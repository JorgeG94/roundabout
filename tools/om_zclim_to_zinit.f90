!! Model-grid z-level T/S climatology -> the `&ocean_zinit_nml source = "file"` layout.
program om_zclim_to_zinit
   !! Converts a z-level ocean climatology that is ALREADY on the model's
   !! horizontal tracer grid — e.g. MOM6's `WOA05_ptemp_salt_annual.nc` style
   !! IC input, `ptemp(time, level, lat, lon)` / `salt(time, level, lat,
   !! lon)`, single time record, land cells and below-bottom levels flagged
   !! by a `_FillValue` attribute — into the layout `rdb_ocean_z_init`
   !! reads: `temp(z, y, x)` / `salt(z, y, x)` / `z_src(z)`, single
   !! precision, no fill values (`tools/om1deg_prepare_inputs.py` builds the
   !! same layout for OM_1deg, but from raw WOA13 on ITS OWN 1-degree grid
   !! with a horizontal nearest-neighbour regrid; this tool does no
   !! horizontal interpolation — the input is already on the target grid,
   !! so only two things need fixing: dropping the leading TIME record and
   !! the fill sentinel).
   !!
   !! Fill handling, per (i, j) column:
   !!  1. a trailing run of fill at the DEEP end (the common case: WOA's
   !!     standard levels stop above the local sea floor) is replaced by the
   !!     deepest valid value above it (constant extrapolation — the same
   !!     rule `interp_column_linear_z` in `rdb_ocean_z_init.F90` applies
   !!     for a model layer deeper than the source column's last level);
   !!  2. a column with NO valid level at all (the source land/sea mask
   !!     disagreeing with the model's, checked for OM4_025's own file: 5
   !!     cells out of 969441 wet) is filled from the nearest column (in a
   !!     growing box) that has data, at every level.
   !! Every column is fully valid on exit; downstream, dry (bathymetry-land)
   !! columns are overwritten by the namelist `land_fill_t`/`land_fill_s`
   !! regardless of what this file carries there, so case 2 only matters for
   !! wet cells.
   !!
   !! Arguments (positional): in.nc out.nc
   !! Options: --temp-var NAME (default ptemp) --salt-var NAME (default salt)
   !!          --level-var NAME (default level)
   use, intrinsic :: iso_fortran_env, only: real32, real64, error_unit, output_unit
   use netcdf, only: nf90_open, nf90_close, nf90_create, nf90_enddef, nf90_noerr, &
                     nf90_nowrite, nf90_clobber, nf90_64bit_offset, nf90_global, &
                     nf90_double, nf90_float, nf90_inq_varid, nf90_inq_dimid, &
                     nf90_inquire_variable, nf90_inquire_dimension, nf90_get_var, &
                     nf90_put_var, nf90_def_dim, nf90_def_var, nf90_put_att, &
                     nf90_get_att, nf90_strerror

   implicit none

   integer, parameter :: sp = real32
   integer, parameter :: dp = real64
   real(dp), parameter :: FILL_DEFAULT = 1.0e19_dp
      !! If the variable carries no `_FillValue` attribute, anything at or
      !! below this (very negative) magnitude is still treated as fill —
      !! WOA-style sentinels are ~1e20/-1e20, never a physical T or S.

   character(len=1024) :: in_path, out_path, arg
   character(len=64) :: temp_var, salt_var, level_var
   integer :: nlon, nlat, nlev, ncid, ncid_out, i

   real(sp), allocatable :: temp(:, :, :), salt(:, :, :)
   real(dp), allocatable :: level(:)
   real(dp) :: fill_t, fill_s

   if (command_argument_count() < 2) then
      write (error_unit, "(a)") "usage: om_zclim_to_zinit in.nc out.nc "// &
         "[--temp-var NAME] [--salt-var NAME] [--level-var NAME]"
      error stop 1
   end if
   call get_command_argument(1, in_path)
   call get_command_argument(2, out_path)
   temp_var = "ptemp"
   salt_var = "salt"
   level_var = "level"
   i = 3
   do while (i <= command_argument_count())
      call get_command_argument(i, arg)
      select case (trim(arg))
      case ("--temp-var")
         i = i + 1
         call get_command_argument(i, temp_var)
      case ("--salt-var")
         i = i + 1
         call get_command_argument(i, salt_var)
      case ("--level-var")
         i = i + 1
         call get_command_argument(i, level_var)
      case default
         write (error_unit, "(a)") "unrecognised option: "//trim(arg)
         error stop 1
      end select
      i = i + 1
   end do

   call check(nf90_open(trim(in_path), nf90_nowrite, ncid), "open "//trim(in_path))
   call read_dims()
   allocate (temp(nlon, nlat, nlev), salt(nlon, nlat, nlev), level(nlev))
   call read_field(trim(temp_var), temp, fill_t)
   call read_field(trim(salt_var), salt, fill_s)
   call read_level()
   call check(nf90_close(ncid), "close "//trim(in_path))

   write (output_unit, "(a,i0,a,i0,a,i0)") "grid ", nlon, " x ", nlat, " x ", nlev
   call fill_column(temp, fill_t, "temp")
   call fill_column(salt, fill_s, "salt")

   call write_out()

contains

   subroutine check(status, what)
      integer, intent(in) :: status
      character(len=*), intent(in) :: what
      if (status /= nf90_noerr) then
         write (error_unit, "(a,a,a,a)") "om_zclim_to_zinit: ", what, ": ", trim(nf90_strerror(status))
         error stop 1
      end if
   end subroutine check

   subroutine read_dims()
      integer :: varid, ndims, n
      integer :: dimids(4)
      call check(nf90_inq_varid(ncid, trim(temp_var), varid), "inq "//trim(temp_var))
      call check(nf90_inquire_variable(ncid, varid, ndims=ndims, dimids=dimids), "inquire "//trim(temp_var))
      if (ndims /= 3 .and. ndims /= 4) then
         write (error_unit, "(a,i0)") "expected a 3-D (level,lat,lon) or 4-D "// &
            "(time,level,lat,lon) source variable, got ndims=", ndims
         error stop 1
      end if
      ! Fortran dimids are fastest-varying first: (lon, lat, level[, time]).
      call check(nf90_inquire_dimension(ncid, dimids(1), len=nlon), "dim1 len")
      call check(nf90_inquire_dimension(ncid, dimids(2), len=nlat), "dim2 len")
      call check(nf90_inquire_dimension(ncid, dimids(3), len=nlev), "dim3 len")
      if (ndims == 4) then
         call check(nf90_inquire_dimension(ncid, dimids(4), len=n), "dim4 len")
         if (n /= 1) then
            write (error_unit, "(a,i0)") "expected a single time record, got ", n
            error stop 1
         end if
      end if
   end subroutine read_dims

   subroutine read_field(vname, arr, fillv)
      character(len=*), intent(in) :: vname
      real(sp), intent(out) :: arr(nlon, nlat, nlev)
      real(dp), intent(out) :: fillv
      integer :: varid, ndims, st
      call check(nf90_inq_varid(ncid, vname, varid), "inq "//vname)
      call check(nf90_inquire_variable(ncid, varid, ndims=ndims), "inquire "//vname)
      if (ndims == 4) then
         call check(nf90_get_var(ncid, varid, arr, start=[1, 1, 1, 1], &
                                 count=[nlon, nlat, nlev, 1]), "get_var "//vname)
      else
         call check(nf90_get_var(ncid, varid, arr, start=[1, 1, 1], &
                                 count=[nlon, nlat, nlev]), "get_var "//vname)
      end if
      st = nf90_get_att(ncid, varid, "_FillValue", fillv)
      if (st /= nf90_noerr) fillv = -FILL_DEFAULT
   end subroutine read_field

   subroutine read_level()
      integer :: varid
      call check(nf90_inq_varid(ncid, trim(level_var), varid), "inq "//trim(level_var))
      call check(nf90_get_var(ncid, varid, level, start=[1], count=[nlev]), "get_var "//trim(level_var))
      do i = 2, nlev
         if (level(i) <= level(i - 1)) then
            write (error_unit, "(a)") "source z-levels not strictly ascending — "// &
               "rdb_ocean_z_init requires monotone z_src"
            error stop 1
         end if
      end do
   end subroutine read_level

   pure function is_fill(v, fillv) result(is_fill_)
      real(sp), intent(in) :: v
      real(dp), intent(in) :: fillv
      logical :: is_fill_
      is_fill_ = (abs(real(v, dp) - fillv) < 1.0e-6_dp*max(1.0_dp, abs(fillv))) .or. &
                 (abs(real(v, dp)) >= FILL_DEFAULT)
   end function is_fill

   subroutine fill_column(arr, fillv, label)
      real(sp), intent(inout) :: arr(nlon, nlat, nlev)
      real(dp), intent(in) :: fillv
      character(len=*), intent(in) :: label
      integer :: ii, jj, k, k0, di, dj, ni, nj, radius
      integer :: n_trailing, n_neighbor_fallback
      logical :: found

      n_trailing = 0
      n_neighbor_fallback = 0

      ! Pass 1: trailing (deep-end) fill -> constant extrapolation from the
      ! deepest valid level above it.
      do jj = 1, nlat
         do ii = 1, nlon
            if (is_fill(arr(ii, jj, 1), fillv)) cycle  ! handled in pass 2
            k0 = 1
            do k = 2, nlev
               if (.not. is_fill(arr(ii, jj, k), fillv)) then
                  k0 = k
               else
                  exit
               end if
            end do
            if (k0 < nlev) then
               do k = k0 + 1, nlev
                  if (is_fill(arr(ii, jj, k), fillv)) then
                     arr(ii, jj, k) = arr(ii, jj, k0)
                     n_trailing = n_trailing + 1
                  end if
               end do
            end if
         end do
      end do

      ! Pass 2: whole-column fill (surface itself missing) -> nearest valid
      ! column, searched as a growing square box (periodic in longitude).
      do jj = 1, nlat
         do ii = 1, nlon
            if (.not. is_fill(arr(ii, jj, 1), fillv)) cycle
            found = .false.
            search: do radius = 1, max(nlon, nlat)
               do dj = -radius, radius
                  nj = jj + dj
                  if (nj < 1 .or. nj > nlat) cycle
                  do di = -radius, radius
                     if (abs(di) /= radius .and. abs(dj) /= radius) cycle
                     ni = modulo(ii - 1 + di, nlon) + 1
                     if (.not. is_fill(arr(ni, nj, 1), fillv)) then
                        arr(ii, jj, :) = arr(ni, nj, :)
                        found = .true.
                        n_neighbor_fallback = n_neighbor_fallback + 1
                        exit search
                     end if
                  end do
               end do
            end do search
            if (.not. found) then
               write (error_unit, "(a,a,a,i0,a,i0)") "om_zclim_to_zinit: no valid ", &
                  label, " column found anywhere for i=", ii, " j=", jj
               error stop 1
            end if
         end do
      end do

      write (output_unit, "(a,a,a,i0,a,i0)") label, ": trailing-fill levels replaced=", &
         "", n_trailing, "  whole-column fallbacks=", n_neighbor_fallback
   end subroutine fill_column

   subroutine write_out()
      integer :: d_x, d_y, d_z, v_x, v_y, v_z, v_t, v_s
      call check(nf90_create(trim(out_path), ior(nf90_clobber, nf90_64bit_offset), ncid_out), &
                 "create "//trim(out_path))
      call check(nf90_def_dim(ncid_out, "z", nlev, d_z), "def z")
      call check(nf90_def_dim(ncid_out, "y", nlat, d_y), "def y")
      call check(nf90_def_dim(ncid_out, "x", nlon, d_x), "def x")
      call check(nf90_def_var(ncid_out, "z_src", nf90_double, [d_z], v_z), "def z_src")
      call check(nf90_put_att(ncid_out, v_z, "units", "m"), "att z_src units")
      call check(nf90_put_att(ncid_out, v_z, "positive", "down"), "att z_src positive")
      call check(nf90_def_var(ncid_out, "temp", nf90_float, [d_x, d_y, d_z], v_t), "def temp")
      call check(nf90_put_att(ncid_out, v_t, "units", "degC"), "att temp units")
      call check(nf90_def_var(ncid_out, "salt", nf90_float, [d_x, d_y, d_z], v_s), "def salt")
      call check(nf90_put_att(ncid_out, v_s, "units", "psu"), "att salt units")
      call check(nf90_put_att(ncid_out, nf90_global, "source", &
                              "tools/om_zclim_to_zinit.f90: "//trim(in_path)// &
                              " -> the ocean_zinit file/z/y/x layout, fill-extended, no horizontal regrid"), &
                 "att global source")
      call check(nf90_enddef(ncid_out), "enddef")
      call check(nf90_put_var(ncid_out, v_z, level), "put z_src")
      call check(nf90_put_var(ncid_out, v_t, temp), "put temp")
      call check(nf90_put_var(ncid_out, v_s, salt), "put salt")
      call check(nf90_close(ncid_out), "close "//trim(out_path))
      write (output_unit, "(a,a)") "wrote ", trim(out_path)
   end subroutine write_out

end program om_zclim_to_zinit
