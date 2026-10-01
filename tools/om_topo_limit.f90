!! Apply MOM6-style MINIMUM_DEPTH/MAXIMUM_DEPTH/MASKING_DEPTH limits, in place.
program om_topo_limit
   !! MOM6's `limit_topography` (`MOM_shared_initialization.F90`): every cell
   !! with `depth > masking_depth` is wet, raised to at least `min_depth` and
   !! capped at `max_depth`; every other cell is land, `depth = 0`. A raw
   !! `topog.nc`-style file (`depth`, `wet`) — unlike the already-limited
   !! `bathy_om1deg.nc` `tools/om1deg_prepare_inputs.py` builds for OM_1deg —
   !! needs this applied before `&ocean_topo_nml topo_config = "file"` reads
   !! it: `roundabout`'s bathymetry reader takes `depth` as given, with no
   !! floor/cap of its own. Modifies `depth`/`wet` IN PLACE (same file, same
   !! shape); every other variable is untouched.
   !!
   !! Arguments (positional): file.nc min_depth max_depth masking_depth
   !! Options: --depth-var NAME (default depth) --wet-var NAME (default wet,
   !!          skipped if absent), --rename-dim OLD:NEW (repeatable — e.g. a
   !!          raw MOM6 `topog.nc`'s `nx`/`ny` to the `x`/`y` roundabout's
   !!          `rdb_bathymetry` reader hardcodes; applied BEFORE the limit,
   !!          renaming is a define-mode op even to a shorter name)
   use, intrinsic :: iso_fortran_env, only: real32, real64, error_unit, output_unit
   use netcdf, only: nf90_open, nf90_close, nf90_write, nf90_noerr, nf90_inq_varid, &
                     nf90_inquire_variable, nf90_get_var, nf90_put_var, nf90_strerror, &
                     nf90_inq_dimid, nf90_redef, nf90_enddef, nf90_rename_dim

   implicit none

   integer, parameter :: sp = real32
   integer, parameter :: dp = real64
   integer, parameter :: MAX_RENAMES = 8

   character(len=1024) :: fpath, arg
   character(len=64) :: depth_var, wet_var
   character(len=64) :: rename_old(MAX_RENAMES), rename_new(MAX_RENAMES)
   integer :: n_renames
   real(dp) :: min_depth, max_depth, masking_depth
   integer :: ncid, varid_d, varid_w, st, ndims, nx, ny, i, nwet, nland, nraised, ncapped
   integer :: dimid, colon
   integer :: dimids(2)

   real(sp), allocatable :: depth(:, :), wet(:, :)
   logical :: has_wet

   if (command_argument_count() < 4) then
      write (error_unit, "(a)") "usage: om_topo_limit file.nc min_depth max_depth masking_depth "// &
         "[--depth-var NAME] [--wet-var NAME] [--rename-dim OLD:NEW ...]"
      error stop 1
   end if
   call get_command_argument(1, fpath)
   min_depth = real_arg(2, "min_depth")
   max_depth = real_arg(3, "max_depth")
   masking_depth = real_arg(4, "masking_depth")
   depth_var = "depth"
   wet_var = "wet"
   n_renames = 0
   i = 5
   do while (i <= command_argument_count())
      call get_command_argument(i, arg)
      select case (trim(arg))
      case ("--depth-var")
         i = i + 1
         call need_value(i, "--depth-var")
         call get_command_argument(i, depth_var)
      case ("--wet-var")
         i = i + 1
         call need_value(i, "--wet-var")
         call get_command_argument(i, wet_var)
      case ("--rename-dim")
         i = i + 1
         call need_value(i, "--rename-dim")
         n_renames = n_renames + 1
         if (n_renames > MAX_RENAMES) then
            write (error_unit, "(a)") "om_topo_limit: too many --rename-dim options"
            error stop 1
         end if
         call get_command_argument(i, arg)
         colon = index(arg, ":")
         if (colon < 2) then
            write (error_unit, "(a)") "om_topo_limit: --rename-dim needs OLD:NEW, got "//trim(arg)
            error stop 1
         end if
         rename_old(n_renames) = arg(1:colon - 1)
         rename_new(n_renames) = arg(colon + 1:)
      case default
         write (error_unit, "(a)") "unrecognised option: "//trim(arg)
         error stop 1
      end select
      i = i + 1
   end do

   call check(nf90_open(trim(fpath), nf90_write, ncid), "open "//trim(fpath))

   if (n_renames > 0) then
      call check(nf90_redef(ncid), "redef")
      do i = 1, n_renames
         call check(nf90_inq_dimid(ncid, trim(rename_old(i)), dimid), "inq_dimid "//trim(rename_old(i)))
         call check(nf90_rename_dim(ncid, dimid, trim(rename_new(i))), &
                    "rename_dim "//trim(rename_old(i))//" -> "//trim(rename_new(i)))
         write (output_unit, "(a,a,a,a)") "  renamed dim ", trim(rename_old(i)), " -> ", trim(rename_new(i))
      end do
      call check(nf90_enddef(ncid), "enddef after rename")
   end if
   call check(nf90_inq_varid(ncid, trim(depth_var), varid_d), "inq "//trim(depth_var))
   call check(nf90_inquire_variable(ncid, varid_d, ndims=ndims, dimids=dimids), "inquire "//trim(depth_var))
   if (ndims /= 2) then
      write (error_unit, "(a)") "om_topo_limit: only a rank-2 depth field is supported"
      error stop 1
   end if
   call get_dim_lens(ncid, dimids, nx, ny)
   allocate (depth(nx, ny), wet(nx, ny))
   call check(nf90_get_var(ncid, varid_d, depth), "get_var "//trim(depth_var))

   st = nf90_inq_varid(ncid, trim(wet_var), varid_w)
   has_wet = (st == nf90_noerr)
   if (has_wet) then
      call check(nf90_get_var(ncid, varid_w, wet), "get_var "//trim(wet_var))
   end if

   nwet = 0
   nland = 0
   nraised = 0
   ncapped = 0
   block
      integer :: ii, jj
      do jj = 1, ny
         do ii = 1, nx
            if (real(depth(ii, jj), dp) > masking_depth) then
               nwet = nwet + 1
               if (real(depth(ii, jj), dp) < min_depth) then
                  depth(ii, jj) = real(min_depth, sp)
                  nraised = nraised + 1
               else if (real(depth(ii, jj), dp) > max_depth) then
                  depth(ii, jj) = real(max_depth, sp)
                  ncapped = ncapped + 1
               end if
               if (has_wet) wet(ii, jj) = 1.0_sp
            else
               nland = nland + 1
               depth(ii, jj) = 0.0_sp
               if (has_wet) wet(ii, jj) = 0.0_sp
            end if
         end do
      end do
   end block

   call check(nf90_put_var(ncid, varid_d, depth), "put_var "//trim(depth_var))
   if (has_wet) call check(nf90_put_var(ncid, varid_w, wet), "put_var "//trim(wet_var))
   call check(nf90_close(ncid), "close "//trim(fpath))

   write (output_unit, "(a,a,a,f0.1,a,f0.1,a,f0.1)") "om_topo_limit: ", trim(fpath), &
      "  min=", min_depth, " max=", max_depth, " masking=", masking_depth
   write (output_unit, "(a,i0,a,i0,a,i0,a,i0)") "  wet=", nwet, " land=", nland, &
      " raised-to-min=", nraised, " capped-to-max=", ncapped

contains

   real(dp) function real_arg(pos, what) result(val)
      !! Positional argument `pos` read as a real; fails loud on a value that
      !! does not parse instead of crashing in the list-directed read.
      integer, intent(in) :: pos
      character(len=*), intent(in) :: what
      character(len=256) :: txt, msg
      integer :: ios
      call get_command_argument(pos, txt)
      read (txt, *, iostat=ios, iomsg=msg) val
      if (ios /= 0) then
         write (error_unit, "(a)") "om_topo_limit: "//what//" is not a number: '"//trim(txt)// &
            "' ("//trim(msg)//")"
         error stop 1
      end if
   end function real_arg

   subroutine need_value(pos, opt)
      !! Fail loud when option `opt` is the last argument (no value at `pos`).
      integer, intent(in) :: pos
      character(len=*), intent(in) :: opt
      if (pos > command_argument_count()) then
         write (error_unit, "(a)") "om_topo_limit: "//opt//" needs a value"
         error stop 1
      end if
   end subroutine need_value

   subroutine check(status, what)
      integer, intent(in) :: status
      character(len=*), intent(in) :: what
      if (status /= nf90_noerr) then
         write (error_unit, "(a,a,a,a)") "om_topo_limit: ", what, ": ", trim(nf90_strerror(status))
         error stop 1
      end if
   end subroutine check

   subroutine get_dim_lens(nc, dimids_in, nx_out, ny_out)
      use netcdf, only: nf90_inquire_dimension
      integer, intent(in) :: nc
      integer, intent(in) :: dimids_in(2)
      integer, intent(out) :: nx_out, ny_out
      ! Fortran dimids are fastest-varying first: (x, y).
      call check(nf90_inquire_dimension(nc, dimids_in(1), len=nx_out), "dim1 len")
      call check(nf90_inquire_dimension(nc, dimids_in(2), len=ny_out), "dim2 len")
   end subroutine get_dim_lens

end program om_topo_limit
