!! Generic NetCDF dimension-range subsetter for the OM_1deg inputs.
program om1deg_subset
   !! Cuts a rectangular sub-range out of one or more named dimensions of a
   !! NetCDF file, copying every dimension, variable, and attribute (global
   !! and per-variable) unchanged apart from the sliced extent. Used to carve
   !! the Southern Ocean domain (`validation_examples/ocean/southern_ocean_1deg/
   !! README.md`) out of the global OM_1deg supergrid, bathymetry, WOA13
   !! initial condition and JRA55-do wind-stress files: every one of those
   !! files stores its fields on a plain `(y, x[, z|time])`-style grid south
   !! of the Arctic tripolar cap, so a single contiguous slice along the `y`
   !! (or supergrid `ny`/`nyp`) index is all four inputs need. Not part of
   !! the model; a data-preparation helper in the pattern of
   !! `tools/om1deg_wind_regrid.f90`.
   !!
   !! Arguments (positional):
   !!   in.nc out.nc [dimname:start0:count ...]
   !!
   !! `start0` is the 0-based first index kept, `count` the number of
   !! indices kept. A dimension named on the command line is truncated to
   !! that range; every other dimension (and every variable that does not
   !! carry a truncated dimension) is copied in full. Handles the `NF90_CHAR`
   !! / `NF90_FLOAT` / `NF90_DOUBLE` variable ranks 1-3 that these four input
   !! files use; anything else is a fail-loud `error stop` naming the
   !! offending variable, not a silent skip.
   use, intrinsic :: iso_fortran_env, only: real32, real64, error_unit, output_unit
   use netcdf, only: nf90_open, nf90_close, nf90_create, nf90_enddef, nf90_noerr, &
                     nf90_nowrite, nf90_clobber, nf90_64bit_offset, nf90_global, &
                     nf90_double, nf90_float, nf90_char, nf90_max_name, &
                     nf90_inquire, nf90_inquire_dimension, nf90_inquire_variable, &
                     nf90_inq_attname, nf90_copy_att, nf90_def_dim, nf90_def_var, &
                     nf90_get_var, nf90_put_var, nf90_put_att, nf90_strerror

   implicit none

   integer, parameter :: dp = real64
   integer, parameter :: sp = real32
   integer, parameter :: MAX_SLICES = 16
   integer, parameter :: MAX_DIMS = 32
   integer, parameter :: MAX_VDIMS = 8

   character(len=1024) :: in_path, out_path, arg
   character(len=64) :: slice_name(MAX_SLICES)
   integer :: slice_start(MAX_SLICES), slice_count(MAX_SLICES)
   integer :: n_slices, nargs, ia, p1, p2

   integer :: ncid_in, ncid_out
   integer :: ndims, nvars, ngatts, unlimdimid
   integer :: d, v, a, vndims, vtype, vnatts
   character(len=nf90_max_name) :: dname, vname, attname
   integer :: dlen
   integer, allocatable :: out_dimid(:), in_dim_start0(:), in_dim_count(:)
   integer, allocatable :: vdimids(:)
   integer :: out_vdimids(MAX_VDIMS)
   integer :: varid_out, k

   nargs = command_argument_count()
   if (nargs < 2) then
      write (error_unit, "(a)") "usage: om1deg_subset in.nc out.nc "// &
         "[dimname:start0:count ...]"
      stop 1
   end if

   call get_command_argument(1, in_path)
   call get_command_argument(2, out_path)

   n_slices = 0
   do ia = 3, nargs
      call get_command_argument(ia, arg)
      p1 = index(arg, ":")
      p2 = index(arg, ":", back=.true.)
      if (p1 == 0 .or. p2 == p1) then
         write (error_unit, "(a)") "om1deg_subset: bad slice spec '"//trim(arg)// &
            "' (want name:start0:count)"
         stop 1
      end if
      n_slices = n_slices + 1
      if (n_slices > MAX_SLICES) then
         write (error_unit, "(a)") "om1deg_subset: too many slice specs"
         stop 1
      end if
      slice_name(n_slices) = arg(1:p1 - 1)
      read (arg(p1 + 1:p2 - 1), *) slice_start(n_slices)
      read (arg(p2 + 1:), *) slice_count(n_slices)
   end do

   call check(nf90_open(trim(in_path), nf90_nowrite, ncid_in), "open "//trim(in_path))
   call check(nf90_inquire(ncid_in, ndims, nvars, ngatts, unlimdimid), "inquire")
   if (ndims > MAX_DIMS) then
      write (error_unit, "(a)") "om1deg_subset: MAX_DIMS too small"
      stop 1
   end if

   allocate (out_dimid(ndims), in_dim_start0(ndims), in_dim_count(ndims))

   call check(nf90_create(trim(out_path), ior(nf90_clobber, nf90_64bit_offset), ncid_out), &
              "create "//trim(out_path))

   ! ---- dimensions: resolve each against the slice specs, define in output
   do d = 1, ndims
      call check(nf90_inquire_dimension(ncid_in, d, dname, dlen), "inquire_dimension")
      in_dim_start0(d) = 0
      in_dim_count(d) = dlen
      do k = 1, n_slices
         if (trim(slice_name(k)) == trim(dname)) then
            if (slice_start(k) < 0 .or. slice_start(k) + slice_count(k) > dlen) then
               write (error_unit, "(a)") "om1deg_subset: slice '"//trim(dname)// &
                  "' out of range for input length"
               stop 1
            end if
            in_dim_start0(d) = slice_start(k)
            in_dim_count(d) = slice_count(k)
            exit
         end if
      end do
      call check(nf90_def_dim(ncid_out, trim(dname), in_dim_count(d), out_dimid(d)), &
                 "def_dim "//trim(dname))
      write (output_unit, "(a,a,a,i0,a,i0,a,i0)") "  dim ", trim(dname), ": in=", dlen, &
         " -> start0=", in_dim_start0(d), " count=", in_dim_count(d)
   end do

   ! ---- global attributes
   do a = 1, ngatts
      call check(nf90_inq_attname(ncid_in, nf90_global, a, attname), "inq_attname global")
      call check(nf90_copy_att(ncid_in, nf90_global, trim(attname), ncid_out, nf90_global), &
                 "copy_att global "//trim(attname))
   end do
   call check(nf90_put_att(ncid_out, nf90_global, "southern_ocean_1deg_subset", &
                           "tools/om1deg_subset.f90: rectangular y-range cut of the "// &
                           "OM_1deg global file "//trim(in_path)), "put_att provenance")

   ! ---- variables: define + copy attributes (data copied in a second pass,
   !      after enddef, so def_var never runs interleaved with get/put_var)
   allocate (vdimids(MAX_VDIMS))
   do v = 1, nvars
      call check(nf90_inquire_variable(ncid_in, v, vname, vtype, vndims, vdimids, vnatts), &
                 "inquire_variable")
      if (vndims > MAX_VDIMS) then
         write (error_unit, "(a)") "om1deg_subset: MAX_VDIMS too small for "//trim(vname)
         stop 1
      end if
      do k = 1, vndims
         out_vdimids(k) = out_dimid(vdimids(k))
      end do
      call check(nf90_def_var(ncid_out, trim(vname), vtype, out_vdimids(1:vndims), varid_out), &
                 "def_var "//trim(vname))
      do a = 1, vnatts
         call check(nf90_inq_attname(ncid_in, v, a, attname), "inq_attname "//trim(vname))
         call check(nf90_copy_att(ncid_in, v, trim(attname), ncid_out, varid_out), &
                    "copy_att "//trim(vname)//" "//trim(attname))
      end do
   end do

   call check(nf90_enddef(ncid_out), "enddef")

   ! ---- data, one variable at a time
   do v = 1, nvars
      call check(nf90_inquire_variable(ncid_in, v, vname, vtype, vndims, vdimids, vnatts), &
                 "inquire_variable (data pass)")
      write (output_unit, "(a,a)") "  copying ", trim(vname)
      select case (vtype)
      case (nf90_double)
         select case (vndims)
         case (1)
            call copy_r8_1d(ncid_in, v, ncid_out, v, vdimids, in_dim_start0, in_dim_count)
         case (2)
            call copy_r8_2d(ncid_in, v, ncid_out, v, vdimids, in_dim_start0, in_dim_count)
         case (3)
            call copy_r8_3d(ncid_in, v, ncid_out, v, vdimids, in_dim_start0, in_dim_count)
         case default
            write (error_unit, "(a,a,a,i0)") "om1deg_subset: unsupported rank for ", &
               trim(vname), " ndims=", vndims
            stop 1
         end select
      case (nf90_float)
         select case (vndims)
         case (1)
            call copy_r4_1d(ncid_in, v, ncid_out, v, vdimids, in_dim_start0, in_dim_count)
         case (2)
            call copy_r4_2d(ncid_in, v, ncid_out, v, vdimids, in_dim_start0, in_dim_count)
         case (3)
            call copy_r4_3d(ncid_in, v, ncid_out, v, vdimids, in_dim_start0, in_dim_count)
         case default
            write (error_unit, "(a,a,a,i0)") "om1deg_subset: unsupported rank for ", &
               trim(vname), " ndims=", vndims
            stop 1
         end select
      case (nf90_char)
         if (vndims /= 1) then
            write (error_unit, "(a,a)") "om1deg_subset: only rank-1 char supported: ", &
               trim(vname)
            stop 1
         end if
         call copy_char_1d(ncid_in, v, ncid_out, v, vdimids, in_dim_start0, in_dim_count)
      case default
         write (error_unit, "(a,a)") "om1deg_subset: unsupported NetCDF type for ", &
            trim(vname)
         stop 1
      end select
   end do

   call check(nf90_close(ncid_in), "close in")
   call check(nf90_close(ncid_out), "close out")
   write (output_unit, "(a)") "om1deg_subset: wrote "//trim(out_path)

contains

   subroutine check(status, what)
      !! Abort with the NetCDF message on any error.
      integer, intent(in) :: status
      character(len=*), intent(in) :: what
      if (status /= nf90_noerr) then
         write (error_unit, "(a)") "om1deg_subset: "//what//": "//trim(nf90_strerror(status))
         stop 1
      end if
   end subroutine check

   subroutine copy_r8_1d(ncin, vin, ncout, vout, dimids, start0, cnt)
      integer, intent(in) :: ncin, vin, ncout, vout
      integer, intent(in) :: dimids(:), start0(:), cnt(:)
      real(dp), allocatable :: buf(:)
      integer :: st(1), ct(1)
      st(1) = start0(dimids(1)) + 1
      ct(1) = cnt(dimids(1))
      allocate (buf(ct(1)))
      call check(nf90_get_var(ncin, vin, buf, start=st, count=ct), "get_var r8_1d")
      call check(nf90_put_var(ncout, vout, buf), "put_var r8_1d")
      deallocate (buf)
   end subroutine copy_r8_1d

   subroutine copy_r8_2d(ncin, vin, ncout, vout, dimids, start0, cnt)
      integer, intent(in) :: ncin, vin, ncout, vout
      integer, intent(in) :: dimids(:), start0(:), cnt(:)
      real(dp), allocatable :: buf(:, :)
      integer :: st(2), ct(2)
      st = [start0(dimids(1)) + 1, start0(dimids(2)) + 1]
      ct = [cnt(dimids(1)), cnt(dimids(2))]
      allocate (buf(ct(1), ct(2)))
      call check(nf90_get_var(ncin, vin, buf, start=st, count=ct), "get_var r8_2d")
      call check(nf90_put_var(ncout, vout, buf), "put_var r8_2d")
      deallocate (buf)
   end subroutine copy_r8_2d

   subroutine copy_r8_3d(ncin, vin, ncout, vout, dimids, start0, cnt)
      integer, intent(in) :: ncin, vin, ncout, vout
      integer, intent(in) :: dimids(:), start0(:), cnt(:)
      real(dp), allocatable :: buf(:, :, :)
      integer :: st(3), ct(3)
      st = [start0(dimids(1)) + 1, start0(dimids(2)) + 1, start0(dimids(3)) + 1]
      ct = [cnt(dimids(1)), cnt(dimids(2)), cnt(dimids(3))]
      allocate (buf(ct(1), ct(2), ct(3)))
      call check(nf90_get_var(ncin, vin, buf, start=st, count=ct), "get_var r8_3d")
      call check(nf90_put_var(ncout, vout, buf), "put_var r8_3d")
      deallocate (buf)
   end subroutine copy_r8_3d

   subroutine copy_r4_1d(ncin, vin, ncout, vout, dimids, start0, cnt)
      integer, intent(in) :: ncin, vin, ncout, vout
      integer, intent(in) :: dimids(:), start0(:), cnt(:)
      real(sp), allocatable :: buf(:)
      integer :: st(1), ct(1)
      st(1) = start0(dimids(1)) + 1
      ct(1) = cnt(dimids(1))
      allocate (buf(ct(1)))
      call check(nf90_get_var(ncin, vin, buf, start=st, count=ct), "get_var r4_1d")
      call check(nf90_put_var(ncout, vout, buf), "put_var r4_1d")
      deallocate (buf)
   end subroutine copy_r4_1d

   subroutine copy_r4_2d(ncin, vin, ncout, vout, dimids, start0, cnt)
      integer, intent(in) :: ncin, vin, ncout, vout
      integer, intent(in) :: dimids(:), start0(:), cnt(:)
      real(sp), allocatable :: buf(:, :)
      integer :: st(2), ct(2)
      st = [start0(dimids(1)) + 1, start0(dimids(2)) + 1]
      ct = [cnt(dimids(1)), cnt(dimids(2))]
      allocate (buf(ct(1), ct(2)))
      call check(nf90_get_var(ncin, vin, buf, start=st, count=ct), "get_var r4_2d")
      call check(nf90_put_var(ncout, vout, buf), "put_var r4_2d")
      deallocate (buf)
   end subroutine copy_r4_2d

   subroutine copy_r4_3d(ncin, vin, ncout, vout, dimids, start0, cnt)
      integer, intent(in) :: ncin, vin, ncout, vout
      integer, intent(in) :: dimids(:), start0(:), cnt(:)
      real(sp), allocatable :: buf(:, :, :)
      integer :: st(3), ct(3)
      st = [start0(dimids(1)) + 1, start0(dimids(2)) + 1, start0(dimids(3)) + 1]
      ct = [cnt(dimids(1)), cnt(dimids(2)), cnt(dimids(3))]
      allocate (buf(ct(1), ct(2), ct(3)))
      call check(nf90_get_var(ncin, vin, buf, start=st, count=ct), "get_var r4_3d")
      call check(nf90_put_var(ncout, vout, buf), "put_var r4_3d")
      deallocate (buf)
   end subroutine copy_r4_3d

   subroutine copy_char_1d(ncin, vin, ncout, vout, dimids, start0, cnt)
      integer, intent(in) :: ncin, vin, ncout, vout
      integer, intent(in) :: dimids(:), start0(:), cnt(:)
      character(len=:), allocatable :: buf
      integer :: st(1), ct(1)
      integer :: n
      st(1) = start0(dimids(1)) + 1
      ct(1) = cnt(dimids(1))
      n = ct(1)
      allocate (character(len=n) :: buf)
      call check(nf90_get_var(ncin, vin, buf, start=st, count=ct), "get_var char_1d")
      call check(nf90_put_var(ncout, vout, buf), "put_var char_1d")
      deallocate (buf)
   end subroutine copy_char_1d

end program om1deg_subset
