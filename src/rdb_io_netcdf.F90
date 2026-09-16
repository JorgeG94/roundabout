!! Thin wrapper around netcdf-fortran (nf90_*) with error checking
module rdb_io_netcdf
   !! Provides NetCDF read/write primitives and an error-checking helper.
   !! Used by bathymetry loader and output writer.
   use, intrinsic :: iso_fortran_env, only: real32
   use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char, c_ptr, c_loc
   use netcdf, only: nf90_noerr, nf90_strerror, nf90_create, nf90_open, &
                     nf90_close, nf90_def_dim, nf90_def_var, nf90_put_att, nf90_enddef, &
                     nf90_put_var, nf90_get_var, nf90_inq_dimid, nf90_inquire_dimension, &
                     nf90_inq_varid, nf90_get_att, nf90_def_var_deflate, &
                     nf90_clobber, nf90_netcdf4, nf90_nowrite, nf90_write, &
                     nf90_unlimited, nf90_double, nf90_float, nf90_int, nf90_global, &
                     nf90_redef, nf90_sync
   use rdb_constants, only: wp
   use pic_logger, only: logger => global_logger
   use pic_strings, only: to_string
   use rdb_error_ring, only: fail
   implicit none
   private

   public :: nc_check
   public :: ensure_directory_exists
   public :: nc_create_file
   public :: nc_open_read
   public :: nc_open_write
   public :: nc_close
   public :: nc_def_dim
   public :: nc_def_var_2d
   public :: nc_def_var_3d
   public :: nc_put_att
   public :: nc_put_att_real
   public :: nc_put_att_int
   public :: nc_put_att_global
   public :: nc_put_att_global_int
   public :: nc_get_att_int
   public :: nc_enddef
   public :: nc_put_var_2d
   public :: nc_put_var_3d_slice
   public :: nc_get_var_2d
   public :: nc_get_var_3d
   public :: nc_get_var_1d
   public :: nc_get_dim_len
   public :: nc_get_varid
   public :: rdb_def_var_1d
   public :: rdb_put_var_1d
   public :: rdb_put_var_2d_slice
   public :: nc_def_var_4d
   public :: nc_put_var_4d_slice
   public :: nc_get_var_slab_3d
   public :: nc_get_att_text
   public :: nc_put_att_real_r4
   public :: nc_put_var_3d_slice_r4
   public :: nc_put_var_4d_slice_r4
   public :: NC_WP
   public :: NC_R4
   public :: output_rank_filename

#ifdef RDB_DOUBLE_PRECISION
   integer, parameter :: NC_WP = nf90_double
      !! NetCDF type matching working precision
#else
   integer, parameter :: NC_WP = nf90_float
      !! NetCDF type matching working precision
#endif

   integer, parameter :: NC_R4 = nf90_float
      !! Explicit 32-bit NetCDF element type.  Used ONLY by the ocean
      !! diagnostic stream when `&ocean_diag_nml output_precision="single"`
      !! asks for a halved-width output file.  Restarts, gauges, coastal
      !! output and console/conservation totals stay on `NC_WP` — a lossy
      !! restart must not be requestable.  The 64-bit counterpart is
      !! `NC_WP` itself — deliberately NOT given an `NC_R8` alias, so
      !! "double" has exactly one spelling and cannot drift away from the
      !! build's working precision.

contains

   subroutine nc_check(status, context, ierr)
      !! Check NetCDF return status and log error if it failed
      integer, intent(in) :: status
         !! Return code from nf90_* call
      character(len=*), intent(in) :: context
         !! Description of what was being attempted
      integer, intent(out), optional :: ierr
         !! Non-zero (the raw `nf90_*` status) on failure when present;
         !! absent behaves as today (`error stop`).  The specific reason
         !! is always logged via `global_logger%error` first, so the
         !! caller-facing code only needs to say THAT a read/write failed.

      if (status /= nf90_noerr) then
         call fail("NetCDF error in "//trim(context)//": "// &
                   trim(nf90_strerror(status)), ierr, status)
         return
      end if

      if (present(ierr)) ierr = nf90_noerr

   end subroutine nc_check

   subroutine nc_create_file(filename, ncid, ierr)
      !! Create a new NetCDF-4 file (overwrites if exists)
      character(len=*), intent(in) :: filename
      integer, intent(out) :: ncid
      integer, intent(out), optional :: ierr
         !! Non-zero (e.g. a missing output directory) when present; absent
         !! behaves as today (`error stop`).

      call nc_check(nf90_create(trim(filename), &
                                ior(nf90_clobber, nf90_netcdf4), ncid), &
                    "creating "//trim(filename), ierr)

   end subroutine nc_create_file

   subroutine nc_open_read(filename, ncid, ierr)
      !! Open an existing NetCDF file for reading
      character(len=*), intent(in) :: filename
      integer, intent(out) :: ncid
      integer, intent(out), optional :: ierr
         !! Non-zero (a missing/unreadable file) when present; absent
         !! behaves as today (`error stop`).

      call nc_check(nf90_open(trim(filename), nf90_nowrite, ncid), &
                    "opening "//trim(filename), ierr)

   end subroutine nc_open_read

   subroutine nc_open_write(filename, ncid, ierr)
      !! Open an existing NetCDF file for writing
      character(len=*), intent(in) :: filename
      integer, intent(out) :: ncid
      integer, intent(out), optional :: ierr
         !! Non-zero on failure when present; absent behaves as today
         !! (`error stop`).

      call nc_check(nf90_open(trim(filename), nf90_write, ncid), &
                    "opening for write "//trim(filename), ierr)

   end subroutine nc_open_write

   subroutine nc_close(ncid, ierr)
      !! Close a NetCDF file
      integer, intent(in) :: ncid
      integer, intent(out), optional :: ierr
         !! Non-zero on failure when present; absent behaves as today
         !! (`error stop`).

      call nc_check(nf90_close(ncid), "closing file", ierr)

   end subroutine nc_close

   subroutine nc_def_dim(ncid, name, length, dimid, ierr)
      !! Define a dimension (use nf90_unlimited for unlimited)
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: name
      integer, intent(in) :: length
      integer, intent(out) :: dimid
      integer, intent(out), optional :: ierr
         !! Non-zero on failure when present; absent behaves as today
         !! (`error stop`).

      call nc_check(nf90_def_dim(ncid, name, length, dimid), &
                    "defining dimension "//trim(name), ierr)

   end subroutine nc_def_dim

   subroutine rdb_def_var_1d(ncid, name, dimid, varid, deflate_level, ierr)
      !! Define a 1D variable with working precision type
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: name
      integer, intent(in) :: dimid
      integer, intent(out) :: varid
      integer, intent(in), optional :: deflate_level
         !! Compression level (0=none, 1-9=deflate). Enables shuffle when > 0.
      integer, intent(out), optional :: ierr
         !! Non-zero on failure when present; absent behaves as today
         !! (`error stop`).

      call nc_check(nf90_def_var(ncid, name, NC_WP, [dimid], varid), &
                    "defining variable "//trim(name), ierr)
      if (present(ierr)) then
         if (ierr /= nf90_noerr) return
      end if

      if (present(deflate_level)) then
         if (deflate_level > 0) then
            call nc_check(nf90_def_var_deflate(ncid, varid, &
                                               shuffle=1, deflate=1, deflate_level=deflate_level), &
                          "enabling compression for "//trim(name), ierr)
         end if
      end if

   end subroutine rdb_def_var_1d

   subroutine rdb_put_var_1d(ncid, varid, data, ierr)
      !! Write a full 1D array
      integer, intent(in) :: ncid
      integer, intent(in) :: varid
      real(wp), intent(in) :: data(:)
      integer, intent(out), optional :: ierr
         !! Non-zero on failure when present; absent behaves as today
         !! (`error stop`).

      call nc_check(nf90_put_var(ncid, varid, data), "writing 1D variable", ierr)

   end subroutine rdb_put_var_1d

   subroutine rdb_put_var_2d_slice(ncid, varid, data, time_index, ierr)
      !! Write a 1D slice into a 2D variable at a given time index
      integer, intent(in) :: ncid
      integer, intent(in) :: varid
      real(wp), intent(in) :: data(:)
      integer, intent(in) :: time_index
      integer, intent(out), optional :: ierr
         !! Non-zero on failure when present; absent behaves as today
         !! (`error stop`).

      call nc_check(nf90_put_var(ncid, varid, data, &
                                 start=[1, time_index], &
                                 count=[size(data), 1]), &
                    "writing 2D slice", ierr)

   end subroutine rdb_put_var_2d_slice

   subroutine nc_def_var_2d(ncid, name, dimids, varid, deflate_level, ierr)
      !! Define a 2D variable with working precision type
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: name
      integer, intent(in) :: dimids(2)
      integer, intent(out) :: varid
      integer, intent(in), optional :: deflate_level
         !! Compression level (0=none, 1-9=deflate). Enables shuffle when > 0.
      integer, intent(out), optional :: ierr
         !! Non-zero on failure when present; absent behaves as today
         !! (`error stop`).

      call nc_check(nf90_def_var(ncid, name, NC_WP, dimids, varid), &
                    "defining variable "//trim(name), ierr)
      if (present(ierr)) then
         if (ierr /= nf90_noerr) return
      end if

      if (present(deflate_level)) then
         if (deflate_level > 0) then
            call nc_check(nf90_def_var_deflate(ncid, varid, &
                                               shuffle=1, deflate=1, deflate_level=deflate_level), &
                          "enabling compression for "//trim(name), ierr)
         end if
      end if

   end subroutine nc_def_var_2d

   subroutine nc_def_var_3d(ncid, name, dimids, varid, deflate_level, xtype, ierr)
      !! Define a 3D variable (x, y, time) with working precision type
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: name
      integer, intent(in) :: dimids(3)
      integer, intent(out) :: varid
      integer, intent(in), optional :: deflate_level
         !! Compression level (0=none, 1-9=deflate). Enables shuffle when > 0.
      integer, intent(in), optional :: xtype
         !! NetCDF element type.  Absent => `NC_WP` (working precision), so
         !! every legacy call site is unchanged.  The ocean diagnostic
         !! stream passes `NC_R4` when single-precision output is requested;
         !! the caller is then responsible for handing `nf90_put_var` real32
         !! data (see `nc_put_var_3d_slice_r4`).
      integer, intent(out), optional :: ierr
         !! Non-zero on failure when present; absent behaves as today
         !! (`error stop`).
      integer :: xt

      xt = NC_WP
      if (present(xtype)) xt = xtype

      call nc_check(nf90_def_var(ncid, name, xt, dimids, varid), &
                    "defining variable "//trim(name), ierr)
      if (present(ierr)) then
         if (ierr /= nf90_noerr) return
      end if

      if (present(deflate_level)) then
         if (deflate_level > 0) then
            call nc_check(nf90_def_var_deflate(ncid, varid, &
                                               shuffle=1, deflate=1, deflate_level=deflate_level), &
                          "enabling compression for "//trim(name), ierr)
         end if
      end if

   end subroutine nc_def_var_3d

   subroutine nc_put_att(ncid, varid, name, value, ierr)
      !! Write a character attribute to a variable
      integer, intent(in) :: ncid
      integer, intent(in) :: varid
      character(len=*), intent(in) :: name
      character(len=*), intent(in) :: value
      integer, intent(out), optional :: ierr
         !! Non-zero on failure when present; absent behaves as today
         !! (`error stop`).

      call nc_check(nf90_put_att(ncid, varid, name, value), &
                    "writing attribute "//trim(name), ierr)

   end subroutine nc_put_att

   subroutine nc_put_att_real(ncid, varid, name, value, ierr)
      !! Write a real-valued attribute to a variable
      integer, intent(in) :: ncid
      integer, intent(in) :: varid
      character(len=*), intent(in) :: name
      real(wp), intent(in) :: value
      integer, intent(out), optional :: ierr
         !! Non-zero on failure when present; absent behaves as today
         !! (`error stop`).

      call nc_check(nf90_put_att(ncid, varid, name, value), &
                    "writing attribute "//trim(name), ierr)

   end subroutine nc_put_att_real

   subroutine nc_put_att_real_r4(ncid, varid, name, value, ierr)
      !! Write a 32-bit real attribute to a variable.  NetCDF rejects a
      !! `_FillValue` whose type differs from the variable's, so a variable
      !! defined with `NC_R4` must have its `_FillValue` / `missing_value`
      !! written through here, not through `nc_put_att_real` (real64).
      integer, intent(in) :: ncid
      integer, intent(in) :: varid
      character(len=*), intent(in) :: name
      real(real32), intent(in) :: value
      integer, intent(out), optional :: ierr
         !! Non-zero on failure when present; absent behaves as today
         !! (`error stop`).

      call nc_check(nf90_put_att(ncid, varid, name, value), &
                    "writing real32 attribute "//trim(name), ierr)

   end subroutine nc_put_att_real_r4

   subroutine nc_put_att_global(ncid, name, value, ierr)
      !! Write a global character attribute
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: name
      character(len=*), intent(in) :: value
      integer, intent(out), optional :: ierr
         !! Non-zero on failure when present; absent behaves as today
         !! (`error stop`).

      call nc_check(nf90_put_att(ncid, nf90_global, name, value), &
                    "writing global attribute "//trim(name), ierr)

   end subroutine nc_put_att_global

   subroutine nc_put_att_int(ncid, varid, name, value, ierr)
      !! Write an integer attribute to a variable
      integer, intent(in) :: ncid
      integer, intent(in) :: varid
      character(len=*), intent(in) :: name
      integer, intent(in) :: value
      integer, intent(out), optional :: ierr
         !! Non-zero on failure when present; absent behaves as today
         !! (`error stop`).

      call nc_check(nf90_put_att(ncid, varid, name, value), &
                    "writing attribute "//trim(name), ierr)

   end subroutine nc_put_att_int

   subroutine nc_put_att_global_int(ncid, name, value, ierr)
      !! Write a global integer attribute
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: name
      integer, intent(in) :: value
      integer, intent(out), optional :: ierr
         !! Non-zero on failure when present; absent behaves as today
         !! (`error stop`).

      call nc_check(nf90_put_att(ncid, nf90_global, name, value), &
                    "writing global attribute "//trim(name), ierr)

   end subroutine nc_put_att_global_int

   subroutine nc_get_att_int(ncid, name, value, ierr)
      !! Read a global integer attribute
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: name
      integer, intent(out) :: value
      integer, intent(out), optional :: ierr
         !! Non-zero on failure when present; absent behaves as today
         !! (`error stop`).

      call nc_check(nf90_get_att(ncid, nf90_global, name, value), &
                    "reading global attribute "//trim(name), ierr)

   end subroutine nc_get_att_int

   subroutine nc_enddef(ncid, ierr)
      !! End define mode, switch to data mode
      integer, intent(in) :: ncid
      integer, intent(out), optional :: ierr
         !! Non-zero on failure when present; absent behaves as today
         !! (`error stop`).

      call nc_check(nf90_enddef(ncid), "ending define mode", ierr)

   end subroutine nc_enddef

   subroutine nc_put_var_2d(ncid, varid, data, ierr)
      !! Write a full 2D array
      integer, intent(in) :: ncid
      integer, intent(in) :: varid
      real(wp), intent(in) :: data(:, :)
      integer, intent(out), optional :: ierr
         !! Non-zero on failure when present; absent behaves as today
         !! (`error stop`).

      call nc_check(nf90_put_var(ncid, varid, data), "writing 2D variable", ierr)

   end subroutine nc_put_var_2d

   subroutine nc_put_var_3d_slice(ncid, varid, data, time_index, ierr)
      !! Write a 2D slice into a 3D variable at a given time index
      integer, intent(in) :: ncid
      integer, intent(in) :: varid
      real(wp), intent(in) :: data(:, :)
      integer, intent(in) :: time_index
      integer, intent(out), optional :: ierr
         !! Non-zero on failure when present; absent behaves as today
         !! (`error stop`).

      call nc_check(nf90_put_var(ncid, varid, data, &
                                 start=[1, 1, time_index], &
                                 count=[size(data, 1), size(data, 2), 1]), &
                    "writing 3D slice", ierr)

   end subroutine nc_put_var_3d_slice

   subroutine nc_put_var_3d_slice_r4(ncid, varid, data, time_index, ierr)
      !! Write a 32-bit 2D slice into a 3D variable at a given time index.
      !! The real32 twin of `nc_put_var_3d_slice`: handing real64 data to
      !! an `NF90_FLOAT` variable is legal (netcdf converts) but writes the
      !! same bytes at twice the host->library traffic, so the diag writer
      !! stages the conversion itself and calls this.
      integer, intent(in) :: ncid
      integer, intent(in) :: varid
      real(real32), intent(in) :: data(:, :)
      integer, intent(in) :: time_index
      integer, intent(out), optional :: ierr
         !! Non-zero on failure when present; absent behaves as today
         !! (`error stop`).

      call nc_check(nf90_put_var(ncid, varid, data, &
                                 start=[1, 1, time_index], &
                                 count=[size(data, 1), size(data, 2), 1]), &
                    "writing 3D slice (real32)", ierr)

   end subroutine nc_put_var_3d_slice_r4

   subroutine nc_put_var_4d_slice_r4(ncid, varid, data, time_index, ierr)
      !! Write a 32-bit 3D slice into a 4D variable at a given time index.
      !! The real32 twin of `nc_put_var_4d_slice`.
      integer, intent(in) :: ncid
      integer, intent(in) :: varid
      real(real32), intent(in) :: data(:, :, :)
      integer, intent(in) :: time_index
      integer, intent(out), optional :: ierr
         !! Non-zero on failure when present; absent behaves as today
         !! (`error stop`).

      call nc_check(nf90_put_var(ncid, varid, data, &
                                 start=[1, 1, 1, time_index], &
                                 count=[size(data, 1), size(data, 2), size(data, 3), 1]), &
                    "writing 4D slice (real32)", ierr)

   end subroutine nc_put_var_4d_slice_r4

   subroutine nc_get_var_2d(ncid, varid, data, ierr)
      !! Read a full 2D array
      integer, intent(in) :: ncid
      integer, intent(in) :: varid
      real(wp), intent(inout) :: data(:, :)
      integer, intent(out), optional :: ierr
         !! Non-zero on failure when present; absent behaves as today
         !! (`error stop`).

      call nc_check(nf90_get_var(ncid, varid, data), "reading 2D variable", ierr)

   end subroutine nc_get_var_2d

   subroutine nc_get_var_3d(ncid, varid, data, ierr)
      !! Read a full 3D array.  Mirrors `nc_get_var_2d`; the target's
      !! shape must match the file variable's Fortran storage order
      !! (the caller handles any C/Fortran dimension reversal, as
      !! `rdb_bathymetry` does for 2D).
      integer, intent(in) :: ncid
      integer, intent(in) :: varid
      real(wp), intent(inout) :: data(:, :, :)
      integer, intent(out), optional :: ierr
         !! Non-zero on failure when present; absent behaves as today
         !! (`error stop`).

      call nc_check(nf90_get_var(ncid, varid, data), "reading 3D variable", ierr)

   end subroutine nc_get_var_3d

   subroutine nc_get_var_1d(ncid, varid, data, ierr)
      !! Read a full 1D array (the read sibling of `rdb_put_var_1d`).
      integer, intent(in) :: ncid
      integer, intent(in) :: varid
      real(wp), intent(inout) :: data(:)
      integer, intent(out), optional :: ierr
         !! Non-zero on failure when present; absent behaves as today
         !! (`error stop`).

      call nc_check(nf90_get_var(ncid, varid, data), "reading 1D variable", ierr)

   end subroutine nc_get_var_1d

   subroutine nc_get_dim_len(ncid, name, length, ierr)
      !! Get the length of a named dimension
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: name
      integer, intent(out) :: length
      integer, intent(out), optional :: ierr
         !! Non-zero (missing dimension) when present; absent behaves as
         !! today (`error stop`).

      integer :: dimid

      call nc_check(nf90_inq_dimid(ncid, name, dimid), &
                    "finding dimension "//trim(name), ierr)
      if (present(ierr)) then
         if (ierr /= nf90_noerr) return
      end if
      call nc_check(nf90_inquire_dimension(ncid, dimid, len=length), &
                    "querying dimension "//trim(name), ierr)

   end subroutine nc_get_dim_len

   subroutine nc_def_var_4d(ncid, name, dimids, varid, deflate_level, xtype, ierr)
      !! Define a 4D variable (x, y, z, time) with working precision type
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: name
      integer, intent(in) :: dimids(4)
      integer, intent(out) :: varid
      integer, intent(in), optional :: deflate_level
         !! Compression level (0=none, 1-9=deflate). Enables shuffle when > 0.
      integer, intent(in), optional :: xtype
         !! NetCDF element type.  Absent => `NC_WP` (working precision).
         !! See `nc_def_var_3d` — same contract, 4D shape.
      integer, intent(out), optional :: ierr
         !! Non-zero on failure when present; absent behaves as today
         !! (`error stop`).
      integer :: xt

      xt = NC_WP
      if (present(xtype)) xt = xtype

      call nc_check(nf90_def_var(ncid, name, xt, dimids, varid), &
                    "defining variable "//trim(name), ierr)
      if (present(ierr)) then
         if (ierr /= nf90_noerr) return
      end if

      if (present(deflate_level)) then
         if (deflate_level > 0) then
            call nc_check(nf90_def_var_deflate(ncid, varid, &
                                               shuffle=1, deflate=1, deflate_level=deflate_level), &
                          "enabling compression for "//trim(name), ierr)
         end if
      end if

   end subroutine nc_def_var_4d

   subroutine nc_put_var_4d_slice(ncid, varid, data, time_index, ierr)
      !! Write a 3D slice into a 4D variable at a given time index
      integer, intent(in) :: ncid
      integer, intent(in) :: varid
      real(wp), intent(in) :: data(:, :, :)
      integer, intent(in) :: time_index
      integer, intent(out), optional :: ierr
         !! Non-zero on failure when present; absent behaves as today
         !! (`error stop`).

      call nc_check(nf90_put_var(ncid, varid, data, &
                                 start=[1, 1, 1, time_index], &
                                 count=[size(data, 1), size(data, 2), size(data, 3), 1]), &
                    "writing 4D slice", ierr)

   end subroutine nc_put_var_4d_slice

   subroutine nc_get_var_slab_3d(ncid, varid, start, count, data, ierr)
      !! Strided slab read: `start`/`count` are the NetCDF start/count
      !! vectors in FORTRAN dimension order, with length matching the
      !! target variable's rank (3 for a 2D-field-plus-time variable, 4
      !! for a 3D-field-plus-time variable).  `data` is always a plain
      !! 3D Fortran array; for the 2D-plus-time case pass a target with
      !! a trailing singleton extent (`data(:,:,1)`) — the netcdf-fortran
      !! generic interface dispatches on the array RANK of `data`, not
      !! the length of `start`/`count`, so one routine covers both
      !! shapes. This is the first strided (non-full-variable) read in
      !! the tree — mirrors `nc_put_var_4d_slice` for the start/count
      !! shape, on the read side.
      integer, intent(in) :: ncid, varid
      integer, intent(in) :: start(:), count(:)
      real(wp), intent(inout) :: data(:, :, :)
      integer, intent(out), optional :: ierr
         !! Non-zero on failure when present; absent behaves as today
         !! (`error stop`).

      call nc_check(nf90_get_var(ncid, varid, data, start=start, count=count), &
                    "reading NetCDF slab", ierr)

   end subroutine nc_get_var_slab_3d

   subroutine nc_get_att_text(ncid, varid, name, value, ok)
      !! Read a character (text) attribute from a variable — used for the
      !! CF `units` attribute on a time axis.  Non-fail-loud: `ok` is
      !! `.false.` (and `value` blank) when the attribute is absent or the
      !! wrong type, so callers can fall back to a default rather than
      !! aborting on a missing/optional attribute.
      integer, intent(in) :: ncid, varid
      character(len=*), intent(in) :: name
      character(len=:), allocatable, intent(out) :: value
         !! Deferred-length (not assumed-length `intent(out)`, which risks
         !! silent truncation — fortitude C072): filled from a fixed local
         !! buffer sized for a CF `units` string.
      logical, intent(out) :: ok

      integer :: status
      character(len=1024) :: buf

      buf = ""
      status = nf90_get_att(ncid, varid, name, buf)
      ok = (status == nf90_noerr)
      if (ok) then
         value = trim(buf)
      else
         value = ""
      end if

   end subroutine nc_get_att_text

   subroutine nc_get_varid(ncid, name, varid, ierr)
      !! Get the variable ID for a named variable
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: name
      integer, intent(out) :: varid
      integer, intent(out), optional :: ierr
         !! Non-zero (no matching variable) when present; absent behaves
         !! as today (`error stop`).

      call nc_check(nf90_inq_varid(ncid, name, varid), &
                    "finding variable "//trim(name), ierr)

   end subroutine nc_get_varid

   function output_rank_filename(base_dir, prefix, rank) result(fname)
      !! Generate a per-rank filename: <base_dir>/<prefix>_rank_NNNNNN.nc
      character(len=*), intent(in) :: base_dir
         !! Output directory
      character(len=*), intent(in) :: prefix
         !! File prefix (e.g. "output" or "restart")
      integer, intent(in) :: rank
         !! MPI rank number
      character(len=512) :: fname

      write (fname, '(A,"/",A,"_rank_",I6.6,".nc")') &
         trim(base_dir), trim(prefix), rank

   end function output_rank_filename

   subroutine ensure_directory_exists(path, ierr)
      !! Create `path` via libc `mkdir(2)` if it does not already exist
      !! (single path component — not a recursive `mkdir -p`; a missing
      !! grandparent still fails, loudly, at the `inquire` re-check
      !! below). Exists to close the headline P7 crash: `&ocean_diag_nml`/
      !! `&output_nml` default `output_dir = "./output"`
      !! (`rdb_config.F90:2316`), so a first-time `create()` against a
      !! fresh checkout with no `./output/` directory would otherwise
      !! reach `nc_create_file` -> `nf90_create` -> ENOENT ->
      !! `error stop`, killing the whole host process.
      !!
      !! Never `error stop`s itself, present `ierr` or not — the caller
      !! decides how to fail (via its OWN `nc_create_file(..., ierr=)`
      !! call on the resulting still-missing directory, or by checking
      !! this routine's `ierr` directly). `ierr` is a plain 0-ok/
      !! non-zero-failed flag (not an `nf90_*` status), because `mkdir`
      !! is not a NetCDF call and has no `nc_check` counterpart.
      character(len=*), intent(in) :: path
      integer, intent(out), optional :: ierr
         !! Non-zero iff `path` still does not exist after the `mkdir`
         !! attempt. Absent is legal — the routine just does its best and
         !! lets the caller's own I/O call fail loud.

      logical :: exists
      integer(c_int) :: rc
      character(kind=c_char), allocatable, target :: c_path(:)
      interface
         function c_mkdir(dirpath, mode) bind(C, name="mkdir") result(r)
            import :: c_ptr, c_int
            implicit none
            type(c_ptr), value, intent(in) :: dirpath
            integer(c_int), value, intent(in) :: mode
            integer(c_int) :: r
         end function c_mkdir
      end interface

      if (present(ierr)) ierr = 0
      if (len_trim(path) == 0) return

      inquire (file=trim(path), exist=exists)
      if (exists) return

      c_path = path_to_c_string(trim(path))
      rc = c_mkdir(c_loc(c_path), int(o'755', c_int))
      if (rc /= 0_c_int) then
         call logger%warning("ensure_directory_exists: mkdir('"//trim(path)// &
                             "') returned "//to_string(int(rc))// &
                             " -- the subsequent NetCDF create/open will "// &
                             "report the exact failure if the directory "// &
                             "is still missing")
      end if

      inquire (file=trim(path), exist=exists)
      if (.not. exists) then
         if (present(ierr)) ierr = -1
      end if

   end subroutine ensure_directory_exists

   pure function path_to_c_string(s) result(c)
      !! Pack a Fortran string into a null-terminated c_char array for
      !! the `mkdir` binding above. Private duplicate of
      !! `rdb_ocean_restart_io`'s `string_to_c` (that one is private to
      !! its own module) — small enough that sharing it is not worth a
      !! cross-module dependency for a two-line helper.
      character(len=*), intent(in) :: s
      character(kind=c_char), allocatable :: c(:)
      integer :: i, n
      n = len_trim(s)
      allocate (c(n + 1))
      do i = 1, n
         c(i) = s(i:i)
      end do
      c(n + 1) = c_null_char
   end function path_to_c_string

end module rdb_io_netcdf
