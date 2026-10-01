!! Serial NetCDF emit hook for the ocean diagnostics registry.
module rdb_ocean_diag_netcdf
   !! Writes the diag manager's per-fire output to a per-rank NetCDF
   !! file. Serial only (one file per process). Each diag var carries its
   !! own `time_<name>` unlimited dim + coord var, so vars firing at
   !! different cadences stay independent. CF-1.8 attributes.
   !!
   !! Usage:
   !!     call open_stream(state%diag, "rdb_diag.nc")
   !!     ! ... step the diag manager as usual ...
   !!     call close_stream(state%diag)
   !!
   !! After `open_stream`, the diag manager emits a slice on every
   !! cadence-fire via the `emit_post_fire` pointer bound here.
   use, intrinsic :: iso_fortran_env, only: int32, real32
   use netcdf, only: nf90_put_var, nf90_unlimited, nf90_global, nf90_sync
   use rdb_constants, only: wp
   use rdb_ocean_diag, only: ocean_diag_t, diag_var_t, DIAG_MISSING_VALUE
   use rdb_io_netcdf, only: nc_check, nc_create_file, nc_close, &
                            nc_def_dim, nc_def_var_3d, nc_def_var_4d, &
                            nc_enddef, nc_put_att, nc_put_att_real, &
                            nc_put_att_real_r4, &
                            nc_put_att_global, nc_put_att_global_int, &
                            nc_put_var_3d_slice, nc_put_var_4d_slice, &
                            nc_put_var_3d_slice_r4, nc_put_var_4d_slice_r4, &
                            rdb_def_var_1d, NC_WP, NC_R4
   use pic_logger, only: logger => global_logger
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_IO
   implicit none
   private

   public :: open_stream, close_stream
   public :: diag_xtype_from_name

contains

   function diag_xtype_from_name(name) result(xtype)
      !! Map the `&ocean_diag_nml output_precision` string onto a NetCDF
      !! element type for the diagnostic DATA variables.
      !!
      !!   * `"double"` (default, or anything empty) => `NC_WP` — the
      !!     working-precision type, i.e. byte-identical to the writer
      !!     before this knob existed.
      !!   * `"single"` => `NC_R4` — halves the bytes written per frame.
      !!
      !! The value is validated by the namelist schema, so an unknown
      !! string should never reach here; it is warned about and falls back
      !! to `NC_WP` (never silently downgrades precision).
      !!
      !! This governs the diagnostic stream ONLY.  Restart files, gauges,
      !! coastal output and the console conservation totals do not consult
      !! it and cannot be made single precision.
      character(len=*), intent(in) :: name
      integer :: xtype

      select case (trim(name))
      case ("", "double")
         xtype = NC_WP
      case ("single")
         xtype = NC_R4
      case default
         call logger%warning("&ocean_diag_nml output_precision = '"//trim(name)// &
                             "' not recognised; writing double precision")
         xtype = NC_WP
      end select
   end function diag_xtype_from_name

   subroutine open_stream(diag, filename, deflate_level, &
                          px, py, i_start, j_start, &
                          nx_global, ny_global, nx_local, ny_local, nghost, &
                          output_precision, ierr)
      !! Create the NetCDF file, define dims + vars for every
      !! currently-registered diag var, and bind the post-fire emit hook.
      !! Must be called AFTER all `register` calls land — vars registered
      !! later get no NetCDF binding and are not written.
      !!
      !! `deflate_level` (0-9, default 0 = uncompressed, bit-identical)
      !! enables gzip compression on data vars; level 1 is usually
      !! optimal (3-5× compression, near-zero CPU).
      !!
      !! The optional decomposition attrs (`px`, `py`, `i_start`, `j_start`,
      !! `nx_global`, `ny_global`, `nx_local`, `ny_local`, `nghost`) are
      !! written as global integer attributes when present, making the file
      !! mergeable by `tools/merge_output.py`.  Absent => no attrs written;
      !! the file is bit-identical to the pre-MPI single-rank output.
      !!
      !! `ierr` (P7 F1 fix): every internal `nc_*` call used to `error
      !! stop` unconditionally — a first-run `create()` against a fresh
      !! checkout with no `./output/` directory killed the whole host
      !! process from `nc_create_file`. When `ierr` is present, any
      !! failure (from `nc_create_file` itself through to the last
      !! attribute write) closes the partially-created file, logs +
      !! rings the specific NetCDF reason (via `nc_check`/`fail`), sets
      !! `ierr = OCEAN_STATUS_ERR_IO`, and returns with `is_open` left
      !! `.false.` — never `error stop`s. Absent `ierr` preserves the
      !! legacy abort behaviour byte-for-byte.
      class(ocean_diag_t), intent(inout) :: diag
      character(len=*), intent(in) :: filename
      integer, intent(in), optional :: deflate_level
      integer, intent(in), optional :: px
         !! Number of ranks in the x direction (process grid).
      integer, intent(in), optional :: py
         !! Number of ranks in the y direction (process grid).
      integer, intent(in), optional :: i_start
         !! 1-based global column index of this rank's first physical cell.
      integer, intent(in), optional :: j_start
         !! 1-based global row index of this rank's first physical cell.
      integer, intent(in), optional :: nx_global
         !! Total number of physical columns in the global domain.
      integer, intent(in), optional :: ny_global
         !! Total number of physical rows in the global domain.
      integer, intent(in), optional :: nx_local
         !! Number of physical columns on this rank (excluding ghost cells).
      integer, intent(in), optional :: ny_local
         !! Number of physical rows on this rank (excluding ghost cells).
      integer, intent(in), optional :: nghost
         !! Number of ghost cells on each side of the local domain.
      character(len=*), intent(in), optional :: output_precision
         !! `"double"` (default / absent) or `"single"` — element type of
         !! the DATA variables.  Absent => `NC_WP` => byte-identical to the
         !! pre-knob writer.  Time coordinate variables stay `NC_WP`
         !! regardless: they are a handful of values per frame and are
         !! routinely differenced.
      integer, intent(out), optional :: ierr
         !! Non-zero (`OCEAN_STATUS_ERR_IO`) on any NetCDF failure when
         !! present; absent behaves as today (`error stop`).
      integer :: i, n1, n2, n3, dl, n3_max
      integer :: local_ierr

      character(len=80) :: time_dim_name, z_dim_name

      if (present(ierr)) ierr = OCEAN_STATUS_OK

      dl = 0
      if (present(deflate_level)) dl = max(0, min(9, deflate_level))

      if (diag%nc_stream%is_open) return
      diag%nc_stream%filename = filename

      diag%nc_stream%xtype = NC_WP
      if (present(output_precision)) then
         diag%nc_stream%xtype = diag_xtype_from_name(output_precision)
      end if

      call nc_create_file(filename, diag%nc_stream%ncid, ierr=local_ierr)
      if (.not. diag_io_ok(local_ierr, ierr)) return
      call nc_put_att_global(diag%nc_stream%ncid, "Conventions", "CF-1.8", ierr=local_ierr)
      if (.not. diag_io_ok(local_ierr, ierr, diag%nc_stream%ncid)) return
      call nc_put_att_global(diag%nc_stream%ncid, "title", &
                             "Roundabout ocean-path diagnostics", ierr=local_ierr)
      if (.not. diag_io_ok(local_ierr, ierr, diag%nc_stream%ncid)) return

      ! Decomposition window attrs — written only in multi-rank runs so
      ! single-rank files stay bit-identical.  The merge tool uses these to
      ! place each rank's physical-cell slice into the global array, trimming
      ! the `nghost` halo on each side before stitching.
      if (present(px)) then
         call nc_put_att_global_int(diag%nc_stream%ncid, "px", px, ierr=local_ierr)
         if (.not. diag_io_ok(local_ierr, ierr, diag%nc_stream%ncid)) return
      end if
      if (present(py)) then
         call nc_put_att_global_int(diag%nc_stream%ncid, "py", py, ierr=local_ierr)
         if (.not. diag_io_ok(local_ierr, ierr, diag%nc_stream%ncid)) return
      end if
      if (present(i_start)) then
         call nc_put_att_global_int(diag%nc_stream%ncid, "i_start", i_start, ierr=local_ierr)
         if (.not. diag_io_ok(local_ierr, ierr, diag%nc_stream%ncid)) return
      end if
      if (present(j_start)) then
         call nc_put_att_global_int(diag%nc_stream%ncid, "j_start", j_start, ierr=local_ierr)
         if (.not. diag_io_ok(local_ierr, ierr, diag%nc_stream%ncid)) return
      end if
      if (present(nx_global)) then
         call nc_put_att_global_int(diag%nc_stream%ncid, "nx_global", nx_global, ierr=local_ierr)
         if (.not. diag_io_ok(local_ierr, ierr, diag%nc_stream%ncid)) return
      end if
      if (present(ny_global)) then
         call nc_put_att_global_int(diag%nc_stream%ncid, "ny_global", ny_global, ierr=local_ierr)
         if (.not. diag_io_ok(local_ierr, ierr, diag%nc_stream%ncid)) return
      end if
      if (present(nx_local)) then
         call nc_put_att_global_int(diag%nc_stream%ncid, "nx_local", nx_local, ierr=local_ierr)
         if (.not. diag_io_ok(local_ierr, ierr, diag%nc_stream%ncid)) return
      end if
      if (present(ny_local)) then
         call nc_put_att_global_int(diag%nc_stream%ncid, "ny_local", ny_local, ierr=local_ierr)
         if (.not. diag_io_ok(local_ierr, ierr, diag%nc_stream%ncid)) return
      end if
      if (present(nghost)) then
         call nc_put_att_global_int(diag%nc_stream%ncid, "nghost", nghost, ierr=local_ierr)
         if (.not. diag_io_ok(local_ierr, ierr, diag%nc_stream%ncid)) return
      end if

      ! Spatial dims from var 1's output buffer (all vars share the
      ! same horizontal extent).
      if (diag%nvars > 0) then
         n1 = size(diag%vars(1)%output_buffer, 1)
         n2 = size(diag%vars(1)%output_buffer, 2)
         call nc_def_dim(diag%nc_stream%ncid, "x", n1, diag%nc_stream%x_dimid, ierr=local_ierr)
         if (.not. diag_io_ok(local_ierr, ierr, diag%nc_stream%ncid)) return
         call nc_def_dim(diag%nc_stream%ncid, "y", n2, diag%nc_stream%y_dimid, ierr=local_ierr)
         if (.not. diag_io_ok(local_ierr, ierr, diag%nc_stream%ncid)) return
      end if

      do i = 1, diag%nvars
         associate (v => diag%vars(i))
            n3 = size(v%output_buffer, 3)

            ! Per-var time dim + coord var (unlimited).
            write (time_dim_name, "(A,A)") "time_", trim(v%name)
            call nc_def_dim(diag%nc_stream%ncid, trim(time_dim_name), &
                            nf90_unlimited, v%nc_time_dimid, ierr=local_ierr)
            if (.not. diag_io_ok(local_ierr, ierr, diag%nc_stream%ncid)) return
            call rdb_def_var_1d(diag%nc_stream%ncid, trim(time_dim_name), &
                                v%nc_time_dimid, v%nc_time_varid, ierr=local_ierr)
            if (.not. diag_io_ok(local_ierr, ierr, diag%nc_stream%ncid)) return
            call nc_put_att(diag%nc_stream%ncid, v%nc_time_varid, &
                            "units", "seconds since simulation start", ierr=local_ierr)
            if (.not. diag_io_ok(local_ierr, ierr, diag%nc_stream%ncid)) return
            call nc_put_att(diag%nc_stream%ncid, v%nc_time_varid, &
                            "long_name", "time", ierr=local_ierr)
            if (.not. diag_io_ok(local_ierr, ierr, diag%nc_stream%ncid)) return

            ! Data var: (x, y, time) for 2D; (x, y, z_<name>, time) for 3D.
            if (n3 == 1) then
               call nc_def_var_3d(diag%nc_stream%ncid, trim(v%name), &
                                  [diag%nc_stream%x_dimid, diag%nc_stream%y_dimid, &
                                   v%nc_time_dimid], v%nc_varid, &
                                  deflate_level=dl, xtype=diag%nc_stream%xtype, ierr=local_ierr)
               if (.not. diag_io_ok(local_ierr, ierr, diag%nc_stream%ncid)) return
            else
               write (z_dim_name, "(A,A)") "z_", trim(v%name)
               call nc_def_dim(diag%nc_stream%ncid, trim(z_dim_name), &
                               n3, v%nc_z_dimid, ierr=local_ierr)
               if (.not. diag_io_ok(local_ierr, ierr, diag%nc_stream%ncid)) return
               call nc_def_var_4d(diag%nc_stream%ncid, trim(v%name), &
                                  [diag%nc_stream%x_dimid, diag%nc_stream%y_dimid, &
                                   v%nc_z_dimid, v%nc_time_dimid], v%nc_varid, &
                                  deflate_level=dl, xtype=diag%nc_stream%xtype, ierr=local_ierr)
               if (.not. diag_io_ok(local_ierr, ierr, diag%nc_stream%ncid)) return
            end if

            if (len_trim(v%units) > 0) then
               call nc_put_att(diag%nc_stream%ncid, v%nc_varid, "units", trim(v%units), ierr=local_ierr)
               if (.not. diag_io_ok(local_ierr, ierr, diag%nc_stream%ncid)) return
            end if
            if (len_trim(v%long_name) > 0) then
               call nc_put_att(diag%nc_stream%ncid, v%nc_varid, "long_name", trim(v%long_name), ierr=local_ierr)
               if (.not. diag_io_ok(local_ierr, ierr, diag%nc_stream%ncid)) return
            end if
            if (len_trim(v%standard_name) > 0) then
               call nc_put_att(diag%nc_stream%ncid, v%nc_varid, &
                               "standard_name", trim(v%standard_name), ierr=local_ierr)
               if (.not. diag_io_ok(local_ierr, ierr, diag%nc_stream%ncid)) return
            end if
            if (v%has_missing) then
               ! Vanished (no-water) target cells carry the sentinel; advertise
               ! it both ways so CF-aware (_FillValue) and legacy
               ! (missing_value) tools mask below-bottom / dry cells.
               !
               ! NetCDF REJECTS a `_FillValue` whose type differs from the
               ! variable's, so the attribute has to follow `xtype`.  The
               ! real32 branch writes `real(DIAG_MISSING_VALUE, real32)` —
               ! the exact same value the data conversion produces for a
               ! sentinel cell, so masking still matches bit-for-bit.
               if (diag%nc_stream%xtype == NC_R4) then
                  call nc_put_att_real_r4(diag%nc_stream%ncid, v%nc_varid, &
                                          "_FillValue", real(DIAG_MISSING_VALUE, real32), ierr=local_ierr)
                  if (.not. diag_io_ok(local_ierr, ierr, diag%nc_stream%ncid)) return
                  call nc_put_att_real_r4(diag%nc_stream%ncid, v%nc_varid, &
                                          "missing_value", real(DIAG_MISSING_VALUE, real32), ierr=local_ierr)
                  if (.not. diag_io_ok(local_ierr, ierr, diag%nc_stream%ncid)) return
               else
                  call nc_put_att_real(diag%nc_stream%ncid, v%nc_varid, &
                                       "_FillValue", DIAG_MISSING_VALUE, ierr=local_ierr)
                  if (.not. diag_io_ok(local_ierr, ierr, diag%nc_stream%ncid)) return
                  call nc_put_att_real(diag%nc_stream%ncid, v%nc_varid, &
                                       "missing_value", DIAG_MISSING_VALUE, ierr=local_ierr)
                  if (.not. diag_io_ok(local_ierr, ierr, diag%nc_stream%ncid)) return
               end if
            end if

            v%nc_time_index = 0
         end associate
      end do

      call nc_enddef(diag%nc_stream%ncid, ierr=local_ierr)
      if (.not. diag_io_ok(local_ierr, ierr, diag%nc_stream%ncid)) return

      ! Single-precision staging buffer.  One shared scratch covering the
      ! largest registered `output_buffer`: every var shares (n1, n2), so
      ! only the vertical extent varies and `n3_max` bounds them all.  This
      ! is why `open_stream` must run after every `register` — the same
      ! precondition the per-var NetCDF binding above already imposes.
      if (diag%nc_stream%xtype == NC_R4 .and. diag%nvars > 0) then
         n3_max = 1
         do i = 1, diag%nvars
            n3_max = max(n3_max, size(diag%vars(i)%output_buffer, 3))
         end do
         if (allocated(diag%nc_stream%stage)) deallocate (diag%nc_stream%stage)
         allocate (diag%nc_stream%stage(n1, n2, n3_max))
      end if

      diag%nc_stream%is_open = .true.
      diag%emit_post_fire => nc_emit_post_fire
   end subroutine open_stream

   function diag_io_ok(local_ierr, ierr, ncid) result(ok)
      !! Translate a raw `nc_check`-style status (0 = ok) from one of the
      !! `nc_*` calls in `open_stream` into the caller's `ierr` contract:
      !! `.true.` on success; on failure, `.false.` with
      !! `ierr = OCEAN_STATUS_ERR_IO` when `ierr` is present (closing
      !! `ncid` first, when given, so a mid-create failure does not leak
      !! the file handle), or `error stop`s with the same generic text
      !! `nc_check` itself would use had `ierr` never been threaded
      !! through — byte-identical legacy behaviour when `ierr` is
      !! omitted. Mirrors `rdb_bathymetry`'s `bathy_io_ok` (P0.1 F1).
      integer, intent(in) :: local_ierr
      integer, intent(out), optional :: ierr
      integer, intent(in), optional :: ncid
      logical :: ok

      integer :: discard_ierr

      ok = (local_ierr == 0)
      if (ok) return

      if (present(ierr)) then
         if (present(ncid)) call nc_close(ncid, ierr=discard_ierr)
         ierr = OCEAN_STATUS_ERR_IO
         return
      end if

      error stop "NetCDF operation failed"
   end function diag_io_ok

   subroutine close_stream(diag)
      !! Flush + close the NetCDF file.  Idempotent.
      class(ocean_diag_t), intent(inout) :: diag
      if (.not. diag%nc_stream%is_open) return
      call nc_close(diag%nc_stream%ncid)
      diag%nc_stream%is_open = .false.
      diag%nc_stream%ncid = -1
      if (allocated(diag%nc_stream%stage)) deallocate (diag%nc_stream%stage)
      nullify (diag%emit_post_fire)
   end subroutine close_stream

   subroutine nc_emit_post_fire(diag, ivar, t)
      !! Post-fire emit hook.  Appends one slice + time value for the
      !! var that just fired.  Bound to `diag%emit_post_fire` by
      !! `open_stream`.
      class(*), intent(inout) :: diag
      integer, intent(in) :: ivar
      real(wp), intent(in) :: t
      integer :: n1, n2, n3, ti

      select type (d => diag)
      class is (ocean_diag_t)
         if (.not. d%nc_stream%is_open) return
         if (ivar < 1 .or. ivar > d%nvars) return
         associate (v => d%vars(ivar))
            if (v%nc_varid < 0) return
            v%nc_time_index = v%nc_time_index + 1
            ti = v%nc_time_index

            call nc_check(nf90_put_var(d%nc_stream%ncid, v%nc_time_varid, &
                                       [t], start=[ti], count=[1]), &
                          "writing diag time coord for "//trim(v%name))

            n3 = size(v%output_buffer, 3)
            if (d%nc_stream%xtype == NC_R4 .and. allocated(d%nc_stream%stage)) then
               ! Convert fp64 -> fp32 in the shared host staging buffer,
               ! then hand real32 data to a real32 variable.  Handing the
               ! fp64 buffer straight to an NF90_FLOAT variable would also
               ! "work" (netcdf converts internally) but doubles the bytes
               ! crossing into the library — the whole point of the knob is
               ! the byte count, so the conversion is ours to do.
               n1 = size(v%output_buffer, 1)
               n2 = size(v%output_buffer, 2)
               d%nc_stream%stage(1:n1, 1:n2, 1:n3) = &
                  real(v%output_buffer(1:n1, 1:n2, 1:n3), real32)
               if (n3 == 1) then
                  call nc_put_var_3d_slice_r4(d%nc_stream%ncid, v%nc_varid, &
                                              d%nc_stream%stage(1:n1, 1:n2, 1), ti)
               else
                  call nc_put_var_4d_slice_r4(d%nc_stream%ncid, v%nc_varid, &
                                              d%nc_stream%stage(1:n1, 1:n2, 1:n3), ti)
               end if
            else if (n3 == 1) then
               call nc_put_var_3d_slice(d%nc_stream%ncid, v%nc_varid, &
                                        v%output_buffer(:, :, 1), ti)
            else
               call nc_put_var_4d_slice(d%nc_stream%ncid, v%nc_varid, &
                                        v%output_buffer, ti)
            end if

            ! Flush the HDF5 metadata cache so frames survive an abort
            ! (e.g. the h-guard ERROR STOP): without this, a killed run
            ! leaves every record dim at 0 and the written data orphaned.
            call nc_check(nf90_sync(d%nc_stream%ncid), &
                          "syncing diag stream after "//trim(v%name))
         end associate
      end select
   end subroutine nc_emit_post_fire

end module rdb_ocean_diag_netcdf
