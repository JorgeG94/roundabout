!! NetCDF read/write of the ocean restart registry.
module rdb_ocean_restart_io
   !! Walks a `restart_registry_t` and writes/reads each registered field
   !! as a FULL local array (interior + ghosts) to/from a per-rank NetCDF
   !! file, alongside the scalar checkpoint metadata (time, step
   !! counters), the decomposition attributes that gate a same-decomp
   !! resume, and the grid/vcoord/tracer metadata that gates a same-schema
   !! resume.
   !!
   !! Durability: the write goes to `<filename>.tmp`, then a POSIX
   !! `rename(2)` atomically swaps it into place — a crash mid-write
   !! leaves the previous checkpoint intact.
   use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char, c_ptr, c_loc
   use rdb_constants, only: wp
   use rdb_decomp, only: decomp_t
   use rdb_ocean_restart, only: restart_registry_t, RESTART_SCHEMA_VERSION
   use rdb_io_netcdf, only: nc_check, nc_create_file, nc_open_read, nc_close, &
                            nc_def_dim, nc_enddef, nc_put_att_global, &
                            nc_put_att_global_int, nc_get_att_int, &
                            NC_WP
   use netcdf, only: nf90_def_var, nf90_put_var, nf90_get_var, nf90_inq_varid, &
                     nf90_double, nf90_int, nf90_noerr
   use pic_logger, only: logger => global_logger
   use pic_strings, only: to_string
   use rdb_error_ring, only: error_ring_push
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_IO
   implicit none
   private

   public :: ocean_restart_write_local
   public :: ocean_restart_read_local
   public :: ocean_restart_check_decomp
   public :: ocean_restart_metadata_t

   type :: ocean_restart_metadata_t
      !! Grid/coordinate/tracer fingerprint validated on read.  Distinct
      !! from the decomposition metadata (px/py/dims) so a file built with
      !! a different vertical resolution, ghost width, vcoord, or tracer
      !! set is rejected loudly rather than silently mis-mapped.
      integer :: nz_ml = 0
         !! Multilayer layer count.
      integer :: nghost = 0
         !! Ghost width.
      integer :: vcoord_type = 0
         !! Vertical-coordinate enum (VCOORD_*).
      character(len=32) :: vcoord_name = ""
         !! Human-readable vcoord tag (for the mismatch message).
      integer :: n_tracers = 0
         !! Registered tracer count.
      integer :: idx_salinity = 0
         !! tracers(:) index of salinity (0 = none).
      integer :: idx_temperature = 0
         !! tracers(:) index of temperature (0 = none).
   end type ocean_restart_metadata_t

contains

   subroutine ocean_restart_write_local(filename, reg, decomp, meta, t, step, &
                                        outer_step_count, ierr)
      !! Write a per-rank ocean restart file from the registry durably
      !! (to `<filename>.tmp` then POSIX-rename).  Device sync
      !! (`!$acc update self`) on the device-mapped arrays must already
      !! have been done by the caller — this is pure host NetCDF I/O.
      !!
      !! `ierr` (P7, same treatment as `rdb_ocean_diag_netcdf`'s
      !! `open_stream`): every internal `nc_*` call, plus the final
      !! atomic rename, used to `error stop` unconditionally. When
      !! `ierr` is present, any failure closes/discards the partial
      !! `.tmp` file, sets `ierr = OCEAN_STATUS_ERR_IO`, and returns
      !! rather than aborting. Absent `ierr` preserves the legacy abort
      !! behaviour byte-for-byte.
      character(len=*), intent(in) :: filename
      type(restart_registry_t), intent(in) :: reg
      type(decomp_t), intent(in) :: decomp
      type(ocean_restart_metadata_t), intent(in) :: meta
      real(wp), intent(in) :: t
      integer, intent(in) :: step
      integer, intent(in) :: outer_step_count
      integer, intent(out), optional :: ierr
         !! Non-zero (`OCEAN_STATUS_ERR_IO`) on any NetCDF or rename
         !! failure when present; absent behaves as today (`error stop`).

      character(len=:), allocatable :: tmpname
      integer :: ncid, e
      integer :: time_varid, step_varid, osc_varid
      integer, allocatable :: varids(:)
      character(len=32) :: dimname_x, dimname_y, dimname_z
      integer :: dx_id, dy_id, dz_id
      integer :: dimids2(2), dimids3(3)
      integer :: rc, local_ierr
      character(len=:), allocatable :: msg

      if (present(ierr)) ierr = OCEAN_STATUS_OK

      tmpname = trim(filename)//".tmp"
      call logger%info("Writing ocean restart file: "//trim(filename)// &
                       " ("//to_string(reg%n)//" fields)")

      call nc_create_file(tmpname, ncid, ierr=local_ierr)
      if (.not. restart_io_ok(local_ierr, ierr)) return
      allocate (varids(reg%n))

      ! One x/y/z dim trio per entry index.  The registry mixes layer
      ! (nk=nz) and interface (nk=nz+1) fields; per-entry dims are simple
      ! and the file is tiny so the duplication costs nothing.
      do e = 1, reg%n
         if (reg%entries(e)%rank == 0) cycle  ! host scalars: no dims
         write (dimname_x, "(A,I0)") "x", e
         write (dimname_y, "(A,I0)") "y", e
         call nc_def_dim(ncid, trim(dimname_x), reg%entries(e)%nx_phys, dx_id, ierr=local_ierr)
         if (.not. restart_io_ok(local_ierr, ierr, ncid)) return
         call nc_def_dim(ncid, trim(dimname_y), reg%entries(e)%ny_phys, dy_id, ierr=local_ierr)
         if (.not. restart_io_ok(local_ierr, ierr, ncid)) return
         if (reg%entries(e)%rank == 2) then
            dimids2 = [dx_id, dy_id]
            call nc_check(nf90_def_var(ncid, trim(reg%entries(e)%tag), NC_WP, &
                                       dimids2, varids(e)), &
                          "defining restart var "//trim(reg%entries(e)%tag), local_ierr)
            if (.not. restart_io_ok(local_ierr, ierr, ncid)) return
         else
            write (dimname_z, "(A,I0)") "z", e
            call nc_def_dim(ncid, trim(dimname_z), reg%entries(e)%nk, dz_id, ierr=local_ierr)
            if (.not. restart_io_ok(local_ierr, ierr, ncid)) return
            dimids3 = [dx_id, dy_id, dz_id]
            call nc_check(nf90_def_var(ncid, trim(reg%entries(e)%tag), NC_WP, &
                                       dimids3, varids(e)), &
                          "defining restart var "//trim(reg%entries(e)%tag), local_ierr)
            if (.not. restart_io_ok(local_ierr, ierr, ncid)) return
         end if
      end do

      ! Host scalars (rank-0 persistent state) as NetCDF scalar vars.
      do e = 1, reg%n
         if (reg%entries(e)%rank /= 0) cycle
         call nc_check(nf90_def_var(ncid, trim(reg%entries(e)%tag), nf90_double, &
                                    varids(e)), &
                       "defining restart scalar "//trim(reg%entries(e)%tag), local_ierr)
         if (.not. restart_io_ok(local_ierr, ierr, ncid)) return
      end do

      ! Scalar metadata (rank-0 globals).
      call nc_check(nf90_def_var(ncid, "time", nf90_double, time_varid), &
                    "defining restart time", local_ierr)
      if (.not. restart_io_ok(local_ierr, ierr, ncid)) return
      call nc_check(nf90_def_var(ncid, "n_steps", nf90_int, step_varid), &
                    "defining restart n_steps", local_ierr)
      if (.not. restart_io_ok(local_ierr, ierr, ncid)) return
      call nc_check(nf90_def_var(ncid, "outer_step_count", nf90_int, osc_varid), &
                    "defining restart outer_step_count", local_ierr)
      if (.not. restart_io_ok(local_ierr, ierr, ncid)) return

      ! Decomposition + schema metadata (validated on read).
      call nc_put_att_global_int(ncid, "schema_version", RESTART_SCHEMA_VERSION, ierr=local_ierr)
      if (.not. restart_io_ok(local_ierr, ierr, ncid)) return
      call nc_put_att_global_int(ncid, "px", decomp%px, ierr=local_ierr)
      if (.not. restart_io_ok(local_ierr, ierr, ncid)) return
      call nc_put_att_global_int(ncid, "py", decomp%py, ierr=local_ierr)
      if (.not. restart_io_ok(local_ierr, ierr, ncid)) return
      call nc_put_att_global_int(ncid, "i_start", decomp%i_start, ierr=local_ierr)
      if (.not. restart_io_ok(local_ierr, ierr, ncid)) return
      call nc_put_att_global_int(ncid, "j_start", decomp%j_start, ierr=local_ierr)
      if (.not. restart_io_ok(local_ierr, ierr, ncid)) return
      call nc_put_att_global_int(ncid, "nx_global", decomp%nx_global, ierr=local_ierr)
      if (.not. restart_io_ok(local_ierr, ierr, ncid)) return
      call nc_put_att_global_int(ncid, "ny_global", decomp%ny_global, ierr=local_ierr)
      if (.not. restart_io_ok(local_ierr, ierr, ncid)) return
      call nc_put_att_global_int(ncid, "nx_local", decomp%nx_local, ierr=local_ierr)
      if (.not. restart_io_ok(local_ierr, ierr, ncid)) return
      call nc_put_att_global_int(ncid, "ny_local", decomp%ny_local, ierr=local_ierr)
      if (.not. restart_io_ok(local_ierr, ierr, ncid)) return

      ! Grid / vcoord / tracer metadata (validated on read — review #7).
      call nc_put_att_global_int(ncid, "nz_ml", meta%nz_ml, ierr=local_ierr)
      if (.not. restart_io_ok(local_ierr, ierr, ncid)) return
      call nc_put_att_global_int(ncid, "nghost", meta%nghost, ierr=local_ierr)
      if (.not. restart_io_ok(local_ierr, ierr, ncid)) return
      call nc_put_att_global_int(ncid, "vcoord_type", meta%vcoord_type, ierr=local_ierr)
      if (.not. restart_io_ok(local_ierr, ierr, ncid)) return
      call nc_put_att_global(ncid, "vcoord_name", trim(meta%vcoord_name), ierr=local_ierr)
      if (.not. restart_io_ok(local_ierr, ierr, ncid)) return
      call nc_put_att_global_int(ncid, "n_tracers", meta%n_tracers, ierr=local_ierr)
      if (.not. restart_io_ok(local_ierr, ierr, ncid)) return
      call nc_put_att_global_int(ncid, "idx_salinity", meta%idx_salinity, ierr=local_ierr)
      if (.not. restart_io_ok(local_ierr, ierr, ncid)) return
      call nc_put_att_global_int(ncid, "idx_temperature", meta%idx_temperature, ierr=local_ierr)
      if (.not. restart_io_ok(local_ierr, ierr, ncid)) return

      call nc_enddef(ncid, ierr=local_ierr)
      if (.not. restart_io_ok(local_ierr, ierr, ncid)) return

      ! Write each FULL local array (interior + ghosts) / scalar.
      do e = 1, reg%n
         associate (en => reg%entries(e))
            select case (en%rank)
            case (0)
               call nc_check(nf90_put_var(ncid, varids(e), en%p0), &
                             "writing restart scalar "//trim(en%tag), local_ierr)
            case (2)
               call nc_check(nf90_put_var(ncid, varids(e), &
                                          en%p2(en%ng + 1:en%ng + en%nx_phys, en%ng + 1:en%ng + en%ny_phys)), &
                             "writing restart "//trim(en%tag), local_ierr)
            case default
               call nc_check(nf90_put_var(ncid, varids(e), &
                                          en%p3(en%ng + 1:en%ng + en%nx_phys, en%ng + 1:en%ng + en%ny_phys, 1:en%nk)), &
                             "writing restart "//trim(en%tag), local_ierr)
            end select
         end associate
         if (.not. restart_io_ok(local_ierr, ierr, ncid)) return
      end do

      call nc_check(nf90_put_var(ncid, time_varid, t), "writing restart time", local_ierr)
      if (.not. restart_io_ok(local_ierr, ierr, ncid)) return
      call nc_check(nf90_put_var(ncid, step_varid, step), "writing restart n_steps", local_ierr)
      if (.not. restart_io_ok(local_ierr, ierr, ncid)) return
      call nc_check(nf90_put_var(ncid, osc_varid, outer_step_count), &
                    "writing restart outer_step_count", local_ierr)
      if (.not. restart_io_ok(local_ierr, ierr, ncid)) return

      call nc_close(ncid, ierr=local_ierr)
      if (.not. restart_io_ok(local_ierr, ierr)) return
      deallocate (varids)

      ! Atomic swap: rename the complete tmp file over the target.  A
      ! crash before this point leaves the previous checkpoint untouched.
      rc = posix_rename(tmpname, trim(filename))
      if (rc /= 0) then
         msg = "Failed to rename restart "//trim(tmpname)//" -> "// &
               trim(filename)//" (rc="//to_string(rc)//")"
         call error_ring_push(msg)
         call logger%error(msg)
         if (present(ierr)) then
            ierr = OCEAN_STATUS_ERR_IO
            return
         end if
         error stop "ocean_restart_write_local: rename failed"
      end if
      call logger%info("Ocean restart written at t = "//to_string(t)//" s")
   end subroutine ocean_restart_write_local

   function restart_io_ok(local_ierr, ierr, ncid) result(ok)
      !! Translate a raw `nc_check`-style status (0 = ok) from one of the
      !! `nc_*` calls in `ocean_restart_write_local` into the caller's
      !! `ierr` contract. Mirrors `rdb_bathymetry`'s `bathy_io_ok` /
      !! `rdb_ocean_diag_netcdf`'s `diag_io_ok` (P0.1 F1 / P7 F1): closes
      !! `ncid` (when given) on failure so a mid-write error does not
      !! leak the file handle, sets `ierr = OCEAN_STATUS_ERR_IO` when
      !! present, or `error stop`s with the generic `nc_check` text when
      !! `ierr` is absent (byte-identical legacy behaviour).
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
   end function restart_io_ok

   integer function posix_rename(oldpath, newpath) result(rc)
      !! Atomically rename `oldpath` -> `newpath` via libc `rename(2)`.
      !! Returns 0 on success, -1 on error.
      character(len=*), intent(in) :: oldpath, newpath
      ! Null-terminated C strings as `target` locals so `c_loc` can hand
      ! libc `rename` a pointer to each (c_ptr by value sidesteps
      ! assumed-size dummies).
      character(kind=c_char), allocatable, target :: c_old(:), c_new(:)
      interface
         function c_rename(old, new) bind(C, name="rename") result(r)
            import :: c_ptr, c_int
            implicit none
            type(c_ptr), value, intent(in) :: old
            type(c_ptr), value, intent(in) :: new
            integer(c_int) :: r
         end function c_rename
      end interface
      c_old = string_to_c(oldpath)
      c_new = string_to_c(newpath)
      rc = int(c_rename(c_loc(c_old), c_loc(c_new)))
   end function posix_rename

   pure function string_to_c(s) result(c)
      !! Pack a Fortran string into a null-terminated c_char array.
      character(len=*), intent(in) :: s
      character(kind=c_char), allocatable :: c(:)
      integer :: i, n
      n = len_trim(s)
      allocate (c(n + 1))
      do i = 1, n
         c(i) = s(i:i)
      end do
      c(n + 1) = c_null_char
   end function string_to_c

   subroutine ocean_restart_check_decomp(filename, decomp, ierr, meta)
      !! Validate schema + decomposition (+ optional grid/vcoord/tracer)
      !! metadata in `filename` against the live `decomp`/`meta`.
      !! `ierr` = 0 on match; 1 = schema, 2 = decomposition,
      !! 3 = grid/vcoord/tracer mismatch.
      character(len=*), intent(in) :: filename
      type(decomp_t), intent(in) :: decomp
      integer, intent(out) :: ierr
      type(ocean_restart_metadata_t), intent(in), optional :: meta

      integer :: ncid, file_schema, fpx, fpy, fnxg, fnyg, fnxl, fnyl
      integer :: fis, fjs
      integer :: fnz, fng, fvc, fntr, fidxs, fidxt
      character(len=:), allocatable :: msg

      ierr = 0
      call nc_open_read(filename, ncid)
      call nc_get_att_int(ncid, "schema_version", file_schema)
      if (file_schema /= RESTART_SCHEMA_VERSION) then
         msg = "Restart schema mismatch: file v"//to_string(file_schema)// &
               " but build expects v"//to_string(RESTART_SCHEMA_VERSION)
         call error_ring_push(msg)
         call logger%error(msg)
         ierr = 1
         call nc_close(ncid)
         return
      end if
      call nc_get_att_int(ncid, "px", fpx)
      call nc_get_att_int(ncid, "py", fpy)
      call nc_get_att_int(ncid, "i_start", fis)
      call nc_get_att_int(ncid, "j_start", fjs)
      call nc_get_att_int(ncid, "nx_global", fnxg)
      call nc_get_att_int(ncid, "ny_global", fnyg)
      call nc_get_att_int(ncid, "nx_local", fnxl)
      call nc_get_att_int(ncid, "ny_local", fnyl)

      if (fpx /= decomp%px .or. fpy /= decomp%py .or. &
          fis /= decomp%i_start .or. fjs /= decomp%j_start .or. &
          fnxg /= decomp%nx_global .or. fnyg /= decomp%ny_global .or. &
          fnxl /= decomp%nx_local .or. fnyl /= decomp%ny_local) then
         msg = "Restart decomposition mismatch: file px="// &
               to_string(fpx)//" py="//to_string(fpy)// &
               " i_start="//to_string(fis)//" j_start="//to_string(fjs)// &
               " nx_local="//to_string(fnxl)//" ny_local="//to_string(fnyl)// &
               " but run has px="//to_string(decomp%px)// &
               " py="//to_string(decomp%py)// &
               " i_start="//to_string(decomp%i_start)// &
               " j_start="//to_string(decomp%j_start)// &
               " nx_local="//to_string(decomp%nx_local)// &
               " ny_local="//to_string(decomp%ny_local)// &
               ".  File: "//trim(filename)// &
               ".  Same-decomp resume only; use the offline "// &
               "redistribute tool to change rank count."
         call error_ring_push(msg)
         call logger%error(msg)
         ierr = 2
         call nc_close(ncid)
         return
      end if

      if (present(meta)) then
         call nc_get_att_int(ncid, "nz_ml", fnz)
         call nc_get_att_int(ncid, "nghost", fng)
         call nc_get_att_int(ncid, "vcoord_type", fvc)
         call nc_get_att_int(ncid, "n_tracers", fntr)
         call nc_get_att_int(ncid, "idx_salinity", fidxs)
         call nc_get_att_int(ncid, "idx_temperature", fidxt)
         if (fnz /= meta%nz_ml .or. fng /= meta%nghost .or. &
             fvc /= meta%vcoord_type .or. fntr /= meta%n_tracers .or. &
             fidxs /= meta%idx_salinity .or. fidxt /= meta%idx_temperature) then
            msg = "Restart grid/vcoord/tracer mismatch: file nz_ml="// &
                  to_string(fnz)//" nghost="//to_string(fng)// &
                  " vcoord_type="//to_string(fvc)// &
                  " n_tracers="//to_string(fntr)// &
                  " idx_S="//to_string(fidxs)//" idx_T="//to_string(fidxt)// &
                  " but run has nz_ml="//to_string(meta%nz_ml)// &
                  " nghost="//to_string(meta%nghost)// &
                  " vcoord_type="//to_string(meta%vcoord_type)// &
                  " ("//trim(meta%vcoord_name)//")"// &
                  " n_tracers="//to_string(meta%n_tracers)// &
                  " idx_S="//to_string(meta%idx_salinity)// &
                  " idx_T="//to_string(meta%idx_temperature)// &
                  ".  File: "//trim(filename)
            call error_ring_push(msg)
            call logger%error(msg)
            ierr = 3
         end if
      end if
      call nc_close(ncid)
   end subroutine ocean_restart_check_decomp

   subroutine ocean_restart_read_local(filename, reg, decomp, meta, t, step, &
                                       outer_step_count, ierr)
      !! Read a per-rank ocean restart into the registry's host arrays
      !! (FULL local extent) + return the scalar metadata.  Validates
      !! decomp + grid/vcoord/tracer metadata first and error-stops on
      !! mismatch unless `ierr` is present.  A REQUIRED field absent from
      !! the file is FATAL; only `optional` entries warn-and-seed.  Must
      !! run BEFORE `ocean_state_enter_data` so the H->D copy carries the
      !! restored values up.
      character(len=*), intent(in) :: filename
      type(restart_registry_t), intent(inout) :: reg
      type(decomp_t), intent(in) :: decomp
      type(ocean_restart_metadata_t), intent(in) :: meta
      real(wp), intent(out) :: t
      integer, intent(out) :: step
      integer, intent(out) :: outer_step_count
      integer, intent(out), optional :: ierr

      integer :: ncid, e, varid, status, check_err
      character(len=:), allocatable :: msg

      t = 0.0_wp
      step = 0
      outer_step_count = 0
      if (present(ierr)) ierr = 0

      call logger%info("Reading ocean restart file: "//trim(filename))

      call ocean_restart_check_decomp(filename, decomp, check_err, meta=meta)
      if (check_err /= 0) then
         if (present(ierr)) then
            ierr = check_err
            return
         end if
         error stop "ocean_restart_read_local: restart metadata mismatch (see log)"
      end if

      call nc_open_read(filename, ncid)
      do e = 1, reg%n
         associate (en => reg%entries(e))
            status = nf90_inq_varid(ncid, trim(en%tag), varid)
            if (status /= nf90_noerr) then
               if (en%optional) then
                  call logger%warning("Restart file lacks optional field "// &
                                      trim(en%tag)//"; leaving it at its seeded value")
                  cycle
               end if
               msg = "Restart file "//trim(filename)// &
                     " is missing REQUIRED field "//trim(en%tag)// &
                     "; refusing to cold-start half the state."
               call error_ring_push(msg)
               call logger%error(msg)
               call nc_close(ncid)
               error stop "ocean_restart_read_local: missing required field"
            end if
            select case (en%rank)
            case (0)
               call nc_check(nf90_get_var(ncid, varid, en%p0), &
                             "reading restart scalar "//trim(en%tag))
            case (2)
               call nc_check(nf90_get_var(ncid, varid, &
                                          en%p2(en%ng + 1:en%ng + en%nx_phys, en%ng + 1:en%ng + en%ny_phys)), &
                             "reading restart "//trim(en%tag))
            case default
               call nc_check(nf90_get_var(ncid, varid, &
                                          en%p3(en%ng + 1:en%ng + en%nx_phys, en%ng + 1:en%ng + en%ny_phys, 1:en%nk)), &
                             "reading restart "//trim(en%tag))
            end select
         end associate
      end do

      call nc_get_varid(ncid, "time", varid)
      call nc_check(nf90_get_var(ncid, varid, t), "reading restart time")
      call nc_get_varid(ncid, "n_steps", varid)
      call nc_check(nf90_get_var(ncid, varid, step), "reading restart n_steps")
      call nc_get_varid(ncid, "outer_step_count", varid)
      call nc_check(nf90_get_var(ncid, varid, outer_step_count), &
                    "reading restart outer_step_count")
      call nc_close(ncid)

      call logger%info("Ocean restart loaded: t = "//to_string(t)//" s, step = "// &
                       to_string(step))
   end subroutine ocean_restart_read_local

   subroutine nc_get_varid(ncid, name, varid)
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: name
      integer, intent(out) :: varid
      call nc_check(nf90_inq_varid(ncid, name, varid), "looking up restart var "//name)
   end subroutine nc_get_varid

end module rdb_ocean_restart_io
