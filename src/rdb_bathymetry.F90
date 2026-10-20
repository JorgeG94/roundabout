!! Bathymetry loading from NetCDF files
module rdb_bathymetry
   !! Reads bathymetry (bottom elevation) from a NetCDF file and populates
   !! the state%barotropic%b array. Expects a 2D variable on a structured grid.
   !! If no file is configured, falls back to flat bottom.
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_io_netcdf, only: nc_check, nc_open_read, nc_close, &
                            nc_get_dim_len, nc_get_varid, nc_get_var_slab_2d
   use netcdf, only: nf90_inquire_variable, nf90_inquire_dimension, nf90_noerr
   use pic_logger, only: logger => global_logger
   use pic_strings, only: to_string
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_IO
   implicit none
   private

   public :: load_bathymetry_into_array
   public :: fill_bathymetry_ghosts_array

contains

   subroutine load_bathymetry_into_array(filename, b, grid, ierr)
      !! Load bathymetry from a NetCDF file into a target 2D array
      !! shaped `(grid%nx_total, grid%ny_total)`.
      !!
      !! The file must contain:
      !!   - Dimensions: x (nx_phys), y (ny_phys)
      !!   - Variable: b(x, y) or elevation(x, y) or depth(x, y)
      !!
      !! Handles NetCDF C/Fortran dimension reversal: files written by
      !! Python/C store b(x,y) in C order, which Fortran reads as b(y,x).
      !! If the direct dimensions don't match but the transposed ones do,
      !! we read transposed and copy correctly.
      !!
      !! Bathymetry lands in the interior; ghost cells are filled by
      !! constant extrapolation from the nearest interior cell via
      !! `fill_bathymetry_ghosts_array`.
      !!
      !! The file holds the WHOLE grid (`grid%nx_global x grid%ny_global`).
      !! A decomposed tile reads only the full-width band of rows its
      !! storage covers (physical rows +- nghost, clipped to the grid) with a
      !! start/count window, extrapolates that band's ghosts exactly as the
      !! whole grid would (at a GLOBAL edge the band edge is the grid edge;
      !! elsewhere the band's rows ARE the tile's ghost rows, read from the
      !! file), and keeps its own columns: every tile is bit-identical to
      !! its slice of the single-rank array, seam ghosts included.
      character(len=*), intent(in) :: filename
      real(wp), intent(inout) :: b(:, :)
      type(hgrid_t), intent(in) :: grid
      integer, intent(out), optional :: ierr
         !! Non-zero on a missing/unreadable file or a dimension mismatch
         !! when present; absent behaves as today (`error stop`).

      integer :: ncid, varid, file_nx, file_ny
      integer :: ng, i, j, ni, nj, jb0, jb1, nb, jo_b
      type(hgrid_t) :: gband
      real(wp), allocatable :: bband(:, :)
      integer :: var_dimids(2)
      integer :: var_ndims, dim1_len, dim2_len
      integer :: local_ierr
      character(len=64) :: dim1_name, dim2_name
      logical :: needs_transpose
      real(wp), allocatable :: b_interior(:, :)

      call logger%info("Loading bathymetry from: "//trim(filename))

      call nc_open_read(filename, ncid, ierr=local_ierr)
      if (.not. bathy_io_ok(local_ierr, ierr)) return

      ! Read grid dimensions (by name — these are the logical sizes)
      call nc_get_dim_len(ncid, "x", file_nx, ierr=local_ierr)
      if (.not. bathy_io_ok(local_ierr, ierr, ncid)) return
      call nc_get_dim_len(ncid, "y", file_ny, ierr=local_ierr)
      if (.not. bathy_io_ok(local_ierr, ierr, ncid)) return

      ! Try variable names in order: b, elevation, depth
      call try_get_bathymetry_var(ncid, varid, ierr=local_ierr)
      if (.not. bathy_io_ok(local_ierr, ierr, ncid)) return

      ! Query the variable's actual dimension order in the file.
      ! Files written by Python/C store b(x,y) in C-order. Fortran's
      ! NetCDF library reverses this, so the variable appears as b(y,x)
      ! in Fortran. We detect this by checking the first dimension name.
      call nc_check(nf90_inquire_variable(ncid, varid, ndims=var_ndims, &
                                          dimids=var_dimids), &
                    "querying bathymetry variable", local_ierr)
      if (.not. bathy_io_ok(local_ierr, ierr, ncid)) return
      call nc_check(nf90_inquire_dimension(ncid, var_dimids(1), &
                                           name=dim1_name, len=dim1_len), &
                    "querying dim 1", local_ierr)
      if (.not. bathy_io_ok(local_ierr, ierr, ncid)) return
      call nc_check(nf90_inquire_dimension(ncid, var_dimids(2), &
                                           name=dim2_name, len=dim2_len), &
                    "querying dim 2", local_ierr)
      if (.not. bathy_io_ok(local_ierr, ierr, ncid)) return

      call logger%info("Bathymetry variable: "//trim(dim1_name)//"="// &
                       to_string(dim1_len)//" x "//trim(dim2_name)//"="// &
                       to_string(dim2_len)//" (Fortran order)")

      ! In Fortran, dim1 is the first array index (varies fastest).
      ! If dim1 is "x" and matches nx, no transpose needed.
      ! If dim1 is "y" (C-order file reversed), we need to transpose.
      needs_transpose = (trim(dim1_name) == "y")

      ! The file describes the WHOLE grid.
      ni = grid%nx_global
      nj = grid%ny_global

      ! Validate dimensions
      if (needs_transpose) then
         if (dim1_len /= nj .or. dim2_len /= ni) then
            call logger%error("Bathymetry grid mismatch: file has "// &
                              to_string(dim2_len)//" x "//to_string(dim1_len)// &
                              " but simulation expects "// &
                              to_string(ni)//" x "// &
                              to_string(nj))
            call nc_close(ncid)
            if (present(ierr)) then
               ierr = OCEAN_STATUS_ERR_IO
               return
            end if
            error stop "Bathymetry grid mismatch"
         end if
      else
         if (dim1_len /= ni .or. dim2_len /= nj) then
            call logger%error("Bathymetry grid mismatch: file has "// &
                              to_string(dim1_len)//" x "//to_string(dim2_len)// &
                              " but simulation expects "// &
                              to_string(ni)//" x "// &
                              to_string(nj))
            call nc_close(ncid)
            if (present(ierr)) then
               ierr = OCEAN_STATUS_ERR_IO
               return
            end if
            error stop "Bathymetry grid mismatch"
         end if
      end if

      ng = grid%nghost
      ! Rows of the whole grid this tile's storage covers (all of them on an
      ! undecomposed grid), full width.
      jb0 = max(1, grid%j_offset_global + 1 - ng)
      jb1 = min(nj, grid%j_offset_global + grid%ny_phys + ng)
      nb = jb1 - jb0 + 1

      ! Read the band into a temporary matching the Fortran storage order.
      if (needs_transpose) then
         allocate (b_interior(nb, ni))
         call nc_get_var_slab_2d(ncid, varid, [jb0, 1], [nb, ni], b_interior, ierr=local_ierr)
      else
         allocate (b_interior(ni, nb))
         call nc_get_var_slab_2d(ncid, varid, [1, jb0], [ni, nb], b_interior, ierr=local_ierr)
      end if
      if (.not. bathy_io_ok(local_ierr, ierr, ncid)) return
      call nc_close(ncid)

      ! Place the band in a whole-width, band-height array with ghosts and
      ! extrapolate ITS ghosts -- exactly the undecomposed construction when
      ! the band is the whole grid.
      call gband%init(ni, nb, ng, grid%dx, grid%dy)
      allocate (bband(gband%nx_total, gband%ny_total))
      bband = 0.0_wp
      if (needs_transpose) then
         ! b_interior is (y, x) in Fortran — transpose to (x, y)
         do j = 1, nb
            do i = 1, ni
               bband(ng + i, ng + j) = b_interior(j, i)
            end do
         end do
      else
         bband(ng + 1:ng + ni, ng + 1:ng + nb) = b_interior
      end if
      deallocate (b_interior)

      ! Fill ghost cells by constant extrapolation from nearest interior cell.
      ! Without this, ghost cells retain b=0 which creates artificial cliffs
      ! against real bathymetry (e.g. b=-50m interior vs b=0 ghost), generating
      ! extreme velocities and tiny CFL timesteps.
      call fill_bathymetry_ghosts_array(bband, gband)

      ! Cut this tile's storage window (ghosts included) out of the band.
      jo_b = grid%j_offset_global - (jb0 - 1)
      do j = 1, grid%ny_total
         do i = 1, grid%nx_total
            b(i, j) = bband(i + grid%i_offset_global, j + jo_b)
         end do
      end do
      deallocate (bband)

      call logger%info("Bathymetry loaded: min = "// &
                       to_string(minval(b(ng + 1:ng + grid%nx_phys, &
                                          ng + 1:ng + grid%ny_phys)))// &
                       " m, max = "// &
                       to_string(maxval(b(ng + 1:ng + grid%nx_phys, &
                                          ng + 1:ng + grid%ny_phys)))// &
                       " m")

      if (present(ierr)) ierr = OCEAN_STATUS_OK
   end subroutine load_bathymetry_into_array

   subroutine try_get_bathymetry_var(ncid, varid, ierr)
      !! Try to find the bathymetry variable by common names
      use netcdf, only: nf90_inq_varid
      integer, intent(in) :: ncid
      integer, intent(out) :: varid
      integer, intent(out), optional :: ierr
         !! Non-zero when none of the recognised variable names is found,
         !! when present; absent behaves as today (leaves `varid`
         !! unset and lets the caller's next NetCDF call fail).

      integer :: status

      ! Try "b" first (our convention)
      status = nf90_inq_varid(ncid, "b", varid)
      if (status == nf90_noerr) then
         if (present(ierr)) ierr = OCEAN_STATUS_OK
         return
      end if

      ! Try "elevation"
      status = nf90_inq_varid(ncid, "elevation", varid)
      if (status == nf90_noerr) then
         if (present(ierr)) ierr = OCEAN_STATUS_OK
         return
      end if

      ! Try "depth"
      status = nf90_inq_varid(ncid, "depth", varid)
      if (status == nf90_noerr) then
         if (present(ierr)) ierr = OCEAN_STATUS_OK
         return
      end if

      ! None found
      call logger%error("No bathymetry variable found (tried: b, elevation, depth)")
      if (present(ierr)) then
         ierr = OCEAN_STATUS_ERR_IO
         return
      end if

   end subroutine try_get_bathymetry_var

   function bathy_io_ok(local_ierr, ierr, ncid) result(ok)
      !! Translate a raw `nc_check`-style status (0 = ok) from one of the
      !! `nc_*` reader calls in `load_bathymetry_into_array` into the
      !! caller's `ierr` contract: `.true.` on success; on failure,
      !! `.false.` with `ierr = OCEAN_STATUS_ERR_IO` when `ierr` is
      !! present (closing `ncid` first, when given, so a mid-read failure
      !! does not leak the file handle), or `error stop`s with the SAME
      !! generic text `nc_check` itself would have used had `ierr` never
      !! been threaded through — keeps the legacy (no `ierr`) behaviour
      !! byte-identical while unblocking the `ierr`-present return path
      !! (F1/F2 of the P0.1 review).
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
   end function bathy_io_ok

   subroutine fill_bathymetry_ghosts_array(b, grid)
      !! Fill ghost cell bathymetry by constant extrapolation from the
      !! nearest interior cell.  This ensures boundary flux computations
      !! see a consistent bottom elevation across the ghost-interior
      !! interface.  Works on a target 2D array directly so both the
      !! coastal and ocean barotropic slots can reuse it.
      real(wp), intent(inout) :: b(:, :)
      type(hgrid_t), intent(in) :: grid

      integer :: ng, i, j, nx, ny

      ng = grid%nghost
      nx = grid%nx_phys
      ny = grid%ny_phys

      ! West and east ghost columns
      do j = 1, grid%ny_total
         do i = 1, ng
            b(i, j) = b(ng + 1, j)                  ! west
            b(ng + nx + i, j) = b(ng + nx, j)       ! east
         end do
      end do

      ! South and north ghost rows (corners already filled above)
      do j = 1, ng
         do i = 1, grid%nx_total
            b(i, j) = b(i, ng + 1)                  ! south
            b(i, ng + ny + j) = b(i, ng + ny)       ! north
         end do
      end do

   end subroutine fill_bathymetry_ghosts_array

end module rdb_bathymetry
