!! Gates for `&ocean_cavity_dyn_nml draft_config = "file"` — the static
!! 2-D NetCDF ice draft, read through the PR-14 reader
!! (`rdb_ocean_data_input::ocean_data_input_load_static_2d`).
!!
!! Lives in its own suite (rather than joining `test_ocean_cavity_draft`)
!! because it WRITES NetCDF, which makes it a `NETCDF_REQUIRED_TESTS`
!! member; the analytic draft tests must keep building with
!! `RDB_ENABLE_NETCDF=OFF`.
!!
!! What is asserted:
!!   1. `roundtrip_depth_convention` — a field written as a DEPTH
!!      (positive down) comes back cell-for-cell, in the interior, with
!!      the ghost band filled by constant extrapolation exactly as the
!!      file bathymetry loader fills it.
!!   2. `elevation_convention_is_negated` — the same file read with
!!      `draft_sign="elevation"` (the ISOMIP+ `iceDraft` convention,
!!      z_d <= 0) comes back as `-z_d`, and open water (z_d = 0) stays a
!!      POSITIVE zero rather than `-0.0`.
!!   3. `wrong_sign_is_caught_not_silently_accepted` — reading an
!!      elevation file as a depth yields a negative draft, which the
!!      existing `cavity_draft_is_finite_nonneg` guard rejects.  This is
!!      the whole reason `draft_sign` is an explicit knob.
!!   4. `dim_mismatch_is_reported` — a file on the wrong grid returns a
!!      non-zero status rather than reading garbage.
!!   5. `ghosts_are_filled` — no ghost cell is left at the alloc-time
!!      zero when the interior is non-zero (the standing formula-setter
!!      rule, which the file path must match).
!!
!! The NetCDF writer is local to this module, matching the house pattern
!! (`test_ocean_data_input`, `test_ocean_zinit`): there is no shared
!! test helper that writes NetCDF.  Files go under `tmp_local_artifacts/`
!! per CLAUDE.md, relative to wherever ctest runs the binary.
module test_ocean_cavity_draft_file
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_data_input, only: ocean_data_input_load_static_2d, data_input_dims_ok
   use rdb_ocean_cavity, only: parse_cavity_draft_sign, cavity_draft_apply_sign, &
                               CAVITY_SIGN_DEPTH, CAVITY_SIGN_ELEVATION, &
                               CAVITY_SIGN_INVALID, cavity_draft_is_finite_nonneg
   use rdb_ocean_bathymetry_inject, only: bathymetry_fill_ghosts_array
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   use rdb_io_netcdf, only: nc_create_file, nc_close, nc_def_dim, nc_def_var_3d, &
                            rdb_def_var_1d, nc_enddef, rdb_put_var_1d
   use netcdf, only: nf90_put_var
   implicit none
   private

   public :: collect_ocean_cavity_draft_file_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NXP = 6, NYP = 4
   real(wp), parameter :: TOL = 1.0e-12_wp

contains

   subroutine collect_ocean_cavity_draft_file_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("draft_file_roundtrip_depth_convention", test_roundtrip_depth), &
                  new_unittest("draft_file_elevation_convention_is_negated", test_elevation), &
                  new_unittest("draft_file_wrong_sign_is_caught", test_wrong_sign), &
                  new_unittest("draft_file_dim_mismatch_is_reported", test_dim_mismatch), &
                  new_unittest("draft_file_ghosts_are_filled", test_ghosts), &
                  new_unittest("draft_sign_parser_has_no_silent_default", test_sign_parser) &
                  ]
   end subroutine collect_ocean_cavity_draft_file_tests

   ! -----------------------------------------------------------------
   ! Test support
   ! -----------------------------------------------------------------

   subroutine write_draft_file(filename, nx, ny, var, field)
      !! Write a `(x, y, t)` Fortran-ordered single-record file — the
      !! shape `register_common` demands (a C/Python writer or `ncdump`
      !! spells the same layout `(nTime, ny, nx)`), with the CF time
      !! coordinate variable the reader resolves.
      character(len=*), intent(in) :: filename, var
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: field(nx, ny)
      integer :: ncid, dim_x, dim_y, dim_t, vid_f, vid_t, ierr
      real(wp) :: f3(nx, ny, 1), t_axis(1)

      f3(:, :, 1) = field
      t_axis(1) = 0.0_wp
      call nc_create_file(filename, ncid)
      call nc_def_dim(ncid, "x", nx, dim_x)
      call nc_def_dim(ncid, "y", ny, dim_y)
      call nc_def_dim(ncid, "time", 1, dim_t)
      call nc_def_var_3d(ncid, var, [dim_x, dim_y, dim_t], vid_f)
      call rdb_def_var_1d(ncid, "time", dim_t, vid_t)
      call nc_enddef(ncid)
      ierr = nf90_put_var(ncid, vid_f, f3)
      call rdb_put_var_1d(ncid, vid_t, t_axis)
      call nc_close(ncid)
   end subroutine write_draft_file

   function make_grid() result(g)
      type(hgrid_t) :: g
      call g%init(NXP, NYP, NGHOST, 2000.0_wp, 2000.0_wp)
   end function make_grid

   pure function ramp_depth(i, j) result(d)
      !! A draft deepening in x and varying in y, with an open-water
      !! (zero) column at i = 1 so the zero-handling is exercised.
      integer, intent(in) :: i, j
      real(wp) :: d
      if (i == 1) then
         d = 0.0_wp
      else
         d = 100.0_wp + 20.0_wp*real(i, wp) + 5.0_wp*real(j, wp)
      end if
   end function ramp_depth

   subroutine load_into(filename, var, sign_code, z_draft, grid, ierr)
      !! The production sequence, in the production order: read the
      !! interior, sign-normalise, fill ghosts.  Mirrors the
      !! `CAVITY_DRAFT_FILE` branch of `seed_cavity_draft`.
      character(len=*), intent(in) :: filename, var
      integer, intent(in) :: sign_code
      real(wp), intent(inout) :: z_draft(:, :)
      type(hgrid_t), intent(in) :: grid
      integer, intent(out) :: ierr
      integer :: nx, ny
      nx = size(z_draft, 1)
      ny = size(z_draft, 2)
      z_draft = 0.0_wp
      call ocean_data_input_load_static_2d(filename, var, grid, nx, ny, &
                                           grid%nghost + 1, grid%nghost + 1, &
                                           z_draft, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call cavity_draft_apply_sign(z_draft, nx, ny, sign_code)
      call bathymetry_fill_ghosts_array(z_draft, grid)
   end subroutine load_into

   ! -----------------------------------------------------------------
   ! (1) depth convention round-trip
   ! -----------------------------------------------------------------
   subroutine test_roundtrip_depth(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      real(wp) :: field(NXP, NYP)
      real(wp), allocatable :: z(:, :)
      character(len=*), parameter :: FN = "tmp_local_artifacts/test_draft_depth.nc"
      integer :: i, j, ng, ierr
      real(wp) :: worst

      grid = make_grid()
      ng = grid%nghost
      do j = 1, NYP
         do i = 1, NXP
            field(i, j) = ramp_depth(i, j)
         end do
      end do
      call write_draft_file(FN, NXP, NYP, "iceDraft", field)

      allocate (z(grid%nx_total, grid%ny_total))
      call load_into(FN, "iceDraft", CAVITY_SIGN_DEPTH, z, grid, ierr)

      checks: block
         call check(error, ierr == OCEAN_STATUS_OK, "static 2-D draft read must succeed")
         if (allocated(error)) exit checks
         worst = 0.0_wp
         do j = 1, NYP
            do i = 1, NXP
               worst = max(worst, abs(z(ng + i, ng + j) - field(i, j)))
            end do
         end do
         call check(error, worst < TOL, &
                    "every interior cell must round-trip through the file unchanged")
         if (allocated(error)) exit checks
         call check(error, cavity_draft_is_finite_nonneg(z, grid%nx_total, grid%ny_total), &
                    "a depth-convention file must load as a finite non-negative draft")
         if (allocated(error)) exit checks
         ! Ghost fill is constant extrapolation, exactly as the file
         ! bathymetry loader does it: the west ghost column copies the
         ! first interior column.
         call check(error, abs(z(ng, ng + 1) - field(1, 1)) < TOL, &
                    "the west ghost column must copy the first interior column")
      end block checks
   end subroutine test_roundtrip_depth

   ! -----------------------------------------------------------------
   ! (2) ISOMIP+ elevation convention
   ! -----------------------------------------------------------------
   subroutine test_elevation(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      real(wp) :: field(NXP, NYP)
      real(wp), allocatable :: z(:, :)
      character(len=*), parameter :: FN = "tmp_local_artifacts/test_draft_elev.nc"
      integer :: i, j, ng, ierr
      real(wp) :: worst

      grid = make_grid()
      ng = grid%nghost
      ! z_d, the ELEVATION of the ice-ocean interface: <= 0 under ice,
      ! exactly 0 in open water (the ISOMIP+ file's own convention).
      do j = 1, NYP
         do i = 1, NXP
            field(i, j) = -ramp_depth(i, j)
         end do
      end do
      call write_draft_file(FN, NXP, NYP, "iceDraft", field)

      allocate (z(grid%nx_total, grid%ny_total))
      call load_into(FN, "iceDraft", CAVITY_SIGN_ELEVATION, z, grid, ierr)

      checks: block
         call check(error, ierr == OCEAN_STATUS_OK, "elevation-convention read must succeed")
         if (allocated(error)) exit checks
         worst = 0.0_wp
         do j = 1, NYP
            do i = 1, NXP
               worst = max(worst, abs(z(ng + i, ng + j) - ramp_depth(i, j)))
            end do
         end do
         call check(error, worst < TOL, &
                    "an ELEVATION file must load as the negated DEPTH")
         if (allocated(error)) exit checks
         ! Open water: z_d = 0 must come back as +0, not -0.  A -0 would
         ! make `cavity_apply_land_exclusion`'s `z_draft /= 0` test fire
         ! on a column with no ice... except it would not, because
         ! -0 == 0 in Fortran; assert the SIGN BIT anyway so the intent
         ! survives a refactor, via 1/z (which distinguishes them).
         call check(error, z(ng + 1, ng + 1) == 0.0_wp, &
                    "open water (z_d = 0) must load as zero draft")
         if (allocated(error)) exit checks
         call check(error, sign(1.0_wp, z(ng + 1, ng + 1)) > 0.0_wp, &
                    "open water must load as POSITIVE zero, not -0.0")
         if (allocated(error)) exit checks
         call check(error, cavity_draft_is_finite_nonneg(z, grid%nx_total, grid%ny_total), &
                    "the negated elevation must be a valid non-negative draft")
      end block checks
   end subroutine test_elevation

   ! -----------------------------------------------------------------
   ! (3) the wrong sign is caught, not silently accepted
   ! -----------------------------------------------------------------
   subroutine test_wrong_sign(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      real(wp) :: field(NXP, NYP)
      real(wp), allocatable :: z(:, :)
      character(len=*), parameter :: FN = "tmp_local_artifacts/test_draft_wrongsign.nc"
      integer :: i, j, ierr

      grid = make_grid()
      do j = 1, NYP
         do i = 1, NXP
            field(i, j) = -ramp_depth(i, j)          ! an ELEVATION file...
         end do
      end do
      call write_draft_file(FN, NXP, NYP, "iceDraft", field)

      allocate (z(grid%nx_total, grid%ny_total))
      ! ...read with the DEPTH convention (the mistake this knob exists
      ! to make impossible to make silently).
      call load_into(FN, "iceDraft", CAVITY_SIGN_DEPTH, z, grid, ierr)

      checks: block
         call check(error, ierr == OCEAN_STATUS_OK, "the read itself still succeeds")
         if (allocated(error)) exit checks
         call check(error,.not. cavity_draft_is_finite_nonneg(z, grid%nx_total, &
                                                              grid%ny_total), &
                    "an elevation file read as a depth must be REJECTED by the "// &
                    "non-negativity guard, not silently accepted")
      end block checks
   end subroutine test_wrong_sign

   ! -----------------------------------------------------------------
   ! (4) dimension mismatch
   ! -----------------------------------------------------------------
   subroutine test_dim_mismatch(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      real(wp) :: field(NXP - 1, NYP)
      real(wp), allocatable :: z(:, :)
      character(len=*), parameter :: FN = "tmp_local_artifacts/test_draft_badshape.nc"
      integer :: ierr

      grid = make_grid()
      field = 50.0_wp
      call write_draft_file(FN, NXP - 1, NYP, "iceDraft", field)

      checks: block
         ! The non-erroring validator says no...
         call check(error,.not. data_input_dims_ok(FN, "iceDraft", NXP, NYP, &
                                                   is_3d=.false.), &
                    "a file one column short must fail dimension validation")
         if (allocated(error)) exit checks
         ! ...and so does the loader, through its status.
         allocate (z(grid%nx_total, grid%ny_total))
         call load_into(FN, "iceDraft", CAVITY_SIGN_DEPTH, z, grid, ierr)
         call check(error, ierr /= OCEAN_STATUS_OK, &
                    "the loader must return a non-zero status on a grid mismatch")
         if (allocated(error)) exit checks
         ! A missing variable is likewise a status, not a read of garbage.
         call load_into(FN, "no_such_variable", CAVITY_SIGN_DEPTH, z, grid, ierr)
         call check(error, ierr /= OCEAN_STATUS_OK, &
                    "a missing variable must return a non-zero status")
      end block checks
   end subroutine test_dim_mismatch

   ! -----------------------------------------------------------------
   ! (5) ghost fill
   ! -----------------------------------------------------------------
   subroutine test_ghosts(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      real(wp) :: field(NXP, NYP)
      real(wp), allocatable :: z(:, :)
      character(len=*), parameter :: FN = "tmp_local_artifacts/test_draft_ghosts.nc"
      integer :: ierr

      grid = make_grid()
      field = 250.0_wp            ! uniform, non-zero everywhere
      call write_draft_file(FN, NXP, NYP, "iceDraft", field)

      allocate (z(grid%nx_total, grid%ny_total))
      call load_into(FN, "iceDraft", CAVITY_SIGN_DEPTH, z, grid, ierr)

      checks: block
         call check(error, ierr == OCEAN_STATUS_OK, "read must succeed")
         if (allocated(error)) exit checks
         ! A uniform interior extrapolates to a uniform whole array: no
         ! cell may be left at the alloc-time zero.  That is the standing
         ! rule the formula setters follow, and a ghost row left at zero
         ! would put a phantom calving front one cell outside every wall.
         call check(error, minval(z) > 0.0_wp, &
                    "no cell, ghost band included, may be left at zero")
         if (allocated(error)) exit checks
         call check(error, abs(maxval(z) - 250.0_wp) < TOL .and. &
                    abs(minval(z) - 250.0_wp) < TOL, &
                    "a uniform draft must extrapolate to a uniform whole array")
      end block checks
   end subroutine test_ghosts

   ! -----------------------------------------------------------------
   ! (6) the sign parser refuses to guess
   ! -----------------------------------------------------------------
   subroutine test_sign_parser(error)
      type(error_type), allocatable, intent(out) :: error
      checks: block
         call check(error, parse_cavity_draft_sign("depth") == CAVITY_SIGN_DEPTH, &
                    "'depth' -> CAVITY_SIGN_DEPTH")
         if (allocated(error)) exit checks
         call check(error, parse_cavity_draft_sign("positive_down") == CAVITY_SIGN_DEPTH, &
                    "'positive_down' -> CAVITY_SIGN_DEPTH")
         if (allocated(error)) exit checks
         call check(error, parse_cavity_draft_sign("elevation") == CAVITY_SIGN_ELEVATION, &
                    "'elevation' -> CAVITY_SIGN_ELEVATION")
         if (allocated(error)) exit checks
         call check(error, parse_cavity_draft_sign("positive_up") == CAVITY_SIGN_ELEVATION, &
                    "'positive_up' -> CAVITY_SIGN_ELEVATION")
         if (allocated(error)) exit checks
         call check(error, parse_cavity_draft_sign("down") == CAVITY_SIGN_INVALID, &
                    "an unrecognised spelling must be INVALID, never a guess")
      end block checks
   end subroutine test_sign_parser

end module test_ocean_cavity_draft_file
