!! Unit tests for the GFS_scale knob on the FV_MOM6 PGF kernel.
!!
!! MOM6's reduced-free-surface-gravity setup (`MOM_PressureForce_FV.F90:1922`)
!! adds a depth-independent Montgomery correction `dM(i,j) =
!! (GFS_scale − 1)·(g/ρ₀)·ρ_surf·η` to the per-layer PGF when
!! `GFS_scale < 1`.  This scales the SURFACE contribution of the slow PGF
!! by `GFS_scale`, leaving the BT substep (which now runs at
!! `g_bt = GFS_scale·g`) to drive η.  The depth-independent shape means
!! the baroclinic part is untouched and the BT mass-flux invariant
!! survives.
!!
!! Tests:
!!   1. `gfs_scale_unity_bit_identical` — gfs_scale=1.0 must produce
!!      PFu/PFv identical to the no-knob path (correction inactive).
!!   2. `gfs_scale_dampens_surface_gradient` — pure-SSH-tilt setup
!!      (NK=1, η sinusoidal, flat bath, uniform ρ).  PFu must scale
!!      linearly with gfs_scale: dpdx(gfs_scale=0.1) = 0.1 · dpdx(gfs_scale=1.0).
!!   3. `gfs_scale_correction_depth_independent` — stratified two-layer
!!      with non-trivial η.  The DIFFERENCE PFu(bed) − PFu(surf) must
!!      be the same with and without the gfs_scale correction (baroclinic
!!      part untouched).
module test_ocean_pgf_gfs_scale
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_ocean_pressure_force, only: ocean_pressure_force_t, &
                                       ocean_pressure_force_compute, &
                                       OPGF_VARIANT_FV_MOM6
   implicit none
   private

   public :: collect_ocean_pgf_gfs_scale_tests

   integer, parameter :: NGHOST = 2

contains

   subroutine collect_ocean_pgf_gfs_scale_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("gfs_scale_unity_bit_identical", &
                               test_gfs_scale_unity), &
                  new_unittest("gfs_scale_dampens_surface_gradient", &
                               test_gfs_scale_surface), &
                  new_unittest("gfs_scale_correction_depth_independent", &
                               test_gfs_scale_baroclinic) &
                  ]
   end subroutine collect_ocean_pgf_gfs_scale_tests

   subroutine run_pgf_with_b(ms, pgf, b_user, dx)
      !! Compute PGF with caller-supplied bathymetry (η = sum h − b).
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_pressure_force_t), intent(inout) :: pgf
      real(wp), intent(in) :: b_user(:, :)
      real(wp), intent(in) :: dx
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      call grid%init(size(ms%h_layer, 1) - 2*NGHOST, &
                     size(ms%h_layer, 2) - 2*NGHOST, &
                     NGHOST, dx, dx)
      call make_cartesian_metrics(metrics, grid)
      call pgf%set_bathymetry(b_user)
      !$acc enter data copyin(ms, pgf)
      call ms%enter_data()
      call pgf%enter_data()
      call ocean_pressure_force_compute(grid, metrics, pgf, ms)
      !$acc update self(pgf%dpdx_face%data, pgf%dpdy_face%data)
      call pgf%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, pgf)
      call destroy_cartesian_metrics(metrics)
   end subroutine run_pgf_with_b

   subroutine test_gfs_scale_unity(error)
      !! gfs_scale = 1.0 must leave the FV_MOM6 PGF bit-identical to
      !! the pre-knob path (no dM correction applied).  Compare two
      !! independent runs with the same setup.
      type(error_type), allocatable, intent(out) :: error
      type(multilayer_state_t) :: ms_a, ms_b
      type(ocean_pressure_force_t) :: pgf_a, pgf_b
      type(hgrid_t) :: grid
      real(wp), parameter :: DX = 5000.0_wp, H_PER = 500.0_wp
      real(wp), allocatable :: b_field(:, :)
      real(wp) :: max_diff
      integer :: i, nx, ny
      checks: block
         call grid%init(8, 6, NGHOST, DX, DX)
         nx = grid%nx_total; ny = grid%ny_total
         ms_a%nz_ml = 2; ms_b%nz_ml = 2
         call ms_a%init(grid); call ms_b%init(grid)
         call pgf_a%init(grid, nz_ml=2); call pgf_b%init(grid, nz_ml=2)
         pgf_a%variant = OPGF_VARIANT_FV_MOM6
         pgf_b%variant = OPGF_VARIANT_FV_MOM6
         pgf_a%rho0 = 1035.0_wp; pgf_a%rho_ref = 1035.0_wp
         pgf_b%rho0 = 1035.0_wp; pgf_b%rho_ref = 1035.0_wp
         pgf_a%gfs_scale = 1.0_wp     ! default
         pgf_b%gfs_scale = 1.0_wp     ! also default

         ! Non-trivial setup: stratified + sloped η.
         allocate (b_field(nx, ny), source=H_PER*2)
         do i = 1, nx
            ms_a%h_layer(i, :, 1) = H_PER + real(i, wp)*0.1_wp   ! bed
            ms_a%h_layer(i, :, 2) = H_PER + real(i, wp)*0.05_wp  ! surf
         end do
         ms_b%h_layer = ms_a%h_layer
         ms_a%rho_layer(:, :, 1) = 1036.0_wp; ms_a%rho_layer(:, :, 2) = 1035.0_wp
         ms_b%rho_layer = ms_a%rho_layer

         call run_pgf_with_b(ms_a, pgf_a, b_field, DX)
         call run_pgf_with_b(ms_b, pgf_b, b_field, DX)

         max_diff = max( &
                    maxval(abs(pgf_a%dpdx_face%data - pgf_b%dpdx_face%data)), &
                    maxval(abs(pgf_a%dpdy_face%data - pgf_b%dpdy_face%data)))
         call check(error, max_diff < 1.0e-14_wp, &
                    "gfs_scale=1 run not bit-identical to pre-knob path")
         deallocate (b_field)
      end block checks
      call pgf_a%destroy(); call pgf_b%destroy()
      call ms_a%destroy(); call ms_b%destroy()
   end subroutine test_gfs_scale_unity

   subroutine test_gfs_scale_surface(error)
      !! Pure-SSH-tilt: NK=1, uniform ρ, sinusoidal η.  PFu must reduce
      !! to -gfs_scale·g·∂η/∂x.  Verify by comparing two runs at
      !! gfs_scale=1.0 and gfs_scale=0.1.
      type(error_type), allocatable, intent(out) :: error
      type(multilayer_state_t) :: ms_a, ms_b
      type(ocean_pressure_force_t) :: pgf_a, pgf_b
      type(hgrid_t) :: grid
      real(wp), parameter :: DX = 1000.0_wp, B0 = 1000.0_wp
      real(wp), parameter :: ETA_AMP = 0.5_wp
      real(wp), parameter :: GFS_SCALE = 0.1_wp
      real(wp), parameter :: PI = 3.14159265358979323846_wp
      real(wp), allocatable :: b_field(:, :)
      integer :: i, j, nx, ny, NX_PHYS
      real(wp) :: eta_i, ratio, expected
      real(wp) :: dpdx_full, dpdx_scaled, dpdx_diff
      checks: block
         NX_PHYS = 16
         call grid%init(NX_PHYS, 4, NGHOST, DX, DX)
         nx = grid%nx_total; ny = grid%ny_total
         ms_a%nz_ml = 1; ms_b%nz_ml = 1
         call ms_a%init(grid); call ms_b%init(grid)
         call pgf_a%init(grid, nz_ml=1); call pgf_b%init(grid, nz_ml=1)
         pgf_a%variant = OPGF_VARIANT_FV_MOM6
         pgf_b%variant = OPGF_VARIANT_FV_MOM6
         pgf_a%rho0 = 1035.0_wp; pgf_a%rho_ref = 1035.0_wp
         pgf_b%rho0 = 1035.0_wp; pgf_b%rho_ref = 1035.0_wp
         pgf_a%gfs_scale = 1.0_wp
         pgf_b%gfs_scale = GFS_SCALE

         ! Uniform-ρ at rho_ref so the only signal is the SSH gradient.
         ms_a%rho_layer = 1035.0_wp
         ms_b%rho_layer = 1035.0_wp

         ! Sinusoidal η — η_i = ETA_AMP·sin(2π·i/NX_PHYS).
         allocate (b_field(nx, ny), source=B0)
         do j = 1, ny
            do i = 1, nx
               eta_i = ETA_AMP*sin(2.0_wp*PI*real(i - NGHOST, wp)/real(NX_PHYS, wp))
               ms_a%h_layer(i, j, 1) = B0 + eta_i
               ms_b%h_layer(i, j, 1) = B0 + eta_i
            end do
         end do

         call run_pgf_with_b(ms_a, pgf_a, b_field, DX)
         call run_pgf_with_b(ms_b, pgf_b, b_field, DX)

         ! Pick a deep-interior u-face away from walls.
         i = NGHOST + NX_PHYS/4; j = NGHOST + 2
         dpdx_full = pgf_a%dpdx_face%data(i, j, 1)
         dpdx_scaled = pgf_b%dpdx_face%data(i, j, 1)

         if (abs(dpdx_full) < 1.0e-10_wp) then
            call check(error, .false., &
                       "gfs_scale surface test: dpdx_full ~ 0, picked a bad face")
            exit checks
         end if

         ratio = dpdx_scaled/dpdx_full
         expected = GFS_SCALE
         call check(error, abs(ratio - expected) < 1.0e-6_wp, &
                    "gfs_scale=0.1: dpdx ratio not 0.1×")
         if (allocated(error)) exit checks

         ! And sanity: the absolute change should be ~ −0.9·g·∂η/∂x.
         dpdx_diff = dpdx_scaled - dpdx_full
         call check(error, abs(dpdx_diff) > 1.0e-8_wp, &
                    "gfs_scale=0.1: correction not applied (no change)")

         deallocate (b_field)
      end block checks
      call pgf_a%destroy(); call pgf_b%destroy()
      call ms_a%destroy(); call ms_b%destroy()
   end subroutine test_gfs_scale_surface

   subroutine test_gfs_scale_baroclinic(error)
      !! The dM correction is depth-independent (no k dependence in its
      !! formula).  So the DIFFERENCE PFu(k=bed) − PFu(k=surf) must be
      !! identical with and without the correction.  This proves the
      !! baroclinic mode is untouched — only the BT (depth-mean) part
      !! is scaled.
      type(error_type), allocatable, intent(out) :: error
      type(multilayer_state_t) :: ms_a, ms_b
      type(ocean_pressure_force_t) :: pgf_a, pgf_b
      type(hgrid_t) :: grid
      real(wp), parameter :: DX = 1000.0_wp, B0 = 1500.0_wp
      real(wp), allocatable :: b_field(:, :)
      integer :: i, j, nx, ny
      real(wp) :: bc_a, bc_b
      real(wp) :: max_bc_diff
      checks: block
         call grid%init(16, 4, NGHOST, DX, DX)
         nx = grid%nx_total; ny = grid%ny_total
         ms_a%nz_ml = 2; ms_b%nz_ml = 2
         call ms_a%init(grid); call ms_b%init(grid)
         call pgf_a%init(grid, nz_ml=2); call pgf_b%init(grid, nz_ml=2)
         pgf_a%variant = OPGF_VARIANT_FV_MOM6
         pgf_b%variant = OPGF_VARIANT_FV_MOM6
         pgf_a%rho0 = 1035.0_wp; pgf_a%rho_ref = 1035.0_wp
         pgf_b%rho0 = 1035.0_wp; pgf_b%rho_ref = 1035.0_wp
         pgf_a%gfs_scale = 1.0_wp
         pgf_b%gfs_scale = 0.1_wp

         ! Stratified: bed layer heavier; bath flat; non-trivial h to
         ! get a baroclinic (interface-tilt) PGF signal.
         allocate (b_field(nx, ny), source=B0)
         do j = 1, ny
            do i = 1, nx
               ms_a%h_layer(i, j, 1) = 1000.0_wp + 5.0_wp*real(i, wp)   ! bed thicker east
               ms_a%h_layer(i, j, 2) = B0 - ms_a%h_layer(i, j, 1)        ! surf thinner east
               ms_b%h_layer(i, j, :) = ms_a%h_layer(i, j, :)
            end do
         end do
         ms_a%rho_layer(:, :, 1) = 1036.0_wp; ms_a%rho_layer(:, :, 2) = 1035.0_wp
         ms_b%rho_layer = ms_a%rho_layer

         call run_pgf_with_b(ms_a, pgf_a, b_field, DX)
         call run_pgf_with_b(ms_b, pgf_b, b_field, DX)

         ! Baroclinic part = PFu(bed) − PFu(surf).  Must agree.
         max_bc_diff = 0.0_wp
         do j = NGHOST + 1, NGHOST + 4
            do i = NGHOST + 2, NGHOST + 14
               bc_a = pgf_a%dpdx_face%data(i, j, 1) - pgf_a%dpdx_face%data(i, j, 2)
               bc_b = pgf_b%dpdx_face%data(i, j, 1) - pgf_b%dpdx_face%data(i, j, 2)
               max_bc_diff = max(max_bc_diff, abs(bc_a - bc_b))
            end do
         end do
         call check(error, max_bc_diff < 1.0e-12_wp, &
                    "gfs_scale changed baroclinic mode (PFu(bed)−PFu(surf))")

         deallocate (b_field)
      end block checks
      call pgf_a%destroy(); call pgf_b%destroy()
      call ms_a%destroy(); call ms_b%destroy()
   end subroutine test_gfs_scale_baroclinic

end module test_ocean_pgf_gfs_scale
