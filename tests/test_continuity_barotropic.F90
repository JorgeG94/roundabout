!! Unit tests for the C-grid barotropic continuity step
!! (rdb_continuity Phase 2b — PPM face reconstruction with the
!! Colella-Woodward monotonic limiter, plus the per-face mass-flux
!! + flux-divergence + forward-Euler update pipeline).  See
!! `docs/ROADMAP_OCEAN.md` Phase 2.
!!
!! Cases:
!!   * Uniform lake at rest — uniform h, zero velocity, h preserved
!!     bit-for-bit (the discrete equivalent of `dh/dt = 0`).
!!   * Non-uniform lake at rest — spatially varying h, zero velocity,
!!     h preserved.  Catches reconstruction bugs that uniform-h
!!     tests miss (a stride error reading the wrong face still
!!     passes uniform-h, but breaks here).
!!   * Gaussian hump propagation — uniform u, hump initialised away
!!     from walls; after a multi-step advection the peak amplitude
!!     must be preserved within 5% of its initial value and the
!!     peak location must have shifted by u*N*dt within ±1 cell.
!!     This is the test that distinguishes PPM from 1st-order
!!     upwind (which would diffuse the peak by ~30-50% over the
!!     same advection distance).
module test_continuity_barotropic
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_barotropic_state, only: barotropic_state_t
   use rdb_continuity, only: continuity_t, &
                             continuity_compute_fluxes_barotropic, &
                             continuity_apply_fluxes_barotropic
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   implicit none
   private

   public :: collect_continuity_barotropic_tests

   integer, parameter :: NX_PHYS = 8
   integer, parameter :: NY_PHYS = 6
   integer, parameter :: NGHOST = 2

contains

   subroutine collect_continuity_barotropic_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("uniform_lake_at_rest", test_uniform_lake_at_rest), &
                  new_unittest("nonuniform_lake_at_rest", test_nonuniform_lake_at_rest), &
                  new_unittest("gaussian_hump_advection", test_gaussian_hump), &
                  new_unittest("closed_basin_mass_conservation", test_closed_basin_mass) &
                  ]
   end subroutine collect_continuity_barotropic_tests

   subroutine make_grid(grid, nx_phys, ny_phys, dx, dy)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      real(wp), intent(in) :: dx, dy
      call grid%init(nx_phys, ny_phys, NGHOST, dx, dy)
   end subroutine make_grid

   subroutine map_in(grid, bs, ct, metrics)
      type(hgrid_t), intent(in) :: grid
      type(barotropic_state_t), intent(inout) :: bs
      type(continuity_t), intent(inout) :: ct
      type(ocean_metrics_t), intent(inout) :: metrics
      call make_cartesian_metrics(metrics, grid)
      !$acc enter data copyin(bs)
      call bs%enter_data()
      !$acc enter data copyin(ct)
      call ct%enter_data()
   end subroutine map_in

   subroutine map_out(bs, ct, metrics)
      type(barotropic_state_t), intent(inout) :: bs
      type(continuity_t), intent(inout) :: ct
      type(ocean_metrics_t), intent(inout) :: metrics
      call ct%exit_data()
      !$acc exit data delete(ct)
      call bs%exit_data()
      !$acc exit data delete(bs)
      call destroy_cartesian_metrics(metrics)
   end subroutine map_out

   subroutine run_step(grid, metrics, ct, bs, dt)
      !! One continuity step on the device.  Assumes the caller
      !! already mapped bs + ct.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(continuity_t), intent(inout) :: ct
      type(barotropic_state_t), intent(inout) :: bs
      real(wp), intent(in) :: dt
      call continuity_compute_fluxes_barotropic(grid, metrics, ct, bs)
      call continuity_apply_fluxes_barotropic(bs, dt)
   end subroutine run_step

   ! -----------------------------------------------------------------
   ! Cases
   ! -----------------------------------------------------------------

   subroutine test_uniform_lake_at_rest(error)
      !! Uniform h, zero velocity -> the discrete continuity step
      !! must reproduce h bit-for-bit.  Constancy preservation —
      !! the most basic property of any conservative scheme.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(barotropic_state_t) :: bs
      type(continuity_t) :: ct
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: H0 = 10.0_wp
      real(wp), parameter :: DT = 0.1_wp
      real(wp) :: max_diff

      call make_grid(grid, NX_PHYS, NY_PHYS, 1.0_wp, 1.0_wp)
      call bs%init(grid)
      call ct%init(grid)
      bs%h = H0
      bs%u_face_x = 0.0_wp
      bs%v_face_y = 0.0_wp

      call map_in(grid, bs, ct, metrics)
      call run_step(grid, metrics, ct, bs, DT)
      call map_out(bs, ct, metrics)

      max_diff = maxval(abs(bs%h - H0))
      call check(error, max_diff < 1.0e-14_wp, &
                 "uniform lake-at-rest: h drift exceeded round-off")

      call ct%destroy()
      call bs%destroy()
   end subroutine test_uniform_lake_at_rest

   subroutine test_nonuniform_lake_at_rest(error)
      !! Spatially varying h with zero velocity must still be
      !! preserved bit-for-bit.  Catches bugs that uniform-h tests
      !! miss — e.g. a kernel that reads h(i+1, j) instead of
      !! h(i-1, j) for the west face under positive u.  Also
      !! exercises the PPM-vs-1st-order-fallback boundary cells
      !! (h(2, :) and h(nx-1, :) etc.) since they have different
      !! initial values.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(barotropic_state_t) :: bs
      type(continuity_t) :: ct
      type(ocean_metrics_t) :: metrics
      real(wp), allocatable :: h_initial(:, :)
      real(wp) :: max_diff
      integer :: i, j, nx, ny

      call make_grid(grid, NX_PHYS, NY_PHYS, 1.0_wp, 1.0_wp)
      call bs%init(grid)
      call ct%init(grid)
      nx = grid%nx_total
      ny = grid%ny_total

      allocate (h_initial(nx, ny))
      do j = 1, ny
         do i = 1, nx
            h_initial(i, j) = 10.0_wp + 0.5_wp*real(i, wp) + 0.3_wp*real(j, wp)
         end do
      end do
      bs%h = h_initial
      bs%u_face_x = 0.0_wp
      bs%v_face_y = 0.0_wp

      call map_in(grid, bs, ct, metrics)
      call run_step(grid, metrics, ct, bs, 0.1_wp)
      call map_out(bs, ct, metrics)

      max_diff = maxval(abs(bs%h - h_initial))
      call check(error, max_diff < 1.0e-14_wp, &
                 "non-uniform lake-at-rest: h drift exceeded round-off")

      deallocate (h_initial)
      call ct%destroy()
      call bs%destroy()
   end subroutine test_nonuniform_lake_at_rest

   subroutine test_gaussian_hump(error)
      !! Advect a Gaussian hump rightward by ~its own width over N
      !! steps, then check (a) the peak amplitude is preserved
      !! within 5% (PPM has minimal diffusion; 1st-order upwind
      !! would lose 30-50%), and (b) the peak has moved by
      !! u*N*dt cells within ±1.
      !!
      !! The kernel uses closed walls + uniform interior u, so the
      !! right wall accumulates mass and the left depletes (the
      !! domain is a "pipe with closed ends, fluid flowing inside").
      !! That accumulation is not what we're testing.  Both peak
      !! amplitude and peak location are measured in an interior
      !! window centred on the expected final-peak location, well
      !! clear of the wall pile-up + drain zones.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(barotropic_state_t) :: bs
      type(continuity_t) :: ct
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: H_BASE = 10.0_wp
      real(wp), parameter :: AMP = 2.0_wp
      real(wp), parameter :: U_CONST = 1.0_wp
      real(wp), parameter :: DT = 0.5_wp
      real(wp), parameter :: DX = 1.0_wp
      real(wp), parameter :: SIGMA = 4.0_wp
      integer, parameter :: NX_HUMP = 64
      integer, parameter :: NY_HUMP = 4
      integer, parameter :: N_STEPS = 8
      integer, parameter :: WINDOW_LO = 8
      integer, parameter :: WINDOW_HI = 36
      real(wp) :: x0, x_expected, dx_expected
      real(wp) :: peak_initial, peak_final, peak_diff
      integer :: i, j, nx, ny, step, j_probe, i_peak_initial, i_peak_final

      call make_grid(grid, NX_HUMP, NY_HUMP, DX, 1.0_wp)
      call bs%init(grid)
      call ct%init(grid)
      nx = grid%nx_total
      ny = grid%ny_total
      j_probe = ny/2

      ! Initial hump centred at x0 = nx/4 (well away from both walls)
      x0 = real(nx, wp)*0.25_wp
      do j = 1, ny
         do i = 1, nx
            bs%h(i, j) = H_BASE + AMP*exp(-((real(i, wp) - x0)/SIGMA)**2)
         end do
      end do
      bs%u_face_x = U_CONST
      bs%v_face_y = 0.0_wp

      peak_initial = maxval(bs%h(WINDOW_LO:WINDOW_HI, j_probe)) - H_BASE
      i_peak_initial = WINDOW_LO + maxloc(bs%h(WINDOW_LO:WINDOW_HI, j_probe), dim=1) - 1

      call map_in(grid, bs, ct, metrics)
      do step = 1, N_STEPS
         call run_step(grid, metrics, ct, bs, DT)
      end do
      call map_out(bs, ct, metrics)

      ! Final peak in the interior window only (excludes wall pile-up
      ! at cell nx and wall drain at cell 1).
      peak_final = maxval(bs%h(WINDOW_LO:WINDOW_HI, j_probe)) - H_BASE
      i_peak_final = WINDOW_LO + maxloc(bs%h(WINDOW_LO:WINDOW_HI, j_probe), dim=1) - 1

      ! Peak amplitude check.  PPM should preserve > 95%; upwind
      ! would drop to ~70% over this advection distance.
      peak_diff = abs(peak_final - peak_initial)/peak_initial
      call check(error, peak_diff < 0.05_wp, &
                 "Gaussian hump peak amplitude diffused more than 5%")
      if (allocated(error)) then
         call ct%destroy(); call bs%destroy(); return
      end if

      ! Peak location check.
      dx_expected = U_CONST*real(N_STEPS, wp)*DT/DX
      x_expected = real(i_peak_initial, wp) + dx_expected
      call check(error, abs(real(i_peak_final, wp) - x_expected) <= 1.0_wp, &
                 "Gaussian hump peak location off by > 1 cell from u*t prediction")

      call ct%destroy()
      call bs%destroy()
   end subroutine test_gaussian_hump

   subroutine test_closed_basin_mass(error)
      !! 100-step closed-basin run with a non-trivial initial h
      !! field and a smooth velocity field that vanishes at the
      !! walls (so the flow stays bounded over the integration
      !! horizon).  Total mass = sum(h)*dx*dy is invariant in exact
      !! arithmetic; in floating-point the drift should sit at
      !! round-off (~N*eps per step, ~100*N*eps over the run).
      !! This is the strongest correctness test for the continuity
      !! kernel — any sign error, stencil shift, wall-handling bug,
      !! or PPM reconstruction asymmetry breaks it.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(barotropic_state_t) :: bs
      type(continuity_t) :: ct
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: H_BASE = 10.0_wp
      real(wp), parameter :: H_AMP = 0.5_wp
      real(wp), parameter :: U_AMP = 0.3_wp
      real(wp), parameter :: DT = 0.1_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      integer, parameter :: N_STEPS = 100
      integer, parameter :: NX_BASIN = 32
      integer, parameter :: NY_BASIN = 16
      real(wp) :: total_initial, total_final, drift
      real(wp) :: h_min, h_max
      integer :: i, j, nx, ny, step

      call make_grid(grid, NX_BASIN, NY_BASIN, 1.0_wp, 1.0_wp)
      call bs%init(grid)
      call ct%init(grid)
      nx = grid%nx_total
      ny = grid%ny_total

      ! Non-trivial h field
      do j = 1, ny
         do i = 1, nx
            bs%h(i, j) = H_BASE + H_AMP* &
                         sin(2.0_wp*PI*real(i, wp)/real(nx, wp))* &
                         cos(2.0_wp*PI*real(j, wp)/real(ny, wp))
         end do
      end do
      ! u_face_x: zero at the walls (i=1, i=nx+1) by construction
      ! so the closed-wall override is consistent with the IC.
      do j = 1, ny
         do i = 1, nx + 1
            bs%u_face_x(i, j) = U_AMP* &
                                sin(PI*real(i - 1, wp)/real(nx, wp))* &
                                sin(2.0_wp*PI*real(j - 1, wp)/real(ny - 1, wp))
         end do
      end do
      ! v_face_y: mirror — zero at j=1 and j=ny+1
      do j = 1, ny + 1
         do i = 1, nx
            bs%v_face_y(i, j) = U_AMP* &
                                sin(PI*real(j - 1, wp)/real(ny, wp))* &
                                sin(2.0_wp*PI*real(i - 1, wp)/real(nx - 1, wp))
         end do
      end do

      total_initial = sum(bs%h)*grid%dx*grid%dy

      call map_in(grid, bs, ct, metrics)
      do step = 1, N_STEPS
         call run_step(grid, metrics, ct, bs, DT)
      end do
      call map_out(bs, ct, metrics)

      total_final = sum(bs%h)*grid%dx*grid%dy
      drift = abs(total_final - total_initial)/abs(total_initial)
      h_min = minval(bs%h)
      h_max = maxval(bs%h)

      ! Mass-conservation floor.  For nx*ny ~ 500 and 100 steps with
      ! ~eps per FP op, expected accumulated relative drift is
      ! O(100*nx*ny*eps) ~ 1e-11.  1e-10 leaves a comfortable margin
      ! over a tight floor.
      call check(error, drift < 1.0e-10_wp, &
                 "total mass drift exceeded 1e-10 round-off floor")
      if (allocated(error)) then
         call ct%destroy(); call bs%destroy(); return
      end if

      ! Stability sanity: h must stay positive and bounded (no NaN).
      call check(error, h_min > 0.0_wp, &
                 "h went negative — kernel unstable")
      if (allocated(error)) then
         call ct%destroy(); call bs%destroy(); return
      end if
      call check(error, h_max < 100.0_wp*H_BASE, &
                 "h grew > 100x H_BASE — kernel unstable")

      call ct%destroy()
      call bs%destroy()
   end subroutine test_closed_basin_mass

end module test_continuity_barotropic
