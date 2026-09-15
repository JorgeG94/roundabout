!! Unit tests for the multilayer C-grid continuity-PPM kernel
!! (rdb_continuity Phase 5a — per-layer lift of the Phase 2
!! barotropic kernel).  Each k-slice runs the same PPM
!! reconstruction + upwind flux + flux-divergence pipeline as the
!! barotropic kernel, parallelised over (k, j, i) in the
!! do-concurrent loops.
!!
!! Cases:
!!   * Per-layer lake-at-rest — uniform h_layer per layer (different
!!     thickness in each), zero velocity: every layer's thickness
!!     must be preserved bit-for-bit, with no cross-layer leakage.
!!   * Stratified Gaussian hump — initialise a Gaussian in one layer
!!     only, with uniform u in that layer only, zero everything in
!!     the other layers.  After advection: the active layer's hump
!!     propagates with > 95% peak preservation (same as barotropic
!!     PPM); inactive layers stay unchanged bit-for-bit.  Probes
!!     that the per-layer kernel doesn't smear vertically.
!!   * Per-layer mass conservation — non-trivial IC in every layer
!!     with smooth wall-vanishing velocities; 100 outer steps;
!!     mass-per-layer drift < 1e-10 (round-off floor).
module test_continuity_multilayer
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_continuity, only: continuity_t, &
                             continuity_compute_fluxes, &
                             continuity_apply_fluxes, &
                             continuity_step_split, &
                             continuity_zonal_flux
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   implicit none
   private

   public :: collect_continuity_multilayer_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 3

contains

   subroutine collect_continuity_multilayer_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("per_layer_lake_at_rest", test_per_layer_lake_at_rest), &
                  new_unittest("active_layer_gaussian_hump", test_active_layer_hump), &
                  new_unittest("per_layer_mass_conservation", test_per_layer_mass), &
                  new_unittest("split_uniform_h_preserved", test_split_uniform), &
                  new_unittest("split_mass_conservation", test_split_mass), &
                  new_unittest("split_zonal_only_matches_unsplit", &
                               test_split_zonal_only), &
                  new_unittest("volcfl_off_bit_identical", &
                               test_volcfl_off_bit_identical), &
                  new_unittest("volcfl_deep_shallow_quiescent", &
                               test_volcfl_deep_shallow), &
                  new_unittest("volcfl_small_cfl_reduces_to_edge", &
                               test_volcfl_cfl_reduction), &
                  new_unittest("renorm_visc_rem_transport_matches", &
                               test_renorm_vr_transport), &
                  new_unittest("renorm_visc_rem_gamma_shares", &
                               test_renorm_vr_shares) &
                  ]
   end subroutine collect_continuity_multilayer_tests

   subroutine make_grid(grid, nx_phys, ny_phys, dx, dy)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      real(wp), intent(in) :: dx, dy
      call grid%init(nx_phys, ny_phys, NGHOST, dx, dy)
   end subroutine make_grid

   subroutine map_in(ms, ct)
      type(multilayer_state_t), intent(inout) :: ms
      type(continuity_t), intent(inout) :: ct
      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(ct)
      call ct%enter_data()
   end subroutine map_in

   subroutine map_out(ms, ct)
      type(multilayer_state_t), intent(inout) :: ms
      type(continuity_t), intent(inout) :: ct
      ! Pull the computed mass fluxes D->H before the device copy is deleted.
      ! exit_data deletes the workspace flux arrays (they are not prognostic
      ! copyout fields), so on -stdpar=gpu the host would otherwise keep its
      ! stale (zero) values — fine on gfortran (host==device), but the
      ! flux-reading vol_cfl subtests need the real device result.
      !$acc update self(ms%mass_flux_x_layer, ms%mass_flux_y_layer)
      call ct%exit_data()
      !$acc exit data delete(ct)
      call ms%exit_data()
      !$acc exit data delete(ms)
   end subroutine map_out

   subroutine run_step(grid, metrics, ct, ms, dt)
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(continuity_t), intent(inout) :: ct
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: dt
      call continuity_compute_fluxes(grid, metrics, ct, ms)
      call continuity_apply_fluxes(ms, dt)
   end subroutine run_step

   ! -----------------------------------------------------------------
   ! Cases
   ! -----------------------------------------------------------------

   subroutine test_per_layer_lake_at_rest(error)
      !! Different thickness in each layer, zero velocity everywhere.
      !! Per-layer thickness must be preserved bit-for-bit with no
      !! cross-layer leakage.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: DT = 0.1_wp
      real(wp) :: h_targets(NZ)
      real(wp) :: max_diff
      integer :: k
      checks: block

         call make_grid(grid, 8, 6, 1.0_wp, 1.0_wp)
         call make_cartesian_metrics(metrics, grid)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ct%init(grid, nz_ml=NZ)
         h_targets = [2.0_wp, 5.0_wp, 3.0_wp]
         do k = 1, NZ
            ms%h_layer(:, :, k) = h_targets(k)
         end do
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp

         call map_in(ms, ct)
         call run_step(grid, metrics, ct, ms, DT)
         call map_out(ms, ct)

         do k = 1, NZ
            max_diff = maxval(abs(ms%h_layer(:, :, k) - h_targets(k)))
            call check(error, max_diff < 1.0e-14_wp, &
                       "per-layer lake-at-rest: thickness drifted")
            if (allocated(error)) exit checks
         end do

      end block checks
      call ct%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_per_layer_lake_at_rest

   subroutine test_active_layer_hump(error)
      !! Gaussian hump in ONE layer with non-zero u; other layers
      !! quiescent.  After multi-step advection: active layer's hump
      !! has propagated with > 95% peak preservation (same physics
      !! the barotropic test demonstrates); inactive layers stay at
      !! their initial uniform thickness bit-for-bit.  Catches
      !! cross-layer index bugs (e.g. a stride that mixes layers
      !! during reconstruction).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: H_BASE = 10.0_wp
      real(wp), parameter :: AMP = 2.0_wp
      real(wp), parameter :: U_CONST = 1.0_wp
      real(wp), parameter :: DT = 0.5_wp
      real(wp), parameter :: SIGMA = 4.0_wp
      integer, parameter :: NX_HUMP = 64
      integer, parameter :: NY_HUMP = 4
      integer, parameter :: N_STEPS = 8
      integer, parameter :: ACTIVE = 2     ! middle layer
      integer, parameter :: WINDOW_LO = 8
      integer, parameter :: WINDOW_HI = 36
      real(wp) :: x0, peak_initial, peak_final, peak_diff, max_inactive
      integer :: i, j, k, nx, ny, step, j_probe
      checks: block

         call make_grid(grid, NX_HUMP, NY_HUMP, 1.0_wp, 1.0_wp)
         call make_cartesian_metrics(metrics, grid)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ct%init(grid, nz_ml=NZ)
         nx = grid%nx_total
         ny = grid%ny_total
         j_probe = ny/2

         ! Quiescent layers everywhere; active layer carries the hump
         ms%h_layer = H_BASE
         x0 = real(nx, wp)*0.25_wp
         do j = 1, ny
            do i = 1, nx
               ms%h_layer(i, j, ACTIVE) = H_BASE + &
                                          AMP*exp(-((real(i, wp) - x0)/SIGMA)**2)
            end do
         end do
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         ms%u_face_x_layer(:, :, ACTIVE) = U_CONST

         peak_initial = maxval(ms%h_layer(WINDOW_LO:WINDOW_HI, j_probe, ACTIVE)) - H_BASE

         call map_in(ms, ct)
         do step = 1, N_STEPS
            call run_step(grid, metrics, ct, ms, DT)
         end do
         call map_out(ms, ct)

         peak_final = maxval(ms%h_layer(WINDOW_LO:WINDOW_HI, j_probe, ACTIVE)) - H_BASE
         peak_diff = abs(peak_final - peak_initial)/peak_initial
         call check(error, peak_diff < 0.05_wp, &
                    "active-layer hump peak diffused > 5%")
         if (allocated(error)) exit checks

         ! Inactive layers must be bit-for-bit unchanged
         do k = 1, NZ
            if (k == ACTIVE) cycle
            max_inactive = maxval(abs(ms%h_layer(:, :, k) - H_BASE))
            call check(error, max_inactive < 1.0e-14_wp, &
                       "inactive layer drifted — cross-layer leakage")
            if (allocated(error)) exit checks
         end do

      end block checks
      call ct%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_active_layer_hump

   subroutine test_per_layer_mass(error)
      !! 100-step closed-basin run on every layer with smooth
      !! wall-vanishing velocity fields; per-layer mass drift must
      !! sit at the 1e-10 round-off floor for all layers.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: H_BASE = 10.0_wp
      real(wp), parameter :: H_AMP = 0.5_wp
      real(wp), parameter :: U_AMP = 0.3_wp
      real(wp), parameter :: DT = 0.05_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      integer, parameter :: N_STEPS = 100
      integer, parameter :: NX_BASIN = 32
      integer, parameter :: NY_BASIN = 16
      real(wp) :: total_initial(NZ), total_final(NZ), drift, layer_phase
      real(wp) :: h_min, h_max
      integer :: i, j, k, nx, ny, step
      checks: block

         call make_grid(grid, NX_BASIN, NY_BASIN, 1.0_wp, 1.0_wp)
         call make_cartesian_metrics(metrics, grid)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ct%init(grid, nz_ml=NZ)
         nx = grid%nx_total
         ny = grid%ny_total

         ! Per-layer IC with a layer-dependent phase so each layer has
         ! a distinct pattern (catches cross-layer mass leaks even
         ! within the closed-basin total).
         do k = 1, NZ
            layer_phase = real(k - 1, wp)*PI/real(NZ, wp)
            do j = 1, ny
               do i = 1, nx
                  ms%h_layer(i, j, k) = H_BASE + H_AMP* &
                                        sin(2.0_wp*PI*real(i, wp)/real(nx, wp) + layer_phase)* &
                                        cos(2.0_wp*PI*real(j, wp)/real(ny, wp))
               end do
            end do
            do j = 1, ny
               do i = 1, nx + 1
                  ms%u_face_x_layer(i, j, k) = U_AMP* &
                                               sin(PI*real(i - 1, wp)/real(nx, wp))* &
                                               sin(2.0_wp*PI*real(j - 1, wp)/real(ny - 1, wp) + layer_phase)
               end do
            end do
            do j = 1, ny + 1
               do i = 1, nx
                  ms%v_face_y_layer(i, j, k) = U_AMP* &
                                               sin(PI*real(j - 1, wp)/real(ny, wp))* &
                                               sin(2.0_wp*PI*real(i - 1, wp)/real(nx - 1, wp) + layer_phase)
               end do
            end do
         end do

         do k = 1, NZ
            total_initial(k) = sum(ms%h_layer(:, :, k))*grid%dx*grid%dy
         end do

         call map_in(ms, ct)
         do step = 1, N_STEPS
            call run_step(grid, metrics, ct, ms, DT)
         end do
         call map_out(ms, ct)

         h_min = minval(ms%h_layer)
         h_max = maxval(ms%h_layer)

         do k = 1, NZ
            total_final(k) = sum(ms%h_layer(:, :, k))*grid%dx*grid%dy
            drift = abs(total_final(k) - total_initial(k))/abs(total_initial(k))
            call check(error, drift < 1.0e-10_wp, &
                       "per-layer mass drift exceeded 1e-10 round-off floor")
            if (allocated(error)) exit checks
         end do
         call check(error, h_min > 0.0_wp, &
                    "h_layer went negative — kernel unstable")
         if (allocated(error)) exit checks
         call check(error, h_max < 100.0_wp*H_BASE, &
                    "h_layer grew > 100x H_BASE — kernel unstable")

      end block checks
      call ct%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_per_layer_mass

   ! -----------------------------------------------------------------
   ! Directionally-split (Lie) form
   ! -----------------------------------------------------------------

   subroutine test_split_uniform(error)
      !! Constancy preservation under the split form.  Uniform h
      !! and uniform face velocities → ∂(uh)/∂x = 0 in the deep
      !! interior, so both substeps leave `h_layer` at H0 bit-for-
      !! bit.  Probed only on cells far from the closed walls
      !! (where boundary-induced flux differences are expected
      !! and equal in both forms).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: H0 = 8.0_wp
      real(wp), parameter :: U0 = 0.4_wp
      real(wp), parameter :: V0 = -0.2_wp
      real(wp), parameter :: DT = 0.05_wp
      real(wp) :: max_dev
      integer :: ig, nx, ny

      call make_grid(grid, 16, 12, 1.0_wp, 1.0_wp)
      call make_cartesian_metrics(metrics, grid)
      ms%nz_ml = NZ
      call ms%init(grid)
      call ct%init(grid, nz_ml=NZ)
      nx = grid%nx_total
      ny = grid%ny_total
      ig = grid%nghost

      ms%h_layer = H0
      ms%u_face_x_layer = U0
      ms%v_face_y_layer = V0

      call map_in(ms, ct)
      call continuity_step_split(grid, metrics, ct, ms, DT)
      call map_out(ms, ct)

      ! Cells > 3 in from each wall see only PPM-interior fluxes
      ! that telescope to zero for uniform h, regardless of u, v.
      max_dev = maxval(abs(ms%h_layer(ig + 4:nx - ig - 3, &
                                      ig + 4:ny - ig - 3, :) - H0))
      call check(error, max_dev < 1.0e-13_wp, &
                 "split form did not preserve uniform h in deep interior")

      call ct%destroy(); call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_split_uniform

   subroutine test_split_mass(error)
      !! 100-step closed-basin run under the directionally-split
      !! form with the same IC pattern as `test_per_layer_mass`.
      !! Total per-layer mass drift must sit at the 1e-10 round-
      !! off floor — the Lie split is conservative by construction
      !! (each substep is a face-flux update on a closed-wall
      !! domain).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: H_BASE = 10.0_wp
      real(wp), parameter :: H_AMP = 0.5_wp
      real(wp), parameter :: U_AMP = 0.3_wp
      real(wp), parameter :: DT = 0.05_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      integer, parameter :: N_STEPS = 100
      integer, parameter :: NX_BASIN = 32
      integer, parameter :: NY_BASIN = 16
      real(wp) :: total_initial(NZ), total_final(NZ), drift, layer_phase
      real(wp) :: h_min
      integer :: i, j, k, nx, ny, step
      checks: block

         call make_grid(grid, NX_BASIN, NY_BASIN, 1.0_wp, 1.0_wp)
         call make_cartesian_metrics(metrics, grid)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ct%init(grid, nz_ml=NZ)
         nx = grid%nx_total
         ny = grid%ny_total

         do k = 1, NZ
            layer_phase = real(k - 1, wp)*PI/real(NZ, wp)
            do j = 1, ny
               do i = 1, nx
                  ms%h_layer(i, j, k) = H_BASE + H_AMP* &
                                        sin(2.0_wp*PI*real(i, wp)/real(nx, wp) + layer_phase)* &
                                        cos(2.0_wp*PI*real(j, wp)/real(ny, wp))
               end do
            end do
            do j = 1, ny
               do i = 1, nx + 1
                  ms%u_face_x_layer(i, j, k) = U_AMP* &
                                               sin(PI*real(i - 1, wp)/real(nx, wp))* &
                                               sin(2.0_wp*PI*real(j - 1, wp)/real(ny - 1, wp) + layer_phase)
               end do
            end do
            do j = 1, ny + 1
               do i = 1, nx
                  ms%v_face_y_layer(i, j, k) = U_AMP* &
                                               sin(PI*real(j - 1, wp)/real(ny, wp))* &
                                               sin(2.0_wp*PI*real(i - 1, wp)/real(nx - 1, wp) + layer_phase)
               end do
            end do
         end do

         do k = 1, NZ
            total_initial(k) = sum(ms%h_layer(:, :, k))*grid%dx*grid%dy
         end do

         call map_in(ms, ct)
         do step = 1, N_STEPS
            call continuity_step_split(grid, metrics, ct, ms, DT)
         end do
         call map_out(ms, ct)

         h_min = minval(ms%h_layer)
         do k = 1, NZ
            total_final(k) = sum(ms%h_layer(:, :, k))*grid%dx*grid%dy
            drift = abs(total_final(k) - total_initial(k))/abs(total_initial(k))
            call check(error, drift < 1.0e-10_wp, &
                       "split-form per-layer mass drift exceeded 1e-10")
            if (allocated(error)) exit checks
         end do
         call check(error, h_min > 0.0_wp, &
                    "split form drove h negative — kernel unstable")

      end block checks
      call ct%destroy(); call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_split_mass

   subroutine test_split_zonal_only(error)
      !! With v_face_y_layer = 0, the meridional substep is a no-op
      !! (mass_flux_y = 0 everywhere → div_y = 0 → h unchanged).
      !! In that degenerate case the split form must produce
      !! bit-identical `h_layer` to the unsplit form, which under
      !! v=0 also only sees x-divergence.  A sanity check that the
      !! lift is structural — no spurious touches to h in either
      !! substep.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms_split, ms_unsplit
      type(continuity_t) :: ct_split, ct_unsplit
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: H_BASE = 10.0_wp
      real(wp), parameter :: U_AMP = 0.3_wp
      real(wp), parameter :: DT = 0.05_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp) :: max_diff
      integer :: i, j, k, nx, ny

      call make_grid(grid, 16, 12, 1.0_wp, 1.0_wp)
      call make_cartesian_metrics(metrics, grid)
      ms_split%nz_ml = NZ
      call ms_split%init(grid)
      ms_unsplit%nz_ml = NZ
      call ms_unsplit%init(grid)
      call ct_split%init(grid, nz_ml=NZ)
      call ct_unsplit%init(grid, nz_ml=NZ)
      nx = grid%nx_total
      ny = grid%ny_total

      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx
               ms_split%h_layer(i, j, k) = H_BASE + &
                                           sin(2.0_wp*PI*real(i, wp)/real(nx, wp))
               ms_unsplit%h_layer(i, j, k) = ms_split%h_layer(i, j, k)
            end do
            do i = 1, nx + 1
               ms_split%u_face_x_layer(i, j, k) = U_AMP* &
                                                  sin(PI*real(i - 1, wp)/real(nx, wp))
               ms_unsplit%u_face_x_layer(i, j, k) = ms_split%u_face_x_layer(i, j, k)
            end do
         end do
      end do
      ms_split%v_face_y_layer = 0.0_wp
      ms_unsplit%v_face_y_layer = 0.0_wp

      call map_in(ms_split, ct_split)
      call continuity_step_split(grid, metrics, ct_split, ms_split, DT)
      call map_out(ms_split, ct_split)

      call map_in(ms_unsplit, ct_unsplit)
      call continuity_compute_fluxes(grid, metrics, ct_unsplit, ms_unsplit)
      call continuity_apply_fluxes(ms_unsplit, DT)
      call map_out(ms_unsplit, ct_unsplit)

      max_diff = maxval(abs(ms_split%h_layer - ms_unsplit%h_layer))
      call check(error, max_diff < 1.0e-13_wp, &
                 "split form with v=0 did not match unsplit form bitwise")

      call ct_split%destroy(); call ct_unsplit%destroy()
      call ms_split%destroy(); call ms_unsplit%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_split_zonal_only

   ! -----------------------------------------------------------------
   ! MOM6 swept-volume (vol_CFL) continuity face thickness
   ! -----------------------------------------------------------------

   subroutine test_volcfl_off_bit_identical(error)
      !! Default-off bit-identity gate.  Run the full split step once
      !! with the default continuity_t (vol_cfl = .false.) and once
      !! with vol_cfl set EXPLICITLY false — the two must be byte-
      !! identical.  Together with the existing lake-at-rest /
      !! Gaussian-hump / mass-conservation cases (all run with the
      !! default vol_cfl = .false.) this guards that the new branch is
      !! never entered when off, so the edge-value flux is unchanged.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms_a, ms_b
      type(continuity_t) :: ct_a, ct_b
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: H_BASE = 10.0_wp
      real(wp), parameter :: U_AMP = 0.3_wp
      real(wp), parameter :: DT = 0.05_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp) :: max_diff
      integer :: i, j, k, nx, ny

      call make_grid(grid, 16, 12, 1.0_wp, 1.0_wp)
      call make_cartesian_metrics(metrics, grid)
      ms_a%nz_ml = NZ; call ms_a%init(grid)
      ms_b%nz_ml = NZ; call ms_b%init(grid)
      call ct_a%init(grid, nz_ml=NZ)
      call ct_b%init(grid, nz_ml=NZ)
      ct_b%vol_cfl = .false.   ! explicit default
      nx = grid%nx_total
      ny = grid%ny_total

      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx
               ms_a%h_layer(i, j, k) = H_BASE + sin(2.0_wp*PI*real(i, wp)/real(nx, wp))
               ms_b%h_layer(i, j, k) = ms_a%h_layer(i, j, k)
            end do
            do i = 1, nx + 1
               ms_a%u_face_x_layer(i, j, k) = U_AMP*sin(PI*real(i - 1, wp)/real(nx, wp))
               ms_b%u_face_x_layer(i, j, k) = ms_a%u_face_x_layer(i, j, k)
            end do
         end do
      end do
      ms_a%v_face_y_layer = 0.3_wp*U_AMP
      ms_b%v_face_y_layer = ms_a%v_face_y_layer

      call map_in(ms_a, ct_a)
      call continuity_step_split(grid, metrics, ct_a, ms_a, DT)
      call map_out(ms_a, ct_a)

      call map_in(ms_b, ct_b)
      call continuity_step_split(grid, metrics, ct_b, ms_b, DT)
      call map_out(ms_b, ct_b)

      max_diff = maxval(abs(ms_a%h_layer - ms_b%h_layer))
      call check(error, max_diff == 0.0_wp, &
                 "vol_cfl=.false. must be byte-identical to the default")

      call ct_a%destroy(); call ct_b%destroy()
      call ms_a%destroy(); call ms_b%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_volcfl_off_bit_identical

   subroutine test_volcfl_deep_shallow(error)
      !! THE regression for the shelf-break near-bed residual.  A
      !! steep 100 m <-> 5000 m thickness step (uniform per column,
      !! quiescent flow) must stay exactly at rest under vol_cfl —
      !! zero velocity gives zero CFL and zero mass flux at EVERY
      !! face regardless of the curvature term — and conserve mass to
      !! the round-off floor.  Flat-bath tests are blind to this
      !! (project memory feedback_zstar_full_needs_deep_shallow_test);
      !! the deep-shallow step exercises the curv3 / dh terms.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: H_DEEP = 5000.0_wp
      real(wp), parameter :: H_SHALLOW = 100.0_wp
      real(wp), parameter :: DT = 50.0_wp
      integer, parameter :: N_STEPS = 50
      integer, parameter :: NX_B = 24
      integer, parameter :: NY_B = 8
      real(wp) :: total_initial, total_final, drift, max_flux, max_dev
      real(wp), allocatable :: h0(:, :, :)
      integer :: i, j, k, nx, ny, step, i_step
      checks: block

         call make_grid(grid, NX_B, NY_B, 2000.0_wp, 2000.0_wp)
         call make_cartesian_metrics(metrics, grid)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ct%init(grid, nz_ml=NZ)
         ct%vol_cfl = .true.
         nx = grid%nx_total
         ny = grid%ny_total
         i_step = nx/2

         ! Steep deep->shallow step at i_step; uniform per column,
         ! split equally across layers.  Quiescent flow.
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx
                  if (i <= i_step) then
                     ms%h_layer(i, j, k) = H_DEEP/real(NZ, wp)
                  else
                     ms%h_layer(i, j, k) = H_SHALLOW/real(NZ, wp)
                  end if
               end do
            end do
         end do
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp

         allocate (h0(nx, ny, NZ))
         h0 = ms%h_layer
         total_initial = sum(ms%h_layer)*grid%dx*grid%dy

         call map_in(ms, ct)
         do step = 1, N_STEPS
            call continuity_step_split(grid, metrics, ct, ms, DT)
         end do
         call map_out(ms, ct)

         ! Must stay at rest bit-for-bit (zero velocity => zero flux).
         max_dev = maxval(abs(ms%h_layer - h0))
         call check(error, max_dev < 1.0e-9_wp, &
                    "vol_cfl deep-shallow quiescent: thickness drifted from rest")
         if (allocated(error)) exit checks

         ! Mass conserved to round-off.
         total_final = sum(ms%h_layer)*grid%dx*grid%dy
         drift = abs(total_final - total_initial)/abs(total_initial)
         call check(error, drift < 1.0e-12_wp, &
                    "vol_cfl deep-shallow quiescent: mass not conserved")
         if (allocated(error)) exit checks

         ! And the diagnosed mass flux is identically zero.
         max_flux = maxval(abs(ms%mass_flux_x_layer)) + maxval(abs(ms%mass_flux_y_layer))
         call check(error, max_flux < 1.0e-9_wp, &
                    "vol_cfl deep-shallow quiescent: nonzero mass flux at rest")

      end block checks
      if (allocated(h0)) deallocate (h0)
      call ct%destroy(); call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_volcfl_deep_shallow

   subroutine test_volcfl_cfl_reduction(error)
      !! CFL-reduction check.  On a smoothly varying h with uniform
      !! u, the vol_cfl zonal mass flux must (a) approach the edge-
      !! value flux as CFL -> 0 (tiny dt), and (b) differ from it by a
      !! resolvable O(CFL) amount at a marginal CFL.  Compares the
      !! single-pass `continuity_zonal_flux` mass flux with vol_cfl
      !! on vs off at two dt values.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms_off, ms_on
      type(continuity_t) :: ct_off, ct_on
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: H_BASE = 100.0_wp
      real(wp), parameter :: H_AMP = 30.0_wp
      real(wp), parameter :: U0 = 1.0_wp
      real(wp), parameter :: DX = 1000.0_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp) :: dt_small, dt_marg, diff_small, diff_marg, flux_scale
      integer :: i, j, k, nx, ny

      call make_grid(grid, 24, 6, DX, DX)
      call make_cartesian_metrics(metrics, grid)
      ms_off%nz_ml = NZ; call ms_off%init(grid)
      ms_on%nz_ml = NZ; call ms_on%init(grid)
      call ct_off%init(grid, nz_ml=NZ); ct_off%vol_cfl = .false.
      call ct_on%init(grid, nz_ml=NZ); ct_on%vol_cfl = .true.
      nx = grid%nx_total
      ny = grid%ny_total

      ! Smooth h ramp with curvature; uniform positive u so the donor
      ! is unambiguous and curv3/dh are non-zero.
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx
               ms_off%h_layer(i, j, k) = H_BASE + &
                                         H_AMP*sin(2.0_wp*PI*real(i, wp)/real(nx, wp))
               ms_on%h_layer(i, j, k) = ms_off%h_layer(i, j, k)
            end do
         end do
      end do
      ms_off%u_face_x_layer = U0
      ms_on%u_face_x_layer = U0
      ms_off%v_face_y_layer = 0.0_wp
      ms_on%v_face_y_layer = 0.0_wp

      ! CFL = U0*dt/dx.  Small: dt -> CFL ~ 1e-4.  Marginal: CFL ~ 0.4.
      dt_small = 1.0e-4_wp*DX/U0
      dt_marg = 0.4_wp*DX/U0

      ! dt = 0 always uses the edge-value flux on the off path; use its
      ! magnitude as the flux scale for relative comparisons.
      call map_in(ms_off, ct_off)
      call continuity_zonal_flux(grid, metrics, ct_off, ms_off, dt_small)
      call map_out(ms_off, ct_off)
      flux_scale = maxval(abs(ms_off%mass_flux_x_layer))

      call map_in(ms_on, ct_on)
      call continuity_zonal_flux(grid, metrics, ct_on, ms_on, dt_small)
      call map_out(ms_on, ct_on)
      diff_small = maxval(abs(ms_on%mass_flux_x_layer - ms_off%mass_flux_x_layer))

      ! Re-run the on-state at the marginal CFL (zonal_flux doesn't mutate h).
      call map_in(ms_on, ct_on)
      call continuity_zonal_flux(grid, metrics, ct_on, ms_on, dt_marg)
      call map_out(ms_on, ct_on)
      diff_marg = maxval(abs(ms_on%mass_flux_x_layer - ms_off%mass_flux_x_layer))

      ! (a) small CFL: the vol_cfl correction is O(CFL)~1e-4 relative to
      ! the flux scale, i.e. it reduces to the edge-value flux.
      call check(error, diff_small < 1.0e-3_wp*flux_scale, &
                 "vol_cfl small-CFL flux did not reduce to edge-value flux")
      if (allocated(error)) return
      ! (b) marginal CFL (~0.4): a resolvable, much larger O(CFL) term —
      ! the correction grows ~ 4000x with CFL going 1e-4 -> 0.4.
      call check(error, diff_marg > 1.0e3_wp*diff_small, &
                 "vol_cfl marginal-CFL flux did not differ by the O(CFL) term")

      call ct_off%destroy(); call ct_on%destroy()
      call ms_off%destroy(); call ms_on%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_volcfl_cfl_reduction

   subroutine setup_renorm_vr_case(grid, metrics, ms, ct, uhbt, vr, u_cor, &
                                   gamma_k, u0, h0, dt)
      !! Shared fixture for the γ-weighted renormaliser tests (SPEC S2b
      !! accept): uniform h + uniform u zonal channel, a per-layer
      !! visc_rem profile γ = gamma_k, and a `uhbt` target set to HALF
      !! the unconstrained transport so `du` is well away from zero.
      !! Runs `continuity_zonal_flux(uhbt, visc_rem, u_cor)` — the same
      !! entry the driver's renorm_visc_rem dispatch takes — and pulls
      !! everything back to the host.  GPU (`mem:separate`): the local
      !! optional arrays are mapped explicitly; directives are inert on
      !! host builds.
      type(hgrid_t), intent(out) :: grid
      type(ocean_metrics_t), intent(out) :: metrics
      type(multilayer_state_t), intent(out) :: ms
      type(continuity_t), intent(out) :: ct
      real(wp), allocatable, intent(out) :: uhbt(:, :), vr(:, :, :), u_cor(:, :, :)
      real(wp), intent(in) :: gamma_k(NZ)
      real(wp), intent(in) :: u0, h0
      real(wp), intent(in) :: dt
      integer :: k, nx, ny

      call make_grid(grid, 12, 4, 1000.0_wp, 1000.0_wp)
      call make_cartesian_metrics(metrics, grid)
      ms%nz_ml = NZ; call ms%init(grid)
      call ct%init(grid, nz_ml=NZ)
      nx = grid%nx_total
      ny = grid%ny_total

      ms%h_layer = h0
      ms%u_face_x_layer = u0
      ms%v_face_y_layer = 0.0_wp

      allocate (uhbt(nx + 1, ny), vr(nx + 1, ny, NZ), u_cor(nx + 1, ny, NZ))
      ! Target: half the unconstrained transport Σ_k u0·h0·dy.
      uhbt = 0.5_wp*u0*h0*real(NZ, wp)*1000.0_wp
      do k = 1, NZ
         vr(:, :, k) = gamma_k(k)
      end do
      u_cor = 0.0_wp

      call map_in(ms, ct)
      !$acc enter data copyin(uhbt, vr, u_cor)
      call continuity_zonal_flux(grid, metrics, ct, ms, dt, uhbt=uhbt, &
                                 visc_rem=vr, u_cor=u_cor)
      !$acc update self(u_cor)
      !$acc exit data delete(uhbt, vr, u_cor)
      call map_out(ms, ct)
   end subroutine setup_renorm_vr_case

   subroutine test_renorm_vr_transport(error)
      !! SPEC S2b accept (1): with γ-weighting on, the renormalised
      !! per-layer fluxes still vertically sum to `uhbt` to RENORM_TOL
      !! at every interior face, AND the returned `u_cor` reproduces
      !! that transport (Σ_k u_cor·h·dy = uhbt) — the property MOM6's
      !! `u_av` relies on.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      real(wp), allocatable :: uhbt(:, :), vr(:, :, :), u_cor(:, :, :)
      real(wp), parameter :: GAMMA(NZ) = [0.1_wp, 0.5_wp, 1.0_wp]
      real(wp), parameter :: U0 = 0.4_wp, H0 = 50.0_wp, DT = 100.0_wp
      real(wp) :: sum_flux, sum_ucor_h, err_flux, err_ucor
      integer :: i, j, k, i0, i1, j0, j1

      call setup_renorm_vr_case(grid, metrics, ms, ct, uhbt, vr, u_cor, &
                                GAMMA, U0, H0, DT)

      ! Interior faces only (strictly inside the physical walls).
      i0 = grid%nghost + 2
      i1 = grid%nghost + grid%nx_phys
      j0 = grid%nghost + 1
      j1 = grid%nghost + grid%ny_phys
      err_flux = 0.0_wp
      err_ucor = 0.0_wp
      do j = j0, j1
         do i = i0, i1
            sum_flux = 0.0_wp
            sum_ucor_h = 0.0_wp
            do k = 1, NZ
               sum_flux = sum_flux + ms%mass_flux_x_layer(i, j, k)
               sum_ucor_h = sum_ucor_h + u_cor(i, j, k)*H0*1000.0_wp
            end do
            err_flux = max(err_flux, abs(sum_flux - uhbt(i, j))/abs(uhbt(i, j)))
            err_ucor = max(err_ucor, abs(sum_ucor_h - uhbt(i, j))/abs(uhbt(i, j)))
         end do
      end do
      call check(error, err_flux < 1.0e-11_wp, &
                 "gamma-weighted renormalised fluxes do not sum to uhbt")
      if (allocated(error)) return
      call check(error, err_ucor < 1.0e-11_wp, &
                 "u_cor transport does not reproduce uhbt")

      call ct%destroy(); call ms%destroy()
      call destroy_cartesian_metrics(metrics)
      deallocate (uhbt, vr, u_cor)
   end subroutine test_renorm_vr_transport

   subroutine test_renorm_vr_shares(error)
      !! SPEC S2b accept (2): the per-layer correction `u_cor − u0` is
      !! proportional to γ_k — the high-drag bottom layer (γ = 0.1)
      !! receives a 10× SMALLER share of `du` than the undamped surface
      !! layer (γ = 1), and the ratios match γ exactly.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      real(wp), allocatable :: uhbt(:, :), vr(:, :, :), u_cor(:, :, :)
      real(wp), parameter :: GAMMA(NZ) = [0.1_wp, 0.5_wp, 1.0_wp]
      real(wp), parameter :: U0 = 0.4_wp, H0 = 50.0_wp, DT = 100.0_wp
      real(wp) :: du_k(NZ), ratio_err
      integer :: i, j, k

      call setup_renorm_vr_case(grid, metrics, ms, ct, uhbt, vr, u_cor, &
                                GAMMA, U0, H0, DT)

      ! One representative interior face.
      i = grid%nghost + 3
      j = grid%nghost + 2
      do k = 1, NZ
         du_k(k) = u_cor(i, j, k) - U0
      end do
      ! The correction must be non-trivial (uhbt = half the transport
      ! ⇒ du < 0) and bottom |du| strictly the smallest.
      call check(error, du_k(NZ) < -1.0e-6_wp, &
                 "surface-layer correction unexpectedly ~0 — inversion inactive?")
      if (allocated(error)) return
      call check(error, abs(du_k(1)) < abs(du_k(NZ)), &
                 "high-drag bottom layer did not receive a smaller du share")
      if (allocated(error)) return
      ! Exact proportionality: du_k/γ_k identical across layers.
      ratio_err = max(abs(du_k(1)/GAMMA(1) - du_k(NZ)/GAMMA(NZ)), &
                      abs(du_k(2)/GAMMA(2) - du_k(NZ)/GAMMA(NZ)))
      call check(error, ratio_err < 1.0e-12_wp*max(1.0_wp, abs(du_k(NZ))), &
                 "per-layer du shares are not proportional to visc_rem")

      call ct%destroy(); call ms%destroy()
      call destroy_cartesian_metrics(metrics)
      deallocate (uhbt, vr, u_cor)
   end subroutine test_renorm_vr_shares

end module test_continuity_multilayer
