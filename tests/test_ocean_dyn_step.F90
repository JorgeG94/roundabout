!! Unit tests for the unsplit SSP-RK2 barotropic driver step
!! (rdb_ocean_dyn Phase 4a).  Couples continuity-PPM + Sadourny
!! Coriolis-adv into one second-order time step.  See
!! `docs/ROADMAP_OCEAN.md` Phase 4.
!!
!! Cases:
!!   * Coupled lake-at-rest — uniform h, zero u, zero v: the full
!!     RK2 driver must reproduce the state bit-for-bit (round-off
!!     floor), the joint test of "both kernels behave correctly on
!!     the trivial fixed point of the system".
!!   * Inertial oscillation through driver — uniform u = U0, v = 0,
!!     uniform h: probes an interior face well away from the wall
!!     stagger asymmetries and asserts the SSP-RK2 magnitude error
!!     is O(dt^4) rather than the O(dt^2) FE error.  Direct
!!     comparison to the analytic cos/sin rotation tolerates only
!!     1e-5 relative on |U| (FE produced ~0.8% drift over the same
!!     run; RK2 should sit at ~2e-7).
!!   * Closed-basin mass conservation through driver — non-trivial
!!     h IC + smooth wall-vanishing velocity field, run 100 outer
!!     steps, total mass drift < 1e-10 (round-off floor).  Same
!!     setup as the Phase 2c test, just routed through the driver.
module test_ocean_dyn_step
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_barotropic_state, only: barotropic_state_t
   use rdb_continuity, only: continuity_t
   use rdb_coriolis_adv, only: coriolis_adv_t
   use rdb_ocean_dyn, only: ocean_dyn_t, ocean_dyn_step_barotropic
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   implicit none
   private

   public :: collect_ocean_dyn_step_tests

   integer, parameter :: NGHOST = 2

contains

   subroutine collect_ocean_dyn_step_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("coupled_lake_at_rest", test_lake_at_rest), &
                  new_unittest("inertial_oscillation_rk2_accuracy", &
                               test_inertial_oscillation_rk2), &
                  new_unittest("closed_basin_mass_conservation", &
                               test_mass_conservation) &
                  ]
   end subroutine collect_ocean_dyn_step_tests

   subroutine make_grid(grid, nx_phys, ny_phys, dx, dy)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      real(wp), intent(in) :: dx, dy
      call grid%init(nx_phys, ny_phys, NGHOST, dx, dy)
   end subroutine make_grid

   subroutine map_in(bs, ct, cor)
      type(barotropic_state_t), intent(inout) :: bs
      type(continuity_t), intent(inout) :: ct
      type(coriolis_adv_t), intent(inout) :: cor
      !$acc enter data copyin(bs)
      call bs%enter_data()
      !$acc enter data copyin(ct)
      call ct%enter_data()
      !$acc enter data copyin(cor)
      call cor%enter_data()
   end subroutine map_in

   subroutine map_out(bs, ct, cor)
      type(barotropic_state_t), intent(inout) :: bs
      type(continuity_t), intent(inout) :: ct
      type(coriolis_adv_t), intent(inout) :: cor
      call cor%exit_data()
      !$acc exit data delete(cor)
      call ct%exit_data()
      !$acc exit data delete(ct)
      call bs%exit_data()
      !$acc exit data delete(bs)
   end subroutine map_out

   ! -----------------------------------------------------------------
   ! Cases
   ! -----------------------------------------------------------------

   subroutine test_lake_at_rest(error)
      !! Uniform h, zero velocity, zero Coriolis -> driver must leave
      !! the state unchanged at round-off.  Joint trivial fixed
      !! point of continuity-PPM (no flux divergence) and Sadourny
      !! Coriolis-adv (no acceleration).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(barotropic_state_t) :: bs
      type(continuity_t) :: ct
      type(coriolis_adv_t) :: cor
      type(ocean_dyn_t) :: dyn
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: H0 = 10.0_wp
      real(wp), parameter :: DT = 0.01_wp
      real(wp) :: max_dh, max_du, max_dv
      checks: block

         call make_grid(grid, 8, 6, 1.0_wp, 1.0_wp)
         call make_cartesian_metrics(metrics, grid)
         call bs%init(grid)
         call ct%init(grid)
         call cor%init(grid)
         call dyn%init(grid)
         bs%h = H0
         bs%u_face_x = 0.0_wp
         bs%v_face_y = 0.0_wp
         bs%coriolis_f = 1.0_wp

         call map_in(bs, ct, cor)
         call ocean_dyn_step_barotropic(grid, metrics, dyn, cor, ct, bs, DT)
         call map_out(bs, ct, cor)

         max_dh = maxval(abs(bs%h - H0))
         max_du = maxval(abs(bs%u_face_x))
         max_dv = maxval(abs(bs%v_face_y))
         call check(error, max_dh < 1.0e-14_wp, &
                    "coupled lake-at-rest: h drift exceeded round-off")
         if (allocated(error)) exit checks
         call check(error, max_du < 1.0e-14_wp, &
                    "coupled lake-at-rest: u drift exceeded round-off")
         if (allocated(error)) exit checks
         call check(error, max_dv < 1.0e-14_wp, &
                    "coupled lake-at-rest: v drift exceeded round-off")
         if (allocated(error)) exit checks

         ! Outer step counter must have advanced
         call check(error, dyn%outer_step_count == 1, &
                    "outer_step_count did not advance")

      end block checks
      call dyn%destroy()
      call cor%destroy()
      call ct%destroy()
      call bs%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_lake_at_rest

   subroutine test_inertial_oscillation_rk2(error)
      !! Inertial oscillation through the SSP-RK2 driver.  IC:
      !! uniform u = U0, v = 0, h = H0 (so the continuity step has
      !! zero divergence on the interior — the probe is far enough
      !! from the walls that the small wall-edge mass perturbations
      !! haven't propagated in over the integration horizon).
      !!
      !! Analytic solution at t = N*dt:
      !!   u(t) =  U0 * cos(f*t),  v(t) = -U0 * sin(f*t)
      !!
      !! Phase 3a's forward-Euler kernel gives ~0.8% magnitude
      !! growth over 157 steps at f*dt = 0.01.  Heun's method
      !! (SSP-RK2) should drop that to ~2e-7.  We demand magnitude
      !! error < 1e-5 (5000x tighter than the FE tolerance from
      !! the Phase 3a test).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(barotropic_state_t) :: bs
      type(continuity_t) :: ct
      type(coriolis_adv_t) :: cor
      type(ocean_dyn_t) :: dyn
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: F_C = 1.0_wp
      real(wp), parameter :: U0 = 1.0_wp
      real(wp), parameter :: H0 = 100.0_wp
      real(wp), parameter :: DT = 0.01_wp
      integer, parameter :: N_STEPS = 157
      real(wp) :: t_end, u_expected, v_expected, mag_expected
      real(wp) :: u_obs, v_obs, mag_obs, mag_err
      integer :: step, i_probe, j_probe, nx, ny
      checks: block

         call make_grid(grid, 16, 16, 1.0_wp, 1.0_wp)
         call make_cartesian_metrics(metrics, grid)
         call bs%init(grid)
         call ct%init(grid)
         call cor%init(grid)
         call dyn%init(grid)
         nx = grid%nx_total
         ny = grid%ny_total

         bs%h = H0
         bs%u_face_x = U0
         bs%v_face_y = 0.0_wp
         bs%coriolis_f = F_C

         call map_in(bs, ct, cor)
         do step = 1, N_STEPS
            call ocean_dyn_step_barotropic(grid, metrics, dyn, cor, ct, bs, DT)
         end do
         call map_out(bs, ct, cor)

         i_probe = nx/2
         j_probe = ny/2
         u_obs = bs%u_face_x(i_probe, j_probe)
         v_obs = bs%v_face_y(i_probe, j_probe)

         t_end = real(N_STEPS, wp)*DT
         u_expected = U0*cos(F_C*t_end)
         v_expected = -U0*sin(F_C*t_end)
         mag_obs = sqrt(u_obs*u_obs + v_obs*v_obs)
         mag_expected = sqrt(u_expected*u_expected + v_expected*v_expected)
         mag_err = abs(mag_obs - mag_expected)/mag_expected

         call check(error, mag_err < 1.0e-5_wp, &
                    "SSP-RK2 inertial oscillation magnitude error > 1e-5")
         if (allocated(error)) exit checks
         call check(error, dyn%outer_step_count == N_STEPS, &
                    "outer_step_count mismatch with N_STEPS")

      end block checks
      call dyn%destroy()
      call cor%destroy()
      call ct%destroy()
      call bs%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_inertial_oscillation_rk2

   subroutine test_mass_conservation(error)
      !! 100-step closed-basin run through the full driver with a
      !! non-trivial initial h and a smooth velocity field that
      !! vanishes at the walls.  Total mass = sum(h)*dx*dy must
      !! drift by less than 1e-10 relative — the round-off floor
      !! the bare continuity kernel demonstrated in Phase 2c, now
      !! demanded of the full Coriolis-advection-coupled driver.
      !!
      !! Sadourny + KE-gradient terms do not change h directly; they
      !! mutate (u, v) which then enters mass_flux_x = u*h_face at
      !! the next stage.  The flux divergence telescopes to wall
      !! fluxes (0) regardless of the velocity field, so total mass
      !! still conserves at round-off.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(barotropic_state_t) :: bs
      type(continuity_t) :: ct
      type(coriolis_adv_t) :: cor
      type(ocean_dyn_t) :: dyn
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: H_BASE = 10.0_wp
      real(wp), parameter :: H_AMP = 0.5_wp
      real(wp), parameter :: U_AMP = 0.3_wp
      real(wp), parameter :: F_C = 1.0_wp
      real(wp), parameter :: DT = 0.05_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      integer, parameter :: N_STEPS = 100
      real(wp) :: total_initial, total_final, drift
      real(wp) :: h_min, h_max
      integer :: i, j, nx, ny, step
      checks: block

         call make_grid(grid, 32, 16, 1.0_wp, 1.0_wp)
         call make_cartesian_metrics(metrics, grid)
         call bs%init(grid)
         call ct%init(grid)
         call cor%init(grid)
         call dyn%init(grid)
         nx = grid%nx_total
         ny = grid%ny_total

         do j = 1, ny
            do i = 1, nx
               bs%h(i, j) = H_BASE + H_AMP* &
                            sin(2.0_wp*PI*real(i, wp)/real(nx, wp))* &
                            cos(2.0_wp*PI*real(j, wp)/real(ny, wp))
            end do
         end do
         do j = 1, ny
            do i = 1, nx + 1
               bs%u_face_x(i, j) = U_AMP* &
                                   sin(PI*real(i - 1, wp)/real(nx, wp))* &
                                   sin(2.0_wp*PI*real(j - 1, wp)/real(ny - 1, wp))
            end do
         end do
         do j = 1, ny + 1
            do i = 1, nx
               bs%v_face_y(i, j) = U_AMP* &
                                   sin(PI*real(j - 1, wp)/real(ny, wp))* &
                                   sin(2.0_wp*PI*real(i - 1, wp)/real(nx - 1, wp))
            end do
         end do
         bs%coriolis_f = F_C

         total_initial = sum(bs%h)*grid%dx*grid%dy

         call map_in(bs, ct, cor)
         do step = 1, N_STEPS
            call ocean_dyn_step_barotropic(grid, metrics, dyn, cor, ct, bs, DT)
         end do
         call map_out(bs, ct, cor)

         total_final = sum(bs%h)*grid%dx*grid%dy
         drift = abs(total_final - total_initial)/abs(total_initial)
         h_min = minval(bs%h)
         h_max = maxval(bs%h)

         call check(error, drift < 1.0e-10_wp, &
                    "driver total mass drift exceeded 1e-10 round-off floor")
         if (allocated(error)) exit checks
         call check(error, h_min > 0.0_wp, &
                    "h went negative through the driver — kernel unstable")
         if (allocated(error)) exit checks
         call check(error, h_max < 100.0_wp*H_BASE, &
                    "h grew > 100x H_BASE through the driver — kernel unstable")

      end block checks
      call dyn%destroy()
      call cor%destroy()
      call ct%destroy()
      call bs%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_mass_conservation

end module test_ocean_dyn_step
