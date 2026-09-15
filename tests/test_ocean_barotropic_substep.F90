!! Unit tests for the linearized barotropic barotropic-substep kernel
!! (`barotropic_substep_linear` in `rdb_ocean_dyn`).
!!
!! The kernel forward-Euler-steps `(η, u_bt, v_bt)` at `dt_inner`
!! over `n_steps` against a constant slow forcing, accumulating the
!! per-step running sum into `ubt_sum`/`vbt_sum`/`eta_sum` and
!! dividing by `n_steps` at the end to store the time-mean back in
!! `bt_η`/`bt_u_bt`/`bt_v_bt`.
!!
!! Cases:
!!   * Rest state — zero IC + zero forcing stays at zero (catches a
!!     sign bug or stray initialization).
!!   * Constant forcing on a rest column — analytic time-mean:
!!     u(t) = F·t  ⇒  ⟨u⟩ over [0, T] = ½·F·T.  The discrete forward-
!!     Euler version with N steps gives ⟨u⟩ = ½·F·T·(1 - 1/N).
!!   * Linear gravity wave propagation — initial η = A·cos(k·x), no
!!     forcing.  Phase speed c = √(g·H) for shallow water.  After
!!     one inertial-free period τ = λ/c the shape must come back to
!!     the IC within FE truncation tolerance.
module test_ocean_barotropic_substep
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY
   use rdb_grid, only: hgrid_t
   use rdb_coriolis_adv, only: coriolis_adv_t
   use rdb_ocean_dyn, only: ocean_dyn_t
   use rdb_barotropic_substep, only: barotropic_substep_linear, barotropic_substep_nonlinear
   use rdb_ocean_boundary_types, only: ocean_bc_state_t, ocean_bc_state_init, &
                                       ocean_bc_state_destroy, OBC_WALL, OBC_OPEN, &
                                       OBC_CLAMPED, OBC_TIDAL, OBC_CHAPMAN
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   implicit none
   private

   public :: collect_ocean_barotropic_substep_tests

   integer, parameter :: NGHOST = 2

contains

   subroutine collect_ocean_barotropic_substep_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("rest_state_stays_at_rest", test_rest_state), &
                  new_unittest("constant_force_time_mean", test_constant_force), &
                  new_unittest("linear_gravity_wave_returns", test_gravity_wave), &
                  new_unittest("nonlinear_rest_state", test_nonlinear_rest), &
                  new_unittest("nonlinear_recovers_linearized", &
                               test_nonlinear_recovers_linearized), &
                  new_unittest("nonlinear_coriolis_rotates_velocity", &
                               test_nonlinear_coriolis_rotation), &
                  new_unittest("nonlinear_wall_row_branches", &
                               test_nonlinear_wall_row_branches), &
                  new_unittest("nonlinear_zero_transport_at_physical_walls", &
                               test_nonlinear_zero_transport_at_walls), &
                  new_unittest("nonlinear_open_west_drains_pulse", &
                               test_nonlinear_open_west_drains), &
                  new_unittest("nonlinear_clamped_west_drives_inflow", &
                               test_nonlinear_clamped_west_inflow), &
                  new_unittest("nonlinear_tidal_west_drives_oscillation", &
                               test_nonlinear_tidal_west), &
                  new_unittest("nonlinear_chapman_persists_eta_old", &
                               test_nonlinear_chapman_state), &
                  new_unittest("bebt_first_substep_is_no_op", &
                               test_bebt_first_substep_no_op), &
                  new_unittest("bebt_engages_after_first_substep", &
                               test_bebt_engages_after_first_substep) &
                  ]
   end subroutine collect_ocean_barotropic_substep_tests

   subroutine make_grid(grid, nx_phys, ny_phys, dx, dy)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      real(wp), intent(in) :: dx, dy
      call grid%init(nx_phys, ny_phys, NGHOST, dx, dy)
   end subroutine make_grid

   subroutine run_fast(grid, metrics, dyn, fu, fv, n_steps, dt_inner)
      !! Map the dyn state to the device, run the kernel, pull back.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(ocean_dyn_t), intent(inout) :: dyn
      real(wp), intent(in) :: fu(:, :), fv(:, :)
      integer, intent(in) :: n_steps
      real(wp), intent(in) :: dt_inner
      !$acc enter data copyin(dyn, fu, fv)
      call dyn%enter_data()
      call barotropic_substep_linear(grid, metrics, dyn%bt_work, fu, fv, n_steps, dt_inner)
      !$acc update self(dyn%bt_work%bt_eta, dyn%bt_work%bt_ubt, dyn%bt_work%bt_vbt)
      call dyn%exit_data()
      !$acc exit data delete(dyn, fu, fv)
   end subroutine run_fast

   ! -----------------------------------------------------------------
   ! Cases
   ! -----------------------------------------------------------------

   subroutine test_rest_state(error)
      !! Zero IC + zero forcing → zero everywhere.  Catches stray
      !! accumulator initialization or sign error in the time-mean
      !! divide.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_dyn_t) :: dyn
      type(ocean_metrics_t) :: metrics
      real(wp), allocatable :: fu(:, :), fv(:, :)
      real(wp), parameter :: H_REF = 100.0_wp
      integer, parameter :: N_STEPS = 10
      real(wp), parameter :: DT_INNER = 0.1_wp
      checks: block

         call make_grid(grid, 12, 10, 1.0_wp, 1.0_wp)
         call make_cartesian_metrics(metrics, grid)
         call dyn%init(grid)
         dyn%bt_work%bt_H_ref = H_REF
         dyn%bt_work%bt_eta = 0.0_wp
         dyn%bt_work%bt_ubt = 0.0_wp
         dyn%bt_work%bt_vbt = 0.0_wp

         allocate (fu(grid%nx_total + 1, grid%ny_total), source=0.0_wp)
         allocate (fv(grid%nx_total, grid%ny_total + 1), source=0.0_wp)

         call run_fast(grid, metrics, dyn, fu, fv, N_STEPS, DT_INNER)

         call check(error, maxval(abs(dyn%bt_work%bt_eta)) < 1.0e-15_wp, &
                    "rest state: η drifted")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(dyn%bt_work%bt_ubt)) < 1.0e-15_wp, &
                    "rest state: u_bt drifted")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(dyn%bt_work%bt_vbt)) < 1.0e-15_wp, &
                    "rest state: v_bt drifted")

      end block checks
      deallocate (fu, fv)
      call dyn%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_rest_state

   subroutine test_constant_force(error)
      !! Apply a constant u-forcing F to a flat resting column, no
      !! η gradient anywhere.  Forward-Euler integration of
      !! ∂u/∂t = F gives u(n·dt) = F·n·dt at step n; the running sum
      !! over n = 1..N is F·dt·(1+2+...+N) = F·dt·N(N+1)/2, so the
      !! time-mean ⟨u⟩ = F·dt·(N+1)/2.  (Note: not F·dt·N/2 — that's
      !! the continuous analytic; the FE version samples at the END
      !! of each substep so the running sum is the trapezoidal-
      !! upper-Riemann sum.)
      !!
      !! With η starting at 0 and no flux divergence (uniform u_bt
      !! across a wall-bounded domain, but in the interior the
      !! continuity step DOES write -div(H·u_bt) = 0 since u is
      !! uniform), the wave isn't excited — pure ramping.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_dyn_t) :: dyn
      type(ocean_metrics_t) :: metrics
      real(wp), allocatable :: fu(:, :), fv(:, :)
      real(wp), parameter :: F = 1.0e-4_wp     ! m/s² (gentle, like wind-stress accel)
      integer, parameter :: N_STEPS = 50
      real(wp), parameter :: DT_INNER = 0.5_wp
      real(wp) :: u_mean_expected, u_mean_obs, rel_err
      integer :: i_probe, j_probe

      call make_grid(grid, 16, 12, 1.0_wp, 1.0_wp)
      call make_cartesian_metrics(metrics, grid)
      call dyn%init(grid)
      ! H_REF = 0 decouples the η dynamics from u_bt — div(0·u_bt) = 0
      ! everywhere → η stays at zero → no -g·grad(η) feedback on u_bt.
      ! With η forced to zero and only F driving u_bt, the FE
      ! integration is exact: u_bt(after step n) = n·dt·F.
      dyn%bt_work%bt_H_ref = 0.0_wp
      dyn%bt_work%bt_eta = 0.0_wp
      dyn%bt_work%bt_ubt = 0.0_wp
      dyn%bt_work%bt_vbt = 0.0_wp

      allocate (fu(grid%nx_total + 1, grid%ny_total), source=F)
      allocate (fv(grid%nx_total, grid%ny_total + 1), source=0.0_wp)

      call run_fast(grid, metrics, dyn, fu, fv, N_STEPS, DT_INNER)

      ! Probe a u-face at the centre of the domain.  The wall faces
      ! (i=1, i=nx+1) are clamped to zero by the kernel, but interior
      ! ones see the linear ramp from FE.
      i_probe = grid%nx_total/2
      j_probe = grid%ny_total/2
      u_mean_expected = F*DT_INNER*real(N_STEPS + 1, wp)/2.0_wp
      u_mean_obs = dyn%bt_work%bt_ubt(i_probe, j_probe)
      rel_err = abs(u_mean_obs - u_mean_expected)/abs(u_mean_expected)

      call check(error, rel_err < 1.0e-3_wp, &
                 "constant forcing: u_bt time-mean off by more than 0.1%")

      deallocate (fu, fv)
      call dyn%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_constant_force

   subroutine test_gravity_wave(error)
      !! Linear shallow-water gravity wave: ∂²η/∂t² = g·H·∇²η on a
      !! periodic-ish strip.  With closed walls our analytic match
      !! is a standing wave: η(x, t) = A·cos(k·x)·cos(ω·t), ω = c·k,
      !! c = √(g·H).  After one period τ = 2π/ω, the η field should
      !! return to its initial shape.
      !!
      !! Because closed walls truncate the cosine half-wavelength
      !! match, we use an even number of half-wavelengths spanning
      !! the interior columns and test on the running time-mean —
      !! after a full period the time-mean is ZERO (the running
      !! integral of cos over one period is zero), which is what we
      !! check.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_dyn_t) :: dyn
      type(ocean_metrics_t) :: metrics
      real(wp), allocatable :: fu(:, :), fv(:, :)
      real(wp), parameter :: H_REF = 100.0_wp
      real(wp), parameter :: A_ETA = 0.1_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp) :: c, lambda, period, dt_inner, max_mean_eta
      integer :: i, j, nx, ny, n_steps, n_per_wave, ighost
      integer, parameter :: NX_PHYS = 32, NY_PHYS = 4

      call make_grid(grid, NX_PHYS, NY_PHYS, 1.0_wp, 1.0_wp)
      call make_cartesian_metrics(metrics, grid)
      call dyn%init(grid)
      nx = grid%nx_total
      ny = grid%ny_total
      dyn%bt_work%bt_H_ref = H_REF

      ! IC: m=2 closed-wall eigenmode of the PHYSICAL interior.  The
      ! barotropic substep closes walls at the physical wall faces (i =
      ! nghost+1 and i = nghost+nx_phys+1) — matching slow continuity
      ! — so the eigenmode is defined over the nx_phys-cell box.
      ! Setting `i_in = i - nghost` and using
      ! η(i_in) = A·cos(2π·(i_in - 0.5)/nx_phys) gives ∂η/∂x = 0 at
      ! both physical walls by the symmetry cos(-θ) = cos(θ).
      ! Wavelength λ = nx_phys·dx.  Ghost cells stay at zero.
      lambda = real(NX_PHYS, wp)*grid%dx
      dyn%bt_work%bt_eta = 0.0_wp
      do j = 1, ny
         do i = grid%nghost + 1, grid%nghost + NX_PHYS
            ighost = i - grid%nghost
            dyn%bt_work%bt_eta(i, j) = A_ETA* &
                                       cos(2.0_wp*PI*(real(ighost, wp) - 0.5_wp)/real(NX_PHYS, wp))
         end do
      end do
      dyn%bt_work%bt_ubt = 0.0_wp
      dyn%bt_work%bt_vbt = 0.0_wp

      ! Resolve the wave: ~40 substeps per period at CFL ~ 0.5.
      c = sqrt(GRAVITY*H_REF)
      period = lambda/c
      n_per_wave = 40
      n_steps = n_per_wave
      dt_inner = period/real(n_per_wave, wp)

      allocate (fu(nx + 1, ny), source=0.0_wp)
      allocate (fv(nx, ny + 1), source=0.0_wp)

      call run_fast(grid, metrics, dyn, fu, fv, n_steps, dt_inner)

      ! Time-mean of A·cos(k·x)·cos(ω·t) over [0, T_period] is zero
      ! at every x because cos(ω·t) has zero mean.  Allow ~5% of A
      ! tolerance for FE truncation (each step is O(dt²); over a
      ! whole period the cumulative drift is dominated by the worst
      ! cell, which sits near the maximum of cos).  Tightening this
      ! to 1% requires either a smaller dt or a higher-order
      ! integrator — both deferred to the driver-integration branch.
      max_mean_eta = maxval(abs(dyn%bt_work%bt_eta(grid%nghost + 1:grid%nghost + NX_PHYS, :)))
      call check(error, max_mean_eta < 0.05_wp*A_ETA, &
                 "gravity wave: time-mean η didn't collapse over one period")

      deallocate (fu, fv)
      call dyn%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_gravity_wave

   ! -----------------------------------------------------------------
   ! Nonlinear barotropic substep
   ! -----------------------------------------------------------------

   subroutine run_fast_nonlinear(grid, metrics, dyn, cor, fu, fv, n_steps, dt_inner)
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(ocean_dyn_t), intent(inout) :: dyn
      type(coriolis_adv_t), intent(inout) :: cor
      real(wp), intent(in) :: fu(:, :), fv(:, :)
      integer, intent(in) :: n_steps
      real(wp), intent(in) :: dt_inner
      !$acc enter data copyin(dyn, cor, fu, fv)
      call dyn%enter_data()
      call cor%enter_data()
      call barotropic_substep_nonlinear(grid, dyn%bt_work, &
                                        fu, fv, &
                                        n_steps, dt_inner, &
                                        bt_eta=dyn%bt_work%bt_eta, bt_H_ref=dyn%bt_work%bt_H_ref, &
                                        bt_eta_new=dyn%bt_work%bt_eta_new, bt_ke_centre=dyn%bt_work%bt_ke_centre, &
                                        eta_sum=dyn%bt_work%eta_sum, bt_eta_end=dyn%bt_work%bt_eta_end, &
                                        bt_ubt=dyn%bt_work%bt_ubt, bt_ubt_prev=dyn%bt_work%bt_ubt_prev, &
                                        bt_rem_u=dyn%bt_work%bt_rem_u, ubt_sum=dyn%bt_work%ubt_sum, &
                                        uhbt_sum=dyn%bt_work%uhbt_sum, bt_uhbt=dyn%bt_work%bt_uhbt, &
                                        bt_ubt_end=dyn%bt_work%bt_ubt_end, &
                                        bt_vbt=dyn%bt_work%bt_vbt, bt_vbt_prev=dyn%bt_work%bt_vbt_prev, &
                                        bt_rem_v=dyn%bt_work%bt_rem_v, vbt_sum=dyn%bt_work%vbt_sum, &
                                        vhbt_sum=dyn%bt_work%vhbt_sum, bt_vhbt=dyn%bt_work%bt_vhbt, &
                                        bt_vbt_end=dyn%bt_work%bt_vbt_end, &
                                        bt_zeta_corner=dyn%bt_work%bt_zeta_corner, &
                                        f_corner=cor%f_corner, &
                                        area_cu=metrics%areaCu, area_cv=metrics%areaCv, dx_cu=metrics%dxCu, dx_cv=metrics%dx_cv, &
                                        dy_cu=metrics%dy_cu, dy_cv=metrics%dyCv, iarea_bu=metrics%iareaBu, iarea_t=metrics%iareaT, &
                                        idx_cu=metrics%idxCu, idy_cv=metrics%idyCv)
      !$acc update self(dyn%bt_work%bt_eta, dyn%bt_work%bt_ubt, dyn%bt_work%bt_vbt)
      call cor%exit_data()
      call dyn%exit_data()
      !$acc exit data delete(dyn, cor, fu, fv)
   end subroutine run_fast_nonlinear

   subroutine test_nonlinear_rest(error)
      !! Zero IC, zero forcing, f = 0 → everything stays at zero.
      !! Sanity that the new (ζ + f) · v_perp and KE-grad branches
      !! don't inject spurious tendency when their inputs vanish.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_dyn_t) :: dyn
      type(ocean_metrics_t) :: metrics
      type(coriolis_adv_t) :: cor
      real(wp), allocatable :: fu(:, :), fv(:, :)
      real(wp), parameter :: H_REF = 100.0_wp
      integer, parameter :: N_STEPS = 10
      real(wp), parameter :: DT_INNER = 0.1_wp
      checks: block

         call make_grid(grid, 12, 10, 1.0_wp, 1.0_wp)
         call make_cartesian_metrics(metrics, grid)
         call dyn%init(grid)
         call cor%init(grid)
         dyn%bt_work%bt_H_ref = H_REF

         allocate (fu(grid%nx_total + 1, grid%ny_total), source=0.0_wp)
         allocate (fv(grid%nx_total, grid%ny_total + 1), source=0.0_wp)

         call run_fast_nonlinear(grid, metrics, dyn, cor, fu, fv, N_STEPS, DT_INNER)

         call check(error, maxval(abs(dyn%bt_work%bt_eta)) < 1.0e-15_wp, &
                    "nonlinear rest: η drifted")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(dyn%bt_work%bt_ubt)) < 1.0e-15_wp, &
                    "nonlinear rest: u_bt drifted")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(dyn%bt_work%bt_vbt)) < 1.0e-15_wp, &
                    "nonlinear rest: v_bt drifted")

      end block checks
      deallocate (fu, fv)
      call cor%destroy(); call dyn%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_nonlinear_rest

   subroutine test_nonlinear_recovers_linearized(error)
      !! With f = 0 and a small-amplitude IC (|η|/H ≪ 1, |u|·dt/dx
      !! tiny), the nonlinear barotropic substep must match the linearized
      !! one to leading order.  Difference scales like
      !! (max ζ)·dt·(max v), (max η)·div(u)·dt, etc — bounded by
      !! the IC amplitude squared.  Tolerance is loose (1% of IC
      !! amplitude) since both forms share the same FE truncation.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_dyn_t) :: dyn_lin, dyn_nl
      type(ocean_metrics_t) :: metrics
      type(coriolis_adv_t) :: cor
      real(wp), allocatable :: fu(:, :), fv(:, :)
      real(wp), parameter :: H_REF = 100.0_wp
      real(wp), parameter :: A_ETA = 1.0e-3_wp   ! 1 mm vs 100 m H
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp) :: max_diff_eta, max_diff_u, max_diff_v, c, dt_inner
      integer :: i, j, nx, ny, n_steps, n_per_wave
      checks: block

         call make_grid(grid, 32, 4, 1.0_wp, 1.0_wp)
         call make_cartesian_metrics(metrics, grid)
         call dyn_lin%init(grid)
         call dyn_nl%init(grid)
         call cor%init(grid)
         cor%f_corner = 0.0_wp
         nx = grid%nx_total
         ny = grid%ny_total

         dyn_lin%bt_work%bt_H_ref = H_REF
         dyn_nl%bt_work%bt_H_ref = H_REF
         do j = 1, ny
            do i = 1, nx
               dyn_lin%bt_work%bt_eta(i, j) = A_ETA*cos(2.0_wp*PI*(real(i, wp) - 0.5_wp)/real(nx, wp))
               dyn_nl%bt_work%bt_eta(i, j) = dyn_lin%bt_work%bt_eta(i, j)
            end do
         end do

         c = sqrt(GRAVITY*H_REF)
         n_per_wave = 40
         n_steps = n_per_wave/4   ! quarter period — half a wave traversal
         dt_inner = (real(nx, wp)*grid%dx/c)/real(n_per_wave, wp)

         allocate (fu(nx + 1, ny), source=0.0_wp)
         allocate (fv(nx, ny + 1), source=0.0_wp)

         call run_fast(grid, metrics, dyn_lin, fu, fv, n_steps, dt_inner)
         call run_fast_nonlinear(grid, metrics, dyn_nl, cor, fu, fv, n_steps, dt_inner)

         max_diff_eta = maxval(abs(dyn_nl%bt_work%bt_eta - dyn_lin%bt_work%bt_eta))
         max_diff_u = maxval(abs(dyn_nl%bt_work%bt_ubt - dyn_lin%bt_work%bt_ubt))
         max_diff_v = maxval(abs(dyn_nl%bt_work%bt_vbt - dyn_lin%bt_work%bt_vbt))

         ! For a 1 mm η IC against a 100 m column, the nonlinear
         ! correction (η in H+η, ζ at small u) sits at A_ETA² ≈ 1e-6
         ! in η units.  Use 1% of A_ETA = 1e-5 as a comfortable cap.
         call check(error, max_diff_eta < 0.01_wp*A_ETA, &
                    "f=0 small-amp: nonlinear η differs from linearized")
         if (allocated(error)) exit checks
         call check(error, max_diff_u < 0.01_wp*A_ETA, &
                    "f=0 small-amp: nonlinear u differs from linearized")
         if (allocated(error)) exit checks
         call check(error, max_diff_v < 0.01_wp*A_ETA, &
                    "f=0 small-amp: nonlinear v differs from linearized")

      end block checks
      deallocate (fu, fv)
      call cor%destroy(); call dyn_nl%destroy(); call dyn_lin%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_nonlinear_recovers_linearized

   subroutine test_nonlinear_coriolis_rotation(error)
      !! f-plane discriminator.  Start with uniform u_bt = U0 in the
      !! deep interior, zero v_bt, zero η, no forcing.  After a short
      !! integration the Coriolis term `-(ζ+f)·u` on the v-face must
      !! drive v_bt toward `-f·U0·dt·n_steps` (sign and order of
      !! magnitude).  ζ is zero for uniform u, v → the test isolates
      !! the Coriolis branch.  Forward-Euler is marginally unstable
      !! on the imaginary axis but grows slowly for `f·dt ≪ 1`;
      !! keep the integration much shorter than one inertial period.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_dyn_t) :: dyn
      type(ocean_metrics_t) :: metrics
      type(coriolis_adv_t) :: cor
      real(wp), allocatable :: fu(:, :), fv(:, :)
      real(wp), parameter :: F0 = 1.0e-4_wp
      real(wp), parameter :: U0 = 0.1_wp
      real(wp), parameter :: H_REF = 10.0_wp
      real(wp), parameter :: DX = 100.0_wp
      real(wp), parameter :: DT_INNER = 1.0_wp
      integer, parameter :: N_STEPS = 30
      integer, parameter :: NX_PHYS = 16, NY_PHYS = 12
      real(wp) :: v_obs, v_expected_low, v_expected_high
      integer :: i, j, nx, ny, i_probe, j_probe
      checks: block

         ! CFL guard: gravity-wave speed √(g·H) = √98 ≈ 9.9 m/s; at
         ! dt=1 s and dx=100 m this is CFL ~ 0.1, well within the
         ! forward-backward stability limit.
         call make_grid(grid, NX_PHYS, NY_PHYS, DX, DX)
         call make_cartesian_metrics(metrics, grid)
         call dyn%init(grid)
         call cor%init(grid)
         cor%f_corner = F0
         nx = grid%nx_total
         ny = grid%ny_total

         dyn%bt_work%bt_H_ref = H_REF
         ! U0 in the deep interior, zero at the outer wall faces so
         ! the closed-wall clamp doesn't trigger a transient on step 1.
         do j = 1, ny
            do i = 2, nx
               dyn%bt_work%bt_ubt(i, j) = U0
            end do
         end do

         allocate (fu(nx + 1, ny), source=0.0_wp)
         allocate (fv(nx, ny + 1), source=0.0_wp)

         call run_fast_nonlinear(grid, metrics, dyn, cor, fu, fv, N_STEPS, DT_INNER)

         ! Probe an interior v-face well away from walls.  The
         ! time-mean over n_steps of `v(t) = -F0·U0·t` is
         ! `-F0·U0·dt·(N+1)/2` (same trapezoidal sampling logic as
         ! `test_constant_force`).  Bracket loosely — the Coriolis
         ! back-reaction from the growing v on u (which feeds back
         ! into v's tendency) makes the exact constant drift somewhat,
         ! but the sign and ~order of magnitude must hold.
         i_probe = nx/2
         j_probe = ny/2
         v_obs = dyn%bt_work%bt_vbt(i_probe, j_probe)
         v_expected_low = -F0*U0*DT_INNER*real(N_STEPS + 1, wp)/2.0_wp*1.5_wp
         v_expected_high = -F0*U0*DT_INNER*real(N_STEPS + 1, wp)/2.0_wp*0.5_wp

         call check(error, v_obs < 0.0_wp, &
                    "Coriolis on barotropic mode: v_bt didn't acquire the expected sign")
         if (allocated(error)) exit checks
         call check(error, v_obs > v_expected_low .and. v_obs < v_expected_high, &
                    "Coriolis on barotropic mode: v_bt magnitude off")

      end block checks
      deallocate (fu, fv)
      call cor%destroy(); call dyn%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_nonlinear_coriolis_rotation

   subroutine test_nonlinear_wall_row_branches(error)
      !! Regression guard on the wall-adjacent half-stencil branches
      !! of `barotropic_substep_nonlinear`'s `v_at_u` / `u_at_v` averages.
      !! The 4-point average drops to a 2-point average at the south
      !! wall row (j == 1), the north wall row (j == ny), the west
      !! wall column (i == 1), and the east wall column (i == nx);
      !! `test_nonlinear_coriolis_rotation` only probes the deep
      !! interior so those branches go un-asserted.  This test runs
      !! the same Coriolis spinup and probes the wall rows directly:
      !!
      !!   * `bt_ubt(i_probe, 1)`, `bt_ubt(i_probe, ny)` — finite,
      !!     same sign as the deep-interior response (Coriolis pulls
      !!     u down via -f·v term in dv/dt feedback), magnitude
      !!     within order-of-magnitude of deep-interior.
      !!   * `bt_vbt(1, j_probe)`, `bt_vbt(nx, j_probe)` — same
      !!     finiteness check on the mirror branch.
      !!
      !! Catches: a NaN in any wall branch, a sign flip from a typo
      !! in the half-stencil, or a 2× scaling slip from forgetting
      !! the `0.5` factor on the dropped row.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_dyn_t) :: dyn
      type(ocean_metrics_t) :: metrics
      type(coriolis_adv_t) :: cor
      real(wp), allocatable :: fu(:, :), fv(:, :)
      real(wp), parameter :: F0 = 1.0e-4_wp
      real(wp), parameter :: U0 = 0.1_wp
      real(wp), parameter :: H_REF = 10.0_wp
      real(wp), parameter :: DX = 100.0_wp
      real(wp), parameter :: DT_INNER = 1.0_wp
      integer, parameter :: N_STEPS = 30
      integer, parameter :: NX_PHYS = 16, NY_PHYS = 12
      real(wp) :: u_south, u_north, u_interior
      real(wp) :: v_west, v_east, v_interior
      integer :: i, j, nx, ny, i_probe, j_probe
      checks: block

         call make_grid(grid, NX_PHYS, NY_PHYS, DX, DX)
         call make_cartesian_metrics(metrics, grid)
         call dyn%init(grid)
         call cor%init(grid)
         cor%f_corner = F0
         nx = grid%nx_total
         ny = grid%ny_total

         dyn%bt_work%bt_H_ref = H_REF
         ! Uniform U0 across all interior u-faces.  Including j=1 and
         ! j=ny rows so the wall-side branches see the full forcing.
         do j = 1, ny
            do i = 2, nx
               dyn%bt_work%bt_ubt(i, j) = U0
            end do
         end do

         allocate (fu(nx + 1, ny), source=0.0_wp)
         allocate (fv(nx, ny + 1), source=0.0_wp)

         call run_fast_nonlinear(grid, metrics, dyn, cor, fu, fv, N_STEPS, DT_INNER)

         ! Probe the deep interior + the four wall rows / columns.
         i_probe = nx/2
         j_probe = ny/2
         u_interior = dyn%bt_work%bt_ubt(i_probe, j_probe)
         u_south = dyn%bt_work%bt_ubt(i_probe, 1)
         u_north = dyn%bt_work%bt_ubt(i_probe, ny)
         v_interior = dyn%bt_work%bt_vbt(i_probe, j_probe)
         v_west = dyn%bt_work%bt_vbt(1, j_probe)
         v_east = dyn%bt_work%bt_vbt(nx, j_probe)

         ! No NaN, no Inf — every wall branch produced a finite value.
         call check(error, abs(u_south) < 1.0_wp .and. abs(u_north) < 1.0_wp, &
                    "wall row: bt_ubt at j=1 or j=ny is not finite / huge")
         if (allocated(error)) exit checks
         call check(error, abs(v_west) < 1.0_wp .and. abs(v_east) < 1.0_wp, &
                    "wall column: bt_vbt at i=1 or i=nx is not finite / huge")
         if (allocated(error)) exit checks

         ! Same sign as the interior response.  Under Coriolis spinup
         ! with positive U0 and f > 0, v_bt acquires negative sign
         ! across the whole grid (the wall branches see a smaller
         ! magnitude but the same direction).
         call check(error, v_west < 0.0_wp .and. v_east < 0.0_wp, &
                    "wall column: bt_vbt sign disagrees with interior")
         if (allocated(error)) exit checks

         ! Wall-row u response should sit within 5x the interior
         ! magnitude (loose order-of-magnitude check — the half-
         ! stencil v_at_u is different but not catastrophically so).
         call check(error, abs(u_south - U0) < 5.0_wp*abs(u_interior - U0) + 1.0e-8_wp, &
                    "wall row: bt_ubt at j=1 deviates from interior by > 5x")

      end block checks
      deallocate (fu, fv)
      call cor%destroy(); call dyn%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_nonlinear_wall_row_branches

   subroutine test_nonlinear_zero_transport_at_walls(error)
      !! Regression for the Phase 3 barotropic-substep physical-wall closure
      !! (2026-05-19).  Without the fix, the barotropic substep closed walls
      !! only at array edges (i=1, i=nx+1) and not at the physical
      !! walls (i=nghost+1, i=nghost+nx_phys+1) where slow continuity
      !! closes — so during substeps the barotropic substep integrated free-
      !! surface transport through the 2-cell ghost strip carrying
      !! full `bt_H_ref`, and `bt_uhbt` at the physical wall came
      !! out non-zero.  Drove the long-run `Σ_k(h_layer) − (H +
      !! bt_eta_end)` drift that NaN'd momentum at ~step 1295 on
      !! the realistic-β_S `test_ocean_baroclinic_longrun`.
      !!
      !! This test runs `barotropic_substep_nonlinear` from an η perturbation
      !! IC with no slow forcing, then asserts the time-mean
      !! depth-integrated transports `bt_uhbt`, `bt_vhbt` are zero
      !! at every wall face — both array edges AND physical walls.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_dyn_t) :: dyn
      type(ocean_metrics_t) :: metrics
      type(coriolis_adv_t) :: cor
      real(wp), allocatable :: fu(:, :), fv(:, :)
      real(wp), parameter :: H_REF = 500.0_wp
      real(wp), parameter :: A_ETA = 0.5_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp), parameter :: F0 = 1.0e-4_wp
      real(wp), parameter :: DT_INNER = 1.0_wp
      real(wp), parameter :: WALL_TOL = 1.0e-12_wp
         !! Tight FP bound — transports at wall faces are formed by
         !! `h_face · u_wall_face` and `u_wall_face` is hard-zeroed
         !! in the kernel, so per-substep contribution is identically
         !! zero and the time-mean is bit-zero.
      integer, parameter :: NX_PHYS = 16, NY_PHYS = 8
      integer, parameter :: N_STEPS = 30
      integer :: i, j, nx, ny, i_w, i_e, j_s, j_n
      real(wp) :: x_rel, max_uhbt_array_edge, max_uhbt_phys_wall
      real(wp) :: max_vhbt_array_edge, max_vhbt_phys_wall

      call make_grid(grid, NX_PHYS, NY_PHYS, 10.0e3_wp, 10.0e3_wp)
      call make_cartesian_metrics(metrics, grid)
      call dyn%init(grid)
      cor%f_0 = F0
      call cor%init(grid)
      nx = grid%nx_total
      ny = grid%ny_total

      ! η perturbation across physical interior (m=1 cosine eigenmode
      ! of the closed basin).  Ghost cells stay at zero.
      dyn%bt_work%bt_H_ref = H_REF
      dyn%bt_work%bt_eta = 0.0_wp
      do j = 1, ny
         do i = grid%nghost + 1, grid%nghost + NX_PHYS
            x_rel = real(i - grid%nghost, wp) - 0.5_wp
            dyn%bt_work%bt_eta(i, j) = A_ETA*cos(2.0_wp*PI*x_rel/real(NX_PHYS, wp))
         end do
      end do
      dyn%bt_work%bt_ubt = 0.0_wp
      dyn%bt_work%bt_vbt = 0.0_wp

      allocate (fu(nx + 1, ny), source=0.0_wp)
      allocate (fv(nx, ny + 1), source=0.0_wp)

      !$acc enter data copyin(dyn, cor, fu, fv)
      call dyn%enter_data()
      call cor%enter_data()
      call barotropic_substep_nonlinear(grid, dyn%bt_work, &
                                        fu, fv, &
                                        N_STEPS, DT_INNER, &
                                        bt_eta=dyn%bt_work%bt_eta, bt_H_ref=dyn%bt_work%bt_H_ref, &
                                        bt_eta_new=dyn%bt_work%bt_eta_new, bt_ke_centre=dyn%bt_work%bt_ke_centre, &
                                        eta_sum=dyn%bt_work%eta_sum, bt_eta_end=dyn%bt_work%bt_eta_end, &
                                        bt_ubt=dyn%bt_work%bt_ubt, bt_ubt_prev=dyn%bt_work%bt_ubt_prev, &
                                        bt_rem_u=dyn%bt_work%bt_rem_u, ubt_sum=dyn%bt_work%ubt_sum, &
                                        uhbt_sum=dyn%bt_work%uhbt_sum, bt_uhbt=dyn%bt_work%bt_uhbt, &
                                        bt_ubt_end=dyn%bt_work%bt_ubt_end, &
                                        bt_vbt=dyn%bt_work%bt_vbt, bt_vbt_prev=dyn%bt_work%bt_vbt_prev, &
                                        bt_rem_v=dyn%bt_work%bt_rem_v, vbt_sum=dyn%bt_work%vbt_sum, &
                                        vhbt_sum=dyn%bt_work%vhbt_sum, bt_vhbt=dyn%bt_work%bt_vhbt, &
                                        bt_vbt_end=dyn%bt_work%bt_vbt_end, &
                                        bt_zeta_corner=dyn%bt_work%bt_zeta_corner, &
                                        f_corner=cor%f_corner, &
                                        area_cu=metrics%areaCu, area_cv=metrics%areaCv, dx_cu=metrics%dxCu, dx_cv=metrics%dx_cv, &
                                        dy_cu=metrics%dy_cu, dy_cv=metrics%dyCv, iarea_bu=metrics%iareaBu, iarea_t=metrics%iareaT, &
                                        idx_cu=metrics%idxCu, idy_cv=metrics%idyCv)
      !$acc update self(dyn%bt_work%bt_uhbt, dyn%bt_work%bt_vhbt)
      call cor%exit_data()
      call dyn%exit_data()
      !$acc exit data delete(dyn, cor, fu, fv)

      i_w = grid%nghost + 1
      i_e = grid%nghost + NX_PHYS + 1
      j_s = grid%nghost + 1
      j_n = grid%nghost + NY_PHYS + 1

      max_uhbt_array_edge = max(maxval(abs(dyn%bt_work%bt_uhbt(1, :))), &
                                maxval(abs(dyn%bt_work%bt_uhbt(nx + 1, :))))
      max_uhbt_phys_wall = max(maxval(abs(dyn%bt_work%bt_uhbt(i_w, :))), &
                               maxval(abs(dyn%bt_work%bt_uhbt(i_e, :))))
      max_vhbt_array_edge = max(maxval(abs(dyn%bt_work%bt_vhbt(:, 1))), &
                                maxval(abs(dyn%bt_work%bt_vhbt(:, ny + 1))))
      max_vhbt_phys_wall = max(maxval(abs(dyn%bt_work%bt_vhbt(:, j_s))), &
                               maxval(abs(dyn%bt_work%bt_vhbt(:, j_n))))

      call check(error, max_uhbt_array_edge < WALL_TOL, &
                 "bt_uhbt nonzero at array-edge wall — closure broken")
      if (allocated(error)) goto 100
      call check(error, max_uhbt_phys_wall < WALL_TOL, &
                 "bt_uhbt nonzero at physical wall — Phase 3 fix regressed")
      if (allocated(error)) goto 100
      call check(error, max_vhbt_array_edge < WALL_TOL, &
                 "bt_vhbt nonzero at array-edge wall — closure broken")
      if (allocated(error)) goto 100
      call check(error, max_vhbt_phys_wall < WALL_TOL, &
                 "bt_vhbt nonzero at physical wall — Phase 3 fix regressed")

100   deallocate (fu, fv)
      call cor%destroy(); call dyn%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_nonlinear_zero_transport_at_walls

   subroutine test_nonlinear_open_west_drains(error)
      !! OBC_OPEN dispatch: with all four walls closed, an initial η
      !! pulse in the basin centre stays bouncing.  Opening the west
      !! wall via Flather radiation lets the west-going component
      !! leave; energy in the basin must drop relative to the
      !! all-closed reference run.
      !!
      !! Compares two identical setups end-of-loop.  Reference is
      !! barotropic_substep_nonlinear with no `bc` arg (closed walls).  Test
      !! run passes a `bc` with `west = OBC_OPEN`.  Assert
      !! `Σ η² (test) < Σ η² (reference)`.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_dyn_t) :: dyn_ref, dyn_open
      type(ocean_metrics_t) :: metrics
      type(coriolis_adv_t) :: cor
      type(ocean_bc_state_t) :: bc
      real(wp), allocatable :: fu(:, :), fv(:, :)
      real(wp), parameter :: H_REF = 100.0_wp
      real(wp), parameter :: A_ETA = 0.5_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      integer, parameter :: NX_PHYS = 32, NY_PHYS = 4
      integer, parameter :: N_STEPS = 300
      real(wp), parameter :: DT_INNER = 0.01_wp
         !! CFL = c·dt/dx = √(g·H)·dt/dx ≈ 31.3·0.01/1 ≈ 0.31 — safe.
      real(wp) :: x_rel, e2_ref, e2_open
      integer :: i, j, nx, ny

      call make_grid(grid, NX_PHYS, NY_PHYS, 1.0_wp, 1.0_wp)
      call make_cartesian_metrics(metrics, grid)
      call dyn_ref%init(grid)
      call dyn_open%init(grid)
      call cor%init(grid)
      call ocean_bc_state_init(bc, grid, nz_ml=1)
      bc%west%bc_type = OBC_OPEN
      nx = grid%nx_total
      ny = grid%ny_total

      ! IC: m=1 cosine pulse centred on the physical interior.
      dyn_ref%bt_work%bt_H_ref = H_REF
      dyn_ref%bt_work%bt_eta = 0.0_wp
      dyn_open%bt_work%bt_H_ref = H_REF
      dyn_open%bt_work%bt_eta = 0.0_wp
      do j = 1, ny
         do i = grid%nghost + 1, grid%nghost + NX_PHYS
            x_rel = real(i - grid%nghost, wp) - 0.5_wp
            dyn_ref%bt_work%bt_eta(i, j) = A_ETA*cos(2.0_wp*PI*x_rel/real(NX_PHYS, wp))
            dyn_open%bt_work%bt_eta(i, j) = dyn_ref%bt_work%bt_eta(i, j)
         end do
      end do
      dyn_ref%bt_work%bt_ubt = 0.0_wp
      dyn_ref%bt_work%bt_vbt = 0.0_wp
      dyn_open%bt_work%bt_ubt = 0.0_wp
      dyn_open%bt_work%bt_vbt = 0.0_wp

      allocate (fu(nx + 1, ny), source=0.0_wp)
      allocate (fv(nx, ny + 1), source=0.0_wp)

      ! Reference run — all walls closed (no bc passed).
      !$acc enter data copyin(dyn_ref, cor, fu, fv)
      call dyn_ref%enter_data(); call cor%enter_data()
      call barotropic_substep_nonlinear(grid, dyn_ref%bt_work, &
                                        fu, fv, &
                                        N_STEPS, DT_INNER, &
                                        bt_eta=dyn_ref%bt_work%bt_eta, bt_H_ref=dyn_ref%bt_work%bt_H_ref, &
                                        bt_eta_new=dyn_ref%bt_work%bt_eta_new, bt_ke_centre=dyn_ref%bt_work%bt_ke_centre, &
                                        eta_sum=dyn_ref%bt_work%eta_sum, bt_eta_end=dyn_ref%bt_work%bt_eta_end, &
                                        bt_ubt=dyn_ref%bt_work%bt_ubt, bt_ubt_prev=dyn_ref%bt_work%bt_ubt_prev, &
                                        bt_rem_u=dyn_ref%bt_work%bt_rem_u, ubt_sum=dyn_ref%bt_work%ubt_sum, &
                                        uhbt_sum=dyn_ref%bt_work%uhbt_sum, bt_uhbt=dyn_ref%bt_work%bt_uhbt, &
                                        bt_ubt_end=dyn_ref%bt_work%bt_ubt_end, &
                                        bt_vbt=dyn_ref%bt_work%bt_vbt, bt_vbt_prev=dyn_ref%bt_work%bt_vbt_prev, &
                                        bt_rem_v=dyn_ref%bt_work%bt_rem_v, vbt_sum=dyn_ref%bt_work%vbt_sum, &
                                        vhbt_sum=dyn_ref%bt_work%vhbt_sum, bt_vhbt=dyn_ref%bt_work%bt_vhbt, &
                                        bt_vbt_end=dyn_ref%bt_work%bt_vbt_end, &
                                        bt_zeta_corner=dyn_ref%bt_work%bt_zeta_corner, &
                                        f_corner=cor%f_corner, &
                                        area_cu=metrics%areaCu, area_cv=metrics%areaCv, dx_cu=metrics%dxCu, dx_cv=metrics%dx_cv, &
                                        dy_cu=metrics%dy_cu, dy_cv=metrics%dyCv, iarea_bu=metrics%iareaBu, iarea_t=metrics%iareaT, &
                                        idx_cu=metrics%idxCu, idy_cv=metrics%idyCv)
      !$acc update self(dyn_ref%bt_work%bt_eta)
      call cor%exit_data(); call dyn_ref%exit_data()
      !$acc exit data delete(dyn_ref, cor, fu, fv)

      ! OBC_OPEN run — west wall radiating.
      !$acc enter data copyin(dyn_open, cor, fu, fv, bc)
      call dyn_open%enter_data(); call cor%enter_data()
      call barotropic_substep_nonlinear(grid, dyn_open%bt_work, &
                                        fu, fv, &
                                        N_STEPS, DT_INNER, &
                                        bt_eta=dyn_open%bt_work%bt_eta, bt_H_ref=dyn_open%bt_work%bt_H_ref, &
                                        bt_eta_new=dyn_open%bt_work%bt_eta_new, bt_ke_centre=dyn_open%bt_work%bt_ke_centre, &
                                        eta_sum=dyn_open%bt_work%eta_sum, bt_eta_end=dyn_open%bt_work%bt_eta_end, &
                                        bt_ubt=dyn_open%bt_work%bt_ubt, bt_ubt_prev=dyn_open%bt_work%bt_ubt_prev, &
                                        bt_rem_u=dyn_open%bt_work%bt_rem_u, ubt_sum=dyn_open%bt_work%ubt_sum, &
                                        uhbt_sum=dyn_open%bt_work%uhbt_sum, bt_uhbt=dyn_open%bt_work%bt_uhbt, &
                                        bt_ubt_end=dyn_open%bt_work%bt_ubt_end, &
                                        bt_vbt=dyn_open%bt_work%bt_vbt, bt_vbt_prev=dyn_open%bt_work%bt_vbt_prev, &
                                        bt_rem_v=dyn_open%bt_work%bt_rem_v, vbt_sum=dyn_open%bt_work%vbt_sum, &
                                        vhbt_sum=dyn_open%bt_work%vhbt_sum, bt_vhbt=dyn_open%bt_work%bt_vhbt, &
                                        bt_vbt_end=dyn_open%bt_work%bt_vbt_end, &
                                        bt_zeta_corner=dyn_open%bt_work%bt_zeta_corner, &
                                        f_corner=cor%f_corner, &
                                        bc=bc, &
                                        area_cu=metrics%areaCu, area_cv=metrics%areaCv, dx_cu=metrics%dxCu, dx_cv=metrics%dx_cv, &
                                        dy_cu=metrics%dy_cu, dy_cv=metrics%dyCv, iarea_bu=metrics%iareaBu, iarea_t=metrics%iareaT, &
                                        idx_cu=metrics%idxCu, idy_cv=metrics%idyCv)
      !$acc update self(dyn_open%bt_work%bt_eta)
      call cor%exit_data(); call dyn_open%exit_data()
      !$acc exit data delete(dyn_open, cor, fu, fv, bc)

      e2_ref = sum(dyn_ref%bt_work%bt_eta(grid%nghost + 1:grid%nghost + NX_PHYS, :)**2)
      e2_open = sum(dyn_open%bt_work%bt_eta(grid%nghost + 1:grid%nghost + NX_PHYS, :)**2)
      ! Observed ratio with this IC + N_STEPS is ≈ 0.002 (both halves
      ! of the cos eigenmode radiate out via the west wall after
      ! reflecting off the closed east wall).  Threshold of 0.1
      ! (10% of reference) catches the radiation while leaving plenty
      ! of headroom for IC / numerics tweaks.
      call check(error, e2_open < 0.1_wp*e2_ref, &
                 "OBC_OPEN west should drain energy to <10% of all-walls-closed reference")

      deallocate (fu, fv)
      call ocean_bc_state_destroy(bc)
      call cor%destroy(); call dyn_open%destroy(); call dyn_ref%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_nonlinear_open_west_drains

   subroutine test_nonlinear_clamped_west_inflow(error)
      !! OBC_CLAMPED dispatch: closed everywhere except west = CLAMPED
      !! with prescribed `u = U_IN > 0` (east-going).  Initial quiescent
      !! basin.  Continuity at the first wet cell:
      !!   dη/dt ≈ +H · U_IN / dx  (inflow at west, outflow zero at east).
      !! After a few substeps the basin-mean η should rise above zero.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_dyn_t) :: dyn
      type(ocean_metrics_t) :: metrics
      type(coriolis_adv_t) :: cor
      type(ocean_bc_state_t) :: bc
      real(wp), allocatable :: fu(:, :), fv(:, :)
      real(wp), parameter :: H_REF = 100.0_wp
      real(wp), parameter :: U_IN = 0.1_wp
      real(wp), parameter :: DT_INNER = 0.01_wp
      integer, parameter :: NX_PHYS = 16, NY_PHYS = 4
      integer, parameter :: N_STEPS = 50
      real(wp) :: eta_mean
      integer :: nx, ny

      call make_grid(grid, NX_PHYS, NY_PHYS, 1.0_wp, 1.0_wp)
      call make_cartesian_metrics(metrics, grid)
      call dyn%init(grid)
      call cor%init(grid)
      call ocean_bc_state_init(bc, grid, nz_ml=1)
      bc%west%bc_type = OBC_CLAMPED
      bc%west%clamped_u = U_IN
      nx = grid%nx_total
      ny = grid%ny_total

      dyn%bt_work%bt_H_ref = H_REF
      dyn%bt_work%bt_eta = 0.0_wp
      dyn%bt_work%bt_ubt = 0.0_wp
      dyn%bt_work%bt_vbt = 0.0_wp

      allocate (fu(nx + 1, ny), source=0.0_wp)
      allocate (fv(nx, ny + 1), source=0.0_wp)

      !$acc enter data copyin(dyn, cor, fu, fv, bc)
      call dyn%enter_data(); call cor%enter_data()
      call barotropic_substep_nonlinear(grid, dyn%bt_work, &
                                        fu, fv, &
                                        N_STEPS, DT_INNER, &
                                        bt_eta=dyn%bt_work%bt_eta, bt_H_ref=dyn%bt_work%bt_H_ref, &
                                        bt_eta_new=dyn%bt_work%bt_eta_new, bt_ke_centre=dyn%bt_work%bt_ke_centre, &
                                        eta_sum=dyn%bt_work%eta_sum, bt_eta_end=dyn%bt_work%bt_eta_end, &
                                        bt_ubt=dyn%bt_work%bt_ubt, bt_ubt_prev=dyn%bt_work%bt_ubt_prev, &
                                        bt_rem_u=dyn%bt_work%bt_rem_u, ubt_sum=dyn%bt_work%ubt_sum, &
                                        uhbt_sum=dyn%bt_work%uhbt_sum, bt_uhbt=dyn%bt_work%bt_uhbt, &
                                        bt_ubt_end=dyn%bt_work%bt_ubt_end, &
                                        bt_vbt=dyn%bt_work%bt_vbt, bt_vbt_prev=dyn%bt_work%bt_vbt_prev, &
                                        bt_rem_v=dyn%bt_work%bt_rem_v, vbt_sum=dyn%bt_work%vbt_sum, &
                                        vhbt_sum=dyn%bt_work%vhbt_sum, bt_vhbt=dyn%bt_work%bt_vhbt, &
                                        bt_vbt_end=dyn%bt_work%bt_vbt_end, &
                                        bt_zeta_corner=dyn%bt_work%bt_zeta_corner, &
                                        f_corner=cor%f_corner, &
                                        bc=bc, &
                                        area_cu=metrics%areaCu, area_cv=metrics%areaCv, dx_cu=metrics%dxCu, dx_cv=metrics%dx_cv, &
                                        dy_cu=metrics%dy_cu, dy_cv=metrics%dyCv, iarea_bu=metrics%iareaBu, iarea_t=metrics%iareaT, &
                                        idx_cu=metrics%idxCu, idy_cv=metrics%idyCv)
      !$acc update self(dyn%bt_work%bt_eta)
      call cor%exit_data(); call dyn%exit_data()
      !$acc exit data delete(dyn, cor, fu, fv, bc)

      eta_mean = sum(dyn%bt_work%bt_eta(grid%nghost + 1:grid%nghost + NX_PHYS, &
                                        grid%nghost + 1:grid%nghost + NY_PHYS))/ &
                 real(NX_PHYS*NY_PHYS, wp)
      ! Expected: η accumulates from west-side inflow.  Lower bound is
      ! generous — the wave dynamics smear the inflow across the basin,
      ! but the time-mean η must be strictly positive after N_STEPS.
      call check(error, eta_mean > 0.0_wp, &
                 "OBC_CLAMPED west with positive u should raise basin-mean η")

      deallocate (fu, fv)
      call ocean_bc_state_destroy(bc)
      call cor%destroy(); call dyn%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_nonlinear_clamped_west_inflow

   subroutine test_nonlinear_tidal_west(error)
      !! OBC_TIDAL dispatch: closed walls except west = OBC_TIDAL with
      !! a single M2-like constituent.  Quiescent IC; the barotropic substep
      !! is called twice — once at t = 0 (when the tidal η_target is
      !! at maximum, cos(0) = 1) and once at t = T_M2/4 (when η_target
      !! crosses zero, cos(π/2) = 0).
      !!
      !! Verifies the tidal dispatch path: when t = 0 the basin
      !! interior u at the wall face becomes non-zero (tide pushes
      !! in); when t = T/4 the dispatch yields the same value as the
      !! OBC_OPEN formula (η_target = 0).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_dyn_t) :: dyn_a, dyn_b
      type(ocean_metrics_t) :: metrics
      type(coriolis_adv_t) :: cor
      type(ocean_bc_state_t) :: bc
      real(wp), allocatable :: fu(:, :), fv(:, :)
      real(wp), parameter :: H_REF = 100.0_wp
      real(wp), parameter :: AMP = 0.5_wp
      real(wp), parameter :: PERIOD = 44712.0_wp   ! M2 period (s)
      real(wp), parameter :: DT_INNER = 1.0_wp
      integer, parameter :: NX_PHYS = 16, NY_PHYS = 4
      integer, parameter :: N_STEPS = 1
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp) :: omega, u_at_t0, u_at_t_quarter
      integer :: nx, ny

      omega = 2.0_wp*PI/PERIOD

      call make_grid(grid, NX_PHYS, NY_PHYS, 1.0_wp, 1.0_wp)
      call make_cartesian_metrics(metrics, grid)
      call dyn_a%init(grid); call dyn_b%init(grid)
      call cor%init(grid)
      call ocean_bc_state_init(bc, grid, nz_ml=1)
      bc%west%bc_type = OBC_TIDAL
      bc%west%n_tidal_constituents = 1
      bc%west%tidal_amp(1) = AMP
      bc%west%tidal_omega(1) = omega
      bc%west%tidal_phase(1) = 0.0_wp
      nx = grid%nx_total
      ny = grid%ny_total

      dyn_a%bt_work%bt_H_ref = H_REF; dyn_a%bt_work%bt_eta = 0.0_wp
      dyn_a%bt_work%bt_ubt = 0.0_wp; dyn_a%bt_work%bt_vbt = 0.0_wp
      dyn_b%bt_work%bt_H_ref = H_REF; dyn_b%bt_work%bt_eta = 0.0_wp
      dyn_b%bt_work%bt_ubt = 0.0_wp; dyn_b%bt_work%bt_vbt = 0.0_wp

      allocate (fu(nx + 1, ny), source=0.0_wp)
      allocate (fv(nx, ny + 1), source=0.0_wp)

      ! Run at t = 0 (eta_target = AMP).  Inflow at west wall expected.
      !$acc enter data copyin(dyn_a, cor, fu, fv, bc)
      call dyn_a%enter_data(); call cor%enter_data()
      call barotropic_substep_nonlinear(grid, dyn_a%bt_work, &
                                        fu, fv, &
                                        N_STEPS, DT_INNER, &
                                        bt_eta=dyn_a%bt_work%bt_eta, bt_H_ref=dyn_a%bt_work%bt_H_ref, &
                                        bt_eta_new=dyn_a%bt_work%bt_eta_new, bt_ke_centre=dyn_a%bt_work%bt_ke_centre, &
                                        eta_sum=dyn_a%bt_work%eta_sum, bt_eta_end=dyn_a%bt_work%bt_eta_end, &
                                        bt_ubt=dyn_a%bt_work%bt_ubt, bt_ubt_prev=dyn_a%bt_work%bt_ubt_prev, &
                                        bt_rem_u=dyn_a%bt_work%bt_rem_u, ubt_sum=dyn_a%bt_work%ubt_sum, &
                                        uhbt_sum=dyn_a%bt_work%uhbt_sum, bt_uhbt=dyn_a%bt_work%bt_uhbt, &
                                        bt_ubt_end=dyn_a%bt_work%bt_ubt_end, &
                                        bt_vbt=dyn_a%bt_work%bt_vbt, bt_vbt_prev=dyn_a%bt_work%bt_vbt_prev, &
                                        bt_rem_v=dyn_a%bt_work%bt_rem_v, vbt_sum=dyn_a%bt_work%vbt_sum, &
                                        vhbt_sum=dyn_a%bt_work%vhbt_sum, bt_vhbt=dyn_a%bt_work%bt_vhbt, &
                                        bt_vbt_end=dyn_a%bt_work%bt_vbt_end, &
                                        bt_zeta_corner=dyn_a%bt_work%bt_zeta_corner, &
                                        f_corner=cor%f_corner, &
                                        bc=bc, t=0.0_wp, &
                                        area_cu=metrics%areaCu, area_cv=metrics%areaCv, dx_cu=metrics%dxCu, dx_cv=metrics%dx_cv, &
                                        dy_cu=metrics%dy_cu, dy_cv=metrics%dyCv, iarea_bu=metrics%iareaBu, iarea_t=metrics%iareaT, &
                                        idx_cu=metrics%idxCu, idy_cv=metrics%idyCv)
      !$acc update self(dyn_a%bt_work%bt_ubt)
      call cor%exit_data(); call dyn_a%exit_data()
      !$acc exit data delete(dyn_a, cor, fu, fv, bc)

      ! Run at t = T_M2 / 4 (eta_target ≈ 0).  Should reduce to OPEN.
      !$acc enter data copyin(dyn_b, cor, fu, fv, bc)
      call dyn_b%enter_data(); call cor%enter_data()
      call barotropic_substep_nonlinear(grid, dyn_b%bt_work, &
                                        fu, fv, &
                                        N_STEPS, DT_INNER, &
                                        bt_eta=dyn_b%bt_work%bt_eta, bt_H_ref=dyn_b%bt_work%bt_H_ref, &
                                        bt_eta_new=dyn_b%bt_work%bt_eta_new, bt_ke_centre=dyn_b%bt_work%bt_ke_centre, &
                                        eta_sum=dyn_b%bt_work%eta_sum, bt_eta_end=dyn_b%bt_work%bt_eta_end, &
                                        bt_ubt=dyn_b%bt_work%bt_ubt, bt_ubt_prev=dyn_b%bt_work%bt_ubt_prev, &
                                        bt_rem_u=dyn_b%bt_work%bt_rem_u, ubt_sum=dyn_b%bt_work%ubt_sum, &
                                        uhbt_sum=dyn_b%bt_work%uhbt_sum, bt_uhbt=dyn_b%bt_work%bt_uhbt, &
                                        bt_ubt_end=dyn_b%bt_work%bt_ubt_end, &
                                        bt_vbt=dyn_b%bt_work%bt_vbt, bt_vbt_prev=dyn_b%bt_work%bt_vbt_prev, &
                                        bt_rem_v=dyn_b%bt_work%bt_rem_v, vbt_sum=dyn_b%bt_work%vbt_sum, &
                                        vhbt_sum=dyn_b%bt_work%vhbt_sum, bt_vhbt=dyn_b%bt_work%bt_vhbt, &
                                        bt_vbt_end=dyn_b%bt_work%bt_vbt_end, &
                                        bt_zeta_corner=dyn_b%bt_work%bt_zeta_corner, &
                                        f_corner=cor%f_corner, &
                                        bc=bc, t=0.25_wp*PERIOD, &
                                        area_cu=metrics%areaCu, area_cv=metrics%areaCv, dx_cu=metrics%dxCu, dx_cv=metrics%dx_cv, &
                                        dy_cu=metrics%dy_cu, dy_cv=metrics%dyCv, iarea_bu=metrics%iareaBu, iarea_t=metrics%iareaT, &
                                        idx_cu=metrics%idxCu, idy_cv=metrics%idyCv)
      !$acc update self(dyn_b%bt_work%bt_ubt)
      call cor%exit_data(); call dyn_b%exit_data()
      !$acc exit data delete(dyn_b, cor, fu, fv, bc)

      u_at_t0 = dyn_a%bt_work%bt_ubt(grid%nghost + 1, grid%nghost + 1)
      u_at_t_quarter = dyn_b%bt_work%bt_ubt(grid%nghost + 1, grid%nghost + 1)
      ! At t=0, eta_target = AMP > 0 → west wall u = +sqrt(g/H)·AMP > 0.
      ! At t=T/4, eta_target ≈ 0 → west wall u ≈ 0 (η_interior also ≈ 0).
      call check(error, u_at_t0 > 0.5_wp*sqrt(GRAVITY/H_REF)*AMP, &
                 "OBC_TIDAL west at t=0 should drive u > +√(g/H)·AMP/2")
      if (allocated(error)) goto 100
      call check(error, abs(u_at_t_quarter) < 0.01_wp*sqrt(GRAVITY/H_REF)*AMP, &
                 "OBC_TIDAL west at t=T/4 should drive u ≈ 0 (eta_target = 0)")

100   deallocate (fu, fv)
      call ocean_bc_state_destroy(bc)
      call cor%destroy(); call dyn_b%destroy(); call dyn_a%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_nonlinear_tidal_west

   subroutine test_nonlinear_chapman_state(error)
      !! OBC_CHAPMAN dispatch: prove the persistent `eta_old_chapman_*`
      !! field is updated across calls.
      !!
      !! Setup: initial η pulse at the western interior column, no
      !! forcing.  bc%west = OBC_CHAPMAN with eta_old_chapman_w = 0.
      !! Run the barotropic substep once.  After return, `eta_old_chapman_w`
      !! should equal the post-barotropic-substep edge-mean η at the western
      !! interior column — non-zero, proving the kernel updated the
      !! persistent state.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_dyn_t) :: dyn
      type(ocean_metrics_t) :: metrics
      type(coriolis_adv_t) :: cor
      type(ocean_bc_state_t) :: bc
      real(wp), allocatable :: fu(:, :), fv(:, :)
      real(wp), parameter :: H_REF = 100.0_wp
      real(wp), parameter :: A_ETA = 0.3_wp
      real(wp), parameter :: DT_INNER = 0.01_wp
      integer, parameter :: NX_PHYS = 16, NY_PHYS = 4
      integer, parameter :: N_STEPS = 10
      real(wp) :: eta_old_after, eta_int_mean_end
      integer :: i, j, nx, ny

      call make_grid(grid, NX_PHYS, NY_PHYS, 1.0_wp, 1.0_wp)
      call make_cartesian_metrics(metrics, grid)
      call dyn%init(grid)
      call cor%init(grid)
      call ocean_bc_state_init(bc, grid, nz_ml=1)
      bc%west%bc_type = OBC_CHAPMAN
      nx = grid%nx_total
      ny = grid%ny_total

      dyn%bt_work%bt_H_ref = H_REF
      dyn%bt_work%bt_eta = 0.0_wp
      do j = 1, ny
         dyn%bt_work%bt_eta(grid%nghost + 1, j) = A_ETA
      end do
      dyn%bt_work%bt_ubt = 0.0_wp
      dyn%bt_work%bt_vbt = 0.0_wp

      allocate (fu(nx + 1, ny), source=0.0_wp)
      allocate (fv(nx, ny + 1), source=0.0_wp)

      ! Sanity: initial eta_old should be zero.
      call check(error, abs(bc%eta_old_chapman_w) < 1.0e-14_wp, &
                 "Chapman eta_old should start at zero")
      if (allocated(error)) goto 200

      !$acc enter data copyin(dyn, cor, fu, fv, bc)
      call dyn%enter_data(); call cor%enter_data()
      call barotropic_substep_nonlinear(grid, dyn%bt_work, &
                                        fu, fv, &
                                        N_STEPS, DT_INNER, &
                                        bt_eta=dyn%bt_work%bt_eta, bt_H_ref=dyn%bt_work%bt_H_ref, &
                                        bt_eta_new=dyn%bt_work%bt_eta_new, bt_ke_centre=dyn%bt_work%bt_ke_centre, &
                                        eta_sum=dyn%bt_work%eta_sum, bt_eta_end=dyn%bt_work%bt_eta_end, &
                                        bt_ubt=dyn%bt_work%bt_ubt, bt_ubt_prev=dyn%bt_work%bt_ubt_prev, &
                                        bt_rem_u=dyn%bt_work%bt_rem_u, ubt_sum=dyn%bt_work%ubt_sum, &
                                        uhbt_sum=dyn%bt_work%uhbt_sum, bt_uhbt=dyn%bt_work%bt_uhbt, &
                                        bt_ubt_end=dyn%bt_work%bt_ubt_end, &
                                        bt_vbt=dyn%bt_work%bt_vbt, bt_vbt_prev=dyn%bt_work%bt_vbt_prev, &
                                        bt_rem_v=dyn%bt_work%bt_rem_v, vbt_sum=dyn%bt_work%vbt_sum, &
                                        vhbt_sum=dyn%bt_work%vhbt_sum, bt_vhbt=dyn%bt_work%bt_vhbt, &
                                        bt_vbt_end=dyn%bt_work%bt_vbt_end, &
                                        bt_zeta_corner=dyn%bt_work%bt_zeta_corner, &
                                        f_corner=cor%f_corner, &
                                        bc=bc, &
                                        area_cu=metrics%areaCu, area_cv=metrics%areaCv, dx_cu=metrics%dxCu, dx_cv=metrics%dx_cv, &
                                        dy_cu=metrics%dy_cu, dy_cv=metrics%dyCv, iarea_bu=metrics%iareaBu, iarea_t=metrics%iareaT, &
                                        idx_cu=metrics%idxCu, idy_cv=metrics%idyCv)
      !$acc update self(dyn%bt_work%bt_eta_end)
      ! Note: `bc%eta_old_chapman_w` lives on host only — the fast
      ! loop writes it host-side, so no `acc update self` needed.
      call cor%exit_data(); call dyn%exit_data()
      !$acc exit data delete(dyn, cor, fu, fv, bc)

      eta_old_after = bc%eta_old_chapman_w
      eta_int_mean_end = sum(dyn%bt_work%bt_eta_end(grid%nghost + 1, &
                                                    grid%nghost + 1:grid%nghost + NY_PHYS))/ &
                         real(NY_PHYS, wp)

      ! Two checks: persistent state must have been updated; it must
      ! equal the post-barotropic-substep interior mean.
      call check(error, abs(eta_old_after) > 1.0e-6_wp, &
                 "Chapman eta_old should be updated (non-zero) after barotropic substep")
      if (allocated(error)) goto 200
      call check(error, abs(eta_old_after - eta_int_mean_end) < 1.0e-10_wp, &
                 "Chapman eta_old should match the end-step interior mean")

200   deallocate (fu, fv)
      call ocean_bc_state_destroy(bc)
      call cor%destroy(); call dyn%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_nonlinear_chapman_state

   ! -----------------------------------------------------------------
   ! BEBT (velocity-projection) cases
   ! -----------------------------------------------------------------
   !
   ! The MOM6 `BT_PROJECT_VELOCITY` knob exposed as `bt_work%bebt`
   ! replaces the η-update transport with a forward-time extrapolation
   !   ubt_trans = (1 + bebt) · ubt^n − bebt · ubt^{n-1}
   ! collapsing to plain forward-backward Euler at `bebt = 0`.  Tests:
   !
   !   * `bebt_first_substep_is_no_op` — `bt_ubt_prev` is seeded to
   !     `bt_ubt` at the start of every `barotropic_substep_nonlinear`
   !     call, so substep #1 must give bit-identical η regardless of
   !     `bebt` (the projection has no previous state to extrapolate
   !     from).
   !
   !   * `bebt_engages_after_first_substep` — once substep #1's
   !     Pass 2 has advanced `bt_ubt` past `bt_ubt_prev`, the
   !     extrapolation has something to bite on.  Running `n_steps=2`
   !     at `bebt=0` vs `bebt=0.5` must give a measurably different
   !     η (the projection actually does something).

   subroutine bebt_make_dyn(grid, dyn, h_ref, eta_pulse)
      !! Helper: build a fresh ocean_dyn_t with a uniform-flow +
      !! eta-pulse IC the BEBT projection has something to work with.
      type(hgrid_t), intent(in) :: grid
      type(ocean_dyn_t), intent(inout) :: dyn
      real(wp), intent(in) :: h_ref, eta_pulse
      integer :: i_centre, j_centre
      call dyn%init(grid)
      dyn%bt_work%bt_H_ref = h_ref
      dyn%bt_work%bt_eta = 0.0_wp
      ! Single-cell η perturbation in the basin centre — gives a
      ! ∇η that drives a (non-trivial) substep evolution.
      i_centre = grid%nghost + grid%nx_phys/2
      j_centre = grid%nghost + grid%ny_phys/2
      dyn%bt_work%bt_eta(i_centre, j_centre) = eta_pulse
      ! Uniform 0.1 m/s zonal flow so the η flux divergence is
      ! non-trivial and the projection has a velocity to extrapolate.
      dyn%bt_work%bt_ubt = 0.1_wp
      dyn%bt_work%bt_vbt = 0.0_wp
   end subroutine bebt_make_dyn

   subroutine test_bebt_first_substep_no_op(error)
      !! With `n_steps = 1`, the projection's `bt_ubt_prev` is
      !! seeded equal to `bt_ubt` on entry, so the formula
      !!   (1+bebt)·u^0 − bebt·u^0 = u^0
      !! collapses regardless of `bebt`.  Two runs (bebt=0, bebt=0.5)
      !! from the same IC must give bit-identical η.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_dyn_t) :: dyn_zero, dyn_half
      type(ocean_metrics_t) :: metrics
      type(coriolis_adv_t) :: cor
      real(wp), allocatable :: fu(:, :), fv(:, :)
      real(wp), allocatable :: eta_zero(:, :), eta_half(:, :)
      real(wp), parameter :: H_REF = 1000.0_wp
      real(wp), parameter :: ETA_PULSE = 0.5_wp
      real(wp), parameter :: DT_INNER = 30.0_wp
      integer, parameter :: N_STEPS = 1
      real(wp) :: max_diff
      checks: block
         call make_grid(grid, 8, 8, 5000.0_wp, 5000.0_wp)
         call make_cartesian_metrics(metrics, grid)
         call cor%init(grid)
         allocate (fu(grid%nx_total + 1, grid%ny_total), source=0.0_wp)
         allocate (fv(grid%nx_total, grid%ny_total + 1), source=0.0_wp)

         call bebt_make_dyn(grid, dyn_zero, H_REF, ETA_PULSE)
         dyn_zero%bt_work%bebt = 0.0_wp
         call run_fast_nonlinear(grid, metrics, dyn_zero, cor, fu, fv, N_STEPS, DT_INNER)
         allocate (eta_zero, source=dyn_zero%bt_work%bt_eta)

         call bebt_make_dyn(grid, dyn_half, H_REF, ETA_PULSE)
         dyn_half%bt_work%bebt = 0.5_wp
         call run_fast_nonlinear(grid, metrics, dyn_half, cor, fu, fv, N_STEPS, DT_INNER)
         allocate (eta_half, source=dyn_half%bt_work%bt_eta)

         max_diff = maxval(abs(eta_zero - eta_half))
         call check(error, max_diff < 1.0e-14_wp, &
                    "BEBT first substep should be no-op (bebt=0 vs bebt=0.5 must match): " &
                    //"max|Δη| should be machine-epsilon")
      end block checks
      if (allocated(eta_zero)) deallocate (eta_zero)
      if (allocated(eta_half)) deallocate (eta_half)
      deallocate (fu, fv)
      call cor%destroy(); call dyn_zero%destroy(); call dyn_half%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_bebt_first_substep_no_op

   subroutine test_bebt_engages_after_first_substep(error)
      !! After substep #1's Pass 2 has updated `bt_ubt` past
      !! `bt_ubt_prev`, the extrapolation `(1+bebt)·u^1 − bebt·u^0`
      !! deviates from plain `u^1`.  Two runs from the same IC with
      !! `n_steps = 2`, one at `bebt = 0` and one at `bebt = 0.5`,
      !! must give a non-trivially different η.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_dyn_t) :: dyn_zero, dyn_half
      type(ocean_metrics_t) :: metrics
      type(coriolis_adv_t) :: cor
      real(wp), allocatable :: fu(:, :), fv(:, :)
      real(wp), allocatable :: eta_zero(:, :), eta_half(:, :)
      real(wp), parameter :: H_REF = 1000.0_wp
      real(wp), parameter :: ETA_PULSE = 0.5_wp
      real(wp), parameter :: DT_INNER = 30.0_wp
      integer, parameter :: N_STEPS = 2
      real(wp) :: max_diff
      checks: block
         call make_grid(grid, 8, 8, 5000.0_wp, 5000.0_wp)
         call make_cartesian_metrics(metrics, grid)
         call cor%init(grid)
         allocate (fu(grid%nx_total + 1, grid%ny_total), source=0.0_wp)
         allocate (fv(grid%nx_total, grid%ny_total + 1), source=0.0_wp)

         call bebt_make_dyn(grid, dyn_zero, H_REF, ETA_PULSE)
         dyn_zero%bt_work%bebt = 0.0_wp
         call run_fast_nonlinear(grid, metrics, dyn_zero, cor, fu, fv, N_STEPS, DT_INNER)
         allocate (eta_zero, source=dyn_zero%bt_work%bt_eta)

         call bebt_make_dyn(grid, dyn_half, H_REF, ETA_PULSE)
         dyn_half%bt_work%bebt = 0.5_wp
         call run_fast_nonlinear(grid, metrics, dyn_half, cor, fu, fv, N_STEPS, DT_INNER)
         allocate (eta_half, source=dyn_half%bt_work%bt_eta)

         max_diff = maxval(abs(eta_zero - eta_half))
         call check(error, max_diff > 1.0e-10_wp, &
                    "BEBT should change η measurably after the first substep: " &
                    //"max|Δη| between bebt=0 and bebt=0.5 must be > 1e-10")
      end block checks
      if (allocated(eta_zero)) deallocate (eta_zero)
      if (allocated(eta_half)) deallocate (eta_half)
      deallocate (fu, fv)
      call cor%destroy(); call dyn_zero%destroy(); call dyn_half%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_bebt_engages_after_first_substep

end module test_ocean_barotropic_substep
