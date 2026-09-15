!! Barotropic linear (Rayleigh) wave drag (Egbert & Ray 2001; Jayne & St
!! Laurent 2001) — the bulk energy sink for the barotropic tide.
!!
!! `&ocean_bt_nml wave_drag` folds a static per-face piston-velocity map
!! `r_H(x,y)` [m/s] into `bt_rem_u/v` (`compute_bt_rem_wave_drag` in
!! `rdb_barotropic_coupling`), composing multiplicatively with the existing
!! `substep_drag` seam.  These tests drive the REAL
!! `barotropic_substep_nonlinear_interior` (not just the coefficient
!! formula) because the sharpest failure modes here are (a) `bt_rem`
!! compounding geometrically across outer steps when `substep_drag` is off
!! (the trap `reset_bt_rem` exists to prevent — see
!! `src/core/ocean/README.md`), and (b) a closure that silently scalarises
!! its spatial map (see `docs/RDB_PHYSICS_PLAN.md` trap #11).
!!
!! Cases (§9 of `PLAN_PR29_wave_drag.md`):
!!   1. `wave_drag_off_bitident`      — r_H≡0 (either enable state) must be
!!      EXACTLY bit-identical to wave_drag disabled, for both
!!      substep_drag=T/F, proving the new dispatch branches are true no-ops
!!      when inert.
!!   2. `wave_drag_plumbing_reaches_state` — `configure_ocean_wave_drag`
!!      actually writes `bt_work%lwd_enable` + `lwd_drag_u/v` (the wire, not
!!      just the ends — CLAUDE.md "trap #3").
!!   3. `free_decay_matches_closed_form` — uniform r_H: closed-form
!!      backward-Euler decay, round-off exact.
!!   4. `spatial_structure_survives`  — two-decade ramp in r_H, each face
!!      decays at its OWN local rate (round-off exact at step 1, before any
!!      spatial feedback can develop).
!!   5. `wave_drag_no_compounding`    — substep_drag off: two independent
!!      fills of `bt_rem` from the SAME state must give the SAME factor, not
!!      its square.
!!   6. `composes_with_substep_drag` — both on: bt_rem is the exact product
!!      of the two independent closed-form factors.
!!   7. `energy_identity_exact`      — the exact discrete backward-Euler
!!      energy identity, KE monotonicity, and first-order dt convergence of
!!      the recovered decay rate to r_H/H.
!!   8. `roughness_proxy_structure`  — the `form="roughness_proxy"` filler
!!      puts drag on a ridge's FLANKS (max |grad b|), ~0 on the flat plain,
!!      0 on land + the ghost ring, and respects the `h2_max` ceiling.
module test_ocean_wave_drag
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_coriolis_adv, only: coriolis_adv_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_dyn, only: ocean_dyn_t
   use rdb_barotropic_substep, only: barotropic_substep_nonlinear_interior
   use rdb_barotropic_coupling, only: compute_bt_rem, reset_bt_rem, &
                                      compute_bt_rem_wave_drag, mask_bt_rem
   use rdb_config, only: config_t
   use rdb_ocean_state, only: ocean_state_t
   use rdb_ocean_setup, only: configure_ocean_wave_drag, wave_drag_roughness_proxy
   implicit none
   private

   public :: collect_ocean_wave_drag_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NX_PHYS = 40, NY_PHYS = 8
   real(wp), parameter :: DX = 1000.0_wp, DY = 1000.0_wp
   real(wp), parameter :: H_REF = 500.0_wp
   integer, parameter :: NZ = 1
   real(wp), parameter :: U0 = 0.05_wp

contains

   subroutine collect_ocean_wave_drag_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("wave_drag_off_bitident", test_off_bitident), &
                  new_unittest("wave_drag_plumbing_reaches_state", test_plumbing), &
                  new_unittest("free_decay_matches_closed_form", test_free_decay), &
                  new_unittest("spatial_structure_survives", test_spatial_structure), &
                  new_unittest("wave_drag_no_compounding", test_no_compounding), &
                  new_unittest("composes_with_substep_drag", test_composes), &
                  new_unittest("energy_identity_exact", test_energy_identity), &
                  new_unittest("roughness_proxy_structure", test_roughness_proxy) &
                  ]
   end subroutine collect_ocean_wave_drag_tests

   ! -----------------------------------------------------------------
   ! Shared harness — flat H_REF, f=0 Cartesian box, single layer.
   ! -----------------------------------------------------------------

   subroutine build_harness(grid, metrics, dyn, cor, ms)
      type(hgrid_t), intent(out) :: grid
      type(ocean_metrics_t), intent(out) :: metrics
      type(ocean_dyn_t), intent(out) :: dyn
      type(coriolis_adv_t), intent(out) :: cor
      type(multilayer_state_t), intent(out) :: ms
      call grid%init(NX_PHYS, NY_PHYS, NGHOST, DX, DY)
      call make_cartesian_metrics(metrics, grid)
      call dyn%init(grid, nz_ml=NZ)
      call cor%init(grid, nz_ml=NZ)
      cor%f_corner = 0.0_wp
      ms%nz_ml = NZ
      call ms%init(grid)
      ms%h_layer = H_REF
      dyn%bt_work%bt_H_ref = H_REF
      dyn%bt_work%bt_ubt = 0.0_wp
      dyn%bt_work%bt_vbt = 0.0_wp
      dyn%bt_work%bt_eta = 0.0_wp
   end subroutine build_harness

   subroutine teardown_harness(metrics, dyn, cor, ms)
      type(ocean_metrics_t), intent(inout) :: metrics
      type(ocean_dyn_t), intent(inout) :: dyn
      type(coriolis_adv_t), intent(inout) :: cor
      type(multilayer_state_t), intent(inout) :: ms
      call dyn%destroy()
      call cor%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine teardown_harness

   subroutine map_harness(dyn, cor, ms)
      !! GPU `mem:separate`: map AFTER the host fills IC + lwd_drag_u/v (the
      !! guarded `copyin` in `barotropic_workstate_enter_data_impl` only
      !! picks up `lwd_drag_u/v` if they are already allocated+filled).
      type(ocean_dyn_t), intent(inout) :: dyn
      type(coriolis_adv_t), intent(inout) :: cor
      type(multilayer_state_t), intent(inout) :: ms
      !$acc enter data copyin(dyn, cor, ms)
      call dyn%enter_data()
      call cor%enter_data()
      call ms%enter_data()
   end subroutine map_harness

   subroutine unmap_harness(dyn, cor, ms)
      type(ocean_dyn_t), intent(inout) :: dyn
      type(coriolis_adv_t), intent(inout) :: cor
      type(multilayer_state_t), intent(inout) :: ms
      call ms%exit_data()
      call cor%exit_data()
      call dyn%exit_data()
      !$acc exit data delete(dyn, cor, ms)
   end subroutine unmap_harness

   subroutine run_wave_drag_dispatch(grid, metrics, dyn, ms, dt_inner, &
                                     bt_substep_drag, r_linear, hbbl)
      !! Mirrors the production 3-branch dispatch in
      !! `rdb_ocean_dyn.F90:run_stage_split` (§5.4 of the plan) byte-for-byte:
      !! `compute_bt_rem` fills `bt_rem` when `substep_drag` is on; otherwise
      !! `reset_bt_rem` resets it to 1 before the wave-drag MULTIPLY so the
      !! multiplicative accumulator never compounds across calls (§7.1's
      !! trap). `mask_bt_rem` always runs last.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(ocean_dyn_t), intent(inout) :: dyn
      type(multilayer_state_t), intent(in) :: ms
      real(wp), intent(in) :: dt_inner, r_linear, hbbl
      logical, intent(in) :: bt_substep_drag
      if (bt_substep_drag) then
         call compute_bt_rem(grid, dyn%bt_work, ms, r_linear, hbbl, dt_inner)
      else if (dyn%bt_work%lwd_enable) then
         call reset_bt_rem(grid, dyn%bt_work)
      end if
      if (dyn%bt_work%lwd_enable) then
         call compute_bt_rem_wave_drag(grid, dyn%bt_work, ms, dt_inner)
      end if
      call mask_bt_rem(grid, metrics, dyn%bt_work)
   end subroutine run_wave_drag_dispatch

   ! -----------------------------------------------------------------
   ! Test 1 — off (or r_H≡0) is bit-identical, both substep_drag states.
   ! -----------------------------------------------------------------

   subroutine test_off_bitident(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: u_off, u_on_zero
      call run_and_probe(error, .false., .false., u_off)
      if (allocated(error)) return
      call run_and_probe(error, .false., .true., u_on_zero)
      if (allocated(error)) return
      call check(error, u_on_zero == u_off, &
                 "substep_drag=F: wave_drag on with r_H=0 not bit-identical to off")
      if (allocated(error)) return
      call run_and_probe(error, .true., .false., u_off)
      if (allocated(error)) return
      call run_and_probe(error, .true., .true., u_on_zero)
      if (allocated(error)) return
      call check(error, u_on_zero == u_off, &
                 "substep_drag=T: wave_drag on with r_H=0 not bit-identical to off")
   end subroutine test_off_bitident

   subroutine run_and_probe(error, bt_substep_drag, lwd_enable, u_probe)
      type(error_type), allocatable, intent(out) :: error
      logical, intent(in) :: bt_substep_drag, lwd_enable
      real(wp), intent(out) :: u_probe
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ocean_dyn_t) :: dyn
      type(coriolis_adv_t) :: cor
      type(multilayer_state_t) :: ms
      real(wp), parameter :: DT_INNER = 30.0_wp
      integer, parameter :: N_STEPS = 5
      integer :: nx, ny, ip, jp

      call build_harness(grid, metrics, dyn, cor, ms)
      nx = grid%nx_total; ny = grid%ny_total
      dyn%bt_work%bt_ubt = U0
      dyn%bt_work%lwd_enable = lwd_enable
      if (lwd_enable) then
         allocate (dyn%bt_work%lwd_drag_u(nx + 1, ny), source=0.0_wp)
         allocate (dyn%bt_work%lwd_drag_v(nx, ny + 1), source=0.0_wp)
      end if

      call map_harness(dyn, cor, ms)
      call run_wave_drag_dispatch(grid, metrics, dyn, ms, DT_INNER, &
                                  bt_substep_drag, 2.5e-5_wp, 10.0_wp)
      call barotropic_substep_nonlinear_interior(grid, metrics, dyn%bt_work, &
                                                 cor%f_corner, N_STEPS, DT_INNER)
      ! `bt_ubt` holds the TIME-MEAN over the substep loop; `bt_ubt_end` is
      ! the end-of-loop snapshot (the actual final decayed value).
      !$acc update self(dyn%bt_work%bt_ubt_end)
      call unmap_harness(dyn, cor, ms)

      ip = grid%nghost + NX_PHYS/2
      jp = grid%nghost + NY_PHYS/2
      u_probe = dyn%bt_work%bt_ubt_end(ip, jp)
      call check(error, u_probe == u_probe, "u_probe is NaN")
      call teardown_harness(metrics, dyn, cor, ms)
   end subroutine run_and_probe

   ! -----------------------------------------------------------------
   ! Test 2 — configure_ocean_wave_drag actually reaches bt_work.
   ! -----------------------------------------------------------------

   subroutine test_plumbing(error)
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(ocean_state_t) :: state
      type(hgrid_t) :: grid
      integer, parameter :: NXP = 8, NYP = 6
      real(wp), parameter :: R_UNIFORM = 1.0e-3_wp, SCALE = 2.0_wp
      real(wp) :: expected

      checks: block
         cfg%sim_type = "ocean"
         cfg%nx = NXP; cfg%ny = NYP; cfg%dx = DX; cfg%dy = DY
         cfg%nz_layers = 1; cfg%nghost = NGHOST
         cfg%ocean%bt%wave_drag = .true.
         cfg%ocean%bt%wave_drag_form = "uniform"
         cfg%ocean%bt%wave_drag_r_uniform = R_UNIFORM
         cfg%ocean%bt%wave_drag_scale = SCALE

         call grid%init(NXP, NYP, NGHOST, DX, DY)
         state%multilayer%nz_ml = 1
         call state%init(grid)

         call configure_ocean_wave_drag(cfg, state, grid, compute_rank=0)

         call check(error, state%dyn%bt_work%lwd_enable, &
                    "configure_ocean_wave_drag did not set lwd_enable")
         if (allocated(error)) exit checks
         call check(error, allocated(state%dyn%bt_work%lwd_drag_u), &
                    "lwd_drag_u not allocated")
         if (allocated(error)) exit checks
         call check(error, allocated(state%dyn%bt_work%lwd_drag_v), &
                    "lwd_drag_v not allocated")
         if (allocated(error)) exit checks

         ! scale=2.0 (!=1) is deliberate: a scale=1 test cannot see a
         ! dropped multiply (§9 of the plan).
         expected = SCALE*R_UNIFORM
         call check(error, abs(maxval(state%dyn%bt_work%lwd_drag_u) - expected) < 1.0e-14_wp, &
                    "wave_drag_scale did not reach lwd_drag_u")
         if (allocated(error)) exit checks
         ! Uniform form -> every INTERIOR face (not the zero-init array edges)
         ! carries the identical value.
         call check(error, abs(minval(state%dyn%bt_work%lwd_drag_u(2:grid%nx_total, :)) &
                               - expected) < 1.0e-14_wp, &
                    "uniform wave_drag_form not uniform at interior faces")
      end block checks
      call state%destroy()
   end subroutine test_plumbing

   ! -----------------------------------------------------------------
   ! Test 3 — uniform r_H: exact closed-form backward-Euler decay.
   ! -----------------------------------------------------------------

   subroutine test_free_decay(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ocean_dyn_t) :: dyn
      type(coriolis_adv_t) :: cor
      type(multilayer_state_t) :: ms
      real(wp), parameter :: R_H = 1.0e-3_wp
      real(wp), parameter :: DT_INNER = 30.0_wp
      integer, parameter :: N_STEPS = 5
      integer :: nx, ny, ip, jp
      real(wp) :: expected, got

      checks: block
         call build_harness(grid, metrics, dyn, cor, ms)
         nx = grid%nx_total; ny = grid%ny_total
         dyn%bt_work%bt_ubt = U0
         dyn%bt_work%lwd_enable = .true.
         allocate (dyn%bt_work%lwd_drag_u(nx + 1, ny), source=R_H)
         allocate (dyn%bt_work%lwd_drag_v(nx, ny + 1), source=R_H)

         call map_harness(dyn, cor, ms)
         call run_wave_drag_dispatch(grid, metrics, dyn, ms, DT_INNER, &
                                     .false., 0.0_wp, 0.0_wp)
         call barotropic_substep_nonlinear_interior(grid, metrics, dyn%bt_work, &
                                                    cor%f_corner, N_STEPS, DT_INNER)
         ! bt_ubt is the substep-loop TIME-MEAN; bt_ubt_end is the
         ! end-of-loop (final decayed) value the closed form describes.
         !$acc update self(dyn%bt_work%bt_ubt_end)
         call unmap_harness(dyn, cor, ms)

         ! Uniform r_H over a uniform IC -> u stays exactly spatially uniform
         ! at every step (no eta/PGF feedback ever develops), so ANY interior
         ! probe matches the pure backward-Euler closed form to round-off.
         ip = grid%nghost + NX_PHYS/2
         jp = grid%nghost + NY_PHYS/2
         expected = U0*(H_REF/(H_REF + R_H*DT_INNER))**N_STEPS
         got = dyn%bt_work%bt_ubt_end(ip, jp)
         call check(error, abs(got - expected) < 1.0e-13_wp*abs(expected), &
                    "free decay does not match the closed-form backward-Euler solution")
      end block checks
      call teardown_harness(metrics, dyn, cor, ms)
   end subroutine test_free_decay

   ! -----------------------------------------------------------------
   ! Test 4 — spatially varying r_H, one substep (no cross-face feedback
   ! has had time to develop yet): each face must decay at ITS OWN local
   ! rate, not a globally-averaged one.
   ! -----------------------------------------------------------------

   subroutine test_spatial_structure(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ocean_dyn_t) :: dyn
      type(coriolis_adv_t) :: cor
      type(multilayer_state_t) :: ms
      real(wp), parameter :: DT_INNER = 30.0_wp
      integer, parameter :: N_STEPS = 1
      integer :: nx, ny, i, j, jp, i_zero
      real(wp) :: expected, got, r_h_i, frac

      checks: block
         call build_harness(grid, metrics, dyn, cor, ms)
         nx = grid%nx_total; ny = grid%ny_total
         dyn%bt_work%bt_ubt = U0
         dyn%bt_work%lwd_enable = .true.
         allocate (dyn%bt_work%lwd_drag_u(nx + 1, ny))
         allocate (dyn%bt_work%lwd_drag_v(nx, ny + 1), source=0.0_wp)
         ! Two-decade ramp in i: 1e-5 .. 1e-3.
         do i = 1, nx + 1
            frac = real(i - 1, wp)/real(nx, wp)
            dyn%bt_work%lwd_drag_u(i, :) = 1.0e-5_wp*10.0_wp**(2.0_wp*frac)
         end do
         ! One exact-zero face, deliberately.
         i_zero = grid%nghost + 5
         dyn%bt_work%lwd_drag_u(i_zero, :) = 0.0_wp

         call map_harness(dyn, cor, ms)
         call run_wave_drag_dispatch(grid, metrics, dyn, ms, DT_INNER, &
                                     .false., 0.0_wp, 0.0_wp)
         call barotropic_substep_nonlinear_interior(grid, metrics, dyn%bt_work, &
                                                    cor%f_corner, N_STEPS, DT_INNER)
         !$acc update self(dyn%bt_work%bt_ubt_end)
         call unmap_harness(dyn, cor, ms)

         jp = grid%nghost + NY_PHYS/2

         ! Representative faces across the ramp (low-, mid-, high-decade).
         do i = grid%nghost + 2, nx - grid%nghost, 10
            frac = real(i - 1, wp)/real(nx, wp)
            r_h_i = 1.0e-5_wp*10.0_wp**(2.0_wp*frac)
            expected = U0*(H_REF/(H_REF + r_h_i*DT_INNER))
            got = dyn%bt_work%bt_ubt_end(i, jp)
            call check(error, abs(got - expected) < 1.0e-13_wp*abs(expected), &
                       "spatially varying r_H: face does not match its own local decay")
            if (allocated(error)) exit checks
         end do

         ! The r_H=0 face must be EXACTLY unchanged.
         call check(error, dyn%bt_work%bt_ubt_end(i_zero, jp) == U0, &
                    "r_H=0 face was perturbed")
      end block checks
      call teardown_harness(metrics, dyn, cor, ms)
   end subroutine test_spatial_structure

   ! -----------------------------------------------------------------
   ! Test 5 — the §7.1 trap: substep_drag off, wave_drag on. Two
   ! independent fills from the SAME re-seeded bt_rem must give the SAME
   ! factor, never its square.
   ! -----------------------------------------------------------------

   subroutine test_no_compounding(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ocean_dyn_t) :: dyn
      type(coriolis_adv_t) :: cor
      type(multilayer_state_t) :: ms
      real(wp), parameter :: R_H = 1.0e-3_wp
      real(wp), parameter :: DT_INNER = 240.0_wp
      integer :: nx, ny, ip, jp
      real(wp) :: factor_a, factor_b

      checks: block
         call build_harness(grid, metrics, dyn, cor, ms)
         nx = grid%nx_total; ny = grid%ny_total
         dyn%bt_work%lwd_enable = .true.
         allocate (dyn%bt_work%lwd_drag_u(nx + 1, ny), source=R_H)
         allocate (dyn%bt_work%lwd_drag_v(nx, ny + 1), source=R_H)

         call map_harness(dyn, cor, ms)
         ip = grid%nghost + NX_PHYS/2
         jp = grid%nghost + NY_PHYS/2

         ! "Outer step" A.
         call run_wave_drag_dispatch(grid, metrics, dyn, ms, DT_INNER, &
                                     .false., 0.0_wp, 0.0_wp)
         !$acc update self(dyn%bt_work%bt_rem_u)
         factor_a = dyn%bt_work%bt_rem_u(ip, jp)

         ! "Outer step" B — re-run the SAME dispatch from the SAME bt_work,
         ! without any external reset.  If `reset_bt_rem` is missing/broken,
         ! this MULTIPLIES into the already-once-multiplied bt_rem_u,
         ! giving factor_a**2 instead of factor_a.
         call run_wave_drag_dispatch(grid, metrics, dyn, ms, DT_INNER, &
                                     .false., 0.0_wp, 0.0_wp)
         !$acc update self(dyn%bt_work%bt_rem_u)
         factor_b = dyn%bt_work%bt_rem_u(ip, jp)
         call unmap_harness(dyn, cor, ms)

         call check(error, factor_a < 1.0_wp - 1.0e-6_wp, &
                    "sanity: factor_a must be a genuine (non-unity) drag factor")
         if (allocated(error)) exit checks
         call check(error, abs(factor_b - factor_a) < 1.0e-14_wp, &
                    "bt_rem compounded across outer steps (reset_bt_rem missing?)")
         if (allocated(error)) exit checks
         call check(error, abs(factor_b - factor_a**2) > 1.0e-6_wp, &
                    "sanity: the compounding failure mode (factor_a**2) must be " &
                    //"distinguishable from the correct answer")
      end block checks
      call teardown_harness(metrics, dyn, cor, ms)
   end subroutine test_no_compounding

   ! -----------------------------------------------------------------
   ! Test 6 — composes multiplicatively with substep_drag.
   ! -----------------------------------------------------------------

   subroutine test_composes(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ocean_dyn_t) :: dyn
      type(coriolis_adv_t) :: cor
      type(multilayer_state_t) :: ms
      real(wp), parameter :: R_H = 1.0e-3_wp
      real(wp), parameter :: R_LINEAR = 2.5e-5_wp, HBBL = 10.0_wp
      real(wp), parameter :: DT_INNER = 240.0_wp
      integer :: nx, ny, ip, jp
      real(wp) :: got, f_wave, f_lin, expected

      checks: block
         call build_harness(grid, metrics, dyn, cor, ms)
         nx = grid%nx_total; ny = grid%ny_total
         dyn%bt_work%lwd_enable = .true.
         allocate (dyn%bt_work%lwd_drag_u(nx + 1, ny), source=R_H)
         allocate (dyn%bt_work%lwd_drag_v(nx, ny + 1), source=R_H)

         call map_harness(dyn, cor, ms)
         call run_wave_drag_dispatch(grid, metrics, dyn, ms, DT_INNER, &
                                     .true., R_LINEAR, HBBL)
         !$acc update self(dyn%bt_work%bt_rem_u)
         call unmap_harness(dyn, cor, ms)

         ip = grid%nghost + NX_PHYS/2
         jp = grid%nghost + NY_PHYS/2
         got = dyn%bt_work%bt_rem_u(ip, jp)

         f_lin = H_REF/(H_REF + R_LINEAR*HBBL*DT_INNER)
         f_wave = H_REF/(H_REF + R_H*DT_INNER)
         expected = f_lin*f_wave

         call check(error, abs(got - expected) < 1.0e-14_wp, &
                    "bt_rem is not the exact product of the two closed-form factors")
         if (allocated(error)) exit checks
         call check(error, abs(got - f_lin) > 1.0e-6_wp .and. abs(got - f_wave) > 1.0e-6_wp, &
                    "bt_rem equals only one factor -- the other was not composed in")
      end block checks
      call teardown_harness(metrics, dyn, cor, ms)
   end subroutine test_composes

   ! -----------------------------------------------------------------
   ! Test 7 — exact discrete energy identity + KE monotonicity + O(dt)
   ! convergence of the recovered decay rate to r_H/H.
   ! -----------------------------------------------------------------

   subroutine test_energy_identity(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: R_H = 1.0e-3_wp
      real(wp), parameter :: DT_A = 200.0_wp
      integer, parameter :: K_STEPS = 4
      real(wp) :: u(0:K_STEPS)
      real(wp) :: lhs, rhs, ratio1, ratio2, rate1, rate2, rate_exact, err1, err2
      integer :: n

      checks: block
         call run_successive_singlesteps(R_H, DT_A, K_STEPS, u)

         ! Exact discrete identity + monotone KE for every successive pair.
         do n = 0, K_STEPS - 1
            lhs = 0.5_wp*H_REF*(u(n + 1)**2 - u(n)**2)
            rhs = -R_H*u(n + 1)**2*DT_A - 0.5_wp*H_REF*(u(n + 1) - u(n))**2
            call check(error, abs(lhs - rhs) < 1.0e-13_wp*H_REF*U0**2, &
                       "discrete energy identity does not close exactly")
            if (allocated(error)) exit checks
            call check(error, abs(u(n + 1)) <= abs(u(n)) + 1.0e-15_wp, &
                       "KE is not monotonically non-increasing")
            if (allocated(error)) exit checks
         end do

         ! O(dt) convergence of the recovered decay rate to the continuous
         ! rate r_H/H (halve dt_bt -> error halves).
         ratio1 = run_single_step_ratio(R_H, DT_A)
         ratio2 = run_single_step_ratio(R_H, DT_A/2.0_wp)
         rate1 = -log(ratio1)/DT_A
         rate2 = -log(ratio2)/(DT_A/2.0_wp)
         rate_exact = R_H/H_REF
         err1 = abs(rate1 - rate_exact)
         err2 = abs(rate2 - rate_exact)
         call check(error, err1 > 0.0_wp, "err1 degenerate (dt too small to resolve)")
         if (allocated(error)) exit checks
         call check(error, abs(err2/err1 - 0.5_wp) < 0.02_wp, &
                    "recovered decay-rate error did not halve with dt_bt (not first order)")
      end block checks
   end subroutine test_energy_identity

   subroutine run_successive_singlesteps(r_h, dt_inner, k_steps, u)
      !! Applies K successive single-substep calls (n_steps=1 each) from a
      !! uniform IC and returns u(0:k_steps) at a fixed interior probe.
      !! Uniform r_H + uniform u keeps the field spatially uniform at every
      !! step (no cross-face feedback), so each individual call is exactly
      !! the pure backward-Euler update -- letting us read off an exact
      !! successive trajectory.
      real(wp), intent(in) :: r_h, dt_inner
      integer, intent(in) :: k_steps
      real(wp), intent(out) :: u(0:k_steps)
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ocean_dyn_t) :: dyn
      type(coriolis_adv_t) :: cor
      type(multilayer_state_t) :: ms
      integer :: nx, ny, ip, jp, n

      call build_harness(grid, metrics, dyn, cor, ms)
      nx = grid%nx_total; ny = grid%ny_total
      dyn%bt_work%bt_ubt = U0
      dyn%bt_work%lwd_enable = .true.
      allocate (dyn%bt_work%lwd_drag_u(nx + 1, ny), source=r_h)
      allocate (dyn%bt_work%lwd_drag_v(nx, ny + 1), source=r_h)

      call map_harness(dyn, cor, ms)
      ip = grid%nghost + NX_PHYS/2
      jp = grid%nghost + NY_PHYS/2
      !$acc update self(dyn%bt_work%bt_ubt)
      u(0) = dyn%bt_work%bt_ubt(ip, jp)
      do n = 1, k_steps
         call run_wave_drag_dispatch(grid, metrics, dyn, ms, dt_inner, &
                                     .false., 0.0_wp, 0.0_wp)
         call barotropic_substep_nonlinear_interior(grid, metrics, dyn%bt_work, &
                                                    cor%f_corner, 1, dt_inner)
         ! n_steps=1 -> bt_ubt (time-mean) == bt_ubt_end (final value)
         ! exactly, AND `bt_ubt` itself already holds it (the substep
         ! kernel's own IC for the next call) -- read bt_ubt_end for
         ! consistency with the multi-step tests, nothing to copy back.
         !$acc update self(dyn%bt_work%bt_ubt_end)
         u(n) = dyn%bt_work%bt_ubt_end(ip, jp)
      end do
      call unmap_harness(dyn, cor, ms)
      call teardown_harness(metrics, dyn, cor, ms)
   end subroutine run_successive_singlesteps

   function run_single_step_ratio(r_h, dt_inner) result(ratio)
      !! u1/u0 after ONE substep of size dt_inner, uniform r_H.
      real(wp), intent(in) :: r_h, dt_inner
      real(wp) :: ratio
      real(wp) :: u(0:1)
      call run_successive_singlesteps(r_h, dt_inner, 1, u)
      ratio = u(1)/u(0)
   end function run_single_step_ratio

   ! -----------------------------------------------------------------
   ! Test 8 — roughness_proxy: drag on the ridge flanks, ~0 on the plain,
   ! 0 on land + the ghost ring, respects h2_max.
   ! -----------------------------------------------------------------

   subroutine test_roughness_proxy(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NXP = 60, NYP = 8
      real(wp), parameter :: KAPPA = 6.2832e-4_wp, N_BOT = 1.0e-3_wp
      ! Below the ridge's peak unclamped <h^2> (~1.2e5 for the profile
      ! below) so the h2_max ceiling genuinely fires at the steepest flank
      ! cells without flattening the whole ridge.
      real(wp), parameter :: H2_MAX = 5.0e4_wp
      real(wp), allocatable :: b(:, :), wet_t(:, :), r_h(:, :)
      integer :: nx, ny, i, j, i_crest, i_land
      real(wp) :: plain_max, flank_val, crest_val, ceiling

      checks: block
         nx = NXP + 2*NGHOST
         ny = NYP + 2*NGHOST
         allocate (b(nx, ny), wet_t(nx, ny), r_h(nx, ny))

         ! Flat abyssal plain (positive-up bottom elevation) + a smooth ridge.
         i_crest = NGHOST + NXP/2
         do j = 1, ny
            do i = 1, nx
               b(i, j) = -3000.0_wp + 2500.0_wp* &
                         exp(-((real(i - i_crest, wp)/6.0_wp)**2))
            end do
         end do
         wet_t = 1.0_wp
         ! A land cell squarely on the (otherwise steep) east flank.
         i_land = i_crest + 6
         wet_t(i_land, NGHOST + 2) = 0.0_wp

         call wave_drag_roughness_proxy(b, wet_t, nx, ny, NGHOST, KAPPA, N_BOT, H2_MAX, r_h)

         ! (a) Ghost ring is exactly zero.
         call check(error, maxval(abs(r_h(1:NGHOST, :))) == 0.0_wp, &
                    "roughness_proxy: west ghost ring not zero")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(r_h(nx - NGHOST + 1:nx, :))) == 0.0_wp, &
                    "roughness_proxy: east ghost ring not zero")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(r_h(:, 1:NGHOST))) == 0.0_wp, &
                    "roughness_proxy: south ghost ring not zero")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(r_h(:, ny - NGHOST + 1:ny))) == 0.0_wp, &
                    "roughness_proxy: north ghost ring not zero")
         if (allocated(error)) exit checks

         ! (b) Land cell forced to exactly zero.
         call check(error, r_h(i_land, NGHOST + 2) == 0.0_wp, &
                    "roughness_proxy: land cell not zero")
         if (allocated(error)) exit checks

         ! (c) Far plain (>=20 cells from the crest, well outside the ridge's
         ! Gaussian footprint) is ~0; the flank (steep gradient, a handful
         ! of cells off the crest) exceeds it by many orders of magnitude;
         ! the crest itself (grad b == 0 by symmetry) is smaller still.
         plain_max = max(maxval(abs(r_h(NGHOST + 1:i_crest - 20, NGHOST + 3))), &
                         maxval(abs(r_h(i_crest + 20:nx - NGHOST, NGHOST + 3))))
         flank_val = r_h(i_crest - 6, NGHOST + 3)   ! west flank, wet
         crest_val = r_h(i_crest, NGHOST + 3)
         call check(error, flank_val > 100.0_wp*max(plain_max, 1.0e-30_wp), &
                    "roughness_proxy: flank drag is not >>plain drag")
         if (allocated(error)) exit checks
         call check(error, flank_val > crest_val, &
                    "roughness_proxy: flank drag must exceed the (near-zero-gradient) crest")
         if (allocated(error)) exit checks

         ! (d) h2_max ceiling respected everywhere, AND genuinely fires at
         ! the steepest flank cell (flank_val sits exactly at the ceiling).
         ceiling = 0.5_wp*KAPPA*H2_MAX*N_BOT
         call check(error, maxval(r_h) <= ceiling + 1.0e-18_wp, &
                    "roughness_proxy: h2_max ceiling violated")
         if (allocated(error)) exit checks
         call check(error, abs(flank_val - ceiling) < 1.0e-12_wp, &
                    "roughness_proxy: h2_max ceiling did not fire at the steep flank")
      end block checks
      deallocate (b, wet_t, r_h)
   end subroutine test_roughness_proxy

end module test_ocean_wave_drag
