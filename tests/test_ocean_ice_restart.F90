!! Bit-exact restart roundtrip gate for the ice->ocean EVP momentum
!! mediation (PR 63).  Closes the F4 hole the old configure-time
!! resume-fold routine (deleted by this PR) documented: the blended stress
!! `ocean_surface_stress_t%tau_x/tau_y` the ocean actually consumes did
!! NOT round-trip bit-exactly, because the old resume path RECONSTRUCTED
!! it from the checkpoint's POST-thermo ice concentration `ci`, while the
!! uninterrupted run's blend used the PRE-thermo/pre-transport `ci` of
!! that same step.  This PR carries the blend's own OUTPUT
!! (`ice%tau_ocn_x/y` + the `ice%tau_ocn_valid` presence scalar) instead
!! of recomputing it, so the resume is formula-agnostic by construction.
!!
!! Harness: `roundtrip_ice`, modelled on `test_ocean_restart.F90`'s
!! `roundtrip` (cold-seeded B, `==` not tolerance, the same
!! write/kill/read/advance-both/compare shape), composed with the ice
!! EVP chain the way the driver orders it:
!!   ocean_dyn_step_split (reads tau_x) -> ice_evp_step -> &
!!   ice_ocean_stress_flux (writes tau_x + mirrors tau_ocn_x/y) -> &
!!   [scripted_thermo]  (stand-in for the freeze/melt/transport chain
!!                        that mutates ci AFTER the blend — the causal
!!                        role F4 names, isolated without dragging the
!!                        whole thermo/transport chain into a restart
!!                        test)
!!
!! GPU discipline (mem:separate, commit 72152870's exact hazard on these
!! exact two arrays): every host re-fill gets `!$acc update device`,
!! every host assertion pulls `!$acc update self` on the COMPONENT
!! arrays (`ice%tau_ocn_x`, `ice%tau_ocn_y`) — NEVER the aggregate
!! `update self(ice)` / `update self(stress)`.
module test_ocean_ice_restart
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_config, only: config_t
   use rdb_grid, only: hgrid_t
   use rdb_decomp, only: decomp_t, decomp_init
   use rdb_ocean_metrics, only: metrics_fill_cartesian, metrics_finalize
   use rdb_ocean_state, only: ocean_state_t, ocean_state_enter_data, &
                              ocean_state_exit_data, ocean_state_seed_from_cfg, &
                              ocean_state_restart_write, ocean_state_restart_read
   use rdb_state, only: register_default_tracers
   use rdb_ocean_dyn, only: ocean_dyn_step_split
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t
   use rdb_ice_column, only: ICE_RHO_ICE
   use rdb_ice_evp, only: ice_evp_step, ice_evp_params_t
   use rdb_ice_ocean_coupler, only: ice_ocean_stress_flux, ice_ocean_stress_resume_apply, &
                                    ice_ocean_stress_cleanup
   implicit none
   private

   public :: collect_ocean_ice_restart_tests

   integer, parameter :: NXP = 8, NYP = 6, NZ = 2, NGHOST = 3
   real(wp), parameter :: DX = 2000.0_wp
   real(wp), parameter :: H_TOTAL = 200.0_wp
   real(wp), parameter :: DT = 900.0_wp
   integer, parameter :: N_INNER = 40
   integer, parameter :: N_STEPS = 4
   real(wp), parameter :: TAU_A_X0 = 0.1_wp
   ! Ice patch: interior columns i in [PATCH_I0, PATCH_I1] start iced
   ! (ci=1); the rest of the interior starts ice-free (ci=0).  The u-face
   ! at PATCH_I1+1 (relative to grid, i.e. straddling the last iced
   ! column and the first ice-free one) then has a_u = 1/2 — the blend is
   ! genuinely non-degenerate there (a full-cover or zero-cover patch
   ! would make (1-a)*tau_a + a*fxoc collapse to a fixed point of the OLD
   ! reconstruct fold too, see the risk list item on this).
   integer, parameter :: PATCH_I0 = NGHOST + 1, PATCH_I1 = NGHOST + 4
   ! scripted_thermo melts THIS column fully out (m_ice -> 0), flipping
   ! its ci from 1 to 0 — the "thermo changed ci after the blend" event
   ! F4 is about.  It sits at the iced edge of the patch so the flip
   ! moves the a_u=1/2 interface, touching faces on both sides.
   integer, parameter :: MELT_I = PATCH_I1, MELT_J = NGHOST + 3

contains

   subroutine collect_ocean_ice_restart_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("restart_bit_exact_ice_evp", test_restart_bit_exact_ice_evp), &
                  new_unittest("restart_fresh_run_is_pure_wind", &
                               test_restart_fresh_run_is_pure_wind), &
                  new_unittest("restart_old_checkpoint_degrades_to_wind", &
                               test_restart_old_checkpoint_degrades_to_wind), &
                  new_unittest("ice_dynamics_off_bitident", test_ice_dynamics_off_bitident) &
                  ]
   end subroutine collect_ocean_ice_restart_tests

   subroutine make_cfg(cfg)
      type(config_t), intent(out) :: cfg
      cfg%sim_type = "ocean"
      cfg%nx = NXP
      cfg%ny = NYP
      cfg%dx = DX
      cfg%dy = DX
      cfg%nz_layers = NZ
      cfg%nghost = NGHOST
      cfg%ocean%topo%max_depth = H_TOTAL
      cfg%initial_temperature = 10.0_wp
      cfg%initial_salinity = 34.0_wp
      cfg%ocean%ic%ic_config = "geostrophic_adjustment"
      cfg%ocean%ic%ga_length_scale = 3.0_wp*DX
      cfg%ocean%ic%ga_eta_amp = 0.02_wp
   end subroutine make_cfg

   subroutine seed_ice_patch(ice)
      !! Host-side ice IC: iced columns [PATCH_I0,PATCH_I1] x full
      !! interior j-range, ice-free elsewhere.  ncat=1 lumped mode:
      !! m_ice > 0 => ci = 1 (`ice_cell_concentration_impl`).
      use rdb_ice_state, only: ocean_sea_ice_t
      type(ocean_sea_ice_t), intent(inout) :: ice
      integer :: i, j
      ice%m_ice = 0.0_wp
      do j = NGHOST + 1, NGHOST + NYP
         do i = PATCH_I0, PATCH_I1
            ice%m_ice(i, j, 1) = 2.0_wp*ICE_RHO_ICE
         end do
      end do
   end subroutine seed_ice_patch

   subroutine build_state_ice_cold(cfg, grid, state, dynamics)
      !! Host-only init/seed/metrics — deliberately STOPS SHORT of
      !! enter_data.  Used for the B (resume) side of a roundtrip: the
      !! restart read (host-only) and the resume-apply host ops must
      !! both run BEFORE the state is mapped, mirroring the driver's
      !! real order (test_ocean_restart.F90's B construction is the
      !! model — it does NOT go through the A-side helper either).
      type(config_t), intent(in) :: cfg
      type(hgrid_t), intent(inout) :: grid
      type(ocean_state_t), intent(inout) :: state
      logical, intent(in) :: dynamics

      call grid%init(NXP, NYP, NGHOST, DX, DX)
      state%multilayer%nz_ml = NZ
      state%ice%enable = .true.
      state%ice%ncat = 1
      state%ice%nk_ice = 1
      state%ice%dynamics = dynamics
      call state%init(grid)
      call metrics_fill_cartesian(state%metrics, grid, grid%dx, grid%dy)
      call metrics_finalize(state%metrics)
      call ocean_state_seed_from_cfg(state, grid, cfg)
      call register_default_tracers( &
         state%multilayer%tracers(state%multilayer%idx_salinity), &
         state%multilayer%tracers(state%multilayer%idx_temperature), cfg)

      ! Uniform wind (the "configure_ocean_forcing" analogue): a full
      ! deterministic array overwrite, same on every build/resume — the
      ! restart read never touches tau_x/tau_y (not registered), so B's
      ! pristine wind must match A's for the compare to mean anything.
      state%surface_stress%tau_x = TAU_A_X0
      state%surface_stress%tau_y = 0.0_wp

      if (dynamics) then
         ! Same STATIC ice patch as A (the patch itself is not restart-
         ! carried — matches test_ocean_restart's wetdry "static fields
         ! agree by construction" discipline).  A resume read overwrites
         ! the dynamics PROGNOSTICS (m_ice is NOT registered/read either
         ! — it is IC state, so B needs the identical seed).
         call seed_ice_patch(state%ice)
      end if
   end subroutine build_state_ice_cold

   subroutine setup_state_ice(cfg, grid, state, sf, dynamics)
      !! Mirror the driver's ocean+ice init path far enough to exercise
      !! the resume seam: init, seed, metrics, ice patch + tau_a
      !! snapshot (driver rdb_driver.F90 :1417-1418 order), enter_data.
      !! Fresh-state path only — a restart-read B must use
      !! build_state_ice_cold instead (enter_data must run AFTER the
      !! read + resume apply, never before).
      type(config_t), intent(in) :: cfg
      type(hgrid_t), intent(inout) :: grid
      type(ocean_state_t), intent(inout) :: state
      type(ocean_surface_flux_t), intent(inout) :: sf
      logical, intent(in) :: dynamics

      call build_state_ice_cold(cfg, grid, state, dynamics)
      if (dynamics) then
         ! Driver order (rdb_driver.F90 :1417-1418): snapshot the
         ! PRISTINE wind into tau_a BEFORE any resume apply.
         state%ice%tau_a_x = state%surface_stress%tau_x
         state%ice%tau_a_y = state%surface_stress%tau_y
      end if

      call sf%init(grid)
      call ocean_state_enter_data(state)
      call sf%enter_data()
   end subroutine setup_state_ice

   function default_par() result(par)
      type(ice_evp_params_t) :: par
      par%evp_sub_steps = 32  ! reduced — gate is bit-exactness, not the rheology golden
   end function default_par

   subroutine advance_ice(grid, state, sf, par, mutate_ci)
      !! One outer step in driver order: ocean reads the PRIOR tau_x,
      !! then the ice EVP substep, then the blend (which mirrors its own
      !! output into ice%tau_ocn_x/y), then — optionally — the scripted
      !! thermo mutation that changes ci AFTER the blend.
      type(hgrid_t), intent(in) :: grid
      type(ocean_state_t), intent(inout) :: state
      type(ocean_surface_flux_t), intent(in) :: sf
      type(ice_evp_params_t), intent(in) :: par
      logical, intent(in) :: mutate_ci

      call ocean_dyn_step_split(grid, state%metrics, state%dyn, state%eos, &
                                state%coriolis_adv, state%continuity, &
                                state%pressure_force, state%hvisc, &
                                state%bdrag, state%surface_stress, &
                                state%vert_advect, state%hdiff_tracer, &
                                state%vdiff, state%vmix, &
                                state%multilayer, DT, N_INNER, sf=sf, &
                                vcoord=state%vcoord, bc=state%bc, &
                                lateral_mix=state%lateral_mix, &
                                epbl=state%epbl, kshear=state%kshear)
      ! Mirrors the driver's gate exactly (rdb_driver.F90): ice_evp_step
      ! AND ice_ocean_stress_flux both run only under dynamics=.true. —
      ! with dynamics=.false. the blend must never touch stress%tau_x.
      if (state%ice%dynamics) then
         call ice_evp_step(grid, state%metrics, state%coriolis_adv%f_corner, state%ice, &
                           state%multilayer, DT, par, .false., .false.)
         call ice_ocean_stress_flux(state%metrics, state%surface_stress, state%ice)
      end if
      if (mutate_ci) call scripted_thermo(state%ice)
   end subroutine advance_ice

   subroutine scripted_thermo(ice)
      !! Deterministic stand-in for the freeze/melt/transport chain that
      !! mutates ci AFTER ice_ocean_stress_flux's blend — melts the
      !! MELT_I/MELT_J column fully out (m_ice -> 0), flipping its ci
      !! from 1 to 0 (ncat=1 lumped: ci = merge(1,0, m_ice>0), a scaling
      !! mutation like `m_ice *= 1.05` would NOT change ci in this mode —
      !! risk list item on the two-mode convention).
      use rdb_ice_state, only: ocean_sea_ice_t
      type(ocean_sea_ice_t), intent(inout) :: ice
      ice%m_ice(MELT_I, MELT_J, 1) = 0.0_wp
      !$acc update device(ice%m_ice)
   end subroutine scripted_thermo

   pure function host_ci(m_ice, i, j) result(ci)
      !! ncat=1 lumped-mode ci at one cell, plain host arithmetic (NOT
      !! the shared do-concurrent `ice_cell_concentration_impl` — this is
      !! a host-only precondition probe, and calling a DC kernel on
      !! host-only arrays under -stdpar=gpu is exactly the hazard this PR
      !! removes from the resume path; do not reintroduce it here).
      real(wp), intent(in) :: m_ice(:, :, :)
      integer, intent(in) :: i, j
      real(wp) :: ci
      if (m_ice(i, j, 1) > 0.0_wp) then
         ci = 1.0_wp
      else
         ci = 0.0_wp
      end if
   end function host_ci

   subroutine teardown_ice(state, sf)
      type(ocean_state_t), intent(inout) :: state
      type(ocean_surface_flux_t), intent(inout) :: sf
      call sf%exit_data()
      call ocean_state_exit_data(state)
      call sf%destroy()
      call state%destroy()
      call ice_ocean_stress_cleanup()
   end subroutine teardown_ice

   logical function arrays_identical_2d(a, b) result(same)
      real(wp), intent(in) :: a(:, :), b(:, :)
      integer :: i0, i1, j0, j1
      i0 = NGHOST + 1
      i1 = NGHOST + NXP
      j0 = NGHOST + 1
      j1 = NGHOST + NYP
      same = all(a(i0:i1, j0:j1) == b(i0:i1, j0:j1))
   end function arrays_identical_2d

   logical function arrays_identical_2d_u(a, b) result(same)
      !! u-face compare: interior faces (NGHOST+1 .. NGHOST+NXP+1).
      real(wp), intent(in) :: a(:, :), b(:, :)
      integer :: i0, i1, j0, j1
      i0 = NGHOST + 1
      i1 = NGHOST + NXP + 1
      j0 = NGHOST + 1
      j1 = NGHOST + NYP
      same = all(a(i0:i1, j0:j1) == b(i0:i1, j0:j1))
   end function arrays_identical_2d_u

   logical function arrays_identical_2d_v(a, b) result(same)
      !! v-face compare: interior faces (NGHOST+1 .. NGHOST+NYP+1).
      real(wp), intent(in) :: a(:, :), b(:, :)
      integer :: i0, i1, j0, j1
      i0 = NGHOST + 1
      i1 = NGHOST + NXP
      j0 = NGHOST + 1
      j1 = NGHOST + NYP + 1
      same = all(a(i0:i1, j0:j1) == b(i0:i1, j0:j1))
   end function arrays_identical_2d_v

   logical function arrays_identical_3d(a, b) result(same)
      real(wp), intent(in) :: a(:, :, :), b(:, :, :)
      integer :: i0, i1, j0, j1
      i0 = NGHOST + 1
      i1 = NGHOST + NXP
      j0 = NGHOST + 1
      j1 = NGHOST + NYP
      same = all(a(i0:i1, j0:j1, :) == b(i0:i1, j0:j1, :))
   end function arrays_identical_3d

   subroutine cleanup_file(fn)
      character(len=*), intent(in) :: fn
      integer :: u, ios
      character(len=256) :: msg
      logical :: exists
      inquire (file=fn, exist=exists)
      if (exists) then
         open (newunit=u, file=fn, status="old", action="readwrite", &
               iostat=ios, iomsg=msg)
         if (ios == 0) close (u, status="delete")
      end if
   end subroutine cleanup_file

   ! =====================================================================
   ! Case 1 — the acceptance gate.
   ! =====================================================================

   subroutine test_restart_bit_exact_ice_evp(error)
      !! RESUME §6 #12's "fully bit-exact resume", verbatim.  N_STEPS
      !! outer steps with the LAST one's scripted_thermo mutating ci
      !! after the blend (the checkpoint therefore holds the POST-thermo
      !! ci, while the write-time blend used PRE-thermo ci — exactly F4).
      !! Write, cold-seed B, read, apply, advance both one more step,
      !! compare bit-for-bit.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(hgrid_t) :: gA, gB
      type(ocean_state_t) :: A, B
      type(ocean_surface_flux_t) :: sfA, sfB
      type(decomp_t) :: decomp
      type(ice_evp_params_t) :: par
      character(len=*), parameter :: FN = "test_ocean_ice_restart_rt.nc"
      real(wp) :: t_read, ci_pre, ci_post
      integer :: step_read, s
      logical :: ok
      real(wp), allocatable :: tau_x_a(:, :), tau_y_a(:, :)
      real(wp), allocatable :: tau_ocn_x_a(:, :), tau_ocn_y_a(:, :)

      call make_cfg(cfg)
      par = default_par()
      call decomp_init(decomp, NXP, NYP, 1, 1, 0)
      call cleanup_file(FN)

      call setup_state_ice(cfg, gA, A, sfA, dynamics=.true.)

      do s = 1, N_STEPS - 1
         call advance_ice(gA, A, sfA, par, mutate_ci=.false.)
      end do
      ! Precondition (iii): ci must change AFTER the blend on the write
      ! step, or the OLD reconstruct fold reproduces tau_x exactly and
      ! this case cannot see the bug (test_ocean_restart's wetdry variant
      ! models exactly this "nontrivial at write time" discipline).
      ci_pre = host_ci(A%ice%m_ice, MELT_I, MELT_J)
      call advance_ice(gA, A, sfA, par, mutate_ci=.true.)
      ci_post = host_ci(A%ice%m_ice, MELT_I, MELT_J)
      call check(error, ci_pre > 0.5_wp .and. ci_post < 0.5_wp, &
                 "precondition: scripted_thermo did not flip ci at the melt column")
      if (allocated(error)) then
         call teardown_ice(A, sfA)
         call cleanup_file(FN)
         return
      end if

      ! Precondition (i): a blend was actually written.
      call check(error, A%ice%tau_ocn_valid > 0.5_wp, &
                 "precondition: tau_ocn_valid not set after ice_ocean_stress_flux")
      if (allocated(error)) then
         call teardown_ice(A, sfA)
         call cleanup_file(FN)
         return
      end if
      ! Precondition (ii): the blend is non-degenerate (some interior
      ! face differs from the pristine wind — proves a in (0,1] fired).
      !$acc update self(A%surface_stress%tau_x, A%surface_stress%tau_y)
      call check(error, &
                 maxval(abs(A%surface_stress%tau_x(NGHOST + 1:NGHOST + NXP + 1, &
                                                   NGHOST + 1:NGHOST + NYP) - TAU_A_X0)) &
                 > 0.0_wp, &
                 "precondition: blend degenerate (tau_x == tau_a_x everywhere)")
      if (allocated(error)) then
         call teardown_ice(A, sfA)
         call cleanup_file(FN)
         return
      end if

      ! Pull the device-computed write-time snapshot to the host BEFORE
      ! taking it — component arrays only (never `update self(A%ice)` /
      ! `update self(A%surface_stress)`, commit 72152870).
      !$acc update self(A%surface_stress%tau_x, A%surface_stress%tau_y)
      !$acc update self(A%ice%tau_ocn_x, A%ice%tau_ocn_y)
      allocate (tau_x_a, source=A%surface_stress%tau_x)
      allocate (tau_y_a, source=A%surface_stress%tau_y)
      allocate (tau_ocn_x_a, source=A%ice%tau_ocn_x)
      allocate (tau_ocn_y_a, source=A%ice%tau_ocn_y)

      call ocean_state_restart_write(A, gA, decomp, FN, real(N_STEPS, wp)*DT, N_STEPS)

      ! State B: cold seed (same static IC as A, INCLUDING the ice patch
      ! — the patch itself is not restart-carried, matching the wetdry
      ! variant's "static fields agree by construction" discipline),
      ! HOST-ONLY (build_state_ice_cold — no enter_data yet), then read
      ! (overwrites the ice-dynamics prognostics), mirror the driver's
      ! configure-time tau_a snapshot + resume apply, THEN enter_data
      ! (mirrors test_ocean_restart.F90's B construction — read runs on
      ! the host BEFORE enter_data).
      call build_state_ice_cold(cfg, gB, B, dynamics=.true.)
      call ocean_state_restart_read(B, gB, decomp, FN, t_read, step_read)
      call check(error, step_read == N_STEPS, "restart step count mismatch")
      if (allocated(error)) return

      ! Direct round-trip: B's freshly-read tau_ocn_x/y must equal A's
      ! write-time snapshot, and tau_ocn_valid must have come back 1.
      call check(error, arrays_identical_2d_u(tau_ocn_x_a, B%ice%tau_ocn_x), &
                 "restored tau_ocn_x /= A's write-time snapshot")
      if (allocated(error)) return
      call check(error, arrays_identical_2d_v(tau_ocn_y_a, B%ice%tau_ocn_y), &
                 "restored tau_ocn_y /= A's write-time snapshot")
      if (allocated(error)) return
      call check(error, B%ice%tau_ocn_valid > 0.5_wp, &
                 "restored tau_ocn_valid /= 1.0")
      if (allocated(error)) return

      ! Mirror the driver's configure sequence: pristine-wind snapshot
      ! FIRST (B's surface_stress%tau_x is still the uniform wind — the
      ! read never touches it, tau_x is not registered), then the apply.
      B%ice%tau_a_x = B%surface_stress%tau_x
      B%ice%tau_a_y = B%surface_stress%tau_y
      call ice_ocean_stress_resume_apply(B%surface_stress, B%ice)

      ! After the apply: bit-identical to A's write-time stress.
      call check(error, arrays_identical_2d_u(tau_x_a, B%surface_stress%tau_x), &
                 "post-apply stress%tau_x /= A's write-time snapshot")
      if (allocated(error)) return
      call check(error, arrays_identical_2d_v(tau_y_a, B%surface_stress%tau_y), &
                 "post-apply stress%tau_y /= A's write-time snapshot")
      if (allocated(error)) return

      ! NOW map — copyin picks up the just-restored + just-applied host
      ! state (h_layer, tracers, ice%str_*/u_ice/v_ice/fxoc/fyoc from the
      ! read; surface_stress%tau_x/y and ice%tau_a_x/y from the apply
      ! above), exactly once.
      call sfB%init(gB)
      call ocean_state_enter_data(B)
      call sfB%enter_data()

      ! Advance both one more step (mutate_ci=.false. — the payload here
      ! is step (N+1), not a second thermo event) and compare bit-for-bit.
      call advance_ice(gA, A, sfA, par, mutate_ci=.false.)
      call advance_ice(gB, B, sfB, par, mutate_ci=.false.)

      !$acc update self(A%surface_stress%tau_x, A%surface_stress%tau_y)
      !$acc update self(B%surface_stress%tau_x, B%surface_stress%tau_y)
      !$acc update self(A%ice%u_ice, A%ice%v_ice, A%ice%str_d, A%ice%str_t, A%ice%str_s)
      !$acc update self(B%ice%u_ice, B%ice%v_ice, B%ice%str_d, B%ice%str_t, B%ice%str_s)
      !$acc update self(A%ice%fxoc, A%ice%fyoc, B%ice%fxoc, B%ice%fyoc)
      !$acc update self(A%multilayer%h_layer, A%multilayer%u_face_x_layer, &
      !$acc&            A%multilayer%v_face_y_layer, A%barotropic%h, &
      !$acc&            A%barotropic%u_face_x, A%barotropic%v_face_y)
      !$acc update self(B%multilayer%h_layer, B%multilayer%u_face_x_layer, &
      !$acc&            B%multilayer%v_face_y_layer, B%barotropic%h, &
      !$acc&            B%barotropic%u_face_x, B%barotropic%v_face_y)

      compare: block
         integer :: it
         ok = arrays_identical_2d_u(A%surface_stress%tau_x, B%surface_stress%tau_x)
         call check(error, ok, "step (N+1): stress%tau_x not bit-identical")
         if (allocated(error)) exit compare
         ok = arrays_identical_2d_v(A%surface_stress%tau_y, B%surface_stress%tau_y)
         call check(error, ok, "step (N+1): stress%tau_y not bit-identical")
         if (allocated(error)) exit compare
         ok = arrays_identical_2d_u(A%ice%u_ice, B%ice%u_ice)
         call check(error, ok, "step (N+1): ice%u_ice not bit-identical")
         if (allocated(error)) exit compare
         ok = arrays_identical_2d_v(A%ice%v_ice, B%ice%v_ice)
         call check(error, ok, "step (N+1): ice%v_ice not bit-identical")
         if (allocated(error)) exit compare
         ok = arrays_identical_2d(A%ice%str_d, B%ice%str_d)
         call check(error, ok, "step (N+1): ice%str_d not bit-identical")
         if (allocated(error)) exit compare
         ok = arrays_identical_2d(A%ice%str_t, B%ice%str_t)
         call check(error, ok, "step (N+1): ice%str_t not bit-identical")
         if (allocated(error)) exit compare
         ok = all(A%ice%str_s(NGHOST + 1:NGHOST + NXP + 1, NGHOST + 1:NGHOST + NYP + 1) == &
                  B%ice%str_s(NGHOST + 1:NGHOST + NXP + 1, NGHOST + 1:NGHOST + NYP + 1))
         call check(error, ok, "step (N+1): ice%str_s not bit-identical")
         if (allocated(error)) exit compare
         ok = arrays_identical_2d_u(A%ice%fxoc, B%ice%fxoc)
         call check(error, ok, "step (N+1): ice%fxoc not bit-identical")
         if (allocated(error)) exit compare
         ok = arrays_identical_2d_v(A%ice%fyoc, B%ice%fyoc)
         call check(error, ok, "step (N+1): ice%fyoc not bit-identical")
         if (allocated(error)) exit compare

         ok = arrays_identical_3d(A%multilayer%h_layer, B%multilayer%h_layer)
         call check(error, ok, "step (N+1): h_layer not bit-identical")
         if (allocated(error)) exit compare
         ok = arrays_identical_3d(A%multilayer%u_face_x_layer, B%multilayer%u_face_x_layer)
         call check(error, ok, "step (N+1): u_face_x_layer not bit-identical")
         if (allocated(error)) exit compare
         ok = arrays_identical_3d(A%multilayer%v_face_y_layer, B%multilayer%v_face_y_layer)
         call check(error, ok, "step (N+1): v_face_y_layer not bit-identical")
         if (allocated(error)) exit compare
         ok = arrays_identical_2d(A%barotropic%h, B%barotropic%h)
         call check(error, ok, "step (N+1): barotropic h not bit-identical")
         if (allocated(error)) exit compare

         do it = 1, size(A%multilayer%tracers)
            !$acc update self(A%multilayer%tracers(it)%hTr, B%multilayer%tracers(it)%hTr)
            ok = arrays_identical_3d(A%multilayer%tracers(it)%hTr, &
                                     B%multilayer%tracers(it)%hTr)
            call check(error, ok, "step (N+1): tracer hTr not bit-identical")
            if (allocated(error)) exit
         end do
      end block compare

      call teardown_ice(A, sfA)
      call teardown_ice(B, sfB)
      call cleanup_file(FN)
   end subroutine test_restart_bit_exact_ice_evp

   ! =====================================================================
   ! Case 2 — fresh run is pure wind.
   ! =====================================================================

   subroutine test_restart_fresh_run_is_pure_wind(error)
      !! No restart file: tau_ocn_valid stays 0.0 (the type default), and
      !! ice_ocean_stress_resume_apply must leave stress%tau_x/y BIT-
      !! IDENTICAL to the pristine wind — guards the catastrophic trap
      !! (§11.2 of the plan): an unconditional copy would hand the ocean
      !! a ZEROED wind stress on step 1 of every fresh ice-dynamics run.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      type(ocean_surface_flux_t) :: sf
      real(wp), allocatable :: tau_x0(:, :), tau_y0(:, :)

      call make_cfg(cfg)
      call setup_state_ice(cfg, grid, state, sf, dynamics=.true.)

      call check(error, state%ice%tau_ocn_valid < 0.5_wp, &
                 "fresh run: tau_ocn_valid must be 0.0 before any blend")
      if (allocated(error)) then
         call teardown_ice(state, sf)
         return
      end if

      allocate (tau_x0, source=state%surface_stress%tau_x)
      allocate (tau_y0, source=state%surface_stress%tau_y)
      call ice_ocean_stress_resume_apply(state%surface_stress, state%ice)
      call check(error, all(state%surface_stress%tau_x == tau_x0), &
                 "fresh run: resume apply must NOT touch stress%tau_x (tau_ocn_valid==0)")
      if (allocated(error)) then
         call teardown_ice(state, sf)
         return
      end if
      call check(error, all(state%surface_stress%tau_y == tau_y0), &
                 "fresh run: resume apply must NOT touch stress%tau_y (tau_ocn_valid==0)")

      call teardown_ice(state, sf)
   end subroutine test_restart_fresh_run_is_pure_wind

   ! =====================================================================
   ! Case 3 — a pre-PR-63 checkpoint degrades to pure wind.
   ! =====================================================================

   subroutine test_restart_old_checkpoint_degrades_to_wind(error)
      !! Write a checkpoint from a state whose ice_ocean_stress_flux was
      !! NEVER called (dynamics=.true. but no outer step ran the blend —
      !! the simplest realisation of "a pre-PR-63 file": tau_ocn_x/y are
      !! registered but hold their init value 0, tau_ocn_valid==0.0).
      !! Reading it into a live B must leave B%ice%tau_ocn_valid==0.0 and
      !! the resume apply must leave the pristine wind — no fatal, no
      !! silently-zeroed stress.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(hgrid_t) :: gA, gB
      type(ocean_state_t) :: A, B
      type(ocean_surface_flux_t) :: sfA, sfB
      type(decomp_t) :: decomp
      character(len=*), parameter :: FN = "test_ocean_ice_restart_old.nc"
      real(wp) :: t_read
      integer :: step_read
      real(wp), allocatable :: tau_x0(:, :), tau_y0(:, :)

      call make_cfg(cfg)
      call decomp_init(decomp, NXP, NYP, 1, 1, 0)
      call cleanup_file(FN)

      call setup_state_ice(cfg, gA, A, sfA, dynamics=.true.)
      ! Deliberately do NOT call advance_ice / ice_ocean_stress_flux — A's
      ! tau_ocn_x/y/valid stay at their init value (0/0/0.0), exactly what
      ! a resume from a checkpoint written before this PR existed would
      ! look like once read back through the new optional=.true. fields.
      call ocean_state_restart_write(A, gA, decomp, FN, 0.0_wp, 0)
      call teardown_ice(A, sfA)

      ! B: host-only cold build (mirrors case 1 — enter_data must run
      ! AFTER the read + resume apply, never before).
      call build_state_ice_cold(cfg, gB, B, dynamics=.true.)
      call ocean_state_restart_read(B, gB, decomp, FN, t_read, step_read)

      call check(error, B%ice%tau_ocn_valid < 0.5_wp, &
                 "old checkpoint: tau_ocn_valid must read back 0.0")
      if (allocated(error)) then
         call B%destroy()
         call cleanup_file(FN)
         return
      end if

      allocate (tau_x0, source=B%surface_stress%tau_x)
      allocate (tau_y0, source=B%surface_stress%tau_y)
      B%ice%tau_a_x = B%surface_stress%tau_x
      B%ice%tau_a_y = B%surface_stress%tau_y
      call ice_ocean_stress_resume_apply(B%surface_stress, B%ice)
      call check(error, all(B%surface_stress%tau_x == tau_x0), &
                 "old checkpoint: resume must leave pristine wind (tau_x)")
      if (.not. allocated(error)) then
         call check(error, all(B%surface_stress%tau_y == tau_y0), &
                    "old checkpoint: resume must leave pristine wind (tau_y)")
      end if

      call sfB%init(gB)
      call ocean_state_enter_data(B)
      call sfB%enter_data()
      call teardown_ice(B, sfB)
      call cleanup_file(FN)
   end subroutine test_restart_old_checkpoint_degrades_to_wind

   ! =====================================================================
   ! Case 4 — dynamics=.false. default-off bit-identity, localised.
   ! =====================================================================

   subroutine test_ice_dynamics_off_bitident(error)
      !! enable=.true., dynamics=.false.: ice_ocean_stress_flux never
      !! runs (driver gate), tau_ocn_x/y stay exactly 0, tau_ocn_valid
      !! stays 0.0, and a full write/read/step-(N+1) round-trip leaves
      !! every ocean prognostic bit-identical.  Localises the intent of
      !! the "default off ⇒ bit-identical" gate (the real proof is the
      !! full `ctest -R rdb` byte-identity).
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(hgrid_t) :: gA, gB
      type(ocean_state_t) :: A, B
      type(ocean_surface_flux_t) :: sfA, sfB
      type(decomp_t) :: decomp
      type(ice_evp_params_t) :: par
      character(len=*), parameter :: FN = "test_ocean_ice_restart_off.nc"
      real(wp) :: t_read
      integer :: step_read, s
      logical :: ok

      call make_cfg(cfg)
      par = default_par()
      call decomp_init(decomp, NXP, NYP, 1, 1, 0)
      call cleanup_file(FN)

      call setup_state_ice(cfg, gA, A, sfA, dynamics=.false.)
      do s = 1, N_STEPS
         call advance_ice(gA, A, sfA, par, mutate_ci=.false.)
      end do

      !$acc update self(A%ice%tau_ocn_x, A%ice%tau_ocn_y)
      call check(error, all(A%ice%tau_ocn_x == 0.0_wp) .and. all(A%ice%tau_ocn_y == 0.0_wp), &
                 "dynamics=.false.: tau_ocn_x/y must stay exactly 0")
      if (allocated(error)) then
         call teardown_ice(A, sfA)
         call cleanup_file(FN)
         return
      end if
      call check(error, A%ice%tau_ocn_valid < 0.5_wp, &
                 "dynamics=.false.: tau_ocn_valid must stay 0.0")
      if (allocated(error)) then
         call teardown_ice(A, sfA)
         call cleanup_file(FN)
         return
      end if

      call ocean_state_restart_write(A, gA, decomp, FN, real(N_STEPS, wp)*DT, N_STEPS)

      call build_state_ice_cold(cfg, gB, B, dynamics=.false.)
      call ocean_state_restart_read(B, gB, decomp, FN, t_read, step_read)
      B%ice%tau_a_x = B%surface_stress%tau_x
      B%ice%tau_a_y = B%surface_stress%tau_y
      call ice_ocean_stress_resume_apply(B%surface_stress, B%ice)
      call sfB%init(gB)
      call ocean_state_enter_data(B)
      call sfB%enter_data()

      call advance_ice(gA, A, sfA, par, mutate_ci=.false.)
      call advance_ice(gB, B, sfB, par, mutate_ci=.false.)

      !$acc update self(A%multilayer%h_layer, B%multilayer%h_layer)
      !$acc update self(A%barotropic%h, B%barotropic%h)
      ok = arrays_identical_3d(A%multilayer%h_layer, B%multilayer%h_layer)
      call check(error, ok, "dynamics=.false.: h_layer not bit-identical")
      if (.not. allocated(error)) then
         ok = arrays_identical_2d(A%barotropic%h, B%barotropic%h)
         call check(error, ok, "dynamics=.false.: barotropic h not bit-identical")
      end if

      call teardown_ice(A, sfA)
      call teardown_ice(B, sfB)
      call cleanup_file(FN)
   end subroutine test_ice_dynamics_off_bitident

end module test_ocean_ice_restart
