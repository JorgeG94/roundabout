!! Windowed tracer advection (`&ocean_vmix_nml dt_tracer_advect_ratio > 1`)
!! under the default `split_scheme = "pred_corr"`, driven through the
!! PRODUCTION engine path (`engine_step` / `engine_step_finalize`).
!!
!! Under `pred_corr` the predictor runs the continuity chain with
!! `TR_MODE_NONE` and only the corrector accumulates the face transports
!! `uhtr/vhtr`.  The accumulator therefore has to carry the corrector's
!! transport at weight 1.0 (the `ssp_rk2` stage weight 0.5 would leave it
!! holding half of what moved `h`), and the drain's post-average budget
!! weight has to follow the scheme (`1/ocean_budget_stage_weight`).
!!
!!   (1) `uniform_tracer_stays_uniform` -- a uniform S and T, a wind-driven
!!       flow that genuinely moves `h`, ratio = 4, three windows (12 steps).  The drain
!!       reconstructs `hprev = h_end + div(uhtr)`; if `uhtr` holds half the
!!       transport, `hprev` is wrong by half the window divergence and the
!!       concentration `hTr/hprev` is biased by `O(dh/h)` -- orders of
!!       magnitude above the 1e-13 gate.
!!   (2) `conservation_over_windows` -- `sum(h*T)` and `sum(h*S)` over a
!!       closed basin with a stratified initial state, to 1e-12.
!!   (3) `ratio4_vs_ratio1` -- the same span with ratio 4 and ratio 1
!!       (both with `dt_therm_ratio = 4`, so only the advection cadence
!!       differs).  The two are different time discretisations of the same
!!       operator: ratio 4 advects over a 4*dt window with the window-mean
!!       transport, a first-order-in-window truncation difference.  The
!!       gate is a fraction of the tracer's own excursion, and the
!!       difference is asserted NON-zero (the windows really differ).
!!   (4) `restart_at_window_boundary` -- checkpoint at a window boundary
!!       (after the driver's flush) and warm restart is bitwise identical
!!       to the uninterrupted run, every restart-registry field.
module test_ocean_tracer_ratio_pc
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_constants, only: wp
   use rdb_config, only: config_t, read_config_from_string, validate_config
   use rdb_ocean_engine, only: ocean_engine_t, engine_setup, engine_enter_data, &
                               engine_step, engine_step_ice, engine_step_finalize, &
                               engine_exit_data, engine_teardown
   use rdb_ocean_state, only: ocean_state_restart_write, ocean_state_build_restart_registry
   use rdb_ocean_restart, only: restart_registry_t
   use rdb_ocean_dyn, only: ocean_dyn_flush_tracer_window, ocean_dyn_t
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   use rdb_comm_env, only: comm_env_init, comm_env_setup_roles
   implicit none
   private

   public :: collect_ocean_tracer_ratio_pc_tests

   real(wp), parameter :: DT = 600.0_wp
   character(len=*), parameter :: FN = "test_ocean_tracer_ratio_pc_rt.nc"
   integer, parameter :: NSTEP = 12
      !! Three 4-step windows: the drain closes a window on steps 4, 8, 12.
   logical :: comm_inited = .false.

contains

   subroutine collect_ocean_tracer_ratio_pc_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("uniform_tracer_stays_uniform", test_uniform), &
                  new_unittest("conservation_over_windows", test_conservation), &
                  new_unittest("ratio4_vs_ratio1", test_ratio_vs_one), &
                  new_unittest("restart_at_window_boundary", test_restart), &
                  new_unittest("window_cadence_predicates", test_cadence), &
                  new_unittest("window_length_exact", test_window_length) &
                  ]
   end subroutine collect_ocean_tracer_ratio_pc_tests

   subroutine ensure_comm()
      if (.not. comm_inited) then
         call comm_env_init()
         call comm_env_setup_roles(.false.)
         comm_inited = .true.
      end if
   end subroutine ensure_comm

   function case_nml(ratio, uniform, vcoord) result(nml)
      !! Closed wind-driven basin with an island (divergent flow, real `h`
      !! motion), `pred_corr` (default), `dt_therm_ratio = 4`.
      integer, intent(in) :: ratio
      logical, intent(in) :: uniform
      character(len=*), intent(in), optional :: vcoord
      character(len=:), allocatable :: nml
      character(len=:), allocatable :: vc
      character(len=*), parameter :: NL = new_line("a")
      character(len=16) :: rs
      character(len=:), allocatable :: tr
      write (rs, '(i0)') ratio
      vc = 'zstar'
      if (present(vcoord)) vc = vcoord
      if (uniform) then
         tr = "&tracer_nml initial_temperature = 12.0, initial_salinity = 35.0 /"//NL
      else
         tr = "&tracer_nml initial_temperature = 12.0, initial_salinity = 35.0, "// &
              "T_init_surface = 20.0, T_init_bottom = 4.0 /"//NL
      end if
      nml = "&sim_nml sim_type = 'ocean' /"//NL// &
            "&time_nml t_end = 86400.0, dt_fixed = 600.0 /"//NL// &
            "&nonhydrostatic_nml nz_layers = 6 /"//NL// &
            tr// &
            "&ocean_bt_nml auto_n_inner = .true. /"//NL// &
            "&ocean_hvisc_nml nu_h = 200.0, lateral_closure = 'smagorinsky', "// &
            "smag_ah = .true. /"//NL// &
            "&ocean_vmix_nml dt_therm_ratio = 4, dt_tracer_advect_ratio = "//trim(rs)//" /"//NL// &
            "&ocean_diag_nml enabled = .false. /"//NL// &
            "&grid_nml nx = 24, ny = 16, nghost = 3, dx = 10000.0, dy = 10000.0 /"//NL// &
            "&physics_nml coriolis_f = 1.0e-4, wind_stress_x = 0.08 /"//NL// &
            "&vcoord_nml vcoord_type = '"//vc//"' /"//NL// &
            "&ocean_topo_nml topo_config = 'island', max_depth = 1000.0, "// &
            "slope_scale = 0.25 /"//NL// &
            "&ocean_bc_nml west = 'wall', east = 'wall', south = 'wall', "// &
            "north = 'wall' /"//NL// &
            "&output_nml output_to_file = .false. /"//NL
   end function case_nml

   subroutine make_engine(engine, cfg, nml, ok, restart_file, t0, step0)
      type(ocean_engine_t), intent(inout) :: engine
      type(config_t), intent(inout) :: cfg
      character(len=*), intent(in) :: nml
      logical, intent(out) :: ok
      character(len=*), intent(in), optional :: restart_file
      real(wp), intent(out), optional :: t0
      integer, intent(out), optional :: step0
      integer :: ierr
      ok = .false.
      call ensure_comm()
      call read_config_from_string(nml, cfg, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call validate_config(cfg, ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      if (present(restart_file)) then
         call engine_setup(engine, cfg, ierr, restart_file=restart_file, &
                           t_restart=t0, step_restart=step0)
      else
         call engine_setup(engine, cfg, ierr)
      end if
      if (ierr /= OCEAN_STATUS_OK) return
      call engine_enter_data(engine, cfg)
      ok = .true.
   end subroutine make_engine

   subroutine advance(engine, cfg, t, nsteps, ok)
      type(ocean_engine_t), intent(inout) :: engine
      type(config_t), intent(inout) :: cfg
      real(wp), intent(inout) :: t
      integer, intent(in) :: nsteps
      logical, intent(out) :: ok
      integer :: n, ierr
      ok = .false.
      do n = 1, nsteps
         call engine_step(engine, DT, t, ierr=ierr)
         if (ierr /= OCEAN_STATUS_OK) return
         call engine_step_ice(engine, cfg, DT, t, ierr=ierr)
         if (ierr /= OCEAN_STATUS_OK) return
         call engine_step_finalize(engine, DT, t, ierr=ierr)
         if (ierr /= OCEAN_STATUS_OK) return
         t = t + DT
      end do
      ok = .true.
   end subroutine advance

   subroutine pull(engine, h, st, tt)
      !! Host copies of h_layer, S and T concentration-content (device ->
      !! host on COMPONENT arrays only), full local extent.
      type(ocean_engine_t), intent(inout), target :: engine
      real(wp), allocatable, intent(out) :: h(:, :, :), st(:, :, :), tt(:, :, :)
      integer :: is, it
      is = engine%state%multilayer%idx_salinity
      it = engine%state%multilayer%idx_temperature
      !$acc update self(engine%state%multilayer%h_layer)
      associate (hs => engine%state%multilayer%tracers(is)%hTr, &
                 ht => engine%state%multilayer%tracers(it)%hTr)
         !$acc update self(hs, ht)
      end associate
      h = engine%state%multilayer%h_layer
      st = engine%state%multilayer%tracers(is)%hTr
      tt = engine%state%multilayer%tracers(it)%hTr
   end subroutine pull

   subroutine test_uniform(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: e
      type(config_t) :: cfg
      real(wp), allocatable :: h(:, :, :), hs(:, :, :), ht(:, :, :), h0(:, :, :)
      real(wp) :: t, err_s, err_t, dh
      integer :: ng, i, j, k
      logical :: ok
      call make_engine(e, cfg, case_nml(4, .true.), ok)
      call check(error, ok, "engine setup failed")
      if (allocated(error)) return
      ng = e%grid%nghost
      call pull(e, h0, hs, ht)
      t = 0.0_wp
      call advance(e, cfg, t, NSTEP - 2, ok)
      ! mid-window sample: proves the flow really moves h
      call pull(e, h, hs, ht)
      dh = 0.0_wp
      do k = 1, size(h, 3)
         do j = ng + 1, size(h, 2) - ng
            do i = ng + 1, size(h, 1) - ng
               if (h(i, j, k) < 1.0e-2_wp) cycle
               dh = max(dh, abs(h(i, j, k) - h0(i, j, k))/max(h0(i, j, k), 1.0e-2_wp))
            end do
         end do
      end do
      if (ok) call advance(e, cfg, t, 2, ok)
      call check(error, ok, "engine failed during the run")
      if (.not. allocated(error)) then
         call pull(e, h, hs, ht)
         err_s = 0.0_wp
         err_t = 0.0_wp
         do k = 1, size(h, 3)
            do j = ng + 1, size(h, 2) - ng
               do i = ng + 1, size(h, 1) - ng
                  if (h(i, j, k) < 1.0e-2_wp) cycle   ! massless / land layers
                  err_s = max(err_s, abs(hs(i, j, k)/h(i, j, k) - 35.0_wp))
                  err_t = max(err_t, abs(ht(i, j, k)/h(i, j, k) - 12.0_wp))
               end do
            end do
         end do
         call check(error, dh > 1.0e-4_wp, "the flow did not move h: test would pass vacuously")
         if (.not. allocated(error)) &
            call check(error, err_s < 1.0e-13_wp*35.0_wp .and. err_t < 1.0e-13_wp*12.0_wp, &
                       "uniform tracer left uniform under pred_corr ratio 4 (accumulator weight)")
         if (.not. allocated(error)) then
            write (*, '(a,es10.2,a,es10.2,a,es10.2)') "  uniform: max|dS|=", err_s, &
               " max|dT|=", err_t, " max rel dh=", dh
         end if
      end if
      call engine_exit_data(e)
      call engine_teardown(e)
   end subroutine test_uniform

   subroutine total_content(engine, s, t_)
      type(ocean_engine_t), intent(inout), target :: engine
      real(wp), intent(out) :: s, t_
      real(wp), allocatable :: h(:, :, :), hs(:, :, :), ht(:, :, :)
      integer :: ng, nx, ny
      call pull(engine, h, hs, ht)
      ng = engine%grid%nghost
      nx = size(h, 1)
      ny = size(h, 2)
      s = sum(hs(ng + 1:nx - ng, ng + 1:ny - ng, :))
      t_ = sum(ht(ng + 1:nx - ng, ng + 1:ny - ng, :))
   end subroutine total_content

   subroutine test_conservation(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: r
      do r = 1, 4, 3
         call conservation_at(error, r)
         if (allocated(error)) return
      end do
   end subroutine test_conservation

   subroutine conservation_at(error, ratio)
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: ratio
      type(ocean_engine_t), target :: e
      type(config_t) :: cfg
      real(wp) :: s0, t0c, s1, t1c, t
      logical :: ok
      call make_engine(e, cfg, case_nml(ratio, .false.), ok)
      call check(error, ok, "engine setup failed")
      if (allocated(error)) return
      call total_content(e, s0, t0c)
      t = 0.0_wp
      call advance(e, cfg, t, NSTEP, ok)
      call check(error, ok, "engine failed during the run")
      if (.not. allocated(error)) then
         call total_content(e, s1, t1c)
         write (*, '(a,i0,a,es10.2,a,es10.2)') "  conservation ratio ", ratio, ": rel d(sum hS)=", abs(s1 - s0)/s0, &
            " rel d(sum hT)=", abs(t1c - t0c)/t0c
         call check(error, abs(s1 - s0) <= 1.0e-12_wp*abs(s0), "sum(h*S) not conserved over 3 windows")
         ! T is stratified (S is uniform, so it cannot see a hold/un-hold
         ! mismatch): the sharper of the two.
         if (.not. allocated(error)) &
            call check(error, abs(t1c - t0c) <= 1.0e-12_wp*abs(t0c), "sum(h*T) not conserved over 3 windows")
      end if
      call engine_exit_data(e)
      call engine_teardown(e)
   end subroutine conservation_at

   subroutine test_ratio_vs_one(error)
      !! Bound: ratio 4 and ratio 1 differ by the time truncation of
      !! advecting with the window-mean transport over 4*dt instead of dt.
      !! That error scales with (window/advective time) times the tracer
      !! excursion, so the gate is a FRACTION of the excursion
      !! max|T - T0| the run itself produced, not an absolute number.
      !! Gate 10%; measured 3% (2.3e-3 of 7.1e-2), identical under ssp_rk2:
      !! it is the splitting difference against the 4*dt vertical mixing.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: e1, e4
      type(config_t) :: c1, c4
      real(wp), allocatable :: h1(:, :, :), s1(:, :, :), t1(:, :, :)
      real(wp), allocatable :: h4(:, :, :), s4(:, :, :), t4(:, :, :), h0(:, :, :), s0(:, :, :), t0a(:, :, :)
      real(wp) :: ta, tb, dmax, excur, c1v, c4v
      integer :: ng, i, j, k
      logical :: ok
      call make_engine(e1, c1, case_nml(1, .false., 'sigma'), ok)
      call check(error, ok, "ratio-1 engine setup failed")
      if (allocated(error)) return
      call make_engine(e4, c4, case_nml(4, .false., 'sigma'), ok)
      call check(error, ok, "ratio-4 engine setup failed")
      if (allocated(error)) return
      ng = e1%grid%nghost
      call pull(e1, h0, s0, t0a)
      ta = 0.0_wp
      tb = 0.0_wp
      call advance(e1, c1, ta, NSTEP, ok)
      call check(error, ok, "ratio-1 run failed")
      call advance(e4, c4, tb, NSTEP, ok)
      call check(error, ok, "ratio-4 run failed")
      if (.not. allocated(error)) then
         call pull(e1, h1, s1, t1)
         call pull(e4, h4, s4, t4)
         dmax = 0.0_wp
         excur = 0.0_wp
         do k = 1, size(h1, 3)
            do j = ng + 1, size(h1, 2) - ng
               do i = ng + 1, size(h1, 1) - ng
                  if (h1(i, j, k) < 1.0e-2_wp .or. h4(i, j, k) < 1.0e-2_wp) cycle
                  c1v = t1(i, j, k)/h1(i, j, k)
                  c4v = t4(i, j, k)/h4(i, j, k)
                  dmax = max(dmax, abs(c1v - c4v))
                  excur = max(excur, abs(c1v - t0a(i, j, k)/h0(i, j, k)))
               end do
            end do
         end do
         write (*, '(a,es10.2,a,es10.2)') "  ratio4-vs-1: max|dT|=", dmax, " excursion=", excur
         call check(error, excur > 1.0e-6_wp, "tracer did not move: comparison is vacuous")
         if (.not. allocated(error)) &
            call check(error, dmax > 0.0_wp, "ratio 4 and ratio 1 identical: windows not exercised")
         if (.not. allocated(error)) &
            call check(error, dmax < 0.1_wp*excur, &
                       "ratio 4 departs from ratio 1 by more than the truncation bound")
      end if
      call engine_exit_data(e1)
      call engine_teardown(e1)
      call engine_exit_data(e4)
      call engine_teardown(e4)
   end subroutine test_ratio_vs_one

   subroutine compare_registries(error, ea, eb, when)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), intent(inout), target :: ea, eb
      character(len=*), intent(in) :: when
      type(restart_registry_t) :: ra, rb
      integer :: e, nbad, nfield
      character(len=:), allocatable :: msg
      call ocean_state_build_restart_registry(ea%state, ea%grid, ra)
      call ocean_state_build_restart_registry(eb%state, eb%grid, rb)
      call check(error, ra%n == rb%n, "registry size differs")
      if (allocated(error)) return
      nbad = 0
      nfield = 0
      msg = ""
      do e = 1, ra%n
         associate (a => ra%entries(e), b => rb%entries(e))
            if (a%device_mapped .and. a%rank == 2) then
               !$acc update self(a%p2)
               !$acc update self(b%p2)
            else if (a%device_mapped .and. a%rank == 3) then
               !$acc update self(a%p3)
               !$acc update self(b%p3)
            end if
            nfield = nfield + 1
            select case (a%rank)
            case (0)
               if (transfer(a%p0, 0_int64) /= transfer(b%p0, 0_int64)) then
                  nbad = nbad + 1
                  msg = msg//" "//trim(a%tag)
               end if
            case (2)
               if (.not. same_bits2(a%p2, b%p2)) then
                  nbad = nbad + 1
                  msg = msg//" "//trim(a%tag)
               end if
            case default
               if (.not. same_bits3(a%p3, b%p3)) then
                  nbad = nbad + 1
                  msg = msg//" "//trim(a%tag)
               end if
            end select
         end associate
      end do
      call check(error, nfield > 20, "too few registry fields compared")
      if (allocated(error)) return
      call check(error, nbad == 0, "restart not bit-identical "//when//":"//msg)
   end subroutine compare_registries

   pure function same_bits2(a, b) result(same)
      real(wp), intent(in) :: a(:, :), b(:, :)
      logical :: same
      same = all(shape(a) == shape(b))
      if (same) same = all(transfer(a, [0_int64]) == transfer(b, [0_int64]))
   end function same_bits2

   pure function same_bits3(a, b) result(same)
      real(wp), intent(in) :: a(:, :, :), b(:, :, :)
      logical :: same
      same = all(shape(a) == shape(b))
      if (same) same = all(transfer(a, [0_int64]) == transfer(b, [0_int64]))
   end function same_bits3

   subroutine test_restart(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: r
      do r = 1, 4, 3
         call restart_at(error, r)
         if (allocated(error)) return
      end do
   end subroutine test_restart

   subroutine restart_at(error, ratio)
      !! Straight 12 steps vs 8 steps + (flush, checkpoint) + warm restart +
      !! 4 steps.  8 = two whole windows, so the checkpoint is window-aligned.
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: ratio
      type(ocean_engine_t), target :: ea, eb, ec
      type(config_t) :: ca, cb, cc
      real(wp) :: ta, tb, tc
      integer :: ierr, step_c
      logical :: ok
      call make_engine(ea, ca, case_nml(ratio, .false.), ok)
      call check(error, ok, "straight engine setup failed")
      if (allocated(error)) return
      ta = 0.0_wp
      call advance(ea, ca, ta, 8 + 4, ok)
      call check(error, ok, "straight run failed")
      if (allocated(error)) return

      call make_engine(eb, cb, case_nml(ratio, .false.), ok)
      call check(error, ok, "writer engine setup failed")
      if (allocated(error)) return
      tb = 0.0_wp
      call advance(eb, cb, tb, 8, ok)
      call check(error, ok, "writer run failed")
      if (allocated(error)) return
      ! the driver's FIX-1 flush: an exact no-op on a window boundary
      call ocean_dyn_flush_tracer_window(eb%grid, eb%state%metrics, eb%state%dyn, &
                                         eb%state%continuity, eb%state%multilayer, &
                                         bc=eb%state%bc)
      call ocean_state_restart_write(eb%state, eb%grid, eb%decomp, FN, tb, 8, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "checkpoint write failed")
      if (allocated(error)) return
      call engine_exit_data(eb)
      call engine_teardown(eb)

      call make_engine(ec, cc, case_nml(ratio, .false.), ok, restart_file=FN, t0=tc, step0=step_c)
      call check(error, ok, "warm restart setup failed")
      if (allocated(error)) return
      call check(error, step_c == 8, "restart step count not restored")
      if (allocated(error)) return
      call advance(ec, cc, tc, 4, ok)
      call check(error, ok, "resumed run failed")
    if (.not. allocated(error)) call compare_registries(error, ea, ec, "after the resumed window, ratio "//char(ichar("0") + ratio))
      call engine_exit_data(ea)
      call engine_teardown(ea)
      call engine_exit_data(ec)
      call engine_teardown(ec)
      call delete_file(FN)
   end subroutine restart_at

   subroutine test_cadence(error)
      !! The predicates, step by step: windows close on steps ratio, 2*ratio,
      !! ... (never earlier), and with equal ratios every remap step is a
      !! drain step (the drain runs first in the step body).
      type(error_type), allocatable, intent(out) :: error
      type(ocean_dyn_t) :: d
      integer :: n, nlast, r
      logical :: ok_drain, ok_order
      do r = 2, 8, 2
         d%dt_tracer_advect_ratio = r
         d%dt_therm_ratio = r
         ok_drain = .true.
         ok_order = .true.
         nlast = 0
         do n = 1, 5*r
            d%outer_step_count = n - 1     ! pre-increment, as in the step body
            if (d%is_tracer_advect_step() .neqv. (mod(n, r) == 0)) ok_drain = .false.
            if (d%is_remap_step() .and. .not. d%is_tracer_advect_step()) ok_order = .false.
            if (d%is_remap_step() .neqv. (mod(n, r) == 0)) ok_order = .false.
         end do
         call check(error, ok_drain, "drain does not close on steps ratio, 2*ratio, ...")
         if (allocated(error)) return
         call check(error, ok_order, "remap step is not a drain step with equal ratios")
         if (allocated(error)) return
      end do
      ! thermo ratio 8 over tracer ratio 4: remap every 8th step, drain every 4th
      d%dt_tracer_advect_ratio = 4
      d%dt_therm_ratio = 8
      d%outer_step_count = 3
      call check(error, d%is_tracer_advect_step() .and. .not. d%is_remap_step(), "step 4 of 4/8")
      if (allocated(error)) return
      d%outer_step_count = 7
      call check(error, d%is_tracer_advect_step() .and. d%is_remap_step(), "step 8 of 4/8")
   end subroutine test_cadence

   subroutine test_window_length(error)
      !! Through the engine: after every step n the open window holds
      !! mod(n, ratio) steps of clock (it resets on steps 4, 8, 12 only --
      !! a remap or thermo step may not cut a window).
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: e
      type(config_t) :: cfg
      real(wp) :: t
      integer :: n
      logical :: ok
      call make_engine(e, cfg, case_nml(4, .false.), ok)
      call check(error, ok, "engine setup failed")
      if (allocated(error)) return
      t = 0.0_wp
      do n = 1, 12
         call advance(e, cfg, t, 1, ok)
         call check(error, ok, "engine step failed")
         if (allocated(error)) exit
         call check(error, abs(e%state%continuity%t_dyn_rel_adv - real(mod(n, 4), wp)*DT) < 1.0e-9_wp, &
                    "window clock after step does not equal mod(step, ratio)*dt")
         if (allocated(error)) exit
      end do
      call engine_exit_data(e)
      call engine_teardown(e)
   end subroutine test_window_length

   subroutine delete_file(fname)
      character(len=*), intent(in) :: fname
      integer :: u, ios
      logical :: exists
      inquire (file=fname, exist=exists)
      if (.not. exists) return
      open (newunit=u, file=fname, status="old", iostat=ios)
      if (ios == 0) close (u, status="delete")
   end subroutine delete_file

end module test_ocean_tracer_ratio_pc
