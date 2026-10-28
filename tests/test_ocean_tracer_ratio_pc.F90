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
   character(len=*), parameter :: KD_NML = "&ocean_tidal_mixing_nml enable = .true., e_uniform = 1.0e-3 /"// &
                                  new_line("a")// &
                                  "&ocean_bdrag_nml form = 'quadratic', cd = 3.0e-3, hbbl = 10.0, bg_vel = 0.1 /"// &
                                  new_line("a")// &
                                  "&ocean_bbl_mixing_nml enable = .true. /"
      !! Tidal + law-of-the-wall BBL mixing: two persistent `kd_int` fields.
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
                  new_unittest("window_length_exact", test_window_length), &
                  new_unittest("ideal_age_surface_reset_cadence", test_age_reset), &
                  new_unittest("diabatic_last_cadence_predicates", test_df_cadence), &
                  new_unittest("diabatic_last_thermo_at_interval_end", test_df_ordering), &
                  new_unittest("diabatic_first_knob_inert_at_ratio1", test_df_ratio1), &
                  new_unittest("diabatic_last_ice_restart", test_df_ice_restart), &
                  new_unittest("diabatic_last_hdiff_dose_once", test_df_hdiff), &
                  new_unittest("diabatic_last_redi_full_dt_subcycling", test_df_redi), &
                  new_unittest("diabatic_last_restart_tidal_bbl_kd", test_df_restart_kd), &
                  new_unittest("diabatic_last_coeffs_at_interval_end", test_df_coeffs) &
                  ]
   end subroutine collect_ocean_tracer_ratio_pc_tests

   subroutine ensure_comm()
      if (.not. comm_inited) then
         call comm_env_init()
         call comm_env_setup_roles(.false.)
         comm_inited = .true.
      end if
   end subroutine ensure_comm

   function case_nml(ratio, uniform, vcoord, dfirst, therm, extra) result(nml)
      !! Closed wind-driven basin with an island (divergent flow, real `h`
      !! motion), `pred_corr` (default), `dt_therm_ratio = 4`.
      integer, intent(in) :: ratio
      logical, intent(in) :: uniform
      character(len=*), intent(in), optional :: vcoord
      logical, intent(in), optional :: dfirst
         !! `&ocean_vmix_nml diabatic_first` (default: leave unset = `.true.`).
      integer, intent(in), optional :: therm
         !! `dt_therm_ratio` (default 4).
      character(len=*), intent(in), optional :: extra
         !! Extra namelist text appended verbatim (e.g. `&ocean_bbl_mixing_nml ... /`).
      character(len=:), allocatable :: nml
      character(len=:), allocatable :: vc, df
      character(len=16) :: ts
      character(len=*), parameter :: NL = new_line("a")
      character(len=16) :: rs
      character(len=:), allocatable :: tr
      write (rs, '(i0)') ratio
      vc = 'zstar'
      if (present(vcoord)) vc = vcoord
      df = ""
      if (present(dfirst)) then
         if (.not. dfirst) df = ", diabatic_first = .false."
      end if
      write (ts, '(i0)') 4
      if (present(therm)) write (ts, '(i0)') therm
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
            "&ocean_vmix_nml dt_therm_ratio = "//trim(ts)//", dt_tracer_advect_ratio = "//trim(rs)//df//" /"//NL// &
            "&ocean_diag_nml enabled = .false. /"//NL// &
            "&grid_nml nx = 24, ny = 16, nghost = 3, dx = 10000.0, dy = 10000.0 /"//NL// &
            "&physics_nml coriolis_f = 1.0e-4, wind_stress_x = 0.08 /"//NL// &
            "&vcoord_nml vcoord_type = '"//vc//"' /"//NL// &
            "&ocean_topo_nml topo_config = 'island', max_depth = 1000.0, "// &
            "slope_scale = 0.25 /"//NL// &
            "&ocean_bc_nml west = 'wall', east = 'wall', south = 'wall', "// &
            "north = 'wall' /"//NL// &
            "&output_nml output_to_file = .false. /"//NL
      if (present(extra)) nml = nml//extra//NL
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

   subroutine advance(engine, cfg, t, nsteps, ok, dtv)
      type(ocean_engine_t), intent(inout) :: engine
      type(config_t), intent(inout) :: cfg
      real(wp), intent(inout) :: t
      integer, intent(in) :: nsteps
      logical, intent(out) :: ok
      real(wp), intent(in), optional :: dtv
         !! Step size (default `DT`).
      integer :: n, ierr
      real(wp) :: dtl
      dtl = DT
      if (present(dtv)) dtl = dtv
      ok = .false.
      do n = 1, nsteps
         call engine_step(engine, dtl, t, ierr=ierr)
         if (ierr /= OCEAN_STATUS_OK) return
         call engine_step_ice(engine, cfg, dtl, t, ierr=ierr)
         if (ierr /= OCEAN_STATUS_OK) return
         call engine_step_finalize(engine, dtl, t, ierr=ierr)
         if (ierr /= OCEAN_STATUS_OK) return
         t = t + dtl
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
         call restart_at(error, r, dfirst=.false.)
         if (allocated(error)) return
      end do
   end subroutine test_restart

   subroutine restart_at(error, ratio, dfirst, extra)
      !! Straight 12 steps vs 8 steps + (flush, checkpoint) + warm restart +
      !! 4 steps.  8 = two whole windows, so the checkpoint is window-aligned.
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: ratio
      logical, intent(in), optional :: dfirst
      character(len=*), intent(in), optional :: extra
      type(ocean_engine_t), target :: ea, eb, ec
      type(config_t) :: ca, cb, cc
      real(wp) :: ta, tb, tc
      integer :: ierr, step_c
      logical :: ok
      call make_engine(ea, ca, case_nml(ratio, .false., dfirst=dfirst, extra=extra), ok)
      call check(error, ok, "straight engine setup failed")
      if (allocated(error)) return
      ta = 0.0_wp
      call advance(ea, ca, ta, 8 + 4, ok)
      call check(error, ok, "straight run failed")
      if (allocated(error)) return

      call make_engine(eb, cb, case_nml(ratio, .false., dfirst=dfirst, extra=extra), ok)
      call check(error, ok, "writer engine setup failed")
      if (allocated(error)) return
      tb = 0.0_wp
      call advance(eb, cb, tb, 8, ok)
      call check(error, ok, "writer run failed")
      if (allocated(error)) return
      if (present(extra)) then
         associate (kt_ => eb%state%vmix_tidal%kd_int, kb_ => eb%state%vdiff%bbl_mix%kd_int)
            !$acc update self(kt_, kb_)
            call check(error, maxval(kt_) > 0.0_wp .and. maxval(kb_) > 0.0_wp, &
                       "tidal/BBL kd_int is zero at the checkpoint: restart test would be vacuous")
         end associate
         if (allocated(error)) return
      end if
      ! the driver's FIX-1 flush: an exact no-op on a window boundary
      call ocean_dyn_flush_tracer_window(eb%grid, eb%state%metrics, eb%state%dyn, &
                                         eb%state%continuity, eb%state%multilayer, &
                                         bc=eb%state%bc)
      call ocean_state_restart_write(eb%state, eb%grid, eb%decomp, FN, tb, 8, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "checkpoint write failed")
      if (allocated(error)) return
      call engine_exit_data(eb)
      call engine_teardown(eb)

      call make_engine(ec, cc, case_nml(ratio, .false., dfirst=dfirst, extra=extra), ok, restart_file=FN, t0=tc, step0=step_c)
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

   subroutine test_age_reset(error)
      !! Ideal age: the surface Dirichlet reset sits on the THERMO step (with the
      !! aging source and the column mixing), for every tracer-advect ratio.
      !! Run one thermo interval (4 steps) at ratio 4 and ratio 1 on the same
      !! stratified wind-driven basin.  After step 1 the surface layer is
      !! exactly zero; after steps 2 and 3 (no column mixing, no reset) it
      !! stays at the ratio-1 value to horizontal-advection accuracy.  Before
      !! the fix ratio 4 reset at step 4, leaving the mixed-up subsurface age
      !! in the surface layer for steps 1-3.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: e4, e1
      type(config_t) :: c4, c1
      character(len=:), allocatable :: nml4, nml1
      real(wp) :: t4, t1, s4, s1, ref
      integer :: n
      logical :: ok
      nml4 = case_nml(4, .false.)//"&ocean_tracers_nml enable_ideal_age = .true. /"//new_line("a")
      nml1 = case_nml(1, .false.)//"&ocean_tracers_nml enable_ideal_age = .true. /"//new_line("a")
      call make_engine(e4, c4, nml4, ok)
      call check(error, ok, "engine (ratio 4) setup failed")
      if (allocated(error)) return
      call make_engine(e1, c1, nml1, ok)
      call check(error, ok, "engine (ratio 1) setup failed")
      if (allocated(error)) return
      t4 = 0.0_wp
      t1 = 0.0_wp
      do n = 1, 3
         call advance(e4, c4, t4, 1, ok)
         if (ok) call advance(e1, c1, t1, 1, ok)
         call check(error, ok, "engine failed during the age run")
         if (allocated(error)) exit
         call surface_age(e4, s4, ref)
         call surface_age(e1, s1, ref)
         call check(error, ref > 0.0_wp, "no interior age: test would pass vacuously")
         if (allocated(error)) exit
         write (*, '(a,i0,a,es10.2,a,es10.2,a,es10.2)') "  age step ", n, ": surf ratio4=", s4, &
            " ratio1=", s1, " interior=", ref
         if (n == 1) call check(error, s4 == 0.0_wp, "ratio 4: surface age not reset on the thermo step")
         if (.not. allocated(error)) &
            call check(error, abs(s4 - s1) <= 1.0e-3_wp*ref, &
                       "ratio 4 surface age departs from ratio 1 within the interval")
         if (allocated(error)) exit
      end do
      call engine_exit_data(e4)
      call engine_teardown(e4)
      call engine_exit_data(e1)
      call engine_teardown(e1)
   end subroutine test_age_reset

   subroutine test_df_cadence(error)
      !! `diabatic_first = .false.`: thermo and remap fire together on steps
      !! N, 2N, ... (once per interval, on a drain step); the GM chain keeps
      !! leading the interval (steps 1, N+1, ...).  `diabatic_first = .true.`
      !! keeps its legacy cadence (thermo leads; remap at the end only when
      !! windowed advection is on).
      type(error_type), allocatable, intent(out) :: error
      type(ocean_dyn_t) :: d
      integer :: n, r
      logical :: ok
      do r = 1, 4, 1
         if (mod(4, r) /= 0) cycle
         d%dt_therm_ratio = 4
         d%dt_tracer_advect_ratio = r
         d%diabatic_first = .false.
         ok = .true.
         do n = 1, 16
            d%outer_step_count = n - 1
            if (d%is_thermo_step() .neqv. (mod(n, 4) == 0)) ok = .false.
            if (d%is_remap_step() .neqv. (mod(n, 4) == 0)) ok = .false.
            if (d%is_gm_refresh_step() .neqv. (mod(n, 4) == 1)) ok = .false.
            if (r > 1 .and. d%is_remap_step() .and. .not. d%is_tracer_advect_step()) ok = .false.
         end do
         call check(error, ok, "diabatic_first=F cadence wrong")
         if (allocated(error)) return
      end do
      ! legacy cadence (default) is untouched
      d%dt_tracer_advect_ratio = 1
      d%diabatic_first = .true.
      ok = .true.
      do n = 1, 16
         d%outer_step_count = n - 1
         if (d%is_thermo_step() .neqv. (mod(n, 4) == 1)) ok = .false.
         if (d%is_remap_step() .neqv. (mod(n, 4) == 1)) ok = .false.
         if (d%is_gm_refresh_step() .neqv. (mod(n, 4) == 1)) ok = .false.
      end do
      call check(error, ok, "diabatic_first=T legacy cadence changed")
   end subroutine test_df_cadence

   subroutine test_df_ordering(error)
      !! Through the engine, ratio 4, one thermo interval = 4 steps.  The
      !! ideal-age source is a thermo-cadence kernel (`+therm_dt*h` below the
      !! surface): with `diabatic_first = .false.` it is exactly 0 after steps
      !! 1-3, appears on step 4 (the interval end, with the drain and the
      !! remap), is unchanged over steps 5-7 and doubles on step 8.  The
      !! default fires it on step 1.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: e, ed
      type(config_t) :: c, cd
      character(len=:), allocatable :: nml
      real(wp) :: t, td, s, a(8), ad1
      integer :: n
      logical :: ok
      nml = case_nml(4, .false., dfirst=.false.)//"&ocean_tracers_nml enable_ideal_age = .true. /"//new_line("a")
      call make_engine(e, c, nml, ok)
      call check(error, ok, "engine setup failed")
      if (allocated(error)) return
      t = 0.0_wp
      do n = 1, 8
         call advance(e, c, t, 1, ok)
         call check(error, ok, "engine step failed")
         if (allocated(error)) exit
         call surface_age(e, s, a(n))
      end do
      if (.not. allocated(error)) then
         write (*, '(a,8es10.2)') "  interior age after steps 1..8: ", a
         call check(error, all(a(1:3) == 0.0_wp), "thermo work ran before the interval end")
         if (.not. allocated(error)) call check(error, a(4) > 0.0_wp, "no thermo work on step 4")
         if (.not. allocated(error)) &
            call check(error, all(abs(a(5:7) - a(4)) <= 0.02_wp*a(4)), "thermo work ran mid-interval")
         if (.not. allocated(error)) &
            call check(error, abs(a(8) - 2.0_wp*a(4)) <= 0.1_wp*a(4), "thermo did not fire exactly once on step 8")
      end if
      if (.not. allocated(error)) then
         ! the deferred thermo ran exactly twice (steps 4 and 8), each time
         ! AFTER the drain had closed the window (window clock 0; it would
         ! read 3*DT if the thermo still ran before the step-4 drain)
         call check(error, e%state%dyn%n_deferred_thermo == 2, "deferred thermo call count")
         if (.not. allocated(error)) &
            call check(error, e%state%dyn%deferred_thermo_clock == 0.0_wp, &
                       "deferred thermo ran before the window drain")
      end if
      call engine_exit_data(e)
      call engine_teardown(e)
      if (allocated(error)) return
      nml = case_nml(4, .false.)//"&ocean_tracers_nml enable_ideal_age = .true. /"//new_line("a")
      call make_engine(ed, cd, nml, ok)
      call check(error, ok, "default engine setup failed")
      if (allocated(error)) return
      td = 0.0_wp
      call advance(ed, cd, td, 1, ok)
      call surface_age(ed, s, ad1)
      call check(error, ad1 > 0.0_wp, "default (diabatic_first) no longer fires thermo on step 1")
      call engine_exit_data(ed)
      call engine_teardown(ed)
   end subroutine test_df_ordering

   function polar_ice_nml(ratio, dfirst) result(nml)
      !! The polar freeze-up case (validation_examples/ocean/sea_ice/
      !! polar_freezeup_thermo.nml, shrunk): frazil nucleates in the first
      !! steps, so the ice state is live across the checkpoint.
      integer, intent(in) :: ratio
      logical, intent(in) :: dfirst
      character(len=:), allocatable :: nml
      character(len=*), parameter :: NL = new_line("a")
      character(len=8) :: rs
      write (rs, '(i0)') ratio
      nml = "&sim_nml sim_type = 'ocean' /"//NL// &
            "&grid_nml nx = 20, ny = 16, dx = 10000.0, dy = 10000.0, nghost = 3 /"//NL// &
            "&time_nml t_end = 86400.0, dt_fixed = 600.0 /"//NL// &
            "&physics_nml coriolis_f = 1.435e-4 /"//NL// &
            "&nonhydrostatic_nml nz_layers = 4 /"//NL// &
            "&tracer_nml initial_salinity = 34.0, initial_temperature = -1.80 /"//NL// &
            "&ocean_ic_nml alpha_T = 1.7e-4, rho_0 = 1027.0 /"//NL// &
            "&ocean_thermo_nml enable_thermodynamics = .true., q_heat = -60.0 /"//NL// &
            "&ocean_topo_nml topo_config = 'flat', max_depth = 50.0 /"//NL// &
            "&ocean_hvisc_nml nu_h = 500.0 /"//NL// &
            "&ocean_vmix_nml kv_ml_invz2 = 0.01, hmix_fixed = 20.0, dt_therm_ratio = 4, "// &
      "dt_tracer_advect_ratio = "//trim(rs)//merge("                          ", ", diabatic_first = .false.", dfirst)//" /"//NL// &
            "&vcoord_nml vcoord_type = 'sigma' /"//NL// &
            "&ocean_bt_nml auto_n_inner = .true. /"//NL// &
            "&ocean_ice_nml enable = .true., ncat = 5, nk_ice = 2, air_temp = -20.0, "// &
            "restore_lambda = 20.0, sw_down = 0.0, transport = .false., dynamics = .false. /"//NL// &
            "&ocean_diag_nml enabled = .false. /"//NL// &
            "&output_nml output_to_file = .false. /"//NL
   end function polar_ice_nml

   subroutine test_df_ice_restart(error)
      !! Ice ON, diabatic_first = .false., windowed advection: 8 steps +
      !! checkpoint at the interval boundary + 4 resumed steps equals 12
      !! straight steps bit for bit (the ice/flux block prepares the fluxes
      !! of the interval-end thermo step after step N-1, so a checkpoint at
      !! step 8 must not strand them).  Frazil must have formed (live ice).
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: ea, eb, ec
      type(config_t) :: ca, cb, cc
      real(wp) :: ta, tb, tc
      integer :: ierr, step_c
      logical :: ok
      call make_engine(ea, ca, polar_ice_nml(4, .false.), ok)
      call check(error, ok, "straight engine setup failed")
      if (allocated(error)) return
      ta = 0.0_wp
      call advance(ea, ca, ta, 12, ok)
      call check(error, ok, "straight run failed")
      if (allocated(error)) return
      call make_engine(eb, cb, polar_ice_nml(4, .false.), ok)
      call check(error, ok, "writer engine setup failed")
      if (allocated(error)) return
      tb = 0.0_wp
      call advance(eb, cb, tb, 8, ok)
      call check(error, ok, "writer run failed")
      if (allocated(error)) return
      call ocean_dyn_flush_tracer_window(eb%grid, eb%state%metrics, eb%state%dyn, &
                                         eb%state%continuity, eb%state%multilayer, &
                                         bc=eb%state%bc)
      call ocean_state_restart_write(eb%state, eb%grid, eb%decomp, FN, tb, 8, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "checkpoint write failed")
      if (allocated(error)) return
      call engine_exit_data(eb)
      call engine_teardown(eb)
      call make_engine(ec, cc, polar_ice_nml(4, .false.), ok, restart_file=FN, t0=tc, step0=step_c)
      call check(error, ok, "warm restart setup failed")
      if (allocated(error)) return
      call advance(ec, cc, tc, 4, ok)
      call check(error, ok, "resumed run failed")
      if (.not. allocated(error)) call compare_registries(error, ea, ec, "ice on, diabatic_first=F")
      call engine_exit_data(ea)
      call engine_teardown(ea)
      call engine_exit_data(ec)
      call engine_teardown(ec)
      call delete_file(FN)
   end subroutine test_df_ice_restart

   subroutine test_df_ratio1(error)
      !! At `dt_therm_ratio = 1` (every step a thermo step) the knob changes
      !! nothing: 6 steps with each setting are bit-identical, with and
      !! without windowed-free ratio-1 tracers.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: ea, eb
      type(config_t) :: ca, cb
      real(wp) :: ta, tb
      logical :: ok
      call make_engine(ea, ca, case_nml(1, .false., therm=1), ok)
      call check(error, ok, "engine (default) setup failed")
      if (allocated(error)) return
      call make_engine(eb, cb, case_nml(1, .false., dfirst=.false., therm=1), ok)
      call check(error, ok, "engine (diabatic_first=F) setup failed")
      if (allocated(error)) return
      ta = 0.0_wp
      tb = 0.0_wp
      call advance(ea, ca, ta, 6, ok)
      if (ok) call advance(eb, cb, tb, 6, ok)
      call check(error, ok, "run failed")
      if (.not. allocated(error)) call compare_registries(error, ea, eb, "across diabatic_first at ratio 1")
      call engine_exit_data(ea)
      call engine_teardown(ea)
      call engine_exit_data(eb)
      call engine_teardown(eb)
   end subroutine test_df_ratio1

   subroutine surface_age(engine, smax, imax)
      !! max |age content| in the surface layer and in the layers below, interior cells.
      type(ocean_engine_t), intent(inout), target :: engine
      real(wp), intent(out) :: smax, imax
      real(wp), allocatable :: a(:, :, :)
      integer :: ng, nz, ia
      ia = engine%state%multilayer%idx_age
      associate (ha => engine%state%multilayer%tracers(ia)%hTr)
         !$acc update self(ha)
      end associate
      a = engine%state%multilayer%tracers(ia)%hTr
      ng = engine%grid%nghost
      nz = size(a, 3)
      smax = maxval(abs(a(ng + 1:size(a, 1) - ng, ng + 1:size(a, 2) - ng, nz)))
      imax = maxval(abs(a(ng + 1:size(a, 1) - ng, ng + 1:size(a, 2) - ng, 1:nz - 1)))
   end subroutine surface_age

   function flat_nml(dtv, therm, ratio, uniform, extra) result(nml)
      !! Closed flat-bottom basin, no wind, no flow: only the tracer physics
      !! acts.  `diabatic_first = .false.` whenever `therm > 1`.
      real(wp), intent(in) :: dtv
      integer, intent(in) :: therm, ratio
      logical, intent(in) :: uniform
      character(len=*), intent(in), optional :: extra
      character(len=:), allocatable :: nml
      character(len=*), parameter :: NL = new_line("a")
      character(len=16) :: ts, rs, ds
      character(len=:), allocatable :: tr, df
      write (ts, '(i0)') therm
      write (rs, '(i0)') ratio
      write (ds, '(f10.2)') dtv
      df = ""
      if (therm > 1) df = ", diabatic_first = .false."
      if (uniform) then
         tr = "&tracer_nml initial_temperature = 12.0, initial_salinity = 35.0 /"//NL
      else
         tr = "&tracer_nml initial_temperature = 12.0, initial_salinity = 35.0, "// &
              "T_init_surface = 20.0, T_init_bottom = 4.0 /"//NL
      end if
      nml = "&sim_nml sim_type = 'ocean' /"//NL// &
            "&time_nml t_end = 86400.0, dt_fixed = "//trim(adjustl(ds))//" /"//NL// &
            "&nonhydrostatic_nml nz_layers = 6 /"//NL// &
            tr// &
            "&ocean_bt_nml auto_n_inner = .true. /"//NL// &
            "&ocean_vmix_nml dt_therm_ratio = "//trim(ts)//", dt_tracer_advect_ratio = "//trim(rs)//df//" /"//NL// &
            "&ocean_diag_nml enabled = .false. /"//NL// &
            "&grid_nml nx = 24, ny = 16, nghost = 3, dx = 10000.0, dy = 10000.0 /"//NL// &
            "&physics_nml coriolis_f = 1.0e-4 /"//NL// &
            "&vcoord_nml vcoord_type = 'zstar' /"//NL// &
            "&ocean_topo_nml topo_config = 'flat', max_depth = 1000.0 /"//NL// &
            "&ocean_tracers_nml enable_ideal_age = .true. /"//NL// &
            "&ocean_bc_nml west = 'wall', east = 'wall', south = 'wall', north = 'wall' /"//NL// &
            "&output_nml output_to_file = .false. /"//NL
      if (present(extra)) nml = nml//extra//NL
   end function flat_nml

   subroutine set_age_step(engine, i_step, a_left, a_right)
      !! Overwrite the age content with a zonal step profile (concentration
      !! `a_left` for i <= i_step, `a_right` beyond) times the layer
      !! thickness, ghosts included, then push it to the device.
      type(ocean_engine_t), intent(inout), target :: engine
      integer, intent(in) :: i_step
      real(wp), intent(in) :: a_left, a_right
      integer :: ia, i
      ia = engine%state%multilayer%idx_age
      do i = 1, size(engine%state%multilayer%h_layer, 1)
         engine%state%multilayer%tracers(ia)%hTr(i, :, :) = &
            merge(a_left, a_right, i <= i_step)*engine%state%multilayer%h_layer(i, :, :)
      end do
      associate (ha => engine%state%multilayer%tracers(ia)%hTr)
         !$acc update device(ha)
      end associate
   end subroutine set_age_step

   subroutine run_age_profile(nml, nsteps, i_step, a_left, a_right, hage, ok, hlay)
      !! Engine run from the stepped age profile; returns the age content.
      character(len=*), intent(in) :: nml
      integer, intent(in) :: nsteps, i_step
      real(wp), intent(in) :: a_left, a_right
      real(wp), allocatable, intent(out) :: hage(:, :, :)
      logical, intent(out) :: ok
      real(wp), allocatable, intent(out), optional :: hlay(:, :, :)
         !! Layer thickness after the run.
      type(ocean_engine_t), target :: e
      type(config_t) :: c
      real(wp) :: t
      integer :: ia
      call make_engine(e, c, nml, ok)
      if (.not. ok) return
      call set_age_step(e, i_step, a_left, a_right)
      t = 0.0_wp
      call advance(e, c, t, nsteps, ok)
      ia = e%state%multilayer%idx_age
      associate (ha => e%state%multilayer%tracers(ia)%hTr)
         !$acc update self(ha)
      end associate
      hage = e%state%multilayer%tracers(ia)%hTr
      if (present(hlay)) then
         !$acc update self(e%state%multilayer%h_layer)
         hlay = e%state%multilayer%h_layer
      end if
      call engine_exit_data(e)
      call engine_teardown(e)
   end subroutine run_age_profile

   subroutine test_df_hdiff(error)
      !! `diabatic_first = .false.`: the tracer Laplacian is applied ONCE per
      !! interval, at the interval end, with the full tracer dt.  A zonal step
      !! of age content on a flat, motionless basin: the hdiff increment (run
      !! with `kappa_h` minus the identical run without) is zero after steps
      !! 1-3 and, after step 4, equals the analytic one-step explicit Laplacian
      !! `D h (A(i+1) - 2 A(i) + A(i-1))`, `D = kappa_h dt_therm / dx^2`, in every
      !! layer below the (reset) surface -- for tracer-advect ratio 4 and 1.
      !! Applying it twice, with dt instead of dt_therm, or not at all, fails.
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: KAP = 4166.6666666666667_wp   ! D = 0.1 at dt_therm = 2400 s
      real(wp), parameter :: AL = 2.0e6_wp, AR = 5.0e5_wp
      integer, parameter :: IC = 15
      real(wp), allocatable :: h0(:, :, :), hk(:, :, :), hz(:, :, :), hl(:, :, :)
      character(len=:), allocatable :: nk
      real(wp) :: dmax, err, d_exp, hh
      integer :: r, nz, j, k, i
      logical :: ok
      nk = "&ocean_hdiff_nml kappa_h = 4166.6666666666667 /"
      do r = 1, 4, 3
         call run_age_profile(flat_nml(DT, 4, r, .true.), 3, IC, AL, AR, h0, ok)
         call check(error, ok, "control run (3 steps) failed")
         if (allocated(error)) return
         call run_age_profile(flat_nml(DT, 4, r, .true., nk), 3, IC, AL, AR, hk, ok)
         call check(error, ok, "kappa_h run (3 steps) failed")
         if (allocated(error)) return
         call check(error, all(h0 == hk), "hdiff acted before the interval end")
         if (allocated(error)) return
         call run_age_profile(flat_nml(DT, 4, r, .true.), 4, IC, AL, AR, hz, ok)
         if (ok) call run_age_profile(flat_nml(DT, 4, r, .true., nk), 4, IC, AL, AR, hk, ok, hl)
         call check(error, ok, "4-step runs failed")
         if (allocated(error)) return
         nz = size(hk, 3)
         dmax = 0.0_wp
         err = 0.0_wp
         j = size(hk, 2)/2
         do k = 1, nz - 1
            do i = 3, size(hk, 1) - 2
               hh = 0.5_wp*(hl(IC, j, k) + hl(IC + 1, j, k))   ! face-mean thickness
               ! analytic increment of the content: D h (A(i+1) - 2A(i) + A(i-1))
               d_exp = 0.0_wp
               if (i == IC) d_exp = 0.1_wp*hh*(AR - AL)
               if (i == IC + 1) d_exp = 0.1_wp*hh*(AL - AR)
               dmax = max(dmax, abs(d_exp))
               err = max(err, abs((hk(i, j, k) - hz(i, j, k)) - d_exp))
            end do
         end do
         write (*, '(a,i0,a,es10.2,a,es10.2)') "  hdiff dose, ratio ", r, ": max|expected|=", dmax, &
            " max err=", err
         call check(error, dmax > 1.0_wp, "no hdiff increment expected: test would pass vacuously")
         if (allocated(error)) return
         call check(error, err <= 1.0e-10_wp*dmax, "hdiff increment differs from one full-dt explicit step")
         if (allocated(error)) return
      end do
   end subroutine test_df_hdiff

   subroutine redi_nitt(nml, nsteps, dtv, n_itt, ok)
      character(len=*), intent(in) :: nml
      integer, intent(in) :: nsteps
      real(wp), intent(in) :: dtv
      integer, intent(out) :: n_itt
      logical, intent(out) :: ok
      type(ocean_engine_t), target :: e
      type(config_t) :: c
      real(wp) :: t
      call make_engine(e, c, nml, ok)
      n_itt = -1
      if (.not. ok) return
      t = 0.0_wp
      call advance(e, c, t, nsteps, ok, dtv)
      n_itt = e%state%redi%n_itt
      call engine_exit_data(e)
      call engine_teardown(e)
   end subroutine redi_nitt

   subroutine test_df_redi(error)
      !! Redi under `diabatic_first = .false.` takes the FULL tracer dt after
      !! the drain, so its CFL-driven subcycle count is the one a single
      !! 2400 s step (dt = 2400, ratio 1) needs -- more than the 600 s one.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), parameter :: RD = "&ocean_redi_nml enable = .true., khtr = 1.2e5, "// &
                                     "check_diffusive_cfl = .true. /"
      integer :: na, nb, nc
      logical :: ok
      call redi_nitt(flat_nml(DT, 4, 4, .false., RD), 4, DT, na, ok)
      call check(error, ok, "df=F ratio 4 run failed")
      if (allocated(error)) return
      call redi_nitt(flat_nml(4.0_wp*DT, 1, 1, .false., RD), 1, 4.0_wp*DT, nb, ok)
      call check(error, ok, "dt=2400 reference run failed")
      if (allocated(error)) return
      call redi_nitt(flat_nml(DT, 1, 1, .false., RD), 1, DT, nc, ok)
      call check(error, ok, "dt=600 reference run failed")
      if (allocated(error)) return
      write (*, '(a,i0,a,i0,a,i0)') "  redi n_itt: df=F ratio4=", na, "  dt2400=", nb, "  dt600=", nc
      call check(error, nb > nc, "reference runs do not discriminate dt (n_itt equal): test is vacuous")
      if (allocated(error)) return
      call check(error, na == nb, "Redi n_itt is not the one the full tracer dt needs")
   end subroutine test_df_redi

   subroutine test_df_restart_kd(error)
      !! `diabatic_first = .false.` + tidal + law-of-the-wall BBL mixing: the
      !! checkpoint at the interval boundary must carry both `kd_int` fields
      !! (their refresh is the interval's last step), so 8 + 4 steps with a
      !! warm restart is bit-identical to 12 straight, ratio 4 and 1.
      type(error_type), allocatable, intent(out) :: error
      integer :: r
      do r = 1, 4, 3
         call restart_at(error, r, dfirst=.false., extra=KD_NML)
         if (allocated(error)) return
      end do
   end subroutine test_df_restart_kd

   subroutine pull_kd_kv(engine, kt, kb, kv)
      !! Host copies of the tidal and BBL `kd_int` and of the assembled `kv`.
      type(ocean_engine_t), intent(inout), target :: engine
      real(wp), allocatable, intent(out) :: kt(:, :, :), kb(:, :, :), kv(:, :, :)
      associate (a_ => engine%state%vmix_tidal%kd_int, b_ => engine%state%vdiff%bbl_mix%kd_int, &
                 c_ => engine%state%vmix%kv)
         !$acc update self(a_, b_, c_)
      end associate
      kt = engine%state%vmix_tidal%kd_int
      kb = engine%state%vdiff%bbl_mix%kd_int
      kv = engine%state%vmix%kv
   end subroutine pull_kd_kv

   subroutine test_df_coeffs(error)
      !! `diabatic_first = .false.`: the mixing coefficients are assembled ONCE
      !! per interval, at its end, after the drain, and carried by the next
      !! interval's dynamics.  Tidal + BBL mixing on a stratified wind-driven
      !! basin, steps 1..8: the slow-solve `kd_int` fields are exactly zero
      !! through steps 1-3 (no in-stage refresh), appear on step 4, are
      !! bit-identical through steps 5-7 and change on step 8; the assembled
      !! `kv` of steps 5-7 is bit-identical to the one assembled at the end of
      !! step 4 (max-abs difference 0) and is rebuilt on step 8.  The deferred
      !! call ran after the drain closed the window (window clock 0).
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: e, ed
      type(config_t) :: c, cd
      real(wp), allocatable :: kt(:, :, :), kb(:, :, :), kv(:, :, :)
      real(wp), allocatable :: kt4(:, :, :), kb4(:, :, :), kv4(:, :, :)
      real(wp) :: t, dkv(8), dkt(8), dkb(8), kvmax4
      integer :: n
      logical :: ok
      call make_engine(e, c, case_nml(4, .false., dfirst=.false., extra=KD_NML), ok)
      call check(error, ok, "engine setup failed")
      if (allocated(error)) return
      t = 0.0_wp
      dkv = 0.0_wp
      dkt = 0.0_wp
      dkb = 0.0_wp
      kvmax4 = 0.0_wp
      do n = 1, 8
         call advance(e, c, t, 1, ok)
         call check(error, ok, "engine step failed")
         if (allocated(error)) exit
         call pull_kd_kv(e, kt, kb, kv)
         if (n == 4) then
            kt4 = kt
            kb4 = kb
            kv4 = kv
            kvmax4 = maxval(abs(kv))
         end if
         if (n >= 5) then
            dkt(n) = maxval(abs(kt - kt4))
            dkb(n) = maxval(abs(kb - kb4))
            dkv(n) = maxval(abs(kv - kv4))
         else
            dkt(n) = maxval(abs(kt))
            dkb(n) = maxval(abs(kb))
         end if
      end do
      if (.not. allocated(error)) then
         write (*, '(a,8es9.1)') "  tidal kd: |kd| at 1-4, then max|kd-kd(4)|: ", dkt
         write (*, '(a,8es9.1)') "  bbl   kd: |kd| at 1-4, then max|kd-kd(4)|: ", dkb
         write (*, '(a,8es9.1)') "  max|kv - kv(4)|                          : ", dkv
         call check(error, all(dkt(1:3) == 0.0_wp) .and. all(dkb(1:3) == 0.0_wp), &
                    "slow solves refreshed before the interval end")
         if (.not. allocated(error)) &
            call check(error, dkt(4) > 0.0_wp .and. dkb(4) > 0.0_wp .and. kvmax4 > 0.0_wp, &
                       "no tidal/BBL diffusivity after the interval-end build")
         if (.not. allocated(error)) &
            call check(error, all(dkt(5:7) == 0.0_wp) .and. all(dkb(5:7) == 0.0_wp), &
                       "kd_int moved inside the interval")
         if (.not. allocated(error)) &
            call check(error, all(dkv(5:7) == 0.0_wp), &
                       "kv used by the next interval differs from the one assembled at step N")
         if (.not. allocated(error)) &
            call check(error, dkt(8) > 0.0_wp .and. dkv(8) > 0.0_wp, &
                       "coefficients not rebuilt at the next interval end")
         if (.not. allocated(error)) &
            call check(error, e%state%dyn%deferred_thermo_clock == 0.0_wp, &
                       "coefficient build ran before the window drain")
      end if
      call engine_exit_data(e)
      call engine_teardown(e)
      if (allocated(error)) return
      ! control: the default ordering rebuilds kv inside the interval
      call make_engine(ed, cd, case_nml(4, .false., extra=KD_NML), ok)
      call check(error, ok, "control engine setup failed")
      if (allocated(error)) return
      t = 0.0_wp
      call advance(ed, cd, t, 1, ok)
      call pull_kd_kv(ed, kt, kb, kv)
      kv4 = kv
      call advance(ed, cd, t, 1, ok)
      call pull_kd_kv(ed, kt, kb, kv)
      call check(error, maxval(abs(kv - kv4)) > 0.0_wp, &
                 "control (diabatic_first) kv never changes: the carried-kv check is vacuous")
      call engine_exit_data(ed)
      call engine_teardown(ed)
   end subroutine test_df_coeffs

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
