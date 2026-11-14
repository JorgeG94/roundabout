!! The CICE-port sea ice (`&ocean_ice_nml scheme = "cice"`) through the PRODUCTION engine path (D11.4):
!! `engine_setup` -> `engine_enter_data` -> `engine_step` -> `engine_step_ice` -> `engine_step_finalize`.
!!
!! Cases:
!!   * `default_untouched`     the default scheme allocates, maps and registers nothing of the new one.
!!   * `refusals`              every configure-time refusal of the first runnable version names its cause,
!!                             and a supported configuration has none (the old "not runnable" text is gone).
!!   * `hooks_installed`       every hook slot of the sequencer holds a production implementation.
!!   * `column_growth_closure` one ocean column's worth of ice under a cold atmosphere: ice grows, no error bits,
!!                             the ice ledger (water, salt, heat) closes every step and the ocean's heat and
!!                             virtual-salt content change by exactly what the ice handed over.
!!   * `column_melt_closure`   warm air and shortwave: ice melts; the same closures.
!!   * `ratio3_closure`        `dt_therm_ratio = 3` (the averager and the melt cap): ocean receives the window mean.
!!   * `smoke_2d`              closed basin, wind, thermodynamics + dynamics + transport + ridging: finite, area <= 1
!!                             after every step, mass and salt conserved against the fluxes, ridging happens.
!!   * `restart_round_trip`    N + checkpoint (mid thermo window) + N == 2N, every registry field bitwise.
!!   * `consumers_follow`      the scheme-neutral aggregates the diagnostics, console totals and the bulk forcing read
!!                             are the new state's.
module test_cice_engine
   use, intrinsic :: iso_fortran_env, only: int32, int64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_config, only: config_t, read_config_from_string
   use rdb_ocean_engine, only: ocean_engine_t, engine_setup, engine_enter_data, &
                               engine_step, engine_step_ice, engine_step_finalize, &
                               engine_exit_data, engine_teardown
   use rdb_ocean_state, only: ocean_state_restart_write, ocean_state_build_restart_registry
   use rdb_ocean_restart, only: restart_registry_t
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   use rdb_comm_env, only: comm_env_init, comm_env_setup_roles
   use rdb_ice_config, only: ice_config_check, ICE_CONFIG_MSG_LEN, ICE_CONFIG_MAX_MSGS
   use rdb_ice_ledger, only: ice_ledger_mirror_t, ICE_LEDGER_WATER, ICE_LEDGER_SALT, ICE_LEDGER_HEAT
   use rdb_ice_cice_hooks, only: ice_cice_bind, ice_cice_hooks_manifest
   use rdb_ice_step_cice, only: ICE_NSTAGE
   use rdb_ocean_console_stats, only: compute_total_tracer
   implicit none
   private

   public :: collect_cice_engine_tests

   real(wp), parameter :: DT = 600.0_wp
   character(len=*), parameter :: FN = "test_cice_engine_rt.nc"
   integer, parameter :: MAXF = 260
   logical :: comm_inited = .false.

   type :: field_t
      character(len=64) :: tag = ""
      real(wp), allocatable :: a(:, :, :)
   end type field_t

   type :: snapshot_t
      integer :: n = 0
      type(field_t) :: f(MAXF)
   end type snapshot_t

   type :: closure_t
      !! Worst closure residuals over a run, each relative to the sum of the absolute terms of its line.
      real(wp) :: water = 0.0_wp, salt = 0.0_wp, heat = 0.0_wp
      real(wp) :: ocean_heat = 0.0_wp, ocean_salt = 0.0_wp
      real(wp) :: vice_first = 0.0_wp, vice_last = 0.0_wp, aice_max = 0.0_wp
      real(wp) :: ridge_max = 0.0_wp
      logical :: finite = .true.
   end type closure_t

contains

   subroutine collect_cice_engine_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("default_untouched", test_default_untouched), &
                  new_unittest("refusals", test_refusals), &
                  new_unittest("hooks_installed", test_hooks_installed), &
                  new_unittest("column_growth_closure", test_column_growth), &
                  new_unittest("column_melt_closure", test_column_melt), &
                  new_unittest("ratio3_closure", test_ratio3), &
                  new_unittest("smoke_2d", test_smoke_2d), &
                  new_unittest("restart_round_trip", test_restart), &
                  new_unittest("consumers_follow", test_consumers), &
                  new_unittest("melt_cap_pending_ratio3", test_melt_cap_pending), &
                  new_unittest("ice_free_freezeup_closure", test_freezeup), &
                  new_unittest("island_and_cold_start", test_island), &
                  new_unittest("wind_handoff_open_and_covered", test_wind_handoff), &
                  new_unittest("sis2_file_conversion", test_sis2_conversion), &
                  new_unittest("sis2_file_wrong_ncat_refused", test_sis2_wrong_ncat) &
                  ]
   end subroutine collect_cice_engine_tests

   ! ------------------------------------------------------------------
   ! namelists and engine plumbing
   ! ------------------------------------------------------------------

   function case_nml(geom, extra) result(nml)
      !! `geom`: `column` (doubly periodic 8x8, uniform, no wind), `basin` (closed 24x16, wind),
      !! `channel` (periodic x, walls y, a seamount, wind).  `extra`: more namelist text appended
      !! (later groups of the same name override nothing: use `set_nml` to change a key).
      character(len=*), intent(in) :: geom
      character(len=*), intent(in), optional :: extra
      character(len=:), allocatable :: nml
      character(len=*), parameter :: NL = new_line("a")
      character(len=:), allocatable :: bc, topo, phys, grid, ic, atm, vmix

      bc = "&ocean_bc_nml west = 'periodic', east = 'periodic', south = 'periodic', north = 'periodic' /"//NL
      topo = "&ocean_topo_nml topo_config = 'flat', max_depth = 100.0 /"//NL
      phys = "&physics_nml coriolis_f = 1.4e-4, wind_stress_x = 0.0, wind_stress_y = 0.0 /"//NL
      grid = "&grid_nml nx = 8, ny = 8, nghost = 3, dx = 20000.0, dy = 20000.0 /"//NL
      ic = "&ocean_ice_ic_nml conc_config = 'uniform', h_ice = 0.5, conc = 1.0, h_snow = 0.0, t_ice = -10.0 /"//NL
      atm = "&ocean_ice_atm_nml uni_tas = 250.0, uni_rlds = 150.0, uni_uas = 3.0, uni_vas = 0.0 /"//NL
      vmix = "&ocean_vmix_nml dt_therm_ratio = 1 /"//NL
      select case (geom)
      case ("column")
         continue
      case ("basin")
         bc = "&ocean_bc_nml west = 'wall', east = 'wall', south = 'wall', north = 'wall' /"//NL
         phys = "&physics_nml coriolis_f = 1.0e-4, wind_stress_x = 0.1, wind_stress_y = 0.03 /"//NL
         grid = "&grid_nml nx = 24, ny = 16, nghost = 3, dx = 20000.0, dy = 20000.0 /"//NL
         ic = "&ocean_ice_ic_nml conc_config = 'uniform', h_ice = 1.0, conc = 0.9, h_snow = 0.1, t_ice = -8.0 /"//NL
         atm = "&ocean_ice_atm_nml uni_tas = 255.0, uni_rlds = 180.0, uni_uas = 12.0, uni_vas = 3.0 /"//NL
         vmix = "&ocean_vmix_nml dt_therm_ratio = 2 /"//NL
      case ("channel")
         bc = "&ocean_bc_nml west = 'periodic', east = 'periodic', south = 'wall', north = 'wall' /"//NL
         phys = "&physics_nml coriolis_f = 1.4e-4, wind_stress_x = 0.1, wind_stress_y = 0.04 /"//NL
         topo = "&ocean_topo_nml topo_config = 'seamount', max_depth = 1000.0, edge_depth = 800.0, "// &
                "slope_scale = 60000.0 /"//NL
         grid = "&grid_nml nx = 24, ny = 16, nghost = 3, dx = 20000.0, dy = 20000.0 /"//NL
         ic = "&ocean_ice_ic_nml conc_config = 'uniform', h_ice = 1.0, conc = 0.7, h_snow = 0.2, t_ice = -8.0 /"//NL
         atm = "&ocean_ice_atm_nml uni_tas = 255.0, uni_rlds = 180.0, uni_uas = 8.0, uni_vas = 2.0 /"//NL
         vmix = "&ocean_vmix_nml dt_therm_ratio = 2 /"//NL
      case default
         error stop "test_cice_engine: unknown geometry"
      end select
      nml = "&sim_nml sim_type = 'ocean' /"//NL// &
            "&time_nml t_end = 86400.0, dt_fixed = 600.0 /"//NL// &
            "&nonhydrostatic_nml nz_layers = 3 /"//NL// &
            "&tracer_nml initial_temperature = -1.0, initial_salinity = 34.0, "// &
            "T_init_surface = -1.0, T_init_bottom = -1.0 /"//NL// &
            "&ocean_bt_nml auto_n_inner = .true. /"//NL// &
            "&ocean_hvisc_nml nu_h = 200.0 /"//NL// &
            "&ocean_diag_nml enabled = .false. /"//NL// &
            grid//phys// &
            "&vcoord_nml vcoord_type = 'sigma' /"//NL// &
            topo// &
            "&ocean_thermo_nml enable_thermodynamics = .true. /"//NL// &
            "&ocean_forcing_nml enable_components = .true. /"//NL// &
            vmix// &
            "&ocean_ice_nml enable = .true., scheme = 'cice', ncat = 5, nilyr = 4, hfreeze = 10.0, "// &
            "dynamics = .true., transport = .true. /"//NL// &
            "&ocean_ice_thermo_nml update_ocn_f = .true. /"//NL// &
            "&ocean_ice_dyn_nml advection = 'upwind', ndte = 40 /"//NL// &
            ic//atm//bc// &
            "&output_nml output_to_file = .false. /"//NL
      if (present(extra)) nml = nml//extra
   end function case_nml

   subroutine ensure_comm()
      if (.not. comm_inited) then
         call comm_env_init()
         call comm_env_setup_roles(.false.)
         comm_inited = .true.
      end if
   end subroutine ensure_comm

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
      if (present(restart_file)) then
         call engine_setup(engine, cfg, ierr, restart_file=restart_file, t_restart=t0, step_restart=step0)
      else
         call engine_setup(engine, cfg, ierr)
      end if
      if (ierr /= OCEAN_STATUS_OK) return
      call engine_enter_data(engine, cfg)
      ok = .true.
   end subroutine make_engine

   subroutine finish_engine(engine)
      type(ocean_engine_t), intent(inout) :: engine
      call engine_exit_data(engine)
      call engine_teardown(engine)
   end subroutine finish_engine

   subroutine one_step(engine, cfg, t, ok)
      type(ocean_engine_t), intent(inout) :: engine
      type(config_t), intent(inout) :: cfg
      real(wp), intent(inout) :: t
      logical, intent(out) :: ok
      integer :: ierr

      ok = .false.
      call engine_step(engine, DT, t, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call engine_step_ice(engine, cfg, DT, t, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call engine_step_finalize(engine, DT, t, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      t = t + DT
      ok = .true.
   end subroutine one_step

   subroutine advance(engine, cfg, t, nsteps, ok)
      type(ocean_engine_t), intent(inout) :: engine
      type(config_t), intent(inout) :: cfg
      real(wp), intent(inout) :: t
      integer, intent(in) :: nsteps
      logical, intent(out) :: ok
      integer :: n

      ok = .false.
      do n = 1, nsteps
         call one_step(engine, cfg, t, ok)
         if (.not. ok) return
      end do
   end subroutine advance

   subroutine pull_cice(engine)
      !! Device -> host of the component arrays the cases read.
      type(ocean_engine_t), intent(inout) :: engine
      !$acc update self(engine%state%ice%cice%aicen, engine%state%ice%cice%vicen, engine%state%ice%cice%vsnon)
      !$acc update self(engine%state%ice%cice%aice, engine%state%ice%cice%vice, engine%state%ice%cice%vsno)
      !$acc update self(engine%state%ice%cice%aice0, engine%state%ice%cice%aice_init, engine%state%ice%cice%ice_err)
      !$acc update self(engine%state%ice%aice, engine%state%ice%vice, engine%state%ice%vsno, engine%state%ice%aice0)
      !$acc update self(engine%state%icec%flux%dardg1dt, engine%state%icec%work%niter)
   end subroutine pull_cice

   function rho0_cp(engine) result(r)
      type(ocean_engine_t), intent(in) :: engine
      real(wp) :: r
      r = engine%state%eos%rho0*engine%state%eos%c_p
   end function rho0_cp

   subroutine ocean_totals(engine, heat, salt)
      !! `sum(rho0 cp h T A)` [J] and `sum(rho0 h S A)` [kg of salinity-mass] over the owned cells.
      type(ocean_engine_t), intent(inout) :: engine
      real(wp), intent(out) :: heat, salt
      associate (ms => engine%state%multilayer)
         heat = rho0_cp(engine)*compute_total_tracer(ms%tracers(ms%idx_temperature)%c, ms%h_layer, &
                                                     engine%state%metrics%areaT, engine%grid%nghost)
         salt = engine%state%eos%rho0*compute_total_tracer(ms%tracers(ms%idx_salinity)%c, ms%h_layer, &
                                                           engine%state%metrics%areaT, engine%grid%nghost)
      end associate
   end subroutine ocean_totals

   subroutine run_closed(engine, cfg, t, nsteps, ok, cl, flw_uni, ratio)
      !! Step with the ice ledger around every ice step; accumulate the ocean-side closure over the
      !! whole run.  The ocean receives each window's hand-over during the NEXT window, so the
      !! comparison holds every window but the last: heat `sum A (fhocn + fswthru) dt`, salt the
      !! per-CELL virtual-salt form `sum A (1000 sum(fsalt dt) - sss_handoff sum(fresh dt))` with the
      !! surface salinity of the hand-off step (`oin%sss_ocn`).  `nsteps` is a multiple of `ratio`.
      type(ocean_engine_t), intent(inout), target :: engine
      type(config_t), intent(inout) :: cfg
      real(wp), intent(inout) :: t
      integer, intent(in) :: nsteps
      logical, intent(out) :: ok
      type(closure_t), intent(out) :: cl
      real(wp), intent(in) :: flw_uni
      integer, intent(in) :: ratio
      type(ice_ledger_mirror_t) :: led
      real(wp) :: h0, s0, h1, s1, sum_h, sum_s, prev_h, prev_s, win_h, dh, ds, hext, scl, efrz, spent, stock
      integer :: n, i, j, ng, is, ie, js, je
      real(wp), allocatable :: area(:, :), acc_fs(:, :), acc_fr(:, :)

      ng = engine%grid%nghost
      is = ng + 1
      ie = ng + engine%grid%nx_phys
      js = ng + 1
      je = ng + engine%grid%ny_phys
      ok = .true.
      cl = closure_t()
      sum_h = 0.0_wp
      sum_s = 0.0_wp
      prev_h = 0.0_wp
      prev_s = 0.0_wp
      win_h = 0.0_wp
      spent = 0.0_wp
      stock = 0.0_wp
      associate (m => engine%state%metrics)
         allocate (area(size(m%areaT, 1), size(m%areaT, 2)), acc_fs(size(m%areaT, 1), size(m%areaT, 2)), &
                   acc_fr(size(m%areaT, 1), size(m%areaT, 2)))
         area = m%areaT
      end associate
      acc_fs = 0.0_wp
      acc_fr = 0.0_wp
      call ocean_totals(engine, h0, s0)
      do n = 1, nsteps
         stock = 0.0_wp
         call led%open(engine%state%ice%cice, engine%state%ice_cgrid%tarea, DT, 1)
         call one_step(engine, cfg, t, ok)
         if (.not. ok) return
         call pull_cice(engine)
         associate (c => engine%state%ice%cice, f => engine%state%icec%flux)
            !$acc update self(f%fresh, f%fsalt, f%fhocn, f%fswthru)
            !$acc update self(engine%state%ice%frazil_heat, engine%state%ice%frzmlt)
            !$acc update self(engine%state%icec%oin%sss_ocn)
            hext = 0.0_wp
            do j = js, je
               do i = is, ie
                  hext = hext + DT*area(i, j)*c%aice_init(i, j)*flw_uni
                  ! bank spent by new ice this step: offered (frzmlt > 0, the bank at the import) minus left
                  efrz = DT*area(i, j)*max(engine%state%ice%frzmlt(i, j), 0.0_wp) &
                         - area(i, j)*engine%state%ice%frazil_heat(i, j)
                  hext = hext - efrz   ! the ice state loses enthalpy by what the bank gave up
                  spent = spent + efrz
                  stock = stock + area(i, j)*engine%state%ice%frazil_heat(i, j)
                  acc_fs(i, j) = acc_fs(i, j) + f%fsalt(i, j)*DT
                  acc_fr(i, j) = acc_fr(i, j) + f%fresh(i, j)*DT
               end do
            end do
            call led%close(c, f, engine%state%ice_cgrid%tarea, 0.0_wp, hext)
            cl%water = max(cl%water, abs(led%res(ICE_LEDGER_WATER))/max(led%scale(ICE_LEDGER_WATER), tiny(1.0_wp)))
            cl%salt = max(cl%salt, abs(led%res(ICE_LEDGER_SALT))/max(led%scale(ICE_LEDGER_SALT), tiny(1.0_wp)))
            cl%heat = max(cl%heat, abs(led%res(ICE_LEDGER_HEAT))/max(led%scale(ICE_LEDGER_HEAT), tiny(1.0_wp)))
            if (.not. ieee_is_finite(led%res(1) + led%res(2) + led%res(3))) cl%finite = .false.
            if (any(c%ice_err /= 0_int32)) cl%finite = .false.
            if (n == 1) cl%vice_first = sum(c%vice(is:ie, js:je))
            cl%vice_last = sum(c%vice(is:ie, js:je))
            cl%aice_max = max(cl%aice_max, maxval(c%aice(is:ie, js:je)))
            cl%ridge_max = max(cl%ridge_max, maxval(f%dardg1dt(is:ie, js:je)))
            win_h = win_h + led%i_fhocn + led%i_fswthru
            if (mod(n, ratio) == 0) then
               prev_h = sum_h
               prev_s = sum_s
               sum_h = sum_h + win_h
               win_h = 0.0_wp
               do j = js, je
                  do i = is, ie
                     sum_s = sum_s + area(i, j)*(1000.0_wp*acc_fs(i, j) &
                                                 - engine%state%icec%oin%sss_ocn(i, j)*acc_fr(i, j))
                  end do
               end do
               acc_fs = 0.0_wp
               acc_fr = 0.0_wp
            end if
         end associate
      end do
      call ocean_totals(engine, h1, s1)
      ! freezing at the clamp RELEASES latent heat into the ocean (T back to the freezing point) and fills
      ! the bank with it: the ocean gained (bank stock + what new ice spent) on top of the hand-over
      write (*, '(a,5(1x,es12.4))') "  ocean heat: dH, bank stock, bank spent, hand-over (applied windows):", h1 - h0, stock, spent, prev_h
      dh = h1 - h0 - (stock + spent)
      ds = s1 - s0
      scl = max(abs(h1 - h0), abs(prev_h), abs(stock + spent), 1.0_wp)   ! the largest term of the line
      cl%ocean_heat = abs(dh - prev_h)/scl
      scl = max(abs(ds), abs(prev_s), 1.0_wp)
      cl%ocean_salt = abs(ds - prev_s)/scl
   end subroutine run_closed

   ! ------------------------------------------------------------------
   ! tests
   ! ------------------------------------------------------------------

   subroutine test_default_untouched(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      type(restart_registry_t) :: reg
      logical :: ok
      character(len=:), allocatable :: nml
      integer :: e
      logical :: found_new

  nml = replace(case_nml("column"), "scheme = 'cice', ncat = 5, nilyr = 4, hfreeze = 10.0, dynamics = .true., transport = .true.", &
                    "ncat = 5, nk_ice = 2")
      nml = replace(nml, "&ocean_ice_dyn_nml advection = 'upwind', ndte = 40 /", "&ocean_ice_dyn_nml ndte = 40 /")
      call make_engine(engine, cfg, nml, ok)
      call check(error, ok, "default-scheme engine setup failed")
      if (allocated(error)) return
      call check(error,.not. engine%state%ice%cice%is_init, "the new state slot is not initialised")
      call check(error,.not. engine%state%icec%is_init, "the run slots are not initialised")
      call check(error, engine%state%icec%bytes() == 0_int64, "the run slots count no bytes")
      call check(error,.not. engine%state%surface_flux%use_ice_components, "no ice flux components")
      call check(error,.not. engine%state%ice%dync%is_init, "the C-grid slot is absent under sis2_c")
      call ocean_state_build_restart_registry(engine%state, engine%grid, reg)
      found_new = .false.
      do e = 1, reg%n
         if (index(reg%entries(e)%tag, "ice_avg_") == 1 .or. index(reg%entries(e)%tag, "ice_strocnT") == 1 .or. &
             index(reg%entries(e)%tag, "ice_state_version") == 1 .or. index(reg%entries(e)%tag, "ice_aicen") == 1) &
            found_new = .true.
      end do
      call check(error,.not. found_new, "no restart tag of the new scheme under sis2")
      call finish_engine(engine)
   end subroutine test_default_untouched

   subroutine test_refusals(error)
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      character(len=ICE_CONFIG_MSG_LEN) :: msgs(ICE_CONFIG_MAX_MSGS)
      integer :: n, ierr, k
      character(len=*), parameter :: NL = new_line("a")
      character(len=64), parameter :: what(15) = [character(len=64) :: &
                                                  "advection", "ndtd", "nghost", "tripolar", "wetdry", "cavity", &
                                         "open edge", "thermodynamics", "components", "hfreeze", "px*py", "dynamics", "transport", &
                                                  "atm forcing", "dataovr wind"]
      character(len=:), allocatable :: nml
      character(len=80) :: key(15)
      logical :: hit

      call ensure_comm()
      nml = case_nml("column")
      call read_config_from_string(nml, cfg, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "the supported configuration parses")
      if (allocated(error)) return
      call ice_config_check(cfg, msgs, n)
      call check(error, n == 0, "the supported configuration has no refusal")
      if (allocated(error)) return
      do k = 1, n
         if (index(msgs(k), "not runnable") > 0) call check(error, .false., "the 'not runnable' refusal is gone")
      end do
      key = [character(len=80) :: "advection='upwind'", "ndtd /= 1", "nghost >= 3", "tripolar north fold", &
             "wetdry", "cavity_dyn", "'wall' or 'periodic'", "enable_thermodynamics", "enable_components", &
             "hfreeze > 0", "px*py > 1", "needs dynamics=.true.", "needs transport=.true.", &
             "atm_forcing_nml", "tau_x/tau_y files"]
      do k = 1, 15
         select case (k)
         case (1)
            nml = replace(case_nml("column"), "advection = 'upwind'", "advection = 'remap'")
         case (2)
            nml = replace(case_nml("column"), "ndte = 40", "ndte = 40, ndtd = 2")
         case (3)
            nml = replace(case_nml("column"), "nghost = 3", "nghost = 2")
         case (4)
            nml = replace(case_nml("column"), "north = 'periodic'", "north = 'tripolar_fold'")
         case (5)
            nml = case_nml("column", "&ocean_wetdry_nml enable = .true. /"//NL)
         case (6)
            nml = case_nml("column", "&ocean_cavity_dyn_nml enable = .true. /"//NL)
         case (7)
            nml = replace(case_nml("column"), "west = 'periodic', east = 'periodic'", "west = 'open', east = 'periodic'")
         case (8)
            nml = replace(case_nml("column"), "&ocean_thermo_nml enable_thermodynamics = .true. /", &
                          "&ocean_thermo_nml enable_thermodynamics = .false. /")
         case (9)
            nml = replace(case_nml("column"), "enable_components = .true.", "enable_components = .false.")
         case (10)
            nml = replace(case_nml("column"), "hfreeze = 10.0", "hfreeze = -1.0")
         case (11)
            nml = case_nml("column", "&mpi_nml px = 2, py = 1 /"//NL)
         case (12)
            nml = replace(case_nml("column"), "hfreeze = 10.0, dynamics = .true.", "hfreeze = 10.0, dynamics = .false.")
         case (13)
            nml = replace(case_nml("column"), "transport = .true.", "transport = .false.")
         case (14)
            nml = case_nml("column", "&ocean_atm_forcing_nml enable = .true. /"//NL)
         case (15)
            nml = case_nml("column", "&ocean_dataovr_nml tau_x_file = 'wind.nc' /"//NL)
         end select
         call read_config_from_string(nml, cfg, ierr=ierr)
         if (ierr /= OCEAN_STATUS_OK) then
            call check(error, .false., "refusal case "//trim(what(k))//": the namelist did not parse")
            return
         end if
         call ice_config_check(cfg, msgs, n)
         hit = mentions(msgs, n, trim(key(k)))
        call check(error, n > 0 .and. hit, "refusal '"//trim(what(k))//"' is reported and names its cause: "//trim(msgs(max(n, 1))))
         if (allocated(error)) return
      end do
   end subroutine test_refusals

   function mentions(msgs, n, key) result(yes)
      character(len=*), intent(in) :: msgs(:), key
      integer, intent(in) :: n
      logical :: yes
      integer :: i
      yes = .false.
      do i = 1, n
         if (index(msgs(i), key) > 0) yes = .true.
      end do
   end function mentions

   subroutine test_hooks_installed(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      logical :: ok, inst(ICE_NSTAGE)

      call make_engine(engine, cfg, case_nml("column"), ok)
      call check(error, ok, "engine setup failed")
      if (allocated(error)) return
      call ice_cice_bind(engine%state, engine%grid)
      call ice_cice_hooks_manifest(inst)
      call check(error, all(inst), "every stage of the manifest has an implementation")
      call finish_engine(engine)
   end subroutine test_hooks_installed

   subroutine report_closure(label, cl)
      character(len=*), intent(in) :: label
      type(closure_t), intent(in) :: cl
      write (*, '(a,a,5(1x,es10.3),a,3(1x,es10.3))') "  closure ", label, cl%water, cl%salt, cl%heat, cl%ocean_heat, &
         cl%ocean_salt, "  vice0/vice1/aice_max", cl%vice_first, cl%vice_last, cl%aice_max
   end subroutine report_closure

   subroutine test_column_growth(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      type(closure_t) :: cl
      real(wp) :: t
      logical :: ok
      character(len=:), allocatable :: nml

      nml = replace(case_nml("column"), "initial_temperature = -1.0", "initial_temperature = -1.84")
      nml = replace(nml, "T_init_surface = -1.0, T_init_bottom = -1.0", "T_init_surface = -1.84, T_init_bottom = -1.84")
      call make_engine(engine, cfg, nml, ok)
      call check(error, ok, "engine setup failed")
      if (allocated(error)) return
      t = 0.0_wp
      call run_closed(engine, cfg, t, 288, ok, cl, 150.0_wp, 1)
      call check(error, ok, "288 ice steps ran (no error bit, no refusal)")
      if (.not. allocated(error)) call report_closure("growth (water salt heat | ocean heat salt)", cl)
      if (.not. allocated(error)) call check(error, cl%finite, "finite, no error bit in any cell")
      if (.not. allocated(error)) call check(error, cl%vice_last > 1.02_wp*cl%vice_first, "the ice grew")
      if (.not. allocated(error)) call check(error, cl%water < 1.0e-9_wp, "ice water ledger closes")
      if (.not. allocated(error)) call check(error, cl%salt < 1.0e-9_wp, "ice salt ledger closes")
      if (.not. allocated(error)) call check(error, cl%heat < 1.0e-5_wp, "ice heat ledger closes (BL99 tolerance)")
      if (.not. allocated(error)) call check(error, cl%ocean_heat < 1.0e-6_wp, "ocean heat = what the ice handed over")
      if (.not. allocated(error)) call check(error, cl%ocean_salt < 1.0e-6_wp, "ocean salt = virtual-salt hand-off (per cell, hand-off salinity)")
      call finish_engine(engine)
   end subroutine test_column_growth

   subroutine test_column_melt(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      type(closure_t) :: cl
      real(wp) :: t
      logical :: ok
      character(len=:), allocatable :: nml

      nml = replace(case_nml("column"), "uni_tas = 250.0, uni_rlds = 150.0", "uni_tas = 275.0, uni_rlds = 300.0")
      nml = replace(nml, "hfreeze = 10.0", "hfreeze = 10.0, sw_down = 250.0")
      nml = replace(nml, "h_ice = 0.5", "h_ice = 0.8")
      call make_engine(engine, cfg, nml, ok)
      call check(error, ok, "engine setup failed")
      if (allocated(error)) return
      t = 0.0_wp
      call run_closed(engine, cfg, t, 288, ok, cl, 300.0_wp, 1)
      call check(error, ok, "288 ice steps ran")
      if (.not. allocated(error)) call report_closure("melt   (water salt heat | ocean heat salt)", cl)
      if (.not. allocated(error)) call check(error, cl%finite, "finite, no error bit in any cell")
      if (.not. allocated(error)) call check(error, cl%vice_last < 0.995_wp*cl%vice_first, "the ice melted")
      if (.not. allocated(error)) call check(error, cl%water < 1.0e-9_wp, "ice water ledger closes")
      if (.not. allocated(error)) call check(error, cl%salt < 1.0e-9_wp, "ice salt ledger closes")
      if (.not. allocated(error)) call check(error, cl%heat < 1.0e-5_wp, "ice heat ledger closes (BL99 tolerance)")
      if (.not. allocated(error)) call check(error, cl%ocean_heat < 1.0e-6_wp, "ocean heat = what the ice handed over")
      if (.not. allocated(error)) call check(error, cl%ocean_salt < 1.0e-9_wp, "ocean salt = virtual-salt hand-off")
      call finish_engine(engine)
   end subroutine test_column_melt

   subroutine test_ratio3(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      type(closure_t) :: cl
      real(wp) :: t
      logical :: ok
      character(len=:), allocatable :: nml

      nml = replace(case_nml("column"), "dt_therm_ratio = 1", "dt_therm_ratio = 3")
      nml = replace(nml, "initial_temperature = -1.0", "initial_temperature = -1.84")
      nml = replace(nml, "T_init_surface = -1.0, T_init_bottom = -1.0", "T_init_surface = -1.84, T_init_bottom = -1.84")
      call make_engine(engine, cfg, nml, ok)
      call check(error, ok, "engine setup failed")
      if (allocated(error)) return
      t = 0.0_wp
      call run_closed(engine, cfg, t, 144, ok, cl, 150.0_wp, 3)
      call check(error, ok, "144 ice steps with dt_therm_ratio = 3 ran")
      if (.not. allocated(error)) call report_closure("ratio3 (water salt heat | ocean heat salt)", cl)
      if (.not. allocated(error)) call check(error, cl%finite, "finite, no error bit in any cell")
      if (.not. allocated(error)) call check(error, cl%water < 1.0e-9_wp .and. cl%salt < 1.0e-9_wp, "ice water and salt close")
      if (.not. allocated(error)) call check(error, cl%ocean_heat < 1.0e-8_wp, "the ocean received the window-mean heat")
      if (.not. allocated(error)) call check(error, cl%ocean_salt < 1.0e-9_wp, "the ocean received the window-mean salt")
      call finish_engine(engine)
   end subroutine test_ratio3

   subroutine test_smoke_2d(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      type(closure_t) :: cl
      real(wp) :: t
      logical :: ok

      call make_engine(engine, cfg, case_nml("basin"), ok)
      call check(error, ok, "engine setup failed")
      if (allocated(error)) return
      t = 0.0_wp
      call run_closed(engine, cfg, t, 300, ok, cl, 180.0_wp, 2)
      call check(error, ok, "300 steps of the closed basin ran")
      if (.not. allocated(error)) call report_closure("basin  (water salt heat | ocean heat salt)", cl)
      if (.not. allocated(error)) write (*, '(a,es10.3)') "  basin max area loss rate to ridging [1/s]: ", cl%ridge_max
      if (.not. allocated(error)) call check(error, cl%finite, "finite, no error bit in any cell")
      if (.not. allocated(error)) call check(error, cl%aice_max <= 1.0_wp + 1.0e-9_wp, "ice area <= 1 after every step")
     if (.not. allocated(error)) call check(error, cl%water < 1.0e-9_wp, "water conserved against the fluxes (transport + ridging)")
      if (.not. allocated(error)) call check(error, cl%salt < 1.0e-9_wp, "salt conserved against the fluxes")
      if (.not. allocated(error)) call check(error, cl%ridge_max > 0.0_wp, "ridging happens under convergence")
      if (.not. allocated(error)) call check(error, cl%ocean_heat < 1.0e-8_wp, "ocean heat = what the ice handed over")
      if (.not. allocated(error)) call check(error, cl%ocean_salt < 1.0e-9_wp, "ocean salt = virtual-salt hand-off")
      call finish_engine(engine)
   end subroutine test_smoke_2d

   subroutine test_melt_cap_pending(error)
      !! Ratio 3, a warm ocean offering a small melt potential (`hfreeze` tiny): the ice draws about half
      !! of it per window.  The hand-over of window 1 reaches the ocean only at the end of window 2, so
      !! the cap of window 2 must count it as already taken: over the five steps before that debit the ice
      !! never takes (much) more than the potential the ocean offered.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      real(wp) :: t, taken, mp, area_sum
      logical :: ok
      character(len=:), allocatable :: nml
      integer :: n, i, j, ng

      nml = replace(case_nml("column"), "dt_therm_ratio = 1", "dt_therm_ratio = 3")
      nml = replace(nml, "hfreeze = 10.0", "hfreeze = 0.04")
      nml = replace(nml, "initial_temperature = -1.0", "initial_temperature = 3.0")
      nml = replace(nml, "T_init_surface = -1.0, T_init_bottom = -1.0", "T_init_surface = 3.0, T_init_bottom = 3.0")
      nml = replace(nml, "h_ice = 0.5", "h_ice = 1.0")
      call make_engine(engine, cfg, nml, ok)
      call check(error, ok, "engine setup failed")
      if (allocated(error)) return
      t = 0.0_wp
      taken = 0.0_wp
      ng = engine%grid%nghost
      area_sum = 0.0_wp
      do n = 1, 5
         call one_step(engine, cfg, t, ok)
         if (.not. ok) exit
         !$acc update self(engine%state%icec%flux%fhocn)
         do j = ng + 1, ng + engine%grid%ny_phys
            do i = ng + 1, ng + engine%grid%nx_phys
               taken = taken - engine%state%icec%flux%fhocn(i, j)*DT*engine%state%metrics%areaT(i, j)
               if (n == 1) area_sum = area_sum + engine%state%metrics%areaT(i, j)
            end do
         end do
      end do
      call check(error, ok, "5 steps ran")
      if (.not. allocated(error)) then
         mp = 1035.0_wp*engine%state%eos%c_p*0.04_wp*(3.0_wp + 0.054_wp*34.0_wp)*area_sum
         write (*, '(a,2(1x,es11.4),a,f7.3)') "  melt cap: taken / potential offered [J]:", taken, mp, "  ratio", taken/mp
         call check(error, taken > 0.3_wp*mp, "the ice does draw a good part of the offered potential")
         if (.not. allocated(error)) call check(error, taken <= 1.02_wp*mp, &
                                                "two windows never take more than the ocean offered")
      end if
      call finish_engine(engine)
   end subroutine test_melt_cap_pending

   subroutine test_freezeup(error)
      !! An initially ICE-FREE, supercooled ocean: frazil, new ice, growth under a cold atmosphere.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      type(closure_t) :: cl
      real(wp) :: t
      logical :: ok
      character(len=:), allocatable :: nml

      nml = replace(case_nml("column"), "conc_config = 'uniform', h_ice = 0.5, conc = 1.0", &
                    "conc_config = 'zero', h_ice = 0.5, conc = 0.0")
      nml = replace(nml, "initial_temperature = -1.0", "initial_temperature = -2.0")
      nml = replace(nml, "T_init_surface = -1.0, T_init_bottom = -1.0", "T_init_surface = -2.0, T_init_bottom = -2.0")
      call make_engine(engine, cfg, nml, ok)
      call check(error, ok, "engine setup failed")
      if (allocated(error)) return
      t = 0.0_wp
      call run_closed(engine, cfg, t, 288, ok, cl, 150.0_wp, 1)
      call check(error, ok, "288 steps from an ice-free ocean ran")
      if (.not. allocated(error)) call report_closure("freezeup (water salt heat | ocean heat salt)", cl)
      if (.not. allocated(error)) call check(error, cl%finite, "finite, no error bit")
      if (.not. allocated(error)) call check(error, cl%vice_first > 0.0_wp, "frazil made ice in the first step")
      if (.not. allocated(error)) call check(error, cl%vice_last > cl%vice_first, "and it grew")
      if (.not. allocated(error)) call check(error, cl%water < 1.0e-9_wp .and. cl%salt < 1.0e-9_wp, "ice water and salt close")
      if (.not. allocated(error)) call check(error, cl%heat < 1.0e-5_wp, "ice heat closes")
      if (.not. allocated(error)) call check(error, cl%ocean_heat < 1.0e-6_wp, "ocean heat closes")
      if (.not. allocated(error)) call check(error, cl%ocean_salt < 1.0e-9_wp, "ocean salt closes")
      call finish_engine(engine)
   end subroutine test_freezeup

   subroutine test_island(error)
      !! A land mask with an island inside the basin: the cold start from `&ocean_ice_ic_nml` seeds wet cells only,
      !! ice never appears on land, and the water and salt ledgers close around it.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      type(closure_t) :: cl
      real(wp) :: t
      logical :: ok
      character(len=:), allocatable :: nml
      integer :: i, j, nland, nbad, ng, nwet

      nml = replace(case_nml("basin"), "topo_config = 'flat', max_depth = 100.0", "topo_config = 'island', max_depth = 100.0, slope_scale = 0.2")
      call make_engine(engine, cfg, nml, ok)
      call check(error, ok, "engine setup failed")
      if (allocated(error)) return
      ng = engine%grid%nghost
      call pull_cice(engine)
      nland = 0
      nwet = 0
      nbad = 0
      associate (c => engine%state%ice%cice, hm => engine%state%ice_cgrid%hm)
         do j = ng + 1, ng + engine%grid%ny_phys
            do i = ng + 1, ng + engine%grid%nx_phys
               if (hm(i, j) < 0.5_wp) then
                  nland = nland + 1
                  if (c%aice(i, j) /= 0.0_wp) nbad = nbad + 1
               else
                  nwet = nwet + 1
               end if
           if (hm(i, j) >= 0.5_wp .and. (abs(c%aice(i, j) - 0.9_wp) > 1.0e-12_wp .or. abs(c%vice(i, j) - 0.9_wp) > 1.0e-12_wp)) then
                  nbad = nbad + 1
               end if
            end do
         end do
      end associate
      write (*, '(a,2(1x,i0))') "  island case: land / wet owned cells:", nland, nwet
      call check(error, nland > 0 .and. nwet > 50, "the island exists and most of the basin is wet")
      if (.not. allocated(error)) call check(error, nbad == 0, "cold start: wet cells = IC concentration and thickness, land empty")
      if (allocated(error)) then
         call finish_engine(engine)
         return
      end if
      t = 0.0_wp
      call run_closed(engine, cfg, t, 200, ok, cl, 180.0_wp, 2)
      call check(error, ok, "200 steps around the island ran")
      if (.not. allocated(error)) call report_closure("island (water salt heat | ocean heat salt)", cl)
      if (.not. allocated(error)) call check(error, cl%finite, "finite, no error bit")
      if (.not. allocated(error)) call check(error, cl%water < 1.0e-9_wp .and. cl%salt < 1.0e-9_wp, "water and salt close")
      call pull_cice(engine)
      nbad = 0
      associate (c => engine%state%ice%cice, hm => engine%state%ice_cgrid%hm)
         do j = ng + 1, ng + engine%grid%ny_phys
            do i = ng + 1, ng + engine%grid%nx_phys
               if (hm(i, j) < 0.5_wp .and. c%aice(i, j) /= 0.0_wp) nbad = nbad + 1
            end do
         end do
      end associate
      if (.not. allocated(error)) call check(error, nbad == 0, "no ice on land after the run")
      call finish_engine(engine)
   end subroutine test_island

   subroutine test_wind_handoff(error)
      !! The ocean's wind stress after the hand-off: an ICE-FREE basin passes the wind through unchanged; an ice-covered one
      !! gets `(1 - a) tau_a - strocn` face by face (`a` the face concentration of the solver): the momentum
      !! identity of the ice + ocean pair against the atmosphere.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: ea, eb
      type(config_t) :: cfg
      character(len=:), allocatable :: nml
      real(wp) :: t, dmax, dx, dy, expect
      logical :: ok
      integer :: i, j, ng

      nml = replace(case_nml("basin"), "conc_config = 'uniform', h_ice = 1.0, conc = 0.9", &
                    "conc_config = 'zero', h_ice = 1.0, conc = 0.0")
      call make_engine(ea, cfg, nml, ok)
      call check(error, ok, "ice-free engine setup failed")
      if (allocated(error)) return
      t = 0.0_wp
      call advance(ea, cfg, t, 3, ok)
      call check(error, ok, "3 ice-free steps")
      if (.not. allocated(error)) then
         !$acc update self(ea%state%surface_stress%tau_x, ea%state%surface_stress%tau_y)
         ng = ea%grid%nghost
         dx = 0.0_wp
         dy = 0.0_wp
         do j = ng + 2, ng + ea%grid%ny_phys - 1
            do i = ng + 2, ng + ea%grid%nx_phys - 1
               dx = max(dx, abs(ea%state%surface_stress%tau_x(i, j) - 0.1_wp))
               dy = max(dy, abs(ea%state%surface_stress%tau_y(i, j) - 0.03_wp))
            end do
         end do
         call check(error, dx <= 1.0e-12_wp .and. dy <= 1.0e-12_wp, "ice-free: the open ocean gets the wind unchanged")
      end if
      call finish_engine(ea)
      if (allocated(error)) return

      call make_engine(eb, cfg, case_nml("basin"), ok)
      call check(error, ok, "ice-covered engine setup failed")
      if (allocated(error)) return
      t = 0.0_wp
      call advance(eb, cfg, t, 3, ok)
      call check(error, ok, "3 ice-covered steps")
      if (.not. allocated(error)) then
         !$acc update self(eb%state%surface_stress%tau_x, eb%state%ice%tau_a_x, eb%state%ice%dync%w%strocnxE)
         !$acc update self(eb%state%ice%dync%w%aiE, eb%state%metrics%wet_u)
         ng = eb%grid%nghost
         dmax = 0.0_wp
         do j = ng + 1, ng + eb%grid%ny_phys
            do i = ng + 2, ng + eb%grid%nx_phys
               expect = ((1.0_wp - eb%state%ice%dync%w%aiE(i - 1, j))*eb%state%ice%tau_a_x(i, j) &
                         - eb%state%ice%dync%w%strocnxE(i - 1, j))*eb%state%metrics%wet_u(i, j)
               dmax = max(dmax, abs(eb%state%surface_stress%tau_x(i, j) - expect))
            end do
         end do
         call check(error, dmax <= 1.0e-12_wp, "covered: tau_ocean = (1 - a) tau_a - strocn at every face")
         call check(error, maxval(abs(eb%state%ice%dync%w%strocnxE)) > 0.0_wp, "and the ice does exert a stress")
      end if
      call finish_engine(eb)
   end subroutine test_wind_handoff

   function sis2_nml(ncat) result(nml)
      integer, intent(in) :: ncat
      character(len=:), allocatable :: nml
      character(len=8) :: sn

      write (sn, "(i0)") ncat
  nml = replace(case_nml("column"), "scheme = 'cice', ncat = 5, nilyr = 4, hfreeze = 10.0, dynamics = .true., transport = .true.", &
                    "ncat = "//trim(sn)//", nk_ice = 2")
      nml = replace(nml, "&ocean_ice_atm_nml uni_tas = 250.0, uni_rlds = 150.0, uni_uas = 3.0, uni_vas = 0.0 /", "")
      nml = replace(nml, "&ocean_ice_thermo_nml update_ocn_f = .true. /", "")
      nml = replace(nml, "&ocean_ice_dyn_nml advection = 'upwind', ndte = 40 /", "")
   end function sis2_nml

   subroutine write_sis2_file(ncat, ok, mass)
      !! A SIS2-scheme run with ice, a few steps, a checkpoint in `FN`; `mass` = its ice + snow mass [kg].
      integer, intent(in) :: ncat
      logical, intent(out) :: ok
      real(wp), intent(out) :: mass
      type(ocean_engine_t), target :: e
      type(config_t) :: cfg
      real(wp) :: t
      integer :: ierr, i, j, c, ng

      mass = 0.0_wp
      call make_engine(e, cfg, sis2_nml(ncat), ok)
      if (.not. ok) return
      t = 0.0_wp
      call advance(e, cfg, t, 4, ok)
      if (.not. ok) return
      call ocean_state_restart_write(e%state, e%grid, e%decomp, FN, t, 4, ierr=ierr)
      ok = ierr == OCEAN_STATUS_OK
      !$acc update self(e%state%ice%part_size, e%state%ice%m_ice, e%state%ice%m_snow)
      ng = e%grid%nghost
      do c = 1, ncat
         do j = ng + 1, ng + e%grid%ny_phys
            do i = ng + 1, ng + e%grid%nx_phys
               mass = mass + e%state%metrics%areaT(i, j)*e%state%ice%part_size(i, j, c)* &
                      (e%state%ice%m_ice(i, j, c) + e%state%ice%m_snow(i, j, c))
            end do
         end do
      end do
      call finish_engine(e)
   end subroutine write_sis2_file

   subroutine test_sis2_conversion(error)
      !! A checkpoint written by the SIS2 scheme starts the new scheme: converted once (mass conserved), then the run proceeds.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: e
      type(config_t) :: cfg
      real(wp) :: m_sis2, m_new, t, t0
      integer :: i, j, ng, step0
      logical :: ok

      call write_sis2_file(5, ok, m_sis2)
      call check(error, ok, "the SIS2 checkpoint was written")
      if (allocated(error)) return
      call make_engine(e, cfg, case_nml("column"), ok, restart_file=FN, t0=t0, step0=step0)
      call check(error, ok, "the new scheme starts from the SIS2 checkpoint")
      if (allocated(error)) return
      call pull_cice(e)
      ng = e%grid%nghost
      m_new = 0.0_wp
      do j = ng + 1, ng + e%grid%ny_phys
         do i = ng + 1, ng + e%grid%nx_phys
            m_new = m_new + e%state%metrics%areaT(i, j)*(e%state%ice%cice%p%c%rhoi*e%state%ice%cice%vice(i, j) &
                                                         + e%state%ice%cice%p%c%rhos*e%state%ice%cice%vsno(i, j))
         end do
      end do
      write (*, '(a,2(1x,es14.6),a,es10.2)') "  SIS2 -> CICE ice+snow mass [kg] in / out:", m_sis2, m_new, "  rel", &
         abs(m_new - m_sis2)/m_sis2
      call check(error, m_sis2 > 0.0_wp .and. abs(m_new - m_sis2) <= 1.0e-6_wp*m_sis2, "mass conserved by the conversion")
      if (.not. allocated(error)) then
         t = t0
         call advance(e, cfg, t, 10, ok)
         call check(error, ok, "the converted run proceeds")
      end if
      call finish_engine(e)
      call delete_file(FN)
   end subroutine test_sis2_conversion

   subroutine test_sis2_wrong_ncat(error)
      !! A SIS2 file written with another category count than the namelist is refused, not converted.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: e
      type(config_t) :: cfg
      real(wp) :: m, t0
      integer :: step0
      logical :: ok

      call write_sis2_file(3, ok, m)
      call check(error, ok, "the 3-category SIS2 checkpoint was written")
      if (allocated(error)) return
      call make_engine(e, cfg, case_nml("column"), ok, restart_file=FN, t0=t0, step0=step0)
      call check(error,.not. ok, "a different category count is refused")
      call delete_file(FN)
   end subroutine test_sis2_wrong_ncat

   subroutine take_snapshot(engine, snap)
      type(ocean_engine_t), intent(inout), target :: engine
      type(snapshot_t), intent(inout) :: snap
      type(restart_registry_t) :: reg
      integer :: e

      call ocean_state_build_restart_registry(engine%state, engine%grid, reg)
      snap%n = 0
      do e = 1, reg%n
         associate (en => reg%entries(e))
            if (en%device_mapped .and. en%rank == 2) then
               !$acc update self(en%p2)
            else if (en%device_mapped .and. en%rank == 3) then
               !$acc update self(en%p3)
            end if
            if (snap%n >= MAXF) error stop "test_cice_engine: raise MAXF"
            snap%n = snap%n + 1
            snap%f(snap%n)%tag = en%tag
            select case (en%rank)
            case (0)
               snap%f(snap%n)%a = reshape([en%p0], [1, 1, 1])
            case (2)
               snap%f(snap%n)%a = reshape(en%p2, [size(en%p2, 1), size(en%p2, 2), 1])
            case default
               snap%f(snap%n)%a = en%p3
            end select
         end associate
      end do
   end subroutine take_snapshot

   subroutine compare(error, sa, sb, when, n_new)
      type(error_type), allocatable, intent(out) :: error
      type(snapshot_t), intent(in) :: sa, sb
      character(len=*), intent(in) :: when
      integer, intent(out) :: n_new
      integer :: e, i, j, k, nbad, nb
      character(len=:), allocatable :: msg

      n_new = 0
      nbad = 0
      msg = ""
      call check(error, sa%n == sb%n, "registry field count differs "//when)
      if (allocated(error)) return
      do e = 1, sa%n
         if (sa%f(e)%tag /= sb%f(e)%tag) then
            call check(error, .false., "registry field order differs: "//trim(sa%f(e)%tag))
            return
         end if
         if (index(sa%f(e)%tag, "ice_avg_") == 1 .or. index(sa%f(e)%tag, "ice_strocnT") == 1 .or. &
             index(sa%f(e)%tag, "ice_aicen") == 1 .or. index(sa%f(e)%tag, "ice_qicen") == 1 .or. &
             index(sa%f(e)%tag, "ice_dync") == 1) n_new = n_new + 1
         if (any(shape(sa%f(e)%a) /= shape(sb%f(e)%a))) then
            nbad = nbad + 1
            msg = msg//" "//trim(sa%f(e)%tag)//"(shape)"
            cycle
         end if
         nb = 0
         do k = 1, size(sa%f(e)%a, 3)
            do j = 1, size(sa%f(e)%a, 2)
               do i = 1, size(sa%f(e)%a, 1)
                  if (transfer(sa%f(e)%a(i, j, k), 0_int64) /= transfer(sb%f(e)%a(i, j, k), 0_int64)) then
                     if (nb == 0) write (*, '(a,a,3(1x,i0),a,es24.16,a,es24.16)') when//": ", &
                        trim(sa%f(e)%tag)//" first at", i, j, k, ": ", sa%f(e)%a(i, j, k), " vs ", &
                        sb%f(e)%a(i, j, k)
                     nb = nb + 1
                  end if
               end do
            end do
         end do
         if (nb > 0) then
            nbad = nbad + 1
            msg = msg//" "//trim(sa%f(e)%tag)
         end if
      end do
      call check(error, nbad == 0, "fields differ "//when//":"//msg)
   end subroutine compare

   subroutine test_restart(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: ea, eb
      type(config_t) :: cfg_a, cfg_b
      type(snapshot_t) :: sa, sb, sa0, sb0
      real(wp) :: t_a, t_b, umax
      integer :: step_b, ierr, nw, na, n_new
      logical :: ok

      nw = 5          ! odd: the checkpoint falls inside a thermodynamics window (ratio 2)
      na = 4
      call make_engine(ea, cfg_a, case_nml("channel"), ok)
      call check(error, ok, "engine A setup failed")
      if (allocated(error)) return
      t_a = 0.0_wp
      call advance(ea, cfg_a, t_a, nw, ok)
      call check(error, ok, "engine A failed before the checkpoint")
      if (allocated(error)) return
      call ocean_state_restart_write(ea%state, ea%grid, ea%decomp, FN, t_a, nw, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "checkpoint write failed")
      if (allocated(error)) return
      call take_snapshot(ea, sa0)
      call advance(ea, cfg_a, t_a, na, ok)
      call check(error, ok, "engine A failed after the checkpoint")
      if (allocated(error)) return
      call take_snapshot(ea, sa)
      !$acc update self(ea%state%ice%u_ice)
      umax = maxval(abs(ea%state%ice%u_ice))
      call finish_engine(ea)
      call check(error, umax > 1.0e-7_wp, "the ice moves in the restart case")
      if (allocated(error)) return

      call make_engine(eb, cfg_b, case_nml("channel"), ok, restart_file=FN, t0=t_b, step0=step_b)
      call check(error, ok, "engine B (warm restart) setup failed")
      if (allocated(error)) return
      call check(error, step_b == nw .and. eb%warm_restart, "warm restart taken at the written step")
      if (allocated(error)) return
      call take_snapshot(eb, sb0)
      call compare(error, sa0, sb0, "at the resume point", n_new)
      if (allocated(error)) then
         call finish_engine(eb)
         return
      end if
      call check(error, n_new >= 15, "the CICE-port tags (state, averager, stress, dynamics) are in the comparison")
      if (allocated(error)) return
      call advance(eb, cfg_b, t_b, na, ok)
      call check(error, ok, "engine B failed after the restart")
      if (allocated(error)) return
      call take_snapshot(eb, sb)
      call finish_engine(eb)
      call delete_file(FN)
      call compare(error, sa, sb, "after the resumed steps", n_new)
   end subroutine test_restart

   subroutine test_consumers(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      real(wp) :: t
      logical :: ok
      integer :: i, j
      real(wp) :: dmax

      call make_engine(engine, cfg, case_nml("basin"), ok)
      call check(error, ok, "engine setup failed")
      if (allocated(error)) return
      t = 0.0_wp
      call advance(engine, cfg, t, 20, ok)
      call check(error, ok, "20 steps ran")
      if (allocated(error)) return
      call pull_cice(engine)
      dmax = 0.0_wp
      associate (c => engine%state%ice%cice, s => engine%state%ice)
         do j = 1, c%ny
            do i = 1, c%nx
               dmax = max(dmax, abs(c%aice(i, j) - s%aice(i, j)), abs(c%vice(i, j) - s%vice(i, j)), &
                          abs(c%vsno(i, j) - s%vsno(i, j)), abs(c%aice0(i, j) - s%aice0(i, j)))
            end do
         end do
         call check(error, dmax == 0.0_wp, "the scheme-neutral aggregates are the new state's, bitwise")
         call check(error, maxval(s%aice) > 0.5_wp, "and carry ice")
      end associate
      call finish_engine(engine)
   end subroutine test_consumers

   function replace(s, a, b) result(r)
      character(len=*), intent(in) :: s, a, b
      character(len=:), allocatable :: r
      integer :: p
      p = index(s, a)
      if (p == 0) then
         r = s
      else
         r = s(1:p - 1)//b//s(p + len(a):)
      end if
   end function replace

   subroutine delete_file(fn_)
      character(len=*), intent(in) :: fn_
      integer :: u, stat
      open (newunit=u, file=fn_, status="old", iostat=stat)
      if (stat == 0) close (u, status="delete")
   end subroutine delete_file

end module test_cice_engine
