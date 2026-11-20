!! CICE-named sea-ice history diagnostics of the CICE-port scheme (`rdb_ocean_diag_ice_cice`), through the
!! PRODUCTION engine path.  The bug that prompted them: under `scheme = "cice"` the CICE names were
!! registered, read the SIS2 scheme's window masses, and wrote zeros that every test accepted.
!!
!! Cases:
!!   * `registered_names`       the cice-only names are registered under `cice` and under no other scheme;
!!                              every registered ice name has a fill.
!!   * `source_mapping`         every id reads ITS OWN source: each flux/step array is stamped with a distinct
!!                              value and the fill must return exactly that value times its unit factor.
!!   * `registered_filled`      growth, melt, flooding and mixed runs: every registered name is non-zero and
!!                              finite somewhere (no allowlist).
!!   * `volume_budget_thermo`   ice at rest (no wind): `dvidtt` equals the independent change of the
!!                              ice volume over the step AND `congel+frazil+snoice-meltt-meltb-meltl` up to
!!                              the ice-sublimation term (bounded by `evap`), cell by cell, growth and melt.
!!   * `volume_budget_dynamics` dynamics on: `dvidtt + dvidtd` is the whole step's volume change; the
!!                              dynamic part is non-zero where the ice moves.
!!   * `units_and_signs`        growth column: `congel` in cm/day matches the thickness change per day;
!!                              melt column: `meltt`, `meltb` positive.
!!   * `tsfc_independent`       `Tsfc` against the area-weighted category temperature computed in the test.
!!   * `observer_bitwise`       the restart registry after N steps is bit-identical with the history
!!                              diagnostics registered and folding every step versus off.
!!   * `cost_report`            only with `RDB_DIAGS_COST=1` (else a no-op pass): ms per ice step on a larger
!!                              grid with the history set registered versus not; prints, asserts nothing.
!!   * `restart_fills_equal`    every fill continues bitwise across a checkpoint (a pure function of
!!                              restart state), documenting the partial first window after a resume.
module test_cice_diags
   use, intrinsic :: iso_fortran_env, only: int64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_config, only: config_t, read_config_from_string, validate_config
   use rdb_ocean_engine, only: ocean_engine_t, engine_setup, engine_enter_data, &
                               engine_step, engine_step_ice, engine_step_finalize, &
                               engine_exit_data, engine_teardown
   use rdb_ocean_state, only: ocean_state_restart_write, ocean_state_build_restart_registry
   use rdb_ocean_restart, only: restart_registry_t
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   use rdb_comm_env, only: comm_env_init, comm_env_setup_roles
   use rdb_ocean_diag_ice_cice, only: cice_fill, CD_CONGEL, CD_FRAZIL, CD_SNOICE, CD_MELTT, CD_MELTB, CD_MELTL, &
                                      CD_DVIDTT, CD_DVIDTD, CD_DAIDTT, CD_DAIDTD, CD_TSFC, CD_FSWABS, CD_FLWUP, &
             CD_FSENS, CD_FLAT, CD_FHOCN, CD_FSWTHRU, CD_FRESH, CD_FSALT, CD_EVAP, CD_ALBSNI, CD_EVAPI, CD_EVAPS, CD_MELTS, CD_FSWDN
   use rdb_test_clock, only: wall_seconds
   implicit none
   private

   public :: collect_cice_diags_tests

   real(wp), parameter :: DT = 600.0_wp
   real(wp), parameter :: CMD = 8.64e6_wp
   real(wp), parameter :: RHOI = 917.0_wp
   integer, parameter :: NID = 25
   integer, parameter :: MAXF = 260
   logical :: comm_inited = .false.

   character(len=12), parameter :: ID_NAME(NID) = [character(len=12) :: &
                          "congel", "frazil", "snoice", "meltt", "meltb", "meltl", "dvidtt", "dvidtd", "daidtt", "daidtd", "Tsfc", &
                      "fswabs_ai", "flwup_ai", "fsens_ai", "flat_ai", "fhocn_ai", "fswthru_ai", "fresh_ai", "fsalt_ai", "evap_ai", &
                                                   "albsni", "evapi_ai", "evaps_ai", "melts", "fswdn_ai"]
      !! The diagnostic of each `CD_*` id, in id order.

   type :: field_t
      character(len=64) :: tag = ""
      real(wp), allocatable :: a(:, :, :)
   end type field_t

   type :: snapshot_t
      integer :: n = 0
      type(field_t) :: f(MAXF)
   end type snapshot_t

contains

   subroutine collect_cice_diags_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("registered_names", test_registered_names), &
                  new_unittest("source_mapping", test_source_mapping), &
                  new_unittest("registered_filled", test_registered_filled), &
                  new_unittest("volume_budget_thermo", test_budget_thermo), &
                  new_unittest("volume_budget_dynamics", test_budget_dynamics), &
                  new_unittest("units_and_signs", test_units_signs), &
                  new_unittest("tsfc_independent", test_tsfc), &
                  new_unittest("observer_bitwise", test_observer), &
                  new_unittest("restart_fills_equal", test_restart_fills), &
                  new_unittest("engine_reuse_odd_run", test_engine_reuse), &
                  new_unittest("melts_snow_budget", test_melts), &
                  new_unittest("cost_report", test_cost_report) &
                  ]
   end subroutine collect_cice_diags_tests

   ! ------------------------------------------------------------------
   ! plumbing
   ! ------------------------------------------------------------------

   function case_nml(kind, outdir, dyn) result(nml)
      !! `kind`: `growth` (cold atmosphere, ice under it), `melt` (warm atmosphere with sunshine, warm
      !! ocean), `mixed` (a west half that grows and an east half that melts, set by `split_atm`).
      !! `outdir` non-empty turns the diagnostic stream ON, in that directory.  `dyn`: dynamics + transport.
      character(len=*), intent(in) :: kind, outdir
      logical, intent(in) :: dyn
      character(len=:), allocatable :: nml
      character(len=*), parameter :: NL = new_line("a")
      character(len=:), allocatable :: diag, atm, ic, tr, dd, phys
      real(wp) :: sw
      character(len=32) :: sws

      if (len_trim(outdir) > 0) then
         diag = "&ocean_diag_nml enabled = .true., filename = 'd.nc', dt_out = 1200.0 /"//NL// &
                "&output_nml output_to_file = .true., output_dir = '"//trim(outdir)//"' /"//NL
      else
         diag = "&ocean_diag_nml enabled = .false. /"//NL//"&output_nml output_to_file = .false. /"//NL
      end if
      ic = "&ocean_ice_ic_nml conc_config = 'uniform', h_ice = 0.5, conc = 0.9, h_snow = 0.05, t_ice = -8.0 /"//NL
      tr = "&tracer_nml initial_temperature = -1.7, initial_salinity = 34.0, "// &
           "T_init_surface = -1.7, T_init_bottom = -1.7 /"//NL
      atm = "&ocean_ice_atm_nml uni_tas = 250.0, uni_rlds = 150.0, uni_uas = 5.0, uni_vas = 0.0 /"//NL
      sw = 0.0_wp
      select case (kind)
      case ("growth")
         continue
      case ("melt")
         ic = "&ocean_ice_ic_nml conc_config = 'uniform', h_ice = 1.0, conc = 0.8, h_snow = 0.0, t_ice = -1.0 /"//NL
         tr = "&tracer_nml initial_temperature = 2.0, initial_salinity = 34.0, "// &
              "T_init_surface = 2.0, T_init_bottom = 2.0 /"//NL
         atm = "&ocean_ice_atm_nml uni_tas = 278.0, uni_rlds = 320.0, uni_uas = 5.0, uni_vas = 0.0 /"//NL
         sw = 250.0_wp
      case ("column")
         ic = "&ocean_ice_ic_nml conc_config = 'uniform', h_ice = 0.5, conc = 1.0, h_snow = 0.0, t_ice = -8.0 /"//NL
      case ("flood")
         ic = "&ocean_ice_ic_nml conc_config = 'uniform', h_ice = 0.6, conc = 0.8, h_snow = 0.4, t_ice = -4.0 /"//NL
      case ("mixed")
         ic = "&ocean_ice_ic_nml conc_config = 'uniform', h_ice = 0.6, conc = 0.8, h_snow = 0.03, t_ice = -4.0 /"//NL
         tr = "&tracer_nml initial_temperature = -1.9, initial_salinity = 34.0, "// &
              "T_init_surface = -1.9, T_init_bottom = -1.9 /"//NL
         atm = "&ocean_ice_atm_nml uni_tas = 250.0, uni_rlds = 150.0, uni_uas = 5.0, uni_vas = 0.0 /"//NL
      case default
         error stop "test_cice_diags: unknown kind"
      end select
      write (sws, '(f8.1)') sw
      dd = "dynamics = .true., transport = .true., snowfall = 0.0, sw_down = "//trim(adjustl(sws))
      ! The CICE-port dynamics, transport and ridging always run (production REFUSES `dynamics = .false.`
      ! under this scheme), so `dyn = .false.` means: no wind, no stress, the ice at rest.
      if (.not. dyn) then
         atm = replace(atm, "uni_uas = 5.0", "uni_uas = 0.0")
         phys = "&physics_nml coriolis_f = 1.0e-4, wind_stress_x = 0.0, wind_stress_y = 0.0 /"//NL
      else
         phys = "&physics_nml coriolis_f = 1.0e-4, wind_stress_x = 0.1, wind_stress_y = 0.03 /"//NL
      end if
      nml = "&sim_nml sim_type = 'ocean' /"//NL// &
            "&time_nml t_end = 86400.0, dt_fixed = 600.0 /"//NL// &
            "&nonhydrostatic_nml nz_layers = 3 /"//NL// &
            tr// &
            "&ocean_bt_nml auto_n_inner = .true. /"//NL// &
            "&ocean_hvisc_nml nu_h = 200.0 /"//NL// &
            diag// &
            "&grid_nml nx = 24, ny = 16, nghost = 3, dx = 20000.0, dy = 20000.0 /"//NL// &
            phys// &
            "&vcoord_nml vcoord_type = 'sigma' /"//NL// &
            "&ocean_topo_nml topo_config = 'flat', max_depth = 100.0 /"//NL// &
            "&ocean_thermo_nml enable_thermodynamics = .true. /"//NL// &
            "&ocean_forcing_nml enable_components = .true. /"//NL// &
            "&ocean_vmix_nml dt_therm_ratio = 2 /"//NL// &
            "&ocean_ice_nml enable = .true., scheme = 'cice', ncat = 5, nilyr = 4, hfreeze = 10.0, "//dd//" /"//NL// &
            "&ocean_ice_thermo_nml update_ocn_f = .true. /"//NL// &
            "&ocean_ice_dyn_nml advection = 'upwind', ndte = 40 /"//NL// &
            ic//atm// &
            "&ocean_bc_nml west = 'wall', east = 'wall', south = 'wall', north = 'wall' /"//NL
   end function case_nml

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

   subroutine warm_east(engine)
      !! The `mixed` kinds: the atmosphere is cold over the west half of the basin and warm and bright over
      !! the east half (per-cell uniform-atmosphere arrays), so growth and melt share one state.
      type(ocean_engine_t), intent(inout) :: engine
      integer :: i0

      i0 = engine%grid%nghost + engine%grid%nx_phys/2
      engine%state%icec%u_tas(i0 + 1:, :) = 278.0_wp
      engine%state%icec%u_rlds(i0 + 1:, :) = 320.0_wp
      engine%state%icec%u_rsds(i0 + 1:, :) = 250.0_wp
      !$acc update device(engine%state%icec%u_tas, engine%state%icec%u_rlds, engine%state%icec%u_rsds)
   end subroutine warm_east

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
      ! a configuration a user could run: production validates it before setup
      call validate_config(cfg, ierr=ierr)
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

   subroutine get_fill(engine, id, a)
      !! Run diagnostic `id` through the fill the DIAGNOSTIC MANAGER registered under its name (the procedure
      !! pointer of `state%diag%vars`, exactly what a run calls every step) into a device-mapped buffer and bring
      !! it back: `a(nx, ny)`.  Calling the kernels by id would skip the `scheme` dispatch inside the registered
      !! fills, which is where the original all-zero bug lived.
      type(ocean_engine_t), intent(inout) :: engine
      integer, intent(in) :: id
      real(wp), allocatable, intent(out) :: a(:, :)
      real(wp), allocatable :: buf(:, :, :)
      integer :: k, kk

      kk = 0
      do k = 1, engine%state%diag%nvars
         if (trim(engine%state%diag%vars(k)%name) == trim(ID_NAME(id))) kk = k
      end do
      if (kk == 0) error stop "test_cice_diags: "//trim(ID_NAME(id))//" is not registered (the engine needs the diagnostic stream)"
      allocate (buf(size(engine%state%metrics%areaT, 1), size(engine%state%metrics%areaT, 2), 1), source=-777.0_wp)
      !$acc enter data copyin(buf)
      call engine%state%diag%vars(kk)%fill(engine%state, buf)
      !$acc update self(buf)
      !$acc exit data delete(buf)
      a = buf(:, :, 1)
   end subroutine get_fill

   subroutine pull(engine)
      !! Device -> host of every array the cases read.
      type(ocean_engine_t), intent(inout) :: engine
      !$acc update self(engine%state%ice%cice%aicen, engine%state%ice%cice%vicen, engine%state%ice%cice%Tsfcn)
      !$acc update self(engine%state%ice%cice%aice, engine%state%ice%cice%vice, engine%state%ice%cice%aice_init)
      !$acc update self(engine%state%ice%cice%vicen_init, engine%state%icec%work%Tf)
      !$acc update self(engine%state%icec%flux%congel, engine%state%icec%flux%frazil_vol, engine%state%icec%flux%snoice)
      !$acc update self(engine%state%icec%flux%meltt, engine%state%icec%flux%meltb, engine%state%icec%flux%meltl)
      !$acc update self(engine%state%icec%flux%evap, engine%state%icec%flux%aice_th, engine%state%icec%flux%vice_th)
   end subroutine pull

   function owned(engine, i, j) result(yes)
      type(ocean_engine_t), intent(in) :: engine
      integer, intent(in) :: i, j
      logical :: yes
      integer :: ng
      ng = engine%grid%nghost
      yes = i > ng .and. i <= ng + engine%grid%nx_phys .and. j > ng .and. j <= ng + engine%grid%ny_phys
   end function owned

   ! ------------------------------------------------------------------
   ! cases
   ! ------------------------------------------------------------------

   subroutine test_registered_names(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      logical :: ok
      integer :: i, k, nfound
      character(len=:), allocatable :: nml

      nml = case_nml("growth", "tmp_cice_diags_reg", .true.)
      call make_engine(engine, cfg, nml, ok)
      call check(error, ok, "engine with the diagnostic stream failed to set up")
      if (allocated(error)) return
      nfound = 0
      do i = 1, NID
         do k = 1, engine%state%diag%nvars
            if (trim(engine%state%diag%vars(k)%name) == trim(ID_NAME(i))) then
               nfound = nfound + 1
               call check(error, associated(engine%state%diag%vars(k)%fill), trim(ID_NAME(i))//" has a fill")
               if (allocated(error)) exit
            end if
         end do
         if (allocated(error)) exit
      end do
      if (.not. allocated(error)) call check(error, nfound == NID, "every one of the 25 CICE names is registered under cice")
      call finish_engine(engine)
      call rm_dir("tmp_cice_diags_o*")
   end subroutine test_registered_names

   subroutine test_source_mapping(error)
      !! Stamp every source with its own value and read each diagnostic: a fill that reads the wrong array
      !! (the SIS2 window masses, a neighbouring flux) cannot return the stamp.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      logical :: ok
      real(wp), allocatable :: a(:, :)
      real(wp) :: want(NID), r
      integer :: id, ng, nx, ny, i0, j0
      character(len=:), allocatable :: nml

      nml = case_nml("growth", "tmp_cice_diags_o1", .true.)
      call make_engine(engine, cfg, nml, ok)
      call check(error, ok, "engine setup failed")
      if (allocated(error)) return
      ng = engine%grid%nghost
      nx = size(engine%state%metrics%areaT, 1)
      ny = size(engine%state%metrics%areaT, 2)
      i0 = ng + 5
      j0 = ng + 5
      associate (f => engine%state%icec%flux)
         f%dt_step = DT
         f%congel = 1.0e-3_wp
         f%frazil_vol = 2.0e-3_wp
         f%snoice = 3.0e-3_wp
         f%meltt = 4.0e-3_wp
         f%meltb = 5.0e-3_wp
         f%meltl = 6.0e-3_wp
         f%fswabs = 70.0_wp
         f%flwout = -80.0_wp
         f%fsens = 90.0_wp
         f%flat = 100.0_wp
         f%fhocn = 110.0_wp
         f%fswthru = 120.0_wp
         f%fresh = 1.3e-5_wp
         f%fsalt = 1.4e-6_wp
         f%evap = -1.5e-6_wp
         f%aice_th = 0.7_wp
         f%vice_th = 0.9_wp
         f%thx(:, :, 1) = -1.2e-6_wp
         f%thx(:, :, 2) = 2.0e-7_wp
         f%thx(:, :, 3) = 5.0e-3_wp
         !$acc update device(f%congel, f%frazil_vol, f%snoice, f%meltt, f%meltb, f%meltl, f%fswabs, f%flwout)
         !$acc update device(f%fsens, f%flat, f%fhocn, f%fswthru, f%fresh, f%fsalt, f%evap, f%aice_th, f%vice_th, f%thx)
      end associate
      associate (c => engine%state%ice%cice)
         c%aice = 0.5_wp
         c%vice = 0.8_wp
         c%aice_init = 0.4_wp
         c%vicen_init = 0.1_wp
         !$acc update device(c%aice, c%vice, c%aice_init, c%vicen_init)
      end associate
      engine%state%icec%work%sw_down = 300.0_wp
      !$acc update device(engine%state%icec%work%sw_down)
      engine%state%icec%rad%alv_ai = 0.6_wp
      engine%state%icec%rad%ali_ai = 0.3_wp
      !$acc update device(engine%state%icec%rad%alv_ai, engine%state%icec%rad%ali_ai)
      r = CMD/DT
      want = 0.0_wp
      want(CD_CONGEL) = 1.0e-3_wp*r
      want(CD_FRAZIL) = 2.0e-3_wp*r
      want(CD_SNOICE) = 3.0e-3_wp*r
      want(CD_MELTT) = 4.0e-3_wp*r
      want(CD_MELTB) = 5.0e-3_wp*r
      want(CD_MELTL) = 6.0e-3_wp*r
      want(CD_DVIDTT) = (0.9_wp - 5*0.1_wp)*r
      want(CD_DVIDTD) = (0.8_wp - 0.9_wp)*r
      want(CD_DAIDTT) = (0.7_wp - 0.4_wp)*r
      want(CD_DAIDTD) = (0.5_wp - 0.7_wp)*r
      want(CD_FSWABS) = 70.0_wp
      want(CD_FLWUP) = -80.0_wp
      want(CD_FSENS) = 90.0_wp
      want(CD_FLAT) = 100.0_wp
      want(CD_FHOCN) = 110.0_wp
      want(CD_FSWTHRU) = 120.0_wp
      want(CD_FRESH) = 1.3e-5_wp*CMD/1000.0_wp
      want(CD_FSALT) = 1.4e-6_wp
      want(CD_EVAP) = -1.5e-6_wp*CMD/1000.0_wp
      want(CD_EVAPI) = -1.2e-6_wp*CMD/1000.0_wp
      want(CD_EVAPS) = 2.0e-7_wp*CMD/1000.0_wp
      want(CD_MELTS) = 5.0e-3_wp*r
      want(CD_FSWDN) = 300.0_wp*0.4_wp
      want(CD_ALBSNI) = 100.0_wp*(0.52_wp*0.6_wp + 0.48_wp*0.3_wp)
      do id = 1, NID
         if (id == CD_TSFC) cycle      ! `tsfc_independent`
         call get_fill(engine, id, a)
         call check(error, abs(a(i0, j0) - want(id)) <= 1.0e-12_wp*max(1.0_wp, abs(want(id))), &
                    trim(ID_NAME(id))//" reads its own source")
         if (allocated(error)) then
            write (*, '(a,a,2es16.8)') "   ", trim(ID_NAME(id)), a(i0, j0), want(id)
            exit
         end if
         call check(error, ieee_is_nan(a(1, 1)), trim(ID_NAME(id))//" is NaN on the ghost ring")
         if (allocated(error)) exit
         call check(error, size(a, 1) == nx .and. size(a, 2) == ny, "shape")
         if (allocated(error)) exit
      end do
      call finish_engine(engine)
      call rm_dir("tmp_cice_diags_o*")
   end subroutine test_source_mapping

   subroutine test_registered_filled(error)
      !! Growth, melt, flooding and mixed runs: every one of the 21 diagnostics is finite on wet cells and non-zero
      !! somewhere in one of them.  No allowlist: the mixed run is supercooled, so even `frazil` fires.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      logical :: ok, seen(NID), fin
      real(wp), allocatable :: a(:, :)
      real(wp) :: t
      integer :: id, kind, i, j, n
      character(len=7), parameter :: KINDS(4) = [character(len=7) :: "growth", "melt", "flood", "mixed"]
      character(len=:), allocatable :: nml

      seen = .false.
      fin = .true.
      do kind = 1, 4
         nml = case_nml(trim(KINDS(kind)), "tmp_cice_diags_o1", .true.)
         call make_engine(engine, cfg, nml, ok)
         call check(error, ok, "engine setup failed")
         if (allocated(error)) return
         if (trim(KINDS(kind)) == "mixed") call warm_east(engine)
         t = 0.0_wp
         call advance(engine, cfg, t, 30, ok)
         call check(error, ok, "30 steps ran")
         if (allocated(error)) then
            call finish_engine(engine)
            return
         end if
         call pull(engine)
         do id = 1, NID
            call get_fill(engine, id, a)
            do j = 1, size(a, 2)
               do i = 1, size(a, 1)
                  if (.not. owned(engine, i, j)) cycle
                  if (.not. ieee_is_finite(a(i, j))) fin = .false.
                  if (abs(a(i, j)) > 0.0_wp .and. ieee_is_finite(a(i, j))) seen(id) = .true.
               end do
            end do
         end do
         call finish_engine(engine)
      end do
      call check(error, fin, "every owned cell of every diagnostic is finite")
      if (allocated(error)) return
      do id = 1, NID
         call check(error, seen(id), trim(ID_NAME(id))//" is registered but never filled (zero everywhere)")
         if (allocated(error)) return
      end do
      n = count(seen)
      write (*, '(a,i0,a,i0,a)') "  registered_filled: ", n, " of ", NID, " diagnostics non-zero"
      call rm_dir("tmp_cice_diags_o*")
   end subroutine test_registered_filled

   subroutine budget_case(error, kind, label)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: kind, label
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      logical :: ok
      real(wp), allocatable :: ei(:, :), dv(:, :), cg(:, :), fz(:, :), sn(:, :), mt(:, :), mb(:, :), ml(:, :), dd(:, :)
      real(wp), allocatable :: v0(:, :), da(:, :), a0(:, :), dad(:, :)
      real(wp) :: t, ident, resid, worst, vmax, sgrow, smelt
      integer :: n, i, j, ng
      logical :: bad_state, bad_ident
      real(wp) :: e1, e2, e3
      character(len=:), allocatable :: nml

      nml = case_nml(kind, "tmp_cice_diags_o1", .false.)
      call make_engine(engine, cfg, nml, ok)
      call check(error, ok, label//": engine setup failed")
      if (allocated(error)) return
      if (kind == "mixed") call warm_east(engine)
      ng = engine%grid%nghost
      t = 0.0_wp
      bad_state = .false.
      bad_ident = .false.
      worst = 0.0_wp
      e1 = 0.0_wp
      e2 = 0.0_wp
      e3 = 0.0_wp
      sgrow = 0.0_wp
      smelt = 0.0_wp
      do n = 1, 24
         call pull(engine)
         v0 = engine%state%ice%cice%vice
         a0 = engine%state%ice%cice%aice
         call one_step(engine, cfg, t, ok)
         call check(error, ok, label//": step ran")
         if (allocated(error)) exit
         call pull(engine)
         call get_fill(engine, CD_DVIDTT, dv)
         call get_fill(engine, CD_DVIDTD, dd)
         call get_fill(engine, CD_DAIDTT, da)
         call get_fill(engine, CD_DAIDTD, dad)
         call get_fill(engine, CD_CONGEL, cg)
         call get_fill(engine, CD_FRAZIL, fz)
         call get_fill(engine, CD_SNOICE, sn)
         call get_fill(engine, CD_MELTT, mt)
         call get_fill(engine, CD_MELTB, mb)
         call get_fill(engine, CD_MELTL, ml)
         call get_fill(engine, CD_EVAPI, ei)
         do j = ng + 1, ng + engine%grid%ny_phys
            do i = ng + 1, ng + engine%grid%nx_phys
               ! independent side: the volume change of the step from the states before and after
               ! (no wind: the ice is at rest, so it is all thermodynamic)
               ! the thermodynamic tendency departs from the whole-step change by exactly the dynamic one
               e1 = max(e1, abs(dv(i, j) - (engine%state%ice%cice%vice(i, j) - v0(i, j))*CMD/DT) - 1.000001_wp*abs(dd(i, j)))
               e2 = max(e2, abs(dd(i, j)))
               e3 = max(e3, abs(da(i, j) - (engine%state%ice%cice%aice(i, j) - a0(i, j))*CMD/DT) - 1.000001_wp*abs(dad(i, j)))
               ! the thermodynamic volume identity, to ROUND-OFF: the ice part of the evaporation (evapi, kg m-2 s-1 =
               ! evapi_ai[cm/day water]*1000/CMD) is the only term besides the five rates and the frazil
               ident = cg(i, j) + fz(i, j) + sn(i, j) - mt(i, j) - mb(i, j) - ml(i, j) &
                       + ei(i, j)*1000.0_wp/engine%state%ice%cice%p%c%rhoi
               resid = dv(i, j) - ident
               worst = max(worst, abs(resid))
               if (abs(resid) > 1.0e-9_wp*max(1.0_wp, abs(ident))) bad_ident = .true.
               sgrow = sgrow + max(cg(i, j), 0.0_wp) + max(fz(i, j), 0.0_wp)
               smelt = smelt + mt(i, j) + mb(i, j)
            end do
         end do
      end do
      bad_state = e1 > 1.0e-7_wp .or. e3 > 1.0e-4_wp
      ! the ice is at rest to drift noise: the dynamic tendency is negligible (1e-7 in the non-mixed kinds)
      bad_state = bad_state .or. e2 > merge(2.0_wp, 1.0e-7_wp, kind == "mixed")
      if (.not. allocated(error)) then
         write (*, '(a,3es10.3)') "  "//label//": (|dvidtt - dV/dt| - |dvidtd|)+, |dvidtd|, |daidtt - dA/dt| max =", e1, e2, e3
         write (*, '(a,a,a,es10.3,a,es10.3,a,es10.3)') "  ", label, ": worst identity residual (round-off)", worst, &
            " cm/day;  growth sum", sgrow, " melt sum", smelt
         call check(error, .not. bad_state, label//": dvidtt is the independent state change up to the (negligible) dynamic tendency; daidtt likewise")
      end if
      if (.not. allocated(error)) call check(error, .not. bad_ident, label//": dvidtt closes the volume identity (rates + frazil + evapi/rhoi) to round-off")
      vmax = 0.0_wp
      if (.not. allocated(error)) then
         if (kind == "mixed") then
            call check(error, sgrow > 0.0_wp .and. smelt > 0.0_wp, label//": growth AND melt both present")
         end if
      end if
      call finish_engine(engine)
      call rm_dir("tmp_cice_diags_o*")
   end subroutine budget_case

   subroutine test_budget_thermo(error)
      type(error_type), allocatable, intent(out) :: error
      call budget_case(error, "growth", "growth")
      if (allocated(error)) return
      call budget_case(error, "melt", "melt")
      if (allocated(error)) return
      call budget_case(error, "mixed", "mixed")
   end subroutine test_budget_thermo

   subroutine test_budget_dynamics(error)
      !! Dynamics, transport and ridging on: the thermodynamic and dynamic tendencies add up to the whole step.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      logical :: ok
      real(wp), allocatable :: dvt(:, :), dvd(:, :), v0(:, :), dat(:, :), dad(:, :), a0(:, :)
      real(wp) :: t, dmax, errmax
      integer :: n, i, j, ng
      character(len=:), allocatable :: nml

      nml = case_nml("growth", "tmp_cice_diags_o1", .true.)
      call make_engine(engine, cfg, nml, ok)
      call check(error, ok, "engine setup failed")
      if (allocated(error)) return
      ng = engine%grid%nghost
      t = 0.0_wp
      dmax = 0.0_wp
      errmax = 0.0_wp
      do n = 1, 30
         call pull(engine)
         v0 = engine%state%ice%cice%vice
         a0 = engine%state%ice%cice%aice
         call one_step(engine, cfg, t, ok)
         call check(error, ok, "step ran")
         if (allocated(error)) exit
         call pull(engine)
         call get_fill(engine, CD_DVIDTT, dvt)
         call get_fill(engine, CD_DVIDTD, dvd)
         call get_fill(engine, CD_DAIDTT, dat)
         call get_fill(engine, CD_DAIDTD, dad)
         do j = ng + 1, ng + engine%grid%ny_phys
            do i = ng + 1, ng + engine%grid%nx_phys
               errmax = max(errmax, abs(dvt(i, j) + dvd(i, j) - (engine%state%ice%cice%vice(i, j) - v0(i, j))*CMD/DT))
               errmax = max(errmax, abs(dat(i, j) + dad(i, j) - (engine%state%ice%cice%aice(i, j) - a0(i, j))*CMD/DT)*1.0e-3_wp)
               dmax = max(dmax, abs(dvd(i, j)))
            end do
         end do
      end do
      if (.not. allocated(error)) then
         write (*, '(a,es10.3,a,es10.3)') "  dynamics case: max|dvidtd| =", dmax, " cm/day, telescoping error", errmax
         call check(error, errmax < 1.0e-7_wp, "dvidtt + dvidtd is the whole step's volume change")
      end if
      if (.not. allocated(error)) call check(error, dmax > 1.0e-3_wp, "the dynamic tendency is non-zero where the ice moves")
      call finish_engine(engine)
      call rm_dir("tmp_cice_diags_o*")
   end subroutine test_budget_dynamics

   subroutine test_units_signs(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      logical :: ok
      real(wp), allocatable :: cg(:, :), mt(:, :), mb(:, :), v0(:, :), ml(:, :), fz(:, :), sn(:, :), dv(:, :), ei(:, :)
      real(wp) :: t, hcm_day, mt_min, mb_min, sub_cm
      integer :: n, i, j, ng, nc
      character(len=:), allocatable :: nml

      ! growth column
      nml = case_nml("column", "tmp_cice_diags_o1", .false.)
      call make_engine(engine, cfg, nml, ok)
      call check(error, ok, "growth engine setup failed")
      if (allocated(error)) return
      ng = engine%grid%nghost
      t = 0.0_wp
      ! an EVEN number of steps in all: the thermo cadence counter is not reset between engines of one process
      call advance(engine, cfg, t, 5, ok)
      call pull(engine)
      v0 = engine%state%ice%cice%vice
      call one_step(engine, cfg, t, ok)
      call pull(engine)
      call get_fill(engine, CD_CONGEL, cg)
      i = ng + 8
      j = ng + 8
      ! per-day thickness change of the cell volume: all of it is congelation (cold air, no snow melt, no new ice)
      call get_fill(engine, CD_MELTT, mt)
      call get_fill(engine, CD_MELTB, mb)
      call get_fill(engine, CD_MELTL, ml)
      call get_fill(engine, CD_FRAZIL, fz)
      call get_fill(engine, CD_SNOICE, sn)
      call get_fill(engine, CD_DVIDTT, dv)
      ! bare ice, aice = 1: the volume changes by congelation, by the (small) basal/top melt and by ice
      ! sublimation `evap/rhoi` (evap < 0)
      hcm_day = (engine%state%ice%cice%vice(i, j) - v0(i, j))*CMD/DT
      call get_fill(engine, CD_EVAPI, ei)
      sub_cm = ei(i, j)*1000.0_wp/engine%state%ice%cice%p%c%rhoi
      call check(error, cg(i, j) > 0.0_wp, "congel > 0 for a growing column")
      if (.not. allocated(error)) then
         write (*, '(a,8es11.4)') "  growth column: congel, dV/dt, dvidtt, sublimation, meltt, meltb, meltl, frazil:", &
            cg(i, j), hcm_day, dv(i, j), sub_cm, mt(i, j), mb(i, j), ml(i, j), fz(i, j)
         call check(error, abs(dv(i, j) - hcm_day) <= 1.0e-9_wp, "dvidtt is the volume change per day")
         if (.not. allocated(error)) call check(error, abs(cg(i, j) - (hcm_day - sub_cm + mt(i, j) + mb(i, j) + ml(i, j) &
                                                                       - fz(i, j) - sn(i, j))) <= 1.0e-9_wp, &
                                                "congel is the volume change per day without the other terms")
      end if
      call finish_engine(engine)
      if (allocated(error)) return
      ! melt column
      nml = case_nml("melt", "tmp_cice_diags_o1", .false.)
      call make_engine(engine, cfg, nml, ok)
      call check(error, ok, "melt engine setup failed")
      if (allocated(error)) return
      t = 0.0_wp
      call advance(engine, cfg, t, 12, ok)
      call pull(engine)
      call get_fill(engine, CD_MELTT, mt)
      call get_fill(engine, CD_MELTB, mb)
      mt_min = huge(1.0_wp)
      mb_min = huge(1.0_wp)
      nc = 0
      do j = ng + 1, ng + engine%grid%ny_phys
         do i = ng + 1, ng + engine%grid%nx_phys
            mt_min = min(mt_min, mt(i, j))
            mb_min = min(mb_min, mb(i, j))
            nc = nc + 1
         end do
      end do
      write (*, '(a,2es11.4,a)') "  melt column: meltt, meltb", mt(ng + 8, ng + 8), mb(ng + 8, ng + 8), " cm/day"
      call check(error, mt_min >= 0.0_wp .and. mb_min >= 0.0_wp, "melt rates are never negative")
      if (.not. allocated(error)) call check(error, mt(ng + 8, ng + 8) > 0.0_wp, "meltt > 0 in warm air with sunshine")
      if (.not. allocated(error)) call check(error, mb(ng + 8, ng + 8) > 0.0_wp, "meltb > 0 over a warm ocean")
      call finish_engine(engine)
      call rm_dir("tmp_cice_diags_o*")
   end subroutine test_units_signs

   subroutine test_tsfc(error)
      !! `Tsfc` is the ICE-AREA mean of the surface temperature, the freezing point only where there is no ice.  The
      !! state is stamped with hand-computed cases (the expectations are literals, not the fill's expression):
      !! a cell of 0.3 melting ice at 0 degC must read 0, not a blend with the open water; two categories of
      !! different area and temperature; an ice-free cell reads the freezing point.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      logical :: ok
      real(wp), allocatable :: a(:, :)
      integer :: ng
      character(len=:), allocatable :: nml

      nml = case_nml("growth", "tmp_cice_diags_o1", .true.)
      call make_engine(engine, cfg, nml, ok)
      call check(error, ok, "engine setup failed")
      if (allocated(error)) return
      ng = engine%grid%nghost
      associate (c => engine%state%ice%cice)
         c%aicen = 0.0_wp
         c%Tsfcn = -5.0_wp
         ! (ng+5, ng+5): 0.3 of melting ice at 0 degC in category 1
         c%aicen(ng + 5, ng + 5, 1) = 0.3_wp
         c%Tsfcn(ng + 5, ng + 5, 1) = 0.0_wp
         ! (ng+6, ng+5): 0.1 at -10 and 0.2 at -2 degC: (-1 - 0.4)/0.3 = -14/3
         c%aicen(ng + 6, ng + 5, 1) = 0.1_wp
         c%Tsfcn(ng + 6, ng + 5, 1) = -10.0_wp
         c%aicen(ng + 6, ng + 5, 2) = 0.2_wp
         c%Tsfcn(ng + 6, ng + 5, 2) = -2.0_wp
         ! (ng+7, ng+5): no ice at all
         !$acc update device(c%aicen, c%Tsfcn)
      end associate
      engine%state%icec%work%Tf = -1.8_wp
      !$acc update device(engine%state%icec%work%Tf)
      call get_fill(engine, CD_TSFC, a)
      call check(error, abs(a(ng + 5, ng + 5) - 0.0_wp) < 1.0e-14_wp, "0.3 of melting ice reads 0 degC (no open-water blend)")
      if (.not. allocated(error)) call check(error, abs(a(ng + 6, ng + 5) + 14.0_wp/3.0_wp) < 1.0e-14_wp, &
                                             "two categories: the ice-area mean -14/3")
      if (.not. allocated(error)) call check(error, abs(a(ng + 7, ng + 5) + 1.8_wp) < 1.0e-14_wp, &
                                             "no ice: the freezing point")
      call finish_engine(engine)
      call rm_dir("tmp_cice_diags_o*")
   end subroutine test_tsfc

   subroutine test_engine_reuse(error)
      !! One engine object, an ODD-length run, then a second engine on the same object: the cadence counter
      !! restarts with the run (it used to leak, and the ice hand-off then refused the window, status 104).
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      logical :: ok
      real(wp) :: t
      character(len=:), allocatable :: nml

      nml = case_nml("growth", "tmp_cice_diags_o1", .true.)
      call make_engine(engine, cfg, nml, ok)
      call check(error, ok, "first engine setup failed")
      if (allocated(error)) return
      t = 0.0_wp
      call advance(engine, cfg, t, 3, ok)
      call check(error, ok, "3 steps ran")
      call finish_engine(engine)
      if (allocated(error)) return
      call make_engine(engine, cfg, nml, ok)
      call check(error, ok, "second engine setup failed")
      if (allocated(error)) return
      call check(error, engine%state%dyn%outer_step_count == 0, "the cadence counter restarts with the run")
      t = 0.0_wp
      if (.not. allocated(error)) then
         call advance(engine, cfg, t, 5, ok)
         call check(error, ok, "5 steps of the reused engine ran (an odd count)")
      end if
      call finish_engine(engine)
      call rm_dir("tmp_cice_diags_o*")
   end subroutine test_engine_reuse

   subroutine test_melts(error)
      !! `melts` against the independent snow-volume change: snow-covered ice under warm air, no snowfall, the
      !! ice at rest.  Over the warm half the cell-mean snow loss per day equals `melts` (the snow
      !! sublimation `evaps` and the snow-ice formation are added back from their own diagnostics).
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      logical :: ok
      real(wp), allocatable :: ms(:, :), es(:, :), sn(:, :), s0(:, :)
      real(wp) :: t, err, tot
      integer :: n, i, j, ng
      character(len=:), allocatable :: nml

      nml = replace(case_nml("mixed", "tmp_cice_diags_o1", .false.), "h_snow = 0.03", "h_snow = 0.15")
      call make_engine(engine, cfg, nml, ok)
      call check(error, ok, "engine setup failed")
      if (allocated(error)) return
      call warm_east(engine)
      ng = engine%grid%nghost
      t = 0.0_wp
      err = 0.0_wp
      tot = 0.0_wp
      do n = 1, 8
         call pull(engine)
         !$acc update self(engine%state%ice%cice%vsno)
         s0 = engine%state%ice%cice%vsno
         call one_step(engine, cfg, t, ok)
         call check(error, ok, "step ran")
         if (allocated(error)) exit
         !$acc update self(engine%state%ice%cice%vsno)
         call get_fill(engine, CD_MELTS, ms)
         call get_fill(engine, CD_EVAPS, es)
         call get_fill(engine, CD_SNOICE, sn)
         do j = ng + 1, ng + engine%grid%ny_phys
            do i = ng + 1 + engine%grid%nx_phys/2 + 2, ng + engine%grid%nx_phys - 2
               ! cm/day of snow: -dVs/dt = melts - evaps_snow/rhos + snoice*rhoi/rhos (snow converted)
               err = max(err, abs(-(engine%state%ice%cice%vsno(i, j) - s0(i, j))*CMD/DT &
                                  - (ms(i, j) - es(i, j)*1000.0_wp/engine%state%ice%cice%p%c%rhos &
                                     + sn(i, j)*engine%state%ice%cice%p%c%rhoi/engine%state%ice%cice%p%c%rhos)))
               tot = tot + ms(i, j)
            end do
         end do
      end do
      if (.not. allocated(error)) then
         write (*, '(a,es10.3,a,es10.3)') "  melts: snow budget error", err, " cm/day (sum melts", tot
         call check(error, tot > 0.0_wp, "snow melts over the warm half")
      end if
      ! the snow budget closes except for the ice motion of the (slowly drifting) ice: bound by the dynamic share
      if (.not. allocated(error)) call check(error, err < 0.1_wp, "melts closes the snow budget (drift noise only)")
      call finish_engine(engine)
      call rm_dir("tmp_cice_diags_o*")
   end subroutine test_melts

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
            if (snap%n >= MAXF) error stop "test_cice_diags: raise MAXF"
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

   subroutine snapshots_equal(sa, sb, nbad, first)
      type(snapshot_t), intent(in) :: sa, sb
      integer, intent(out) :: nbad
      character(len=*), intent(out) :: first
      integer :: e, i, j, k

      nbad = 0
      first = ""
      if (sa%n /= sb%n) then
         nbad = -1
         first = "field count"
         return
      end if
      do e = 1, sa%n
         if (sa%f(e)%tag /= sb%f(e)%tag .or. any(shape(sa%f(e)%a) /= shape(sb%f(e)%a))) then
            nbad = nbad + 1
            if (len_trim(first) == 0) first = trim(sa%f(e)%tag)//" (tag or shape)"
            cycle
         end if
         do k = 1, size(sa%f(e)%a, 3)
            do j = 1, size(sa%f(e)%a, 2)
               do i = 1, size(sa%f(e)%a, 1)
                  if (transfer(sa%f(e)%a(i, j, k), 0_int64) /= transfer(sb%f(e)%a(i, j, k), 0_int64)) then
                     nbad = nbad + 1
                     if (len_trim(first) == 0) first = trim(sa%f(e)%tag)
                  end if
               end do
            end do
         end do
      end do
   end subroutine snapshots_equal

   subroutine test_observer(error)
      !! Engine A: no diagnostic stream.  Engine B: the stream on, every CICE name registered and folded
      !! every step by the diagnostic manager.  The restart registry must agree to the last bit.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: ea, eb
      type(config_t) :: cfg_a, cfg_b
      type(snapshot_t) :: sa, sb
      real(wp) :: t_a, t_b
      logical :: ok
      integer :: nbad
      character(len=64) :: first
      character(len=:), allocatable :: nml_a, nml_b

      nml_a = case_nml("mixed", "", .true.)
      nml_b = case_nml("mixed", "tmp_cice_diags_obs", .true.)
      call make_engine(ea, cfg_a, nml_a, ok)
      call check(error, ok, "engine A setup failed")
      if (allocated(error)) return
      call make_engine(eb, cfg_b, nml_b, ok)
      call check(error, ok, "engine B setup failed")
      if (allocated(error)) then
         call finish_engine(ea)
         return
      end if
      call warm_east(ea)
      call warm_east(eb)
      call check(error, eb%diag_enabled .and. eb%state%diag%nvars > 20, "engine B registered the diagnostics")
      t_a = 0.0_wp
      t_b = 0.0_wp
      if (.not. allocated(error)) call advance(ea, cfg_a, t_a, 24, ok)
      if (.not. allocated(error)) call advance(eb, cfg_b, t_b, 24, ok)
      if (.not. allocated(error)) then
         call take_snapshot(ea, sa)
         call take_snapshot(eb, sb)
         call snapshots_equal(sa, sb, nbad, first)
         write (*, '(a,i0,a)') "  observer: ", sa%n, " registry fields compared"
         call check(error, nbad == 0, "the restart registry differs with the history diagnostics on (first: "//trim(first)//")")
         if (.not. allocated(error)) call check(error, sa%n > 30, "a real registry was compared")
      end if
      call finish_engine(ea)
      call finish_engine(eb)
      call rm_dir("tmp_cice_diags_obs")
   end subroutine test_observer

   subroutine test_restart_fills(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: ea, eb
      type(config_t) :: cfg_a, cfg_b
      real(wp) :: t_a, t_b
      real(wp), allocatable :: fa(:, :), fb(:, :)
      logical :: ok
      integer :: ierr, step_b, id, i, j, ng, nbad
      character(len=*), parameter :: FN = "test_cice_diags_rt.nc"
      character(len=:), allocatable :: nml, nml_b

      nml = case_nml("mixed", "tmp_cice_diags_o1", .true.)
      nml_b = replace(nml, "tmp_cice_diags_o1", "tmp_cice_diags_o2")
      call make_engine(ea, cfg_a, nml, ok)
      call check(error, ok, "engine A setup failed")
      if (allocated(error)) return
      call warm_east(ea)
      t_a = 0.0_wp
      call advance(ea, cfg_a, t_a, 9, ok)      ! odd: the checkpoint falls inside a thermodynamics window
      call check(error, ok, "engine A failed before the checkpoint")
      if (allocated(error)) return
      call ocean_state_restart_write(ea%state, ea%grid, ea%decomp, FN, t_a, 9, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "checkpoint write failed")
      if (allocated(error)) return
      call advance(ea, cfg_a, t_a, 5, ok)
      call check(error, ok, "engine A failed after the checkpoint")
      if (allocated(error)) return
      call make_engine(eb, cfg_b, nml_b, ok, restart_file=FN, t0=t_b, step0=step_b)
      call check(error, ok, "engine B (warm restart) setup failed")
      if (allocated(error)) return
      call warm_east(eb)
      call advance(eb, cfg_b, t_b, 5, ok)
      call check(error, ok, "engine B failed after the restart")
      nbad = 0
      ng = ea%grid%nghost
      if (.not. allocated(error)) then
         call pull(ea)
         call pull(eb)
         do id = 1, NID
            call get_fill(ea, id, fa)
            call get_fill(eb, id, fb)
            do j = ng + 1, ng + ea%grid%ny_phys
               do i = ng + 1, ng + ea%grid%nx_phys
                  if (transfer(fa(i, j), 0_int64) /= transfer(fb(i, j), 0_int64)) nbad = nbad + 1
               end do
            end do
         end do
         call check(error, nbad == 0, "every fill continues bitwise across a checkpoint")
      end if
      call finish_engine(ea)
      call finish_engine(eb)
      call delete_file(FN)
      call rm_dir("tmp_cice_diags_o*")
   end subroutine test_restart_fills

   subroutine test_cost_report(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=16) :: v
      integer :: stat, n, run
      real(wp) :: ms(2)
      character(len=:), allocatable :: nml
      character(len=*), parameter :: GRID_SMALL = "nx = 24, ny = 16"
      character(len=*), parameter :: GRID_BIG = "nx = 512, ny = 512"

      call get_environment_variable("RDB_DIAGS_COST", v, status=stat)
      if (stat /= 0) return
      n = 40
      do run = 1, 2
         if (run == 1) nml = replace(case_nml("mixed", "", .true.), GRID_SMALL, GRID_BIG)
         if (run == 2) nml = replace(case_nml("mixed", "tmp_cice_diags_cost", .true.), GRID_SMALL, GRID_BIG)
         if (run == 2) nml = replace(nml, "dt_out = 1200.0", "dt_out = 864000.0")
         call time_steps(nml, n, ms(run), error)
         if (allocated(error)) return
      end do
      write (*, '(a,2f10.3,a,f9.3,a)') "  cost (512x512, 5 cat): ms/step without, with the history set:", ms(1), ms(2), &
         "  difference", ms(2) - ms(1), " ms/step"
      call rm_dir("tmp_cice_diags_cost")
   end subroutine test_cost_report

   subroutine time_steps(nml, n, ms_per_step, error)
      character(len=*), intent(in) :: nml
      integer, intent(in) :: n
      real(wp), intent(out) :: ms_per_step
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      real(wp) :: t
      real(wp) :: c0, c1
      logical :: ok

      ms_per_step = 0.0_wp
      call make_engine(engine, cfg, nml, ok)
      call check(error, ok, "cost engine setup failed")
      if (allocated(error)) return
      t = 0.0_wp
      call advance(engine, cfg, t, 6, ok)
      c0 = wall_seconds()
      call advance(engine, cfg, t, n, ok)
      c1 = wall_seconds()
      call check(error, ok, "cost steps ran")
      ms_per_step = 1.0e3_wp*(c1 - c0)/real(n, wp)
      call finish_engine(engine)
   end subroutine time_steps

   subroutine delete_file(fn_)
      character(len=*), intent(in) :: fn_
      integer :: u, stat
      open (newunit=u, file=fn_, status="old", iostat=stat)
      if (stat == 0) close (u, status="delete")
   end subroutine delete_file

   subroutine rm_dir(d)
      character(len=*), intent(in) :: d
      call execute_command_line("rm -rf "//trim(d))
   end subroutine rm_dir

end module test_cice_diags
