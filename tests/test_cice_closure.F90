!! Closure of the CICE-port sea ice (`[ocean.ice] scheme = "cice"`): the named-term conservation ledger
!! (`rdb_ice_ledger_terms`) closes stage by stage, what the ice hands over is what the ocean receives, and
!! ocean + ice together conserve salt, water and heat against the atmosphere that was delivered (D11.5).
!! Everything runs through the production engine path (`engine_step` -> `engine_step_ice` -> `engine_step_finalize`).
!!
!! Cases:
!!   * A `column_growth`, `column_melt`, `column_snow`, `column_newice`, `column_lateral`: a doubly periodic uniform
!!     box (every cell is the same column) in each thermodynamic regime; every line (water, salt, heat) of every
!!     stage closes to round-off (heat in `therm1`: the vertical solver's tolerance), the regime's own named term is
!!     active, and the new-ice energy is `-rhoi Lfresh V`, which equals minus the bank the ocean paid (the SIGN).
!!   * B `interface_uniform`, `interface_prescribed`: the window means of what the ice handed over equal what the
!!     ocean's assembler holds, cell by cell (heat, shortwave, virtual salt), read from the ocean's arrays.
!!   * C `global_uniform`, `global_prescribed_mixed`: a CLOSED basin, ocean + ice, freezing and melting and moving
!!     ice: total salt, water and heat against the delivered atmosphere, with the known non-delivered terms as
!!     named, quantified residuals.
!!   * D `dynamics_only`, `dynamics_only_snow`: thermodynamics off; EVP + transport + ridging conserve the
!!     area-weighted volume, snow, salt and energy; ocean and atmosphere lines are zero (with snow, the ridged
!!     half of it is the one named exchange).
!!   * E `restart_ledger`: N + checkpoint + N == 2N, the ledger accumulators bitwise; `restart_ledger_off_on` and
!!     `restart_ledger_on_off`: a checkpoint written with the ledger off is read with it on and the reverse.
!!   * F `zap_booked`: the ITD clean-up's zaps (a category below `puny`, the round-off excess `aice > 1`, a near-empty
!!     cell) are the `cleanup` stage's ocean-ward flux: each is non-zero, has the flux slot's signs (mass and salt
!!     positive, heat negative), equals the amount an independent host expectation removes, and the stage closes.
!!     `clip_booked`: the remap transport's clip record (content dropped where the new area is 0) is a named source of
!!     the `transport` stage: with it booked the stage closes, with it withheld the residual is the whole drop.
!!     `snow_reset_reported`: the energy the profile initialisation adds by resetting snow thinner than `hs_min` is
!!     reported (`tab%snow_reset`), equals the independent host value, and accounts for the `therm1` heat residual
!!     to the solver's bound without being subtracted from it.
!!   * `ledger_on_equals_off`: the ledger is an OBSERVER: the whole restart registry (excluding `ice_led_*`)
!!     is bitwise the same after N steps with it on and off.  `ledger_nan_sticky`: a non-finite step poisons the
!!     accumulators and the console for good.
!!
!! What the tests can and cannot see.  The ledger's `atm_rad` term is built from the ice's OWN `fswabs`,
!! `flwout`, `fsens` (and `evap`), so a sign error applied consistently there would close every line; only `flw`,
!! `fsnow`, `frain` are forcing inputs.  `column_melt` therefore also ties `fswabs` to the forcing's downward
!! shortwave and the albedo range independently.  The vertical solver's error (`ICE_FERRMAX`, a per-category flux
!! tolerance) is a BOUND on the heat closure, never a subtracted term.
module test_cice_closure
   use, intrinsic :: iso_fortran_env, only: int32, int64, output_unit
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_config, only: config_t, read_config_from_toml_string
   use rdb_ocean_engine, only: ocean_engine_t, engine_setup, engine_enter_data, &
                               engine_step, engine_step_ice, engine_step_finalize, &
                               engine_exit_data, engine_teardown
   use rdb_ocean_state, only: ocean_state_restart_write, ocean_state_build_restart_registry
   use rdb_ocean_restart, only: restart_registry_t
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   use rdb_comm_env, only: comm_env_init, comm_env_setup_roles
   use rdb_ice_errors, only: ICE_ERR_FATAL_MASK, ICE_ERR_THERMO_CONSERV
   use rdb_ice_thickness, only: ICE_FERRMAX
   use rdb_efp, only: efp_t
   use rdb_ice_remap_flux, only: remap_ledger_cells_impl
   use rdb_ice_transport_pack, only: transport_pack
   use rdb_ice_flux_state, only: ice_flux_zero
   use rdb_ice_ledger_terms, only: ice_lt_ix_state, ICE_LT_NCLIP, ICE_LT_THERM1, ICE_LT_REMAP, ICE_LT_NEWICE, ICE_LT_LATERAL, ICE_LT_CLEANUP, &
                                   ICE_LT_DYNAMICS, ICE_LT_TRANSPORT, ICE_LT_RIDGING, ICE_LT_NSTAGE, &
                                   ICE_LT_WATER, ICE_LT_SALT, ICE_LT_HEAT, ice_lt_stage_name, ice_lt_line_name
   use rdb_ocean_console_stats, only: compute_total_tracer
   use rdb_ice_cice_hooks, only: ice_cice_ledger_report
   use rdb_ice_vprofile, only: ice_init_vertical_profile
   use test_cice_synth_atm, only: write_synth_atm, remove_synth_atm
   implicit none
   private

   public :: collect_cice_closure_tests

   real(wp), parameter :: DT = 600.0_wp
   character(len=*), parameter :: FN_ATM = "test_cice_closure_atm.nc"
   character(len=*), parameter :: FN_LST = "test_cice_closure_atm.txt"
   character(len=*), parameter :: FN_RT = "test_cice_closure_rt.nc"
   integer, parameter :: MAXF = 300
   logical :: comm_inited = .false.

   type :: field_t
      character(len=64) :: tag = ""
      real(wp), allocatable :: a(:, :, :)
   end type field_t

   type :: snapshot_t
      integer :: n = 0
      type(field_t) :: f(MAXF)
   end type snapshot_t

   type :: stage_run_t
      !! What a run's per-step tables add up to.
      real(wp) :: rel(3, ICE_LT_NSTAGE) = 0.0_wp
         !! Worst `|res|/scale` per line and stage over the steps.
      real(wp) :: dstate(3, ICE_LT_NSTAGE) = 0.0_wp
      real(wp) :: fout(3, ICE_LT_NSTAGE) = 0.0_wp
      real(wp) :: absfout(3, ICE_LT_NSTAGE) = 0.0_wp
      real(wp) :: snow = 0.0_wp, rain = 0.0_wp, evap = 0.0_wp, lw_in = 0.0_wp, snow_enth = 0.0_wp
      real(wp) :: atm_rad = 0.0_wp, atm_lat = 0.0_wp, frz_lat = 0.0_wp
      real(wp) :: fsw(ICE_LT_NSTAGE) = 0.0_wp
      real(wp) :: rsum(3) = 0.0_wp
         !! Signed residual of each line summed over the steps and stages (the named `R`).
      real(wp) :: efrz_indep = 0.0_wp
         !! The bank energy the ocean paid for new ice, from the ocean's side only.
      real(wp) :: frz_mismatch = 0.0_wp
         !! Worst `|frz_lat + E_frz|/max(|frz_lat|, ...)` per step.
      real(wp) :: newice_state_mismatch = 0.0_wp
      real(wp) :: hb0 = 0.0_wp, hb_delta = 0.0_wp, stock_now = 0.0_wp
         !! The ocean's frazil-clamp accumulator [K m summed with the cell areas] at the start and its change; the
         !! bank standing [J].
      real(wp) :: clamp_mismatch = 0.0_wp
         !! `|clamp heat - bank standing - new-ice latent heat|`, relative.
      real(wp) :: asm_mismatch = 0.0_wp
         !! Worst relative difference between the assembler's open-water heat and the bulk forcing's.
      real(wp) :: fswabs_sum = 0.0_wp, sw_in_sum = 0.0_wp
         !! `dt sum A fswabs` and `dt sum A aice_init sw_down`.
      real(wp) :: bl99_ratio = 0.0_wp
         !! Worst per-step `|res(heat, therm1)|` over `1.1 ICE_FERRMAX dt sum(A aice_init)`: the vertical solver's
         !! tolerance (`ICE_FERRMAX` [W m-2], Bitz and Lipscomb 1999) per unit ice-covered area.
      logical :: finite = .true.
      integer :: nsteps = 0
   end type stage_run_t

contains

   subroutine collect_cice_closure_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      character(len=64) :: only
      logical, allocatable :: keep(:)
      integer :: k
      testsuite = [ &
                  new_unittest("column_growth", test_column_growth), &
                  new_unittest("column_melt", test_column_melt), &
                  new_unittest("column_snow", test_column_snow), &
                  new_unittest("column_newice", test_column_newice), &
                  new_unittest("column_lateral", test_column_lateral), &
                  new_unittest("interface_uniform", test_interface_uniform), &
                  new_unittest("interface_prescribed", test_interface_prescribed), &
                  new_unittest("global_uniform", test_global_uniform), &
                  new_unittest("global_prescribed_mixed", test_global_prescribed), &
                  new_unittest("global_prescribed_frazil", test_global_prescribed_frazil), &
                  new_unittest("dynamics_only", test_dynamics_only), &
                  new_unittest("dynamics_only_snow", test_dynamics_only_snow), &
                  new_unittest("restart_ledger", test_restart_ledger), &
                  new_unittest("ledger_off_inert", test_ledger_off_inert), &
                  new_unittest("ledger_on_equals_off", test_ledger_on_equals_off), &
                  new_unittest("ledger_nan_sticky", test_ledger_nan_sticky), &
                  new_unittest("zap_booked", test_zap_booked), &
                  new_unittest("clip_booked", test_clip_booked), &
                  new_unittest("snow_reset_reported", test_snow_reset_reported), &
                  new_unittest("restart_ledger_off_on", test_restart_off_on), &
                  new_unittest("restart_ledger_on_off", test_restart_on_off), &
                  new_unittest("cice_coupled_water_closure", test_mass_water_closure), &
                  new_unittest("cice_coupled_water_closure_prescribed", test_mass_water_closure_prescribed), &
                  new_unittest("cice_salt_closure_mass", test_mass_salt_closure), &
                  new_unittest("cice_heat_closure_mass", test_mass_heat_closure), &
                  new_unittest("cice_closure_mass_ratio8_deferred", test_mass_closure_ratio8), &
                  new_unittest("rain_split_sums_to_total", test_rain_split), &
                  new_unittest("ice_load_no_flow", test_ice_load_no_flow), &
                  new_unittest("ice_load_drift", test_ice_load_drift), &
                  new_unittest("psurf_ice_only_seam", test_psurf_ice_only), &
                  new_unittest("psurf_atm_load_ice", test_psurf_atm_load_ice), &
                  new_unittest("psurf_atm_load_ib", test_psurf_atm_load_ib), &
                  new_unittest("ice_load_restart_exact", test_ice_load_restart) &
                  ]
      ! Dev filter: `RDB_CICE_CLOSURE_ONLY=<substring>` keeps the cases whose name contains it.
      call get_environment_variable("RDB_CICE_CLOSURE_ONLY", only)
      if (len_trim(only) > 0) then
         keep = [(index(testsuite(k)%name, trim(only)) > 0, k=1, size(testsuite))]
         testsuite = pack(testsuite, keep)
      end if
   end subroutine collect_cice_closure_tests

   ! ------------------------------------------------------------------
   ! inputs and engine plumbing
   ! ------------------------------------------------------------------

   pure function has_nan(line) result(yes)
      !! The text holds a NaN (any spelling the runtime writes).
      character(len=*), intent(in) :: line
      logical :: yes
      integer :: k

      yes = .false.
      do k = 1, len(line) - 2
         if (line(k:k + 2) == "NaN" .or. line(k:k + 2) == "NAN" .or. line(k:k + 2) == "nan") yes = .true.
      end do
   end function has_nan

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

   function uni_nml(geom) result(nml)
      !! The uniform-atmosphere cases, ledger on.  `column`: doubly periodic 8x8, uniform, no wind.
      !! `basin`: closed 24x16 with wind and an f-plane.
      character(len=*), intent(in) :: geom
      character(len=:), allocatable :: nml
      character(len=*), parameter :: NL = new_line("a")
      character(len=:), allocatable :: bc, phys, grid, ic, atm

      bc = "[ocean.bc]"//NL// &
           'west = "periodic"'//NL// &
           'east = "periodic"'//NL// &
           'south = "periodic"'//NL// &
           'north = "periodic"'//NL
      phys = "[physics]"//NL// &
             "coriolis_f = 1.4e-4"//NL// &
             "wind_stress_x = 0.0"//NL// &
             "wind_stress_y = 0.0"//NL
      grid = "[grid]"//NL// &
             "nx = 8"//NL// &
             "ny = 8"//NL// &
             "nghost = 3"//NL// &
             "dx = 20000.0"//NL// &
             "dy = 20000.0"//NL
      ic = "[ocean.ice_ic]"//NL// &
           'conc_config = "uniform"'//NL// &
           "h_ice = 0.5"//NL// &
           "conc = 1.0"//NL// &
           "h_snow = 0.0"//NL// &
           "t_ice = -10.0"//NL
      atm = "[ocean.ice_atm]"//NL// &
            "uni_tas = 250.0"//NL// &
            "uni_rlds = 150.0"//NL// &
            "uni_uas = 3.0"//NL// &
            "uni_vas = 0.0"//NL
      select case (geom)
      case ("column")
         continue
      case ("basin")
         bc = "[ocean.bc]"//NL// &
              'west = "wall"'//NL// &
              'east = "wall"'//NL// &
              'south = "wall"'//NL// &
              'north = "wall"'//NL
         phys = "[physics]"//NL// &
                "coriolis_f = 1.0e-4"//NL// &
                "wind_stress_x = 0.1"//NL// &
                "wind_stress_y = 0.03"//NL
         grid = "[grid]"//NL// &
                "nx = 24"//NL// &
                "ny = 16"//NL// &
                "nghost = 3"//NL// &
                "dx = 20000.0"//NL// &
                "dy = 20000.0"//NL
         ic = "[ocean.ice_ic]"//NL// &
              'conc_config = "uniform"'//NL// &
              "h_ice = 1.0"//NL// &
              "conc = 0.9"//NL// &
              "h_snow = 0.1"//NL// &
              "t_ice = -8.0"//NL
         atm = "[ocean.ice_atm]"//NL// &
               "uni_tas = 255.0"//NL// &
               "uni_rlds = 180.0"//NL// &
               "uni_uas = 12.0"//NL// &
               "uni_vas = 3.0"//NL
      case default
         error stop "test_cice_closure: unknown geometry"
      end select
      nml = "[sim]"//NL// &
            'sim_type = "ocean"'//NL// &
            "[time]"//NL// &
            "t_end = 86400.0"//NL// &
            "dt_fixed = 600.0"//NL// &
            "[nonhydrostatic]"//NL// &
            "nz_layers = 3"//NL// &
            "[tracer]"//NL// &
            "initial_temperature = -1.0"//NL// &
            "initial_salinity = 34.0"//NL// &
            "T_init_surface = -1.0"//NL// &
            "T_init_bottom = -1.0"//NL// &
            "[ocean.bt]"//NL// &
            "auto_n_inner = true"//NL// &
            "[ocean.hvisc]"//NL// &
            "nu_h = 200.0"//NL// &
            "[ocean.diag]"//NL// &
            "enabled = false"//NL// &
            grid//NL// &
            phys//NL// &
            "[vcoord]"//NL// &
            'vcoord_type = "sigma"'//NL// &
            "[ocean.topo]"//NL// &
            'topo_config = "flat"'//NL// &
            "max_depth = 100.0"//NL// &
            "[ocean.thermo]"//NL// &
            "enable_thermodynamics = true"//NL// &
            "[ocean.forcing]"//NL// &
            "enable_components = true"//NL// &
            "[ocean.vmix]"//NL// &
            "dt_therm_ratio = 1"//NL// &
            "[ocean.ice]"//NL// &
            "enable = true"//NL// &
            'scheme = "cice"'//NL// &
            "ncat = 5"//NL// &
            "nilyr = 4"//NL// &
            "hfreeze = 10.0"//NL// &
            "dynamics = true"//NL// &
            "transport = true"//NL// &
            "ledger = true"//NL// &
            "[ocean.ice_thermo]"//NL// &
            "update_ocn_f = true"//NL// &
            "[ocean.ice_dyn]"//NL// &
            'advection = "upwind"'//NL// &
            "ndte = 40"//NL// &
            ic//NL// &
            atm//NL// &
            bc//NL// &
            "[output]"//NL// &
            "output_to_file = false"//NL
   end function uni_nml

   function pres_nml(ratio) result(nml)
      !! A closed spherical high-latitude basin (12 x 10 cells of 0.25 degree) under the synthetic prescribed
      !! atmosphere (bulk forcing on), the ice 0.6 covered; ledger on.
      integer, intent(in) :: ratio
      character(len=:), allocatable :: nml
      character(len=*), parameter :: NL = new_line("a")
      character(len=32) :: sratio

      write (sratio, '(i0)') ratio
      nml = "[sim]"//NL// &
            'sim_type = "ocean"'//NL// &
            "[time]"//NL// &
            "t_end = 86400.0"//NL// &
            "dt_fixed = 600.0"//NL// &
            'start_date = "2000-01-01"'//NL// &
            "[nonhydrostatic]"//NL// &
            "nz_layers = 3"//NL// &
            "[ocean.diag]"//NL// &
            "enabled = false"//NL// &
            "[tracer]"//NL// &
            "initial_temperature = -1.0"//NL// &
            "initial_salinity = 34.0"//NL// &
            "T_init_surface = -1.0"//NL// &
            "T_init_bottom = -1.0"//NL// &
            "[ocean.bt]"//NL// &
            "auto_n_inner = true"//NL// &
            "[ocean.hvisc]"//NL// &
            "nu_h = 200.0"//NL// &
            "[grid]"//NL// &
            "nx = 12"//NL// &
            "ny = 10"//NL// &
            "nghost = 3"//NL// &
            "dx = 0.25"//NL// &
            "dy = 0.25"//NL// &
            "[ocean.grid]"//NL// &
            'grid_config = "spherical"'//NL// &
            "lon_west = 20.0"//NL// &
            "lat_south = 70.0"//NL// &
            "rad_earth = 6.371e6"//NL// &
            'coriolis_scheme = "planetary"'//NL// &
            "[physics]"//NL// &
            "wind_stress_x = 0.0"//NL// &
            "wind_stress_y = 0.0"//NL// &
            "[vcoord]"//NL// &
            'vcoord_type = "sigma"'//NL// &
            "[ocean.topo]"//NL// &
            'topo_config = "flat"'//NL// &
            "max_depth = 100.0"//NL// &
            "[ocean.thermo]"//NL// &
            "enable_thermodynamics = true"//NL// &
            "[ocean.forcing]"//NL// &
            "enable_components = true"//NL// &
            "[ocean.vmix]"//NL// &
            "dt_therm_ratio = "//trim(sratio)//NL// &
            "[ocean.ice]"//NL// &
            "enable = true"//NL// &
            'scheme = "cice"'//NL// &
            "ncat = 5"//NL// &
            "nilyr = 4"//NL// &
            "hfreeze = 10.0"//NL// &
            "dynamics = true"//NL// &
            "transport = true"//NL// &
            "ledger = true"//NL// &
            "[ocean.ice_thermo]"//NL// &
            "update_ocn_f = true"//NL// &
            "[ocean.ice_dyn]"//NL// &
            'advection = "upwind"'//NL// &
            "ndte = 40"//NL// &
            "[ocean.ice_ic]"//NL// &
            'conc_config = "uniform"'//NL// &
            "h_ice = 1.0"//NL// &
            "conc = 0.6"//NL// &
            "h_snow = 0.1"//NL// &
            "t_ice = -8.0"//NL// &
            "[ocean.atm_forcing]"//NL// &
            "enable = true"//NL// &
            'sw_interp = "linear"'//NL// &
            'fluxes_regrid = "bilinear"'//NL// &
            'time_mode = "interannual"'//NL// &
            "wind_relative = false"//NL// &
            'filelist_tas = "'//FN_LST//'"'//NL// &
            'filelist_huss = "'//FN_LST//'"'//NL// &
            'filelist_uas = "'//FN_LST//'"'//NL// &
            'filelist_vas = "'//FN_LST//'"'//NL// &
            'filelist_psl = "'//FN_LST//'"'//NL// &
            'filelist_rlds = "'//FN_LST//'"'//NL// &
            'filelist_rsds = "'//FN_LST//'"'//NL// &
            'filelist_prra = "'//FN_LST//'"'//NL// &
            'filelist_prsn = "'//FN_LST//'"'//NL// &
            "[ocean.bulk_flux]"//NL// &
            "enable = true"//NL// &
            "[ocean.bc]"//NL// &
            'west = "wall"'//NL// &
            'east = "wall"'//NL// &
            'south = "wall"'//NL// &
            'north = "wall"'//NL// &
            "[output]"//NL// &
            "output_to_file = false"//NL
   end function pres_nml

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
      call read_config_from_toml_string(nml, cfg, ierr=ierr)
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

   subroutine one_step(engine, cfg, t, ok, was_thermo)
      type(ocean_engine_t), intent(inout) :: engine
      type(config_t), intent(inout) :: cfg
      real(wp), intent(inout) :: t
      logical, intent(out) :: ok
      logical, intent(out), optional :: was_thermo
      integer :: ierr

      ok = .false.
      call engine_step(engine, DT, t, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      if (present(was_thermo)) was_thermo = engine%state%dyn%is_thermo_step()
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
      type(ocean_engine_t), intent(inout) :: engine

      !$acc update self(engine%state%ice%cice%aicen, engine%state%ice%cice%vicen, engine%state%ice%cice%vsnon)
      !$acc update self(engine%state%ice%cice%aice, engine%state%ice%cice%vice, engine%state%ice%cice%vsno)
      !$acc update self(engine%state%ice%cice%aice0, engine%state%ice%cice%aice_init, engine%state%ice%cice%ice_err)
      !$acc update self(engine%state%ice%aice, engine%state%ice%vice, engine%state%ice%vsno, engine%state%ice%aice0)
      !$acc update self(engine%state%icec%flux%dardg1dt, engine%state%icec%work%niter)
   end subroutine pull_cice

   subroutine pull_ice_flux(engine)
      type(ocean_engine_t), intent(inout) :: engine

      !$acc update self(engine%state%icec%flux%fresh, engine%state%icec%flux%fsalt)
      !$acc update self(engine%state%icec%flux%fhocn, engine%state%icec%flux%fswthru)
      !$acc update self(engine%state%icec%work%frain, engine%state%icec%work%fsnow, engine%state%icec%work%flw)
      !$acc update self(engine%state%icec%oin%sss_ocn)
      !$acc update self(engine%state%ice%frazil_heat, engine%state%ice%frzmlt, engine%state%ice%heat_budget_frazil)
      !$acc update self(engine%state%icec%flux%fswabs, engine%state%icec%work%sw_down)
   end subroutine pull_ice_flux

   subroutine pull_sf(engine)
      type(ocean_engine_t), intent(inout) :: engine

      associate (sf => engine%state%surface_flux)
         !$acc update self(sf%heat_ice, sf%sw_ice, sf%salt_ice, sf%Q_heat, sf%Q_salt)
         !$acc update self(sf%q_sw, sf%q_lw, sf%q_lat, sf%q_sens, sf%q_fus, sf%evap, sf%lprec, sf%fprec)
         !$acc update self(sf%heat_content_lprec, sf%heat_content_fprec)
         !$acc update self(sf%fw_net, sf%seaice_melt, sf%Q_salt_tracer, sf%p_surf)
      end associate
   end subroutine pull_sf

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

   function ocean_volume_mass(engine) result(m)
      !! `rho0 sum(h A)` [kg] over the owned cells: the Boussinesq column mass (`rho0` converts the volume the
      !! layers carry to the mass the real-mass step spends).
      type(ocean_engine_t), intent(inout) :: engine
      real(wp) :: m
      integer :: i, j, k, ng

      !$acc update self(engine%state%multilayer%h_layer)
      ng = engine%grid%nghost
      m = 0.0_wp
      do k = 1, size(engine%state%multilayer%h_layer, 3)
         do j = ng + 1, ng + engine%grid%ny_phys
            do i = ng + 1, ng + engine%grid%nx_phys
               m = m + engine%state%metrics%areaT(i, j)*engine%state%multilayer%h_layer(i, j, k)
            end do
         end do
      end do
      m = m*engine%state%eos%rho0
   end function ocean_volume_mass

   subroutine delete_file(fn_)
      character(len=*), intent(in) :: fn_
      integer :: u, stat
      character(len=128) :: msg

      open (newunit=u, file=fn_, status="old", action="readwrite", iostat=stat, iomsg=msg)
      if (stat == 0) close (u, status="delete")
   end subroutine delete_file

   ! ------------------------------------------------------------------
   ! A: the column regimes
   ! ------------------------------------------------------------------

   function profile_repair(engine, d_ice, d_snow) result(delta)
      !! The energy the vertical thermodynamics' profile initialisation ADDS to the state before the step starts
      !! (its round-off repair of a layer at or above its melting point, and of snow too thin to carry a
      !! temperature), `sum A a_n (einit_n - sum_k q_k h_k)`.  The thermodynamics' own energy check starts from the
      !! repaired energy, so this part is invisible to it and appears in the ledger's `therm1` residual.  Host call on
      !! the pulled start-of-step state.
      type(ocean_engine_t), intent(in) :: engine
      real(wp), intent(out) :: d_ice, d_snow
      real(wp) :: delta
      integer :: i, j, n, k, ng, nl, ns
      real(wp) :: q0(8), qs0(8)
      real(wp), allocatable :: zq(:), zqs(:), zt(:), zts(:), tm(:), zs(:)
      real(wp) :: hin, hilyr, hsn, hslyr, einit, raw
      integer(int32) :: ie

      delta = 0.0_wp
      d_ice = 0.0_wp
      d_snow = 0.0_wp
      associate (c => engine%state%ice%cice)
         nl = c%nilyr
         ns = c%nslyr
         ng = c%nghost
         allocate (zq(nl), zt(nl), tm(nl), zs(nl), zqs(ns), zts(ns))
         do j = ng + 1, ng + c%ny_phys
            do i = ng + 1, ng + c%nx_phys
               do n = 1, c%ncat
                  if (c%aicen(i, j, n) <= c%p%c%puny) cycle
                  raw = 0.0_wp
                  do k = 1, nl
                     zq(k) = c%qicen(i, j, n, k)
                     q0(k) = zq(k)
                     zs(k) = c%sicen(i, j, n, k)
                     raw = raw + zq(k)*c%vicen(i, j, n)/real(nl, wp)
                  end do
                  do k = 1, ns
                     zqs(k) = c%qsnon(i, j, n, k)
                     qs0(k) = zqs(k)
                     raw = raw + zqs(k)*c%vsnon(i, j, n)/real(ns, wp)
                  end do
                  ie = 0_int32
                  call ice_init_vertical_profile(c%aicen(i, j, n), c%vicen(i, j, n), c%vsnon(i, j, n), nl, ns, c%p%c, &
                                                 c%p%thermo%l_brine, zs, zq, zqs, zt, zts, tm, hin, hilyr, hsn, hslyr, einit, ie)
                  if (ie /= 0_int32) cycle
                  do k = 1, nl
                     d_ice = d_ice + engine%state%ice_cgrid%tarea(i, j)*(zq(k) - q0(k))*c%vicen(i, j, n)/real(nl, wp)
                  end do
                  do k = 1, ns
                     d_snow = d_snow + engine%state%ice_cgrid%tarea(i, j)*(zqs(k) - qs0(k))*c%vsnon(i, j, n)/real(ns, wp)
                  end do
                  delta = delta + engine%state%ice_cgrid%tarea(i, j)*(c%aicen(i, j, n)*einit - raw)
               end do
            end do
         end do
      end associate
   end function profile_repair

   subroutine pull_cat_state(engine)
      type(ocean_engine_t), intent(inout) :: engine

      !$acc update self(engine%state%ice%cice%aicen, engine%state%ice%cice%vicen, engine%state%ice%cice%vsnon)
      !$acc update self(engine%state%ice%cice%qicen, engine%state%ice%cice%sicen, engine%state%ice%cice%qsnon)
   end subroutine pull_cat_state

   subroutine init_clamp(engine, run)
      !! The ocean's frazil-clamp accumulator at the start of a run (its change is the heat the clamp added).
      type(ocean_engine_t), intent(inout) :: engine
      type(stage_run_t), intent(inout) :: run
      integer :: i, j, ng

      !$acc update self(engine%state%ice%heat_budget_frazil)
      ng = engine%grid%nghost
      run%hb0 = 0.0_wp
      do j = ng + 1, ng + engine%grid%ny_phys
         do i = ng + 1, ng + engine%grid%nx_phys
            run%hb0 = run%hb0 + engine%state%ice_cgrid%tarea(i, j)*engine%state%ice%heat_budget_frazil(i, j, 1)
         end do
      end do
   end subroutine init_clamp

   subroutine accumulate_tab(engine, run)
      !! Fold the ledger's last step into `run`, and read the bank the ocean paid from the ocean's side.
      type(ocean_engine_t), intent(inout) :: engine
      type(stage_run_t), intent(inout) :: run
      integer :: q, s, i, j, ng, is, ie, js, je
      real(wp) :: efrz, r, a_dt, hb, stock

      associate (L => engine%state%icec%ledger)
         if (.not. L%valid .or. .not. L%finite) then
            run%finite = .false.
            return
         end if
         do s = 1, ICE_LT_NSTAGE
            do q = 1, 3
               r = abs(L%tab%res(q, s))/max(L%tab%scale(q), tiny(1.0_wp))
               run%rel(q, s) = max(run%rel(q, s), r)
               run%dstate(q, s) = run%dstate(q, s) + L%tab%dstate(q, s)
               run%fout(q, s) = run%fout(q, s) + L%tab%fout(q, s)
               run%absfout(q, s) = run%absfout(q, s) + abs(L%tab%fout(q, s))
               run%rsum(q) = run%rsum(q) + L%tab%res(q, s)
            end do
            run%fsw(s) = run%fsw(s) + L%tab%fsw(s)
         end do
         run%snow = run%snow + L%tab%snow
         run%rain = run%rain + L%tab%rain
         run%evap = run%evap + L%tab%evap
         run%lw_in = run%lw_in + L%tab%lw_in
         run%snow_enth = run%snow_enth + L%tab%snow_enth
         run%atm_rad = run%atm_rad + L%tab%atm_rad
         run%atm_lat = run%atm_lat + L%tab%atm_lat
         run%frz_lat = run%frz_lat + L%tab%frz_lat
         run%nsteps = run%nsteps + 1
         ! the bank the ocean paid this step: offered (frzmlt > 0, the bank at the import) minus what is left
         ng = engine%grid%nghost
         is = ng + 1
         ie = ng + engine%grid%nx_phys
         js = ng + 1
         je = ng + engine%grid%ny_phys
         efrz = 0.0_wp
         a_dt = 0.0_wp
         hb = 0.0_wp
         stock = 0.0_wp
         do j = js, je
            do i = is, ie
               hb = hb + engine%state%ice_cgrid%tarea(i, j)*engine%state%ice%heat_budget_frazil(i, j, 1)
               stock = stock + engine%state%ice_cgrid%tarea(i, j)*engine%state%ice%frazil_heat(i, j)
               run%fswabs_sum = run%fswabs_sum + DT*engine%state%ice_cgrid%tarea(i, j)*engine%state%icec%flux%fswabs(i, j)
               run%sw_in_sum = run%sw_in_sum + DT*engine%state%ice_cgrid%tarea(i, j)*engine%state%ice%cice%aice_init(i, j)* &
                               engine%state%icec%work%sw_down(i, j)
               a_dt = a_dt + DT*engine%state%ice_cgrid%tarea(i, j)*engine%state%ice%cice%aice_init(i, j)
               efrz = efrz + DT*engine%state%ice_cgrid%tarea(i, j)*max(engine%state%ice%frzmlt(i, j), 0.0_wp) &
                      - engine%state%ice_cgrid%tarea(i, j)*engine%state%ice%frazil_heat(i, j)
            end do
         end do
         run%efrz_indep = run%efrz_indep + efrz
         run%hb_delta = rho0_cp(engine)*(hb - run%hb0)
         run%stock_now = stock
         ! the ocean's clamp heat (its own accumulator) = the bank standing + what new ice spent (the ice's latent heat)
         r = max(abs(run%hb_delta), abs(run%frz_lat), 1.0_wp)
         run%clamp_mismatch = abs(run%hb_delta - run%stock_now + run%frz_lat)/r
       run%bl99_ratio = max(run%bl99_ratio, abs(L%tab%res(ICE_LT_HEAT, ICE_LT_THERM1))/(1.1_wp*ICE_FERRMAX*max(a_dt, tiny(1.0_wp))))
         if (abs(efrz) > 0.0_wp .or. abs(L%tab%frz_lat) > 0.0_wp) then
            run%frz_mismatch = max(run%frz_mismatch, abs(L%tab%frz_lat + efrz)/max(abs(efrz), abs(L%tab%frz_lat)))
            run%newice_state_mismatch = max(run%newice_state_mismatch, &
                                            abs(L%tab%dstate(3, ICE_LT_NEWICE) - L%tab%frz_lat)/max(abs(L%tab%frz_lat), &
                                                                                                    tiny(1.0_wp)))
         end if
      end associate
   end subroutine accumulate_tab

   subroutine run_regime(engine, cfg, nsteps, ok, run)
      type(ocean_engine_t), intent(inout) :: engine
      type(config_t), intent(inout) :: cfg
      integer, intent(in) :: nsteps
      logical, intent(out) :: ok
      type(stage_run_t), intent(out) :: run
      real(wp) :: t
      integer :: n

      t = 0.0_wp
      ok = .true.
      call init_clamp(engine, run)
      do n = 1, nsteps
         call one_step(engine, cfg, t, ok)
         if (.not. ok) return
         call pull_cice(engine)
         call pull_ice_flux(engine)
         if (any(iand(engine%state%ice%cice%ice_err, ICE_ERR_FATAL_MASK) /= 0_int32)) run%finite = .false.
         call accumulate_tab(engine, run)
      end do
   end subroutine run_regime

   subroutine report_regime(label, run)
      character(len=*), intent(in) :: label
      type(stage_run_t), intent(in) :: run
      integer :: s

      write (output_unit, "(a,a,i0,a)") "  ", label, run%nsteps, " steps: worst |res|/scale per stage (water, salt, heat)"
      do s = 1, ICE_LT_NSTAGE
         write (output_unit, "(4x,a10,3es11.3)") ice_lt_stage_name(s), run%rel(1, s), run%rel(2, s), run%rel(3, s)
      end do
      write (output_unit, "(a,es10.3)") "    heat in therm1 over the solver tolerance, worst per-step ratio (<= 1):", run%bl99_ratio
      write (output_unit, "(a,3es12.4)") "    sums of dstate (M, S, E) over stage newice:", run%dstate(:, ICE_LT_NEWICE)
      write (output_unit, "(a,es12.4,a,es12.4,a,es11.3)") "    new-ice energy (ledger):", run%frz_lat, "; bank the ocean paid:", &
         -run%efrz_indep, "; worst relative mismatch:", run%frz_mismatch
   end subroutine report_regime

   subroutine check_stage_closure(error, run, tol_mass, tol_heat)
      !! Every stage closes: mass and salt to `tol_mass`, heat to `tol_heat`, except heat in `therm1`, which closes
      !! to the vertical solver's own tolerance (`run%bl99_ratio <= 1`: `1.1 ICE_FERRMAX` per unit ice area).
      type(error_type), allocatable, intent(out) :: error
      type(stage_run_t), intent(in) :: run
      real(wp), intent(in) :: tol_mass, tol_heat
      integer :: s

      call check(error, run%finite, "finite ledger, no fatal ice error bit")
      if (allocated(error)) return
      do s = 1, ICE_LT_NSTAGE
         call check(error, run%rel(ICE_LT_WATER, s) < tol_mass, "water closes in stage "//trim(ice_lt_stage_name(s)))
         if (allocated(error)) return
         call check(error, run%rel(ICE_LT_SALT, s) < tol_mass, "salt closes in stage "//trim(ice_lt_stage_name(s)))
         if (allocated(error)) return
         if (s /= ICE_LT_THERM1) then
            call check(error, run%rel(ICE_LT_HEAT, s) < tol_heat, "heat closes in stage "//trim(ice_lt_stage_name(s)))
            if (allocated(error)) return
         end if
      end do
      call check(error, run%bl99_ratio <= 1.0_wp, "heat in therm1 is within the vertical solver's tolerance")
   end subroutine check_stage_closure

   subroutine test_column_growth(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      type(stage_run_t) :: run
      logical :: ok
      character(len=:), allocatable :: nml

      nml = replace(uni_nml("column"), "initial_temperature = -1.0", "initial_temperature = -1.84")
      nml = replace(nml, "T_init_surface = -1.0"//new_line('a')//"T_init_bottom = -1.0", "T_init_surface = -1.84"//new_line('a')//"T_init_bottom = -1.84")
      call make_engine(engine, cfg, nml, ok)
      call check(error, ok, "engine setup failed")
      if (allocated(error)) return
      call run_regime(engine, cfg, 72, ok, run)
      call check(error, ok, "72 steps ran")
      if (.not. allocated(error)) call report_regime("growth: ", run)
      if (.not. allocated(error)) call check_stage_closure(error, run, 1.0e-11_wp, 1.0e-11_wp)
      if (.not. allocated(error)) call check(error, run%dstate(ICE_LT_WATER, ICE_LT_THERM1) > 0.0_wp, "the ice grew in therm1")
      if (.not. allocated(error)) call check(error, run%atm_rad < 0.0_wp, "the ice lost heat to the atmosphere")
      if (.not. allocated(error)) call check(error, run%rain == 0.0_wp, "no rain through the ice (uniform atmosphere)")
      if (.not. allocated(error)) call check(error, run%fswabs_sum == 0.0_wp, "no shortwave forcing: nothing absorbed")
      call finish_engine(engine)
   end subroutine test_column_growth

   subroutine test_column_melt(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      type(stage_run_t) :: run
      logical :: ok
      character(len=:), allocatable :: nml

      nml = replace(uni_nml("column"), "uni_tas = 250.0"//new_line('a')//"uni_rlds = 150.0", "uni_tas = 275.0"//new_line('a')//"uni_rlds = 300.0")
      nml = replace(nml, "hfreeze = 10.0", "hfreeze = 10.0"//new_line('a')//"sw_down = 250.0")
      nml = replace(nml, "h_ice = 0.5", "h_ice = 0.8")
      call make_engine(engine, cfg, nml, ok)
      call check(error, ok, "engine setup failed")
      if (allocated(error)) return
      call run_regime(engine, cfg, 144, ok, run)
      call check(error, ok, "144 steps ran")
      if (.not. allocated(error)) call report_regime("melt: ", run)
      if (.not. allocated(error)) call check_stage_closure(error, run, 1.0e-11_wp, 1.0e-11_wp)
      if (.not. allocated(error)) call check(error, run%dstate(ICE_LT_WATER, ICE_LT_THERM1) < 0.0_wp, "the ice melted in therm1")
      if (.not. allocated(error)) call check(error, run%fout(ICE_LT_WATER, ICE_LT_THERM1) > 0.0_wp, "meltwater went to the ocean")
      ! the ice's absorbed shortwave against the forcing's downward shortwave and the albedo range (0.02 = 1 - the
      ! snow albedo 0.98, 0.64 = 1 - the bare-ice near-infrared albedo 0.36): independent of the ledger's own terms
      if (.not. allocated(error)) write (output_unit, "(a,es10.3)") "    absorbed / incident shortwave on the ice:", run%fswabs_sum/run%sw_in_sum
      if (.not. allocated(error)) call check(error, run%sw_in_sum > 0.0_wp .and. run%fswabs_sum >= 0.02_wp*run%sw_in_sum &
                                             .and. run%fswabs_sum <= 0.64_wp*run%sw_in_sum, &
                                             "fswabs is the forcing's shortwave times the absorbed fraction (albedo bounds)")
      call finish_engine(engine)
   end subroutine test_column_melt

   subroutine test_column_snow(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      type(stage_run_t) :: run
      logical :: ok
      character(len=:), allocatable :: nml

      nml = replace(uni_nml("column"), "hfreeze = 10.0", "hfreeze = 10.0"//new_line('a')//"snowfall = 2.0e-4")
      call make_engine(engine, cfg, nml, ok)
      call check(error, ok, "engine setup failed")
      if (allocated(error)) return
      call run_regime(engine, cfg, 72, ok, run)
      call check(error, ok, "72 steps ran")
      if (.not. allocated(error)) call report_regime("snow: ", run)
      if (.not. allocated(error)) call check_stage_closure(error, run, 1.0e-11_wp, 1.0e-11_wp)
      if (.not. allocated(error)) call check(error, run%snow > 0.0_wp, "snow fell on the ice")
      if (.not. allocated(error)) call check(error, run%snow_enth > 0.0_wp, "the falling snow carried (negative) enthalpy")
      call finish_engine(engine)
   end subroutine test_column_snow

   subroutine test_column_newice(error)
      !! An ice-free supercooled ocean: the frazil clamp banks heat, new ice forms.  Pins the SIGN of the new-ice
      !! energy: the stored energy FALLS by exactly the bank the ocean paid (read from the ocean's arrays).
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      type(stage_run_t) :: run
      logical :: ok
      character(len=:), allocatable :: nml

      nml = replace(uni_nml("column"), 'conc_config = "uniform"'//new_line('a')//"h_ice = 0.5"//new_line('a')//"conc = 1.0", &
                    'conc_config = "zero"'//new_line('a')//"h_ice = 0.5"//new_line('a')//"conc = 0.0")
      nml = replace(nml, "initial_temperature = -1.0", "initial_temperature = -2.0")
      nml = replace(nml, "T_init_surface = -1.0"//new_line('a')//"T_init_bottom = -1.0", "T_init_surface = -2.0"//new_line('a')//"T_init_bottom = -2.0")
      call make_engine(engine, cfg, nml, ok)
      call check(error, ok, "engine setup failed")
      if (allocated(error)) return
      call run_regime(engine, cfg, 96, ok, run)
      call check(error, ok, "96 steps ran")
      if (.not. allocated(error)) call report_regime("new ice in open water: ", run)
      if (.not. allocated(error)) call check_stage_closure(error, run, 1.0e-11_wp, 1.0e-11_wp)
      if (.not. allocated(error)) call check(error, run%frz_lat < 0.0_wp, "new ice SINKS the stored energy (the sign)")
      if (.not. allocated(error)) call check(error, run%efrz_indep > 0.0_wp, "the ocean paid a bank")
      if (.not. allocated(error)) call check(error, run%frz_mismatch < 1.0e-12_wp, &
                                             "the new-ice energy is minus the bank the ocean paid")
      if (.not. allocated(error)) write (output_unit, "(a,es10.3)") "    ocean clamp heat vs bank + new-ice latent, relative:", run%clamp_mismatch
      if (.not. allocated(error)) call check(error, run%clamp_mismatch < 1.0e-10_wp, &
                                             "the ocean's own frazil-clamp heat = the bank standing + the new ice's latent heat")
      if (.not. allocated(error)) call check(error, run%newice_state_mismatch < 1.0e-12_wp, &
                                             "and the state's energy change in the new-ice stage is that term")
      if (.not. allocated(error)) call check(error, run%dstate(ICE_LT_WATER, ICE_LT_NEWICE) > 0.0_wp, "new ice has mass")
      if (.not. allocated(error)) call check(error, run%fout(ICE_LT_WATER, ICE_LT_NEWICE) < 0.0_wp, &
                                             "and the ocean is charged for it (update_ocn_f)")
      call finish_engine(engine)
   end subroutine test_column_newice

   subroutine test_column_lateral(error)
      !! Ice over a warm ocean: the melt potential melts the base and the floe edges.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      type(stage_run_t) :: run
      logical :: ok
      character(len=:), allocatable :: nml

      nml = replace(uni_nml("column"), "conc = 1.0"//new_line('a')//"h_snow = 0.0", "conc = 0.6"//new_line('a')//"h_snow = 0.0")
      nml = replace(nml, "initial_temperature = -1.0", "initial_temperature = 0.5")
      nml = replace(nml, "T_init_surface = -1.0"//new_line('a')//"T_init_bottom = -1.0", "T_init_surface = 0.5"//new_line('a')//"T_init_bottom = 0.5")
      call make_engine(engine, cfg, nml, ok)
      call check(error, ok, "engine setup failed")
      if (allocated(error)) return
      call run_regime(engine, cfg, 72, ok, run)
      call check(error, ok, "72 steps ran")
      if (.not. allocated(error)) call report_regime("lateral + basal melt: ", run)
      if (.not. allocated(error)) call check_stage_closure(error, run, 1.0e-11_wp, 1.0e-11_wp)
      if (.not. allocated(error)) call check(error, run%absfout(ICE_LT_HEAT, ICE_LT_LATERAL) > 0.0_wp, &
                                             "the lateral melt booked heat to the ocean")
      if (.not. allocated(error)) call check(error, run%dstate(ICE_LT_WATER, ICE_LT_LATERAL) < 0.0_wp, &
                                             "and removed mass from the ice")
      call finish_engine(engine)
   end subroutine test_column_lateral

   ! ------------------------------------------------------------------
   ! B: the interface
   ! ------------------------------------------------------------------

   subroutine interface_case(error, nml, ratio, nsteps, label, prescribed)
      !! Window means of the ice's `fhocn`, `fswthru`, `fsalt`, `fresh` against the ocean's `heat_ice`, `sw_ice`,
      !! `salt_ice` at every hand-over, cell by cell, and the assembled `Q_heat`, `Q_salt` of the cells with no
      !! other source (uniform atmosphere).
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: nml, label
      integer, intent(in) :: ratio, nsteps
      logical, intent(in) :: prescribed
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      real(wp), allocatable :: a_h(:, :), a_sw(:, :), a_fs(:, :), a_fr(:, :)
      real(wp) :: t, e_h, e_sw, e_s, e_w, ws_h, ws_sw, ws_s, ws_w, exp_h, exp_sw, exp_s, exp_fr, scl
      integer :: n, i, j, ng, is, ie, js, je, nwin, nx, ny
      logical :: ok, thermo
      character(len=:), allocatable :: nm

      nm = nml
      if (prescribed) call write_synth_atm("mixed", FN_ATM, FN_LST)
      call make_engine(engine, cfg, nm, ok)
      call check(error, ok, "engine setup failed")
      if (allocated(error)) then
         if (prescribed) call remove_synth_atm(FN_ATM, FN_LST)
         return
      end if
      ng = engine%grid%nghost
      is = ng + 1
      ie = ng + engine%grid%nx_phys
      js = ng + 1
      je = ng + engine%grid%ny_phys
      nx = size(engine%state%metrics%areaT, 1)
      ny = size(engine%state%metrics%areaT, 2)
      allocate (a_h(nx, ny), a_sw(nx, ny), a_fs(nx, ny), a_fr(nx, ny))
      a_h = 0.0_wp
      a_sw = 0.0_wp
      a_fs = 0.0_wp
      a_fr = 0.0_wp
      ws_h = 0.0_wp
      ws_sw = 0.0_wp
      ws_s = 0.0_wp
      ws_w = 0.0_wp
      nwin = 0
      t = 0.0_wp
      do n = 1, nsteps
         call one_step(engine, cfg, t, ok, thermo)
         call check(error, ok, "step ran")
         if (allocated(error)) exit
         call pull_cice(engine)
         call pull_ice_flux(engine)
         call pull_sf(engine)
         associate (f => engine%state%icec%flux, w => engine%state%icec%work, sf => engine%state%surface_flux, &
                    c => engine%state%ice%cice)
            !$acc update self(f%fresh, f%fsalt, f%fhocn, f%fswthru)
            do j = js, je
               do i = is, ie
                  a_h(i, j) = a_h(i, j) + f%fhocn(i, j)
                  a_sw(i, j) = a_sw(i, j) + f%fswthru(i, j)
                  a_fs(i, j) = a_fs(i, j) + f%fsalt(i, j)
                  ! the share of the rain that falls through the ice is removed from the hand-over (prescribed atmosphere)
                  if (prescribed) then
                     a_fr(i, j) = a_fr(i, j) + f%fresh(i, j) - w%frain(i, j)*c%aice_init(i, j)
                  else
                     a_fr(i, j) = a_fr(i, j) + f%fresh(i, j)
                  end if
               end do
            end do
            if (thermo) then
               nwin = nwin + 1
               e_h = 0.0_wp
               e_sw = 0.0_wp
               e_s = 0.0_wp
               e_w = 0.0_wp
               exp_h = 0.0_wp
               exp_sw = 0.0_wp
               exp_s = 0.0_wp
               exp_fr = 0.0_wp
               do j = js, je
                  do i = is, ie
                     exp_h = a_h(i, j)/real(ratio, wp)
                     exp_sw = a_sw(i, j)/real(ratio, wp)
                     exp_s = 1000.0_wp*a_fs(i, j)/real(ratio, wp) - engine%state%icec%oin%sss_ocn(i, j)*a_fr(i, j)/real(ratio, wp)
                     exp_fr = a_fr(i, j)/real(ratio, wp)
                     scl = max(abs(exp_h), 1.0e-3_wp)
                     e_h = max(e_h, abs(sf%heat_ice(i, j) - exp_h)/scl)
                     scl = max(abs(exp_sw), 1.0e-3_wp)
                     e_sw = max(e_sw, abs(sf%sw_ice(i, j) - exp_sw)/scl)
                     scl = max(abs(exp_s), 1.0e-9_wp)
                     e_s = max(e_s, abs(sf%salt_ice(i, j) - exp_s)/scl)
                     if (.not. prescribed) then
                        ! nothing else feeds the assembler: Q_heat and Q_salt are the ice's
                        scl = max(abs(exp_h + exp_sw), 1.0e-3_wp)
                        e_w = max(e_w, abs(sf%Q_heat(i, j) - (sf%heat_ice(i, j) + sf%sw_ice(i, j)))/scl, &
                                  abs(sf%Q_salt(i, j) - sf%salt_ice(i, j))/max(abs(exp_s), 1.0e-9_wp))
                     end if
                  end do
               end do
               ws_h = max(ws_h, e_h)
               ws_sw = max(ws_sw, e_sw)
               ws_s = max(ws_s, e_s)
               ws_w = max(ws_w, e_w)
               a_h = 0.0_wp
               a_sw = 0.0_wp
               a_fs = 0.0_wp
               a_fr = 0.0_wp
            end if
         end associate
         if (any(iand(engine%state%ice%cice%ice_err, ICE_ERR_FATAL_MASK) /= 0_int32)) then
            call check(error, .false., "fatal ice error bit")
            exit
         end if
      end do
      if (.not. allocated(error)) then
   write (output_unit, "(a,a,i0,a,4es11.3)") "  ", label, nwin, " hand-overs, worst relative mismatch heat, sw, salt, assembled:", &
            ws_h, ws_sw, ws_s, ws_w
         call check(error, nwin >= 3, "enough hand-overs to compare")
         if (.not. allocated(error)) call check(error, ws_h < 1.0e-13_wp, "ice heat handed over = window mean of fhocn")
         if (.not. allocated(error)) call check(error, ws_sw < 1.0e-13_wp, "ice shortwave handed over = window mean of fswthru")
         if (.not. allocated(error)) call check(error, ws_s < 1.0e-13_wp, &
                                                "ice salt handed over = 1000 fsalt - S_surf fresh (window means)")
         if (.not. allocated(error) .and. .not. prescribed) &
            call check(error, ws_w < 1.0e-13_wp, "the assembled Q_heat, Q_salt are the ice's (no other source)")
      end if
      call finish_engine(engine)
      if (prescribed) call remove_synth_atm(FN_ATM, FN_LST)
   end subroutine interface_case

   subroutine test_interface_uniform(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: nml

      nml = replace(uni_nml("basin"), "dt_therm_ratio = 1", "dt_therm_ratio = 2")
      nml = replace(nml, "uni_tas = 255.0"//new_line('a')//"uni_rlds = 180.0", "uni_tas = 268.0"//new_line('a')//"uni_rlds = 230.0")
      nml = replace(nml, "hfreeze = 10.0", "hfreeze = 10.0"//new_line('a')//"sw_down = 120.0")
      call interface_case(error, nml, 2, 40, "uniform atmosphere, ratio 2: ", .false.)
   end subroutine test_interface_uniform

   subroutine test_interface_prescribed(error)
      type(error_type), allocatable, intent(out) :: error

      call interface_case(error, pres_nml(3), 3, 42, "prescribed atmosphere, ratio 3: ", .true.)
   end subroutine test_interface_prescribed

   ! ------------------------------------------------------------------
   ! C: the coupled global closure
   ! ------------------------------------------------------------------

   subroutine global_case(error, nml, prescribed, nsteps, label, mass, what)
      !! ocean + ice on a closed basin.  The totals come from the STATE (ocean tracers, the ice's EFP sums); the
      !! flows from the flux arrays; the open-water atmosphere from the BULK-FORCING outputs (`q_sw`, `q_lw`,
      !! `q_lat`, `q_sens`, `q_fus` (latent heat of melting the open-water snow), the precipitation enthalpies, `cp T evap`), not from the assembler; the frazil clamp
      !! from the ocean's own accumulator.  The ledger enters only through its SOURCE terms (the atmosphere
      !! exchange of the ice and the new-ice latent heat); its residual `R` is neither subtracted nor used: the
      !! vertical solver's tolerance is a BOUND on the closure.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: nml, label
      logical, intent(in) :: prescribed
      integer, intent(in) :: nsteps
      logical, intent(in), optional :: mass
         !! Real fresh-water mass (`freshwater = 'mass'` with the ice loading): the water, salt and heat closures
         !! of the REAL-MASS hand-off (the ocean's volume moves) replace the virtual-salt ones.
      character(len=*), intent(in), optional :: what
         !! With `mass`: which closure this call asserts (`water`, `salt`, `heat`, `all`); every term is always
         !! computed and reported.
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      type(stage_run_t) :: run
      logical :: is_mass
      character(len=8) :: which
      real(wp) :: moc0, moc1, oc_atm_n, oc_atm_pend, oc_atm_app, fresh_pend, fresh_app, hmc_n, hmc_pend, hmc_app
      real(wp) :: wdt, r_oc_now, r_oc_worst, melt_n, salt_n, salt_pend, salt_app_all, salt_book_all
      logical :: hand, deliver
      real(wp) :: r_oc, sc_oc, r_cpl, sc_cpl, r_heat_m, hmc_conv_closure, melt_mis, fwnet_mis, f_last, rho0
      real(wp) :: repair_ice, repair_snow, di, ds, repair_sum, t, h0, so0, h1, so1, m0, si0, e0, m1, si1, e1, cp_w, sst
      real(wp) :: open_forc, open_asm, open_forc_app, open_forc_pend, ice_asm, ice_app, ice_pend, ice_asm_last
      real(wp) :: icefl, icefl_all, q_n, qs_n, fsalt_n, fresh_n, fresh_eff_n, dil_n, dil_app, dil_pend
      real(wp) :: fsalt_last, fsalt_app_all, fresh_all, fresh_eff_sum, a_ice_dt, sum_ext_h, w_in
      real(wp) :: expected, r_heat, sc_heat, bound_bl99, r_salt, sc_salt, r_water, sc_water
      real(wp) :: clamp_ocean, stock, dil_ocean, h_mass_conv, tsum, asum, gap_rain_enth, gap_open_water
      real(wp) :: precip_nominal, rain_ice_sum, flux_rain_equiv, wmis, hmis, vmax, ridge_max, hand_mis
      character(len=200) :: line
      integer :: n, i, j, ng, is, ie, js, je, nz, nbits
      logical :: ok
      real(wp), allocatable :: area(:, :)

      nbits = 0
      is_mass = .false.
      if (present(mass)) is_mass = mass
      which = "all"
      if (present(what)) which = what
      moc0 = 0.0_wp
      oc_atm_pend = 0.0_wp
      oc_atm_app = 0.0_wp
      fresh_pend = 0.0_wp
      fresh_app = 0.0_wp
      hmc_pend = 0.0_wp
      hmc_app = 0.0_wp
      melt_mis = 0.0_wp
      fwnet_mis = 0.0_wp
      if (prescribed) call write_synth_atm("mixed", FN_ATM, FN_LST)
      call make_engine(engine, cfg, nml, ok)
      call check(error, ok, "engine setup failed")
      if (allocated(error)) then
         if (prescribed) call remove_synth_atm(FN_ATM, FN_LST)
         return
      end if
      ng = engine%grid%nghost
      is = ng + 1
      ie = ng + engine%grid%nx_phys
      js = ng + 1
      je = ng + engine%grid%ny_phys
      nz = size(engine%state%multilayer%h_layer, 3)
      allocate (area(size(engine%state%metrics%areaT, 1), size(engine%state%metrics%areaT, 2)))
      area = engine%state%metrics%areaT
      cp_w = engine%state%eos%c_p
      call ocean_totals(engine, h0, so0)
      rho0 = engine%state%eos%rho0
      deliver = .false.
      r_oc_worst = 0.0_wp
      salt_pend = 0.0_wp
      salt_app_all = 0.0_wp
      salt_book_all = 0.0_wp
      moc0 = ocean_volume_mass(engine)
      call engine%state%icec%ledger%totals(engine%state%ice%cice, m0, si0, e0)
      run = stage_run_t()
      call init_clamp(engine, run)
      t = 0.0_wp
      sum_ext_h = 0.0_wp
      open_forc_app = 0.0_wp
      open_forc_pend = 0.0_wp
      ice_app = 0.0_wp
      ice_pend = 0.0_wp
      icefl_all = 0.0_wp
      dil_app = 0.0_wp
      dil_pend = 0.0_wp
      fsalt_app_all = 0.0_wp
      fsalt_last = 0.0_wp
      fresh_all = 0.0_wp
      fresh_eff_sum = 0.0_wp
      w_in = 0.0_wp
      a_ice_dt = 0.0_wp
      ice_asm_last = 0.0_wp
      gap_rain_enth = 0.0_wp
      gap_open_water = 0.0_wp
      precip_nominal = 0.0_wp
      rain_ice_sum = 0.0_wp
      wmis = 0.0_wp
      hmis = 0.0_wp
      vmax = 0.0_wp
      ridge_max = 0.0_wp
      repair_sum = 0.0_wp
      repair_ice = 0.0_wp
      repair_snow = 0.0_wp
      do n = 1, nsteps
         call pull_cat_state(engine)
         repair_sum = repair_sum + profile_repair(engine, di, ds)
         repair_ice = repair_ice + di
         repair_snow = repair_snow + ds
         call one_step(engine, cfg, t, ok, hand)
         call check(error, ok, "step ran")
         if (allocated(error)) exit
         ! `hand`: this step's ice step handed the ocean its window (the thermo cadence, counter advanced); the
         ! ocean delivers it in the NEXT step (`deliver`).  `wdt` is the window those hand-over fluxes cover and is
         ! zero between hand-offs, so the sums below count each hand-over once whatever `dt_therm_ratio` is.
         wdt = merge(engine%state%dyn%therm_dt(DT), 0.0_wp, hand)
         call pull_cice(engine)
         call pull_ice_flux(engine)
         call pull_sf(engine)
         !$acc update self(engine%state%multilayer%fw_dhT)
         !$acc update self(engine%state%multilayer%tracers(engine%state%multilayer%idx_temperature)%c)
         !$acc update self(engine%state%ice%u_ice, engine%state%ice%v_ice)
         call accumulate_tab(engine, run)
         !$acc update self(engine%state%icec%work%ierr_n)
         nbits = nbits + count(iand(engine%state%icec%work%ierr_n(is:ie, js:je, :), ICE_ERR_THERMO_CONSERV) /= 0_int32)
         vmax = max(vmax, maxval(abs(engine%state%ice%u_ice)))
         ridge_max = max(ridge_max, maxval(engine%state%icec%flux%dardg1dt(is:ie, js:je)))
         if (any(iand(engine%state%ice%cice%ice_err, ICE_ERR_FATAL_MASK) /= 0_int32)) then
            call check(error, .false., "fatal ice error bit")
            exit
         end if
         associate (c => engine%state%ice%cice, f => engine%state%icec%flux, wk => engine%state%icec%work, &
                    sf => engine%state%surface_flux, ms => engine%state%multilayer, L => engine%state%icec%ledger)
            icefl = 0.0_wp
            ice_asm = 0.0_wp
            q_n = 0.0_wp
            qs_n = 0.0_wp
            open_forc = 0.0_wp
            fsalt_n = 0.0_wp
            fresh_n = 0.0_wp
            fresh_eff_n = 0.0_wp
            dil_n = 0.0_wp
            oc_atm_n = 0.0_wp
            hmc_n = 0.0_wp
            melt_n = 0.0_wp
            salt_n = 0.0_wp
            do j = js, je
               do i = is, ie
                  sst = ms%tracers(ms%idx_temperature)%c(i, j, nz)
                  icefl = icefl + DT*area(i, j)*(f%fhocn(i, j) + f%fswthru(i, j))
                  ice_asm = ice_asm + wdt*area(i, j)*(sf%heat_ice(i, j) + sf%sw_ice(i, j))
                  q_n = q_n + wdt*area(i, j)*sf%Q_heat(i, j)
                  qs_n = qs_n + wdt*area(i, j)*sf%Q_salt(i, j)
                  ! the open water's heat from the bulk-forcing outputs (the assembler's formula, rebuilt from them)
                  if (is_mass) then
                     ! real mass: the enthalpy of the open water's mass fluxes is the layer-temperature content the
                     ! mass step books (below), never a `Q_heat` term
                     open_forc = open_forc + wdt*area(i, j)*(sf%q_sw(i, j) + sf%q_lw(i, j) + sf%q_lat(i, j) &
                                                             + sf%q_sens(i, j) + sf%q_fus(i, j))
                     oc_atm_n = oc_atm_n + wdt*area(i, j)*(sf%lprec(i, j) + sf%fprec(i, j) + sf%evap(i, j))
                     ! what the ocean was handed of the ice: the assembled window means (`seaice_melt`, `salt_ice`)
                     melt_n = melt_n + wdt*area(i, j)*sf%seaice_melt(i, j)
                     salt_n = salt_n + wdt*area(i, j)*sf%salt_ice(i, j)/1000.0_wp
                     if (deliver) hmc_n = hmc_n + rho0_cp(engine)*area(i, j)*ms%fw_dhT(i, j)
            if (hand .and. engine%state%dyn%dt_therm_ratio <= 1) melt_mis = max(melt_mis, abs(sf%seaice_melt(i, j) - f%fresh(i, j)))
                     if (hand) fwnet_mis = max(fwnet_mis, abs(sf%fw_net(i, j) - (sf%lprec(i, j) + sf%fprec(i, j) + sf%evap(i, j) &
                                                                                 + sf%seaice_melt(i, j) + sf%vprec(i, j) &
                                                                                 + sf%lrunoff(i, j) + sf%frunoff(i, j))))
                  else
                     open_forc = open_forc + wdt*area(i, j)*(sf%q_sw(i, j) + sf%q_lw(i, j) + sf%q_lat(i, j) &
                                                             + sf%q_sens(i, j) + sf%q_fus(i, j) &
                                                             + sf%heat_content_lprec(i, j) &
                                                             + sf%heat_content_fprec(i, j) + cp_w*sst*sf%evap(i, j))
                  end if
                  fsalt_n = fsalt_n + DT*area(i, j)*f%fsalt(i, j)
                  fresh_n = fresh_n + DT*area(i, j)*f%fresh(i, j)
                  a_ice_dt = a_ice_dt + DT*area(i, j)*c%aice_init(i, j)
                  if (prescribed) then
                     ! the rain that falls through the ice is not handed to the ocean (it only passes the ice's ledger)
                     fresh_eff_n = fresh_eff_n + DT*area(i, j)*(f%fresh(i, j) - wk%frain(i, j)*c%aice_init(i, j))
                     dil_n = dil_n + DT*area(i, j)*engine%state%icec%oin%sss_ocn(i, j)* &
                             (f%fresh(i, j) - wk%frain(i, j)*c%aice_init(i, j))
                     rain_ice_sum = rain_ice_sum + DT*area(i, j)*wk%frain(i, j)*c%aice_init(i, j)
                     ! the precipitation the open water receives is the forcing's times the open fraction
                     wmis = max(wmis, abs(sf%lprec(i, j) - wk%frain(i, j)*(1.0_wp - c%aice_init(i, j))) &
                                /max(abs(wk%frain(i, j)), 1.0e-30_wp), &
                                abs(sf%fprec(i, j) - wk%fsnow(i, j)*(1.0_wp - c%aice_init(i, j))) &
                                /max(abs(wk%fsnow(i, j)), 1.0e-30_wp))
                     ! the enthalpy of the rain on the ice-covered fraction, from its own formula, never delivered
                     gap_rain_enth = gap_rain_enth + cp_w*sst*DT*area(i, j)*wk%frain(i, j)*c%aice_init(i, j)
                     gap_open_water = gap_open_water + DT*area(i, j)*(1.0_wp - c%aice_init(i, j))* &
                                      (wk%fsnow(i, j) + wk%frain(i, j))
                     precip_nominal = precip_nominal + DT*area(i, j)*(wk%fsnow(i, j) + wk%frain(i, j))
                     if (abs(sf%lprec(i, j)) > 0.0_wp) &
                        hmis = max(hmis, abs(sf%heat_content_lprec(i, j) - cp_w*sst*sf%lprec(i, j))/abs(cp_w*sst*sf%lprec(i, j)))
                  else
                     fresh_eff_n = fresh_eff_n + DT*area(i, j)*f%fresh(i, j)
                     dil_n = dil_n + DT*area(i, j)*engine%state%icec%oin%sss_ocn(i, j)*f%fresh(i, j)
                     gap_open_water = gap_open_water + DT*area(i, j)*(1.0_wp - c%aice_init(i, j))*wk%fsnow(i, j)
                     precip_nominal = precip_nominal + DT*area(i, j)*wk%fsnow(i, j)
                  end if
               end do
            end do
            ! the assembler's open-water share against the rebuilt one, this step
            open_asm = q_n - ice_asm
            if (hand .and. (abs(open_forc) > 0.0_wp .or. abs(open_asm) > 0.0_wp)) &
               run%asm_mismatch = max(run%asm_mismatch, abs(open_asm - open_forc)/max(abs(open_forc), abs(open_asm)))
            ! the ocean applies the flux assembled by the PREVIOUS step
            oc_atm_app = oc_atm_app + oc_atm_pend
            oc_atm_pend = oc_atm_n
            fresh_app = fresh_app + fresh_pend
            fresh_pend = merge(melt_n, fresh_n, is_mass)
            salt_app_all = salt_app_all + salt_pend
            salt_pend = salt_n
            salt_book_all = salt_book_all + fsalt_n
            hmc_app = hmc_app + hmc_n
            deliver = hand
            if (is_mass) then
               r_oc_now = (ocean_volume_mass(engine) - moc0) - (fresh_app + oc_atm_app)
               r_oc_worst = max(r_oc_worst, abs(r_oc_now))
            end if
            open_forc_app = open_forc_app + open_forc_pend
            ice_app = ice_app + ice_pend
            dil_app = dil_app + dil_pend
            fsalt_app_all = fsalt_app_all + fsalt_last
            open_forc_pend = open_forc
            ice_pend = ice_asm
            if (is_mass) dil_n = 0.0_wp
            dil_pend = dil_n
            fsalt_last = fsalt_n
            ice_asm_last = ice_asm
            icefl_all = icefl_all + icefl
            fresh_all = fresh_all + fresh_n
            fresh_eff_sum = fresh_eff_sum + fresh_eff_n
            sum_ext_h = sum_ext_h + L%tab%lw_in - L%tab%snow_enth + L%tab%atm_rad + L%tab%atm_lat
            w_in = w_in + L%tab%snow + L%tab%evap
         end associate
      end do
      if (.not. allocated(error)) then
         call ocean_totals(engine, h1, so1)
         call engine%state%icec%ledger%totals(engine%state%ice%cice, m1, si1, e1)
         ! real mass: the salt in transit = booked by the ice (every step) - the windows the ocean has been handed and applied
         if (is_mass) fsalt_last = salt_book_all - salt_app_all
         clamp_ocean = run%hb_delta
         stock = run%stock_now
         ! ---- heat: the closure.  dH_ocean + dE_ice = what the ice received from the atmosphere + what the open
         !      water received (applied windows) + the ice hand-over the ocean has applied - what the ice handed
         !      over (all steps) + the ice's new-ice latent heat + the heat the frazil clamp added to the ocean.
         expected = sum_ext_h + open_forc_app + ice_app - icefl_all + run%frz_lat + clamp_ocean
         r_heat = (h1 - h0) + (e1 - e0) - expected
         bound_bl99 = 1.1_wp*ICE_FERRMAX*a_ice_dt
         sc_heat = max(abs(h1 - h0), abs(e1 - e0), abs(sum_ext_h), abs(open_forc_app), abs(icefl_all), 1.0_wp)
         hand_mis = abs(ice_app + ice_asm_last - icefl_all)/max(abs(icefl_all), 1.0_wp)
         ! ---- the heat the model does not carry, from its own formulas (reported, see the module header)
         tsum = 0.0_wp
         asum = 0.0_wp
         do j = js, je
            do i = is, ie
               tsum = tsum + area(i, j)*engine%state%multilayer%tracers(engine%state%multilayer%idx_temperature)%c(i, j, nz)
               asum = asum + area(i, j)
            end do
         end do
         ! water that leaves the ocean at its own temperature carries cp T dM against the 0 degC reference; the
         ! ocean volume never changes (virtual salt), so it is booked nowhere: zero by configuration
         h_mass_conv = cp_w*(tsum/asum)*(-(m1 - m0))
         moc1 = ocean_volume_mass(engine)
         f_last = fresh_pend
         ! ---- real mass.  dM_ocean = what the ice (fresh, window mean) and the open-water atmosphere (lprec + fprec +
         !      evap) handed over and the ocean has APPLIED (the previous step's assembly), read from the ice flux and the
         !      bulk-forcing arrays, never from `fw_net`.  Coupled: ocean + ice + snow together change only by the
         !      atmosphere's exchange with ALL fractions (the open water's, the ice's snow/sublimation, the rain
         !      through the ice) minus what the ice booked and the ocean has not yet applied.
         r_oc = (moc1 - moc0) - (fresh_app + oc_atm_app)
         sc_oc = max(abs(moc1 - moc0), abs(fresh_app), abs(oc_atm_app), 1.0e-30_wp)
         r_cpl = (moc1 - moc0) + (m1 - m0) - (oc_atm_app + w_in + run%rain - f_last)
         sc_cpl = max(abs(moc1 - moc0), abs(m1 - m0), abs(fresh_all), 1.0e-30_wp)
         ! heat: the same closure with the mass flux's enthalpy as the layer-temperature content the mass step books
         ! (`cp T dM`, ice reference 0 degC liquid).  Without it the residual IS that named term.
         expected = sum_ext_h + open_forc_app + ice_app - icefl_all + run%frz_lat + clamp_ocean + hmc_app
         r_heat_m = (h1 - h0) + (e1 - e0) - expected
         hmc_conv_closure = r_heat_m + hmc_app
         ! ---- salt (kg): 1e-3 dSo + dSi + fsalt in transit + the applied virtual-salt dilution = 0
         r_salt = (si1 - si0) + (so1 - so0)/1000.0_wp + fsalt_last + dil_app/1000.0_wp
         sc_salt = max(abs(si1 - si0), abs(so1 - so0)/1000.0_wp, abs(dil_app)/1000.0_wp, 1.0e-30_wp)
         ! the water the ocean was handed, read from its S/h arrays: 1000 sum(fsalt) - d(salt content) is the dilution
         dil_ocean = 1000.0_wp*fsalt_app_all - (so1 - so0)
         ! ---- water.  Under the virtual-salt hand-off the ocean has NO water budget (its volume never changes); the
         !      water the ice gives it is felt only as the dilution above.  What is conserved: ice + snow mass change =
         !      the atmosphere's mass exchange with the ice (snow, condensation, rain through the ice) minus what the
         !      ice booked as `fresh` (state totals against the flux arrays), and the ocean's receipt of that
         !      `fresh` is the dilution it shows in its S/h arrays (salt line).
         r_water = (m1 - m0) + fresh_all - w_in - run%rain
         sc_water = max(abs(m1 - m0), abs(fresh_all), abs(w_in), 1.0e-30_wp)
         flux_rain_equiv = abs(gap_rain_enth)/(max(real(nsteps, wp)*DT, 1.0_wp)*sum(area(is:ie, js:je)))
         write (output_unit, "(a,a)") "  ", label
         write (output_unit, "(a,6es12.4)") "    heat terms: dH_ocean, dE_ice, atm->ice, atm->open (applied), ice hand-over applied, handed over:", &
            h1 - h0, e1 - e0, sum_ext_h, open_forc_app, ice_app, icefl_all
      write (output_unit, "(a,3es12.4)") "    heat terms: new-ice latent (ice), frazil clamp (ocean accumulator), bank standing:", &
            run%frz_lat, clamp_ocean, stock
     write (output_unit, "(a,3es12.4)") "    heat: closure residual, its bound (BL99), ledger R (info only):", r_heat, bound_bl99, &
            run%rsum(ICE_LT_HEAT)
         write (output_unit, "(a,8es10.2)") "    ledger heat residual per stage (cumulative):", engine%state%icec%ledger%cum_res(ICE_LT_HEAT, :)
         write (output_unit, "(a,2es12.4)") "    therm1 heat residual vs the thin-snow enthalpy reset (energy added before the step, from the pulled state):", &
            engine%state%icec%ledger%cum_res(ICE_LT_HEAT, ICE_LT_THERM1), repair_sum
         write (output_unit, "(a,2es12.4)") "    of which ice layers, snow:", repair_ice, repair_snow
         write (output_unit, "(a,es10.3,a,i0)") "    worst per-step therm1 heat residual over its bound:", run%bl99_ratio, &
            "; categories flagged by the vertical thermodynamics' own energy check:", nbits
         write (output_unit, "(a,es10.3,a,es10.3,a,es10.3)") "    heat: |residual - snow reset|/bound", abs(r_heat - repair_sum)/bound_bl99, &
            "; |residual|/scale", abs(r_heat)/sc_heat, "; hand-over mismatch", hand_mis
      write (output_unit, "(a,es12.4,a,es10.3)") "    H_mass_conv (cp T dM, booked nowhere, zero by configuration):", h_mass_conv, &
            " J, relative to dH", abs(h_mass_conv)/max(abs(h1 - h0), 1.0_wp)
  write (output_unit, "(a,3es12.4,a,es10.3)") "    salt: dSi, dSo/1000, applied dilution/1000:", si1 - si0, (so1 - so0)/1000.0_wp, &
            dil_app/1000.0_wp, "; relative residual", abs(r_salt)/sc_salt
 write (output_unit, "(a,2es12.4)") "    water the ocean was handed (dilution), from its S arrays vs from the ice flux (kg psu):", &
            dil_ocean, dil_app
         write (output_unit, "(a,3es12.4,a,es10.3)") "    water: dM, booked fresh, sources:", m1 - m0, fresh_all, w_in, &
            "; relative residual", abs(r_water)/sc_water
         write (output_unit, "(a,es12.4,1x,es10.3,a)") "    diagnostic, not delivered: rain enthalpy on the ice-covered fraction", &
            gap_rain_enth, flux_rain_equiv, " W/m2 equivalent"
         write (output_unit, "(a,2es12.4)") "    diagnostic, not delivered: precipitation on the open fraction, nominal (kg):", &
            gap_open_water, precip_nominal
   write (output_unit, "(a,3es10.3)") "    assembler vs rebuilt open-water heat, precip weights, rain enthalpy (worst relative):", &
            run%asm_mismatch, wmis, hmis
         write (output_unit, "(a,es10.3,a,es10.3)") "    ice speed max", vmax, ", ridging rate max", ridge_max
         line = engine%state%icec%ledger%console_line()
         write (output_unit, "(a,a)") "    console: ", trim(line)
         call ice_cice_ledger_report(engine%state, 0)
         call check(error, index(line, "Ice ledger: ") == 1, "the console line is the ledger's")
         if (is_mass) then
            write (output_unit, "(a)") "    --- real fresh-water mass ---"
            write (output_unit, "(a,3es12.4,a,es10.3)") "    ocean: dM, applied ice fresh, applied atmosphere (kg):", moc1 - moc0, &
               fresh_app, oc_atm_app, "; relative residual", abs(r_oc)/sc_oc
  write (output_unit, "(a,4es12.4,a,es10.3)") "    coupled: dM_ocean, dM_ice+snow, atmosphere (all fractions), undelivered (kg):", &
               moc1 - moc0, m1 - m0, oc_atm_app + w_in + run%rain, f_last, "; relative residual", abs(r_cpl)/sc_cpl
            write (output_unit, "(a,es12.4,a,es12.4)") "    NAMED ledger term: ice salt mass (kg) dSi =", si1 - si0, &
               "; its share of the ice mass change", abs(si1 - si0)/max(abs(m1 - m0), 1.0e-30_wp)
        write (output_unit, "(a,es12.4,a,es10.3)") "    salt (ocean + ice, kg): residual", r_salt, "; relative", abs(r_salt)/sc_salt
            write (output_unit, "(a,es12.4,a,es12.4,a,es10.3)") "    NAMED ledger term H_mass_conv (cp T dM, J) =", hmc_app, &
               "; closure residual without it", hmc_conv_closure, "; with it (BL99 bound)", bound_bl99
            write (output_unit, "(a,es12.4,a,es12.4)") "    heat: residual with the term", r_heat_m, "; snow reset", repair_sum
    write (output_unit, "(a,es10.3,a,es10.3)") "    heat: |residual - snow reset| / bound", abs(r_heat_m - repair_sum)/bound_bl99, &
               "; worst per-step |ocean mass - delivered|/column mass", r_oc_worst/moc0
            write (output_unit, "(a,2es10.3)") "    seaice_melt vs fresh, fw_net vs the component sum (worst abs, kg m-2 s-1):", melt_mis, fwnet_mis
            if (.not. allocated(error)) call check(error, run%finite, "finite")
  if (.not. allocated(error)) call check(error, melt_mis <= 1.0e-12_wp*max(maxval(abs(engine%state%icec%flux%fresh)), 1.0e-30_wp), &
                                                   "the adapter hands the ice's fresh over as seaice_melt")
        if (.not. allocated(error)) call check(error, fwnet_mis <= 1.0e-12_wp*max(maxval(abs(sf_fw_net_ref(engine))), 1.0e-30_wp), &
                                                   "fw_net is the sum of the mass components")
            if (.not. allocated(error)) call check(error, abs(r_water)/sc_water < 1.0e-13_wp, &
                                                   "ice mass change = the atmosphere's exchange - the booked fresh")
            if (.not. allocated(error) .and. (which == "water" .or. which == "all")) then
               call check(error, abs(r_oc) < 1.0e-13_wp*moc0, "the ocean's column mass moved by exactly what the ice and the atmosphere handed it (to round-off of the column mass)")
               if (.not. allocated(error)) call check(error, r_oc_worst < 1.0e-13_wp*moc0, &
       "the ocean's column mass followed the delivered flux arrays at EVERY step, the first interval after the cold start included")
               if (.not. allocated(error) .and. engine%state%dyn%dt_therm_ratio <= 1) &
                  call check(error, abs(r_cpl) < 1.0e-13_wp*moc0, &
                             "ocean + ice + snow mass = the atmosphere over all fractions (to round-off of the column mass)")
               if (.not. allocated(error)) call check(error, abs(si1 - si0) > 0.0_wp, "the ice salt mass is a live named term")
            end if
            if (.not. allocated(error) .and. (which == "salt" .or. which == "all")) then
               call check(error, abs(r_salt) < 1.0e-13_wp*abs(so0)/1000.0_wp, "ocean salt + ice salt is exact: no dilution term, no volume for the ice's salt")
            end if
            if (.not. allocated(error) .and. (which == "heat" .or. which == "all")) then
               call check(error, abs(r_heat_m - repair_sum) <= bound_bl99 + 1.0e-11_wp*sc_heat, &
        "ocean + ice heat closes with the named H_mass_conv term at the temperatures the mass step used, within the solver's bound")
               if (.not. allocated(error)) call check(error, abs(hmc_app) > 0.0_wp, "H_mass_conv is a live named term")
               if (.not. allocated(error)) call check(error, run%asm_mismatch < 1.0e-12_wp, &
                                                      "the assembler's open-water heat is the bulk forcing's non-mass part")
            end if
            call finish_engine(engine)
            if (prescribed) call remove_synth_atm(FN_ATM, FN_LST)
            return
         end if
         if (.not. allocated(error)) call check(error, run%finite, "finite")
         if (.not. allocated(error)) call check(error, abs(r_salt)/sc_salt < 1.0e-10_wp, &
                                                "ocean salt + ice salt + virtual-salt dilution: closed")
         if (.not. allocated(error)) call check(error, abs(dil_ocean - dil_app)/max(abs(dil_app), 1.0_wp) < 1.0e-10_wp, &
                                                "the water the ocean was handed, read from its own S arrays, is the ice's fresh")
         if (.not. allocated(error)) call check(error, abs(r_water)/sc_water < 1.0e-13_wp, &
                                                "ice mass change = the atmosphere's exchange - the booked fresh")
         if (.not. allocated(error)) call check(error, hand_mis < 1.0e-12_wp, &
                                                "the ice hand-over: applied + in transit == what the ice booked")
         if (.not. allocated(error)) call check(error, abs(r_heat - repair_sum) <= bound_bl99 + 1.0e-11_wp*sc_heat, &
                                          "ocean + ice heat closes, up to the thin-snow reset, within the solver's tolerance bound")
         if (.not. allocated(error)) call check(error, run%asm_mismatch < 1.0e-12_wp, &
                                                "the assembler's open-water heat is the bulk forcing's")
         if (.not. allocated(error)) call check(error, run%clamp_mismatch < 1.0e-10_wp, &
                                                "the ocean's frazil clamp heat = the bank standing + the new-ice latent heat")
         if (.not. allocated(error)) call check(error, run%frz_mismatch < 1.0e-12_wp, "new-ice energy = minus the bank paid")
         if (.not. allocated(error)) call check(error, vmax > 1.0e-6_wp, "the ice moves")
         if (.not. allocated(error) .and. prescribed) then
            call check(error, wmis < 1.0e-12_wp, "the open water receives the forcing's precipitation times aice0")
            if (.not. allocated(error)) call check(error, hmis < 5.0e-2_wp, "and its enthalpy cp T lprec (the SST of the bulk forcing is one step older: ~1 %)")
         end if
      end if
      call finish_engine(engine)
      if (prescribed) call remove_synth_atm(FN_ATM, FN_LST)
   end subroutine global_case

   subroutine test_global_uniform(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: nml

      ! cold air over a warmish ocean: surface growth, basal and lateral melt, wind-driven ridging in a closed basin
      nml = replace(uni_nml("basin"), "initial_temperature = -1.0", "initial_temperature = -0.6")
      nml = replace(nml, "T_init_surface = -1.0"//new_line('a')//"T_init_bottom = -1.0", "T_init_surface = -0.6"//new_line('a')//"T_init_bottom = -0.6")
      nml = replace(nml, "hfreeze = 10.0", "hfreeze = 10.0"//new_line('a')//"snowfall = 3.0e-5")
      call global_case(error, nml, .false., 144, "uniform atmosphere, closed basin, wind: ")
   end subroutine test_global_uniform

   subroutine test_global_prescribed(error)
      type(error_type), allocatable, intent(out) :: error

      call global_case(error, pres_nml(1), .true., 96, "prescribed mixed atmosphere, closed basin: ")
   end subroutine test_global_prescribed

   subroutine test_global_prescribed_frazil(error)
      !! The same basin with an ocean at the freezing point: the cold half makes frazil, so the ocean's clamp, the new-ice
      !! latent heat and the bank are all active.  Here the heat in `therm1` exceeds the per-ice-area solver bound (1.2 of it
      !! cumulatively, 3.1 in the worst step) although the vertical thermodynamics' own per-category check flags nothing:
      !! 96 % of it (1.78e12 of 1.85e12 J) is the profile initialisation resetting snow thinner than `hs_min` to the melting
      !! enthalpy `-rhos Lfresh` BEFORE the check's `einit` is formed (upstream behaviour, ported verbatim), an energy
      !! creation outside the check.  The test computes that term independently (the ice-core routine on the pulled
      !! start-of-step state) and bounds the rest by the solver tolerance.
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: nml

      nml = replace(pres_nml(1), "initial_temperature = -1.0", "initial_temperature = -1.84")
      nml = replace(nml, "T_init_surface = -1.0"//new_line('a')//"T_init_bottom = -1.0", "T_init_surface = -1.84"//new_line('a')//"T_init_bottom = -1.84")
      call global_case(error, nml, .true., 96, "prescribed mixed atmosphere, supercooled ocean (frazil): ")
   end subroutine test_global_prescribed_frazil

   ! ------------------------------------------------------------------
   ! G: real fresh-water mass (`[ocean.forcing] freshwater = "mass"`) with the ice loading
   ! ------------------------------------------------------------------

   function mass_nml(base) result(nml)
      !! `base` with the real-mass hand-off and the ice weight on the surface pressure.
      character(len=*), intent(in) :: base
      character(len=:), allocatable :: nml
      character(len=*), parameter :: NL = new_line("a")

      nml = replace(base, "[ocean.forcing]"//new_line('a')// &
                    "enable_components = true"//new_line('a'), &
                    "[ocean.forcing]"//new_line('a')// &
                    "enable_components = true"//new_line('a')// &
                    'freshwater = "mass"'//new_line('a')// &
                    "[ocean.psurf]"//new_line('a')// &
                    "enable = true"//new_line('a')// &
                    "ice_load = true"//new_line('a'))
   end function mass_nml

   function sf_fw_net_ref(engine) result(a)
      type(ocean_engine_t), intent(in) :: engine
      real(wp), allocatable :: a(:, :)

      a = engine%state%surface_flux%fw_net
   end function sf_fw_net_ref

   function mass_uniform_nml() result(nml)
      !! The uniform-atmosphere closed basin of `test_global_uniform` under real mass: cold air over a
      !! slightly warm ocean, so the ice grows on top and melts below, the wind moves and ridges it.
      character(len=:), allocatable :: nml

      nml = replace(uni_nml("basin"), "initial_temperature = -1.0", "initial_temperature = -0.6")
      nml = replace(nml, "T_init_surface = -1.0"//new_line('a')//"T_init_bottom = -1.0", "T_init_surface = -0.6"//new_line('a')//"T_init_bottom = -0.6")
      nml = replace(nml, "hfreeze = 10.0", "hfreeze = 10.0"//new_line('a')//"snowfall = 3.0e-5")
      nml = mass_nml(nml)
   end function mass_uniform_nml

   subroutine test_mass_water_closure(error)
      !! Ocean + ice + snow mass changes only by the atmosphere's exchange over ALL fractions (relative 1e-12), read
      !! from the ice's and the bulk forcing's own arrays; the ice's salt mass is the one named ledger term.
      type(error_type), allocatable, intent(out) :: error

      call global_case(error, mass_uniform_nml(), .false., 144, "real mass, uniform atmosphere, closed basin: ", &
                       mass=.true., what="water")
   end subroutine test_mass_water_closure

   subroutine test_mass_water_closure_prescribed(error)
      !! The same under the prescribed mixed atmosphere: rain and snow on both fractions, evaporation of the open water.
      type(error_type), allocatable, intent(out) :: error

      call global_case(error, mass_nml(pres_nml(1)), .true., 96, "real mass, prescribed mixed atmosphere: ", &
                       mass=.true., what="water")
   end subroutine test_mass_water_closure_prescribed

   subroutine test_mass_salt_closure(error)
      !! Ocean salt + ice salt is exact: the ice's salt moves as concentration, no volume, and no dilution term.
      type(error_type), allocatable, intent(out) :: error

      call global_case(error, mass_uniform_nml(), .false., 144, "real mass salt, uniform atmosphere: ", &
                       mass=.true., what="salt")
   end subroutine test_mass_salt_closure

   subroutine test_mass_heat_closure(error)
      !! Heat closes up to the named `H_mass_conv = cp T dM` term (the ocean's mass leaves at its own temperature,
      !! the ice's reference is 0 degC liquid), reported, within the vertical solver's bound.
      type(error_type), allocatable, intent(out) :: error

      call global_case(error, mass_uniform_nml(), .false., 144, "real mass heat, uniform atmosphere: ", &
                       mass=.true., what="heat")
   end subroutine test_mass_heat_closure

   function mass_ratio8_nml() result(nml)
      !! `mass_uniform_nml` at `dt_therm_ratio = 8` in the deferred order (`diabatic_first = .false.`).
      character(len=:), allocatable :: nml

      nml = replace(mass_uniform_nml(), "[ocean.vmix]"//new_line('a')// &
                    "dt_therm_ratio = 1"//new_line('a'), &
                    "[ocean.vmix]"//new_line('a')// &
                    "dt_therm_ratio = 8"//new_line('a')// &
                    "diabatic_first = false"//new_line('a'))
   end function mass_ratio8_nml

   subroutine test_mass_closure_ratio8(error)
      !! Water, salt and heat of ocean + ice at `dt_therm_ratio = 8`, deferred order, from the cold start on (18 thermo
      !! intervals; the first one included, where the ocean's mass follows the delivered flux arrays at every step).
      !! The ocean's atmosphere is the hand-off step's sample, held for the interval, so the check is against the
      !! delivered amounts (windows booked once, each delivered in the step after its hand-off).
      type(error_type), allocatable, intent(out) :: error

      call global_case(error, mass_ratio8_nml(), .false., 144, "real mass, ratio 8 deferred, uniform atmosphere: ", &
                       mass=.true., what="all")
   end subroutine test_mass_closure_ratio8

   subroutine test_rain_split(error)
      !! The open fraction's rain (`lprec`, weight `aice0`) and the covered fraction's (through the ice, in `fresh`,
      !! weight `aice_init`) are the whole precipitation exactly once: cell by cell the two weights add to 1 on wet
      !! cells and the delivered masses add to the nominal precipitation; the hook keeps the rain-through term.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      real(wp) :: t, open_rain, ice_rain, nominal_rain, open_snow, ice_snow, nominal_snow, wsum_mis, a
      real(wp) :: booked_rain, melt_mis
      real(wp), allocatable :: aice0_prev(:, :)
      integer :: n, i, j, ng, is, ie, js, je
      logical :: ok

      call write_synth_atm("mixed", FN_ATM, FN_LST)
      call make_engine(engine, cfg, mass_nml(pres_nml(1)), ok)
      call check(error, ok, "engine setup failed")
      if (allocated(error)) then
         call remove_synth_atm(FN_ATM, FN_LST)
         return
      end if
      ng = engine%grid%nghost
      is = ng + 1
      ie = ng + engine%grid%nx_phys
      js = ng + 1
      je = ng + engine%grid%ny_phys
      t = 0.0_wp
      open_rain = 0.0_wp
      ice_rain = 0.0_wp
      nominal_rain = 0.0_wp
      open_snow = 0.0_wp
      ice_snow = 0.0_wp
      nominal_snow = 0.0_wp
      wsum_mis = 0.0_wp
      booked_rain = 0.0_wp
      melt_mis = 0.0_wp
      call pull_cice(engine)
      aice0_prev = engine%state%ice%aice0
      do n = 1, 48
         call one_step(engine, cfg, t, ok)
         call check(error, ok, "step ran")
         if (allocated(error)) exit
         call pull_cice(engine)
         call pull_ice_flux(engine)
         call pull_sf(engine)
         do j = js, je
            do i = is, ie
               a = DT*engine%state%metrics%areaT(i, j)
               ! the bulk forcing weights by the aggregate of the PREVIOUS step's end, the ice by its start state
               if (engine%state%multilayer%wet_mask(i, j) > 0.5_wp) &
                  wsum_mis = max(wsum_mis, abs(aice0_prev(i, j) + engine%state%ice%cice%aice_init(i, j) - 1.0_wp))
               open_rain = open_rain + a*engine%state%surface_flux%lprec(i, j)
               open_snow = open_snow + a*engine%state%surface_flux%fprec(i, j)
               ice_rain = ice_rain + a*engine%state%icec%work%frain(i, j)*engine%state%ice%cice%aice_init(i, j)
               ice_snow = ice_snow + a*engine%state%icec%work%fsnow(i, j)*engine%state%ice%cice%aice_init(i, j)
               nominal_rain = nominal_rain + a*engine%state%icec%work%frain(i, j)
               nominal_snow = nominal_snow + a*engine%state%icec%work%fsnow(i, j)
               melt_mis = max(melt_mis, abs(engine%state%surface_flux%seaice_melt(i, j) - engine%state%icec%flux%fresh(i, j)))
            end do
         end do
         booked_rain = booked_rain + engine%state%icec%ledger%tab%rain
         aice0_prev = engine%state%ice%aice0
      end do
      if (.not. allocated(error)) then
  write (output_unit, "(a,3es13.5)") "  rain split: open (lprec), through the ice, nominal (kg):", open_rain, ice_rain, nominal_rain
       write (output_unit, "(a,3es13.5)") "  snow split: open (fprec), on the ice, nominal (kg):", open_snow, ice_snow, nominal_snow
       write (output_unit, "(a,es10.3,a,es13.5)") "  worst |aice0 + aice_init - 1|:", wsum_mis, "; ledger rain booked:", booked_rain
         call check(error, nominal_rain > 0.0_wp .and. nominal_snow > 0.0_wp, "the prescribed atmosphere precipitates")
         if (.not. allocated(error)) call check(error, wsum_mis < 1.0e-12_wp, "aice0 + aice_init = 1 on wet cells")
         if (.not. allocated(error)) call check(error, abs(open_rain + ice_rain - nominal_rain) < 1.0e-12_wp*nominal_rain, &
                                                "rain: open share + covered share = the whole precipitation")
         if (.not. allocated(error)) call check(error, abs(open_snow + ice_snow - nominal_snow) < 1.0e-12_wp*nominal_snow, &
                                                "snow: open share + covered share = the whole snowfall")
         if (.not. allocated(error)) call check(error, abs(booked_rain - ice_rain) < 1.0e-12_wp*nominal_rain, &
                                                "the rain through the ice is booked once, by the ice (kept in fresh)")
  if (.not. allocated(error)) call check(error, melt_mis <= 1.0e-12_wp*max(maxval(abs(engine%state%icec%flux%fresh)), 1.0e-30_wp), &
                                                "the hand-off keeps it: seaice_melt is fresh, rain through the ice included")
      end if
      call finish_engine(engine)
      call remove_synth_atm(FN_ATM, FN_LST)
   end subroutine test_rain_split

   function no_flow_nml(load, ratio, deferred, baroclinic) result(nml)
      !! A closed spherical basin of 12 x 10 cells, 100 m deep, three sigma layers, at rest, no wind; ice thin and
      !! complete in the northern half (`latitudes`), cold still air, an ocean at the freezing point: the ice
      !! grows, the open half does not.  `load` = the ice weight on the surface pressure (real mass); off, it is
      !! real mass WITHOUT the load: the configuration `validate_config` refuses (the negative control).
      !! `ratio`/`deferred`: `dt_therm_ratio` and `diabatic_first = .false.` (the deferred order).
      !! `baroclinic` (default off): the density follows T and S.  Off, the linear EOS has no T/S dependence, which
      !! isolates the BAROTROPIC balance of the load (the brine the freezing rejects would otherwise drive a real,
      !! slow baroclinic flow under the ice).
      logical, intent(in) :: load
      integer, intent(in), optional :: ratio
      logical, intent(in), optional :: deferred, baroclinic
      character(len=:), allocatable :: nml
      character(len=*), parameter :: NL = new_line("a")
      character(len=16) :: sratio
      character(len=:), allocatable :: vmix
      integer :: n
      logical :: defer, bcl

      n = 1
      if (present(ratio)) n = ratio
      defer = .false.
      if (present(deferred)) defer = deferred
      bcl = .false.
      if (present(baroclinic)) bcl = baroclinic
      write (sratio, '(i0)') n
      vmix = "[ocean.vmix]"//NL//"dt_therm_ratio = "//trim(sratio)//NL
      if (defer) vmix = vmix//"diabatic_first = false"//NL
      nml = "[sim]"//NL// &
            'sim_type = "ocean"'//NL// &
            "[time]"//NL// &
            "t_end = 86400.0"//NL// &
            "dt_fixed = 600.0"//NL// &
            "[nonhydrostatic]"//NL// &
            "nz_layers = 3"//NL// &
            "[ocean.diag]"//NL// &
            "enabled = false"//NL// &
            "[tracer]"//NL// &
            "initial_temperature = -1.80"//NL// &
            "initial_salinity = 34.0"//NL// &
            "T_init_surface = -1.80"//NL// &
            "T_init_bottom = -1.80"//NL// &
            "[ocean.bt]"//NL// &
            "auto_n_inner = true"//NL// &
            "[ocean.hvisc]"//NL// &
            "nu_h = 200.0"//NL// &
            "[grid]"//NL// &
            "nx = 12"//NL// &
            "ny = 10"//NL// &
            "nghost = 3"//NL// &
            "dx = 0.25"//NL// &
            "dy = 0.25"//NL// &
            "[ocean.grid]"//NL// &
            'grid_config = "spherical"'//NL// &
            "lon_west = 20.0"//NL// &
            "lat_south = 70.0"//NL// &
            "rad_earth = 6.371e6"//NL// &
            'coriolis_scheme = "planetary"'//NL// &
            "[physics]"//NL// &
            "wind_stress_x = 0.0"//NL// &
            "wind_stress_y = 0.0"//NL// &
            "[vcoord]"//NL// &
            'vcoord_type = "sigma"'//NL// &
            "[ocean.topo]"//NL// &
            'topo_config = "flat"'//NL// &
            "max_depth = 100.0"//NL// &
            "[ocean.thermo]"//NL// &
            "enable_thermodynamics = true"//NL// &
            "[ocean.forcing]"//NL// &
            "enable_components = true"//NL// &
            'freshwater = "mass"'//NL// &
            vmix//NL// &
            "[ocean.ice]"//NL// &
            "enable = true"//NL// &
            'scheme = "cice"'//NL// &
            "ncat = 5"//NL// &
            "nilyr = 4"//NL// &
            "hfreeze = 10.0"//NL// &
            "dynamics = true"//NL// &
            "transport = true"//NL// &
            "[ocean.ice_thermo]"//NL// &
            "update_ocn_f = true"//NL// &
            "[ocean.ice_dyn]"//NL// &
            'advection = "upwind"'//NL// &
            "ndte = 40"//NL// &
            "[ocean.ice_ridge]"//NL// &
            "kstrength = 0"//NL// &
            "Pstar = 1.0"//NL// &
            "[ocean.ice_ic]"//NL// &
            'conc_config = "latitudes"'//NL// &
            "arctic_edge = 71.2"//NL// &
            "h_ice = 0.05"//NL// &
            "conc = 1.0"//NL// &
            "h_snow = 0.0"//NL// &
            "t_ice = -8.0"//NL// &
            "[ocean.ice_atm]"//NL// &
            "uni_tas = 240.0"//NL// &
            "uni_rlds = 120.0"//NL// &
            "uni_uas = 0.0"//NL// &
            "uni_vas = 0.0"//NL// &
            "[ocean.bc]"//NL// &
            'west = "wall"'//NL// &
            'east = "wall"'//NL// &
            'south = "wall"'//NL// &
            'north = "wall"'//NL// &
            "[output]"//NL// &
            "output_to_file = false"//NL
      if (.not. bcl) nml = nml//NL// &
                           "[ocean.ic]"//NL// &
                           "alpha_T = 0.0"//NL// &
                           "beta_S = 0.0"//NL
      if (load) then
         ! `ice_load_balance_ic` is on by default: the cold start is lowered by the ice's weight, no hand IC
         nml = nml//NL// &
               "[ocean.psurf]"//NL// &
               "enable = true"//NL// &
               "ice_load = true"//NL
      end if
   end function no_flow_nml

   subroutine no_flow_run(error, load, nsteps, umax, dh_ice, depression, u_first, ratio, deferred, baroclinic, d_atm)
      !! `umax`: the largest face velocity over the whole run.  `u_first`: over the steps up to and including the
      !! first one that DELIVERS the ice's water (step `ratio` deferred, `ratio + 1` in-stage, 2 at ratio 1).
      !! `d_atm`: the largest sea-level step (m) the ICE's EXCHANGE WITH THE ATMOSPHERE made over the run,
      !! `sum_steps |evap| dt / rho_0` on the worst column: sublimation and condensation change the ice's weight
      !! (`fresh` contains `evap`, so the ocean's water does not follow) and the ocean answers that real change.
      type(error_type), allocatable, intent(out) :: error
      logical, intent(in) :: load
      integer, intent(in) :: nsteps
      real(wp), intent(out) :: umax, dh_ice, depression, u_first
      integer, intent(in), optional :: ratio
      logical, intent(in), optional :: deferred, baroclinic
      real(wp), intent(out), optional :: d_atm
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      real(wp) :: t, v0, v1
      real(wp), allocatable :: cum(:, :)
      integer :: n, n_first, nr
      logical :: ok, defer

      umax = 0.0_wp
      u_first = 0.0_wp
      dh_ice = 0.0_wp
      depression = 0.0_wp
      nr = 1
      if (present(ratio)) nr = ratio
      defer = .false.
      if (present(deferred)) defer = deferred
      n_first = 2
      if (nr > 1) n_first = merge(nr, nr + 1, defer)
      call make_engine(engine, cfg, no_flow_nml(load, ratio, deferred, baroclinic), ok)
      call check(error, ok, "engine setup failed")
      if (allocated(error)) return
      if (load) then
         !$acc update self(engine%state%surface_flux%p_ice_load)
         depression = maxval(engine%state%surface_flux%p_ice_load)/(engine%state%eos%rho0*engine%state%dyn%bt_work%g_bt)
      end if
      call pull_cice(engine)
      v0 = maxval(engine%state%ice%vice)
      allocate (cum(size(engine%state%metrics%areaT, 1), size(engine%state%metrics%areaT, 2)), source=0.0_wp)
      t = 0.0_wp
      do n = 1, nsteps
         call one_step(engine, cfg, t, ok)
         call check(error, ok, "step ran")
         if (allocated(error)) exit
         !$acc update self(engine%state%multilayer%u_face_x_layer, engine%state%multilayer%v_face_y_layer)
         !$acc update self(engine%state%icec%flux%evap)
         umax = max(umax, maxval(abs(engine%state%multilayer%u_face_x_layer)), &
                    maxval(abs(engine%state%multilayer%v_face_y_layer)))
         cum = cum + DT*abs(engine%state%icec%flux%evap)
         if (n == n_first) u_first = umax
      end do
      call pull_cice(engine)
      v1 = maxval(engine%state%ice%vice)
      dh_ice = v1 - v0
      if (present(d_atm)) d_atm = maxval(cum)/engine%state%eos%rho0
      call finish_engine(engine)
   end subroutine no_flow_run

   subroutine test_ice_load_no_flow(error)
      !! Freeze half a closed basin at rest.  With the ice's weight on the surface pressure the water the ice takes
      !! (`rho_i dV/rho_0`) and the load (`g rho_i dV`) change on the SAME outer step (the ice's exchange with the
      !! atmosphere excepted, see (2)), so the ocean stays at rest.
      !! (1) THE TIMING: through the first delivery of the ice's water max|u| is round-off, at `dt_therm_ratio = 1`,
      !! at ratio 8 in the deferred order (`diabatic_first = .false.`, where the delivery is one step before the
      !! interval's end) and at ratio 2 in the in-stage order, from the cold start the `ice_load_balance_ic` lowered
      !! (no hand IC).  The imbalance a mistimed load would leave is the ice's whole growth, 250 times the next
      !! term.  (2) THE WHOLE RUN: what remains afterwards is not the delivery but the ice's own exchange with the
      !! ATMOSPHERE (sublimation / condensation: `fresh` contains `evap`, so the ocean's water does not follow the
      !! weight that the atmosphere takes), a real sea-level step of `d_atm` that the ocean answers at gravity-wave
      !! speed `sqrt(g/H)`: the run's max|u| is held to `sqrt(g/H) d_atm` (barotropic: density independent of T, S).
      !! NEGATIVE CONTROL: the same real-mass run without the load (`validate_config` refuses that; the test builds
      !! it directly) removes the water from under the ice only, a sea-surface step of the growth at the ice edge,
      !! and the flow it drives is orders of magnitude larger.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: u_off, g_off, dep_off, uf_off
      real(wp) :: u1, g1, d1, uf1, a1, u8, g8, d8, uf8, a8, u2, g2, d2, uf2, a2, c
      real(wp), parameter :: ROUNDOFF = 1.0e-11_wp

      call no_flow_run(error, .true., 72, u1, g1, d1, uf1, ratio=1, d_atm=a1)
      if (allocated(error)) return
      call no_flow_run(error, .true., 72, u8, g8, d8, uf8, ratio=8, deferred=.true., d_atm=a8)
      if (allocated(error)) return
      call no_flow_run(error, .true., 72, u2, g2, d2, uf2, ratio=2, d_atm=a2)
      if (allocated(error)) return
      call no_flow_run(error, .false., 72, u_off, g_off, dep_off, uf_off)
      if (allocated(error)) return
      c = sqrt(9.8_wp/100.0_wp)
      write (output_unit, "(a)") "  ice_load_no_flow: max|u| over the run (m/s), through the first delivery, ice growth (m), IC depression (m), atmosphere sea-level step (m)"
      write (output_unit, "(a,4es11.3,f8.4)") "    ratio 1:               ", u1, uf1, g1, a1, d1
      write (output_unit, "(a,4es11.3,f8.4)") "    ratio 8, deferred:     ", u8, uf8, g8, a8, d8
      write (output_unit, "(a,4es11.3,f8.4)") "    ratio 2, in-stage:     ", u2, uf2, g2, a2, d2
      write (output_unit, "(a,3es11.3)") "    negative control (no load):", u_off, uf_off, g_off
      call check(error, g1 > 1.0e-3_wp .and. g8 > 1.0e-3_wp .and. g2 > 1.0e-3_wp .and. g_off > 1.0e-3_wp, &
                 "the ice grows (the test freezes water)")
      if (allocated(error)) return
      call check(error, u_off > 1.0e-4_wp, "negative control: without the load the freezing drives a flow")
      if (allocated(error)) return
      call check(error, d1 > 1.0e-3_wp, "the IC balance lowered the start under the ice")
      if (allocated(error)) return
      call check(error, uf1 < ROUNDOFF, "ratio 1: through the first delivery the loaded ocean is at rest to round-off")
      if (allocated(error)) return
   call check(error, uf8 < ROUNDOFF, "ratio 8, deferred order: through the first delivery the loaded ocean is at rest to round-off")
      if (allocated(error)) return
   call check(error, uf2 < ROUNDOFF, "ratio 2, in-stage order: through the first delivery the loaded ocean is at rest to round-off")
      if (allocated(error)) return
      call check(error, u1 < c*a1 .and. u8 < c*a8 .and. u2 < c*a2, &
                 "over the whole run the flow is below the gravity-wave answer to the atmosphere's own sea-level step")
      if (allocated(error)) return
      call check(error, u_off > 100.0_wp*max(u1, u8, u2), "the load removes the flow (control / loaded >= 100)")
   end subroutine test_ice_load_no_flow

   subroutine test_ice_load_drift(error)
      !! WHERE the residual flow of the loaded, freezing basin comes from (it was 3.9e-6 m/s and attributed to the
      !! drifting ice).  The ice drifts at about 6e-6 m/s here, and mass conservation of a load of equivalent water
      !! depth `d` moving at `v` over a column of depth `H` bounds the flow it stirs by `u ~ (d/H) v`, ~5e-9 m/s: the
      !! drift is NOT the residual.  The barotropic basin's residual equals the run with a different density law, i.e. it
      !! is neither the drift nor the brine's baroclinic flow (both are tested below as a bound), it is the sea-level
      !! step of the ice's sublimation (`ice_load_no_flow`).
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: u_bt, g_bt_, dep_bt, uf_bt, u_bcl, g_bcl, dep_bcl, uf_bcl, a_bt, a_bcl
      real(wp), parameter :: ICE_SPEED_OBS = 6.0e-6_wp
         !! the fastest ice of this basin over the run, measured (it follows the ocean's own flow)

      call no_flow_run(error, .true., 72, u_bt, g_bt_, dep_bt, uf_bt, d_atm=a_bt)
      if (allocated(error)) return
      call no_flow_run(error, .true., 72, u_bcl, g_bcl, dep_bcl, uf_bcl, baroclinic=.true., d_atm=a_bcl)
      if (allocated(error)) return
      write (output_unit, "(a,es10.3,a,es10.3,a,es10.3)") "  ice_load_drift: barotropic max|u| =", u_bt, " m/s, baroclinic (T/S density) ", u_bcl, &
         " m/s; the drift's bound (d/H) v:", (dep_bt/100.0_wp)*ICE_SPEED_OBS
      call check(error, u_bcl < 1.5_wp*u_bt + 1.0e-6_wp, "the brine's baroclinic flow does not exceed the barotropic residual by more than 1e-6 m/s")
      if (allocated(error)) return
      call check(error, (dep_bt/100.0_wp)*ICE_SPEED_OBS < 1.0e-3_wp*u_bt, &
                 "the drift's bound (d/H) v is three orders below the residual: the drift is not the residual")
   end subroutine test_ice_load_drift

   subroutine psurf_seam_case(error, mass, atm_load)
      !! The assembled `p_surf` on a periodic east-west channel under the bulk forcing with a latitude-dependent
      !! sea-level pressure.  `mass`: real mass with the ice weight (else virtual fresh water, no ice load);
      !! `atm_load`: the atmosphere's pressure joins `p_surf`.
      type(error_type), allocatable, intent(out) :: error
      logical, intent(in) :: mass, atm_load
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      character(len=:), allocatable :: nml
      real(wp) :: t, lat, psl, seam, diffuse, rho0g, ib
      integer :: n, i, j, ng, nxp, nyp, g
      logical :: ok

      if (mass) then
         nml = mass_nml(pres_nml(1))
         if (atm_load) nml = replace(nml, "ice_load = true", "ice_load = true"//new_line('a')//"atm_load = true")
      else
         nml = pres_nml(1)
         nml = replace(nml, "[ocean.forcing]"//new_line('a')// &
                       "enable_components = true"//new_line('a'), &
                       "[ocean.forcing]"//new_line("a")// &
                       "enable_components = true"//new_line("a")// &
                       "[ocean.psurf]"//new_line("a")// &
                       "enable = true"//new_line("a")// &
                       "atm_load = true"//new_line("a"))
      end if
      nml = replace(nml, 'west = "wall"'//new_line('a')//'east = "wall"', 'west = "periodic"'//new_line('a')//'east = "periodic"')
      call write_synth_atm("psl_lat", FN_ATM, FN_LST)
      call make_engine(engine, cfg, nml, ok)
      call check(error, ok, "engine setup failed")
      if (allocated(error)) then
         call remove_synth_atm(FN_ATM, FN_LST)
         return
      end if
      ng = engine%grid%nghost
      nxp = engine%grid%nx_phys
      nyp = engine%grid%ny_phys
      t = 0.0_wp
      do n = 1, 4
         call one_step(engine, cfg, t, ok)
         call check(error, ok, "step ran")
         if (allocated(error)) exit
      end do
      if (.not. allocated(error)) then
         call pull_sf(engine)
         !$acc update self(engine%state%surface_flux%p_surf_atm, engine%state%surface_flux%p_surf)
         !$acc update self(engine%state%p_surf%eta_ib)
         if (mass) then
            !$acc update self(engine%state%surface_flux%p_ice_load)
         end if
         rho0g = engine%state%eos%rho0*engine%state%dyn%bt_work%g_bt
         seam = 0.0_wp
         diffuse = 0.0_wp
         ib = 0.0_wp
         associate (sf => engine%state%surface_flux)
            call check(error, maxval(sf%p_surf_atm) > 7.0e4_wp, "the bulk forcing wrote the sea-level pressure")
            do j = ng + 1, ng + nyp
               lat = 70.0_wp + (real(j - ng, wp) - 0.5_wp)*0.25_wp
               psl = 79000.0_wp + 300.0_wp*lat
               do i = 1, size(sf%p_surf, 1)
                  if (atm_load) then
                     if (mass) then
                        diffuse = max(diffuse, abs(sf%p_surf(i, j) - psl - sf%p_ice_load(i, j)))
                     else
                        diffuse = max(diffuse, abs(sf%p_surf(i, j) - psl))
                     end if
                     if (.not. mass) ib = max(ib, abs(engine%state%p_surf%eta_ib(i, j) + psl/rho0g))
                  else
                     ! the ice weight is the ONLY addition to the configure seed (0 Pa): the atmosphere's
                     ! pressure never reaches p_surf, owned cells or ghosts
                     diffuse = max(diffuse, abs(sf%p_surf(i, j) - sf%p_ice_load(i, j)))
                  end if
               end do
               do g = 1, ng
                  seam = max(seam, abs(sf%p_surf(g, j) - sf%p_surf(g + nxp, j)), &
                             abs(sf%p_surf(ng + nxp + g, j) - sf%p_surf(ng + g, j)))
               end do
            end do
            write (output_unit, "(a,l1,a,l1,a,3es11.3)") "  psurf seam: mass ", mass, ", atm_load ", atm_load, &
               "; worst |p_surf - expected|, periodic-seam jump (Pa), |eta_ib + psl/(rho0 g)| (m):", diffuse, seam, ib
            if (.not. allocated(error)) call check(error, diffuse < 1.0e-4_wp, &
                                                   "p_surf is the expected sum of its parts, ghost columns included")
            if (.not. allocated(error)) call check(error, seam == 0.0_wp, &
                                                   "p_surf has no discontinuity across the periodic seam")
            if (.not. allocated(error) .and. atm_load .and. .not. mass) &
               call check(error, ib < 1.0e-9_wp, "the inverse-barometer response is -psl/(rho0 g_bt) on every column")
         end associate
      end if
      call finish_engine(engine)
      call remove_synth_atm(FN_ATM, FN_LST)
   end subroutine psurf_seam_case

   subroutine test_psurf_ice_only(error)
      !! `ice_load` adds ONLY the ice's weight: the bulk forcing's sea-level pressure (1e5 Pa) must not reach `p_surf`,
      !! and the periodic seam is continuous (it was a ~10 m inverse-barometer step before: owned cells held `psl`,
      !! ghosts held the seed).
      type(error_type), allocatable, intent(out) :: error

      call psurf_seam_case(error, .true., .false.)
   end subroutine test_psurf_ice_only

   subroutine test_psurf_atm_load_ice(error)
      !! `atm_load` + `ice_load`: `p_surf = psl + ice weight`, with the atmosphere's seam ghosts exchanged.
      type(error_type), allocatable, intent(out) :: error

      call psurf_seam_case(error, .true., .true.)
   end subroutine test_psurf_atm_load_ice

   subroutine test_psurf_atm_load_ib(error)
      !! `atm_load` alone: the analytic inverse barometer `eta_ib = -psl/(rho0 g_bt)` on every column, ghosts included.
      type(error_type), allocatable, intent(out) :: error

      call psurf_seam_case(error, .false., .true.)
   end subroutine test_psurf_atm_load_ib

   subroutine test_ice_load_restart(error)
      !! N steps + checkpoint + N steps == 2N, bitwise, with the ice loading on (the carried load `sf_p_ice_load`).
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: ea, eb
      type(config_t) :: cfg_a, cfg_b
      real(wp) :: ta, tb, t0
      integer :: step0, n, ierr
      logical :: ok
      character(len=:), allocatable :: nml

      nml = no_flow_nml(.true.)
      nml = replace(nml, "uni_tas = 240.0", "uni_tas = 235.0")
      call make_engine(ea, cfg_a, nml, ok)
      call check(error, ok, "engine A setup failed")
      if (allocated(error)) return
      ta = 0.0_wp
      do n = 1, 6
         call one_step(ea, cfg_a, ta, ok)
         call check(error, ok, "A step ran")
         if (allocated(error)) then
            call finish_engine(ea)
            return
         end if
      end do
      call ocean_state_restart_write(ea%state, ea%grid, ea%decomp, FN_RT, ta, 6, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "checkpoint write failed")
      if (allocated(error)) then
         call finish_engine(ea)
         return
      end if
      do n = 1, 6
         call one_step(ea, cfg_a, ta, ok)
      end do
      call pull_sf(ea)
      call make_engine(eb, cfg_b, nml, ok, restart_file=FN_RT, t0=t0, step0=step0)
      call check(error, ok, "engine B (restart) setup failed")
      if (allocated(error)) then
         call finish_engine(ea)
         return
      end if
      tb = t0
      do n = 1, 6
         call one_step(eb, cfg_b, tb, ok)
      end do
      call pull_sf(eb)
      !$acc update self(ea%state%multilayer%h_layer, eb%state%multilayer%h_layer)
      call check(error, all(ea%state%surface_flux%p_surf == eb%state%surface_flux%p_surf), "p_surf is bitwise after the resume")
      if (.not. allocated(error)) call check(error, all(ea%state%surface_flux%p_ice_load == eb%state%surface_flux%p_ice_load), &
                                             "the carried ice load is bitwise after the resume")
      if (.not. allocated(error)) call check(error, all(ea%state%multilayer%h_layer == eb%state%multilayer%h_layer), &
                                             "h_layer is bitwise after the resume")
      if (.not. allocated(error)) call check(error, maxval(ea%state%surface_flux%p_ice_load) > 0.0_wp, "the load is live")
      call finish_engine(eb)
      call finish_engine(ea)
      call delete_file(FN_RT)
   end subroutine test_ice_load_restart

   ! ------------------------------------------------------------------
   ! D: dynamics-only invariance
   ! ------------------------------------------------------------------

   subroutine dynamics_case(error, with_snow)
      type(error_type), allocatable, intent(out) :: error
      logical, intent(in) :: with_snow
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      type(stage_run_t) :: run
      character(len=:), allocatable :: nml
      real(wp) :: t, m0, s0, e0, m1, s1, e1, vmax, ridge_max, worst_dyn, worst_tr, worst_rdg, fl_abs, fresh_rdg, worst_res
      integer :: n, ng, is, ie, js, je
      logical :: ok

      nml = uni_nml("basin")
      if (.not. with_snow) nml = replace(nml, "h_snow = 0.1", "h_snow = 0.0")
      nml = replace(nml, "wind_stress_x = 0.1"//new_line('a')//"wind_stress_y = 0.03", "wind_stress_x = 0.2"//new_line('a')//"wind_stress_y = 0.06")
      call make_engine(engine, cfg, nml, ok)
      call check(error, ok, "engine setup failed")
      if (allocated(error)) return
      engine%state%icec%skip_thermo = .true.
      ng = engine%grid%nghost
      is = ng + 1
      ie = ng + engine%grid%nx_phys
      js = ng + 1
      je = ng + engine%grid%ny_phys
      call engine%state%icec%ledger%totals(engine%state%ice%cice, m0, s0, e0)
      run = stage_run_t()
      t = 0.0_wp
      vmax = 0.0_wp
      ridge_max = 0.0_wp
      worst_dyn = 0.0_wp
      worst_tr = 0.0_wp
      worst_rdg = 0.0_wp
      worst_res = 0.0_wp
      fl_abs = 0.0_wp
      fresh_rdg = 0.0_wp
      do n = 1, 60
         call one_step(engine, cfg, t, ok)
         call check(error, ok, "step ran")
         if (allocated(error)) exit
         call pull_cice(engine)
         call pull_ice_flux(engine)
         !$acc update self(engine%state%ice%u_ice, engine%state%ice%v_ice)
         call accumulate_tab(engine, run)
         vmax = max(vmax, maxval(abs(engine%state%ice%u_ice)))
         ridge_max = max(ridge_max, maxval(engine%state%icec%flux%dardg1dt(is:ie, js:je)))
         associate (L => engine%state%icec%ledger)
            worst_dyn = max(worst_dyn, maxval(abs(L%tab%dstate(:, ICE_LT_DYNAMICS))))
            worst_tr = max(worst_tr, abs(L%tab%dstate(1, ICE_LT_TRANSPORT))/max(m0, 1.0_wp), &
                           abs(L%tab%dstate(2, ICE_LT_TRANSPORT))/max(s0, 1.0e-30_wp), &
                           abs(L%tab%dstate(3, ICE_LT_TRANSPORT))/max(abs(e0), 1.0_wp))
            worst_rdg = max(worst_rdg, abs(L%tab%res(1, ICE_LT_RIDGING))/max(m0, 1.0_wp))
            worst_res = max(worst_res, maxval(abs(L%tab%res(1, :)))/max(m0, 1.0_wp), &
                            maxval(abs(L%tab%res(2, :)))/max(s0, 1.0e-30_wp), maxval(abs(L%tab%res(3, :)))/max(abs(e0), 1.0_wp))
            fl_abs = fl_abs + sum(abs(L%tab%fout(:, ICE_LT_DYNAMICS))) + sum(abs(L%tab%fout(:, ICE_LT_THERM1))) &
                     + abs(L%tab%snow) + abs(L%tab%lw_in) + abs(L%tab%atm_rad)
            fresh_rdg = fresh_rdg + L%tab%fout(1, ICE_LT_RIDGING)
         end associate
      end do
      if (.not. allocated(error)) then
         call engine%state%icec%ledger%totals(engine%state%ice%cice, m1, s1, e1)
         write (output_unit, "(a,l1,a)") "  dynamics only (snow = ", with_snow, "):"
        write (output_unit, "(a,3es11.3)") "    relative change of mass, salt, energy over 60 steps:", (m1 - m0)/m0, (s1 - s0)/s0, &
            (e1 - e0)/max(abs(e0), 1.0_wp)
         write (output_unit, "(a,es11.3,a,es11.3,a,es11.3)") "    worst per-step dynamics change", worst_dyn, &
            ", transport change (relative)", worst_tr, ", ridging mass residual (relative)", worst_rdg
         write (output_unit, "(a,es11.3,a,es11.3,a,es10.3)") "    ridged snow booked to the ocean (kg)", fresh_rdg, &
            "; ice speed max", vmax, "; ridging rate max"//" ", ridge_max
         write (output_unit, "(a,3es11.3,a,3es11.3)") "    ocean booking |water kg, salt kg, heat J| transport:", run%absfout(:, ICE_LT_TRANSPORT), &
            "; ridging:", run%absfout(:, ICE_LT_RIDGING)
         write (output_unit, "(a,es11.3)") "    worst stage residual relative to the totals:", worst_res
         call check(error, run%finite, "finite")
         if (.not. allocated(error)) call check(error, vmax > 1.0e-6_wp .and. ridge_max > 0.0_wp, "the ice moves and ridges")
         if (.not. allocated(error)) call check(error, fl_abs == 0.0_wp, &
                                                "no atmosphere or ocean line is booked outside the ridging stage")
         if (.not. allocated(error)) call check(error, worst_tr < 1.0e-14_wp, &
                                                "EVP + transport (one probe spans both) conserve mass, salt and energy")
         if (.not. allocated(error)) call check(error, abs(m1 - m0 + sum(run%fout(1, :)))/m0 < 1.0e-14_wp, &
                                                "mass: the ice changes only by what the ocean was booked (the ridged snow)")
         if (.not. allocated(error)) call check(error, abs(s1 - s0)/s0 < 1.0e-14_wp, "salt conserved")
         if (.not. allocated(error)) call check(error, abs(e1 - e0 + sum(run%fout(3, :)))/max(abs(e0), 1.0_wp) < 1.0e-14_wp, &
                                                "energy: the ice changes only by the ridged snow's enthalpy booked to the ocean")
         if (.not. allocated(error)) call check(error, worst_res < 1.0e-14_wp, &
                                                "every stage line closes, relative to the totals of the conserved quantity")
         if (.not. allocated(error) .and. .not. with_snow) then
            call check(error, sum(run%absfout(:, ICE_LT_RIDGING)) + sum(run%absfout(:, ICE_LT_TRANSPORT)) == 0.0_wp, &
                       "without snow neither the transport nor the ridging books anything to the ocean")
            if (.not. allocated(error)) call check(error, abs(m1 - m0)/m0 < 1.0e-14_wp, "the area-weighted volume is conserved")
         end if
         if (.not. allocated(error) .and. with_snow) call check(error, run%fout(1, ICE_LT_RIDGING) > 0.0_wp, &
                                                                "the ridged snow reached the ocean's ledger")
      end if
      call finish_engine(engine)
   end subroutine dynamics_case

   subroutine test_dynamics_only(error)
      type(error_type), allocatable, intent(out) :: error

      call dynamics_case(error, .false.)
   end subroutine test_dynamics_only

   subroutine test_dynamics_only_snow(error)
      type(error_type), allocatable, intent(out) :: error

      call dynamics_case(error, .true.)
   end subroutine test_dynamics_only_snow

   ! ------------------------------------------------------------------
   ! the knob off, and the default scheme: nothing runs, nothing is registered
   ! ------------------------------------------------------------------

   subroutine test_ledger_off_inert(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      type(restart_registry_t) :: reg
      character(len=:), allocatable :: nml
      real(wp) :: t
      integer :: e, n_led, w
      logical :: ok

      do w = 1, 2
         if (w == 1) then
            ! the CICE scheme with the ledger knob off
            nml = replace(uni_nml("column"), new_line('a')//"ledger = true", "")
         else
            ! the default scheme
            nml = replace(uni_nml("column"), 'scheme = "cice"'//new_line('a')//"ncat = 5"//new_line('a')//"nilyr = 4"//new_line('a')//"hfreeze = 10.0"//new_line('a')//"dynamics = true"//new_line('a')//&
                          "transport = true"//new_line('a')//"ledger = true", "ncat = 5"//new_line('a')//"nk_ice = 2")
            nml = replace(nml, "[ocean.ice_dyn]"//new_line('a')// &
                          'advection = "upwind"'//new_line('a')// &
                          "ndte = 40"//new_line('a'), "[ocean.ice_dyn]"//new_line('a')// &
                          "ndte = 40"//new_line('a'))
         end if
         call make_engine(engine, cfg, nml, ok)
         call check(error, ok, "engine setup failed")
         if (allocated(error)) return
         call check(error,.not. engine%state%icec%ledger%on .and. .not. engine%state%icec%ledger%is_init, &
                    "the ledger slot is inert")
         if (.not. allocated(error)) call check(error, engine%state%icec%ledger%bytes() == 0_int64, "and counts no bytes")
         if (.not. allocated(error)) then
            t = 0.0_wp
            call advance(engine, cfg, t, 4, ok)
            call check(error, ok, "4 steps ran")
         end if
         if (.not. allocated(error)) then
            call ice_cice_ledger_report(engine%state, 0)    ! must print nothing and touch nothing
            call ocean_state_build_restart_registry(engine%state, engine%grid, reg)
            n_led = 0
            do e = 1, reg%n
               if (index(reg%entries(e)%tag, "ice_led_") == 1) n_led = n_led + 1
            end do
            call check(error, n_led == 0, "no restart tag of the ledger")
         end if
         call finish_engine(engine)
         if (allocated(error)) return
      end do
   end subroutine test_ledger_off_inert

   ! ------------------------------------------------------------------
   ! E: restart
   ! ------------------------------------------------------------------

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
            if (snap%n >= MAXF) error stop "test_cice_closure: raise MAXF"
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

   subroutine compare(error, sa, sb, when, n_led)
      type(error_type), allocatable, intent(out) :: error
      type(snapshot_t), intent(in) :: sa, sb
      character(len=*), intent(in) :: when
      integer, intent(out) :: n_led
      integer :: e, i, j, k, nbad
      character(len=:), allocatable :: msg

      n_led = 0
      nbad = 0
      msg = ""
      call check(error, sa%n == sb%n, "registry field count differs "//when)
      if (allocated(error)) return
      do e = 1, sa%n
         if (sa%f(e)%tag /= sb%f(e)%tag) then
            call check(error, .false., "registry field order differs: "//trim(sa%f(e)%tag))
            return
         end if
         if (index(sa%f(e)%tag, "ice_led_") == 1) n_led = n_led + 1
         if (any(shape(sa%f(e)%a) /= shape(sb%f(e)%a))) then
            nbad = nbad + 1
            msg = msg//" "//trim(sa%f(e)%tag)//"(shape)"
            cycle
         end if
         do k = 1, size(sa%f(e)%a, 3)
            do j = 1, size(sa%f(e)%a, 2)
               do i = 1, size(sa%f(e)%a, 1)
                  if (transfer(sa%f(e)%a(i, j, k), 0_int64) /= transfer(sb%f(e)%a(i, j, k), 0_int64)) then
                     nbad = nbad + 1
                     msg = msg//" "//trim(sa%f(e)%tag)
                     exit
                  end if
               end do
            end do
         end do
      end do
      call check(error, nbad == 0, "fields differ "//when//":"//msg)
   end subroutine compare

   subroutine strip_ledger(snap)
      !! Remove the ledger's registry entries (`ice_led_*`) from a snapshot.
      type(snapshot_t), intent(inout) :: snap
      integer :: e, m

      m = 0
      do e = 1, snap%n
         if (index(snap%f(e)%tag, "ice_led_") == 1) cycle
         m = m + 1
         if (m /= e) snap%f(m) = snap%f(e)
      end do
      snap%n = m
   end subroutine strip_ledger

   subroutine test_ledger_on_equals_off(error)
      !! The ledger is an OBSERVER: with it on the step runs therm2 as four sweeps and probes between them, and
      !! the whole restart registry (everything but `ice_led_*`) is bitwise what the default run holds.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: ea, eb
      type(config_t) :: cfg_a, cfg_b
      type(snapshot_t) :: sa, sb
      character(len=:), allocatable :: nml_on, nml_off
      real(wp) :: t_a, t_b
      integer :: n_led
      logical :: ok

      nml_on = replace(uni_nml("basin"), "dt_therm_ratio = 1", "dt_therm_ratio = 2")
      nml_off = replace(nml_on, new_line('a')//"ledger = true", "")
      call make_engine(ea, cfg_a, nml_off, ok)
      call check(error, ok, "engine (ledger off) setup failed")
      if (allocated(error)) return
      t_a = 0.0_wp
      call advance(ea, cfg_a, t_a, 24, ok)
      call check(error, ok, "24 steps, ledger off")
      if (.not. allocated(error)) then
         call take_snapshot(ea, sa)
      end if
      call finish_engine(ea)
      if (allocated(error)) return
      call make_engine(eb, cfg_b, nml_on, ok)
      call check(error, ok, "engine (ledger on) setup failed")
      if (allocated(error)) return
      t_b = 0.0_wp
      call advance(eb, cfg_b, t_b, 24, ok)
      call check(error, ok, "24 steps, ledger on")
      if (.not. allocated(error)) then
         call take_snapshot(eb, sb)
         call strip_ledger(sb)
         call compare(error, sa, sb, "ledger on vs off", n_led)
      end if
      call finish_engine(eb)
   end subroutine test_ledger_on_equals_off

   subroutine test_ledger_nan_sticky(error)
      !! A non-finite step must not vanish from the console: the accumulators turn NaN for good, the counter counts.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      real(wp) :: t, n0
      character(len=200) :: line
      logical :: ok

      call make_engine(engine, cfg, uni_nml("column"), ok)
      call check(error, ok, "engine setup failed")
      if (allocated(error)) return
      t = 0.0_wp
      call advance(engine, cfg, t, 3, ok)
      call check(error, ok, "3 steps ran")
      if (.not. allocated(error)) then
         associate (L => engine%state%icec%ledger)
            line = L%console_line()
            call check(error,.not. has_nan(line), "finite before the poisoning")
            n0 = L%nsteps
            ! a poisoned summand in the new-ice stage's mass total (what a NaN in the state decomposes to)
            L%lv = efp_t()
            L%probed = .false.
            L%probed(ICE_LT_NEWICE) = .true.
            L%lv(ice_lt_ix_state(1, ICE_LT_NEWICE))%poison = 1_int64
            L%is_open = .true.
            call L%end_step()
            line = L%console_line()
            write (output_unit, "(a,a)") "  poisoned: ", trim(line)
            if (.not. allocated(error)) call check(error, has_nan(line), "the console shows NaN")
            if (.not. allocated(error)) call check(error,.not. ieee_is_finite(L%cum_rel(ICE_LT_WATER)), "cum_rel is NaN")
            if (.not. allocated(error)) call check(error,.not. any(ieee_is_finite(L%worst)), "worst is NaN")
            if (.not. allocated(error)) call check(error, L%nsteps == n0 + 1.0_wp, "the step is counted")
            ! a later finite step does not wash it out
            if (.not. allocated(error)) then
               call one_step(engine, cfg, t, ok)
               line = L%console_line()
               call check(error, ok .and. has_nan(line), "still NaN after a finite step")
            end if
         end associate
      end if
      call finish_engine(engine)
   end subroutine test_ledger_nan_sticky

   subroutine restart_cross(error, on_first)
      !! A checkpoint written with the ledger `on_first` is read with it the other way: both must read, run and (ledger on)
      !! start or resume the accumulators sensibly.
      type(error_type), allocatable, intent(out) :: error
      logical, intent(in) :: on_first
      type(ocean_engine_t), target :: ea, eb
      type(config_t) :: cfg_a, cfg_b
      character(len=:), allocatable :: nml_on, nml_off
      real(wp) :: t_a, t_b
      integer :: ierr, step_b
      logical :: ok

      nml_on = replace(uni_nml("basin"), "dt_therm_ratio = 1", "dt_therm_ratio = 2")
      nml_off = replace(nml_on, new_line('a')//"ledger = true", "")
      if (on_first) then
         call make_engine(ea, cfg_a, nml_on, ok)
      else
         call make_engine(ea, cfg_a, nml_off, ok)
      end if
      call check(error, ok, "engine A setup failed")
      if (allocated(error)) return
      t_a = 0.0_wp
      call advance(ea, cfg_a, t_a, 5, ok)
      call check(error, ok, "engine A ran")
      if (.not. allocated(error)) then
         call ocean_state_restart_write(ea%state, ea%grid, ea%decomp, FN_RT, t_a, 5, ierr=ierr)
         call check(error, ierr == OCEAN_STATUS_OK, "checkpoint write failed")
      end if
      call finish_engine(ea)
      if (allocated(error)) return
      if (on_first) then
         call make_engine(eb, cfg_b, nml_off, ok, restart_file=FN_RT, t0=t_b, step0=step_b)
      else
         call make_engine(eb, cfg_b, nml_on, ok, restart_file=FN_RT, t0=t_b, step0=step_b)
      end if
      call check(error, ok .and. eb%warm_restart, "engine B reads the checkpoint")
      if (allocated(error)) then
         call delete_file(FN_RT)
         return
      end if
      call advance(eb, cfg_b, t_b, 4, ok)
      call check(error, ok, "engine B ran 4 steps")
      if (.not. allocated(error) .and. .not. on_first) then
         call check(error, eb%state%icec%ledger%nsteps == 4.0_wp, "off -> on: the accumulators start at zero and count 4 steps")
         if (.not. allocated(error)) call check(error, eb%state%icec%ledger%cum_rel(ICE_LT_WATER) < 1.0e-11_wp .and. &
                                                eb%state%icec%ledger%cum_rel(ICE_LT_SALT) < 1.0e-11_wp .and. &
                                                eb%state%icec%ledger%cum_rel(ICE_LT_HEAT) < 1.0e-6_wp, "and close")
      end if
      call finish_engine(eb)
      call delete_file(FN_RT)
   end subroutine restart_cross

   subroutine test_restart_off_on(error)
      type(error_type), allocatable, intent(out) :: error

      call restart_cross(error, .false.)
   end subroutine test_restart_off_on

   subroutine test_restart_on_off(error)
      type(error_type), allocatable, intent(out) :: error

      call restart_cross(error, .true.)
   end subroutine test_restart_on_off

   subroutine test_restart_ledger(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: ea, eb
      type(config_t) :: cfg_a, cfg_b
      type(snapshot_t) :: sa, sb, sa0, sb0
      real(wp) :: t_a, t_b, cum_a(3), cum_b(3), cres_a(3, ICE_LT_NSTAGE), wa(3), extra_a(3 + ICE_LT_NCLIP + 1), heat_rel_net
      integer :: step_b, ierr, nw, na, n_led, q, s
      logical :: ok, same
      character(len=:), allocatable :: nml

      nw = 7
      na = 6
      ! snow thinner than `hs_min` everywhere, so the thin-snow reset term is non-zero and the restart has something to carry
      ! (its energy is part of the therm1 heat residual: the closure below subtracts the REPORTED term)
      nml = replace(uni_nml("basin"), "dt_therm_ratio = 1", "dt_therm_ratio = 2")
      nml = replace(nml, "h_snow = 0.1", "h_snow = 0.00005")
      call make_engine(ea, cfg_a, nml, ok)
      call check(error, ok, "engine A setup failed")
      if (allocated(error)) return
      t_a = 0.0_wp
      call advance(ea, cfg_a, t_a, nw, ok)
      call check(error, ok, "engine A failed before the checkpoint")
      if (allocated(error)) return
      call ocean_state_restart_write(ea%state, ea%grid, ea%decomp, FN_RT, t_a, nw, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "checkpoint write failed")
      if (allocated(error)) return
      call take_snapshot(ea, sa0)
      call advance(ea, cfg_a, t_a, na, ok)
      call check(error, ok, "engine A failed after the checkpoint")
      if (allocated(error)) return
      call take_snapshot(ea, sa)
      cum_a = [ea%state%icec%ledger%cum_scale(1), ea%state%icec%ledger%cum_scale(2), ea%state%icec%ledger%cum_scale(3)]
      cres_a = ea%state%icec%ledger%cum_res
      wa = ea%state%icec%ledger%worst
      extra_a = [ea%state%icec%ledger%cum_zap, ea%state%icec%ledger%cum_clip, ea%state%icec%ledger%cum_snow_reset]
      call finish_engine(ea)

      call make_engine(eb, cfg_b, nml, ok, restart_file=FN_RT, t0=t_b, step0=step_b)
      call check(error, ok, "engine B (warm restart) setup failed")
      if (allocated(error)) return
      call check(error, step_b == nw .and. eb%warm_restart, "warm restart taken at the written step")
      if (allocated(error)) return
      call take_snapshot(eb, sb0)
      call compare(error, sa0, sb0, "at the resume point", n_led)
      if (allocated(error)) then
         call finish_engine(eb)
         return
      end if
      call check(error, n_led >= 30, "the ledger accumulators are registry entries and were read back")
      if (allocated(error)) then
         call finish_engine(eb)
         return
      end if
      call check(error, eb%state%icec%ledger%nsteps == real(nw, wp), "the step counter was restored")
      call advance(eb, cfg_b, t_b, na, ok)
      call check(error, ok, "engine B failed after the restart")
      if (allocated(error)) then
         call finish_engine(eb)
         return
      end if
      call take_snapshot(eb, sb)
      cum_b = [eb%state%icec%ledger%cum_scale(1), eb%state%icec%ledger%cum_scale(2), eb%state%icec%ledger%cum_scale(3)]
      write (output_unit, "(a,i0,a,2es13.5)") "  ledger entries in the registry: ", n_led, "; cumulative heat scale A, B:", &
         cum_a(3), cum_b(3)
      call check(error, all([(transfer(cum_a(q), 0_int64) == transfer(cum_b(q), 0_int64), q=1, 3)]), &
                 "cumulative scales: N + restart + N == 2N, bitwise")
      same = .true.
      do s = 1, ICE_LT_NSTAGE
         do q = 1, 3
            if (transfer(cres_a(q, s), 0_int64) /= transfer(eb%state%icec%ledger%cum_res(q, s), 0_int64)) same = .false.
         end do
      end do
      if (.not. allocated(error)) call check(error, same, "cumulative residuals: N + restart + N == 2N, bitwise")
      if (.not. allocated(error)) call check(error, all([(transfer(wa(q), 0_int64) == &
                                                          transfer(eb%state%icec%ledger%worst(q), 0_int64), q=1, 3)]), &
                                             "worst per-step residuals (all three lines) bitwise")
      if (.not. allocated(error)) call check(error, all(transfer(extra_a, 0_int64, size(extra_a)) == &
                                                        transfer([eb%state%icec%ledger%cum_zap, eb%state%icec%ledger%cum_clip, &
                                                                  eb%state%icec%ledger%cum_snow_reset], 0_int64, size(extra_a))), &
                                             "cumulative zap, clip record and thin-snow reset: N + restart + N == 2N, bitwise")
      if (.not. allocated(error)) call check(error, any(extra_a /= 0.0_wp), &
                                             "the restored accumulators are not all zero (the check can fail)")
      ! the heat line net of the reported thin-snow reset (an energy the profile initialisation adds in `therm1`)
      heat_rel_net = 0.0_wp
      do s = 1, ICE_LT_NSTAGE
         heat_rel_net = max(heat_rel_net, abs(eb%state%icec%ledger%cum_res(ICE_LT_HEAT, s) &
                                              - merge(eb%state%icec%ledger%cum_snow_reset, 0.0_wp, s == ICE_LT_THERM1)))
      end do
      heat_rel_net = heat_rel_net/max(eb%state%icec%ledger%cum_scale(ICE_LT_HEAT), tiny(1.0_wp))
      if (.not. allocated(error)) call check(error, eb%state%icec%ledger%cum_rel(ICE_LT_WATER) < 1.0e-11_wp .and. &
                                             eb%state%icec%ledger%cum_rel(ICE_LT_SALT) < 1.0e-11_wp .and. &
                                             heat_rel_net < 1.0e-6_wp, &
                                             "the closure holds after the restart")
      call finish_engine(eb)
      call delete_file(FN_RT)
      if (.not. allocated(error)) call compare(error, sa, sb, "after the resumed steps", n_led)
   end subroutine test_restart_ledger

   ! ------------------------------------------------------------------
   ! F: the zaps, the transport clip record, the thin-snow reset
   ! ------------------------------------------------------------------

   subroutine push_cat_state(engine)
      !! Host category state -> device (what a test edited).
      type(ocean_engine_t), intent(inout) :: engine

      !$acc update device(engine%state%ice%cice%aicen, engine%state%ice%cice%vicen, engine%state%ice%cice%vsnon)
      !$acc update device(engine%state%ice%cice%qicen, engine%state%ice%cice%qsnon, engine%state%ice%cice%sicen)
      !$acc update device(engine%state%ice%cice%aice, engine%state%ice%cice%vice, engine%state%ice%cice%vsno)
      !$acc update device(engine%state%ice%cice%aice0)
   end subroutine push_cat_state

   subroutine pull_cat_all(engine)
      !! Device -> host category state including the enthalpies.
      type(ocean_engine_t), intent(inout) :: engine

      call pull_cice(engine)
      !$acc update self(engine%state%ice%cice%qicen, engine%state%ice%cice%qsnon, engine%state%ice%cice%sicen)
   end subroutine pull_cat_all

   subroutine test_zap_booked(error)
      !! Two of the three ITD clean-up zaps, forced in an otherwise sparse ice cover, are the `cleanup` stage's ocean-ward
      !! flux with Icepack's signs and the amounts an independent host expectation removes, and the stage closes.
      !! Step 1: (A) a category of `0.5 puny` of area beside a healthy one, and (C) a near-empty cell (one category
      !! of `0.3 puny`): thermodynamics leaves such categories alone, so the amounts match the host expectation to the frazil
      !! leftovers.  The cover elsewhere is removed so the totals the ledger differences are small enough that the zaps are
      !! far above the sums' round-off.  (The round-off excess `1 < aice < 1 + puny` cannot be forced through the engine: the
      !! stages before the clean-up reduce it; its kernel closure is pinned by `test_cice_itd_cleanup` `aice_excess_scaling`.)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      character(len=:), allocatable :: nml
      real(wp) :: t, puny, expm, exps, exph, hh, zap(3)
      integer :: ng, ia, ja, ib, jb, ic, jc, n, k
      logical :: ok
      real(wp), parameter :: TOL = 1.0e-6_wp

      ! a supercooled ocean: no lateral melt, so the removals of the clean-up are the only thing that takes a thin
      ! category away (new ice is made in category 1: the forced categories are the thick one, 5)
      nml = replace(uni_nml("basin"), "initial_temperature = -1.0", "initial_temperature = -1.84")
      nml = replace(nml, "T_init_surface = -1.0"//new_line('a')//"T_init_bottom = -1.0", "T_init_surface = -1.84"//new_line('a')//"T_init_bottom = -1.84")
      call make_engine(engine, cfg, nml, ok)
      call check(error, ok, "engine setup failed")
      if (allocated(error)) return
      t = 0.0_wp
      call advance(engine, cfg, t, 2, ok)
      call check(error, ok, "2 steps ran")
      if (allocated(error)) return
      call pull_cat_all(engine)
      ng = engine%grid%nghost
      ia = ng + 6
      ja = ng + 6
      ib = ng + 12
      jb = ng + 6
      ic = ng + 18
      jc = ng + 6
      expm = 0.0_wp
      exps = 0.0_wp
      exph = 0.0_wp
      associate (c => engine%state%ice%cice)
         puny = c%p%c%puny
         n = maxloc(c%aicen(ia, ja, :), 1)
         call wipe_cover_except(c, [ia, ib, ic], [ja, jb, jc])
         hh = 5.0_wp
         k = c%ncat
         ! (A) a second category of 0.5 puny beside the healthy one, with the healthy one's enthalpies
         c%aicen(ia, ja, k) = 0.5_wp*puny
         c%vicen(ia, ja, k) = 0.5_wp*puny*hh
         c%vsnon(ia, ja, k) = 0.5_wp*puny*0.1_wp
         c%qicen(ia, ja, k, :) = c%qicen(ia, ja, n, :)
         c%qsnon(ia, ja, k, :) = c%qsnon(ia, ja, n, :)
         call zap_expect(c, ia, ja, k, engine%state%ice_cgrid%tarea(ia, ja), 1.0_wp, expm, exps, exph)
         ! (C) a near-empty cell: ONE category of 0.3 puny (its enthalpies from the template)
         c%aicen(ic, jc, :) = 0.0_wp
         c%vicen(ic, jc, :) = 0.0_wp
         c%vsnon(ic, jc, :) = 0.0_wp
         c%aicen(ic, jc, k) = 0.3_wp*puny
         c%vicen(ic, jc, k) = 0.3_wp*puny*hh
         c%vsnon(ic, jc, k) = 0.3_wp*puny*0.1_wp
         c%qicen(ic, jc, k, :) = c%qicen(ia, ja, n, :)
         c%qsnon(ic, jc, k, :) = c%qsnon(ia, ja, n, :)
         call zap_expect(c, ic, jc, k, engine%state%ice_cgrid%tarea(ic, jc), 1.0_wp, expm, exps, exph)
         call refresh_aggregates(c)
      end associate
      call push_cat_state(engine)
      call one_step(engine, cfg, t, ok)
      call check(error, ok, "the step with the forced zaps ran")
      if (allocated(error)) then
         call finish_engine(engine)
         return
      end if
      associate (L => engine%state%icec%ledger)
         do k = ICE_LT_REMAP, ICE_LT_CLEANUP
            write (output_unit, "(a,a9,3es13.5)") "  ocean-ward flux of stage ", ice_lt_stage_name(k), L%tab%fout(:, k)
         end do
         zap = L%tab%fout(:, ICE_LT_CLEANUP)
         write (output_unit, "(a,3es13.5)") "  independent host expectation of the zaps:                    ", expm, exps, exph
         call check(error, expm > 0.0_wp .and. exps > 0.0_wp .and. exph < 0.0_wp, "the expectation has the flux signs")
         if (.not. allocated(error)) call check(error, L%tab%fout(1, ICE_LT_CLEANUP) > 0.0_wp .and. &
                                                L%tab%fout(2, ICE_LT_CLEANUP) > 0.0_wp .and. &
                                                L%tab%fout(3, ICE_LT_CLEANUP) < 0.0_wp, &
                                                "the zap is non-zero with mass and salt positive, heat negative")
         ! (the booked zap holds the forced removals exactly plus the thin leftovers of the frazil of the supercooled ocean
         ! that the same clean-up removes, a few percent: a lower bound at round-off, an upper bound at 5 %)
         if (.not. allocated(error)) call check(error, zap(1) >= (1.0_wp - TOL)*expm .and. zap(1) <= 1.05_wp*expm, &
                                                "zapped mass is the host expectation (+ the frazil leftovers)")
         if (.not. allocated(error)) call check(error, zap(2) >= (1.0_wp - TOL)*exps .and. zap(2) <= 1.05_wp*exps, &
                                                "zapped salt is the host expectation (+ the frazil leftovers)")
         if (.not. allocated(error)) call check(error, abs(zap(3)) >= (1.0_wp - TOL)*abs(exph) .and. &
                                                abs(zap(3)) <= 1.05_wp*abs(exph), &
                                                "zapped heat is the host expectation (+ the frazil leftovers)")
         ! the stage closes: the state lost what the flux slot handed over (the residual is far below the zap)
         if (.not. allocated(error)) call check(error, abs(L%tab%res(1, ICE_LT_CLEANUP)) <= 1.0e-3_wp*expm, &
                                                "cleanup stage: water closes against the zap")
         if (.not. allocated(error)) call check(error, abs(L%tab%res(2, ICE_LT_CLEANUP)) <= 1.0e-3_wp*exps, &
                                                "cleanup stage: salt closes against the zap")
         if (.not. allocated(error)) call check(error, abs(L%tab%res(3, ICE_LT_CLEANUP)) <= 1.0e-3_wp*abs(exph), &
                                                "cleanup stage: heat closes against the zap")
         if (.not. allocated(error)) call check(error, L%tab%dstate(1, ICE_LT_CLEANUP) < 0.0_wp, &
                                                "the state lost the mass the zap handed over")
         if (.not. allocated(error)) call check(error, L%step_rel(ICE_LT_WATER, 0) < 1.0e-12_wp .and. &
                                                L%step_rel(ICE_LT_SALT, 0) < 1.0e-12_wp, &
                                                "the whole step closes in water and salt")
         if (.not. allocated(error)) call check(error, L%cum_zap(1) >= L%tab%fout(1, ICE_LT_CLEANUP), &
                                                "the cumulative zap term has it")
      end associate
      if (allocated(error)) then
         call finish_engine(engine)
         return
      end if
      call finish_engine(engine)
   end subroutine test_zap_booked

   subroutine wipe_cover_except(c, ik, jk)
      !! Remove the ice everywhere but the cells `(ik, jk)`; keep their state.
      use rdb_ice_cice_state, only: ice_cice_state_t
      type(ice_cice_state_t), intent(inout) :: c
      integer, intent(in) :: ik(:), jk(:)
      integer :: i, j, m
      logical :: keep

      do j = 1, c%ny
         do i = 1, c%nx
            keep = .false.
            do m = 1, size(ik)
               if (i == ik(m) .and. j == jk(m)) keep = .true.
            end do
            if (keep) cycle
            c%aicen(i, j, :) = 0.0_wp
            c%vicen(i, j, :) = 0.0_wp
            c%vsnon(i, j, :) = 0.0_wp
         end do
      end do
   end subroutine wipe_cover_except

   subroutine refresh_aggregates(c)
      !! The host cell aggregates from the category arrays (a test edit does not leave them stale).
      use rdb_ice_cice_state, only: ice_cice_state_t
      type(ice_cice_state_t), intent(inout) :: c

      c%aice = sum(c%aicen, dim=3)
      c%vice = sum(c%vicen, dim=3)
      c%vsno = sum(c%vsnon, dim=3)
      c%aice0 = 1.0_wp - c%aice
   end subroutine refresh_aggregates

   subroutine zap_expect(c, i, j, n, ta, frac, expm, exps, exph)
      !! ADD to the expected zap totals the removal of the fraction `frac` of category `n` of cell `(i, j)`
      !! (mass, salt, enthalpy with the flux slot's signs), from the host arrays alone.
      use rdb_ice_cice_state, only: ice_cice_state_t
      type(ice_cice_state_t), intent(in) :: c
      integer, intent(in) :: i, j, n
      real(wp), intent(in) :: ta, frac
      real(wp), intent(inout) :: expm, exps, exph
      integer :: k

      expm = expm + ta*frac*(c%p%c%rhoi*c%vicen(i, j, n) + c%p%c%rhos*c%vsnon(i, j, n))
      exps = exps + ta*frac*c%p%c%rhoi*c%vicen(i, j, n)*c%p%thermo%ice_ref_salinity*1.0e-3_wp
      do k = 1, c%nilyr
         exph = exph + ta*frac*c%qicen(i, j, n, k)*c%vicen(i, j, n)/real(c%nilyr, wp)
      end do
      do k = 1, c%nslyr
         exph = exph + ta*frac*c%qsnon(i, j, n, k)*c%vsnon(i, j, n)/real(c%nslyr, wp)
      end do
   end subroutine zap_expect

   subroutine test_clip_booked(error)
      !! The remap transport's clip record is a named source of the `transport` stage, and the wiring between the remap's
      !! per-cell record and the ledger is exercised end to end.  A category of one cell is PACKED by the transport's own
      !! `transport_pack` (so the tracers are in the remap's frame), its area flux is set to take `0.5 puny` more than the
      !! cell holds with no tracer flux (so the update's guards drop its whole content), and `remap_ledger_cells_impl` -- the
      !! kernel of the production step -- writes the clip planes; the state loses the category between `begin` and the
      !! `transport` probe; the ledger is handed `fw%clip` itself.  The stage must close to round-off, and the booked terms
      !! must equal the independently computed removal (a factor `nilyr` in the energy frame would fail both).  The same
      !! sequence with the record withheld (the negative control) leaves the whole drop as the stage's residual.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      real(wp) :: t, ta, v, sn, de_i, de_s, rhoi, rhos, lf, sref, w_exp, s_exp, h_exp, res_w(2), res_h(2), puny, a0
      real(wp), allocatable :: tr0(:)
      integer :: ng, i, j, n, k, pass, nqs, nq2
      integer(int32) :: ierr_p
      logical :: ok

      call make_engine(engine, cfg, replace(uni_nml("basin"), 'advection = "upwind"', 'advection = "remap"'), ok)
      call check(error, ok, "engine setup failed")
      if (allocated(error)) return
      t = 0.0_wp
      call advance(engine, cfg, t, 2, ok)
      call check(error, ok, "2 steps ran")
      if (allocated(error)) return
      associate (c => engine%state%ice%cice, L => engine%state%icec%ledger, w => engine%state%icec%work, &
                 tw => engine%state%icec%tws)
         call check(error, allocated(tw%fw%clip), "the clip record is allocated under the ledger with the remap")
         if (allocated(error)) then
            call finish_engine(engine)
            return
         end if
         ng = engine%grid%nghost
         i = ng + 7
         j = ng + 5
         rhoi = c%p%c%rhoi
         rhos = c%p%c%rhos
         lf = c%p%c%Lfresh
         sref = c%p%thermo%ice_ref_salinity
         puny = c%p%c%puny
         nq2 = tw%a%tab%nt_qsno + 1
         nqs = tw%a%tab%nt_qsno + 2
         do pass = 1, 2
            call pull_cat_all(engine)
            n = maxloc(c%aicen(i, j, :), 1)
            ta = engine%state%ice_cgrid%tarea(i, j)
            ! the packed state by the transport's own routine, then ONLY this category of this cell kept
            call transport_pack(tw%a, c%p%c, c%nilyr, c%nslyr, 1, c%nx, 1, c%ny, c%aice0, c%aicen, c%vicen, c%vsnon, &
                                c%Tsfcn, c%qicen, c%sicen, c%qsnon, ierr_p)
            !$acc update self(tw%a%aim, tw%a%trm)
            a0 = tw%a%aim(i, j, n)
            tr0 = tw%a%trm(i, j, :, n)
            tw%a%aim = 0.0_wp
            tw%a%trm = 0.0_wp
            tw%a%aim(i, j, n) = a0
            tw%a%trm(i, j, :, n) = tr0
            tw%fw%mflx_e = 0.0_wp
            tw%fw%mflx_n = 0.0_wp
            tw%fw%mtflx_e = 0.0_wp
            tw%fw%mtflx_n = 0.0_wp
            tw%fw%mflx_e(i, j, 1) = (tw%a%aim(i, j, n) + 0.5_wp*puny)*ta
            tw%fw%clip = 0.0_wp
            !$acc update device(tw%a%aim, tw%a%trm, tw%fw%mflx_e, tw%fw%mflx_n, tw%fw%mtflx_e, tw%fw%mtflx_n, tw%fw%clip)
            call remap_ledger_cells_impl(c%nx, c%ny, c%ncat, tw%a%tab%ntrace, tw%fw%nbatch, ng + 1, ng + engine%grid%nx_phys, &
                                         ng + 1, ng + engine%grid%ny_phys, n, n, puny, tw%a%tab%tracer_type(1:tw%a%tab%ntrace), &
                                         tw%a%tab%depend(1:tw%a%tab%ntrace), 4, nq2, nqs, engine%state%ice_cgrid%tarea, &
                                         tw%a%aim, tw%a%trm, tw%fw%mflx_e, tw%fw%mflx_n, tw%fw%mtflx_e, tw%fw%mtflx_n, tw%fw%clip)
            call ice_flux_zero(engine%state%icec%flux)     ! (the sequencer zeroes the accumulators before `begin`)
            call L%begin(c, DT, w%fsnow, w%frain, w%flw, .false.)
            ! the removal the remap's guards stand for: the whole category (ice, snow; the enthalpies stay behind)
            v = c%vicen(i, j, n)
            sn = c%vsnon(i, j, n)
            de_i = 0.0_wp
            do k = 1, c%nilyr
               de_i = de_i + c%qicen(i, j, n, k)*v/real(c%nilyr, wp)
            end do
            de_s = 0.0_wp
            do k = 1, c%nslyr
               de_s = de_s + c%qsnon(i, j, n, k)*sn/real(c%nslyr, wp)
            end do
            c%aicen(i, j, n) = 0.0_wp
            c%vicen(i, j, n) = 0.0_wp
            c%vsnon(i, j, n) = 0.0_wp
            call refresh_aggregates(c)
            call push_cat_state(engine)
            if (pass == 1) call L%book_clip(tw%fw%clip)
            call L%probe(ICE_LT_TRANSPORT, c, engine%state%icec%flux)
            call L%end_step()
            w_exp = -(rhoi*v + rhos*sn)*ta
            s_exp = -rhoi*v*sref*1.0e-3_wp*ta
            h_exp = -(de_i + de_s)*ta
            res_w(pass) = L%tab%res(1, ICE_LT_TRANSPORT)
            res_h(pass) = L%tab%res(3, ICE_LT_TRANSPORT)
            if (pass == 1) then
            write (output_unit, "(a,3es13.5)") "  booked clip source (water, salt, heat):", L%tab%clip_w, L%tab%clip_s, L%tab%clip_h
            write (output_unit, "(a,3es13.5)") "  state change in the stage (water, salt, heat):", L%tab%dstate(:, ICE_LT_TRANSPORT)
               call check(error, L%tab%n_clip == 1.0_wp, "one clipped cell is reported")
               if (.not. allocated(error)) call check(error, abs(L%tab%clip_w - w_exp) <= 1.0e-12_wp*abs(w_exp), &
                                                      "booked water equals the removal")
               if (.not. allocated(error)) call check(error, abs(L%tab%clip_s - s_exp) <= 1.0e-12_wp*abs(s_exp), &
                                                      "booked salt equals the removal")
               if (.not. allocated(error)) call check(error, abs(L%tab%clip_h - h_exp) <= 1.0e-9_wp*abs(h_exp), &
                                                      "booked heat equals the removal (the remap's energy frame is the ledger's)")
               if (.not. allocated(error)) call check(error, abs(L%tab%dstate(1, ICE_LT_TRANSPORT) - w_exp) <= &
                                                      1.0e-6_wp*abs(w_exp), "the state actually lost that mass")
               if (.not. allocated(error)) call check(error, abs(L%tab%res(1, ICE_LT_TRANSPORT)) <= 1.0e-6_wp*abs(w_exp) .and. &
                                                      abs(L%tab%res(3, ICE_LT_TRANSPORT)) <= 1.0e-6_wp*abs(h_exp), &
                                                      "the transport stage closes with the clip record booked")
            end if
         end do
         ! the negative control: the same removal with the record withheld is the stage's whole residual
         if (.not. allocated(error)) then
   write (output_unit, "(a,2es13.5,a,es13.5)") "  transport-stage residual (water, heat) booked vs withheld:", res_w(1), res_h(1), &
               " vs ", res_w(2)
            call check(error, abs(res_w(2)) > 0.5_wp*abs(w_exp) .and. abs(res_h(2)) > 0.5_wp*abs(h_exp), &
                       "withheld, the drop shows as the residual (the test can fail)")
         end if
      end associate
      call finish_engine(engine)
   end subroutine test_clip_booked

   subroutine test_snow_reset_reported(error)
      !! The thin-snow enthalpy reset is reported as its own diagnostic term.  A category's snow is made thinner than
      !! `hs_min` with an enthalpy colder than the melting reference in two cells; the term equals the independent host
      !! sum `A (-rhos Lfresh - q) v/nslyr`, is part of (not subtracted from) the `therm1` heat residual, and
      !! accounts for it to within the solver's bound.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      real(wp) :: t, ta, expect, hs, rhos, lf, resid, bound
      integer :: ng, ii(4), jj(4), m, n, k
      logical :: ok

      call make_engine(engine, cfg, uni_nml("basin"), ok)
      call check(error, ok, "engine setup failed")
      if (allocated(error)) return
      t = 0.0_wp
      call advance(engine, cfg, t, 2, ok)
      call check(error, ok, "2 steps ran")
      if (allocated(error)) return
      call pull_cat_all(engine)
      associate (c => engine%state%ice%cice)
         ng = engine%grid%nghost
         ii = [ng + 5, ng + 15, ng + 9, ng + 19]
         jj = [ng + 5, ng + 9, ng + 12, ng + 4]
         rhos = c%p%c%rhos
         lf = c%p%c%Lfresh
         expect = 0.0_wp
         do m = 1, 4
            n = maxloc(c%aicen(ii(m), jj(m), :), 1)
            hs = 0.5_wp*c%p%c%hs_min
            c%vsnon(ii(m), jj(m), n) = hs*c%aicen(ii(m), jj(m), n)
            c%qsnon(ii(m), jj(m), n, :) = -1.5_wp*rhos*lf
            ta = engine%state%ice_cgrid%tarea(ii(m), jj(m))
            do k = 1, c%nslyr
               expect = expect + ta*((-rhos*lf) + 1.5_wp*rhos*lf)*c%vsnon(ii(m), jj(m), n)/real(c%nslyr, wp)
            end do
         end do
         call refresh_aggregates(c)
      end associate
      call push_cat_state(engine)
      call one_step(engine, cfg, t, ok)
      call check(error, ok, "the step ran")
      if (allocated(error)) then
         call finish_engine(engine)
         return
      end if
      associate (L => engine%state%icec%ledger)
         resid = L%tab%res(3, ICE_LT_THERM1)
         bound = L%tab%heat_bound
         write (output_unit, "(a,3es13.5)") "  thin-snow reset (reported, host expectation) and therm1 heat residual:", &
            L%tab%snow_reset, expect, resid
         call check(error, expect > 0.0_wp, "the forced reset adds energy")
         if (.not. allocated(error)) call check(error, abs(L%tab%snow_reset - expect) <= 1.0e-9_wp*expect, &
                                                "the reported reset equals the independent host value")
         if (.not. allocated(error)) call check(error, abs(resid - L%tab%snow_reset) <= bound, &
                                                "the reset accounts for the therm1 heat residual to the solver bound")
         if (.not. allocated(error)) call check(error, L%tab%snow_reset > 20.0_wp*bound, &
                                                "the reset is far above the solver bound (the attribution is not vacuous)")
         if (.not. allocated(error)) call check(error, L%cum_snow_reset >= L%tab%snow_reset, "the cumulative term has it")
      end associate
      call finish_engine(engine)
   end subroutine test_snow_reset_reported

end module test_cice_closure
