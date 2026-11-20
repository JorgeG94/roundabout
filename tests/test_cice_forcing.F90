!! The CICE-port sea ice (`[ocean.ice] scheme = "cice"`) under a PRESCRIBED, TIME-VARYING atmosphere
!! (`[ocean.atm_forcing] enable` + the bulk forcing), through the production engine path (D11.4b).
!! The atmosphere is a small synthetic NetCDF set written in the test: nine fields, each affine in
!! (lon, lat, time) so the bilinear regrid and the time blend reproduce the closed form exactly.
!!
!! Cases:
!!   * `sis2_inert`            the default scheme allocates, maps and registers nothing of the new wiring, and the
!!                             stress it writes is the weighted one of the bulk forcing (optionally dumps two restart
!!                             files for the cross-build byte-identity check, `RDB_CICE_FORCING_DUMP`).
!!   * `refusals`              the uniform and the prescribed atmosphere together are refused; the prescribed
!!                             atmosphere alone is accepted; a file wind stress stays refused.
!!   * `unweighted_stress`     the bulk forcing exposes the unweighted wind stress of the step: it equals the
!!                             closed-form bulk formula, and it moves from one step to the next.
!!   * `stress_blend_current_wind` ice-free faces carry the wind stress unchanged, ice-covered faces
!!                             `(1 - a) tau + ice-ocean stress` with the wind of the CURRENT step, at the top of the
!!                             step and after the ice step; the configure-time snapshot is not used.
!!   * `consumed_stress_independent` the stress the ocean dynamics consumes against an independent host build.
!!   * `face_fraction_current` the wind share is 1 - the face mean of the CURRENT concentration (host check).
!!   * `ice_inputs_from_forcing` the ice's atmosphere inputs are the forcing's of this step (state, rotation, units).
!!   * `flux_partition`        every surface flux of a partly ice-covered cell: the open-water share (weight `aice0`)
!!                             and the ice share through the hand-over add up to the cell.
!!   * `closure_column`, `closure_basin`, `closure_basin_ratio2` the ice ledger and the ocean budgets with a
!!                             time-varying atmosphere.
!!   * `restart_tags_required` the ice-ocean stress tags are required under the prescribed atmosphere.
!!   * `restart_time_varying`  N + checkpoint (mid thermo window) + N == 2N, every registry field bitwise.
module test_cice_forcing
   use, intrinsic :: iso_fortran_env, only: int32, int64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
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
   use rdb_ice_errors, only: ICE_ERR_FATAL_MASK
   use rdb_ice_config, only: ice_config_check, ICE_CONFIG_MSG_LEN, ICE_CONFIG_MAX_MSGS
   use rdb_ice_ledger, only: ice_ledger_mirror_t, ICE_LEDGER_WATER, ICE_LEDGER_SALT, ICE_LEDGER_HEAT
   use rdb_ocean_console_stats, only: compute_total_tracer
   use rdb_ocean_bulk_flux, only: bulk_flux_column, bulk_flux_longwave_up, BULK_T0_KELVIN
   use test_cice_synth_atm, only: write_synth_atm, remove_synth_atm
   implicit none
   private

   public :: collect_cice_forcing_tests

   real(wp), parameter :: DT = 600.0_wp
   character(len=*), parameter :: FN_ATM = "test_cice_forcing_atm.nc"
   character(len=*), parameter :: FN_LST = "test_cice_forcing_atm.txt"
   character(len=*), parameter :: FN_RT = "test_cice_forcing_rt.nc"
   integer, parameter :: MAXF = 300
   real(wp), parameter :: LON_W = 20.0_wp, LAT_S = 70.0_wp, DLL = 0.25_wp
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
      real(wp) :: water = 0.0_wp, salt = 0.0_wp, heat = 0.0_wp
      real(wp) :: ocean_heat = 0.0_wp, ocean_salt = 0.0_wp
      real(wp) :: assembled_heat = 0.0_wp, assembled_salt = 0.0_wp, ice_share = 0.0_wp
      real(wp) :: vice_first = 0.0_wp, vice_last = 0.0_wp, aice_max = 0.0_wp
      logical :: finite = .true.
   end type closure_t

contains

   subroutine collect_cice_forcing_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("sis2_inert", test_sis2_inert), &
                  new_unittest("refusals", test_refusals), &
                  new_unittest("unweighted_stress", test_unweighted_stress), &
                  new_unittest("stress_blend_current_wind", test_stress_blend), &
                  new_unittest("consumed_stress_independent", test_consumed_stress), &
                  new_unittest("face_fraction_current", test_face_fraction_current), &
                  new_unittest("ice_inputs_from_forcing", test_ice_inputs), &
                  new_unittest("flux_partition", test_flux_partition), &
                  new_unittest("closure_column", test_closure_column), &
                  new_unittest("closure_basin", test_closure_basin), &
                  new_unittest("closure_basin_ratio2", test_closure_basin_ratio2), &
                  new_unittest("restart_time_varying", test_restart), &
                  new_unittest("restart_tags_required", test_restart_tags_required) &
                  ]
   end subroutine collect_cice_forcing_tests

   ! ------------------------------------------------------------------
   ! the synthetic atmosphere
   ! ------------------------------------------------------------------

   subroutine write_atm_files(scn)
      character(len=*), intent(in) :: scn
      call write_synth_atm(scn, FN_ATM, FN_LST)
   end subroutine write_atm_files

   subroutine remove_atm_files()
      call remove_synth_atm(FN_ATM, FN_LST)
   end subroutine remove_atm_files

   ! ------------------------------------------------------------------
   ! inputs and engine plumbing
   ! ------------------------------------------------------------------

   function case_nml(geom, conc, ratio, prescribed, extra) result(nml)
      !! A closed spherical high-latitude basin (`column` 8x8, `basin` 12x10 cells of 0.25 degree), the
      !! CICE-port ice under it.  `conc`: the initial concentration (0 = no ice).  `prescribed`: the
      !! prescribed atmosphere with the bulk forcing, else the uniform atmosphere.
      character(len=*), intent(in) :: geom
      real(wp), intent(in) :: conc
      integer, intent(in) :: ratio
      logical, intent(in) :: prescribed
      character(len=*), intent(in), optional :: extra
      character(len=:), allocatable :: nml
      character(len=*), parameter :: NL = new_line("a")
      character(len=:), allocatable :: grid, ic, atm
      character(len=32) :: sconc, sratio, snx, sny

      select case (geom)
      case ("column")
         snx = "8"
         sny = "8"
      case ("basin")
         snx = "12"
         sny = "10"
      case default
         error stop "test_cice_forcing: unknown geometry"
      end select
      write (sconc, '(f6.3)') conc
      write (sratio, '(i0)') ratio
      grid = "[grid]"//NL// &
             "nx = "//trim(snx)//NL// &
             "ny = "//trim(sny)//NL// &
             "nghost = 3"//NL// &
             "dx = 0.25"//NL// &
             "dy = 0.25"//NL// &
             "[ocean.grid]"//NL// &
             'grid_config = "spherical"'//NL// &
             "lon_west = 20.0"//NL// &
             "lat_south = 70.0"//NL// &
             "rad_earth = 6.371e6"//NL// &
             'coriolis_scheme = "planetary"'//NL
      if (conc > 0.0_wp) then
         ic = "[ocean.ice_ic]"//NL// &
              'conc_config = "uniform"'//NL// &
              "h_ice = 1.0"//NL// &
              "conc = "//trim(adjustl(sconc))//NL// &
              "h_snow = 0.1"//NL// &
              "t_ice = -8.0"//NL
      else
         ic = "[ocean.ice_ic]"//NL// &
              'conc_config = "zero"'//NL
      end if
      if (prescribed) then
         atm = "[ocean.atm_forcing]"//NL// &
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
               "enable = true"//NL
      else
         atm = "[ocean.ice_atm]"//NL// &
               "uni_tas = 255.0"//NL// &
               "uni_rlds = 170.0"//NL// &
               "uni_uas = 3.0"//NL// &
               "uni_vas = 0.0"//NL
      end if
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
            grid//NL// &
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
            "[ocean.ice_thermo]"//NL// &
            "update_ocn_f = true"//NL// &
            "[ocean.ice_dyn]"//NL// &
            'advection = "upwind"'//NL// &
            "ndte = 40"//NL// &
            ic//NL// &
            atm//NL// &
            "[ocean.bc]"//NL// &
            'west = "wall"'//NL// &
            'east = "wall"'//NL// &
            'south = "wall"'//NL// &
            'north = "wall"'//NL// &
            "[output]"//NL// &
            "output_to_file = false"//NL
      if (present(extra)) nml = nml//extra
   end function case_nml

   function sis2_nml(with_ice) result(nml)
      !! The default scheme under the same prescribed atmosphere (bulk forcing on).
      logical, intent(in) :: with_ice
      character(len=:), allocatable :: nml
      character(len=*), parameter :: NL = new_line("a")

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
            "[ocean.atm_forcing]"//NL// &
            "enable = true"//NL// &
            'sw_interp = "linear"'//NL// &
            'fluxes_regrid = "bilinear"'//NL// &
            'time_mode = "interannual"'//NL// &
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
      if (with_ice) then
         nml = nml//NL// &
               "[ocean.ice]"//NL// &
               "enable = true"//NL// &
               "dynamics = true"//NL// &
               "transport = true"//NL// &
               "[ocean.ice_ic]"//NL// &
               'conc_config = "uniform"'//NL// &
               "h_ice = 1.0"//NL// &
               "conc = 0.6"//NL// &
               "h_snow = 0.1"//NL// &
               "t_ice = -8.0"//NL
      end if
   end function sis2_nml

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

   subroutine pull_stress(engine)
      !! Device -> host of the stress hand-off fields.
      type(ocean_engine_t), intent(inout) :: engine
      !$acc update self(engine%state%surface_stress%tau_x, engine%state%surface_stress%tau_y)
      !$acc update self(engine%state%ice%tau_a_x, engine%state%ice%tau_a_y)
      !$acc update self(engine%state%ice%dync%w%aiE, engine%state%ice%dync%w%aiN)
      !$acc update self(engine%state%ice%dync%w%strocnxE, engine%state%ice%dync%w%strocnyN)
      !$acc update self(engine%state%ice%dync%strairxE, engine%state%ice%dync%strairyN)
   end subroutine pull_stress

   subroutine pull_atm(engine)
      type(ocean_engine_t), intent(inout) :: engine
      !$acc update self(engine%state%atm%tas, engine%state%atm%huss, engine%state%atm%uas, engine%state%atm%vas)
      !$acc update self(engine%state%atm%psl, engine%state%atm%rlds, engine%state%atm%rsds)
      !$acc update self(engine%state%atm%prra, engine%state%atm%prsn)
   end subroutine pull_atm

   subroutine pull_cice(engine)
      type(ocean_engine_t), intent(inout) :: engine
      !$acc update self(engine%state%ice%cice%aicen, engine%state%ice%cice%vicen, engine%state%ice%cice%vsnon)
      !$acc update self(engine%state%ice%cice%aice, engine%state%ice%cice%vice, engine%state%ice%cice%vsno)
      !$acc update self(engine%state%ice%cice%aice0, engine%state%ice%cice%aice_init, engine%state%ice%cice%ice_err)
      !$acc update self(engine%state%ice%aice, engine%state%ice%vice, engine%state%ice%vsno, engine%state%ice%aice0)
   end subroutine pull_cice

   subroutine pull_sf(engine)
      type(ocean_engine_t), intent(inout) :: engine
      associate (sf => engine%state%surface_flux)
         !$acc update self(sf%q_sw, sf%q_lw, sf%q_sens, sf%q_lat, sf%evap, sf%lprec, sf%fprec)
         !$acc update self(sf%heat_content_lprec, sf%heat_content_fprec, sf%heat_content_massin, sf%heat_content_massout)
         !$acc update self(sf%heat_ice, sf%sw_ice, sf%salt_ice, sf%Q_heat, sf%Q_salt)
      end associate
   end subroutine pull_sf

   subroutine pull_ice_flux(engine)
      type(ocean_engine_t), intent(inout) :: engine
      !$acc update self(engine%state%icec%flux%fresh, engine%state%icec%flux%fsalt)
      !$acc update self(engine%state%icec%flux%fhocn, engine%state%icec%flux%fswthru)
      !$acc update self(engine%state%icec%work%frain, engine%state%icec%work%fsnow, engine%state%icec%work%flw)
      !$acc update self(engine%state%icec%work%uatm, engine%state%icec%work%vatm, engine%state%icec%work%sw_down)
      !$acc update self(engine%state%icec%work%potT, engine%state%icec%work%rhoa, engine%state%icec%work%Qa)
      !$acc update self(engine%state%icec%oin%sss_ocn)
      !$acc update self(engine%state%ice%frazil_heat, engine%state%ice%frzmlt)
   end subroutine pull_ice_flux

   function rho0_cp(engine) result(r)
      type(ocean_engine_t), intent(in) :: engine
      real(wp) :: r
      r = engine%state%eos%rho0*engine%state%eos%c_p
   end function rho0_cp

   subroutine ocean_totals(engine, heat, salt)
      type(ocean_engine_t), intent(inout) :: engine
      real(wp), intent(out) :: heat, salt
      associate (ms => engine%state%multilayer)
         heat = rho0_cp(engine)*compute_total_tracer(ms%tracers(ms%idx_temperature)%c, ms%h_layer, &
                                                     engine%state%metrics%areaT, engine%grid%nghost)
         salt = engine%state%eos%rho0*compute_total_tracer(ms%tracers(ms%idx_salinity)%c, ms%h_layer, &
                                                           engine%state%metrics%areaT, engine%grid%nghost)
      end associate
   end subroutine ocean_totals

   ! ------------------------------------------------------------------
   ! default scheme inert / refusals
   ! ------------------------------------------------------------------

   subroutine test_sis2_inert(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      character(len=512) :: dumpdir
      character(len=:), allocatable :: nml
      integer :: with, ierr, ilen, nsteps
      real(wp) :: t
      logical :: ok

      nsteps = 6
      call get_environment_variable("RDB_CICE_FORCING_DUMP", dumpdir, length=ilen)
      call write_atm_files("cold")
      do with = 0, 1
         nml = sis2_nml(with == 1)
         call make_engine(engine, cfg, nml, ok)
         call check(error, ok, "sis2 engine setup failed")
         if (allocated(error)) then
            call remove_atm_files()
            return
         end if
         call check(error,.not. engine%state%bulk_forcing%cice_wind, "no CICE wind under sis2")
         call check(error,.not. engine%state%icec%is_init, "no CICE run slot under sis2")
         call check(error,.not. engine%state%ice%cice%is_init, "no CICE state under sis2")
         if (.not. allocated(error)) then
            t = 0.0_wp
            call advance(engine, cfg, t, nsteps, ok)
            call check(error, ok, "sis2 steps ran")
         end if
         if (.not. allocated(error)) then
            if (ilen > 0) then
               call ocean_state_restart_write(engine%state, engine%grid, engine%decomp, &
                                              trim(dumpdir)//"/sis2_bulk_"//trim(merge("ice  ", "noice", with == 1))//".nc", &
                                              t, nsteps, ierr=ierr)
               call check(error, ierr == OCEAN_STATUS_OK, "dump write")
            end if
         end if
         call finish_engine(engine)
         if (allocated(error)) then
            call remove_atm_files()
            return
         end if
      end do
      call remove_atm_files()
   end subroutine test_sis2_inert

   function mentions(msgs, n, key) result(yes)
      character(len=*), intent(in) :: msgs(:), key
      integer, intent(in) :: n
      logical :: yes
      integer :: k

      yes = .false.
      do k = 1, n
         if (index(msgs(k), key) > 0) yes = .true.
      end do
   end function mentions

   subroutine test_refusals(error)
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      character(len=ICE_CONFIG_MSG_LEN) :: msgs(ICE_CONFIG_MAX_MSGS)
      integer :: n, ierr
      character(len=*), parameter :: NL = new_line("a")

      call ensure_comm()
      ! prescribed atmosphere alone: accepted by the ice config check
      call read_config_from_toml_string(case_nml("column", 0.6_wp, 1, .true.), cfg, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "the prescribed-atmosphere input parses")
      if (allocated(error)) return
      call ice_config_check(cfg, msgs, n)
      call check(error, n == 0, "the prescribed atmosphere alone has no refusal: "//trim(msgs(max(n, 1))))
      if (allocated(error)) return
      ! both atmospheres: refused
      call read_config_from_toml_string(case_nml("column", 0.6_wp, 1, .true., &
                                                 "[ocean.ice_atm]"//NL// &
                                                 "uni_tas = 250.0"//NL), cfg, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "the doubled-atmosphere input parses")
      if (allocated(error)) return
      call ice_config_check(cfg, msgs, n)
      call check(error, n > 0 .and. mentions(msgs, n, "is ambiguous"), "uni_* with a prescribed atmosphere is refused")
      if (allocated(error)) return
      call read_config_from_toml_string(replace(case_nml("column", 0.6_wp, 1, .true.), "hfreeze = 10.0", &
                                                "hfreeze = 10.0"//new_line('a')//"snowfall = 1.0e-5"), cfg, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "the snowfall input parses")
      if (allocated(error)) return
      call ice_config_check(cfg, msgs, n)
      call check(error, n > 0 .and. mentions(msgs, n, "is ambiguous"), "snowfall with a prescribed atmosphere is refused")
      if (allocated(error)) return
      ! a file wind stress stays refused
      call read_config_from_toml_string(case_nml("column", 0.6_wp, 1, .false., &
                                                 "[ocean.dataovr]"//NL// &
                                                 'tau_x_file = "wind.nc"'//NL), cfg, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "the file-wind input parses")
      if (allocated(error)) return
      call ice_config_check(cfg, msgs, n)
      call check(error, n > 0 .and. mentions(msgs, n, "tau_x/tau_y files"), "a file wind stress stays refused")
   end subroutine test_refusals

   ! ------------------------------------------------------------------
   ! stress
   ! ------------------------------------------------------------------

   subroutine expected_wind_stress(engine, cfg, i, j, taux, tauy)
      !! The bulk formula at physical cell `(i, j)` from the pulled atmosphere and surface temperature
      !! (`wind_relative = .false.`).
      type(ocean_engine_t), intent(inout), target :: engine
      type(config_t), intent(in) :: cfg
      integer, intent(in) :: i, j
      real(wp), intent(out) :: taux, tauy
      real(wp) :: ca, sa, ug, vg, sens, lat, evp, cd, ch, ce, ustar, sst
      integer :: kt

      associate (a => engine%state%atm, ms => engine%state%multilayer, m => engine%state%metrics)
         ca = cos(m%angle_dx(i, j))
         sa = sin(m%angle_dx(i, j))
         ug = ca*a%uas(i, j) + sa*a%vas(i, j)
         vg = -sa*a%uas(i, j) + ca*a%vas(i, j)
         kt = ms%nz_ml
         sst = ms%tracers(ms%idx_temperature)%c(i, j, kt)
         call bulk_flux_column(ug, vg, 0.0_wp, 0.0_wp, a%tas(i, j), a%huss(i, j), a%psl(i, j), sst, &
                               cfg%ocean%bulk_flux%z_wind, cfg%ocean%bulk_flux%z_ta, cfg%ocean%bulk_flux%n_iter, &
                               taux, tauy, sens, lat, evp, cd, ch, ce, ustar)
      end associate
   end subroutine expected_wind_stress

   subroutine test_unweighted_stress(error)
      !! Under the new scheme the bulk forcing keeps the UNWEIGHTED stress: at partly ice-covered faces it is the
      !! closed-form bulk stress (not `aice0` times it), and it moves from step to step with the wind.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      real(wp), allocatable :: prev(:, :)
      real(wp) :: t, tx, ty, e, worst, moved
      integer :: n, i, j, ng, is, ie, js, je, ierr
      logical :: ok

      call write_atm_files("cold")
      call make_engine(engine, cfg, case_nml("basin", 0.6_wp, 1, .true.), ok)
      call check(error, ok, "engine setup failed")
      if (allocated(error)) then
         call remove_atm_files()
         return
      end if
      call check(error, engine%state%bulk_forcing%cice_wind, "the CICE wind path is on")
      if (.not. allocated(error)) then
         ng = engine%grid%nghost
         is = ng + 1
         ie = ng + engine%grid%nx_phys
         js = ng + 1
         je = ng + engine%grid%ny_phys
         allocate (prev(size(engine%state%ice%tau_a_x, 1), size(engine%state%ice%tau_a_x, 2)))
         prev = 0.0_wp
         t = 0.0_wp
         do n = 1, 4
            call engine_step(engine, DT, t, ierr=ierr)
            call check(error, ierr == OCEAN_STATUS_OK, "engine_step")
            if (allocated(error)) exit
            call pull_stress(engine)
            call pull_atm(engine)
            !$acc update self(engine%state%multilayer%tracers(engine%state%multilayer%idx_temperature)%c)
            !$acc update self(engine%state%icec%flux%fresh)
            worst = 0.0_wp
            do j = js, je
               do i = is + 1, ie
                  ! the x face between cells i-1 and i: the mean of the two centred stresses
                  call expected_wind_stress(engine, cfg, i - 1, j, tx, ty)
                  e = tx
                  call expected_wind_stress(engine, cfg, i, j, tx, ty)
                  e = 0.5_wp*(e + tx)
                  if (.not. ieee_is_finite(engine%state%ice%tau_a_x(i, j))) then
                     call check(error, .false., "tauw_x is finite")
                     exit
                  end if
                  worst = max(worst, abs(engine%state%ice%tau_a_x(i, j) - e)/max(abs(e), 1.0e-6_wp))
               end do
            end do
            if (allocated(error)) exit
            write (*, '(a,i0,a,es10.3)') "  step ", n, ": unweighted tau_x vs the closed form, worst relative ", worst
            ! the exposure is the stress a run without ice would carry: nothing multiplied it by the open fraction
            call check(error, worst < 1.0e-4_wp, "the exposed stress is the unweighted bulk stress (the surface temperature drifts 1e-4 K a step)")
            if (allocated(error)) exit
            if (n > 1) then
               moved = maxval(abs(engine%state%ice%tau_a_x(is + 1:ie, js:je) - prev(is + 1:ie, js:je)))
               call check(error, moved > 1.0e-5_wp, "the exposed stress follows the wind from step to step")
               if (allocated(error)) exit
            end if
            prev = engine%state%ice%tau_a_x
            call engine_step_ice(engine, cfg, DT, t, ierr=ierr)
            call check(error, ierr == OCEAN_STATUS_OK, "engine_step_ice")
            if (allocated(error)) exit
            call engine_step_finalize(engine, DT, t, ierr=ierr)
            call check(error, ierr == OCEAN_STATUS_OK, "engine_step_finalize")
            if (allocated(error)) exit
            t = t + DT
         end do
      end if
      call finish_engine(engine)
      call remove_atm_files()
   end subroutine test_unweighted_stress

   subroutine face_blend_x(engine, i, j, tauw, tau)
      !! The face blend at x face `i`, row `j`: `((1 - aiE) tauw - strocnxE) * wet_u` with the ice fields of the
      !! pulled state; `i - 1` is the `E` index of the face.
      type(ocean_engine_t), intent(in) :: engine
      integer, intent(in) :: i, j
      real(wp), intent(in) :: tauw
      real(wp), intent(out) :: tau

      tau = ((1.0_wp - engine%state%ice%dync%w%aiE(i - 1, j))*tauw - engine%state%ice%dync%w%strocnxE(i - 1, j)) &
            *engine%state%metrics%wet_u(i, j)
   end subroutine face_blend_x

   subroutine face_blend_y(engine, i, j, tauw, tau)
      type(ocean_engine_t), intent(in) :: engine
      integer, intent(in) :: i, j
      real(wp), intent(in) :: tauw
      real(wp), intent(out) :: tau

      tau = ((1.0_wp - engine%state%ice%dync%w%aiN(i, j - 1))*tauw - engine%state%ice%dync%w%strocnyN(i, j - 1)) &
            *engine%state%metrics%wet_v(i, j)
   end subroutine face_blend_y

   subroutine blend_residual(engine, wind_x, wind_y, worst_x, worst_y, scale)
      !! Worst `|ss%tau - blend(wind)|` over the interior faces, and the largest stress.
      type(ocean_engine_t), intent(in) :: engine
      real(wp), intent(in) :: wind_x(:, :), wind_y(:, :)
      real(wp), intent(out) :: worst_x, worst_y, scale
      real(wp) :: tau
      integer :: i, j, ng, is, ie, js, je

      ng = engine%grid%nghost
      is = ng + 1
      ie = ng + engine%grid%nx_phys
      js = ng + 1
      je = ng + engine%grid%ny_phys
      worst_x = 0.0_wp
      worst_y = 0.0_wp
      scale = 0.0_wp
      do j = js, je
         do i = is + 1, ie
            call face_blend_x(engine, i, j, wind_x(i, j), tau)
            worst_x = max(worst_x, abs(engine%state%surface_stress%tau_x(i, j) - tau))
            scale = max(scale, abs(tau))
         end do
      end do
      do j = js + 1, je
         do i = is, ie
            call face_blend_y(engine, i, j, wind_y(i, j), tau)
            worst_y = max(worst_y, abs(engine%state%surface_stress%tau_y(i, j) - tau))
            scale = max(scale, abs(tau))
         end do
      end do
   end subroutine blend_residual

   subroutine test_stress_blend(error)
      !! (a) A time-varying wind reaches the ocean through the blend with the CURRENT step's wind: ice-free
      !! faces unchanged, ice-covered faces `(1 - a) tau + ice-ocean stress`.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: engine, ef
      type(config_t) :: cfg, cff
      real(wp), allocatable :: wind_prev_x(:, :), wind_prev_y(:, :)
      real(wp) :: t, tf, wx, wy, sc, wx_old, wy_old, sc_old, dmax, amax
      integer :: n, ierr, ng, is, ie, js, je, i, j
      logical :: ok

      call write_atm_files("cold")
      call make_engine(engine, cfg, case_nml("basin", 0.6_wp, 1, .true.), ok)
      call check(error, ok, "covered engine setup failed")
      if (allocated(error)) then
         call remove_atm_files()
         return
      end if
      ng = engine%grid%nghost
      is = ng + 1
      ie = ng + engine%grid%nx_phys
      js = ng + 1
      je = ng + engine%grid%ny_phys
      t = 0.0_wp
      do n = 1, 6
         call engine_step(engine, DT, t, ierr=ierr)
         call check(error, ierr == OCEAN_STATUS_OK, "engine_step")
         if (allocated(error)) exit
         call pull_stress(engine)
         ! top of the step: the blend of the wind of THIS step with the ice state the last ice step left
         call blend_residual(engine, engine%state%ice%tau_a_x, engine%state%ice%tau_a_y, wx, wy, sc)
         if (n > 1) then
            call blend_residual(engine, wind_prev_x, wind_prev_y, wx_old, wy_old, sc_old)
            write (*, '(a,i0,a,2es10.3,a,2es10.3)') "  step ", n, " (top): residual vs the current wind", wx, wy, &
               "; vs the previous step's wind", wx_old, wy_old
            call check(error, max(wx, wy) <= 1.0e-12_wp*max(sc, 1.0e-3_wp), "ocean stress = blend with the CURRENT wind (top)")
            if (allocated(error)) exit
            call check(error, max(wx_old, wy_old) > 1.0e-9_wp, "and it is not the blend with the previous step's wind")
            if (allocated(error)) exit
         else
            call check(error, max(wx, wy) <= 1.0e-12_wp*max(sc, 1.0e-3_wp), "ocean stress = blend (first step)")
            if (allocated(error)) exit
         end if
         if (allocated(wind_prev_x)) deallocate (wind_prev_x, wind_prev_y)
         wind_prev_x = engine%state%ice%tau_a_x
         wind_prev_y = engine%state%ice%tau_a_y
         call engine_step_ice(engine, cfg, DT, t, ierr=ierr)
         call check(error, ierr == OCEAN_STATUS_OK, "engine_step_ice")
         if (allocated(error)) exit
         call pull_stress(engine)
         ! after the ice step: the hand-off blend, with the wind the ice just felt and the NEW ice state
         call blend_residual(engine, engine%state%ice%tau_a_x, engine%state%ice%tau_a_y, wx, wy, sc)
         call check(error, max(wx, wy) <= 1.0e-12_wp*max(sc, 1.0e-3_wp), "ocean stress = blend with the current wind (after)")
         if (allocated(error)) exit
         call engine_step_finalize(engine, DT, t, ierr=ierr)
         call check(error, ierr == OCEAN_STATUS_OK, "engine_step_finalize")
         if (allocated(error)) exit
         t = t + DT
      end do
      if (.not. allocated(error)) then
         amax = maxval(engine%state%ice%dync%w%aiE(is:ie, js:je))
         call check(error, amax > 0.3_wp, "the faces are partly ice covered")
      end if
      call finish_engine(engine)
      if (allocated(error)) then
         call remove_atm_files()
         return
      end if

      ! ice-free: the ocean stress IS the unweighted wind stress
      call remove_atm_files()
      call write_atm_files("warm")
      call make_engine(ef, cff, case_nml("basin", 0.0_wp, 1, .true.), ok)
      call check(error, ok, "ice-free engine setup failed")
      if (allocated(error)) then
         call remove_atm_files()
         return
      end if
      tf = 0.0_wp
      do n = 1, 4
         call engine_step(ef, DT, tf, ierr=ierr)
         call check(error, ierr == OCEAN_STATUS_OK, "ice-free engine_step")
         if (allocated(error)) exit
         call pull_stress(ef)
         dmax = 0.0_wp
         sc = 0.0_wp
         do j = js, je
            do i = is + 1, ie
               dmax = max(dmax, abs(ef%state%surface_stress%tau_x(i, j) - ef%state%ice%tau_a_x(i, j)))
               sc = max(sc, abs(ef%state%ice%tau_a_x(i, j)))
            end do
         end do
         write (*, '(a,i0,a,es10.3,a,es10.3)') "  ice-free step ", n, ": |tau - wind stress| worst", dmax, " of", sc
         call check(error, sc > 1.0e-4_wp, "there is a wind stress")
         if (allocated(error)) exit
         call check(error, dmax <= 1.0e-14_wp*sc, "an ice-free face carries the wind stress unchanged")
         if (allocated(error)) exit
         call engine_step_ice(ef, cff, DT, tf, ierr=ierr)
         call engine_step_finalize(ef, DT, tf, ierr=ierr)
         tf = tf + DT
      end do
      call finish_engine(ef)
      call remove_atm_files()
   end subroutine test_stress_blend

   subroutine face_mean_host(engine, i0, j0, i1, j1, q, r)
      !! Area-weighted, wet-masked two-point mean of the T field `q` between cells `(i0, j0)` and `(i1, j1)`.
      type(ocean_engine_t), intent(in) :: engine
      integer, intent(in) :: i0, j0, i1, j1
      real(wp), intent(in) :: q(:, :)
      real(wp), intent(out) :: r
      real(wp) :: w1, w2

      w1 = engine%state%ice_cgrid%hm(i0, j0)*engine%state%ice_cgrid%tarea(i0, j0)
      w2 = engine%state%ice_cgrid%hm(i1, j1)*engine%state%ice_cgrid%tarea(i1, j1)
      r = 0.0_wp
      if (w1 + w2 /= 0.0_wp) r = (w1*q(i0, j0) + w2*q(i1, j1))/(w1 + w2)
   end subroutine face_mean_host

   subroutine test_face_fraction_current(error)
      !! The wind share of the ocean stress at the top of a step is `1 - (face mean of the CURRENT aice)`,
      !! computed here on the host from the pulled aggregate; it is the complement of the open-water weight
      !! (`aice0`) every other flux uses, and it differs from the concentration the dynamics left (one step older).
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      real(wp), allocatable :: ai_dyn(:, :)
      real(wp) :: t, ax, ay, worst, comp, lag
      integer :: n, i, j, ierr, ng, is, ie, js, je
      logical :: ok

      call write_atm_files("cold")
      call make_engine(engine, cfg, case_nml("basin", 0.6_wp, 1, .true.), ok)
      call check(error, ok, "engine setup failed")
      if (allocated(error)) then
         call remove_atm_files()
         return
      end if
      ng = engine%grid%nghost
      is = ng + 1
      ie = ng + engine%grid%nx_phys
      js = ng + 1
      je = ng + engine%grid%ny_phys
      allocate (ai_dyn(size(engine%state%ice%dync%w%aiE, 1), size(engine%state%ice%dync%w%aiE, 2)))
      ai_dyn = 0.0_wp
      t = 0.0_wp
      lag = 0.0_wp
      do n = 1, 6
         call engine_step(engine, DT, t, ierr=ierr)
         call check(error, ierr == OCEAN_STATUS_OK, "engine_step")
         if (allocated(error)) exit
         call pull_stress(engine)
         call pull_cice(engine)
         worst = 0.0_wp
         comp = 0.0_wp
         do j = js, je
            do i = is + 1, ie
               call face_mean_host(engine, i - 1, j, i, j, engine%state%ice%aice, ax)
               worst = max(worst, abs(engine%state%ice%dync%w%aiE(i - 1, j) - ax))
               call face_mean_host(engine, i - 1, j, i, j, engine%state%ice%aice0, ay)
               comp = max(comp, abs((1.0_wp - engine%state%ice%dync%w%aiE(i - 1, j)) - ay))
               if (n > 1) lag = max(lag, abs(engine%state%ice%dync%w%aiE(i - 1, j) - ai_dyn(i - 1, j)))
            end do
         end do
         write (*, '(a,i0,a,es10.3,a,es10.3)') "  step ", n, ": aiE vs the host face mean of the current aice", worst, &
            "; (1 - aiE) vs the face mean of aice0 ", comp
         call check(error, worst <= 1.0e-14_wp, "the stress weight is the face mean of the CURRENT concentration")
         if (allocated(error)) exit
         call check(error, comp <= 1.0e-14_wp, "wind share + open-water weight add up to one")
         if (allocated(error)) exit
         call engine_step_ice(engine, cfg, DT, t, ierr=ierr)
         call check(error, ierr == OCEAN_STATUS_OK, "engine_step_ice")
         if (allocated(error)) exit
         call pull_stress(engine)
         ai_dyn = engine%state%ice%dync%w%aiE      ! the concentration the dynamics left (the older one)
         call engine_step_finalize(engine, DT, t, ierr=ierr)
         t = t + DT
      end do
      if (.not. allocated(error)) then
         write (*, '(a,es10.3)') "  largest face-concentration change vs the dynamics' one: ", lag
         call check(error, lag > 1.0e-4_wp, "the concentration changes over a step (the lag the rebuild removes is real)")
      end if
      call finish_engine(engine)
      call remove_atm_files()
   end subroutine test_face_fraction_current

   subroutine test_consumed_stress(error)
      !! The stress the ocean dynamics CONSUMES (captured after `engine_step`, before the ice step overwrites
      !! it with the hand-off blend) against quantities built independently on the host: the closed-form bulk
      !! wind stress (not the exposed array), the host face mean of the current concentration, and the
      !! previous ice step's ice-ocean stress.  PROVES: the consumed stress is `(1 - a_cur) tau_wind - strocn_prev`
      !! per face from those inputs.  DOES NOT PROVE: the ice-side momentum balance (strair + strocn against the
      !! ice acceleration and internal stress needs the dynamics' internals); that rests on `test_cice_dyn_c_*`.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      real(wp) :: t, ax, tl, tr, e, worst, sc, wy_worst
      integer :: n, i, j, ierr, ng, is, ie, js, je
      logical :: ok

      call write_atm_files("cold")
      call make_engine(engine, cfg, case_nml("basin", 0.6_wp, 1, .true.), ok)
      call check(error, ok, "engine setup failed")
      if (allocated(error)) then
         call remove_atm_files()
         return
      end if
      ng = engine%grid%nghost
      is = ng + 1
      ie = ng + engine%grid%nx_phys
      js = ng + 1
      je = ng + engine%grid%ny_phys
      t = 0.0_wp
      do n = 1, 6
         call engine_step(engine, DT, t, ierr=ierr)
         call check(error, ierr == OCEAN_STATUS_OK, "engine_step")
         if (allocated(error)) exit
         call pull_stress(engine)
         call pull_cice(engine)
         call pull_atm(engine)
         !$acc update self(engine%state%multilayer%tracers(engine%state%multilayer%idx_temperature)%c)
         worst = 0.0_wp
         wy_worst = 0.0_wp
         sc = 0.0_wp
         do j = js, je
            do i = is + 1, ie
               call expected_wind_stress(engine, cfg, i - 1, j, tl, e)
               call expected_wind_stress(engine, cfg, i, j, tr, e)
               call face_mean_host(engine, i - 1, j, i, j, engine%state%ice%aice, ax)
               e = (1.0_wp - ax)*0.5_wp*(tl + tr) - engine%state%ice%dync%w%strocnxE(i - 1, j)
               worst = max(worst, abs(engine%state%surface_stress%tau_x(i, j) - e))
               sc = max(sc, abs(e))
            end do
         end do
         do j = js + 1, je
            do i = is, ie
               call expected_wind_stress(engine, cfg, i, j - 1, e, tl)
               call expected_wind_stress(engine, cfg, i, j, e, tr)
               call face_mean_host(engine, i, j - 1, i, j, engine%state%ice%aice, ax)
               e = (1.0_wp - ax)*0.5_wp*(tl + tr) - engine%state%ice%dync%w%strocnyN(i, j - 1)
               wy_worst = max(wy_worst, abs(engine%state%surface_stress%tau_y(i, j) - e))
               sc = max(sc, abs(e))
            end do
         end do
         write (*, '(a,i0,a,2es10.3,a,es10.3)') "  step ", n, ": consumed stress vs the independent value (x, y)", worst, &
            wy_worst, " of", sc
         call check(error, sc > 1.0e-3_wp, "there is stress")
         if (allocated(error)) exit
         call check(error, max(worst, wy_worst) <= 1.0e-4_wp*sc, "the consumed stress is the independent blend")
         if (allocated(error)) exit
         call engine_step_ice(engine, cfg, DT, t, ierr=ierr)
         call engine_step_finalize(engine, DT, t, ierr=ierr)
         t = t + DT
      end do
      call finish_engine(engine)
      call remove_atm_files()
   end subroutine test_consumed_stress

   ! ------------------------------------------------------------------
   ! atmosphere inputs and flux partition
   ! ------------------------------------------------------------------

   subroutine test_ice_inputs(error)
      !! (2) the ice's atmosphere inputs are the prescribed forcing of this step, in the units and rotation the
      !! bulk forcing uses; the downward shortwave for the ice radiation is the forcing's `rsds`.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      real(wp) :: t, ca, sa, e_u, e_v, e_pt, e_rho, worst, sw_err
      integer :: n, i, j, ierr, ng, is, ie, js, je
      logical :: ok

      call write_atm_files("cold")
      call make_engine(engine, cfg, case_nml("basin", 0.6_wp, 1, .true.), ok)
      call check(error, ok, "engine setup failed")
      if (allocated(error)) then
         call remove_atm_files()
         return
      end if
      ng = engine%grid%nghost
      is = ng + 1
      ie = ng + engine%grid%nx_phys
      js = ng + 1
      je = ng + engine%grid%ny_phys
      t = 0.0_wp
      do n = 1, 4
         call engine_step(engine, DT, t, ierr=ierr)
         call engine_step_ice(engine, cfg, DT, t, ierr=ierr)
         call check(error, ierr == OCEAN_STATUS_OK, "engine_step_ice")
         if (allocated(error)) exit
         call pull_atm(engine)
         call pull_ice_flux(engine)
         worst = 0.0_wp
         sw_err = 0.0_wp
         associate (a => engine%state%atm, w => engine%state%icec%work, m => engine%state%metrics)
            do j = js, je
               do i = is, ie
                  ca = cos(m%angle_dx(i, j))
                  sa = sin(m%angle_dx(i, j))
                  e_u = ca*a%uas(i, j) + sa*a%vas(i, j)
                  e_v = -sa*a%uas(i, j) + ca*a%vas(i, j)
                  e_pt = a%tas(i, j)*(1.0e5_wp/a%psl(i, j))**0.286_wp
                  e_rho = a%psl(i, j)/(287.058_wp*(1.0_wp + 0.608_wp*a%huss(i, j))*a%tas(i, j))
                  worst = max(worst, abs(w%uatm(i, j) - e_u)/max(abs(e_u), 1.0_wp), &
                              abs(w%vatm(i, j) - e_v)/max(abs(e_v), 1.0_wp), &
                              abs(w%potT(i, j) - e_pt)/e_pt, abs(w%rhoa(i, j) - e_rho)/e_rho, &
                              abs(w%Qa(i, j) - a%huss(i, j))/1.0e-3_wp, &
                              abs(w%flw(i, j) - a%rlds(i, j))/a%rlds(i, j), &
                              abs(w%fsnow(i, j) - a%prsn(i, j))/a%prsn(i, j), &
                              abs(w%frain(i, j) - a%prra(i, j))/a%prra(i, j))
                  sw_err = max(sw_err, abs(w%sw_down(i, j) - a%rsds(i, j))/a%rsds(i, j))
               end do
            end do
            call check(error, maxval(abs(a%uas(is:ie, js:je))) > 1.0_wp, "there is a wind")
            call check(error, maxval(abs(a%vas(is:ie, js:je) - a%vas(is, js))) > 1.0e-3_wp, "the wind varies in space")
         end associate
         write (*, '(a,i0,a,es10.3,a,es10.3)') "  step ", n, ": ice inputs vs the forcing, worst relative ", worst, &
            "; sw_down ", sw_err
         call check(error, worst < 1.0e-12_wp, "the ice's atmosphere inputs are the forcing's")
         if (allocated(error)) exit
         call check(error, sw_err == 0.0_wp, "the ice's downward shortwave is the forcing's rsds")
         if (allocated(error)) exit
         call engine_step_finalize(engine, DT, t, ierr=ierr)
         t = t + DT
      end do
      call finish_engine(engine)
      call remove_atm_files()
   end subroutine test_ice_inputs

   subroutine test_flux_partition(error)
      !! (3) per flux, the open-water share carries weight `aice0` and the ice share reaches the ocean through the
      !! hand-over only; the two add up to one cell.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      real(wp) :: t, ow, rsds, sst, ca, sa, ug, vg, tx, ty, sens, lat, evp, cd, ch, ce, ustar, lw_up
      real(wp) :: e_sw, e_lw, e_sens, e_lat, e_evap, e_lp, e_fp, e_ice_salt, esw_ice
      real(wp) :: w_sw, w_lw, w_sens, w_lat, w_evap, w_lp, w_fp, w_part, w_salt, w_h, w_isw
      real(wp) :: a0(2), albedo
      integer :: n, i, j, ierr, ng, is, ie, js, je, kt
      logical :: ok

      call write_atm_files("cold")
      call make_engine(engine, cfg, case_nml("basin", 0.6_wp, 1, .true.), ok)
      call check(error, ok, "engine setup failed")
      if (allocated(error)) then
         call remove_atm_files()
         return
      end if
      albedo = cfg%ocean%atm_forcing%ocean_albedo
      ng = engine%grid%nghost
      is = ng + 1
      ie = ng + engine%grid%nx_phys
      js = ng + 1
      je = ng + engine%grid%ny_phys
      t = 0.0_wp
      do n = 1, 5
         call engine_step(engine, DT, t, ierr=ierr)
         call check(error, ierr == OCEAN_STATUS_OK, "engine_step")
         if (allocated(error)) exit
         call pull_atm(engine)
         call pull_cice(engine)
         call pull_sf(engine)
         !$acc update self(engine%state%multilayer%tracers(engine%state%multilayer%idx_temperature)%c)
         w_sw = 0.0_wp
         w_lw = 0.0_wp
         w_sens = 0.0_wp
         w_lat = 0.0_wp
         w_evap = 0.0_wp
         w_lp = 0.0_wp
         w_fp = 0.0_wp
         w_part = 0.0_wp
         associate (a => engine%state%atm, m => engine%state%metrics, sf => engine%state%surface_flux, &
                    ms => engine%state%multilayer, c => engine%state%ice%cice)
            kt = ms%nz_ml
            do j = js, je
               do i = is, ie
                  ow = engine%state%ice%aice0(i, j)
                  w_part = max(w_part, abs(ow + c%aice(i, j) - 1.0_wp))
                  sst = ms%tracers(ms%idx_temperature)%c(i, j, kt)
                  ca = cos(m%angle_dx(i, j))
                  sa = sin(m%angle_dx(i, j))
                  ug = ca*a%uas(i, j) + sa*a%vas(i, j)
                  vg = -sa*a%uas(i, j) + ca*a%vas(i, j)
                  call bulk_flux_column(ug, vg, 0.0_wp, 0.0_wp, a%tas(i, j), a%huss(i, j), a%psl(i, j), sst, &
                                        cfg%ocean%bulk_flux%z_wind, cfg%ocean%bulk_flux%z_ta, &
                                        cfg%ocean%bulk_flux%n_iter, tx, ty, sens, lat, evp, cd, ch, ce, ustar)
                  lw_up = bulk_flux_longwave_up(sst + BULK_T0_KELVIN, cfg%ocean%bulk_flux%emissivity)
                  rsds = a%rsds(i, j)
                  e_sw = ow*(1.0_wp - albedo)*rsds
                  e_lw = ow*(a%rlds(i, j) - lw_up)
                  e_sens = ow*sens
                  e_lat = ow*lat
                  e_evap = -ow*evp
                  e_lp = ow*a%prra(i, j)
                  e_fp = ow*a%prsn(i, j)
                  w_sw = max(w_sw, abs(sf%q_sw(i, j) - e_sw)/rsds)
                  w_lw = max(w_lw, abs(sf%q_lw(i, j) - e_lw)/abs(a%rlds(i, j)))
                  w_sens = max(w_sens, abs(sf%q_sens(i, j) - e_sens)/(abs(sens) + 1.0_wp))
                  w_lat = max(w_lat, abs(sf%q_lat(i, j) - e_lat)/(abs(lat) + 1.0_wp))
                  w_evap = max(w_evap, abs(sf%evap(i, j) - e_evap)/(abs(evp) + 1.0e-9_wp))
                  w_lp = max(w_lp, abs(sf%lprec(i, j) - e_lp)/a%prra(i, j))
                  w_fp = max(w_fp, abs(sf%fprec(i, j) - e_fp)/a%prsn(i, j))
               end do
            end do
            a0 = [minval(engine%state%ice%aice0(is:ie, js:je)), maxval(engine%state%ice%aice0(is:ie, js:je))]
         end associate
         write (*, '(a,i0,a,2f7.3,a,es9.2)') "  step ", n, " open fraction range", a0, "; aice0 + aice - 1 worst", w_part
         write (*, '(a,7es10.2)') "    open-water share vs closed form (sw, lw, sens, lat, evap, lprec, fprec):", &
            w_sw, w_lw, w_sens, w_lat, w_evap, w_lp, w_fp
         call check(error, a0(1) > 0.0_wp .and. a0(2) < 1.0_wp, "the cells are partly ice covered")
         if (allocated(error)) exit
         call check(error, w_part < 1.0e-14_wp, "open + ice fractions add up to one")
         if (allocated(error)) exit
         call check(error, w_sw < 1.0e-12_wp, "open-water shortwave = aice0 (1 - albedo) SW")
         if (allocated(error)) exit
         call check(error, w_lw < 1.0e-3_wp, "open-water longwave = aice0 (LW down - emission; the surface temperature of the closed form is one ocean step newer)")
         if (allocated(error)) exit
         call check(error, w_sens < 1.0e-3_wp, "open-water sensible heat = aice0 x the bulk value")
         if (allocated(error)) exit
         call check(error, w_lat < 1.0e-3_wp, "open-water latent heat = aice0 x the bulk value")
         if (allocated(error)) exit
         call check(error, w_evap < 1.0e-3_wp, "evaporation = aice0 x the bulk value")
         if (allocated(error)) exit
         call check(error, w_lp < 1.0e-14_wp .and. w_fp < 1.0e-14_wp, "rain and snow on the open fraction = aice0 x the forcing")
         if (allocated(error)) exit
         call engine_step_ice(engine, cfg, DT, t, ierr=ierr)
         call check(error, ierr == OCEAN_STATUS_OK, "engine_step_ice")
         if (allocated(error)) exit
         call pull_sf(engine)
         call pull_cice(engine)
         call pull_ice_flux(engine)
         w_salt = 0.0_wp
         w_h = 0.0_wp
         w_isw = 0.0_wp
         esw_ice = 0.0_wp
         associate (sf => engine%state%surface_flux, f => engine%state%icec%flux, w => engine%state%icec%work, &
                    c => engine%state%ice%cice, o => engine%state%icec%oin, a => engine%state%atm)
            do j = js, je
               do i = is, ie
                  ! ratio 1: the window is this step; the ice share is the ice's own hand-over
                  w_h = max(w_h, abs(sf%heat_ice(i, j) - f%fhocn(i, j))/(abs(f%fhocn(i, j)) + 1.0_wp))
                  w_isw = max(w_isw, abs(sf%sw_ice(i, j) - f%fswthru(i, j))/(abs(f%fswthru(i, j)) + 1.0_wp))
                  ! rain through the ice is NOT diluted: the salt hand-over is the ice's own water and salt
                  e_ice_salt = 1000.0_wp*f%fsalt(i, j) - o%sss_ocn(i, j)*(f%fresh(i, j) - w%frain(i, j)*c%aice_init(i, j))
                  w_salt = max(w_salt, abs(sf%salt_ice(i, j) - e_ice_salt)/(abs(e_ice_salt) + 1.0e-6_wp))
                  ! shortwave: open water + transmitted through the ice stays within the incident shortwave
                  esw_ice = max(esw_ice, (sf%q_sw(i, j) + sf%sw_ice(i, j) - a%rsds(i, j))/a%rsds(i, j))
               end do
            end do
         end associate
         write (*, '(a,4es10.2)') "    ice share: heat, sw, salt vs the hand-over; (open + transmitted - incident)/incident max:", &
            w_h, w_isw, w_salt, esw_ice
         call check(error, w_h < 1.0e-13_wp, "ice heat = fhocn")
         if (allocated(error)) exit
         call check(error, w_isw < 1.0e-13_wp, "ice shortwave = fswthru (transmitted through the ice)")
         if (allocated(error)) exit
         call check(error, w_salt < 1.0e-9_wp, "ice salt = 1000 fsalt - S (fresh - rain through the ice)")
         if (allocated(error)) exit
         call check(error, esw_ice <= 1.0e-12_wp, "open-water + transmitted shortwave never exceeds the incident")
         if (allocated(error)) exit
         call engine_step_finalize(engine, DT, t, ierr=ierr)
         t = t + DT
      end do
      call finish_engine(engine)
      call remove_atm_files()
   end subroutine test_flux_partition

   ! ------------------------------------------------------------------
   ! closures
   ! ------------------------------------------------------------------

   subroutine run_closed(engine, cfg, t, nsteps, ok, cl, ratio)
      !! Step with the ice ledger around every ice step.  The external sources of the ledger are the
      !! per-cell inputs the ice read this step (`work%flw, fsnow, frain`).  The ocean side: the heat and salt
      !! the assembler applied (`Q_heat`, `Q_salt` of every step, `ratio = 1`) against the change of the ocean
      !! content, and the ice share of those against the ice's hand-over.
      type(ocean_engine_t), intent(inout), target :: engine
      type(config_t), intent(inout) :: cfg
      real(wp), intent(inout) :: t
      integer, intent(in) :: nsteps
      logical, intent(out) :: ok
      type(closure_t), intent(out) :: cl
      integer, intent(in) :: ratio
      type(ice_ledger_mirror_t) :: led
      real(wp) :: h0, s0, h1, s1, hext, wext, efrz, spent, stock, ap_h, ap_s, ice_asm, ice_led, lfresh, scl
      real(wp) :: dh, ds, pend_h, pend_s
      integer :: n, i, j, ng, is, ie, js, je
      real(wp), allocatable :: area(:, :)

      ng = engine%grid%nghost
      is = ng + 1
      ie = ng + engine%grid%nx_phys
      js = ng + 1
      je = ng + engine%grid%ny_phys
      ok = .true.
      cl = closure_t()
      spent = 0.0_wp
      stock = 0.0_wp
      ap_h = 0.0_wp
      ap_s = 0.0_wp
      pend_h = 0.0_wp
      pend_s = 0.0_wp
      ice_asm = 0.0_wp
      ice_led = 0.0_wp
      lfresh = engine%state%ice%cice%p%c%Lfresh
      allocate (area(size(engine%state%metrics%areaT, 1), size(engine%state%metrics%areaT, 2)))
      area = engine%state%metrics%areaT
      call ocean_totals(engine, h0, s0)
      do n = 1, nsteps
         stock = 0.0_wp
         call led%open(engine%state%ice%cice, engine%state%ice_cgrid%tarea, DT, 1)
         call one_step(engine, cfg, t, ok)
         if (.not. ok) return
         call pull_cice(engine)
         call pull_ice_flux(engine)
         call pull_sf(engine)
         associate (c => engine%state%ice%cice, f => engine%state%icec%flux, w => engine%state%icec%work, &
                    sf => engine%state%surface_flux)
            hext = 0.0_wp
            wext = 0.0_wp
            do j = js, je
               do i = is, ie
                  hext = hext + DT*area(i, j)*c%aice_init(i, j)*(w%flw(i, j) - lfresh*w%fsnow(i, j))
                  wext = wext + DT*area(i, j)*(w%fsnow(i, j) + w%frain(i, j))*c%aice_init(i, j)
                  efrz = DT*area(i, j)*max(engine%state%ice%frzmlt(i, j), 0.0_wp) &
                         - area(i, j)*engine%state%ice%frazil_heat(i, j)
                  hext = hext - efrz
                  spent = spent + efrz
                  stock = stock + area(i, j)*engine%state%ice%frazil_heat(i, j)
                  if (ratio == 1) then
                     ice_asm = ice_asm + area(i, j)*DT*(sf%heat_ice(i, j) + sf%sw_ice(i, j))
                  end if
               end do
            end do
            if (ratio == 1) then
               ! the ocean applies the flux assembled by the PREVIOUS step (assembled at the end of a step,
               ! applied by the next one)
               ap_h = ap_h + pend_h
               ap_s = ap_s + pend_s
               pend_h = 0.0_wp
               pend_s = 0.0_wp
               do j = js, je
                  do i = is, ie
                     pend_h = pend_h + area(i, j)*DT*sf%Q_heat(i, j)
                     pend_s = pend_s + area(i, j)*DT*sf%Q_salt(i, j)
                  end do
               end do
            end if
            call led%close(c, f, engine%state%ice_cgrid%tarea, wext, hext)
            ice_led = ice_led + led%i_fhocn + led%i_fswthru
            cl%water = max(cl%water, abs(led%res(ICE_LEDGER_WATER))/max(led%scale(ICE_LEDGER_WATER), tiny(1.0_wp)))
            cl%salt = max(cl%salt, abs(led%res(ICE_LEDGER_SALT))/max(led%scale(ICE_LEDGER_SALT), tiny(1.0_wp)))
            cl%heat = max(cl%heat, abs(led%res(ICE_LEDGER_HEAT))/max(led%scale(ICE_LEDGER_HEAT), tiny(1.0_wp)))
            if (.not. ieee_is_finite(led%res(1) + led%res(2) + led%res(3))) cl%finite = .false.
            if (any(iand(c%ice_err, ICE_ERR_FATAL_MASK) /= 0_int32)) then
               cl%finite = .false.
               write (*, '(a,i0,a,i0)') "  ice error bits at step ", n, ": ", int(maxval(iand(c%ice_err, ICE_ERR_FATAL_MASK)))
            end if
            if (n == 1) cl%vice_first = sum(c%vice(is:ie, js:je))
            cl%vice_last = sum(c%vice(is:ie, js:je))
            cl%aice_max = max(cl%aice_max, maxval(c%aice(is:ie, js:je)))
         end associate
      end do
      call ocean_totals(engine, h1, s1)
      if (ratio == 1) then
         dh = h1 - h0 - (stock + spent)
         ds = s1 - s0
         scl = max(abs(h1 - h0), abs(ap_h), abs(stock + spent), 1.0_wp)
         cl%ocean_heat = abs(dh - ap_h)/scl
         scl = max(abs(ds), abs(ap_s), 1.0_wp)
         cl%ocean_salt = abs(ds - ap_s)/scl
  write (*, '(a,4es12.4)') "  ocean heat: dH, bank (stock+spent), applied Q_heat, ice share:", h1 - h0, stock + spent, ap_h, ice_asm
         ! the ice share the assembler holds against the ice's own hand-over (this step's)
         scl = max(abs(ice_asm), abs(ice_led), 1.0_wp)
         cl%ice_share = abs(ice_asm - ice_led)/scl
      end if
   end subroutine run_closed

   subroutine report(label, cl)
      character(len=*), intent(in) :: label
      type(closure_t), intent(in) :: cl
      write (*, '(a,a)') "  ", label
      write (*, '(a,3es11.3)') "    ice ledger water, salt, heat (worst relative):", cl%water, cl%salt, cl%heat
      write (*, '(a,2es11.3,a,es11.3)') "    ocean heat, salt vs the applied fluxes:", cl%ocean_heat, cl%ocean_salt, &
         "; ice share held vs handed over:", cl%ice_share
      write (*, '(a,3es11.3)') "    vice first/last, aice max:", cl%vice_first, cl%vice_last, cl%aice_max
   end subroutine report

   subroutine closure_case(error, geom, conc, ratio, nsteps, label)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: geom, label
      real(wp), intent(in) :: conc
      integer, intent(in) :: ratio, nsteps
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      type(closure_t) :: cl
      real(wp) :: t
      logical :: ok

      call write_atm_files("cold")
      call make_engine(engine, cfg, case_nml(geom, conc, ratio, .true.), ok)
      call check(error, ok, "engine setup failed")
      if (allocated(error)) then
         call remove_atm_files()
         return
      end if
      t = 0.0_wp
      call run_closed(engine, cfg, t, nsteps, ok, cl, ratio)
      call check(error, ok, "the run completed")
      if (.not. allocated(error)) then
         call report(label, cl)
         call check(error, cl%finite, "finite, no fatal ice error bit (the counted ridging-pass bit is allowed)")
         if (.not. allocated(error)) call check(error, cl%water < 3.0e-11_wp, "ice water ledger closes with the forcing")
         if (.not. allocated(error)) call check(error, cl%salt < 3.0e-11_wp, "ice salt ledger closes with the forcing")
         if (.not. allocated(error)) call check(error, cl%heat < 3.0e-9_wp, "ice heat ledger closes with the forcing")
         if (.not. allocated(error) .and. ratio == 1) then
            call check(error, cl%ocean_heat < 1.0e-11_wp, "the ocean heat change is the applied heat flux")
         if (.not. allocated(error)) call check(error, cl%ocean_salt < 1.0e-11_wp, "the ocean salt change is the applied salt flux")
   if (.not. allocated(error)) call check(error, cl%ice_share < 1.0e-9_wp, "the ice share held by the ocean is the ice's hand-over")
         end if
         if (.not. allocated(error)) call check(error, cl%aice_max > 0.3_wp .and. cl%aice_max <= 1.0_wp + 1.0e-12_wp, &
                                                "the pack exists and area <= 1")
      end if
      call finish_engine(engine)
      call remove_atm_files()
   end subroutine closure_case

   subroutine test_closure_column(error)
      type(error_type), allocatable, intent(out) :: error
      call closure_case(error, "column", 0.9_wp, 1, 12, "column, concentration 0.9, ratio 1")
   end subroutine test_closure_column

   subroutine test_closure_basin(error)
      type(error_type), allocatable, intent(out) :: error
      call closure_case(error, "basin", 0.6_wp, 1, 18, "basin, concentration 0.6, ratio 1")
   end subroutine test_closure_basin

   subroutine test_closure_basin_ratio2(error)
      type(error_type), allocatable, intent(out) :: error
      call closure_case(error, "basin", 0.6_wp, 2, 18, "basin, concentration 0.6, ratio 2")
   end subroutine test_closure_basin_ratio2

   ! ------------------------------------------------------------------
   ! restart
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
            if (snap%n >= MAXF) error stop "test_cice_forcing: raise MAXF"
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

   subroutine compare(error, sa, sb, when, n_face)
      type(error_type), allocatable, intent(out) :: error
      type(snapshot_t), intent(in) :: sa, sb
      character(len=*), intent(in) :: when
      integer, intent(out) :: n_face
      integer :: e, i, j, k, nbad, nb
      character(len=:), allocatable :: msg

      n_face = 0
      nbad = 0
      msg = ""
      call check(error, sa%n == sb%n, "registry field count differs "//when)
      if (allocated(error)) return
      do e = 1, sa%n
         if (sa%f(e)%tag /= sb%f(e)%tag) then
            call check(error, .false., "registry field order differs: "//trim(sa%f(e)%tag))
            return
         end if
         if (index(sa%f(e)%tag, "ice_dync_strocn") == 1) n_face = n_face + 1
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
      real(wp) :: t_a, t_b, umax, wmax
      integer :: step_b, ierr, nw, na, n_face
      logical :: ok
      character(len=:), allocatable :: nml_pen

      nw = 5          ! odd: the checkpoint falls inside a thermodynamics window (ratio 2)
      ! Shortwave penetration ON through `sw_source = 'q_sw'`: the through-ice shortwave `sw_ice`
      ! is then read EVERY step via `sw_net` but refreshed only at the thermo hand-off, so the
      ! mid-window checkpoint is bitwise only if the restart carries it (`sf_sw_ice`).
      nml_pen = replace(case_nml("basin", 0.6_wp, 2, .true.), &
                        "[ocean.thermo]"//new_line('a')// &
                        "enable_thermodynamics = true"//new_line('a'), &
                        "[ocean.thermo]"//new_line('a')// &
                        "enable_thermodynamics = true"//new_line('a')// &
                        "sw_pen_frac = 0.4"//new_line('a')// &
                        'sw_source = "q_sw"'//new_line('a'))
      na = 4
      call write_atm_files("cold")
      call make_engine(ea, cfg_a, nml_pen, ok)
      call check(error, ok, "engine A setup failed")
      if (allocated(error)) then
         call remove_atm_files()
         return
      end if
      t_a = 0.0_wp
      call advance(ea, cfg_a, t_a, nw, ok)
      call check(error, ok, "engine A failed before the checkpoint")
      if (allocated(error)) then
         call finish_engine(ea)
         call remove_atm_files()
         return
      end if
      call ocean_state_restart_write(ea%state, ea%grid, ea%decomp, FN_RT, t_a, nw, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "checkpoint write failed")
      if (allocated(error)) then
         call finish_engine(ea)
         call remove_atm_files()
         return
      end if
      call take_snapshot(ea, sa0)
      call advance(ea, cfg_a, t_a, na, ok)
      call check(error, ok, "engine A failed after the checkpoint")
      if (.not. allocated(error)) then
         call take_snapshot(ea, sa)
         !$acc update self(ea%state%ice%u_ice)
         !$acc update self(ea%state%ice%tau_a_x)
         umax = maxval(abs(ea%state%ice%u_ice))
         wmax = maxval(abs(ea%state%ice%tau_a_x))
         call check(error, umax > 1.0e-7_wp, "the ice moves in the restart case")
         if (.not. allocated(error)) call check(error, wmax > 1.0e-4_wp, "and there is a wind stress")
      end if
      call finish_engine(ea)
      if (allocated(error)) then
         call remove_atm_files()
         return
      end if

      call make_engine(eb, cfg_b, nml_pen, ok, restart_file=FN_RT, t0=t_b, step0=step_b)
      call check(error, ok, "engine B (warm restart) setup failed")
      if (allocated(error)) then
         call remove_atm_files()
         return
      end if
      call check(error, step_b == nw .and. eb%warm_restart, "warm restart taken at the written step")
      if (.not. allocated(error)) then
         call take_snapshot(eb, sb0)
         call compare(error, sa0, sb0, "at the resume point", n_face)
      end if
      if (.not. allocated(error)) call check(error, n_face >= 2, "the ice-ocean stress tags are in the comparison")
      if (.not. allocated(error)) then
         call advance(eb, cfg_b, t_b, na, ok)
         call check(error, ok, "engine B failed after the restart")
      end if
      if (.not. allocated(error)) then
         call take_snapshot(eb, sb)
         call compare(error, sa, sb, "after the resumed steps", n_face)
      end if
      call finish_engine(eb)
      call delete_file(FN_RT)
      call remove_atm_files()
   end subroutine test_restart

   subroutine tag_flags(engine, req_x, req_y, found)
      !! Whether the two ice-ocean stress tags are registered, and REQUIRED (not optional).
      type(ocean_engine_t), intent(inout), target :: engine
      logical, intent(out) :: req_x, req_y, found
      type(restart_registry_t) :: reg
      integer :: e

      req_x = .false.
      req_y = .false.
      found = .false.
      call ocean_state_build_restart_registry(engine%state, engine%grid, reg)
      do e = 1, reg%n
         if (trim(reg%entries(e)%tag) == "ice_dync_strocnxE") then
            req_x = .not. reg%entries(e)%optional
            found = .true.
         end if
         if (trim(reg%entries(e)%tag) == "ice_dync_strocnyN") req_y = .not. reg%entries(e)%optional
      end do
   end subroutine tag_flags

   subroutine test_restart_tags_required(error)
      !! Under the prescribed atmosphere a checkpoint WITHOUT the ice-ocean stress would resume with the bare
      !! wind under the ice: the tags are required there (the registry's missing-required-field path is fatal),
      !! optional under the uniform atmosphere.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: ea, eb
      type(config_t) :: cfg_a, cfg_b
      logical :: ok, rx, ry, found

      call write_atm_files("cold")
      call make_engine(ea, cfg_a, case_nml("column", 0.6_wp, 1, .true.), ok)
      call check(error, ok, "prescribed engine setup failed")
      if (.not. allocated(error)) then
         call tag_flags(ea, rx, ry, found)
         call check(error, found, "the stress tags are registered under the prescribed atmosphere")
         if (.not. allocated(error)) call check(error, rx .and. ry, "and REQUIRED there")
         call finish_engine(ea)
      end if
      call remove_atm_files()
      if (allocated(error)) return
      call make_engine(eb, cfg_b, case_nml("column", 0.6_wp, 1, .false.), ok)
      call check(error, ok, "uniform engine setup failed")
      if (allocated(error)) return
      call tag_flags(eb, rx, ry, found)
      call check(error, found, "the stress tags are registered under the uniform atmosphere")
      if (.not. allocated(error)) call check(error,.not. rx .and. .not. ry, "and optional there")
      call finish_engine(eb)
   end subroutine test_restart_tags_required

   function replace(s, a, b) result(r)
      character(len=*), intent(in) :: s, a, b
      character(len=:), allocatable :: r
      integer :: k

      k = index(s, a)
      if (k == 0) then
         r = s
      else
         r = s(:k - 1)//b//s(k + len(a):)
      end if
   end function replace

   subroutine delete_file(fn_)
      character(len=*), intent(in) :: fn_
      integer :: u, ios

      open (newunit=u, file=fn_, status="old", iostat=ios)
      if (ios == 0) close (u, status="delete")
   end subroutine delete_file

end module test_cice_forcing
