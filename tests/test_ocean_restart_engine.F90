!! Bit-exact warm restart through the PRODUCTION engine path.
!!
!! `test_ocean_restart` builds its two states by hand and never runs
!! `engine_setup`, so it cannot see what the configure chain does to a
!! state AFTER the restart read.  Three defects lived exactly there, and
!! all three made a 1/4-degree Southern Ocean resume diverge from the run
!! that wrote the checkpoint (2026-10-01):
!!   1. the `pred_corr` predictor's carried viscous tendency
!!      (`hvisc%du_visc`/`dv_visc`) was not checkpointed, and the scratch
!!      attach zero-filled it even once it was;
!!   2. `configure_ocean_land_mask` re-seeded every land column to uniform
!!      `H_VANISHED`, while the running model (the ALE remap regrids land
!!      columns like any other) had carried a different layout;
!!   3. the init-time periodic wrap, the host halo exchange and the
!!      device warm-up exchange re-derived the checkpointed ghosts, and at
!!      a step boundary the duplicated seam faces of `u` are not a pure
!!      function of the owned interior.
!!
!! Gate: engine A steps N_WRITE steps, writes a checkpoint, steps
!! N_AFTER more.  Engine B resumes from the checkpoint and steps N_AFTER.
!! EVERY restart-registry field -- FULL local arrays, ghosts included,
!! and the registered scalars -- must be bitwise identical, and so must
!! the carried surface forcing the next step would read (`tau_x/tau_y`,
!! `stress_mag`, `Q_heat/Q_salt`).  Each step is the driver's sequence:
!! `engine_step` -> `engine_step_ice` -> `engine_step_finalize`.
!!
!! Variants:
!!   * island_periodic_zfixed -- land (an island) under `z_fixed` (so the
!!     remap reshapes the land columns), a periodic axis (seam ghosts), and
!!     lateral viscosity under the default `pred_corr` (the carried
!!     tendency);
!!   * sea_ice_periodic{,_ssp_rk2} -- sea ice with thermo + ITD, EVP
!!     dynamics (CFL clip, `project_ci`), the ice->ocean stress blend, the
!!     category transport, snowfall and penetrating shortwave, in a cooled
!!     re-entrant channel.  `dt_therm_ratio = 2` with N_WRITE odd puts the
!!     checkpoint in the MIDDLE of a thermo window: the ice fluxes the
!!     ocean integrates on the resumed step come from the window that
!!     closed before the checkpoint, so the restart must carry them.  The
!!     periodic edge puts ice across the wrap seam, whose ghosts the
!!     checkpoint holds and setup must not re-derive (the cold-start ice
!!     exchange is skipped on a warm restart);
!!   * sea_ice_components -- the same under `&ocean_forcing_nml
!!     enable_components`, where `Q_heat`/`Q_salt` are ASSEMBLED at the
!!     end of each thermo window and carried, unchanged, through the steps
!!     until the next one;
!!   * gm_meke_varmix_therm2 -- Gent-McWilliams + MEKE + VarMix over a
!!     stratified seamount under `dt_therm_ratio = 2`, N_WRITE odd.  The GM
!!     operator runs AFTER the dynamics of every outer step and leaves
!!     `gm_src`, which MEKE reads at the top of the next thermo step; on
!!     the step between two thermo refreshes it also reads the slopes and
!!     the VarMix(+MEKE) KhTh the last refresh left.  All of it is carried
!!     across the checkpoint, so all of it must be restart-registered.
!!   * seamount_periodic_zstar_visc_rem_kshear_vertex -- compat row c069:
!!     the vertex kappa-shear closure's corner carrier `kshear%kd_corner`
!!     is only (re)allocated by `configure_ocean_kappa_shear`, which runs
!!     AFTER the restart read in `engine_setup` -- so unlike every other
!!     carried field above, it cannot be put in the registry the READ
!!     sees unless something allocates it EARLIER.  A warm resume's first
!!     stage read it at the cold `0` instead of the checkpointed value.
module test_ocean_restart_engine
   use, intrinsic :: iso_fortran_env, only: int64
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_config, only: config_t, read_config_from_string, validate_config
   use rdb_ocean_engine, only: ocean_engine_t, engine_setup, engine_enter_data, &
                               engine_step, engine_step_ice, engine_step_finalize, &
                               engine_exit_data, engine_teardown
   use rdb_ocean_state, only: ocean_state_restart_write, ocean_state_build_restart_registry, &
                              ocean_state_restart_write_drop_field
   use rdb_ocean_restart, only: restart_registry_t
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   use rdb_comm_env, only: comm_env_init, comm_env_setup_roles
   use rdb_io_netcdf, only: nc_create_file, nc_close, nc_def_dim, nc_def_var_2d, &
                            nc_def_var_3d, nc_enddef, rdb_def_var_1d, rdb_put_var_1d, nc_put_att
   use netcdf, only: nf90_put_var
   implicit none
   private

   public :: collect_ocean_restart_engine_tests

   real(wp), parameter :: DT = 600.0_wp
   character(len=*), parameter :: FN = "test_ocean_restart_engine_rt.nc"
   character(len=*), parameter :: ATM_FILE = "rsw_atm.nc"
   character(len=*), parameter :: ATM_LIST = "rsw_atm.list"
   character(len=*), parameter :: ATM_SW0_FILE = "rsw_atm_sw0.nc"
   character(len=*), parameter :: ATM_SW0_LIST = "rsw_atm_sw0.list"
   integer, parameter :: MAXF = 160
   character(len=*), parameter :: NL = new_line("a")

   logical :: comm_inited = .false.

   type :: field_t
      character(len=64) :: tag = ""
      real(wp), allocatable :: a(:, :, :)
   end type field_t

   type :: snapshot_t
      !! Host copies of every registry field (full local extent; a scalar
      !! is a 1x1x1 array) plus the carried surface forcing.
      integer :: n = 0
      type(field_t) :: f(MAXF)
   end type snapshot_t

contains

   subroutine collect_ocean_restart_engine_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("restart_engine_bit_exact_island_periodic_zfixed", &
                               test_engine_bit_exact), &
                  new_unittest("restart_engine_bit_exact_sea_ice_periodic", &
                               test_engine_bit_exact_ice), &
                  new_unittest("restart_engine_bit_exact_sea_ice_periodic_ssp_rk2", &
                               test_engine_bit_exact_ice_ssp), &
                  new_unittest("restart_engine_bit_exact_sea_ice_components", &
                               test_engine_bit_exact_ice_components), &
                  new_unittest("restart_engine_bit_exact_sea_ice_components_bulk", &
                               test_engine_bit_exact_ice_components_bulk), &
                  new_unittest("engine_ice_enabled_icefree_shortwave_matches_ice_disabled", &
                               test_engine_icefree_sw), &
                  new_unittest("restart_engine_bit_exact_gm_meke_varmix_therm2", &
                               test_engine_bit_exact_gm), &
                  new_unittest("restart_engine_bit_exact_visc_rem_chain", &
                               test_engine_bit_exact_visc_rem), &
                  new_unittest("restart_engine_bit_exact_visc_rem_kshear_vertex", &
                               test_engine_bit_exact_visc_rem_kshear_vertex), &
                  new_unittest("restart_engine_bit_exact_bkgnd_profile", &
                               test_engine_bit_exact_bkgnd_profile), &
                  new_unittest("restart_engine_bkgnd_full_reseed_on_restart", &
                               test_engine_bkgnd_full_reseed_on_restart), &
                  new_unittest("restart_engine_cold_seeds_kv_without_checkpoint", &
                               test_engine_cold_seeds_kv_without_checkpoint), &
                  new_unittest("restart_engine_cavity_old_checkpoint_loads", &
                               test_engine_cavity_old_checkpoint_loads), &
                  new_unittest("restart_engine_bit_exact_cavity_melt_drag", &
                               test_engine_bit_exact_cavity_melt_drag), &
                  new_unittest("restart_engine_bit_exact_cavity_melt_only", &
                               test_engine_bit_exact_cavity_melt_only), &
                  new_unittest("restart_engine_bit_exact_cavity_therm2_window_end", &
                               test_engine_bit_exact_cavity_therm2_end), &
                  new_unittest("restart_engine_bit_exact_cavity_therm2_mid_window", &
                               test_engine_bit_exact_cavity_therm2_mid), &
                  new_unittest("restart_engine_bit_exact_cavity_melt_only_therm2", &
                               test_engine_bit_exact_cavity_melt_only_therm2), &
                  new_unittest("restart_engine_cavity_ustar_not_vacuous", &
                               test_engine_cavity_ustar_not_vacuous), &
                  new_unittest("restart_engine_bit_exact_cavity_mass", &
                               test_engine_bit_exact_cavity_mass), &
                  new_unittest("restart_engine_bit_exact_cavity_mass_uniform_open", &
                               test_engine_bit_exact_cavity_mass_uo), &
                  new_unittest("restart_engine_bit_exact_cavity_ustar", &
                               test_engine_bit_exact_cavity_ustar), &
                  new_unittest("restart_engine_bit_exact_cavity_ssp_rk2", &
                               test_engine_bit_exact_cavity_ssp), &
                  new_unittest("restart_engine_bit_exact_cavity_tdrag_sigma", &
                               test_engine_bit_exact_cavity_tdrag_sigma) &
                  ]
   end subroutine collect_ocean_restart_engine_tests

   subroutine ensure_comm()
      !! See `test_ocean_budget_periodic_sponge_serial::ensure_comm`.
      if (.not. comm_inited) then
         call comm_env_init()
         call comm_env_setup_roles(.false.)
         comm_inited = .true.
      end if
   end subroutine ensure_comm

   function case_nml(label) result(nml)
      !! The namelist of one variant (see the module docstring).
      character(len=*), intent(in) :: label
      character(len=:), allocatable :: nml
      character(len=*), parameter :: NL = new_line("a")
      if (index(label, "cavity_") == 1) then
         nml = cavity_nml(label)
         return
      end if
      select case (label)
      case ("island_periodic_zfixed")
         ! Periodic channel with an island, `z_fixed` (tanh) layers, wind,
         ! Smagorinsky + a Laplacian floor, KPP (default), `pred_corr`
         ! (default).
         nml = "&sim_nml sim_type = 'ocean' /"//NL// &
               "&time_nml t_end = 86400.0, dt_fixed = 600.0 /"//NL// &
               "&nonhydrostatic_nml nz_layers = 6 /"//NL// &
               "&tracer_nml initial_temperature = 12.0, initial_salinity = 35.0, "// &
               "T_init_surface = 20.0, T_init_bottom = 4.0 /"//NL// &
               "&ocean_bt_nml auto_n_inner = .true. /"//NL// &
               "&ocean_hvisc_nml nu_h = 200.0, lateral_closure = 'smagorinsky', "// &
               "smag_ah = .true. /"//NL// &
               "&ocean_diag_nml enabled = .false. /"//NL// &
               "&grid_nml nx = 24, ny = 16, nghost = 3, dx = 10000.0, dy = 10000.0 /"//NL// &
               "&physics_nml coriolis_f = 1.0e-4, wind_stress_x = 0.08 /"//NL// &
               "&vcoord_nml vcoord_type = 'z_fixed', z_fixed_profile = 'tanh', "// &
               "z_fixed_dz_top = 20.0, z_fixed_tanh_center = 0.5, "// &
               "z_fixed_tanh_width = 0.25 /"//NL// &
               "&ocean_topo_nml topo_config = 'island', max_depth = 1000.0, "// &
               "slope_scale = 0.25 /"//NL// &
               "&ocean_bc_nml west = 'periodic', east = 'periodic', south = 'wall', "// &
               "north = 'wall' /"//NL// &
               "&output_nml output_to_file = .false. /"//NL
      case ("island_periodic_zfixed_visc_rem")
         ! PR-1: the SAME island/periodic/z_fixed/pred_corr case as
         ! `island_periodic_zfixed`, but with the visc_rem chain live
         ! (`implicit_drag` + `visc_rem_chain`) so
         ! `bt_work%visc_rem_u/v` holds NON-trivial values (not just the
         ! `source=1.0` init) at the checkpoint — the restart test that
         ! closes compat row `restart_visc_rem`.  Linear bed drag with
         ! `hbbl = 0` (no HBBL band — `implicit_drag` + `hbbl > 0` stays
         ! refused at configure).
         nml = "&sim_nml sim_type = 'ocean' /"//NL// &
               "&time_nml t_end = 86400.0, dt_fixed = 600.0 /"//NL// &
               "&nonhydrostatic_nml nz_layers = 6 /"//NL// &
               "&tracer_nml initial_temperature = 12.0, initial_salinity = 35.0, "// &
               "T_init_surface = 20.0, T_init_bottom = 4.0 /"//NL// &
               "&ocean_bt_nml auto_n_inner = .true., visc_rem_chain = .true. /"//NL// &
               "&ocean_hvisc_nml nu_h = 200.0, lateral_closure = 'smagorinsky', "// &
               "smag_ah = .true. /"//NL// &
               "&ocean_bdrag_nml form = 'linear', r = 2.0e-4, hbbl = 0.0, bg_vel = 0.1 /"//NL// &
               "&ocean_vdiff_nml implicit_drag = .true. /"//NL// &
               "&ocean_diag_nml enabled = .false. /"//NL// &
               "&grid_nml nx = 24, ny = 16, nghost = 3, dx = 10000.0, dy = 10000.0 /"//NL// &
               "&physics_nml coriolis_f = 1.0e-4, wind_stress_x = 0.08 /"//NL// &
               "&vcoord_nml vcoord_type = 'z_fixed', z_fixed_profile = 'tanh', "// &
               "z_fixed_dz_top = 20.0, z_fixed_tanh_center = 0.5, "// &
               "z_fixed_tanh_width = 0.25 /"//NL// &
               "&ocean_topo_nml topo_config = 'island', max_depth = 1000.0, "// &
               "slope_scale = 0.25 /"//NL// &
               "&ocean_bc_nml west = 'periodic', east = 'periodic', south = 'wall', "// &
               "north = 'wall' /"//NL// &
               "&output_nml output_to_file = .false. /"//NL
      case ("seamount_periodic_zstar_visc_rem_kshear_vertex")
         ! c069 restart fix regression.  A wind-driven channel over a
         ! seamount (real depth variation -> real vertical shear, needed
         ! so the vertex kappa-shear corner solve produces a non-trivial
         ! `kd_corner` -- a flat-bottom / weakly-sheared column leaves it
         ! at exactly 0 whether warm-started or not, which would make this
         ! case pass for the wrong reason) on `zstar` (MOM6 z*), with the
         ! kappa-shear VERTEX closure live (`&ocean_kappa_shear_nml
         ! enable = .true., at_vertex = .true.`) and the visc_rem chain
         ! (`&ocean_bt_nml visc_rem_chain = .true.`).
         ! `wind_stress_x =
         ! 2.0` is deliberately strong -- it is what pushes the corner
         ! solve's converged `kd_corner` well above `kappa_trunc`
         ! (~1e-9) by step 12, so the bug (below) is not masked by both
         ! sides happening to floor to the same 0.
         !
         ! `visc_rem_precompute` -- called at the START of every stage,
         ! BEFORE that stage's own `vmix_apply_in_stage` refreshes the
         ! kappa-shear corner carrier `kshear%kd_corner` -- reads
         ! `kd_corner` as the vertex form's corner Kv source on EVERY
         ! stage, including the FIRST stage of a warm-restarted run.
         ! `kd_corner` is only ever (re)allocated by
         ! `configure_ocean_kappa_shear`'s `init_vertex` call, which runs
         ! well AFTER the restart read in `engine_setup` -- so without
         ! the early, knob-gated `init_vertex` call added ahead of the
         ! restart read (`rdb_ocean_engine.F90`), a warm resume's first
         ! stage read `kd_corner` at the cold `0` `init_vertex` seeds it
         ! with instead of the checkpointed value, perturbing
         ! `bt_visc_rem_u/v` at the `kappa_trunc` scale and, from there,
         ! every carried field (compat row c069). This is the MINIMAL
         ! reproduction found by bisecting c069's 16-axis cell: only
         ! `vcoord = zstar` + `kappa_shear(at_vertex)` + the visc_rem
         ! chain are load-bearing; geometry/tripolar, GM/MEKE, KPP,
         ! Bryan-Lewis and pseudo-salt are NOT (the real c069 cell fails
         ! identically with or without them, and passes identically with
         ! or without them once this fix is in).
         nml = "&sim_nml sim_type = 'ocean' /"//NL// &
               "&time_nml t_end = 86400.0, dt_fixed = 600.0 /"//NL// &
               "&nonhydrostatic_nml nz_layers = 10 /"//NL// &
               "&tracer_nml initial_temperature = 9.0, initial_salinity = 34.8, "// &
               "T_init_surface = 16.0, T_init_bottom = 3.0, "// &
               "S_init_surface = 34.6, S_init_bottom = 35.0 /"//NL// &
               "&ocean_bt_nml auto_n_inner = .true., split_scheme = 'pred_corr', "// &
               "visc_rem_chain = .true. /"//NL// &
               "&ocean_hvisc_nml nu_h = 10.0 /"//NL// &
               "&ocean_bdrag_nml form = 'quadratic', cd = 0.003, hbbl = 0.0, "// &
               "bg_vel = 0.05 /"//NL// &
               "&ocean_vdiff_nml implicit_drag = .true. /"//NL// &
               "&ocean_kappa_shear_nml enable = .true., at_vertex = .true. /"//NL// &
               "&ocean_vmix_nml use_closure = .true., use_kpp = .false. /"//NL// &
               "&ocean_thermo_nml enable_thermodynamics = .true., q_heat = -60.0 /"//NL// &
               "&ocean_diag_nml enabled = .false. /"//NL// &
               "&grid_nml nx = 24, ny = 16, nghost = 4, dx = 20000.0, dy = 20000.0 /"//NL// &
               "&physics_nml coriolis_f = 1.0e-4, wind_stress_x = 2.0 /"//NL// &
               "&vcoord_nml vcoord_type = 'zstar' /"//NL// &
               "&ocean_topo_nml topo_config = 'seamount', max_depth = 2000.0, "// &
               "edge_depth = 1600.0, slope_scale = 80000.0, "// &
               "coriolis_beta = 2.0e-11 /"//NL// &
               "&ocean_bc_nml west = 'periodic', east = 'periodic', south = 'wall', "// &
               "north = 'wall' /"//NL// &
               "&output_nml output_to_file = .false. /"//NL
      case ("island_periodic_zfixed_bkgnd")
         ! PR-2 review: the SAME island/periodic/z_fixed/pred_corr/
         ! visc_rem-chain case as `island_periodic_zfixed_visc_rem`
         ! (closed-BC: south/north = 'wall'), PLUS `&ocean_vmix_nml
         ! bkgnd_profile` (Bryan-Lewis) -- `kd_bg` is set ONLY by
         ! `vmix_seed_backgrounds`, so a warm restart that wrongly skips
         ! the FULL seed (the bug this case was added to catch: the
         ! first `reseed_arrays` fix skipped kt/ks/kd_bg too, not just
         ! kv) loses the background diffusivity entirely post-resume.
         nml = "&sim_nml sim_type = 'ocean' /"//NL// &
               "&time_nml t_end = 86400.0, dt_fixed = 600.0 /"//NL// &
               "&nonhydrostatic_nml nz_layers = 6 /"//NL// &
               "&tracer_nml initial_temperature = 12.0, initial_salinity = 35.0, "// &
               "T_init_surface = 20.0, T_init_bottom = 4.0 /"//NL// &
               "&ocean_bt_nml auto_n_inner = .true., visc_rem_chain = .true. /"//NL// &
               "&ocean_hvisc_nml nu_h = 200.0, lateral_closure = 'smagorinsky', "// &
               "smag_ah = .true. /"//NL// &
               "&ocean_bdrag_nml form = 'linear', r = 2.0e-4, hbbl = 0.0, bg_vel = 0.1 /"//NL// &
               "&ocean_vdiff_nml implicit_drag = .true. /"//NL// &
               "&ocean_vmix_nml bkgnd_profile = .true., bkgnd_kd_sfc = 2.0e-5, "// &
               "bkgnd_kd_deep = 3.0e-4, bkgnd_z0 = 1500.0, "// &
               "pp81_kappa_bg = 4.0e-5 /"//NL// &
               "&ocean_diag_nml enabled = .false. /"//NL// &
               "&grid_nml nx = 24, ny = 16, nghost = 3, dx = 10000.0, dy = 10000.0 /"//NL// &
               "&physics_nml coriolis_f = 1.0e-4, wind_stress_x = 0.08 /"//NL// &
               "&vcoord_nml vcoord_type = 'z_fixed', z_fixed_profile = 'tanh', "// &
               "z_fixed_dz_top = 20.0, z_fixed_tanh_center = 0.5, "// &
               "z_fixed_tanh_width = 0.25 /"//NL// &
               "&ocean_topo_nml topo_config = 'island', max_depth = 1000.0, "// &
               "slope_scale = 0.25 /"//NL// &
               "&ocean_bc_nml west = 'periodic', east = 'periodic', south = 'wall', "// &
               "north = 'wall' /"//NL// &
               "&output_nml output_to_file = .false. /"//NL
      case ("sea_ice_periodic", "sea_ice_periodic_ssp_rk2", "sea_ice_components")
         ! The decomposition harness's sea-ice channel (cold, cooled,
         ! oblique wind over a seamount, 0.7 conc of 1 m ice), with the
         ! category transport, snowfall and a shortwave on top.  The
         ! column sits just above the liquidus so the ice survives the
         ! basal flux (see that harness's comment), the EVP stresses are
         ! live, and `cfl_trunc` is tight enough that the final clip fires.
         nml = "&sim_nml sim_type = 'ocean' /"//NL// &
               "&time_nml t_end = 86400.0, dt_fixed = 600.0 /"//NL// &
               "&nonhydrostatic_nml nz_layers = 4 /"//NL// &
               "&tracer_nml initial_temperature = -1.83, initial_salinity = 34.0, "// &
               "T_init_surface = -1.82, T_init_bottom = -1.835 /"//NL// &
               "&ocean_bt_nml auto_n_inner = .true., split_scheme = '"// &
               trim(merge("ssp_rk2  ", "pred_corr", label == "sea_ice_periodic_ssp_rk2"))//"' /"//NL// &
               "&ocean_hvisc_nml nu_h = 200.0, lateral_closure = 'smagorinsky', "// &
               "smag_ah = .true. /"//NL// &
               "&ocean_diag_nml enabled = .false. /"//NL// &
               "&grid_nml nx = 24, ny = 16, nghost = 3, dx = 20000.0, dy = 20000.0 /"//NL// &
               "&physics_nml coriolis_f = 1.4e-4, wind_stress_x = 0.1, wind_stress_y = 0.04 /"//NL// &
               "&vcoord_nml vcoord_type = 'sigma' /"//NL// &
               "&ocean_topo_nml topo_config = 'seamount', max_depth = 1000.0, "// &
               "edge_depth = 800.0, slope_scale = 60000.0 /"//NL// &
               "&ocean_thermo_nml enable_thermodynamics = .true., q_heat = -100.0 /"//NL// &
               "&ocean_vmix_nml dt_therm_ratio = 2 /"//NL// &
               "&ocean_ice_nml enable = .true., ncat = 5, nk_ice = 2, air_temp = -20.0, "// &
               "restore_lambda = 20.0, sw_down = 40.0, snowfall = 2.0e-5, "// &
               "dynamics = .true., evp_sub_steps = 30, cfl_trunc = 0.01, "// &
               "project_ci = .true., transport = .true., adv_substeps = 2 /"//NL// &
               "&ocean_ice_ic_nml conc_config = 'uniform', h_ice = 1.0, conc = 0.7 /"//NL// &
               "&ocean_bc_nml west = 'periodic', east = 'periodic', south = 'wall', "// &
               "north = 'wall' /"//NL// &
               "&output_nml output_to_file = .false. /"//NL
         if (label == "sea_ice_components") then
            nml = nml//"&ocean_forcing_nml enable_components = .true. /"//NL
         end if
      case ("sea_ice_components_bulk")
         ! Ice + components + the atmosphere-driven bulk forcing: the path on
         ! which the ice couplers used to overwrite the bulk consumer's
         ! open-water `q_sw` (the decomposition harness's `bulk_flux_ice`,
         ! on one rank).  The transmitted shortwave, the open-water
         ! shortwave and the brine/heat components are all live, and the
         ! checkpoint lands mid thermo window.
         nml = "&sim_nml sim_type = 'ocean' /"//NL// &
               "&time_nml t_end = 86400.0, dt_fixed = 600.0, start_date = '2000-06-21' /"//NL// &
               "&nonhydrostatic_nml nz_layers = 4 /"//NL// &
               "&ocean_bt_nml auto_n_inner = .true., split_scheme = 'pred_corr' /"//NL// &
               "&ocean_hvisc_nml nu_h = 200.0, lateral_closure = 'smagorinsky', "// &
               "smag_ah = .true. /"//NL// &
               "&ocean_diag_nml enabled = .false. /"//NL// &
               "&ocean_eos_nml eos = 'wright' /"//NL// &
               "&ocean_thermo_nml enable_thermodynamics = .true. /"//NL// &
               "&ocean_forcing_nml enable_components = .true. /"//NL// &
               "&ocean_bulk_flux_nml enable = .true., n_iter = 5, z_wind = 10.0, "// &
               "z_ta = 10.0 /"//NL// &
               "&tracer_nml initial_temperature = -1.83, initial_salinity = 34.0, "// &
               "T_init_surface = -1.82, T_init_bottom = -1.835 /"//NL// &
               atm_group(ATM_LIST)// &
               "&ocean_vmix_nml dt_therm_ratio = 2 /"//NL// &
               "&ocean_ice_nml enable = .true., ncat = 5, nk_ice = 2, "// &
               "restore_lambda = 20.0, dynamics = .true., evp_sub_steps = 30, "// &
               "cfl_trunc = 0.01, project_ci = .true., transport = .true., "// &
               "adv_substeps = 2 /"//NL// &
               "&ocean_ice_ic_nml conc_config = 'uniform', h_ice = 1.0, conc = 0.7 /"//NL// &
               "&grid_nml nx = 24, ny = 16, nghost = 3, dx = 1.0, dy = 1.0 /"//NL// &
               "&ocean_grid_nml grid_config = 'spherical', lon_west = 0.0, "// &
               "lat_south = 20.0, rad_earth = 6.371e6, coriolis_scheme = 'planetary' /"//NL// &
               "&vcoord_nml vcoord_type = 'sigma' /"//NL// &
               "&ocean_topo_nml topo_config = 'spoon', max_depth = 3000.0, "// &
               "edge_depth = 300.0, slope_scale = 300000.0 /"//NL// &
               "&ocean_bc_nml west = 'periodic', east = 'periodic', south = 'wall', "// &
               "north = 'wall' /"//NL// &
               "&output_nml output_to_file = .false. /"//NL
      case ("gm_meke_varmix_therm2")
         ! Stratified seamount on z* with GM + MEKE + VarMix (the compat
         ! matrix's `gm_meke` overlay), thermo every second step.
         nml = "&sim_nml sim_type = 'ocean' /"//NL// &
               "&time_nml t_end = 86400.0, dt_fixed = 600.0 /"//NL// &
               "&nonhydrostatic_nml nz_layers = 6 /"//NL// &
               "&tracer_nml initial_temperature = 12.0, initial_salinity = 35.0, "// &
               "T_init_surface = 20.0, T_init_bottom = 4.0 /"//NL// &
               "&ocean_bt_nml auto_n_inner = .true. /"//NL// &
               "&ocean_hvisc_nml nu_h = 200.0 /"//NL// &
               "&ocean_diag_nml enabled = .false. /"//NL// &
               "&grid_nml nx = 24, ny = 16, nghost = 3, dx = 10000.0, dy = 10000.0 /"//NL// &
               "&physics_nml coriolis_f = 1.0e-4, wind_stress_x = 0.08 /"//NL// &
               "&vcoord_nml vcoord_type = 'zstar' /"//NL// &
               "&ocean_topo_nml topo_config = 'seamount', max_depth = 1000.0, "// &
               "edge_depth = 800.0, slope_scale = 40000.0 /"//NL// &
               "&ocean_thermo_nml enable_thermodynamics = .true. /"//NL// &
               "&ocean_vmix_nml dt_therm_ratio = 2 /"//NL// &
               "&ocean_slopes_nml enable = .true. /"//NL// &
               "&ocean_gm_nml enable = .true., khth = 500.0 /"//NL// &
               "&ocean_meke_nml enable = .true., gmcoeff = 0.15, khcoeff = 1.0, "// &
               "damping = 1.0e-6 /"//NL// &
               "&ocean_wavespeed_nml enable = .true. /"//NL// &
               "&ocean_varmix_nml enable = .true., use_visbeck = .true., "// &
               "khth_slope_cff = 0.1, khtr_slope_cff = 0.1, visbeck_l_scale = 3.0e4, "// &
               "khth_max = 2000.0, khtr_max = 2000.0 /"//NL// &
               "&ocean_bc_nml west = 'periodic', east = 'periodic', south = 'wall', "// &
               "north = 'wall' /"//NL// &
               "&output_nml output_to_file = .false. /"//NL
      case default
         error stop "test_ocean_restart_engine: unknown case"
      end select
   end function case_nml

   function cavity_nml(label) result(nml)
      !! Ice-shelf cavity with basal melt (and optionally top drag) over a
      !! seamount, closed basin, `enable_components`: the engine-driven
      !! cavity namelists of the cavity engine tests.  The label selects the
      !! variant by substring: `zfixed` = z_fixed (no KPP, 8 fixed-z layers,
      !! top drag on) else sigma with the default KPP and 4 layers;
      !! `tdrag` = top drag on a sigma cavity; `ustar` = sigma, 16 thin
      !! layers and a large tidal friction velocity (`u_tide = 0.2`), so the
      !! under-ice KPP boundary layer is active; `therm2` = `dt_therm_ratio
      !! = 2`; `ssp` = `ssp_rk2` (else `pred_corr`); `mass` =
      !! `freshwater = 'mass'`; `_uo` = `volume_compensation =
      !! 'uniform_open_ocean'`.
      character(len=*), intent(in) :: label
      character(len=:), allocatable :: nml
      character(len=*), parameter :: NL = new_line("a")
      character(len=:), allocatable :: vc, nzl, cav, vm, td, bt, melt_extra, vk
      logical :: zfixed, top_drag
      zfixed = index(label, "zfixed") > 0
      top_drag = zfixed .or. index(label, "tdrag") > 0
      vk = ""
      melt_extra = ""
      if (zfixed) then
         nzl = "8"
         vc = "&vcoord_nml vcoord_type = 'z_fixed', zfixed_closed_faces = .true. /"//NL
         cav = ", h_min_cavity = 250.0 /"//NL
         vk = "use_kpp = .false., pp81_nu0 = 0.0, pp81_nu_bg = 1.0e-3, pp81_kappa_bg = 5.0e-5"
      else
         nzl = "4"
         if (index(label, "ustar") > 0) then
            nzl = "16"
            melt_extra = ", u_tide = 0.2"
         end if
         vc = "&vcoord_nml vcoord_type = 'sigma' /"//NL
         cav = ", h_min_cavity = 100.0 /"//NL
      end if
      if (index(label, "therm2") > 0) then
         if (len(vk) > 0) vk = vk//", "
         vk = vk//"dt_therm_ratio = 2"
      end if
      vm = ""
      if (len(vk) > 0) vm = "&ocean_vmix_nml "//vk//" /"//NL
      if (index(label, "mass") > 0) melt_extra = melt_extra//", freshwater = 'mass'"
      if (index(label, "_uo") > 0) melt_extra = melt_extra//", volume_compensation = 'uniform_open_ocean'"
      bt = "&ocean_bt_nml auto_n_inner = .true., split_scheme = 'pred_corr' /"//NL
      if (index(label, "ssp") > 0) bt = "&ocean_bt_nml auto_n_inner = .true., split_scheme = 'ssp_rk2' /"//NL
      td = ""
      if (top_drag) td = "&ocean_tdrag_nml enable = .true., form = 'quadratic', cd = 2.5e-3 /"//NL
      nml = "&sim_nml sim_type = 'ocean' /"//NL// &
            "&time_nml t_end = 86400.0, dt_fixed = 600.0 /"//NL// &
            "&nonhydrostatic_nml nz_layers = "//nzl//" /"//NL// &
            "&tracer_nml initial_salinity = 34.5, T_init_surface = 3.0, T_init_bottom = 1.0 /"//NL// &
            bt// &
            "&ocean_hvisc_nml nu_h = 200.0, lateral_closure = 'smagorinsky', smag_ah = .true. /"//NL// &
            "&ocean_diag_nml enabled = .false. /"//NL// &
            "&grid_nml nx = 26, ny = 18, nghost = 3, dx = 20000.0, dy = 20000.0 /"//NL// &
            "&physics_nml coriolis_f = -1.409e-4 /"//NL// &
            vc//vm// &
            "&ocean_topo_nml topo_config = 'seamount', max_depth = 1000.0, edge_depth = 400.0, "// &
            "slope_scale = 80000.0, wind_config = 'constant', taux_magnitude = 0.05 /"//NL// &
            "&ocean_cavity_dyn_nml enable = .true., draft_config = 'linear', "// &
            "draft_depth = 500.0, draft_slope = -4.0e-4, draft_x0 = 0.0, "// &
            "draft_x1 = 380000.0, draft_y1 = 260000.0"//cav// &
            "&ocean_forcing_nml enable_components = .true. /"//NL// &
            "&ocean_cavity_melt_nml enable = .true., exchange_law = 'const_gamma', "// &
            "gamma_t = 2.2e-2, gamma_s = -1.0, cdrag_top = 2.5e-3, "// &
            "ice_conduction = 'insulating', s_ice = 0.0, far_field_depth = 10.0"// &
            melt_extra//" /"//NL// &
            td// &
            "&ocean_pgf_nml form = 'fv_mom6', p_top_in_bc = .true. /"//NL// &
            "&ocean_eos_nml eos = 'linear', tfreeze_set = 'isomip' /"//NL// &
            "&ocean_bc_nml west = 'wall', east = 'wall', south = 'wall', north = 'wall' /"//NL// &
            "&output_nml output_to_file = .false. /"//NL
   end function cavity_nml

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
      !! The driver's step sequence (the ice block is a no-op unless
      !! `&ocean_ice_nml enable`).
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

   subroutine take_snapshot(engine, snap)
      !! Every restart-registry field, device -> host first (COMPONENT
      !! arrays only, never the aggregate derived type -- see CLAUDE.md),
      !! full local extent, plus the carried surface forcing.
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
            select case (en%rank)
            case (0)
               call push(snap, en%tag, reshape([en%p0], [1, 1, 1]))
            case (2)
               call push(snap, en%tag, reshape(en%p2, [size(en%p2, 1), size(en%p2, 2), 1]))
            case default
               call push(snap, en%tag, en%p3)
            end select
         end associate
      end do
      !$acc update self(engine%state%surface_stress%tau_x, engine%state%surface_stress%tau_y)
      !$acc update self(engine%state%surface_stress%stress_mag)
      !$acc update self(engine%state%surface_flux%Q_heat, engine%state%surface_flux%Q_salt)
      call push2(snap, "surface_stress_tau_x", engine%state%surface_stress%tau_x)
      call push2(snap, "surface_stress_tau_y", engine%state%surface_stress%tau_y)
      call push2(snap, "surface_stress_mag", engine%state%surface_stress%stress_mag)
      call push2(snap, "surface_flux_Q_heat", engine%state%surface_flux%Q_heat)
      call push2(snap, "surface_flux_Q_salt", engine%state%surface_flux%Q_salt)
   end subroutine take_snapshot

   subroutine push(snap, tag, a)
      type(snapshot_t), intent(inout) :: snap
      character(len=*), intent(in) :: tag
      real(wp), intent(in) :: a(:, :, :)
      if (snap%n >= MAXF) error stop "test_ocean_restart_engine: raise MAXF"
      snap%n = snap%n + 1
      snap%f(snap%n)%tag = tag
      snap%f(snap%n)%a = a
   end subroutine push

   subroutine push2(snap, tag, a)
      type(snapshot_t), intent(inout) :: snap
      character(len=*), intent(in) :: tag
      real(wp), intent(in) :: a(:, :)
      call push(snap, tag, reshape(a, [size(a, 1), size(a, 2), 1]))
   end subroutine push2

   subroutine compare(error, sa, sb, n_ice, when)
      !! Bitwise, every element of every field (ghosts included).  Reports
      !! every mismatching field (first cell each) before failing.
      type(error_type), allocatable, intent(out) :: error
      type(snapshot_t), intent(in) :: sa, sb
      character(len=*), intent(in) :: when
      integer, intent(out) :: n_ice
         !! How many `ice_*` registry fields were compared (a guard that
         !! the ice variants really reach the ice state).
      integer :: e, i, j, k, nbad_fields, nb
      character(len=:), allocatable :: msg
      character(len=200) :: line

      n_ice = 0
      nbad_fields = 0
      msg = ""
      call check(error, sa%n == sb%n, "registry field count differs "//when)
      if (allocated(error)) return
      do e = 1, sa%n
         if (sa%f(e)%tag /= sb%f(e)%tag) then
            call check(error, .false., "registry field order differs: "//trim(sa%f(e)%tag)// &
                       " vs "//trim(sb%f(e)%tag))
            return
         end if
         if (index(sa%f(e)%tag, "ice_") == 1) n_ice = n_ice + 1
         if (any(shape(sa%f(e)%a) /= shape(sb%f(e)%a))) then
            nbad_fields = nbad_fields + 1
            msg = msg//" "//trim(sa%f(e)%tag)//"(shape)"
            cycle
         end if
         nb = 0
         line = ""
         do k = 1, size(sa%f(e)%a, 3)
            do j = 1, size(sa%f(e)%a, 2)
               do i = 1, size(sa%f(e)%a, 1)
                  if (transfer(sa%f(e)%a(i, j, k), 0_int64) /= &
                      transfer(sb%f(e)%a(i, j, k), 0_int64)) then
                     if (nb == 0) write (line, '(a,3(1x,i0),a,es24.16,a,es24.16)') &
                        when//": "//trim(sa%f(e)%tag)//" first at", i, j, k, ": ", &
                        sa%f(e)%a(i, j, k), " vs ", sb%f(e)%a(i, j, k)
                     nb = nb + 1
                  end if
               end do
            end do
         end do
         if (nb > 0) then
            write (*, '(a,a,i0,a)') trim(line), "  (", nb, " cells)"
            nbad_fields = nbad_fields + 1
            msg = msg//" "//trim(sa%f(e)%tag)
         end if
      end do
      call check(error, nbad_fields == 0, "fields differ "//when//":"//msg)
   end subroutine compare

   pure function itoa(n) result(str)
      integer, intent(in) :: n
      character(len=:), allocatable :: str
      character(len=16) :: buf
      write (buf, '(i0)') n
      str = trim(buf)
   end function itoa

   subroutine run_round_trip(error, label, n_write, n_after, n_ice, resume_point_only)
      !! Straight N_WRITE + N_AFTER vs N_WRITE + checkpoint + warm restart
      !! through `engine_setup` + N_AFTER, compared on every field -- twice:
      !! at the RESUME POINT (the writer's state right after the checkpoint
      !! vs the reader's right after `engine_enter_data`, i.e. exactly what
      !! the first resumed step reads -- no setup pass may have re-derived
      !! any of it), and after the N_AFTER steps.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: label
      integer, intent(in) :: n_write, n_after
      integer, intent(out) :: n_ice
      logical, intent(in), optional :: resume_point_only
         !! `.true.` = only the AT-RESUME-POINT snapshot comparison runs
         !! (checkpoint plumbing); the post-restart N_AFTER steps are
         !! skipped, so a case whose answer is known to drift for a
         !! reason OTHER than the checkpoint/restore wiring can still
         !! assert the registry round-trip is exact without asserting
         !! something this test cannot (yet) make true.  No current case
         !! uses `.true.` -- the `restart_visc_rem` gap that used to need
         !! it is closed (see `test_engine_bit_exact_visc_rem`); kept as
         !! general test infra for a future compat gap of the same shape.
         !! Default `.false.` = the full round trip (every case today).
      type(ocean_engine_t), target :: ea, eb
      type(config_t) :: cfg_a, cfg_b
      type(snapshot_t) :: sa, sb, sa0, sb0
      real(wp) :: t_a, t_b
      integer :: step_b, ierr
      logical :: ok, point_only

      n_ice = 0
      point_only = .false.
      if (present(resume_point_only)) point_only = resume_point_only
      ! ---- A: n_write steps, checkpoint, n_after more ----
      call make_engine(ea, cfg_a, case_nml(label), ok)
      call check(error, ok, "engine A setup failed")
      if (allocated(error)) return
      t_a = 0.0_wp
      call advance(ea, cfg_a, t_a, n_write, ok)
      call check(error, ok, "engine A failed before the checkpoint")
      if (allocated(error)) return
      call ocean_state_restart_write(ea%state, ea%grid, ea%decomp, FN, t_a, n_write, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "checkpoint write failed")
      if (allocated(error)) return
      call take_snapshot(ea, sa0)
      call advance(ea, cfg_a, t_a, n_after, ok)
      call check(error, ok, "engine A failed after the checkpoint")
      if (allocated(error)) return
      call take_snapshot(ea, sa)
      call engine_exit_data(ea)
      call engine_teardown(ea)

      ! ---- B: resume from the checkpoint, n_after steps ----
      call make_engine(eb, cfg_b, case_nml(label), ok, restart_file=FN, t0=t_b, step0=step_b)
      call check(error, ok, "engine B (warm restart) setup failed")
      if (allocated(error)) return
      call check(error, step_b == n_write, "restart step count not restored")
      if (allocated(error)) return
      call check(error, eb%warm_restart, "engine B did not take the warm-restart path")
      if (allocated(error)) return
      call take_snapshot(eb, sb0)
      call compare(error, sa0, sb0, n_ice, trim(label)//" nw="//itoa(n_write)//": at the resume point")
      if (allocated(error) .or. point_only) then
         call engine_exit_data(eb)
         call engine_teardown(eb)
         call delete_file(FN)
         return
      end if
      call advance(eb, cfg_b, t_b, n_after, ok)
      call check(error, ok, "engine B failed after the restart")
      if (allocated(error)) return
      call take_snapshot(eb, sb)
      call engine_exit_data(eb)
      call engine_teardown(eb)
      call delete_file(FN)

      call compare(error, sa, sb, n_ice, trim(label)//" nw="//itoa(n_write)//": after the resumed steps")
   end subroutine run_round_trip

   subroutine test_engine_bit_exact(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: n_ice
      call run_round_trip(error, "island_periodic_zfixed", 6, 3, n_ice)
   end subroutine test_engine_bit_exact

   subroutine test_engine_bit_exact_visc_rem(error)
      !! `visc_rem_chain` live (D1 follow-up: was `correction_visc_rem`,
      !! retired 2026-10 -- `visc_rem_chain` gives the same NON-trivial
      !! `visc_rem_u/v` at checkpoint via the self-sufficient forcing/
      !! renorm/bt_rem_from producer, plus broader restart coverage of
      !! those three).  PR-2 (bt-rem-from-av-rem) closes
      !! compat row `restart_visc_rem` for real: checkpointing
      !! `bt_work%visc_rem_u/v` (PR-1) was necessary but NOT sufficient,
      !! because `visc_rem_precompute` (the pre-substep producer, gated
      !! on `is_pc .or. forcing_visc_rem .or. renorm_visc_rem`) builds its
      !! remnant matrix from the PREVIOUS stage's `vmix%kv` -- a one-stage
      !! lag, by design, mirroring MOM6's `vertvisc_coef` -- and `vmix%kv`
      !! itself was never checkpointed.  Registering `vmix_kv`
      !! (`ocean_state_build_restart_registry`,
      !! `src/core/ocean/state/ocean_state/rdb_ocean_state_restart.F90`) closes the gap: this
      !! now runs the FULL round trip (N_AFTER steps post-restart, not
      !! just the at-resume-point snapshot) and is bit-exact end to end.
      type(error_type), allocatable, intent(out) :: error
      integer :: n_ice
      call run_round_trip(error, "island_periodic_zfixed_visc_rem", 6, 3, n_ice)
   end subroutine test_engine_bit_exact_visc_rem

   subroutine test_engine_bit_exact_visc_rem_kshear_vertex(error)
      !! c069 restart fix: the visc_rem chain + kappa-shear VERTEX form
      !! together on `zstar` -- see `case_nml`'s
      !! `seamount_periodic_zstar_visc_rem_kshear_vertex` docstring for the
      !! mechanism. Fails before the `kd_corner` early-allocation +
      !! restart-registration fix (`rdb_ocean_engine.F90`,
      !! `ocean_state_build_restart_registry` in
      !! `ocean_state/rdb_ocean_state_restart.F90`); bit-exact end to end after it.
      type(error_type), allocatable, intent(out) :: error
      integer :: n_ice
      call run_round_trip(error, "seamount_periodic_zstar_visc_rem_kshear_vertex", 12, 12, n_ice)
   end subroutine test_engine_bit_exact_visc_rem_kshear_vertex

   subroutine test_engine_bit_exact_bkgnd_profile(error)
      !! PR-2 review fix: `vmix_seed_backgrounds`'s FULL round trip under
      !! a `kd_bg`-dependent config.  The first `reseed_arrays` version
      !! of the warm-restart fix returned before seeding `kt`/`ks`/
      !! `kd_bg` and the `kv`/`ks` boundary zero -- `kd_bg` is set ONLY
      !! in that routine, so a Bryan-Lewis background (`&ocean_vmix_nml
      !! bkgnd_profile`) silently lost its depth-varying floor on every
      !! warm restart of a `visc_rem_chain`/`implicit_drag` (closed
      !! south/north wall) case.  Full round trip (not resume-point-only)
      !! on `island_periodic_zfixed_bkgnd`.
      type(error_type), allocatable, intent(out) :: error
      integer :: n_ice
      call run_round_trip(error, "island_periodic_zfixed_bkgnd", 6, 3, n_ice)
   end subroutine test_engine_bit_exact_bkgnd_profile

   subroutine test_engine_bkgnd_full_reseed_on_restart(error)
      !! Narrower companion to `test_engine_bit_exact_bkgnd_profile`:
      !! asserts the array values THEMSELVES, right at the resume point
      !! (before any post-restart step runs `vmix_assemble`, which would
      !! otherwise recompute `kd_bg` from scratch every stage and mask a
      !! stale seed within one step -- confirmed empirically: the
      !! bitwise round-trip test above does NOT, by itself, distinguish
      !! a correct full reseed from the original bug, because
      !! `vmix_bkgnd_fill_impl` self-heals `kd_bg` on the very first
      !! `vmix_assemble` call either way).  `island_periodic_zfixed_bkgnd`
      !! sets `pp81_kappa_bg = 4.0e-5`, distinct from the 1.0e-5 type
      !! default `ocean_vmix_init`'s COLD seed would leave `kt`/`ks`/
      !! `kd_bg` at if the warm-restart config-copy wrongly skipped them
      !! (the original bug) -- so reading back the type default here
      !! instead of 4.0e-5 is the direct, trajectory-independent
      !! signature of that bug.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: ea, eb
      type(config_t) :: cfg_a, cfg_b
      real(wp) :: t_a, t_b
      integer :: step_b, ierr
      logical :: ok
      integer :: ip, jp

      checks: block
         call make_engine(ea, cfg_a, case_nml("island_periodic_zfixed_bkgnd"), ok)
         call check(error, ok, "engine A setup failed")
         if (allocated(error)) exit checks
         t_a = 0.0_wp
         call advance(ea, cfg_a, t_a, 6, ok)
         call check(error, ok, "engine A failed before the checkpoint")
         if (allocated(error)) exit checks
         call ocean_state_restart_write(ea%state, ea%grid, ea%decomp, FN, t_a, 6, ierr=ierr)
         call check(error, ierr == OCEAN_STATUS_OK, "checkpoint write failed")
         if (allocated(error)) exit checks
         call engine_exit_data(ea)
         call engine_teardown(ea)

         call make_engine(eb, cfg_b, case_nml("island_periodic_zfixed_bkgnd"), ok, &
                          restart_file=FN, t0=t_b, step0=step_b)
         call check(error, ok, "engine B (warm restart) setup failed")
         if (allocated(error)) exit checks
         call check(error, eb%state%vmix%kv_from_restart, &
                    "kv_from_restart must be .true. -- vmix_kv WAS in this checkpoint")
         if (allocated(error)) exit checks

         ip = ea%grid%nghost + 3; jp = ea%grid%nghost + 3
         call check(error, eb%state%vmix%kd_bg(ip, jp, 3) == 4.0e-5_wp, &
                    "kd_bg at the resume point was not fully reseeded to the "// &
                    "configured pp81_kappa_bg (4.0e-5) -- the original bug read "// &
                    "back the 1.0e-5 type default instead")
         if (allocated(error)) exit checks
         call check(error, eb%state%vmix%kt(ip, jp, 3) == 4.0e-5_wp, &
                    "kt at the resume point was not fully reseeded to the "// &
                    "configured pp81_kappa_bg")
         if (allocated(error)) exit checks
         call check(error, eb%state%vmix%ks(ip, jp, 3) == 4.0e-5_wp, &
                    "ks at the resume point was not fully reseeded to the "// &
                    "configured pp81_kappa_bg")
      end block checks
      call engine_exit_data(eb)
      call engine_teardown(eb)
      call delete_file(FN)
   end subroutine test_engine_bkgnd_full_reseed_on_restart

   subroutine test_engine_cold_seeds_kv_without_checkpoint(error)
      !! PR-2 review: a restart from a checkpoint that does NOT carry
      !! `vmix_kv` (simulating an older checkpoint schema, before
      !! `vmix_kv` was registered) must cold-seed `kv` from the
      !! configured `pp81_nu_bg` background -- not crash, not leave `kv`
      !! at whatever the allocator handed back, and `kv_from_restart`
      !! must read back `.false.` (the signal `configure_ocean_lateral`
      !! uses to decide whether to skip the reseed).  Uses
      !! `ocean_state_restart_write_drop_field` (test-only,
      !! `RDB_ENABLE_TESTING`) to write a checkpoint with `vmix_kv`
      !! stashed under a mangled tag, so the read genuinely cannot find
      !! it -- not a hand-set flag standing in for the file contents.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: ea, eb
      type(config_t) :: cfg_a, cfg_b
      real(wp) :: t_a, t_b
      integer :: step_b, ierr
      logical :: ok
      real(wp) :: kv_interior, kv_bed, expect_bg
      integer :: ip, jp

      checks: block
         call make_engine(ea, cfg_a, case_nml("island_periodic_zfixed_bkgnd"), ok)
         call check(error, ok, "engine A setup failed")
         if (allocated(error)) exit checks
         t_a = 0.0_wp
         call advance(ea, cfg_a, t_a, 6, ok)
         call check(error, ok, "engine A failed before the checkpoint")
         if (allocated(error)) exit checks
         call ocean_state_restart_write_drop_field(ea%state, ea%grid, ea%decomp, FN, t_a, &
                                                   6, "vmix_kv", ierr=ierr)
         call check(error, ierr == OCEAN_STATUS_OK, "checkpoint write (drop vmix_kv) failed")
         if (allocated(error)) exit checks
         call engine_exit_data(ea)
         call engine_teardown(ea)

         call make_engine(eb, cfg_b, case_nml("island_periodic_zfixed_bkgnd"), ok, &
                          restart_file=FN, t0=t_b, step0=step_b)
         call check(error, ok, "engine B (warm restart without vmix_kv) setup failed")
         if (allocated(error)) exit checks
         call check(error,.not. eb%state%vmix%kv_from_restart, &
                    "kv_from_restart must be .false. when vmix_kv was absent from the file")
         if (allocated(error)) exit checks

         ! kv must be the COLD-SEEDED background (pp81_nu_bg), interior
         ! AND the closed-BC boundary interfaces zeroed -- i.e. the full
         ! array reseed ran, not a half-applied one.
         ip = ea%grid%nghost + 3; jp = ea%grid%nghost + 3
         expect_bg = eb%state%vmix%pp81_nu_bg
         kv_interior = eb%state%vmix%kv(ip, jp, 3)
         kv_bed = eb%state%vmix%kv(ip, jp, 1)
         call check(error, kv_interior == expect_bg, &
                    "kv was not cold-seeded to pp81_nu_bg without the checkpoint field")
         if (allocated(error)) exit checks
         call check(error, kv_bed == 0.0_wp, &
                    "kv's closed-BC boundary interface was not zeroed on cold-seed")
         if (allocated(error)) exit checks

         ! "still runs": advance a few more steps with no error.
         call advance(eb, cfg_b, t_b, 3, ok)
         call check(error, ok, "engine B failed to advance after a cold-seeded kv resume")
      end block checks
      call engine_exit_data(eb)
      call engine_teardown(eb)
      call delete_file(FN)
   end subroutine test_engine_cold_seeds_kv_without_checkpoint

   subroutine test_engine_cavity_old_checkpoint_loads(error)
      !! OLD-CHECKPOINT path of the cavity restart fields.  The six
      !! `cavity_*` / `shelf_stress` entries are optional on read, so a
      !! checkpoint written before they were registered must still load
      !! (with the "Restart file lacks optional field" warning) and the
      !! absent field keeps its zero seed.  `drop_field` drops ONE tag per
      !! call; the one dropped here is `shelf_stress`, the one that feeds
      !! the melt-only under-ice `u_*` of the first resumed step (the case
      !! is the melt-only `cavity_sigma_ustar`, where that stress is
      !! non-zero).  Asserts: setup succeeds, the run is a warm restart,
      !! the dropped field is at its zero seed while a carried sibling
      !! (`cavity_ustar`) was restored, and it steps and stays finite.
      !! NOT bitwise by design (the first resumed step has no under-ice
      !! u_*); the bitwise gate is the `restart_engine_bit_exact_cavity_*`
      !! set.
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: ea, eb
      type(config_t) :: cfg_a, cfg_b
      type(snapshot_t) :: sb
      real(wp) :: t_a, t_b
      integer :: step_b, ierr, f
      logical :: ok, ea_live, eb_live

      ea_live = .false.
      eb_live = .false.

      checks: block
         call make_engine(ea, cfg_a, case_nml("cavity_sigma_ustar"), ok)
         ea_live = ok
         call check(error, ok, "engine A setup failed")
         if (allocated(error)) exit checks
         t_a = 0.0_wp
         call advance(ea, cfg_a, t_a, 5, ok)
         call check(error, ok, "engine A failed before the checkpoint")
         if (allocated(error)) exit checks
         !$acc update self(ea%state%surface_stress%stress_shelf)
         call check(error, maxval(abs(ea%state%surface_stress%stress_shelf)) > 0.0_wp, &
                    "vacuous: stress_shelf is zero in the written run")
         if (allocated(error)) exit checks
         call ocean_state_restart_write_drop_field(ea%state, ea%grid, ea%decomp, FN, t_a, &
                                                   5, "shelf_stress", ierr=ierr)
         call check(error, ierr == OCEAN_STATUS_OK, "checkpoint write (drop shelf_stress) failed")
         if (allocated(error)) exit checks
         call engine_exit_data(ea)
         call engine_teardown(ea)
         ea_live = .false.

         call make_engine(eb, cfg_b, case_nml("cavity_sigma_ustar"), ok, &
                          restart_file=FN, t0=t_b, step0=step_b)
         eb_live = ok
         call check(error, ok, "engine B (old-style cavity checkpoint) setup failed")
         if (allocated(error)) exit checks
         call check(error, eb%warm_restart, "engine B did not take the warm-restart path")
         if (allocated(error)) exit checks
         call check(error, step_b == 5, "restart step count not restored")
         if (allocated(error)) exit checks
         !$acc update self(eb%state%surface_stress%stress_shelf, eb%state%cavity_flux%ustar)
         call check(error, maxval(abs(eb%state%surface_stress%stress_shelf)) == 0.0_wp, &
                    "the dropped shelf_stress did not keep its zero seed")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(eb%state%cavity_flux%ustar)) > 0.0_wp, &
                    "the carried cavity_ustar was not restored")
         if (allocated(error)) exit checks
         call advance(eb, cfg_b, t_b, 4, ok)
         call check(error, ok, "engine B failed to step after the old-style resume")
         if (allocated(error)) exit checks
         call take_snapshot(eb, sb)
         do f = 1, sb%n
            call check(error, all(ieee_is_finite(sb%f(f)%a)), &
                       "non-finite "//trim(sb%f(f)%tag)//" after the old-style resume")
            if (allocated(error)) exit checks
         end do
      end block checks
      ! Tear down whatever is still mapped on EVERY path: a check that
      ! fails above `exit checks` must not leave an engine's device
      ! mappings (and the halo buffers) behind for the next case.
      if (ea_live) then
         call engine_exit_data(ea)
         call engine_teardown(ea)
      end if
      if (eb_live) then
         call engine_exit_data(eb)
         call engine_teardown(eb)
      end if
      call delete_file(FN)
   end subroutine test_engine_cavity_old_checkpoint_loads

   subroutine test_engine_bit_exact_ice(error)
      type(error_type), allocatable, intent(out) :: error
      call ice_round_trip(error, "sea_ice_periodic")
   end subroutine test_engine_bit_exact_ice

   subroutine test_engine_bit_exact_ice_ssp(error)
      type(error_type), allocatable, intent(out) :: error
      call ice_round_trip(error, "sea_ice_periodic_ssp_rk2")
   end subroutine test_engine_bit_exact_ice_ssp

   subroutine test_engine_bit_exact_ice_components(error)
      type(error_type), allocatable, intent(out) :: error
      call ice_round_trip(error, "sea_ice_components")
   end subroutine test_engine_bit_exact_ice_components

   subroutine test_engine_bit_exact_gm(error)
      !! N_WRITE = 5 (odd) with `dt_therm_ratio = 2`: the checkpoint lands
      !! mid thermo window, so the resumed step runs the GM operator on the
      !! carried slopes / KhTh, and the next thermo step's MEKE reads the
      !! carried `gm_src`.
      type(error_type), allocatable, intent(out) :: error
      integer :: n_ice
      call run_round_trip(error, "gm_meke_varmix_therm2", 5, 4, n_ice)
   end subroutine test_engine_bit_exact_gm

   subroutine test_engine_bit_exact_cavity_melt_drag(error)
      !! Basal melt + top drag, `dt_therm_ratio = 1`: the first resumed step
      !! must integrate the carried melt heat/salt flux.
      type(error_type), allocatable, intent(out) :: error
      integer :: n_ice
      call run_round_trip(error, "cavity_zfixed", 5, 4, n_ice)
   end subroutine test_engine_bit_exact_cavity_melt_drag

   subroutine test_engine_bit_exact_cavity_melt_only(error)
      !! Melt only (no top drag, KPP on): the under-ice `u_*` of the first
      !! resumed step comes from the melt slot's lagged `ustar`.
      type(error_type), allocatable, intent(out) :: error
      integer :: n_ice
      call run_round_trip(error, "cavity_sigma_melt_only", 5, 4, n_ice)
   end subroutine test_engine_bit_exact_cavity_melt_only

   subroutine test_engine_bit_exact_cavity_therm2_end(error)
      !! `dt_therm_ratio = 2`, checkpoint at a WINDOW BOUNDARY (even N_WRITE):
      !! the resumed step is a thermo step.
      type(error_type), allocatable, intent(out) :: error
      integer :: n_ice
      call run_round_trip(error, "cavity_zfixed_therm2", 4, 4, n_ice)
   end subroutine test_engine_bit_exact_cavity_therm2_end

   subroutine test_engine_bit_exact_cavity_therm2_mid(error)
      !! `dt_therm_ratio = 2`, checkpoint MID window (odd N_WRITE).
      type(error_type), allocatable, intent(out) :: error
      integer :: n_ice
      call run_round_trip(error, "cavity_zfixed_therm2", 5, 4, n_ice)
   end subroutine test_engine_bit_exact_cavity_therm2_mid

   subroutine test_engine_bit_exact_cavity_melt_only_therm2(error)
      !! Melt only + KPP at `dt_therm_ratio = 2`, checkpoint at a window
      !! boundary and mid window: the non-thermo resumed step reads the
      !! lagged melt `ustar` through `stress_shelf`.
      type(error_type), allocatable, intent(out) :: error
      integer :: n_ice
      call run_round_trip(error, "cavity_sigma_melt_only_therm2", 4, 4, n_ice)
      if (allocated(error)) return
      call run_round_trip(error, "cavity_sigma_melt_only_therm2", 5, 4, n_ice)
   end subroutine test_engine_bit_exact_cavity_melt_only_therm2

   subroutine test_engine_bit_exact_cavity_mass(error)
      !! `freshwater = 'mass'`: the in-stage mass step of the first resumed
      !! thermo step reads the previous step's `melt`/`active`/`s_far`/`t_b`.
      !! Ratio 1, then ratio 2 at a window boundary and mid window.
      type(error_type), allocatable, intent(out) :: error
      integer :: n_ice
      call run_round_trip(error, "cavity_zfixed_mass", 5, 4, n_ice)
      if (allocated(error)) return
      call run_round_trip(error, "cavity_zfixed_mass_therm2", 4, 4, n_ice)
      if (allocated(error)) return
      call run_round_trip(error, "cavity_zfixed_mass_therm2", 5, 4, n_ice)
   end subroutine test_engine_bit_exact_cavity_mass

   subroutine test_engine_bit_exact_cavity_mass_uo(error)
      !! `freshwater = 'mass'` + `volume_compensation = 'uniform_open_ocean'`.
      type(error_type), allocatable, intent(out) :: error
      integer :: n_ice
      call run_round_trip(error, "cavity_zfixed_mass_uo_therm2", 4, 4, n_ice)
      if (allocated(error)) return
      call run_round_trip(error, "cavity_zfixed_mass_uo_therm2", 5, 4, n_ice)
   end subroutine test_engine_bit_exact_cavity_mass_uo

   subroutine test_engine_bit_exact_cavity_ustar(error)
      !! Melt-only under-ice `u_*` made REAL: sigma, 16 thin layers,
      !! `u_tide = 0.2`, so the under-ice KPP boundary layer reaches below
      !! the top layer and `stress_shelf` (lagged melt `ustar`) matters.
      type(error_type), allocatable, intent(out) :: error
      integer :: n_ice
      call run_round_trip(error, "cavity_sigma_ustar", 5, 4, n_ice)
      if (allocated(error)) return
      call run_round_trip(error, "cavity_sigma_ustar_therm2", 4, 4, n_ice)
      if (allocated(error)) return
      call run_round_trip(error, "cavity_sigma_ustar_therm2", 5, 4, n_ice)
   end subroutine test_engine_bit_exact_cavity_ustar

   subroutine test_engine_cavity_ustar_not_vacuous(error)
      !! Non-vacuity guard of the `cavity_sigma_ustar` cases: (1) under-ice
      !! KPP is live (boundary-layer depth exceeds the top-layer thickness on
      !! a covered column) and the lagged `stress_shelf` is non-zero; (2)
      !! zeroing `stress_shelf` ONCE in an otherwise straight run (exactly
      !! what a resume without it does) changes the answer.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: ea, ec
      type(config_t) :: cfg_a, cfg_c
      type(snapshot_t) :: sa, sc
      real(wp) :: t_a, t_c, ratio, smax
      integer :: i, j, e
      logical :: ok, differs

      call make_engine(ea, cfg_a, case_nml("cavity_sigma_ustar"), ok)
      call check(error, ok, "engine A setup failed")
      if (allocated(error)) return
      call make_engine(ec, cfg_c, case_nml("cavity_sigma_ustar"), ok)
      call check(error, ok, "engine C setup failed")
      if (allocated(error)) return
      t_a = 0.0_wp
      t_c = 0.0_wp
      call advance(ea, cfg_a, t_a, 5, ok)
      call advance(ec, cfg_c, t_c, 5, ok)
      call check(error, ok, "advance to the perturbation point failed")
      if (allocated(error)) return

      !$acc update self(ea%state%vmix%bl_depth, ea%state%multilayer%h_layer)
      !$acc update self(ea%state%surface_stress%stress_shelf)
      ratio = 0.0_wp
      smax = maxval(ea%state%surface_stress%stress_shelf)
      do j = 1, ea%grid%ny_total
         do i = 1, ea%grid%nx_total
            if (ea%state%metrics%cover_frac(i, j) > 0.5_wp .and. &
                ea%state%multilayer%h_layer(i, j, ea%state%multilayer%nz_ml) > 0.0_wp) then
               ratio = max(ratio, ea%state%vmix%bl_depth(i, j)/ &
                           ea%state%multilayer%h_layer(i, j, ea%state%multilayer%nz_ml))
            end if
         end do
      end do
      write (*, '(a,es12.4,a,es12.4)') "ustar guard: max(bl_depth/h_top) under ice = ", &
         ratio, ";  max stress_shelf = ", smax
      call check(error, ratio > 1.0_wp, "under-ice KPP boundary layer does not exceed the top layer")
      if (allocated(error)) return
      call check(error, smax > 0.0_wp, "stress_shelf is zero in the straight run")
      if (allocated(error)) return

      ec%state%surface_stress%stress_shelf = 0.0_wp
      !$acc update device(ec%state%surface_stress%stress_shelf)
      call advance(ea, cfg_a, t_a, 4, ok)
      call advance(ec, cfg_c, t_c, 4, ok)
      call check(error, ok, "advance after the perturbation failed")
      if (allocated(error)) return
      call take_snapshot(ea, sa)
      call take_snapshot(ec, sc)
      differs = .false.
      do e = 1, sa%n
         if (any(shape(sa%f(e)%a) /= shape(sc%f(e)%a))) cycle
         if (any(transfer(sa%f(e)%a, [0_int64]) /= transfer(sc%f(e)%a, [0_int64]))) then
            differs = .true.
         end if
      end do
      call check(error, differs, "zeroing stress_shelf did not change the answer: the case is vacuous")
      call engine_exit_data(ea)
      call engine_teardown(ea)
      call engine_exit_data(ec)
      call engine_teardown(ec)
   end subroutine test_engine_cavity_ustar_not_vacuous

   subroutine test_engine_bit_exact_cavity_ssp(error)
      !! `ssp_rk2` applies the surface flux in both stages.
      type(error_type), allocatable, intent(out) :: error
      integer :: n_ice
      call run_round_trip(error, "cavity_sigma_melt_only_ssp", 5, 4, n_ice)
      if (allocated(error)) return
      call run_round_trip(error, "cavity_sigma_ustar_therm2_ssp", 5, 4, n_ice)
      if (allocated(error)) return
      call run_round_trip(error, "cavity_zfixed_therm2_ssp", 5, 4, n_ice)
      if (allocated(error)) return
      call run_round_trip(error, "cavity_zfixed_mass_therm2_ssp", 5, 4, n_ice)
   end subroutine test_engine_bit_exact_cavity_ssp

   subroutine test_engine_bit_exact_cavity_tdrag_sigma(error)
      !! Top drag on SIGMA (KPP live), so the axes are not confounded.
      type(error_type), allocatable, intent(out) :: error
      integer :: n_ice
      call run_round_trip(error, "cavity_sigma_tdrag_ustar", 5, 4, n_ice)
      if (allocated(error)) return
      call run_round_trip(error, "cavity_sigma_tdrag_ustar_therm2", 5, 4, n_ice)
   end subroutine test_engine_bit_exact_cavity_tdrag_sigma

   function atm_group(list) result(g)
      !! `&ocean_atm_forcing_nml` pointing all nine fields at one filelist.
      character(len=*), intent(in) :: list
      character(len=:), allocatable :: g
      g = "&ocean_atm_forcing_nml enable = .true., time_mode = 'interannual', "// &
          "sw_interp = 'linear', filelist_tas = '"//list//"', filelist_huss = '"//list// &
          "', filelist_uas = '"//list//"', filelist_vas = '"//list//"', filelist_psl = '"// &
          list//"', filelist_rlds = '"//list//"', filelist_rsds = '"//list// &
          "', filelist_prra = '"//list//"', filelist_prsn = '"//list//"' /"//NL
   end function atm_group

   subroutine write_atm_file(path, listname, tas0, huss0, rsds_const)
      !! One global lon-lat file (10 degrees, two time records, all nine CMOR
      !! names) of spatially varying prescribed atmospheric state, plus the
      !! one-line filelist the namelist points at.  `rsds_const >= 0` makes
      !! the downwelling shortwave spatially/temporally UNIFORM at that value
      !! (the shortwave probe); negative keeps the varying field.
      character(len=*), intent(in) :: path, listname
      real(wp), intent(in) :: tas0, huss0, rsds_const
      integer, parameter :: NLON = 36, NLAT = 18
      real(wp), parameter :: PI = 3.14159265358979323846_wp
      character(len=4), parameter :: VN(9) = [character(len=4) :: "tas", "huss", "uas", &
                                                                                       "vas", "psl", "rlds", "rsds", "prra", "prsn"]
      real(wp) :: lon(NLON), lat(NLAT), lonb(2, NLON), latb(2, NLAT), tax(2), f(NLON, NLAT, 2)
      real(wp) :: sx, sy, dlo, dla
      integer :: ncid, dlon_, dlat_, dbnd, dtim, vlon, vlat, vlonb, vlatb, vtim, vf(9)
      integer :: v, i, j, k, ierr, u

      dlo = 360.0_wp/real(NLON, wp)
      dla = 180.0_wp/real(NLAT, wp)
      do i = 1, NLON
         lonb(1, i) = dlo*real(i - 1, wp)
         lonb(2, i) = dlo*real(i, wp)
         lon(i) = 0.5_wp*(lonb(1, i) + lonb(2, i))
      end do
      do j = 1, NLAT
         latb(1, j) = -90.0_wp + dla*real(j - 1, wp)
         latb(2, j) = -90.0_wp + dla*real(j, wp)
         lat(j) = 0.5_wp*(latb(1, j) + latb(2, j))
      end do
      tax = [0.0_wp, 10.0_wp]
      call nc_create_file(path, ncid)
      call nc_def_dim(ncid, "lon", NLON, dlon_)
      call nc_def_dim(ncid, "lat", NLAT, dlat_)
      call nc_def_dim(ncid, "bnds", 2, dbnd)
      call nc_def_dim(ncid, "time", 2, dtim)
      call rdb_def_var_1d(ncid, "lon", dlon_, vlon)
      call rdb_def_var_1d(ncid, "lat", dlat_, vlat)
      call nc_def_var_2d(ncid, "lon_bnds", [dbnd, dlon_], vlonb)
      call nc_def_var_2d(ncid, "lat_bnds", [dbnd, dlat_], vlatb)
      call rdb_def_var_1d(ncid, "time", dtim, vtim)
      call nc_put_att(ncid, vtim, "units", "days since 2000-06-20 00:00:00")
      do v = 1, 9
         call nc_def_var_3d(ncid, trim(VN(v)), [dlon_, dlat_, dtim], vf(v))
      end do
      call nc_enddef(ncid)
      call rdb_put_var_1d(ncid, vlon, lon)
      call rdb_put_var_1d(ncid, vlat, lat)
      ierr = nf90_put_var(ncid, vlonb, lonb)
      ierr = nf90_put_var(ncid, vlatb, latb)
      call rdb_put_var_1d(ncid, vtim, tax)
      do v = 1, 9
         do k = 1, 2
            do j = 1, NLAT
               do i = 1, NLON
                  sx = sin(2.0_wp*PI*lon(i)/37.0_wp + 0.3_wp)
                  sy = cos(2.0_wp*PI*lat(j)/29.0_wp - 0.5_wp)
                  select case (v)
                  case (1)
                     f(i, j, k) = tas0 + 4.0_wp*sx + 0.15_wp*(lat(j) - 20.0_wp) + 0.05_wp*tax(k)
                  case (2)
                     f(i, j, k) = huss0*(1.0_wp + 0.2_wp*sy + 0.1_wp*sx)
                  case (3)
                     f(i, j, k) = 6.0_wp + 4.0_wp*sx + 0.05_wp*(lat(j) - 20.0_wp)
                  case (4)
                     f(i, j, k) = -3.0_wp + 3.0_wp*sy + 0.04_wp*lon(i) + 0.02_wp*tax(k)
                  case (5)
                     f(i, j, k) = 1.01e5_wp - 300.0_wp*sx + 5.0_wp*lat(j)
                  case (6)
                     f(i, j, k) = 330.0_wp + 15.0_wp*sy + 20.0_wp*sx
                  case (7)
                     if (rsds_const >= 0.0_wp) then
                        f(i, j, k) = rsds_const
                     else
                        f(i, j, k) = 180.0_wp + 40.0_wp*sx*sy + 0.5_wp*lat(j)
                     end if
                  case (8)
                     f(i, j, k) = 2.0e-5_wp*(1.0_wp + 0.5_wp*sx)
                  case default
                     f(i, j, k) = 1.0e-6_wp*(1.0_wp + 0.5_wp*sy)
                  end select
               end do
            end do
         end do
         ierr = nf90_put_var(ncid, vf(v), f)
      end do
      call nc_close(ncid)
      open (newunit=u, file=listname, status='replace', action='write')
      write (u, '(a)') trim(path)
      close (u)
   end subroutine write_atm_file

   subroutine test_engine_bit_exact_ice_components_bulk(error)
      !! Restart round trip of the ice + components + bulk-forcing path,
      !! checkpoint mid thermo window (N_WRITE odd, `dt_therm_ratio = 2`).
      type(error_type), allocatable, intent(out) :: error
      integer :: n_ice
      call write_atm_file(ATM_FILE, ATM_LIST, 255.0_wp, 6.0e-4_wp, -1.0_wp)
      call ice_round_trip(error, "sea_ice_components_bulk")
      call delete_file(ATM_FILE)
      call delete_file(ATM_LIST)
      n_ice = 0
   end subroutine test_engine_bit_exact_ice_components_bulk

   function probe_nml(ratio, ice_on, list) result(nml)
      !! The 6x6 spherical flat-100 m, 4-layer, 15 degC warm basin: the
      !! atmosphere-driven bulk forcing + components, ice enabled with NO ice
      !! present (default `conc_config`) or disabled.
      integer, intent(in) :: ratio
      logical, intent(in) :: ice_on
      character(len=*), intent(in) :: list
      character(len=:), allocatable :: nml
      character(len=2) :: sr
      write (sr, '(i0)') ratio
      nml = "&sim_nml sim_type = 'ocean' /"//NL// &
            "&grid_nml nx = 6, ny = 6, nghost = 3, dx = 1.0, dy = 1.0 /"//NL// &
            "&time_nml t_end = 4800.0, dt_fixed = 600.0, start_date = '2000-06-21' /"//NL// &
            "&nonhydrostatic_nml nz_layers = 4 /"//NL// &
            "&ocean_bt_nml auto_n_inner = .true., split_scheme = 'pred_corr' /"//NL// &
            "&ocean_diag_nml enabled = .false. /"//NL// &
            "&ocean_eos_nml eos = 'wright' /"//NL// &
            "&ocean_thermo_nml enable_thermodynamics = .true. /"//NL// &
            "&ocean_forcing_nml enable_components = .true. /"//NL// &
            "&ocean_bulk_flux_nml enable = .true., n_iter = 5, z_wind = 10.0, z_ta = 10.0 /"//NL// &
            "&tracer_nml initial_temperature = 15.0, initial_salinity = 35.0 /"//NL// &
            atm_group(list)// &
            "&ocean_vmix_nml dt_therm_ratio = "//trim(sr)//" /"//NL// &
            "&ocean_ice_nml enable = "//trim(merge(".true. ", ".false.", ice_on))//" /"//NL// &
            "&ocean_grid_nml grid_config = 'spherical', lon_west = 0.0, lat_south = 20.0, "// &
            "rad_earth = 6.371e6, coriolis_scheme = 'planetary' /"//NL// &
            "&vcoord_nml vcoord_type = 'sigma' /"//NL// &
            "&ocean_topo_nml topo_config = 'flat', max_depth = 100.0 /"//NL// &
            "&ocean_bc_nml west = 'wall', east = 'wall', south = 'wall', north = 'wall' /"//NL// &
            "&output_nml output_to_file = .false. /"//NL
   end function probe_nml

   subroutine run_probe(ratio, ice_on, list, nsteps, t_top, ok)
      !! Run the probe `nsteps` steps; return the surface-layer temperature.
      integer, intent(in) :: ratio, nsteps
      logical, intent(in) :: ice_on
      character(len=*), intent(in) :: list
      real(wp), allocatable, intent(out) :: t_top(:, :)
      logical, intent(out) :: ok
      type(ocean_engine_t), target :: e
      type(config_t) :: cfg
      real(wp) :: t
      integer :: nz, it

      call make_engine(e, cfg, probe_nml(ratio, ice_on, list), ok)
      if (.not. ok) return
      t = 0.0_wp
      call advance(e, cfg, t, nsteps, ok)
      if (ok) then
         nz = e%state%multilayer%nz_ml
         it = e%state%multilayer%idx_temperature
         !$acc update self(e%state%multilayer%tracers(it)%c)
         t_top = e%state%multilayer%tracers(it)%c(:, :, nz)
      end if
      call engine_exit_data(e)
      call engine_teardown(e)
   end subroutine run_probe

   subroutine test_engine_icefree_sw(error)
      !! The regression for "ice enabled => the ocean receives no shortwave".
      !! Surface-layer warming caused by the shortwave (T at rsds = 400 minus
      !! T at rsds = 0) must be the SAME with ice enabled (no ice present) as
      !! with ice disabled, at `dt_therm_ratio` 1 and 3 -- and nonzero.
      !! Before the fix the ice-enabled warming was exactly 0.
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: t4_on(:, :), t0_on(:, :), t4_off(:, :), t0_off(:, :)
      real(wp) :: d_on, d_off
      integer :: ratios(2), k, ip, jp
      logical :: ok

      call write_atm_file(ATM_FILE, ATM_LIST, 288.0_wp, 6.0e-3_wp, 400.0_wp)
      call write_atm_file(ATM_SW0_FILE, ATM_SW0_LIST, 288.0_wp, 6.0e-3_wp, 0.0_wp)
      ratios = [1, 3]
      do k = 1, 2
         call run_probe(ratios(k), .true., ATM_LIST, 6, t4_on, ok)
         call check(error, ok, "probe (ice on, sw 400) failed")
         if (allocated(error)) exit
         call run_probe(ratios(k), .true., ATM_SW0_LIST, 6, t0_on, ok)
         call check(error, ok, "probe (ice on, sw 0) failed")
         if (allocated(error)) exit
         call run_probe(ratios(k), .false., ATM_LIST, 6, t4_off, ok)
         call check(error, ok, "probe (ice off, sw 400) failed")
         if (allocated(error)) exit
         call run_probe(ratios(k), .false., ATM_SW0_LIST, 6, t0_off, ok)
         call check(error, ok, "probe (ice off, sw 0) failed")
         if (allocated(error)) exit
         ip = 3 + 3
         jp = 3 + 3
         d_on = t4_on(ip, jp) - t0_on(ip, jp)
         d_off = t4_off(ip, jp) - t0_off(ip, jp)
         call check(error, d_off > 1.0e-6_wp, "precondition: shortwave warms the ice-disabled run")
         if (allocated(error)) exit
         call check(error, abs(d_on - d_off) <= 1.0e-9_wp*d_off, &
                    "ice-enabled shortwave warming equals the ice-disabled value (ratio "// &
                    itoa(ratios(k))//")")
         if (allocated(error)) exit
      end do
      call delete_file(ATM_FILE)
      call delete_file(ATM_LIST)
      call delete_file(ATM_SW0_FILE)
      call delete_file(ATM_SW0_LIST)
   end subroutine test_engine_icefree_sw

   subroutine ice_round_trip(error, label)
      !! N_WRITE = 5 (odd) with `dt_therm_ratio = 2`: the checkpoint lands
      !! mid thermo window; N_AFTER = 4 then closes two more windows.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: label
      integer :: n_ice
      call run_round_trip(error, label, 5, 4, n_ice)
      if (allocated(error)) return
      call check(error, n_ice >= 20, "the ice registry fields were not compared")
   end subroutine ice_round_trip

   subroutine delete_file(fname)
      character(len=*), intent(in) :: fname
      integer :: u, ios
      logical :: exists
      inquire (file=fname, exist=exists)
      if (.not. exists) return
      open (newunit=u, file=fname, status="old", iostat=ios)
      if (ios == 0) close (u, status="delete")
   end subroutine delete_file

end module test_ocean_restart_engine
