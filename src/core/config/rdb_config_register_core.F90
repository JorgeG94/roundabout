submodule(rdb_config) rdb_config_register_core
   !! Strict-schema registration: dispatcher + core/time/grid/misc groups.
   implicit none

contains

   ! ==================================================================
   ! Strict-schema construction (rdb_nml_schema).  Lives here (not in
   ! rdb_config_schema) because read_config builds + parses the schema
   ! itself — a module-level `use rdb_config_schema` from rdb_config
   ! would be a circular dependency.  rdb_config_schema re-exports
   ! build_rdb_schema so the production entry points keep their
   ! original import.
   ! ==================================================================

   module subroutine build_rdb_schema(cfg, schema)
      !! Register the validated groups + the still-external groups onto
      !! `schema`, capturing defaults from `cfg`.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(out) :: schema

      call register_sim(cfg, schema)
      call register_grid(cfg, schema)
      call register_time(cfg, schema)
      call register_mpi(cfg, schema)
      call register_logging(cfg, schema)
      call register_vcoord(cfg, schema)
      call register_physics(cfg, schema)
      call register_output(cfg, schema)
      call register_boundary(cfg, schema)
      call register_tracer(cfg, schema)
      call register_initial_condition(cfg, schema)
      call register_nonhydrostatic(cfg, schema)
      call register_kappa_shear(cfg, schema)
      call register_slopes(cfg, schema)
      call register_gm(cfg, schema)
      call register_redi(cfg, schema)
      call register_varmix(cfg, schema)
      call register_meke(cfg, schema)
      call register_tidal_mixing(cfg, schema)
      call register_conv(cfg, schema)
      call register_porous(cfg, schema)
      call register_ddiff(cfg, schema)
      call register_ocean_tides(cfg, schema)
      call register_ocean_psurf(cfg, schema)
      call register_ocean_cavity_dyn(cfg, schema)
      call register_ocean_cavity_melt(cfg, schema)
      call register_epbl(cfg, schema)
      call register_wavespeed(cfg, schema)
      call register_foxkemper(cfg, schema)
      call register_ocean_grid(cfg, schema)
      call register_ocean_coriolis(cfg, schema)
      call register_ocean_thermo(cfg, schema)
      call register_ocean_forcing(cfg, schema)
      call register_ocean_ice(cfg, schema)
      call register_ocean_ice_ic(cfg, schema)
      call register_ocean_restore(cfg, schema)
      call register_ocean_geothermal(cfg, schema)
      call register_ocean_sponge(cfg, schema)
      call register_ocean_tracers(cfg, schema)
      call register_ocean_bt(cfg, schema)
      call register_ocean_debug(cfg, schema)
      call register_ocean_mpi(cfg, schema)
      call register_ocean_wetdry(cfg, schema)
      call register_ocean_pgf(cfg, schema)
      call register_ocean_eos(cfg, schema)
      call register_ocean_bdrag(cfg, schema)
      call register_ocean_tdrag(cfg, schema)
      call register_ocean_hdiff(cfg, schema)
      call register_ocean_hvisc(cfg, schema)
      call register_ocean_vmix(cfg, schema)
      call register_ocean_vdiff(cfg, schema)
      call register_ocean_continuity(cfg, schema)
      call register_ocean_isopycnal(cfg, schema)
      call register_ocean_topo(cfg, schema)
      call register_ocean_ic(cfg, schema)
      call register_ocean_zinit(cfg, schema)
      call register_ocean_data(cfg, schema)
      call register_ocean_dataovr(cfg, schema)
      call register_ocean_diag(cfg, schema)
      call register_ocean_bc(cfg, schema)
   end subroutine build_rdb_schema

   module subroutine register_sim(cfg, schema)
      !! `&sim_nml`: the simulation regime selector.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      character(len=:), pointer :: ps

      g%name = "sim"
      g%doc = "Simulation-regime selector."
      ps => cfg%sim_type
      call g%add(nml_enum("sim_type", ps, &
                          "Simulation regime (ocean is the only regime this build ships)", &
                          allowed=[character(len=8) :: "ocean"]))
      call schema%add_group(g)
   end subroutine register_sim

   module subroutine register_grid(cfg, schema)
      !! `&grid_nml`: structured-grid geometry.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      integer, pointer :: pi
      real(wp), pointer :: pr
      character(len=:), pointer :: ps

      g%name = "grid"
      g%doc = "Structured-grid geometry."
      pi => cfg%nx
      call g%add(nml_int("nx", pi, "Number of physical cells in x", min=1))
      pi => cfg%ny
      call g%add(nml_int("ny", pi, "Number of physical cells in y", min=1))
      pr => cfg%dx
      call g%add(nml_real("dx", pr, "Cell size in x", units="m"))
      pr => cfg%dy
      call g%add(nml_real("dy", pr, "Cell size in y", units="m"))
      pi => cfg%nghost
      call g%add(nml_int("nghost", pi, "Ghost cells on each side", min=1))
      call schema%add_group(g)
   end subroutine register_grid

   module subroutine register_time(cfg, schema)
      !! `&time_nml`: time-integration controls.  `t_end` is interpreted
      !! in `time_unit` and converted to seconds by the post-parse
      !! cascade in read_config; `dt_fixed`/`dt_max` stay in seconds.
      !! `time_unit` is a plain string (the cascade validates its set).
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      integer, pointer :: pi
      real(wp), pointer :: pr
      character(len=:), pointer :: ps

      g%name = "time"
      g%doc = "Time-integration controls."
      pr => cfg%t_end
      call g%add(nml_real("t_end", pr, "Simulation end time (in time_unit)"))
      pr => cfg%cfl
      call g%add(nml_real("cfl", pr, "CFL number for adaptive timestep"))
      pr => cfg%dt_max
      call g%add(nml_real("dt_max", pr, "Maximum allowable timestep", units="s", &
                          dead_on_ocean_path="accepted and type/range-validated but read "// &
                          "nowhere in src/ -- the ocean path's adaptive-timestep ceiling "// &
                          "comes from the barotropic gravity-wave CFL (auto_n_inner), not "// &
                          "this knob (found by the P4 dead-knob sweep, 2026-09-10)."))
      pr => cfg%dt_fixed
      call g%add(nml_real("dt_fixed", pr, "Fixed timestep (0 = adaptive CFL)", units="s"))
      pi => cfg%cfl_interval
      call g%add(nml_int("cfl_interval", pi, "Recompute CFL timestep every N steps", min=1, &
                         dead_on_ocean_path="accepted and range-validated but read nowhere "// &
                         "in src/ -- there is no adaptive-CFL recompute cadence on the ocean "// &
                         "path to interval-gate (found by the P4 dead-knob sweep, 2026-09-10)."))
      ps => cfg%time_unit
      call g%add(nml_string("time_unit", ps, &
                            "Unit for the long-time fields: s/min/hr/day/year"))
      ps => cfg%start_date
      call g%add(nml_string("start_date", ps, &
                            "Proleptic-Gregorian date t=0 corresponds to "// &
                            "('' = plain relative clock, bit-identical)"))
      ps => cfg%calendar
      call g%add(nml_enum("calendar", ps, "Calendar backing start_date", &
                          allowed=[character(len=19) :: "gregorian", "standard", &
                                   "proleptic_gregorian"]))
      call schema%add_group(g)
   end subroutine register_time

   module subroutine register_mpi(cfg, schema)
      !! `&mpi_nml`: MPI domain decomposition.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      integer, pointer :: pi

      g%name = "mpi"
      g%doc = "MPI domain decomposition."
      pi => cfg%px
      call g%add(nml_int("px", pi, "MPI processes in x", min=1))
      pi => cfg%py
      call g%add(nml_int("py", pi, "MPI processes in y", min=1))
      call schema%add_group(g)
   end subroutine register_mpi

   module subroutine register_logging(cfg, schema)
      !! `&logging_nml`: logger verbosity + status cadence.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      real(wp), pointer :: pr
      character(len=:), pointer :: ps

      g%name = "logging"
      g%doc = "Logger verbosity + status cadence."
      ps => cfg%log_level
      call g%add(nml_enum("log_level", ps, "Log verbosity", &
                          allowed=[character(len=11) :: "debug", "verbose", "info", &
                                   "performance", "warning", "error"]))
      pr => cfg%status_interval
      call g%add(nml_real("status_interval", pr, &
                          "Status-print cadence (in time_unit; 0 = every 100 steps)"))
      call schema%add_group(g)
   end subroutine register_logging

   module subroutine register_vcoord(cfg, schema)
      !! `&vcoord_nml`: vertical-coordinate + ALE-remap controls.
      !! `vcoord_type` / `remap_method` / `zstar_stretching` are enums
      !! with the canonical sets accepted by the rdb_vcoord parsers.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      integer, pointer :: pi
      real(wp), pointer :: pr
      real(wp), pointer :: pra(:)
      character(len=:), pointer :: ps
      logical, pointer :: pl

      g%name = "vcoord"
      g%doc = "Vertical-coordinate + ALE-remap controls."
      ps => cfg%vcoord_type
      call g%add(nml_enum("vcoord_type", ps, "Vertical coordinate type", &
                          allowed=[character(len=12) :: "sigma", "zsigma", &
                                   "zstar", "zstar_full", "zstar_sigma", "z", &
                                   "eulerian_z", "isopycnal", "lagrangian", "gprime", &
                                   "z_fixed", "rho", "hycom"]))
      ps => cfg%thickness_config
      call g%add(nml_enum("thickness_config", ps, &
                          "Initial layer-thickness profile (ocean path)", &
                          allowed=[character(len=9) :: "sigma", "uniform_z"]))
      ps => cfg%remap_method
      call g%add(nml_enum("remap_method", ps, "Vertical remapping method", &
                          allowed=[character(len=6) :: "pcm", "plm", "ppm", "ppm_h4", "pqm"]))
      pr => cfg%zstar_h_surf_target
      call g%add(nml_real("zstar_h_surf_target", pr, &
                          "z*-full: target surface-layer thickness (0 = auto)", units="m"))
      ps => cfg%zstar_stretching
      call g%add(nml_enum("zstar_stretching", ps, &
                          "z*-full: surface-concentration stretching", &
                          allowed=[character(len=7) :: "log", "uniform"], &
                          dead_on_ocean_path="registered but never copied onto "// &
                          "ocean_state%vcoord -- the ocean path always runs "// &
                          "STRETCH_UNIFORM regardless of this setting "// &
                          "(rdb_ocean_vcoord.F90:179). See docs/ocean_python_api_plan.md P4.5."))
      pr => cfg%zstar_h_min
      call g%add(nml_real("zstar_h_min", pr, "z*-full: vanishing-layer floor", units="m"))
      pi => cfg%zstar_n_surf
      call g%add(nml_int("zstar_n_surf", pi, &
                         "z*-full: number of fine near-surface layers (0 = auto)", min=0, &
                         dead_on_ocean_path="registered but never copied onto "// &
                         "ocean_state%vcoord -- the ocean path always runs n_surf=0 "// &
                         "regardless of this setting (rdb_ocean_vcoord.F90:182). "// &
                         "See docs/ocean_python_api_plan.md P4.5."))
      pr => cfg%rho_ref_pressure
      call g%add(nml_real("rho_ref_pressure", pr, &
                          "rho-coord: reference pressure for potential density", units="Pa"))
      pr => cfg%rho_target_light
      call g%add(nml_real("rho_target_light", pr, &
                          "rho-coord: lightest (surface) target density", units="kg/m^3"))
      pr => cfg%rho_target_dense
      call g%add(nml_real("rho_target_dense", pr, &
                          "rho-coord: densest (bed) target density", units="kg/m^3"))
      ps => cfg%rho_target_profile
      call g%add(nml_enum("rho_target_profile", ps, &
                          "rho/hycom target densities: uniform (rho_target_light.."// &
                          "rho_target_dense linspace) or list (rho_target_list)", &
                          allowed=[character(len=7) :: "uniform", "list"]))
      pra => cfg%rho_target_list
      call g%add(nml_real_array("rho_target_list", pra, &
                                "rho_target_profile='list': interface densities, "// &
                                "lightest first (exactly nz_layers+1 entries)", &
                                units="kg/m^3"))
      pr => cfg%regrid_time_scale
      call g%add(nml_real("regrid_time_scale", pr, &
                          "ALE regrid grid time-filter timescale (0 = jump to target)", &
                          units="s", min=0.0_wp))
      pl => cfg%remap_vel_conserve_ke
      call g%add(nml_logical("remap_vel_conserve_ke", pl, &
                             "ALE velocity remap: KE-conserving baroclinic-anomaly rescale"))
      pl => cfg%remap_boundary_extrap
      call g%add(nml_logical("remap_boundary_extrap", pl, &
                             "ALE remap: linear-exact one-sided reconstruction in the "// &
                             "k=1/k=nz boundary cells (MOM6 BOUNDARY_EXTRAPOLATION)"))
      pl => cfg%remap_nonuniform_weights
      call g%add(nml_logical("remap_nonuniform_weights", pl, &
                             "ALE remap: non-uniform-grid PLM slope + PPM edge weights "// &
                             "(Colella-Woodward 1984 eqs 1.6-1.8) instead of the "// &
                             "equal-thickness specialisations"))
      pl => cfg%remap_check_preconditions
      call g%add(nml_logical("remap_check_preconditions", pl, &
                             "ALE remap: fail loud when a column violates the overlap "// &
                             "sweep's preconditions (non-negative thicknesses, "// &
                             "matching column totals)"))
      pl => cfg%check_vanished_content
      call g%add(nml_logical("check_vanished_content", pl, &
                             "I1' tripwire: fail loud if any layer at or below "// &
                             "H_VANISHED does not hold its donor live layer's "// &
                             "concentration (debug/validation)"))
      pl => cfg%zfixed_closed_faces
      call g%add(nml_logical("zfixed_closed_faces", pl, &
                             "z_fixed / zstar / zstar_full partial steps: close every face whose "// &
                             "layer is an inert filler on either side (z-level wall, "// &
                             "free-slip)"))
      ps => cfg%z_fixed_profile
      call g%add(nml_enum("z_fixed_profile", ps, &
                          "z_fixed / zstar levels / hycom z* floor nominal "// &
                          "layer-thickness profile: uniform (max_depth/nz), list "// &
                          "(z_fixed_dz) or tanh stretching", &
                          allowed=[character(len=7) :: "uniform", "list", "tanh"]))
      pra => cfg%z_fixed_dz
      call g%add(nml_real_array("z_fixed_dz", pra, &
                                "z_fixed_profile='list': nominal layer thicknesses, "// &
                                "surface first (exactly nz_layers entries). "// &
                                "tools/vgrid_to_dz.py converts a MOM6 ocean_vgrid.nc "// &
                                "(or a plain interface-depth list) into a ready-to-paste "// &
                                "block of this knob -- e.g. ACCESS-OM3's 75-level z* "// &
                                "grid, used with vcoord_type='zstar'", units="m"))
      pr => cfg%z_fixed_dz_top
      call g%add(nml_real("z_fixed_dz_top", pr, &
                          "z_fixed_profile='tanh': surface-layer nominal thickness", &
                          units="m"))
      pr => cfg%z_fixed_tanh_center
      call g%add(nml_real("z_fixed_tanh_center", pr, &
                          "z_fixed_profile='tanh': transition centre, fraction of the "// &
                          "layer-index span (0 = surface, 1 = bed)", min=0.0_wp, max=1.0_wp))
      pr => cfg%z_fixed_tanh_width
      call g%add(nml_real("z_fixed_tanh_width", pr, &
                          "z_fixed_profile='tanh': transition width, fraction of the "// &
                          "layer-index span"))
      call schema%add_group(g)
   end subroutine register_vcoord

   module subroutine register_physics(cfg, schema)
      !! `&physics_nml`: barotropic physics + bottom drag + wind stress.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      real(wp), pointer :: pr
      logical, pointer :: pl

      g%name = "physics"
      g%doc = "Barotropic physics: bottom drag, wind stress, Coriolis."
      pr => cfg%manning_n
      call g%add(nml_real("manning_n", pr, "Manning roughness coefficient", &
                          dead_on_ocean_path="coastal-legacy A-grid bottom-drag knob; the A-grid "// &
                          "god state that read it is gone (see rdb_state.F90's module docstring). "// &
                          "The ocean path's drag is &ocean_bdrag_nml."))
      pr => cfg%wind_stress_x
      call g%add(nml_real("wind_stress_x", pr, "Surface wind stress in x", units="Pa"))
      pr => cfg%wind_stress_y
      call g%add(nml_real("wind_stress_y", pr, "Surface wind stress in y", units="Pa"))
      pr => cfg%coriolis_f
      call g%add(nml_real("coriolis_f", pr, "Coriolis parameter f", units="1/s"))
      call schema%add_group(g)
   end subroutine register_physics

   module subroutine register_output(cfg, schema)
      !! `&output_nml`: file output + I/O + forcing/restart/gauge files.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      real(wp), pointer :: pr
      integer, pointer :: pi
      logical, pointer :: pl
      character(len=:), pointer :: ps

      g%name = "output"
      g%doc = "File output, I/O server, forcing/restart/gauge files."
      pl => cfg%output_to_file
      call g%add(nml_logical("output_to_file", pl, "Enable file output (NetCDF snapshots)", &
                             dead_on_ocean_path="accepted and validated but read nowhere in "// &
                             "src/ -- NetCDF output is gated by RDB_ENABLE_NETCDF at build "// &
                             "time and by the diag/restart registries at run time, not by this "// &
                             "flag (found by the P4 dead-knob sweep, 2026-09-10)."))
      ps => cfg%output_dir
      call g%add(nml_string("output_dir", ps, "Directory for output files"))
      pr => cfg%restart_interval
      call g%add(nml_real("restart_interval", pr, "Time between restart writes (0 = none)", units="s"))
      pl => cfg%compress_output
      call g%add(nml_logical("compress_output", pl, "Deflate-compress NetCDF output"))
      pi => cfg%compress_level
      call g%add(nml_int("compress_level", pi, "Deflate level (1=fast, 9=max)", min=1, max=9))
      pl => cfg%use_io_server
      call g%add(nml_logical("use_io_server", pl, "Dedicate one MPI rank per node as I/O server"))
      ps => cfg%bathymetry_file
      call g%add(nml_string("bathymetry_file", ps, "NetCDF bathymetry file (empty = flat)"))
      ps => cfg%restart_file
      call g%add(nml_string("restart_file", ps, "Restart file for warm start (empty = cold)"))
      call schema%add_group(g)
   end subroutine register_output

   module subroutine register_boundary(cfg, schema)
      !! `&boundary_nml`: per-side BC types, tidal forcing, inflow/
      !! discharge/clamped/sponge/nesting parameters.  The four `bc_*`
      !! keys are enums over the bc_type_from_string accepted set.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      real(wp), pointer :: pr
      real(wp), pointer :: pa(:)
      integer, pointer :: pi
      character(len=:), pointer :: ps
      character(len=9), parameter :: bc_allowed(*) = &
                                     [character(len=9) :: "wall", "open", "tidal", "nested", "inflow", &
                                                           "discharge", "clamped", "sponge", "chapman"]

      g%name = "boundary"
      g%doc = "Boundary conditions, tidal forcing, inflow/discharge/sponge/nesting."
      ps => cfg%bc_west
      call g%add(nml_enum("bc_west", ps, "West boundary type", allowed=bc_allowed, &
                          dead_on_ocean_path="coastal-legacy: reaches only warn_unknown_bc "// &
                          "and nothing else. The ocean path's edge selector is &ocean_bc_nml "// &
                          "west/east/south/north."))
      ps => cfg%bc_east
      call g%add(nml_enum("bc_east", ps, "East boundary type", allowed=bc_allowed, &
                          dead_on_ocean_path="coastal-legacy: reaches only warn_unknown_bc "// &
                          "and nothing else. The ocean path's edge selector is &ocean_bc_nml "// &
                          "west/east/south/north."))
      ps => cfg%bc_south
      call g%add(nml_enum("bc_south", ps, "South boundary type", allowed=bc_allowed, &
                          dead_on_ocean_path="coastal-legacy: reaches only warn_unknown_bc "// &
                          "and nothing else. The ocean path's edge selector is &ocean_bc_nml "// &
                          "west/east/south/north."))
      ps => cfg%bc_north
      call g%add(nml_enum("bc_north", ps, "North boundary type", allowed=bc_allowed, &
                          dead_on_ocean_path="coastal-legacy: reaches only warn_unknown_bc "// &
                          "and nothing else. The ocean path's edge selector is &ocean_bc_nml "// &
                          "west/east/south/north."))
      pi => cfg%n_tidal_constituents
      call g%add(nml_int("n_tidal_constituents", pi, &
                         "Active tidal constituents (0 = legacy single)", min=0))
      pa => cfg%tidal_amp
      call g%add(nml_real_array("tidal_amp", pa, "Constituent amplitudes", units="m"))
      pa => cfg%tidal_phase
      call g%add(nml_real_array("tidal_phase", pa, "Constituent phases", units="rad"))
      pa => cfg%tidal_omega
      call g%add(nml_real_array("tidal_omega", pa, "Constituent angular frequencies", units="rad/s"))
      pr => cfg%inflow_salinity
      call g%add(nml_real("inflow_salinity", pr, "Inflow salinity (<0 = zero-gradient)", units="PSU", &
                          dead_on_ocean_path="stored onto tracer_t%tr_inflow via "// &
                          "register_default_tracers but tr_inflow is read nowhere on the "// &
                          "ocean path (found by the P4 dead-knob sweep, 2026-09-10)."))
      pr => cfg%inflow_temperature
      call g%add(nml_real("inflow_temperature", pr, "Inflow temperature (<0 = zero-gradient)", units="degC", &
                          dead_on_ocean_path="stored onto tracer_t%tr_inflow via "// &
                          "register_default_tracers but tr_inflow is read nowhere on the "// &
                          "ocean path (found by the P4 dead-knob sweep, 2026-09-10)."))
      pi => cfg%sponge_width
      call g%add(nml_int("sponge_width", pi, "Sponge layer width in cells", min=0))
      pr => cfg%sponge_strength
      call g%add(nml_real("sponge_strength", pr, "Sponge relaxation rate", units="1/s"))
      call schema%add_group(g)
   end subroutine register_boundary

   module subroutine register_tracer(cfg, schema)
      !! `&tracer_nml`: salinity + temperature IC/EOS/bounds + sediment.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      real(wp), pointer :: pr
      logical, pointer :: pl
      character(len=:), pointer :: ps

      g%name = "tracer"
      g%doc = "Salinity + temperature IC/EOS/bounds, sediment transport."
      pr => cfg%initial_salinity
      call g%add(nml_real("initial_salinity", pr, "Initial salinity (uniform IC)", units="PSU"))
      pr => cfg%S_ref
      call g%add(nml_real("S_ref", pr, "EOS reference salinity", units="PSU", &
                          dead_on_ocean_path="RETIRED -- stored onto tracer_t%eos_ref via "// &
                          "register_default_tracers but eos_ref is read nowhere on the ocean "// &
                          "path. The live ocean-path spelling is &ocean_ic_nml S_ref; setting "// &
                          "THIS one to anything other than its default is a fail-loud "// &
                          "configure error (validate_config)."))
      pr => cfg%beta_S
      call g%add(nml_real("beta_S", pr, "Haline contraction coefficient", units="kg/m^3/PSU", &
                          dead_on_ocean_path="RETIRED -- stored onto tracer_t%eos_coeff via "// &
                          "register_default_tracers but eos_coeff is read nowhere on the "// &
                          "ocean path. The live ocean-path spelling is &ocean_ic_nml beta_S; "// &
                          "setting THIS one to anything other than its default is a fail-loud "// &
                          "configure error (validate_config)."))
      pr => cfg%S_min
      call g%add(nml_real("S_min", pr, "Lower physical bound for salinity", units="PSU", &
                          dead_on_ocean_path="stored onto tracer_t%tr_min via "// &
                          "register_default_tracers but tr_min is read nowhere on the ocean "// &
                          "path (found by the P4 dead-knob sweep, 2026-09-10)."))
      pr => cfg%S_max
      call g%add(nml_real("S_max", pr, "Upper physical bound for salinity", units="PSU", &
                          dead_on_ocean_path="stored onto tracer_t%tr_max via "// &
                          "register_default_tracers but tr_max is read nowhere on the ocean "// &
                          "path (found by the P4 dead-knob sweep, 2026-09-10)."))
      pr => cfg%kappa_S_bg
      call g%add(nml_real("kappa_S_bg", pr, "Background vertical salinity diffusivity", units="m^2/s", &
                          dead_on_ocean_path="stored onto tracer_t%kappa_bg via "// &
                          "register_default_tracers but kappa_bg is read nowhere on the "// &
                          "ocean path -- the ocean path's background salt diffusivity is "// &
                          "&ocean_vmix_nml ks_bg (found by the P4 dead-knob sweep, 2026-09-10)."))
      pr => cfg%S_init_surface
      call g%add(nml_real("S_init_surface", pr, &
                          "Initial surface salinity at k=nz (linear-in-layer stratified IC; "// &
                          "needs S_init_bottom non-zero too)", units="PSU", min=0.0_wp))
      pr => cfg%S_init_bottom
      call g%add(nml_real("S_init_bottom", pr, &
                          "Initial bed salinity at k=1 (linear-in-layer stratified IC; "// &
                          "needs S_init_surface non-zero too)", units="PSU", min=0.0_wp))
      pr => cfg%initial_temperature
      call g%add(nml_real("initial_temperature", pr, "Initial temperature (uniform IC)", units="degC"))
      pr => cfg%T_ref
      call g%add(nml_real("T_ref", pr, "EOS reference temperature", units="degC", &
                          dead_on_ocean_path="RETIRED -- stored onto tracer_t%eos_ref via "// &
                          "register_default_tracers but eos_ref is read nowhere on the ocean "// &
                          "path. The live ocean-path spelling is &ocean_ic_nml T_ref; setting "// &
                          "THIS one to anything other than its default is a fail-loud "// &
                          "configure error (validate_config)."))
      pr => cfg%alpha_T
      call g%add(nml_real("alpha_T", pr, "Thermal expansion coefficient", units="kg/m^3/degC", &
                          dead_on_ocean_path="RETIRED -- the coastal-legacy alpha_T (D2.5): "// &
                          "stored onto tracer_t%eos_coeff via register_default_tracers but "// &
                          "eos_coeff is read nowhere on the ocean path. The live ocean-path "// &
                          "spelling is &ocean_ic_nml alpha_T; setting THIS one to anything "// &
                          "other than its default is a fail-loud configure error "// &
                          "(validate_config)."))
      pr => cfg%T_min
      call g%add(nml_real("T_min", pr, "Lower physical bound for temperature", units="degC", &
                          dead_on_ocean_path="stored onto tracer_t%tr_min via "// &
                          "register_default_tracers but tr_min is read nowhere on the ocean "// &
                          "path (found by the P4 dead-knob sweep, 2026-09-10)."))
      pr => cfg%T_max
      call g%add(nml_real("T_max", pr, "Upper physical bound for temperature", units="degC", &
                          dead_on_ocean_path="stored onto tracer_t%tr_max via "// &
                          "register_default_tracers but tr_max is read nowhere on the ocean "// &
                          "path (found by the P4 dead-knob sweep, 2026-09-10)."))
      pr => cfg%kappa_T_bg
      call g%add(nml_real("kappa_T_bg", pr, "Background vertical temperature diffusivity", units="m^2/s", &
                          dead_on_ocean_path="stored onto tracer_t%kappa_bg via "// &
                          "register_default_tracers but kappa_bg is read nowhere on the "// &
                          "ocean path -- the ocean path's background heat diffusivity is "// &
                          "&ocean_vmix_nml kt_bg (found by the P4 dead-knob sweep, 2026-09-10)."))
      pr => cfg%T_init_surface
      call g%add(nml_real("T_init_surface", pr, "Initial surface temperature (stratified IC)", units="degC"))
      pr => cfg%T_init_bottom
      call g%add(nml_real("T_init_bottom", pr, "Initial bed temperature (stratified IC)", units="degC"))
      call schema%add_group(g)
   end subroutine register_tracer

   module subroutine register_initial_condition(cfg, schema)
      !! `&initial_condition_nml`: coastal IC selector + parameters.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      real(wp), pointer :: pr
      character(len=:), pointer :: ps

      g%name = "initial_condition"
      g%doc = "Coastal initial-condition selector + parameters."
      pr => cfg%h0
      call g%add(nml_real("h0", pr, "Background depth (gaussian_hump IC)", units="m"))
      call schema%add_group(g)
   end subroutine register_initial_condition

   module subroutine register_nonhydrostatic(cfg, schema)
      !! `&nonhydrostatic_nml`: NH/multilayer switches, CG-Poisson
      !! controls, PP81 vmix, KPP knobs, k-eps, Smagorinsky, BPG, mode
      !! split.  `keps_stability` and `bpg_method` are enums.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      real(wp), pointer :: pr
      integer, pointer :: pi
      logical, pointer :: pl
      character(len=:), pointer :: ps

      g%name = "nonhydrostatic"
      g%doc = "Non-hydrostatic / multilayer: CG-Poisson, PP81/KPP/k-eps vmix, Smagorinsky, BPG, mode split."
      pi => cfg%nz_layers
      call g%add(nml_int("nz_layers", pi, "Number of sigma layers", min=1))
      pl => cfg%use_multilayer
      call g%add(nml_logical("use_multilayer", pl, "Enable coupled hydrostatic vertical layers"))
      pr => cfg%rho_0
      call g%add(nml_real("rho_0", pr, "Reference density for EOS", units="kg/m^3"))
      pr => cfg%kpp_ri_crit
      call g%add(nml_real("kpp_ri_crit", pr, "Critical bulk Richardson number for KPP BL-depth"))
      pr => cfg%kpp_cs_nonlocal
      call g%add(nml_real("kpp_cs_nonlocal", pr, "KPP non-local (counter-gradient) transport coefficient"))
      pr => cfg%kpp_c_vt2
      call g%add(nml_real("kpp_c_vt2", pr, "KPP V_t^2 unresolved-turbulence coefficient (0 = off)"))
      pr => cfg%hdiff_kappa
      call g%add(nml_real("hdiff_kappa", pr, &
                          "Horizontal tracer diffusion coefficient (COASTAL path only; "// &
                          "ocean's along-coordinate equivalent is &ocean_hdiff_nml kappa_h)", &
                          units="m^2/s", &
                          dead_on_ocean_path="coastal-legacy, per this key's own doc string: "// &
                          "stored onto tracer_t%hdiff_kappa via register_default_tracers but "// &
                          "hdiff_kappa is read nowhere on the ocean path. Use "// &
                          "&ocean_hdiff_nml kappa_h."))
      call schema%add_group(g)
   end subroutine register_nonhydrostatic

   module subroutine register_ocean_grid(cfg, schema)
      !! `&ocean_grid_nml`: horizontal-grid generator + geometry
      !! (curvilinear-grid stream).  `grid_config` / `coriolis_scheme`
      !! enums mirror `parse_grid_config` / `parse_coriolis_scheme` in
      !! `rdb_ocean_metrics`.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      character(len=:), pointer :: ps
      real(wp), pointer :: pr

      g%name = "ocean_grid"
      g%doc = "Horizontal-grid generator + geometry (ocean curvilinear stream)."
      ps => cfg%ocean%grid%grid_config
      call g%add(nml_enum("grid_config", ps, &
                          "Horizontal grid generator", &
                          allowed=[character(len=12) :: "cartesian", "spherical", &
                                   "supergrid", "tripolar"]))
      pr => cfg%ocean%grid%lon_west
      call g%add(nml_real("lon_west", pr, &
                          "West edge of the domain (spherical)", units="degrees_east"))
      pr => cfg%ocean%grid%lat_south
      call g%add(nml_real("lat_south", pr, &
                          "South edge of the domain (spherical)", units="degrees_north", &
                          min=-90.0_wp, max=90.0_wp))
      pr => cfg%ocean%grid%rad_earth
      call g%add(nml_real("rad_earth", pr, "Earth radius for the spherical metric", &
                          units="m", min=0.0_wp))
      ps => cfg%ocean%grid%supergrid_file
      call g%add(nml_string("supergrid_file", ps, &
                            "Path to the MOM6 supergrid NetCDF (supergrid grid_config)"))
      ps => cfg%ocean%grid%coriolis_scheme
      call g%add(nml_enum("coriolis_scheme", ps, &
                          "Coriolis source for the metrics f-fill", &
                          allowed=[character(len=10) :: "beta_plane", "planetary"]))
      pr => cfg%ocean%grid%omega
      call g%add(nml_real("omega", pr, "Planetary rotation rate (planetary scheme)", &
                          units="rad/s"))
      pr => cfg%ocean%grid%phi_join
      call g%add(nml_real("phi_join", pr, &
                          "Join latitude for the tripolar bipolar cap", &
                          units="degrees_north", min=-90.0_wp, max=90.0_wp))
      pr => cfg%ocean%grid%lon_pole
      call g%add(nml_real("lon_pole", pr, &
                          "Longitude of the first tripolar cap pole (partner +180)", &
                          units="degrees_east"))
      ps => cfg%ocean%grid%axis_units
      call g%add(nml_enum("axis_units", ps, &
                          "Units of the Cartesian domain extent (MOM6 AXIS_UNITS)", &
                          allowed=[character(len=8) :: "meters", "degrees", "km"]))
      pr => cfg%ocean%grid%len_lon
      call g%add(nml_real("len_lon", pr, &
                          "Total x-extent of the Cartesian domain in axis_units "// &
                          "(MOM6 LENLON; derives dx when > 0)"))
      pr => cfg%ocean%grid%len_lat
      call g%add(nml_real("len_lat", pr, &
                          "Total y-extent of the Cartesian domain in axis_units "// &
                          "(MOM6 LENLAT; derives dy when > 0)"))
      call schema%add_group(g)
   end subroutine register_ocean_grid

   module subroutine register_ocean_debug(cfg, schema)
      !! `&ocean_debug_nml`: forensic probes (all heavy, all default off).
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      logical, pointer :: pl
      integer, pointer :: pi

      g%name = "ocean_debug"
      g%doc = "Forensic probes for the ocean dyn-core (all heavy, all default off)."
      pl => cfg%ocean%debug%budget
      call g%add(nml_logical("budget", pl, "Per-stage BT power-budget probe (BUDGET rows)"))
      pl => cfg%ocean%debug%ke_attr
      call g%add(nml_logical("ke_attr", pl, &
                             "Per-segment layer-KE attribution meter (KE_ATTR rows)"))
      pi => cfg%ocean%debug%ke_attr_start_step
      call g%add(nml_int("ke_attr_start_step", pi, &
                         "First outer step the KE meter samples (0 = from start)", min=0))
      pi => cfg%ocean%debug%ke_attr_end_step
      call g%add(nml_int("ke_attr_end_step", pi, &
                         "Last outer step the KE meter samples (0 = unbounded)", min=0))
      pl => cfg%ocean%debug%chksum
      call g%add(nml_logical("chksum", pl, &
                             "MOM6-style per-phase field checksums (CHKSUM rows) + "// &
                             "HOTFACE argmax-face anatomy at the tendency seams"))
      pi => cfg%ocean%debug%chksum_start_step
      call g%add(nml_int("chksum_start_step", pi, &
                         "First outer step the chksum probe samples (0 = from start)", min=0))
      pi => cfg%ocean%debug%chksum_end_step
      call g%add(nml_int("chksum_end_step", pi, &
                         "Last outer step the chksum probe samples (0 = unbounded)", min=0))
      pl => cfg%ocean%debug%chksum_interior
      call g%add(nml_logical("chksum_interior", pl, &
                             "Reduce over physical cells only, making `bits` a valid "// &
                             "1-rank-vs-N-rank gate"))
      call schema%add_group(g)
   end subroutine register_ocean_debug

   module subroutine register_ocean_mpi(cfg, schema)
      !! `&ocean_mpi_nml`: multi-rank MPI debug / tuning controls.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      logical, pointer :: pl

      g%name = "ocean_mpi"
      g%doc = "Multi-rank MPI debug/tuning controls."
      pl => cfg%ocean%mpi%poison_ghosts
      call g%add(nml_logical("poison_ghosts", pl, &
                             "Sentinel-NaN the exchange-covered ghost bands at each outer-step "// &
                             "start; unexchanged-ghost consumption becomes a loud NaN. "// &
                             "Default off = bit-identical."))
      call schema%add_group(g)
   end subroutine register_ocean_mpi

   module subroutine register_ocean_topo(cfg, schema)
      !! `&ocean_topo_nml`: basin geometry + surface forcing + Coriolis tilt.
      !! `topo_config` / `wind_config` enums mirror the dispatch
      !! select-cases in rdb_ocean_state / rdb_ocean_setup.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      real(wp), pointer :: pr
      character(len=:), pointer :: ps

      g%name = "ocean_topo"
      g%doc = "Basin geometry + surface forcing + Coriolis tilt."
      ps => cfg%ocean%topo%topo_config
      call g%add(nml_enum("topo_config", ps, "Bathymetry profile selector", &
                          allowed=[character(len=12) :: "flat", "spoon", "seamount", &
                                   "neverworld2", "island", "double_drake", &
                                   "isomip_plus", "file"]))
      pr => cfg%ocean%topo%max_depth
      call g%add(nml_real("max_depth", pr, "Basin maximum depth", units="m"))
      pr => cfg%ocean%topo%edge_depth
      call g%add(nml_real("edge_depth", pr, "Spoon-bathymetry edge depth", units="m"))
      pr => cfg%ocean%topo%slope_scale
      call g%add(nml_real("slope_scale", pr, "Spoon-bathymetry exponential decay scale", units="m"))
      pr => cfg%ocean%topo%nl_continent_amp
      call g%add(nml_real("nl_continent_amp", pr, &
                          "Neverworld2 continent amplitude (1=full continents, 0=aquaplanet+channel)", &
                          units="nondim", min=0.0_wp))
      pr => cfg%ocean%topo%nl_roughness_amp
      call g%add(nml_real("nl_roughness_amp", pr, &
                          "Neverworld2 bathymetry roughness amplitude", units="nondim", min=0.0_wp))
      pr => cfg%ocean%topo%nl_min_depth
      call g%add(nml_real("nl_min_depth", pr, &
                          "Neverworld2 minimum-depth floor (MOM6 MINIMUM_DEPTH analogue)", &
                          units="m", min=0.0_wp))
      ps => cfg%ocean%topo%wind_config
      call g%add(nml_enum("wind_config", ps, "Surface wind-stress dispatch", &
                          allowed=[character(len=11) :: "constant", "2gyre", "neverworld2"]))
      pr => cfg%ocean%topo%taux_magnitude
      call g%add(nml_real("taux_magnitude", pr, "Peak zonal wind stress for wind_config=2gyre/neverworld2", &
                          units="Pa"))
      pr => cfg%ocean%topo%coriolis_beta
      call g%add(nml_real("coriolis_beta", pr, "Meridional gradient of f (0 = f-plane)", &
                          units="1/(s m)"))
      pr => cfg%ocean%topo%coriolis_y_ref
      call g%add(nml_real("coriolis_y_ref", pr, "Reference y where f = coriolis_f under beta-plane", &
                          units="m"))
      pr => cfg%ocean%topo%x_origin
      call g%add(nml_real("x_origin", pr, &
                          "Absolute x of the domain west edge, for topo_config='isomip_plus' "// &
                          "(ISOMIP+ ocean box starts at the MISMIP+ x = 320 km)", units="m"))
      call schema%add_group(g)
   end subroutine register_ocean_topo

   module subroutine register_ocean_ic(cfg, schema)
      !! `&ocean_ic_nml`: initial-condition overlay + EOS reference state.
      !! `ic_config` enum mirrors the IC dispatch in rdb_ocean_state
      !! (default "" keeps the analytical T(z) IC).
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      real(wp), pointer :: pr
      real(wp), pointer :: pra(:)
      integer, pointer :: pi
      character(len=:), pointer :: ps

      g%name = "ocean_ic"
      g%doc = "Initial-condition overlay + EOS reference state."
      ps => cfg%ocean%ic%ic_config
      call g%add(nml_enum("ic_config", ps, "IC overlay tag", &
                          allowed=[character(len=22) :: "", "eady", &
                                   "geostrophic_adjustment", "baroclinic_jet"]))
      pr => cfg%ocean%ic%alpha_T
      call g%add(nml_real("alpha_T", pr, &
                          "Linear-EOS thermal-expansion coefficient (DIMENSIONAL: "// &
                          "multiply a fractional 1/degC coefficient by rho_0)", &
                          units="kg/m^3/degC"))
      pr => cfg%ocean%ic%beta_S
      call g%add(nml_real("beta_S", pr, &
                          "Linear-EOS haline contraction coefficient (DIMENSIONAL: "// &
                          "multiply a fractional 1/PSU coefficient by rho_0)", &
                          units="kg/m^3/PSU", min=0.0_wp))
      pr => cfg%ocean%ic%T_ref
      call g%add(nml_real("T_ref", pr, "Linear-EOS reference temperature", &
                          units="degC", min=-273.15_wp))
      pr => cfg%ocean%ic%S_ref
      call g%add(nml_real("S_ref", pr, "Linear-EOS reference salinity", &
                          units="PSU", min=0.0_wp))
      pr => cfg%ocean%ic%rho_0
      call g%add(nml_real("rho_0", pr, "Reference density for the linear EOS / Boussinesq PGF", &
                          units="kg/m^3"))
      pra => cfg%ocean%ic%layer_rho_init
      call g%add(nml_real_array("layer_rho_init", pra, &
                                "Per-layer initial density (k=1 bed -> k=nz surface; -1 = EOS init)", &
                                units="kg/m^3"))
      pr => cfg%ocean%ic%rho_lightest
      call g%add(nml_real("rho_lightest", pr, &
                          "Linear density-range IC: surface (lightest) layer density; -1 = off", &
                          units="kg/m^3"))
      pr => cfg%ocean%ic%rho_range
      call g%add(nml_real("rho_range", pr, &
                          "Linear density-range IC: top-to-bottom density contrast (MOM6 DENSITY_RANGE)", &
                          units="kg/m^3"))
      pr => cfg%ocean%ic%eady_dT_dy
      call g%add(nml_real("eady_dT_dy", pr, "Eady IC meridional T gradient", units="degC/m"))
      pr => cfg%ocean%ic%eady_dT_dz
      call g%add(nml_real("eady_dT_dz", pr, "Eady IC vertical stratification", units="degC/m"))
      pr => cfg%ocean%ic%eady_T_ref
      call g%add(nml_real("eady_T_ref", pr, "Eady IC reference temperature at z=0,y=y_mid", &
                          units="degC"))
      pr => cfg%ocean%ic%eady_pert_amp
      call g%add(nml_real("eady_pert_amp", pr, "Eady IC symmetry-breaking perturbation amplitude", &
                          units="degC"))
      pi => cfg%ocean%ic%eady_pert_seed
      call g%add(nml_int("eady_pert_seed", pi, "RNG seed for the Eady IC perturbation"))
      pr => cfg%ocean%ic%ga_eta_amp
      call g%add(nml_real("ga_eta_amp", pr, "Geostrophic-adjustment IC: SSH-bump amplitude", units="m"))
      pr => cfg%ocean%ic%ga_length_scale
      call g%add(nml_real("ga_length_scale", pr, "Geostrophic-adjustment IC: SSH-bump e-folding scale", &
                          units="m"))
      pr => cfg%ocean%ic%ga_x_center
      call g%add(nml_real("ga_x_center", pr, "Geostrophic-adjustment IC: bump x-centre (<0 = auto)", &
                          units="m"))
      pr => cfg%ocean%ic%ga_y_center
      call g%add(nml_real("ga_y_center", pr, "Geostrophic-adjustment IC: bump y-centre (<0 = auto)", &
                          units="m"))
      pr => cfg%ocean%ic%jet_half_width
      call g%add(nml_real("jet_half_width", pr, &
                          "Baroclinic-jet IC: tanh jet half-width L", units="m"))
      pr => cfg%ocean%ic%interface_amp
      call g%add(nml_real("interface_amp", pr, &
                          "Baroclinic-jet IC: interface displacement amplitude", units="m"))
      pr => cfg%ocean%ic%pert_amp_frac
      call g%add(nml_real("pert_amp_frac", pr, &
                          "Baroclinic-jet IC: meander amplitude as a fraction of interface_amp"))
      pi => cfg%ocean%ic%pert_nx
      call g%add(nml_int("pert_nx", pi, &
                         "Baroclinic-jet IC: zonal perturbation wavenumber (integer)"))
      pr => cfg%ocean%ic%upper_layer_rest
      call g%add(nml_real("upper_layer_rest", pr, &
                          "Baroclinic-jet IC: upper (surface) layer rest thickness H1", &
                          units="m"))
      call schema%add_group(g)
   end subroutine register_ocean_ic

   module subroutine register_ocean_zinit(cfg, schema)
      !! `&ocean_zinit_nml`: z-level T/S initial-condition overlay (A2).
      !! Default `enable = .false.` is a no-op (analytical IC unchanged).
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      logical, pointer :: pl
      real(wp), pointer :: pr
      character(len=:), pointer :: ps

      g%name = "ocean_zinit"
      g%doc = "Z-level T/S initial-condition overlay (A2)."
      pl => cfg%ocean%zinit%enable
      call g%add(nml_logical("enable", pl, &
                             "Master switch (default off; requires RDB_ENABLE_NETCDF=ON)"))
      ps => cfg%ocean%zinit%source
      call g%add(nml_enum("source", ps, &
                          "Where the T(z)/S(z) profile comes from: a "// &
                          "pre-regridded NetCDF, or the analytic affine "// &
                          "lin_* profile (no file)", &
                          allowed=[character(len=15) :: "file", "linear"]))
      ps => cfg%ocean%zinit%file
      call g%add(nml_string("file", ps, "Path to the model-grid T/S NetCDF"))
      ps => cfg%ocean%zinit%t_var
      call g%add(nml_string("t_var", ps, &
                            "Temperature variable-name override (blank tries temp/T/temperature)"))
      ps => cfg%ocean%zinit%s_var
      call g%add(nml_string("s_var", ps, &
                            "Salinity variable-name override (blank tries salt/S/salinity)"))
      ps => cfg%ocean%zinit%z_var
      call g%add(nml_string("z_var", ps, &
                            "Source-axis variable-name override (blank tries z_src/z/depth/lev)"))
      pr => cfg%ocean%zinit%land_fill_t
      call g%add(nml_real("land_fill_t", pr, "Fallback temperature for dry columns", units="degC"))
      pr => cfg%ocean%zinit%land_fill_s
      call g%add(nml_real("land_fill_s", pr, "Fallback salinity for dry columns", units="PSU"))
      pr => cfg%ocean%zinit%lin_t_ref
      call g%add(nml_real("lin_t_ref", pr, &
                          "source='linear': temperature at the z = 0 datum", units="degC"))
      pr => cfg%ocean%zinit%lin_dt_dz
      call g%add(nml_real("lin_dt_dz", pr, &
                          "source='linear': dT/dz, z positive UP (stable > 0)", &
                          units="degC/m"))
      pr => cfg%ocean%zinit%lin_s_ref
      call g%add(nml_real("lin_s_ref", pr, &
                          "source='linear': salinity at the z = 0 datum", units="PSU"))
      pr => cfg%ocean%zinit%lin_ds_dz
      call g%add(nml_real("lin_ds_dz", pr, &
                          "source='linear': dS/dz, z positive UP (stable < 0)", &
                          units="PSU/m"))
      call schema%add_group(g)
   end subroutine register_ocean_zinit

   module subroutine register_ocean_data(cfg, schema)
      !! `&ocean_data_nml`: the shared time-varying NetCDF input reader
      !! (PR-14).  Two knobs, no per-field entries — see
      !! `ocean_data_config_t`.  Registration of individual (file,
      !! variable, destination) triples is programmatic, via each
      !! consumer's own namelist group calling
      !! `ocean_data_input_register_2d/_3d` at setup.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      integer, pointer :: pi
      logical, pointer :: pl

      g%name = "ocean_data"
      g%doc = "Shared time-varying NetCDF input reader (PR-14): registry sizing only."
      pi => cfg%ocean%data%max_fields
      call g%add(nml_int("max_fields", pi, &
                         "Size of the reader's field registry", min=1))
      pl => cfg%ocean%data%verbose
      call g%add(nml_logical("verbose", pl, &
                             "Log every bracket advance (field, records, weight, time)"))
      call schema%add_group(g)
   end subroutine register_ocean_data

   module subroutine register_ocean_dataovr(cfg, schema)
      !! `&ocean_dataovr_nml`: file-backed surface forcing (PR-15).  Flat
      !! per-tag knobs (`<tag>_file`, `<tag>_var`, `<tag>_scale`,
      !! `<tag>_add`) rather than a MOM6-style parallel-array data table
      !! — the schema engine has no string-array key type, and the set of
      !! recognised tags is fixed by the slots they feed, so a flat
      !! layout is both expressible today and self-documenting in
      !! `docs/generated_nml_knobs.md`.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      logical, pointer :: pl
      real(wp), pointer :: pr
      character(len=:), pointer :: ps

      g%name = "ocean_dataovr"
      g%doc = "File-backed surface forcing (PR-15): wind stress, heat, freshwater."
      pl => cfg%ocean%dataovr%enable
      call g%add(nml_logical("enable", pl, &
                             "Master switch (default off; requires RDB_ENABLE_NETCDF=ON)"))
      ps => cfg%ocean%dataovr%time_mode
      call g%add(nml_enum("time_mode", ps, &
                          "Shared time-axis mode for every active tag", &
                          allowed=[character(len=6) :: "linear", "cyclic", "static"]))
      pr => cfg%ocean%dataovr%cycle_period
      call g%add(nml_real("cycle_period", pr, &
                          "Climatology period (required > 0 when time_mode=cyclic)", &
                          units="s", min=0.0_wp))
      pr => cfg%ocean%dataovr%t_offset
      call g%add(nml_real("t_offset", pr, &
                          "Added to the model time before the file-axis lookup", units="s"))
      pl => cfg%ocean%dataovr%oor_clamp
      call g%add(nml_logical("oor_clamp", pl, &
                             "Clamp an out-of-range query to the end record instead of aborting"))

      call add_dataovr_entry(g, cfg%ocean%dataovr%tau_x, "tau_x", &
                             "zonal wind stress", "Pa")
      call add_dataovr_entry(g, cfg%ocean%dataovr%tau_y, "tau_y", &
                             "meridional wind stress", "Pa")
      call add_dataovr_entry(g, cfg%ocean%dataovr%heat, "heat", &
                             "net surface heat flux (positive down)", "W m-2")
      call add_dataovr_entry(g, cfg%ocean%dataovr%evap, "evap", &
                             "evaporative mass flux (<= 0)", "kg m-2 s-1")
      call add_dataovr_entry(g, cfg%ocean%dataovr%lprec, "lprec", &
                             "liquid precipitation (>= 0)", "kg m-2 s-1")
      call add_dataovr_entry(g, cfg%ocean%dataovr%salt, "salt", &
                             "surface salt flux (positive salinifies)", "kg m-2 s-1")

      call schema%add_group(g)
   end subroutine register_ocean_dataovr

   module subroutine add_dataovr_entry(g, e, tag, what, units)
      !! Register the four flat knobs of one `&ocean_dataovr_nml` tag.
      !! `e` is `target, intent(in)` for the same reason `cfg` is in every
      !! `register_ocean_*` above: the schema stores a pointer to the
      !! target and assigns through it at parse time, and the actual
      !! argument is always a component of the `target` `cfg`, so the
      !! association outlives this call.
      type(nml_group_t), intent(inout) :: g
      type(dataovr_entry_config_t), target, intent(in) :: e
      character(len=*), intent(in) :: tag, what, units
      real(wp), pointer :: pr
      character(len=:), pointer :: ps

      ps => e%file
      call g%add(nml_string(tag//"_file", ps, &
                            "Path to the "//what//" NetCDF (blank ⇒ tag not file-driven)"))
      ps => e%var
      call g%add(nml_string(tag//"_var", ps, &
                            "Variable name inside "//tag//"_file (required when it is set)"))
      pr => e%scale
      call g%add(nml_real(tag//"_scale", pr, &
                          "Multiplier applied to the "//what//" slab at read", units=units))
      pr => e%add_offset
      call g%add(nml_real(tag//"_add", pr, &
                          "Offset added to the "//what//" slab after scale", units=units))
   end subroutine add_dataovr_entry

   module subroutine register_ocean_diag(cfg, schema)
      !! `&ocean_diag_nml`: ocean diag-manager output controls.  `dt_out`
      !! is interpreted in `time_unit` and converted by the post-parse
      !! cascade in read_config.  `vgrid` enum is {layer, z_fixed}.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      logical, pointer :: pl
      real(wp), pointer :: pr
      real(wp), pointer :: pra(:)
      integer, pointer :: pi
      character(len=:), pointer :: ps

      g%name = "ocean_diag"
      g%doc = "Ocean diag-manager output controls."
      pl => cfg%ocean%diag%enabled
      call g%add(nml_logical("enabled", pl, "Enable the per-step diag-manager hook"))
      ps => cfg%ocean%diag%filename
      call g%add(nml_string("filename", ps, "Output file basename (per-rank suffix appended)"))
      pr => cfg%ocean%diag%dt_out
      call g%add(nml_real("dt_out", pr, "Diag-fire cadence (in time_unit)"))
      ps => cfg%ocean%diag%vgrid
      call g%add(nml_enum("vgrid", ps, "Default output vertical grid for layer-shaped diags", &
                          allowed=[character(len=7) :: "layer", "z_fixed", "sigma", "zstar", "density"]))
      pra => cfg%ocean%diag%z_levels
      call g%add(nml_real_array("z_levels", pra, "Output z-levels when vgrid=z_fixed (positive down)", &
                                units="m"))
      pi => cfg%ocean%diag%n_z_levels
      call g%add(nml_int("n_z_levels", pi, "Number of z_levels entries used (0 = none)", min=0))
      pra => cfg%ocean%diag%sigma_levels
      call g%add(nml_real_array("sigma_levels", pra, &
                                "Output sigma fractions (0..1) for sigma output; empty = auto-uniform"))
      pi => cfg%ocean%diag%n_sigma_levels
      call g%add(nml_int("n_sigma_levels", pi, "Number of sigma_levels entries used (0 = auto)", min=0))
      pra => cfg%ocean%diag%zstar_levels
      call g%add(nml_real_array("zstar_levels", pra, &
                                "Output z* reference depths for zstar output; empty = auto-uniform", &
                                units="m"))
      pi => cfg%ocean%diag%n_zstar_levels
      call g%add(nml_int("n_zstar_levels", pi, "Number of zstar_levels entries used (0 = auto)", min=0))
      ! NVHPC 26.5 codegen workaround (verified by bisection, see PR-9 report):
      ! passing a 4th `nml_real_array(...)` function-result temporary
      ! directly to `g%add(...)` inside THIS subroutine corrupts unrelated
      ! heap state (manifests as a SIGSEGV deep in an unrelated later
      ! namelist read, e.g. the old hand-rolled `read_ocean_bc_nml`, since
      ! removed by the P4.5 schema migration — `register_ocean_bc` uses
      ! the same workaround for its twelve real_array knobs).  z_levels/sigma_levels/
      ! zstar_levels (the first 3 real_array calls here) are unaffected;
      ! materialising the result into an explicit local variable before
      ! `g%add` avoids the miscompile.  Do not "simplify" this back to the
      ! one-line call-expression form.
      block
         type(nml_real_array_key_t) :: rho_key
         pra => cfg%ocean%diag%rho_levels
         rho_key = nml_real_array("rho_levels", pra, &
                                  "Output potential-density bin edges when vgrid=density; "// &
                                  "strictly increasing, light->dense", &
                                  units="kg/m^3")
         call g%add(rho_key)
      end block
      pi => cfg%ocean%diag%n_rho_levels
      call g%add(nml_int("n_rho_levels", pi, "Number of rho_levels entries used (0 = none)", min=0))
      pl => cfg%ocean%diag%mask_vanished_layers
      call g%add(nml_logical("mask_vanished_layers", pl, &
                             "Mask below-bottom/pinched remap cells to missing_value (default off)"))
      pl => cfg%ocean%diag%reproducing_sums
      call g%add(nml_logical("reproducing_sums", pl, &
                             "Order-invariant EFP console totals + budget terms: identical "// &
                             "console on every rank count, exact drift residual "// &
                             "(.false. = the pre-v0.1.0 FP sums)"))
      ps => cfg%ocean%diag%output_precision
      call g%add(nml_enum("output_precision", ps, &
                          "Element width of the diag NetCDF data vars; "// &
                          "'single' halves the bytes written (restarts/gauges/"// &
                          "console totals stay double regardless)", &
                          allowed=[character(len=6) :: "double", "single"]))
      ps => cfg%ocean%diag%diag_remap_scheme
      call g%add(nml_enum("diag_remap_scheme", ps, &
                          "Reconstruction for the conservative diagnostic vertical remap", &
                          allowed=[character(len=6) :: "pcm", "plm", "ppm", "ppm_h4", "pqm"]))
      ps => cfg%ocean%diag%diags
      call g%add(nml_string("diags", ps, &
                            "Unified diag list modifying the default set: name[:off][:cadence][:op] "// &
                            "(e.g. 'vorticity_z:1d  KE:off  temperature:6h:mean'); "// &
                            "default empty = canonical set"))
      call schema%add_group(g)
   end subroutine register_ocean_diag

end submodule rdb_config_register_core
