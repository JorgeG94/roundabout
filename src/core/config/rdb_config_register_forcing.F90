submodule(rdb_config) rdb_config_register_forcing
   !! Strict-schema registration: forcing/boundary/ice/cavity groups.
   implicit none

contains

   module subroutine register_ocean_tides(cfg, schema)
      !! `&ocean_tides` (C1 equilibrium astronomical body-force tide).
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      logical, pointer :: pl
      real(wp), pointer :: pr
      character(len=:), pointer :: ps

      g%name = "ocean_tides"
      g%doc = "Equilibrium (astronomical) body-force tide (C1) + scalar SAL (C2)."

      pl => cfg%ocean%tides%enable
      call g%add(nml_logical("enable", pl, &
                             "Master switch (requires a non-cartesian grid)"))
      pl => cfg%ocean%tides%use_sal
      call g%add(nml_logical("use_sal", pl, &
                             "Apply scalar self-attraction & loading (C2)"))
      pr => cfg%ocean%tides%beta_sal
      call g%add(nml_real("beta_sal", pr, &
                          "Scalar SAL factor beta (~0.085-0.12)"))
      pl => cfg%ocean%tides%add_nodal
      call g%add(nml_logical("add_nodal", pl, &
                             "Apply the 18.6-yr nodal f/u corrections"))
      ps => cfg%ocean%tides%constituents
      call g%add(nml_string("constituents", ps, &
                            "Active constituent list (e.g. 'M2 S2 K1 O1')"))
      ps => cfg%ocean%tides%ref_date
      call g%add(nml_string("ref_date", ps, &
                            "Astronomical reference date YYYY-MM-DD"))
      ps => cfg%ocean%tides%nodal_ref_date
      call g%add(nml_string("nodal_ref_date", ps, &
                            "Nodal reference date ('' => ref_date)"))

      call schema%add_group(g)
   end subroutine register_ocean_tides

   module subroutine register_ocean_psurf(cfg, schema)
      !! `&ocean_psurf` (PR-17 atmospheric surface-pressure loading /
      !! inverse barometer).
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      logical, pointer :: pl
      real(wp), pointer :: pr

      g%name = "ocean_psurf"
      g%doc = "Atmospheric surface-pressure loading / inverse barometer (PR-17)."

      pl => cfg%ocean%psurf%enable
      call g%add(nml_logical("enable", pl, &
                             "Master switch (split-solver only; requires "// &
                             "&ocean_forcing_nml enable_components=.true.)"))
      pl => cfg%ocean%psurf%in_eos
      call g%add(nml_logical("in_eos", pl, &
                             "Also feed the surface load to the equation of "// &
                             "state as the top-of-column pressure p_top "// &
                             "(requires enable=.true.; refused with the "// &
                             "unported pressure builders)"))
      pr => cfg%ocean%psurf%p_surf_const
      call g%add(nml_real("p_surf_const", pr, &
                          "Uniform atmospheric surface pressure (Pa) seeded "// &
                          "into p_surf_atm (uniform => provably inert)"))

      call schema%add_group(g)
   end subroutine register_ocean_psurf

   module subroutine register_ocean_cavity_dyn(cfg, schema)
      !! `&ocean_cavity_dyn` (P5.1 static ice-shelf cavity geometry: the
      !! prescribed draft + the barotropic datum that absorbs it).
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      logical, pointer :: pl
      real(wp), pointer :: pr
      character(len=:), pointer :: ps

      g%name = "ocean_cavity_dyn"
      g%doc = "Static ice-shelf cavity geometry: prescribed draft + "// &
              "barotropic datum bt_H_ref = b - z_draft."

      pl => cfg%ocean%cavity_dyn%enable
      call g%add(nml_logical("enable", pl, &
                             "Master switch (single-rank, split solver, "// &
                             "fv_mom6 PGF, sigma/z_fixed only)"))
      ps => cfg%ocean%cavity_dyn%draft_config
      call g%add(nml_enum("draft_config", ps, &
                          "Draft source: analytic shape, or 'file' (static 2-D "// &
                          "NetCDF on the model grid, single rank)", &
                          allowed=[character(len=15) :: "none", "flat", &
                                   "linear", "file"]))
      ps => cfg%ocean%cavity_dyn%draft_source
      call g%add(nml_enum("draft_source", ps, &
                          "Whether the formula gives the ice-base DEPTH or "// &
                          "an ice THICKNESS ('in_situ' isostasy is deferred)", &
                          allowed=[character(len=15) :: "draft", "thickness", &
                                   "in_situ"]))
      pr => cfg%ocean%cavity_dyn%draft_depth
      call g%add(nml_real("draft_depth", pr, &
                          "Draft amplitude (ice thickness under "// &
                          "draft_source='thickness')", units="m"))
      pr => cfg%ocean%cavity_dyn%draft_slope
      call g%add(nml_real("draft_slope", pr, &
                          "d(draft)/dx for draft_config='linear' "// &
                          "(dimensionless; converted to grid units)"))
      pr => cfg%ocean%cavity_dyn%draft_x0
      call g%add(nml_real("draft_x0", pr, &
                          "Western edge of the shelf box, and the anchor of "// &
                          "the 'linear' profile (+/-1e30 => no limit)", units="m"))
      pr => cfg%ocean%cavity_dyn%draft_x1
      call g%add(nml_real("draft_x1", pr, &
                          "Eastern edge of the shelf box = the calving front "// &
                          "(+/-1e30 => no limit)", units="m"))
      pr => cfg%ocean%cavity_dyn%draft_y0
      call g%add(nml_real("draft_y0", pr, &
                          "Southern edge of the shelf box (+/-1e30 => no limit)", &
                          units="m"))
      pr => cfg%ocean%cavity_dyn%draft_y1
      call g%add(nml_real("draft_y1", pr, &
                          "Northern edge of the shelf box (+/-1e30 => no limit)", &
                          units="m"))
      ps => cfg%ocean%cavity_dyn%draft_file
      call g%add(nml_string("draft_file", ps, &
                            "draft_config='file': NetCDF path (variable must be "// &
                            "(x,y,t) Fortran order, on the model grid; record 1 read)"))
      ps => cfg%ocean%cavity_dyn%draft_var
      call g%add(nml_string("draft_var", ps, &
                            "draft_config='file': 2-D variable name (ISOMIP+ ships "// &
                            "'iceDraft')"))
      ps => cfg%ocean%cavity_dyn%draft_sign
      call g%add(nml_enum("draft_sign", ps, &
                          "draft_config='file': sign convention of the file values "// &
                          "(ISOMIP+ iceDraft is an ELEVATION)", &
                          allowed=[character(len=14) :: "depth", "positive_down", &
                                   "elevation", "positive_up"]))
      pr => cfg%ocean%cavity_dyn%h_min_cavity
      call g%add(nml_real("h_min_cavity", pr, &
                          "Grounding cutoff: b - z_draft below this is LAND "// &
                          "(never a thin film under grounded ice)", units="m"))
      pr => cfg%ocean%cavity_dyn%grounded_max_frac
      call g%add(nml_real("grounded_max_frac", pr, &
                          "Fail loud if more than this fraction of the "// &
                          "interior columns ground"))
      pr => cfg%ocean%cavity_dyn%rho_ice
      call g%add(nml_real("rho_ice", pr, &
                          "Ice density, consulted only by "// &
                          "draft_source='thickness'", units="kg/m^3"))
      pl => cfg%ocean%cavity_dyn%trim_ic_for_p_surf
      call g%add(nml_logical("trim_ic_for_p_surf", pl, &
                             "Trim the initial column top so the displaced "// &
                             "water's weight equals the ice load (MOM6 "// &
                             "TRIM_IC_FOR_P_SURF; linear EOS + zinit "// &
                             "source='linear' only)"))

      call schema%add_group(g)
   end subroutine register_ocean_cavity_dyn

   module subroutine register_ocean_cavity_melt(cfg, schema)
      !! `&ocean_cavity_melt` (P2b ice-shelf basal-melt thermodynamics:
      !! the three-equation interface, its exchange law and the
      !! far-field sampling depth).  The `exchange_law` and
      !! `ice_conduction` enums MIRROR `parse_cavity_exchange_law` /
      !! `parse_cavity_ice_mode` in `rdb_ocean_cavity_melt` — the two
      !! lists move together, and `validate_config` re-checks them
      !! belt-and-braces so a RESERVED law is refused by name rather
      !! than silently falling through the kernel's dispatch.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      logical, pointer :: pl
      real(wp), pointer :: pr
      character(len=:), pointer :: ps

      g%name = "ocean_cavity_melt"
      g%doc = "Ice-shelf basal-melt thermodynamics: the three-equation "// &
              "interface, its exchange law, and the far-field sampling depth."

      pl => cfg%ocean%cavity_melt%enable
      call g%add(nml_logical("enable", pl, &
                             "Master switch (requires &ocean_cavity_dyn_nml, "// &
                             "tfreeze_set='isomip' and the surface-flux "// &
                             "component set)"))
      ps => cfg%ocean%cavity_melt%exchange_law
      call g%add(nml_enum("exchange_law", ps, &
                          "Turbulent exchange law; laws other than "// &
                          "const_gamma/hj99/yung25 are RESERVED and refused "// &
                          "at configure", &
                          allowed=[character(len=12) :: "const_gamma", "hj99", &
                                   "yung25", "jenkins91", "rosevear22", "vt19", &
                                   "mk18", "burchard22", "jenkins21"]))
      pr => cfg%ocean%cavity_melt%gamma_t
      call g%add(nml_real("gamma_t", pr, &
                          "Dimensionless heat-transfer coefficient Gamma_T "// &
                          "(ISOMIP+ starting guess; tune per coordinate)", &
                          min=0.0_wp))
      pr => cfg%ocean%cavity_melt%gamma_s
      call g%add(nml_real("gamma_s", pr, &
                          "Dimensionless salt-transfer coefficient Gamma_S "// &
                          "(negative = unset = gamma_t/35)"))
      pr => cfg%ocean%cavity_melt%cdrag_top
      call g%add(nml_real("cdrag_top", pr, &
                          "Top drag coefficient for the MELT friction velocity "// &
                          "(no momentum drag yet - that is Phase 4)", &
                          min=0.0_wp))
      pr => cfg%ocean%cavity_melt%u_tide
      call g%add(nml_real("u_tide", pr, &
                          "RMS tidal velocity in the melt u* only, never the drag", &
                          units="m/s", min=0.0_wp))
      pr => cfg%ocean%cavity_melt%ustar_min
      call g%add(nml_real("ustar_min", pr, &
                          "Friction-velocity floor (Yung et al. 2025 eq. 14)", &
                          units="m/s", min=0.0_wp))
      ps => cfg%ocean%cavity_melt%ice_conduction
      call g%add(nml_enum("ice_conduction", ps, &
                          "Ice-side conduction; 'diffusive' is RESERVED and "// &
                          "refused (it changes the melt/freeze branch logic)", &
                          allowed=[character(len=10) :: "insulating", "adv_diff", &
                                   "diffusive"]))
      pr => cfg%ocean%cavity_melt%t_ice
      call g%add(nml_real("t_ice", pr, &
                          "Ice interior temperature; read by ice_conduction="// &
                          "'adv_diff' only", units="degC"))
      pr => cfg%ocean%cavity_melt%s_ice
      call g%add(nml_real("s_ice", pr, &
                          "Ice salinity; must stay strictly below the far-field "// &
                          "salinity", units="g/kg", min=0.0_wp))
      pr => cfg%ocean%cavity_melt%far_field_depth
      call g%add(nml_real("far_field_depth", pr, &
                          "Thickness below the ice base the far-field T/S/u are "// &
                          "averaged over (METRES, not layers)", units="m", &
                          min=0.0_wp))
      ps => cfg%ocean%cavity_melt%freshwater
      call g%add(nml_enum("freshwater", ps, &
                          "Meltwater delivery: 'virtual' (default, fixed "// &
                          "column mass) or 'mass' (real Boussinesq volume on "// &
                          "the top layer)", &
                          allowed=[character(len=8) :: "virtual", "mass"]))
      ps => cfg%ocean%cavity_melt%volume_compensation
      call g%add(nml_enum("volume_compensation", ps, &
                          "Sea-level compensation for freshwater='mass': "// &
                          "'none' (default) or 'uniform_open_ocean' (remove the "// &
                          "melt volume again over uncovered wet cells)", &
                          allowed=[character(len=20) :: "none", "uniform_open_ocean"]))

      call schema%add_group(g)
   end subroutine register_ocean_cavity_melt

   module subroutine register_ocean_forcing(cfg, schema)
      !! `&ocean_forcing_nml`: surface-flux component-set gate (PR-12).
      !! One knob, deliberately — components are filled by fillers, not
      !! by scalar namelist knobs (that would be a knob-per-component,
      !! all dead in v1).
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      logical, pointer :: pl

      g%name = "ocean_forcing"
      g%doc = "Surface-flux component-set gate (PR-12): allocate the "// &
              "q_sw/evap/heat_content_*/... set and derive Q_heat/Q_salt "// &
              "from it every thermo step."
      pl => cfg%ocean%forcing%enable_components
      call g%add(nml_logical("enable_components", pl, &
                             "Allocate the surface-flux component set and run "// &
                             "the assembler (default off => byte-identical)"))
      call schema%add_group(g)
   end subroutine register_ocean_forcing

   module subroutine register_ocean_ice(cfg, schema)
      !! `&ocean_ice_nml`: sea-ice model switch + category/layer counts
      !! (SIS2 port scaffold) + the PR-3c v1 restoring atmospheric-forcing
      !! scalars.  Default off ⇒ byte-identical.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      logical, pointer :: pl
      integer, pointer :: pi
      real(wp), pointer :: pr
      real(wp), pointer :: pra(:)

      g%name = "ocean_ice"
      g%doc = "Sea-ice model (SIS2 port) master switch + category/layer counts."
      pl => cfg%ocean%ice%enable
      call g%add(nml_logical("enable", pl, &
                             "Master switch for the sea-ice slot (default off => byte-identical)"))
      pi => cfg%ocean%ice%ncat
      call g%add(nml_int("ncat", pi, &
                         "Number of ice thickness categories (cat 0 = open water)", min=1))
      pra => cfg%ocean%ice%hlim
      call g%add(nml_real_array("hlim", pra, &
                                "ITD category lower thickness edges; unset => SIS2 default table", &
                                units="m"))
      pi => cfg%ocean%ice%nk_ice
      call g%add(nml_int("nk_ice", pi, &
                         "Vertical ice layers per category (2 = Winton two-layer)", min=1))
      pr => cfg%ocean%ice%air_temp
      call g%add(nml_real("air_temp", pr, &
                          "Prescribed slab-atmosphere air temperature for the v1 restoring filler", &
                          units="degC"))
      pr => cfg%ocean%ice%restore_lambda
      call g%add(nml_real("restore_lambda", pr, &
                          "Surface-flux restoring coefficient dSF/dT (0 = passive column)", &
                          units="W/m^2/K", min=0.0_wp))
      pr => cfg%ocean%ice%sw_down
      call g%add(nml_real("sw_down", pr, &
                          "Downwelling shortwave into the ice top for the v1 filler", &
                          units="W/m^2", min=0.0_wp))
      pr => cfg%ocean%ice%snowfall
      call g%add(nml_real("snowfall", pr, &
                          "Uniform frozen-precipitation rate onto the ice top for the v1 filler", &
                          units="kg/m^2/s", min=0.0_wp))
      pl => cfg%ocean%ice%snow_ice
      call g%add(nml_logical("snow_ice", pl, &
                             "Enable Archimedes snow-ice flooding conversion "// &
                             "(SIS2 SN2IC; default off => byte-identical)"))
      pl => cfg%ocean%ice%transport
      call g%add(nml_logical("transport", pl, &
                             "Enable horizontal category ice/snow transport "// &
                             "(default off => byte-identical)"))
      pi => cfg%ocean%ice%adv_substeps
      call g%add(nml_int("adv_substeps", pi, &
                         "Advective sub-iterations per transport call (SIS2 NSTEPS_ADV)", min=1))
      pr => cfg%ocean%ice%roll_factor
      call g%add(nml_real("roll_factor", pr, &
                          "Thin-ice rolling floor factor (SIS2 SEA_ICE_ROLL_FACTOR); 0 disables rolling", &
                          min=0.0_wp))
      pl => cfg%ocean%ice%dynamics
      call g%add(nml_logical("dynamics", pl, &
                             "Enable C-grid EVP ice dynamics (default off => byte-identical)"))
      pl => cfg%ocean%ice%a_face_stress
      call g%add(nml_logical("a_face_stress", pl, &
                             "Weight EVP wind stress + ice-ocean drag by face ice "// &
                             "concentration (momentum-conserving; default off => byte-identical)"))
      pr => cfg%ocean%ice%p0
      call g%add(nml_real("p0", pr, &
                          "Ice-strength pressure constant (SIS2 ICE_STRENGTH_PSTAR)", &
                          units="Pa", min=0.0_wp))
      pr => cfg%ocean%ice%c0
      call g%add(nml_real("c0", pr, &
                          "Ice-strength exponent constant (SIS2 ICE_STRENGTH_CSTAR)", &
                          min=0.0_wp))
      pr => cfg%ocean%ice%ec
      call g%add(nml_real("ec", pr, &
                          "Yield-curve axis ratio (SIS2 ICE_YIELD_ELLIPTICITY); 0 = cavitating fluid", &
                          min=0.0_wp))
      pr => cfg%ocean%ice%cdw
      call g%add(nml_real("cdw", pr, &
                          "Ice-ocean drag coefficient (SIS2 ICE_CDRAG_WATER)", min=0.0_wp))
      pr => cfg%ocean%ice%rho_ocean
      call g%add(nml_real("rho_ocean", pr, &
                          "Ice-drag reference density (SIS2 RHO_OCEAN)", &
                          units="kg/m^3", min=0.0_wp))
      pi => cfg%ocean%ice%evp_sub_steps
      call g%add(nml_int("evp_sub_steps", pi, &
                         "EVP subcycles per slow step (SIS2 NSTEPS_DYN)", min=1))
      pr => cfg%ocean%ice%del_sh_min_scale
      call g%add(nml_real("del_sh_min_scale", pr, &
                          "Viscosity-floor scale (SIS2 ICE_DEL_SH_MIN_SCALE)", min=0.0_wp))
      pr => cfg%ocean%ice%tdamp
      call g%add(nml_real("tdamp", pr, &
                          "Elastic damping timescale rule (SIS2 ICE_TDAMP_ELASTIC): "// &
                          ">0 seconds, ==0 auto (0.2*dt_slow), <0 fraction of dt_slow"))
      pr => cfg%ocean%ice%cfl_trunc
      call g%add(nml_real("cfl_trunc", pr, &
                          "Transport-CFL ceiling on the final ice velocity "// &
                          "(SIS2 CFL_TRUNCATE, default 0.5 there); 0 disables", min=0.0_wp))
      pl => cfg%ocean%ice%cfl_trunc_dyn_its
      call g%add(nml_logical("cfl_trunc_dyn_its", pl, &
                             "Also clip the ice velocity every EVP subcycle "// &
                             "(SIS2 CFL_TRUNC_DYN_ITS)"))
      pl => cfg%ocean%ice%project_ci
      call g%add(nml_logical("project_ci", pl, &
                             "Project ice concentration forward within the EVP "// &
                             "subcycle loop and recompute the ice strength "// &
                             "(SIS2 PROJECT_ICE_CONCENTRATION)"))
      call schema%add_group(g)
   end subroutine register_ocean_ice

   module subroutine register_ocean_ice_ic(cfg, schema)
      !! `&ocean_ice_ic_nml`: sea-ice ANALYTIC initial-condition path (PR
      !! 24). Default `conc_config="zero"` ⇒ byte-identical. `"file"` is
      !! deliberately NOT in the `allowed=` list — file-backed ICs are
      !! PR-14 (v1.1), out of scope here.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      character(len=:), pointer :: ps
      real(wp), pointer :: pr

      g%name = "ocean_ice_ic"
      g%doc = "Sea-ice analytic initial-condition path (PR 24)."
      ps => cfg%ocean%ice_ic%conc_config
      call g%add(nml_enum("conc_config", ps, &
                          "Sea-ice initial concentration configuration", &
                          allowed=[character(len=12) :: "zero", "uniform", "latitudes"]))
      pr => cfg%ocean%ice_ic%conc
      call g%add(nml_real("conc", pr, &
                          "Uniform-mode concentration ('uniform' only)", &
                          min=0.0_wp, max=1.0_wp))
      pr => cfg%ocean%ice_ic%h_ice
      call g%add(nml_real("h_ice", pr, &
                          "Ice thickness where seeded (SIS2 ICE_INIT_MASS as a thickness)", &
                          units="m", min=0.0_wp))
      pr => cfg%ocean%ice_ic%h_snow
      call g%add(nml_real("h_snow", pr, &
                          "Snow thickness where seeded (SIS2 SNOW_INIT_MASS as a thickness)", &
                          units="m", min=0.0_wp))
      pr => cfg%ocean%ice_ic%t_ice
      call g%add(nml_real("t_ice", pr, &
                          "Ice/snow temperature fed through the exact enthalpy inversion "// &
                          "(SIS2 ICE_TEMPERATURE_IC)", units="degC"))
      pr => cfg%ocean%ice_ic%s_ice
      call g%add(nml_real("s_ice", pr, &
                          "Ice bulk salinity (SIS2 ICE_SALINITY_IC)", &
                          units="PSU", min=0.0_wp))
      pr => cfg%ocean%ice_ic%arctic_edge
      call g%add(nml_real("arctic_edge", pr, &
                          "'latitudes' Arctic ice edge (SIS2 ARCTIC_ICE_EDGE_IC)", &
                          units="degrees_north", min=-91.0_wp, max=91.0_wp))
      pr => cfg%ocean%ice_ic%antarctic_edge
      call g%add(nml_real("antarctic_edge", pr, &
                          "'latitudes' Antarctic ice edge (SIS2 ANTARCTIC_ICE_EDGE_IC)", &
                          units="degrees_north", min=-91.0_wp, max=91.0_wp))
      call schema%add_group(g)
   end subroutine register_ocean_ice_ic

   module subroutine register_ocean_restore(cfg, schema)
      !! `&ocean_restore_nml`: surface buoyancy restoring (MOM6
      !! `RESTOREBUOY`).  Piston-velocity relaxation of top-layer T / S
      !! toward scalar targets.  Default OFF ⇒ bit-identical.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      logical, pointer :: pl
      real(wp), pointer :: pr

      g%name = "ocean_restore"
      g%doc = "Surface buoyancy restoring (MOM6 RESTOREBUOY): "// &
              "piston-velocity relaxation of top-layer T / S toward targets."
      pl => cfg%ocean%restore%enable_restore_temp
      call g%add(nml_logical("enable_restore_temp", pl, &
                             "Master switch for SST restoring"))
      pl => cfg%ocean%restore%enable_restore_salt
      call g%add(nml_logical("enable_restore_salt", pl, &
                             "Master switch for SSS restoring"))
      pr => cfg%ocean%restore%piston_t
      call g%add(nml_real("piston_t", pr, &
                          "SST piston velocity (MOM6 FLUXCONST_T)", &
                          units="m/day", min=0.0_wp))
      pr => cfg%ocean%restore%piston_s
      call g%add(nml_real("piston_s", pr, &
                          "SSS piston velocity (MOM6 FLUXCONST_S)", &
                          units="m/day", min=0.0_wp))
      pr => cfg%ocean%restore%restore_sst
      call g%add(nml_real("restore_sst", pr, &
                          "Scalar target SST", units="degC"))
      pr => cfg%ocean%restore%restore_sss
      call g%add(nml_real("restore_sss", pr, &
                          "Scalar target SSS", units="PSU"))
      call schema%add_group(g)
   end subroutine register_ocean_restore

   module subroutine register_ocean_geothermal(cfg, schema)
      !! `&ocean_geothermal_nml`: geothermal bottom-heat-flux switch +
      !! scalar flux.  Bed-side analogue of the surface heat flux.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      logical, pointer :: pl
      real(wp), pointer :: pr

      g%name = "ocean_geothermal"
      g%doc = "Geothermal bottom heat flux (bed-side analogue of the surface heat flux)."
      pl => cfg%ocean%geothermal%enable
      call g%add(nml_logical("enable", pl, &
                             "Apply a constant geothermal bottom heat flux to the bed layer"))
      pr => cfg%ocean%geothermal%q_geo
      call g%add(nml_real("q_geo", pr, &
                          "Constant bottom heat flux (positive into the ocean from below)", &
                          units="W/m^2"))
      call schema%add_group(g)
   end subroutine register_ocean_geothermal

   module subroutine register_ocean_sponge(cfg, schema)
      !! `&ocean_sponge_nml`: the map-driven sponge (PR-23). Default
      !! `enable = .false.` ⇒ the legacy `&ocean_bc_nml` band kernels run
      !! unchanged ⇒ bit-identical. No `damp_max` / `idamp_file` / `*_var`
      !! keys in v1 — those validate and do nothing until PR-23b adds their
      !! consumer (CLAUDE.md trap class "dead knobs").
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      logical, pointer :: pl
      real(wp), pointer :: pr
      integer, pointer :: pi
      character(len=:), pointer :: ps

      g%name = "ocean_sponge"
      g%doc = "Map-driven sponge: per-cell Idamp + 3-D reference state (PR-23)."
      pl => cfg%ocean%sponge%enable
      call g%add(nml_logical("enable", pl, &
                             "Master switch (default off; legacy band kernels run when off)"))
      ps => cfg%ocean%sponge%damp_source
      call g%add(nml_enum("damp_source", ps, &
                          "How idamp_h/u/v are filled", &
                          allowed=[character(len=4) :: "band", "file"]))
      ps => cfg%ocean%sponge%target_source
      call g%add(nml_enum("target_source", ps, &
                          "Reference-state source", &
                          allowed=[character(len=8) :: "ic", "linear_z", "file"]))
      ps => cfg%ocean%sponge%ramp
      call g%add(nml_enum("ramp", ps, &
                          "Band ramp shape from the sponge wall inward "// &
                          "('linear' is ISOMIP+ Eq. 20)", &
                          allowed=[character(len=6) :: "cosine", "linear"]))
      pr => cfg%ocean%sponge%lin_t_ref
      call g%add(nml_real("lin_t_ref", pr, &
                          "target_source='linear_z': T at the z = 0 datum", units="degC"))
      pr => cfg%ocean%sponge%lin_dt_dz
      call g%add(nml_real("lin_dt_dz", pr, &
                          "target_source='linear_z': dT/dz, z positive UP "// &
                          "(stable => > 0)", units="degC/m"))
      pr => cfg%ocean%sponge%lin_s_ref
      call g%add(nml_real("lin_s_ref", pr, &
                          "target_source='linear_z': S at the z = 0 datum", units="PSU"))
      pr => cfg%ocean%sponge%lin_ds_dz
      call g%add(nml_real("lin_ds_dz", pr, &
                          "target_source='linear_z': dS/dz, z positive UP "// &
                          "(stable => < 0)", units="PSU/m"))
      pl => cfg%ocean%sponge%relax_uv
      call g%add(nml_logical("relax_uv", pl, &
                             "Relax u/v toward u_ref/v_ref"))
      pl => cfg%ocean%sponge%relax_tracers
      call g%add(nml_logical("relax_tracers", pl, &
                             "Relax every registered tracer toward its 3-D reference field"))
      pl => cfg%ocean%sponge%relax_h
      call g%add(nml_logical("relax_h", pl, &
                             "Interior-interface thickness damping (NOT IMPLEMENTED in v1 — "// &
                             "deferred to PR-23b, requires vcoord_type='lagrangian')"))
      pi => cfg%ocean%sponge%west_width
      call g%add(nml_int("west_width", pi, &
                         "Per-edge band-width override, cells (<0 => inherit "// &
                         "&ocean_bc_nml sponge_width)"))
      pi => cfg%ocean%sponge%east_width
      call g%add(nml_int("east_width", pi, "as west_width"))
      pi => cfg%ocean%sponge%south_width
      call g%add(nml_int("south_width", pi, "as west_width"))
      pi => cfg%ocean%sponge%north_width
      call g%add(nml_int("north_width", pi, "as west_width"))
      pr => cfg%ocean%sponge%west_strength
      call g%add(nml_real("west_strength", pr, &
                          "Per-edge peak relaxation-rate override (<0 => inherit "// &
                          "&ocean_bc_nml sponge_strength)", units="1/s"))
      pr => cfg%ocean%sponge%east_strength
      call g%add(nml_real("east_strength", pr, "as west_strength", units="1/s"))
      pr => cfg%ocean%sponge%south_strength
      call g%add(nml_real("south_strength", pr, "as west_strength", units="1/s"))
      pr => cfg%ocean%sponge%north_strength
      call g%add(nml_real("north_strength", pr, "as west_strength", units="1/s"))
      call schema%add_group(g)
   end subroutine register_ocean_sponge

   module subroutine register_ocean_tracers(cfg, schema)
      !! `&ocean_tracers_nml`: prognostic-tracer registry switches.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      logical, pointer :: pl
      real(wp), pointer :: pr

      g%name = "ocean_tracers"
      g%doc = "Prognostic-tracer registry switches."
      pl => cfg%ocean%tracers%enable_ideal_age
      call g%add(nml_logical("enable_ideal_age", pl, &
                             "Register a passive ideal-age tracer"))
      pr => cfg%ocean%tracers%ideal_age_young_val
      call g%add(nml_real("ideal_age_young_val", pr, &
                          "Surface-band ideal-age Dirichlet value (0 = today's hard-coded reset)", &
                          units="s"))
      pr => cfg%ocean%tracers%ideal_age_sfc_growth_rate
      call g%add(nml_real("ideal_age_sfc_growth_rate", pr, &
                          "Exponential growth rate of the ideal-age surface value (0 = constant)", &
                          units="1/s"))
      pl => cfg%ocean%tracers%enable_pseudo_salt
      call g%add(nml_logical("enable_pseudo_salt", pl, &
                             "Register the pseudo-salt verification tracer "// &
                             "(diagnostic; seeded to S, given salinity's boundary fluxes)"))
      call schema%add_group(g)
   end subroutine register_ocean_tracers

   module subroutine register_ocean_bc(cfg, schema)
      !! `&ocean_bc_nml`: open-boundary condition config (P4.5 migration
      !! off the hand-rolled `read_ocean_bc_nml`).  All defaults reproduce
      !! a closed-wall run (bit-identical to nmls that omit this block).
      !! The `west`/`east`/`south`/`north` edge-type enum lists every arm
      !! of `ocean_bc_type_from_string` (rdb_ocean_boundary_types.F90),
      !! not just the subset named in the config-type docstring — that
      !! function is what the parsed string is actually fed to, and the
      !! native reader accepted (and `validate_config` did not reject) any
      !! of its arms on any edge, so the schema must not narrow that.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      real(wp), pointer :: pr
      real(wp), pointer :: pra(:)
      integer, pointer :: pi
      logical, pointer :: pl
      character(len=:), pointer :: ps
      type(nml_real_array_key_t) :: rr_key
         !! Local materialisation for every `nml_real_array(...)` result
         !! before `g%add` — NVHPC 26.5 corrupts unrelated heap state when
         !! a 4th-or-later real_array function-result temporary is passed
         !! directly to `g%add(...)` in the same subroutine (bisected in
         !! `register_ocean_diag`'s comment above; this group registers
         !! twelve, so every one uses the workaround).
      character(len=*), parameter :: edge_allowed(*) = &
                                     [character(len=13) :: "wall", "open", "tidal", "nested", "inflow", &
                                                           "discharge", "clamped", "sponge", "chapman", "periodic", "tripolar_fold"]

      g%name = "ocean_bc"
      g%doc = "Open-boundary condition config: per-edge BC type, clamped/inflow "// &
              "Dirichlet values, sponge band, Orlanski radiation, full Flather, "// &
              "per-edge tidal constituents, tidal-OBC nodal correction, "// &
              "wall-velocity masking."

      ps => cfg%ocean%bc%west
      call g%add(nml_enum("west", ps, "West edge boundary type", allowed=edge_allowed))
      ps => cfg%ocean%bc%east
      call g%add(nml_enum("east", ps, "East edge boundary type", allowed=edge_allowed))
      ps => cfg%ocean%bc%south
      call g%add(nml_enum("south", ps, "South edge boundary type", allowed=edge_allowed))
      ps => cfg%ocean%bc%north
      call g%add(nml_enum("north", ps, &
                          "North edge boundary type ('tripolar_fold' is meaningful "// &
                          "here only — Murray bipolar cap; requires "// &
                          "grid_config='tripolar' + periodic west/east)", &
                          allowed=edge_allowed))

      pr => cfg%ocean%bc%west_clamped_eta
      call g%add(nml_real("west_clamped_eta", pr, "Clamped SSH, west edge", units="m"))
      pr => cfg%ocean%bc%east_clamped_eta
      call g%add(nml_real("east_clamped_eta", pr, "Clamped SSH, east edge", units="m"))
      pr => cfg%ocean%bc%south_clamped_eta
      call g%add(nml_real("south_clamped_eta", pr, "Clamped SSH, south edge", units="m"))
      pr => cfg%ocean%bc%north_clamped_eta
      call g%add(nml_real("north_clamped_eta", pr, "Clamped SSH, north edge", units="m"))
      pr => cfg%ocean%bc%west_clamped_u
      call g%add(nml_real("west_clamped_u", pr, "Clamped normal velocity, west edge", units="m/s"))
      pr => cfg%ocean%bc%east_clamped_u
      call g%add(nml_real("east_clamped_u", pr, "Clamped normal velocity, east edge", units="m/s"))
      pr => cfg%ocean%bc%south_clamped_v
      call g%add(nml_real("south_clamped_v", pr, "Clamped normal velocity, south edge", units="m/s"))
      pr => cfg%ocean%bc%north_clamped_v
      call g%add(nml_real("north_clamped_v", pr, "Clamped normal velocity, north edge", units="m/s"))

      pr => cfg%ocean%bc%west_inflow_S
      call g%add(nml_real("west_inflow_S", pr, "Inflow salinity, west edge", units="PSU"))
      pr => cfg%ocean%bc%west_inflow_T
      call g%add(nml_real("west_inflow_T", pr, "Inflow temperature, west edge", units="degC"))
      pr => cfg%ocean%bc%east_inflow_S
      call g%add(nml_real("east_inflow_S", pr, "Inflow salinity, east edge", units="PSU"))
      pr => cfg%ocean%bc%east_inflow_T
      call g%add(nml_real("east_inflow_T", pr, "Inflow temperature, east edge", units="degC"))
      pr => cfg%ocean%bc%south_inflow_S
      call g%add(nml_real("south_inflow_S", pr, "Inflow salinity, south edge", units="PSU"))
      pr => cfg%ocean%bc%south_inflow_T
      call g%add(nml_real("south_inflow_T", pr, "Inflow temperature, south edge", units="degC"))
      pr => cfg%ocean%bc%north_inflow_S
      call g%add(nml_real("north_inflow_S", pr, "Inflow salinity, north edge", units="PSU"))
      pr => cfg%ocean%bc%north_inflow_T
      call g%add(nml_real("north_inflow_T", pr, "Inflow temperature, north edge", units="degC"))

      pi => cfg%ocean%bc%sponge_width
      call g%add(nml_int("sponge_width", pi, "Sponge band width in cells (0 disables)"))
      pr => cfg%ocean%bc%sponge_strength
      call g%add(nml_real("sponge_strength", pr, &
                          "Peak sponge relaxation rate at the outer face of the band", &
                          units="1/s"))
      pl => cfg%ocean%bc%sponge_relax_tracers
      call g%add(nml_logical("sponge_relax_tracers", pl, &
                             "Extend the sponge relaxation to h_layer + tracers"))

      pr => cfg%ocean%bc%res_lscale_out
      call g%add(nml_real("res_lscale_out", pr, &
                          "Outflow reservoir length scale (0 = disabled)", units="m"))
      pr => cfg%ocean%bc%res_lscale_in
      call g%add(nml_real("res_lscale_in", pr, &
                          "Inflow reservoir length scale (0 = instantaneous inflow)", units="m"))

      ps => cfg%ocean%bc%radiation_scheme
      call g%add(nml_enum("radiation_scheme", ps, &
                          "'anomaly' = v1 BT-mean + zero-gradient anomaly (default); "// &
                          "'orlanski' = per-layer implicit-upwind radiation (Orlanski 1976)", &
                          allowed=[character(len=8) :: "anomaly", "orlanski"]))
      pr => cfg%ocean%bc%orlanski_rx_max
      call g%add(nml_real("orlanski_rx_max", pr, &
                          "Clamp on the Orlanski nondimensional phase speed", &
                          units="grid cells/step"))
      pr => cfg%ocean%bc%orlanski_gamma
      call g%add(nml_real("orlanski_gamma", pr, &
                          "Running-mean weight (0=full running mean, 1=instant)"))
      pr => cfg%ocean%bc%nudge_tau_in
      call g%add(nml_real("nudge_tau_in", pr, "Inflow nudging timescale (0 = off)", units="s"))
      pr => cfg%ocean%bc%nudge_tau_out
      call g%add(nml_real("nudge_tau_out", pr, "Outflow nudging timescale (0 = off)", units="s"))

      ps => cfg%ocean%bc%flather_form
      call g%add(nml_enum("flather_form", ps, &
                          "'legacy' = v1 Flather (u_ext=0, no interior vel, default); "// &
                          "'full' = half-characteristic form (Flather 1976)", &
                          allowed=[character(len=6) :: "legacy", "full"]))
      pr => cfg%ocean%bc%west_ext_u
      call g%add(nml_real("west_ext_u", pr, "Exterior barotropic u, west edge", units="m/s"))
      pr => cfg%ocean%bc%east_ext_u
      call g%add(nml_real("east_ext_u", pr, "Exterior barotropic u, east edge", units="m/s"))
      pr => cfg%ocean%bc%south_ext_v
      call g%add(nml_real("south_ext_v", pr, "Exterior barotropic v, south edge", units="m/s"))
      pr => cfg%ocean%bc%north_ext_v
      call g%add(nml_real("north_ext_v", pr, "Exterior barotropic v, north edge", units="m/s"))

      pi => cfg%ocean%bc%west_n_tidal
      call g%add(nml_int("west_n_tidal", pi, "Active tidal constituents, west edge", &
                         min=0, max=OBC_MAX_TIDAL_CFG))
      pra => cfg%ocean%bc%west_tidal_amp
      rr_key = nml_real_array("west_tidal_amp", pra, "Constituent amplitudes, west edge", &
                              units="m")
      call g%add(rr_key)
      pra => cfg%ocean%bc%west_tidal_phase
      rr_key = nml_real_array("west_tidal_phase", pra, "Constituent phases, west edge", &
                              units="rad")
      call g%add(rr_key)
      pra => cfg%ocean%bc%west_tidal_omega
      rr_key = nml_real_array("west_tidal_omega", pra, &
                              "Constituent angular frequencies, west edge", units="rad/s")
      call g%add(rr_key)

      pi => cfg%ocean%bc%east_n_tidal
      call g%add(nml_int("east_n_tidal", pi, "Active tidal constituents, east edge", &
                         min=0, max=OBC_MAX_TIDAL_CFG))
      pra => cfg%ocean%bc%east_tidal_amp
      rr_key = nml_real_array("east_tidal_amp", pra, "Constituent amplitudes, east edge", &
                              units="m")
      call g%add(rr_key)
      pra => cfg%ocean%bc%east_tidal_phase
      rr_key = nml_real_array("east_tidal_phase", pra, "Constituent phases, east edge", &
                              units="rad")
      call g%add(rr_key)
      pra => cfg%ocean%bc%east_tidal_omega
      rr_key = nml_real_array("east_tidal_omega", pra, &
                              "Constituent angular frequencies, east edge", units="rad/s")
      call g%add(rr_key)

      pi => cfg%ocean%bc%south_n_tidal
      call g%add(nml_int("south_n_tidal", pi, "Active tidal constituents, south edge", &
                         min=0, max=OBC_MAX_TIDAL_CFG))
      pra => cfg%ocean%bc%south_tidal_amp
      rr_key = nml_real_array("south_tidal_amp", pra, "Constituent amplitudes, south edge", &
                              units="m")
      call g%add(rr_key)
      pra => cfg%ocean%bc%south_tidal_phase
      rr_key = nml_real_array("south_tidal_phase", pra, "Constituent phases, south edge", &
                              units="rad")
      call g%add(rr_key)
      pra => cfg%ocean%bc%south_tidal_omega
      rr_key = nml_real_array("south_tidal_omega", pra, &
                              "Constituent angular frequencies, south edge", units="rad/s")
      call g%add(rr_key)

      pi => cfg%ocean%bc%north_n_tidal
      call g%add(nml_int("north_n_tidal", pi, "Active tidal constituents, north edge", &
                         min=0, max=OBC_MAX_TIDAL_CFG))
      pra => cfg%ocean%bc%north_tidal_amp
      rr_key = nml_real_array("north_tidal_amp", pra, "Constituent amplitudes, north edge", &
                              units="m")
      call g%add(rr_key)
      pra => cfg%ocean%bc%north_tidal_phase
      rr_key = nml_real_array("north_tidal_phase", pra, "Constituent phases, north edge", &
                              units="rad")
      call g%add(rr_key)
      pra => cfg%ocean%bc%north_tidal_omega
      rr_key = nml_real_array("north_tidal_omega", pra, &
                              "Constituent angular frequencies, north edge", units="rad/s")
      call g%add(rr_key)

      pl => cfg%ocean%bc%obc_tidal_nodal
      call g%add(nml_logical("obc_tidal_nodal", pl, &
                             "Apply the 18.6-yr nodal factor + equilibrium/nodal phase to "// &
                             "the open-boundary tidal elevation forcing (shares the "// &
                             "&ocean_tides_nml astro generator); when true *_tidal_phase "// &
                             "becomes a Greenwich phase lag"))
      pl => cfg%ocean%bc%mask_wall_velocity
      call g%add(nml_logical("mask_wall_velocity", pl, &
                             "Zero the T-cell wet-mask in the ghost cells beyond every "// &
                             "solid WALL edge at setup, so wall-normal C-grid face "// &
                             "velocities are masked rather than left as garbage "// &
                             "(default ON; set .false. to reproduce a pre-fix legacy "// &
                             "closed-basin baseline)"))

      call schema%add_group(g)
   end subroutine register_ocean_bc

   module subroutine register_ocean_atm_forcing(cfg, schema)
      !! `&ocean_atm_forcing_nml` (PR-C2c): JRA55-do-style prescribed
      !! atmospheric state.  See `ocean_atm_forcing_config_t` for the
      !! per-knob docs.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      logical, pointer :: pl
      integer, pointer :: pi
      real(wp), pointer :: pr
      character(len=:), pointer :: ps

      g%name = "ocean_atm_forcing"
      g%doc = "JRA55-do-style prescribed atmospheric state driving "// &
              "bulk-flux ocean/ice surface forcing (PR-C2c)."

      pl => cfg%ocean%atm_forcing%enable
      call g%add(nml_logical("enable", pl, "Master switch"))
      ps => cfg%ocean%atm_forcing%source
      call g%add(nml_enum("source", ps, "Atmospheric-state provider", &
                          allowed=[character(len=16) :: "jra55do"]))

      ps => cfg%ocean%atm_forcing%filelist_tas
      call g%add(nml_string("filelist_tas", ps, &
                            "tools/jra_filelist.py-style file list for tas "// &
                            "(air temperature); blank = fail loud when enable=.true."))
      ps => cfg%ocean%atm_forcing%filelist_huss
      call g%add(nml_string("filelist_huss", ps, "File list for huss (specific humidity)"))
      ps => cfg%ocean%atm_forcing%filelist_uas
      call g%add(nml_string("filelist_uas", ps, "File list for uas (eastward 10 m wind)"))
      ps => cfg%ocean%atm_forcing%filelist_vas
      call g%add(nml_string("filelist_vas", ps, "File list for vas (northward 10 m wind)"))
      ps => cfg%ocean%atm_forcing%filelist_psl
      call g%add(nml_string("filelist_psl", ps, "File list for psl (sea-level pressure)"))
      ps => cfg%ocean%atm_forcing%filelist_rlds
      call g%add(nml_string("filelist_rlds", ps, "File list for rlds (downwelling longwave)"))
      ps => cfg%ocean%atm_forcing%filelist_rsds
      call g%add(nml_string("filelist_rsds", ps, &
                            "File list for rsds (downwelling shortwave, interval mean)"))
      ps => cfg%ocean%atm_forcing%filelist_prra
      call g%add(nml_string("filelist_prra", ps, "File list for prra (liquid precip rate)"))
      ps => cfg%ocean%atm_forcing%filelist_prsn
      call g%add(nml_string("filelist_prsn", ps, "File list for prsn (frozen precip rate)"))

      ps => cfg%ocean%atm_forcing%var_tas
      call g%add(nml_string("var_tas", ps, "Variable name inside the tas filelist"))
      ps => cfg%ocean%atm_forcing%var_huss
      call g%add(nml_string("var_huss", ps, "Variable name inside the huss filelist"))
      ps => cfg%ocean%atm_forcing%var_uas
      call g%add(nml_string("var_uas", ps, "Variable name inside the uas filelist"))
      ps => cfg%ocean%atm_forcing%var_vas
      call g%add(nml_string("var_vas", ps, "Variable name inside the vas filelist"))
      ps => cfg%ocean%atm_forcing%var_psl
      call g%add(nml_string("var_psl", ps, "Variable name inside the psl filelist"))
      ps => cfg%ocean%atm_forcing%var_rlds
      call g%add(nml_string("var_rlds", ps, "Variable name inside the rlds filelist"))
      ps => cfg%ocean%atm_forcing%var_rsds
      call g%add(nml_string("var_rsds", ps, "Variable name inside the rsds filelist"))
      ps => cfg%ocean%atm_forcing%var_prra
      call g%add(nml_string("var_prra", ps, "Variable name inside the prra filelist"))
      ps => cfg%ocean%atm_forcing%var_prsn
      call g%add(nml_string("var_prsn", ps, "Variable name inside the prsn filelist"))

      ps => cfg%ocean%atm_forcing%time_mode
      call g%add(nml_enum("time_mode", ps, &
                          "interannual (default) or repeat_year (RYF)", &
                          allowed=[character(len=16) :: "interannual", "repeat_year"]))
      ps => cfg%ocean%atm_forcing%ryf_start_date
      call g%add(nml_string("ryf_start_date", ps, &
                            "RYF window start date (time_mode=repeat_year only)"))
      pr => cfg%ocean%atm_forcing%ryf_period_days
      call g%add(nml_real("ryf_period_days", pr, &
                          "RYF window length (time_mode=repeat_year only)", units="days"))

      ps => cfg%ocean%atm_forcing%states_regrid
      call g%add(nml_enum("states_regrid", ps, &
                          "Horizontal regrid for tas/huss/uas/vas/psl; only "// &
                          "'bilinear' is implemented -- 'patch' fails loud "// &
                          "naming the xESMF gap", &
                          allowed=[character(len=16) :: "bilinear", "patch"]))
      ps => cfg%ocean%atm_forcing%fluxes_regrid
      call g%add(nml_enum("fluxes_regrid", ps, &
                          "Horizontal regrid for rlds/rsds/prra/prsn", &
                          allowed=[character(len=16) :: "conservative", "bilinear"]))
      ps => cfg%ocean%atm_forcing%sw_interp
      call g%add(nml_enum("sw_interp", ps, &
                          "coszen (default, held interval + diurnal cosine-zenith "// &
                          "factor) or linear (plain time-blend)", &
                          allowed=[character(len=16) :: "coszen", "linear"]))
      pr => cfg%ocean%atm_forcing%sw_offset_s
      call g%add(nml_real("sw_offset_s", pr, &
                          "Time offset applied to the rsds query time before the "// &
                          "bracket search (sw_interp=coszen)", units="s"))
      pi => cfg%ocean%atm_forcing%orb_year
      call g%add(nml_int("orb_year", pi, &
                         "Orbital-parameter year (reserved; v1 coszen uses fixed "// &
                         "Spencer (1971) coefficients regardless)"))
      pl => cfg%ocean%atm_forcing%wind_relative
      call g%add(nml_logical("wind_relative", pl, &
                             "Subtract the ocean surface current from the "// &
                             "atmospheric wind before the bulk flux"))
      pr => cfg%ocean%atm_forcing%ocean_albedo
      call g%add(nml_real("ocean_albedo", pr, &
                          "Ocean surface shortwave albedo (Large & Yeager 2009)"))
      pl => cfg%ocean%atm_forcing%verbose
      call g%add(nml_logical("verbose", pl, &
                             "Log every bracket advance of the 9 underlying readers"))

      call schema%add_group(g)
   end subroutine register_ocean_atm_forcing

   module subroutine register_ocean_bulk_flux(cfg, schema)
      !! `&ocean_bulk_flux_nml`: Monin-Obukhov bulk-flux kernel knobs
      !! (`bulk_flux_config_t`, `rdb_ocean_bulk_flux`).  `enable` is
      !! implied `.true.` whenever `&ocean_atm_forcing_nml enable` is --
      !! `validate_config` fails loud on an explicit `.false.` conflict.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      logical, pointer :: pl
      integer, pointer :: pi
      real(wp), pointer :: pr

      g%name = "ocean_bulk_flux"
      g%doc = "Monin-Obukhov bulk-flux kernel (Large & Yeager 2004/2009)."

      pl => cfg%ocean%bulk_flux%enable
      call g%add(nml_logical("enable", pl, &
                             "Master switch (implied by &ocean_atm_forcing_nml "// &
                             "enable; an explicit .false. with atm_forcing on "// &
                             "is a conflict)"))
      pi => cfg%ocean%bulk_flux%n_iter
      call g%add(nml_int("n_iter", pi, "Monin-Obukhov iteration count", min=1))
      pr => cfg%ocean%bulk_flux%z_wind
      call g%add(nml_real("z_wind", pr, "Wind reference height", units="m"))
      pr => cfg%ocean%bulk_flux%z_ta
      call g%add(nml_real("z_ta", pr, &
                          "Air temperature / specific humidity reference height", &
                          units="m"))
      pr => cfg%ocean%bulk_flux%emissivity
      call g%add(nml_real("emissivity", pr, "Longwave emissivity of seawater"))

      call schema%add_group(g)
   end subroutine register_ocean_bulk_flux

end submodule rdb_config_register_forcing
