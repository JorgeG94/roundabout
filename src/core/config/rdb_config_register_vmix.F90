submodule(rdb_config) rdb_config_register_vmix
   !! Strict-schema registration: vertical-mixing + EOS/PGF/drag groups.
   implicit none

contains

   module subroutine register_kappa_shear(cfg, schema)
      !! `&ocean_kappa_shear` (JHL08 shear-driven interior mixing).
      !! Exclusive `> 0` bounds (ri_crit, kappa_0, tol_err, prandtl_turb)
      !! are left to the configure-time strict checks.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      logical, pointer :: pl
      real(wp), pointer :: pr
      integer, pointer :: pi

      g%name = "ocean_kappa_shear"
      g%doc = "JHL08 shear-driven interior turbulence (kappa-shear)."

      pl => cfg%ocean%kshear%enable
      call g%add(nml_logical("enable", pl, &
                             "Master switch (requires use_closure + thermodynamics)"))
      pr => cfg%ocean%kshear%ri_crit
      call g%add(nml_real("ri_crit", pr, "Critical Richardson number", &
                          units="nondim"))
      pr => cfg%ocean%kshear%shearmix_rate
      call g%add(nml_real("shearmix_rate", pr, "Shear source-rate coefficient"))
      pr => cfg%ocean%kshear%fri_curvature
      call g%add(nml_real("fri_curvature", pr, "Ri-function curvature in the shear source"))
      pr => cfg%ocean%kshear%c_n
      call g%add(nml_real("c_n", pr, "TKE decay-rate coefficient vs stratification N"))
      pr => cfg%ocean%kshear%c_s
      call g%add(nml_real("c_s", pr, "TKE decay-rate coefficient vs shear S"))
      pr => cfg%ocean%kshear%lambda
      call g%add(nml_real("lambda", pr, "Buoyancy mixing-length-scale coefficient"))
      pr => cfg%ocean%kshear%lz_rescale
      call g%add(nml_real("lz_rescale", pr, "Boundary-distance length-scale rescale factor"))
      pr => cfg%ocean%kshear%kappa_0
      call g%add(nml_real("kappa_0", pr, "Background diffusivity (pre-step kappa)", &
                          units="m^2/s"))
      pr => cfg%ocean%kshear%kappa_seed
      call g%add(nml_real("kappa_seed", pr, "Iteration seed diffusivity", units="m^2/s"))
      pr => cfg%ocean%kshear%kappa_trunc
      call g%add(nml_real("kappa_trunc", pr, "Diffusivity truncated to 0 below this", &
                          units="m^2/s"))
      pr => cfg%ocean%kshear%tke_bg
      call g%add(nml_real("tke_bg", pr, "Background TKE (Q denominator floor)", &
                          units="m^2/s^2"))
      pr => cfg%ocean%kshear%tol_err
      call g%add(nml_real("tol_err", pr, "Picard convergence tolerance"))
      pi => cfg%ocean%kshear%max_inner_it
      call g%add(nml_int("max_inner_it", pi, "Inner Picard iteration cap", min=1))
      pi => cfg%ocean%kshear%max_substep_it
      call g%add(nml_int("max_substep_it", pi, "Outer adaptive-substep iteration cap", min=1))
      pr => cfg%ocean%kshear%src_max_chg
      call g%add(nml_real("src_max_chg", pr, "Adaptive-dt source-change tolerance band"))
      pr => cfg%ocean%kshear%prandtl_turb
      call g%add(nml_real("prandtl_turb", pr, "Kv = prandtl_turb * Kd into the momentum solve"))
      pr => cfg%ocean%kshear%vel_underflow
      call g%add(nml_real("vel_underflow", pr, "Velocity snap-to-zero magnitude", &
                          units="m/s"))
      pl => cfg%ocean%kshear%massless_merge
      call g%add(nml_logical("massless_merge", pl, &
                             "Merge vanished (<H_VANISHED) layers onto the "// &
                             "massive sub-grid before the column solve "// &
                             "(default off; identity columns bypass)"))
      pl => cfg%ocean%kshear%at_vertex
      call g%add(nml_logical("at_vertex", pl, &
                             "Solve the JHL08 columns at C-grid corners "// &
                             "(vorticity points) from the native face "// &
                             "velocities, averaging corner Kd back to tracer "// &
                             "points (MOM6 VERTEX_SHEAR; v1 = Kd only)"))
      pl => cfg%ocean%kshear%vertex_geometric_mean
      call g%add(nml_logical("vertex_geometric_mean", pl, &
                             "Geometric (vs arithmetic) mean in the "// &
                             "corner->centre Kd average (MOM6 "// &
                             "VERTEX_SHEAR_GEOMETRIC_MEAN)"))
      pr => cfg%ocean%kshear%vertex_geomean_kdmin
      call g%add(nml_real("vertex_geomean_kdmin", pr, &
                          "Floor applied to each corner Kd before the "// &
                          "geometric mean (inert unless vertex_geometric_mean; "// &
                          "OM5 configs use 1e-9)", &
                          units="m^2/s", min=0.0_wp))

      call schema%add_group(g)
   end subroutine register_kappa_shear

   module subroutine register_tidal_mixing(cfg, schema)
      !! `&ocean_tidal_mixing` (St-Laurent/Simmons internal-tide mixing).
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      logical, pointer :: pl
      real(wp), pointer :: pr
      character(len=:), pointer :: ps

      g%name = "ocean_tidal_mixing"
      g%doc = "St-Laurent/Simmons internal-tide interior diapycnal mixing."

      pl => cfg%ocean%tidal_mixing%enable
      call g%add(nml_logical("enable", pl, &
                             "Master switch (requires use_closure + thermodynamics)"))
      pr => cfg%ocean%tidal_mixing%gamma
      call g%add(nml_real("gamma", pr, "Local-dissipation fraction q (GAMMA_ITIDES)", &
                          units="nondim"))
      pr => cfg%ocean%tidal_mixing%mu
      call g%add(nml_real("mu", pr, "Mixing efficiency Gamma_mix (MU_ITIDES)", &
                          units="nondim"))
      pr => cfg%ocean%tidal_mixing%zeta
      call g%add(nml_real("zeta", pr, "Bottom decay scale (INT_TIDE_DECAY_SCALE)", &
                          units="m"))
      pr => cfg%ocean%tidal_mixing%kd_max
      call g%add(nml_real("kd_max", pr, "Per-layer physical Kd cap (<0 => no cap)", &
                          units="m^2/s"))
      pr => cfg%ocean%tidal_mixing%prandtl_tidal
      call g%add(nml_real("prandtl_tidal", pr, "Kv = prandtl_tidal * Kd"))
      pr => cfg%ocean%tidal_mixing%min_zbot
      call g%add(nml_real("min_zbot", pr, "Mask off where column depth H < min_zbot", &
                          units="m"))
      pr => cfg%ocean%tidal_mixing%e_uniform
      call g%add(nml_real("e_uniform", pr, "Uniform bottom internal-tide energy input E", &
                          units="W m-2"))
      pl => cfg%ocean%tidal_mixing%e_compute
      call g%add(nml_logical("e_compute", pl, &
                             "State-dependent E = min(TKE_coef*N_bot, e_max) (v1.1)"))
      pr => cfg%ocean%tidal_mixing%kappa_itides
      call g%add(nml_real("kappa_itides", pr, "Topographic wavenumber (v1.1 E recompute)", &
                          units="m^-1"))
      pr => cfg%ocean%tidal_mixing%kappa_h2
      call g%add(nml_real("kappa_h2", pr, "KAPPA_H2_FACTOR (v1.1 E recompute)"))
      pr => cfg%ocean%tidal_mixing%utide
      call g%add(nml_real("utide", pr, "RMS barotropic tidal velocity (v1.1 E recompute)", &
                          units="m/s"))
      pr => cfg%ocean%tidal_mixing%h2_rough
      call g%add(nml_real("h2_rough", pr, "Sub-grid topographic roughness variance <h^2>", &
                          units="m^2"))
      pr => cfg%ocean%tidal_mixing%frac_rough
      call g%add(nml_real("frac_rough", pr, "Roughness clamp <h^2> <= (frac_rough*H)^2"))
      pr => cfg%ocean%tidal_mixing%e_max
      call g%add(nml_real("e_max", pr, "TKE_itide_max cap on E", units="W m-2"))
      ps => cfg%ocean%tidal_mixing%tideamp_source
      call g%add(nml_enum("tideamp_source", ps, &
                          "Tidal speed source: uniform (utide) | file (static 2-D field; needs e_compute)", &
                          allowed=[character(len=7) :: "uniform", "file"]))
      ps => cfg%ocean%tidal_mixing%tideamp_file
      call g%add(nml_string("tideamp_file", ps, "NetCDF file with the tidal-amplitude field"))
      ps => cfg%ocean%tidal_mixing%tideamp_var
      call g%add(nml_string("tideamp_var", ps, "Variable name in tideamp_file (m/s)"))
      ps => cfg%ocean%tidal_mixing%h2_source
      call g%add(nml_enum("h2_source", ps, &
                          "Roughness-variance source: uniform (h2_rough) | file (static 2-D field; needs e_compute)", &
                          allowed=[character(len=7) :: "uniform", "file"]))
      ps => cfg%ocean%tidal_mixing%h2_file
      call g%add(nml_string("h2_file", ps, "NetCDF file with the roughness-variance field"))
      ps => cfg%ocean%tidal_mixing%h2_var
      call g%add(nml_string("h2_var", ps, "Variable name in h2_file (m^2)"))
      pl => cfg%ocean%tidal_mixing%n2_bot_rough_avg
      call g%add(nml_logical("n2_bot_rough_avg", pl, &
                             "Near-bed N^2 averaged over the roughness height and imposed in that span (OM3: true)"))
      ps => cfg%ocean%tidal_mixing%energy_limiter
      call g%add(nml_enum("energy_limiter", ps, &
                     "kd_max handling: clamp (clip Kd, discard power) | reinject (cap layer TKE, pass excess up) (OM3: reinject)", &
                          allowed=[character(len=8) :: "clamp", "reinject"]))
      ps => cfg%ocean%tidal_mixing%bed_share
      call g%add(nml_enum("bed_share", ps, &
                          "Bed-layer power share: discard | lump into the first interior layer (OM3: lump)", &
                          allowed=[character(len=7) :: "discard", "lump"]))
      ps => cfg%ocean%tidal_mixing%n2_layer_form
      call g%add(nml_enum("n2_layer_form", ps, &
                          "Layer N^2: legacy (interface value) | mean_interfaces (mean of the two bounding jumps; always under polzin) (OM3: mean_interfaces)", &
                          allowed=[character(len=15) :: "legacy", "mean_interfaces"]))
      ps => cfg%ocean%tidal_mixing%profile
      call g%add(nml_enum("profile", ps, &
                  "Vertical structure: stlaurent (exponential) | polzin (WKB-stretched algebraic; needs e_compute) (OM3: polzin)", &
                          allowed=[character(len=9) :: "stlaurent", "polzin"]))
      pr => cfg%ocean%tidal_mixing%nu_polzin
      call g%add(nml_real("nu_polzin", pr, "Polzin decay-scale constant NU_POLZIN"))
      pr => cfg%ocean%tidal_mixing%nbotref_polzin
      call g%add(nml_real("nbotref_polzin", pr, "Polzin reference bottom N (NBOTREF_POLZIN)", &
                          units="s^-1"))
      pr => cfg%ocean%tidal_mixing%polzin_decay_scale_factor
      call g%add(nml_real("polzin_decay_scale_factor", pr, "Multiplier on the Polzin decay scale"))
      pr => cfg%ocean%tidal_mixing%polzin_scale_max_factor
      call g%add(nml_real("polzin_scale_max_factor", pr, "Polzin stretched scale capped at this multiple of the depth"))
      pr => cfg%ocean%tidal_mixing%polzin_min_decay_scale
      call g%add(nml_real("polzin_min_decay_scale", pr, "Floor on the unstretched Polzin decay scale", &
                          units="m"))

      call schema%add_group(g)
   end subroutine register_tidal_mixing

   module subroutine register_bbl_mixing(cfg, schema)
      !! `&ocean_bbl_mixing` (law-of-the-wall BBL diapycnal mixing).
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      logical, pointer :: pl
      real(wp), pointer :: pr

      g%name = "ocean_bbl_mixing"
      g%doc = "Law-of-the-wall bottom-boundary-layer diapycnal mixing (tracer diffusivity)."

      pl => cfg%ocean%bbl_mixing%enable
      call g%add(nml_logical("enable", pl, &
                             "Master switch (needs the per-face viscous BBL + thermodynamics)"))
      pr => cfg%ocean%bbl_mixing%effic
      call g%add(nml_real("effic", pr, "Fraction of the BBL mean-KE loss spent on mixing (BBL_EFFIC; OM3: 0.01)", &
                          units="nondim"))
      pr => cfg%ocean%bbl_mixing%max_decay
      call g%add(nml_real("max_decay", pr, "Upper bound on the upward power decay length (BBL_MIXING_MAX_DECAY)", &
                          units="m"))
      pl => cfg%ocean%bbl_mixing%use_omega
      call g%add(nml_logical("use_omega", pl, "Floor N^2 at Omega^2 in the energy bound (LOTW_BBL_USE_OMEGA)"))
      pr => cfg%ocean%bbl_mixing%von_karman
      call g%add(nml_real("von_karman", pr, "Von Karman constant (VON_KARMAN_BBL)"))
      pr => cfg%ocean%bbl_mixing%kd_max
      call g%add(nml_real("kd_max", pr, "Cap on the BBL diffusivity and unstratified fill value (KD_MAX)", &
                          units="m^2/s"))

      call schema%add_group(g)
   end subroutine register_bbl_mixing

   module subroutine register_conv(cfg, schema)
      !! `&ocean_conv` (Brunt-Vaisala-triggered convective adjustment,
      !! CVMix_conv-style).
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      logical, pointer :: pl
      real(wp), pointer :: pr

      g%name = "ocean_conv"
      g%doc = "Brunt-Vaisala-triggered convective adjustment (interior closure contributor)."

      pl => cfg%ocean%conv%enable
      call g%add(nml_logical("enable", pl, &
                             "Master switch (requires use_closure + thermodynamics)"))
      pr => cfg%ocean%conv%kd_conv
      call g%add(nml_real("kd_conv", pr, "Convective tracer diffusivity (KD_CONV)", &
                          units="m^2/s"))
      pr => cfg%ocean%conv%prandtl_conv
      call g%add(nml_real("prandtl_conv", pr, "Kv_conv = prandtl_conv * kd_conv"))
      pr => cfg%ocean%conv%n2_thresh
      call g%add(nml_real("n2_thresh", pr, "Trigger threshold on N^2 (BV_SQR_CONV)", &
                          units="s^-2"))

      call schema%add_group(g)
   end subroutine register_conv

   module subroutine register_ddiff(cfg, schema)
      !! `&ocean_ddiff` (double diffusion: salt fingering + diffusive
      !! convection, CVMix_ddiff-style).
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      logical, pointer :: pl
      real(wp), pointer :: pr

      g%name = "ocean_ddiff"
      g%doc = "Double diffusion: salt fingering + diffusive convection (folded into the heat/salt split)."

      pl => cfg%ocean%ddiff%enable
      call g%add(nml_logical("enable", pl, &
                             "Master switch (requires thermodynamics)"))
      pr => cfg%ocean%ddiff%strat_param_max
      call g%add(nml_real("strat_param_max", pr, &
                          "R_rho salt-fingering cutoff (STRAT_PARAM_MAX)"))
      pr => cfg%ocean%ddiff%kappa_ddiff_s
      call g%add(nml_real("kappa_ddiff_s", pr, &
                          "Leading salt-fingering diffusivity K_f (KAPPA_DDIFF_S)", &
                          units="m^2/s"))
      pr => cfg%ocean%ddiff%ddiff_exp1
      call g%add(nml_real("ddiff_exp1", pr, "Inner fingering exponent (DDIFF_EXP1)"))
      pr => cfg%ocean%ddiff%ddiff_exp2
      call g%add(nml_real("ddiff_exp2", pr, "Outer fingering exponent (DDIFF_EXP2)"))
      pr => cfg%ocean%ddiff%param1
      call g%add(nml_real("param1", pr, "MC76 convection exterior coeff (KAPPA_DDIFF_PARAM1)"))
      pr => cfg%ocean%ddiff%param2
      call g%add(nml_real("param2", pr, "MC76 convection middle coeff (KAPPA_DDIFF_PARAM2)"))
      pr => cfg%ocean%ddiff%param3
      call g%add(nml_real("param3", pr, "MC76 convection interior coeff (KAPPA_DDIFF_PARAM3)"))
      pr => cfg%ocean%ddiff%mol_diff
      call g%add(nml_real("mol_diff", pr, &
                          "Molecular diffusivity scaling convection (MOL_DIFF)", &
                          units="m^2/s"))
      pl => cfg%ocean%ddiff%use_k90
      call g%add(nml_logical("use_k90", pl, &
                             "Convection form: MC76 (default) vs Kelley-90"))

      call schema%add_group(g)
   end subroutine register_ddiff

   module subroutine register_epbl(cfg, schema)
      !! `&ocean_epbl` (Reichl & Hallberg 2018 energetics-based PBL).
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      logical, pointer :: pl
      real(wp), pointer :: pr
      integer, pointer :: pi
      character(len=:), pointer :: ps

      g%name = "ocean_epbl"
      g%doc = "RH18 energetics-based planetary boundary layer."

      pl => cfg%ocean%epbl%enable
      call g%add(nml_logical("enable", pl, &
                             "Master switch (replaces the KPP overlay; requires use_closure)"))
      ps => cfg%ocean%epbl%mstar_scheme
      call g%add(nml_enum("mstar_scheme", ps, "Surface TKE mstar scheme", &
                          allowed=[character(len=8) :: "constant", "om4", "rh18"]))
      pr => cfg%ocean%epbl%mstar
      call g%add(nml_real("mstar", pr, "Constant-scheme mstar"))
      pr => cfg%ocean%epbl%mstar_cap
      call g%add(nml_real("mstar_cap", pr, "Cap for OM4/RH18 mstar; off when < 0"))
      pr => cfg%ocean%epbl%mstar_coef1
      call g%add(nml_real("mstar_coef1", pr, "OM4 stabilizing coefficient"))
      pr => cfg%ocean%epbl%c_ek
      call g%add(nml_real("c_ek", pr, "OM4 Ekman coefficient"))
      pr => cfg%ocean%epbl%mstar_conv_adj
      call g%add(nml_real("mstar_conv_adj", pr, "Convective mstar reduction in [0,1]"))
      pr => cfg%ocean%epbl%rh18_cn1
      call g%add(nml_real("rh18_cn1", pr, "RH18 mstar fit coefficient cn1"))
      pr => cfg%ocean%epbl%rh18_cn2
      call g%add(nml_real("rh18_cn2", pr, "RH18 mstar fit coefficient cn2"))
      pr => cfg%ocean%epbl%rh18_cn3
      call g%add(nml_real("rh18_cn3", pr, "RH18 mstar fit coefficient cn3"))
      pr => cfg%ocean%epbl%rh18_cs1
      call g%add(nml_real("rh18_cs1", pr, "RH18 mstar fit coefficient cs1"))
      pr => cfg%ocean%epbl%rh18_cs2
      call g%add(nml_real("rh18_cs2", pr, "RH18 mstar fit coefficient cs2"))
      pr => cfg%ocean%epbl%nstar
      call g%add(nml_real("nstar", pr, "Convective PE -> TKE efficiency"))
      pr => cfg%ocean%epbl%tke_decay
      call g%add(nml_real("tke_decay", pr, "Ekman-depth / TKE-decay-scale ratio"))
      pr => cfg%ocean%epbl%wstar_ustar_coef
      call g%add(nml_real("wstar_ustar_coef", pr, "Convective weight in the velocity scale"))
      ps => cfg%ocean%epbl%vel_scale_scheme
      call g%add(nml_enum("vel_scale_scheme", ps, "Velocity-scale scheme", &
                          allowed=[character(len=9) :: "cube_root", "rh18"]))
      pr => cfg%ocean%epbl%vstar_scale_fac
      call g%add(nml_real("vstar_scale_fac", pr, "Overall vstar multiplier"))
      pr => cfg%ocean%epbl%vstar_surf_fac
      call g%add(nml_real("vstar_surf_fac", pr, "RH18 mechanical surface vstar factor"))
      pr => cfg%ocean%epbl%von_karman
      call g%add(nml_real("von_karman", pr, "von Karman kappa in Kd = vstar*kappa*mixlen"))
      pr => cfg%ocean%epbl%ekman_scale_coef
      call g%add(nml_real("ekman_scale_coef", pr, "Rotational mixing-length rolloff"))
      pr => cfg%ocean%epbl%min_mix_len
      call g%add(nml_real("min_mix_len", pr, "Mixing-length floor", units="m"))
      pr => cfg%ocean%epbl%mixlen_exponent
      call g%add(nml_real("mixlen_exponent", pr, "Shape-function exponent"))
      pr => cfg%ocean%epbl%translay_scale
      call g%add(nml_real("translay_scale", pr, "Transition-layer shape floor (in [0,1) when iterating)"))
      pl => cfg%ocean%epbl%mld_iteration
      call g%add(nml_logical("mld_iteration", pl, "Self-consistent MLD root-find"))
      pr => cfg%ocean%epbl%mld_tol
      call g%add(nml_real("mld_tol", pr, "MLD convergence tolerance", units="m"))
      pi => cfg%ocean%epbl%mld_max_its
      call g%add(nml_int("mld_max_its", pi, "Max MLD iterations", min=1))
      pl => cfg%ocean%epbl%mld_bisection
      call g%add(nml_logical("mld_bisection", pl, "Bisection instead of false position"))
      pl => cfg%ocean%epbl%mld_use_prev_guess
      call g%add(nml_logical("mld_use_prev_guess", pl, "Seed from the previous step's MLD"))
      pr => cfg%ocean%epbl%omega
      call g%add(nml_real("omega", pr, "Earth rotation rate", units="1/s"))
      pr => cfg%ocean%epbl%omega_frac
      call g%add(nml_real("omega_frac", pr, "Blend |f| with 2*Omega"))
      pr => cfg%ocean%epbl%prandtl
      call g%add(nml_real("prandtl", pr, "Kv = prandtl*Kd into the momentum solve"))
      ps => cfg%ocean%epbl%combine
      call g%add(nml_enum("combine", ps, "Combine vs interior closure kv/kt: add | max | max_shear (bg + max(epbl, kappa-shear))", &
                          allowed=[character(len=9) :: "add", "max", "max_shear"]))
      pl => cfg%ocean%epbl%tke_diags
      call g%add(nml_logical("tke_diags", pl, "Compute per-column TKE budget diagnostics"))
      pl => cfg%ocean%epbl%use_lt
      call g%add(nml_logical("use_lt", pl, "Langmuir-turbulence enhancement (LF17 wind-only)"))
      ps => cfg%ocean%epbl%lt_scheme
      call g%add(nml_enum("lt_scheme", ps, "Langmuir enhancement scheme", &
                          allowed=[character(len=8) :: "rescale", "additive"]))
      pr => cfg%ocean%epbl%lt_enhance_coef
      call g%add(nml_real("lt_enhance_coef", pr, "Langmuir enhancement coefficient"))
      pr => cfg%ocean%epbl%lt_enhance_exp
      call g%add(nml_real("lt_enhance_exp", pr, "Langmuir-number exponent"))
      pr => cfg%ocean%epbl%lt_max_enhance
      call g%add(nml_real("lt_max_enhance", pr, "Cap on the multiplicative enhancement"))
      pr => cfg%ocean%epbl%la_frac_hbl
      call g%add(nml_real("la_frac_hbl", pr, "Stokes SL-average depth fraction"))
      pr => cfg%ocean%epbl%lt_lac1
      call g%add(nml_real("lt_lac1", pr, "Stability-modified La coefficient 1"))
      pr => cfg%ocean%epbl%lt_lac2
      call g%add(nml_real("lt_lac2", pr, "Stability-modified La coefficient 2"))
      pr => cfg%ocean%epbl%lt_lac3
      call g%add(nml_real("lt_lac3", pr, "Stability-modified La coefficient 3"))
      pr => cfg%ocean%epbl%lt_lac4
      call g%add(nml_real("lt_lac4", pr, "Stability-modified La coefficient 4"))
      pr => cfg%ocean%epbl%lt_lac5
      call g%add(nml_real("lt_lac5", pr, "Stability-modified La coefficient 5"))

      call schema%add_group(g)
   end subroutine register_epbl

   module subroutine register_wavespeed(cfg, schema)
      !! `&ocean_wavespeed` (B1 first-baroclinic wave speed + Rd).
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      logical, pointer :: pl
      real(wp), pointer :: pr
      integer, pointer :: pi

      g%name = "ocean_wavespeed"
      g%doc = "First-baroclinic wave speed + Rossby deformation radius (diagnostic)."

      pl => cfg%ocean%wavespeed%enable
      call g%add(nml_logical("enable", pl, &
                             "Master switch (diagnostic; default off)"))
      pr => cfg%ocean%wavespeed%mono_n2
      call g%add(nml_real("mono_n2", pr, &
                          "N2-monotonising depth; NOT IMPLEMENTED (refused together with use_ebt); < 0 = off"))
      pl => cfg%ocean%wavespeed%use_ebt
      call g%add(nml_logical("use_ebt", pl, &
   "Also compute the equivalent-barotropic speed cg1_ebt + structure ebt_struct (consumed by GM with ocean_gm_nml use_ebt_struct)"))
      pi => cfg%ocean%wavespeed%n_wavespeed
      call g%add(nml_int("n_wavespeed", pi, "Recompute cadence (every N steps)", min=1))

      call schema%add_group(g)
   end subroutine register_wavespeed

   module subroutine register_ocean_thermo(cfg, schema)
      !! `&ocean_thermo_nml`: thermodynamics switch + scalar surface fluxes.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      logical, pointer :: pl
      real(wp), pointer :: pr
      real(wp), pointer :: pra(:)
      integer, pointer :: pi
      character(len=:), pointer :: ps

      g%name = "ocean_thermo"
      g%doc = "Thermodynamics master switch + scalar surface fluxes."
      ps => cfg%ocean%thermo%sw_source
      call g%add(nml_enum("sw_source", ps, &
                          "Shortwave irradiance source: net_heat (legacy) or q_sw", &
                          allowed=[character(len=8) :: "net_heat", "q_sw"]))
      ps => cfg%ocean%thermo%kpp_sw_method
      call g%add(nml_enum("kpp_sw_method", ps, &
                          "KPP shortwave-in-BL method: all_sw | mxl_sw | lv1_sw", &
                          allowed=[character(len=8) :: "all_sw", "mxl_sw", "lv1_sw"]))
      pl => cfg%ocean%thermo%epbl_sw_ctke
      call g%add(nml_logical("epbl_sw_ctke", pl, &
                             "Charge the EPBL TKE ledger for penetrating shortwave"))
      pl => cfg%ocean%thermo%enable_thermodynamics
      call g%add(nml_logical("enable_thermodynamics", pl, &
                             "Run EOS + tracer advection + vertical mixing"))
      pr => cfg%ocean%thermo%q_heat
      call g%add(nml_real("q_heat", pr, "Net surface heat flux (positive down)", &
                          units="W/m^2"))
      pr => cfg%ocean%thermo%q_salt
      call g%add(nml_real("q_salt", pr, "Net surface salt flux (positive salinifies)", &
                          units="kg/m^2/s"))
      pr => cfg%ocean%thermo%sw_pen_frac
      call g%add(nml_real("sw_pen_frac", pr, &
                          "Penetrating fraction of q_heat (0 = off, all at surface)"))
      pr => cfg%ocean%thermo%sw_band_ratio
      call g%add(nml_real("sw_band_ratio", pr, &
                          "Two-band shortwave band-1 weight R (Jerlov type I)"))
      pr => cfg%ocean%thermo%sw_zeta1
      call g%add(nml_real("sw_zeta1", pr, &
                          "Shortwave band-1 e-folding depth", units="m"))
      pr => cfg%ocean%thermo%sw_zeta2
      call g%add(nml_real("sw_zeta2", pr, &
                          "Shortwave band-2 e-folding depth", units="m"))
      ps => cfg%ocean%thermo%sw_scheme
      call g%add(nml_enum("sw_scheme", ps, &
                          "Penetrating-SW scheme: default (two-band) | manizza05 "// &
                          "(chlorophyll-dependent N-band, Manizza et al. 2005)", &
                          allowed=[character(len=9) :: "default", "manizza05"]))
      pi => cfg%ocean%thermo%sw_nbands
      call g%add(nml_int("sw_nbands", pi, &
                         "manizza05: number of bands (3..6); NIR split evenly over bands 3..n", &
                         min=3, max=MAX_SW_BANDS))
      pr => cfg%ocean%thermo%sw_blue_frac
      call g%add(nml_real("sw_blue_frac", pr, &
                          "manizza05: blue share of the visible (0.42 of total) shortwave", &
                          min=0.0_wp, max=1.0_wp))
      ps => cfg%ocean%thermo%sw_chl_file
      call g%add(nml_string("sw_chl_file", ps, &
                            "manizza05: NetCDF chlorophyll climatology (mg/m^3), monthly cyclic"))
      ps => cfg%ocean%thermo%sw_chl_var
      call g%add(nml_string("sw_chl_var", ps, &
                            "manizza05: variable name of the chlorophyll field"))
      pra => cfg%ocean%thermo%sw_opacity_coef
      call g%add(nml_real_array("sw_opacity_coef", pra, &
                                "manizza05: (a_n, b_n) pairs, kappa_n = a_n + b_n*chl**p_n", &
                                units="1/m"))
      pra => cfg%ocean%thermo%sw_chl_power
      call g%add(nml_real_array("sw_chl_power", pra, &
                                "manizza05: chlorophyll exponents p_n"))
      pr => cfg%ocean%thermo%sw_opacity_land
      call g%add(nml_real("sw_opacity_land", pr, &
                          "manizza05: opacity of land columns", units="1/m", min=0.0_wp))
      call schema%add_group(g)
   end subroutine register_ocean_thermo

   module subroutine register_ocean_wetdry(cfg, schema)
      !! `&ocean_wetdry_nml`: dynamic wetting/drying for the BT substep.
      !! Design + validation numbers: docs/ocean_wetdry_plan.md.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      logical, pointer :: pl
      real(wp), pointer :: pr

      g%name = "ocean_wetdry"
      g%doc = "Dynamic wetting/drying for the split-explicit barotropic substep."
      pl => cfg%ocean%wetdry%enable
      call g%add(nml_logical("enable", pl, &
                             "Dynamic wet/dry: upwind BT face thickness + positive-definite "// &
                             "outflow limiter + bed-blocking momentum gate"))
      pr => cfg%ocean%wetdry%dry_depth
      call g%add(nml_real("dry_depth", pr, &
                          "Total-depth dry threshold", units="m", min=0.0_wp))
      pr => cfg%ocean%wetdry%rewet_depth
      call g%add(nml_real("rewet_depth", pr, &
                          "Hysteresis re-wet threshold (must exceed dry_depth)", &
                          units="m", min=0.0_wp))
      pr => cfg%ocean%wetdry%land_margin
      call g%add(nml_real("land_margin", pr, &
                          "Static-land headroom above rest MSL (intertidal cutoff)", &
                          units="m", min=0.0_wp))
      call schema%add_group(g)
   end subroutine register_ocean_wetdry

   module subroutine register_ocean_pgf(cfg, schema)
      !! `&ocean_pgf_nml`: pressure-gradient-force kernel selector + knobs.
      !! `form` enum mirrors `parse_opgf_variant` in rdb_ocean_pressure_force.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      real(wp), pointer :: pr
      logical, pointer :: pl
      integer, pointer :: pi
      character(len=:), pointer :: ps

      g%name = "ocean_pgf"
      g%doc = "Pressure-gradient-force kernel selector + knobs."
      ps => cfg%ocean%pgf%form
      ! "montgomery" is a long-spelling alias for "mont": `parse_opgf_variant`
      ! accepts it, so the schema must too or the alias is unreachable through
      ! the namelist (the enum gate rejects before the parser is ever called).
      ! Both canonicalise to OPGF_VARIANT_MONT; the only literal comparison on
      ! this key anywhere is against "gprime", so the extra spelling is inert.
      call g%add(nml_enum("form", ps, "PGF kernel variant", &
                          allowed=[character(len=10) :: "mont", "montgomery", "fv_lite", &
                                   "fv_wright", "gprime", "fv_mom6"]))
      pr => cfg%ocean%pgf%gprime_gfs
      call g%add(nml_real("gprime_gfs", pr, "Free-surface gravity for the gprime PGF", &
                          units="m/s^2"))
      pr => cfg%ocean%pgf%gprime_gint
      call g%add(nml_real("gprime_gint", pr, "Internal-interface reduced gravity (gprime PGF)", &
                          units="m/s^2"))
      pr => cfg%ocean%pgf%gfs_scale
      call g%add(nml_real("gfs_scale", pr, "Free-surface gravity scaling (FV_MOM6, MOM6 GFS_scale)"))
      pr => cfg%ocean%pgf%maxvel
      call g%add(nml_real("maxvel", pr, "Velocity-truncation clamp (0 = disabled)", units="m/s"))
      pr => cfg%ocean%pgf%cfl_trunc
      call g%add(nml_real("cfl_trunc", pr, "Advective-CFL velocity truncation threshold (0 = disabled)"))
      pl => cfg%ocean%pgf%mass_weight
      call g%add(nml_logical("mass_weight", pl, &
                             "FV_MOM6 shelf-break hWght mass-weighting at unequal-depth faces"))
      pl => cfg%ocean%pgf%reconstruct_for_pressure
      call g%add(nml_logical("reconstruct_for_pressure", pl, &
                             "FV_MOM6 in-layer PLM/PPM T/S reconstruction for the density integral"))
      pi => cfg%ocean%pgf%recon_scheme
      call g%add(nml_int("recon_scheme", pi, &
                         "In-layer reconstruction scheme: 1=PLM, 2=PPM", min=1, max=2))
      pl => cfg%ocean%pgf%insitu_density
      call g%add(nml_logical("insitu_density", pl, &
                             "FV_MOM6 constant-by-layer density at its in-situ pressure "// &
                             "(MOM6 PCM parity); .false. = legacy potential density at "// &
                             "&ocean_eos_nml p_ref. Inert for the linear EOS"))
      pl => cfg%ocean%pgf%p_top_in_bc
      call g%add(nml_logical("p_top_in_bc", pl, &
                             "FV_MOM6: add the top-of-column load ms%p_top to "// &
                             "the pressure-stack surface boundary condition "// &
                             "pa(nz+1) = rho_ref*g*eta + p_top"))
      call schema%add_group(g)
   end subroutine register_ocean_pgf

   module subroutine register_ocean_eos(cfg, schema)
      !! `&ocean_eos_nml`: equation-of-state variant selector, the
      !! freezing-point (liquidus) coefficient set, and the
      !! potential-density reference pressure.
      !! `eos` enum mirrors `parse_eos_variant` in rdb_eos;
      !! `tfreeze_set` mirrors `parse_tfreeze_set` in the same module.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      character(len=:), pointer :: ps
      real(wp), pointer :: pr

      g%name = "ocean_eos"
      g%doc = "Equation-of-state variant selector, liquidus set + reference pressure."
      ps => cfg%ocean%eos%eos
      call g%add(nml_enum("eos", ps, "Equation-of-state variant", &
                          allowed=[character(len=10) :: "linear", "wright", &
                                   "roquet_spv", "roquet_rho", "teos10"]))
      ps => cfg%ocean%eos%tfreeze_set
      call g%add(nml_enum("tfreeze_set", ps, &
                          "Named liquidus coefficient set for eos_freezing_point "// &
                          "(T_f = l1*S + l2 + l3*p): 'seaice' = SIS2/MOM6 "// &
                          "(-0.054, 0, -7.53e-8), 'isomip' = ISOMIP+ "// &
                          "(-0.0573, 0.0832, -7.53e-8)", &
                          allowed=[character(len=6) :: "seaice", "isomip"]))
      ps => cfg%ocean%eos%ts_convention
      call g%add(nml_enum("ts_convention", ps, &
                          "Meaning of the T/S tracers: 'pot_prac' = potential "// &
                          "temperature + practical salinity (default); 'cons_abs' = "// &
                          "TEOS-10 Conservative Temperature + Absolute Salinity "// &
                          "(roquet_rho/roquet_spv only)", &
                          allowed=[character(len=8) :: "pot_prac", "cons_abs"]))
      ps => cfg%ocean%eos%tfreeze_form
      call g%add(nml_enum("tfreeze_form", ps, &
                          "Liquidus form: 'linear' (tfreeze_set coefficients, default) "// &
                          "or 'teos_poly' (TEOS-10 CT_freezing_poly(SA,p); needs "// &
                          "ts_convention='cons_abs')", &
                          allowed=[character(len=9) :: "linear", "teos_poly"]))
      pr => cfg%ocean%eos%tfreeze_sat_frac
      call g%add(nml_real("tfreeze_sat_frac", pr, &
                          "Air-saturation fraction of the TEOS-10 polynomial liquidus "// &
                          "(0 = no dissolved air; ACCESS-OM3 parity)", &
                          min=0.0_wp, max=1.0_wp))
      pr => cfg%ocean%eos%c_p
      call g%add(nml_real("c_p", pr, &
                          "Seawater heat capacity for ocean heat-flux <-> temperature "// &
                          "conversion (MOM6 C_P); 3992 (OM3) by default under both "// &
                          "ts_conventions, TEOS-10 cp0 = 3991.86795711963", &
                          units="J kg-1 K-1", min=3000.0_wp, max=4500.0_wp))
      pr => cfg%ocean%eos%p_ref
      call g%add(nml_real("p_ref", pr, &
                          "Reference pressure for the potential density "// &
                          "ms%rho_layer (horizontally uniform by design; "// &
                          "0 => surface density)", &
                          units="Pa", min=0.0_wp))
      call schema%add_group(g)
   end subroutine register_ocean_eos

   module subroutine register_ocean_bdrag(cfg, schema)
      !! `&ocean_bdrag_nml`: bottom-drag selector + coefficients.
      !! `form` enum mirrors `parse_bdrag_variant` in rdb_ocean_bottom_drag.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      real(wp), pointer :: pr
      character(len=:), pointer :: ps
      logical, pointer :: pl

      g%name = "ocean_bdrag"
      g%doc = "Bottom-drag selector + coefficients."
      ps => cfg%ocean%bdrag%form
      call g%add(nml_enum("form", ps, "Bottom-drag variant", &
                          allowed=[character(len=9) :: "quadratic", "linear"]))
      pr => cfg%ocean%bdrag%cd
      call g%add(nml_real("cd", pr, "Quadratic drag coefficient (0 disables)"))
      pr => cfg%ocean%bdrag%r
      call g%add(nml_real("r", pr, "Linear Rayleigh coefficient (0 disables)", units="1/s"))
      pr => cfg%ocean%bdrag%hbbl
      call g%add(nml_real("hbbl", pr, "BBL thickness for distributed drag (0 = bed-only)", &
                          units="m"))
      pr => cfg%ocean%bdrag%bg_vel
      call g%add(nml_real("bg_vel", pr, "Background velocity floor for distributed drag", &
                          units="m/s"))
      pr => cfg%ocean%bdrag%bbl_thick_min
      call g%add(nml_real("bbl_thick_min", pr, "Minimum effective BBL thickness", units="m"))
      pr => cfg%ocean%bdrag%bed_factor
      call g%add(nml_real("bed_factor", pr, "Multiplier on the bed-layer drag tendency only"))
      pl => cfg%ocean%bdrag%channel_drag
      call g%add(nml_logical("channel_drag", pl, &
                             "Enable per-layer lateral side-wall (channel) Rayleigh drag"))
      pr => cfg%ocean%bdrag%cdrag_side
      call g%add(nml_real("cdrag_side", pr, &
                          "Side-wall drag coefficient (0 disables channel drag)"))
      pl => cfg%ocean%bdrag%implicit
      call g%add(nml_logical("implicit", pl, &
                             "Backward-Euler implicit drag (stable for thin bottom layers)"))
      call schema%add_group(g)
   end subroutine register_ocean_bdrag

   module subroutine register_ocean_tdrag(cfg, schema)
      !! `&ocean_tdrag_nml`: ice-shelf TOP-drag selector + coefficients.
      !! `form` enum mirrors `parse_tdrag_variant` in rdb_ocean_top_drag.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      real(wp), pointer :: pr
      character(len=:), pointer :: ps
      logical, pointer :: pl

      g%name = "ocean_tdrag"
      g%doc = "Ice-shelf top-drag selector + coefficients (mirror of &ocean_bdrag_nml)."
      pl => cfg%ocean%tdrag%enable
      call g%add(nml_logical("enable", pl, &
                             "Enable the ice-shelf top drag (requires "// &
                             "&ocean_cavity_dyn_nml enable)"))
      ps => cfg%ocean%tdrag%form
      call g%add(nml_enum("form", ps, "Top-drag variant", &
                          allowed=[character(len=9) :: "quadratic", "linear"]))
      pr => cfg%ocean%tdrag%cd
      call g%add(nml_real("cd", pr, &
                          "Quadratic top-drag coefficient (0 disables); must equal "// &
                          "&ocean_cavity_melt_nml cdrag_top when melt is on", &
                          min=0.0_wp))
      pr => cfg%ocean%tdrag%r
      call g%add(nml_real("r", pr, "Linear Rayleigh top-drag coefficient (0 disables)", &
                          units="1/s", min=0.0_wp))
      pr => cfg%ocean%tdrag%htbl
      call g%add(nml_real("htbl", pr, &
                          "Top-boundary-layer thickness for distributed drag "// &
                          "(0 = layer-nz only)", units="m", min=0.0_wp))
      pr => cfg%ocean%tdrag%bg_vel
      call g%add(nml_real("bg_vel", pr, &
                          "Background velocity floor in the quadratic top-drag speed", &
                          units="m/s", min=0.0_wp))
      pr => cfg%ocean%tdrag%tbl_thick_min
      call g%add(nml_real("tbl_thick_min", pr, &
                          "Minimum effective TBL thickness (0 = fall back to h_min)", &
                          units="m", min=0.0_wp))
      pl => cfg%ocean%tdrag%implicit
      call g%add(nml_logical("implicit", pl, &
                             "Backward-Euler top drag inside the drag kernel "// &
                             "(stable for thin top layers)"))
      call schema%add_group(g)
   end subroutine register_ocean_tdrag

   module subroutine register_ocean_vmix(cfg, schema)
      !! `&ocean_vmix_nml`: vertical-mixing module switches + knobs.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      logical, pointer :: pl
      real(wp), pointer :: pr
      integer, pointer :: pi
      character(len=:), pointer :: ps

      g%name = "ocean_vmix"
      g%doc = "Vertical-mixing module switches + knobs."
      pl => cfg%ocean%vmix%use_closure
call g%add(nml_logical("use_closure", pl, "Master switch for the interior closure (PP81 or interior_closure) and its contributors"))
      ps => cfg%ocean%vmix%interior_closure
      call g%add(nml_enum("interior_closure", ps, &
                          "Interior closure under use_closure: pp81 (default) | none (background + other schemes only)", &
                          allowed=[character(len=8) :: "pp81", "none"]))
      pl => cfg%ocean%vmix%use_kpp
      call g%add(nml_logical("use_kpp", pl, "KPP surface-boundary-layer overlay (needs use_closure)"))
      pl => cfg%ocean%vmix%direct_stress
      call g%add(nml_logical("direct_stress", pl, "Spread wind stress over hmix_stress (MOM6 DIRECT_STRESS)"))
      pr => cfg%ocean%vmix%hmix_stress
      call g%add(nml_real("hmix_stress", pr, "Surface-slab thickness for direct_stress", units="m"))
      pr => cfg%ocean%vmix%kv_ml_invz2
      call g%add(nml_real("kv_ml_invz2", pr, "Near-surface 1/(z hmix)^2 viscosity (MOM6 KV_ML_INVZ2)", &
                          units="m^2/s"))
      pr => cfg%ocean%vmix%hmix_fixed
      call g%add(nml_real("hmix_fixed", pr, "Mixed-layer thickness for the KV_ML_INVZ2 profile", units="m"))
      pl => cfg%ocean%vmix%harmonic_visc
      call g%add(nml_logical("harmonic_visc", pl, "Harmonic-mean face thickness in vdiff (MOM6 HARMONIC_VISC)"))
      pi => cfg%ocean%vmix%dt_therm_ratio
      call g%add(nml_int("dt_therm_ratio", pi, "Thermo step runs at ratio*dt_dyn (MOM6 DT_THERM)", min=1))
      pl => cfg%ocean%vmix%diabatic_first
      call g%add(nml_logical("diabatic_first", pl, &
                             "Thermo step + ALE remap lead (T) or close (F, OM3) each dt_therm interval (MOM6 DIABATIC_FIRST)"))
      pi => cfg%ocean%vmix%dt_tracer_advect_ratio
      call g%add(nml_int("dt_tracer_advect_ratio", pi, &
                         "Horizontal tracer advect runs at ratio*dt_dyn (MOM6 DT_TRACER_ADVECT)", min=1))
      ps => cfg%ocean%vmix%tracer_advection_scheme
      call g%add(nml_enum("tracer_advection_scheme", ps, &
                          "Horizontal tracer-advection scheme (ocean, every-step and windowed drain): ppm|ppm_h3", &
                          allowed=[character(len=8) :: "ppm", "ppm_h3"]))
      ps => cfg%ocean%vmix%tracer_recon
      call g%add(nml_enum("tracer_recon", ps, &
                          "Windowed tracer-advect drain face reconstruction (ocean): ppm|weno5|weno7|weno9", &
                          allowed=[character(len=8) :: "ppm", "weno5", "weno7", "weno9"]))
      pr => cfg%ocean%vmix%kv_max
      call g%add(nml_real("kv_max", pr, "Assembly ceiling on kv; huge=off (MOM6 Kd_max momentum)", &
                          units="m^2/s"))
      pr => cfg%ocean%vmix%kd_max
      call g%add(nml_real("kd_max", pr, "Assembly ceiling on kt/ks; huge=off (MOM6 Kd_max)", &
                          units="m^2/s"))
      pi => cfg%ocean%vmix%kd_smooth_iterations
      call g%add(nml_int("kd_smooth_iterations", pi, "1-2-1 smoothing passes on kv/kt (MOM6 Kd_smooth)", min=0))
      pl => cfg%ocean%vmix%vmix_guard
      call g%add(nml_logical("vmix_guard", pl, "Debug-gated negative/NaN diffusivity guard (assembly)"))
      ! C7 Bryan-Lewis + Henyey depth-varying background.
      pl => cfg%ocean%vmix%bkgnd_profile
      call g%add(nml_logical("bkgnd_profile", pl, "Bryan-Lewis depth-varying background diffusivity"))
      pr => cfg%ocean%vmix%bkgnd_kd_sfc
      call g%add(nml_real("bkgnd_kd_sfc", pr, "Bryan-Lewis surface-asymptote background Kd", units="m^2/s"))
      pr => cfg%ocean%vmix%bkgnd_kd_deep
      call g%add(nml_real("bkgnd_kd_deep", pr, "Bryan-Lewis deep-asymptote background Kd", units="m^2/s"))
      pr => cfg%ocean%vmix%bkgnd_z0
      call g%add(nml_real("bkgnd_z0", pr, "Bryan-Lewis transition-centre depth", units="m"))
      pr => cfg%ocean%vmix%bkgnd_delta
      call g%add(nml_real("bkgnd_delta", pr, "Bryan-Lewis transition half-width", &
                          units="m", min=1.0e-6_wp))
      pr => cfg%ocean%vmix%bkgnd_prandtl
call g%add(nml_real("bkgnd_prandtl", pr, "Background Prandtl number: kv floor = prandtl*Kd_bg + kv_bg; -1 (exact) off, 0 refused", &
                          min=-1.0_wp))
      pl => cfg%ocean%vmix%bkgnd_henyey
      call g%add(nml_logical("bkgnd_henyey", pl, &
                             "Henyey IGW latitude factor on the scalar background "// &
                             "(excludes bkgnd_profile; needs a non-cartesian grid)"))
      pr => cfg%ocean%vmix%bkgnd_kd_min
      call g%add(nml_real("bkgnd_kd_min", pr, &
                          "Minimum background Kd under the Henyey scaling "// &
                          "(negative = 0.01*kt_bg)", units="m^2/s"))
      pr => cfg%ocean%vmix%bkgnd_henyey_n0_2omega
      call g%add(nml_real("bkgnd_henyey_n0_2omega", pr, &
                          "Henyey N0/(2*Omega) reference stratification ratio", &
                          min=1.0_wp))
      pr => cfg%ocean%vmix%bkgnd_henyey_max_lat
      call g%add(nml_real("bkgnd_henyey_max_lat", pr, &
                          "Latitude poleward of which the Henyey factor floors", &
                          units="degN", min=0.0_wp))
      ! PP81 interior closure + KPP BL-depth constants — route to the
      ! OCEAN vmix slot (see rdb_ocean_setup.F90:configure_ocean_lateral).
      ! Do NOT confuse with the similarly-named &nonhydrostatic_nml
      ! kpp_ri_crit/kpp_cs_nonlocal/kpp_c_vt2, which drove the coastal KPP
      ! and are VESTIGIAL since that path was carved out — nothing reads
      ! them.  These &ocean_vmix_nml keys are the live ones.
      pr => cfg%ocean%vmix%pp81_nu0
      call g%add(nml_real("pp81_nu0", pr, "PP81 Richardson-dependent viscosity scale", &
                          units="m^2/s", min=0.0_wp))
      pr => cfg%ocean%vmix%pp81_nu_bg
      call g%add(nml_real("pp81_nu_bg", pr, "PP81 background viscosity (also seeds vmix%kv_bg)", &
                          units="m^2/s", min=0.0_wp))
      pr => cfg%ocean%vmix%pp81_kappa_bg
      call g%add(nml_real("pp81_kappa_bg", pr, &
                          "PP81 background diffusivity (also seeds vmix%kt_bg/ks_bg)", &
                          units="m^2/s", min=0.0_wp))
      pr => cfg%ocean%vmix%pp81_alpha
      call g%add(nml_real("pp81_alpha", pr, &
                          "PP81 Richardson-number scaling coefficient (paper value 5)", &
                          min=1.0e-12_wp))
      pr => cfg%ocean%vmix%shear2_floor
      call g%add(nml_real("shear2_floor", pr, &
                          "Floor on |du/dz|^2+|dv/dz|^2 in the PP81 Ri denominator", &
                          units="1/s^2", min=1.0e-30_wp))
      pr => cfg%ocean%vmix%kpp_ri_crit
      call g%add(nml_real("kpp_ri_crit", pr, "Critical bulk Richardson number for KPP BL-depth", &
                          min=1.0e-12_wp))
      pr => cfg%ocean%vmix%kpp_cs_nonlocal
      call g%add(nml_real("kpp_cs_nonlocal", pr, &
                          "KPP non-local (counter-gradient) transport coefficient C_s", &
                          min=0.0_wp))
      pr => cfg%ocean%vmix%kpp_c_vt2
      call g%add(nml_real("kpp_c_vt2", pr, &
                          "KPP unresolved-turbulence V_t^2 coefficient (0 disables V_t^2)", &
                          min=0.0_wp))
      ! E4: source of the alpha/beta pair the KPP B_0 and the
      ! double-diffusion density ratio use.  The `allowed=` list must stay
      ! in lockstep with `parse_buoyancy_coeffs` (rdb_ocean_vmix.F90).
      ps => cfg%ocean%vmix%buoyancy_coeffs
      call g%add(nml_enum("buoyancy_coeffs", ps, &
                          "Source of alpha/beta for KPP B_0 + double diffusion: "// &
                          "constant (scalar &ocean_ic_nml pair) | eos (active EOS derivatives)", &
                          allowed=[character(len=8) :: "constant", "eos"]))
      call schema%add_group(g)
   end subroutine register_ocean_vmix

   module subroutine register_ocean_vdiff(cfg, schema)
      !! `&ocean_vdiff_nml`: backward-Euler vertical-friction folding knobs.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      logical, pointer :: pl
      real(wp), pointer :: pr

      g%name = "ocean_vdiff"
      g%doc = "Backward-Euler vertical-friction solver knobs."
      pl => cfg%ocean%vdiff%implicit_stress
      call g%add(nml_logical("implicit_stress", pl, &
                             "Fold wind stress into the vdiff surface (k=nz) RHS"))
      pl => cfg%ocean%vdiff%accel_visc_rem
      call g%add(nml_logical("accel_visc_rem", pl, &
                             "RETIRED (refused when set)", &
                             dead_on_ocean_path="RETIRED -- PR-3's audit found no MOM6 "// &
                             "state-update equivalent: btstep_layer_accel applies the "// &
                             "depth-mean barotropic acceleration uniformly across every "// &
                             "layer, no visc_rem weight. Setting it .true. is a "// &
                             "fail-loud configure error (validate_config). The real MOM6 "// &
                             "visc_rem*velocity mechanisms are ocean_bt_nml "// &
                             "renorm_visc_rem and rescale_strong_drag."))
      pl => cfg%ocean%vdiff%implicit_drag
      call g%add(nml_logical("implicit_drag", pl, &
                             "Fold bottom drag into the vdiff bed (k=1) diagonal"))
      pl => cfg%ocean%vdiff%implicit_top_drag
      call g%add(nml_logical("implicit_top_drag", pl, &
                             "Fold the ice-shelf top drag into the vdiff surface "// &
                             "(k=nz) diagonal, masking the wind RHS under cover"))
      pl => cfg%ocean%vdiff%hvel_mom6
      call g%add(nml_logical("hvel_mom6", pl, &
                             "MOM6 HARMONIC_VISC parity: harmonic momentum face "// &
                             "thickness with the near-bed upwind blend, and arithmetic "// &
                             "h_shear. Suppresses grounded-sliver momentum as MOM6 does"))
      pl => cfg%ocean%vdiff%hvel_harmonic
      call g%add(nml_logical("hvel_harmonic", pl, &
                             "MOM6 HARMONIC_VISC for the hvel_mom6 face thickness: .false. = "// &
                             "MOM6 default (arithmetic + z_clear near-bed harmonic blend), "// &
                             ".true. = harmonic + near-bed upwind-arithmetic blend"))
      pr => cfg%ocean%vdiff%hbbl_visc
      call g%add(nml_real("hbbl_visc", pr, &
                          "Bottom-layer scale for the hvel_mom6 botfn blend without the BBL "// &
                          "glue, and the glue's HBBL when ocean_bdrag hbbl = 0 (MOM6 HBBL)", &
                          units="m"))
      pl => cfg%ocean%vdiff%bbl_glue
      call g%add(nml_logical("bbl_glue", pl, &
                             "MOM6 BOTTOMDRAGLAW: per-face set_viscous_BBL kv_bbl/bbl_thick "// &
                             "from the ocean_bdrag law (quadratic or linear); kv_bbl botfn "// &
                             "glue at near-bed interfaces + piston bed row replace the bed "// &
                             "drag apply. Requires hvel_mom6"))
      pr => cfg%ocean%vdiff%bbl_piston
      call g%add(nml_real("bbl_piston", pr, &
                          "Historical constant BBL piston u* (hand-built slots only; a "// &
                          "configured glue takes its drag from ocean_bdrag)", &
                          units="m/s", min=0.0_wp))
      pl => cfg%ocean%vdiff%hvel_upwind
      call g%add(nml_logical("hvel_upwind", pl, &
                             "Near-bed upwind blend in the hvel_mom6 thickness build "// &
                             "(.false. = pure harmonic; the u-sign blend flip-flops on "// &
                             "roundoff at rest and collapses the BBL glue, PGF_BUG.md par.9.8)"))
      call schema%add_group(g)
   end subroutine register_ocean_vdiff

   module subroutine register_ocean_continuity(cfg, schema)
      !! `&ocean_continuity_nml`: continuity-PPM positivity controls.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      real(wp), pointer :: pr
      logical, pointer :: pl

      g%name = "ocean_continuity"
      g%doc = "Continuity-PPM positivity controls."
      pr => cfg%ocean%continuity%h_min
      call g%add(nml_real("h_min", pr, "PPM positivity-limiter thickness floor", units="m"))
      pl => cfg%ocean%continuity%ppm_limit_pos
      call g%add(nml_logical("ppm_limit_pos", pl, "Positivity-preserving PPM face limiter (MOM6 PPM_limit_pos)"))
      pl => cfg%ocean%continuity%vol_cfl
      call g%add(nml_logical("vol_cfl", pl, "Swept-volume continuity-PPM face thickness (MOM6 vol_CFL)"))
      pl => cfg%ocean%continuity%positive_definite
      call g%add(nml_logical("positive_definite", pl, &
                             "Positive-definite split continuity (h>=h_lim, zero mass created; "// &
                             "reconstruction floor + per-donor outflux limiter)"))
      pl => cfg%ocean%continuity%renorm_consistent_flux
      call g%add(nml_logical("renorm_consistent_flux", pl, &
                             "uhbt renormalisation: continuous flux for a layer whose "// &
                             "upwind donor flips under the correction (MOM6 "// &
                             "zonal_flux_adjust); fixes the wrong-sign transport at "// &
                             "thickness jumps"))
      call schema%add_group(g)
   end subroutine register_ocean_continuity

end submodule rdb_config_register_vmix
