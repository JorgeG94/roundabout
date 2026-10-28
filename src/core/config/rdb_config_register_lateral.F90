submodule(rdb_config) rdb_config_register_lateral
   !! Strict-schema registration: lateral-closure + dynamics groups.
   implicit none

contains

   module subroutine register_slopes(cfg, schema)
      !! `&ocean_slopes` (Griffies 1998 isopycnal-slope diagnostics).
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      logical, pointer :: pl
      real(wp), pointer :: pr

      g%name = "ocean_slopes"
      g%doc = "Isopycnal (neutral) slope diagnostics (Griffies 1998)."

      pl => cfg%ocean%slopes%enable
      call g%add(nml_logical("enable", pl, &
                             "Master switch (diagnostic; default off ⇒ no-op)"))
      pr => cfg%ocean%slopes%kd_smooth
      call g%add(nml_real("kd_smooth", pr, &
                          "Vert-fill smoothing diffusivity (× dt fills massless layers)", &
                          units="m^2/s"))
      pr => cfg%ocean%slopes%min_dz_for_n2
      call g%add(nml_real("min_dz_for_n2", pr, &
                          "Minimum layer thickness floored in the N²/drdz denominator", &
                          units="m"))

      call schema%add_group(g)
   end subroutine register_slopes

   module subroutine register_gm(cfg, schema)
      !! `&ocean_gm` (Gent-McWilliams thickness diffusion, capability [2]).
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      logical, pointer :: pl
      real(wp), pointer :: pr
      character(len=:), pointer :: ps

      g%name = "ocean_gm"
      g%doc = "Gent-McWilliams thickness diffusion (eddy bolus transport)."

      pl => cfg%ocean%gm%enable
      call g%add(nml_logical("enable", pl, &
                             "Master switch (requires &ocean_slopes_nml enable; default off ⇒ no-op)"))
      pr => cfg%ocean%gm%khth
      call g%add(nml_real("khth", pr, &
                          "Thickness diffusivity KhTh (constant-fill; production 1e2-1e3)", &
                          units="m^2/s"))
      pr => cfg%ocean%gm%khth_max_cfl
      call g%add(nml_real("khth_max_cfl", pr, &
                          "Fraction of the diffusive CFL the face KH may use", units="nondim"))
      pr => cfg%ocean%gm%khth_slope_max
      call g%add(nml_real("khth_slope_max", pr, &
                          "Slope magnitude above which the safe-streamfunction blend takes over", &
                          units="nondim"))
      pl => cfg%ocean%gm%use_ebt_struct
      call g%add(nml_logical("use_ebt_struct", pl, &
                             "EBT vertical structure of KhTh (KHTH_USE_EBT_STRUCT; needs wavespeed use_ebt; default off)"))
      pl => cfg%ocean%gm%use_fgnv
      call g%add(nml_logical("use_fgnv", pl, &
                    "Solve the FGNV (2010) vertical BVP for the GM streamfunction (KHTH_USE_FGNV_STREAMFUNCTION; needs wavespeed)"))
      pr => cfg%ocean%gm%fgnv_filter_scale
      call g%add(nml_real("fgnv_filter_scale", pr, &
                          "FGNV smoothing-term scale: c^2 = scale*c^2, psi_unlim amplified by (1+scale) (FGNV_FILTER_SCALE, > 0)", &
                          units="nondim"))
      pr => cfg%ocean%gm%fgnv_c_min
      call g%add(nml_real("fgnv_c_min", pr, &
                          "FGNV_C_MIN, accepted for parity but not applied (>= 0)", units="m/s"))
      pr => cfg%ocean%gm%fgnv_strat_floor
      call g%add(nml_real("fgnv_strat_floor", pr, &
                          "FGNV N^2 floor as a fraction of Omega: N^2 >= (floor*Omega)^2 (FGNV_STRAT_FLOOR, > 0)", &
                          units="nondim"))
      ps => cfg%ocean%gm%gm_src_form
      call g%add(nml_enum("gm_src_form", ps, &
                          "MEKE's GM energy source: 'local' = rho0/4 KH S^2 N^2 h (default); "// &
                          "'work' = PE release by the final bolus transports (MEKE_GM_SRC_ALT=False)", &
                          allowed=[character(len=15) :: "local", "work"]))

      call schema%add_group(g)
   end subroutine register_gm

   module subroutine register_redi(cfg, schema)
      !! `&ocean_redi` (continuous neutral / along-isopycnal tracer
      !! diffusion, capability [3]).
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      logical, pointer :: pl
      real(wp), pointer :: pr

      g%name = "ocean_redi"
      g%doc = "Redi continuous neutral (along-isopycnal) tracer diffusion."

      pl => cfg%ocean%redi%enable
      call g%add(nml_logical("enable", pl, &
                             "Master switch (default off ⇒ no-op; augments hdiff_tracer)"))
      pl => cfg%ocean%redi%continuous
      call g%add(nml_logical("continuous", pl, &
                             "Continuous variant (discontinuous deferred R4; .false. rejected)"))
      pr => cfg%ocean%redi%khtr
      call g%add(nml_real("khtr", pr, &
                          "Redi neutral diffusivity KhTr (production 1e2-1e3)", units="m^2/s"))
      pl => cfg%ocean%redi%check_diffusive_cfl
      call g%add(nml_logical("check_diffusive_cfl", pl, &
                             "Subcycle Redi so the explicit-diffusion CFL stays <= 1 per pass (default off)"))
      pr => cfg%ocean%redi%max_diffusive_cfl
      call g%add(nml_real("max_diffusive_cfl", pr, &
                          "Redi: cap face coefficients at this cell CFL (<=0 off; OM3 2.0)", units="nondim"))

      call schema%add_group(g)
   end subroutine register_redi

   module subroutine register_varmix(cfg, schema)
      !! `&ocean_varmix` (spatially-varying GM/Redi coefficients, [4]).
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      logical, pointer :: pl
      real(wp), pointer :: pr
      integer, pointer :: pi

      g%name = "ocean_varmix"
      g%doc = "Spatially-varying GM/Redi lateral-diffusivity coefficients."

      pl => cfg%ocean%varmix%enable
      call g%add(nml_logical("enable", pl, &
                             "Master switch (requires slopes + wavespeed; default off ⇒ GM uses const khth)"))
      pl => cfg%ocean%varmix%use_visbeck
      call g%add(nml_logical("use_visbeck", pl, &
                             "Add the Visbeck/Eady khth_slope_cff·L²·SN baroclinicity term"))
      pl => cfg%ocean%varmix%resoln_scaled_khth
      call g%add(nml_logical("resoln_scaled_khth", pl, &
                             "Scale the assembled KhTh by the resolution function"))
      pl => cfg%ocean%varmix%resoln_scaled_khtr
      call g%add(nml_logical("resoln_scaled_khtr", pl, &
                             "Scale the assembled KhTr by the resolution function"))
      pl => cfg%ocean%varmix%gill_equatorial_ld
      call g%add(nml_logical("gill_equatorial_ld", pl, &
                             "Gill (1982) equatorial-Ld convention (factor 2 in beta_dx2)"))
      pl => cfg%ocean%varmix%interpolate_res_fn
      call g%add(nml_logical("interpolate_res_fn", pl, &
                             "Interpolate centre Res_fn to faces (else interpolate cg1, MOM6 default)"))
      pi => cfg%ocean%varmix%kh_res_fn_power
      call g%add(nml_int("kh_res_fn_power", pi, "Resolution-function power p (integer >= 1; 2 is the production form)"))
      pr => cfg%ocean%varmix%kh_res_scale_coef
      call g%add(nml_real("kh_res_scale_coef", pr, &
                          "Resolution-function alpha ((alpha·cg1)^p denominator coef)", units="nondim"))
      pr => cfg%ocean%varmix%khth
      call g%add(nml_real("khth", pr, "Background thickness diffusivity KhTh", units="m^2/s"))
      pr => cfg%ocean%varmix%khtr
      call g%add(nml_real("khtr", pr, "Background tracer diffusivity KhTr (future Redi)", units="m^2/s"))
      pr => cfg%ocean%varmix%khth_slope_cff
      call g%add(nml_real("khth_slope_cff", pr, "Visbeck coefficient for the KhTh chain", units="nondim"))
      pr => cfg%ocean%varmix%khtr_slope_cff
      call g%add(nml_real("khtr_slope_cff", pr, "Visbeck coefficient for the KhTr chain", units="nondim"))
      pr => cfg%ocean%varmix%khth_min
      call g%add(nml_real("khth_min", pr, "Lower clamp on KhTh", units="m^2/s"))
      pr => cfg%ocean%varmix%khth_max
      call g%add(nml_real("khth_max", pr, "Upper clamp on KhTh (<= 0 ⇒ no cap)", units="m^2/s"))
      pr => cfg%ocean%varmix%khtr_min
      call g%add(nml_real("khtr_min", pr, "Lower clamp on KhTr", units="m^2/s"))
      pr => cfg%ocean%varmix%khtr_max
      call g%add(nml_real("khtr_max", pr, "Upper clamp on KhTr (<= 0 ⇒ no cap)", units="m^2/s"))
      pr => cfg%ocean%varmix%visbeck_l_scale
      call g%add(nml_real("visbeck_l_scale", pr, &
                          "Visbeck length scale L (m); if < 0, |L|²·areaCu", units="m"))
      pr => cfg%ocean%varmix%visbeck_max_slope
      call g%add(nml_real("visbeck_max_slope", pr, &
                          "S² limiter scale (<= 0 ⇒ no limit)", units="nondim"))

      call schema%add_group(g)
   end subroutine register_varmix

   module subroutine register_meke(cfg, schema)
      !! `&ocean_meke` (prognostic mesoscale eddy kinetic energy, [5]).
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      logical, pointer :: pl
      real(wp), pointer :: pr

      g%name = "ocean_meke"
      g%doc = "Prognostic mesoscale eddy kinetic energy (GM<->eddy loop)."

      pl => cfg%ocean%meke%enable
      call g%add(nml_logical("enable", pl, &
                             "Master switch (requires &ocean_gm_nml enable; default off ⇒ no-op)"))
      pr => cfg%ocean%meke%gmcoeff
      call g%add(nml_real("gmcoeff", pr, "PE->MEKE conversion efficiency (< 0 ⇒ off)", units="nondim"))
      pr => cfg%ocean%meke%frcoeff
      call g%add(nml_real("frcoeff", pr, "Frictional mean->eddy conversion (< 0 ⇒ off; >=0 sources hvisc KE dissipation)", &
                          units="nondim"))
      pr => cfg%ocean%meke%bgsrc
      call g%add(nml_real("bgsrc", pr, "Background energy source", units="m^2/s^3"))
      pr => cfg%ocean%meke%damping
      call g%add(nml_real("damping", pr, "Linear MEKE dissipation rate", units="1/s"))
      pr => cfg%ocean%meke%kh
      call g%add(nml_real("kh", pr, "Background lateral diffusion of MEKE (< 0 ⇒ off)", units="m^2/s"))
      pr => cfg%ocean%meke%k4
      call g%add(nml_real("k4", pr, "Background biharmonic diffusion of MEKE (< 0 ⇒ off)", units="m^4/s"))
      pr => cfg%ocean%meke%khcoeff
      call g%add(nml_real("khcoeff", pr, "MEKE->Kh scaling (<= 0 ⇒ closure off)", units="nondim"))
      pr => cfg%ocean%meke%cd_scale
      call g%add(nml_real("cd_scale", pr, "Bottom/column eddy-velocity ratio", units="nondim"))
      pr => cfg%ocean%meke%cb
      call g%add(nml_real("cb", pr, "Coefficient in gamma_bot (bottomFac2)", units="nondim"))
      pr => cfg%ocean%meke%ct
      call g%add(nml_real("ct", pr, "Coefficient in gamma_bt (barotrFac2)", units="nondim"))
      pr => cfg%ocean%meke%min_gamma2
      call g%add(nml_real("min_gamma2", pr, "Floor on gamma_b^2/gamma_t^2", units="nondim"))
      pr => cfg%ocean%meke%uscale
      call g%add(nml_real("uscale", pr, "Background eddy velocity scale for bottom drag", units="m/s"))
      pr => cfg%ocean%meke%dtscale
      call g%add(nml_real("dtscale", pr, "Time-stepping acceleration factor", units="nondim"))
      pr => cfg%ocean%meke%khth_fac
      call g%add(nml_real("khth_fac", pr, "Geom-mean kh -> VarMix KhTh factor (0 ⇒ inert)", units="nondim"))
      pr => cfg%ocean%meke%khtr_fac
      call g%add(nml_real("khtr_fac", pr, "Geom-mean kh -> VarMix KhTr factor (0 ⇒ inert)", units="nondim"))
      pl => cfg%ocean%meke%backscatter
      call g%add(nml_logical("backscatter", pl, &
                             "Enable MEKE -> momentum harmonic backscatter (negative viscosity); default off ⇒ bit-identical"))
      pr => cfg%ocean%meke%backscatter_visc_coeff_ku
      call g%add(nml_real("backscatter_visc_coeff_ku", pr, &
                          "MEKE_VISCOSITY_COEFF_KU: harmonic backscatter efficiency Ku=coeff*sqrt(2*gt2*E)*Lmix (0 ⇒ inert)", &
                          units="nondim"))
      pr => cfg%ocean%meke%khmeke_fac
      call g%add(nml_real("khmeke_fac", pr, "meke%kh -> MEKE self-diffusivity factor", units="nondim"))
      pr => cfg%ocean%meke%advection_factor
      call g%add(nml_real("advection_factor", pr, "Barotropic-transport advection scaling (0 ⇒ off)", &
                          units="nondim"))
      pr => cfg%ocean%meke%cdrag
      call g%add(nml_real("cdrag", pr, "Bottom drag coefficient for MEKE", units="nondim"))
      pl => cfg%ocean%meke%use_bbl_drag
      call g%add(nml_logical("use_bbl_drag", pl, &
                             "Add resolved |u_bed|^2 to the MEKE bottom-drag rate (default off ⇒ bit-identical)"))
      pr => cfg%ocean%meke%alpha_deform
      call g%add(nml_real("alpha_deform", pr, "Weight on deformation length scale", units="nondim"))
      pr => cfg%ocean%meke%alpha_rhines
      call g%add(nml_real("alpha_rhines", pr, "Weight on Rhines length scale (v1 default 0 ⇒ inert)", &
                          units="nondim"))
      pr => cfg%ocean%meke%alpha_eady
      call g%add(nml_real("alpha_eady", pr, "Weight on Eady length scale (needs VarMix SN)", units="nondim"))
      pr => cfg%ocean%meke%alpha_frict
      call g%add(nml_real("alpha_frict", pr, "Weight on frictional-arrest length scale", units="nondim"))
      pr => cfg%ocean%meke%alpha_grid
      call g%add(nml_real("alpha_grid", pr, "Weight on grid length scale", units="nondim"))

      call schema%add_group(g)
   end subroutine register_meke

   module subroutine register_porous(cfg, schema)
      !! `&ocean_porous` (Adcroft 2013 porous barriers: subgrid
      !! sill/strait narrowing of the C-grid transport face widths).
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      logical, pointer :: pl
      real(wp), pointer :: pr
      character(len=:), pointer :: ps

      g%name = "ocean_porous"
      g%doc = "Porous barriers: subgrid sill/strait blocking of the C-grid face widths."

      pl => cfg%ocean%porous%enable
      call g%add(nml_logical("enable", pl, &
                             "Master switch (ocean multilayer path only; "// &
                             "fails loud with bt_halo > 0 or wetdry enable)"))
      ps => cfg%ocean%porous%source
      call g%add(nml_enum("source", ps, &
                          "Along-face bathymetry-statistics source "// &
                          "('resolved' is a wet-gated proxy; 'file' is deferred)", &
                          allowed=[character(len=15) :: "resolved", "file"]))
      ps => cfg%ocean%porous%eta_interp
      call g%add(nml_enum("eta_interp", ps, &
                          "Interface height at the velocity point "// &
                          "(PORBAR_ETA_INTERP); 'max' is the LEAST blocking, "// &
                          "'min' the most", &
                          allowed=[character(len=15) :: "max", "min", &
                                   "arithmetic", "harmonic"]))
      pr => cfg%ocean%porous%masking_depth
      call g%add(nml_real("masking_depth", pr, &
                          "Faces shallower than this stay fully open "// &
                          "(PORBAR_MASKING_DEPTH, positive below the surface)", &
                          units="m"))

      call schema%add_group(g)
   end subroutine register_porous

   module subroutine register_foxkemper(cfg, schema)
      !! `&ocean_foxkemper` (Fox-Kemper et al. 2008/2011 mixed-layer-eddy
      !! restratification, capability B5).
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      logical, pointer :: pl
      real(wp), pointer :: pr

      g%name = "ocean_foxkemper"
      g%doc = "Fox-Kemper mixed-layer-eddy restratification (B5)."

      pl => cfg%ocean%foxkemper%enable
      call g%add(nml_logical("enable", pl, &
                             "Master switch (off => bit-identity)"))
      pr => cfg%ocean%foxkemper%ce
      call g%add(nml_real("ce", pr, "FK08 coefficient Ce (0.06-0.08)"))
      pr => cfg%ocean%foxkemper%f_floor
      call g%add(nml_real("f_floor", pr, "|f| regularisation floor", units="1/s"))
      pr => cfg%ocean%foxkemper%mld_decay_time
      call g%add(nml_real("mld_decay_time", pr, &
                          "Running-mean MLD filter time-scale (0 = off, instantaneous MLD)", &
                          units="s"))
      pr => cfg%ocean%foxkemper%tail_dh
      call g%add(nml_real("tail_dh", pr, "mu cubic-tail extension (0 = exact mu)"))
      pl => cfg%ocean%foxkemper%use_mom_mixrate
      call g%add(nml_logical("use_mom_mixrate", pl, &
                             "FK11 momentum-mixrate timescale (PRODUCTION-RECOMMENDED) vs bare Ce/|f|"))
      pl => cfg%ocean%foxkemper%resolution_taper
      call g%add(nml_logical("resolution_taper", pl, &
                             "B2 res_fn hook (hard error if on without B2)"))
      pl => cfg%ocean%foxkemper%use_bodner
      call g%add(nml_logical("use_bodner", pl, &
                             "Bodner 2023 frontogenesis-arrest MLE (overrides ce/mixrate)"))
      pr => cfg%ocean%foxkemper%cr
      call g%add(nml_real("cr", pr, "Bodner 2023 efficiency coefficient Cr (0 = off)"))
      pr => cfg%ocean%foxkemper%bodner_mstar
      call g%add(nml_real("bodner_mstar", pr, "Bodner mechanical (u*) weight in w'u'"))
      pr => cfg%ocean%foxkemper%bodner_nstar
      call g%add(nml_real("bodner_nstar", pr, "Bodner convective (w*) weight in w'u'"))
      pr => cfg%ocean%foxkemper%min_wstar2
      call g%add(nml_real("min_wstar2", pr, "Floor on w'u' (1/0 armour)", units="m^2/s^2"))

      call schema%add_group(g)
   end subroutine register_foxkemper

   module subroutine register_ocean_coriolis(cfg, schema)
      !! `&ocean_coriolis_nml`: Coriolis-advection scheme selector.
      !! `form` enum mirrors `parse_pv_variant` in rdb_coriolis_adv
      !! (canonical names; short aliases hk/energy are not advertised).
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      character(len=:), pointer :: ps
      logical, pointer :: pl

      g%name = "ocean_coriolis"
      g%doc = "Coriolis-advection scheme selector."
      ps => cfg%ocean%coriolis%form
      call g%add(nml_enum("form", ps, "Coriolis-advection variant", &
                          allowed=[character(len=15) :: "sadourny", "sadourny_hk", &
                                   "sadourny_energy"]))
      ps => cfg%ocean%coriolis%pv_adv_scheme
      call g%add(nml_enum("pv_adv_scheme", ps, &
                          "PV face interpolation (Sadourny path): centered (default) "// &
                          "or weno3/weno5/weno7 (WENO-Z); weno5/weno7 need nghost>=4/5", &
                          allowed=[character(len=8) :: "centered", "weno3", &
                                   "weno5", "weno7"]))
      pl => cfg%ocean%coriolis%use_state_fluxes
      call g%add(nml_logical("use_state_fluxes", pl, &
                             "mom6-corrector CorAdv consumes continuity's renormalised "// &
                             "mass fluxes (MOM6 mass-consistent uh/vh)"))
      pl => cfg%ocean%coriolis%bound_coriolis
      call g%add(nml_logical("bound_coriolis", pl, &
                             "clamp the Coriolis accel to the (f+zeta)*v velocity-form "// &
                             "range (MOM6 BOUND_CORIOLIS; sadourny_energy, or sadourny "// &
                             "with centred pv_adv_scheme)"))
      ps => cfg%ocean%coriolis%corner_h
      call g%add(nml_enum("corner_h", ps, "PV corner-thickness construction (energy scheme)", &
                          allowed=[character(len=9) :: "cell_mean", "mom6_area"]))
      call schema%add_group(g)
   end subroutine register_ocean_coriolis

   module subroutine register_ocean_bt(cfg, schema)
      !! `&ocean_bt_nml`: split-explicit barotropic substep controls.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      integer, pointer :: pi
      real(wp), pointer :: pr
      logical, pointer :: pl
      character(len=:), pointer :: ps

      g%name = "ocean_bt"
      g%doc = "Split-explicit barotropic substep controls."
      pi => cfg%ocean%bt%n_inner
      call g%add(nml_int("n_inner", pi, &
                         "Barotropic substeps per outer step (0 = unsplit)", min=0))
      pl => cfg%ocean%bt%auto_n_inner
      call g%add(nml_logical("auto_n_inner", pl, &
                             "Derive n_inner from the gravity-wave CFL at setup"))
      pr => cfg%ocean%bt%cfl_bt_safety
      call g%add(nml_real("cfl_bt_safety", pr, &
                          "Safety fraction on the BT CFL when auto_n_inner"))
      pr => cfg%ocean%bt%dtbt
      call g%add(nml_real("dtbt", pr, &
                      "OM3-style BT step (MOM6 DTBT): <0 = fraction of max stable, 0 = -0.98, >0 = seconds; replaces auto_n_inner"))
      pr => cfg%ocean%bt%dtbt_reset_period
      call g%add(nml_real("dtbt_reset_period", pr, &
                          "Seconds between dtbt recomputes: 0 = every step, >0 = period, <0 = init only"))
      pr => cfg%ocean%bt%bebt
      call g%add(nml_real("bebt", pr, &
                          "Forward-velocity-projection weight (MOM6 BEBT; default 0.1 = MOM6)"))
      pl => cfg%ocean%bt%use_cont_type
      call g%add(nml_logical("use_cont_type", pl, &
                             "Use the BT_cont flux-bounded closure (MOM6 USE_BT_CONT_TYPE)"))
      pl => cfg%ocean%bt%cont_corr_bounds
      call g%add(nml_logical("cont_corr_bounds", pl, &
                             "Use BT_cont flux limits for the eta-correction bound"))
      pl => cfg%ocean%bt%upstream_h_face
      call g%add(nml_logical("upstream_h_face", pl, &
                             "Use per-face upstream-PPM column-sum thickness in the BT chain"))
      pl => cfg%ocean%bt%correction_h_weighted
      call g%add(nml_logical("correction_h_weighted", pl, &
                             "RETIRED h-weighted BT-corrector fold (refused when set)", &
                             dead_on_ocean_path="RETIRED -- the h-weighted barotropic-"// &
                             "correction fold was energy-non-conserving (a positive "// &
                             "0.5*D^2*H*(kappa-1) source plus shear feedback) and MOM6 "// &
                             "has no such fold; setting it .true. is a fail-loud "// &
                             "configure error (validate_config). Drag-aware damping "// &
                             "is visc_rem_chain (it does not re-weight this fold -- "// &
                             "correction_visc_rem, which used to, is itself retired)."))
      pl => cfg%ocean%bt%visc_rem_chain
      call g%add(nml_logical("visc_rem_chain", pl, &
                             "PR-3 (D1): ONE switch for exactly MOM6's visc_rem/av_rem/"// &
                             "bt_rem set -- equivalent to switching on the visc_rem "// &
                             "producer plus forcing_visc_rem + renorm_visc_rem + "// &
                             "bt_rem_from_visc_rem all at once (never a superset); the "// &
                             "BT-correction fold stays UNIFORM (MOM6 accel_layer_u never "// &
                             "weights it; correction_visc_rem, which used to, is retired). "// &
                             "strong_drag/rescale_strong_drag stay separate keys (their "// &
                             "own MOM6 params). Still requires ocean_vdiff_nml "// &
                             "implicit_drag or bbl_glue, else visc_rem is inert "// &
                             "(=1, warned)."))
      pl => cfg%ocean%bt%correction_visc_rem
      call g%add(nml_logical("correction_visc_rem", pl, &
                             "RETIRED (refused when set)", &
                             dead_on_ocean_path="RETIRED -- MOM6's accel_layer_u applies "// &
                             "the BT-correction acceleration UNIFORMLY across every layer, "// &
                             "then the SAME implicit "// &
                             "friction the glue uses distributes it -- never twice. This "// &
                             "fold re-weighted it a second time by visc_rem/<visc_rem>_h, "// &
                             "an unbounded ratio that NaNs the 1-degree Southern Ocean z* "// &
                             "open-step case under bbl_glue at step ~40. Setting it "// &
                             ".true. is a fail-loud configure error (validate_config). "// &
                             "Use visc_rem_chain instead."))
      pl => cfg%ocean%bt%bt_rem_from_visc_rem
      call g%add(nml_logical("bt_rem_from_visc_rem", pl, &
                             "bt_rem_u/v = mask*av_rem**(1/n_inner), av_rem the frhat-"// &
                             "weighted depth mean of visc_rem (MOM6's barotropic "// &
                             "solver); self-sufficient (the producer runs whenever "// &
                             "this is on), mutually exclusive with substep_drag and "// &
                             "bt_halo > 0. An equivalent subset of visc_rem_chain, kept "// &
                             "for granular testing -- prefer visc_rem_chain."))
      pl => cfg%ocean%bt%strong_drag
      call g%add(nml_logical("strong_drag", pl, &
                             "MOM6 BT_STRONG_DRAG: rational-approximation bt_rem form "// &
                             "n_inner*av_rem/(1+(n_inner-1)*av_rem) instead of the plain "// &
                             "power; requires bt_rem_from_visc_rem (or visc_rem_chain)"))
      pl => cfg%ocean%bt%rescale_strong_drag
      call g%add(nml_logical("rescale_strong_drag", pl, &
                             "MOM6 RESCALE_STRONG_DRAG: rescale the BT-correction "// &
                             "increment by min(bt_rem**n_inner/av_rem, 1.0); requires "// &
                             "strong_drag"))
      ps => cfg%ocean%bt%split_scheme
      call g%add(nml_enum("split_scheme", ps, &
                          "Outer split-explicit time scheme: pred_corr (DEFAULT; MOM6 "// &
                          "predictor-corrector, slow tendencies on the u_av/h_av step "// &
                          "time-means, forward-backward gravity-wave pairing; lifts the "// &
                          "internal-wave dt ceiling) or ssp_rk2 (EXPERIMENTAL; two-stage "// &
                          "SSP average, widest envelope — the only scheme wired through "// &
                          "eulerian_z, wet/dry and dt_tracer_advect_ratio>1 — but it "// &
                          "spuriously grows internal gravity waves out of a stratified "// &
                          "REST state, En 2.992E-05 vs 1.739E-09 at day 25 on "// &
                          "resting_stratified_channel.nml; a (omega*dt)^4 noise floor, "// &
                          "so forced viscous runs sit decades above it and quiescent or "// &
                          "long spin-up runs do not)", &
                          allowed=[character(len=9) :: "pred_corr", "ssp_rk2"], &
                          retired=[character(len=9) :: "mom6_pc", "split_rk2"], &
                          retired_hint="the MOM6 predictor-corrector is now spelled "// &
                          "'pred_corr' and is the DEFAULT, so this key can also simply "// &
                          'be deleted; write split_scheme = "pred_corr" to keep it '// &
                          "pinned"))
      pr => cfg%ocean%bt%pc_be
      call g%add(nml_real("pc_be", pr, &
                          "pred_corr predictor fraction BE (MOM6 BE, 0.6 reference)", &
                          min=0.0_wp, max=1.0_wp))
      pl => cfg%ocean%bt%renorm_visc_rem
      call g%add(nml_logical("renorm_visc_rem", pl, &
                             "gamma-weighted continuity transport-matching inversion "// &
                             "(MOM6 u_cor = u + du*visc_rem); self-sufficient (the "// &
                             "producer runs whenever this is on). An equivalent subset "// &
                             "of visc_rem_chain, kept for granular testing -- prefer "// &
                             "visc_rem_chain."))
      pl => cfg%ocean%bt%forcing_visc_rem
      call g%add(nml_logical("forcing_visc_rem", pl, &
                             "MOM6 wt_u parity: h*visc_rem-weight the BT forcing "// &
                             "depth-mean so friction-damped (glued) layers do not "// &
                             "force the fast loop; self-sufficient (the producer runs "// &
                             "whenever this is on). An equivalent subset of "// &
                             "visc_rem_chain, kept for granular testing -- prefer "// &
                             "visc_rem_chain."))
      ps => cfg%ocean%bt%frhat_scheme
      call g%add(nml_enum("frhat_scheme", ps, &
                          "Per-layer barotropic face-thickness closure for every depth "// &
                          "mean that reads a face thickness (derive_bt_from_layers, "// &
                          "face_depth_mean_u/v, face_depth_mean_rem_u/v's wt_u, "// &
                          "apply_bt_correction's open/visc_rem folds, and av_rem via "// &
                          "face_depth_mean_u/v): arithmetic (DEFAULT, bit-identical -- "// &
                          "the plain two-abutting-cell mean) or hybrid (MOM6 btcalc "// &
                          "HVEL_SCHEME=HYBRID -- arithmetic "// &
                          "mean above the shallower column's bed, harmonic mean below "// &
                          "it, suppressing a thin partial-bed layer's weight next to a "// &
                          "thick abutting one)", &
                          allowed=[character(len=10) :: "arithmetic", "hybrid"]))
      pl => cfg%ocean%bt%correction_bc_pgf
      call g%add(nml_logical("correction_bc_pgf", pl, &
                             "Per-layer baroclinic-PGF retro-correction for the eta change "// &
                             "(requires &ocean_pgf_nml form='fv_mom6')"))
      pl => cfg%ocean%bt%bc_pgf_forcing
      call g%add(nml_logical("bc_pgf_forcing", pl, &
                             "Force the BT substep with the depth mean of the full slow layer "// &
                             "PGF, shedding only the free-surface term the PGF itself carries "// &
                             "(MOM6 BT_force/eta_PF); .false. = legacy split that discarded the "// &
                             "depth-mean baroclinic PGF (no JEBAR)"))
      pl => cfg%ocean%bt%substep_drag
      call g%add(nml_logical("substep_drag", pl, &
                             "Apply a per-substep BT velocity damping factor"))
      pl => cfg%ocean%bt%substep_zeta_ke
      call g%add(nml_logical("substep_zeta_ke", pl, &
                             "Integrate live zeta_bt + KE-gradient in the BT fast loop "// &
                             "(.false. = MOM6 parity: planetary Coriolis only, "// &
                             "zeta/KE frozen in the slow forcing)"))
      pl => cfg%ocean%bt%wave_drag
      call g%add(nml_logical("wave_drag", pl, &
                             "Barotropic linear wave drag master switch (MOM6 BT_LINEAR_WAVE_DRAG)"))
      ps => cfg%ocean%bt%wave_drag_form
      call g%add(nml_enum("wave_drag_form", ps, "Wave-drag r_H filler", &
                          allowed=[character(len=15) :: "uniform", "roughness_proxy", "file"]))
      pr => cfg%ocean%bt%wave_drag_scale
      call g%add(nml_real("wave_drag_scale", pr, &
                          "Global tuning multiplier on r_H (MOM6 BT_WAVE_DRAG_SCALE)"))
      pr => cfg%ocean%bt%wave_drag_r_uniform
      call g%add(nml_real("wave_drag_r_uniform", pr, &
                          "Piston velocity r_H for wave_drag_form='uniform'", units="m/s"))
      pr => cfg%ocean%bt%wave_drag_kappa
      call g%add(nml_real("wave_drag_kappa", pr, &
                          "Topographic wavenumber for wave_drag_form='roughness_proxy'", &
                          units="1/m"))
      pr => cfg%ocean%bt%wave_drag_n_bot
      call g%add(nml_real("wave_drag_n_bot", pr, &
                          "Reference bottom N for wave_drag_form='roughness_proxy'", units="1/s"))
      pr => cfg%ocean%bt%wave_drag_h2_max
      call g%add(nml_real("wave_drag_h2_max", pr, &
                          "Ceiling on <h^2> proxy for wave_drag_form='roughness_proxy'", &
                          units="m^2"))
      ps => cfg%ocean%bt%wave_drag_file
      call g%add(nml_string("wave_drag_file", ps, &
                            "Reserved for PR-14 (MOM6 BT_WAVE_DRAG_FILE); unused today"))
      ps => cfg%ocean%bt%wave_drag_var
      call g%add(nml_string("wave_drag_var", ps, &
                            "Reserved for PR-14 (MOM6 BT_WAVE_DRAG_VAR); unused today"))
      pi => cfg%ocean%bt%bt_halo
      call g%add(nml_int("bt_halo", pi, &
                         "Wide-halo BT march-in width "// &
                         "(-1 = auto: resolves to 0, the march-in is opt-in; "// &
                         "0 = explicit off, bit-identical)", &
                         min=BT_HALO_AUTO_SENTINEL))
      call schema%add_group(g)
   end subroutine register_ocean_bt

   module subroutine register_ocean_hdiff(cfg, schema)
      !! `&ocean_hdiff_nml`: along-coordinate tracer Laplacian
      !! (`rdb_ocean_hdiff_tracer`).  Not the neutral/isopycnal path —
      !! that is `&ocean_redi_nml`, already reachable.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      real(wp), pointer :: pr

      g%name = "ocean_hdiff"
      g%doc = "Along-coordinate (not neutral) constant-coefficient tracer diffusion."
      pr => cfg%ocean%hdiff%kappa_h
      call g%add(nml_real("kappa_h", pr, &
                          "Horizontal tracer diffusivity, along the model coordinate "// &
                          "(0 disables; ocean path only — coastal's equivalent knob is "// &
                          "top-level hdiff_kappa)", &
                          units="m^2/s", min=0.0_wp))
      call schema%add_group(g)
   end subroutine register_ocean_hdiff

   module subroutine register_ocean_hvisc(cfg, schema)
      !! `&ocean_hvisc_nml`: lateral-viscosity closure + coefficients.
      !! `lateral_closure` enum mirrors `parse_lateral_closure` in
      !! rdb_ocean_lateral_mix.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      real(wp), pointer :: pr
      real(wp), pointer :: pa(:)
      integer, pointer :: pi
      logical, pointer :: pl
      character(len=:), pointer :: ps

      g%name = "ocean_hvisc"
      g%doc = "Lateral-viscosity closure + coefficients."
      ps => cfg%ocean%hvisc%lateral_closure
      call g%add(nml_enum("lateral_closure", ps, "Lateral-mixing closure tag", &
                          allowed=[character(len=12) :: "none", "leith", "smagorinsky", &
                                   "smag", "biharmonic", "leith_biharm"]))
      pr => cfg%ocean%hvisc%c_smag
      call g%add(nml_real("c_smag", pr, "Smagorinsky coefficient"))
      pr => cfg%ocean%hvisc%c_leith
      call g%add(nml_real("c_leith", pr, "Leith coefficient"))
      pr => cfg%ocean%hvisc%kh_vel_scale
      call g%add(nml_real("kh_vel_scale", pr, "Velocity scale for the resolution viscosity floor", &
                          units="m/s"))
      pr => cfg%ocean%hvisc%ah_bg
      call g%add(nml_real("ah_bg", pr, "Background harmonic viscosity (<0 = derive)", units="m^2/s"))
      pr => cfg%ocean%hvisc%ah_max
      call g%add(nml_real("ah_max", pr, "Cap on per-face harmonic viscosity", units="m^2/s"))
      pl => cfg%ocean%hvisc%smag_ah
      call g%add(nml_logical("smag_ah", pl, "Flow-aware biharmonic viscosity (MOM6 SMAGORINSKY_AH)"))
      pr => cfg%ocean%hvisc%smag_bi_const
      call g%add(nml_real("smag_bi_const", pr, "Biharmonic Smagorinsky constant (MOM6 SMAG_BI_CONST)"))
      pr => cfg%ocean%hvisc%c_leith_bi
      call g%add(nml_real("c_leith_bi", pr, "Biharmonic Leith constant (MOM6 LEITH_BI_CONST)"))
      pr => cfg%ocean%hvisc%nu_4_bg
      call g%add(nml_real("nu_4_bg", pr, "Background biharmonic viscosity floor", units="m^4/s"))
      pr => cfg%ocean%hvisc%nu_4_max
      call g%add(nml_real("nu_4_max", pr, "Cap on per-face biharmonic viscosity", units="m^4/s"))
      pr => cfg%ocean%hvisc%nu_h
      call g%add(nml_real("nu_h", pr, "Constant horizontal eddy viscosity", units="m^2/s"))
      pr => cfg%ocean%hvisc%nu_4
      call g%add(nml_real("nu_4", pr, "Constant biharmonic eddy viscosity", units="m^4/s"))
      pl => cfg%ocean%hvisc%no_slip
      call g%add(nml_logical("no_slip", pl, &
                             "Coastal lateral BC: .false.=free-slip (×wet_q), .true.=no-slip (×(2-wet_q))"))
      pl => cfg%ocean%hvisc%stress_tensor
      call g%add(nml_logical("stress_tensor", pl, &
                             "MOM6 thickness-weighted stress-div operator + per-cell CFL + coast-mask"))
      pr => cfg%ocean%hvisc%bound_coef
      call g%add(nml_real("bound_coef", pr, "Per-cell viscosity-CFL safety coefficient (MOM6 HORVISC_BOUND_COEF)"))
      pl => cfg%ocean%hvisc%bound_kh
      call g%add(nml_logical("bound_kh", pl, &
                             "Per-face harmonic viscosity CFL clamp on the velocity-Laplacian paths (MOM6 BOUND_KH)"))
      pl => cfg%ocean%hvisc%resoln_scaled_visc
      call g%add(nml_logical("resoln_scaled_visc", pl, &
                  "Scale dynamic LAPLACIAN viscosity (Leith/Smag_KH A_h, not biharmonic) by VarMix Res_fn (MOM6 RESOLN_SCALED_KH)"))
      pr => cfg%ocean%hvisc%kh_vel_scale_live
      call g%add(nml_real("kh_vel_scale_live", pr, &
                          "Live velocity-scale viscosity Kh=U*dx*|u| (0=off; distinct from kh_vel_scale floor)", &
                          units="m/s"))
      pr => cfg%ocean%hvisc%kh_sin_lat
      call g%add(nml_real("kh_sin_lat", pr, &
                          "Latitude Laplacian viscosity floor kh_sin_lat*|sin(lat)|^kh_pwr_of_sine (0=off; MOM6 KH_SIN_LAT)", &
                          units="m^2/s"))
      pr => cfg%ocean%hvisc%kh_pwr_of_sine
      call g%add(nml_real("kh_pwr_of_sine", pr, "Exponent of |sin(lat)| in the kh_sin_lat floor (MOM6 KH_PWR_OF_SINE)"))
      pr => cfg%ocean%hvisc%ah_vel_scale
      call g%add(nml_real("ah_vel_scale", pr, &
                          "Biharmonic velocity-scale floor Ah>=U*Delta^3 under smag_ah (0=off; MOM6 AH_VEL_SCALE)", &
                          units="m/s"))
      pl => cfg%ocean%hvisc%bound_coriolis_biharm
      call g%add(nml_logical("bound_coriolis_biharm", pl, &
                             "Coriolis bound on the Smagorinsky biharmonic (needs smag_ah; MOM6 BOUND_CORIOLIS_BIHARM)"))
      pr => cfg%ocean%hvisc%bound_coriolis_vel
      call g%add(nml_real("bound_coriolis_vel", pr, &
                          "Velocity scale of the bound_coriolis_biharm term (MOM6 BOUND_CORIOLIS_VEL)", &
                          units="m/s"))
      pr => cfg%ocean%hvisc%kh_aniso
      call g%add(nml_real("kh_aniso", pr, &
                          "Anisotropic Laplacian viscosity magnitude (MOM6 KH_ANISO; stress_tensor path only)", &
                          units="m^2/s"))
      pi => cfg%ocean%hvisc%aniso_mode
      call g%add(nml_int("aniso_mode", pi, &
                         "Anisotropy direction mode (0=grid-relative aniso_dir; MOM6 ANISOTROPIC_MODE)", min=0))
      pa => cfg%ocean%hvisc%aniso_dir
      call g%add(nml_real_array("aniso_dir", pa, &
                                "Anisotropy direction (n1,n2) in grid i,j components (MOM6 ANISO_GRID_DIR)"))
      call schema%add_group(g)
   end subroutine register_ocean_hvisc

   module subroutine register_ocean_isopycnal(cfg, schema)
      !! `&ocean_isopycnal_nml`: grounding-stability controls for the
      !! Lagrangian vertical coordinate.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      real(wp), pointer :: pr
      logical, pointer :: pl

      g%name = "ocean_isopycnal"
      g%doc = "Lagrangian grounding-stability controls."
      pr => cfg%ocean%isopycnal%angstrom_h
      call g%add(nml_real("angstrom_h", pr, &
                          "Minimum-thickness floor on the Lagrangian continuity h-update (MOM6 Angstrom_H analogue)", &
                          units="m"))
      pl => cfg%ocean%isopycnal%reset_vanished_u
      call g%add(nml_logical("reset_vanished_u", pl, &
                             "Zero face velocity when layer vanished on both adjacent cells"))
      pl => cfg%ocean%isopycnal%cfl_ignore_vanished
      call g%add(nml_logical("cfl_ignore_vanished", pl, &
                             "Exclude vanished layers from MaxCFL / panic / CFL truncation"))
      pl => cfg%ocean%isopycnal%pgf_skip_nonoverlap
      call g%add(nml_logical("pgf_skip_nonoverlap", pl, &
                             "Zero the face PGF where a grounded layer's z-extents do not "// &
                             "overlap across the face AND the layer is at the floor "// &
                             "(<= 2*max(angstrom_h, H_VANISHED)) on one side "// &
                             "(VCOORD_LAGRANGIAN only, mont / "// &
                             "fv_lite / fv_wright / fv_mom6; default ON — kills the spurious "// &
                             "at-rest grounded-layer pressure gradient)"))
      pl => cfg%ocean%isopycnal%conservative_floor
      call g%add(nml_logical("conservative_floor", pl, &
                             "Conservative min-thickness borrow (replaces the injecting angstrom_h floor)"))
      pl => cfg%ocean%isopycnal%check_h_positive
      call g%add(nml_logical("check_h_positive", pl, &
                             "DEBUG: abort on the first negative h_layer, naming the stage + (i,j,k)"))
      call schema%add_group(g)
   end subroutine register_ocean_isopycnal

end submodule rdb_config_register_lateral
