### &sim_nml

Simulation-regime selector.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `sim_type` | `"ocean"` |  | Simulation regime (ocean is the only regime this build ships) |

### &grid_nml

Structured-grid geometry.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `nx` | `200` |  | Number of physical cells in x |
| `ny` | `1` |  | Number of physical cells in y |
| `dx` | `0.1000000000E+01` | m | Cell size in x |
| `dy` | `0.1000000000E+01` | m | Cell size in y |
| `nghost` | `3` |  | Ghost cells on each side |

### &time_nml

Time-integration controls.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `t_end` | `0.1000000000E+01` |  | Simulation end time (in time_unit) |
| `cfl` | `0.4500000000E+00` |  | CFL number for adaptive timestep |
| `dt_max` | `0.1000000000E+11` | s | Maximum allowable timestep |
| `dt_fixed` | `0.0000000000E+00` | s | Fixed timestep (0 = adaptive CFL) |
| `cfl_interval` | `1` |  | Recompute CFL timestep every N steps |
| `time_unit` | `"s"` |  | Unit for the long-time fields: s/min/hr/day/year |

### &mpi_nml

MPI domain decomposition.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `px` | `1` |  | MPI processes in x |
| `py` | `1` |  | MPI processes in y |

### &logging_nml

Logger verbosity + status cadence.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `log_level` | `"info"` |  | Log verbosity |
| `status_interval` | `0.0000000000E+00` |  | Status-print cadence (in time_unit; 0 = every 100 steps) |

### &vcoord_nml

Vertical-coordinate + ALE-remap controls.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `vcoord_type` | `"sigma"` |  | Vertical coordinate type |
| `thickness_config` | `"sigma"` |  | Initial layer-thickness profile (ocean path) |
| `remap_method` | `"ppm"` |  | Vertical remapping method |
| `zstar_h_surf_target` | `0.0000000000E+00` | m | z*-full: target surface-layer thickness (0 = auto) |
| `zstar_stretching` | `"log"` |  | z*-full: surface-concentration stretching |
| `zstar_h_min` | `0.1000000000E-03` | m | z*-full: vanishing-layer floor |
| `zstar_n_surf` | `0` |  | z*-full: number of fine near-surface layers (0 = auto) |
| `rho_ref_pressure` | `0.2000000000E+08` | Pa | rho-coord: reference pressure for potential density |
| `rho_target_light` | `0.1020000000E+04` | kg/m^3 | rho-coord: lightest (surface) target density |
| `rho_target_dense` | `0.1030000000E+04` | kg/m^3 | rho-coord: densest (bed) target density |
| `regrid_time_scale` | `0.0000000000E+00` | s | ALE regrid grid time-filter timescale (0 = jump to target) |
| `remap_vel_conserve_ke` | `.false.` |  | ALE velocity remap: KE-conserving baroclinic-anomaly rescale |
| `remap_boundary_extrap` | `.false.` |  | ALE remap: linear-exact one-sided reconstruction in the k=1/k=nz boundary cells (MOM6 BOUNDARY_EXTRAPOLATION) |
| `remap_nonuniform_weights` | `.false.` |  | ALE remap: non-uniform-grid PLM slope + PPM edge weights (Colella-Woodward 1984 eqs 1.6-1.8) instead of the equal-thickness specialisations |
| `remap_check_preconditions` | `.false.` |  | ALE remap: fail loud when a column violates the overlap sweep's preconditions (non-negative thicknesses, matching column totals) |

### &physics_nml

Barotropic physics: bottom drag, wind stress, Coriolis.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `manning_n` | `0.0000000000E+00` |  | Manning roughness coefficient |
| `wind_stress_x` | `0.0000000000E+00` | Pa | Surface wind stress in x |
| `wind_stress_y` | `0.0000000000E+00` | Pa | Surface wind stress in y |
| `coriolis_f` | `0.0000000000E+00` | 1/s | Coriolis parameter f |

### &output_nml

File output, I/O server, forcing/restart/gauge files.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `output_to_file` | `.false.` |  | Enable file output (NetCDF snapshots) |
| `output_dir` | `"./output"` |  | Directory for output files |
| `restart_interval` | `0.0000000000E+00` | s | Time between restart writes (0 = none) |
| `compress_output` | `.false.` |  | Deflate-compress NetCDF output |
| `compress_level` | `1` |  | Deflate level (1=fast, 9=max) |
| `use_io_server` | `.false.` |  | Dedicate one MPI rank per node as I/O server |
| `bathymetry_file` | `""` |  | NetCDF bathymetry file (empty = flat) |
| `restart_file` | `""` |  | Restart file for warm start (empty = cold) |

### &boundary_nml

Boundary conditions, tidal forcing, inflow/discharge/sponge/nesting.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `bc_west` | `"wall"` |  | West boundary type |
| `bc_east` | `"wall"` |  | East boundary type |
| `bc_south` | `"wall"` |  | South boundary type |
| `bc_north` | `"wall"` |  | North boundary type |
| `n_tidal_constituents` | `0` |  | Active tidal constituents (0 = legacy single) |
| `tidal_amp` | `0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00` | m | Constituent amplitudes |
| `tidal_phase` | `0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00` | rad | Constituent phases |
| `tidal_omega` | `0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00` | rad/s | Constituent angular frequencies |
| `inflow_salinity` | `-0.1000000000E+01` | PSU | Inflow salinity (<0 = zero-gradient) |
| `inflow_temperature` | `-0.9990000000E+03` | degC | Inflow temperature (<0 = zero-gradient) |
| `sponge_width` | `0` |  | Sponge layer width in cells |
| `sponge_strength` | `0.0000000000E+00` | 1/s | Sponge relaxation rate |

### &tracer_nml

Salinity + temperature IC/EOS/bounds, sediment transport.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `initial_salinity` | `0.3500000000E+02` | PSU | Initial salinity (uniform IC) |
| `S_ref` | `0.0000000000E+00` | PSU | EOS reference salinity |
| `beta_S` | `0.7800000000E+00` | kg/m^3/PSU | Haline contraction coefficient |
| `S_min` | `0.0000000000E+00` | PSU | Lower physical bound for salinity |
| `S_max` | `0.4000000000E+02` | PSU | Upper physical bound for salinity |
| `kappa_S_bg` | `0.1000000000E-04` | m^2/s | Background vertical salinity diffusivity |
| `S_init_surface` | `0.0000000000E+00` | PSU | Initial surface salinity at k=nz (linear-in-layer stratified IC; needs S_init_bottom non-zero too) |
| `S_init_bottom` | `0.0000000000E+00` | PSU | Initial bed salinity at k=1 (linear-in-layer stratified IC; needs S_init_surface non-zero too) |
| `initial_temperature` | `0.1500000000E+02` | degC | Initial temperature (uniform IC) |
| `T_ref` | `0.1500000000E+02` | degC | EOS reference temperature |
| `alpha_T` | `0.1700000000E+00` | kg/m^3/degC | Thermal expansion coefficient |
| `T_min` | `-0.2000000000E+01` | degC | Lower physical bound for temperature |
| `T_max` | `0.4000000000E+02` | degC | Upper physical bound for temperature |
| `kappa_T_bg` | `0.1000000000E-04` | m^2/s | Background vertical temperature diffusivity |
| `T_init_surface` | `0.0000000000E+00` | degC | Initial surface temperature (stratified IC) |
| `T_init_bottom` | `0.0000000000E+00` | degC | Initial bed temperature (stratified IC) |

### &initial_condition_nml

Coastal initial-condition selector + parameters.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `h0` | `0.1000000000E+01` | m | Background depth (gaussian_hump IC) |

### &nonhydrostatic_nml

Non-hydrostatic / multilayer: CG-Poisson, PP81/KPP/k-eps vmix, Smagorinsky, BPG, mode split.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `nz_layers` | `2` |  | Number of sigma layers |
| `use_multilayer` | `.false.` |  | Enable coupled hydrostatic vertical layers |
| `rho_0` | `0.1000000000E+04` | kg/m^3 | Reference density for EOS |
| `kpp_ri_crit` | `0.3000000000E+00` |  | Critical bulk Richardson number for KPP BL-depth |
| `kpp_cs_nonlocal` | `0.6300000000E+01` |  | KPP non-local (counter-gradient) transport coefficient |
| `kpp_c_vt2` | `0.1800000000E+01` |  | KPP V_t^2 unresolved-turbulence coefficient (0 = off) |
| `hdiff_kappa` | `0.0000000000E+00` | m^2/s | Horizontal tracer diffusion coefficient (COASTAL path only; ocean's along-coordinate equivalent is &ocean_hdiff_nml kappa_h) |

### &ocean_kappa_shear_nml

JHL08 shear-driven interior turbulence (kappa-shear).

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `enable` | `.false.` |  | Master switch (requires use_closure + thermodynamics) |
| `ri_crit` | `0.2500000000E+00` | nondim | Critical Richardson number |
| `shearmix_rate` | `0.8900000000E-01` |  | Shear source-rate coefficient |
| `fri_curvature` | `-0.9700000000E+00` |  | Ri-function curvature in the shear source |
| `c_n` | `0.2400000000E+00` |  | TKE decay-rate coefficient vs stratification N |
| `c_s` | `0.1400000000E+00` |  | TKE decay-rate coefficient vs shear S |
| `lambda` | `0.8200000000E+00` |  | Buoyancy mixing-length-scale coefficient |
| `lz_rescale` | `0.1000000000E+01` |  | Boundary-distance length-scale rescale factor |
| `kappa_0` | `0.1000000000E-06` | m^2/s | Background diffusivity (pre-step kappa) |
| `kappa_seed` | `0.1000000000E+01` | m^2/s | Iteration seed diffusivity |
| `kappa_trunc` | `0.1000000000E-08` | m^2/s | Diffusivity truncated to 0 below this |
| `tke_bg` | `0.0000000000E+00` | m^2/s^2 | Background TKE (Q denominator floor) |
| `tol_err` | `0.1000000000E+00` |  | Picard convergence tolerance |
| `max_inner_it` | `50` |  | Inner Picard iteration cap |
| `max_substep_it` | `13` |  | Outer adaptive-substep iteration cap |
| `src_max_chg` | `0.1000000000E+02` |  | Adaptive-dt source-change tolerance band |
| `prandtl_turb` | `0.1000000000E+01` |  | Kv = prandtl_turb * Kd into the momentum solve |
| `vel_underflow` | `0.0000000000E+00` | m/s | Velocity snap-to-zero magnitude |
| `massless_merge` | `.false.` |  | Merge vanished (<H_VANISHED) layers onto the massive sub-grid before the column solve (default off; identity columns bypass) |
| `at_vertex` | `.false.` |  | Solve the JHL08 columns at C-grid corners (vorticity points) from the native face velocities, averaging corner Kd back to tracer points (MOM6 VERTEX_SHEAR; v1 = Kd only) |
| `vertex_geometric_mean` | `.false.` |  | Geometric (vs arithmetic) mean in the corner->centre Kd average (MOM6 VERTEX_SHEAR_GEOMETRIC_MEAN) |
| `vertex_geomean_kdmin` | `0.0000000000E+00` | m^2/s | Floor applied to each corner Kd before the geometric mean (inert unless vertex_geometric_mean; OM5 configs use 1e-9) |

### &ocean_slopes_nml

Isopycnal (neutral) slope diagnostics (Griffies 1998).

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `enable` | `.false.` |  | Master switch (diagnostic; default off ⇒ no-op) |
| `kd_smooth` | `0.1000000000E-05` | m^2/s | Vert-fill smoothing diffusivity (× dt fills massless layers) |
| `min_dz_for_n2` | `0.1000000000E+01` | m | Minimum layer thickness floored in the N²/drdz denominator |

### &ocean_gm_nml

Gent-McWilliams thickness diffusion (eddy bolus transport).

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `enable` | `.false.` |  | Master switch (requires &ocean_slopes_nml enable; default off ⇒ no-op) |
| `khth` | `0.0000000000E+00` | m^2/s | Thickness diffusivity KhTh (constant-fill; production 1e2-1e3) |
| `khth_max_cfl` | `0.1000000000E+00` | nondim | Fraction of the diffusive CFL the face KH may use |
| `khth_slope_max` | `0.1000000000E-01` | nondim | Slope magnitude above which the safe-streamfunction blend takes over |

### &ocean_redi_nml

Redi continuous neutral (along-isopycnal) tracer diffusion.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `enable` | `.false.` |  | Master switch (default off ⇒ no-op; augments hdiff_tracer) |
| `continuous` | `.true.` |  | Continuous variant (discontinuous deferred R4; .false. rejected) |
| `khtr` | `0.0000000000E+00` | m^2/s | Redi neutral diffusivity KhTr (production 1e2-1e3) |

### &ocean_varmix_nml

Spatially-varying GM/Redi lateral-diffusivity coefficients.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `enable` | `.false.` |  | Master switch (requires slopes + wavespeed; default off ⇒ GM uses const khth) |
| `use_visbeck` | `.false.` |  | Add the Visbeck/Eady khth_slope_cff·L²·SN baroclinicity term |
| `resoln_scaled_khth` | `.false.` |  | Scale the assembled KhTh by the resolution function |
| `resoln_scaled_khtr` | `.false.` |  | Scale the assembled KhTr by the resolution function |
| `gill_equatorial_ld` | `.true.` |  | Gill (1982) equatorial-Ld convention (factor 2 in beta_dx2) |
| `interpolate_res_fn` | `.false.` |  | Interpolate centre Res_fn to faces (else interpolate cg1, MOM6 default) |
| `kh_res_fn_power` | `2` |  | Resolution-function power p (even) |
| `kh_res_scale_coef` | `0.1000000000E+01` | nondim | Resolution-function alpha ((alpha·cg1)^p denominator coef) |
| `khth` | `0.0000000000E+00` | m^2/s | Background thickness diffusivity KhTh |
| `khtr` | `0.0000000000E+00` | m^2/s | Background tracer diffusivity KhTr (future Redi) |
| `khth_slope_cff` | `0.0000000000E+00` | nondim | Visbeck coefficient for the KhTh chain |
| `khtr_slope_cff` | `0.0000000000E+00` | nondim | Visbeck coefficient for the KhTr chain |
| `khth_min` | `0.0000000000E+00` | m^2/s | Lower clamp on KhTh |
| `khth_max` | `0.0000000000E+00` | m^2/s | Upper clamp on KhTh (<= 0 ⇒ no cap) |
| `khtr_min` | `0.0000000000E+00` | m^2/s | Lower clamp on KhTr |
| `khtr_max` | `0.0000000000E+00` | m^2/s | Upper clamp on KhTr (<= 0 ⇒ no cap) |
| `visbeck_l_scale` | `0.0000000000E+00` | m | Visbeck length scale L (m); if < 0, |L|²·areaCu |
| `visbeck_max_slope` | `0.0000000000E+00` | nondim | S² limiter scale (<= 0 ⇒ no limit) |

### &ocean_meke_nml

Prognostic mesoscale eddy kinetic energy (GM<->eddy loop).

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `enable` | `.false.` |  | Master switch (requires &ocean_gm_nml enable; default off ⇒ no-op) |
| `gmcoeff` | `-0.1000000000E+01` | nondim | PE->MEKE conversion efficiency (< 0 ⇒ off) |
| `frcoeff` | `-0.1000000000E+01` | nondim | Frictional mean->eddy conversion (< 0 ⇒ off; >=0 sources hvisc KE dissipation) |
| `bgsrc` | `0.0000000000E+00` | m^2/s^3 | Background energy source |
| `damping` | `0.0000000000E+00` | 1/s | Linear MEKE dissipation rate |
| `kh` | `-0.1000000000E+01` | m^2/s | Background lateral diffusion of MEKE (< 0 ⇒ off) |
| `k4` | `-0.1000000000E+01` | m^4/s | Background biharmonic diffusion of MEKE (< 0 ⇒ off) |
| `khcoeff` | `0.1000000000E+01` | nondim | MEKE->Kh scaling (<= 0 ⇒ closure off) |
| `cd_scale` | `0.0000000000E+00` | nondim | Bottom/column eddy-velocity ratio |
| `cb` | `0.2500000000E+02` | nondim | Coefficient in gamma_bot (bottomFac2) |
| `ct` | `0.5000000000E+02` | nondim | Coefficient in gamma_bt (barotrFac2) |
| `min_gamma2` | `0.1000000000E-03` | nondim | Floor on gamma_b^2/gamma_t^2 |
| `uscale` | `0.0000000000E+00` | m/s | Background eddy velocity scale for bottom drag |
| `dtscale` | `0.1000000000E+01` | nondim | Time-stepping acceleration factor |
| `khth_fac` | `0.0000000000E+00` | nondim | Geom-mean kh -> VarMix KhTh factor (0 ⇒ inert) |
| `khtr_fac` | `0.0000000000E+00` | nondim | Geom-mean kh -> VarMix KhTr factor (0 ⇒ inert) |
| `backscatter` | `.false.` |  | Enable MEKE -> momentum harmonic backscatter (negative viscosity); default off ⇒ bit-identical |
| `backscatter_visc_coeff_ku` | `0.0000000000E+00` | nondim | MEKE_VISCOSITY_COEFF_KU: harmonic backscatter efficiency Ku=coeff*sqrt(2*gt2*E)*Lmix (0 ⇒ inert) |
| `khmeke_fac` | `0.0000000000E+00` | nondim | meke%kh -> MEKE self-diffusivity factor |
| `advection_factor` | `0.0000000000E+00` | nondim | Barotropic-transport advection scaling (0 ⇒ off) |
| `cdrag` | `0.2500000000E-02` | nondim | Bottom drag coefficient for MEKE |
| `use_bbl_drag` | `.false.` |  | Add resolved |u_bed|^2 to the MEKE bottom-drag rate (default off ⇒ bit-identical) |
| `alpha_deform` | `0.0000000000E+00` | nondim | Weight on deformation length scale |
| `alpha_rhines` | `0.0000000000E+00` | nondim | Weight on Rhines length scale (v1 default 0 ⇒ inert) |
| `alpha_eady` | `0.0000000000E+00` | nondim | Weight on Eady length scale (needs VarMix SN) |
| `alpha_frict` | `0.0000000000E+00` | nondim | Weight on frictional-arrest length scale |
| `alpha_grid` | `0.0000000000E+00` | nondim | Weight on grid length scale |

### &ocean_tidal_mixing_nml

St-Laurent/Simmons internal-tide interior diapycnal mixing.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `enable` | `.false.` |  | Master switch (requires use_closure + thermodynamics) |
| `gamma` | `0.3333000000E+00` | nondim | Local-dissipation fraction q (GAMMA_ITIDES) |
| `mu` | `0.2000000000E+00` | nondim | Mixing efficiency Gamma_mix (MU_ITIDES) |
| `zeta` | `0.5000000000E+03` | m | Bottom decay scale (INT_TIDE_DECAY_SCALE) |
| `kd_max` | `0.1000000000E-01` | m^2/s | Per-layer physical Kd cap (<0 => no cap) |
| `prandtl_tidal` | `0.1000000000E+01` |  | Kv = prandtl_tidal * Kd |
| `min_zbot` | `0.0000000000E+00` | m | Mask off where column depth H < min_zbot |
| `e_uniform` | `0.0000000000E+00` | W m-2 | Uniform bottom internal-tide energy input E |
| `e_compute` | `.false.` |  | State-dependent E = min(TKE_coef*N_bot, e_max) (v1.1) |
| `kappa_itides` | `0.6283200000E-03` | m^-1 | Topographic wavenumber (v1.1 E recompute) |
| `kappa_h2` | `0.1000000000E+01` |  | KAPPA_H2_FACTOR (v1.1 E recompute) |
| `utide` | `0.0000000000E+00` | m/s | RMS barotropic tidal velocity (v1.1 E recompute) |
| `h2_rough` | `0.0000000000E+00` | m^2 | Sub-grid topographic roughness variance <h^2> |
| `frac_rough` | `0.1000000000E+00` |  | Roughness clamp <h^2> <= (frac_rough*H)^2 |
| `e_max` | `0.1000000000E+04` | W m-2 | TKE_itide_max cap on E |

### &ocean_conv_nml

Brunt-Vaisala-triggered convective adjustment (interior closure contributor).

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `enable` | `.false.` |  | Master switch (requires use_closure + thermodynamics) |
| `kd_conv` | `0.1000000000E+01` | m^2/s | Convective tracer diffusivity (KD_CONV) |
| `prandtl_conv` | `0.1000000000E+01` |  | Kv_conv = prandtl_conv * kd_conv |
| `n2_thresh` | `0.0000000000E+00` | s^-2 | Trigger threshold on N^2 (BV_SQR_CONV) |

### &ocean_porous_nml

Porous barriers: subgrid sill/strait blocking of the C-grid face widths.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `enable` | `.false.` |  | Master switch (ocean multilayer path only; fails loud with bt_halo > 0 or wetdry enable) |
| `source` | `"resolved"` |  | Along-face bathymetry-statistics source ('resolved' is a wet-gated proxy; 'file' is deferred) |
| `eta_interp` | `"max"` |  | Interface height at the velocity point (PORBAR_ETA_INTERP); 'max' is the LEAST blocking, 'min' the most |
| `masking_depth` | `0.0000000000E+00` | m | Faces shallower than this stay fully open (PORBAR_MASKING_DEPTH, positive below the surface) |

### &ocean_ddiff_nml

Double diffusion: salt fingering + diffusive convection (folded into the heat/salt split).

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `enable` | `.false.` |  | Master switch (requires thermodynamics) |
| `strat_param_max` | `0.2550000000E+01` |  | R_rho salt-fingering cutoff (STRAT_PARAM_MAX) |
| `kappa_ddiff_s` | `0.1000000000E-03` | m^2/s | Leading salt-fingering diffusivity K_f (KAPPA_DDIFF_S) |
| `ddiff_exp1` | `0.1000000000E+01` |  | Inner fingering exponent (DDIFF_EXP1) |
| `ddiff_exp2` | `0.3000000000E+01` |  | Outer fingering exponent (DDIFF_EXP2) |
| `param1` | `0.9090000000E+00` |  | MC76 convection exterior coeff (KAPPA_DDIFF_PARAM1) |
| `param2` | `0.4600000000E+01` |  | MC76 convection middle coeff (KAPPA_DDIFF_PARAM2) |
| `param3` | `-0.5400000000E+00` |  | MC76 convection interior coeff (KAPPA_DDIFF_PARAM3) |
| `mol_diff` | `0.1500000000E-05` | m^2/s | Molecular diffusivity scaling convection (MOL_DIFF) |
| `use_k90` | `.false.` |  | Convection form: MC76 (default) vs Kelley-90 |

### &ocean_tides_nml

Equilibrium (astronomical) body-force tide (C1) + scalar SAL (C2).

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `enable` | `.false.` |  | Master switch (requires a non-cartesian grid) |
| `use_sal` | `.false.` |  | Apply scalar self-attraction & loading (C2) |
| `beta_sal` | `0.0000000000E+00` |  | Scalar SAL factor beta (~0.085-0.12) |
| `add_nodal` | `.false.` |  | Apply the 18.6-yr nodal f/u corrections |
| `constituents` | `"M2 S2 N2 K2 K1 O1 P1 Q1"` |  | Active constituent list (e.g. 'M2 S2 K1 O1') |
| `ref_date` | `"1900-01-01"` |  | Astronomical reference date YYYY-MM-DD |
| `nodal_ref_date` | `""` |  | Nodal reference date ('' => ref_date) |

### &ocean_psurf_nml

Atmospheric surface-pressure loading / inverse barometer (PR-17).

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `enable` | `.false.` |  | Master switch (split-solver only; requires &ocean_forcing_nml enable_components=.true.) |
| `in_eos` | `.false.` |  | Also feed the surface load to the equation of state as the top-of-column pressure p_top (requires enable=.true.; refused with the unported pressure builders) |
| `p_surf_const` | `0.0000000000E+00` |  | Uniform atmospheric surface pressure (Pa) seeded into p_surf_atm (uniform => provably inert) |

### &ocean_epbl_nml

RH18 energetics-based planetary boundary layer.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `enable` | `.false.` |  | Master switch (replaces the KPP overlay; requires use_closure) |
| `mstar_scheme` | `"om4"` |  | Surface TKE mstar scheme |
| `mstar` | `0.1200000000E+01` |  | Constant-scheme mstar |
| `mstar_cap` | `-0.1000000000E+01` |  | Cap for OM4/RH18 mstar; off when < 0 |
| `mstar_coef1` | `0.3000000000E+00` |  | OM4 stabilizing coefficient |
| `c_ek` | `0.8500000000E-01` |  | OM4 Ekman coefficient |
| `mstar_conv_adj` | `0.0000000000E+00` |  | Convective mstar reduction in [0,1] |
| `rh18_cn1` | `0.2750000000E+00` |  | RH18 mstar fit coefficient cn1 |
| `rh18_cn2` | `0.8000000000E+01` |  | RH18 mstar fit coefficient cn2 |
| `rh18_cn3` | `-0.5000000000E+01` |  | RH18 mstar fit coefficient cn3 |
| `rh18_cs1` | `0.2000000000E+00` |  | RH18 mstar fit coefficient cs1 |
| `rh18_cs2` | `0.4000000000E+00` |  | RH18 mstar fit coefficient cs2 |
| `nstar` | `0.2000000000E+00` |  | Convective PE -> TKE efficiency |
| `tke_decay` | `0.2500000000E+01` |  | Ekman-depth / TKE-decay-scale ratio |
| `wstar_ustar_coef` | `0.1000000000E+01` |  | Convective weight in the velocity scale |
| `vel_scale_scheme` | `"cube_root"` |  | Velocity-scale scheme |
| `vstar_scale_fac` | `0.1000000000E+01` |  | Overall vstar multiplier |
| `vstar_surf_fac` | `0.1200000000E+01` |  | RH18 mechanical surface vstar factor |
| `von_karman` | `0.4100000000E+00` |  | von Karman kappa in Kd = vstar*kappa*mixlen |
| `ekman_scale_coef` | `0.1000000000E+01` |  | Rotational mixing-length rolloff |
| `min_mix_len` | `0.0000000000E+00` | m | Mixing-length floor |
| `mixlen_exponent` | `0.2000000000E+01` |  | Shape-function exponent |
| `translay_scale` | `0.1000000000E+00` |  | Transition-layer shape floor (in [0,1) when iterating) |
| `mld_iteration` | `.true.` |  | Self-consistent MLD root-find |
| `mld_tol` | `0.1000000000E+01` | m | MLD convergence tolerance |
| `mld_max_its` | `20` |  | Max MLD iterations |
| `mld_bisection` | `.false.` |  | Bisection instead of false position |
| `mld_use_prev_guess` | `.false.` |  | Seed from the previous step's MLD |
| `omega` | `0.7292100000E-04` | 1/s | Earth rotation rate |
| `omega_frac` | `0.0000000000E+00` |  | Blend |f| with 2*Omega |
| `prandtl` | `0.1000000000E+01` |  | Kv = prandtl*Kd into the momentum solve |
| `combine` | `"add"` |  | Combine vs interior closure kv/kt |
| `tke_diags` | `.false.` |  | Compute per-column TKE budget diagnostics |
| `use_lt` | `.false.` |  | Langmuir-turbulence enhancement (LF17 wind-only) |
| `lt_scheme` | `"rescale"` |  | Langmuir enhancement scheme |
| `lt_enhance_coef` | `0.4470000000E+00` |  | Langmuir enhancement coefficient |
| `lt_enhance_exp` | `-0.1330000000E+01` |  | Langmuir-number exponent |
| `lt_max_enhance` | `0.5000000000E+01` |  | Cap on the multiplicative enhancement |
| `la_frac_hbl` | `0.4000000000E-01` |  | Stokes SL-average depth fraction |
| `lt_lac1` | `-0.8700000000E+00` |  | Stability-modified La coefficient 1 |
| `lt_lac2` | `0.0000000000E+00` |  | Stability-modified La coefficient 2 |
| `lt_lac3` | `0.0000000000E+00` |  | Stability-modified La coefficient 3 |
| `lt_lac4` | `0.9500000000E+00` |  | Stability-modified La coefficient 4 |
| `lt_lac5` | `0.9500000000E+00` |  | Stability-modified La coefficient 5 |

### &ocean_wavespeed_nml

First-baroclinic wave speed + Rossby deformation radius (diagnostic).

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `enable` | `.false.` |  | Master switch (diagnostic; default off) |
| `mono_n2` | `-0.1000000000E+01` |  | DEFERRED N2-monotonising depth (EBT path); < 0 = off |
| `use_ebt` | `.false.` |  | DEFERRED equivalent-barotropic variant |
| `n_wavespeed` | `1` |  | Recompute cadence (every N steps) |

### &ocean_foxkemper_nml

Fox-Kemper mixed-layer-eddy restratification (B5).

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `enable` | `.false.` |  | Master switch (off => bit-identity) |
| `ce` | `0.6250000000E-01` |  | FK08 coefficient Ce (0.06-0.08) |
| `f_floor` | `0.1000000000E-04` | 1/s | |f| regularisation floor |
| `mld_decay_time` | `0.0000000000E+00` | s | Running-mean MLD filter time-scale (0 = off, instantaneous MLD) |
| `tail_dh` | `0.0000000000E+00` |  | mu cubic-tail extension (0 = exact mu) |
| `use_mom_mixrate` | `.false.` |  | FK11 momentum-mixrate timescale (PRODUCTION-RECOMMENDED) vs bare Ce/|f| |
| `resolution_taper` | `.false.` |  | B2 res_fn hook (hard error if on without B2) |
| `use_bodner` | `.false.` |  | Bodner 2023 frontogenesis-arrest MLE (overrides ce/mixrate) |
| `cr` | `0.0000000000E+00` |  | Bodner 2023 efficiency coefficient Cr (0 = off) |
| `bodner_mstar` | `0.5000000000E+00` |  | Bodner mechanical (u*) weight in w'u' |
| `bodner_nstar` | `0.6600000000E-01` |  | Bodner convective (w*) weight in w'u' |
| `min_wstar2` | `0.1000000000E-23` | m^2/s^2 | Floor on w'u' (1/0 armour) |

### &ocean_grid_nml

Horizontal-grid generator + geometry (ocean curvilinear stream).

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `grid_config` | `"cartesian"` |  | Horizontal grid generator |
| `lon_west` | `0.0000000000E+00` | degrees_east | West edge of the domain (spherical) |
| `lat_south` | `0.0000000000E+00` | degrees_north | South edge of the domain (spherical) |
| `rad_earth` | `0.6378000000E+07` | m | Earth radius for the spherical metric |
| `supergrid_file` | `""` |  | Path to the MOM6 supergrid NetCDF (supergrid grid_config) |
| `coriolis_scheme` | `"beta_plane"` |  | Coriolis source for the metrics f-fill |
| `omega` | `0.7292100000E-04` | rad/s | Planetary rotation rate (planetary scheme) |
| `phi_join` | `0.6500000000E+02` | degrees_north | Join latitude for the tripolar bipolar cap |
| `lon_pole` | `0.1000000000E+03` | degrees_east | Longitude of the first tripolar cap pole (partner +180) |
| `axis_units` | `"meters"` |  | Units of the Cartesian domain extent (MOM6 AXIS_UNITS) |
| `len_lon` | `0.0000000000E+00` |  | Total x-extent of the Cartesian domain in axis_units (MOM6 LENLON; derives dx when > 0) |
| `len_lat` | `0.0000000000E+00` |  | Total y-extent of the Cartesian domain in axis_units (MOM6 LENLAT; derives dy when > 0) |

### &ocean_coriolis_nml

Coriolis-advection scheme selector.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `form` | `"sadourny"` |  | Coriolis-advection variant |
| `pv_adv_scheme` | `"centered"` |  | PV face interpolation (Sadourny path): centered (default) or weno3/weno5/weno7 (WENO-Z); weno5/weno7 need nghost>=3/4 |
| `use_state_fluxes` | `.false.` |  | mom6-corrector CorAdv consumes continuity's renormalised mass fluxes (MOM6 mass-consistent uh/vh) |
| `bound_coriolis` | `.false.` |  | clamp the energy-scheme Coriolis accel to the (f+zeta)*v velocity-form range (MOM6 BOUND_CORIOLIS; sadourny_energy only) |
| `corner_h` | `"cell_mean"` |  | PV corner-thickness construction (energy scheme) |

### &ocean_thermo_nml

Thermodynamics master switch + scalar surface fluxes.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `sw_source` | `"net_heat"` |  | Shortwave irradiance source: net_heat (legacy) or q_sw |
| `kpp_sw_method` | `"mxl_sw"` |  | KPP shortwave-in-BL method: all_sw | mxl_sw | lv1_sw |
| `epbl_sw_ctke` | `.true.` |  | Charge the EPBL TKE ledger for penetrating shortwave |
| `enable_thermodynamics` | `.true.` |  | Run EOS + tracer advection + vertical mixing |
| `q_heat` | `0.0000000000E+00` | W/m^2 | Net surface heat flux (positive down) |
| `q_salt` | `0.0000000000E+00` | kg/m^2/s | Net surface salt flux (positive salinifies) |
| `sw_pen_frac` | `0.0000000000E+00` |  | Penetrating fraction of q_heat (0 = off, all at surface) |
| `sw_band_ratio` | `0.5800000000E+00` |  | Two-band shortwave band-1 weight R (Jerlov type I) |
| `sw_zeta1` | `0.3500000000E+00` | m | Shortwave band-1 e-folding depth |
| `sw_zeta2` | `0.2300000000E+02` | m | Shortwave band-2 e-folding depth |

### &ocean_forcing_nml

Surface-flux component-set gate (PR-12): allocate the q_sw/evap/heat_content_*/... set and derive Q_heat/Q_salt from it every thermo step.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `enable_components` | `.false.` |  | Allocate the surface-flux component set and run the assembler (default off => byte-identical) |

### &ocean_ice_nml

Sea-ice model (SIS2 port) master switch + category/layer counts.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `enable` | `.false.` |  | Master switch for the sea-ice slot (default off => byte-identical) |
| `ncat` | `5` |  | Number of ice thickness categories (cat 0 = open water) |
| `hlim` | `-0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01` | m | ITD category lower thickness edges; unset => SIS2 default table |
| `nk_ice` | `2` |  | Vertical ice layers per category (2 = Winton two-layer) |
| `air_temp` | `0.0000000000E+00` | degC | Prescribed slab-atmosphere air temperature for the v1 restoring filler |
| `restore_lambda` | `0.0000000000E+00` | W/m^2/K | Surface-flux restoring coefficient dSF/dT (0 = passive column) |
| `sw_down` | `0.0000000000E+00` | W/m^2 | Downwelling shortwave into the ice top for the v1 filler |
| `snowfall` | `0.0000000000E+00` | kg/m^2/s | Uniform frozen-precipitation rate onto the ice top for the v1 filler |
| `snow_ice` | `.false.` |  | Enable Archimedes snow-ice flooding conversion (SIS2 SN2IC; default off => byte-identical) |
| `transport` | `.false.` |  | Enable horizontal category ice/snow transport (default off => byte-identical) |
| `adv_substeps` | `1` |  | Advective sub-iterations per transport call (SIS2 NSTEPS_ADV) |
| `roll_factor` | `0.1000000000E+01` |  | Thin-ice rolling floor factor (SIS2 SEA_ICE_ROLL_FACTOR); 0 disables rolling |
| `dynamics` | `.false.` |  | Enable C-grid EVP ice dynamics (default off => byte-identical) |
| `a_face_stress` | `.false.` |  | Weight EVP wind stress + ice-ocean drag by face ice concentration (momentum-conserving; default off => byte-identical) |
| `p0` | `0.2750000000E+05` | Pa | Ice-strength pressure constant (SIS2 ICE_STRENGTH_PSTAR) |
| `c0` | `0.2000000000E+02` |  | Ice-strength exponent constant (SIS2 ICE_STRENGTH_CSTAR) |
| `ec` | `0.2000000000E+01` |  | Yield-curve axis ratio (SIS2 ICE_YIELD_ELLIPTICITY); 0 = cavitating fluid |
| `cdw` | `0.3240000000E-02` |  | Ice-ocean drag coefficient (SIS2 ICE_CDRAG_WATER) |
| `rho_ocean` | `0.1030000000E+04` | kg/m^3 | Ice-drag reference density (SIS2 RHO_OCEAN) |
| `evp_sub_steps` | `432` |  | EVP subcycles per slow step (SIS2 NSTEPS_DYN) |
| `del_sh_min_scale` | `0.2000000000E+01` |  | Viscosity-floor scale (SIS2 ICE_DEL_SH_MIN_SCALE) |
| `tdamp` | `-0.2000000000E+00` |  | Elastic damping timescale rule (SIS2 ICE_TDAMP_ELASTIC): >0 seconds, ==0 auto (0.2*dt_slow), <0 fraction of dt_slow |
| `cfl_trunc` | `0.0000000000E+00` |  | Transport-CFL ceiling on the final ice velocity (SIS2 CFL_TRUNCATE, default 0.5 there); 0 disables |
| `cfl_trunc_dyn_its` | `.false.` |  | Also clip the ice velocity every EVP subcycle (SIS2 CFL_TRUNC_DYN_ITS) |
| `project_ci` | `.false.` |  | Project ice concentration forward within the EVP subcycle loop and recompute the ice strength (SIS2 PROJECT_ICE_CONCENTRATION) |

### &ocean_ice_ic_nml

Sea-ice analytic initial-condition path (PR 24).

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `conc_config` | `"zero"` |  | Sea-ice initial concentration configuration |
| `conc` | `0.0000000000E+00` |  | Uniform-mode concentration ('uniform' only) |
| `h_ice` | `0.0000000000E+00` | m | Ice thickness where seeded (SIS2 ICE_INIT_MASS as a thickness) |
| `h_snow` | `0.0000000000E+00` | m | Snow thickness where seeded (SIS2 SNOW_INIT_MASS as a thickness) |
| `t_ice` | `-0.4000000000E+01` | degC | Ice/snow temperature fed through the exact enthalpy inversion (SIS2 ICE_TEMPERATURE_IC) |
| `s_ice` | `0.4000000000E+01` | PSU | Ice bulk salinity (SIS2 ICE_SALINITY_IC) |
| `arctic_edge` | `0.9100000000E+02` | degrees_north | 'latitudes' Arctic ice edge (SIS2 ARCTIC_ICE_EDGE_IC) |
| `antarctic_edge` | `-0.9100000000E+02` | degrees_north | 'latitudes' Antarctic ice edge (SIS2 ANTARCTIC_ICE_EDGE_IC) |

### &ocean_restore_nml

Surface buoyancy restoring (MOM6 RESTOREBUOY): piston-velocity relaxation of top-layer T / S toward targets.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `enable_restore_temp` | `.false.` |  | Master switch for SST restoring |
| `enable_restore_salt` | `.false.` |  | Master switch for SSS restoring |
| `piston_t` | `0.0000000000E+00` | m/day | SST piston velocity (MOM6 FLUXCONST_T) |
| `piston_s` | `0.0000000000E+00` | m/day | SSS piston velocity (MOM6 FLUXCONST_S) |
| `restore_sst` | `0.0000000000E+00` | degC | Scalar target SST |
| `restore_sss` | `0.0000000000E+00` | PSU | Scalar target SSS |

### &ocean_geothermal_nml

Geothermal bottom heat flux (bed-side analogue of the surface heat flux).

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `enable` | `.false.` |  | Apply a constant geothermal bottom heat flux to the bed layer |
| `q_geo` | `0.0000000000E+00` | W/m^2 | Constant bottom heat flux (positive into the ocean from below) |

### &ocean_sponge_nml

Map-driven sponge: per-cell Idamp + 3-D reference state (PR-23).

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `enable` | `.false.` |  | Master switch (default off; legacy band kernels run when off) |
| `damp_source` | `"band"` |  | How idamp_h/u/v are filled |
| `target_source` | `"ic"` |  | Reference-state source |
| `relax_uv` | `.true.` |  | Relax u/v toward u_ref/v_ref |
| `relax_tracers` | `.true.` |  | Relax every registered tracer toward its 3-D reference field |
| `relax_h` | `.false.` |  | Interior-interface thickness damping (NOT IMPLEMENTED in v1 — deferred to PR-23b, requires vcoord_type='lagrangian') |
| `west_width` | `-1` |  | Per-edge band-width override, cells (<0 => inherit &ocean_bc_nml sponge_width) |
| `east_width` | `-1` |  | as west_width |
| `south_width` | `-1` |  | as west_width |
| `north_width` | `-1` |  | as west_width |
| `west_strength` | `-0.1000000000E+01` | 1/s | Per-edge peak relaxation-rate override (<0 => inherit &ocean_bc_nml sponge_strength) |
| `east_strength` | `-0.1000000000E+01` | 1/s | as west_strength |
| `south_strength` | `-0.1000000000E+01` | 1/s | as west_strength |
| `north_strength` | `-0.1000000000E+01` | 1/s | as west_strength |

### &ocean_tracers_nml

Prognostic-tracer registry switches.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `enable_ideal_age` | `.false.` |  | Register a passive ideal-age tracer |
| `ideal_age_young_val` | `0.0000000000E+00` | s | Surface-band ideal-age Dirichlet value (0 = today's hard-coded reset) |
| `ideal_age_sfc_growth_rate` | `0.0000000000E+00` | 1/s | Exponential growth rate of the ideal-age surface value (0 = constant) |
| `enable_pseudo_salt` | `.false.` |  | Register the pseudo-salt verification tracer (diagnostic; seeded to S, given salinity's boundary fluxes) |

### &ocean_bt_nml

Split-explicit barotropic substep controls.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `n_inner` | `0` |  | Barotropic substeps per outer step (0 = unsplit) |
| `auto_n_inner` | `.false.` |  | Derive n_inner from the gravity-wave CFL at setup |
| `cfl_bt_safety` | `0.6500000000E+00` |  | Safety fraction on the BT CFL when auto_n_inner |
| `bebt` | `0.0000000000E+00` |  | Forward-velocity-projection weight (MOM6 BEBT) |
| `use_cont_type` | `.false.` |  | Use the BT_cont flux-bounded closure (MOM6 USE_BT_CONT_TYPE) |
| `cont_corr_bounds` | `.false.` |  | Use BT_cont flux limits for the eta-correction bound |
| `upstream_h_face` | `.false.` |  | Use per-face upstream-PPM column-sum thickness in the BT chain |
| `correction_h_weighted` | `.false.` |  | h-weight the post-substep BT corrector (MOM6 frhatu) |
| `correction_visc_rem` | `.false.` |  | h*visc_rem joint corrector weight (requires correction_h_weighted; visc_rem is produced by vdiff and is inert, =1, without ocean_vdiff_nml implicit_drag) |
| `split_scheme` | `"pred_corr"` |  | Outer split-explicit time scheme: pred_corr (DEFAULT; MOM6 predictor-corrector, slow tendencies on the u_av/h_av step time-means, forward-backward gravity-wave pairing; lifts the internal-wave dt ceiling) or ssp_rk2 (EXPERIMENTAL; two-stage SSP average, widest envelope — the only scheme wired through eulerian_z, wet/dry and dt_tracer_advect_ratio>1 — but it spuriously grows internal gravity waves out of a stratified REST state, En 2.992E-05 vs 1.739E-09 at day 25 on resting_stratified_channel.nml; a (omega*dt)^4 noise floor, so forced viscous runs sit decades above it and quiescent or long spin-up runs do not) |
| `pc_be` | `0.6000000000E+00` |  | pred_corr predictor fraction BE (MOM6 BE, 0.6 reference) |
| `renorm_visc_rem` | `.false.` |  | gamma-weighted continuity transport-matching inversion (MOM6 u_cor = u + du*visc_rem; requires correction_visc_rem) |
| `forcing_visc_rem` | `.false.` |  | MOM6 wt_u parity: h*visc_rem-weight the BT forcing depth-mean so friction-damped (glued) layers do not force the fast loop (requires correction_visc_rem) |
| `correction_bc_pgf` | `.false.` |  | Per-layer baroclinic-PGF retro-correction for the eta change |
| `substep_drag` | `.false.` |  | Apply a per-substep BT velocity damping factor |
| `substep_zeta_ke` | `.true.` |  | Integrate live zeta_bt + KE-gradient in the BT fast loop (.false. = MOM6 parity: planetary Coriolis only, zeta/KE frozen in the slow forcing) |
| `wave_drag` | `.false.` |  | Barotropic linear wave drag master switch (MOM6 BT_LINEAR_WAVE_DRAG) |
| `wave_drag_form` | `"uniform"` |  | Wave-drag r_H filler |
| `wave_drag_scale` | `0.1000000000E+01` |  | Global tuning multiplier on r_H (MOM6 BT_WAVE_DRAG_SCALE) |
| `wave_drag_r_uniform` | `0.0000000000E+00` | m/s | Piston velocity r_H for wave_drag_form='uniform' |
| `wave_drag_kappa` | `0.6283200000E-03` | 1/m | Topographic wavenumber for wave_drag_form='roughness_proxy' |
| `wave_drag_n_bot` | `0.1000000000E-02` | 1/s | Reference bottom N for wave_drag_form='roughness_proxy' |
| `wave_drag_h2_max` | `0.2500000000E+05` | m^2 | Ceiling on <h^2> proxy for wave_drag_form='roughness_proxy' |
| `wave_drag_file` | `""` |  | Reserved for PR-14 (MOM6 BT_WAVE_DRAG_FILE); unused today |
| `wave_drag_var` | `"rH"` |  | Reserved for PR-14 (MOM6 BT_WAVE_DRAG_VAR); unused today |
| `bt_halo` | `-1` |  | Wide-halo BT march-in width (-1 = auto: 8 under multi-rank if compatible, else 0; 0 = explicit off, bit-identical) |

### &ocean_debug_nml

Forensic probes for the ocean dyn-core (all heavy, all default off).

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `budget` | `.false.` |  | Per-stage BT power-budget probe (BUDGET rows) |
| `ke_attr` | `.false.` |  | Per-segment layer-KE attribution meter (KE_ATTR rows) |
| `ke_attr_start_step` | `0` |  | First outer step the KE meter samples (0 = from start) |
| `ke_attr_end_step` | `0` |  | Last outer step the KE meter samples (0 = unbounded) |
| `chksum` | `.false.` |  | MOM6-style per-phase field checksums (CHKSUM rows) + HOTFACE argmax-face anatomy at the tendency seams |
| `chksum_start_step` | `0` |  | First outer step the chksum probe samples (0 = from start) |
| `chksum_end_step` | `0` |  | Last outer step the chksum probe samples (0 = unbounded) |
| `chksum_interior` | `.false.` |  | Reduce over physical cells only, making `bits` a valid 1-rank-vs-N-rank gate |

### &ocean_mpi_nml

Multi-rank MPI debug/tuning controls.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `poison_ghosts` | `.false.` |  | Sentinel-NaN the exchange-covered ghost bands at each outer-step start; unexchanged-ghost consumption becomes a loud NaN. Default off = bit-identical. |

### &ocean_wetdry_nml

Dynamic wetting/drying for the split-explicit barotropic substep.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `enable` | `.false.` |  | Dynamic wet/dry: upwind BT face thickness + positive-definite outflow limiter + bed-blocking momentum gate |
| `dry_depth` | `0.5000000000E-01` | m | Total-depth dry threshold |
| `rewet_depth` | `0.1000000000E+00` | m | Hysteresis re-wet threshold (must exceed dry_depth) |
| `land_margin` | `0.5000000000E+01` | m | Static-land headroom above rest MSL (intertidal cutoff) |

### &ocean_pgf_nml

Pressure-gradient-force kernel selector + knobs.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `form` | `"mont"` |  | PGF kernel variant |
| `gprime_gfs` | `0.9810000000E+01` | m/s^2 | Free-surface gravity for the gprime PGF |
| `gprime_gint` | `0.9800000000E-02` | m/s^2 | Internal-interface reduced gravity (gprime PGF) |
| `gfs_scale` | `0.1000000000E+01` |  | Free-surface gravity scaling (FV_MOM6, MOM6 GFS_scale) |
| `maxvel` | `0.0000000000E+00` | m/s | Velocity-truncation clamp (0 = disabled) |
| `cfl_trunc` | `0.0000000000E+00` |  | Advective-CFL velocity truncation threshold (0 = disabled) |
| `mass_weight` | `.false.` |  | FV_MOM6 shelf-break hWght mass-weighting at unequal-depth faces |
| `reconstruct_for_pressure` | `.false.` |  | FV_MOM6 in-layer PLM/PPM T/S reconstruction for the density integral |
| `recon_scheme` | `1` |  | In-layer reconstruction scheme: 1=PLM, 2=PPM |

### &ocean_eos_nml

Equation-of-state variant selector, liquidus set + reference pressure.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `eos` | `"linear"` |  | Equation-of-state variant |
| `tfreeze_set` | `"seaice"` |  | Named liquidus coefficient set for eos_freezing_point (T_f = l1*S + l2 + l3*p): 'seaice' = SIS2/MOM6 (-0.054, 0, -7.53e-8), 'isomip' = ISOMIP+ (-0.0573, 0.0832, -7.53e-8) |
| `p_ref` | `0.0000000000E+00` | Pa | Reference pressure for the potential density ms%rho_layer (horizontally uniform by design; 0 => surface density) |

### &ocean_bdrag_nml

Bottom-drag selector + coefficients.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `form` | `"quadratic"` |  | Bottom-drag variant |
| `cd` | `0.0000000000E+00` |  | Quadratic drag coefficient (0 disables) |
| `r` | `0.0000000000E+00` | 1/s | Linear Rayleigh coefficient (0 disables) |
| `hbbl` | `0.0000000000E+00` | m | BBL thickness for distributed drag (0 = bed-only) |
| `bg_vel` | `0.0000000000E+00` | m/s | Background velocity floor for distributed drag |
| `bbl_thick_min` | `0.0000000000E+00` | m | Minimum effective BBL thickness |
| `bed_factor` | `0.1000000000E+01` |  | Multiplier on the bed-layer drag tendency only |
| `channel_drag` | `.false.` |  | Enable per-layer lateral side-wall (channel) Rayleigh drag |
| `cdrag_side` | `0.0000000000E+00` |  | Side-wall drag coefficient (0 disables channel drag) |
| `implicit` | `.false.` |  | Backward-Euler implicit drag (stable for thin bottom layers) |

### &ocean_hdiff_nml

Along-coordinate (not neutral) constant-coefficient tracer diffusion.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `kappa_h` | `0.0000000000E+00` | m^2/s | Horizontal tracer diffusivity, along the model coordinate (0 disables; ocean path only — coastal's equivalent knob is top-level hdiff_kappa) |

### &ocean_hvisc_nml

Lateral-viscosity closure + coefficients.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `lateral_closure` | `"none"` |  | Lateral-mixing closure tag |
| `c_smag` | `0.1500000000E+00` |  | Smagorinsky coefficient |
| `c_leith` | `0.1000000000E+01` |  | Leith coefficient |
| `kh_vel_scale` | `0.0000000000E+00` | m/s | Velocity scale for the resolution viscosity floor |
| `ah_bg` | `-0.1000000000E+01` | m^2/s | Background harmonic viscosity (<0 = derive) |
| `ah_max` | `0.1000000000E+05` | m^2/s | Cap on per-face harmonic viscosity |
| `smag_ah` | `.false.` |  | Flow-aware biharmonic viscosity (MOM6 SMAGORINSKY_AH) |
| `smag_bi_const` | `0.6000000000E-01` |  | Biharmonic Smagorinsky constant (MOM6 SMAG_BI_CONST) |
| `c_leith_bi` | `0.0000000000E+00` |  | Biharmonic Leith constant (MOM6 LEITH_BI_CONST) |
| `nu_4_bg` | `0.0000000000E+00` | m^4/s | Background biharmonic viscosity floor |
| `nu_4_max` | `0.1000000000E+13` | m^4/s | Cap on per-face biharmonic viscosity |
| `nu_h` | `0.0000000000E+00` | m^2/s | Constant horizontal eddy viscosity |
| `nu_4` | `0.0000000000E+00` | m^4/s | Constant biharmonic eddy viscosity |
| `no_slip` | `.false.` |  | Coastal lateral BC: .false.=free-slip (×wet_q), .true.=no-slip (×(2-wet_q)) |
| `stress_tensor` | `.false.` |  | MOM6 thickness-weighted stress-div operator + per-cell CFL + coast-mask |
| `bound_coef` | `0.8000000000E+00` |  | Per-cell viscosity-CFL safety coefficient (MOM6 HORVISC_BOUND_COEF) |
| `bound_kh` | `.false.` |  | Per-face harmonic viscosity CFL clamp on the velocity-Laplacian paths (MOM6 BOUND_KH) |
| `resoln_scaled_visc` | `.false.` |  | Scale dynamic LAPLACIAN viscosity (Leith/Smag_KH A_h, not biharmonic) by VarMix Res_fn (MOM6 RESOLN_SCALED_KH) |
| `kh_vel_scale_live` | `0.0000000000E+00` | m/s | Live velocity-scale viscosity Kh=U*dx*|u| (0=off; distinct from kh_vel_scale floor) |
| `kh_aniso` | `0.0000000000E+00` | m^2/s | Anisotropic Laplacian viscosity magnitude (MOM6 KH_ANISO; stress_tensor path only) |
| `aniso_mode` | `0` |  | Anisotropy direction mode (0=grid-relative aniso_dir; MOM6 ANISOTROPIC_MODE) |
| `aniso_dir` | `0.1000000000E+01, 0.0000000000E+00` |  | Anisotropy direction (n1,n2) in grid i,j components (MOM6 ANISO_GRID_DIR) |

### &ocean_vmix_nml

Vertical-mixing module switches + knobs.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `use_closure` | `.true.` |  | Master switch for the interior closure (PP81) |
| `use_kpp` | `.true.` |  | KPP surface-boundary-layer overlay (needs use_closure) |
| `direct_stress` | `.false.` |  | Spread wind stress over hmix_stress (MOM6 DIRECT_STRESS) |
| `hmix_stress` | `0.2000000000E+02` | m | Surface-slab thickness for direct_stress |
| `kv_ml_invz2` | `0.0000000000E+00` | m^2/s | Near-surface 1/(z hmix)^2 viscosity (MOM6 KV_ML_INVZ2) |
| `hmix_fixed` | `0.2000000000E+02` | m | Mixed-layer thickness for the KV_ML_INVZ2 profile |
| `harmonic_visc` | `.false.` |  | Harmonic-mean face thickness in vdiff (MOM6 HARMONIC_VISC) |
| `dt_therm_ratio` | `1` |  | Thermo step runs at ratio*dt_dyn (MOM6 DT_THERM) |
| `dt_tracer_advect_ratio` | `1` |  | Horizontal tracer advect runs at ratio*dt_dyn (MOM6 DT_TRACER_ADVECT) |
| `tracer_recon` | `"ppm"` |  | Windowed tracer-advect drain face reconstruction (ocean): ppm|weno5|weno7|weno9 |
| `kv_max` | `0.1797693135+309` | m^2/s | Assembly ceiling on kv; huge=off (MOM6 Kd_max momentum) |
| `kd_max` | `0.1797693135+309` | m^2/s | Assembly ceiling on kt/ks; huge=off (MOM6 Kd_max) |
| `kd_smooth_iterations` | `0` |  | 1-2-1 smoothing passes on kv/kt (MOM6 Kd_smooth) |
| `vmix_guard` | `.false.` |  | Debug-gated negative/NaN diffusivity guard (assembly) |
| `bkgnd_profile` | `.false.` |  | Bryan-Lewis depth-varying background diffusivity |
| `bkgnd_kd_sfc` | `0.1000000000E-04` | m^2/s | Bryan-Lewis surface-asymptote background Kd |
| `bkgnd_kd_deep` | `0.1300000000E-03` | m^2/s | Bryan-Lewis deep-asymptote background Kd |
| `bkgnd_z0` | `0.2500000000E+04` | m | Bryan-Lewis transition-centre depth |
| `bkgnd_delta` | `0.2220000000E+03` | m | Bryan-Lewis transition half-width |
| `bkgnd_prandtl` | `0.1000000000E+01` |  | Background Prandtl number Kv_bg=prandtl*Kd_bg |
| `bkgnd_henyey` | `.false.` |  | Henyey IGW latitude factor on the scalar background (excludes bkgnd_profile; needs a non-cartesian grid) |
| `bkgnd_kd_min` | `-0.1000000000E+01` | m^2/s | Minimum background Kd under the Henyey scaling (negative = 0.01*kt_bg) |
| `bkgnd_henyey_n0_2omega` | `0.2000000000E+02` |  | Henyey N0/(2*Omega) reference stratification ratio |
| `bkgnd_henyey_max_lat` | `0.9500000000E+02` | degN | Latitude poleward of which the Henyey factor floors |
| `pp81_nu0` | `0.1000000000E-01` | m^2/s | PP81 Richardson-dependent viscosity scale |
| `pp81_nu_bg` | `0.1000000000E-03` | m^2/s | PP81 background viscosity (also seeds vmix%kv_bg) |
| `pp81_kappa_bg` | `0.1000000000E-04` | m^2/s | PP81 background diffusivity (also seeds vmix%kt_bg/ks_bg) |
| `pp81_alpha` | `0.5000000000E+01` |  | PP81 Richardson-number scaling coefficient (paper value 5) |
| `shear2_floor` | `0.1000000000E-09` | 1/s^2 | Floor on |du/dz|^2+|dv/dz|^2 in the PP81 Ri denominator |
| `kpp_ri_crit` | `0.3000000000E+00` |  | Critical bulk Richardson number for KPP BL-depth |
| `kpp_cs_nonlocal` | `0.6300000000E+01` |  | KPP non-local (counter-gradient) transport coefficient C_s |
| `kpp_c_vt2` | `0.1800000000E+01` |  | KPP unresolved-turbulence V_t^2 coefficient (0 disables V_t^2) |

### &ocean_vdiff_nml

Backward-Euler vertical-friction solver knobs.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `implicit_stress` | `.false.` |  | Fold wind stress into the vdiff surface (k=nz) RHS |
| `accel_visc_rem` | `.false.` |  | MOM6 parity: attenuate the slow explicit accelerations by the per-layer viscous remnant (u = u0 + visc_rem*(u-u0) after the applies); requires ocean_bt_nml correction_visc_rem |
| `implicit_drag` | `.false.` |  | Fold bottom drag into the vdiff bed (k=1) diagonal |
| `hvel_mom6` | `.false.` |  | MOM6 HARMONIC_VISC parity: harmonic momentum face thickness with the near-bed upwind blend, and arithmetic h_shear. Suppresses grounded-sliver momentum as MOM6 does |
| `hbbl_visc` | `0.1000000000E+02` |  | Bottom-layer scale for the hvel_mom6 botfn blend (MOM6 HBBL) |
| `bbl_glue` | `.false.` |  | MOM6 bottomdraglaw coupling parity: kv_bbl botfn glue at near-bed interfaces + piston bed drag. Absorbs the spurious grounded-layer PGF as MOM6 does (PGF_BUG.md par.9). Requires hvel_mom6 + implicit_drag + linear bottom drag |
| `bbl_piston` | `0.3000000000E-03` | m/s | BBL drag piston velocity u* for bbl_glue (MOM6 CDRAG*DRAG_BG_VEL) |
| `hvel_upwind` | `.true.` |  | Near-bed upwind blend in the hvel_mom6 thickness build (.false. = pure harmonic; the u-sign blend flip-flops on roundoff at rest and collapses the BBL glue, PGF_BUG.md par.9.8) |

### &ocean_continuity_nml

Continuity-PPM positivity controls.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `h_min` | `0.1000000000E-05` | m | PPM positivity-limiter thickness floor |
| `ppm_limit_pos` | `.false.` |  | Positivity-preserving PPM face limiter (MOM6 PPM_limit_pos) |
| `vol_cfl` | `.false.` |  | Swept-volume continuity-PPM face thickness (MOM6 vol_CFL) |
| `positive_definite` | `.false.` |  | Positive-definite split continuity (h>=h_lim, zero mass created; reconstruction floor + per-donor outflux limiter) |

### &ocean_isopycnal_nml

Lagrangian grounding-stability controls.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `angstrom_h` | `0.0000000000E+00` | m | Minimum-thickness floor on the Lagrangian continuity h-update (MOM6 Angstrom_H analogue) |
| `reset_vanished_u` | `.false.` |  | Zero face velocity when layer vanished on both adjacent cells |
| `cfl_ignore_vanished` | `.false.` |  | Exclude vanished layers from MaxCFL / panic / CFL truncation |
| `pgf_skip_nonoverlap` | `.true.` |  | Zero the face PGF where a grounded layer's z-extents do not overlap across the face (VCOORD_LAGRANGIAN only, mont / fv_lite / fv_wright / fv_mom6; default ON — kills the spurious at-rest grounded-layer pressure gradient) |
| `conservative_floor` | `.false.` |  | Conservative min-thickness borrow (replaces the injecting angstrom_h floor) |
| `check_h_positive` | `.false.` |  | DEBUG: abort on the first negative h_layer, naming the stage + (i,j,k) |

### &ocean_topo_nml

Basin geometry + surface forcing + Coriolis tilt.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `topo_config` | `"flat"` |  | Bathymetry profile selector |
| `max_depth` | `0.2000000000E+04` | m | Basin maximum depth |
| `edge_depth` | `0.1000000000E+03` | m | Spoon-bathymetry edge depth |
| `slope_scale` | `0.4000000000E+06` | m | Spoon-bathymetry exponential decay scale |
| `nl_continent_amp` | `0.1000000000E+01` | nondim | Neverworld2 continent amplitude (1=full continents, 0=aquaplanet+channel) |
| `nl_roughness_amp` | `0.5000000000E-01` | nondim | Neverworld2 bathymetry roughness amplitude |
| `nl_min_depth` | `0.5000000000E+03` | m | Neverworld2 minimum-depth floor (MOM6 MINIMUM_DEPTH analogue) |
| `wind_config` | `"constant"` |  | Surface wind-stress dispatch |
| `taux_magnitude` | `0.1000000000E+00` | Pa | Peak zonal wind stress for wind_config=2gyre/neverworld2 |
| `coriolis_beta` | `0.0000000000E+00` | 1/(s m) | Meridional gradient of f (0 = f-plane) |
| `coriolis_y_ref` | `0.0000000000E+00` | m | Reference y where f = coriolis_f under beta-plane |

### &ocean_ic_nml

Initial-condition overlay + EOS reference state.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `ic_config` | `""` |  | IC overlay tag |
| `alpha_T` | `0.1700000000E-03` | kg/m^3/degC | Linear-EOS thermal-expansion coefficient (DIMENSIONAL: multiply a fractional 1/degC coefficient by rho_0) |
| `beta_S` | `0.7600000000E-03` | kg/m^3/PSU | Linear-EOS haline contraction coefficient (DIMENSIONAL: multiply a fractional 1/PSU coefficient by rho_0) |
| `T_ref` | `0.1000000000E+02` | degC | Linear-EOS reference temperature |
| `S_ref` | `0.3500000000E+02` | PSU | Linear-EOS reference salinity |
| `rho_0` | `0.1035000000E+04` | kg/m^3 | Reference density for the linear EOS / Boussinesq PGF |
| `layer_rho_init` | `-0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01` | kg/m^3 | Per-layer initial density (k=1 bed -> k=nz surface; -1 = EOS init) |
| `rho_lightest` | `-0.1000000000E+01` | kg/m^3 | Linear density-range IC: surface (lightest) layer density; -1 = off |
| `rho_range` | `0.2000000000E+01` | kg/m^3 | Linear density-range IC: top-to-bottom density contrast (MOM6 DENSITY_RANGE) |
| `eady_dT_dy` | `-0.2000000000E-04` | degC/m | Eady IC meridional T gradient |
| `eady_dT_dz` | `0.1000000000E-01` | degC/m | Eady IC vertical stratification |
| `eady_T_ref` | `0.1000000000E+02` | degC | Eady IC reference temperature at z=0,y=y_mid |
| `eady_pert_amp` | `0.1000000000E-02` | degC | Eady IC symmetry-breaking perturbation amplitude |
| `eady_pert_seed` | `12345` |  | RNG seed for the Eady IC perturbation |
| `ga_eta_amp` | `0.1000000000E+01` | m | Geostrophic-adjustment IC: SSH-bump amplitude |
| `ga_length_scale` | `0.5000000000E+05` | m | Geostrophic-adjustment IC: SSH-bump e-folding scale |
| `ga_x_center` | `-0.1000000000E+01` | m | Geostrophic-adjustment IC: bump x-centre (<0 = auto) |
| `ga_y_center` | `-0.1000000000E+01` | m | Geostrophic-adjustment IC: bump y-centre (<0 = auto) |
| `jet_half_width` | `0.4000000000E+05` | m | Baroclinic-jet IC: tanh jet half-width L |
| `interface_amp` | `0.2000000000E+03` | m | Baroclinic-jet IC: interface displacement amplitude |
| `pert_amp_frac` | `0.2000000000E+00` |  | Baroclinic-jet IC: meander amplitude as a fraction of interface_amp |
| `pert_nx` | `3` |  | Baroclinic-jet IC: zonal perturbation wavenumber (integer) |
| `upper_layer_rest` | `0.5000000000E+03` | m | Baroclinic-jet IC: upper (surface) layer rest thickness H1 |

### &ocean_zinit_nml

Z-level T/S initial-condition overlay (A2).

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `enable` | `.false.` |  | Master switch (default off; requires RDB_ENABLE_NETCDF=ON) |
| `file` | `""` |  | Path to the model-grid T/S NetCDF |
| `t_var` | `""` |  | Temperature variable-name override (blank tries temp/T/temperature) |
| `s_var` | `""` |  | Salinity variable-name override (blank tries salt/S/salinity) |
| `z_var` | `""` |  | Source-axis variable-name override (blank tries z_src/z/depth/lev) |
| `land_fill_t` | `0.1000000000E+02` | degC | Fallback temperature for dry columns |
| `land_fill_s` | `0.3500000000E+02` | PSU | Fallback salinity for dry columns |

### &ocean_data_nml

Shared time-varying NetCDF input reader (PR-14): registry sizing only.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `max_fields` | `16` |  | Size of the reader's field registry |
| `verbose` | `.false.` |  | Log every bracket advance (field, records, weight, time) |

### &ocean_dataovr_nml

File-backed surface forcing (PR-15): wind stress, heat, freshwater.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `enable` | `.false.` |  | Master switch (default off; requires RDB_ENABLE_NETCDF=ON) |
| `time_mode` | `"linear"` |  | Shared time-axis mode for every active tag |
| `cycle_period` | `0.0000000000E+00` | s | Climatology period (required > 0 when time_mode=cyclic) |
| `t_offset` | `0.0000000000E+00` | s | Added to the model time before the file-axis lookup |
| `oor_clamp` | `.false.` |  | Clamp an out-of-range query to the end record instead of aborting |
| `tau_x_file` | `""` |  | Path to the zonal wind stress NetCDF (blank ⇒ tag not file-driven) |
| `tau_x_var` | `""` |  | Variable name inside tau_x_file (required when it is set) |
| `tau_x_scale` | `0.1000000000E+01` | Pa | Multiplier applied to the zonal wind stress slab at read |
| `tau_x_add` | `0.0000000000E+00` | Pa | Offset added to the zonal wind stress slab after scale |
| `tau_y_file` | `""` |  | Path to the meridional wind stress NetCDF (blank ⇒ tag not file-driven) |
| `tau_y_var` | `""` |  | Variable name inside tau_y_file (required when it is set) |
| `tau_y_scale` | `0.1000000000E+01` | Pa | Multiplier applied to the meridional wind stress slab at read |
| `tau_y_add` | `0.0000000000E+00` | Pa | Offset added to the meridional wind stress slab after scale |
| `heat_file` | `""` |  | Path to the net surface heat flux (positive down) NetCDF (blank ⇒ tag not file-driven) |
| `heat_var` | `""` |  | Variable name inside heat_file (required when it is set) |
| `heat_scale` | `0.1000000000E+01` | W m-2 | Multiplier applied to the net surface heat flux (positive down) slab at read |
| `heat_add` | `0.0000000000E+00` | W m-2 | Offset added to the net surface heat flux (positive down) slab after scale |
| `evap_file` | `""` |  | Path to the evaporative mass flux (<= 0) NetCDF (blank ⇒ tag not file-driven) |
| `evap_var` | `""` |  | Variable name inside evap_file (required when it is set) |
| `evap_scale` | `0.1000000000E+01` | kg m-2 s-1 | Multiplier applied to the evaporative mass flux (<= 0) slab at read |
| `evap_add` | `0.0000000000E+00` | kg m-2 s-1 | Offset added to the evaporative mass flux (<= 0) slab after scale |
| `lprec_file` | `""` |  | Path to the liquid precipitation (>= 0) NetCDF (blank ⇒ tag not file-driven) |
| `lprec_var` | `""` |  | Variable name inside lprec_file (required when it is set) |
| `lprec_scale` | `0.1000000000E+01` | kg m-2 s-1 | Multiplier applied to the liquid precipitation (>= 0) slab at read |
| `lprec_add` | `0.0000000000E+00` | kg m-2 s-1 | Offset added to the liquid precipitation (>= 0) slab after scale |
| `salt_file` | `""` |  | Path to the surface salt flux (positive salinifies) NetCDF (blank ⇒ tag not file-driven) |
| `salt_var` | `""` |  | Variable name inside salt_file (required when it is set) |
| `salt_scale` | `0.1000000000E+01` | kg m-2 s-1 | Multiplier applied to the surface salt flux (positive salinifies) slab at read |
| `salt_add` | `0.0000000000E+00` | kg m-2 s-1 | Offset added to the surface salt flux (positive salinifies) slab after scale |

### &ocean_diag_nml

Ocean diag-manager output controls.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `enabled` | `.true.` |  | Enable the per-step diag-manager hook |
| `filename` | `"ocean_diag"` |  | Output file basename (per-rank suffix appended) |
| `dt_out` | `0.3600000000E+04` |  | Diag-fire cadence (in time_unit) |
| `vgrid` | `"layer"` |  | Default output vertical grid for layer-shaped diags |
| `z_levels` | `-0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01` | m | Output z-levels when vgrid=z_fixed (positive down) |
| `n_z_levels` | `0` |  | Number of z_levels entries used (0 = none) |
| `sigma_levels` | `-0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01` |  | Output sigma fractions (0..1) for sigma output; empty = auto-uniform |
| `n_sigma_levels` | `0` |  | Number of sigma_levels entries used (0 = auto) |
| `zstar_levels` | `-0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01` | m | Output z* reference depths for zstar output; empty = auto-uniform |
| `n_zstar_levels` | `0` |  | Number of zstar_levels entries used (0 = auto) |
| `rho_levels` | `-0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01, -0.1000000000E+01` | kg/m^3 | Output potential-density bin edges when vgrid=density; strictly increasing, light->dense |
| `n_rho_levels` | `0` |  | Number of rho_levels entries used (0 = none) |
| `mask_vanished_layers` | `.false.` |  | Mask below-bottom/pinched remap cells to missing_value (default off) |
| `reproducing_sums` | `.false.` |  | Order-invariant EFP console totals: rank-count reproducible, exact drift residual (default off = byte-identical) |
| `output_precision` | `"double"` |  | Element width of the diag NetCDF data vars; 'single' halves the bytes written (restarts/gauges/console totals stay double regardless) |
| `diag_remap_scheme` | `"ppm"` |  | Reconstruction for the conservative diagnostic vertical remap |
| `diags` | `""` |  | Unified diag list modifying the default set: name[:off][:cadence][:op] (e.g. 'vorticity_z:1d  KE:off  temperature:6h:mean'); default empty = canonical set |

### &ocean_bc_nml

Open-boundary condition config: per-edge BC type, clamped/inflow Dirichlet values, sponge band, Orlanski radiation, full Flather, per-edge tidal constituents, tidal-OBC nodal correction, wall-velocity masking.

| Knob | Default | Units | Description |
|------|---------|-------|-------------|
| `west` | `"wall"` |  | West edge boundary type |
| `east` | `"wall"` |  | East edge boundary type |
| `south` | `"wall"` |  | South edge boundary type |
| `north` | `"wall"` |  | North edge boundary type ('tripolar_fold' is meaningful here only — Murray bipolar cap; requires grid_config='tripolar' + periodic west/east) |
| `west_clamped_eta` | `0.0000000000E+00` | m | Clamped SSH, west edge |
| `east_clamped_eta` | `0.0000000000E+00` | m | Clamped SSH, east edge |
| `south_clamped_eta` | `0.0000000000E+00` | m | Clamped SSH, south edge |
| `north_clamped_eta` | `0.0000000000E+00` | m | Clamped SSH, north edge |
| `west_clamped_u` | `0.0000000000E+00` | m/s | Clamped normal velocity, west edge |
| `east_clamped_u` | `0.0000000000E+00` | m/s | Clamped normal velocity, east edge |
| `south_clamped_v` | `0.0000000000E+00` | m/s | Clamped normal velocity, south edge |
| `north_clamped_v` | `0.0000000000E+00` | m/s | Clamped normal velocity, north edge |
| `west_inflow_S` | `0.3500000000E+02` | PSU | Inflow salinity, west edge |
| `west_inflow_T` | `0.1000000000E+02` | degC | Inflow temperature, west edge |
| `east_inflow_S` | `0.3500000000E+02` | PSU | Inflow salinity, east edge |
| `east_inflow_T` | `0.1000000000E+02` | degC | Inflow temperature, east edge |
| `south_inflow_S` | `0.3500000000E+02` | PSU | Inflow salinity, south edge |
| `south_inflow_T` | `0.1000000000E+02` | degC | Inflow temperature, south edge |
| `north_inflow_S` | `0.3500000000E+02` | PSU | Inflow salinity, north edge |
| `north_inflow_T` | `0.1000000000E+02` | degC | Inflow temperature, north edge |
| `sponge_width` | `0` |  | Sponge band width in cells (0 disables) |
| `sponge_strength` | `0.0000000000E+00` | 1/s | Peak sponge relaxation rate at the outer face of the band |
| `sponge_relax_tracers` | `.false.` |  | Extend the sponge relaxation to h_layer + tracers |
| `res_lscale_out` | `0.0000000000E+00` | m | Outflow reservoir length scale (0 = disabled) |
| `res_lscale_in` | `0.0000000000E+00` | m | Inflow reservoir length scale (0 = instantaneous inflow) |
| `radiation_scheme` | `"anomaly"` |  | 'anomaly' = v1 BT-mean + zero-gradient anomaly (default); 'orlanski' = per-layer implicit-upwind radiation (Orlanski 1976) |
| `orlanski_rx_max` | `0.1000000000E+02` | grid cells/step | Clamp on the Orlanski nondimensional phase speed |
| `orlanski_gamma` | `0.1000000000E+01` |  | Running-mean weight (0=full running mean, 1=instant) |
| `nudge_tau_in` | `0.0000000000E+00` | s | Inflow nudging timescale (0 = off) |
| `nudge_tau_out` | `0.0000000000E+00` | s | Outflow nudging timescale (0 = off) |
| `flather_form` | `"legacy"` |  | 'legacy' = v1 Flather (u_ext=0, no interior vel, default); 'full' = half-characteristic form (Flather 1976) |
| `west_ext_u` | `0.0000000000E+00` | m/s | Exterior barotropic u, west edge |
| `east_ext_u` | `0.0000000000E+00` | m/s | Exterior barotropic u, east edge |
| `south_ext_v` | `0.0000000000E+00` | m/s | Exterior barotropic v, south edge |
| `north_ext_v` | `0.0000000000E+00` | m/s | Exterior barotropic v, north edge |
| `west_n_tidal` | `0` |  | Active tidal constituents, west edge |
| `west_tidal_amp` | `0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00` | m | Constituent amplitudes, west edge |
| `west_tidal_phase` | `0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00` | rad | Constituent phases, west edge |
| `west_tidal_omega` | `0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00` | rad/s | Constituent angular frequencies, west edge |
| `east_n_tidal` | `0` |  | Active tidal constituents, east edge |
| `east_tidal_amp` | `0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00` | m | Constituent amplitudes, east edge |
| `east_tidal_phase` | `0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00` | rad | Constituent phases, east edge |
| `east_tidal_omega` | `0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00` | rad/s | Constituent angular frequencies, east edge |
| `south_n_tidal` | `0` |  | Active tidal constituents, south edge |
| `south_tidal_amp` | `0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00` | m | Constituent amplitudes, south edge |
| `south_tidal_phase` | `0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00` | rad | Constituent phases, south edge |
| `south_tidal_omega` | `0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00` | rad/s | Constituent angular frequencies, south edge |
| `north_n_tidal` | `0` |  | Active tidal constituents, north edge |
| `north_tidal_amp` | `0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00` | m | Constituent amplitudes, north edge |
| `north_tidal_phase` | `0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00` | rad | Constituent phases, north edge |
| `north_tidal_omega` | `0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00, 0.0000000000E+00` | rad/s | Constituent angular frequencies, north edge |
| `obc_tidal_nodal` | `.false.` |  | Apply the 18.6-yr nodal factor + equilibrium/nodal phase to the open-boundary tidal elevation forcing (shares the &ocean_tides_nml astro generator); when true *_tidal_phase becomes a Greenwich phase lag |
| `mask_wall_velocity` | `.true.` |  | Zero the T-cell wet-mask in the ghost cells beyond every solid WALL edge at setup, so wall-normal C-grid face velocities are masked rather than left as garbage (default ON; set .false. to reproduce a pre-fix legacy closed-basin baseline) |
