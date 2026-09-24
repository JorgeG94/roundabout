!! Runtime configuration via Fortran namelists
module rdb_config
   !! Reads simulation parameters from a namelist input file
#ifdef LFORTRAN_PASSING
   use rdb_constants, only: wp, DEG2RAD, &
                            nz_stack_required, nz_stack_is_sufficient
#else
   use rdb_constants, only: NZ_STACK_MAX, wp, DEG2RAD, &
                            nz_stack_required, nz_stack_is_sufficient
#endif
   use pic_logger, only: logger => global_logger
   use rdb_error_ring, only: fail, error_ring_push
   use pic_strings, only: to_string
   use pic_ascii, only: to_lower
   use rdb_nml_schema, only: nml_schema_t, nml_group_t, &
                             nml_real, nml_int, nml_logical, nml_string, &
                             nml_enum, nml_real_array, nml_real_array_key_t
   use rdb_ice_enthalpy, only: ice_t_freeze
   use rdb_ice_init, only: ice_ic_parse_conc_config, ICE_IC_CONC_INVALID
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_CONFIG_PARSE, &
                               OCEAN_STATUS_ERR_CONFIG_VALIDATE
   implicit none
   private

#ifdef LFORTRAN_PASSING
   integer, parameter :: NZ_STACK_MAX = 64
      !! LFortran 0.64 workaround: module-local copy of the rdb_constants value
      !! (an imported parameter used as an explicit-shape dummy bound inside a
      !! PURE call becomes an impure getter under LFortran). Keep in sync (=64).
#endif

   ! ------------------------------------------------------------------
   ! Retired coastal-legacy `&tracer_nml` linear-EOS quartet.
   !
   ! These four keys fed `tracer_t%eos_coeff` / `tracer_t%eos_ref` via
   ! `register_default_tracers` — a path the A-grid coastal solvers read
   ! and the C-grid ocean path never did.  They are now RETIRED rather
   ! than merely annotated dead: `validate_config` fails loud whenever a
   ! namelist moves one off its historical default, naming the live
   ! `&ocean_ic_nml` replacement.  The defaults are named here so the
   ! type declaration and the guard can never drift apart.
   ! ------------------------------------------------------------------
   real(wp), parameter :: LEGACY_TRACER_S_REF = 0.0_wp
      !! Historical `&tracer_nml S_ref` default (PSU).
   real(wp), parameter :: LEGACY_TRACER_BETA_S = 0.78_wp
      !! Historical `&tracer_nml beta_S` default (kg/m^3 per PSU).
   real(wp), parameter :: LEGACY_TRACER_T_REF = 15.0_wp
      !! Historical `&tracer_nml T_ref` default (degC).
   real(wp), parameter :: LEGACY_TRACER_ALPHA_T = 0.17_wp
      !! Historical `&tracer_nml alpha_T` default (kg/m^3 per degC).

   public :: config_t
   public :: ocean_config_t, ocean_grid_config_t, &
             ocean_coriolis_config_t, ocean_thermo_config_t, &
             ocean_ice_config_t, ocean_ice_ic_config_t, &
             ocean_restore_config_t, &
             ocean_geothermal_config_t, &
             ocean_sponge_config_t, &
             ocean_tracers_config_t, &
             ocean_bt_config_t, ocean_pgf_config_t, ocean_eos_config_t, &
             ocean_bdrag_config_t, &
             ocean_tdrag_config_t, &
             ocean_hdiff_config_t, &
             ocean_hvisc_config_t, ocean_vmix_config_t, &
             ocean_vdiff_config_t, &
             ocean_epbl_config_t, ocean_wave_speed_config_t, &
             ocean_foxkemper_config_t, &
             ocean_kappa_shear_config_t, &
             ocean_gm_config_t, ocean_varmix_config_t, &
             ocean_meke_config_t, &
             ocean_tidal_mixing_config_t, ocean_conv_config_t, ocean_tides_config_t, &
             ocean_porous_config_t, &
             ocean_psurf_config_t, &
             ocean_cavity_dyn_config_t, ocean_cavity_melt_config_t, &
             ocean_continuity_config_t, ocean_isopycnal_config_t, &
             ocean_topo_config_t, &
             ocean_ic_config_t, ocean_zinit_config_t, &
             ocean_data_config_t, &
             ocean_dataovr_config_t, dataovr_entry_config_t, &
             ocean_diag_config_t, ocean_bc_config_t, &
             ocean_wetdry_config_t, ocean_mpi_config_t
   public :: read_config
   public :: read_config_from_string
   public :: validate_config
   public :: diag_density_levels_ok
   public :: dataovr_entry_is_valid, dataovr_time_is_valid
   public :: dataovr_freshwater_needs_components, dataovr_any_tag_set
   public :: build_rdb_schema
   public :: resolve_bt_halo
   public :: bt_halo_auto_exclusion
   public :: p_top_has_producer
   public :: substep_drag_ignores_bdrag_form
   public :: cavity_draft_is_uniform
   public :: MAX_TIDAL_CONSTITUENTS
   public :: MAX_OCEAN_DIAG_Z_LEVELS
   public :: MAX_OCEAN_LAYER_RHO_INIT
   public :: MAX_ICE_HLIM_VALS
   public :: BT_HALO_AUTO_SENTINEL
   public :: BT_HALO_AUTO_WIDTH
   public :: ice_hlim_count
   public :: ice_hlim_spec_is_valid

   integer, parameter :: BT_HALO_AUTO_SENTINEL = -1
      !! `&ocean_bt_nml bt_halo` default: AUTO.  Resolved at configure (in the
      !! driver, where `compute_size` is known) to `BT_HALO_AUTO_WIDTH` under a
      !! compatible multi-rank run, else 0.  See `resolve_bt_halo`.
   integer, parameter :: BT_HALO_AUTO_WIDTH = 8
      !! Wide-halo BT march-in width chosen when `bt_halo` auto-resolves ON
      !! (the validated production width).

   integer, parameter :: MAX_TIDAL_CONSTITUENTS = 10
      !! Maximum number of tidal constituents
   integer, parameter :: MAX_OCEAN_DIAG_Z_LEVELS = 64
      !! Maximum number of z-levels for ocean-diag z_fixed output vgrid
   integer, parameter :: MAX_OCEAN_LAYER_RHO_INIT = 64
      !! Maximum number of per-layer density init entries.
   integer, parameter :: MAX_ICE_HLIM_VALS = 16
      !! Maximum number of `&ocean_ice_nml hlim` entries a user may STATE
      !! (PR-58). Caps how many edges may be listed, not `ncat` (unbounded
      !! elsewhere) — production ITDs are 5-10 categories, so 16 is a
      !! generous soft cap; exceeding it is not possible for `ncat<=15`,
      !! and for larger `ncat` the remaining edges simply extrapolate
      !! (`ice_hlim_spec_is_valid`).  Deliberately NOT matched to
      !! `MAX_OCEAN_LAYER_RHO_INIT` (64) — that cap is dimensioned by
      !! `nz`, this one by `ncat+1`.

   ! Nested ocean config — one sub-type per concern, nested under
   ! `config_t%ocean`.

   type :: ocean_grid_config_t
      !! `&ocean_grid_nml`: horizontal-grid generator + geometry knobs.
      character(len=16) :: grid_config = "cartesian"
         !! Grid generator: "cartesian" (default, uniform), "spherical"
         !! (lon-lat sector), "supergrid" (MOM6 mosaic NetCDF), "tripolar".
      real(wp) :: lon_west = 0.0_wp
         !! West edge of the physical domain (degrees) — spherical only.
      real(wp) :: lat_south = 0.0_wp
         !! South edge of the physical domain (degrees) — spherical only.
      real(wp) :: rad_earth = 6.378e6_wp
         !! Earth radius for the spherical metric (m).
      character(len=256) :: supergrid_file = ""
         !! Path to the MOM6 supergrid NetCDF — "supergrid" only.
      character(len=16) :: coriolis_scheme = "beta_plane"
         !! Coriolis source for the metrics f-fill: "beta_plane" (default,
         !! f_0 + beta*y) or "planetary" (2*omega*sin(geolat)).
      real(wp) :: omega = 7.2921e-5_wp
         !! Planetary rotation rate (rad/s) — used by "planetary".
      real(wp) :: phi_join = 65.0_wp
         !! Join latitude (degrees) for "tripolar": lon-lat below, bipolar
         !! Arctic cap above.
      real(wp) :: lon_pole = 100.0_wp
         !! Longitude (degrees) of the first tripolar cap pole; partner at
         !! `lon_pole + 180`.
      character(len=16) :: axis_units = "meters"
         !! Units of the Cartesian domain extent (`grid_config="cartesian"`),
         !! mirroring MOM6 AXIS_UNITS: "meters" (default), "degrees", "km".
         !! With "degrees"/"km" the domain is sized by `len_lon`/`len_lat`
         !! below and `&grid_nml dx`/`dy` are DERIVED (ignored as inputs).
         !! "meters" ⇒ `&grid_nml dx`/`dy` used verbatim (bit-identical).
      real(wp) :: len_lon = 0.0_wp
         !! Total x-extent of the Cartesian domain in `axis_units` (MOM6
         !! LENLON).  <= 0 (default) ⇒ unset.  > 0 ⇒ a uniform
         !! `dx = to_metres(len_lon)/nx` is derived, where for "degrees"
         !! `to_metres(v) = rad_earth·v·π/180` (arc length, no cos(lat) —
         !! MOM6 cartesian) and for "km" `to_metres(v) = 1000·v`.
      real(wp) :: len_lat = 0.0_wp
         !! Total y-extent of the Cartesian domain in `axis_units` (MOM6
         !! LENLAT).  Paired with `len_lon`; > 0 ⇒ `dy = to_metres(len_lat)/ny`.
      !! NOTE: spherical cell sizes reuse `&grid_nml dx`/`dy` as dlon/dlat
      !! in DEGREES when `grid_config="spherical"`.
   end type ocean_grid_config_t

   type :: ocean_coriolis_config_t
      !! `&ocean_coriolis_nml`.
      character(len=16) :: form = "sadourny"
         !! Coriolis-advection variant: "sadourny" (default, enstrophy form),
         !! "sadourny_energy" (energy-conserving transport form q·vh),
         !! "sadourny_hk" (Hollingsworth-Källén).
      character(len=16) :: pv_adv_scheme = "centered"
         !! PV face-interpolation scheme, orthogonal to `form` (Sadourny path
         !! only): "centered" (default, 2-point corner average ⇒ bit-identical)
         !! or "weno3"/"weno5"/"weno7" — upwind-biased WENO-Z reconstruction of
         !! the corner absolute vorticity onto the velocity faces (MOM6
         !! WENOVI{3,5,7}TH), sharpening submesoscale PV fronts without the
         !! centred form's global dissipation.  weno5/weno7 (radius-3/4
         !! stencils) require `nghost >= 3`/`4` (fail-loud otherwise).
      logical :: use_state_fluxes = .false.
         !! MOM6 mass-consistent CorAdCalc: the mom6 corrector's
         !! Coriolis-adv consumes the predictor continuity's renormalised
         !! transports (`ms%mass_flux_*_layer` — the uh/vh from the same
         !! solve that produced `u_av`) instead of recomputing `u·h_face`
         !! from the evaluation-state prognostics.  Requires
         !! `form="sadourny_energy"` + `&ocean_bt_nml split_scheme="pred_corr"`
         !! (fail-loud).  Default `.false.` = bit-identical.
      character(len=16) :: corner_h = "cell_mean"
         !! PV corner-thickness construction for the energy scheme.
         !! `"cell_mean"` (default, bit-identical) — Roundabout's wet-area-weighted
         !! 4-cell mean `h_corner = Σ(wet·area·h)/Σ(wet·area)`, floored at
         !! `H_MIN_PV`.  `"mom6_area"` — MOM6's exact form:
         !! `q = abs_vort·Area_q/(hArea_q +
         !! vol_neglect)` with `hArea_q = Σ(mask·area·h)`, `Area_q =
         !! Σ(mask·area)`.  **The two are algebraically IDENTICAL above the
         !! floor** (same numerator + denominator); they differ ONLY in the
         !! vanishing-thickness guard — Roundabout CAPS `q` at `abs_vort/H_MIN_PV`,
         !! MOM6's `vol_neglect` is pure 1/0 armor (no cap), so at
         !! sub-`H_MIN_PV` corners MOM6's `q` is LARGER.  For realistic layers
         !! (`h ≳ 1e-12`) the knob is a round-off no-op.  Energy scheme only
         !! (fail-loud otherwise).
      logical :: bound_coriolis = .false.
         !! MOM6 `BOUND_CORIOLIS`: clamp the
         !! energy-scheme Coriolis acceleration into the range of the four
         !! neighbouring VELOCITY-form estimates `(f+ζ)·v` BEFORE the KE-gradient
         !! subtraction.  Kills the thin-layer PV blow-up: `q=(f+ζ)/h_corner`
         !! is huge at a vanishing edge and `q·vh` can inject an acceleration
         !! the local velocity field cannot support; the `(f+ζ)·v` bound caps
         !! it.  Energy scheme only (`form="sadourny_energy"`, fail-loud
         !! otherwise) — MOM6's double-gyre control run enables it.  Default
         !! `.false.` ⇒ untaken branch ⇒ bit-identical.
   end type ocean_coriolis_config_t

   type :: ocean_thermo_config_t
      !! Thermodynamics master switches + scalar surface fluxes.
      logical :: enable_thermodynamics = .true.
         !! When `.true.` (default), the dyn step runs EOS + tracer
         !! advection + vertical mixing.  `.false.` = adiabatic.
      real(wp) :: q_heat = 0.0_wp
         !! Scalar net surface heat flux (W/m^2, positive down).
      real(wp) :: q_salt = 0.0_wp
         !! Scalar net surface salt flux (kg salt/m^2/s, positive =
         !! surface salinifies).
      real(wp) :: sw_pen_frac = 0.0_wp
         !! Penetrating fraction of `q_heat` as a two-band exponential
         !! (Paulson & Simpson 1977).  0 = off.
      real(wp) :: sw_band_ratio = 0.58_wp
         !! Band-1 weight `R` of the two-band shortwave decay (Jerlov type I).
      real(wp) :: sw_zeta1 = 0.35_wp
         !! Band-1 (red/near-IR) e-folding depth (m).
      real(wp) :: sw_zeta2 = 23.0_wp
         !! Band-2 (blue/green) e-folding depth (m).
      character(len=32) :: sw_source = "net_heat"
         !! (PR-21) Irradiance source for shortwave penetration + the
         !! boundary-layer SW coupling.  `"net_heat"` (default,
         !! bit-identical): `I0 = sw_pen_frac*q_heat`.  `"q_sw"`:
         !! `I0 = sw_pen_frac*q_sw` (the dedicated >= 0 shortwave
         !! component) — kills the night-time negative-`I0` hazard;
         !! REQUIRES `&ocean_forcing_nml enable_components=.true.`
         !! (validate_config aborts otherwise).
      character(len=32) :: kpp_sw_method = "mxl_sw"
         !! (PR-21) KPP boundary-layer shortwave method (MOM6
         !! `KPP_SHORTWAVE_METHOD`): `all_sw` | `mxl_sw` (default) |
         !! `lv1_sw`.  Only the SW absorbed inside the boundary layer
         !! stabilises `B_0`.  Inert at `sw_pen_frac = 0` ⇒ bit-identical.
      logical :: epbl_sw_ctke = .true.
         !! (PR-21) Charge the EPBL TKE ledger for penetrating shortwave
         !! (the in-layer PE-cost `Phi(tau)` shape).  Inert at
         !! `sw_pen_frac = 0` ⇒ bit-identical.
   end type ocean_thermo_config_t

   type :: ocean_forcing_config_t
      !! `&ocean_forcing_nml`: PR-12 surface-flux component-set gate.
      !! Default OFF ⇒ bit-identical (no component array allocated, no
      !! assembler kernel runs; `ocean_surface_flux_t%Q_heat`/`Q_salt`
      !! are filled exactly as `&ocean_thermo_nml` always has).
      logical :: enable_components = .false.
         !! Allocate the surface-flux component set (`q_sw`, `evap`,
         !! `heat_content_*`, ...) on `ocean_surface_flux_t` and run
         !! `ocean_surface_flux_assemble` every thermo step to derive
         !! `Q_heat`/`Q_salt` from it.  Individual components carry no
         !! knobs of their own — they are filled by fillers (readers,
         !! the ice coupler, ...), not by this namelist.
   end type ocean_forcing_config_t

   type :: ocean_ice_config_t
      !! `&ocean_ice_nml`: sea-ice model master switch + category/layer
      !! counts (SIS2 port, `PLAN_SEA_ICE.md`).  PR 0 scaffold — the knobs
      !! size the gated `ocean_sea_ice_t` slot.  PR 3c adds the v1 scalar
      !! restoring atmospheric-forcing filler (`air_temp`/`restore_lambda`/
      !! `sw_down`) that drives the column live in the driver.  Default OFF
      !! ⇒ byte-identical.
      logical :: enable = .false.
         !! Master switch.  Off (default) ⇒ the sea-ice slot is never
         !! initialised, mapped, or stepped ⇒ byte-identical.
      integer :: ncat = 5
         !! Number of ice thickness categories (category 0 = open water;
         !! SIS2 default 5).
      real(wp) :: hlim(MAX_ICE_HLIM_VALS) = -1.0_wp
         !! PR-58: ITD category lower thickness edges (m), overriding the
         !! hardcoded SIS2 default table. `-1.0` sentinel (all entries) =
         !! unset.  Semantics when set: `hlim(1..n)` where `n` is the
         !! leading non-sentinel run, `2 <= n <= ncat+1`, strictly
         !! increasing, `hlim(1) > 0`.  Fills `h_lim(1..min(ncat+1,n))`;
         !! the remainder extrapolates by constant width,
         !! `h_lim(k) = 2*h_lim(k-1) - h_lim(k-2)`.  Unset => the
         !! hardcoded SIS2 default table, byte-identical (default-off
         !! house rule).  See `ice_hlim_count`/`ice_hlim_spec_is_valid`.
      integer :: nk_ice = 2
         !! Vertical ice layers per category (2 = Winton two-layer, v1).
      real(wp) :: air_temp = 0.0_wp
         !! Prescribed slab-atmosphere air temperature (degC) for the v1
         !! restoring atmospheric-forcing filler. `SF(T) = restore_lambda*(T - air_temp)`.
      real(wp) :: restore_lambda = 0.0_wp
         !! Surface-flux restoring coefficient (W/m^2/K), i.e. dSF/dT of the
         !! linearized SEB. 0 (default) => a passive column (dsf_dt=0, sf_0=0):
         !! no thermostat, so no restoring drive. A production melt season uses
         !! ~20. NOTE: with restore_lambda=0 the seam fields are (0,0,sw_down)
         !! => the column feels only shortwave + basal flux.
      real(wp) :: sw_down = 0.0_wp
         !! Downwelling shortwave into the ice top (W/m^2) for the v1 filler.
      real(wp) :: snowfall = 0.0_wp
         !! Uniform frozen-precipitation rate onto the ice top (kg/m^2/s)
         !! for the v1 restoring atmospheric-forcing filler (PR 26). Spread
         !! onto the `atm_fprec` seam exactly like `sw_down` -> `atm_sw_dn`.
         !! 0 (default) => byte-identical (`has_snowfall=.false.`, the new
         !! kernels never fire).

      ! ---- PR-27: snow-ice flooding ----
      logical :: snow_ice = .false.
         !! Archimedes freeboard snow-ice conversion (SIS2
         !! `ice_resize_SIS2` final block; diagnosed as SN2IC). When the
         !! snow load submerges the snow-ice interface, convert the
         !! submerged snow mass to ice in the top layer. Default OFF =>
         !! byte-identical.

      ! ---- PR 4b: category ice transport + compress_ice ----
      logical :: transport = .false.
         !! Enable horizontal category ice/snow transport
         !! (`rdb_ice_transport%ice_transport_step`) — SIS2's default
         !! (velocity, non-merged) `ice_cat_transport` path. Default OFF ⇒
         !! byte-identical for every existing nml/test (house rule).
      integer :: adv_substeps = 1
         !! Advective sub-iterations per transport call (SIS2 `NSTEPS_ADV`,
         !! default 1). `dt_adv = dt_therm / adv_substeps`.
      real(wp) :: roll_factor = 1.0_wp
         !! SIS2 `SEA_ICE_ROLL_FACTOR` — thin-ice rolling floor in
         !! `cell_ave_state_to_ice_state`. 0 disables rolling.

      ! ---- PR 5: C-grid EVP ice dynamics ----
      logical :: dynamics = .false.
         !! EVP rheology master switch (`rdb_ice_evp%ice_evp_step`).
         !! Default off => byte-identical (the transport sampler keeps
         !! filling `u_ice`/`v_ice` from the ocean surface layer).
      logical :: a_face_stress = .false.
         !! PR 62: weight the atmospheric stress AND the ice-ocean drag in
         !! the EVP momentum balance by the face ice concentration
         !! `a_u = 0.5*(ci(i-1,j) + ci(i,j))` (same expression/edge
         !! convention as `ice_ocean_stress_flux_impl`'s `a_u`), giving the
         !! textbook `m du/dt = grad.sigma + a*(tau_a - tau_w)` (Hibler 1979
         !! eq. 1) and an EXACTLY closing ice<->ocean momentum budget at
         !! every fractional cover, not just steady free drift. Both terms
         !! must be weighted together -- weighting the wind alone converts
         !! today's leak (zero at steady free drift) into a permanent one.
         !! Default off => byte-identical.
      real(wp) :: p0 = 2.75e4_wp
         !! SIS2 `ICE_STRENGTH_PSTAR` — ice-strength pressure constant [Pa].
      real(wp) :: c0 = 20.0_wp
         !! SIS2 `ICE_STRENGTH_CSTAR` — ice-strength exponent constant [nondim].
      real(wp) :: ec = 2.0_wp
         !! SIS2 `ICE_YIELD_ELLIPTICITY` — yield-curve axis ratio [nondim].
         !! 0 => cavitating-fluid rheology (str_t/str_s stay exactly 0).
      real(wp) :: cdw = 3.24e-3_wp
         !! SIS2 `ICE_CDRAG_WATER` — ice-ocean drag coefficient [nondim].
      real(wp) :: rho_ocean = 1030.0_wp
         !! SIS2 `RHO_OCEAN` — ice-drag reference density [kg/m^3].
         !! Deliberately independent of the ocean's `rho0` (usually 1035).
      integer :: evp_sub_steps = 432
         !! SIS2 `NSTEPS_DYN` — EVP subcycles per slow (outer) step.
      real(wp) :: del_sh_min_scale = 2.0_wp
         !! SIS2 `ICE_DEL_SH_MIN_SCALE` — viscosity-floor scale [nondim].
      real(wp) :: tdamp = -0.2_wp
         !! SIS2 `ICE_TDAMP_ELASTIC`. `>0` => seconds; `==0` =>
         !! `max(0.2*dt_slow, 3*dt)`; `<0` => `max(-tdamp*dt_slow, 3*dt)`.

      ! ---- PR 36: EVP velocity CFL truncation + PROJECT_ICE_CONCENTRATION ----
      real(wp) :: cfl_trunc = 0.0_wp
         !! SIS2 `CFL_TRUNCATE` (SIS2 default **0.5**) — transport-CFL
         !! ceiling on the FINAL ice velocity the EVP call hands to
         !! transport: `|u| <= 0.95*cfl_trunc*areaT(donor)/(dt_transport*
         !! dy_cu)`. `0` (Roundabout default, house bit-identity rule) disables
         !! the clip entirely => byte-identical. Recommended production
         !! value 0.5 ("instability can occur past 0.5", SIS2).
      logical :: cfl_trunc_dyn_its = .false.
         !! SIS2 `CFL_TRUNC_DYN_ITS` (SIS2 default `.false.`, matches).
         !! Also clip `u_ice`/`v_ice` to the EXACT (no 0.95 back-off, no
         !! count) bound at the bottom of every EVP subcycle, not just the
         !! final velocity. Requires `cfl_trunc > 0`.
      logical :: project_ci = .false.
         !! SIS2 `PROJECT_ICE_CONCENTRATION` (SIS2 default **`.true.`**).
         !! Project the ice concentration forward along the current
         !! divergence each subcycle (`ci_proj = ci*exp(-t_cum*sh_dd)`) and
         !! recompute `pres_mice` from it, stiffening the rheology under
         !! convergence (and weakening it under divergence) within the
         !! call. `.false.` (Roundabout default, house bit-identity rule) holds
         !! `pres_mice` at its pre-loop value for all subcycles => byte-identical.
   end type ocean_ice_config_t

   type :: ocean_ice_ic_config_t
      !! `&ocean_ice_ic_nml`: sea-ice ANALYTIC initial-condition path (PR
      !! 24). Default `conc_config="zero"` early-returns before touching
      !! a single ice array => byte-identical to today for every existing
      !! nml/test. v1 ships `"uniform"` (scalar) and `"latitudes"` (SIS2
      !! polar-cap analytic form) only — file-backed ICs are deliberately
      !! NOT in the allowed list (`"file"` is v1.1, PR-14).
      character(len=32) :: conc_config = "zero"
         !! `"zero"` (default, no-op) / `"uniform"` / `"latitudes"`.
      real(wp) :: conc = 0.0_wp
         !! Uniform-mode concentration [0,1] (nondim). `"uniform"` only.
      real(wp) :: h_ice = 0.0_wp
         !! Ice thickness (m) where seeded — SIS2 `ICE_INIT_MASS` as a
         !! thickness.
      real(wp) :: h_snow = 0.0_wp
         !! Snow thickness (m) where seeded — SIS2 `SNOW_INIT_MASS` as a
         !! thickness. Snow is a one-way street (no snowfall source
         !! anywhere in the tree) — this can only ever melt.
      real(wp) :: t_ice = -4.0_wp
         !! Ice/snow temperature (degC) fed through the exact
         !! `ice_enth_from_ts` inversion — SIS2 `ICE_TEMPERATURE_IC`
         !! default.
      real(wp) :: s_ice = 4.0_wp
         !! Ice bulk salinity (PSU) — SIS2 `ICE_SALINITY_IC` default,
         !! matching `ICE_BULK_SALINITY` (`ocean_sea_ice_init`'s own
         !! `sal_ice` source value).
      real(wp) :: arctic_edge = 91.0_wp
         !! `"latitudes"` Arctic ice edge (degrees_north) — SIS2
         !! `ARCTIC_ICE_EDGE_IC` default (no cell qualifies => no ice).
      real(wp) :: antarctic_edge = -91.0_wp
         !! `"latitudes"` Antarctic ice edge (degrees_north) — SIS2
         !! `ANTARCTIC_ICE_EDGE_IC` default (no cell qualifies => no ice).
   end type ocean_ice_ic_config_t

   type :: ocean_restore_config_t
      !! `&ocean_restore_nml`: surface buoyancy restoring of top-layer T/S
      !! toward scalar targets via a piston velocity.  Default OFF = bit-identical.
      logical :: enable_restore_temp = .false.
         !! Master switch for SST restoring (effective when `piston_t /= 0`).
      logical :: enable_restore_salt = .false.
         !! Master switch for SSS restoring (`piston_s /= 0`).
      real(wp) :: piston_t = 0.0_wp
         !! SST piston velocity (m/day).
      real(wp) :: piston_s = 0.0_wp
         !! SSS piston velocity (m/day).
      real(wp) :: restore_sst = 0.0_wp
         !! Scalar target SST (degC).
      real(wp) :: restore_sss = 0.0_wp
         !! Scalar target SSS (PSU).
   end type ocean_restore_config_t

   type :: ocean_geothermal_config_t
      !! `&ocean_geothermal_nml`: bottom-heat-flux knobs.  Default OFF = bit-identical.
      logical :: enable = .false.
         !! Master switch.  Default `.false.`.
      real(wp) :: q_geo = 0.0_wp
         !! Scalar constant bottom heat flux (W/m^2, positive up).  Typical ~0.05-0.1.
   end type ocean_geothermal_config_t

   type :: ocean_sponge_config_t
      !! `&ocean_sponge_nml`: the map-driven sponge (PR-23). Default OFF
      !! (`enable = .false.`) ⇒ the legacy `&ocean_bc_nml` band kernels run
      !! unchanged ⇒ bit-identical. See `rdb_ocean_sponge.F90`'s module
      !! docstring for the physics; `docs/plans/PLAN_PR23_real_sponge.md`
      !! for the design record.
      logical :: enable = .false.
         !! Master switch. Off ⇒ legacy band kernels ⇒ byte-identical.
      character(len=16) :: damp_source = "band"
         !! How `idamp_h`/`idamp_u`/`idamp_v` are filled. `"band"` (default,
         !! implemented): cosine ramp from every `OBC_SPONGE`-tagged edge,
         !! summed at overlaps. `"file"` recognised but aborts at
         !! `validate_config` in v1 (PR-23b, needs the PR-14 reader).
      character(len=16) :: target_source = "ic"
         !! Reference-state source.
         !!
         !! `"ic"` (default, implemented): snapshot the seeded initial
         !! condition (reachable as a "nudge toward a parent climatology"
         !! via `&ocean_zinit_nml`, no new reader).
         !!
         !! `"linear_z"`: an ANALYTIC affine geopotential profile,
         !! `T(z) = lin_t_ref + lin_dt_dz*z` and the salinity twin,
         !! re-evaluated on the LIVE layer geometry once per outer step.
         !! Independent of the initial condition, which is what makes an
         !! ISOMIP+ Ocean1 / Ocean2 (restore to a different water mass
         !! than you start from) expressible. Every tracer that is NOT
         !! temperature or salinity keeps the `"ic"` snapshot, and so do
         !! `u_ref`/`v_ref`.
         !!
         !! `"file"` recognised but aborts at `validate_config` (PR-23b).
      character(len=16) :: ramp = "cosine"
         !! Shape of the `damp_source="band"` ramp from the sponge-tagged
         !! wall (`d = 0`) inward.
         !!
         !! `"cosine"` (default, bit-identical): `0.5*(1 + cos(pi*d/band))`
         !! — the legacy band kernel's own ramp.
         !!
         !! `"linear"`: `(band - d - 0.5)/band`, which is the CELL-CENTRE
         !! evaluation of a rate that rises linearly from zero at the
         !! interior edge of the band to full `sponge_strength` at the
         !! wall. That is ISOMIP+ Eq. (20),
         !! `gamma(x) = gamma0*max(0, (x - x_r0)/(x_r1 - x_r0))`
         !! (Asay-Davis et al. 2016), with `gamma0 = sponge_strength` and
         !! `band = (x_r1 - x_r0)/dx` cells — no separate `tau_boundary`
         !! or x-range knob is needed, because the existing strength is
         !! already `1/tau` in 1/s and the existing width already names
         !! the range.
      real(wp) :: lin_t_ref = 0.0_wp
         !! `target_source="linear_z"`: potential temperature (degC) at
         !! the `z = 0` datum. Same convention as `&ocean_zinit_nml
         !! lin_t_ref`.
      real(wp) :: lin_dt_dz = 0.0_wp
         !! `target_source="linear_z"`: dT/dz (degC/m) with **z positive
         !! UP**, so a stable column has `lin_dt_dz > 0`. Same convention
         !! as `&ocean_zinit_nml lin_dt_dz`.
      real(wp) :: lin_s_ref = 35.0_wp
         !! `target_source="linear_z"`: salinity (PSU) at the `z = 0`
         !! datum. Same convention as `&ocean_zinit_nml lin_s_ref`.
      real(wp) :: lin_ds_dz = 0.0_wp
         !! `target_source="linear_z"`: dS/dz (PSU/m), z positive UP, so a
         !! stable column has `lin_ds_dz < 0`. Same convention as
         !! `&ocean_zinit_nml lin_ds_dz`.
      logical :: relax_uv = .true.
         !! Relax `u`/`v` toward `u_ref`/`v_ref`. Default `.true.` matches
         !! today's legacy path (momentum is the one thing the legacy
         !! sponge always damps); MOM6's `SPONGE_UV` defaults `.false.`.
      logical :: relax_tracers = .true.
         !! Relax every registered tracer toward its 3-D reference field.
      logical :: relax_h = .false.
         !! Interior-interface thickness damping. NOT IMPLEMENTED in PR-23
         !! v1 — `enable=.true., relax_h=.true.` aborts fail-loud at
         !! `validate_config` (deferred to PR-23b; see the plan §14 Q1).
      integer :: west_width = -1
         !! Per-edge sponge-band width override (cells). `< 0` ⇒ inherit
         !! `&ocean_bc_nml sponge_width`.
      integer :: east_width = -1
      integer :: south_width = -1
      integer :: north_width = -1
      real(wp) :: west_strength = -1.0_wp
         !! Per-edge peak relaxation-rate override (1/s). `< 0` ⇒ inherit
         !! `&ocean_bc_nml sponge_strength`.
      real(wp) :: east_strength = -1.0_wp
      real(wp) :: south_strength = -1.0_wp
      real(wp) :: north_strength = -1.0_wp
   end type ocean_sponge_config_t

   type :: ocean_tracers_config_t
      !! Prognostic-tracer registry switches.
      logical :: enable_ideal_age = .false.
         !! When `.true.`, register a passive "ideal age" tracer.  Default off.
      real(wp) :: ideal_age_young_val = 0.0_wp
         !! Surface-band age value (s).  Default 0 = today's hard-coded
         !! reset ⇒ bit-identical.  MOM6 CS%young_val.
      real(wp) :: ideal_age_sfc_growth_rate = 0.0_wp
         !! Exponential growth rate of the surface value (1/s).  0 (default)
         !! ⇒ young_val is constant ⇒ no exp() ⇒ bit-identical.  MOM6
         !! CS%growth_rate (which is per-year; 1/30 yr^-1 = 1.057e-9 s^-1).
      logical :: enable_pseudo_salt = .false.
         !! When `.true.`, register the pseudo-salt verification tracer
         !! (diagnostic; seeded to S, given salinity's boundary fluxes;
         !! typically run to measure passive-vs-active transport-path
         !! error).  Default off.  Incompatible with SSS piston restoring
         !! and sea-ice (both are un-mirrored salinity sources — see
         !! `validate_config`).
   end type ocean_tracers_config_t

   type :: ocean_bt_config_t
      integer :: n_inner = 0
         !! Mode-split barotropic substeps per outer step.  `0` (default)
         !! routes the unsplit `ocean_dyn_step`; `>= 1` routes the split
         !! solver with that many fast substeps (production ~30-100).
         !! Overwritten when `auto_n_inner = .true.`.
      logical :: auto_n_inner = .false.
         !! If `.true.`, derive `n_inner` at setup from the gravity-wave
         !! CFL.  Default `.false.`.
      real(wp) :: cfl_bt_safety = 0.65_wp
         !! Safety fraction on the shallow-water CFL bound when
         !! `auto_n_inner = .true.` (typical 0.65-0.7).
      real(wp) :: bebt = 0.1_wp
         !! BT continuity-flux velocity projection weight (MOM6 `BEBT`).
         !! Default `0.1` = MOM6's default (`MOM_barotropic.F90` get_param
         !! "BEBT", default=0.1): damps the barotropic gravity waves at
         !! `|λ|² = 1 − b·a²` per substep (`a = c·dt_bt·k_eff`), which is
         !! what damps the barotropic grid-scale mode under `pred_corr`.
         !! `0.0` = pure forward-backward Euler (neutral, the pre-2026-09-22
         !! default).  The FB stability limit tightens to
         !! `a ≤ 2/√(1+2·bebt)` (MOM6 `set_dtbt`'s `1+2·BEBT` factor).
      logical :: use_cont_type = .false.
         !! When `.true.`, the BT substep uses a piecewise-cubic
         !! flux-bounded closure instead of `uh = u·h_face`.  Default `.false.`.
      logical :: cont_corr_bounds = .false.
         !! When `.true.` (and `use_cont_type`), the η-correction bound
         !! uses the BT_cont face-by-face flux limits.  Default `.false.`.
      logical :: upstream_h_face = .false.
         !! When `.true.`, the BT chain uses per-face upstream-PPM
         !! column-sum thickness instead of centred `h_face`.  Default `.false.`.
      logical :: correction_h_weighted = .false.
         !! When `.true.`, the post-substep BT corrector distributes
         !! per-layer Δu by `h_face(k) / ⟨h⟩_h` instead of uniformly.
      logical :: correction_visc_rem = .false.
         !! Joint h·visc_rem weight for the BT corrector.  Requires
         !! `correction_h_weighted = .true.` (checked in `validate_config`;
         !! the uniform-Δu branch never reads `visc_rem`, so without it the
         !! knob would validate and silently do nothing).  `visc_rem` (the
         !! bt_work%visc_rem_u/v seam) is produced by
         !! `vdiff_apply_momentum` from the SAME factorized momentum
         !! tridiagonal with RHS ≡ 1 — see `rdb_ocean_vdiff.F90`.  Its rows
         !! sum to exactly 1 unless `&ocean_vdiff_nml implicit_drag=.true.`
         !! breaks the bed row, so with `implicit_drag=.false.` visc_rem ≡ 1
         !! identically and this knob is mathematically inert (warned, not
         !! an error — a legal, merely no-op configuration).
      character(len=16) :: split_scheme = "pred_corr"
         !! Outer time-integration scheme for the split-explicit ocean path
         !! (SPEC §4 S3/S4).  Both schemes are supported and both are under
         !! test; they differ in what they cost you.
         !!
         !! `"pred_corr"` (DEFAULT since 2026-09-14) = the MOM6
         !! predictor-corrector.  Off-centred predictor at `pc_be·dt`, slow
         !! tendencies (Coriolis-advection, horizontal viscosity) evaluated
         !! on the `u_av`/`h_av` step time-means, ONE prognostic update in
         !! the corrector, forward-backward gravity-wave pairing — neutrally
         !! stable to `ω·dt = 2`, which lifts the internal-wave `dt` ceiling
         !! and removes the resting-state energy growth described below.  It
         !! preserves the Eady benchmark's genuine baroclinic mode and makes
         !! it slightly stronger (max|v| ×2074 over 60 days vs ×1480).  The
         !! per-step cost is close to ssp_rk2's: `coriolis_coast` (48²×4, 20
         !! simulated days, one V100, quiescent machine) runs 70.3 s under
         !! pred_corr against 65.9 s under ssp_rk2 — about 7 % slower, not a
         !! different order.
         !!
         !! `validate_config` refuses `pred_corr` FAIL-LOUD, never silently,
         !! outside its v1 envelope: `eulerian_z`, `&ocean_wetdry_nml
         !! enable`, `dt_tracer_advect_ratio > 1`.  Six shipped namelists
         !! pin `ssp_rk2` for those reasons.
         !!
         !! `"ssp_rk2"` (**EXPERIMENTAL**) = two identical stages + SSP
         !! average.  Still fully supported and fully tested — experimental
         !! is a label on the ANSWER, not a deprecation of the code path,
         !! and the stability suite runs an `ssp_rk2` twin of every case
         !! whose namelist does not pin a scheme.  It has the widest
         !! envelope: every vcoord including `eulerian_z`, wet/dry, and
         !! `dt_tracer_advect_ratio > 1` — and it is the only scheme wired
         !! through the windowed tracer-advection path, which is why some
         !! namelists must pin it.
         !!
         !! **What it costs you, measured.**  The two-stage average
         !! amplifies an internal gravity wave by `√(1 + (ω·dt)⁴/4)` per
         !! step, so it MANUFACTURES energy from a motionless stratified
         !! state.  On
         !! `validation_examples/ocean/eady/resting_stratified_channel.nml`
         !! (flat bed, periodic, stably stratified, at rest, seeded with
         !! ±0.5 mK of noise, NO energy source of any kind) En reaches
         !! **2.992E-05 m²/s² by day 25** — 7.7 mm/s of current out of
         !! nothing — still climbing on a 2.5-day e-folding, where
         !! `pred_corr` on the identical file holds **1.739E-09** (17 000×
         !! less, 83-day e-folding).  The cause is the OUTER time splitting
         !! and nothing else: the Coriolis form, the ALE remap and the PGF
         !! form were each substituted and each moved the answer by < 0.1 %
         !! — all three are EXONERATED; removing the stratification dropped
         !! En 119×; removing the lateral viscosity RAISED it (viscosity
         !! damps the mode, it is not its source).
         !!
         !! **How to read that number.**  The error is a `(ω·dt)⁴` noise
         !! floor, set by how hard `dt` is pushed against the internal-wave
         !! period.  A forced, energetic, viscous run sits decades above the
         !! floor and never notices it; a quiescent, weakly-damped or long
         !! spin-up run does not, and there the manufactured energy is the
         !! signal.  The suite carries this as a scoped XFAIL on
         !! `resting_stratified_channel__ssp_rk2`, not as institutional
         !! memory.
      real(wp) :: pc_be = 0.6_wp
         !! pred_corr predictor fraction (MOM6 `BE`, 0.6 in the control run):
         !! the predictor's provisional velocity advances to `BE·dt`; only
         !! the corrector takes the full step.  Unused under ssp_rk2.
      logical :: renorm_visc_rem = .false.
         !! MOM6 continuity-inversion parity (SPEC §4 S2b): pass
         !! `visc_rem_u/v` into the slow-continuity transport-matching
         !! renormaliser, switching it to the γ-weighted `du` + Jacobian
         !! form (`u_cor = u + du·γ_k`)
         !! — a heavily-frictioned layer receives a smaller share of the
         !! barotropic correction than an undamped one, and the same
         !! weighting lands in the mass fluxes.  Behaviour change when
         !! on (the flux expression tree gains the γ factors).  Requires
         !! `correction_visc_rem = .true.` (the visc_rem producer);
         !! checked in `validate_config`.
      logical :: forcing_visc_rem = .false.
         !! MOM6 `wt_u` parity for the BT FORCING assembly:
         !! weight each layer's
         !! contribution to `F_bt_u/v` (and to the PGF-projection
         !! subtraction, which must use the same weights) by
         !! `h_face·visc_rem(k)` instead of `h_face`.  Layers the implicit
         !! friction will immediately damp — grounded sliver stacks under
         !! the vdiff BBL glue — then contribute nothing to the fast loop,
         !! which otherwise integrates the spurious grounded-layer PGF's
         !! depth-mean ballistically (PGF_BUG.md §9).  Requires
         !! `correction_visc_rem = .true.` (the visc_rem producer);
         !! checked in `validate_config`.
      logical :: correction_bc_pgf = .false.
         !! Adds a per-layer baroclinic-PGF retro-correction for the η
         !! change during the BT substep.  Requires `pgf%form = "fv_mom6"`.
      logical :: substep_drag = .false.
         !! Multiplies the per-face BT velocity update by a damping factor
         !! every inner step.
      logical :: substep_zeta_ke = .true.
         !! When `.true.` (default, bit-identical) the BT fast loop
         !! integrates its LIVE relative vorticity + KE gradient
         !! (`(ζ_bt+f)·v − ∇KE_bt`) every substep.  `.false.` = MOM6
         !! parity: the fast loop carries planetary Coriolis only
         !! (`f·v_at_u`, live velocities); ζ_bt-advection and ∇KE stay
         !! FROZEN inside the depth-mean forcing `F_bt_*_fast`, exactly
         !! like MOM6's `q = f/D` weights (`btstep_find_Cor` — planetary
         !! only, no ζ, no KE term).  The
         !! `subtract_fast_cor_ref` reference reduces to its `f·v̄` part
         !! so the τ=0 cancellation stays exact.  Motivation: the live
         !! nonlinear terms host an exponential ~4Δx wall-trapped BT
         !! mode on shelf rims (wall-biased KE stencil: wall faces
         !! contribute zero KE, so −∇KE points into walls and corners)
         !! that the `bound_kh`-clamped viscosity (λ·dt ≤ bound_coef/8)
         !! cannot damp at large dt — the 1024²×100 dt=600 SE-corner
         !! column evacuation (h-guard step 106).
      logical :: wave_drag = .false.
         !! Master switch for the barotropic linear (Rayleigh) wave drag
         !! — the bulk energy sink for the barotropic tide (Egbert & Ray
         !! 2001; Jayne & St Laurent 2001), MOM6 `BT_LINEAR_WAVE_DRAG`.
         !! Default `.false.` => bit-identical.  Composes multiplicatively
         !! with `substep_drag` inside `bt_rem_u/v`.
      character(len=32) :: wave_drag_form = "uniform"
         !! `r_H` filler: "uniform" (global scalar), "roughness_proxy"
         !! (resolved-bathymetry-variance proxy for Jayne & St Laurent
         !! `<h^2>`), or "file" (reserved for PR-14's NetCDF map reader —
         !! fail-loud not-implemented today).
      real(wp) :: wave_drag_scale = 1.0_wp
         !! Global tuning multiplier on `r_H` (MOM6 `BT_WAVE_DRAG_SCALE`).
      real(wp) :: wave_drag_r_uniform = 0.0_wp
         !! Piston velocity `r_H` (m/s) for `wave_drag_form = "uniform"`.
      real(wp) :: wave_drag_kappa = 6.2832e-4_wp
         !! Topographic wavenumber `kappa` (1/m) for
         !! `wave_drag_form = "roughness_proxy"`.  Same default as
         !! `&ocean_tidal_mixing_nml kappa_itides` — keep them equal.
      real(wp) :: wave_drag_n_bot = 1.0e-3_wp
         !! Reference near-bottom buoyancy frequency `N_bot` (1/s) for
         !! `wave_drag_form = "roughness_proxy"`.
      real(wp) :: wave_drag_h2_max = 2.5e4_wp
         !! Ceiling on the resolved-bathymetry `<h^2>` proxy (m^2; ⇔
         !! h_rms <= 158 m) for `wave_drag_form = "roughness_proxy"`.
      character(len=256) :: wave_drag_file = ""
         !! Reserved for PR-14 (MOM6 `BT_WAVE_DRAG_FILE`); unused today.
      character(len=64) :: wave_drag_var = "rH"
         !! Reserved for PR-14 (MOM6 `BT_WAVE_DRAG_VAR`); unused today.
      integer :: bt_halo = BT_HALO_AUTO_SENTINEL
         !! Wide-halo march-in width.  `-1` (default) = AUTO: resolved at
         !! configure to `BT_HALO_AUTO_WIDTH` (8) IFF the run is multi-rank
         !! (`compute_size > 1`) AND no march-in exclusion is active, else 0.
         !! `0` = explicit off (per-substep grouped exchange, v1 bit-identical).
         !! Even positive value: widen the BT ghost band to `nghost + bt_halo`
         !! and fire one grouped exchange every `bt_halo/2` substeps.  Odd
         !! values are rounded DOWN to even with a warning.  Fail-loud
         !! exclusions (apply when bt_halo is set EXPLICITLY > 0; AUTO instead
         !! silently resolves to 0): wet/dry enable, use_cont_type,
         !! upstream_h_face, tides enable, psurf enable, porous enable,
         !! cavity_dyn enable, supergrid/tripolar grid_config.  The set of record is
         !! `bt_halo_auto_exclusion`; keep it and the `validate_config`
         !! checks in lockstep.
   end type ocean_bt_config_t

   type :: ocean_debug_config_t
      !! `&ocean_debug_nml`: forensic probes for the ocean dyn-core.
      !! Every knob defaults off ⇒ bit-identical (untaken branches); all
      !! are HEAVY when on (device waits / serialised apply chains) —
      !! investigation tools, not production diagnostics.  Output is
      !! greppable fixed-format rows (BUDGET / KE_ATTR / CHKSUM /
      !! HOTFACE); reader: `tools/read_chksum.py`.
      logical :: budget = .false.
         !! Per-stage BT power-budget probe (BUDGET rows).
      logical :: ke_attr = .false.
         !! Per-segment layer-KE attribution meter (KE_ATTR rows) —
         !! serialises the velocity-apply chain while sampling.
      integer :: ke_attr_start_step = 0
         !! First outer step the KE meter samples (0 = from start).
      integer :: ke_attr_end_step = 0
         !! Last outer step the KE meter samples (0 = no upper bound).
      logical :: chksum = .false.
         !! MOM6-style per-phase field checksums (CHKSUM rows) + the
         !! HOTFACE argmax-face anatomy at the tendency seams — the
         !! first-diverging-operator / non-finite-minting attribution
         !! probe.  Device waits + reductions per phase seam.
      integer :: chksum_start_step = 0
         !! First outer step the chksum probe samples (0 = from start).
      integer :: chksum_end_step = 0
         !! Last outer step the chksum probe samples (0 = no upper bound).
      logical :: chksum_interior = .false.
         !! Restrict every chksum reduction to PHYSICAL cells (drop the
         !! ghost ring).  Default `.false.` = legacy whole-array output.
         !!
         !! Turn this ON to compare a 1-rank run against an N-rank run.
         !! With ghosts included the reduced value set genuinely differs
         !! between decompositions even for CORRECT code, so the
         !! otherwise decomposition-invariant `bits` column is unusable
         !! for that comparison.  Interior-only makes the value set
         !! identical by construction, and `bits` is an exact integer
         !! (associative, commutative) reduction — so a correct run
         !! yields the SAME `bits` on any layout, and the first stage
         !! whose `bits` disagree is where a decomposition bug lives.
   end type ocean_debug_config_t

   type :: ocean_mpi_config_t
      !! `&ocean_mpi_nml`: multi-rank MPI debug and tuning controls.
      logical :: poison_ghosts = .false.
         !! Debug knob — sentinel-NaN the exchange-covered ghost bands at
         !! each outer-step start; any unexchanged-ghost consumption becomes
         !! a loud NaN.  Default `.false.` = bit-identical (untaken branch).
   end type ocean_mpi_config_t

   type :: ocean_wetdry_config_t
      !! `&ocean_wetdry_nml`: dynamic wetting/drying for the split-explicit
      !! barotropic substep (docs/ocean_wetdry_plan.md).  Default off ⇒
      !! byte-identical everywhere (the wd_* workspaces stay unallocated and
      !! the substep runs its unmodified centred-face path).
      logical :: enable = .false.
         !! Master switch.  When `.true.` the BT substep switches to upwind
         !! face thickness + the positive-definite outflow limiter + the
         !! bed-blocking momentum gate; requires sigma / zstar(-lite)
         !! vertical coordinates, `&ocean_continuity_nml ppm_limit_pos`,
         !! and single-rank (all fail-loud at configure).
      real(wp) :: dry_depth = 0.05_wp
         !! Total-depth dry threshold (m): a column with `D < dry_depth`
         !! is dynamically dry (masked).  Must be > 0.
      real(wp) :: rewet_depth = 0.10_wp
         !! Hysteresis re-wet threshold (m): a dry column re-wets only when
         !! `D > rewet_depth`.  Must be > `dry_depth` (kills wet/dry
         !! flip-flop chatter at the front).
      real(wp) :: land_margin = 5.0_wp
         !! Static-land headroom (m) above rest MSL, used ONLY when enable=.true.
         !! A column is STATIC LAND (metrics zeroed at configure) iff its bed is
         !! above the highest credible water level: `b < -land_margin` (b positive-
         !! down, so -land_margin is `land_margin` metres above rest MSL). Columns
         !! with `-land_margin <= b < LAND_DEPTH_THRESHOLD` are INTERTIDAL: they
         !! seed wet_mask=1 (real metrics) and the dynamic wd_wet_dyn gate handles
         !! their wetting/drying. Set to the expected tidal amplitude + surge. When
         !! enable=.false. this knob is inert (the LAND_DEPTH_THRESHOLD seed runs).
   end type ocean_wetdry_config_t
   type :: ocean_pgf_config_t
      character(len=16) :: form = "mont"
         !! Pressure-gradient kernel variant.  "mont" (default) is the
         !! Boussinesq Montgomery-potential form — general purpose (valid over
         !! sloping bathymetry and every vcoord), algebraically identical to
         !! "fv_lite" on columns of equal layer thickness, and exact at rest in
         !! isopycnal columns; it is also the cheapest, since it builds no
         !! pressure stack.  "fv_lite" and "fv_wright" are the FV sigma-aware
         !! kernels.  "fv_lite" carries the SAME physics content as "mont"
         !! (layer-mean rho, no in-layer quadrature) and differs from it only
         !! where layer thicknesses are UNEQUAL across a face — there the two
         !! are different discretisations of the same term, neither exact.
         !! "gprime" is the reduced-gravity layered PGF (NK = 2 only).
         !! "fv_mom6" is the faithful FV-Bouss port (used by `gfs_scale`).
      real(wp) :: gprime_gfs = 9.81_wp
         !! Free-surface gravity (m/s²) for the gprime PGF.
      real(wp) :: gprime_gint = 0.0098_wp
         !! Internal-interface reduced gravity (m/s²) for the gprime PGF.
      real(wp) :: gfs_scale = 1.0_wp
         !! Free-surface gravity scaling for the FV_MOM6 PGF.  Default 1.0
         !! = no reduction.  When < 1, applies a Montgomery `dM` correction
         !! and evolves η at the reduced gravity.  Active only under
         !! `form = "fv_mom6"`.
      real(wp) :: maxvel = 0.0_wp
         !! Velocity-truncation clamp (m/s): `u = sign(u)·min(|u|, maxvel)`
         !! after each outer dyn step.  Default 0 = disabled.
      real(wp) :: cfl_trunc = 0.0_wp
         !! Advective-CFL velocity truncation threshold (nondimensional).
         !! When > 0, faces with `|u|·dt/dx > cfl_trunc` are clipped to
         !! `0.9·cfl_trunc·dx/dt` (runs before the `maxvel` cap).  Default
         !! 0.0 = disabled (bit-identical).
      logical :: mass_weight = .false.
         !! FV_MOM6 shelf-break mass-weighting.  When `.true.`, biases the
         !! layer density in the horizontal pressure integral toward the
         !! thinner column at unequal-depth faces, cancelling the spurious
         !! bottom-layer PGF.  Reduces to bit-identical on aligned columns.
         !! Default `.false.`.  Active only under `form = "fv_mom6"`.
      logical :: reconstruct_for_pressure = .false.
         !! FV_MOM6 in-layer T/S reconstruction (Adcroft, Hallberg &
         !! Harrison 2008; White, Adcroft & Hallberg 2009).  `.false.`
         !! (default) uses layer-mean (PCM) density ⇒ bit-identical.
         !! `.true.` builds each layer's pressure anomaly via a 5-point
         !! Boole quadrature of a monotone PLM/PPM T/S profile.  Active only
         !! under `form = "fv_mom6"` (fail-loud at configure otherwise).
      integer :: recon_scheme = 1
         !! In-layer reconstruction scheme: 1 = PLM, 2 = PPM.  Mirrors
         !! MOM6 `Recon_Scheme`.  Only consulted when
         !! `reconstruct_for_pressure = .true.`.
      logical :: p_top_in_bc = .false.
         !! Add the top-of-column load `multilayer_state_t%p_top` (Pa) to
         !! the FV_MOM6 pressure-stack surface boundary condition:
         !! `pa(nz+1) = rho_ref*g*eta_geo + p_top`.  Default `.false.` =>
         !! `pa(nz+1) = rho_ref*g*eta_geo`, bit-identical to every run
         !! before this knob existed.
         !!
         !! WHY IT IS NOT A DOUBLE COUNT.  A depth-uniform `p_top`
         !! perturbs EVERY layer's `PFu` by the same `-(1/rho_0)*grad
         !! p_top` (the theorem in `compute_fv_mom6_impl`'s docstring).
         !! The split solver subtracts the depth mean of the layer PGF
         !! from the barotropic forcing and then folds the barotropic
         !! solution back over the layers, so that uniform piece cancels
         !! identically and the load's barotropic response is carried by
         !! the `eta_forcing` seam alone — the two seams are orthogonal,
         !! not additive.  On the UNSPLIT driver (`n_inner = 0`) there is
         !! neither a depth-mean replacement nor a seam, so this term is
         !! the load's ONLY path into the momentum: a correction, not a
         !! duplicate.
         !!
         !! What it buys where the load is LARGE (an ice-shelf draft,
         !! `5e6 Pa`): `pa` is built as an anomaly about `rho_ref*g*z`,
         !! and with the load cancelled inside `pa(nz+1)` against
         !! `rho_ref*g*eta_geo` the whole stack stays `O(1e4 Pa)` instead
         !! of `O(5e6 Pa)`, which shrinks the `h_neglect` face-divisor
         !! leak by the same factor.
         !!
         !! FV_MOM6 ONLY (both the PCM and the `reconstruct_for_pressure`
         !! branch) — `validate_config` refuses it for any other `form`,
         !! which carries no `pa` stack to inject into.  Orthogonal to
         !! `&ocean_psurf_nml in_eos`: that knob puts the same `p_top`
         !! into the EOS's IN-SITU pressure ARGUMENTS, this one into the
         !! PGF's pressure BOUNDARY CONDITION.  Either, both or neither.
         !! `p_top` itself is produced today only by the `&ocean_psurf_nml`
         !! seam, so with `enable = .false.` it is the zero array and this
         !! knob is inert — a warning says so rather than leaving it
         !! silent.
   end type ocean_pgf_config_t
   type :: ocean_eos_config_t
      character(len=16) :: eos = "linear"
         !! Equation-of-state variant: "linear" (default, two-tracer),
         !! "wright" (Wright 1997 rational), "roquet_spv" (Roquet et al.
         !! 2015 specific-volume polynomial).  "roquet_spv" is incompatible
         !! with the "fv_wright" PGF (fails loud at configure).
      character(len=16) :: tfreeze_set = "seaice"
         !! Named seawater freezing-point (liquidus) coefficient set for
         !! `eos_freezing_point`, which evaluates the linear form
         !!
         !!   T_f = λ1·S + λ2 + λ3·p
         !!
         !! for every EOS variant (MOM6 keeps `TFREEZE_FORM = "LINEAR"`
         !! as its default under any density branch).
         !!
         !!   * `"seaice"` (DEFAULT ⇒ bit-identical) — the SIS2/MOM6
         !!     sea-ice liquidus, λ = (−0.054 °C/PSU, 0 °C,
         !!     −7.53e-8 °C/Pa).  `T_f(35 PSU, 0 Pa) = −1.89 °C`.  This is
         !!     the set the shipped sea-ice column model was ported and
         !!     tested against.
         !!   * `"isomip"` — the ISOMIP+ protocol liquidus (Asay-Davis et
         !!     al. 2016, GMD 9, Table 4 p. 2483; consumed in their
         !!     eq. (25) p. 2485), λ = (−0.0573 °C/PSU, 0.0832 °C,
         !!     −7.53e-8 °C/Pa).  Required for an ice-shelf-cavity run
         !!     claiming ISOMIP+ compliance, and the value the ice-shelf
         !!     literature is unanimous on.
         !!
         !! WHY IT MATTERS.  At S = 34.5 the two sets differ by ~0.03 °C
         !! — a few percent of a typical Antarctic thermal driving, and
         !! enough to flip the SIGN of a basal melt rate over a 0.03 °C
         !! band of ocean temperature.  A mistyped value is therefore a
         !! fail-loud `validate_config` error, never a silent fallback.
         !!
         !! NAMED SETS ONLY, on purpose: λ1/λ2/λ3 are a fitted triple and
         !! there is no free-form coefficient knob, so a configuration
         !! cannot mix λ1 from one source with λ2 from another.  A
         !! NONLINEAR liquidus (MOM6 `MILLERO_78`, a TEOS-10 polynomial)
         !! is a different functional form and would arrive as its own
         !! `form` selector at the documented seam in
         !! `eos_freezing_point`, not as another member of this list.
         !!
         !! SCOPE: this is the OCEAN-side liquidus only — what
         !! `eos_freezing_point` returns, i.e. the sea-surface freezing
         !! temperature the frazil, frazil-uptake and basal-flux kernels
         !! work against.  The SIS2 ice model's INTERNAL brine-pocket
         !! liquidus slope (`ICE_DTF_DS`, `rdb_ice_enthalpy`) is baked
         !! into its closed-form enthalpy<->temperature map and is NOT
         !! switched here; under `"isomip"` the two therefore disagree by
         !! ~0.03 °C.  The ISOMIP+ set is for ice-shelf-cavity work, where
         !! the sea-ice column model is normally off.
      real(wp) :: p_ref = 0.0_wp
         !! Reference pressure (Pa, `>= 0`) at which the model's POTENTIAL
         !! density `ms%rho_layer` is evaluated.  Default 0 (surface
         !! density, σ₀) ⇒ bit-identical to every run before this knob
         !! existed.  Read by the "wright" and "roquet_spv" variants;
         !! "linear" has no pressure dependence and ignores it.
         !!
         !! HORIZONTALLY UNIFORM BY DESIGN — it is a scalar and must stay
         !! one.  `rho_layer` is differenced ALONG a layer (the Montgomery
         !! PGF, the FV-lite / FV-MOM6-PCM integrands) and VERTICALLY (the
         !! vmix N² builders), so a reference pressure varying with (i,j)
         !! would give two columns of identical water at the same depth
         !! densities differing by `∂ρ/∂p·Δp` — a spurious along-layer
         !! gradient and hence a spurious pressure gradient force.  A
         !! spatially varying surface load goes to the IN-SITU builders via
         !! `&ocean_psurf_nml in_eos`, never here.
         !!
         !! What it buys: the thermobaric state at which the effective
         !! α/β are evaluated.  Near the freezing point at ice-shelf-cavity
         !! pressures that matters, so a cavity or abyssal study is better
         !! referenced to a representative depth (2.0e7 Pa ≈ 2000 dbar, the
         !! usual σ₂ choice) than to the surface.
         !!
         !! Distinct from `&vcoord_nml rho_ref_pressure`, which references
         !! the RHO / HYCOM target-density COORDINATE and the density-space
         !! diagnostic remap.  For a density-coordinate run the two should
         !! normally be set to the SAME value so coordinate and dynamics
         !! agree on what "density" means; they are kept independent
         !! because a diagnostic remap to a different reference is a
         !! legitimate request.
   end type ocean_eos_config_t
   type :: ocean_bdrag_config_t
      character(len=16) :: form = "quadratic"
         !! Bottom-drag variant: "quadratic" (default, log-layer
         !! du/dt = -C_d·|U|·u/h) or "linear" (Rayleigh du/dt = -r·u).
      real(wp) :: cd = 0.0_wp
         !! Quadratic drag coefficient (dimensionless).  Typical 2.5e-3.
         !! Zero disables the quadratic branch.
      real(wp) :: r = 0.0_wp
         !! Linear Rayleigh coefficient (1/s).  Active when `form = "linear"`.
      real(wp) :: hbbl = 0.0_wp
         !! Bottom-boundary-layer thickness (m) over which drag is
         !! distributed.  Zero (default) = bed-only mode.
      real(wp) :: bg_vel = 0.0_wp
         !! Background velocity floor (m/s) for the distributed quadratic
         !! form.  Typical 0.1 m/s.
      real(wp) :: bbl_thick_min = 0.0_wp
         !! Minimum effective BBL thickness (m).  Typical 0.1 m.
      real(wp) :: bed_factor = 1.0_wp
         !! Multiplier on the bed-layer (k=1) drag tendency only.  Default
         !! `1.0` = bit-identical.
      logical :: channel_drag = .false.
         !! Lateral side-wall (channel) drag: a per-layer Rayleigh sink at
         !! every layer with a partially blocked cross-stream perimeter.
         !! Default `.false.`.  Flat-bottom/all-wet ⇒ exact no-op.
      real(wp) :: cdrag_side = 0.0_wp
         !! Side-wall drag coefficient (dimensionless) for the channel drag.
         !! Zero (default) disables the branch.
      logical :: implicit = .false.
         !! Backward-Euler (implicit) bottom drag: `u^{n+1}=u/(1+dt·λ)`,
         !! unconditionally stable.  Default `.false.` = explicit
         !! forward-Euler (bit-identical) but unstable on thin shelf layers.
         !! Recommend `.true.` for shallow-shelf runs.
   end type ocean_bdrag_config_t

   type :: ocean_tdrag_config_t
      !! `&ocean_tdrag_nml` — ICE-SHELF TOP drag, the mirror of
      !! `&ocean_bdrag_nml` at `k = nz`.  Requires
      !! `&ocean_cavity_dyn_nml enable` (without a draft there is no ice
      !! base and `cover_frac` is identically zero, so the kernel would
      !! be a no-op with a cost).  Default `enable = .false.` ⇒ the slot
      !! is a placeholder, no kernel runs, every path is bit-identical.
      logical :: enable = .false.
         !! Master switch.  Requires `&ocean_cavity_dyn_nml enable`.
      character(len=16) :: form = "quadratic"
         !! Top-drag variant: "quadratic" (default, `du/dt =
         !! -C_d*|U|*u/h`, the ISOMIP+ prescription) or "linear"
         !! (Rayleigh `du/dt = -r*u`).  Enum mirrors
         !! `parse_tdrag_variant` in `rdb_ocean_top_drag`.
      real(wp) :: cd = 0.0_wp
         !! Quadratic drag coefficient (dimensionless).  ISOMIP+ value
         !! 2.5e-3 (Asay-Davis et al. 2016 Table 4).  Zero disables the
         !! quadratic branch.  When `&ocean_cavity_melt_nml enable`, this
         !! is THE `C_d` for both momentum and the melt friction
         !! velocity — see the `cdrag_top` agreement rule in
         !! `validate_config`.
      real(wp) :: r = 0.0_wp
         !! Linear Rayleigh coefficient (1/s).  Active when
         !! `form = "linear"`.
      real(wp) :: htbl = 0.0_wp
         !! Top-boundary-layer thickness (m) the stress is distributed
         !! over (the mirror of `&ocean_bdrag_nml hbbl`).  Zero (default)
         !! = layer-`nz`-only.  Positive values keep the explicit drag
         !! rate finite where a sigma coordinate thins the top layer near
         !! a grounding line.
      real(wp) :: bg_vel = 0.0_wp
         !! Background velocity floor (m/s) in the quadratic speed.
         !! Zero (default) ⇒ the layer-only quadratic branch is the exact
         !! algebraic mirror of the bottom drag's.
      real(wp) :: tbl_thick_min = 0.0_wp
         !! Minimum effective TBL thickness (m) in the `stress/h_tbl`
         !! denominator.  Zero (default) falls back to the kernel's
         !! `h_min` (1e-3 m), matching the bottom drag.
      logical :: implicit = .false.
         !! Backward-Euler top drag inside the drag kernel:
         !! `u^{n+1} = u/(1 + dt*lambda)`, unconditionally stable for any
         !! top-layer thickness.  Default `.false.` = explicit forward
         !! Euler (bit-identical to the pre-knob path, but conditionally
         !! unstable when `lambda*dt > 1`).  Mutually exclusive with
         !! `&ocean_vdiff_nml implicit_top_drag`.
   end type ocean_tdrag_config_t

   type :: ocean_hdiff_config_t
      !! `&ocean_hdiff_nml` — along-coordinate tracer Laplacian
      !! (`rdb_ocean_hdiff_tracer`).  Diffusion along the MODEL
      !! coordinate, not neutral surfaces (that is `&ocean_redi_nml`,
      !! a separate already-reachable capability).
      real(wp) :: kappa_h = 0.0_wp
         !! Constant horizontal tracer diffusivity (m^2/s).  Zero
         !! (default) is a no-op — the kernel short-circuits without
         !! touching hTr, so this knob is bit-identical when unset.
         !! Stability bound (explicit forward-Euler), checked at
         !! configure using the coarsest `dx`/`dy` (Cartesian metres):
         !!   kappa_h * dt_therm * (1/dx^2 + 1/dy^2) <= 0.5
   end type ocean_hdiff_config_t
   type :: ocean_hvisc_config_t
      character(len=16) :: lateral_closure = "none"
         !! Lateral-mixing closure tag.  "none" (default) keeps the scalar
         !! `nu_h`.  "leith" = Leith vorticity-gradient scaling;
         !! "smagorinsky"/"smag" = Smagorinsky strain-rate scaling.  Either
         !! gives a per-face viscosity floored at `ah_bg`, capped at `ah_max`.
      real(wp) :: c_smag = 0.15_wp
         !! Dimensionless Smagorinsky coefficient.  Default 0.15 (range 0.1–0.2).
      real(wp) :: c_leith = 1.0_wp
         !! Dimensionless Leith coefficient.  Typical range 1.0–2.0.
      real(wp) :: kh_vel_scale = 0.0_wp
         !! Velocity scale (m/s) for the grid-resolution viscosity floor
         !! `nu_min = kh_vel_scale · dx`.  Typical 3e-3 m/s.
      real(wp) :: ah_bg = -1.0_wp
         !! Background harmonic viscosity (m²/s) for the lateral closure.
         !! Negative (default) = derive from nu_h and kh_vel_scale.
      real(wp) :: ah_max = 1.0e4_wp
         !! Hard cap on the per-face harmonic viscosity (m²/s).
      real(wp) :: kh_vel_scale_live = 0.0_wp
         !! Live velocity-scale viscosity coefficient (m/s).  When positive,
         !! a state-dependent `A_vel = kh_vel_scale_live·dx·|u|` is
         !! max-combined into the per-face harmonic viscosity every step.
         !! Default 0 ⇒ off, bit-identical.
      real(wp) :: kh_aniso = 0.0_wp
         !! Anisotropic Laplacian viscosity magnitude (m²/s; Smith &
         !! McWilliams 2003).  When positive AND `stress_tensor = .true.`,
         !! adds direction-tensor terms aligned with `aniso_dir`.  Default 0
         !! ⇒ isotropic, bit-identical.  Stress-tensor path only.
      integer :: aniso_mode = 0
         !! Anisotropy-direction mode.  Only mode 0 (grid-relative
         !! `(n1,n2) = aniso_dir`) is implemented; other values rejected at
         !! configure (no silent fall-back).
      real(wp) :: aniso_dir(2) = [1.0_wp, 0.0_wp]
         !! Anisotropy direction vector `(n1,n2)` in grid i,j components.
         !! Default `(1,0)` = grid-i.  Normalised internally.
      logical :: smag_ah = .false.
         !! Smagorinsky biharmonic.  When `.true.` the biharmonic viscosity
         !! becomes flow-aware `nu4_face = C_b · L⁴ · |D|`, clamped to
         !! `[nu_4_bg, nu_4_max]`; the scalar `nu_4` is bypassed.
      real(wp) :: smag_bi_const = 0.06_wp
         !! Nondimensional biharmonic Smagorinsky constant.  Typical 0.015–0.06.
      real(wp) :: c_leith_bi = 0.0_wp
         !! Nondimensional biharmonic Leith constant.  Drives
         !! `lateral_closure="leith_biharm"` (per-face
         !! `nu4_face = C_lb · grid_sp⁶ · inv_PI6 · |∇²ζ|`, clamped to
         !! `[nu_4_bg, nu_4_max]`).  Default 0.0 = no-op (startup warning if
         !! `leith_biharm` selected at 0.0).
      real(wp) :: nu_4_bg = 0.0_wp
         !! Background biharmonic viscosity floor (m⁴/s) for the flow-aware path.
      real(wp) :: nu_4_max = 1.0e12_wp
         !! Static hard cap on per-face biharmonic viscosity (m⁴/s).  The
         !! forward-Euler biharmonic-CFL bound is enforced automatically per
         !! cell (scaled by `bound_coef`), so this is just an extra ceiling.
      real(wp) :: nu_h = 0.0_wp
         !! Constant horizontal eddy viscosity for momentum (m²/s).  Default
         !! 0 = bit-identical; wind-driven runs need this set.
      real(wp) :: nu_4 = 0.0_wp
         !! Constant biharmonic eddy viscosity for momentum (m⁴/s).  Default 0.
      logical :: no_slip = .false.
         !! Lateral coast BC (static land mask).  `.false.` (default) =
         !! FREE-SLIP (corner vorticity/shear strain multiplied by `wet_q`);
         !! `.true.` = NO-SLIP (factor `2 - wet_q`).  All-wet ⇒ bit-identical.
      logical :: stress_tensor = .false.
         !! Horizontal-viscosity operator.  `.false.` (default) =
         !! velocity-Laplacian `A·∇²u` with global caps, bit-identical.
         !! `.true.` = thickness-weighted stress-divergence form
         !! (momentum-conserving) with per-cell CFL limiter and coast
         !! masking; reduces to the Laplacian on uniform-grid/h/all-wet.
      real(wp) :: bound_coef = 0.8_wp
         !! CFL safety coefficient for the per-cell viscosity limiters
         !! (harmonic when `stress_tensor`, biharmonic always when `nu_4 > 0`
         !! or Smag_AH).  Default 0.8; must be ≤ 1.
      logical :: bound_kh = .false.
         !! MOM6 `BOUND_KH` analogue for the velocity-Laplacian paths:
         !! per-face clamp of the harmonic viscosity to
         !! `bound_coef·0.125/(dt·(1/dx²+1/dy²))` (~¼ of the forward-
         !! Euler stability limit).  Prevents an over-large `nu_h` (or
         !! flow-aware coefficient) from putting the frozen depth-mean
         !! viscous forcing on the barotropic mode into the phase-
         !! reversed ANTI-damping regime at grid scale, which pumps an
         !! exponential rim-trapped barotropic mode via the split-RK2
         !! Δu corrector (the 600² Lagrangian double-gyre blow-up).
         !! Recommended `.true.` for any config where
         !! `nu·dt·(1/dx²+1/dy²) ≳ 0.1`.  Default `.false.` ⇒
         !! bit-identical.
      logical :: resoln_scaled_visc = .false.
         !! Resolution-scaled Laplacian viscosity (Hallberg 2013).  When
         !! `.true.` the dynamic Laplacian viscosity is multiplied by the
         !! VarMix resolution function `Res_fn ∈ [0,1]` before the `ah_max`
         !! clamp (biharmonic `nu_4` never scaled).  Requires
         !! `&ocean_varmix_nml enable=.true.`.  Default `.false.` ⇒ bit-identical.
   end type ocean_hvisc_config_t
   type :: ocean_vmix_config_t
      logical :: use_closure = .true.
         !! Master switch for the vmix (vertical-mixing) module.  When
         !! `.true.` (default), the chosen interior closure (PP81)
         !! populates `vmix%kv` / `vmix%kt` from the local Richardson
         !! number and vdiff reads them.  When `.false.`, vdiff uses
         !! its scalar `K_v_*` defaults and vmix kernels never run.
      logical :: use_kpp = .true.
         !! KPP surface-boundary-layer overlay.  Requires `use_closure`.
         !! Adds the cubic-profile kv overlay plus non-local γ_T/γ_S transport.
      logical :: direct_stress = .false.
         !! When `.true.` the wind stress is distributed across the top
         !! `hmix_stress` metres rather than the surface-most layer.
      real(wp) :: hmix_stress = 20.0_wp
         !! Thickness (m) of the surface slab for `direct_stress`.  Default 20 m.
      real(wp) :: kv_ml_invz2 = 0.0_wp
         !! Extra near-surface vertical viscosity (m²/s) with a `1/(z·hmix)²`
         !! profile in the top `hmix_fixed` metres.  Zero (default) = no change.
      real(wp) :: hmix_fixed = 20.0_wp
         !! Mixed-layer thickness (m) for the `kv_ml_invz2` profile.  Default 20 m.
      logical :: harmonic_visc = .false.
         !! When `.true.` vdiff uses the harmonic mean of adjacent layer
         !! thicknesses in the face-thickness denominator.  Default `.false.`
         !! = arithmetic mean.
      integer :: dt_therm_ratio = 1
         !! Thermodynamic + tracer-advection step runs at `ratio · dt_dynamic`.
         !! Default 1.
      integer :: dt_tracer_advect_ratio = 1
         !! Horizontal tracer advection runs every `ratio · dt_dynamic` via
         !! accumulated mass fluxes drained in one windowed advect.  1
         !! (default) ⇒ every step, bit-identical.  `dt_therm_ratio` must be
         !! an integer multiple of this (configure-time check).
      real(wp) :: kv_max = huge(1.0_wp)
         !! Ceiling on momentum viscosity kv (m^2/s).  Default `huge` = no clip.
      real(wp) :: kd_max = huge(1.0_wp)
         !! Ceiling on tracer diffusivity kt/ks (m^2/s).  Default `huge` = no clip.
      integer :: kd_smooth_iterations = 0
         !! 1-2-1 horizontal smoothing passes on kv/kt at interfaces.  Default 0 = off.
      logical :: vmix_guard = .false.
         !! Debug-gated negative/NaN diffusivity guard.  Default off.
      logical :: bkgnd_profile = .false.
         !! Enable the Bryan & Lewis (1979) depth-varying background
         !! diffusivity instead of the constant interior background.
         !! Default off (bit-identical).
      real(wp) :: bkgnd_kd_sfc = 1.0e-5_wp
         !! Surface-asymptote background tracer diffusivity (m^2/s).
      real(wp) :: bkgnd_kd_deep = 1.3e-4_wp
         !! Deep-asymptote background tracer diffusivity (m^2/s).
      real(wp) :: bkgnd_z0 = 2500.0_wp
         !! Bryan-Lewis transition-centre depth (m, positive down).
      real(wp) :: bkgnd_delta = 222.0_wp
         !! Bryan-Lewis transition half-width (m).
      real(wp) :: bkgnd_prandtl = 1.0_wp
         !! Background Prandtl number: Kv_bg = bkgnd_prandtl * Kd_bg.
      logical :: bkgnd_henyey = .false.
         !! Henyey, Wright & Flatte (1986) JGR 91:8487 latitude-dependent
         !! internal-wave factor (constant-`N0` simplification of Harrison &
         !! Hallberg 2008, JPO 38:1894), scaling the SCALAR background tracer
         !! diffusivities `kt_bg`/`ks_bg` and floored at `bkgnd_kd_min`.
         !! MUTUALLY EXCLUSIVE with `bkgnd_profile` (Bryan-Lewis) and
         !! requires a non-cartesian `grid_config` — both fail-loud at
         !! configure.  Default off ⇒ bit-identical.
      real(wp) :: bkgnd_kd_min = -1.0_wp
         !! Minimum background tracer diffusivity (m^2/s) under the Henyey
         !! scaling — `max(Kd_min, Kd*L(phi))`.  Negative = unset ⇒ resolved
         !! to 0.01*kt_bg (the reference `KD_MIN` default).
      real(wp) :: bkgnd_henyey_n0_2omega = 20.0_wp
         !! Ratio of the assumed reference buoyancy frequency N0 to twice
         !! the planetary rotation rate (nondim).
      real(wp) :: bkgnd_henyey_max_lat = 95.0_wp
         !! Latitude (degN) poleward of which the Henyey factor is reset to
         !! its minimum floor; compared against |latitude| so both
         !! hemispheres clamp.  > 90 so inert for any real latitude out of
         !! the box.
      character(len=8) :: tracer_recon = "ppm"
         !! Face-reconstruction scheme for the WINDOWED horizontal
         !! tracer-advection drain (`dt_tracer_advect_ratio > 1` path):
         !! ppm (default, CW-PPM, bit-identical) | weno5 | weno7 | weno9
         !! (WENO-Z swept-average, rung-adaptive; reuses the coastal
         !! reconstruction ladder).  Per-rung nghost minimum (weno5→3,
         !! weno7→4, weno9→5) is fail-loud at configure.  This is the
         !! SEPARATE ocean knob — independent of the coastal `tracer_recon`.

      ! ---- PP81 interior closure + KPP BL-depth constants ----
      ! Each default equals the current `ocean_vmix_t` field default
      ! (rdb_ocean_vmix.F90) so an nml that never mentions these keys is
      ! bit-identical.  `kpp_*` names (not bare `ri_crit`/`cs_nonlocal`/
      ! `c_vt2`) avoid a collision with `&ocean_kshear_nml ri_crit` (a
      ! DIFFERENT physical quantity, the JHL08 critical Ri).  These route
      ! to the OCEAN vmix slot — do not confuse with the similarly-named
      ! `&nonhydrostatic_nml kpp_*` keys, which are the live COASTAL path.
      real(wp) :: pp81_nu0 = 1.0e-2_wp
         !! PP81 Richardson-dependent viscosity scale (m^2/s).
      real(wp) :: pp81_nu_bg = 1.0e-4_wp
         !! PP81 background viscosity (m^2/s).  Also seeds `vmix%kv_bg`
         !! (the assembly floor) — see `vmix_seed_backgrounds`.
      real(wp) :: pp81_kappa_bg = 1.0e-5_wp
         !! PP81 background diffusivity (m^2/s).  Also seeds
         !! `vmix%kt_bg`/`vmix%ks_bg`.
      real(wp) :: pp81_alpha = 5.0_wp
         !! PP81 Richardson-number scaling coefficient (paper value 5;
         !! implementations vary 4-10).
      real(wp) :: shear2_floor = 1.0e-10_wp
         !! Floor on |du/dz|^2 + |dv/dz|^2 in the PP81 Ri denominator (1/s^2).
      real(wp) :: kpp_ri_crit = 0.3_wp
         !! Critical bulk Richardson number for the KPP BL-depth sweep
         !! (LMD94 §3; MOM6 KPP_BULK_RI default 0.3).
      real(wp) :: kpp_cs_nonlocal = 6.3_wp
         !! KPP non-local (counter-gradient) transport coefficient C_s
         !! (LMD94 eq 20, limit value).
      real(wp) :: kpp_c_vt2 = 1.8_wp
         !! KPP unresolved-turbulence coefficient for the V_t^2 term in
         !! the bulk-Ri denominator (LMD94 eq 23).  0 disables V_t^2.

      ! ---- Source of the thermal-expansion / haline-contraction pair ----
      character(len=8) :: buoyancy_coeffs = "constant"
         !! Where the vmix closures take their α (thermal expansion) and
         !! β (haline contraction) from — `"constant"` (DEFAULT) or
         !! `"eos"`.  Parsed by `parse_buoyancy_coeffs`
         !! (`rdb_ocean_vmix.F90`); keep the `allowed=` list and that
         !! routine's `case` arms in lockstep.
         !!
         !!   * `"constant"` — the scalar `&ocean_ic_nml alpha_T` /
         !!     `beta_S` off the EOS handle, whatever the active EOS is.
         !!     Historical behaviour ⇒ **bit-identical** for every shipped
         !!     namelist.  For `eos = "linear"` these ARE the true
         !!     coefficients, so the setting is physically exact there.
         !!   * `"eos"` — `eos_buoyancy_coeffs` evaluated per column /
         !!     per interface from the ACTIVE equation of state
         !!     (`α = −∂ρ/∂T`, `β = +∂ρ/∂S`, both analytic).  Under
         !!     `eos = "linear"` it returns those same handle members
         !!     bit-for-bit, so the two settings are byte-identical
         !!     there — the knob only bites under a NONLINEAR EOS
         !!     (`wright` / `roquet`).
         !!
         !! Routes THREE consumers: the KPP `B_0` surface buoyancy flux
         !! in both passes of `vmix_kpp_overlay_impl` (and hence the
         !! non-local γ gate, which switches on `B_0 < 0`), and the
         !! double-diffusion density ratio `R_ρ = α·ΔT / β·ΔS` in
         !! `vmix_split_ddiff_*_impl`.  Every OTHER α/β-like quantity on
         !! the ocean path already tracks the active EOS: EPBL,
         !! kappa-shear, tidal mixing, the isopycnal slopes and Redi go
         !! through `eos_specvol_derivs`, while PP81's N², the convective
         !! trigger, the wave speed and MLE difference `ms%rho_layer`
         !! itself.
         !!
         !! Why it matters: seawater's thermal expansion collapses toward
         !! zero near the freezing point and roughly doubles by 1000 dbar,
         !! so under an ice shelf a constant α mis-sizes the melt-driven
         !! buoyancy flux that sets the boundary layer.  `validate_config`
         !! WARNS (does not refuse) on cavity melt × nonlinear EOS ×
         !! `"constant"`.
   end type ocean_vmix_config_t
   type :: ocean_vdiff_config_t
      !! Backward-Euler vertical-friction solver knobs (`&ocean_vdiff_nml`).
      !! Both default `.false.` ⇒ the vdiff tridiagonal is built exactly as
      !! today and the explicit wind-stress / bottom-drag applies stay live
      !! ⇒ bit-identical to prior production.
      logical :: implicit_stress = .false.
         !! Fold the surface wind stress into the vdiff surface (k=nz) RHS
         !! row (Neumann top-BC) instead of the explicit pre-solve add.
      logical :: accel_visc_rem = .false.
         !! MOM6 `MOM_dynamics_split_RK2` parity: attenuate the slow
         !! EXPLICIT accelerations by the per-layer viscous remnant —
         !! `u = u_entry + visc_rem·(u_applied − u_entry)` after the
         !! CorAdv/PGF/hvisc/drag applies, before the BT correction — so
         !! friction-dominated near-massless layers cannot receive a
         !! full-strength dt·F kick (MOM6 `u = u_init + dt·visc_rem_u·
         !! (CAu + PFu + diffu)`).  Requires `&ocean_bt_nml
         !! correction_visc_rem` (the visc_rem producer; fail-loud at
         !! configure).  Split path only (v1).  Default off ⇒ bit-identical.
      logical :: implicit_drag = .false.
      logical :: hvel_mom6 = .false.
      real(wp) :: hbbl_visc = 10.0_wp
         !! Fold the bottom drag into the vdiff bed (k=1) diagonal (stress
         !! bottom-BC) instead of the explicit pre-solve add.  Mutually
         !! exclusive with `&ocean_bdrag_nml implicit` (split-apply) and
         !! incompatible with HBBL-distributed drag (`hbbl > 0`); both fail
         !! loud at configure.
      logical :: implicit_top_drag = .false.
         !! Fold the ICE-SHELF TOP drag into the vdiff `k = nz` DIAGONAL
         !! (`&ocean_tdrag_nml`'s mirror of `implicit_drag`) instead of
         !! the explicit pre-solve add.  The wind stress already owns
         !! that row's RHS; a drag is a diagonal term, so the two
         !! compose — but on a face the ice covers, the wind RHS is
         !! MASKED OFF here (there is no atmosphere under a shelf).
         !! That masking is now BELT AND BRACES and kept deliberately:
         !! the cover mask zeroes the `tau` pair at its source, so the
         !! factor multiplies zero, and it stays so the fold is correct
         !! STANDALONE if a later forcing path ever writes `tau` after
         !! the configure-time mask.
         !! Requires `&ocean_tdrag_nml enable`; mutually exclusive with
         !! `&ocean_tdrag_nml implicit` (both would damp the top layer)
         !! and with `htbl > 0` (the fold is one `k = nz` rate and
         !! cannot represent a distributed band).  All fail loud at
         !! configure.  Default `.false.` ⇒ bit-identical.
      logical :: bbl_glue = .false.
         !! MOM6 `bottomdraglaw` coupling parity: raise the momentum-solve
         !! interface viscosity to `kv_bbl` within botfn reach of the bed
         !! (harmonic-z bookkeeping) and replace the bed-row `dt·λ` Rayleigh
         !! fold with the piston `kv_bbl/(min(hvel₁/2, bbl_thick))` — the
         !! absorber that keeps MOM6's layered rest state at rest under the
         !! same spurious grounded-layer PGF we compute (PGF_BUG.md §9).
         !! Requires `hvel_mom6`, `implicit_drag`, and
         !! `&ocean_bdrag_nml form="linear"`; all enforced at configure.
      real(wp) :: bbl_piston = 3.0e-4_wp
         !! BBL drag piston velocity u* (m/s) for `bbl_glue` — MOM6
         !! `CDRAG·DRAG_BG_VEL` (0.003·0.1 at reference defaults);
         !! `kv_bbl = bbl_piston·hbbl_visc`.
      logical :: hvel_upwind = .true.
         !! Near-bed upwind (arithmetic-donor) blend in the `hvel_mom6`
         !! face-thickness build.  Default `.true.` = MOM6 parity /
         !! bit-identical hvel_mom6 behaviour.  `.false.` = pure harmonic
         !! hvel — the blend's u-sign test flip-flops on roundoff
         !! velocities at (near-)rest and collapses the BBL glue at
         !! whichever faces flip (PGF_BUG.md §9.8).
   end type ocean_vdiff_config_t
   type :: ocean_epbl_config_t
      !! Energetics-based planetary boundary layer (`&ocean_epbl_nml`).
      !! All defaults preserve bit-identity: `enable = .false.` keeps
      !! the existing PP81 + KPP path untouched.  Knob table + MOM6
      !! name mapping: `docs/generated_nml_knobs.md`.
      logical :: enable = .false.
         !! Master switch.  Requires `vmix%use_closure`; replaces the
         !! KPP overlay (configure logs the override).
      character(len=16) :: mstar_scheme = "om4"
         !! "constant" / "om4" / "rh18" (MOM6 EPBL_MSTAR_SCHEME).
      real(wp) :: mstar = 1.2_wp
         !! Constant-scheme mstar (MOM6 MSTAR).
      real(wp) :: mstar_cap = -1.0_wp
         !! Cap for OM4/RH18; off when < 0 (MOM6 MSTAR_CAP).
      real(wp) :: mstar_coef1 = 0.3_wp
         !! OM4 stabilizing coefficient (MOM6 MSTAR2_COEF1).
      real(wp) :: c_ek = 0.085_wp
         !! OM4 Ekman coefficient (MOM6 MSTAR2_COEF2).
      real(wp) :: mstar_conv_adj = 0.0_wp
         !! Convective mstar reduction in [0,1] (MOM6 MSTAR_CONV_ADJ).
      real(wp) :: rh18_cn1 = 0.275_wp
      real(wp) :: rh18_cn2 = 8.0_wp
      real(wp) :: rh18_cn3 = -5.0_wp
      real(wp) :: rh18_cs1 = 0.2_wp
      real(wp) :: rh18_cs2 = 0.4_wp
         !! RH18 mstar fits (MOM6 RH18_MSTAR_CN1..CS2).
      real(wp) :: nstar = 0.2_wp
         !! Convective PE -> TKE efficiency (MOM6 NSTAR).
      real(wp) :: tke_decay = 2.5_wp
         !! Ekman-depth / TKE-decay-scale ratio (MOM6 TKE_DECAY).
      real(wp) :: wstar_ustar_coef = 1.0_wp
         !! Convective weight in the velocity scale (MOM6 WSTAR_USTAR_COEF).
      character(len=16) :: vel_scale_scheme = "cube_root"
         !! "cube_root" / "rh18" (MOM6 EPBL_VEL_SCALE_SCHEME).
      real(wp) :: vstar_scale_fac = 1.0_wp
         !! Overall vstar multiplier (MOM6 EPBL_VEL_SCALE_FACTOR).
      real(wp) :: vstar_surf_fac = 1.2_wp
         !! RH18 mechanical surface vstar factor (MOM6 VSTAR_SURF_FAC).
      real(wp) :: von_karman = 0.41_wp
         !! kappa in Kd = vstar*kappa*mixlen (MOM6 VON_KARMAN_CONST).
      real(wp) :: ekman_scale_coef = 1.0_wp
         !! Rotational mixing-length rolloff (MOM6 EKMAN_SCALE_COEF).
      real(wp) :: min_mix_len = 0.0_wp
         !! Mixing-length floor, m (MOM6 EPBL_MIN_MIX_LEN).
      real(wp) :: mixlen_exponent = 2.0_wp
         !! Shape-function exponent (MOM6 MIX_LEN_EXPONENT).
      real(wp) :: translay_scale = 0.1_wp
         !! Transition-layer shape floor, must be in [0,1) when
         !! iterating (MOM6 EPBL_TRANSITION_SCALE).
      logical :: mld_iteration = .true.
         !! Self-consistent MLD root-find (MOM6 USE_MLD_ITERATION).
      real(wp) :: mld_tol = 1.0_wp
         !! MLD convergence tolerance, m (MOM6 EPBL_MLD_TOLERANCE).
      integer :: mld_max_its = 20
         !! Max MLD iterations (MOM6 EPBL_MLD_MAX_ITS).
      logical :: mld_bisection = .false.
         !! Bisection instead of false position (MOM6 EPBL_MLD_BISECTION).
      logical :: mld_use_prev_guess = .false.
         !! Seed from the previous step's MLD (MOM6 MLD_ITERATION_GUESS).
      real(wp) :: omega = 7.2921e-5_wp
         !! Earth rotation rate, 1/s (MOM6 OMEGA).
      real(wp) :: omega_frac = 0.0_wp
         !! Blend |f| with 2*Omega (MOM6 ML_OMEGA_FRAC).
      real(wp) :: prandtl = 1.0_wp
         !! Kv = prandtl*Kd into the momentum solve (MOM6 EPBL_PRANDTL).
      character(len=8) :: combine = "add"
         !! "add" / "max" vs the interior closure's kv/kt
         !! (MOM6 EPBL_IS_ADDITIVE).
      logical :: tke_diags = .false.
         !! Compute the per-column TKE budget diagnostic terms.
      logical :: use_lt = .false.
         !! Langmuir-turbulence enhancement (LF17 wind-only path,
         !! no wave model needed).  MOM6 EPBL_LT / USE_LA_LI2016.
      character(len=16) :: lt_scheme = "rescale"
         !! "rescale" (multiplicative) / "additive"
         !! (MOM6 EPBL_LANGMUIR_SCHEME).
      real(wp) :: lt_enhance_coef = 0.447_wp
         !! Enhancement coefficient (MOM6 LT_ENHANCE_COEF).
      real(wp) :: lt_enhance_exp = -1.33_wp
         !! Langmuir-number exponent (MOM6 LT_ENHANCE_EXP).
      real(wp) :: lt_max_enhance = 5.0_wp
         !! Cap on the multiplicative enhancement.
      real(wp) :: la_frac_hbl = 0.04_wp
         !! Stokes SL-average depth fraction (MOM6 LA_DEPTH_RATIO).
      real(wp) :: lt_lac1 = -0.87_wp
      real(wp) :: lt_lac2 = 0.0_wp
      real(wp) :: lt_lac3 = 0.0_wp
      real(wp) :: lt_lac4 = 0.95_wp
      real(wp) :: lt_lac5 = 0.95_wp
         !! Stability-modified La coefficients (MOM6 LT_MOD_LAC1..5);
         !! set all to 0 for the unmodified La.
   end type ocean_epbl_config_t
   type :: ocean_wave_speed_config_t
      !! `&ocean_wavespeed_nml`: first-baroclinic wave speed + deformation
      !! radius.  Diagnostic, default off — feeds GM/Redi/MEKE resolution scaling.
      logical :: enable = .false.
         !! Master switch.  Default off — bit-identity preserved.
      real(wp) :: mono_n2 = -1.0_wp
         !! DEFERRED: N^2-monotonising depth.  `< 0` = off.
      logical :: use_ebt = .false.
         !! DEFERRED: equivalent-barotropic / pressure-Neumann variant.
      integer :: n_wavespeed = 1
         !! Recompute cadence (every N steps); slow diagnostic.
   end type ocean_wave_speed_config_t

   type :: ocean_foxkemper_config_t
      !! `&ocean_foxkemper_nml`: Fox-Kemper mixed-layer-eddy restratification.
      !! Default off ⇒ bit-identical.  Reads `epbl%mld`.
      logical :: enable = .false.
         !! Master switch (off => bit-identity).
      real(wp) :: ce = 0.0625_wp
         !! FK08 coefficient Ce (0.06-0.08).
      real(wp) :: f_floor = 1.0e-5_wp
         !! |f| regularisation floor (1/s).
      real(wp) :: mld_decay_time = 0.0_wp
         !! Running-mean MLD filter time-scale (s) used when the MLD is
         !! retreating.  0 (default) => filter off => instantaneous EPBL MLD
         !! (bit-identical).  Positive damps step-to-step MLD swings.
      real(wp) :: tail_dh = 0.0_wp
         !! mu cubic-tail extension (smoother); default 0 = exact mu.
      logical :: use_mom_mixrate = .false.
         !! FK11 momentum-mixrate timescale vs bare Ce/|f|.  Default off.
      logical :: resolution_taper = .false.
         !! B2 res_fn double-counting hook (hard error if on without B2).
      logical :: use_bodner = .false.
         !! Bodner et al. (2023) frontogenesis-arrest MLE variant: replaces
         !! the `Ce/|f|` timescale with `Cr·Δs·|f|·h/w'u'`, the frontal-arrest
         !! length set by boundary-layer turbulence.  Overrides `ce`/
         !! `use_mom_mixrate`.  Default off ⇒ classic Fox-Kemper.
      real(wp) :: cr = 0.0_wp
         !! Bodner (2023) efficiency coefficient `Cr` (nondim).  Default 0
         !! (MOM6 default) ⇒ no Bodner transport until set (~0.02-0.08).
      real(wp) :: bodner_mstar = 0.5_wp
         !! Bodner mechanical (u*) weight in w'u' = (mstar·u*³+nstar·w*³)^{2/3}.
      real(wp) :: bodner_nstar = 0.066_wp
         !! Bodner convective (w*) weight in the same w'u' estimate.
      real(wp) :: min_wstar2 = 1.0e-24_wp
         !! Floor on w'u' (m²/s²) — pure 1/0 armour when u* and the surface
         !! buoyancy flux both vanish.
   end type ocean_foxkemper_config_t

   type :: ocean_kappa_shear_config_t
      !! `&ocean_kappa_shear_nml`: shear-driven interior turbulence (Jackson,
      !! Hallberg & Legg 2008).  Coexists with KPP/EPBL/PP81 — kappa is ADDED
      !! to the interior diffusivities.  Default off ⇒ bit-identical.
      logical :: enable = .false.
         !! Master switch.  Requires `vmix%use_closure` + thermodynamics
         !! (validated at configure).
      real(wp) :: ri_crit = 0.25_wp
         !! Critical Richardson number (JHL08 Ri_c).
      real(wp) :: shearmix_rate = 0.089_wp
         !! Shear-source-rate coefficient (JHL08).
      real(wp) :: fri_curvature = -0.97_wp
         !! Ri-function curvature in the shear source (JHL08).
      real(wp) :: c_n = 0.24_wp
         !! TKE decay-rate coefficient vs stratification N (JHL08).
      real(wp) :: c_s = 0.14_wp
         !! TKE decay-rate coefficient vs shear S (JHL08).
      real(wp) :: lambda = 0.82_wp
         !! Buoyancy mixing-length-scale coefficient (JHL08).
      real(wp) :: lz_rescale = 1.0_wp
         !! Boundary-distance length-scale rescale factor (JHL08).
      real(wp) :: kappa_0 = 1.0e-7_wp
         !! Background diffusivity, m^2/s; also the pre-step kappa (JHL08).
      real(wp) :: kappa_seed = 1.0_wp
         !! Iteration seed diffusivity, m^2/s (JHL08).
      real(wp) :: kappa_trunc = 1.0e-9_wp
         !! Diffusivity below this is truncated to 0, m^2/s (JHL08).
      real(wp) :: tke_bg = 0.0_wp
         !! Background TKE, m^2/s^2; Q is a denominator, floored (JHL08).
      real(wp) :: tol_err = 0.1_wp
         !! Picard convergence tolerance (JHL08).
      integer :: max_inner_it = 50
         !! Inner Picard iteration cap (JHL08).
      integer :: max_substep_it = 13
         !! Outer adaptive-substep iteration cap (JHL08).
      real(wp) :: src_max_chg = 10.0_wp
         !! Adaptive-dt source-change tolerance band (JHL08).
      real(wp) :: prandtl_turb = 1.0_wp
         !! Kv = prandtl_turb * Kd into the momentum solve (JHL08).
      real(wp) :: vel_underflow = 0.0_wp
         !! Velocity snap-to-zero magnitude, m/s, in the projection (JHL08).
      logical :: massless_merge = .false.
         !! D4: merge vanished (< H_VANISHED) layers onto the massive
         !! sub-grid before the column solve (vs the blunt gather floor).
         !! Default off (bit-identical); identity columns bypass the merge.
      logical :: at_vertex = .false.
         !! Solve the JHL08 columns at C-grid CORNERS (vorticity points)
         !! from the native face velocities, then average the corner Kd
         !! back to tracer points (MOM6 VERTEX_SHEAR — the OM5-class
         !! production setting).  Default off ⇒ bit-identical.  v1 scope
         !! is Kd only: Kv stays `prandtl_turb * kd_int` at cell centres
         !! (the corner->face viscosity seam is a follow-up PR).
      logical :: vertex_geometric_mean = .false.
         !! Geometric (vs arithmetic) mean in the corner->centre Kd
         !! average (MOM6 VERTEX_SHEAR_GEOMETRIC_MEAN).  Pair with
         !! `vertex_geomean_kdmin` — with a 0 floor the geometric mean
         !! hard-zeros Kd along every shear-zone edge.
      real(wp) :: vertex_geomean_kdmin = 0.0_wp
         !! Floor, m^2/s, applied to each corner Kd BEFORE the geometric
         !! mean (MOM6 VERTEX_SHEAR_GEOMETRIC_MEAN_KDMIN; inert unless
         !! `vertex_geometric_mean`).  OM5 configs use 1e-9.
   end type ocean_kappa_shear_config_t
   type :: ocean_slopes_config_t
      !! `&ocean_slopes_nml`: isopycnal (neutral) slope diagnostics (Griffies
      !! 1998) at C-grid interfaces.  Foundational gate for GM/Redi/VarMix.
      !! Default off ⇒ no-op.
      logical :: enable = .false.
         !! Master switch.  Default off ⇒ `ocean_slopes_compute` no-ops.
      real(wp) :: kd_smooth = 1.0e-6_wp
         !! Vert-fill smoothing diffusivity, m^2/s (× dt fills massless
         !! layers before the gradients are formed).
      real(wp) :: min_dz_for_n2 = 1.0_wp
         !! Minimum layer thickness, m, floored in the N² / drdz
         !! denominator so vanished layers don't spike the slope.
   end type ocean_slopes_config_t
   type :: ocean_gm_config_t
      !! `&ocean_gm_nml`: Gent-McWilliams thickness diffusion — eddy-induced
      !! bolus transport folded into continuity.  Reads the slopes slot.
      !! Default off ⇒ no bolus flux.
      logical :: enable = .false.
         !! Master switch.  Default off ⇒ `gm_compute_transports` no-ops.
         !! Requires `&ocean_slopes_nml enable` (configure-time error).
      real(wp) :: khth = 0.0_wp
         !! Thickness diffusivity KhTh, m^2/s (constant-fill for v1;
         !! production 1e2-1e3).  0 ⇒ no bolus transport even when enabled.
      real(wp) :: khth_max_cfl = 0.1_wp
         !! Fraction of the diffusive CFL the face KH may use.
      real(wp) :: khth_slope_max = 0.01_wp
         !! Slope magnitude above which the safe-streamfunction blend takes
         !! over (MOM6 `slope_max`).
   end type ocean_gm_config_t
   type :: ocean_redi_config_t
      !! `&ocean_redi_nml`: Redi continuous neutral (along-isopycnal) tracer
      !! diffusion — the rotated diffusion tensor.  Default off ⇒ no neutral flux.
      logical :: enable = .false.
         !! Master switch.  Default off ⇒ no-op ⇒ bit-identical.
      logical :: continuous = .true.
         !! Continuous variant (closed-form linear neutral surfaces, the
         !! production path).  `.false.` (discontinuous) DEFERRED, rejected
         !! at configure.
      real(wp) :: khtr = 0.0_wp
         !! Redi neutral diffusivity KhTr (m^2/s; production 1e2-1e3).  0 ⇒
         !! no neutral flux even when enabled.
   end type ocean_redi_config_t
   type :: ocean_varmix_config_t
      !! `&ocean_varmix_nml`: spatially-varying GM/Redi lateral-diffusivity
      !! coefficients = (background + Visbeck) × resolution_fn, clamped.
      !! Reads the slopes + wavespeed slots.  Default off ⇒ GM keeps constant khth.
      logical :: enable = .false.
         !! Master switch.  Default off.  Requires `&ocean_slopes_nml enable`
         !! + `&ocean_wavespeed_nml enable` (configure-time errors).
      logical :: use_visbeck = .false.
         !! Add the Visbeck/Eady `khth_slope_cff·L²·SN` baroclinicity term.
      logical :: resoln_scaled_khth = .false.
         !! Scale the assembled KhTh by the resolution function.
      logical :: resoln_scaled_khtr = .false.
         !! Scale the assembled KhTr by the resolution function.
      logical :: gill_equatorial_ld = .true.
         !! Gill (1982) equatorial-Ld convention (factor 2 in beta_dx2).
      logical :: interpolate_res_fn = .false.
         !! Interpolate the centre Res_fn to faces (`.true.`) vs interpolate
         !! cg1 to faces then recompute Res_fn (`.false.`, MOM6 default).
      integer :: kh_res_fn_power = 2
         !! Resolution-function power p (even).
      real(wp) :: kh_res_scale_coef = 1.0_wp
         !! Resolution-function alpha (the `(alpha·cg1)^p` denom coefficient).
      real(wp) :: khth = 0.0_wp
         !! Background thickness diffusivity KhTh (m²/s) the Visbeck term +
         !! Res_fn scale.
      real(wp) :: khtr = 0.0_wp
         !! Background tracer diffusivity KhTr (m²/s) for the future Redi.
      real(wp) :: khth_slope_cff = 0.0_wp
         !! Visbeck coefficient α_s for the KhTh chain.
      real(wp) :: khtr_slope_cff = 0.0_wp
         !! Visbeck coefficient for the KhTr chain.
      real(wp) :: khth_min = 0.0_wp
         !! Lower clamp on the assembled KhTh (m²/s).
      real(wp) :: khth_max = 0.0_wp
         !! Upper clamp on KhTh (m²/s); ≤ 0 ⇒ no cap.
      real(wp) :: khtr_min = 0.0_wp
         !! Lower clamp on KhTr (m²/s).
      real(wp) :: khtr_max = 0.0_wp
         !! Upper clamp on KhTr (m²/s); ≤ 0 ⇒ no cap.
      real(wp) :: visbeck_l_scale = 0.0_wp
         !! Visbeck length scale L (m); if < 0, |L|²·areaCu is used.
      real(wp) :: visbeck_max_slope = 0.0_wp
         !! S² limiter scale; ≤ 0 ⇒ no S² limit.
   end type ocean_varmix_config_t
   type :: ocean_meke_config_t
      !! `&ocean_meke_nml`: mesoscale eddy kinetic energy — a 2D prognostic
      !! eddy-energy field fed back as a thickness/tracer diffusivity into
      !! VarMix's KhTh/KhTr.  Reads `gm%gm_src`.  Default off ⇒ no-op.
      logical :: enable = .false.
         !! Master switch.  Default off.  Requires `&ocean_gm_nml
         !! enable=.true.` (configure-time error).
      real(wp) :: gmcoeff = -1.0_wp
         !! Efficiency of PE->MEKE conversion (nondim); < 0 ⇒ GM source off.
      real(wp) :: frcoeff = -1.0_wp
         !! Frictional mean->eddy conversion efficiency (nondim); < 0 ⇒ off
         !! (default).  When >= 0 the lateral-viscosity KE dissipation rate
         !! is sourced into MEKE.
      real(wp) :: bgsrc = 0.0_wp
         !! Background energy source (m^2/s^3).
      real(wp) :: damping = 0.0_wp
         !! Local depth-independent linear MEKE dissipation rate (1/s).
      real(wp) :: kh = -1.0_wp
         !! Background lateral diffusion of MEKE (m^2/s); < 0 ⇒ off.
      real(wp) :: k4 = -1.0_wp
         !! Background biharmonic diffusion of MEKE (m^4/s); < 0 ⇒ off.
      real(wp) :: khcoeff = 1.0_wp
         !! Scaling converting MEKE into Kh (nondim); <= 0 ⇒ closure off.
      real(wp) :: cd_scale = 0.0_wp
         !! Bottom/column eddy-velocity ratio (nondim); enters bottomFac2.
      real(wp) :: cb = 25.0_wp
         !! Coefficient in the gamma_bot (bottomFac2) expression (nondim).
      real(wp) :: ct = 50.0_wp
         !! Coefficient in the gamma_bt (barotrFac2) expression (nondim).
      real(wp) :: min_gamma2 = 1.0e-4_wp
         !! Floor on gamma_b^2 / gamma_t^2 (nondim).
      real(wp) :: uscale = 0.0_wp
         !! Background eddy velocity scale for bottom drag (m/s).
      real(wp) :: dtscale = 1.0_wp
         !! Scale factor accelerating MEKE time-stepping (nondim).
      real(wp) :: khth_fac = 0.0_wp
         !! Factor on the geom-mean kh added into VarMix KhTh (nondim).
         !! 0 (default) ⇒ feedback inert ⇒ bit-identical seam.
      real(wp) :: khtr_fac = 0.0_wp
         !! Factor on the geom-mean kh added into VarMix KhTr (nondim).
      logical :: backscatter = .false.
         !! Enable the MEKE → momentum harmonic backscatter (negative
         !! viscosity) energy return (Gap 2, v1).  Default off ⇒ bit-identical.
         !! Requires a flow-aware lateral closure (`lateral_closure /= none`)
         !! for the returned energy to reach the momentum tendency.
      real(wp) :: backscatter_visc_coeff_ku = 0.0_wp
         !! MOM6 `MEKE_VISCOSITY_COEFF_KU` — harmonic backscatter efficiency
         !! `Ku = coeff·sqrt(2·gamma_t²·E)·Lmix` (m²/s).  Subtracted (face-
         !! averaged, CFL-floored) from the resolved harmonic viscosity.
         !! 0 (default) ⇒ inert.  HARMONIC only in v1 (biharmonic `Au` +
         !! EBT/SQG vertical structure deferred).
      real(wp) :: khmeke_fac = 0.0_wp
         !! Factor relating meke%kh to MEKE's own lateral diffusivity (nondim).
      real(wp) :: advection_factor = 0.0_wp
         !! Barotropic-transport advection scaling (nondim); 0 ⇒ off (v1).
      real(wp) :: cdrag = 2.5e-3_wp
         !! Bottom drag coefficient for MEKE (nondim); enters drag_rate +
         !! Lfrict.  Set from `&ocean_bdrag_nml cdrag_side` at configure if
         !! that is > 0, else this default.
      logical :: use_bbl_drag = .false.
         !! Add the resolved bottom-boundary-layer eddy velocity to the MEKE
         !! bottom-drag rate: `drag_rate = rho0·i_mass·sqrt(cdrag²·(2·bf2·E +
         !! |u_bed|² + uscale²))` (MOM6 `drag_rate_visc` term, here the
         !! bed-layer speed).  Default `.false.` ⇒ the `|u_bed|²` term is 0 ⇒
         !! bit-identical to the prior MEKE drag.
      real(wp) :: alpha_deform = 0.0_wp
         !! Weight on the deformation length scale (nondim).
      real(wp) :: alpha_rhines = 0.0_wp
         !! Weight on the Rhines length scale (nondim); v1 default 0 ⇒ inert.
      real(wp) :: alpha_eady = 0.0_wp
         !! Weight on the Eady length scale (needs VarMix SN) (nondim).
      real(wp) :: alpha_frict = 0.0_wp
         !! Weight on the frictional-arrest length scale (nondim).
      real(wp) :: alpha_grid = 0.0_wp
         !! Weight on the grid length scale (nondim).
   end type ocean_meke_config_t
   type :: ocean_tidal_mixing_config_t
      !! St-Laurent/Simmons internal-tide interior mixing
      !! (`&ocean_tidal_mixing_nml`).  Bottom-intensified diapycnal
      !! diffusivity from the local dissipation of internal-tide energy
      !! over rough topography (Jayne & St Laurent 2001; St Laurent et
      !! al. 2002; Simmons et al. 2004).  An INTERIOR closure: its Kd is
      !! ADDED to the surface PBL schemes (KPP or EPBL) and to
      !! PP81/background/kappa-shear.  `enable = .false.` keeps the
      !! existing path bit-identical.  Knob table:
      !! `docs/generated_nml_knobs.md`.
      logical :: enable = .false.
         !! Master switch.  Requires `vmix%use_closure` + thermodynamics
         !! (validated at configure).
      real(wp) :: gamma = 0.3333_wp
         !! Local-dissipation fraction q (GAMMA_ITIDES).
      real(wp) :: mu = 0.2_wp
         !! Mixing efficiency Gamma_mix (MU_ITIDES).
      real(wp) :: zeta = 500.0_wp
         !! Bottom decay scale, m (INT_TIDE_DECAY_SCALE).
      real(wp) :: kd_max = 1.0e-2_wp
         !! Per-layer physical Kd cap, m^2/s; < 0 => no cap.
      real(wp) :: prandtl_tidal = 1.0_wp
         !! Kv = prandtl_tidal * Kd into the momentum solve.
      real(wp) :: min_zbot = 0.0_wp
         !! Mask off where column depth H < min_zbot, m.
      real(wp) :: e_uniform = 0.0_wp
         !! Uniform bottom energy input E, W m-2 (v1 prescribed field;
         !! default 0 => inert even when enabled).
      logical :: e_compute = .false.
         !! v1.1 state-dependent E = min(TKE_coef*N_bot, e_max).
      real(wp) :: kappa_itides = 6.2832e-4_wp
         !! Topographic wavenumber, m^-1 (v1.1 E recompute).
      real(wp) :: kappa_h2 = 1.0_wp
         !! KAPPA_H2_FACTOR (v1.1 E recompute).
      real(wp) :: utide = 0.0_wp
         !! RMS barotropic tidal velocity, m/s (v1.1 E recompute).
      real(wp) :: h2_rough = 0.0_wp
         !! Sub-grid topographic roughness variance <h^2>, m^2.
      real(wp) :: frac_rough = 0.1_wp
         !! Roughness clamp <h^2> <= (frac_rough*H)^2.
      real(wp) :: e_max = 1.0e3_wp
         !! TKE_itide_max cap on E, W m-2.
   end type ocean_tidal_mixing_config_t
   type :: ocean_porous_config_t
      !! Porous barriers (`&ocean_porous_nml`): subgrid sill/strait
      !! blocking of the C-grid face widths via the Adcroft (2013)
      !! three-parameter fit.  The layer-averaged OPEN-AREA fraction
      !! multiplies `dy_cu` / `dx_cv` in the layer mass transport
      !! (continuity-PPM + Coriolis advection), so a deep sill blocks the
      !! bottom layers while the surface layers stay fully open.
      !! `enable = .false.` keeps the existing path bit-identical.
      !! Knob table: `docs/generated_nml_knobs.md`.
      logical :: enable = .false.
         !! Master switch.  Requires `sim_type='ocean'` + a multilayer
         !! run (validated at configure).
      character(len=32) :: source = "resolved"
         !! Where the along-face min/max/mean topographic heights come
         !! from.  `"resolved"`: three samples of the RESOLVED
         !! bathymetry along the face (two corners + midpoint) — a
         !! documented PROXY, not true subgrid data.  `"file"`: an
         !! offline subgrid-bathymetry file (MOM6 `topog_edge.nc`),
         !! NOT implemented — fails loud pending the file-forcing
         !! backend.
      character(len=32) :: eta_interp = "max"
         !! Interface-height-at-velocity-point rule (MOM6
         !! `PORBAR_ETA_INTERP`): `"max"` (the higher, i.e. shallower, of
         !! the two adjacent interfaces — the default, and the LEAST
         !! blocking, since the open width increases with interface
         !! height), `"min"` (the most blocking), `"arithmetic"`,
         !! `"harmonic"`.
      real(wp) :: masking_depth = 0.0_wp
         !! Faces whose mean along-face depth is SHALLOWER than this
         !! (m, positive below the sea surface) are left fully open
         !! (MOM6 `PORBAR_MASKING_DEPTH`).  0 ⇒ apply everywhere the
         !! face is below sea level.
   end type ocean_porous_config_t

   type :: ocean_conv_config_t
      !! Brunt-Vaisala-triggered convective adjustment (`&ocean_conv_nml`,
      !! CVMix_conv-style).  Where the interior N^2 < n2_thresh
      !! (dense-over-light), raises kt -> max(kt, kd_conv) and
      !! kv -> max(kv, prandtl_conv*kd_conv), strictly below the active
      !! surface boundary layer (KPP/EPBL own that).  An INTERIOR
      !! closure CONTRIBUTOR -- it does not clip/floor on its own; it
      !! feeds `vmix_assemble`.  `enable = .false.` keeps the existing
      !! path bit-identical.  Knob table: `docs/generated_nml_knobs.md`.
      logical :: enable = .false.
         !! Master switch.  Requires `vmix%use_closure` +
         !! `thermo%enable_thermodynamics` (validated at configure).
      real(wp) :: kd_conv = 1.0_wp
         !! Convective tracer diffusivity (MOM6 `KD_CONV`), m^2/s.
      real(wp) :: prandtl_conv = 1.0_wp
         !! Kv_conv = prandtl_conv * kd_conv (MOM6 `PRANDTL_CONV`),
         !! nondim.
      real(wp) :: n2_thresh = 0.0_wp
         !! Trigger threshold on N^2 (MOM6 `BV_SQR_CONV`), s^-2.
   end type ocean_conv_config_t
   type :: ocean_ddiff_config_t
      !! Double diffusion (`&ocean_ddiff_nml`, CVMix_ddiff-style): salt
      !! fingering (Large et al. 1994) + diffusive convection
      !! (Marmorino-Caldwell 1976 / Kelley 1990).  Folded INTO the
      !! heat/salt split (not a kv/kt contributor) to give an asymmetric
      !! ks-vs-kt divergence.  Defaults are the CVMix defaults; `enable =
      !! .false.` keeps the path bit-identical.  Knob table:
      !! `docs/generated_nml_knobs.md`.
      logical :: enable = .false.
         !! Master switch (MOM6 `USE_CVMIX_DDIFF`).  Requires
         !! `thermo%enable_thermodynamics` (validated at configure).
      real(wp) :: strat_param_max = 2.55_wp
         !! R_rho salt-fingering cutoff (CVMix `STRAT_PARAM_MAX`), nondim.
      real(wp) :: kappa_ddiff_s = 1.0e-4_wp
         !! Leading salt-fingering salinity diffusivity K_f (CVMix
         !! `KAPPA_DDIFF_S`), m^2/s.  K_T = 0.7*K_S.
      real(wp) :: ddiff_exp1 = 1.0_wp
         !! Inner (bracket) fingering exponent (CVMix `DDIFF_EXP1`).
      real(wp) :: ddiff_exp2 = 3.0_wp
         !! Outer fingering exponent (CVMix `DDIFF_EXP2`); exp1=1,exp2=3
         !! is the Large et al. cubic.
      real(wp) :: param1 = 0.909_wp
         !! MC76 convection exterior coeff (CVMix `KAPPA_DDIFF_PARAM1`).
      real(wp) :: param2 = 4.6_wp
         !! MC76 convection middle coeff (CVMix `KAPPA_DDIFF_PARAM2`).
      real(wp) :: param3 = -0.54_wp
         !! MC76 convection interior coeff (CVMix `KAPPA_DDIFF_PARAM3`).
      real(wp) :: mol_diff = 1.5e-6_wp
         !! Molecular diffusivity scaling the convection branch (CVMix
         !! `MOL_DIFF`), m^2/s -- the molecular value, not a background eddy.
      logical :: use_k90 = .false.
         !! Convection form: `.false.` = Marmorino-Caldwell 1976 (default),
         !! `.true.` = Kelley 1990.
   end type ocean_ddiff_config_t
   type :: ocean_tides_config_t
      !! Equilibrium (astronomical) body-force tide (`&ocean_tides_nml`,
      !! capability C1).  Drives `-g grad(eta - eta_eq)` in the barotropic
      !! momentum solve from a set of harmonic constituents.  `enable =
      !! .false.` keeps the existing path bit-identical.  Requires a
      !! non-cartesian grid (needs lat/lon).
      logical :: enable = .false.
         !! Master switch (default off => bit-identical).
      logical :: use_sal = .false.
         !! Apply scalar self-attraction & loading (C2). Default off =>
         !! bit-identical.
      real(wp) :: beta_sal = 0.0_wp
         !! Scalar SAL factor beta (~0.085-0.12); `eta_sal = beta_sal*eta`.
      logical :: add_nodal = .false.
         !! Apply the 18.6-yr nodal f/u corrections (fixed at nodal_ref_date).
      character(len=64) :: constituents = "M2 S2 N2 K2 K1 O1 P1 Q1"
         !! Whitespace/comma-separated active constituent list.
      character(len=16) :: ref_date = "1900-01-01"
         !! Astronomical reference date "YYYY-MM-DD"; model t=0 == ref_date.
      character(len=16) :: nodal_ref_date = ""
         !! Nodal reference date; "" => ref_date.
   end type ocean_tides_config_t
   type :: ocean_psurf_config_t
      !! Atmospheric surface-pressure loading / inverse barometer
      !! (`&ocean_psurf_nml`, PR-17).  Folds `eta_ib = -p_surf/(rho0 g)`
      !! into the barotropic `eta_forcing` seam so an atmospheric high
      !! depresses SSH (~1 cm/hPa; Wunsch & Stammer 1997).  `enable =
      !! .false.` keeps the existing path bit-identical.  Split-solver only
      !! and mutually exclusive with `&ocean_bt_nml bt_halo > 0`.  Requires
      !! `&ocean_forcing_nml enable_components=.true.` (the `p_surf` field
      !! is allocated only with the component set).  ρ₀ is NOT a knob here —
      !! it is taken from `ocean_state%eos%rho0` (the single ρ₀ of record).
      logical :: enable = .false.
         !! Master switch (default off => bit-identical).
      logical :: in_eos = .false.
         !! Include the surface load in the EOS's **IN-SITU** pressure
         !! arguments, as the top-of-column pressure `p_top`
         !! (`multilayer_state_t%p_top`): a hydrostatic pressure that used
         !! to start at 0 Pa at the free surface starts at `p_top(i,j)`
         !! instead.  Default `.false.` => `p_top` stays the zero array and
         !! every EOS evaluation is bit-identical.  This is the E3 seam for
         !! ice-shelf cavities, where 1e6-2e7 Pa of ice load makes the p=0
         !! assumption a systematic ~4-5 kg/m^3 density error under a
         !! nonlinear EOS.
         !!
         !! IN-SITU ONLY.  It deliberately does NOT touch `ms%rho_layer`,
         !! which is a POTENTIAL density at the horizontally uniform
         !! `&ocean_eos_nml p_ref` — offsetting a potential-density
         !! reference per column would manufacture an along-layer density
         !! gradient (see that knob's docstring).  The N² builders inherit
         !! `rho_layer` and are therefore unaffected and self-consistent.
         !!
         !! v1 reaches exactly ONE consumer: the FV_WRIGHT Picard column
         !! sweep (`&ocean_pgf_nml form="fv_wright"`), the only in-situ EOS
         !! pressure in the dyn core today.  With any other PGF form and
         !! none of the closures below enabled the knob is INERT — a
         !! rank-0 warning says so rather than leaving it silent.  The
         !! closures that build their OWN surface-relative hydrostatic
         !! pressure have NOT been ported, so `validate_config` REFUSES
         !! `in_eos = .true.` together with any of them (EPBL,
         !! kappa-shear, tidal mixing, Redi, isopycnal slopes, sea ice,
         !! and the PGF in-layer reconstruction) rather than run a
         !! silently inconsistent pressure.  Requires `enable = .true.`
         !! (the load itself comes from `sf%p_surf`).
      real(wp) :: p_surf_const = 0.0_wp
         !! Uniform atmospheric surface pressure (Pa) seeded into
         !! `sf%p_surf_atm` at configure.  A UNIFORM load is provably inert
         !! (gauge invariance — only grad(p_surf) is physical), so a
         !! non-zero value with no file/override path emits a rank-0
         !! warning.  v1 fill path (the analogue of `&ocean_thermo_nml
         !! q_heat`).  The NetCDF forcing reader itself ships
         !! (`&ocean_dataovr_nml`); `p_surf` is simply not one of its six
         !! tags (tau_x, tau_y, heat, salt, evap, lprec) yet, so a
         !! file-driven load needs a `register_tag` entry, not new reader
         !! machinery.
   end type ocean_psurf_config_t
   type :: ocean_cavity_dyn_config_t
      !! Static ice-shelf cavity GEOMETRY (`&ocean_cavity_dyn_nml`,
      !! Phase 5.1).  A prescribed, time-constant ice draft `z_draft(i,j)`
      !! (m, positive DOWN — the depth of the ice base below `z = 0`) is
      !! laid over the bed and absorbed into the barotropic DATUM:
      !!
      !!     bt_H_ref = b - z_draft      (was: bt_H_ref = b)
      !!
      !! so the column starts with `bt_eta = sum(h_layer) - bt_H_ref = 0`
      !! under the shelf and every consumer of the water-column thickness
      !! `D = bt_H_ref + bt_eta` is correct without its own cavity branch.
      !! This is Losch (2008) §2.1's convention verbatim ("the
      !! 'sea-surface height' eta is the deviation from the 'reference'
      !! ice-shelf draft h"), not a divergence from it.
      !!
      !! The isostatic load `p_ice_ref = rho_ref*GRAVITY*z_draft` (Pa) is
      !! built at configure and assembled into the top-of-column pressure
      !!
      !!     ms%p_top = metrics%p_ice_ref + sf%p_surf
      !!
      !! whose consumers are the FV_MOM6 surface BC
      !! (`&ocean_pgf_nml p_top_in_bc`, REQUIRED unless the draft is
      !! uniform) and the in-situ EOS pressure (`&ocean_psurf_nml
      !! in_eos`).  The load is deliberately NOT added to `sf%p_surf`:
      !! the datum `bt_H_ref = b - z_draft` already carries its whole
      !! barotropic effect, and `eta_ib` is built from the assembled
      !! `sf%p_surf`, so only the load ANOMALY belongs on that seam.
      !!
      !! `enable = .false.` (default) keeps `z_draft` at its `(1,1)`
      !! placeholder, `bt_H_ref = b`, and every path bit-identical.
      !! Knob table: `docs/generated_nml_knobs.md`.
      logical :: enable = .false.
         !! Master switch.  Requires the ocean multilayer path, the split
         !! solver, `&ocean_pgf_nml form="fv_mom6"`, `vcoord_type` in
         !! {sigma, zstar} and a single rank; mutually exclusive with
         !! wet/dry, porous barriers, sea ice, `bt_halo > 0`, tidal SAL
         !! and `gfs_scale /= 1` (every one of those fails loud at
         !! configure, naming the knob and the reason).
      character(len=32) :: draft_config = "none"
         !! Analytic draft shape.  `"none"` (default): `z_draft = 0`
         !! everywhere — the identity, even with `enable = .true.`.
         !! `"flat"`: uniform `draft_depth` inside the shelf box
         !! `[draft_x0, draft_x1] x [draft_y0, draft_y1]`, 0 outside (the
         !! open ocean beyond the calving front at `draft_x1`).
         !! `"linear"`: `z_draft = draft_depth + draft_slope*(x - draft_x0)`
         !! inside the same box, clipped at 0 below.  `"file"`: a static
         !! 2-D NetCDF draft — NOT implemented (fails loud); the
         !! MPI-correct static-2-D reader is a later slice, and it is the
         !! only route to an ISOMIP+ draft, which has no analytic form
         !! (Asay-Davis et al. 2016 §3.1.1).
      character(len=32) :: draft_source = "draft"
         !! What the draft is prescribed FROM.  `"draft"` (default): the
         !! geometry above IS the ice-base depth.  `"thickness"`: the
         !! formula gives an ice THICKNESS, converted by the
         !! Boussinesq-isostatic (flotation) relation
         !! `z_draft = rho_ice*h_ice/rho_0`.  `"in_situ"` (true isostasy,
         !! `p_ice = g*integral(rho_hat)`) needs a per-column root find,
         !! does not admit exact discrete rest, and is deliberately NOT
         !! implemented (fails loud).
      real(wp) :: draft_depth = 0.0_wp
         !! Draft amplitude (m, positive down) — the uniform value for
         !! `"flat"`, the value at `draft_x0` for `"linear"`.  Under
         !! `draft_source = "thickness"` it is an ice THICKNESS instead.
      real(wp) :: draft_slope = 0.0_wp
         !! `"linear"` only: d(draft)/dx, dimensionless (m of draft per m
         !! of x).  Positive deepens the ice base toward +x.  Converted
         !! from metres to GRID units at the dispatch (metres on a
         !! Cartesian grid, degrees on spherical/curvilinear), the same
         !! way `&ocean_topo_nml slope_scale` is.
      real(wp) :: draft_x0 = -1.0e30_wp
         !! Western edge of the shelf box (m, GLOBAL physical coordinate),
         !! and the ANCHOR of the `"linear"` profile (`draft_depth` is the
         !! draft AT `draft_x0`), which is why `"linear"` requires a
         !! finite value here.  The default +/-1e30 on all four bounds is
         !! the "no limit on this side" sentinel: the shelf then covers
         !! the whole domain INCLUDING the ghost band, which is what a
         !! shelf that reaches a wall needs (a box stopping at x = 0 puts
         !! a phantom calving front one cell outside the west wall).
      real(wp) :: draft_x1 = 1.0e30_wp
         !! Eastern edge of the shelf box = the CALVING FRONT (m): beyond
         !! it the draft is 0 (open ocean).  Default: no eastern limit.
      real(wp) :: draft_y0 = -1.0e30_wp
         !! Southern edge of the shelf box (m).  Default: no limit.
      real(wp) :: draft_y1 = 1.0e30_wp
         !! Northern edge of the shelf box (m).  Default: no limit.
      character(len=256) :: draft_file = ""
         !! `draft_config="file"`: path to the NetCDF carrying the static
         !! ice draft.  The variable named by `draft_var` must be rank 3
         !! in FORTRAN storage order `(x, y, t)` — which is how a C or
         !! Python writer (and `ncdump`) spells `(nTime, ny, nx)` — with a
         !! time coordinate variable; RECORD 1 is read and the field is
         !! never re-read.  There is NO horizontal interpolation: the file
         !! must already be on the model grid (`nx x ny` physical cells),
         !! exactly as `bathymetry_file` and `&ocean_zinit_nml file`
         !! require.  Single rank only (fail-loud otherwise).
      character(len=64) :: draft_var = "iceDraft"
         !! `draft_config="file"`: name of the 2-D variable to read.  The
         !! default is the ISOMIP+ geometry file's own spelling
         !! (Asay-Davis et al. 2016 Sect. 3.3).
      character(len=16) :: draft_sign = "depth"
         !! `draft_config="file"`: the SIGN CONVENTION of the file values.
         !! There is no default that guesses from the data — the two
         !! conventions differ by the whole ice load, and a field that is
         !! partly open water (zeros) is indistinguishable by inspection.
         !!
         !!   * `"depth"` / `"positive_down"` (default) — the values ARE
         !!     the ice-base depth, `>= 0`, Roundabout's own convention.
         !!   * `"elevation"` / `"positive_up"` — the values are the
         !!     ice-base ELEVATION `z_d`, `<= 0` under a floating shelf.
         !!     Negated on load.  **This is what the ISOMIP+ geometry file
         !!     needs**: its `iceDraft` is "the elevation of the
         !!     ice-ocean interface (z_d)".
         !!
         !! A file in the wrong convention produces a negative depth and
         !! is caught fail-loud by the existing non-negativity check, not
         !! silently accepted.
      real(wp) :: h_min_cavity = 10.0_wp
         !! GROUNDING cutoff (m): a column whose water thickness
         !! `b - z_draft` is below this is LAND — it goes through the same
         !! `seed_wet_mask_impl` the bathymetry uses, so the static
         !! metric-zeroing land mask and the finite land-state hold for
         !! free.  Never a thin film of water under grounded ice.
         !! ISOMIP+ §3.1.5 leaves the choice to the modeller and notes
         !! ~40 m (two cells) for z-level models; sigma is less
         !! restricted, hence 10 m.
      real(wp) :: grounded_max_frac = 0.5_wp
         !! Sanity bound: if more than this fraction of the interior
         !! columns ground, configure fails loud rather than silently
         !! running a domain that is mostly land.
      real(wp) :: rho_ice = 918.0_wp
         !! Ice density (kg/m^3), consulted ONLY by
         !! `draft_source = "thickness"`.
   end type ocean_cavity_dyn_config_t

   type :: ocean_cavity_melt_config_t
      !! Ice-shelf basal-melt THERMODYNAMICS (`&ocean_cavity_melt_nml`,
      !! Phase 2b).  The three-equation interface of Holland & Jenkins
      !! (1999), solved once per thermo step on every ice-covered column
      !! and delivered to the ocean as two OWNED surface-flux components
      !! (`heat_cavity`, `salt_cavity`).
      !!
      !! GEOMETRY IS A DIFFERENT GROUP.  The draft, the cover mask and
      !! the barotropic datum are `&ocean_cavity_dyn_nml`'s, and this
      !! group REQUIRES it: without a draft there is no interface to melt
      !! and `cover_frac` is identically zero.  The split is the
      !! per-concern sub-namelist convention, and it is also the honest
      !! one — a datum-only cavity run is a legitimate configuration.
      !!
      !! VIRTUAL OR REAL MASS — `freshwater` picks (Phase 3).
      !! `"virtual"` (default, bit-identical) delivers the meltwater as a
      !! virtual salt flux at fixed column mass; `"mass"` adds the real
      !! Boussinesq volume to the top layer and lets the dilution happen
      !! by itself.  See that knob, and
      !! `docs/CAPABILITIES_AND_LIMITATIONS.md`.
      !!
      !! `enable = .false.` (default) ⇒ no slot arrays, no kernel, no
      !! component written, byte-identical.  Knob table:
      !! `docs/generated_nml_knobs.md`.
      logical :: enable = .false.
         !! Master switch.  Requires `&ocean_cavity_dyn_nml enable`,
         !! `&ocean_eos_nml tfreeze_set="isomip"` and
         !! `&ocean_forcing_nml enable_components`; mutually exclusive
         !! with atmospheric surface forcing (wind stress, surface
         !! restoring, shortwave penetration, the uniform scalar
         !! `q_heat`/`q_salt`) until the per-cell cover mask lands, and
         !! with sea ice.  Every one of those fails loud at configure,
         !! naming the knob and the follow-up.  `&ocean_psurf_nml`
         !! composes freely: the liquidus reads the ASSEMBLED
         !! `ms%p_top = p_ice_ref + sf%p_surf`, so an atmospheric load
         !! under the shelf depresses the freezing point with no extra
         !! wiring.
      character(len=32) :: exchange_law = "const_gamma"
         !! Turbulent exchange-velocity law.  `"const_gamma"` (default)
         !! is `gamma = Gamma*u*` — Jenkins, Nicholls & Corr (2010)
         !! eqs. (1),(2),(5) p. 2300 and the ISOMIP+ form.  `"hj99"` is
         !! Holland & Jenkins (1999) eqs. (14)-(18) p. 1792 (needs a
         !! non-zero Coriolis parameter under the cover).  `"yung25"` is
         !! Yung et al. (2025) "StratFeedback" eqs. (7)-(8) p. 5832.
         !! Every other name the kernel's enum reserves
         !! (`jenkins91`, `rosevear22`, `vt19`, `mk18`, `burchard22`,
         !! `jenkins21`) PARSES but is refused at configure naming
         !! `CAVITY_MELT_NOT_IMPLEMENTED` — a reserved law and a typo
         !! must stay distinguishable.
      real(wp) :: gamma_t = 2.2e-2_wp
         !! Dimensionless heat-transfer coefficient `Gamma_T` of
         !! `gamma_t = Gamma_T*u*`.  ISOMIP+ starting guess, Asay-Davis
         !! et al. (2016) §3.2.1 p. 2487 — **a starting guess, not a
         !! constant of nature**: the protocol has participants tune it,
         !! and Yung et al. (2026) Table 2 p. 2058 shows the twelve
         !! submissions spanning 0.011 to 0.2.  Re-derive it per vertical
         !! coordinate; a single value across a coordinate sweep makes
         !! the sweep measure its own tuning.
      real(wp) :: gamma_s = -1.0_wp
         !! Dimensionless salt-transfer coefficient `Gamma_S`.
         !! **Negative = unset ⇒ resolved to `gamma_t/35`**, the ISOMIP+
         !! ratio (Asay-Davis et al. (2016) Table 4 p. 2483, after
         !! Jenkins, Nicholls & Corr (2010) p. 2309: the ratio "should
         !! lie somewhere in the range 35-70.  Adopting a value at the
         !! lower end of this range...").  Zero and positive values are
         !! taken literally, so `gamma_s = 0` is refused as a range
         !! error rather than silently re-triggering the default.
      real(wp) :: cdrag_top = 2.5e-3_wp
         !! Top drag coefficient `C_D,top` entering the MELT friction
         !! velocity `u*^2 = C_D (U^2 + u_tide^2)`.  ISOMIP+ Table 4
         !! p. 2483.  The least constrained number in the subject: the
         !! literature spans 1.5e-3 (Holland & Jenkins 1999 Table 1) to
         !! 9.7e-3 (Jenkins et al. 2010 Table 2).  **This knob does not
         !! yet drive any momentum drag** — top drag is Phase 4; here it
         !! only scales `u*` for the exchange velocities.
      real(wp) :: u_tide = 1.0e-2_wp
         !! RMS tidal velocity (m/s) in the melt friction velocity —
         !! ISOMIP+ Table 4 p. 2483 and eq. (27) p. 2485, after Jenkins,
         !! Nicholls & Corr (2010) eq. (10) p. 2309.  TRAP, and the
         !! protocol says it outright (p. 2486): "The computation of top
         !! and bottom drag do not incorporate utidal" — it belongs to
         !! the melt `u*` ONLY.
      real(wp) :: ustar_min = 1.0e-4_wp
         !! Friction-velocity floor (m/s) — Yung et al. (2025) eq. (14)
         !! p. 5836, value from their Table 2 p. 5838.  It exists because
         !! "a friction velocity of zero (perhaps created by initialising
         !! the model at rest) will result in identically zero melt ...
         !! which would be inconsistent with the presence of heat
         !! available for melting".
      character(len=32) :: ice_conduction = "insulating"
         !! Ice-side heat conduction.  `"insulating"` (default) is
         !! `q_ice = 0`, which the ISOMIP+ protocol PRESCRIBES
         !! (Asay-Davis et al. (2016) Table 4 p. 2483 sets `kappa_i = 0`
         !! and p. 2485 instructs participants not to use the H&J99
         !! advection-diffusion scheme); `t_ice` is then unread.
         !! `"adv_diff"` is Holland & Jenkins (1999) eq. (31) p. 1794 in
         !! its melting asymptote, which collapses to
         !! `q_ice = m_mass*c_i*(T_b - T_ice)` and is zero on freezing.
         !! `"diffusive"` is RESERVED and refused: it changes the
         !! melt/freeze BRANCH logic, not just a coefficient.
      real(wp) :: t_ice = -25.0_wp
         !! Ice interior temperature (degC), read by
         !! `ice_conduction="adv_diff"` only — Holland & Jenkins (1999)
         !! Table 1 p. 1790 uses `T_S ~ -25`.  Under `"insulating"` the
         !! kernel substitutes exactly zero, so a stale value here cannot
         !! leak into an insulating run.
      real(wp) :: s_ice = 0.0_wp
         !! Ice salinity (g/kg), `>= 0`.  Zero is the ISOMIP+ value
         !! (Asay-Davis et al. (2016) Table 4 p. 2483).  It must stay
         !! strictly below the far-field salinity — that inequality is
         !! what the three-equation root bracketing rests on — and a
         !! column violating it is COUNTED and given zero melt, not
         !! guessed at.
      real(wp) :: far_field_depth = 10.0_wp
         !! Thickness (m) below the ice base over which the far-field
         !! `(T, S, u, v)` are thickness-averaged, with a partial last
         !! layer.  **Metres, deliberately, never "layer nz".**  The melt
         !! rate is roughly linear in the thermal driving it is handed,
         !! and how far from the ice that was sampled is the dominant
         !! resolution artefact in the subject (Gwyther et al. 2020;
         !! Burchard et al. (2022) Table 2 p. 15 — the all-bulk error
         !! GROWS under refinement; Yung et al. (2026) p. 2074).  No
         !! protocol prescribes a value; 10 m is this repository's
         !! default and it must be held FIXED across any
         !! vertical-coordinate comparison, or the comparison measures
         !! the sampling depth instead.
      character(len=16) :: freshwater = "virtual"
         !! How the meltwater reaches the ocean.
         !!
         !! `"virtual"` (**default ⇒ bit-identical**) — no mass moves.
         !! The dilution is emulated at FIXED column mass by the exact
         !! fixed-mass equivalent salt flux `-m*(S_far - s_ice)`
         !! (derivation in `rdb_ocean_cavity_flux`'s module docstring).
         !!
         !! `"mass"` — the meltwater is a REAL Boussinesq VOLUME source
         !! on the top layer, `dh = m*dt/rho_0`, and the salinity falls
         !! by dilution on its own.  The virtual salt flux is then NOT
         !! also applied to the tracer (that would double-count); it is
         !! RETAINED in the assembled `Q_salt` solely as the surface
         !! buoyancy forcing KPP/EPBL read for `B_0`, and removed again
         !! from salinity (and from the pseudo-salt mirror) by the same
         !! in-stage kernel that adds the volume.  The top layer's heat
         !! additionally gains the enthalpy of the added water,
         !! `m*c_w*T_b`.
         !!
         !! `"mass"` is refused (fail loud, naming the follow-up)
         !! together with dynamic wet/dry and with a windowed
         !! tracer-advection ratio > 1.
      character(len=32) :: volume_compensation = "none"
         !! What to do with the volume `freshwater="mass"` adds to a
         !! CLOSED domain.  Requires `freshwater="mass"`.
         !!
         !! `"none"` (default) — nothing; the domain fills up.  Correct
         !! for a short run and for a domain with an open boundary that
         !! can pass the volume out.
         !!
         !! `"uniform_open_ocean"` — each thermo step the
         !! domain-integrated melt volume is removed again, spread
         !! UNIFORMLY (per unit area) over the wet cells the ice does
         !! NOT cover, each parcel carrying that cell's own T and S so
         !! no concentration there is changed.  Tracked as a mass, salt
         !! and heat SINK in all three console budgets.  This is the
         !! sea-level compensation ISOMIP+ Sect. 3.1.3 allows for the
         !! closed Ocean3/4 domains; Ocean0-2 have a restoring sponge
         !! that does not remove volume, so without it their cavity
         !! fills at metres per year.
   end type ocean_cavity_melt_config_t

   type :: ocean_continuity_config_t
      real(wp) :: h_min = 1.0e-6_wp
         !! Floor used by the PPM positivity limiter (MOM6
         !! `GV%Angstrom_H`).  When `ppm_limit_pos = .true.` and a
         !! cell's centred thickness is at or below this floor, the
         !! PPM face reconstruction is forced to a constant (= upwind
         !! for that cell); above the floor the parabola's interior
         !! minimum is forced to `h_min` if it would dip below.
         !! Default 1e-6 m preserves "pure positivity" semantics.
         !! MOM6 typically uses 1–4 m (≈ `2·Angstrom_H`).  No-op when
         !! `ppm_limit_pos = .false.`.
      logical :: ppm_limit_pos = .false.
         !! MOM6 `PPM_limit_pos` analogue.  When `.true.`, the PPM
         !! face-thickness reconstruction in continuity adds a
         !! positivity-preserving limiter that prevents the parabola
         !! from dipping below `h_min` inside a cell.  At thin /
         !! vanishing layers (shelf-break, seamount, wet/dry boundary)
         !! the reconstruction collapses to a constant — equivalent
         !! to upwind for that cell.  Bounds mass flux through
         !! near-vanishing layers by the actual layer thickness.
      logical :: vol_cfl = .false.
         !! MOM6 `vol_CFL` analogue.  When `.false.` (default) the
         !! continuity-PPM face thickness is the downwind PPM EDGE
         !! value (the CFL→0 limit), bit-identical to the pre-knob
         !! behaviour.  When `.true.` the donor-side face thickness is
         !! the swept-volume integral of the reconstructed parabola
         !! (adds the missing O(CFL) term).  Fixes the dt-sensitive
         !! near-bed mass residual at steep shelf breaks.
      logical :: positive_definite = .false.
         !! Positive-definite continuity (MOM6-prevention + Roundabout-conservation).
         !! When `.true.` the split layer continuity guarantees every layer
         !! stays `>= h_lim` after each direction pass — the PPM reconstruction
         !! edges are floored at `2·h_lim` (P1, MOM6's positive-definite
         !! PPM) and (P2, later) the per-donor outfluxes are scaled down
         !! (never thickness-inflated ⇒ zero mass created).  `h_lim` is derived
         !! at setup: `angstrom_h` on VCOORD_LAGRANGIAN, else 0 (⇒ the floor is
         !! inert on non-Lagrangian coords).  Default `.false.` ⇒ untaken
         !! branches only ⇒ bit-identical.  Fail-loud composed with
         !! `&ocean_wetdry_nml enable` (that module owns its own barotropic
         !! limiter; composition deferred).
      logical :: renorm_consistent_flux = .true.
         !! Continuous flux model in the `uhbt`/`vhbt` Newton
         !! renormalisation (see `continuity_t%renorm_consistent_flux`):
         !! a layer whose upwind donor flips under the barotropic
         !! correction carries `(u0 + du)·h_face(new donor)` instead of the
         !! historical `flux0 + du·h_face(new donor)`, which jumps by
         !! `u0·(h_new − h_old)` at the flip and leaves the solve with no
         !! root when `uhbt` falls in the gap — a wrong-sign layer
         !! transport, an O(η) mismatch between the layer and barotropic
         !! free surfaces at every thickness jump (a sigma layer over a
         !! bathymetric step), and an exponentially pumped barotropic
         !! grid-scale mode under `pred_corr` (finding B of the
         !! vertical-coordinate stability matrix).  Newton is bracketed by
         !! bisection (MOM6 `zonal_flux_adjust`).  Default `.true.`
         !! (MOM6 behaviour, maintainer decision 2026-09-22); it is
         !! bit-identical to the historical model on every face where no
         !! donor flips.  `.false.` restores the historical discontinuous
         !! model.  Ignored by the wet/dry single-step form.
   end type ocean_continuity_config_t
   type :: ocean_isopycnal_config_t
      !! `&ocean_isopycnal_nml`: grounding-stability controls for the
      !! Lagrangian (isopycnal-class, remap-free) vertical coordinate.
      !! Every knob here engages ONLY under sim_type='ocean' +
      !! VCOORD_LAGRANGIAN, so the whole group is bit-identical on every other
      !! vertical coordinate.  All are default-OFF except
      !! `pgf_skip_nonoverlap`, which is a correctness fix and defaults ON.
      real(wp) :: angstrom_h = 0.0_wp
         !! Phase 1: minimum-thickness floor (m) on the Lagrangian
         !! continuity h-update — MOM6 GV%Angstrom_H analogue used as the
         !! max(h_new, Angstrom_H) clamp. 0.0 ⇒ off ⇒ bit-identical
         !! (the floor loop is skipped). Physical floor (D4 taxonomy), NOT
         !! H_DIV_EPS (too small to bound 1/h) and NOT hardwired H_VANISHED
         !! (that is the skip/merge marker). Recommended run value 1e-3..1e-2 m.
         !! NOT strictly conservative: injects <= angstrom_h*areaT per floored
         !! layer-cell (zero on non-grounding steps). Thermo-on caveat (R7):
         !! the floor lifts h but not the companion hTr, so on a floored layer
         !! the implied Tr=hTr/h shifts — harmless for the adiabatic
         !! (enable_thermodynamics=.false.) isopycnal config; thermo-on
         !! isopycnal correctness is OUT OF SCOPE for v1.
      logical :: reset_vanished_u = .false.
         !! Phase 2: zero the face velocity of a layer vanished on BOTH
         !! adjacent cells (a massless layer carries no independent momentum).
      logical :: cfl_ignore_vanished = .false.
         !! Phase 3: exclude vanished layers from the console MaxCFL / panic /
         !! CFL truncation so a thin-layer velocity spike cannot falsely abort.
      logical :: conservative_floor = .false.
         !! Conservative minimum-thickness mode.  When `.true.` (requires
         !! `angstrom_h > 0` and VCOORD_LAGRANGIAN — fail-loud otherwise) the
         !! non-conservative `max(h, angstrom_h)` injection is REPLACED by a
         !! per-column conservative borrow: a sub-floor layer's deficit is taken
         !! from the surplus layers of the same column, so total thickness,
         !! momentum, and every tracer mass are preserved to round-off and the
         !! result is independent of the floor value.  Fires only where a layer
         !! is below floor; a strict no-op elsewhere (isopycnals not pinned).
         !! Default `.false.` ⇒ the legacy injecting floor ⇒ bit-identical.
      logical :: pgf_skip_nonoverlap = .true.
         !! Grounded-layer PGF gate.  **Default ON** — unlike every other knob
         !! in this group this one is a CORRECTNESS fix, not an opt-in, and it
         !! is still bit-identical everywhere except VCOORD_LAGRANGIAN (the
         !! only coordinate it is applied under).
         !!
         !! The FV/Montgomery face PGF is a two-point Jacobian,
         !! `-(1/ρ0)·[Δp_centre + g·ρ_layer·Δz_centre]/dx`.  Those two terms
         !! cancel AT REST only while the two abutting layer centres sit in a
         !! common z-interval whose ambient density IS `ρ_layer`.  In an
         !! isopycnal column a layer grounded against sloping topography is
         !! squeezed onto the `angstrom_h` floor on the shallow side while
         !! staying massive on the deep side; the two centres are then hundreds
         !! of metres apart, the interval between them spans OTHER density
         !! classes, and the residual `g·(ρ_layer − ρ̄_ambient)·∂z/∂x` is a
         !! spurious acceleration on a motionless ocean.  Nothing bounds it, so
         !! it spins a basin-scale current up out of nothing (the
         !! `seamount_conservative_floor` quiescent case reached 6.7 cm/s).
         !!
         !! When `.true.` the face PGF is ZEROED wherever the layer's z-extents
         !! in the two abutting columns do not overlap — i.e. exactly where the
         !! layer has wedged out against the bed and there is no common depth to
         !! difference the pressure across.  Overlapping (physically shared)
         !! layers are untouched, so the isopycnal interior is unchanged; a
         !! wedged-out layer can still be re-wetted by continuity's own upwind
         !! flux and by the barotropic correction.  Set `.false.` to recover the
         !! pre-fix behaviour for comparison.
         !!
         !! Covers the `mont`, `fv_lite`, `fv_wright` and `fv_mom6` PGF forms.
         !! The first three build the layer-centre depth in their face passes
         !! anyway; `fv_mom6` (layer-integrated FV-Bouss) does not, so the
         !! buffer is allocated and filled there for this gate ALONE, which is
         !! why the decision has to be taken before the PGF slot is
         !! initialised.  The defect is the geometry, not the quadrature:
         !! `fv_mom6` reproduces the same at-rest spin-up (En 2.17e-03 at day 2
         !! on `seamount_conservative_floor`, → 1.9e-27 gated), and so does
         !! `mont` (En 2.242e-03 ungated, → 7.65e-31 gated) even though its
         !! Montgomery recursion is a third kind of face expression again —
         !! a grounded layer puts the two columns' interface heights hundreds
         !! of metres apart, so `M` stops being horizontally uniform at rest.
         !! `gprime` differences interface positions directly and has nothing
         !! to gate — it is the one form that warns.
      logical :: check_h_positive = .false.
         !! DEBUG GUARD (vcoord-agnostic despite living in this group): abort on
         !! the FIRST negative `h_layer`, naming the pipeline stage that wrote
         !! it plus the offending `(i,j,k)` and that column's full profile.
         !! Without it a negative thickness is only visible later and elsewhere
         !! — as a NaN in the console stats, or as an hourly diagnostic minimum
         !! — by which point the producing kernel is unidentifiable.
         !! Costs one device-side min-reduction over `h_layer` per instrumented
         !! stage (~6% of an outer step at 600x600x50); no H<-D traffic on the
         !! healthy path. Default `.false.` ⇒ not called ⇒ bit-identical.
   end type ocean_isopycnal_config_t
   type :: ocean_topo_config_t
      !! Basin geometry + surface-forcing setup.  Despite the name,
      !! this carries forcing knobs (wind_config / taux_magnitude)
      !! and Coriolis tilt (beta / y_ref) alongside the bathymetry
      !! profile, because all of them live in the same logical
      !! "set up the box and what's pushing it" slot at init — the
      !! `&ocean_topo_nml` block is read once before the IC is
      !! seeded.  This grouping is settled (reviewed 2026-05-27);
      !! the `form` knob in `&ocean_coriolis_nml` is the *scheme*
      !! selector (Sadourny vs HK), distinct from the `beta`/`y_ref`
      !! *parameters* here.
      character(len=16) :: topo_config = "flat"
         !! Bathymetry profile selector: "flat", "spoon", etc.
      real(wp) :: max_depth = 2000.0_wp
         !! Basin maximum depth (m).  Used by spoon bathymetry and
         !! as the uniform depth for `topo_config="flat"`.
      real(wp) :: edge_depth = 100.0_wp
         !! Edge depth (m) for spoon bathymetry.  MOM6 default 100 m.
      real(wp) :: slope_scale = 400000.0_wp
         !! Exponential decay scale (m) for spoon bathymetry.  MOM6
         !! default TOPOG_SLOPE_SCALE = 400 km.
      real(wp) :: nl_continent_amp = 1.0_wp
         !! Continent amplitude for topo_config="neverworld2": scales
         !! the continent/ridge terms of the Neverworld2 basin (1.0 =
         !! full continents, 0.0 = aquaplanet with the southern channel
         !! only).  MOM6 NL_CONTINENT_AMP default 1.0.
      real(wp) :: nl_roughness_amp = 0.05_wp
         !! Roughness amplitude for topo_config="neverworld2": amplitude
         !! of the wavy bathymetry signal.  MOM6 NL_ROUGHNESS_AMP default
         !! 0.05.
      real(wp) :: nl_min_depth = 500.0_wp
         !! Minimum-depth floor (m) for topo_config="neverworld2"
         !! (MOM6 MINIMUM_DEPTH analogue).  The fractional-depth formula
         !! drives the depth to 0 (land) at the continent/wall cells; the
         !! C-grid dyn-core has no robust wet/dry path yet, so true land +
         !! the resulting sub-metre surface layers blow up.  Flooring to
         !! `nl_min_depth` keeps the basin all-wet (continents become
         !! shallow steering shelves) and stable.  Lower = stronger
         !! topographic steering but thinner layers (nl_min_depth/nz).
      character(len=16) :: wind_config = "constant"
         !! Surface wind-stress dispatch: "constant" (uniform tau_x/y
         !! from `wind_stress_x/y`) or "2gyre" (MOM6 sinusoidal).
      real(wp) :: taux_magnitude = 0.1_wp
         !! Peak zonal wind stress (Pa) for wind_config = "2gyre".
         !! MOM6 TAUX_MAGNITUDE default.
      real(wp) :: coriolis_beta = 0.0_wp
         !! Meridional gradient of f (s⁻¹ m⁻¹).  0 = f-plane (uses
         !! `cfg%coriolis_f` as a uniform value); non-zero engages
         !! `set_beta_plane` with f₀ = coriolis_f.
      real(wp) :: coriolis_y_ref = 0.0_wp
         !! Reference y-coordinate (m) at which f = coriolis_f under
         !! the beta-plane: f(y) = coriolis_f + beta*(y - y_ref).
      real(wp) :: x_origin = 0.0_wp
         !! Absolute x-coordinate (m) of the domain's WEST edge, for the
         !! formula bathymetries whose published formula is written in an
         !! absolute coordinate frame rather than a domain-relative one.
         !! Read ONLY by `topo_config="isomip_plus"`, whose bedrock
         !! polynomial is a function of the MISMIP+ `x` that runs from the
         !! ice divide at 0, while the ISOMIP+ *ocean* box starts at
         !! `x = 320 km` (Asay-Davis et al. 2016, Table 3 `x0`).  Default
         !! `0` ⇒ the model's own `x = 0` west edge ⇒ bit-identical for
         !! every other `topo_config`.
   end type ocean_topo_config_t
   type :: ocean_ic_config_t
      !! Initial-condition overlay + EOS reference state.  Drives
      !! the analytical configurations (eady-front, geostrophic-
      !! adjustment, …) plus the linear-EOS reference values that
      !! the rest of the ocean dyn step reads.
      character(len=32) :: ic_config = ""
         !! IC overlay tag.  "" (default) keeps the analytical T(z)
         !! + zero-velocity IC.  "eady" runs the Eady-front overlay;
         !! "geostrophic_adjustment" drops a Gaussian SSH bump on an
         !! f-plane; "baroclinic_jet" runs the two-layer reduced-gravity
         !! baroclinic-instability jet (SIM_DETAILS.md §5).
      real(wp) :: alpha_T = 1.7e-4_wp
         !! Linear-EOS thermal-expansion coefficient, **DIMENSIONAL**
         !! (kg/m³ per °C) — see `beta_S` for the conversion from the
         !! fractional 1/°C coefficient most protocols quote.
      real(wp) :: beta_S = 7.6e-4_wp
         !! Linear-EOS haline contraction coefficient, **DIMENSIONAL**
         !! (kg/m³ per PSU).
         !!
         !! UNITS TRAP.  Roundabout's linear EOS is written as the
         !! DENSITY-ANOMALY form
         !!
         !!   rho = rho_0 + beta_S·(S − S_ref) − alpha_T·(T − T_ref)
         !!
         !! so `alpha_T`/`beta_S` carry kg/m³ per unit T/S.  Most
         !! protocols (ISOMIP+, Asay-Davis et al. 2016 among them) quote
         !! the FRACTIONAL coefficients of the equivalent form
         !!
         !!   rho = rho_0·(1 − alpha·(T − T_ref) + beta·(S − S_ref))
         !!
         !! with alpha in 1/°C and beta in 1/PSU.  Convert by
         !! multiplying through by `rho_0`:
         !!
         !!   alpha_T = rho_0 · alpha      beta_S = rho_0 · beta
         !!
         !! e.g. ISOMIP+ (alpha = 3.733e-5 1/°C, beta = 7.843e-4 1/PSU,
         !! rho_0 = 1027.51) becomes `alpha_T = 3.8356948e-2`,
         !! `beta_S = 8.0587609e-1`.  Feeding the fractional numbers
         !! straight in under-states the density response ~1000×, which
         !! looks like a plausible but far too weakly stratified run.
      real(wp) :: T_ref = 10.0_wp
         !! Linear-EOS reference temperature (°C) — the T at which the
         !! thermal anomaly term vanishes.
      real(wp) :: S_ref = 35.0_wp
         !! Linear-EOS reference salinity (PSU) — the S at which the
         !! haline anomaly term vanishes.
      real(wp) :: rho_0 = 1035.0_wp
         !! Reference density (kg/m³) for the linear EOS and Boussinesq
         !! PGF — the density at `(T_ref, S_ref)`.
      real(wp) :: layer_rho_init(MAX_OCEAN_LAYER_RHO_INIT) = -1.0_wp
         !! Per-layer initial density (kg/m³), `k=1` bed → `k=nz`
         !! surface, mirroring MOM6's `COORD_CONFIG="gprime"` IC.  Any
         !! positive value bypasses the EOS-derived init and writes the
         !! user-specified densities straight into `ms%rho_layer` at the
         !! end of `ocean_state_seed_from_cfg`.  Count of non-sentinel
         !! (>= 0) entries must equal `nz_layers`.  Default `-1.0`
         !! sentinel ⇒ existing init path (bit-identical).
      real(wp) :: rho_lightest = -1.0_wp
         !! Linear density-range IC — MOM6 `COORD_CONFIG="linear"`
         !! analogue.  When `>= 0`, generates `nz_layers` linearly-spaced
         !! layer densities spanning `[rho_lightest, rho_lightest +
         !! rho_range]` and writes them into `ms%rho_layer` (bed `k=1`
         !! heaviest → surface `k=nz` lightest), reusing the
         !! `layer_rho_init` write path — so the column scales by just
         !! bumping `nz_layers`, no hand-listed densities.  Mutually
         !! exclusive with `layer_rho_init`.  Default `-1.0` sentinel ⇒
         !! off (bit-identical).
      real(wp) :: rho_range = 2.0_wp
         !! Linear density-range IC: total top-to-bottom density contrast
         !! (kg/m³), MOM6 `DENSITY_RANGE` (default 2.0).  Read only when
         !! `rho_lightest >= 0`.
      real(wp) :: eady_dT_dy = -2.0e-5_wp
         !! Meridional gradient of T (°C/m) for the Eady IC.
      real(wp) :: eady_dT_dz = 0.01_wp
         !! Vertical stratification (°C/m) for the Eady IC.
      real(wp) :: eady_T_ref = 10.0_wp
         !! Reference temperature (°C) for the Eady IC at z=0, y=y_mid.
      real(wp) :: eady_pert_amp = 1.0e-3_wp
         !! Amplitude (°C) of the random temperature perturbation
         !! seeded into the Eady IC to break symmetry.
      integer :: eady_pert_seed = 12345
         !! RNG seed for the Eady IC perturbation.
      real(wp) :: ga_eta_amp = 1.0_wp
         !! Geostrophic-adjustment IC: amplitude (m) of the Gaussian
         !! SSH bump.
      real(wp) :: ga_length_scale = 50000.0_wp
         !! Geostrophic-adjustment IC: e-folding scale (m) of the
         !! Gaussian SSH bump.
      real(wp) :: ga_x_center = -1.0_wp
         !! x-coordinate (m) of the bump centre.  Negative ⇒ auto
         !! (basin midpoint).
      real(wp) :: ga_y_center = -1.0_wp
         !! y-coordinate (m) of the bump centre.  Negative ⇒ auto
         !! (basin midpoint).
      real(wp) :: jet_half_width = 40000.0_wp
         !! Baroclinic-jet IC (`ic_config="baroclinic_jet"`): tanh jet
         !! half-width L (m).  Spec "tuned" default 40 km; the Julia-
         !! faithful case uses 100 km.  Read only for the baroclinic-jet IC.
      real(wp) :: interface_amp = 200.0_wp
         !! Baroclinic-jet IC: interface displacement amplitude Δξ (m).
         !! Spec "tuned" default 200 m; faithful uses 100 m.
      real(wp) :: pert_amp_frac = 0.2_wp
         !! Baroclinic-jet IC: front-localised meander amplitude as a
         !! fraction of Δξ (`A_pert = pert_amp_frac·Δξ`).  Spec default 0.2.
      integer :: pert_nx = 3
         !! Baroclinic-jet IC: zonal perturbation wavenumber (integer number
         !! of wavelengths across the periodic domain).  Integer so the
         !! analytic x-continuation wraps exactly (`k_x·L_x = 2π·pert_nx`).
         !! Spec default 3.
      real(wp) :: upper_layer_rest = 500.0_wp
         !! Baroclinic-jet IC: rest thickness of the upper (surface, k=2)
         !! layer H₁ (m).  The lower (bed, k=1) rest thickness is
         !! `max_depth − upper_layer_rest`.  Spec §4 default 500 m
         !! (H₂ = 1500 m over the 2000 m column).
   end type ocean_ic_config_t
   type :: ocean_zinit_config_t
      !! Z-level T/S initial-condition overlay (`&ocean_zinit_nml`,
      !! capability A2).  Reads pre-regridded T/S from a model-grid
      !! NetCDF and interpolates linearly in depth onto the seeded
      !! layer-centre depths, overwriting any analytical T/S.  All
      !! defaults preserve bit-identity: `enable = .false.` is a no-op.
      !! Design + MOM6 divergences: `local_archive/specs/a2_zinit_spec.md`.
      logical :: enable = .false.
         !! Master switch.  Default `.false.` keeps the analytical IC.
         !! Requires `RDB_ENABLE_NETCDF=ON` at build time (the overlay
         !! lives in the NetCDF-gated `rdb_ocean_z_init`, including the
         !! `source = "linear"` path, which opens nothing).
      character(len=32) :: source = "file"
         !! Where the T(z)/S(z) profile comes from.
         !!
         !! `"file"` (default) reads the pre-regridded NetCDF named by
         !! `file` and interpolates linearly in depth.
         !!
         !! `"linear"` evaluates the ANALYTIC affine profiles
         !! `T(z) = lin_t_ref + lin_dt_dz*z`, `S(z) = lin_s_ref +
         !! lin_ds_dz*z` at every layer centre's true geopotential depth
         !! — no file, no interpolation, defined on ghosts and land too.
         !! This is the profile an idealised sloping-lid cavity (or any
         !! ISOMIP+-style case) needs: the `&tracer_nml T_init_surface` /
         !! `T_init_bottom` family is linear in LAYER INDEX, so under a
         !! terrain-following coordinate with a tilted lid it tilts the
         !! isopycnals with the coordinate and the column is NOT at rest.
      character(len=256) :: file = ""
         !! Path to the model-grid T/S NetCDF (dims x/y/z; vars
         !! temp/salt/z_src with the documented name-fallbacks).
      character(len=32) :: t_var = ""
         !! Temperature variable-name override; "" tries temp/T/temperature.
      character(len=32) :: s_var = ""
         !! Salinity variable-name override; "" tries salt/S/salinity.
      character(len=32) :: z_var = ""
         !! Source-axis variable-name override; "" tries z_src/z/depth/lev.
      real(wp) :: land_fill_t = 10.0_wp
         !! Fallback temperature (°C) written to dry (wet_mask <= 0) columns.
      real(wp) :: land_fill_s = 35.0_wp
         !! Fallback salinity (PSU) written to dry (wet_mask <= 0) columns.
         !! `source = "file"` only — the analytic path has a value
         !! everywhere and uses it.
      real(wp) :: lin_t_ref = 0.0_wp
         !! `source = "linear"`: temperature (degC) at the `z = 0` datum.
      real(wp) :: lin_dt_dz = 0.0_wp
         !! `source = "linear"`: `dT/dz` (degC/m) with **z positive UP**
         !! — the same convention as `&ocean_ic_nml eady_dT_dz`.  A
         !! thermally STABLE column has `lin_dt_dz > 0` (warm on top).
         !! `0` (the default) is a uniform column.
      real(wp) :: lin_s_ref = 35.0_wp
         !! `source = "linear"`: salinity (PSU) at the `z = 0` datum.
      real(wp) :: lin_ds_dz = 0.0_wp
         !! `source = "linear"`: `dS/dz` (PSU/m) with **z positive UP**.
         !! Note the polarity is the INVERSE of temperature's: salty
         !! water belongs at the bed, so a halinely STABLE column has
         !! `lin_ds_dz < 0`.  `0` (the default) is a uniform column.
   end type ocean_zinit_config_t
   type :: ocean_data_config_t
      !! `&ocean_data_nml`: the shared time-varying NetCDF input reader
      !! (`rdb_ocean_data_input`, PR-14).  Deliberately two knobs — this
      !! group carries NO per-field entries.  Registration is
      !! programmatic: each consumer (surface forcing, OBC segments,
      !! sponge targets, tidal-mixing maps, ...) calls
      !! `ocean_data_input_register_2d/_3d` from its OWN namelist group
      !! at setup and gets back an opaque `id`.  A field never appears
      !! here.  Zero registered fields on every shipped namelist today
      !! ⇒ `ocean_data_input_update_all` is a no-op ⇒ bit-identical.
      integer :: max_fields = 16
         !! Size of the reader's field registry.  A registration beyond
         !! this is a fail-loud abort (`max_fields >= 1` is enforced by
         !! the schema `min=1`).
      logical :: verbose = .false.
         !! Log every bracket advance (field, record pair, weight, model
         !! time) via `logger%info`.  Default off — quiet per-step path.
   end type ocean_data_config_t
   type :: dataovr_entry_config_t
      !! One file-driven surface-forcing tag in `&ocean_dataovr_nml`.
      !! A blank `file` means "this tag is not file-driven" — the slot
      !! keeps whatever the configure-time scalar/formula path seeded.
      character(len=256) :: file = ""
         !! Path to the pre-regridded `field(x, y, time)` NetCDF.  Blank
         !! ⇒ tag inactive.
      character(len=64) :: var = ""
         !! Variable name inside `file`.  Required when `file` is set —
         !! a blank `var` with a non-blank `file` is a fail-loud setup
         !! abort, never a guessed name.
      real(wp) :: scale = 1.0_wp
         !! Multiplied into the slab once at read (unit conversion).
      real(wp) :: add_offset = 0.0_wp
         !! Added after `scale`, once at read.
   end type dataovr_entry_config_t
   type :: ocean_dataovr_config_t
      !! `&ocean_dataovr_nml`: file-backed surface forcing (PR-15, the A3
      !! keystone consumer).  Binds time-varying NetCDF fields, through
      !! the shared PR-14 reader (`&ocean_data_nml`), onto the surface
      !! stress and surface-flux slots.  `enable = .false.` registers
      !! nothing ⇒ `ocean_data_input_update_all` stays a no-op ⇒
      !! bit-identical.
      !!
      !! The heat/freshwater tags write the surface-flux COMPONENT
      !! fields, not `Q_heat`/`Q_salt` directly, because
      !! `ocean_surface_flux_assemble` fully overwrites those two from
      !! the components every thermo step — a direct write would be
      !! silently clobbered.  They therefore require
      !! `&ocean_forcing_nml enable_components = .true.` (fail-loud
      !! otherwise).  Sign conventions are inherited from the component
      !! slots: `evap <= 0`, `lprec >= 0`, both kg/m^2/s.
      logical :: enable = .false.
         !! Master switch.  Default off ⇒ no registration, bit-identical.
         !! Requires `RDB_ENABLE_NETCDF=ON` at build time.
      character(len=16) :: time_mode = "linear"
         !! Shared time-axis mode for every active tag: `linear`
         !! (interannual; out-of-range is an abort unless `oor_clamp`),
         !! `cyclic` (climatology wrapped through `cycle_period`), or
         !! `static` (record 1, read once).  Per-tag overrides are a
         !! documented v2 gap, not an oversight.
      real(wp) :: cycle_period = 0.0_wp
         !! Climatology period (s).  Required > 0 when
         !! `time_mode = "cyclic"`; ignored otherwise.
      real(wp) :: t_offset = 0.0_wp
         !! Added to the model time before the file-axis lookup (s) —
         !! shifts the model epoch onto the file's epoch.
      logical :: oor_clamp = .false.
         !! `.true.` clamps a query outside the file time axis to the
         !! end record (one warning); default `.false.` aborts.
      type(dataovr_entry_config_t) :: tau_x
         !! Zonal wind stress (Pa) -> `surface_stress%tau_x` (east faces).
      type(dataovr_entry_config_t) :: tau_y
         !! Meridional wind stress (Pa) -> `surface_stress%tau_y` (north faces).
      type(dataovr_entry_config_t) :: heat
         !! Net surface heat flux (W/m^2, positive down) ->
         !! `surface_flux%heat_added`.
      type(dataovr_entry_config_t) :: evap
         !! Evaporative mass flux (kg/m^2/s, <= 0) -> `surface_flux%evap`.
      type(dataovr_entry_config_t) :: lprec
         !! Liquid precipitation (kg/m^2/s, >= 0) -> `surface_flux%lprec`.
      type(dataovr_entry_config_t) :: salt
         !! Surface salt flux (kg salt/m^2/s, positive salinifies) ->
         !! `surface_flux%salt_flux`.
   end type ocean_dataovr_config_t
   type :: ocean_diag_config_t
      logical :: enabled = .true.
         !! Enable per-step diag-manager hook in driver_run_ocean.
      character(len=256) :: filename = "ocean_diag"
         !! Output file basename — `output_rank_filename` appends the
         !! per-rank suffix.  Final path:
         !! `<output_dir>/<filename>_rank_NNNNNN.nc`.
      real(wp) :: dt_out = 3600.0_wp
         !! Cadence (s) at which the diag manager fires every variable's
         !! time op (INSTANT writes, MEAN flushes the accumulator, etc.).
      character(len=16) :: vgrid = "layer"
         !! Default output vertical grid for layer-shaped diags: "layer"
         !! (native `nz_ml` layers, default), "z_fixed" (conservative remap
         !! to `z_levels`), "sigma" (terrain-following, `sigma_levels`),
         !! "zstar" (SSH-tracking, `zstar_levels`), or "density" (isopycnal
         !! bins, `rho_levels` — requires `n_rho_levels > 0`, strictly
         !! increasing).  Per-diagnostic overrides via the `diags` list
         !! `:layer/z/sigma/zstar/density` attribute.
      real(wp) :: z_levels(MAX_OCEAN_DIAG_Z_LEVELS) = -1.0_wp
         !! Output z-levels (m, positive downward) when vgrid = "z_fixed".
      integer :: n_z_levels = 0
         !! Number of entries used in z_levels (0 = none).
      real(wp) :: sigma_levels(MAX_OCEAN_DIAG_Z_LEVELS) = -1.0_wp
         !! Output sigma levels — cumulative fractions (0..1, shallow->deep)
         !! when sigma output is used.  Empty => auto-generate `nz_ml`
         !! uniform fractions.
      integer :: n_sigma_levels = 0
         !! Number of entries used in sigma_levels (0 = auto-uniform).
      real(wp) :: zstar_levels(MAX_OCEAN_DIAG_Z_LEVELS) = -1.0_wp
         !! Output z* reference interface depths (m, positive-down,
         !! shallow->deep; deepest = H_ref) when z* output is used.  Empty =>
         !! auto-generate `nz_ml` uniform depths (== sigma; supply a
         !! non-uniform reference for z* to differ).
      integer :: n_zstar_levels = 0
         !! Number of entries used in zstar_levels (0 = auto-uniform).
      real(wp) :: rho_levels(MAX_OCEAN_DIAG_Z_LEVELS) = -1.0_wp
         !! Output target potential-density bin edges (kg/m^3, strictly
         !! increasing, light->dense) when vgrid = "density".  No auto-fill —
         !! `validate_config` requires `n_rho_levels > 0` whenever the
         !! density vgrid is selected (globally or per-diagnostic).
      integer :: n_rho_levels = 0
         !! Number of entries used in rho_levels (0 = none).
      logical :: mask_vanished_layers = .false.
         !! When `.true.`, remapped (non-layer) diagnostics fill target cells
         !! that overlap no water (below-bottom / pinched-out in a shallow
         !! column) with a missing sentinel and tag the NetCDF variable with
         !! `_FillValue` / `missing_value`.  Default `.false.` => those cells
         !! read 0 (bit-identical to the legacy writer).
      logical :: reproducing_sums = .false.
         !! PR-32: when `.true.`, the ocean console-conservation totals
         !! (Mass/KE/Salt/Heat + sea-ice area) use order-invariant
         !! extended-fixed-point (EFP) summation (`rdb_efp` +
         !! `halo_allreduce_efp_list`) instead of plain FP
         !! `!$acc parallel loop reduction(+:acc)` + `MPI_SUM` — the
         !! printed totals become bit-identical across rank counts and
         !! reduction orders, and the `Error` residual is formed via
         !! `efp_real_diff` (a fixed-point difference) rather than a
         !! double subtraction of two already-quantised totals.  Default
         !! `.false.` => the console block is byte-identical to the
         !! pre-PR-32 output.  See
         !! `docs/CAPABILITIES_AND_LIMITATIONS.md`'s conservation-contract
         !! section for the achievable guarantee + the `EFP_MAX_RANKS =
         !! 131072` envelope.
      character(len=16) :: output_precision = "double"
         !! Element width of the DIAGNOSTIC NetCDF data variables:
         !! "double" (default => byte-identical to the pre-knob writer) or
         !! "single" (fp32 => ~half the bytes per frame; diagnostic output
         !! is write-bandwidth bound, so this is the main lever on its
         !! cost).  fp32 carries ~7 decimal digits — ~1e-5 degC, ~1e-5 PSU,
         !! ~1e-9 m/s at 1 m/s — orders of magnitude below the model's own
         !! discretisation error.
         !!
         !! Scope, deliberately narrow: this knob reaches the diagnostic
         !! stream ONLY.  Restart files, console conservation totals
         !! and checksums
         !! are unconditionally working precision — a restart that does not
         !! round-trip exactly makes a resumed run a different run, so
         !! there is no way to ask for a lossy one.  The diag TIME
         !! coordinate variables also stay double.
         !!
         !! Caveat: a regression baseline that compares diagnostic NetCDF
         !! byte-for-byte will not match a "single" file; opting in means
         !! regenerating those baselines deliberately.
      character(len=16) :: diag_remap_scheme = "ppm"
         !! In-cell reconstruction for the conservative diagnostic vertical
         !! remap (z_fixed / density vgrids): "pcm", "plm", "ppm" (default),
         !! "ppm_h4", or "pqm".  All conserve the column integral; PPM is the
         !! accurate default (PCM stair-steps).  Inert when every diagnostic
         !! is on the native "layer" grid.
      character(len=512) :: diags = ""
         !! Unified diagnostic selection list — one authoritative spec that
         !! modifies the canonical default set.  Whitespace/comma-separated
         !! entries, each `name[:attr]...` with self-identifying colon
         !! attributes (order-free):
         !!
         !!   * `off`               — turn this diagnostic off (skip it)
         !!   * a cadence `1h`/`6h`/`30m`/`1d` — output interval override
         !!   * an op `instant`/`mean`/`max`/`min` — time-reduction override
         !!
         !! A `name` matching a canonical default applies its attributes (or
         !! drops it with `:off`); a `name` from the derived-diagnostic
         !! catalog (`rdb_ocean_diag_derived`) is added.  Example:
         !! `"vorticity_z:1d  KE:off  temperature:6h:mean"`.  Default empty => the
         !! canonical set unchanged => bit-identical output for existing
         !! namelists.  An unknown name fails loud at setup.  (Per-diagnostic
         !! vertical-coordinate attributes — `z*`/`sigma` — land with the
         !! multi-coordinate remap; today every diag uses the global `vgrid`.)
   end type ocean_diag_config_t

   integer, parameter, public :: OBC_MAX_TIDAL_CFG = 8
      !! Maximum tidal constituents per edge in the &ocean_bc_nml namelist.
      !! Mirrors OBC_MAX_TIDAL_CONSTITUENTS in rdb_ocean_boundary_types.

   type :: ocean_bc_config_t
      !! Boundary-condition configuration read from `&ocean_bc_nml`.
      !! All defaults reproduce a closed-wall run — existing nmls that
      !! omit this block are bit-identical to prior behaviour.
      character(len=16) :: west = "wall"
         !! BC type string for the west  edge.  Recognised: "wall" (default),
         !! "open", "tidal", "clamped", "sponge", "chapman", "periodic".
         !! "tripolar_fold" is accepted on the NORTH edge only (Murray
         !! bipolar cap; requires grid_config="tripolar" + periodic w/e).
      character(len=16) :: east = "wall"
      character(len=16) :: south = "wall"
      character(len=16) :: north = "wall"

      ! Clamped Dirichlet values — per-edge, all default 0.
      real(wp) :: west_clamped_eta = 0.0_wp
      real(wp) :: east_clamped_eta = 0.0_wp
      real(wp) :: south_clamped_eta = 0.0_wp
      real(wp) :: north_clamped_eta = 0.0_wp
      real(wp) :: west_clamped_u = 0.0_wp
      real(wp) :: east_clamped_u = 0.0_wp
      real(wp) :: south_clamped_v = 0.0_wp
      real(wp) :: north_clamped_v = 0.0_wp

      ! Inflow tracer concentrations (per-edge; used as clamped_tracer
      ! for salinity index 1 and temperature index 2 by default).
      real(wp) :: west_inflow_S = 35.0_wp
      real(wp) :: west_inflow_T = 10.0_wp
      real(wp) :: east_inflow_S = 35.0_wp
      real(wp) :: east_inflow_T = 10.0_wp
      real(wp) :: south_inflow_S = 35.0_wp
      real(wp) :: south_inflow_T = 10.0_wp
      real(wp) :: north_inflow_S = 35.0_wp
      real(wp) :: north_inflow_T = 10.0_wp

      ! Sponge band — single set shared by whichever edge is OBC_SPONGE.
      integer  :: sponge_width = 0
         !! Number of cells in the sponge band.  0 disables the sponge.
      real(wp) :: sponge_strength = 0.0_wp
         !! Peak relaxation rate (1/s) at the outer face of the band.
      logical  :: sponge_relax_tracers = .false.
         !! When .true., extend the sponge relaxation to h_layer + tracers.
         !! Default .false. — preserves bit-identity of existing sponge runs.

      ! Open-edge tracer reservoir (§1, v2).
      real(wp) :: res_lscale_out = 0.0_wp
         !! Outflow reservoir length scale (m).  0 = feature disabled (v1
         !! sign-switch path unchanged).  Non-zero engages the implicit
         !! Marchesiello et al. 2001 reservoir on all open-ish edges.
      real(wp) :: res_lscale_in = 0.0_wp
         !! Inflow reservoir length scale (m).  0 = instantaneous inflow
         !! (T_data applied on inflow unconditionally — same as v1 for inflow
         !! direction when feature is active).

      ! Per-layer Orlanski radiation (§2, v2).
      character(len=16) :: radiation_scheme = "anomaly"
         !! "anomaly" (default) = v1 BT-mean + zero-gradient anomaly.
         !! "orlanski" = per-layer implicit-upwind radiation (Orlanski 1976).
      real(wp) :: orlanski_rx_max = 10.0_wp
         !! Clamp on the Orlanski nondimensional phase speed (grid cells / step).
      real(wp) :: orlanski_gamma = 1.0_wp
         !! Running-mean weight (0=full running mean, 1=instant; default=1).
      real(wp) :: nudge_tau_in = 0.0_wp
         !! Inflow nudging timescale (s, Marchesiello et al. 2001).  0 = off.
      real(wp) :: nudge_tau_out = 0.0_wp
         !! Outflow nudging timescale (s).  0 = off.

      ! Full Flather with exterior velocity (§4, v2).
      character(len=16) :: flather_form = "legacy"
         !! "legacy" (default) = v1 Flather (u_ext=0, no interior vel);
         !! "full"   = half-characteristic form (Flather 1976).
      real(wp) :: west_ext_u = 0.0_wp  !! Exterior barotropic u, west  (m/s).
      real(wp) :: east_ext_u = 0.0_wp  !! Exterior barotropic u, east  (m/s).
      real(wp) :: south_ext_v = 0.0_wp  !! Exterior barotropic v, south (m/s).
      real(wp) :: north_ext_v = 0.0_wp  !! Exterior barotropic v, north (m/s).

      ! Tidal constituents (west edge only in v1 — extend as needed).
      integer  :: west_n_tidal = 0
      real(wp) :: west_tidal_amp(OBC_MAX_TIDAL_CFG) = 0.0_wp
      real(wp) :: west_tidal_phase(OBC_MAX_TIDAL_CFG) = 0.0_wp
      real(wp) :: west_tidal_omega(OBC_MAX_TIDAL_CFG) = 0.0_wp
      integer  :: east_n_tidal = 0
      real(wp) :: east_tidal_amp(OBC_MAX_TIDAL_CFG) = 0.0_wp
      real(wp) :: east_tidal_phase(OBC_MAX_TIDAL_CFG) = 0.0_wp
      real(wp) :: east_tidal_omega(OBC_MAX_TIDAL_CFG) = 0.0_wp
      integer  :: south_n_tidal = 0
      real(wp) :: south_tidal_amp(OBC_MAX_TIDAL_CFG) = 0.0_wp
      real(wp) :: south_tidal_phase(OBC_MAX_TIDAL_CFG) = 0.0_wp
      real(wp) :: south_tidal_omega(OBC_MAX_TIDAL_CFG) = 0.0_wp
      integer  :: north_n_tidal = 0
      real(wp) :: north_tidal_amp(OBC_MAX_TIDAL_CFG) = 0.0_wp
      real(wp) :: north_tidal_phase(OBC_MAX_TIDAL_CFG) = 0.0_wp
      real(wp) :: north_tidal_omega(OBC_MAX_TIDAL_CFG) = 0.0_wp

      ! Tidal-OBC nodal/astronomical correction (capability C3).
      logical  :: obc_tidal_nodal = .false.
         !! Apply the 18.6-yr nodal factor `f_c` + equilibrium/nodal phase
         !! `(V_c + u_c)` to the open-boundary tidal elevation forcing, keeping
         !! the interior body tide and the boundary tide phase-consistent (same
         !! generator, same shared `&ocean_tides_nml` reference epoch).  Default
         !! .false. ⇒ legacy static-phase OBC sum (bit-identical).  When .true.
         !! the `*_tidal_phase` becomes a Greenwich phase LAG (subtracted).

      ! Solid-wall velocity masking (MOM6-faithful wall-ghost land fill).
      logical  :: mask_wall_velocity = .true.
         !! Zero the T-cell wet-mask in the ghost cells beyond every solid
         !! WALL edge at setup, so the derived C-grid face masks
         !! (`wet_u`/`wet_v`/`wet_q`) are 0 at the wall face and the existing
         !! per-stage `mask_layer_velocities` clears the wall-normal velocity
         !! (MOM6's mask-in-the-update; the halo beyond a wall is land).
         !! DEFAULT ON (2026-07): a solid wall must carry zero normal velocity;
         !! leaving it garbage (flux-masked but nonzero) was a bug. Set
         !! `.false.` to recover the pre-fix legacy answer for a closed-basin
         !! run (e.g. to reproduce an old double_gyre baseline exactly).
         !! Without it a flat all-wet channel leaves the wall-normal face
         !! velocity as unmasked garbage — the wall flux is masked so mass
         !! conserves, but the raw v/u at the wall drifts to the maxvel clamp
         !! and paints a spurious vorticity band in diagnostics.  Only WALL
         !! edges are touched: periodic edges keep their wrapped (wet) ghosts,
         !! open/OBC edges keep the interior value (OBC override).  Default
         !! .false. ⇒ wall ghosts untouched ⇒ bit-identical to legacy runs
         !! (this also changes the continuity PPM mirror-h at the wall, which
         !! is why it must be opt-in).
   end type ocean_bc_config_t

   type :: ocean_config_t
      !! All ocean-path knobs, nested by concern.  Parallel to the
      !! state composition pattern (`ocean_state%dyn`, `%vcoord`, etc.).
      !! Read-side: see the per-sub-nml blocks in `read_config`.
      !! Call-side: `cfg%ocean%bt%use_cont_type`,
      !! `cfg%ocean%pgf%form`, etc.
      type(ocean_grid_config_t)       :: grid
      type(ocean_coriolis_config_t)   :: coriolis
      type(ocean_thermo_config_t)     :: thermo
      type(ocean_forcing_config_t)    :: forcing
         !! Surface-flux component-set gate (`&ocean_forcing_nml`, PR-12).
         !! Default off ⇒ byte-identical.
      type(ocean_ice_config_t)        :: ice
         !! Sea-ice slot config (`&ocean_ice_nml`).  Default off ⇒
         !! byte-identical (PR 0 scaffold).
      type(ocean_ice_ic_config_t)     :: ice_ic
         !! Sea-ice initial-condition config (`&ocean_ice_ic_nml`, PR 24).
         !! Default `conc_config="zero"` ⇒ byte-identical.
      type(ocean_restore_config_t)    :: restore
      type(ocean_geothermal_config_t) :: geothermal
      type(ocean_sponge_config_t)     :: sponge
      type(ocean_tracers_config_t)    :: tracers
      type(ocean_bt_config_t)         :: bt
      type(ocean_debug_config_t)      :: debug
         !! Forensic probes (`&ocean_debug_nml`).  Default all-off ⇒
         !! bit-identical.
      type(ocean_pgf_config_t)        :: pgf
      type(ocean_eos_config_t)        :: eos
      type(ocean_bdrag_config_t)      :: bdrag
      type(ocean_tdrag_config_t)      :: tdrag
         !! Ice-shelf TOP drag (`&ocean_tdrag_nml`).  Default off ⇒
         !! bit-identical.
      type(ocean_hdiff_config_t)      :: hdiff
         !! Along-coordinate tracer Laplacian (`&ocean_hdiff_nml`).
         !! Default `kappa_h = 0.0` ⇒ bit-identical.
      type(ocean_hvisc_config_t)      :: hvisc
      type(ocean_vmix_config_t)       :: vmix
      type(ocean_vdiff_config_t)      :: vdiff
      type(ocean_epbl_config_t)       :: epbl
      type(ocean_wave_speed_config_t) :: wavespeed
      type(ocean_foxkemper_config_t)  :: foxkemper
      type(ocean_kappa_shear_config_t) :: kshear
      type(ocean_slopes_config_t)     :: slopes
      type(ocean_gm_config_t)         :: gm
      type(ocean_redi_config_t)       :: redi
      type(ocean_varmix_config_t)     :: varmix
      type(ocean_meke_config_t)       :: meke
      type(ocean_tidal_mixing_config_t) :: tidal_mixing
      type(ocean_conv_config_t)       :: conv
      type(ocean_porous_config_t)     :: porous
      type(ocean_ddiff_config_t)      :: ddiff
      type(ocean_tides_config_t)      :: tides
      type(ocean_psurf_config_t)      :: psurf
      type(ocean_cavity_dyn_config_t) :: cavity_dyn
         !! Static ice-shelf cavity geometry (`&ocean_cavity_dyn_nml`).
         !! Default-OFF ⇒ bit-identical.
      type(ocean_cavity_melt_config_t) :: cavity_melt
         !! Ice-shelf basal-melt thermodynamics
         !! (`&ocean_cavity_melt_nml`).  Default-OFF ⇒ bit-identical.
      type(ocean_continuity_config_t)  :: continuity
      type(ocean_isopycnal_config_t)   :: isopycnal
         !! Lagrangian grounding-stability controls (`&ocean_isopycnal_nml`).
         !! Default-OFF ⇒ bit-identical.
      type(ocean_topo_config_t)        :: topo
      type(ocean_ic_config_t)         :: ic
      type(ocean_zinit_config_t)      :: zinit
      type(ocean_data_config_t)       :: data
         !! Shared time-varying NetCDF input reader (`&ocean_data_nml`,
         !! PR-14).  Two knobs only — see `ocean_data_config_t`.
      type(ocean_dataovr_config_t)    :: dataovr
         !! File-backed surface forcing (`&ocean_dataovr_nml`, PR-15) —
         !! the first production consumer of `data`.  Default-OFF ⇒
         !! bit-identical.
      type(ocean_diag_config_t)       :: diag
      type(ocean_bc_config_t)         :: bc
         !! Open-boundary condition config (`&ocean_bc_nml`).
         !! All defaults = wall — existing nmls without this block are
         !! bit-identical to prior behaviour.
      type(ocean_wetdry_config_t)     :: wetdry
         !! Dynamic wetting/drying (`&ocean_wetdry_nml`).  Default off ⇒
         !! byte-identical.
      type(ocean_mpi_config_t)        :: mpi
         !! Multi-rank MPI debug / tuning controls (`&ocean_mpi_nml`).
         !! Default off ⇒ bit-identical.
   end type ocean_config_t

   type :: config_t
      !! Container for all runtime simulation parameters

      ! Simulation regime.  "ocean" (Arakawa C-grid + continuity-PPM) is
      ! the only regime this build ships — the coastal A-grid path was
      ! carved out into its own repository.  Retained as a knob so an
      ! existing namelist keeps parsing; validate_config rejects any other
      ! value.  See docs/ROADMAP_OCEAN.md.
      character(len=16) :: sim_type = "ocean"
         !! Simulation regime: "coastal" or "ocean"

      ! Grid parameters
      integer :: nx = 200
         !! Number of physical cells in x-direction
      integer :: ny = 1
         !! Number of physical cells in y-direction
      real(wp) :: dx = 1.0_wp
         !! Cell size in x-direction (m)
      real(wp) :: dy = 1.0_wp
         !! Cell size in y-direction (m)
      integer :: nghost = 3
         !! Number of ghost cells on each side.  Default 3 (not 2): the
         !! ocean dyn-core's continuity/tracer PPM reconstruction is a
         !! 5-point stencil, so a rank-seam face's ghost-side donor needs
         !! two neighbours beyond itself — full-order reconstruction AT the
         !! seam requires three ghost columns.  At nghost=2 the PPM
         !! local-array-edge fallback fires at the seam and degrades the
         !! face to first order, producing a seam-inconsistent tracer mass
         !! flux: it conserves a uniform tracer (constant field cancels)
         !! but LEAKS a structured one (advecting T with rank-A-out /=
         !! rank-B-in on the shared seam face drifts hT), invisible to the
         !! mass/salt closure gates.  Measured on seamount_bench_full
         !! (np2 x-split): heat closure grew to +5e-11 by day 2 at nghost=2,
         !! stays at round-off (~5e-14, no growth) at nghost=3.  Single-rank
         !! runs are unaffected by the wider halo (same interior); the extra
         !! ghost column only costs a little memory.

      ! Time parameters
      real(wp) :: t_end = 1.0_wp
         !! Simulation end time, interpreted in `time_unit` (default s).
      real(wp) :: cfl = 0.45_wp
         !! CFL number for adaptive timestep
      real(wp) :: dt_max = 1.0e10_wp
         !! Maximum allowable timestep (s) — always in seconds, NOT
         !! affected by `time_unit`.  Timesteps are naturally short
         !! and CFL-bounded; writing them in `"day"` would force
         !! awkward fractional values.
      real(wp) :: dt_fixed = 0.0_wp
         !! Fixed timestep (s), 0 = adaptive CFL.  Same convention
         !! as `dt_max` — always seconds.
      integer :: cfl_interval = 1
         !! Recompute CFL timestep every N steps (1 = every step)
      character(len=8) :: time_unit = "s"
         !! Unit applied to the long-time fields: `t_end`,
         !! `status_interval` (logging cadence), `ocean_diag_dt_out`
         !! (diag-write cadence).  `dt_fixed` and `dt_max` stay in
         !! seconds — see their docstrings.  Read at namelist parse,
         !! then multiplied through so the rest of the code keeps
         !! using seconds internally.  Valid values: `"s"`, `"min"`,
         !! `"hr"`, `"day"`, `"year"`.  Year is 365.25 days (Julian,
         !! matches MOM6 / UDUNITS convention).  Default `"s"`
         !! preserves backward compatibility — every existing nml is
         !! bit-identical.

      ! Physics parameters
      real(wp) :: manning_n = 0.0_wp
         !! Manning roughness coefficient
      real(wp) :: wind_stress_x = 0.0_wp
         !! Surface wind stress in x-direction (Pa)
      real(wp) :: wind_stress_y = 0.0_wp
         !! Surface wind stress in y-direction (Pa)
      real(wp) :: coriolis_f = 0.0_wp
         !! Coriolis parameter f (1/s), typically 2*Omega*sin(lat)

      ! Output parameters
      logical :: output_to_file = .false.
         !! Enable file output (NetCDF snapshots)
      character(len=256) :: output_dir = "./output"
         !! Directory for output files
      real(wp) :: restart_interval = 0.0_wp
         !! Time between restart file writes (s), 0 = no restarts
      logical :: compress_output = .false.
         !! Enable deflate compression for NetCDF output (NetCDF4/HDF5)
      integer :: compress_level = 1
         !! Deflate compression level (1=fast, 9=max). 1 is usually optimal.
      logical :: use_io_server = .false.
         !! Dedicate one MPI rank per node as I/O server (MPI only)

      ! Output variable selection (multilayer)

      ! I/O file parameters
      character(len=256) :: bathymetry_file = ""
         !! Path to NetCDF bathymetry file (empty = flat bottom)
      character(len=256) :: restart_file = ""
         !! Path to restart file for warm start (empty = cold start)

      ! Nesting parameters

      ! Boundary parameters
      character(len=16) :: bc_west = "wall"
         !! West boundary: "wall", "open", "tidal", "nested"
      character(len=16) :: bc_east = "wall"
         !! East boundary: "wall", "open", "tidal", "nested"
      character(len=16) :: bc_south = "wall"
         !! South boundary: "wall", "open", "tidal", "nested"
      character(len=16) :: bc_north = "wall"
         !! North boundary: "wall", "open", "tidal", "nested"

      ! Tidal forcing parameters

      ! Multi-constituent tidal forcing
      integer :: n_tidal_constituents = 0
         !! Number of active tidal constituents (0 = use legacy single)
      real(wp) :: tidal_amp(MAX_TIDAL_CONSTITUENTS) = 0.0_wp
         !! Constituent amplitudes (m)
      real(wp) :: tidal_phase(MAX_TIDAL_CONSTITUENTS) = 0.0_wp
         !! Constituent phases (radians)
      real(wp) :: tidal_omega(MAX_TIDAL_CONSTITUENTS) = 0.0_wp
         !! Constituent angular frequencies (rad/s)

      ! Inflow/discharge/clamped BC parameters
      real(wp) :: inflow_salinity = -1.0_wp
         !! Prescribed salinity for inflow BC (PSU); <0 = zero-gradient
      real(wp) :: inflow_temperature = -999.0_wp
         !! Prescribed temperature for inflow BC (degC); <0 = zero-gradient

      ! Sponge layer parameters
      integer :: sponge_width = 0
         !! Sponge layer width in cells
      real(wp) :: sponge_strength = 0.0_wp
         !! Sponge relaxation rate (1/s)

      ! Metadata parameters

      ! MPI domain decomposition
      integer :: px = 1
         !! Number of MPI processes in x-direction
      integer :: py = 1
         !! Number of MPI processes in y-direction

      ! Non-hydrostatic pressure correction
      integer :: nz_layers = 2
         !! Number of sigma layers for NH solver

      ! Coupled hydrostatic multi-layer
      logical :: use_multilayer = .false.
         !! Enable coupled hydrostatic vertical layers (alternative to NH)
      real(wp) :: rho_0 = 1000.0_wp
         !! Reference density for EOS (kg/m^3)
      real(wp) :: kpp_ri_crit = 0.3_wp
         !! Critical bulk Richardson number for the KPP BL-depth sweep.
      real(wp) :: kpp_cs_nonlocal = 6.3_wp
         !! Non-local (counter-gradient) transport coefficient C_s
         !! (LMD94 eq 20, limit value).
      real(wp) :: kpp_c_vt2 = 1.8_wp
         !! Unresolved-turbulence coefficient for the V_t^2 term in the
         !! bulk-Ri denominator (LMD94 eq 23).  Set 0 to disable V_t^2.
      real(wp) :: hdiff_kappa = 0.0_wp
         !! Horizontal diffusion coefficient for tracers (m^2/s).
         !! VESTIGIAL: its only consumers were the coastal solvers, which
         !! left with the carve-out, so nothing reads it today.  The
         !! along-coordinate equivalent is the separate
         !! `&ocean_hdiff_nml kappa_h` knob (`rdb_ocean_hdiff_tracer`).
         !! Smooths sharp density fronts to reduce BPG overshoot.
         !! Typical values: 1-10 m^2/s depending on grid resolution.

      ! Tracer parameters (salinity + temperature)
      real(wp) :: initial_salinity = 35.0_wp
         !! Initial salinity for all layers (PSU).  Used as a uniform IC
         !! unless `S_init_surface` and `S_init_bottom` are BOTH set
         !! non-zero, in which case the seed builds a linear S(z) profile
         !! from `S_init_bottom` at k=1 (bed) to `S_init_surface` at
         !! k=nz (surface).  NOTE the stable polarity is the inverse of
         !! temperature: salty/dense water belongs at the bed, so a stable
         !! haline column has `S_init_bottom > S_init_surface` (fresh
         !! surface, e.g. a river-plume column).
      real(wp) :: S_init_surface = 0.0_wp
         !! Initial salinity at the surface layer (k=nz_ml, PSU).
         !! Activates linear S(z) stratification when set non-zero
         !! together with `S_init_bottom`.
      real(wp) :: S_init_bottom = 0.0_wp
         !! Initial salinity at the bed (k=1, PSU).
      real(wp) :: S_ref = LEGACY_TRACER_S_REF
         !! RETIRED coastal-legacy EOS reference salinity (PSU).  The live
         !! ocean-path spelling is `&ocean_ic_nml S_ref`; moving this one
         !! off its default is a fail-loud configure error.
      real(wp) :: beta_S = LEGACY_TRACER_BETA_S
         !! RETIRED coastal-legacy haline contraction coefficient
         !! (kg/m^3 per PSU).  Live spelling: `&ocean_ic_nml beta_S`.
      real(wp) :: S_min = 0.0_wp
         !! Lower physical bound for salinity (PSU)
      real(wp) :: S_max = 40.0_wp
         !! Upper physical bound for salinity (PSU)
      real(wp) :: kappa_S_bg = 1.0e-5_wp
         !! Background vertical salinity diffusivity (m^2/s)
      real(wp) :: initial_temperature = 15.0_wp
         !! Initial potential temperature for all layers (degC).  Used
         !! as a uniform IC unless `T_init_surface` and `T_init_bottom`
         !! are BOTH set non-zero, in which case the seed builds a
         !! linear T(z) profile from `T_init_bottom` at k=1 (bed) to
         !! `T_init_surface` at k=nz (surface).  Stratified IC is the
         !! prerequisite for baroclinic-instability-driven eddies.
      real(wp) :: T_init_surface = 0.0_wp
         !! Initial temperature at the surface layer (k=nz_ml, degC).
         !! Activates linear T(z) stratification when set non-zero
         !! together with `T_init_bottom`.
      real(wp) :: T_init_bottom = 0.0_wp
         !! Initial temperature at the bed (k=1, degC).
      real(wp) :: T_ref = LEGACY_TRACER_T_REF
         !! RETIRED coastal-legacy EOS reference temperature (degC).  The
         !! live ocean-path spelling is `&ocean_ic_nml T_ref`; moving this
         !! one off its default is a fail-loud configure error.
      real(wp) :: alpha_T = LEGACY_TRACER_ALPHA_T
         !! RETIRED coastal-legacy thermal expansion coefficient
         !! (kg/m^3 per degC).  Live spelling: `&ocean_ic_nml alpha_T`.
      real(wp) :: T_min = -2.0_wp
         !! Lower physical bound for temperature (degC); seawater freezing
      real(wp) :: T_max = 40.0_wp
         !! Upper physical bound for temperature (degC)
      real(wp) :: kappa_T_bg = 1.0e-5_wp
         !! Background vertical temperature diffusivity (m^2/s)
      ! Sediment transport (single class, gated; default off = bit-identical).

      ! Initial condition selection
      real(wp) :: h0 = 1.0_wp
         !! Background depth for gaussian_hump IC (m)

      ! Mode splitting parameters

      ! Vertical coordinate parameters
      character(len=16) :: vcoord_type = "sigma"
         !! Vertical coordinate type: "sigma", "zsigma", "zstar", etc.
      character(len=16) :: thickness_config = "sigma"
         !! INITIAL layer-thickness profile for the ocean path (`sim_type =
         !! 'ocean'`); ignored by the coastal regimes.  Distinct from
         !! `vcoord_type`, which selects the RUNNING coordinate: this knob only
         !! decides what `h_layer` is seeded to at t=0.
         !!
         !!   "sigma"     (default) — even split of the LOCAL depth,
         !!                `h_layer(i,j,k) = b(i,j)/nz_ml`.  Bit-identical to
         !!                the pre-knob behaviour for every existing nml.
         !!   "uniform_z" — MOM6 `THICKNESS_CONFIG="uniform"` port: uniform z
         !!                interfaces laid over the GLOBAL `ocean_max_depth`,
         !!                clipped bottom-up to the local bathymetry, with
         !!                sub-floor layers collapsed to the isopycnal
         !!                `angstrom_h`.  Gives FLAT resting isopycnals under a
         !!                horizontally-uniform density stack — the layered /
         !!                isopycnal (VCOORD_LAGRANGIAN) resting state.  With
         !!                "sigma" the density interfaces instead follow the
         !!                bathymetry, which puts the full `rho_range` contrast
         !!                across every shelf break at t=0 and releases that
         !!                available potential energy as a slumping gravity
         !!                current.
      character(len=16) :: remap_method = "ppm"
         !! Vertical remapping method: "pcm", "plm", "ppm", "ppm_h4", "pqm".  Default PPM
         !! (piecewise parabolic, Colella & Woodward 1984) — the state-of-the-art
         !! conservative vertical-remap stencil for the coastal ALE coords. "pqm" is
         !! piecewise quartic (White & Adcroft 2008); it falls back to PPM for nz < 5.
      ! Full MOM6 z* parameters (VCOORD_ZSTAR_FULL)
      real(wp) :: zstar_h_surf_target = 0.0_wp
         !! Target physical thickness of the surface layer (m).
         !! 0 => auto: fall back to the lite (uniform sigma) profile
         !! per column (effectively zstar-lite per-column).  Set > 0 to
         !! anchor the surface layer at a fixed thickness regardless of H.
      character(len=16) :: zstar_stretching = "log"
         !! Surface-concentration stretching: "log" | "uniform"
      real(wp) :: zstar_h_min = 1.0e-4_wp
         !! Vanishing-layer floor (m).  Layers that would land below the
         !! local bed get clipped to this thickness rather than going to zero.
      integer :: zstar_n_surf = 0
         !! Number of "fine" near-surface layers using stretching.
         !! 0 => auto (use max(1, nz/3))
      ! Isopycnal coordinate (VCOORD_RHO) parameters
      real(wp) :: rho_ref_pressure = 2.0e7_wp
         !! Reference pressure (Pa, default 2e7 = 2000 dbar) for the
         !! potential density that defines the VCOORD_RHO coordinate.
      real(wp) :: rho_target_light = 1020.0_wp
         !! Lightest target interface potential density (kg/m³) — the
         !! surface interface for VCOORD_RHO.  With `rho_target_dense`
         !! builds the uniform light→dense target-density linspace
         !! (MOM6 `ALE_COORDINATE_CONFIG=UNIFORM` analogue for the RHO
         !! mode).  A stretched (RFNC1-style) profile is a P3 follow-up.
      real(wp) :: rho_target_dense = 1030.0_wp
         !! Densest target interface potential density (kg/m³) — the bed
         !! interface for VCOORD_RHO.  See `rho_target_light`.
      ! ALE-regrid refinements (ocean path)
      real(wp) :: regrid_time_scale = 0.0_wp
         !! Grid time-filter timescale τ (s).  The ALE regrid relaxes the
         !! coordinate a fraction dt/(τ+dt) toward the target each thermo
         !! step instead of jumping (White & Adcroft 2008), damping the
         !! σ/z* per-step grid-motion shock.  0 (default) = jump to target
         !! = bit-identical.
      logical :: remap_vel_conserve_ke = .false.
         !! KE-conserving rescale of the remapped layer velocities: scale
         !! the baroclinic anomaly per column so column KE is preserved
         !! (Adcroft & Hallberg 2006), capped 1.25×; barotropic mean
         !! untouched.  .false. (default) = momentum-only = bit-identical.
      logical :: remap_boundary_extrap = .false.
         !! Linear-exact one-sided reconstruction in the ALE remap's two
         !! boundary cells (MOM6 `BOUNDARY_EXTRAPOLATION`) instead of the
         !! PCM flatten, which leaves PLM/PPM/PPM_H4/PQM first-order at
         !! `k=1` and `k=nz`.  .false. (default) = bit-identical.
      logical :: remap_nonuniform_weights = .false.
         !! Non-uniform-grid reconstruction weights in the ALE remap's PLM
         !! slope and PPM edge estimate (Colella & Woodward 1984 eqs
         !! 1.6-1.8) instead of their equal-thickness specialisations, which
         !! are linear-exact only on a UNIFORM source column.  PPM_H4 and
         !! PQM already carry thickness-weighted stencils and are unchanged.
         !! .false. (default) = bit-identical.
      logical :: remap_check_preconditions = .false.
         !! Assert, once per ALE remap (thermo cadence), that every column
         !! satisfies what the overlap sweep has always assumed: non-negative
         !! source AND target thicknesses, and matching column totals.  Both
         !! are caller obligations and neither was ever checked; a violation
         !! silently CREATES or DELETES tracer mass.  Diagnostic knob —
         !! .false. (default) = the check never runs = bit-identical.
      logical :: check_vanished_content = .false.
         !! **I1′ tripwire** — assert `h_layer <= H_VANISHED ⇒ hTr =
         !! h_layer·c_live` (the donor live layer's concentration; `hTr = 0`
         !! in a column with no live layer) for every registered tracer, once
         !! per outer step, immediately after the enforcement point that
         !! establishes it (`multilayer_state_t%enforce_vanished_content`).
         !! A violation logs the offending cell count and the worst
         !! `|hTr − h·c_live|` and `error stop`s.
         !!
         !! HEAVY only in the sense that it adds two device reductions per
         !! tracer per step; the healthy path does no H←D copy.  Default
         !! `.false.`.  Turn it ON for any configuration whose coordinate
         !! actually vanishes layers (`z_fixed`, `zstar_full`, wet/dry) —
         !! the stability suite does.  Inert on a family with no fillers.
         !!
         !! See `src/core/ocean/README.md`, "The vanished-layer content
         !! rule".
      logical :: zfixed_closed_faces = .false.
         !! **Partial-step z-level face closure** under
         !! `vcoord_type = "z_fixed"` (Adcroft, Hill & Marshall 1997;
         !! Losch 2008 for the ice-shelf cavity).  A layer whose nominal
         !! geopotential range lies inside the bed — or inside the ice
         !! draft — carries an inert FILLER of thickness `zstar_h_min`
         !! (`<= H_VANISHED`).  A velocity face where layer `k` is a
         !! filler on EITHER side is not a thin passage, it is a WALL for
         !! that layer: no normal velocity, no mass / tracer flux, and
         !! FREE-SLIP on the tangential component.  Leaving it open makes
         !! the FV pressure gradient integrate across a staircase step of
         !! height `Δz_step`, which drives
         !! `|ρ′|·g·Δz_step/(ρ₀·dx)` out of a resting stratified state —
         !! independent of the filler thickness, so no `h`-gate reaches
         !! it.
         !!
         !! ON builds a STATIC 0/1 per-layer face mask
         !! (`ocean_metrics_t%open_u` / `open_v`) once at configure from
         !! the `z_fixed` target at `η = 0`, and composes it
         !! multiplicatively with the land metrics and the porous-barrier
         !! open-area fraction:
         !! `dy_eff(i,j,k) = dy_cu(i,j)·por_face_area_u(i,j,k)·open_u(i,j,k)`.
         !!
         !! Default `.false.` ⇒ the mask arrays stay at their `(1,1,1)`
         !! placeholder, no kernel branch is taken, byte-identical.
         !! Refused on any coordinate but `z_fixed`, and without a
         !! resolved `z_fixed_h_ref` (there would be no fillers to close).

      ! Logging parameters
      character(len=16) :: log_level = "info"
         !! Log verbosity: "debug", "verbose", "info", "performance", "warning", "error"
      real(wp) :: status_interval = 0.0_wp
         !! Print status every N seconds of simulation time (0 = default every 100 steps)

      ! Ocean idealised setup (sim_type='ocean' only).  Drives the
      ! double-gyre and other analytical configurations.  Defaults
      ! reproduce a flat-bottom basin with no analytical wind / Coriolis
      ! pattern so existing tests are unaffected.
      ! `topo_config`        migrated to `cfg%ocean%topo%topo_config`.
      ! `ic_config`        migrated to `cfg%ocean%ic%ic_config`.
      ! `ocean_alpha_T`    migrated to `cfg%ocean%ic%alpha_T`.
      ! `ocean_rho_0`      migrated to `cfg%ocean%ic%rho_0`.
      ! `ocean_layer_rho_init` migrated to `cfg%ocean%ic%layer_rho_init`.
      ! `eady_dT_dy`       migrated to `cfg%ocean%ic%eady_dT_dy`.
      ! `eady_dT_dz`       migrated to `cfg%ocean%ic%eady_dT_dz`.
      ! `eady_T_ref`       migrated to `cfg%ocean%ic%eady_T_ref`.
      ! `eady_pert_amp`    migrated to `cfg%ocean%ic%eady_pert_amp`.
      ! `eady_pert_seed`   migrated to `cfg%ocean%ic%eady_pert_seed`.
      ! `ga_eta_amp`       migrated to `cfg%ocean%ic%ga_eta_amp`.
      ! `ga_length_scale`  migrated to `cfg%ocean%ic%ga_length_scale`.
      ! `ga_x_center`      migrated to `cfg%ocean%ic%ga_x_center`.
      ! `ga_y_center`      migrated to `cfg%ocean%ic%ga_y_center`.
      ! `wind_config`     migrated to `cfg%ocean%topo%wind_config`.
      ! `ocean_max_depth` migrated to `cfg%ocean%topo%max_depth`.
      ! `ocean_use_bt_cont_type`     migrated to `cfg%ocean%bt%use_cont_type`.
      ! `ocean_bt_upstream_h_face`   migrated to `cfg%ocean%bt%upstream_h_face`.
      ! `ocean_bt_cont_corr_bounds`  migrated to `cfg%ocean%bt%cont_corr_bounds`.
      ! `ocean_edge_depth`  migrated to `cfg%ocean%topo%edge_depth`.
      ! `ocean_slope_scale` migrated to `cfg%ocean%topo%slope_scale`.
      ! `taux_magnitude`    migrated to `cfg%ocean%topo%taux_magnitude`.
      ! `coriolis_beta`     migrated to `cfg%ocean%topo%coriolis_beta`.
      ! `coriolis_y_ref`    migrated to `cfg%ocean%topo%coriolis_y_ref`.
      ! `ocean_coriolis_form` migrated to `cfg%ocean%coriolis%form`
      ! (PR-B1 first knob-group migration, proof of pattern).
      ! `ocean_bdrag_form`           migrated to `cfg%ocean%bdrag%form`.
      ! `ocean_bdrag_cd`             migrated to `cfg%ocean%bdrag%cd`.
      ! `ocean_bdrag_r`              migrated to `cfg%ocean%bdrag%r`.
      ! `ocean_bdrag_hbbl`           migrated to `cfg%ocean%bdrag%hbbl`.
      ! `ocean_bdrag_bg_vel`         migrated to `cfg%ocean%bdrag%bg_vel`.
      ! `ocean_bdrag_bbl_thick_min`  migrated to `cfg%ocean%bdrag%bbl_thick_min`.
      ! `ocean_bdrag_bed_factor`     migrated to `cfg%ocean%bdrag%bed_factor`.
      ! `ocean_continuity_h_min` migrated to `cfg%ocean%continuity%h_min`.
      ! `ocean_continuity_ppm_limit_pos` migrated to
      ! `cfg%ocean%continuity%ppm_limit_pos`.
      ! `ocean_debug_bt_budget`         migrated to `cfg%ocean%debug%budget`.
      ! `ocean_bt_correction_bc_pgf`    migrated to `cfg%ocean%bt%correction_bc_pgf`.
      ! `ocean_bt_substep_drag`         migrated to `cfg%ocean%bt%substep_drag`.
      ! `ocean_bt_correction_visc_rem`  migrated to `cfg%ocean%bt%correction_visc_rem`.
      ! `ocean_enable_thermodynamics` migrated to
      ! `cfg%ocean%thermo%enable_thermodynamics`.
      ! `ocean_enable_ideal_age` migrated to
      ! `cfg%ocean%tracers%enable_ideal_age`.
         !! Use for circulation-pathway diagnostics + spurious
         !! diapycnal-mixing validation against MOM6.

      ! ---- Horizontal viscosity ----
      ! `ocean_lateral_closure` migrated to `cfg%ocean%hvisc%lateral_closure`.
      ! `ocean_c_smag`          migrated to `cfg%ocean%hvisc%c_smag`.
      ! `ocean_c_leith`         migrated to `cfg%ocean%hvisc%c_leith`.
      ! `ocean_kh_vel_scale`    migrated to `cfg%ocean%hvisc%kh_vel_scale`.
      ! `ocean_ah_bg`           migrated to `cfg%ocean%hvisc%ah_bg`.
      ! `ocean_ah_max`          migrated to `cfg%ocean%hvisc%ah_max`.
      ! `ocean_smag_ah`         migrated to `cfg%ocean%hvisc%smag_ah`.
      ! `ocean_smag_bi_const`   migrated to `cfg%ocean%hvisc%smag_bi_const`.
      ! `ocean_nu_4_bg`         migrated to `cfg%ocean%hvisc%nu_4_bg`.
      ! `ocean_nu_4_max`        migrated to `cfg%ocean%hvisc%nu_4_max`.
      ! `ocean_direct_stress`  migrated to `cfg%ocean%vmix%direct_stress`.
      ! `ocean_hmix_stress`    migrated to `cfg%ocean%vmix%hmix_stress`.
      ! `ocean_kv_ml_invz2`    migrated to `cfg%ocean%vmix%kv_ml_invz2`.
      ! `ocean_hmix_fixed`     migrated to `cfg%ocean%vmix%hmix_fixed`.
      ! `ocean_harmonic_visc`  migrated to `cfg%ocean%vmix%harmonic_visc`.

      ! `ocean_dt_therm_ratio` migrated to `cfg%ocean%vmix%dt_therm_ratio`.

      ! `ocean_pgf_form`     migrated to `cfg%ocean%pgf%form`.
      ! `ocean_gprime_gfs`   migrated to `cfg%ocean%pgf%gprime_gfs`.
      ! `ocean_gprime_gint`  migrated to `cfg%ocean%pgf%gprime_gint`.
      ! `ocean_gfs_scale`    migrated to `cfg%ocean%pgf%gfs_scale`.
      ! `ocean_maxvel`       migrated to `cfg%ocean%pgf%maxvel`.
      ! `ocean_nu_h` migrated to `cfg%ocean%hvisc%nu_h`.
      ! `ocean_nu_4` migrated to `cfg%ocean%hvisc%nu_4`.
      ! `n_inner`         migrated to `cfg%ocean%bt%n_inner`.
      ! `auto_n_inner`    migrated to `cfg%ocean%bt%auto_n_inner`.
      ! `cfl_bt_safety`   migrated to `cfg%ocean%bt%cfl_bt_safety`.
      ! `ocean_bebt`      migrated to `cfg%ocean%bt%bebt`.
      ! `ocean_bt_correction_h_weighted` migrated to
      ! `cfg%ocean%bt%correction_h_weighted`.
      ! `vmix_use_closure` migrated to `cfg%ocean%vmix%use_closure`.
      ! `vmix_use_kpp`     migrated to `cfg%ocean%vmix%use_kpp`.

      ! Ocean diag manager knobs migrated to `cfg%ocean%diag` —
      ! `enabled`, `filename`, `dt_out`, `vgrid`, `z_levels`, `n_z_levels`.

      type(ocean_config_t) :: ocean
         !! Nested ocean knobs.  Per-concern sub-types live under here
         !! (coriolis, thermo, bt, pgf, bdrag, hvisc, vmix, continuity,
         !! topo, ic, diag).  PR-B1 migrates knobs from this struct's
         !! flat fields into the nested form one group at a time —
         !! callers go from `cfg%ocean_X` → `cfg%ocean%group%X`.
         !! PR-B2 will split the giant `&ocean_setup_nml` to match.
   end type config_t

contains

   subroutine read_config(filename, cfg, schema, ierr)
      !! Read simulation configuration from a namelist file.
      !!
      !! Groups still read natively (`read(nml=)`) are seeded from the
      !! `cfg` defaults, read, and copied back as before.  Groups that
      !! have migrated to the strict schema (`rdb_nml_schema`) are
      !! validated and applied here, AFTER the native reads, by building
      !! a local schema (defaults captured from `cfg`'s pristine field
      !! initialisers) and strict-parsing the file — so `error stop` on
      !! any unknown key / range / enum violation fires for every caller
      !! (driver, FFI, unit tests), not just the production entry points.
      !!
      !! When `schema` is present the built schema is returned for the
      !! caller's parameter-doc dumps.  Its key pointers alias into
      !! `cfg`; both entry points keep their `cfg` alive (and `target`),
      !! so returning it is safe.
      character(len=*), intent(in) :: filename
         !! Path to the namelist input file
      type(config_t), target, intent(out) :: cfg
         !! Populated configuration
      type(nml_schema_t), intent(out), optional :: schema
         !! Built + parsed schema (aliases into `cfg`); for doc dumps.
      integer, intent(out), optional :: ierr
         !! Non-zero on a strict-parse failure (unknown group/key,
         !! range/enum violation) when present; absent behaves as today
         !! (`error stop`).

      ! Dispatch so the schema is built DIRECTLY into the caller's
      ! variable when requested — never copied.  Intrinsic assignment of
      ! nml_schema_t deep-copies polymorphic allocatable key boxes, which
      ! NVHPC miscompiles (heap corruption in the FFI tests); building in
      ! place sidesteps it and is cheaper anyway.
      if (present(schema)) then
         call read_config_impl(filename, cfg, schema, ierr=ierr)
      else
         block
            type(nml_schema_t) :: sch
            call read_config_impl(filename, cfg, sch, ierr=ierr)
         end block
      end if

   end subroutine read_config

   subroutine read_config_impl(filename, cfg, sch, ierr)
      !! Body of [[read_config]]: build the strict schema (capturing
      !! defaults from the pristine `cfg` field initialisers), parse the
      !! file (validate + apply straight into `cfg`), run the post-parse
      !! time-unit cascade.
      character(len=*), intent(in) :: filename
         !! Path to the namelist input file
      type(config_t), target, intent(out) :: cfg
         !! Populated configuration
      type(nml_schema_t), intent(out) :: sch
         !! Built + parsed schema (aliases into `cfg`)
      integer, intent(out), optional :: ierr
         !! Non-zero on a strict-parse failure when present; absent
         !! behaves as today (`error stop`).

      ! `cfg` arrives at its pristine field-initialiser defaults.  Every
      ! namelist group is now on the strict schema: build it (capturing
      ! those defaults), then parse — which validates AND applies every
      ! override straight into `cfg` via the schema's key pointers.  No
      ! native readers, no flat shadow locals, no copy-back remain.

      ! Build the strict schema so its key constructors capture defaults
      ! from the pristine `cfg` field initialisers.
      call build_rdb_schema(cfg, sch)

      ! A missing namelist file is not an error — the pristine `cfg`
      ! defaults stand (production main.F90 fails loudly on a missing
      ! file earlier; tests + the FFI rely on this defaults-on-absent
      ! behaviour).  Guard the strict parse so it never error-stops on
      ! "cannot open file".
      block
         logical :: file_exists
         inquire (file=trim(filename), exist=file_exists)
         if (.not. file_exists) then
            call logger%debug("Namelist file not found, keeping defaults: "//trim(filename))
            if (present(ierr)) ierr = OCEAN_STATUS_OK
            return
         end if
      end block

      ! Strict schema parse: validate + APPLY every group straight into
      ! `cfg`.  When `ierr` is present, route through the schema's own
      ! non-aborting `status`/`errors` pair (already implemented in
      ! `rdb_nml_schema`) instead of the bare (strict, error-stopping)
      ! form — every collected message is still logged, just without the
      ! `error stop`.  Absent `ierr` ⇒ unchanged bare call, for EVERY
      ! existing caller (driver, benchmarks, unit tests).
      if (present(ierr)) then
         block
            integer :: parse_status, e
            character(len=:), allocatable :: parse_errors(:)
            call sch%parse(trim(filename), status=parse_status, errors=parse_errors)
            if (parse_status /= 0) then
               if (allocated(parse_errors)) then
                  do e = 1, size(parse_errors)
                     call error_ring_push(trim(parse_errors(e)))
                     call logger%error(trim(parse_errors(e)))
                  end do
               end if
               ierr = OCEAN_STATUS_ERR_CONFIG_PARSE
               return
            end if
         end block
      else
         call sch%parse(trim(filename))
      end if

      call apply_time_unit_cascade(cfg)
      call apply_cartesian_degrees(cfg)
      if (present(ierr)) ierr = OCEAN_STATUS_OK

   end subroutine read_config_impl

   subroutine read_config_from_string(text, cfg, schema, ierr)
      !! In-memory sibling of [[read_config]]: build + strict-parse +
      !! apply a namelist held entirely in the `text` buffer (newline-
      !! separated), with NO filesystem touch — no temp file, no `chdir`,
      !! no cleanup.  Same validation behaviour as the file path (the
      !! schema parser and the native &ocean_bc_nml read both run from
      !! the in-memory line array); same `error stop` default and same
      !! optional non-aborting `ierr` as [[read_config]].
      character(len=*), intent(in) :: text
         !! Whole namelist as one string (records separated by '\n').
      type(config_t), target, intent(out) :: cfg
      type(nml_schema_t), intent(out), optional :: schema
      integer, intent(out), optional :: ierr
         !! Non-zero on a strict-parse failure when present; absent
         !! behaves as today (`error stop`). The natural FFI/create()
         !! door: a bad namelist string becomes a returned status
         !! instead of aborting the host process.

      if (present(schema)) then
         call read_config_from_string_impl(text, cfg, schema, ierr=ierr)
      else
         block
            type(nml_schema_t) :: sch
            call read_config_from_string_impl(text, cfg, sch, ierr=ierr)
         end block
      end if
   end subroutine read_config_from_string

   subroutine read_config_from_string_impl(text, cfg, sch, ierr)
      !! Body of [[read_config_from_string]] — mirrors read_config_impl
      !! with the line array replacing the file unit.
      character(len=*), intent(in) :: text
      type(config_t), target, intent(out) :: cfg
      type(nml_schema_t), intent(out) :: sch
      integer, intent(out), optional :: ierr
         !! Non-zero on a strict-parse failure when present; absent
         !! behaves as today (`error stop`).

      character(len=:), allocatable :: lines(:)
      integer :: n_lines

      call build_rdb_schema(cfg, sch)

      ! Empty buffer → keep pristine defaults (mirrors the missing-file
      ! branch of read_config_impl).
      if (len_trim(text) == 0) then
         if (present(ierr)) ierr = OCEAN_STATUS_OK
         return
      end if

      call split_to_lines(text, lines, n_lines)

      ! Strict schema parse from the same lines: validate + APPLY every
      ! migrated group straight into `cfg`.  Same status/errors routing
      ! as read_config_impl — see its comment for the rationale.
      if (present(ierr)) then
         block
            integer :: parse_status, e
            character(len=:), allocatable :: parse_errors(:)
            call sch%parse_lines(lines, n_lines, status=parse_status, errors=parse_errors)
            if (parse_status /= 0) then
               if (allocated(parse_errors)) then
                  do e = 1, size(parse_errors)
                     call error_ring_push(trim(parse_errors(e)))
                     call logger%error(trim(parse_errors(e)))
                  end do
               end if
               ierr = OCEAN_STATUS_ERR_CONFIG_PARSE
               return
            end if
         end block
      else
         call sch%parse_lines(lines, n_lines)
      end if

      call apply_time_unit_cascade(cfg)
      call apply_cartesian_degrees(cfg)
      if (present(ierr)) ierr = OCEAN_STATUS_OK
   end subroutine read_config_from_string_impl

   subroutine split_to_lines(text, lines, n_lines)
      !! Split a newline-separated buffer into a fixed-len character array
      !! (one record per line, trailing CR stripped) for internal-file
      !! reads + the schema line walker.
      character(len=*), intent(in) :: text
      character(len=:), allocatable, intent(out) :: lines(:)
      integer, intent(out) :: n_lines

      integer :: i, start, maxlen, cur, nl
      character(len=1) :: c

      ! First pass: count records + longest line.
      nl = 1
      maxlen = 0
      cur = 0
      do i = 1, len(text)
         c = text(i:i)
         if (c == achar(10)) then
            nl = nl + 1
            maxlen = max(maxlen, cur)
            cur = 0
         else if (c /= achar(13)) then
            cur = cur + 1
         end if
      end do
      maxlen = max(maxlen, cur, 1)

      allocate (character(len=maxlen) :: lines(nl))
      ! Blank-fill via a SECTION, not `lines = ""`.  `lines` has a DEFERRED
      ! length, so whole-variable intrinsic assignment from a zero-length
      ! scalar re-allocates it with len = 0 (F2018 10.2.1.3p3) — ifx does
      ! exactly that, and every record then parses as empty.  A section
      ! designator is not an allocatable variable, so it blank-pads in place.
      lines(:) = ""

      ! Second pass: copy each record.
      n_lines = 0
      start = 1
      cur = 0
      do i = 1, len(text)
         c = text(i:i)
         if (c == achar(10)) then
            n_lines = n_lines + 1
            if (cur > 0) lines(n_lines) = text(start:start + cur - 1)
            start = i + 1
            cur = 0
         else if (c /= achar(13)) then
            ! CR is skipped (a trailing CR before LF is excluded because
            ! `cur` counts only content chars and `start` marks the first).
            if (cur == 0) start = i
            cur = cur + 1
         end if
      end do
      ! Trailing record (no final newline).
      n_lines = n_lines + 1
      if (cur > 0) lines(n_lines) = text(start:start + cur - 1)
   end subroutine split_to_lines

   subroutine apply_time_unit_cascade(cfg)
      !! Post-parse time-unit fixup shared by the file + string config
      !! paths.  `t_end`, `status_interval`, `ocean_diag%dt_out` are given
      !! in `time_unit` from `&time_nml`; multiply through to SI seconds.
      !! `dt_fixed`/`dt_max` stay in seconds; default `time_unit="s"`
      !! gives a 1× factor (bit-identical to before).
      type(config_t), intent(inout) :: cfg

      real(wp) :: time_factor

      select case (trim(cfg%time_unit))
      case ("s", "sec", "second", "seconds")
         time_factor = 1.0_wp
      case ("min", "minute", "minutes")
         time_factor = 60.0_wp
      case ("hr", "hour", "hours", "h")
         time_factor = 3600.0_wp
      case ("day", "days", "d")
         time_factor = 86400.0_wp
      case ("year", "years", "yr")
         ! Julian year (365.25 d) — matches MOM6 / UDUNITS.
         time_factor = 365.25_wp*86400.0_wp
      case default
         call logger%error("Invalid time_unit = '"//trim(cfg%time_unit)// &
                           "': must be one of s/min/hr/day/year")
         time_factor = 1.0_wp  ! fallthrough so validation reports it
      end select
      cfg%t_end = cfg%t_end*time_factor
      cfg%status_interval = cfg%status_interval*time_factor
      cfg%ocean%diag%dt_out = cfg%ocean%diag%dt_out*time_factor
   end subroutine apply_time_unit_cascade

   subroutine apply_cartesian_degrees(cfg)
      !! Post-parse Cartesian domain sizing (MOM6 GRID_CONFIG="cartesian" +
      !! AXIS_UNITS / LENLON / LENLAT).  When `&ocean_grid_nml len_lon`/`len_lat`
      !! are set (> 0) on an ocean Cartesian grid, DERIVE the uniform
      !! `&grid_nml dx`/`dy` (metres) from the domain extent, the physical cell
      !! count (`nx`/`ny` = MOM6 NIGLOBAL/NJGLOBAL) and `axis_units`, exactly as
      !! MOM6's `set_grid_metrics_cartesian`:
      !!   degrees : dx = rad_earth · len_lon · π/180 / nx  (arc length, no cos(lat))
      !!   km      : dx = 1000 · len_lon / nx
      !!   meters  : dx =        len_lon / nx
      !! and likewise dy from len_lat / ny.  This runs BEFORE grid init, so the
      !! metrics, barotropic CFL and Coriolis all see the derived metres.
      !! Default (len_lon <= 0) leaves `dx`/`dy` untouched ⇒ bit-identical.
      type(config_t), intent(inout) :: cfg

      real(wp) :: dx_scale, dy_scale

      ! Opt-in: only engage when the extent knobs are set.
      if (cfg%ocean%grid%len_lon <= 0.0_wp .and. cfg%ocean%grid%len_lat <= 0.0_wp) return

      if (trim(cfg%sim_type) /= "ocean") then
         call logger%error("len_lon/len_lat are &ocean_grid_nml knobs; "// &
                           "set sim_type='ocean' to use them")
         return
      end if
      if (trim(cfg%ocean%grid%grid_config) /= "cartesian") then
         call logger%error("len_lon/len_lat (MOM6 LENLON/LENLAT) apply to "// &
                           "grid_config='cartesian' only (got '"// &
                           trim(cfg%ocean%grid%grid_config)//"'); the spherical/"// &
                           "tripolar paths already read &grid_nml dx/dy as degrees")
         return
      end if
      if (cfg%ocean%grid%len_lon <= 0.0_wp .or. cfg%ocean%grid%len_lat <= 0.0_wp) then
         call logger%error("len_lon AND len_lat must both be > 0 to size the "// &
                           "Cartesian domain (got len_lon="// &
                           to_string(cfg%ocean%grid%len_lon)//", len_lat="// &
                           to_string(cfg%ocean%grid%len_lat)//")")
         return
      end if

      ! Per-unit metres-per-degree/km/metre scale (applied to len/N below).
      select case (trim(cfg%ocean%grid%axis_units))
      case ("meters", "meter", "m")
         dx_scale = 1.0_wp
         dy_scale = 1.0_wp
      case ("km", "kilometer", "kilometers")
         dx_scale = 1000.0_wp
         dy_scale = 1000.0_wp
      case ("degrees", "degree", "deg")
         ! Arc length on a sphere of radius rad_earth; NO cos(lat) — MOM6's
         ! Cartesian grid is a flat plane, unlike its spherical grid.
         dx_scale = cfg%ocean%grid%rad_earth*DEG2RAD
         dy_scale = cfg%ocean%grid%rad_earth*DEG2RAD
      case default
         call logger%error("Invalid axis_units = '"// &
                           trim(cfg%ocean%grid%axis_units)// &
                           "': must be meters/degrees/km")
         return
      end select

      cfg%dx = dx_scale*cfg%ocean%grid%len_lon/real(cfg%nx, wp)
      cfg%dy = dy_scale*cfg%ocean%grid%len_lat/real(cfg%ny, wp)

      call logger%info("Cartesian domain sized from extent (MOM6-style): len_lon="// &
                       to_string(cfg%ocean%grid%len_lon)//", len_lat="// &
                       to_string(cfg%ocean%grid%len_lat)//" "// &
                       trim(cfg%ocean%grid%axis_units)//" over "// &
                       to_string(cfg%nx)//"x"//to_string(cfg%ny)//" cells => dx="// &
                       to_string(cfg%dx)//" m, dy="//to_string(cfg%dy)//" m")
   end subroutine apply_cartesian_degrees

   subroutine validate_config(cfg, ierr)
      !! Validate configuration parameters after reading
      !!
      !! Logs errors for invalid values that would crash the solver,
      !! and warnings for suspicious but non-fatal settings.
      !!
      !! `ierr`, when present, returns `OCEAN_STATUS_ERR_CONFIG_VALIDATE`
      !! instead of `error stop`-ing on the first accumulated failure
      !! (every individual check still logs via `global_logger%error`
      !! unchanged); absent behaves as today.
      use rdb_ocean_lateral_mix, only: parse_lateral_closure, &
                                       lateral_closure_is_implemented, &
                                       lateral_closure_conflicts_smag_ah, &
                                       has_biharmonic_backstop, &
                                       leith_biharm_is_inert
      use rdb_ocean_horizontal_viscosity, only: aniso_mode_is_implemented
      use rdb_vcoord, only: parse_vcoord_type, vcoord_h_min_is_coherent
      use rdb_constants, only: VCOORD_SIGMA, VCOORD_ZSTAR, VCOORD_EULERIAN_Z, &
                               VCOORD_ZSIGMA, VCOORD_LAGRANGIAN, VCOORD_ZSTAR_SIGMA, &
                               VCOORD_ZSTAR_FULL, VCOORD_Z_FIXED, H_VANISHED
      use rdb_ocean_boundary_types, only: ocean_bc_type_from_string, OBC_PERIODIC, OBC_WALL, &
                                          OBC_INVALID
      use rdb_coriolis_adv, only: parse_pv_variant, pv_variant_is_implemented, &
                                  parse_pv_adv_scheme, pv_adv_scheme_is_implemented, &
                                  pv_adv_required_nghost
      use rdb_ocean_bottom_drag, only: parse_bdrag_variant, bdrag_variant_is_implemented
      use rdb_ocean_pressure_force, only: parse_opgf_variant, gprime_nz_is_supported
      use rdb_ocean_tidal_mixing, only: tidal_mixing_is_inert
      use rdb_ocean_pseudo_salt, only: pseudo_salt_conflicts_restore, &
                                       pseudo_salt_conflicts_ice, &
                                       pseudo_salt_needs_thermo_warning
      use rdb_ocean_surface_flux, only: sw_source_is_implemented
      use rdb_ocean_vmix, only: kpp_sw_method_is_implemented, &
                                bkgnd_henyey_conflicts_profile
      use rdb_eos, only: parse_tfreeze_set, TFREEZE_SET_INVALID
      use rdb_ocean_top_drag, only: parse_tdrag_variant, tdrag_variant_is_implemented, &
                                    TDRAG_LINEAR, TDRAG_QUADRATIC
      use rdb_ocean_cavity_melt, only: parse_cavity_exchange_law, parse_cavity_ice_mode, &
                                       CAVITY_LAW_INVALID, CAVITY_LAW_CONST_GAMMA, &
                                       CAVITY_LAW_HJ99, CAVITY_LAW_YUNG25, &
                                       CAVITY_ICE_INVALID, CAVITY_ICE_INSULATING, &
                                       CAVITY_ICE_ADV_DIFF, &
                                       parse_cavity_freshwater, parse_cavity_volume_comp, &
                                       CAVITY_FW_INVALID, CAVITY_FW_VIRTUAL, CAVITY_FW_MASS, &
                                       CAVITY_VC_INVALID, CAVITY_VC_NONE, CAVITY_VC_UNIFORM_OPEN
      use rdb_ocean_cavity, only: parse_cavity_draft_sign, CAVITY_SIGN_INVALID
      type(config_t), intent(in) :: cfg
      integer, intent(out), optional :: ierr
         !! Non-zero on any cross-knob semantic validation failure when
         !! present; absent behaves as today (`error stop`).

      logical :: has_error
      integer :: lateral_closure_code

      has_error = .false.

      ! Simulation regime
      if (trim(cfg%sim_type) /= "ocean") then
         call logger%error("Invalid sim_type = '"//trim(cfg%sim_type)// &
                           "': must be 'ocean' (the coastal regime was split out into "// &
                           "its own repository)")
         has_error = .true.
      end if

      ! Grid parameters
      if (cfg%nx < 1) then
         call logger%error("Invalid nx = "//to_string(cfg%nx)//": must be >= 1")
         has_error = .true.
      end if
      if (cfg%ny < 1) then
         call logger%error("Invalid ny = "//to_string(cfg%ny)//": must be >= 1")
         has_error = .true.
      end if
      if (cfg%nghost < 1) then
         call logger%error("Invalid nghost = "//to_string(cfg%nghost)//": must be >= 1")
         has_error = .true.
      end if
      if (cfg%dx <= 0.0_wp) then
         call logger%error("Invalid dx = "//to_string(cfg%dx)//": must be > 0")
         has_error = .true.
      end if
      if (cfg%dy <= 0.0_wp) then
         call logger%error("Invalid dy = "//to_string(cfg%dy)//": must be > 0")
         has_error = .true.
      end if

      ! Time parameters
      if (cfg%cfl <= 0.0_wp .or. cfg%cfl > 1.0_wp) then
         call logger%error("Invalid cfl = "//to_string(cfg%cfl)//": must be in (0, 1]")
         has_error = .true.
      end if
      if (cfg%t_end <= 0.0_wp) then
         call logger%error("Invalid t_end = "//to_string(cfg%t_end)//": must be > 0")
         has_error = .true.
      end if

      ! Multilayer parameters
      if (cfg%use_multilayer) then
         if (cfg%nz_layers < 1) then
            call logger%error("Invalid nz_layers = "//to_string(cfg%nz_layers)// &
                              ": must be >= 1 for multilayer")
            has_error = .true.
         end if
      end if

      ! Per-column stack-workspace bound.  FAIL LOUD, both regimes.
      !
      ! Every layered kernel (BPG, ALE remap, kappa-shear, Redi, the
      ! diag-remap fills) carries fixed-size `NZ_STACK_MAX` thread-local
      ! column arrays.  Overrunning them corrupts thread-local storage
      ! SILENTLY — wrong answers, no crash — so this must refuse the run,
      ! not warn.  It replaces a warning-only check that additionally
      ! never fired on the ocean path: no ocean namelist sets
      ! `use_multilayer` (it defaults .false.), so the ocean regime had no
      ! nz guard at all beyond the wavespeed- and Redi-specific ones.
      !
      ! The requirement is `nz + 1` — see `nz_stack_required`.
      if (cfg%use_multilayer .or. trim(cfg%sim_type) == "ocean") then
         if (.not. nz_stack_is_sufficient(cfg%nz_layers)) then
            call logger%error("nz_layers = "//to_string(cfg%nz_layers)// &
                              " needs NZ_STACK_MAX >= "// &
                              to_string(nz_stack_required(cfg%nz_layers))// &
                              " but this binary was compiled with NZ_STACK_MAX = "// &
                              to_string(NZ_STACK_MAX)// &
                              ". Per-column stack kernels (BPG, ALE remap, kappa-shear, "// &
                              "Redi, diag remap) would overrun thread-local storage and "// &
                              "silently produce wrong answers. Rebuild with "// &
                              "-DRDB_NZ_STACK_MAX="// &
                              to_string(nz_stack_required(cfg%nz_layers))//" (or larger).")
            has_error = .true.
         end if
      end if

      ! Ocean windowed-drain tracer_recon guard (Q6).  A SEPARATE knob from
      ! the coastal tracer_recon above — it lives in `&ocean_vmix_nml` and
      ! only takes effect on sim_type='ocean' (ignored on coastal, like all
      ! ocean knobs).  Reuses the coastal parse + per-rung nghost helpers
      ! but NOT the coastal support-status predicate (which rejects ocean by
      ! design).  Fail loud on pqm/unknown and on an insufficient nghost for
      ! the requested rung (weno5→3, weno7→4, weno9→5).  The default nghost
      ! is 3 (bumped from 2 — see the grid config default), so ppm/weno5 are
      ! satisfied by default and only weno7/9 need an explicit extra column.
      if (trim(cfg%sim_type) == "ocean") then
         block
            use rdb_recon_weno, only: parse_tracer_recon, &
                                      tracer_recon_required_nghost, &
                                      TRACER_RECON_PPM
            integer :: ocn_recon_code
            ocn_recon_code = parse_tracer_recon(trim(cfg%ocean%vmix%tracer_recon))
            if (ocn_recon_code == -2) then
               call logger%error("&ocean_vmix_nml tracer_recon='pqm' is not implemented; "// &
                                 "use ppm, weno5, weno7 or weno9")
               has_error = .true.
            else if (ocn_recon_code < 0) then
               call logger%error("&ocean_vmix_nml tracer_recon='"// &
                                 trim(cfg%ocean%vmix%tracer_recon)// &
                                 "' not recognised; valid: ppm, weno5, weno7, weno9")
               has_error = .true.
            else if (ocn_recon_code /= TRACER_RECON_PPM) then
               if (cfg%nghost < tracer_recon_required_nghost(ocn_recon_code)) then
                  call logger%error("&ocean_vmix_nml tracer_recon='"// &
                                    trim(cfg%ocean%vmix%tracer_recon)// &
                                    "' requires nghost >= "// &
                                    to_string(tracer_recon_required_nghost(ocn_recon_code))// &
                                    " but nghost = "//to_string(cfg%nghost))
                  has_error = .true.
               end if
            end if
         end block
      end if

      ! (PR-21) Shortwave source + boundary-layer coupling selectors.
      ! Fail-loud on unknown strings (PR-6 idiom); q_sw source requires
      ! the PR-12 component set (makes the host-side source branch total);
      ! the net-heat source with penetrating SW warns about the
      ! night-time negative-I0 hazard (§2.4 of the plan).
      if (trim(cfg%sim_type) == "ocean") then
         if (.not. sw_source_is_implemented(trim(cfg%ocean%thermo%sw_source))) then
            call logger%error("&ocean_thermo_nml sw_source='"// &
                              trim(cfg%ocean%thermo%sw_source)// &
                              "' not recognised; valid: net_heat, q_sw")
            has_error = .true.
         end if
         if (.not. kpp_sw_method_is_implemented(trim(cfg%ocean%thermo%kpp_sw_method))) then
            call logger%error("&ocean_thermo_nml kpp_sw_method='"// &
                              trim(cfg%ocean%thermo%kpp_sw_method)// &
                              "' not recognised; valid: all_sw, mxl_sw, lv1_sw")
            has_error = .true.
         end if
         if (trim(cfg%ocean%thermo%sw_source) == "q_sw" .and. &
             .not. cfg%ocean%forcing%enable_components) then
            call logger%error("&ocean_thermo_nml sw_source='q_sw' requires "// &
                              "&ocean_forcing_nml enable_components=.true. "// &
                              "(q_sw is allocated only with the component set)")
            has_error = .true.
         end if
         if (trim(cfg%ocean%thermo%sw_source) == "net_heat" .and. &
             cfg%ocean%thermo%sw_pen_frac > 0.0_wp) then
            call logger%warning("shortwave penetration is reading the NET heat flux "// &
                                "(&ocean_thermo_nml sw_source='net_heat', sw_pen_frac>0): "// &
                                "I0 < 0 wherever Q_heat < 0, which drives unphysical "// &
                                "negative irradiance down the two-band profile. Set "// &
                                "sw_source='q_sw' with &ocean_forcing_nml "// &
                                "enable_components=.true.")
         end if
      end if

      ! Sponge width
      if (cfg%sponge_width < 0) then
         call logger%error("Invalid sponge_width = "//to_string(cfg%sponge_width)// &
                           ": must be >= 0")
         has_error = .true.
      end if

      ! Warn about unknown BC strings (they silently default to WALL)
      call warn_unknown_bc(cfg%bc_west, "bc_west")
      call warn_unknown_bc(cfg%bc_east, "bc_east")
      call warn_unknown_bc(cfg%bc_south, "bc_south")
      call warn_unknown_bc(cfg%bc_north, "bc_north")

      ! `thickness_config` is consumed only by the ocean state builder's IC
      ! seed.  The coastal regimes lay their layers on a separate path and
      ! would silently ignore a non-default value, so fail loud instead.
      if (trim(cfg%thickness_config) /= "sigma" .and. &
          trim(cfg%sim_type) /= "ocean") then
         call logger%error("thickness_config = '"//trim(cfg%thickness_config)// &
                           "' is an ocean-path knob (sim_type = 'ocean'); the coastal "// &
                           "regimes seed their layers on a separate path and would "// &
                           "silently ignore it")
         has_error = .true.
      end if

      ! Retired coastal-legacy `&tracer_nml` linear-EOS quartet.  These
      ! four keys reach `tracer_t%eos_coeff`/`eos_ref` and nothing else;
      ! the C-grid ocean path's linear EOS reads `&ocean_ic_nml` alone.
      ! A knob that validates and silently does nothing is the bug, so
      ! moving any of them off its historical default is fatal and names
      ! the live replacement rather than being quietly ignored.
      if (cfg%alpha_T /= LEGACY_TRACER_ALPHA_T) then
         call logger%error("&tracer_nml alpha_T is RETIRED on the ocean path: it only ever "// &
                           "reached tracer_t%eos_coeff, which no ocean kernel reads, so "// &
                           "setting it changed nothing. Use &ocean_ic_nml alpha_T "// &
                           "(kg/m^3 per degC) instead, or delete the key.")
         has_error = .true.
      end if
      if (cfg%beta_S /= LEGACY_TRACER_BETA_S) then
         call logger%error("&tracer_nml beta_S is RETIRED on the ocean path: it only ever "// &
                           "reached tracer_t%eos_coeff, which no ocean kernel reads, so "// &
                           "setting it changed nothing. Use &ocean_ic_nml beta_S "// &
                           "(kg/m^3 per PSU) instead, or delete the key.")
         has_error = .true.
      end if
      if (cfg%T_ref /= LEGACY_TRACER_T_REF) then
         call logger%error("&tracer_nml T_ref is RETIRED on the ocean path: it only ever "// &
                           "reached tracer_t%eos_ref, which no ocean kernel reads, so "// &
                           "setting it changed nothing. Use &ocean_ic_nml T_ref (degC) "// &
                           "instead, or delete the key.")
         has_error = .true.
      end if
      if (cfg%S_ref /= LEGACY_TRACER_S_REF) then
         call logger%error("&tracer_nml S_ref is RETIRED on the ocean path: it only ever "// &
                           "reached tracer_t%eos_ref, which no ocean kernel reads, so "// &
                           "setting it changed nothing. Use &ocean_ic_nml S_ref (PSU) "// &
                           "instead, or delete the key.")
         has_error = .true.
      end if

      ! Linear-EOS stratified salinity IC: both ends must be set (the
      ! seed gates on `/= 0` for BOTH, mirroring T_init_surface/bottom),
      ! so exactly one non-zero is a silently-uniform column — the same
      ! class of bug the retired knobs above were.
      if ((cfg%S_init_surface /= 0.0_wp) .neqv. (cfg%S_init_bottom /= 0.0_wp)) then
         call logger%error("&tracer_nml S_init_surface / S_init_bottom must BOTH be "// &
                           "non-zero to build the linear S(z) IC (they gate together, "// &
                           "exactly as T_init_surface/T_init_bottom do); one alone is "// &
                           "ignored and the column stays uniform at initial_salinity.")
         has_error = .true.
      end if

      ! Freezing-point (liquidus) coefficient set.  Belt-and-braces on top
      ! of the `nml_enum allowed=` list, in the style of the
      ! `conc_config` check above — and NOT optional: `parse_tfreeze_set`
      ! deliberately has no default fallback, so a string that reaches
      ! `eos_apply_tfreeze_set` unrecognised would silently leave the
      ! sea-ice set in place.  At S = 34.5 the two shipped sets differ by
      ! ~0.03 degC, which is enough to flip the sign of an ice-shelf
      ! basal melt rate — a mistyped liquidus must stop the run.
      if (parse_tfreeze_set(cfg%ocean%eos%tfreeze_set) == TFREEZE_SET_INVALID) then
         call logger%error("&ocean_eos_nml tfreeze_set = '"// &
                           trim(cfg%ocean%eos%tfreeze_set)// &
                           "' is not recognised (allowed: seaice, isomip). "// &
                           "'seaice' is the SIS2/MOM6 sea-ice liquidus "// &
                           "(-0.054*S - 7.53e-8*p); 'isomip' is the ISOMIP+ "// &
                           "ice-shelf-cavity set (-0.0573*S + 0.0832 - 7.53e-8*p, "// &
                           "Asay-Davis et al. 2016 Table 4).")
         has_error = .true.
      end if

      ! Ocean diag manager
      if (trim(cfg%sim_type) == "ocean") then
         if (trim(cfg%ocean%diag%vgrid) /= "layer" .and. &
             trim(cfg%ocean%diag%vgrid) /= "z_fixed" .and. &
             trim(cfg%ocean%diag%vgrid) /= "sigma" .and. &
             trim(cfg%ocean%diag%vgrid) /= "zstar" .and. &
             trim(cfg%ocean%diag%vgrid) /= "density") then
            call logger%error("Invalid ocean_diag_vgrid = '"// &
                              trim(cfg%ocean%diag%vgrid)// &
                              "': must be 'layer', 'z_fixed', 'sigma', 'zstar', or 'density'")
            has_error = .true.
         end if
         if (trim(cfg%ocean%diag%vgrid) == "z_fixed" .and. &
             cfg%ocean%diag%n_z_levels <= 0) then
            call logger%error("ocean_diag_vgrid = 'z_fixed' requires "// &
                              "ocean_diag_n_z_levels > 0")
            has_error = .true.
         end if
         if (cfg%ocean%diag%n_z_levels < 0 .or. &
             cfg%ocean%diag%n_z_levels > MAX_OCEAN_DIAG_Z_LEVELS) then
            call logger%error("ocean_diag_n_z_levels = "// &
                              to_string(cfg%ocean%diag%n_z_levels)// &
                              " out of range [0, "// &
                              to_string(MAX_OCEAN_DIAG_Z_LEVELS)//"]")
            has_error = .true.
         end if
         ! Density bins have no auto-fill (unlike sigma/z*), so the global
         ! vgrid='density' selection OR any per-diagnostic ':density'/':rho'
         ! attribute in `diags` requires a valid rho_levels list.
         if (.not. diag_density_levels_ok(trim(cfg%ocean%diag%vgrid), &
                                          trim(cfg%ocean%diag%diags), &
                                          cfg%ocean%diag%n_rho_levels, &
                                          cfg%ocean%diag%rho_levels)) then
            call logger%error("&ocean_diag_nml vgrid='density' (or a ':density'/':rho' "// &
                              "entry in diags) requires n_rho_levels > 0, "// &
                              "rho_levels(1:n_rho_levels) strictly increasing, and "// &
                              "n_rho_levels <= "//to_string(MAX_OCEAN_DIAG_Z_LEVELS))
            has_error = .true.
         end if
         if (cfg%ocean%diag%enabled .and. cfg%ocean%diag%dt_out <= 0.0_wp) then
            call logger%error("ocean_diag_dt_out must be > 0 when "// &
                              "ocean_diag_enabled = .true.")
            has_error = .true.
         end if
         if (cfg%ocean%diag%enabled .and. cfg%ocean%diag%dt_out > cfg%t_end) then
            call logger%warning("ocean_diag_dt_out ("// &
                                to_string(cfg%ocean%diag%dt_out)//" s) > t_end ("// &
                                to_string(cfg%t_end)//" s): the diagnostic will "// &
                                "fire at most once — is dt_out in your time_unit "// &
                                "(e.g. days), not seconds?")
         end if
         if (trim(cfg%ocean%topo%topo_config) /= "flat" .and. &
             trim(cfg%ocean%topo%topo_config) /= "spoon" .and. &
             trim(cfg%ocean%topo%topo_config) /= "seamount" .and. &
             trim(cfg%ocean%topo%topo_config) /= "neverworld2" .and. &
             trim(cfg%ocean%topo%topo_config) /= "island" .and. &
             trim(cfg%ocean%topo%topo_config) /= "double_drake" .and. &
             trim(cfg%ocean%topo%topo_config) /= "isomip_plus" .and. &
             trim(cfg%ocean%topo%topo_config) /= "file") then
            call logger%error("Invalid topo_config = '"//trim(cfg%ocean%topo%topo_config)// &
                              "': must be 'flat', 'spoon', 'seamount', 'neverworld2', "// &
                              "'island', 'double_drake', 'isomip_plus', or 'file'")
            has_error = .true.
         end if
         if (trim(cfg%ocean%ic%ic_config) /= "" .and. &
             trim(cfg%ocean%ic%ic_config) /= "eady" .and. &
             trim(cfg%ocean%ic%ic_config) /= "geostrophic_adjustment" .and. &
             trim(cfg%ocean%ic%ic_config) /= "baroclinic_jet") then
            call logger%error("Invalid ic_config = '"//trim(cfg%ocean%ic%ic_config)// &
                              "': must be '', 'eady', 'geostrophic_adjustment', or "// &
                              "'baroclinic_jet'")
            has_error = .true.
         end if
         if (trim(cfg%ocean%topo%topo_config) == "file" .and. &
             len_trim(cfg%bathymetry_file) == 0) then
            call logger%error("topo_config = 'file' requires bathymetry_file "// &
                              "to be set in &output_nml")
            has_error = .true.
         end if
         if (trim(cfg%ocean%topo%wind_config) /= "constant" .and. &
             trim(cfg%ocean%topo%wind_config) /= "2gyre" .and. &
             trim(cfg%ocean%topo%wind_config) /= "neverworld2") then
            call logger%error("Invalid wind_config = '"//trim(cfg%ocean%topo%wind_config)// &
                              "': must be 'constant', '2gyre', or 'neverworld2'")
            has_error = .true.
         end if
         if (cfg%ocean%topo%max_depth <= 0.0_wp) then
            call logger%error("ocean_max_depth must be > 0")
            has_error = .true.
         end if
         if (trim(cfg%ocean%topo%topo_config) == "spoon") then
            if (cfg%ocean%topo%edge_depth <= 0.0_wp) then
               call logger%error("ocean_edge_depth must be > 0 for spoon bathymetry")
               has_error = .true.
            end if
            if (cfg%ocean%topo%edge_depth >= cfg%ocean%topo%max_depth) then
               call logger%error("ocean_edge_depth must be < ocean_max_depth")
               has_error = .true.
            end if
            if (cfg%ocean%topo%slope_scale <= 0.0_wp) then
               call logger%error("ocean_slope_scale must be > 0 for spoon bathymetry")
               has_error = .true.
            end if
         end if

         ! `uniform_z` seeds collapsed layers at the isopycnal `angstrom_h`
         ! floor, whereas the wet/dry path pins the emerged-column invariant
         ! `Sum h_layer = nz*2*H_VANISHED` and floors `bt_h` to exactly that
         ! same value.  Mixing the two would seed a bt_h that disagrees with
         ! Sum h_layer on every emerged column, so refuse the combination.
         if (trim(cfg%thickness_config) == "uniform_z" .and. &
             cfg%ocean%wetdry%enable) then
            call logger%error("thickness_config = 'uniform_z' is incompatible with "// &
                              "&ocean_wetdry_nml enable = .true.: the emerged-column "// &
                              "seed invariant (bt_h = nz*2*H_VANISHED) assumes the "// &
                              "'sigma' even split")
            has_error = .true.
         end if

         ! `VCOORD_ZSIGMA` is NOT a working coordinate on the ocean path.
         ! Its deep branch reads `z_ref_global` as a table of absolute
         ! reference depths in METRES
         ! (`z_top_k = min(z_ref_global(nz-k), column_total)`,
         ! `rdb_ocean_vcoord :: ocean_vcoord_compute_target_h_impl`), but the
         ! ONLY writer of that array anywhere in `src/` is the DIMENSIONLESS
         ! `z_ref_global(k) = k/nz` init in `ocean_vcoord_init` — nothing on
         ! the namelist path or the Python path ever replaces it with metres.
         ! So every z-level interval is `1/nz` metres, every layer collapses,
         ! and the whole column is dumped into `target_h(:,:,1)` (the BED
         ! layer) by the deficit line.  `Sum target_h = H + eta` still holds,
         ! which is exactly why no conservation test ever caught it; the
         ! placement is measured interface-by-interface in
         ! `test_ocean_vcoord_interface_depths ::
         ! documents_zsigma_dimensionless_zref_collapse` (nine 0.1 m layers
         ! in the top 90 cm of a 1000 m column).
         !
         ! Refuse it rather than silently running a broken coordinate.  No
         ! shipped namelist selects it.  Follow-up: fill `z_ref_global` in
         ! metres (the family is the natural seat for a sigma-near-the-top /
         ! z-below HYBRID), then delete this refusal and the `documents_*`
         ! test with it.  `VCOORD_ZSTAR_SIGMA` consumes the same table
         ! FRACTIONALLY (rescaled by `z_ref_global(nz)`) so it is unaffected
         ! by the units — but note that with a uniform table its deep branch
         ! is numerically indistinguishable from SIGMA.
         if (parse_vcoord_type(cfg%vcoord_type, default_code=VCOORD_EULERIAN_Z) &
             == VCOORD_ZSIGMA) then
            call logger%error("&vcoord_nml vcoord_type = 'zsigma' is refused on the "// &
                              "ocean path: its deep branch reads `z_ref_global` as "// &
                              "absolute depths in METRES, but the only writer of that "// &
                              "table is the dimensionless `k/nz` init, so every "// &
                              "z-level interval is 1/nz metres and the whole column "// &
                              "collapses into the bed layer (the column sum is still "// &
                              "exact, which is why it looked healthy).  Use 'sigma', "// &
                              "'zstar', 'zstar_sigma' or 'zstar_full'; ZSIGMA returns "// &
                              "when `z_ref_global` is filled in metres.")
            has_error = .true.
         end if

         ! `&vcoord_nml zstar_h_min` carries TWO different contracts, picked
         ! by the coordinate family rather than by the value (see
         ! `rdb_vcoord :: vcoord_h_min_role`).  On the GEOMETRIC families
         ! (ZSTAR_FULL / Z_FIXED) it is the anti-zero thickness of below-bed
         ! FILLER layers that are meant to read as vanished downstream, so it
         ! belongs at or below the D4 skip/merge marker `H_VANISHED`; on the
         ! DENSITY families (RHO / HYCOM) the collapsed layers carry real
         ! tracer mass and that path deliberately floors at
         ! `max(zstar_h_min, 2*H_VANISHED)` instead.  Nothing else pins the
         ! knob, so an INERT-role run with `zstar_h_min > H_VANISHED` silently
         ! promotes its below-bed filler to LIVE layers (real EOS density off
         ! ghost T/S, a PGF column entry, a remap-drain concentration, a vdiff
         ! interface) while the coordinate still treats them as throwaway.
         !
         ! This was a WARNING until the rigid-top work made it blocking.  The
         ! reason it is now an ERROR: a coordinate that vanishes layers
         ! against the TOP of the column (an ice-shelf cavity) puts its
         ! fillers where the surface fluxes, the pressure gradient, the melt
         ! sampler and the tracer budgets all read — so a filler that is not
         ! skipped is not a cosmetic slip, it is a conservation hole.  The
         ! only configuration in the tree that sat in the band was the repo's
         ! own Python worked example (`python/tests/test_worked_example.py`,
         ! `ZStarFull(h_min=1.0e-3)`), corrected in the same change.
         !
         ! NOTE the boundary this must NOT move: five shipped namelists set
         ! `zstar_h_min = 1.5e-4`, H_VANISHED EXACTLY, which is legal —
         ! `vcoord_h_min_is_coherent` is a strict `>` and every downstream
         ! vanish test is a strict `> H_VANISHED`, so a layer on the marker
         ! reads as vanished.  Relaxing either to `>=` would refuse the
         ! canonical double-gyre reference; `test_ocean_vcoord_hygiene ::
         ! h_min_on_the_marker_is_accepted` guards that.
         !
         ! A non-positive floor was already refused and still is — it defeats
         ! the knob's single documented purpose.
         block
            integer :: hmin_vcoord_code
            hmin_vcoord_code = parse_vcoord_type(cfg%vcoord_type, &
                                                 default_code=VCOORD_EULERIAN_Z)
            if (.not. vcoord_h_min_is_coherent(hmin_vcoord_code, cfg%zstar_h_min)) then
               if (cfg%zstar_h_min <= 0.0_wp) then
                  call logger%error("&vcoord_nml zstar_h_min = "// &
                                    to_string(cfg%zstar_h_min)//" must be > 0: it exists "// &
                                    "so a vanishing layer's target thickness is never "// &
                                    "exactly zero (kernels that divide by h_layer)")
               else
                  call logger%error("&vcoord_nml zstar_h_min = "// &
                                    to_string(cfg%zstar_h_min)//" m exceeds H_VANISHED = "// &
                                    to_string(H_VANISHED)//" m under vcoord_type = '"// &
                                    trim(cfg%vcoord_type)//"': that family uses the knob "// &
                                    "as an anti-zero floor for filler layers that are "// &
                                    "MEANT to stay vanished — above H_VANISHED they "// &
                                    "become dynamically live (EOS/PGF/remap-drain/vdiff) "// &
                                    "while the coordinate still treats them as "// &
                                    "throwaway.  Use a value <= "// &
                                    to_string(H_VANISHED)//"; for a genuinely live "// &
                                    "minimum layer thickness use &ocean_isopycnal_nml "// &
                                    "angstrom_h (the D4 floor knob); the rho/hycom "// &
                                    "regrid has its own keep-alive floor")
               end if
               has_error = .true.
            end if
            ! NOTE the boundary, deliberately NOT warned about at runtime:
            ! the five shipped namelists set `zstar_h_min = 1.5e-4`, which is
            ! H_VANISHED EXACTLY.  That is legal — every downstream vanish
            ! test is a strict `> H_VANISHED`, so a layer sitting on the
            ! marker still reads as vanished — but it carries zero margin,
            ! and any gate relaxed to `>=` would change those runs' answers.
            ! A per-run warning here would fire on the canonical double-gyre
            ! reference forever while recommending nothing an operator can do
            ! without changing answers, so the fact lives in the docs and at
            ! the gate (`rdb_ocean_remap :: H_FLOOR`) instead.
         end block

         ! GM thickness diffusion consumes the stored isopycnal slope, so
         ! the slopes slot MUST be enabled.  Loud invariant (the kernel is
         ! pure/device and cannot fail loud); a silent no-op would hide a
         ! misconfigured GM run.
         if (cfg%ocean%gm%enable .and. .not. cfg%ocean%slopes%enable) then
            call logger%error("&ocean_gm_nml enable=.true. requires "// &
                              "&ocean_slopes_nml enable=.true. (GM reads the stored slope)")
            has_error = .true.
         end if
         ! khth_slope_max enters as 1/slope_max^2 in the safe-streamfunction
         ! limiter; a zero/negative value divides by zero.
         if (cfg%ocean%gm%enable .and. cfg%ocean%gm%khth_slope_max <= 0.0_wp) then
            call logger%error("&ocean_gm_nml khth_slope_max must be > 0 "// &
                              "(enters the safe-streamfunction limiter as 1/slope_max^2)")
            has_error = .true.
         end if

         ! Redi (capability [3]) ships the CONTINUOUS variant only; the
         ! discontinuous (regula-falsi) path is deferred (R4).  Loud
         ! invariant — the device kernel cannot fail loud.
         if (cfg%ocean%redi%enable .and. .not. cfg%ocean%redi%continuous) then
            call logger%error("&ocean_redi_nml continuous=.false. (discontinuous "// &
                              "variant) is not implemented yet (deferred R4)")
            has_error = .true.
         end if

         ! VarMix (capability [4]) needs the stored slope + N² (slopes slot)
         ! for the Eady term AND the first-mode cg1 (wavespeed slot) for the
         ! resolution function.  Loud invariant — the kernel is pure/device.
         if (cfg%ocean%varmix%enable .and. .not. cfg%ocean%slopes%enable) then
            call logger%error("&ocean_varmix_nml enable=.true. requires "// &
                              "&ocean_slopes_nml enable=.true. (VarMix reads the stored slope + N^2)")
            has_error = .true.
         end if
         if (cfg%ocean%varmix%enable .and. .not. cfg%ocean%wavespeed%enable) then
            call logger%error("&ocean_varmix_nml enable=.true. requires "// &
                              "&ocean_wavespeed_nml enable=.true. (the resolution function needs cg1)")
            has_error = .true.
         end if

         ! Resolution-scaled momentum viscosity (Gap 1) reads the VarMix
         ! resolution function `Res_fn`, so VarMix MUST be enabled.  Loud
         ! invariant — the lateral-mix kernel is pure/device.
         if (cfg%ocean%hvisc%resoln_scaled_visc .and. .not. cfg%ocean%varmix%enable) then
            call logger%error("&ocean_hvisc_nml resoln_scaled_visc=.true. requires "// &
                              "&ocean_varmix_nml enable=.true. (the resolution function lives in VarMix)")
            has_error = .true.
         end if

         ! MEKE (capability [5]) sources its eddy energy from the GM PE
         ! release (`gm%gm_src`), so GM MUST be enabled.  Loud invariant —
         ! the kernel is pure/device and cannot fail loud.  (The VarMix
         ! feedback seam is OPTIONAL: with VarMix off, MEKE still evolves E
         ! but `meke%kh` has no face accumulator to feed — documented, not an
         ! error.)
         if (cfg%ocean%meke%enable .and. .not. cfg%ocean%gm%enable) then
            call logger%error("&ocean_meke_nml enable=.true. requires "// &
                              "&ocean_gm_nml enable=.true. (MEKE sources from gm_src)")
            has_error = .true.
         end if

         ! ---- Lateral-mixing closure: fail loud on unimplemented tags ----
         ! `parse_lateral_closure` maps the &ocean_hvisc_nml string to an
         ! LMIX_* code; a mistyped/garbage string yields LMIX_INVALID and a
         ! tag with no dispatcher path is not implemented.  Either case must
         ! ABORT here — never silently fall through to background-only
         ! viscosity (the original silent-wrong-answer foot-gun).
         lateral_closure_code = parse_lateral_closure(cfg%ocean%hvisc%lateral_closure)
         if (.not. lateral_closure_is_implemented(lateral_closure_code)) then
            call logger%error("Invalid/unimplemented &ocean_hvisc_nml lateral_closure = '"// &
                              trim(cfg%ocean%hvisc%lateral_closure)//"': must be one of "// &
                              "'none', 'leith', 'smagorinsky'/'smag', 'biharmonic', "// &
                              "'leith_biharm' — a closure with no kernel must not "// &
                              "silently disable lateral viscosity")
            has_error = .true.
         end if
         ! `leith_biharm` and `smag_ah` are BOTH flow-aware biharmonic
         ! closures that fill the same `nu4_face_*` arrays; the dispatcher
         ! runs the closure first, then `compute_smag_ah` unconditionally,
         ! so an enabled `smag_ah` would silently OVERWRITE the
         ! Leith-biharmonic fill.  Fail loud rather than let one closure
         ! silently win (MOM6 max-combines them; we do not).
         if (lateral_closure_conflicts_smag_ah(lateral_closure_code, &
                                               cfg%ocean%hvisc%smag_ah)) then
            call logger%error("&ocean_hvisc_nml lateral_closure='leith_biharm' and "// &
                              "smag_ah=.true. are both flow-aware biharmonic closures "// &
                              "filling nu4_face — smag_ah would overwrite the "// &
                              "Leith-biharmonic fill; enable only one")
            has_error = .true.
         end if
         ! Only anisotropy mode 0 (grid-relative `aniso_dir`) has an
         ! implemented direction tensor; a requested-but-unimplemented mode
         ! must abort, not silently fall back to the grid-i default.
         if (.not. aniso_mode_is_implemented(cfg%ocean%hvisc%aniso_mode)) then
            call logger%error("&ocean_hvisc_nml aniso_mode must be 0 "// &
                              "(grid-relative aniso_dir) — other modes are not "// &
                              "implemented and must not silently fall back to grid-i")
            has_error = .true.
         end if
         ! MEKE backscatter injects a NEGATIVE harmonic viscosity that
         ! amplifies grid-scale modes; only a POSITIVE biharmonic
         ! (nu_4 / smag_ah / leith_biharm) can dissipate them.  Abort if
         ! backscatter is on without a biharmonic backstop configured —
         ! the CFL floor bounds the growth rate, it does not make a
         ! negative harmonic operator stable on its own.
         if (cfg%ocean%meke%backscatter .and. &
             .not. has_biharmonic_backstop(cfg%ocean%hvisc%nu_4, &
                                           cfg%ocean%hvisc%smag_ah, &
                                           cfg%ocean%hvisc%smag_bi_const, &
                                           lateral_closure_code, &
                                           cfg%ocean%hvisc%c_leith_bi, &
                                           cfg%ocean%hvisc%nu_4_bg)) then
            call logger%error("&ocean_meke_nml backscatter=.true. requires a "// &
                              "biharmonic backstop with a NON-ZERO coefficient "// &
                              "(&ocean_hvisc_nml nu_4>0 with no flow-aware "// &
                              "closure selected, or smag_ah=.true. with "// &
                              "smag_bi_const>0, or lateral_closure='leith_biharm' "// &
                              "with c_leith_bi>0 — or nu_4_bg>0 as a floor in "// &
                              "either flow-aware case) — the negative harmonic "// &
                              "backscatter is unstable without one")
            has_error = .true.
         end if

         ! ---- PR-6 fail-loud dispatch pack ----
         ! Each of these string→enum / config→kernel seams previously
         ! accepted a value the schema (or a hand-written reader) advertised,
         ! then silently ran DIFFERENT physics than the name promised.  Each
         ! guard consumes a `pure` predicate next to its dispatcher and
         ! aborts once via `has_error` — same idiom as the lateral-closure
         ! guard above.  Every guard is inert on valid config (bit-identical).
         ! Coriolis form: a typo (→ PV_VARIANT_INVALID) or the
         ! reserved-but-unwired `al81` (→ PV_VARIANT_AL81) must NOT silently
         ! run enstrophy-only Sadourny — a different conservation law.
         if (.not. pv_variant_is_implemented( &
             parse_pv_variant(cfg%ocean%coriolis%form))) then
            call logger%error("&ocean_coriolis_nml form = '"// &
                              trim(cfg%ocean%coriolis%form)//"' is not implemented — "// &
                              "must be 'sadourny', 'sadourny_hk' or 'sadourny_energy' "// &
                              "('al81' is reserved but its Arakawa-Lamb kernel is not "// &
                              "yet wired)")
            has_error = .true.
         end if
         ! PV face-interpolation scheme: a typo (→ PV_ADV_INVALID) must not
         ! silently fall back to centred.
         if (.not. pv_adv_scheme_is_implemented( &
             parse_pv_adv_scheme(cfg%ocean%coriolis%pv_adv_scheme))) then
            call logger%error("&ocean_coriolis_nml pv_adv_scheme = '"// &
                              trim(cfg%ocean%coriolis%pv_adv_scheme)//"' is not recognised — "// &
                              "must be 'centered', 'weno3', 'weno5' or 'weno7'")
            has_error = .true.
         end if
         ! weno5/weno7 stencil radius (3/4) needs a wider halo — mirror the
         ! tracer-WENO ladder's per-rung nghost gate (fail-loud, not a silent
         ! near-boundary order collapse over the whole domain).
         if (cfg%nghost < pv_adv_required_nghost( &
             parse_pv_adv_scheme(cfg%ocean%coriolis%pv_adv_scheme))) then
            call logger%error("&ocean_coriolis_nml pv_adv_scheme='"// &
                              trim(cfg%ocean%coriolis%pv_adv_scheme)//"' requires nghost >= "// &
                              to_string(pv_adv_required_nghost( &
                                        parse_pv_adv_scheme(cfg%ocean%coriolis%pv_adv_scheme)))// &
                              " but nghost = "//to_string(cfg%nghost))
            has_error = .true.
         end if
         ! WENO PV reconstruction is wired only into the Sadourny enstrophy
         ! path; pairing it with the hk/energy forms would silently run the
         ! centred interpolation those kernels hard-code.
         if (trim(cfg%ocean%coriolis%pv_adv_scheme) /= "centered" .and. &
             trim(cfg%ocean%coriolis%form) /= "sadourny") then
            call logger%error("&ocean_coriolis_nml pv_adv_scheme='"// &
                              trim(cfg%ocean%coriolis%pv_adv_scheme)// &
                              "' is only wired into form='sadourny'; got form='"// &
                              trim(cfg%ocean%coriolis%form)//"'")
            has_error = .true.
         end if
         ! Mass-consistent CorAdCalc (use_state_fluxes): only the
         ! sadourny_energy transport form consumes uh/vh, and only the
         ! mom6 corrector has a predictor continuity solve to be
         ! consistent WITH — any other combination would silently run
         ! the recompute path while the namelist promised otherwise.
         if (cfg%ocean%coriolis%use_state_fluxes) then
            if (trim(cfg%ocean%coriolis%form) /= "sadourny_energy") then
               call logger%error("&ocean_coriolis_nml use_state_fluxes requires "// &
                                 "form='sadourny_energy' (the transport form is what "// &
                                 "consumes uh/vh); got form='"// &
                                 trim(cfg%ocean%coriolis%form)//"'")
               has_error = .true.
            end if
            if (trim(cfg%ocean%bt%split_scheme) /= "pred_corr") then
               call logger%error("&ocean_coriolis_nml use_state_fluxes requires "// &
                                 "&ocean_bt_nml split_scheme='pred_corr' (the corrector "// &
                                 "consumes the predictor chain's fluxes); got "// &
                                 "split_scheme='"//trim(cfg%ocean%bt%split_scheme)//"'")
               has_error = .true.
            end if
         end if
         ! BOUND_CORIOLIS clamps the ENERGY-scheme CAu (q·vh) to the
         ! (f+ζ)·v velocity-form range; it is wired into the energy impl
         ! only (MOM6 also applies it to the enstrophy scheme, but Roundabout's
         ! enstrophy path is a separate kernel — v1 covers the audit target).
         ! Silently ignoring it on another form would be a positivity-request
         ! foot-gun, so fail loud.
         if (cfg%ocean%coriolis%bound_coriolis .and. &
             trim(cfg%ocean%coriolis%form) /= "sadourny_energy") then
            call logger%error("&ocean_coriolis_nml bound_coriolis is implemented for "// &
                              "form='sadourny_energy' only (the energy-scheme q·vh blow-up "// &
                              "it cures); got form='"//trim(cfg%ocean%coriolis%form)//"'")
            has_error = .true.
         end if
         ! corner_h="mom6_area" is wired into the energy impl's Pass 2 only
         ! (MOM6 shares the q construction across schemes, but Roundabout's
         ! enstrophy/hk are separate kernels — v1 covers the audit target).
         if (trim(cfg%ocean%coriolis%corner_h) == "mom6_area" .and. &
             trim(cfg%ocean%coriolis%form) /= "sadourny_energy") then
            call logger%error("&ocean_coriolis_nml corner_h='mom6_area' is implemented for "// &
                              "form='sadourny_energy' only; got form='"// &
                              trim(cfg%ocean%coriolis%form)//"'")
            has_error = .true.
         end if
         ! Bottom-drag form: linear (τ=ρ·r·u) and quadratic (τ=ρ·C_d·|u|·u)
         ! obey different physics; a typo must abort, not silently pick the
         ! quadratic default.
         if (.not. bdrag_variant_is_implemented( &
             parse_bdrag_variant(cfg%ocean%bdrag%form))) then
            call logger%error("&ocean_bdrag_nml form = '"// &
                              trim(cfg%ocean%bdrag%form)//"' is not implemented — "// &
                              "must be 'linear'/'rayleigh' or 'quadratic'/'cd'")
            has_error = .true.
         end if
         ! gprime PGF is a 2-layer reduced-gravity form — it hard-writes
         ! ONLY the top+bottom layer, so nz/=2 silently zeros the PGF on the
         ! other layers.  Config-time check against nz_layers (not nz_ml).
         if (.not. gprime_nz_is_supported( &
             parse_opgf_variant(cfg%ocean%pgf%form), cfg%nz_layers)) then
            call logger%error("&ocean_pgf_nml form = 'gprime' requires "// &
                              "&nonhydrostatic_nml nz_layers = 2 (the reduced-gravity "// &
                              "form writes only the top + bottom layer; use "// &
                              "the default form='mont' for a general-nz PGF)")
            has_error = .true.
         end if
         ! Leith-biharmonic with c_leith_bi<=0 is a provable no-op (ν₄ is
         ! linear in c_leith_bi) — the user asked for biharmonic dissipation
         ! and got none.  Promoted from a configure-time warning to an abort.
         if (leith_biharm_is_inert(lateral_closure_code, cfg%ocean%hvisc%c_leith_bi)) then
            call logger%error("&ocean_hvisc_nml lateral_closure='leith_biharm' with "// &
                              "c_leith_bi <= 0 is inert (ν₄ is linear in c_leith_bi) — "// &
                              "set c_leith_bi > 0 or choose a different closure")
            has_error = .true.
         end if
         ! Tidal mixing enabled with no energy source (e_uniform<=0 and
         ! e_compute off) makes Kd ≡ 0 — the whole flux sweep runs for
         ! nothing.
         if (tidal_mixing_is_inert(cfg%ocean%tidal_mixing%enable, &
                                   cfg%ocean%tidal_mixing%e_uniform, &
                                   cfg%ocean%tidal_mixing%e_compute)) then
            call logger%error("&ocean_tidal_mixing_nml enable=.true. with e_uniform <= 0 "// &
                              "and e_compute=.false. supplies no energy — Kd is "// &
                              "identically zero; set e_uniform > 0 or e_compute=.true.")
            has_error = .true.
         end if
         ! OBC edge strings bypass the schema (external group), so a typo
         ! silently closed the boundary to a wall.  Reject INVALID per edge,
         ! naming which edge (the guard runs before configure_ocean_bc, so an
         ! invalid tag never reaches the setup path).
         if (ocean_bc_type_from_string(cfg%ocean%bc%west) == OBC_INVALID) then
            call logger%error("&ocean_bc_nml west = '"//trim(cfg%ocean%bc%west)// &
                              "' is not a recognised boundary type (wall/open/tidal/"// &
                              "nested/inflow/discharge/clamped/sponge/chapman/periodic/"// &
                              "tripolar_fold)")
            has_error = .true.
         end if
         if (ocean_bc_type_from_string(cfg%ocean%bc%east) == OBC_INVALID) then
            call logger%error("&ocean_bc_nml east = '"//trim(cfg%ocean%bc%east)// &
                              "' is not a recognised boundary type")
            has_error = .true.
         end if
         if (ocean_bc_type_from_string(cfg%ocean%bc%south) == OBC_INVALID) then
            call logger%error("&ocean_bc_nml south = '"//trim(cfg%ocean%bc%south)// &
                              "' is not a recognised boundary type")
            has_error = .true.
         end if
         if (ocean_bc_type_from_string(cfg%ocean%bc%north) == OBC_INVALID) then
            call logger%error("&ocean_bc_nml north = '"//trim(cfg%ocean%bc%north)// &
                              "' is not a recognised boundary type")
            has_error = .true.
         end if

         ! Along-coordinate tracer Laplacian (`&ocean_hdiff_nml kappa_h`):
         ! the explicit forward-Euler stability bound
         ! kappa_h*dt_therm*(1/dx^2+1/dy^2) <= 0.5 used to be checked HERE
         ! with the nominal `cfg%dx`/`cfg%dy` — broken on non-Cartesian
         ! grids (DEGREES there, not metres) and blind to the true minimum
         ! cell even on Cartesian.  MOVED to `ocean_stability_audit`
         ! (`rdb_ocean_stability_audit.F90`), which runs AFTER the real
         ! per-cell `ocean_metrics_t` arrays are built and works on every
         ! grid type.  See that module's docstring, or
         ! docs/CAPABILITIES_AND_LIMITATIONS.md, for the full check list.
      end if

      ! ---- Implicit vdiff stress/drag folding (mutual exclusions) ----
      ! `implicit_drag` (fold into the vdiff bed diagonal) is mutually
      ! exclusive with the legacy `&ocean_bdrag_nml implicit` split-apply
      ! path (both would damp the bed velocity ⇒ double drag), and the
      ! bed-only fold cannot represent HBBL-distributed drag (`hbbl > 0`)
      ! — that needs a per-layer rate (follow-up PR).  Fail loud.
      if (cfg%ocean%vdiff%implicit_drag .and. cfg%ocean%bdrag%implicit) then
         call logger%error("ocean_vdiff implicit_drag is mutually exclusive with "// &
                           "ocean_bdrag implicit (split-apply): set only one")
         has_error = .true.
      end if
      if (cfg%ocean%vdiff%implicit_drag .and. cfg%ocean%bdrag%hbbl > 0.0_wp) then
         call logger%error("ocean_vdiff implicit_drag does not yet support "// &
                           "HBBL-distributed drag (ocean_bdrag hbbl > 0); use the "// &
                           "split-apply path (ocean_bdrag implicit) for HBBL")
         has_error = .true.
      end if
      ! `implicit_top_drag` is the `k = nz` twin of the rule above, plus
      ! one of its own: the surface row is the row the WIND stress owns
      ! as a Neumann RHS, so the fold both adds a diagonal term and
      ! masks that RHS off on the faces the ice covers.  That is only
      ! meaningful with a top drag configured.
      if (cfg%ocean%vdiff%implicit_top_drag) then
         if (.not. cfg%ocean%tdrag%enable) then
            call logger%error("&ocean_vdiff_nml implicit_top_drag=.true. requires "// &
                              "&ocean_tdrag_nml enable=.true.  The fold consumes the "// &
                              "top-drag slot's lambda_top_u/v and its face cover "// &
                              "masks; with the slot disabled those are placeholder "// &
                              "arrays and the knob would silently do nothing except "// &
                              "look like a top drag was configured.")
            has_error = .true.
         end if
         if (cfg%ocean%tdrag%implicit) then
            call logger%error("&ocean_vdiff_nml implicit_top_drag is mutually "// &
                              "exclusive with &ocean_tdrag_nml implicit: both damp "// &
                              "the top layer, so running both is a DOUBLE COUNT, not "// &
                              "a stronger drag.  Pick one — the vdiff fold if the "// &
                              "column also has real vertical viscosity to couple "// &
                              "against, the in-kernel backward-Euler form otherwise.")
            has_error = .true.
         end if
         if (cfg%ocean%tdrag%htbl > 0.0_wp) then
            call logger%error("&ocean_vdiff_nml implicit_top_drag does not support "// &
                              "the HTBL-distributed top drag (&ocean_tdrag_nml htbl "// &
                              "> 0): the fold is a SINGLE k = nz Rayleigh rate on the "// &
                              "diagonal and cannot represent a band spread over "// &
                              "several layers (the mirror of the implicit_drag/HBBL "// &
                              "restriction).  Use &ocean_tdrag_nml implicit for a "// &
                              "distributed top drag.")
            has_error = .true.
         end if
      end if
      ! `bbl_glue` (MOM6 bottomdraglaw coupling parity, PGF_BUG.md §9)
      ! needs the harmonic-z bookkeeping that only the hvel_mom6 path
      ! builds, replaces the implicit-drag bed fold (so that fold must be
      ! on for there to be exactly one bed sink), and is a CONSTANT-piston
      ! parity — only the linear drag form maps onto it (the quadratic
      ! |U|-dependent piston composition is deferred).  Fail loud on each.
      if (cfg%ocean%vdiff%bbl_glue) then
         if (.not. cfg%ocean%vdiff%hvel_mom6) then
            call logger%error("ocean_vdiff bbl_glue requires hvel_mom6=.true. — the "// &
                              "botfn glue reads the harmonic height-above-bed stack "// &
                              "that only the hvel_mom6 path accumulates")
            has_error = .true.
         end if
         if (.not. cfg%ocean%vdiff%implicit_drag) then
            call logger%error("ocean_vdiff bbl_glue requires implicit_drag=.true. — "// &
                              "the glue's piston bed row REPLACES the implicit-drag "// &
                              "fold; without it the bed would carry no sink at all")
            has_error = .true.
         end if
         if (trim(cfg%ocean%bdrag%form) /= "linear") then
            call logger%error("ocean_vdiff bbl_glue requires ocean_bdrag form='linear' "// &
                              "(constant-piston parity; quadratic composition deferred)")
            has_error = .true.
         end if
      end if
      ! `implicit_stress` injects the wind stress at the surface row (k=nz)
      ! only — it cannot represent the DIRECT_STRESS distribution of stress
      ! over the top `hmix_stress` metres (`vmix%direct_stress`).  Enabling
      ! both would silently drop the distribution (the fold wins).  Symmetric
      ! to the implicit_drag/HBBL guard above.  Fail loud.
      if (cfg%ocean%vdiff%implicit_stress .and. cfg%ocean%vmix%direct_stress) then
         call logger%error("ocean_vdiff implicit_stress is incompatible with "// &
                           "distributed wind stress (ocean_vmix direct_stress): the "// &
                           "surface-row fold cannot spread stress over hmix_stress")
         has_error = .true.
      end if

      ! `correction_visc_rem` (PR-19) only has a reader in the h-weighted
      ! branch of apply_bt_correction — the uniform branch never touches
      ! `visc_rem_u/v`.  Without `correction_h_weighted`, the knob would
      ! validate, log, and silently do nothing.  Fail loud rather than
      ! ship a no-op configuration.
      if (cfg%ocean%bt%correction_visc_rem .and. .not. cfg%ocean%bt%correction_h_weighted) then
         call logger%error("&ocean_bt_nml correction_visc_rem=.true. requires "// &
                           "correction_h_weighted=.true. — the uniform-Δu BT-corrector "// &
                           "branch never reads visc_rem, so the knob would silently "// &
                           "do nothing")
         has_error = .true.
      end if
      ! The vdiff operator's row sums are exactly 1 (a no-flux-top/no-flux-
      ! bottom viscous operator cannot remove a uniform acceleration)
      ! UNLESS the implicit-drag fold breaks the k=1 row sum.  So without
      ! `implicit_drag`, the remnant producer still runs but returns
      ! gamma ≡ 1 identically, and the corrector reduces to the plain
      ! h-weighted path.  Legal and mathematically correct — merely inert.
      ! Warn (not error): same "enabled but inert" precedent as the tidal-
      ! mixing e_uniform=0 warning.
      ! `forcing_visc_rem` reads the same visc_rem arrays the corrector
      ! weighting does — without the producer knob they stay ≡ 1 and the
      ! forcing weighting silently degenerates to the plain h-mean.
      if (cfg%ocean%bt%forcing_visc_rem .and. .not. cfg%ocean%bt%correction_visc_rem) then
         call logger%error("&ocean_bt_nml forcing_visc_rem=.true. requires "// &
                           "correction_visc_rem=.true. — that knob is the visc_rem "// &
                           "producer; without it the weights are identically 1")
         has_error = .true.
      end if
      ! pred_corr (SPEC S4) envelope: any ALE / Lagrangian-within-step
      ! vcoord (lagrangian, sigma, zstar, zsigma, zstar_sigma,
      ! zstar_full) — MOM6's own model: the dynamics step is always
      ! Lagrangian, regridding is orthogonal, and the pc corrector flows
      ! into the same thermo-cadence ALE remap the ssp path uses.  The
      ! LEGACY pure Eulerian-z path stays excluded: its per-stage
      ! vertical advection + BT-fold h-rescale would compose differently
      ! inside the restructured loop (the fold's rescale lands BEFORE the
      ! deferred continuity — a double eta application) — untested, fail
      ! loud.  Also still excluded: dynamic wet/dry (per-stage masking
      ! composes with the ssp stages only) and windowed tracer advection
      ! (the predictor's TR_MODE_NONE + window accounting is unexercised).
      if (trim(cfg%ocean%bt%split_scheme) == "pred_corr") then
         if (parse_vcoord_type(cfg%vcoord_type, &
                               default_code=VCOORD_EULERIAN_Z) == VCOORD_EULERIAN_Z) then
            call logger%error("&ocean_bt_nml split_scheme='pred_corr' requires an ALE "// &
                              "vertical coordinate (lagrangian/sigma/zstar/zsigma/"// &
                              "zstar_sigma/zstar_full) — the legacy 'eulerian_z' "// &
                              "per-stage vertical-advection + h-rescale path is not "// &
                              "wired through the restructured pc loop; got '"// &
                              trim(cfg%vcoord_type)//"'. pred_corr is the DEFAULT, so "// &
                              "set split_scheme='ssp_rk2' to keep this configuration.")
            has_error = .true.
         end if
         if (cfg%ocean%wetdry%enable) then
            call logger%error("&ocean_bt_nml split_scheme='pred_corr' is incompatible "// &
                              "with &ocean_wetdry_nml enable (v1 envelope). pred_corr "// &
                              "is the DEFAULT, so set split_scheme='ssp_rk2' to keep "// &
                              "this configuration.")
            has_error = .true.
         end if
         if (cfg%ocean%vmix%dt_tracer_advect_ratio > 1) then
            call logger%error("&ocean_bt_nml split_scheme='pred_corr' requires "// &
                              "dt_tracer_advect_ratio = 1 (v1 envelope). pred_corr is "// &
                              "the DEFAULT, so set split_scheme='ssp_rk2' to keep this "// &
                              "configuration.")
            has_error = .true.
         end if
         ! z-level closed faces need their SOLID WALLS to be walls.  With
         ! `mask_wall_velocity = .false.` a wall face keeps a free per-layer
         ! velocity that carries no mass (continuity zeroes the wall flux)
         ! and whose depth mean the BT fold resets to zero every stage —
         ! but whose BAROCLINIC part nothing restores under pred_corr: the
         ! Coriolis reads `u_av`, which the transport renormaliser never
         ! writes at a skipped wall face, so the wall velocity never sees
         ! its own rotation and integrates the layer Coriolis of its
         ! interior neighbour without bound.  The closed-face mask is what
         ! makes that neighbour baroclinic (a closed bed layer beside open
         ! ones).  Measured on the rotating ledge basin of
         ! tests/test_ocean_zfixed_cor_ref.F90: KE+PE x1783 over 4000
         ! outer steps with the physical interior at x1.24 — the whole
         ! growth sits on the unmasked wall faces — against x0.74 with
         ! the walls masked.  The mask is the DEFAULT; refused rather than
         ! silently overridden.
         if (cfg%zfixed_closed_faces .and. .not. cfg%ocean%bc%mask_wall_velocity) then
            call logger%error("&ocean_bt_nml split_scheme='pred_corr' with "// &
                              "&vcoord_nml zfixed_closed_faces=.true. requires "// &
                              "&ocean_bc_nml mask_wall_velocity=.true. (the default): "// &
                              "an unmasked solid-wall face carries a baroclinic layer "// &
                              "velocity that pred_corr's u_av-evaluated Coriolis never "// &
                              "rotates, so it grows without bound. Set "// &
                              "mask_wall_velocity=.true., or split_scheme='ssp_rk2'.")
            has_error = .true.
         end if
      end if
      ! `accel_visc_rem` (the per-layer slow-apply attenuation) reads the
      ! same producer arrays — without the producer knob they stay ≡ 1
      ! and the reweight is a silent no-op.
      if (cfg%ocean%vdiff%accel_visc_rem .and. .not. cfg%ocean%bt%correction_visc_rem) then
         call logger%error("&ocean_vdiff_nml accel_visc_rem=.true. requires "// &
                           "&ocean_bt_nml correction_visc_rem=.true. — that knob is the "// &
                           "visc_rem producer; without it the weights are identically 1 "// &
                           "and the reweight silently does nothing")
         has_error = .true.
      end if
      ! `renorm_visc_rem` (SPEC S2b) reads the same producer arrays.
      if (cfg%ocean%bt%renorm_visc_rem .and. .not. cfg%ocean%bt%correction_visc_rem) then
         call logger%error("&ocean_bt_nml renorm_visc_rem=.true. requires "// &
                           "correction_visc_rem=.true. — that knob is the visc_rem "// &
                           "producer; without it the gamma weights are identically 1")
         has_error = .true.
      end if
      if (cfg%ocean%bt%correction_visc_rem .and. .not. cfg%ocean%vdiff%implicit_drag) then
         call logger%warning("&ocean_bt_nml correction_visc_rem=.true. but "// &
                             "&ocean_vdiff_nml implicit_drag=.false.: the vdiff operator "// &
                             "then carries no drag, so visc_rem = 1 identically and the "// &
                             "BT corrector reduces to the plain h-weighted path")
      end if
      if (substep_drag_ignores_bdrag_form(cfg)) then
         call logger%warning("&ocean_bt_nml substep_drag=.true. with "// &
                             "&ocean_bdrag_nml form='"// &
                             trim(adjustl(cfg%ocean%bdrag%form))//"': the "// &
                             "barotropic substep damping is built from the "// &
                             "LINEAR coefficient &ocean_bdrag_nml r (times hbbl) "// &
                             "only, so it does NOT follow this bottom drag — "// &
                             "with r = 0 (the default) substep_drag is a no-op, "// &
                             "otherwise it damps the barotropic mode with a "// &
                             "linear drag the slow step never applies.  Use "// &
                             "form='linear', or drop substep_drag")
      end if

      ! Equilibrium tide (C1) requires lat/lon — meaningless on a
      ! cartesian grid (geolatT/geolonT stay 0).  Fail loud rather than
      ! silently forcing a flat basin.
      if (cfg%ocean%tides%enable .and. &
          trim(cfg%ocean%grid%grid_config) == "cartesian") then
         call logger%error("ocean_tides_nml: enable=.true. requires a "// &
                           "non-cartesian grid_config (spherical/tripolar/"// &
                           "supergrid) — the equilibrium tide needs lat/lon")
         has_error = .true.
      end if
      if (cfg%ocean%tides%enable .and. len_trim(cfg%ocean%tides%constituents) == 0) then
         call logger%error("ocean_tides_nml: enable=.true. but constituents "// &
                           "list is empty")
         has_error = .true.
      end if

      ! C7 Henyey latitude factor: copy of the equilibrium-tide cartesian
      ! guard above.  `geolatT` is identically zero on a cartesian grid (only
      ! the spherical/supergrid/tripolar fills populate it), so EVERY column
      ! would take the equatorial factor L(0 deg) = 0 and the background
      ! would collapse to a uniform `bkgnd_kd_min` everywhere while the
      ! console cheerfully logged "Henyey IGW latitude factor ON".  The
      ! `kd_min` floor makes that less catastrophic than the unfloored
      ! zeroing it used to be, but it is still a latitude parameterisation
      ! on a grid with no meaningful latitude — strictly worse than leaving
      ! the knob off, so refuse at configure.
      if (cfg%ocean%vmix%bkgnd_henyey .and. &
          trim(cfg%ocean%grid%grid_config) == "cartesian") then
         call logger%error("&ocean_vmix_nml bkgnd_henyey=.true. requires a "// &
                           "non-cartesian &ocean_grid_nml grid_config (spherical/"// &
                           "tripolar/supergrid) — geolatT is meaningless on cartesian, "// &
                           "so every column would take the equatorial L(0 deg)=0 and "// &
                           "the background would flatten to a uniform bkgnd_kd_min")
         has_error = .true.
      end if
      ! Bryan-Lewis and Henyey are MUTUALLY EXCLUSIVE background schemes,
      ! matching the reference code, which FATALs when a second background
      ! scheme is selected.  Henyey scales the SCALAR background; letting it
      ! also multiply the Bryan-Lewis deep asymptote would suppress abyssal /
      ! internal-tide mixing at the equator, which is not what the Henyey
      ! scaling describes.
      if (bkgnd_henyey_conflicts_profile(cfg%ocean%vmix%bkgnd_henyey, &
                                         cfg%ocean%vmix%bkgnd_profile)) then
         call logger%error("&ocean_vmix_nml bkgnd_henyey and bkgnd_profile are "// &
                           "mutually exclusive background schemes — pick ONE "// &
                           "(Bryan-Lewis depth profile, or the Henyey latitude "// &
                           "factor on the scalar background)")
         has_error = .true.
      end if
      ! Surface-pressure loading / inverse barometer (PR-17).  The seam it
      ! folds into (`eta_forcing`) lives in the barotropic substep, so the
      ! feature is split-solver only; the `p_surf` field it reads is
      ! allocated only with the PR-12 component set.  Every violation fails
      ! loud rather than silently reading an unallocated array (a silent
      ! stale device read on the mem:separate GPU build).
      if (cfg%ocean%psurf%enable) then
         if (.not. cfg%ocean%forcing%enable_components) then
            call logger%error("&ocean_psurf_nml enable=.true. requires "// &
                              "&ocean_forcing_nml enable_components=.true. "// &
                              "(p_surf is allocated only with the component set)")
            has_error = .true.
         end if
         if (cfg%ocean%bt%n_inner < 1 .and. .not. cfg%ocean%bt%auto_n_inner) then
            call logger%error("&ocean_psurf_nml enable=.true. requires the "// &
                              "split-explicit solver (&ocean_bt_nml n_inner >= 1 "// &
                              "or auto_n_inner=.true.) — the eta_forcing seam "// &
                              "lives in the barotropic substep")
            has_error = .true.
         end if
         ! A UNIFORM surface pressure has no gradient => provably inert
         ! (gauge invariance).  A user enabling p_surf with only a uniform
         ! constant and no file/override path gets nothing — say so rather
         ! than repeat the &ocean_tidal_mixing_nml e_uniform silent-no-op.
         if (cfg%ocean%psurf%p_surf_const /= 0.0_wp) then
            if (cfg%ocean%psurf%in_eos) then
               ! With in_eos the gauge argument does NOT apply to the EOS:
               ! it is nonlinear in pressure, so a spatially UNIFORM load
               ! still changes the IN-SITU densities it reaches.  Say so
               ! instead of the (then wrong) "provably inert" warning.
               call logger%info("&ocean_psurf_nml in_eos=.true.: a uniform "// &
                                "p_surf_const is inert on the barotropic seam "// &
                                "(only grad(p_surf) is physical there) but NOT "// &
                                "in the IN-SITU EOS pressure — it compresses "// &
                                "the water. The potential density ms%rho_layer "// &
                                "is unaffected either way (uniform p_ref).")
            else
               call logger%warning("&ocean_psurf_nml enable=.true. with a uniform "// &
                                   "p_surf_const and no file/override path: only "// &
                                   "grad(p_surf) is physical, so a spatially "// &
                                   "constant load is inert (gauge invariance). "// &
                                   "File-driven p_surf lands in PR-14/PR-15.")
            end if
         end if
      end if
      ! ---- Top-of-column pressure in the IN-SITU EOS arguments (E3) ----
      ! `in_eos` offsets the EOS's IN-SITU pressure arguments by
      ! `ms%p_top`.  v1 ports exactly one: the FV_WRIGHT Picard column
      ! sweep.  It does NOT touch `ms%rho_layer`, which is a potential
      ! density at the horizontally uniform `&ocean_eos_nml p_ref` and must
      ! stay that way.  Every OTHER in-situ consumer builds its own
      ! surface-relative hydrostatic pressure starting at 0 Pa at the free
      ! surface and has NOT been ported; running them against a loaded
      ! column would mix two incompatible pressure conventions inside one
      ! time step with no symptom.  Refuse fail-loud rather than be
      ! silently inconsistent — each line below is a named follow-up, not a
      ! permanent limit.
      if (cfg%ocean%psurf%in_eos) then
         ! Inert-configuration warning, not a refusal: with no in-situ
         ! consumer selected the knob legitimately does nothing, and the
         ! house rule (cf. &ocean_tidal_mixing_nml e_uniform) is to SAY so
         ! rather than let a user believe a cavity load reached the EOS.
         ! E4 adds a THIRD ported consumer: `buoyancy_coeffs="eos"` seeds
         ! the KPP B_0 coefficients at `p_top` and the double-diffusion
         ! interface stack from it, so the knob is no longer inert when
         ! that is selected.
         if (trim(adjustl(cfg%ocean%pgf%form)) /= "fv_wright" .and. &
             .not. cfg%ocean%epbl%enable .and. &
             trim(adjustl(cfg%ocean%vmix%buoyancy_coeffs)) /= "eos") then
            call logger%warning("&ocean_psurf_nml in_eos=.true. is INERT for "// &
                                "&ocean_pgf_nml form='"// &
                                trim(adjustl(cfg%ocean%pgf%form))//"' with no "// &
                                "other ported in-situ consumer: the ported ones "// &
                                "are the FV_WRIGHT Picard column sweep, the "// &
                                "EPBL column stack (&ocean_epbl_nml) and the "// &
                                "EOS-derived buoyancy coefficients "// &
                                "(&ocean_vmix_nml buoyancy_coeffs='eos'). The "// &
                                "other PGF forms read the POTENTIAL density "// &
                                "ms%rho_layer, which is referenced to the "// &
                                "uniform &ocean_eos_nml p_ref BY DESIGN and is "// &
                                "not offset by the load.")
         end if
         if (.not. cfg%ocean%psurf%enable) then
            call logger%error("&ocean_psurf_nml in_eos=.true. requires "// &
                              "enable=.true. (p_top is filled from the "// &
                              "assembled sf%p_surf the seam owns)")
            has_error = .true.
         end if
         ! EPBL was refused here until Phase 4b.  It is PORTED now:
         ! `epbl_column_kernel` seeds its stack at `ms%p_top(i,j)` when
         ! `in_eos`, which moves BOTH consumers of that stack (the
         ! in-situ `eos_specvol_derivs` argument and the PE weight
         ! `dmass*p_mid*dsv`) together.  Gate:
         ! `test_ocean_bl_under_ice` (`epbl_p_top_*`).
         if (cfg%ocean%kshear%enable) then
            call logger%error("&ocean_psurf_nml in_eos=.true. is not supported "// &
                              "with &ocean_kappa_shear_nml enable=.true. — "// &
                              "ks_solve_column builds its own interface pressure "// &
                              "from 0 Pa at the surface; not yet ported to p_top")
            has_error = .true.
         end if
         if (cfg%ocean%tidal_mixing%enable) then
            call logger%error("&ocean_psurf_nml in_eos=.true. is not supported "// &
                              "with &ocean_tidal_mixing_nml enable=.true. — "// &
                              "tidal_mixing_column_kernel builds its own interface "// &
                              "pressure from 0 Pa at the surface; not yet ported")
            has_error = .true.
         end if
         if (cfg%ocean%redi%enable) then
            call logger%error("&ocean_psurf_nml in_eos=.true. is not supported "// &
                              "with &ocean_redi_nml enable=.true. — "// &
                              "redi_build_column seeds Pint(1)=0 at the surface; "// &
                              "not yet ported to p_top")
            has_error = .true.
         end if
         if (cfg%ocean%slopes%enable) then
            call logger%error("&ocean_psurf_nml in_eos=.true. is not supported "// &
                              "with &ocean_slopes_nml enable=.true. — "// &
                              "pressure_above_x sums g*rho0*h down from 0 Pa at "// &
                              "the surface; not yet ported to p_top")
            has_error = .true.
         end if
         if (cfg%ocean%pgf%reconstruct_for_pressure) then
            call logger%error("&ocean_psurf_nml in_eos=.true. is not supported "// &
                              "with &ocean_pgf_nml reconstruct_for_pressure=.true. "// &
                              "— boole_dpa_intz_layer builds the EOS pressure as "// &
                              "p = -g*rho0*z from the surface-relative interface "// &
                              "height; not yet ported to p_top")
            has_error = .true.
         end if
         if (cfg%ocean%ice%enable) then
            call logger%error("&ocean_psurf_nml in_eos=.true. is not supported "// &
                              "with &ocean_ice_nml enable=.true. — the freezing "// &
                              "point is evaluated at p=0 in the frazil / basal-flux "// &
                              "kernels, and the liquidus pressure depression is "// &
                              "exactly what an ice-shelf load changes; not yet ported")
            has_error = .true.
         end if
      end if
      ! ---- Z-level T/S initial-condition overlay (`&ocean_zinit_nml`) ----
      if (cfg%ocean%zinit%enable) then
         select case (trim(adjustl(cfg%ocean%zinit%source)))
         case ("file")
            if (len_trim(cfg%ocean%zinit%file) == 0) then
               call logger%error("&ocean_zinit_nml enable=.true. with "// &
                                 "source='file' requires a non-blank file= path.")
               has_error = .true.
            end if
         case ("linear")
            ! An analytic profile plus a file path is ambiguous: the file
            ! would be silently ignored.  Say so rather than pick one.
            if (len_trim(cfg%ocean%zinit%file) > 0) then
               call logger%error("&ocean_zinit_nml source='linear' takes the "// &
                                 "analytic lin_* profile and opens NOTHING, so the "// &
                                 "file='"//trim(adjustl(cfg%ocean%zinit%file))// &
                                 "' you also set would be silently ignored.  Pick "// &
                                 "one: drop file=, or set source='file'.")
               has_error = .true.
            end if
         case default
            call logger%error("&ocean_zinit_nml source='"// &
                              trim(adjustl(cfg%ocean%zinit%source))// &
                              "' is not a known profile source; expected 'file' "// &
                              "(pre-regridded NetCDF) or 'linear' (analytic "// &
                              "affine T(z)/S(z)).")
            has_error = .true.
         end select
      end if
      ! ---- Top-of-column load in the PGF surface boundary condition (P5.0) ----
      ! `pa(nz+1) = rho_ref*g*eta_geo + ms%p_top`.  Only the FV_MOM6 family
      ! builds a `pa` stack at all — `mont` hard-zeroes `M(nz)` and
      ! `fv_lite`/`fv_wright` seed `p_edge(nz+1) = 0` — so there is
      ! literally no boundary condition to inject anywhere else.  Refuse
      ! rather than accept a knob that would silently do nothing on a form
      ! the user believes is carrying an ice load.
      if (cfg%ocean%pgf%p_top_in_bc) then
         if (trim(adjustl(cfg%ocean%pgf%form)) /= "fv_mom6") then
            call logger%error("&ocean_pgf_nml p_top_in_bc=.true. requires "// &
                              "form='fv_mom6' (got '"// &
                              trim(adjustl(cfg%ocean%pgf%form))//"'). Only the "// &
                              "FV_MOM6 family builds the pa(nz+1) pressure-stack "// &
                              "boundary condition the load is injected into; mont "// &
                              "hard-zeroes M(nz) and fv_lite/fv_wright seed "// &
                              "p_edge(nz+1)=0.")
            has_error = .true.
         end if
         ! Inert-configuration warning, not a refusal (the house rule, cf.
         ! &ocean_psurf_nml in_eos and &ocean_tidal_mixing_nml e_uniform).
         ! `ms%p_top = metrics%p_ice_ref + sf%p_surf` has TWO producers,
         ! and the warning must name both or it is a lie: the
         ! &ocean_psurf_nml seam (the atmospheric half, `sf%p_surf`) and
         ! the &ocean_cavity_dyn_nml ice-shelf load (the static half,
         ! `p_ice_ref = rho_ref*g*z_draft`, assembled in
         ! `configure_ocean_cavity`).  Under a cavity `p_top` carries the
         ! ice load and `p_top_in_bc` is not merely live — it is REQUIRED
         ! for a varying draft, refused above and again at configure.  The
         ! warning fires only when NEITHER producer is on, which is the
         ! one case in which `p_top` really is the zero array it was
         ! allocated as.
         if (.not. p_top_has_producer(cfg)) then
            call logger%warning("&ocean_pgf_nml p_top_in_bc=.true. is INERT "// &
                                "without a producer for ms%p_top: enable "// &
                                "&ocean_psurf_nml (the atmospheric surface-pressure "// &
                                "seam) or &ocean_cavity_dyn_nml (the static "// &
                                "ice-shelf load), else p_top is the zero array and "// &
                                "pa(nz+1) is unchanged.")
         end if
      end if
      ! ---- Static ice-shelf cavity geometry (&ocean_cavity_dyn_nml, P5.1) ----
      ! The draft is absorbed into the barotropic DATUM (bt_H_ref =
      ! b - z_draft), so every consumer of the water-column thickness
      ! D = bt_H_ref + bt_eta is correct with no cavity branch of its own.
      ! What that buys is paid for by a narrow envelope, and EVERY
      ! restriction below fails loud naming the knob and the reason: a
      ! cavity that silently runs outside it looks plausible and is wrong
      ! (a coordinate anchored at z = 0 under 500 m of ice, a second
      ! un-reconciled surface load, a wide-halo BT clone with no draft).
      if (cfg%ocean%cavity_dyn%enable) then
         if (trim(cfg%sim_type) /= "ocean") then
            call logger%error("&ocean_cavity_dyn_nml enable=.true. requires "// &
                              "sim_type='ocean'")
            has_error = .true.
         end if
         ! --- geometry source envelope ---
         select case (trim(adjustl(cfg%ocean%cavity_dyn%draft_config)))
         case ("none", "flat", "linear")
            continue
         case ("file")
            ! Ships single-rank, through the PR-14 static-2-D reader
            ! (`ocean_data_input_load_static_2d`), which DOES apply the
            ! global offset — but the cavity as a whole is single-rank
            ! fenced below, and the loader re-asserts it.
            if (len_trim(cfg%ocean%cavity_dyn%draft_file) == 0) then
               call logger%error("&ocean_cavity_dyn_nml draft_config='file' requires "// &
                                 "draft_file")
               has_error = .true.
            end if
            if (len_trim(cfg%ocean%cavity_dyn%draft_var) == 0) then
               call logger%error("&ocean_cavity_dyn_nml draft_config='file' requires "// &
                                 "draft_var (the 2-D variable name; ISOMIP+ ships "// &
                                 "'iceDraft')")
               has_error = .true.
            end if
            if (parse_cavity_draft_sign(cfg%ocean%cavity_dyn%draft_sign) == &
                CAVITY_SIGN_INVALID) then
               call logger%error("&ocean_cavity_dyn_nml draft_sign='"// &
                                 trim(adjustl(cfg%ocean%cavity_dyn%draft_sign))// &
                                 "' is not recognised (depth|positive_down|"// &
                                 "elevation|positive_up).  There is no default that "// &
                                 "guesses from the data: the ISOMIP+ file carries an "// &
                                 "ELEVATION (z_d <= 0) and a depth file carries "// &
                                 "z_draft >= 0, and the two differ by the whole load.")
               has_error = .true.
            end if
         case default
            call logger%error("&ocean_cavity_dyn_nml draft_config='"// &
                              trim(adjustl(cfg%ocean%cavity_dyn%draft_config))// &
                              "' is not recognised (none|flat|linear|file)")
            has_error = .true.
         end select
         select case (trim(adjustl(cfg%ocean%cavity_dyn%draft_source)))
         case ("draft", "thickness")
            continue
         case ("in_situ")
            call logger%error("&ocean_cavity_dyn_nml draft_source='in_situ' (true "// &
                              "isostasy, p_ice = g*int(rho)) is not implemented: it "// &
                              "needs a per-column root find and does NOT admit exact "// &
                              "discrete rest in the split solver.  Use 'draft' "// &
                              "(the Boussinesq-isostatic flotation load ISOMIP+ "// &
                              "prescribes) or 'thickness'.")
            has_error = .true.
         case default
            call logger%error("&ocean_cavity_dyn_nml draft_source='"// &
                              trim(adjustl(cfg%ocean%cavity_dyn%draft_source))// &
                              "' is not recognised (draft|thickness|in_situ)")
            has_error = .true.
         end select
         ! `"linear"` measures its profile from `draft_x0`, so an
         ! unbounded (sentinel) anchor would make `draft_depth` meaningless
         ! and the whole shelf depth an artefact of 1e30*slope.
         if (trim(adjustl(cfg%ocean%cavity_dyn%draft_config)) == "linear" .and. &
             abs(cfg%ocean%cavity_dyn%draft_x0) >= 1.0e29_wp) then
            call logger%error("&ocean_cavity_dyn_nml draft_config='linear' requires "// &
                              "a finite draft_x0: it is the ANCHOR of the profile "// &
                              "(draft_depth is the draft AT draft_x0), not just the "// &
                              "western edge of the box.")
            has_error = .true.
         end if
         if (cfg%ocean%cavity_dyn%draft_depth < 0.0_wp) then
            call logger%error("&ocean_cavity_dyn_nml draft_depth must be >= 0 "// &
                              "(it is a DEPTH below z = 0, positive down)")
            has_error = .true.
         end if
         if (cfg%ocean%cavity_dyn%h_min_cavity <= 0.0_wp) then
            call logger%error("&ocean_cavity_dyn_nml h_min_cavity must be > 0 "// &
                              "(the grounding cutoff; 0 would admit a zero-thickness "// &
                              "water column under the ice)")
            has_error = .true.
         end if
         if (cfg%ocean%cavity_dyn%grounded_max_frac <= 0.0_wp .or. &
             cfg%ocean%cavity_dyn%grounded_max_frac > 1.0_wp) then
            call logger%error("&ocean_cavity_dyn_nml grounded_max_frac must be in "// &
                              "(0, 1] (the fraction of interior columns allowed to "// &
                              "ground)")
            has_error = .true.
         end if
         if (trim(adjustl(cfg%ocean%cavity_dyn%draft_source)) == "thickness" .and. &
             cfg%ocean%cavity_dyn%rho_ice <= 0.0_wp) then
            call logger%error("&ocean_cavity_dyn_nml draft_source='thickness' "// &
                              "requires rho_ice > 0")
            has_error = .true.
         end if
         ! --- the one atmospheric-forcing path the cover mask does NOT
         !     reach (P2c) ---
         ! Every static forcing field is masked: the wind pair and the
         ! scalar q_heat/q_salt once at configure, the component bands
         ! every thermo step in the assembler.  The FILE-DRIVEN override
         ! is the exception: `ocean_data_forcing_apply` rewrites `tau_x`/
         ! `tau_y` (and, with `heat_to_component=.false.`, `Q_heat`
         ! itself) from the next time bracket with no access to
         ! `metrics%cover_frac` — it is handed `ss`, `sf`, `grid` and
         ! `bc`, and nothing else.  Re-masking per bracket means
         ! threading the metrics slot through the reader, which is a
         ! separate change.  Refused rather than half-wired: a cavity run
         ! whose wind is silently restored to its unmasked file value on
         ! the first bracket read looks entirely plausible and is wrong.
         if (cfg%ocean%dataovr%enable) then
            call logger%error("&ocean_cavity_dyn_nml enable=.true. is mutually "// &
                              "exclusive with &ocean_dataovr_nml enable=.true.  The "// &
                              "ice-cover mask on the atmospheric forcing is applied "// &
                              "to the wind pair at configure and to the surface-flux "// &
                              "components in the assembler; the data-override reader "// &
                              "rewrites tau_x/tau_y (and Q_heat, unless "// &
                              "heat_to_component=.true.) per time bracket without "// &
                              "the cover, which would restore the unmasked "// &
                              "atmosphere under the shelf.  Follow-up: thread "// &
                              "cover_frac through ocean_data_forcing_apply and the "// &
                              "ocean_seam_refresh_surface_stress seam.")
            has_error = .true.
         end if
         ! --- pressure-gradient envelope ---
         if (trim(adjustl(cfg%ocean%pgf%form)) /= "fv_mom6") then
            call logger%error("&ocean_cavity_dyn_nml enable=.true. requires "// &
                              "&ocean_pgf_nml form='fv_mom6' (got '"// &
                              trim(adjustl(cfg%ocean%pgf%form))//"'). Only the "// &
                              "FV_MOM6 family builds the pa(nz+1) pressure-stack "// &
                              "boundary condition the ice load is injected into; "// &
                              "mont hard-zeroes M(nz) and fv_lite/fv_wright seed "// &
                              "p_edge(nz+1)=0.")
            has_error = .true.
         end if
         ! P5.2 — the load must have a consumer once it has a GRADIENT.
         ! `p_top_in_bc` is the only route by which the isostatic load
         ! rho_ref*g*z_draft reaches the FV_MOM6 pa(nz+1) surface BC;
         ! without it a varying draft leaves the pressure stack ~5e6 Pa
         ! off its anomaly scale, the unsplit driver feels a raw
         ! g*grad(z_draft), and `correction_h_weighted` turns the
         ! uncancelled depth-uniform force into a real per-layer shear.
         ! REFUSED rather than auto-enabled: an answer-changing knob that
         ! a second namelist group switches on behind the user's back is
         ! exactly the class of silent coupling this file exists to
         ! prevent.  A UNIFORM draft is exempt — a load with no gradient
         ! is bit-identically inert in the top BC (the theorem in
         ! `compute_fv_mom6_impl`'s docstring) — which is what keeps the
         ! flat-lid datum-equivalence gate expressible.  `draft_config`
         ! "none" is uniform (identically zero); "flat" is uniform only
         ! when no box bound clips it, else the calving front is a step.
         if (.not. cfg%ocean%pgf%p_top_in_bc) then
            block
               character(len=:), allocatable :: dcfg
               dcfg = trim(adjustl(cfg%ocean%cavity_dyn%draft_config))
               if (.not. cavity_draft_is_uniform(cfg)) then
                  call logger%error("&ocean_cavity_dyn_nml enable=.true. with "// &
                                    "draft_config='"//dcfg//"' requires "// &
                                    "&ocean_pgf_nml p_top_in_bc=.true.  That knob is "// &
                                    "the ONLY route by which the isostatic load "// &
                                    "rho_ref*g*z_draft reaches the FV_MOM6 pa(nz+1) "// &
                                    "surface boundary condition; without it a draft "// &
                                    "that VARIES leaves the pressure stack ~5e6 Pa "// &
                                    "off its anomaly scale and the column out of "// &
                                    "hydrostatic balance.  Only a draft that is "// &
                                    "uniform over the whole domain is exempt (a load "// &
                                    "with no gradient is provably inert there).")
                  has_error = .true.
               end if
            end block
         end if
         if (cfg%ocean%pgf%gfs_scale /= 1.0_wp) then
            call logger%error("&ocean_cavity_dyn_nml enable=.true. requires "// &
                              "&ocean_pgf_nml gfs_scale=1: the datum "// &
                              "(bt_H_ref, which the barotropic substep feels "// &
                              "through g_bt = gfs_scale*GRAVITY) and the load "// &
                              "(rho_ref*GRAVITY*z_draft, which the PGF feels "// &
                              "through GRAVITY) would then sit on two different "// &
                              "gravities and drift apart.")
            has_error = .true.
         end if
         ! --- vertical-coordinate envelope ---
         ! sigma (and zstar-lite, which shares its ocean branch) rescales
         ! the live column and so follows the draft for free; `z_fixed`
         ! is the FIRST family taught about the ice base explicitly
         ! (P6.2) — it reads `vcoord%z_top`, keeps its nominal interface
         ! depths GEOPOTENTIAL, vanishes the layers that outcrop into the
         ! ice to the inert filler and cuts the first live layer at the
         ! draft (Yung, Hallberg, Adcroft & Morrison 2026, JAMES 18,
         ! e2025MS005645, Fig. 1b).  Its own envelope is fenced below.
         !
         ! The refusal used to be a two-value whitelist whose message
         ! enumerated SIX families and named neither `lagrangian` nor
         ! `zstar_sigma` — both of which it refused.  A refusal that does
         ! not name what it refused, or gives a reason that is not the
         ! real one, sends the operator to fix the wrong thing.  So the
         ! accept test stays a whitelist (nothing else has been validated
         ! under a shelf) but the message now carries the offending
         ! family's OWN reason, and the reasons are not all "anchors at
         ! z = 0": three of the seven refused families are geometrically
         ! datum-safe and are refused for want of validation, which is a
         ! different, and recoverable, kind of no.  Measured per family in
         ! `test_ocean_vcoord_interface_depths`.
         block
            integer :: cav_vcoord_code
            real(wp) :: cav_h_nominal
            character(len=:), allocatable :: cav_reason
            cav_vcoord_code = parse_vcoord_type(cfg%vcoord_type, &
                                                default_code=VCOORD_EULERIAN_Z)
            if (.not. (cav_vcoord_code == VCOORD_SIGMA .or. &
                       cav_vcoord_code == VCOORD_ZSTAR .or. &
                       cav_vcoord_code == VCOORD_Z_FIXED)) then
               select case (cav_vcoord_code)
               case (VCOORD_LAGRANGIAN)
                  cav_reason = "'lagrangian' is geometrically datum-FREE (the target "// &
                               "IS the live h_layer and the remap is a no-op), so the "// &
                               "placement objection does not apply to it.  It is "// &
                               "refused in v1 only because no cavity run has been "// &
                               "validated on an isopycnal coordinate; it is the "// &
                               "intended isopycnal control leg of the coordinate study"
               case (VCOORD_EULERIAN_Z)
                  cav_reason = "'eulerian_z' is a stretched SIGMA with the free "// &
                               "surface DROPPED (it is not a geopotential coordinate "// &
                               "despite the name), and the dropped eta — not a z "// &
                               "anchor — is why it cannot carry an ice-base "// &
                               "displacement"
               case (VCOORD_ZSIGMA)
                  cav_reason = "'zsigma' is refused on EVERY path, cavity or not: its "// &
                               "deep branch reads z_ref_global as metres while the "// &
                               "only writer fills it with the dimensionless k/nz"
               case (VCOORD_ZSTAR_SIGMA)
                  cav_reason = "'zstar_sigma' is a purely FRACTIONAL rescale of the "// &
                               "global reference table, hence datum-invariant and "// &
                               "geometrically safe under a draft.  It is refused in "// &
                               "v1 only because it follows BOTH boundaries (so it "// &
                               "solves nothing a cavity needs) and no cavity run has "// &
                               "been validated on it"
               case (VCOORD_ZSTAR_FULL)
                  cav_reason = "'zstar_full' builds its per-column reference table "// &
                               "from the TRUE bed while the target walk sees only the "// &
                               "live thickness, so under a draft it keeps the table's "// &
                               "SHALLOW entries: the fine near-surface band lands "// &
                               "against the ice base and the deep water is carried by "// &
                               "shallow-ocean spacing.  It does not merely anchor at "// &
                               "z = 0, it inverts which half of the column is resolved"
               case default
                  cav_reason = "'"//trim(cfg%vcoord_type)//"' is a DENSITY-space "// &
                               "coordinate with no geometric anchor, so the placement "// &
                               "objection does not apply.  It is refused in v1 because "// &
                               "cavity x isopycnal is unvalidated; note also that "// &
                               "hycom's z* nominal-floor band accumulates from the "// &
                               "column top, so under a shelf that band is "// &
                               "draft-FOLLOWING and must not be quoted as z-like"
               end select
               call logger%error("&ocean_cavity_dyn_nml enable=.true. accepts "// &
                                 "vcoord_type='sigma', 'zstar' (zstar-lite) or "// &
                                 "'z_fixed' ONLY — the two families that rescale the "// &
                                 "live column and so follow the ice base for free, "// &
                                 "plus the one that has been TAUGHT the ice base "// &
                                 "(z_fixed reads vcoord%z_top, keeps its nominal "// &
                                 "interface depths geopotential, and vanishes the "// &
                                 "layers that outcrop into the ice).  Got '"// &
                                 trim(cfg%vcoord_type)//"': "//cav_reason//".")
               has_error = .true.
            end if

            if (cav_vcoord_code == VCOORD_Z_FIXED) then
               ! ---- The staircase, and why this is a WARNING ----
               !
               ! A z-like coordinate removes the interior sigma tilt but
               ! replaces it with the ice-base STAIRCASE: where the draft
               ! crosses a nominal level the two columns' filler counts
               ! differ by one, the interface offset across that face
               ! jumps by up to `h_nominal`, and the FV_MOM6 acceleration
               ! in a vanished layer is h-INDEPENDENT.  The corrections
               ! that arrest it — top-side mass weighting (MWIPG) and the
               ! interior reference interface, Yung, Hallberg, Adcroft &
               ! Morrison (2026), JAMES 18, e2025MS005645, §3.3.2 and
               ! §3.2 — are NOT implemented.  A UNIFORM draft has no
               ! staircase (every column vanishes the same layers and
               ! cuts at the same depth, so every interface offset is
               ! identically zero) and is validated bit-zero to 30 days;
               ! a VARYING draft is not validated at all, and measured on
               ! `cavity_sloping_lid_rest_zfixed.nml` it is 47x the sigma
               ! leg's resting pressure-gradient residual at step 1 and
               ! ends in a non-finite state on day 18 (gfortran) or a
               ! saturated En = 3.8E-04 (nvfortran GPU).  That is the
               ! knob-OFF state.  `&vcoord_nml zfixed_closed_faces`
               ! (Adcroft, Hill & Marshall 1997 partial steps) closes the
               ! staircase faces, and with it the same file completes 30
               ! days and ISOMIP+ Ocean0 runs 30 days — but the residual on
               ! the faces left open is still ~2 decades above the sigma
               ! leg, so the warning stays, reworded for each state.
               !
               ! WARNING, not a refusal, deliberately: those corrections
               ! are the next slices and they need this configuration to
               ! be runnable to be developed and measured against.  What
               ! the user must not do is walk into it silently.
               if (.not. cavity_draft_is_uniform(cfg)) then
                  if (cfg%zfixed_closed_faces) then
                     call logger%warning("&vcoord_nml vcoord_type='z_fixed' under "// &
                                         "&ocean_cavity_dyn_nml with a draft that VARIES "// &
                                         "(draft_config='"// &
                                         trim(adjustl(cfg%ocean%cavity_dyn%draft_config))// &
                                         "') is EXPERIMENTAL.  &vcoord_nml "// &
                                         "zfixed_closed_faces=.true. closes the ice-base "// &
                                         "STAIRCASE faces, which is what makes this "// &
                                         "configuration runnable: the sloping-lid rest "// &
                                         "case completes 30 days and ISOMIP+ Ocean0 runs "// &
                                         "30 days.  It does not remove the staircase "// &
                                         "residual on the faces that stay open — the "// &
                                         "top-side mass weighting (MWIPG) and the interior "// &
                                         "reference interface, Yung, Hallberg, Adcroft & "// &
                                         "Morrison (2026), JAMES 18, e2025MS005645, "// &
                                         "sections 3.3.2 and 3.2, are not implemented — so "// &
                                         "the resting sloping-lid case still carries about "// &
                                         "two decades more spurious energy than its sigma "// &
                                         "leg.  Only a UNIFORM draft is bit-zero at rest.")
                  else
                     call logger%warning("&vcoord_nml vcoord_type='z_fixed' under "// &
                                         "&ocean_cavity_dyn_nml with a draft that VARIES "// &
                                         "(draft_config='"// &
                                         trim(adjustl(cfg%ocean%cavity_dyn%draft_config))// &
                                         "') and &vcoord_nml zfixed_closed_faces=.false. "// &
                                         "is NOT VALIDATED.  Only a UNIFORM draft is: "// &
                                         "there every column vanishes the same layers and "// &
                                         "cuts at the same depth, so the answer is "// &
                                         "bit-zero at rest.  Where the draft crosses a "// &
                                         "nominal level the filler count changes column "// &
                                         "to column and the FV pressure gradient across "// &
                                         "the open ice-base STAIRCASE drives a spurious "// &
                                         "flow: measured at rest, 47x the sigma leg's "// &
                                         "step-1 residual, and the run does not survive "// &
                                         "30 days.  Set zfixed_closed_faces=.true. (the "// &
                                         "partial-step face closure, with which a varying "// &
                                         "draft does survive 30 days), or use "// &
                                         "vcoord_type='sigma' for a sloping lid.")
                  end if
               end if
               ! ---- z_fixed × cavity, v1 envelope ----
               !
               ! Under a quasi-geopotential coordinate the ice base cuts
               ! the nominal stack, so on an ICE-COVERED column `k = nz`
               ! is an inert filler (`h <= H_VANISHED`), NOT the
               ! ice-adjacent live layer.  This is a state no consumer in
               ! the tree has ever seen: no family on this branch
               ! vanishes a layer against the TOP.  The shared
               ! `k_top(i,j)` (the first live layer, counting down) that
               ! fixes them is the NEXT slice — P6.3 for the tracer/flux
               ! consumers, P6.4 for momentum and the boundary-layer
               ! schemes — so until it lands the only thing this
               ! combination may run is ADIABATIC DYNAMICS.
               !
               ! What is NOT refused, and why: OPEN-OCEAN columns have
               ! `z_top = 0`, so `k = nz` is live there and behaves
               ! exactly as today; and on a COVERED column the binary
               ! cover mask (`cover_frac`) already zeroes every
               ! ATMOSPHERIC input at source — wind stress (and with it
               ! `stress_mag`, hence u*), the scalar and assembled
               ! surface heat/salt fluxes, the shortwave deposit and both
               ! restoring increments.  Those therefore compose with a
               ! vanished `k = nz` and are left alone.  What follows is
               ! everything that acts ON the ice-adjacent layer itself,
               ! and so is not masked by anything.
               ! NOTE `&ocean_cavity_melt_nml enable` and
               ! `&ocean_tdrag_nml enable` used to be refused HERE, and
               ! are not any more: both are routed through the shared
               ! first-live-layer index `ms%k_top` (and its two face
               ! twins) — P6.3/P6.4.  The melt heat/salt deposit and the
               ! `freshwater="mass"` column source land on
               ! `k_top(i,j)`; the ice-ocean drag's band walk starts at
               ! `k_top_u/v`, gates on `H_VANISHED` instead of on zero,
               ! and captures its implicit-fold rate on the same row the
               ! vdiff diagonal adds it to.  `k_top ≡ nz` off a rigid
               ! top, so nothing else moved.
               if (cfg%ocean%vmix%use_kpp) then
                  call logger%error("&ocean_vmix_nml use_kpp=.true. is refused with "// &
                                    "vcoord_type='z_fixed' under a cavity: KPP's "// &
                                    "surface reference column is k = nz "// &
                                    "(b_ref = -g*rho_layer(nz)/rho_0, "// &
                                    "d_centre_ref = h(nz)/2), which on a covered "// &
                                    "column is the filler — a spurious buoyancy jump "// &
                                    "at the very first interface and a reference "// &
                                    "depth of ~1e-4 m.  Set use_kpp=.false. (an "// &
                                    "adiabatic cavity run is the v1 envelope).  "// &
                                    "Needs the shared k_top: follow-up P6.4.")
                  has_error = .true.
               end if
               if (cfg%ocean%epbl%enable) then
                  call logger%error("&ocean_epbl_nml enable=.true. is refused with "// &
                                    "vcoord_type='z_fixed' under a cavity: the EPBL "// &
                                    "column captures the surface at k = nz and gates "// &
                                    "on `h > 0` rather than on H_VANISHED, so a "// &
                                    "filler top layer reaches the specific-volume "// &
                                    "derivatives as an absurd concentration.  Needs "// &
                                    "the shared k_top: follow-up P6.4.")
                  has_error = .true.
               end if
               if (cfg%ocean%tracers%enable_ideal_age) then
                  call logger%error("&ocean_tracers_nml enable_ideal_age=.true. is "// &
                                    "refused with vcoord_type='z_fixed' under a "// &
                                    "cavity: the young-band reset writes "// &
                                    "hTr_age(:,:,nz) = young*h(:,:,nz), which on a "// &
                                    "covered column ventilates a filler (and there "// &
                                    "is no ventilation under a shelf at all).  "// &
                                    "Follow-up P6.3.")
                  has_error = .true.
               end if
               if (cfg%ocean%gm%enable) then
                  call logger%error("&ocean_gm_nml enable=.true. is refused with "// &
                                    "vcoord_type='z_fixed' under a cavity: the "// &
                                    "non-divergence closure dumps the residual "// &
                                    "streamfunction transport into k = nz.  GM x "// &
                                    "cavity is unvalidated on any coordinate; "// &
                                    "revisited with the coordinate study.")
                  has_error = .true.
               end if
               if (cfg%ocean%redi%enable .or. cfg%ocean%slopes%enable) then
                  call logger%error("&ocean_redi_nml / &ocean_slopes_nml enable="// &
                                    ".true. is refused with vcoord_type='z_fixed' "// &
                                    "under a cavity: the isopycnal-slope surface "// &
                                    "fill is built from h(:,:,nz), the filler on a "// &
                                    "covered column.  Unvalidated with a cavity on "// &
                                    "any coordinate; revisited with the coordinate "// &
                                    "study.")
                  has_error = .true.
               end if
               ! ---- kappa_h: refused ONLY with the face mask off ----
               !
               ! The objection is unchanged where it still applies: the
               ! along-coordinate tracer diffusion gates its T = hTr/h
               ! division on `h > 0` (1/0 armour) and not on H_VANISHED,
               ! weights the face flux by the ARITHMETIC mean thickness
               ! — a filler beside a 20 m cell is weighted by 10 m — and
               ! carries no mass-availability limiter, so it can drive a
               ! vanished cell's hTr strongly negative in one step.
               !
               ! But `&vcoord_nml zfixed_closed_faces` already answers
               ! it, from the other end: `ocean_hdiff_tracer_step`
               ! multiplies every face flux by `metrics%open_u/open_v`,
               ! which is exactly zero wherever the layer is a filler on
               ! EITHER side.  A filler cell then has all four of its
               ! own-layer faces closed, so its hTr divergence is
               ! identically zero and the garbage concentration the
               ! `h > 0` gate computes for it never leaves the cell.
               ! Flux-zero, not flux-limited, so the scheme stays
               ! conservative by construction.  That is the P6.5 gate,
               ! reached through the mask instead of through a new
               ! threshold — and it is a strictly stronger statement,
               ! because it also closes the partial⇄filler face the
               ! threshold alone would leave open on the thick side.
               if (cfg%ocean%hdiff%kappa_h /= 0.0_wp .and. &
                   .not. cfg%zfixed_closed_faces) then
                  call logger%error("&ocean_hdiff_nml kappa_h /= 0 is refused with "// &
                                    "vcoord_type='z_fixed' under a cavity UNLESS "// &
                                    "&vcoord_nml zfixed_closed_faces=.true.: the "// &
                                    "along-coordinate tracer diffusion gates its "// &
                                    "T = hTr/h division on `h > 0` (1/0 armour), "// &
                                    "not on H_VANISHED, weights the face flux by the "// &
                                    "ARITHMETIC mean thickness — so a filler beside "// &
                                    "a 20 m cell is weighted by 10 m — and has no "// &
                                    "mass-availability limiter, so it can drive a "// &
                                    "vanished cell's hTr strongly negative in one "// &
                                    "step.  With the partial-step face mask on, "// &
                                    "every face touching a filler is CLOSED and the "// &
                                    "flux is exactly zero, which answers all three.  "// &
                                    "Set zfixed_closed_faces=.true., or kappa_h=0.")
                  has_error = .true.
               end if
               if (cfg%ocean%kshear%enable) then
                  call logger%error("&ocean_kappa_shear_nml enable=.true. is refused "// &
                                    "with vcoord_type='z_fixed' under a cavity: the "// &
                                    "JHL08 column solve closes its SURFACE row on "// &
                                    "k = nz (u_c/v_c/t_c/s_c(nz)), which on a "// &
                                    "covered column is an inert filler inside the "// &
                                    "ice draft, and its own massless-merge helper is "// &
                                    "off by default.  Not covered by the shared "// &
                                    "k_top (P6.3/P6.4), which routes the FORCING "// &
                                    "sites; a column solver wants the compacted "// &
                                    "column rdb_massless already builds.")
                  has_error = .true.
               end if
               if (cfg%ocean%tidal_mixing%enable) then
                  call logger%error("&ocean_tidal_mixing_nml enable=.true. is "// &
                                    "refused with vcoord_type='z_fixed' under a "// &
                                    "cavity: the N^2 column sets its top boundary "// &
                                    "at k = nz and forms the k = nz-1 interface "// &
                                    "spacing as 0.5*(h(nz-1) + h(nz)), which on a "// &
                                    "covered column HALVES that spacing against a "// &
                                    "filler and inflates the buoyancy frequency at "// &
                                    "the first live interface.  Unvalidated under a "// &
                                    "rigid top.")
                  has_error = .true.
               end if
               if (cfg%regrid_time_scale > 0.0_wp) then
                  call logger%error("&vcoord_nml regrid_time_scale > 0 is refused "// &
                                    "with vcoord_type='z_fixed' under a cavity: the "// &
                                    "grid time-filter is a per-layer convex blend "// &
                                    "that does not respect the vanish marker, so a "// &
                                    "filler relaxing toward h_min from a live "// &
                                    "thickness passes THROUGH H_VANISHED and the "// &
                                    "layer oscillates live/dead on successive steps "// &
                                    "— the remap drain deleting its content on the "// &
                                    "step it reads dead.  Excluding the fillers from "// &
                                    "the blend is its own change.")
                  has_error = .true.
               end if
               ! ISOMIP+ (Asay-Davis et al. 2016 §3.1.5): "the minimum
               ! thickness is likely to be approximately two grid cells
               ! (~40 m if z levels are equally spaced)".  Under a
               ! terrain-following coordinate a just-afloat column still
               ! carries all nz layers; under a z-like one a column
               ! thinner than 2*h_nominal carries ONE partial live layer
               ! and nz-1 fillers, and every column-walking consumer then
               ! operates on a single cell.  `h_nominal` is
               ! `max_depth/nz` — Z_FIXED's spacing is written from
               ! `&ocean_topo_nml max_depth`, deliberately with no second
               ! spelling.
               cav_h_nominal = 0.0_wp
               if (cfg%nz_layers > 0) then
                  cav_h_nominal = cfg%ocean%topo%max_depth/real(cfg%nz_layers, wp)
               end if
               if (cav_h_nominal > 0.0_wp .and. &
                   cfg%ocean%cavity_dyn%h_min_cavity < 2.0_wp*cav_h_nominal) then
                  call logger%error("&ocean_cavity_dyn_nml h_min_cavity = "// &
                                    to_string(cfg%ocean%cavity_dyn%h_min_cavity)// &
                                    " m is below 2*h_nominal = "// &
                                    to_string(2.0_wp*cav_h_nominal)//" m "// &
                                    "(h_nominal = &ocean_topo_nml max_depth / "// &
                                    "nz_layers = "//to_string(cav_h_nominal)// &
                                    " m), which vcoord_type='z_fixed' requires under "// &
                                    "a cavity: a thinner water column carries a "// &
                                    "single partial live layer and nz-1 inert "// &
                                    "fillers.  This is ISOMIP+'s own rule "// &
                                    "(Asay-Davis et al. 2016, GMD 9, 2471, "// &
                                    "§3.1.5).  NOTE it moves the grounding line "// &
                                    "relative to a sigma leg — a confound the "// &
                                    "coordinate study must control for.")
                  has_error = .true.
               end if
            end if
         end block
         if (trim(cfg%thickness_config) == "uniform_z") then
            call logger%error("&ocean_cavity_dyn_nml enable=.true. is mutually "// &
                              "exclusive with thickness_config='uniform_z': that "// &
                              "seed lays uniform z interfaces from z = 0 down, so "// &
                              "under a draft it would fill the ice with water.  Use "// &
                              "the default sigma-style split.")
            has_error = .true.
         end if
         ! NOTE: `&ocean_zinit_nml enable` used to be refused here — the
         ! overlay measured depth from the COLUMN TOP, which under a
         ! draft is `z_draft` metres below `z = 0`, so a geopotential
         ! profile landed systematically too shallow.  `build_z_ctr` now
         ! takes the column-top depth and the seed passes
         ! `metrics%z_draft`, so the two compose and the refusal is gone.
         ! --- solver envelope ---
         if (cfg%ocean%bt%n_inner < 1 .and. .not. cfg%ocean%bt%auto_n_inner) then
            call logger%error("&ocean_cavity_dyn_nml enable=.true. requires the "// &
                              "SPLIT solver (&ocean_bt_nml n_inner >= 1, or "// &
                              "auto_n_inner): the unsplit driver has neither the "// &
                              "barotropic correction nor the eta_forcing seam, so it "// &
                              "carries no barotropic response to a surface load and "// &
                              "the cavity is unvalidated there.")
            has_error = .true.
         end if
         if (cfg%ocean%bt%bt_halo > 0) then
            call logger%error("&ocean_cavity_dyn_nml enable=.true. is mutually "// &
                              "exclusive with &ocean_bt_nml bt_halo > 0: the "// &
                              "wide-halo BT clone rebuilds its own metrics from the "// &
                              "grid formula and carries no z_draft, so its reference "// &
                              "depth would be the bed and its solve would ignore the "// &
                              "ice.")
            has_error = .true.
         end if
         if (cfg%px*cfg%py > 1) then
            call logger%error("&ocean_cavity_dyn_nml enable=.true. is single-rank in "// &
                              "v1 (px*py = "//to_string(cfg%px*cfg%py)//").  The "// &
                              "draft halo itself is two lines, but the grounding "// &
                              "statistics are global reductions the v1 configure does "// &
                              "not take, and draft_config='file' inherits the "// &
                              "local-nx reader.  Lifting the fence is its own PR.")
            has_error = .true.
         end if
         ! --- mutually exclusive capabilities ---
         if (cfg%ocean%wetdry%enable) then
            call logger%error("&ocean_cavity_dyn_nml enable=.true. is mutually "// &
                              "exclusive with &ocean_wetdry_nml enable: both decide "// &
                              "whether a column can carry water, from different "// &
                              "thresholds, and wet/dry also forces split_scheme="// &
                              "'ssp_rk2'.")
            has_error = .true.
         end if
         if (cfg%ocean%porous%enable) then
            call logger%error("&ocean_cavity_dyn_nml enable=.true. is mutually "// &
                              "exclusive with &ocean_porous_nml enable: both narrow "// &
                              "the same faces from static geometry and the "// &
                              "combination is unvalidated.")
            has_error = .true.
         end if
         if (cfg%ocean%ice%enable) then
            call logger%error("&ocean_cavity_dyn_nml enable=.true. is mutually "// &
                              "exclusive with &ocean_ice_nml enable: sea ice is a "// &
                              "SECOND surface load on the same column, and the two "// &
                              "are not reconciled in v1.")
            has_error = .true.
         end if
         if (cfg%ocean%tides%enable .and. cfg%ocean%tides%use_sal) then
            call logger%error("&ocean_cavity_dyn_nml enable=.true. is mutually "// &
                              "exclusive with &ocean_tides_nml use_sal: the scalar "// &
                              "SAL elevation is beta_sal*bt_eta, and under the cavity "// &
                              "datum bt_eta is the departure from the LOADED "// &
                              "equilibrium, not the sea-surface elevation the SAL "// &
                              "response is defined on.  The body tide itself "// &
                              "(enable alone) is datum-independent and allowed.")
            has_error = .true.
         end if
         ! --- honest-inertness warnings (the house rule: warn, do not refuse) ---
         if (trim(adjustl(cfg%ocean%cavity_dyn%draft_config)) == "none") then
            call logger%warning("&ocean_cavity_dyn_nml enable=.true. with "// &
                                "draft_config='none': z_draft is identically zero, so "// &
                                "bt_H_ref = b and the run is bit-identical to a "// &
                                "cavity-free one.")
         end if
         ! Design Q4: WARN, do not refuse, on cavity x in_eos=.false.  The
         ! load is depth-uniform in `pa`, so the rest and equivalence
         ! gates pass either way; what is wrong without `in_eos` is the
         ! THERMOBARICITY, and only a pressure-dependent EOS has any.  A
         ! linear EOS is pressure-blind, so there the knob is honestly
         ! inert and there is nothing to warn about.
         if (.not. cfg%ocean%psurf%in_eos .and. &
             trim(adjustl(cfg%ocean%eos%eos)) /= "linear") then
            call logger%warning("&ocean_cavity_dyn_nml enable=.true. with a "// &
                                "NONLINEAR equation of state (&ocean_eos_nml eos='"// &
                                trim(adjustl(cfg%ocean%eos%eos))//"') but without "// &
                                "&ocean_psurf_nml in_eos=.true.: the in-situ EOS "// &
                                "pressure still starts at 0 Pa at the ice base, so "// &
                                "the EOS ignores up to ~5e6 Pa of ice load (wrong "// &
                                "thermobaricity, wrong freezing point).  Harmless "// &
                                "for the rest and equivalence gates; not for a "// &
                                "production cavity.")
         end if
      end if
      ! ---- Ice-shelf basal melt, v1 scope (Phase 2b) ----
      ! Every restriction fails loud.  One of them exists because a piece
      ! of the coupling is NOT in this slice: there is no per-cell cover
      ! mask on the atmospheric forcing yet, and leaving that silently
      ! half-wired would run an atmosphere through several hundred metres
      ! of solid ice.  The interface pressure is NOT such a gap — the
      ! cavity assembles `ms%p_top = p_ice_ref + sf%p_surf` (P5.2) and
      ! the melt liquidus is its third consumer.
      if (cfg%ocean%cavity_melt%enable) then
         if (.not. cfg%ocean%cavity_dyn%enable) then
            call logger%error("&ocean_cavity_melt_nml enable=.true. requires "// &
                              "&ocean_cavity_dyn_nml enable=.true.  The melt "// &
                              "interface stands on that group's geometry: without a "// &
                              "draft there is no ice base, cover_frac is identically "// &
                              "zero and the kernel would never be called.")
            has_error = .true.
         end if
         if (trim(adjustl(cfg%ocean%eos%tfreeze_set)) /= "isomip") then
            call logger%error("&ocean_cavity_melt_nml enable=.true. requires "// &
                              "&ocean_eos_nml tfreeze_set='isomip' (got '"// &
                              trim(adjustl(cfg%ocean%eos%tfreeze_set))//"').  The "// &
                              "two shipped liquidi differ by ~0.03 degC at S = 34.5, "// &
                              "which is enough to flip the SIGN of the basal melt "// &
                              "rate over a 0.03 degC band of far-field temperature — "// &
                              "the sea-ice set is not a defensible default for a "// &
                              "cavity.")
            has_error = .true.
         end if
         if (.not. cfg%ocean%forcing%enable_components) then
            call logger%error("&ocean_cavity_melt_nml enable=.true. requires "// &
                              "&ocean_forcing_nml enable_components=.true.  The melt "// &
                              "fluxes are delivered as the OWNED components "// &
                              "heat_cavity/salt_cavity, which are allocated only "// &
                              "with the component set, and only "// &
                              "ocean_surface_flux_assemble folds them into "// &
                              "Q_heat/Q_salt.  Writing Q_heat/Q_salt directly is "// &
                              "refused by the fill contract — the sea-ice coupler "// &
                              "full-overwrites those.")
            has_error = .true.
         end if
         block
            integer :: melt_law_code, melt_ice_code
            melt_law_code = parse_cavity_exchange_law(cfg%ocean%cavity_melt%exchange_law)
            select case (melt_law_code)
            case (CAVITY_LAW_CONST_GAMMA, CAVITY_LAW_HJ99, CAVITY_LAW_YUNG25)
               continue
            case (CAVITY_LAW_INVALID)
               call logger%error("&ocean_cavity_melt_nml exchange_law = '"// &
                                 trim(adjustl(cfg%ocean%cavity_melt%exchange_law))// &
                                 "' is not recognised (shipped: const_gamma, hj99, "// &
                                 "yung25).")
               has_error = .true.
            case default
               call logger%error("&ocean_cavity_melt_nml exchange_law = '"// &
                                 trim(adjustl(cfg%ocean%cavity_melt%exchange_law))// &
                                 "' is RESERVED, not implemented "// &
                                 "(CAVITY_MELT_NOT_IMPLEMENTED).  Its enum value is "// &
                                 "nailed down so adding it later is not a "// &
                                 "renumbering, and the Python prototype has it, but "// &
                                 "no Fortran physics ships.  Shipped: const_gamma, "// &
                                 "hj99, yung25.")
               has_error = .true.
            end select
            ! Holland & Jenkins (1999) eq. (15) p. 1792 takes
            ! `ln(u* xi_N eta*^2 / (|f| h_nu))` and eq. (18) divides by
            ! `f L_O`: the law does not exist on the equator.  The
            ! per-column check over the covered cells is taken at
            ! configure (`configure_ocean_cavity_melt`), where f_centre
            ! exists; here we can only catch the whole-domain f = 0 case,
            ! and catching it early is worth the duplication.
            if (melt_law_code == CAVITY_LAW_HJ99) then
               if (cfg%coriolis_f == 0.0_wp .and. &
                   cfg%ocean%topo%coriolis_beta == 0.0_wp) then
                  call logger%error("&ocean_cavity_melt_nml exchange_law='hj99' on an "// &
                                    "f = 0 grid (coriolis_f = 0, coriolis_beta = 0). "// &
                                    "Holland & Jenkins (1999) eq. (15) takes "// &
                                    "ln(.../|f| h_nu) and eq. (18) divides by f*L_O, "// &
                                    "so the law has no value there — the same "// &
                                    "fail-loud stance &ocean_vmix_nml bkgnd_henyey "// &
                                    "takes on a cartesian grid.  Use "// &
                                    "exchange_law='const_gamma' or set a Coriolis "// &
                                    "parameter.")
                  has_error = .true.
               end if
            end if
            melt_ice_code = parse_cavity_ice_mode(cfg%ocean%cavity_melt%ice_conduction)
            select case (melt_ice_code)
            case (CAVITY_ICE_INSULATING, CAVITY_ICE_ADV_DIFF)
               continue
            case (CAVITY_ICE_INVALID)
               call logger%error("&ocean_cavity_melt_nml ice_conduction = '"// &
                                 trim(adjustl(cfg%ocean%cavity_melt%ice_conduction))// &
                                 "' is not recognised (shipped: insulating, "// &
                                 "adv_diff).")
               has_error = .true.
            case default
               call logger%error("&ocean_cavity_melt_nml ice_conduction = '"// &
                                 trim(adjustl(cfg%ocean%cavity_melt%ice_conduction))// &
                                 "' is RESERVED, not implemented "// &
                                 "(CAVITY_MELT_NOT_IMPLEMENTED).  The steady "// &
                                 "diffusive form makes q_ice independent of m_mass, "// &
                                 "so sign(m) = sign(T*) no longer holds and the "// &
                                 "pre-solve melt/freeze branch has to be revisited — "// &
                                 "it is not a coefficient change.")
               has_error = .true.
            end select
         end block
         if (cfg%ocean%cavity_melt%gamma_t <= 0.0_wp) then
            call logger%error("&ocean_cavity_melt_nml gamma_t must be > 0 (it "// &
                              "multiplies u* to give the heat exchange velocity; "// &
                              "zero is identically zero melt, which is what "// &
                              "enable=.false. is for).")
            has_error = .true.
         end if
         if (cfg%ocean%cavity_melt%gamma_s >= 0.0_wp .and. &
             cfg%ocean%cavity_melt%gamma_s <= 0.0_wp) then
            call logger%error("&ocean_cavity_melt_nml gamma_s = 0 is refused: the "// &
                              "three-equation form divides by gamma_s.  Leave it "// &
                              "NEGATIVE (the default) to take the ISOMIP+ "// &
                              "gamma_t/35, or set a positive value.")
            has_error = .true.
         end if
         if (cfg%ocean%cavity_melt%far_field_depth <= 0.0_wp) then
            call logger%error("&ocean_cavity_melt_nml far_field_depth must be > 0 m "// &
                              "(it is the thickness the far-field T/S/u are averaged "// &
                              "over; zero would sample nothing).")
            has_error = .true.
         end if
         ! --- Phase 3: real freshwater MASS, and its sea-level partner ---
         ! `freshwater="mass"` moves the meltwater as a REAL Boussinesq
         ! volume on the top layer.  Two envelope holes are refused BY
         ! NAME rather than half-wired, because each would put the mass
         ! source and the machinery that owns `h_layer(:,:,nz)` out of
         ! step with one another:
         !
         !   * dynamic wet/dry re-decides every step which columns carry
         !     water, and its positive-definite outflow limiter is the
         !     other writer of a top-layer thickness source.  Composing
         !     the two needs the limiter to SEE the melt volume;
         !   * a windowed tracer-advection ratio > 1 freezes `hTr` for
         !     `ratio` steps while `h` keeps moving, so the dilution the
         !     mass form relies on would be applied to a tracer load that
         !     is deliberately stale.
         block
            integer :: fw_code, vc_code
            fw_code = parse_cavity_freshwater(cfg%ocean%cavity_melt%freshwater)
            vc_code = parse_cavity_volume_comp(cfg%ocean%cavity_melt%volume_compensation)
            if (fw_code == CAVITY_FW_INVALID) then
               call logger%error("&ocean_cavity_melt_nml freshwater = '"// &
                                 trim(adjustl(cfg%ocean%cavity_melt%freshwater))// &
                                 "' is not recognised (shipped: virtual, mass).")
               has_error = .true.
            end if
            if (vc_code == CAVITY_VC_INVALID) then
               call logger%error("&ocean_cavity_melt_nml volume_compensation = '"// &
                                 trim(adjustl(cfg%ocean%cavity_melt%volume_compensation))// &
                                 "' is not recognised (shipped: none, "// &
                                 "uniform_open_ocean).")
               has_error = .true.
            end if
            if (fw_code == CAVITY_FW_MASS) then
               if (cfg%ocean%wetdry%enable) then
                  call logger%error("&ocean_cavity_melt_nml freshwater='mass' is "// &
                                    "refused with &ocean_wetdry_nml enable=.true.  "// &
                                    "Wet/dry owns the other top-layer thickness "// &
                                    "source (its positive-definite outflow limiter) "// &
                                    "and re-decides per step which columns hold "// &
                                    "water; composing the two needs the limiter to "// &
                                    "see the melt volume.  Follow-up: "// &
                                    "'cavity real freshwater under wet/dry'.")
                  has_error = .true.
               end if
               if (cfg%ocean%vmix%dt_tracer_advect_ratio > 1) then
                  call logger%error("&ocean_cavity_melt_nml freshwater='mass' is "// &
                                    "refused with &ocean_vmix_nml "// &
                                    "dt_tracer_advect_ratio > 1.  The windowed drain "// &
                                    "holds hTr fixed for the window while h keeps "// &
                                    "moving, so the dilution the mass form relies on "// &
                                    "would act on a deliberately stale tracer load.  "// &
                                    "Follow-up: 'cavity real freshwater under the "// &
                                    "windowed tracer-advect drain'.")
                  has_error = .true.
               end if
            end if
            if (vc_code == CAVITY_VC_UNIFORM_OPEN .and. fw_code /= CAVITY_FW_MASS) then
               call logger%error("&ocean_cavity_melt_nml volume_compensation="// &
                                 "'uniform_open_ocean' requires freshwater='mass'.  "// &
                                 "The virtual form adds no volume, so there is "// &
                                 "nothing to compensate and the sink would be a "// &
                                 "pure, unexplained mass loss.")
               has_error = .true.
            end if
         end block
         ! --- the cover mask SHIPS (P2c) ---
         ! Wind stress, surface restoring, shortwave penetration and the
         ! uniform scalar q_heat/q_salt were each refused here while
         ! there was no per-cell ice-COVER MASK on the atmospheric
         ! forcing.  There is one now:
         !
         !   * the wind-stress PAIR is masked face-wise at configure
         !     (`ocean_surface_stress_apply_cover`, called from
         !     `configure_ocean_cavity`), which also silences the
         !     implicit vdiff stress fold and the MLE front sampler —
         !     both read `ss%tau_x` raw — and refreshes `stress_mag`, so
         !     KPP/EPBL `u*` is zero-wind under cover;
         !   * `q_heat`/`q_salt`, the radiative/turbulent bands, the
         !     mass-flux enthalpies and `salt_flux` are masked in
         !     `ocean_surface_flux_assemble` (cover-aware twin), which is
         !     the one place they are still separable from the cavity's
         !     own `heat_cavity`/`salt_cavity` — and which makes the
         !     `Q_heat`/`Q_salt` KPP/EPBL read for `B_0` the MASKED
         !     values;
         !   * shortwave penetration and surface restoring carry the
         !     factor themselves (they do not route through `Q_*`).
         !
         ! `&ocean_psurf_nml enable` was never refused here: the cavity
         ! is the `ms%p_top` producer (P5.2: `p_ice_ref + sf%p_surf`), so
         ! the seam composes with the ice load rather than clobbering it.
         if (cfg%ocean%ice%enable) then
            call logger%error("&ocean_cavity_melt_nml enable=.true. is mutually "// &
                              "exclusive with &ocean_ice_nml enable=.true.  Sea ice "// &
                              "full-overwrites heat_added/salt_flux and runs its own "// &
                              "instant-relaxation basal flux with a DIFFERENT "// &
                              "liquidus set; two interface thermodynamics in one "// &
                              "column is not a configuration.")
            has_error = .true.
         end if
         ! (E4) Constant α/β under a NONLINEAR EOS in a cavity — a WARNING,
         ! deliberately, not a refusal.  ISOMIP+ prescribes the LINEAR EOS
         ! (Asay-Davis et al. 2016 Table 4), and there the constants ARE
         ! that EOS's exact coefficients, so every shipped cavity namelist
         ! is unaffected and must keep running untouched.  Under Wright or
         ! Roquet the pair is a constant stand-in for a coefficient that
         ! collapses toward zero at the freezing point and roughly doubles
         ! by 1000 dbar — exactly the corner a cavity sits in.
         if (trim(adjustl(cfg%ocean%vmix%buoyancy_coeffs)) == "constant" .and. &
             trim(adjustl(cfg%ocean%eos%eos)) /= "linear") then
            call logger%warning("&ocean_cavity_melt_nml enable=.true. with a "// &
                                "NONLINEAR EOS (&ocean_eos_nml eos='"// &
                                trim(adjustl(cfg%ocean%eos%eos))//"') but "// &
                                "&ocean_vmix_nml buoyancy_coeffs='constant'.  The "// &
                                "KPP surface buoyancy flux B_0 and the "// &
                                "double-diffusion density ratio will size the "// &
                                "melt-driven buoyancy with the SCALAR "// &
                                "&ocean_ic_nml alpha_T/beta_S, not with the "// &
                                "derivatives of the density this run actually "// &
                                "integrates.  Near the freezing point thermal "// &
                                "expansion is several times smaller than at 10 degC "// &
                                "and grows strongly with pressure, so the boundary "// &
                                "layer under the shelf can be mis-sized (and in the "// &
                                "cold-fresh corner mis-signed).  Set "// &
                                "buoyancy_coeffs='eos' unless you are reproducing a "// &
                                "constant-coefficient reference.")
         end if
      end if
      ! ---- Ice-shelf TOP drag (&ocean_tdrag_nml, Phase 4a) ----
      ! Default off ⇒ this whole block is skipped and every path is
      ! bit-identical.  Every restriction below fails loud: a top drag
      ! that is silently inert (no cover, no coefficient) is worse than
      ! no top drag, because a cavity circulation would then be
      ! frictionless at the ice and LOOK like it was damped.
      if (cfg%ocean%tdrag%enable) then
         if (.not. cfg%ocean%cavity_dyn%enable) then
            call logger%error("&ocean_tdrag_nml enable=.true. requires "// &
                              "&ocean_cavity_dyn_nml enable=.true.  The top drag "// &
                              "stands on that group's geometry: without a draft "// &
                              "there is no ice base, cover_frac is identically zero, "// &
                              "and every face mask would be zero — an inert kernel "// &
                              "with a cost.  An OPEN surface's momentum boundary "// &
                              "condition is the wind stress (&ocean_topo_nml "// &
                              "wind_config), not a drag.")
            has_error = .true.
         end if
         block
            integer :: tdrag_code
            tdrag_code = parse_tdrag_variant(cfg%ocean%tdrag%form)
            if (.not. tdrag_variant_is_implemented(tdrag_code)) then
               call logger%error("&ocean_tdrag_nml form = '"// &
                                 trim(adjustl(cfg%ocean%tdrag%form))// &
                                 "' is not recognised (quadratic|linear).  The two "// &
                                 "forms take coefficients of DIFFERENT dimensions "// &
                                 "(cd dimensionless, r in 1/s), so a typo cannot be "// &
                                 "defaulted.")
               has_error = .true.
            else if (tdrag_code == TDRAG_QUADRATIC .and. cfg%ocean%tdrag%cd <= 0.0_wp) then
               call logger%error("&ocean_tdrag_nml form='quadratic' requires cd > 0 "// &
                                 "(ISOMIP+ prescribes 2.5e-3).  cd = 0 is an "// &
                                 "identically zero drag, which is what enable=.false. "// &
                                 "is for.")
               has_error = .true.
            else if (tdrag_code == TDRAG_LINEAR .and. cfg%ocean%tdrag%r <= 0.0_wp) then
               call logger%error("&ocean_tdrag_nml form='linear' requires r > 0 "// &
                                 "(1/s).  r = 0 is an identically zero drag, which "// &
                                 "is what enable=.false. is for.")
               has_error = .true.
            end if
            ! ---- ONE drag coefficient for momentum and melt ----
            ! MOM6 carries two independent top-drag coefficients (one in
            ! the momentum BC, one in the melt u*); we deliberately do
            ! not.  A cavity in which the ice base takes momentum out of
            ! the flow at one C_d and reports a friction velocity built
            ! on another is not a closure, it is two closures sharing a
            ! boundary.  So: when both groups are on, the melt slot TAKES
            ! its C_d from this group (`configure_ocean_cavity_melt`), and
            ! a user who set both to DIFFERENT values is told rather than
            ! silently overridden.
            if (cfg%ocean%cavity_melt%enable .and. tdrag_code == TDRAG_QUADRATIC) then
               if (cfg%ocean%cavity_melt%cdrag_top /= cfg%ocean%tdrag%cd) then
                  call logger%error("&ocean_tdrag_nml cd = "// &
                                    to_string(cfg%ocean%tdrag%cd)//" and "// &
                                    "&ocean_cavity_melt_nml cdrag_top = "// &
                                    to_string(cfg%ocean%cavity_melt%cdrag_top)// &
                                    " disagree.  There is ONE ice-base drag "// &
                                    "coefficient in this model: the same C_d sets "// &
                                    "the momentum sink and the melt friction "// &
                                    "velocity u* = sqrt(C_d*(U^2 + u_tide^2)).  Set "// &
                                    "them equal (or leave cdrag_top at its default "// &
                                    "and set &ocean_tdrag_nml cd alone).")
                  has_error = .true.
               end if
            end if
            if (cfg%ocean%cavity_melt%enable .and. tdrag_code == TDRAG_LINEAR) then
               call logger%warning("&ocean_tdrag_nml form='linear' with "// &
                                   "&ocean_cavity_melt_nml enable=.true.: the melt "// &
                                   "friction velocity is quadratic by construction "// &
                                   "(u* = sqrt(C_d*(U^2 + u_tide^2))), so it keeps "// &
                                   "its own cdrag_top and the momentum sink uses r. "// &
                                   "The two boundary conditions are then NOT the "// &
                                   "same closure — intended for analytic work only.")
            end if
         end block
      end if
      ! ---- Dynamic wetting/drying v1 scope (docs/ocean_wetdry_plan.md §6) ----
      ! Every restriction fails loud: silently running wet/dry outside its
      ! validated envelope is the coastal ZSTAR_FULL salt-leak foot-gun class.
      if (cfg%ocean%wetdry%enable) then
         if (.not. (cfg%ocean%wetdry%rewet_depth > cfg%ocean%wetdry%dry_depth &
                    .and. cfg%ocean%wetdry%dry_depth > 0.0_wp)) then
            call logger%error("&ocean_wetdry_nml requires rewet_depth > dry_depth > 0 "// &
                              "(hysteresis band; equal thresholds flip-flop the front)")
            has_error = .true.
         end if
         ! v1 vertical-coordinate restriction: sigma / zstar-lite only.  On
         ! sigma the layer partition scales all layers to zero TOGETHER as
         ! D -> 0, so a dry column is just "all layers vanished" (a state the
         ! H_VANISHED gates already handle).  ZSTAR_FULL bed layers pinch
         ! independently of the surface — the argument fails and the
         ! documented 1-2%/cycle intertidal salt leak would compound.
         ! Parse (not string-compare): the vcoord string has aliases
         ! ("z-star", "SIGMA", ...); same ocean default as rdb_ocean_vcoord.
         block
            integer :: wd_vcoord_code
            wd_vcoord_code = parse_vcoord_type(cfg%vcoord_type, &
                                               default_code=VCOORD_EULERIAN_Z)
            if (.not. (wd_vcoord_code == VCOORD_SIGMA .or. &
                       wd_vcoord_code == VCOORD_ZSTAR)) then
               call logger%error("&ocean_wetdry_nml enable=.true. supports "// &
                                 "vcoord_type='sigma' or 'zstar' (zstar-lite) only — "// &
                                 "got '"//trim(cfg%vcoord_type)//"' (ZSTAR_FULL bed "// &
                                 "layers pinch independently; other coords unvalidated)")
               has_error = .true.
            end if
         end block
         ! Positive-definite continuity composes its own per-layer outflux
         ! limiter; wet/dry owns the barotropic drying limiter.  The two
         ! positivity machineries have not been reconciled (v1) — fail loud
         ! rather than run both over the same faces.
         if (cfg%ocean%continuity%positive_definite) then
            call logger%error("&ocean_wetdry_nml enable=.true. is mutually "// &
                              "exclusive with &ocean_continuity_nml "// &
                              "positive_definite=.true. (wet/dry owns its own "// &
                              "barotropic limiter; composition deferred)")
            has_error = .true.
         end if
         ! Layer-side positivity is the continuity-PPM limiter's job — the BT
         ! limiter bounds the column, ppm_limit_pos bounds each layer.
         if (.not. cfg%ocean%continuity%ppm_limit_pos) then
            call logger%error("&ocean_wetdry_nml enable=.true. requires "// &
                              "&ocean_continuity_nml ppm_limit_pos=.true. "// &
                              "(per-layer positivity under drying transports)")
            has_error = .true.
         end if
         ! The BT_cont / upstream-h Pass-1 flux branches are not composed
         ! with the wet/dry upwind+limiter branch in v1 (each replaces the
         ! same face-thickness logic).
         if (cfg%ocean%bt%use_cont_type .or. cfg%ocean%bt%upstream_h_face) then
            call logger%error("&ocean_wetdry_nml enable=.true. is mutually "// &
                              "exclusive with &ocean_bt_nml use_cont_type / "// &
                              "upstream_h_face (all three replace the BT Pass-1 "// &
                              "face thickness; composition deferred)")
            has_error = .true.
         end if
         ! v1 single-rank: the theta/mask ghost exchange lands with the
         ! C-grid MPI halo work.
         if (cfg%px*cfg%py > 1) then
            call logger%error("&ocean_wetdry_nml enable=.true. is single-rank in "// &
                              "v1 (px*py = 1); the wet-mask/limiter halo exchange "// &
                              "is deferred to the C-grid MPI work")
            has_error = .true.
         end if
         ! Wet/dry needs the split solver: the limiter lives in the BT substep.
         if (cfg%ocean%bt%n_inner < 1 .and. .not. cfg%ocean%bt%auto_n_inner) then
            call logger%error("&ocean_wetdry_nml enable=.true. requires the "// &
                              "split-explicit solver (&ocean_bt_nml n_inner >= 1 "// &
                              "or auto_n_inner=.true.)")
            has_error = .true.
         end if
         ! Surface-flux composition: only the MAIN heat/salt deposit is
         ! dynamic-mask gated in v1.  sw_pen removes the surface deposit
         ! and re-adds a distributed profile — on a masked (dry) column
         ! it would remove heat that was never added; the restoring
         ! piston rate lambda = piston/h_top and geothermal's bed deposit
         ! blow up on a ~0-thickness sliver.  All three mis-compose:
         ! fail loud, defer the gated variants.
         if (cfg%ocean%thermo%sw_pen_frac > 0.0_wp) then
            call logger%error("&ocean_wetdry_nml enable=.true. does not yet "// &
                              "compose with shortwave penetration "// &
                              "(&ocean_thermo_nml sw_pen_frac > 0): the additive "// &
                              "correction would withdraw a surface deposit the "// &
                              "dry-column mask never made")
            has_error = .true.
         end if
         if (cfg%ocean%restore%enable_restore_temp .or. &
             cfg%ocean%restore%enable_restore_salt) then
            call logger%error("&ocean_wetdry_nml enable=.true. does not yet "// &
                              "compose with surface restoring (&ocean_restore_nml): "// &
                              "the piston rate lambda = piston/h_top diverges on a "// &
                              "drying column")
            has_error = .true.
         end if
         if (cfg%ocean%geothermal%enable) then
            call logger%error("&ocean_wetdry_nml enable=.true. does not yet "// &
                              "compose with geothermal heating "// &
                              "(&ocean_geothermal_nml): the bed deposit blows up "// &
                              "a dried column's vanished bed layer")
            has_error = .true.
         end if
      end if
      ! ---- Pseudo-salt: fail-loud exclusions ----
      ! Pseudo-salt's contract is "receives EXACTLY the operators salinity
      ! receives" (PR-28).  SSS piston restoring and sea-ice frazil/basal
      ! salt exchange are salinity sources this PR does NOT mirror into
      ! pseudo-salt — an un-mirrored source would make the deviation
      ! diagnostic (pseudo_salt - S) measure that source instead of the
      ! passive-transport-path error it exists to isolate.  Abort rather
      ! than silently ship a lying diagnostic.
      if (pseudo_salt_conflicts_restore(cfg%ocean%tracers%enable_pseudo_salt, &
                                        cfg%ocean%restore%enable_restore_salt)) then
         call logger%error("&ocean_tracers_nml enable_pseudo_salt=.true. does not "// &
                           "compose with &ocean_restore_nml enable_restore_salt=.true. "// &
                           "(SSS piston restoring): the deviation diagnostic would "// &
                           "measure the un-mirrored restoring term, not the "// &
                           "passive-transport-path error")
         has_error = .true.
      end if
      if (pseudo_salt_conflicts_ice(cfg%ocean%tracers%enable_pseudo_salt, &
                                    cfg%ocean%ice%enable)) then
         call logger%error("&ocean_tracers_nml enable_pseudo_salt=.true. does not "// &
                           "compose with &ocean_ice_nml enable=.true.: sea-ice "// &
                           "frazil/basal salt exchange is an un-mirrored salinity "// &
                           "source that would corrupt the deviation diagnostic")
         has_error = .true.
      end if
      if (pseudo_salt_needs_thermo_warning(cfg%ocean%tracers%enable_pseudo_salt, &
                                           cfg%ocean%thermo%enable_thermodynamics)) then
         call logger%warning("&ocean_tracers_nml enable_pseudo_salt=.true. with "// &
                             "&ocean_thermo_nml enable_thermodynamics=.false.: "// &
                             "pseudo-salt will never receive the surface salt flux "// &
                             "or KPP non-local mirrors, so it measures pure "// &
                             "advection/diffusion transport only — a legitimate "// &
                             "but narrower use of the diagnostic")
      end if
      ! ---- Sea-ice: whole-slot envelope (applies to enable=.true.) ----
      ! The sea-ice slot carries NO cross-rank halo exchange anywhere — not
      ! in the column thermodynamics, the ITD transport, the EVP dynamics,
      ! nor the ocean coupler.  It is single-rank by design in v1 (the
      ! C-grid MPI halo is deferred).  Guard at the `enable` level so the
      ! fail-loud covers thermo-only configurations too, not just the
      ! transport/dynamics sub-features guarded below.
      if (cfg%ocean%ice%enable) then
         if (cfg%px*cfg%py > 1) then
            call logger%error("&ocean_ice_nml enable=.true. is single-rank in v1 "// &
                              "(px*py = 1); the sea-ice path carries no cross-rank "// &
                              "halo exchange — MPI support is deferred to the C-grid "// &
                              "MPI work")
            has_error = .true.
         end if
         ! Sea ice is tuned for global, ice-covered basins — not for
         ! regional/coastal runs.  It runs every outer step and is NOT
         ! performance-tuned; expect a large throughput hit.  Warn loudly
         ! so it is never mistaken for a free rider on a coastal run.
         call logger%warning("=====================================================")
         call logger%warning("&ocean_ice_nml enable=.true.: SEA ICE IS ON.")
         call logger%warning("Sea ice is tuned for GLOBAL, ice-covered simulations "// &
                             "and is NOT performance-tuned. Expect it to be VERY "// &
                             "SLOW on regional/coastal domains where most of the "// &
                             "grid is ice-free.")
         call logger%warning("=====================================================")
      end if
      ! ---- Sea-ice PR 4b: category ice transport envelope ----
      ! `transport=.true.` is a request to run `ice_transport_step` every
      ! thermo window; the v1 envelope excludes configurations the kernel
      ! does not (yet) handle — fail loud rather than silently produce a
      ! wrong/unstable answer.
      if (cfg%ocean%ice%transport) then
         if (.not. cfg%ocean%ice%enable) then
            call logger%error("&ocean_ice_nml transport=.true. requires "// &
                              "enable=.true. (the ice slot must be live)")
            has_error = .true.
         end if
         if (cfg%ocean%ice%ncat <= 1) then
            call logger%error("&ocean_ice_nml transport=.true. requires ncat > 1 "// &
                              "(ncat==1 is the frozen legacy lumped mode — no ITD, "// &
                              "nothing to transport by category)")
            has_error = .true.
         end if
         ! The driver only reaches `ice_transport_step` inside the
         ! `enable_thermodynamics .and. is_thermo_step()` ice block, so
         ! `transport=.true.` with thermo off is a WORDLESS no-op — fail
         ! loud (house rule) instead.
         if (.not. cfg%ocean%thermo%enable_thermodynamics) then
            call logger%error("&ocean_ice_nml transport=.true. requires "// &
                              "&ocean_thermo_nml enable_thermodynamics=.true. "// &
                              "(the transport call fires only on the thermo cadence)")
            has_error = .true.
         end if
         ! v1 envelope: closed/land boundaries only.  Ice ghost cells are
         ! never wrapped, so a periodic ocean edge would silently transport
         ! ice into/out of an unmirrored ghost band.
         if (ocean_bc_type_from_string(cfg%ocean%bc%west) == OBC_PERIODIC .or. &
             ocean_bc_type_from_string(cfg%ocean%bc%east) == OBC_PERIODIC .or. &
             ocean_bc_type_from_string(cfg%ocean%bc%south) == OBC_PERIODIC .or. &
             ocean_bc_type_from_string(cfg%ocean%bc%north) == OBC_PERIODIC) then
            call logger%error("&ocean_ice_nml transport=.true. is incompatible with "// &
                              "any PERIODIC &ocean_bc_nml edge (v1 envelope: closed/"// &
                              "land boundaries only — ice ghost cells are never wrapped)")
            has_error = .true.
         end if
         ! Single-rank only in v1 (mirrors the wetdry / semi_implicit
         ! precedent above): the category face-flux workspace carries no
         ! cross-rank halo exchange yet.
         if (cfg%px*cfg%py > 1) then
            call logger%error("&ocean_ice_nml transport=.true. is single-rank in v1 "// &
                              "(px*py = 1); the per-category face-flux halo exchange "// &
                              "is deferred to the C-grid MPI work")
            has_error = .true.
         end if
         if (cfg%ocean%ice%adv_substeps < 1) then
            call logger%error("&ocean_ice_nml adv_substeps must be >= 1")
            has_error = .true.
         end if
      end if

      ! ---- Sea-ice PR 5: C-grid EVP dynamics envelope ----
      ! `dynamics=.true.` requests `ice_evp_step` every outer step; the v1
      ! envelope excludes configurations the kernel does not (yet) handle —
      ! fail loud rather than silently produce a wrong/unstable answer.
      ! Unlike transport, EVP DOES allow periodic edges (the ghost-wrap
      ! machinery in `rdb_ice_evp` mirrors the ocean's own periodic-wrap
      ! contract), so this is a separate (looser on periodicity, otherwise
      ! similar) envelope, not a re-use of the transport block above.
      if (cfg%ocean%ice%dynamics) then
         if (.not. cfg%ocean%ice%enable) then
            call logger%error("&ocean_ice_nml dynamics=.true. requires "// &
                              "enable=.true. (the ice slot must be live)")
            has_error = .true.
         end if
         ! Single-rank only in v1: the EVP substep loop's ghost-wrap
         ! machinery carries no cross-rank halo exchange yet.
         if (cfg%px*cfg%py > 1) then
            call logger%error("&ocean_ice_nml dynamics=.true. is single-rank in v1 "// &
                              "(px*py = 1); the EVP substep halo exchange is deferred "// &
                              "to the C-grid MPI work")
            has_error = .true.
         end if
         ! v1 envelope: every edge must be WALL or PERIODIC — no OBC/tidal/
         ! sponge/clamped/Chapman edges for ice (the momentum solve has no
         ! notion of an open ice boundary yet).
         if (.not. any(ocean_bc_type_from_string(cfg%ocean%bc%west) == &
                       [OBC_WALL, OBC_PERIODIC]) .or. &
             .not. any(ocean_bc_type_from_string(cfg%ocean%bc%east) == &
                       [OBC_WALL, OBC_PERIODIC]) .or. &
             .not. any(ocean_bc_type_from_string(cfg%ocean%bc%south) == &
                       [OBC_WALL, OBC_PERIODIC]) .or. &
             .not. any(ocean_bc_type_from_string(cfg%ocean%bc%north) == &
                       [OBC_WALL, OBC_PERIODIC])) then
            call logger%error("&ocean_ice_nml dynamics=.true. requires every "// &
                              "&ocean_bc_nml edge to be 'wall' or 'periodic' (v1: no "// &
                              "OBC/tidal/sponge/clamped/Chapman edges for ice dynamics)")
            has_error = .true.
         end if
         ! No tripolar north fold: the EVP ghost-wrap contract is plain
         ! periodic/wall only.
         if (trim(cfg%ocean%bc%north) == "tripolar_fold") then
            call logger%error("&ocean_ice_nml dynamics=.true. is incompatible with "// &
                              "north='tripolar_fold' (v1: plain periodic/wall ghost "// &
                              "policy only)")
            has_error = .true.
         end if
         if (cfg%ocean%ice%evp_sub_steps < 1) then
            call logger%error("&ocean_ice_nml evp_sub_steps must be >= 1")
            has_error = .true.
         end if
         ! Physical-parameter positivity guards — fail loud on a typo/garbage
         ! value that would silently produce a wrong or NaN rheology. `ec`
         ! allows 0 (the documented cavitating-fluid mode: `i_ec2` is forced
         ! to 0 in `rdb_ice_evp`), so it is guarded `>= 0`, not `> 0`.
         ! `tdamp` is legitimately sign-free (`< 0` selects the
         ! fraction-of-slow-step form) — no guard.
         if (cfg%ocean%ice%p0 <= 0.0_wp) then
            call logger%error("&ocean_ice_nml p0 (ice-strength P*) must be > 0")
            has_error = .true.
         end if
         if (cfg%ocean%ice%c0 <= 0.0_wp) then
            call logger%error("&ocean_ice_nml c0 (ice-strength C*) must be > 0")
            has_error = .true.
         end if
         if (cfg%ocean%ice%ec < 0.0_wp) then
            call logger%error("&ocean_ice_nml ec (yield ellipticity) must be >= 0 "// &
                              "(0 = cavitating-fluid rheology)")
            has_error = .true.
         end if
         if (cfg%ocean%ice%cdw <= 0.0_wp) then
            call logger%error("&ocean_ice_nml cdw (ice-ocean drag coefficient) must be > 0")
            has_error = .true.
         end if
         if (cfg%ocean%ice%rho_ocean <= 0.0_wp) then
            call logger%error("&ocean_ice_nml rho_ocean (ice-drag reference density) "// &
                              "must be > 0")
            has_error = .true.
         end if
         if (cfg%ocean%ice%del_sh_min_scale <= 0.0_wp) then
            call logger%error("&ocean_ice_nml del_sh_min_scale (viscosity-floor scale) "// &
                              "must be > 0")
            has_error = .true.
         end if
         ! Drift-only mode (dynamics without transport) has no CFL backstop
         ! on u_ice: transport's positivity abort is the only ice-velocity
         ! guard. Warn (NOT error — the envelope is unchanged); this
         ! combination is unvalidated.
         if (.not. cfg%ocean%ice%transport) then
            call logger%warning("&ocean_ice_nml dynamics=.true. without "// &
                                "transport=.true.: drift-only mode has no CFL backstop "// &
                                "on u_ice (transport's positivity abort is the only "// &
                                "ice-velocity guard) and is unvalidated")
         end if
         ! PR 36: the per-iteration clip needs a ceiling to clip to — SIS2's
         ! own triple gate (`do_trunc_its = CFL_check_its .and. CFL_trunc>0
         ! .and. dt_slow>0`) makes this combination a silent no-op; fail
         ! loud instead of silently reproducing SIS2's silence.
         if (cfg%ocean%ice%cfl_trunc_dyn_its .and. cfg%ocean%ice%cfl_trunc <= 0.0_wp) then
            call logger%error("&ocean_ice_nml cfl_trunc_dyn_its=.true. requires "// &
                              "cfl_trunc > 0 (the per-iteration check needs a ceiling)")
            has_error = .true.
         end if
         ! SIS2 permits cfl_trunc > 1 but calls it unwise ("instability can
         ! occur past 0.5") — warn, do not fail loud.
         if (cfg%ocean%ice%cfl_trunc > 1.0_wp) then
            call logger%warning("&ocean_ice_nml cfl_trunc > 1.0 is permitted but unwise "// &
                                "(SIS2: 'instability can occur past 0.5')")
         end if
      end if

      ! ---- Sea-ice PR 36: cfl_trunc / cfl_trunc_dyn_its / project_ci ----
      ! require dynamics; a knob that only has meaning inside a disabled
      ! feature is a wordless no-op — fail loud instead (same pattern as
      ! a_face_stress below). Deliberately OUTSIDE the `if (dynamics)`
      ! envelope above: it must fire precisely when `dynamics` is false.
      if (cfg%ocean%ice%cfl_trunc > 0.0_wp .and. .not. cfg%ocean%ice%dynamics) then
         call logger%error("&ocean_ice_nml cfl_trunc > 0 requires dynamics=.true. "// &
                           "(the knob only reaches the EVP velocity solve)")
         has_error = .true.
      end if
      if (cfg%ocean%ice%cfl_trunc_dyn_its .and. .not. cfg%ocean%ice%dynamics) then
         call logger%error("&ocean_ice_nml cfl_trunc_dyn_its=.true. requires "// &
                           "dynamics=.true. (the knob only reaches the EVP velocity solve)")
         has_error = .true.
      end if
      if (cfg%ocean%ice%project_ci .and. .not. cfg%ocean%ice%dynamics) then
         call logger%error("&ocean_ice_nml project_ci=.true. requires dynamics=.true. "// &
                           "(the knob only reaches the EVP subcycle loop)")
         has_error = .true.
      end if

      ! ---- Map-driven sponge (&ocean_sponge_nml, PR-23) ----
      if (cfg%ocean%sponge%enable) then
         if (.not. sponge_source_is_implemented(cfg%ocean%sponge%damp_source)) then
            call logger%error("&ocean_sponge_nml enable=.true. with unimplemented "// &
                              "damp_source = '"//trim(cfg%ocean%sponge%damp_source)// &
                              "': must be 'band' ('file' is PR-23b, needs the PR-14 reader)")
            has_error = .true.
         end if
         if (.not. sponge_target_is_implemented(cfg%ocean%sponge%target_source)) then
            call logger%error("&ocean_sponge_nml enable=.true. with unimplemented "// &
                              "target_source = '"//trim(cfg%ocean%sponge%target_source)// &
                              "': must be 'ic' or 'linear_z' "// &
                              "('file' is PR-23b, needs the PR-14 reader)")
            has_error = .true.
         end if
         if (.not. sponge_ramp_is_valid(cfg%ocean%sponge%ramp)) then
            call logger%error("&ocean_sponge_nml ramp = '"//trim(cfg%ocean%sponge%ramp)// &
                              "': must be 'cosine' (default, the legacy band shape) "// &
                              "or 'linear' (ISOMIP+ Eq. 20)")
            has_error = .true.
         end if
         ! A `linear_z` target with both gradients AND both references left
         ! at their defaults would relax toward T = 0 degC / S = 35 PSU
         ! everywhere — almost certainly not what was meant, and silent.
         if (trim(cfg%ocean%sponge%target_source) == "linear_z" .and. &
             cfg%ocean%sponge%lin_t_ref == 0.0_wp .and. &
             cfg%ocean%sponge%lin_dt_dz == 0.0_wp .and. &
             cfg%ocean%sponge%lin_ds_dz == 0.0_wp) then
            call logger%warning("&ocean_sponge_nml target_source='linear_z' with "// &
                                "lin_t_ref = lin_dt_dz = lin_ds_dz = 0: the sponge will "// &
                                "relax toward a uniform T = 0 degC column. Set the "// &
                                "lin_* profile knobs.")
         end if
         ! relax_h (interior-interface thickness damping) is NOT implemented
         ! in PR-23 v1 — deferred to PR-23b alongside the file targets (the
         ! plan's §14 Q1 recommendation). Abort unconditionally rather than
         ! silently ignore the knob.
         if (cfg%ocean%sponge%relax_h) then
            call logger%error("&ocean_sponge_nml relax_h=.true. is not implemented "// &
                              "in PR-23 v1 (interior-interface thickness damping is "// &
                              "deferred to PR-23b; see docs/plans/PLAN_PR23_real_sponge.md "// &
                              "§14 Q1) — leave relax_h=.false.")
            has_error = .true.
         end if
         ! The sponge is dispatched only from run_stage_split (the split
         ! barotropic solver); the unsplit run_stage never receives bc or sp.
         ! With auto_n_inner=.true. n_inner resolves later at setup, so only
         ! a HARD n_inner < 1 with auto off is an error here (mirrors the
         ! &ocean_wetdry_nml n_inner precedent above).
         if (cfg%ocean%bt%n_inner < 1 .and. .not. cfg%ocean%bt%auto_n_inner) then
            call logger%error("&ocean_sponge_nml enable=.true. requires the "// &
                              "split-explicit solver (&ocean_bt_nml n_inner >= 1 "// &
                              "or auto_n_inner=.true.) — the sponge is dispatched only "// &
                              "from the split barotropic path")
            has_error = .true.
         end if
         ! A "band" idamp source with no OBC_SPONGE-tagged edge would build
         ! an all-zero map — the sponge would silently do nothing.
         if (trim(cfg%ocean%sponge%damp_source) == "band") then
            if (trim(cfg%ocean%bc%west) /= "sponge" .and. &
                trim(cfg%ocean%bc%east) /= "sponge" .and. &
                trim(cfg%ocean%bc%south) /= "sponge" .and. &
                trim(cfg%ocean%bc%north) /= "sponge") then
               call logger%error("&ocean_sponge_nml enable=.true., damp_source='band' "// &
                                 "requires at least one &ocean_bc_nml edge = 'sponge' "// &
                                 "(otherwise idamp_h/u/v would build all-zero)")
               has_error = .true.
            end if
         end if
         ! Per-edge width/strength overrides: sentinel is < 0 (inherit); a
         ! non-sentinel width must be non-negative, and a positive width
         ! paired with a zero strength is a silent no-op sponge.
         if (cfg%ocean%sponge%west_width < -1 .or. cfg%ocean%sponge%east_width < -1 .or. &
             cfg%ocean%sponge%south_width < -1 .or. cfg%ocean%sponge%north_width < -1) then
            call logger%error("&ocean_sponge_nml *_width overrides must be >= -1 "// &
                              "(-1 = inherit &ocean_bc_nml sponge_width)")
            has_error = .true.
         end if
         if ((cfg%ocean%sponge%west_width > 0 .and. cfg%ocean%sponge%west_strength == 0.0_wp) .or. &
             (cfg%ocean%sponge%east_width > 0 .and. cfg%ocean%sponge%east_strength == 0.0_wp) .or. &
             (cfg%ocean%sponge%south_width > 0 .and. cfg%ocean%sponge%south_strength == 0.0_wp) .or. &
             (cfg%ocean%sponge%north_width > 0 .and. cfg%ocean%sponge%north_strength == 0.0_wp)) then
            call logger%error("&ocean_sponge_nml a positive *_width override with "// &
                              "*_strength == 0.0 would build an all-zero band on that "// &
                              "edge — set a nonzero *_strength or leave *_width at -1")
            has_error = .true.
         end if
      end if
      ! ---- Barotropic linear wave drag (&ocean_bt_nml wave_drag) ----
      ! `wave_drag_form` is an `nml_enum` (allowed = uniform/roughness_proxy/
      ! file), so a garbage tag is already rejected at nml-apply time; the
      ! one form-specific case left to catch here is "file", which is
      ! REGISTERED but has no reader yet (PR-14) — fail loud rather than
      ! silently no-op (the audit's trap #12: a named-but-dead option).
      if (trim(cfg%ocean%bt%wave_drag_form) == "file") then
         call logger%error("&ocean_bt_nml wave_drag_form='file' is not "// &
                           "implemented — the NetCDF map reader is PR-14; "// &
                           "use wave_drag_form='uniform' or 'roughness_proxy'")
         has_error = .true.
      end if
      if (cfg%ocean%bt%wave_drag_scale < 0.0_wp) then
         call logger%error("&ocean_bt_nml wave_drag_scale must be >= 0 "// &
                           "(negative r_H makes bt_rem > 1, exponentially "// &
                           "growing the barotropic mode)")
         has_error = .true.
      end if
      if (cfg%ocean%bt%wave_drag_r_uniform < 0.0_wp) then
         call logger%error("&ocean_bt_nml wave_drag_r_uniform must be >= 0")
         has_error = .true.
      end if
      if (cfg%ocean%bt%wave_drag_kappa <= 0.0_wp) then
         call logger%error("&ocean_bt_nml wave_drag_kappa must be > 0")
         has_error = .true.
      end if
      if (cfg%ocean%bt%wave_drag_n_bot < 0.0_wp) then
         call logger%error("&ocean_bt_nml wave_drag_n_bot must be >= 0")
         has_error = .true.
      end if
      if (cfg%ocean%bt%wave_drag_h2_max <= 0.0_wp) then
         call logger%error("&ocean_bt_nml wave_drag_h2_max must be > 0")
         has_error = .true.
      end if
      if (cfg%ocean%bt%wave_drag) then
         if (trim(cfg%ocean%bt%wave_drag_form) == "uniform" .and. &
             cfg%ocean%bt%wave_drag_r_uniform <= 0.0_wp) then
            call logger%warning("&ocean_bt_nml wave_drag=.true. with "// &
                                "wave_drag_form='uniform' and wave_drag_r_uniform <= 0 "// &
                                "=> inert (r_H == 0 everywhere)")
         end if
         if (trim(cfg%ocean%bt%wave_drag_form) == "roughness_proxy") then
            call logger%warning("&ocean_bt_nml wave_drag_form='roughness_proxy' derives "// &
                                "<h^2> from RESOLVED bathymetry variance — it is a "// &
                                "placeholder for a real subgrid roughness map "// &
                                "(PR-14/PR-30), not a substitute")
         end if
      end if
      ! ---- Sea-ice PR-58: ITD category-bound override envelope ----
      ! `hlim=...` is consumed entirely inside `ice%init` (called from
      ! `ocean_state_init_from_config`, BEFORE the thermo cadence), so —
      ! unlike transport/dynamics — there is no `enable_thermodynamics`
      ! clause to copy here: a wordless no-op would require the value to
      ! be read on a cadence it never reaches, and `hlim` has no such gate.
      if (ice_hlim_count(cfg%ocean%ice%hlim) > 0) then
         if (.not. cfg%ocean%ice%enable) then
            call logger%error("&ocean_ice_nml hlim=... requires enable=.true. "// &
                              "(the ice slot must be live for the override to apply)")
            has_error = .true.
         end if
         block
            logical :: hlim_ok
            character(len=:), allocatable :: reason
            call ice_hlim_spec_is_valid(cfg%ocean%ice%hlim, cfg%ocean%ice%ncat, hlim_ok, reason)
            if (.not. hlim_ok) then
               call logger%error("&ocean_ice_nml hlim: "//reason)
               has_error = .true.
            end if
         end block
      end if
      ! ---- Sea-ice PR 24: analytic initial-condition envelope ----
      ! `conc_config /= "zero"` is a request to seed a live ice pack
      ! before the first step; the envelope excludes configurations the
      ! seeder cannot handle consistently — fail loud rather than
      ! silently seed a meaningless or self-melting pack.
      if (trim(cfg%ocean%ice_ic%conc_config) /= "zero") then
         ! V1: mirrors "transport=.true. requires enable=.true." above.
         if (.not. cfg%ocean%ice%enable) then
            call logger%error("&ocean_ice_ic_nml conc_config /= 'zero' requires "// &
                              "&ocean_ice_nml enable=.true. (the ice slot must be live)")
            has_error = .true.
         end if
         ! V2: an IC with no thickness is a wordless no-op; h_ice > 0
         ! also keeps the pack above h_lim(1) = 1e-10 m so the ITD cat-1
         ! compress cannot fire on the very first thermo window.
         if (cfg%ocean%ice_ic%h_ice <= 0.0_wp) then
            call logger%error("&ocean_ice_ic_nml conc_config /= 'zero' requires "// &
                              "h_ice > 0 (an IC with no thickness is a wordless no-op, "// &
                              "and h_ice must clear the category-1 floor)")
            has_error = .true.
         end if
         ! V3: mirrors the "wordless no-op => fail loud" argument used
         ! for transport above.
         if (trim(cfg%ocean%ice_ic%conc_config) == "uniform" .and. &
             cfg%ocean%ice_ic%conc <= 0.0_wp) then
            call logger%error("&ocean_ice_ic_nml conc_config='uniform' requires "// &
                              "conc > 0 (conc <= 0 is a wordless open-water no-op)")
            has_error = .true.
         end if
         ! V4: copy of the equilibrium-tide cartesian-grid guard above —
         ! geolatT/geolonT stay 0 on a cartesian grid (only
         ! metrics_init_spherical fills them), so a latitudes IC would
         ! silently ice the whole basin or none of it.
         if (trim(cfg%ocean%ice_ic%conc_config) == "latitudes" .and. &
             trim(cfg%ocean%grid%grid_config) == "cartesian") then
            call logger%error("&ocean_ice_ic_nml conc_config='latitudes' requires a "// &
                              "non-cartesian &ocean_grid_nml grid_config (spherical/"// &
                              "tripolar/supergrid) — geolatT is meaningless on cartesian")
            has_error = .true.
         end if
         ! V5: the legacy lumped mode (ncat==1) has no partial-cover
         ! representation — ice_cell_concentration_impl returns a BINARY
         ! ci, so a fractional conc would be silently reinterpreted.
         if (cfg%ocean%ice%ncat == 1 .and. &
             trim(cfg%ocean%ice_ic%conc_config) == "uniform" .and. &
             cfg%ocean%ice_ic%conc /= 1.0_wp) then
            call logger%error("&ocean_ice_ic_nml conc_config='uniform' with "// &
                              "&ocean_ice_nml ncat==1 requires conc == 1.0 (the legacy "// &
                              "lumped mode has no partial-cover representation — "// &
                              "ice_cell_concentration_impl returns a BINARY ci)")
            has_error = .true.
         end if
         ! V6: at/above the liquidus, ice_enth_from_ts returns the
         ! LIQUID-WATER branch — the pack would carry m_ice > 0 whose
         ! enthalpy says "this is water" and melt out on the first
         ! thermo step, silently.
         if (cfg%ocean%ice_ic%t_ice >= ice_t_freeze(cfg%ocean%ice_ic%s_ice)) then
            call logger%error("&ocean_ice_ic_nml conc_config /= 'zero' requires "// &
                              "t_ice < the liquidus at s_ice (t_ice >= t_freeze(s_ice) "// &
                              "silently yields LIQUID-WATER enthalpy — the pack would "// &
                              "melt out on the first thermo step)")
            has_error = .true.
         end if
         ! V7: orphan snow. The column entry gate is on m_ice, and
         ! ice_adjust_categories' massless-cleanup step would silently
         ! absorb snow-with-no-ice into the "massless" taxonomy.
         if (cfg%ocean%ice_ic%h_snow > 0.0_wp .and. cfg%ocean%ice_ic%h_ice <= 0.0_wp) then
            call logger%error("&ocean_ice_ic_nml h_snow > 0 requires h_ice > 0 "// &
                              "(orphan snow with no ice is silently absorbed as "// &
                              "'massless' by the ITD restore)")
            has_error = .true.
         end if
         ! V8: belt-and-braces on top of nml_enum allowed= — catches an
         ! unrecognised conc_config string even if the schema layer is
         ! ever bypassed.
         if (ice_ic_parse_conc_config(cfg%ocean%ice_ic%conc_config) == ICE_IC_CONC_INVALID) then
            call logger%error("&ocean_ice_ic_nml conc_config='"// &
                              trim(cfg%ocean%ice_ic%conc_config)// &
                              "' is not recognised (allowed: zero, uniform, latitudes)")
            has_error = .true.
         end if
         ! V9: the IC is skipped on resume (a restart REPLACES it, never
         ! merges) — warn, don't abort, so the legitimate "same nml, cold
         ! then warm" workflow still works.
         if (len_trim(cfg%restart_file) > 0) then
            call logger%warning("&ocean_ice_ic_nml conc_config /= 'zero' with "// &
                                "restart_file set: the ice IC is SKIPPED on a warm "// &
                                "restart (the restart replaces it, it does not merge "// &
                                "with it) — this namelist's ice IC knobs have no effect "// &
                                "on this run")
         end if
      end if
      ! ---- Sea-ice PR 26: snowfall source-term envelope ----
      ! `snowfall > 0` is a request for `ice_atm_forcing_restoring` +
      ! `ice_thermo_driver_step`/`ice_snowfall_ocean_share` to fire every
      ! thermo window; both live inside the `enable_thermodynamics .and.
      ! is_thermo_step()` ice block in the driver, so a config with the
      ! ice slot or thermodynamics off would leave `snowfall` a wordless
      ! no-op — fail loud instead (house rule, mirrors the transport
      ! envelope above).
      if (cfg%ocean%ice%snowfall > 0.0_wp) then
         if (.not. cfg%ocean%ice%enable) then
            call logger%error("&ocean_ice_nml snowfall > 0 requires "// &
                              "enable=.true. (the ice slot must be live)")
            has_error = .true.
         end if
         if (.not. cfg%ocean%thermo%enable_thermodynamics) then
            call logger%error("&ocean_ice_nml snowfall > 0 requires "// &
                              "&ocean_thermo_nml enable_thermodynamics=.true. "// &
                              "(the atmospheric-forcing fill and the column "// &
                              "snow add both fire only on the thermo cadence)")
            has_error = .true.
         end if
      end if

      ! ---- Sea-ice PR 62: a_face_stress requires dynamics ----
      ! `a_face_stress` only reaches the EVP momentum solve — with
      ! `dynamics=.false.` it would be a wordless no-op.  Fail loud instead
      ! (house rule).  Deliberately OUTSIDE the `if (dynamics)` envelope
      ! above: it must fire precisely when `dynamics` is false.
      if (cfg%ocean%ice%a_face_stress .and. .not. cfg%ocean%ice%dynamics) then
         call logger%error("&ocean_ice_nml a_face_stress=.true. requires dynamics=.true. "// &
                           "(the knob only reaches the EVP momentum solve)")
         has_error = .true.
      end if

      ! ---- Sea-ice PR 27: Archimedes snow-ice flooding envelope ----
      ! `snow_ice=.true.` is a request for `ice_snow_ice_flood` to fire
      ! inside `ice_column_step`'s resize stage every thermo window; that
      ! call lives inside the `enable_thermodynamics .and. is_thermo_
      ! step()` ice block in the driver, so a config with the ice slot or
      ! thermodynamics off would leave `snow_ice` a wordless no-op — fail
      ! loud instead (house rule, mirrors the transport/snowfall envelopes
      ! above).
      if (cfg%ocean%ice%snow_ice) then
         if (.not. cfg%ocean%ice%enable) then
            call logger%error("&ocean_ice_nml snow_ice=.true. requires "// &
                              "enable=.true. (the ice slot must be live)")
            has_error = .true.
         end if
         if (.not. cfg%ocean%thermo%enable_thermodynamics) then
            call logger%error("&ocean_ice_nml snow_ice=.true. requires "// &
                              "&ocean_thermo_nml enable_thermodynamics=.true. "// &
                              "(the column runs only on the thermo cadence)")
            has_error = .true.
         end if
      end if

      ! ---- BT wide-halo march-in (bt_halo > 0) exclusions ----
      if (cfg%ocean%bt%bt_halo > 0) then
         if (cfg%ocean%wetdry%enable) then
            call logger%error("&ocean_bt_nml bt_halo > 0 is mutually exclusive "// &
                              "with &ocean_wetdry_nml enable (wd arrays not widened in v1)")
            has_error = .true.
         end if
         if (cfg%ocean%bt%use_cont_type) then
            call logger%error("&ocean_bt_nml bt_halo > 0 is mutually exclusive "// &
                              "with use_cont_type (BTCL arrays not widened in v1)")
            has_error = .true.
         end if
         if (cfg%ocean%bt%upstream_h_face) then
            call logger%error("&ocean_bt_nml bt_halo > 0 is mutually exclusive "// &
                              "with upstream_h_face (upstream h-face arrays not widened in v1)")
            has_error = .true.
         end if
         if (cfg%ocean%tides%enable) then
            call logger%error("&ocean_bt_nml bt_halo > 0 is mutually exclusive "// &
                              "with &ocean_tides_nml enable (wide eta_forcing copy deferred)")
            has_error = .true.
         end if
         if (cfg%ocean%psurf%enable) then
            call logger%error("&ocean_bt_nml bt_halo > 0 is mutually exclusive "// &
                              "with &ocean_psurf_nml enable (wide eta_forcing copy deferred)")
            has_error = .true.
         end if
         if (cfg%ocean%porous%enable) then
            ! `bt_wide` carries its OWN `metrics_w`, re-filled from the
            ! grid formula on the widened grid.  Nothing fills its porous
            ! statistics, so `bt_wide_substep` would hand the fast loop the
            ! UN-narrowed `metrics_w%dy_cu` / `dx_cv` — silently reverting
            ! the porous-aware barotropic solve and giving an answer that
            ! is neither the porous one nor the baseline.  Fail loud until
            ! the wide shadow carries `dy_cu_bt` / `dx_cv_bt` too.
            call logger%error("&ocean_bt_nml bt_halo > 0 is mutually exclusive "// &
                              "with &ocean_porous_nml enable (the wide-halo "// &
                              "metrics shadow carries no porous statistics, so the "// &
                              "BT solve would silently transport on un-narrowed widths)")
            has_error = .true.
         end if
         if (cfg%ocean%cavity_dyn%enable) then
            ! Same failure mode as porous, one level up: `metrics_w` is
            ! re-filled from the grid formula and carries no `z_draft`, so
            ! the wide fast loop would take the BED as its reference depth
            ! and solve a column that is twice as deep as the cavity's.
            ! (The cavity block above refuses this from its own side too —
            ! the knob that is "wrong" depends on which one the user meant.)
            call logger%error("&ocean_bt_nml bt_halo > 0 is mutually exclusive "// &
                              "with &ocean_cavity_dyn_nml enable (the wide-halo "// &
                              "metrics shadow carries no ice draft, so the BT solve "// &
                              "would reference the bed instead of the ice base)")
            has_error = .true.
         end if
         if (trim(cfg%ocean%grid%grid_config) == "supergrid" .or. &
             trim(cfg%ocean%grid%grid_config) == "tripolar") then
            call logger%error("&ocean_bt_nml bt_halo > 0 is mutually exclusive "// &
                              "with grid_config='"//trim(cfg%ocean%grid%grid_config)// &
                              "' (file-based metrics cannot be re-filled on a wide grid; "// &
                              "cartesian/spherical formula fills supported)")
            has_error = .true.
         end if
      end if
      ! ---- Porous barriers (&ocean_porous_nml) exclusions ----
      if (cfg%ocean%porous%enable) then
         if (cfg%ocean%wetdry%enable) then
            ! A drying column drives every layer's face thickness below
            ! H_VANISHED, where the open fractions all go to zero and the
            ! barotropic width goes with them.  That is self-consistent but
            ! wholly untested against the wet/dry outflow limiter, and the
            ! two schemes each own part of the same "can this face carry
            ! transport" decision.  Refuse the combination rather than ship
            ! an unvalidated interaction.
            call logger%error("&ocean_porous_nml enable is mutually exclusive with "// &
                              "&ocean_wetdry_nml enable (the vanishing-column "// &
                              "interaction between the open fractions and the "// &
                              "wet/dry outflow limiter is unvalidated)")
            has_error = .true.
         end if
      end if

      ! ---- Ideal-age tracer knob sanity (PR-7, fail-quiet class) ----
      ! Neither combination is unsafe, so both are warnings, not aborts —
      ! and the default (0/0) must not trip either guard.
      if (cfg%ocean%tracers%ideal_age_sfc_growth_rate /= 0.0_wp .and. &
          cfg%ocean%tracers%ideal_age_young_val == 0.0_wp) then
         call logger%warning("&ocean_tracers_nml ideal_age_sfc_growth_rate is nonzero but "// &
                             "ideal_age_young_val = 0: young_val*exp(...) is identically 0, "// &
                             "so the growth-rate knob is silently inert (MOM6 seeds its "// &
                             "vintage tracer with young_val = 1e-20 for exactly this reason)")
      end if
      if (.not. cfg%ocean%tracers%enable_ideal_age .and. &
          (cfg%ocean%tracers%ideal_age_young_val /= 0.0_wp .or. &
           cfg%ocean%tracers%ideal_age_sfc_growth_rate /= 0.0_wp)) then
         call logger%warning("&ocean_tracers_nml ideal_age_young_val / "// &
                             "ideal_age_sfc_growth_rate are set but enable_ideal_age = "// &
                             ".false.: no age tracer is registered, so both knobs are "// &
                             "silently inert")
      end if

      ! ---- Surface-flux component set (PR-12) ----
      ! The component set only feeds the thermo path (the assembler
      ! derives Q_heat/Q_salt, which apply_tracers stamps onto the
      ! tracer hTr); with thermo off the whole thing is inert. Fail loud
      ! rather than let a user carry the (real, gated-off) device memory
      ! for nothing.
      if (forcing_components_need_thermo(cfg%ocean%forcing%enable_components, &
                                         cfg%ocean%thermo%enable_thermodynamics)) then
         call logger%error("&ocean_forcing_nml enable_components=.true. requires "// &
                           "&ocean_thermo_nml enable_thermodynamics=.true. "// &
                           "(the component set feeds the tracer path; it is inert "// &
                           "under enable_thermodynamics=.false.)")
         has_error = .true.
      end if

      ! ---- &ocean_dataovr_nml (PR-15 file-backed surface forcing) ----
      !
      ! Validated HERE, not at registration, on the "a broken config must
      ! not run at all" principle: `ocean_data_forcing_configure` runs
      ! deep in the driver's setup phase, AFTER the grid, bathymetry,
      ! initial condition and a dozen slots are built.  Catching a
      ! namelist typo there means the user pays for all of that first.
      ! Everything below is decidable from `cfg` alone, so it costs
      ! nothing to decide it before any work happens.
      !
      ! The fail-loud guards inside `rdb_ocean_data_forcing` STAY: they
      ! are the gate for a programmatic caller (a future Python driver
      ! registers fields directly and never parses a namelist), and this
      ! block is unreachable for that path.
      if (cfg%ocean%dataovr%enable) then
#ifdef RDB_NO_NETCDF
         call logger%error("&ocean_dataovr_nml enable=.true. requires a NetCDF build "// &
                           "(RDB_ENABLE_NETCDF=ON); the reader is compiled out here")
         has_error = .true.
#endif
         ! D1: a named file with no variable name.  No guessed fallbacks
         ! for forcing — an unintended variable is worse than a stop.
         if (.not. dataovr_entry_is_valid(cfg%ocean%dataovr%tau_x)) then
            call logger%error("&ocean_dataovr_nml tau_x_file is set but tau_x_var is blank")
            has_error = .true.
         end if
         if (.not. dataovr_entry_is_valid(cfg%ocean%dataovr%tau_y)) then
            call logger%error("&ocean_dataovr_nml tau_y_file is set but tau_y_var is blank")
            has_error = .true.
         end if
         if (.not. dataovr_entry_is_valid(cfg%ocean%dataovr%heat)) then
            call logger%error("&ocean_dataovr_nml heat_file is set but heat_var is blank")
            has_error = .true.
         end if
         if (.not. dataovr_entry_is_valid(cfg%ocean%dataovr%evap)) then
            call logger%error("&ocean_dataovr_nml evap_file is set but evap_var is blank")
            has_error = .true.
         end if
         if (.not. dataovr_entry_is_valid(cfg%ocean%dataovr%lprec)) then
            call logger%error("&ocean_dataovr_nml lprec_file is set but lprec_var is blank")
            has_error = .true.
         end if
         if (.not. dataovr_entry_is_valid(cfg%ocean%dataovr%salt)) then
            call logger%error("&ocean_dataovr_nml salt_file is set but salt_var is blank")
            has_error = .true.
         end if

         ! D2: cyclic climatology with no period is undecidable, not
         ! defaultable — there is no sane fallback for "one year".
         if (.not. dataovr_time_is_valid(cfg%ocean%dataovr)) then
            call logger%error("&ocean_dataovr_nml time_mode='cyclic' requires cycle_period > 0")
            has_error = .true.
         end if

         ! D3: the freshwater tags exist ONLY in the surface-flux
         ! component set; without it there is no slot to write into.
         if (dataovr_freshwater_needs_components(cfg%ocean%dataovr, &
                                                 cfg%ocean%forcing%enable_components)) then
            call logger%error("&ocean_dataovr_nml evap/lprec file forcing requires "// &
                              "&ocean_forcing_nml enable_components=.true. "// &
                              "(there is no non-component freshwater slot to write into)")
            has_error = .true.
         end if

         ! D4: a named file that does not exist.  A typo'd path is the
         ! single most common way to get this group wrong, and without
         ! this check it surfaces only once registration opens the file
         ! — deep in setup, as a bare "NetCDF operation failed" that
         ! does not even name the path.  `inquire` needs no NetCDF, so
         ! this works in every build.  It is a genuine
         ! time-of-check/time-of-use race (the file could vanish before
         ! registration), which is fine: registration still fails loud,
         ! this only moves the COMMON case earlier and makes it legible.
         call check_dataovr_file(cfg%ocean%dataovr%tau_x, "tau_x", has_error)
         call check_dataovr_file(cfg%ocean%dataovr%tau_y, "tau_y", has_error)
         call check_dataovr_file(cfg%ocean%dataovr%heat, "heat", has_error)
         call check_dataovr_file(cfg%ocean%dataovr%evap, "evap", has_error)
         call check_dataovr_file(cfg%ocean%dataovr%lprec, "lprec", has_error)
         call check_dataovr_file(cfg%ocean%dataovr%salt, "salt", has_error)

         ! D5: enabled but nothing driven is a wordless no-op — the
         ! house "silent no-op => fail loud" rule (cf. V2/V3 above).
         if (.not. dataovr_any_tag_set(cfg%ocean%dataovr)) then
            call logger%error("&ocean_dataovr_nml enable=.true. but no <tag>_file is set — "// &
                              "either name a forcing file or set enable=.false.")
            has_error = .true.
         end if
      end if

      if (has_error) then
         ! NOTE: this is a terminal ROLLUP over many independent checks above
         ! (most of which log their own specific reason via `logger%error`
         ! but do not individually push to the error ring — instrumenting
         ! all ~160 of them is out of scope for this pass). This message is
         ! therefore GENERIC by construction; a caller wanting the specific
         ! reason should also look at ring index 1 (the most recently
         ! instrumented check, if any fired) rather than only index 0.
         call fail("Configuration validation failed — see errors above", ierr, &
                   OCEAN_STATUS_ERR_CONFIG_VALIDATE)
         return
      end if
      if (present(ierr)) ierr = OCEAN_STATUS_OK

   end subroutine validate_config

   ! ======================================================================
   ! &ocean_dataovr_nml predicates (PR-15).
   !
   ! Pure and non-aborting by design: a `pure` function cannot `error
   ! stop`, which is exactly what makes these testable from both sides
   ! without a subprocess/death-test harness — the house idiom, see
   ! `ice_hlim_spec_is_valid` / `zinit_dims_ok` / `data_input_dims_ok`.
   ! `validate_config` turns a `.false.` into the abort.
   ! ======================================================================

   subroutine check_dataovr_file(e, tag, has_error)
      !! Existence check for one tag's file.  NOT `pure` — `inquire` is
      !! an I/O statement, which is exactly why this is a subroutine
      !! setting `has_error` rather than another predicate: the pure
      !! ones above stay pure and independently testable, and this
      !! keeps the filesystem dependency in one visible place.
      !! Silent for a tag that names no file.
      type(dataovr_entry_config_t), intent(in) :: e
      character(len=*), intent(in) :: tag
      logical, intent(inout) :: has_error
      logical :: exists

      if (len_trim(e%file) == 0) return
      inquire (file=trim(e%file), exist=exists)
      if (.not. exists) then
         call logger%error("&ocean_dataovr_nml "//tag//"_file not found: '"// &
                           trim(e%file)//"'")
         has_error = .true.
      end if
   end subroutine check_dataovr_file

   pure function dataovr_entry_is_valid(e) result(ok)
      !! `.true.` unless a tag names a file without naming the variable
      !! inside it.  A blank `file` (tag not file-driven) is valid, and
      !! a `var` with no `file` is harmless — the tag simply never
      !! registers — so only the one combination is rejected.
      type(dataovr_entry_config_t), intent(in) :: e
      logical :: ok
      ok = .not. (len_trim(e%file) > 0 .and. len_trim(e%var) == 0)
   end function dataovr_entry_is_valid

   pure function dataovr_time_is_valid(cfg) result(ok)
      !! `.true.` unless `time_mode='cyclic'` was requested without a
      !! positive `cycle_period`.  Checked only for the cyclic mode;
      !! `cycle_period` is ignored by `linear`/`static`.
      type(ocean_dataovr_config_t), intent(in) :: cfg
      logical :: ok
      ok = .true.
      if (trim(cfg%time_mode) == "cyclic") ok = (cfg%cycle_period > 0.0_wp)
   end function dataovr_time_is_valid

   pure function dataovr_freshwater_needs_components(cfg, use_components) result(bad)
      !! `.true.` when `evap`/`lprec` are file-driven but the
      !! surface-flux component set they write into is switched off.
      !! Note the polarity: this reports the PROBLEM, not validity —
      !! named for how it reads at the call site.
      type(ocean_dataovr_config_t), intent(in) :: cfg
      logical, intent(in) :: use_components
      logical :: bad
      bad = (len_trim(cfg%evap%file) > 0 .or. len_trim(cfg%lprec%file) > 0) &
            .and. (.not. use_components)
   end function dataovr_freshwater_needs_components

   pure function dataovr_any_tag_set(cfg) result(any_set)
      !! `.true.` when at least one tag names a file, i.e. enabling the
      !! group would actually do something.
      type(ocean_dataovr_config_t), intent(in) :: cfg
      logical :: any_set
      any_set = len_trim(cfg%tau_x%file) > 0 .or. len_trim(cfg%tau_y%file) > 0 .or. &
                len_trim(cfg%heat%file) > 0 .or. len_trim(cfg%evap%file) > 0 .or. &
                len_trim(cfg%lprec%file) > 0 .or. len_trim(cfg%salt%file) > 0
   end function dataovr_any_tag_set

   pure function ice_hlim_count(hlim) result(n)
      !! Length of the LEADING run of non-sentinel (>= 0) entries in
      !! `&ocean_ice_nml hlim` (PR-58). Contiguity matters here (unlike
      !! `apply_layer_rho_init`'s bare `count(>= 0)`): a value AFTER the
      !! first sentinel is silently dropped by this count, so
      !! `ice_hlim_spec_is_valid` checks contiguity explicitly rather than
      !! trusting the count alone.
      real(wp), intent(in) :: hlim(:)
      integer :: n
      integer :: k

      n = 0
      do k = 1, size(hlim)
         if (hlim(k) < 0.0_wp) exit
         n = n + 1
      end do
   end function ice_hlim_count

   pure subroutine ice_hlim_spec_is_valid(hlim, ncat, ok, reason)
      !! Shape/monotonicity predicate for a supplied `&ocean_ice_nml hlim`
      !! list (PR-58). A SUBROUTINE, not a function — a pure FUNCTION may
      !! not carry an `intent(out)` dummy (`reason`); pure SUBROUTINES
      !! can. `n = ice_hlim_count(hlim)`; `ok` iff ALL of:
      !!   (a) n >= 2                          -- one edge gives no width
      !!                                          to extrapolate (SIS2
      !!                                          divergence D1: Roundabout
      !!                                          fails loud here instead
      !!                                          of SIS2's silent
      !!                                          fallback to the default
      !!                                          table).
      !!   (b) n <= ncat + 1                   -- more edges than the ITD
      !!                                          has bins.
      !!   (c) all(hlim(n+1:) < 0)             -- contiguity: no value
      !!                                          after the first
      !!                                          sentinel.
      !!   (d) hlim(1) > 0                     -- mh_lim(1) > 0 is a live
      !!                                          gate (`rdb_ice_itd`
      !!                                          `ice_adjust_categories`).
      !!   (e) hlim(k+1) > hlim(k), k = 1..n-1  -- strictly increasing.
      !! Checking (d)+(e) on the SUPPLIED list is sufficient for the FULL
      !! `h_lim(1..ncat+1)`: the constant-width extrapolation preserves
      !! both properties by induction
      !! (`h_lim(k)-h_lim(k-1) = h_lim(k-1)-h_lim(k-2)`), so no redundant
      !! post-extrapolation check is needed.
      real(wp), intent(in) :: hlim(:)
      integer, intent(in) :: ncat
      logical, intent(out) :: ok
      character(len=:), allocatable, intent(out) :: reason
      integer :: n, k

      n = ice_hlim_count(hlim)
      ok = .false.
      reason = ""

      if (n < 2) then
         reason = "hlim must supply >= 2 entries (a single edge has no "// &
                  "width to extrapolate); got "//to_string(n)
         return
      end if
      if (n > ncat + 1) then
         reason = "hlim supplies "//to_string(n)//" entries, more than "// &
                  "ncat+1 ("//to_string(ncat + 1)//")"
         return
      end if
      if (.not. all(hlim(n + 1:) < 0.0_wp)) then
         reason = "hlim is not a contiguous leading run: a value follows "// &
                  "the first sentinel (entry "//to_string(n + 1)//")"
         return
      end if
      if (hlim(1) <= 0.0_wp) then
         reason = "hlim(1) must be > 0"
         return
      end if
      do k = 1, n - 1
         if (.not. (hlim(k + 1) > hlim(k))) then
            reason = "hlim must be strictly increasing (violated at entry "// &
                     to_string(k + 1)//")"
            return
         end if
      end do

      ok = .true.
      reason = "ok"
   end subroutine ice_hlim_spec_is_valid

   pure function resolve_bt_halo(requested, compute_size, exclusion_active) result(width)
      !! Resolve the `&ocean_bt_nml bt_halo` sentinel to a concrete march-in
      !! width.
      !!   requested == BT_HALO_AUTO_SENTINEL (-1, the default) => AUTO:
      !!     BT_HALO_AUTO_WIDTH (8) IFF compute_size > 1 AND no march-in
      !!     exclusion is active; else 0 (serial, or an incompatible feature).
      !!   requested >= 0 (user set it explicitly) => returned UNCHANGED; the
      !!     fail-loud exclusion checks in validate_config police an explicit
      !!     bt_halo > 0 against an incompatible feature (user asked for the
      !!     impossible), so the auto-resolution never overrides an explicit 0
      !!     or an explicit width.
      integer, intent(in) :: requested
         !! The namelist value: BT_HALO_AUTO_SENTINEL (-1) for auto, else >= 0.
      integer, intent(in) :: compute_size
         !! Number of compute ranks (1 = serial).
      logical, intent(in) :: exclusion_active
         !! .true. iff any march-in exclusion feature is on (see
         !! `bt_halo_auto_exclusion`).  Only consulted for the auto sentinel.
      integer :: width
      if (requested >= 0) then
         width = requested
      else if (compute_size > 1 .and. .not. exclusion_active) then
         width = BT_HALO_AUTO_WIDTH
      else
         width = 0
      end if
   end function resolve_bt_halo

   pure subroutine bt_halo_auto_exclusion(cfg, excluded, reason)
      !! Single source of truth for the BT march-in exclusion set (mirrors the
      !! explicit-`bt_halo > 0` fail-loud checks in `validate_config`).  Reports
      !! whether ANY exclusion is active and names the first one (for the AUTO
      !! resolution log).  When none is active, `reason` is "".
      type(config_t), intent(in) :: cfg
      logical, intent(out) :: excluded
      character(len=:), allocatable, intent(out) :: reason
      reason = ""
      if (cfg%ocean%wetdry%enable) then
         reason = "wetdry enable"
      else if (cfg%ocean%bt%use_cont_type) then
         reason = "use_cont_type"
      else if (cfg%ocean%bt%upstream_h_face) then
         reason = "upstream_h_face"
      else if (cfg%ocean%tides%enable) then
         reason = "tides enable"
      else if (cfg%ocean%psurf%enable) then
         ! Surface-pressure loading rides the SAME `eta_forcing` seam as
         ! the body tide, and the wide BT clone carries no copy of it:
         ! `run_stage_split` `error stop`s on `bt_halo > 0` with an
         ! `eta_forcing` actual.  Must exclude alongside tides.
         reason = "psurf enable"
      else if (cfg%ocean%porous%enable) then
         reason = "porous enable"
      else if (cfg%ocean%cavity_dyn%enable) then
         ! The wide BT clone rebuilds `metrics_w` from the grid formula and
         ! carries no `z_draft`, so its reference depth would be the BED —
         ! it would solve a 1000 m ocean where the cavity has 500 m of
         ! water under 500 m of ice.  `validate_config` refuses an EXPLICIT
         ! `bt_halo > 0`; AUTO must resolve to 0 here or a multi-rank
         ! cavity run manufactures a width the user never asked for and
         ! then trips that abort.
         reason = "cavity_dyn enable"
      else if (trim(cfg%ocean%grid%grid_config) == "supergrid" .or. &
               trim(cfg%ocean%grid%grid_config) == "tripolar") then
         reason = "grid_config='"//trim(cfg%ocean%grid%grid_config)//"'"
      end if
      excluded = len_trim(reason) > 0
   end subroutine bt_halo_auto_exclusion

   pure function p_top_has_producer(cfg) result(has)
      !! Is there anything in this configuration that WRITES
      !! `multilayer_state_t%p_top`?
      !!
      !! `p_top = metrics%p_ice_ref + sf%p_surf` has exactly two
      !! producers, and both halves count:
      !!
      !!   * `&ocean_psurf_nml enable` — the atmospheric surface-pressure
      !!     seam, which fills `sf%p_surf` and refreshes `p_top` once per
      !!     outer step in `ocean_dyn_step_split`;
      !!   * `&ocean_cavity_dyn_nml enable` — the STATIC ice-shelf load
      !!     `p_ice_ref = rho_ref*GRAVITY*z_draft`, assembled into
      !!     `p_top` by `configure_ocean_cavity`.  Static is not the same
      !!     as absent: the draft never changes, so the configure-time
      !!     assembly is the final value and there is nothing to refresh.
      !!
      !! Used by the `&ocean_pgf_nml p_top_in_bc` inert-knob warning,
      !! which must name both — a cavity run is precisely the case where
      !! `p_top_in_bc` is not merely live but REQUIRED (for a draft that
      !! varies), so warning that it is inert there was telling the user
      !! the opposite of the truth.
      type(config_t), intent(in) :: cfg
      logical :: has
      has = cfg%ocean%psurf%enable .or. cfg%ocean%cavity_dyn%enable
   end function p_top_has_producer

   pure function substep_drag_ignores_bdrag_form(cfg) result(ignores)
      !! Is `&ocean_bt_nml substep_drag` blind to the configured bottom
      !! drag?
      !!
      !! `compute_bt_rem` builds the barotropic substep damping
      !! `Htot/(Htot + r·hbbl·dt_inner)` from the LINEAR coefficient
      !! `&ocean_bdrag_nml r` alone.  Under any `form` other than
      !! `"linear"` the slow bottom drag is not that operator, so the
      !! knob either does nothing at all (`r = 0`, the default — the
      !! factor is identically 1) or damps the barotropic mode with a
      !! linear drag the slow step never applies.  Neither is what the
      !! user asked for; the configure WARNING that names both knobs
      !! reads this predicate.  A warning, not a refusal: the knob is
      !! harmless-if-useless in the default case, and a deliberate
      !! BT-only linear sponge is a legitimate (if unusual) request.
      type(config_t), intent(in) :: cfg
      logical :: ignores
      ignores = cfg%ocean%bt%substep_drag .and. &
                trim(adjustl(cfg%ocean%bdrag%form)) /= "linear"
   end function substep_drag_ignores_bdrag_form

   pure function cavity_draft_is_uniform(cfg) result(uniform)
      !! Is the configured ice-shelf draft UNIFORM over the whole array?
      !!
      !! The one predicate, so the two rules that turn on it cannot drift
      !! apart: the `&ocean_pgf_nml p_top_in_bc` REFUSAL (a load with a
      !! gradient must have a consumer — a uniform load is bit-identically
      !! inert in the FV_MOM6 top BC, which is the theorem in
      !! `compute_fv_mom6_impl`'s docstring), and the `z_fixed` x cavity
      !! staircase WARNING (only a uniform draft is validated on a
      !! quasi-geopotential coordinate).
      !!
      !! `draft_config = "none"` is uniform because it is identically
      !! zero.  `"flat"` is uniform only when NO box bound clips it —
      !! a clipped flat draft has a calving front, which is a step, and a
      !! step is the largest gradient in the domain.  `"linear"` and
      !! `"file"` are never assumed uniform: this is a NAMELIST-level
      !! predicate and cannot see the filled array (the configure-time
      !! twin in `configure_ocean_cavity` tests `maxval /= minval` on the
      !! field itself, which is the stricter check and runs later).
      type(config_t), intent(in) :: cfg
      logical :: uniform
      character(len=:), allocatable :: dcfg
      dcfg = trim(adjustl(cfg%ocean%cavity_dyn%draft_config))
      uniform = (dcfg == "none")
      if (dcfg == "flat") then
         uniform = abs(cfg%ocean%cavity_dyn%draft_x0) >= 1.0e29_wp .and. &
                   abs(cfg%ocean%cavity_dyn%draft_x1) >= 1.0e29_wp .and. &
                   abs(cfg%ocean%cavity_dyn%draft_y0) >= 1.0e29_wp .and. &
                   abs(cfg%ocean%cavity_dyn%draft_y1) >= 1.0e29_wp
      end if
   end function cavity_draft_is_uniform

   subroutine warn_unknown_bc(bc_str, param_name)
      !! Warn if a BC string does not match any known type
      character(len=*), intent(in) :: bc_str
      character(len=*), intent(in) :: param_name

      select case (trim(bc_str))
      case ("wall", "open", "tidal", "nested", "inflow", &
            "discharge", "clamped", "sponge", "chapman")
         continue
      case default
         call logger%warning("Unknown boundary type '"//trim(bc_str)// &
                             "' for "//param_name//", defaulting to wall")
      end select

   end subroutine warn_unknown_bc

   pure logical function forcing_components_need_thermo(enable_components, &
                                                        enable_thermodynamics) result(bad)
      !! `&ocean_forcing_nml enable_components=.true.` with
      !! `&ocean_thermo_nml enable_thermodynamics=.false.` is a
      !! configuration that validates but does nothing (the component
      !! set only feeds the assembler->apply_tracers thermo path) —
      !! named predicate per the `lateral_closure_is_implemented` idiom
      !! so `validate_config` reads as one line.
      logical, intent(in) :: enable_components, enable_thermodynamics
      bad = enable_components .and. .not. enable_thermodynamics
   end function forcing_components_need_thermo

   pure logical function diag_density_levels_ok(vgrid, diags, n_rho_levels, rho_levels) result(ok)
      !! .true. iff a density-coordinate diagnostic selection (global
      !! `vgrid='density'`, or a per-diagnostic `:density`/`:rho` attribute
      !! anywhere in `diags`) has a usable `rho_levels` axis: at least one
      !! entry, strictly increasing (`invert_density_targets` assumes a
      !! monotone light->dense target list), and within the declared
      !! `MAX_OCEAN_DIAG_Z_LEVELS` array bound.  Unlike sigma/z*, density
      !! bins have no auto-fill, so an unset axis must abort rather than
      !! silently size to zero.  When density is not requested at all this
      !! is unconditionally `.true.` (no constraint).
      character(len=*), intent(in) :: vgrid
      character(len=*), intent(in) :: diags
      integer, intent(in) :: n_rho_levels
      real(wp), intent(in) :: rho_levels(:)
      logical :: density_requested
      integer :: k

      density_requested = (trim(vgrid) == "density") .or. &
                          (index(diags, ":density") > 0) .or. &
                          (index(diags, ":rho") > 0)
      if (.not. density_requested) then
         ok = .true.
         return
      end if

      ok = .true.
      if (n_rho_levels <= 0) then
         ok = .false.
         return
      end if
      if (n_rho_levels > size(rho_levels)) then
         ok = .false.
         return
      end if
      do k = 2, n_rho_levels
         if (rho_levels(k) <= rho_levels(k - 1)) then
            ok = .false.
            return
         end if
      end do
   end function diag_density_levels_ok
   pure logical function sponge_source_is_implemented(tag) result(ok)
      !! .true. iff `&ocean_sponge_nml damp_source` names a source that is
      !! actually filled. `"file"` is a recognised name (PR-23b, needs the
      !! PR-14 reader) but has no kernel yet — a source with no filler
      !! must abort, not silently build an all-zero map (the
      !! `lateral_closure_is_implemented` idiom).
      !!
      !! Kept under its original name because it is the DAMP axis; the
      !! TARGET axis has its own predicate below. They were one function
      !! until `"linear_z"` arrived, at which point a shared list would
      !! have accepted `damp_source="linear_z"` — a spelling with no
      !! meaning on that axis — as valid.
      character(len=*), intent(in) :: tag
      select case (trim(tag))
      case ("band")
         ok = .true.
      case default
         ok = .false.
      end select
   end function sponge_source_is_implemented

   pure logical function sponge_target_is_implemented(tag) result(ok)
      !! .true. iff `&ocean_sponge_nml target_source` names a reference
      !! state that is actually filled: `"ic"` (snapshot of the seeded
      !! initial condition) or `"linear_z"` (analytic affine geopotential
      !! profile, re-evaluated on the live layer geometry). `"file"` is
      !! recognised but has no reader yet (PR-23b).
      character(len=*), intent(in) :: tag
      select case (trim(tag))
      case ("ic", "linear_z")
         ok = .true.
      case default
         ok = .false.
      end select
   end function sponge_target_is_implemented

   pure logical function sponge_ramp_is_valid(tag) result(ok)
      !! .true. iff `&ocean_sponge_nml ramp` names a band shape the
      !! `damp_source="band"` filler implements.
      character(len=*), intent(in) :: tag
      select case (trim(tag))
      case ("cosine", "linear")
         ok = .true.
      case default
         ok = .false.
      end select
   end function sponge_ramp_is_valid

   ! ==================================================================
   ! Strict-schema construction (rdb_nml_schema).  Lives here (not in
   ! rdb_config_schema) because read_config builds + parses the schema
   ! itself — a module-level `use rdb_config_schema` from rdb_config
   ! would be a circular dependency.  rdb_config_schema re-exports
   ! build_rdb_schema so the production entry points keep their
   ! original import.
   ! ==================================================================

   subroutine build_rdb_schema(cfg, schema)
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

   subroutine register_sim(cfg, schema)
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

   subroutine register_grid(cfg, schema)
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

   subroutine register_time(cfg, schema)
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
      call schema%add_group(g)
   end subroutine register_time

   subroutine register_mpi(cfg, schema)
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

   subroutine register_logging(cfg, schema)
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

   subroutine register_vcoord(cfg, schema)
      !! `&vcoord_nml`: vertical-coordinate + ALE-remap controls.
      !! `vcoord_type` / `remap_method` / `zstar_stretching` are enums
      !! with the canonical sets accepted by the rdb_vcoord parsers.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      integer, pointer :: pi
      real(wp), pointer :: pr
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
                             "z_fixed partial steps: close every face whose layer is "// &
                             "an inert filler on either side (z-level wall, free-slip)"))
      call schema%add_group(g)
   end subroutine register_vcoord

   subroutine register_kappa_shear(cfg, schema)
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

   subroutine register_slopes(cfg, schema)
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

   subroutine register_gm(cfg, schema)
      !! `&ocean_gm` (Gent-McWilliams thickness diffusion, capability [2]).
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      logical, pointer :: pl
      real(wp), pointer :: pr

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

      call schema%add_group(g)
   end subroutine register_gm

   subroutine register_redi(cfg, schema)
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

      call schema%add_group(g)
   end subroutine register_redi

   subroutine register_varmix(cfg, schema)
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
      call g%add(nml_int("kh_res_fn_power", pi, "Resolution-function power p (even)"))
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

   subroutine register_meke(cfg, schema)
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

   subroutine register_tidal_mixing(cfg, schema)
      !! `&ocean_tidal_mixing` (St-Laurent/Simmons internal-tide mixing).
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      logical, pointer :: pl
      real(wp), pointer :: pr

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

      call schema%add_group(g)
   end subroutine register_tidal_mixing

   subroutine register_conv(cfg, schema)
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

   subroutine register_porous(cfg, schema)
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

   subroutine register_ddiff(cfg, schema)
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

   subroutine register_ocean_tides(cfg, schema)
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

   subroutine register_ocean_psurf(cfg, schema)
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

   subroutine register_ocean_cavity_dyn(cfg, schema)
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
                             "fv_mom6 PGF, sigma/zstar only)"))
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

      call schema%add_group(g)
   end subroutine register_ocean_cavity_dyn

   subroutine register_ocean_cavity_melt(cfg, schema)
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

   subroutine register_epbl(cfg, schema)
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
      call g%add(nml_enum("combine", ps, "Combine vs interior closure kv/kt", &
                          allowed=[character(len=3) :: "add", "max"]))
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

   subroutine register_wavespeed(cfg, schema)
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
                          "DEFERRED N2-monotonising depth (EBT path); < 0 = off"))
      pl => cfg%ocean%wavespeed%use_ebt
      call g%add(nml_logical("use_ebt", pl, &
                             "DEFERRED equivalent-barotropic variant"))
      pi => cfg%ocean%wavespeed%n_wavespeed
      call g%add(nml_int("n_wavespeed", pi, "Recompute cadence (every N steps)", min=1))

      call schema%add_group(g)
   end subroutine register_wavespeed

   subroutine register_foxkemper(cfg, schema)
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

   subroutine register_physics(cfg, schema)
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

   subroutine register_output(cfg, schema)
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

   subroutine register_boundary(cfg, schema)
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

   subroutine register_tracer(cfg, schema)
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

   subroutine register_initial_condition(cfg, schema)
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

   subroutine register_nonhydrostatic(cfg, schema)
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

   subroutine register_ocean_grid(cfg, schema)
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

   subroutine register_ocean_coriolis(cfg, schema)
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
                          "or weno3/weno5/weno7 (WENO-Z); weno5/weno7 need nghost>=3/4", &
                          allowed=[character(len=8) :: "centered", "weno3", &
                                   "weno5", "weno7"]))
      pl => cfg%ocean%coriolis%use_state_fluxes
      call g%add(nml_logical("use_state_fluxes", pl, &
                             "mom6-corrector CorAdv consumes continuity's renormalised "// &
                             "mass fluxes (MOM6 mass-consistent uh/vh)"))
      pl => cfg%ocean%coriolis%bound_coriolis
      call g%add(nml_logical("bound_coriolis", pl, &
                             "clamp the energy-scheme Coriolis accel to the (f+zeta)*v "// &
                             "velocity-form range (MOM6 BOUND_CORIOLIS; sadourny_energy only)"))
      ps => cfg%ocean%coriolis%corner_h
      call g%add(nml_enum("corner_h", ps, "PV corner-thickness construction (energy scheme)", &
                          allowed=[character(len=9) :: "cell_mean", "mom6_area"]))
      call schema%add_group(g)
   end subroutine register_ocean_coriolis

   subroutine register_ocean_thermo(cfg, schema)
      !! `&ocean_thermo_nml`: thermodynamics switch + scalar surface fluxes.
      type(config_t), target, intent(in) :: cfg
      type(nml_schema_t), intent(inout) :: schema
      type(nml_group_t) :: g
      logical, pointer :: pl
      real(wp), pointer :: pr
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
      call schema%add_group(g)
   end subroutine register_ocean_thermo

   subroutine register_ocean_forcing(cfg, schema)
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

   subroutine register_ocean_ice(cfg, schema)
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

   subroutine register_ocean_ice_ic(cfg, schema)
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

   subroutine register_ocean_restore(cfg, schema)
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

   subroutine register_ocean_geothermal(cfg, schema)
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

   subroutine register_ocean_sponge(cfg, schema)
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

   subroutine register_ocean_tracers(cfg, schema)
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

   subroutine register_ocean_bt(cfg, schema)
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
                             "h-weight the post-substep BT corrector (MOM6 frhatu)"))
      pl => cfg%ocean%bt%correction_visc_rem
      call g%add(nml_logical("correction_visc_rem", pl, &
                             "h*visc_rem joint corrector weight (requires "// &
                             "correction_h_weighted; visc_rem is produced by vdiff and "// &
                             "is inert, =1, without ocean_vdiff_nml implicit_drag)"))
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
                             "(MOM6 u_cor = u + du*visc_rem; requires correction_visc_rem)"))
      pl => cfg%ocean%bt%forcing_visc_rem
      call g%add(nml_logical("forcing_visc_rem", pl, &
                             "MOM6 wt_u parity: h*visc_rem-weight the BT forcing "// &
                             "depth-mean so friction-damped (glued) layers do not "// &
                             "force the fast loop (requires correction_visc_rem)"))
      pl => cfg%ocean%bt%correction_bc_pgf
      call g%add(nml_logical("correction_bc_pgf", pl, &
                             "Per-layer baroclinic-PGF retro-correction for the eta change"))
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
                         "(-1 = auto: 8 under multi-rank if compatible, else 0; "// &
                         "0 = explicit off, bit-identical)", &
                         min=BT_HALO_AUTO_SENTINEL))
      call schema%add_group(g)
   end subroutine register_ocean_bt

   subroutine register_ocean_debug(cfg, schema)
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

   subroutine register_ocean_mpi(cfg, schema)
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

   subroutine register_ocean_wetdry(cfg, schema)
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

   subroutine register_ocean_pgf(cfg, schema)
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
      pl => cfg%ocean%pgf%p_top_in_bc
      call g%add(nml_logical("p_top_in_bc", pl, &
                             "FV_MOM6: add the top-of-column load ms%p_top to "// &
                             "the pressure-stack surface boundary condition "// &
                             "pa(nz+1) = rho_ref*g*eta + p_top"))
      call schema%add_group(g)
   end subroutine register_ocean_pgf

   subroutine register_ocean_eos(cfg, schema)
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
                                   "roquet_spv", "teos10"]))
      ps => cfg%ocean%eos%tfreeze_set
      call g%add(nml_enum("tfreeze_set", ps, &
                          "Named liquidus coefficient set for eos_freezing_point "// &
                          "(T_f = l1*S + l2 + l3*p): 'seaice' = SIS2/MOM6 "// &
                          "(-0.054, 0, -7.53e-8), 'isomip' = ISOMIP+ "// &
                          "(-0.0573, 0.0832, -7.53e-8)", &
                          allowed=[character(len=6) :: "seaice", "isomip"]))
      pr => cfg%ocean%eos%p_ref
      call g%add(nml_real("p_ref", pr, &
                          "Reference pressure for the potential density "// &
                          "ms%rho_layer (horizontally uniform by design; "// &
                          "0 => surface density)", &
                          units="Pa", min=0.0_wp))
      call schema%add_group(g)
   end subroutine register_ocean_eos

   subroutine register_ocean_bdrag(cfg, schema)
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

   subroutine register_ocean_tdrag(cfg, schema)
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

   subroutine register_ocean_hdiff(cfg, schema)
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

   subroutine register_ocean_hvisc(cfg, schema)
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

   subroutine register_ocean_vmix(cfg, schema)
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
      call g%add(nml_logical("use_closure", pl, "Master switch for the interior closure (PP81)"))
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
      pi => cfg%ocean%vmix%dt_tracer_advect_ratio
      call g%add(nml_int("dt_tracer_advect_ratio", pi, &
                         "Horizontal tracer advect runs at ratio*dt_dyn (MOM6 DT_TRACER_ADVECT)", min=1))
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
      call g%add(nml_real("bkgnd_prandtl", pr, "Background Prandtl number Kv_bg=prandtl*Kd_bg", &
                          min=0.0_wp))
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

   subroutine register_ocean_vdiff(cfg, schema)
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
                             "MOM6 parity: attenuate the slow explicit accelerations by "// &
                             "the per-layer viscous remnant (u = u0 + visc_rem*(u-u0) "// &
                             "after the applies); requires ocean_bt_nml correction_visc_rem"))
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
      pr => cfg%ocean%vdiff%hbbl_visc
      call g%add(nml_real("hbbl_visc", pr, &
                          "Bottom-layer scale for the hvel_mom6 botfn blend (MOM6 HBBL)"))
      pl => cfg%ocean%vdiff%bbl_glue
      call g%add(nml_logical("bbl_glue", pl, &
                             "MOM6 bottomdraglaw coupling parity: kv_bbl botfn glue at "// &
                             "near-bed interfaces + piston bed drag. Absorbs the spurious "// &
                             "grounded-layer PGF as MOM6 does (PGF_BUG.md par.9). Requires "// &
                             "hvel_mom6 + implicit_drag + linear bottom drag"))
      pr => cfg%ocean%vdiff%bbl_piston
      call g%add(nml_real("bbl_piston", pr, &
                          "BBL drag piston velocity u* for bbl_glue (MOM6 CDRAG*DRAG_BG_VEL)", &
                          units="m/s", min=0.0_wp))
      pl => cfg%ocean%vdiff%hvel_upwind
      call g%add(nml_logical("hvel_upwind", pl, &
                             "Near-bed upwind blend in the hvel_mom6 thickness build "// &
                             "(.false. = pure harmonic; the u-sign blend flip-flops on "// &
                             "roundoff at rest and collapses the BBL glue, PGF_BUG.md par.9.8)"))
      call schema%add_group(g)
   end subroutine register_ocean_vdiff

   subroutine register_ocean_continuity(cfg, schema)
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

   subroutine register_ocean_isopycnal(cfg, schema)
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
                             "overlap across the face (VCOORD_LAGRANGIAN only, mont / "// &
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

   subroutine register_ocean_topo(cfg, schema)
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

   subroutine register_ocean_ic(cfg, schema)
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

   subroutine register_ocean_zinit(cfg, schema)
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

   subroutine register_ocean_data(cfg, schema)
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

   subroutine register_ocean_dataovr(cfg, schema)
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

   subroutine add_dataovr_entry(g, e, tag, what, units)
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

   subroutine register_ocean_diag(cfg, schema)
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
                             "Order-invariant EFP console totals: rank-count reproducible, "// &
                             "exact drift residual (default off = byte-identical)"))
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

   subroutine register_ocean_bc(cfg, schema)
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

end module rdb_config
