!! Runtime configuration via Fortran namelists
module rdb_config
   !! Reads simulation parameters from a namelist input file
#ifdef LFORTRAN_PASSING
   use rdb_constants, only: wp, DEG2RAD, SEAWATER_CP_DEFAULT, &
                            nz_stack_required, nz_stack_is_sufficient
#else
   use rdb_constants, only: NZ_STACK_MAX, wp, DEG2RAD, SEAWATER_CP_DEFAULT, &
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
   use rdb_ocean_bulk_flux, only: bulk_flux_config_t, bulk_flux_config_n_iter_ok
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
   real(wp), parameter, public :: BT_DTBT_UNSET = -1.0e30_wp
      !! Sentinel for `&ocean_bt_nml dtbt` "not given".
   real(wp), parameter :: LEGACY_TRACER_S_REF = 0.0_wp
      !! Historical `&tracer_nml S_ref` default (PSU).
   real(wp), parameter :: LEGACY_TRACER_BETA_S = 0.78_wp
      !! Historical `&tracer_nml beta_S` default (kg/m^3 per PSU).
   real(wp), parameter :: LEGACY_TRACER_T_REF = 15.0_wp
      !! Historical `&tracer_nml T_ref` default (degC).
   real(wp), parameter :: LEGACY_TRACER_ALPHA_T = 0.17_wp
      !! Historical `&tracer_nml alpha_T` default (kg/m^3 per degC).

   public :: MAX_SW_BANDS
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
             ocean_atm_forcing_config_t, &
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
   public :: bt_dtbt_is_set, bt_split_derived
   public :: bt_halo_auto_exclusion
   public :: p_top_has_producer
   public :: substep_drag_ignores_bdrag_form
   public :: bbl_glue_is_effective
   public :: atm_forcing_prereq_code, atm_forcing_prereq_message
   public :: ATM_PREREQ_OK, ATM_PREREQ_START_DATE, ATM_PREREQ_COMPONENTS, &
             ATM_PREREQ_CARTESIAN, ATM_PREREQ_FILELIST, ATM_PREREQ_SW_SOURCE, &
             ATM_PREREQ_DATAOVR, ATM_PREREQ_BULK_FLUX, ATM_PREREQ_N_ITER, &
             ATM_PREREQ_ORB_YEAR
   integer, parameter :: ATM_PREREQ_OK = 0
      !! `atm_forcing_prereq_code`: every prerequisite holds.
   integer, parameter :: ATM_PREREQ_START_DATE = 1
      !! `&time_nml start_date` is blank.
   integer, parameter :: ATM_PREREQ_COMPONENTS = 2
      !! `&ocean_forcing_nml enable_components` is off.
   integer, parameter :: ATM_PREREQ_CARTESIAN = 3
      !! `grid_config = "cartesian"` (no geographic lon/lat).
   integer, parameter :: ATM_PREREQ_FILELIST = 4
      !! A `filelist_<field>` path is blank.
   integer, parameter :: ATM_PREREQ_SW_SOURCE = 5
      !! `sw_pen_frac > 0` without `sw_source = "q_sw"`.
   integer, parameter :: ATM_PREREQ_DATAOVR = 6
      !! A `&ocean_dataovr_nml` wind/heat/evap/lprec file is set.
   integer, parameter :: ATM_PREREQ_BULK_FLUX = 7
      !! `&ocean_bulk_flux_nml enable` is not set.
   integer, parameter :: ATM_PREREQ_N_ITER = 8
      !! `&ocean_bulk_flux_nml n_iter < 1`.
   integer, parameter :: ATM_PREREQ_ORB_YEAR = 9
      !! `orb_year` is not the default (not wired into the coszen kernel).
   public :: ocean_bt_forcing_visc_rem_on
   public :: ocean_bt_renorm_visc_rem_on
   public :: ocean_bt_rem_from_visc_rem_on
   public :: ocean_bt_visc_rem_producer_on
   public :: cavity_draft_is_uniform
   public :: zfixed_cavity_nu_h_below_envelope
   public :: ZFIXED_CAVITY_NU_H_MIN
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
      !! driver, where `compute_size` is known).  AUTO resolves to 0 — the
      !! march-in is OPT-IN (set `bt_halo` explicitly) because it is not yet
      !! bit-reproducible against the serial run.  See `resolve_bt_halo`.
   integer, parameter :: BT_HALO_AUTO_WIDTH = 8
      !! The recommended EXPLICIT march-in width (the validated production
      !! width).  AUTO no longer selects it — see `resolve_bt_halo`.

   real(wp), parameter :: ZFIXED_CAVITY_NU_H_MIN = 2.0_wp
      !! Lower envelope of `&ocean_hvisc_nml nu_h` (m²/s) under
      !! `vcoord_type = "z_fixed"` x `&ocean_cavity_dyn_nml enable`.  Below
      !! it `validate_config` WARNS (see
      !! `zfixed_cavity_nu_h_below_envelope` for the measurement).

   integer, parameter :: MAX_TIDAL_CONSTITUENTS = 10
      !! Maximum number of tidal constituents
   integer, parameter :: MAX_OCEAN_DIAG_Z_LEVELS = 128
      !! Maximum number of z-levels for ocean-diag z_fixed output vgrid.
      !! [C1] Raised 64 -> 128 for ACCESS-OM3 parity: the 75-level OM3
      !! z* grid needs every model interface available as a z_fixed/
      !! sigma/zstar/density diag level.  Sizes four fixed `config_t`
      !! arrays only (`z_levels`, `sigma_levels`, `zstar_levels`,
      !! `rho_levels`, all host-side namelist storage, never mapped to
      !! device); the diag manager's actual output buffers are
      !! allocated to the USED count (`n_z_levels`, ...), not this bound,
      !! and are independently capped by `NZ_STACK_MAX` (128 by default,
      !! `rdb_constants`) which already covered 128 before this change.
      !! Existing namelists are unaffected: `n_z_levels` (etc.) they set
      !! stays <= 64 <= 128, so this is bit-identical for them.
   integer, parameter :: MAX_OCEAN_LAYER_RHO_INIT = 64
      !! Maximum number of per-layer density init entries.
   integer, parameter :: MAX_Z_FIXED_DZ = 128
      !! Maximum number of `&vcoord_nml z_fixed_dz` entries (nominal
      !! `z_fixed` layer thicknesses).
   integer, parameter :: MAX_SW_BANDS = 6
      !! Maximum number of penetrating-shortwave bands under
      !! `&ocean_thermo_nml sw_scheme = "manizza05"`.
   real(wp), parameter :: CFL_BT_SAFETY_DEFAULT = 0.65_wp
      !! Default `&ocean_bt_nml cfl_bt_safety`; `validate_config` compares against it.
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
         !! Coriolis acceleration into the range of the four
         !! neighbouring VELOCITY-form estimates `(f+ζ)·v` BEFORE the KE-gradient
         !! subtraction.  Kills the thin-layer PV blow-up: `q=(f+ζ)/h_corner`
         !! is huge at a vanishing edge and `q·vh` can inject an acceleration
         !! the local velocity field cannot support; the `(f+ζ)·v` bound caps
         !! it.  `form="sadourny_energy"`, or `form="sadourny"` (enstrophy) with
         !! `pv_adv_scheme="centered"`; `sadourny_hk` and the WENO PV variants
         !! fail loud.  MOM6's double-gyre control run enables it.  Default
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
      character(len=16) :: sw_scheme = "default"
         !! Penetrating-shortwave scheme.  `"default"` (bit-identical): the
         !! two-band Paulson & Simpson (1977) decay above (`sw_pen_frac`,
         !! `sw_band_ratio`, `sw_zeta1/2`).  `"manizza05"`: the
         !! chlorophyll-dependent N-band scheme of Manizza et al. (2005,
         !! GRL 32, L05603) -- the WHOLE shortwave is split into bands at
         !! the surface and each band decays as `exp(-kappa_n z)` with a
         !! per-column `kappa_n(chl)`.  Requires `sw_chl_file`; `sw_pen_frac`
         !! must stay 0 (or 1) because the scheme has no separate
         !! penetrating fraction; `validate_config` refuses anything else.
         !! The deposition kernel, the KPP `B_0` coupling and the EPBL SW
         !! TKE ledger all consume the same per-column band opacities.
      integer :: sw_nbands = 3
         !! Number of bands under `"manizza05"` (3 .. MAX_SW_BANDS).  Band
         !! 1 = blue, 2 = red, the near-infrared share is split evenly over
         !! bands 3..n.  Bands beyond 3 need their own `sw_opacity_coef` /
         !! `sw_chl_power` entries.
      real(wp) :: sw_blue_frac = 0.5_wp
         !! Blue share of the visible band (MOM6 `BLUE_FRAC_SW`).  The
         !! visible band is `0.42` of the total shortwave.
      character(len=256) :: sw_chl_file = ""
         !! NetCDF climatology of surface chlorophyll (mg m-3) on the model
         !! horizontal grid, read monthly-cyclic over a 365-day year whose
         !! day 0 is model `t = 0`.  Required under `"manizza05"`.
      character(len=64) :: sw_chl_var = "CHL"
         !! Variable name of the chlorophyll field in `sw_chl_file`.
      real(wp) :: sw_opacity_coef(2*MAX_SW_BANDS) = [0.0232_wp, 0.074_wp, 0.225_wp, &
                                                     0.037_wp, 2.86_wp, 0.0_wp, &
                                                     0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                                                     0.0_wp, 0.0_wp]
         !! Manizza (2005) opacity coefficients `(a_n, b_n)` pairs, m-1:
         !! `kappa_n = a_n + b_n * chl**p_n`.  Defaults are the three-band
         !! values (MOM6 `OPACITY_VALUES_MANIZZA`).
      real(wp) :: sw_chl_power(MAX_SW_BANDS) = [0.674_wp, 0.629_wp, 0.0_wp, &
                                                0.0_wp, 0.0_wp, 0.0_wp]
         !! Chlorophyll exponents `p_n` (MOM6 `CHOROPHYLL_POWER_MANIZZA`).
      real(wp) :: sw_opacity_land = 10.0_wp
         !! Opacity assigned to land columns, m-1 (`OPACITY_LAND_VALUE`).
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
         !! Scalar target SST (degC), in the MODEL convention: potential
         !! temperature under `ts_convention = "pot_prac"`, Conservative
         !! Temperature under `"cons_abs"`.  Never converted.
      real(wp) :: restore_sss = 0.0_wp
         !! Scalar target SSS (PSU), in the MODEL convention: practical
         !! salinity under `ts_convention = "pot_prac"`, Absolute Salinity
         !! under `"cons_abs"` (OM3 restores toward a file variable that is
         !! already Absolute Salinity).  Never converted.
      character(len=256) :: sss_file = ""
         !! File-based SSS restoring (OM3 `SALT_RESTORE_FILE`): NetCDF
         !! `field(x, y, time)` on the model grid, in the convention named by
         !! `sss_convention`.  Blank (default) keeps the scalar target
         !! `restore_sss` and the content-injection path, bit-identical.
         !! Requires `enable_restore_salt` and `piston_s > 0`.  The restoring
         !! acts as a SALT FLUX `F = rho0 * w * dS` (`w = piston_s/86400`,
         !! `dS` clamped to `max_delta_sss`) in the `Q_salt` unit
         !! convention (S in g/kg times kg/m^2/s), deposited in the top layer
         !! and mirrored into the salt budget contributor.
      character(len=64) :: sss_var = "salt"
         !! Variable name in `sss_file` (OM3 `SALT_RESTORE_VARIABLE`).
      logical :: sss_time_cyclic = .true.
         !! `.true.` ⇒ the file is a climatology wrapped with period
         !! `sss_cycle_period` (reader `DATA_TIME_CYCLIC`); `.false.` ⇒
         !! linear in time, out-of-range model time is fatal.
      real(wp) :: sss_cycle_period = 31536000.0_wp
         !! Climatology period (s) for `sss_time_cyclic` (default 365 days).
      real(wp) :: max_delta_sss = 999.0_wp
         !! Cap on `|S_target - S|` (OM3 `MAX_DELTA_SRESTORE`, g/kg).  The
         !! default 999 is "no clamp".  Needs `sss_file`.
      logical :: sss_adjust_net_zero = .false.
         !! Remove the global net salt input (OM3
         !! `ADJUST_NET_SRESTORE_TO_ZERO`).  Needs `sss_file`.
      logical :: sss_adjust_by_scaling = .false.
         !! With `sss_adjust_net_zero`: remove the net by SCALING the larger
         !! of the positive / negative area integrals (OM3
         !! `ADJUST_NET_SRESTORE_BY_SCALING`; zero contour unmoved).
         !! `.false.` ⇒ subtract the area-mean flux over wet cells instead.
      logical :: sss_as_sflux = .false.
         !! OM3 `SRESTORE_AS_SFLUX`.  The only restoring form this model has
         !! is the salt-flux form (no virtual-precipitation form), so with
         !! `sss_file` this MUST be `.true.` (fail-loud otherwise); without a
         !! file the scalar path is already a salt-flux form and either
         !! value is accepted.
      character(len=16) :: sss_convention = "model"
         !! Salinity convention of the FILE: `"model"` (default; used as-is,
         !! OM3's `asalt` under `cons_abs`) or `"practical"` (WOA-style SP,
         !! converted per cell with `sr_from_sp` under `cons_abs`; a no-op
         !! under `pot_prac`).
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
         !! CFL, evaluated per WET cell (local depth with local cell size,
         !! MOM6 `set_dtbt`; land never limits it).  Default `.false.`.
      real(wp) :: cfl_bt_safety = CFL_BT_SAFETY_DEFAULT
         !! Safety fraction on the shallow-water CFL bound when
         !! `auto_n_inner = .true.` (typical 0.65-0.7).
      real(wp) :: dtbt = BT_DTBT_UNSET
         !! OM3-style barotropic step definition (MOM6 `DTBT`).  Unset
         !! (`BT_DTBT_UNSET`, the default) leaves `n_inner` to `n_inner` /
         !! `auto_n_inner` + `cfl_bt_safety`.  When set it REPLACES them:
         !! `< 0` is the fraction `|dtbt|` of the largest stable step
         !! `dt_max` (domain minimum over wet cells of
         !! `1/dt_max^2 = 1/2 (1 + 2 bebt) [ (1/A_T) sum_faces g D_f dy_f/dx_f
         !! + sum_corners f^2 ]`), `0` means `-0.98`, and `> 0` is a fixed
         !! substep in seconds.  `n_inner = ceiling(dt / dtbt_used)` in every
         !! case.  Fails loud together with `auto_n_inner = .true.` or a
         !! non-default `cfl_bt_safety`.
      real(wp) :: dtbt_reset_period = -1.0_wp
         !! Seconds between recomputations of `dt_max` (and so `n_inner`)
         !! from the live layer thicknesses, when `dtbt` is set: `0`
         !! recomputes every dynamics step, `> 0` every that many seconds
         !! (a mid-period restart re-derives at its first step, so it is not
         !! bit-identical to the uninterrupted run), `< 0` (default) only at
         !! initialisation.  Ignored (fails loud if `>= 0`) when `dtbt` is
         !! unset.
      real(wp) :: bebt = 0.1_wp
         !! BT continuity-flux velocity projection weight (MOM6 `BEBT`).
         !! Default `0.1` = MOM6's own default for that parameter: damps
         !! the barotropic gravity waves at
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
         !! RETIRED (2026-10-02) — setting it `.true.` is a fail-loud
         !! `validate_config` error.  It distributed the barotropic
         !! increment by `h_face(k)/⟨h⟩_h`, which is not energy-conserving:
         !! beyond the barotropic `ΔKE` it adds a positive-definite source
         !! `½Δ²·H·(κ−1)`, `κ = Σh³Σh/(Σh²)²`, plus a shear feedback that
         !! grew the stretched-`z_fixed` 1-degree Southern Ocean to a
         !! non-finite state on day 16.  MOM6 has no h-weighted fold (its
         !! barotropic acceleration is the same in every layer).  The key
         !! stays registered only so the refusal can say why; drag-aware
         !! weighting is `correction_visc_rem`.
      logical :: correction_visc_rem = .false.
         !! RETIRED (2026-10, D1 follow-up) — setting `.true.` is a
         !! fail-loud `validate_config` error.  This used to weight the
         !! BT-corrector fold by `visc_rem(k)/⟨visc_rem⟩_h` instead of
         !! uniformly, but MOM6's barotropic solver gives every layer the
         !! SAME `u_accel_bt` via `accel_layer_u` — no `visc_rem` weight —
         !! folded into `up` BEFORE
         !! `vertvisc` distributes it via the SAME implicit friction the
         !! BBL glue uses, so MOM6 never damps it twice.  This fold did,
         !! and the SECOND, unbounded `visc_rem_k/⟨visc_rem⟩_h` ratio is
         !! what NaNs the 1-degree Southern Ocean z* open-step case under
         !! `bbl_glue` at step ~40 (see `validate_config`'s refusal
         !! message for the measured isolation). `visc_rem_chain` is the
         !! replacement and uses the UNIFORM fold, matching MOM6. The key
         !! stays registered only so the refusal can say why; the
         !! underlying kernel dispatch (`apply_bt_correction`'s
         !! `use_visc_rem`) and its direct unit tests are untouched.
      logical :: visc_rem_chain = .true.
         !! PR-3 (visc_rem audit + unification, D1 — revised 2026-10 once
         !! MOM6 settled the BT-correction fold question): ONE switch for
         !! exactly MOM6's `vertvisc_remnant`/`av_rem`/`bt_rem` set —
         !! DEFAULT ON since PR-4 (the flip, 2026-10), TOGETHER with
         !! `&ocean_vdiff_nml hvel_mom6`/`bbl_glue` — MOM6 gives the
         !! barotropic solver the same friction the layered vertical
         !! solve applies; without it the BBL glue's bed piston never
         !! reaches the fast mode (NaN at day 253 on the 1-degree
         !! Southern Ocean closed-face case under glue-only). `.false.`
         !! restores the pre-PR-4 linear-piston `substep_drag` behaviour
         !! for `bt_rem` (mutually exclusive with this switch, D2 below) —
         !! equivalent to switching on the visc_rem PRODUCER (decoupled
         !! from the retired `correction_visc_rem`, runs whenever any
         !! consumer below needs it) plus `forcing_visc_rem` (MOM6
         !! `wt_u`), `renorm_visc_rem` (MOM6 continuity `u_cor = u +
         !! du*visc_rem`, which IS how MOM6 ties `visc_rem` to a velocity
         !! correction) and `bt_rem_from_visc_rem` (MOM6's barotropic
         !! `av_rem`/`bt_rem` viscous-remnant depth mean) all at once.  The
         !! barotropic-correction fold (`apply_bt_correction`) stays
         !! UNIFORM under this switch — MOM6's `accel_layer_u` never
         !! weights it by `visc_rem`, so neither does this chain; see
         !! `correction_visc_rem`'s docstring for why that fold is
         !! retired, not folded in here.  Same caveat as before: without
         !! `&ocean_vdiff_nml implicit_drag=.true.` or `bbl_glue=.true.`,
         !! `visc_rem` stays ≡ 1 and the whole chain is a legal, warned,
         !! no-op.  The three `*_visc_rem` consumer knobs
         !! (`forcing_visc_rem`/`renorm_visc_rem`/`bt_rem_from_visc_rem`)
         !! stay individually settable (never retired) for the existing
         !! fine-grained tests, each an equivalent SUBSET of this switch,
         !! never a superset, and each now SELF-SUFFICIENT (no longer
         !! "requires" a separate producer knob — the producer runs
         !! whenever any one of them is on).  `strong_drag`/
         !! `rescale_strong_drag` stay separate namelist keys per D1
         !! (MOM6 has its own `BT_STRONG_DRAG`/`RESCALE_STRONG_DRAG`
         !! parameters for them) and require this switch (or
         !! `bt_rem_from_visc_rem`) on, same as before. `accel_visc_rem`
         !! (`&ocean_vdiff_nml`) is likewise NOT part of this chain and is
         !! RETIRED (see its own docstring) — no MOM6 state-update
         !! equivalent. Default `.false.` ⇒ no answer change (PR-4 is the
         !! default flip).  D2: mutually exclusive with `substep_drag`
         !! (checked via `bt_rem_from_visc_rem`'s existing requirement,
         !! which this switch satisfies identically to setting it
         !! directly).
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
         !! on (the flux expression tree gains the γ factors).
         !! Self-sufficient (D1 follow-up) — the visc_rem producer runs
         !! whenever this is on, decoupled from the retired
         !! `correction_visc_rem`.
      logical :: forcing_visc_rem = .false.
         !! MOM6 `wt_u` parity for the BT FORCING assembly:
         !! weight each layer's
         !! contribution to `F_bt_u/v` (and to the PGF-projection
         !! subtraction, which must use the same weights) by
         !! `h_face·visc_rem(k)` instead of `h_face`.  Layers the implicit
         !! friction will immediately damp — grounded sliver stacks under
         !! the vdiff BBL glue — then contribute nothing to the fast loop,
         !! which otherwise integrates the spurious grounded-layer PGF's
         !! depth-mean ballistically (PGF_BUG.md §9).
         !! Self-sufficient (D1 follow-up) — the visc_rem producer runs
         !! whenever this is on, decoupled from the retired
         !! `correction_visc_rem`.
      logical :: bt_rem_from_visc_rem = .false.
         !! PR-2 (bt-rem-from-av-rem): `bt_rem_u/v` built from the SAME
         !! viscous remnant the layered momentum solve uses, instead of
         !! the linear-piston `substep_drag` law or the static `1.0`
         !! no-op. Matches MOM6's barotropic solver:
         !! `av_rem = Σ_k frhat_k·visc_rem_k` (the visc_rem depth mean,
         !! `frhat_k` the face layer fraction `derive_bt_from_layers`
         !! already uses), then `bt_rem = mask·av_rem**(1/n_inner)`
         !! (zero where `mask·av_rem <= 0`), built ONCE per barotropic
         !! call — after the visc_rem producer, before the inner
         !! substeps — into the existing `bt_work%bt_rem_u/v` multiplier
         !! the substep loop already reads.  This is the fix for the
         !! bbl_glue day-253 instability (the barotropic solver was
         !! seeing weak explicit drag while the layers were strongly
         !! glued): the fast mode now feels the SAME friction.  Default
         !! `.false.` ⇒ no answer change.  Self-sufficient (D1
         !! follow-up) — the visc_rem producer runs whenever this is on,
         !! decoupled from the retired `correction_visc_rem`; mutually
         !! exclusive with
         !! `substep_drag` (bed drag would be double-counted — once
         !! inside visc_rem via the glue/implicit_drag fold, once again
         !! via the linear piston) and with `bt_halo > 0` (the wide-halo
         !! BT clone's `metrics_w`/halo-widened arrays carry no
         !! `av_rem`/`visc_rem` ghost width yet — same posture as
         !! porous). Composes with `wave_drag` (multiplied in after).
      logical :: strong_drag = .false.
         !! MOM6 `BT_STRONG_DRAG` (default `.false.`):
         !! replace the plain power form with the rational approximation
         !! `bt_rem = mask·n_inner·av_rem/(1 + (n_inner-1)·av_rem)`, which
         !! damps LESS aggressively per substep for a given `av_rem` —
         !! recommended only if the plain `av_rem**(1/n_inner)` form still
         !! leaves the acceptance-gate run unstable (D3, a measured
         !! deviation to report, not a default).  Requires
         !! `bt_rem_from_visc_rem = .true.`; inert otherwise (checked in
         !! `validate_config`).
      logical :: rescale_strong_drag = .false.
         !! MOM6 `RESCALE_STRONG_DRAG`: under `strong_drag`,
         !! `bt_rem**n_inner /= av_rem` exactly (the rational form is only
         !! an approximation), so the barotropic-correction `Δu` is
         !! rescaled by `min(bt_rem**n_inner/av_rem, 1.0)` before it is
         !! folded into the layers, keeping the correction consistent
         !! with the TRUE depth-mean remnant.  Requires `strong_drag =
         !! .true.`; inert (and refused) otherwise.
      character(len=16) :: frhat_scheme = "hybrid"
         !! `rdb_barotropic_coupling::frhat_h_face_step`'s per-layer
         !! face-thickness closure for EVERY barotropic depth mean that
         !! reads a layer's face thickness: `derive_bt_from_layers`,
         !! `face_depth_mean_u/v`, `face_depth_mean_rem_u/v`
         !! (`forcing_visc_rem`'s `wt_u`), `apply_bt_correction`'s
         !! open/visc_rem folds, and — through `face_depth_mean_u/v` —
         !! `compute_bt_rem_from_visc_rem`'s `av_rem`.
         !!
         !! `"hybrid"` — DEFAULT since PR-4 (the flip, 2026-10), TOGETHER
         !! with `visc_rem_chain`/`hvel_mom6`/`bbl_glue`: MOM6's `btcalc`
         !! `HVEL_SCHEME=HYBRID` (the MOM6 default) — arithmetic mean above the
         !! shallower of the two abutting columns' bed depths, harmonic
         !! mean below it (with a linear blend across the transition), so
         !! a thin partial-bed layer's weight is suppressed next to a
         !! thick abutting layer instead of inflated — closing the
         !! `docs/visc_rem_bt_rem_plan.md` Section 7 gap now that `av_rem`
         !! feeds the barotropic fast loop by default too. Does NOT change
         !! `derive_bt_from_layers`'s `use_upstream_h_face` branch (a
         !! roundabout-only alternative with no MOM6 frhat counterpart) or
         !! `apply_bt_correction`'s uniform fold (no depth mean involved
         !! there at all — matches MOM6 `accel_layer_u`).
         !!
         !! `"arithmetic"` — the pre-port plain two-abutting-cell mean,
         !! `h_face(k) = 0.5*(h_L(k)+h_R(k))`. `.false.` on all four of
         !! `hvel_mom6`/`bbl_glue`/`visc_rem_chain`/`frhat_scheme` together
         !! restores the pre-flip behaviour bit for bit.
      logical :: correction_bc_pgf = .false.
         !! Adds a per-layer baroclinic-PGF retro-correction for the η
         !! change during the BT substep.  Requires `pgf%form = "fv_mom6"`.
      logical :: bc_pgf_forcing = .true.
         !! MOM6 split (`BT_force` + `eta_PF`): the barotropic substep is
         !! forced by the depth mean of the FULL slow PGF — baroclinic
         !! part included — minus only the free-surface term that PGF
         !! itself carries at the η it was evaluated on (none for the
         !! surface-relative MONT/FV_LITE/FV_WRIGHT forms), which the
         !! substep's own `-g_bt·∇η` replaces.  `.false.` restores the legacy split, which
         !! subtracted the WHOLE depth-mean PGF and so never let the
         !! barotropic mode feel the baroclinic bottom-pressure gradient
         !! (JEBAR; ~0 Sv through Drake Passage on the global 1° case
         !! against MOM6's ~160 Sv).
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
         !! Wide-halo march-in width.  `-1` (default) = AUTO, which resolves
         !! to 0: the march-in is OPT-IN.  A decomposed run with the march-in
         !! is not bit-identical to the serial run over variable bathymetry or
         !! with open boundaries (measured by `test_ocean_decomp_bitid_mpi`;
         !! flat-bottom closed / periodic / spherical cases are), so the
         !! default keeps every decomposition bit-reproducible.  `8`
         !! (`BT_HALO_AUTO_WIDTH`) is the recommended explicit width.
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
         !!
         !! Cost: 5 EOS evaluations per layer + 15 per face per layer (T and
         !! S vary through the layer, so no closed form and no hoisting),
         !! each inlined into a kernel that runs one GPU thread per CELL
         !! (edges, layer and face integrals alike).  Global 1-degree PPM,
         !! 5 days, one V100: `ocean_pgf` 3.5 s under `wright` (was 13.7 s,
         !! then 6.0 s), 6.6 s under `roquet_spv` (was 14.6 s, then
         !! 12.4 s) — against 1.9 / 2.5 s for the PCM in-situ default.
      integer :: recon_scheme = 1
         !! In-layer reconstruction scheme: 1 = PLM, 2 = PPM.  Mirrors
         !! MOM6 `Recon_Scheme`.  Only consulted when
         !! `reconstruct_for_pressure = .true.`.
      logical :: insitu_density = .true.
         !! FV_MOM6 constant-by-layer (PCM) density at its IN-SITU pressure
         !! `p = -g*rho0*z` (MOM6 `int_density_dz_generic_pcm`; MOM6
         !! parity, default).  `.false.` = the legacy integral of the
         !! POTENTIAL density `ms%rho_layer` at the uniform
         !! `&ocean_eos_nml p_ref`, which loses the pressure dependence of
         !! the horizontal density gradient away from `p_ref` (on the global
         !! 1-degree spin-up: Drake Passage ~80 Sv against MOM6's ~155 Sv).
         !! Consulted only by `form = "fv_mom6"` with
         !! `reconstruct_for_pressure = .false.` and a PRESSURE-DEPENDENT
         !! EOS (`wright`, `roquet_spv`); for `linear` in-situ and
         !! potential density coincide and the legacy path runs,
         !! bit-identical.
         !!
         !! Cost (global 1-degree, 5 days, one V100; `.false.`: `ocean_pgf`
         !! 1.09 s, time loop 54.2 s):
         !!   * `wright` — the layer integral is ANALYTIC (MOM6
         !!     `int_density_dz_wright`): one polynomial evaluation per layer
         !!     and per cross-face sub-column.  `ocean_pgf` 1.86 s, time loop
         !!     55.0 s (+1.5 %).
         !!   * `roquet_spv` — no closed form (the integrand `1/SV(p)` is
         !!     rational in depth); the 5-point Boole rule of MOM6
         !!     `int_density_dz_generic_pcm`, with the (T, S) part of the EOS
         !!     evaluated ONCE per sub-column (`roquet_pcm_dpa_intz`, the SpV
         !!     value inlined from `rdb_roquet_spv.inc`) and only the
         !!     pressure Horner per Boole point.  `ocean_pgf` 2.55 s (was
         !!     10.73 s, then 3.64 s with the (T, S) part an out-of-line
         !!     call), time loop 55.8 s (+3 %, was 63.9 s): 1.36x Wright.
         !! Both run their integrals one GPU thread per CELL, not per column.
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
         !! Under the MOM6 split (`&ocean_bt_nml bc_pgf_forcing`,
         !! default) the depth mean of the layer PGF forces the barotropic
         !! mode, so the `p_surf` part of `p_top` is shed from that
         !! forcing as `g*grad(eta_ib)` and the `eta_forcing` seam carries
         !! it once; the static `p_ice_ref` part cancels inside `pa(nz+1)`
         !! against the datum-shifted `eta_geo`.  (The legacy split
         !! subtracted the whole depth mean, so there the uniform piece
         !! cancelled identically.)  On the UNSPLIT driver (`n_inner = 0`)
         !! there is no seam, so this term is the load's ONLY path into
         !! the momentum.
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
         !! 2015 specific-volume polynomial, polyTEOS10_75t, ΔS=24),
         !! "roquet_rho" (Roquet et al. 2015 IN-SITU DENSITY polynomial,
         !! polyTEOS10_55t, ΔS=32 — an independent fit of the same paper;
         !! returns ρ directly with no 1/SV division, so it is cheaper
         !! per point, but is NOT bit-identical to "roquet_spv" — the two
         !! agree only to the paper's stated rms accuracy, ~1e-3 kg/m³
         !! class).  Both "roquet_spv" and "roquet_rho" are incompatible
         !! with the "fv_wright" PGF (fails loud at configure).  "teos10"
         !! parses but is refused at configure (not implemented; the two
         !! Roquet branches above ARE the TEOS-10 polynomial fits).
         !!
         !! Per-EOS cost.  One point evaluation: linear ~5 flops; Wright
         !! ~20 flops + 1 division; Roquet ~130 flops, 2 sqrt, 2 divisions
         !! (value only — `roquet_spv_value`; the derivatives cost as much
         !! again and are computed only where a consumer uses them).  Where
         !! it shows (global 1-degree, 5 days, one V100, FV-MOM6 in-situ
         !! default): `ocean_pgf` 1.09 s linear, 1.88 s Wright (analytic),
         !! 2.55 s Roquet (factored Boole); time loop 54.3 / 55.1 / 55.8 s.
         !! On the CPU (gfortran) the gap is wider: `rho_layer` is ~7x
         !! dearer under Roquet than Wright, and `benchmark_ale` (Roquet,
         !! 8 steps) spends 4.9 s in `ocean_pgf` (33.1 s with the generic
         !! Boole rule, 6.5 s factored but unvectorised).
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
      character(len=16) :: ts_convention = "pot_prac"
         !! Meaning of the prognostic T/S tracers.
         !!
         !!   * `"pot_prac"` (DEFAULT ⇒ bit-identical) — potential
         !!     temperature (degC) and practical salinity; the Roquet et al.
         !!     (2015) variants convert to SR/CT internally.
         !!   * `"cons_abs"` — Conservative Temperature and Absolute
         !!     Salinity (the TEOS-10 pair, IOC/SCOR/IAPSO 2010).  The
         !!     Roquet points take T as CT and S as SA directly (no SR
         !!     factor, unit Jacobian); the heat-to-temperature conversion
         !!     still uses `c_p` (it does NOT switch to TEOS-10 cp0 by
         !!     itself -- OM3 runs CT/SA with `C_P = 3992`); the tracers
         !!     are labelled `sea_water_conservative_temperature` /
         !!     `sea_water_absolute_salinity`.
         !!
         !! `cons_abs` REQUIRES `eos` in {roquet_rho, roquet_spv}.  It is
         !! supported with `&ocean_ice_nml enable`: the ice-ocean boundary
         !! converts once (theta = pt_from_ct(SA, CT), SP = sp_from_sr(SA),
         !! theta freezing point), so the ice stays theta/SP while the ocean
         !! heat budget stays CT/SA.  It is refused (fail-loud) with
         !! `&ocean_cavity_melt_nml enable` (its liquidus is undecided).
         !!
         !! The other consumers convert or are in the model convention:
         !!   * a file `&ocean_zinit_nml` holding theta / SP
         !!     (`input_convention = "pot_prac"`) is converted at read time,
         !!     SA first then CT from SA;
         !!   * `&ocean_bulk_flux_nml` evaluates its formulae on
         !!     `theta = pt_from_ct(SA, CT)` of the surface layer (the only
         !!     surface-side conversion; heat contents of the precipitation
         !!     streams keep the model T with `c_p`);
         !!   * the `&ocean_restore_nml` `restore_sst` / `restore_sss`
         !!     targets are read IN the model convention (CT / SA; OM3
         !!     restores toward a file variable that is already Absolute
         !!     Salinity, and does not restore temperature);
         !!   * the derived diagnostics `temperature_pt` / `salinity_sp`
         !!     give theta / SP (CMOR `thetao` / `so`); `temperature` /
         !!     `salinity` are CT / SA (`bigthetao` / `absso`).
         !! Scalar/analytic initial conditions are interpreted IN the model
         !! convention and are NOT converted.
      character(len=16) :: tfreeze_form = "linear"
         !! Liquidus functional form for `eos_freezing_point`.
         !!
         !!   * `"linear"` (DEFAULT ⇒ bit-identical) — the `tfreeze_set`
         !!     coefficients.
         !!   * `"teos_poly"` — the TEOS-10 polynomial `CT_f(SA, p,
         !!     sat_frac)` (manual section 3.33; McDougall & Barker 2011).
         !!     S is read as SA, the result is a Conservative Temperature
         !!     (about -1.92 degC at SA = 35.16504, p = 0).  Requires
         !!     `ts_convention = "cons_abs"`.
      real(wp) :: tfreeze_sat_frac = 0.0_wp
         !! Air-saturation fraction of the TEOS-10 polynomial liquidus
         !! (`tfreeze_form = "teos_poly"`; 0 = no dissolved air, 1 = fully
         !! air-saturated).  Ignored by the linear form.  The default 0 is
         !! ACCESS-OM3 parity: OM3 runs `TFREEZE_FORM = "TEOS_POLY"`, the
         !! polynomial without the dissolved-air term.
      real(wp) :: c_p = SEAWATER_CP_DEFAULT
         !! Seawater heat capacity (J/(kg K)) used for every ocean-side
         !! heat-flux <-> temperature-tendency conversion (MOM6 `C_P`).
         !! Default 3992 under BOTH `ts_convention`s: it is the historical
         !! value and the one ACCESS-OM3 sets (its CT/SA run keeps
         !! `C_P = 3992`).  The exact TEOS-10 value for Conservative
         !! Temperature is cp0 = 3991.86795711963 (IOC/SCOR/IAPSO 2010,
         !! eq. 3.3.3); set it explicitly to get Q/(rho0*cp0) as the CT
         !! tendency.  Sea-ice thermodynamics keep their own constant.
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
      real(wp) :: kh_sin_lat = 0.0_wp
         !! Latitude-dependent Laplacian viscosity floor amplitude (m²/s):
         !! the background Laplacian becomes
         !! `max(existing, kh_sin_lat·|sin φ|^kh_pwr_of_sine)` (MOM6
         !! `KH_SIN_LAT`).  φ is the face-centre latitude (mean of the two
         !! `geolatBu` corners of the face).  It is a FLOOR under the dynamic
         !! closures (max-combined into the per-face `ah_face_*`, like
         !! `kh_vel_scale_live`).  Default 0 ⇒ off, bit-identical.  Fails
         !! loud on a `cartesian` grid (`geolat ≡ 0`).
      real(wp) :: kh_pwr_of_sine = 4.0_wp
         !! Exponent of `|sin φ|` in the `kh_sin_lat` floor (MOM6
         !! `KH_PWR_OF_SINE`).  Inert while `kh_sin_lat = 0`.
      real(wp) :: ah_vel_scale = 0.0_wp
         !! Biharmonic velocity-scale floor (m/s): the background that the
         !! Smagorinsky biharmonic is `max`ed against becomes
         !! `max(nu_4_bg, ah_vel_scale·(Δ²)^(3/2))`, `Δ²` the harmonic mean
         !! of the squared T-cell spacings (MOM6 `AH_VEL_SCALE`).  Requires
         !! `smag_ah`.  Default 0 ⇒ off, bit-identical.
      logical :: bound_coriolis_biharm = .false.
         !! Coriolis bound on the Smagorinsky biharmonic (MOM6
         !! `BOUND_CORIOLIS_BIHARM`): `A = |S|·(C_b·Δ⁴ + Δ⁶·f_max·|S|/(5V²))`
         !! with `V = bound_coriolis_vel` and `f_max` the largest `|f|` of the
         !! T cell's four corners.  Requires `smag_ah`.  Default off.
      real(wp) :: bound_coriolis_vel = 6.0_wp
         !! Velocity scale `V` (m/s, > 0) of the `bound_coriolis_biharm`
         !! term (MOM6 `BOUND_CORIOLIS_VEL`, whose default is `MAXVEL`).
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
      character(len=8) :: interior_closure = "pp81"
         !! Interior mixing closure under `use_closure`: `"pp81"` (default;
         !! Pacanowski & Philander 1981 Richardson-number closure,
         !! bit-identical) or `"none"` (no shear-driven interior closure --
         !! the interior `kv`/`kt` are the background floor (scalar, Henyey
         !! et al. 1986 / Harrison & Hallberg 2008 latitude factor, or
         !! Bryan-Lewis) plus whatever ePBL, kappa-shear, tides, convection
         !! and double diffusion contribute).  `"none"` is how an
         !! ACCESS-OM3-like run (`USE_PP81 = False`) avoids double-counting
         !! shear mixing.  Requires `use_closure`.  `large94`/`cvmix` are
         !! reserved and fail loud.
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
      real(wp) :: bkgnd_prandtl = -1.0_wp
         !! Background Prandtl number (MOM6 `PRANDTL_BKGND`).  Default -1 =
         !! off (independent `kv_bg` floor, bit-identical).  `> 0`: the
         !! interior momentum floor is `bkgnd_prandtl * kt_background +
         !! kv_bg` with `kt_background` the tracer background in use
         !! (scalar `kt_bg`, or the Henyey `max(bkgnd_kd_min, kt_bg*L(phi))`).
         !! ACCESS-OM3 equivalent: `bkgnd_prandtl = 1`, `pp81_nu_bg = 0`.
         !! Under `bkgnd_profile` the floor is `bkgnd_prandtl * kd_bg` and an
         !! unset value resolves to 1 (historical behaviour).
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
         !! (the assembly floor) — see `vmix_seed_backgrounds` -- and so
         !! is the additive `kv_bg` term of the `bkgnd_prandtl` floor even
         !! when `interior_closure = "none"` (set 0 for OM3's `KV = 0`).
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
      !! The implicit folds default `.false.` (explicit wind-stress apply);
      !! MOM6's face treatment (`hvel_mom6` + `bbl_glue`) defaults ON.
      logical :: implicit_stress = .false.
         !! Fold the surface wind stress into the vdiff surface (k=nz) RHS
         !! row (Neumann top-BC) instead of the explicit pre-solve add.
      logical :: accel_visc_rem = .false.
         !! RETIRED (PR-3, 2026-10): setting `.true.` is a fail-loud
         !! `validate_config` error.  This attenuated the slow EXPLICIT
         !! accelerations by the per-layer viscous remnant —
         !! `u = u_entry + visc_rem·(u_applied − u_entry)` after the
         !! CorAdv/PGF/hvisc/drag applies, before the BT correction — but
         !! PR-3's audit found NO MOM6 state-update equivalent:
         !! `btstep_layer_accel` and the split-explicit RK2 corrector's
         !! `up`/`vp` update both apply the depth-mean barotropic acceleration
         !! `u_accel_bt` UNIFORMLY across every layer, with no `visc_rem`
         !! weight anywhere in that path.  The real MOM6 mechanisms that
         !! multiply a velocity correction by `visc_rem` are
         !! `&ocean_bt_nml renorm_visc_rem` (continuity `u_cor = u +
         !! du·visc_rem`) and `rescale_strong_drag`.  The underlying
         !! kernels (`accel_visc_rem_snapshot`/`accel_visc_rem_reweight`
         !! in `rdb_ocean_dyn.F90`) and their direct unit tests
         !! (`tests/test_ocean_accel_visc_rem.F90`) are untouched — only
         !! the configure-time path to reach them from a namelist is
         !! refused.
      logical :: implicit_drag = .false.
         !! Fold the bottom drag into the vdiff bed-row diagonal (stress
         !! bottom-BC; the row is the face's first LIVE layer `k_bot_u/v`,
         !! `k = 1` off `z_fixed`) instead of the explicit pre-solve add.  Mutually
         !! exclusive with `&ocean_bdrag_nml implicit` (split-apply) and,
         !! without `bbl_glue`, incompatible with HBBL-distributed drag
         !! (`hbbl > 0`); both fail loud at configure.  Under `bbl_glue`
         !! the glue's piston replaces the fold.
      logical :: hvel_mom6 = .true.
         !! MOM6 momentum face thickness `hvel` for the vertical-friction
         !! solve (`vertvisc_coef`), with MOM6's arithmetic `h_shear` and
         !! the height-above-bed stack the BBL glue reads; the branch is
         !! `hvel_harmonic`.  DEFAULT ON (2026-10, with `bbl_glue`: MOM6's
         !! face treatment everywhere).  `.false.` ⇒ the historical
         !! arithmetic `h_u` / `face_thick` pair.
      real(wp) :: hbbl_visc = 10.0_wp
         !! MOM6 `HBBL` (m) for the `hvel_mom6` botfn blend without the
         !! BBL glue (under the glue each face's `bbl_thick` normalises the
         !! height above the bed), and the glue's HBBL when
         !! `&ocean_bdrag_nml hbbl = 0`.
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
      logical :: hvel_harmonic = .false.
         !! Which MOM6 face-thickness branch `hvel_mom6` builds (MOM6
         !! `HARMONIC_VISC`).  `.false.` (MOM6's default): arithmetic face
         !! thickness, blended to harmonic near the bed for thin -> thick
         !! flow, height above the bed `max(zh, z_clear)` — every face layer
         !! below the shallower bed of a step sits inside the BBL.  `.true.`:
         !! the harmonic branch (harmonic thickness, upwind-arithmetic blend
         !! gated by `hvel_upwind`) — the historical `hvel_mom6`.
      logical :: bbl_glue = .true.
         !! DEFAULT ON (2026-10; `.false.` restores the explicit / folded
         !! bed drag).  MOM6 `BOTTOMDRAGLAW`: the bottom drag lives in the vertical
         !! viscosity.  MOM6 `set_viscous_BBL` computes, per face and once
         !! per outer step, the BBL viscosity `kv_bbl = sqrt(CDRAG)·u*·
         !! bbl_thick` (`u*` from the HBBL-mean speed under the quadratic
         !! law, `sqrt(CDRAG)·DRAG_BG_VEL` under the linear one) and the
         !! rotation/stratification-limited `bbl_thick` (KW99); the
         !! momentum solve then raises the interface viscosity by
         !! `(kv_bbl − KV)·botfn` near the bed and takes the bed row as the
         !! piston `kv_bbl/(min(hvel₁/2, bbl_thick))`, which REPLACES the
         !! explicit / folded bed drag.  The drag law and its parameters
         !! come from `&ocean_bdrag_nml` (`form`, `cd` or `r`, `hbbl` —
         !! falling back to `hbbl_visc` when 0 — `bg_vel`,
         !! `bbl_thick_min`); with no drag configured the glue is inert.
         !! Requires `hvel_mom6` (fail-loud at configure).
      real(wp) :: bbl_piston = 3.0e-4_wp
         !! Historical constant BBL piston velocity (m/s), `kv_bbl =
         !! bbl_piston·hbbl_visc`: read only by a hand-built vdiff slot
         !! without the per-face BBL (unit tests).  A configured run's glue
         !! takes its drag from `&ocean_bdrag_nml`.
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
      character(len=10) :: combine = "add"
         !! "add" / "max" vs the interior closure's kv/kt
         !! (MOM6 EPBL_IS_ADDITIVE), or "max_shear": kappa-shear is merged
         !! first and ePBL adds only what exceeds it, i.e.
         !! bg + max(ePBL, kappa-shear) + others (ACCESS-OM3,
         !! EPBL_IS_ADDITIVE=False).  Identical to "add" when kappa-shear
         !! is off; refused with use_kpp.
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
         !! N^2-monotonising depth.  `< 0` = off.  NOT IMPLEMENTED:
         !! `validate_config` refuses `>= 0` together with `use_ebt`.
      logical :: use_ebt = .false.
         !! Also compute the equivalent-barotropic speed `cg1_ebt` and
         !! vertical structure `ebt_struct` (`cg1` stays first-baroclinic).
         !! GM consumes the structure with `&ocean_gm_nml use_ebt_struct`.
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
      logical :: use_ebt_struct = .false.
         !! MOM6 `KHTH_USE_EBT_STRUCT`: multiply the 2-D face KhTh (base +
         !! Visbeck + MEKE, after the CFL clamp) by the equivalent-barotropic
         !! mode structure at each interface (Hallberg 2013; Killworth &
         !! Blundell 2003).  Requires `&ocean_wavespeed_nml enable` +
         !! `use_ebt` and GM enabled (configure-time errors).  Default off
         !! ⇒ bit-identical.
      logical :: use_fgnv = .false.
         !! Solve the Ferrari, Griffies, Nurser & Vallis (2010) vertical
         !! boundary-value problem for the GM streamfunction
         !! (MOM6 `KHTH_USE_FGNV_STREAMFUNCTION`).  Default off ⇒ the local
         !! `-KhTh*S` streamfunction, bit-identical.  Requires GM and
         !! `&ocean_wavespeed_nml enable` (the first-baroclinic `cg1`).
      real(wp) :: fgnv_filter_scale = 1.0_wp
         !! `FGNV_FILTER_SCALE`, > 0.  Scales the SMOOTHING TERM:
         !! `c^2 = fgnv_filter_scale * c_face^2`, and `psi_unlim` is amplified
         !! by `(1 + fgnv_filter_scale)` before the solve, so the
         !! effective filter speed is `sqrt(fgnv_filter_scale) * cg1`
         !! (OM3: 0.1 ⇒ `0.316 cg1`).  Inert unless `use_fgnv`.
      real(wp) :: fgnv_c_min = 0.0_wp
         !! `FGNV_C_MIN`, m/s, >= 0.  Accepted for parity; the reference
         !! formulation does not apply it, so it is INERT here (OM3 sets 0).
      real(wp) :: fgnv_strat_floor = 1.0e-15_wp
         !! `FGNV_STRAT_FLOOR`, > 0.  `N^2` floor as a fraction of the
         !! rotation rate: `N^2 >= (fgnv_strat_floor * Omega)^2`, keeping the
         !! tridiagonal system non-degenerate.
      character(len=16) :: gm_src_form = "local"
         !! Form of the PE-release source `gm_src` that feeds MEKE
         !! (inverse of MOM6 `MEKE_GM_SRC_ALT`).  `"local"` (default,
         !! bit-identical): `rho0/4 Sum KH S^2 N^2 h` from the unlimited local
         !! slope.  `"work"` (OM3, `MEKE_GM_SRC_ALT = False`): the rate of
         !! potential-energy release by the FINAL bolus transports (FGNV,
         !! EBT structure, limiters included), per face then averaged to
         !! cells; needs `&ocean_gm_nml enable`.
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
      logical :: check_diffusive_cfl = .false.
         !! Subcycle the neutral flux so the explicit-diffusion stability
         !! limit is respected: `n_itt = max(1, ceiling(max_i CFL_i - 4*epsilon))` with
         !! `CFL_i = 2*sum_f(KhTr_f*dt*dy_f/dx_f)/area_i` (the sum over the
         !! cell's four faces, after the `max_diffusive_cfl` cap).  The
         !! global max is taken over all ranks, so `n_itt` is
         !! decomposition-independent.  Default off.  The MOM6 name
         !! `CHECK_DIFFUSIVE_CFL`; the OM3 value is `True`.
      real(wp) :: max_diffusive_cfl = -1.0_wp
         !! `<= 0` ⇒ off.  Otherwise each face coefficient is capped at
         !! `(max_diffusive_cfl/8)*min(area_L, area_R)` so no cell exceeds
         !! this `CFL_i`; without `check_diffusive_cfl` it also sets
         !! `n_itt = ceiling(max_diffusive_cfl)`.  The MOM6 name
         !! `MAX_TR_DIFFUSION_CFL`; the OM3 value is 2.0.
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
         !! Resolution-function power p (integer >= 1).
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
      character(len=8) :: tideamp_source = "uniform"
         !! RMS barotropic tidal speed source: `"uniform"` (the scalar
         !! `utide`, bit-identical default) | `"file"` (the static 2-D
         !! `tideamp_var` of `tideamp_file`; needs `e_compute`).  File
         !! fields are masked and roughness-clamped ONCE at setup with the
         !! resting column depth D = sum of `h_layer`.
      character(len=256) :: tideamp_file = ""
         !! NetCDF file holding the tidal-amplitude field (`tideamp_source
         !! = "file"`).
      character(len=64) :: tideamp_var = "tideamp"
         !! Variable name in `tideamp_file`, m/s (OM3: `tideamp`).
      character(len=8) :: h2_source = "uniform"
         !! Sub-grid roughness variance source: `"uniform"` (the scalar
         !! `h2_rough`) | `"file"` (the static 2-D `h2_var` of `h2_file`;
         !! needs `e_compute`).
      character(len=256) :: h2_file = ""
         !! NetCDF file holding the roughness-variance field (`h2_source =
         !! "file"`).
      character(len=64) :: h2_var = "h2"
         !! Variable name in `h2_file`, m^2 (OM3: `h2`).
      logical :: n2_bot_rough_avg = .false.
         !! Near-bed stratification averaged over the roughness height
         !! `sqrt(<h^2>)` (Melet et al. 2013): the bed-most interfaces up
         !! to that height give one `N_b^2`, which also REPLACES the layer
         !! `N^2` inside the span.  Default off = the raw two-deepest-layer
         !! gradient (bit-identical).  OM3 parity value: `.true.`.
      character(len=8) :: energy_limiter = "clamp"
         !! Per-layer cap `kd_max` handling: `"clamp"` (default: clip the
         !! layer Kd, DISCARD the clipped power) | `"reinject"` (cap the
         !! layer TKE at `rho0*h*(N^2+Omega^2)*kd_max`, offer the excess to
         !! the layer above; energy-conserving).  OM3 parity value:
         !! `"reinject"` with `kd_max = 0.1`.
      character(len=8) :: bed_share = "discard"
         !! Fate of the bed layer's share of the dissipated power:
         !! `"discard"` (default: owned by the bottom boundary layer, not
         !! deposited) | `"lump"` (carried into the first interior layer).
         !! OM3 parity value: `"lump"`.
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
      logical :: trim_ic_for_p_surf = .false.
         !! Trim the INITIAL column under the ice so it is at rest (MOM6
         !! `TRIM_IC_FOR_P_SURF`, `trim_for_ice`).  The load
         !! `p_ice_ref = rho_ref*g*z_draft` is the displaced weight at the
         !! REFERENCE density; a stratified column's displaced water weighs
         !! `g*int_{-z_draft}^{0} rho dz`, and the difference is a depth-
         !! uniform bottom-pressure gradient the MOM6 barotropic split
         !! (`&ocean_bt_nml bc_pgf_forcing`) adjusts to.  With the knob on,
         !! the load is kept (the ice MASS is what is prescribed) and each
         !! loaded column's initial top is moved to the depth `s` where
         !! `g*int_{-s}^{0} rho dz = p_ice_ref`, i.e. an initial
         !! `eta = z_draft - s` (a few cm under ISOMIP+ COLD), with T/S
         !! then evaluated at the trimmed layer centres.  Closed form, exact at
         !! the discrete FV interfaces: requires `&ocean_eos_nml
         !! eos="linear"` and `&ocean_zinit_nml enable, source="linear"`
         !! (the analytic profile is what defines `rho` above the ice
         !! base); anything else fails loud.  Default off => bit-identical.
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
         !! in the two abutting columns do not overlap AND the layer is at the
         !! floor on at least one side (`<= 2·max(angstrom_h, H_VANISHED)`,
         !! `nonoverlap_vanish_tol_for`) — i.e. exactly where the layer has
         !! wedged out against the bed and there is no common depth to
         !! difference the pressure across.  Non-overlap ALONE is not grounding:
         !! a sigma-seeded stack over a bathymetric step puts layers that are
         !! massive on both sides in disjoint z-intervals, and zeroing their PGF
         !! while continuity keeps moving their mass breaks the PGF-work / PE
         !! exchange (the compat-matrix staircase grew to the CFL panic at step
         !! 178; `test_ocean_lagrangian_staircase`).  Overlapping (physically shared)
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
      character(len=16) :: input_convention = "pot_prac"
         !! What the `source = "file"` T/S variables hold.
         !!
         !!   * `"pot_prac"` (DEFAULT => bit-identical) -- potential
         !!     temperature (degC) and practical salinity, as World Ocean
         !!     Atlas `ptemp` / `salt` do.  Under `&ocean_eos_nml
         !!     ts_convention = "cons_abs"` they are converted on the host,
         !!     after the depth interpolation and before the device map:
         !!     SA = sr_from_sp(SP) FIRST, then CT = ct_from_pt(SA, theta)
         !!     with the already-converted SA (IOC/SCOR/IAPSO 2010;
         !!     McDougall & Barker 2011).  Reference Salinity is used as
         !!     Absolute Salinity (no delta-SA atlas), the ACCESS-OM3
         !!     behaviour.
         !!   * `"cons_abs"` -- the file already holds CT / SA; read as is.
         !!     Refused (fail-loud) under a `pot_prac` model: there is no
         !!     CT -> theta initial-condition path.
         !!
         !! Dry-column `land_fill_t/s` and the `source = "linear"` profile
         !! are always in the MODEL convention and are never converted.
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

   type :: ocean_atm_forcing_config_t
      !! `&ocean_atm_forcing_nml`: JRA55-do-style prescribed atmospheric
      !! state (PR-C2c), read through nine `register_2d_filelist`
      !! instances and driving the bulk-flux-computed ocean/ice surface
      !! fluxes.  Default `enable = .false.` registers nothing ⇒
      !! bit-identical.  Field tags throughout are the 9 CMOR/JRA55-do
      !! names: `tas, huss, uas, vas, psl` (instantaneous states),
      !! `rlds, rsds, prra, prsn` (interval-mean fluxes).
      logical :: enable = .false.
         !! Master switch.
      character(len=16) :: source = "jra55do"
         !! Provider tag (currently only `"jra55do"` is implemented).
      character(len=256) :: filelist_tas = ""
         !! Path to a `tools/jra_filelist.py`-style one-file-per-line
         !! list for `tas` (air temperature). Blank = fail loud when
         !! `enable = .true.`.
      character(len=256) :: filelist_huss = ""
         !! Filelist for `huss` (specific humidity).
      character(len=256) :: filelist_uas = ""
         !! Filelist for `uas` (eastward 10 m wind).
      character(len=256) :: filelist_vas = ""
         !! Filelist for `vas` (northward 10 m wind).
      character(len=256) :: filelist_psl = ""
         !! Filelist for `psl` (sea-level pressure).
      character(len=256) :: filelist_rlds = ""
         !! Filelist for `rlds` (downwelling longwave).
      character(len=256) :: filelist_rsds = ""
         !! Filelist for `rsds` (downwelling shortwave, interval mean).
      character(len=256) :: filelist_prra = ""
         !! Filelist for `prra` (liquid precipitation rate).
      character(len=256) :: filelist_prsn = ""
         !! Filelist for `prsn` (frozen precipitation rate).
      character(len=64) :: var_tas = "tas"
         !! Variable name inside the `tas` filelist's files.
      character(len=64) :: var_huss = "huss"
         !! Variable name inside the `huss` filelist's files.
      character(len=64) :: var_uas = "uas"
         !! Variable name inside the `uas` filelist's files.
      character(len=64) :: var_vas = "vas"
         !! Variable name inside the `vas` filelist's files.
      character(len=64) :: var_psl = "psl"
         !! Variable name inside the `psl` filelist's files.
      character(len=64) :: var_rlds = "rlds"
         !! Variable name inside the `rlds` filelist's files.
      character(len=64) :: var_rsds = "rsds"
         !! Variable name inside the `rsds` filelist's files.
      character(len=64) :: var_prra = "prra"
         !! Variable name inside the `prra` filelist's files.
      character(len=64) :: var_prsn = "prsn"
         !! Variable name inside the `prsn` filelist's files.
      character(len=16) :: time_mode = "interannual"
         !! `"interannual"` (default) or `"repeat_year"` (RYF), passed
         !! straight to `register_2d_filelist`.
      character(len=32) :: ryf_start_date = "1990-05-01"
         !! RYF window start date (`time_mode = "repeat_year"` only).
      real(wp) :: ryf_period_days = 365.0_wp
         !! RYF window length in days (`time_mode = "repeat_year"` only).
      character(len=16) :: states_regrid = "bilinear"
         !! Horizontal regrid for the instantaneous states (`tas, huss,
         !! uas, vas, psl`). Only `"bilinear"` is implemented; a
         !! patch-recovery option is not offered.
      character(len=16) :: fluxes_regrid = "conservative"
         !! Horizontal regrid for the interval-mean fluxes (`rlds, rsds,
         !! prra, prsn`): `"conservative"` (default) or `"bilinear"`.
      character(len=16) :: sw_interp = "coszen"
         !! `"coszen"` (default): hold the shortwave record constant
         !! over its interval and modulate by a diurnal cosine-zenith
         !! factor. `"linear"`: plain time-blend like every other flux
         !! (disables the hold + the coszen kernel).
      real(wp) :: sw_offset_s = -5400.0_wp
         !! Time offset (s) applied to the `rsds` query time before the
         !! bracket search -- moves an interval-centred timestamp onto
         !! the interval start for the hold (`sw_interp = "coszen"`).
      integer :: orb_year = 1990
         !! Orbital-parameter year, reserved for a future Berger (1978)
         !! upgrade. The v1 coszen kernel uses the fixed Spencer (1971)
         !! Fourier-fit coefficients, so a non-default value is refused
         !! at configure rather than silently ignored.
      logical :: wind_relative = .true.
         !! Subtract the ocean surface current from the atmospheric wind
         !! before computing the bulk-flux relative wind speed.
      real(wp) :: ocean_albedo = 0.066_wp
         !! Ocean surface shortwave albedo (Large & Yeager 2009).
      logical :: verbose = .false.
         !! Log every bracket advance of the 9 underlying readers.
   end type ocean_atm_forcing_config_t

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
      logical :: reproducing_sums = .true.
         !! PR-32: when `.true.` (the DEFAULT since v0.1.0), the ocean
         !! console-conservation totals (Mass/KE/Salt/Heat + sea-ice area)
         !! and the salt/heat closed-budget `out`/`src` terms use
         !! order-invariant extended-fixed-point (EFP) summation
         !! (`rdb_efp` + `halo_allreduce_efp_list`) instead of plain FP
         !! `!$acc parallel loop reduction(+:acc)` + `MPI_SUM` — the
         !! printed console is identical on every rank count (1 included)
         !! and reduction order, and the `Error` residual is formed via
         !! `efp_real_diff` (a fixed-point difference) rather than a
         !! double subtraction of two already-quantised totals.  It runs
         !! only at the status cadence and changes diagnostic TEXT only,
         !! never the trajectory.  `.false.` restores the pre-PR-32 FP
         !! console (whose last digits depend on the decomposition).  See
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
      type(ocean_atm_forcing_config_t) :: atm_forcing
         !! JRA55-do-style prescribed atmospheric state
         !! (`&ocean_atm_forcing_nml`, PR-C2c).  Default-OFF ⇒
         !! bit-identical.
      type(bulk_flux_config_t)        :: bulk_flux
         !! Bulk-flux kernel configuration (`&ocean_bulk_flux_nml`).
         !! Default-OFF ⇒ bit-identical. `enable` is REQUIRED (set it
         !! explicitly) whenever `atm_forcing%enable` is on --
         !! `validate_config` refuses atm forcing without it.
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
      character(len=32) :: start_date = ""
         !! `&time_nml`: proleptic-Gregorian calendar date the model
         !! clock's `t=0` corresponds to ("YYYY-MM-DD" or "YYYY-MM-DD
         !! hh:mm:ss"). Empty (the default) leaves `t` a plain relative
         !! clock -- bit-identical to every namelist shipped before this
         !! knob existed. Validated at `validate_config` (parsed via
         !! `rdb_calendar::parse_date`); consumed by the absolute-date
         !! forcing readers of later OM3 waves, not by anything in this
         !! wave.
      character(len=24) :: calendar = "gregorian"
         !! `&time_nml`: calendar backing `start_date`. Only proleptic
         !! Gregorian is implemented
         !! (`rdb_calendar::calendar_name_is_implemented` --
         !! gregorian/standard/proleptic_gregorian all mean the same day
         !! count here). `noleap`/`360_day` are not offered: the enum
         !! refuses them rather than silently aliasing them to Gregorian.
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
         !! surface interface for VCOORD_RHO / VCOORD_HYCOM.  With
         !! `rho_target_dense` builds the uniform light→dense
         !! target-density linspace (MOM6 `ALE_COORDINATE_CONFIG=UNIFORM`
         !! analogue) under `rho_target_profile = "uniform"`.
      real(wp) :: rho_target_dense = 1030.0_wp
         !! Densest target interface potential density (kg/m³) — the bed
         !! interface for VCOORD_RHO / VCOORD_HYCOM.  See `rho_target_light`.
      character(len=16) :: rho_target_profile = "uniform"
         !! Target-density profile of `vcoord_type = "rho" | "hycom"`:
         !! `"uniform"` (default — the `rho_target_light`→`rho_target_dense`
         !! linspace, byte-identical) or `"list"` (the `nz_layers+1`
         !! interface densities in `rho_target_list`, light first — MOM6
         !! target densities from a list/file, e.g. the
         !! `HYBRID:file,sigma2,dz` coordinate of OM4).  A uniform list is a
         !! poor fit wherever most of the volume sits in a narrow density
         !! range (the Southern Ocean: half its volume inside ~2 linspace
         !! layers).  Refused on any other coordinate.
      real(wp) :: rho_target_list(MAX_Z_FIXED_DZ + 1) = -1.0_wp
         !! `rho_target_profile = "list"`: interface potential densities
         !! (kg/m³) at `rho_ref_pressure`, LIGHTEST (surface) FIRST, strictly
         !! increasing.  Exactly `nz_layers+1` leading positive entries; the
         !! rest unset (`<= 0`, default `-1`).
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
         !! actually vanishes layers (`z_fixed`, `zstar`, `zstar_full`, wet/dry) —
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
         !! Accepted on `z_fixed`, `zstar` (MOM6 z*) and `zstar_full` —
         !! the three GEOMETRIC families whose bed-side fillers sit at
         !! fixed reference depths (under `zstar`: the `z_fixed` fillers,
         !! decided at `η = 0` and exactly static under the dilation;
         !! under `zstar_full`: every layer below a column's partial cell
         !! when the column is shallower than the `zstar_h_surf_target`
         !! fine zone).  The mask is built from the coordinate's target at
         !! `η = 0`, and under `zstar` / `zstar_full` the IC is then seeded
         !! on that target.  Refused on every other coordinate, and without
         !! fillers to close (`z_fixed` / `zstar` without a resolved
         !! `z_fixed_h_ref`; `zstar_full` with `zstar_h_surf_target <= 0`).
         !! The name is historical.
      character(len=16) :: z_fixed_profile = "uniform"
         !! Nominal layer-thickness profile of `vcoord_type = "z_fixed"`
         !! and of `vcoord_type = "zstar"` (MOM6 z*, which dilates this
         !! same nominal profile per column by `(H + η)/H`):
         !! `"uniform"` (default — `max_depth/nz_layers` everywhere,
         !! byte-identical), `"list"` (the thicknesses in `z_fixed_dz`,
         !! surface first) or `"tanh"` (a hyperbolic-tangent stretching
         !! from `z_fixed_dz_top` at the surface, scaled to sum to
         !! `&ocean_topo_nml max_depth`; see `rdb_vcoord ::
         !! z_fixed_nominal_dz`).  The z_fixed target builder, the
         !! closed-face mask, `k_top`, the cavity partial-top rule and the
         !! bed partial-cell rule all read the profile.
         !!
         !! It is the z* COORDINATE RESOLUTION, and `vcoord_type = "hycom"`
         !! reads it too, as its z* nominal floor (MOM6 HYCOM1 floors its
         !! interfaces at the `coordinateResolution` that
         !! `ALE_COORDINATE_CONFIG` / `HYBRID:` defines): interface `k` is
         !! kept at least `Σ dz·(H+η)/H` deep.  "uniform" there means
         !! `max_depth/nz_layers` METRES, not `1/nz_layers` of the column.
         !! The `z_fixed_*` names are kept for every reader rather than
         !! introducing a parallel table that could disagree.  `vcoord_type =
         !! "zstar"` (MOM6 z*) reads it as its nominal levels as well.  Refused
         !! on any other coordinate.
      real(wp) :: z_fixed_dz(MAX_Z_FIXED_DZ) = -1.0_wp
         !! `z_fixed_profile = "list"`: nominal layer thicknesses (m),
         !! SURFACE FIRST (MOM6 `ALE_COORDINATE_CONFIG = "PARAM:..."` /
         !! `vgrid` order).  Exactly `nz_layers` leading positive entries;
         !! the rest unset (`<= 0`, default `-1`).  The profile's total is
         !! its sum; a column deeper than that puts the excess in the bed
         !! layer.
      real(wp) :: z_fixed_dz_top = 2.0_wp
         !! `z_fixed_profile = "tanh"`: surface-layer nominal thickness (m).
         !! Must satisfy `nz_layers*z_fixed_dz_top < max_depth`.
      real(wp) :: z_fixed_tanh_center = 0.5_wp
         !! `z_fixed_profile = "tanh"`: transition centre as a fraction of
         !! the layer-index span, `[0, 1]` (0 = surface, 1 = bed).
      real(wp) :: z_fixed_tanh_width = 0.25_wp
         !! `z_fixed_profile = "tanh"`: transition width as a fraction of
         !! the layer-index span (`> 0`; small = abrupt).

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

   interface
      module subroutine read_config(filename, cfg, schema, ierr)
            !! See submodule [[rdb_config_read]].
         implicit none
         character(len=*), intent(in) :: filename
         type(config_t), target, intent(out) :: cfg
         type(nml_schema_t), intent(out), optional :: schema
         integer, intent(out), optional :: ierr
      end subroutine read_config

      module subroutine read_config_impl(filename, cfg, sch, ierr)
            !! See submodule [[rdb_config_read]].
         implicit none
         character(len=*), intent(in) :: filename
         type(config_t), target, intent(out) :: cfg
         type(nml_schema_t), intent(out) :: sch
         integer, intent(out), optional :: ierr
      end subroutine read_config_impl

      module subroutine read_config_from_string(text, cfg, schema, ierr)
            !! See submodule [[rdb_config_read]].
         implicit none
         character(len=*), intent(in) :: text
         type(config_t), target, intent(out) :: cfg
         type(nml_schema_t), intent(out), optional :: schema
         integer, intent(out), optional :: ierr
      end subroutine read_config_from_string

      module subroutine read_config_from_string_impl(text, cfg, sch, ierr)
            !! See submodule [[rdb_config_read]].
         implicit none
         character(len=*), intent(in) :: text
         type(config_t), target, intent(out) :: cfg
         type(nml_schema_t), intent(out) :: sch
         integer, intent(out), optional :: ierr
      end subroutine read_config_from_string_impl

      module subroutine split_to_lines(text, lines, n_lines)
            !! See submodule [[rdb_config_read]].
         implicit none
         character(len=*), intent(in) :: text
         character(len=:), allocatable, intent(out) :: lines(:)
         integer, intent(out) :: n_lines
      end subroutine split_to_lines

      module subroutine apply_time_unit_cascade(cfg)
            !! See submodule [[rdb_config_read]].
         implicit none
         type(config_t), intent(inout) :: cfg
      end subroutine apply_time_unit_cascade

      module subroutine apply_cartesian_degrees(cfg)
            !! See submodule [[rdb_config_read]].
         implicit none
         type(config_t), intent(inout) :: cfg
      end subroutine apply_cartesian_degrees

      module subroutine validate_config(cfg, ierr)
            !! See submodule [[rdb_config_validate]].
         implicit none
         type(config_t), intent(in) :: cfg
         integer, intent(out), optional :: ierr
      end subroutine validate_config

      module subroutine check_dataovr_file(e, tag, has_error)
            !! See submodule [[rdb_config_validate]].
         implicit none
         type(dataovr_entry_config_t), intent(in) :: e
         character(len=*), intent(in) :: tag
         logical, intent(inout) :: has_error
      end subroutine check_dataovr_file

      pure module function dataovr_entry_is_valid(e) result(ok)
            !! See submodule [[rdb_config_validate]].
         implicit none
         type(dataovr_entry_config_t), intent(in) :: e
         logical :: ok
      end function dataovr_entry_is_valid

      pure module function dataovr_time_is_valid(cfg) result(ok)
            !! See submodule [[rdb_config_validate]].
         implicit none
         type(ocean_dataovr_config_t), intent(in) :: cfg
         logical :: ok
      end function dataovr_time_is_valid

      pure module function dataovr_freshwater_needs_components(cfg, use_components) result(bad)
            !! See submodule [[rdb_config_validate]].
         implicit none
         type(ocean_dataovr_config_t), intent(in) :: cfg
         logical, intent(in) :: use_components
         logical :: bad
      end function dataovr_freshwater_needs_components

      pure module function atm_forcing_prereq_code(cfg) result(code)
            !! See submodule [[rdb_config_validate]].
         implicit none
         type(config_t), intent(in) :: cfg
         integer :: code
      end function atm_forcing_prereq_code

      pure module function atm_forcing_prereq_message(code) result(msg)
            !! See submodule [[rdb_config_validate]].
         implicit none
         integer, intent(in) :: code
         character(len=:), allocatable :: msg
      end function atm_forcing_prereq_message

      pure module function dataovr_any_tag_set(cfg) result(any_set)
            !! See submodule [[rdb_config_validate]].
         implicit none
         type(ocean_dataovr_config_t), intent(in) :: cfg
         logical :: any_set
      end function dataovr_any_tag_set

      pure module function ice_hlim_count(hlim) result(n)
            !! See submodule [[rdb_config_validate]].
         implicit none
         real(wp), intent(in) :: hlim(:)
         integer :: n
      end function ice_hlim_count

      pure module subroutine ice_hlim_spec_is_valid(hlim, ncat, ok, reason)
            !! See submodule [[rdb_config_validate]].
         implicit none
         real(wp), intent(in) :: hlim(:)
         integer, intent(in) :: ncat
         logical, intent(out) :: ok
         character(len=:), allocatable, intent(out) :: reason
      end subroutine ice_hlim_spec_is_valid

      pure module function resolve_bt_halo(requested, compute_size, exclusion_active) result(width)
            !! See submodule [[rdb_config_validate]].
         implicit none
         integer, intent(in) :: requested
         integer, intent(in) :: compute_size
         logical, intent(in) :: exclusion_active
         integer :: width
      end function resolve_bt_halo

      pure module subroutine bt_halo_auto_exclusion(cfg, excluded, reason)
            !! See submodule [[rdb_config_validate]].
         implicit none
         type(config_t), intent(in) :: cfg
         logical, intent(out) :: excluded
         character(len=:), allocatable, intent(out) :: reason
      end subroutine bt_halo_auto_exclusion

      pure module function p_top_has_producer(cfg) result(has)
            !! See submodule [[rdb_config_validate]].
         implicit none
         type(config_t), intent(in) :: cfg
         logical :: has
      end function p_top_has_producer

      pure module function bbl_glue_is_effective(cfg) result(on)
            !! See submodule [[rdb_config_validate]].
         implicit none
         type(config_t), intent(in) :: cfg
         logical :: on
      end function bbl_glue_is_effective

      pure module function substep_drag_ignores_bdrag_form(cfg) result(ignores)
            !! See submodule [[rdb_config_validate]].
         implicit none
         type(config_t), intent(in) :: cfg
         logical :: ignores
      end function substep_drag_ignores_bdrag_form

      pure module function ocean_bt_forcing_visc_rem_on(cfg) result(on)
            !! See submodule [[rdb_config_validate]].
         implicit none
         type(config_t), intent(in) :: cfg
         logical :: on
      end function ocean_bt_forcing_visc_rem_on

      pure module function ocean_bt_renorm_visc_rem_on(cfg) result(on)
            !! See submodule [[rdb_config_validate]].
         implicit none
         type(config_t), intent(in) :: cfg
         logical :: on
      end function ocean_bt_renorm_visc_rem_on

      pure module function ocean_bt_rem_from_visc_rem_on(cfg) result(on)
            !! See submodule [[rdb_config_validate]].
         implicit none
         type(config_t), intent(in) :: cfg
         logical :: on
      end function ocean_bt_rem_from_visc_rem_on

      pure module function ocean_bt_visc_rem_producer_on(cfg) result(on)
            !! See submodule [[rdb_config_validate]].
         implicit none
         type(config_t), intent(in) :: cfg
         logical :: on
      end function ocean_bt_visc_rem_producer_on

      pure module function cavity_draft_is_uniform(cfg) result(uniform)
            !! See submodule [[rdb_config_validate]].
         implicit none
         type(config_t), intent(in) :: cfg
         logical :: uniform
      end function cavity_draft_is_uniform

      pure module function zfixed_cavity_nu_h_below_envelope(cfg) result(below)
            !! See submodule [[rdb_config_validate]].
         implicit none
         type(config_t), intent(in) :: cfg
         logical :: below
      end function zfixed_cavity_nu_h_below_envelope

      module subroutine warn_unknown_bc(bc_str, param_name)
            !! See submodule [[rdb_config_validate]].
         implicit none
         character(len=*), intent(in) :: bc_str
         character(len=*), intent(in) :: param_name
      end subroutine warn_unknown_bc

      pure logical module function forcing_components_need_thermo(enable_components, &
                                                                  enable_thermodynamics) result(bad)
            !! See submodule [[rdb_config_validate]].
         implicit none
         logical, intent(in) :: enable_components, enable_thermodynamics
      end function forcing_components_need_thermo

      pure logical module function diag_density_levels_ok(vgrid, diags, n_rho_levels, rho_levels) result(ok)
            !! See submodule [[rdb_config_validate]].
         implicit none
         character(len=*), intent(in) :: vgrid
         character(len=*), intent(in) :: diags
         integer, intent(in) :: n_rho_levels
         real(wp), intent(in) :: rho_levels(:)
      end function diag_density_levels_ok

      pure logical module function sponge_source_is_implemented(tag) result(ok)
            !! See submodule [[rdb_config_validate]].
         implicit none
         character(len=*), intent(in) :: tag
      end function sponge_source_is_implemented

      pure logical module function sponge_target_is_implemented(tag) result(ok)
            !! See submodule [[rdb_config_validate]].
         implicit none
         character(len=*), intent(in) :: tag
      end function sponge_target_is_implemented

      pure logical module function sponge_ramp_is_valid(tag) result(ok)
            !! See submodule [[rdb_config_validate]].
         implicit none
         character(len=*), intent(in) :: tag
      end function sponge_ramp_is_valid

      module subroutine build_rdb_schema(cfg, schema)
            !! See submodule [[rdb_config_register_core]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(out) :: schema
      end subroutine build_rdb_schema

      module subroutine register_sim(cfg, schema)
            !! See submodule [[rdb_config_register_core]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_sim

      module subroutine register_grid(cfg, schema)
            !! See submodule [[rdb_config_register_core]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_grid

      module subroutine register_time(cfg, schema)
            !! See submodule [[rdb_config_register_core]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_time

      module subroutine register_mpi(cfg, schema)
            !! See submodule [[rdb_config_register_core]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_mpi

      module subroutine register_logging(cfg, schema)
            !! See submodule [[rdb_config_register_core]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_logging

      module subroutine register_vcoord(cfg, schema)
            !! See submodule [[rdb_config_register_core]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_vcoord

      module subroutine register_kappa_shear(cfg, schema)
            !! See submodule [[rdb_config_register_vmix]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_kappa_shear

      module subroutine register_slopes(cfg, schema)
            !! See submodule [[rdb_config_register_lateral]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_slopes

      module subroutine register_gm(cfg, schema)
            !! See submodule [[rdb_config_register_lateral]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_gm

      module subroutine register_redi(cfg, schema)
            !! See submodule [[rdb_config_register_lateral]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_redi

      module subroutine register_varmix(cfg, schema)
            !! See submodule [[rdb_config_register_lateral]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_varmix

      module subroutine register_meke(cfg, schema)
            !! See submodule [[rdb_config_register_lateral]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_meke

      module subroutine register_tidal_mixing(cfg, schema)
            !! See submodule [[rdb_config_register_vmix]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_tidal_mixing

      module subroutine register_conv(cfg, schema)
            !! See submodule [[rdb_config_register_vmix]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_conv

      module subroutine register_porous(cfg, schema)
            !! See submodule [[rdb_config_register_lateral]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_porous

      module subroutine register_ddiff(cfg, schema)
            !! See submodule [[rdb_config_register_vmix]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_ddiff

      module subroutine register_ocean_tides(cfg, schema)
            !! See submodule [[rdb_config_register_forcing]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_ocean_tides

      module subroutine register_ocean_psurf(cfg, schema)
            !! See submodule [[rdb_config_register_forcing]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_ocean_psurf

      module subroutine register_ocean_cavity_dyn(cfg, schema)
            !! See submodule [[rdb_config_register_forcing]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_ocean_cavity_dyn

      module subroutine register_ocean_cavity_melt(cfg, schema)
            !! See submodule [[rdb_config_register_forcing]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_ocean_cavity_melt

      module subroutine register_epbl(cfg, schema)
            !! See submodule [[rdb_config_register_vmix]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_epbl

      module subroutine register_wavespeed(cfg, schema)
            !! See submodule [[rdb_config_register_vmix]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_wavespeed

      module subroutine register_foxkemper(cfg, schema)
            !! See submodule [[rdb_config_register_lateral]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_foxkemper

      module subroutine register_physics(cfg, schema)
            !! See submodule [[rdb_config_register_core]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_physics

      module subroutine register_output(cfg, schema)
            !! See submodule [[rdb_config_register_core]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_output

      module subroutine register_boundary(cfg, schema)
            !! See submodule [[rdb_config_register_core]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_boundary

      module subroutine register_tracer(cfg, schema)
            !! See submodule [[rdb_config_register_core]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_tracer

      module subroutine register_initial_condition(cfg, schema)
            !! See submodule [[rdb_config_register_core]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_initial_condition

      module subroutine register_nonhydrostatic(cfg, schema)
            !! See submodule [[rdb_config_register_core]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_nonhydrostatic

      module subroutine register_ocean_grid(cfg, schema)
            !! See submodule [[rdb_config_register_core]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_ocean_grid

      module subroutine register_ocean_coriolis(cfg, schema)
            !! See submodule [[rdb_config_register_lateral]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_ocean_coriolis

      module subroutine register_ocean_thermo(cfg, schema)
            !! See submodule [[rdb_config_register_vmix]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_ocean_thermo

      module subroutine register_ocean_forcing(cfg, schema)
            !! See submodule [[rdb_config_register_forcing]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_ocean_forcing

      module subroutine register_ocean_atm_forcing(cfg, schema)
            !! See submodule [[rdb_config_register_forcing]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_ocean_atm_forcing

      module subroutine register_ocean_bulk_flux(cfg, schema)
            !! See submodule [[rdb_config_register_forcing]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_ocean_bulk_flux

      module subroutine register_ocean_ice(cfg, schema)
            !! See submodule [[rdb_config_register_forcing]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_ocean_ice

      module subroutine register_ocean_ice_ic(cfg, schema)
            !! See submodule [[rdb_config_register_forcing]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_ocean_ice_ic

      module subroutine register_ocean_restore(cfg, schema)
            !! See submodule [[rdb_config_register_forcing]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_ocean_restore

      module subroutine register_ocean_geothermal(cfg, schema)
            !! See submodule [[rdb_config_register_forcing]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_ocean_geothermal

      module subroutine register_ocean_sponge(cfg, schema)
            !! See submodule [[rdb_config_register_forcing]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_ocean_sponge

      module subroutine register_ocean_tracers(cfg, schema)
            !! See submodule [[rdb_config_register_forcing]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_ocean_tracers

      module subroutine register_ocean_bt(cfg, schema)
            !! See submodule [[rdb_config_register_lateral]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_ocean_bt

      module subroutine register_ocean_debug(cfg, schema)
            !! See submodule [[rdb_config_register_core]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_ocean_debug

      module subroutine register_ocean_mpi(cfg, schema)
            !! See submodule [[rdb_config_register_core]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_ocean_mpi

      module subroutine register_ocean_wetdry(cfg, schema)
            !! See submodule [[rdb_config_register_vmix]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_ocean_wetdry

      module subroutine register_ocean_pgf(cfg, schema)
            !! See submodule [[rdb_config_register_vmix]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_ocean_pgf

      module subroutine register_ocean_eos(cfg, schema)
            !! See submodule [[rdb_config_register_vmix]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_ocean_eos

      module subroutine register_ocean_bdrag(cfg, schema)
            !! See submodule [[rdb_config_register_vmix]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_ocean_bdrag

      module subroutine register_ocean_tdrag(cfg, schema)
            !! See submodule [[rdb_config_register_vmix]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_ocean_tdrag

      module subroutine register_ocean_hdiff(cfg, schema)
            !! See submodule [[rdb_config_register_lateral]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_ocean_hdiff

      module subroutine register_ocean_hvisc(cfg, schema)
            !! See submodule [[rdb_config_register_lateral]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_ocean_hvisc

      module subroutine register_ocean_vmix(cfg, schema)
            !! See submodule [[rdb_config_register_vmix]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_ocean_vmix

      module subroutine register_ocean_vdiff(cfg, schema)
            !! See submodule [[rdb_config_register_vmix]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_ocean_vdiff

      module subroutine register_ocean_continuity(cfg, schema)
            !! See submodule [[rdb_config_register_vmix]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_ocean_continuity

      module subroutine register_ocean_isopycnal(cfg, schema)
            !! See submodule [[rdb_config_register_lateral]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_ocean_isopycnal

      module subroutine register_ocean_topo(cfg, schema)
            !! See submodule [[rdb_config_register_core]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_ocean_topo

      module subroutine register_ocean_ic(cfg, schema)
            !! See submodule [[rdb_config_register_core]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_ocean_ic

      module subroutine register_ocean_zinit(cfg, schema)
            !! See submodule [[rdb_config_register_core]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_ocean_zinit

      module subroutine register_ocean_data(cfg, schema)
            !! See submodule [[rdb_config_register_core]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_ocean_data

      module subroutine register_ocean_dataovr(cfg, schema)
            !! See submodule [[rdb_config_register_core]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_ocean_dataovr

      module subroutine add_dataovr_entry(g, e, tag, what, units)
            !! See submodule [[rdb_config_register_core]].
         implicit none
         type(nml_group_t), intent(inout) :: g
         type(dataovr_entry_config_t), target, intent(in) :: e
         character(len=*), intent(in) :: tag, what, units
      end subroutine add_dataovr_entry

      module subroutine register_ocean_diag(cfg, schema)
            !! See submodule [[rdb_config_register_core]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_ocean_diag

      module subroutine register_ocean_bc(cfg, schema)
            !! See submodule [[rdb_config_register_forcing]].
         implicit none
         type(config_t), target, intent(in) :: cfg
         type(nml_schema_t), intent(inout) :: schema
      end subroutine register_ocean_bc
   end interface

contains

   pure function bt_dtbt_is_set(bt) result(is_set)
      !! `.true.` when `&ocean_bt_nml dtbt` was given (differs from the
      !! `BT_DTBT_UNSET` sentinel).
      type(ocean_bt_config_t), intent(in) :: bt
      logical :: is_set
      is_set = bt%dtbt > BT_DTBT_UNSET*0.5_wp
   end function bt_dtbt_is_set

   pure function bt_split_derived(bt) result(derived)
      !! `.true.` when `n_inner` is resolved at setup rather than read from
      !! the namelist: `auto_n_inner` or an explicit `dtbt`.
      type(ocean_bt_config_t), intent(in) :: bt
      logical :: derived
      derived = bt%auto_n_inner .or. bt_dtbt_is_set(bt)
   end function bt_split_derived

end module rdb_config
