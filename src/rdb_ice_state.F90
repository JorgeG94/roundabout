!! Sea-ice slot scaffold (SIS2 port, PR 0/1/3a).
module rdb_ice_state
   !! Gated `ocean_sea_ice_t` slot for the SIS2-derived sea-ice model
   !! (`PLAN_SEA_ICE.md`).  PR 0 was plumbing only; PR 1 added the first
   !! prognostic array — the ocean-side frazil supercooling accumulator
   !! `frazil_heat` (filled by `rdb_ice_frazil`) plus its heat-budget
   !! contributor array.  PR 2 added the enthalpy library
   !! (`rdb_ice_enthalpy`).  PR 3a lands the Winton column prognostics:
   !! `part_size`, `m_ice`, `m_snow`, `enth_ice`, `enth_snow`, `sal_ice`
   !! — driven only by the new `rdb_ice_column` test suite (nothing in
   !! the driver calls the column kernels yet).  A run with
   !! `&ocean_ice_nml enable = .false.` (the default) never calls into
   !! this module beyond the parent gates, so existing configurations
   !! stay byte-identical.
   !!
   !! PR 3b lands the ocean<->ice coupling seam: `rdb_ice_frazil_uptake`
   !! spends the frazil bank as new category-1 ice (`m_frozen_diag`) and
   !! `rdb_ice_ocean_coupler` refreshes the surface `Q_salt` field from
   !! the resulting brine-rejection rate (`salt_flux_diag`) — both new
   !! fields below.
   !!
   !! PR 3c lands the coupleable atmospheric-forcing seam
   !! (`atm_sf0`/`atm_dsfdt`/`atm_sw_dn`, filled by `rdb_ice_atm_forcing`),
   !! the ocean->ice basal-flux field (`fb`, filled by
   !! `rdb_ice_basal_flux`, which also samples the `sst_seam`/`ssurf_seam`/
   !! `tfw_seam` scratch reused by `rdb_ice_thermo_driver`), the column's
   !! per-category scratch outputs (`tsurf_out`/`h2o_ocn_to_ice`/
   !! `h2o_ice_to_ocn`/`heat_to_ocn`/`sw_thru`), and the melt-side
   !! ocean-coupling diags `heat_flux_diag`/`m_melt_diag` (parallel to
   !! `salt_flux_diag`/`m_frozen_diag`) that `rdb_ice_ocean_coupler`'s
   !! `ice_ocean_heat_flux` and the driver's resume fold consume.
   !!
   !! PR 4a lands the multi-category ITD: `h_lim`/`mh_lim` (SIS2 category
   !! thickness/mass bounds, `ice_itd_category_bounds`, THIS module) and the
   !! `fb_part_sum` scratch (the `rdb_ice_thermo_driver` fb-cover snapshot).
   !! Mirrors `ocean_meke_t` (the default-off, gated- allocation,
   !! restart-persistent precedent). PR-58 adds an optional `hlim_cfg`
   !! override of the bin edges (`ice_itd_category_bounds`'s `hlim_vals`
   !! dummy) — see `hlim_cfg` below.
   !!
   !! PR 4b lands horizontal category transport (`rdb_ice_transport`):
   !! `u_ice`/`v_ice` (C-grid face velocities, persistent, allocated
   !! whenever the ice slot is live — PR 5's EVP dynamics will fill them
   !! instead of the v1 ocean-surface-velocity sampler) and the
   !! transport-only scratch workspace (`mca_ice`/`mca_snow`,
   !! `uh_ice`/`vh_ice`/`uh_snow`/`vh_snow`, `htot_work`/`hl_x_work`/
   !! `hr_x_work`/`hl_y_work`/`hr_y_work`/`uhtot_work`/`vhtot_work`)
   !! allocated ONLY when
   !! `&ocean_ice_nml transport = .true.` (the `this%transport` flag,
   !! set from config before `init`, mirrors `enable`/`ncat`). Later
   !! rungs: PR 5 C-grid EVP dynamics (replaces the u_ice/v_ice filler).
   !!
   !! PR 5 lands C-grid EVP rheology (`rdb_ice_evp`): `dynamics` (master
   !! flag, set from config before `init` like `transport`) gates the
   !! driver's `ice_evp_step` call (which now WRITES `u_ice`/`v_ice`
   !! instead of the transport sampler) and turns on the stress fields
   !! `str_d`/`str_t` (T-cells) + `str_s` (corners) — the H&D (1997)
   !! divergence/tension/shear stress tensor components, persistent +
   !! restart-carried. `tau_a_x`/`tau_a_y` are a configure-time snapshot
   !! of the wind stress felt by the ice (no ice-specific bulk drag law,
   !! D7); `fxoc`/`fyoc` are the subcycle-averaged ice->ocean stress that
   !! `ice_ocean_stress_flux` (`rdb_ice_ocean_coupler`) blends into
   !! `ocean_surface_stress_t%tau_x/tau_y` by concentration. All seven
   !! fields are allocated UNCONDITIONALLY whenever the ice slot is live
   !! (same contract as `u_ice`/`v_ice` — cheap 2-D arrays), zero-init.
   !!
   !! PR 63 closes the resume hole PR 5 documented (F4): `fxoc`/`fyoc`
   !! restart-carry the subcycle-averaged drag, but `ice_ocean_stress_flux`
   !! blends it against the PRE-thermo `ci` of that same step, and the
   !! checkpoint only ever sees the POST-thermo `ci` — so reconstructing
   !! the blend at resume from the checkpointed masses gives a different
   !! answer whenever thermo/transport changed `ci` after the blend.
   !! `tau_ocn_x`/`tau_ocn_y` mirror the exact blended value
   !! `ice_ocean_stress_flux` last wrote into
   !! `ocean_surface_stress_t%tau_x/tau_y` — restart-carried, so the
   !! resume path COPIES it instead of recomputing it (formula-agnostic:
   !! the carry stays exact under any future change to the blend
   !! formula). `tau_ocn_valid` is the host-only "was a blend ever
   !! written" flag (0.0 fresh run / pre-PR-63 checkpoint, 1.0 otherwise)
   !! — SIS2's `query_initialized('stress_mag')` fallback idiom
   !! (`ice_type.F90:264`, `ice_model.F90:2410`), needed because an
   !! unconditional copy on a fresh run would hand the ocean a zeroed
   !! wind stress on step 1. See `ice_ocean_stress_resume_apply`
   !! (`rdb_ice_ocean_coupler`).
   !!
   !! **Two-mode convention (`part_size`/`m_ice`/`m_snow`), load-bearing.**
   !! `ncat == 1` is the LEGACY LUMPED mode landed by PR 3b/3c:
   !! `part_size` is never maintained (`part_size(:,:,0)=1` forever,
   !! category 1 unused), and `m_ice`/`m_snow`/`enth_ice`/`enth_snow`/
   !! `sal_ice` are per unit CELL area. All existing code paths
   !! (`ice_frazil_uptake_impl`, the ncat=1 `ice_thermo_driver_reduce_impl`)
   !! run byte-for-byte unchanged under this mode — bit-identity with PR 3c
   !! is guaranteed BY CONSTRUCTION (dispatch happens only at the three
   !! shim points in `rdb_ice_frazil_uptake`/`rdb_ice_thermo_driver`/
   !! `rdb_ice_itd`, never inside a kernel).  `ncat > 1` is the SIS2 ITD
   !! mode (PR 4a): `part_size` is LIVE (Σ_cat part_size = 1, category 0 =
   !! open water), and `m_ice`/`m_snow` become per unit ICE-COVERED-area
   !! intensive quantities (SIS2 `mH_ice`/`mH_snow`) — the frazil uptake
   !! annexes open water into an occupied/thinnest category
   !! (`ice_frazil_uptake_multicat_impl`), the thermo-driver reduce
   !! part-weights the per-cell diags, and `ice_adjust_categories`
   !! (`rdb_ice_itd`) restores the ITD (whole-category "move all of it"
   !! shift, `SIS_transport.F90:611-891`) once per thermo window.
   !! `enth_ice`/`enth_snow`/`sal_ice` stay per-MASS intensive in BOTH
   !! modes (unchanged meaning). Restart files written by PR-3c at
   !! ncat>1 carry the OLD lumped meaning — the feature branch is
   !! unreleased, so this break is accepted (not carried forward).
   !!
   !! PR 26 lands the snowfall source term: `has_snowfall` (host-side
   !! dispatch gate, latched from `&ocean_ice_nml snowfall /= 0` before
   !! `init`, mirrors `transport`/`dynamics`), the fourth atmospheric-seam
   !! field `atm_fprec` (kg/m^2/s, filled by `rdb_ice_atm_forcing`
   !! alongside `atm_sw_dn`, consumed by `rdb_ice_column`'s
   !! `ice_snow_accumulate`), and the ocean-delivery pair
   !! `snow_part_ocn`/`fprec_ocn_diag` (the PRE-column ice-free-cover
   !! snapshot and the resulting ocean-bound frozen-precipitation share,
   !! filled by `rdb_ice_thermo_driver`/`rdb_ice_snow` — same
   !! PRE-column-snapshot contract as `fb_part_sum`). All three arrays
   !! are allocated UNCONDITIONALLY whenever the ice slot is live (cheap
   !! 2-D arrays, same contract as the rest of the `atm_*` seam);
   !! `has_snowfall` is what makes `snowfall=0` byte-identical by
   !! construction. `fprec_ocn_diag` is the binding seam PR-16 consumes
   !! as a `net_massin` source (see `rdb_ice_snow`'s module docstring).
   !!
   !! PR 27 lands the Archimedes freeboard snow-ice flood: `snow_ice`
   !! (host-side gate, latched from `&ocean_ice_nml snow_ice` BEFORE
   !! `init`, same convention as `has_snowfall`/`transport`/`dynamics`)
   !! and `snow_to_ice` (per-category kg/m² output, SIS2 `SN2IC`,
   !! allocated UNCONDITIONALLY whenever the ice slot is live, same
   !! contract as `sw_thru` — filled by `ice_thermo_columns`, zeroed
   !! unconditionally every thermo step, `copyin`/`delete` mapped, NOT
   !! restart-carried, no downstream consumer yet). `snow_ice=.false.`
   !! (default) is what makes the flood byte-identical by construction.
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_mem_report, only: arr_bytes
   use rdb_ice_column, only: ICE_BULK_SALINITY, ICE_NK_MAX, ICE_RHO_ICE
   use rdb_ice_enthalpy, only: ICE_CP_BRINE, ICE_CP_ICE
   use, intrinsic :: iso_fortran_env, only: int64
   implicit none
   private

   public :: ocean_sea_ice_t
   public :: evp_workspace_t
   public :: ice_itd_category_bounds
   public :: ice_cell_concentration_impl

   type :: evp_workspace_t
      !! Persistent EVP-rheology scratch (SIS2 C-grid dynamics), hung off
      !! `ocean_sea_ice_t` as the `evp_ws` slot — eagerly allocated in the
      !! ice state's `init` and GPU-mapped via `enter_data`, so it matches
      !! the rest of `src/core/ocean/` (which carries ZERO module-level
      !! `save` allocatables — every kernel puts its scratch on a state
      !! DT).  Replaces the retired `rdb_ice_evp` module-level `save`
      !! workspace + its lazy `evp_workspace_ensure`/`ice_evp_cleanup`
      !! lifecycle.  Allocated / mapped ONLY when `ocean_sea_ice_t%dynamics`
      !! is on (the workspace is EVP-dynamics-only).
      !!
      !! `rdb_ice_evp`'s `ice_evp_dynamics` shim unpacks these components
      !! into the explicit-shape flat-impl args of `ice_evp_dynamics_impl`
      !! — NEVER read `ws%component` inside a `do concurrent` (that triggers
      !! per-launch descriptor copies).
      !!
      !! F3 (dummy-argument aliasing safety): `mis_in_w`/`mice_in_w`/
      !! `ci_in_w` are the `ice_evp_step` gather-INPUT buffers, DELIBERATELY
      !! separate from `mis_w`/`mice_w`/`ci_w` (which `ice_evp_dynamics`
      !! fills from its `intent(in)` gathered dummies) — the concentration
      !! gather never writes the same arrays the core reads-then-fills.
      integer :: nx = 0
         !! Cached T-cell x extent (nx_total).
      integer :: ny = 0
         !! Cached T-cell y extent (ny_total).
      real(wp), allocatable :: mis_w(:, :), mice_w(:, :), ci_w(:, :)
      real(wp), allocatable :: mis_in_w(:, :), mice_in_w(:, :), ci_in_w(:, :)
      real(wp), allocatable :: pres_mice_w(:, :), del_sh_min_pr_w(:, :)
      real(wp), allocatable :: sh_dd_w(:, :), sh_dt_w(:, :)
      real(wp), allocatable :: zeta_w(:, :), del_sh_w(:, :)
      real(wp), allocatable :: mask_t_w(:, :)
      real(wp), allocatable :: mi_u_w(:, :), mask_u_w(:, :), u_tmp_w(:, :)
      real(wp), allocatable :: mi_v_w(:, :), mask_v_w(:, :)
      real(wp), allocatable :: sh_ds_w(:, :), mi_ratio_a_q_w(:, :)
      real(wp), allocatable :: q_w(:, :), mask_q_w(:, :)
      real(wp), allocatable :: a_u_w(:, :), a_v_w(:, :)
         !! PR 62: face ice concentration, `a_u_w(i,j) = 0.5*(ci_w(i-1,j) +
         !! ci_w(i,j))` (same expression/edge convention as
         !! `evp_mi_face_impl`'s `mi_u`, and as `ice_ocean_stress_flux_impl`'s
         !! `a_u`). Filled by `evp_mi_face_impl(ci_w, ...)` ONLY when
         !! `par%a_face_stress`; read ONLY under the same gate — uninitialised
         !! device memory is unreachable when the knob is off.
   contains
      procedure, non_overridable :: init => evp_workspace_init
      procedure, non_overridable :: enter_data => evp_workspace_enter_data
      procedure, non_overridable :: exit_data => evp_workspace_exit_data
      procedure, non_overridable :: destroy => evp_workspace_destroy
      procedure, non_overridable :: bytes => evp_workspace_bytes
   end type evp_workspace_t

   type :: ocean_sea_ice_t
      !! Sea-ice state slot.  All fields default to the inert
      !! (`enable=.false.`) configuration so an ocean run that never sets
      !! `&ocean_ice_nml` is byte-identical to a build without this slot.
      logical :: is_init = .false.
         !! True between `init` and `destroy`; gate on this (never on
         !! `allocated`, which misses the GPU mapping).
      logical :: enable = .false.
         !! Master switch (`&ocean_ice_nml enable`).  Off (default) ⇒ the
         !! parent state never inits / maps / steps this slot ⇒ no-op.
      integer :: ncat = 5
         !! Number of ice thickness categories (`&ocean_ice_nml ncat`).
         !! Category 0 is open water; SIS2 default 5.  Sizes the
         !! per-category arrays when PR 3+ adds them.
      integer :: nk_ice = 2
         !! Vertical ice layers per category (`&ocean_ice_nml nk_ice`).
         !! 2 = Winton two-layer (v1).  Bottom-up once arrays land
         !! (k=1 ice bottom — SIS2's top-down index is flipped at this
         !! state layer, per the plan).
      logical :: transport = .false.
         !! Horizontal category transport switch (`&ocean_ice_nml
         !! transport`), set from config BEFORE `init` (same convention as
         !! `enable`/`ncat`).  Gates allocation of the transport-only
         !! workspace below (memory: the per-category face fluxes are
         !! ~4·ncat·N reals — do not pay it by default).  `u_ice`/`v_ice`
         !! are unconditional (allocated whenever the slot is live) since
         !! PR 5's EVP dynamics will need them regardless of whether
         !! category transport is on.
      logical :: dynamics = .false.
         !! EVP rheology master switch (`&ocean_ice_nml dynamics`), set
         !! from config BEFORE `init` (same convention as `enable`/
         !! `transport`). Gates the driver's `ice_evp_step` call and the
         !! transport sampler skip; does NOT gate any allocation here —
         !! `str_d`/`str_t`/`str_s`/`tau_a_*`/`fxoc`/`fyoc`/`tau_ocn_*`
         !! are unconditional (cheap 2-D arrays), same contract as
         !! `u_ice`/`v_ice`.
      logical :: has_snowfall = .false.
         !! Host-side dispatch gate for the PR 26 snowfall source term,
         !! latched from `(cfg%ocean%ice%snowfall /= 0.0_wp)` BEFORE `init`
         !! (same convention as `transport`/`dynamics`).  Gates the
         !! pre-column `snow_part_ocn` fill and the `ice_snowfall_ocean_
         !! share` contributor call; does NOT gate allocation of `atm_fprec`/
         !! `snow_part_ocn`/`fprec_ocn_diag` (cheap 2-D arrays, same
         !! contract as `atm_sw_dn`/`fb_part_sum`) — the `has_sw`/`has_heat`
         !! idiom that makes `snowfall=0` byte-identical by construction,
         !! not by IEEE argument.
      real(wp) :: tau_ocn_valid = 0.0_wp
         !! HOST-ONLY (PR 63). Exactly 0.0 or 1.0. `1.0` iff
         !! `tau_ocn_x`/`tau_ocn_y` hold a blend written by this run's
         !! `ice_ocean_stress_flux` or restored from a checkpoint that
         !! carried one; `0.0` on a fresh run or a resume from a
         !! pre-PR-63 checkpoint. Set by `ice_ocean_stress_flux`, carried
         !! by `register_scalar` (which forces `device_mapped=.false.`),
         !! read only by `ice_ocean_stress_resume_apply`. NEVER
         !! device-mapped and NEVER read in a kernel — the SIS2
         !! `query_initialized('stress_mag')` idiom
         !! (`ice_model.F90:2410`), one host scalar instead of a
         !! per-field registry query.
      logical :: snow_ice = .false.
         !! Archimedes freeboard snow-ice flood master switch
         !! (`&ocean_ice_nml snow_ice`, PR 27), latched from config
         !! BEFORE `init` (same convention as `transport`/`dynamics`).
         !! Read by `ice_thermo_driver_step` and passed to
         !! `ice_thermo_columns` as `do_snow_ice`; does NOT gate
         !! allocation of `snow_to_ice` (cheap per-category array, same
         !! contract as `sw_thru`) — `.false.` (default) is what makes
         !! `snow_ice=.false.` byte-identical by construction.

      real(wp), allocatable :: u_ice(:, :)
         !! Ice velocity at C-grid u-faces (m/s), shape (nx_total+1,
         !! ny_total).  v1 INTERIM: sampled from the ocean surface layer
         !! every transport call (`ice_transport_step` Phase 0); PR 5's
         !! EVP dynamics will write this field instead of the sampler, same
         !! faces.  Persistent (unconditional allocation whenever the ice
         !! slot is live).
      real(wp), allocatable :: v_ice(:, :)
         !! Ice velocity at C-grid v-faces (m/s), shape (nx_total,
         !! ny_total+1).  Same interim/persistent contract as `u_ice`.

      ! ---- PR 5: C-grid EVP ice dynamics ----
      real(wp), allocatable :: str_d(:, :)
         !! Divergence stress tensor component [Pa*m], T-cells, shape
         !! (nx_total, ny_total). SIS2 `CS%str_d`. Persistent, restart-
         !! carried (optional). Relaxes to `-0.5*pres_mice*mice` at zero
         !! strain rate (the elliptical yield curve's pressure intercept —
         !! ice-free cells have `mice=0` so this is exactly 0 there, an
         !! emergent property, never a Dirichlet ice-edge special case).
      real(wp), allocatable :: str_t(:, :)
         !! Tension stress tensor component [Pa*m], T-cells, shape
         !! (nx_total, ny_total). SIS2 `CS%str_t`. Persistent, restart-
         !! carried (optional).
      real(wp), allocatable :: str_s(:, :)
         !! Shearing stress tensor component (cross term) [Pa*m],
         !! corners, shape (nx_total+1, ny_total+1). SIS2 `CS%str_s`.
         !! Persistent, restart-carried (optional).
      real(wp), allocatable :: tau_a_x(:, :)
         !! Atmospheric (wind) stress felt by the ice, snapshotted onto
         !! u-faces [Pa], shape (nx_total+1, ny_total). Configure-time
         !! copy of the ocean's wind-stress field (D7: one atmospheric
         !! stress field, no ice-specific bulk drag law) — NOT restart-
         !! carried (rebuilt at configure from the wind-stress config on
         !! every run, resumed or fresh).
      real(wp), allocatable :: tau_a_y(:, :)
         !! Ditto, v-faces, shape (nx_total, ny_total+1).
      real(wp), allocatable :: fxoc(:, :)
         !! Subcycle-averaged ice->ocean stress [Pa], u-faces, shape
         !! (nx_total+1, ny_total). SIS2 `fxoc`. Persistent, restart-
         !! carried (optional): `ice_ocean_stress_flux` blends this into
         !! `ocean_surface_stress_t%tau_x` every outer step. PR 63:
         !! this registration is no longer what makes the resume
         !! bit-exact — `tau_ocn_x`/`tau_ocn_y` below carry the blend's
         !! OUTPUT directly, so the resume path never needs to
         !! reconstruct it from `fxoc` + the checkpointed `ci`. Kept
         !! restart-carried for diagnostic continuity across a resume
         !! (like `m_frozen_diag`/`m_melt_diag`) — its only cross-step
         !! reader inside a call is `ice_evp_dynamics`'s own zeroing at
         !! entry, so it carries no state INTO the EVP solve either.
      real(wp), allocatable :: fyoc(:, :)
         !! Ditto, v-faces, shape (nx_total, ny_total+1). SIS2 `fyoc`.
      real(wp), allocatable :: tau_ocn_x(:, :)
         !! PR 63. Mirror of the exact value `ice_ocean_stress_flux` last
         !! wrote into `ocean_surface_stress_t%tau_x` [Pa], u-faces, shape
         !! (nx_total+1, ny_total). PERSISTENT, restart-carried
         !! (optional). Formula-agnostic by construction — it is a COPY
         !! of the blend's output, never a recomputation, so it stays
         !! exact under any future change to the blend formula (SIS2
         !! `Ice%flux_u`, `ice_type.F90:242`). Valid to read on the host
         !! only after `!$acc update self` of the COMPONENT array (never
         !! the aggregate `ice` — commit `72152870`). See
         !! `ice_ocean_stress_resume_apply` (`rdb_ice_ocean_coupler`).
      real(wp), allocatable :: tau_ocn_y(:, :)
         !! Ditto, v-faces, shape (nx_total, ny_total+1). SIS2
         !! `Ice%flux_v`.

      type(evp_workspace_t) :: evp_ws
         !! Persistent EVP-dynamics scratch (`rdb_ice_evp`'s subcycle
         !! kernels).  Allocated + GPU-mapped ONLY when `dynamics` is on
         !! (see `evp_workspace_t`) — the successor to the retired
         !! `rdb_ice_evp` module-level `save` workspace.  Rides
         !! `ocean_state_enter_data` via this slot's `enter_data`.

      ! ---- PR 4b: category transport workspace (allocated iff `transport`) ----
      real(wp), allocatable :: mca_ice(:, :, :)
         !! Cell-averaged ice mass per category (CAS space, kg/m² of CELL
         !! area — SIS2 `CAS%m_ice`), shape (nx_total, ny_total, ncat).
         !! SCRATCH: filled at Phase 1 (IST->CAS), consumed/updated through
         !! Phase 2, discarded after Phase 3 (CAS->IST).
      real(wp), allocatable :: mca_snow(:, :, :)
         !! Cell-averaged snow mass per category (CAS space, kg/m² of CELL
         !! area), shape (nx_total, ny_total, ncat).  Same lifecycle as
         !! `mca_ice`.
      real(wp), allocatable :: uh_ice(:, :, :)
         !! Per-category zonal ice-mass face transport (kg/s), shape
         !! (nx_total+1, ny_total, ncat).  SCRATCH, recomputed every pass.
      real(wp), allocatable :: vh_ice(:, :, :)
         !! Per-category meridional ice-mass face transport (kg/s), shape
         !! (nx_total, ny_total+1, ncat).  SCRATCH.
      real(wp), allocatable :: uh_snow(:, :, :)
         !! Per-category zonal snow-mass face transport (kg/s), shape
         !! (nx_total+1, ny_total, ncat).  SCRATCH.
      real(wp), allocatable :: vh_snow(:, :, :)
         !! Per-category meridional snow-mass face transport (kg/s), shape
         !! (nx_total, ny_total+1, ncat).  SCRATCH.
      real(wp), allocatable :: htot_work(:, :)
         !! Category-SUMMED mass, 2-D scratch (kg/m² of CELL area), shape
         !! (nx_total, ny_total).  Reused for the ice solve then the snow
         !! solve, per direction (SIS2 `zonal_mass_flux`'s `htot`).
      real(wp), allocatable :: hl_x_work(:, :)
         !! PPM face-edge reconstruction of the category-SUMMED zonal mass,
         !! FACE-shaped `(nx_total+1, ny_total)` (mirrors `rdb_continuity`'s
         !! `h_face_left_x`).  `hl_x_work(i,j)` is the value AT east face
         !! `i` extrapolated from the LEFT cell `(i-1)` (i.e. cell i-1's
         !! own right/downwind edge); `hr_x_work(i,j)` is AT face `i` from
         !! the RIGHT cell `i` (cell i's own left edge).  FACE shape is
         !! load-bearing: the array-edge fallback writes index `nx+1`
         !! (out of bounds for a `(nx,ny)` array — the pre-fix bug).
         !! SCRATCH, reused per medium.
      real(wp), allocatable :: hr_x_work(:, :)
         !! Right-cell zonal face-edge companion to `hl_x_work`,
         !! `(nx_total+1, ny_total)`.  SCRATCH.
      real(wp), allocatable :: hl_y_work(:, :)
         !! North-face PPM edge, `(nx_total, ny_total+1)` — the meridional
         !! twin of `hl_x_work` (`h_face_left_y` convention). SCRATCH.
      real(wp), allocatable :: hr_y_work(:, :)
         !! North-face right-cell PPM edge, `(nx_total, ny_total+1)`.
         !! SCRATCH.
      real(wp), allocatable :: uhtot_work(:, :)
         !! Category-SUMMED zonal face transport (kg/s), shape
         !! (nx_total+1, ny_total).  SCRATCH, reused per medium.
      real(wp), allocatable :: vhtot_work(:, :)
         !! Category-SUMMED meridional face transport (kg/s), shape
         !! (nx_total, ny_total+1).  SCRATCH, reused per medium.
      real(wp), allocatable :: tr_flux_x_work(:, :, :)
         !! Per-category zonal DONOR-VALUE buffer (`val(donor of face
         !! I)`), shape (nx_total+1, ny_total, ncat).  GPU-race-free
         !! gather/scatter split (mirrors `rdb_continuity`'s
         !! `Tr_face_left_x` pattern): a first `do concurrent` over FACES
         !! reads `val` only at each face's OWN donor cell and writes here
         !! (never the same cell two faces disagree on); the cell-update
         !! kernel then reads only this face buffer (to recover `val_e`/
         !! `val_w` directly — NOT flux-weighted, so no back-division is
         !! needed) plus its OWN cell's prior mass/value — no kernel ever
         !! reads a NEIGHBOUR cell's `val` while another iteration writes
         !! that neighbour.  Reused across every riding field (`m_ice`,
         !! each `enth_ice`/`sal_ice` layer, `enth_snow`) one at a time
         !! within a pass — SCRATCH.
      real(wp), allocatable :: tr_flux_y_work(:, :, :)
         !! Meridional twin of `tr_flux_x_work`, shape (nx_total,
         !! ny_total+1, ncat).

      real(wp), allocatable :: frazil_heat(:, :)
         !! Ocean-side frazil supercooling BANK (J/m², T-cells, shape
         !! (nx_total, ny_total)).  `ice_frazil_accumulate` clamps the
         !! surface layer at the freezing point and deposits the removed
         !! heat deficit `ρ·Cp·h·(T_f − T)⁺` here; the ice thermodynamics
         !! (PR 3) will SPEND it as new-ice formation.  Persistent state
         !! (accumulates across steps, restart-carried) — distinct from
         !! the drain-and-zero budget accumulator below.
      real(wp), allocatable :: heat_budget_frazil(:, :, :)
         !! hTr (K·m) change per cell per outer step attributed to the
         !! frazil clamp — the `budgets` HEAT contributor mirroring the
         !! `heat_budget_surface` convention (positive = heat added to
         !! the ocean; the clamp WARMS the surface layer up to T_f).
         !! Shape (nx_total, ny_total, 1): the clamp only ever touches
         !! the surface layer, and the drain integrates over k anyway.
         !! Drained-and-zeroed by `ocean_budgets_drain_contributors`.

      ! ---- PR 3a: Winton column prognostics ----
      real(wp), allocatable :: part_size(:, :, :)
         !! Fractional ice-category area, shape (nx_total, ny_total,
         !! 0:ncat). Category 0 = OPEN WATER; Σ_cat part_size = 1.
         !! TWO-MODE CONVENTION (module docstring): `ncat==1` — legacy
         !! lumped mode, deliberately UNTOUCHED by every kernel
         !! (`part_size(:,:,0)=1` forever, category 1 unused, thermo runs
         !! per unit CELL area). `ncat>1` — SIS2 ITD mode (PR 4a), LIVE:
         !! `ice_frazil_uptake_multicat_impl` annexes open water into an
         !! occupied category, and `ice_adjust_categories` restores the
         !! per-cell partition each thermo window.
      real(wp), allocatable :: m_ice(:, :, :)
         !! Total ice mass per category (kg/m²), shape (nx_total,
         !! ny_total, ncat). TWO-MODE CONVENTION: `ncat==1` — per unit
         !! CELL area (legacy). `ncat>1` — per unit ICE-COVERED area
         !! (SIS2 `mH_ice`, intensive quantity DIVIDED by `part_size`,
         !! not multiplied). Layer masses are always `m_ice/nk_ice` in
         !! EITHER mode — the column step always ends with
         !! `ice_rebalance_layers`, so equal-mass layering is a state
         !! invariant regardless of the area convention.
      real(wp), allocatable :: m_snow(:, :, :)
         !! Snow mass per category (kg/m²), shape (nx_total, ny_total,
         !! ncat). Same two-mode area convention as `m_ice` (per-CELL at
         !! ncat==1, per-ICE-area at ncat>1).
      real(wp), allocatable :: enth_ice(:, :, :, :)
         !! Ice specific enthalpy (J/kg), shape (nx_total, ny_total,
         !! ncat, nk_ice). BOTTOM-UP: k=1 = ice BOTTOM (ocean side),
         !! k=nk_ice = ice TOP (atm/snow side) — the opposite of the
         !! SIS2/`rdb_ice_column` internal top-down convention; the
         !! flip happens at the `ice_column_step` gather/scatter
         !! boundary (`rdb_ice_column`'s module docstring, TRAP #2).
      real(wp), allocatable :: enth_snow(:, :, :, :)
         !! Snow specific enthalpy (J/kg), shape (nx_total, ny_total,
         !! ncat, 1).
      real(wp), allocatable :: sal_ice(:, :, :, :)
         !! Ice bulk salinity (PSU), shape (nx_total, ny_total, ncat,
         !! nk_ice), same BOTTOM-UP k convention as `enth_ice`.

      ! ---- PR 3b: frazil->ice uptake + brine-rejection coupling ----
      real(wp), allocatable :: m_frozen_diag(:, :)
         !! New frazil-ice mass formed in the LAST thermo window (kg per
         !! m² of CELL area — v1 category convention, see
         !! `rdb_ice_frazil_uptake`), shape (nx_total, ny_total).
         !! OVERWRITTEN (not accumulated) each uptake; zero where no
         !! freezing occurred. Restart-carried for post-resume
         !! diagnostic continuity.
      real(wp), allocatable :: salt_flux_diag(:, :)
         !! Ice -> ocean virtual salt flux from the last uptake
         !! (PSU·kg/m²/s, POSITIVE SALINIFIES — the surface-flux Q_salt
         !! convention), shape (nx_total, ny_total). PERSISTENT +
         !! restart-carried: `rdb_ice_ocean_coupler` refills `Q_salt`
         !! from this field each thermo step, and the driver's
         !! configure-time resume fold re-applies it so the first
         !! post-resume window sees the same flux the uninterrupted run
         !! would have. PR 3c: the column driver ADDS its net-melt
         !! contribution here on top of what the frazil uptake (which
         !! runs first each window) wrote — see `rdb_ice_thermo_driver`.

      ! ---- PR 3c: coupleable atmospheric-forcing seam ----
      real(wp), allocatable :: atm_sf0(:, :)
         !! Net upward surface flux at T_surf=0 (W/m^2), the SEB intercept
         !! `sf_0` fed to ice_thermo_columns. Filled every thermo step by
         !! ice_atm_forcing_restoring (v1) or the future bulk-flux subsystem.
         !! SCRATCH (recomputed each step) — NOT restart-carried.
      real(wp), allocatable :: atm_dsfdt(:, :)
         !! SEB slope dSF/dT (W/m^2/K) — `dsf_dt`. Same seam contract.
      real(wp), allocatable :: atm_sw_dn(:, :)
         !! Downwelling shortwave into the ice top (W/m^2) — `sw_dn`. Same.
      real(wp), allocatable :: atm_fprec(:, :)
         !! PR 26: frozen-precipitation rate onto the ice top (kg/m^2/s),
         !! >= 0 — the seam `ice_atm_forcing_restoring` fills from
         !! `&ocean_ice_nml snowfall` exactly like `atm_sw_dn` <- `sw_down`.
         !! Same SCRATCH contract (refilled every thermo step, NOT restart-
         !! carried). PR 55 fills this from data instead of a uniform
         !! scalar; the shape/units/lifecycle are fixed by this PR (see
         !! `rdb_ice_snow`'s module docstring for the binding seam spec).

      ! ---- PR 3c: ocean -> ice basal flux + the sample seam it shares
      !       with the column driver ----
      real(wp), allocatable :: fb(:, :)
         !! Ocean -> ice-base heat flux (W/m^2), filled by
         !! `rdb_ice_basal_flux%ice_compute_basal_flux` from the
         !! above-freezing SST. SCRATCH (recomputed each step) — NOT
         !! restart-carried. Same lifecycle as the `atm_*` seam.
      real(wp), allocatable :: sst_seam(:, :)
         !! Sampled sea-surface temperature (degC) at the one-step-lagged
         !! surface, filled by `ice_compute_basal_flux` and reused by
         !! `rdb_ice_thermo_driver` — one sample, shared by `fb` and the
         !! column's ocean-side inputs. SCRATCH.
      real(wp), allocatable :: ssurf_seam(:, :)
         !! Sampled sea-surface salinity (PSU), same sharing contract as
         !! `sst_seam`. SCRATCH.
      real(wp), allocatable :: tfw_seam(:, :)
         !! Sampled seawater freezing temperature at the surface (degC),
         !! same sharing contract as `sst_seam`. SCRATCH.

      ! ---- PR 3c: column per-category scratch outputs ----
      real(wp), allocatable :: tsurf_out(:, :, :)
         !! Column skin temperature (degC), shape (nx_total, ny_total,
         !! ncat) — diagnostic, recomputed every thermo step. SCRATCH.
      real(wp), allocatable :: h2o_ocn_to_ice(:, :, :)
         !! Mass FROZEN from the ocean onto the ice base (kg/m²), shape
         !! (nx_total, ny_total, ncat), >= 0. SCRATCH.
      real(wp), allocatable :: h2o_ice_to_ocn(:, :, :)
         !! Meltwater mass to the ocean (kg/m²), shape (nx_total,
         !! ny_total, ncat), >= 0. SCRATCH.
      real(wp), allocatable :: heat_to_ocn(:, :, :)
         !! Leftover melt energy drained back to the ocean (J/m²), shape
         !! (nx_total, ny_total, ncat), >= 0. SCRATCH.
      real(wp), allocatable :: sw_thru(:, :, :)
         !! Shortwave transmitted through the ice to the ocean (W/m²),
         !! shape (nx_total, ny_total, ncat), >= 0. Per-category column
         !! output; PR 31 reduces it to the per-cell `sw_thru_diag` below
         !! (which IS coupled to the ocean surface-flux slot). SCRATCH
         !! (recomputed every thermo step, `delete`-mapped).
      real(wp), allocatable :: snow_to_ice(:, :, :)
         !! PR 27: Archimedes freeboard snow-ice conversion this thermo
         !! window (kg/m² per category), shape (nx_total, ny_total,
         !! ncat), >= 0 — SIS2 `SN2IC`. Filled by `ice_thermo_columns`
         !! (zeroed unconditionally, written only when `snow_ice` is on
         !! and the column floods). Same lifecycle as `sw_thru`: filled,
         !! not yet consumed by any diagnostic/coupling seam
         !! (a future diagnostics PR wires it in — `sw_thru` precedent).
         !! NOT restart-carried (window scratch).

      ! ---- PR 3c: melt-side ocean coupling diags ----
      real(wp), allocatable :: heat_flux_diag(:, :)
         !! Ice -> ocean net surface heat flux from the last thermo window
         !! (W/m^2, POSITIVE DOWN into the ocean — the Q_heat convention),
         !! shape (nx_total, ny_total). Q_heat = heat_to_ocn/dt_therm - fb
         !! (see `rdb_ice_thermo_driver`). PERSISTENT + restart-carried
         !! (mirrors salt_flux_diag): `rdb_ice_ocean_coupler` refills
         !! Q_heat from this each thermo step, and the driver resume fold
         !! re-applies it.  This is the NON-shortwave heat share only — the
         !! penetrating shortwave `sw_thru_diag` below is a SEPARATE field
         !! (PR 31), never folded in here (that would double-count against
         !! the `q_sw` component / the components-off Q_heat fold — see
         !! `ice_ocean_sw_flux` in `rdb_ice_ocean_coupler`).
      real(wp), allocatable :: sw_thru_diag(:, :)
         !! PR 31: ice -> ocean shortwave transmitted through the ice to
         !! the water below, reduced to per unit CELL area from the
         !! per-category `sw_thru` (W/m^2, >= 0, POSITIVE DOWN), shape
         !! (nx_total, ny_total).  Refilled every thermo step by
         !! `ice_thermo_driver_reduce{,_multicat}_impl` (area-weighted at
         !! ncat>1, lumped at ncat==1 — the same convention as
         !! `heat_flux_diag`), consumed by `ice_ocean_sw_flux`
         !! (`rdb_ice_ocean_coupler`): with `&ocean_forcing_nml
         !! enable_components` on it fills the `q_sw` surface-flux
         !! component (assembler sums it into `Q_heat`); off, it is added
         !! directly into `Q_heat`.  PERSISTENT + restart-carried (mirrors
         !! `heat_flux_diag`): the coupler rebuilds the ocean SW every
         !! thermo step from this and the driver's resume fold re-applies
         !! it, so a resume does not cold-start the ice-SW contribution for
         !! one thermo window.  Distinct from the per-category `sw_thru`
         !! above (SCRATCH, recomputed, `delete`-mapped).
      real(wp), allocatable :: m_melt_diag(:, :)
         !! Net meltwater mass to the ocean in the last window (kg/m^2 of
         !! cell, = h2o_ice_to_ocn - h2o_ocn_to_ice summed over cats;
         !! POSITIVE = net melt). OVERWRITTEN each window. Restart-carried
         !! for post-resume diagnostic continuity (like m_frozen_diag).
      real(wp), allocatable :: fprec_ocn_diag(:, :)
         !! PR 26: frozen precipitation delivered directly to the ocean
         !! (kg/m^2/s, per unit CELL area, >= 0) — the share of
         !! `atm_fprec` that lands where there is no ice
         !! (`snow_part_ocn`-weighted; open water at ncat>1, ice-free cells
         !! at ncat==1). Zero on land and on fully ice-covered cells.
         !! SCRATCH: zeroed + rewritten every thermo window by
         !! `ice_snowfall_ocean_share`, NOT restart-carried. Binding seam
         !! spec (consumed by PR-16 as a `net_massin` source): see
         !! `rdb_ice_snow`'s module docstring.

      ! ---- PR 4a: multi-category ITD ----
      real(wp), allocatable :: h_lim(:)
         !! Category lower thickness limits (m), shape (1:ncat+1) — SIS2
         !! `cat_thick_lim` (`SIS_state_initialization.F90:45-78`).
         !! `h_lim(c)` is the lower bound of category `c`; `h_lim(ncat+1)`
         !! is stored but never used as an upper cap (category `ncat` is
         !! unbounded above — SIS2 keeps it too). Computed at `init` by
         !! `ice_itd_category_bounds` from `ncat` and (PR-58) the optional
         !! `hlim_cfg` override below — NOT
         !! restart-carried (cheap to recompute, and doing so avoids a
         !! stale-`ncat` mismatch after a restart that changes `ncat`).
      real(wp), allocatable :: mh_lim(:)
         !! Category lower MASS limits (kg/m²), shape (1:ncat+1),
         !! `= ICE_RHO_ICE*h_lim` — SIS2 `mH_cat_bound`
         !! (`SIS_state_initialization.F90:75-77`). `ice_adjust_categories`
         !! compares `m_ice` against THIS array (mass space, like SIS2),
         !! never against `h_lim` directly. Same non-restart-carried
         !! lifecycle as `h_lim`.
      real(wp), allocatable :: hlim_cfg(:)
         !! PR-58: HOST-ONLY config latch: the user's `&ocean_ice_nml
         !! hlim` list, trimmed to its supplied length, set by
         !! `ocean_state_init_from_config` BEFORE `init`. Unallocated =>
         !! the SIS2 default table. NOT device-mapped (no kernel reads
         !! it — `h_lim`/`mh_lim` are the mapped artefacts, computed at
         !! `init` from this before `enter_data`) and NOT
         !! restart-carried (same lifecycle as `h_lim`).
      real(wp), allocatable :: fb_part_sum(:, :)
         !! Pre-column fb-charged ice-cover fraction snapshot (shape
         !! (nx_total, ny_total)) — `Σ_c part_size(c)` over the categories
         !! that pass the column's OWN entry gate
         !! (`m_ice(i,j,c) > ICE_RHO_ICE*H_VANISHED`), evaluated
         !! POST-uptake / PRE-column by `rdb_ice_thermo_driver` so the
         !! weight matches exactly the categories the column charges `fb`
         !! to before `ice_thermo_columns` mutates `m_ice`. SCRATCH
         !! (recomputed every thermo step at ncat>1 only) — zero-init,
         !! NOT restart-carried, same lifecycle as the `atm_*`/`fb`
         !! seam fields.
      real(wp), allocatable :: snow_part_ocn(:, :)
         !! PR 26: pre-column ice-free-cover fraction snapshot (shape
         !! (nx_total, ny_total)), `1 - Σ_c part_size(c)` over the
         !! categories that pass the column's OWN entry gate
         !! (`m_ice(i,j,c) > ICE_RHO_ICE*H_VANISHED`), evaluated PRE-column
         !! by `ice_snow_part_ocn_fill_impl` for the same reason as
         !! `fb_part_sum` — the weight must match the cover snow actually
         !! caught, before `ice_thermo_columns` mutates `m_ice`. At ncat==1
         !! this is 0 or 1 (binary concentration). SCRATCH (recomputed
         !! every thermo step, only when `has_snowfall`) — zero-init, NOT
         !! restart-carried, same lifecycle as `fb_part_sum`.

      ! ---- Cached extents ----
      integer :: nx_total = 0
      integer :: ny_total = 0
   contains
      procedure, non_overridable :: init => ocean_sea_ice_init
      procedure, non_overridable :: destroy => ocean_sea_ice_destroy
      procedure, non_overridable :: enter_data => ocean_sea_ice_enter_data
      procedure, non_overridable :: exit_data => ocean_sea_ice_exit_data
      procedure, non_overridable :: bytes => ocean_sea_ice_bytes
   end type ocean_sea_ice_t

contains

   subroutine ocean_sea_ice_init(this, grid)
      !! Cache the grid extents, allocate the frazil (PR 1) + Winton
      !! column (PR 3a) prognostics, and mark the slot live.  Allocation
      !! is gated on `enable` at the parent call site (memory Rule 2),
      !! so an ice-off run never carries these arrays.
      class(ocean_sea_ice_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid

      if (this%nk_ice > ICE_NK_MAX .or. this%nk_ice < 1) then
         error stop "rdb_ice_state: nk_ice must be in [1, ICE_NK_MAX]"
      end if
      if (this%ncat < 1) then
         error stop "rdb_ice_state: ncat must be >= 1"
      end if
      if (ICE_CP_BRINE /= ICE_CP_ICE) then
         error stop "rdb_ice_state: rdb_ice_column assumes ICE_CP_BRINE == ICE_CP_ICE " &
            //"(the Newton/false-position branches are not ported)"
      end if
      ! PR-58: defence-in-depth — `init` is reachable from tests (and any
      ! other caller) without going through `validate_config`, which
      ! already polices `hlim_cfg`'s shape at configure time.
      if (allocated(this%hlim_cfg)) then
         if (size(this%hlim_cfg) < 2 .or. size(this%hlim_cfg) > this%ncat + 1) then
            error stop "rdb_ice_state: hlim_cfg must hold 2..ncat+1 entries"
         end if
      end if

      this%nx_total = grid%nx_total
      this%ny_total = grid%ny_total
      allocate (this%frazil_heat(this%nx_total, this%ny_total), source=0.0_wp)
      allocate (this%heat_budget_frazil(this%nx_total, this%ny_total, 1), source=0.0_wp)

      allocate (this%part_size(this%nx_total, this%ny_total, 0:this%ncat), source=0.0_wp)
      this%part_size(:, :, 0) = 1.0_wp
      allocate (this%m_ice(this%nx_total, this%ny_total, this%ncat), source=0.0_wp)
      allocate (this%m_snow(this%nx_total, this%ny_total, this%ncat), source=0.0_wp)
      allocate (this%enth_ice(this%nx_total, this%ny_total, this%ncat, this%nk_ice), &
                source=0.0_wp)
      allocate (this%enth_snow(this%nx_total, this%ny_total, this%ncat, 1), source=0.0_wp)
      allocate (this%sal_ice(this%nx_total, this%ny_total, this%ncat, this%nk_ice), &
                source=ICE_BULK_SALINITY)

      allocate (this%m_frozen_diag(this%nx_total, this%ny_total), source=0.0_wp)
      allocate (this%salt_flux_diag(this%nx_total, this%ny_total), source=0.0_wp)

      allocate (this%atm_sf0(this%nx_total, this%ny_total), source=0.0_wp)
      allocate (this%atm_dsfdt(this%nx_total, this%ny_total), source=0.0_wp)
      allocate (this%atm_sw_dn(this%nx_total, this%ny_total), source=0.0_wp)
      allocate (this%atm_fprec(this%nx_total, this%ny_total), source=0.0_wp)

      allocate (this%fb(this%nx_total, this%ny_total), source=0.0_wp)
      allocate (this%sst_seam(this%nx_total, this%ny_total), source=0.0_wp)
      allocate (this%ssurf_seam(this%nx_total, this%ny_total), source=0.0_wp)
      allocate (this%tfw_seam(this%nx_total, this%ny_total), source=0.0_wp)

      allocate (this%tsurf_out(this%nx_total, this%ny_total, this%ncat), source=0.0_wp)
      allocate (this%h2o_ocn_to_ice(this%nx_total, this%ny_total, this%ncat), source=0.0_wp)
      allocate (this%h2o_ice_to_ocn(this%nx_total, this%ny_total, this%ncat), source=0.0_wp)
      allocate (this%heat_to_ocn(this%nx_total, this%ny_total, this%ncat), source=0.0_wp)
      allocate (this%sw_thru(this%nx_total, this%ny_total, this%ncat), source=0.0_wp)
      allocate (this%snow_to_ice(this%nx_total, this%ny_total, this%ncat), source=0.0_wp)

      allocate (this%heat_flux_diag(this%nx_total, this%ny_total), source=0.0_wp)
      allocate (this%sw_thru_diag(this%nx_total, this%ny_total), source=0.0_wp)
      allocate (this%m_melt_diag(this%nx_total, this%ny_total), source=0.0_wp)
      allocate (this%fprec_ocn_diag(this%nx_total, this%ny_total), source=0.0_wp)

      allocate (this%h_lim(this%ncat + 1), source=0.0_wp)
      allocate (this%mh_lim(this%ncat + 1), source=0.0_wp)
      if (allocated(this%hlim_cfg)) then
         call ice_itd_category_bounds(this%ncat, this%h_lim, this%mh_lim, hlim_vals=this%hlim_cfg)
      else
         call ice_itd_category_bounds(this%ncat, this%h_lim, this%mh_lim)
      end if
      allocate (this%fb_part_sum(this%nx_total, this%ny_total), source=0.0_wp)
      allocate (this%snow_part_ocn(this%nx_total, this%ny_total), source=0.0_wp)

      ! PR 4b: u_ice/v_ice are unconditional (PR 5 EVP needs them
      ! regardless of whether category transport is on).
      allocate (this%u_ice(this%nx_total + 1, this%ny_total), source=0.0_wp)
      allocate (this%v_ice(this%nx_total, this%ny_total + 1), source=0.0_wp)

      ! PR 5: EVP stress + tau fields, unconditional (cheap 2-D arrays;
      ! same allocation contract as u_ice/v_ice — allocated whenever the
      ! ice slot is live, not gated on `dynamics`).
      allocate (this%str_d(this%nx_total, this%ny_total), source=0.0_wp)
      allocate (this%str_t(this%nx_total, this%ny_total), source=0.0_wp)
      allocate (this%str_s(this%nx_total + 1, this%ny_total + 1), source=0.0_wp)
      allocate (this%tau_a_x(this%nx_total + 1, this%ny_total), source=0.0_wp)
      allocate (this%tau_a_y(this%nx_total, this%ny_total + 1), source=0.0_wp)
      allocate (this%fxoc(this%nx_total + 1, this%ny_total), source=0.0_wp)
      allocate (this%fyoc(this%nx_total, this%ny_total + 1), source=0.0_wp)

      ! PR 63: tau_ocn_x/y are unconditional (cheap 2-D arrays, same
      ! contract as fxoc/fyoc). tau_ocn_valid needs no allocation — its
      ! default 0.0_wp on the type already means "no blend written yet".
      allocate (this%tau_ocn_x(this%nx_total + 1, this%ny_total), source=0.0_wp)
      allocate (this%tau_ocn_y(this%nx_total, this%ny_total + 1), source=0.0_wp)

      ! PR 5: EVP subcycle scratch, gated on `this%dynamics` (the workspace
      ! is dynamics-only — same as the retired module `evp_workspace_ensure`
      ! gating).  Rides `ocean_state_enter_data` via this%evp_ws%enter_data.
      if (this%dynamics) then
         call this%evp_ws%init(this%nx_total, this%ny_total)
      end if

      ! PR 4b: transport workspace, gated on `this%transport` (memory
      ! Rule 2 — do not pay the per-category face-flux footprint by
      ! default).
      if (this%transport) then
         allocate (this%mca_ice(this%nx_total, this%ny_total, this%ncat), source=0.0_wp)
         allocate (this%mca_snow(this%nx_total, this%ny_total, this%ncat), source=0.0_wp)
         allocate (this%uh_ice(this%nx_total + 1, this%ny_total, this%ncat), source=0.0_wp)
         allocate (this%vh_ice(this%nx_total, this%ny_total + 1, this%ncat), source=0.0_wp)
         allocate (this%uh_snow(this%nx_total + 1, this%ny_total, this%ncat), source=0.0_wp)
         allocate (this%vh_snow(this%nx_total, this%ny_total + 1, this%ncat), source=0.0_wp)
         allocate (this%htot_work(this%nx_total, this%ny_total), source=0.0_wp)
         allocate (this%hl_x_work(this%nx_total + 1, this%ny_total), source=0.0_wp)
         allocate (this%hr_x_work(this%nx_total + 1, this%ny_total), source=0.0_wp)
         allocate (this%hl_y_work(this%nx_total, this%ny_total + 1), source=0.0_wp)
         allocate (this%hr_y_work(this%nx_total, this%ny_total + 1), source=0.0_wp)
         allocate (this%uhtot_work(this%nx_total + 1, this%ny_total), source=0.0_wp)
         allocate (this%vhtot_work(this%nx_total, this%ny_total + 1), source=0.0_wp)
         allocate (this%tr_flux_x_work(this%nx_total + 1, this%ny_total, this%ncat), source=0.0_wp)
         allocate (this%tr_flux_y_work(this%nx_total, this%ny_total + 1, this%ncat), source=0.0_wp)
      end if

      this%is_init = .true.
   end subroutine ocean_sea_ice_init

   subroutine ocean_sea_ice_destroy(this)
      !! Reverse of `init`.  Clears `is_init` first (use-after-destroy
      !! guard), then releases allocations.
      class(ocean_sea_ice_t), intent(inout) :: this
      this%is_init = .false.
      if (allocated(this%frazil_heat)) deallocate (this%frazil_heat)
      if (allocated(this%heat_budget_frazil)) deallocate (this%heat_budget_frazil)
      if (allocated(this%part_size)) deallocate (this%part_size)
      if (allocated(this%m_ice)) deallocate (this%m_ice)
      if (allocated(this%m_snow)) deallocate (this%m_snow)
      if (allocated(this%enth_ice)) deallocate (this%enth_ice)
      if (allocated(this%enth_snow)) deallocate (this%enth_snow)
      if (allocated(this%sal_ice)) deallocate (this%sal_ice)
      if (allocated(this%m_frozen_diag)) deallocate (this%m_frozen_diag)
      if (allocated(this%salt_flux_diag)) deallocate (this%salt_flux_diag)
      if (allocated(this%atm_sf0)) deallocate (this%atm_sf0)
      if (allocated(this%atm_dsfdt)) deallocate (this%atm_dsfdt)
      if (allocated(this%atm_sw_dn)) deallocate (this%atm_sw_dn)
      if (allocated(this%atm_fprec)) deallocate (this%atm_fprec)
      if (allocated(this%fb)) deallocate (this%fb)
      if (allocated(this%sst_seam)) deallocate (this%sst_seam)
      if (allocated(this%ssurf_seam)) deallocate (this%ssurf_seam)
      if (allocated(this%tfw_seam)) deallocate (this%tfw_seam)
      if (allocated(this%tsurf_out)) deallocate (this%tsurf_out)
      if (allocated(this%h2o_ocn_to_ice)) deallocate (this%h2o_ocn_to_ice)
      if (allocated(this%h2o_ice_to_ocn)) deallocate (this%h2o_ice_to_ocn)
      if (allocated(this%heat_to_ocn)) deallocate (this%heat_to_ocn)
      if (allocated(this%sw_thru)) deallocate (this%sw_thru)
      if (allocated(this%snow_to_ice)) deallocate (this%snow_to_ice)
      if (allocated(this%heat_flux_diag)) deallocate (this%heat_flux_diag)
      if (allocated(this%sw_thru_diag)) deallocate (this%sw_thru_diag)
      if (allocated(this%m_melt_diag)) deallocate (this%m_melt_diag)
      if (allocated(this%fprec_ocn_diag)) deallocate (this%fprec_ocn_diag)
      if (allocated(this%h_lim)) deallocate (this%h_lim)
      if (allocated(this%mh_lim)) deallocate (this%mh_lim)
      if (allocated(this%hlim_cfg)) deallocate (this%hlim_cfg)
      if (allocated(this%fb_part_sum)) deallocate (this%fb_part_sum)
      if (allocated(this%snow_part_ocn)) deallocate (this%snow_part_ocn)
      if (allocated(this%u_ice)) deallocate (this%u_ice)
      if (allocated(this%v_ice)) deallocate (this%v_ice)
      if (allocated(this%str_d)) deallocate (this%str_d)
      if (allocated(this%str_t)) deallocate (this%str_t)
      if (allocated(this%str_s)) deallocate (this%str_s)
      if (allocated(this%tau_a_x)) deallocate (this%tau_a_x)
      if (allocated(this%tau_a_y)) deallocate (this%tau_a_y)
      if (allocated(this%fxoc)) deallocate (this%fxoc)
      if (allocated(this%fyoc)) deallocate (this%fyoc)
      if (allocated(this%tau_ocn_x)) deallocate (this%tau_ocn_x)
      if (allocated(this%tau_ocn_y)) deallocate (this%tau_ocn_y)
      this%tau_ocn_valid = 0.0_wp
      if (this%dynamics) call this%evp_ws%destroy()
      if (allocated(this%mca_ice)) deallocate (this%mca_ice)
      if (allocated(this%mca_snow)) deallocate (this%mca_snow)
      if (allocated(this%uh_ice)) deallocate (this%uh_ice)
      if (allocated(this%vh_ice)) deallocate (this%vh_ice)
      if (allocated(this%uh_snow)) deallocate (this%uh_snow)
      if (allocated(this%vh_snow)) deallocate (this%vh_snow)
      if (allocated(this%htot_work)) deallocate (this%htot_work)
      if (allocated(this%hl_x_work)) deallocate (this%hl_x_work)
      if (allocated(this%hr_x_work)) deallocate (this%hr_x_work)
      if (allocated(this%hl_y_work)) deallocate (this%hl_y_work)
      if (allocated(this%hr_y_work)) deallocate (this%hr_y_work)
      if (allocated(this%uhtot_work)) deallocate (this%uhtot_work)
      if (allocated(this%vhtot_work)) deallocate (this%vhtot_work)
      if (allocated(this%tr_flux_x_work)) deallocate (this%tr_flux_x_work)
      if (allocated(this%tr_flux_y_work)) deallocate (this%tr_flux_y_work)
      this%nx_total = 0
      this%ny_total = 0
   end subroutine ocean_sea_ice_destroy

   subroutine ocean_sea_ice_enter_data(this)
      ! Poly TBP delegating to a `type(...)`-arg `_impl` (AMD-crash rule).
      class(ocean_sea_ice_t), intent(inout) :: this
      select type (this)
      type is (ocean_sea_ice_t)
         call ocean_sea_ice_enter_data_impl(this)
      end select
   end subroutine ocean_sea_ice_enter_data

   subroutine ocean_sea_ice_enter_data_impl(this)
      !! Attach the slot's allocatables to the device.  `copyin` (not
      !! `create`) throughout: every array is host-initialised (0, or
      !! `part_size`/`sal_ice`'s non-zero defaults) at init, and a
      !! restart read seeds the persistent ones host-side before the
      !! mapping.  PR 3c's seam + scratch fields are also `copyin`: they
      !! are host-zeroed at init and (re)filled device-side every thermo
      !! step, so a plain `copyin` is correct (and cheap — one-time).
      !! PR 4b: `u_ice`/`v_ice` are `copyin` (host-zeroed at init, and the
      !! v1 sampler / PR-5 EVP both WRITE them device-side every call —
      !! `copyin` establishes presence, same contract as the `atm_*` seam).
      !! PR 5: `str_d`/`str_t`/`str_s` are `copyin` (host-zeroed at init or
      !! restart-seeded; `ice_evp_dynamics` reads them before the first
      !! write on every call — `limit_stresses` runs first). `tau_a_x`/
      !! `tau_a_y` are `copyin` (host-snapshotted at configure, then read
      !! device-side every substep — never written on-device). `fxoc`/
      !! `fyoc` are `copyin` (host-zeroed at init or restart-seeded; the
      !! EVP core zeros them itself at call entry before accumulating).
      !! PR 63: `tau_ocn_x`/`tau_ocn_y` are `copyin` (host-zeroed at init
      !! or restart-seeded; `ice_tau_mirror_impl` WRITES them device-side
      !! every call to `ice_ocean_stress_flux` — same reasoning as
      !! `fxoc`/`fyoc`). `tau_ocn_valid` is NEVER mapped — it is host-only
      !! scalar state, set on the host by `ice_ocean_stress_flux` and read
      !! on the host by `ice_ocean_stress_resume_apply`; `register_scalar`
      !! forces `device_mapped=.false.` on the restart side for the same
      !! reason.
      !! The transport workspace is `create` (pure scratch, recomputed
      !! from scratch every pass — never read before written), mapped
      !! only `if (this%transport)` (memory Rule 2, same gate as the
      !! host-side allocation).
      type(ocean_sea_ice_t), intent(inout) :: this
      if (.not. this%is_init) return
      !$acc enter data copyin(this%frazil_heat, this%heat_budget_frazil)
      !$acc enter data copyin(this%part_size, this%m_ice, this%m_snow)
      !$acc enter data copyin(this%enth_ice, this%enth_snow, this%sal_ice)
      !$acc enter data copyin(this%m_frozen_diag, this%salt_flux_diag)
      !$acc enter data copyin(this%atm_sf0, this%atm_dsfdt, this%atm_sw_dn, this%atm_fprec)
      !$acc enter data copyin(this%fb, this%sst_seam, this%ssurf_seam, this%tfw_seam)
      !$acc enter data copyin(this%tsurf_out, this%h2o_ocn_to_ice, this%h2o_ice_to_ocn)
      !$acc enter data copyin(this%heat_to_ocn, this%sw_thru, this%snow_to_ice)
      !$acc enter data copyin(this%heat_flux_diag, this%sw_thru_diag, this%m_melt_diag, &
      !$acc&                   this%fprec_ocn_diag)
      !$acc enter data copyin(this%h_lim, this%mh_lim, this%fb_part_sum, this%snow_part_ocn)
      !$acc enter data copyin(this%u_ice, this%v_ice)
      !$acc enter data copyin(this%str_d, this%str_t, this%str_s)
      !$acc enter data copyin(this%tau_a_x, this%tau_a_y)
      !$acc enter data copyin(this%fxoc, this%fyoc)
      !$acc enter data copyin(this%tau_ocn_x, this%tau_ocn_y)
      ! tau_ocn_valid: NEVER mapped (host-only, register_scalar contract).
      if (this%dynamics) call this%evp_ws%enter_data()
      if (this%transport) then
         !$acc enter data create(this%mca_ice, this%mca_snow)
         !$acc enter data create(this%uh_ice, this%vh_ice, this%uh_snow, this%vh_snow)
         !$acc enter data create(this%htot_work, this%hl_x_work, this%hr_x_work)
         !$acc enter data create(this%hl_y_work, this%hr_y_work)
         !$acc enter data create(this%uhtot_work, this%vhtot_work)
         !$acc enter data create(this%tr_flux_x_work, this%tr_flux_y_work)
      end if
   end subroutine ocean_sea_ice_enter_data_impl

   subroutine ocean_sea_ice_exit_data(this)
      class(ocean_sea_ice_t), intent(inout) :: this
      select type (this)
      type is (ocean_sea_ice_t)
         call ocean_sea_ice_exit_data_impl(this)
      end select
   end subroutine ocean_sea_ice_exit_data

   subroutine ocean_sea_ice_exit_data_impl(this)
      !! Reverse of `enter_data_impl` — all six PR-3a fields are
      !! prognostic state (post-run host inspection + restart write),
      !! same as the persistent frazil bank; the drain-and-zero budget
      !! scratch is dropped. PR 3b's `m_frozen_diag`/`salt_flux_diag`
      !! are both PERSISTENT (restart-carried) — `copyout`, not `delete`.
      !! PR 3c: `heat_flux_diag`/`m_melt_diag` are likewise PERSISTENT
      !! (`copyout`); PR 31's `sw_thru_diag` is PERSISTENT too (`copyout`,
      !! restart-carried — NOT its per-category parent `sw_thru`, which
      !! stays `delete`); the `atm_*` seam, `fb`/`sst_seam`/`ssurf_seam`/
      !! `tfw_seam`, and the per-category column scratch
      !! (`tsurf_out`/`h2o_ocn_to_ice`/`h2o_ice_to_ocn`/`heat_to_ocn`/
      !! `sw_thru`) are all recomputed every thermo step — `delete`. PR 4a:
      !! `h_lim`/`mh_lim` are `delete` (recomputed from `ncat` at `init`,
      !! never mutated after); `fb_part_sum` is `delete` (recomputed every
      !! thermo step, same lifecycle as the other scratch seams). PR 4b:
      !! `u_ice`/`v_ice` are `copyout` (post-run host inspection, same
      !! contract as the melt/frazil diags); the transport workspace is
      !! `delete` (pure scratch), gated the same `if (this%transport)` as
      !! the enter_data map.  PR 5: `str_d`/`str_t`/`str_s`/`fxoc`/`fyoc`
      !! are `copyout` (PERSISTENT, restart-carried); `tau_a_x`/`tau_a_y`
      !! are `delete` (configure-time snapshot, not restart-carried —
      !! rebuilt fresh on every run). PR 26: `atm_fprec` is `delete`
      !! (same lifecycle as the rest of the `atm_*` seam); `fprec_ocn_diag`/
      !! `snow_part_ocn` are `delete` (SCRATCH, recomputed every thermo
      !! step, same lifecycle as `fb_part_sum`). PR 63: `tau_ocn_x`/`tau_ocn_y`
      !! are `copyout` — PERSISTENT, restart-carried, and the restart write
      !! reads the host copy (opposite lifecycle to `tau_a_x`/`tau_a_y`,
      !! do not copy that pattern). `tau_ocn_valid` is never mapped, so
      !! there is nothing to unmap.
      type(ocean_sea_ice_t), intent(inout) :: this
      if (.not. this%is_init) return
      !$acc exit data copyout(this%frazil_heat)
      !$acc exit data delete(this%heat_budget_frazil)
      !$acc exit data copyout(this%part_size, this%m_ice, this%m_snow)
      !$acc exit data copyout(this%enth_ice, this%enth_snow, this%sal_ice)
      !$acc exit data copyout(this%m_frozen_diag, this%salt_flux_diag)
      !$acc exit data delete(this%atm_sf0, this%atm_dsfdt, this%atm_sw_dn, this%atm_fprec)
      !$acc exit data delete(this%fb, this%sst_seam, this%ssurf_seam, this%tfw_seam)
      !$acc exit data delete(this%tsurf_out, this%h2o_ocn_to_ice, this%h2o_ice_to_ocn)
      !$acc exit data delete(this%heat_to_ocn, this%sw_thru, this%snow_to_ice)
      !$acc exit data copyout(this%heat_flux_diag, this%sw_thru_diag, this%m_melt_diag)
      !$acc exit data delete(this%fprec_ocn_diag)
      !$acc exit data delete(this%h_lim, this%mh_lim, this%fb_part_sum, this%snow_part_ocn)
      !$acc exit data copyout(this%u_ice, this%v_ice)
      !$acc exit data copyout(this%str_d, this%str_t, this%str_s)
      !$acc exit data delete(this%tau_a_x, this%tau_a_y)
      !$acc exit data copyout(this%fxoc, this%fyoc)
      !$acc exit data copyout(this%tau_ocn_x, this%tau_ocn_y)
      if (this%dynamics) call this%evp_ws%exit_data()
      if (this%transport) then
         !$acc exit data delete(this%mca_ice, this%mca_snow)
         !$acc exit data delete(this%uh_ice, this%vh_ice, this%uh_snow, this%vh_snow)
         !$acc exit data delete(this%htot_work, this%hl_x_work, this%hr_x_work)
         !$acc exit data delete(this%hl_y_work, this%hr_y_work)
         !$acc exit data delete(this%uhtot_work, this%vhtot_work)
         !$acc exit data delete(this%tr_flux_x_work, this%tr_flux_y_work)
      end if
   end subroutine ocean_sea_ice_exit_data_impl

   pure function ocean_sea_ice_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the sea-ice slot (one
      !! `arr_bytes` term per array, mirroring `ocean_meke_bytes`).
      class(ocean_sea_ice_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = 0_int64
      if (.not. this%is_init) return
      nbytes = arr_bytes(this%frazil_heat) &
               + arr_bytes(this%heat_budget_frazil) &
               + arr_bytes(this%part_size) &
               + arr_bytes(this%m_ice) &
               + arr_bytes(this%m_snow) &
               + arr_bytes(this%enth_ice) &
               + arr_bytes(this%enth_snow) &
               + arr_bytes(this%sal_ice) &
               + arr_bytes(this%m_frozen_diag) &
               + arr_bytes(this%salt_flux_diag) &
               + arr_bytes(this%atm_sf0) &
               + arr_bytes(this%atm_dsfdt) &
               + arr_bytes(this%atm_sw_dn) &
               + arr_bytes(this%atm_fprec) &
               + arr_bytes(this%fb) &
               + arr_bytes(this%sst_seam) &
               + arr_bytes(this%ssurf_seam) &
               + arr_bytes(this%tfw_seam) &
               + arr_bytes(this%tsurf_out) &
               + arr_bytes(this%h2o_ocn_to_ice) &
               + arr_bytes(this%h2o_ice_to_ocn) &
               + arr_bytes(this%heat_to_ocn) &
               + arr_bytes(this%sw_thru) &
               + arr_bytes(this%snow_to_ice) &
               + arr_bytes(this%heat_flux_diag) &
               + arr_bytes(this%sw_thru_diag) &
               + arr_bytes(this%m_melt_diag) &
               + arr_bytes(this%fprec_ocn_diag) &
               + arr_bytes(this%h_lim) &
               + arr_bytes(this%mh_lim) &
               + arr_bytes(this%hlim_cfg) &
               + arr_bytes(this%fb_part_sum) &
               + arr_bytes(this%snow_part_ocn) &
               + arr_bytes(this%u_ice) &
               + arr_bytes(this%v_ice) &
               + arr_bytes(this%str_d) &
               + arr_bytes(this%str_t) &
               + arr_bytes(this%str_s) &
               + arr_bytes(this%tau_a_x) &
               + arr_bytes(this%tau_a_y) &
               + arr_bytes(this%fxoc) &
               + arr_bytes(this%fyoc) &
               + arr_bytes(this%tau_ocn_x) &
               + arr_bytes(this%tau_ocn_y) &
               + arr_bytes(this%mca_ice) &
               + arr_bytes(this%mca_snow) &
               + arr_bytes(this%uh_ice) &
               + arr_bytes(this%vh_ice) &
               + arr_bytes(this%uh_snow) &
               + arr_bytes(this%vh_snow) &
               + arr_bytes(this%htot_work) &
               + arr_bytes(this%hl_x_work) &
               + arr_bytes(this%hr_x_work) &
               + arr_bytes(this%hl_y_work) &
               + arr_bytes(this%hr_y_work) &
               + arr_bytes(this%uhtot_work) &
               + arr_bytes(this%vhtot_work) &
               + arr_bytes(this%tr_flux_x_work) &
               + arr_bytes(this%tr_flux_y_work)
      ! EVP rheology scratch — allocated + device-mapped only when
      ! `dynamics` is on, so the term self-zeroes on the EVP-off path.
      nbytes = nbytes + this%evp_ws%bytes()
   end function ocean_sea_ice_bytes

   ! =====================================================================
   ! EVP subcycle workspace lifecycle (successor to the retired
   ! `rdb_ice_evp` module-level `save` scratch + `evp_workspace_ensure`).
   ! =====================================================================

   pure function evp_workspace_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the EVP scratch (0 when
      !! unallocated, i.e. whenever `&ocean_ice_nml dynamics` is off).
      !!
      !! Every array here is device-mapped by `evp_workspace_enter_data_impl`
      !! and was the single omission in `ocean_sea_ice_bytes` (~154 MB at
      !! 1000x800).  One `arr_bytes` term per array — add one here when a
      !! new allocatable joins the type.
      class(evp_workspace_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = arr_bytes(this%mis_w) &
               + arr_bytes(this%mice_w) &
               + arr_bytes(this%ci_w) &
               + arr_bytes(this%mis_in_w) &
               + arr_bytes(this%mice_in_w) &
               + arr_bytes(this%ci_in_w) &
               + arr_bytes(this%pres_mice_w) &
               + arr_bytes(this%del_sh_min_pr_w) &
               + arr_bytes(this%sh_dd_w) &
               + arr_bytes(this%sh_dt_w) &
               + arr_bytes(this%zeta_w) &
               + arr_bytes(this%del_sh_w) &
               + arr_bytes(this%mask_t_w) &
               + arr_bytes(this%mi_u_w) &
               + arr_bytes(this%mask_u_w) &
               + arr_bytes(this%u_tmp_w) &
               + arr_bytes(this%mi_v_w) &
               + arr_bytes(this%mask_v_w) &
               + arr_bytes(this%a_u_w) &
               + arr_bytes(this%a_v_w) &
               + arr_bytes(this%sh_ds_w) &
               + arr_bytes(this%mi_ratio_a_q_w) &
               + arr_bytes(this%q_w) &
               + arr_bytes(this%mask_q_w)
   end function evp_workspace_bytes

   subroutine evp_workspace_init(this, nx, ny)
      !! Allocate the EVP scratch for an `(nx,ny)` T-cell domain, zero-init.
      !! Called ONCE per workspace lifetime — from `ocean_sea_ice_init`
      !! (gated on `dynamics`) and from the EVP test harness.  `nx`/`ny`
      !! declared before the arrays that use them is unnecessary here (all
      !! allocatables), but the extents are cached for `destroy`.
      class(evp_workspace_t), intent(inout) :: this
      integer, intent(in) :: nx, ny

      this%nx = nx
      this%ny = ny

      allocate (this%mis_w(nx, ny), source=0.0_wp)
      allocate (this%mice_w(nx, ny), source=0.0_wp)
      allocate (this%ci_w(nx, ny), source=0.0_wp)
      allocate (this%mis_in_w(nx, ny), source=0.0_wp)
      allocate (this%mice_in_w(nx, ny), source=0.0_wp)
      allocate (this%ci_in_w(nx, ny), source=0.0_wp)
      allocate (this%pres_mice_w(nx, ny), source=0.0_wp)
      allocate (this%del_sh_min_pr_w(nx, ny), source=0.0_wp)
      allocate (this%sh_dd_w(nx, ny), source=0.0_wp)
      allocate (this%sh_dt_w(nx, ny), source=0.0_wp)
      allocate (this%zeta_w(nx, ny), source=0.0_wp)
      allocate (this%del_sh_w(nx, ny), source=0.0_wp)
      allocate (this%mask_t_w(nx, ny), source=0.0_wp)

      allocate (this%mi_u_w(nx + 1, ny), source=0.0_wp)
      allocate (this%mask_u_w(nx + 1, ny), source=0.0_wp)
      allocate (this%u_tmp_w(nx + 1, ny), source=0.0_wp)
      allocate (this%mi_v_w(nx, ny + 1), source=0.0_wp)
      allocate (this%mask_v_w(nx, ny + 1), source=0.0_wp)
      allocate (this%a_u_w(nx + 1, ny), source=0.0_wp)
      allocate (this%a_v_w(nx, ny + 1), source=0.0_wp)

      allocate (this%sh_ds_w(nx + 1, ny + 1), source=0.0_wp)
      allocate (this%mi_ratio_a_q_w(nx + 1, ny + 1), source=0.0_wp)
      allocate (this%q_w(nx + 1, ny + 1), source=0.0_wp)
      allocate (this%mask_q_w(nx + 1, ny + 1), source=0.0_wp)
   end subroutine evp_workspace_init

   subroutine evp_workspace_enter_data(this)
      ! Poly TBP delegating to a `type(...)`-arg `_impl` (AMD-crash rule,
      ! mirrors `ocean_sea_ice_enter_data`).
      class(evp_workspace_t), intent(inout) :: this
      select type (this)
      type is (evp_workspace_t)
         call evp_workspace_enter_data_impl(this)
      end select
   end subroutine evp_workspace_enter_data

   subroutine evp_workspace_enter_data_impl(this)
      !! Attach the scratch to the device with `create` — every array is
      !! pure scratch, fully (re)written before it is read on every
      !! `ice_evp_dynamics` call (masks rebuilt, cell fields filled, the
      !! subcycle kernels overwrite the rest), so `create` (not `copyin`)
      !! is correct.  One-level components of the `type(...)` dummy — same
      !! shape as `ocean_sea_ice_enter_data_impl`, no associate-leaf needed.
      type(evp_workspace_t), intent(inout) :: this
      if (.not. allocated(this%mis_w)) return
      !$acc enter data create(this%mis_w, this%mice_w, this%ci_w)
      !$acc enter data create(this%mis_in_w, this%mice_in_w, this%ci_in_w)
      !$acc enter data create(this%pres_mice_w, this%del_sh_min_pr_w)
      !$acc enter data create(this%sh_dd_w, this%sh_dt_w, this%zeta_w, this%del_sh_w)
      !$acc enter data create(this%mask_t_w)
      !$acc enter data create(this%mi_u_w, this%mask_u_w, this%u_tmp_w)
      !$acc enter data create(this%mi_v_w, this%mask_v_w)
      !$acc enter data create(this%a_u_w, this%a_v_w)
      !$acc enter data create(this%sh_ds_w, this%mi_ratio_a_q_w, this%q_w, this%mask_q_w)
   end subroutine evp_workspace_enter_data_impl

   subroutine evp_workspace_exit_data(this)
      class(evp_workspace_t), intent(inout) :: this
      select type (this)
      type is (evp_workspace_t)
         call evp_workspace_exit_data_impl(this)
      end select
   end subroutine evp_workspace_exit_data

   subroutine evp_workspace_exit_data_impl(this)
      !! Reverse of `enter_data_impl` — pure scratch, `delete` throughout.
      type(evp_workspace_t), intent(inout) :: this
      if (.not. allocated(this%mis_w)) return
      !$acc exit data delete(this%mis_w, this%mice_w, this%ci_w)
      !$acc exit data delete(this%mis_in_w, this%mice_in_w, this%ci_in_w)
      !$acc exit data delete(this%pres_mice_w, this%del_sh_min_pr_w)
      !$acc exit data delete(this%sh_dd_w, this%sh_dt_w, this%zeta_w, this%del_sh_w)
      !$acc exit data delete(this%mask_t_w)
      !$acc exit data delete(this%mi_u_w, this%mask_u_w, this%u_tmp_w)
      !$acc exit data delete(this%mi_v_w, this%mask_v_w)
      !$acc exit data delete(this%a_u_w, this%a_v_w)
      !$acc exit data delete(this%sh_ds_w, this%mi_ratio_a_q_w, this%q_w, this%mask_q_w)
   end subroutine evp_workspace_exit_data_impl

   subroutine evp_workspace_destroy(this)
      !! Reverse of `init` (host deallocation).  Idempotent — safe on an
      !! already-clean workspace.  Device unmapping is `exit_data`'s job
      !! (call it first).
      class(evp_workspace_t), intent(inout) :: this
      if (allocated(this%mis_w)) deallocate (this%mis_w)
      if (allocated(this%mice_w)) deallocate (this%mice_w)
      if (allocated(this%ci_w)) deallocate (this%ci_w)
      if (allocated(this%mis_in_w)) deallocate (this%mis_in_w)
      if (allocated(this%mice_in_w)) deallocate (this%mice_in_w)
      if (allocated(this%ci_in_w)) deallocate (this%ci_in_w)
      if (allocated(this%pres_mice_w)) deallocate (this%pres_mice_w)
      if (allocated(this%del_sh_min_pr_w)) deallocate (this%del_sh_min_pr_w)
      if (allocated(this%sh_dd_w)) deallocate (this%sh_dd_w)
      if (allocated(this%sh_dt_w)) deallocate (this%sh_dt_w)
      if (allocated(this%zeta_w)) deallocate (this%zeta_w)
      if (allocated(this%del_sh_w)) deallocate (this%del_sh_w)
      if (allocated(this%mask_t_w)) deallocate (this%mask_t_w)
      if (allocated(this%mi_u_w)) deallocate (this%mi_u_w)
      if (allocated(this%mask_u_w)) deallocate (this%mask_u_w)
      if (allocated(this%u_tmp_w)) deallocate (this%u_tmp_w)
      if (allocated(this%mi_v_w)) deallocate (this%mi_v_w)
      if (allocated(this%mask_v_w)) deallocate (this%mask_v_w)
      if (allocated(this%a_u_w)) deallocate (this%a_u_w)
      if (allocated(this%a_v_w)) deallocate (this%a_v_w)
      if (allocated(this%sh_ds_w)) deallocate (this%sh_ds_w)
      if (allocated(this%mi_ratio_a_q_w)) deallocate (this%mi_ratio_a_q_w)
      if (allocated(this%q_w)) deallocate (this%q_w)
      if (allocated(this%mask_q_w)) deallocate (this%mask_q_w)
      this%nx = 0
      this%ny = 0
   end subroutine evp_workspace_destroy

   pure subroutine ice_itd_category_bounds(ncat, h_lim, mh_lim, hlim_vals)
      !! SIS2 `initialize_ice_categories`
      !! (`SIS_state_initialization.F90:45-78`): absent `hlim_vals` fills
      !! the first `min(ncat+1, 8)` entries from the default
      !! lower-thickness-limit table `HLIM_DFLT_TABLE`; a supplied
      !! `hlim_vals` fills the first `min(ncat+1, size(hlim_vals))`
      !! entries from IT instead (PR-58 — the SIS2 `hLim_vals` optional
      !! dummy Roundabout's port originally declined to carry). Either way,
      !! the remainder extrapolates by constant width,
      !! `h_lim(k) = 2*h_lim(k-1) - h_lim(k-2)`, resuming at ONE PAST
      !! however many entries were actually supplied — not at a fixed
      !! index 9. `mh_lim = ICE_RHO_ICE*h_lim` (SIS2 `mH_cat_bound`,
      !! `SIS_state_initialization.F90:75-77`). `h_lim(c)`/`mh_lim(c)` is
      !! the LOWER bound of category `c` (1-based); index `ncat+1` is the
      !! lower bound that WOULD start category `ncat+1` — stored (SIS2
      !! keeps it) but never used as an upper cap on category `ncat`
      !! (`ice_adjust_categories`'s upward pass stops at `c = ncat-1`).
      !! No fixed-size local arrays — `ncat` is unbounded here (unlike
      !! the per-column `ICE_NK_MAX`-capped layer arrays).
      !!
      !! Thickness-distribution discretization: Thorndike, Rothrock,
      !! Maykut & Colony (1975), JGR 80, 4501-4513; Bitz, Holland, Weaver
      !! & Eby (2001), JGR 106, 2441-2463.
      !!
      !! Documented divergences from SIS2's `initialize_ice_categories`
      !! (PR-58 plan §4):
      !!   D1 — SIS2 silently falls back to the default table for a
      !!        1-element `hlim_vals` (`size(hLim_vals) > 1` gate,
      !!        `:60`). Roundabout's caller (`rdb_config%validate_config`)
      !!        FAILS LOUD instead via `ice_hlim_spec_is_valid` (house
      !!        rule: a wordless no-op is worse than an error). This
      !!        routine itself has no opinion — the size floor is
      !!        enforced by the caller, not here.
      !!   D2 — SIS2's namelist CALLER truncates to `CatIce` entries
      !!        (`ice_model.F90:2110`), so its top edge can never be
      !!        stated via the namelist. This routine (and Roundabout's
      !!        `&ocean_ice_nml hlim`) accept the procedure's own native
      !!        range, `2..ncat+1`.
      integer, intent(in) :: ncat
         !! Number of ice thickness categories (declared first —
         !! decl-order).
      real(wp), intent(out) :: h_lim(ncat + 1)
         !! Category lower thickness limits (m), `h_lim(1:ncat+1)`.
      real(wp), intent(out) :: mh_lim(ncat + 1)
         !! Category lower mass limits (kg/m²), `= ICE_RHO_ICE*h_lim`.
      real(wp), intent(in), optional :: hlim_vals(:)
         !! PR-58: optional user-supplied lower edges (m), `2..ncat+1`
         !! entries, strictly increasing, `hlim_vals(1) > 0` — SIS2
         !! `initialize_ice_categories`'s `hLim_vals`
         !! (`SIS_state_initialization.F90:50`). Assumed-shape is correct
         !! here: this is a host-side, once-per-run `init` routine with
         !! no `do concurrent`, not a per-step kernel. Absent => the
         !! pre-PR-58 behaviour, bit-for-bit (the default-table branch
         !! below is unmodified). Validated by the caller
         !! (`ice_hlim_spec_is_valid`) before it ever reaches here — this
         !! routine trusts its shape.

      integer, parameter :: N_HLIM_DFLT = 8
         !! Number of entries in the SIS2 default table (declared before
         !! the array that uses it — decl-order).
      real(wp), parameter :: HLIM_DFLT_TABLE(N_HLIM_DFLT) = &
                             [1.0e-10_wp, 0.1_wp, 0.3_wp, 0.7_wp, 1.1_wp, 1.5_wp, 2.0_wp, 2.5_wp]
      integer :: n_given, k
         !! `n_given`: how many entries were actually filled from the
         !! source (table OR `hlim_vals`) before extrapolation takes
         !! over. Renamed from the pre-PR-58 `n_dflt` (PR-58 plan §11.1)
         !! because after this change it no longer means "how many came
         !! from the table" — it means "how many came from wherever we
         !! got them", and the extrapolation loop below is only correct
         !! if it resumes at `n_given + 1` regardless of source.

      if (present(hlim_vals)) then
         n_given = min(ncat + 1, size(hlim_vals))
         do k = 1, n_given
            h_lim(k) = hlim_vals(k)
         end do
      else
         n_given = min(ncat + 1, N_HLIM_DFLT)
         do k = 1, n_given
            h_lim(k) = HLIM_DFLT_TABLE(k)
         end do
      end if
      do k = n_given + 1, ncat + 1
         h_lim(k) = 2.0_wp*h_lim(k - 1) - h_lim(k - 2)
      end do
      do k = 1, ncat + 1
         mh_lim(k) = ICE_RHO_ICE*h_lim(k)
      end do
   end subroutine ice_itd_category_bounds

   pure subroutine ice_cell_concentration_impl(wet_T, part_size, m_ice, m_snow, &
                                               mis, mice, ci, ncat, nx, ny)
      !! Two-mode per-cell concentration/mass gather (module docstring
      !! convention), shared by `rdb_ice_evp` (`ice_evp_step`) and
      !! `rdb_ice_ocean_coupler` (`ice_ocean_stress_flux`) so neither
      !! module depends on the other.
      !!
      !! `ncat > 1` (ITD live): `mis = Σ_c part_size(c)*(m_ice(c)+
      !! m_snow(c))`, `mice = Σ_c part_size(c)*m_ice(c)`,
      !! `ci = min(1, Σ_{c>=1} part_size(c))`.
      !! `ncat == 1` (legacy lumped; per-CELL masses, D5 — SIS2 has no
      !! lumped mode): `mis = m_ice(1)+m_snow(1)`, `mice = m_ice(1)`,
      !! `ci = merge(1, 0, m_ice(1) > 0)`.
      !! All × `wet_T` (land ⇒ 0).  Decl-order: integer dims before the
      !! explicit-shape arrays that use them.
      integer, intent(in) :: ncat, nx, ny
      real(wp), intent(in) :: wet_T(nx, ny)
      real(wp), intent(in) :: part_size(nx, ny, 0:ncat)
      real(wp), intent(in) :: m_ice(nx, ny, ncat)
      real(wp), intent(in) :: m_snow(nx, ny, ncat)
      real(wp), intent(out) :: mis(nx, ny)
      real(wp), intent(out) :: mice(nx, ny)
      real(wp), intent(out) :: ci(nx, ny)
      integer :: i, j, c
      real(wp) :: mis_sum, mice_sum, ci_sum

      if (ncat == 1) then
         do concurrent(j=1:ny, i=1:nx)
            if (wet_T(i, j) > 0.5_wp .and. m_ice(i, j, 1) > 0.0_wp) then
               mis(i, j) = m_ice(i, j, 1) + m_snow(i, j, 1)
               mice(i, j) = m_ice(i, j, 1)
               ci(i, j) = 1.0_wp
            else
               mis(i, j) = 0.0_wp
               mice(i, j) = 0.0_wp
               ci(i, j) = 0.0_wp
            end if
         end do
      else
         do concurrent(j=1:ny, i=1:nx) local(c, mis_sum, mice_sum, ci_sum)
            if (wet_T(i, j) > 0.5_wp) then
               mis_sum = 0.0_wp
               mice_sum = 0.0_wp
               ci_sum = 0.0_wp
               do c = 1, ncat
                  mis_sum = mis_sum + part_size(i, j, c)*(m_ice(i, j, c) + m_snow(i, j, c))
                  mice_sum = mice_sum + part_size(i, j, c)*m_ice(i, j, c)
                  ci_sum = ci_sum + part_size(i, j, c)
               end do
               mis(i, j) = mis_sum
               mice(i, j) = mice_sum
               ci(i, j) = min(1.0_wp, ci_sum)
            else
               mis(i, j) = 0.0_wp
               mice(i, j) = 0.0_wp
               ci(i, j) = 0.0_wp
            end if
         end do
      end if
   end subroutine ice_cell_concentration_impl

end module rdb_ice_state
