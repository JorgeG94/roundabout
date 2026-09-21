!! Ocean vertical mixing parameterisation state.
module rdb_ocean_vmix
   !! Holds the closure state for the ocean dynamical core's vertical
   !! mixing kernel.  Produces 3D `kv`, `kt` (and `ks` once the
   !! double-diffusion path lands) diffusivity fields at layer
   !! interfaces; the existing `rdb_ocean_vdiff` Thomas solver
   !! consumes them directly via the `kv_source` argument.
   !!
   !! Two closure paths planned:
   !!   * PP81 (Pacanowski-Philander 1981) — Richardson-number
   !!     stability function; cheap, no boundary-layer detection,
   !!     reasonable for stratified-interior dynamics.  Live as of
   !!     this branch.
   !!   * KPP (Large, McWilliams & Doney 1994) — adds a Ri-bulk
   !!     boundary-layer detector + shape functions + non-local
   !!     transport term.  Future work; the BL state slots
   !!     (`bl_depth`, `gamma_t`, `gamma_s`) live here for that.
   !!
   !! Interface convention: `kv(:, :, k)` lives at the bottom interface
   !! of layer k (between layers k-1 and k); `kv(:, :, 1)` is the bed
   !! (forced to zero, closed BC), `kv(:, :, nz+1)` is the surface
   !! (forced to zero, closed top).  Same convention the
   !! `rdb_ocean_vdiff` impl assumes.
   !!
   !! **Convective adjustment** (`vmix_apply_convection`) is a Brunt-
   !! Vaisala-triggered CONTRIBUTOR: where the interior N^2 < n2_thresh
   !! (dense-over-light), it raises `kt -> max(kt, kd_conv)` and
   !! `kv -> max(kv, prandtl_conv*kd_conv)`, applied ONLY below the
   !! active surface boundary-layer depth (KPP `bl_depth` / EPBL `mld`
   !! own the BL response).  `kd_conv` defaults to 1 m^2/s -- five orders
   !! of magnitude above PP81's `Ri<0` clip ceiling (~1.01e-2 m^2/s) --
   !! and is admissible ONLY because `vdiff_apply_tracers` /
   !! `vdiff_apply_momentum` are backward-Euler (Thomas) tridiagonal
   !! solves: unconditionally stable, so `kappa*dt/dz^2 >> 1` does not
   !! blow up.  An explicit vertical diffusion could not use this
   !! coefficient.  Reference: Brunt-Vaisala convective trigger as
   !! implemented in CVMix (Griffies et al., `CVMix_convection`), wired
   !! in MOM6 as `USE_CVMix_CONVECTION` / `KD_CONV` / `PRANDTL_CONV` /
   !! `BV_SQR_CONV`; the underlying "represent convection as a very
   !! large diapycnal diffusivity" idea traces to Cox (1984) and
   !! Marotzke (1991).  Documented divergences from MOM6's
   !! `MOM_CVMix_conv`:
   !!   * D1 -- the N^2 trigger differences `ms%rho_layer` (a potential
   !!     density at the single global `eos%p_ref`, `rdb_eos.F90`
   !!     `eos_wright_impl`), not a locally-referenced (interface-
   !!     pressure) density as MOM6 evaluates.  Neglects thermobaricity;
   !!     identical to what PP81 already assumes for the same
   !!     expression.
   !!   * D2 -- uses `max(kt, kd_conv)` / `max(kv, prandtl_conv*kd_conv)`
   !!     rather than MOM6's additive `Kd = Kd + kd_col`.  The two differ
   !!     by at most the resolved interior Kd (<=1.01e-2, <=1% of
   !!     `kd_conv`) -- physically immaterial -- but `max()` is an exact,
   !!     idempotent floor, which matters because this runs every RK2
   !!     stage.  SAFE ONLY because `vmix_compute_pp81` ASSIGNS (not
   !!     accumulates) `kv`/`kt` every stage when `interior_closure ==
   !!     VMIX_INTERIOR_PP81` (the only reachable tag today): each stage
   !!     starts clean so `max()` cannot ratchet.  A future interior
   !!     closure that does not rewrite kv/kt every stage would let this
   !!     `max()` ratchet monotonically -- document this dependency for
   !!     whoever wires `VMIX_INTERIOR_LARGE94` / `VMIX_INTERIOR_CVMIX`.
   !!   * D3 -- MOM6 calls `calculate_CVMix_conv` BEFORE
   !!     `energetic_PBL_get_MLD` fills `BLD` for that step,
   !!     so CVMix_conv can mask against a
   !!     stale/zero BLD under EPBL (MOM6 emits a warning about this).
   !!     Roundabout does not need the warning: `vmix_apply_in_stage` runs
   !!     `epbl_compute` -> `epbl_merge_into_kv_kt` -> convection in the
   !!     SAME stage, so `epbl%mld` is current when convection reads it.
   !!     An improvement, not a gap.
   !!   * D4 -- writes `kv` and `kt` ONLY, never `ks`.  Unlike `kv`/`kt`,
   !!     `ks` is not rewritten by any per-stage contributor today (PP81
   !!     touches only `kv`/`kt`), so `ks = max(ks, kd_conv)` would
   !!     ratchet monotonically and never relax.  Salt convects for free
   !!     once the `ks` <- `kt` split (a separate PR) lands downstream of
   !!     this call, because that split copies the already-raised `kt`.
   !!     Until then `ks` stays the constant `pp81_kappa_bg` background
   !!     and salt genuinely does not convect -- a pre-existing
   !!     limitation this module does not worsen.
   !!
   !! `kd_max` interaction: `vmix_assemble` is the single downstream
   !! ceiling gate (contract above).  With the default `kd_max = huge`
   !! the convective 1 m^2/s passes through untouched; a caller who sets
   !! `kd_max < kd_conv` will silently cap the convective value -- that
   !! is the gate doing its job, not a convection bug.
   use rdb_constants, only: wp, GRAVITY, PI, H_VANISHED, DEG2RAD
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t, SEAWATER_CP, &
                                     sw_transmission
   use rdb_eos, only: eos_t
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_mem_report, only: arr_bytes
   implicit none
   private

   public :: ocean_vmix_t
   public :: vmix_compute_pp81
   public :: vmix_apply_kpp_overlay
   public :: vmix_apply_nonlocal_tendencies
   public :: vmix_add_kv_ml_invz2
   public :: vmix_split_kd_heat_salt
   public :: vmix_assemble
   public :: vmix_bkgnd_fill_impl
   public :: vmix_assemble_clip_henyey_impl
   public :: henyey_lat_factor_impl
   public :: vmix_resolve_kd_min
   public :: bkgnd_henyey_conflicts_profile
   public :: vmix_apply_convection
   public :: vmix_interior_closure_is_implemented
   public :: parse_kpp_sw_method
   public :: kpp_sw_method_is_implemented
   public :: kpp_surface_buoyancy_flux

   ! KPP shortwave-in-boundary-layer methods (MOM6 KPP_SHORTWAVE_METHOD).
   ! The surface buoyancy flux `B_0` is charged only for the SW ABSORBED
   ! inside the boundary layer; the fraction that leaks below `h_b`
   ! cannot stabilise it.  `Q_bl = Q_heat - I0·T(depth)` (Large,
   ! McWilliams & Doney 1994, App. B), with `T` the shared two-band
   ! `sw_transmission`.
   integer, parameter, public :: KPP_SW_ALL = 1
      !! `all_sw`: charge B_0 with the full net heat flux (legacy).
   integer, parameter, public :: KPP_SW_MXL = 2
      !! `mxl_sw` (default): subtract the SW that leaks below `h_b`.
   integer, parameter, public :: KPP_SW_LV1 = 3
      !! `lv1_sw`: subtract the SW that leaks below the top model layer.

   ! Interior-mixing closure tags.  Only `PP81` is implemented;
   ! `LARGE94` and `CVMIX` are reserved tags for future work.  KPP
   ! ships as an overlay on top of PP81 (`vmix_apply_kpp_overlay_*`),
   ! not as a separate interior tag.
   integer, parameter, public :: VMIX_INTERIOR_PP81 = 1
      !! Pacanowski-Philander 1981 — implemented; the default.
   integer, parameter, public :: VMIX_INTERIOR_LARGE94 = 2
      !! Large et al. 1994 interior closure (future; not yet wired).
   integer, parameter, public :: VMIX_INTERIOR_CVMIX = 3
      !! CVMix-compatible (future).

   real(wp), parameter, public :: HENYEY_MIN_SINLAT = 1.0e-10_wp
      !! Floor on |sin(latitude)| used ONLY inside the `N0_2Omega /
      !! |sin(latitude)|` ratio of the Henyey factor
      !! (`henyey_lat_factor_impl`) to avoid the 1/0 singularity exactly at
      !! the equator — Henyey, Wright & Flatte (1986) JGR 91:8487; the
      !! constant-`N0` simplification implemented here is Harrison &
      !! Hallberg (2008) JPO 38:1894.  A fixed constant, not a namelist
      !! knob.  The OUTER multiplication in the factor uses the TRUE
      !! (unfloored) |sin(latitude)|, so the assembled factor still goes
      !! smoothly to zero at the equator rather than blowing up.  PUBLIC so
      !! the unit tests can DERIVE the poleward-clamp oracle from it
      !! instead of hard-coding a magic number that would silently rot if
      !! this value ever changed.

   real(wp), parameter, public :: HENYEY_KD_MIN_FRAC = 0.01_wp
      !! Fraction of the scalar background tracer diffusivity used as the
      !! DEFAULT minimum diffusivity under the Henyey latitude scaling —
      !! MOM6 `KD_MIN`, whose documented default is `0.01*KD`.  Applied by
      !! `vmix_resolve_kd_min` when `bkgnd_kd_min` is left at its negative
      !! "unset" sentinel.  PUBLIC so the unit tests derive the floor
      !! oracle from it rather than hard-coding 1e-7.

   type :: ocean_vmix_t
      logical :: is_init = .false.
         !! True between `init` and `destroy`.  Prefer this to
         !! `allocated(...)` — tracks GPU device attachment too.

      ! ---- Scheme selection ----
      logical :: use_closure = .false.
         !! Master switch.  When `.false.` the driver skips the
         !! compute step and `rdb_ocean_vdiff` falls back to its
         !! scalar `K_v_*` constants.  Default off so all existing
         !! tests stay on the scalar path.
      logical :: use_kpp = .false.
         !! Enable KPP surface-boundary-layer scheme (future).
      integer :: interior_closure = VMIX_INTERIOR_PP81
         !! Interior mixing closure tag.

      ! ---- PP81 constants ----
      ! Default values from Pacanowski & Philander (1981) JPO 11:1443.
      real(wp) :: pp81_nu0 = 1.0e-2_wp
         !! Numerator viscosity at Ri=0 (m^2/s).
      real(wp) :: pp81_nu_bg = 1.0e-4_wp
         !! Background interior viscosity (m^2/s).
      real(wp) :: pp81_kappa_bg = 1.0e-5_wp
         !! Background tracer diffusivity (m^2/s).
      real(wp) :: pp81_alpha = 5.0_wp
         !! Richardson-number multiplier in (1 + α*Ri).  Original
         !! paper uses 5; some implementations use 4 or 10.
      real(wp) :: rho0 = 1035.0_wp
         !! Boussinesq reference density (kg/m^3).  Used in the
         !! N² = -g/ρ_0 * dρ/dz expression, in the friction velocity
         !! u_* = √(|τ|/ρ_0), and in the kinematic surface fluxes
         !! q_T = Q_heat/(ρ_0·cp) / q_S = Q_salt/ρ_0 that set the KPP
         !! surface buoyancy flux B_0.
         !!
         !! ASSIGNED FROM CONFIG by `configure_ocean_reference_density`,
         !! which copies the single rho0 of record (`&ocean_ic_nml rho_0`
         !! -> `eos%rho0`).  The literal here is only the pre-configure
         !! type default.  UNLIKE the other reference densities this one
         !! IS read on-device (`this%rho0` inside the `do concurrent`
         !! bodies of `vmix_compute_pp81` / `vmix_kpp_overlay_impl` /
         !! `vmix_convective_impl`), so under `mem:separate` it is only
         !! correct because the configure pass runs strictly BEFORE
         !! `ocean_state_enter_data`'s `copyin` — same as the `pp81_*`
         !! scalars beside it.  A later write needs `!$acc update device`.
      real(wp) :: shear2_floor = 1.0e-10_wp
         !! Lower bound on |∂u/∂z|² to avoid Ri = N²/0 blow-up in
         !! quiescent water columns.

      ! ---- EOS handle for the KPP buoyancy flux ----
      ! Shared flat-POD copy of `ocean_state%eos`, set once at
      ! configure.  The KPP convective-velocity scale B_0 needs the
      ! surface α (thermal expansion) and β (haline contraction);
      ! it reads `eos%alpha_T` / `eos%beta_S` so the closure tracks
      ! the dyn-core EOS coefficients (previously a private pair that
      ! was NEVER refreshed from the eos slot — the latent staleness
      ! bug this centralization fixes).  Maps onto the device with
      ! the parent `this`.  Both are DIMENSIONAL (kg/m^3 per degC /
      ! psu) — see `kpp_surface_buoyancy_flux` for the `1/rho_0` that
      ! turns them into a buoyancy flux.
      type(eos_t) :: eos

      ! ---- KPP-specific ----
      real(wp) :: ri_crit = 0.3_wp
         !! Critical bulk Richardson number for KPP BL depth.

      ! ---- KV_ML_INVZ2 surface-band viscosity (MOM6) ----
      real(wp) :: kv_ml_invz2 = 0.0_wp
         !! Extra kinematic viscosity (m²/s) inside the surface band
         !! of thickness `hmix_fixed`.  Profile is `(hmix_fixed / z)²`
         !! where z is the distance from the surface to the interface.
         !! Mirrors MOM6's KV_ML_INVZ2 + HMIX_FIXED combo — provides
         !! the dissipation that lets the wind drive thin surface
         !! layers without exciting a numerical eigenmode.  Zero
         !! (default) means no addition; the kv field comes purely
         !! from PP81 / KPP.
      real(wp) :: hmix_fixed = 20.0_wp
         !! Surface-band thickness (m) over which `kv_ml_invz2`
         !! profile is active.  MOM6 production default 20 m.
      real(wp) :: cs_nonlocal = 6.3_wp
         !! Non-local transport coefficient C_s (LMD94 eq 20).
         !! The full form `6.32·(1 - 0.5·exp(-σ/0.1))` is asymptotic;
         !! at σ > 0.1 it sits within 1% of the limit value 6.32, and
         !! the non-local flux is most important in that range.  Using
         !! the limit value keeps the kernel branch-free and matches
         !! MOM6's `KPP_Cstar` default.
      real(wp) :: c_vt2 = 1.8_wp
         !! Unresolved-turbulence coefficient for the V_t² term in the
         !! bulk-Ri denominator (LMD94 eq 23).  Folded form of
         !! `(C_v · √(-β_T) · √(c_s · ε)) / κ`.  Default 1.8 matches
         !! the standard LMD94 tuning.  Set to 0 to disable V_t² —
         !! recovers the shear-only bulk-Ri sweep bit-identically;
         !! useful as a discriminator in unit tests.
         !! `bl_depth(:, :)` from the previous step seeds w_s for V_t²,
         !! so on the very first call (h_b_lagged = 0) the w_*
         !! contribution evaluates to zero and the algorithm
         !! self-bootstraps.
      integer :: kpp_sw_method = KPP_SW_MXL
         !! Shortwave-in-boundary-layer method for `B_0` (MOM6
         !! `KPP_SHORTWAVE_METHOD`).  Only bites when penetrating SW is
         !! active (`sf%has_sw`); inert at `sw_pen_frac = 0` ⇒ the
         !! default `mxl_sw` is bit-identical to the legacy path there.

      ! ---- Assembly stage: backgrounds / ceilings / smoothing / guard ----
      ! The single downstream gate (`vmix_assemble`) that every interior
      ! and overlay closure feeds into before vdiff consumes kv/kt/ks.
      ! Defaults reproduce the pre-assembly numerics bit-for-bit: the
      ! backgrounds match the constant `pp81_*_bg` that PP81 already adds
      ! (so the floor is a no-op for the closure path), the ceilings are
      ! `huge` (no clip), smoothing is off, and the guard is off.
      real(wp) :: kv_bg = 1.0e-4_wp
         !! Background floor on momentum viscosity (m^2/s).  Default
         !! matches `pp81_nu_bg` so the floor never raises a PP81 value
         !! (PP81 writes `pp81_nu_bg + nu0·factor ≥ pp81_nu_bg`).
      real(wp) :: kt_bg = 1.0e-5_wp
         !! Background floor on temperature diffusivity (m^2/s).
         !! Default matches `pp81_kappa_bg`.
      real(wp) :: ks_bg = 1.0e-5_wp
         !! Background floor on salinity diffusivity (m^2/s).
         !! Default matches `pp81_kappa_bg`.
      real(wp) :: kv_max = huge(1.0_wp)
         !! Ceiling on momentum viscosity (m^2/s).  Default `huge` =
         !! no clip (bit-identical).  MOM6 `Kd_max` momentum analogue.
      real(wp) :: kd_max = huge(1.0_wp)
         !! Ceiling on tracer diffusivity (kt, ks) (m^2/s).  Default
         !! `huge` = no clip.  MOM6 `Kd_max`.
      integer :: kd_smooth_iterations = 0
         !! Number of 1-2-1 horizontal smoothing passes applied to
         !! kv/kt at each interior interface.  Default 0 = off (no
         !! smoothing kernel runs, bit-identical).  MOM6 `Kd_smooth`.
      logical :: vmix_guard = .false.
         !! Debug-gated negative/NaN guard.  When `.true.` the assembly
         !! scans kv/kt/ks for a negative or NaN value and `error stop`s
         !! (or returns a non-zero status via the testable path).
         !! Default off (cheap reduction skipped) so production runs pay
         !! nothing.

      ! ---- C7 background mixing: Bryan-Lewis XOR Henyey ----
      ! Two MUTUALLY EXCLUSIVE alternatives to the scalar background floor,
      ! matching the reference code's one-background-scheme rule:
      !   * `bkgnd_profile` — Bryan & Lewis (1979) depth profile, replacing
      !     the scalar floor with a per-interface `kd_bg` field.
      !   * `bkgnd_henyey`  — Henyey (1986) latitude factor on the SCALAR
      !     `kt_bg`/`ks_bg`, floored at `bkgnd_kd_min`.
      ! Enabling both fails loud at configure.  Both default off ⇒ the
      ! scalar path is used verbatim (bit-identical).
      !
      ! When `bkgnd_profile` is on, the SCALAR `kt_bg`/`ks_bg` additive
      ! floor in `vmix_assemble` is replaced by a per-interface
      ! Bryan & Lewis (1979) JGR 84:2503 depth profile
      !   Kd_bg(z) = Kd_sfc + (Kd_deep - Kd_sfc)*(0.5 + atan((|z|-z0)/Delta)/PI)
      ! evaluated from the CURRENT column interface depths (cumulative
      ! h_layer from the surface down) so it is correct under any vcoord.
      ! Momentum follows MOM6: Kv_bg = bkgnd_prandtl * Kd_bg.  Default off
      ! ⇒ the scalar background path is used verbatim (bit-identical).
      logical :: bkgnd_profile = .false.
         !! Master switch for the depth-varying Bryan-Lewis background.
         !! Default `.false.` ⇒ `vmix_assemble` uses the scalar
         !! kv_bg/kt_bg/ks_bg floor exactly as before (bit-identical).
      real(wp) :: bkgnd_kd_sfc = 1.0e-5_wp
         !! Surface-asymptote background tracer diffusivity (m^2/s).
         !! MOM6 BRYAN_LEWIS_C2-side default scale.
      real(wp) :: bkgnd_kd_deep = 1.3e-4_wp
         !! Deep-asymptote background tracer diffusivity (m^2/s).
      real(wp) :: bkgnd_z0 = 2500.0_wp
         !! Transition-centre depth (m, positive down) where the profile
         !! reaches the (sfc+deep)/2 midpoint.
      real(wp) :: bkgnd_delta = 222.0_wp
         !! Transition half-width (m): atan argument is (|z|-z0)/Delta.
      real(wp) :: bkgnd_prandtl = 1.0_wp
         !! Background Prandtl number — Kv_bg = bkgnd_prandtl * Kd_bg
         !! (MOM6 ties the background viscosity to the background tracer
         !! diffusivity through a Prandtl factor).  Default 1.0.
      logical :: bkgnd_henyey = .false.
         !! Henyey, Wright & Flatte (1986) JGR 91:8487 latitude-dependent
         !! internal-wave factor, scaling the SCALAR background tracer
         !! diffusivities `kt_bg`/`ks_bg` by a horizontal-only factor
         !! `L(phi)` computed from `geolatT` (see `henyey_lat_factor_impl`
         !! for the exact form), with the result floored at `bkgnd_kd_min`:
         !!
         !!   kt_floor(i,j) = max(bkgnd_kd_min, kt_bg * L(phi))
         !!
         !! The implemented variant is the SIMPLIFIED one of Harrison &
         !! Hallberg (2008) JPO 38:1894, which assumes the in-situ
         !! stratification equals a constant reference `N0` rather than the
         !! evolving column N — so the factor depends only on latitude +
         !! `bkgnd_henyey_n0_2omega` / `bkgnd_henyey_max_lat` and needs no
         !! per-step recompute.
         !!
         !! MUTUALLY EXCLUSIVE with `bkgnd_profile` (Bryan-Lewis), matching
         !! the reference formulation, which selects ONE background scheme
         !! and FATALs when a second is requested.  Enabling both fails loud
         !! at configure.  This is also the cheaper arrangement: no
         !! per-interface `kd_bg` field is filled at all, the latitude
         !! factor is a per-column scalar folded straight into the existing
         !! floor/ceiling clip.
         !!
         !! Also requires a non-cartesian `grid_config` (validated fail-loud
         !! at configure): `geolatT` is identically zero on a cartesian grid,
         !! so every column would take the equatorial factor `L(0 deg) = 0`
         !! and the background would collapse to a uniform `bkgnd_kd_min`
         !! everywhere — a latitude parameterisation on a grid with no
         !! meaningful latitude.  Default `.false.` ⇒ the scalar background
         !! floor path is used verbatim (bit-identical).
         !!
         !! SCOPE: the factor scales the two TRACER background floors only;
         !! the momentum floor `kv_bg` is left alone.  Roundabout's scalar
         !! background path deliberately carries `kv_bg` as an INDEPENDENT
         !! momentum floor rather than `prandtl * kt_bg` (the shipped
         !! defaults 1e-4 / 1e-5 imply Pr = 10), so there is no single
         !! background `Kd` for the reference code's `Kv_bkgnd =
         !! PRANDTL_BKGND * Kd` tie to reproduce here.  Henyey scales the
         !! diapycnal DIFFUSIVITY, which is what `kt_bg`/`ks_bg` are.
      real(wp) :: bkgnd_kd_min = -1.0_wp
         !! Minimum background tracer diffusivity (m^2/s) under the Henyey
         !! latitude scaling — MOM6 `KD_MIN`, applied as
         !! `max(Kd_min, Kd * L(phi))`.  Without it `kd_bg` would collapse
         !! toward zero at the equator (`L(0 deg) = 0` exactly) and at the
         !! poleward `max_lat` clamp, instead of the documented behaviour
         !! "the Henyey profile is returned to the MINIMUM diffusivity".
         !!
         !! NEGATIVE = unset sentinel ⇒ resolved to
         !! `HENYEY_KD_MIN_FRAC * kt_bg` (MOM6's `0.01*KD` default) by
         !! `vmix_resolve_kd_min`, which `vmix_assemble` calls on every
         !! entry so a directly-constructed slot (unit tests) gets the same
         !! default as the configure path.  Read only when `bkgnd_henyey`
         !! is on.
      real(wp) :: bkgnd_henyey_n0_2omega = 20.0_wp
         !! Ratio of the assumed reference buoyancy frequency `N0` to twice
         !! the planetary rotation rate (nondim).  Physically
         !! `N0 >> 2*Omega` always, so this stays well above 1
         !! (configure-time `min=1` guard keeps the internal `acosh`
         !! argument in-domain even under a misconfigured value).
      real(wp) :: bkgnd_henyey_max_lat = 95.0_wp
         !! Latitude (degN) poleward of which the factor is reset to its
         !! equatorial (near-zero) floor.  Deliberately > 90 by default so
         !! the clamp is INERT for any real latitude out of the box; lower
         !! it to activate the optional poleward cutoff.

      ! ---- Convective adjustment (Brunt-Vaisala trigger, CVMix_conv-style) ----
      ! `&ocean_conv_nml`.  A CONTRIBUTOR that raises kv/kt via max() where
      ! the interior N^2 < n2_thresh, applied strictly below the active
      ! surface boundary layer.  Default off => bit-identical.  See the
      ! module docstring for the D1-D4 divergences from MOM6's
      ! MOM_CVMix_conv and the kd_max interaction.
      logical :: conv_enable = .false.
         !! Master switch.  Requires `use_closure` + thermodynamics
         !! (validated at configure).  Default off => bit-identical.
      real(wp) :: conv_kd = 1.0_wp
         !! Convective tracer diffusivity (m^2/s).  MOM6 `KD_CONV`
         !! default 1.0 -- ~1e5x the background, admissible only because
         !! vdiff is an unconditionally-stable backward-Euler solve.
      real(wp) :: conv_prandtl = 1.0_wp
         !! Kv_conv = conv_prandtl * kd_conv.  MOM6 `PRANDTL_CONV`
         !! default 1.0.
      real(wp) :: conv_n2_thresh = 0.0_wp
         !! Trigger threshold on N^2 (s^-2).  MOM6 `BV_SQR_CONV` default
         !! 0.0.  Strict `<` so an exactly-neutral interface (N^2 == 0)
         !! does not trigger.

      ! ---- Double diffusion (salt fingering + diffusive convection) ----
      ! `&ocean_ddiff_nml`.  NOT a kv/kt contributor -- it is folded INTO
      ! `vmix_split_kd_heat_salt` (a pre-split ks write would be clobbered
      ! by `ks := kt`), producing an ASYMMETRIC ks-vs-kt divergence:
      !   ks = kt_pre + kd_extra_s ;  kt = kt_pre + kd_extra_t
      ! Interior interfaces only, branched on the SIGNED alpha*dT / beta*dS
      ! (never a pre-divided R_rho -- dodges R_rho<=0, ->inf, 0/0).  The
      ! CVMix (Large et al. 1994 / Marmorino-Caldwell 1976 / Kelley 1990)
      ! closed forms; constants are the CVMix defaults.  Default off =>
      ! bit-identical.  v1 uses the constant linear-EOS alpha_T/beta_S off
      ! `this%eos` (bitwise-exact for the linear EOS); nonlinear
      ! per-interface derivatives are a deferred refinement.
      logical :: ddiff_enable = .false.
         !! Master switch (MOM6 `USE_CVMIX_DDIFF`).  Requires
         !! thermodynamics (validated at configure).  Default off.
      real(wp) :: ddiff_strat_param_max = 2.55_wp
         !! R_rho salt-fingering cutoff (CVMix `STRAT_PARAM_MAX`).  Above
         !! it fingering diffusivity is zero.
      real(wp) :: ddiff_kappa_s = 1.0e-4_wp
         !! Leading salt-fingering salinity diffusivity K_f (m^2/s,
         !! CVMix `KAPPA_DDIFF_S`).  K_T = 0.7*K_S (0.7 hard-wired).
      real(wp) :: ddiff_exp1 = 1.0_wp
         !! Inner (bracket) exponent of the fingering clamped form
         !! (CVMix `DDIFF_EXP1`).
      real(wp) :: ddiff_exp2 = 3.0_wp
         !! Outer exponent of the fingering clamped form (CVMix
         !! `DDIFF_EXP2`); `exp1=1, exp2=3` is the Large et al. cubic.
      real(wp) :: ddiff_param1 = 0.909_wp
         !! MC76 diffusive-convection exterior coeff (CVMix
         !! `KAPPA_DDIFF_PARAM1`).
      real(wp) :: ddiff_param2 = 4.6_wp
         !! MC76 middle coeff (CVMix `KAPPA_DDIFF_PARAM2`).
      real(wp) :: ddiff_param3 = -0.54_wp
         !! MC76 interior coeff (CVMix `KAPPA_DDIFF_PARAM3`).
      real(wp) :: ddiff_mol_diff = 1.5e-6_wp
         !! Molecular diffusivity scaling the convection branch (m^2/s,
         !! CVMix `MOL_DIFF`) -- the *molecular* value, NOT a background
         !! eddy diffusivity.
      logical :: ddiff_use_k90 = .false.
         !! Diffusive-convection form: `.false.` = Marmorino-Caldwell 1976
         !! (MC76, default), `.true.` = Kelley 1990 (K90).

      ! ---- Diagnostic / prognostic 2D fields (KPP) ----
      real(wp), allocatable :: bl_depth(:, :)
         !! KPP boundary-layer depth (m, positive down).
      real(wp), allocatable :: b0(:, :)
         !! Surface buoyancy flux `B_0` (m^2/s^3) the KPP overlay's
         !! convective scale was built from, persisted per column on the
         !! SECOND pass (the one that uses the freshly-diagnosed
         !! `bl_depth`).  Sign convention matches EPBL's `epbl%b0`:
         !! `> 0` stabilizing (heating / freshening), `< 0` destabilizing
         !! (cooling / salting), and `w_*^3 = max(0, -b0)*bl_depth`.
         !! Diagnostic only — nothing reads it back into the closure; it
         !! exists so the two boundary-layer schemes' surface forcing can
         !! be compared directly (`test_ocean_buoyancy_flux`).  Zero until
         !! the first KPP overlay call.

      ! ---- 3D mixing coefficients on layer interfaces ----
      ! Shape (nx, ny, nz_ml+1); k=1 at the bed, k=nz_ml+1 at the
      ! free surface.  Consumed by the vertical-diffusion solve.
      real(wp), allocatable :: kv(:, :, :)
         !! Momentum (u, v) eddy viscosity at interfaces.
      real(wp), allocatable :: kt(:, :, :)
         !! Temperature eddy diffusivity.
      real(wp), allocatable :: ks(:, :, :)
         !! Salinity eddy diffusivity.

      ! ---- KPP non-local (counter-gradient) transport ----
      ! Interface-located downward flux in tracer-units · m/s.  Only
      ! non-zero inside the BL under destabilizing surface forcing
      ! (B_0 < 0).  Divergence ∂γ/∂z appears as a source in the
      ! tracer-vdiff RHS.  Allocated alongside `kv` / `kt` when
      ! `nz_ml` is supplied; same shape `(nx, ny, nz_ml + 1)`.
      real(wp), allocatable :: gamma_t(:, :, :)
         !! Non-local temperature flux (°C·m/s) at interfaces.
      real(wp), allocatable :: gamma_s(:, :, :)
         !! Non-local salinity flux (PSU·m/s) at interfaces.

      ! ---- Assembly smoothing scratch ----
      ! Double-buffer for the 1-2-1 horizontal smoothing kernel — a
      ! 1-2-1 pass cannot run in place (neighbours must read pre-pass
      ! values).  Allocated in `init` and mapped in `enter_data`
      ! alongside kv/kt so the smoothing kernel never lazy-attaches a
      ! component to an already-mapped parent (that collides with the
      ! present table).  Only read when `kd_smooth_iterations > 0`.
      real(wp), allocatable :: smooth_scratch(:, :, :)
         !! Scratch (nx, ny, nz+1) for one smoothing pass.

      ! ---- C7 Bryan-Lewis background floor field ----
      ! Per-interface tracer-background diffusivity (m^2/s), shape
      ! (nx, ny, nz+1).  Filled by `vmix_bkgnd_fill` (called from
      ! `vmix_assemble`) each stage from the current column thicknesses
      ! when `bkgnd_profile` is on; otherwise never read.  Allocated in
      ! `init` so the enter_data orchestrator maps it unconditionally
      ! (avoids a lazy device attach on an already-mapped parent).
      real(wp), allocatable :: kd_bg(:, :, :)
         !! Bryan-Lewis depth-profile tracer background (m^2/s) at
         !! interfaces.  kv background floor = bkgnd_prandtl * kd_bg.
   contains
      procedure, non_overridable :: init => ocean_vmix_init
      procedure, non_overridable :: destroy => ocean_vmix_destroy
      procedure, non_overridable :: enter_data => ocean_vmix_enter_data
      procedure, non_overridable :: exit_data => ocean_vmix_exit_data
      procedure, non_overridable :: bytes => ocean_vmix_bytes
      procedure, non_overridable :: seed_backgrounds => vmix_seed_backgrounds
   end type ocean_vmix_t

contains

   pure function vmix_interior_closure_is_implemented(code) result(ok)
      !! `.true.` only for `VMIX_INTERIOR_PP81` — the one interior
      !! closure with a real kernel today.  `VMIX_INTERIOR_LARGE94` and
      !! `VMIX_INTERIOR_CVMIX` are declared-but-unimplemented reservations
      !! (no kernel); selecting one would leave `kv`/`kt` stale (no
      !! interior mixing at all).  `interior_closure` has no namelist key
      !! yet, so this is defence-in-depth (PR-6) for the next code/config
      !! consumer that sets it — the predicate is the gate.
      integer, intent(in) :: code
      logical :: ok
      ok = (code == VMIX_INTERIOR_PP81)
   end function vmix_interior_closure_is_implemented

   subroutine ocean_vmix_init(this, grid, nz_ml)
      !! Allocate the kv / kt / ks diffusivity fields at layer
      !! interfaces.  Default values: kv = kv_bg (= pp81_nu_bg),
      !! kt = ks = kt_bg (= pp81_kappa_bg) at interior interfaces;
      !! boundary interfaces k=1 and k=nz+1 are zeroed (closed BC).
      !! BL fields (`bl_depth`, `gamma_*`) seed at zero.  `ks` is
      !! allocated so the assembly stage (`vmix_assemble`) can
      !! floor/clip it; `vmix_split_kd_heat_salt` derives its live
      !! value from `kt` every stage, and `vdiff_apply_tracers`
      !! consumes it for salinity + every passive tracer.
      !!
      !! `smooth_scratch` is NOT allocated here — it is allocated lazily
      !! in `configure_ocean_vmix` only when `kd_smooth_iterations > 0`
      !! (before `ocean_state_enter_data`).  The `enter_data` path
      !! already guards on `allocated(smooth_scratch)`.
      class(ocean_vmix_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      integer, intent(in), optional :: nz_ml
      integer :: nx, ny, nz

      nx = grid%nx_total
      ny = grid%ny_total
      nz = 1
      if (present(nz_ml)) nz = nz_ml

      allocate (this%kv(nx, ny, nz + 1))
      allocate (this%kt(nx, ny, nz + 1))
      allocate (this%ks(nx, ny, nz + 1))
      allocate (this%kd_bg(nx, ny, nz + 1))
      ! kv_bg/kt_bg/ks_bg + kv/kt/ks/kd_bg all seed from pp81_* here
      ! (structural invariant, see `vmix_seed_backgrounds`).  The SAME
      ! call re-derives them after `&ocean_vmix_nml pp81_*` reaches the
      ! slot at configure — see `rdb_ocean_setup.F90:configure_ocean_lateral`.
      call this%seed_backgrounds()
      allocate (this%bl_depth(nx, ny), source=0.0_wp)
      allocate (this%b0(nx, ny), source=0.0_wp)
      allocate (this%gamma_t(nx, ny, nz + 1), source=0.0_wp)
      allocate (this%gamma_s(nx, ny, nz + 1), source=0.0_wp)

      this%is_init = .true.
   end subroutine ocean_vmix_init

   pure subroutine vmix_seed_backgrounds(this)
      !! Seed `kv_bg`/`kt_bg`/`ks_bg` and the `kv`/`kt`/`ks`/`kd_bg` arrays
      !! from the current `pp81_nu_bg`/`pp81_kappa_bg` fields, then zero
      !! the closed-BC boundary interfaces on `kv`/`ks`.  Extracted out of
      !! `ocean_vmix_init` (§2E structural invariant: `kv_bg == pp81_nu_bg`,
      !! `kt_bg == ks_bg == pp81_kappa_bg`, so `vmix_assemble`'s background
      !! floor is a no-op for the shipped PP81 closure) so BOTH the initial
      !! seed and the `&ocean_vmix_nml pp81_*` config-copy re-derive
      !! consistently — a bare field copy without this call would leave
      !! the assembly floor clamped against the OLD (type-default)
      !! background even after a user sets a new one.  Requires
      !! `kv`/`kt`/`ks`/`kd_bg` already allocated (true after `init`; the
      !! config-copy call runs strictly after `init_from_config`).
      class(ocean_vmix_t), intent(inout) :: this
      integer :: nz1

      this%kv_bg = this%pp81_nu_bg
      this%kt_bg = this%pp81_kappa_bg
      this%ks_bg = this%pp81_kappa_bg

      this%kv = this%pp81_nu_bg
      this%kt = this%pp81_kappa_bg
      this%ks = this%pp81_kappa_bg
      this%kd_bg = this%pp81_kappa_bg

      ! Zero boundary interfaces (closed BC contract).
      nz1 = size(this%kv, 3)
      this%kv(:, :, 1) = 0.0_wp
      this%kv(:, :, nz1) = 0.0_wp
      this%ks(:, :, 1) = 0.0_wp
      this%ks(:, :, nz1) = 0.0_wp
   end subroutine vmix_seed_backgrounds

   subroutine ocean_vmix_destroy(this)
      class(ocean_vmix_t), intent(inout) :: this
      this%is_init = .false.
      if (allocated(this%bl_depth)) deallocate (this%bl_depth)
      if (allocated(this%b0)) deallocate (this%b0)
      if (allocated(this%kv)) deallocate (this%kv)
      if (allocated(this%kt)) deallocate (this%kt)
      if (allocated(this%ks)) deallocate (this%ks)
      if (allocated(this%gamma_t)) deallocate (this%gamma_t)
      if (allocated(this%gamma_s)) deallocate (this%gamma_s)
      if (allocated(this%smooth_scratch)) deallocate (this%smooth_scratch)
      if (allocated(this%kd_bg)) deallocate (this%kd_bg)
   end subroutine ocean_vmix_destroy

   subroutine ocean_vmix_enter_data(this)
      class(ocean_vmix_t), intent(inout) :: this
      select type (this)
      type is (ocean_vmix_t)
         call ocean_vmix_enter_data_impl(this)
      end select
   end subroutine ocean_vmix_enter_data

   subroutine ocean_vmix_enter_data_impl(this)
      type(ocean_vmix_t), intent(inout) :: this
      if (allocated(this%kv)) then
         !$acc enter data copyin(this%kv)
      end if
      if (allocated(this%kt)) then
         !$acc enter data copyin(this%kt)
      end if
      if (allocated(this%ks)) then
         !$acc enter data copyin(this%ks)
      end if
      if (allocated(this%smooth_scratch)) then
         !$acc enter data copyin(this%smooth_scratch)
      end if
      if (allocated(this%bl_depth)) then
         !$acc enter data copyin(this%bl_depth)
      end if
      if (allocated(this%b0)) then
         !$acc enter data copyin(this%b0)
      end if
      if (allocated(this%gamma_t)) then
         !$acc enter data copyin(this%gamma_t)
      end if
      if (allocated(this%gamma_s)) then
         !$acc enter data copyin(this%gamma_s)
      end if
      if (allocated(this%kd_bg)) then
         !$acc enter data copyin(this%kd_bg)
      end if
   end subroutine ocean_vmix_enter_data_impl

   subroutine ocean_vmix_exit_data(this)
      class(ocean_vmix_t), intent(inout) :: this
      select type (this)
      type is (ocean_vmix_t)
         call ocean_vmix_exit_data_impl(this)
      end select
   end subroutine ocean_vmix_exit_data

   subroutine ocean_vmix_exit_data_impl(this)
      type(ocean_vmix_t), intent(inout) :: this
      if (allocated(this%kv)) then
         !$acc exit data delete(this%kv)
      end if
      if (allocated(this%kt)) then
         !$acc exit data delete(this%kt)
      end if
      if (allocated(this%ks)) then
         !$acc exit data delete(this%ks)
      end if
      if (allocated(this%smooth_scratch)) then
         !$acc exit data delete(this%smooth_scratch)
      end if
      if (allocated(this%bl_depth)) then
         !$acc exit data delete(this%bl_depth)
      end if
      if (allocated(this%b0)) then
         !$acc exit data delete(this%b0)
      end if
      if (allocated(this%gamma_t)) then
         !$acc exit data delete(this%gamma_t)
      end if
      if (allocated(this%gamma_s)) then
         !$acc exit data delete(this%gamma_s)
      end if
      if (allocated(this%kd_bg)) then
         !$acc exit data delete(this%kd_bg)
      end if
   end subroutine ocean_vmix_exit_data_impl

   pure subroutine vmix_compute_pp81(grid, this, ms)
      !! Pacanowski-Philander (1981) Richardson-number closure.
      !! Inlined kernel — keeps the `do concurrent` body adjacent to
      !! its derived-type accesses.  We tried the outer-shim +
      !! `_impl` pattern but NVHPC's stdpar codegen produced more
      !! descriptor-marshalling memcpys at the shim boundary than the
      !! direct-access form generates inside the kernel.  Direct
      !! `ms%foo(i,j,k)` access is what other working hot-path
      !! kernels (continuity, coriolis_adv, the barotropic substep) use; the
      !! shim pattern was an experiment that didn't help here.
      type(hgrid_t), intent(in) :: grid
      type(ocean_vmix_t), intent(inout) :: this
      type(multilayer_state_t), intent(in) :: ms

      integer :: i, j, k, nx, ny, nz
      real(wp) :: u_km1, v_km1, u_k, v_k
      real(wp) :: du_dz, dv_dz, shear2, dz_face
      real(wp) :: n2, ri, denom, ri_factor

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml

      ! Boundary interfaces: keep them at zero (closed top + bottom).
      do concurrent(j=1:ny, i=1:nx)
         this%kv(i, j, 1) = 0.0_wp
         this%kt(i, j, 1) = 0.0_wp
         this%kv(i, j, nz + 1) = 0.0_wp
         this%kt(i, j, nz + 1) = 0.0_wp
      end do

      ! Interior interfaces k = 2..nz.  Richardson stability factor
      ! in one pass per (i, j, k).  Loop spans the full (i, j) extent
      ! including walls — face arrays are (nx+1, ny, *) and
      ! (nx, ny+1, *), so i+1 / j+1 reach valid storage.  At closed
      ! walls the face velocities are zero ⇒ shear² floors ⇒ Ri huge
      ! ⇒ same kv as adjacent interior.  A wall-only fallback that
      ! sets kv = pp81_nu_bg breaks horizontal symmetry by O(nu0/denom²)
      ! and seeds a baroclinic-noise instability in stratified
      ! closed basins.
      do concurrent(k=2:nz, j=1:ny, i=1:nx) &
         local(u_km1, v_km1, u_k, v_k, du_dz, dv_dz, &
               shear2, dz_face, n2, ri, denom, ri_factor)
         u_km1 = 0.5_wp*(ms%u_face_x_layer(i, j, k - 1) + ms%u_face_x_layer(i + 1, j, k - 1))
         v_km1 = 0.5_wp*(ms%v_face_y_layer(i, j, k - 1) + ms%v_face_y_layer(i, j + 1, k - 1))
         u_k = 0.5_wp*(ms%u_face_x_layer(i, j, k) + ms%u_face_x_layer(i + 1, j, k))
         v_k = 0.5_wp*(ms%v_face_y_layer(i, j, k) + ms%v_face_y_layer(i, j + 1, k))

         dz_face = 0.5_wp*(ms%h_layer(i, j, k - 1) + ms%h_layer(i, j, k))
         if (dz_face <= 0.0_wp) then
            this%kv(i, j, k) = this%pp81_nu_bg
            this%kt(i, j, k) = this%pp81_kappa_bg
            cycle
         end if

         du_dz = (u_k - u_km1)/dz_face
         dv_dz = (v_k - v_km1)/dz_face
         shear2 = max(du_dz*du_dz + dv_dz*dv_dz, this%shear2_floor)

         n2 = -GRAVITY*(ms%rho_layer(i, j, k) - ms%rho_layer(i, j, k - 1))/ &
              (this%rho0*dz_face)
         ri = n2/shear2
         denom = 1.0_wp + this%pp81_alpha*max(ri, 0.0_wp)
         ri_factor = 1.0_wp/(denom*denom)
         this%kv(i, j, k) = this%pp81_nu_bg + this%pp81_nu0*ri_factor
         this%kt(i, j, k) = this%pp81_kappa_bg + this%pp81_nu0*ri_factor/denom
      end do
   end subroutine vmix_compute_pp81

   pure subroutine vmix_apply_kpp_overlay(grid, this, ms, ss, sf)
      !! Thin host-side shim over `vmix_kpp_overlay_impl` — same
      !! signature as before this PR, so the call site in
      !! `rdb_ocean_dyn.F90` is unchanged.  Selects the shortwave
      !! irradiance source HOST-SIDE (`sf%sw_from_qsw`): the PR-12 `q_sw`
      !! component is allocated only under `use_components`, so it is only
      !! ever passed on the branch guarded by the host flag
      !! (validate_config forces `enable_components` when
      !! `sw_source="q_sw"`, making this total).  `sf%has_sw` is passed as
      !! the `sw_active` gate — false ⇒ the `_impl`'s `B_0` reduces to the
      !! unmodified legacy source line, bit-for-bit.
      type(hgrid_t), intent(in) :: grid
      type(ocean_vmix_t), intent(inout) :: this
      type(multilayer_state_t), intent(in) :: ms
      type(ocean_surface_stress_t), intent(in) :: ss
      type(ocean_surface_flux_t), intent(in) :: sf
      if (sf%sw_from_qsw) then
         call vmix_kpp_overlay_impl(grid, this, ms, ss, sf, &
                                    grid%nx_total, grid%ny_total, sf%q_sw, &
                                    sf%has_sw, sf%sw_pen_frac, sf%sw_band_ratio, &
                                    sf%sw_zeta1, sf%sw_zeta2, this%kpp_sw_method)
      else
         call vmix_kpp_overlay_impl(grid, this, ms, ss, sf, &
                                    grid%nx_total, grid%ny_total, sf%Q_heat, &
                                    sf%has_sw, sf%sw_pen_frac, sf%sw_band_ratio, &
                                    sf%sw_zeta1, sf%sw_zeta2, this%kpp_sw_method)
      end if
   end subroutine vmix_apply_kpp_overlay

   pure subroutine vmix_kpp_overlay_impl(grid, this, ms, ss, sf, &
                                         nx_arg, ny_arg, sw_src, sw_active, &
                                         sw_pen_frac, sw_R, sw_zeta1, sw_zeta2, &
                                         sw_method)
      !! KPP boundary-layer overlay on top of the interior closure
      !! already in `this%kv` / `this%kt`.  Phase 1 was shear-driven
      !! only; Phase 2 added the convective velocity scale `w_*` and
      !! γ_T/γ_S non-local transport; Phase 3 (this revision) adds
      !! the V_t² unresolved-turbulence term in the bulk-Ri
      !! denominator (LMD94 eq 23).  Surface BL is now feature-
      !! complete except for Langmuir / Stokes enhancement.
      !!
      !! Two passes per column:
      !!
      !!   1. Bulk-Ri sweep top-down (from the top layer centre) to
      !!      find the boundary-layer depth h_b.  Reference values
      !!      are u_ref, v_ref, B_ref at the top layer centre.  At
      !!      every subsequent layer centre k:
      !!         Ri_bulk(k) = (B_ref - B(k)) * (z_k - z_ref) /
      !!                      max(|V(k) - V_ref|² + V_t²(d_k),
      !!                          shear²_floor)
      !!      with B(k) = -g·ρ(k)/ρ_0.  When Ri_bulk first crosses
      !!      `ri_crit` (default 0.3), linearly interpolate between
      !!      the previous and current layer centres to recover the
      !!      crossing depth.  If we walk to the bed without
      !!      crossing, h_b = total column depth.
      !!
      !!      V_t²(d) = (c_vt2 / ri_crit) · d · N(d) · w_s_col
      !!      adds an unresolved-turbulence contribution that
      !!      sharpens the BL depth diagnosis when grid-scale shear
      !!      is weak (LMD94 eq 23).  N(d) = √(max(0, ΔB/Δd)) is the
      !!      local buoyancy frequency.  w_s_col is computed once
      !!      per column using `this%bl_depth(i, j)` from the
      !!      previous step (lagged h_b) → self-bootstrapping on
      !!      the first call.  Setting `c_vt2 = 0` disables V_t²
      !!      bit-identically.
      !!
      !!   2. KPP overlay: at every interface k with depth from the
      !!      surface d(k) < h_b, compute
      !!         σ = d(k) / h_b
      !!         G(σ) = σ · (1 - σ)²
      !!         w_s = √(u_*² + w_*²)
      !!         kv_kpp = h_b · w_s · G(σ)
      !!      and `max` it into `this%kv(:, :, k)` and
      !!      `this%kt(:, :, k)`.  Outside the BL the interior
      !!      values (PP81 etc.) carry through unchanged.
      !!
      !!   u_* = √(|τ|/ρ_0) is the shear-driven scale.
      !!   w_* = max(0, -B_0·h_b)^(1/3) is the convective scale —
      !!         non-zero only when the surface buoyancy flux is
      !!         destabilizing (B_0 < 0).
      !!   B_0 = (g/ρ_0)·(α_T·F_T - β_S·F_S) where F_T = Q_heat/
      !!         (ρ_0·cp) and F_S = Q_salt/ρ_0 are the kinematic
      !!         surface heat / salt fluxes.  B_0 > 0 = stabilizing
      !!         (heating / freshening), B_0 < 0 = destabilizing
      !!         (cooling / salting).  The `1/ρ_0` is load-bearing:
      !!         `α_T` / `β_S` here are the DIMENSIONAL linear-EOS
      !!         sensitivities (kg m^-3 per degC / psu), not the
      !!         fractional ones — see `kpp_surface_buoyancy_flux`,
      !!         which owns the expression and the convention.  B_0 is
      !!         persisted per column into `this%b0` on the second pass
      !!         (diagnostic; the same quantity as `epbl%b0`).
      !!
      !! `sf` is REQUIRED (A7): B_0 is computed PER COLUMN from the 2D
      !! `sf%Q_heat / Q_salt` fields inside the kernel, so file-driven
      !! spatially-varying forcing (data-override, A3) feeds KPP with no
      !! further change.  Zero-flux behaviour = zero-filled fields.
      !!
      !! Penetrating shortwave (PR-21): when `sw_active`, `B_0` is charged
      !! only for the SW ABSORBED INSIDE the boundary layer.  The heat the
      !! BL actually feels is `Q_bl = Q_heat - I0·T(depth)` (Large,
      !! McWilliams & Doney 1994, App. B), where `I0 = sw_pen_frac·sw_src`
      !! and `T` is the shared two-band `sw_transmission` (Paulson &
      !! Simpson 1977).  `sw_method` selects the reference depth:
      !! `all_sw` ⇒ no correction (legacy), `mxl_sw` ⇒ `T(h_b)`,
      !! `lv1_sw` ⇒ `T(h_layer(nz))`.  `Q_heat`/`Q_salt` stay on `sf`;
      !! the only new array dummy is the explicit-shape `sw_src`.
      !! NOTE (Roundabout KPP fidelity, roadmap trap #4 / PR-11): the overlay
      !! has no Monin-Obukhov length, so `B_0 > 0` ⇒ `w_* ≡ 0` and the SW
      !! method is a numerical no-op under net heating; it bites only in
      !! the sunny-but-net-cooling regime (`B_0 <= 0`, `q_sw > 0`).
      type(hgrid_t), intent(in) :: grid
      type(ocean_vmix_t), intent(inout) :: this
      type(multilayer_state_t), intent(in) :: ms
      type(ocean_surface_stress_t), intent(in) :: ss
      type(ocean_surface_flux_t), intent(in) :: sf
      integer, intent(in) :: nx_arg, ny_arg
         !! Grid extents — declared before `sw_src` (decl-order hook).
      real(wp), intent(in) :: sw_src(nx_arg, ny_arg)
         !! Caller-selected irradiance source (`sf%Q_heat` or `sf%q_sw`),
         !! explicit-shape (per-RK2-stage kernel: no assumed-shape waiver).
      logical, intent(in) :: sw_active
         !! Host-side `sf%has_sw` gate — false ⇒ `B_0` uses the unmodified
         !! legacy source line (bit-identity).
      real(wp), intent(in) :: sw_pen_frac, sw_R, sw_zeta1, sw_zeta2
         !! Two-band SW parameters (from `sf`), by value.
      integer, intent(in) :: sw_method
         !! `KPP_SW_ALL` | `KPP_SW_MXL` | `KPP_SW_LV1`.

      integer :: i, j, k, nx, ny, nz
      real(wp) :: tau_mag, u_star
      real(wp) :: B_0, wstar3, w_star, w_s
      real(wp) :: u_ref, v_ref, b_ref, d_centre_ref
      real(wp) :: u_k, v_k, b_k, d_centre_k, d_running
      real(wp) :: shear2, ri_bulk, ri_prev, d_prev
      real(wp) :: delta_d, frac, denom, h_b
      real(wp) :: sigma, g_shape, kv_kpp, d_face_k
      real(wp) :: q_T_kin, q_S_kin, gamma_factor
      real(wp) :: h_b_lagged, wstar3_lagged, w_s_col
      real(wp) :: delta_b, n_brunt, vt2
      real(wp) :: q_bl
      logical :: crossing_found
      logical :: destabilizing

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml

      ! Surface buoyancy flux B_0 is computed PER COLUMN inside the
      ! BL-depth loop from the 2D flux fields (A7): identical arithmetic
      ! per column to the old host scalar when the fields are constant
      ! (bitwise), spatially varying when data-override fills them.

      ! Zero γ before re-populating.
      do concurrent(k=1:nz + 1, j=1:ny, i=1:nx)
         this%gamma_t(i, j, k) = 0.0_wp
         this%gamma_s(i, j, k) = 0.0_wp
      end do

      ! ---- Pass 1: per-column BL depth ----
      do concurrent(j=1:ny, i=1:nx) &
         local(k, u_ref, v_ref, b_ref, d_centre_ref, &
               u_k, v_k, b_k, d_centre_k, &
               shear2, ri_bulk, ri_prev, d_prev, &
               delta_d, frac, denom, h_b, crossing_found, &
               tau_mag, u_star, &
               h_b_lagged, wstar3_lagged, w_s_col, &
               delta_b, n_brunt, vt2, B_0, destabilizing, q_T_kin, q_S_kin, q_bl)
         u_ref = 0.5_wp*(ms%u_face_x_layer(i, j, nz) + ms%u_face_x_layer(i + 1, j, nz))
         v_ref = 0.5_wp*(ms%v_face_y_layer(i, j, nz) + ms%v_face_y_layer(i, j + 1, nz))
         b_ref = -GRAVITY*ms%rho_layer(i, j, nz)/this%rho0
         d_centre_ref = 0.5_wp*ms%h_layer(i, j, nz)

         ! PR-12 dedup: |tau| at cell centres is now a shared field
         ! (ocean_surface_stress_set_derived, same 3-line FP op order as
         ! the inline computation this replaces — bit-identical, §7.5).
         tau_mag = ss%stress_mag(i, j)
         u_star = sqrt(tau_mag/this%rho0)
         h_b_lagged = this%bl_depth(i, j)
         ! B_0 charges only the SW absorbed inside the (lagged) BL depth:
         !   Q_bl = Q_heat - I0·T(h_b)  (MXL_SW).  Gated on sw_active so
         !   sw_pen_frac=0 keeps the unmodified legacy source line.
         if (sw_active) then
            q_bl = sf%Q_heat(i, j)
            if (sw_method == KPP_SW_MXL) then
               q_bl = q_bl - sw_pen_frac*sw_src(i, j)* &
                      sw_transmission(h_b_lagged, sw_R, sw_zeta1, sw_zeta2)
            else if (sw_method == KPP_SW_LV1) then
               q_bl = q_bl - sw_pen_frac*sw_src(i, j)* &
                      sw_transmission(ms%h_layer(i, j, nz), sw_R, sw_zeta1, sw_zeta2)
            end if
            q_T_kin = q_bl/(this%rho0*sf%cp)
         else
            q_T_kin = sf%Q_heat(i, j)/(this%rho0*sf%cp)
         end if
         q_S_kin = sf%Q_salt(i, j)/this%rho0
         B_0 = kpp_surface_buoyancy_flux(this%eos%alpha_T, this%eos%beta_S, &
                                         this%rho0, q_T_kin, q_S_kin)
         destabilizing = (B_0 < 0.0_wp)
         wstar3_lagged = max(0.0_wp, -B_0)*h_b_lagged
         w_s_col = sqrt(u_star*u_star + wstar3_lagged**(2.0_wp/3.0_wp))

         ri_prev = 0.0_wp
         d_prev = d_centre_ref
         h_b = 0.0_wp
         crossing_found = .false.

         d_centre_k = d_centre_ref
         do k = nz - 1, 1, -1
            if (.not. crossing_found) then
               d_centre_k = d_centre_k + 0.5_wp*(ms%h_layer(i, j, k + 1) + ms%h_layer(i, j, k))
               u_k = 0.5_wp*(ms%u_face_x_layer(i, j, k) + ms%u_face_x_layer(i + 1, j, k))
               v_k = 0.5_wp*(ms%v_face_y_layer(i, j, k) + ms%v_face_y_layer(i, j + 1, k))
               b_k = -GRAVITY*ms%rho_layer(i, j, k)/this%rho0

               delta_d = d_centre_k - d_centre_ref
               delta_b = b_ref - b_k

               n_brunt = sqrt(max(0.0_wp, delta_b/delta_d))
               vt2 = (this%c_vt2/this%ri_crit)*d_centre_k*n_brunt*w_s_col

               shear2 = max((u_k - u_ref)*(u_k - u_ref) + (v_k - v_ref)*(v_k - v_ref) + vt2, &
                            this%shear2_floor)
               ri_bulk = delta_b*delta_d/shear2

               if (ri_bulk >= this%ri_crit) then
                  denom = ri_bulk - ri_prev
                  if (abs(denom) > 1.0e-12_wp) then
                     frac = (this%ri_crit - ri_prev)/denom
                  else
                     frac = 0.5_wp
                  end if
                  h_b = d_prev + frac*(d_centre_k - d_prev)
                  crossing_found = .true.
               end if

               ri_prev = ri_bulk
               d_prev = d_centre_k
            end if
         end do

         if (.not. crossing_found) then
            h_b = 0.0_wp
            do k = 1, nz
               h_b = h_b + ms%h_layer(i, j, k)
            end do
         end if

         this%bl_depth(i, j) = h_b
      end do

      ! ---- Pass 2: overlay kv_kpp inside the BL ----
      do concurrent(j=1:ny, i=1:nx) &
         local(k, d_running, d_face_k, sigma, g_shape, kv_kpp, h_b, &
               tau_mag, u_star, &
               wstar3, w_star, w_s, gamma_factor, B_0, destabilizing, &
               q_T_kin, q_S_kin, q_bl)
         ! PR-12 dedup — see Pass 1's comment above.
         tau_mag = ss%stress_mag(i, j)
         u_star = sqrt(tau_mag/this%rho0)
         h_b = this%bl_depth(i, j)
         ! B_0 charges only the SW absorbed inside the (current) BL depth
         !   Q_bl = Q_heat - I0·T(h_b)  (MXL_SW).  See Pass 1.  Gated on
         !   sw_active so sw_pen_frac=0 keeps the legacy source line.
         if (sw_active) then
            q_bl = sf%Q_heat(i, j)
            if (sw_method == KPP_SW_MXL) then
               q_bl = q_bl - sw_pen_frac*sw_src(i, j)* &
                      sw_transmission(h_b, sw_R, sw_zeta1, sw_zeta2)
            else if (sw_method == KPP_SW_LV1) then
               q_bl = q_bl - sw_pen_frac*sw_src(i, j)* &
                      sw_transmission(ms%h_layer(i, j, nz), sw_R, sw_zeta1, sw_zeta2)
            end if
            q_T_kin = q_bl/(this%rho0*sf%cp)
         else
            q_T_kin = sf%Q_heat(i, j)/(this%rho0*sf%cp)
         end if
         q_S_kin = sf%Q_salt(i, j)/this%rho0
         B_0 = kpp_surface_buoyancy_flux(this%eos%alpha_T, this%eos%beta_S, &
                                         this%rho0, q_T_kin, q_S_kin)
         this%b0(i, j) = B_0
         destabilizing = (B_0 < 0.0_wp)
         wstar3 = max(0.0_wp, -B_0)*h_b
         w_star = wstar3**(1.0_wp/3.0_wp)
         w_s = sqrt(u_star*u_star + w_star*w_star)
         if (h_b > 0.0_wp .and. w_s > 0.0_wp) then
            d_running = 0.0_wp
            do k = nz, 1, -1
               d_running = d_running + ms%h_layer(i, j, k)
               d_face_k = d_running
               if (k > 1 .and. d_face_k < h_b) then
                  sigma = d_face_k/h_b
                  g_shape = sigma*(1.0_wp - sigma)*(1.0_wp - sigma)
                  kv_kpp = h_b*w_s*g_shape
                  this%kv(i, j, k) = max(this%kv(i, j, k), kv_kpp)
                  this%kt(i, j, k) = max(this%kt(i, j, k), kv_kpp)
                  if (destabilizing) then
                     gamma_factor = this%cs_nonlocal*g_shape
                     this%gamma_t(i, j, k) = gamma_factor*q_T_kin
                     this%gamma_s(i, j, k) = gamma_factor*q_S_kin
                  end if
               end if
            end do
         end if
      end do
   end subroutine vmix_kpp_overlay_impl

   pure function kpp_surface_buoyancy_flux(alpha_T, beta_S, rho0, q_T_kin, q_S_kin) &
      result(b0)
      !! Surface buoyancy flux for the KPP overlay's convective scale:
      !!
      !!   `B_0 = (g/ρ_0)·(α_T·F_T − β_S·F_S)`   [m^2/s^3]
      !!
      !! **Units convention.** `alpha_T` / `beta_S` are the LINEAR-EOS
      !! DIMENSIONAL sensitivities of the density ANOMALY form used
      !! throughout this tree (`rdb_eos`):
      !!
      !!   `rho = rho_0 + beta_S·(S − S_ref) − alpha_T·(T − T_ref)`
      !!
      !! so `alpha_T = -∂ρ/∂T` [kg m^-3 K^-1] and
      !! `beta_S = ∂ρ/∂S` [kg m^-3 psu^-1] — NOT the fractional
      !! `(1/ρ)·∂ρ/∂T` coefficients (~2e-4 K^-1) that many texts write
      !! as α.  Buoyancy is `b = −g·ρ'/ρ_0`, so converting a dimensional
      !! `α_T` into a buoyancy flux costs a `1/ρ_0` — dropping it makes
      !! `B_0` a factor `ρ_0` (~1035) too large.  Fed a FRACTIONAL
      !! coefficient `α_frac = α_T/ρ_0` the caller must therefore pass
      !! `rho0 = 1`, and the two spellings agree exactly.
      !!
      !! `q_T_kin` / `q_S_kin` are the KINEMATIC surface fluxes
      !! `Q_heat/(ρ_0·cp)` [K m/s] and `Q_salt/ρ_0` [psu m/s] (the
      !! caller charges `q_T_kin` for penetrating shortwave first, so
      !! this stays a pure algebraic kernel).  `b0 > 0` is stabilizing
      !! (heating / freshening), `b0 < 0` destabilizing.
      !!
      !! Same quantity as EPBL's `b0 = g·ρ_0·(dSV/dT·q_T + dSV/dS·q_S)`:
      !! the linear EOS has `dSV/dT = +α_T/ρ_0²`, `dSV/dS = −β_S/ρ_0²`,
      !! so the two reduce to the identical expression (they agree to
      !! round-off, not bitwise — the FP op orders differ).
      !$acc routine seq
      real(wp), intent(in) :: alpha_T, beta_S
         !! Dimensional linear-EOS sensitivities (kg/m^3 per degC / psu).
      real(wp), intent(in) :: rho0
         !! Boussinesq reference density (kg/m^3).
      real(wp), intent(in) :: q_T_kin, q_S_kin
         !! Kinematic surface heat / salt fluxes (K m/s, psu m/s).
      real(wp) :: b0

      b0 = (GRAVITY/rho0)*(alpha_T*q_T_kin - beta_S*q_S_kin)
   end function kpp_surface_buoyancy_flux

   pure function parse_kpp_sw_method(name) result(tag)
      !! Map the `&ocean_thermo_nml kpp_sw_method` string to a
      !! `KPP_SW_*` tag; `-1` for an unrecognised string (fail-loud at
      !! `validate_config`).
      character(len=*), intent(in) :: name
      integer :: tag
      select case (trim(name))
      case ("all_sw")
         tag = KPP_SW_ALL
      case ("mxl_sw")
         tag = KPP_SW_MXL
      case ("lv1_sw")
         tag = KPP_SW_LV1
      case default
         tag = -1
      end select
   end function parse_kpp_sw_method

   pure function kpp_sw_method_is_implemented(name) result(ok)
      !! Fail-loud predicate for `kpp_sw_method` — `validate_config`
      !! aborts on any string this rejects.
      character(len=*), intent(in) :: name
      logical :: ok
      ok = (parse_kpp_sw_method(name) > 0)
   end function kpp_sw_method_is_implemented

   pure subroutine vmix_apply_nonlocal_tendencies(grid, this, ms, dt)
      !! Apply the KPP non-local (counter-gradient) tracer tendency
      !! computed by `vmix_apply_kpp_overlay`.  Updates
      !! `hT`, `hS` from the divergence of `gamma_t` / `gamma_s` at
      !! interfaces:
      !!
      !!   hT(k) += dt · (γ_T(k+1) - γ_T(k))
      !!
      !! (interface k+1 sits above layer k, interface k below — the
      !! ROMS-style k=1-bed convention).  γ at the bed (k=1) and
      !! free surface (k=nz+1) is identically zero, so the surface
      !! and bed layers see only the lower / upper neighbour's γ.
      !!
      !! No-op when γ is all zero (stabilizing flux or KPP
      !! disabled).  Touches only the layers with index_salinity /
      !! index_temperature in the tracer registry — extra tracers
      !! are diffusive only under KPP.
      type(hgrid_t), intent(in) :: grid
      type(ocean_vmix_t), intent(in) :: this
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: dt

      integer :: i, j, k, nx, ny, nz, it_T, it_S, it_ps

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml
      it_T = ms%idx_temperature
      it_S = ms%idx_salinity
      it_ps = ms%idx_pseudo_salt

      if (it_T > 0) then
         do concurrent(k=1:nz, j=1:ny, i=1:nx)
            ms%tracers(it_T)%hTr(i, j, k) = ms%tracers(it_T)%hTr(i, j, k) + &
                                            dt*(this%gamma_t(i, j, k + 1) - &
                                                this%gamma_t(i, j, k))
         end do
      end if
      if (it_S > 0) then
         do concurrent(k=1:nz, j=1:ny, i=1:nx)
            ms%tracers(it_S)%hTr(i, j, k) = ms%tracers(it_S)%hTr(i, j, k) + &
                                            dt*(this%gamma_s(i, j, k + 1) - &
                                                this%gamma_s(i, j, k))
         end do
      end if
      ! Pseudo-salt mirror: identical gamma_s expression, identical dt
      ! (§5.4.2) — `it_ps` is a HOST scalar hoisted above, exactly like
      ! `it_T`/`it_S`, so this is legal inside `do concurrent`.
      if (it_ps > 0) then
         do concurrent(k=1:nz, j=1:ny, i=1:nx)
            ms%tracers(it_ps)%hTr(i, j, k) = ms%tracers(it_ps)%hTr(i, j, k) + &
                                             dt*(this%gamma_s(i, j, k + 1) - &
                                                 this%gamma_s(i, j, k))
         end do
      end if
   end subroutine vmix_apply_nonlocal_tendencies

   pure subroutine vmix_add_kv_ml_invz2(grid, this, ms)
      !! Augment `this%kv` with an extra near-surface viscosity:
      !!
      !!     kv_extra(z) = kv_ml_invz2 · (hmix_fixed / max(z, dz_min))²
      !!
      !! where `z` is the depth (m) of the layer interface below the
      !! free surface.  Active only for interfaces whose `z <
      !! hmix_fixed`.  Adds to the existing `kv` (which `kt` does not
      !! see — momentum-only).
      !!
      !! Interface convention: `kv(:, :, k)` lives at the bottom of
      !! layer k.  Surface-most interior interface is at `k = nz`
      !! (between layer nz and the wall at nz+1).  Surface boundary
      !! `k = nz+1` stays at zero (closed top), and we don't touch it.
      !!
      !! When `kv_ml_invz2 <= 0` the routine is a no-op so existing
      !! configs are bit-identical.
      type(hgrid_t), intent(in) :: grid
      type(ocean_vmix_t), intent(inout) :: this
      type(multilayer_state_t), intent(in) :: ms

      integer :: i, j, k, nx, ny, nz
      real(wp) :: kv_extra, hmix, z, dz_min_safe, hmix_sq

      if (this%kv_ml_invz2 <= 0.0_wp) return
      if (this%hmix_fixed <= 0.0_wp) return

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml
      hmix = this%hmix_fixed
      hmix_sq = hmix*hmix

      ! Floor on `z` to prevent the 1/z² profile from diverging right
      ! under the surface.  Half a target layer thickness is a safe,
      ! grid-resolved minimum.
      dz_min_safe = 0.5_wp*hmix/real(max(nz, 1), wp)

      ! Interior interfaces: k = nz (surface-most) down to k = 2.
      ! The depth of interface k below the surface is the cumulative
      ! thickness of layers nz, nz-1, ..., k.
      do concurrent(j=1:ny, i=1:nx) local(k, z, kv_extra)
         z = 0.0_wp
         do k = nz, 2, -1
            z = z + ms%h_layer(i, j, k)
            if (z >= hmix) exit
            kv_extra = this%kv_ml_invz2*hmix_sq/(max(z, dz_min_safe)**2)
            this%kv(i, j, k) = this%kv(i, j, k) + kv_extra
         end do
      end do
   end subroutine vmix_add_kv_ml_invz2

   pure subroutine vmix_apply_convection(grid, this, ms, bl_depth)
      !! Brunt-Vaisala-triggered convective adjustment -- a CONTRIBUTOR
      !! into kv/kt (see the module docstring for the full physics and
      !! the D1-D4 divergences from MOM6's `MOM_CVMix_conv`).  Early-
      !! returns when `conv_enable` is off (bit-identical).  `bl_depth`
      !! (m, positive down) is the caller's active surface-boundary-layer
      !! depth: pass `this%bl_depth` under KPP / no BL scheme, `epbl%mld`
      !! under EPBL -- both are permanently-zero, device-resident fields
      !! when their owning scheme is off, so either is a free "no BL
      !! scheme -> mask nothing" default (`z_int(k) >= 0` for every
      !! interior interface).
      type(hgrid_t), intent(in) :: grid
      type(ocean_vmix_t), intent(inout) :: this
      type(multilayer_state_t), intent(in) :: ms
      real(wp), intent(in) :: bl_depth(:, :)

      integer :: nx, ny, nz

      if (.not. this%conv_enable) return

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml

      call vmix_convection_impl(nx, ny, nz + 1, this%kv, this%kt, &
                                ms%rho_layer, ms%h_layer, bl_depth, &
                                this%conv_kd, this%conv_prandtl, &
                                this%conv_n2_thresh, this%rho0)
   end subroutine vmix_apply_convection

   pure subroutine vmix_convection_impl(nx, ny, nzp1, kv, kt, rho_layer, h_layer, &
                                        bl_depth, kd_conv, prandtl_conv, n2_thresh, rho0)
      !! Flat, explicit-shape kernel (the "vmix incident" module --
      !! assumed-shape dummies here produced 1.4M per-launch descriptor-
      !! walk memcpys; explicit-shape only, no exceptions).  Interior
      !! interfaces k = 2..nzp1-1 only; boundary interfaces k=1 (bed) and
      !! k=nzp1 (surface) stay at the closed-BC zero and are never
      !! touched (same contract as PP81 / `vmix_add_kv_ml_invz2`).
      !!
      !! Depth walk copied verbatim from `vmix_add_kv_ml_invz2`: `z`
      !! accumulates the cumulative thickness of layers nz, nz-1, ...
      !! down to interface k -- the depth of interface k below the free
      !! surface.  Masked out (left at the upstream PP81 value) above
      !! the active BL (`z < bl_depth`); the BL owns its own convective
      !! response (KPP w_* / EPBL convective TKE) and double-counting
      !! would over-mix and defeat that scheme's own energetics.
      !!
      !! N^2 expression is bit-for-bit PP81's (`vmix_compute_pp81`):
      !! `n2 = -g*(rho_k - rho_km1)/(rho0*dz_face)`, positive when
      !! stable.  `dz_face` is gated on `H_VANISHED` (dynamic-vanish
      !! skip), NOT `H_DIV_EPS`/bare `<= 0` -- a pinched ZSTAR_FULL
      !! interface must be skipped, not divided by (a `0/0` there would
      !! otherwise be silently laundered to background by
      !! `vmix_assemble`'s default un-guarded clip).
      integer, intent(in) :: nx, ny, nzp1
      real(wp), intent(inout) :: kv(nx, ny, nzp1)
      real(wp), intent(inout) :: kt(nx, ny, nzp1)
      real(wp), intent(in) :: rho_layer(nx, ny, nzp1 - 1)
      real(wp), intent(in) :: h_layer(nx, ny, nzp1 - 1)
      real(wp), intent(in) :: bl_depth(nx, ny)
      real(wp), intent(in) :: kd_conv, prandtl_conv, n2_thresh, rho0

      integer :: i, j, k, nz
      real(wp) :: dz_face, n2, z

      nz = nzp1 - 1

      do concurrent(j=1:ny, i=1:nx) local(k, z, dz_face, n2)
         z = 0.0_wp
         do k = nz, 2, -1
            z = z + h_layer(i, j, k)
            dz_face = 0.5_wp*(h_layer(i, j, k - 1) + h_layer(i, j, k))
            if (dz_face <= H_VANISHED) cycle
            if (z < bl_depth(i, j)) cycle
            n2 = -GRAVITY*(rho_layer(i, j, k) - rho_layer(i, j, k - 1))/(rho0*dz_face)
            if (n2 < n2_thresh) then
               kt(i, j, k) = max(kt(i, j, k), kd_conv)
               kv(i, j, k) = max(kv(i, j, k), prandtl_conv*kd_conv)
            end if
         end do
      end do
   end subroutine vmix_convection_impl
   pure subroutine vmix_split_kd_heat_salt(grid, this, ms)
      !! Analogue of MOM6's heat/salt diffusivity split.  Derives the
      !! per-tracer diffusivities from the assembled interior/boundary-
      !! layer diffusivity:
      !!     Kd_heat = Kd_int + Kd_extra_T  ->  kt
      !!     Kd_salt = Kd_int + Kd_extra_S  ->  ks
      !! `kt` holds Kd_int on entry (every contributor writes it).  No
      !! double-diffusion contributor exists yet, so Kd_extra_{T,S} = 0
      !! and the split reduces to `ks := kt` => bit-identical.  PR-33
      !! extends this to
      !!     ks = kt + kd_extra_s ;  kt = kt + kd_extra_t
      !! BOTH computed from the SAME pre-split kt -- read kt into a
      !! local before writing it, or the kt update poisons the ks
      !! update.
      !!
      !! `vmix_split_kd_heat_salt` is not a contributor -- it is the
      !! last statement before `vmix_assemble`, always.  Any PR adding
      !! a kv/kt contributor inserts it ABOVE this split.
      !!
      !! Runs over the full k = 1..nz+1 range (not 2..nz): `kt`'s
      !! boundary interfaces are zeroed every stage by
      !! `vmix_compute_pp81`, so an all-k copy reproduces `ks`'s
      !! init-time boundary zeros exactly and keeps the closed-BC
      !! invariant true by construction.
      type(hgrid_t), intent(in) :: grid
      type(ocean_vmix_t), intent(inout) :: this
      type(multilayer_state_t), intent(in) :: ms

      integer :: nx, ny, nz

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml

      if (this%ddiff_enable) then
         ! Outer shim: dereference the registry tracers (idx 1/2) on the
         ! host and hand the flat top-level arrays to the do-concurrent
         ! impl (array-of-derived-type device indirection rule).
         call vmix_split_ddiff_impl(nx, ny, nz + 1, this%kt, this%ks, &
                                    ms%tracers(ms%idx_temperature)%hTr, &
                                    ms%tracers(ms%idx_salinity)%hTr, ms%h_layer, &
                                    this%eos%alpha_T, this%eos%beta_S, &
                                    this%ddiff_strat_param_max, this%ddiff_kappa_s, &
                                    this%ddiff_exp1, this%ddiff_exp2, &
                                    this%ddiff_param1, this%ddiff_param2, &
                                    this%ddiff_param3, this%ddiff_mol_diff, &
                                    this%ddiff_use_k90)
      else
         call vmix_split_kd_heat_salt_impl(nx, ny, nz + 1, this%kt, this%ks)
      end if
   end subroutine vmix_split_kd_heat_salt

   pure subroutine vmix_split_kd_heat_salt_impl(nx, ny, nzp1, kt, ks)
      !! Explicit-shape args so the do concurrent stays descriptor-walk
      !! free (assumed-shape dummies in a do concurrent make NVHPC walk
      !! descriptors per launch -- this runs every stage).
      integer, intent(in) :: nx, ny, nzp1
      real(wp), intent(in) :: kt(nx, ny, nzp1)
      real(wp), intent(inout) :: ks(nx, ny, nzp1)

      integer :: i, j, k

      do concurrent(k=1:nzp1, j=1:ny, i=1:nx)
         ks(i, j, k) = kt(i, j, k)
      end do
   end subroutine vmix_split_kd_heat_salt_impl

   pure subroutine vmix_split_ddiff_impl(nx, ny, nzp1, kt, ks, temp_h, salt_h, &
                                         h_layer, alpha_T, beta_S, strat_param_max, &
                                         kappa_s, exp1, exp2, param1, param2, param3, &
                                         mol_diff, use_k90)
      !! Double-diffusion split.  Replaces `ks := kt` with the asymmetric
      !!   ks = kt_pre + kd_extra_s ;  kt = kt_pre + kd_extra_t
      !! both from the SAME pre-split kt (read into `kt_pre` before either
      !! write).  CVMix `cvmix_coeffs_ddiff` algebra (Large et al. 1994
      !! fingering; Marmorino-Caldwell 1976 / Kelley 1990 convection).
      !!
      !! Interior interfaces k=2..nz only (boundary interfaces k=1/nzp1
      !! carry no fingering and stay at kt's closed-BC zero, reproducing
      !! `ks`'s boundary zeros -- same all-k structure as the plain split).
      !! Interface k sits at the bottom of layer k, between layer k-1
      !! (lower) and layer k (upper, toward-surface): bottom-up convention.
      !!
      !! Branch on the SIGNED alpha*dT / beta*dS (never a pre-divided
      !! R_rho) so R_rho<=0, R_rho->inf and 0/0 never arise; R_rho is
      !! formed only inside a branch whose denominator is provably nonzero.
      !! No floor/ceiling here -- `vmix_assemble` is the single gate.
      !! Explicit-shape dummies (per the "vmix incident" descriptor-walk
      !! rule); linear-EOS constant alpha_T/beta_S passed by value.
      integer, intent(in) :: nx, ny, nzp1
      real(wp), intent(inout) :: kt(nx, ny, nzp1)
      real(wp), intent(inout) :: ks(nx, ny, nzp1)
      real(wp), intent(in) :: temp_h(nx, ny, nzp1 - 1)
         !! Temperature thickness-integral hTr = T*h_layer (degC*m).
      real(wp), intent(in) :: salt_h(nx, ny, nzp1 - 1)
         !! Salinity thickness-integral hTr = S*h_layer (PSU*m).
      real(wp), intent(in) :: h_layer(nx, ny, nzp1 - 1)
      real(wp), intent(in) :: alpha_T, beta_S
         !! |d rho/dT| and d rho/dS (kg/m^3 per degC / per PSU).
      real(wp), intent(in) :: strat_param_max, kappa_s, exp1, exp2
      real(wp), intent(in) :: param1, param2, param3, mol_diff
      logical, intent(in) :: use_k90

      integer :: i, j, k, nz
      real(wp) :: kt_pre, adT, bdS, rrho, ddiff, kd_t, kd_s, hu, hl

      nz = nzp1 - 1

      do concurrent(k=1:nzp1, j=1:ny, i=1:nx) &
         local(kt_pre, adT, bdS, rrho, ddiff, kd_t, kd_s, hu, hl)
         kt_pre = kt(i, j, k)
         kd_t = 0.0_wp
         kd_s = 0.0_wp

         if (k >= 2 .and. k <= nz) then
            hu = h_layer(i, j, k)      ! upper layer (toward surface)
            hl = h_layer(i, j, k - 1)  ! lower layer
            if (hu > H_VANISHED .and. hl > H_VANISHED) then
               ! alpha*dT and beta*dS across interface k (upper - lower)
               adT = alpha_T*(temp_h(i, j, k)/hu - temp_h(i, j, k - 1)/hl)
               bdS = beta_S*(salt_h(i, j, k)/hu - salt_h(i, j, k - 1)/hl)

               if (adT >= bdS .and. bdS > 0.0_wp) then
                  ! ---- salt fingering (R_rho >= 1) ----
                  rrho = adT/bdS
                  if (rrho < strat_param_max) then
                     ddiff = (1.0_wp - ((rrho - 1.0_wp)/ &
                                        (strat_param_max - 1.0_wp))**exp1)**exp2
                     kd_s = kappa_s*ddiff
                  end if
                  kd_t = 0.7_wp*kd_s
               else if (adT >= bdS .and. adT < 0.0_wp) then
                  ! ---- diffusive convection (0 < R_rho < 1) ----
                  rrho = adT/bdS
                  if (use_k90) then
                     ddiff = mol_diff*8.7_wp*rrho**1.1_wp
                  else
                     ddiff = mol_diff*param1* &
                             exp(param2*exp(param3*(1.0_wp/rrho - 1.0_wp)))
                  end if
                  kd_t = ddiff
                  if (rrho < 0.5_wp) then
                     kd_s = 0.15_wp*rrho*ddiff
                  else
                     kd_s = (1.85_wp*rrho - 0.85_wp)*ddiff
                  end if
               end if
            end if
         end if

         ks(i, j, k) = kt_pre + kd_s
         kt(i, j, k) = kt_pre + kd_t
      end do
   end subroutine vmix_split_ddiff_impl

   subroutine vmix_assemble(grid, this, ms, geolat, status)
      !! The single downstream gate of the vmix diffusivity assembly —
      !! the `set_diffusivity`-style stage that every interior and
      !! overlay closure feeds into before vdiff consumes kv/kt/ks.
      !!
      !! CONTRACT.  Every interior / overlay closure CONTRIBUTES into
      !! `kv` / `kt` / `ks` upstream of this call: PP81 writes the
      !! interior background + Richardson term, KPP / EPBL overlay (max
      !! or add), kappa-shear merges additively, KV_ML_INVZ2 augments
      !! the surface band.  `vmix_assemble` is the SINGLE place the
      !! assembled fields are gated.  Future Area-C contributors (tidal
      !! mixing, double-diffusion, geothermal-adjacent floors, depth-
      !! varying background profiles) plug in as additional upstream
      !! contributors — they do NOT add their own floor/clip; they rely
      !! on this stage.  Applied here, in order, on interior interfaces
      !! k = 2..nz (the boundary interfaces k=1/bed and k=nz+1/surface
      !! stay at the closed-BC zero, untouched):
      !!   (a) optional debug-gated negative/NaN guard (`vmix_guard`,
      !!       default off) — runs FIRST, on the raw closure output before
      !!       any floor or clip.  The guard must see raw values because the
      !!       subsequent clip would launder NaN/negatives: NVHPC -O2
      !!       evaluates min(max(NaN, bg), huge) to bg, converting NaN
      !!       diffusivities into silent plausible mixing.  With `status`
      !!       present the routine returns a non-zero code (testable path);
      !!       without it a tripped guard `error stop`s.
      !!   (b) background floors — `kv ≥ kv_bg`, `kt ≥ kt_bg`,
      !!       `ks ≥ ks_bg`.  Defaults match `pp81_nu_bg` / `pp81_kappa_bg`
      !!       (structural invariant set in `init`) so the floor is a no-op
      !!       for the shipped closure path (bit-identical).  This is the
      !!       single place the constant background is enforced going forward.
      !!       Two MUTUALLY EXCLUSIVE opt-in variants replace the scalar
      !!       floor: `bkgnd_profile` (Bryan-Lewis per-interface `kd_bg`
      !!       field, kv floor = `bkgnd_prandtl·kd_bg`) and `bkgnd_henyey`
      !!       (Henyey latitude factor on the scalar tracer floors,
      !!       `max(bkgnd_kd_min, kt_bg·L(phi))`).  Configure refuses both
      !!       at once, so this dispatch is a three-way `if/else if/else`.
      !!   (c) ceilings — `kv ≤ kv_max`, `kt,ks ≤ kd_max`.  Defaults
      !!       `huge(1.0)` ⇒ no clip ⇒ bit-identical.
      !!   (d) optional 1-2-1 horizontal smoothing of kv/kt/ks at interfaces
      !!       (`kd_smooth_iterations`, default 0 = off).  Wet-mask aware:
      !!       contributions from dry neighbours are excluded and the stencil
      !!       weight is renormalised over wet cells only; a dry centre column
      !!       is left unchanged.  NOTE: under future MPI, multi-pass
      !!       smoothing requires nghost >= kd_smooth_iterations and a halo
      !!       refresh between passes; with nghost = 2 only 1-2 passes are
      !!       safe without the refresh.
      !!
      !! `ks` is derived from `kt` by `vmix_split_kd_heat_salt`, called
      !! upstream of this gate (after the last kv/kt contributor), and
      !! consumed by `vdiff_apply_tracers` (salinity + every passive
      !! tracer).  It equals `kt` until a double-diffusion contributor
      !! lands (PR-33).
      type(hgrid_t), intent(in) :: grid
      type(ocean_vmix_t), intent(inout) :: this
      type(multilayer_state_t), intent(in) :: ms
      real(wp), intent(in), optional :: geolat(grid%nx_total, grid%ny_total)
         !! T-point geographic latitude (degN, `ocean_metrics_t%geolatT` —
         !! the T stagger specifically; `geolatBu` is corner-shaped
         !! `(nx+1, ny+1)` and is NOT interchangeable here).  Explicit-shape
         !! so the value flows into `vmix_assemble_clip_henyey_impl`'s
         !! `do concurrent` without a descriptor walk.  Only read when
         !! `bkgnd_henyey`.  Required in that case —
         !! absent then `error stop`s (a caller forgot to thread `metrics`
         !! through).  Every production call site (`vmix_apply_in_stage`)
         !! has `metrics` in scope and always passes it; test harnesses that
         !! never enable `bkgnd_henyey` may omit it.
      integer, intent(out), optional :: status
         !! 0 = ok; 1 = guard tripped (negative or NaN K).  Only
         !! written when `vmix_guard` is on.  When absent and the guard
         !! trips, the routine `error stop`s instead.

      integer :: nx, ny, nz, it
      integer :: bad_count

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml

      if (present(status)) status = 0

      ! (a) optional negative/NaN guard — BEFORE clip so laundering
      ! (NaN→bg via min/max) cannot suppress a real closure error.
      if (this%vmix_guard) then
         call vmix_guard_impl(nx, ny, nz + 1, this%kv, this%kt, this%ks, bad_count)
         if (bad_count > 0) then
            if (present(status)) then
               status = 1
            else
               error stop "vmix_assemble: negative or NaN diffusivity detected"
            end if
         end if
      end if

      ! (b)+(c) floors + ceilings on interior interfaces.  One DC over
      ! (i,j,k) with explicit-shape args via the _impl helper (assumed-
      ! shape dummies in a do concurrent make NVHPC walk descriptors
      ! per launch).  Two floor paths:
      !   * default (bkgnd_profile off): the SCALAR kv_bg/kt_bg/ks_bg
      !     floor exactly as before — bit-identical.
      !   * C7 Bryan-Lewis on: fill the per-interface kd_bg depth profile
      !     from the current column thicknesses (correct under any vcoord)
      !     and floor with that field; kv floor = bkgnd_prandtl * kd_bg.
      if (this%bkgnd_profile) then
         call vmix_bkgnd_fill_impl(nx, ny, nz + 1, this%kd_bg, ms%h_layer, &
                                   this%bkgnd_kd_sfc, this%bkgnd_kd_deep, &
                                   this%bkgnd_z0, this%bkgnd_delta)
         call vmix_assemble_clip_profile_impl(nx, ny, nz + 1, this%kv, this%kt, this%ks, &
                                              this%kd_bg, this%bkgnd_prandtl, &
                                              this%kv_max, this%kd_max)
      else if (this%bkgnd_henyey) then
         ! C7 Henyey: the latitude factor scales the SCALAR kt_bg/ks_bg
         ! floors, floored at `bkgnd_kd_min` — see the `bkgnd_henyey`
         ! docstring + `henyey_lat_factor_impl` for the exact form.  This
         ! branch is unreachable together with `bkgnd_profile`: the two are
         ! mutually exclusive (configure refuses both), so the factor never
         ! touches the Bryan-Lewis deep asymptote.
         !
         ! `geolat` is required here — a caller that turned this knob on but
         ! never threaded `metrics%geolatT` through is a wiring bug, not a
         ! silent no-op.  Deliberately a SEPARATE `_impl` from the plain
         ! scalar clip (rather than one routine with optional Henyey args):
         ! with no `present()` inside the kernel there is no way for a
         ! partially-supplied argument list to degrade into a silent no-op,
         ! and the default path's routine body is literally unchanged from
         ! before Henyey existed.
         if (.not. present(geolat)) then
            error stop "vmix_assemble: bkgnd_henyey requires the geolat "// &
               "argument (thread ocean_metrics_t%geolatT through the caller)"
         end if
         call vmix_assemble_clip_henyey_impl(nx, ny, nz + 1, this%kv, this%kt, this%ks, &
                                             this%kv_bg, this%kt_bg, this%ks_bg, &
                                             this%kv_max, this%kd_max, &
                                             this%bkgnd_henyey_n0_2omega, &
                                             this%bkgnd_henyey_max_lat, &
                                             vmix_resolve_kd_min(this%bkgnd_kd_min, &
                                                                 this%kt_bg), &
                                             geolat)
      else
         call vmix_assemble_clip_impl(nx, ny, nz + 1, this%kv, this%kt, this%ks, &
                                      this%kv_bg, this%kt_bg, this%ks_bg, &
                                      this%kv_max, this%kd_max)
      end if

      ! (d) optional 1-2-1 horizontal smoothing of kv / kt.  scratch is a
      ! persistent slot buffer mapped in enter_data — no lazy attach.
      ! Wet-mask aware: dry neighbours are excluded, stencil renormalised.
      if (this%kd_smooth_iterations > 0) then
         do it = 1, this%kd_smooth_iterations
            call vmix_smooth_121_impl(nx, ny, nz + 1, this%kv, this%smooth_scratch, ms%wet_mask)
            call vmix_smooth_121_impl(nx, ny, nz + 1, this%kt, this%smooth_scratch, ms%wet_mask)
            call vmix_smooth_121_impl(nx, ny, nz + 1, this%ks, this%smooth_scratch, ms%wet_mask)
         end do
      end if
   end subroutine vmix_assemble

   pure subroutine vmix_assemble_clip_impl(nx, ny, nzp1, kv, kt, ks, &
                                           kv_bg, kt_bg, ks_bg, kv_max, kd_max)
      !! Floor + ceiling on interior interfaces k = 2..nzp1-1.  Explicit-
      !! shape args so the do concurrent stays descriptor-walk free.
      integer, intent(in) :: nx, ny, nzp1
      real(wp), intent(inout) :: kv(nx, ny, nzp1)
      real(wp), intent(inout) :: kt(nx, ny, nzp1)
      real(wp), intent(inout) :: ks(nx, ny, nzp1)
      real(wp), intent(in) :: kv_bg, kt_bg, ks_bg, kv_max, kd_max

      integer :: i, j, k

      do concurrent(k=2:nzp1 - 1, j=1:ny, i=1:nx)
         kv(i, j, k) = min(max(kv(i, j, k), kv_bg), kv_max)
         kt(i, j, k) = min(max(kt(i, j, k), kt_bg), kd_max)
         ks(i, j, k) = min(max(ks(i, j, k), ks_bg), kd_max)
      end do
   end subroutine vmix_assemble_clip_impl

   pure subroutine vmix_bkgnd_fill_impl(nx, ny, nzp1, kd_bg, h_layer, &
                                        kd_sfc, kd_deep, z0, delta)
      !! Fill the per-interface Bryan & Lewis (1979) JGR 84:2503 background
      !! tracer-diffusivity profile from the CURRENT column interface depths.
      !!
      !!   Kd_bg(z) = Kd_sfc + (Kd_deep - Kd_sfc)·[½ + atan((|z| - z0)/Δ)/π]
      !!
      !! z is the interface depth below the free surface, accumulated from
      !! the current `h_layer` (surface k=nzp1 → bed k=1) so the profile is
      !! correct under any vcoord (z*, sigma, …) and tracks the moving free
      !! surface.  Boundary interfaces k=1 (bed) and k=nzp1 (surface) are
      !! left untouched (the assembly clip only floors k=2..nzp1-1).
      !!
      !! There is no Henyey-scaled variant of this fill: Bryan-Lewis and the
      !! Henyey latitude factor are mutually exclusive (configure refuses
      !! both), matching the reference formulation.  The Henyey path scales
      !! the SCALAR background instead — see `vmix_assemble_clip_henyey_impl`.
      !!
      !! Explicit-shape args so the do concurrent stays descriptor-walk free.
      integer, intent(in) :: nx, ny, nzp1
      real(wp), intent(inout) :: kd_bg(nx, ny, nzp1)
      real(wp), intent(in) :: h_layer(nx, ny, nzp1 - 1)
      real(wp), intent(in) :: kd_sfc, kd_deep, z0, delta
         !! Surface / deep asymptotes (m^2/s), transition centre depth (m),
         !! transition half-width (m).

      integer :: i, j, k
      real(wp) :: depth, inv_pi

      inv_pi = 1.0_wp/PI

      ! Per column: walk the interfaces top-down (surface → bed),
      ! accumulating depth.  Interface k sits at the bottom of layer k;
      ! its depth below the surface is the sum of the thicknesses of
      ! layers k+1 .. nz (nz = nzp1-1) — i.e. depth grows as k decreases.
      do concurrent(j=1:ny, i=1:nx) local(depth, k)
         depth = 0.0_wp
         ! Interface k is the bottom of layer k; its depth below the
         ! surface is the cumulative thickness of layers k..nz.  Descending
         ! from nz, adding h_layer(k) before evaluating interface k gives
         ! exactly that running sum.
         do k = nzp1 - 1, 2, -1
            depth = depth + h_layer(i, j, k)
            kd_bg(i, j, k) = kd_sfc + (kd_deep - kd_sfc)* &
                             (0.5_wp + atan((depth - z0)/delta)*inv_pi)
         end do
      end do
   end subroutine vmix_bkgnd_fill_impl

   pure function henyey_lat_factor_impl(lat_deg, n0_2omega, max_lat) result(fac)
      !$acc routine seq
      !! Henyey, Wright & Flatte (1986) JGR 91:8487 latitude dependence of
      !! the internal-wave-driven mixing rate, in the SIMPLIFIED constant-`N0`
      !! form of Harrison & Hallberg (2008) JPO 38:1894 — the in-situ column
      !! stratification is replaced by a fixed reference `N0`, so the factor
      !! collapses to a pure function of latitude:
      !!
      !!   L(φ) = |sin φ| · acosh(N0_2Ω / max(|sin φ|, HENYEY_MIN_SINLAT))
      !!          / [ sin 30° · acosh(N0_2Ω / sin 30°) ]
      !!
      !! `N0_2Ω = N0/(2Ω)`.  `L(30°) = 1` to round-off by construction (the
      !! denominator IS the numerator at 30°), which is what makes 30°
      !! useless as a test latitude for "was the factor applied at all" —
      !! a test that cannot distinguish `×L` from `×1` there.  Use 45°/90°.
      !!
      !! Equator singularity: `|sin φ|` is floored to `HENYEY_MIN_SINLAT`
      !! ONLY inside the ratio that feeds `acosh` (that is the 1/0 guard);
      !! the OUTER multiplication keeps the TRUE (unfloored) `|sin φ|`, so
      !! `L → 0` smoothly at the equator instead of blowing up, and
      !! `L(0°) = 0` EXACTLY.
      !!
      !! Poleward of `max_lat` (degN, compared against `|lat_deg|` so BOTH
      !! hemispheres clamp) `|sin φ|` is reset to `HENYEY_MIN_SINLAT`
      !! everywhere in the expression, collapsing `L` to a tiny positive
      !! floor (~1.2e-9 at the default `n0_2omega`) — NOT to the exact zero
      !! the equator produces.  Inert at the default `max_lat = 95` (> 90,
      !! never trips for a real latitude).
      !!
      !! This function returns the RAW factor `L(phi)` — it applies no
      !! minimum-diffusivity floor, so `L(0°)` really is exactly 0 and the
      !! poleward clamp really does return ~1.2e-9.  The floor lives one
      !! level up, in `vmix_assemble_clip_henyey_impl`, as
      !! `max(kd_min, kt_bg·L(phi))` — matching the reference code, which
      !! wraps the scaled diffusivity in `max(Kd_min, Kd·L)` with `Kd_min`
      !! defaulting to `0.01·Kd` (`HENYEY_KD_MIN_FRAC`).  Keeping the floor
      !! out of the factor is what lets the unit tests assert the closed
      !! form of `L` directly against an `acosh` oracle.
      !!
      !! `pure` + `!$acc routine seq`: called per column from the
      !! `do concurrent` in `vmix_assemble_clip_henyey_impl`, and directly
      !! from the unit tests as the analytical oracle's subject.
      real(wp), intent(in) :: lat_deg
         !! T-point latitude, degrees north (may be negative).
      real(wp), intent(in) :: n0_2omega
         !! `N0/(2Ω)` reference-stratification ratio (nondim, ≥ 1).
      real(wp), intent(in) :: max_lat
         !! Poleward cutoff latitude (degN, compared against `|lat_deg|`).
      real(wp) :: fac

      real(wp) :: abs_sinlat

      ! ONE write site per local (a `local()` variable reassigned across a
      ! branch is the recorded gfortran-15.1 corruption shape) — `merge`
      ! evaluates both arms, and both are finite for every input: |sin| is
      ! bounded by 1 and n0_2omega ≥ 1, so the acosh argument is always ≥ 1.
      abs_sinlat = merge(HENYEY_MIN_SINLAT, abs(sin(lat_deg*DEG2RAD)), &
                         abs(lat_deg) > max_lat)
      ! sin(30 deg) = 0.5 exactly; the denominator normalises L(30 deg) = 1.
      fac = abs_sinlat*acosh(n0_2omega/max(HENYEY_MIN_SINLAT, abs_sinlat)) &
            /(0.5_wp*acosh(2.0_wp*n0_2omega))
   end function henyey_lat_factor_impl

   pure function bkgnd_henyey_conflicts_profile(bkgnd_henyey, bkgnd_profile) result(conflict)
      !! `.true.` when both background schemes are selected at once.
      !!
      !! Bryan-Lewis and Henyey are MUTUALLY EXCLUSIVE, matching the
      !! reference code, which calls its `check_bkgnd_scheme` guard once
      !! from each scheme's block and FATALs when a scheme was already
      !! selected.  Henyey scales the SCALAR background; composing it with
      !! the Bryan-Lewis depth profile would also scale the deep asymptote,
      !! which stands for abyssal / internal-tide mixing over rough
      !! topography — not the open-ocean internal-wave continuum Henyey
      !! describes, and not equatorially suppressed in nature.
      !!
      !! A `pure` predicate rather than an inline `.and.` in
      !! `validate_config` so the RULE is unit-testable: the guard sites
      !! themselves are `error stop` / logger paths a unit test cannot
      !! reach without a subprocess harness, which left the rule as a
      !! mutation-survivable hole.  Same shape as
      !! `pseudo_salt_conflicts_restore` / `lateral_closure_conflicts_smag_ah`.
      logical, intent(in) :: bkgnd_henyey
         !! `&ocean_vmix_nml bkgnd_henyey`.
      logical, intent(in) :: bkgnd_profile
         !! `&ocean_vmix_nml bkgnd_profile` (Bryan-Lewis).
      logical :: conflict

      conflict = bkgnd_henyey .and. bkgnd_profile
   end function bkgnd_henyey_conflicts_profile

   pure function vmix_resolve_kd_min(kd_min, kt_bg) result(kd_min_eff)
      !! Resolve the `bkgnd_kd_min` "unset" sentinel to the reference
      !! default `HENYEY_KD_MIN_FRAC * kt_bg` (MOM6 `KD_MIN`, default
      !! `0.01*KD`).  A NEGATIVE `kd_min` means "not set by the user";
      !! zero and positive values are taken literally, so `kd_min = 0`
      !! is a legal way to ask for no floor at all.
      !!
      !! Called from `vmix_assemble` on every entry rather than resolved
      !! once at configure, so a directly-constructed slot (every unit test
      !! that never goes through `configure_ocean_vmix`) gets exactly the
      !! same default as a production run.  It is a two-flop scalar.
      real(wp), intent(in) :: kd_min
         !! The raw `bkgnd_kd_min` knob; negative = unset sentinel.
      real(wp), intent(in) :: kt_bg
         !! Scalar background tracer diffusivity the default is a fraction of.
      real(wp) :: kd_min_eff

      kd_min_eff = merge(HENYEY_KD_MIN_FRAC*kt_bg, kd_min, kd_min < 0.0_wp)
   end function vmix_resolve_kd_min

   pure subroutine vmix_assemble_clip_henyey_impl(nx, ny, nzp1, kv, kt, ks, &
                                                  kv_bg, kt_bg, ks_bg, kv_max, kd_max, &
                                                  n0_2omega, henyey_max_lat, kd_min, &
                                                  geolat)
      !! `vmix_assemble_clip_impl` with the Henyey latitude factor
      !! (`henyey_lat_factor_impl`) scaling the two SCALAR TRACER background
      !! floors, each then floored at the minimum diffusivity:
      !!
      !!   kt floor = max(kd_min, kt_bg · L(φ))
      !!   ks floor = max(kd_min, ks_bg · L(φ))
      !!
      !! matching the reference `Kd_sfc = max(Kd_min, Kd · L)`.  Without the
      !! `kd_min` floor the background would collapse toward zero at the
      !! equator (`L(0°) = 0` exactly) and at the poleward clamp.
      !!
      !! The momentum floor stays the plain scalar `kv_bg` — see the
      !! `bkgnd_henyey` docstring for why (rdb's `kv_bg` is an
      !! independent momentum floor, not `prandtl · kt_bg`).
      !!
      !! Bryan-Lewis composition is impossible here by construction: this
      !! routine takes no depth profile, and `vmix_assemble` reaches it only
      !! on the `bkgnd_profile == .false.` arm.
      !!
      !! Interior interfaces k = 2..nzp1-1.  The latitude factor is a column
      !! constant, so it is computed once per (i,j) and reused down the
      !! column.  No optional arguments: the Henyey parameter group is
      !! supplied in full or this routine is not the one you call.
      !!
      !! Explicit-shape args so the do concurrent stays descriptor-walk free.
      integer, intent(in) :: nx, ny, nzp1
      real(wp), intent(inout) :: kv(nx, ny, nzp1)
      real(wp), intent(inout) :: kt(nx, ny, nzp1)
      real(wp), intent(inout) :: ks(nx, ny, nzp1)
      real(wp), intent(in) :: kv_bg, kt_bg, ks_bg, kv_max, kd_max
      real(wp), intent(in) :: n0_2omega, henyey_max_lat
         !! `N0/(2Ω)` ratio (nondim) and the poleward cutoff latitude (degN).
      real(wp), intent(in) :: kd_min
         !! Minimum background tracer diffusivity (m^2/s), ALREADY resolved
         !! through `vmix_resolve_kd_min` — this kernel never sees the
         !! negative sentinel.
      real(wp), intent(in) :: geolat(nx, ny)
         !! T-point latitude (degN) — `ocean_metrics_t%geolatT`, the T
         !! stagger.  The corner field `geolatBu` is `(nx+1, ny+1)` and is
         !! NOT a substitute.

      integer :: i, j, k
      real(wp) :: kt_floor, ks_floor

      do concurrent(j=1:ny, i=1:nx) local(kt_floor, ks_floor, k)
         ! ONE write site per local (a `local()` variable reassigned across
         ! a branch is the recorded gfortran-15.1 corruption shape).
         kt_floor = max(kd_min, kt_bg* &
                        henyey_lat_factor_impl(geolat(i, j), n0_2omega, henyey_max_lat))
         ks_floor = max(kd_min, ks_bg* &
                        henyey_lat_factor_impl(geolat(i, j), n0_2omega, henyey_max_lat))
         do k = 2, nzp1 - 1
            kv(i, j, k) = min(max(kv(i, j, k), kv_bg), kv_max)
            kt(i, j, k) = min(max(kt(i, j, k), kt_floor), kd_max)
            ks(i, j, k) = min(max(ks(i, j, k), ks_floor), kd_max)
         end do
      end do
   end subroutine vmix_assemble_clip_henyey_impl

   pure subroutine vmix_assemble_clip_profile_impl(nx, ny, nzp1, kv, kt, ks, &
                                                   kd_bg, prandtl, kv_max, kd_max)
      !! C7 floor + ceiling: same as `vmix_assemble_clip_impl` but the
      !! tracer floor is the per-interface Bryan-Lewis `kd_bg` field and the
      !! momentum floor is `prandtl·kd_bg` (MOM6 background Prandtl tie).
      !! Interior interfaces k = 2..nzp1-1.
      integer, intent(in) :: nx, ny, nzp1
      real(wp), intent(inout) :: kv(nx, ny, nzp1)
      real(wp), intent(inout) :: kt(nx, ny, nzp1)
      real(wp), intent(inout) :: ks(nx, ny, nzp1)
      real(wp), intent(in) :: kd_bg(nx, ny, nzp1)
      real(wp), intent(in) :: prandtl, kv_max, kd_max

      integer :: i, j, k

      do concurrent(k=2:nzp1 - 1, j=1:ny, i=1:nx)
         kv(i, j, k) = min(max(kv(i, j, k), prandtl*kd_bg(i, j, k)), kv_max)
         kt(i, j, k) = min(max(kt(i, j, k), kd_bg(i, j, k)), kd_max)
         ks(i, j, k) = min(max(ks(i, j, k), kd_bg(i, j, k)), kd_max)
      end do
   end subroutine vmix_assemble_clip_profile_impl

   pure subroutine vmix_smooth_121_impl(nx, ny, nzp1, fld, scratch, wet_mask)
      !! One in-plane 1-2-1 horizontal smoothing pass on interior
      !! interfaces k = 2..nzp1-1.  Wet-mask aware: contributions from
      !! dry neighbours (wet_mask == 0) are excluded and the 9-point
      !! stencil weight is renormalised over the wet cells only.  A dry
      !! centre column (wet_mask(i,j) == 0) is left unchanged — no
      !! leakage into or out of dry cells.
      !!
      !! Reads `fld` into `scratch`, then writes the smoothed result back
      !! to `fld`.  Interior (i,j) only (i = 2..nx-1, j = 2..ny-1) —
      !! edge columns are left unchanged so closed-wall ghost rows don't
      !! bleed into the smoothed field.
      integer, intent(in) :: nx, ny, nzp1
      real(wp), intent(inout) :: fld(nx, ny, nzp1)
      real(wp), intent(inout) :: scratch(nx, ny, nzp1)
      real(wp), intent(in) :: wet_mask(nx, ny)
         !! 1 = wet, 0 = dry.  Shape (nx, ny).

      integer :: i, j, k
      real(wp) :: w_c, w_e, w_n, w_w, w_s_pt, w_ne, w_nw, w_se, w_sw, wsum

      do concurrent(k=1:nzp1, j=1:ny, i=1:nx)
         scratch(i, j, k) = fld(i, j, k)
      end do
      do concurrent(k=2:nzp1 - 1, j=2:ny - 1, i=2:nx - 1) &
         local(w_c, w_e, w_n, w_w, w_s_pt, w_ne, w_nw, w_se, w_sw, wsum)
         ! Standard 9-point [[1,2,1],[2,4,2],[1,2,1]]/16 weights, masked.
         w_c = 4.0_wp*wet_mask(i, j)
         w_e = 2.0_wp*wet_mask(i + 1, j)
         w_w = 2.0_wp*wet_mask(i - 1, j)
         w_n = 2.0_wp*wet_mask(i, j + 1)
         w_s_pt = 2.0_wp*wet_mask(i, j - 1)
         w_ne = wet_mask(i + 1, j + 1)
         w_nw = wet_mask(i - 1, j + 1)
         w_se = wet_mask(i + 1, j - 1)
         w_sw = wet_mask(i - 1, j - 1)
         wsum = w_c + w_e + w_w + w_n + w_s_pt + w_ne + w_nw + w_se + w_sw
         if (wsum <= 0.0_wp .or. wet_mask(i, j) == 0.0_wp) then
            ! Dry centre or fully isolated cell: leave unchanged.
            fld(i, j, k) = scratch(i, j, k)
         else
            fld(i, j, k) = (w_c*scratch(i, j, k) &
                            + w_e*scratch(i + 1, j, k) + w_w*scratch(i - 1, j, k) &
                            + w_n*scratch(i, j + 1, k) + w_s_pt*scratch(i, j - 1, k) &
                            + w_ne*scratch(i + 1, j + 1, k) + w_nw*scratch(i - 1, j + 1, k) &
                            + w_se*scratch(i + 1, j - 1, k) + w_sw*scratch(i - 1, j - 1, k)) &
                           /wsum
         end if
      end do
   end subroutine vmix_smooth_121_impl

   pure subroutine vmix_guard_impl(nx, ny, nzp1, kv, kt, ks, bad_count)
      !! Count negative or NaN diffusivities across interior interfaces.
      !! A reduction over fresh scratch — uses `!$acc parallel loop
      !! reduction` (the project rule for device reductions; bare
      !! `count`/`sum` over device data can silently return 0).
      integer, intent(in) :: nx, ny, nzp1
      real(wp), intent(in) :: kv(nx, ny, nzp1)
      real(wp), intent(in) :: kt(nx, ny, nzp1)
      real(wp), intent(in) :: ks(nx, ny, nzp1)
      integer, intent(out) :: bad_count

      integer :: i, j, k

      bad_count = 0
      do concurrent(k=2:nzp1 - 1, j=1:ny, i=1:nx) reduce(+:bad_count)
         if (kv(i, j, k) < 0.0_wp .or. kv(i, j, k) /= kv(i, j, k)) bad_count = bad_count + 1
         if (kt(i, j, k) < 0.0_wp .or. kt(i, j, k) /= kt(i, j, k)) bad_count = bad_count + 1
         if (ks(i, j, k) < 0.0_wp .or. ks(i, j, k) /= ks(i, j, k)) bad_count = bad_count + 1
      end do
   end subroutine vmix_guard_impl

   pure function ocean_vmix_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the vertical mixing slot
      !! (0 when unallocated). One arr_bytes term per array — add a
      !! term here when a new allocatable joins the type.
      class(ocean_vmix_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = arr_bytes(this%bl_depth) &
               + arr_bytes(this%b0) &
               + arr_bytes(this%kv) &
               + arr_bytes(this%kt) &
               + arr_bytes(this%ks) &
               + arr_bytes(this%gamma_t) &
               + arr_bytes(this%gamma_s) &
               + arr_bytes(this%smooth_scratch) &
               + arr_bytes(this%kd_bg)
   end function ocean_vmix_bytes

end module rdb_ocean_vmix
