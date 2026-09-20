# Roundabout — Closure & Scheme Matrix

**Single source of truth for "which closure / scheme is enabled in which
regime."** When you add, port, or remove a closure, **update the relevant cell
in the same PR** — `tools/check_closure_matrix.py` (pre-commit/CI) fails the
build if a namelist closure knob has no row here, or a listed test is missing.

Cell convention: a cell holds the **namelist knob** that enables the scheme in
that regime (so the matrix doubles as "how to turn it on"), or `default` (active
unless another is chosen), or `—` (not available in that regime).

This build ships the **ocean** regime only (`sim_type='ocean'`, Arakawa C-grid);
the coastal A-grid / unstructured path was split out into its own repository.
Ocean knobs live in `&ocean_<group>_nml` (the `ocean_` prefix is dropped from
each key); the shared vertical/tracer config lives in `&nonhydrostatic_nml` and
`&vcoord_nml`.

Detailed prose + limitations: [`CAPABILITIES_AND_LIMITATIONS.md`](CAPABILITIES_AND_LIMITATIONS.md).
Physics/namelist reference: [`REFERENCE.md`](REFERENCE.md). Acronyms:
[`codebase/CONCEPTS.md`](codebase/CONCEPTS.md). When this file and the code
disagree, **the code is authority and this file is a bug** — fix it.

Last reconciled with the code: **2026-09-14** (outer time-split section added — `&ocean_bt_nml split_scheme` had no row at all, and the fast-loop Coriolis-reference seam it selects had none either; earlier: 2026-08-23 (coastal-regime carve-out — the coastal-S / coastal-U columns, the coastal-only closure rows and the `&physics_nml` / `&nonhydrostatic_nml` coastal knob tables were removed with the coastal path; earlier: 2026-08-04 vertical-mixing section; background-mixing row added — Bryan-Lewis had shipped with no matrix row, Henyey lands in the same PR that fixes the gap; row rewritten when Henyey was brought to MOM6 parity: Bryan-Lewis XOR Henyey instead of Henyey-requires-Bryan-Lewis, `bkgnd_kd_min` floor added)).

---

## Vertical mixing (3D / multilayer only — no-op at nz=1 and in 2D)

| Closure | ocean | Numerics | Primary test |
|---|---|---|---|
| PP81 (Pacanowski–Philander Ri) | `default` | implicit (backward-Euler tridiagonal) | `test_ocean_pp81` |
| KPP (Large et al. 1994) | `use_kpp` (Phase 1) | overlay on PP81; surface buoyancy flux `B_0 = (g/ρ₀)·(α_T·F_T − β_S·F_S)` with `F_T = Q_heat/(ρ₀·c_p)`, `F_S = Q_salt/ρ₀` — the `1/ρ₀` is load-bearing because `α_T`/`β_S` are the DIMENSIONAL linear-EOS sensitivities (kg m⁻³ per °C / PSU), and it is the same number EPBL forms as `g·ρ₀·(dSV/dT·F_T + dSV/dS·F_S)`; persisted per column as the diagnostic `vmix%b0` (mirror of `epbl%b0`); `B_0` shortwave method `&ocean_thermo_nml kpp_sw_method` (all_sw/mxl_sw/lv1_sw, PR-21) | `test_ocean_kpp*`, `test_ocean_buoyancy_flux`, `test_ocean_sw_bl_coupling` |
| EPBL (Reichl & Hallberg 2018 energetics) | `&ocean_epbl_nml enable` | implicit; closed-form energy solve; replaces the KPP overlay, combines with PP81 (`add`/`max`); penetrating-SW TKE ledger `&ocean_thermo_nml epbl_sw_ctke` (PR-21) | `test_ocean_epbl`, `test_ocean_sw_bl_coupling` |
| kappa-shear (Jackson, Hallberg & Legg 2008 prognostic shear) | `&ocean_kappa_shear_nml enable` | per-column coupled (κ, TKE) Picard solve + adaptive substepping; INTERIOR closure — additive with KPP/EPBL/PP81, thermo-cadence compute, every-stage merge; opt-in VERTEX form (`at_vertex`, MOM6 VERTEX_SHEAR / OM5) solves at C-grid corners from native face velocities and averages corner Kd back to tracer points (arithmetic or `vertex_geometric_mean` + `vertex_geomean_kdmin` floor); Kv routed corner→face (`prandtl_turb·kd_corner` via `vdiff_apply_momentum kv_corner_source`, MOM6 `Kv_shear_Bu` — cell-centred kv merge suppressed in vertex mode) | `test_ocean_kappa_shear` |
| tidal mixing (St-Laurent/Simmons internal-tide, Jayne & St Laurent 2001 / St Laurent et al. 2002 / Simmons et al. 2004) | `&ocean_tidal_mixing_nml enable` | bottom-intensified `Kd = q·μ·E·F(z)/(ρ·(N²+Ω²))`, exp decay from bed (scale ζ), conservative flux bookkeeping; INTERIOR closure — additive with KPP/EPBL/PP81/kappa-shear, thermo-cadence compute, every-stage merge; E prescribed (v1) or `E=½ρ₀κ⟨h²⟩U²·N_bot` (Jayne & St Laurent 2001, `e_compute`) | `test_ocean_tidal_mixing` |
| convective adjustment (Brunt-Väisälä trigger, CVMix `CVMix_convection` / Cox 1984 / Marotzke 1991) | `&ocean_conv_nml enable` | `N² < n2_thresh` (interior interface, dense-over-light) ⇒ `kt = max(kt, kd_conv)`, `kv = max(kv, prandtl_conv·kd_conv)`; INTERIOR closure applied BELOW the active KPP/EPBL boundary layer; CONTRIBUTOR (max floor), runs every stage before `vmix_assemble`; `kd_conv` default 1.0 m²/s (~100× PP81's own Ri<0 ceiling) admissible only because vdiff is backward-Euler; writes `kv`/`kt` only, never `ks` (see CAPABILITIES_AND_LIMITATIONS.md) | `test_ocean_convection` |
| background mixing (Bryan & Lewis 1979 depth profile **XOR** Henyey, Wright & Flatte 1986 latitude factor) | `&ocean_vmix_nml bkgnd_profile` **xor** `bkgnd_henyey` | Two MUTUALLY EXCLUSIVE replacements for the SCALAR `kt_bg`/`ks_bg`/`kv_bg` assembly floor (enabling both fails loud at configure, matching MOM6's one-background-scheme rule). **`bkgnd_profile`**: per-interface `kd_bg(z) = kd_sfc + (kd_deep−kd_sfc)·[½+atan((\|z\|−z0)/Δ)/π]` recomputed each stage from the live column thickness (correct under any vcoord); `kv` floor = `bkgnd_prandtl·kd_bg`. **`bkgnd_henyey`** (Harrison & Hallberg 2008 constant-`N0` simplification; requires a non-cartesian `grid_config`, fail-loud at configure — `geolatT ≡ 0` on cartesian so every column would take the equatorial `L(0°)=0`): scales the SCALAR TRACER floors by a latitude-only factor `L(φ)` and floors the result, `max(bkgnd_kd_min, kt_bg·L(φ))` (MOM6 `KD_MIN`; negative `bkgnd_kd_min` ⇒ `0.01·kt_bg`). Knobs `bkgnd_henyey_n0_2omega`/`bkgnd_henyey_max_lat`, defaults 20/95°, clamp compared against `\|φ\|` so both hemispheres; equator singularity floored via a fixed `1e-10` `\|sin φ\|` guard inside the `acosh` ratio only, so `L→0` smoothly and the `Kd_min` floor is what the background lands on. `kv_bg` is NOT latitude-scaled (rdb's momentum floor is independent of `kt_bg`, not `prandtl·Kd`). FLOOR contributor into `vmix_assemble`'s (b) stage, replacing (not combining with) the scalar background when on. Both default off ⇒ bit-identical | `test_ocean_bkgnd_mixing` |

Notes: ocean vmix is **implicit** (backward-Euler tridiagonal). The
two-equation k-ε closure and the Galperin/Canuto stability functions were
coastal-path closures and left with that split; the GLS-generic family (k-ω /
MY2.5) is **not** built.

**Ocean diffusivity assembly (`vmix_assemble`, C11).** On the ocean path every
interior / overlay closure (PP81, KPP, EPBL, kappa-shear, tidal mixing,
KV_ML_INVZ2, convective adjustment) **contributes** into `kv`/`kt`. `ks` (the salt diffusivity) is not
a contributor target — it is **derived** from `kt` by `vmix_split_kd_heat_salt`
(MOM6 `Kd_salt = Kd_int + Kd_extra_S`, PR-20), which runs after the last
contributor and immediately before `vmix_assemble`. With no double-diffusion
contributor shipped yet, `Kd_extra_S ≡ 0` so the split reduces to `ks := kt`
⇒ bit-identical. `vmix_assemble` is the single downstream gate run once per
stage before vdiff. It applies, in order: background floors
(`kv_bg`/`kt_bg`/`ks_bg`, defaulting to `pp81_*_bg` so the floor is a no-op for
the shipped path), ceilings (`kv_max`/`kd_max`, default `huge` = off), optional
1-2-1 horizontal smoothing of kv/kt/ks (`kd_smooth_iterations`, default 0), and
an optional negative/NaN guard (`vmix_guard`, default off). All defaults are
bit-identical to the pre-assembly chain. `vdiff_apply_tracers` takes
`kt_source` (temperature) and `ks_source` (salinity + every passive tracer,
MOM6's `Kd_salt` convention) — `ks ≡ kt` until a double-diffusion contributor
lands (PR-33). Future Area-C contributors (tidal mixing, double-diffusion,
geothermal-adjacent floors, background profiles) plug in upstream of the
split. Tests: `test_ocean_vmix_assembly`, `test_ocean_vdiff`.

## Horizontal viscosity

| Closure | ocean | Test |
|---|---|---|
| Smagorinsky Laplacian | `smag` (Smag_KH) | — (exercised by ML tests) |
| Biharmonic ∇⁴ (constant ν₄) | `nu_4` | `test_ocean_hvisc` |
| Smagorinsky_AH (flow-aware biharmonic) | `smag_ah` | `test_ocean_smag_ah` |
| Leith (vorticity-gradient ν_h) | `lateral_closure="leith"` | `test_ocean_leith` |
| Leith-biharmonic (∇²ζ-scaled ν₄, MOM6 `LEITH_AH`) | `lateral_closure="leith_biharm"` (`c_leith_bi`) | `test_ocean_hvisc_leith_biharm` |
| Resolution-scaled viscosity (Hallberg 2013) | `&ocean_hvisc_nml resoln_scaled_visc` (needs `&ocean_varmix_nml enable`) | `test_ocean_hvisc_resoln` |
| Constant ν_h floor | `nu_h` | — |
| Live velocity-scale ν_h (Kh = U·dx·\|u\|) | `kh_vel_scale_live` | `test_ocean_hvisc_aniso` |
| Anisotropic ν_h (Smith & McWilliams 2003) | `kh_aniso` + `aniso_dir` (stress-tensor path) | `test_ocean_hvisc_aniso` |
| Tensor-strain Laplacian (free-slip walls) | `form="tensor"` | — |
| Fox-Kemper ML-eddy restratification (B5) | `&ocean_foxkemper_nml enable` | `test_ocean_foxkemper` |
| Gent-McWilliams thickness diffusion ([2]) | `&ocean_gm_nml enable` (needs `&ocean_slopes_nml enable`) | `test_ocean_gm` |
| Redi neutral (along-isopycnal) tracer diffusion ([3]) | `&ocean_redi_nml enable` (continuous variant) | `test_ocean_redi` |
| MEKE prognostic eddy energy ([5]) | `&ocean_meke_nml enable` (needs `&ocean_gm_nml enable`) | `test_ocean_meke` |
| MEKE harmonic backscatter ([5], negative-ν momentum return) | `&ocean_meke_nml backscatter` + `backscatter_visc_coeff_ku` (needs a flow-aware closure active) | `test_ocean_meke_backscatter` |

Notes: explicit biharmonic is CFL-capped (`ν₄·dt·((π/dx)²+(π/dy)²)² ≤ 2`); MOM6
production reaches higher ν₄ via implicit time-stepping we don't have. The ocean
velocity-Laplacian paths (scalar `nu_h` and the flow-aware per-face closures) take
an optional per-face harmonic clamp `&ocean_hvisc_nml bound_kh` (MOM6 `BOUND_KH`,
ceiling `bound_coef·0.125/(dt·(1/dx²+1/dy²))`, `test_ocean_hvisc_kh_bound`) —
load-bearing on the split-explicit ocean path, where an over-large `ν_h·dt/dx²`
anti-damps grid-scale barotropic gravity modes through the frozen `F_bt` forcing
(see `src/core/ocean/README.md`); the Lagrangian double-gyre configs run it with
`bound_coef = 0.15`. The ocean tensor form is free-slip. The ocean `lateral_closure` knob is **fail-loud** — a
tag with no dispatcher kernel (or a mistyped/garbage string) aborts at configure
(`validate_config` → `lateral_closure_is_implemented`) instead of silently
falling back to background-only viscosity. `leith_biharm` fills the per-face
`nu4_face_*` (same biharmonic apply as `smag_ah`); the harmonic Laplacian stays
on the scalar `nu_h` underneath it. The biharmonic add-on (`nu_4` / `smag_ah`
/ `leith_biharm`) composes with **all three** harmonic dispatch arms,
including the `stress_tensor` path — an early `return` used to silently
disable it whenever `stress_tensor=.true.`, which is now fixed; `kh_aniso`
(only consulted on the `stress_tensor` path) is likewise no longer mutually
exclusive with the biharmonic — both are independent linear operators that
superpose.

Fox-Kemper (B5, `rdb_ocean_mle`) is a lateral **restratification** closure, not
a viscosity: submesoscale mixed-layer eddies slump lateral buoyancy fronts via
an overturning streamfunction `Ψ = Ce·(H_ml²/|f|)·∇b̄·μ(z)`.  It injects
ML-confined per-layer mass transports (`uhml`/`vhml`) into the continuity mass
fluxes BEFORE the divergence — never touching velocities — so it is conservative
by construction (`Σ_k a(k)=0`, a closed overturning cell).  `H_ml` is taken from
`epbl%mld` (EPBL is the enabler; `enable` requires `ocean_epbl_nml enable`).  Runs
at THERMO cadence, once per outer step.  Two timescale forms (`use_mom_mixrate`):
the bare `Ce/max(|f|,f_floor)` FK08 floor form (analytic-gate default), and the
FK11 momentum-mixrate form (`use_mom_mixrate=.true.`, **production-recommended** —
it suppresses restratification under vigorous mixing).  `resolution_taper` is a
B2 hook (hard config error until B2 lands); slow-filtered-MLD second transport is
deferred (instantaneous MLD only).  Default off ⇒ bit-identical.

Gent-McWilliams ([2], `rdb_ocean_gm`) is the interior eddy-induced **bolus**
thickness diffusion (GM90/Griffies98): the skew-flux streamfunction
`Ψ = -KhTh·dy·S` on the stored isopycnal slope `S` (`&ocean_slopes_nml`, a loud
configure prerequisite) is turned into per-layer thickness transports
(`uhD`/`vhD`) by the MOM6 `uhtot` column recurrence with a safe-streamfunction
slope limiter + a mass-availability limiter (keeps `h ≥ H_VANISHED` without a
post-hoc clamp).  Folded into the continuity mass fluxes BEFORE the divergence
like Fox-Kemper — conservative by construction (`Σ_k uhD = 0`).  KhTh is a 2D
face field (constant-fill for v1; CFL-clamped via `khth_max_cfl`); the VarMix /
MEKE seam will make it spatially varying (`+=`).  `gm_src` carries the
`-¼·Σ_k ρ₀·KH·S²·N²·h` PE release for the future MEKE coupling.  Runs at THERMO
cadence (owns the slope refresh in the split driver).  Deferred: bottom-blocking,
FGNV/EBT/int_slope.  Default off ⇒ bit-identical.
**MOM6 divergence (recorded, not deferred — PR-8):** MOM6's
thickness-diffusion (ALE mode) always runs its top layer through a
linear return-flow closure (`nk_linear = max(GV%nkml,1) = 1`, not a namelist
parameter).  Roundabout's `gm_column_x`/`_y` never did this — the field that would
have selected it (`nk_linear`) was permanently 0, and PR-8 found the branch it
gated was wired to the wrong end of the column (bed-most, not surface-most, per
Roundabout's bottom-up convention) — so it was deleted as dead rather than wired.
Wiring a correct surface region is a live physics change (moves every GM
answer) left to a future PR + analytical test; see the `!! DIVERGENCE
(MOM6 nk_linear):` docstring on `gm_column_x`.

MEKE ([5], `rdb_ocean_meke`) closes the GM↔eddy-energy loop: a 2D prognostic
eddy-kinetic-energy field `E(i,j)` (Jansen 2015 / Eden-Greatbatch 2008 / Marshall
2012) sourced by GM's `gm_src` PE release (`&ocean_gm_nml enable` is a loud
prerequisite — MEKE needs it), damped by an implicit backward-Euler bottom drag
(`drag_rate = ρ₀·i_mass·√(cdrag²·(2·γb²·E + u_bbl² + uscale²))`, the ρ₀ factor
being MOM6's `GV%H_to_RZ`; fixed in PR-5 — the drag rate was previously missing
this factor and was ~ρ₀ too weak, a dimensional bug not a knob),
and transported by a harmonic-mass Laplacian (+ optional biharmonic).  Strang
split per thermo step: explicit source bump → drag half → diffusion → drag half
(the two half-drags collapse to one full drag when neither `meke_kh≥0` nor
`meke_k4≥0`).  The derived diffusivity `kh = khcoeff·√(2·γt²·E)·Lmix` (harmonic
sum of deformation / frictional / Rhines / Eady / grid scales, each
alpha-gated) is added as the geometric mean `khth_fac·√(kh_i·kh_{i+1})` into the
VarMix face KhTh/KhTr accumulator BEFORE GM's CFL clamp (so GM consumes the
MEKE-augmented KhTh next; one-step gm_src lag).  `meke` is restart-persistent.
`khth_fac=khtr_fac=0` (default) ⇒ feedback inert; VarMix off ⇒ E still evolves
but the feedback has no face accumulator.  The Rhines scale is live
(`alpha_rhines>0`; `β=|∇f|` from the MEKE slot's own `f_centre`, filled at setup
from the Coriolis path and scaled by `idxT`/`idyT`), and upwind barotropic
advection is live (`advection_factor>0`; mass-weighted `baroHu` from
`mass_flux_*_layer`, conservative).  Remaining upstream gaps: per-face BBL
`drag_visc` (needs `Kv_bbl`/`bbl_thick` on the bottom-drag slot; drag carried by
`cdrag`) and the frictional `mom_src` source (needs a per-cell lateral-dissipation
output on the hvisc slot).  Default off ⇒ bit-identical.

MEKE harmonic backscatter (`&ocean_meke_nml backscatter`) closes the
eddy-energy→momentum loop in the other direction: `meke_step` fills
`Ku = backscatter_visc_coeff_ku·√(2·E)·Lmix` (MOM6 `MEKE_VISCOSITY_COEFF_KU`;
plain `√(2·E)` — no `γt²`, unlike kh; harmonic only in v1, vertical structure
`BS_struct=1`), and `meke_backscatter_apply`
subtracts a face-average of `Ku` from the resolved per-face harmonic viscosity
(`lateral_mix%ah_face_x/y`) so the NET coefficient can go negative — a
negative-viscosity energy return into the resolved flow.  The net is floored at
the forward-Euler viscous-CFL lower bound `−0.8·0.5/(dt·(idx²+idy²))` per face
(MOM6 `BACKSCATTER_UNDERBOUND`), which bounds the negative mode's growth rate but
does NOT stabilise it alone — a positive biharmonic backstop (`nu_4`/`smag_ah`/
`leith_biharm`) is mandatory, enforced fail-loud at configure.
Runs after `ocean_lateral_mix_compute` (reads the prior thermo step's `Ku`), needs
a flow-aware closure active for `ah_face_*` to be consumed.  Default off ⇒
bit-identical.  Deferred: biharmonic `Au`, EBT/SQG vertical structure.

## Bottom drag

| Form | ocean | Notes |
|---|---|---|
| Quadratic log-layer (Cd) | `default` | MOM6/ROMS `Cd ≈ 2.5e-3` |
| Linear (Rayleigh) | `form="linear"` | Rayleigh drag rate `r` (1/s) |
| HBBL-distributed | `hbbl` | spreads stress over bottom `hbbl` m |
| Implicit-fold (backward-Euler vdiff bed diagonal) | `&ocean_vdiff_nml implicit_drag` | folds bottom drag into the vdiff bed (`k=1`) diagonal as a stress bottom-BC instead of the explicit pre-solve add; thin-layer (z*/ZSTAR_FULL pinch-out) CFL-stable. Mutually exclusive with `&ocean_bdrag_nml implicit` (split-apply) and HBBL (`hbbl>0`), fail-loud at configure. Test: `test_ocean_vdiff_implicit_stress_drag` |

## Porous barriers (subgrid topography; ocean only)

Represent a subgrid sill/strait by reducing the **open** fraction of a C-grid
face below its full width, per layer, so a deep sill blocks the bottom layers
while the surface layers stay fully open.

| Scheme | ocean | Numerics | Primary test |
|---|---|---|---|
| Porous barriers (Adcroft 2013 three-parameter fit) | `&ocean_porous_nml enable` | per-face `d_min`/`d_max`/`d_avg` along-face topographic heights ⇒ monotone open-width profile `w(η)`; the layer-averaged OPEN-AREA fraction is the exact difference of the profile's vertical integral over the layer, `min(1, (A(η_hi)−A(η_lo))/(η_hi−η_lo))`. A DEGENERATE face (`d_max ≤ d_min`) and a VANISHED layer (`dz ≤ H_VANISHED`) are special-cased — see the caveat below. Recomputed on the device ONCE PER OUTER STEP and held across both RK2 stages (MOM6's cadence), then MULTIPLIED into the per-layer transports of continuity-PPM (both the fused and the direction-split forms, both directions, before the barotropic renormalisation, which carries the open fraction in BOTH its `sum_h` denominator and its per-layer increment) and of the Coriolis/advection **transport** forms (`sadourny_energy`, `sadourny_hk`). The default velocity-form `sadourny` has no `u·h·dy_cu` transport to narrow, so porosity reaches it only through continuity. The BAROTROPIC substep transports on `dy_cu_bt`/`dx_cv_bt` — the same widths scaled by the COLUMN-INTEGRATED open fraction, which is identically the thickness-weighted mean of the per-layer fractions (an exact telescoping identity, both being the same integral of `w` over the column). Without that the barotropic solve would be porous-blind and the per-layer renormalisation to `uhbt` would hand the blocked transport straight back, leaving the barrier a vertical redistribution that never reduces net flow. **This is NOT `BT_cont` parity**: MOM6's production barotropic face area is `Σ_k (dy_Cu·por_k)·h_marginal_k·visc_rem_k` (PPM marginal thickness AND `visc_rem`, neither of which appears here); its `Σ_k h_k·(dy_Cu·por_k)` form is the open-boundary-segment branch only, and `set_local_BT_cont_types` carries no `por` at all. NOT applied to: `areaCu`/`areaCv` (so BT Coriolis + KE use un-narrowed areas against narrowed transports, as in MOM6), GM/Redi/MEKE/hdiff, or sea ice. **Fail-loud exclusions** (`validate_config`): `&ocean_bt_nml bt_halo > 0` (the wide-halo BT clone re-fills its own `metrics_w` from the grid formula and nothing gives it the porous statistics, so the wide fast loop would silently transport on UN-narrowed widths) and `&ocean_wetdry_nml enable` (the vanishing-column interaction with the wet/dry outflow limiter is unvalidated). Single-rank. Default off ⇒ no kernel launch at all ⇒ byte-identical | `test_ocean_porous` |

**Subgrid-data caveat.** The fit needs min/max/mean of a bathymetry finer than
the model grid; MOM6 reads that from an offline `topog_edge.nc`. Roundabout has no
such file plumbing yet, so `source="resolved"` (the only implemented source)
samples the **resolved** bathymetry at three along-face points (the two face
corners and the midpoint) — a documented PROXY that captures along-face slope
but is blind to genuine subgrid structure. `source="file"` fails loud pending
the same file-forcing backend `&ocean_bt_nml wave_drag_form="file"` waits on.

Three properties of the proxy are load-bearing and easy to get wrong:

* **Uniform ALONG the face ⇒ inert, not a wall.** All three samples are
  cell-centre averages, so a ridge or shelf break running *parallel* to the
  face collapses them onto one value. The bare fit's limit there is a STEP at
  the two-cell mean height — a hard wall on every layer below it, on a face the
  grid resolves as open. The kernel therefore treats a degenerate statistic
  (`d_max ≤ d_min`) as FULLY OPEN and blocks nothing. Flat bathymetry is the
  same degenerate case, which is what makes the resolved source a literal no-op
  on a flat basin. Faces that *cross* the structure still narrow normally.
* **Land cells are gated out.** A corner sample averages four cells; a land
  elevation in that stencil would pull `d_max` up and block a wet–wet face the
  grid fully resolves (measured: ~37% spurious blockage from one 4000 m land
  diagonal beside a 4000 m column). A corner whose stencil is not entirely wet
  falls back to the two-cell face midpoint.
* **`eta_interp="max"` blocks the LEAST, not the most.** `w` is monotone
  increasing in the interface height, so the rule returning the *higher*
  (shallower) of the two adjacent interfaces leaves the most of the face open.
  `"min"` is the most blocking. Combined with the point above, `"max"` on
  smooth resolved bathymetry is close to inert.

A layer whose face thickness is at or below `H_VANISHED` (1.5e-4 m) gets
`por = 0`; MOM6 uses `Angstrom_Z` (1e-10 m), 1.5 million times smaller, so
layers between the two thresholds get a real fraction there and zero here, and
the column-integrated fraction is the thickness-weighted mean over the
NON-vanished layers only. A wholly vanished column blocks every layer AND
zeroes the barotropic width, so the two modes agree.

Single-rank only: the open-area fields have no halo exchange (MOM6
`pass_vector`s them), so a multi-rank C-grid run would see stale seam columns.

## Wind-stress application (multilayer)

| Mode | ocean | Notes |
|---|---|---|
| Surface-concentrated (top layer) | `default` (DIRECT_STRESS) | generates Ekman/ML shear; required for solver Kato-Phillips √t |
| Implicit-fold (backward-Euler vdiff surface row) | `&ocean_vdiff_nml implicit_stress` | folds wind stress into the vdiff surface (`k=nz`) RHS row as a Neumann top-BC instead of the explicit pre-solve add; thin-layer CFL-stable. Incompatible with `&ocean_vmix_nml direct_stress` (distributed stress), fail-loud at configure. Test: `test_ocean_vdiff_implicit_stress_drag` |

## Baroclinic pressure gradient / PGF

| Variant | ocean | Knob |
|---|---|---|
| Montgomery / FV-lite / FV-Wright / gprime / FV-MOM6 | `form=` | `&ocean_pgf_nml`; **default `"mont"`** |
| Montgomery potential (Boussinesq `M = p/ρ0 + (g·ρ/ρ0)·z`, vertical recursion seeded at the free surface + ONE horizontal difference, plus the `z_eff`-weighted horizontal-density term) | `form="mont"` (**default**) | `&ocean_pgf_nml`; general-purpose — valid over sloping bathymetry and every vcoord, EXACT at rest in isopycnal (`VCOORD_LAGRANGIAN`) columns wherever a layer is present on both sides of a face, and algebraically identical to `fv_lite` on aligned columns. Never forms a pressure stack, so it avoids `fv_lite`'s ~1e4 cancellation when differencing a ~2e7 Pa column. Where layer thicknesses are UNEQUAL across a face (σ/z*σ over a slope) `mont` and `fv_lite` are DIFFERENT discretisations of the same term — neither is exact there; measured on the 34 shipped namelists that do not pin `form` the end-of-run energy differs by ≤ 5e-5 relative, and on quiescent σ-seamount cases `mont` is 3-5 orders QUIETER. Tests: `test_ocean_pgf_mont`, `test_ocean_pgf_fv` |
| FV-MOM6 in-layer T/S reconstruction (PLM/PPM Boole density integral, BOTH the per-column vertical integral and the 5-point cross-face horizontal one) | `reconstruct_for_pressure=.true.` (+ `recon_scheme=1\|2`) | `&ocean_pgf_nml`; FV_MOM6 only, default off ⇒ PCM bit-identical. **EXACT AT REST**: with a linear EOS and T/S linear in z the face acceleration is round-off (~2e-15 m/s², `C·ε·g·H/dx`) for ANY layer geometry, tilted σ layers included — the defining property of the analytic-FV PGF (Adcroft, Hallberg & Harrison 2008; Yung et al. 2026 §2.4, Fig. 5b→5c). Two pieces are load-bearing and were both missing before 2026-09-20: the boundary layers (k=1, k=nz) take a LINEAR-EXACT one-sided edge pair instead of a PCM flatten (a flatten left the full σ truncation error in exactly the layers next to the tilted boundary), and `intx_dpa`/`inty_dpa` are a 5-point Boole quadrature over sub-columns at the INTERPOLATED interface height instead of the two-column trapezoid `½(dpa_L+dpa_R)` (the trapezoid left the curvature residual `g·(−dρ/dz)·Δe²/12` at every tilted interface — the σ "second-kind" pressure-gradient error, Haney 1991). Measured on a 48×6×15 σ column over a 226→709 m linear bed: PCM 2.59e-8, interior-only reconstruction 9.43e-9, both pieces 9.6e-16 m/s². Costs 25 EOS evaluations per face per layer. Tests: `test_ocean_pgf_sigma_rest`, `test_ocean_pgf_reconstruct` |
| Reference densities — the Boussinesq divisor `ρ₀` (`du/dt = −(1/ρ₀)∂p/∂x`, every variant) and the FV-MOM6 anomaly baseline `ρ_ref` (`pa(top) = ρ_ref·g·η`, layer anomaly `(ρ_k−ρ_ref)·g·h`; also `compute_pbce`'s `g·ρ_ref/ρ₀`) | not separately settable — both take `&ocean_ic_nml rho_0` | Kept as two members (the roles differ; interchanging them is a known MOM6 bug class) but sourced from the ONE configured ρ₀ via `eos%rho0` in `configure_ocean_pgf` — the same scalar the EOS, EPBL, kappa-shear, tidal mixing, wave speed, the `η_ib` surface-pressure seam, GM/MEKE/Redi/MLE and the isopycnal slopes take. Before that wiring the slot kept a hard-coded 1035 while the EOS followed the namelist, so `rho_0 /= 1035` ran the EOS and the pressure gradient on two different reference densities, silently. Default `rho_0 = 1035` ⇒ bit-identical. Test: `test_ocean_pgf_rho_ref` |
| Grounded-layer gate (zero the face PGF where the layer's z-extents do not overlap across the face) | `pgf_skip_nonoverlap=.true.` | `&ocean_isopycnal_nml`; **default ON**, applied under `VCOORD_LAGRANGIAN` only ⇒ every other vcoord bit-identical. Covers `mont`/`fv_lite`/`fv_wright`/`fv_mom6` (`gprime` N/A — it warns); on `fv_mom6` the `z_centre` buffer is allocated for this gate alone, so the flag is latched before the PGF slot is initialised. Without it a grounded isopycnal layer leaves `g·(ρ_layer−ρ̄_ambient)·∂z/∂x` of PGF **at rest** — `test_ocean_pgf_grounded` |

## Coriolis

| Scheme | ocean | Knob |
|---|---|---|
| Sadourny PV-flux (enstrophy, velocity form) | `form="sadourny"` | ocean default |
| Sadourny energy-conserving (transport form q·vh) | `form="sadourny_energy"` | MOM6 SADOURNY75_ENERGY; faithful |
| Sadourny + Hollingsworth-Källén guard | `form="sadourny_hk"` | Arakawa-Hsu PV stencil |
| BOUND_CORIOLIS velocity-form clamp (energy scheme) | `bound_coriolis` (default off; `form="sadourny_energy"` only, fail-loud) | MOM6 `BOUND_CORIOLIS`; clamps CAu/CAv into the `(f+ζ)·v` range before the KE-grad subtraction. **Inert on all-wet columns** (Roundabout's cell-mean corner-h ⇒ PV flux already convex-bounded; MOM6's `hArea_q` telescopes identically ⇒ MOM6's clamp is inert too); engages only at land-masked / floored corners. H200 rim A/B: NULL |
| PV corner-thickness construction (energy scheme) | `corner_h` (default `"cell_mean"`; `"mom6_area"` `form="sadourny_energy"` only, fail-loud) | MOM6's area-weighted PV form. **Algebraically identical to `cell_mean` above the `H_MIN_PV` floor** (same `Σarea·h`/`Σarea`); differs only in the vanishing-thickness guard (MOM6 `vol_neglect` = pure 1/0 armor, no cap). Round-off no-op for `h≳1e-12`; predicted inert for the dt=800 disease |

## Outer time-split scheme (`&ocean_bt_nml split_scheme`)

Selects the OUTER (baroclinic) integrator that wraps the barotropic fast
loop. Both schemes ship and both are under test; the stability suite runs an
`ssp_rk2` twin of every case whose namelist does not pin a scheme.

| Scheme | ocean | What it does | Primary test |
|---|---|---|---|
| `pred_corr` | **`default`** | MOM6 predictor-corrector: `pc_be`-off-centred predictor, slow tendencies on the `u_av`/`h_av` step time-means, ONE prognostic update in the corrector, forward-backward gravity-wave pairing. Neutrally stable to `ω·dt = 2`. `validate_config` refuses it FAIL-LOUD outside its v1 envelope (`eulerian_z`, wet/dry, `dt_tracer_advect_ratio > 1`) — those configurations must pin `ssp_rk2` | `test_ocean_pred_corr`, `test_ocean_cor_ref_seiche` |
| `ssp_rk2` | `split_scheme = "ssp_rk2"` — **EXPERIMENTAL** | Two identical stages + SSP average. Widest envelope: the only scheme that takes `eulerian_z`, `&ocean_wetdry_nml enable` and `dt_tracer_advect_ratio > 1`. **Amplifies internal gravity waves by `√(1+(ω·dt)⁴/4)` per step** — see the note below. Fully supported and fully tested; *experimental* labels the ANSWER, not the code path | `test_ocean_dyn_split`, `test_ocean_analytical` |
| Fast-loop Coriolis reference (`Cor_ref`, MOM6 `ubt_Cor`) | `default` (scheme-derived) | `subtract_fast_cor_ref` removes the barotropic Coriolis/advection already frozen into `F_bt` so the substep does not integrate it twice. The reference velocity `bt_work%cor_ref_u/v` is built by `set_cor_ref_velocity` from the SAME velocity the slow `cor%pv_flux_*` was evaluated on — the stage-entry `bt_ubt/bt_vbt` under `ssp_rk2`, the depth mean of `u_av/v_av` under `pred_corr`, weighted exactly as the forcing depth-mean was (h, or h·visc_rem under `forcing_visc_rem`). Mismatch them and the residual `f × (v̄_ref − v̄_slow)` forces every substep and pumps the basin's gravest Poincaré seiche | `test_ocean_bt_cor_ref`, `test_ocean_cor_ref_seiche` |

**Why `ssp_rk2` is labelled EXPERIMENTAL — the number.** On
`validation_examples/ocean/eady/resting_stratified_channel.nml` — flat bed,
periodic, stably stratified, at rest, ±0.5 mK seed, **no energy source** —
`ssp_rk2` manufactures **En = 2.992E-05 m²/s²** by day 25 (7.7 mm/s rms) and
is still climbing on a 2.5-day e-folding; `pred_corr` on the identical file
holds **1.739E-09** (17 000×, 83-day e-folding). The outer split is the
cause: the Coriolis form, the ALE remap and the PGF form were each
substituted and each moved the answer < 0.1 %; `dT/dz = 0` dropped En 119×;
removing the lateral viscosity RAISED it. It is a `(ω·dt)⁴` noise floor, so
a forced, energetic, viscous run sits decades above it and never notices,
while a quiescent or long-spin-up one does not. The suite carries it as a
scoped XFAIL on `resting_stratified_channel__ssp_rk2` — **do not close it by
raising `en_rest_max`**.  `ssp_rk2` stays supported and stays on the scheme
axis; the label tells you what choosing it costs, it does not deprecate it.

**`pred_corr` became the default on 2026-09-14.** Its GPU path is clean:
the seven failures that retargeting the dyn-core unit suite onto it used to
leave on the **NVHPC GPU build only** (`periodic`, `obc_baroclinic`,
`dyn_split`, `ice_restart`, `wetdry` non-finite within 1-5 steps, plus
`restart` bit-exactness and `p_surf` gauge invariance) had one cause —
`scratch_3d_buffer_t` attached uninitialised device memory and the pred_corr
predictor reads `du_visc`/`dv_visc` without recomputing them — and the suite
is 187/187 on both toolchains. It still does NOT eliminate the resting-state
growth, only slows it ~36× (83-day e-folding), which is why
`resting_stratified_channel` remains a scoped XFAIL on the settle gate under
the default too. Its per-step cost is ~7 % over `ssp_rk2`. Full detail in
[`CAPABILITIES_AND_LIMITATIONS.md`](CAPABILITIES_AND_LIMITATIONS.md) and the
`split_scheme` docstring in `rdb_config.F90`.

## Vertical coordinates (dispatch via `vcoord_type`)

| Coord | ocean |
|---|---|
| Sigma | `"sigma"` |
| Z-sigma hybrid | `"zsigma"` |
| Z-star (lite) | `"zstar"` |
| Z-star/sigma hybrid | `"zstar_sigma"` |
| Z-star full (per-column) | `"zstar_full"` |
| Eulerian-Z | `"eulerian_z"` |
| Z-fixed (gprime) | `"z_fixed"` / `"gprime"` |
| Isopycnal (rho) | `"rho"` (validation-grade; collapses weakly-stratified columns — HYCOM hybrid is the production follow-on) |
| Hybrid z*/isopycnal (HYCOM) | `"hycom"` (isopycnal interior + z* surface-resolution floor; the production isopycnal coord — reuses the `rho` inversion + `dsig` z* band) |

### ALE remap reconstruction (dispatch via `remap_method`)

The conservative integrate-and-redistribute is shared across all coords; only
the per-column reconstruction order changes. All methods conserve column
integrals to roundoff and are monotone (no new extrema).

| Method | `remap_method` | Order | Notes |
|---|---|---|---|
| Piecewise constant | `"pcm"` | 0th | donor cell, diffusive |
| Piecewise linear | `"plm"` | 1st | minmod-limited |
| Piecewise parabolic | `"ppm"` (`default`) | 2nd | Colella & Woodward 1984; uniform-grid `(7/12,−1/12)` edge estimate |
| PPM, non-uniform H4 edges | `"ppm_h4"` | 2nd (parabola) + 4th edges | thickness-weighted (White & Adcroft 2008) edge estimate — holds 4th order on non-uniform ALE layers where PPM degrades to 2nd; cuts spurious diapycnal mixing (Ilicak 2012). Reuses the PPM limiter + parabola + redistribute verbatim. Boundary edges: PCM-outermost + 3-cell H3 (cubic-fit upgrade deferred to PQM). |
| Piecewise quartic (PQM) | `"pqm"` | 4th–5th | White & Adcroft (2008) PQM_IH4IH3: implicit-h4 edge values + implicit-h3 edge slopes (per-column tridiagonal solves) → degree-4 reconstruction + W&A monotonicity limiter; prototype shows ~5th-order convergence (vs PPM ~2nd) on smooth profiles, ~10–20× lower remap error. Conservative + monotone. `N<5` falls back to PPM. ~98 regs (no spill); cadence-bounded. Opt-in (PPM stays default). |

#### Boundary-cell closure (`remap_boundary_extrap`) — orthogonal to the order above

Every reconstruction above PCM needs a stencil the outermost cells do not
have. By default — matching MOM6 `BOUNDARY_EXTRAPOLATION = False` — `k=1`
and `k=nz` collapse to PCM, so **the remap is first-order in the two cells
next to the bed and the surface whichever method is selected**: PLM, PPM,
PPM_H4 and PQM share the closure and remap a linear-in-z profile with the
same O(h) error there, while their interiors are already exact for it.

| Knob | Default | Effect |
|---|---|---|
| `&vcoord_nml remap_boundary_extrap` | `.false.` (bit-identical) | `.true.` ⇒ the boundary cells take the linear-exact one-sided edge pair `q ± dq_up·h_self/(h_self+h_nbr)` (`boundary_half_jump`, the remap-side twin of the FV PGF's `boundary_edges_linear`), making the whole column exact for a profile linear in z. Inert for `"pcm"` — there is no reconstruction to close. |

**Why it matters.** A stratified ocean at rest has a tracer profile linear
in z, and under a terrain-following coordinate the ALE remap runs on it
every thermo step. The first-order boundary closure therefore injects a
spurious diapycnal tracer flux into those two layers on every step; over a
slope it differs between neighbouring columns, which is a horizontal
density gradient, which is a pressure-gradient force — and with rotation it
feeds a growing grid mode trapped in exactly those layers. Measured on an
undamped 48×6×15 σ-over-slope rest case (bed 226 → 709 m, f-plane 75 °S,
exact FV pressure gradient, zero viscosity/drag/mixing):

| closure | `σ_En` (days 10–45) | En at day 45 |
|---|---|---|
| default (PCM flatten) — `ppm`, `plm`, `ppm_h4`, `pqm` | 0.353–0.355 /day | 1.6–1.7E-20 |
| `remap_boundary_extrap = .true.` — `ppm` / `plm` / `pqm` | 0.042–0.047 /day | 3.9–4.7E-25 |
| no remap at all (`vcoord_type = "lagrangian"`) | 0.046 /day | 1.5E-25 |
| `remap_method = "pcm"` (knob inert) | 0.875 /day | 7.9E-13 |

i.e. with the knob on, the mode's growth rate falls to the no-remap floor.
Gate: `test_remap_boundary_extrap` (four exactness cases FAIL with the
default closure, by construction).

#### Non-uniform-grid weights (`remap_nonuniform_weights`) — also orthogonal

The boundary closure above buys linear exactness only for a **uniform**
source column. PLM's `0.5·minmod(Δq_l, Δq_r)` slope and PPM's
`(7/12, −1/12)` edge estimate are the **equal-thickness specialisations** of
Colella & Woodward (1984) eqs (1.6)–(1.8), so on a **stretched** source
column — which is what every geometric family but σ-on-a-flat-bed hands the
remap — they carry an O(Δh/h) error through the whole **interior**, not just
at the two boundary cells. `"ppm_h4"` and `"pqm"` already carry
thickness-weighted stencils and are unaffected.

| Knob | Default | Effect |
|---|---|---|
| `&vcoord_nml remap_nonuniform_weights` | `.false.` (bit-identical) | `.true.` ⇒ PPM takes the CW84 (1.6) edge value on the true stencil thicknesses (with the unlimited (1.7) jump, so it reduces **exactly** to `(7/12, −1/12)` on a uniform column), and PLM takes the CW84 (1.7)+(1.8) h-weighted slope (MOM6 `PLM_slope_cw`). Inert for `"pcm"`; reaches `"ppm_h4"`/`"pqm"` only through their small-`nz` fallbacks. |

Measured on random stretched columns (`nz = 12`, `h` uniform on [0.5, 6.5] m,
`dq/dz = 0.25`, `remap_boundary_extrap = .true.`), max |q_new − q_exact| on a
profile linear in z:

| method | knob off | knob on |
|---|---|---|
| `plm` | 1.5E-01 | 2.8E-14 |
| `ppm` | 1.7E-01 | 4.3E-14 |
| `ppm_h4` | 5.0E-14 | 5.0E-14 (already exact) |
| `pqm` | 5.2E-14 | 5.2E-14 (already exact) |

**One caveat worth reading before turning it on.** PPM reduces exactly on a
uniform column, PLM does **not**: CW84 (1.7) there is the *centred*
difference under the (1.8) bound, where the shipped kernel uses the strictly
more diffusive minmod. So the knob swaps PLM's limiter as well as its
weighting and moves the `plm` answer even with no stretching — still
monotone, still conservative, now linear-exact. Gate:
`test_remap_nonuniform` (`uniform_source_plm_swaps_limiter` pins exactly
that, and `knob_off_is_first_order` pins the defect the knob removes).

## Tracer advection & reconstruction schemes

Tracer transport has three distinct reconstruction jobs — horizontal tracer
advection, vertical (in-z) advection, and the conservative ALE vertical remap.
They are *separate code paths* even when they share a paper, so this section
lists the primitives, then the regime × backend matrix, then the knobs.

### (a) Reconstruction primitives

All PPM rows cite Colella & Woodward 1984; MUSCL is Barth-Jespersen-limited
Green-Gauss (van Leer-class). "Order" is the formal spatial order of the
reconstruction.

| Primitive | file:line | Order | One-line | Paper |
|---|---|---|---|---|
| PCM remap (`remap_column_pcm`) | `src/ALE/rdb_remap_column.F90:73` | 0th (donor) | donor cell, diffusive, monotone | — |
| PLM remap (`remap_column_plm`) | `:131` | 1st | minmod-limited linear | van Leer 1979 |
| PPM remap (`remap_column_ppm`) | `:227` | 2nd | piecewise parabolic, C-W limiter | Colella & Woodward 1984 |
| PQM remap (`remap_column_pqm`) | `src/ALE/rdb_remap_column.F90` | 4th–5th | piecewise quartic, implicit h4 values + h3 slopes + W&A limiter; `N<5`→PPM | White & Adcroft 2008 |
| ALE remap dispatcher (`remap_column`) | `:48` | — | `select case(method)`; **unknown tag → PLM** (`:68-69`) | — |
| Ocean continuity PPM slope (`ppm_limited_slope`) | `src/core/ocean/kernels/continuity_ppm/rdb_continuity.F90:1645` | 2nd | van-Leer MC slope for the parabola | C-W 1984 |
| Ocean continuity PPM cell limiter (`ppm_cell_limiter`) | `:1669` | 2nd | C-W eq 1.10 monotonic limiter | C-W 1984 |
| Ocean positivity limiter (`ppm_limit_pos`) | `:1696` | — | positivity *modifier* on the PPM face values | — |
| Ocean positive-definite outflux limiter (`pd_limit_zonal_impl`/`pd_limit_meridional_impl`) | `src/core/ocean/kernels/continuity_ppm/rdb_continuity.F90` | — | per-donor θ scaling of outgoing layer mass flux ⇒ `h ≥ h_lim`, **zero mass created** (contrast: MOM6's injecting `max(h,Angstrom)` clamp is NOT ported); + a `2·h_lim` PPM edge floor (`:467/:514` parity) | — |

There are **two separate PPM implementations** for two different jobs — ALE
conservative remap (`rdb_remap_column.F90:227`) and ocean continuity/transport
face-reconstruction (`rdb_continuity.F90`). Separate code, same paper, **not
duplication** — each is shaped for its own loop/stencil/conservation
constraint.

### (b) Regime × backend matrix

| Horizontal tracer advect | Vertical (in-z) advect | ALE vertical remap |
|---|---|---|
| **PPM, HARDCODED (no knob)** (`tracer_advect_zonal_one_impl:1314`, `tracer_advect_meridional_one_impl:1384`; 2-cell near-wall band drops to 1st-order PCM/donor) | **1st-order UPWIND** (`src/core/ocean/kernels/vertical_advection/rdb_ocean_vertical_advection.F90:247` `tracer_advect_vertical_one_impl`, docstring "first-order upwind-in-z") | PCM/PLM/PPM via `remap_method` (default `"ppm"`); the parsed method is threaded through `rdb_ocean_setup.F90:500` → `vcoord%remap_method` → `rdb_ocean_dyn.F90:1061` |

### (c) Knobs

| Knob | file:line (decl) | Default | Selects | Scope |
|---|---|---|---|---|
| `remap_method` | `src/core/rdb_config.F90:2030` (enum `:4097`, parse `src/ALE/rdb_vcoord.F90:628`) | `"ppm"` (pcm/plm/ppm/ppm_h4/pqm) | vertical ALE remap order | the ocean setup parses it into `vcoord%remap_method` (`rdb_ocean_setup.F90:838`) and the ocean remap caller passes it (`rdb_ocean_dyn.F90:1061`) |
| `ppm_limit_pos` | `src/core/rdb_config.F90:569` (parse `:2967-2968`) | `.false.` | positivity-limiter **modifier** on the ocean continuity PPM face values — not a scheme-family selector | ocean continuity |
| `positive_definite` | `src/core/rdb_config.F90` (`register_ocean_continuity`) | `.false.` | positive-definite split continuity **modifier**: `2·h_lim` PPM edge floor + per-donor θ outflux limiter ⇒ every layer `h ≥ h_lim` with zero mass created; `h_lim = angstrom_h` on VCOORD_LAGRANGIAN else 0. Mutually exclusive with `&ocean_wetdry_nml enable` | ocean continuity |

### (d) Honesty notes (the bits a reader will get wrong)

1. **`remap_method` drives the ocean ALE remap.** The ocean ALE remap honors it:
   `parse_remap_method(cfg%remap_method)` → `vcoord%remap_method`
   (`rdb_ocean_setup.F90:500`), passed through to the remap call
   (`rdb_ocean_dyn.F90:1061`). The `m = REMAP_PPM` in `rdb_ocean_remap.F90:121`
   is only the *fallback default for the optional `method` argument* — the
   caller always supplies one. (The *horizontal* tracer advect is the only
   hardcoded-PPM path; the ALE remap is not.)

## Equation of state

| EOS | ocean | Knob |
|---|---|---|
| Linear two-tracer | `default` | `&ocean_eos_nml eos="linear"` |
| Wright (1997) nonlinear | available | `&ocean_eos_nml eos="wright"` |
| Roquet et al. (2015) SpV (TEOS-10-class) | available | `&ocean_eos_nml eos="roquet_spv"` (not with `fv_wright` PGF) |

Linear-EOS reference state — `ρ = ρ_0 + β_S·(S−S_ref) − α_T·(T−T_ref)`:

| Parameter | Knob | Default | Units | Test |
|---|---|---|---|---|
| `α_T` | `&ocean_ic_nml alpha_T` | `1.7e-4` | kg/m³ per °C | `test_ocean_linear_eos_knobs` |
| `β_S` | `&ocean_ic_nml beta_S` | `7.6e-4` | kg/m³ per PSU | `test_ocean_linear_eos_knobs` |
| `T_ref` | `&ocean_ic_nml T_ref` | `10.0` | °C | `test_ocean_linear_eos_knobs` |
| `S_ref` | `&ocean_ic_nml S_ref` | `35.0` | PSU | `test_ocean_linear_eos_knobs` |
| `ρ_0` | `&ocean_ic_nml rho_0` | `1035.0` | kg/m³ | `test_ocean_linear_eos_knobs` |

`α_T`/`β_S` are **DIMENSIONAL** (kg/m³ per unit), not the fractional
1/°C, 1/PSU coefficients protocols usually quote — multiply those by
`ρ_0` first (`alpha_T = ρ_0·α`). See `docs/REFERENCE.md` §Equation of
state. The coastal-legacy `&tracer_nml alpha_T/beta_S/T_ref/S_ref`
spellings are RETIRED and fail loud at configure.

## Initial conditions (analytical seeds)

| IC | ocean | Knob | Test |
|---|---|---|---|
| Uniform T / S | `default` | `&tracer_nml initial_temperature / initial_salinity` | `test_ocean_salinity_ic` |
| Linear-in-layer T(z) (bed `k=1` → surface `k=nz`) | available | `&tracer_nml T_init_bottom` + `T_init_surface` (both non-zero) | `test_ocean_eady_ic` |
| Linear-in-layer S(z) (bed `k=1` → surface `k=nz`) | available | `&tracer_nml S_init_bottom` + `S_init_surface` (both non-zero; one alone fails loud) | `test_ocean_salinity_ic` |

## Boundary tracer fluxes (ocean path)

| Flux | ocean | Knob | Test |
|---|---|---|---|
| Surface heat / salt (top layer `k=nz`) | scalar | `&ocean_thermo_nml q_heat / q_salt` | `test_ocean_surface_flux` |
| Surface-flux component set (`q_sw/q_lw/q_lat/q_sens/heat_added`, mass fluxes + `heat_content_*` enthalpy companions, `salt_flux`, `p_surf_atm`/`p_surf`; `Q_heat`/`Q_salt` become derived views assembled from the const + components) | per-component 2D | `&ocean_forcing_nml enable_components` (default off ⇒ no array allocated, `Q_heat`/`Q_salt` unchanged) | `test_ocean_surface_forcing_type` |
| Geothermal bottom heat (lowest massive layer, `k=1`) | scalar | `&ocean_geothermal_nml enable / q_geo` | `test_ocean_geothermal` |
| Surface buoyancy restoring (top layer `k=nz` T/S → scalar targets, piston velocity; non-conservative) | scalar | `&ocean_restore_nml enable_restore_temp / enable_restore_salt` | `test_ocean_restore` |
| Tidal body forcing (equilibrium tide `η_eq`; momentum feels `−g∇(η−η_eq)` folded into the barotropic PGF; equilibrium-argument astronomy, 10 constituents, needs lat/lon grid) | barotropic | `&ocean_tides_nml enable` | `test_ocean_tides_astronomy`, `test_ocean_tidal_forcing` |
| Barotropic linear wave drag (Egbert & Ray 2001; Jayne & St Laurent 2001) — static per-face piston velocity `r_H` MULTIPLIED into `bt_rem_u/v` inside the BT substep; `form="uniform"` (global scalar) or `form="roughness_proxy"` (resolved-bathymetry-variance placeholder for the real subgrid `⟨h²⟩`, pending PR-14); composes with `substep_drag` | barotropic | `&ocean_bt_nml wave_drag` | `test_ocean_wave_drag` |
| Atmospheric surface-pressure loading / inverse barometer (Wunsch & Stammer 1997) — `η_ib = −p_surf/(ρ₀·g_bt)` folded into the barotropic `eta_forcing` seam so momentum feels `−(1/ρ₀)∇p_surf` (~1 cm/hPa; a high depresses SSH); composes additively with the tide on the same seam; split-solver only, needs `enable_components` (reads `sf%p_surf`); ρ₀ from `eos%rho0` | barotropic | `&ocean_psurf_nml enable` | `test_ocean_p_surf`, `test_ocean_p_surf_bitident` |
| Reference density of the forcing terms — the `dt/(ρ₀·cp)` / `dt/ρ₀` divisors on EVERY surface heat/salt source (sea-ice coupling included), the wind-stress `τ/(ρ₀·h_top)`, the KPP `N²` / `u* = √(\|τ\|/ρ₀)` / kinematic fluxes behind `B_0`, and the geothermal `dt·Q_geo/(ρ₀·cp)` | not separately settable — all take `&ocean_ic_nml rho_0` | Fanned out from the ONE configured ρ₀ (`eos%rho0`) by `configure_ocean_reference_density`, joining the PGF and the twelve other slots that already took it. Before that wiring each slot kept a hard-coded 1035 that nothing assigned, so `rho_0 /= 1035` ran the EOS on the configured density and every forcing term on 1035, silently. `vmix%rho0` is the only one read on-device (inside `do concurrent`), and is correct because configure precedes `enter_data`. The implicit stress fold (`&ocean_vdiff_nml implicit_stress`) takes ρ₀ as an argument and now FAILS LOUD if it is omitted rather than defaulting to 1035. Default `rho_0 = 1035` ⇒ bit-identical. Test: `test_ocean_forcing_rho_ref` |
| Top-of-column pressure in the EOS's IN-SITU argument (E3, ice-shelf-cavity prerequisite) — the assembled `sf%p_surf` is copied once per outer step into `ms%p_top` (Pa, always allocated + zero-filled + device-mapped), and the ported in-situ builder evaluates at `p_top + hydrostatic` instead of starting at 0 Pa at the free surface. v1 ports exactly one consumer: the FV-Wright Picard column sweep (`form="fv_wright"`); with any other PGF form it is a documented no-op (warned at configure). It deliberately does NOT touch `ms%rho_layer` — a POTENTIAL density at the horizontally uniform `&ocean_eos_nml p_ref`, differenced along layers and vertically, so a per-column reference would fabricate an along-layer density gradient (N² therefore unchanged and self-consistent). EOS ARGUMENT only — the PGF top boundary condition (`pa(nz+1)`, `p_edge(nz+1)`) is unchanged. `validate_config` REFUSES `in_eos` with the unported in-situ builders (EPBL, kappa-shear, tidal mixing, Redi, isopycnal slopes, PGF in-layer reconstruction, sea ice) | EOS in-situ pressure | `&ocean_psurf_nml in_eos` (default off ⇒ byte-identical) | `test_ocean_eos_p_top` |
| Potential-density reference pressure — the pressure `ms%rho_layer` is referenced to (Wright / Roquet; the linear EOS ignores it). Previously a declared-but-never-assigned `eos%p_ref` permanently stuck at 0; now a namelist knob. HORIZONTALLY UNIFORM BY DESIGN (a per-column reference would fabricate an along-layer density gradient). Selects the thermobaric state at which the effective α/β are evaluated — relevant near the freezing point at cavity pressures. Distinct from `&vcoord_nml rho_ref_pressure` (the RHO/HYCOM target-density coordinate + density-space diag remap); normally set both to the same value for a density-coordinate run | EOS | `&ocean_eos_nml p_ref` (default 0.0 ⇒ byte-identical) | `test_ocean_eos_p_top` |
| Freezing-point (liquidus) coefficient SET — `eos_freezing_point` evaluates the linear `T_f = λ1·S + λ2 + λ3·p` for every EOS variant (MOM6 `TFREEZE_FORM="LINEAR"` parity), with the triple carried on the `eos_t` handle (`tfr_s`/`tfr_0`/`tfr_p`) instead of hard-coded, so the ice slot's handle (`engine%state%eos`, read by frazil / frazil uptake / basal flux) and the vmix / EPBL / kappa-shear / tidal-mixing copies all inherit it. `"seaice"` = SIS2/MOM6 (−0.054, 0, −7.53e-8); `"isomip"` = ISOMIP+ (−0.0573, 0.0832, −7.53e-8), Asay-Davis et al. 2016 Table 4 p. 2483 / eq. (25) p. 2485. The two are ~0.031 °C apart at S = 34.5 — enough to flip the SIGN of an ice-shelf basal melt rate, hence a named-set knob rather than a silent choice. NAMED SETS ONLY (no free-form λ knobs: the triple is a fit); an unrecognised name is a fail-loud `validate_config` error; a NONLINEAR form (Millero 1978, TEOS-10 polynomial) would arrive at the documented dispatch seam in `eos_freezing_point`, not as another set. Assigned in `configure_ocean_drag`, the earliest configure stage — before every handle copy and before `enter_data` | EOS / sea ice | `&ocean_eos_nml tfreeze_set` (default `"seaice"` ⇒ bit-identical at every shipped call site, which all pass `p = 0`) | `test_ocean_freezing_point` |

## Sea ice (ocean path only — `&ocean_ice_nml`, default off ⇒ byte-identical)

SIS2 port on the ocean C-grid (`src/core/ice/`, 14 modules). The `enable`
master switch gates the whole subsystem — off (default) ⇒ the slot is never
initialised, mapped, or stepped ⇒ byte-identical. **A high-quality sea-ice
dynamical core and column model, not yet a sea-ice model**: the rows below
that carry a knob are landed and independently tested; the rows marked `—`
are absences documented on purpose (see `docs/CAPABILITIES_AND_LIMITATIONS.md`
§"Sea ice" for the physical consequence of each). Envelope: the whole slot is
fail-loud single-rank (`enable=.true.` alone already requires `px*py=1`,
`rdb_config.F90` — not just `transport`/`dynamics`, which carry their own
identical guards on top); `transport` additionally forbids any periodic edge;
`dynamics` (EVP) allows wall/periodic edges only, no tripolar fold, no
OBC/tidal/sponge/clamped/Chapman edge. The authoritative knob set lives in the
`ocean_ice_config_t` / `&ocean_ice_nml` block of `src/core/rdb_config.F90`.

| Capability | ocean | Test |
|---|---|---|
| Winton (2000) two-layer column thermodynamics + enthalpy | `&ocean_ice_nml enable` (`nk_ice=2` Winton two-layer) | `test_ocean_ice_column`, `test_ocean_ice_enthalpy` |
| Multi-category ITD (thickness-space category restore) | `enable` + `ncat>1` | `test_ocean_ice_itd` |
| Category ice/snow transport + `compress_ice` (category-summed PPM) | `transport` (needs `ncat>1` + non-periodic edges; `adv_substeps` advective sub-iterations) | `test_ocean_ice_transport` |
| C-grid EVP dynamics (elastic-viscous-plastic momentum) | `dynamics` (`evp_sub_steps`; `p0`/`c0`/`ec`/`cdw`/`rho_ocean` strength+drag; `del_sh_min_scale`; `tdamp`; `a_face_stress` momentum-conserving wind+drag weighting; `cfl_trunc`/`cfl_trunc_dyn_its` transport-CFL velocity ceiling; `project_ci` in-loop concentration projection) | `test_ocean_ice_evp` |
| Ridging / rafting | **`—` (not available)** — `compress_ice` is **area compaction** = SIS2's own `DO_RIDGING=.false.` fallback, not a participation/redistribution scheme; convergent-regime ITD (and the EVP strength that reads it) is biased. SIS2's own ridging option is an Icepack wrapper (`ice_ridge.F90`), default off | — |
| Snowfall source (`m_snow` accumulation; PR 26) | `&ocean_ice_nml snowfall` (requires `enable` + `&ocean_thermo_nml enable_thermodynamics`) | `test_ocean_ice_snowfall` |
| Snow-ice flooding (Archimedes freeboard, SIS2 `SN2IC`; PR 27) | `&ocean_ice_nml snow_ice` (requires `enable` + `&ocean_thermo_nml enable_thermodynamics`) | `test_ocean_ice_column` (`snow_ice_*` cases) |
| Melt ponds | **`—` (not available)** — `m_pond = 0.0_wp` dead local (`rdb_ice_column.F90`) | — |
| Lateral melt / floe size | **`—` (not available)** — surface + basal melt only; MIZ retreat biased | — |
| Frazil (bank + uptake) | `enable` (no knob of its own) — **surface-layer only** (`k=nz`); MOM6/SIS2 check the full column; Boussinesq virtual salt flux (ocean water mass unchanged) | `test_ocean_frazil`, `test_ocean_ice_coupling` |
| Ice atmospheric forcing | `&ocean_ice_nml air_temp` / `restore_lambda` / `sw_down` — **v1 stub**: scalar, uniform, slab-restoring; no bulk formulae, no 2-D fields, no file input, no time dependence | `test_ocean_ice_driver_column` |
| Ice → ocean coupling | `enable` — salt ✓ closed (virtual); heat ✓ closed (incl. transmitted shortwave `sw_thru`, PR 31 — `ice_ocean_sw_flux` delivers it to `q_sw`/`Q_heat`); momentum ✗ **not conserved at fractional cover** (D7); freshwater/mass ✗ **not coupled at all** — ice carries no weight (no dynamic sea-surface loading) | `test_ocean_ice_coupling`, `test_ocean_ice_driver_column` |
| Ice initial condition (analytic seeding: seeds a live pack before the first step — no frazil growth required; mass binned into the ITD via `ice%mh_lim`, `enth_ice`/`enth_snow` set from `t_ice`/`s_ice` via the exact `ice_enth_from_ts` inversion) | `&ocean_ice_ic_nml conc_config` = `"zero"` (default, no-op) / `"uniform"` (scalar `h_ice`/`conc`/`h_snow`/`t_ice`/`s_ice`) / `"latitudes"` (SIS2 polar-cap 0/1 step off `geolatT`, needs a non-cartesian grid). v1 analytic-only — file-backed ICs deferred to PR-14 | `test_ocean_ice_init` |

Notes: `transport` is fenced against periodic edges (fail-loud at configure);
EVP `dynamics` is not — it runs under wall or periodic edges. With
`dynamics=.false.` ice velocity falls back to the ocean-surface-layer sampler;
with it on, `rdb_ice_evp` writes `ice%u_ice`/`v_ice` from the momentum solve.
`tdamp < 0` is the SIS2 special case `|tdamp|·dt_slow` for the elastic damping
timescale. `snowfall` spreads a uniform frozen-precipitation rate onto
`ice%atm_fprec` (v1 slab filler, same seam as `air_temp`/`sw_down`); the
ice-free share (open water at `ncat>1`, ice-free cells at `ncat==1`, and any
category failing the column's own entry gate) is delivered to the ocean as a
latent-heat sink + virtual freshening via `rdb_ice_snow%ice_snowfall_ocean_
share` — a DELIBERATE divergence from SIS2, which instead orphans that share
as snow on a zero-ice category (a state Roundabout's `rdb_ice_transport` treats as
a fail-loud violation). `cfl_trunc` (SIS2 `CFL_TRUNCATE`, default 0.5 there,
`0` here ⇒ byte-identical) clips the final ice velocity to the transport-CFL
bound, demoting `ice_transport_step`'s conservation/positivity abort to a
backstop rather than the only defence; `project_ci` (SIS2
`PROJECT_ICE_CONCENTRATION`, default `.true.` there, `.false.` here) projects
`ci` forward within the EVP subcycle loop so `pres_mice` stiffens under
convergence instead of holding the pre-loop value for the whole call — inert at
a saturated (`ci=1`) jam by construction. `snow_ice` converts submerged snow
mass to the top ice layer in ONE non-iterative Archimedes step
(`ice_snow_ice_flood`, `rdb_ice_mass`), exactly conserving column
mass/enthalpy/salt; it exchanges NOTHING with the ocean (SIS2's own
formulation — `snow_to_ice` is diagnostic-only). DIVERGENCE FROM PHYSICAL
REALITY (inherited from SIS2, not fixed here): the converted ice carries the
SNOW's enthalpy and ZERO salinity, not the flooding SEAWATER's — snow enthalpy
is far more negative than near-freezing seawater, so the new ice is too COLD
and too FRESH. True seawater flooding (drawing ocean mass + its
enthalpy/salinity into the pore space) is out of scope. Deferred: multi-rank,
ridging / mechanical
redistribution beyond `compress_ice`, real (non-virtual) freshwater mass for
the snowfall ocean share (PR-16). The `—`-in-the-ocean-column rows above are
documented absences, not omissions — an undocumented gap and a documented one
look identical from outside, which is why each has a row.

---

## Diagnostics (no feedback on the dynamics)

| Diagnostic | ocean | Numerics | Primary test |
|---|---|---|---|
| First-baroclinic wave speed + Rossby radius (B1; Chelton et al. 1998) | `&ocean_wavespeed_nml enable` (default off) | per-column rigid-lid Sturm–Liouville eigensolve from `rho_layer`; backtracking convective merge + fixed-budget Sturm-count bisection; called once per outer step (thermo-cadence-gated, further gated by `n_wavespeed`) in `ocean_dyn_step_split`, BEFORE `varmix_compute`/`run_meke_step` — feeds B2 (GM/Redi/MEKE resolution scaling). `rd_over_dx = rd/metrics%dxT` (metres, not `grid%dx` — degrees on spherical/supergrid/tripolar); `f_centre`/`beta_centre` come from `metrics_fill_coriolis` (planetary or beta-plane), not a hard-coded beta-plane. `ocean_dyn_step` (the unsplit path) does not call it — trap #5, PR-3. | `test_ocean_wave_speed` |

---

## Tunable parameters (key dials per scheme)

The *curated* dial set — defaults are the namelist defaults (`0` / `.false.` =
off/inert). The exhaustive per-field list (units, CF metadata) is the `config_t`
docstrings (FORD) + [`REFERENCE.md`](REFERENCE.md).

**Porous barriers** — `&ocean_porous_nml` (ocean only)

| Knob | Default | Meaning |
|---|---|---|
| `enable` | `.false.` | Master switch. Off ⇒ the open-area fields stay at their `(1,1,1)` placeholder and no porous kernel is launched (byte-identical). Fails loud with `&ocean_bt_nml bt_halo > 0` and with `&ocean_wetdry_nml enable`. |
| `source` | `"resolved"` | Along-face statistics source. `"resolved"` = wet-gated corner/midpoint samples of the RESOLVED bathymetry (a proxy). `"file"` = offline subgrid file — fails loud, deferred. |
| `eta_interp` | `"max"` | Interface height at the velocity point (MOM6 `PORBAR_ETA_INTERP`): `max` (the higher/shallower interface — the LEAST blocking, since `w` increases with height), `min` (the most blocking), `arithmetic`, `harmonic`. |
| `masking_depth` | `0.0` m | Faces whose mean along-face depth is SHALLOWER than this stay fully open (MOM6 `PORBAR_MASKING_DEPTH`, positive below the surface). |

**Implicit stress/drag fold** — `&ocean_vdiff_nml` (ocean path only; folds the
surface wind stress + bottom drag into the backward-Euler vertical-friction
tridiagonal as BCs instead of explicit pre-solve adds — thin-layer
(z*/ZSTAR_FULL pinch-out) CFL-robust). Both default off ⇒ bit-identical.
| Knob | Default | Meaning |
|---|---|---|
| `implicit_stress` | `.false.` | Wind stress → vdiff surface (`k=nz`) RHS row (Neumann top-BC). Incompatible with `&ocean_vmix_nml direct_stress`. |
| `implicit_drag` | `.false.` | Bottom drag → vdiff bed (`k=1`) diagonal (stress bottom-BC). Mutually exclusive with `&ocean_bdrag_nml implicit` and HBBL (`hbbl>0`). |

**Barotropic linear wave drag** — `&ocean_bt_nml` (ocean path only; bulk
energy sink for the barotropic tide, Egbert & Ray 2001 / Jayne & St Laurent
2001). `wave_drag=.false.` (default) ⇒ bit-identical. `form="uniform"` needs
`wave_drag_r_uniform>0` (a warning fires if not). `form="roughness_proxy"` is
a documented PLACEHOLDER — resolved-bathymetry `⟨grad b⟩²` variance, not the
real subgrid `⟨h²⟩` (pending PR-14's file reader / PR-30). `form="file"` is
registered but fails loud at configure (not implemented — PR-14).
| Knob | Default | Meaning |
|---|---|---|
| `wave_drag` | `.false.` | Master switch. |
| `wave_drag_form` | `"uniform"` | `"uniform"` \| `"roughness_proxy"` \| `"file"` (fail-loud, unimplemented). |
| `wave_drag_scale` | `1.0` | Global tuning multiplier on `r_H`. |
| `wave_drag_r_uniform` | `0.0` | Piston velocity `r_H` (m/s) for `form="uniform"`. |
| `wave_drag_kappa` | `6.2832e-4` | Topographic wavenumber (1/m) for `form="roughness_proxy"` (matches `&ocean_tidal_mixing_nml kappa_itides`). |
| `wave_drag_n_bot` | `1.0e-3` | Reference bottom `N` (1/s) for `form="roughness_proxy"`. |
| `wave_drag_h2_max` | `2.5e4` | Ceiling on the `⟨h²⟩` proxy (m²) for `form="roughness_proxy"`. |

**PP81 + KPP constants** — `&ocean_vmix_nml` (routes to
`ocean_state%vmix`). `kv_bg`/`kt_bg`/`ks_bg` (the assembly floor,
`vmix_assemble`) are re-derived from `pp81_nu_bg`/`pp81_kappa_bg` at
configure time (`vmix_seed_backgrounds`) so the floor tracks a user-set
background rather than the type default.
| Knob | Default | Meaning |
|---|---|---|
| `pp81_nu0` | `1e-2` | PP81 Richardson-dependent viscosity scale (m²/s). |
| `pp81_nu_bg` | `1e-4` | PP81 background viscosity (m²/s); also seeds `vmix%kv_bg`. |
| `pp81_kappa_bg` | `1e-5` | PP81 background diffusivity (m²/s); also seeds `vmix%kt_bg`/`ks_bg`. |
| `pp81_alpha` | `5.0` | Ri scaling in `ν = ν_bg + ν₀/(1+α·Ri)²` (paper value 5; implementations vary 4-10). |
| `shear2_floor` | `1e-10` | Floor on `\|∂u/∂z\|²+\|∂v/∂z\|²` in the PP81 Ri denominator (1/s²). |
| `kpp_ri_crit` | `0.3` | Critical bulk Richardson number for the KPP BL-depth sweep (LMD94 §3). |
| `kpp_cs_nonlocal` | `6.3` | Non-local (counter-gradient) transport coefficient `C_s` (LMD94 eq 20). |
| `kpp_c_vt2` | `1.8` | Unresolved-shear `V_t²` coefficient (LMD94 eq 23); `0` disables `V_t²`. |

**Along-coordinate tracer Laplacian (ocean path)** — `&ocean_hdiff_nml`
(`rdb_ocean_hdiff_tracer`). Diffusion along the model coordinate, NOT along
neutral surfaces (`&ocean_redi_nml` is the separate neutral/isopycnal path).
Default `kappa_h = 0.0` keeps the kernel's short-circuit intact ⇒
bit-identical. Configure-time guard aborts if
`kappa_h·dt_therm·(1/dx²+1/dy²) > 0.5` (Cartesian-metres check; skipped on
spherical/tripolar grids and when `dt_fixed` is adaptive).
| Knob | Default | Meaning |
|---|---|---|
| `kappa_h` | `0.0` | Horizontal tracer diffusivity (m²/s); `0` = no-op. |

**EPBL** — `&ocean_epbl_nml` (ocean path only; replaces the KPP overlay, PP81
interior + background continue underneath).  Full knob table + the MOM6 name
mapping: [`docs/generated_nml_knobs.md`](generated_nml_knobs.md); the
load-bearing dials:
| Knob | Default | Meaning |
|---|---|---|
| `enable` | `.false.` | Master switch (requires `use_closure` + thermodynamics; logs + disables `use_kpp`). |
| `mstar_scheme` | `"om4"` | Mechanical-TKE efficiency: `constant` / `om4` / `rh18`. |
| `mstar` | `1.2` | Constant-scheme mstar. |
| `nstar` | `0.2` | Convective-PE → TKE efficiency. |
| `tke_decay` | `2.5` | Ekman-depth / TKE-decay-scale ratio. |
| `mld_iteration` | `.true.` | Self-consistent MLD root-find (false position, `mld_tol = 1 m`). |
| `translay_scale` | `0.1` | Transition-layer mixing-length floor. |
| `prandtl` | `1.0` | `Kv = prandtl · Kd` into the momentum solve. |
| `combine` | `"add"` | Fold into PP81 kv/kt additively or by `max`. |
| `tke_diags` | `.false.` | Per-column TKE budget terms (close to round-off; asserted by `test_ocean_epbl`). |
| `use_lt` | `.false.` | Langmuir enhancement of mstar — LF17 wind-only statistical waves (one La per column from u* + MLD; no wave model). `lt_scheme = "rescale"\|"additive"`, Reichl & Li (2019) coefficients, `lt_lac1..5` stability modification. |

**kappa-shear** — `&ocean_kappa_shear_nml` (ocean path only; interior
closure, coexists with KPP/EPBL — its κ is ADDED to `kt` and
`prandtl_turb·κ` to `kv` every stage; the column solve runs at thermo
cadence).  Full knob table + the upstream-name mapping:
[`docs/generated_nml_knobs.md`](generated_nml_knobs.md); the
load-bearing dials:
| Knob | Default | Meaning |
|---|---|---|
| `enable` | `.false.` | Master switch (requires `use_closure` + thermodynamics; no exclusions). |
| `ri_crit` | `0.25` | Critical Richardson number for the shear source. |
| `shearmix_rate` | `0.089` | Source-rate coefficient (JHL08 calibration). |
| `kappa_0` | `1.0e-7` | Background κ floor / well-posedness smoothing (m²/s). |
| `tol_err` | `0.1` | Picard convergence tolerance (also scales the adaptive-dt bands). |
| `max_inner_it` / `max_substep_it` | `50` / `13` | Inner Picard / outer substep caps. |
| `prandtl_turb` | `1.0` | `Kv = prandtl_turb · Kd` on the momentum side of the merge. |
| `massless_merge` | `.false.` | D4: fold vanished (`< H_VANISHED`) layers onto the column's massive sub-grid before the solve (vs the blunt gather floor); identity columns bypass, so on-path is bit-identical on healthy envelopes. |
| `at_vertex` | `.false.` | MOM6 `VERTEX_SHEAR` (OM5-class production form): solve the columns at C-grid CORNERS from the native face velocities (no centre-average shear damping), then average corner Kd back to tracer points.  Kv is routed corner→face (`prandtl_turb·kd_corner` into `vdiff_apply_momentum`'s `kv_corner_source`, MOM6 `Kv_shear_Bu`) with the cell-centred kv merge suppressed; `tke_int` zeroed in vertex mode. |
| `vertex_geometric_mean` | `.false.` | Geometric (vs arithmetic) corner→centre Kd mean (MOM6 `VERTEX_SHEAR_GEOMETRIC_MEAN`). |
| `vertex_geomean_kdmin` | `0.0` | Floor (m²/s) on each corner Kd BEFORE the geometric mean (MOM6 `VERTEX_SHEAR_GEOMETRIC_MEAN_KDMIN`; inert unless geometric).  With 0 the geometric mean hard-zeros Kd at every shear-zone edge; OM5 configs use `1e-9`. |

**Vmix assembly gate** — `&ocean_vmix_nml` (ocean path only; the single
downstream floor/clip/smooth/guard stage `vmix_assemble`, run once per stage
after every closure has contributed into kv/kt and `vmix_split_kd_heat_salt`
has derived ks from kt).  All defaults are bit-identical to the pre-assembly
chain.
| Knob | Default | Meaning |
|---|---|---|
| `kv_max` | `huge` (off) | Ceiling on momentum viscosity kv (MOM6 `Kd_max` momentum). |
| `kd_max` | `huge` (off) | Ceiling on tracer diffusivity kt/ks (MOM6 `Kd_max`). |
| `kd_smooth_iterations` | `0` | 1-2-1 horizontal smoothing passes on kv/kt/ks (MOM6 `Kd_smooth`). |
| `vmix_guard` | `.false.` | Debug-gated negative/NaN diffusivity guard (error-stop, or status-returning testable path). |

**Geothermal bottom heat flux** — `&ocean_geothermal_nml` (ocean path only;
bed-side analogue of the surface heat flux, deposited into the lowest massive
layer)
| Knob | Default | Meaning |
|---|---|---|
| `enable` | `.false.` | Master switch (default off ⇒ bit-identical to no geothermal). |
| `q_geo` | `0.0` | Constant bottom heat flux (W/m², positive into the ocean; typical ~0.05–0.1). |

**Surface buoyancy restoring** — `&ocean_restore_nml` (ocean path only;
relaxes the surface top layer T/S toward scalar targets at rate
`λ = piston/h_top`; piston-velocity form, MOM6 `RESTOREBUOY`-inspired;
non-conservative — increment mirrored into the surface heat/salt budget)
| Knob | Default | Meaning |
|---|---|---|
| `enable_restore_temp` | `.false.` | Master switch for SST restoring (effective only when `piston_t /= 0`). |
| `enable_restore_salt` | `.false.` | Master switch for SSS restoring (effective only when `piston_s /= 0`). |
| `piston_t` | `0.0` | SST piston velocity (m/day, MOM6 `FLUXCONST_T`; converted to m/s at seed). |
| `piston_s` | `0.0` | SSS piston velocity (m/day, MOM6 `FLUXCONST_S`). |
| `restore_sst` | `0.0` | Scalar target SST (degC). |
| `restore_sss` | `0.0` | Scalar target SSS (PSU). |

**Ocean** dials live in the per-concern `&ocean_<group>_nml` sub-namelists
(`vmix`, `hvisc`, `bdrag`, `pgf`, `coriolis`, …) — see the working envelope in
[`CAPABILITIES_AND_LIMITATIONS.md`](CAPABILITIES_AND_LIMITATIONS.md).

<!-- TODO: windowed tracer advect (dt_tracer_advect_ratio) — add a row to a
     tracer-advection cadence section when the tracer-advection matrix branch
     (docs/tracer-advection-matrix) merges; the scheme lives in
     continuity_t (rdb_continuity.F90) and is gated by &ocean_vmix_nml. -->

## Keeping this current

1. **Per-capability PR step.** The CLAUDE.md "adding a new capability" workflow
   gains a step: *namelist knob → kernel → unit test → **closure-matrix row** →
   pre-commit + fortitude → ctest → commit.* One row per scheme is cheap.
2. **Automated gate.** `tools/check_closure_matrix.py` parses the closure knobs
   from `src/core/rdb_config.F90` namelists and the ocean sub-namelists, asserts
   each appears somewhere in this file, and checks that every test named here
   exists in `tests/CMakeLists.txt`. Wired into pre-commit + CI so drift fails
   the build instead of rotting silently (this matrix was born because the prose
   in `CAPABILITIES_AND_LIMITATIONS.md` drifted — it still claims PP81 is the
   only coastal mixing scheme, three closures ago).
