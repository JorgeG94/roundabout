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

Last reconciled with the code: **2026-09-22** (`&vcoord_nml zfixed_closed_faces` now REFUSES `&ocean_bt_nml correction_bc_pgf` / `substep_drag` / `wave_drag`, the three barotropic paths still on full-column weights; `substep_drag` under a non-linear `&ocean_bdrag_nml form` WARNS; `apply_bt_correction` takes `metrics` as a REQUIRED argument.  Before that, 2026-09-21: the shared FIRST-LIVE-LAYER index `multilayer_state_t%k_top` — the `Z-fixed (gprime)` vcoord row's cavity envelope rewritten: melt and ice-shelf top drag are no longer refused under `z_fixed`, and `&ocean_hdiff_nml kappa_h /= 0` is accepted whenever `&vcoord_nml zfixed_closed_faces` is on, because the face mask already zeroes every hdiff flux touching a filler; `&ocean_kappa_shear_nml enable` and `&ocean_tidal_mixing_nml enable` newly REFUSED there, both found unfenced in the consumer survey; earlier: 2026-09-20 (ice-shelf cavity REAL FRESHWATER MASS — `&ocean_cavity_melt_nml freshwater`/`volume_compensation` added, the `Virtual salt flux` row narrowed to the default and two new rows added for the mass form and its sea-level sink, and the `Budgets` row extended with the tracked MASS source `ms%mass_src`; earlier the same day: ice-shelf cavity DIAGNOSTICS — the `Diagnostics` row flipped from `not implemented` to thirteen catalog entries with the NaN-outside-the-cavity convention and the fail-loud prerequisite gate; earlier the same day: ice-shelf cover MASK on the atmospheric forcing — the `Cover mask on ATMOSPHERIC forcing` row flipped from `not implemented — REFUSED` to shipped, and the four refusals it stood on (wind stress, restoring, shortwave penetration, scalar `q_heat`/`q_salt`) became acceptances; earlier the same day: ice-shelf cavity basal melt COUPLED — `&ocean_cavity_melt_nml` fills the previously-empty knob column of the kernel-only section, and the far-field sampling / interface-pressure / flux-delivery / virtual-salt / budget / status rows were added for the coupling itself, with the melt liquidus becoming the THIRD consumer of `ms%p_top`; earlier the same day: the ice-shelf cavity LOAD partition wired end to end — `ms%p_top = p_ice_ref + sf%p_surf`, the `p_top_in_bc` requirement for a varying draft, and the measured resting-cavity residual; earlier the same day again: the basal-melt section existed as kernel-only, with every ocean cell reading `—`; earlier: 2026-09-14 (outer time-split section added — `&ocean_bt_nml split_scheme` had no row at all, and the fast-loop Coriolis-reference seam it selects had none either; earlier: 2026-08-23 (coastal-regime carve-out — the coastal-S / coastal-U columns, the coastal-only closure rows and the `&physics_nml` / `&nonhydrostatic_nml` coastal knob tables were removed with the coastal path; earlier: 2026-08-04 vertical-mixing section; background-mixing row added — Bryan-Lewis had shipped with no matrix row, Henyey lands in the same PR that fixes the gap; row rewritten when Henyey was brought to MOM6 parity: Bryan-Lewis XOR Henyey instead of Henyey-requires-Bryan-Lewis, `bkgnd_kd_min` floor added)).)

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

## Ice-shelf top drag (`&ocean_tdrag_nml`; ocean only, default off)

Kernel + slot: `src/parameterizations/vertical/rdb_ocean_top_drag.F90`
(`ocean_top_drag_t`), computed and applied beside the bottom drag in
`run_stage` / `run_stage_split` and added into the barotropic slow forcing by
`add_top_drag_into_F_slow`.  The mirror of the bottom drag about the middle of
the column: the sink lives at `k = nz`, masked by ice cover on FACES.  Requires
`&ocean_cavity_dyn_nml enable` (fail-loud) — without a draft every face mask is
zero.  Default off ⇒ bit-identical.

**The face cover rule is OR:** `cover_u(i,j) = max(cover_frac(i−1,j),
cover_frac(i,j))`, so the CALVING-FRONT face feels the drag.  AND would leave a
free-slip band exactly where the cavity outflow leaves.  **One `C_d`:** with
`&ocean_cavity_melt_nml` on, the melt `u*` takes its coefficient from
`&ocean_tdrag_nml cd` and a disagreeing `cdrag_top` is refused at configure —
MOM6 carries two independent coefficients, this model deliberately does not.

| Form | ocean | Notes |
|---|---|---|
| Quadratic log-layer (Cd) | `default` when enabled | ISOMIP+ `C_d = 2.5e-3` (Asay-Davis et al. 2016, Table 4, which prescribes the same quadratic law top and bottom); `du/dt = −C_d·\|U\|·u/h_nz`.  Test: `test_ocean_top_drag` |
| Linear (Rayleigh) | `form="linear"` | rate `r` (1/s); the form with a closed-form spin-down.  Test: `test_ocean_top_drag` |
| HTBL-distributed | `htbl` | spreads the stress over the top `htbl` m — the `hbbl` mode reflected; keeps the explicit rate finite where sigma thins the top layer near a grounding line (Killworth & Edwards 1999).  Distributed linear rate `r·htbl/h_nz` asserted in closed form |
| Implicit (backward-Euler in the drag kernel) | `implicit` | tendency formed as `−λ·u/(1 + dt·λ)` so the ordinary apply gives `u/(1 + dt·λ)`; unconditionally stable for any top-layer thickness.  Test asserts it stays monotone and same-signed at `dt·λ = 2.5`, where the explicit form amplifies and flips sign |
| Implicit (fold into the vdiff `k=nz` diagonal) | `&ocean_vdiff_nml implicit_top_drag` | the OTHER implicit form, and a different thing: `dt·λ_top` on the surface diagonal of the backward-Euler vertical-friction tridiagonal, so the drag is solved TOGETHER with the interior shear rather than ahead of it.  The wind stress already owns that row's RHS — a drag is a diagonal term and a stress is an RHS term, so they compose — and on a covered face the wind RHS is scaled by `(1 − cover)`: no atmosphere under a shelf.  The explicit apply is gated off.  Layer-`nz` only (`htbl > 0` refused, the mirror of the `implicit_drag`/HBBL restriction) and mutually exclusive with `implicit` above (double count).  Test: `test_ocean_top_drag` (`u/(1 + dt·λ)` at `dt·λ = 8` on a 5 cm pinched top layer; open face takes the full wind, covered face takes none) |
| Top-drag stress magnitude | (always, when enabled) | `stress_top(nx,ny)` = cell-centred `\|τ_top\|` (N/m²), device-resident.  Copied INLINE by `run_stage` / `run_stage_split` into `surface_stress%stress_shelf` in the same stage that computes it, and read there by KPP and EPBL as `u_*² = (stress_mag + stress_shelf)/ρ₀` — the total upper-boundary momentum flux, with disjoint supports (`stress_mag ≡ 0` under cover, `stress_top ≡ 0` off it).  No lag.  Tests: `test_ocean_bl_under_ice` |
| Under-ice `u_*` for the boundary-layer schemes | `&ocean_tdrag_nml enable`, else `&ocean_cavity_melt_nml enable` | KPP (`rdb_ocean_vmix`, two sites) and EPBL (`rdb_ocean_epbl`) take `u_* = √((stress_mag + stress_shelf)/ρ₀)`.  `stress_shelf` is ALWAYS allocated and zero without a cavity ⇒ one code path, no optional dummy, `x + 0.0 ≡ x` ⇒ bit-identical.  Preferred source is the top drag (in-stage); with the top drag OFF and melt ON, `engine_step_finalize` fills it from `ρ₀·u_*²` with the melt slot's own `u_*` — the SAME `C_d` under the one-coefficient rule, but at thermo cadence, so **that** path lags one outer step.  Tests: `test_ocean_bl_under_ice` |

## Porous barriers (subgrid topography; ocean only)

Represent a subgrid sill/strait by reducing the **open** fraction of a C-grid
face below its full width, per layer, so a deep sill blocks the bottom layers
while the surface layers stay fully open.

| Scheme | ocean | Numerics | Primary test |
|---|---|---|---|
| Porous barriers (Adcroft 2013 three-parameter fit) | `&ocean_porous_nml enable` | per-face `d_min`/`d_max`/`d_avg` along-face topographic heights ⇒ monotone open-width profile `w(η)`; the layer-averaged OPEN-AREA fraction is the exact difference of the profile's vertical integral over the layer, `min(1, (A(η_hi)−A(η_lo))/(η_hi−η_lo))`. A DEGENERATE face (`d_max ≤ d_min`) and a VANISHED layer (`dz ≤ H_VANISHED`) are special-cased — see the caveat below. Recomputed on the device ONCE PER OUTER STEP and held across both RK2 stages (MOM6's cadence), then MULTIPLIED into the per-layer transports of continuity-PPM (both the fused and the direction-split forms, both directions, before the barotropic renormalisation, which carries the open fraction in BOTH its `sum_h` denominator and its per-layer increment) and of the Coriolis/advection **transport** forms (`sadourny_energy`, `sadourny_hk`). The default velocity-form `sadourny` has no `u·h·dy_cu` transport to narrow, so porosity reaches it only through continuity. The BAROTROPIC substep transports on `dy_cu_bt`/`dx_cv_bt` — the same widths scaled by the COLUMN-INTEGRATED open fraction, which is identically the thickness-weighted mean of the per-layer fractions (an exact telescoping identity, both being the same integral of `w` over the column). Without that the barotropic solve would be porous-blind and the per-layer renormalisation to `uhbt` would hand the blocked transport straight back, leaving the barrier a vertical redistribution that never reduces net flow. **This is NOT `BT_cont` parity**: MOM6's production barotropic face area is `Σ_k (dy_Cu·por_k)·h_marginal_k·visc_rem_k` (PPM marginal thickness AND `visc_rem`, neither of which appears here); its `Σ_k h_k·(dy_Cu·por_k)` form is the open-boundary-segment branch only, and `set_local_BT_cont_types` carries no `por` at all. NOT applied to: `areaCu`/`areaCv` (so BT Coriolis + KE use un-narrowed areas against narrowed transports, as in MOM6), GM/Redi/MEKE/hdiff, or sea ice. **Fail-loud exclusions** (`validate_config`): `&ocean_bt_nml bt_halo > 0` (the wide-halo BT clone re-fills its own `metrics_w` from the grid formula and nothing gives it the porous statistics, so the wide fast loop would silently transport on UN-narrowed widths) and `&ocean_wetdry_nml enable` (the vanishing-column interaction with the wet/dry outflow limiter is unvalidated). Single-rank. Default off ⇒ no kernel launch at all ⇒ byte-identical | `test_ocean_porous` |

## Partial-step z-level face closure (`vcoord_type = "z_fixed"`; ocean only)

The other per-layer face gate, and **independent** of the porous one: porous
barriers narrow a face continuously for topography the grid cannot resolve;
this closes a face completely for a layer the coordinate has already declared
dead.  Effective width is the product, `dy_eff(I,j,k) = dy_cu(I,j) ·
por_face_area_u(I,j,k) · open_u(I,j,k)`.

| Scheme | ocean | Numerics | Primary test |
|---|---|---|---|
| Partial-step z-level face closure (Adcroft, Hill & Marshall 1997; Losch 2008 §2.1) | `&vcoord_nml zfixed_closed_faces` | Under `z_fixed` a layer whose nominal geopotential range lies inside the bed or the ice draft is an inert FILLER of thickness `zstar_h_min` (`≤ H_VANISHED`).  A velocity face at which that layer is a filler on EITHER side is a **WALL** for that layer, not a thin passage: no normal velocity, no mass or tracer flux, FREE-SLIP.  Left open, the FV pressure gradient integrates across a staircase step of height up to `h_nominal` and drives `\|ρ′\|·g·Δz_step/(ρ₀·dx)` out of a resting stratified state — thickness-independent, so no `h`-gate reaches it.  The mask `metrics%open_u/open_v` (`(nx+1,ny,nz)`/`(nx,ny+1,nz)`, 0/1) is built ONCE at configure by the `pure` builder `ocean_vcoord_closed_face_masks` from `ocean_vcoord_z_fixed_target` at `η = 0` — the same kernel the ALE regrid and the IC seed use, so "live" has ONE definition — and is STATIC (bed and draft are static; `η` is absorbed by the first LIVE layer).  Array-edge faces stay open, as in the porous kernel; a periodic seam is an INTERIOR index and is covered because the target is built from ghost-filled `bt_H_ref`/`z_top`.  Consumers: the four continuity-PPM flux sites; the barotropic renormaliser's weight `wk = dy_cu·por·open` (so `uhbt` is distributed over the OPEN layers only and `Σ_k mass_flux = uhbt` stays exact); the transport-Coriolis mass fluxes; `mask_layer_velocities` after every stage; `&ocean_hdiff_nml kappa_h`; the HARMONIC velocity-Laplacian viscosity (each neighbour difference gated by the NEIGHBOUR's open flag, the tendency by the face's own — otherwise a zeroed closed-face velocity acts as Dirichlet-0, i.e. NO-slip drag); the momentum vertical-friction tridiagonal (face column `min(h_L,h_R)`, coupling CUT across any interface touching a closed layer — the decoupling the tracer matrix already had); and the ALE face-velocity remap (same `min` face column, closed layers dropped from source AND target).  Tracer advection rides the already-masked `mass_flux_*_layer` and needs nothing, and `mask_layer_velocities` is re-asserted after the ALE remap (the remap is the last velocity writer of the step).  **Barotropic consistency** is part of the closure, not an extra: `derive_bt_from_layers`, `face_depth_mean_u/v` (+ the `visc_rem` twins), `set_cor_ref_velocity` (the `pred_corr` fast-loop Coriolis reference, MOM6 `ubt_Cor`) and `apply_bt_correction` all weight by `h_face·open`, so `ubt` is the OPEN-column depth mean and a CLOSED layer receives no barotropic increment BY CONSTRUCTION (not by a mask cleaning up after a uniform fold).  `metrics` is a REQUIRED argument of all four, not an optional one: it WAS optional, `set_cor_ref_velocity` omitted it, and the full-column mean it silently took returned `φ·ū_open` instead of `ū_open` — leaving `Δa_u = +f·(1−φ_v)·v̄` forcing EVERY barotropic substep proportionally to the barotropic velocity itself; the OPEN layers keep whichever distribution `use_h_weighted` selected, because forcing the h-weighted form concentrates the whole increment into the thickest open layer and was measured to blow the sloping-lid cavity up on day 2.  `dy_cu_bt`/`dx_cv_bt` carry `dy_cu·(Σ_k h_face·por·open)/(Σ_k h_face)`, recomputed from the LIVE `h` every outer step at the porous cadence and SUPERSEDING the porous write (the combined fraction already contains it).  `bt_H_ref` is untouched — a closed face removes transport capacity, not water.  **Fail-loud exclusions**: any coordinate but `z_fixed`; unresolved `z_fixed_h_ref`; empty `bt_H_ref`; `&ocean_wetdry_nml enable`; `&ocean_bt_nml bt_halo > 0`; GM / Redi / Fox-Kemper MLE; `&ocean_hvisc_nml nu_4 > 0` or `stress_tensor` (the free-slip closure covers the harmonic kernels only); and the three barotropic paths still on FULL-column weights — `&ocean_bt_nml correction_bc_pgf` (`compute_pbce`/`compute_gtot_faces`/the bc-PGF `du_bc` block: depth-mean-zero on the full column, not the open one), `substep_drag` (`compute_bt_rem`) and `wave_drag` (`compute_bt_rem_wave_drag`) (both damp on the full-column `Htot`).  Single-rank.  Default off ⇒ the masks stay at their `(1,1,1)` placeholder, no branch is taken ⇒ byte-identical | `test_ocean_zfixed_closed_faces`, `test_ocean_zfixed_bt_seiche`, `test_ocean_zfixed_cor_ref` |

**What it does and does not buy.**  It removes the catastrophic, GEOMETRIC part
of the `z_fixed` resting-state defect: `cavity_sloping_lid_rest_zfixed.nml`
goes from NaN during day 18 to completing 30 days at `En = 2.332E-06 m²/s²`
with budgets at `6E-13`; the ISOMIP+ ice-free control SATURATES at
`1.005E-08` by day 2, 56 000× below its knob-off twin; and the ISOMIP+ cavity
leg goes from NaN during day 3 to `8.842E-09` at day 2 (max|u| ≈ 3.7 mm/s
against the throwaway spike's 0.26 m/s) and runs on to day 4.  It does **not** give a resting state — the classical
partial-step PGF error at an OPEN face between a PARTIAL and a FULL cell
survives, and that is Yung et al. (2026) §3.2, a separate slice.  The bit-zero
`cavity_flat_lid_rest_zfixed` gate is preserved exactly with the knob on.

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

## Ice-shelf cavity geometry + load (static draft, barotropic datum, top-of-column pressure; ocean only)

A prescribed, time-constant ice draft `z_draft(i,j)` (m, positive **down**)
laid over the bed and absorbed into the **barotropic datum**, so a cavity
column starts from its loaded equilibrium rather than from `z = 0`; its
isostatic load is carried into the **pressure** through
`multilayer_state_t%p_top` and into the **barotropic mode** through the
datum — each exactly once, and never onto the `eta_forcing` seam.

| Piece | Form | Knob | Numerics / envelope |
|---|---|---|---|
| Static ice draft `z_draft` | `"none"` (identity) / `"flat"` (uniform inside the shelf box, 0 beyond the calving front `draft_x1`) / `"linear"` (`d0 + s·(x−x0)`, clipped at 0) / `"file"` (**fails loud** — the MPI-correct static-2-D reader is a later slice, and it is the only route to an ISOMIP+ draft, which has no analytic form) | `&ocean_cavity_dyn_nml enable` + `draft_config` | Filled on the host immediately after the bathymetry, FULL array including ghosts (formula evaluated at the ghost position), then carried through the bathymetry's own periodic/fold re-wrap + halo sequence. Box corners are METRES converted to GRID units (degrees on spherical/curvilinear) at the dispatch, like `&ocean_topo_nml slope_scale`; `±1e30` is the "no limit on this side" sentinel, which is what a shelf that reaches a wall needs. `draft_source="thickness"` converts an ice thickness by flotation, `ρ_ice·h/ρ₀`; `"in_situ"` isostasy fails loud. Tests: `test_ocean_cavity_draft` |
| **Datum** `bt_H_ref = b − z_draft` | the whole dynamical effect of a STATIC load | (implied by `enable`) | `bt_eta = Σh_layer − bt_H_ref` is then the deviation from the **loaded** equilibrium — zero at rest under the shelf — which is Losch (2008) §2.1's own convention. Every consumer of the water-column thickness `D = bt_H_ref + bt_eta` (BT continuity face thickness, Chapman phase speed, ALE `remap_h_ref`, the BT↔layer rescale) is then correct with NO cavity branch of its own. The geopotential stack is a separate datum and stays absolute (`e_face(1) = −b` from the true bed, so the column top lands at `−z_draft + η` by itself). Counted-once invariant: `ρ_ref·g·z_draft + (bt_H_ref − b)·ρ_ref·g ≡ 0`, asserted at configure. Gate: `test_ocean_cavity_equivalence` — a flat lid over a deep bed evolves like a shallower ocean |
| **Grounding** → LAND | `b − z_draft < h_min_cavity` ⇒ `wet_mask = 0` | `h_min_cavity` (default 10 m), `grounded_max_frac` | Routed through the SAME `seed_wet_mask_impl` the bathymetry uses, so the static metric-zeroing land mask and the finite land-state hold follow for free. **Never a thin film of water under grounded ice.** More than `grounded_max_frac` of the interior columns grounded ⇒ fail loud. No ice over land: the draft is forced to 0 wherever `b < LAND_DEPTH_THRESHOLD` (count logged), which keeps the counted-once invariant exact on every column |
| Isostatic load `p_ice_ref = (ρ_ref·GRAVITY)·z_draft` | Boussinesq-isostatic (flotation) load, ISOMIP+ §3.1.1's own `p = −ρ_sw·g·z_d` form | (implied by `enable`) | Formed as the SAME single product the FV-MOM6 surface BC forms, so `pa(nz+1) = ρ_ref·g·(−z_draft) + p_ice_ref` cancels to bit-zero at rest (to the rounding of one product under FMA contraction). Not the TRUE overburden of a stratified column (`g∫ρ̂`), which differs by `−g∫(ρ̂−ρ₀)` and leaves a depth-uniform `N²·z_draft·∇z_draft` residual in the raw PGF — chosen anyway because it is the load that makes the DISCRETE BAROTROPIC state exactly at rest, and the residual is annihilated by the split's depth-mean replacement. `draft_source="in_situ"` (the true solve) fails loud |
| **Load partition** `ms%p_top = p_ice_ref + sf%p_surf` | static ice load + atmospheric/anomaly load | consumers: `&ocean_pgf_nml p_top_in_bc`, `&ocean_psurf_nml in_eos` | Assembled on the host in `configure_ocean_cavity` (FINAL for a cavity without the psurf seam — the draft is static), and re-assembled once per outer step as an INLINE `do concurrent` in `ocean_dyn_step_split` when the psurf seam makes `sf%p_surf` live. Whole array, ghosts included ⇒ no halo of its own. `p_ice_ref` is **never** added to `sf%p_surf`: `eta_ib = −p_surf/(ρ₀·g_bt)` is built from that total and the datum already carries the same load, so only the ANOMALY reaches the `eta_forcing` seam — which for the isostatic default IS `sf%p_surf`. ⇒ an inverse-barometer run with no cavity is bit-identical, and a cavity with `p_surf = 0` sends the seam nothing (bit-identical to having no seam at all). Gates: `test_ocean_cavity_load`, `test_ocean_cavity_equivalence` |
| **`p_top_in_bc` REQUIRED for a varying draft** | fail-loud refusal, not auto-enable | `&ocean_pgf_nml p_top_in_bc` | A draft that VARIES must have its load in the `pa(nz+1)` BC: otherwise the stack sits ~5e6 Pa off its anomaly scale, the unsplit driver feels a raw `g·∇z_draft ≈ 2e-2 m/s²`, and `&ocean_bt_nml correction_h_weighted` redistributes the uncancelled uniform force as a real per-layer shear. Refused at configure (on the FILLED array) and in `validate_config` (on the namelist shape). A UNIFORM draft is EXEMPT — a load with no gradient is bit-identically inert in the top BC — which is what keeps the flat-lid equivalence gate expressible loaded and unloaded |
| **Rest residual under a sloping draft** | `PFu(k) − ⟨PFu⟩_h = N²·D³·(3σ_k²−1)/(12·dx·H̄)`, `D = s·dx` | — | The sigma-coordinate PGF truncation: second order in `dx`, CUBIC in the draft slope, independent of `nz`. Measured end to end (flat bed, flat isopycnals, `f = 0`, no forcing, 12000 s, `s = 1e-3`, `N² = 1e-5`): **`max\|u\| = 7.08e-8 m/s`** against the derived `a_peak·t = 1.04e-7`; with `N² = 0` it vanishes and the run sits at `2.6e-13 m/s`. This is the baseline for the sloping-coordinate PGF corrections (a later phase). Gates: `test_ocean_cavity_load` (unit) and `validation_examples/ocean/ice_shelf_cavity/` (end to end, 48x6x15 @ 2 km, 30 days, ZERO viscosity/drag/mixing): flat lid **`0.000E+00`**, uniform density **`1.34E-20`**, stratified sloping lid **`3.04E-08`** at day 30 after a 20-day plateau at the derived `1.0-2.1E-09`, then a ~3.2-day e-folding SATURATING at `3.89E-06` by day 60. dt-independent, only delayed by viscosity, absent without the lid slope or without stratification. Scoped XFAIL on `energy:rest-settles`; gates at `REST_1MM_S`, NOT the design's `REST_100UM_S` -- that is Phase 6's acceptance criterion (Yung et al. 2026, JAMES 18, e2025MS005645) |

**Envelope (every row fails loud at configure, naming the knob and the
reason):** `&ocean_pgf_nml form="fv_mom6"` and `gfs_scale = 1`; `vcoord_type`
in {`sigma`, `zstar`} — the two families that rescale the live column and so
follow the ice base for free; the refusal message now carries the OFFENDING
family's own reason, and they are not all "anchors at `z = 0`": `lagrangian`,
`zstar_sigma`, `rho` and `hycom` are geometrically datum-safe and are refused
for want of VALIDATION, while `z_fixed` hangs its stack from the column top,
`zstar_full` inverts which half of the column is resolved, `eulerian_z` drops
the free surface and `zsigma` is refused everywhere (see the vertical-coordinate
table below) — and `thickness_config /= "uniform_z"`; the SPLIT
solver (`n_inner ≥ 1`); single rank; `bt_halo = 0` (also in
`bt_halo_auto_exclusion`, so AUTO resolves to 0 instead of manufacturing a
width); and mutually exclusive with wet/dry, porous barriers, sea ice and
`&ocean_tides_nml use_sal`. `&ocean_zinit_nml enable` COMPOSES (its refusal
was lifted with P5.3 — the overlay measures depth from `z = 0`, not from the
column top). Default off ⇒ byte-identical.

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
| Top-of-column load in the FV-MOM6 pressure-stack surface BC (P5.0, ice-shelf-cavity prerequisite) — `pa(nz+1) = ρ_ref·g·η_geo + ms%p_top` instead of `ρ_ref·g·η_geo`. Reaches BOTH the PCM and the `reconstruct_for_pressure` branch (same Pass-1 seed, same face assembly). NOT a double count of the barotropic `eta_forcing` seam: a depth-uniform `p_top` perturbs EVERY layer's `PFu` by the same `−(1/ρ₀)∇p_top` (theorem in the `compute_fv_mom6_impl` docstring), and the split solver replaces the depth mean of the layer PGF with the barotropic solution (`F_bt_u_fast = F_bt_u − ⟨PFu⟩_h`, then `apply_bt_correction`), so that piece is annihilated — the two seams are orthogonal. On the UNSPLIT driver there is neither replacement nor seam, so the term is the load's only path into the momentum. What it buys where the load is large: `pa` is an anomaly stack about `ρ_ref·g·z`, so cancelling the load inside `pa(nz+1)` keeps the stack `O(1e4 Pa)` instead of `O(5e6 Pa)` and shrinks the `h_neglect` face-divisor leak by the same factor. `ms%p_top` has TWO producers that SUM — the `&ocean_cavity_dyn_nml` static ice load and the `&ocean_psurf_nml` seam (per-step refresh gate is the DISJUNCTION `in_eos .or. p_top_in_bc`, so neither consumer can read a stale value; a cavity without psurf needs no refresh, its configure seed being final); with neither producer the knob is inert and says so. REQUIRED, fail-loud, for a cavity whose draft varies. `compute_pbce` is deliberately load-blind (`pbce = ∂p_k/∂η`, and a static load has no `∂/∂η`). Orthogonal to `&ocean_psurf_nml in_eos` — that is the EOS pressure ARGUMENT, this is the PGF pressure BOUNDARY CONDITION | `p_top_in_bc=.true.` | `&ocean_pgf_nml`; FV_MOM6 only (fail-loud at configure otherwise — mont hard-zeroes `M(nz)`, fv_lite/fv_wright seed `p_edge(nz+1)=0`), default off ⇒ byte-identical. Test: `test_ocean_pgf_p_top_bc` |
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
| Z-sigma hybrid | `"zsigma"` — **REFUSED at configure**: the deep branch reads `z_ref_global` as metres while its only writer fills it with the dimensionless `k/nz`, so every z-level interval is `1/nz` m and the whole column collapses into the bed layer (with `Σ = H + η` still exact). Returns when the table is filled in metres — `test_ocean_vcoord_hygiene :: zsigma_is_refused`, `test_ocean_vcoord_interface_depths :: documents_zsigma_dimensionless_zref_collapse` |
| Z-star (lite) | `"zstar"` |
| Z-star/sigma hybrid | `"zstar_sigma"` |
| Z-star full (per-column) | `"zstar_full"` |
| Eulerian-Z | `"eulerian_z"` |
| Z-fixed (gprime) | `"z_fixed"` / `"gprime"` — nominal interfaces at `(nz-k)*h_nominal` below `z = 0` with `h_nominal = &ocean_topo_nml max_depth / nz_layers`; bed-side layers vanish to `zstar_h_min`. The ONLY z-like family legal under an ice-shelf cavity: there it reads `vcoord%z_top = metrics%z_draft`, vanishes the layers that outcrop into the ice to the same inert filler, and cuts the first live layer at the ice base into a partial top cell (min `0.1*h_nominal`, else the sliver merges into the layer below). `z_top = 0` ⇒ bit-identical. Under a cavity the top-side consumers read `multilayer_state_t%k_top(i,j)` (+ the `min`-rule face twins `k_top_u`/`k_top_v`), the first LIVE layer counting down, filled once at configure from the same `η = 0` target the closed-face mask uses and falling back to `nz` wherever nothing vanishes against the top: melt heat/salt/mass, top drag (and hence `stress_shelf` and the under-ice `u_*`), the implicit stress/drag folds and the momentum surface row, `mld_density`. Cavity envelope: melt and top drag are ACCEPTED; `kappa_h /= 0` is accepted iff `&vcoord_nml zfixed_closed_faces` (the face mask zeroes every flux touching a filler); KPP, EPBL, ideal age, GM/Redi/slopes, kappa-shear, tidal mixing and `regrid_time_scale > 0` still fail loud. |
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
| Barotropic transport renormalisation (`renormalise_zonal_flux_to_uhbt` / `renormalise_meridional_flux_to_vhbt`) | `src/core/ocean/kernels/continuity_ppm/rdb_continuity.F90` | — | Newton on a uniform `du` so `Σ_k uh_k = uhbt`, donor re-picked per iteration, CFL-bracketed. `&ocean_continuity_nml renorm_consistent_flux` (default off ⇒ byte-identical): a layer whose donor FLIPS carries `(u0+du)·h_face(new donor)` instead of the historical `flux0 + du·h_face(new donor)` — which jumps by `u0·(h_new−h_old)` at the flip and has NO root when `uhbt` lands in the gap (a wrong-sign layer transport at a thickness jump, e.g. sigma over a bathymetric step; the layer `η` then leaves the barotropic `η_end`) — and bisects when Newton leaves its bracket (MOM6 `zonal_flux_adjust`). Unflipped layers keep the historical expression, so a face with no flip is bit-identical either way | `test_continuity_multilayer` (`renorm_donor_flip_*`, `renorm_consistent_flux_no_flip_bit_identical`), `test_ocean_step_bt_mode` |

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

### Where α and β come from (E4)

`ρ_0`, `α_T` and `β_S` above are the LINEAR EOS's own coefficients.
Under `eos = "wright"` / `"roquet_spv"` they are not that EOS's
coefficients at all — the true α collapses toward zero at the freezing
point and grows strongly with pressure. `eos_buoyancy_coeffs` (and its
sign twin `eos_density_derivs`) evaluate `α = −∂ρ/∂T` and `β = +∂ρ/∂S`
from the ACTIVE EOS in closed form, per variant, as `pure elemental`
`!$acc routine seq` point routines. Measured (Wright, S = 34.5):

| state | α (kg/m³/°C) | vs the 10 °C value |
|---|---|---|
| −1.9 °C, 0 Pa | `2.619e-2` | ÷ 6.5 |
| −1.9 °C, 5e6 Pa (500 dbar) | `4.262e-2` | ÷ 4.0 |
| −1.9 °C, 1e7 Pa (1000 dbar) | `5.871e-2` | ÷ 2.9 |
| 10 °C (S = 35), 0 Pa | `1.711e-1` | — |

Which consumers take which:

| Consumer | Source of α/β | Pressure | Knob |
|---|---|---|---|
| KPP `B_0` (both passes of `vmix_kpp_overlay_impl`; also gates the non-local γ via `B_0 < 0`) | selectable | `ms%p_top` under `&ocean_psurf_nml in_eos`, else `eos%p_ref` — a SURFACE flux has no depth of its own | `&ocean_vmix_nml buoyancy_coeffs` |
| Double-diffusion density ratio `R_ρ = α·ΔT / β·ΔS` (`&ocean_ddiff_nml`) | selectable | true in-situ interface pressure, seeded from `ms%p_top` and accumulated `g·ρ₀·h` downward (EPBL's stack shape) | `&ocean_vmix_nml buoyancy_coeffs` |
| EPBL, kappa-shear, tidal mixing, isopycnal slopes, Redi | ALWAYS the active EOS (`eos_specvol_derivs`) | each builds its own in-situ stack; only EPBL seeds from `p_top` | — |
| PP81 N², convective-adjustment N², KPP bulk-Ri buoyancy, wave speed, MLE `b_ml`, `compute_pbce` | ALWAYS the active EOS — they difference `ms%rho_layer` itself | `eos%p_ref` (baked into `rho_layer`) | — |
| `ocean_cavity_const_t%alpha_T`/`beta_S` (hj99 / yung25 melt buoyancy) | its OWN pair, deliberately | n/a | none — ISOMIP+ Table 4 calibration constants, and **FRACTIONAL** (1/°C, 1/PSU), not the dimensional `eos_t` pair; left at their defaults by `configure_ocean_cavity` so ISOMIP+ comparability is not silently broken |
| Sea-ice coupler (`rdb_ice_ocean_coupler`) | reads no α/β and no density derivative at all | n/a | — |

| Source of α/β for KPP `B_0` + double diffusion | ocean | Knob | Test |
|---|---|---|---|
| `constant` — the scalar `&ocean_ic_nml alpha_T`/`beta_S` pair, whatever the active EOS | `default` ⇒ bit-identical | `&ocean_vmix_nml buoyancy_coeffs="constant"` | `test_ocean_eos_buoyancy` |
| `eos` — `eos_buoyancy_coeffs` from the ACTIVE EOS, per column (KPP) / per interface (ddiff). BYTE-IDENTICAL to `constant` under `eos="linear"` (that branch returns the handle members themselves, not `−ρ²·dSV/dX`), so the knob only bites under a nonlinear EOS | available | `&ocean_vmix_nml buoyancy_coeffs="eos"` | `test_ocean_eos_buoyancy` |

`validate_config` **warns** (does not refuse) on cavity melt × nonlinear
EOS × `buoyancy_coeffs="constant"`: ISOMIP+ prescribes the linear EOS, so
every shipped cavity namelist is unaffected and keeps running untouched.

**Known, out of scope, and NOT introduced here:** the KPP `B_0` line is
`g·(α·q_T − β·q_S)` with a DIMENSIONAL α, i.e. it is missing the `1/ρ₀`
that turns `∂ρ` into a buoyancy — EPBL's twin of the same quantity has it
(`b0 = g·ρ₀·(dSV/dT·q_T + dSV/dS·q_S)`, and `dSV/dT = α/ρ₀²`). The two
therefore disagree by a factor ρ₀ for the same physical flux. E4 changes
only WHERE α comes from, in the same units, so the discrepancy is
untouched and the default stays bit-identical; fixing it moves every KPP
answer and belongs in its own PR.

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
| Barotropic linear wave drag (Egbert & Ray 2001; Jayne & St Laurent 2001) — static per-face piston velocity `r_H` MULTIPLIED into `bt_rem_u/v` inside the BT substep; `form="uniform"` (global scalar) or `form="roughness_proxy"` (resolved-bathymetry-variance placeholder for the real subgrid `⟨h²⟩`, pending PR-14); composes with `substep_drag`.  Refused (with `substep_drag` and `correction_bc_pgf`) under `&vcoord_nml zfixed_closed_faces`: `Htot` is the FULL-column face depth.  `substep_drag` reads only `&ocean_bdrag_nml r`, so under `form /= "linear"` it is a no-op (`r = 0`) or a drag the slow step never applies — configure WARNING | barotropic | `&ocean_bt_nml wave_drag` | `test_ocean_wave_drag`, `test_ocean_zfixed_closed_faces` (refusals), `test_ocean_bt_substep_drag` (`warns_when_bdrag_not_linear`) |
| Atmospheric surface-pressure loading / inverse barometer (Wunsch & Stammer 1997) — `η_ib = −p_surf/(ρ₀·g_bt)` folded into the barotropic `eta_forcing` seam so momentum feels `−(1/ρ₀)∇p_surf` (~1 cm/hPa; a high depresses SSH); composes additively with the tide on the same seam; split-solver only, needs `enable_components` (reads `sf%p_surf`); ρ₀ from `eos%rho0` | barotropic | `&ocean_psurf_nml enable` | `test_ocean_p_surf`, `test_ocean_p_surf_bitident` |
| Reference density of the forcing terms — the `dt/(ρ₀·cp)` / `dt/ρ₀` divisors on EVERY surface heat/salt source (sea-ice coupling included), the wind-stress `τ/(ρ₀·h_top)`, the KPP `N²` / `u* = √(\|τ\|/ρ₀)` / kinematic fluxes behind `B_0`, and the geothermal `dt·Q_geo/(ρ₀·cp)` | not separately settable — all take `&ocean_ic_nml rho_0` | Fanned out from the ONE configured ρ₀ (`eos%rho0`) by `configure_ocean_reference_density`, joining the PGF and the twelve other slots that already took it. Before that wiring each slot kept a hard-coded 1035 that nothing assigned, so `rho_0 /= 1035` ran the EOS on the configured density and every forcing term on 1035, silently. `vmix%rho0` is the only one read on-device (inside `do concurrent`), and is correct because configure precedes `enter_data`. The implicit stress fold (`&ocean_vdiff_nml implicit_stress`) takes ρ₀ as an argument and now FAILS LOUD if it is omitted rather than defaulting to 1035. Default `rho_0 = 1035` ⇒ bit-identical. Test: `test_ocean_forcing_rho_ref` |
| Top-of-column pressure in the EOS's IN-SITU argument (E3, ice-shelf-cavity prerequisite) — the assembled `sf%p_surf` is copied once per outer step into `ms%p_top` (Pa, always allocated + zero-filled + device-mapped), and the ported in-situ builder evaluates at `p_top + hydrostatic` instead of starting at 0 Pa at the free surface. Ported consumers: the FV-Wright Picard column sweep (`form="fv_wright"`) and, since Phase 4b, the EPBL column stack (`&ocean_epbl_nml`, gated by the host scalar `epbl%in_eos` — the seed moves BOTH consumers of that stack, the in-situ `eos_specvol_derivs` argument and the PE weight `dmass·p_mid·dsv`, because they are the same pressure; a uniform load is gauge-neutral under a linear EOS, measured exactly zero, and moves the answer under Wright through `α(p)`/`β(p)` alone). With no ported consumer selected it is a documented no-op (warned at configure). It deliberately does NOT touch `ms%rho_layer` — a POTENTIAL density at the horizontally uniform `&ocean_eos_nml p_ref`, differenced along layers and vertically, so a per-column reference would fabricate an along-layer density gradient (N² therefore unchanged and self-consistent). EOS ARGUMENT only — this knob leaves the PGF top boundary condition alone; injecting the same `ms%p_top` there is the separate, orthogonal `&ocean_pgf_nml p_top_in_bc` (FV-MOM6 `pa(nz+1)`, see the Pressure Gradient section; FV-Wright's `p_edge(nz+1) = 0` is still untouched). `validate_config` REFUSES `in_eos` with the STILL-unported in-situ builders (kappa-shear, tidal mixing, Redi, isopycnal slopes, PGF in-layer reconstruction, sea ice). NOTE the `in_eos` gate is not redundant with `p_top` being zero: a cavity fills `p_top` with the ice load either way, so the gate is what keeps a cavity run bit-identical | EOS in-situ pressure | `&ocean_psurf_nml in_eos` (default off ⇒ byte-identical) | `test_ocean_eos_p_top`, `test_ocean_bl_under_ice` |
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

## Ice-shelf cavity basal melt (`&ocean_cavity_melt_nml`; ocean only, default off)

Kernel: `src/parameterizations/vertical/rdb_ocean_cavity_melt.F90`.
Coupling: `src/parameterizations/vertical/rdb_ocean_cavity_flux.F90` (the
`ocean_cavity_flux_t` slot, the far-field sampler and the per-thermo-step
driver), called from `engine_step_finalize` immediately BEFORE
`ocean_surface_flux_assemble`. Requires `&ocean_cavity_dyn_nml` (the draft,
the cover mask and the barotropic datum), `&ocean_eos_nml
tfreeze_set="isomip"` and `&ocean_forcing_nml enable_components`; every one of
those is a fail-loud configure requirement. Default off ⇒ bit-identical.

**v1 is VIRTUAL SALT, thermodynamics only.** The meltwater carries no mass, so
it adds no volume and no direct buoyancy (Phase 3); the ice-base MOMENTUM
sink is a separate opt-in (`&ocean_tdrag_nml`, above — off by default, and
sharing this group's `C_d` when on).  KPP and EPBL DO now feel the ice
(Phase 4b): their `u_*` is `√((stress_mag + stress_shelf)/ρ₀)`, where the
cover mask drives `stress_mag` to EXACTLY zero under a shelf and
`stress_shelf` carries the ice-ocean stress in its place.  Wind stress, surface
restoring, shortwave penetration and the uniform `&ocean_thermo_nml
q_heat/q_salt` are no longer refused: they are MASKED (see the cover-mask row
below).  `&ocean_psurf_nml` composes freely: the liquidus reads the assembled
`p_top = p_ice_ref + sf%p_surf`, so an atmospheric load under the shelf
depresses the freezing point with no extra wiring. See
`docs/CAPABILITIES_AND_LIMITATIONS.md`.

| Capability | ocean | Numerics | Primary test |
|---|---|---|---|
| Three-equation interface solve (liquidus + heat + salt) | `&ocean_cavity_melt_nml enable` | closed-form quadratic in `S_b` with the LARGER root (derived, no paper states it) and the cancellation-safe `q = −½(B + sign(B)√Δ)` form, because `Γ_S = Γ_T/35` puts production runs in the `\|B\| ≫ √(4AC)` corner; melt/freeze branch decided from `T*` BEFORE the solve, so it cannot disagree with its own answer; liquidus read off the `eos_t` handle (`&ocean_eos_nml tfreeze_set`), never duplicated | `test_ocean_cavity_melt` |
| Two-equation variant (`S_b = S_w`) | `&ocean_cavity_melt_nml enable` | the `γ_S → ∞` limit; IDENTICAL to the shipped sea-ice `rdb_ice_basal_flux` relax-to-`T_f` law at `γ_T·dt = h` (asserted, not asserted-to-be-similar) | `test_ocean_cavity_melt` |
| Exchange law `const_gamma` (`γ = Γ·u*`, Jenkins et al. 2010 / ISOMIP+) | `&ocean_cavity_melt_nml enable` | explicit in the melt rate; `Γ_T = 2.2e−2`, `Γ_S = Γ_T/35` are the ISOMIP+ STARTING GUESS, tuned per model (0.011–0.2 across the twelve ISOMIP+ submissions) — a knob, not a constant | `test_ocean_cavity_melt` |
| Exchange law `hj99` (Holland & Jenkins 1999 eqs. 14–18) | `&ocean_cavity_melt_nml enable` | turbulent + molecular sublayer with the McPhee `η*` stability parameter ⇒ IMPLICIT in the melt rate; outer bisection on `ln(L⁺)` (not Newton: the map kinks where the buoyancy flux changes sign); fail-loud at `f = 0` (the law has `\|f\|` in a logarithm and divides by it) | `test_ocean_cavity_melt` |
| Exchange law `yung25` (Yung et al. 2025 StratFeedback, eqs. 7–8) | `&ocean_cavity_melt_nml enable` | two power laws in `L⁺` capped at the Vreugdenhil & Taylor (2019) maxima ⇒ implicit, same bisection; its neutral limit is its OWN cap (0.012 / 3.9e−4), NOT the configured `Γ_T` | `test_ocean_cavity_melt` |
| Exchange laws `jenkins91` / `rosevear22` / `vt19` / `mk18` / `burchard22` / `jenkins21` | **RESERVED** (refused at configure) | enum values are nailed down and the dispatcher returns `CAVITY_MELT_NOT_IMPLEMENTED` — never a silent fallback to the default law. Each is implemented in the Python prototype `ice_shelf_melt/melt.py`, with its published ambiguities recorded there | `test_ocean_cavity_melt` (`reserved_laws_refuse`) |
| Ice conduction: insulating / H&J99 advective-diffusive | `&ocean_cavity_melt_nml enable` | insulating is the ISOMIP+ prescription (`κ_i = 0`); adv-diff collapses to `L_eff = L_f + c_i(T_b − T_ice)` via H&J99 eq. (31), zeroed on freezing as that paper prescribes, so the closed form survives and neither `H_I` nor `κ_I` is needed. The purely diffusive option is RESERVED (it moves the zero-melt point, so it changes the branch logic too) | `test_ocean_cavity_melt` |
| Friction velocity `u*` (tidal floor) | `&ocean_cavity_melt_nml enable` | `u* = max(√(C_d(u²+v²+u_tide²)), u*_min)` — ISOMIP+ eq. (27) + the Yung et al. (2025) floor; the tidal term belongs to the MELT `u*` only, never to the momentum drag | `test_ocean_cavity_melt` |
| **Far-field sampling** over `far_field_depth` METRES below the ice base | `&ocean_cavity_melt_nml far_field_depth` (10 m) | Thickness-weighted mean of `T`, `S` and the CELL-CENTRED velocity over the layers spanning that depth, with a PARTIAL last layer; vanished layers (`h ≤ H_VANISHED`) skipped, not clamped; `far_field_depth ≥ column` uses the whole column. **Metres, never "layer nz"** — the melt rate is roughly linear in the thermal driving it is handed, and the sampling distance is the dominant resolution artefact in the subject (Gwyther et al. 2020; Burchard et al. (2022) Table 2 p. 15, where the all-bulk error GROWS under refinement; Yung et al. (2026) p. 2074), so "layer nz" would make melt a function of the vertical coordinate's cell thickness — exactly what a coordinate study must hold fixed. Hold it FIXED across a coordinate sweep | `test_ocean_cavity_flux` |
| **Interface pressure** = `multilayer_state_t%p_top` | (implied by `enable`) | ONE pressure, THREE consumers: the FV-MOM6 pressure-stack surface BC, the in-situ EOS, and — with this slice — the melt liquidus. NOT `eos%p_ref` (a horizontally-uniform potential-density reference) and NOT the `eta_forcing` seam (a gradient). The melt path is a pure CONSUMER: the sole producer is `configure_ocean_cavity`'s `p_ice_ref + sf%p_surf` (P5.2), re-assembled per outer step when the psurf seam is live, so `configure_ocean_cavity_melt` writes nothing and instead ASSERTS `p_top >= p_ice_ref` on every cell (`cavity_count_unloaded_p_top`, exact because `p_surf >= 0`) — which catches a skipped or reordered producer rather than trusting the call order | `test_ocean_cavity_flux` (`ice_pump_deeper_melts_more`) |
| **Flux delivery** — the two OWNED components `heat_cavity` / `salt_cavity` | (implied by `enable`) | `heat_cavity = −q_ocean` (W/m², positive-down convention: `q_ocean > 0` warms the INTERFACE, so the ocean cools) and `salt_cavity = −m_mass·(S_far − s_ice)` (positive salinifies ⇒ melting freshens). NEVER `heat_added`/`salt_flux`, which the sea-ice coupler full-overwrites. The assembler folds both into `Q_heat`/`Q_salt` with a plain `+`; the filler latches `has_heat`/`has_salt` host-side and full-overwrites its own components each thermo step | `test_ocean_cavity_flux` (`flux_signs_cool_and_freshen`) |
| **Virtual salt flux** (`freshwater="virtual"`, the DEFAULT) | `&ocean_cavity_melt_nml freshwater` | `salt_cavity = −m_mass·(S_far − s_ice)` is the EXACT fixed-mass equivalent of adding mass `m_mass` at salinity `s_ice`: from `d(M·S)/dt = m·S_i` and `dM/dt = m`, `M·dS/dt = −m·(S − S_i)`. The heat twin of that dilution term, `−m·c_w·(T_w − T_b)`, is deliberately dropped — it is `m/(ρ_w·γ_T) ≈ 1e−3` of `q_ocean` because `T_w − T_b` is hundredths of a degree, whereas `S_far − s_ice ≈ 34 g/kg` makes the salt term the entire meltwater buoyancy signal. **Meltwater adds no VOLUME** in this mode — first-order for cavity circulation | `test_ocean_cavity_flux` (`virtual_salt_is_dilution_identity`) |
| **REAL freshwater MASS** (`freshwater="mass"`) | `&ocean_cavity_melt_nml freshwater` (default `"virtual"` ⇒ bit-identical) | `dh = m·dt/ρ₀` added to `h_layer(:,:,nz)`, `d(hS) = dh·s_ice`, `d(hT) = dh·T_b`, so the salinity falls by DILUTION and the closed form is `S(t) = S₀h₀/(h₀ + m t/ρ₀)`. **ρ₀, not the freshwater density**: the model is Boussinesq, its mass total IS `ρ·volume`, so a second density would put the tracked source and the tracked total on different scales. The heat reference is the model's implicit **0 °C** (`total_heat = Σ hT·area·ρ`), so in `hTr` space the enthalpy is simply `dh·T_b` — no heat capacity enters and the melt law's `c_w` never has to match the model's `c_p`. The virtual flux is NOT also applied to the tracer (that is the double count) but IS still assembled into `Q_salt`, because it is the entire meltwater buoyancy signal KPP/EPBL read for `B_0`; `ocean_cavity_mass_step` removes it again from salinity AND from the pseudo-salt mirror in the same stage, as the exact negation of the stamp. Passive tracers need nothing — `h` grows, `hC` does not. Applied IN-STAGE right after the surface-flux apply, so it spends the same `melt`, at the same stage weight, and `derive_bt_from_layers` puts the volume into `bt_eta` at the next stage with no barotropic forcing term. Refused with wet/dry and with `dt_tracer_advect_ratio > 1`, each naming its follow-up. Freezing that would drive the top layer through the vanish marker is clamped at `H_CAVITY_FLOOR = 2·H_VANISHED` — STRICTLY above it, because every thin-layer gate in the tree tests the marker with a strict `>`, so a layer left ON it reads as vanished and the next ALE remap's `c = hTr/h` guard returns 0, deleting the heat and salt the clamp was protecting — then counted, and FATAL.  A column already at or below the floor has the withdrawal REFUSED (`dh = 0`) rather than pinned up to it, so the clamp can never become a deposit | `test_ocean_cavity_freshwater` |
| **Sea-level compensation** (`volume_compensation="uniform_open_ocean"`) | `&ocean_cavity_melt_nml volume_compensation` (default `"none"`); requires `freshwater="mass"` | Each thermo step the domain-integrated melt volume is removed again, uniformly per unit area over the wet cells the ice does NOT cover, each parcel carrying that cell's own `T` and `S` (implemented as a top-layer RATIO `h_new/h_old` published in `comp_scale` and applied to EVERY tracer, so no concentration anywhere changes) and TRACKED as a sink in all three budgets. `"none"` is correct for a short run or an open boundary; for a CLOSED ISOMIP+ Ocean0 box it is metres per year of sea level (~30 m/yr over ~1e10 m² into ~4e10 m² of open surface). ISOMIP+ Sect. 3.1.3 leaves the restored Ocean0-2 uncompensated (their only seam is a restoring band that moves no volume) and allows compensation for the closed Ocean3/4. No open ocean at all ⇒ the withdrawal is exactly zero and the volume honestly stays in | `test_ocean_cavity_freshwater` |
| **Budgets** | (implied by `enable`) | Both components ride `Q_heat`/`Q_salt`, so they are integrated by the ordinary `ocean_surface_flux_apply_tracers` and land in the EXISTING `ms%heat_budget_surface` / `ms%salt_budget_surface` contributors — which the console already folds with the correct `ocean_budget_stage_weight`. No separate frazil-style accumulator and no full-weight/half-weight decision: the source is applied by the same per-stage kernel as every other surface flux, so it carries the same weight by construction. Under `freshwater="mass"` the same rule extends to MASS: the new host scalar `multilayer_state_t%mass_src` is accumulated with the SAME per-stage weight as `mass_out` and printed as the console's `src` on the `Mass` line, so `(M − M₀) + mass_out − mass_src` stays at round-off while the total grows; salt and heat still need no new accumulator because every real-mass increment is mirrored into the existing surface contributors | `test_ocean_cavity_flux` (`budget_closes_to_roundoff`), `test_ocean_cavity_freshwater`, `test_ocean_cavity_grounded_budget` |
| **Status accounting** | (implied by `enable`) | The per-column `CAVITY_MELT_*` plane is reduced ON DEVICE (`do concurrent … reduce`, never a host `count()` over a device-resident array) into four counters; `NONFINITE_INPUT`/`NONFINITE_STATE` on a covered column is **FATAL** (the safe state would hide an already-corrupt column behind a plausible run), while `NOT_CONVERGED` / `NO_PHYSICAL_ROOT` / `BAD_INPUT` / `NO_CORIOLIS` are COUNTED, warned, and given the kernel's zero-melt safe state. int64 running totals drain like `continuity_t%n_limited_total`. The fail-loud decision is a `pure` predicate (`cavity_melt_status_is_fatal`) so the suite tests it without provoking `error stop` | `test_ocean_cavity_flux` (`status_nonfinite_is_fatal`) |
| **Cover mask on ATMOSPHERIC forcing** | follows `&ocean_cavity_dyn_nml enable` (no knob of its own) | Under `cover_frac = 1` there is no atmosphere: every atmospheric term is scaled by `1 - cover_frac`, the cavity's own `heat_cavity`/`salt_cavity` are NOT.  Three sites, each chosen by what its consumers read.  (1) **The `tau` PAIR, at its source**, once from `configure_ocean_cavity` and again through `ocean_surface_stress_set_derived`'s optional `cover_frac` — masking the derived views would miss the implicit vdiff stress fold (`tau_u=ss%tau_x`) and the MLE front sampler, which read `tau` raw; the same call refreshes `stress_mag`, so KPP/EPBL `u*` is EXACTLY zero under cover.  FACE RULE: closed when EITHER neighbour is covered (`1 - max(cover_L, cover_R)`) — the only choice that leaves no wind on a covered cell; its price is a one-face front transition where the first OPEN cell keeps one live face and reads `\|tau\|/2`.  (2) **`ocean_surface_flux_assemble`** for `Q_heat_const`/`Q_salt_const`, `q_sw`/`q_lw`/`q_lat`/`q_sens`/`heat_added`, both `heat_content_mass*` and `salt_flux` — NOT at apply time, because `Q_heat`/`Q_salt` are what KPP/EPBL read for `B_0` and because the assembler is the last place the atmospheric and cavity bands are separable.  Components off ⇒ no assembler ⇒ the static scalar fill is masked once at configure (`ocean_surface_flux_apply_cover_const`).  (3) **SW penetration + surface restoring**, which route around `Q_*`, carry their own optional `cover_frac`.  Every path is an OPTIONAL argument (absent ⇒ byte-identical) and every mask is idempotent.  v1 cover is BINARY.  The one unreached path is `&ocean_dataovr_nml` file forcing (rewrites `tau`/`Q_heat` per bracket with no cover) — refused on the geometry group | `test_ocean_cavity_flux` (`cavity_cover_*`) |
| **Diagnostics** (13 derived-catalog entries) | `&ocean_cavity_melt_nml enable` (11) / `&ocean_cavity_dyn_nml enable` (2) | `melt` (kg m⁻² s⁻¹, + = melting), `melt_m_per_yr` (ISOMIP+ reporting unit: `/ρ_fw = 1000` × 365-day year, factor **31 536** exactly — Asay-Davis et al. 2016 §3.3; deliberately NOT the model ρ₀ = 1035, which is 3.5 % off every published figure), `thermal_driving` (`T_far − T_f(S_far, p_top)`), `haline_driving` (`S_far − S_b`), `tbdry`, `sbdry`, `tfreeze_ib` (the FAR FIELD's in-situ liquidus at the ice base — its difference from `tbdry` is the three-equation correction, which is why both ship), `exch_vel_t`/`exch_vel_s` (the CONVERGED γ_T/γ_S, carried on the slot by `cavity_melt_point_gamma`: under `hj99`/`yung25` they are implicit in the interface state, so a diagnostic that re-derived them from `u*` would report the NEUTRAL values), `ustar_shelf`, `cavity_melt_status`; plus the geometry pair `z_draft` and `water_column` (= `bt_H_ref + bt_eta`, the datum identity, gated on the geometry group alone). MISSING-VALUE CONVENTION: IEEE NaN outside the cover (masked by `cavity_flux%active` — covered AND wet AND the sample found mass — not by `cover_frac`), never zero, so a domain mean is a CAVITY mean; the console reports `missing=n/total` and the NetCDF variable carries `_FillValue`. A diagnostic requested without its prerequisite knob FAILS LOUD at configure rather than registering a plane of missing values | `test_ocean_cavity_diags` |

Notes: every solver returns a `CAVITY_MELT_*` status instead of `error stop` —
a `pure` `!$acc routine seq` procedure can neither log nor abort, and one bad
column must not take the run down. On any non-OK status the outputs are the
documented SAFE STATE (`m_mass` EXACTLY zero, `S_b = S_w`,
`T_b = T_f(S_w, p_b)`) so a failed column contributes nothing to a budget
rather than a plausible wrong number. Every clamp is `ieee_is_finite`-guarded
BEFORE the min/max, per the repo's NaN-laundering hazard.  The data-parallel
seam is `cavity_melt_columns` (a `do concurrent` over columns, inside the
module): on the GPU build nvlink cannot resolve a `!$acc routine seq` symbol
out of `librdb_core.so` into a `do concurrent` compiled in another translation
unit, so the coupling PR must extend that routine rather than write its own
column loop in the engine.

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
| `implicit_top_drag` | `.false.` | Ice-shelf top drag → vdiff surface (`k=nz`) diagonal, AND the wind-stress RHS on that row scaled by `(1 − cover)` so a covered face takes the drag and no wind. Requires `&ocean_tdrag_nml enable`; mutually exclusive with `&ocean_tdrag_nml implicit` and with `htbl>0`. Test: `test_ocean_top_drag`. |

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
