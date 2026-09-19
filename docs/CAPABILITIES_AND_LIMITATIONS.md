# Roundabout — Capabilities and Limitations

A snapshot of what the solver can and can't do today. Sourced from `ROADMAP_OCEAN.md`, `CLAUDE.md`, and the current state of the codebase. Every "can't" item below is something the codebase doesn't do yet — not a hard "won't"; consult the roadmap for project-scope tags and validation criteria for the planned items.

> **Scope note.** The `sim_type='coastal'` regime — the A-grid HLL/HLLC path, the unstructured triangular (KNP central-upwind) backend, the semi-implicit (Casulli) free surface, the non-hydrostatic extension, and the C/Python FFI — was split out into a **separate repository**. This tree is ocean-only, and nothing below describes those paths.

*Last reconciled against the codebase: 2026-08-23 (coastal carve-out).*

For the full feature checklist (planned, validated, in flight, deferred), see `ROADMAP_OCEAN.md`. For namelist / config keys + physics, see `REFERENCE.md`; for the per-closure shipped surface, `CLOSURE_MATRIX.md`. For acronyms (CWC, PGF, PP81, ALE, EPBL, MEKE, …) and concept definitions (vanishing layer, conservative soft clamp, mode split, …), see `codebase/CONCEPTS.md`.

**Regime status**: one path, `sim_type='ocean'` — Tier-1 operational as of
2026-05-24. Arakawa C-grid + continuity-PPM + split-explicit RK2 + Wright
EOS + Smagorinsky KH/AH + KPP + ALE remap. The MOM6-reference double-gyre
setup runs stable to day 580 on a single V100. Most physics slots are
wired; deferred work is the Hollingsworth-Källén Coriolis correction,
full KPP (V_t² + non-local), non-hydrostatic on the C-grid, and parts of
the MPI surface. See the [**Ocean path**](#ocean-path-sim_typeocean)
section for the detailed shipped surface — it is the authority for
everything below. Production-ready for regional hydrostatic ocean
configurations within those limits.

---

## Time-step constraint — read this first

> **The ocean path is not gravity-wave-CFL-bound.** `sim_type='ocean'` is **split-explicit**: the barotropic (gravity-wave) mode is sub-cycled in a fast forward-backward loop (`&ocean_bt_nml auto_n_inner` derives the substep count from the gravity-wave CFL each step), so the **outer/baroclinic `Δt` is bounded by advective/baroclinic CFL, not by `dx/√(gH)`** — O(minutes), e.g. `dt=1200 s` in the MOM6-reference double gyre.

The driver requires a fixed outer step (`&time_nml dt_fixed > 0`); there is no adaptive-CFL helper on this path.

**The outer split scheme is a real choice with a measured price on each side** (`&ocean_bt_nml split_scheme`, both shipped, both under test).

* **`"pred_corr"` (DEFAULT since 2026-09-14)** — the MOM6 predictor-corrector: `pc_be`-off-centred predictor, slow tendencies on the `u_av`/`h_av` step time-means, one prognostic update in the corrector, forward-backward gravity-wave pairing. Neutrally stable to `ω·dt = 2`: it removes the resting-state growth described below, lifts the internal-wave `dt` ceiling, and leaves the Eady benchmark's physical baroclinic mode intact and slightly stronger (max|v| ×2074 over 60 days vs ×1480). Its per-step cost is the same order as `ssp_rk2`'s — `coriolis_coast`, 48²×4, 20 simulated days on one V100: 70.3 s under pred_corr vs 65.9 s under ssp_rk2, about 7 % more. `validate_config` refuses it **fail-loud** outside its v1 envelope (`eulerian_z`, `&ocean_wetdry_nml enable`, `dt_tracer_advect_ratio > 1`); six shipped namelists pin `ssp_rk2` for one of those reasons. **Its own residual limit:** it does not ELIMINATE the resting-state growth, it slows it ~36× (83-day e-folding against ssp_rk2's 2.5-day one), so `resting_stratified_channel` stays a scoped XFAIL on the *settle* gate under the default — it passes the magnitude gate with margin.
* **`"ssp_rk2"` (EXPERIMENTAL)** — two identical stages + SSP average; widest envelope (every vcoord including `eulerian_z`, `&ocean_wetdry_nml enable`, `dt_tracer_advect_ratio > 1`) and the only scheme wired through windowed tracer advection. It remains **fully supported and fully tested** — the stability suite runs an `ssp_rk2` twin of every case whose namelist does not pin a scheme — and *experimental* is a label on the ANSWER, not a deprecation of the code path. **The defect, quantified:** the two-stage average amplifies an internal gravity wave by `√(1+(ω·dt)⁴/4)` per step, so it manufactures energy out of a motionless stratified state. On `validation_examples/ocean/eady/resting_stratified_channel.nml` — flat bed, periodic, stably stratified, at rest, ±0.5 mK seed, **no energy source of any kind** — En reaches **2.992E-05 m²/s² at day 25** (7.7 mm/s rms) and is still climbing on a 2.5-day e-folding; `pred_corr` on the identical file holds **1.739E-09** (**17 000×** less, 83-day e-folding). The outer time splitting is the cause and everything else was **exonerated** by direct substitution: the Coriolis form (`sadourny_hk`, `sadourny_energy`), the ALE remap (`vcoord = sigma`) and the PGF form (`fv_lite`) each moved the answer by < 0.1 %, `eady_dT_dz = 0` dropped En 119×, and REMOVING the lateral viscosity RAISED it — viscosity damps this mode, it does not cause it. **What it means in practice:** the growth is a `(ω·dt)⁴` noise floor, so what decides whether it matters is how hard `dt` is pushed against the internal-wave period. A forced, energetic, viscous configuration runs decades above that floor and never notices it; a quiescent, weakly-damped, or long-spin-up one does not, and there the manufactured energy IS the signal. Carried as a scoped XFAIL on `resting_stratified_channel__ssp_rk2`.
* **How the default came to move (2026-09-14).** Three things blocked the predictor-corrector and all three are closed. (1) the **land-masked-coastline NaN** — `coriolis_coast` went En 6.98e-07 (day 8) → 2.78e-04 (day 12) → 1.57e-01 (day 16) → NaN (day 18); the cause was the fast-loop Coriolis reference (`subtract_fast_cor_ref`) being evaluated on the stage-entry `u^n` while the slow Coriolis folded into `F_bt` was evaluated on `u_av`, so the uncancelled `f × (v̄_av − v̄^n)` forced every barotropic substep and pumped the basin's gravest Poincaré seiche. With `set_cor_ref_velocity` building the reference from the same velocity (MOM6 `ubt_Cor`), the case holds En 2.0e-07…1.1e-06 for 20 days and exits 0. (2) the **`nz = 1` GPU fault** — the vertical-friction velocity (and tracer) tridiagonals indexed 0 at a single layer; both now carry a single-layer path that puts the wind stress and the bottom drag in the one row. (3) the **seven GPU-only unit-test failures** that stood after those two — five non-finite within 1-5 steps (`periodic`, `obc_baroclinic`, `dyn_split`, `ice_restart`, `wetdry`) plus `restart` bit-exactness and uniform-`p_surf` gauge invariance, green on gfortran throughout. All seven had ONE cause, in the GPU data path rather than the scheme: `scratch_3d_buffer_t%enter_data` attached its payload with `!$acc enter data create`, which copies no host value, while the `pred_corr` PREDICTOR deliberately reads `hv%du_visc`/`dv_visc` without recomputing them (MOM6's `diffu(u[n-1])` reuse) — and at step 1 there is no previous producer, so the host read the zero `init` promised and the device read the allocator's leftovers. `enter_data` now device-zeroes the payload, gated by `rdb_test_scratch_3d_device`, and the suite is **187/187 on both toolchains**. The flip re-baselined every shipped answer and every golden for the ~54 of 64 namelists that do not pin a scheme; the six that pin `ssp_rk2` and the four that pin `pred_corr` are unchanged.

Several individual operators are additionally unconditionally stable on their sub-problem: vertical viscosity/diffusion is a backward-Euler tridiagonal per column, bottom drag and wind stress can be folded into that same solve (`&ocean_vdiff_nml implicit_drag` / `implicit_stress`), and Coriolis is the closed-form rotation.

This is a **sub-cycled fast mode, not a semi-implicit free surface**. Roundabout does not treat the gravity-wave term implicitly in the SCHISM/Casulli sense; the barotropic substeps resolve it explicitly, just at their own smaller `Δt`.

---

## Solver core

**Structured Arakawa C-grid, hydrostatic, Boussinesq.** Layer thickness is prognostic (continuity-PPM transport, no Poisson constraint); momentum is vector-invariant with a Sadourny PV-flux Coriolis-advection operator and a finite-volume pressure gradient; the vertical coordinate is ALE (advance Lagrangian, remap conservatively). Horizontal grids: `cartesian`, `spherical`, `supergrid` (MOM6 mosaic), `tripolar` (Murray 1996 cap + single-rank north fold).

The full operator-by-operator surface, with knobs and limits, is in the [Ocean path](#ocean-path-sim_typeocean) section below.

---

## Physics gaps (cross-cutting)

> Everything **shipped** is enumerated in the [Ocean path](#ocean-path-sim_typeocean) section. This list is the standing gaps.

- **GLS, k-ω / MY2.5** turbulence closures. Shipped instead: PP81 interior, KPP boundary layer (default on), EPBL (`&ocean_epbl_nml`), kappa-shear (`&ocean_kappa_shear_nml`), tidal mixing, convective adjustment, double diffusion, and two mutually exclusive background schemes.
- **TEOS-10 GSW in-situ EOS branch.** Wright (1997) and a Roquet et al. (2015) polynomial ship; a full GSW branch does not.
- **Sediment transport, biogeochemistry, vegetation drag.** No source terms beyond the surface/bottom flux set.
- **Online bulk aerodynamic flux formulae** (Zeng / COARE). Forcing must be pre-computed to model-grid `τ`, `Q_net`, `E−P` — there is no in-model conversion from `(U10, T_air, q_air, SST, SLP)`.
- **Sea ice** is a dynamical core + column model, not a complete ice model — no ridging, melt ponds or lateral melt; see the [Sea ice](#sea-ice) subsection and "What Roundabout is not".

---

## Forcing (ocean path)

### Shipped
- Wind stress (programmatic / analytic profiles; `&ocean_surface_*`).
- **Time-varying NetCDF surface forcing** (`&ocean_dataovr_nml`, default off) — wind stress, heat, salt, evaporation and liquid precipitation bound to pre-regridded model-grid files through the shared PR-14 reader. See the full entry in the ocean-physics list below for the tag table, the `enable_components` destination rule, and the v1 limits (shared time mode across tags, single-rank).
- Surface heat + salt flux, shortwave penetration, surface buoyancy restoring.
- Equilibrium body-force tide + scalar SAL + boundary-tide nodal correction.
- **Atmospheric surface-pressure loading / inverse barometer** (`&ocean_psurf_nml enable`, default off ⇒ byte-identical). `η_ib = −p_surf/(ρ₀·g_bt)` is folded into the barotropic `eta_forcing` seam, so the momentum feels `−(1/ρ₀)∇p_surf` — an atmospheric high depresses SSH ~1 cm/hPa. Composes additively with the equilibrium tide + scalar SAL on the same seam; ρ₀ from `eos%rho0`. **v1 fill is a uniform scalar** (`p_surf_const`, seeded into `sf%p_surf_atm`) which is provably inert (only ∇p_surf is physical) — a genuine load needs `p_surf` wired up as an `&ocean_dataovr_nml` tag (the reader ships; `p_surf` is not yet one of its tags). Requires `&ocean_forcing_nml enable_components=.true.` (reads `sf%p_surf`), the split solver (`n_inner ≥ 1`), and is mutually exclusive with `&ocean_bt_nml bt_halo > 0` — an explicit width fails loud at configure, and the `bt_halo` AUTO default (the multi-rank `-1` sentinel) resolves to 0 under psurf rather than to the wide-halo march-in.

### Not yet shipped
- **File-driven `p_surf`** (reanalysis MSLP) — the NetCDF forcing reader itself
  ships (`&ocean_dataovr_nml`, above); `p_surf` is simply not one of its six
  tags (`tau_x`, `tau_y`, `heat`, `salt`, `evap`, `lprec`). Adding it is a
  `register_tag` entry plus a destination slot, not new reader machinery.
- **Sea-ice mass loading into `p_surf`** — deferred to the ice-loading PR (the `eta_ib` seam and the `sf%p_surf` overwrite convention are in place for it).
- **`p_surf` in the EOS pressure argument** / ice-shelf-cavity surface-pressure curvature corrections / `MAX_P_SURF` load cap — out of scope.

---

## I/O

### Shipped
- Per-rank NetCDF output (CF-1.8) driven by the diagnostics manager
  (`&ocean_diag_nml`), offline merge via `tools/merge_output.py`
- Diagnostics registry + cadence dispatch, per-diag output vcoord remap
  (z / sigma / density bins), MEAN/MAX/MIN/INTEGRAL time-ops
- Optional deflate compression
- Single-precision **diagnostic** output (`&ocean_diag_nml output_precision =
  "single"`, default `"double"` ⇒ byte-identical) — halves the bytes per frame
  on the write-bandwidth-bound diag stream.  Deliberately scoped to that
  one stream: restarts and the console conservation totals / checksums are
  unconditionally working precision, so a lossy restart is not requestable.
  Opting in means regenerating any regression baseline that byte-compares
  diagnostic NetCDF.
- Restart / warm start, bit-exact round-trip (including the wet/dry
  hysteresis registry)

### Not yet shipped
- MPI I/O server hand-off for the diag manager — `&output_nml use_io_server`
  is accepted but warns and is ignored on this path; the emit is serial
  per-rank NetCDF.
- Gauge / station time-series output
- Lagrangian drifters
- Parallel collective NetCDF
- In-situ visualization

---

## Performance and scaling

### Shipped
- Single source `do concurrent` + OpenACC, runs on NVHPC GPU/multicore, gfortran, ifx (OpenMP-target build also exercised).
- CUDA-aware MPI for halo exchange (`RDB_CUDA_AWARE_MPI`).
- Per-slot up-front memory reporting (`rdb_mem_report`) reconciled against the measured device mapping, so a gated-off closure provably costs nothing.
- Gated allocation for default-off closures (EPBL, kappa-shear, tidal mixing, GM, Redi, MLE, VarMix, MEKE, isopycnal slopes) — a plain run pays none of their multi-GB footprint.

### Headline numbers (single V100)

| Configuration | Wall time |
|---|---|
| Tasman 2 km (eddy-resolving, KPP + ALE + sponges) | ~34 s / simulated day |
| Double-gyre MOM6 ref (44 × 40 × 2, dt=1200) | ~16 s / 30 simulated days |

### Not yet shipped
- Mixed-precision tracers (fp32/bf16)
- Dynamic load balancing for shifting wet fractions
- Topology-aware (NVLink/NVSwitch) MPI tuning
- Fault tolerance (ULFM, mid-run checkpoint-restart on rank failure)
- GPU CI runner (GPU tests run locally today)

---

## Validation

Analytical / exact-solution tests shipped: vertical diffusion vs erfc, gravity-wave phase speed, geostrophic adjustment, lake-at-rest over a seamount, quiescent island-at-rest, Thacker moving shoreline (wet/dry), Eady baroclinic growth rate, Nansen free drift and Stefan melt (sea ice). ~177 `rdb_*` ctest entries over 180 test sources — 162 labelled `ocean`, 18 labelled `core` (`ctest -L ocean` / `-L core`); regime labels are mandatory per row and enforced by a pre-commit hook. MPI consistency tests for the C-grid halo (`tests/mpi/`). CPU↔GPU bit-comparison via the gcc + nvhpc CI matrix. A golden-output regression harness lives in `tests/regression/`.

Canonical benchmark configs live under `validation_examples/ocean/` — `seamount/`, `geostrophic_adjustment/`, `eady/`, `eddy_test/`, `double_gyre/`, `acc_channel/`, `neverworld2/`, `baroclinic_channel/`, `island_at_rest/`, `flow_past_island/`, `double_drake/`, `sea_ice/`, `tides/`, `sponge_demo/`, `ideal_age/`, `epbl_mld/`, `data_forcing/`, and others. Each has a README with expected behaviour, analytical scales (when applicable), and diagnostic interpretation of failure modes. `validation_examples/test_cases/tasman_validation/` carries the regional real-bathymetry case.

Not shipped: formalised observational comparison for a real domain, performance regression CI.

---

## Data assimilation

- **Shipped:** nothing beyond the namelist/restart surface and the sponge/nudging boundary machinery (`&ocean_sponge_nml`, per-edge asymmetric nudging), which is enough for relaxation toward a reference state.
- **Not shipped:** in-model observation operators, ensemble support, 4D-Var adjoint. The in-process state-setter API went with the FFI carve-out.

---

## Machine learning

Nothing built today, and no in-process API surface for a surrogate to plug into since the C/Python FFI was carved out with the coastal repo. A surrogate swap-in would need that seam rebuilt first.

---

## Ocean path (`sim_type='ocean'`)

The dynamical core of this tree, designed for regional /
basin-scale hydrostatic ocean. Operational as of 2026-05-24; the
MOM6-reference double-gyre setup runs stable to day 580 on a single
V100. Specific shipped surface below.

### Architecture

Arakawa C-grid (`u_face_x`, `v_face_y`, `eta` centred), split-explicit
RK2:

1. Slow tendencies — Coriolis advection, PGF, hvisc, vmix (implicit
   tridiagonal), bottom drag, surface stress, tracer hdiff/vdiff,
   vertical advection.
2. Barotropic fast loop — N forward-backward-Euler substeps on
   `(η, u_bt, v_bt)` + nonlinear `ζ + KE` terms; `auto_n_inner=.true.`
   derives N from the gravity-wave CFL each step.
3. BT correction — distributes the barotropic Δu back into the
   layers; optionally h-weighted (`&ocean_bt_nml correction_h_weighted`)
   and further biased by the vdiff-produced viscous remnant γ
   (`correction_visc_rem`, requires `correction_h_weighted`) — γ is
   the momentum tridiagonal's own sensitivity to a uniform barotropic
   acceleration, so it biases the corrector's Δu against layers inside
   a frictional bottom boundary layer.  γ ≡ 1 (the joint weight
   reduces to plain h-weighting) unless `&ocean_vdiff_nml
   implicit_drag = .true.` also folds bottom drag into the vdiff
   operator.
4. Stage 2 averages.

Continuity is a transport equation (`∂h/∂t = -∇·(hu)`) solved with
**continuity-PPM** — no Poisson constraint, no FFT projection.

### Shipped — dyn-core operators

- **PGF**: Montgomery potential (`form="mont"`, `OPGF_VARIANT_MONT` —
  **the default** — the Boussinesq `M = p/ρ0 + (g·ρ/ρ0)·z` recursion plus
  one horizontal difference; general-purpose, valid over sloping
  bathymetry and every vcoord, and algebraically `fv_lite` on aligned
  columns), z-corrected FV (`fv_lite`), FV + Wright in-situ Picard re-eval
  (`fv_wright`), reduced-gravity 2-layer (`gprime`).  `gprime` is
  **NK=2-only** and now configure-enforced (PR-6): `form="gprime"` with
  `nz_layers /= 2` is a fail-loud abort, not a silent zero-PGF on the
  extra layers — use the default `form="mont"` for a general-nz PGF.
  `mont` and `fv_lite` carry the SAME physics content (layer-mean ρ, no
  in-layer quadrature) and coincide algebraically wherever the two columns
  meeting at a face have equal layer thicknesses.  They are genuinely
  different discretisations where thicknesses are UNEQUAL — σ / z*σ over a
  slope — and NEITHER is exact there; that residual is the open
  sigma-coordinate PGF error, not a property of one form.
- **EOS**: Wright (1997) nonlinear rational fit (production); linear
  EOS available.
- **Coriolis**: Sadourny PV-flux (`sadourny`); Sadourny +
  Hollingsworth-Källén guard (`sadourny_hk`, enums live, kernel
  refinement deferred).
- **Lateral closures**: Leith (vorticity-gradient ν_h per face);
  Smagorinsky_KH + Smagorinsky_AH (flow-aware biharmonic); constant
  `nu_h` / `nu_4` floors.
- **Barotropic linear (Rayleigh) wave drag** (`&ocean_bt_nml wave_drag`,
  Egbert & Ray 2001; Jayne & St Laurent 2001): the bulk energy sink for
  the barotropic tide — a static per-face piston velocity `r_H(x,y)`
  MULTIPLIED into the BT-substep damping factor `bt_rem_u/v`, composing
  with `substep_drag`. Default off ⇒ bit-identical. `form="uniform"` or
  `form="roughness_proxy"` (a documented resolved-bathymetry-variance
  placeholder for the real subgrid `⟨h²⟩`); `form="file"` fails loud at
  configure pending PR-14's NetCDF map reader.
- **Fox-Kemper mixed-layer-eddy restratification** (FK08/FK11, B5,
  `&ocean_foxkemper_nml enable`): submesoscale ML eddies slump lateral
  buoyancy fronts via an overturning streamfunction
  `Ψ = Ce·(H_ml²/|f|)·∇b̄·μ(z)`.  Injects ML-confined per-layer mass
  transports into the continuity mass fluxes BEFORE the divergence (never
  touches velocities) ⇒ conservative by construction (`Σ_k a(k)=0`, a
  closed overturning cell — verified to round-off).  `H_ml` from
  `epbl%mld` (EPBL is the enabler; requires `ocean_epbl_nml enable`).
  Surface transport is westward toward the dense column (restratifying).
  Two timescale forms: bare `Ce/max(|f|,f_floor)` (analytic-gate default)
  and the FK11 momentum-mixrate form (`use_mom_mixrate`, **production-
  recommended** — suppresses restratification under vigorous mixing).
  Runs at thermo cadence, once per outer step; bandwidth-bound 2D + per-
  layer μ scatter.  Default off ⇒ bit-identical.  Deferred:
  `resolution_taper` (B2 res_fn double-counting hook, hard error until
  B2), slow-filtered-MLD second transport (instantaneous MLD only),
  μ cubic tail (`tail_dh`).  Tests: `test_ocean_foxkemper`.
- **Vertical mixing**: PP81 (Richardson) interior + KPP Phase 1
  (shear-driven bulk-Ri sweep + shape function `G(σ) = σ(1-σ)²`) +
  KV_ML_INVZ2 (MOM6 inverse-z² surface band over `HMIX_FIXED`) +
  DT_THERM cadence skip + HARMONIC_VISC face-thickness option.
- **Convective adjustment** (Brunt-Väisälä trigger, CVMix
  `CVMix_convection`-style, `&ocean_conv_nml enable`, default off ⇒
  bit-identical): where the interior `N² < n2_thresh` (dense-over-light),
  raises `kt -> max(kt, kd_conv)` / `kv -> max(kv, prandtl_conv·kd_conv)`,
  strictly below the active KPP/EPBL boundary-layer depth (the BL owns
  its own convective response — `w_*` / convective TKE).  `kd_conv`
  defaults to 1.0 m²/s, ~100× PP81's own `Ri<0` ceiling (~1.01e-2
  m²/s), admissible only because `vdiff_apply_tracers` /
  `vdiff_apply_momentum` are unconditionally-stable backward-Euler
  solves.  A CONTRIBUTOR (`max()` floor, not an additive increment):
  feeds `vmix_assemble` like every other interior closure.  Two v1
  limits vs. MOM6's `MOM_CVMix_conv`: **D1** the N² trigger uses a
  single global `p_ref` potential density (`ms%rho_layer`), not a
  locally-referenced interface-pressure density, so thermobaricity is
  neglected (same assumption PP81 already makes); **D2** uses `max()`
  rather than MOM6's additive `Kd += kd_col` (differ by at most the
  resolved interior Kd, ≤1% of `kd_conv` — physically immaterial, and
  `max()` is exact/idempotent under the every-RK2-stage call cadence).
  Writes `kv`/`kt` only, never `ks` — `ks` is not rewritten by any
  per-stage contributor today, so a `max()` floor on it would ratchet
  monotonically and never relax; salt convects for free once a future
  `ks <- kt` split runs downstream of this call (until then salt
  genuinely does not convect, a pre-existing limitation this closure
  does not worsen).  Tests: `test_ocean_convection`.
- **EPBL** (Reichl & Hallberg 2018 energetics-based PBL,
  `&ocean_epbl_nml enable`): prognostic TKE budget per column —
  wind (`mstar·ρ₀u*³dt`) + convective (`nstar`-weighted) energy
  spent interface-by-interface via a closed-form energy solve;
  the MLD is where the energy runs out (false-position root-find,
  previous-step seed).  mstar schemes `constant`/`om4`/`rh18`;
  gravity-wave column-height correction; per-column TKE-budget
  diagnostics that close to round-off.  Replaces the KPP overlay
  when on (mutually exclusive); combines with PP81 `add`/`max`;
  `Kv = prandtl·Kd`.  Runs at thermo cadence; scalar-carry GPU
  sweep (no per-column work arrays).  **Langmuir enhancement
  shipped** (`use_lt`): LF17 wind-only statistical waves — COARE 3.5
  u*→U10 + Phillips-spectrum surface-layer-averaged Stokes drift,
  one scalar La per column, Reichl & Li (2019) mstar enhancement
  (rescale/additive) with the Li et al. (2016) stability-modified
  La; zero wave-model coupling, zero Stokes arrays.  Design + MOM6
  knob mapping: `docs/generated_nml_knobs.md`.
  **Penetrating shortwave now charged** (PR-21, `&ocean_thermo_nml
  epbl_sw_ctke`, default on): the per-layer `cTKE` ledger pays the
  in-layer PE cost of the exponentially-distributed SW absorption
  (`Phi(h/zeta)·`skin-cost, MOM6 `absorbRemainingSW` shape), filled in
  the prep sweep and drained in the interface sweep; charging the truth
  (rather than dumping all SW at the skin) DEEPENS the midday MLD.  Inert
  at `sw_pen_frac=0` ⇒ bit-identical.
- **kappa-shear** (Jackson, Hallberg & Legg 2008 prognostic shear
  turbulence, `&ocean_kappa_shear_nml enable`): INTERIOR closure —
  per-column coupled (κ, TKE) steady-state equations solved by Picard
  iteration with adaptive time-substepping; the diffusivity diffuses
  vertically with a stratification/rotation/boundary-limited decay
  length, so resolved shear layers entrain correctly at marginal Ri and
  the mixing is self-limiting (relaxes the column toward Ri ≳ Ri_c).
  Coexists with KPP/EPBL (additive: `kt += κ`, `kv += prandtl_turb·κ`
  every stage; column solve at thermo cadence).  Tracer-point path,
  Picard-only, Boussinesq; no massless-layer merging yet (vanishing
  ZSTAR_FULL layers are floored in the gather).  GPU: one
  `do concurrent` column kernel, all-local fixed-size arrays (the
  MRE-measured L1 layout; 255-register/12.5%-occupancy bound).  Kernel
  matches the validated single-column prototype to ~1e-13.  Knob mapping:
  `docs/generated_nml_knobs.md`.
- **First-baroclinic wave speed + Rossby radius** (B1, Chelton et al.
  1998, `&ocean_wavespeed_nml enable`; **diagnostic, default off**):
  per-column `cg1` (m/s) + first-mode Rossby deformation radius `Rd`
  (m) + `Rd/dx`, called once per outer step (thermo-cadence-gated,
  further gated by `n_wavespeed`) in `ocean_dyn_step_split`, BEFORE
  both `varmix_compute` call sites and `run_meke_step` — feeds the B2
  GM/Redi/MEKE resolution-aware scaling (`&ocean_hvisc_nml
  resoln_scaled_visc`, `&ocean_varmix_nml resoln_scaled_khth/khtr`,
  MEKE's `Ldeform`).  `ocean_dyn_step` (the unsplit path) does not call
  it.  Rigid-lid Sturm–Liouville eigensolve discretised
  from `rho_layer` (no T/S re-eval): per-interface reduced gravity with
  a static-instability floor, a backtracking convective layer-merge
  preconditioner (welds inverted/degenerate layers so no zero-`gprime`
  row survives the active range — closes the partial-inversion blow-up),
  then a fixed-budget (40-iteration, no early-exit → warp-divergence-free)
  Sturm-count bisection for the largest `c²`.  Equatorial `Rd` uses the
  smooth `cg1/sqrt(f² + 2β·cg1)` blend, where `f` (`f_centre`) and `β`
  (`beta_centre = |∇f_centre|`) are both static fields filled at
  configure from `metrics_fill_coriolis` (planetary dispatch on
  spherical/tripolar, bit-identical beta-plane elsewhere) — not a
  hard-coded beta-plane / namelist-scalar `β`.  `rd_over_dx = Rd /
  metrics%dxT` (metres — NOT `grid%dx`, which is degrees on
  spherical/supergrid/tripolar; bit-identical to the old expression on
  Cartesian, where `dxT ≡ grid%dx`).  GPU: outer-shim + explicit-shape
  flat-impl kernel, 5 fixed-size column arrays, divergence-free
  fixed-budget bisection (no data-dependent early-exit) — occupancy
  comfortably above the kappa-shear cliff per the dev-time MRE.
  Matches the verified Python reference prototype to ~1e-6; golden
  nz=4 column `cg1 = 3.266407 m/s`.  Deferred: the equivalent-barotropic
  (`use_ebt`) and N²-monotonising (`mono_n2`) EBT refinements (knobs
  reserved, default off); the `|∇|f||`-vs-`|∇f|` equator kink shared
  with VarMix/MEKE (house idiom, not fixed here).
- **Bottom drag**: `linear` (Rayleigh rate, 1/s) and `quadratic`
  (log-layer, MOM6/ROMS default `Cd ≈ 2.5e-3`). Both have an
  HBBL-distributed mode that spreads the stress across the bottom
  `hbbl` metres rather than the bed-most layer alone.
- **Time-varying NetCDF input reader** (PR-14, `rdb_ocean_data_input`):
  the shared, target-agnostic reader every forced-hindcast/regional-
  nesting capability builds on. A consumer registers a `(file,
  variable, destination)` triple via `ocean_data_input_register_2d/_3d`
  (or the OBC-segment variants) and gets back an opaque `id`; the
  driver's one-line `ocean_data_input_update_all` hook refreshes every
  registered field's time bracket each outer step and blends linearly
  in time ON-DEVICE into the caller's own, whole, mapped array —
  `linear`/`cyclic` (explicit-period climatology)/`static` time modes,
  out-of-range default abort (opt-in clamp per field),
  `fill_static_host{,_3d}` for setup-time static fills before a
  destination is device-mapped. Files must be pre-regridded to the
  model horizontal grid (no in-core horizontal interpolation, no
  vertical remap of a source z axis, no calendar). **Still missing**
  are file-backed OBC segment data and NetCDF climatology restoring,
  both of which this reader unblocks but does not itself deliver.
- **File-backed surface forcing** (PR-15, `&ocean_dataovr_nml`,
  `rdb_ocean_data_forcing`) — the first consumer of the reader above,
  and what makes an atmospherically forced hindcast possible. Flat
  per-tag knobs (`<tag>_file`, `<tag>_var`, `<tag>_scale`,
  `<tag>_add`) bind a time-varying NetCDF field onto a forcing slot;
  a blank `<tag>_file` leaves that slot on its configure-time
  scalar/formula value, and `enable = .false.` (default) registers
  nothing ⇒ bit-identical. Tags: `tau_x`/`tau_y` (→ the C-grid face
  stress, with `stress_mag` re-derived after every blend so KPP/EPBL
  see the current wind), `heat`, `salt`, `evap`, `lprec`. The heat and
  salt destinations depend on `&ocean_forcing_nml enable_components`:
  off ⇒ `Q_heat`/`Q_salt` direct, on ⇒ the `heat_added`/`salt_flux`
  components (because `ocean_surface_flux_assemble` rebuilds
  `Q_heat`/`Q_salt` from the component set every thermo step). `evap`
  and `lprec` exist only in the component set and fail loud without
  it. Time-axis mode (`linear`/`cyclic`/`static`), `cycle_period`,
  `t_offset` and the out-of-range policy are **shared by every tag** in
  v1 — a per-tag time mode (interannual winds alongside a cyclic SST
  climatology) is the known follow-up. Stress ghost cells are filled by halo
  exchange + periodic wrap + tripolar fold, never by extrapolation, and
  the same refresh runs for the analytic `wind_config` seeds — so the
  wind seam is correct under decomposition (the flux tags need no ghost
  fill: they are applied column-locally).  The stress file is C-grid
  staggered and must carry `nx_phys+1` x-face values / `ny_phys+1`
  y-face values. Online bulk formulae are out of scope: v1 consumes
  offline-preprocessed, model-grid fluxes. The sea-ice `tau_a`
  atmospheric-stress snapshot is still taken once at configure, so it
  does **not** follow a time-varying wind.
- **Dynamic wetting/drying** (`&ocean_wetdry_nml enable`, default off ⇒
  byte-identical; spec + prototype numbers in `docs/ocean_wetdry_plan.md`):
  per-substep hysteresis cell wet mask (`dry_depth`/`rewet_depth`) +
  UPWIND BT face thickness + per-cell positive-definite outflow limiter
  (`θ = min(1, avail/outflow)`, faces scaled by `min(θ_L, θ_R)` ⇒ total
  depth `D ≥ 0` unconditionally, NO thin-film mass injection — conserves
  to round-off) + bed-blocking momentum gate (a face into a dry cell is
  a wall unless the wet surface overtops the dry bed) + `FROUDE_CAP`
  thin-face runaway-velocity guard.  Composes multiplicatively on top of
  the static land masks; layer velocities at blocked faces reset via
  `mask_layer_velocities`; surface heat/salt flux masked on dry columns.
  Envelope (all fail-loud at configure): sigma / zstar-lite vcoord,
  single-rank, split solver + `ppm_limit_pos` required; mutually
  exclusive with BT_cont / upstream-h / sw_pen / restoring / geothermal.
  Thacker moving-shoreline gates: period < 2%, shoreline < 4 cells
  (measured −0.9% / 3.6 cells), ~12%/period front dissipation on the
  frictionless runup (documented limit).  **v2 intertidal seed**
  (`land_margin`, default 5 m, consulted only when `enable=.true.`):
  the static land seed becomes `b < -land_margin` (bed above the
  highest credible water level) instead of the 2 m rest-depth test, so
  truly intertidal terrain — bed between LWL and HWL, including
  rest-dry columns — keeps real metrics and floods/dries dynamically.
  The seed floors each layer to `2·H_VANISHED` (never 0, and above the
  strict `h_old > H_VANISHED` remap-drain gate so seeded S/T survives
  the first regrid), so a starts-dry column sits at a *vanished*
  `D ≈ nz·2·H_VANISHED > 0` with bed-blocked faces until overtopped —
  not literal `D = 0`.  Caveat: formula-topo land at `b ≈ 0` (island /
  double_drake), `neverworld2` continents, and make_bathy `b = 1 m`
  land-flag cells are all intertidal/wet under the knob-on criterion —
  true static land needs `b < -land_margin`.  `configure_ocean_wetdry`
  emits a WARNING via a condition-based post-seed scan (counts interior
  cells with `b ∈ [-land_margin, LAND_DEPTH_THRESHOLD)`), so this is
  caught for name-less topos and file bathymetry too.  **Zero-depth (D≈0)
  hardening** (Finding A, `docs/ocean_wetdry_zero_depth_plan.md`): the
  EOS `h > 0` gates are tightened to `h > H_VANISHED` (a mid-drain layer
  in `(0, H_VANISHED]` returns the reference density, not a corrupted
  `hS/h`), the seed + `bt_h` floor and the drain-surviving `2·H_VANISHED`
  value are in place, and the controlled vanished-column full-driver test
  (`wetdry_vanished_column_flood`) confirms a seeded-vanished interior
  column with a density gradient floods cleanly, no NaN, tracers conserved.
  DEFERRED: a proof at *exactly* `D = 0` (all layers at literal 0, no
  floor) — the shipped path never produces it (seed + limiter keep
  `h ≥ 2·H_VANISHED` / `h ≥ 0`), so this is a defence-in-depth corner, not
  an operational gap.  Restart round-trip is
  BIT-EXACT (`wd_wet_dyn` hysteresis state registry-carried; transient
  wd_* recomputed) — gated by `restart_bit_exact_wetdry`.
  `test_ocean_wetdry` (+ `wetdry_starts_dry`, `wetdry_seed_criterion`,
  `wetdry_emerged_beach_multilayer_step`, `wetdry_emerged_tracer_gradient`,
  `wetdry_vanished_column_flood`)
  + `test_ocean_wetdry_driver` (full-driver tracer conservation through a
  dry→wet→dry band cycle: `Σ hS`/`Σ hT` drift 7.5e-16).  **Remaining v1
  gap** (unchanged by v2): the surface-flux dry gate is *binary* — a
  column held wet in the hysteresis band (`dry_depth < D < rewet_depth`)
  receives the full unthrottled heat/salt flux over its mm-scale top
  layer (no thickness-aware throttle), quantified in
  `test_ocean_wetdry_driver`.
- **Surface stress**: 2D `tau_x/tau_y` on the top layer; MOM6
  `DIRECT_STRESS` distributes it over `hmix_stress` metres.
- **Surface heat / salt flux** (`&ocean_thermo_nml q_heat / q_salt`):
  scalar `Q_heat` / `Q_salt` stamped into the top layer (`k=nz`).
- **Surface-flux component set** (`&ocean_forcing_nml enable_components`,
  default off ⇒ no array allocated, `Q_heat`/`Q_salt` unchanged): grows
  the scalar `Q_heat`/`Q_salt` fill into a MOM6-shaped component set on
  `ocean_surface_flux_t` — `q_sw`/`q_lw`/`q_lat`/`q_sens`/`heat_added`,
  the mass-flux set (`evap`/`lprec`/`fprec`/`vprec`/`lrunoff`/`frunoff`/
  `seaice_melt`), a `heat_content_*` enthalpy companion per mass flux,
  `salt_flux`, and `p_surf_atm`/`p_surf`. `Q_heat`/`Q_salt` become
  **derived views** rebuilt every thermo step by
  `ocean_surface_flux_assemble` from the const scalar + the components.
  **v1 limits, stated plainly**: mass fluxes are stored and carry
  enthalpy/salt bookkeeping but do **NOT** change column mass (no real
  freshwater — that is a named follow-up); `salt_flux` stays **virtual**
  (no column-mass change, same as the pre-PR-12 `Q_salt`); `p_surf`/
  `p_surf_atm` are stored with **no consumer** yet (the inverse-barometer
  PGF fold is a same-release-cycle follow-up); `q_sw` is now a selectable
  irradiance source for shortwave penetration + the boundary-layer SW
  coupling (`&ocean_thermo_nml sw_source="q_sw"`, PR-21 — requires
  `enable_components`; default `"net_heat"` stays bit-identical). The only v1 filler is the
  sea-ice coupler (`rdb_ice_ocean_coupler`), which writes `salt_flux`/
  `heat_added` instead of overwriting `Q_salt`/`Q_heat` directly when
  components are on — bit-identical to the legacy full-overwrite path.
- **Geothermal bottom heat flux** (`&ocean_geothermal_nml enable /
  q_geo`, default off): bed-side analogue of the surface heat flux —
  a constant `Q_geo` (W/m², positive into the ocean from below;
  Davies/Pollack global mean ~0.05–0.1) deposited into the lowest
  *massive* layer (`k=1` in the common case, falling to the first
  layer above a pinched ZSTAR_FULL bed). Heat only, thermo cadence,
  `heat_budget_geothermal` accounting contributor.
- **Vertical advection** (Eulerian): w diagnosed from continuity;
  tracer + h_layer conservation gated by remap.
- **ALE remap**: PPM stencil (`remap_method="ppm"`, default); momentum,
  h_layer, all tracers; drift ≤ 1e-15 / 10 remaps. Optional **PPM_H4**
  (`remap_method="ppm_h4"`) swaps the interior edge estimate for the
  thickness-weighted non-uniform 4th-order stencil (White & Adcroft 2008),
  which stays 4th-order on non-uniform ALE layers where plain PPM degrades to
  2nd — cutting the spurious diapycnal mixing injected per remap (Ilicak 2012).
  Reuses the PPM limiter + parabola + conservative redistribute verbatim
  (conservation/monotonicity unchanged). Boundary edges use PCM-outermost +
  3-cell H3; the 4×4 cubic-fit boundary is a documented upgrade. Optional
  **PQM** (`remap_method="pqm"`) — White & Adcroft (2008) `PQM_IH4IH3`:
  implicit-h4 edge values + implicit-h3 edge slopes (per-column tridiagonal
  solves) + monotonicity limiter, ~5th-order convergence on smooth profiles
  vs PPM's 2nd. **Reachable** (both `&vcoord_nml remap_method` and
  `&ocean_diag_nml diag_remap_scheme`); `N<5` silently falls back to PPM
  (documented design, not a configure-time error — the canonical
  `double_gyre_mom6.nml` runs `nz=2`, so selecting `pqm` there is PPM).
  PQM's own boundary cells still collapse to PCM (the 4×4 cubic-fit upgrade
  above applies to PQM too — deferred, a kernel change not a wiring one).
  The remap (and the Fox-Kemper fold) are **gated on the
  thermo cadence** (`dt_therm_ratio`, MOM6 DT_THERM): every outer step at the
  default `ratio=1` (bit-identical), once per thermo interval at `ratio>1`
  (Lagrangian-then-remap, conservative; the Lagrangian state is valid for the
  fixed-z diag output and bit-exact for restart via the saved `outer_step_count`).
  Horizontal tracer advection is NOT yet on the coarser cadence — the deferred
  flux-accumulation phase decouples `DT_TRACER_ADVECT` next.
- **Along-coordinate tracer Laplacian** (`&ocean_hdiff_nml kappa_h`,
  default `0.0` = off): constant-coefficient conservative curvilinear
  flux-form diffusion on `T = hTr/h`, thickness-weighted faces. Diffuses
  along the MODEL coordinate, not neutral surfaces — on a σ/z* grid over a
  slope that has a diapycnal component (accepted; `&ocean_redi_nml` is the
  separate neutral path). Physical-edge wall closure (not array-edge) —
  correct under both single-rank all-wall and MPI-decomposed runs.
  Configure-time guard (`rdb_ocean_stability_audit.F90`, part of the
  stability audit below) aborts on `κ_h·dt_therm·(1/dx_min²+1/dy_min²) >
  0.5`, using the REAL per-cell metric minimum on any grid type (fixed
  from an earlier Cartesian-only, nominal-`dx`/`dy` check that silently
  skipped spherical/tripolar grids).
- **Configure-time stability audit** (`rdb_ocean_stability_audit.F90`,
  runs once per solver creation, after the real metric arrays are built)
  — named-number, named-fix diagnostics for the failure modes that used
  to surface only as a bare `[nan-catch]` count mid-run: (1) viscous CFL
  `nu_h·dt/dx_min² > 0.125` (ERROR unless `bound_kh`/`stress_tensor` is
  enabled, in which case informational only — the runtime clamp already
  protects it); (2) the `kappa_h` diffusive number above (ERROR); (3)
  Munk sidewall boundary-layer resolution `delta_M >= 2` cells (WARNING);
  (4) `&ocean_hvisc_nml ah_max < nu_h` silently clamping the configured
  viscosity (WARNING). All four use the ACTUAL minimum grid cell
  (`ocean_stability_min_cell`), never nominal `&grid_nml dx`/`dy`
  (degrees on spherical/tripolar). The companion runtime diagnostic
  (`apply_velocity_truncation`'s NaN-catch, `rdb_ocean_dyn.F90`) reports
  the `(i,j,k)` of the first non-finite face plus the local cell size and
  viscous CFL there, gated behind the existing catch so it costs nothing
  on a healthy run.
- **Tracer packages** — `multilayer_state_t%register_passive_tracer(grid,
  name, units, long_name, idx)` (6-arg; `idx=0` on refusal, `registry_locked`
  after `enter_data`) grows the S+T[+age] registry at setup; every
  passive-transport kernel (advection, ALE remap, vertical exchange,
  vertical/horizontal diffusion, Redi, sponge, halo/fold, OBC ghost +
  reservoirs, restart) already loops the registry, so a newly-registered
  tracer rides them for free — only a named EOS/surface-flux coupling and a
  diag `fill_<name>` (the ocean diag surface is NOT registry-driven) need
  hand-wiring. Budget attribution (`heat_budget_*`/`salt_budget_*`) is a
  registry property (`tracer_t%budget_id`), not an index coincidence — a
  passive tracer defaults to no budget slot. Shipped packages: **ideal age**
  (`&ocean_tracers_nml enable_ideal_age`, MOM6 `USE_IDEAL_AGE_TRACER`
  analogue — 1 s/s interior aging, surface reset) and **pseudo-salt**
  (`&ocean_tracers_nml enable_pseudo_salt`, Shao 2016 verification tracer —
  seeded to S, given exactly S's surface salt flux + KPP/EPBL non-local
  mirror; the deviation `pseudo_salt − S` measures the passive-vs-active
  transport-path error; fail-loud excluded from SSS restoring + sea-ice,
  both un-mirrored salinity sources). Both default off ⇒ bit-identical.
  Deferred: a namelist-declarable arbitrary dye/CFC/BGC package (needs a
  generic registry-driven diag-fill path, not just the registration API).

### Shipped — coordinates + BCs + I/O

- **All ten VCOORD_* families** dispatch through the same ALE remap:
  `LAGRANGIAN`, `EULERIAN_Z`, `SIGMA`, `ZSIGMA`, `ZSTAR`, `ZSTAR_FULL`,
  `ZSTAR_SIGMA`, `Z_FIXED`, `RHO`, `HYCOM` (enum of record:
  `src/core/rdb_constants.F90`). Eight are `select case` branches of
  `ocean_vcoord_compute_target_h`; `RHO` and `HYCOM` need per-layer T/S plus
  the EOS and so enter through the sibling `compute_target_h_rho`, dispatched
  from the same remap driver (`src/ALE/rdb_ocean_remap.F90`). `LAGRANGIAN`
  early-returns — the target IS the live `h_layer`, so the remap is a no-op.
  On the ocean path `ZSTAR` shares the `SIGMA` branch: in this barotropic
  `(H, η)` form the two target formulas are the same expression.
  `RHO` is validation-grade alone (weakly-stratified columns collapse);
  `HYCOM` is the production hybrid.
- **OBC dispatch wired end-to-end** (2026-06-10, `&ocean_bc_nml` →
  per-edge tags → driver → kernels). Shipped types: WALL (default),
  OPEN (Flather + per-layer zero-gradient baroclinic anomaly), TIDAL
  (multi-constituent η), CLAMPED (Dirichlet η/u/v + per-tracer inflow
  values), CHAPMAN (scalar edge-mean Orlanski), SPONGE (momentum decay
  + optional tracer relaxation via `sponge_relax_tracers`), and
  **PERIODIC** (per-axis ghost-wrap reentrant boundary; requires
  `nghost >= 3`, paired edges; bit-exact seam — a circularly shifted
  IC reproduces the shifted solution bit-for-bit). NESTED behaves as
  OPEN pending the parent-data orchestrator; INFLOW/DISCHARGE are
  tag-reserved only. Open-edge tracer ghosts are upwind-aware
  (inflow uses the per-edge boundary value, outflow zero-gradient),
  or — with `res_lscale_in/out` — evolve through per-face **tracer
  reservoirs** (implicit relaxation between interior and boundary
  data at flow-dependent rates; no sign-switch chatter). Radiating
  edges optionally use **per-layer Orlanski radiation**
  (`radiation_scheme="orlanski"`, running-mean phase speed,
  `rx_max` clip) with asymmetric inflow/outflow **nudging**
  (`nudge_tau_in/out`), and the barotropic Flather supports the
  **full half-characteristic form** with exterior velocity
  (`flather_form="full"`, per-edge `*_ext_u/v`). Boundary-corner
  relative vorticity is zeroed at every non-periodic edge.
- **Horizontal grids** (`&ocean_grid_nml grid_config`): `cartesian`
  (default), `spherical` lon-lat sector, `supergrid` (MOM6 mosaic
  reader), and **`tripolar`** — Murray (1996) bipolar Arctic cap above
  `phi_join` + lon-lat below, closed by a single-rank **north fold**
  (`north="tripolar_fold"`, requires periodic west/east). The fold
  reverses-i + sign-flips vector normals + antisymmetrically projects
  the duplicated v/corner seam row; metric + `f_corner` ghosts are
  folded once at configure. Kernels read full 2D metric arrays only.
  Vector→geographic output rotation at the seam is deferred (output
  shows model-frame velocities in the cap).
- **Diag manager**: registry + cadence + procedure-pointer fill
  dispatch. Default device-resident fills for `SSH, T, S, u_centre,
  v_centre, KE`; CONSERVATIVE vertical remap onto fixed-z
  (`DIAG_VGRID_Z_FIXED`), sigma, z*, or density bins (`DIAG_VGRID_DENSITY`,
  **namelist-reachable**: `&ocean_diag_nml vgrid="density"` + `rho_levels`/
  `n_rho_levels`, or per-diagnostic `name:density`/`name:rho` in `diags`;
  configure-time guard requires `n_rho_levels > 0` and strictly increasing —
  density bins have no auto-fill, unlike sigma/z*) —
  donor-cell overlap integral via the shared `remap_column` kernel,
  intensive + extensive variants, density targets reusing the RHO
  `invert_density_targets` solve; time-ops INSTANT / MEAN / MAX / MIN /
  **INTEGRAL** (`name:integral` — cumulative `Σ(sample·dt)` over the cadence
  window, undivided; the budget-closure operator, distinct from MEAN);
  serial per-rank NetCDF emit (CF-1.8). `&ocean_diag_nml diags` requesting
  a canonical diagnostic whose feature gate is closed (e.g. `age` with
  `enable_ideal_age=.false.`, `ice_conc` with the ice off) now
  `logger%warning`s naming the diagnostic and the likely gate, instead of
  silently dropping it — `name:off` stays a silent no-op either way.
- **Restart / checkpoint** (bit-exact, MPI-native registry):
  per-slot restart registry (`restart_registry_t`) — every
  prognostic-owning slot registers its arrays at state-build time, the
  manager (`ocean_state_restart_write` / `_read`) walks the registry
  with the `!$acc update self` device-pull discipline (device-mapped
  entries only), writes each field as a full local array (interior +
  ghosts), and writes a per-rank `restart_rank_NNNNNN.nc` durably
  (`.tmp` + POSIX-rename) with decomposition + grid/vcoord/tracer +
  schema metadata, validated on read. Checkpoints barotropic
  (`h, u_face_x, v_face_y`), multilayer (`h_layer, u/v_face_x/y_layer`,
  every `tracers(:)%hTr`, `rho_layer` — carried across dt_therm-skipped
  steps), KPP lagged `bl_depth`, EPBL `mld`+`kd_int` and kappa-shear
  `kd_int`+`tke_int` (merged every stage, refreshed only at thermo
  cadence — so required for bit-exactness), the live BC persistent state
  (Chapman `eta_old_chapman_*` scalars + tracer reservoirs `tres_*` when
  allocated), and `outer_step_count` (dt_therm alignment). A missing
  REQUIRED field on read is fatal. Resume requires the SAME
  decomposition + matching grid/vcoord/tracer metadata (errors loudly on
  mismatch — fatal when the driver omits `ierr`); cross-rank
  redistribution is a future offline tool. Driver-wired at
  `cfg%restart_interval` cadence + clean end; warm-start from
  `cfg%restart_file`. Bit-exact step-(N+1) roundtrip gate:
  `test_ocean_restart` (closed-wall + periodic-x seam + physics-rich
  KPP/heat-flux/dt_therm + decomp-mismatch negatives via both the check
  path and the production read wrapper), green on GPU. Sea-ice EVP
  dynamics is covered too: the ice->ocean momentum mediation
  (`surface_stress%tau_x/y`) is carried on the ice slot
  (`ice_tau_ocn_x/y` + the `ice_tau_ocn_valid` presence scalar, both
  `optional=.true.`) and COPIED back at configure by
  `ice_ocean_stress_resume_apply` — formula-agnostic by construction,
  replacing an earlier reconstruct-from-checkpointed-concentration fold
  that was NOT bit-exact whenever a checkpoint step's thermo/transport
  changed ice concentration after the blend. Gate:
  `restart_bit_exact_ice_evp` (`tests/test_ocean_ice_restart.F90`),
  alongside `restart_fresh_run_is_pure_wind` (no restart file ⇒ pure
  wind, not a zeroed stress) and `restart_old_checkpoint_degrades_to_wind`
  (a pre-this-capability checkpoint resumes with a documented
  one-window cold-start of the τ mediation, not a fatal). Diag
  MEAN/MAX/MIN windows reset on restart by design. The dead
  `ocean_obc_t` scaffold is NOT checkpointed (never enabled /
  device-mapped); OBC registers its own live state when it lands.
- **Initial conditions**: analytical ICs (uniform / linear-T(z) /
  Eady-front / geostrophic-adjustment / per-layer `gprime` density)
  plus **z-level T/S from NetCDF** (`&ocean_zinit_nml enable`,
  capability A2). The file must be **pre-regridded to the model
  horizontal grid** (`nx_phys × ny_phys`, bathymetry-loader
  precedent); A2 does the in-core vertical step — linear-in-depth
  interpolation of `temp`/`salt` onto the seeded layer-centre depths
  with constant extrapolation beyond the source z-range, written as
  `hTr = value · h_layer`. Dry columns (`wet_mask ≤ 0`) get the
  namelist `land_fill_t/_s` constants. Default off (`enable=.false.`)
  ⇒ analytical IC bit-identical. Source `temp` is taken as the
  prognostic T directly. **Not yet shipped:** on-the-fly horizontal
  regrid + nearest-wet land flood-fill (MOM6 `horiz_interp_and_extrap`
  half — deferred to a Python preprocessor + v2), and conservative
  (cell-integral-preserving) vertical remap (v1 uses point linear
  interp at layer centres; the first ALE remap re-grids anyway).

### Sea ice

**A high-quality sea-ice dynamical core and column model, not yet a
sea-ice model.** SIS2 port under `src/core/ice/` (14 modules, ~6,900
lines; ~6,400 lines of tests), gated by `&ocean_ice_nml enable` (default
off ⇒ the slot is never initialised, mapped, or stepped ⇒
byte-identical). The physics that is there can be trusted; the physics
that is missing is most of the ice mass budget. Limits first:

- **No ridging.** `compress_ice` (`rdb_ice_transport.F90`) is a
  thinnest-first cascade that returns `part_size(0) ≥ 0` by in-place
  compaction/promotion — mass- and area-conserving, but it thickens ice
  **in its own category** at zero energetic cost. That is not a ridging
  parameterisation: no participation function (which thin ice deforms),
  no redistribution function (where the deformed mass goes), no work
  done against gravity. **It is SIS2's own `DO_RIDGING=.false.`
  fallback** — SIS2's `compress_ice`/`else` branch (`SIS_transport.F90`)
  is the identical routine and role, and SIS2 itself calls it *"a
  minimalist version of a sea-ice ridging scheme."* Roundabout matches
  SIS2's shipped default and is missing SIS2's option, which is itself
  an Icepack wrapper (`ice_ridge.F90`), not native SIS2 code — porting
  it means vendoring Icepack, not a small follow-up. Consequence: under
  convergence, ice thickens uniformly instead of building a thick
  ridged tail — the ITD is wrong in a convergent regime (too little ice
  in the thick categories), and because the EVP ice strength
  `pres_mice = (p0/ρ_ice)·exp(−c0·max(1−ci,0))` reads that distribution,
  the strength is biased with it. This bites hardest exactly where EVP
  matters most: convergent, near-shore, and Antarctic-shear regimes.
- **No snowfall source.** `m_snow` is allocated, transported, melted,
  and rebalanced — and can only ever *decrease* (zero occurrences of
  `snowfall`/`fprec`/`precip` under `src/core/ice/`). Consequence: the
  snow machinery is present but effectively unreachable outside a
  restart — the CSIM4 snow-vs-bare-ice albedo blend
  (`rdb_ice_optics.F90`), the snow-covered melt point, and the snow/ice
  interface conductivity never fire on a cold-start run.
- **No snow-ice flooding.** No freeboard or submergence calculation
  anywhere in `src/core/ice/` — the mass module docstring
  (`rdb_ice_mass.F90`) states freeboard/flooding mass paths are
  "deliberately NOT ported." Consequence: where snow load depresses the
  freeboard below sea level, flooding (snow → snow-ice conversion) is
  unrepresented; in the Antarctic this is a large fraction of total ice
  mass.
- **No melt ponds.** `m_pond = 0.0_wp` is a hardwired dead local in
  `rdb_ice_column.F90`, kept only so the surrounding branch matches
  SIS2's shape — there is no pond state on `ocean_sea_ice_t`.
  Consequence: summer melt is biased low and the melt-albedo feedback
  is absent; the CSIM4 melting-albedo ramp is a crude stand-in, not a
  representation of ponds.
- **No lateral melt / floe-size effects.** Surface and basal melt only
  (zero occurrences of `lateral_melt`/`floe` under `src/core/ice/`) —
  MIZ retreat is biased.
- **No ice initial-condition path.** No IC namelist group, no reader —
  ice can enter a run only by frazil growth from an ice-free ocean, or
  by restart.
- **Coupling is partial.** Salt is closed (virtual brine-rejection
  flux, Boussinesq — the ocean's water mass is not reduced). Heat is
  closed, including the transmitted shortwave `sw_thru` (shortwave
  penetrating the ice into the ocean below): PR 31 reduces the
  per-category `sw_thru` to a per-cell field and `ice_ocean_sw_flux`
  (`rdb_ice_ocean_coupler.F90`) delivers it to the ocean — into the
  `q_sw` surface-flux component when `&ocean_forcing_nml
  enable_components` is on (assembled into `Q_heat`), else added
  directly into `Q_heat`. The v1 divergences from SIS2 that remain:
  the ocean's own background shortwave is not lead-fraction weighted by
  ice cover (the configure-time `q_heat`/`q_sw` scalars stay ice-blind),
  and the ocean shortwave is single-broadband (no `VIS_DIF`-style
  spectral band assignment). Momentum reaches the ocean as the
  concentration-weighted blend of the wind snapshot and the EVP
  ice-ocean drag (`ice_ocean_stress_flux`), and that blend also
  refreshes the derived cell-centred `|tau|`
  (`ocean_surface_stress_refresh_mag`), so the boundary-layer schemes'
  friction velocity `u_* = sqrt(|tau|/rho0)` follows the ice-mediated
  stress instead of the configure-time wind — under full ice cover in a
  windless run the pre-fix `u_*` was identically zero
  (`test_ocean_ice_stress_mag`). Momentum is still **not conserved at
  fractional ice cover**: the ice feels the full wind-stress snapshot
  rather than a concentration-weighted bulk drag law (documented
  divergence D7, `rdb_ice_evp.F90`). There is **no freshwater/mass
  coupling at all** (`rdb_ice_frazil_uptake.F90`: "Freshwater/mass
  coupling is PR-3c+ territory") — ice carries no weight, so it applies
  no dynamic sea-surface loading.
- **Frazil is surface-layer-only.** The supercooling clamp/bank
  (`rdb_ice_frazil`) checks only the top layer (`k=nz`); MOM6/SIS2
  check the full water column for supercooled water.

**Shipped** (each independently validated):

- **Winton (2000) two-layer column thermodynamics + enthalpy**
  (`&ocean_ice_nml enable`, `nk_ice=2`) — brine-salinity-dependent,
  energy-conserving, Stefan-verified to <1% relative at days
  1/5/10/30 over a 32-day hourly integration.
  `test_ocean_ice_column`, `test_ocean_ice_enthalpy`.
- **Multi-category ITD restore** (`ncat>1`) — whole-category shift, not
  Lipscomb (2001) remap (SIS2 doesn't ship that either — its "remap" is
  a `SIS_transport.F90` TODO comment, never code). `test_ocean_ice_itd`.
- **Category ice/snow transport + `compress_ice`** (`&ocean_ice_nml
  transport`) — category-summed PPM (reuses `rdb_continuity`) with a
  proportionate per-category flux split and PCM tracer riding.
  `test_ocean_ice_transport`.
- **C-grid EVP rheology** (`&ocean_ice_nml dynamics`) — full elliptical
  yield curve, replacement pressure with a grid/Tdamp-scaled floor,
  device-resident subcycled loop, correct C-grid staggering; golden-
  tested against the analytic Nansen free-drift solution
  `|u| = √(τ/(ρ_ocean·c_dw))` to 15 digits. With `dynamics=.false.` ice
  velocity falls back to sampling the ocean surface layer.
  `test_ocean_ice_evp`.
- **Frazil bank + uptake** and **ice → ocean coupling** (salt closed;
  heat closed, incl. transmitted shortwave `sw_thru`, PR 31) —
  `test_ocean_frazil`, `test_ocean_ice_coupling`,
  `test_ocean_ice_driver_column`, `test_ocean_ice_diags`.

**Envelope.** The whole slot is fail-loud single-rank: `enable=.true.`
alone already requires `px*py=1` (`rdb_config.F90`, covers the
thermo-only path, not just transport/dynamics), and `transport` /
`dynamics` carry their own identical guards on top. `transport`
additionally forbids **any** periodic edge; `dynamics` (EVP) allows
wall or periodic edges but forbids `north="tripolar_fold"` and any
OBC/tidal/sponge/clamped/Chapman edge — there is no open-boundary
support for ice at all. There is no ice halo-exchange code anywhere in
`src/core/ice/`. Combined with the absent IC path above, there is today
no reachable configuration that is both realistic (frazil-grown or
restarted ice, on a mixed real edge set) and multi-rank.

### Deferred (planned, not yet shipped)

- **Nonlinear collapse of the Eady front** (`validation_examples/ocean/eady/`,
  2026-09-13) — past day ~65 the `dT_dy = -2e-5` front's eddies pass 1 m/s,
  hit the 2 m/s `maxvel` clamp and NaN-catch by day ~75 at dt = 600 s
  (`n_inner` 20 and 30, clamp removed, `nu_h = 100` all fail); dt = 300 s
  survives to day 100 but with `En` above the front's available potential
  energy. The shipped 60-day window is the validated linear phase.
- **Resting-state growth at low Laplacian viscosity** (same dir, open) — a
  motionless, stably stratified periodic channel seeded with ±0.5 mK
  white-noise T grows `En` from 0 to 3e-5 m²/s² in 25 days at `nu_h = 0`
  (Smagorinsky + `nu_4 = 5e8` on; broadband in x, smooth bed/surface-
  intensified vertical structure, e-folding 2.5 days, and the resting T
  range shrinks 0.5 K at each end). `nu_h = 20` lowers the level 230× at
  day 25 but not the rate; `nu_h = 100` holds it at 1.6e-11. Not yet
  attributed to a kernel.

- **Hollingsworth-Källén Coriolis correction** — Arakawa-Hsu 4-corner
  weighted-PV stencil; the PV form is current production.
- **Full KPP Phase 2+** — `V_t²` unresolved-turbulence term in the
  bulk-Ri sweep + full non-local γ_T/γ_S transport + bottom BL
  extension.
- **Smagorinsky and biharmonic lateral kernels** declared with enums
  (`LMIX_SMAGORINSKY`, `LMIX_LEITH_BIHARM`, `LMIX_BIHARMONIC`);
  Smag_KH + Smag_AH are live; Leith biharmonic is the next item.
- **Non-hydrostatic on the C-grid** — hydrostatic only today; a
  non-hydrostatic pressure correction on the C-grid layout is
  unimplemented.
- **Tides** — equilibrium + SAL ship (`&ocean_tides_nml`); OBC-tide
  nodal correction ships (`&ocean_bc_nml obc_tidal_nodal`). Internal-tide
  drag on the barotropic mode ships too, but as a SEPARATE capability —
  `&ocean_bt_nml wave_drag` (barotropic linear/Rayleigh wave drag, Egbert
  & Ray 2001 / Jayne & St Laurent 2001), independent of `&ocean_tides_nml`
  and Cartesian-compatible. Its `r_H(x,y)` map is `form="uniform"` (a
  global scalar) or `form="roughness_proxy"` (a resolved-bathymetry-
  variance PLACEHOLDER for the real subgrid `⟨h²⟩`, not a substitute for
  it); `form="file"` (a real subgrid-roughness map) is registered but
  fails loud at configure — the NetCDF reader is PR-14.
- **MPI I/O server hand-off** for diag manager (serial NetCDF is
  the current emit path).
- **Per-layer Orlanski phase-speed radiation** at open edges (v1
  ships Flather mean + zero-gradient baroclinic anomaly) and
  file-backed / tidal-table boundary data sources (the polymorphic
  `update(t, bc)` call site is wired in the driver; only the
  constant backend ships).

### Validation surface

`validation_examples/ocean/` carries the canonical benchmarks:

| Subdir | What it tests |
|---|---|
| `seamount/` | Quiescent IC over Gaussian seamount; zero motion forever; caught three latent ocean bugs in May-21. |
| `geostrophic_adjustment/` | Rossby SSH bump on f-plane; IG-wave radiation at √(gH); 2π/f oscillation at centre. |
| `eady/` | Baroclinic instability of a thermal-wind front in a periodic channel; the kx = 2 (125 km) channel mode grows at 1.93e-6 1/s vs 1.98e-6 from Eady (1949) theory (re-baselined 2026-09-13: `dT_dy = -2e-5`, `nu_h = 20`, 60 days). Linear-growth benchmark only — see the nml header for why the run stops before the front's nonlinear collapse. |
| `eddy_test/` | Stratified wind-driven β-plane double-gyre; WBC + eddy shedding. |
| `double_gyre/` | MOM6 ocean_only/double_gyre reproduction; stable 580+ days. |
| `neverworld2/` | Idealized single-basin (Marques et al. 2022): 60°×140° spherical sector, re-entrant southern (Drake) channel, MOM6-inspired continents + banded zonal wind. v1 = geometry bring-up (wind-only, linear-T EOS, zstar+KPP); SST-restore + thermocline IC are fast-follows. |
| `baroclinic_channel/` | Eady/Phillips baroclinic instability in a periodic-x channel; the clean, fast, well-resolved eddy-generation test (growth → roll-up → saturation, conservation, bounded CFL). 2-layer (~3 min) + 15-layer demo. |
| `island_at_rest/` | Quiescent basin with an interior land block; stays exactly at rest (zero motion, mass-exact) — the land-mask Tier-1.5 trap. |
| `flow_past_island/` | Wind-driven flow onto an island; verifies no-normal-flow (velocity ≡ 0 inside land) + mass conservation. |
| `double_drake/` | Ferreira et al. 2010 two-continent world (seam-straddling meridional walls + reentrant southern channel); land-masking showcase. |
| `sea_ice/` | Polar freeze-up: `polar_freezeup_thermo.nml` (column + frazil + brine) and `polar_freezeup_dynamics.nml` (EVP drift) — physics-sanity runs, not a performance benchmark. |

### Working envelope (proven production-stable)

```
&ocean_hvisc_nml      nu_h = 10000.0,  smag_ah = .true. /
&ocean_bdrag_nml      form = "linear", hbbl = 10.0, bg_vel = 0.1, r = 2.5e-5 /
&ocean_coriolis_nml   form = "sadourny" /
&ocean_pgf_nml        form = "gprime", maxvel = 6.0 /   ! 2-layer reduced gravity + clamp
                                                        ! (optional cfl_trunc = 0.5 adds an
                                                        !  advective-CFL velocity truncation
                                                        !  ahead of the absolute maxvel cap)
&ocean_bt_nml         auto_n_inner = .true. /
&vcoord_nml           vcoord_type = "sigma" /
&time_nml             dt_fixed = 1200.0 /
```

(Ocean knobs split into per-concern sub-namelists as of the
namelist-UX refactor — `&ocean_<group>_nml` with the `ocean_`
prefix dropped from each key.  Run `tools/nml_split.py <file>` to
migrate an old `&ocean_setup_nml` namelist.)

This is the regime that hits day 580 on the MOM6-ref double-gyre.
Outside this envelope (lower `nu_h`, alternative drag form, alternative
PGF on real bathymetry) needs case-by-case validation.

### Headline ocean performance

| Config | Wall time (single V100) |
|---|---|
| Tasman 2 km (eddy-resolving, KPP + ALE + sponges) | **~34 s / simulated day** |
| Double-gyre MOM6 ref (44 × 40 × 2, dt=1200) | ~16 s / 30 simulated days |

---

## Conservation contract (console Mass / Salt / Heat `Error`)

The periodic `[stats]` console block prints a `Mass : <total>  Error <residual>` line for each conserved quantity.  The `Error` meaning depends on the regime and the run configuration.

### What the `Error` measures

For a conserved quantity Q the residual is

```
Error = (Q_total - Q_ref) + out_Q - src_Q
```

where `Q_ref` is the value latched on the first status report, `out_Q` is the cumulative net outflux through open boundaries (positive = left the domain), and `src_Q` is the cumulative surface source (positive = added to the domain).  On a conservative, closed-domain, unforced run all three terms are zero and `Error` is the true numerical-leak residual **subject to the summation-order floor described below**.  With open boundaries or surface forcing the budget correction keeps `Error` a round-off residual rather than a spurious flag for the intended physics.

When a budget term is not instrumented (see limitations below) the console falls back to raw drift `(Q_total - Q_ref) / Q_ref`, byte-identical to the pre-budget behaviour.

### The summation-order floor, and the `reproducing_sums` escape hatch (PR-32)

By default every `Q_total` above is a plain floating-point `!$acc parallel loop reduction(+:acc)` device reduction, combined across ranks by `MPI_SUM` on doubles.  Neither is associative or order-deterministic: a changed rank count, a changed domain decomposition, or a changed GPU reduction-tree shape can all change the last few bits of `Q_total` — and `Q_ref` is latched as a `real(wp)` (double), so even a perfectly exact sum feeding an unchanged `(Q_total - Q_ref)` subtraction would still be quantised to ~1 ulp of `Q_total` (for a basin-scale `Q_total ~ 1e21`, that ulp is ~1.3e5 in absolute units — precisely the band the `Error` column is asked to resolve).  Two consequences:

- **`Total mass` / `Total KE` / `Total salt` / `Total heat` are NOT rank-count-reproducible** with `reproducing_sums = .false.` (the default) — the printed totals drift at round-off (~1e-14 to 1e-15 relative) when the decomposition changes, even for bit-identical physics.
- **The `Error` floor is set by summation noise, not physics**, at that same ~1e-14 to 1e-15 relative band.

`&ocean_diag_nml reproducing_sums = .true.` routes `Total mass` / `Total KE` / `Total salt` / `Total heat` (and the sea-ice area totals) through an Extended-Fixed-Point (EFP) reproducing sum (Hallberg & Adcroft 2014; `src/framework/rdb_efp.F90` + the comm-facade `halo_allreduce_efp_list`, ONE collective in place of the default path's several separate `MPI_SUM` calls) and forms the `Error` residual by differencing in FIXED POINT (`efp_real_diff`) rather than subtracting two already-quantised doubles.  With it on:

- The printed totals become **bit-identical across rank counts and reduction/decomposition orders** (EFP's decomposition is a pure function of each summand's value; integer bin addition is exact and order-invariant, modulo the bounds below).
- The residual floor drops to the EFP quantum (`2^-3P` per summand, aggregated `<= N * 2^-3P` over `N` cells) — around `1e-54` relative against a `~1e21` total, i.e. the floor becomes irrelevant to any physically-meaningful drift rather than merely "exact".
- **Rank envelope**: `EFP_MAX_RANKS = 2^(53 - 36) = 131072`.  The cross-rank combine transports the six fixed-point bins as exactly-representable `real64` values (there is no `integer(int64)` MPI allreduce in `pic_mpi_lib`), which is exact only while every partial sum a rank count could form stays `<= 2^53`; both this bound and the per-summand bin-1 bound are enforced fail-loud (`error stop`), never a silent fallback.  131072 ranks is far beyond any Roundabout run.
- Default is `.false.` ⇒ the console block (and every latched `console_stats_t` field) is byte-identical to the pre-PR-32 output — this knob changes diagnostic TEXT only, never the prognostic trajectory.
- Scope: the primary totals feeding `Error` (Mass/KE/Salt/Heat + sea-ice area) go through EFP; the salt/heat closed-budget `out`/`src` terms (the boundary-outflux and surface-source corrections above) remain on the pre-existing FP `halo_allreduce_sum` path in both branches — their individual magnitudes are typically much smaller than the total, so the summation-order floor is proportionately less consequential there. `compute_max_cfl` is untouched by this knob in either branch: a global max is already exact and order-invariant in floating point, so there is nothing to fix.

### Instrumented paths

| Regime | Mass | Salt | Heat |
|---|---|---|---|
| Ocean (`sim_type='ocean'`, split RK2) | closed residual | closed residual | closed residual |
| Ocean (`sim_type='ocean'`, unsplit) | closed residual | closed residual | closed residual |

Windowed tracer advection (`&ocean_vmix_nml dt_tracer_advect_ratio > 1`) is
instrumented too, since 2026-09-12: `continuity_tracer_drain`'s sub-cycle and
BOTH halves of the per-stage concentration hold accumulate into the same
`*_budget_horiz_adv` arrays, so the closed residual holds at every report, not
only on window boundaries.  It used to fall back to raw drift, which cannot
subtract a source — a conservative run with a surface heat flux then reported
the heat the flux legitimately added as a ~5e-5 "leak" (four shipped
`acc_channel` namelists; now ~3e-14).

### Known limitations (fall back to raw drift)

- **Redi + open boundary** (`&ocean_redi_nml enable = .true.` with any open-edge BC): the neutral-diffusion flux exits through the open boundary without being counted in `out`; reverts to raw drift.

---

## What Roundabout is not

- **Not a semi-implicit free-surface model.** The split-explicit barotropic sub-cycling gives O(minute) outer steps, but that is sub-cycling the fast mode, not treating the gravity-wave term implicitly in the SCHISM/Casulli sense. There is no implicit η-solve on this tree.
- **Not a wave model, and not non-hydrostatic.** The dynamical core is hydrostatic; there is no phase-resolving surface-gravity-wave capability, no wave-maker / absorbing BCs, and no external wave-model (SWAN/WW3) coupling.
- **Not a coastal / estuarine shock-capturing solver.** Dam-break, hydraulic-jump, tidal-bore and inundation work lives in the separate coastal repository — the Riemann/KNP machinery is not in this tree.
- **Not a BGC model.** No biology source terms, no FABM coupling.
- **Not a sediment / morphodynamic model.** Bathymetry is fixed.
- **A high-quality sea-ice dynamical core and column model, not yet a
  complete sea-ice model.** Winton column thermodynamics, the
  multi-category ITD, category transport + `compress_ice`, C-grid
  EVP dynamics, an analytic ice initial-condition path, a snowfall
  source term, and Archimedes snow-ice flooding are landed (default off,
  `&ocean_ice_nml`) and independently validated (Nansen free-drift to 15
  digits, Stefan melt to <1%) — see the "Sea ice" subsection above for
  the full limits list. The physics that is there can be trusted; the
  physics that is missing is most of the ice mass budget: no ridging
  (`compress_ice` is SIS2's own `DO_RIDGING=.false.` fallback, not a
  participation/redistribution scheme), no melt ponds, and no lateral
  melt. Snow-ice flooding (`snow_ice`) carries a known SIS2-inherited
  approximation — the converted ice takes the SNOW's enthalpy and zero
  salinity rather than the flooding SEAWATER's, so it is too cold and
  too fresh; true seawater-based flooding (ocean mass/heat/salt sink
  from refreezing pore water) is a materially larger, deferred closure.
  Coupling is partial (transmitted shortwave `sw_thru` is now coupled
  to the ocean — PR 31; ice<->ocean
  momentum conserved at fractional cover only with the opt-in
  `&ocean_ice_nml a_face_stress` — the default-off legacy form leaks
  `(1-a)(tau_a-fxoc)` per face, no freshwater/mass exchange), and the
  whole slot is fail-loud single-rank.
- **Not a fully equipped basin-scale climate ocean model — yet.** The
  ocean path (`sim_type='ocean'`) is now operational for regional /
  basin hydrostatic configurations (see the "Ocean path" section
  above), but is missing tides, SAL correction,
  and the Hollingsworth-Källén Coriolis guard that production
  climate-scale models rely on. Eddy-resolving regional ocean (Tasman
  2 km, MOM6-ref double-gyre) works today; ice-coupled global climate
  doesn't.

---

## Compiler / platform support

- **NVHPC** (GPU + multicore stdpar) — primary target, exercised in CI.
- **gfortran** — exercised in CI for portability.
- **ifx** — buildable; some OpenMP-target map-clause edge cases observed in vendor tooling, not load-bearing for the canonical NVHPC build.
- **AMD / Intel GPUs** — unsupported today (no HIP / SYCL backend; would be a portability research project).

---

## How to read this alongside the roadmaps

The roadmap (`ROADMAP_OCEAN.md`) carries the full work checklist with project-scope tags (`[E]/[M]/[H]/[V]`) and `Validate:` lines per planned item. This document summarises the *current* state; if there's a discrepancy, the codebase is the authority and the roadmaps are the next-best source.
