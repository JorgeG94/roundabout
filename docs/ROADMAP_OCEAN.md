# Roundabout Ocean dyn-core — roadmap

Long-horizon backlog for the ocean path (`sim_type='ocean'`). The executed
phase-map / design-of-record (what the `Phase N` anchors in the ocean code cite)
is the [**Executed phase map**](#executed-phase-map-design-of-record) section at
the bottom of this file. Cross-cutting numerics and engineering (WENO, SSP-RK3,
IMEX, mixed precision, DA, ML, HPC/CS research) are not otherwise recorded here.

Last reconciled with the code: **2026-09-15**.

Scope tags: **`[E]`** honours / summer · **`[M]`** masters / 1st-PhD-chapter · **`[H]`** full PhD · **`[V]`** research programme. Each item names a *Validate:* target. `(tentative)` = judgement call whether to do it at all.

## Where the dyn-core sits

A **MOM6-architecture generalized-coordinate (ALE) dynamical core**, run in its
**z\*/sigma mode**, GPU-native, aimed at the **eddy-resolving regional** niche:

- Hydrostatic Boussinesq **Arakawa C-grid**; layer thickness is prognostic —
  continuity is a transport equation (continuity-PPM), not a free-surface Poisson solve.
- **Split-explicit** — slow baroclinic + fast forward-backward barotropic subcycle,
  CFL-derived `n_inner`. Two outer schemes: `pred_corr` (default; predictor-corrector,
  tendencies on the barotropic-window means) and `ssp_rk2` (Heun, experimental).
- **Vector-invariant momentum** — Sadourny PV-flux Coriolis + KE-gradient; FV pressure gradient — Boussinesq Montgomery potential (default), z-corrected,
  Wright-Picard, MOM6-analytic, gprime.
- Wright nonlinear EOS; Leith + Smagorinsky lateral; PP81 + KPP vertical; PPM ALE remap.

On the coordinate axis it now spans both the *geometric* GVC families
(sigma / z\* / hybrids + Lagrangian / Eulerian endpoints) **and** the *isopycnal*
interior — `VCOORD_RHO` (density-space regrid) and `VCOORD_HYCOM` (Bleck hybrid)
ship, alongside the full mesoscale-eddy closure stack (isopycnal slopes → GM
thickness diffusion, Redi neutral diffusion, VarMix coefficients, prognostic
MEKE). GM remains *opt-in / off by design* at the 1 km submesoscale-resolving
target (where you resolve the mesoscale), but neutral (Redi) diffusion is useful
even eddy-resolving to control z\* spurious diapycnal mixing. On the model map it
sits with regional eddy-resolving z\* models (MITgcm / ROMS / CROCO / NEMO
regional), GPU-native being the differentiator.

**Status:** Tier-1 operational — double-gyre stable to day 580; 2 km Tasman
eddy-resolving at ~34 s / simulated-day on one V100. Geometric + isopycnal/hybrid
coordinates, mesoscale-eddy closures, PERIODIC (reentrant) boundaries, and
windowed tracer advection all shipped, alongside the sea-ice train
(`&ocean_ice_nml`). The C-grid MPI multi-rank halo ships and is gated by a
decomposition-invariance suite (1 rank vs N, same answer) plus the
`tests/mpi/` set; OBC dispatch is the remaining regional-enabling gap
(tidal forcing ships; its validation campaign does not). Several individual
features remain single-rank — porous barriers, wet/dry, sea ice, the
tripolar north fold and the windowed tracer-advection drain. The user guide is a
Sphinx site under `rdb_docs/` (theory, discretisation, coordinates, PGF,
closures, running simulations, the Python driver, validation); per-procedure
reference comes from the FORD `!!` docstrings and is built separately.

## The fork

Two directions the dyn-core can grow:

- **Lane 1 — finish the regional tool:** OBCs and tidal validation (§1). PBCs and
  the C-grid MPI halo are done. Lets you run a *real* regional domain at scale; the GPU-native angle is
  the differentiator. **This is the remaining near-term priority.**
- **Lane 2 — isopycnal / water-mass fidelity:** *largely delivered* — the RHO +
  HYCOM coordinates and the neutral-diffusion / GM / VarMix / MEKE stack now ship
  (§2). Remaining lane-2 items are partial cells and the NH-on-ocean / 1 km-global
  reach (§3–4).

## 1. Regional-enabling (lane 1)

- [~] `[M]` **OBCs** — finish NESTED / CLAMPED / SPONGE dispatch, a **file-based nesting
  data source** (only `constant` exists today), tracer OBC values, validated clean
  radiation. *Very important* — without them you're limited to closed basins.
  (v0.1.0 milestone; detail in the checklist.) *Validate:* an outgoing Rossby / gravity
  wave radiates cleanly; one-way nest reproduces the parent at the boundary.
- [x] `[M]` **PBCs** — periodic / reentrant boundaries. **Shipped:** per-axis ghost-wrap
  (`nghost >= 3`, paired edges, bit-exact seam — shifted IC ⇒ shifted solution
  bit-for-bit, `test_ocean_periodic`); Eady runs as a true periodic-x channel
  (`validation_examples/ocean/eady/`). Unlocked the idealized channels (Eady; ACC pending a run).
- [~] `[M]` **Tides** — **the forcing ships**: the equilibrium body force + scalar
  self-attraction & loading (`&ocean_tides_nml enable` / `use_sal` / `beta_sal`),
  the 18.6-yr nodal correction (`add_nodal`, also reaching OBC_TIDAL via
  `&ocean_bc_nml obc_tidal_nodal`), and barotropic linear wave drag
  (`&ocean_bt_nml wave_drag`) all run off the shared astro generator
  (`rdb_ocean_tide_astro`); configs in `validation_examples/ocean/tides/`.
  Remaining: the *validation* campaign. *Validate:* Australia-wide barotropic
  vs TPXO M2/S2/K1/O1; SAL phase shift ~10–15%.
- [x] `[M]` **C-grid MPI multi-rank halo** — **shipped:** staggered face/centre/corner
  halos + decomp-gated walls, exercised by `tests/mpi/` (`test_halo_ocean_mpi`,
  `test_ocean_dyn_mpi`, `test_ocean_bt_cfl_mpi`) and by the
  decomposition-invariance gate (`tests/regression/decomp_invariance.py`), which
  exists because two real bugs were invisible at one rank *by construction*.
  Remaining: strong-scaling characterisation, and the features still flagged
  single-rank (porous barriers, wet/dry, sea ice, tripolar fold, windowed
  tracer-advect drain).
- [ ] `[E]` **HK Coriolis** — `sadourny_hk` ships opt-in; PV / Sadourny is the validated
  default. Validate HK, then promote or leave opt-in. *Validate:* Hollingsworth instability
  suppressed on a pathology test.
- [ ] `[E]` **S&H dsig-stretch option** — if ROMS Song-Haidvogel stretching is needed,
  expose it as a `stretch_type="song_haidvogel"` option feeding the existing sigma path
  (`dsig`-builder ~30 lines); the dead `VCOORD_SCOORD` enum has been removed (D7).
  *Validate:* generated `dsig` matches the S&H analytic formula for given `(θ_s, θ_b, h_c)`.

## 2. Coordinate & mixing physics (the isopycnal fork — lane 2)

- [x] `[M]` **True isopycnal (RHO) vertical coordinate** — **shipped:** `VCOORD_RHO`
  locates each target-density surface per column (PPM density reconstruction +
  Newton root-find) and remaps the interfaces onto it; `VCOORD_HYCOM` adds the
  Bleck hybrid (isopycnal interior + z\*-like surface floor). All ten VCOORD
  families dispatch through the same ALE remap path. (The `VCOORD_LAGRANGIAN`
  endpoint still drifts off-isopycnal under diabatic forcing by design — RHO is
  the re-pinning path.)
- [x] `[H]` **Neutral / Redi diffusion + Gent-McWilliams** — **shipped:** isopycnal
  slopes (Griffies 1998) → GM thickness diffusion (bolus), Redi neutral diffusion
  (continuous sweep + GPU flux), VarMix coefficients (Res_fn + Visbeck/Eady), and
  prognostic MEKE closing the GM↔MEKE loop. **Caveat unchanged:** GM is a
  *coarse*-resolution parameterization — at the 1–2 km eddy-resolving target run it
  *off*; neutral diffusion stays useful to control z\* spurious diapycnal mixing.
- [ ] `[H]` **Partial cells** (cut cells at bathymetry) — only with a z-coordinate. *Validate:*
  internal-tide over a shelf break vs terrain-following ref; no small-cell dt penalty.

## 3. Non-hydrostatic on the C-grid

- [ ] `[H]` **NH-on-ocean** — a Stelling-Zijlema-class non-hydrostatic pressure
  correction (CG-Poisson) on the C-grid ocean layout. This was scoped as an *adaptation
  of the proven coastal sigma-coordinate NH solver*; that implementation (and the
  elliptic-solver machinery it shared with the semi-implicit free surface) left with
  the coastal carve-out, so the item stands but the reference implementation is no
  longer in this tree — build it from the published method, or against the coastal
  repository. *Validate:* solitary / internal-wave propagation against an analytic
  dispersion relation (the coastal NH solver is no longer available as the in-tree
  reference).

## 4. Scale

- [ ] `[V]` **1 km global** — the original 1 km-submesoscale target (Phase 7+ in the phase
  map): multi-GPU throughput + load balance at basin / global scale on the shipped MPI halo. *Validate:*
  submesoscale-resolving global run at sustained throughput; conservation over a multi-month integration.

## 5. Sea ice

No longer tentative — the ladder was built. `src/core/ice/` ships the Winton
enthalpy column thermodynamics, multi-category ITD, ITD-aware transport, C-grid
EVP dynamics, frazil, and the ice-ocean salt/heat/stress couplers, all under
`&ocean_ice_nml` (default off ⇒ byte-identical). Remaining: ridging, the
ice-shell / shelf-cavity interface kernel, and a polar validation campaign.

## 6. Ocean validation

- [ ] `[M]` **MOM6 cross-comparison campaign** — beyond the double-gyre: baroclinic instability
  (Phillips / Eady), overflow / DOME, seamount / internal-tide. *Validate:* published MOM6
  reference within tolerance per case.
- [~] `[M]` **Eady / ACC channel** — PBC shipped; Eady runs as a periodic-x channel
  (`validation_examples/ocean/eady/`) and, since 2026-09-13, reproduces the linear-theory
  growth rate of its channel mode (1.93e-6 vs 1.98e-6 1/s; gated in
  `tests/regression/stability_manifest.py`). Remaining: an ACC transport run; the front's
  nonlinear collapse past day ~65 is not integrable at dt = 600 s (open, see the nml header).

## 7. Time integration

- [x] `[M]` **Outer split scheme** — **shipped:** `&ocean_bt_nml split_scheme`
  selects `pred_corr` (default; predictor-corrector, tendencies on the
  barotropic-window means `u_av`/`h_av`, forward-backward gravity-wave pairing)
  or `ssp_rk2` (Heun, experimental). The Coriolis reference the fast loop
  subtracts is evaluated on the same velocity the slow tendency used, under the
  same weights, so the two cancel exactly (MOM6 `ubt_Cor`); evaluating it at
  stage entry instead pumped a barotropic seiche. Both vertical-friction
  tridiagonals gained an `nz = 1` path, and `scratch_3d` buffers are attached
  device-zeroed rather than as allocator leftovers.
- [ ] `[E]` **`ssp_rk2` residual growth** — Heun amplifies oscillatory modes at
  `(ω·Δt)⁴`, so a stratified rest state manufactures internal-wave energy. Gated
  by `validation_examples/ocean/eady/resting_stratified_channel.nml`, which every
  scheme must pass on magnitude; the *settle* criterion is a scoped expected
  failure. Either retire `ssp_rk2` or give the outer loop a neutral / dissipative
  option (forward-backward is neutral; an RK3 outer step damps at the same order).
  *Validate:* a motionless stratified periodic channel does not grow.
- [ ] `[M]` **Cross-rank bitwise reproducibility** — `En`/`MaxCFL` diverge across
  rank counts while mass stays bit-identical; non-reproducible FP reductions are
  the suspect. EFP reproducing sums ship (`reproducing_sums`, default off).
  *Validate:* identical time series at 104 vs 208 ranks.

## Prioritization (ocean)

0. **Keep the gates honest** — the two-tier stability suite
   (`tests/regression/stability.py`) asserts on each run's own time series rather
   than on goldens, runs every unpinned case under both outer schemes, and is the
   reason several silent defects surfaced at all. New physics lands with a case
   that fails without it.
1. **Regional-enabling first** — OBC dispatch (PBC and the C-grid MPI halo done; tidal forcing done, tidal validation open). Unlocks real regional domains; the distinctive GPU-native value.
2. ~~Fork decision~~ — *resolved*: the RHO/HYCOM coordinates + neutral/GM/VarMix/MEKE stack shipped. Remaining coordinate/mixing work is partial cells (§2) only.
3. **NH-on-ocean + 1 km-global** — larger, later.

## Executed phase map (design-of-record)

The `Phase N` anchors the ocean source cites (~152 comments). Tier-1 dyn-core is operational — the MOM6 double-gyre runs stable to day 580 on one V100.

- **Phase -1** — Directory refactor to MOM6-style organise-by-concern (`core/{coastal,ocean}` + `ALE/`, `tracer/`, `parameterizations/`, `pressure_force/`, `equation_of_state/`, `framework/`); `sim_type` dispatcher. **shipped** (`core/coastal/` and the `sim_type` dispatcher's second branch have since been carved out; `sim_type` is pinned to `"ocean"`)
- **Phase 0** — Ocean scaffolding: `rdb_ocean_solver` stub + empty `core/ocean/` dirs (real stepping lives in `rdb_ocean_dyn.F90`). **shipped** (the `rdb_ocean_solver` Phase-0 abort stub was deleted as dead code — zero `use` sites, PR-8 — once `rdb_ocean_dyn` / `driver_run_ocean` made it vestigial)
- **Phase 0d** — Ocean god-state + per-kernel state types declared, all component `init`s no-op. **shipped**
- **Phase 0e** — Ocean slot shells (vmix / diag / vcoord / tides / budgets) declared, arrays unallocated. **shipped** (the units-registry / debug-probe / data-override / river scaffolding declared in this phase was deleted as dead code — no kernel, no knob, PR-8)
- **Phase 0f** — Framework infra shells: NetCDF-reader hook, safe-math, checksums interface decls. **shipped** (`rdb_checksums` deleted as dead code — a size-only stub masquerading as a hash, PR-8; `rdb_safe_math` kept as a tested library with no production consumer, its build-time dispatch option `RDB_BITWISE_REPRO` removed, PR-8)
- **Phase 1** — C-grid barotropic state (`barotropic_state_t`, edge velocities) + `scratch_3d_buffer_t`; device round-trip test. **shipped**
- **Phase 2** — Continuity-PPM kernel, barotropic only — per-face mass fluxes, `∂h/∂t+∇·(hu)=0`. **shipped**
- **Phase 2a** — First-order upwind continuity baseline (`PPM_VARIANT_UPWIND`). **shipped**
- **Phase 2b** — PPM face reconstruction (Colella-Woodward). **shipped**
- **Phase 3** — PV-conserving Coriolis-advection (vector-invariant, KE-gradient). **shipped**
- **Phase 3a** — Plain 4-point C-grid Coriolis force (`du/dt=f·v`). **shipped**
- **Phase 3b** — Sadourny (1975) energy-conserving form (corner ζ, centre KE). **shipped**
- **Phase 3c** — Hollingsworth-Källén correction (`sadourny_hk`, Arakawa-Hsu 6-pass). **shipped**
- **Phase 4** — Split-explicit dynamics driver (`ocean_dyn_step_split`, BT-to-BC feedback). **shipped**
- **Phase 4a** — Unsplit SSP-RK2 (Heun) reference path, no fast/slow split. **shipped**
- **Phase 4b** — Split-explicit fast-barotropic substep loop + CFL-derived `n_inner`. **shipped**
- **Phase 5** — Wire tracer + base physics (PP81, lateral mix, PGF) through the ocean path. **shipped**
- **Phase 5a** — Leith closure (+ Smagorinsky_KH / biharmonic Smagorinsky_AH overlay). **shipped**
- **Phase 5b** — KPP boundary-layer vertical mixing (+ V_t², non-local γ_T/γ_S). **shipped**
- **Phase 5c** — Wright (1997) nonlinear EOS. **shipped**
- **Phase 5d** — FV pressure force (MONT / FV_LITE / FV_WRIGHT / GPRIME / FV_MOM6 variants). **shipped** (MONT became a true Boussinesq Montgomery potential — geopotential recursion + thickness-weighted density term — and the default, 2026-09-12)
- **Phase 5e** — Tidal interior body-force + SAL + nodal correction (`rdb_ocean_tides.F90` + `rdb_ocean_tide_astro.F90`; OBC_TIDAL ships; barotropic linear wave drag ships). **shipped** (the TPXO validation campaign is still open — see §1)
- **Phase 5f** — Non-hydrostatic on the C-grid ocean layout (no kernels; the coastal NH solver that was to be the reference left with the coastal carve-out). **TODO**
- **Phase 5g** — Vertical coord + ALE remap on the ocean path (six VCOORD families, PPM remap, centre + face). **shipped**
- **Phase 6** — Diagnostics manager + ocean validation suite. **partial**
- **Phase 6a** — Diag registry + cadence + log sink (SSH / T / S / u / v / KE default fills). **shipped**
- **Phase 6b** — Windowed flux-accumulated tracer advection (`dt_tracer_advect_ratio` knob, `&ocean_vmix_nml`): `uhtr`/`vhtr` face-transport accumulators + `t_dyn_rel_adv` in `continuity_t`; swept-average CW-PPM drain with fixed-budget CFL sub-cycling; mandatory flush before output/restart/remap; MOM6 `DT_TRACER_ADVECT` cadence (Adcroft & Hallberg 2006). Single-rank; multi-rank halo deferred to E1. **shipped** (default cadence stays
`dt_tracer_advect_ratio=1` fused per-stage; the windowed drain is the opt-in
`ratio>1` path — see `docs/howto/tracer_advect_cadence.md`).
- **Phase 6c** — Land/ocean mask + layer→fixed-z diag vremap + time-mean ops; budgets path. **shipped**
- **Phase 6d** — Serial NetCDF emit (per-rank, CF-1.8). **shipped** (the MPI I/O-server rank was coastal-path machinery and left with the carve-out; the ocean path emits per-rank directly)
- **Phase 7** — Multi-GPU C-grid MPI halo (face/centre/corner); single-GPU perf characterised. **shipped** (halo + decomp-invariance gate; multi-GPU strong scaling not yet characterised)
- **Phase 8** — Deferred v2: TEOS-10, stochastic, rivers (**MEKE + internal-tide /
  tidal interior mixing + sea ice now shipped**; A↔C nesting is moot — there is
  no A-grid path in this repository any more). **partial**

Dyn-core is Tier-1 operational in closed / sponged / periodic basins, single- and multi-rank; the forward backlog (OBC dispatch, tidal validation, isopycnal fork, NH-on-ocean, 1 km-global) lives in the sections above.
