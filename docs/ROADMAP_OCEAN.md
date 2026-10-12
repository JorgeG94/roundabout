# Roundabout Ocean dyn-core — roadmap

Last reconciled with the code: **2026-10-12** (beta at `635f35d92`).

## Where we are today

- The ocean path (`sim_type='ocean'`) is **Tier-1 operational**. The MOM6 double-gyre runs stable to day 580. A 2 km Tasman case runs at about 34 s per simulated day on one V100.
- It is a MOM6-style ALE core on an Arakawa C-grid. Ten vertical-coordinate families ship, including isopycnal (`RHO`) and hybrid (`HYCOM`).
- The mesoscale stack (GM, Redi, VarMix, MEKE), sea ice, periodic boundaries, tides (forcing), and the multi-rank C-grid halo all ship. A decomposed run is bit-identical to a serial run.
- The biggest gaps are **file-backed open-boundary data (nesting)**, **tide validation against TPXO**, and **non-hydrostatic physics**.
- Current push: **ACCESS-OM3 parity**. C1, C7 part 1, C3, C2a and C2b are merged. C2c (the atmospheric-state provider) is in progress.

## How to read this file

1. **[Done](#1-done)** is what already ships, one line each, with evidence.
2. **[In progress and next](#2-in-progress-and-next)** is the work that is open or near-term.
3. **[Backlog and deferred](#3-backlog-and-deferred)** is later work, grouped by theme.
4. **[Appendices](#appendices)** hold the rationale, measurements and the full phase map. Source comments cite the `Phase N` labels in Appendix A.

Scope tags: **`[E]`** honours / summer · **`[M]`** masters / 1st-PhD-chapter · **`[H]`** full PhD · **`[V]`** research programme.
`(tentative)` means it is a judgement call whether to do it at all.
Each open item names a *Validate:* target.

Where this file and the code disagree, the code wins. Per-feature detail lives in [`CAPABILITIES_AND_LIMITATIONS.md`](CAPABILITIES_AND_LIMITATIONS.md) and [`CLOSURE_MATRIX.md`](CLOSURE_MATRIX.md).
Cross-cutting numerics and engineering (WENO, SSP-RK3, IMEX, mixed precision, DA, ML, HPC/CS research) are not otherwise recorded here.

---

## 1. Done

### Dynamics, time stepping and grids

| Item | What shipped |
|---|---|
| C-grid core | Continuity-PPM, PV-conserving Coriolis, split-explicit barotropic loop (Phases 1-4b). |
| Outer split scheme | `&ocean_bt_nml split_scheme`: `pred_corr` (default) or `ssp_rk2` (experimental). See [Appendix B](#appendix-b-detail-for-done-and-open-items). |
| Periodic boundaries (PBCs) | Per-axis ghost-wrap, bit-exact seam (`test_ocean_periodic`). Eady runs as a periodic-x channel. |
| Tripolar north fold | Bit-identical on every split, including `px > 1` (`b35c5f63e`, `e8d614aa1`, `42adac5ff`, `6d91130b2`, `40236aba0`, `11ed42542`). |
| Hollingsworth-Källén Coriolis | `sadourny_hk` ships opt-in (Phase 3c). Promotion is still open, see §2. |
| Vertical coordinates | All ten `VCOORD_*` families dispatch through one ALE remap path (the `VCOORD_*` enum in `src/core/rdb_constants.F90`). |
| True isopycnal coordinate | `VCOORD_RHO` and `VCOORD_HYCOM` (Bleck hybrid). |
| Partial bed cells on `z_fixed` | Partial bed cell, bed fillers and closed-face mask ship (`tests/test_ocean_zfixed_closed_faces.F90`, `test_ocean_zfixed_k_bot.F90`). General cut cells are still open, see §3. |

### Physics

| Item | What shipped |
|---|---|
| Mesoscale closures | Isopycnal slopes, GM thickness diffusion, Redi neutral diffusion, VarMix, prognostic MEKE. |
| Tides (forcing) | Equilibrium body force, scalar SAL, nodal correction (also for OBC_TIDAL), barotropic wave drag. The TPXO validation is still open, see §2. |
| Sea ice | `src/core/ice/`: Winton thermodynamics, multi-category ITD, ITD-aware transport, C-grid EVP, frazil, ice-ocean couplers. All under `&ocean_ice_nml` (default off). |
| Sea ice on many ranks | Thermo, EVP and stress blend on the ocean decomposition (`84ddefc87`); category transport across ranks and periodic seams (`313a8c0c1`); restart round trip (`9df1c5b90`). |
| Sea ice under the tripolar fold | One rank only: category state (`cc75df88f`) and EVP dynamics (A6, `e1ae86fa7`). More than one rank is still refused (A7, #143). |
| Windowed tracer advection | `dt_tracer_advect_ratio`; single-rank. |

### Parallelism and reproducibility

| Item | What shipped |
|---|---|
| C-grid MPI halo | Face, centre and corner halos. Gated by `tests/mpi/` and the decomposition-invariance suite. |
| Cross-rank bitwise reproducibility | `tests/mpi/test_ocean_decomp_bitid_mpi.F90` compares every restart-registry field bitwise on every `px x py` split. The console uses EFP reproducing sums, now the default (`reproducing_sums` in `ocean_diag_config_t`, `src/core/rdb_config.F90`). |
| Porous barriers on many ranks | Covered by the decomposition bit-identity test (the `periodic_channel_zstar` case in `tests/mpi/test_ocean_decomp_bitid_mpi.F90`). |

### Boundaries, validation and ACCESS-OM3 parity

| Item | What shipped |
|---|---|
| OBC dispatch | WALL, OPEN (Flather), TIDAL, CLAMPED, CHAPMAN, SPONGE, PERIODIC, Orlanski radiation, asymmetric nudging. Nesting data is still open, see §2. |
| Eady channel | Linear growth rate reproduced (1.93e-6 vs 1.98e-6 1/s), gated in `tests/regression/stability_manifest.py`. An ACC run is still open. |
| Global 1-degree + sea ice | One year on one V100 (`d6e46b6d7`, 7.45 s per simulated day). |
| **C1** 75-level vertical grid | ACCESS-OM3's exact grid runs (`69310e028`, `tools/vgrid_to_dz.py`, `tests/test_ocean_vgrid_om3.F90`). |
| **C7 part 1** OM3 diagnostics | `umo`/`vmo`, `tauuo`/`tauvo`, `Kd_interface` (`89f78f07c`, `src/core/ocean/diag/rdb_ocean_diag_derived.F90`). |
| **C3** Large-Yeager bulk flux | Standalone column kernel plus a 2-D driver (`babb12d6e`, review fixes `e963d3469`; `src/parameterizations/vertical/rdb_ocean_bulk_flux.F90`). **Not wired in:** no namelist, no engine call site. |
| **C2a** calendar and multi-file time axis | `rdb_calendar`, `&time_nml start_date`/`calendar`, `register_2d_filelist` (`3d4ff94d7`, review fixes `d3921d57a`; `src/core/rdb_calendar.F90`). |
| **C2b** lat-lon regridding weights | Bilinear and first-order conservative weights onto spherical, supergrid or tripolar grids (`b1b4b94b3`, `c9006a540`; `src/core/ocean/io/rdb_regrid_latlon.F90`). **Not wired in:** no reader, no namelist. Conservative scratch scaling is #202. |

### Code health (October 2026)

| Item | Commit |
|---|---|
| `rdb_config` split into submodules under `src/core/config/` | `22ee5f5ed` |
| `rdb_ocean_setup` split into submodules under `src/core/ocean/state/setup/` | `65fdf0d48` |
| `rdb_ocean_state` host-only code moved to `src/core/ocean/state/ocean_state/` | `cdcbdcd16` |
| Halo pack/unpack loops converted to `do concurrent` | `22c1d4ba9` |
| Reductions converted to `do concurrent reduce(...)` | `4c8434801` |
| Closure-matrix check also runs on the config submodules | `a7b0b9c79` |

---

## 2. In progress and next

| Item | Goal | Status |
|---|---|---|
| `[M]` **C2c** atmospheric-state provider | JRA55-do reader, regrid and bulk-flux wiring, in three chunks. | C2c-1 (config and reader window) is open PR [#205](https://github.com/JorgeG94/roundabout/pull/205). Blocked on #202 for the conservative regrid at 1440x1080. |
| **C3 wiring** | Connect the bulk-flux kernel to a namelist, an atmospheric-state reader and an engine call site. | Not started. The C3 commit (`babb12d6e`) names C2b/C2c as the work packages for the reader and regridder. |
| `[M]` **OBC nesting data** | A file-backed data source, so a one-way nest reproduces its parent. | Only the `constant` backend exists (`ocean_boundary_data_constant_t` in `src/core/ocean/boundary/rdb_ocean_boundary_data.F90`, held as `bc_source` in `ocean_engine_t`). `OBC_NESTED` is tag-reserved and behaves as OPEN (`OBC_NESTED` in `src/core/ocean/boundary/rdb_ocean_boundary_types.F90`). |
| `[M]` **Tide validation** | Australia-wide barotropic run against TPXO M2/S2/K1/O1. | Forcing ships. No TPXO comparison in `validation_examples/ocean/tides/`. |
| `[M]` **MOM6 cross-comparison** | Phillips/Eady, overflow/DOME, seamount/internal tide against published MOM6. | Not started beyond the double-gyre. |
| `[M]` **ACC channel** | An ACC transport run on the periodic channel. | Eady part done; ACC run pending (`validation_examples/ocean/acc_channel/` exists). |
| `[E]` **HK Coriolis** | Validate `sadourny_hk` on a Hollingsworth-instability test, then promote or leave opt-in. | Opt-in; PV/Sadourny is the validated default. |
| `[E]` **`ssp_rk2` residual growth** | Retire `ssp_rk2` or give the outer loop a neutral option. | Scoped XFAIL on the *settle* gate. |

<details>
<summary>Validate targets for these items</summary>

- **OBCs:** an outgoing Rossby or gravity wave radiates cleanly; a one-way nest reproduces the parent at the boundary. Also open: tracer OBC values from file and validated clean radiation (v0.1.0 milestone).
- **Tides:** Australia-wide barotropic run vs TPXO M2/S2/K1/O1; SAL phase shift about 10-15 %.
- **MOM6 cross-comparison:** published MOM6 reference within tolerance per case.
- **HK Coriolis:** Hollingsworth instability suppressed on a pathology test.
- **`ssp_rk2`:** a motionless stratified periodic channel does not grow.

</details>

---

## 3. Backlog and deferred

### Regional enabling

| Item | Note |
|---|---|
| `[E]` **Song-Haidvogel `dsig` stretch** | Expose ROMS stretching as `stretch_type="song_haidvogel"` feeding the sigma path (`dsig` builder, about 30 lines). The dead `VCOORD_SCOORD` enum is gone (D7). Not in `src/` today. *Validate:* generated `dsig` matches the S&H formula for given `(θ_s, θ_b, h_c)`. |
| `[H]` **Per-layer Orlanski radiation** | Per-edge depth-uniform Orlanski ships. Per-layer phase speed does not. |
| `[M]` **Remaining single-rank features** | Ice-shelf cavity, wet/dry, Chapman edges, `dt_tracer_advect_ratio > 1`, sea ice under the tripolar fold, in-memory geometry injection. Each fails loud on more than one rank. |
| `[M]` **Restarts across a rank-count change** | Restarts resume on the same decomposition only. |
| `[M]` **Strong-scaling characterisation** | Multi-GPU strong scaling is not yet measured. |
| `[M]` **Large-rank-count reproducibility re-run** | The decomposition bit-identity test covers 1, 2 and 4 ranks. A 104 vs 208 rank comparison of the console time series (`En`, `MaxCFL`) has not been re-run on this tree. *Validate:* identical time series at 104 and 208 ranks. |

### Coordinates and mixing

| Item | Note |
|---|---|
| `[H]` **Cut cells at bathymetry** | Partial bed cells exist on `z_fixed` only. General cut cells would also serve other z-coordinates. *Validate:* internal tide over a shelf break vs a terrain-following reference; no small-cell `dt` penalty. |
| `[M]` **Hollingsworth-Källén promotion** | See §2. |
| `[E]` **Retire `zstar_full`** | Slated once z\* is proven on the 1-degree Southern Ocean and the global case. |

### Non-hydrostatic and scale

| Item | Note |
|---|---|
| `[H]` **NH-on-ocean** (Phase 5f) | A Stelling-Zijlema-class pressure correction (CG-Poisson) on the C-grid. No kernels in the tree; build from the published method. *Validate:* solitary or internal-wave propagation against an analytic dispersion relation. |
| `[V]` **1 km global** (Phase 7+) | Multi-GPU throughput and load balance at basin or global scale. *Validate:* a submesoscale-resolving global run at sustained throughput, with conservation over a multi-month integration. The 1-degree global case already runs on one GPU. |

### Sea ice

| Item | Note |
|---|---|
| `[M]` **Ridging** | `compress_ice` is SIS2's `DO_RIDGING=.false.` fallback. |
| `[M]` **Ice-shell / shelf-cavity interface kernel** | Not started. |
| `[M]` **Polar validation campaign** | Not started. |

### Physics deferred (Phase 8)

| Item | Note |
|---|---|
| TEOS-10 | `EOS_VARIANT_TEOS10` is declared but marked "future" (`EOS_VARIANT_TEOS10` in `src/equation_of_state/rdb_eos.F90`). |
| Stochastic parameterisations | Nothing in `src/`. |
| Rivers | Runoff fields exist in `rdb_ocean_surface_flux.F90`. The `ocean_river_t` slot listed in `src/core/ocean/README.md` has no source file. |

### Performance and infrastructure

See "Performance and scaling" in [`CAPABILITIES_AND_LIMITATIONS.md`](CAPABILITIES_AND_LIMITATIONS.md): mixed precision, dynamic load balancing, NVLink-aware MPI tuning, a GPU CI runner.

---

## Priorities

0. **Keep the gates honest.** The two-tier stability suite (`tests/regression/stability.py`) asserts on each run's own time series, not on goldens. It runs every unpinned case under both outer schemes. New physics lands with a case that fails without it.
1. **Regional enabling.** Nesting data and tide validation. This unlocks real regional domains, which is the distinctive GPU-native value.
2. **ACCESS-OM3 parity.** Finish C2c, which wires C3.
3. **NH-on-ocean and 1 km global.** Larger and later.

The old "fork decision" between regional tool and isopycnal fidelity is resolved: the RHO/HYCOM coordinates and the neutral/GM/VarMix/MEKE stack shipped.

---

## Appendices

### Appendix A: Executed phase map (design-of-record)

The `Phase N` anchors the ocean source cites (about 152 comments). Tier-1 dyn-core is operational. The MOM6 double-gyre runs stable to day 580 on one V100.

<details>
<summary>All phases (click to expand)</summary>

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
- **Phase 5e** — Tidal interior body-force + SAL + nodal correction (`rdb_ocean_tides.F90` + `rdb_ocean_tide_astro.F90`; OBC_TIDAL ships; barotropic linear wave drag ships). **shipped** (the TPXO validation campaign is still open — see §2)
- **Phase 5f** — Non-hydrostatic on the C-grid ocean layout (no kernels; the coastal NH solver that was to be the reference left with the coastal carve-out). **TODO**
- **Phase 5g** — Vertical coord + ALE remap on the ocean path (originally six VCOORD families; ten ship today — PPM remap, centre + face). **shipped**
- **Phase 6** — Diagnostics manager + ocean validation suite. **partial**
- **Phase 6a** — Diag registry + cadence + log sink (SSH / T / S / u / v / KE default fills). **shipped**
- **Phase 6b** — Windowed flux-accumulated tracer advection (`dt_tracer_advect_ratio` knob, `&ocean_vmix_nml`): `uhtr`/`vhtr` face-transport accumulators + `t_dyn_rel_adv` in `continuity_t`; swept-average CW-PPM drain with fixed-budget CFL sub-cycling; mandatory flush before output/restart/remap; MOM6 `DT_TRACER_ADVECT` cadence (Adcroft & Hallberg 2006). Single-rank; multi-rank halo deferred to E1. **shipped** (default cadence stays `dt_tracer_advect_ratio=1` fused per-stage; the windowed drain is the opt-in `ratio>1` path — see `docs/howto/tracer_advect_cadence.md`).
- **Phase 6c** — Land/ocean mask + layer→fixed-z diag vremap + time-mean ops; budgets path. **shipped**
- **Phase 6d** — Serial NetCDF emit (per-rank, CF-1.8). **shipped** (the MPI I/O-server rank was coastal-path machinery and left with the carve-out; the ocean path emits per-rank directly)
- **Phase 7** — Multi-GPU C-grid MPI halo (face/centre/corner); single-GPU perf characterised. **shipped** (halo + decomp-invariance gate; multi-GPU strong scaling not yet characterised)
- **Phase 8** — Deferred v2: TEOS-10, stochastic, rivers (**MEKE + internal-tide / tidal interior mixing + sea ice now shipped**; A↔C nesting is moot — there is no A-grid path in this repository any more). **partial**

Dyn-core is Tier-1 operational in closed / sponged / periodic basins, single- and multi-rank. The forward backlog is in §2 and §3.

</details>

### Appendix B: Detail for done and open items

<details>
<summary>What the dyn-core is (architecture summary)</summary>

A **MOM6-architecture generalized-coordinate (ALE) dynamical core**, run in its **z\*/sigma mode**, GPU-native, aimed at the **eddy-resolving regional** niche.

- Hydrostatic Boussinesq **Arakawa C-grid**; layer thickness is prognostic. Continuity is a transport equation (continuity-PPM), not a free-surface Poisson solve.
- **Split-explicit**: slow baroclinic step plus a fast forward-backward barotropic subcycle, CFL-derived `n_inner`. Two outer schemes: `pred_corr` (default; predictor-corrector, tendencies on the barotropic-window means) and `ssp_rk2` (Heun, experimental).
- **Vector-invariant momentum**: Sadourny PV-flux Coriolis plus KE-gradient. FV pressure gradient: Boussinesq Montgomery potential (default), z-corrected, Wright-Picard, MOM6-analytic, gprime.
- Wright nonlinear EOS; Leith + Smagorinsky lateral; PP81 + KPP vertical; PPM ALE remap.

On the coordinate axis it spans both the *geometric* GVC families (sigma / z\* / hybrids + Lagrangian / Eulerian endpoints) **and** the *isopycnal* interior: `VCOORD_RHO` (density-space regrid) and `VCOORD_HYCOM` (Bleck hybrid). The full mesoscale-eddy closure stack ships alongside (isopycnal slopes, GM thickness diffusion, Redi neutral diffusion, VarMix coefficients, prognostic MEKE).

GM remains *opt-in / off by design* at the 1 km submesoscale-resolving target, where you resolve the mesoscale. Neutral (Redi) diffusion is useful even eddy-resolving, to control z\* spurious diapycnal mixing. On the model map it sits with regional eddy-resolving z\* models (MITgcm / ROMS / CROCO / NEMO regional); GPU-native is the differentiator.

</details>

<details>
<summary>Status paragraph (long form)</summary>

Tier-1 operational: double-gyre stable to day 580; 2 km Tasman eddy-resolving at about 34 s per simulated day on one V100. Geometric + isopycnal/hybrid coordinates, mesoscale-eddy closures, PERIODIC (reentrant) boundaries, and windowed tracer advection all shipped, alongside the sea-ice train (`&ocean_ice_nml`).

The C-grid MPI multi-rank halo ships and is gated by a decomposition-invariance suite (1 rank vs N, same answer) plus the `tests/mpi/` set. OBC dispatch is wired; file-backed nesting data is the remaining regional-enabling gap (tidal forcing ships; its validation campaign does not).

Several individual features remain single-rank: the ice-shelf cavity, wet/dry, sea ice under the tripolar fold and the windowed tracer-advection drain. (Porous barriers, sea ice and the tripolar north fold with `px > 1` have since gone multi-rank; see §1.)

The user guide is a Sphinx site under `rdb_docs/` (theory, discretisation, coordinates, PGF, closures, running simulations, the Python driver, validation). Per-procedure reference comes from the FORD `!!` docstrings and is built separately.

</details>

<details>
<summary>The fork (history): regional tool vs isopycnal fidelity</summary>

Two directions the dyn-core could grow:

- **Lane 1, finish the regional tool:** OBCs and tidal validation. PBCs and the C-grid MPI halo are done. This lets you run a *real* regional domain at scale. **Still the remaining near-term priority.**
- **Lane 2, isopycnal / water-mass fidelity:** *largely delivered*. The RHO + HYCOM coordinates and the neutral-diffusion / GM / VarMix / MEKE stack ship. Remaining lane-2 items are cut cells and the NH-on-ocean / 1 km-global reach.

</details>

<details>
<summary>Periodic boundaries (PBCs), detail</summary>

Per-axis ghost-wrap (`nghost >= 3`, paired edges, bit-exact seam: shifted IC gives shifted solution bit-for-bit, `test_ocean_periodic`). Eady runs as a true periodic-x channel (`validation_examples/ocean/eady/`). This unlocked the idealized channels (Eady; ACC pending a run).

</details>

<details>
<summary>Tides, detail</summary>

The forcing ships: the equilibrium body force plus scalar self-attraction and loading (`&ocean_tides_nml enable` / `use_sal` / `beta_sal`), the 18.6-yr nodal correction (`add_nodal`, also reaching OBC_TIDAL via `&ocean_bc_nml obc_tidal_nodal`), and barotropic linear wave drag (`&ocean_bt_nml wave_drag`). All run off the shared astro generator (`rdb_ocean_tide_astro`); configs are in `validation_examples/ocean/tides/`.

Remaining: the *validation* campaign. *Validate:* Australia-wide barotropic run vs TPXO M2/S2/K1/O1; SAL phase shift about 10-15 %.

</details>

<details>
<summary>C-grid MPI halo, detail</summary>

Staggered face/centre/corner halos and decomposition-gated walls, exercised by `tests/mpi/` (`test_halo_ocean_mpi`, `test_ocean_dyn_mpi`, `test_ocean_bt_cfl_mpi`) and by the decomposition-invariance gate (`tests/regression/decomp_invariance.py`). That gate exists because two real bugs were invisible at one rank *by construction*.

Remaining: strong-scaling characterisation, and the features still flagged single-rank (see §3).

</details>

<details>
<summary>Isopycnal coordinate and mesoscale closures, detail</summary>

- **RHO:** `VCOORD_RHO` locates each target-density surface per column (PPM density reconstruction plus Newton root-find) and remaps the interfaces onto it. `VCOORD_HYCOM` adds the Bleck hybrid (isopycnal interior plus a z\*-like surface floor). The `VCOORD_LAGRANGIAN` endpoint still drifts off-isopycnal under diabatic forcing by design; RHO is the re-pinning path.
- **GM / Redi / VarMix / MEKE:** isopycnal slopes (Griffies 1998) feed GM thickness diffusion (bolus) and Redi neutral diffusion (continuous sweep plus GPU flux). VarMix gives Res_fn and Visbeck/Eady coefficients. Prognostic MEKE closes the GM-MEKE loop. **Caveat:** GM is a *coarse*-resolution parameterization. At the 1-2 km eddy-resolving target run it *off*; neutral diffusion stays useful against z\* spurious diapycnal mixing.

</details>

<details>
<summary>Non-hydrostatic on the C-grid, history</summary>

This was scoped as an *adaptation of the proven coastal sigma-coordinate NH solver*. That implementation, and the elliptic-solver machinery it shared with the semi-implicit free surface, left with the coastal carve-out. The item stands, but the reference implementation is no longer in this tree. Build it from the published method, or against the coastal repository. The coastal NH solver is no longer available as the in-tree reference for validation.

</details>

<details>
<summary>Sea ice, detail</summary>

No longer tentative; the ladder was built. `src/core/ice/` ships the Winton enthalpy column thermodynamics, multi-category ITD, ITD-aware transport, C-grid EVP dynamics, frazil, and the ice-ocean salt/heat/stress couplers, all under `&ocean_ice_nml` (default off, byte-identical). Remaining: ridging, the ice-shell / shelf-cavity interface kernel, and a polar validation campaign.

</details>

<details>
<summary>Eady / ACC channel, detail</summary>

PBC shipped. Eady runs as a periodic-x channel (`validation_examples/ocean/eady/`) and, since 2026-09-13, reproduces the linear-theory growth rate of its channel mode (1.93e-6 vs 1.98e-6 1/s; gated in `tests/regression/stability_manifest.py`). Remaining: an ACC transport run. The front's nonlinear collapse past day about 65 is not integrable at dt = 600 s (open; see the nml header).

</details>

<details>
<summary>Outer split scheme (§1) and the <code>ssp_rk2</code> residual (§2), detail</summary>

- **Shipped:** `&ocean_bt_nml split_scheme` selects `pred_corr` (default; predictor-corrector, tendencies on the barotropic-window means `u_av`/`h_av`, forward-backward gravity-wave pairing) or `ssp_rk2` (Heun, experimental).
- The Coriolis reference the fast loop subtracts is evaluated on the same velocity the slow tendency used, under the same weights, so the two cancel exactly (MOM6 `ubt_Cor`). Evaluating it at stage entry instead pumped a barotropic seiche.
- Both vertical-friction tridiagonals gained an `nz = 1` path. `scratch_3d` buffers are attached holding their host payload (zero from `init`) rather than as allocator leftovers.
- **Open (`ssp_rk2`):** Heun amplifies oscillatory modes at `(ω·Δt)⁴`, so a stratified rest state manufactures internal-wave energy. It is gated by `validation_examples/ocean/eady/resting_stratified_channel.nml`, which every scheme must pass on magnitude. The *settle* criterion is a scoped expected failure. Options: retire `ssp_rk2`, or give the outer loop a neutral or dissipative option (forward-backward is neutral; an RK3 outer step damps at the same order). *Validate:* a motionless stratified periodic channel does not grow.

</details>

<details>
<summary>Cross-rank reproducibility, history</summary>

The original item read: "`En`/`MaxCFL` diverge across rank counts while mass stays bit-identical; non-reproducible FP reductions are the suspect. EFP reproducing sums ship (`reproducing_sums`, default off). *Validate:* identical time series at 104 vs 208 ranks."

Since then `reproducing_sums` became the default (`ocean_diag_config_t` in `src/core/rdb_config.F90`), and `tests/mpi/test_ocean_decomp_bitid_mpi.F90` gates bitwise equality of the full restart registry across every `px x py` factorisation at 1, 2 and 4 ranks. The 104 vs 208 rank comparison is tracked as an open item in §3 (Regional enabling).

</details>

<details>
<summary>Status changes in this revision (2026-10-12)</summary>

- Cross-rank bitwise reproducibility: open to **done**.
- Tripolar fold and sea ice: removed from the single-rank list (sea ice under the fold stays single-rank).
- Porous barriers: removed from the single-rank list.
- OBCs: relabelled from "finish NESTED/CLAMPED/SPONGE dispatch" to "file-backed nesting data". CLAMPED, SPONGE and per-tracer inflow values ship.
- Partial cells: `z_fixed` partial bed cells recorded as shipped; general cut cells stay open.
- Phase 5g: "six VCOORD families" corrected to ten.
- Added: ACCESS-OM3 parity items C1, C7 part 1, C3, C2a, C2b, and the October 2026 code-health refactors.

</details>
