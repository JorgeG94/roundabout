# Roundabout — GPU-Native Ocean Solver

I had to learn about ocean dynamics so I wrote a code to do it. It runs on CPUs and GPUs
with `do concurrent` being the way to access parallelism. MPI is availbale too.

It is also a way to find bugs on Intel, AMD, and NVIDIA and check for LFortran portability.

> **Where things are.** Acronyms → [`docs/codebase/CONCEPTS.md`](docs/codebase/CONCEPTS.md) · Which closure/scheme is enabled + the tunable knobs → [`docs/CLOSURE_MATRIX.md`](docs/CLOSURE_MATRIX.md) · Design contract (god-state slot map) → [`src/core/ocean/README.md`](src/core/ocean/README.md) · Capabilities + limits → [`docs/CAPABILITIES_AND_LIMITATIONS.md`](docs/CAPABILITIES_AND_LIMITATIONS.md) · Physics + namelist → [`docs/REFERENCE.md`](docs/REFERENCE.md) · Style + GPU (`do concurrent` + OpenACC) programming → [`FORTRAN_STYLE.md`](FORTRAN_STYLE.md) · Repo map + run lifecycle → [`docs/codebase/INDEX.md`](docs/codebase/INDEX.md).
>
> **Per-module / per-procedure docs come from the source `!!` docstrings — rendered by FORD (built on PR to `main`). Don't hand-write parallel markdown describing what a module/kernel does; that's the copy that rots.**

## Build

```bash
module load cmake nvhpc
cmake -B build -S . && cmake --build build
cd build && ctest --output-on-failure
```

Key CMake options: `RDB_ENABLE_GPU` (default **OFF** — opt-in; enable with `-DRDB_ENABLE_GPU=ON` on the NVHPC toolchain. An explicit ON with a compiler lacking a GPU path is a configure-time error, never a silent CPU downgrade; the value prints on the `GPU offload:` configure line), `RDB_ENABLE_MPI` (default OFF), `RDB_CUDA_AWARE_MPI`, `RDB_GPU_ARCH` (cc70/80/90), `RDB_ENABLE_DOUBLE` (default ON). Full table in `cmake/options.cmake`.

`./rdb --validate-only case.nml` runs every configure-time check — `validate_config` (all failures, not just the first) and the whole `engine_setup` sequence — and exits **0** (accepted) or **3** (refused, each reason logged) without mapping the device, stepping or writing anything. It is how the compatibility matrix (`tests/regression/compat_matrix.py`) asks the model what it refuses.

## Local development conventions (agents + contributors)

Two standing rules for anything run against this repo — humans and AI agents alike:

- **Never `pip install` (or otherwise install packages).** Use only what the
  toolchain already provides (`mac_env.sh` / the loaded module environment). If a
  Python helper needs a package that isn't present, rewrite it against the
  standard library — a missing dependency is a signal to simplify, not to install.
- **All local/throwaway artifacts go in `tmp_local_artifacts/` at the repo root** —
  prototypes, scratch scripts, generated data, plots, logs, profiler/benchmark
  dumps, intermediate command output. It is git-ignored, so nothing leaks into a
  commit, and every worktree + every agent writes to the same known place. Do
  **NOT** write to `/tmp`, `/private/tmp`, the session scratchpad, or anywhere
  else in the tree — `tmp_local_artifacts/` is the *only* place. (Committed
  tooling under `tools/` is the standing exception — that's real code, not an
  artifact.)
- **Build directories are named `build_<toolchain>`** (`build_serial`,
  `build_cc90`, ...) — `.gitignore` covers `build*/` and `*_build/`; a build
  tree outside those patterns is how a stray `git add -A` once poisoned the
  object store.  Never name a build dir anything else.
- **Local builds are non-MPI.** Configure with `RDB_ENABLE_MPI=OFF` (the
  default) — do NOT enable MPI for local build/test. The mac dev box is
  single-rank and the sea-ice path is single-rank by design; MPI is a
  cluster-only build concern.

## Pre-commit workflow

Before each commit, in order:

1. **`pre-commit run --all`** — formatting (whitespace, trailing newlines, end-of-file).
2. **Fortitude lint** — Fortran static checks (`fortitude check` or via the pre-commit hook).
3. **`ctest`** — full test suite. **NEVER pass `-j N` when the GPU build is active** — every worker shares one GPU, so parallel test execution causes spurious failures + hangs. Building with `-j` is fine; only `ctest -j N` is the problem. Use `ctest -R rdb` (the rdb-only subset, 177 tests, ~35 s) for fast regression checks — note the FULL suite additionally carries ~26 pic/test-drive dependency self-tests, so "N/N" reports should say which scope they ran.
4. **Validate the docs — always, before merging anything.** FORD owns per-module / per-procedure docs from the source `!!` docstrings (don't hand-write those). But the hand-maintained **synthesis** docs rot silently and a stale cell is a bug: re-read [`docs/CLOSURE_MATRIX.md`](docs/CLOSURE_MATRIX.md) (drift-checked by `tools/check_closure_matrix.py`), the design-contract README (`src/core/ocean/README.md`), the `docs/howto/` extension guides, and [`docs/CAPABILITIES_AND_LIMITATIONS.md`](docs/CAPABILITIES_AND_LIMITATIONS.md). If your change adds, ports, removes, or renames a closure / scheme / knob / extension seam, update the relevant doc **in the same PR** — and when a doc and the code disagree, the code is authority and the doc is the bug to fix.

When adding a new capability (typical pattern for ocean-physics work):

1. **Namelist knob** with default = off (preserves bit-identity for existing nmls + tests).
2. **Kernel implementation** gated on the knob.
3. **Unit test** covering the new code path (analytical tests are high-leverage — every one added to the ocean core has caught a real bug).
4. **`pre-commit run --all` + fortitude** — pass before committing.
5. **`ctest`** (no `-j`) green on the GPU build.
6. **One commit per capability; one PR per capability.** Smaller PRs review faster and bisect cleanly.

**Compatibility matrix.** A change that adds, lifts or moves a configure-time refusal, adds an axis value (a closure, scheme or coordinate), or fixes a composition crash runs `python3 tests/regression/compat_matrix.py run --build-dir build_gfortran --jobs 4` (≈1 min on CPU; pairwise over every closure / coordinate / edge axis on a 24x16x10 domain). A cell that passes then runs the LEGS: ENERGY (En against the PASS population), RESTART (12 + engine-path restart + 12 == 24, bitwise) and, with `--mpi-build-dir`, DECOMP (1 rank == 2x2 == 4x1, owned cells, bitwise); the GPU leg adds CROSS_BACKEND against the CPU report (`--cells-from` + `--reference`; `tests/regression/compat_nightly.sh` runs both). `--diff-base origin/main` runs only the cells a diff touches. An unexplained refusal, crash, non-finite, budget miss or leg failure is a FAIL; a `KNOWN_GAP` row in `tests/regression/compat_expect.py` that now passes is an XPASS and also fails — delete the row in the same PR. See `tests/regression/README.md`, "The pairwise compatibility matrix".

## Project Layout

```
src/
  core/
    rdb_{constants,config,decomp,profiler,array_utils}.F90   infra
    rdb_grid, rdb_tracer, rdb_{barotropic,multilayer}_state,
    rdb_barotropic_workstate                                 shared state types
    rdb_state      register_default_tracers only — the A-grid god state is gone
    ocean/         the C-grid dyn-core — read src/core/ocean/README.md first
      state/ dynamics/split_rk2/ vcoord/ boundary/ forcing/ diag/ io/
      kernels/{barotropic,continuity_ppm,coriolis_adv,vertical_advection}/
    ice/           sea ice — thermo/, itd/, transport/, dynamics/, coupling/, state/
  # shared physics
  equation_of_state/   rdb_eos (+ rdb_ocean_eos_compute shim in core/ocean/state/)
  pressure_force/      rdb_ocean_pgf_reconstruct + rdb_ocean_pressure_force
  tracer/              rdb_recon_weno + rdb_ocean_{ideal_age,pseudo_salt}
  ALE/                            rdb_remap_column + rdb_ocean_remap + rdb_vcoord
  shared_module_utilities/  *.inc — small `pure`/`!$acc routine seq` helper BODIES
                   `#include`d into several modules' `contains` so each gets a
                   local, inlinable copy of ONE source text (see its README)
  parameterizations/
    vertical/      vmix, vdiff, epbl, kappa_shear, tidal_mixing, bottom_drag,
                   surface_flux + surface_stress, wave_speed, geothermal
    lateral/       horizontal_viscosity, hdiff_tracer, lateral_mix, meke,
                   gm, redi, varmix, mle, isopycnal_slopes
  framework/     god-state helpers: scratch_3d, safe_math, efp, mem_report
  driver/        rdb_driver — driver_run (the compute-rank lifecycle)
  io/            banner, bathymetry, io_netcdf (+ output_rank_filename)
  comm/          ONE comm implementation (no stub twin) — always over pic_mpi_lib; see comm/README.md
app/main.F90     Thin entry — calls driver_run
tests/           test_*.F90 (unit + analytical, flat) + mpi/
benchmarks/      bench_ocean (namelist-driven, NetCDF-free dyn-core throughput)
validation_examples/ocean/    canonical / analytical benchmarks
```

## Architecture

### One regime (`sim_type='ocean'`)

`sim_type='ocean'` is the only accepted value (`validate_config` fails loud on
anything else) — Arakawa C-grid + continuity-PPM + PV-conserving Coriolis +
split-explicit RK2, targeting 1 km submesoscale-resolving regional and global
hydrostatic configurations. **Tier-1 dyn-core operational** (2026-05-24) — the
canonical MOM6 double-gyre reference
(`validation_examples/ocean/double_gyre/double_gyre_mom6.nml`) runs stable to
day 580 with the production envelope (`nu_h=10000` + Smagorinsky + distributed
linear drag, `VCOORD_SIGMA`, `dt=1200s`). **Read `src/core/ocean/README.md`
first** — it's the slot map + design contract (4 rules, the `is_init` /
`scratch_3d_buffer_t` / OpenMP-only / no-direct-MPI conventions) + "how to pick
up a slot" recipe. The executed phase-map + long-horizon backlog live in
`docs/ROADMAP_OCEAN.md`.

`ocean_state_t` composes the ocean slot types + the `framework/` modules
(scratch buffer, safe-math). The dynamical-core slots are filled (continuity,
Coriolis-adv, EOS, PGF, vmix, hvisc/hdiff, bottom drag, ALE remap, diag
manager); deferred work is mostly Phase-5b+ refinements (Hollingsworth-Källén
correction, NH-on-ocean, OBC dispatch). See the "Ocean dyn-core" subsection
below for the current shipped surface.

The `no-mpi-in-rdb` pre-commit hook enforces "all MPI calls go through
`pic_mpi_lib`".

### Backend

Structured Cartesian only. SoA `(nx,ny[,nz])` for coalesced GPU access. The
unstructured (triangular, KNP central-upwind) backend left with the coastal
split — there is no `*_unstr` code path, no mesh reader, and no Hilbert
partitioner in this tree.

### GPU Parallelism

All data-parallel loops use `do concurrent`. GPU offload is NVHPC `-stdpar=gpu` plus a thin layer of OpenACC directives for what `do concurrent` alone can't express: reductions use `!$acc parallel loop reduction(...)` (inert comment on non-OpenACC compilers), and `!$acc enter/exit data` / `!$acc update` is used directly — no `#ifdef` wrapping. State arrays stay device-resident; only diagnostic output triggers D→H transfers. `RDB_GPU_OFFLOAD` is reserved for CUDA-aware MPI selection + multi-GPU build-config guard. There is **no full-OpenMP port** — `do concurrent` is never rewritten to an `!$omp` worksharing construct. For compilers without an OpenACC path, `tools/acc_to_omp.py` auto-translates the `!$acc` directives to `!$omp` equivalents (`do concurrent` untouched); the result is regenerated and CI-verified (build + `ctest -R rdb` under gfortran `-fopenmp`) on every push to `main` onto branch `auto/dc-openmp` (`.github/workflows/sync-dc-openmp.yml`).

### Vertical Layer Convention (load-bearing)

The layer stack is **bottom-up (ROMS-style)**: `k=1` is the bed, `k=nz` is the surface layer. Surface forcings land at `(:,:,nz)` (`k_top` under a `z_fixed` rigid top); bed forcings at `(:,:,1)` — precisely, at the first LIVE layer counting up, `multilayer_state_t%k_bot` (`≡ 1` off `z_fixed`). Any new layered kernel must follow this — PGF, rho_ref, per-layer drag, vcoord output, surface tracer fluxes, vertical advection/diffusion all assume it. Regression gates: `test_ocean_conservation_salt_heat`, `test_ocean_sw_penetration` (positive Q warms `k=nz` and preserves stratification), `test_ocean_vcoord*`.

### Ocean dyn-core (`sim_type='ocean'`, structured C-grid)

The dynamical core, designed for basin- and global-scale hydrostatic regional
ocean. Lives at
`src/core/ocean/` + `src/parameterizations/vertical/`.

**Prognostics** — Arakawa C-grid `(η, u_face_x, v_face_y, h_layer,
hu/hv_layer, S, T, ...)`. Continuity is a transport equation
(`∂h_total/∂t = -∇·(hu)`) solved with **continuity-PPM** — the layer
thickness IS the prognostic variable, not a constraint enforced by a
Poisson solver. Wright (1997) nonlinear EOS is the production path;
linear EOS available for validation.

**Time integration** — split-explicit, outer scheme selected by
`&ocean_bt_nml split_scheme`. Both schemes ship, both are under test, and
they differ in what they cost you.

- **`"pred_corr"` — the DEFAULT** (since 2026-09-14), and the MOM6
  predictor-corrector. Off-centred predictor at `pc_be·dt`, slow tendencies on
  the `u_av`/`h_av` step time-means, ONE prognostic update in the corrector,
  forward-backward gravity-wave pairing — neutrally stable to `ω·dt = 2`,
  which lifts the internal-wave `dt` ceiling and removes the resting-state
  growth below. It preserves the Eady benchmark's physical mode (slightly
  stronger: max|v| ×2074 over 60 days vs ×1480), and costs ~7 % more wall time
  per step than `ssp_rk2` (`coriolis_coast`, 20 simulated days: 70.3 s vs
  65.9 s on one V100) — not a different order. `validate_config` refuses it
  **fail-loud** outside its v1 envelope (`eulerian_z`, wet/dry,
  `dt_tracer_advect_ratio > 1`); six shipped namelists pin `ssp_rk2` for that
  reason. It shipped as `"mom6_pc"`, then briefly as `"split_rk2"` — which
  collided with `dynamics/split_rk2/`, the directory holding the outer-loop
  machinery of BOTH schemes. Both dead spellings now fail loud naming
  `pred_corr`. **Residual limit:** it does not eliminate the
  resting-state growth, it slows it ~36× (83-day e-folding), so
  `resting_stratified_channel` stays a scoped XFAIL on the *settle* gate.
- **`"ssp_rk2"` — EXPERIMENTAL.** Two identical stages + SSP average. Widest
  envelope: every vcoord including `eulerian_z`, wet/dry, and
  `dt_tracer_advect_ratio > 1`; the only scheme wired through the windowed
  tracer-advection path. **Fully supported and fully tested** — experimental
  labels the ANSWER, not the code path. **What it costs, measured:** the
  two-stage average amplifies an internal gravity wave by `√(1+(ω·dt)⁴/4)` per
  step, so it *manufactures* energy from a motionless stratified state. On
  `validation_examples/ocean/eady/resting_stratified_channel.nml` (flat bed,
  periodic, stably stratified, at rest, ±0.5 mK seed, **no energy source**)
  En reaches **2.992E-05 m²/s²** by day 25 — 7.7 mm/s of current out of
  nothing — still climbing on a 2.5-day e-folding, where `pred_corr` holds
  **1.739E-09** (17 000×, 83-day e-folding). The outer split is the cause and
  nothing else: the Coriolis form, the ALE remap and the PGF form were each
  substituted and each moved the answer < 0.1 % — all three **exonerated**;
  `dT/dz = 0` dropped En 119×; removing the lateral viscosity *raised* it. It
  is a `(ω·dt)⁴` noise floor — a forced, energetic, viscous run sits decades
  above it and never notices; a quiescent, weakly-damped or long-spin-up run
  does not, and there the manufactured energy IS the signal. Carried as a
  scoped XFAIL on `resting_stratified_channel__ssp_rk2`, not as institutional
  memory.

The stability suite runs an `ssp_rk2` twin of every case whose namelist does
not pin a scheme, so neither branch can rot. Both schemes wrap the same
nonlinear barotropic fast loop (forward-backward Euler substeps for η +
barotropic u/v + ζ + KE) + BT correction back into the layers.
**MOM6-parity defaults (2026-09-22):** the fast loop damps its gravity waves
with MOM6's `&ocean_bt_nml bebt = 0.1` (was 0 — pure, neutral forward-backward,
under which nothing damped the barotropic grid-scale mode `pred_corr` leaves
alone), and the continuity's `uhbt` renormalisation uses the continuous
donor-flip flux (`&ocean_continuity_nml renorm_consistent_flux = .true.`,
MOM6 `zonal_flux_adjust`; bit-identical wherever no donor flips) — together
they close FINDING B of the vertical-coordinate matrix. Namelists that set
`bebt` explicitly (mostly 0.2) keep their value.
`auto_n_inner=.true.`
derives `n_inner` from the gravity-wave CFL **once at configure time**
(`configure_ocean_bt_split` in `rdb_ocean_setup.F90`), writing the
resolved value back into `cfg%ocean%bt%n_inner` — which is why `cfg` is
`intent(inout)` all the way up to `main`. It is NOT re-derived per step;
CFL truncation is a counter (`dyn%ntrunc_total`), not a controller.
Under `ssp_rk2` the two stages are averaged; under `pred_corr` the
corrector's single full-dt update IS the step (averaging forward-backward
stages annihilates the internal wave), and the salt/heat console budget weight
follows the scheme (`ocean_budget_stage_weight`: 0.5 vs 1.0). The fast loop's
Coriolis/advection reference (`bt_work%cor_ref_u/v`, MOM6 `ubt_Cor`) is
likewise scheme-dependent and MUST match the velocity the slow Coriolis was
evaluated on — see `set_cor_ref_velocity`, and
`tests/test_ocean_cor_ref_seiche.F90`, the rotating closed-basin non-growth
test that guards it.

**Physics shipped** — the per-closure detail (knobs, formulae, per-path test
names) lives in [`docs/CLOSURE_MATRIX.md`](docs/CLOSURE_MATRIX.md) (drift-checked
by `tools/check_closure_matrix.py`); the knob→default map in
`docs/generated_nml_knobs.md`; narrative treatment in `rdb_docs/`. In brief:
- **PGF** (`&ocean_pgf_nml form=`) — Montgomery (**default**) / FV-lite / FV-Wright / gprime (2-layer) / FV-MOM6; optional in-layer PLM/PPM T/S reconstruction (Boole quadrature).
- **Coriolis** (`&ocean_coriolis_nml form=`) — Sadourny enstrophy PV-flux (default, velocity form), energy transport form (`sadourny_energy`), HK correction (`sadourny_hk`). Orthogonal `pv_adv_scheme=` selects the corner-vorticity→face interpolation in the Sadourny path: `centered` (default, 2-pt average ⇒ bit-identical) or `weno3`/`weno5`/`weno7` (upwind-biased WENO-Z reconstruction, MOM6 WENOVI{3,5,7}TH — sharpens submesoscale PV fronts; scale-selective dissipation that controls grid-scale vorticity noise where explicit viscosity is low). weno5/weno7 (radius 3/4 stencils) require `nghost≥3`/`4` (fail-loud via `pv_adv_required_nghost`).
- **Lateral closures** (`&ocean_hvisc_nml`, `&ocean_meke_nml`) — Leith, Smagorinsky KH+AH, Leith-biharmonic, constant floors + per-cell CFL clamps, resolution-scaled visc, anisotropic + live velocity-scale, MEKE backscatter. `stress_tensor` composes with the biharmonic add-on (`nu_4`/`smag_ah`/`leith_biharm`), matching MOM6 — it no longer disables it, and `kh_aniso` + biharmonic is a legal combination. Dispatcher is fail-loud; biharmonic backstop with a NON-ZERO dissipation coefficient (not just a flow-aware closure selected) mandatory under backscatter.
- **Vertical mixing** — PP81 interior + KPP boundary overlay (default on); every closure contributes into `kv`/`kt` — `ks` is derived from `kt` by `vmix_split_kd_heat_salt` (last statement before the gate; MOM6 `Kd_salt = Kd_int + Kd_extra_S` / `Kd_heat = Kd_int + Kd_extra_T` — `ks ≡ kt` unless `&ocean_ddiff_nml` double diffusion is on, which folds the asymmetric salt-fingering / diffusive-convection `Kd_extra` into the split) — with `vmix_assemble` the single floors/ceilings/smoothing gate over `kv`/`kt`/`ks`. `vdiff_apply_tracers` takes `kt_source` (temperature) and `ks_source` (salinity + every passive tracer).
- **EPBL** (`&ocean_epbl_nml`) — Reichl-Hallberg energetics PBL + Langmuir; mutually exclusive with KPP.
- **kappa-shear** (`&ocean_kappa_shear_nml`) — JHL08 prognostic interior shear turbulence, additive interior merge; opt-in vertex form (`at_vertex`, MOM6 VERTEX_SHEAR/OM5 — corner solve on native face velocities, corner Kd averaged back to centres, arithmetic or geometric+`kdmin`; Kv routed corner→face into the momentum vdiff (MOM6 `Kv_shear_Bu`), cell-centred kv merge suppressed in vertex mode).
- **tidal mixing** (`&ocean_tidal_mixing_nml`) — St-Laurent/Simmons bottom-intensified diapycnal diffusivity.
- **convective adjustment** (`&ocean_conv_nml`) — Brunt-Väisälä trigger (`N² < n2_thresh`), `kt`/`kv` raised to `kd_conv`/`prandtl_conv·kd_conv` below the active KPP/EPBL boundary layer; `max()` contributor, writes `kv`/`kt` only (never `ks`).
- **background mixing** (`&ocean_vmix_nml bkgnd_profile` **xor** `bkgnd_henyey`) — two MUTUALLY EXCLUSIVE replacements for the scalar `kt_bg`/`ks_bg`/`kv_bg` floor (enabling both fails loud at configure, matching MOM6's one-background-scheme rule): Bryan-Lewis (1979) depth-varying `kd_bg` field, or the Henyey (1986) latitude factor — constant-`N0` simplification of Harrison & Hallberg (2008), latitude-only, no per-step N dependency — scaling the SCALAR tracer floors, `max(bkgnd_kd_min, kt_bg·L(φ))` (MOM6 `KD_MIN`, negative ⇒ `0.01·kt_bg`; fail-loud on a cartesian grid, where `geolatT ≡ 0` makes every column equatorial); both default off ⇒ bit-identical.
- **tides** — body forcing (`&ocean_tides_nml`, equilibrium `η_eq`) + scalar SAL + boundary-tide nodal correction (`&ocean_bc_nml obc_tidal_nodal`), all off the shared astro generator.
- **Barotropic linear wave drag** (`&ocean_bt_nml wave_drag`, Egbert & Ray 2001; Jayne & St Laurent 2001) — static per-face piston velocity `r_H` MULTIPLIED into `bt_rem_u/v` inside the BT substep (`form="uniform"` or a resolved-bathymetry-variance `"roughness_proxy"` placeholder for the real subgrid `⟨h²⟩`; `"file"` fails loud pending PR-14); composes with `substep_drag`.
- **Porous barriers** (`&ocean_porous_nml`, Adcroft 2013) — subgrid sill/strait blocking: a per-layer OPEN-AREA fraction from the three-parameter along-face `d_min`/`d_max`/`d_avg` fit narrows `dy_cu`/`dx_cv` in the continuity-PPM and transport-Coriolis mass fluxes, and the column-integrated fraction narrows the BAROTROPIC substep widths (`dy_cu_bt`/`dx_cv_bt`) so the BT solve isn't porous-blind; recomputed once per outer step (MOM6 cadence). Multi-rank bit-identical (`test_ocean_decomp_bitid_mpi`); fails loud with `&ocean_bt_nml bt_halo > 0` (the wide-halo BT clone carries no porous stats) and `&ocean_wetdry_nml enable`. `source="resolved"` is a documented, wet-gated RESOLVED-bathymetry proxy: it only sees variation ALONG the face, so a degenerate statistic (bathymetry uniform along the face) is left fully open rather than walled at the two-cell mean depth (true subgrid needs an offline `topog_edge.nc`; `"file"` fails loud). `eta_interp="max"` (default) blocks the LEAST. CHANNEL_DRAG form drag is a separate follow-up.
- **Bottom drag** (`&ocean_bdrag_nml`) — linear / quadratic (default), each with an HBBL-distributed mode; every form sits on the face's first live layer `k_bot_u/v` (`k=1` off `z_fixed`).
- **Ice-shelf top drag** (`&ocean_tdrag_nml`, default off) — the mirror of the bottom drag at `k=nz` on ice-covered FACES (a face is under ice if EITHER abutting cell is, so the calving-front face is dragged): quadratic (ISOMIP+ `C_d=2.5e-3`) / linear, HTBL-distributed mode, optional backward-Euler form. Requires `&ocean_cavity_dyn_nml`; shares ONE `C_d` with `&ocean_cavity_melt_nml cdrag_top` (disagreement fails loud). Reaches the barotropic mode through `F_slow` like the bottom drag. Publishes `stress_top`, which the RK2 stage drivers copy inline into `surface_stress%stress_shelf` — the under-ice `u_*` source BOTH boundary-layer schemes now read (Phase 4b: `u_*^2 = (stress_mag + stress_shelf)/ρ₀`; same stage, no lag). With melt on but this group off, `engine_step_finalize` fills `stress_shelf` from the melt slot's own `u_*` instead (`ρ₀·u_*²`, one thermo step lagged).
- **Implicit stress/drag fold** (`&ocean_vdiff_nml`) — folds wind stress + bottom drag (bed-row `k_bot` diagonal, `k=1` off `z_fixed`) + ice-shelf top drag (`implicit_top_drag`, `k=nz` diagonal, and masks the wind RHS under cover) into the backward-Euler vertical-friction solve (kills thin-layer CFL blow-up).
- **Shortwave penetration** (`&ocean_thermo_nml sw_pen_frac`; `sw_source=net_heat|q_sw`; two-band `sw_transmission` shared by the deposition kernel + both BL schemes) — with MOM6 `KPP_SHORTWAVE_METHOD` `B_0` coupling (`kpp_sw_method`) and the EPBL penetrating-SW TKE ledger (`epbl_sw_ctke`, `Phi(tau)` in-layer PE cost); all inert at `sw_pen_frac=0` — + **surface buoyancy restoring** (`&ocean_restore_nml`).
- **ALE remap** — PPM, momentum + tracers, thermo-cadence gated (`is_thermo_step()`).
- **Windowed tracer advection** (`&ocean_vmix_nml dt_tracer_advect_ratio`) — CW-PPM drain decoupled from per-step dynamics; single-rank.
- **Diagnostics** — registry + cadence dispatch (device-side fill, H←D on the cadence fire); unified `&ocean_diag_nml diags` selection knob; per-diag output vcoord remap (+ `ice_conc`/`ice_thick` when `&ocean_ice_nml enable`; mean ice conc/thick also on the console status line; derived catalog also carries `ice_speed`/`ice_u`/`ice_v`).
- **Boundaries** (`&ocean_bc_nml`) — WALL/OPEN(Flather)/TIDAL/CLAMPED/CHAPMAN/SPONGE/PERIODIC + Orlanski radiation + asymmetric nudging.
- **Interior land masking** (static free-slip walls via metric-zeroing) + **dynamic wet/dry** (`&ocean_wetdry_nml`; sigma only, single-rank, positive-definite outflow limiter).
- **Tracer registry** — `multilayer_state_t%register_passive_tracer(grid, name, units, long_name, idx)` grows the S+T[+age] registry at setup (6-arg; `idx=0` on refusal; `registry_locked` after `enter_data` closes the mem:separate device-map foot-gun); every passive-transport kernel already loops it, so a new tracer rides for free. Budget attribution is `tracer_t%budget_id` (`NONE`/`HEAT`/`SALT`), not an index comparison. Shipped packages: ideal age + **pseudo-salt** (`&ocean_tracers_nml enable_pseudo_salt`, Shao 2016 verification tracer — seeded to S, given S's surface salt flux + KPP nonlocal mirror; deviation measures the passive-vs-active transport-path error). Recipe: `docs/howto/add_passive_tracer.md`.

Most new-physics knobs default off ⇒ bit-identical (KPP + ALE-remap are the on-by-default exceptions); each path is guarded by an analytical test.

**Horizontal grids** — `&ocean_grid_nml grid_config`: `cartesian`
(default, bit-identity), `spherical` lon-lat sector, `supergrid`
(MOM6 mosaic reader), and `tripolar` — Murray (1996) bipolar Arctic
cap above `phi_join` + ordinary lon-lat below, closed by a north fold
(`north="tripolar_fold"`, requires periodic west/east) that is
single-rank IN X: north-south splits (`px = 1`, any `py`) run
bit-identical to one rank — only the north-edge rank folds
(`bc%north_fold` is rank-local) — while `px > 1` is refused at configure
until the distributed fold exchange lands.
Kernels consume full 2D metric arrays only (`ocean_metrics_t` slot);
the fold exchange (`rdb_ocean_fold` + `rdb_ocean_fold_apply`)
reverses-i and sign-flips vector normals, projecting the
duplicated-DOF fold-line row antisymmetric -- storage row `ng+nj+1`, because
roundabout's `v` sits on the SOUTH face of its cell and corners at the SW
(MOM6's north-face rule is one row off here; see `rdb_ocean_fold` header).

**Vertical coords** — all ten `VCOORD_*` families dispatch through the
same ALE remap path; see the next subsection.

**Working envelope** (matches `validation_examples/ocean/double_gyre/double_gyre_mom6.nml`).
Ocean knobs live in per-concern sub-namelists — `&ocean_<group>_nml`
(coriolis, thermo, bt, pgf, bdrag, hvisc, vmix, continuity, topo, ic,
diag) with the `ocean_` prefix dropped from each key.  `tools/nml_split.py`
migrates a legacy `&ocean_setup_nml`:

```
&ocean_hvisc_nml      nu_h = 10000.0,  smag_ah = .true. /
&ocean_bdrag_nml      form = "linear", r = 2.5e-5, hbbl = 10.0, bg_vel = 0.1 /
&ocean_coriolis_nml   form = "sadourny" /
&ocean_pgf_nml        form = "gprime" /     ! 2-layer reduced gravity
&vcoord_nml           vcoord_type = "sigma" /
&time_nml             dt_fixed = 1200.0 /
```

580-day spinup stable; ~16 s wallclock for 30 days on a single V100.

**Deferred** — Hollingsworth-Källén Coriolis correction (PV form is
production, HK as a guard option),
non-hydrostatic on the C-grid, per-layer Orlanski
phase-speed radiation + file-backed boundary-data backends (Flather +
zero-gradient anomaly and the constant backend ship today), MPI
per-feature multi-rank support for the
single-rank closures (wet/dry, sea ice, the ice-shelf cavity, Chapman
edges, the tripolar fold split east-west (`px > 1`; north-south splits
ship), windowed tracer-advect drain), restarts across a rank-count change,
a bit-identical barotropic march-in — the C-grid MPI halo itself ships,
bit-identical to one rank.

### Vertical Coordinates

Dispatched via `vcoord_type` (namelist) → `parse_vcoord_type` → `VCOORD_*` enum.
**Ten** families ship; the enum is defined in `src/core/rdb_constants.F90` (the
authority) and the target-grid build lives in
`src/core/ocean/vcoord/rdb_ocean_vcoord.F90`. Eight are `select case` branches
of `ocean_vcoord_compute_target_h`; the two density-space coords need per-layer
T/S + the EOS and so come in through the sibling
`compute_target_h_rho` wrapper, dispatched from the same ALE remap driver
(`src/ALE/rdb_ocean_remap.F90`). The `compute_target_h` `case default` is
`error stop` — fail-loud, no silent fallback.

**The placement gate is `tests/test_ocean_vcoord_interface_depths.F90`** — it asserts the absolute geopotential DEPTH of every target interface, not just that thicknesses are positive and sum to the column (those pass even with the stack in the wrong half). Its third geometry is a **displaced column top** (a rigid lid / ice base): sigma-like families divide the live column proportionally and are correct; `Z_FIXED` and `ZSTAR_FULL` hang their stacks from the column top and are wrong by up to 500 m, pinned as `documents_*` cases the rigid-top slice must flip. Per-family table: `docs/CAPABILITIES_AND_LIMITATIONS.md`.

- `VCOORD_LAGRANGIAN` (-1) — pure Lagrangian / isopycnal; `target_h` is the live `h_layer` and the remap is a no-op (early return, no kernel launch). Parsed from `lagrangian` / `isopycnal`.
- `VCOORD_EULERIAN_Z` (0) — `H · dsig(k)`, η ignored. The ocean path's "leave the IC layers alone" default.
- `VCOORD_SIGMA` (1) — terrain-following, `(H + η) · dsig(k)`.
- `VCOORD_ZSIGMA` (2) — smoothstep blend sigma→fixed z-levels. **Refused at configure**: the deep branch reads `z_ref_global` as absolute depths in METRES, but its only writer is the dimensionless `z_ref_global(k) = k/nz` init, so every z-level interval is `1/nz` m and the whole column collapses into the bed layer — with `Σ target_h = H + η` still exact, which is why no conservation test caught it. Returns when the table is filled in metres (it is the natural seat for a sigma-near-the-top / z-below hybrid).
- `VCOORD_ZSTAR` (4) — **MOM6 z\*** (`REGRIDDING_COORDINATE_MODE = "Z*"`, MOM6 `build_zstar_column`; `ocean_vcoord_zstar_target`): the `z_fixed` nominal profile (`z_fixed_h_ref = max_depth`, uniform or `&vcoord_nml z_fixed_profile = "list"|"tanh"` — the same `z_fixed_zi` table, not a parallel one) DILATED per column by the free-surface stretching, over `z_fixed`'s partial bed cell and `zstar_h_min` bed fillers. Every layer's liveness is decided at `η = 0` (the `z_fixed` bed rule) and the live interfaces sit at `s·Z_k` below the free surface, `s = (H+η−n_f·h_min)/(H−n_f·h_min)` — MOM6's `(H+η)/H` with the filler stack kept out of the dilation — so the live/filler pattern is EXACTLY static for `η` of either sign (dilation keeps ratios; MOM6 coord_zlike.F90:109-143) and the `η = 0` target IS `z_fixed`'s bit for bit. Hence `zfixed_closed_faces` applies unchanged (static mask from the `η = 0` target, every consumer reused, `h_layer` seeded on the target under the knob or `&ocean_zinit_nml`; gate `test_ocean_zstar_closed_faces`). A partial bed cell scales with `s`, so it stays above `H_VANISHED` for any `η` that does not all but dry the column. Refused under an ice-shelf cavity (no rigid-top branch) and with wet/dry; `zstar_h_min` carries the INERT contract. **Until 2026-10 `"zstar"` was literally the `VCOORD_SIGMA` branch** (`case (VCOORD_SIGMA, VCOORD_ZSTAR)`), so every namelist that selects it changed answers with this coordinate (accepted as a bug fix). Aliases `zstar_lite`/`z-star` still parse to it.
- `VCOORD_ZSTAR_FULL` (5) — per-column `z_ref(0:nz)` from local bathymetry; surface layer anchored at `zstar_h_surf_target` regardless of H. Bed-side layers can vanish (`zstar_h_min`, which must stay ≤ `H_VANISHED` on this family and on `Z_FIXED` — a larger value fails loud at configure, because a filler above the marker is dynamically live while the coordinate still treats it as throwaway). Operators that divide by `h_layer` gate on `H_VANISHED = 1.5e-4 m`. The fine zone is `nz/3` layers of `zstar_h_surf_target` (`zstar_n_surf`/`zstar_stretching` are dead on the ocean path), the rest is terrain-following; a column shallower than the fine zone ends in a partial cell over fillers — the `z_fixed` staircase — so `&vcoord_nml zfixed_closed_faces` (name historical) is accepted here too: static mask from the `η = 0` target, every consumer reused, `h_layer` seeded on the target. It closes the filler faces only; the coarse zone's open faces keep the sigma PGF error. **Slated for retirement** in favour of `zstar` (MOM6 z\*) once that is proven on the 1-degree Southern Ocean and the global case — its coarse zone is terrain-following (50× `z_fixed`'s energy on the 1-degree SO, sigma PGF error at rx0 0.99), and its `η < 0` branch clips from the bed instead of dilating. **Caveat:** intertidal domains still leak 1-2% salt/cycle from wet/dry destruction. Recommend sigma for those.
- `VCOORD_ZSTAR_SIGMA` (6) — sigma in shallow water, z*-lite in deep. Conserves by construction. Consumes `z_ref_global` FRACTIONALLY (rescaled by `z_ref_global(nz)`), so it is immune to the units defect above — but with the uniform default table its deep branch is numerically indistinguishable from `VCOORD_SIGMA`, so the shipped namelists that select it are running sigma.
- `VCOORD_Z_FIXED` (7) — fixed-z interfaces from `z_fixed_h_ref` (uniform `max_depth/nz`, the default) or a STRETCHED nominal profile (`&vcoord_nml z_fixed_profile = "list"` with `z_fixed_dz`, surface first, or `"tanh"` from `z_fixed_dz_top`; `rdb_vcoord :: z_fixed_nominal_dz` → `vcoord%z_fixed_zi`/`z_fixed_dz`, read by the target builder and hence the closed-face mask, `k_top`, the partial top/bed cells and the on-target IC seed (cavity, or `&ocean_zinit_nml`) — gate `test_ocean_zfixed_stretched`), bed-side layers vanishing to `zstar_h_min` (which must stay ≤ `H_VANISHED` on this family — see `ZSTAR_FULL`) in shallow water (MOM6 `COORD_CONFIG="gprime"` layering as a per-step ALE target). Falls back to uniform sigma when the knob is unset. Parsed from `z_fixed` / `z_levels` / `gprime`. **The only z-like family legal under an ice-shelf cavity**: it reads `vcoord%z_top` (= `metrics%z_draft`, else 0), keeps its nominal interface depths geopotential, vanishes the layers that outcrop into the ice to the same inert filler and cuts the first live layer at the ice base into a partial top cell (Yung et al. 2026 Fig. 1b). `z_top = 0` ⇒ bit-identical. Under a cavity the top-side consumers read `multilayer_state_t%k_top` (the first LIVE layer counting down, `min`-rule face twins, filled once at configure from the same `η = 0` target the closed-face mask uses, `≡ nz` off a rigid top ⇒ bit-identical): melt heat/salt/mass, top drag (hence `stress_shelf` and the under-ice `u_*`), the implicit stress/drag folds and the momentum surface row, `mld_density`. The bed side has the mirror index on EVERY `z_fixed` run (not only under a cavity — every column shallower than the nominal stack has bed fillers): `multilayer_state_t%k_bot` (first LIVE layer counting UP, `max`-rule face twins `k_bot_u/v`, same `η = 0` target, built by `configure_ocean_k_bot`, face-halo-exchanged, static / not restart state, `≡ 1` elsewhere ⇒ bit-identical), read by bottom drag (bed-only, HBBL band, implicit-fold rate), the vdiff bed row (identity rows below it), geothermal, the tidal-mixing bed anchor, MEKE's bed speed and the bed-reaching shortwave; gate `test_ocean_zfixed_k_bot`. `implicit_drag` + `hbbl > 0` stays refused (the fold is one 2-D rate; a band needs a per-layer one). `&vcoord_nml zfixed_closed_faces = .false.` (the default) is REFUSED at configure over a stepped bed (any wet face whose two columns' bed lies in different nominal layers — the open-staircase PGF blows up in a few steps); flat / step-free beds stay accepted. Melt and top drag are accepted; `kappa_h /= 0` is accepted iff `&vcoord_nml zfixed_closed_faces`; KPP / EPBL / ideal age / GM / Redi / kappa-shear / tidal mixing / `regrid_time_scale` still fail loud; see `docs/CAPABILITIES_AND_LIMITATIONS.md`.
- `VCOORD_RHO` (8) — isopycnal: interfaces placed on prescribed potential-density surfaces `rho_target(0:nz)` (`rho_target_profile`: the `rho_target_light`→`dense` linspace, or `"list"` = `rho_target_list`) by inverting a PPM reconstruction of the column density. Validation-grade alone (weakly-stratified columns collapse). A column thinner than `nz·max(zstar_h_min, 2·H_VANISHED)` (land, dry) is not regridded (both density families).
- `VCOORD_HYCOM` (9) — hybrid z*/isopycnal (Bleck 2002, MOM6 `coord_hycom`): the same density inversion plus a bottom-up density monotonize before it and a z* nominal-floor sweep after it. Fixed-resolution near-surface z* band with an isopycnal interior — the production GVC coordinate. The floor is the z* coordinate resolution **in metres**, `z(k) ≥ Σ dz·(H+η)/H` clamped to the bed, read from the SAME `&vcoord_nml z_fixed_profile` table `z_fixed` uses (`"uniform"` = `max_depth/nz` m; `"list"`/`"tanh"` = `z_fixed_dz` / the tanh ramp) — MOM6 HYCOM1's `coordinateResolution`. It was `Σ dsig·(H+η)` with `dsig ≡ 1/nz` until 2026-10-02: a column fraction, i.e. a SIGMA floor (95 % of the 1° Southern Ocean's interfaces). Reuses `rho_target` / `rho_ref_pressure`; non-uniform targets via `rho_target_profile = "list"` + `rho_target_list` (both families).

(3 is unused — the enum has a gap, not a missing family.)

### Boundary Conditions

See the "Boundaries" bullet of the Ocean dyn-core section — `&ocean_bc_nml` selects WALL / OPEN (Flather) / TIDAL / CLAMPED / CHAPMAN / SPONGE / PERIODIC per edge, with Orlanski radiation and asymmetric nudging; the types live in `rdb_ocean_boundary_types` and the kernels in `src/core/ocean/boundary/`. Interior land masking (static free-slip walls via metric-zeroing) and dynamic wet/dry (`&ocean_wetdry_nml`) are separate.

### MPI

`rdb_decomp` splits the domain into rectangular subdomains. Halos via `mpi_f08` Isend/Irecv on `compute_comm`. Two-phase comm_env init: `comm_env_init()` → read config → `comm_env_setup_roles(.false.)` (there is no dedicated I/O-server rank any more — that was a coastal-path feature). Each GPU rank binds via `acc_set_device_num(gpu_rank = mod(node_rank, n_devices))`.

**Multi-GPU one node: pin `CUDA_VISIBLE_DEVICES` per rank BEFORE `MPI_Init`.** With `RDB_CUDA_AWARE_MPI=ON` (GPU-direct halo), UCX/hpcx creates a CUDA primary context at `MPI_Init` — which runs *before* `acc_set_device_num` — so with all GPUs visible every rank also lands a context on device 0 (`{0},{0,1},{0,2},{0,3}` in `nvidia-smi`, binding diagnostic still looks correct). Fix: each rank must see only its own GPU. The launcher (or the script, before the first CUDA call) sets `CUDA_VISIBLE_DEVICES` from the local-rank env var (`OMPI_COMM_WORLD_LOCAL_RANK`, …); the `mod(node_rank, n_devices)` clamp then binds the single visible device 0. With pinning, plain `mpirun -np N ./rdb …` works — no `bash -c` wrapper. Host-staged MPI (`RDB_CUDA_AWARE_MPI=OFF`) doesn't hit this (UCX never touches CUDA).

Halo coverage: `rdb_ocean_halo` + `ocean_halo_exchange_ml_state`. Tracer halo uses the outer-shim + flat-impl pattern.

Tripolar fold under MPI: every seam site runs exchange → periodic wrap → fold, and the fold is applied only by the rank that owns the north edge (`bc%north_fold` = tag `.and.` `has_north`, re-derived in `ocean_bc_state_set_edges`), which with `px = 1` holds the whole fold row, so the local fold kernels are exact. `px > 1` needs a distributed fold exchange (the mirror point `(ni+1-i, nj+1-d)` lives on another rank) and is refused at configure.

**A decomposed run IS the serial run — bitwise.** `tests/mpi/test_ocean_decomp_bitid_mpi` (ctest at 1/2/4 ranks) compares every restart-registry field, owned cells, bitwise, on every `px x py` factorisation (2x1, 1x2, 4x1, 2x2, 1x4) of seven configurations x both split schemes (walls + island, periodic z* + porous, tidal/Flather OBC, spherical, Orlanski/clamped/sponge OBC, the closure set, and the windowed file readers). Any mismatch is a bug. The recurring bug class it finds: a pass that writes a tile's PHYSICAL span only (an OBC fill, a sponge, a wall re-close) or keys a closure on the global edge TAG without `has_*` — at an MPI seam the neighbour's ghost copy goes stale or the seam is treated as a wall. Gate new edge/ghost writers on `has_*`, and refresh seam ghosts (collective, gated on rank-uniform tags) after any physical-span-only writer that something reads before the next exchange. `&mpi_nml px/py` unset on N ranks is auto-factored (north-south under a tripolar fold). The supergrid, bathymetry and z-level IC readers are windowed per rank (start/count). The console uses the EFP reproducing sums by default (`&ocean_diag_nml reproducing_sums`), so it prints the same digits on any rank count. The barotropic march-in (`&ocean_bt_nml bt_halo > 0`) is OPT-IN — AUTO resolves to 0 — because it is not bit-identical over variable bathymetry / with OBCs.

**Single-rank features fail loud at configure on >1 rank** (in `engine_setup`, keyed on the actual rank count): ice-shelf cavity, wet/dry, sea ice, Chapman OBC edges, `dt_tracer_advect_ratio > 1`, the tripolar fold with `px > 1` (or a north tile shorter than `nghost+1`), in-memory (staged) geometry injection, `nghost < 3`, an explicit `bt_halo` wider than the smallest tile. Restarts resume on the same decomposition only. Full list + the known non-bit-identical paths: `docs/CAPABILITIES_AND_LIMITATIONS.md`, *MPI (domain decomposition)*.

### I/O

Per-rank diag + restart files (`<prefix>_rank_NNNNNN.nc`) — no gather. Optional deflate (3-5x). Offline merge via `tools/merge_output.py`. Diagnostics run off the registry + cadence dispatch (device-side fill, H←D on the cadence fire) with per-diag output-vcoord remap.

## Gotchas

- **Vertical indexing**: see "Vertical Layer Convention" above. Surface = `nz`, bed = `1`. Don't flip.
- **Thin-layer constants of record (ocean path)**: `H_VANISHED` (1.5e-4, dynamic-vanish — skip/merge, don't clamp) and `H_DIV_EPS` (1e-20, pure 1/0 armour) live in `rdb_constants` with documented roles. New kernels pick the right one (D4 taxonomy in the docstrings); the massless-merge helper is `rdb_massless`.
- **The vanished-layer content rule is STRUCTURAL — don't re-derive it per kernel.** Tracers are CONTENT (`hTr = h·c`), so the tree carries one invariant: **I1′ — `h <= H_VANISHED` ⇒ `hTr = h·c_live` for every registered tracer**, `c_live` the concentration of the filler's DONOR (the nearest live layer above; the topmost live layer, `k_top`, for fillers inside an ice shelf; `hTr = 0` in a column with no live layer). It is restored by POOLING each filler run with its donor and mixing the pool to one concentration — content moves within the column ⇒ **nothing is recorded in any budget**. It replaced I1 (`hTr = 0`), which kept a filler's mass but not its concentration, so the thickness continuity moved out of a filler arrived as fresh 0 °C water and a uniform tracer did not stay uniform (gate: `tests/test_ocean_vanished_constancy.F90`). Three pieces: ONE definition in `src/shared_module_utilities/rdb_vanished_layer.inc` (`rdb_vl_is_live` / `rdb_vl_conc` / `rdb_vl_column_conc` — fillers read off their donor — / `rdb_vl_holds_live_conc` / `rdb_vl_merge_content`), `#include`d into each consuming module's `contains`; ONE enforcement point, `multilayer_state_t%enforce_vanished_content`, at the tail of `ocean_dyn_step_split` (a textual no-op on a column with no filler ⇒ sigma / z*-lite / `eulerian_z` are bit-identical; host twin `enforce_vanished_content_host` for the pre-`enter_data` seed); and two gates — `&vcoord_nml check_vanished_content` (fail-loud tripwire, default off, ON in every shipped namelist that vanishes layers) plus the diff-aware `vanished-layer` pre-commit hook. Substituting something OTHER than `c_live` on a vanished layer is legal and several consumers do it deliberately (EOS → reference T/S; sponge → nearest massive layer; diagnostics → NaN; melt far-field → skip) — declare it with `! vanished-ok: <reason>`. Contract + the full consumer table: [`src/core/ocean/README.md`](src/core/ocean/README.md) §"The vanished-layer content rule". The one-sided version of the guard (read side floored, write side not) is what deleted salt/heat un-budgeted on day 16 of the ISOMIP+ z_fixed case.
- **Tests with NetCDF/HDF5**: run with `OMP_NUM_THREADS=1` (via CTest `ENVIRONMENT`) — non-thread-safe HDF5 + `-stdpar=multicore` crashes otherwise.
- **Writing GPU tests (even with no GPU on the dev box)**: the GPU build is `-gpu=...,mem:separate` (**no** managed/unified memory), so a `do concurrent` / `!$acc` kernel gets **NO** implicit host↔device copies — *every* array a kernel touches must be device-present, or it silently reads/writes stale host memory (symptoms: budget/accumulator stays `0`, or division by an unmapped scratch array yields `NaN`, with no crash). A test that only ever runs on the multicore/host build will pass while being GPU-broken. Rules of thumb: **(1)** before calling any kernel, map every state object AND its scratch companion — e.g. mapping `ms` but forgetting `ct` (`call ct%enter_data()`) leaves `continuity`'s `h_face` buffers off-device → NaN. **(2)** Arrays mapped `create` by `enter_data` (fluxes the production step recomputes on-device, e.g. `mass_flux_*_layer`) do **not** carry host values you set before the map — push them with `!$acc update device(...)` after `enter_data`. The ONE exception is `scratch_3d_buffer_t`: its `enter_data` attaches with `copyin` of the HOST payload (zero from `init`, or what a warm-restart read restored), so a buffer whose producer is SKIPPED for a step reads as the zero `init` promised on both toolchains (gated by `rdb_test_scratch_3d_device`; without it MOM6's `diffu(u[n-1])` reuse in the `pred_corr` predictor read the allocator's leftovers at step 1) and a checkpointed buffer (`hvisc_du_visc`/`dv_visc`) keeps its restored value (an earlier device zero-fill wiped it — `restart_bit_exact_viscous`). **(3)** Host-set inputs read back after a kernel need `!$acc update self(...)`; local scratch passed straight to a `*_one_impl` needs its own `!$acc enter data copyin/create ... exit data`. **(4)** SETUP-TIME host editors must run BEFORE the map, not after: `metrics_apply_land_mask` is plain host loops with no `!$acc update` of its own (by design — the host edit is what the copyin captures), so masking a metrics slot that `make_cartesian_metrics` has already mapped leaves the DEVICE all-wet. The island then does not exist on the GPU while every host-side assertion, reading the masked host copy, still agrees with itself — pass the mask through `make_cartesian_metrics(metrics, grid, wet_mask=...)` and the order is unskippable. **(5)** An `ieee_get_flag` assertion only sees what THIS thread's FPU did: a kernel on the device raises no host flag, so such a gate is vacuous on the offload build and must be skipped there (`KERNEL_IN_HOST_FPENV` in `tests/test_continuity_multilayer.F90`), never left to "pass". All these directives are inert no-ops on host/multicore builds, so add them unconditionally. Canonical template: `test_open_boundary_out_closes` in `tests/test_ocean_conservation_salt_heat.F90` (its comment block documents the `mem:separate` contract). **Verify on the actual GPU build** (`source gadi_nvhpc.sh` on Gadi) — a green multicore run proves nothing about device data motion.
- **Never `!$acc update self` / `copyin` a WHOLE derived type with allocatable components**: the aggregate D→H copy overwrites the host component descriptors with DEVICE addresses, so the next host read of a component (`obj%arr(...)`) segfaults. Always update/map the COMPONENT arrays (`!$acc update self(obj%arr_x, obj%arr_y)`), never the aggregate `obj`. (Bit us in `test_ocean_ice_evp` — the host `check` read `stress%tau_x` right after `!$acc update self(stress)`.)
- **Persistent kernel workspaces**: never local-allocate scratch in per-step kernels on `-stdpar=gpu` (stalls device on entry/exit). Use module-level allocatables, lazy-allocate via `*_workspace_ensure`, release via `*_cleanup` in the exit-data path. Keep in-loop names short — but NOT with `associate` over derived-type components (see the ifx gotcha below); spell the component out or pass it as an explicit-shape `_impl` dummy.
- **Outer-shim + flat-impl call sites** (for the array-of-derived-types device indirection rule — see memory): the tracer-registry loops (`ocean_halo_exchange_ml_state`, the per-tracer remap/vdiff/hdiff drivers). Top-level allocatables reach the device directly; only the registry indirection needs the shim.
- **Sharing a helper BODY across kernel modules — use an include, not a hand-duplicated `*_impl` copy.** When several modules must agree on a small `pure` (+ `!$acc routine seq`) leaf helper that is called from inside a `do concurrent`, put the body in `src/shared_module_utilities/<concern>.inc` and `#include` it in each consumer's `contains`. Every module gets a module-LOCAL `private` copy the compiler inlines readily, and the tree has ONE source text to edit. See that directory's README for the admission test (rule-not-routine, small + leaf, called inside a kernel, one concern per file) and the naming rules. **Measured** (nvfortran 26.5, `-O2 -stdpar=gpu -gpu=cc70,mem:separate`, V100): a per-cell accessor written as a raw inline expression, as an included module-local helper, and as a **cross-module** `!$acc routine seq` helper produces device PTX with **zero `call` instructions in all three** — so on this compiler version the cross-module form inlines too, and the include's value is the single source of truth, not a speed-up. This REPLACES the older advice to hand-duplicate helper bodies as `*_impl` copies; existing duplicated copies are historical, not a pattern to extend. (Note: passing an array *element* to any helper makes nvfortran report the enclosing array as a whole-array `implicit copy` rather than a sliced `copyin` — free under `mem:separate` when the array is already device-resident, not free otherwise.)
- **Hoist neighbour-column builds out of per-layer loops**: any kernel shaped `for cell { for k { for nbr { build_column(nbr); query_at_z(k); }}}` should build the neighbour column once per cell, then query per-k. Done for Jacobian BPG kernels (2-2.4× per kernel).
- **Formula bathymetry setters must fill ghost rows**: `set_bathymetry_*` (spoon/seamount/…) that leave `h_layer=0` at ghost rows make EOS fall back to ρ=ρ₀ → a spurious density jump at wall-adjacent faces → ~12-hr e-fold blowup. The file-loader path fills ghosts; formula paths must too.
- **Formula-bathymetry length scales are in GRID units**: `set_bathymetry_{spoon,seamount}` work in grid coordinates — metres on Cartesian, **degrees** on spherical/curvilinear. The metres `&ocean_topo_nml slope_scale` knob is converted at the dispatch via `topo_length_to_grid_units`; without it a metres scale against a degrees position makes `exp(−r²/L²)≈1` everywhere and the basin collapses to a flat `peak_depth` (the spherical-seamount bug).
- **`if/else` clamps LAUNDER NaN under `-fast` (host AND device)**: nvfortran's
  relaxed-FP default (+`-fast`, no `-Kieee`) lowers `if (u > hi) u = hi; else if
  (u < lo) u = lo` to a NaN-blind min/max select — a NaN input comes OUT as
  `lo`/`hi` (reproducer: `tmp_local_artifacts/nan_clamp_repro.f90`; `-Kieee`
  preserves NaN but costs ~10-20%). Consequence: a clamp can silently convert
  corruption into plausible extreme values (bit us as ±maxvel phantom
  velocities — LOGBOOK Phase 10.17). Any clamp that can ever see non-finite
  data must be `ieee_is_finite`-guarded (see `apply_maxvel_clamp` /
  the truncation NaN-catch for the pattern); comparisons with NaN are FALSE, so
  `if (abs(u) > thresh)` guards also silently skip NaN — pair them with the
  catch, never rely on them to bound corrupted data.
- **EFP/fixed-point reductions must propagate non-finite input** — same hazard
  class as the NaN-blind clamp above, in the reduction layer instead of a
  clamp. The device-callable EFP decompose (`rdb_ocean_console_stats::
  efp_decompose_impl`) originally had NO NaN handling at all: a NaN summand
  hit `int(NaN, int64)` (compiler-undefined) and the fixed-point bins it fed
  reconstructed as a silent, plausible-looking finite number — observed as
  the console printing `En 0.000E+00` / `Salt 0.000` / `Temp 0.000` on an
  already-diverged (`Mass NaN`) run instead of aborting, which let twelve
  tier-1 stability-matrix cells slip past their `completed;finite` XFAIL pin
  once `reproducing_sums=.true.` became the default (a0755ada1). Fixed by
  giving `efp_t` (`rdb_efp`) a `poison` counter: `efp_decompose`/
  `efp_decompose_impl` set it (via `ieee_is_finite`, never a bare `>=`/`<`
  comparison) for a NaN, `+-Inf`, or bin-1-overflowing summand,
  `efp_plus`/`efp_minus` add it forward, `efp_to_real` returns
  `ieee_value(x, ieee_quiet_nan)` whenever it is nonzero, and
  `efp_to_transport`/`efp_from_transport` carry it through the SAME
  `MPI_SUM` collective as the bins (`EFP_TRANSPORT_WIDTH`) so one poisoned
  rank makes every rank read NaN identically. The plain-FP path
  (`reproducing_sums=.false.`) already propagated NaN correctly through
  ordinary IEEE `+`/`MPI_SUM` — nothing to fix there. Tests:
  `tests/test_efp.F90` (`test_efp_poison_propagates_nan`/`_inf`,
  `test_efp_poison_cross_rank_transport`, `test_efp_finite_path_bit_identical`),
  `tests/test_ocean_console_stats_efp.F90`
  (`test_efp_console_kernel_nan_propagates`), `tests/mpi/test_efp_mpi.F90`.
- **Assumed-shape dummies in `do concurrent` kernels**: NVHPC walks the descriptor with per-launch memcpys; use explicit-shape args (`arr(nx,ny,nz)`); cadence-bounded code may waive with `! assumed-shape-ok: <reason>`; enforced diff-aware by pre-commit (`dc-assumed-shape`).
- **A host-gated call that hands a state array to an EXTERNAL subroutine still costs, even when never taken**: passing e.g. `ms%mass_flux_x_layer` as an `intent(inout)` actual makes nvfortran treat the array as escaping and pessimises *every* `do concurrent` in the calling routine. Measured at **+4.8%** total solver time (600×600×50 ocean, all of it in `ocean_continuity`, 8.97 → 10.58 s) for a `call porous_narrow_3d(...)` sitting behind `if (metrics%use_porous)` with the knob OFF. Fix: write the guarded pass INLINE as a `do concurrent` in the same routine. A same-module helper does **not** help — it is the call, not the module boundary. (Not the extra dummy arguments either: removing them recovered nothing.) Suspect this whenever an "inert, host-gated" addition shows up in a profile.
- **Never wrap a `do concurrent` kernel in `associate` over a derived-type component**: under `-qopenmp` (how ifx maps `do concurrent` onto threads) ifx 2025.0/2026.0 evaluates an ASSOCIATE name whose selector is an allocatable component of a DT dummy as **zero** inside the loop body — silently, when the body also has a branch chain calling an inlinable `pure` module function. It killed the whole Coriolis term in `coriolis_adv_compute_tendencies_sadourny` (`f_corner` read as 0 ⇒ dead `dv/dt`, then NaN in the long runs) while gfortran and nvfortran were correct. NVHPC has the sibling failure with `!$omp target` inside `associate` over mapped workspaces. Spell the components out. Reproducer + writeup on the project wiki.
- **Never assign a scalar to a whole deferred-length `character` ARRAY**: `character(len=:), allocatable :: lines(:)` then `lines = ""` RE-ALLOCATES `lines` with `len=0` (F2018 10.2.1.3p3 — the length type parameter of the variable and of the expression differ). ifx implements this correctly; gfortran and nvfortran do not, which is how `split_to_lines` shipped for months silently discarding every in-memory namelist on ifx (so `rdb_ocean_create(namelist_text)` ran pure defaults with no error). Blank-fill through a SECTION — `lines(:) = ""` — which is not an allocatable variable and so pads in place everywhere.
- **Declare integer dims before the explicit-shape arrays that use them**: ifx warns #8586 ("non-standard out-of-order declaration") when an array's bound references a dummy integer declared on a later line. gfortran/nvfortran accept it silently. Fix: move `integer, intent(in) :: nx, ny, nz` (or equivalent) above the first `real(...) :: arr(nx, ny, nz)` in every procedure spec part; enforced by the `decl-order` pre-commit hook.
- **Never instantiate a default-initialised derived type in device code**: under `-acc=gpu` (nvfortran 23.9-26.5, any arch/`-O`), an `!$acc routine seq` procedure that builds a type WITH default initialisers (function result, `intent(out)` dummy, or a local) emits a device template `common._st_<hex(filename)>__<module>_NN`, and once `2*len(file stem) + len(module name)` reaches ~59 (locals: somewhat longer) `fort2` segfaults. Symptom: `fort2 TERMINATED by signal 11`, no diagnostic, and it VANISHES when you copy the file to a short name, so bisection lies — always bisect under the ORIGINAL filename. Pass the scalars instead (`rdb_ocean_cavity_melt`'s `*_f` variants thread the per-column `f_cor` that way). `intent(in)` bundles are fine.

## Coding Conventions

- File prefix: `rdb_` for modules, `test_` for tests (`tests/test_<name>.F90`, MPI tests in `tests/mpi/`). Extensions: `.F90`.
- Types: `_t` suffix. Constants: UPPERCASE. Everything else: snake_case.
- `use` always with `only:`. `implicit none`, `private` by default, explicit `public`. FORD docs: `!!` after declarations.
- Intrinsic modules need the modifier: `use, intrinsic :: iso_fortran_env, only: ...` (fortitude rule C122). Same for `iso_c_binding`.
- No `print *` — use `pic_logger` (`global_logger`).
- **Default new procedures to `pure`** — make a new function/subroutine `pure` unless it genuinely needs a side effect (I/O, logging, MPI, mutating module state, lazy `allocate`). Pure is the norm here; non-pure is the justified exception. See `FORTRAN_STYLE.md` §Pure by Default.
- Working precision: `wp` from `rdb_constants` (currently `real64`).
- MPI portability: `#ifdef USE_LEGACY_MPI` selects `use mpi` vs `use mpi_f08`.
- Dependencies: `pic` (types, strings, timers, logger), `pic-mpi`, `test-drive`, NetCDF-Fortran.

## Performance (V100)

| Config | Grid | GPUs | Wall time |
|---|---|---|---|
| Tasman 2 km (eddy-resolving, KPP + ALE + sponges) | regional | 1 | **~34 s / simulated day** |
| Double-gyre MOM6 ref (`double_gyre_mom6.nml`) | 44 × 40 × 2, dt=1200 | 1 | ~16 s / 30 days |
