# Physics-coverage regression suite — ocean (run-clean gate + golden compare)

A stdlib-only Python harness that runs a curated set of **bathymetry-free ocean**
namelist cases for a handful of timesteps each and asserts every one **runs
clean** (P0 — exit 0, no non-finite / crash markers) **and does not drift**
against a committed golden final-state summary (P1 — high-tolerance, not
bitwise).

- **P0 (`run_regression.py`):** runner skeleton + curated manifest +
  run-clean/NaN gate, CPU + GPU. No goldens.
- **P1 (`compare.py`, here):** golden field-summary generation +
  high-tolerance compare. Layered ON the P0 runner (reuses its case-execution
  machinery). Distinct mode from the pure run-clean gate.
- **P2 (`coverage.py`):** gcov corpus-coverage mode + per-closure report.
- **Entry point (`run_all.py`):** one command that invokes everything —
  compare (run-clean + drift) on CPU and/or GPU, plus coverage — returning a
  CI-ready exit code. Standalone (not a make/CTest target); CI can call it later.
- **Scope:** ocean only — the coastal regime was carved out into its own
  repository and is not a later phase here.

See `docs/regression_suite_plan.md` for the full design.

## What the suite runs

**35 bathymetry-free ocean cases**, each run for a handful of outer timesteps
(the runner rewrites `&time_nml` so `t_end = n_steps * dt_fixed`). Every case is
formula bathymetry (no DEM/file load) and finishes in O(seconds). The single
source of truth is `manifest.py` (each case carries the `tags` = closures it
exercises); this table is the human summary.

### Core dyn-core & analytical
| case | what it exercises |
|------|-------------------|
| `double_gyre_mom6` | canonical MOM6 double-gyre (NK=2): split-RK2, gprime reduced-gravity PGF, Sadourny PV Coriolis, linear drag, z* |
| `double_gyre_linear_nk10` | multi-layer (NK=10): FV-lite PGF, Sadourny-energy Coriolis, Smagorinsky viscosity |
| `geostrophic_adjustment` | barotropic (NK=1) free-surface + Coriolis geostrophic balance |
| `seamount` | quiescent stratified seamount: nonlinear EOS, rest-state PGF error, ALE remap |
| `benchmark_ale` | PPM ALE remap across many layers (NK=75) |

### Baroclinic / mesoscale
| case | what it exercises |
|------|-------------------|
| `eddy_test` | mesoscale eddy spin-down (NK=15) |
| `baroclinic_2layer` | baroclinic-instability re-entrant channel (NK=2) |
| `eady` | Eady linear baroclinic growth / thermal-wind balance |

### Vertical mixing, tracers & tides
| case | what it exercises |
|------|-------------------|
| `epbl_basin` | EPBL energetics boundary-layer mixing (Reichl-Hallberg) |
| `kpp_basin` | KPP boundary-layer mixing (the other vmix path) |
| `ideal_age_demo` | passive-tracer registry + transport (ideal age) |
| `body_tide_basin` | equilibrium body-tide forcing + scalar SAL |

### Boundaries, land masking & sea ice
| case | what it exercises |
|------|-------------------|
| `island_at_rest` | interior land masking (free-slip walls), quiescent — should stay at rest |
| `flow_past_island` | land masking under a driven flow (island wake) |
| `coriolis_coast` | Coriolis against a solid wall boundary |
| `sponge_real_demo` | SPONGE BC nudging / relaxation |
| `polar_freezeup_dynamics` | sea-ice EVP dynamics + ice-ocean coupling |

### P4 gap-fillers (seamount formula topo — turn on closures the base corpus missed)
| case | what it exercises |
|------|-------------------|
| `seamount_gm_redi_meke` | mesoscale lateral stack: isopycnal slopes + Gent-McWilliams + Redi + MEKE |
| `seamount_meke_backscatter` | MEKE backscatter (negative viscosity) + length-scale weights + BBL drag (+ Smagorinsky backstop) |
| `seamount_tidal_mixing` | St-Laurent/Simmons internal-tide bottom-intensified diapycnal mixing |
| `seamount_pgf_reconstruct` | FV-MOM6 PGF with in-layer **PLM** T/S density-integral reconstruction |
| `seamount_pgf_ppm` | FV-MOM6 PGF with in-layer **PPM** reconstruction (edges + Boole quadrature) |
| `seamount_conservative_floor` | Lagrangian grounding + conservative minimum-thickness borrow |
| `seamount_obc_baroclinic` | baroclinic open boundary: Orlanski radiation + Marchesiello nudging + tracer reservoirs (M2-tidal west / Flather east) |

### WENO tracer reconstruction (windowed drain)
| case | what it exercises |
|------|-------------------|
| `eady_weno5` | Eady front + windowed tracer-advect drain (`dt_tracer_advect_ratio=2`) with **WENO5-Z** face reconstruction (`tracer_recon="weno5"`) — the weno5 + PLM swept-average face helpers + rung-degradation ladder in `rdb_recon_weno.F90` |
| `eady_weno7` | same base + drain, **WENO7-Z** (`tracer_recon="weno7"`, `nghost=4`) — reaches the higher-order cubic-candidate stencil (`weno7_face_swept`) weno5 can't |

### Outer split scheme (`&ocean_bt_nml split_scheme`)
The base corpus runs the default `split_scheme="pred_corr"` (the MOM6
predictor-corrector; it became the default on 2026-09-14), so the
`SPLIT_SCHEME_PRED_CORR` branches in `rdb_ocean_dyn.F90` — predictor stage,
corrector on the time-mean `u_av/h_av`, `restore_state`, the `is_pc`
conditionals — are the ones the whole corpus exercises. The `ssp_rk2`
branches are covered by `eady_weno5` and `eady_weno7`, which pin `"ssp_rk2"`
in their own namelists because the `pred_corr` v1 envelope REFUSES their
windowed tracer-advect configuration, so both branches of the dispatcher stay
golden-covered.

**Skipped on update — how a golden goes stale.** `compare.py --update-golden`
writes a golden only for a case that PASSES the run-clean gate; a gate-failing
case keeps its old file. Five goldens (`double_gyre_mom6`,
`double_gyre_dataovr`, `flow_past_island`, `coriolis_coast`, `island_at_rest`)
sat stale that way: the NaN missing-data sentinel in the T/S diagnostics made
their console `mean` print `NaN`, the gate rejected them, and every
regeneration skipped them — so three still held `ssp_rk2`-era answers under an
unpinned namelist after the 2026-09-14 default flip, and `coriolis_coast` held
its pre-wind-fix (En = 0) answer. The console reduction now skips non-finite
cells (`missing=<n>/<total>` on the line), all five were regenerated, and an
update that skips any case now names the stale goldens and exits non-zero.

The two cases below PIN `pred_corr` explicitly. They predate the flip — they
existed to keep the then-non-default branches warm — and they stay because
the pin is what holds their goldens still across any future default move,
and because each is a distinct configuration in its own right (a
non-Lagrangian ALE vcoord under the pc loop, #395).

The STABILITY suite carries the broader axis: it re-runs every case whose
namelist does not pin a scheme under a forced `ssp_rk2` as `<case>__ssp_rk2`
(43 twins, ~2.5 min of extra tier-2 wall on 4 CPU workers), which is what
keeps the experimental scheme from rotting and is where its one open defect
(resting-state internal-gravity-wave growth) is recorded as a scoped XFAIL.
See `stability_manifest.py`, "The OUTER SPLIT-SCHEME axis".

| case | what it exercises |
|------|-------------------|
| `seamount_pred_corr` | quiescent stratified-free seamount (zstar_sigma, NK=15), `pred_corr` PINNED — rest-preservation guard: the pc loop holds the rest state to machine precision on a non-Lagrangian ALE coord |
| `double_gyre_pred_corr` | active wind-driven double-gyre (fv_lite PGF, NK=10, zstar), `pred_corr` PINNED — the corrector runs on **real** tendencies, a meaningful golden drift target |

**Coverage** (measured when the corpus was 28 cases; not re-measured at 35): **~53% of ocean-closure source lines** (gcov).
The two WENO cases also lift the shared tracer module
`src/tracer/rdb_recon_weno.F90` (weno5/7/9 + WENO-Z
swept-average faces, reused by the ocean windowed drain) from **~3% to ~28%**
(the residual is the untested `weno9` rung). Known gaps (future cases):
tripolar fold/bipolar (needs a tripolar grid), I/O / restart / budget
diagnostics (low physics value). Run
`coverage.py` for the live per-closure gap list.

## Files

| file | purpose |
|------|---------|
| `manifest.py` | the curated case list — nml path, `n_steps`, `timeout_s`, `backends`, physics `tags`, optional per-case `tol`, optional per-field `atol` + `atol_reason` |
| `run_regression.py` | the P0 orchestrator (CPU-serial \| GPU work-queue), the run-clean/NaN gate, the report |
| `compare.py` | the P1 golden-summary compare (`--update-golden` \| `--compare`) |
| `golden/<case>.json` | committed per-case final-state field summaries (tiny) |
| `coverage.py` | P2 gcov coverage mode — measures ocean-closure coverage of the corpus + emits the gap list |
| `run_all.py` | **single entry point** for the golden suite — invokes compare (CPU/GPU) + coverage, aggregates to one pass/fail + CI-ready exit code |
| `stability.py` | **the two-tier stability suite** (see below) — parses the model's own console time series and asserts on PHYSICS, not on a golden |
| `stability_manifest.py` | its case list: every tracked ocean namelist (73 base cases), with per-case run length, physics assertions, tier-2 downscale spec and known-failure markers |
| `downscale.py` | the dimensionless-number rules a tier-2 twin must satisfy, plus the standalone checker that validates every twin |
| `README.md` | this file |

## Build the CPU app

The suite runs the normal `rdb` executable. Load a gfortran + NetCDF toolchain
however your site does it (`module load`, Spack — see `environments/` —, conda),
then build just the app (not the whole test suite):

```bash
# load your gfortran / NetCDF toolchain first — and only that one:
# stacking a second toolchain in the same shell puts two NetCDF builds on the
# link line and the NetCDF tests fail in confusing ways.
cmake -B build_gcc -S . -DRDB_ENABLE_GPU=OFF -DRDB_ENABLE_MPI=OFF
cmake --build build_gcc --target rdb -j
```

This produces `build_gcc/rdb`, which the runner auto-discovers.

## Run

```bash
# CPU, serial — the P0 bar
python3 tests/regression/run_regression.py --backend cpu \
        --build-dir build_gcc --out tests/regression/last_cpu.json

# GPU — one case per GPU at a time, farmed across visible devices
python3 tests/regression/run_regression.py --backend gpu \
        --build-dir build_gpu --gpus 0,1,2,3 --out last_gpu.json
```

`--help` documents every flag (`--backend`, `--build-dir`, `--binary`, `--out`,
`--gpus`, `--cases`, `--budget-min`, `--keep-scratch`, `--scratch-root`).

The runner exits **nonzero** if any case fails the gate or the budget is blown,
so it drops straight into CI.

## What each case run does

For every manifest case the runner:

1. Copies the committed nml into a per-case scratch dir (default
   `tmp_local_artifacts/regression/<case>/`, which is git-ignored).
2. Patches `&time_nml` so the run is `n_steps` **outer** steps — it reads
   `dt_fixed` and sets `t_end = n_steps * dt_fixed` seconds, normalizing
   `time_unit` to `"second"`. The committed nml is never mutated.
3. Redirects `&output_nml output_dir` to the scratch dir **and** runs the
   process with `cwd=<scratch>`, so no run ever writes into the repo tree.
4. Runs `rdb <temp.nml>` with the per-case `timeout_s`, capturing exit code
   + combined stdout/stderr.

**Gate (P0):** PASS iff exit code `0` **and** no bad markers appear in the
output — `nan`, `inf`/`infinity` (word-boundary matched so benign banner text
does not false-trip), `error stop`, `not finite`, `abort`, `panic`,
`segmentation`, `floating point exception`, `backtrace`. Scratch is removed on
success and kept on failure for debugging (`--keep-scratch` keeps it always).

## Backends

- **`--backend cpu`** runs the cases serially (the gfortran CPU build).
- **`--backend gpu`** farms cases across the visible GPUs with a work queue —
  **exactly one case per GPU at a time** (a GPU is never shared, same rule as
  `ctest` on this repo), pinning each case with `CUDA_VISIBLE_DEVICES=<id>`.
  `--gpus 0,1,2,3` sets the device list explicitly; otherwise it auto-detects
  via `nvidia-smi -L` (falling back to a single device `0` if absent).

## Budget

Each backend must finish in **under 30 minutes** (`--budget-min`, default 30).
The runner tracks cumulative wallclock and **fails loud** the moment a backend
exceeds it — a silent 45-minute suite is a bug. Per-case `timeout_s` catches an
individual hang. In practice the CPU pass is ~30 s total (the two multi-layer
cases, `eddy_test` ~14 s and `benchmark_ale` ~10 s, dominate).

## Golden-summary compare (P1)

`compare.py` adds high-tolerance **drift** detection on top of the P0 gate. It
reuses the P0 case-execution machinery (temp-nml patch + output isolation +
run-clean/NaN gate) and then compares a tiny per-case **final-state summary** to
a committed golden.

**What a golden is** — for each case, at the final step: the `min/max/mean` of
every prognostic the run prints as a `[diag]` line (SSH/η, u, v, T, S, KE, any
active tracer, ice_*), plus the scalar diagnostics on the `[stats]` line (`En`,
total `Mass`, and `Salt`/`Temp` when thermodynamics are reported). These are
**parsed from the model's own console output** (the last occurrence of each
field) — no new instrumentation. Stored as `golden/<case>.json` (a few hundred
bytes each; committed).

**Guaranteeing a final-state print** — the committed nmls emit diags/stats on an
hours-or-days cadence, so a short regression run would otherwise only print at
`t=0`. The runner's `patch_namelist` (`force_emit=True`, default) rewrites
`&logging_nml status_interval` and `&ocean_diag_nml dt_out` down to half an
outer step (seconds, after the `time_unit` normalization) so **both fire every
step** — the final step always emits a usable summary. (Harmless to the P0 gate:
it only adds log lines / more NaN sampling.)

```bash
# (Re)generate goldens from a trusted CPU build, then commit them.
# (gfortran / NetCDF toolchain loaded in this shell)
python3 tests/regression/compare.py --update-golden --backend cpu --build-dir build_gcc

# Compare (default mode) — CPU:
python3 tests/regression/compare.py --compare --backend cpu --build-dir build_gcc
# Compare — GPU, farmed across devices (NVHPC toolchain, in a FRESH shell —
# never stack it on top of the gfortran one):
python3 tests/regression/compare.py --compare --backend gpu --build-dir build_gpu --gpus 0,1,2,3
```

**Tolerance model** — each summary value is compared numpy-`isclose` style:
`|golden − run| ≤ atol_field + rtol·scale`, where `scale` is the *field's*
magnitude (diag field: `max|min|,|max|` over golden **and** run; stats scalar:
`max|golden|,|run|`). An active field is bounded by `rtol·scale`, and a
spurious spin-up lifts the scale and **fails loudly**.

**There is no global absolute floor.** `atol_field` is `0` unless the case's
manifest entry names that field in an explicit `atol` map (`{"SSH": 1e-9,
"stats:En": 1e-18, ...}`) with an `atol_reason`. Until 2026-09-23 a global
`atol = 1e-7` applied to every value of every case, and it passed drifts it
should have reported: on the v0.1.0 integration branch `seamount` and
`seamount_pred_corr` printed **PASS** next to a 0.70 / 0.98 relative drift,
because their whole fields sit at 1e-11 m / 1e-13 m/s and any change below
1e-7 — a 10⁴-fold growth of that noise included — was inside the floor; the
same floor hid `ideal_age_demo` SSH (8e-7 m) moving 9 %, and a
`seamount_pred_corr` golden that `main` itself no longer reproduces. A
**physically-zero field carrying roundoff noise** (a resting case's `u` /
`SSH` / `KE` / `En`) therefore declares its floor in `manifest.py`
(`RESTING_NOISE_ATOL`: ≥100× the measured noise, so a noise change passes
but a spin-up past 1 nm / 1 nm/s fails), and the report names every value
that passed *only* on it (`ok; N value(s) inside manifest atol only`, plus a
`passed on its manifest atol only:` list). The report's `|d|/allow` column is
the verdict quantity (≤ 1 on every PASS row); `rel drift` is `|d|/scale`.
`compare.py` fails a case whose `atol` names a field its golden does not
carry, and asserts the verdict can never be PASS while a checked value is
outside its allowance.

```bash
python3 tests/regression/compare.py --self-test   # no model; ALL PASS expected
```

The self-test replays the measured seamount pair (it must FAIL without the
manifest floor and PASS with it, and a 1000× spin-up must still FAIL), the
`ideal_age_demo` drift, NaN / absent-field / exact-zero cases, the verdict
invariant, the `atol` validation, and checks every manifest `atol` against its
committed golden. It is registered in ctest as `rdb_regression_compare_selftest`.

**Chosen defaults (measured, not guessed)** — `rtol = 1e-3`. The
worst REAL-field (non-roundoff) CPU(gfortran)↔GPU(nvfortran) relative drift
across the 22 well-behaved cases is **1.4e-6** (`ideal_age_demo` SSH:min);
`rtol=1e-3` is ~700× that — loose enough to absorb FMA / reduction-order /
transcendental-library divergence (GPU is deterministic run-to-run here) yet
tight enough that a real physics regression fails. Goldens are generated on the
**CPU** build; `--compare --backend gpu` then passes against the same goldens
within these tolerances (measured with the old global floor; the GPU leg has not
been re-measured since the floor became per-field — a GPU FAIL on a field of a
resting case is the cue to add it to that case's `atol`, not to widen `rtol`).

**Per-case `tol`** — two early-transient baroclinic-instability cases
(`baroclinic_2layer`, `eady`) carry a manifest `tol=0.25`: their tiny
cross-channel `v` / free-surface `SSH` **extrema** diverge 9–15% CPU↔GPU at the
noise floor, while their integrated `En/Mass/Salt/Temp` agree within the global
`1e-3`. Per the design, such cases are **flagged with a per-case tol (documented
inline in the manifest), not hidden by widening the global bound**.

Exit is nonzero if any case fails the gate, drifts past tolerance, is missing a
golden, or the budget is blown — CI-ready.

## Coverage (P2)

`coverage.py` measures, with **gcov**, what fraction of the **ocean
physics/closure source** the corpus exercises, and surfaces the biggest
UNCOVERED closure files as the gap list that drives adding cases (P4). It is
ocean-only and stdlib-only Python.

```bash
# gfortran + gcov + NetCDF toolchain loaded in this shell
python3 tests/regression/coverage.py    # build_cov → run corpus → report
```

End to end it:

1. **Builds a coverage binary** in `build_cov/` with the existing
   `-DRDB_ENABLE_COVERAGE=ON` option (gfortran `-O0 -g --coverage`, GPU + MPI
   off). Reuses an existing `build_cov/rdb` unless `--rebuild`. **Do not add
   `-DCMAKE_BUILD_TYPE=Debug`** — the coverage flags already force `-O0 -g`, and
   Debug's `-fcheck=bounds` trips a latent out-of-bounds in `profiler_report`
   (`rdb_profiler.F90`) at exit, crashing every otherwise-clean case.
2. **Runs the ocean manifest corpus SERIALLY** against it (reusing
   `run_regression.run_case`'s temp-nml patch + output isolation) so the `.gcda`
   counters accumulate. Counters are zeroed once up front.
3. **Captures coverage** — always a stdlib gcov fallback (exact per-line counts
   parsed from the `.gcov` files), plus the lcov/genhtml HTML report + a per-file
   cross-check *if* lcov is on PATH (not required — the env is flaky).
4. **Reports** overall line coverage (whole tree + ocean-closure subset), a
   per-closure table sorted ascending, and an explicit **GAPS** list (closures
   at/under `--gap-threshold`, default 10%). Writes `build_cov/coverage_summary.json`.

`--help` documents every flag (`--build-dir`, `--rebuild`, `--jobs`, `--gcov`,
`--cases`, `--no-lcov`, `--gap-threshold`, `--budget-min`, `--skip-run`, `--out`,
`--scratch-root`). Exit is nonzero only on a hard failure (budget blown / no
coverage data) — **low coverage is the report's point, not a failure**.

## Adding a case

Append an entry to `CASES` in `manifest.py` — path to a **bathymetry-free** ocean
nml, a small `n_steps`, a `timeout_s`, the valid `backends`, and `tags` naming
the physics it targets. Keep it fast (O(seconds)); drop or shorten anything that
can't be. Then generate its golden with
`compare.py --update-golden --backend cpu` and commit `golden/<case>.json`. If
`--compare --backend gpu` shows a genuinely large-but-benign CPU↔GPU spread for
that case (confirm the integrated `En/Mass` still agree), add a documented
per-case `tol` rather than loosening the global tolerance.

---

# The two-tier stability suite (`stability.py`)

Everything above is the **golden-drift** suite: it runs each case for 6–10 outer
timesteps and compares the final state to a committed answer. That is a good
gate for what it is, and it has two structural blind spots that let real
defects ship:

1. **6–10 steps is 1.7 simulated hours** at `dt = 600`. The defects found on
   2026-09-11 manifest at hour 3, hour 7, step 131 and "days". The suite stops
   before the physics exists.
2. **A golden enshrines whatever the run did.** A golden captured from a broken
   run passes forever. `eady` passed its golden for months *because* the golden was the
   over-damped answer.

`stability.py` is the answer to both. It runs each namelist long enough that
its physics manifests and asserts on the model's **own console time series** —
`[stats]` (En, MaxCFL, Mass, Salt, Temp) and the per-step closed-budget
residuals — against statements about *physics*. There is no stored answer, so
there is nothing to enshrine.

```bash
# Tier 2 — the CI gate. gfortran CPU, downscaled twins, ~2 min on 6 workers.
# (gfortran / NetCDF toolchain loaded in this shell)
python3 tests/regression/stability.py --tier 2 --build-dir build_gcc --jobs 6

# Tier 1 — full scale, GPU, local / nightly. ~10 min on 4 GPUs.
# (NVHPC toolchain, in a FRESH shell — one toolchain per shell)
python3 tests/regression/stability.py --tier 1 --backend gpu \
        --build-dir build --gpus 0,1,2,3

# One case, every assertion printed (not just the failures):
python3 tests/regression/stability.py --tier 2 --cases eady -v

# Who tests the test: replay the MEASURED 2026-09-11 failure signatures
# through the assertions and check each is classified correctly. Runs no
# model; a second or two. Also in CI.
python3 tests/regression/stability.py --self-test
```

The self-test is worth reading before trusting the suite. It pins, among other
things, that the 2026-09-11 instability trips `energy:no-fast-growth` **and
that an absolute energy ceiling and a whole-run growth fit both miss it
entirely** — which is the whole reason the windowed rate exists.

## The two tiers

| | **Tier 1** | **Tier 2** |
|---|---|---|
| purpose | full-scale physics gate | **the CI gate** |
| where | local box / self-hosted GPU runner / nightly | standard GitHub-hosted runner |
| toolchain | NVHPC `nvfortran`, `RDB_ENABLE_GPU=ON` | **gfortran**, `RDB_ENABLE_GPU=OFF` |
| hardware | 1–4 GPUs | 2–4 CPU cores, ~7 GB RAM, no NVIDIA hardware |
| cases | **all 62** tracked ocean namelists | **49** (13 are tier-1 only — see below) |
| run length | per case, sized to when the physics appears (`eady` 60 d, `seamount` 30 d, `baroclinic_2layer` 40 d) | 40–140 simulated hours |
| wallclock | **21 min across 4 GPUs (1270 s)** — `baroclinic_15layer` alone is 1030 s of it, because its instability does not appear until day 30 | **155 s on 4 workers; ~575 s serial** |
| drops | — | the assertions that genuinely need a long integration (a baroclinic growth rate cannot be measured in 40 simulated hours). They are reported as `SKIP … TIER 1 ONLY`, never weakened to pass. |

Tier 1 is **not** wired into GitHub Actions: a hosted runner has no GPU. Run it
on a GPU box before merging anything that touches ocean physics, and from a
nightly job. Tier 2 runs on every push and pull request —
`.github/workflows/ocean-stability.yml`.

Budget note: tier 2 is **155 s wall on 4 workers** on this dev box (575 s of
CPU in total, longest single case 28 s). A hosted runner is roughly 2× slower
with 2–4 cores, so expect **5–10 minutes** for the sweep, plus the gfortran
build (cached) and `ctest -R rdb` (186 tests, 65 s measured on gfortran).
Comfortably inside a 20-minute job.

Current state of both tiers on `fix/eos-every-dynamics-step`:

```
TIER 1   62 cases, 1269 s (4x V100)   54 PASS   0 FAIL   8 XFAIL   0 XPASS
TIER 2   49 cases,  155 s (4 workers) 44 PASS   0 FAIL   5 XFAIL   0 XPASS
```

## What it asserts

Every assertion prints **what was expected, what was observed, and what it
means** — a contributor should be able to act on a red run without reading this
file.

| assertion | what it says |
|---|---|
| `completed` | reached the requested step count, exit 0, no abort marker. Uses the model's own `Total steps:` line, not the `[stats]` cadence. |
| `finite` | no NaN/Inf in any conserved scalar, or in any field **extremum**. (Field `mean` is excluded on purpose — see *the NaN that is not a bug*.) |
| `conserve:Mass` / `:Salt` / `:Heat` | the model's own closed-budget residual stays at roundoff (1e-11 closed, 1e-9 with an open boundary or surface flux; a healthy run sits at 1e-14–1e-15). |
| `energy:rest` | a **quiescent** case must stay at rest — peak `En` under a velocity bar. |
| `energy:rest-settles` | spurious motion must **equilibrate**: the final `En` must be below 95 % of the run's peak. Sharper than the magnitude bar — a spurious pressure gradient that settles is a tolerable discretisation error; one still at its maximum when the clock runs out is not. |
| `energy:no-growth` | **an unforced/adiabatic run's kinetic energy MUST NOT GROW.** This is the gate that catches what NaN and budget checks cannot: on 2026-09-11 an instability grew `En` 30–60× while Mass, Salt and Temp stayed exact to every printed digit and nothing crashed. |
| `energy:no-fast-growth` | an unforced case may grow at its **physical** rate and no faster — the fastest rate over any *window* of the run, not a whole-run fit. (A whole-run fit cannot tell them apart: `eady`'s instability went 32× in two days then sat flat for 23, which averages to a perfectly physical-looking `sigma`.) |
| `energy:destabilises` | an instability case must actually **grow somewhere**. A baroclinic case that only ever decays has never demonstrated its own physics. |
| `energy:growth-rate` | tier 1 only — the measured growth rate must match the **analytical** expectation the case exists to demonstrate. Fitted on `En` by default, or — where `sigma_expected` names a `field` — on that `[diag]` field's extrema from `fit_from_day` on (`eady` reads `max|v|`, which is zero in the basic state, because `En` there is 4.3e-3 m²/s² of jet that drifts *down* while the mode climbs three decades underneath). |
| `energy:responds` | a **forced** case must actually move. A forcing path that is silently a no-op passes every NaN, budget and golden check ever written — doing nothing is exactly conservative. This is the only thing that sees it, and it caught one (`coriolis_coast`). |
| `energy:bounded` / `energy:saturating` | a forced case stays inside a physical envelope, and is saturating by the end of a spin-up-length run. |
| `cfl:bounded` / `cfl:no-runaway` | MaxCFL under its ceiling, and a monotone climb already past 10 % of the ceiling must not project to a trip within one more run length. |
| `claim:<name>` | the namelist header's **own** testable statement. Never weakened to pass and never satisfied by editing the namelist: when a header and the code disagree, that *is* the finding. |

### Statuses

`PASS` · `FAIL` (a regression — fix it) · `XFAIL` (a **known, documented**
defect: the reason and a pointer print with it, so red is informative rather
than flaky) · `XPASS` (a case marked known-failing that now passes — drop its
`known_failure` marker).

### The NaN that is not a bug

`fill_tracer_impl` writes an IEEE quiet NaN into any cell whose layer has
vanished, deliberately, so a pinched-out or below-bottom cell reads as *missing*
rather than as a plausible 0 °C / 0 PSU. Land, halo and vanished bed cells are
all such cells, so **every tracer-concentration diagnostic carries NaN** and the
whole-array `mean` on the `[diag]` line is NaN for any case with land, a halo or
a vanishing coordinate. The printed `min`/`max` survive because `minval`/`maxval`
(and the GPU `reduce(min:)`) are NaN-blind — they launder it.

`run_regression.py`'s bare-text `\bnan\b` scan cannot make that distinction and
fails the flagship `double_gyre_mom6` on sight. This suite parses the numbers
and decides per field: extrema **are** gated, the sentinel mean is reported as
an `ARTIFACT` line. That the mean is unusable is itself a defect worth fixing
(the reduction should skip the sentinel) — it is surfaced, not hidden.

# The vertical-coordinate REST MATRIX (`vcoord_matrix.py`)

Everything above runs each shipped namelist under **whatever vertical
coordinate that namelist chose**. Nothing ran ONE problem under ALL of them, so
nothing could answer the question a user actually has:

> **Which vertical coordinate is trustworthy on which geometry?**

That gap is not theoretical. In two days on the cavity branch it hid, at the
same time: `VCOORD_ZSIGMA` collapsing its whole column into the bed layer (with
`Σ target_h = H + η` still exact, which is why no conservation test caught it);
`VCOORD_ZSTAR_SIGMA` being numerically indistinguishable from `VCOORD_SIGMA` in
all 21 shipped namelists that select it; `VCOORD_ZSTAR_FULL` resolving the
*wrong half* of the column under a lid; a rest-state growth mode over **any**
sloping boundary that no shipped case saw because every sloping case carries
viscosity or runs short; and sigma failing outright at `rx0 = 0.73`.

Every one of those is invisible to a per-namelist suite and obvious in a
family × geometry table. **The table is the product.**

## How it is built

One canonical template — [`vcoord_templates/rest_matrix.nml.in`](vcoord_templates/rest_matrix.nml.in)
— plus one parameter dict per **problem** and one per **family**, in
[`vcoord_matrix.py`](vcoord_matrix.py). `build_matrix()` emits one namelist per
cell into `tmp_local_artifacts/vcoord_matrix/` at manifest-import time and
registers the rows through the ordinary `_case(...)`. Dozens of
near-identical namelists are deliberately **not** checked in: they rot, a
reader cannot tell which cell of the matrix a given file is, and a change to
the shared problem then has to be applied by hand N times. The emitted files
are plain text, carry a generated header naming their cell, and are left on
disk — any cell reproduces by hand with
`./rdb tmp_local_artifacts/vcoord_matrix/<cell>.nml`.

```bash
# Just the matrix, tier 2 (CPU):
python3 tests/regression/stability.py --tier 2 --build-dir build_gcc \
        --jobs 6 --tags vcoord_matrix

# The full matrix, tier 1 (30 simulated days per cell; 10 for the seamount
# and rx0 problems), on EACH toolchain --
# one leg at a time with --tags vcoord_matrix_inviscid / _viscous if wanted:
python3 tests/regression/stability.py --tier 1 --backend gpu \
        --build-dir build_cc70 --gpus 0 --tags vcoord_matrix \
        --scratch-root tmp_local_artifacts/st_gpu --out tmp_local_artifacts/vcm_t1_gpu.json
python3 tests/regression/stability.py --tier 1 --backend cpu \
        --build-dir build_gcc --jobs 4 --tags vcoord_matrix \
        --scratch-root tmp_local_artifacts/st_gcc --out tmp_local_artifacts/vcm_t1_gcc.json

# ...then pin the markers and envelopes from them (never by hand):
python3 tests/regression/vcoord_matrix_pin.py \
        --t1 gfortran=tmp_local_artifacts/vcm_t1_gcc.json \
        --t1 nvfortran=tmp_local_artifacts/vcm_t1_gpu.json \
        --t2 gfortran=… --t2 nvfortran=… --provenance "…"
```

A run prints the family × problem table at the end and writes it next to the
ordinary report as `<out>_vcoord_matrix.json`, one record per cell — so drift
in the table can be diffed against last week's instead of being read out of a
terminal.

## The problem

**One problem, at rest.** A motionless, stably stratified ocean on an f-plane
(48 × 6 × 15 at 2 km, `dt = 600 s`, `f = −1.409e-4`), with **no forcing and no
vertical mixing**, and with the isopycnals laid **flat in geopotential `z`** by
`&ocean_zinit_nml source="linear"` — *not* by a layer-index profile, which
would tilt them with a terrain-following coordinate and hand the run real
available potential energy to convert. Over a slope that distinction is the
whole experiment.

Such a state is the global minimum of potential energy under rigid boundaries,
so APE ≡ 0 and there is **no energy source of any kind**. Every joule of
kinetic energy a cell develops was manufactured by the discretisation.

### The fixed configuration — what v0.1.0 recommends

Every cell of both legs runs the configuration the 2026-09-21 rerun found
necessary (`design/vcoord_ale_audit.md` in the prototypes repo), and **this is
the configuration v0.1.0 recommends for any run over sloping topography**:

| knob | why |
|---|---|
| `&ocean_pgf_nml form = "fv_mom6", reconstruct_for_pressure = .true.` | the exact finite-volume PGF: the sigma rest-state seed falls from 2.6e-8 to 9.6e-16 m s⁻² and the slope plateaus fall 12–17 decades to machine zero |
| `&vcoord_nml remap_boundary_extrap = .true.` | linear-exact boundary cells in the remap (the first-order boundary cell was the rest-state amplifier: growth 0.355 → 0.047 /day) |
| `&vcoord_nml remap_nonuniform_weights = .true.` | PLM/PPM linear-exact on a stretched column (1.5e-01 → 3e-14) |
| `&vcoord_nml remap_check_preconditions = .true.` | fail loud the first step the remap is handed a negative layer or a mismatched column total, instead of silently creating or deleting tracer |
| `&vcoord_nml zfixed_closed_faces = .true.` (`z_fixed` only) | closed staircase faces with an open-column barotropic mode — without it the staircase PGF residual is unbounded |

`stability.py --self-test` asserts that every emitted cell carries all of
them.

The FINDING B fix, `&ocean_continuity_nml renorm_consistent_flux = .true.`,
and MOM6's fast-loop damping `&ocean_bt_nml bebt = 0.1` are the model
DEFAULTS since 2026-09-22 (MOM6 parity), so the template does not list them —
every cell runs them. Re-pinned on both toolchains with them in: every
viscous cell of `sigma`, `zstar`, `zstar_sigma` and `zstar_full` passes, rx0
0.1–0.8 included.

### Two legs

The legs share the template and differ **only** in the dissipation (the
self-test diffs every viscous namelist against its inviscid twin with the
`&ocean_hvisc_nml` / `&ocean_bdrag_nml` groups removed):

* **INVISCID** (`vcm_*`, 115 cells) — `nu_h = 0`, no drag. **The hard probe.**
  A viscosity that removes a spurious mode removes the measurement with it
  (`nu_h = 10 m² s⁻¹` costs the rx0 mode twelve decades), so this leg keeps it
  out and reports what the discretisation does on its own. Its
  terrain-following growth is **documented expected behaviour**: MOM6 shows the
  same mode (below). Its failures are XFAIL with that reason.
* **VISCOUS** (`vcmv_*`, 90 cells — the runnable ones; a refusal is a property
  of the configuration, and the `N² = 0` control belongs to the inviscid
  question) — MOM6's shipped seamount closure, translated. **This leg is the
  PASS/FAIL gate**: a cell inside its family's documented rx0 envelope must
  pass every assertion; `stability.py --self-test` refuses a marker on any
  in-envelope viscous cell.

### The viscous closure: MOM6's seamount, translated

MOM6's own seamount rest test (`ocean_only/seamount`, dev/gfdl `d74a11f9c`,
resolved `MOM_parameter_doc.all`) ships `LAPLACIAN = True`, `KH = 1000 m² s⁻¹`,
`KH_VEL_SCALE = 0.003 m s⁻¹`, `BOUND_KH = True`, `BIHARMONIC = False`,
`BOTTOMDRAGLAW = LINEAR_DRAG = True`, `CDRAG = 0.002`, `DRAG_BG_VEL = 0.05 m s⁻¹`,
`HBBL = 10 m`, `KV = 1e-4 m² s⁻¹`, on a 5 km grid at `dt = 900 s`.

| MOM6 | roundabout (viscous leg) | the mapping, and why |
|---|---|---|
| `LAPLACIAN`, `KH` | `nu_h = 160 m² s⁻¹`, `stress_tensor = .true.` | **Definition:** both are the `ν` of `∂u/∂t = ν∇²u`. MOM6's `diffu = (1/h)∇·(h·Kh·S)` with tension `S_xx = u_x − v_y` and shear `S_xy = v_x + u_y` (`MOM_hor_visc.F90`) reduces to `Kh∇²u` on a uniform grid with uniform `h` — the cross terms cancel — and so does roundabout's scalar path (`hvisc_compute_scalar_impl`, 5-point velocity Laplacian). They differ where `h` varies between neighbours: MOM6 thickness-weights and conserves momentum, the scalar path does neither. roundabout's `stress_tensor = .true.` **is** MOM6's operator (thickness-weighted stress, `÷(h_u + h_neglect)`, per-cell BOUND_KH clamp, coast masks), so the leg selects it. **Magnitude:** translated, not copied. What a Laplacian does to a grid-scale mode is damp it at `ν·k²_grid ∝ ν/Δx²`, and the modes this matrix exists for are 2–3 Δx wide; MOM6's `KH/Δx² = 1000/5000² = 4.0e-5 s⁻¹`, i.e. `160 m² s⁻¹` at 2 km. Copying `1000 m² s⁻¹` would be 6.25× MOM6's grid-scale damping — and MOM6 itself would not run it here: its BOUND_KH ceiling `0.1·Δx²/dt` is 667 m² s⁻¹ at (2 km, 600 s). |
| `KH_VEL_SCALE` | not set | `KH_VEL_SCALE·Δx` = 6 m² s⁻¹ at 2 km (15 at 5 km) is below `KH` on both grids, so it never binds. (roundabout's `kh_vel_scale` only seeds `ah_bg` for the flow-aware closures and is inert on the constant-`nu_h` path.) |
| `LINEAR_DRAG`, `CDRAG·DRAG_BG_VEL`, `HBBL` | `form = "linear"`, `r = 1.0e-5 s⁻¹`, `hbbl = 10 m` | MOM6's linear drag is a bottom **stress** `τ/ρ₀ = CDRAG·DRAG_BG_VEL·u_bbl = 1.0e-4 m s⁻¹ · u_bbl`, grid-independent, over the bottom `HBBL`. roundabout's distributed linear form applies `∂u_k/∂t = −r·u_k·h_in_bbl,k/h_k`, whose column integral is `r·hbbl·u_bbl`; `r = 1e-4/10` reproduces the stress exactly. |
| `KV = 1e-4` | not carried | its damping of the first baroclinic mode, `KV·(π/H)² ≈ 1e-9 s⁻¹`, is four decades under the rates measured here, and `use_closure = .false.` keeps `KD = 0`, which the tracer-extrema gate needs. |

`stress_tensor` is refused under `zfixed_closed_faces`, so the `z_fixed` cells
of the viscous leg carry the scalar operator (with the closed-face free-slip
masks). The self-test pins `nu_h = 160` and `r = 1e-5` to their derivation.

### The MOM6 baseline

The inviscid leg's terrain-following growth is **not a roundabout defect**.
MOM6 ocean_only (GNU build, dev/gfdl `d74a11f9c`, nothing modified), on its own
seamount at rest in **sigma** coordinates with every dissipation off and
`f = 1e-4 s⁻¹`:

| rx0 | dissipation | `f` [s⁻¹] | e-folding of KE [d] | R² | `En_KE` day 30 [m² s⁻²] |
|---|---|---|---|---|---|
| 0.000 | inviscid | 1e-4 | *no growth* | – | **0** (exactly) |
| 0.103 | inviscid | 1e-4 | **0.81** | 1.0000 | 2.16e-11 |
| 0.412 | inviscid | 1e-4 | **0.85** | 1.0000 | 2.00e-12 |
| 0.764 | inviscid | 1e-4 | **0.85** | 0.9999 | 2.90e-12 |
| 0.412 / 0.764 | `KH = 10` | 1e-4 | 0.95 | 1.0000 | 5.6e-14 / 4.0e-14 |
| 0.412 / 0.764 | inviscid | 0 | *no growth* | – | 5e-24 / 1e-22 |
| any | **shipped closure** | 1e-4 | *no growth* | – | ~1e-26 |

A textbook exponential out of round-off, rate independent of rx0, zero on a
flat bed, absent without rotation — the same mode the roundabout forensics
isolated (baroclinic, 2–3 Δx, pinned to the steepest face, dt-independent,
seeded by a 1.1e-16 m s⁻² PGF residual that is below one ulp of the
hydrostatic pressure). An isotropic 40 × 40 seamount (rx0 0.363) e-folds at
0.77 d, so it is not a slice artefact. Two differences are recorded rather than
explained away: `KH = 10` only shrinks the amplitude in MOM6 (35–500×) but
kills the mode in roundabout; and MOM6's z* and native-layer coordinates are
far worse on this test (O(1 m s⁻¹) within a day, which no MOM6 PGF or remap
option rescues at rx0 0.76). **Compare rates and controls, not amplitudes** —
an amplitude seeded by round-off is compiler-dependent, let alone
model-dependent. Recipes, JSON series and the full table live in the
prototypes repository (`mom6_baselines/README.md`); MOM6 is never a test
dependency.

### Geometries

| problem | geometry | `Δe` per face | bar |
|---|---|---|---|
| `flat` | flat bed, free surface | 0 | `REST_1UM_S` |
| `slope` | constant gradient, 1000 → 500 m over 48 cells (`topo_config="file"`) | 10.6 m | `REST_1MM_S` |
| `seamount_gentle` | Gaussian seamount, peak 300 m, `L = 40 km` = 20 cells | ~30 m | `REST_SEAMOUNT` |
| `seamount_steep` | the same, `L = 15 km` = 7.5 cells | ~75 m | `REST_SEAMOUNT` |
| `rx0_010` … `rx0_080` | two-level shelf/trough, ONE face at `rx0 = 0.1/0.2/0.4/0.6/0.8` | `2·H̄·rx0` = 150 … 1200 m | `REST_1MM_S` |
| `lid_flat` | flat ice lid over a flat bed (cavity) | 0 | `REST_1UM_S` |
| `lid_slope` | linearly sloping ice lid + calving front (cavity) | 13.8 m | `REST_1MM_S` |

`rx0 = |H_a − H_b| / (H_a + H_b)` is the terrain-following **stiffness**
(Beckmann & Haidvogel 1993, *J. Phys. Oceanogr.* **23**, 1736–1753, §2c, who
write it `r = |Δh|/(2h̄)`; Haney 1991, *J. Phys. Oceanogr.* **21**, 610–619
states the same thing as hydrostatic consistency). `rdb_ocean_stability_audit`
bounds it at **0.2** and warns above; the ladder walks *past* the bound so the
matrix can say where each family stops being trustworthy rather than asserting
the bound and hoping. The motivating failure — the ISOMIP+ trough sidewall that
takes the shipped case non-finite at day 3.2 — sits at **0.73**, between the
0.6 and 0.8 rungs. `tools/make_bathy_nc.py` dials it exactly:
`H_deep = H̄(1+rx0)`, `H_shallow = H̄(1−rx0)`.

### Axes

* **coordinate family** — all ten: `lagrangian`, `eulerian_z`, `sigma`,
  `zstar`, `zstar_sigma`, `zstar_full`, `z_fixed`, `rho`, `hycom`, `zsigma`.
  A family **refused by design** gets a row asserting the *refusal*, so an
  accidental un-refusal turns the cell `NOT-REFUSED` and says so.
  `eulerian_z` pins `split_scheme="ssp_rk2"` because the `pred_corr` v1
  envelope refuses it fail-loud — which also keeps that row out of the
  automatic scheme axis.
* **stratification** — `linear` (uniform `N² = 5.77e-6 s⁻²` from a linear
  `S(z)`) and `unstrat` (`N² = 0`, the control: the truncation
  `G(K) ∝ ρ₀N²Δe³` is not *small* but **absent** however steep the geometry).
* **EOS** — `linear` (ISOMIP+ Table 4) and `wright` on the clean slope.
* **outer split** — every non-pinning cell gets its `__ssp_rk2` twin from the
  existing scheme axis, for free.

## The metrics, with formulas

Everything comes from the model's own console output — no new Python
dependencies, no NetCDF reader.

| metric | assertion | formula / source |
|---|---|---|
| spurious energy **level** | `energy:rest` | peak `En` from `[stats]`, reported as `\|u\|_rms = sqrt(2·En)` |
| spurious energy **rate** | `energy:rest-growth-rate` | least squares of `log En` over the **tail** (last 60 %) of the run; `En ~ exp(2σt)`, so the reported amplitude rate is `slope/2`. **R² is reported, not gated** — a plateau has a near-zero slope and a meaningless R², an instability has both, and letting the reader see the pair is honest where picking an R² threshold would not be. Tier 1 only: the bar is a 20-day e-folding and a 3.3-day twin cannot fit one. |
| **truncation estimate** | printed in the row's note | `a_peak = N²·Δe³/(6·Δx·H̄)`, derived in [`validation_examples/ocean/ice_shelf_cavity/README.md`](../../validation_examples/ocean/ice_shelf_cavity/README.md) and verified there against a four-decade slope ladder to within 33 %. A rotating run balances it geostrophically at `U = a_peak/\|f\|`, so `En_plateau ≈ ½U²`. It is an **estimate** (one dominant face step, geostrophic balance) and is labelled as one; the bar never derives from it. |
| tracer bounds | `tracer:no-new-extrema` | `[diag]` `temperature` / `salinity` extrema may not leave the **first sample's** range by more than max(1e-6, 1.5 print quanta). The console line is `ES13.5`, so salinity near 34 is resolved to 1e-4 PSU: a one-last-digit flip is rounding of two printed numbers, not a measurement (it was the one gfortran-vs-V100 flip on `slope × zstar_full`). At rest, with no flux and no mixing, a new extremum is spurious diapycnal mixing from the coordinate's own regrid. |
| remap guard | `remap:preconditions` | the run never tripped `&vcoord_nml remap_check_preconditions` (a negative layer or a mismatched column total handed to the remap). Evaluated even when the run aborted — the guard's job is to abort. On nvfortran the guard's line can be split by the interleaved `ERROR STOP`; the parser records the firing either way. |
| thickness positivity | `thickness:positive` | `min h_layer` over the run, from the `h_layer` derived diagnostic (`&ocean_diag_nml diags = "h_layer"`), strictly positive. |
| budget residuals | `conserve:{Mass,Salt,Heat}` | the model's own closed-budget `Error`, 1e-11 relative. |
| clamp / truncation counters | `counters:no-truncation` | the solver's `CFL truncations: N total` line. A rest case that truncates is not resting, whatever its energy says. |
| stiffness | recorded in the JSON | the configure-time audit's own `rx0 = …` line. |

### The growth-rate bar, derived

`REST_SIGMA_MAX` is **0.05 per day of `En`** — a 20-day e-folding — expressed
in 1/s of amplitude. Derived, not tuned:

* the sloping-boundary sigma instability the matrix exists to find runs at
  `σ_En = 0.333/day` (3.0-day e-folding) on the measured cavity case, and never
  below `0.094/day` anywhere on its slope ladder → the bar sits **1.9× under
  the slowest measured instance of the defect**;
* the residual creep the *default* outer split leaves behind is `0.012/day`
  (83-day e-folding) → the bar sits **4× above what a healthy-but-imperfect run
  does**.

One decade between those two is what makes a single bar workable.
`stability.py --self-test` replays both measured series through the fit and
checks it fails the first and passes the second — who tests the test.

## THE BASELINE TABLE

**Tier 1 — the real horizon.** 205 cells (115 inviscid + 90 viscous), 30
simulated days each, `pred_corr`, `RDB_ENABLE_MPI=OFF`, measured on **both**
toolchains with the MOM6-parity defaults of 2026-09-22 (`bebt = 0.1`,
`renorm_consistent_flux = .true.`) AND the vanished-layer content rule I1′
(PR #50) together: gfortran 15.1.0 Release on the CPU (8 workers, 26 min)
and nvfortran 26.5 on four V100s, cc70 (one worker per device, 56 min — the
GPU is launch-bound at 48 × 6 × 15). `ok` = every
assertion passed (on `flat`/`lid_flat` that means `En = 0.000E+00` at every
sample, exactly); a number is the peak `En` in m² s⁻² of a cell that
completes but fails a gate (`leak` = the salt/heat budget left round-off);
`✗ dN` = aborted on day N (remap guard or non-finite); `refused` = rejected
at configure by design, and the row asserts the refusal. `a / b` =
gfortran / nvfortran where they differ; otherwise both.

**Both toolchains agree on PASS/FAIL on every one of the 205 cells** (the
two `slope × hycom × wright` cells used to differ — FINDING C, GPU only — and
agree since its fix, e8a1ab68d); where they differ it is in the number, the
day or the abort mode, never the verdict.

### inviscid leg

| problem (rx0) | lagrangian | eulerian_z | sigma | zstar | zstar_sigma | zstar_full | z_fixed | rho | hycom | zsigma |
|---|---|---|---|---|---|---|---|---|---|---|
| `flat` (0) | ok | ok | ok | ok | ok | ok | ok | ok | ok | refused |
| `slope` (0.00709) | ok | 1.1e-03 | ok | ok | ok | ok | 4.9e-07 | 1.4e-06 / 1.3e-06 | 1.7e-06 / 1.5e-06 | refused |
| `lid_flat` (0) | refused | refused | ok | ok | refused | refused | ok | refused | refused | refused |
| `lid_slope` (0.0138) | refused | refused | ok | ok | refused | refused | 1.5e-04 | refused | refused | refused |
| `seamount_gentle` (0.03) | ok | 2.8e-04 / 2.7e-04 | ok | ok | ok | ok | 1.2e-05 | ✗ d21 / ✗ d19 | ✗ d19 / ✗ d28 | refused |
| `seamount_steep` (0.078) | ✗ d4 | ✗ d18 / ✗ d15 | ok | ok | ok | ok | 1.2e-04 | ✗ d7 | ✗ d6 | refused |
| `rx0_010` (0.1) | ✗ d12 | 7.5e-04 | 2.1e-12 / 1.9e-12 | 2.1e-12 / 1.9e-12 | ok | ok | 8.1e-08 | ✗ d4 / ✗ d3 | ✗ d9 | refused |
| `rx0_020` (0.2) | ✗ d6 / ✗ d7 | ✗ d23 | 7.4e-04 / 7.2e-04 | 7.4e-04 / 7.2e-04 | 9.6e-04 / 9.7e-04 | 6.3e-05 / 6.8e-05 | 3.0e-07 | ✗ d1 | ✗ d1 | refused |
| `rx0_040` (0.4) | ✗ d5 / ✗ d4 | ✗ d11 | 3.7e-03 / 3.4e-01 | 3.7e-03 / 3.4e-01 | 5.5e-03 / 1.7e-02 | 3.0e-04 / 3.2e-04 | 1.8e-08 | ✗ d0 | ✗ d1 / ✗ d0 | refused |
| `rx0_060` (0.6) | ✗ d4 | ✗ d7 | 3.1e-07 / 1.6e-07 | 3.1e-07 / 1.6e-07 | 3.5e-07 / 5.3e-07 | 2.6e-09 / 1.8e-08 | 2.3e-07 | ✗ d0 | ✗ d0 | refused |
| `rx0_080` (0.8) | ✗ d5 | ✗ d4 | ✗ d17 | ✗ d17 | ✗ d17 | ✗ d14 / ✗ d13 | ok | ✗ d0 | ✗ d0 | refused |

### viscous leg

| problem (rx0) | lagrangian | eulerian_z | sigma | zstar | zstar_sigma | zstar_full | z_fixed | rho | hycom |
|---|---|---|---|---|---|---|---|---|---|
| `flat` (0) | ok | ok | ok | ok | ok | ok | ok | ok | ok |
| `slope` (0.00709) | ok | ok | ok | ok | ok | ok | 4.4e-08 | ✗ d0 | ✗ d0 |
| `lid_flat` (0) | - | - | ok | ok | - | - | ok | - | - |
| `lid_slope` (0.0138) | - | - | ok | ok | - | - | 5.3e-06 | - | - |
| `seamount_gentle` (0.03) | ok | ok | ok | ok | ok | ok | 1.2e-07 | ✗ d0 | ✗ d0 |
| `seamount_steep` (0.078) | ✗ d16 | ok | ok | ok | ok | ok | 8.0e-08 | ✗ d0 | ✗ d0 |
| `rx0_010` (0.1) | 1.2e-09 / 4.7e-09 | ok | ok | ok | ok | ok | ok | ✗ d0 | ✗ d0 |
| `rx0_020` (0.2) | ✗ d15 | ok | ok | ok | ok | ok | 1.4e-08 | ✗ d0 | ✗ d0 |
| `rx0_040` (0.4) | ✗ d9 | ok | ok | ok | ok | ok | 2.4e-10 | ✗ d0 | ✗ d0 |
| `rx0_060` (0.6) | ✗ d7 | ok | ok | ok | ok | ok | 1.2e-08 | ✗ d0 | ✗ d0 |
| `rx0_080` (0.8) | ✗ d7 | ok | ok | ok | ok | ok | 7.5e-09 | ✗ d0 | ✗ d0 |

The EOS and stratification controls (not in the grids above): `slope × sigma
× wright` passes both legs on both toolchains (inviscid En 7.9e-23, viscous
5.2e-24); `slope × z_fixed × wright` behaves exactly like its linear twin
(4.9e-07 / 4.4e-08, salinity overshoot only — no leak under I1′);
`slope × hycom × wright` completes the inviscid leg on both toolchains (En
2.6e-06 / 2.9e-06, failing the rest gates like its linear twin; with the
parity defaults but WITHOUT I1′ it aborted on the remap guard at day 29.5 /
24.0, and it died at step 0 on the GPU before the FINDING C fix) and aborts
at step 0 in the viscous leg (FINDING A). The `N² = 0` controls:
`lid_slope × sigma` holds 6.5e-20 / 2.1e-19 (machine zero, both toolchains);
`rx0_080 × sigma` unstratified still reaches the CFL wall (day 22 on both)
through the slower barotropic residual the forensics measured.

**Counts, tier 1** (identical on both toolchains, against the re-pinned
markers): inviscid **33 PASS / 82 XFAIL** (23 of them the by-design refusals)
/ 0 FAIL / 0 XPASS; viscous **58 PASS / 32 XFAIL** / 0 FAIL / 0 XPASS (was
41 / 49 with I1′ alone: the 17 FINDING B cells now pass on both toolchains).
**Tier 2** (the CI slice, 3.33 days, 120 rows with the `__ssp_rk2` twins):
**90 PASS / 30 XFAIL** / 0 FAIL / 0 XPASS on both toolchains.

**How to read it.**

1. **The fixed configuration works where it was aimed.** Every
   terrain-following family rests at machine zero on the slope, both
   seamounts and the sloping ice lid in BOTH legs (`En ≤ 9.1e-21` in the 2026-09-23 re-pin), where the
   pre-fix table carried 3e-09 … 2.7e-04 plateaus. The inviscid growth that
   remains lives only on the single-face `rx0` ladder of the INVISCID leg:
   sigma/zstar at rung 0.1 reach 2e-12 but still fit a growth rate over the
   20-day bar, rungs 0.2–0.6 now complete 30 days with it, and 0.8 reaches
   the CFL wall — the mode MOM6 shares. The VISCOUS leg rests on every rung
   (FINDING B, fixed by the 2026-09-22 defaults).
2. **`sigma` ≡ `zstar`** to every printed digit, in both legs, on both
   toolchains; `zstar_sigma` now differs from them on the ladder (it is the
   only place its z* branch engages).
3. **`z_fixed` is clean only without fillers** (`flat`, `lid_flat`). Wherever
   layers vanish it used to leak salt and heat at 1e-7 … 1e-6 relative; the
   vanished-layer content rule I1′ closed that to round-off in every cell.
   What remains is `tracer:no-new-extrema` (salinity overshoot from the
   regrid) and, under the sloping lid, the staircase residual — which under
   I1′ is 5–12× lower but approaches its bounded ceiling slowly and
   monotonically, so it also fails `energy:rest-settles` at 30 days
   (viscous: 5.24e-06 d30 → 5.86e-06 d60 → 6.19e-06 d90).
4. **`rho` / `hycom` are clean on a flat bed only.** Inviscid they abort on
   every sloping geometry within days (the negative layer the rest-state
   mode writes); viscous they abort at step 0 (FINDING A).
5. **`eulerian_z` (which pins `ssp_rk2`) rests to the top of the viscous
   ladder** — `rx0 = 0.8` with the closure (0.6 before the 2026-09-22
   defaults) — while its
   inviscid leg is the worst geometric family on a slope (1.1e-03).
6. **`lagrangian`** inherits the sigma-shaped initial column and never
   regrids; over a step it collapses a layer in both legs.

## The rx0 ENVELOPES — what v0.1.0 claims

A family's **envelope** is the largest geometry `rx0 = max|H_a − H_b| / (H_a +
H_b)` (over every wet-wet face; exact for the ladder, from the Gaussian formula
for the seamounts, `Δe/(2H̄)` for the slope and the lids — `geometry_rx0` in
`vcoord_matrix.py`) at which **every viscous cell of the family at or below it
passes every assertion on every toolchain measured**. It is derived by
`vcoord_matrix_pin.py`, pinned in `vcoord_matrix_measured.py`, and gated: the
self-test refuses a marker on any viscous cell inside its envelope. The
measured geometries are rx0 = 0, 0.0071 (slope), 0.0138 (sloping lid), 0.030
(gentle seamount), 0.078 (steep seamount), 0.1, 0.2, 0.4, 0.6, 0.8.

| family | viscous envelope (gate) | first viscous failure | inviscid: all-assertions pass up to | inviscid: completes 30 d at |
|---|---|---|---|---|
| `sigma`, `zstar` | **rx0 ≤ 0.8** (top of the ladder) | none measured | 0.078 | ≤ 0.6 |
| `zstar_sigma` | **rx0 ≤ 0.8** (top of the ladder) | none measured | 0.1 | ≤ 0.6 |
| `zstar_full` | **rx0 ≤ 0.8** (top of the ladder) | none measured | 0.1 | ≤ 0.6 |
| `eulerian_z` (ssp_rk2) | **rx0 ≤ 0.8** (top of the ladder) | none measured | flat only | ≤ 0.03, and 0.1 |
| `lagrangian` | **rx0 ≤ 0.03** | seamount_steep (collapse) | 0.03 | ≤ 0.03 |
| `z_fixed` (closed faces) | **flat only** | slope (salinity overshoot) | flat only | every geometry |
| `rho`, `hycom` | **flat only** | slope (FINDING A) | flat only | ≤ 0.0071 |

Measured 2026-09-23 with the MOM6-parity defaults (`bebt = 0.1`,
`renorm_consistent_flux = .true.`), gfortran 15.1 CPU and nvfortran 26.5
V100, identical viscous verdicts on both. The classical Beckmann–Haidvogel
bound is 0.2; the terrain-following viscous envelopes now reach the top of
the measured ladder (0.8). Before those defaults they stopped at 0.078 / 0.1,
capped by FINDING B. The INVISCID leg is unchanged in its all-gates envelope
(the MOM6-shared rest-state mode, not FINDING B), but now COMPLETES 30 days up
to rx0 0.6 in every terrain-following family (0.4 used to abort); 0.8 still
aborts on the remap guard.

## Per-toolchain tolerance

A marker pins **which assertions** a cell fails, never a number: its assertion
list is the UNION over gfortran and nvfortran and over both tiers (the
tolerance band), and its reason quotes every toolchain's own number. A cell
whose number drifts inside its band stays XFAIL on both toolchains; a cell
that fails a NEW assertion on either turns FAIL; a cell that passes outright
on one toolchain and not the other is marked `toolchain_dependent` and
reports PASS (not XPASS) where it passes. The bars are physics bars and are
never read from a run.

What differed on the previous measurement (five cells between gfortran and
V100) and on this one:

* **the abort mode, not the verdict** — which of `finite` /
  `remap:preconditions` reports a blow-up that both toolchains have
  (`rx0_080 × {zstar_sigma, sigma/unstrat}`, `seamount_steep × lagrangian`);
  the union covers both;
* **`tracer:no-new-extrema` on `slope × zstar_full`** — one last printed digit
  of salinity (3.46700E+01 vs 3.46701E+01). That was a harness artefact: the
  gate's 1e-6 slack is below the 1e-4 resolution of the console `[diag]`
  line it reads. The gate now allows 1.5 print quanta (`diag_print_quantum`
  in `stability.py`; self-tested both ways), which removes the flip and
  leaves every real overshoot (z_fixed: 7e-4 … 0.14 PSU) failing;
* **FINDING C** — the only genuine verdict difference, GPU only; FIXED
  (e8a1ab68d), the toolchains now agree on every cell.

## FINDINGS — what the two legs turned up

Three sets of cells FAIL where one would expect them to rest — MOM6's
shipped closure rests its seamount at every steepness it can reach
(rx0 ≤ 0.76), and roundabout's own inviscid leg *completes* several of these
cells. They were localised to a first-order term, **not** tuned away: the
envelope table below is capped by them, and each marker carries its finding
by name (`finding_visc_pred_corr`, `finding_stress_density` in
`vcoord_matrix.py`; FINDING C's `finding_gpu_wright_density` was retired with
its fix). Nothing in the Fortran was changed by this PR.

### FINDING B — `pred_corr` × Laplacian viscosity destabilises a stepped terrain-following column

**Where.** `sigma`/`zstar` at every ladder rung 0.1–0.8 (the inviscid leg
completes 0.1, 0.2 and 0.6), `zstar_sigma` at 0.2/0.4/0.8, `zstar_full` at
0.2–0.6. Smooth geometries (slope, both seamounts, the sloping lid) rest at
round-off (`En ≤ 4.5e-23`) in every terrain-following family.

**First failure** (`vcmv_rx0_010_sigma`, gfortran, field output every 3 h; the V100 aborts the same cell on day 19):
En creeps at round-off (`|u| ~ 3e-12 m s⁻¹`, all of it in the **top layer at
the step face**, `i = 24–25`, baroclinic: column mean 40× under the deviation)
until ≈ day 11.5 (outer step ≈ 1660); then a **BAROTROPIC, grid-scale (2 Δx),
domain-wide** mode takes over (column mean ≈ max|u|, deviation 40× smaller)
and grows ×1e4 in 30 h (amplitude rate ≈ 8.5e-5 s⁻¹), until continuity writes
a negative layer (−255 m at step 2087) and the remap guard stops the run.
On `rx0_060` with the scalar operator the same sequence starts at day 8.5.

**The term.** One knob at a time against the gate configuration
(`vcmv_rx0_010_sigma`, 30 days, gfortran):

| substitution | outcome |
|---|---|
| (gate: `nu_h = 160`, stress tensor, linear drag, `pred_corr`, `dt = 600`) | **aborts day 14** |
| no viscosity and no drag (= the inviscid cell) | rests, En 1.9e-13 |
| no drag (viscosity alone) | aborts day 15 |
| `split_scheme = "ssp_rk2"` | **rests**, En 1.1e-24 |
| `f = 0` | **rests**, En 2.2e-23 |
| `dt = 300 s` | **rests**, En 5.2e-24 |
| `dt = 360 s` (MOM6's dt/Δx) / `450 s` | aborts day 20 / 18.5 |
| `nu_h = 10` | **rests**, En 9.3e-20 |
| `nu_h = 40` | aborts day 17 |
| `coriolis form = "sadourny_energy"` (MOM6's seamount form) | aborts day 14.5 |

and, on `rx0_060` with the scalar operator and no drag (14.6-day probes;
"onset" = first sample to jump ×100): `nu_h` 40/80/160 onset d6.7/d6.5/d8.8;
`bound_kh` no effect (160 is under its 333 m² s⁻¹ ceiling);
`remap_boundary_extrap = .false.` **earlier** (d3.8); PCM remap much earlier
(d2.1); `pc_be = 1` earlier (d5.0); `remap_nonuniform_weights = .false.` /
`correction_h_weighted` later (d10.8 / d10.2); y-periodic earlier (d2.3, so
the walls are not the cause); `nghost = 3` bit-identical; a frozen regrid
(`regrid_time_scale = 10 d`), `N² = 0`, `f = 0`, `ssp_rk2` and `dt = 300/450`
clean to the probe horizon.

**Localised (`fix/pred-corr-viscous-bt-mode`).** The table above is a
table of SEEDS, not of stability conditions. Seeding the same cell with a
1e-6 m s⁻¹ barotropic velocity pattern makes EVERY row of it blow up within
2–3 days — `nu_h = 0` (no drag either), `nu_h = 10`, `f = 0`, `dt = 300`,
`dt = 360`, `N² = 0`, ONE layer — except `ssp_rk2` and `bebt ≥ 0.05`, which
decay it. A flat bed seeded the same way is exactly neutral. The viscosity
is not in the mechanism, and neither is the frozen viscous forcing: under
`pred_corr` the viscous term acts on the barotropic mode only through the
time-mean `u_av`, so a grid-scale gravity mode is damped per step by the
factor `1 − λ·dt·sinc²(ω·dt/2)` (`λ = ν·k²`), which can never exceed one —
at 2 km / 600 s, `λ·dt = 0.096` but `sinc² = 1.3e-3`, i.e. the Laplacian
barely touches that mode at all (measured on a seeded flat-bed channel:
`KE+PE` 0.71 with `nu_h = 160` vs 0.83 with none after 1000 steps). The anti-damping
the `bound_kh` docstring records needs the forcing on `u^n` (`ssp_rk2`,
factor `1 − λ·dt·sin(ω·dt)/(ω·dt)`), and even there it is `≤ λ/(ω·dt) =
0.019 λ = 3e-6 s⁻¹`, about thirty times under the measured growth.

The defect is in the barotropic transport renormalisation
(`renormalise_zonal_flux_to_uhbt`, continuity). Its Newton solve re-picks
each layer's upwind donor at the corrected velocity but keeps the OLD
donor's `u0·h_old` in the flux, `flux0 + du·h_new`; at a flip that model
jumps by `u0·(h_new − h_old)·dy`. Over the step both PPM edges are
flattened to the cell thicknesses (825 | 675 m), so when the corrector's
layer velocity (the END-of-step barotropic velocity) and the target `uhbt`
(the TIME-MEAN transport) have opposite signs and
`|ū_mean| < |u_end|·Δh/h_new` (Δh/h = 0.18 at rx0 0.1, while the time mean
of a 2Δx gravity mode is ~`sinc(ω·dt/2) ≈ 4 %` of its end value — true about
every other step), the solve has no root, cycles for its 8 iterations and
hands continuity a layer transport of the WRONG SIGN (measured −3.8e-2 vs
+2.9e-3 m³ s⁻¹). The layer free surface then leaves the barotropic `η_end` by
O(η) in the two step columns — a spurious η dipole at every such step that
pumps the undamped barotropic grid-scale mode (8 % per step in amplitude at
rx0 = 0.1 and 21 % at 0.3 on one layer, at `f = 0`, inviscid, no remap;
8.5e-5 s⁻¹ in the matrix cell). MOM6's `zonal_flux_adjust` recomputes each
layer's flux at `u + du` with its own donor (continuous) and brackets Newton
with bisection, so it cannot do this.

**Fix:** `&ocean_continuity_nml renorm_consistent_flux = .true.` — the
DEFAULT since 2026-09-22 (MOM6 parity), together with MOM6's `&ocean_bt_nml
bebt = 0.1`; faces where no donor flips are byte-identical to the historical
model. Gates: `test_continuity_multilayer`
(`renorm_donor_flip_lands_on_uhbt_x`/`_vhbt_y`,
`renorm_consistent_flux_no_flip_bit_identical`) and `test_ocean_step_bt_mode`
(a 48-cell rx0 0.1 rotating channel under `pred_corr`: layer-vs-barotropic
`η` 7.2e-4 m and `KE+PE` ×670 in 150 steps without the knob, 2.0e-13 m and
×1.08 with it). `vcmv_rx0_010_sigma` with the knob rests 30 days (En
1.4e-24); seeded as above it stays neutral. The re-pin with both defaults on
(2026-09-23, both toolchains) turned every FINDING B marker XPASS; none is
left.

### FINDING A — the stress-divergence viscosity drives a density-space column negative

**Where.** `rho` and `hycom`, every non-flat geometry, both outer schemes, at
outer steps 4–8 (−5e-5 m … −0.3 m). Flat geometries are clean.

**The term.** Against `vcmv_slope_rho`:

| substitution | outcome |
|---|---|
| gate (`stress_tensor = .true.`, drag) | remap guard, step 7 |
| no drag | still fails (remap guard, step 7 — identical) |
| scalar operator (`stress_tensor = .false.`), with drag, 30 d | **rests**, En 2.0e-23 (hycom: 7.0e-25) |
| no viscosity | clean over the 0.4 d probed (the inviscid cell completes 30 d, En 2.3e-06) |

So it is the thickness-weighted operator on the thin, keep-alive layers a
density coordinate carries where its targets outcrop. The one MOM6 thin-layer
safeguard roundabout's port lacks is `hrat_min = min(1, h_min/(h + h_neglect))`
scaling the BOUND_KH ceiling (`MOM_hor_visc.F90`) — **a hypothesis**, not
measured. **What a user can do today:** the scalar operator.

### FINDING C — GPU only: the density-space target builder faults under the Wright EOS (FIXED, `e8a1ab68d`)

**Status: FIXED** by e8a1ab68d (`fix(vcoord): rho/hycom x Wright illegal
address on the GPU build`): compute-sanitizer put the fault at a HOST address
read through an `associate` name passed to the out-of-module EOS routine, and
the target builders no longer wrap a `do concurrent` in `associate`.
Re-measured 2026-09-23: every `rho`/`hycom` × Wright cell now runs on the GPU
and reaches the same verdict as on gfortran. The record below is kept as the
history of the finding.

**Where.** nvfortran 26.5 / V100 (cc70), every `rho` or `hycom` run with
`&ocean_eos_nml eos = "wright"`: the matrix's `slope × hycom × wright` in both
legs, and — probed on purpose — `flat × hycom × wright` (a flat bed at rest)
and `slope × rho × wright`. gfortran completes the same inviscid cell to day
30 (En 2.2e-06).

**First failure.** The first regrid (no `[stats]` line is ever printed):
`CUDA_ERROR_ILLEGAL_ADDRESS` in `ocean_vcoord_compute_target_h_rho_impl`
(`src/core/ocean/vcoord/rdb_ocean_vcoord.F90`, the column `do concurrent`).
The linear EOS on the same cells is clean on the GPU; a gfortran
`-fcheck=all` Debug build runs `flat × hycom × wright` with no bounds
violation. The Wright coefficients are `parameter`s, so the point EOS call is
not the suspect; the kernel wraps its `do concurrent` in an `associate` over
`this%` components — the NVHPC mapping hazard CLAUDE.md records — which is a
hypothesis, not a localisation. No `ctest` case runs `rho`/`hycom` × Wright
on the device, which is why the V100 suite is green. Related, and the same
kernel: when FINDING A hands it a negative layer the GPU faults there too,
before the remap guard can report — the non-monotone-column refusal the audit
proposed (fix item 2) would turn both into fail-loud messages.

## Adding a family or a problem

* **A new coordinate family** — append it to `FAMILIES` in `vcoord_matrix.py`
  with its `vcoord_type`, its status (`run` or `refused`), any extra
  `&vcoord_nml` lines it needs, and a one-paragraph note saying what the row
  *measures*. If it is refused under a cavity, it also belongs outside
  `CAVITY_ACCEPTED`. [`docs/howto/add_vertical_coordinate.md`](../../docs/howto/add_vertical_coordinate.md)
  makes this a required step: a family that is not in the matrix has no
  envelope statement. Both legs pick it up automatically.
* **A new problem** — append to `PROBLEMS` with its `class`, its `bar`
  constant, its worst per-face `Δe`, its `topo` configuration, and `tier2` /
  `tier2_viscous` only if it belongs in the CI slice of that leg.
* **Then re-measure, on both toolchains, and pin by tool.** Sweep tier 1 and
  tier 2 on gfortran and on nvfortran (`--tags vcoord_matrix --out …`,
  `--scratch-root` distinct per toolchain if they run at once), then
  `python3 tests/regression/vcoord_matrix_pin.py --t1 gfortran=… --t1
  nvfortran=… --t2 gfortran=… --t2 nvfortran=… --provenance "…"`. The
  generated `vcoord_matrix_measured.py` is never hand-edited.
* **Never** close a cell by widening `en_rest_max`, by shortening the run, or
  by putting viscosity into the INVISCID leg — and a viscous cell that fails
  where you expected it to rest is a FINDING to localise, not a marker.

## Tiers, cost, and what runs in CI

| | tier 2 (CPU, **the CI gate**) | tier 1 (local only) |
|---|---|---|
| inviscid leg | `flat`, `slope`, `rx0_060` × every family | all 115 cells |
| viscous leg | `flat`, `slope`, `seamount_steep` (the steepest seamount; the rx0 0.1–0.8 ladder is tier-1 only) × every family | all 90 cells |
| length | 480 steps = 3.33 simulated days | 4320 steps = 30 simulated days; 1440 = 10 days for the `seamount_*` and `rx0_*` problems (maintainer decision 2026-09-25: under the MOM6 barotropic split the inviscid sloping-boundary mode MOM6 shares reaches the CFL wall inside 30 days on the steeper rungs, so 30 days only measured when it aborts; `T1_STEPS_SEAMOUNT` in `vcoord_matrix.py`) |
| rows | 120 (incl. the `__ssp_rk2` twins) | 205 (twins are tier-1-skipped by the curated scheme axis) |
| measured wall | 443 s on 2 workers of a contended 4-core box (gfortran); 1081 s on one V100 | 3438 s on 3 CPU workers (gfortran); 11 910 s on one V100 (nvfortran 26.5, launch-bound on this box) |
| the rate gate | `SKIP … TIER 1 ONLY` — the bar is a 20-day e-folding | evaluated |

**The decision.** `.github/workflows/ocean-stability.yml` runs on a hosted
4-core gfortran runner under a 45-minute job budget, at two cadences
(maintainer decision, 2026-09-22): every pull request runs the tier-2 suite
WITHOUT the matrix slice (`--exclude-tags vcoord_matrix`, 107 cases), and a
NIGHTLY scheduled run ON MAIN (plus `workflow_dispatch`) runs the full
corpus, matrix slice included (227 cases). The viscous slice is chosen to
carry the gate's IN-ENVELOPE geometries — the steep seamount is the steepest
of them that fits the tier-2 budget (the rx0 0.1–0.8 ladder, in-envelope
since the 2026-09-22 defaults, is tier 1) — so a regression that breaks the
v0.1.0 claim reddens the nightly run on main within a day. The full matrix is
tier 1: it needs 30 simulated days per cell to see a growth rate (the viscous
leg's FINDING B, before its fix, fired between day 2 and day 28), which
neither fits the budget nor runs on a
hosted runner's hardware for the GPU half; it is re-run locally, on BOTH
toolchains, whenever a vertical-coordinate, remap, PGF, viscosity or
split-step change lands, and the pinned table is regenerated from it.

## Not done, and not stubbed

* **The NONLINEAR (exponential thermocline) stratification axis.**
  `&ocean_zinit_nml source="linear"` is affine in `z` by construction and
  `source="file"` aborts at `validate_config` (PR-23b), so there is no way to
  lay a non-affine rest state today without a new analytic profile or the file
  reader. The axis is absent rather than faked.
* **The ADJUSTMENT problems and the RPE (spurious-mixing) diagnostic** — lock
  exchange, internal seiche, overflow, and Winters et al. (1995,
  *J. Fluid Mech.* **289**, 115–128) sorted-density background potential
  energy as the per-family mixing metric (Ilicak, Adcroft, Griffies &
  Hallberg 2012, *Ocean Modelling* **45–46**, 37–58; Petersen, Jacobsen,
  Ringler, Hecht & Maltrud 2015, *Ocean Modelling* **86**, 93–113). Not
  started. The rest matrix measures what a coordinate does to a fluid that
  should not move; RPE measures what it does to one that does. They are
  different questions and the second is a separate piece of work — it needs a
  cadence-bounded host-side sort and a hypsometric fill, i.e. a Fortran
  diagnostic, not a post-processor.

## Downscaling rules (`downscale.py`)

**A tier-2 twin is built by preserving the dimensionless numbers, never by
shrinking the grid.** Naive shrinking changes every resolved scale at once and
is exactly how the 2026-09-11 defects hid. The rules below are **checked for
every twin before it runs**; a violating twin fails the suite as a manifest bug
rather than quietly testing different physics.

| rule | statement | why |
|---|---|---|
| **R1** Munk layer | `nu_h >= beta * (2*dy)**3` | under two cells the western boundary current is unresolved and you get grid-scale wall noise that looks like physics. At −45 S with `dy = 0.5°` that floor is 2.2e4. |
| **R2** `ah_max` clamps `nu_h` | `ah_max >= nu_h` | `ah_max` is a magnitude **ceiling** on the assembled viscosity. Raising `nu_h` to satisfy R1 without raising `ah_max` is **inert** — the clamp throws the increase away. |
| **R3** viscous CFL | `nu_h*dt/dx² <= 0.125`, or `bound_kh = .true.` | explicit Laplacian friction is forward-Euler. `ah_max` does **not** help here: it bounds magnitude, not CFL. `bound_kh` is the per-cell CFL clamp. |
| **R4** deformation radius | `Rd/dx >= 4` for an eddying case | below ~2 cells per `Rd` no eddies form at **any** run length. |
| **R5** available potential energy | an eddying twin must declare its APE source | a uniform-density IC produces no eddies at any resolution — `acc_channel` as shipped never would have. |
| **R6** `nghost` vs scheme | weno5 → 3, weno7 → 4, periodic-x → 3 | a downscale must not quietly drop `nghost` to save cells. |
| **R7** domain / `Rd` | `min(L) >= 8*Rd` for an eddying case | a smaller box measures the box, not the turbulence. |

Two shapes of downscale follow from those rules:

- **A periodic channel is SHORTENED** — fewer cells, *same* `dx`. Every
  dimensionless number is preserved exactly; only the number of wavelengths the
  box holds goes down, which is what R7 bounds. This is the preferred recipe
  (`baroclinic_2layer`, `baroclinic_15layer`, `bc_inst_tuned_512`).
- **A closed basin can only be COARSENED** — the basin *is* the case — which
  moves `Rd/dx` and the Munk width in cells by the same factor. Acceptable for
  a quiescent case with no boundary current and no eddy field (the seamount
  family, 96×96 @ 4 km → 48×48 @ 8 km, a geometry the repo already ships as
  `seamount_conservative_floor`). **Fatal for an eddy-resolving case**, which
  is why several are tier-1 only.

**Where a case cannot be downscaled without changing its physics it stays
tier-1 only, and says which rule blocks it.** That is a valid outcome, not a
gap. Thirteen cases are in that position:

- `acc_channel_eddy`, `eddy_test`, `eddy_test_quick` — eddy-resolving; the only
  affordable downscale halves `Rd/dx` and breaks **R4**.
- `acc_channel_kitchensink_xl`, `seamount_bench_full` — their own coarse twin is
  already in tier 2.
- the eight `bench_scaling/*` benchmarks — the problem **size** is what they
  measure, so a downscaled twin measures nothing. They still get the tier-1
  stability gate.

## Gate E6 — ISOMIP+ Ocean0 melt-off under `z_fixed`, 180 days (tier 1 only)

Two tier-1 rows run `isomip_plus/ocean0_idealised_zfixed.nml` with basal melt
and top drag OFF (run-time `t1_overrides`, the committed file is untouched)
for 180 simulated days — the horizon is load-bearing, because the numerical
calving-front regime does not appear until day ~105:

```bash
# NVHPC GPU build, one V100 (~18 min per leg, sequential on one device)
CUDA_VISIBLE_DEVICES=2 python3 tests/regression/stability.py --tier 1 \
    --backend gpu --build-dir build_cc70 --gpus 2 -v \
    --cases isomip_plus_ocean0_zfixed_meltoff,isomip_plus_ocean0_zfixed_meltoff_nu30
```

| row | gate | measured (V100, 2026-09-24) |
|---|---|---|
| `isomip_plus_ocean0_zfixed_meltoff` (`nu_h = 6`) | `En(30 d) < 1E-07` | `5.632E-08` |
| | 5-day log-rate d25–30 < d15–20 | `0.0507` < `0.0664 /day` |
| | `En(180 d) < 1E-05` (and peak, `energy:rest`) | `3.147E-06` (peak `3.804E-06`, d173.5) |
| | `En(180)/En(150) < 2` | `1.39` |
| `isomip_plus_ocean0_zfixed_meltoff_nu30` | peak `En < 5E-07` (`REST_1MM_S`) | `2.121E-07` |

Measured on the v0.1.0 defaults (`bebt = 0.1`, `renorm_consistent_flux`, the
I1′ vanished-layer rule); within 3 % of `cavity_rest_growth_diagnosis.md` §Q
through day 90; regime 2 arrives on the same schedule (`~d105`) and
saturates at the same level (peak `3.804E-06` vs `3.704E-06`), but ~15 days
later (growth stops by `~d165`, not `~d150`). Zero `[nan-catch]`, `MaxCFL ≤ 0.0093`, budgets
`≤ 8.8E-12` over the 180 days on both rows.
The `nu_h = 30` row carries a scoped tier-1 XFAIL on `energy:rest-settles`
only (regime 1 is an `En ∝ t` boundary current that never settles; its last
sample is its peak). The `nu_h = 6` row carries NONE: regime 2 oscillates
about its saturation level and ends at 83 % of its peak, so the settle gate
passes — it reads the phase of that oscillation at day 180 (it missed by a
hair, 95.3 %, before `bebt = 0.1`), which is why a trip there is a thing to
localise, not to re-mark.
No tier-2 twin: see the `t2_reason`s.

## Adding a case to the stability suite

Add an entry to `stability_manifest.py`:

1. Pick the **regime** — `rest` / `adiabatic` / `baroclinic` / `forced`. This
   chooses the energy assertion, and it is the decision that matters: "energy
   must not grow" is only a correct statement for a case with nothing driving
   it.
2. Set `tier1` `n_steps` **from when the case's physics appears**, not from a
   global number. If it is an instability case, that is several e-folding times.
3. Set `tier2` `n_steps` to fit a few seconds of gfortran CPU, plus
   `overrides` + `dimensionless` if the twin changes the grid — then run
   `python3 tests/regression/downscale.py` to check it.
4. If the case's header makes a testable claim, encode it as a `claim`. **Do
   not** edit the namelist to make a number match; if it fails, that is the
   finding.
5. If it fails today for a known reason, add `known_failure` with a one-line
   reason and a pointer. Do **not** silence it by weakening an assertion.

## The known-failing namelists (as of 2026-09-13)

Six cases fail today (nine on 2026-09-12; the three `eady` cases were fixed). Each is marked `known_failure` in
`stability_manifest.py` with a reason and a pointer, so CI red stays
informative — **none of them is silenced by weakening an assertion or editing a
namelist.** They are grouped by the defect they expose.

### 1. `coriolis_coast` — the wind is not wired up (`energy:responds`)

The file sets `&ocean_topo_nml taux_magnitude = 0.05` under
`wind_config = "constant"`, where that knob is **inert** — `rdb_config`
documents `taux_magnitude` as *"peak zonal wind stress for
`wind_config = 2gyre/neverworld2`"* — and leaves `&physics_nml
wind_stress_x = 0.0`, which is the knob the constant path actually reads. The
case therefore runs its full 10 days at `En = 0.000E+00`, and has never tested
the coastal Coriolis energy budget its header describes ("an eastward wind
0.05 N/m² spins up a gyre that wraps the island and builds relative vorticity at
the free-slip corners"). Verified by re-running with `&physics_nml
wind_stress_x = 0.05`: `En` reaches 1.3e-6 and `MaxCFL` 2.6e-4 within half a day.

This is the class of defect nothing else can see: a forcing path that is
silently a no-op is perfectly conservative, perfectly finite, and matches its
golden exactly.

### 2–4. `eady`, `eady_weno5`, `eady_weno7` — over-damped — **FIXED 2026-09-13**

`eady` shipped `dT_dy = -5e-6` (a front 4× weaker than the one its header's
`Ri ~ 156` / `sigma ~ 2.5e-6` / "1200×" described) with `nu_h = 100`, and
**decayed everywhere** (sigma −2.8e-7 1/s over 25 days). The two WENO variants
are the same basin and decayed identically.

Re-baselined to the front the header always described — `dT_dy = -2e-5`,
`nu_h = 20`, 60 days — after a `nu_h × dT_dy` sweep with modal spectra
(`validation_examples/ocean/eady/eady.nml` header carries the table). Two
things the sweep overturned:

* the earlier "`nu_h = 10` recovers Eady growth for the weak front to within
  31 %" was **not the Eady mode**: at `nu_h ≤ 10` every zonal wavenumber
  kx = 1..8 grew at the *same* rate, faster for a *weaker* front, and the
  growth is still there with no front at all (`dT_dy = 0`, resting stratified
  channel, `nu_h = 0`: `max|v|` 0 → 2 cm/s in 25 days; 1.5 mm/s at
  `nu_h = 20`). That is a numerical growth of the resting state (open item;
  documented in the nml header) and the reason any growth in this basin must
  be checked for wavenumber selectivity before it is called Eady;
* at `dT_dy = -2e-5` the kx = 2 (125 km) mode is unambiguous — it grows
  ~22000× while kx ≥ 5 decay — at 1.93e-6 1/s against 1.98e-6 from channel
  theory, insensitive to `nu_h` from 10 to 100.

The gates now read `max|v|` from the `[diag]` line (`sigma_expected` with
`"field": "v"`, fit from day 30; band [1.5, 2.5] e-6) plus a `field_ratio`
claim on the header's "~1500×" (bar 300×), because `En` — 4.3e-3 m²/s² of
basic-state jet, drifting down 6 % — cannot see the mode until day ~50. All
three pass both tiers. `sigma_fast_max` is pinned at 6e-6 (the default
derivation from the wider band would land at 1.25e-5 and miss the 1.0e-5
thermo-cadence instability the WENO cases guard against).

Not fixed, and outside the shipped window: the front's nonlinear collapse past
day ~65 hits the 2 m/s clamp and NaN-catches by day ~75 at dt = 600 s
(`n_inner` 20 and 30, `maxvel = 10`, `nu_h = 100` all fail; dt = 300 s survives
to day 100 but with `En` above the front's available potential energy). The
60-day `t_end` is that boundary, stated.

### 5–8. `acc_channel_{weno5,kitchensink,sw_penetration,kitchensink_xl}` — heat leak (`conserve:Heat`)

The model's own closed Heat budget residual reaches `|Error| ~ 5e-5` relative —
10⁵× the 1e-9 bound the rest of the corpus holds (a healthy case sits at
1e-15). The correlation is exact:

| case | `dt_therm_ratio` | `dt_tracer_advect_ratio` | surface heat flux | Heat `Error` |
|---|---|---|---|---|
| `acc_channel` | 1 | – | −40 W/m² | 1e-15 ✅ |
| `acc_channel_quiescent` | 1 | – | 0 | 1e-15 ✅ |
| `eady_weno5` / `eady_weno7` | 2 | 2 | **none** | 1e-15 ✅ |
| `acc_channel_weno5` | 2 | 2 | −40 | **5e-5** ❌ |
| `acc_channel_kitchensink` | 2 | 2 | −40 | **5.6e-5** ❌ |
| `acc_channel_sw_penetration` | 2 | 2 | +40 | **3.9e-5** ❌ |
| `acc_channel_kitchensink_xl` | 2 | 2 | −40 | **5e-5** ❌ |

So the leak is in the windowed-drain / thermo-cadence path's handling of a
**surface flux** — it needs both the cadence lag and a flux to appear — and it
**survives** the 2026-09-11 concentration-rescale fix, which addressed the
energy instability on the same path. See
`docs/spec_windowed_tracer_advect_hardening.md`.

### 9. `seamount_conservative_floor` — a quiescent seamount that would not stay quiescent — **FIXED**

**Symptom.** Peak `En` 2.24e-3, i.e. **6.7 cm/s of spurious current** on a
motionless ocean: 45× the 1 cm/s seamount-test bar and ~60× its own siblings.
Worse than the magnitude, it **ended at its peak** — the spurious mode never
equilibrated (30 simulated days reached 1.2 cm/s and were still climbing).

**Root cause — the grounded-layer PGF Jacobian, not the coordinate.** The face
pressure gradient is a two-point Jacobian,

```
PGF_x = -(1/rho0) * [ (p_c^R - p_c^L)/dx  +  g*rho_layer*(z_c^R - z_c^L)/dx ]
```

whose two terms cancel at rest only while the two abutting layer centres lie in
a common z-interval whose **ambient** density is that layer's own `rho_layer`.
Under `vcoord_type="lagrangian"` a layer grounded against sloping topography is
squeezed onto the `angstrom_h` floor on the shallow side of a face while
remaining massive on the deep side — the two centres are then *hundreds of
metres* apart in z, the interval between them is filled with other density
classes, and what survives the cancellation is
`g*(rho_layer - rho_ambient)*dz/dx` of acceleration on a resting ocean.

Measured at t = 0 on the shipped namelist: the PGF is **exactly zero**
(1e-16 … 1e-18) on every face whose layer is massive on both sides, and
non-zero on precisely the faces that touch a vanished cell, with the amplitude
falling off linearly in the number of density classes between the layer and the
bed — the analytic signature of the term above. Nothing bounded the resulting
face velocity, which then leaked into the interior through the barotropic
correction (a *uniform* per-layer increment, which is why the surface layer —
whose own PGF is identically zero — ended up moving at 11 cm/s).

**What it was not.** The coordinate itself is innocent, and so is the PGF form:
`fv_mom6` reproduced the same defect (4.16e-4 at day 1, 2.17e-3 at day 2, and
it is gated now too — see below); running the same namelist under
`sigma` / `zstar_sigma` / `zstar_full` was 100× *worse* (that IC is only
meaningful in an isopycnal coordinate); and at `angstrom_h = 1e-10`, where the
floored layers carry no mass at all, the spurious energy was unchanged — so it
was never a thin-layer bookkeeping artefact.

**Fix.** `&ocean_isopycnal_nml pgf_skip_nonoverlap` (**default ON**,
`VCOORD_LAGRANGIAN` only, `fv_lite` / `fv_wright` / `fv_mom6`) zeroes the face
PGF wherever the layer's z-extents in
the two abutting columns do not overlap — exactly where the layer has wedged
out against the bed and there is no common depth to difference the pressure
across. Overlapping layers are untouched, so every other vertical coordinate is
bit-identical and a wedged-out layer can still be re-wetted by continuity's own
upwind flux and by the barotropic correction.

**After:** `En` 2.24e-3 → **1.39e-26** at day 2 and 1.58e-26 at day 30 — machine
zero, held for 30 simulated days. The case now carries the `REST_1UM_S`
"must be bit-zero" bar instead of the 1 cm/s seamount tolerance. Pinned by
`tests/test_ocean_pgf_grounded.F90`.

The same namelist run with `form="fv_mom6"` goes 2.17e-3 → **1.91e-27** at
day 2 on the same gate. `form="mont"` is still NOT gated and still NaNs on this
case at day 1 — that form differences layer-centre pressure with no
z-correction whatsoever, which is structurally wrong on any coordinate whose
layer centres are not z-aligned, so the grounding gate is not what it is
missing. `form="fv_wright"` runs but settles at 6.3e-6 rather than machine
zero: it rebuilds the pressure stack from the *Wright in-situ* density of the
S/T fields (uniform here) instead of the isopycnal `rho_layer` the coordinate
was seeded from, so the seeded state is not an exact discrete rest state for
it. Neither of those is a defect of this gate.

### Near-misses worth watching (currently passing)

- the whole 96×96 `zstar_sigma` seamount family sits at **1.08 mm/s** of
  spurious velocity — 9× under the 1 cm/s bar, but three orders of magnitude
  above `seamount_flat`'s 1e-23. The number is printed on every run so drift
  toward the bar is visible in review.
- `baroclinic_15layer` reaches `En = 0.344` (|u|_rms 0.83 m/s) by day 50,
  against a 2.0 ceiling.
- the tracer-concentration diag `mean` is NaN on every case with land, a halo
  or a vanishing coordinate — a real (if benign) defect in the reduction, not a
  blow-up. See *the NaN that is not a bug*.
