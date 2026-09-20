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
| `manifest.py` | the curated case list — nml path, `n_steps`, `timeout_s`, `backends`, physics `tags`, optional per-case `tol` |
| `run_regression.py` | the P0 orchestrator (CPU-serial \| GPU work-queue), the run-clean/NaN gate, the report |
| `compare.py` | the P1 golden-summary compare (`--update-golden` \| `--compare`) |
| `golden/<case>.json` | committed per-case final-state field summaries (tiny) |
| `coverage.py` | P2 gcov coverage mode — measures ocean-closure coverage of the corpus + emits the gap list |
| `run_all.py` | **single entry point** for the golden suite — invokes compare (CPU/GPU) + coverage, aggregates to one pass/fail + CI-ready exit code |
| `stability.py` | **the two-tier stability suite** (see below) — parses the model's own console time series and asserts on PHYSICS, not on a golden |
| `stability_manifest.py` | its case list: all 68 tracked ocean namelists, with per-case run length, physics assertions, tier-2 downscale spec and known-failure markers |
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
`|golden − run| ≤ atol + rtol·scale`, where `scale` is the *field's* magnitude
(diag field: `max|min|,|max|` over golden **and** run; stats scalar:
`max|golden|,|run|`). This means a **physically-zero field carrying roundoff
noise** (e.g. a quiescent seamount's velocities at ~1e-12, which differ ~6%
*relatively* CPU↔GPU) is absorbed by `atol`, while an active field is bounded by
`rtol·scale` and a spurious spin-up lifts the scale and **fails loudly**.

**Chosen defaults (measured, not guessed)** — `rtol = 1e-3`, `atol = 1e-7`. The
worst REAL-field (non-roundoff) CPU(gfortran)↔GPU(nvfortran) relative drift
across the 22 well-behaved cases is **1.4e-6** (`ideal_age_demo` SSH:min);
`rtol=1e-3` is ~700× that — loose enough to absorb FMA / reduction-order /
transcendental-library divergence (GPU is deterministic run-to-run here) yet
tight enough that a real physics regression fails. Goldens are generated on the
**CPU** build; `--compare --backend gpu` then passes against the same goldens
within these tolerances.

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
