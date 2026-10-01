---
name: pr-review
description: Review a roundabout branch, PR, or working-tree diff against the project's Fortran house style and its GPU / MPI / compiler-portability rules — the conventions the pre-commit hooks and CI cannot check. Use when asked to review a PR, review a diff, check Fortran style, or check whether a change is up to standards before pushing. Also the skill the `/claude-review` GitHub workflow runs.
---

# roundabout PR review

Review changed Fortran against `FORTRAN_STYLE.md` and the rules in `CLAUDE.md`
(its *Gotchas* section is the list of things that have already cost real
debugging time). The point of this skill is the part **no hook catches** — read
what CI and the hooks already report first, so you never spend a finding on it,
then read the diff for the rules below.

`CLAUDE.md` wins where it and `FORTRAN_STYLE.md` disagree. The style guide was
inherited from another project, and these parts of it are **stale** — do not
raise findings from them:

| `FORTRAN_STYLE.md` says | roundabout actually does |
|---|---|
| kinds from `pic_types` (`dp`) | `wp` from `rdb_constants`, `_wp` literals |
| files `rdb_<name>.f90`, tests `test_rdb_<name>.f90` | every source file is `.F90`; tests are `tests/test_<name>.F90` |
| errors via `error_t` / `create_error` from `rdb_error` | no such module: optional `ierr` returning `OCEAN_STATUS_*` codes (`rdb_ocean_status`), errors logged with `global_logger%error`, and `error stop` (fail loud) when no `ierr` is passed |
| units Bohr / Hartree, `to_bohr()` | SI: m, s, kg, Pa, °C, psu |
| use `associate` for long expressions | **never** `associate` over a derived-type component around a `do concurrent` (see §3a) |
| `pic_blas` for BLAS | not used anywhere; there is no BLAS in this code |

## 1. Establish the target

Unless one was named, review the branch against `main`:

```bash
git diff origin/main...HEAD --stat     # scope
git diff origin/main...HEAD            # the change itself
git log --oneline origin/main..HEAD    # one capability per commit?
```

For a GitHub PR: `gh pr view <N>` (title, body, base — a stacked PR's base is
another PR's branch, so diff against that base, not `main`), `gh pr diff <N>`,
`gh pr checks <N>`.

## 2. What is automated — do not review by hand what these report

**pre-commit** (`.pre-commit-config.yaml`; locally
`pre-commit run --from-ref origin/main --to-ref HEAD`):

| hook | owns |
|---|---|
| `fprettify` | layout — never comment on indentation or line breaks |
| `fortitude` (`fortitude.toml`, line length 178) | `implicit none`, `use ... only`, `private`, forbidden constructs. It **excludes `tests/`, `benchmarks/` and `tools/`**, so check those by eye |
| `no-mpi-in-rdb` | direct MPI anywhere (use `pic_mpi_lib`) |
| `dc-assumed-shape` | assumed-shape dummies in `do concurrent` kernels (waiver `! assumed-shape-ok: <reason>`) |
| `dc-intrinsic-shadow`, `dc-transfer` | intrinsic-shadowing locals and `transfer()` in `do concurrent` |
| `decl-order` | integer dims declared before the explicit-shape arrays that use them (ifx #8586) |
| `vanished-layer` | hand-rolled `hTr/h` or `H_VANISHED` tests (waiver `! vanished-ok: <reason>`) |
| `openmp-portability` | new `!$acc` shapes that break ifx's OpenMP path |
| `check-closure-matrix` | knob / `docs/CLOSURE_MATRIX.md` / test drift |
| `test-regime-labels` | every rdb ctest entry carries an `ocean`/`core` label |
| `check-bytes-accounting` | every array counted in its type's `bytes()` |
| `commit-msg-no-session-links` | Claude session links in commit messages |

The hooks that are diff-aware against `origin/main` only see *new*
offenders. Judge any **new waiver** comment in the diff by hand: the reason
must be specific and true.

**CI** (`.github/workflows/`): `portability.yml` builds and runs `ctest` on
gcc 15, intel 2025, ifx + Intel MPI, nvidia-hpc, nvfortran + HPC-X, gfortran
+ OpenMPI, plus **flang-new and LFortran with NetCDF off**; `ocean-stability.yml`
runs the gfortran tier-2 stability sweep; `pre-commit.yml` runs the hooks.
Read them with `gh pr checks <N>` and quote failures verbatim.

Report hook and CI hits as **"CI will fail on this"**, separate from your own
findings. When running inside the GitHub workflow, the Fortran tools are not
installed — read `gh pr checks` instead of running them, and do not build.

## 3. Review by hand — what escapes the hooks

Ordered roughly by how much damage each does. Cite `file:line` for every
finding and say what to do.

### 3a. GPU correctness — the GPU build is `mem:separate`, no unified memory

A kernel gets **no** implicit host↔device copies. A test that only runs on
the CPU build passes while the GPU build is wrong — so these are Must fix even
when every CI leg is green.

- **Every array a kernel touches must be device-present.** A new array on a
  state type needs its `!$acc enter data` (and `exit data`) next to its
  siblings'. In tests, map every state object **and** its scratch companion.
- Arrays mapped `create` carry no host values: host-set inputs need
  `!$acc update device(...)` after `enter_data`; host reads after a kernel
  need `!$acc update self(...)`.
- **Setup-time host editors run before the map**, not after (the
  `metrics_apply_land_mask` trap: the device stays all-wet while every
  host-side check passes).
- **Never `!$acc update self` / `copyin` a whole derived type with
  allocatable components** — update the component arrays.
- **Never instantiate a default-initialised derived type in device code**
  (`!$acc routine seq` result, `intent(out)` dummy or local): nvfortran's
  `fort2` segfaults with no diagnostic. Pass scalars; `intent(in)` bundles
  are fine.
- **No `associate` over a derived-type component around a `do concurrent`**:
  ifx reads it as zero inside the loop. Spell the component out or pass an
  explicit-shape `_impl` dummy.
- **No local `allocate` in a per-step kernel**: use a module-level workspace
  with `*_workspace_ensure` and a `*_cleanup` in the exit-data path.
- **A host-gated call that passes a state array to another subroutine still
  costs** (+4.8 % measured, with the knob off): write the guarded pass inline
  as a `do concurrent` in the same routine.
- An `ieee_get_flag` assertion is vacuous on the device build; it must be
  skipped there, never left to "pass".

### 3b. Numerics that fail silently

- **Clamps launder NaN under `-fast`**: an `if/else` clamp becomes a
  NaN-blind min/max. Any clamp that can see non-finite data needs an
  `ieee_is_finite` guard; `if (abs(u) > thresh)` silently skips NaN too.
- **Fixed-point (EFP) reductions must propagate non-finite input** through the
  `poison` counter — never convert a NaN to an integer bin.
- **The vertical stack is bottom-up**: `k = 1` is the bed, `k = nz` the
  surface. Surface forcing lands at `nz`, bed forcing at `1`. A new layered
  kernel that assumes the opposite is Must fix.
- **Thin layers**: dynamic vanishing tests use `H_VANISHED`, pure division
  armour uses `H_DIV_EPS` (both `rdb_constants`). Vanished-layer content goes
  through the `rdb_vl_*` helpers in
  `src/shared_module_utilities/rdb_vanished_layer.inc`, never a re-derived
  rule; the enforcement point is `multilayer_state_t%enforce_vanished_content`.
- **Formula bathymetry** must fill ghost rows, and its length scales are in
  grid units (degrees on spherical grids) — convert via
  `topo_length_to_grid_units`.

### 3c. MPI — a decomposed run must be bitwise the serial run

- A new pass that writes a tile's physical span only (an OBC fill, a sponge,
  a wall re-close), or keys on a global edge tag, must gate on `has_*` and
  refresh seam ghosts afterwards if anything reads them before the next
  exchange. Keying a closure on the edge tag alone turns an MPI seam into a
  wall.
- A new feature either joins `tests/mpi/test_ocean_decomp_bitid_mpi` or fails
  loud in `engine_setup` on more than one rank. A feature that silently runs
  multi-rank without either is Must fix.
- All MPI goes through `pic_mpi_lib` (the hook checks the imports; check new
  collectives are reached by every rank).

### 3d. Compiler portability (gcc, ifx, nvfortran, flang-new, LFortran)

- **NetCDF-only code** — sources, test cases and their `tester.F90` /
  `tests/CMakeLists.txt` registration — sits under `#ifndef RDB_NO_NETCDF`,
  and a whole NetCDF-only test goes in `NETCDF_REQUIRED_TESTS`. The flang-new
  and LFortran legs build with NetCDF off and fail otherwise.
- **A new importer of `NZ_STACK_MAX`** needs the module-local
  `#ifdef LFORTRAN_PASSING` parameter like its siblings, and the file count
  in `.github/workflows/portability.yml` must follow.
- Never assign a scalar to a whole deferred-length `character` array
  (`lines = ""` re-allocates with `len=0` on ifx); use `lines(:) = ""`.
- Intrinsic modules need `use, intrinsic ::`.

### 3e. Knobs — the plumbing that fails silently

A new namelist knob (`docs/howto/add_namelist_knob.md`) needs **all** of:
the field and its default in `src/core/rdb_config.F90`; registration in that
group's `register_ocean_<group>`; validation if any combination is illegal
(fail loud); propagation in `src/core/ocean/state/rdb_ocean_setup.F90`; the
kernel gated on it; a `!!` docstring on the field;
`docs/generated_nml_knobs.md` regenerated; and the Python ledger —
`python/rdb/_config_generated.py` (descriptors and `N_KNOBS`) plus the count
in `python/tests/test_config.py` — **regenerated with
`tools/gen_python_config.py`, never text-merged**. After a rebase, check the
three counts still agree.

- **Default off ⇒ bit-identical** for every existing namelist and test (KPP
  and ALE remap are the on-by-default exceptions). A new default-on knob
  changes answers — see 3f.
- Check the knob does not already exist under another name.

### 3f. When answers move

- A changed `tests/regression/golden/*.json` needs the reason in the PR body,
  the maintainer's authorisation, and its own commit. Unexplained golden
  churn is Must fix.
- The same goes for the stability pins (`tests/regression/stability_manifest.py`,
  `vcoord_matrix_measured.py`) and the global reference year
  (`validation_examples/ocean/global_1deg/reference_daily.csv` plus its README
  table). A re-pin that widens an envelope must say why.
- **A performance claim needs a repeated measurement**, ideally NVTX ranges
  per kernel; single-shot before/after numbers are a Should fix.
- A validation README that quotes results should name the build (commit or
  PR stack) that produced them, so they can be reproduced.

### 3g. Interface quality and style

- **`intent` on every dummy argument**, no exceptions.
- **Six arguments or fewer on a public procedure**; group the rest in a
  derived type. Private kernels and `_impl` helpers with explicit-shape
  arrays get latitude.
- **`private` by default** with an explicit `public ::` list; a new public
  entity needs a reason to be public. `use ..., only:` everywhere, without
  pulling in far more than the file uses.
- **`pure` by default.** A new non-pure procedure needs a real side effect
  (I/O, logging, MPI, module state, a lazy allocate).
- **No naked output** (`print *`, `write(*,*)`, `write(6,*)`): use
  `pic_logger`'s `global_logger` with a real level. Watch for leftover debug
  writes.
- **No literal kinds** (`real(8)`, `1.0d0`): `wp` and `_wp`.
- **No magic numbers**: a named `parameter` for thresholds and tolerances.
  Physical constants live in `rdb_constants` (`GRAVITY`, `RHO_WATER`, ...);
  a re-declared or truncated copy is a finding — the `g_bt = 9.81` vs
  `GRAVITY = 9.80665` mismatch once produced a spurious barotropic force.
- **`destroy`** on any new type with allocatable components; `allocatable`
  over `pointer`; `character(len=:), allocatable` over fixed lengths; `block`
  to scope temporaries; nesting deeper than three wants an early return or a
  helper.
- **A helper body shared across kernel modules** goes in
  `src/shared_module_utilities/<concern>.inc` and is `#include`d — not a
  hand-duplicated `*_impl` copy, and not a second inline copy of a formula an
  include already holds.
- **No emojis** in `.F90` files, comments and strings included.
- **Cite the paper, not another codebase** (`!! Wright (1997)`, not
  `!! ported from MOM_EOS_Wright.F90`). GPL code is off-limits — flag anything
  that reads as ported from a GPL model.

### 3h. Naming

| thing | rule |
|---|---|
| source file | `rdb_<name>.F90`, module of the same name, one per file |
| test | `tests/test_<name>.F90`; MPI tests in `tests/mpi/` |
| derived type | `snake_case` with a `_t` suffix |
| variables, procedures | `snake_case`, descriptive; single letters only for loop indices |
| `parameter` | `UPPER_SNAKE_CASE` |

### 3i. Documentation

- `!!` FORD docstrings, **below** the declaration they document, on every new
  public type, procedure, type component and knob field. Say what it is and
  its units, not what the code obviously does. Module/procedure docs come
  from these docstrings — a hand-written markdown copy is a finding.
- The **synthesis docs** must move in the same PR when a closure, scheme,
  knob or extension seam is added, removed or renamed:
  `docs/CLOSURE_MATRIX.md`, `docs/CAPABILITIES_AND_LIMITATIONS.md`,
  `src/core/ocean/README.md`, the relevant `docs/howto/` guide, and
  `CLAUDE.md` when the change creates a new gotcha. A doc that now disagrees
  with the code is a bug.

## 4. Coverage

- Does a new module, public procedure or bug fix come with a test? A
  bug fix's test must fail without the fix — say so if the PR does not show
  it. Analytical tests are preferred: every one added to the ocean core has
  caught a real bug.
- Is the test registered in **both** `tests/CMakeLists.txt` (`RDB_TESTS`,
  with its regime label) **and** `tests/tester.F90` (the fpm aggregate)? MPI
  tests go under `RDB_ENABLE_MPI` with their rank counts.
- Does a GPU-touching test map everything it hands to a kernel (3a)?
- Does the PR say which `ctest` scope it ran (`ctest -R rdb` vs the full
  suite) and on which toolchains? A GPU-path change needs a GPU run.

## 5. Report

Group as:

1. **CI will fail** — verbatim hook or check output, with `file:line`.
2. **Must fix** — correctness, GPU data motion, MPI bit-identity,
   portability, or a house rule broken.
3. **Should fix** — style and interface-quality findings.
4. **Consider** — suggestions the author may reasonably decline.

Every finding cites `file.F90:line` and says what to do, not just what is
wrong. If a section found nothing, say so in one line rather than omitting it,
so the author knows it was looked at. Do not invent findings to fill a
section. Do not approve or request changes; the maintainer decides.

**In the GitHub workflow:** post the whole report as one PR comment with
`gh pr comment <N> --body-file <file>`, and put each **Must fix** finding
inline on its line with the inline-comment tool as well. Post nothing else.
