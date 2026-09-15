# Extending Roundabout — contributor how-to guides

Recipe-style walkthroughs for the common extension points, written so they work
whether you're editing by hand or driving an LLM coding assistant (Claude Code,
etc.). Each guide traces one *seam*: the dispatch mechanism, the ordered list of
files to touch, the test to copy, and the GPU gotchas that bite.

These are **process docs**. They deliberately do *not* restate what each module
or procedure does — that lives in the source `!!` docstrings and is rendered by
FORD on every PR to `main`. When a guide says "the kernel computes X", go read
the kernel; the guide's job is to tell you *where the seam is and what to edit*.

> Symbol names (modules, procedures, enums, namelist groups) are stable and are
> what these guides anchor on. Line numbers drift — treat any you see in git
> history or agent output as "approximately here", and grep for the symbol.

## Read these first (orientation, not duplicated here)

| You want… | Read |
|---|---|
| Repo map + run lifecycle | [`docs/codebase/INDEX.md`](../codebase/INDEX.md) |
| Acronyms / concepts | [`docs/codebase/CONCEPTS.md`](../codebase/CONCEPTS.md) |
| Which closure/scheme is on + the knobs | [`docs/CLOSURE_MATRIX.md`](../CLOSURE_MATRIX.md) |
| Ocean design contract (god-state slot map, the 4 rules) | [`src/core/ocean/README.md`](../../src/core/ocean/README.md) |
| Every namelist group + knob + default (code-derived) | [`docs/generated_nml_knobs.md`](../generated_nml_knobs.md) |
| Physics + namelist reference | [`docs/REFERENCE.md`](../REFERENCE.md) |
| Style + GPU programming (`do concurrent` + OpenACC) | [`FORTRAN_STYLE.md`](../../FORTRAN_STYLE.md) |
| Capabilities + limits | [`docs/CAPABILITIES_AND_LIMITATIONS.md`](../CAPABILITIES_AND_LIMITATIONS.md) |

## The guides

| Guide | Has a dispatch seam already? |
|---|---|
| [Add a namelist knob](add_namelist_knob.md) | n/a — the config plumbing itself |
| [Migrate a namelist group to the schema](migrate_namelist_to_schema.md) | n/a — the config plumbing itself |
| [Add a passive tracer](add_passive_tracer.md) | yes — the tracer registry "just works" |
| [Add a vertical coordinate](add_vertical_coordinate.md) | yes — `VCOORD_*` enum + target-thickness generator |
| [Add an equation of state](add_equation_of_state.md) | yes — `EOS_VARIANT_*` + `&ocean_eos_nml eos` |
| [Add a diagnostic](add_diagnostic.md) | yes — the diag registry + `fill_*` |
| [Add a turbulence closure](add_closure.md) | yes — the vmix contributor chain / `LMIX_*` select-case |
| [Add a boundary condition](add_boundary_condition.md) | yes — the `OBC_*` enum + per-edge dispatch |

## Two-minute orientation

- **One regime.** `sim_type` is pinned to `'ocean'` (the only enum value; a
  different value is a fail-loud `validate_config` error). The A-grid coastal
  path and the unstructured triangular backend used to live here and have been
  split into their own repository — if you find a doc, comment or branch that
  says "mirror this on the unstructured side", it predates that split.
- **The dyn-core** is an Arakawa C-grid, continuity-PPM, split-explicit RK2
  hydrostatic ocean model under `src/core/ocean/`, with shared physics in
  `src/{equation_of_state,pressure_force,tracer,ALE,parameterizations}/` and an
  optional sea-ice subsystem under `src/core/ice/`.
- **State composition**: `ocean_state_t`
  (`src/core/ocean/state/rdb_ocean_state.F90`) is a thin composition of slot
  types, each owning its allocatables + bound `init`/`destroy`/`enter_data`.
  Read [`src/core/ocean/README.md`](../../src/core/ocean/README.md) — it is the
  slot map, the 4 rules, and the "how to pick up a slot" recipe.
- **Test regimes**: every row in `tests/CMakeLists.txt` carries a `<regime>`
  label — `ocean` (C-grid dyn-core + sea-ice), `core` (shared infrastructure),
  or `ocean+core` (both). Those three are the *only* legal values; configure
  fails loud on anything else, and the `test-regime-labels` pre-commit hook
  enforces it too. Run one with `ctest -L ocean` / `ctest -L core`.

## Anatomy of a capability (the standard workflow)

Every guide is a specialisation of this. From the top-level `CLAUDE.md`:

1. **Namelist knob, default = off.** Preserves bit-identity for existing namelists
   and tests. This is non-negotiable — every existing run must stay byte-for-byte.
2. **Kernel implementation, gated on the knob.**
3. **Unit test covering the new path.** Analytical tests are the highest-leverage
   thing you can write here — *every* analytical test added to the ocean core has
   caught a real bug. Default to "add a test" over "audit the code".
4. **`pre-commit run --all` + fortitude** — pass before committing.
5. **`ctest` green on the GPU build** (never `-j N` — every worker shares one GPU;
   `ctest -R rdb` is the rdb-only regression subset, ~177 tests).
6. **One capability per commit, one PR per capability.** Smaller PRs review faster
   and bisect cleanly. (Exception: several branches that each just add a namelist +
   IC + test → bundle into one `feat/<area>` PR.)

## Cross-cutting GPU rules (you *will* hit these)

The solver is GPU-native via `do concurrent` + OpenACC. State arrays stay
device-resident; only diagnostic output triggers D→H copies. These rules recur in
every guide — they're collected here so the guides can just point at them. Full
detail in [`FORTRAN_STYLE.md`](../../FORTRAN_STYLE.md) and the "Gotchas" section of
[`CLAUDE.md`](../../CLAUDE.md).

- **`do concurrent` loop order is `(k, j, i)` / `(j, i)` — j outer, i inner.**
  The reverse kills NVHPC GPU coalescing.
- **Flat-impl kernels take explicit-shape arrays** `arr(nx, ny, nz)`, never
  assumed-shape `arr(:,:,:)` — assumed-shape makes NVHPC emit a descriptor-walk
  memcpy per launch (millions of them in a real run). Enforced diff-aware by the
  `dc-assumed-shape` pre-commit hook.
- **Never `allocate` scratch inside a per-step kernel on `-stdpar=gpu`.** It stalls
  the device on entry/exit and bridges device inputs with implicit memcpys. Two
  sanctioned patterns: module-level allocatable + lazy `*_workspace_ensure` +
  `*_cleanup`, or a field on the slot mapped once in `enter_data`. Per-column
  scratch uses fixed-size `NZ_STACK_MAX` stack arrays (a `local()` clause requires
  fixed size; a dummy-sized automatic array crashes with
  `CUDA_ERROR_ILLEGAL_ADDRESS`).
- **A loop-invariant `select case` goes *inside* one `do concurrent`, not split
  into one loop per case.** Splitting roughly doubles GPU launches (~1.8× slower);
  a uniform branch inside one loop is ~free.
- **A new allocatable on an ocean slot must be wired into the
  `ocean_state_enter_data` orchestrator.** A silent omission is a 150–1500× memcpy
  explosion, not a crash — so it's easy to miss.
- **The GPU build is `mem:separate`** (no managed/unified memory): a kernel gets
  **no** implicit host↔device copies, so every array it touches must be
  device-present or it silently reads/writes stale host memory. A test that only
  runs on the multicore/host build proves nothing about device data motion. The
  canonical template is `test_open_boundary_out_closes` in
  `tests/test_ocean_conservation_salt_heat.F90`.
- **Never `!$acc update self` / `copyin` a WHOLE derived type with allocatable
  components** — the aggregate D→H copy overwrites the host descriptors with
  device addresses. Update the component arrays.
- **`pure` is the default** for new procedures. A `pure` helper called from inside
  a `do concurrent` GPU loop needs `!$acc routine seq`. NVHPC won't inline such
  helpers across module boundaries in hot kernels — duplicate the body as an
  `_impl` copy only when a profile justifies it (don't preemptively).
- **Bottom-up layer convention: `k=1` is the bed, `k=nz` is the surface.** Surface
  forcings land at `(:,:,nz)`; bed forcings at `(:,:,1)`. Any new per-layer kernel
  must follow it. Regression gates that guard it: `test_ocean_surface_flux`
  (positive Q warms `k=nz`), `test_ocean_geothermal` (bed heat at `k=1`),
  `test_ocean_conservation_salt_heat`.
- **Default-off = bit-identity** (rule 1 above, restated because it's the one most
  often forgotten): the default value must reproduce existing behaviour
  (`0.0` / `.false.` / a sentinel / a `1.0` multiplier), and the kernel must be
  gated so an un-set knob changes nothing.

## Working with an LLM coding assistant here

This repo is set up to be navigated by an agent. A few things make that work well:

- **The top-level [`CLAUDE.md`](../../CLAUDE.md) is the agent's primary brief** — it
  carries the architecture, the layer convention, the GPU rules, and a "Gotchas"
  section that is essentially a list of mistakes an LLM would otherwise make. Point
  the assistant at it first.
- **Give the assistant the relevant guide from this directory** as the task seam,
  plus the design contract README (`src/core/ocean/README.md`).
- **The assistant cannot see GPU correctness from the source.** NVHPC stdpar +
  OpenACC has sharp edges (descriptor-walk memcpys, illegal-address crashes from
  automatic arrays, silent approximate-math from intrinsic-shadowing locals). The
  ground truth is: build on the NVHPC GPU target and run `ctest` (no `-j`). Treat a
  green CPU/gfortran build as necessary-not-sufficient.
- **Prototype algorithm questions in Python first.** For a per-layer / per-cell
  numerical question, a tiny Python prototype on a representative setup answers it
  in seconds versus a ~9-minute GPU rebuild. Keep such Python prototypes and
  minimal `.F90` reproducers local — they're a development aid, not shipped
  repo artefacts.
- **Verify, then claim.** If `ctest` fails, say so with the output; don't report a
  refactor as done until it's green on the GPU build. The analytical tests are the
  arbiter.
