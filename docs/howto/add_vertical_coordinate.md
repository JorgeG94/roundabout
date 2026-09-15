# How to add a new vertical coordinate

Adds a new `VCOORD_*` family. There is one production implementation — the ocean
C-grid coordinate type `ocean_vcoord_t` in
`src/core/ocean/vcoord/rdb_ocean_vcoord.F90` — plus two shared pieces it builds
on:

- **The enum**, in `src/core/rdb_constants.F90` (`VCOORD_LAGRANGIAN`,
  `VCOORD_EULERIAN_Z`, `VCOORD_SIGMA`, `VCOORD_ZSIGMA`, `VCOORD_ZSTAR`,
  `VCOORD_ZSTAR_FULL`, `VCOORD_ZSTAR_SIGMA`, `VCOORD_Z_FIXED`, `VCOORD_RHO`,
  `VCOORD_HYCOM`).
- **The string parser** `parse_vcoord_type` in `src/ALE/rdb_vcoord.F90`, wrapped
  by `parse_ocean_vcoord_type` (same parser, default `VCOORD_EULERIAN_Z`).

> `rdb_vcoord.F90` also still carries `vcoord_t` and the column-level
> `vcoord_target_dz_column` / `zstar_full_build_column` generators. Only
> `parse_vcoord_type`, `parse_remap_method` and `STRETCH_*` have production
> consumers today; `vcoord_target_dz_column` survives as a column-level
> reference implementation exercised by `tests/test_vcoord_target.F90`. Adding
> an arm there is optional — mirroring a new coordinate into it buys you a cheap
> column-level unit test, nothing more.

## The seam — key insight

The remap *reconstruction* (`remap_column` in `src/ALE/rdb_remap_column.F90`)
does **not** switch on the coordinate — it switches only on the reconstruction
method (PCM/PLM/PPM/PPM_H4/PQM). A coordinate is defined *entirely* by its
**target-thickness generator**: the routine that, given total depth `H` and `η`,
produces the target layer thicknesses the conservative remap maps onto. That's
where the `select case (coord_type)` lives:

- `ocean_vcoord_compute_target_h(this, total_h, eta)` →
  `ocean_vcoord_compute_target_h_impl` in `rdb_ocean_vcoord.F90` — fills
  `target_h(i,j,k)` directly in 2D/3D `do concurrent` loops (everything inlined;
  no cross-module call). Its `case default` is `error stop`, so an unhandled
  coordinate crashes rather than silently falling back.
- Density-space coordinates (`VCOORD_RHO`, `VCOORD_HYCOM`) take a separate
  entry, `ocean_vcoord_compute_target_h_rho(..., T, S, eos, hybrid)`, because
  they need the EOS and an interface inversion (`invert_density_targets`).
- The per-column-reference coordinate (`VCOORD_ZSTAR_FULL`) is special-cased
  with its own builder, `ocean_vcoord_build_zref_full`, run once at setup from
  local bathymetry.

## Recipe (ordered) — add `VCOORD_FOO`

1. **Enum** — add `integer, parameter :: VCOORD_FOO = <n>` in
   `src/core/rdb_constants.F90` and to its `public ::` list.
2. **Parser** — add a `case ("foo", "FOO", …)` → `VCOORD_FOO` in
   `parse_vcoord_type` (`src/ALE/rdb_vcoord.F90`). Re-export `VCOORD_FOO` from
   `rdb_ocean_vcoord.F90` if callers outside need the tag.
3. **`ocean_vcoord_t`** (if FOO needs new params): add the knob fields; allocate
   any new array in `ocean_vcoord_init` and release it in
   `ocean_vcoord_destroy`; map it in `ocean_vcoord_enter_data_impl` /
   `_exit_data_impl`. **A new ocean array must also be reached from the
   `ocean_state_enter_data` orchestrator** or you get a 150–1500× memcpy
   explosion.
4. **Target generator** — add `case (VCOORD_FOO)` to
   `ocean_vcoord_compute_target_h_impl`, filling `target_h`. Bind every
   derived-type component you touch through the existing `associate` block (ifx
   ICEs on a `this%component` reference inside an offloaded `do concurrent`).
   If FOO is density-space, extend `ocean_vcoord_compute_target_h_rho_impl`
   instead and add the dispatch in `rdb_ocean_remap.F90` that routes to it.
5. **Remap early-return** — if FOO is a no-op coordinate (like
   `VCOORD_LAGRANGIAN` / `VCOORD_EULERIAN_Z`), add it to the early-return guards
   in `src/ALE/rdb_ocean_remap.F90` (there is one per remap entry point).
6. **Config + namelist** — `vcoord_type` already exists in `&vcoord_nml`; add
   the new string to its enum registration (`register_vcoord` in
   `src/core/rdb_config.F90`) and any new knob fields per
   [add_namelist_knob.md](add_namelist_knob.md).
7. **Driver wiring** — `src/driver/rdb_driver.F90` already routes
   `ocean_state%vcoord%coord_type = parse_ocean_vcoord_type(cfg%vcoord_type)`;
   copy any new knob from `cfg%` onto `ocean_state%vcoord%…` in the same block.
   It runs **before** the IC seed, because the seed calls
   `vcoord%build_zref_full(b)` — keep that order if FOO needs a per-column
   build at setup.
8. **Initial thickness** — `&vcoord_nml thickness_config` seeds `h_layer` at
   setup independently of the coordinate. If FOO needs a matching seed profile,
   add it there too, or a run starts far from its own target and spends the
   first steps remapping.
9. **Docs + tests** — add a row to the "Vertical coordinates" table in
   `docs/CLOSURE_MATRIX.md` and the tests below.

## Tests — the pattern to copy

- `tests/test_vcoord_parse.F90` — string/alias → enum, including the
  default-fallback cases and the ocean wrapper's `VCOORD_EULERIAN_Z` default.
  **Add your aliases here.**
- `tests/test_ocean_vcoord.F90` — the canonical unit test for a new coordinate:
  per-case target-thickness asserts plus a `*_conservation_per_column` case for
  every family (`zsigma_conservation_per_column`,
  `zstar_sigma_conservation_per_column`, `zstar_full_conservation_per_column`),
  the `build_zref` cases, and the lagrangian no-op.
- `tests/test_ocean_vcoord_rho.F90`, `tests/test_ocean_vcoord_hycom.F90` — the
  density-space path + interface inversion.
- `tests/test_ocean_remap.F90` — remap-level: `eulerian_z_is_noop`, sigma
  redistribution, tracer conservation, `zstar_full_tracer_conservation`,
  per-face momentum conservation, budget telescoping.
- `tests/test_ocean_remap_e2e.F90` — solver-level: `e2e_lake_at_rest_sigma` /
  `e2e_lake_at_rest_zstar_full` (a resting stratified lake must stay at rest —
  the sharpest detector of a generator that doesn't conserve), mass/tracer
  conservation, remap cadence.
- `tests/test_vcoord_target.F90`, `tests/test_vcoord_zstar_full.F90` — the
  column-level reference generators.

Invariants to assert: (a) `sum_k target_h == H + η` per column; (b) layers stay
≥ floor; (c) at solver level, tracer/mass conservation across steps, and a
lake-at-rest run stays at rest. **A deep-shallow bathymetry + flow test is
required** to catch over-allocation bugs — flat-bottom / static-stratification
tests are blind to them.

## Gotchas

- **Conservation is load-bearing.** Every branch must guarantee
  `sum_k target_h == H + η`; the conservative remap depends on it. Decide where
  any clip deficit goes (ZSIGMA clips the deep-branch intervals to the local
  column total; ZSTAR_FULL trims the surface).
- **Bottom-up ROMS ordering** (`k=1` bed, `k=nz` surface). Builders that
  construct z-levels top-down must *reverse* into that order — don't flip.
- **Vanishing layers.** ZSTAR_FULL and Z_FIXED clip bed-side layers to
  `zstar_h_min`; operators that divide by `h_layer` gate on `H_VANISHED`
  (dynamic vanish — skip/merge) or `H_DIV_EPS` (pure 1/0 armour), which have
  documented, distinct roles in `rdb_constants`. Pick the right one. A negative
  interface NaNs the CFL — defend degenerate columns.
- **`associate` every component used in an offloaded loop.** ifx's
  `do concurrent` → OpenMP-target lowering ICEs on a derived-type allocatable
  component referenced directly inside the loop body; the whole generator is
  wrapped in one `associate` for that reason.
- **`case default` is `error stop`.** Adding the enum value without adding the
  arm turns FOO into a crash, which is the intended behaviour — but it means
  step 4 is not optional.
- **A no-op coordinate needs the early return too.** Skipping step 5 makes
  the remap run against a `target_h` nobody filled.
- **vcoord choice is a first-order validation lever.** The MOM6 double-gyre gap
  was ultimately sigma-vs-isopycnal — validate FOO against a reference run that
  uses the matching coordinate.

## Source pointers

- `src/core/rdb_constants.F90` — the `VCOORD_*` enum (+ `H_VANISHED` /
  `H_DIV_EPS`).
- `src/core/ocean/vcoord/rdb_ocean_vcoord.F90` — `ocean_vcoord_t`,
  `parse_ocean_vcoord_type`, `ocean_vcoord_compute_target_h[_rho]`,
  `ocean_vcoord_build_zref_full`, `invert_density_targets`.
- `src/ALE/rdb_vcoord.F90` — `parse_vcoord_type`, `parse_remap_method`,
  `STRETCH_*`, and the column-level reference generators.
- `src/ALE/rdb_ocean_remap.F90` — the remap driver + the no-op early returns;
  `src/ALE/rdb_remap_column.F90` — the reconstruction.
- `src/driver/rdb_driver.F90` — the `parse_ocean_vcoord_type` wiring (before the
  IC seed).
- The "Vertical Coordinates" section of the top-level [`CLAUDE.md`](../../CLAUDE.md).
