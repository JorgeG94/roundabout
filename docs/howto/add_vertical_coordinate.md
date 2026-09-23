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
  `ocean_vcoord_compute_target_h_impl` in `rdb_ocean_vcoord.F90` — a host
  dispatcher that hands the `ocean_vcoord_t` components, as explicit-shape
  and scalar dummies, to the flat `ocean_vcoord_geometric_target` kernel
  (`VCOORD_Z_FIXED` to `ocean_vcoord_z_fixed_target`), whose 2D/3D
  `do concurrent` loops fill `target_h(i,j,k)` (no cross-module call). The
  kernel's `case default` is `error stop`, so an unhandled coordinate crashes
  rather than silently falling back.
- Density-space coordinates (`VCOORD_RHO`, `VCOORD_HYCOM`) take a separate
  entry, `ocean_vcoord_compute_target_h_rho(..., T, S, eos, hybrid)`, because
  they need the EOS and an interface inversion (`invert_density_targets`);
  its kernel is the flat `ocean_vcoord_rho_target`, one same-module
  `!$acc routine seq` `ocean_vcoord_rho_target_column` call per column.
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
   `ocean_vcoord_geometric_target`, filling `target_h`. Pass every
   `ocean_vcoord_t` component you need into that kernel as an explicit-shape
   array or scalar dummy (scalars `intent(in), value`, like the existing
   knobs) from the `ocean_vcoord_compute_target_h_impl`
   dispatcher — never reference `this%component` inside the `do concurrent`
   (ifx ICEs) and never wrap it in `associate` over the components (ifx reads
   the names as zero; nvfortran `-stdpar=gpu` handed a by-reference scalar
   associate-name to a device callee as a HOST address — the RHO × Wright
   `CUDA_ERROR_ILLEGAL_ADDRESS`, gated by `test_ocean_vcoord_wright_device`).
   If FOO is density-space, extend `ocean_vcoord_rho_target_column` (reached from
   `ocean_vcoord_compute_target_h_rho_impl`) instead and add the dispatch in
   `rdb_ocean_remap.F90` that routes to it.
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
- **`tests/test_ocean_vcoord_interface_depths.F90` — the interface-DEPTH gate,
  and it is now part of the recipe: a new family MUST add its rows here.**

Invariants to assert: (a) `sum_k target_h == H + η` per column; (b) layers stay
≥ floor; (c) at solver level, tracer/mass conservation across steps, and a
lake-at-rest run stays at rest. **A deep-shallow bathymetry + flow test is
required** to catch over-allocation bugs — flat-bottom / static-stratification
tests are blind to them.

**(d) The absolute geopotential DEPTH of every target interface.** (a)–(c) are
all satisfied *identically* by a stack laid in the wrong half of the column, so
on their own they cannot see a coordinate that puts its fine resolution 500 m
too deep, nor one whose whole column has collapsed into the bed layer — which
is exactly how `VCOORD_ZSIGMA` ran for months against a dimensionless
`z_ref_global` while every sum test stayed green. Assert

    e(K) = −b + Σ_{k ≤ K} target_h(k),   e(0) = −b (the bed),

against a hand-derived analytic table, in three geometries: flat bed at
`η = 0`; a **sloping** bed (so the vanishing branch runs); and the same columns
with the **column top displaced** to `z = −z_top` — a rigid lid, i.e. an ice
shelf. The third is the one that separates the families: a sigma-like family
must divide the live column proportionally, while a z-like family must vanish
the layers whose nominal range lies above `−z_top` against the **top** and keep
the live layers at their open-ocean depths. No family does the latter today
(the builder is handed a column *thickness* and nothing else, so it cannot know
where the column starts); the z-like case-3 rows are therefore carried as
`documents_*` cases that assert the CURRENT placement and are flagged in-line
as the assertions the rigid-top slice must flip.

## Gotchas

- **Conservation is load-bearing.** Every branch must guarantee
  `sum_k target_h == H + η`; the conservative remap depends on it. Decide where
  any clip deficit goes (ZSIGMA clips the deep-branch intervals to the local
  column total; ZSTAR_FULL trims the surface).
- **Bottom-up ROMS ordering** (`k=1` bed, `k=nz` surface). Builders that
  construct z-levels top-down must *reverse* into that order — don't flip.
- **The column top is NOT always `z = 0`.** A family that measures absolute
  depth must take it from `ocean_vcoord_t%z_top(i,j)` — the geopotential depth
  of the top of the WATER column, `metrics%z_draft` under an ice-shelf cavity
  and `0` everywhere else, filled once at configure by
  `configure_ocean_cavity`. It is allocated UNCONDITIONALLY at
  `(nx_total, ny_total)` with `source = 0.0_wp`, precisely so a kernel can take
  it as an explicit-shape dummy without ever meeting `metrics%z_draft`'s
  `(1,1)` placeholder. `z_top ≡ 0` must reproduce the pre-cavity arithmetic
  BIT-for-bit — write the branch so that is true by construction (subtracting
  an exact zero is bit-identical in IEEE round-to-nearest, and stays so under
  FMA contraction), and pin it with an `==` assertion. `VCOORD_Z_FIXED` is the
  worked example: `ocean_vcoord_z_fixed_target`, gated by
  `test_ocean_vcoord_interface_depths` (the per-family interface-DEPTH
  table, which asserts `e(K) = -max((nz-K)*h_nominal, z_top)` under a lid)
  and `test_ocean_vcoord_zfixed_cavity` (the partial cell, the sliver
  merge and the `==` bit-identity). Until your family does this,
  leave it in `validate_config`'s cavity refusal list with a reason.
- **Vanishing layers, at BOTH ends.** ZSTAR_FULL and Z_FIXED clip bed-side
  layers to `zstar_h_min`; under a rigid top Z_FIXED also clips TOP-side
  layers to the same filler and cuts the first live layer into a partial cell.
  A partial cell needs a minimum thickness or it ships slivers: the bed's
  threshold is `zstar_h_min` and the top's is `0.1*h_nominal`
  (`Z_FIXED_TOP_PARTIAL_FRAC`, MITgcm's `hFacMin` / Losch 2008 §2.1), with a
  sub-threshold cut MERGED into the neighbour away from the boundary. Operators
  that divide by `h_layer` gate on `H_VANISHED`
  (dynamic vanish — skip/merge) or `H_DIV_EPS` (pure 1/0 armour), which have
  documented, distinct roles in `rdb_constants`. Pick the right one. A negative
  interface NaNs the CFL — defend degenerate columns.
- **Decide which side of `H_VANISHED` your floor lands on, and say so.**
  The two shipped answers are both correct and they are opposites, so
  `rdb_vcoord :: vcoord_h_min_role` names them and `vcoord_h_min_is_coherent`
  (consumed by `validate_config`) checks the choice. A geometric family (ZSTAR_FULL / Z_FIXED) produces filler layers
  that lie *below the bed*: they hold no water, the floor is only there so
  `target_h` is never exactly zero, and it must stay **≤ `H_VANISHED`** so
  every h-dividing kernel skips them. Thinner is better there — each filler
  interface carries the full topographic slope, so the spurious rest-state PGF
  transport scales *with* the floor. A density family (RHO / HYCOM) collapses
  *real* layers that carry tracer mass anywhere in the column: those must
  survive the remap drain (`h_old > H_VANISHED`, strict), so that path floors
  at **`max(zstar_h_min, 2·H_VANISHED)`** instead. Never hoist one into the
  other — and never reach for `H_VANISHED` as a positivity floor (D4 forbids
  it); the live-minimum-thickness knob is `&ocean_isopycnal_nml angstrom_h`.
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
