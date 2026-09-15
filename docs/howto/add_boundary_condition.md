# How to add a new boundary condition

There are eleven OBC types: `OBC_WALL=1`, `OBC_OPEN=2` (Flather), `OBC_TIDAL=3`,
`OBC_NESTED=4`, `OBC_INFLOW=5`, `OBC_DISCHARGE=6`, `OBC_CLAMPED=7`,
`OBC_SPONGE=8`, `OBC_CHAPMAN=9` (Orlanski), `OBC_PERIODIC=10`,
`OBC_TRIPOLAR_FOLD=11`, plus the `OBC_INVALID=-1` parse sentinel. Adding one
touches the shared types, the barotropic substep dispatch, the baroclinic
per-layer apply, and the config.

## The seam — module split + dispatch

- **Types/constants**: `src/core/ocean/boundary/rdb_ocean_boundary_types.F90`
  defines `OBC_*`, `ocean_bc_face_tag_t` (the per-edge descriptor: tag, clamped
  η/u/v + `clamped_tracer(:)`, tidal constituents + nodal factors, sponge width /
  strength), `ocean_bc_state_t` (the four edge tags + the periodic / north-fold
  flags, the `has_west/east/south/north` physical-edge gates under MPI, the
  Chapman persistent `eta_old_*`, the Orlanski `rx_*` / `u_prev_*` running state,
  the OBC tracer reservoirs `tres_*`, and the exterior/data-backend arrays),
  `ocean_bc_state_init` / `_destroy` / `_enter_data` / `_exit_data` /
  `_set_edges`, and `ocean_bc_type_from_string`.
- **`ocean_bc_type_from_string` is fail-loud**: it is case-insensitive and an
  unrecognised name returns `OBC_INVALID`, which `validate_config` rejects by
  edge name. A typo must never silently close the boundary to a wall.
- **Barotropic dispatch**:
  `src/core/ocean/kernels/barotropic/rdb_barotropic_substep.F90`.
  Inside the forward-backward substep loop there is one
  `select case (bc_<side>)` per edge (west/east/south/north), each guarded by
  the matching `has_<side>` flag, with arms for `OBC_OPEN`,
  `OBC_TIDAL`/`OBC_CHAPMAN` (they share the Flather algebra with an η target),
  `OBC_CLAMPED`, `OBC_PERIODIC` (leave the value; the ghost wrap overwrites it),
  and `case default` = wall (zero normal transport).
- **Baroclinic / per-layer dispatch**:
  `src/core/ocean/boundary/rdb_ocean_obc_baroclinic.F90` —
  `ocean_obc_apply_baroclinic` (per-layer normal velocity: Orlanski radiation,
  clamped, asymmetric in/out nudging, with the `snapshot_u_prev_*` running
  state), `ocean_obc_fill_ghosts`, `ocean_obc_refill_ghost_ssh`, and
  `ocean_obc_update_reservoirs` (the tracer reservoirs behind an open edge).
  Each is written as per-edge `pure` helpers (`apply_orlanski_west`,
  `apply_clamped_zonal_east`, `apply_nudge_meridional_north`, …).
- **Mass flux**: `src/core/ocean/kernels/continuity_ppm/rdb_continuity.F90`
  reads the edge tags directly to zero the wall-face mass flux and to decide
  whether the column renormalisation may skip closed edges.
- **Separate kernels, not BC arms**: sponge relaxation
  (`rdb_ocean_sponge.F90`), periodic ghost wrap (`rdb_ocean_periodic.F90`), and
  the tripolar north fold (`rdb_ocean_fold.F90` + `rdb_ocean_fold_apply.F90`)
  are their own passes — the `OBC_SPONGE` / `OBC_PERIODIC` / `OBC_TRIPOLAR_FOLD`
  tags select them, they are not implemented inside the substep `select case`.
- **Boundary data backends**: `rdb_ocean_boundary_data.F90` defines
  `ocean_boundary_data_source_t` with `ocean_boundary_data_constant_t` as the
  shipped implementation; file-backed backends are the documented follow-up.

## Recipe (ordered) — add `OBC_FOO`

1. `rdb_ocean_boundary_types.F90` — add `integer, parameter, public :: OBC_FOO = <n>`;
   add a `case ("foo")` in `ocean_bc_type_from_string`. If FOO needs per-edge
   config scalars, add fields to `ocean_bc_face_tag_t`; if it needs persistent
   ghost state, mirror the Chapman `eta_old_*` / Orlanski `u_prev_*` pattern —
   allocate it in `ocean_bc_state_init` and map it in
   `ocean_bc_state_enter_data` (with the matching `exit_data` + `destroy`).
2. `rdb_barotropic_substep.F90` — add a `case (OBC_FOO)` to **all four** edge
   `select case` blocks. Keep the ghost-only discipline: write the boundary face
   / ghost cell, never an interior value.
3. `rdb_ocean_obc_baroclinic.F90` — if FOO needs a per-layer treatment, add the
   arm in `ocean_obc_apply_baroclinic` plus the four `pure` per-edge helpers,
   and decide whether it participates in `ocean_obc_fill_ghosts` and the tracer
   reservoirs.
4. `rdb_continuity.F90` — decide whether FOO is closed to mass flux (the
   `OBC_WALL` zeroing) and whether the renormalisation may skip it.
5. `rdb_config.F90` — add any new `&ocean_bc_nml` keys to `read_ocean_bc_nml`
   (declare on `ocean_bc_config_t`, add to the `namelist /ocean_bc_nml/`
   statement, seed the local from `cfg%`, copy back) per
   [add_namelist_knob.md](add_namelist_knob.md). `&ocean_bc_nml` is the one
   group still on the native reader — it is registered with
   `schema%add_external_group("ocean_bc")` so the strict parser accepts it. Add
   any range/consistency check to `validate_config`.
6. `src/core/ocean/state/rdb_ocean_setup.F90` — copy the parsed config onto the
   edge tags in `configure_ocean_bc`.
7. **GPU map** any new *allocatable* persistent ghost array in
   `ocean_bc_state_enter_data` **and** confirm the slot is reached from the
   `ocean_state_enter_data` orchestrator. A fixed-size component on an
   already-mapped struct rides along automatically.

**Test:** add cases to `tests/test_ocean_boundary.F90` (parse + defaults +
device-mapping lifecycle) and, if FOO has per-layer behaviour, to
`tests/test_ocean_obc_baroclinic.F90`.

## Tests — the pattern to copy

- `tests/test_ocean_boundary.F90` — `bc_state_defaults` (a default state is a
  closed-wall run), `bc_state_destroy`, `type_string_known` (**add your string
  here**), `type_string_fallback` (an unknown name must yield `OBC_INVALID`, not
  a wall), `constant_data_source`.
- `tests/test_ocean_obc_baroclinic.F90` — the physics asserts: outflow
  depth-mean, inflow tracer rise, radiating gravity-wave decay, the
  `wall_bc_regression` bit-identity guard, the reservoir relax rates + their
  bit-identity test, Orlanski `rx` unit and wave drain.
- `tests/test_ocean_obc_eta_ghost.F90` — corner ghost fill for η and tracers,
  and the `refill_ghost_ssh` no-op case.
- `tests/test_ocean_obc_tide_nodal.F90` — the boundary-tide nodal correction.
- `tests/test_ocean_conservation_salt_heat.F90` — `test_open_boundary_out_closes`
  is also the canonical `mem:separate` device-mapping template.

Pattern: a tiny grid, set one edge to your BC, run the pass, assert the
closed-form boundary values — plus a bit-identity test proving that a run
without your BC is byte-for-byte unchanged.

## Gotchas

- **Wall is the default fall-through** in the barotropic `select case` — an
  unhandled tag closes the edge. Keep that invariant for FOO, and make sure the
  intended behaviour is a real `case`, not the default.
- **`OBC_INVALID` is a parse sentinel, never a live edge tag.**
  `validate_config` rejects it (naming the offending edge) before
  `configure_ocean_bc` consumes any parse result.
- **`OBC_SPONGE` selects a separate relaxation kernel** (`ocean_sponge_apply` /
  `ocean_sponge_apply_tracers`) run after the BC pass; its ghost treatment falls
  through to wall. Same shape for `OBC_PERIODIC` (ghost wrap) and
  `OBC_TRIPOLAR_FOLD` (the fold exchange, which reverses i and sign-flips vector
  normals).
- **`has_west/east/south/north`** gate every edge block: under MPI only the ranks
  that own a physical domain edge apply the BC. Never dispatch on the tag alone.
- **Per-edge tracer arrays are sized at `ocean_bc_state_init`** from
  `size(multilayer%tracers)` — a tracer registered after that point is silently
  skipped at every clamped edge. See [add_passive_tracer.md](add_passive_tracer.md).
- **GPU mapping**: a new *allocatable* persistent ghost array must be added to
  `ocean_bc_state_enter_data` / `_exit_data`. On the `mem:separate` build a
  missed map is a silent stale-host read, not a crash.
- **Periodic + fold have validators** — `ocean_bc_validate_periodic` and
  `ocean_bc_validate_fold` enforce the legal edge combinations (the north fold
  requires periodic west/east). Extend them if FOO constrains the edge set.

## Source pointers

- `src/core/ocean/boundary/rdb_ocean_boundary_types.F90` — constants,
  `ocean_bc_face_tag_t`, `ocean_bc_state_t`, `ocean_bc_type_from_string`,
  `ocean_bc_state_init` / `_enter_data` / `_set_edges`, the periodic + fold
  validators, `obc_tide_nodal_fill`.
- `src/core/ocean/kernels/barotropic/rdb_barotropic_substep.F90` —
  the per-edge barotropic dispatch.
- `src/core/ocean/boundary/rdb_ocean_obc_baroclinic.F90` — per-layer apply,
  ghost fill, SSH refill, tracer reservoirs.
- `src/core/ocean/boundary/{rdb_ocean_sponge,rdb_ocean_periodic,rdb_ocean_fold,rdb_ocean_fold_apply,rdb_ocean_boundary_data}.F90`.
- `src/core/ocean/kernels/continuity_ppm/rdb_continuity.F90` — wall
  mass-flux zeroing.
- `src/core/rdb_config.F90` — `read_ocean_bc_nml`, `validate_config`;
  `src/core/ocean/state/rdb_ocean_setup.F90` — `configure_ocean_bc`.
- The "Boundaries" bullet of the top-level [`CLAUDE.md`](../../CLAUDE.md).
