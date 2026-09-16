# How to add a new diagnostic

Scope: the **ocean** path (`sim_type='ocean'`). The diag manager lives in
`src/core/ocean/diag/`.

The ocean diagnostics are a *registry* with cadence dispatch and GPU-resident
`fill_*` kernels: each diagnostic computes on-device into its `output_buffer`, and
**only the cadence fire pulls host←device**. (Per-step `update self` was measured at
3× overhead — don't reintroduce it.)

## The seam

- `ocean_diag_t` holds `vars(:) : diag_var_t` (capacity-doubling), the fixed-z
  output levels, a NetCDF stream, and an `emit_post_fire` proc pointer.
- `diag_var_t` carries name/units/standard_name, a `fill` proc pointer (and an
  optional `remap`), a time-op (`DIAG_OP_INSTANT/MEAN/MAX/MIN/INTEGRAL`), a `dt_out`
  cadence, and the buffers (`output_buffer`, `layer_buffer`, `accumulator`).
- `ocean_diag_register(...)` grows the registry and allocates the buffers; for a
  `LAYER`-vgrid var the `fill` writes `output_buffer` directly, otherwise it writes
  `layer_buffer` and a `remap` produces `output_buffer` on the fixed-z grid.
- `ocean_diag_step(...)` accumulates, decides `fire`, runs fill/fold/finalise
  **on-device**, and only then does the single `!$acc update self(v%output_buffer)
  if_present` before the host-side min/max/sum log line and the NetCDF emit.

The `fill_*` pattern is outer-shim + flat-impl: the public `fill_<x>` casts
`class(*) → ocean_state_t`, dereferences any tracer-registry indirection host-side,
and forwards explicit-shape device arrays to a `pure fill_<x>_impl` running the work
as `do concurrent`. (See `fill_ke` / `fill_ke_impl` for the cleanest example.)

## Recipe (ordered) — a new derived field

1. **Write the fill** in `rdb_ocean_diag_derived.F90` (opt-in, name-keyed) or
   `rdb_ocean_diag_fills.F90` (always-on default):
   - a shim `fill_<name>` (host `select type`) + a `pure fill_<name>_impl`
     (`do concurrent`), modelled on an existing derived fill;
   - explicit-shape args, `local()` + fixed-size `NZ_STACK_MAX` for any column
     scratch;
   - add it to the module `public ::` list;
   - for a derived field, add a `CATALOG` entry in `ensure_catalog_initialised`
     (name/units/long_name/standard_name/`is_layered`) and bump the count.
2. **CMake** — no edit if you reuse an existing diag module (all listed already); a
   brand-new module file must be added to `src/CMakeLists.txt`.
3. **Register it** — either it flows through `register_derived` by name (driver/test
   opt-in), or for an always-on field add a `call state%diag%register(...)` in
   `register_default_diags`.
4. **Cadence / time-op** — pass `time_op=` and `dt_out=` to `register`. For fixed-z
   output, call `set_output_z_levels` *before* register and pass
   `output_vgrid=DIAG_VGRID_Z_FIXED, remap=remap_layer_to_z` (intensive) or
   `remap=remap_layer_to_z_extensive` (extensive). For density-space output, call
   `set_output_density_levels` *before* register and pass
   `output_vgrid=DIAG_VGRID_DENSITY, remap=remap_layer_to_density` (or
   `remap_layer_to_density_extensive`); `z_out` then carries monotone target
   densities and the manager threads `rho_out` to the remap automatically.
5. **Nothing else.** `ocean_diag_step` and the NetCDF emit are generic over the
   registry — no step/emit edits.

## The ordering gotcha (read this)

Buffers are allocated at *register* time; `ocean_diag_enter_data` walks `vars(:)`
and `!$acc enter data copyin`s each buffer, and is itself called from
`ocean_state_enter_data`. The driver calls `configure_ocean_diag` (register +
open_stream) **before** `ocean_state_enter_data`. So:

> **You must register before `enter_data` runs.** You do *not* add a line to the
> orchestrator for a diag buffer — `diag%enter_data` already loops every registered
> var. But register *after* `enter_data` and your buffer is never device-attached →
> the `do concurrent` in your `_impl` reads a fresh empty descriptor → garbage /
> 150–1500× memcpy explosion.

## Tests

- `tests/test_ocean_diag.F90` — registry growth (capacity doubling preserves names),
  cadence (no fill below cadence; fill + reset on fire), the default six-var set
  (nvars, first = SSH, last = KE), `fill_ssh`/`fill_ke` value checks, the
  conservative z-remap (cell-average consistency / whole-column clip / e2e), all
  four time-ops (mean/max/min/integral + `instant_no_accum`), and the derived
  catalog (vorticity of rigid rotation = 2Ω, ke_total, transport, mld).
- `tests/test_ocean_diag_remap.F90` — conservative z-remap (intensive consistency +
  integral preservation; extensive column-sum conservation), density-space remap
  (bin value match + lightest→surface ordering; extensive conservation), and the
  on-device (enter_data) round-trip.
- `tests/test_ocean_diag_netcdf.F90` — 2D/3D var written per fire; close unbinds the
  emit hook.
- `tests/test_ocean_driver_diag.F90` — end-to-end through the driver (default set +
  expected frame count), including the device path.

Pattern: a closed-form field on a synthetic state (rigid rotation, uniform flow,
step profile) with an exact-value assert, plus a cadence test (fires at `dt_out`,
not before) and a time-op test.

## Gotchas

- **GPU-resident fill, single D→H at fire.** Never pull the buffer inside your
  `_impl`.
- **Register before `ocean_state_enter_data`** (above).
- **`output_buffer` lifecycle** is owned by the manager (alloc in register, attach
  in enter_data, detach in exit_data reverse order, dealloc in destroy). The I/O
  server gets a non-owning view.
- **NVHPC stdpar landmines in `_impl`**: explicit-shape (not assumed-shape) args;
  fixed-size stack arrays under `local()` (dummy-sized automatic arrays crash with
  `CUDA_ERROR_ILLEGAL_ADDRESS`); avoid local names that shadow intrinsics (use
  `scale_factor`, not `scale`); one `do concurrent` per time-op case, not
  case-inside-loop.
- **The vertical remap is CONSERVATIVE** (donor-cell overlap integral via
  `remap_column`, reconstruction set by `&ocean_diag_nml diag_remap_scheme` =
  pcm/plm/ppm/ppm_h4, default **ppm**; all schemes conserve the column
  integral) and comes in intensive / extensive variants for both
  z (`remap_layer_to_z{,_extensive}`) and density
  (`remap_layer_to_density{,_extensive}`). Pick the variant matching the var's
  `is_extensive`: intensive remaps the value as a cell-average; extensive divides
  in / multiplies out by thickness so the column integral is preserved. `z_out`
  carries target INTERFACE depths (implicit 0 surface) clipped to the column total;
  density targets reuse the RHO `invert_density_targets` bracket+Newton solve.
- **MPI I/O-server hand-off is pending.** The current writer emits one
  `output_rank_NNNNNN.nc` per process with no gather -- merge offline with
  `tools/merge_output.py`. The manager has the design hooks; a gathering
  I/O-server path is a separate, unscheduled piece of work.

## Source pointers

- `src/core/ocean/diag/rdb_ocean_diag.F90` — manager: `ocean_diag_t`, `diag_var_t`,
  `register`, `step`, enter/exit_data, the enums.
- `src/core/ocean/diag/rdb_ocean_diag_fills.F90` — the six default fills +
  `register_default_diags`.
- `src/core/ocean/diag/rdb_ocean_diag_derived.F90` — the name-keyed derived catalog.
- `src/core/ocean/diag/rdb_ocean_diag_netcdf.F90` — the serial NetCDF emit.
- `src/driver/rdb_driver.F90` — `configure_ocean_diag` (registers + opens stream
  *before* `ocean_state_enter_data`).
