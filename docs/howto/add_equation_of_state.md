# How to add a new equation of state

One EOS module: `src/equation_of_state/rdb_eos.F90`. It already has a
variant enum + `select case` dispatch (linear, Wright 1997 nonlinear, Roquet
specific-volume polynomial), a namelist selector, and a device-callability gate.
Adding a variant is filling in an arm plus wiring the parse.

## The seam

`ocean_eos_compute(eos, ms, active)`
(`src/core/ocean/state/rdb_ocean_eos_compute.F90`) is the outer shim: it pulls
`ms%tracers(idx)%hTr` off the registry on the *host* and hands flat arrays to
`eos_compute_arrays` in `rdb_eos.F90`, which does
`select case (eos%variant)` over `EOS_VARIANT_LINEAR` / `EOS_VARIANT_WRIGHT_97`
/ `EOS_VARIANT_ROQUET_SPV` (`EOS_VARIANT_TEOS10` is reserved and **fails loud**
in `eos_validate` — it is not device-callable). Each arm calls a `pure *_impl`
(`eos_linear_impl`, `eos_wright_impl`, `eos_roquet_spv_impl`) that is a single
`do concurrent(k,j,i)` with a vanishing-layer fallback to `rho_0`. The variant
lives on `eos_t%variant`.

**Config wiring** — `&ocean_eos_nml eos` (default `"linear"`) is parsed by
`parse_eos_variant` and assigned in `configure_ocean_*`
(`src/core/ocean/state/rdb_ocean_setup.F90`):
`ocean_state%eos%variant = parse_eos_variant(cfg%ocean%eos%eos)`, followed by
`eos_validate`. Note `parse_eos_variant` falls back to `EOS_VARIANT_LINEAR` on
an unrecognised string (bit-identity default); the *fail-loud* gate is
`eos_validate`, over the device-callable set.

**Point-wise entry points** — `eos_density_point`, `eos_specvol_derivs` and
`eos_freezing_point` are the `!$acc routine seq` handles other kernels (vcoord
`VCOORD_RHO` inversion, sea ice, diagnostics) call per point. A new variant must
extend those too, or those consumers silently keep the old formula. A device
point routine cannot `error stop` — that is exactly why `eos_validate` gates the
variant set on the host.

**Pressure coupling.** A pressure-dependent EOS couples to the pressure-gradient
force. `ocean_eos_compute` at `p_ref = 0` produces only a *seed* surface density;
the real in-situ density is produced by `eos_wright_pgf_column_sweep_impl` during
PGF Pass 1, gated on the *PGF's own* variant (`OPGF_VARIANT_FV_WRIGHT`) — a
one-Picard-step top-down column sweep that re-evaluates ρ at each layer's
half-pressure. EOS runs before the PGF every slow tendency. So a new
pressure-dependent EOS needs *two* arms: an EOS variant and a matching PGF
variant + column sweep.

## Recipe — a new variant

1. **Enum** — add `EOS_VARIANT_<NAME>` as a public `integer, parameter` in
   `rdb_eos.F90`.
2. **Kernel impl** — add `pure subroutine eos_<name>_impl(...)` modelled on
   `eos_linear_impl`: explicit-shape 3D arrays, one `do concurrent(k,j,i)`,
   vanishing-layer fallback to `rho_0`. If it's pressure-dependent, also add a
   public `eos_<name>_pgf_column_sweep_impl` modelled on the Wright one.
3. **Dispatch arm** — add `case (EOS_VARIANT_<NAME>)` inside the *existing*
   `select case (eos%variant)` in `eos_compute_arrays` (keep it one select-case
   — see gotchas).
4. **Point routines** — extend `eos_density_point` / `eos_specvol_derivs` (and
   `eos_freezing_point` if the variant changes it), and add the new tag to the
   device-callable arm of `eos_validate`.
5. **New coefficients** on `eos_t` if needed (scalar defaults; any allocatable
   scratch must be released in `eos_destroy` *and* wired into
   `ocean_state_enter_data`).
6. **Namelist** — add the string to the `&ocean_eos_nml eos` enum registration
   in `register_ocean_eos` (`src/core/rdb_config.F90`) and the matching
   `case` in `parse_eos_variant` — the two must stay in sync (the registration's
   FORD comment says so explicitly). See
   [add_namelist_knob.md](add_namelist_knob.md).
7. **PGF arm** (only if pressure-dependent) — add an `OPGF_VARIANT_*` and the
   column-sweep call in
   `src/pressure_force/rdb_ocean_pressure_force.F90`, and extend
   `parse_opgf_variant`.
8. **Test** — copy `tests/test_ocean_wright_eos.F90` (nonlinear) or
   `tests/test_ocean_eos_roquet.F90` (which also covers the fail-loud and
   default-off cases); register the row in `tests/CMakeLists.txt` with regime
   `ocean`.

## Tests — the pattern to copy

Analytical point-value asserts, one per branch of the formula:

- `tests/test_ocean_eos.F90` (linear): `eos_at_reference` — ρ at reference S/T
  equals `rho_0` to ~1e-12; `eos_pure_salinity_anomaly` /
  `eos_pure_temperature_anomaly` equal `±β·dS` / `−α·dT`; `eos_combined_per_layer`
  with distinct `dS_k`/`dT_k` + non-uniform `h` (guards cross-layer index bugs).
  The harness runs `enter_data → compute → exit_data` so the GPU path is
  exercised.
- `tests/test_ocean_wright_eos.F90` (nonlinear): `wright_reference_seawater` (a
  known reference value within a loose tolerance), `wright_pure_water`,
  `wright_monotonic_TS` (warmer → lighter, saltier → denser), `wright_cabbeling`
  (a blend of two equal-density parcels is denser than the linear mean — the
  discriminator a linear EOS *fails*), and `wright_vanishing_layer`
  (`h ≤ 0 → rho_0`).
- `tests/test_ocean_eos_roquet.F90`: published check values, chain-rule
  derivatives, **`roquet_default_off_bit_identity`**, **`roquet_fail_loud_fv_wright`**
  (an unsupported EOS×PGF pairing must abort, not silently mismatch),
  fresh-water finiteness, and `roquet_3d_fill_matches_point` — the 3D fill and
  the point routine must agree.
- `tests/test_ocean_eos_handle.F90`: the point-routine handles (`eos_density_point`,
  `eos_specvol_derivs`) locked against reference values per variant, consumer
  consistency, and `eos_handle_fail_loud_variant_set`.

## Gotchas

- **One `select case` inside one `do concurrent`** — splitting variant dispatch into
  one DC-launch per case is ~1.8× slower; a uniform branch is ~free.
- **`pure` + `!$acc routine seq`** for any helper called from the device loop; NVHPC
  won't inline cross-module `acc routine seq` helpers in hot kernels (the Wright
  coefficients are inlined as module `parameter`s for exactly this reason).
- **The 3D fill and the point routine must not drift.** They are separate code
  paths over the same formula; `roquet_3d_fill_matches_point` is the test that
  keeps them honest — write the equivalent for your variant.
- **Pressure dependence lives on the PGF integration path**, not just the surface
  `rho_layer`. A new pressure-dependent EOS needs an analogous column sweep + a
  matching `OPGF_VARIANT_*`, and the unsupported pairings must fail loud.
- **Reference-density consistency** — the FV-PGF subtracts a reference density
  so the integrand is the small anomaly; a new EOS must remain consistent with
  the reference state it feeds.
- **Vanishing-layer fallback is mandatory** — every impl guards `h ≤ 0` and falls
  back to `rho_0` / reference T,S (required for ZSTAR_FULL and the wet/dry paths
  where bed-side layers pinch out).
- **New allocatable EOS scratch** must be freed in `eos_destroy` *and* wired into
  `ocean_state_enter_data` (silent omission → 150–1500× memcpy explosion).
- **T/S convention**: `eos_t` carries the `TS_POT_PRAC` / `TS_CONS_ABS` tag. A new
  variant must state which convention its coefficients are fitted in — mixing
  potential/practical with conservative/absolute is a silent bias, not an error.

## Source pointers

- `src/equation_of_state/rdb_eos.F90` — `eos_t`, `EOS_VARIANT_*`,
  `parse_eos_variant`, `eos_validate`, `eos_compute_arrays` + the `*_impl`
  kernels, the Wright PGF column sweep, `eos_density_point`,
  `eos_specvol_derivs`, `eos_freezing_point`.
- `src/core/ocean/state/rdb_ocean_eos_compute.F90` — the registry outer shim.
- `src/pressure_force/rdb_ocean_pressure_force.F90` — `OPGF_VARIANT_*`,
  `parse_opgf_variant`, the in-situ re-eval call.
- `src/core/ocean/dynamics/split_rk2/rdb_ocean_dyn.F90` — EOS-before-PGF call order.
- `src/core/rdb_config.F90` — `register_ocean_eos` (`&ocean_eos_nml`);
  `src/core/ocean/state/rdb_ocean_setup.F90` — the `parse_eos_variant` assignment.
- The "Equation of state" section of [`docs/CLOSURE_MATRIX.md`](../CLOSURE_MATRIX.md).
