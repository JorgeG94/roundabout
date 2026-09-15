# How to add a new turbulence closure / parameterization

Covers vertical mixing (PP81 interior, KPP, EPBL, kappa-shear, tidal mixing,
convective adjustment) and lateral closures (Leith, Smagorinsky). Modules live
under `src/parameterizations/{vertical,lateral}/`.

Worked references: EPBL, kappa-shear and tidal mixing were all added against the
seam described here and are the templates to copy.

## The seam — the vertical-mixing contributor chain

The ocean path does **not** pick one closure with an `if / else if / else`.
Every closure is a *contributor* that accumulates into the interface-shaped
diffusivity fields `vmix%kv` (momentum) / `vmix%kt` (heat) on
`ocean_vmix_t`, and one downstream gate applies the floors and ceilings. The
chain is assembled in `vmix_apply_in_stage`
(`src/core/ocean/dynamics/split_rk2/rdb_ocean_dyn.F90`), called once per RK2
stage from `run_stage`:

```
PP81 interior (vmix_compute_pp81)
  → KPP overlay (vmix_apply_kpp_overlay)  XOR  EPBL merge (epbl_merge_into_kv_kt)
  → kappa-shear additive merge (kappa_shear_merge_into_kv_kt)
  → tidal-mixing additive merge (tidal_mixing_merge_into_kt)
  → KV_ML_INVZ2 surface band (vmix_add_kv_ml_invz2)
  → convective adjustment (vmix_apply_convection)
  → vmix_split_kd_heat_salt   ! derives ks from kt — NOT a contributor, stays last
  → vmix_assemble             ! the SINGLE gate: backgrounds, kv_max/kd_max,
                              ! optional smoothing + guard; applies to kv, kt, ks
  → vdiff (vdiff_apply_momentum / vdiff_apply_tracers)
```

The interior closure itself is selected by `vmix%interior_closure`
(`VMIX_INTERIOR_PP81` = the only implemented one; `VMIX_INTERIOR_LARGE94` /
`VMIX_INTERIOR_CVMIX` are reserved-but-unwired and **fail loud**, they are not
silent fall-throughs). It has no namelist key today — `&ocean_vmix_nml
use_closure` is the master on/off switch.

### Interface contract a vertical closure must satisfy

- **Inputs** (read from the state slots): per-layer `h_layer`, `hu_layer`,
  `hv_layer`, `rho_layer` (shear + stratification/Ri) from
  `multilayer_state_t`; surface stress (`ocean_surface_stress_t`) and
  surface buoyancy/heat/salt fluxes (`ocean_surface_flux_t`) for any
  boundary-layer scheme.
- **Outputs**: `kv` / `kt` at **interfaces** (shape `(nx, ny, nz+1)`), which
  `vdiff` then applies implicitly (backward-Euler tridiagonal). Do **not**
  apply your own explicit diffusive flux and do **not** clamp — clamping is
  `vmix_assemble`'s job, and doing it upstream double-applies the floor.
- **`ks` is derived, never written by a closure.** `vmix_split_kd_heat_salt`
  produces `ks` from `kt` (with the `&ocean_ddiff_nml` double-diffusion
  asymmetry folded in when enabled) as the last statement before the gate. A new
  closure writes `kv`/`kt` only.
- **Persistent state** lives on your own slot type (e.g. `ocean_epbl_t`,
  `ocean_kappa_shear_t`), which owns its allocatables + bound
  `init`/`destroy`/`enter_data`/`exit_data`, and is composed onto
  `ocean_state_t`.
- **Thermo cadence.** Expensive prognostic closures refresh their diffusivity
  only on a thermo step (`dyn%enable_thermodynamics .and. dyn%is_thermo_step()
  .and. stage == 1`) and merge the cached field every stage — see the
  `kappa_shear_compute` / `kappa_shear_merge_into_kv_kt` split.

## Recipe (ordered) — a new vertical closure `<name>`

1. **Config knobs** — a new `&ocean_<name>_nml` group with `enable = .false.`
   per [add_namelist_knob.md](add_namelist_knob.md). Default off ⇒ bit-identical.
2. **Slot type + fields** — `ocean_<name>_t` in
   `src/parameterizations/vertical/rdb_ocean_<name>.F90`, with the
   interface-shaped (`nz+1`) persistent fields it needs and bound
   `init`/`destroy`/`enter_data`/`exit_data`.
3. **Compose onto `ocean_state_t`** (`src/core/ocean/state/rdb_ocean_state.F90`)
   and **wire it into the `ocean_state_enter_data` orchestrator** — a missed
   array there is a 150–1500× memcpy explosion, not a crash.
4. **Kernel** — `pure` subroutines, explicit-shape args `arr(nx, ny, nz)`,
   `do concurrent(j, i)` over columns (j outer, i inner), bottom-up `k`.
   Split it as `<name>_compute` (refresh, thermo cadence) +
   `<name>_merge_into_kv_kt` (per-stage accumulate) if it is expensive.
   Register the file in `src/CMakeLists.txt`.
5. **Wire the chain** — add the `present(<name>) .and. <name>%enable` block in
   `vmix_apply_in_stage`, in the right place in the ordering above (before
   `vmix_split_kd_heat_salt`, always before `vmix_assemble`), and thread the
   optional slot argument through `ocean_dyn_step` → `run_stage` →
   `vmix_apply_in_stage`. Decide and document precedence vs KPP/EPBL (KPP and
   EPBL are mutually exclusive — that exclusion is enforced fail-loud at
   configure).
6. **Configure step** — a `configure_ocean_<name>(cfg, ocean_state, ...)` in
   `src/core/ocean/state/rdb_ocean_setup.F90` copies `cfg%ocean%<name>%*` onto
   the slot and does the cross-knob validation the schema cannot express
   (mutual exclusions, "requires X" checks). Fail loud; never silently downgrade.
7. **Test + matrix row** — add `tests/test_ocean_<name>.F90` to
   `tests/CMakeLists.txt` with regime `ocean` (see the ctest row format in
   [add_namelist_knob.md](add_namelist_knob.md)), and a row to
   `docs/CLOSURE_MATRIX.md` (the drift hook fails the build if a cited test is
   unregistered, and the curated-knob list there covers `&ocean_hvisc_nml`,
   `&ocean_bdrag_nml` and `&ocean_kappa_shear_nml` by name).

## Lateral closures

Under `src/parameterizations/lateral/`:

- `rdb_ocean_lateral_mix.F90` — the flow-aware closures, dispatched by
  `select case (this%closure)` over `LMIX_NONE` / `LMIX_LEITH` /
  `LMIX_SMAGORINSKY` / `LMIX_BIHARMONIC` / `LMIX_LEITH_BIHARM`, parsed from
  `&ocean_hvisc_nml lateral_closure` by `parse_lateral_closure`. The
  dispatcher is **fail-loud**: `lateral_closure_is_implemented` /
  `LMIX_INVALID` mean a typo aborts rather than silently selecting `LMIX_NONE`.
  Harmonic closures fill per-face `ah_face_x/y`; the biharmonic one fills
  `nu4_face_x/y`.
- `rdb_ocean_horizontal_viscosity.F90` — the apply step (stress tensor,
  per-cell CFL clamps, the biharmonic add-on, anisotropy), plus
  `rdb_ocean_meke.F90` for the backscatter budget.

A new lateral closure follows the same shape: add the `LMIX_*` tag, extend
`parse_lateral_closure` **and** `lateral_closure_is_implemented` together (they
must stay in sync — the second is what keeps the fail-loud honest), add the
per-face coefficient kernel, and gate it. Remember the standing rule: a
biharmonic backstop with a NON-ZERO dissipation coefficient is mandatory under
MEKE backscatter.

## Tests — the pattern to copy

- `tests/test_ocean_kpp.F90` — analytic shape-function check;
  **`kpp_no_wind_is_no_op`** and **`kpp_vt2_no_forcing_invariant`** (the
  default-off / off-forcing bit-identity invariant — the canonical proof that a
  new overlay collapses to the base closure when inactive); BL depth in a
  stratified vs unstratified column. Also `tests/test_ocean_kpp_convective.F90`,
  `tests/test_ocean_kpp_nonlocal.F90`.
- `tests/test_ocean_vmix_assembly.F90` — the gate itself: defaults pass through
  untouched, `kd_max` clamps, floors apply, 1-2-1 smoothing stencil, the
  negative/NaN guard, `ks` derived from `kt`. **Any new contributor must leave
  these green.**
- `tests/test_ocean_epbl.F90`, `tests/test_ocean_kappa_shear.F90`,
  `tests/test_ocean_tidal_mixing.F90` — per-closure analytic + merge tests.
- `tests/test_ocean_convection.F90` — the model for a `max()`-style
  contributor: `conv_disabled_bit_identical`, `conv_stable_column_bit_identical`,
  then the physics asserts.
- Layer-convention gates every closure must keep green:
  `test_ocean_surface_flux` (positive Q warms `k=nz`), `test_ocean_geothermal`
  (bed heat at `k=1`), `test_ocean_conservation_salt_heat`.

Implementations here are clean-room from the published equations — **not**
derived from MOM6, GOTM or any other codebase's source. Keep new closures the
same way: cite the paper, not another codebase. Validate against an analytic
limit (decay law, log layer, similarity profile) where one exists.

## Gotchas

- **No local-allocate scratch in a per-step kernel on `-stdpar=gpu`.** Use either a
  module-level allocatable + lazy `*_workspace_ensure` + `*_cleanup`, or fields on
  the slot mapped once in `enter_data`. Per-column scratch is fixed-size
  `NZ_STACK_MAX` (a `local()` clause requires fixed size; a dummy-sized automatic
  array crashes with `CUDA_ERROR_ILLEGAL_ADDRESS`).
- **`do concurrent(j, i)` ordering** (j outer) for coalescing; **explicit-shape**
  kernel args, never assumed-shape.
- **`!$acc routine seq`** on cross-module `pure` helpers called from the device loop;
  for hot inner loops NVHPC won't inline them — keep helpers in-module or duplicate
  as `_impl`.
- **Bottom-up `k`** (`k=1` bed, `k=nz` surface); interface `k` sits between layers
  `k` and `k+1`; surface forcing at the top interface.
- **Don't clamp in your own kernel** — `vmix_assemble` is the single floors /
  ceilings / smoothing gate over `kv`, `kt` and `ks`. Contributing a
  pre-floored field silently changes the assembled result.
- **Don't write `ks`** — it is derived from `kt` by `vmix_split_kd_heat_salt`,
  which must stay the last statement before the gate.
- **A host-gated call that hands a state array to an EXTERNAL subroutine costs
  even when never taken** (measured +4.8% on `ocean_continuity`): nvfortran treats
  the array as escaping and pessimises every `do concurrent` in the calling
  routine. Write a guarded pass inline as a `do concurrent` in the same routine.
- **The CLOSURE_MATRIX can lag the code.** Treat the chain in
  `vmix_apply_in_stage` as the source of truth for what's wired; update the matrix
  row when you add a closure.

## Source pointers

- `src/core/ocean/dynamics/split_rk2/rdb_ocean_dyn.F90` — `vmix_apply_in_stage`
  (the contributor chain + the fail-loud interior-closure guard), `run_stage`.
- `src/parameterizations/vertical/rdb_ocean_vmix.F90` — `ocean_vmix_t`,
  `vmix_compute_pp81`, `vmix_apply_kpp_overlay`, `vmix_add_kv_ml_invz2`,
  `vmix_apply_convection`, `vmix_split_kd_heat_salt`, `vmix_assemble`,
  `VMIX_INTERIOR_*`.
- `src/parameterizations/vertical/rdb_ocean_vdiff.F90` — the implicit
  solve that consumes `kv`/`kt`/`ks`.
- `src/parameterizations/vertical/{rdb_ocean_epbl,rdb_ocean_kappa_shear,rdb_ocean_tidal_mixing}.F90`
  — worked contributor examples.
- `src/parameterizations/lateral/{rdb_ocean_lateral_mix,rdb_ocean_horizontal_viscosity,rdb_ocean_meke}.F90`.
- `src/core/ocean/state/rdb_ocean_setup.F90` — `configure_ocean_vmix`,
  `configure_ocean_epbl`, `configure_ocean_kappa_shear`,
  `configure_ocean_tidal_mixing`, `configure_ocean_conv`, `configure_ocean_lateral`.
- `docs/CLOSURE_MATRIX.md` — which closure is on + the knobs;
  `docs/CLOSURE_MATRIX.md` — the per-closure detail (drift-checked).
