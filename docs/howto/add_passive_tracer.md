# How to add a new passive tracer

Scope: the ocean C-grid path (`multilayer_state_t`) — the only regime this
build ships.

This is the easiest extension in the codebase, by design. The multilayer tracer
pipeline iterates a *registry*: append a tracer with `register_passive_tracer(...)`
and every transport kernel — advection, ALE remap, vertical exchange,
vertical/horizontal diffusion, Redi, sponge, halo/fold, OBC ghost + reservoirs,
restart, conservative clamp, RK2 averaging — picks it up automatically. You write
zero kernel code.

## The seam

`multilayer_state_t` (`src/core/rdb_multilayer_state.F90`)
owns `tracers(:)`, an allocatable array of `tracer_t`
(`src/core/rdb_tracer.F90`). Salinity and temperature are the first
two entries (`idx_salinity`, `idx_temperature`), optionally followed by ideal
age (`idx_age`) and pseudo-salt (`idx_pseudo_salt`); a new passive tracer
appends after them. Every passive-transport kernel loops
`do it = 1, size(tracers)`, so a new entry "just works". What stays S/T-specific
by design: the EOS and the surface fluxes (a passive tracer has `eos_coeff = 0`,
so it has zero density feedback by construction).

`register_passive_tracer` has a deliberately **narrow 6-argument** signature —
`(grid, name, units, long_name, idx)` on the bound type — with no
`tr_min`/`tr_max`/`kappa_bg`/`hdiff_kappa` (those are dead on this path; see
`FORTRAN_STYLE.md` §Limit Procedure Arguments) — and a **lock** that closes the
`enter_data` foot-gun with a message instead of a silent GPU fault.

## Recipe (ordered)

1. **Call `register_passive_tracer` between `multilayer%init` and
   `ocean_bc_state_init`** — inside `ocean_state_init`
   (`src/core/ocean/state/rdb_ocean_state.F90`), right after
   `call this%multilayer%init(grid, ...)`. **This ordering is load-bearing**:
   `ocean_bc_state_init` sizes the per-edge `bc%*%clamped_tracer(:)` arrays
   and the OBC tracer reservoirs from `size(this%multilayer%tracers)` at that
   moment — register after it and the new tracer is silently skipped at every
   `OBC_CLAMPED` edge (a bug that looks like clean physics: `D ≈ 0` in the
   interior, drifting at the boundary).
   ```fortran
   call ms%register_passive_tracer(grid, "dye", "1", "Dye tracer", idx)
   if (idx <= 0) error stop "registration refused"   ! callers MUST check
   ```
   `idx` is `intent(out)`: the new slot number, or **`0` on refusal**
   (registry not initialised, or already locked by `enter_data`). The call
   forces `eos_coeff = 0`, `eos_ref = 0`, `budget_id = TRACER_BUDGET_NONE`,
   `standard_name = ""`; set any other public component (the pipeline opt-outs
   `do_horizontal_advection` / `do_vertical_exchange` /
   `do_vertical_diffusion` / `do_horizontal_diffusion` / `do_clamp`,
   `standard_name`) on `tracers(idx)` afterwards. `idx_salinity` /
   `idx_temperature` / `idx_age` are unchanged.
2. **`registry_locked`** — set `.true.` by `enter_data`, cleared by
   `exit_data`. A registration attempt while locked refuses (`idx = 0`,
   logged) instead of leaving an unmapped `hTr` on the `mem:separate` GPU
   build.
3. **Seed `hTr`** — manual: after layer thicknesses / any analytical or z-file
   IC land, `ms%tracers(idx)%hTr = ms%h_layer * tr_init` (or copy from another
   tracer, as `ocean_pseudo_salt_seed` does from salinity — see below).
   Must run before `enter_data`.
4. **Budget attribution is a registry property, not an index coincidence.**
   `tracer_t%budget_id` (`TRACER_BUDGET_NONE` default / `_HEAT` / `_SALT`)
   drives the ocean budget-dispatch blocks
   (`select case (ms%tracers(it)%budget_id)`) that fill
   `heat_budget_*`/`salt_budget_*`. A passive tracer stays at the default
   `NONE` — it must NOT opt into a budget slot unless it genuinely is a
   heat/salt-equivalent flow (pseudo-salt deliberately stays `NONE`; only S/T
   get `_SALT`/`_HEAT` in `multilayer_state_init`).
5. **Diagnostic output is NOT registry-driven**: `fill_salinity` / `fill_age` /
   … each name their index explicitly in `rdb_ocean_diag_fills.F90`. A named
   tracer package needs its own `fill_<name>` + a `register_one_canonical(...)`
   call guarded on `idx_<name> > 0` (mirror `fill_age` / `idx_age`) to show up in
   NetCDF output. This is the one step the registry does not do for you.
6. **Knob + test** — gate the package on a `&ocean_tracers_nml` knob defaulting
   to `.false.` per [add_namelist_knob.md](add_namelist_knob.md), and add a
   `tests/test_ocean_<name>.F90` row to `tests/CMakeLists.txt` with regime
   `ocean`.

## Why the order matters

`register_passive_tracer` grows `tracers(:)` via `move_alloc` + reallocation
(deep-copying the existing entries, which reallocates their `hTr`/`hTr0`). The
slot's `enter_data` snapshots the array and maps each `hTr`/`hTr0` to the device
in a two-step map (`copyin` the array descriptor, then a per-element loop
mapping each `hTr`). Register *after* `enter_data` and the new slot's `hTr`
would never be mapped — hence `registry_locked` refusing rather than letting you
hit `CUDA_ERROR_ILLEGAL_ADDRESS` on the first kernel.

## What you must NOT do

- Don't write `eos_coeff` — it must stay `0`, or the EOS would read density off
  an unbacked tracer.
- Don't assume `idx_salinity` / `idx_temperature` / `idx_age` shift — they don't.
- Don't register after `enter_data` (it will refuse), or after
  `ocean_bc_state_init` (it will silently under-size the boundary arrays).
- Don't ignore the returned `idx` — `0` means refused.

## Worked example — pseudo-salt

`src/tracer/rdb_ocean_pseudo_salt.F90` is the reference consumer: a
verification tracer seeded to salinity and given exactly salinity's surface salt
flux + KPP/EPBL non-local mirror (the two operators the registry does NOT deliver
automatically — everything else, it rides for free).
`&ocean_tracers_nml enable_pseudo_salt` (default `.false.`) gates it end to end:
`ocean_pseudo_salt_register` (setup, between `multilayer%init` and
`ocean_bc_state_init`), `ocean_pseudo_salt_seed` (host, after every S write,
before `enter_data`), and the two mirrored operators self-gate inside
`ocean_surface_flux_apply_tracers` / `vmix_apply_nonlocal_tendencies` on
`ms%idx_pseudo_salt > 0` — no new dyn-step call site.
`src/tracer/rdb_ocean_ideal_age.F90` is the second example (a source
term instead of a flux mirror).

## Test — the pattern to copy

`tests/test_ocean_pseudo_salt.F90` is the pattern for a new tracer package,
including the `idx = 0`-refused / locked-registry / GPU device-indirection
canary tests.

The canonical *transport* check is a **tracking** test, not a sum check:
register a "dye" with zero density coupling, seed it equal to salinity, run a
short integration, and assert `maxval(|hTr_dye - hTr_S|)` stays at round-off.
The dye must track salinity bit-for-bit through advection + remap + vertical
exchange + clamp + RK2. A sum check would falsely pass if a kernel silently
skipped the tracer; the tracking check won't.

## Gotchas

- **The `bc%n_tracers` ordering trap** — see step 1; this is the single most
  likely way to get a passive tracer wrong.
- **`registry_locked`, not prose** — the registry refuses a post-`enter_data`
  registration instead of only documenting the hazard.
- **No automatic diagnostic output** — step 5; forgetting the `fill_<name>`
  means the tracer transports correctly but never appears in NetCDF output.
- **Outer-shim + flat-impl device indirection.** A kernel can't dereference
  `ms%tracers(it)%hTr(...)` from inside a device loop (the array-of-derived-types
  descriptor isn't reachable on the device). The pattern: an outer shim
  dereferences the registry on the *host* and hands a flat array to a
  `pure *_impl` worker whose dummy is a plain `arr(nx,ny,nz)`. Every registry
  kernel already does this — you only need to know it if you write a *new*
  registry-iterating kernel (you don't, for a plain passive tracer).

## Source pointers

- `src/core/rdb_tracer.F90` — `tracer_t`, `TRACER_BUDGET_*`, the
  pipeline opt-out flags.
- `src/core/rdb_multilayer_state.F90` —
  `register_passive_tracer`, `registry_locked`, the index handles.
- `src/core/ocean/state/rdb_ocean_state.F90` — `ocean_state_init`: the
  `multilayer%init` → register → `ocean_bc_state_init` ordering.
- `src/tracer/{rdb_ocean_pseudo_salt,rdb_ocean_ideal_age}.F90` — the
  two shipped packages.
- `src/core/ocean/diag/rdb_ocean_diag_fills.F90` — `fill_*` +
  `register_one_canonical`.
- The "Tracer registry" bullet of the top-level [`CLAUDE.md`](../../CLAUDE.md).
