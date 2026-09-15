# How to add a new namelist knob / config option

A knob is the entry point for almost every capability (workflow step 1). Get the
plumbing right once and the rest of the guides reuse it.

All config lives in `src/core/rdb_config.F90`. There are two public config types:

- **`config_t`** — the flat parent: grid, time, physics, tracer, vcoord,
  boundary, MPI, output, logging. Each field is declared *with its default inline*
  (`real(wp) :: nu_h = 0.0_wp`).
- **`ocean_config_t`** — the nested ocean knobs, one sub-type per concern
  (`ocean_coriolis_config_t`, `ocean_hvisc_config_t`, `ocean_vmix_config_t`, …),
  composed onto `config_t` as `type(ocean_config_t) :: ocean`. Ocean knobs live
  under per-concern namelist groups `&ocean_<group>_nml` with the `ocean_` prefix
  dropped from each key (so `nu_h` lives under `&ocean_hvisc_nml`).

Parsing goes through the **strict schema** (`rdb_nml_schema`, type
`nml_schema_t`): one registration line per knob carries the name, a pointer into
the config, optional bounds/enum set, units, and the doc string. That single
list is what the parser, the parameter dumps and `render_markdown` all read
from, and unknown keys / type / range / enum violations fail loud with
did-you-mean suggestions. `&ocean_bc_nml` is the one remaining group on the old
native `read(nml=)` reader (`read_ocean_bc_nml`), registered with
`schema%add_external_group("ocean_bc")` so the strict parser tolerates it —
migrating it is the known follow-up, see
[migrate_namelist_to_schema.md](migrate_namelist_to_schema.md).

## The seam

The default lives on the *type*; the schema registration captures
`default = current target value` at registration time and writes the parsed
value back through the pointer. So an absent namelist block is a no-op — that's
what gives you bit-identity for free.

```
type field (default)  ->  registered in the group's register_* routine
   ->  parsed + validated by nml_schema_t  ->  written back through the pointer
   ->  propagate cfg -> ocean_state (configure_ocean_*)
   ->  consume in kernel (gated)
```

## Recipe — a new `&ocean_hvisc_nml` key `foo`

1. **Field + default** on the right sub-type in `rdb_config.F90`
   (here inside `ocean_hvisc_config_t`): `real(wp) :: foo = 0.0_wp  !! <FORD doc>`.
   The default must reproduce current behaviour. **There is no `default=`
   argument on the registration** — the field initialiser *is* the default.
2. **Register it** in that group's `register_ocean_<group>` subroutine (here
   `register_ocean_hvisc`, called from `build_rdb_schema`):
   ```fortran
   pr => cfg%ocean%hvisc%foo
   call g%add(nml_real("foo", pr, "What foo does", units="m^2/s", min=0.0_wp))
   ```
   Use `nml_real` / `nml_int` / `nml_logical` / `nml_string` / `nml_enum` /
   `nml_real_array`. A string knob with a fixed accepted set **must** be
   `nml_enum` with the `allowed=` list — and that list must stay in sync with
   the matching `parse_*` routine in the physics module (the existing
   registrations say so in their FORD comments; `lateral_closure` ↔
   `parse_lateral_closure`, `eos` ↔ `parse_eos_variant`, …).
3. **Validate (optional).** Anything the engine can't express — exclusive
   bounds, cross-knob "requires X" / mutual-exclusion rules — goes in
   `validate_config` (`rdb_config.F90`) or in the group's `configure_*` routine.
   Set `has_error` / log and let the routine `error stop`. Never silently
   downgrade a bad combination.
4. **Propagate to state.** In `src/core/ocean/state/rdb_ocean_setup.F90`, in the
   relevant `configure_ocean_*` block: `ocean_state%hvisc%foo = cfg%ocean%hvisc%foo`.
   (Diagnostics are the exception — `configure_ocean_diag` lives in
   `src/driver/rdb_driver.F90`.)
5. **Consume in the kernel, gated on the knob.**
6. **FORD `!!` docstring on the field** — that *is* the documentation. Don't write
   parallel markdown describing the knob.
7. **Regenerate `docs/generated_nml_knobs.md`** — from a build directory,
   `./rdb_nml_doc ../docs/generated_nml_knobs.md`. The file is code-derived
   (defaults come straight from the field initialisers), so it never drifts as
   long as you regenerate after touching a group. Never hand-edit it.

## Default-off / bit-identity

The default must be the value that reproduces existing behaviour. Established
conventions in the tree: `0.0` (e.g. `nu_h`, `nu_4`), `.false.` (every `enable`
master switch), a negative sentinel ("unset, derive it" — `ah_bg = -1.0`,
`bkgnd_kd_min = -1.0`), a `1.0` multiplier, or `huge` as "no ceiling"
(`kv_max` / `kd_max`). The kernel is gated so that the default changes nothing —
verify by running `ctest` and confirming existing tests stay byte-for-byte.

## Adding the test row

Every test suite is one `|`-separated row in `RDB_TESTS`
(`tests/CMakeLists.txt`):

```
"<test-name>|<module-name>|<collect-fn>|<suite-display>|<regime>"
```

`<regime>` is the ctest LABEL and **must** be one of `ocean` (C-grid dyn-core +
sea-ice), `core` (shared infrastructure), or `ocean+core` (both labels).
Configure fails loud on anything else. Also append the module to
`tests/tester.F90` so the `fpm test` path stays in sync.

## The CLOSURE_MATRIX drift hook

`tools/check_closure_matrix.py` (pre-commit + CI) fires when you touch
`src/core/rdb_config.F90`, `docs/CLOSURE_MATRIX.md`, or `tests/CMakeLists.txt`.
It enforces:

1. **Knob coverage** — every auto-discovered `logical :: vmix_use_* =` toggle,
   plus a curated list (`ocean_hvisc_nml`, `ocean_bdrag_nml`,
   `ocean_kappa_shear_nml`) and the sea-ice selectors (`ocean_ice_nml`,
   `nk_ice`, `adv_substeps`, `evp_sub_steps`, `del_sh_min_scale`, `tdamp`), must
   appear as a string in `docs/CLOSURE_MATRIX.md`. Curated names must still
   exist in the config (rename guard).
2. **Test coverage** — every `test_*` token cited in the matrix must be
   registered in `tests/CMakeLists.txt`.
3. **The ALE-remap default** (`remap_method = "ppm"`) in the config must match
   the default the matrix states.

So: if you add a closure selector to the curated list, or touch one of the
guarded defaults, you **must** update `docs/CLOSURE_MATRIX.md` in the same
change or the hook fails the commit. A plain numeric knob that is not a closure
selector isn't enforced by the hook — but the matrix is still the single source
of truth for "how do I turn this on", so keep it current by convention.

## Tests

There isn't a dedicated "config" test per capability — the knob is exercised by
the feature's own analytical test (the kernel gated on it). What you *do*
verify: existing tests stay bit-identical with the default, and the feature test
flips the knob on. `tests/test_config.F90` covers the reader end to end and
`tests/test_config_schema.F90` covers the schema engine — extend the latter if
your knob introduces a new *shape* (first array knob in a group, a new enum
form, …).

## Gotchas

- **The field initialiser is the default.** Passing a "default" anywhere else
  (a local, the configure routine) means the parameter dump and the docs lie.
- **Enum registrations and `parse_*` routines must move together.** An allowed
  string with no `case` arm parses to the fallback and silently does the wrong
  thing; a `case` arm with no allowed string is rejected by the parser.
- **`&ocean_bc_nml` is the odd one out** — it is still a native reader, so a new
  key there means four edits (local declaration, the `namelist /ocean_bc_nml/`
  statement, seed from `cfg%`, copy back), not one registration line.
- **The closure-matrix hook can block a commit that compiles fine** — update
  `docs/CLOSURE_MATRIX.md` in the same change for a curated / vmix knob.
- **There is no C/Python FFI in this build.** Knobs reach the solver only
  through the namelist file — don't go looking for an API surface to wire them
  into.
- **`tools/nml_split.py`** migrates a legacy flat `&ocean_setup_nml` into the
  per-group `&ocean_<group>_nml` layout. You don't need it to *add* a knob; you'd
  touch its remap dicts only if you *rename* one (and then also the migration-ledger
  comment block in `rdb_config.F90`).

## Source pointers

- `src/core/rdb_config.F90` — `config_t`, `ocean_config_t` + sub-types,
  `read_config`, `build_rdb_schema` + the `register_*` routines,
  `read_ocean_bc_nml`, `validate_config`.
- `src/core/rdb_nml_schema.F90` — `nml_schema_t`, the `nml_*` constructors,
  `render_markdown`.
- `src/core/rdb_config_schema.F90` — the re-export shim (`build_rdb_schema`,
  `nml_dirname`).
- `src/core/ocean/state/rdb_ocean_setup.F90` — the `configure_ocean_*`
  propagation + cross-knob validation; `src/driver/rdb_driver.F90` —
  `configure_ocean_diag`.
- `app/nml_doc.F90` (`rdb_nml_doc`) → `docs/generated_nml_knobs.md`.
- `tools/check_closure_matrix.py`, `tools/nml_split.py`, `docs/CLOSURE_MATRIX.md`.
