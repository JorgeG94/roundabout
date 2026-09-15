# Migrating a namelist group to the strict schema (`rdb_nml_schema`)

The schema engine (`src/core/rdb_nml_schema.F90`, type `nml_schema_t`) is the
strict, self-documenting replacement for the soft native `read(nml=)` readers
in `rdb_config.F90`. One registration line per knob carries the name, the
pointer into the config, optional bounds/enum set, units, and the doc string —
and is the single list the parser, the parameter dumps
(`rdb_parameter_doc.{all,short}`, MOM6-style), and `render_markdown` all
read from. Migrated groups fail loudly (`error stop`, all errors collected,
did-you-mean suggestions) on unknown keys, type/range/enum violations, and
unsupported syntax. Destined for `pic` once trialed; the engine depends only
on `iso_fortran_env` + `pic_logger`.

**Status: every group is migrated except `&ocean_bc_nml`.** That one still has a
native reader (`read_ocean_bc_nml` in `rdb_config.F90`) and is registered with
`schema%add_external_group("ocean_bc")` so the strict parser accepts the block
without validating it. Migrating it is the outstanding job this guide is for.
The worked examples to copy are the pilot groups `register_kappa_shear` and
`register_epbl`, and — for a group with an array knob and a big enum —
`register_ocean_hvisc`. All the `register_*` routines live in `rdb_config.F90`
and are called from `build_rdb_schema`; `rdb_config_schema.F90` is only a
re-export shim (it also owns `nml_dirname`), kept separate to avoid a circular
dependency.

## Contract (load-bearing)

1. **Defaults come from the config-type field initializers, never from the
   registration.** Constructors capture `default = current target value` at
   registration time — so register after the `config_t` exists with its
   defaults, before parsing. There is no `default=` argument on purpose.
2. **Group names are registered bare** (`"ocean_epbl"`); files spell
   `&ocean_epbl_nml`.
3. Validation that the engine can't express (exclusive bounds, cross-knob
   requirements involving non-config state) stays where it is
   (`validate_config` in `rdb_config.F90`, and the `configure_ocean_*`
   routines in `src/core/ocean/state/rdb_ocean_setup.F90`) — migrating a group
   must never weaken an existing check.
4. **A string knob with a fixed accepted set becomes `nml_enum`**, and its
   `allowed=` list must match the physics module's `parse_*` routine
   arm-for-arm. Extract the list from that routine; don't invent it.

## Recipe (one group per commit)

1. In `build_rdb_schema` (`rdb_config.F90`): remove the group's
   `add_external_group` line; add a `register_<group>` subroutine registering
   **every** field of its config type — `nml_real/int/logical/string/enum/real_array(name,
   cfg%...%field, doc [, units] [, min/max] [, allowed])` — and call it from
   `build_rdb_schema`.
2. In `rdb_config.F90`: delete the group's `read_*_nml` subroutine, its
   `namelist /…/` statement, shadow locals, and the call site in
   `read_config`. Keep the config type + initializers.
3. Run the **full** suite (gfortran + GPU builds): every repo `.nml` that
   reaches the driver is now strictly validated for this group — fix stale
   keys in the files (or the registration, if the config type says the key
   is real).
4. Extend `tests/test_config_schema.F90` if the group introduces a new
   shape (first array knob, etc.). Its existing cases are the template:
   `happy_path_kshear_epbl`, `typo_key_suggestion`, `out_of_range`,
   `doc_short_non_default`, `legacy_group_external`.
5. **Regenerate `docs/generated_nml_knobs.md`** (below) — a newly migrated group
   appears there for the first time, since external groups contribute nothing to
   the render.

Parse + dumps are already wired at the production entry point (`app/main.F90`,
after the output-dir bootstrap, rank 0 only) — new groups need no wiring there.

## Regenerating the knob reference

`docs/generated_nml_knobs.md` is the full, code-derived table of every
namelist group + knob (name, default, units, description), rendered by the
`rdb_nml_doc` CLI tool (`app/nml_doc.F90`). It builds a default
`config_t`, registers the schema, and calls `render_markdown`. Because the
defaults come straight from the `config_t` field initialisers, the file
never drifts from the code as long as you regenerate it after touching a
group. From a build directory:

```
./rdb_nml_doc ../docs/generated_nml_knobs.md   # arg 1 = output path (default ./nml_knobs.md)
```

The tool is host-only (no NetCDF / GPU) and is built alongside `rdb`.
Never hand-edit the generated file.
