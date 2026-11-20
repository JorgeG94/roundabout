# Registering an input group in the strict schema (`rdb_nml_schema`)

The schema engine (`src/core/rdb_nml_schema.F90`, type `nml_schema_t`) is the
strict, self-documenting description of the input. One registration line per
knob carries the name, the pointer into the config, optional bounds/enum set,
units, and the doc string, and is the single list the TOML reader
(`src/core/config/rdb_toml_read.F90`), the parameter dumps
(`rdb_parameter_doc.{all,short}`, MOM6-style), and `render_markdown` all read
from. The input fails loudly (`error stop`, all errors collected, did-you-mean
suggestions) on unknown tables and keys and on type/range/enum violations. The
module keeps its historical `nml_*` names; a group `ocean_<x>` is the TOML table
`[ocean.<x>]`, any other group `<x>` the table `[<x>]`.

Every group is registered. The worked examples to copy are the pilot groups
`register_kappa_shear` and `register_epbl`, and — for a group with an array knob
and a big enum — `register_ocean_hvisc`. All the `register_*` routines are
declared as `module subroutine` interfaces in `rdb_config.F90` and implemented in
the `src/core/config/rdb_config_register_{core,vmix,lateral,forcing}.F90`
submodules, called from `build_rdb_schema`; `rdb_config_schema.F90` is only a
re-export shim, kept separate to avoid a circular dependency.

## Contract (load-bearing)

1. **Defaults come from the config-type field initializers, never from the
   registration.** Constructors capture `default = current target value` at
   registration time — so register after the `config_t` exists with its
   defaults, before parsing. There is no `default=` argument on purpose.
2. **Group names are registered bare** (`"ocean_epbl"`); the file spells
   the table `[ocean.epbl]`.
3. Validation that the engine can't express (exclusive bounds, cross-knob
   requirements involving non-config state) stays where it is
   (`validate_config` in `src/core/config/rdb_config_validate.F90`, and the `configure_ocean_*`
   routines in `src/core/ocean/state/rdb_ocean_setup.F90` and its
   `setup/rdb_ocean_setup_<group>.F90` submodules) — migrating a group
   must never weaken an existing check.
4. **A string knob with a fixed accepted set becomes `nml_enum`**, and its
   `allowed=` list must match the physics module's `parse_*` routine
   arm-for-arm. Extract the list from that routine; don't invent it.

## Recipe (one group per commit)

1. In `build_rdb_schema` (`src/core/config/rdb_config_register_core.F90`): add a
   `register_<group>` subroutine registering **every** field of its config type —
   `nml_real/int/logical/string/enum/real_array(name, cfg%...%field, doc [, units]
   [, min/max] [, allowed])` — and call it from `build_rdb_schema`.
2. Keep the config type + initializers in `rdb_config.F90`; the schema writes the
   parsed value back through the registered pointer, so there is no reader code
   to write.
3. Run the **full** suite (gfortran + GPU builds): every shipped `.toml` that
   reaches the driver is strictly validated for this group — fix stale keys in
   the files (or the registration, if the config type says the key is real).
4. Extend `tests/test_config_schema.F90` (and `tests/test_config_schema_toml.F90`,
   the engine's own behaviours) if the group introduces a new shape (first array
   knob, etc.). Their existing cases are the template.
5. **Regenerate `docs/generated_nml_knobs.md`** (below) — the new group appears
   there.

Parse + dumps are already wired at the production entry point (`app/main.F90`,
after the output-dir bootstrap, rank 0 only) — new groups need no wiring there.

## Regenerating the knob reference

`docs/generated_nml_knobs.md` is the full, code-derived table of every
input group + knob (name, default, units, description), rendered by the
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
