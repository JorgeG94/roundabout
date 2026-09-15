#!/usr/bin/env python3
"""Generate `python/rdb/_config_generated.py` from the live nml_schema_t.

Stdlib only (project rule: never `pip install`; no numpy anywhere in this
pipeline). Reads the JSON emitted by `rdb_nml_json` (built from
`app/nml_json.F90` -> `schema%render_json`, `src/core/rdb_nml_schema.F90`)
and writes ONE Python class per namelist group, with one descriptor
attribute per knob (see `python/rdb/_knob.py` for the descriptor
machinery: `Real`/`Int`/`Bool`/`Str`/`Enum`/`RealArray`).

This script parses no Fortran source -- the JSON is the only input, and it
was produced by introspecting the actual `nml_schema_t` object, not by
scraping text. That is the whole point of the pipeline (see
`tmp_local_artifacts/python_ffi_scope/06_python_surface_design.md` D5):
the generated Python layer cannot drift from the Fortran schema because it
is a direct rendering of it.

`&ocean_bc_nml` is emitted here like every other group (P4.5 migrated it
onto the schema; there are no more externally-registered groups).

Usage:

    ./build_shared/rdb_nml_json tmp_local_artifacts/schema.json
    python3 tools/gen_python_config.py tmp_local_artifacts/schema.json \\
            python/rdb/_config_generated.py

`python/tests/test_config_drift.py` runs exactly this pipeline into a temp
file and diffs it against the checked-in `_config_generated.py` -- that is
the anti-drift gate (there is no pre-commit hook for it, mirroring
`docs/generated_nml_knobs.md`'s situation before this phase; the ctest
suite is the gate here instead).
"""

from __future__ import annotations

import json
import keyword
import sys
from pathlib import Path

_KIND_TO_CLASS = {
    "real": "Real",
    "int": "Int",
    "logical": "Bool",
    "string": "Str",
    "enum": "Enum",
    "real_array": "RealArray",
}

_HEADER = '''"""Generated typed config layer -- DO NOT EDIT.

Emitted by `tools/gen_python_config.py` from the live `nml_schema_t`
(`src/core/rdb_config.F90` -> `build_rdb_schema`, walked via
`schema%render_json`, `src/core/rdb_nml_schema.F90`). One class per
namelist group; one descriptor attribute per knob, carrying the schema's
own doc/units/default/bounds/enum-set (see `python/rdb/_knob.py`).

ABSENCE IS THE DEFAULT, NOT THE FORTRAN VALUE: every knob starts unset,
and only explicitly-assigned knobs serialise via `to_namelist()`
(`python/rdb/_config.py`). This is what gives bit-identity to existing
namelists -- a config of all defaults serialises to nothing.

Regenerate with:

    ./build_shared/rdb_nml_json tmp_local_artifacts/schema.json
    python3 tools/gen_python_config.py tmp_local_artifacts/schema.json \\
            python/rdb/_config_generated.py

`python/tests/test_config_drift.py` is the anti-drift gate: it runs this
pipeline into a temp file and diffs it against this checked-in file,
failing with the regeneration command above if they differ.
"""

from __future__ import annotations

from ._knob import Bool, Enum, Group, Int, Real, RealArray, Str

'''


def pascal_case(name: str) -> str:
    return "".join(part.capitalize() for part in name.split("_") if part)


def py_attr(name: str) -> str:
    """A namelist key name is already a valid-ish Python identifier
    (letters/digits/underscore); guard only against a Python keyword
    collision (none currently exist in the schema, but a future knob
    could be named e.g. `class`)."""
    if keyword.iskeyword(name):
        return name + "_"
    return name


def emit_kwarg(name, value) -> str:
    return f"{name}={value!r}"


def emit_key(key: dict) -> list[str]:
    """Render one `attr = KindClass(...)` descriptor assignment as a list
    of source lines."""
    kind = key["kind"]
    cls = _KIND_TO_CLASS.get(kind)
    if cls is None:
        raise ValueError(f"unknown key kind {kind!r} for key {key['name']!r}")

    attr = py_attr(key["name"])
    args = [repr(key["name"])]
    kwargs = [
        emit_kwarg("doc", key["doc"]),
        emit_kwarg("units", key["units"]),
        emit_kwarg("required", key["required"]),
    ]

    if kind == "real":
        kwargs.append(emit_kwarg("default", key["default"]))
        if key["has_min"]:
            kwargs.append(emit_kwarg("has_min", True))
            kwargs.append(emit_kwarg("vmin", key["vmin"]))
        if key["has_max"]:
            kwargs.append(emit_kwarg("has_max", True))
            kwargs.append(emit_kwarg("vmax", key["vmax"]))
    elif kind == "int":
        kwargs.append(emit_kwarg("default", key["default"]))
        if key["has_min"]:
            kwargs.append(emit_kwarg("has_min", True))
            kwargs.append(emit_kwarg("vmin", key["vmin"]))
        if key["has_max"]:
            kwargs.append(emit_kwarg("has_max", True))
            kwargs.append(emit_kwarg("vmax", key["vmax"]))
    elif kind == "logical":
        kwargs.append(emit_kwarg("default", key["default"]))
    elif kind == "string":
        kwargs.append(emit_kwarg("default", key["default"]))
        kwargs.append(emit_kwarg("max_len", key["max_len"]))
    elif kind == "enum":
        kwargs.append(emit_kwarg("default", key["default"]))
        kwargs.append(emit_kwarg("allowed", tuple(key["allowed"])))
    elif kind == "real_array":
        kwargs.append(emit_kwarg("default", tuple(key["default"])))
        kwargs.append(emit_kwarg("size", key["size"]))

    dead = key.get("dead_on_ocean_path") or ""
    if dead:
        kwargs.append(emit_kwarg("dead_on_ocean_path", dead))

    all_args = ", ".join(args + kwargs)
    line = f"    {attr} = {cls}({all_args})"
    # Keep generated lines from becoming unreadably long; wrap at commas
    # when over ~95 columns. Simple, deterministic, no external formatter
    # (stdlib-only rule).
    if len(line) <= 95:
        return [line]
    wrapped = [f"    {attr} = {cls}("]
    for i, a in enumerate(args + kwargs):
        comma = "," if i < len(args) + len(kwargs) - 1 else ","
        wrapped.append(f"        {a}{comma}")
    wrapped.append("    )")
    return wrapped


def emit_group(group: dict) -> tuple[str, list[str]]:
    """Return (attr_name, source_lines) for one group class."""
    gname = group["name"]
    cls_name = pascal_case(gname)
    lines = [f"class {cls_name}(Group):"]
    doc = group["doc"] or f"`&{gname}_nml`."
    lines.append(f'    """`&{gname}_nml` -- {doc}"""')
    lines.append("")
    lines.append(f"    _nml_name = {gname!r}")
    keys = group.get("keys", [])
    if not keys:
        lines.append("")
    for key in keys:
        lines.append("")
        lines.extend(emit_key(key))
    lines.append("")
    return gname, lines


def generate(schema: dict) -> str:
    groups = schema["groups"]
    out_lines = [_HEADER.rstrip("\n")]
    group_attrs = []
    for group in groups:
        attr, lines = emit_group(group)
        group_attrs.append((attr, pascal_case(attr)))
        out_lines.append("")
        out_lines.extend(lines)

    out_lines.append("")
    out_lines.append(
        "#: schema group name -> generated Group subclass, in the order "
        "build_rdb_schema")
    out_lines.append(
        "#: registers them. Combined with OceanBc (python/rdb/_config_bc.py)")
    out_lines.append("#: by python/rdb/_config.py to build the full Config object.")
    out_lines.append("GENERATED_GROUPS = {")
    for attr, cls in group_attrs:
        out_lines.append(f"    {attr!r}: {cls},")
    out_lines.append("}")
    out_lines.append("")

    ngroups = len(groups)
    nknobs = sum(len(g.get("keys", [])) for g in groups)
    out_lines.append(f"N_GROUPS = {ngroups}")
    out_lines.append(f"N_KNOBS = {nknobs}")

    # Exactly one trailing newline, no blank line before EOF -- matches
    # what pre-commit's end-of-file-fixer enforces on the checked-in file,
    # so a regeneration never drifts on whitespace alone.
    return "\n".join(out_lines) + "\n"


def main(argv: list[str]) -> int:
    if len(argv) != 3:
        sys.stderr.write(
            "usage: gen_python_config.py <schema.json> <_config_generated.py>\n")
        return 2
    schema_path = Path(argv[1])
    out_path = Path(argv[2])

    with schema_path.open("r", encoding="utf-8") as f:
        schema = json.load(f)

    text = generate(schema)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text, encoding="utf-8")
    ngroups = len(schema["groups"])
    nknobs = sum(len(g.get("keys", [])) for g in schema["groups"])
    sys.stderr.write(
        f"wrote {out_path} ({ngroups} groups, {nknobs} knobs)\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
