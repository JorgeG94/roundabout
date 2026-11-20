"""Read a TOML config FILE into a typed `Config`, and merge configs.

The reader is the standard library's ``tomllib`` (Python >= 3.11) -- no
parser of our own, no dependency. What this module adds is only the
mapping from TOML tables to the typed groups, and it takes that mapping
from the generated classes (``Group._table``, emitted by
``tools/gen_python_config.py``), never from a second hand-kept list.

Values reach the `Config` through the SAME descriptor ``__set__`` used
everywhere else (`rdb._knob._KnobDescriptor`), so range/enum/type
validation is not re-implemented here: a malformed value still raises the
identical `ConfigParseError` shape. Key spelling is exact-case (the TOML
front end is), and a TOML integer is accepted for a real knob.
"""

from __future__ import annotations

import sys
from pathlib import Path

from ._errors import ConfigConflictError, ConfigParseError


def _name_to_attr():
    from ._config import ALL_GROUPS
    return {cls._group: attr for attr, cls in ALL_GROUPS.items()}


def _table_to_attr():
    from ._config import ALL_GROUPS
    return {cls._table: attr for attr, cls in ALL_GROUPS.items()}


def _loads(text: str) -> dict:
    try:
        import tomllib
    except ImportError as exc:  # pragma: no cover - python < 3.11
        raise ImportError(
            "reading a TOML file needs the standard-library `tomllib`, "
            "which Python 3.11+ ships; this interpreter is "
            + ".".join(map(str, sys.version_info[:3]))) from exc
    try:
        return tomllib.loads(text)
    except tomllib.TOMLDecodeError as exc:
        raise ConfigParseError(f"from_toml: invalid TOML: {exc}") from exc


def _flatten(doc: dict) -> dict:
    """``{"ocean": {"bt": {...}}, "time": {...}}`` -> ``{"ocean.bt":
    {...}, "time": {...}}``: one entry per leaf table. Only the leading
    ``ocean`` nests; any other table is a group of its own, and a key
    outside every table is refused by name."""
    out = {}
    for key, val in doc.items():
        if not isinstance(val, dict):
            raise ConfigParseError(
                f"from_toml: key {key!r} sits outside any table")
        if key == "ocean":
            for sub, body in val.items():
                if not isinstance(body, dict):
                    raise ConfigParseError(
                        f"from_toml: key 'ocean.{sub}' is not a table")
                out[f"ocean.{sub}"] = body
        else:
            out[key] = val
    return out


def config_from_toml_text(text: str):
    """Parse `text` into a fresh, typed `Config` -- every knob the file
    sets lands through the SAME descriptor `__set__` a Python script would
    use, so `Model.from_toml(path).config` is the "same typed config" a
    hand-written `Config()` would produce, not a second representation.
    """
    from ._config import Config

    table_to_attr = _table_to_attr()
    cfg = Config()
    for table, keys in _flatten(_loads(text)).items():
        attr = table_to_attr.get(table)
        if attr is None:
            raise ConfigParseError(
                f"from_toml: unknown table '[{table}]' -- not one of the "
                f"{len(table_to_attr)} tables this package's schema knows")
        inst = getattr(cfg, attr)
        cls = type(inst)
        for key, value in keys.items():
            desc = getattr(cls, key, None)
            if desc is None or not hasattr(desc, "kind"):
                raise ConfigParseError(
                    f"from_toml: unknown key {key!r} in [{table}]")
            if desc.kind == "real_array" and not isinstance(value, list):
                value = [value]
            setattr(inst, key, value)
    return cfg


#: Str knobs holding a filesystem path -- resolved against the config
#: file's own directory (no `os.chdir()`, ever). Named explicitly rather
#: than by a `_file` suffix heuristic so a future path-shaped knob must be
#: added here deliberately.
_PATH_KNOBS = {
    ("ocean_grid", "supergrid_file"),
    ("output", "bathymetry_file"),
    ("output", "restart_file"),
    ("ocean_bt", "wave_drag_file"),
    ("ocean_zinit", "file"),
    ("ocean_dataovr", "tau_x_file"),
    ("ocean_dataovr", "tau_y_file"),
    ("ocean_dataovr", "heat_file"),
    ("ocean_dataovr", "evap_file"),
    ("ocean_dataovr", "lprec_file"),
    ("ocean_dataovr", "salt_file"),
}


def resolve_relative_file_paths(cfg, base_dir):
    """Rewrite every explicitly-set path-shaped knob that is a RELATIVE
    path to be relative to `base_dir` (the config file's own directory)
    instead of the process cwd."""
    for group_attr, key in _PATH_KNOBS:
        inst = getattr(cfg, group_attr)
        if key not in inst._explicit_names():
            continue
        value = getattr(inst, key)
        if not value or Path(value).is_absolute():
            continue
        setattr(inst, key, str((Path(base_dir) / value).resolve()))
    return cfg


def merge_explicit(dst_cfg, src_cfg):
    """Merge every explicitly-set knob of `src_cfg` into `dst_cfg`.
    Raises `ConfigConflictError` (naming both values) if the same knob is
    explicit on both with DIFFERENT values -- a conflict between a file
    and a curated object raises rather than silently last-writer-wins.
    """
    from ._knob import MISSING

    name_to_attr = _name_to_attr()
    for group, key, value in src_cfg.explicit_knobs():
        attr = name_to_attr[group]
        dst_inst = getattr(dst_cfg, attr)
        current = getattr(dst_inst, key)
        if current is not MISSING and current != value:
            table = type(dst_inst)._table
            raise ConfigConflictError(
                f"from_toml: [{table}] {key} is set to {current!r} by the "
                f"config file AND to {value!r} by a curated keyword "
                f"argument -- refusing to silently pick one. Edit the "
                f"file, drop the curated override, or set "
                f"model.config.{attr}.{key} explicitly after loading.")
        setattr(dst_inst, key, value)
    return dst_cfg
