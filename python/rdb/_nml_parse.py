"""P6 -- a small, deliberately bounded Fortran-namelist READER, used only by
`Model.from_namelist`.

This is NOT a general Fortran namelist parser (no `$`-continuation, no
repeated/merged groups, no `arr(2:4) = ...` index slicing, no `n*value`
repeat-count shorthand). It handles exactly the shape both `Config.
to_namelist()` emits and every namelist in this repo's `validation_examples/`
+ test fixtures uses: `&group_nml key = value, key2 = value2 /` blocks,
one or more per file, `!` end-of-line comments, quoted strings (`'` or
`"`), `.true.`/`.false.`, and comma-separated lists for array knobs.

Values are handed back as PLAIN PYTHON OBJECTS (str/bool/int/float or a
list of numbers) and assigned onto a `Config` through the SAME descriptor
`__set__` used everywhere else (`rdb._knob._KnobDescriptor`), so range/
enum/type validation is not re-implemented here -- a malformed value still
raises the identical `ConfigParseError` shape.
"""

from __future__ import annotations

import re
from pathlib import Path

from ._errors import ConfigConflictError, ConfigParseError

_GROUP_START_RE = re.compile(r"&([A-Za-z_]\w*)")
_KEY_RE = re.compile(r"([A-Za-z_]\w*)\s*=")


def _strip_comments(text: str) -> str:
    """Remove `!...` to end-of-line, respecting quotes."""
    out = []
    in_quote = None
    i = 0
    n = len(text)
    while i < n:
        c = text[i]
        if in_quote:
            out.append(c)
            if c == in_quote:
                in_quote = None
        elif c in ("'", '"'):
            in_quote = c
            out.append(c)
        elif c == "!":
            j = text.find("\n", i)
            if j == -1:
                break
            i = j
            continue
        else:
            out.append(c)
        i += 1
    return "".join(out)


def _mask_quotes(text: str) -> str:
    """Same length as `text`; quoted spans replaced with 'Q' so a `=`/`,`
    inside a quoted string never confuses the key/value splitter."""
    out = list(text)
    in_quote = None
    for i, c in enumerate(text):
        if in_quote:
            out[i] = "Q" if c != in_quote else c
            if c == in_quote:
                in_quote = None
        elif c in ("'", '"'):
            in_quote = c
    return "".join(out)


def _split_top_commas(text: str) -> list:
    masked = _mask_quotes(text)
    parts = []
    depth = 0
    start = 0
    for i, c in enumerate(masked):
        if c == "(":
            depth += 1
        elif c == ")":
            depth -= 1
        elif c == "," and depth == 0:
            parts.append(text[start:i])
            start = i + 1
    parts.append(text[start:])
    return [p.strip() for p in parts if p.strip() != ""]


def _literal(token: str):
    token = token.strip()
    if len(token) >= 2 and token[0] == token[-1] and token[0] in ("'", '"'):
        return token[1:-1]
    low = token.lower()
    if low in (".true.", "t", ".t."):
        return True
    if low in (".false.", "f", ".f."):
        return False
    try:
        return int(token)
    except ValueError:
        pass
    try:
        return float(token.replace("d", "e").replace("D", "e"))
    except ValueError:
        raise ConfigParseError(f"cannot parse namelist literal {token!r}")


def _parse_group_body(body: str) -> dict:
    body = _strip_comments(body)
    masked = _mask_quotes(body)
    matches = list(_KEY_RE.finditer(masked))
    out = {}
    for idx, m in enumerate(matches):
        key = m.group(1)
        val_start = m.end()
        val_end = matches[idx + 1].start() if idx + 1 < len(matches) \
            else len(body)
        raw = body[val_start:val_end].strip()
        if raw.endswith(","):
            raw = raw[:-1].strip()
        tokens = _split_top_commas(raw)
        if not tokens:
            continue
        values = [_literal(t) for t in tokens]
        out[key] = values[0] if len(values) == 1 else values
    return out


def parse_namelist_text(text: str) -> dict:
    """``&group_nml key=val, key2=val2 / ...`` -> ``{group_name:
    {key: value, ...}, ...}``, group name WITHOUT the leading `&` or the
    trailing `_nml`. Later blocks for the same group MERGE (later keys
    win) rather than replace -- matches how a repeated namelist block
    behaves under Fortran's own reader."""
    out = {}
    pos = 0
    n = len(text)
    while True:
        m = _GROUP_START_RE.search(text, pos)
        if m is None:
            break
        name = m.group(1)
        if name.lower().endswith("_nml"):
            name = name[:-4]
        # Scan forward, QUOTE-AWARE, for the terminating `/` -- a naive
        # regex up to the first `/` breaks on any quoted string value
        # containing one (e.g. bathymetry_file = "sub/bathy.nc").
        body_start = m.end()
        i = body_start
        in_quote = None
        while i < n:
            c = text[i]
            if in_quote:
                if c == in_quote:
                    in_quote = None
            elif c == "!":
                # Skip the whole comment. MUST come before the quote check:
                # this scanner runs on RAW text, so an ordinary English
                # apostrophe inside a comment ("the setter's x/y lengths")
                # would otherwise open a phantom quote and swallow every
                # `/` group terminator until the next apostrophe. That made
                # Model.from_namelist() fail on real files in the tree,
                # including acc_channel.nml and double_gyre_mom6.nml.
                while i < n and text[i] != "\n":
                    i += 1
                continue
            elif c in ("'", '"'):
                in_quote = c
            elif c == "/":
                break
            i += 1
        body = text[body_start:i]
        pos = i + 1
        parsed = _parse_group_body(body)
        out.setdefault(name.lower(), {}).update(parsed)
    return out


def _name_to_attr():
    from ._config import ALL_GROUPS
    return {cls._nml_name.lower(): attr for attr, cls in ALL_GROUPS.items()}


def config_from_namelist_text(text: str):
    """Parse `text` into a fresh, typed `Config` -- every knob the file
    sets lands through the SAME descriptor `__set__` a Python script would
    use, so `Model.from_namelist(path).config` is the "same typed config"
    a hand-written `Config()` would produce (D7 P6 / D2's migration
    story), not a second representation.
    """
    from ._config import Config

    name_to_attr = _name_to_attr()
    parsed = parse_namelist_text(text)
    cfg = Config()
    for group_name, keys in parsed.items():
        attr = name_to_attr.get(group_name)
        if attr is None:
            raise ConfigParseError(
                f"from_namelist: unknown namelist group '&{group_name}_nml' "
                f"-- not one of the {len(name_to_attr)} groups this "
                f"package's schema knows")
        inst = getattr(cfg, attr)
        cls = type(inst)
        for key, value in keys.items():
            desc = getattr(cls, key, None)
            if desc is None or not hasattr(desc, "kind"):
                raise ConfigParseError(
                    f"from_namelist: unknown key {key!r} in "
                    f"&{group_name}_nml")
            if desc.kind == "real_array" and not isinstance(value, list):
                value = [value]
            setattr(inst, key, value)
    return cfg


#: Str knobs holding a filesystem path -- resolved against the namelist's
#: own directory (D2's `Model.from_namelist` migration story: no
#: `os.chdir()`, ever). Named explicitly rather than by a `_file` suffix
#: heuristic so a future path-shaped knob must be added here deliberately.
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
    path to be relative to `base_dir` (the namelist file's own directory)
    instead of the process cwd. Retires the `os.chdir()` the recovered
    wrapper did around every file create."""
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
    explicit on both with DIFFERENT values -- "a conflict between a file
    and a curated object raises rather than silently last-writer-wins."
    """
    from ._knob import MISSING

    name_to_attr = _name_to_attr()
    for group_nml, key, value in src_cfg.explicit_knobs():
        attr = name_to_attr[group_nml]
        dst_inst = getattr(dst_cfg, attr)
        current = getattr(dst_inst, key)
        if current is not MISSING and current != value:
            raise ConfigConflictError(
                f"from_namelist: &{group_nml}_nml {key} is set to "
                f"{current!r} by the namelist file AND to {value!r} by a "
                f"curated keyword argument -- refusing to silently pick "
                f"one. Edit the file, drop the curated override, or set "
                f"model.config.{attr}.{key} explicitly after loading.")
        setattr(dst_inst, key, value)
    return dst_cfg
