"""Render and line-patch roundabout TOML inputs (stdlib only, no `tomllib`).

The regression drivers build and edit model inputs as plain text.  Some of
them run under python3.10 (the interpreter that has netCDF4), which has no
`tomllib`, so this module parses nothing: it WRITES TOML and PATCHES the
line-per-key files the repo ships (see `tools/nml_to_toml.py` for the shape).

Groups are named the way the schema and the drivers always have, with the
`_nml` suffix (`time_nml`, `ocean_bt_nml`); the table path follows the one
mapping rule (`group_to_table`, mirroring the Fortran function of the same
name): only the LEADING `ocean_` nests.

    ocean_bt_nml  ->  [ocean.bt]        time_nml  ->  [time]

Values are Python values: bool -> true/false, int, float, str -> basic
string, list/tuple -> array.

Patching works on a file written one key per line (arrays on one line); it
raises on a multi-line array rather than guess.
"""

import re

OCEAN_PREFIX = "ocean_"

_HEADER = re.compile(r"^\s*\[\s*([A-Za-z0-9_.\-]+)\s*\]\s*(?:#.*)?$")
_KEY = re.compile(r"^\s*([A-Za-z0-9_\-]+)\s*=")


def group_to_table(group):
    """`ocean_hvisc_nml`/`ocean_hvisc` -> `ocean.hvisc`; `time_nml` -> `time`."""
    g = group[:-4] if group.endswith("_nml") else group
    if len(g) > len(OCEAN_PREFIX) and g.startswith(OCEAN_PREFIX):
        return "ocean." + g[len(OCEAN_PREFIX):]
    return g


def fmt_value(v):
    """TOML text of a Python value."""
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        if v != v:
            return "nan"
        if v in (float("inf"), float("-inf")):
            return "inf" if v > 0 else "-inf"
        return repr(v)
    if isinstance(v, str):
        esc = v.replace("\\", "\\\\").replace('"', '\\"')
        esc = esc.replace("\n", "\\n").replace("\t", "\\t")
        return '"' + esc + '"'
    if isinstance(v, (list, tuple)):
        return "[" + ", ".join(fmt_value(x) for x in v) + "]"
    raise TypeError("cannot render {!r}".format(v))


def render(groups, header=""):
    """`{group: {key: value}}` -> TOML text, groups and keys sorted."""
    out = []
    for line in header.splitlines():
        out.append("# " + line if line else "#")
    tables = sorted((group_to_table(g), g) for g in groups)
    for table, g in tables:
        out.append("[{}]".format(table))
        for k in sorted(groups[g]):
            out.append("{} = {}".format(k, fmt_value(groups[g][k])))
        out.append("")
    return "\n".join(out)


def _strip_comment(line):
    """The line without its trailing `# comment` (strings respected)."""
    quote = None
    i = 0
    while i < len(line):
        c = line[i]
        if quote:
            if c == "\\" and quote == '"':
                i += 1
            elif c == quote:
                quote = None
        elif c in "'\"":
            quote = c
        elif c == "#":
            return line[:i]
        i += 1
    return line


def _scan(lines):
    """Per line: (table or None, key or None).  `table` is the table the
    line belongs to; `key` is set on assignment lines only."""
    table = None
    info = []
    for n, line in enumerate(lines):
        code = _strip_comment(line)
        m = _HEADER.match(line)
        if m:
            table = m.group(1)
            info.append((table, None))
            continue
        k = _KEY.match(code)
        if k:
            rhs = code[k.end():]
            if rhs.count("[") != rhs.count("]"):
                raise ValueError("line {}: multi-line array is not patchable"
                                 .format(n + 1))
            info.append((table, k.group(1)))
        else:
            info.append((table, None))
    return info


def get(text, group, key):
    """Raw value text of `key` in `group`, or None."""
    table = group_to_table(group)
    lines = text.splitlines()
    for line, (t, k) in zip(lines, _scan(lines)):
        if t == table and k == key:
            code = _strip_comment(line)
            return code[_KEY.match(code).end():].strip()
    return None


def get_float(text, group, key):
    raw = get(text, group, key)
    return None if raw is None else float(raw)


def has_group(text, group):
    table = group_to_table(group)
    return any(_HEADER.match(l) and _HEADER.match(l).group(1) == table
               for l in text.splitlines())


def patch(text, overrides):
    """Apply `{group: {key: value}}` to `text`.

    A key already in its table is replaced in place (its trailing comment is
    dropped); a key missing from an existing table is inserted after the
    table's last assignment; a table the file lacks is appended.
    """
    lines = text.splitlines()
    for group in sorted(overrides):
        table = group_to_table(group)
        kv = overrides[group]
        info = _scan(lines)
        present = {k for (t, k) in info if t == table and k}
        in_table = [i for i, (t, _k) in enumerate(info) if t == table]
        has_header = any(_HEADER.match(lines[i]) for i in in_table)
        for i in in_table:
            k = info[i][1]
            if k in kv:
                lines[i] = "{} = {}".format(k, fmt_value(kv[k]))
        new = [k for k in sorted(kv) if k not in present]
        if not new:
            continue
        add = ["{} = {}".format(k, fmt_value(kv[k])) for k in new]
        if has_header:
            last = max(i for i in in_table if info[i][1] is not None
                       or _HEADER.match(lines[i]))
            lines[last + 1:last + 1] = add
        else:
            while lines and not lines[-1].strip():
                lines.pop()
            lines += ([""] if lines else []) + ["[{}]".format(table)] + add
    return "\n".join(lines) + "\n"
