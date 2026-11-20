#!/usr/bin/env python3
"""Convert a roundabout namelist (.nml) to the TOML input format, text-level.

TEMPORARY: this tool exists for the namelist -> TOML migration and is deleted
at the end of it.  Stdlib only.

    tools/nml_to_toml.py --nml-json build_gfortran/rdb_nml_json in.nml [out.toml]
    tools/nml_to_toml.py --schema schema.json in.nml [out.toml]

The conversion keeps the file's shape so the result reads as hand-written:
comments (`!` -> `#`), blank lines, key order and the column of trailing
comments survive; one key per line; groups become tables.

Mapping (mirrors `group_to_table` in `src/core/config/rdb_toml_read.F90`):

    &ocean_<x>_nml  ->  [ocean.<x>]      (only the LEADING `ocean_` nests)
    &<x>_nml        ->  [<x>]

Keys are emitted in the REGISTERED spelling taken from the schema JSON
(`rdb_nml_json`), because the namelist reader is case-insensitive and TOML is
not.  Values are typed from the schema: logical -> true/false, string/enum ->
basic string, real_array -> `[a, b]` (always an array, even for one element),
real/int -> a TOML number resolving to the same value (original digits kept
where already valid TOML; `1.`, `.5`, `1d3`, `1.e5` are rewritten).

Anything the converter cannot translate faithfully (an unknown group or key,
a repeat count `3*1.0`, a kind suffix, a continuation line) raises an error;
nothing is guessed.  Repeated keys inside one group keep the LAST value, as
the namelist reader does.
"""

import argparse
import json
import re
import subprocess
import sys
import tempfile

OCEAN_PREFIX = "ocean_"


class ConvertError(Exception):
    pass


def group_to_table(group):
    """Table path of a bare group name (`ocean_hvisc` -> `ocean.hvisc`)."""
    if len(group) > len(OCEAN_PREFIX) and group.startswith(OCEAN_PREFIX):
        return "ocean." + group[len(OCEAN_PREFIX):]
    return group


# --------------------------------------------------------------------------
# Schema
# --------------------------------------------------------------------------
class Schema:
    def __init__(self, doc):
        # lower-cased bare group -> (canonical group, {lower key: (name, kind)})
        self.groups = {}
        for g in doc["groups"]:
            keys = {k["name"].lower(): (k["name"], k["kind"]) for k in g["keys"]}
            self.groups[g["name"].lower()] = (g["name"], keys)

    @classmethod
    def from_file(cls, path):
        with open(path) as fh:
            return cls(json.load(fh))

    @classmethod
    def from_binary(cls, binary):
        with tempfile.TemporaryDirectory() as td:
            out = td + "/schema.json"
            subprocess.run([binary, out], check=True, capture_output=True)
            return cls.from_file(out)

    def lookup_group(self, name):
        bare = name.lower()
        if bare.endswith("_nml"):
            bare = bare[:-4]
        if bare not in self.groups:
            raise ConvertError("unknown group &{}".format(name))
        return self.groups[bare]


# --------------------------------------------------------------------------
# Lexing helpers
# --------------------------------------------------------------------------
_TOKEN = re.compile(
    r"""\s*(?:
        (?P<str>'(?:[^']|'')*'|"(?:[^"]|"")*")
      | (?P<comma>,)
      | (?P<eq>=)
      | (?P<slash>/)
      | (?P<word>[^\s,='"/]+)
    )""", re.X)


def split_comment(line):
    """(code, comment) with `comment` including no leading `!`, or None."""
    quote = None
    i = 0
    while i < len(line):
        c = line[i]
        if quote:
            if c == quote:
                if i + 1 < len(line) and line[i + 1] == quote:
                    i += 1  # doubled quote escape
                else:
                    quote = None
        elif c in "'\"":
            quote = c
        elif c == "!":
            return line[:i], line[i + 1:]
        i += 1
    return line, None


def tokenize(code):
    toks = []
    pos = 0
    code = code.rstrip()
    while pos < len(code):
        m = _TOKEN.match(code, pos)
        if not m:
            raise ConvertError("cannot tokenize {!r}".format(code[pos:]))
        kind = m.lastgroup
        toks.append((kind, m.group(kind), m.start(kind), m.end(kind)))
        pos = m.end()
    return toks


# --------------------------------------------------------------------------
# Value formatting
# --------------------------------------------------------------------------
_NUM = re.compile(r"^([+-]?)(\d*)(?:(\.)(\d*))?(?:[eEdD]([+-]?\d+))?$")


def fmt_number(tok, want_int):
    m = _NUM.match(tok)
    if not m or not (m.group(2) or m.group(4)):
        raise ConvertError("not a number: {!r}".format(tok))
    sign, ip, dot, fp, ex = m.groups()
    ip = ip or ""
    fp = fp or ""
    if want_int:
        if dot or ex:
            raise ConvertError("integer key given {!r}".format(tok))
        return sign.replace("+", "") + (ip.lstrip("0") or "0")
    # TOML: no leading zeros in the integer part, a digit on both sides of
    # the point.  Keep the digits otherwise.
    ip = ip.lstrip("0") or "0"
    if not dot and not ex:
        return sign.replace("+", "") + ip  # a bare integer for a real key
    out = sign.replace("+", "") + ip
    if dot:
        out += "." + (fp or "0")
    if ex is not None:
        out += "e" + ex
    return out


def fmt_bool(tok):
    t = tok.lower().strip(".")
    if t in ("t", "true"):
        return "true"
    if t in ("f", "false"):
        return "false"
    raise ConvertError("not a logical: {!r}".format(tok))


def fmt_string(tok):
    if tok[0] in "'\"":
        q = tok[0]
        body = tok[1:-1].replace(q + q, q)
    else:
        raise ConvertError("unquoted string value {!r}".format(tok))
    esc = body.replace("\\", "\\\\").replace('"', '\\"')
    return '"' + esc + '"'


_PLACEHOLDER = re.compile(r"^\{\{\w+\}\}$")


def fmt_values(kind, vals):
    """vals: list of (tokkind, text).  Returns the TOML value text.

    A bare `{{NAME}}` token (a text-template placeholder, as in
    `tests/regression/vcoord_templates/`) passes through verbatim.
    """
    if len(vals) == 1 and vals[0][0] == "word" and _PLACEHOLDER.match(vals[0][1]):
        return vals[0][1]
    if kind == "real_array":
        return "[" + ", ".join(fmt_number(t, False) for k, t in vals) + "]"
    if len(vals) != 1:
        raise ConvertError("{} key given {} values".format(kind, len(vals)))
    k, t = vals[0]
    if "*" in t and k == "word":
        raise ConvertError("repeat count not supported: {!r}".format(t))
    if kind == "logical":
        return fmt_bool(t)
    if kind in ("string", "enum"):
        return fmt_string(t) if k == "str" else fmt_string(t)
    if kind in ("real", "int"):
        if k != "word":
            raise ConvertError("number expected, got {!r}".format(t))
        return fmt_number(t, kind == "int")
    raise ConvertError("unknown key kind {!r}".format(kind))


# --------------------------------------------------------------------------
# Conversion
# --------------------------------------------------------------------------
def assignments(toks):
    """Split a token list into ([(key, [values...])], ended_by_slash)."""
    out = []
    ended = False
    i = 0
    n = len(toks)
    cur = None
    while i < n:
        kind, text, _s, _e = toks[i]
        if kind == "slash":
            ended = True
            if i != n - 1:
                raise ConvertError("text after the group terminator `/`")
            break
        if kind == "word" and i + 1 < n and toks[i + 1][0] == "eq":
            cur = (text, toks[i], [])
            out.append(cur)
            i += 2
            continue
        if cur is None:
            raise ConvertError("value without a key: {!r}".format(text))
        if kind == "comma":
            i += 1
            continue
        if kind in ("word", "str"):
            cur[2].append((kind, text))
            i += 1
            continue
        raise ConvertError("unexpected token {!r}".format(text))
    return out, ended


def group_dedents(lines):
    """Per line: the indent of the first assignment of the group the line is
    in (0 outside a group, or when the header line carries the keys)."""
    res = [0] * len(lines)
    i = 0
    while i < len(lines):
        code = split_comment(lines[i])[0].strip()
        if not code.startswith("&"):
            i += 1
            continue
        j, first = i, None
        if "=" in code:
            first = 0
        while True:
            c = split_comment(lines[j])[0]
            if first is None and "=" in c:
                first = len(c) - len(c.lstrip())
            if c.strip().endswith("/") or j + 1 >= len(lines):
                break
            j += 1
        for k in range(i, j + 1):
            res[k] = first or 0
        i = j + 1
    return res


def convert(text, schema):
    out = []
    group = None          # (canonical bare group, keymap) or None
    seen_tables = set()
    seen_keys = set()

    all_lines = text.splitlines()
    dedents = group_dedents(all_lines)
    for ln, raw in enumerate(all_lines, 1):
        try:
            dedent = dedents[ln - 1]
            code, comment = split_comment(raw)
            stripped = code.strip()
            indent = len(code) - len(code.lstrip())
            cut = min(indent, dedent) if group else 0
            if not stripped:
                out.append("" if comment is None else
                           " " * (indent - cut) + "#" + comment)
                continue
            if _PLACEHOLDER.match(stripped):
                out.append(" " * (indent - cut) + stripped)
                continue
            if stripped.startswith("&") and group is None:
                toks = tokenize(stripped[1:])
                if not toks or toks[0][0] != "word":
                    raise ConvertError("bad group header")
                canon, keys = schema.lookup_group(toks[0][1])
                table = group_to_table(canon)
                if table in seen_tables:
                    raise ConvertError("group &{} appears twice".format(canon))
                seen_tables.add(table)
                seen_keys = set()
                group = (canon, keys)
                rest = toks[1:]
                if not rest:
                    out.append("[{}]".format(table) +
                               ("" if comment is None else "  #" + comment))
                    continue
                out.append("[{}]".format(table))
                body, ended = assignments(rest)
                emit_keys(out, group, body, comment, None, "", 0, 0, seen_keys)
                if ended:
                    group = None
                continue
            if group is None:
                raise ConvertError("assignment outside a group")
            if stripped == "/":
                group = None
                if comment is not None:
                    out.append("#" + comment)
                continue
            toks = tokenize(stripped)
            body, ended = assignments(toks)
            cut = min(indent, dedent)
            emit_keys(out, group, body, comment, toks, code, indent, cut, seen_keys)
            if ended:
                group = None
        except ConvertError as exc:
            raise ConvertError("line {}: {} [{}]".format(ln, exc, raw.strip()))
    if group is not None:
        raise ConvertError("unterminated group")
    res = []
    for line in out:  # the removed `/` lines leave doubled blank lines
        if line == "" and res and res[-1] == "":
            continue
        res.append(line)
    while res and res[-1] == "":
        res.pop()
    return "\n".join(res) + "\n"


def emit_keys(out, group, body, comment, toks, code, indent, cut, seen):
    """Append one TOML line per assignment; the comment rides the last one."""
    canon, keys = group
    for idx, (key, keytok, vals) in enumerate(body):
        low = key.lower()
        if low not in keys:
            raise ConvertError("unknown key {} in &{}".format(key, canon))
        if low in seen:
            raise ConvertError("key {} repeated in &{}".format(key, canon))
        seen.add(low)
        name, kind = keys[low]
        if not vals:
            raise ConvertError("key {} has no value".format(key))
        value = fmt_values(kind, vals)
        if len(body) == 1 and toks is not None:
            # keep the original gaps around `=` so the columns still line up
            c = code.strip()
            eq = toks[1]
            between = c[keytok[3]:eq[2]]
            after = c[eq[3]:]
            gap = after[:len(after) - len(after.lstrip())] or " "
            pad = max(1, len(between) + len(key) - len(name)) if between else 1
            line = "{}{}={}{}".format(name, " " * pad, gap, value)
        else:
            line = "{} = {}".format(name, value)
        if idx == len(body) - 1 and comment is not None:
            trailing = len(code) - len(code.rstrip())
            line += " " * max(2, trailing - cut) + "#" + comment
        out.append(line)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--schema", help="schema JSON written by rdb_nml_json")
    src.add_argument("--nml-json", help="path to the rdb_nml_json executable")
    ap.add_argument("input")
    ap.add_argument("output", nargs="?", help="default: stdout")
    args = ap.parse_args(argv)
    schema = Schema.from_file(args.schema) if args.schema \
        else Schema.from_binary(args.nml_json)
    with open(args.input) as fh:
        text = fh.read()
    try:
        res = convert(text, schema)
    except ConvertError as exc:
        sys.exit("{}: {}".format(args.input, exc))
    if args.output:
        with open(args.output, "w") as fh:
            fh.write(res)
    else:
        sys.stdout.write(res)
    return 0


if __name__ == "__main__":
    sys.exit(main())
