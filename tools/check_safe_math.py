#!/usr/bin/env python3
"""Keep bare transcendentals out of the model: they go through `rdb_safe_math`.

The intrinsic `exp(x)` is not a function of `x` alone: the compiler may route it
through a vector math library (glibc `libmvec`, NVHPC `libpgmath`) in a loop's
vector body and through scalar libm in its remainder, so the same cell rounds
differently on a different tile width or rank count.  `rdb_safe_math` is written
in IEEE `+ - * /` and `sqrt` only and is position independent (see the module
header, `src/framework/rdb_safe_math.F90`, and `plans/multirank/SAFE_MATH_PLAN.md`).

Flagged in `src/**` (`.F90`, `.f90`, `.inc`) and `app/**`:

  * a call of  exp log log10 sin cos tan asin acos atan atan2 sinh cosh tanh
    asinh acosh atanh erf erfc gamma log_gamma;
  * `**` with a non-integer exponent: a real literal (`0.5`, `1.0_wp`, `2e0`),
    or a parenthesised exponent that holds a real literal, a division, a call or
    a name that does not look integer-typed (see INT_NAME below).  NOT flagged:
    `**2`, `**n`, `**(n-1)`, `**name` (a bare identifier is assumed integer; the
    type is not resolved -- the nm check, tools/check_no_vector_math.py, is the
    backstop for what a regex cannot see).

Not flagged (false-positive classes): comments, string literals, continuation
text of a comment, names that merely contain a function name (`cosine_x`,
`explicit`, `safe_exp`, `rdb_sm_exp`), type-bound / component references
(`obj%exp(...)`), a name the file declares as a variable / dummy / procedure
(`real :: exp(n)`, `function log(`), `parameter` initialisers (constant folded
by the compiler), and the library itself (`src/framework/rdb_safe_math.F90`,
`src/shared_module_utilities/rdb_safe_math_*.inc`).

Waiver, on the SAME physical line:    ! safe-math-ok: <reason>
The reason is required; a bare `! safe-math-ok:` is itself an error.

Usage:
  tools/check_safe_math.py --diff origin/beta     # default mode: only lines changed
                                                  # vs the merge base are errors
  tools/check_safe_math.py --all                  # whole-tree inventory, exit 0
  tools/check_safe_math.py --all --strict         # ... exit 1 if any site exists
  tools/check_safe_math.py [--diff REF] FILE...   # explicit files
Exit status: 0 clean, 1 findings, 2 usage / git error.
"""

from __future__ import annotations

import argparse
import collections
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

INTRINSICS = (
    "exp log log10 sin cos tan asin acos atan atan2 sinh cosh tanh "
    "asinh acosh atanh erf erfc gamma log_gamma"
).split()
SOURCE_SUFFIXES = (".F90", ".f90", ".inc")
DEFAULT_ROOTS = ("src", "app")
EXEMPT_RE = re.compile(
    r"^src/(framework/rdb_safe_math\.F90|shared_module_utilities/rdb_safe_math_[^/]*\.inc)$"
)

WAIVER_RE = re.compile(r"!\s*safe-math-ok:(.*)$", re.I)
CALL_RE = re.compile(
    r"(?<![\w%])(" + "|".join(sorted(INTRINSICS, key=len, reverse=True)) + r")\s*\(",
    re.I,
)
# A real literal: 1.0, .5, 1., 2e0, 1d-3, anything with a kind suffix after a dot/exponent.
REAL_LIT_RE = re.compile(
    r"(?<![\w.])(?:\d+\.\d*(?![a-zA-Z_]*\()|\.\d+|\d+[eEdD][+-]?\d+)(?:_\w+)?", re.I
)
# Names that read as integers by the old implicit-typing convention (i-n).  The code
# base is `implicit none`, so this is only a heuristic for a parenthesised exponent.
INT_NAME = re.compile(r"^[i-nI-N]\w*$")
PROC_RE = re.compile(
    r"^\s*(?!end\b)(?:[\w()=*,:\s]*?\b)?(function|subroutine)\s+(\w+)", re.I
)
END_PROC_RE = re.compile(r"^\s*end\s+(function|subroutine)\b", re.I)


# --------------------------------------------------------------------------
# Lexing helpers
# --------------------------------------------------------------------------
def mask_code(line: str) -> tuple[str, str]:
    """Return (code, comment): `code` has string contents blanked (quotes kept)
    and the trailing `!` comment removed; `comment` is the removed text."""
    out: list[str] = []
    quote = ""
    for i, ch in enumerate(line):
        if quote:
            out.append(ch if ch == quote else " ")
            if ch == quote:
                quote = ""
        elif ch in "'\"":
            quote = ch
            out.append(ch)
        elif ch == "!":
            return "".join(out), line[i:]
        else:
            out.append(ch)
    return "".join(out), ""


def declared_shadows(codes: list[str]) -> set[str]:
    """Names of intrinsics that this file declares as variables / dummies /
    procedures, so that `name(` is not an intrinsic call."""
    shadows: set[str] = set()
    names = set(INTRINSICS)
    for code in codes:
        low = code.lower()
        m = re.search(r"\b(?:function|subroutine|interface)\s+(\w+)", low)
        if m and m.group(1) in names:
            shadows.add(m.group(1))
        if "::" not in low:
            continue
        rhs = low.split("::", 1)[1]
        depth, start = 0, 0
        pieces = []
        for i, ch in enumerate(rhs):
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1
            elif ch == "," and depth == 0:
                pieces.append(rhs[start:i])
                start = i + 1
        pieces.append(rhs[start:])
        for p in pieces:
            m = re.match(r"\s*&?\s*(\w+)", p)
            if m and m.group(1) in names:
                shadows.add(m.group(1))
    return shadows


def power_exponent(code: str, pos: int) -> str | None:
    """Exponent text following the `**` that ends at `pos`, or None."""
    n = len(code)
    i = pos
    while i < n and code[i] == " ":
        i += 1
    if i < n and code[i] in "+-":
        i += 1
        while i < n and code[i] == " ":
            i += 1
    if i >= n:
        return None
    if code[i] == "(":
        depth, j = 0, i
        while j < n:
            if code[j] == "(":
                depth += 1
            elif code[j] == ")":
                depth -= 1
                if depth == 0:
                    return code[i : j + 1]
            j += 1
        return code[i:]  # continues on the next line
    m = re.match(r"[\w.]+(?:\([^()]*\))?", code[i:])
    return m.group(0) if m else None


def exponent_is_real(exp_text: str) -> bool:
    t = exp_text.strip()
    if t.startswith("("):
        inner = t[1:-1] if t.endswith(")") else t[1:]
        if REAL_LIT_RE.search(inner) or "/" in inner:
            return True
        if re.search(r"\b(?:real|dble|float|sqrt|exp|log)\s*\(", inner, re.I):
            return True
        for name in re.findall(r"[A-Za-z_]\w*", inner):
            if re.fullmatch(r"\d*[eEdD]\d*|\d+_\w+", name):
                continue
            if not INT_NAME.match(name):
                return True
        return False
    if REAL_LIT_RE.fullmatch(t) or REAL_LIT_RE.match(t):
        return True
    if re.match(r"(?:real|dble|float)\s*\(", t, re.I):
        return True
    return False  # integer literal, bare identifier, array element: assumed integer


# --------------------------------------------------------------------------
# Scanner
# --------------------------------------------------------------------------
class Hit:
    __slots__ = ("path", "line", "kind", "what", "proc", "text")

    def __init__(self, path, line, kind, what, proc, text):
        self.path, self.line, self.kind = path, line, kind
        self.what, self.proc, self.text = what, proc, text


def scan_text(text: str, path: str = "<mem>") -> list[Hit]:
    """All un-waived findings in one source text (waiver errors included, kind 'waiver')."""
    raw = text.splitlines()
    split = [mask_code(r) for r in raw]
    shadows = declared_shadows([c for c, _ in split])
    hits: list[Hit] = []
    stack: list[str] = []
    for n, (code, comment) in enumerate(split, 1):
        stripped = code.lstrip()
        if stripped.startswith("#"):  # cpp directive
            continue
        # procedure context (for the inventory)
        if END_PROC_RE.match(code):
            if stack:
                stack.pop()
        else:
            pm = PROC_RE.match(code)
            if pm and "::" not in code.split(pm.group(1), 1)[0]:
                stack.append(pm.group(2).lower())
        proc = stack[-1] if stack else "<module>"
        found: list[tuple[str, str]] = []
        is_param = bool(re.search(r"\bparameter\b", code, re.I)) and "::" in code
        if not is_param:
            for m in CALL_RE.finditer(code):
                name = m.group(1).lower()
                if name in shadows:
                    continue
                found.append(("call", name))
        if not is_param:
            for m in re.finditer(r"\*\*", code):
                ex = power_exponent(code, m.end())
                if ex is not None and exponent_is_real(ex):
                    found.append(("pow", "**" + ex.strip()))
        if not found:
            continue
        w = WAIVER_RE.search(comment)
        if w:
            if not w.group(1).strip():
                hits.append(Hit(path, n, "waiver", "safe-math-ok needs a reason", proc, raw[n - 1]))
            continue
        for kind, what in found:
            hits.append(Hit(path, n, kind, what, proc, raw[n - 1]))
    return hits


def exempt(rel: str) -> bool:
    return bool(EXEMPT_RE.match(rel))


def walk(paths: list[Path]) -> list[Path]:
    files: list[Path] = []
    for p in paths:
        if p.is_file():
            files.append(p)
        elif p.is_dir():
            files += [f for f in sorted(p.rglob("*")) if f.suffix in SOURCE_SUFFIXES]
    return files


def rel_of(f: Path) -> str:
    try:
        return f.resolve().relative_to(ROOT).as_posix()
    except ValueError:
        return f.as_posix()


def scan_files(paths: list[Path]) -> list[Hit]:
    hits: list[Hit] = []
    for f in walk(paths):
        rel = rel_of(f)
        if exempt(rel) or f.suffix not in SOURCE_SUFFIXES:
            continue
        hits += scan_text(f.read_text(errors="replace"), rel)
    return hits


# --------------------------------------------------------------------------
# Diff awareness
# --------------------------------------------------------------------------
def parse_unified0(diff: str) -> dict[str, set[int]]:
    result: dict[str, set[int]] = {}
    cur: str | None = None
    for line in diff.splitlines():
        if line.startswith("+++ "):
            cur = line[6:] if line.startswith("+++ b/") else None
            if cur is not None:
                result.setdefault(cur, set())
        elif line.startswith("@@") and cur is not None:
            m = re.match(r"@@ -\d+(?:,\d+)? \+(\d+)(?:,(\d+))? @@", line)
            if m:
                start, cnt = int(m.group(1)), int(m.group(2) if m.group(2) is not None else 1)
                result[cur].update(range(start, start + cnt))
    return result


def changed_lines(ref: str) -> dict[str, set[int]]:
    """Lines added or changed in the working tree relative to merge-base(ref, HEAD)."""
    base = subprocess.run(
        ["git", "merge-base", ref, "HEAD"], capture_output=True, text=True, cwd=ROOT
    )
    if base.returncode != 0:
        raise RuntimeError(f"no merge base with '{ref}': {base.stderr.strip()}")
    out = subprocess.run(
        ["git", "diff", "--unified=0", "--no-color", base.stdout.strip(), "--"],
        capture_output=True, text=True, check=True, cwd=ROOT,
    ).stdout
    return parse_unified0(out)


# --------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------
def inventory(hits: list[Hit]) -> str:
    sites = [h for h in hits if h.kind != "waiver"]
    lines: list[str] = []
    by_dir: collections.Counter = collections.Counter()
    by_file: collections.Counter = collections.Counter()
    by_func: collections.Counter = collections.Counter()
    by_proc: collections.Counter = collections.Counter()
    for h in sites:
        parts = h.path.split("/")
        by_dir["/".join(parts[:3]) if parts[0] == "src" and len(parts) > 3 else "/".join(parts[:-1])] += 1
        by_file[h.path] += 1
        by_func["real **" if h.kind == "pow" else h.what] += 1
        by_proc[(h.path, h.proc)] += 1
    lines.append(f"safe-math inventory: {len(sites)} bare sites in {len(by_file)} files")
    lines.append("\n-- per function --")
    lines += [f"{c:6d}  {k}" for k, c in by_func.most_common()]
    lines.append("\n-- per directory --")
    lines += [f"{c:6d}  {k}" for k, c in sorted(by_dir.items(), key=lambda kv: (-kv[1], kv[0]))]
    lines.append("\n-- per file --")
    lines += [f"{c:6d}  {k}" for k, c in sorted(by_file.items(), key=lambda kv: (-kv[1], kv[0]))]
    lines.append("\n-- per procedure (file :: procedure) --")
    lines += [f"{c:6d}  {f} :: {p}" for (f, p), c in sorted(by_proc.items(), key=lambda kv: (-kv[1], kv[0]))]
    return "\n".join(lines)


def format_hit(h: Hit) -> str:
    if h.kind == "waiver":
        return f"{h.path}:{h.line}: {h.what}"
    if h.kind == "pow":
        return f"{h.path}:{h.line}: real-exponent power `{h.what}` -- use rdb_safe_math pow/cbrt  [{h.text.strip()}]"
    return f"{h.path}:{h.line}: bare `{h.what}(` -- use rdb_safe_math `safe_{h.what}`  [{h.text.strip()}]"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("paths", nargs="*", help="files / dirs (default: src app)")
    ap.add_argument("--diff", metavar="REF", help="only lines changed vs merge-base(REF, HEAD) are errors")
    ap.add_argument("--all", action="store_true", help="whole-tree inventory; exit 0 unless --strict")
    ap.add_argument("--strict", action="store_true", help="with --all: exit 1 if any site exists")
    args = ap.parse_args(argv)
    paths = [ROOT / p if not Path(p).is_absolute() else Path(p) for p in (args.paths or DEFAULT_ROOTS)]
    hits = scan_files(paths)

    if args.all:
        print(inventory(hits))
        bad = [h for h in hits if h.kind == "waiver"]
        for h in bad:
            print(format_hit(h))
        return 1 if (args.strict and hits) or bad else 0

    if not args.diff:
        ap.error("give --diff REF (hook mode) or --all (inventory)")
    try:
        changed = changed_lines(args.diff)
    except (RuntimeError, subprocess.CalledProcessError) as e:
        print(f"check_safe_math: cannot compute diff against '{args.diff}': {e}", file=sys.stderr)
        return 2
    new = [h for h in hits if h.line in changed.get(h.path, ())]
    for h in new:
        print(format_hit(h))
    if new:
        print(
            f"\ncheck_safe_math: {len(new)} new bare transcendental(s). Use rdb_safe_math, or waive on the "
            "line with `! safe-math-ok: <reason>`.  Inventory: tools/check_safe_math.py --all",
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
