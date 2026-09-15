#!/usr/bin/env python3
"""Add explicit `map()` clauses for OPTIONAL array dummies read inside a
`!$omp target` region.

Why
---
LLVM Flang emits an *unconditional implicit* map for an `optional`,
EXPLICIT-SHAPE array dummy referenced inside a `target` region.  The map size
comes from the declared bounds — which are computable whether or not the
argument is present — but the base address of an absent argument is NULL, so
libomptarget calls `hsa_amd_memory_lock(0x0, <size>)` and the run dies with

    PluginInterface error: Failure to copy data from host to device.
    Pointers: host = 0x0, device = 0x..., size = 1438240 ...
    HSA_STATUS_ERROR_INVALID_ARGUMENT

Flang gets this right for descriptor-passed dummies (assumed-shape /
allocatable / pointer), where it null-checks the box before building the
bounds; an explicit-shape dummy has no box to check.  Naming the argument in
an explicit `map()` clause takes the guarded path, and the OpenMP rule that a
map of an absent optional dummy is ignored then applies.

This runs as a POST-PASS over the already-translated tree, so it covers target
regions emitted by BOTH translators — `dc_to_omp.py` (from `do concurrent`)
and `acc_to_omp.py` (from `!$acc parallel loop`, which carry `reduction`).
`main` is untouched: OpenACC resolves these references through the host
present-table and never hits the bug.

Map type follows the dummy's intent: `to` for `intent(in)`, `tofrom`
otherwise.  For an argument that is already device-resident (the normal case
here — state arrays are attached by `enter_data`) the clause is a present-table
lookup and moves no data either way.

Idempotent: a directive that already maps the name is left alone.

Usage:
    python tools/omp_optional_map.py [--write] [-j N] src/ app/ benchmarks/ tests/
"""

from __future__ import annotations

import argparse
import re
import sys
from functools import partial
from pathlib import Path

from _parallel import pmap, resolve_jobs
from _progress import track

# An *opening* `!$omp target` worksharing directive. Excludes `!$omp end
# target`, the declarative `!$omp declare target`, and the standalone data
# constructs (`target enter/exit data`, `target update`) — none of which
# enclose a region whose body could reference a dummy.
OMP_TARGET_OPEN_RE = re.compile(
    r"^\s*!\$omp\s+target\b(?!\s+(?:enter|exit)\s+data\b)(?!\s+update\b)",
    re.IGNORECASE,
)
OMP_SENTINEL_RE = re.compile(r"^\s*!\$omp\b", re.IGNORECASE)
DO_RE = re.compile(r"^\s*(?:[a-zA-Z_][a-zA-Z0-9_]*\s*:\s*)?do\b", re.IGNORECASE)
END_DO_RE = re.compile(r"^\s*end\s*do\b", re.IGNORECASE)
PROC_OPEN_RE = re.compile(
    r"^\s*[^!]*\b(?:subroutine|function)\b\s+[a-zA-Z]\w*", re.IGNORECASE
)
PROC_END_RE = re.compile(r"^\s*end\s*(?:subroutine|function)\b", re.IGNORECASE)
INTENT_IN_RE = re.compile(r"\bintent\s*\(\s*in\s*\)", re.IGNORECASE)
OPTIONAL_RE = re.compile(r"\boptional\b", re.IGNORECASE)

# Free-form source limit; wrap the appended clause onto an OpenMP continuation
# rather than push a directive past it.
MAX_LINE = 132


def strip_strings(line: str) -> str:
    out, quote = [], None
    for ch in line:
        if quote:
            out.append(" " if ch != quote else ch)
            if ch == quote:
                quote = None
        else:
            if ch in ("'", '"'):
                quote = ch
            out.append(ch)
    return "".join(out)


def code_part(line: str) -> str:
    """Line with any trailing comment removed. OpenMP sentinels are kept: an
    `!$omp` directive IS code for our purposes."""
    s = strip_strings(line)
    if OMP_SENTINEL_RE.match(s):
        return line
    bang = s.find("!")
    return line[:bang] if bang >= 0 else line


def split_top_level(s: str, sep: str = ",") -> list[str]:
    parts, depth, cur = [], 0, ""
    for ch in s:
        if ch in "([":
            depth += 1
        elif ch in ")]":
            depth -= 1
        if ch == sep and depth == 0:
            parts.append(cur)
            cur = ""
        else:
            cur += ch
    parts.append(cur)
    return parts


def logical_lines(lines: list[str], lo: int, hi: int):
    """Yield (start_index, joined_code) for each continued statement in
    [lo, hi), so a declaration split across `&` continuations is seen whole."""
    i = lo
    while i < hi:
        c = code_part(lines[i]).rstrip()
        if OMP_SENTINEL_RE.match(c):
            i += 1
            continue
        start = i
        joined = c
        while joined.endswith("&") and i + 1 < hi:
            i += 1
            nxt = code_part(lines[i]).strip()
            joined = joined[:-1].rstrip() + " " + nxt.lstrip("&").strip()
        yield start, joined
        i += 1


def parse_optional_array_dummies(lines: list[str], lo: int, hi: int) -> dict[str, str]:
    """Map lowercase name -> map-type for every OPTIONAL ARRAY dummy declared
    in the procedure spanning [lo, hi)."""
    found: dict[str, str] = {}
    for _, stmt in logical_lines(lines, lo, hi):
        if "::" not in stmt:
            continue
        head, _, tail = stmt.partition("::")
        if not OPTIONAL_RE.search(head):
            continue
        maptype = "to" if INTENT_IN_RE.search(head) else "tofrom"
        for item in split_top_level(tail):
            item = item.strip()
            if "(" not in item:
                continue  # scalar optional: passed by value/ref, no map emitted
            name, _, shape = item.partition("(")
            name = name.strip()
            shape = shape.strip()
            # Only EXPLICIT-SHAPE dummies are affected. Assumed-shape `(:,:)`
            # and assumed-size `(*)` are passed as descriptors, and flang
            # already null-checks the box before building the map bounds.
            if shape.startswith(":") or shape.startswith("*"):
                continue
            if name:
                found[name.lower()] = maptype
    return found


def parse_procedures(lines: list[str]) -> list[dict]:
    procs: list[dict] = []
    stack: list[int] = []
    for idx, raw in enumerate(lines):
        c = code_part(raw)
        if OMP_SENTINEL_RE.match(c):
            continue
        if PROC_END_RE.match(c):
            if stack:
                procs.append({"open": stack.pop(), "end": idx})
        elif PROC_OPEN_RE.match(c):
            stack.append(idx)
    return procs


def directive_last_line(lines: list[str], first: int) -> int:
    """Index of the final physical line of the directive starting at `first`
    (follows `&` continuations)."""
    i = first
    while code_part(lines[i]).rstrip().endswith("&") and i + 1 < len(lines):
        i += 1
    return i


def region_end(lines: list[str], dir_last: int) -> int:
    """Index of the last line of the region body: the `end do` matching the
    first `do` after the directive. Covers a collapsed nest (matching the
    OUTERMOST do spans the whole nest) and works whether or not the translator
    emitted a closing `!$omp end ...`."""
    i = dir_last + 1
    while i < len(lines):
        c = code_part(lines[i])
        if OMP_SENTINEL_RE.match(c) or not c.strip():
            i += 1
            continue
        if DO_RE.match(c) and not END_DO_RE.match(c):
            break
        return dir_last  # no loop found; nothing to scan
    depth = 0
    while i < len(lines):
        c = code_part(lines[i])
        if OMP_SENTINEL_RE.match(c):
            i += 1
            continue
        if END_DO_RE.match(c):
            depth -= 1
            if depth == 0:
                return i
        elif DO_RE.match(c):
            depth += 1
        i += 1
    return len(lines) - 1


def append_clause(lines: list[str], dir_first: int, dir_last: int, clause: str) -> None:
    """Append `clause` to the directive, continuing onto a new sentinel line if
    it would push past the free-form line limit."""
    line = lines[dir_last]
    if len(line.rstrip()) + 1 + len(clause) <= MAX_LINE:
        lines[dir_last] = line.rstrip() + " " + clause
        return
    indent = re.match(r"\s*", lines[dir_first]).group(0)
    lines[dir_last] = line.rstrip() + " &"
    lines.insert(dir_last + 1, f"{indent}!$omp {clause}")


def transform_lines(lines: list[str]) -> tuple[list[str], int]:
    out = list(lines)
    added = 0
    # Rebuild the procedure table each pass: append_clause can insert a
    # continuation line and shift every index after it.
    idx = 0
    while True:
        procs = parse_procedures(out)
        target = None
        for i in range(idx, len(out)):
            if OMP_TARGET_OPEN_RE.match(code_part(out[i])):
                target = i
                break
        if target is None:
            break
        enclosing = None
        for p in procs:
            if p["open"] < target < p["end"]:
                if enclosing is None or (p["end"] - p["open"]) < (
                    enclosing["end"] - enclosing["open"]
                ):
                    enclosing = p
        idx = target + 1
        if enclosing is None:
            continue
        opts = parse_optional_array_dummies(out, enclosing["open"], enclosing["end"])
        if not opts:
            continue
        dir_last = directive_last_line(out, target)
        directive = " ".join(code_part(l) for l in out[target:dir_last + 1])
        body_end = region_end(out, dir_last)
        body = "\n".join(code_part(l) for l in out[target:body_end + 1])
        wanted: dict[str, list[str]] = {}
        for name, maptype in opts.items():
            if not re.search(rf"\b{re.escape(name)}\b", body, re.IGNORECASE):
                continue
            if re.search(rf"map\s*\([^)]*\b{re.escape(name)}\b", directive, re.IGNORECASE):
                continue  # already mapped — idempotent re-run
            wanted.setdefault(maptype, []).append(name)
        for maptype in ("to", "tofrom"):
            names = sorted(wanted.get(maptype, []))
            if not names:
                continue
            append_clause(out, target, directive_last_line(out, target),
                          f"map({maptype}: {', '.join(names)})")
            added += 1
    return out, added


def _convert(path: Path) -> tuple[Path, list[str], int, bool]:
    text = path.read_text()
    lines, n = transform_lines(text.splitlines())
    return path, lines, n, text.endswith("\n")


def _write(unit, write: bool = False):
    path, lines, n, trailing_nl = unit
    if n and write:
        path.write_text("\n".join(lines) + ("\n" if trailing_nl else ""))
    return path, n


def walk(roots: list[Path]) -> list[Path]:
    files: list[Path] = []
    for r in roots:
        if r.is_file():
            files.append(r)
        else:
            files.extend(sorted(r.rglob("*.F90")))
    return files


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("roots", nargs="+", type=Path)
    ap.add_argument("--write", action="store_true", help="rewrite files in place")
    ap.add_argument("-j", "--jobs", type=int, default=1)
    args = ap.parse_args()

    files = walk(args.roots)
    jobs = resolve_jobs(args.jobs)
    units = list(track(pmap(_convert, files, jobs), total=len(files),
                       label="omp_optional_map"))
    total = 0
    touched = 0
    for path, n in pmap(partial(_write, write=args.write), units, jobs):
        if n:
            touched += 1
            total += n
            print(f"  {path}: +{n} map clause(s)")
    verb = "added" if args.write else "would add"
    print(f"{verb} {total} map clause(s) across {touched} file(s) "
          f"({len(files)} scanned)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
