#!/usr/bin/env python3
"""Lint: the vanished-layer rule must not be re-invented per kernel.

Invariant I1' — `h <= H_VANISHED  =>  hTr = h*c_live`, the donor live
layer's concentration — and its companion "the concentration of layer k"
have exactly ONE definition in this tree:
`src/shared_module_utilities/rdb_vanished_layer.inc` (`rdb_vl_is_live`,
`rdb_vl_conc`, `rdb_vl_column_conc`, `rdb_vl_holds_live_conc`,
`rdb_vl_merge_content`).  Every hand-rolled copy of that
rule is a place it can drift, and the drift is silent: the day-16 z_fixed
salt/heat break was a guard that was right on the read side and absent on
the write side of the SAME routine.

So this lint flags, on NEWLY ADDED lines only:

  [vanished-raw-divide]
      a raw `<something>hTr.../<something>h...` division — a concentration
      recovered from a thickness without going through `rdb_vl_conc` (or a
      documented, different substitution).

  [vanished-threshold]
      a comparison against `H_VANISHED` outside the sanctioned modules —
      a second spelling of "is this layer live?".

Neither is wrong per se.  Plenty of consumers legitimately substitute
something OTHER than zero on a vanished layer (the EOS wants reference T/S
so `rho = rho_0`; the sponge wants the nearest massive layer; the
diagnostics want NaN; the melt far-field sampler wants to skip the layer
entirely).  Those are different, deliberate rules, and they are exactly the
places a reader must be told which rule is in force.  Hence the waiver:

    ! vanished-ok: <reason>

on the offending line, or on any of the three lines above it.  State WHICH
substitution is intended and why `rdb_vl_conc` is not it.

Sanctioned modules (the rule's own home + the places that define the
families of vanished layers) are exempt entirely -- see SANCTIONED below.

Diff-aware: `--diff REF` reports only lines added relative to REF, so
pre-existing sites on main stay non-blocking and only new ones fail.  This
mirrors `tools/dc_assumed_shape_lint.py`.

Usage:

    python3 tools/vanished_layer_lint.py src/ app/ --diff origin/main

Stdlib only.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path
from typing import Iterable

SOURCE_SUFFIXES = {".f", ".F", ".f90", ".F90"}

WAIVER_TOKEN = "vanished-ok"
WAIVER_LOOKBACK = 4  # the line itself plus three above it

# Modules that OWN the rule or own the geometry that creates vanished
# layers.  They are expected to spell `H_VANISHED` out and to divide by a
# thickness; asking them for a waiver on every line would be noise.
SANCTIONED = {
    # the rule itself
    "src/shared_module_utilities/rdb_vanished_layer.inc",
    # the constant of record
    "src/core/rdb_constants.F90",
    # the coordinate: it MAKES the fillers, and owns the h_min contract
    "src/ALE/rdb_vcoord.F90",
    "src/core/ocean/vcoord/rdb_ocean_vcoord.F90",
    # the remap: the rule's first consumer, and the k_top / closed-face
    # producers that translate "vanished" into an index or a mask
    "src/ALE/rdb_ocean_remap.F90",
    # the state type that owns the enforcement point + the I1' scan
    "src/core/rdb_multilayer_state.F90",
}

# A division whose numerator mentions a tracer-CONTENT name and whose
# denominator mentions a THICKNESS name.  Deliberately narrow: it is meant
# to catch `hTr(i,j,k)/h_layer(i,j,k)` and its spellings, not every `/`.
_CONTENT = r"(?:hTr|htr|hT\b|hS\b|temp_h|salt_h|thtr|shtr|htr_col|hTr_col)"
_THICK = r"(?:h_layer|h_old|h_new|h_col|h_sd|hh|dz|h\b|hk\b|he\b|hc\b)"
_DIVIDE_RE = re.compile(
    rf"{_CONTENT}[A-Za-z0-9_]*\s*(?:\([^()]*\))?\s*/\s*"
    rf"(?:max\s*\(\s*)?{_THICK}[A-Za-z0-9_]*"
)

_THRESHOLD_RE = re.compile(r"[<>=]=?\s*H_VANISHED|H_VANISHED\s*[<>=]=?")


def walk_files(paths: list[Path]) -> Iterable[Path]:
    for p in paths:
        if p.is_file() and p.suffix in SOURCE_SUFFIXES:
            yield p
        elif p.is_dir():
            for f in sorted(p.rglob("*")):
                if f.suffix in SOURCE_SUFFIXES:
                    yield f


def changed_lines(ref: str) -> dict[Path, set[int]]:
    """Map path -> set of line numbers ADDED relative to `ref`."""
    out = subprocess.run(
        ["git", "diff", "--unified=0", "--no-color", ref, "--"],
        capture_output=True, text=True, check=True,
    ).stdout
    result: dict[Path, set[int]] = {}
    current: Path | None = None
    for line in out.splitlines():
        if line.startswith("+++ b/"):
            current = Path(line[6:])
            result.setdefault(current, set())
        elif line.startswith("@@") and current is not None:
            m = re.match(r"@@ -\d+(?:,\d+)? \+(\d+)(?:,(\d+))? @@", line)
            if m:
                start = int(m.group(1))
                n = int(m.group(2) or "1")
                for ln in range(start, start + n):
                    result[current].add(ln)
    return result


def has_waiver(lines: list[str], lineno: int) -> bool:
    """Waiver on the line itself or any of the three comment lines above."""
    for look in range(WAIVER_LOOKBACK):
        idx = lineno - 1 - look
        if idx < 0:
            break
        text = lines[idx]
        if WAIVER_TOKEN in text:
            return True
        if look > 0:
            stripped = text.strip()
            if stripped and not stripped.startswith("!"):
                break
    return False


def strip_comment(line: str) -> str:
    """Return only the executable text: trailing `!` comment dropped and the
    CONTENTS of every string literal blanked out.

    Blanking the literals matters: a fail-loud message that quotes the rule
    ("h <= H_VANISHED => hTr = h*c_live") is documentation, not a second spelling of
    the test, and must not trip the lint."""
    out = []
    quote = ""
    for ch in line:
        if quote:
            if ch == quote:
                quote = ""
                out.append(ch)
            else:
                out.append(" ")
        elif ch in "'\"":
            quote = ch
            out.append(ch)
        elif ch == "!":
            break
        else:
            out.append(ch)
    return "".join(out)


def lint_file(path: Path, diff_filter: set[int] | None) -> list[tuple[int, str, str]]:
    """Return [(lineno, rule_id, source_text)] for this file."""
    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return []
    hits: list[tuple[int, str, str]] = []
    for i, raw in enumerate(lines):
        lineno = i + 1
        if diff_filter is not None and lineno not in diff_filter:
            continue
        code = strip_comment(raw)
        if not code.strip():
            continue
        rule = None
        if _DIVIDE_RE.search(code):
            rule = "vanished-raw-divide"
        elif _THRESHOLD_RE.search(code):
            rule = "vanished-threshold"
        if rule is None:
            continue
        if has_waiver(lines, lineno):
            continue
        hits.append((lineno, rule, raw))
    return hits


def lint(paths: list[Path], diff_lines: dict[Path, set[int]] | None):
    results = []
    for f in walk_files(paths):
        key = f.as_posix()
        if key in SANCTIONED:
            continue
        filt = None
        if diff_lines is not None:
            if f not in diff_lines:
                continue
            filt = diff_lines[f]
            if not filt:
                continue
        for lineno, rule, text in lint_file(f, filt):
            results.append((f, lineno, rule, text))
    return results


ADVICE = {
    "vanished-raw-divide": (
        "a tracer content divided by a thickness is a CONCENTRATION, and this "
        "tree has one definition of it",
        "use `rdb_vl_conc(hTr, h)` from "
        "`src/shared_module_utilities/rdb_vanished_layer.inc` "
        "(`#include` it in this module's `contains`)",
    ),
    "vanished-threshold": (
        "a comparison against H_VANISHED is a second spelling of "
        "\"is this layer live?\"",
        "use `rdb_vl_is_live(h)` from "
        "`src/shared_module_utilities/rdb_vanished_layer.inc`",
    ),
}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("paths", nargs="*", default=["src", "app"], type=Path,
                    help="files or directories to scan (default: src app)")
    ap.add_argument("--diff", metavar="REF",
                    help="only report findings on lines added vs REF")
    args = ap.parse_args()

    diff_lines: dict[Path, set[int]] | None = None
    if args.diff:
        try:
            diff_lines = changed_lines(args.diff)
        except subprocess.CalledProcessError as exc:
            print(f"vanished_layer_lint: cannot resolve diff base '{args.diff}' "
                  f"(exit {exc.returncode}) — skipping cleanly.")
            return 0

    hits = lint([Path(p) for p in args.paths], diff_lines)

    if not hits:
        print("vanished_layer_lint: clean")
        return 0

    for path, lineno, rule, text in hits:
        why, fix = ADVICE[rule]
        print(f"ERROR  {path}:{lineno}  [{rule}]  {why}")
        print(f"       {text.strip()}")
        print(f"  FIX: {fix}.")
        print(f"  WAIVER: if this site deliberately substitutes something ELSE "
              f"on a vanished layer (reference T/S for the EOS, the nearest "
              f"massive layer, NaN, a skip), say so with "
              f"`! {WAIVER_TOKEN}: <reason>` on or above this line.")
        print()

    print(f"vanished_layer_lint: {len(hits)} error(s)")
    return 1


if __name__ == "__main__":
    sys.exit(main())
