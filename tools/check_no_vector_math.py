#!/usr/bin/env python3
"""Fail when a project object file references a compiler vector-math entry point.

A loop that calls the intrinsic `exp`/`sin`/`pow`... may be compiled to a call of a
VECTOR math routine in its vector body and the scalar libm routine in its
remainder, so a cell rounds differently depending on its position in the tile
(`plans/multirank/SAFE_MATH_PLAN.md`, section 0.4 / 3b).  `tools/check_safe_math.py`
finds the call sites it can read; this tool is the backstop on the compiled
objects: it runs `nm -u` over the model's own objects (`src/` and `app/` sources: everything
under the build directory except `_deps/`, CMake's probe directories, tests and benchmarks) and looks for
the undefined symbols below.

Symbol patterns (regex, matched against each undefined symbol):
  glibc libmvec (gfortran, x86_64):  ^_ZGV[bcde][MN]\\d+[vl]+_\\w+   e.g. _ZGVdN4v_exp, _ZGVbN2vv_pow
  NVHPC libpgmath vector entries:    ^__f[ds]_\\w+   (__fd_exp_4_avx2, __fd_sin_2_avx2, __fd_exp_4_mn;
                                     fd = double, fs = single)
                                     ^__p[ds]_\\w+   (packed/array forms, double / single)
                                     ^__gv[ds]_\\w+  (generic-vector forms, double / single)
                                     ^__mth_i_\\w+   (scalar runtime math, e.g. __mth_i_dpowd)
The glibc pattern was taken from a real gfortran 15.1 build of this tree
(`-O3 -march=native -ffp-contract=off`, glibc's math-vector-fortran.h pre-include);
the NVHPC patterns are from the plan's inventory notes (`__fd_*`, `__fs_*`,
`__mth_i_*`) plus the sibling `__pd_/__ps_/__gvd_/__gvs_` families, which no object
here has been observed to use: they are listed so they cannot slip in unnoticed.
The allow-list below was seeded from a GNU build only.

Allow-list (`tools/vector_math_allowlist.txt`): one object per line, the path
relative to its `<target>.dir/`, optionally followed by `  # symbols`.  The check
FAILS when
  * an object not in the list references a pattern (a new offender), or
  * a listed object that exists in the build no longer references any (STALE: delete
    the line -- the same XPASS rule as tests/regression/compat_expect.py).
A listed object missing from the build dir is only a warning (partial build).

Usage:
  tools/check_no_vector_math.py --build-dir build_gfortran
  tools/check_no_vector_math.py --build-dir build_gfortran --write-allowlist   # re-seed
Exit status: 0 ok, 1 new offender or stale entry, 2 usage / no objects found.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_ALLOWLIST = Path(__file__).resolve().parent / "vector_math_allowlist.txt"

PATTERNS = [
    re.compile(r"^_ZGV[a-zA-Z]+\d*[a-zA-Z]+_\w+$"),
    re.compile(r"^__f[ds]_\w+$"),
    re.compile(r"^__p[ds]_\w+$"),
    re.compile(r"^__gv[ds]_\w+$"),
    re.compile(r"^__mth_i_\w+$"),
]
SKIP_DIRS = ("_deps", "CMakeTmp")
MODEL_PREFIXES = ("src/", "app/")  # tests/ and benchmarks/ objects are not the model


def object_key(obj: Path, build: Path) -> str:
    """Path of the object relative to its `<target>.dir/` (or to the build dir)."""
    rel = obj.relative_to(build).as_posix()
    m = re.search(r"(?:^|/)[^/]+\.dir/(.*)$", rel)
    return m.group(1) if m else rel


def project_objects(build: Path) -> list[Path]:
    out = []
    for o in sorted(build.rglob("*.o")):
        parts = o.relative_to(build).parts
        if any(p in SKIP_DIRS or p.startswith("CompilerId") for p in parts):
            continue
        if object_key(o, build).startswith(MODEL_PREFIXES):
            out.append(o)
    return out


def vector_symbols(obj: Path) -> list[str]:
    r = subprocess.run(["nm", "-u", str(obj)], capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"nm -u {obj}: {r.stderr.strip()}")
    syms = set()
    for line in r.stdout.splitlines():
        parts = line.split()
        if parts and any(p.match(parts[-1]) for p in PATTERNS):
            syms.add(parts[-1])
    return sorted(syms)


def parse_allowlist(text: str) -> set[str]:
    keys = set()
    for line in text.splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            keys.add(line)
    return keys


def compare(found: dict[str, list[str]], allowed: set[str], present: set[str]):
    """-> (new, stale, missing): offenders not listed, listed-but-clean, listed-but-not-built."""
    new = sorted(k for k in found if k not in allowed)
    stale = sorted(k for k in allowed if k in present and k not in found)
    missing = sorted(k for k in allowed if k not in present)
    return new, stale, missing


def format_allowlist(found: dict[str, list[str]]) -> str:
    head = (
        "# Objects that reference compiler vector-math entry points (see tools/check_no_vector_math.py).\n"
        "# Each line is a burn-down item: remove it when the object stops referencing them (the check fails\n"
        "# on a stale line).  Seeded from a GNU build; regenerate with --write-allowlist.\n"
    )
    return head + "".join(f"{k}  # {' '.join(v)}\n" for k, v in sorted(found.items()))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--build-dir", required=True, type=Path)
    ap.add_argument("--allowlist", type=Path, default=DEFAULT_ALLOWLIST)
    ap.add_argument("--write-allowlist", action="store_true", help="overwrite the allow-list from this build")
    args = ap.parse_args(argv)
    build = args.build_dir.resolve()
    objs = project_objects(build) if build.is_dir() else []
    if not objs:
        print(f"check_no_vector_math: no object files under {build}", file=sys.stderr)
        return 2
    found: dict[str, list[str]] = {}
    present: set[str] = set()
    for o in objs:
        key = object_key(o, build)
        present.add(key)
        syms = vector_symbols(o)
        if syms:
            found[key] = syms
    if args.write_allowlist:
        args.allowlist.write_text(format_allowlist(found))
        print(f"check_no_vector_math: wrote {len(found)} objects to {args.allowlist}")
        return 0
    allowed = parse_allowlist(args.allowlist.read_text()) if args.allowlist.exists() else set()
    new, stale, missing = compare(found, allowed, present)
    for k in missing:
        print(f"warning: allow-listed object not in this build: {k}")
    for k in new:
        print(f"NEW  {k}: {' '.join(found[k])}")
    for k in stale:
        print(f"STALE {k}: no longer references vector math -- delete it from {args.allowlist.name}")
    if new or stale:
        print(
            "\ncheck_no_vector_math: FAIL.  Route the call through rdb_safe_math (the adoption PRs), or, for a\n"
            "deliberate new offender, add it to tools/vector_math_allowlist.txt.",
            file=sys.stderr,
        )
        return 1
    print(f"check_no_vector_math: OK ({len(objs)} objects, {len(found)} allow-listed offenders)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
