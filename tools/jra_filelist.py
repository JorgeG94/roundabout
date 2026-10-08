#!/usr/bin/env python3
"""Write the ordered per-variable file list `register_2d_filelist`
(`rdb_ocean_data_input`, OM3 wave-1 PR-C2a) needs.

Fortran has no glob, and a JRA55-do-style forcing dataset is "one file
per variable per year" (`tas_..._gr_19580101...-19581231...nc`, one such
file per variable per calendar year, 1958..2023 + a partial 2024-01
file). `register_2d_filelist` wants one path per line, in TIME ORDER,
for a single variable -- this script builds that list by walking a
directory and matching a variable-name prefix, stdlib only (no
`glob`-equivalent dependency beyond `os`/`re`, which this box already
has; never `pip install`).

Sort key: the date-range token the JRA55-do naming convention embeds
right before the `.nc` extension (`..._YYYYMMDDhhmm-YYYYMMDDhhmm.nc`) --
sorting on the token's START date puts the files in time order even if
the directory listing itself is not alphabetical for some other reason
(it usually is, for this naming convention, but sort explicitly rather
than assume).

Usage:
    jra_filelist.py DATA_DIR tas [--pattern 'tas_*_gr_*.nc'] [-o out.txt]

With no `--pattern`, `tools/jra_filelist.py DATA_DIR VAR` matches any
file under `DATA_DIR` whose name starts with `VAR_` and ends `.nc`.
Prints one path per line to stdout, or to `-o FILE` when given -- the
file `register_2d_filelist`'s `files(:)` argument reads, one path per
line (see `rdb_ocean_data_input.F90`).
"""
from __future__ import annotations

import argparse
import os
import re
import sys

# JRA55-do / input4MIPs date-range token: two YYYYMMDDhhmm stamps
# separated by a dash, immediately before ".nc".
_DATE_RANGE_RE = re.compile(r"_(\d{8,12})-(\d{8,12})\.nc$")


def _sort_key(path: str) -> str:
    """Date-range start stamp if the filename carries one, else the
    bare filename (so a non-conforming name still sorts deterministically
    rather than raising)."""
    m = _DATE_RANGE_RE.search(os.path.basename(path))
    if m:
        return m.group(1)
    return os.path.basename(path)


def find_files(data_dir: str, var: str, pattern: str | None) -> list[str]:
    """Walk `data_dir` (recursively) collecting files that match either
    an explicit shell-glob-style `pattern` (fnmatch, stdlib) or, with no
    pattern given, any file named `<var>_*.nc`."""
    import fnmatch

    glob_pat = pattern if pattern else f"{var}_*.nc"
    matches = []
    for root, _dirs, files in os.walk(data_dir):
        for name in files:
            if fnmatch.fnmatch(name, glob_pat):
                matches.append(os.path.join(root, name))
    matches.sort(key=_sort_key)
    return matches


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("data_dir", help="directory to walk (recursively)")
    ap.add_argument("var", help="variable name (JRA55-do CMOR name, e.g. 'tas')")
    ap.add_argument("--pattern", default=None,
                    help="explicit fnmatch pattern (default: '<var>_*.nc')")
    ap.add_argument("-o", "--output", default=None,
                    help="write the list here (default: stdout)")
    args = ap.parse_args(argv)

    files = find_files(args.data_dir, args.var, args.pattern)
    if not files:
        print(f"jra_filelist: no files matched under {args.data_dir!r} "
              f"for var={args.var!r} pattern={args.pattern!r}", file=sys.stderr)
        return 1

    text = "\n".join(files) + "\n"
    if args.output:
        with open(args.output, "w") as f:
            f.write(text)
    else:
        sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
