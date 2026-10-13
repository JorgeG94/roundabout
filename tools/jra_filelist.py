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
    jra_filelist.py --streams datm.streams.xml [--stream-name NAME] [-o out.txt]

With no `--pattern`, `tools/jra_filelist.py DATA_DIR VAR` matches any
file under `DATA_DIR` whose name starts with `VAR_` and ends `.nc`.
Prints one path per line to stdout, or to `-o FILE` when given -- the
file `register_2d_filelist`'s `files(:)` argument reads, one path per
line (see `rdb_ocean_data_input.F90`).

`--streams XML` is the second mode (OM3 C2c-3b): instead of walking a
directory, it reads a DATM data-model stream descriptor -- the
`datm.streams.xml` every OM3/CESM-family coupled run already carries,
one `<stream_info name="...">` block per named stream, each owning a
`<datafiles><file>...</file>...</datafiles>` list already in the time
order the stream reader expects. `--stream-name` selects one block when
the file defines more than one (the usual case: a JRA55-do deployment
typically splits the nine fields across a handful of named streams);
omitted, every stream's `<file>` entries are concatenated in document
order. Same stdout/`-o` output contract as the directory-walk mode.
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


def find_files_streams(xml_path: str, stream_name: str | None) -> list[str]:
    """Extract `<file>` entries from a DATM `datm.streams.xml`-style
    descriptor, stdlib `xml.etree.ElementTree` only (no new dependency).
    Walks every `<stream_info>` element in document order; `name` is
    read from either a `name=` attribute or a `<name>`/`<streamname>`
    child (the two spellings this format is seen with), and a block is
    skipped when `stream_name` is given and does not match. Every
    `<file>` text under that block's `<datafiles>` is appended in
    document order (the order the stream reader consumes them in) --
    concatenated across every surviving stream when more than one
    matches (or when `stream_name` is omitted)."""
    import xml.etree.ElementTree as ET

    tree = ET.parse(xml_path)
    root = tree.getroot()
    files: list[str] = []
    for stream in root.iter("stream_info"):
        name = stream.get("name")
        if name is None:
            name = stream.findtext("name") or stream.findtext("streamname") or ""
        if stream_name and name != stream_name:
            continue
        datafiles = stream.find("datafiles")
        if datafiles is None:
            continue
        for f in datafiles.findall("file"):
            if f.text and f.text.strip():
                files.append(f.text.strip())
    return files


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("data_dir", nargs="?", default=None,
                    help="directory to walk (recursively); unused with --streams")
    ap.add_argument("var", nargs="?", default=None,
                    help="variable name (JRA55-do CMOR name, e.g. 'tas'); unused with --streams")
    ap.add_argument("--pattern", default=None,
                    help="explicit fnmatch pattern (default: '<var>_*.nc')")
    ap.add_argument("--streams", default=None, metavar="XML",
                    help="parse a datm.streams.xml instead of walking a directory")
    ap.add_argument("--stream-name", default=None,
                    help="--streams only: the <stream_info name=...> to extract "
                        "(default: concatenate every stream in the file)")
    ap.add_argument("-o", "--output", default=None,
                    help="write the list here (default: stdout)")
    args = ap.parse_args(argv)

    if args.streams:
        files = find_files_streams(args.streams, args.stream_name)
        if not files:
            print(f"jra_filelist: no <file> entries found in {args.streams!r} "
                  f"(stream_name={args.stream_name!r})", file=sys.stderr)
            return 1
    else:
        if not args.data_dir or not args.var:
            ap.error("data_dir and var are required unless --streams is given")
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
