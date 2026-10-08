#!/usr/bin/env python3
"""[C1] Convert a MOM6 `ocean_vgrid.nc` into a roundabout `&vcoord_nml
z_fixed_dz` list, so the ocean dyn-core can run on ACCESS-OM3's exact
vertical grid (75 levels, z* coordinate).

Usage
-----
    python3 vgrid_to_dz.py ocean_vgrid.nc
    python3 vgrid_to_dz.py --text-interfaces interfaces.txt
    python3 vgrid_to_dz.py --self-test

Input formats
-------------
1. A MOM6 "supergrid" vertical-grid file (what `ocean_vgrid.nc` on Gadi is):
   a NetCDF file with a 1-D variable `zeta` holding `2*NK + 1` interface
   depths (m, positive down, `zeta[0] = 0` at the surface) -- the supergrid
   doubles the resolution of the NK model-grid interfaces, so the model
   interfaces are the EVEN-index entries, `zeta[::2]`.  Requires the
   `netCDF4` package (present in this project's python3.10; never
   `pip install` -- see CLAUDE.md).  If the file instead carries a plain
   `dz` variable (already model-grid layer thicknesses, NK values), it is
   read directly.

2. Stdlib-only fallback, no NetCDF reader needed: `--text-interfaces FILE`
   reads a plain text file of `NK + 1` interface depths (m, positive down,
   starting at 0), whitespace/comma/newline separated -- i.e. the content
   of `zeta[::2]` typed or dumped by any other tool.

Both paths are "top-down": the first thickness is the surface layer, the
last is the bed. 75 is ACCESS-OM3's layer count; this tool does not
hardcode it -- it emits exactly as many `z_fixed_dz` entries as the input
has layers.

Output ordering -- READ BEFORE PASTING
---------------------------------------
roundabout's `&vcoord_nml z_fixed_dz` is documented and enforced
SURFACE-FIRST (`src/core/rdb_config.F90`: "z_fixed_profile = 'list' needs
exactly nz_layers ... leading positive z_fixed_dz entries (surface first,
no gaps)"). The internal bottom-up table the kernels actually index
(k=1 = bed, k=nz = surface, per CLAUDE.md's Vertical Layer Convention) is
built by FLIPPING this surface-first namelist input --
`ocean_vcoord_set_z_fixed_profile` in
`src/core/ocean/vcoord/rdb_ocean_vcoord.F90` does the flip.

MOM6's `zeta` (and `dz`) are ALREADY top-down / surface-first, so this
tool does NOT reverse anything it reads: the thicknesses it prints are
used verbatim as `z_fixed_dz`. Pair the printed block with
`vcoord_type = "zstar"` (ACCESS-OM3's coordinate) and
`z_fixed_profile = "list"`.

Self-test
---------
`--self-test` builds a synthetic vgrid (both via `netCDF4`, when
importable, and via the stdlib-only text-interfaces path), converts it,
and checks the recovered thicknesses against the known answer. Exits 0 on
pass, 1 on failure; no model is run.
"""

from __future__ import annotations

import argparse
import re
import sys


def interfaces_to_dz(interfaces):
    """NK+1 monotone interface depths (top-down, positive down) -> NK
    surface-first layer thicknesses."""
    dz = [interfaces[k + 1] - interfaces[k] for k in range(len(interfaces) - 1)]
    for k, d in enumerate(dz):
        if d <= 0.0:
            raise ValueError(
                "interface depths must be strictly increasing (non-positive "
                "thickness at layer {} = {})".format(k + 1, d)
            )
    return dz


def read_vgrid_netcdf(path):
    """Read a MOM6 `ocean_vgrid.nc`: prefer a plain `dz` variable (already
    model-grid, surface-first thicknesses) and fall back to deriving
    thicknesses from the supergrid `zeta` interface array."""
    try:
        import netCDF4
    except ImportError as exc:
        raise RuntimeError(
            "netCDF4 is not importable -- use --text-interfaces instead, "
            "or run this on a Python with netCDF4 (never `pip install` "
            "in this repo)"
        ) from exc

    with netCDF4.Dataset(path, "r") as ds:
        if "dz" in ds.variables:
            dz = [float(v) for v in ds.variables["dz"][:].ravel()]
            return dz
        if "zeta" not in ds.variables:
            raise ValueError(
                "'{}' has neither a 'dz' nor a 'zeta' variable".format(path)
            )
        zeta = [float(v) for v in ds.variables["zeta"][:].ravel()]

    if len(zeta) % 2 == 0:
        raise ValueError(
            "'zeta' has {} values; a supergrid interface array must have "
            "2*NK + 1 (odd) entries".format(len(zeta))
        )
    model_interfaces = zeta[::2]
    return interfaces_to_dz(model_interfaces)


def read_text_interfaces(path):
    """Stdlib-only fallback: a text file of NK+1 interface depths,
    whitespace/comma/newline separated."""
    with open(path, "r") as fh:
        text = fh.read()
    tokens = [t for t in re.split(r"[\s,]+", text.strip()) if t]
    interfaces = [float(t) for t in tokens]
    if len(interfaces) < 2:
        raise ValueError(
            "'{}' must list at least 2 interface depths (1 layer)".format(path)
        )
    return interfaces_to_dz(interfaces)


def format_nml_block(dz, line_width=88):
    """Render a ready-to-paste `&vcoord_nml z_fixed_profile = "list",
    z_fixed_dz = ...` block, surface first, matching roundabout's
    namelist formatting (see docs/generated_nml_knobs.md)."""
    values = ["{:.6f}".format(d) for d in dz]
    lines = ['&vcoord_nml', '   vcoord_type     = "zstar",', '   z_fixed_profile = "list",']
    prefix = "   z_fixed_dz      = "
    cur = prefix
    for i, v in enumerate(values):
        tok = v + ("," if i < len(values) - 1 else "")
        if len(cur) + len(tok) + 1 > line_width and cur != prefix:
            lines.append(cur)
            cur = " " * len(prefix)
        cur += tok + " "
    lines.append(cur.rstrip())
    lines.append("/")
    header = (
        "! {} layer thicknesses (m), surface first, sum = {:.3f} m\n"
        "! pair with &ocean_topo_nml max_depth >= {:.3f}, "
        "&nonhydrostatic_nml nz_layers = {}"
    ).format(len(dz), sum(dz), sum(dz), len(dz))
    return header + "\n" + "\n".join(lines) + "\n"


def convert(path=None, text_interfaces=None):
    if text_interfaces is not None:
        dz = read_text_interfaces(text_interfaces)
    elif path is not None:
        dz = read_vgrid_netcdf(path)
    else:
        raise ValueError("need either a vgrid NetCDF path or --text-interfaces")
    return dz


# ---------------------------------------------------------------------------
# Self-test: synthetic vgrid, known answer, both input paths. Stdlib only
# except for the netCDF4-backed sub-test, which is skipped (not failed) when
# netCDF4 is not importable.
# ---------------------------------------------------------------------------


def _expected_dz(nk):
    """Synthetic OM3-like profile: fine near the surface, coarsening with
    depth -- same shape used by the Fortran regression test
    (tests/test_ocean_vgrid_om3.F90)."""
    return [10.0 + 2.0 * k for k in range(nk)]


def _close(a, b, tol=1.0e-6):
    return abs(a - b) <= tol * max(1.0, abs(b))


def self_test():
    fails = []

    def check(label, ok, detail=""):
        print("  {:<64} {}".format(label, "ok" if ok else "FAILED " + detail))
        if not ok:
            fails.append(label)

    nk = 75
    expected = _expected_dz(nk)
    model_interfaces = [0.0]
    for d in expected:
        model_interfaces.append(model_interfaces[-1] + d)

    # (1) interfaces_to_dz recovers the exact thicknesses.
    got = interfaces_to_dz(model_interfaces)
    check(
        "(1) interfaces_to_dz recovers {} known thicknesses".format(nk),
        len(got) == nk and all(_close(g, e) for g, e in zip(got, expected)),
    )

    # (2) stdlib-only text-interfaces path, written to a temp file.
    import os
    import tempfile

    fd, tmp_path = tempfile.mkstemp(suffix=".txt")
    try:
        with os.fdopen(fd, "w") as fh:
            fh.write(",\n".join("{:.6f}".format(v) for v in model_interfaces))
        got_text = read_text_interfaces(tmp_path)
        check(
            "(2) --text-interfaces path (stdlib only, no netCDF4)",
            len(got_text) == nk and all(_close(g, e) for g, e in zip(got_text, expected)),
        )
    finally:
        os.remove(tmp_path)

    # (3) format_nml_block: sanity on structure/content, surface-first order
    #     preserved verbatim (no flip -- that's the kernel's job).
    block = format_nml_block(expected)
    check(
        '(3a) block selects vcoord_type = "zstar" + z_fixed_profile = "list"',
        'vcoord_type     = "zstar"' in block and 'z_fixed_profile = "list"' in block,
    )
    check(
        "(3b) first printed value is the surface layer (10.0 m), not the bed",
        "10.000000" in block.split("z_fixed_dz")[1].split("\n")[0],
    )
    check(
        "(3c) nz_layers hint matches the input count ({})".format(nk),
        "nz_layers = {}".format(nk) in block,
    )

    # (4) netCDF4-backed path: build a real ocean_vgrid.nc with a supergrid
    #     `zeta` (2*NK+1 values, odd-index entries interpolated -- their
    #     exact value is irrelevant since only zeta[::2] is read) and
    #     separately a `dz`-only file. Skipped (not failed) with no netCDF4.
    try:
        import netCDF4
        import numpy as np
    except ImportError:
        print("  (4) netCDF4 path ................................ SKIPPED "
              "(netCDF4/numpy not importable)")
        netCDF4 = None

    if netCDF4 is not None:
        fd, nc_path = tempfile.mkstemp(suffix=".nc")
        os.close(fd)
        try:
            supergrid = [0.0] * (2 * nk + 1)
            for k, zi in enumerate(model_interfaces):
                supergrid[2 * k] = zi
            for k in range(nk):
                supergrid[2 * k + 1] = 0.5 * (supergrid[2 * k] + supergrid[2 * k + 2])
            with netCDF4.Dataset(nc_path, "w") as ds:
                ds.createDimension("nkp1x2", len(supergrid))
                var = ds.createVariable("zeta", "f8", ("nkp1x2",))
                var[:] = np.array(supergrid, dtype="f8")
            got_nc = read_vgrid_netcdf(nc_path)
            check(
                "(4a) netCDF4 'zeta' (supergrid) path recovers the thicknesses",
                len(got_nc) == nk and all(_close(g, e) for g, e in zip(got_nc, expected)),
            )
        finally:
            os.remove(nc_path)

        fd, nc_path2 = tempfile.mkstemp(suffix=".nc")
        os.close(fd)
        try:
            with netCDF4.Dataset(nc_path2, "w") as ds:
                ds.createDimension("nk", nk)
                var = ds.createVariable("dz", "f8", ("nk",))
                var[:] = np.array(expected, dtype="f8")
            got_dz = read_vgrid_netcdf(nc_path2)
            check(
                "(4b) netCDF4 plain 'dz' variable path is taken verbatim",
                len(got_dz) == nk and all(_close(g, e) for g, e in zip(got_dz, expected)),
            )
        finally:
            os.remove(nc_path2)

    # (5) a non-monotone interface list is refused.
    bad = list(model_interfaces)
    bad[3] = bad[2] - 1.0
    refused = False
    try:
        interfaces_to_dz(bad)
    except ValueError:
        refused = True
    check("(5) non-monotone interfaces refused loudly", refused)

    print(
        "\nvgrid_to_dz.py self-test: {}".format(
            "ALL PASS" if not fails else "{} FAILED".format(len(fails))
        )
    )
    return 1 if fails else 0


def build_parser():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("vgrid", nargs="?", help="MOM6 ocean_vgrid.nc (requires netCDF4)")
    p.add_argument(
        "--text-interfaces",
        metavar="FILE",
        help="stdlib-only fallback: text file of NK+1 interface depths",
    )
    p.add_argument(
        "--self-test",
        action="store_true",
        help="run the built-in self-test (synthetic vgrid, known answer); no model is run",
    )
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.self_test:
        print("vgrid_to_dz.py self-test -- no model is run\n")
        return self_test()
    if not args.vgrid and not args.text_interfaces:
        build_parser().print_help()
        return 2
    try:
        dz = convert(path=args.vgrid, text_interfaces=args.text_interfaces)
    except (ValueError, RuntimeError) as exc:
        print("error: {}".format(exc), file=sys.stderr)
        return 1
    sys.stdout.write(format_nml_block(dz))
    return 0


if __name__ == "__main__":
    sys.exit(main())
