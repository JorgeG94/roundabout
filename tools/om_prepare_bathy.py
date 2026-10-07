#!/usr/bin/env python3
"""Arbitrary MOM6 supergrid + topog -> MOM6-limited bathymetry, generic grid size.

Generalises the bathymetry-limiting half of ``tools/om1deg_prepare_inputs.py``
(and duplicates, in Python, the rule ``tools/om_topo_limit.f90`` applies in
Fortran) to ANY MOM6 ``ocean_hgrid.nc`` + ``topog.nc`` pair, not just
OM_1deg's. It exists because ``om1deg_prepare_inputs.py`` is deliberately
standard-library-only (a hand-rolled NetCDF-3 classic reader/writer — see its
own docstring), and real ACCESS-OM3 (vk83) grids such as the 25 km
JRA55-do config are shipped as NetCDF-4/HDF5, which that reader cannot open
at all. This tool uses ``netCDF4`` + ``numpy`` instead (NOT the standard
library) specifically so it can open either container format transparently.
Both are available on Gadi via::

    module use /g/data/xp65/public/modules
    module load conda/analysis3

and locally via ``python3.10`` (see project memory: "python3.10 has
netCDF4").  Per CLAUDE.md, never ``pip install`` either package — if they
are missing, that is a signal this tool is being run in the wrong
environment, not a reason to install.

Rule applied (MOM6 ``limit_topography``, ``MASKING_DEPTH`` taken as the
land/sea cutoff): every cell with ``depth > masking_depth`` is wet, raised
to at least ``min_depth`` and capped at ``max_depth``; every other cell is
land, ``depth = 0``. Modifies depth/wet in a FRESH output file (never the
source in place — the source may be a read-only /g/data mount). The
bathymetry variable is auto-detected the same way ``rdb_bathymetry``'s
reader falls back (``depth``, then ``elevation``, then ``b``); a ``wet``
mask, if present in the source, is rewritten consistently with the limited
depth, and written if absent.

``--hgrid`` is an optional cross-check, not a transform: if given, this
tool reads the supergrid's node dims (``nxp``/``nyp``, MOM6 standard) and
FAILS LOUD if the implied T-point grid ``((nxp-1)/2, (nyp-1)/2)`` does not
match ``--topog``'s own horizontal dims — the same class of silent
grid-size mismatch the formula-bathymetry gotchas in CLAUDE.md warn about,
caught here before it ever reaches the model.

Usage::

    python3 tools/om_prepare_bathy.py --topog topog.nc --out bathy_limited.nc \\
        --min-depth 9.5 --max-depth 6500.0 --masking-depth 0.0 \\
        [--hgrid ocean_hgrid.nc]

Exit status 0 on success; non-zero with a message on stderr on any failure
(missing input, dimension mismatch, no bathymetry variable found).
"""

from __future__ import annotations

import argparse
import sys

try:
    import numpy as np
except ImportError as exc:  # pragma: no cover - environment problem, not a bug
    sys.stderr.write(
        "om_prepare_bathy.py needs numpy (see this file's docstring for how "
        f"to get it without pip install): {exc}\n"
    )
    sys.exit(1)

try:
    import netCDF4
except ImportError as exc:  # pragma: no cover - environment problem, not a bug
    sys.stderr.write(
        "om_prepare_bathy.py needs netCDF4 (see this file's docstring for "
        f"how to get it without pip install): {exc}\n"
    )
    sys.exit(1)

BATHY_VAR_FALLBACKS = ("depth", "elevation", "b")


def _find_bathy_var(ds):
    for name in BATHY_VAR_FALLBACKS:
        if name in ds.variables:
            return name
    raise SystemExit(
        f"no bathymetry variable found (tried: {', '.join(BATHY_VAR_FALLBACKS)})"
    )


def _check_hgrid(hgrid_path, nx, ny):
    """Fail loud if the supergrid's implied T-point grid != (nx, ny)."""
    with netCDF4.Dataset(hgrid_path, "r") as ds:
        # MOM6 supergrid: node dims nxp/nyp = 2*nx+1, 2*ny+1.
        dim_names = set(ds.dimensions)
        if "nxp" in dim_names and "nyp" in dim_names:
            nxp = len(ds.dimensions["nxp"])
            nyp = len(ds.dimensions["nyp"])
        else:
            # Fall back to the x/y coordinate variable's own shape.
            if "x" not in ds.variables:
                raise SystemExit(
                    f"{hgrid_path}: no nxp/nyp dims and no 'x' variable to "
                    "infer the supergrid shape from"
                )
            nyp, nxp = ds.variables["x"].shape
        if (nxp - 1) % 2 != 0 or (nyp - 1) % 2 != 0:
            raise SystemExit(
                f"{hgrid_path}: supergrid node dims ({nxp}, {nyp}) are not "
                "odd -- not a MOM6 supergrid"
            )
        hgrid_nx, hgrid_ny = (nxp - 1) // 2, (nyp - 1) // 2
        if (hgrid_nx, hgrid_ny) != (nx, ny):
            raise SystemExit(
                f"{hgrid_path} implies a ({hgrid_nx} x {hgrid_ny}) T-point "
                f"grid but --topog is ({nx} x {ny}) -- grid-size mismatch, "
                "refusing to guess (CLAUDE.md: formula-bathymetry-length-"
                "scale class of bug)"
            )


def limit_topography(depth, wet_in, min_depth, max_depth, masking_depth):
    """MOM6 ``limit_topography``: wet = depth > masking_depth; clamp in [min,max]."""
    wet = depth > masking_depth
    out = np.where(wet, np.clip(depth, min_depth, max_depth), 0.0)
    n_wet = int(wet.sum())
    n_land = int((~wet).sum())
    n_raised = int(np.count_nonzero(wet & (depth < min_depth)))
    n_capped = int(np.count_nonzero(wet & (depth > max_depth)))
    if wet_in is not None:
        # Report (not enforce) disagreement with a source wet mask -- the
        # depth-derived mask above is authoritative, matching om_topo_limit.f90.
        n_disagree = int(np.count_nonzero(wet != (wet_in > 0)))
        if n_disagree:
            sys.stderr.write(
                f"note: {n_disagree} cells where the recomputed wet mask "
                "disagrees with the source 'wet' variable (depth-derived "
                "mask wins)\n"
            )
    return out, wet, n_wet, n_land, n_raised, n_capped


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--topog", required=True, help="source topog.nc (depth[, wet])")
    ap.add_argument("--out", required=True, help="output path, fresh file")
    ap.add_argument("--hgrid", default=None, help="optional supergrid cross-check")
    ap.add_argument("--min-depth", type=float, default=9.5, dest="min_depth")
    ap.add_argument("--max-depth", type=float, default=6500.0, dest="max_depth")
    ap.add_argument("--masking-depth", type=float, default=0.0, dest="masking_depth")
    ap.add_argument(
        "--depth-var",
        default=None,
        help="override bathymetry variable name (default: autodetect depth/elevation/b)",
    )
    args = ap.parse_args()

    with netCDF4.Dataset(args.topog, "r") as src:
        depth_var = args.depth_var or _find_bathy_var(src)
        depth = np.array(src.variables[depth_var][:], dtype=np.float64)
        ny, nx = depth.shape
        wet_in = np.array(src.variables["wet"][:]) if "wet" in src.variables else None

    if args.hgrid:
        _check_hgrid(args.hgrid, nx, ny)

    limited, wet, n_wet, n_land, n_raised, n_capped = limit_topography(
        depth, wet_in, args.min_depth, args.max_depth, args.masking_depth
    )

    with netCDF4.Dataset(args.out, "w", format="NETCDF4") as dst:
        dst.createDimension("y", ny)
        dst.createDimension("x", nx)
        dvar = dst.createVariable("depth", "f8", ("y", "x"))
        dvar[:] = limited
        dvar.units = "m"
        dvar.long_name = "bathymetry depth, positive down, MOM6-limited"
        wvar = dst.createVariable("wet", "i4", ("y", "x"))
        wvar[:] = wet.astype(np.int32)
        wvar.long_name = "1 = ocean, 0 = land (depth-derived)"
        dst.source_file = args.topog
        dst.min_depth = args.min_depth
        dst.max_depth = args.max_depth
        dst.masking_depth = args.masking_depth

    sys.stdout.write(
        f"{args.out}: {nx} x {ny}, wet={n_wet} land={n_land} "
        f"raised={n_raised} capped={n_capped} "
        f"(min={args.min_depth} max={args.max_depth} masking={args.masking_depth})\n"
    )


if __name__ == "__main__":
    main()
