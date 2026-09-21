#!/usr/bin/env python3
"""Generate an analytic bathymetry NetCDF for `&ocean_topo_nml topo_config="file"`.

Why this exists.  The formula bathymetries the solver ships (`flat`,
`spoon`, `seamount`, `island`, `double_drake`, `neverworld2`,
`isomip_plus`) do not include the two geometries the vertical-coordinate
stability matrix needs, and both are defined by ONE number that has to be
dialled exactly:

  * **`slope`** -- a gentle constant-gradient bed.  Every terrain-following
    coordinate's pressure-gradient truncation scales as the CUBE of the
    per-face interface offset (`a_peak = N^2 * de^3 / (6 dx Hbar)`), so a
    bed whose gradient is CONSTANT is the only geometry on which that
    scaling can be read off cleanly: every wet-wet face carries the same
    `de` and the same `rx0`.
  * **`rx0_ladder`** -- a two-level shelf/trough cross-section whose
    terrain-following STIFFNESS

        rx0 = |H_a - H_b| / (H_a + H_b)

    (Beckmann & Haidvogel 1993, J. Phys. Oceanogr. 23, 1736-1753, section
    2c, who write it as `r = |dh| / (2 hbar)`; the same statement as
    Haney's 1991 hydrostatic-consistency condition) takes an EXACT
    prescribed value at the step.  `rdb_ocean_stability_audit` bounds that
    number at 0.2 and warns above it; the ladder walks it from 0.1 to 0.8
    so the matrix can say where each coordinate family stops being
    trustworthy instead of asserting the bound and hoping.

    Given a target `rx0` and a mean depth `Hbar`, the two levels follow
    from the definition in closed form:

        H_deep    = Hbar * (1 + rx0)
        H_shallow = Hbar * (1 - rx0)

    which reproduces `rx0` exactly at the single face joining them.

Writes NetCDF-3 *classic* (CDF-1) by hand through `make_forcing_nc.build_cdf1`
-- this box has no `netCDF4`, no `numpy` and no `scipy`, and the project
forbids installing packages.

The loader (`src/io/rdb_bathymetry.F90`) reads `b(x, y)` as bottom depth
POSITIVE-DOWN, sized `nx_phys x ny_phys`, and fills the ghost rows itself
by constant extrapolation.

Dimension order: netCDF stores dimension ids in C order (first varies
SLOWEST) while the netCDF-Fortran API takes them in Fortran order, so a
variable Fortran must see as `b(x, y)` is written with the dimid list
`[y, x]` and data laid out with x varying fastest.

Usage:
    make_bathy_nc.py OUT.nc --nx 48 --ny 6 --profile slope \
                     --deep 1000 --shallow 500
    make_bathy_nc.py OUT.nc --nx 48 --ny 6 --profile rx0_ladder \
                     --hbar 750 --rx0 0.6
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_forcing_nc import build_cdf1  # noqa: E402


def profile_slope(nx, ny, deep, shallow):
    """Constant-gradient bed, `deep` at i = 1 and `shallow` at i = nx.

    Every interior face carries the SAME depth step `(deep - shallow)/(nx-1)`
    and therefore the same interface offset and the same `rx0` (which drifts
    only through the `1/(H_a + H_b)` denominator).  That is what makes this
    the clean geometry for a truncation-scaling measurement.
    """
    if nx < 2:
        raise ValueError("slope needs nx >= 2")
    step = (deep - shallow) / float(nx - 1)
    col = [deep - step * i for i in range(nx)]
    return [col[i] for _ in range(ny) for i in range(nx)], step


def profile_rx0_ladder(nx, ny, hbar, rx0):
    """Two-level shelf/trough with ONE step of exactly the requested `rx0`.

    The deep half occupies `i < nx/2`, the shallow half the rest, so there
    is a single wet-wet face at the join carrying
    `rx0 = |H_d - H_s| / (H_d + H_s)`; every other face is flat and carries
    `rx0 = 0`.  Isolating the stiffness on ONE face is deliberate: a ladder
    that ramps would convolve the step's response with the ramp's.
    """
    if not 0.0 < rx0 < 1.0:
        raise ValueError("rx0 must be in (0, 1)")
    h_deep = hbar * (1.0 + rx0)
    h_shallow = hbar * (1.0 - rx0)
    half = nx // 2
    col = [h_deep if i < half else h_shallow for i in range(nx)]
    return ([col[i] for _ in range(ny) for i in range(nx)],
            (h_deep, h_shallow, abs(h_deep - h_shallow) / (h_deep + h_shallow)))


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("out")
    p.add_argument("--nx", type=int, required=True)
    p.add_argument("--ny", type=int, required=True)
    p.add_argument("--profile", choices=("slope", "rx0_ladder"), required=True)
    p.add_argument("--deep", type=float, default=1000.0)
    p.add_argument("--shallow", type=float, default=500.0)
    p.add_argument("--hbar", type=float, default=750.0)
    p.add_argument("--rx0", type=float, default=0.2)
    a = p.parse_args()

    if a.profile == "slope":
        data, step = profile_slope(a.nx, a.ny, a.deep, a.shallow)
        note = "constant-gradient bed, depth step {:.4f} m per face".format(step)
    else:
        data, (hd, hs, r) = profile_rx0_ladder(a.nx, a.ny, a.hbar, a.rx0)
        note = ("two-level shelf/trough, H_deep {:.3f} m, H_shallow {:.3f} m, "
                "rx0 {:.6f} at the single step".format(hd, hs, r))

    dims = [("y", a.ny), ("x", a.nx)]
    variables = [{
        "name": "b",
        "dimids": [0, 1],          # C order: y slowest, x fastest
        "attrs": {"units": "m", "long_name": "bottom depth, positive down",
                  "note": note},
        "data": data,
    }]
    with open(a.out, "wb") as fh:
        fh.write(build_cdf1(dims, variables))
    sys.stdout.write("{}: {} ({} x {})\n".format(a.out, note, a.nx, a.ny))
    return 0


if __name__ == "__main__":
    sys.exit(main())
