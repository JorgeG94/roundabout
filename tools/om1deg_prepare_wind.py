#!/usr/bin/env python3
"""JRA55-do 10 m wind -> OM_1deg C-grid wind-stress file for ``global_1deg_wind.nml``.

Writes ``wind_jra55do_<year>_<avg>h.nc`` (default: daily means,
``$RDB_DATA_DIR/OM_1deg/wind_jra55do_1958_24h.nc``) in the layout the
``&ocean_dataovr_nml`` ``tau_x``/``tau_y`` tags read:

* ``taux(xf = 361, y = 320, time)`` on the u faces (the WEST face of T-cell
  ``i``, supergrid node ``(2i-1, 2j)``), ``tauy(x = 360, yf = 321, time)`` on
  the v faces (the SOUTH face, node ``(2i, 2j-1)``) — one face more than
  cells in the staggered direction, Fortran order, float32;
* ``time`` in ``days since <year>-01-01`` at the bin centres (0.5, 1.5, ...
  for daily bins), so the namelist needs ``t_offset = 0`` and
  ``cycle_period`` = the year's length.

Method, per 3-hourly JRA55-do record (``uas``/``vas``, point values on the
TL319 640x320 lon/Gaussian-lat grid):

1. bilinear interpolation of the east/north wind to each face;
2. quadratic bulk stress ``tau = rho_air * C_d(|U|) * |U| * (u_E, v_N)``,
   ``rho_air = 1.22 kg m-3``; ``C_d`` is the Large & Yeager (2004, NCAR Tech.
   Note TN-460, eq. 6a) neutral 10 m drag coefficient
   ``(2.7/U + 0.142 + 0.0764 U) * 1e-3`` with ``U >= 0.5 m/s`` (the NCAR/FMS
   floor), or ``--cd CONST``. No stability correction (the JRA 10 m wind is
   used as if it were the neutral wind), and the ocean surface velocity is
   ignored (absolute, not relative, wind);
3. rotation onto the grid axes by the grid angle ``a`` at the face:
   ``tau_i = cos(a) tau_E + sin(a) tau_N``, ``tau_j = -sin(a) tau_E + cos(a)
   tau_N`` (the ``ocean_metrics_t%angle_dx`` convention, counter-clockwise
   from east). ``a`` is computed from the supergrid node positions, NOT read
   from the mosaic's ``angle_dx``: in OM_1deg's Arctic cap that variable is
   the heading in degree space, without the ``cos(lat)`` metric (10.4 against
   the true 32.7 degrees at 73 N — what MOM6's ``GRID_ROTATION_ANGLE_BUGS``
   is about), and it is one-sided on the periodic seam column and the fold
   row. The two agree south of the cap, where both are zero;
4. trapezoidal mean over ``--avg-hours`` bins (a record on a bin edge
   counts half to each side; the ``padded`` JRA file carries the next year's
   00:00 record, so every daily bin has exactly 8 weights).

The stress is formed at 3-hourly resolution BEFORE averaging: averaging the
wind first and then applying the quadratic law would drop the synoptic
variance and underestimate the mean stress.

The number crunching (2920 records x 2 x 115 k faces, reading 3.4 GB of
compressed NetCDF-4) is done by a small Fortran helper,
``tools/om1deg_wind_regrid.f90``, which this script compiles with the system
Fortran compiler and ``nf-config`` (NetCDF-Fortran) and then runs — pure
Python over 670 M face evaluations would take hours; the helper takes about
a minute. Nothing to install.

    export RDB_DATA_DIR=...     # holds OM_1deg/ and JRA55do/
    python3 tools/fetch_om1deg.py --jra-wind
    python3 tools/om1deg_prepare_wind.py [--year 1958] [--avg-hours 24] [--cd ly04]
"""

import argparse
import calendar
import os
import shlex
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
JRA_TAG = "input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-4-0_gr"


def jra_file(jra_dir, var, year):
    """The padded GFDL file if present, else the plain input4MIPs name."""
    stem = f"{var}_{JRA_TAG}_{year}01010000-{year}12312100"
    for name in (stem + ".padded.nc", stem + ".nc"):
        path = os.path.join(jra_dir, name)
        if os.path.isfile(path):
            return path
    raise SystemExit(f"no JRA55-do {var} file for {year} in {jra_dir} "
                     "(run tools/fetch_om1deg.py --jra-wind)")


def nf_config(flag):
    exe = shutil.which("nf-config")
    if exe is None:
        raise SystemExit("nf-config (NetCDF-Fortran) not on PATH — load the NetCDF-Fortran "
                         "module of the compiler you build roundabout with")
    out = subprocess.run([exe, flag], check=True, capture_output=True, text=True).stdout
    return shlex.split(out)


def build_helper(build_dir, fc):
    """Compile tools/om1deg_wind_regrid.f90 into `build_dir`; return the binary path."""
    src = os.path.join(HERE, "om1deg_wind_regrid.f90")
    exe = os.path.join(build_dir, "om1deg_wind_regrid")
    cmd = [fc, "-O2", *nf_config("--fflags"), "-J", build_dir, src, "-o", exe,
           *nf_config("--flibs")]
    print("building helper:", " ".join(cmd))
    subprocess.run(cmd, check=True)
    return exe


def main():
    root = os.environ.get("RDB_DATA_DIR")
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--year", type=int, default=1958)
    ap.add_argument("--avg-hours", type=int, default=24,
                    help="averaging bin in hours, a divisor of 24 (default 24 = daily)")
    ap.add_argument("--cd", default="ly04",
                    help="'ly04' (Large & Yeager 2004 neutral 10 m) or a constant, e.g. 1.2e-3")
    ap.add_argument("--om1deg", default=None, help="OM_1deg directory (default $RDB_DATA_DIR/OM_1deg)")
    ap.add_argument("--jra", default=None, help="JRA55-do directory (default $RDB_DATA_DIR/JRA55do)")
    ap.add_argument("--out", default=None,
                    help="output file (default <om1deg>/wind_jra55do_<year>_<avg>h.nc)")
    ap.add_argument("--fc", default=os.environ.get("FC", "gfortran"), help="Fortran compiler")
    ap.add_argument("--build-dir", default=None,
                    help="where the helper is compiled (default: a temporary directory)")
    args = ap.parse_args()

    if (args.om1deg is None or args.jra is None) and not root:
        ap.error("set RDB_DATA_DIR or pass both --om1deg and --jra")
    om1deg = os.path.abspath(args.om1deg or os.path.join(root, "OM_1deg"))
    jra = os.path.abspath(args.jra or os.path.join(root, "JRA55do"))
    if args.cd != "ly04":
        float(args.cd)  # fail early on a typo
    out = os.path.abspath(args.out or os.path.join(
        om1deg, f"wind_jra55do_{args.year}_{args.avg_hours}h.nc"))
    hgrid = os.path.join(om1deg, "ocean_hgrid.nc")
    uas, vas = jra_file(jra, "uas", args.year), jra_file(jra, "vas", args.year)

    tmp = None
    build_dir = args.build_dir
    if build_dir is None:
        tmp = tempfile.TemporaryDirectory(prefix="om1deg_wind_")
        build_dir = tmp.name
    os.makedirs(build_dir, exist_ok=True)
    try:
        exe = build_helper(build_dir, args.fc)
        subprocess.run([exe, hgrid, uas, vas, out, str(args.year), args.cd,
                        str(args.avg_hours)], check=True)
    finally:
        if tmp is not None:
            tmp.cleanup()

    ndays = 366 if calendar.isleap(args.year) else 365
    print(f"\n{out}: {os.path.getsize(out) / 1e6:.1f} MB")
    print("namelist (&ocean_dataovr_nml): time_mode = \"cyclic\", "
          f"cycle_period = {ndays * 86400.0:.1f}, t_offset = 0.0, "
          "tau_x_var = \"taux\", tau_y_var = \"tauy\"")
    return 0


if __name__ == "__main__":
    sys.exit(main())
