#!/usr/bin/env python3
"""Vertical tilt and surface-vorticity signature of an isolated vortex.

Reads the single-rank NetCDF diagnostics of a `balanced_vortex` run
(`vorticity_z`, `temperature`, optionally `SSH`; merge multi-rank output with
`tools/merge_output.py` first) and writes a text table:

  * per snapshot, the vortex CENTRE of the surface layer and of a deep
    reference layer (a same-sign-weighted centroid within a window of
    `--window` x R about the previous surface centre);
  * the TILT: the horizontal offset of the upper (surface) centre from the
    deep centre, per unit depth difference.  `tilt_deg` is the angle from the
    vertical, atan(|offset| / dz) (`slope` is the bare |offset| / dz, m per m: the tilt
    angle saturates toward 90 degrees once the offset exceeds the depth
    difference, so read `slope`); `azimuth_deg` is the direction the UPPER
    part is displaced toward, measured counter-clockwise from east (90 = north);
  * the SURFACE-VORTICITY signature: the extremum of the vortex-sign vorticity
    and its same-sign area integral (circulation, m^2/s) within the window.

Two centre definitions are reported (`--method` picks which one drives the
tilt columns):

  T      weights = anomaly of the layer temperature against the layer's
         far-field (outer rim) mean, vortex-sign only.  Defined at every
         depth the temperature anomaly reaches, so it is the robust choice
         for the deep layers.  (default)
  zeta   weights = same-sign relative vorticity.  At depth the vorticity
         signal is weak and may change sign; use for the surface layer and
         the near-surface layers.

Layer index 0 is the bed (the model's bottom-up convention).  Coordinates
are cartesian: cell (i, j) (0-based, ghosts included) sits at
x = (i - nghost + 1/2) dx.  Standard library + numpy + netCDF4 only.

Example:
    python3.10 analyze_vortex_tilt.py out_vortex_beta/vortex_beta_rank_000000.nc \\
        --dx 5000 --nghost 4 --depth 4000 --radius 50000 --tilt-depth 1500
"""
import argparse
import math
import sys

import numpy as np
import netCDF4 as nc


def load(path, name):
    ds = nc.Dataset(path)
    if name not in ds.variables:
        return None, None
    var = ds.variables[name]
    tname = [d for d in var.dimensions if d.startswith("time")][0]
    t = np.array(ds.variables[tname][:], dtype=float)
    arr = np.ma.filled(var[:], np.nan).astype(float)
    return t, arr


def centroid(weights, x, y, mask):
    """Weighted centroid and sum of weights over `mask` (2-d arrays)."""
    w = np.where(mask, weights, 0.0)
    sw = w.sum()
    if not np.isfinite(sw) or sw <= 0.0:
        return float("nan"), float("nan"), 0.0
    return (w * x).sum() / sw, (w * y).sum() / sw, sw


def rim_mean(field2d, ng, width=3):
    """Mean over a ring `width` cells wide just inside the physical edge."""
    ny, nx = field2d.shape
    phys = field2d[ng:ny - ng, ng:nx - ng]
    ring = np.ones(phys.shape, dtype=bool)
    ring[width:-width, width:-width] = False
    return float(np.nanmean(phys[ring]))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("file", help="single-rank diagnostics NetCDF file")
    ap.add_argument("--dx", type=float, required=True, help="grid spacing (m)")
    ap.add_argument("--dy", type=float, default=None, help="y spacing (m); default dx")
    ap.add_argument("--nghost", type=int, default=4)
    ap.add_argument("--depth", type=float, required=True, help="total depth H (m)")
    ap.add_argument("--radius", type=float, required=True, help="vortex radius R (m)")
    ap.add_argument("--window", type=float, default=2.5, help="window half-width, in R")
    ap.add_argument("--tilt-depth", type=float, default=1500.0,
                    help="depth (m) of the deep reference layer for the tilt")
    ap.add_argument("--method", choices=["T", "zeta"], default="T")
    ap.add_argument("--f0", type=float, default=1.0e-4)
    ap.add_argument("--cold", action="store_true",
                    help="cold-core / cyclone (default: infer from SSH, else warm)")
    ap.add_argument("--profile-day", type=float, default=None,
                    help="also print the centre-vs-depth profile at this day "
                         "(default: the last snapshot)")
    ap.add_argument("--png", default=None, help="write a PNG here if matplotlib exists")
    a = ap.parse_args()
    dy = a.dy or a.dx
    ng = a.nghost

    tz, zeta = load(a.file, "vorticity_z")
    tT, temp = load(a.file, "temperature")
    tS, ssh = load(a.file, "SSH")
    if zeta is None or temp is None:
        sys.exit("need vorticity_z and temperature in " + a.file)
    nt, nz, ny, nx = temp.shape
    dz = a.depth / nz
    depth_k = (nz - 1 - np.arange(nz) + 0.5) * dz     # k = 0 is the bed
    kd = int(np.argmin(np.abs(depth_k - a.tilt_depth)))
    ks = nz - 1
    xs = (np.arange(nx) - ng + 0.5) * a.dx
    ys = (np.arange(ny) - ng + 0.5) * dy
    X, Y = np.meshgrid(xs, ys)
    phys = np.zeros((ny, nx), dtype=bool)
    phys[ng:ny - ng, ng:nx - ng] = True
    cell_area = a.dx * dy

    # vortex sign: warm core (+1) = anticyclone
    if a.cold:
        warm = False
    elif ssh is not None:
        warm = float(np.nanmax(ssh[0][phys])) > -float(np.nanmin(ssh[0][phys]))
    else:
        warm = True
    s_T = 1.0 if warm else -1.0          # sign of the T anomaly
    s_z = -1.0 if warm else 1.0          # sign of the vorticity (f > 0)

    cx = 0.5 * (xs[ng] + xs[nx - ng - 1])      # domain centre = IC centre
    cy = 0.5 * (ys[ng] + ys[ny - ng - 1])
    win = a.window * a.radius

    def layer_centre(t_i, k, c0, method):
        m = phys & ((X - c0[0]) ** 2 + (Y - c0[1]) ** 2 <= win ** 2)
        if method == "zeta":
            w = np.maximum(0.0, s_z * zeta[t_i, k])
        else:
            anom = temp[t_i, k] - rim_mean(temp[t_i, k], ng)
            w = np.maximum(0.0, s_T * anom)
            peak = np.nanmax(np.where(m, w, 0.0))
            w = np.where(w >= 0.1 * peak, w, 0.0)      # ignore the faint tail
        return centroid(w, X, Y, m)

    print("# vortex: %s core (%s), R = %.0f km, window = %.0f km, method = %s"
          % ("warm" if warm else "cold", "anticyclone" if warm else "cyclone",
             a.radius / 1e3, win / 1e3, a.method))
    print("# surface layer k=%d (depth %.0f m); deep reference layer k=%d (depth %.0f m); "
          "dz_ref = %.0f m" % (ks, depth_k[ks], kd, depth_k[kd], depth_k[kd] - depth_k[ks]))
    print("# azimuth: direction of the UPPER centre relative to the deep centre, "
          "deg CCW from east (90 = north)")
    cols = ("day", "xs_km", "ys_km", "xd_km", "yd_km", "off_km", "tilt_deg",
            "azim_deg", "slope", "zeta_ext", "zeta/f0", "circ_m2s", "circ/circ0")
    print("# " + " ".join("%11s" % c for c in cols))

    rows = []
    c_prev = (cx, cy)
    circ0 = None
    for ti in range(nt):
        day = tz[ti] / 86400.0
        # track the surface centre by vorticity (stable, sharp)
        xs_c, ys_c, _ = layer_centre(ti, ks, c_prev, "zeta")
        if not math.isfinite(xs_c):
            xs_c, ys_c = c_prev
        c_prev = (xs_c, ys_c)
        cs = layer_centre(ti, ks, c_prev, a.method)
        cd = layer_centre(ti, kd, c_prev, a.method)
        mwin = phys & ((X - c_prev[0]) ** 2 + (Y - c_prev[1]) ** 2 <= win ** 2)
        zsur = s_z * zeta[ti, ks]
        zext = float(np.nanmax(np.where(mwin, zsur, -np.inf)))
        circ = float(np.nansum(np.where(mwin & (zsur > 0.0), zsur, 0.0)) * cell_area)
        if circ0 is None:
            circ0 = circ
        offx = cs[0] - cd[0]
        offy = cs[1] - cd[1]
        off = math.hypot(offx, offy) if math.isfinite(offx) and math.isfinite(offy) else float("nan")
        dzr = depth_k[kd] - depth_k[ks]
        dzr = abs(dzr) if dzr != 0 else float("nan")
        tilt = math.degrees(math.atan2(off, dzr)) if math.isfinite(off) else float("nan")
        azim = math.degrees(math.atan2(offy, offx)) % 360.0 if off > 0 else float("nan")
        slope = off / dzr if math.isfinite(off) else float("nan")   # m per m
        row = (day, cs[0] / 1e3, cs[1] / 1e3, cd[0] / 1e3, cd[1] / 1e3, off / 1e3,
               tilt, azim, slope, s_z * zext, s_z * zext / a.f0, circ, circ / circ0)
        rows.append(row)
        print("  " + " ".join("%11.4f" % v if abs(v) < 1e5 and (abs(v) >= 1e-3 or v == 0) else "%11.4e" % v
                              for v in row))

    # centre vs depth at one snapshot
    pi = nt - 1
    if a.profile_day is not None:
        pi = int(np.argmin(np.abs(tz / 86400.0 - a.profile_day)))
    print("# centre vs depth at day %.2f (method %s)" % (tz[pi] / 86400.0, a.method))
    print("# %4s %9s %9s %9s %12s" % ("k", "depth_m", "x_km", "y_km", "weight"))
    c0 = (rows[pi][1] * 1e3, rows[pi][2] * 1e3)
    for k in range(nz - 1, -1, -1):
        cxk, cyk, sw = layer_centre(pi, k, c0, a.method)
        print("  %4d %9.0f %9.2f %9.2f %12.4e" % (k, depth_k[k], cxk / 1e3, cyk / 1e3, sw))

    if a.png:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            print("# matplotlib not available: no PNG written")
            return
        r = np.array(rows)
        fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
        ax[0].plot(r[:, 1], r[:, 2], "o-", label="surface")
        ax[0].plot(r[:, 3], r[:, 4], "s-", label="deep (%.0f m)" % depth_k[kd])
        ax[0].set_xlabel("x (km)")
        ax[0].set_ylabel("y (km)")
        ax[0].legend()
        ax[0].set_aspect("equal", "datalim")
        ax[1].plot(r[:, 0], r[:, 8], "o-")
        ax[1].set_xlabel("day")
        ax[1].set_ylabel("tilt slope: offset / dz (m/m)")
        ax[2].plot(r[:, 0], r[:, 10], "o-", label="peak zeta / f0")
        ax[2].plot(r[:, 0], r[:, 12], "s-", label="circulation / initial")
        ax[2].set_xlabel("day")
        ax[2].legend()
        fig.tight_layout()
        fig.savefig(a.png, dpi=120)
        print("# wrote " + a.png)


if __name__ == "__main__":
    main()
