#!/usr/bin/env python3
"""M2 tidal validation — least-squares harmonic fit of ocean SSH output.

Fits a pure M2 sinusoid (period 12.4206 h) to the SSH time series in a Roundabout
ocean diagnostics NetCDF and reports, per probe cell, the M2 amplitude and the
fraction of variance the single constituent explains (R^2).  A short record
cannot resolve M2 by FFT (bin spacing is too coarse), so the harmonic fit is
the robust check: R^2 -> 1 means the signal IS M2.

Usage:
    python analyze_m2.py <output_rank_000000.nc> [--spinup-days N]

Expected (see README.md):
    body_tide_basin      : wall probe R^2 ~ 0.98, interior amp ~ 1 cm (equilibrium)
    seamount_tidal_obc   : every probe R^2 > 0.999, amp ~ 40-65 cm, W->E phase lag
"""
import argparse
import numpy as np
import netCDF4 as nc

M2_HOURS = 12.4206
OMEGA_M2 = 2.0 * np.pi / (M2_HOURS * 3600.0)  # rad/s


def fit_m2(t, series):
    """Return (amplitude_m, R^2, phase_rad) of a least-squares M2 fit."""
    g = np.column_stack([np.ones_like(t), np.cos(OMEGA_M2 * t), np.sin(OMEGA_M2 * t)])
    coef, *_ = np.linalg.lstsq(g, series, rcond=None)
    pred = g @ coef
    ss_tot = np.sum((series - series.mean()) ** 2)
    r2 = 1.0 - np.sum((series - pred) ** 2) / ss_tot if ss_tot > 0 else 0.0
    return np.hypot(coef[1], coef[2]), r2, np.arctan2(coef[2], coef[1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ncfile")
    ap.add_argument("--spinup-days", type=float, default=1.5,
                    help="discard this many days of spin-up before fitting")
    args = ap.parse_args()

    ds = nc.Dataset(args.ncfile)
    t = ds.variables["time_SSH"][:].astype(float)          # seconds since start
    ssh = ds.variables["SSH"][:].astype(float)             # (time, y, x)
    keep = t >= t[0] + args.spinup_days * 86400.0
    t, ssh = t[keep], ssh[keep]
    ny, nx = ssh.shape[1:]
    span_d = (t[-1] - t[0]) / 86400.0
    print(f"{args.ncfile}")
    print(f"  grid {ny}x{nx}, fit window {span_d:.1f} d, {len(t)} samples "
          f"(spin-up {args.spinup_days:.1f} d dropped)")

    # Report a mid-row west / centre / east transect + the max-variance cell.
    jm = ny // 2
    probes = {"west":   (jm, 1),
              "centre": (jm, nx // 2),
              "east":   (jm, nx - 2)}
    var = ssh.var(axis=0)
    pj, pi = np.unravel_index(np.argmax(var), var.shape)
    probes["max-var"] = (pj, pi)

    worst_r2 = 1.0
    for label, (j, i) in probes.items():
        amp, r2, ph = fit_m2(t, ssh[:, j, i])
        worst_r2 = min(worst_r2, r2)
        print(f"  {label:8s} (y={j:3d},x={i:3d}): M2 amp = {amp*100:7.2f} cm  "
              f"R^2 = {r2:.4f}  phase = {np.degrees(ph):7.1f} deg")

    verdict = ("M2 CONFIRMED" if worst_r2 > 0.95
               else "M2 clear" if worst_r2 > 0.90
               else "SIGNAL NOISY")
    print(f"  VERDICT: {verdict} (min R^2 = {worst_r2:.4f})")
    return 0 if worst_r2 > 0.90 else 1


if __name__ == "__main__":
    raise SystemExit(main())
