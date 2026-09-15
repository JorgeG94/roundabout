#!/usr/bin/env python3
"""Surface relative-vorticity (zeta/f) movie for the baroclinic channel.

Works for both baroclinic_2layer.nml and baroclinic_15layer.nml output
(reads dx, f, nghost from args; surface = top z-level of the diag file).

Usage:
    python animate_baroclinic.py NETCDF [OUT.mp4] [--dx 2000] [--f 1e-4]
                                 [--ng 3] [--fps 30]
Defaults: OUT=baroclinic.mp4, dx=2000 m, f=1e-4 s^-1, ng=3, fps=30.
Falls back to .gif if ffmpeg is unavailable.
"""
import argparse
import numpy as np
import netCDF4
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation


def main():
    p = argparse.ArgumentParser()
    p.add_argument("netcdf")
    p.add_argument("out", nargs="?", default="baroclinic.mp4")
    p.add_argument("--dx", type=float, default=2000.0)
    p.add_argument("--f", type=float, default=1.0e-4)
    p.add_argument("--ng", type=int, default=3)
    p.add_argument("--fps", type=int, default=30)
    a = p.parse_args()

    d = netCDF4.Dataset(a.netcdf)
    u, v = d.variables["u"], d.variables["v"]   # (t, z, y, x)
    ng, k = a.ng, -1                            # surface = top z-level
    nt = u.shape[0]
    tvar = d.variables.get("time_u")
    days = (np.asarray(tvar[:]) / 86400.0) if tvar is not None else np.arange(nt)

    def zf(fr):
        U = np.asarray(u[fr, k, ng:-ng, ng:-ng])
        V = np.asarray(v[fr, k, ng:-ng, ng:-ng])
        return (np.gradient(V, a.dx, axis=1) - np.gradient(U, a.dx, axis=0)) / a.f

    ny, nx = zf(0).shape
    vmax = np.percentile(np.abs(zf(int(0.85 * nt))), 99.0)
    fig, ax = plt.subplots(figsize=(11, 5.5 * ny / nx + 1))
    im = ax.imshow(zf(0), origin="lower", cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                   aspect="auto", extent=[0, nx * a.dx / 1e3, 0, ny * a.dx / 1e3])
    ax.set_xlabel("x [km] (periodic)")
    ax.set_ylabel("y [km]")
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cb.set_label(r"surface $\zeta / f$")
    ttl = ax.set_title("")
    fig.tight_layout()

    def upd(fr):
        im.set_data(zf(fr))
        ttl.set_text(f"Baroclinic channel — surface $\\zeta/f$ — day {days[fr]:.1f}")
        return im, ttl

    ani = animation.FuncAnimation(fig, upd, frames=nt, blit=False)
    out = a.out
    try:
        ani.save(out, writer=animation.FFMpegWriter(fps=a.fps, bitrate=6000), dpi=130)
    except Exception:
        out = out.rsplit(".", 1)[0] + ".gif"
        ani.save(out, writer="pillow", fps=a.fps)
    print("saved", out, f"({nt} frames)")


if __name__ == "__main__":
    main()
