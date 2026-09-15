#!/usr/bin/env python3
"""3-panel movie of the idealized-ACC channel showcase run.

Panels:
  1. SSH (m) + surface flow quivers  -> the zonal jet + standing meander
     downstream of the seamount ridge.
  2. MLD_EPBL (m)                     -> the EPBL active mixing-layer depth.
  3. Kd_KSHEAR (m^2/s) at the upper interface -> the kappa-shear interior
     diffusivity (log scale; lives in the wind-driven surface shear here).

Spherical channel: the x axis is longitude (deg E from lon_west), the y
axis is latitude (deg N).  The ridge crest sits mid-channel; the periodic
seam is at the left/right edges.

Usage:
    python animate_acc_channel.py
    python animate_acc_channel.py [NETCDF] [OUTPUT]
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import netCDF4
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import LogNorm

p = argparse.ArgumentParser()
p.add_argument("netcdf", nargs="?",
               default="./out_acc_channel/acc_channel_rank_000000.nc")
p.add_argument("out", nargs="?", default="acc_channel_3panel.mp4")
p.add_argument("--fps", type=int, default=8)
p.add_argument("--dpi", type=int, default=160)
p.add_argument("--nghost", type=int, default=3)
p.add_argument("--lon-west", type=float, default=0.0)
p.add_argument("--lat-south", type=float, default=-52.0)
p.add_argument("--dlon", type=float, default=0.5)
p.add_argument("--dlat", type=float, default=0.5)
p.add_argument("--quiver-skip", type=int, default=2,
               help="show 1-in-N velocity arrows over the SSH panel")
args = p.parse_args()

ds = netCDF4.Dataset(args.netcdf)
NG = args.nghost
t = ds.variables["time_SSH"][:]
SSH = ds.variables["SSH"][:, NG:-NG, NG:-NG]
u = ds.variables["u"][:, :, NG:-NG, NG:-NG]
v = ds.variables["v"][:, :, NG:-NG, NG:-NG]
mld = ds.variables["MLD_EPBL"][:, NG:-NG, NG:-NG]
kks = ds.variables["Kd_KSHEAR"][:, :, NG:-NG, NG:-NG]
ds.close()

nt = len(t)
nz = u.shape[1]
ny, nx = SSH.shape[-2:]
ks = nz - 1   # surface layer (k=1 bed -> k=nz surface)
kif = nz - 1  # upper interface for Kd

print(f"Loaded: {nt} snapshots, grid {ny}x{nx}x{nz}, run {t[-1]/86400:.0f} days")


def vmax(field, pct=99.5):
    return np.percentile(np.abs(field), pct)


ssh_vmax = max(vmax(SSH), 1e-3)
spd = np.sqrt(u[:, ks]**2 + v[:, ks]**2)
mld_max = max(np.percentile(mld, 99.5), 10.0)
mld_min = max(np.percentile(mld[mld > 0], 1) if np.any(mld > 0) else 1.0, 1.0)
kd_surf = kks[:, kif]
kd_floor = 1e-6
kd_cap = max(np.percentile(kd_surf[kd_surf > kd_floor], 99.5)
             if np.any(kd_surf > kd_floor) else 1e-3, 1e-4)

print(f"  SSH vmax = +/-{ssh_vmax:.3f} m")
print(f"  MLD range = [{mld_min:.0f}, {mld_max:.0f}] m")
print(f"  Kd_KSHEAR cap = {kd_cap:.2e} m^2/s")

# Geographic cell-centre coordinates (degrees).
lon = args.lon_west + (np.arange(nx) + 0.5) * args.dlon
lat = args.lat_south + (np.arange(ny) + 0.5) * args.dlat
Xc, Yc = np.meshgrid(lon, lat)

fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.4), constrained_layout=True)
fig.patch.set_facecolor("white")

# --- SSH + flow panel ---
ax = axes[0]
im_ssh = ax.pcolormesh(Xc, Yc, SSH[0], cmap="RdBu_r",
                       vmin=-ssh_vmax, vmax=ssh_vmax, shading="gouraud")
cb = plt.colorbar(im_ssh, ax=ax, fraction=0.046, pad=0.04)
cb.set_label("SSH (m)")
qsk = args.quiver_skip
qx = Xc[::qsk, ::qsk]; qy = Yc[::qsk, ::qsk]
quiver = ax.quiver(qx, qy, u[0, ks][::qsk, ::qsk], v[0, ks][::qsk, ::qsk],
                   color="k", scale=max(vmax(spd, 99) * 15, 1.0),
                   width=0.004, alpha=0.8)
ax.set_title("SSH + surface flow", fontsize=11)
ax.set_xlabel("lon (deg E)"); ax.set_ylabel("lat (deg N)")

# --- MLD panel ---
ax = axes[1]
im_mld = ax.pcolormesh(Xc, Yc, mld[0], cmap="viridis_r",
                       vmin=mld_min, vmax=mld_max, shading="gouraud")
cb = plt.colorbar(im_mld, ax=ax, fraction=0.046, pad=0.04)
cb.set_label("MLD_EPBL (m)")
ax.set_title("EPBL mixing-layer depth", fontsize=11)
ax.set_xlabel("lon (deg E)")

# --- Kd_KSHEAR panel (log scale) ---
ax = axes[2]
kd0 = np.clip(kd_surf[0], kd_floor, kd_cap)
im_kd = ax.pcolormesh(Xc, Yc, kd0, cmap="magma",
                      norm=LogNorm(vmin=kd_floor, vmax=kd_cap),
                      shading="gouraud")
cb = plt.colorbar(im_kd, ax=ax, fraction=0.046, pad=0.04)
cb.set_label("Kd_KSHEAR (m$^2$/s)")
ax.set_title("kappa-shear diffusivity (surface interface)", fontsize=11)
ax.set_xlabel("lon (deg E)")

day_text = fig.suptitle(f"Roundabout idealized ACC channel - day {t[0]/86400:.1f}",
                        fontsize=13, fontweight="bold")


def update(idx):
    im_ssh.set_array(SSH[idx].ravel())
    quiver.set_UVC(u[idx, ks][::qsk, ::qsk], v[idx, ks][::qsk, ::qsk])
    im_mld.set_array(mld[idx].ravel())
    im_kd.set_array(np.clip(kd_surf[idx], kd_floor, kd_cap).ravel())
    day_text.set_text(f"Roundabout idealized ACC channel - day {t[idx]/86400:5.1f}")
    return im_ssh, quiver, im_mld, im_kd, day_text


print(f"Rendering {nt} frames at {args.fps} fps to {args.out}...")
anim = animation.FuncAnimation(fig, update, frames=nt, interval=1000/args.fps,
                               blit=False)

try:
    writer = animation.FFMpegWriter(fps=args.fps, bitrate=2000,
                                    extra_args=["-crf", "18",
                                                "-preset", "slow",
                                                "-pix_fmt", "yuv420p"])
    anim.save(args.out, writer=writer, dpi=args.dpi)
    print(f"  -> wrote {args.out}")
except (FileNotFoundError, RuntimeError) as e:
    print(f"  ffmpeg failed ({e}); falling back to .gif")
    gif = Path(args.out).with_suffix(".gif")
    anim.save(gif, writer=animation.PillowWriter(fps=args.fps), dpi=120)
    print(f"  -> wrote {gif}")
