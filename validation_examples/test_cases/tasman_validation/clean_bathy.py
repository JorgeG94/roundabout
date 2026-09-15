#!/usr/bin/env python
"""Aggressive C-grid topology cleanup for a Roundabout ocean bathymetry file.

make_bathy.py already removes the worst degenerate cells (wet cells with
>=3 land neighbours).  Real coastlines additionally leave SHALLOW,
SEMI-ENCLOSED embayments connected to the open ocean only through a
1-cell-wide neck — the static land mask gives these almost no flux path,
and the wind-forced bottom layer there goes unstable (the Tasmania day-2
NaN was localised to exactly such a 100 m pocket).

This post-processor produces a SECOND, more aggressively cleaned bathy:
  1. iterate to convergence:
       a. land-fill any wet cell with >=3 of 4 edge-neighbours land,
       b. land-fill any 1-cell-wide channel (land on N&S or on E&W),
  2. keep only the LARGEST 4-connected wet region (the open ocean);
     fill every smaller disconnected basin / lake.
Surviving ocean-cell DEPTHS are never modified — topology only.

Usage:  ./clean_bathy.py in_bathy.nc out_bathy.nc
"""
import sys
import numpy as np
import netCDF4 as nc
from scipy import ndimage

LAND_DEPTH_THRESHOLD = 2.0
LAND_DEPTH = 0.0


def degenerate_count(wet):
    land = ~wet
    nland = np.zeros(wet.shape, dtype=int)
    nland[1:, :] += land[:-1, :]; nland[:-1, :] += land[1:, :]
    nland[:, 1:] += land[:, :-1]; nland[:, :-1] += land[:, 1:]
    nland[0, :] += 1; nland[-1, :] += 1; nland[:, 0] += 1; nland[:, -1] += 1
    ew = np.zeros(wet.shape, bool)
    ew[1:-1, :] = wet[1:-1, :] & land[:-2, :] & land[2:, :]
    ns = np.zeros(wet.shape, bool)
    ns[:, 1:-1] = wet[:, 1:-1] & land[:, :-2] & land[:, 2:]
    return int((wet & (nland >= 3)).sum()), int(ew.sum() + ns.sum())


def clean(b):
    wet = b >= LAND_DEPTH_THRESHOLD
    for _ in range(200):
        land = ~wet
        nland = np.zeros(wet.shape, dtype=int)
        nland[1:, :] += land[:-1, :]; nland[:-1, :] += land[1:, :]
        nland[:, 1:] += land[:, :-1]; nland[:, :-1] += land[:, 1:]
        nland[0, :] += 1; nland[-1, :] += 1; nland[:, 0] += 1; nland[:, -1] += 1
        deg = wet & (nland >= 3)
        # 1-cell-wide channels (land on opposite sides)
        ch = np.zeros(wet.shape, bool)
        ch[1:-1, :] |= wet[1:-1, :] & land[:-2, :] & land[2:, :]
        ch[:, 1:-1] |= wet[:, 1:-1] & land[:, :-2] & land[:, 2:]
        kill = deg | ch
        if not kill.any():
            break
        wet &= ~kill
    # largest connected wet region only (4-connectivity)
    lbl, n = ndimage.label(wet, structure=np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]]))
    if n > 1:
        sizes = ndimage.sum(np.ones_like(lbl), lbl, index=np.arange(1, n + 1))
        keep = 1 + int(np.argmax(sizes))
        wet = lbl == keep
    out = np.where(wet, b, LAND_DEPTH)
    return out, wet


def main():
    src, dst = sys.argv[1], sys.argv[2]
    d = nc.Dataset(src)
    b = np.array(d.variables["b"][:])
    wet0 = b >= LAND_DEPTH_THRESHOLD
    d3a, dchan = degenerate_count(wet0)
    out, wet = clean(b)
    d3b, dchan2 = degenerate_count(wet)
    print(f"in : land={(~wet0).mean()*100:.1f}%  3-nbr-degenerate={d3a}  1-cell-channels={dchan}")
    print(f"out: land={(~wet).mean()*100:.1f}%  3-nbr-degenerate={d3b}  1-cell-channels={dchan2}  "
          f"filled={int((wet0 & ~wet).sum())} cells")
    w = nc.Dataset(dst, "w")
    w.createDimension("y", b.shape[0]); w.createDimension("x", b.shape[1])
    v = w.createVariable("b", "f8", ("y", "x")); v[:] = out
    v.convention = "b = depth positive-down (Roundabout); topology-cleaned"
    w.close()
    print(f"wrote {dst}")


if __name__ == "__main__":
    main()
