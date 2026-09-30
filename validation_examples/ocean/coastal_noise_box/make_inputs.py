#!/usr/bin/env python3
"""Write bathy.nc + ic.nc for coastal_noise_box.nml (stdlib only, NetCDF-3).

Geometry (1/4 degree, 160 x 80 cells, lon 0..40 E, lat 66..46 S):
  * south coast: land south of lat_c(lon) = -62.5 + 2.5 sin(2 pi lon / 40),
    a staircase coast at every angle (one period, so the periodic seam is
    seamless).  --profile shelf (default): depth 400 m at the coast rising
    linearly to 4000 m over 300 km.  --profile wall: flat 4000 m, vertical
    staircase walls.  --profile rough: 20 m at the coast over 150 km, times a
    deterministic +-25 % grid-scale roughness (the rough-bathymetry control).
  * island: land within 120 km of (20 E, 55 S), same slope rule within 250 km.
IC: T(lat, d) = 2 + 8 exp(-d/800) + 2 tanh((lat + 54)/1.5) exp(-d/1000) degC,
S = 35 psu, on a 23-level source axis (zinit interpolates to the layers).

usage: make_inputs.py [--outdir DIR] [--profile shelf|wall|rough] [--res 0.25]
"""
import argparse
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "tools"))
from make_forcing_nc import build_cdf1  # noqa: E402

LON0, LAT0, LX, LY, H = 0.0, -66.0, 40.0, 20.0, 4000.0


def depth(lon, lat, i, j, profile):
    latc = LAT0 + 3.5 + 2.5 * math.sin(2 * math.pi * (lon - LON0) / LX)
    dist = (lat - latc) * 111.2
    noise = ((i * 73856093) ^ (j * 19349663)) % 1000 / 1000.0
    ilon, ilat = LON0 + 0.5 * LX, LAT0 + 11.0
    r = math.hypot((lon - ilon) * 111.2 * math.cos(math.radians(ilat)), (lat - ilat) * 111.2)
    if dist <= 0 or r < 120.0:
        return 0.0
    if profile == "wall":
        return H
    if profile == "rough":
        d = min(H, 20.0 + (H - 20.0) * dist / 150.0)
        if r < 220.0:
            d = min(d, 20.0 + (H - 20.0) * (r - 120.0) / 100.0)
        return min(H, d * (0.75 + 0.5 * noise))
    d = min(H, 400.0 + (H - 400.0) * dist / 300.0)
    if r < 250.0:
        d = min(d, 400.0 + (H - 400.0) * (r - 120.0) / 130.0)
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--profile", default="shelf", choices=["shelf", "wall", "rough"])
    ap.add_argument("--res", type=float, default=0.25)
    a = ap.parse_args()
    nx, ny = int(round(LX / a.res)), int(round(LY / a.res))
    b = []
    for j in range(ny):
        for i in range(nx):
            b.append(depth(LON0 + (i + 0.5) * a.res, LAT0 + (j + 0.5) * a.res, i, j, a.profile))
    with open(os.path.join(a.outdir, "bathy.nc"), "wb") as fh:
        fh.write(build_cdf1([("y", ny), ("x", nx)], [
            {"name": "b", "dimids": [0, 1], "data": b,
             "attrs": {"units": "m", "long_name": "bottom depth, positive down"}}]))
    z = [0.0, 5, 10, 20, 30, 50, 75, 100, 150, 200, 300, 400, 500, 700, 900,
         1200, 1500, 2000, 2500, 3000, 3500, 4000, 4500]
    T, S = [], []
    for d in z:
        for j in range(ny):
            lat = LAT0 + (j + 0.5) * a.res
            t = 2 + 8 * math.exp(-d / 800.0) + 2 * math.tanh((lat + 54.0) / 1.5) * math.exp(-d / 1000.0)
            T.extend([t] * nx)
            S.extend([35.0] * nx)
    with open(os.path.join(a.outdir, "ic.nc"), "wb") as fh:
        fh.write(build_cdf1([("z_src", len(z)), ("y", ny), ("x", nx)], [
            {"name": "z_src", "dimids": [0], "attrs": {"units": "m", "positive": "down"}, "data": z},
            {"name": "temp", "dimids": [0, 1, 2], "attrs": {"units": "degC"}, "data": T},
            {"name": "salt", "dimids": [0, 1, 2], "attrs": {"units": "psu"}, "data": S}]))
    print("wrote bathy.nc + ic.nc (%d x %d, profile=%s)" % (nx, ny, a.profile))


main()
