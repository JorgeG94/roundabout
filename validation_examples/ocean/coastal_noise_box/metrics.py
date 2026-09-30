#!/usr/bin/env python3
"""Coastal-noise metrics for coastal_noise_box (surface-layer corner zeta).

usage: metrics.py RUNDIR [--res 0.25] [--sponge-rows N] [--no-periodic-x]
       [--k LAYER]
Needs `ncdump` on PATH (load the NetCDF toolchain first); stdlib only.

Reads RUNDIR/output/restart_rank_000000.nc via ncdump (stdlib only),
rebuilds the dynamics' circulation-form corner relative vorticity on the
chosen layer (default the surface, k=nz), masked free-slip by wet_q, and
prints:
  rms_int   rms zeta at corners >= 8 cells from land, outside the sponge
  rms_coast rms zeta at wet corners with a land T-cell within 2 cells
  rms_spg   rms zeta in the north sponge rows (not coastal)
  ratio     rms_coast / rms_int
  nyq_*     2dx-ness: <(d2x^2 + d2y^2)/2> / <zeta^2>, d2 = (z+ - 2z + z-)/4
            (1.0 = pure +-alternating, 0.375 = white noise, ~0 = smooth)
  prof      rms zeta by distance-to-land 1..8
"""
import math
import os
import re
import subprocess
import sys



def ncvars(path, names):
    out = subprocess.run(["ncdump", "-v", ",".join(names), path],
                         capture_output=True, text=True, check=True).stdout
    data = out.split("data:", 1)[1]
    res = {}
    for name in names:
        m = re.search(r"\n\s*%s\s*=\s*(.*?);" % re.escape(name), data, re.S)
        txt = m.group(1).replace("\n", " ")
        res[name] = [float(t) if t.strip() != "_" else float("nan")
                     for t in txt.split(",")]
    return res


def dims(path):
    out = subprocess.run(["ncdump", "-h", path],
                         capture_output=True, text=True, check=True).stdout
    d = dict((m.group(1), int(m.group(2)))
             for m in re.finditer(r"\s(\w+) = (\d+) ;", out))
    return d


def main():
    a = sys.argv[1:]
    run = a[0]
    opt = {"--sponge-rows": -1, "--lon0": 0.0, "--lat0": -66.0, "--res": 0.25,
           "--ng": 3, "--k": -1, "--cart": 0.0}
    flags = {"--periodic-x"}
    i = 1
    while i < len(a):
        if a[i] == "--no-periodic-x":
            flags.discard("--periodic-x"); i += 1
        else:
            opt[a[i]] = float(a[i + 1]); i += 2
    opt["--dlon"] = opt["--dlat"] = opt["--res"]
    if opt["--sponge-rows"] < 0:  # the 5-degree band + 1 degree margin
        opt["--sponge-rows"] = round(6.0 / opt["--res"])
    path = os.path.join(run, "output", "restart_rank_000000.nc")
    d = dims(path)
    v = ncvars(path, ["ml_u_face_x_layer", "ml_v_face_y_layer", "ml_h_layer"])
    NX, NY, NZ = d["x4"], d["y4"], d["z4"]
    ng = int(opt["--ng"])
    k = NZ - 1 if opt["--k"] < 0 else int(opt["--k"]) - 1
    U = v["ml_u_face_x_layer"]; V = v["ml_v_face_y_layer"]; H = v["ml_h_layer"]

    def u(i, j):  # i in 0..NX, j in 0..NY-1  (face i = west face of T(i))
        return U[(k * NY + j) * (NX + 1) + i]

    def vv(i, j):
        return V[(k * (NY + 1) + j) * NX + i]

    def htot(i, j):
        return sum(H[(kk * NY + j) * NX + i] for kk in range(NZ))

    wet = [[1 if htot(i, j) > 2.0 else 0 for i in range(NX)] for j in range(NY)]
    # physical-domain view: ghost cells beyond a wall are land
    px = "--periodic-x" in flags
    py = "--periodic-y" in flags
    ni, nj = NX - 2 * ng, NY - 2 * ng

    def wetT(i, j):  # physical indices 0..ni-1, 0..nj-1
        if px:
            i %= ni
        if py:
            j %= nj
        if i < 0 or i >= ni or j < 0 or j >= nj:
            return 0
        return wet[j + ng][i + ng]

    R = 6.371e6
    dlon = math.radians(opt["--dlon"]); dlat = math.radians(opt["--dlat"])
    cart = opt["--cart"]

    def lat_T(j):
        return math.radians(opt["--lat0"] + (j + 0.5) * opt["--dlat"])

    def lat_q(j):
        return math.radians(opt["--lat0"] + j * opt["--dlat"])

    def dxCu(j):
        return cart if cart else R * math.cos(lat_T(j)) * dlon

    def dyCv():
        return cart if cart else R * dlat

    def areaBu(j):
        return cart * cart if cart else R * R * math.cos(lat_q(j)) * dlon * dlat

    # corner (i,j) = SW corner of T(i,j), physical i in 0..ni-1 (+1), j 0..nj
    Z = {}
    for j in range(1, nj):
        for i in range(0, ni + (0 if px else 1)):
            if not px and (i == 0 or i == ni):
                continue
            wq = wetT(i - 1, j - 1) * wetT(i, j - 1) * wetT(i - 1, j) * wetT(i, j)
            if not wq:
                continue
            I = i + ng; J = j + ng
            circ = ((vv(I, J) - vv(I - 1, J)) * dyCv()
                    - (u(I, J) * dxCu(j) - u(I, J - 1) * dxCu(j - 1)))
            Z[(i % ni if px else i, j)] = circ / areaBu(j)

    def landdist(i, j, maxd=10):
        for dd in range(1, maxd + 1):
            for jj in range(j - dd, j + dd):
                for ii in range(i - dd, i + dd):
                    if not wetT(ii, jj):
                        return dd
        return maxd + 1

    spg = int(opt["--sponge-rows"])
    reg = {"int": [], "coast": [], "spg": []}
    prof = {}
    dist = {}
    for (i, j), z in Z.items():
        dd = landdist(i, j)
        dist[(i, j)] = dd
        prof.setdefault(min(dd, 9), []).append(z)
        insp = spg > 0 and j >= nj - spg
        if dd <= 2:
            reg["coast"].append((i, j))
        elif insp:
            reg["spg"].append((i, j))
        elif dd >= 8:
            reg["int"].append((i, j))

    def get(i, j):
        if px:
            i %= ni
        return Z.get((i, j))

    def rms(pts):
        return math.sqrt(sum(Z[p] ** 2 for p in pts) / max(1, len(pts)))

    def nyq(pts):
        num = 0.0; den = 0.0
        for (i, j) in pts:
            z = Z[(i, j)]
            terms = []
            a1, b1 = get(i - 1, j), get(i + 1, j)
            if a1 is not None and b1 is not None:
                terms.append(((a1 - 2 * z + b1) / 4) ** 2)
            a2, b2 = get(i, j - 1), get(i, j + 1)
            if a2 is not None and b2 is not None:
                terms.append(((a2 - 2 * z + b2) / 4) ** 2)
            if terms:
                num += sum(terms) / len(terms); den += z * z
        return num / den if den > 0 else float("nan")

    ri, rc, rs = rms(reg["int"]), rms(reg["coast"]), rms(reg["spg"])
    print("%-28s rms_int %.3e rms_coast %.3e rms_spg %.3e ratio %.2f spg/int %.2f "
          "nyq_int %.2f nyq_coast %.2f nyq_spg %.2f  n=%d/%d/%d" % (
              os.path.basename(run.rstrip("/")), ri, rc, rs, rc / ri if ri else 0,
              rs / ri if ri else 0,
              nyq(reg["int"]), nyq(reg["coast"]), nyq(reg["spg"]) if reg["spg"] else 0,
              len(reg["int"]), len(reg["coast"]), len(reg["spg"])))
    print("   prof(dist 1..9+): " + " ".join(
        "%d:%.2e" % (dd, math.sqrt(sum(z * z for z in prof[dd]) / len(prof[dd])))
        for dd in sorted(prof)))


main()
