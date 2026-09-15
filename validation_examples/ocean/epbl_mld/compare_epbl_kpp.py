#!/usr/bin/env python3
"""EPBL-vs-KPP mixed-layer-depth comparison for the epbl_mld basin runs.

Reads the per-rank diag NetCDF from up to three runs (kpp_basin,
epbl_basin, epbl_lt_basin), computes an apples-to-apples
density-threshold MLD from the T/S output for every run, overlays the
EPBL scheme's own MLD_EPBL diagnostic where available, and writes
`epbl_mld_comparison.png` + a summary table to stdout.

MLD criterion: depth where sigma(z) exceeds the surface value by
DRHO_CRIT = 0.03 kg/m3 (de Boyer Montegut et al. 2004 convention),
linearly interpolated between layer centres.  Density from the same
linear EOS the runs use: rho = rho0 - alpha_T*(T - T_ref) (salinity is
uniform in this setup so the haline term drops).

Output-layout assumptions (verified against the diag manager):
  * file <out_dir>/<name>_rank_000000.nc, serial run;
  * `temperature` shaped (time, z, y, x) or (time, x, y, z) depending
    on the writer — handled by locating the z axis by length;
  * vertical index follows the model convention k=1 = bed, so the
    LAST z index is the surface layer;
  * uniform 20 m sigma layers on the flat 600 m bottom.

Run from the repo root after the three model runs:
    python3 validation_examples/ocean/epbl_mld/compare_epbl_kpp.py
Optionally pass --base-dir to point at where the out_* directories
were created (defaults to the current working directory).
"""

import argparse
import glob
import os
import sys

import numpy as np

try:
    from netCDF4 import Dataset
except ImportError:
    sys.exit("compare_epbl_kpp.py needs the netCDF4 python module")

ALPHA_T = 0.2        # kg/m3/K — matches &ocean_ic_nml alpha_T
DRHO_CRIT = 0.03     # kg/m3 MLD threshold
MAX_DEPTH = 600.0    # m, flat bottom
NZ = 30
DZ = MAX_DEPTH / NZ

RUNS = [
    ("KPP", "out_kpp_basin", "kpp_basin"),
    ("EPBL", "out_epbl_basin", "epbl_basin"),
    ("EPBL+LT", "out_epbl_lt_basin", "epbl_lt_basin"),
]


def find_file(base_dir, out_dir, name):
    pats = [
        os.path.join(base_dir, out_dir, f"{name}_rank_000000.nc"),
        os.path.join(base_dir, out_dir, f"{name}*rank*0.nc"),
        os.path.join(base_dir, out_dir, "*.nc"),
    ]
    for pat in pats:
        hits = sorted(glob.glob(pat))
        if hits:
            return hits[0]
    return None


def get_var_tzyx(ds, name):
    """Return variable as (time, z, y, x) with z index 0 = bed."""
    if name not in ds.variables:
        return None
    v = ds.variables[name][:]
    arr = np.ma.filled(v, np.nan)
    if arr.ndim == 3:  # (time, y, x) 2D field
        return arr
    if arr.ndim != 4:
        raise ValueError(f"{name}: unexpected rank {arr.ndim}")
    # time is axis 0 by construction; find the z axis among 1..3 by
    # matching NZ (grid is 44x40 + ghosts, never 30).
    zax = [a for a in (1, 2, 3) if arr.shape[a] == NZ]
    if len(zax) != 1:
        raise ValueError(f"{name}: cannot identify z axis in shape {arr.shape}")
    arr = np.moveaxis(arr, zax[0], 1)
    return arr


def mld_from_T(T_tzyx):
    """Density-threshold MLD (m) per (time, y, x).

    z index 0 = bed, -1 = surface.  Scans downward from the surface.
    """
    nt, nz, ny, nx = T_tzyx.shape
    # layer-centre depth below surface for z index k (bed-up):
    z_centre = (nz - np.arange(nz) - 0.5) * DZ
    rho = -ALPHA_T * T_tzyx  # constant offsets cancel in the threshold
    rho_sfc = rho[:, -1:, :, :]
    excess = rho - rho_sfc  # >= 0 below the mixed layer
    mld = np.full((nt, ny, nx), MAX_DEPTH)
    # walk downward: surface index nz-1, deeper = smaller index
    for t in range(nt):
        ex = excess[t]
        found = np.zeros((ny, nx), dtype=bool)
        for k in range(nz - 2, -1, -1):
            newly = (~found) & (ex[k] >= DRHO_CRIT)
            if newly.any():
                # linear interp between layer centres k+1 (above) and k
                e_above = ex[k + 1][newly]
                e_here = ex[k][newly]
                frac = (DRHO_CRIT - e_above) / np.maximum(e_here - e_above, 1e-30)
                mld[t][newly] = z_centre[k + 1] + frac * (z_centre[k] - z_centre[k + 1])
                found |= newly
        mld[t][~found] = MAX_DEPTH
    return mld


def interior(a2d, ng=2):
    return a2d[..., ng:-ng, ng:-ng]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base-dir", default=".",
                    help="directory containing the out_* run directories")
    ap.add_argument("--out", default="epbl_mld_comparison.png")
    args = ap.parse_args()

    results = {}
    for label, out_dir, name in RUNS:
        path = find_file(args.base_dir, out_dir, name)
        if path is None:
            print(f"[skip] {label}: no NetCDF under {out_dir}/")
            continue
        ds = Dataset(path)
        T = get_var_tzyx(ds, "temperature")
        if T is None:
            print(f"[skip] {label}: no 'temperature' in {path}")
            continue
        entry = {"mld": mld_from_T(T), "path": path}
        mld_diag = get_var_tzyx(ds, "MLD_EPBL")
        if mld_diag is not None:
            entry["mld_diag"] = mld_diag
        results[label] = entry
        ds.close()
        print(f"[ok]   {label}: {path}  frames={entry['mld'].shape[0]}")

    if not results:
        sys.exit("no runs found — run the three namelists first")

    # ---- summary table ----
    print("\n  run        day-1 mean   final mean   final max   (interior MLD, m)")
    for label, e in results.items():
        m = interior(e["mld"])
        line = (f"  {label:<9}  {np.nanmean(m[0]):10.1f}  "
                f"{np.nanmean(m[-1]):10.1f}  {np.nanmax(m[-1]):9.1f}")
        if "mld_diag" in e:
            d = interior(e["mld_diag"])
            line += f"   | MLD_EPBL final mean {np.nanmean(d[-1]):.1f}"
        print(line)

    # ---- plots ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n_maps = len(results)
    fig = plt.figure(figsize=(5 + 4 * n_maps, 8))

    # (1) basin-mean MLD time series
    ax = fig.add_subplot(2, 1, 1)
    for label, e in results.items():
        m = interior(e["mld"])
        days = np.arange(1, m.shape[0] + 1)
        ax.plot(days, np.nanmean(m, axis=(1, 2)), label=f"{label} (Δρ criterion)")
        if "mld_diag" in e:
            d = interior(e["mld_diag"])
            ax.plot(days, np.nanmean(d, axis=(1, 2)), "--",
                    label=f"{label} (MLD_EPBL diag)")
    ax.set_xlabel("day")
    ax.set_ylabel("basin-mean MLD (m)")
    ax.invert_yaxis()
    ax.legend()
    ax.grid(alpha=0.3)
    ax.set_title("EPBL vs KPP — basin-mean mixed-layer depth")

    # (2) final-frame MLD maps
    vmax = max(np.nanmax(interior(e["mld"][-1])) for e in results.values())
    for i, (label, e) in enumerate(results.items()):
        axm = fig.add_subplot(2, n_maps, n_maps + 1 + i)
        pc = axm.pcolormesh(interior(e["mld"][-1]), vmin=0, vmax=vmax,
                            cmap="viridis_r")
        axm.set_title(f"{label} — final MLD (m)")
        fig.colorbar(pc, ax=axm, shrink=0.85)

    fig.tight_layout()
    fig.savefig(args.out, dpi=140)
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
