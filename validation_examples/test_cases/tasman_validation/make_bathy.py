#!/usr/bin/env python
"""
Build a Cartesian bathymetry NetCDF for the Tasman calibration run.

Pulls a GEBCO sub-area (lat/lon raster), projects to an equirectangular
Cartesian grid centred mid-domain, interpolates onto the model's
`(nx_phys, ny_phys)` resolution, and writes a NetCDF the Roundabout ocean
driver can ingest via `topo_config = "file"`.

Conventions:
  * Roundabout's `b` is bottom depth POSITIVE-DOWN (ocean = positive).
    GEBCO ships `elevation` positive-up, so we flip the sign and clip
    above zero to drop land cells onto a no-flow value.
  * Output file has dims `(x, y)` with variable `b(x, y)` (Roundabout's
    loader also accepts `elevation` / `depth` names).

Workflow:
  1. Download a GEBCO sub-area from https://download.gebco.net/
     (web form: pick the lat/lon box, get NetCDF back).  Default box
     below is the EAC corner — adjust if you want a bigger domain.
  2. Run this script: `python make_bathy.py --gebco gebco.nc --dx 5e3`
  3. Resulting `tasman_bathy_5km.nc` is what `tasman_5km.nml` references.

Land handling:
  * Cells with GEBCO elevation > 0 (subaerial) get `b = LAND_DEPTH` (0 m),
    which is below Roundabout's `LAND_DEPTH_THRESHOLD` (2 m), so the static
    land mask seeds `wet_mask = 0` there and walls them as free-slip
    continents (metric-zeroing + wet_q/wet_T) — no min-depth pinching.
  * Cells shallower than `MIN_OCEAN_DEPTH` get clamped to `MIN_OCEAN_DEPTH`
    to avoid divide-by-tiny issues in shallow shelf regions.  This clamp
    is applied to OCEAN cells ONLY, BEFORE land is stamped — otherwise the
    clamp would bump land up past the threshold and drown the coastline.

Dependencies: numpy, xarray, netCDF4 (for writing).  No pyproj
required for the equirectangular projection — the trig fits in 6 lines.
"""

import argparse
import numpy as np
import xarray as xr


# ---------------------------------------------------------------------
# Defaults — tweak via CLI args.  Box matches the GEBCO sub-area download
# `gebco_2026_n-32.5_s-41.5_w152.3_e163.7.nc` (southern Tasman / Bass
# Strait region, ~1010 × 1000 km centred on -37°S 158°E).
# ---------------------------------------------------------------------
DEFAULT_GEBCO = "gebco_2026_n-32.5_s-41.5_w152.3_e163.7.nc"
DEFAULT_LON_MIN = 152.3
DEFAULT_LON_MAX = 163.7
DEFAULT_LAT_MIN = -41.5
DEFAULT_LAT_MAX = -32.5
DEFAULT_LAT0 = -37.0   # projection latitude reference (mid-domain)
DEFAULT_DX_M = 2.0e3   # 2 km nominal resolution — eddy-resolving

LAND_DEPTH_THRESHOLD = 2.0  # m, Roundabout's wet_mask cutoff (rdb_constants)
LAND_DEPTH = 0.0          # m, subaerial cells; < LAND_DEPTH_THRESHOLD (2 m)
                          # so Roundabout's static land mask walls them as
                          # free-slip continents (no min-depth pinching).
MIN_OCEAN_DEPTH = 100.0   # m, shallowest ocean cell
                          # ↑ At 30 sigma layers, MIN_OCEAN_DEPTH / nz controls
                          #   the worst-case top-layer thickness.  100 m / 30
                          #   ≈ 3.3 m surface layer — wind stress τ/(ρ·h) stays
                          #   below ~3e-5 m/s² peak, stable under any plausible
                          #   τ_x.  Drop to 5 m and you can blow up in 6 h.
EARTH_RADIUS = 6371.0e3   # m


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gebco", default=DEFAULT_GEBCO,
                        help="Path to GEBCO NetCDF (lat/lon raster).")
    parser.add_argument("--out", default=None,
                        help="Output Cartesian bathymetry NetCDF "
                             "(default: tasman_bathy_<dx>km.nc).")
    parser.add_argument("--lon-min", type=float, default=DEFAULT_LON_MIN)
    parser.add_argument("--lon-max", type=float, default=DEFAULT_LON_MAX)
    parser.add_argument("--lat-min", type=float, default=DEFAULT_LAT_MIN)
    parser.add_argument("--lat-max", type=float, default=DEFAULT_LAT_MAX)
    parser.add_argument("--lat0", type=float, default=DEFAULT_LAT0,
                        help="Projection-centre latitude (deg).")
    parser.add_argument("--dx", type=float, default=DEFAULT_DX_M,
                        help="Cartesian resolution (m).")
    args = parser.parse_args()

    out_path = args.out
    if out_path is None:
        out_path = f"tasman_bathy_{int(args.dx/1e3)}km.nc"

    # --- Load + subset GEBCO ---
    ds = xr.open_dataset(args.gebco)
    # GEBCO uses 'lat', 'lon'.  Some sub-area cuts use 'latitude'/'longitude'.
    lat_name = "lat" if "lat" in ds.coords else "latitude"
    lon_name = "lon" if "lon" in ds.coords else "longitude"
    elev_name = "elevation" if "elevation" in ds.data_vars else "z"

    sub = ds.sel({lat_name: slice(args.lat_min, args.lat_max),
                  lon_name: slice(args.lon_min, args.lon_max)})
    lat = sub[lat_name].values
    lon = sub[lon_name].values
    elev = sub[elev_name].values  # shape (nlat, nlon), positive-up

    print(f"GEBCO subset: {len(lon)} × {len(lat)} cells, "
          f"elev range [{elev.min():.0f}, {elev.max():.0f}] m")

    # --- Equirectangular projection centred at lat0 ---
    cos_lat0 = np.cos(np.deg2rad(args.lat0))
    x_gebco = EARTH_RADIUS * cos_lat0 * np.deg2rad(lon - args.lon_min)
    y_gebco = EARTH_RADIUS * np.deg2rad(lat - args.lat_min)

    # --- Target Cartesian grid ---
    x_max = EARTH_RADIUS * cos_lat0 * np.deg2rad(args.lon_max - args.lon_min)
    y_max = EARTH_RADIUS * np.deg2rad(args.lat_max - args.lat_min)
    nx = int(round(x_max / args.dx))
    ny = int(round(y_max / args.dx))
    x_model = (np.arange(nx) + 0.5) * args.dx   # cell centres
    y_model = (np.arange(ny) + 0.5) * args.dx

    print(f"Cartesian grid: {nx} × {ny} cells at {args.dx/1e3} km, "
          f"covers {x_max/1e3:.0f} × {y_max/1e3:.0f} km")

    # --- Bilinear interpolation onto model grid ---
    # xarray's interp does bilinear by default.  Stage a temporary DataArray
    # on the projected coords so we don't fight lat/lon → Cartesian roundoff.
    da = xr.DataArray(
        elev,
        coords={"y": y_gebco, "x": x_gebco},
        dims=("y", "x"),
    )
    elev_model = da.interp(x=x_model, y=y_model, method="linear").values

    # --- Convert elevation (+up) → depth (+down) + handle land ---
    # Order matters: clamp the shallow-ocean floor FIRST (ocean cells only),
    # THEN stamp land.  Doing it the other way round bumps the land cells up
    # through MIN_OCEAN_DEPTH and drowns the coastline — Roundabout's static land
    # mask seeds `wet_mask` from `b >= LAND_DEPTH_THRESHOLD` (2 m), so land
    # MUST end up below that to be walled.
    depth = -elev_model            # ocean: positive; land: negative/zero
    is_ocean = depth > 0.0
    depth = np.where(is_ocean & (depth < MIN_OCEAN_DEPTH), MIN_OCEAN_DEPTH, depth)
    depth = np.where(is_ocean, depth, LAND_DEPTH)   # land last; LAND_DEPTH < 2 m ⇒ masked

    # --- C-grid topology cleanup ("remove lakes / edit topography") ---
    # A real coastline leaves numerically-degenerate single-cell features that
    # Roundabout's static land mask cannot represent: a wet cell with >=3 land
    # neighbours has all/most of its C-grid face metrics zeroed, so it has no
    # valid flux path -> continuity/PGF divide by a vanishing volume -> Inf
    # (vcoord-independent NaN on day ~2).  Iteratively convert any wet cell
    # with >=3 of its 4 edge-neighbours land to land, until none remain.  This
    # removes isolated cells, fills 1-cell bays, and erodes 1-cell channels /
    # peninsulas.  It is TOPOLOGY cleanup, NOT depth smoothing (depths of
    # surviving ocean cells are untouched) — the standard MOM6/ROMS step.
    wet = depth >= LAND_DEPTH_THRESHOLD
    for _ in range(50):
        land = ~wet
        nland = np.zeros(wet.shape, dtype=int)
        nland[1:, :] += land[:-1, :]; nland[:-1, :] += land[1:, :]
        nland[:, 1:] += land[:, :-1]; nland[:, :-1] += land[:, 1:]
        # domain-edge faces count as land (no neighbour across the wall)
        nland[0, :] += 1; nland[-1, :] += 1; nland[:, 0] += 1; nland[:, -1] += 1
        degenerate = wet & (nland >= 3)
        if not degenerate.any():
            break
        wet &= ~degenerate
        depth = np.where(degenerate, LAND_DEPTH, depth)
    n_removed = int((depth < LAND_DEPTH_THRESHOLD).sum() - (~is_ocean).sum())

    print(f"Model bathy: depth range [{depth.min():.1f}, "
          f"{depth.max():.1f}] m, "
          f"{(depth < LAND_DEPTH_THRESHOLD).sum()} land cells "
          f"({n_removed} added by topology cleanup), "
          f"{(depth >= LAND_DEPTH_THRESHOLD).sum()} ocean cells")

    # --- Write NetCDF ---
    # Roundabout's loader expects dims named `x`, `y` with variable `b`.
    # Storage order in Fortran is `b(x, y)` after the C/Fortran flip the
    # loader handles automatically.
    out_ds = xr.Dataset(
        {"b": (("y", "x"), depth.astype(np.float64))},
        coords={"x": x_model, "y": y_model},
        attrs={
            "title": "Tasman calibration bathymetry",
            "source": f"GEBCO subset {args.lon_min}-{args.lon_max}E "
                      f"{args.lat_min}-{args.lat_max}S",
            "projection": f"equirectangular, lat0={args.lat0}",
            "dx_m": args.dx,
            "land_depth_m": LAND_DEPTH,
            "min_ocean_depth_m": MIN_OCEAN_DEPTH,
            "convention": "b = depth positive-down (Roundabout)",
        },
    )
    out_ds.to_netcdf(out_path)
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
