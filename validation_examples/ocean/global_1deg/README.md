# Global 1° tripolar ocean, unforced — roundabout's first global run

`global_1deg_unforced.nml` runs the MOM6 **OM_1deg** geometry (360 × 320
tracer points, bipolar Arctic cap, periodic in x, north tripolar fold) from
a **WOA13 January** ocean at rest, with **no forcing** — no wind, no heat or
freshwater flux, no sea ice — for one simulated year.

What it tests is the dynamical core on a real global ocean: that it holds a
realistic stratified state over real topography for a year without NaN,
with bounded energy, and with mass, salt and heat closed to round-off. It
is the reduced form of the v0.1.0 exit criterion E8 (the forced version
needs bulk formulae and an atmospheric-state reader that do not exist yet).
It is **not** a climate: with nothing driving it, the flow is the
geostrophic adjustment of the WOA density field plus what the grid and the
closures make of it.

| | |
|---|---|
| Grid | OM_1deg `ocean_hgrid.nc` supergrid (`grid_config = "supergrid"`), `west/east = "periodic"`, `north = "tripolar_fold"`, `south = "wall"` |
| Bathymetry | OM_1deg `topog.nc` with MOM6's own OM_1deg limits: wet cells raised to `MINIMUM_DEPTH` 9.5 m, capped at `MAXIMUM_DEPTH` 6500 m (`MASKING_DEPTH = 0`) |
| Vertical | 50 `z_fixed` levels, tanh-stretched from 2 m at the surface to 258 m at depth (`z_fixed_profile = "tanh"`), partial-step faces closed (`zfixed_closed_faces`) |
| Initial state | WOA13 decav January potential temperature + salinity, at rest, η = 0 |
| Time step | `dt = 1800 s` (OM_1deg `DT`), `pred_corr`, barotropic `n_inner` from the per-wet-cell CFL (32) |
| Physics | Wright EOS, FV-MOM6 PGF, Sadourny energy-conserving Coriolis, Smagorinsky Laplacian (0.15, floor 2000 m²/s) + biharmonic (0.06) with `bound_kh`, quadratic bottom drag distributed over `hbbl = 10 m` (OM_1deg `HBBL`, `DRAG_BG_VEL`), KPP + PP81, `maxvel = 6 m/s` (OM_1deg `MAXVEL`) |
| Output | daily: SSH, and the top-10 m means of u, v, T, S (one conservative `z_fixed` output level), single precision |
| Cost | 11 s per simulated day on one V100 (a year in 67 min), 10.4 GB of device memory |

## 1. Inputs (not in the repository)

```bash
export RDB_DATA_DIR=/somewhere/with/4GB            # re-downloadable data only
python3 tools/fetch_om1deg.py                      # -> $RDB_DATA_DIR/OM_1deg/
python3 tools/om1deg_prepare_inputs.py             # bathy_om1deg.nc + ic_woa13_jan.nc
```

`fetch_om1deg.py` streams `OM_1deg.tgz` (9 MB) and `obs.woa13.tgz`
(1.06 GB, only the two January T/S files are kept) from GFDL's public
MOM6-testing ftp tree, checks sizes and SHA-256, and is idempotent.
`om1deg_prepare_inputs.py` writes the model-grid bathymetry and the WOA13
initial condition (nearest neighbour horizontally, missing values filled
level by level from the nearest WOA cell with data at the same depth; the
model interpolates linearly in depth). Both are standard-library Python;
nothing to install. Takes about 30 s.

## 2. The reference run — the `rdb` executable

The namelist reads its three input files from `./INPUT/`:

```bash
mkdir run && cd run
ln -s "$RDB_DATA_DIR/OM_1deg" INPUT
cp /path/to/roundabout/validation_examples/ocean/global_1deg/global_1deg_unforced.nml .
CUDA_VISIBLE_DEVICES=0 /path/to/build/rdb global_1deg_unforced.nml > run.log
```

(or put absolute paths in `supergrid_file`, `bathymetry_file` and
`&ocean_zinit_nml file`). The console prints, every simulated day,
`[stats] Day … En … MaxCFL …` and the mass / salt / heat budget lines
(`Error` = the relative closure residual, `out` = the tracked boundary
flux, which is round-off in this closed domain).

## 3. The same run through the Python interface

`run_global_1deg.py` builds the same configuration knob by knob with the
typed `rdb.Config` API (and checks it against the namelist before it
starts), drives it with `rdb.Model(config).step(...)`, reads the
diagnostic manager's in-memory buffers every day
(`model.diagnostic("SSH")`, `"u"`, `"v"` — no NetCDF round trip) and
renders the movie frames as it goes:

```bash
export RDB_DATA_DIR=... RDB_LIB=/path/to/build/librdb_core.so CUDA_VISIBLE_DEVICES=0
python3 validation_examples/ocean/global_1deg/run_global_1deg.py --days 365 --out global_1deg_py
# -> global_1deg_py/global_1deg.mp4, .gif, daily.txt, frames/, and the diag file
```

The Python-driven run is the same computation as the executable's: its
diagnostic file is bit-identical to the reference run's (checked on the
year run, every record of SSH, u, v, T and S).

## 4. The movie from the executable's output

```bash
module load netcdf-c   # nccopy flattens the NetCDF-4 diag file for the stdlib reader
python3 validation_examples/ocean/global_1deg/global_movie.py \
    run/output/global_1deg_rank_000000.nc movie_frames --data "$RDB_DATA_DIR/OM_1deg"
```

`global_movie.py` is standard-library Python: a nearest-cell lookup from
the tripolar grid onto a regular 1/3° lon-lat image (the Arctic cap needs
no special case), anchor-table colour maps (inferno for speed, a diverging
map for SSH), a built-in bitmap font, PPM frames, then `ffmpeg` to MP4
(H.264, yuv420p) and GIF.

## 5. What the year shows

Measured on one V100 (nvfortran 26.5, `-gpu=cc70,mem:separate`), 2026-09-24:
365 days, 17 520 steps, 4013 s wall (11 s per simulated day), 10.4 GB of
device memory.

| day | En (m²/s²) | MaxCFL | Mass Error | Salt Error | Heat Error | max 3-D \|u\|,\|v\| (m/s) and where |
|---:|---:|---:|---:|---:|---:|---|
| 1 | 2.532e-04 | 0.031 | -1.80e-14 | -4.4e-16 | -2.3e-16 | 2.00, Strait of Gibraltar, ~550 m |
| 10 | 3.556e-04 | 0.112 | -1.80e-13 | -3.9e-15 | -3.2e-15 | 2.97, Celebes Sea trench, 4800 m |
| 30 | 5.416e-04 | 0.160 | -5.41e-13 | -1.2e-14 | -1.0e-14 | 4.23, Celebes Sea trench |
| 60 | 5.777e-04 | 0.156 | -1.08e-12 | -2.4e-14 | -2.1e-14 | 4.00, Celebes Sea trench |
| 90 | 5.747e-04 | 0.135 | -1.62e-12 | -3.7e-14 | -3.1e-14 | 3.46, Celebes Sea trench |
| 180 | 5.375e-04 | 0.091 | -3.25e-12 | -7.4e-14 | -6.1e-14 | 2.38, Celebes Sea trench |
| 240 | 4.996e-04 | 0.069 | -4.33e-12 | -1.0e-13 | -8.0e-14 | 1.75, Celebes Sea trench |
| 300 | 4.631e-04 | 0.059 | -5.41e-12 | -1.3e-13 | -9.8e-14 | 1.45, Celebes Sea trench |
| 365 | 4.345e-04 | 0.039 | -6.58e-12 | -1.5e-13 | -1.2e-13 | 1.12, Strait of Gibraltar |

* **No NaN, no CFL truncation, no positive-definite-limiter event, and the
  6 m/s `maxvel` clamp never reached** (3-D maximum over the year 4.33 m/s,
  day 38).
* **Energy is bounded.** En rises from rest over the first two months as the
  WOA density field adjusts geostrophically, peaks at 5.79e-04 m²/s² on
  day 65, and then decays slowly (4.35e-04 at day 365): with no forcing the
  closures spin it down.
* **Budgets close to round-off, linearly.** The relative residuals grow at
  a constant rate — mass −1.80e-14 per day, salt −4.1e-16, heat −3.4e-16 —
  i.e. round-off accumulating, not a leak; the tracked boundary fluxes
  (`out`) stay at round-off size in this closed domain (mass ≤ 1e2 kg of
  1.4e21, heat ≤ 7e10 J of 5.0e21).
* The fastest water all year is the Celebes Sea overflow (§6) and the
  Gibraltar exchange (1–2 m/s through a one-cell, 600 m channel — the
  Mediterranean outflow, the right order of magnitude). The surface
  (top-10 m daily mean) never exceeds 0.50 m/s: the equatorial current
  system and the Antarctic Circumpolar Current's fronts carry it, and the
  movie shows them spinning up and slowly decaying.
* **The Florida Straits jet is gone.** With uniform 130 m layers (the
  scoping probe) a 9.5 m cell next to 587 m cells in the Straits carried
  6 m/s by day 2 — with bed-only or HBBL drag alike. With the stretched
  profile the same cells are 2–5 m layers like their neighbours: the
  Straits stayed below 0.7 m/s over the first 12 days of the probe, and
  over the year they never carry the domain maximum.
* **The Python-driven run is bit-identical to the executable's**: every
  value of every record of SSH, T, S, u and v in the two diagnostic files
  (43.5 M values per field, ghosts included) is equal, and so are En and
  the mass drift each day.

## 6. Known limits of this configuration

* **Unforced, and 1°.** The surface circulation is the adjustment of the
  WOA January density field, slowly spinning down; there is no wind-driven
  gyre and no seasonal cycle. The movie is a picture of the dynamical core
  holding a global ocean, not of the ocean.
* **Viscosity is below what the grid wants in the western boundary
  layers.** The configure step says so: the Munk layer
  `(nu_h/β)^(1/3)` spans 0.4 cells with the 2000 m²/s floor. MOM6's OM_1deg
  uses a 2-D background (`KH_BG_2D`) file for this, which roundabout does
  not read yet.
* **One-cell trenches and sills.** At 1° some deep trenches are one cell
  wide and some sills sit much deeper than the real ones. Where a real sill
  keeps a basin warm (the Celebes Sea, 3.3 °C at depth behind a sill that
  the 1° topography puts at ~2700 m) the colder Pacific water overflows
  into it from the start and runs down a one-cell trench as a narrow bed
  jet (4.3 m/s at 4800 m by day 38, saturated by the bottom drag, then
  decaying to 0.9 m/s by day 365). It is the model
  geometry disagreeing with the observed state, not an instability, and it
  is the fastest flow in the run. With bed-only drag (see the next point)
  the same jet reaches the 6 m/s clamp by day 19.
* **`z_fixed` + bed-only bottom drag gives no bottom drag.** The bed-only
  mode (`hbbl = 0`) drags layer `k = 1`, which is an inert filler in every
  column shallower than the deepest level; configure warns. This namelist
  uses `hbbl = 10` (MOM6's value), which reaches the live bottom layer.
* **No channel list, no `CHANNEL_DRAG` file, no GM/MEKE, KPP instead of
  EPBL, no neutral diffusion, no sea ice.** None of them is needed for a
  stable year; all of them matter for a realistic one.
