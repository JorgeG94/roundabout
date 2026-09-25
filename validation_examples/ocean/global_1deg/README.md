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

**Wind-forced variant:** `global_1deg_wind.nml` adds one year of JRA55-do
wind stress through `&ocean_dataovr_nml` and changes nothing else — §6.

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

## 1. Reproduce it — one command

Given a roundabout build (the `rdb` executable; for `--python` also
`librdb_core.so`, i.e. configured with `-DRDB_BUILD_SHARED=ON`):

```bash
validation_examples/ocean/global_1deg/reproduce.sh \
    --data-dir /somewhere/with/4GB --build-dir /path/to/build --gpu 0
```

That fetches the inputs, prepares them, runs the year, renders the movie and
checks the run against the committed reference. `--quick` runs 10 days
instead and checks them against the reference's first 10 days — about
two and a half minutes on one V100 including the movie, the way to check a
new machine or toolchain.

| option | |
|---|---|
| `--data-dir DIR` | data root; the inputs land in `DIR/OM_1deg` (default `$RDB_DATA_DIR`, else it stops and says so). About 4 GB, all re-downloadable. |
| `--build-dir DIR` | the build tree holding `rdb` and `librdb_core.so` (default `$RDB_BUILD_DIR`) |
| `--run-dir DIR` | where the run happens (default `./global_1deg_run`) |
| `--days N` / `--quick` | simulated days (default 365) / 10 days |
| `--python` | drive the run through the Python interface (`run_global_1deg.py`) instead of the executable |
| `--movie` / `--no-movie` | require / skip the movie (default: made when `ffmpeg` — and, for the executable path, NetCDF-C's `nccopy` — is on `PATH`, otherwise skipped with a note) |
| `--gpu N` | `CUDA_VISIBLE_DEVICES=N` |
| `--skip-fetch` | the data are already in place |
| `--no-check` | do not compare against the reference |

Both paths make the movie: the executable path renders it afterwards from
the diagnostic file (`global_movie.py`), the Python path in-process as the
run goes. It prints where everything went: the run directory holds
`global_1deg_unforced.nml` (with `t_end` set from `--days`), the `INPUT`
symlink to the data, `run.log` (the console), `output/` (the diagnostic
file; with `--python` also `daily.txt`, the frames and the movie),
`movie/` (the executable path's frames and movie) and `check.txt`.

**On a machine without NetCDF:** this case needs a build with
NetCDF-Fortran (and NetCDF-C) — the bathymetry and initial-condition
readers and the diagnostic file use it. Where the system has neither,
`tools/build_netcdf.sh` builds both from source against an existing HDF5;
see the docs on building NetCDF from source
(`docs/howto/deploy_without_netcdf.md`, both arriving with the
`feat/netcdf-bootstrap` branch). The scripts here are standard-library
Python and never need it (only the executable path's movie uses
NetCDF-C's `nccopy`).

## 2. The check against the reference

`reference_daily.csv` is the executable's year: En, MaxCFL and the mass /
salt / heat `Error` of every day, measured on one V100 with nvfortran 26.5
(`-gpu=cc70,mem:separate`, double precision); the header records the
commit. `check_against_reference.py RUN.log [--days N | --quick]` parses a
new run's console and prints a day-by-day table and PASS/FAIL (exit status
0/1). The executable's console carries everything; the Python driver's
`[py] day` line carries En and the mass drift only — the Python API does
not expose MaxCFL or the salt and heat totals — so a `--python` run is
checked on those two (the two runs are bit-identical, so the executable's
check covers the rest). The two kinds of number are checked differently,
because another toolchain does not produce the same bits and round-off
grows:

* **Energy — a relative band that widens with time.** En is a global
  integral of a smooth, large-scale adjustment, so it stays close across
  toolchains long after individual features decorrelate, but it does
  decorrelate:

  | days | En band | why |
  |---:|---:|---|
  | 1–10 (`--quick`) | 0.5 % | deterministic geostrophic adjustment; the console prints En to 4 digits (rounding alone is up to 0.2 %) |
  | 11–30 | 2 % | spin-up; the one-cell features (Celebes overflow, Gibraltar) begin to decorrelate |
  | 31–90 | 5 % | En peaks (day 65); the eddying part of the flow has decorrelated |
  | 91–365 | 10 % | slow spin-down of a decorrelated flow, still pinned by the initial state and the closures |

  MaxCFL, a pointwise maximum and far more sensitive to where one fast cell
  sits, gets twice En's band and must stay below 0.5 (the reference peaks
  at 0.16). On the reference toolchain the run is deterministic and matches
  every printed digit (all 365 days re-run at the reference commit: identical)
  — the bands exist for the others. gfortran 15.1 on the CPU (serial) also
  matches En and MaxCFL to every printed digit over the 10 quick days
  (a 5550 s run). The bands after day 10 are not yet measured across toolchains
  (a CPU year is days of wall time); they follow from how the flow evolves.
* **Budgets — round-off in this run, not equality with the reference.** The
  `Error` of a closed domain is accumulated round-off, which is
  toolchain-specific. Each series must stay under the envelope
  `|Error(d)| ≤ floor + rate · d` and must not accelerate (the mean daily
  increment of the second half of the run at most 4× that of the first
  half, plus 1e-12/day of jitter). `rate` is the per-step round-off that
  accumulates: 1e-13/day for mass and 1e-14/day for salt and heat — a few
  machine epsilons per step at 48 steps a day, 5× and 25× above the
  reference's own rates. `floor` = 1e-11 is the round-off of the global
  totals themselves, double sums over 5.8 M cells: the GPU's tree
  reductions are nearly exact, a serial CPU loop is not — gfortran 15.1
  measures a flat 1.3–1.5e-12 (mass) and a ±2e-13 jitter (salt, heat) from
  day 1 on. A real leak is 1e-10 or more in a day.

## 3. What `reproduce.sh` does, step by step

In order: fetch the inputs (idempotent, SHA-256 checked) → prepare the
model-grid inputs (skipped while `bathy_om1deg.nc` and `ic_woa13_jan.nc`
are newer than the files they are made from) → set up the run directory
(`INPUT` symlink, namelist copy with `t_end` from `--days`) → run → movie →
`check_against_reference.py`. Each step can be run by hand:

### Inputs (not in the repository)

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

### The reference run — the `rdb` executable

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

### The same run through the Python interface

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

### The movie from the executable's output

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

## 4. What the year shows

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
* The fastest water all year is the Celebes Sea overflow (§5) and the
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

## 5. Known limits of this configuration

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

## 6. Wind-forced: `global_1deg_wind.nml`

The same configuration with **one** change: a year of JRA55-do wind stress,
read from a file through `&ocean_dataovr_nml` and cycled. Grid, bathymetry,
vertical grid, WOA13 initial state, physics, time step and output are the
unforced case's (the diagnostics add the depth-integrated transport
`transport_x`/`transport_y`, which is output only). The whole run is driven
from the namelist by the `rdb` executable; no Python is involved.

**What is forced and what is not.** Wind stress only. There is **no heat
flux, no freshwater flux and no sea ice**: temperature and salinity change
only by transport and mixing, with no SST feedback and no restoring, so the
heat and salt budgets must close to round-off exactly as in the unforced run.
The stress is computed offline from the absolute wind — the ocean surface
velocity is ignored.

### 6.1 Inputs

```bash
export RDB_DATA_DIR=/somewhere/with/8GB
python3 tools/fetch_om1deg.py                 # OM_1deg grid + WOA13 (§3)
python3 tools/om1deg_prepare_inputs.py        # bathymetry + initial condition
python3 tools/fetch_om1deg.py --jra-wind      # -> $RDB_DATA_DIR/JRA55do/ (3.4 GB)
python3 tools/om1deg_prepare_wind.py          # -> $RDB_DATA_DIR/OM_1deg/wind_jra55do_1958_24h.nc
```

**The wind data.** JRA55-do **v1.4.0** (Tsujino et al. 2018, *Ocean
Modelling* 130, 79–139), `uas` and `vas`, the 10 m eastward/northward wind as
3-hourly point values on the TL319 Gaussian grid (640 × 320, ~0.56°), year
**1958**. `--jra-wind` streams GFDL's `reanalysis.tgz` (12.4 GB, the OM_1deg
forcing set, `ftp://ftp.gfdl.noaa.gov/perm/Alistair.Adcroft/MOM6-testing/`)
and keeps only the two wind members (1.87 + 1.89 GB, NetCDF-4, SHA-256
checked). GFDL's files are the `padded` ones: 2921 records, 1958-01-01 00:00
to 1959-01-01 00:00. That tarball holds 1958 only; another year has to come
from input4MIPs/ESGF (`source_id = MRI-JRA55-do-1-4-0`; the preprocessing
takes `--year`). **Licence** (the files' own `license` attribute): Creative
Commons Attribution-[NonCommercial-]ShareAlike 4.0, with the input4MIPs terms
of use (https://pcmdi.llnl.gov/CMIP6/TermsOfUse). Cite Tsujino et al. (2018),
and do not redistribute the data or the derived stress file under other terms.

**The preprocessing** (`tools/om1deg_prepare_wind.py`). For every 3-hourly
record, in this order:

1. **Bilinear interpolation** of `uas`, `vas` onto the model's C-grid faces.
   `taux` goes on the u faces (supergrid node `(2i-1, 2j)`, the west face of
   T-cell `i`, 361 × 320) and `tauy` on the v faces (node `(2i, 2j-1)`, the
   south face, 360 × 321). Longitude is periodic; beyond the outermost
   Gaussian row (±89.57°) the row value is used.
2. **Bulk stress** in geographic axes, `τ = ρ_air C_d(|U|) |U| (u_E, v_N)`,
   with `ρ_air = 1.22 kg m⁻³` and the Large & Yeager (2004, NCAR TN-460,
   eq. 6a) neutral 10 m drag coefficient
   `C_d = (2.7/U + 0.142 + 0.0764 U)·10⁻³`, with `U` floored at 0.5 m/s (the
   NCAR/FMS floor). `--cd 1.2e-3` selects a constant instead. Two
   simplifications: the JRA 10 m wind is used as if it were the neutral wind
   (no stability iteration), and the wind is absolute (the ocean velocity is
   not subtracted).
3. **Rotation onto the grid axes**, `τ_i = cos a τ_E + sin a τ_N`,
   `τ_j = −sin a τ_E + cos a τ_N` — the `ocean_metrics_t%angle_dx` convention,
   `a` measured counter-clockwise from east. The angle `a` is computed from
   the supergrid node positions: it is the direction of the chord between
   the face's two i-neighbours, projected on the local tangent plane in 3-D,
   so it stays accurate next to the pole. It is **not** read from the
   mosaic's `angle_dx`, which has two problems in OM_1deg's
   `ocean_hgrid.nc`:
   * **In the Arctic cap it has no `cos(lat)` factor.** It is the heading in
     degree space, `atan2(Δlat, Δlon)` — 10.4° against the true 32.7° at
     73° N, and more than 10° wrong at 69 000
     supergrid nodes. This is the error that
     MOM6's `GRID_ROTATION_ANGLE_BUGS = False` exists to avoid.
   * **Its edge nodes are one-sided.** Column 1 differs from its periodic
     twin, column 721, by up to 89°, and the fold row is not antisymmetric
     about the fold.

   South of the cap the two angles agree exactly (both are 0). With the
   computed angle, the written file is exactly periodic
   (`taux(1) ≡ taux(361)`) and exactly antisymmetric across the fold
   (`tauy(i, 321) ≡ −tauy(361−i, 321)`).

   The model's own `metrics%angle_dx` is read from the same mosaic variable.
   No kernel reads it yet, but the future in-model regridder must not trust
   it in the cap.
4. **Daily means with trapezoidal weights.** A record on a day boundary
   counts half to each day, so every day has exactly 8 weights and is
   centred on 12:00. The stress is formed at 3-hourly resolution *before*
   averaging: averaging the wind first would drop the synoptic variance from
   the quadratic law and underestimate the mean stress. `--avg-hours 6`
   would give 6-hourly records (4 × the file size).

The arithmetic (2920 records × 231 k faces, reading 3.4 GB of compressed
NetCDF-4) is done by a small Fortran helper, `tools/om1deg_wind_regrid.f90`.
The Python script compiles it with the system Fortran compiler and
`nf-config`, then runs it: **46 s** in all, 1 GB of memory. Pure Python
would take hours over the same 670 M face evaluations. Nothing needs to be
installed.

The output is `wind_jra55do_1958_24h.nc`: NetCDF classic (64-bit offset),
float32, **337 MB**. It holds `taux(xf, y, time)` and `tauy(x, yf, time)`,
and `time` is "days since 1958-01-01", 0.5 … 364.5. The layout is exactly
what the `tau_x`/`tau_y` tags read: Fortran order, and one face more than
there are cells in the staggered direction (`validation_examples/ocean/data_forcing/README.md`).
In the namelist:

```
&ocean_dataovr_nml
   enable       = .true.
   time_mode    = "cyclic"
   cycle_period = 31536000.0   ! 365 days
   t_offset     = 0.0          ! model t = 0 is 1958-01-01 00:00
   tau_x_file   = "INPUT/wind_jra55do_1958_24h.nc",  tau_x_var = "taux"
   tau_y_file   = "INPUT/wind_jra55do_1958_24h.nc",  tau_y_var = "tauy"
/
```

The reader interpolates linearly between the daily records and wraps across
the year boundary: day 0.0 is a 50/50 blend of Dec 31 and Jan 1.

The forcing it produces (annual mean over wet cells, zonal means in 5°
bands) is the familiar one:

* **Southern Hemisphere westerlies:** `τ_E` peaks at **+0.12 N m⁻²** at
  50–55° S.
* **Trade winds:** −0.04 to −0.05 N m⁻² at 10–25° N and S.
* **Northern Hemisphere westerlies:** +0.05 N m⁻² at 35–50° N.
* **Antarctic coastal easterlies:** −0.03 N m⁻².
* **Wind-stress curl:** negative over the northern subtropics (20–35° N,
  about −6·10⁻⁸ N m⁻³) and positive over the subpolar north (50–60° N). The
  southern subtropics are positive and the Southern Ocean south of 55° S is
  negative. These are the signs that drive the subtropical and subpolar
  gyres.
* **Extremes:** the largest daily mean on an ocean face is 5.5 N m⁻², on
  the Antarctic coast at 52° E (katabatic outflow). The largest 3-hourly
  value on any face, land included, is 13 N m⁻².

### 6.2 The run

```bash
mkdir run && cd run && ln -s "$RDB_DATA_DIR/OM_1deg" INPUT
cp /path/to/roundabout/validation_examples/ocean/global_1deg/global_1deg_wind.nml .
CUDA_VISIBLE_DEVICES=1 /path/to/build/rdb global_1deg_wind.nml > run.log
```

Measured on one V100, 2026-09-24/25: nvfortran 26.5,
`-DRDB_ENABLE_GPU=ON -DRDB_GPU_ARCH=cc70`, `mem:separate`. The year took
365 days and 17 520 steps in **3993 s of wall time (10.9 s per simulated
day)**, using 10.4 GB of device memory. The diagnostic file is 1.2 GB.

`reproduce.sh` and `check_against_reference.py` (§1–§2) cover the unforced
case only. This run's daily series is not a committed reference yet.

| day | En (m²/s²) | MaxCFL | Mass Error | Salt Error | Heat Error | max top-10 m speed (m/s), where | Drake (Sv) | Gulf Stream box | Kuroshio box | eq. Pacific u (m/s) |
|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|
| 1 | 2.691e-04 | 0.051 | -1.80e-14 | -6.4e-16 | -1.7e-16 | 0.81, Hudson Strait | 0.2 | 0.20 | 0.26 | -0.007 |
| 10 | 4.238e-04 | 0.112 | -1.80e-13 | -3.9e-15 | -3.0e-15 | 1.29, Irish shelf | 1.7 | 0.76 | 0.65 | -0.016 |
| 30 | 6.326e-04 | 0.168 | -5.41e-13 | -1.2e-14 | -8.8e-15 | 1.32, Taiwan Strait | -0.4 | 0.51 | 0.97 | -0.042 |
| 60 | 6.692e-04 | 0.161 | -1.08e-12 | -2.4e-14 | -1.8e-14 | 0.99, Tierra del Fuego coast | -1.5 | 0.99 | 0.55 | -0.084 |
| 90 | 6.635e-04 | 0.142 | -1.62e-12 | -3.6e-14 | -2.7e-14 | 1.17, Chukchi Sea, Alaska coast | -5.0 | 0.49 | 0.51 | -0.003 |
| 120 | 6.615e-04 | 0.129 | -2.16e-12 | -4.8e-14 | -3.5e-14 | 1.02, Hudson Bay east coast | -0.4 | 0.27 | 0.37 | -0.108 |
| 180 | 6.992e-04 | 0.092 | -3.25e-12 | -7.2e-14 | -5.2e-14 | 1.63, Sri Lanka coast (SW monsoon) | 5.6 | 0.36 | 0.51 | -0.290 |
| 240 | 6.907e-04 | 0.071 | -4.33e-12 | -9.8e-14 | -6.8e-14 | 1.37, Brazil coast 22° S | 4.8 | 0.62 | 0.27 | -0.327 |
| 300 | 6.346e-04 | 0.079 | -5.41e-12 | -1.2e-13 | -8.4e-14 | 2.04, Tierra del Fuego coast | 17.5 | 0.91 | 0.92 | -0.141 |
| 365 | 6.248e-04 | 0.050 | -6.58e-12 | -1.5e-13 | -1.0e-13 | 1.07, Kamchatka coast | 1.7 | 0.81 | 0.49 | -0.127 |

Columns:

* **Errors** are the relative closure residuals printed on the console.
* **The speed columns** are the fastest daily-mean water in the top 10 m
  (from the diagnostic file). The "box" columns take that maximum over
  80–60° W × 28–42° N (Gulf Stream) and 120–150° E × 24–40° N (Kuroshio).
* **Drake** is `Σ transport_x · dy` along 67.5° W, from Antarctica to South
  America.
* **eq. Pacific u** is the mean top-10 m zonal velocity over 160° E–100° W,
  2° S–2° N.

`python_prototypes/global_1deg/wind_analysis.py` produces all of these (one
row per day in `wind_daily.txt` there).

* **Stable.** The year ran without NaN, without a CFL truncation and
  without a positive-definite-limiter event. Nothing was logged as a
  warning. The `maxvel` clamp (6 m/s) was never reached: MaxCFL peaked at
  0.173 on day 38.
* **Energy is bounded.** `En` rises over the first two months, to 6.3e-04
  by day 30 (the unforced run: 5.4e-04). It reaches its maximum of
  **7.16e-04 m²/s² on day 183** (austral winter westerlies) and ends the
  year at 6.25e-04. It follows the seasonal cycle of the wind, not a trend.
* **Budgets close to round-off, linearly.** The residual grows at the same
  per-day rate as in the unforced run:

  | | per day | day 365 |
  |---|---:|---:|
  | mass | −1.81e-14 | −6.58e-12 |
  | salt | −4.1e-16 | −1.5e-13 |
  | heat | −2.8e-16 | −1.0e-13 |

  With no surface fluxes, the wind changes nothing here.

  The tracked boundary term `out` does change. It should be zero in this
  closed domain. It stays tiny: mass ≤ 82 kg of 1.4e21, salt ≤ 6.0e10 of
  4.8e22 (1e-12 relative), heat ≤ 1.2e12 J of 5.0e21 (2.5e-10 relative).
  But it is larger than the unforced run's: heat 1.2e12 J against
  ≤ 7e10 J over the unforced year (17×); at day 30, −8.5e9 J against
  ≤ 2.1e9 J. The printed heat total moves in its tenth digit, consistently.
  So some flux in this closed domain, larger when the flow is stronger, is
  booked as boundary flow; the tripolar fold seam is the first suspect. It
  is not a leak in the budget sense, since `Error` stays at round-off, but
  it is worth one look.
* **The fastest surface water** (top-10 m daily mean) is **2.70 m/s on
  day 338**, on the Labrador coast (56.5° W, 53.8° N). Other near-2 m/s
  maxima sit at single coastal cells:
  * the Antarctic coast at 85.5° E, 66.5° S (katabatic wind);
  * southern Brazil, 50.5° W, 30.5° S (2.44 m/s, day 362);
  * the Tierra del Fuego coast (2.04 m/s, day 300).

  These are wind-driven shelf jets in 9.5–50 m cells. The median over the
  year of the daily maximum is 1.31 m/s.

  The namelist output has no 3-D velocity, so the 3-D maximum is not
  reported. MaxCFL bounds it.

### 6.3 Physics sanity

* **The wind drives the western boundary currents.** Over the same boxes,
  the unforced run's surface maximum stays at 0.21–0.25 m/s (Gulf Stream)
  and 0.29–0.36 m/s (Kuroshio) all year. Forced, the maxima are:

  | | typical | peak |
  |---|---:|---:|
  | Gulf Stream | 0.5–1.0 m/s | 1.85 m/s on day 80 |
  | Kuroshio | 0.5–1.0 m/s | 1.33 m/s on day 361 |

  The Brazil–Malvinas confluence and the Agulhas show the same
  strengthening. That is the right order for a 1° model in its first year.
* **The subtropical gyres are there, but only as a barotropic spin-up.**
  Over the last 30 days, the SSH of the subtropical box minus the subpolar
  box is:

  | | model | observed (dynamic topography) |
  |---|---:|---:|
  | North Atlantic | +3.5 cm | ≈ 1 m |
  | North Pacific | +6.3 cm | ≈ 1 m |

  A year is enough for the barotropic Sverdrup response (Rossby-wave time
  scales of weeks) but not for the baroclinic one (decades).
* **Equatorial currents.** The top-10 m flow over the equatorial Pacific is
  westward all year: −0.1 to −0.33 m/s, strongest in austral winter–spring
  (days 180–240) with the trades. That is the South Equatorial Current. The
  Equatorial Undercurrent is below the 10 m diagnostic and not measured.
* **Drake Passage transport does NOT grow toward the observed ~130–170 Sv.**
  It fluctuates around zero all year: −5 Sv on day 90, +5.6 Sv on day 180,
  a maximum of 19.8 Sv on day 261, +1.7 Sv on day 365. The last 30 days
  average −0.3 Sv, and sections 4° either side agree to about 1 Sv. The
  unforced twin, re-run for 30 days with the same transport diagnostic,
  sits at −2.5 Sv, so the wind adds only a few Sv.

  The rest of the section is consistent with this:
  * The SSH step across the passage is 8–11 cm, where the observed step is
    about 1.2 m.
  * The top-10 m flow is eastward at about 0.1 m/s, but it is compensated
    by westward flow at depth.

  From rest with `η = 0`, the geostrophic adjustment of the WOA density puts
  the thermal-wind shear in the water column with ~zero depth-integrated
  transport. The wind's momentum input into the barotropic mode is then
  apparently taken out by topographic form stress within days.

  The expectation for this experiment was growth toward O(100) Sv within
  the year. This run does not show it. Nothing here says whether that is a
  roundabout defect or a consequence of this protocol (rest, η = 0, no GM,
  no buoyancy forcing). **It is the open question of this run.** The
  discriminating experiment is MOM6 OM_1deg on the same protocol: rest,
  WOA13 January, this stress file, no fluxes.

### 6.4 The movie

The movie is rendered from the diagnostic file of the namelist run. It needs
`nccopy` on PATH for the NetCDF-4 → NetCDF-3 flattening.

```bash
python3 validation_examples/ocean/global_1deg/global_movie.py \
    run/output/global_1deg_wind_rank_000000.nc movie_frames \
    --data "$RDB_DATA_DIR/OM_1deg" --interp bilinear --speed-max 1.0 \
    --stem global_1deg_wind --title "ROUNDABOUT  GLOBAL 1 DEG  JRA55-DO WIND 1958"
```

`--interp bilinear` blends the four model cells around each lon-lat pixel,
in the model's index space. The inversion is on the tangent plane, so the
Arctic cap needs no special case. Land corners get zero weight, so
coastlines do not bleed. The result is visibly smoother than the default
nearest-cell map, at about 1.3 s per frame (the year in about 8 min).

The 60 fps version is `ffmpeg minterpolate`; see
`python_prototypes/global_1deg/make_movie.py --fps60 mci`.

**One thing in the movie is not the wind.** The zonal stripes in the surface
speed, a few model rows apart in the Southern Ocean and elsewhere, are in the
unforced run too, row for row. For example, at 100° W, 60–70° S, rows 31, 34,
39 and 44 carry −0.01 m/s against +0.03 m/s in their neighbours on day 1 of
both runs. They come from the initial condition: `om1deg_prepare_inputs.py`
maps WOA13 onto the model grid by nearest neighbour, so neighbouring model
rows that fall on the same or on non-adjacent WOA rows give a staircase in
density, and the staircase has thermal-wind jets. Horizontal bilinear
interpolation in the IC preparation would remove them. That changes the
unforced reference, so it is left for its own change.
