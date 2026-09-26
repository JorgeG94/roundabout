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
| Cost | 12.7 s per simulated day on one V100 (a year in 77 min), 10.4 GB of device memory |

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
  | 1–10 (`--quick`) | 0.5 % | deterministic adjustment: the barotropic mode answers the baroclinic pressure field within days (En peaks on day 3), then settles; the console prints En to 4 digits (rounding alone is up to 0.2 %) |
  | 11–30 | 2 % | spin-up; the one-cell straits and shelves that carry the fastest water begin to decorrelate |
  | 31–90 | 5 % | En plateau (~5.1e-04, a secondary maximum on day 63); the eddying part of the flow has decorrelated |
  | 91–365 | 10 % | slow spin-down of a decorrelated flow, still pinned by the initial state and the closures |

  MaxCFL, a pointwise maximum and far more sensitive to where one fast cell
  sits, gets twice En's band and must stay below 0.5 (the reference peaks
  at 0.243 on day 28). On the reference toolchain the run is deterministic
  and matches every printed digit (the previous reference year, before the
  barotropic-split fix, re-ran all 365 days identically) — the bands exist
  for the others. On that previous reference, gfortran 15.1 on the CPU
  (serial) also matched En and MaxCFL to every printed digit over the 10
  quick days (a 5550 s run); not yet repeated on this one. The bands after
  day 10 are not yet measured across toolchains
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
previous reference year's run, every record of SSH, u, v, T and S).

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

Measured on one V100 (nvfortran 26.5, `-gpu=cc70,mem:separate`), 2026-09-25,
code `f995eeef2` (the MOM6 barotropic split, `&ocean_bt_nml bc_pgf_forcing`,
plus the FV_MOM6 in-situ density): 365 days, 17 520 steps, 4637 s wall
(12.7 s per simulated day), 10.4 GB of device memory.

| day | En (m²/s²) | MaxCFL | Mass Error | Salt Error | Heat Error | max top-10 m speed (m/s) and where |
|---:|---:|---:|---:|---:|---:|---|
| 1 | 5.764e-04 | 0.046 | -1.80e-14 | -5.6e-16 | -2.4e-16 | 0.59, Cape Hatteras shelf, 30 m cell |
| 10 | 5.166e-04 | 0.238 | -1.80e-13 | -3.7e-15 | -3.3e-15 | 0.85, Taiwan Strait, 17 m cell |
| 30 | 5.067e-04 | 0.239 | -5.41e-13 | -1.2e-14 | -9.8e-15 | 0.80, North Carolina shelf, 16 m cell |
| 60 | 5.097e-04 | 0.224 | -1.08e-12 | -2.4e-14 | -2.1e-14 | 0.70, Taiwan Strait |
| 90 | 5.100e-04 | 0.220 | -1.62e-12 | -3.7e-14 | -3.2e-14 | 0.65, Taiwan Strait |
| 180 | 4.614e-04 | 0.146 | -3.25e-12 | -7.3e-14 | -6.4e-14 | 0.67, North Carolina shelf |
| 240 | 4.390e-04 | 0.112 | -4.33e-12 | -9.7e-14 | -8.4e-14 | 0.67, Bering Strait, 42 m cell |
| 300 | 4.172e-04 | 0.087 | -5.41e-12 | -1.2e-13 | -1.0e-13 | 0.69, Bering Strait |
| 365 | 3.941e-04 | 0.066 | -6.58e-12 | -1.5e-13 | -1.2e-13 | 0.71, Bering Strait |

(Speeds are the daily means of the diagnostic file's one top-10 m level;
the file carries no deeper velocity.)

* **No NaN, no CFL truncation, no positive-definite-limiter event, no
  NaN-catch**: each of these is logged when it fires, and the console has
  none. The 6 m/s `maxvel` clamp keeps no counter, and the 3-D maximum
  speed was not re-measured for this year (the diagnostic file holds only
  the top 10 m; the previous reference year's 3-D maximum was 4.33 m/s, in
  the Celebes trench, day 38). MaxCFL stays at or below 0.243 all year.
* **Why the curve differs from the previous reference.** That year was
  measured before the barotropic mode was driven by the depth mean of the
  baroclinic pressure gradient (JEBAR, MOM6's barotropic split). Now the
  barotropic mode answers the WOA pressure field at once: En reaches
  5.76e-04 on day 1 (was 2.53e-04) and peaks on day 3 instead of rising
  over two months; the faster, stronger barotropic flow lifts the MaxCFL
  peak from 0.16 to 0.243 (day 28).
* **Energy is bounded.** En peaks at 5.924e-04 m²/s² on day 3, settles to a
  plateau of 5.0–5.1e-04 through day ~100 (a secondary maximum of 5.13e-04
  on day 63), and then decays slowly (3.94e-04 at day 365): with no forcing
  the closures spin it down.
* **Budgets close to round-off, linearly.** The relative residuals grow at
  a constant rate — mass −1.80e-14 per day, salt −4.0e-16, heat −3.4e-16 —
  i.e. round-off accumulating, not a leak; the tracked boundary fluxes
  (`out`) stay at round-off size in this closed domain (mass ≤ 1e2 kg of
  1.4e21, heat ≤ 6e10 J of 5.0e21).
* **The fastest surface water is in one-cell shallow straits and shelves.**
  The top-10 m daily mean peaks at 0.92 m/s on day 2 on the Yucatán shelf
  (a 10 m cell); after that the domain maximum, 0.6–0.9 m/s, sits in the
  Taiwan Strait (days 4–115), on the North Carolina shelf south of Cape
  Hatteras (on and off, days 17–197) and in the Bering Strait (from day
  116 on) — cells 16–42 m deep that carry the barotropic flow. The previous
  reference's surface never exceeded 0.50 m/s. Elsewhere the equatorial
  current system and the Antarctic Circumpolar Current's fronts carry the
  surface flow (the Drake Passage box peaks at 0.49 m/s), and the movie
  shows them spinning up and slowly decaying. The Drake Passage transport
  is not measured: the diagnostic file carries no depth-integrated transport.
* **The Florida Straits jet is gone.** With uniform 130 m layers (the
  scoping probe) a 9.5 m cell next to 587 m cells in the Straits carried
  6 m/s by day 2 — with bed-only or HBBL drag alike. With the stretched
  profile the same cells are 2–5 m layers like their neighbours: the
  Straits stayed below 0.7 m/s over the first 12 days of the probe, and
  over this year their top-10 m speed peaks at 0.52 m/s (day 43) and never
  carries the surface maximum.
* **The Python-driven run is bit-identical to the executable's** (checked on
  the previous reference year; not re-run for this one): every value of
  every record of SSH, T, S, u and v in the two diagnostic files (43.5 M
  values per field, ghosts included) is equal, and so are En and the mass
  drift each day.

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
  jet. In the previous reference year (before the barotropic-split fix; the
  3-D field was not re-measured on this one) it reached 4.3 m/s at 4800 m by
  day 38, saturated by the bottom drag, then decayed to 0.9 m/s by day 365.
  It is the model geometry disagreeing with the observed state, not an
  instability, and it was the fastest flow in that run. With bed-only drag (see the next point)
  the same jet reaches the 6 m/s clamp by day 19.
* **`z_fixed` + bed-only bottom drag gives no bottom drag.** The bed-only
  mode (`hbbl = 0`) drags layer `k = 1`, which is an inert filler in every
  column shallower than the deepest level; configure warns. This namelist
  uses `hbbl = 10` (MOM6's value), which reaches the live bottom layer.
* **No channel list, no `CHANNEL_DRAG` file, no GM/MEKE, KPP instead of
  EPBL, no neutral diffusion, no sea ice.** None of them is needed for a
  stable year; all of them matter for a realistic one.
