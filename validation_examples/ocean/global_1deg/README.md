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
| Cost | 11.1 s per simulated day on one V100 (a year in ~68 min), 10.4 GB of device memory |

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
  | 31–90 | 5 % | En plateau (~5.0–5.2e-04, a secondary maximum on day 71); the eddying part of the flow has decorrelated |
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

**On several ranks** (MPI build, `-DRDB_ENABLE_MPI=ON`; several GPUs also need
`-DRDB_CUDA_AWARE_MPI=ON`): the same namelist runs unchanged — with
`&mpi_nml` unset the engine splits the tripolar grid north-south (`px = 1`,
`py = ranks`; the fold is single-rank in x) and every rank reads only its
band of the three input files.  Pin one GPU per rank BEFORE `MPI_Init`:

```bash
cat > pin.sh <<'SH'
#!/bin/bash
export CUDA_VISIBLE_DEVICES=$OMPI_COMM_WORLD_LOCAL_RANK
exec "$@"
SH
chmod +x pin.sh
mpirun -np 4 ./pin.sh /path/to/build/rdb global_1deg_unforced.nml > run.log
```

The state and the console are bit-identical to the 1-rank run (the console
uses the reproducing sums by default).  Measured, 2 days: 12.9 / 8.1 / 5.8 s
per simulated day on 1 / 2 / 4 V100s (`docs/CAPABILITIES_AND_LIMITATIONS.md`,
*MPI (domain decomposition)*, has the CPU numbers).

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

Measured on one V100 (nvfortran 26.5, `-gpu=cc70,mem:separate`), 2026-10-01,
on branch `fix/biharmonic-zfixed-closed-faces` (213e4e012) — the flow-aware
Smagorinsky biharmonic now treats `z_fixed` closed face-layers as free-slip
walls (it used to read the closed layers' zero velocities as open water,
i.e. a spurious drag at every staircase step): 365 days, 17 520 steps,
6771 s wall on a GPU shared with other jobs (the previous year took 4064 s
on an idle card; the change itself adds no measurable cost). The top-10 m speed
column is measured directly from this run's own diagnostic file (a small
stdlib NetCDF scan, the same approach as
`python_prototypes/global_1deg/wind_analysis.py`).

| day | En (m²/s²) | MaxCFL | Mass Error | Salt Error | Heat Error | max top-10 m speed (m/s) and where |
|---:|---:|---:|---:|---:|---:|---|
| 1 | 5.779e-04 | 0.049 | -1.80e-14 | -5.2e-16 | -4.2e-16 | 0.58, Cape Hatteras shelf |
| 10 | 5.307e-04 | 0.241 | -1.80e-13 | -3.7e-15 | -2.7e-15 | 0.86, Taiwan Strait |
| 30 | 5.259e-04 | 0.235 | -5.41e-13 | -1.2e-14 | -8.6e-15 | 0.83, North Carolina shelf |
| 60 | 5.330e-04 | 0.229 | -1.08e-12 | -2.3e-14 | -1.9e-14 | 0.72, Taiwan Strait |
| 90 | 5.352e-04 | 0.230 | -1.62e-12 | -3.6e-14 | -2.8e-14 | 0.67, Taiwan Strait |
| 180 | 5.160e-04 | 0.177 | -3.25e-12 | -7.1e-14 | -5.7e-14 | 0.70, North Carolina shelf |
| 240 | 4.841e-04 | 0.128 | -4.33e-12 | -9.5e-14 | -7.5e-14 | 0.71, Bering Strait |
| 300 | 4.579e-04 | 0.107 | -5.41e-12 | -1.2e-13 | -9.3e-14 | 0.72, Bering Strait |
| 365 | 4.329e-04 | 0.091 | -6.58e-12 | -1.4e-13 | -1.1e-13 | 0.73, Bering Strait |

(Speeds are the daily means of the diagnostic file's one top-10 m level;
the file carries no deeper velocity.)

* **No NaN, no CFL truncation, no positive-definite-limiter event, no
  NaN-catch**: each of these is logged when it fires, and the console has
  none. The 6 m/s `maxvel` clamp keeps no counter, and the 3-D maximum
  speed was not re-measured for this year (the diagnostic file holds only
  the top 10 m; an earlier reference year's 3-D maximum was 4.33 m/s, in
  the Celebes trench, day 38 — see §5). MaxCFL stays at or below 0.24728
  all year (day 9).
* **Why the curve differs from the previous reference.** The previous
  reference (`fix/free-slip-visc-sponge-seed`, 2026-09-30) ran the same
  physics with the flow-aware Smagorinsky biharmonic reading `z_fixed` closed
  face-layers as open water: their velocities are held at zero, so every
  staircase step acted as a no-slip wall for the biharmonic — a spurious
  drag. With closed faces treated as free-slip walls, the circulation keeps
  more energy: En is higher by up to 10.7 % and MaxCFL by up to 37 %,
  relative, over the year (day 365: 4.329e-04 vs 3.946e-04, MaxCFL 0.091 vs
  0.066). The early spin-up shape (En peaking on day 3) is unchanged; the
  difference grows through the year, mostly in the slow spin-down.
* **Energy is bounded.** En peaks at 6.042e-04 m²/s² on day 3, settles to a
  plateau of 5.2–5.4e-04 through day ~100 (a secondary maximum of 5.388e-04
  on day 99), and then decays slowly (4.329e-04 at day 365): with no forcing
  the closures spin it down.
* **Budgets close to round-off, linearly.** The relative residuals grow at
  a constant rate — mass −1.80e-14 per day, salt −4.0e-16, heat −3.1e-16 —
  i.e. round-off accumulating, not a leak; the tracked boundary fluxes
  (`out`) stay at round-off size in this closed domain (mass at most 21 kg
  of 1.4e21, heat at most 3.8e5 J of 5.0e21).
* **The fastest surface water is in one-cell shallow straits and shelves.**
  (Details below measured on the previous reference year; this year's table
  rows put the daily maximum at the same sites within 0.03 m/s.)
  The top-10 m daily mean peaks at 0.94 m/s on day 5 in the Taiwan Strait (a
  17 m cell), just above the 0.92 m/s reached on day 2 on the Yucatán shelf
  (9.5 m); after that the domain maximum, 0.6–0.9 m/s, sits mostly in the
  Taiwan Strait (days 4–111, with brief exceptions), on the North Carolina
  shelf south of Cape Hatteras (on and off, days 16–184) and in the Bering
  Strait (from day 112 on) — cells 15.5–42 m deep that carry the barotropic
  flow. Elsewhere the equatorial current system and the Antarctic
  Circumpolar Current's fronts carry the surface flow (a Drake Passage box,
  72–63° W × 70–52° S, peaks at about 0.51 m/s, day 327), and the movie
  shows them spinning up and slowly decaying. The Drake Passage transport
  is not measured: the diagnostic file carries no depth-integrated transport.
* **The Florida Straits jet is gone.** (Previous reference year.) With uniform 130 m layers (the
  scoping probe) a 9.5 m cell next to 587 m cells in the Straits carried
  6 m/s by day 2 — with bed-only or HBBL drag alike. With the stretched
  profile the same cells are 2–5 m layers like their neighbours: the
  Straits stayed below 0.7 m/s over the first 12 days of the probe, and
  over this year their top-10 m speed peaks at 0.53 m/s (day 46) and never
  carries the surface maximum.
* **The Python-driven run is bit-identical to the executable's** (checked on
  an earlier reference year; not re-run for this one): every value of
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

Measured on one V100, 2026-09-30: nvfortran 26.5,
`-DRDB_ENABLE_GPU=ON -DRDB_GPU_ARCH=cc70`, `mem:separate`, on branch
`fix/free-slip-visc-sponge-seed` — the same two fixes as the unforced year
(§4): the velocity-form Laplacian/biharmonic is now free-slip at land, and
`z_fixed` layers are now seeded on their target when zinit sets T/S. The
year took 365 days and 17 520 steps in
**4067 s of wall time (11.1 s per simulated day)**, using 10.4 GB of device
memory. The diagnostic file is 1.2 GB.

`reproduce.sh` and `check_against_reference.py` (§1–§2) cover the unforced
case only. This run's daily series is not a committed reference yet.

| day | En (m²/s²) | MaxCFL | Mass Error | Salt Error | Heat Error | max top-10 m speed (m/s), where | Drake (Sv) | Drake SSH step (m) | Gulf Stream box | Kuroshio box | eq. Pacific u (m/s) |
|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|
| 1 | 5.983e-04 | 0.068 | -1.80e-14 | -5.2e-16 | -4.2e-16 | 0.98, Hudson Strait | -175.7 | 0.37 | 0.57 | 0.42 | -0.029 |
| 10 | 5.785e-04 | 0.238 | -1.80e-13 | -3.7e-15 | -2.7e-15 | 1.26, North Sea (1.5° W, 55.5° N) | 156.9 | 1.67 | 0.43 | 0.49 | -0.079 |
| 30 | 5.655e-04 | 0.241 | -5.41e-13 | -1.1e-14 | -8.4e-15 | 1.12, Gulf of Maine | 147.4 | 1.51 | 0.61 | 0.45 | -0.091 |
| 60 | 5.527e-04 | 0.225 | -1.08e-12 | -2.3e-14 | -1.8e-14 | 1.38, Gulf Stream (80.5° W, 31.3° N) | 151.1 | 1.44 | 1.38 | 0.47 | -0.084 |
| 90 | 5.590e-04 | 0.219 | -1.62e-12 | -3.5e-14 | -2.7e-14 | 1.20, Chukchi Sea, Alaska coast | 148.8 | 1.44 | 0.45 | 0.55 | -0.018 |
| 120 | 5.597e-04 | 0.201 | -2.16e-12 | -4.7e-14 | -3.6e-14 | 1.09, Hudson Bay east coast | 152.0 | 1.47 | 0.92 | 0.63 | -0.114 |
| 180 | 6.067e-04 | 0.152 | -3.25e-12 | -7.1e-14 | -5.4e-14 | 1.75, Sri Lanka coast (SW monsoon) | 158.8 | 1.51 | 0.66 | 0.80 | -0.269 |
| 240 | 6.081e-04 | 0.119 | -4.33e-12 | -9.6e-14 | -7.2e-14 | 1.50, Brazil coast 22° S | 147.7 | 1.45 | 0.56 | 0.51 | -0.329 |
| 300 | 5.563e-04 | 0.098 | -5.41e-12 | -1.2e-13 | -8.9e-14 | 2.06, Tierra del Fuego coast | 159.4 | 1.51 | 0.77 | 0.74 | -0.133 |
| 365 | 5.418e-04 | 0.085 | -6.58e-12 | -1.5e-13 | -1.1e-13 | 1.12, Kamchatka coast | 152.9 | 1.51 | 0.77 | 0.45 | -0.125 |

Columns:

* **Errors** are the relative closure residuals printed on the console.
* **The speed columns** are the fastest daily-mean water in the top 10 m
  (from the diagnostic file). The "box" columns take that maximum over
  80–60° W × 28–42° N (Gulf Stream) and 120–150° E × 24–40° N (Kuroshio).
* **Drake** is `Σ transport_x · dy` along 67.5° W, from Antarctica to South
  America. **Drake SSH step** is SSH at the northernmost wet cell of that
  section minus SSH at the southernmost.
* **eq. Pacific u** is the mean top-10 m zonal velocity over 160° E–100° W,
  2° S–2° N.

`python_prototypes/global_1deg/wind_analysis.py` produces all of these (one
row per day in `wind_gfix_daily.txt` there).

* **Stable.** The year ran without NaN, without a CFL truncation and
  without a positive-definite-limiter event. Nothing was logged as a
  warning. MaxCFL peaked at 0.2437 on day 28, while the barotropic mode
  adjusts; the `maxvel` clamp (6 m/s) is far above anything the diagnostics
  show.
* **Energy is bounded.** `En` peaks at **6.343e-04 m²/s² on day 3**, during
  the barotropic adjustment to the WOA density field, and then stays at
  5.3–6.2e-04 for the rest of the year, higher in austral winter (6.1e-04 on
  days 180–240) with the westerlies. It ends at 5.42e-04. The unforced year
  over the same days: 5.98e-04 on day 3, 3.95e-04 on day 365.
* **Budgets close to round-off, linearly.** The residual grows at the same
  per-day rate as in the unforced run:

  | | per day | day 365 |
  |---|---:|---:|
  | mass | −1.80e-14 | −6.58e-12 |
  | salt | −4.1e-16 | −1.5e-13 |
  | heat | −3.0e-16 | −1.1e-13 |

  With no surface fluxes, the wind changes nothing here.

  The tracked boundary term `out` should be zero in this closed domain. It
  stays tiny, at most mass 137 kg of 1.4e21, salt 9.9e10 of 4.8e22 (2e-12
  relative), heat 3.5e11 J of 5.0e21 (7e-11 relative). It is larger than
  the unforced year's (mass at most 11.4 kg, salt at most 5.1e6, heat at
  most 2.5e5 J, §4): some flux in this closed domain, larger when the flow
  is stronger, is booked as boundary flow. **Measured, not the tripolar
  fold**: a single-rank 2-day
  re-run of both namelists (`gfortran`/nvfortran cc70 agree; `state`/`En`
  0.000% against `reference_daily.csv`) gives day-2 `out` of mass 1.45 kg /
  salt -62.6 / heat -3.84 (unforced) vs mass 0.067 kg / salt 602 / heat
  -11.1 (wind-forced) — both already at round-off-telescoping SCALE
  (relative to the ~1e21-1e22 totals, 1e-20 to 1e-23) after only 2 of 365
  days, growing with total solver work (more arithmetic ⇒ more accumulated
  round-off) rather than with anything fold-specific: the sign flips
  between runs and between quantities, salt grew but mass shrank going
  from unforced to wind-forced, which a systematic per-step fold leak
  would not do. Earlier wording blamed the fold as "the first suspect";
  that was speculation made before this measurement, not evidence for it.
  It is not a leak in the budget sense, since `Error` stays at round-off.
* **The fastest surface water** (top-10 m daily mean) is **3.02 m/s on day
  338**, at a single cell in the Labrador Sea (56.5° W, 53.8° N) — essentially
  tied with the Antarctic-coast katabatic maximum, 3.01 m/s on day 71 at
  86.5° E, 66.5° S; each lasts a day. Other near-2 m/s maxima sit at single
  coastal cells: the Brazil–Malvinas box reaches 2.34 m/s (day 265) and the
  Tierra del Fuego coast 2.38 m/s (day 288). These are wind-driven shelf
  jets in 9.5–50 m cells. The median over the year of the daily maximum is
  1.40 m/s.

  The namelist output has no 3-D velocity, so the 3-D maximum is not
  reported. MaxCFL bounds it.

### 6.3 Physics sanity

The MOM6 twin of this protocol — MOM6 `dev/gfdl` `d74a11f9c`, same grid,
bathymetry, vertical grid, WOA13 state and stress file, run 90 days — is
recorded in `python_prototypes/mom6_baselines/global_1deg_wind/`. It is the
yardstick below.

* **Drake Passage carries ~150 Sv, as in MOM6.**

  | | roundabout | MOM6 |
  |---|---:|---:|
  | day 1 (the adjustment surge) | −175.7 Sv | −188.0 Sv |
  | day 10 | 156.9 Sv | 157.9 Sv |
  | mean, days 10–90 | 151.9 Sv (140.3–163.1) | 161 Sv (143–171) |
  | mean, days 91–365 | 150.8 Sv (140.1–164.7) | not run |
  | day 365 | 152.9 Sv | not run |
  | SSH step, days 10–365 | 1.33–1.67 m (mean 1.48 m) | 1.56–1.59 m (days 10–90) |

  The observed transport is about 130–170 Sv; the observed step about
  1.2–1.5 m. Sections 4° either side (i = 229, 237) agree to about 1.5 Sv.
  Like MOM6, the transport comes from the WOA13 density field, not from the
  wind: starting from rest with `η = 0`, the depth-integrated baroclinic
  pressure gradient over the topography forces the barotropic mode, which
  answers within days by external gravity-wave adjustment — the −175 Sv
  surge on day 1, then about 150 Sv by day 10 with a 1.5 m SSH step across
  the passage. The wind's share is a few Sv: MOM6's no-wind control carries
  6–9 Sv less on days 10–15.

  An earlier version of this run gave **~0 Sv** all year and an SSH step of
  a few cm. That was a roundabout defect, and this comparison is what found
  it: the split-explicit scheme discarded the depth mean of the baroclinic
  PGF, so the barotropic mode never felt JEBAR. MOM6 forces the barotropic
  mode with that depth mean minus only its free-surface part
  (`BT_force`/`eta_PF`); roundabout now does the same
  (`&ocean_bt_nml bc_pgf_forcing`, default on). That fix alone gave
  73–83 Sv. The FV_MOM6 PGF then evaluated its constant-by-layer density at
  the surface reference pressure rather than in situ; with the in-situ
  density (`&ocean_pgf_nml insitu_density`, default on) the transport
  reaches the ~150 Sv here.
* **The subtropical gyres match MOM6.** The SSH of the subtropical box
  minus the subpolar box (N. Atlantic 70–40° W × 25–35° N against
  50–30° W × 50–60° N; N. Pacific 150° E–170° W × 20–32° N against
  160° E–170° W × 45–55° N):

  | | roundabout, days 61–90 | MOM6, days 61–90 | roundabout, days 336–365 | observed (dynamic topography) |
  |---|---:|---:|---:|---:|
  | North Atlantic | +1.031 m | +1.031 m | +0.956 m | ≈ 1 m |
  | North Pacific | +0.961 m | +0.967 m | +0.979 m | ≈ 1 m |

  As with Drake, most of this is the initial density field, which the
  barotropic mode now adjusts to; the wind's Sverdrup response is a small
  part in one year.
* **The western boundary currents are too fast.** Maximum top-10 m speed
  in the boxes:

  | | roundabout, median / peak | roundabout max, days 1–90 | MOM6 max, days 1–90 | unforced roundabout |
  |---|---:|---:|---:|---:|
  | Gulf Stream | 0.70 / 1.95 m/s (day 80) | 1.95 | 0.71 | 0.38–0.81 |
  | Kuroshio | 0.54 / 1.24 m/s (days 205, 361) | 1.09 | 0.74 | 0.37–0.54 |

  So the peaks are still 1.5–2.7× MOM6's — 2.74× (Gulf Stream) and 1.47×
  (Kuroshio), against 2.72× and 1.47× on the previous reference. **Neither
  the barotropic fix nor the free-slip coastal viscosity changed them
  materially**: the free-slip Laplacian/biharmonic was the standing
  candidate explanation (a no-slip wall drains momentum from the boundary
  layer that free-slip lets accumulate), and it moved the Gulf Stream peak
  by +0.8 % and the Kuroshio peak by −0.4 %, both within the day-to-day
  noise of which single cell holds the box maximum. The first suspect
  remains the initial state: `om1deg_prepare_inputs.py` maps WOA13 onto the
  model grid by nearest neighbour, and the resulting density staircase
  carries thermal-wind jets (§6.4), where MOM6 regrids WOA horizontally. A
  second is lateral viscosity near the boundary, where OM_1deg adds a 2-D
  background viscosity file that this configuration does not read.
  **This is still the open question of this run.**
* **Equatorial currents.** The top-10 m flow over the equatorial Pacific is
  westward all year, −0.008 to −0.34 m/s, strongest in austral
  winter–spring (days 180–240) with the trades: the South Equatorial
  Current. On day 30 roundabout has −0.09 m/s where MOM6 has +0.12 m/s;
  worth a look together with the boundary currents. The Equatorial
  Undercurrent is below the 10 m diagnostic and not measured.

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
