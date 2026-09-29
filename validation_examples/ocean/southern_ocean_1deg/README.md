# Southern Ocean 1 degree, wind-forced — a single-GPU cut of OM_1deg

`southern_ocean_1deg_wind.nml` runs a **regional cut** of the MOM6
**OM_1deg** tripolar geometry — everything from the Antarctic coast up to
30.02 S, periodic east-west, no Arctic cap — from a **WOA13 January** ocean
at rest, forced by one (cycled) year of **JRA55-do 1958** wind stress, for
**two simulated years**. It is the same physics, vertical grid, time step
and closures as
[`../global_1deg/global_1deg_wind.nml`](../global_1deg/global_1deg_wind.nml),
with **no physics change**: the lateral viscosity is the global run's
(Smagorinsky + `nu_h = 2000` m²/s floor, `ah_max` 10000). The only things
that change are the domain (Southern Ocean instead of global) and the north
edge (a sponge instead of the tripolar fold). An earlier version of this
namelist raised the viscosity to 150000 m²/s to survive a day-573 blow-up;
that blow-up was a model defect in the sponge edge, now fixed — see
"Stability" (§5). The
point is a single-GPU, few-hour run of the Antarctic Circumpolar Current
that a laptop-class movie script can turn into a "watch it spin" video —
see `southern_movie.py` below.

| | |
|---|---|
| Grid | OM_1deg `ocean_hgrid.nc` supergrid, cut to tracer rows 0-96 (360 x 97), `west/east = "periodic"`, `south = "wall"` (Antarctica), `north = "sponge"` |
| Domain | 77.88 S (south wall, the same row the global grid uses) to 30.02 S (north sponge edge) |
| Bathymetry | Same y-cut of OM_1deg `topog.nc`, same MOM6 limits (`MINIMUM_DEPTH` 9.5 m, `MAXIMUM_DEPTH` 6500 m, `MASKING_DEPTH` 0) |
| Vertical | 50 `z_fixed` levels, tanh-stretched 2 m (surface) to ~255 m (bed), partial-step faces closed — identical to the global run |
| Initial state | WOA13 decav January T/S, same y-cut, at rest, eta = 0 |
| Forcing | JRA55-do 1958 daily-mean wind stress, same y-cut, cycled annually |
| Time step | `dt = 1800 s`, `pred_corr`, barotropic `n_inner` from the per-wet-cell CFL |
| North BC | `&ocean_bc_nml north = "sponge"`: OBC_SPONGE (wall at the outer face + a 5-cell/~5 degree interior relaxation band), `&ocean_sponge_nml target_source = "ic"` — relax T/S/u/v toward the WOA13 snapshot |
| Output | daily: SSH, top-10 m u/v/T (T doubling as SST), surface `vorticity_z`, depth-integrated `transport_x`/`transport_y`; single precision, 2-D only |
| Length | 730 days (2 years) |

## 1. Why this domain, and why a sponge (not an open/Flather edge)

**South.** The OM_1deg grid is a plain lon-lat grid south of the Arctic
tripolar cap (checked directly against `ocean_hgrid.nc`: at every column,
row `j=0`'s latitude is 77.8796 S and it is exactly uniform across the row
— see "How the inputs were cut" below), and row `j=0` is already the wall
the global runs (`global_1deg_unforced.nml`, `global_1deg_wind.nml`) use
for `south = "wall"`. Cutting at `j=0` reuses that wall unchanged; no new
land mask.

**North, 30.02 S.** The supergrid has a node row (a cell EDGE, not a
straddled cell centre) at exactly 30.0234 S for every column: supergrid
node row 194 (0-indexed), the north face of T-row 96. Keeping tracer rows
`j=0..96` (97 rows) puts the domain's north edge exactly there — a clean
cut, chosen for landing almost exactly on the requested "about 30 S", not
for any dynamical reason.

**Sponge, not an open (Flather) edge.** `&ocean_bc_nml` OPEN/Flather needs
a boundary DATA series (`data_eta`, and per-tracer inflow values) to
radiate against; this cut has none — there is no parent-model boundary
archive for this line, only the same static WOA13 climatology the interior
is seeded from. `OBC_SPONGE` is the honest use of what's actually
available: the outer face is a wall (no barotropic through-flow leaks out
uncontrolled), and a `~5` cell (`sponge_width = 5`, roughly 5 degrees at
this resolution) interior band relaxes momentum toward `u_ref`/`v_ref`
(zero — the WOA13 state starts at rest) and every tracer toward the WOA13
snapshot, `target_source = "ic"`, on a 10-day timescale
(`sponge_strength = 1.1574074e-4 = 1/(10 day)`, the same magnitude the
`isomip_plus` examples use). This *is* the model's existing OBC/sponge
machinery (`docs/CLOSURE_MATRIX.md`, `docs/REFERENCE.md` §4) — no new
physics, per the brief. SSH is not separately restored: the real sponge's
interior-thickness relaxation (`&ocean_sponge_nml relax_h`) is not
implemented in this tree (v1), so the wall face is what pins the
barotropic mode at the edge; T/S/u/v restoring is what's live.

**No tripolar fold.** The cut stops 30 degrees of latitude south of the
Arctic cap, so `north = "tripolar_fold"` was never in play; the fold
exchange code is untouched.

## 2. The supergrid reader accepts a non-global subset

`metrics_fill_from_supergrid` (`src/core/ocean/state/rdb_ocean_metrics.F90`)
reads `ni = grid%nx_global`, `nj = grid%ny_global` from `&grid_nml nx/ny`
and only checks that the supergrid file's `nxp/nyp/nx/ny` dimensions equal
`2*ni+1`/`2*nj+1`/`2*ni`/`2*nj` — it never assumes the file is a *global*
grid. The topology cross-checks are both satisfied by this cut without any
source change:

* **Tripolar-fold detection** (`supergrid_top_row_folds`): compares every
  node of the file's TOP row against its mirror `nxp+1-m` as unit vectors
  on the sphere (both longitude and latitude). Our cut's top row (30.02 S)
  is an ordinary constant-latitude row with genuinely different longitudes
  at `m` and at its mirror, so the test correctly returns `.false.` —
  matching `&ocean_bc_nml north = "sponge"` (not `tripolar_fold`).
* **Periodic east-west**: checked as `y(1,:) == y(nxp,:)`, which holds
  because the cut keeps the full longitude circle (only `y` is cut).

So no source change was needed for this configuration; `nx/src` gate
scope in "Gates" below is empty.

## 3. Inputs — how they were cut

All four OM_1deg files this configuration needs share the same convention
south of the cap: a plain `(y, x[, z|time])`-style grid, `y` index 0 at the
south wall. A single contiguous slice along `y` (or the supergrid's
`ny`/`nyp`) is everything all four inputs need — there is no
reprojection, no re-triangulation, nothing i-dependent. `tools/
om1deg_subset.f90` is a small generic NetCDF dimension-range subsetter
(same `nf90_*` calling style as `tools/om1deg_wind_regrid.f90`): given
`dimname:start0:count` specs on the command line, it copies every
dimension, variable and attribute of the input file unchanged except for
the named dimensions, which are truncated. One tool, four invocations:

```bash
source /home/jorge/nci/cdx/roundabout/nvhpc_env.sh   # or gcc_env.sh — any NetCDF-Fortran toolchain
gfortran -O2 $(nf-config --fflags) -o om1deg_subset tools/om1deg_subset.f90 $(nf-config --flibs)

mkdir -p /home/jorge/nci/cdx/data/OM_1deg_southern
D=/home/jorge/nci/cdx/data/OM_1deg
O=/home/jorge/nci/cdx/data/OM_1deg_southern

./om1deg_subset $D/ocean_hgrid.nc            $O/ocean_hgrid.nc            nyp:0:195 ny:0:194
./om1deg_subset $D/bathy_om1deg.nc           $O/bathy_om1deg.nc           y:0:97
./om1deg_subset $D/ic_woa13_jan.nc           $O/ic_woa13_jan.nc           y:0:97
./om1deg_subset $D/wind_jra55do_1958_24h.nc  $O/wind_jra55do_1958_24h.nc  y:0:97 yf:0:98
```

(`ocean_hgrid.nc`'s supergrid dims double-plus-one: 97 T-rows -> `nyp = 2*97+1 = 195`
node rows, `ny = 2*97 = 194` cell-row segments. The bathymetry, IC and wind
files are already on the model T/C-grid, so they take the T-row count `97`
— and the wind file's meridional FACE dimension `yf` takes `97+1 = 98`, one
more than the T-row count, same as the global file's `yf = ny+1`
convention.) Verified: the cut `depth` field's min/max are 0/6500 m,
matching the global bathymetry's own limits (land is still 0, the deep
cap is still 6500); every cut file's global attributes are the source
file's plus one provenance note. Output sizes: `ocean_hgrid.nc` 6.7 MB,
`bathy_om1deg.nc` 0.28 MB, `ic_woa13_jan.nc` 27 MB, `wind_jra55do_1958_24h.nc`
98 MB (versus 337 MB uncut) — about 3.4 GB less to move around than the
full global set for a domain that only needs the JRA55-do wind file at
regional scale.

## 4. Run it

```bash
source /home/jorge/nci/cdx/roundabout/nvhpc_env.sh
cmake -B build_cc70 -S . -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=nvfortran \
    -DRDB_ENABLE_GPU=ON -DRDB_ENABLE_MPI=OFF
cmake --build build_cc70 -j 8

mkdir -p /home/jorge/nci/cdx/data/runs/southern_ocean_1deg
cd /home/jorge/nci/cdx/data/runs/southern_ocean_1deg
ln -s /home/jorge/nci/cdx/data/OM_1deg_southern INPUT
cp /path/to/roundabout/validation_examples/ocean/southern_ocean_1deg/southern_ocean_1deg_wind.nml .
CUDA_VISIBLE_DEVICES=0 OMP_NUM_THREADS=1 /path/to/build_cc70/rdb southern_ocean_1deg_wind.nml \
    > run.log
```

The console prints `[stats] Day ... En ... MaxCFL ...` and the mass/salt/heat
budget lines daily, exactly like the global runs; with the sponge active,
the relaxation shows up as a `src` term on the salt and heat budget lines
(the sponge is a tracked source). The mass line's `out` must stay at
round-off: the sponge edge is a closed wall, and a growing `out` there is
the signature of the defect described in "Stability" just below.

## 5. Stability: the day-573 blow-up was a leaking sponge edge

The first attempt at this namelist (global viscosity, `nu_h = 2000`) ran
for 573 days and then blew up: at outer step 27528 (day 573.5) the
barotropic correction produced ~4000 non-finite faces
(`[nan-catch] ... non-finite BT-correction faces`) and the run aborted on
the vanished-layer invariant (`ERROR STOP I1' violated`). It reproduces
bit for bit, and from a day-540 restart.

**Where.** Not in the ACC, not at the Antarctic coast, not in a western
boundary layer: in the **north ghost row**, beyond the sponge edge at
30 S. The console's SSH minimum (it includes the ghost cells) had been
falling linearly the whole run — −4.9 m at day 61, −8.3 m at day 361,
−10.6 m at day 573, about 1.1 cm/day — at ghost column `i = 49`
(108.5 E). In the last output before the NaN it jumped to −12.2 m while
every interior field stayed smooth. Under `z_fixed` the top four nominal
layers are 2.0 + 2.8 + 3.8 + 4.9 m thick, so a −10 to −13 m surface
empties them in the ghost column; the fillers sat one face away from the
interior and the barotropic correction divided by them. Meanwhile the
interior row next to the edge rose from +0.01 m to +0.59 m, and the
console's `Mass ... out` line grew linearly to −7.4e16 kg (a closed
domain prints round-off there): **about 0.9 Sv of water was flowing
into the domain through the "wall"** and the ghost row, which nothing
refills, was draining.

**Why.** `OBC_SPONGE` is documented as a WALL at its outer face, with an
interior relaxation band. The barotropic substep did close it (its
`select case` default zeroes `vbt`), but every closure in the slow path —
the layer mass-flux zeroing in the continuity, the `uhbt`/`vhbt` wall
reconciliation, the lateral tracer-diffusion walls — compared the raw
tag with `OBC_WALL` and so treated the sponge face as OPEN. The layer
velocity on that face was carried as a flux, drawn from the ghost row
beyond the edge — a reservoir outside the domain and outside every
budget, with nothing to refill it. The 150000 m²/s viscosity only
slowed the leak (0.1 Sv, the ghost SSH at −2.7 m after two years); it
did not close it. The fix is in the model
(`ocean_bc_outer_face_tag`, which maps SPONGE to WALL for those closures;
gate: `test_ocean_dyn_split` `split_obc_sponge_outer_face_is_a_wall`,
which requires a zero-width sponge edge to be bit-identical to a wall).
With it the mass `out` term is round-off (21 kg over two years) and the
ghost SSH stays at 0.

**The Munk-layer warning was a coincidence.** `configure` prints

```
Munk sidewall boundary layer under-resolved: delta_M = (nu_h/beta)^(1/3)
= 46624 m spans only 0.49 cells (dx = 95957 m, beta = 1.97e-10 1/(m*s)
near j=100); ...
```

It is a static worst case over every cell, wet or not, taken at the
largest `beta` — here the northernmost row at 30 S, the sponge row, not
the Antarctic coast. The global run prints the same warning (0.4 cells)
and runs cleanly, and with the leak fixed this domain runs two years at
the same viscosity.

**The viscosity is the global run's.** At the ACC jets the Smagorinsky
term is far below the floor: the surface vorticity in the jet band
(61–45 S) has an rms of 7.7e-7 1/s (99th percentile 2.5e-6), so
`(0.15 Δ)² |D|` with `Δ ≈ 65 km` is 70–250 m²/s and the Laplacian
viscosity there is the 2000 m²/s floor. MOM6's OM_1deg reads
`KH_background_2d.nc` instead, which is 5000 m²/s between about 65 S and
51 S and 0 north of that; roundabout has no reader for it (the global
README names this as a known limit).

**The two-year run** (V100, `CUDA_VISIBLE_DEVICES`-pinned, 2026-09-29):
730 days, no NaN, no truncation, no limiter events; relative mass/salt/heat
closure residuals −1.5e-11 / −2.4e-13 / −1.1e-13 at day 730.

| | this namelist (fixed model) | old 150000 m²/s workaround |
|---|---:|---:|
| Drake Passage, days 30–365 | 139.7 Sv (121.7–153.8) | 134.3 Sv (115.8–149.7) |
| Drake Passage, days 366–730 | 133.2 Sv (116.4–149.3) | 126.1 Sv (109.6–141.5) |
| Drake SSH step, year 1 / year 2 | 1.32 / 1.25 m | 1.30 / 1.19 m |
| circumpolar-mean jet latitude (top 10 m) | 51.7 S | 50.6 S |
| En, last 30 days | 1.39e-3 m²/s² | 6.8e-4 m²/s² |

For comparison, the global wind-forced run gives 148–160 Sv over its
first year. The remaining gap is not the viscosity (it is the same), and
not the sponge's momentum relaxation (`relax_uv = .false.` moves the days
30–180 mean by +0.8 Sv, 142.1 → 142.9). The transport is set by the
WOA13 density field and declines slowly here as that field relaxes with
no buoyancy forcing; the rest of the gap is the domain itself (a walled
edge at 30 S instead of the rest of the ocean) and is not pursued here.

## 6. Diagnostics to report

`../global_1deg/wind_analysis.py`'s section-transport recipe (Drake
Passage = `sum(transport_x * dyt)` over the wet cells of the 67.5 W column
between 70 S and 52 S, divided by 1e6 for Sv) applies unchanged to this
cut — same longitude convention, same `transport_x` diagnostic. An adapted
copy lives in `python_prototypes/southern_ocean_1deg/` (outside this repo)
alongside the run's own analysis and the movie outputs, per the project's
"large media and run-specific analysis don't belong in the model repo"
convention; only the reproducible SCRIPTS (the subsetter, the namelist,
the movie renderer here) are committed to roundabout.

## 7. The movie — `southern_movie.py`

South-polar-stereographic, centred on Antarctica, extending to the domain's
own north edge (`--lat-edge`, default -30, matching the sponge). Main
panel: surface relative vorticity normalised by the local Coriolis
parameter (`zeta / |f|`, the diverging `BALANCE` colour map, symmetric
limits from the run's own 98th-percentile |zeta/f|); a boxed speed inset
(top-left, `INFERNO`) shows the top-10 m current; a day label; and — when
the diagnostic file carries `transport_x` (it does here) — a Drake Passage
transport time strip along the bottom, computed live from the same file
(no separate table needed). Standard library only: reuses
`../global_1deg/global_movie.py`'s `Canvas`/font/colour tables/`encode`
(`ffmpeg`) and `tools/om1deg_prepare_inputs.py`'s NetCDF-3 reader.

```bash
module load netcdf-c   # nccopy flattens the diagnostic file
python3 southern_movie.py \
    RUN/output/southern_ocean_1deg_wind_rank_000000.nc movie_frames \
    --data /home/jorge/nci/cdx/data/OM_1deg_southern --fps60
```

Produces `movie_frames/southern_ocean_1deg_wind.mp4` (+`.gif`), the 60 fps
`minterpolate`-smoothed twin (`_smooth60.mp4`), and a last-frame PNG.

**Known simplification**: the tripolar (here, plain lat-lon) -> image remap
is nearest-cell (flood-filled owner per pixel), not the bilinear-in-
index-space remap `global_movie.BilinearLatLonMap` does for the
equirectangular global movies — that construction is written for a lon-lat
rectangle target; a bilinear polar version needs its own pixel ->
(lat, lon) inverse and is a reasonable follow-up, not done here.

## 8. Known limits

* **1 degree does not resolve mesoscale eddies.** The ACC's transport comes
  through as jets and fronts in the vorticity field (and in the SSH step
  across the passage), not as the eddy field a 1/10 degree or finer
  configuration would show — the same caveat the global 1 degree README
  states for its own basin currents. The movie is "watch the ACC spin up
  and the fronts sharpen", not "watch eddies shed".
* **No sea ice.** The Antarctic coastal boundary is a bare land wall with
  no ice-ocean thermodynamics; near-coastal dynamics (polynyas, coastal
  currents under ice) are absent.
* **Sponge, not a live open boundary.** The north edge relaxes toward a
  FIXED WOA13 January climatology, not toward a time-varying parent-model
  state — there is no seasonal cycle imported from outside the domain, and
  the 10-day relaxation timescale is a modelling choice, not a measured
  quantity. Genuine boundary-forced variability (parent-model eddies,
  Rossby waves entering from the subtropics) cannot appear.
* **The wind is 1958 JRA55-do, cycled.** Like the global wind run, there is
  no heat or freshwater flux and no SST feedback — only wind stress. Any
  seasonal cycle in the diagnostics is the wind's seasonal cycle (plus
  whatever the sponge's fixed target contributes), not a full air-sea
  coupled seasonal cycle.
* **Same western/southern-boundary caveats as the global run** (background
  viscosity below what the grid wants near sharp topography, one-cell
  sills, WOA13 nearest-neighbour thermal-wind striping in the IC) — see
  `../global_1deg/README.md` §5, all inherited unchanged since the physics
  and vertical grid are identical.

## 9. Files

| | |
|---|---|
| `southern_ocean_1deg_wind.nml` | the namelist (§ above) |
| `southern_movie.py` | the polar movie renderer |
| `tools/om1deg_subset.f90` (repo root `tools/`) | the generic NetCDF subsetter used to cut the four inputs |
