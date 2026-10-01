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
see `southern_movie.py` below. **Current results** (two full years, five
fixes landed since the first clean run — see §5-7): Drake Passage 170.5 /
167.5 Sv (year 1 / year 2), closed mass/salt/heat budgets to 10⁻¹¹-10⁻¹³,
and a side-by-side comparison against the 1/4-degree companion run.

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

## 5. Stability: the day-573 blow-up was a leaking sponge edge (history, fixed)

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

**Superseded (history, not current numbers).** The first clean two-year run
after this fix (V100, `CUDA_VISIBLE_DEVICES`-pinned, 2026-09-29) reported
Drake Passage transport of 139.7 Sv (year 1) / 133.2 Sv (year 2) and a
circumpolar-mean jet at 51.7 S / 50.6 S, against an old `nu_h = 150000`
workaround run at 134.3 / 126.1 Sv. That pair of runs predated four further
fixes that landed afterward — the periodic east-west seam's ghost refresh,
the velocity-form Laplacian/biharmonic going free-slip at land, `z_fixed`
being seeded on its own target under `&ocean_zinit_nml`, and `vorticity_z`
being written as the dynamics' own corner ζ with a proper `_FillValue` —
each of which touches this domain (the periodic seam runs the full
longitude circle here; the free-slip fix changes every land-adjacent
viscous stress; the `zinit`/`z_fixed` fix changes the WOA13 seeding this
domain depends on entirely, since there is no file-backed open boundary to
fall back on). Those old numbers are **not reproduced below** — see
"Results: the two-year run" (§6) for the current, fully-fixed-code numbers,
which land noticeably closer to the 1/4-degree run's.

## 6. Results: the two-year run (all five fixes, current numbers)

**The build.** This branch combined with the five fixes the results depend
on, all now on `main`: the sponge outer-face wall (#82), the periodic-seam
ghost refresh (#84), the free-slip viscosity and the `z_fixed` + zinit
on-target seed (#88), and the `vorticity_z` diagnostic (#86). This branch
alone does not reproduce these numbers; `main` with this PR merged does.

**The run.** `southern_ocean_1deg_wind.nml` on that build, one V100,
single rank, `CUDA_VISIBLE_DEVICES`-pinned, host-staged halo path (no MPI):
730 days, 35040 steps, **3116.4 s wall — 4.3 s/simulated day, ~55 simulated
years per wallclock day (SYPD)**. Grepped the full `run.log` for
`nan`/`truncat`/`error stop`/`non-finite`: **zero hits** — no NaN, no CFL
truncation, no limiter event, no error-stop over the full two years.
Profiler (`Profiler Report: Compute`, % of the 3116.4 s): the barotropic
solver dominates as usual (`ocean_barotropic_solver` 11.5%,
`ocean_continuity` 7.7%, `ocean_ale_remap` 6.5%, `ocean_F_slow_assembly`
4.7%, `ocean_vdiff_apply` 4.3%) — a single-rank run, so unlike the 1/4-degree
run's profile (Sec. 7 there) there is no `ocean_comms_bt`/`ocean_comms_ml`
network cost, just the host-staged halo bookkeeping.

**Budgets at day 730** (relative `Error`; `src` is the sponge relaxation
term, which the north band is expected to carry since it is a tracked
source, not a leak — a nonzero `out` on a closed wall face would be the
leaking-sponge defect Sec. 5 fixed):

| | value | relative Error | `out` | `src` |
|---|---:|---:|---:|---:|
| Mass | 4.277525218E+20 kg | −1.454E−11 | −2.084 kg | 0 (wall, no source) |
| Salt | 1.481850327E+22 | −3.068E−13 | −1.268E+06 | +1.260E+18 (sponge) |
| Heat | 1.047293361E+21 | −2.934E−13 | +1.235E+05 | +1.235E+19 (sponge) |

`En` (day 730) **7.107E-04 m²/s²**, `MaxCFL` **0.0259**. Mass/salt/heat are
closed to 10⁻¹¹-10⁻¹³ relative error over two years — the same order the
1/4-degree run's fixed sponge edge achieves (−2.870E−11/+3.575E−13/
+1.062E−12, `../southern_ocean_025/README.md` Sec. 7) and consistent with
the OBC_SPONGE outer-face fix (Sec. 5) holding at both resolutions.

**Diagnostics.** Reproduced with the SAME script the 1/4-degree run uses,
generalised (2026-10-01) to take the grid/run/diagnostic paths and the
analysis bands as CLI arguments instead of being forked per grid — one
script, `python_prototypes/southern_ocean_025/southern_025_analysis.py`,
now drives both READMEs. It auto-detects and strips the ghost-cell halo
(`ng = (shape − physical_size) / 2` on each horizontal axis, the same
technique `southern_movie.py` already used), so it reads this run's raw,
un-merged `*_rank_000000.nc` (single rank, `ng = nghost = 3`) exactly as it
reads the 1/4-degree run's already-ghost-stripped `merged.nc` (`ng = 0`).
Regression-checked against the published 1/4-degree numbers before use:
reran it against that run's `merged.nc` and reproduced every number in
`../southern_ocean_025/README.md` Sec. 7 (budgets) and Sec. 8 (Drake, SSH,
jet, EKE, vorticity) bit-for-bit (175.6/170.0 Sv, +1.552/+1.503 m, 55.92 S,
the EKE and vorticity tables) before trusting it on the 1-degree grid.

```bash
python3.10 python_prototypes/southern_ocean_025/southern_025_analysis.py \
    --run RUN_DIR --diag RUN_DIR/output/southern_ocean_1deg_wind_rank_000000.nc \
    --data /home/jorge/nci/cdx/data/OM_1deg_southern --bathy bathy_om1deg.nc \
    --out python_prototypes/southern_ocean_1deg/final --stem southern_1deg_2yr \
    --sponge-rows 6
```

(`--sponge-rows 6` is the north-edge analysis band, deliberately a little
wider than the model's own `sponge_width = 5` cells, matching how the
1/4-degree analysis uses 24 analysis rows against that run's `sponge_width
= 20`.) Output: `python_prototypes/southern_ocean_1deg/final/
southern_1deg_2yr_daily.txt` (one row/day) and `southern_1deg_full_output.txt`
(full stdout, committed — see "Commit" below).

**Drake Passage transport** (`sum(transport_x * dyt)`, wet cells of the
67.5 W column, 70-52 S — the section lands on the SAME nominal longitude as
the 1/4-degree run's section, both cuts sharing OM's node-latitude/longitude
convention):

| | mean | range | seasonal cycle |
|---|---:|---:|---|
| Year 1 (days 30-365) | **170.5 Sv** | 153.1-186.0 Sv | trough ~166 Sv around bin 7 (late June/July), peak ~176 Sv bin 9 (Sep) |
| Year 2 (days 366-730) | **167.5 Sv** | 150.6-183.5 Sv | trough ~163 Sv bin 7, peak ~173 Sv bin 9 |

Early spin-up (days 1-10) already ranges 30.7-185.0 Sv and is inside the
year-1 band by day 3-4 — the same fast barotropic-adjustment-to-the-IC
mechanism `../global_1deg/README.md` Sec. 6.3 and the 1/4-degree README
Sec. 8 both document (the transport comes from the WOA13 density field
answering within days, not from wind spin-up or eddy growth over weeks).
Section-placement check: columns 4 cells either side of i=233 (67.5 W) give
165.6/165.3 Sv on the last day vs 165.3 Sv at the home column — sub-Sv,
confirming the number is not sensitive to the exact column.

**This is the headline change from the stale Sec. 5 numbers**: 139.7/133.2
Sv (pre the four later fixes) vs **170.5/167.5 Sv** now — the fixed-code
run sits much closer to the 1/4-degree run's 175.6/170.0 Sv (Sec. 7 below).
Which of the four fixes moved it most is NOT isolated here — this run
carries all four together and no ablation run was done, so attributing the
~31 Sv change to any one of them individually would be speculation beyond
what was measured.

**SSH step across the ACC** (same section, northernmost minus southernmost
wet cell): Year 1 mean **+1.549 m** (range +1.370 to +1.703), Year 2 mean
**+1.520 m** (range +1.350 to +1.636) — slightly declining between years,
the WOA13 density field relaxing slowly under the sponge with no buoyancy
forcing, the same mechanism the 1/4-degree run's SSH step shows with WOA05.

**Jet latitude.** Circumpolar-mean (all longitudes) top-10 m zonal speed,
time-meaned per year, restricted to rows that are >=90% wet and outside the
6-row north sponge analysis band: a single dominant jet core at **53.16 S**
in both years (zonal-mean `u` 0.082 m/s year 1, 0.083 m/s year 2) — no
second statistically distinct jet core passes the 20%-prominence /
2-degree-separation test in either year, same as the 1/4-degree run. This
sits ~2.8 degrees equatorward of the 1/4-degree run's 55.92 S — smaller
than the ~4-5 degree gap the stale (pre-fix) 1-degree numbers showed
against the SAME 1/4-degree reference, consistent with (not conclusive
proof of) the fixed code also narrowing the resolution-driven jet-position
gap, alongside the un-isolated IC difference (Sec. "1° vs 1/4°" below).

**Surface EKE** (`0.5*<u'^2+v'^2>` about each year's own time mean,
area-weighted over the supergrid's T-cell area):

| | domain mean | ACC band mean (61-45 S) | max |
|---|---:|---:|---:|
| Year 1 | 1.0481E-03 m²/s² | 1.0912E-03 m²/s² | 1.8648E-01 m²/s² |
| Year 2 | 9.9639E-04 m²/s² | 1.0474E-03 m²/s² | 1.8187E-01 m²/s² |

Unlike the 1/4-degree run's EKE (which rises ~6-10% year-over-year as its
resolved eddy field equilibrates), this run's EKE is flat-to-slightly
*declining* between years — consistent with this 1-degree grid not
resolving the first baroclinic Rossby radius here. The MAX column is the
tell: this run's max EKE (1.86E-01) is itself ~2x the 1/4-degree run's
(8.63E-02) even though its DOMAIN-MEAN EKE is lower — the variance
concentrates in a handful of large, smoothed, topographically-forced
jets/fronts rather than spreading across a broad field of resolved eddies
(see the movie, Sec. 8).

**Vorticity quality** — surface relative-vorticity rms, interior vs.
coastal ring (wet cells within 2 cells, Chebyshev distance, of any land
cell) vs. north sponge band (top 6 rows), masking the diagnostic's
`vorticity_z` land sentinel (`_FillValue = 1e20`):

| day | interior rms (1/s) | coastal rms | coastal/interior | sponge rms | sponge/interior |
|---:|---:|---:|---:|---:|---:|
| 30 | 3.891E-07 | 1.399E-06 | 3.60x | 2.747E-07 | 0.71x |
| 180 | 6.526E-07 | 1.277E-06 | 1.96x | 3.077E-07 | 0.47x |
| 365 | 6.923E-07 | 1.001E-06 | 1.45x | 3.023E-07 | 0.44x |
| 730 | 6.964E-07 | 9.749E-07 | 1.40x | 3.343E-07 | 0.48x |

The interior rms climbs ~1.8x from day 30 to day 365 (spin-up) then holds
flat through year 2, same qualitative shape as the 1/4-degree run. The
coastal/interior and sponge/interior RATIOS track the 1/4-degree run's
closely (day 730: 1.40x/0.48x here vs 1.48x/0.49x there — Sec. 7 below),
even though the absolute interior rms is ~2.8x smaller here (6.96E-07 vs
1.97E-06 1/s) — a coarser grid damps the resolved vorticity gradients
themselves, but the RELATIVE quality of the coastal ring and the sponge
band (both diagnostics of numerical hygiene, not of the physical eddy
field) is essentially resolution-independent at this grid-to-grid ratio.

## 7. 1° vs 1/4°, side by side

Same fixed code (the Sec. 6 build), same physics
family, same vertical profile, same wind forcing and north-sponge design;
different horizontal resolution AND different initial condition (the one
variable this task does not isolate — see below).

| | 1° (this README) | 1/4° (`../southern_ocean_025/README.md`) |
|---|---:|---:|
| Drake Passage, year 1 / year 2 | 170.5 / 167.5 Sv | 175.6 / 170.0 Sv |
| SSH step, year 1 / year 2 | 1.549 / 1.520 m | 1.552 / 1.503 m |
| Jet latitude | 53.16 S | 55.92 S |
| EKE, ACC band, year 1 / year 2 | 1.09 / 1.05E-03 m²/s² | 1.52 / 1.67E-03 m²/s² |
| ζ rms, day 730: interior | 6.96E-07 1/s | 1.97E-06 1/s |
| ζ rms, day 730: coastal/interior | 1.40x | 1.48x |
| ζ rms, day 730: sponge/interior | 0.48x | 0.49x |
| Cost | 4.3 s/day, 1 V100, 1 rank | 29.71 s/day, 4 V100s, 4 ranks |
| Day-730 relative Error (mass/salt/heat) | −1.5E−11 / −3.1E−13 / −2.9E−13 | −2.9E−11 / +3.6E−13 / +1.1E−12 |

**What the numbers say, measured, not speculated beyond it:**

* **Drake transport and SSH step are now close** (170.5 vs 175.6 Sv year 1,
  a 3% gap; SSH step within 2 mm) — a dramatic narrowing from the stale
  pre-fix 1-degree numbers' 20%+ gap against the same 1/4-degree reference
  (Sec. 6 above). Both runs' transport is set within days by their own
  initial density field (both READMEs document the same early-spin-up
  plateau independently), not by eddy spin-up, so the remaining 3-5 Sv gap
  is consistent with either resolution or the IC difference below — this
  task does not run a same-grid IC swap to isolate which.
* **Jet latitude differs by ~2.8 degrees** (53.16 S vs 55.92 S), smaller
  than the stale pre-fix gap (~4-5 degrees) but not closed. A single-front
  1-degree jet sitting equatorward of an eddying 1/4-degree jet's position
  is the same qualitative pattern the stale numbers showed, just less
  pronounced.
* **EKE is NOT comparable as "the same physics at lower resolution"**: the
  1-degree run does not resolve the first baroclinic Rossby radius here (an
  established limit, Sec. 9 below), so its surface kinetic-energy variance
  reflects jet/front smoothing and topographic forcing, not a damped
  version of the 1/4-degree run's resolved eddy field. The 1/4-degree run's
  EKE is HIGHER on average (1.52-1.67E-03 vs 1.09-1.05E-03) but its MAX is
  LOWER (8.63E-02 vs 1.86E-01) — the 1-degree run concentrates its
  variance in a few strong, coarse jets rather than spreading it across a
  broad eddy field.
* **Vorticity hygiene (coastal ring, sponge band) is resolution-independent
  at this ratio**: both ratios agree to within 0.08 of each other even
  though the absolute vorticity scale differs by ~2.8x. This is evidence
  the coastal-ring and sponge-band elevations are a property of the
  numerics (the staircase coastline, the relaxation band), not an artefact
  that scales with how well eddies are resolved.
* **The initial condition is NOT the same field at two resolutions** — this
  is the one honestly unmeasured confound. The 1-degree run seeds from
  WOA13 decav January T/S, nearest-neighbour onto the OM_1deg model grid
  (Sec. "Inputs" above); the 1/4-degree run seeds from WOA05 annual T/S
  ALREADY on the OM4_025 model grid (no horizontal interpolation, only a
  vertical fill-gap repair — `../southern_ocean_025/README.md` Sec. 2).
  Different climatology (WOA13 vs WOA05), different season (January vs
  annual), different interpolation (nearest-neighbour vs none needed). The
  1/4-degree README's own Sec. 8 attributes part of its higher transport to
  this sharper, un-smoothed IC rather than to its resolved eddy field (EKE
  is near-equilibrated by year 1, not still growing). Disentangling
  resolution from IC sharpness would need a same-grid IC swap (e.g. run the
  1-degree domain from a WOA05-on-OM4_025-regridded-to-OM_1deg field, or
  vice versa); this task measures what is, not what a controlled ablation
  would show.

## 8. The movie — `southern_movie.py`

South-polar-stereographic, centred on Antarctica, extending to the domain's
own north edge (`--lat-edge`, default -30, matching the sponge). Two
independently selectable panels, `--main FIELD` and `--inset FIELD`
(`FIELD` one of `vorticity` / `speed` / `ssh`; **default `--main speed
--inset ssh`**):

* `vorticity` — surface relative vorticity normalised by the local
  Coriolis parameter (`zeta / |f|`, diverging `BALANCE` colour map,
  symmetric limits from the run's own 98th-percentile |zeta/f|).
* `speed` — top-10 m current speed (sequential `INFERNO` colour map,
  `[0, --speed-max]`; `--speed-max` defaults to the run's own
  99.5th-percentile |u|).
* `ssh` — sea-surface-height anomaly about the AREA-WEIGHTED domain mean
  (recomputed every frame; diverging `BALANCE` colour map, symmetric
  `+/- --ssh-max`; `--ssh-max` defaults to the run's own 99.5th-percentile
  |anomaly|).

All colour limits are sampled once across the run (~100 frames) and held
fixed for every frame. The main panel is the full-size disc on the right;
the inset is a proper-size (300 px) second disc in its own column to the
left — never overlapping the main disc — with its own box, colorbar and
label. Every panel's colorbar/label states the field name and its units.
A day label, the title, and — when the diagnostic file carries
`transport_x` (it does here) — a Drake Passage transport time strip along
the bottom, computed live from the same file (no separate table needed).
Standard library only: reuses `../global_1deg/global_movie.py`'s
`Canvas`/font/colour tables/`encode` (`ffmpeg`) and
`tools/om1deg_prepare_inputs.py`'s NetCDF-3 reader.

```bash
module load netcdf-c   # nccopy flattens the diagnostic file
python3 southern_movie.py \
    RUN/output/southern_ocean_1deg_wind_rank_000000.nc movie_frames \
    --data /home/jorge/nci/cdx/data/OM_1deg_southern --fps60
# select the panels explicitly, e.g. the original vorticity + speed view:
python3 southern_movie.py RUN/output/....nc movie_frames \
    --data /home/jorge/nci/cdx/data/OM_1deg_southern \
    --main vorticity --inset speed --fps60
```

Produces `movie_frames/southern_ocean_1deg_wind.mp4` (+`.gif`), the 60 fps
`minterpolate`-smoothed twin (`_smooth60.mp4`), and a last-frame PNG.

**Rendered** (this run, default `--main speed --inset ssh`, all 730 daily
frames, `python_prototypes/southern_ocean_1deg/final/`):
`southern_ocean_1deg_2yr.mp4` (20.9 MB, 24 fps) + `_smooth60.mp4` (60 fps
`minterpolate`) + a GIF. Checked the first, middle and last frames directly:

* **Day 1** — the speed panel is almost dark (near-zero everywhere, INFERNO
  scale maxing at 0.28 m/s), with structure only as a couple of bright
  filaments near topographic pinch points on the NE side of the disc and a
  faint patchy ring elsewhere — the same geostrophic-adjustment-to-the-IC
  signature the 1/4-degree movie shows at day 1, just coarser. The SSH
  inset already shows a clear large-scale dipole (red/high on one side,
  blue/low on the other) from the initial density field, not from the wind
  (which has barely acted for one day). Drake strip reads `+30.7 Sv`,
  matching Sec. 6's day-1 table entry.
* **Day 365** — a continuous, braided bright band of elevated speed wraps
  the entire circumpolar band — the ACC core as a persistent, structured
  jet with filamented fine structure, not a single smooth ring, but also
  not the broad field of discrete eddies the 1/4-degree movie shows at the
  same day. The SSH inset shows a strong, smooth low (blue) centred over
  Antarctica ringed by a high (red) further out — the ACC's dynamic-height
  signature. Drake strip: `+170.6 Sv`, matching Sec. 6's year-1 mean.
* **Day 730** — visually almost indistinguishable from day 365: the same
  jet band, same rough intensity, same filament pattern, not visibly more
  developed — the frame-by-frame confirmation of Sec. 6's EKE table (flat
  to slightly declining year-over-year, unlike the 1/4-degree run's still-
  rising EKE). Drake strip: `+166.7 Sv`, matching Sec. 6's year-2 mean.
* **Antarctica itself is solid, uniform grey in every frame** — the
  `PolarSouthMap` pole-void fix (below) and the vorticity land-sentinel fix
  (movie-side, "treat the vorticity land sentinel as missing") both hold at this resolution: no false
  open-water wedge through the pole.

**Known simplification**: the tripolar (here, plain lat-lon) -> image remap
is nearest-cell (flood-filled owner per pixel), not the bilinear-in-
index-space remap `global_movie.BilinearLatLonMap` does for the
equirectangular global movies — that construction is written for a lon-lat
rectangle target; a bilinear polar version needs its own pixel ->
(lat, lon) inverse and is a reasonable follow-up, not done here. The grid's
i-index IS smoothly periodic in longitude (checked directly against
`ocean_hgrid.nc`) and this flood fill is correctly non-periodic in PIXEL
space, since a polar disc image has no left/right wrap.

**Fixed: the pole-void seam.** This is a regional cut that stops at the
domain's south WALL (`lat ~ -77.8`, well short of the geographic pole), so
the disc's centre — poleward of that wall — has no model data at all, at
any longitude. The nearest-cell flood fill used to fill that gap with
whichever boundary-ring cell's pixel-space breadth-first search reached it
first, which near the projection singularity (`rho -> 0`) is not the same
as shortest physical distance: an open-ocean cell near the Ross Sea sector
(`lon ~ -180`) used to win a wide fan of pixels reaching the disc centre,
cutting a false open-water wedge through what should be — and, at the
model's own south wall, IS — solid Antarctica. `PolarSouthMap` now masks
every pixel closer to the projection centre than the southernmost real
grid row's projected radius to LAND explicitly, instead of leaving it to
flood-fill happenstance.

## 9. Known limits

* **1 degree does not resolve mesoscale eddies.** The ACC's transport comes
  through as jets and fronts in the vorticity field (and in the SSH step
  across the passage), not as the eddy field a 1/10 degree or finer
  configuration would show — the same caveat the global 1 degree README
  states for its own basin currents. The movie is "watch the ACC spin up
  and the fronts sharpen", not "watch eddies shed". Measured, not just
  asserted: Sec. 7's EKE comparison shows this run's domain-mean EKE lower
  but its MAX roughly double the 1/4-degree run's (1.86E-01 vs 8.63E-02
  m²/s²) — the variance concentrates in a handful of coarse, smoothed jets
  rather than spreading across a resolved eddy field.
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

## 10. Files

| | |
|---|---|
| `southern_ocean_1deg_wind.nml` | the namelist (§ above) |
| `southern_movie.py` | the polar movie renderer |
| `tools/om1deg_subset.f90` (repo root `tools/`) | the generic NetCDF subsetter used to cut the four inputs |

## 11. Commit

The reproducible pieces — this namelist and `southern_movie.py` — are
committed to `roundabout` (no `tools/` changes this round; the subsetter was
already in the tree). The run-specific analysis script
(`python_prototypes/southern_ocean_025/southern_025_analysis.py`, now
generalised to drive both this README and the 1/4-degree one from CLI
arguments rather than being forked per grid) and this run's per-day table
(`python_prototypes/southern_ocean_1deg/final/southern_1deg_2yr_daily.txt`)
and full text output (`southern_1deg_full_output.txt`) are committed to
`python_prototypes/` (a separate repository), per the same "large media and
run-specific analysis don't belong in the model repo" convention the
1/4-degree README follows. The rendered movie (MP4s + a last-frame PNG,
`python_prototypes/southern_ocean_1deg/final/`) is left on disk but NOT
committed — large media stays out of git in that repo too. The raw
diagnostic file (881 MB, single rank, no merge needed) stays in the run
directory (`/home/jorge/nci/cdx/data/runs/southern_1deg_2yr_final/`),
outside both repositories.
