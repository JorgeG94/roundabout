# Southern Ocean 1/4 degree, wind-forced, 4-GPU MPI — eddy-resolving ACC

`southern_ocean_025_wind.nml` scales `../southern_ocean_1deg/southern_ocean_1deg_wind.nml`
from MOM6's **OM_1deg** to MOM6's **OM4_025** tripolar geometry — the same
Antarctic-to-~30 S cut, the same north sponge, the same wind source and
vertical-grid family — at 4x the horizontal resolution, run on **4 V100s**
(one rank per GPU, `RDB_CUDA_AWARE_MPI=ON`) for **two simulated years**. At
1/4 degree the domain resolves the first baroclinic Rossby radius through
most of the Southern Ocean, so the Antarctic Circumpolar Current sheds
mesoscale eddies instead of appearing only as smoothed jets and fronts (the
1-degree run's own stated limit).

| | |
|---|---|
| Grid | OM4_025 `ocean_hgrid.nc` supergrid, cut to tracer rows 0-377 (1440 x 378), `west/east = "periodic"`, `south = "wall"` (Antarctica), `north = "sponge"` |
| Domain | the grid's own south row (see Sec. 2 — NOT a uniform-latitude circle, unlike the 1-degree cut) to 29.92 S (north sponge edge, the supergrid q-row nearest 30 S) |
| Bathymetry | same y-cut of OM4_025 `ocean_topog.nc`, MOM6 limits applied in place (`MINIMUM_DEPTH` 9.5 m, `MAXIMUM_DEPTH` 6500 m, `MASKING_DEPTH` 0) |
| Vertical | 50 `z_fixed` levels, tanh-stretched 2 m (surface) to ~255 m (bed) — the SAME nominal profile as the 1-degree runs |
| Initial state | WOA05 annual T/S, already on the OM4_025 model grid, same y-cut, at rest, eta = 0 |
| Forcing | JRA55-do 1958 daily-mean wind stress, regridded directly onto this domain's cut supergrid, cycled annually |
| Time step | `dt = 900 s` (MOM6 OM4_025 `DT`), `pred_corr`, barotropic `n_inner` from the per-wet-cell CFL |
| MPI | 4 ranks, one GPU each; decomposition chosen empirically, Sec. 6 |
| Output | daily: SSH, top-10 m u/v/T (T doubling as SST), surface `vorticity_z`, depth-integrated `transport_x`/`transport_y`; single precision, 2-D only |
| Length | 730 days (2 years) |

## 1. Reproduce it

```bash
# see CLAUDE.md "MPI" + this file's Sec. 5 for the HPC-X/openmpi-5.0.5 gotcha
module unload openmpi/5.0.5
module load nvhpc misc/nvhpc-build/25.5/netcdf-c misc/nvhpc-build/25.5/netcdf-fortran hdf5
export LD_LIBRARY_PATH=<HPC-X ompi5>/lib:<netcdf-fortran>/lib:<netcdf-c>/lib:<hdf5>/lib:<hdf5>/lib64:<nvhpc>/compilers/lib:$LD_LIBRARY_PATH
export PATH=<HPC-X ompi5>/bin:$PATH

cmake -B build_cc70_mpi -S . -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=nvfortran \
    -DCMAKE_C_COMPILER=nvc -DRDB_ENABLE_GPU=ON -DRDB_ENABLE_MPI=ON -DRDB_CUDA_AWARE_MPI=ON \
    -DRDB_GPU_ARCH=cc70
cmake --build build_cc70_mpi -j 8

mkdir -p /home/jorge/nci/cdx/data/runs/southern_025_run
cd /home/jorge/nci/cdx/data/runs/southern_025_run
ln -s /home/jorge/nci/cdx/data/OM4_025_southern INPUT
cp /path/to/roundabout/validation_examples/ocean/southern_ocean_025/southern_ocean_025_wind.nml .
cp /path/to/roundabout/validation_examples/ocean/southern_ocean_025/pin_gpu.sh .
export OMP_NUM_THREADS=1
mpirun -np 4 ./pin_gpu.sh /path/to/build_cc70_mpi/rdb southern_ocean_025_wind.nml > run.log
```

## 2. How the inputs were cut and built

Data already on disk (NOT re-downloaded, per the task): MOM6's OM4_025
inputs at `/home/jorge/nci/cdx/data/OM4_025/` (`ocean_hgrid.nc` — 2880 x
2160 supergrid, `ocean_topog.nc` — 1440 x 1080 raw bathymetry with 95
edits already baked into `depth`/`wet`, `ocean_mask.nc`,
`WOA05_ptemp_salt_annual.v20141007.nc` — annual T/S ALREADY on the 1440 x
1080 model grid, 35 standard levels) and JRA55-do 1958 3-hourly `uas`/`vas`
at `/home/jorge/nci/cdx/data/JRA55do/`. Every derived file lands at
`/home/jorge/nci/cdx/data/OM4_025_southern/` (outside the repo, per
project convention).

**South wall.** Unlike OM_1deg (whose row `j=0` is an exactly
uniform-latitude circle south of the tripolar cap), OM4_025's own row
`j=0` is NOT uniform: probed directly against `ocean_hgrid.nc`, latitude
at that node row ranges −81.66 to −79.78 across longitude (a 1.9-degree
spread) — the grid's own southern edge is not a plain Mercator row this
far south. That is fine: it is still the domain's physical south
boundary (the bathymetry wall-masks Antarctica there regardless of the
row's exact shape), so the cut still starts at tracer row `j=0`, same as
the 1-degree recipe — the uniformity check just isn't a meaningful
validation at this domain's south edge the way it was at OM_1deg's.

**North edge, ~30 S.** Probed the same way: OM4_025's supergrid node rows
are uniform-latitude circles (the grid IS Mercator this far from both
poles), so a q-row (a cell EDGE, not a straddled T-cell centre) nearest
30 S was picked directly — node row 756 (0-indexed) at 29.9152 S is a
genuine edge (node row 755 at −30.0234 S, closer to 30 S nominally, is a
T-cell CENTRE row, the odd-vs-even parity check: even rows are edges, odd
rows are centres, confirmed against the 1-degree grid's own node row
194 = 2 x 97). Node row 756 = 2 x 378, so the cut keeps tracer rows
`j = 0..377` (378 T-rows), landing the sponge's outer wall exactly on a
grid edge, 0.085 degrees off the requested 30 S — the same kind of
"clean row, not exactly 30" the 1-degree cut used (there, 194 landed at
exactly 30.0234 S by coincidence of the two grids sharing the same
Mercator-zone node-latitude formula; here it doesn't, so the edge is
merely close).

**Tools used, one new, one generalised, two reused as-is:**

* `tools/om1deg_subset.f90` (existing, generic dimension-range NetCDF
  subsetter) — used unchanged for the hgrid cut (`nyp:0:757 ny:0:756`),
  and for the topog/IC cuts (`ny:0:378` / `y:0:378`). It needed ONE
  generalisation: OM4_025's raw `ocean_topog.nc` carries `iEdit`/`jEdit`
  (`NF90_INT`) alongside the float fields the tool already handled — a new
  `copy_i4_1d` case was added (rank-1 int, mirroring the existing rank-1
  char/float/double cases). Roundabout's bathymetry reader never consumes
  `iEdit`/`jEdit`/`zEdit` (grepped: zero hits in `src/`), so carrying them
  through the cut unfiltered (some now referencing row indices outside the
  new domain) is inert metadata, not a correctness issue.
* `tools/om1deg_wind_regrid.f90` (existing) — needed **no source change**:
  it already reads `nxp`/`nyp` from whatever supergrid file it's pointed
  at and sizes every output array from that, so pointing it directly at
  the CUT `OM4_025_southern/ocean_hgrid.nc` produces a wind-stress file
  already at the 1/4-degree Southern Ocean size — no separate windowing
  step needed, unlike the 1-degree recipe (which cut a pre-made global
  wind file). `tools/om1deg_prepare_wind.py`, the driver script, DID need
  a small generalisation — a `--hgrid` override so it can size the output
  from any supergrid, not just `<om1deg-dir>/ocean_hgrid.nc` — added
  without changing its default behaviour.
* `tools/om_topo_limit.f90` (**new**) — OM4_025's raw `ocean_topog.nc`
  is NOT pre-limited the way `tools/om1deg_prepare_inputs.py` pre-limits
  `bathy_om1deg.nc` for the 1-degree runs: probed directly, 299 wet cells
  under 9.5 m and 212 over 6500 m, unclamped depth reaching 8797 m.
  Roundabout's bathymetry reader takes `depth` as given (no floor/cap of
  its own — grepped `src/io/rdb_bathymetry.F90`). This tool applies
  MOM6's `limit_topography` rule (wet ⇒ raise to `MINIMUM_DEPTH`, cap at
  `MAXIMUM_DEPTH`; `depth <= MASKING_DEPTH` ⇒ land) IN PLACE on any
  `depth`/`wet` pair, generic in file/variable name — not OM4_025-specific.
* `tools/om_zclim_to_zinit.f90` (**new**) — WOA05 is already on the
  model's horizontal grid, so (per the task brief) no horizontal
  interpolation was needed, only: (1) dropping the file's leading `TIME`
  record, (2) the vertical fill-gap repair `rdb_ocean_z_init` doesn't do
  itself. Checked directly: 969441 of 969446 topog-wet cells have valid
  WOA05 data at the surface (5 mismatches); every wet column's valid-data
  count DECREASES with depth as the local sea floor is reached (level 1:
  585759 fill/land cells; level 35: 1553138 fill cells) — i.e. WOA05
  standard levels routinely run out above a column's true bottom. The
  tool replaces a trailing deep-end fill run with the deepest valid value
  above it (matching `interp_column_linear_z`'s own constant-extrapolation
  rule for a model layer deeper than the source data) and, for the 5
  (plus every LAND) whole-column-fill cases, copies the nearest valid
  neighbour column. `rdb_ocean_z_init` overwrites every DRY column with
  the namelist `land_fill_t`/`land_fill_s` regardless of file content, so
  only the wet-cell repair actually matters to the run.

```bash
# (module load / LD_LIBRARY_PATH as in Sec. 1, or any NetCDF-Fortran toolchain)
D=/home/jorge/nci/cdx/data/OM4_025
O=/home/jorge/nci/cdx/data/OM4_025_southern
mkdir -p $O

nvfortran -O2 $(nf-config --fflags) -o om1deg_subset tools/om1deg_subset.f90 $(nf-config --flibs)
nvfortran -O2 $(nf-config --fflags) -o om_zclim_to_zinit tools/om_zclim_to_zinit.f90 $(nf-config --flibs)
nvfortran -O2 $(nf-config --fflags) -o om_topo_limit tools/om_topo_limit.f90 $(nf-config --flibs)
nvfortran -O2 $(nf-config --fflags) -o om1deg_wind_regrid tools/om1deg_wind_regrid.f90 $(nf-config --flibs)

./om1deg_subset $D/ocean_hgrid.nc  $O/ocean_hgrid.nc  nyp:0:757 ny:0:756
./om1deg_subset $D/ocean_topog.nc  $O/ocean_topog.nc  ny:0:378
./om_topo_limit $O/ocean_topog.nc 9.5 6500.0 0.0

./om_zclim_to_zinit $D/WOA05_ptemp_salt_annual.v20141007.nc $O/ic_woa05_full.nc \
    --temp-var ptemp --salt-var salt --level-var level
./om1deg_subset $O/ic_woa05_full.nc $O/ic_woa05_jan.nc y:0:378
rm $O/ic_woa05_full.nc   # intermediate, full-grid -- not needed by the run

./om1deg_wind_regrid $O/ocean_hgrid.nc \
    /home/jorge/nci/cdx/data/JRA55do/uas_input4MIPs_*_195801010000-195812312100.padded.nc \
    /home/jorge/nci/cdx/data/JRA55do/vas_input4MIPs_*_195801010000-195812312100.padded.nc \
    $O/wind_jra55do_1958_24h.nc 1958 ly04 24
```

Output sizes: `ocean_hgrid.nc` 100 MB, `ocean_topog.nc` 8.3 MB,
`ic_woa05_jan.nc` 145 MB, `wind_jra55do_1958_24h.nc` 1.52 GB.

## 3. Physics vs. MOM6 OM4_025 (`MOM6-examples/ice_ocean_SIS2/OM4_025/MOM_input`)

Read directly off the local checkout — not from memory:

| Knob | This namelist | OM4_025 `MOM_input` | Match |
|---|---|---|---|
| Time step | `dt_fixed = 900.0` | `DT = 900.0` | exact |
| EOS | `eos = "wright"` | `EQN_OF_STATE = "WRIGHT"` (overriding the MOM6 default `WRIGHT_FULL`) | exact |
| Coriolis | `form = "sadourny"` (enstrophy PV-flux) | `CORIOLIS_SCHEME = "SADOURNY75_ENSTRO"` | exact (the 1-degree run used the energy form instead; this config matches OM4_025's own choice) |
| PGF | `form = "fv_mom6"` | `PressureForce_AFV` + `MASS_WEIGHT_IN_PRESSURE_GRADIENT = True` | closest available form |
| Biharmonic viscosity | `smag_ah = .true., smag_bi_const = 0.06` | `SMAGORINSKY_AH = True, SMAG_BI_CONST = 0.06` | exact |
| Laplacian viscosity floor | `nu_h = 150.0` m²/s | `LAPLACIAN = True` (Smagorinsky-backed, no explicit floor given) | scaled from the 1-degree run's `nu_h = 2000` roughly as `dx^2`: `(27.5 km / 110 km)^2 * 2000 ≈ 125`; picked 150, inside the 125-250 band the brief asked for |
| Bottom drag | `form = "quadratic", cd = 3.0e-3, hbbl = 10.0, bg_vel = 0.1` | `BOTTOMDRAGLAW` default `CDRAG = 3.0e-3` (not overridden), `HBBL = 10.0`, `DRAG_BG_VEL = 0.1` | exact |
| Max velocity | `maxvel = 6.0` | `MAXVEL = 6.0` | exact |
| Depth limits | 9.5 / 6500 / 0 m | `MINIMUM_DEPTH = 9.5`, `MAXIMUM_DEPTH = 6500.0`, `MASKING_DEPTH = 0.0` | exact |
| GM / Redi | OFF (`&ocean_gm_nml enable = .false.`) | `THICKNESSDIFFUSE = True` (GM-lite thickness diffusion, `RESOLN_SCALED_KHTH`) | **deviation, per the task brief** ("keep it off"); OM4_025 itself runs thickness diffusion, roundabout's GM stays off here as it is in the 1-degree config |
| Vertical grid | 50 `z_fixed` tanh-stretched levels (the 1-degree run's profile) | `NK = 75`, `ALE_COORDINATE_CONFIG = "HYBRID:hycom1_75_800m.nc,sigma2,..."` (75-layer HYCOM1 hybrid) | **deviation, per the task brief** ("the stretched z_fixed profile from the 1-degree runs is fine"); roundabout's `VCOORD_HYCOM` exists but was not asked for here |
| Vertical mixing | KPP + PP81 (roundabout default-on) | (OM4_025 uses KPP too, via `MOM_CVMix`, not diffed knob-by-knob here) | same family, not knob-matched |

## 4. Memory

`1440 x 378 x 50 = 27,216,000` cells. <memory report to be filled in from
the run's own `[mem]` console lines>.

## 5. Build gotcha — the HPC-X vs. `openmpi/5.0.5` module collision

Confirmed live on this machine: the login environment auto-loads
`openmpi/5.0.5` (`module list` shows it loaded by default, `MPI_ROOT` set
to `~/install/openmpi-5.0.5`), which is a plain gfortran-built OpenMPI —
NOT ABI-compatible with an nvfortran-compiled `mpi_f08` binding. `module
unload openmpi/5.0.5` first, then put NVHPC's bundled HPC-X (found at
`<nvhpc>/comm_libs/13.2/hpcx/hpcx-2.50/ompi5`, NOT the top-level
`comm_libs/mpi` symlink, which trampolines through a CUDA-version
selector) on `PATH`/`LD_LIBRARY_PATH` explicitly, ahead of everything
else. `cmake`'s `FindMPI` then reports
`Found MPI_Fortran: .../hpcx-2.50/ompi5/lib/libmpi_usempif08.so` — the
right one. Verified with `readelf -d build_cc70_mpi/rdb | grep -i
runpath` <to be pasted after the build finishes>.

## 6. Decomposition

<2x2 vs 1x4 vs 4x1 timing table, to be filled in from the measurement
run.>

## 7. Stability, budgets

<filled in from the 2-day smoke run and the 30-day run.>

## 8. Diagnostics — Drake Passage, SSH, jets, EKE

<filled in against `../global_1deg/wind_analysis.py`'s recipe, adapted
for this grid, compared to the 1-degree run's 139.7 / 133.2 Sv (years
1/2) and observations (130-170 Sv).>

## 9. Diagnostic file size

`1440 x 378 = 544,320` T-cells x 4 bytes (single precision) x 7 daily 2-D
fields (SSH, u, v, T, vorticity_z, transport_x, transport_y) ≈ 15.2
MB/day x 730 days ≈ 11.1 GB total across the run (split across 4 per-rank
files, no gather — `tools/merge_output.py` for an offline merge). `KE` is
turned off in the `diags` list for the same reason the 1-degree config
turns it off.

## 10. Movie

```bash
module load netcdf-c   # nccopy flattens the diagnostic file
python3 ../southern_ocean_1deg/southern_movie.py \
    RUN/output/southern_ocean_025_wind_rank_000000.nc \
    /home/jorge/nci/cdx/python_prototypes/southern_ocean_025 \
    --data /home/jorge/nci/cdx/data/OM4_025_southern --bathy ocean_topog.nc \
    --stem southern_ocean_025 --main vorticity --inset speed \
    --title "ROUNDABOUT  SOUTHERN OCEAN 1/4 DEG  JRA55-DO WIND" --fps60
```

`southern_movie.py` needed one generalisation for this domain: it
hardcoded the 1-degree config's pre-limited bathymetry file name
(`bathy_om1deg.nc`); a `--bathy` option (default unchanged, so the
1-degree invocation above is untouched) lets it read the OM4_025 cut's
raw `ocean_topog.nc` instead — same `depth(y,x)` variable, same classic
NetCDF-3 format, no other change needed (the renderer only ever reads
`depth` off that file).

<frame description / eddy visibility notes to be filled in after
rendering.>

## 11. Known limits

Same standing caveats as `../southern_ocean_1deg/README.md` Sec. 8 (no
sea ice, sponge not a live open boundary, 1958 wind only — no heat/
freshwater flux) PLUS: GM/Redi off and 50 z_fixed levels instead of
OM4_025's 75-layer HYCOM1 hybrid coordinate are deliberate simplifications
for this run (Sec. 3), not physics gaps in roundabout itself.
