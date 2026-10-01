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
| Vertical | 50 `z_fixed` levels, tanh-stretched 2 m (surface) to ~258 m (bed) — the SAME nominal profile as the 1-degree runs |
| Initial state | WOA05 annual T/S, already on the OM4_025 model grid, same y-cut, at rest, eta = 0 |
| Forcing | JRA55-do 1958 daily-mean wind stress, regridded directly onto this domain's cut supergrid, cycled annually |
| Time step | `dt = 900 s` (MOM6 OM4_025 `DT`), `pred_corr`, barotropic `n_inner` from the per-wet-cell CFL |
| MPI | 4 ranks, one GPU each; decomposition chosen empirically, Sec. 6 |
| Output | daily: SSH, top-10 m u/v/T/S (T doubling as SST), surface `vorticity_z`, depth-integrated `transport_x`/`transport_y`; single precision, 2-D only |
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

`1440 x 378 x 50 = 27,216,000` cells. From the run's own console lines
(`run.log`, setup banner + the per-GPU allocator report):

```
Memory:   state arrays counted ~ 6.50 GB (ocean state, exact allocatable footprint)
    barotropic (C-grid): 12.89 MB
    layers + tracers: 1.94 GB
    metrics + dyn-core + closures: 4.55 GB
Setup time: 6.997071900000000 s
    enter_data  : 0.8957532000000000 s

Memory budget (ocean state): host allocations across setup ~ 6.62 GB (RSS growth; advisory, not a device estimate)
  device 0 free before mapping: 31.40 GB of 31.72 GB (host-growth upper bound: 6.62 GB)
device memory in use: 6.99 GB of 31.72 GB (ocean state mapped; includes CUDA context ~300-500 MB + other processes)
  state mapped: ~6.67 GB (device-measured; shared-device activity can perturb)
```

At the end of the 730-day run: `device memory in use: 11.71 GB (ocean end
of run)`, grown ~4.73 GB after setup (lazy kernel workspaces + the one-off
CUDA kernel/runtime reservation). Each rank's V100 carries 32 GB, so this
run used well under half a card per GPU — 4-GPU MPI here is bandwidth- and
communication-bound (Sec. 6-7 below), not memory-bound.

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
right one. The 2-year production run (Sec. 7) is itself the empirical
confirmation that the fix holds: 4 ranks over `RDB_CUDA_AWARE_MPI=ON`
GPU-direct halos, 70080 steps, no ABI-mismatch crash, no hang — the
`build_cc70_mpi` tree used for it was not kept around after the run (and
the task brief asked not to rebuild/rerun), so there is no fresh
`readelf -d ... | grep runpath` capture to paste here; the console does
confirm the GPU-direct path was active: `[gpu-bind] halo path: GPU-direct
(CUDA-aware MPI)`.

## 6. Decomposition

Chosen from a 48-step, compute-only timing test (comment baked into
`southern_ocean_025_wind.nml`'s `&mpi_nml`): **4x1 (east-west strips) 8.17 s
/ 1x4 (north-south strips) 9.21 s / 2x2 9.37 s** — 4x1 wins by ~11-13%.
Confirmed by the production run's own startup line:

```
Domain decomposition: 4 ranks as 4 x 1 (px x py) (y-strips) over 1440 x 378 cells
  rank-0 subdomain: 360 x 378 interior cells
```

4x1 keeps every rank's halo a full-height north-south strip (378 rows),
so east-west halo exchange (the periodic seam) stays the only
inter-rank communication that crosses more than one tile boundary; 1x4
and 2x2 instead cut the latitude range into shorter strips/blocks and pay
more small messages per step for the same physical boundary length at
this aspect ratio (1440 wide x 378 tall). The profiler confirms comms is a
real cost at this scale, not in the noise: `ocean_comms_bt` is 13.0% of
total compute time and `ocean_comms_ml` another 2.2% (Sec. 7's profiler
table) — second only to the barotropic solver itself.

## 7. Stability, budgets

**The build.** The 2-year results below (Secs. 7-8) came from this branch
combined with the five fixes they depend on, all now on `main`: the sponge
outer-face wall (#82), the periodic-seam ghost refresh (#84), the free-slip
viscosity and the `z_fixed` + zinit on-target seed (#88), and the
`vorticity_z` diagnostic (#86). This branch alone does not reproduce them;
`main` with this PR (and #83 under it) merged does.

**730 days, 70080 steps, zero NaN / truncation / limiter / error-stop
events** (grepped the full `run.log`). Wall time: 21689.56 s compute /
21696.58 s total (**~6.03 h** on 4 V100s) — **29.71 s/simulated day**,
i.e. **~7.97 simulated years per wallclock day** (SYPD). Throughput:
0.440 Mcells/s per GPU, 1.759 Mcells/s total (4 ranks).

**Energy.** `En` rises from the rest-state `0` through the wind spin-up,
crosses `1e-3 m²/s²` around day ~240 (~8 months, matching the brief),
peaks at **1.163E-03 m²/s²** on day 626, and settles into a
1.05-1.16E-03 m²/s² equilibrated band for the remainder of year 2 — no
drift, no blow-up precursor. `MaxCFL` never exceeds **0.290** (day 136,
an early-spin-up transient) and typically sits 0.05-0.13; day-730 values
are `En 1.079E-03 m²/s²`, `MaxCFL 0.0698`.

**Budgets at day 730** (relative `Error`; `src` is the sponge relaxation
term, which the North BC's interior band is expected to carry since it is
a tracked source, not a leak — see `../southern_ocean_1deg/README.md`
Sec. 5 for why a nonzero `out` there would be the leaking-sponge defect
this model fixed):

| | value | relative Error | `out` | `src` |
|---|---:|---:|---:|---:|
| Mass | 4.301436089E+20 kg | −2.870E−11 | +4.826 kg | 0 (wall, no source) |
| Salt | 1.490118669E+22 | +3.575E−13 | −1.365E+05 | +1.160E+18 (sponge) |
| Heat | 1.052467857E+21 | +1.062E−12 | −6.310E+04 | +1.126E+19 (sponge) |

Mass/salt/heat are closed to 10⁻¹¹-10⁻¹³ relative error over two years on
4 GPUs with GPU-direct MPI halos — the same order the 1-degree run's fixed
sponge-edge model achieves (Sec. 5 there: −1.5E−11/−2.4E−13/−1.1E−13 at
day 730), carrying over cleanly to 4x the resolution and 4 ranks.

**Profiler** (`Profiler Report: Compute`, % of the 21689.56 s compute
total): `ocean_barotropic_solver` 17.0%, `ocean_comms_bt` 13.0%,
`ocean_continuity` 6.8%, `ocean_ale_remap` 5.7%, `ocean_F_slow_assembly`
3.9%, `ocean_vdiff_apply` 3.7%, `ocean_pgf` 2.9%, `ocean_hvisc` 2.7%,
`ocean_comms_ml` 2.2%, `ocean_coriolis_adv` 1.6%, `ocean_vmix_compute`
1.4% — the fast barotropic loop (solver + its own halo comms) is 30% of
total time, same qualitative shape as the single-GPU 1-degree run, with
the BT halo exchange now a first-class cost of being on 4 ranks.

## 8. Diagnostics — Drake Passage, SSH, jets, EKE

Recipe: `python_prototypes/southern_ocean_025/southern_025_analysis.py`
(committed there, not here — see "Commit" below), a 1/4-degree-grid
adaptation of `../global_1deg/wind_analysis.py`'s section-transport
definition (`Σ transport_x · dy` along the wet cells of the 67.5 W column,
70-52 S) — for this cut the nearest column is `i=930` at 67.625 W
(0.125 degree off, same physical cut as the 1-degree run's i-column
within the grids' own resolution), rows covering 69.95-52.10 S, 126 of
150 cells wet.

**Drake Passage transport.**

| | mean | range | cf. 1-degree (old, Sec. 5 there) |
|---|---:|---:|---:|
| Year 1 (days 30-365) | **175.6 Sv** | 159.5-189.4 Sv | 139.7 Sv (121.7-153.8) |
| Year 2 (days 366-730) | **170.0 Sv** | 154.6-182.8 Sv | 133.2 Sv (116.4-149.3) |

Both years show a mild seasonal cycle, ~165-170 Sv trough around the
coarse ~monthly bin covering late June/July, ~177-180 Sv peaks around
bin 2 (Feb) and bin 9 (Sep) — full bin table in
`southern_025_full_output.txt`. The 1-degree numbers above are the
previously-published run on the *old* code (pre the sponge-outer-face and
velocity-form-Laplacian fixes this PR's run carries); the task brief notes
the 1-degree re-measurement on the same fixed code is still in progress,
so this table is the best available like-for-like reference, not a final
one.

**The Drake finding: this is the initial density field, not eddy
spin-up, and not the section's exact longitude.** Three checks:

1. **Section placement is not the explanation.** This cut's nearest-Drake
   column sits at 67.625 W vs the 1-degree cut's 67.5 W — 0.125 degree
   apart, and `southern_025_analysis.py` confirms columns 4 cells either
   side (67.375 W / 67.875 W) agree to within ~1 Sv of the day-730 value
   (169.4 / 168.5 Sv vs 168.9 Sv at i=930) — an order of magnitude too
   small to explain a 30-40 Sv gap.
2. **The transport is set within days, not months.** Early spin-up
   (day 1-10): 47.6, 173.4, 188.2, 179.1, 165.2, 155.1, 165.1, 181.7,
   177.8, 185.2 Sv — the run is already in its year-1 range (159.5-189.4
   Sv) by day 2-3, matching the same mechanism `../global_1deg/README.md`
   Sec. 6.3 documents for the global run ("the transport comes from the
   WOA13 density field, not from the wind: starting from rest with η=0,
   the depth-integrated baroclinic pressure gradient ... forces the
   barotropic mode, which answers within days"). Mesoscale eddies take
   weeks-months to spin up from a resting state; a day-2 plateau is too
   fast to be eddy-driven.
3. **Eddy kinetic energy is real but not growing.** Surface EKE (Sec.
   below) is **1.26E-03 m²/s² (domain) / 1.52E-03 (ACC band)** in year 1
   and **1.34E-03 / 1.67E-03** in year 2 — a ~6-10% year-over-year rise,
   not the multiplicative growth an eddy field still spinning up would
   show. The eddy field is present and resolved (this cut sheds eddies,
   Sec. 10) but it is a secondary, largely-already-equilibrated
   contributor, not the cause of the 1/4-vs-1-degree gap.

The remaining, measured difference between the initial conditions is the
one piece this task's two READMEs both document independently: this run's
IC is WOA05 potential temperature/salinity **already on the OM4_025 model
grid** (no horizontal interpolation, only the vertical fill-gap repair,
Sec. 2 above), while the 1-degree cut's IC is WOA13 **nearest-neighbour**
onto the OM_1deg grid. A sharper, un-smoothed density front across Drake
Passage in the source gives a stronger initial thermal-wind shear and
hence a stronger barotropic adjustment transport — consistent with (1)
and (2) above, though isolating its magnitude from the resolution
difference itself would need a same-grid IC swap this task doesn't run.
What the measurements rule out is section placement and eddy spin-up as
the explanation; what they're consistent with, without fully isolating
it, is the sharper model-grid IC. Both 170-176 Sv-mean values sit at or
slightly above the widely-cited observational range (130-170 Sv) and
MOM6's own global ~148-160 Sv — plausible for an eddy-resolving regional
cut whose transport is set largely by its own (undegraded) density field
plus a real, resolved eddy field, but on the high side of the literature,
consistent with the no-heat/freshwater-flux, single-year-cycled-wind
limitations below.

**SSH step across the ACC** (same section, northernmost minus southernmost
wet cell): Year 1 mean **+1.552 m** (range +1.365 to +1.658), Year 2 mean
**+1.503 m** (range +1.358 to +1.603) — slightly declining between years,
the same slow relaxation-of-the-IC pattern the 1-degree README attributes
to "the WOA13 density field ... declin[ing] slowly ... as that field
relaxes with no buoyancy forcing" (Sec. 5 there); same mechanism here with
WOA05.

**Jet latitude.** Circumpolar-mean (all longitudes) top-10 m zonal speed,
time-meaned per year, restricted to rows that are >=90% wet (excludes the
ragged coastal fringe) and outside the north sponge band: a single
dominant jet core at **55.92 S** in both years (zonal-mean `u` 0.086 m/s
year 1, 0.085 m/s year 2) — no second statistically distinct jet core
passes a 20%-prominence / 2-degree-separation test in either year. This
sits ~4-5 degrees poleward of the 1-degree run's 51.7 S (year 1) / 50.6 S
(year 2) old-code jet — consistent with an eddy-resolving run's ACC core
sitting closer to its "true" (observationally, multi-fronted but
net-equatorward-of-this) position than a 1-degree run's single smoothed
front, though a definitive attribution again needs the 1-degree
re-measurement on the fixed code.

**Surface EKE** (`0.5*<u'^2+v'^2>` about each year's own time mean,
area-weighted over the supergrid's T-cell area):

| | domain mean | ACC band mean (61-45 S) | max |
|---|---:|---:|---:|
| Year 1 | 1.2553E-03 m²/s² | 1.5156E-03 m²/s² | 8.63E-02 m²/s² |
| Year 2 | 1.3409E-03 m²/s² | 1.6716E-03 m²/s² | 8.26E-02 m²/s² |

A real, resolved mesoscale eddy field (max EKE two orders of magnitude
above the domain mean, concentrated in the ACC band) that is
near-equilibrated by year 1 — the 1-degree run, which does not resolve
the first baroclinic Rossby radius here, has no comparable number to
report (its ACC "spins up and fronts sharpen" rather than sheds eddies,
per its own README Sec. 8).

**Vorticity quality** — surface relative-vorticity rms, interior vs.
coastal ring (wet cells within 2 cells, Chebyshev distance, of any land
cell) vs. north sponge band (top 24 rows), masking the diagnostic's
`vorticity_z` land sentinel (`_FillValue = 1e20`):

| day | interior rms (1/s) | coastal rms | coastal/interior | sponge rms | sponge/interior |
|---:|---:|---:|---:|---:|---:|
| 30 | 1.155E-06 | 4.243E-06 | 3.67x | 5.712E-07 | 0.49x |
| 180 | 1.811E-06 | 4.236E-06 | 2.34x | 7.798E-07 | 0.43x |
| 365 | 1.962E-06 | 3.055E-06 | 1.56x | 8.988E-07 | 0.46x |
| 730 | 1.966E-06 | 2.919E-06 | 1.48x | 9.652E-07 | 0.49x |

(Day-180 interior/sponge figures match the earlier 180-day preview
exactly — 1.81E-06 and 0.43x — confirming the ring/band definitions are
stable across the preview and this full run.) The interior rms climbs
~1.7x from day 30 to day 365 as the eddy field spins up and then holds
flat through year 2 (1.962E-06 -> 1.966E-06), tracking the EKE
equilibration above. The coastal ring stays elevated above the interior
throughout (1.5-3.7x) — Sec. 11's staircase-coastline limitation — while
its *ratio* to the (growing) interior falls over time simply because the
denominator grows faster than the coastal numerator. The sponge band
stays suppressed (0.43-0.49x the interior) the whole run, as expected
from a 24-row relaxation-toward-rest band.

## 9. Diagnostic file size

`1440 x 378 = 544,320` T-cells x 4 bytes (single precision) x 8 daily 2-D
fields (SSH, u, v, T, S, vorticity_z, transport_x, transport_y — the
header table above undercounted at 7, missing `salinity`, which the core
diagnostic set always carries alongside `temperature`/`u`/`v`/`SSH`;
corrected here). Actual, measured: **13.13 GB** across the 4 per-rank
output files (`ls -la output/*_wind_rank_*.nc`, 3,283,561,976 bytes each)
over 730 days, i.e. ~18.0 MB/day total — close to the original per-day
estimate, the gap being exactly that extra `salinity` field.
`tools/merge_output.py`
--grid ocean_diag merges the 4 per-rank files (ghost-stripped) into one
global `merged.nc`, 12.7 GB (smaller than the raw sum because the 3-cell
ghost halo on every rank's tile edge is dropped). `KE` is turned off in
the `diags` list for the same file-size reason the 1-degree config turns
it off. Restart files (one 30-day-cadence snapshot kept per rank, the
final one at day 730): 511.1 MB x 4 = 2.04 GB.

## 10. Movie

```bash
module load nvhpc misc/nvhpc-build/25.5/netcdf-c misc/nvhpc-build/25.5/netcdf-fortran hdf5   # nccopy flattens the diagnostic file
python3 ../southern_ocean_1deg/southern_movie.py \
    RUN/merged.nc \
    /home/jorge/nci/cdx/python_prototypes/southern_ocean_025 \
    --data /home/jorge/nci/cdx/data/OM4_025_southern --bathy ocean_topog.nc \
    --stem southern_ocean_025_2yr --main vorticity --inset speed --fps60 --every 2 \
    --title "ROUNDABOUT  SOUTHERN OCEAN 1/4 DEG  JRA55-DO WIND"
```

`--every 2` (365 of the 730 daily frames) was the practical choice here,
not a size one: at this canvas size a 730-frame render takes roughly 2
hours single-threaded (flood-fill polar remap + font rendering per frame,
no caching across frames), while 365 frames takes about 25 minutes and
both `_2yr.mp4` (24.3 MB) and `_2yr_smooth60.mp4` (30.7 MB) land well
under the ~60 MB target either way — `--every 1` would very likely still
fit the budget, just cost 4x the wall time for a video whose 24/60 fps
playback looks the same to the eye at this frame density.

`southern_movie.py` needed one generalisation for this domain: it
hardcoded the 1-degree config's pre-limited bathymetry file name
(`bathy_om1deg.nc`); a `--bathy` option (default unchanged, so the
1-degree invocation above is untouched) lets it read the OM4_025 cut's
raw `ocean_topog.nc` instead — same `depth(y,x)` variable, same classic
NetCDF-3 format, no other change needed (the renderer only ever reads
`depth` off that file). This PR also carries
`feat(southern-movie): selectable speed/ssh/vorticity panels, fix pole
seam` and `fix(southern-movie): treat the vorticity land sentinel as
missing` from the stacked 1-degree PR — the second fix matters here
specifically: this run's `vorticity_z` diagnostic uses `_FillValue=1e20`
over land (Sec. 8), and without that fix the renderer's colour scale
would be blown out by the sentinel.

**Rendered**: `southern_ocean_025_2yr.mp4` (24.3 MB, 24 fps, 365 frames at
`--every 2`) + `_smooth60.mp4` (30.7 MB, 60 fps `minterpolate`) + a GIF +
last-frame PNG, from the full `merged.nc` (both are comfortably under the
~60 MB target without needing a coarser `--every`). Checked the first,
middle and last frames directly:

* **Day 1** — the vorticity panel is almost featureless (near-zero,
  pale), with structure only as thin alternating-sign filaments at a
  handful of topographic pinch points (the Drake Passage mouth, south of
  Africa, south of Australia/Tasmania) — the geostrophic adjustment of
  Sec. 8's "Drake finding" concentrating where the bathymetry forces it,
  before wind spin-up or eddies. The speed inset matches: a faint,
  patchy high-speed ring at the same handful of sites, the rest of the
  domain near zero. The Drake strip reads `+43.6 Sv` at day 1 (same order
  as the analysis script's `+47.6 Sv`, read at a slightly different
  instant within day 1).
* **Day 365** — a dense field of alternating-sign vorticity filaments and
  mesoscale eddies wraps the entire circumpolar band, most intense
  downstream of Drake Passage and south of Africa/Australia — the
  classic standing-meander ACC eddy-shedding sites. The speed inset shows
  a continuous, structured high-speed ring all the way around Antarctica,
  qualitatively thicker and richer than day 1's isolated filaments. Drake
  strip: `+174.0 Sv`, matching Sec. 8's year-1 table.
* **Day 729** (the last `--every 2` sample) — visually as rich and
  turbulent as day 365, not visibly more energetic — the frame-by-frame
  confirmation of Sec. 8's EKE table (only a 6-10% year-over-year rise):
  the eddy field is fully developed well within year 1 and holds, rather
  than continuing to intensify through year 2. Drake strip: `+168.9 Sv`.
* **Antarctica itself is solid, uniform grey in every frame** — the
  `PolarSouthMap` pole-void fix (inherited from the stacked 1-degree PR)
  and the vorticity land-sentinel fix (this PR's own stacked fix,
  `e3552f9a6`) both hold at this resolution: no false open-water wedge
  through the pole, and no `_FillValue=1e20` sentinel blowing out the
  colour scale.

## 11. Known limits

Same standing caveats as `../southern_ocean_1deg/README.md` Sec. 8 (no
sea ice, sponge not a live open boundary, 1958 wind only — no heat/
freshwater flux, no seasonal cycle beyond the wind's own) PLUS, specific
to this 1/4-degree 4-GPU configuration:

* **GM/Redi off.** OM4_025 itself runs `THICKNESSDIFFUSE=True`
  (GM-lite); this config keeps GM off, matching the 1-degree run's
  choice, per the task brief (Sec. 3). At 1/4 degree the domain partly
  resolves its own eddy field instead (Sec. 8's EKE), which is the
  physical role GM parameterises at coarser resolution — the two are not
  simply additive, so this is a real simplification, not a
  resolution-makes-it-moot non-issue.
* **50 `z_fixed` levels, not OM4_025's 75-layer HYCOM1 hybrid.** Per the
  task brief, this run keeps the 1-degree run's vertical profile rather
  than standing up `VCOORD_HYCOM` for this config (Sec. 3); it is an
  available roundabout vertical coordinate, just not exercised here.
* **Coastal roughness: a 2-cell-wide ring of elevated vorticity around
  every coastline** (Sec. 8's vorticity table — 1.5-3.7x the interior
  rms) from the staircase representation of the coast on a Cartesian
  index grid plus unsmoothed OM4_025 bathymetry. The run's own startup
  diagnostics quantify it directly: `z_fixed closed faces: 1663 LIVE
  cells have all four own-layer faces closed` (isolated-water one-cell
  spikes) and `u closed 3832471/7008000, v closed 3831357/7008900`
  (~54.7% of all layer-faces closed by the partial-step mask, dominated
  by the coastline's cell count at this resolution). Topography
  smoothing (not implemented for this run) is the standard remedy; it was
  out of scope here (the task brief asked for the OM4_025 bathymetry
  as-is, limited only by `tools/om_topo_limit.f90`'s MOM6 depth clamps,
  Sec. 2).
* **4-GPU MPI is beta.** This is the first multi-GPU production run in
  this validation suite (the 1-degree run is single-GPU); the
  `RDB_CUDA_AWARE_MPI=ON` GPU-direct halo path, the `px/py` decomposition
  choice (Sec. 6) and the BT-halo comms cost (Sec. 7's profiler table,
  13% of compute time in `ocean_comms_bt` alone) are all new load on this
  configuration relative to the single-GPU 1-degree run, even though this
  run itself completed cleanly with zero NaN/truncation/limiter events
  over two full years (Sec. 7).

## 12. Commit

The reproducible pieces — this namelist, `pin_gpu.sh`, the two `tools/`
generalisations (`om1deg_subset.f90`'s `copy_i4_1d` case,
`om1deg_prepare_wind.py`'s `--hgrid` override), the two new `tools/`
programs (`om_topo_limit.f90`, `om_zclim_to_zinit.f90`) and
`southern_movie.py`'s `--bathy` option + land-sentinel fix — are committed
to `roundabout`. The run-specific analysis script
(`southern_025_analysis.py`), its per-day table (`southern_025_daily.txt`)
and full text output (`southern_025_full_output.txt`) are committed to
`python_prototypes/southern_ocean_025/` (a separate repository). The
rendered movie (MP4s + a last-frame PNG) is left on disk in that same
directory but NOT committed — large media stays out of git per that
repo's own convention (see the `southern_ocean_1deg` commit this PR's run
follows, and `../southern_ocean_1deg/README.md` Sec. 6). `merged.nc`
itself (12.7 GB) stays in the run directory
(`/home/jorge/nci/cdx/data/runs/southern_025_2yr_final/`), outside both
repositories.
