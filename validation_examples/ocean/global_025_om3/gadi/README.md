# Gadi kit — global 0.25° OM3 ocean-only, 2 nodes x 4 H200 GPUs

A first scaling/feasibility test for `../global_025_om3_wind.nml`: the
ACCESS-OM3 25 km JRA55-do IAF grid (`release-MC_25km_jra_iaf`,
github.com/ACCESS-NRI/access-om3-configs — see
`tmp_local_artifacts/global025/OM3_25KM_INPUTS.md` in the repo for every
source path + citation), ocean-only (no sea ice), wind-stress forced, run
on **2 Gadi nodes x 4 H200 GPUs** (8 MPI ranks, one per GPU) for 10
simulated days.

**This kit was built with no Gadi access.** Every module name, queue
name, path and decomposition choice below is either sourced from the
public ACCESS-OM3 config repo, reused from `../../southern_ocean_025/`'s
proven recipe, or marked `UNKNOWN`/"check:" — see the assumption table at
the bottom before running anything.

## 1. Prepare inputs

```bash
RDB_DATA_DIR=/scratch/<project>/<user>/rdb_data \
    ./prepare_inputs.sh
```

Copies the real OM3 25 km grid/bathymetry/vertical-grid/IC from `/g/data`
(paths in `tmp_local_artifacts/global025/OM3_25KM_INPUTS.md`) into
`$RDB_DATA_DIR/OM3_025/`, then:

- applies MOM6's `limit_topography` bathymetry rule via the new, generic
  `tools/om_prepare_bathy.py` (`netCDF4` + `numpy` — handles OM3's
  NetCDF-4/HDF5 files directly, unlike `tools/om1deg_prepare_inputs.py`'s
  classic-only reader);
- repairs the WOA T/S IC's vertical fill gaps via the EXISTING, already-generic
  `tools/om_zclim_to_zinit.f90` (no source change needed — it already
  auto-detects a 3-D or 4-D source variable and takes
  `--temp-var`/`--salt-var`/`--level-var` overrides);
- regrids JRA55-do 1958 daily-mean wind stress onto OM3's own supergrid via
  the existing `tools/om1deg_prepare_wind.py` + `tools/om1deg_wind_regrid.f90`.

**Fails loud** before copying anything if any `/g/data` source path is
missing, listing every missing path at once.

## 2. Build

```bash
./build.sh
```

`cmake` configure + build with `RDB_ENABLE_GPU=ON -DRDB_GPU_ARCH=cc90
-DRDB_ENABLE_MPI=ON -DRDB_CUDA_AWARE_MPI=ON`, build dir `build_nvhpc_cc90_mpi`
(repo convention: `build_<toolchain>`). Module names inside are the
`southern_ocean_025` recipe's names (nvhpc 25.5 + its netcdf-c/-fortran +
hdf5) — **run `module avail` on Gadi and fix these** if the H200 nodes
sit on a different NVHPC/CUDA stack (likely, since H200 implies a newer
CUDA/driver baseline than the V100 nodes `southern_ocean_025` was built
on).

## 3. Run

```bash
RDB_DATA_DIR=/scratch/<project>/<user>/rdb_data qsub run_2node.pbs
```

8 ranks (`px=4` within a node x `py=2` across the 2 nodes), one GPU each,
`CUDA_VISIBLE_DEVICES` pinned per rank via `pin_gpu.sh` BEFORE `MPI_Init`
(CLAUDE.md "Multi-GPU one node" — required under `RDB_CUDA_AWARE_MPI=ON`
or every rank on a node lands a CUDA context on device 0).
`OMP_NUM_THREADS=1`, walltime 1 h, `restart_interval=0` (no restarts — this
is a feasibility check, not a spin-up).

## What to look at

- **`s/day`** (wall-clock seconds per simulated day) on the console status
  line (`&logging_nml status_interval = 1.0` day) — the scaling number
  this test exists to produce. Compare against the single-GPU
  `double_gyre_mom6.nml` and the 4-GPU `southern_ocean_025` numbers in the
  main `CLAUDE.md` performance table to sanity-check it's in a plausible
  range for 8x the ranks and ~8x the OM4_025 cut's domain.
- **`MaxCFL`** on the same console line — should stay bounded (the
  barotropic substep count `auto_n_inner` derives from it); a climbing
  `MaxCFL` with no corresponding growth in the budgets below is the first
  sign of an instability, not yet a blowup.
- **Salt/heat budgets** (`console` Mass/Salt/Temp/En lines) — with no heat
  or freshwater flux (wind-only forcing), these should close to round-off
  over 10 days, same invariant `global_1deg_wind.nml`'s README documents.
  `En` (kinetic + potential energy) climbing fast from literally nothing
  in the first few days is expected (wind spin-up); climbing *unbounded*
  past day 2-3 is not.
- **The Munk sidewall warning** the local OM4_025 smoke-test validate run
  printed (`tmp_local_artifacts/global025/om4_smoke/validate_out.log` in
  the repo) — "Raise `&ocean_hvisc_nml nu_h` to at least ~3950 m²/s" — is
  about the WESTERN-BOUNDARY viscous layer being under-resolved at
  `nu_h=150`, not refused, but worth checking this OM3-grid run for the
  same symptom (grid-scale noise hugging western boundaries) before
  trusting anything past this feasibility check.

## Assumptions the maintainer MUST verify on Gadi

| # | Assumption | Where | Risk if wrong |
|---|---|---|---|
| 1 | The 5 `/g/data` source paths in `prepare_inputs.sh` still exist and are readable (dated `2026.06.11`/`2026.03.16` snapshots — the config branch may have moved on) | `prepare_inputs.sh` STEP 1 | Hard fail, loud, before any copy — low risk, but re-check the branch for a newer path if it fails |
| 2 | `MINIMUM_DEPTH`/`MASKING_DEPTH` for OM3 (not found in `MOM_input`/`MOM_override` by this survey) — reused OM4_025's 9.5 m / 0.0 m | `OM3_25KM_INPUTS.md`, `prepare_inputs.sh` `OM3_MIN_DEPTH`/`OM3_MASKING_DEPTH` env overrides | Wrong wet/land mask at shallow shelves only — bathymetry still physically sane either way |
| 3 | `woa23_ts_01_mom.nc` variable names (`ptemp`/`salt`/`level` assumed) | `prepare_inputs.sh` `OM3_TEMP_VAR`/`OM3_SALT_VAR`/`OM3_LEVEL_VAR` env overrides | `om_zclim_to_zinit` fails loud (`nf90_inq_varid` miss) — easy to fix, run `ncdump -h` first |
| 4 | `woa23_ts_01_mom.nc` / `ocean_hgrid.nc` / `topog.nc` are NetCDF-4/HDF5 (assumed, not confirmed) | why `tools/om_prepare_bathy.py` uses `netCDF4`+`numpy` instead of extending the classic-only `om1deg_prepare_inputs.py` | If actually classic, everything still works (netCDF4 reads both) |
| 5 | The `qv56` JRA55-do `atmos/` directory's `uas`/`vas` filenames match the glob `tools/om1deg_wind_regrid.f90` expects (`*_input4MIPs_*_YYYYMMDDHHMM-YYYYMMDDHHMM.padded.nc`) | `prepare_inputs.sh` STEP 4 | **Most likely failure point** — the `.padded.nc` suffix implies a preprocessing step was applied to the LOCAL JRA55do copy `southern_ocean_025` used, which the raw `qv56` input4MIPs archive may not have; if the glob misses, run whatever `tools/fetch_om1deg.py --jra-wind` does for that padding against the `qv56` files first, or adjust the glob |
| 6 | `CORIOLIS_SCHEME`/`EQN_OF_STATE`/`BOTTOMDRAGLAW`/`SMAGORINSKY` constants for OM3 (none fetched — reused from `southern_ocean_025`'s OM4_025 values) | `../global_025_om3_wind.nml`, every block flagged `UNKNOWN for OM3` | Physics plausible either way (both are production MOM6 global configs); answers will differ from a true OM3 match |
| 7 | `REENTRANT_X=True` for OM3 (not found explicitly — inferred by analogy with every other global MOM6 tripolar config) | `OM3_25KM_INPUTS.md` | If actually false this is not a real global config — would fail loud immediately (NaN at the date line) |
| 8 | H200 PBS queue name, `ncpus`/`mem` per GPU, project code | `run_2node.pbs` (every `check:` line) | Job submission rejected by PBS — fix before first `qsub`, not a physics risk |
| 9 | `px=4 x py=2` decomposition — chosen to respect the tripolar-fold row/column minimums and map onto 2 nodes x 4 GPUs, but NOT through the empirical timing test `southern_ocean_025/README.md` Sec. 6 used | `../global_025_om3_wind.nml` `&mpi_nml`, `run_2node.pbs` | Correctness is fine either way (any respecting factorisation is bit-identical); only wall-time is at stake |
| 10 | NVHPC/CUDA module versions for the H200 nodes (reused nvhpc 25.5 from the V100-era `southern_ocean_025` recipe) | `build.sh` | H200 likely needs a newer CUDA baseline — `module avail` and fix before building |
