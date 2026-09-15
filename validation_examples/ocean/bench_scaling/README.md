# March-in scaling matrix — run kit

Binary: `build_gpu_mpi/rdb` (nvhpc **26.3** — do NOT use 26.5: its
hpcx UCC CUDA plugin has an nvml symbol error at multi-rank AND 26.5
carries the known BT descriptor-memcpy perf regression).
Env: `module load nvhpc/26.3 misc/nvhpc-build/25.5/netcdf-c/
misc/nvhpc-build/25.5/netcdf-fortran/ hdf5`, `OMP_NUM_THREADS=1`.

Four configs, identical apart from `&ocean_bt_nml bt_halo` (verified
by diff): `seamount_bt{0,8}.nml` (512x256x30, 2 d, dt=300, diags off)
and `dg50_bt{0,8}.nml` (300x300x50 double-gyre, 2 d, dt=1200, diags
off).

## Command matrix (from the repo root; 3 runs per cell, take median)

```
cd /home/jorge/nci/cdx/rdb_marchin_wt
OMP_NUM_THREADS=1 mpirun -np 1 bench_scaling/gpu_pin.sh build_gpu_mpi/rdb bench_scaling/seamount_bt0.nml
OMP_NUM_THREADS=1 mpirun -np 2 bench_scaling/gpu_pin.sh build_gpu_mpi/rdb bench_scaling/seamount_bt0.nml
OMP_NUM_THREADS=1 mpirun -np 4 bench_scaling/gpu_pin.sh build_gpu_mpi/rdb bench_scaling/seamount_bt0.nml
# ... same three for seamount_bt8.nml, dg50_bt0.nml, dg50_bt8.nml
```

`gpu_pin.sh` sets `CUDA_VISIBLE_DEVICES` from
`OMPI_COMM_WORLD_LOCAL_RANK` BEFORE MPI_Init (the UCX primary-context
gotcha).  Sanity: during a run `nvidia-smi` must show exactly ONE pid
per device — no {0},{0,1} doubling on device 0.

## What to record per run (all on stdout, rank 0)

1. `Total wall time: ... s` (and `Solver time` if printed).
2. The two comm profiler lines from the final report ("Profiler
   Report: Compute"): `ocean_comms_bt` and `ocean_comms_ml` — columns
   are seconds / call-count / % of time_loop.  The headline metric is
   (ocean_comms_bt + ocean_comms_ml) % of wall.
3. (Optional) the last `halo exch:` counter line — deterministic
   exchange/message counts; at bt_halo=8 `bt_group` should be ~8x
   lower and `bt_u` = 0 vs the bt_halo=0 leg.

## Preliminary numbers (pre-handover, this machine, pinned, quiet)

seamount bt_halo=0: np1 42.6/42.6/42.7 s; np2 41.6/41.9/41.9 s
(3 clean samples each; np2 profiler: ocean_comms_bt 16.0 s = 38.4%,
ocean_comms_ml 1.8 s = 4.3% of the time loop — i.e. ~43% of np2 wall
is comms at bt_halo=0 even after the 3D-batching win).  vs the O4
baseline (same config, pre-arc): np1 43.3 / np2 72.9 / np4 96.6 s.
