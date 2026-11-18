# Benchmarks

`bench_ocean` drives the ocean dyn-core through the engine lifecycle (see its header).

## `bench_ice_kernels` (ice only, no ocean)

Links only `rdb_ice_core_objs` + `rdb_core_objs`. Build with `-DRDB_BUILD_BENCHMARKS=ON`:
`cmake --build <dir> --target bench_ice_kernels`; run `bench_ice_kernels [nx ny [nwarm nrep [nvrows [tol]]]]`
(defaults 512 512 10 50 96). Synthetic seeded state (~30 % ice-free cells), SoA arrays, ncat 5, nilyr 7, nslyr 1,
device-resident. Each loop is one `do concurrent` launch, median of `nrep` after `nwarm` warm-ups, `!$acc wait`
before every timer read. Before timing, each loop is verified once against a host evaluation of the same cell
kernel on the first `nvrows` rows (bitwise with gfortran, 1e-12 on offload / nvfortran builds); a failure exits 1.
`rdb_test_cice_bench_smoke` runs it on a 48 x 32 grid so the loops cannot rot. The atmosphere kernel and the
strength driver are listed as a STUB until those PRs are integrated.

Results (2026-11-13, V100, nvfortran 26.5 `-stdpar=gpu`; gfortran 15.1 `-O3`, one thread; ms per launch;
last column = ms per 25 km step, 1440 x 1080, scaled from the 512 x 512 run for the CPU rows):

| loop | GPU 1440x1080 | GPU 512x512 | gfortran 1 thread 512x512 | nvfortran multicore 1 / 8 threads 512x512 |
|---|---|---|---|---|
| BL99, 6 passes, per (n,j,i) | 425.5 | 71.6 | 2892 | 10113 / 2027 |
| shortwave per (n,j,i) | 10.87 | 1.97 | 241 | 324 / 66 |
| aggregate, packed gather, per (j,i) | 11.02 | 1.99 | 781 | 685 / 88 |
| aggregate, in place, per (j,i) | 9.64 | 1.68 | 792 | 609 / 80 |
| `ice_lateral_melt_step` | 3.18 | 0.60 | 284 | 235 / 32 |
| `ice_err_reduce`, clean | 0.083 | 0.049 | 0.24 | 0.30 / 0.04 |
| surface heat flux + derivative per (n,j,i) | 8.82 | 1.53 | 227 | 153 / 31 |

The BL99 loop is about 40x more expensive per cell-category than its siblings (54.7 ns vs 1.1 - 1.4 ns on the
GPU). `-Minfo=accel -gpu=ptxinfo` for the benchmark file shows why: the BL99 kernel uses 152 registers and a
1232-byte per-thread stack frame (the whole column working set lives in local memory), against 66 registers and
184 bytes for the in-place aggregate; the packed aggregate needs 1648 bytes of stack. To reproduce, compile
`bench_ice_kernels.F90` by hand with the flags of `benchmarks/CMakeFiles/bench_ice_kernels.dir/flags.make` plus
`-Minfo=accel -gpu=ptxinfo`.

Note (G.2): the BL99 row times the unchecked `*_core` routines, as `ice_temperature_changes` calls them, so a
before/after across that change includes the removed per-routine finiteness scans.

Update (2026-11-13, profile built once before the 6 passes, extrapolation to 1440 x 1152, V100, 1440 x 1152 run):
BL99 takes 307.7 ms with the random ~30 %-ice-free mask (70 % active, 52.97 ns per active cell-category) and
47.2 ms with `--band` (a coherent band, first 20 % of rows: 28.51 ns per active cell-category). The aggregate
loops are unchanged (packed 13.2 ms, in place 9.3 ms), shortwave 11.4 / 1.7 ms, surface flux 8.9 / 1.3 ms.
