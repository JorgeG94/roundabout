# Barotropic fast loop: launches per substep, and what was fused

Branch `perf/bt-fastloop-fuse` on `b9430f924`.  Bench: `validation_examples/ocean/bench_scaling/seamount_bt0.nml`
(512 x 256 x 30, `nghost = 3`, 576 steps x 2 stages, `n_inner = 65`, so 74 880 substeps per run, walls, bebt = 0.1,
`pred_corr`, `bt_halo = 0`), one V100 (np 1) and two V100s (np 2, HPC-X CUDA-aware MPI, 2 x 1 split), cc70
`-stdpar=gpu -acc=gpu`.  Kernel accounting from `NVCOMPILER_ACC_TIME=1` (the timer serialises, so the totals of those
runs are inflated; only the per-kernel device times and launch counts are used from them).

## 1. Measurement, before (`b9430f924`)

Launches per substep, default path (`bt_halo = 0`, one rank, `pred_corr`, walls), device time per launch from the
accounting run:

| # | source (base) | sweep | device us |
|---|---|---|---|
| 1 | 1189 | Pass 1: eta continuity + KE + zeta (grid = 262: NOT collapsed, one gang per row) | 61 |
| 2 | 1296 | eta swap + zeta wall closure | 5 |
| 3 | 1430 | BEBT stash `ubt_prev <- ubt` | 5 |
| 4 | 1433 | BEBT stash `vbt_prev <- vbt` | 5 |
| 5 | 1439 | Pass 2b: u update | 20 |
| 6 | 1501 | u wall / OBC BC sweep (1-D over j) | 3 |
| 7 | 1725 | Pass 2c: v update | 18 |
| 8 | 1774 | v wall / OBC BC sweep (1-D over i) | 3 |
| 9 | 1951 | time-mean accumulation | 14 |

9 launches, 134 us of device time; plus **two implicit host-to-device copyins per substep** (the `grid` aggregate,
read inside both async regions, was `present_or_copyin`-ed each time: "data copyin transfers: 74880" at both region
lines).  The substep is not launch-latency bound on one rank on this bench: the 61 us Pass 1 (34 % of the loop) was
occupancy-bound, because the `find_uhbt` / `find_vhbt` calls in its body kept the compiler from collapsing the
(j, i) nest (262 blocks of 128 threads, about 20 % occupancy).

np 2 adds, per rank per substep (2 x 1 split, one neighbour): the mid-substep `ocean_halo_face_x(bt_ubt)` (1 pack +
1 unpack launch) and `ocean_halo_packed_wide` (3 pack + 3 unpack launches: one per array), i.e. 8 comms launches
against 9 compute launches, each preceded by an `!$acc wait(1)` that drains the compute queue.  `ocean_comms_bt` was
13.2-14.1 s of the 20-25 s `ocean_barotropic_solver` at np 2 (150 912 calls): the 2-GPU scaling wall is the exchange,
not the kernels.

Other branches (not on the bench path), launches ADDED per substep: wet/dry replaces Pass 1 by 6 sweeps; periodic x:
+1 eta wrap, +2 u wrap, +1 v wrap (+1 per `per_y` equivalent); `per_y`: +1 eta, +1 u, +2 v; tripolar fold (single
px): +1 eta, +1 u, +2 v (and under `px > 1` two host fold exchanges per substep); each open / tidal / Chapman /
clamped edge: +1 eta ghost fill; `bt_halo > 0`: in-loop wide exchanges instead of the per-substep one.  Tides, porous
barriers, wave drag, BT_cont and the upstream-h face add arithmetic to existing sweeps, no launches.

## 2. What changed

1. **No per-substep aggregate copyin.**  The loop-invariant `grid%nghost / nx_phys / ny_phys` and the `bt_work`
   scalars the kernels read (`use_bt_cont_type`, `use_upstream_h_face`, `wd_dry_depth`, `wd_rewet_depth`) are hoisted
   into locals; the kernels no longer name the `grid` aggregate, so the two H2D copyins per substep are gone (the
   "copyin transfers" lines disappear from the accounting).  Bitwise identical.
2. **Pass 1 collapses.**  Two copies of the sweep, differing only in the BT_cont branch.  The default copy has no
   `find_*` call and collapses (grid 262 -> 1061 blocks): 61 -> 42 us.  The flux select keeps the three-way merge of
   the single-sweep form (dropping to a two-way merge changed the GPU answers at round-off through FMA contraction;
   the three-way form is bitwise the old one on the bench).
3. **BEBT stash folded** (2 launches): Pass 2b / 2c copy their own face before overwriting it (own-cell read then
   write, no neighbour reads); the outer-edge faces neither sweep owns are stashed in the eta-swap sweep.
4. **Closed-wall BC folded** (2 launches): per edge a mode (none / zero / clamped / needs-the-sweep, `bc_fuse_mode`).
   When no edge needs the neighbour-reading radiation legs (Flather, tidal, Chapman) and the array edge is not
   marched, the stand-alone u and v BC sweeps are skipped: the wall faces are branch-free selects at the end of
   Pass 2b / 2c (after the wet/dry gate, as before), the u outer-edge columns zero in the eta-swap sweep (Pass 2b
   neither reads nor writes them), the v outer-edge rows zero in the sum sweep (Pass 2b READS them, so they cannot go
   earlier; any periodic wrap or fold between 2c and the sum would copy the unzeroed rows, so periodic / fold keep the
   stand-alone sweep).
5. **Packed exchange in one kernel per side** (`pk_{x,y}_{pack,unpack}3` in `rdb_ocean_halo.F90`): the
   `(eta, ubt, vbt)` call with no optional array packs / unpacks all three in one launch with the SAME buffer layout
   (messages, counts and the wide-halo message-count assertions unchanged).  6 -> 2 launches per exchange.
6. Cleanups requested on the base: `ocean_halo_packed_wide` returns before `pk_ensure` when the halo is not
   initialised; the `wide_halo_guards` text is right at `bt_halo = 0`; stale comments / unused import; stage 7 text.

## 3. Launches per substep after

| rank layout | compute | comms | total | before |
|---|---|---|---|---|
| np 1, walls (bench) | 5 (Pass 1, eta swap, 2b, 2c, sum) | 0 | 5 | 9 |
| np 2 x 1, walls (bench) | 5 | 4 (mid 1 + 1, packed 1 + 1) | 9 | 17 |
| periodic / fold / OBC-radiation edge | 5 + the same wraps / fills as before (the BC sweep(s) stay) | | | |

Device time on the bench: 134 -> 106 us per substep (42 + 6 + 22 + 21 + 15).

## 4. What is left, and why

* Pass 1 -> swap: Pass 1 reads `eta` at neighbours, the swap writes it; a neighbour read of a field written earlier in
  the same fused loop is exactly the unfusable case.  Removing the swap would need ping-pong buffers threaded through
  every reader (wraps, OBC fills, the exchange, the diagnostics): not a launch-level change.
* Pass 2b -> 2c: 2c reads `ubt` at neighbour faces after 2b wrote them (forward-backward ordering, and at np > 1 the
  mid-substep exchange sits between them).  Kept.
* The sum sweep: reading `ubt` after the mid-substep exchange refreshed seam ghosts, and `vbt` after its outer rows are
  zeroed; folding it into the next substep's Pass 1 would sum post-exchange ghost values (different seam-ghost sums).
  Kept.
* The two exchanges per substep (mid `ubt`, end-of-substep triple), each an `!$acc wait(1)` + 1 pack + MPI + 1 unpack:
  the v update reads the u it just wrote, and the fold sits between the two on `px > 1`.  Kept.
* Everything under periodic / fold / radiation OBC / wet-dry / march-in keeps its sweeps: those paths read neighbours
  across the BC or wrap, and none is on the benchmark's critical path.

## 5. Timings (single runs, shared host and shared GPUs: read as indicative)

`ocean_barotropic_solver` (the fast loop plus its comms) and total wall, seconds, base = `b9430f924` binary, new = this
branch, same flags.  The V100s and the host were shared with other jobs (1-min host load 11-30), so absolute numbers
move by a factor of two between sessions; only same-session pairs are comparable.

| session | np | base solver / wall | new solver / wall |
|---|---|---|---|
| early, intermediate binaries | 1 | 14.2-18.3 / 48-54 | 10.1 / 41 (BC folded, no comm change on 1 rank) |
| early | 2 | 24.0-25.0 / 44-45 | 16.8-16.9 / 36.5-36.7 (`ocean_comms_bt` 14.1 -> 8.9) |
| late (GPU shared, load 15-30) | 1 | 34.0 / 89.6 (x2) | 25.3 / 80.4 (x2) |
| late | 2 | 61.8 / 102.8 (x2) | 45.7 / 87.1 (x1) |

Fast-loop time falls by about 26-30 % at np 1 and np 2; total wall by 10-20 %.  The np 2 fast loop still costs about
as much as np 1 (the exchange: `acc wait` + MPI latency + 2 launches, twice per substep), so the two-GPU scaling is
improved in absolute terms but not fixed: what is left there is MPI latency and the forced queue drain, not launches.

## 6. Verification

* GPU restart fields (26 at np 1, 2 x 26 at np 2, every field, every cell, ghosts included): bitwise identical to the
  base binary after every step of the series above.
* gfortran serial `ctest -R rdb`: 290 / 290; GPU serial `ctest -R rdb` (no `-j`, card 2): 289 / 289.
* `test_ocean_decomp_bitid_mpi`, whole and unfiltered, gfortran + OpenMPI: np 2 and np 4 all checks passed.
