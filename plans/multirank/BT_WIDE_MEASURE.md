# Barotropic wide-halo march-in: measurement (stage 1)

Goal of the stage: make `&ocean_bt_nml bt_halo > 0` bitwise equal to the ONE-rank `bt_halo = 0` run
(owned cells of every restart-registry field), so it can become the default.  Plan:
`plans/multirank/BT_WIDE_HALO_PLAN.md` (PR-A + PR-B).  Everything below was measured on beta
`fb1de45b8` plus this branch; "pre-fix" means the branch with only the new test cases.

## Setup

Cases (`tests/mpi/test_ocean_decomp_bitid_mpi.F90`): `bth_<geometry>_w<W>[_b0]`, geometry `seamount`
(closed walls, variable bathymetry, sigma), `channel` (the same seamount, periodic west/east),
`spherical` (lon-lat sector, 0.2 degree, seamount, Wright EOS); `W` in 2, 4, 8; `_b0` = `bebt = 0`
(the default is 0.1); both split schemes; 48 x 32 x 4, `nghost = 3`, `dt = 900 s`, `n_inner = 20`
(fixed: 19 / 9 / 4 in-loop wide exchanges per BT solve at W = 2 / 4 / 8; the case asserts the `bt_group`
counter equals `2*steps*(2 + cycles)` and the mid-substep `bt_u` counter is 0, so it cannot pass with the
march-in off).  Reference: ONE rank at `bt_halo = 0`; decomposed runs use `bt_halo = W` on every
factorisation of 2 and 4 ranks whose smallest tile is at least `nghost + W` wide (W = 8 is skipped on
1x4: 8 rows < 11).  Compare: every restart-registry field, owned cells, bitwise.
Pre-fix runs use `RDB_BITID_BTH_NOCHAIN=1` (`visc_rem_chain = .false.`: the default chain was refused).
Compiler: gfortran 15.1 (FMA contraction off, project flags), OpenMPI 5.0.5.

## Result before the fix (gfortran CPU)

Summary over 72 (np 2) and 96 (np 4) decomposed runs per step count (3 geometries x 3 W x 2 bebt x 2 schemes x layouts):

* every `_b0` run (`bebt = 0`): IDENTICAL, all W, all geometries, both schemes, np 2 and 4, 1 and 48 steps.
* every `W = 2` run: IDENTICAL (with `bebt = 0.1` too: `num_cycles = 1`, nothing to go stale).
* `W = 4, 8` with `bebt = 0.1` on `seamount` and `channel`: FAIL in every run, from the FIRST outer step,
  mismatching cells reach the seam itself (distance 0).  Fields: `ml_h/u/v_layer`, `tracer_c_*`, `*_av_layer`,
  `bt_eta_end`, `vmix_bl_depth` (and `vmix_kv` once).
* the first spherical case (1 degree, spoon basin) was IDENTICAL even without the fix (not
  investigated further: at 100 km spacing the stale-`prev` amplitude is evidently too small to flip a bit).  It was replaced by a 0.2 degree seamount sector that does show it
  (table at the end of this section).

### 1 outer step (2 BT solves), np 2
| case | scheme | layout | status | fields | cells | max abs | cells from seam (min..max) | bt_eta_end cells (max abs) | n_inner |
|---|---|---|---|---|---|---|---|---|---|
| bth_channel_w4 | pred_corr | 1x2 | FAIL | 5 | 3543 | 3.07e-19 | 0..8 | 345 (3.1e-19) | 20 |
| bth_channel_w4 | pred_corr | 2x1 | FAIL | 5 | 9711 | 1.26e-17 | 0..24 | 733 (1.3e-17) | 20 |
| bth_channel_w4 | ssp_rk2 | 1x2 | FAIL | 5 | 9295 | 4.23e-18 | 0..13 | 893 (4.2e-18) | 20 |
| bth_channel_w4 | ssp_rk2 | 2x1 | FAIL | 5 | 21328 | 7.17e-17 | 0..24 | 1402 (7.2e-17) | 20 |
| bth_channel_w8 | pred_corr | 1x2 | FAIL | 5 | 2736 | 1.42e-19 | 0..8 | 311 (1.4e-19) | 20 |
| bth_channel_w8 | pred_corr | 2x1 | FAIL | 5 | 8562 | 6.21e-18 | 0..24 | 664 (6.2e-18) | 20 |
| bth_channel_w8 | ssp_rk2 | 1x2 | FAIL | 5 | 7318 | 6.44e-18 | 0..12 | 770 (6.4e-18) | 20 |
| bth_channel_w8 | ssp_rk2 | 2x1 | FAIL | 5 | 19822 | 1.07e-16 | 0..24 | 1291 (1.1e-16) | 20 |
| bth_seamount_w4 | pred_corr | 1x2 | FAIL | 5 | 3711 | 3.16e-19 | 0..9 | 331 (3.2e-19) | 20 |
| bth_seamount_w4 | pred_corr | 2x1 | FAIL | 5 | 4963 | 1.25e-17 | 0..7 | 387 (1.3e-17) | 20 |
| bth_seamount_w4 | ssp_rk2 | 1x2 | FAIL | 5 | 8843 | 4.34e-18 | 0..11 | 758 (4.3e-18) | 20 |
| bth_seamount_w4 | ssp_rk2 | 2x1 | FAIL | 5 | 10745 | 7.19e-17 | 0..13 | 712 (7.2e-17) | 20 |
| bth_seamount_w8 | pred_corr | 1x2 | FAIL | 5 | 2946 | 2.02e-19 | 0..9 | 282 (2.0e-19) | 20 |
| bth_seamount_w8 | pred_corr | 2x1 | FAIL | 5 | 4247 | 6.05e-18 | 0..8 | 321 (6.1e-18) | 20 |
| bth_seamount_w8 | ssp_rk2 | 1x2 | FAIL | 5 | 6785 | 6.30e-18 | 0..11 | 613 (6.3e-18) | 20 |
| bth_seamount_w8 | ssp_rk2 | 2x1 | FAIL | 5 | 9885 | 1.06e-16 | 0..13 | 655 (1.1e-16) | 20 |

### 1 outer step, np 4
| case | scheme | layout | status | fields | cells | max abs | cells from seam (min..max) | bt_eta_end cells (max abs) | n_inner |
|---|---|---|---|---|---|---|---|---|---|
| bth_channel_w4 | pred_corr | 1x4 | FAIL | 5 | 12357 | 4.35e-18 | 0..7 | 1125 (4.3e-18) | 20 |
| bth_channel_w4 | pred_corr | 2x2 | FAIL | 5 | 11833 | 1.26e-17 | 0..15 | 912 (1.3e-17) | 20 |
| bth_channel_w4 | pred_corr | 4x1 | FAIL | 5 | 19178 | 1.26e-17 | 0..12 | 1399 (1.3e-17) | 20 |
| bth_channel_w4 | ssp_rk2 | 1x4 | FAIL | 5 | 20179 | 5.13e-18 | 0..7 | 1526 (5.1e-18) | 20 |
| bth_channel_w4 | ssp_rk2 | 2x2 | FAIL | 5 | 23006 | 7.34e-17 | 0..15 | 1477 (7.3e-17) | 20 |
| bth_channel_w4 | ssp_rk2 | 4x1 | FAIL | 5 | 26639 | 7.35e-17 | 0..12 | 1536 (7.3e-17) | 20 |
| bth_channel_w8 | pred_corr | 2x2 | FAIL | 5 | 10329 | 6.21e-18 | 0..15 | 854 (6.2e-18) | 20 |
| bth_channel_w8 | pred_corr | 4x1 | FAIL | 5 | 17354 | 8.30e-18 | 0..12 | 1324 (8.3e-18) | 20 |
| bth_channel_w8 | ssp_rk2 | 2x2 | FAIL | 5 | 21458 | 1.07e-16 | 0..15 | 1400 (1.1e-16) | 20 |
| bth_channel_w8 | ssp_rk2 | 4x1 | FAIL | 5 | 26626 | 1.31e-16 | 0..12 | 1536 (1.3e-16) | 20 |
| bth_seamount_w4 | pred_corr | 1x4 | FAIL | 5 | 11732 | 4.35e-18 | 0..7 | 998 (4.3e-18) | 20 |
| bth_seamount_w4 | pred_corr | 2x2 | FAIL | 5 | 8001 | 1.25e-17 | 0..7 | 640 (1.3e-17) | 20 |
| bth_seamount_w4 | pred_corr | 4x1 | FAIL | 5 | 14685 | 1.26e-17 | 0..7 | 1093 (1.3e-17) | 20 |
| bth_seamount_w4 | ssp_rk2 | 1x4 | FAIL | 5 | 19775 | 5.31e-18 | 0..7 | 1415 (5.3e-18) | 20 |
| bth_seamount_w4 | ssp_rk2 | 2x2 | FAIL | 5 | 15878 | 7.18e-17 | 0..11 | 1051 (7.2e-17) | 20 |
| bth_seamount_w4 | ssp_rk2 | 4x1 | FAIL | 5 | 24237 | 7.50e-17 | 0..11 | 1435 (7.5e-17) | 20 |
| bth_seamount_w8 | pred_corr | 2x2 | FAIL | 5 | 6751 | 6.05e-18 | 0..9 | 556 (6.1e-18) | 20 |
| bth_seamount_w8 | pred_corr | 4x1 | FAIL | 5 | 13237 | 8.19e-18 | 0..7 | 1006 (8.2e-18) | 20 |
| bth_seamount_w8 | ssp_rk2 | 2x2 | FAIL | 5 | 14080 | 1.06e-16 | 0..13 | 960 (1.1e-16) | 20 |
| bth_seamount_w8 | ssp_rk2 | 4x1 | FAIL | 5 | 23194 | 1.75e-16 | 0..10 | 1363 (1.8e-16) | 20 |

### 48 outer steps, np 2
| case | scheme | layout | status | fields | cells | max abs | cells from seam (min..max) | bt_eta_end cells (max abs) | n_inner |
|---|---|---|---|---|---|---|---|---|---|
| bth_channel_w4 | pred_corr | 1x2 | FAIL | 11 | 42526 | 3.64e-12 | 0..15 | 1536 (5.1e-13) | 20 |
| bth_channel_w4 | pred_corr | 2x1 | FAIL | 10 | 42373 | 3.52e-12 | 0..24 | 1536 (5.4e-13) | 20 |
| bth_channel_w4 | ssp_rk2 | 1x2 | FAIL | 10 | 43440 | 1.48e-12 | 0..15 | 1536 (4.5e-13) | 20 |
| bth_channel_w4 | ssp_rk2 | 2x1 | FAIL | 11 | 43189 | 1.70e-12 | 0..24 | 1536 (4.6e-13) | 20 |
| bth_channel_w8 | pred_corr | 1x2 | FAIL | 11 | 42613 | 3.52e-12 | 0..15 | 1536 (7.2e-13) | 20 |
| bth_channel_w8 | pred_corr | 2x1 | FAIL | 11 | 42377 | 2.62e-12 | 0..24 | 1536 (5.6e-13) | 20 |
| bth_channel_w8 | ssp_rk2 | 1x2 | FAIL | 11 | 43273 | 1.48e-12 | 0..15 | 1536 (4.7e-13) | 20 |
| bth_channel_w8 | ssp_rk2 | 2x1 | FAIL | 11 | 43522 | 1.59e-12 | 0..24 | 1536 (5.5e-13) | 20 |
| bth_seamount_w4 | pred_corr | 1x2 | FAIL | 11 | 41722 | 5.16e-13 | 0..15 | 1536 (5.2e-13) | 20 |
| bth_seamount_w4 | pred_corr | 2x1 | FAIL | 12 | 41504 | 2.10e-08 | 0..23 | 1536 (6.1e-13) | 20 |
| bth_seamount_w4 | ssp_rk2 | 1x2 | FAIL | 10 | 42997 | 1.48e-12 | 0..15 | 1536 (5.7e-13) | 20 |
| bth_seamount_w4 | ssp_rk2 | 2x1 | FAIL | 10 | 42959 | 8.53e-13 | 0..23 | 1536 (5.1e-13) | 20 |
| bth_seamount_w8 | pred_corr | 1x2 | FAIL | 12 | 41970 | 2.10e-08 | 0..15 | 1536 (5.4e-13) | 20 |
| bth_seamount_w8 | pred_corr | 2x1 | FAIL | 11 | 41872 | 5.04e-13 | 0..23 | 1536 (5.0e-13) | 20 |
| bth_seamount_w8 | ssp_rk2 | 1x2 | FAIL | 10 | 42945 | 1.70e-12 | 0..15 | 1536 (4.6e-13) | 20 |
| bth_seamount_w8 | ssp_rk2 | 2x1 | FAIL | 10 | 43001 | 1.48e-12 | 0..23 | 1536 (5.2e-13) | 20 |

### 48 outer steps, np 4
| case | scheme | layout | status | fields | cells | max abs | cells from seam (min..max) | bt_eta_end cells (max abs) | n_inner |
|---|---|---|---|---|---|---|---|---|---|
| bth_channel_w4 | pred_corr | 1x4 | FAIL | 11 | 43386 | 3.64e-12 | 0..7 | 1536 (5.9e-13) | 20 |
| bth_channel_w4 | pred_corr | 2x2 | FAIL | 11 | 42776 | 2.62e-12 | 0..15 | 1536 (6.5e-13) | 20 |
| bth_channel_w4 | pred_corr | 4x1 | FAIL | 11 | 43193 | 3.52e-12 | 0..12 | 1536 (6.0e-13) | 20 |
| bth_channel_w4 | ssp_rk2 | 1x4 | FAIL | 11 | 44225 | 1.54e-12 | 0..7 | 1536 (4.5e-13) | 20 |
| bth_channel_w4 | ssp_rk2 | 2x2 | FAIL | 10 | 43627 | 1.48e-12 | 0..15 | 1536 (4.4e-13) | 20 |
| bth_channel_w4 | ssp_rk2 | 4x1 | FAIL | 10 | 43875 | 1.48e-12 | 0..12 | 1536 (4.4e-13) | 20 |
| bth_channel_w8 | pred_corr | 2x2 | FAIL | 10 | 43127 | 3.64e-12 | 0..15 | 1536 (5.7e-13) | 20 |
| bth_channel_w8 | pred_corr | 4x1 | FAIL | 11 | 43114 | 3.64e-12 | 0..12 | 1536 (6.3e-13) | 20 |
| bth_channel_w8 | ssp_rk2 | 2x2 | FAIL | 10 | 43828 | 3.13e-12 | 0..15 | 1536 (4.5e-13) | 20 |
| bth_channel_w8 | ssp_rk2 | 4x1 | FAIL | 10 | 43877 | 1.65e-12 | 0..12 | 1536 (4.9e-13) | 20 |
| bth_seamount_w4 | pred_corr | 1x4 | FAIL | 12 | 42779 | 2.10e-08 | 0..7 | 1536 (5.4e-13) | 20 |
| bth_seamount_w4 | pred_corr | 2x2 | FAIL | 10 | 41952 | 4.93e-13 | 0..15 | 1536 (4.9e-13) | 20 |
| bth_seamount_w4 | pred_corr | 4x1 | FAIL | 11 | 42154 | 5.12e-13 | 0..11 | 1536 (5.1e-13) | 20 |
| bth_seamount_w4 | ssp_rk2 | 1x4 | FAIL | 10 | 43635 | 1.48e-12 | 0..7 | 1536 (4.8e-13) | 20 |
| bth_seamount_w4 | ssp_rk2 | 2x2 | FAIL | 10 | 43251 | 1.65e-12 | 0..15 | 1536 (4.9e-13) | 20 |
| bth_seamount_w4 | ssp_rk2 | 4x1 | FAIL | 10 | 43655 | 1.82e-12 | 0..11 | 1536 (4.6e-13) | 20 |
| bth_seamount_w8 | pred_corr | 2x2 | FAIL | 12 | 42292 | 2.10e-08 | 0..15 | 1536 (5.2e-13) | 20 |
| bth_seamount_w8 | pred_corr | 4x1 | FAIL | 11 | 42579 | 6.12e-13 | 0..11 | 1536 (6.1e-13) | 20 |
| bth_seamount_w8 | ssp_rk2 | 2x2 | FAIL | 10 | 43678 | 1.76e-12 | 0..15 | 1536 (5.4e-13) | 20 |
| bth_seamount_w8 | ssp_rk2 | 4x1 | FAIL | 10 | 43607 | 1.76e-12 | 0..11 | 1536 (5.0e-13) | 20 |

### Spherical 0.2 degree seamount, 48 steps, np 2, fix disabled
| case | scheme | layout | status | fields | cells | max abs | cells from seam (min..max) | bt_eta_end cells (max abs) | n_inner |
|---|---|---|---|---|---|---|---|---|---|
| bth_spherical_w4 | pred_corr | 1x2 | FAIL | 12 | 50528 | 1.40e-11 | 0..15 | 1536 (7.9e-13) | 20 |
| bth_spherical_w4 | pred_corr | 2x1 | FAIL | 14 | 54577 | 1.81e-09 | 0..23 | 1536 (1.4e-09) | 20 |
| bth_spherical_w4 | ssp_rk2 | 1x2 | FAIL | 14 | 50987 | 3.49e-11 | 0..15 | 1536 (7.2e-13) | 20 |
| bth_spherical_w4 | ssp_rk2 | 2x1 | FAIL | 14 | 50903 | 3.90e-11 | 0..23 | 1536 (6.7e-13) | 20 |
| bth_spherical_w8 | pred_corr | 1x2 | FAIL | 12 | 50758 | 1.68e-11 | 0..15 | 1536 (8.1e-13) | 20 |
| bth_spherical_w8 | pred_corr | 2x1 | FAIL | 12 | 50659 | 1.52e-11 | 0..23 | 1536 (9.5e-13) | 20 |
| bth_spherical_w8 | ssp_rk2 | 1x2 | FAIL | 14 | 51001 | 3.65e-11 | 0..15 | 1536 (6.8e-13) | 20 |
| bth_spherical_w8 | ssp_rk2 | 2x1 | FAIL | 14 | 50887 | 3.02e-11 | 0..23 | 1536 (6.4e-13) | 20 |

Magnitudes: `max abs` is the largest absolute difference over the field's owned cells (fields are O(1e-3)
to O(1e3), so 1e-17 .. 1e-12 absolute is 1e-10 .. 1e-14 relative).  `vmix_bl_depth` reached 2.1e-8 at
48 steps (a threshold-detecting field).

## What the cause was

RC1 of the plan, confirmed.  The in-loop wide exchange refreshed `eta`, `ubt`, `vbt` only.  With
`bebt > 0` the continuity (Pass 1) also reads `bt_ubt_prev` / `bt_vbt_prev` in the ghost band; they are
stashed from the locally marched `u`, so after an exchange they are stale in the band, the invalid front
restarts `num_cycles - 1` cells deeper every cycle and reaches the owned cells after about three cycles.
Evidence it is exactly that and nothing else: `bebt = 0` twin clean everywhere; `W = 2` (`num_cycles = 1`)
clean; adding `ubt_prev` / `vbt_prev` to the in-loop exchange (two wide singles after the group, only when
`bebt > 0`) makes every case IDENTICAL at np 2 and 4, both schemes, 1 and 48 steps, gfortran CPU, and at
np 2 on the V100.  The exchange point is NOT moved inside the substep.

No further cause showed in closed, periodic or spherical geometry.  Not found is not proven absent: the
cases are 48 x 32 grids, 48 steps, W <= 8.

## Reach (the poisoned-ring proof)

After `k` substeps without an exchange the outermost `k` eta layers, `k + 1` u / v layers and `k` layers
of `ubt_prev` / `vbt_prev` are invalid (one cell per substep).  `RDB_BT_INHALO_POISON=2` (the existing
in-halo hook, now also live under `bt_halo > 0`) overwrites exactly those layers with NaN after every
substep that does not exchange; owned cells stay bitwise equal to the one-rank run (`*_btwide_poison`
ctest legs; `test_ocean_dyn_mpi` poisoned-ring leg).  `+8` extra layers (past the band) makes the
result differ (non-vacuity).  Without the `prev` exchange the poisoned run is NaN in owned cells.

## Configurations measured but not lifted or changed

| configuration (bth_* / bthx_* case) | W | `bebt` | result at np 2 and 4 (gfortran), also V100 np 2 |
|---|---|---|---|
| Flather open east edge, `bthx_flather` | 4, 8 | 0.1, 0 | IDENTICAL, also with the poison hook |
| interior land block straddling the seam, `bthx_island` | 4, 8 | 0.1, 0 | IDENTICAL, also with the poison hook |
| tidal west edge + Flather east, seamount (dev case, since removed) | 4, 8 | 0.1, 0 | FAIL in every run: 14 fields, ~67 000 cells, max abs 3.9, `bt_eta_end` 0.2: independent of `bebt`; RC2 (no model time reaches the wide fast loop, the tidal eta target stays at t = 0) |

Because the tidal edge is accepted today and is O(1) wrong, `validate_config` now refuses
`bt_halo > 0` with a tidal `&ocean_bc_nml` edge (unit test `tidal_edge_bt_halo_refused`); it was the only
silently wrong configuration the measurement found.  Everything else that was refused stays refused
(tides, psurf, porous, closed faces, cavity, wet/dry, Chapman, supergrid/tripolar, `use_cont_type`,
`upstream_h_face`).  `visc_rem_chain` is no longer refused (the fast loop reads only `bt_rem` and the
force, both wide-exchanged each stage); the cases above run at the default chain.  `inhalo_march` with
`bt_halo > 0` is refused by its own rule (test `inhalo_march_refusals`, case 1); before the chain lift that
case was refused for the chain's reason instead, so it did not test the inhalo rule.

## Timing (V100, nvfortran 26.5, HPC-X CUDA-aware, one rank per GPU, GPUs 2 and 3, median of 3)

`seamount_*` namelists of `validation_examples/ocean/bench_scaling` (2 days at dt = 300 s; extreme: 1 day,
60 layers, EPBL, kappa-shear).  `base_*` = the unmodified beta binary; `new_*` = this branch.  `bt8off` =
the namelist as shipped before (`visc_rem_chain = .false.`), `bt8on` = the shipped namelist now
(the pin removed).  The profiler buckets are exclusive; `ocean_comms_bt` is rank 0.

| run | wall s (3 runs) | median wall | median ocean_barotropic_solver | median ocean_comms_bt (calls) | median ocean_comms_ml | final halo exch line |
|---|---|---|---|---|---|---|
| base_bt0_np1 | 43.8, 43.8, 43.7 | 43.8 | 13.98 | 4.28 (76032) | 0.02 | `halo exch: total=87558 bt_group=76033 bt_u=0 ml=3458 c2d=2 c3d=8064 fx2=1 fx3=0 fy2=0 fy3=0 msgs=0` |
| base_bt8off_np2 | 31.9, 32.2, 31.8 | 31.9 | 14.07 | 3.37 (20736) | 1.66 | `halo exch: total=46088 bt_group=20738 bt_u=0 ml=3458 c2d=1155 c3d=9216 fx2=3457 fx3=2304 fy2=3456 fy3=2304 msgs=104854` |
| base_extoff_np2 | 21.8, 22.0, 21.8 | 21.8 | 7.04 | 1.69 (10368) | 1.03 | `halo exch: total=23048 bt_group=10370 bt_u=0 ml=1730 c2d=579 c3d=4608 fx2=1729 fx3=1152 fy2=1728 fy3=1152 msgs=52438` |
| new_bt0_np1 | 44.0, 43.7, 43.6 | 43.7 | 13.93 | 4.26 (76032) | 0.02 | `halo exch: total=87558 bt_group=76033 bt_u=0 ml=3458 c2d=2 c3d=8064 fx2=1 fx3=0 fy2=0 fy3=0 msgs=0` |
| new_bt0_np2 | 45.8, 45.2, 45.4 | 45.4 | 26.28 | 15.72 (150912) | 1.62 | `halo exch: total=170503 bt_group=76033 bt_u=74880 ml=3458 c2d=3 c3d=9216 fx2=1 fx3=3456 fy2=0 fy3=3456 msgs=339859` |
| new_bt8off_np1 | 41.7, 41.6, 41.8 | 41.7 | 13.74 | 2.71 (20736) | 0.02 | `halo exch: total=77191 bt_group=20738 bt_u=0 ml=3458 c2d=1154 c3d=8064 fx2=21889 fx3=0 fy2=21888 fy3=0 msgs=0` |
| new_bt8off_np2 | 33.2, 33.5, 33.5 | 33.5 | 15.66 | 4.80 (20736) | 1.64 | `halo exch: total=82952 bt_group=20738 bt_u=0 ml=3458 c2d=1155 c3d=9216 fx2=21889 fx3=2304 fy2=21888 fy3=2304 msgs=141718` |
| new_bt8on_np1 | 43.5, 43.6, 43.6 | 43.6 | 13.76 | 3.19 (20736) | 0.02 | `halo exch: total=77191 bt_group=20738 bt_u=0 ml=3458 c2d=1154 c3d=8064 fx2=21889 fx3=0 fy2=21888 fy3=0 msgs=0` |
| new_bt8on_np2 | 35.2, 35.0, 34.7 | 35.0 | 15.74 | 4.86 (20736) | 1.63 | `halo exch: total=85256 bt_group=20738 bt_u=0 ml=3458 c2d=1155 c3d=9216 fx2=21889 fx3=3456 fy2=21888 fy3=3456 msgs=144022` |
| new_exton_np1 | 33.6, 33.6, 33.6 | 33.6 | 6.78 | 1.50 (10368) | 0.01 | `halo exch: total=38599 bt_group=10370 bt_u=0 ml=1730 c2d=578 c3d=4032 fx2=10945 fx3=0 fy2=10944 fy3=0 msgs=0` |
| new_exton_np2 | 23.8, 23.6, 24.0 | 23.8 | 7.90 | 2.48 (10368) | 1.02 | `halo exch: total=42632 bt_group=10370 bt_u=0 ml=1730 c2d=579 c3d=4608 fx2=10945 fx3=1728 fy2=10944 fy3=1728 msgs=72022` |

Reading:

* One GPU, `bt_halo = 0`: `ocean_barotropic_solver` 13.98 s (base, runs 13.92 / 13.98 / 14.02) vs 13.93 s
  (13.93 / 13.93 / 14.22): unchanged; the printed `[stats]` / budget lines of the two binaries are
  identical (bit-identical answers).  The 1-GPU bucket of the wide path itself (np1 `bt_halo = 8`) is 13.8 s.
* `bt_halo = 0` at np 2 (45.4 s) is SLOWER than np 1 (43.7 s): 15.7 s in `ocean_comms_bt`, 150 912 calls
  (every substep exchanges twice: group + mid-substep u).  `bt_halo = 8` at np 2 (35.0 s) is 1.30x faster
  than `bt_halo = 0` at np 2 and 1.25x faster than np 1: `ocean_comms_bt` 4.9 s over 20 736 calls.
* The `prev` exchange costs: the same case without it (`base_bt8off_np2`, wrong answers at `bebt = 0.1`)
  is 31.9 s, with it 33.5 s: +5 %.  Exchanges per BT solve go from `10 + 3*floor((n-1)/nc) + 3` to
  `10 + 5*floor((n-1)/nc) + 3` (the plan's formula): `fx2` 3 457 -> 21 889 over the run.  The packed group
  (stage F) removes this extra cost.
* `visc_rem_chain` on vs off at `bt_halo = 8`, np 2: 35.0 s vs 33.5 s (the chain's own producer and BBL
  glue work, not communication).
* `seamount_extreme_bt8` (n_inner = 65), np 2: 23.8 s vs 33.6 s at np 1: 1.41x.

## Open items for the next stages

* RC2: tidal edges (model time into `bt_wide_substep`, a wide `eta_forcing`): measured above.
* RC3 (static fields from their owners, land -0/+0): `bthx_island` is bitwise, so the plan's -0/+0
  concern did not show on a static land block; the formula refill still blocks supergrid / porous /
  closed faces / cavity.
* RC4 (vector-lane round-off): not seen on gfortran (contraction off) at tile widths 16 .. 24; nvfortran
  CPU not measured.
* RC5 (tripolar): stage 4 below.

## Stage 2: static fields from their owners (fix/bt-wide-halo-static)

Before (stage-1 binary, same `bthx_*` cases, gfortran 15.1 CPU, `bth` family at np 1 and 2; the carried
`hvisc_du_visc`/`hvisc_dv_visc` are now compared for the `bth*` cases):

| case | np 1 (1x1, `bt_halo = W` vs 0) | np 2 (2x1, 1x2) |
|---|---|---|
| `bthx_coast_w4/w8/w4_b0` (land peninsulas touching the Flather east edge) | 41 000 .. 58 000 mismatching owned values (52 876 .. 57 884 with the carried tendencies compared), every field | 53 260 .. 58 400 |
| `bthx_ffull_w4/w8` (`flather_form = 'full'`) | 4 (the carried `hvisc_dv_visc`, owned north-face row) | 4 |
| `bthx_tidal`, `bthx_tides`, `bthx_psurf`, `bthx_super`, `bthx_porous`, `bthx_zfixed` | refused at configure | refused |

After: every case above IDENTICAL, W = 4 and 8, both split schemes, `bebt` 0.1 and 0, at np 1 (1x1), 2 (2x1, 1x2)
and 4; the whole `bth*` family 84 / 252 / 302 comparisons at np 1 / 2 / 4 (34 np 4 layouts skipped by the
tile-fit rule: `nghost + W` wider than the tile), also with `RDB_BT_INHALO_POISON=2` and `RDB_BT_RING_POISON=2`.
The reviewer's 31 namelists (`review_btwide/runs`) all reproduce bitwise against one rank at `bt_halo = 0`,
except the four `wdragrp*` (`wave_drag_form = 'roughness_proxy'`), which also differ at `bt_halo = 0`
(`wdragrp0`, `wdragrp210`): a decomposition non-identity of the roughness proxy that has nothing to do with
the march-in.

## Stage 3: periodic axis meeting an open / clamped edge (merged `fix/periodic-edge-decomp-bitid`)

Reviewer findings on 20acd150b, both rooted in `bt_halo = 0` beta:

* periodic x + open/clamped south: the y-pass of the open-edge eta ghost fill skipped the periodic x-ghost
  columns, leaving the corner at the wrap of Pass-1 ghost-row garbage (one rank) vs the exchange's
  (decomposed).  Fixed on beta (fill over every column); the one-rank march-in difference at W = 4 had the same
  cause (the wide kernel runs the same fill).  One-rank answers moved for periodic-x + open/tidal/clamped
  south/north: 6 h seamount, clamped south, h 1.0e-7 m, u 3.6e-6, v 2.2e-6 m/s.
* periodic y with a beta plane or lon-lat / supergrid metrics: the duplicated seam row carries two values of f
  and of the metrics, which copy survives depends on the split (1x2: u 1e-4).  Not well posed (f jumps across
  the wrap): REFUSED at configure on beta (`south/north = periodic` needs `grid_config = cartesian`,
  `coriolis_beta = 0`).  `exchange_corner_host` therefore never sees a y-varying duplicated row.

New wide cases `bthx_perxo`, `bthx_perxc` (W = 4, 8), `bthx_peryo`, `bthx_perxy` (W = 4): bitwise at 1 rank
(`bt_halo = W` against 0) and 2 / 4 ranks.

## Stage 4: packed wide exchange, cadence at the proven reach, `H_ref` once

(Measured on `perf/bt-wide-halo-packed` BEFORE the merge with the static-fields stages above: the tables below are that branch's numbers, stacked on stage 1, `8e5691b86`.  The merged result is stage 5, at the end of this file; stage 4's refusal list (tides, psurf, tidal edges, porous barriers, closed faces, supergrid) was overtaken by stages 2 and 3.)  Same bit-for-bit contract against the one-rank `bt_halo = 0` run.

## What changed

* `ocean_halo_packed_wide`: all arrays of an exchange in ONE message per neighbour per direction pass
  (one X round + one Y round).  Used for the in-loop refresh (eta, ubt, vbt, + `ubt_prev`/`vbt_prev` at
  `bebt > 0`), the per-solve entry set (eta, ubt, vbt, rem_u/v, force_u/v) and the exit set (time means,
  `*_end` snapshots, uhbt/vhbt at normal width).  Buffers reserved by `ocean_halo_reserve`
  (`ocean_halo_pk_grow_count()` stays 0).
* In-loop cadence `bt_wide_cycle_len(ng_w) = ng_w - 1 - WIDE_CYCLE_MARGIN` (`WIDE_CYCLE_MARGIN = 1`), was `bt_halo/2`.
* `ubt_prev`/`vbt_prev` dropped from the entry exchange (the loop reseeds them wholesale when `bebt > 0`,
  never reads them otherwise); static `H_ref` copied + exchanged once.

## Messages-per-solve formula (checked by the counters)

Per BT solve at `bt_halo = W > 0`, `ng_w = nghost + W`, `cad = ng_w - 2`:

    packed exchanges  P = 2 + floor((n_inner - 1) / cad)      (1 entry + in-loop + 1 exit)
    Isends per rank   = P * (active X neighbours + active Y neighbours)   [np 2 wall: P]
    bt_group = bt_u = per-array face singles = 0;  H_ref: one centre_2d per RUN

Before (stage 1): `10 + 5*floor((n-1)/(W/2)) + 3` two-pass array exchanges per solve (bebt > 0); at n = 65, W = 8:
93 rounds; now `2 + floor(64/9) = 9` rounds.  Asserted by `check_bth_counters` (`pk == 2*steps*P`, `bt_group == 0`,
`bt_u == 0`, `c2d == 1`) in every `bth_*` / `bthx_*` case, and by the `test_ocean_dyn_mpi` canary
(counters and np 2 rank-0 Isend count).  Observed: seamount_bt8 (n = 65) np 2: `pk=10368 = 576*2*9`;
Isends 144 022 -> 47 255 per run.

## Where the cadence limit really is (fail-without evidence)

The poison hook's thousands digit stretches the cadence by that many substeps (still poisoning every layer
invalidated so far).  `bt_wide_cycle_len = ng_w - 2` is the new cadence, the proven reach is `ng_w - 1`:

| cadence | result |
|---|---|
| new (`ng_w - 2`), NaN-poisoned | all `bth_*` / `bthx_*` bitwise, CPU np 2/4, GPU np 2 |
| new + 1 (`ng_w - 1`, the proven limit), poisoned | all bitwise (CPU np 2, 4; GPU np 2; TD np 2) |
| new + 2 (`ng_w`, one substep past the reach) | breaks: NaN reaches the owned cells (ctest `*_cycle_break`, WILL_FAIL; GPU leg exit 1; `test_ocean_dyn_mpi`: owned result differs) |

So the margin is exactly one substep, measured, not assumed.  It stays for what 48 x 32 closed / periodic /
spherical cases cannot show.

## Timing (V100, nvfortran 26.5, HPC-X CUDA-aware, one rank per GPU, median of 3)

`final_*` = stage-1 binary, `s2_*` = this branch.  `bt8` = `seamount_bt8` (default chain), `ext` =
`seamount_extreme_bt8`, `nk0` / `nk8` = `double_gyre_nk100_512x256` (nghost = 3, output off, `t_end = 45000 s`,
`bt_halo` 0 / 8), `bt4` / `bt16` = `seamount_bt8` with `bt_halo = 4 / 16`.  EVERY run held both GPU locks,
with `nvidia-smi` and the 1-min load average recorded inside the lock (`.smi` / `.load` files).  The first
matrix (taken while the coordinator's own GPU suites overlapped the cards) and a second (machine load 30-50 from
other projects) were DISCARDED; this third set interleaves `final`/`s2` runs, and the pairs `bt0_np1`,
`ext_np2`, `nk8_np2` were re-taken once more (second set of 3 each, replacing the first).  The host was still
busy (load 4 to 35 during the matrix), so absolute walls are 3 to 8 percent above the quiet-machine numbers of
stage 1 (e.g. `bt0_np1` 43.7 s then, 45.7 s now); compare pairs, not rows from different sessions.

| run | wall s (3 runs) | median wall | median `ocean_barotropic_solver` s | median `ocean_comms_bt` s (calls) | median `ocean_comms_ml` s | final `halo exch:` line |
|---|---|---|---|---|---|---|
| final_bt0_np1 | 45.6, 46.7, 46.6 | 46.6 | 15.59 | 3.31 (76032) | 0.04 | `halo exch: total=87558 bt_group=76033 bt_u=0 ml=3458 c2d=2 c3d=8064 fx2=1 fx3=0 fy2=0 fy3=0 msgs=0` |
| final_bt0_np2 | 45.5, 47.1, 45.6 | 45.6 | 26.48 | 16.03 (150912) | 1.64 | `halo exch: total=170503 bt_group=76033 bt_u=74880 ml=3458 c2d=3 c3d=9216 fx2=1 fx3=3456 fy2=0 fy3=3456 msgs=339859` |
| final_bt8_np1 | 43.8, 43.8, 43.8 | 43.8 | 13.80 | 3.08 (20736) | 0.03 | `halo exch: total=77191 bt_group=20738 bt_u=0 ml=3458 c2d=1154 c3d=8064 fx2=21889 fx3=0 fy2=21888 fy3=0 msgs=0` |
| final_bt8_np2 | 35.4, 35.1, 34.6 | 35.1 | 15.87 | 5.20 (20736) | 1.64 | `halo exch: total=85256 bt_group=20738 bt_u=0 ml=3458 c2d=1155 c3d=9216 fx2=21889 fx3=3456 fy2=21888 fy3=3456 msgs=144022` |
| final_ext_np1 | 33.6, 33.5, 33.6 | 33.6 | 6.75 | 1.21 (10368) | 0.01 | `halo exch: total=38599 bt_group=10370 bt_u=0 ml=1730 c2d=578 c3d=4032 fx2=10945 fx3=0 fy2=10944 fy3=0 msgs=0` |
| final_ext_np2 | 27.4, 24.1, 23.8 | 24.1 | 8.06 | 2.73 (10368) | 1.03 | `halo exch: total=42632 bt_group=10370 bt_u=0 ml=1730 c2d=579 c3d=4608 fx2=10945 fx3=1728 fy2=10944 fy3=1728 msgs=72022` |
| final_nk8_np2 | 50.3, 50.5, 50.5 | 50.5 | 2.09 | 0.73 (2400) | 2.32 | `halo exch: total=28808 bt_group=2402 bt_u=0 ml=3602 c2d=1203 c3d=7200 fx2=3601 fx3=3600 fy2=3600 fy3=3600 msgs=48020` |
| s2_bt0_np1 | 46.6, 45.7, 44.8 | 45.7 | 14.98 | 3.61 (76032) | 0.03 | `halo exch: total=87558 bt_group=76033 bt_u=0 pk=0 ml=3458 c2d=2 c3d=8064 fx2=1 fx3=0 fy2=0 fy3=0 msgs=0` |
| s2_bt0_np2 | 45.8, 45.7, 46.3 | 45.8 | 26.65 | 16.00 (150912) | 1.65 | `halo exch: total=170503 bt_group=76033 bt_u=74880 pk=0 ml=3458 c2d=3 c3d=9216 fx2=1 fx3=3456 fy2=0 fy3=3456 msgs=339859` |
| s2_bt16_np2 | 30.7, 31.2, 30.9 | 30.9 | 11.76 | 1.12 (5760) | 1.63 | `halo exch: total=25353 bt_group=2 bt_u=0 pk=5760 ml=3458 c2d=4 c3d=9216 fx2=1 fx3=3456 fy2=0 fy3=3456 msgs=42647` |
| s2_bt4_np2 | 32.1, 32.2, 32.5 | 32.2 | 13.05 | 2.52 (16128) | 1.65 | `halo exch: total=35721 bt_group=2 bt_u=0 pk=16128 ml=3458 c2d=4 c3d=9216 fx2=1 fx3=3456 fy2=0 fy3=3456 msgs=53015` |
| s2_bt8_np1 | 43.1, 43.4, 43.4 | 43.4 | 13.39 | 2.31 (10368) | 0.02 | `halo exch: total=21896 bt_group=2 bt_u=0 pk=10368 ml=3458 c2d=3 c3d=8064 fx2=1 fx3=0 fy2=0 fy3=0 msgs=0` |
| s2_bt8_np2 | 32.3, 32.6, 31.7 | 32.3 | 12.90 | 1.79 (10368) | 1.69 | `halo exch: total=29961 bt_group=2 bt_u=0 pk=10368 ml=3458 c2d=4 c3d=9216 fx2=1 fx3=3456 fy2=0 fy3=3456 msgs=47255` |
| s2_ext_np1 | 33.6, 33.5, 33.5 | 33.5 | 6.68 | 1.30 (5184) | 0.01 | `halo exch: total=10952 bt_group=2 bt_u=0 pk=5184 ml=1730 c2d=3 c3d=4032 fx2=1 fx3=0 fy2=0 fy3=0 msgs=0` |
| s2_ext_np2 | 23.1, 22.2, 22.3 | 22.3 | 6.30 | 0.89 (5184) | 1.04 | `halo exch: total=14985 bt_group=2 bt_u=0 pk=5184 ml=1730 c2d=4 c3d=4608 fx2=1 fx3=1728 fy2=0 fy3=1728 msgs=23639` |
| s2_nk0_np1 | 86.2, 86.2, 86.1 | 86.2 | 1.11 | 0.26 (6000) | 0.03 | `halo exch: total=15606 bt_group=6001 bt_u=0 pk=0 ml=3602 c2d=2 c3d=6000 fx2=1 fx3=0 fy2=0 fy3=0 msgs=0` |
| s2_nk0_np2 | 51.0, 52.4, 51.8 | 51.8 | 2.47 | 1.53 (10800) | 2.59 | `halo exch: total=28807 bt_group=6001 bt_u=4800 pk=0 ml=3602 c2d=3 c3d=7200 fx2=1 fx3=3600 fy2=0 fy3=3600 msgs=55217` |
| s2_nk8_np1 | 89.8, 88.3, 87.4 | 88.3 | 1.74 | 0.02 (2400) | 0.04 | `halo exch: total=12008 bt_group=2 bt_u=0 pk=2400 ml=3602 c2d=3 c3d=6000 fx2=1 fx3=0 fy2=0 fy3=0 msgs=0` |
| s2_nk8_np2 | 49.8, 50.0, 63.9 | 50.0 | 1.85 | 0.55 (2400) | 2.28 | `halo exch: total=20409 bt_group=2 bt_u=0 pk=2400 ml=3602 c2d=4 c3d=7200 fx2=1 fx3=3600 fy2=0 fy3=3600 msgs=34821` |

Reading:

* `bt_halo = 0` is untouched (`bt0`: `ocean_barotropic_solver` 15.6 s / 15.0 s on np 1, identical counters,
  identical `[stats]` lines as before; only the printed `pk=0` field is new).
* `bt_halo = 8` np 2: 35.1 -> 32.3 s (-8 %), `ocean_comms_bt` 5.20 -> 1.79 s, 20 736 -> 10 368 calls.
  extreme: 24.1 -> 22.3 s (`ocean_comms_bt` 2.73 -> 0.89 s).  `nk100` (ML-comm bound): 50.5 -> 50.0 s.
* Width: np 2 `seamount_bt8` namelist at `bt_halo` 4 / 8 / 16: 32.2 / 32.3 / 30.9 s, `ocean_comms_bt`
  2.52 / 1.79 / 1.12 s, `ocean_barotropic_solver` 13.05 / 12.90 / 11.76 s: the extra wide-band compute is not
  visible at these widths, so the optimum is at least 16 on this 256 x 256 tile (not searched past 16).
* np 1 -> np 2 speedup (1 GPU per rank): `bt0` 1.00x (45.7 -> 45.8 s), `bt8` 1.25x -> **1.34x**
  (43.8/35.1 -> 43.4/32.3), `ext` 1.39x -> **1.50x**, `nk100` 1.66x (`bt_halo = 0`, 86.2 -> 51.8 s) and
  **1.77x** (`bt_halo = 8`, 88.3 -> 49.8 s).

## Where the np 2 wall goes now (`seamount_bt8`, `bt_halo = 8`, median run, 32.3 s)

| part | s | share |
|---|---|---|
| `ocean_comms_bt` | 1.8 | 5.5 % |
| `ocean_comms_ml` | 1.7 | 5.2 % |
| `ocean_barotropic_solver` (the fast loop and its wide copies) | 12.9 | 40 % |
| everything else | 15.9 | 49 % |

The next bottleneck is the fast loop itself, not communication: from np 1 to np 2 the other buckets scale
(`ocean_continuity` 7.2 -> 4.3 s, `ocean_ale_remap` 5.3 -> 2.9, `ocean_vdiff_apply` 3.8 -> 2.1) while
`ocean_barotropic_solver` stays at 13.4 -> 13.3 s: its kernels run on 2-D fields of 131 k cells (65 k per rank) and
are launch-latency bound, so halving the tile does not halve them.  Fewer, fatter fast-loop kernels (fusing the
Pass 1 / 2b / 2c / accumulation launches) is the lever, then the ML exchanges (`ocean_comms_ml`, 3 458 calls, flat at
about 1.7 s).

## bt_halo = 0: can the packed group serve it?

Feasible and answer-preserving, not done here: the per-substep group exchange (eta, ubt, vbt: three sequential
two-pass rounds) is exactly a `ocean_halo_packed_wide(oh_nghost, ...)` call (same strips, same order), cutting
the per-substep rounds from 4 (3 + the mid-substep u exchange) to 2.  The mid-substep u exchange sits between
Pass 2b and 2c and cannot join it.

## Stage 5: the two lines merged (`perf/bt-wide-halo-packed-static`)

Stage 2 / 3 (static fields from their owners, tides, tidal / Flather edges, porous, closed faces, supergrid,
periodic-edge mixes) were proven with the old cadence (`bt_halo/2`) and one exchange per array; stage 4 (packed
exchange, cadence `ng_w - 2`) was proven on the stage-1 refusal set.  The merge carries both, and the questions
the two authors could not see were answered by running the poison / cadence legs over EVERY `bth_*` and `bthx_*`
case (48 x 32, `n_inner = 40`, gfortran 15.1 + OpenMPI 5.0.5, CPU):

| leg | np 1 (1x1 `bt_halo = W` vs 0) | np 2 | np 4 (38 layouts skipped by the tile-fit rule) |
|---|---|---|---|
| plain | 96 identical | 288 | 346 |
| `RDB_BT_INHALO_POISON=2` + `RDB_BT_RING_POISON=2` | 96 | 288 | 346 |
| cadence + 1 (`RDB_BT_INHALO_POISON=1002`, the proven limit) | 96 | 288 | 346 |
| whole unfiltered `test_ocean_decomp_bitid_mpi` | - | 770 | 1004 |

GPU (V100, nvfortran 26.5, HPC-X CUDA-aware, GPUs 2 and 3, `RDB_ENABLE_GPU=ON`): the `bth` family plain and with
`RDB_BT_INHALO_POISON=2` + `RDB_BT_RING_POISON=2`: np 1 96 / 96 identical, np 2 288 / 288 identical, no failures; cadence + 1 at np 2: 288 identical.

Cadence + 2 (`RDB_BT_INHALO_POISON=2002`) BREAKS (NaN reaches the owned cells) on all 18 case families run one by
one at np 2 (`bth_seamount_w4`, `_b0`, `bth_channel_w4`, `bth_spherical_w4`, and `bthx_` flather, island, coast,
ffull, tidal, tides, psurf, super, porous, zfixed, perxo, perxc, peryo, perxy; each is a ctest `WILL_FAIL` leg).
So the reach argument (one cell per substep, `cadence <= ng_w - 1` exact) holds with the open / tidal / Flather
edge fills, the periodic wraps, the closed-face and porous widths, the supergrid metrics and the wide
`eta_forcing`, and it is tight there too: the one substep of margin is the whole of the slack.

What the merge had to change:

* **Counter formula.**  Packed exchanges per solve are unchanged (`2 + floor((n_inner-1)/(ng_w-2))`, asserted
  `pk == 2*steps*P`).  The static stage adds traffic that is not packed and is now asserted too:
  centre-2D = 1 (`H_ref`, once) + 2*steps under a body tide / psurf (the wide `eta_forcing`, once per solve)
  + steps under psurf (its own per-step seam refresh, present at `bt_halo = 0`); x-face-2D = y-face-2D = 2*steps
  under porous barriers / closed faces (the two recomputed BT widths follow the packed set in one face exchange
  each: the packed call takes at most three x / y arrays and the entry set already uses them), 0 otherwise;
  `bt_group = bt_u = 0` everywhere.  Periodic cases (a rank its own neighbour at one rank) satisfy the same formula.
* **Ring poison.**  `RDB_BT_RING_POISON` poisoned `ubt_prev` / `vbt_prev` and `H_ref` rings every solve.  `H_ref`
  is exchanged once, so poisoning it again on later solves garbled a seam ring nothing refills: it is poisoned
  with the first solve only.  `ubt_prev` / `vbt_prev` left the entry exchange: at `bebt = 0` the kernel still
  evaluates `(1+bebt)*ubt - bebt*ubt_prev`, so a NaN ring was `0*NaN` in a cell with no physical dependency
  (NaN in the owned cells at step 1 of `bth_seamount_w4_b0`, np 2); they are no longer ring-poisoned.  The
  physical rings of every other wide array, `ubt`, `vbt`, `eta`, rem, force, the widths and `eta_forcing`, stay
  poisoned and bitwise.
* `H_ref` once is sound here: the only thing that moves it is the ice-shelf cavity, which `bt_halo > 0` refuses.
* The wider exit set (time means, `*_end`, uhbt / vhbt) is bitwise with periodic wraps and open-edge ghosts: a
  physical edge is never filled by the exchange (only seams and periodic wraps are), exactly as the one-rank
  `bt_halo = 0` path leaves it, and every `bthx_*` case compares the following steps too.
* New ctest legs: `*_btwide_cycle_break_<class>` (14 classes).

Timing: the stage 4 tables above were measured on the packed branch before the merge; not re-measured.

## Stage 4b: the tripolar north fold (fix/bt-wide-halo-fold; proven with the OLD cadence `bt_halo/2` and per-array exchanges, BEFORE the merge with the packed stage; see stage 6 for the merged result)

Mechanism.  The wide fast loop already carried the fold inline (`bc%north_fold`, the fold index
constants built from `grid%nghost`, which on the wide grid is `nghost + bt_halo`), so at `px = 1` the
wide band is folded by the local kernels over all `nghost + W` rows with no further change.  The
fold is applied EVERY substep, at the same two points per substep as the `bt_halo = 0` loop (eta after
Pass 1, ubt after Pass 2b, vbt after Pass 2c), not at the wide exchange only: the fold-line row (`v`,
storage row `ng+nj+1`) is projected from its mirror each substep (west half := minus the east half), and
the mirror point sits in the same tile (`px = 1`) or another one (`px > 1`).  A band refreshed only at
the wide exchange would evolve locally between exchanges, and the locally evolved mirror image is not
bitwise the owner's value (the reflected stencil sums in the other order), so the owned fold-line row
would differ from the one-rank run from the second substep.  At `px > 1` the pair of fold exchanges per
substep (`bt_mid`: eta + ubt, `bt_late`: vbt with the fold-line projection) now runs on a second routing
context of `rdb_ocean_fold_exchange`, built from `fold_plan_build(ni, px, nghost + bt_halo, rx)`: the
same owner-routed mechanism (every value from the tile that owns the mirror point, ghost columns
included), `FOLD_CTX_NORMAL` for every state-level fold and `FOLD_CTX_WIDE` for the march-in.  The plan,
its device index lists and its buffers are per context (the pack / unpack kernels take them as
explicit-shape dummies).  Only the north rank row pays the extra fold exchanges; the seams keep the
`bt_halo/2` cadence.  Static wide fields beyond the fold are not folded: the outer `bt_halo` rows are clamped copies never read (ring poison), the first `nghost` rows are the model's own folded ghost rows and ARE read (corrected in stage 6); every dynamic row there is refolded each substep.  `bt_halo > 0` with the fold is refused when the north tile
is no taller than `nghost + bt_halo` rows (or, with `px > 1`, a tile no wider).

New cases `bth_trip_w{2,4,8}`, `bth_trip_w{4,8}_b0` (Murray cap above 74 N, double-drake land walls
reaching the fold line, meridional + zonal stress driving flow into the fold row) and
`bthx_tripsea_w{4,8}`, `bthx_tripsea_w4_b0` (a seamount in the cap): 48 x 32, nghost 3, `n_inner = 20`,
both split schemes.

| configuration | before | after |
|---|---|---|
| `bt_halo > 0` + tripolar fold, any layout | REFUSED at configure | accepted |
| refusal lifted only (local fold, `px = 1`: 1x1, 1x2, 1x4) | bitwise (the inline fold already ran over the wide band) | bitwise |
| refusal lifted only (`px > 1`: 2x1, 2x2, 4x1; plan built for `nghost`, asked for `nghost + W` rows) | FAIL: ~77 000 mismatching owned values, every field (`ml_h_layer` up to 2 m) | bitwise |
| wide loop with the fold switched off (sabotage, the non-vacuity control) | FAIL at 1x1, 2x1, 1x2: ~76 000 mismatching values | not applicable |

After: every `bth_trip_*` / `bthx_tripsea_*` case IDENTICAL at 1 rank (1x1), 2 ranks (2x1, 1x2) and 4
ranks (4x1, 2x2, 1x4), both split schemes, W = 2, 4, 8 where the smallest tile fits (`nghost + W + 1`;
the 1x4 layout at W = 8 is skipped), also with `RDB_BT_INHALO_POISON=2`, `RDB_BT_RING_POISON=1`,
`RDB_BT_RING_POISON=2` and both.


## Stage 6: the fold on the packed cadence (`feat/bt-wide-halo-all`)

Merge of `fix/bt-wide-halo-fold` (stage 4b: fold every substep, second owner-routed plan `FOLD_CTX_WIDE`) into the
packed static stage (stage 5).  gfortran 15.1 + OpenMPI 5.0.5 CPU, V100 where stated; 48 x 32, `n_inner = 40`, `nghost = 3`.

Interactions, and the evidence:

* **Reach across the fold.**  The fold runs every substep at the same two points as `bt_halo = 0` (eta after Pass 1,
  u after Pass 2b, v after Pass 2c), after the substep's periodic wrap, and the packed in-loop exchange only touches
  seams (the north fold edge has no Y neighbour), so the order is exchange -> wrap -> fold.  Every row beyond the
  fold is refreshed from owned mirror values every substep, so the fold edge consumes no reach.  Evidence: cadence + 1
  (`RDB_BT_INHALO_POISON=1002`) is bitwise over every `bth_trip_*` case at np 2 (30 identical) and np 4 (36), and
  cadence + 2 (`2002`) breaks on `bth_trip_w4` at np 2 and np 4 (NaN-catch, 286 / 572 hits) and on `bthx_tripsea`,
  `bthx_trip_{tides,porous,zfixed}` and `bth_trip_w2` at np 2 and 4 (run one by one), like every other family.  The
  in-loop poison does not poison the fold edge itself: those rows are expected valid.
* **Entry / exit sets, `H_ref` once.**  Packed entry and exit sets are exchanged at the seams only; the north rows
  beyond the fold are the model's own folded ghost rows (offset-copied; `H_ref` once) or are refolded.  The
  fold-line `v` update reads `H_ref` one row north of the fold line: with the north rows of the outer ring poisoned to
  depth `W` + 2 (`RDB_BT_RING_POISON=122`, new hundreds digit = north rows only) the owned cells stay bitwise at
  np 1, 2 and 4; at depth `W` + 3 = all `nghost` model rows (`132`) the run breaks (NaN-catch at every np), so the
  innermost model ghost row IS read, as at `bt_halo = 0`.  Rows 2 and 3 are not read in this configuration.  The exit
  set (time means, `*_end`, uhbt / vhbt ghost rows) is NOT compared directly against the reference (the comparison is
  owned cells); it is verified through the following steps (the pred_corr corrector reads them) and the write / resume
  leg (below), not row by row.
* **Counters.**  `oh_*` counters are unchanged by the fold (it is not a halo exchange): `pk = 2*steps*(2 +
  floor((n_inner-1)/(ng_w-2)))`, `bt_group = bt_u = 0`, `H_ref` centre-2D once, as for every `bth` case.  New: a
  per-context fold-group counter (`ocean_fold_exchange_calls`): on a rank that folds through the wide plan
  (north rank row, `px > 1`) exactly `2*steps*2*n_inner` groups (eta + ubt, then vbt with the projection, per substep,
  per solve), on every other rank 0, independent of the cadence.  Isend / Irecv counts inside a group depend on the plan's
  peers and are not asserted.
* **Ring poison, north edge.**  `RDB_BT_RING_POISON=2` (all edges) holds over the fold cases; see the `122` / `132` result above.
* **New cases.**  `bthx_trip_{tides,porous,zfixed}_w4` (fold + body tide, porous barriers, closed faces) are accepted and
  bitwise at np 2 and 4 (their per-solve `eta_forcing` / width exchanges are asserted by the counters); `bthx_tripsea_w4`
  also runs a checkpoint / resume leg with `bt_halo = 4` (2x1 at np 2; 4x1 and 2x2 at np 4), bitwise against the straight
  run, ghosts included.
* **Guard.**  The substep `error stop`s if the march-in meets a distributed fold without the wide plan.

Gate counts on the final code (identical = IDENTICAL lines; skips are tile-fit skips):

| leg | np 1 | np 2 | np 4 |
|---|---|---|---|
| `bth` plain | 118 | 356 | 432 (44 skipped) |
| `RDB_BT_INHALO_POISON=2` + `RDB_BT_RING_POISON=2` | 118 | 356 | 432 |
| cadence + 1 (`1002`) | - | 356 | 432 |
| whole unfiltered `test_ocean_decomp_bitid_mpi` | - | 838 | 1090 |

`test_ocean_dyn_mpi` np 1, 2 pass; ctest `-R mpi` (excluding the plain whole-suite legs): 59 / 59 (includes the btwide poison,
limit, break (WILL_FAIL), fold-ring legs); serial gfortran `ctest -R rdb -j 6`: 290 / 290; compat slice vs `origin/beta`: 50 PASS, 0 FAIL,
0 XPASS; `tools/check_closure_matrix.py` OK.  GPU (V100, GPUs 2 and 3, np 1 and 2): `bth` plain, poison and cadence + 1 all identical (np 1: 118 / 118 / 118;
np 2: 356 / 356 / 356), and the `tripolar` / `drain_ratio_fold` cases (np 1: 6, np 2: 18 identical).

Not verified: nvfortran CPU; GPU np 4; GPU cycle-break and fold-ring legs; GPU `test_ocean_dyn_mpi`; ice-on-fold or other
accepted combinations (sponge, Orlanski, nudged edges, wave drag) with the fold; timing of the merged result (stage 4 numbers
were measured without the fold); `bt_halo` odd values across the fold beyond the rounding rule.  Nothing new is refused: the
fold stays accepted with every tested combination; AUTO still resolves to 0.


## Stage 7: AUTO resolves to a wide halo (`feat/bt-wide-halo-auto`)

`&ocean_bt_nml bt_halo = -1` (the default) now resolves on more than one rank to a wide halo wherever a
bitwise `bth*` case proves it, and to 0 everywhere else.

### Measurement behind the width

V100 x 2 (GPUs 2 and 3, one rank each, HPC-X CUDA-aware, nvfortran 26.5), binary built from this branch, explicit
`bt_halo = W`.  SINGLE runs for the W scan (not medians), and the host was busy (1-min load 5 to 21: my own builds and
a CPU MPI test run, plus other projects) with the locks held and `nvidia-smi` clean on 2 and 3.  Read the
wall column with +-10 % in mind, and the `ocean_barotropic_solver` column (a profiler bucket, steadier) as the signal.
Benches: `seamount_bt8` (512 x 256 x 30, 2 d, `n_inner` = 65), `dg50` (300 x 300 x 50, 2 d, `n_inner` = 48), and the
same `dg50` with `dt_fixed` cut so that `n_inner` is small: `dgm` (dt 400, `n_inner` = 16, 0.67 d) and `dgs` (dt 150,
`n_inner` = 6, 0.25 d; these two runs are only 6 to 8 s long, setup is a visible part of them).

| bench (`n_inner`) | W | wall s | `ocean_barotropic_solver` s | `ocean_comms_bt` s |
|---|---|---|---|---|
| seamount (65) | 0 | 46.3 | 27.1 | 16.7 |
| | 2 | 34.9 | 15.6 | 4.7 |
| | 4 | 36.2 | 16.2 | 4.3 |
| | 8 | 37.5 | 16.9 | 3.3 |
| | 16 | 38.8 | 17.8 | 2.0 |
| dg50 (48) | 0 | 13.6 | 7.5 | 4.0 |
| | 2 | 8.5 | 3.1 | 1.0 |
| | 4 | 8.0 | 2.6 | 0.52 |
| | 8 | 7.9 | 2.5 | 0.41 |
| | 16 | 9.5 | 3.6 | 0.28 |
| dgm (16) | 0 | 6.9 | 1.68 | 0.96 |
| | 2 | 7.7 | 1.86 | 1.05 |
| | 4 | 7.6 | 1.73 | 0.91 |
| | 8 | 7.3 | 1.56 | 0.76 |
| | 16 | 7.2 | 1.46 | 0.66 |
| dgs (6) | 0 | 5.9 | 0.71 | 0.39 |
| | 2 | 5.8 | 0.59 | 0.16 |
| | 4 | 5.8 | 0.58 | 0.16 |
| | 8 | 5.8 | 0.56 | 0.14 |
| | 16 | 5.7 | 0.55 | 0.12 |

Reading.  The first wide step (W = 2) takes most of the gain wherever the solve has many substeps (the in-loop
exchange count drops from 3 + 1 rounds per substep to one packed one per cycle, and the mid-substep u exchange goes
away).  Past W = 4 the curve is flat to rising: the redundant band compute eats what the saved exchanges win (seamount:
+4 % solver per doubling; dg50 at 16: worse than 8).  W = 2 LOSES to `bt_halo = 0` on the mid case (`dgm`, 16 substeps,
cycle of 3 substeps), so a narrow cycle only pays when the solve is long or the tile forces it.  W = 8 is within noise
of the best on every bench, so the target is 8.  Wider than a solve is useless (the cycle `nghost + W - 2` past
`n_inner` buys no further exchange), so AUTO caps W at `n_inner - nghost + 2`; at `n_inner` = 6 that is 5, so 4 (the
rung below), which measured the same as 8.  `n_inner` below `BT_HALO_AUTO_MIN_NINNER = 4` resolves to 0: the cap leaves
no rung above 2 and nothing below `n_inner` = 6 was measured, so the guard is conservative, not measured.

The rule (`resolve_bt_halo`): one rank 0; excluded 0; else the widest rung of {8, 4, 2} not above the `n_inner` cap
that fits the smallest tile (`nghost + W <= min(nx/px, ny/py)`, one more across the tripolar fold), shrinking and never
refusing; 0 if none fits.  AUTO is resolved after `configure_ocean_bt_split` because `auto_n_inner` derives `n_inner`
there; every input is rank-uniform.  Rungs are only widths a `bth*` case proves (2, 4, 8): W = 16 has no case and is not
selected.

### AUTO vs `bt_halo = 0` on the same binary, 2 GPUs (two runs each, interleaved, load 10 to 20)

| bench | resolved | `bt_halo = 0` wall s | AUTO wall s | solver s (0 / AUTO) |
|---|---|---|---|---|
| seamount (`n_inner` 65) | 8 | 46.9, 47.2 | 32.0, 32.9 | 27.6 / 12.7, 13.4 |
| dg50 (48) | 8 | 10.1, 10.2 | 7.9, 7.9 | 4.7 / 2.5 |
| dgm (16) | 8 (cap 15 -> 8) | 6.9, 7.1 | 6.2, 6.3 | 1.7 / 1.0 |
| dgs (6) | 4 | 5.9, 5.9 | 5.8, 5.8 | 0.70 / 0.57 |

So AUTO buys 30 % of wall on the long-solve seamount bench, 22 % on `dg50`, 10 % on `dgm`, 2 % on the 6-substep case
where the run is short and dominated by everything but the barotropic solve.  On a realistic case with a 3D
mixing-heavy step (extreme, `nk100`: ML-comm bound) the stage-5 measurements stand: -5 % to ~0 %.  All of this is np 2;
np 1 is unchanged by construction (AUTO = 0 there).

### Exclusions (AUTO = 0), by choice (b) unless stated

Single list `bt_halo_auto_exclusion` (`src/core/config/rdb_config_validate.F90`); an explicit `bt_halo` is honoured.
The first five also refuse an explicit `bt_halo > 0` as before.

| feature | choice | why |
|---|---|---|
| wet/dry, `use_cont_type`, `upstream_h_face`, Chapman edge, cavity dynamics | refused today (a) n/a | not widened |
| `inhalo_march`, `dtbt_reset_period >= 0` | refused today | exclusive / `n_inner` re-derived mid-run |
| sea ice | (b) | no `bth*` case; ice runs its own seam exchanges and the EVP / transport halo reach |
| edges clamped, sponge, nested, inflow, discharge | (b) | only wall, open (Flather), tidal, periodic and the fold have cases (`bthx_perxc` pairs a clamped south edge with periodic x only) |
| Orlanski radiation, open-edge nudging, tracer reservoirs | (b) | the wide band evolves under array-edge closures the OBC target does not see |
| wave drag | (b) | multiplies `bt_rem` inside the substep; no case |
| map-driven sponge (`&ocean_sponge_nml`) | (b) | no case |
| ice-shelf cavity melt | (b) | cavity dyn is already refused; melt without it is not exercised |

Each is a candidate for (a) later: add a `bth`/`bthx` case at np 2 and 4, then drop its line from the exclusion list;
the unset-`bt_halo` cases of the test already assert that the excluded ones resolve to 0.

### Test gate

Every case of `test_ocean_decomp_bitid_mpi` that leaves `bt_halo` unset runs AUTO on its decomposed legs against the
one-rank (AUTO = 0) reference.  The test records the resolved width per layout (`case ...: AUTO bt_halo = N`), checks
the packed-exchange counter against it, asserts that sea ice / wet/dry / cavity / dtbt / wave drag / chapman / sponge
cases resolved to 0 and that island, spherical, tripolar, periodic channel, closures and eulz_periodic resolved wide
wherever the tile holds `nghost + 2` (+1 under the fold), prints the totals and fails below 10 wide cases.  Restart:
`island_basin`, `spherical` and `tripolar` write/resume under AUTO on both legs; `island_basin` and `tripolar` also write
at `bt_halo = 0` and resume under AUTO and the reverse (owned cells against the one-rank reference: the ghosts of the
two marches legitimately differ).  The non-trivial widths this exercises include 6 / 4 / 2 shrunk by the odd 26 x 18 tiles.
