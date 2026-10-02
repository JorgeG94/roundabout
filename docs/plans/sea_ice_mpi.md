# Plan: MPI-aware sea ice

Status: design only, no code changed. Target: v0.1.0 global 0.25° tripolar
(1440×1080×75) with sea ice on several GPUs, comparable to ACCESS-OM3.
The contract is CLAUDE.md §MPI: **a decomposed run is the serial run, bitwise**.

Approach: reuse what already exists. Exchanges go through the ocean halo layer
(`src/comm/rdb_ocean_halo.F90`). Every sum goes through the existing EFP path.
Restarts use the restart registry, and tests use the `decomp_bitid` harness and
`test_ocean_restart_engine`. No new comm framework is needed. The ice uses the
ocean's decomposition and the ocean's `nghost`, and exchanges once per EVP
subcycle.

Citations are `path:line` against `origin/main` at `11b2d134d`.

---

## 1. Current state

### Refusal sites (all must be lifted or narrowed)

| Site | What it refuses |
|---|---|
| `src/core/ocean/rdb_ocean_engine.F90:405-411` | `&ocean_ice_nml enable` on >1 actual rank (engine gate, keyed on `csize`) |
| `src/core/rdb_config.F90:6653-6660` | `enable` with `px*py > 1` |
| `src/core/rdb_config.F90:6714-6720` | `transport` with `px*py > 1` |
| `src/core/rdb_config.F90:6742-6748` | `dynamics` with `px*py > 1` |
| `src/core/rdb_config.F90:6700-6710` | `transport` with ANY periodic edge ("ice ghost cells are never wrapped"). **Blocks every global run, even on one rank.** |
| `src/core/rdb_config.F90:6767-6772` | `dynamics` with `north='tripolar_fold'` |
| `tests/mpi/test_ocean_decomp_bitid_mpi.F90:755` | fence test asserting the refusal |

### Wiring into the ocean step

- Callers run `engine_step` → `engine_step_ice` → `engine_step_finalize`
  (`src/driver/rdb_driver.F90:361-373`, `src/api/rdb_ocean_api.F90:594`).
- `engine_step_ice` (`rdb_ocean_engine.F90:1364-1545`) does the following:
  - frazil accumulate on every step (`:1420`);
  - EVP on every outer step, followed by the ice→ocean stress blend (`:1435-1457`);
  - at thermo cadence (`:1499-1500`): transport (`:1506-1519`), then the
    forcing → basal → frazil uptake → column → snowfall → brine/heat/SW
    couplers → ITD chain (`:1521-1543`).
- Setup:
  - The slot is allocated at `rdb_ocean_state.F90:487`, with its flags at `:579-596`.
  - EVP params and the `tau_a` snapshot are taken at `rdb_ocean_engine.F90:766-779`.
  - The `Q_salt`/`Q_heat` resume fold is at `:648-661`.
  - The IC is at `:1016-1033`. It is already skipped on a warm restart (`:1017`).
  - Device mapping: `rdb_ocean_state.F90:891` and `:926`.
- The ice uses the ocean's tiles: arrays are sized from the same `grid%nx_total/ny_total`
  (`src/core/ice/state/rdb_ice_state.F90:686-687`).

### Ghost policy today

The ice code has no exchange anywhere.
- EVP wraps its own copies only on periodic axes (`rdb_ice_evp.F90:513-518, 545-550, 568-571, 631-634`).
- Every ghost outside the physical domain is pinned to land (`evp_build_masks_impl`, `rdb_ice_evp.F90:726-734`).
- Transport zeroes all CAS ghosts (`rdb_ice_transport.F90:289-298`).
- The IC seeds ghosts by formula once (`rdb_ice_init.F90:302-308`). After that, nothing refreshes the ice ghost state.

---

## 2. Kernel inventory

Ice grid: **C-grid**. `u` sits on the WEST face of cell `(i,j)`, `v` on the
SOUTH face, and corners are SW. Stress `str_d`/`str_t` sit at T-cells and
`str_s` at corners (`rdb_ice_evp.F90:12-17`). Below, "phys" means the loop
covers physical cells/faces only and "all" means the full local array, ghosts
included.

| Kernel (file:line) | Computes | Stencil / ghost reads | Range | Per subcycle? |
|---|---|---|---|---|
| `ice_cell_concentration_impl` `rdb_ice_state.F90:1324-1382` | gather `mis/mice/ci` from categories | point | all (reads IST ghosts) | no (per call) |
| `evp_build_masks_impl` `rdb_ice_evp.F90:701-761` | `mask_t/u/v/q` | ±1; zeroes non-periodic ghosts | all | no |
| `evp_pres_mice_impl` `:825-846`, `evp_mi_face_impl` `:1020-1053` | strength, face mass | point / ±1 | all | no |
| `evp_q_and_mi_ratio_impl` `:1055-1114` | corner `q`, `mi_ratio_A_q` | 4 cells around corner | all | no |
| `ice_limit_stresses` `:1120-1177` | clamp stresses | point / 4 cells | all | no (once per call) |
| `evp_sh_ds_impl` `:1224-1258` | shear strain at corners | ±1 face | all | yes |
| `evp_sh_dd_dt_impl` `:1260-1277` | div/tension at cells | faces i,i+1 / j,j+1 | all | yes |
| `evp_zeta_impl` `:1279-1304` | viscosity | 4 corners | all | yes |
| `evp_stress_relax_impl` `:1306-1322`, `evp_str_s_relax_impl` `:1324-1363` | elastic relax of `str_d/t`, `str_s` | point / 4 cells | all | yes |
| `evp_u_momentum_impl` `:1401-1553`, `evp_v_momentum_impl` `:1555-1681` | `ui/vi`, accumulate `fxoc/fyoc` | stresses ±1, `vi`/`vo`/`mi_v` at i-1,j+1 | **phys faces incl. both edge faces** (`:1457-1462`, `:1591-1594`) | yes |
| `evp_truncate_*` `:875-1018` | CFL clip (+ count) | donor `areaT` ±1 | phys | optional / final |
| `evp_average_stress_impl` `:1683-1701` | `fxoc/fyoc` mean | point | all | no |
| transport `ice_cat_flux_{x,y}_impl` `rdb_ice_transport.F90:421-570, 572-693` | PPM on summed mass, proportionate category split | **±3 cells per face** (5-pt PPM on each donor) | all rows/cols | per adv substep |
| `ice_gather_flux_*` / `ice_ride_update_*` `:717-1045` | PCM riding of `m_ice`, `enth_ice`, `sal_ice`, `enth_snow` | donor ±1 | all | per adv substep |
| `ice_mass_update_{x,y}_impl` `:1047-1081` | `mca` update | ±1 face | x: phys i, all j; y: phys j, all i | per adv substep |
| CAS↔IST, compress `:269-299, 1127-1265` | category state | point | phys | per transport call |
| ITD `ice_adjust_categories_impl` `rdb_ice_itd.F90:104-140` | re-bin categories | column | phys | thermo |
| column thermo `ice_thermo_columns` `rdb_ice_column.F90:812-933`; driver fills/reduces `rdb_ice_thermo_driver.F90:147-401` | Winton column, per-cell coupling diags | column | phys | thermo |
| atm forcing `rdb_ice_atm_forcing.F90:36-51` | restoring seam | point | all | thermo |
| frazil `rdb_ice_frazil.F90:76-130`; uptake `rdb_ice_frazil_uptake.F90:139-336`; basal `rdb_ice_basal_flux.F90:69-135`; snowfall `rdb_ice_snow.F90:99-159` | ocean-surface sampling, frazil bank, basal flux | column (reads ocean `k=nz` T/S) | phys | every step / thermo |
| brine/heat/SW couplers `rdb_ice_ocean_coupler.F90:93-288` | write `Q_salt/Q_heat/q_sw` | point | all | thermo |
| stress blend `rdb_ice_ocean_coupler.F90:318-434` | `tau=(1-a)tau_a+a·fxoc`, `a_u=½(ci(i-1)+ci(i))` | ±1 cell **incl. ghost `ci`** | all faces | every step |

The thermodynamic part of the ice (column, ITD, frazil, basal, snowfall,
couplers) is column-local and needs no halo. Only three parts reach across
cells: EVP, transport, and the stress blend.

---

## 3. Halo needs → minimal exchange set

All exchanges call the existing primitives. These are `ocean_halo_centre`
(2D/3D), `ocean_halo_face_x`, and `ocean_halo_face_y` (`rdb_ocean_halo.F90:143-152`).
Each one does a two-pass X-then-Y exchange, so corner ghosts come out valid
(`:5-15`). Each one also follows the D1 rule: the west/south rank owns the
duplicated seam face (`:17-39`). Each call moves device buffers with CUDA-aware
MPI (`device_resident` defaults to true, `:545-552`) and sends one message per
direction for all levels (`:561-568`). The same order as the ocean's seam sites
applies (`rdb_ocean_halo_state.F90:101-106, 135-159`): exchange, then the local
periodic wrap only on an axis the halo does not own
(`ocean_halo_is_decomposed_x/_y`), then the fold.

| # | When | Fields | Primitive | Replaces |
|---|---|---|---|---|
| X1 | end of the thermo block (after ITD, `rdb_ocean_engine.F90:1542`), and once at cold-start setup (host, `device_resident=.false.`, **skipped on warm restart**) | IST: `part_size` (nz=ncat+1), `m_ice`, `m_snow` (ncat), `enth_ice`, `sal_ice` (nz=ncat·nk), `enth_snow` (ncat) | `ocean_halo_centre` 3D, 6 calls. 4-D arrays pass as contiguous `nz=ncat*nk`. | nothing (ghosts are stale today) |
| X2 | each EVP subcycle, top of loop | `ui`, `vi` | `face_x` 2D + `face_y` 2D | wraps `rdb_ice_evp.F90:568-571` |
| X3 | after the final CFL clip | `ui`, `vi` | same | wraps `:631-634` |
| X4 | top of each transport adv substep | `mca_ice`, `mca_snow`, `m_ice`, `enth_ice`, `sal_ice`, `enth_snow` | `ocean_halo_centre` 3D | zeroing at `rdb_ice_transport.F90:289-292` |
| X5 | after the stress blend | `tau_x`, `tau_y` (+ `stress_mag` re-derived) | existing `ocean_seam_refresh_surface_stress` (`rdb_ocean_halo_state.F90:81-163`) | `ocean_surface_stress_refresh_mag` at `rdb_ice_ocean_coupler.F90:388` |

Why this set is enough:

- **EVP inputs.** X1 makes the IST ghosts valid, so `ice_cell_concentration_impl`
  (which runs over the full array) yields valid ghost `mis/mice/ci`. `mask_t` must
  then follow `wet_T` in seam ghosts instead of being zeroed. Gate the zeroing at
  `rdb_ice_evp.F90:726-734` on `has_*` (non-periodic physical edge only), which is
  CLAUDE.md's "gate edge/ghost writers on `has_*`" rule. `ice_evp_step` therefore
  needs `bc` and not just `periodic_x/y` (`rdb_ice_evp.F90:288-289`).
- **Stresses need no exchange.** Every subcycle kernel except momentum runs over
  the full array, and the stress updates are point-local in the stresses. Given
  exchanged `ui/vi`, the ghost stresses are recomputed locally. Only the
  outermost ring is wrong, and with `nghost=3` no physical face reads it:
  momentum at face `ng+1` reads cells `ng, ng+1` and corner `ng+1`. The serial
  code already relies on this: it wraps the stresses once per call (`:545-550`)
  and never inside the loop. As a result, no corner-field exchange primitive is
  needed. (`evp_wrap_corner_impl`, `:1183-1222`, stays for non-decomposed
  periodic axes.)
- **Transport.** The x-pass flux and mass update run on all `j`
  (`rdb_ice_transport.F90:519, 1060`), so after X4 (which has valid corners) the
  x-pass refreshes the ghost rows that the y-pass reads. One exchange per substep
  is therefore enough at `nghost=3`. After transport, CAS→IST, compress and ITD
  are physical-only, and X1 refreshes their ghosts.
- **Coupler.** `a_u` at the westmost owned face reads ghost `ci`
  (`rdb_ice_ocean_coupler.F90:414-422`), which X1 makes valid. `fxoc` ghost faces
  are 0, because the momentum loop writes physical faces only. X5 overwrites the
  ghost `tau` with the owner's value, which the ocean's KPP/EPBL/MLE ghost reads
  need (`rdb_ocean_halo_state.F90:85-99`).

### Required kernel fixes, not exchanges

- **F1 — transport wall zeroing.** `rdb_ice_transport.F90:545-548` and `:665-668`
  zero faces `ng+1` and `ng+nphys+1` unconditionally, which turns every seam into
  a wall. Gate them on `has_west/east/south/north` and on non-periodic edges.
- **F2 — PPM edge fallback overwrites a value the seam face uses.** The fallback
  at `:490-502` overwrites `hr_x_work(3,:)` (y twin: `hr_y_work(:,3)`, `:621-633`)
  after the PPM loop (`:470-487`) wrote it. With `nghost=3`, face 4 is the
  west/south seam face, and its `u>0` flux reads `hr_x_work(3)` (`:522-525`).
  The result: a first-order value at the seam, PPM in serial, so the run is not
  bit-identical. Fix: do not overwrite `hr_*_work(3)`. This is bit-identical
  serially, because face `ng+1` is a zeroed wall on one rank.
- **F3 — local wraps** in EVP must skip decomposed axes. This is the same
  `skip_x/skip_y` idiom as `rdb_ocean_dyn.F90:3627-3629`.

**Uncertain:** whether the ocean surface-layer `u/v` ghosts (`uo/vo`, which EVP
reads at `i-1, j+1`, `rdb_ice_evp.F90:1482-1487`) are fresh after `engine_step`.
The last ML exchange is at the stage end (`rdb_ocean_dyn.F90:5082`), but ALE remap
and finalisation come after it and I did not trace them. PR 1's bit-id case
settles this. If the ghosts turn out stale, exchange the contiguous slice
`u_face_x_layer(:,:,nz)` / `v_face_y_layer(:,:,nz)` in place before the EVP call:
two more exchanges per outer step.

---

## 4. Halo strategy: one exchange per EVP subcycle

- **Knob.** `&ocean_ice_nml evp_sub_steps`, default **432** (`rdb_config.F90:377`).
  EVP runs on **every outer step** (`rdb_ocean_engine.F90:1427-1434`), not at
  thermo cadence. Default `nghost = 3` (`rdb_config.F90:3188`).
- **Cost of X2.** There are 2 calls per subcycle. Each is two-pass with up to 4
  neighbour messages, so 8 isend+8 irecv and 4 waitall per subcycle, which is
  ≈3456 messages and 1728 waits per outer step.
  - Example tile: 1440×1080 on 4×4 ranks gives 360×270 cells, `ng=3`. Each
    message is ≈(ng+1)·276 ≈ 1.1k doubles (E/W) or ≈(ng+1)·367 ≈ 1.5k doubles
    (N/S), so ≈80 KB per subcycle and ≈35 MB per outer step per rank.
  - The cost is latency-bound. At an assumed ~10 µs per wait, that is ≈17 ms per
    outer step. Compare ≈10 kernel launches per subcycle × 432 ≈ 4300 launches
    per step, which the GPU pays regardless. These numbers are estimates and
    have not been measured.
- **Wide halo** (exchange every *k* subcycles). Each subcycle uses about 2 rings
  (u → strain → zeta → `str_s` → momentum), so this needs ≥2k ghost rings. It
  would also need the momentum loop (`rdb_ice_evp.F90:1457-1462`) extended into
  the halo, and a wider array set (the BT-clone pattern; primitives exist at
  `rdb_ocean_halo.F90:1673-1760`).
- **Recommendation.** Do the per-subcycle exchange on the ocean's `nghost`. It
  needs no array reshaping and is bit-identical by construction. Consider a wide
  halo only if profiling the 0.25° case shows the X2 waits are a significant
  share of step time. It is not part of this plan.

---

## 5. Global reductions in the ice path

| Reduction (file:line) | Today | Steers? | Required |
|---|---|---|---|
| transport `vmax` max-speed, early return on `vmax == 0` (`rdb_ice_transport.F90:183-185`, kernel `:240-262`) | rank-local | **yes**. One rank skips the CAS↔IST round trip (which perturbs last bits, `:180-182`), another does not, and the next X4 deadlocks. | `halo_allreduce_max` (`src/comm/rdb_halo.F90:798`). A max is exact, so no EFP is needed. |
| transport validity `ok` (`:359-360, 413-414`, kernel `:1087-1121`) and compress `ok` (`:206-209`, kernel `:1256`) | rank-local | **yes** (early `return` + driver `error stop`, `rdb_ocean_engine.F90:1514-1519`). A rank-local abort strands the others in X4. | Allreduce min of a 0/1 flag (`halo_allreduce_min`, `rdb_halo.F90:777`) **before** any further exchange. All ranks then abort together. |
| EVP `n_trunc` count (`rdb_ice_evp.F90:985, 1003`) | rank-local; only rank 0 logs its own count (`rdb_ocean_engine.F90:1448`) | no (warning only) | Exact sum (`halo_allreduce_sum` of an integer-valued real). Only when `cfl_trunc > 0`. |
| console ice area / conc / thickness (`rdb_ocean_console_stats.F90:453-456, 476-480` EFP; `:510-535` FP) | **already global + EFP** | no | nothing |
| frazil heat budget (`rdb_ocean_console_stats.F90:664-673` EFP; `:690-694` FP) | **already global + EFP** | no | nothing |
| diag fills `ice_conc`/`ice_thick` (`rdb_ocean_diag_fills.F90:337-357`), `ice_speed/u/v` (`rdb_ocean_diag_derived.F90:873-920`) | per-cell, per-rank files | no | nothing |

No other reductions exist: ITD, column, frazil and couplers contain no
reductions across cells. Every sum is already on EFP. The two control-steering
reductions are a max and a logical, both of which are exact, so they only need
to be **rank-uniform**.

---

## 6. Restart

### Registered today

All entries are optional, with full local arrays and ghosts (`ng=0` in the
registry). They are registered at `rdb_ocean_state.F90:2095-2302`:

- `ice_frazil_heat`;
- `ice_m_frozen_diag`, `ice_salt_flux_diag`, `ice_m_melt_diag`, `ice_heat_flux_diag`, `ice_sw_thru_diag`;
- `ice_u_ice`, `ice_v_ice`, `ice_str_d`, `ice_str_t`, `ice_str_s`, `ice_fxoc`, `ice_fyoc`;
- `ice_tau_ocn_x/y` + scalar `ice_tau_ocn_valid`;
- `ice_part_size`, `ice_m_ice`, `ice_m_snow`;
- `ice_enth_ice_k*`, `ice_sal_ice_k*`, `ice_enth_snow_k*` (per-k slices).

### Rules from `67f559c78` / `8e1931f20`, applied to ice

1. **Every piece of state carried across steps is registered.** I checked this.
   - The EVP workspace is rebuilt at the start of every call (`rdb_ice_evp.F90:508-559`, before the subcycle loop at `:561`).
   - The transport workspace is rebuilt from IST on every call (`rdb_ice_transport.F90:189-190`).
   - The coupler scratch is rebuilt every call (`rdb_ice_ocean_coupler.F90:359-364`).
   - The `atm_*`/`fb`/seam fields are refilled every thermo step (`rdb_ice_state.F90:484-508`).
   - `tau_a_x/y` is a deliberate configure-time snapshot (`rdb_ocean_engine.F90:776-777`).
   - `heat_budget_frazil` (unregistered) is a console-only accumulator. The ocean's console drift restarts from the resume point (`rdb_driver.F90:283-290`).
   - **Conclusion: no new carried scratch is needed.**
2. **Persistent arrays attach by `copyin`.** All ice persistent fields already do
   (`rdb_ice_state.F90:899-914`), and only scratch uses `create` (`:916-923, 1166-1174`).
3. **No setup pass may re-derive checkpointed ghosts.**
   - The IC is already skipped (`rdb_ocean_engine.F90:1017`).
   - The new cold-start X1 must sit under `.not. did_restart`, like the ocean's wrap and exchange (`rdb_ocean_engine.F90:848-853, 909-913`, from `8e1931f20`).
   - It must not be added to `engine_enter_data`'s warm-up, which is skipped on warm restart (`:1243-1246`).
   - `ice_ocean_stress_resume_apply` copies full arrays, ghosts included (`rdb_ice_ocean_coupler.F90:458-501`), so it is already compliant.
4. **Per-rank files round-trip bitwise.** No new mechanism is needed: X1–X5 run
   inside the step, so the straight run and the resumed run perform the same
   exchanges.

---

## 7. Tripolar fold

The ice does not apply the fold anywhere today: `rdb_config.F90:6767-6772`
refuses it, and `rdb_ice_evp.F90:93` notes "D6". The fold primitives exist in
`src/core/ocean/boundary/rdb_ocean_fold.F90:104-107`: centre, u-face, v-face,
and corner with a `negate` flag (`:305-312`). The distributed version
(`px > 1`) is the sibling plan `docs/plans/tripolar_fold_px_gt_1.md`. The ice
needs it for the field classes below.

| Class | Ice fields | Fold call |
|---|---|---|
| scalar, centre | IST (`part_size`, `m_ice`, `m_snow`, `enth_*`, `sal_ice`), CAS (`mca_*`), `frazil_heat` | `fold_north_centre` |
| vector, faces | `u_ice`/`v_ice`, `fxoc`/`fyoc`, `tau_a_*`, `tau_ocn_*` | `fold_north_u_face` / `fold_north_v_face` (sign flip + duplicated-row projection) |
| 2nd-rank tensor components (invariant under the fold's 180° rotation: two sign flips) | `str_d`, `str_t` (centre), `str_s` (corner) | `fold_north_centre`; `fold_north_corner(negate=.false.)` |

Where the fold goes in the ice step: after X1, X2, X3 and X4, each time after
the exchange and wrap; inside X5 (already folded, `rdb_ocean_halo_state.F90:155-159`);
and **after the stress update in every subcycle**. Across the fold, recomputing
the stresses in the mirrored orientation is not bit-exact, so the redundant-compute
argument from §3 does not carry over. That adds 3 rank-local fold kernels per
subcycle when `px = 1`.

**Uncertain, needs its own test:** the transport's y-flux through the duplicated
fold `v`-row. Mass conservation across the fold has not been analysed here.

---

## 8. GPU device residency

- **Device-resident:** all persistent ice arrays (`copyin`, `rdb_ice_state.F90:899-914`),
  the transport workspace (`:917-923`), the EVP workspace (`:916`, `:1166-1174`),
  and the coupler scratch (`rdb_ice_ocean_coupler.F90:297-300`).
  `ocean_sea_ice_t%enter_data` rides `ocean_state_enter_data` (`rdb_ocean_state.F90:891`).
- **Host-only by design:** `ice_init_apply` (`rdb_ice_init.F90:265-269`, before
  mapping), `ice_ocean_stress_resume_apply` (`rdb_ice_ocean_coupler.F90:478-490`,
  before mapping), and `tau_ocn_valid` (never mapped, `rdb_ice_state.F90:244-255`).
- **Per-step host syncs:** these are only the reduction scalars. That is
  `n_trunc` once per step, `vmax` once per transport call, and the two `ok`
  checks per pass plus compress. The allreduces in §5 add no extra D→H, because
  these scalars already land on the host.
- **New exchanges:** they use `device_resident=.true.` on arrays that are already
  mapped. The cold-start X1 uses `.false.` and runs before `enter_data`. Per
  CLAUDE.md, update and map component arrays only, never `!$acc update` the
  aggregate `ice`. `part_size` (lower bound 0) and the 4-D `enth_ice`/`sal_ice`
  are contiguous, so they can pass to the explicit-shape halo dummies by
  sequence association.

---

## 9. Tests

1. **Decomposition bit-identity** (`tests/mpi/test_ocean_decomp_bitid_mpi.F90`).
   - Add an ice case to `CASES` (`:120-124`) and `case_nml` (`:164-318`).
     Suggested setup:
     - 26×18, `nghost=3`, re-entrant channel (`west/east='periodic'`, walls N/S);
     - `&ocean_ice_nml enable, ncat=5, transport, dynamics, evp_sub_steps=30`;
     - `&ocean_ice_ic_nml conc_config='uniform'`, wind, thermo on, `dt_therm_ratio = 2` so both cadences run.
     - Ice must actually cross seams, and the periodic edge exercises F1/F3 and X5.
   - The harness already runs serial, then 2×1, 1×2 at np=2 and 4×1, 1×4, 2×2 at np=4 (`:11-23`).
   - Required harness edits:
     - `run_one` must call `engine_step_ice` between step and finalize (`:346-349`);
     - raise `MAXF` (`:94`), because the ice adds about 23 registry fields;
     - drop the ice entry from `check_single_rank_fences` (`:755`) or keep it as a tripolar-only fence.
   - Gate it on the **GPU** build. Commit `26cdf275a` found that the CPU
     vectoriser's body/remainder FMA split is position-dependent, and 432
     iterated EVP subcycles could amplify that. If a CPU build fails only in EVP
     fields, treat it as that issue, not as a halo bug: check `-Mnofma` on the
     EVP module before exempting anything.
2. **Engine-path restart** (`tests/test_ocean_restart_engine.F90`).
   - Add an ice variant: `advance` (`:127-142`) calls `engine_step_ice`.
   - `take_snapshot` (`:144-167`) also compares every ice registry field, full extent and ghosts included.
   - Use `N_WRITE` odd with `dt_therm_ratio=2`, so the checkpoint lands mid-thermo-window.
   - For per-rank files: give the bit-id ice case a write/resume leg (checkpoint at step 24, resume, compare full arrays against the straight decomposed run). This reuses the harness without adding a new MPI test.
3. **Unit tests.** No new halo unit test is needed, since the primitives are
   already covered by `tests/mpi/test_halo_ocean_mpi.F90`. Add one serial
   regression for F2: a transport flux at face `ng+1` with `u>0` and a curved
   profile should equal the interior PPM flux.

---

## 10. Open decisions (real choices only)

1. **Serial periodic answer change.** Today, on one rank, the coupler reads
   stale IST ghosts across a periodic seam (`rdb_ice_evp.F90:130-137`, the
   documented "PRE-EXISTING" D7 note). X1 and X5 fix that, which changes serial
   answers for existing periodic + ice-dynamics configs. *Options:* accept the
   change, or keep the stale path behind a knob. *Recommend:* accept it and note
   it in the PR. Bit-identity across decompositions cannot hold otherwise.
2. **F2: fix the fallback, or require `nghost ≥ 4` for multi-rank transport.**
   *Recommend:* fix it. It is serially bit-identical and keeps one halo width for
   ocean and ice.
3. **Where X1 sits.** Options: the end of the thermo block, once per window, or
   the top of each consumer (EVP every step, coupler, transport). *Recommend:*
   the end of the thermo block. The IST changes only there, so this is the
   fewest exchanges.

These are settled rather than open: the ice shares the ocean decomposition
(same arrays, same `ocean_halo_init` at `rdb_ocean_engine.F90:894`, and column
coupling is point-local); the halo width is the ocean `nghost`; and each IST
array goes out as one `centre_3d` call carrying all its categories, so there is
no cross-array packing.

---

## 11. Slicing

| PR | Scope | Files | Tests | Still refused after | Size |
|---|---|---|---|---|---|
| **1. Dynamics + thermo + coupler on >1 rank** | X1, X2, X3, X5; F3; `has_*`-gated `mask_t`; `bc` into `ice_evp_step`; `n_trunc` allreduce; lift `rdb_ocean_engine.F90:405-411`, `rdb_config.F90:6653-6660, 6742-6748`; add a transport-only multi-rank fence | `rdb_ocean_engine.F90`, `rdb_ice_evp.F90`, `rdb_ice_ocean_coupler.F90`, `rdb_config.F90`, `docs/CAPABILITIES_AND_LIMITATIONS.md` | bit-id ice case with `transport=.false.`; resolve the §3 `uo/vo` question here | transport on >1 rank and on periodic grids; tripolar + dynamics | L |
| **2. Transport on >1 rank and periodic** | X4; F1; F2; rank-uniform `vmax` / `ok`; lift `rdb_config.F90:6700-6710` (periodic) and `:6714-6720` (px·py) | `rdb_ice_transport.F90`, `rdb_ice_itd.F90` (interface only, if needed), `rdb_config.F90` | bit-id case with `transport=.true.`; F2 regression | tripolar + transport | M (merge into PR 1 if it stays small) |
| **3. Restart round trip** | cold-start-only X1 at setup; ice variant of the engine restart test; write/resume leg in the bit-id ice case | `rdb_ocean_engine.F90`, `tests/test_ocean_restart_engine.F90`, `tests/mpi/test_ocean_decomp_bitid_mpi.F90` | listed in §9.2 | tripolar | S |
| **4. Ice on the tripolar fold** (after the fold plan lands) | fold calls from §7 after X1–X4 and after the subcycle stress update; lift `rdb_config.F90:6767-6772`; transport across the fold row | `rdb_ice_evp.F90`, `rdb_ice_transport.F90`, `rdb_ocean_engine.F90`, `rdb_config.F90` | tripolar ice case in the bit-id harness (`px=1` first, then `px>1` once the fold plan lands); conservation across the fold | — | M |

Each PR updates `docs/CAPABILITIES_AND_LIMITATIONS.md` (the MPI single-rank
list) and the CLAUDE.md "Single-rank features" sentence in the same change.
