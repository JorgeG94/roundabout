# Distributed tripolar north fold (`px > 1`): design plan

Status: plan only, no code. Target: v0.1.0 global 0.25° tripolar (OM4_025-sized,
1440 × 1080 × 75) on many GPUs, **bitwise equal to the serial run** (the
`test_ocean_decomp_bitid_mpi` contract). Citations are `path:line` against
`origin/main` at the time of writing.

## 0. Decisions in one screen

- **Owner-routed exchange:** each north-ghost / fold-row value comes from the rank
  that OWNS the mirror point (T/v: the cell's tile; u/corner: the D1 face owner),
  never from a halo copy, so the fold needs no preceding X pass (§2.3, §5).
- **Sign at unpack**; **projection stays a copy** (east half authoritative), so no
  floating-point arithmetic has to agree across ranks (§3, §6).
- **px = 1 keeps the local kernels**; px > 1 always uses the exchange, the
  self-mirror tile via a local copy. Separate fold call, not a third halo pass (v1).
- **Barotropic loop:** fold at the serial program points, 2-3 per substep (§5.3).

## 1. Current state

### 1.1 Staggering and index maps

Storage index = `nghost + physical` (`rdb_ocean_fold.F90:26-27`). T at centre, u on
the WEST face, v on the SOUTH face, corner at SW (`rdb_ocean_fold.F90:28-31`). The
fold line is `y = nj`; a point `(x, y)` north of it maps to `(ni - x, 2nj - y)`
through a 180° rotation: true vector components negate, scalars and vorticity copy
(`rdb_ocean_fold.F90:38-43`). Storage maps (`rdb_ocean_fold.F90:46-51`):

| Stagger | i-map (storage)  | j-map (storage)   | fold-line row | kernel |
|---------|------------------|-------------------|---------------|--------|
| T       | `2ng+ni+1 - i`   | `2ng+2nj+1 - j`   | none          | `fold_north_centre_{2d,3d}` `:134-168` |
| u (Cu)  | `2ng+ni+2 - i`   | `2ng+2nj+1 - j`   | none          | `fold_north_u_face_*` `:177-211`, negates `:191,:209` |
| v (Cv)  | `2ng+ni+1 - i`   | `2ng+2nj+2 - j`   | `ng+nj+1`     | `fold_north_v_face_*` `:224-295` |
| corner  | `2ng+ni+2 - i`   | `2ng+2nj+2 - j`   | `ng+nj+1`     | `fold_north_corner_2d` `:305-345`, `negate` arg |

T/u rows `j ≥ ng+nj+1` are pure halo fills (`rdb_ocean_fold.F90:53-54,145-148`).
v/corner have a row ON the fold line, storage `ng+nj+1` = north face of the last
T-row (`rdb_ocean_fold.F90:55-62`); halo fill starts at `j_fold+1`
(`:237,:241-243`). In every kernel `nx_phys`/`ny_phys` are the TILE's extents
(`:143-145`), which is why the kernels are exact only when the tile holds the
whole row (`rdb_ocean_fold.F90:10-22`).

### 1.2 Duplicated fold-line DOF

On row `ng+nj+1`, slots `i` and `i'` are the same physical face/vertex with opposite
orientation (`rdb_ocean_fold.F90:75-84`). The projection, over every storage column
including periodic ghosts, computes `p = modulo(i-i_lo, ni)+1`, `pm = ni+1-p` (v)
or `pm = modulo(ni+1-p, ni)+1` (corner); `p < pm` takes `sgn·fld(ng+pm)`,
`p == pm` is zeroed for vectors, `p > pm` is the read-only source
(`rdb_ocean_fold.F90:250-258,286-294,336-344`). Self-conjugate slots: v column
`(ni+1)/2` for odd ni only; corners `p = 1` and `p = ni/2+1` (the bipoles)
(`rdb_ocean_fold.F90:85-90`). Same rule for the scalar metric copy
`metrics_fold_north_cv_scalar` (`rdb_ocean_metrics.F90:2042-2062`) and for the final
continuity mass flux (`rdb_continuity.F90:2769-2781`).

### 1.3 Gating and refusals

- `bc%north_fold` is rank-local: `= (tag == OBC_TRIPOLAR_FOLD) .and. has_north`
  (`rdb_ocean_boundary_types.F90:147-157`), set in `ocean_bc_validate_fold`
  (`:385-451`, assignment `:449`) and re-derived in `ocean_bc_state_set_edges`
  (`:589-612`, assignment `:611`), which the engine calls after `configure_ocean_bc`
  (`rdb_ocean_engine.F90:797-800`). Validate rules: north edge only, periodic w/e,
  `nghost ≥ 3` (`rdb_ocean_boundary_types.F90:406-446`).
- **px > 1 refusal**: `engine_setup` only (`rdb_ocean_engine.F90:347-364`).
  `validate_config` has no px-fold check (grep of `rdb_config.F90` for fold/px finds
  none). Auto-factor forces `px = 1, py = csize` for a folded grid
  (`rdb_ocean_engine.F90:330-336`).
- **North tile shorter than `nghost+1`**: `rdb_ocean_engine.F90:365-378` (checks
  `ny/py`, the smallest tile because the remainder goes to the first rows,
  `rdb_decomp.F90:93-101`).
- Test pin: `check_px_refused` asserts px = 2 is refused
  (`tests/mpi/test_ocean_tripolar_fold_mpi.F90:430-444`, called `:88`).
- Related exclusions that stay: `bt_halo > 0` with tripolar
  (`rdb_config.F90:7218-7224`, AUTO → 0 at `:7606-7608`); sea ice is single-rank
  (`rdb_config.F90:6653-6660`) and refuses the fold with dynamics (`:6765-6771`).
- Composition rule today: "periodic-x wrap FIRST, then the fold"
  (`rdb_ocean_fold.F90:93-94`, `rdb_ocean_fold_apply.F90:7-9`), and at step-time
  sites "exchange → periodic wrap → fold" (`rdb_ocean_halo_state.F90:101-106`).

## 2. Rank-pair mapping for px > 1

### 2.1 Tiles

`rank = ry·px + rx` (`rdb_decomp.F90:78-80,173-185`). Widths: `base = ni/px`,
`rem = mod(ni,px)`; ranks `rx < rem` get `base+1`, i.e. the remainder goes to the
WEST tiles (`rdb_decomp.F90:82-91`). Only the north rank row (`ry = py-1`,
`has_north`, `:107`) participates. Tile `rx` owns cells `[a_rx, b_rx]`.

### 2.2 Mirror in global physical indices

From §1.1, subtracting `ng` (all column indices reduced `modulo(·-1, ni)+1`):

| Stagger | column `c` → mirror | rows read (north tile's own storage) | rows written |
|---|---|---|---|
| T, v | `ni+1-c` | T: `ng+nyl+1-d`, d=1..ng | `ng+nyl+d` |
| u, corner | `ni+2-c` | v/corner: `ng+nyl+1-d`, d=1..ng, plus fold row `ng+nyl+1` (source) | `ng+nyl+1+d`, plus fold row (west half) |

(`nyl` = north tile's `ny_phys`.) Rows are identical on sender and receiver because
both are north-row tiles of the same height.

### 2.3 Owner routing

- T / v column `c`: owner = tile with `c ∈ [a, b]`.
- u / corner face `c` (west face of cell `c`, `c = 1 ≡ ni+1`): mirror face `ni+2-c`
  is the EAST face of cell `ni+1-c`. Under D1 the east seam face belongs to the
  cell's own tile (`rdb_ocean_halo.F90:17-24,28-33`), so
  **owner(u-mirror of face c) = owner(T-mirror of cell c)**. The sender reads its
  owned faces `[ng+2, ng+w+1]` (local), never the D1 copy at `ng+1`.
- Receiver column sets: T/v storage cells `[a-ng, b+ng]`; u/corner storage faces
  `[a-ng, b+ng+1]` — one extra face whose mirror owner is one cell further west.
  That is the whole "u differs by one" offset; the rank set can grow by one.

### 2.4 Cases (ng = 3)

| ni, px | tiles | mirror of physical cells | with ghosts |
|---|---|---|---|
| 30, 4 (rem 2) | [1,8] [9,16] [17,23] [24,30] | tile0 → [23,30]: rx2+rx3 (straddle); tile3 → [1,7]: rx0 | tile0 storage [-2,11] → senders rx0 (self), rx2, rx3 |
| 30, 3 (rem 0) | [1,10] [11,20] [21,30] | tile1 → itself (self-mirror) | tile1 → rx0, rx1, rx2 |
| 32, 3 (rem 2) | [1,11] [12,22] [23,32] | tile1 → [11,21]: rx0 + rx1 (self-mirror straddles) | |
| 10, 3 (rem 1) | [1,4] [5,7] [8,10] | tile0 → [7,10]: rx1 + rx2 | |

Rules: physical mirror of tile `rx` is exactly tile `px-1-rx` **iff `rem == 0`**
(else tile 0 is wider than tile `px-1` and the range straddles). Ghost columns
always reach the mirror tile's neighbours, and the receiver itself whenever its
ghost band wraps (every px ≥ 2). Peers per receiver ≤
`ceil((w_max+2ng+1)/w_min) + 1`; 3 is typical. Fold-row projection pairs are a
subset of the ghost-row pairs (same mirror columns).

The plan (per-peer index lists) is a pure function of `(ni, px, ng, stagger)`;
every rank computes every tile from `decomp_init` (`rdb_decomp.F90:53-109`), no
handshake. `base = 0` (px > ni) is refused.

## 3. Message layout

- **One message per ordered (sender, receiver) pair per call**; self-pair is a
  local gather/scatter through the same buffers (no MPI).
- **Columns:** per peer, a precomputed list `(src_local_col(e), dst_local_col(e))`,
  `e = 1..n_peer`, separately for the x-families {T,v} and {u,corner}. Lists are
  index lists, not ranges (a peer's set can wrap modulo ni). Built once, mapped to
  the device once.
- **Rows:** T, u: `ng` rows. v, corner: `ng+1` rows (the fold row is always sent;
  the receiver uses it only where `p < pm`). Packing the fold row unconditionally
  keeps senders stagger-agnostic; the waste is half a row.
- **Order:** field-major, then layer, then row, then entry:
  `off_f + ((L-1)·nrow_f + (d-1))·n_peer + e` — the existing halo's layer-outer
  convention (`rdb_ocean_halo.F90:565-566,615-635`).
- **Sign:** applied at unpack (`dst = s_f · buf`), `s_f = -1` for true vectors.
  Self-conjugate fold-row slots are set to 0 locally for vectors, no message.
  Reason: one pack kernel for every field; the sign is caller metadata (the same
  `negate` notion as `fold_north_corner_2d`, `rdb_ocean_fold.F90:311-313`); IEEE
  negation is exact, so placement cannot change bits.
- **Batching:** one message per peer per *group* (table below); per-field
  pack/unpack calls into one buffer at offsets (outer-shim for the tracer registry,
  as `ocean_halo_exchange_ml_state` does, `rdb_ocean_halo_state.F90:68-75`).

| Group | Fields | Sites |
|---|---|---|
| ml_state | h (T), u (u,−), v (v,−), every hTr (T) | §7 rows D1-D4, E1 |
| centre_3d | h, every hTr | C1, D5 |
| time_means | u_av, v_av, h_av | D3 |
| bt_mid | bt_eta (T), bt_ubt (u,−) | B1-B2 |
| bt_late | bt_vbt (v,−) | B3 |
| stress | tau_x (u,−), tau_y (v,−) | H1 |
| single | any one 2D/3D field | C2, E2-E5, S* |

- **3D:** all `nz` layers in the one message, as the batched halo (O3,
  `rdb_ocean_halo.F90:84-89`).
- **Tags:** one new tag `TAG_OC_FOLD = 15` (halo uses 11-14,
  `rdb_ocean_halo.F90:160-169`; `rdb_halo` 1-4). One message per ordered pair per
  call and every call completes with `waitall` before return, so `(source, tag)`
  is unambiguous, px = 2 included.
- **Size (OM4_025, px = 8, nz = 75, ng = 3):** ~186 columns × (3+3+4+3+3) rows ×
  75 ≈ 2.2·10⁵ values ≈ 1.8 MB per receiver for ml_state with T,S — small next to
  the regular halo.

## 4. Integration with the halo layer

- **New module** `src/comm/rdb_ocean_fold_exchange.F90` next to
  `rdb_ocean_halo.F90`, MPI only through `pic_mpi_lib`
  (`isend/irecv/waitall`, `HALO_ISEND_N`/`HALO_IRECV_N` as in
  `rdb_ocean_halo.F90:91-96,106-108`). The `no-mpi-in-rdb` hook enforces this
  (`.pre-commit-config.yaml:45-50`, `tools/no_mpi_in_rdb.sh`). Note:
  `src/comm/README.md` does not exist on main; CLAUDE.md is the only statement of
  the rule.
- **Pure plan module** `src/core/ocean/boundary/rdb_ocean_fold_plan.F90` (index
  math only, `pure`, no MPI) so the routing is unit-testable serially.
- **Init:** built by `ocean_halo_init` (or a sibling called next to it,
  `rdb_ocean_engine.F90:892-899`) from the same `decomp_t`, gated on
  `fold .and. px > 1`. Persistent buffers sized from the plan with the
  exit-delete/alloc/enter-create grow pattern of `ocean_halo_buffers_ensure_nz`
  (`rdb_ocean_halo.F90:450-506`); the tracer count is fixed once the registry
  locks.
- **Device vs host:** mirror the halo exactly — `device_resident` switch,
  `!$acc parallel loop` pack/unpack, `!$acc host_data use_device` around the MPI
  calls (`rdb_ocean_halo.F90:612-674`). Index lists are device-resident.
  **Observation (verify):** in `rdb_ocean_halo.F90` the device path always passes
  device addresses; `RDB_CUDA_AWARE_MPI` only selects a log line
  (`rdb_comm_env.F90:255-259`; define added at `cmake/dependencies.cmake:154-156`).
  I found no host-staged branch in the ocean halo. The fold exchange should inherit
  whatever the halo does rather than invent its own staging.
- **Dispatch:** `ocean_fold_wrap_state`, `ocean_fold_wrap_centre_3d_state`,
  `ocean_fold_wrap_eta_2d` (`rdb_ocean_fold_apply.F90:39-121`) keep their
  signatures and pick local kernel (px = 1) or exchange (px > 1). Add thin
  dispatchers (`ocean_fold_u/v/centre/corner`, 2D/3D, `negate`) for the sites that
  call raw kernels today (§7). Collectivity: every north-row rank has
  `north_fold = .true.` and reaches the same sites, so the exchange is collective
  over exactly that row; other ranks skip, as today.
- **BT group exchanges** stay as they are (`ocean_halo_bt_group_2d`,
  `rdb_ocean_halo.F90:1573-1592`; called `rdb_barotropic_substep.F90:1676`); the
  fold is added beside them (§5.3).

## 5. Ordering

### 5.1 Rule for px > 1

At every seam site: **halo exchange (X pass incl. periodic MPI wrap, then Y pass)
→ local periodic wrap on undecomposed axes only (y, if any) → distributed fold.**
The halo runs X then Y and the Y pass never touches the north ghosts of the north
tile (`need_n` false there, `rdb_ocean_halo.F90:1768`). Because the fold is
owner-routed, it writes every storage column of the north ghost rows (corner ghosts
north × x-ghost included) directly; it does not read x-ghosts. A later X pass copies
a neighbour's already-folded values into those corners, which are the same bits.

### 5.2 Why this is bitwise the serial result

Every fold write is `±(an owned value)`; the serial fold writes `±` the same owned
value, read directly or through a periodic-wrap copy, which is a bitwise copy.
Copies and exact negations commute, so the order relative to the X pass does not
change bits. The one assumption already made by the MPI halo: the two copies of the
periodic seam u face are bit-equal (`rdb_ocean_halo.F90:33-40`).

### 5.3 Barotropic fast loop (needs care)

The serial folds are inline inside `!$acc kernels async(1)` regions, interleaved
with the passes: η after the η update (`rdb_barotropic_substep.F90:1098-1106`),
ubt after Pass 2b (`:1397-1404`), vbt after Pass 2c and BEFORE the time-mean
accumulators (`:1620-1642`, accumulators `:1644-1657`). Regions: `:799-1420`,
`:1433-1658`. For px > 1 the inline blocks are skipped
(`do_fold .and. .not. ocean_halo_is_decomposed_x()`) and replaced by:

1. **bt_mid** (η + ubt) right after the mid-substep u exchange
   (`:1421-1432`, which already does `!$acc wait(1)`). This moves the η fold from
   before Pass 2b to after it. Correct only if Pass 2b (`:1174`ff) reads no η north
   ghost that feeds an owned value. **Not verified** — if it does, split the first
   region at `:1106` and fold η there (a third exchange).
2. **bt_late** (vbt) after `:1642`, by splitting the second region so the
   accumulators sum projected values (the fold-row v is an owned face).
3. The end-of-substep group exchange (`:1659-1680`) stays where it is.

`bt_halo > 0` stays excluded for tripolar (`rdb_config.F90:7218-7224`).

### 5.4 Engine init

Today: wrap → fold (`rdb_ocean_engine.F90:830-891`) → `ocean_halo_init` (`:894`)
→ host halo (`:901-916`). The distributed fold needs the plan, so for px > 1 the
init-time folds of `b`, the prognostics (cold start only, `:840-853`), `z_draft`,
`cover_frac` and `bt_H_ref` move after `:916`, host mode. The px = 1 path is left
untouched.

## 6. Ownership of the fold-row DOFs

- The east-half slot (`p > pm`) is authoritative and its owner never modifies it;
  the west-half owner receives it and negates (vector) or copies (scalar); `p == pm`
  is zeroed locally (vector) or left (scalar). This is the existing serial rule
  (§1.2). With two ranks each holding half the row, they do not "agree" by
  computing the same expression: one value moves, exactly.
- West-half slots in a receiver's x-ghost columns are overwritten by the fold too
  (the X pass may have delivered a pre-projection neighbour value). East-half ghost
  slots keep the X-pass value, which is the owner's unmodified value.
- No source of the fold is ever a destination (sources: rows below the fold line
  and the east half of the fold row), so pack-all → MPI → unpack-all has no hazard
  even in place.
- If anyone later switches the projection to an average, the operand order must be
  canonical (east operand first) on every rank; flag that in the kernel header.

## 7. Inventory of fold call sites (`grep fold_north_|ocean_fold_wrap|north_fold`)

| # | file:line | field(s) | kind | host/device | px > 1 action |
|---|---|---|---|---|---|
| D1 | `rdb_ocean_dyn.F90:3627-3632` post-continuity | h, u, v, hTr | T, u−, v− | device | dispatcher |
| D2 | `rdb_ocean_dyn.F90:4228-4234` stage entry | h, u, v, hTr | T, u−, v− | device | dispatcher |
| D3 | `rdb_ocean_dyn.F90:4283-4298` pred_corr means | u_av, v_av, h_av | u−, v−, T | device | time_means group |
| D4 | `rdb_ocean_dyn.F90:5073-5088` stage end | h, u, v, hTr | T, u−, v− | device | dispatcher |
| D5 | `rdb_ocean_dyn.F90:5092-5112` `refresh_tracer_ghosts` (callers `:3491,:4949,:5066`) | h, hTr (via centre_3d_state) | T | device | dispatcher |
| C1 | `rdb_continuity.F90:2717-2721` mid Lie split | h, hTr | T | device | dispatcher |
| C2 | `rdb_continuity.F90:2769-2781` | mass_flux_y_layer | v− (projection) | device | single; no halo precedes it — owner routing required |
| C3 | `rdb_continuity.F90:3645-3894` drain (`:3908-3954`) | uhtr, vhtr, tr_work, pal/par/pa6, … | T/u−/v− | device | none: `dt_tracer_advect_ratio > 1` is single-rank (`rdb_ocean_engine.F90:379-385`) |
| B1 | `rdb_barotropic_substep.F90:1098-1106` | bt_eta | T | device (inline) | bt_mid |
| B2 | `rdb_barotropic_substep.F90:1397-1404` | bt_ubt | u− | device (inline) | bt_mid |
| B3 | `rdb_barotropic_substep.F90:1620-1642` | bt_vbt | v− + projection | device (inline) | bt_late |
| H1 | `rdb_ocean_halo_state.F90:152-160` (configure `rdb_ocean_setup.F90:661`, `rdb_ocean_data_forcing.F90:228`) | tau_x, tau_y | u−, v− | both | stress group |
| E1 | `rdb_ocean_engine.F90:853` | ml state | T, u−, v− | host | after `:916` |
| E2 | `rdb_ocean_engine.F90:851,869,876,887` | b, z_draft, cover_frac, bt_H_ref | T | host | after `:916` |
| S1 | `rdb_ocean_state.F90:1761-1799` `seed_wrap_static_2d` (calls `:1305,:1332-1333`) | b, z_draft, cover_frac | T | host | already skips px > 1 (`:1788-1789`); see §11 |
| S2 | `rdb_ocean_metrics.F90:821` in `metrics_apply_land_mask` (`north_fold` from `rdb_ocean_setup.F90:356-366`) | wm (wet-mask copy) | T | host | pass `north_fold=.false.` + a pre-folded copy |
| S3 | `rdb_ocean_setup.F90:772-799` `fill_f_corner_seam_ghosts` | f_corner | corner, copy | host | returns early on px > 1 (`:787`); needs corner fold for beta-plane (§11) |
| M1 | `rdb_ocean_metrics.F90:1897-2062` (via `:1375,:1388`, `metrics_fill_tripolar` `:1672-1733`, engine staged `rdb_ocean_engine.F90:710-713`) | all metrics, Cu/Cv/Bu scalars | all | host | none: built on the whole grid / full-width band then windowed (`rdb_ocean_metrics.F90:1694-1709,1258-1280`); staged path is single-rank (`rdb_ocean_engine.F90:695-700`) |
| A1 | `rdb_ocean_api.F90:1538-1559` | b, bt_H_ref | T | host | dispatcher if the API runs multi-rank (unverified) |
| L1 | `rdb_ocean_mle.F90:575-584` | logic only (`n_seam`) | — | — | none |

Not fold sites (grep-verified): diagnostics (`src/core/ocean/diag/`), restart/io
(`src/core/ocean/io/`), ALE (`src/ALE/`), sea ice (`src/core/ice/`; EVP declares
no fold, `rdb_ice_evp.F90:93`).

## 8. Tests

### 8.1 Unit: `tests/mpi/test_ocean_fold_exchange_mpi.F90` (new)

- Synthetic index-encoded fields, exact in `wp`:
  `val = fid·10⁸ + gi·10⁵ + gj·10² + k` (gi, gj global physical, reduced periodic).
  Fill owned cells only, poison ghosts with NaN, run the exchange, check every
  north-ghost and fold-row slot against the analytic mirror of §2.2 with sign
  (`−val` for vectors, `0` at self-conjugate vector slots, untouched east-half fold
  slots). Also compare the whole storage window against a reference: the global
  array on every rank, periodic-wrapped then folded by the existing local kernels
  — the pattern of `tests/mpi/test_halo_ocean_mpi.F90:10-17`.
- Matrix: T/u/v/corner × 2D/3D × negate on/off × groups; ni even and odd
  (v self-conjugate column); `rem = 0` and `rem > 0`. Branch on `nprocs` as the halo
  test does (`tests/mpi/test_halo_ocean_mpi.F90:3-22`): np 1 (dispatcher must hit
  the local kernels bit-identically), np 2 (2x1), np 3 (3x1, self-mirror middle
  tile, ni = 30 and 32), np 4 (4x1 with ni = 30 → 8/8/7/7, and 2x2). Host mode
  (`device_resident=.false.`) like the halo test; device path covered end-to-end
  by 8.2 on the GPU build.
- Serial unit test of the plan module (no MPI): brute-force owner search vs the
  closed-form lists over a sweep of `(ni, px, ng)`.
- CMake: add to `MPI_TESTS` (`tests/CMakeLists.txt:580-596`) plus 1/3/4-rank legs
  in the style of `:687-695`.

### 8.2 Engine: `tests/mpi/test_ocean_decomp_bitid_mpi.F90`

- Add a `tripolar` case to `CASES` (`:120-124`), namelist copied from the tripolar
  test (`tests/mpi/test_ocean_tripolar_fold_mpi.F90:107-137`) with `nx = 30,
  dx = 12.0` (not `NX_G = 26`: `360/26` is inexact; `compare` is size-agnostic,
  `:458-509`), cap ≥ 73N (`test_ocean_tripolar_fold_mpi.F90:31-36`), both schemes.
  `run_case` already loops every factorisation (`:549-551`): np 2 → 2x1, 1x2;
  np 4 → 4x1 (8/8/7/7, uneven), 2x2, 1x4. ctest already runs 1/2/4 ranks
  (`tests/CMakeLists.txt:589-596,687-695`). This test compares OWNED windows only.
- Extend `test_ocean_tripolar_fold_mpi` (it compares whole storage, ghosts
  included, `:281-341`): replace `check_px_refused` (`:430-444`) by decomposed runs
  at px = 2 (and px = 4 / 2x2 on np 4); add a 3-rank leg next to
  `tests/CMakeLists.txt:674-682` for odd px.

## 9. Open design decisions

| # | Question | Options | Recommendation |
|---|---|---|---|
| 1 | Sign flip | send / receive | Receive (§3). |
| 2 | Source of values | owner (D1) / any held copy after X pass | Owner: independent of X pass; needed for C2. |
| 3 | One generic exchange vs per-stagger | generic + stagger enum / 4 routines | Generic engine, stagger enum picks list family + row count; thin typed wrappers. |
| 4 | Keep local kernels | px = 1 only / also self-mirror / drop | px = 1 only. px > 1 always uses the exchange (self-pair = local gather), one code path per regime. |
| 5 | Pack location | host / device | Both, same `device_resident` switch as the halo. |
| 6 | Batch with regular halo | 3rd pass inside primitives / separate call | Separate in v1 (C2 and BT have no halo at the fold point; the primitives cannot tell scalar from vector face fields). Revisit for performance. |
| 7 | Fold row in v/corner messages | always / separate message | Always (`ng+1` rows). |
| 8 | BT fold points | 3 exchanges at serial points / 2 (η moved to mid) | 2 if the Pass-2b η read audit passes, else 3. |
| 9 | Auto-factor for tripolar | keep `px = 1, py = N` / perimeter-minimising | Keep px = 1 default until px > 1 is measured; px > 1 explicit. |
| 10 | Testing while refused | optional `engine_setup` argument (test-only) / one big PR | Test-only optional dummy argument `allow_distributed_fold`, removed by the last PR. Not a namelist knob. |
| 11 | Minimum tile width | none / `w ≥ ng+1` | Refuse `w_min < ng+1`, symmetric with the height rule; I found no width guard in the regular halo. |
| 12 | Whole-grid metric build per rank (analytic tripolar) | accept / tile-local construction | Accept for v0.1.0 (transient, host, setup only; `rdb_ocean_metrics.F90:1694-1709`). A file supergrid already reads a full-width band. |
| 13 | Is px > 1 needed for v0.1.0? | px = 1 strips (supported now) / px > 1 | Out of scope to decide here. Strips already run: 1080/16 = 67 rows per rank, but each rank's halo is 2 × 1440 columns vs 2 × (360+270) for 4x4 tiles. |

## 10. PR slicing (one PR each, in order)

| PR | Scope | Files | Tests | Still refused after | Size |
|---|---|---|---|---|---|
| 1 | Fold plan: pure index math, per-peer lists, stagger families, self-conjugate masks | new `rdb_ocean_fold_plan.F90`; CMake | serial plan unit test | px > 1 (engine) | S |
| 2 | Exchange engine: buffers, host/device pack, tags, groups, dispatch by px | new `src/comm/rdb_ocean_fold_exchange.F90`; halo init hook | `test_ocean_fold_exchange_mpi` at np 1/2/3/4 | px > 1 | M-L |
| 3 | Baroclinic + continuity sites D1-D5, C1-C2, H1; fold_apply dispatchers; engine init reorder (E1-E2) for px > 1; test-only setup bypass (decision 10) | `rdb_ocean_fold_apply.F90`, `rdb_ocean_dyn.F90`, `rdb_continuity.F90`, `rdb_ocean_halo_state.F90`, `rdb_ocean_engine.F90` | tripolar bitid with `n_inner = 0` if that mode runs a tripolar case (unverified), else component tests on a decomposed state | px > 1 | M |
| 4 | BT fast loop B1-B3: skip inline folds on px > 1, bt_mid/bt_late exchanges, region splits; Pass-2b η audit | `rdb_barotropic_substep.F90` | `test_ocean_tripolar_fold_mpi` at px = 2 via bypass | px > 1 | M |
| 5 | Setup sites S2, S3 (corner stagger, host), A1 audit, seed/consumer gap (§11) | `rdb_ocean_setup.F90`, `rdb_ocean_metrics.F90`, `rdb_ocean_state.F90`, `rdb_ocean_api.F90` | beta-plane tripolar case in the tripolar MPI test | px > 1 | M |
| 6 | Lift refusal: delete `rdb_ocean_engine.F90:355-364`, add width rule (decision 11), remove bypass; invert `check_px_refused`; bitid tripolar case; docs | engine, both MPI tests, `tests/CMakeLists.txt`, `rdb_ocean_fold.F90` header `:8-22`, `docs/CAPABILITIES_AND_LIMITATIONS.md:301,318,1009-1029`, CLAUDE.md "Horizontal grids" / "Tripolar fold under MPI" / single-rank list | full bitid 1/2/4 + tripolar 1/2/3/4 | north tile `< ng+1` rows, tile width `< ng+1`, `bt_halo > 0`, sea ice | M |
| 7 | (optional) Overlap / batching with the regular halo, OM4_025 timing on GPUs | comm modules | unchanged gates | — | M |

Sea ice, wet/dry and the other single-rank features (CLAUDE.md "Single-rank features
fail loud") are unaffected and stay refused on > 1 rank throughout.

## 11. Facts I could not verify (resolve in the PR that touches them)

- **Pass 2b η reads** (§5.3): whether any owned value between
  `rdb_barotropic_substep.F90:1106` and `:1420` reads η north-ghost rows. Decides
  2 vs 3 BT exchanges.
- **Seed-time gap:** `seed_wrap_static_2d` skips the fold on px > 1
  (`rdb_ocean_state.F90:1788-1789`) while setup-time consumers read those ghosts
  before the engine's halo pass (`rdb_ocean_state.F90:1293-1305`). The same gap
  exists for x-ghosts on decomposed periodic runs, which pass bitid; whether
  north ghosts are equally harmless is not shown.
- **Configure-time stress refresh** (`rdb_ocean_setup.F90:661`, via
  `configure_ocean_forcing` at `rdb_ocean_engine.F90:734`) runs before
  `ocean_halo_init` (`:894`); what the halo primitives do then on > 1 rank was not
  traced. The fold there needs the plan; either init the plan earlier (it depends
  only on `decomp` and `nghost`) or rely on the first step's refresh.
- **Beta-plane f on tripolar:** `metrics_fill_coriolis` evaluates the beta-plane
  formula in the north ghost rows (`rdb_ocean_metrics.F90:2433-2441`); planetary f
  is read from windowed, already-folded `geolatBu` (`:2421-2426`). With the early
  return at `rdb_ocean_setup.F90:787`, beta-plane f_corner ghost rows are unfolded
  on px > 1. Whether beta-plane tripolar is a supported configuration is unknown.
- **CUDA-aware vs host-staged** halo path (§4).
- **API multi-rank** (`rdb_ocean_api.F90:1538-1559`): whether this entry point is
  reachable on > 1 rank.
- **`refresh_tracer_ghosts` re-folds `h_layer` too** (`ocean_fold_wrap_centre_3d_state`
  folds h, `rdb_ocean_fold_apply.F90:99-100`), despite its docstring saying thickness
  is left alone (`rdb_ocean_dyn.F90:5093-5097`). It is idempotent, but in the
  distributed path it costs an extra field per message. Keep identical semantics
  in the port; fix separately if wanted.
- **`n_inner = 0` tripolar** (PR 3 test plan): not checked whether that mode is a
  valid tripolar configuration.
