# Ocean dynamical core — scaffold map

This directory hosts the ocean-regime path (`sim_type='ocean'`) — the
only regime this repository builds. The goal is **MOM6-style "organise
by concern" with a single composed god state** so a contributor can pick
up one slot and ship without spelunking the rest of the tree.

*Reconciled with the tree 2026-08-23 (coastal carve-out): the god-state
access path is `ocean_state`, not `state%ocean`; references to the
coastal A-grid / unstructured path were removed. The four rules, the
slot map, and the `is_init` / `scratch_3d_buffer_t` / OpenMP-only /
no-direct-MPI conventions are unchanged.*

Slots marked **✓** below have real allocations + bound `enter_data` /
`exit_data` + kernel(s) + tests landed (on the corresponding
`feat/ocean-*` feature branch, pending merge to main). Everything else
is still a Phase 0e shell: the type is declared, the component arrays
are named, `init / destroy` are no-ops that compile and run, and the
slot is wired into `ocean_state_t`. The "who reads / who writes" column
is the hand-off contract for the contributor implementing each shell.

**Implementation status (2026-05-16)** — highlights pending merge:
- Tier-1 prognostic stack: barotropic + multilayer C-grid state, continuity-PPM, Sadourny Coriolis-adv (enstrophy default / energy-conserving `sadourny_energy` / Arakawa-Hsu `sadourny_hk`), split-RK2 driver, Wright EOS, FV-PGF, vertical advection, hvisc, hdiff-tracer, vmix (PP81 + KPP shear/convective/non-local), bottom drag, 2D wind stress.
- Phase 5a Leith closure (`feat/ocean-leith`) — `ah_face_*` from vorticity-gradient magnitude with `ah_max` clip; hvisc kernel optionally reads per-face coefficients.
- Phase 5g vertical coord + ALE remap (`feat/ocean-vcoord-*`, six stacked) — `ocean_vcoord_t` slot with `compute_target_h` for EULERIAN_Z / SIGMA / ZSTAR / ZSIGMA / ZSTAR_SIGMA / ZSTAR_FULL / Z_FIXED; centre + face PPM remap orchestrator wired into the split driver as optional `vcoord=`. The isopycnal `VCOORD_RHO` (P2, validation-grade) is the one state-dependent coord — a separate `compute_target_h_rho(this, total_h, eta, T, S, eos)` TBP (the shared `compute_target_h` stays `pure (total_h, eta)`); the remap orchestrator passes `eos=` only for RHO and builds T/S concentrations into `remap_conc_{t,s}` scratch.
- Phase 6a–6d diagnostics manager (`feat/ocean-diag-*`, four stacked) — registry + cadence + procedure-pointer fills (SSH, T, S, u, v, KE) + layer→z linear remap + MEAN/MAX/MIN time-ops + serial NetCDF emit.  Driver wire-up still pending.
- Implicit stress/drag fold (`feat/ocean-implicit-stress-drag`, `&ocean_vdiff_nml`, default off ⇒ bit-identical) — folds the surface wind stress (`implicit_stress` → vdiff `k=nz` RHS row, Neumann top-BC) and bottom drag (`implicit_drag` → vdiff `k=1` diagonal, stress bottom-BC) into the backward-Euler vertical-friction tridiagonal, replacing the explicit pre-solve adds; thin-layer (z*/ZSTAR_FULL pinch-out) quadratic-drag/wind CFL-robust.  No new state slot — reuses the bottom-drag slot's `lambda_bot_u/v` + the existing vdiff tridiagonal workspaces; lives in `rdb_ocean_vdiff.F90` + `rdb_ocean_bottom_drag.F90`.  Configure-time fail-loud exclusions: `implicit_drag` vs `&ocean_bdrag_nml implicit` + HBBL (`hbbl>0`); `implicit_stress` vs `&ocean_vmix_nml direct_stress`.  Test: `test_ocean_vdiff_implicit_stress_drag`.
- Ice-shelf top drag (`feat/cavity-top-drag`, `&ocean_tdrag_nml`, default off ⇒ bit-identical) — a quadratic / linear momentum sink at `k = nz` on ice-covered FACES, the mirror of the bottom drag; its own slot `ocean_top_drag_t` in `rdb_ocean_top_drag.F90`.  **Where it sits in the step:** the compute runs with the other slow velocity tendencies (`run_stage` / `run_stage_split`, immediately after `ocean_channel_drag_compute_tendencies` and before the surface stress); the tendency is added into `bt_work%F_slow_u/v` by `add_top_drag_into_F_slow` BEFORE the depth mean that drives the barotropic substep — the same route the bottom drag takes, and the reason it damps the fast mode at all; and the explicit apply runs in the batched velocity-apply chain on OpenACC queue 1, next to the bottom-drag apply, skipped when the rate is folded into the vdiff `k = nz` diagonal instead.  Requires `&ocean_cavity_dyn_nml enable`; ONE `C_d` shared with `&ocean_cavity_melt_nml cdrag_top` (disagreement is fail-loud).  **Two implicit forms, and they are different things:** `&ocean_tdrag_nml implicit` is backward-Euler inside the drag kernel (compatible with `htbl`), while `&ocean_vdiff_nml implicit_top_drag` folds `dt·λ_top` into the vdiff `k = nz` DIAGONAL — the mirror of `implicit_drag` at the other end of the column, layer-`nz` only, and the path on which the wind-stress RHS is masked by `(1 − cover)` (no atmosphere under a shelf).  Mutually exclusive, fail-loud.  Test: `test_ocean_top_drag`.
- Tidal body forcing (`feat/ocean-tides-body`, `&ocean_tides_nml enable`, default off ⇒ bit-identical) — fills the previously-dormant `ocean_tides_t` slot (C1 of the tides chain). Astronomy generator `rdb_ocean_tide_astro` (mean longitudes → equilibrium argument `V_c` with ±π/2 diurnal signs → nodal `f/u`; 10-constituent catalog verified 1:1 vs `MOM_tidal_forcing.F90`). Slot carries `cos_struct/sin_struct(nx,ny,3)` (built at init from `geolatT/geolonT`) + `eta_eq(nx,ny)`; `tides_update_eta_eq` = per-outer-step host scalar update + a `do concurrent` `_impl` fill (2 mul + 1 add per constituent per cell). Consumed in `barotropic_substep_{linear,nonlinear}` via an optional explicit-shape `eta_forcing` folded into the four η-gradient PGF sites (`∇(bt_eta − eta_eq)`), held static across the inner substep loop. `enter_data`/`exit_data` TBPs wired into the orchestrator (select-type→`_impl`). Fail-loud on cartesian grid. Tests: `test_ocean_tides_astronomy`, `test_ocean_tidal_forcing`, `test_ocean_tides_disabled_bitident`.
- Scalar SAL (`feat/ocean-tides-sal`, C2, `&ocean_tides_nml use_sal`/`beta_sal`, default off ⇒ bit-identical) — self-attraction & loading `η_sal = β·η` folded into the C1 seam: `tides_update_eta_sal` fills `eta_sal = beta_sal·bt_eta` (lagged one outer step) and the combined `eta_forcing = eta_eq + eta_sal`; the barotropic momentum then feels the effective-gravity `−g(1−β)∇η`. New `eta_forcing`/`eta_sal` slot arrays (mapped in enter_data). Accad-Pekeris load form (body tide full strength) — documented divergence from MOM6's `(1−β)(η−η_eq)`. Test: `test_ocean_tides_sal`. C4 internal-tide drag fields still declared-but-unmapped — next PR on this seam.
- OBC boundary tides (`feat/ocean-obc-tides`, C3, `&ocean_bc_nml obc_tidal_nodal`, default off ⇒ bit-identical) — folds the nodal factor `f_c` + equilibrium/nodal phase `(V_c+u_c)` into the OBC v2 per-edge tidal `eta_target` sum in `barotropic_substep_nonlinear`, reusing C1's `rdb_ocean_tide_astro` generator and the shared `&ocean_tides_nml` reference epoch (interior body tide ↔ boundary tide stay phase-consistent). New fixed-size members `tidal_fnodal`/`tidal_arg` on `ocean_bc_face_tag_t` (ride the parent `copyin(this)` — no separate map) + scalar `tidal_nodal` on `ocean_bc_state_t`; `pure obc_match_constituent` (nearest catalog ω, rel-tol 1e-4, 0 ⇒ fail-loud) + `pure obc_tide_nodal_fill`; host precompute in `configure_obc_edge_nodal` at configure. Nodal-ON uses the MOM6 `OBC_TIDE` lag convention (`cos(ω·t + V_c+u_c − φ_c)`); nodal-OFF is byte-identical to the legacy `cos(ω·t + φ_c)`. Test: `test_ocean_obc_tide_nodal`. C4 fields still declared-but-unmapped.
- Reachability pack (PR-9) — five finished-but-unreachable capabilities given their missing namelist spellings (no new physics, all default off/unchanged ⇒ bit-identical): **PQM** vertical remap (`&vcoord_nml remap_method="pqm"`, `&ocean_diag_nml diag_remap_scheme="pqm"`; `remap_column_pqm` was already dispatched, `nz<5` falls back to PPM); **density-coordinate diagnostics** (`&ocean_diag_nml vgrid="density"` + `rho_levels`/`n_rho_levels`, or per-diagnostic `:density`/`:rho` in `diags`; configure-time guard requires a strictly-increasing list — density bins have no sigma/z*-style auto-fill); **`DIAG_OP_INTEGRAL`** (`:integral` in `diags` — Σ(sample·dt) undivided, the budget-closure operator, was already dispatched in `fold_sample_unmasked_impl`/`ocean_diag_step`); **`&ocean_hdiff_nml kappa_h`** (the along-coordinate tracer Laplacian `tracer_hdiff` was called every RK2 stage but always short-circuited on `kappa_h<=0` with no way to set it; landed alongside a physical-edge (not array-edge) wall-closure fix — `tracer_hdiff_one_impl` now zeroes the flux at `nghost+1`/`nghost+nx_phys+1`, matching `continuity_zonal_flux`'s `has_*`/`OBC_WALL` seam gate); and **`&ocean_vmix_nml pp81_*`/`kpp_*`** (PP81 interior + KPP BL-depth constants — the `&nonhydrostatic_nml kpp_*` keys existed but routed to the (since carved-out) coastal state only; the ocean copy re-derives `kv_bg`/`kt_bg`/`ks_bg` from the new `pp81_nu_bg`/`pp81_kappa_bg` via `vmix_seed_backgrounds`, extracted out of `ocean_vmix_init` so both the initial seed and the config-copy stay in lockstep). See `docs/CLOSURE_MATRIX.md` for the physics + knob tables.

## God state (ish)

`ocean_state : ocean_state_t` (file: `core/ocean/state/rdb_ocean_state.F90`).

Declared and initialised in `driver_run_ocean` (`src/driver/rdb_driver.F90`)
and passed down by argument.  `sim_type` is pinned to `'ocean'` — the coastal
A-grid / unstructured regimes and the `state_t` wrapper that used to compose
this slot as `state%ocean` were split into their own repository, so the ocean
god-state is now the top-level state object rather than a component of one.

## Design contract

`ocean_state_t` is the single source of truth for the ocean path, but
it is deliberately **not a global**.  Four rules keep it from sliding
into a god-object:

1. **Handle, not global.** `ocean_state` is passed as a subroutine
   argument.  No module ever does `use rdb_ocean_state, only:
   the_one_state` or stashes it in a `save :: ...` variable.
   Searchable property: `grep -r "ocean_state%" src/` should only
   match argument call sites; `grep -r "save :: .*state" src/` should
   return nothing.

2. **Each slot owns its memory.** Every allocation is attached on the
   slot's `init`, every GPU mapping in its `enter_data`, every
   release in `destroy` / `exit_data`.  Memory ownership is per-slot,
   not centralised — a contributor adding a new slot does not edit
   `rdb_ocean_state`'s init except to add the one `call
   this%newslot%init(grid)` line.  **Default-off closures are gated,
   though:** EPBL, kappa-shear, tidal-mixing, GM, Redi, MLE, VarMix,
   MEKE and the isopycnal slopes only allocate when their `enable` flag
   is set, so a plain run pays none of their (multi-GB at scale)
   footprint.  Their `enable` flags are latched in `init_from_config`
   **before** `init` so the `if (this%slot%enable) call
   this%slot%init(...)` gate can read them; the matching
   `enter_data`/`exit_data` calls in the parent walk are gated the same
   way (those bodies carry no internal `allocated` guard).  The runtime
   is safe with a gated-off slot because its kernels already
   `if (.not. enable) return` and `destroy` is `allocated`-guarded.
   `ocean_state%bytes()` (reported before `enter_data`, reconciled
   against the measured device mapping — see `rdb_mem_report`) sums an
   `arr_bytes` term per array, so a gated-off slot counts 0 and the
   upfront memory report tracks the gating directly.

   **Gating WITHIN a live slot** (a scheme-variant scratch buffer that
   only one `variant` can reach — `ocean_pressure_force_t`'s FV_MOM6 /
   reconstruction stack, `continuity_t%windowed_advection`) follows the
   same shape: latch the discriminator in `init_from_config` before
   `init`, and let `bytes()` ride `arr_bytes` so the count stays
   conditional in the same way the allocation is.  Two extra
   obligations, because an intra-slot gate has no `if (.not. enable)
   return` backstop:
   - **Prove unreachability, per buffer, at the call site.** Under
     `-gpu=...,mem:separate` an unallocated array reaching a `do
     concurrent` kernel does NOT reliably fault — nvfortran emits an
     implicit `copyin` for explicit-shape dummies, so you get stale host
     data and a plausible wrong answer. Structural proofs only: a
     dispatch that `return`s before the buffer is named, or a consumer
     that already `error stop`s on the wrong variant.
   - **Default the gate OFF (allocate everything).** Direct
     `slot%init(...)` call sites in tests and benchmarks set the
     discriminator *after* `init`, so a gate that defaults on would
     decide on a stale value. `init_from_config` turns it on; anything
     else keeps the historical behaviour. Where the discriminator is
     re-derived later (`configure_ocean_*`), compare it against the
     latched value and `error stop` on a mismatch.

3. **Read-mostly across slot boundaries.** Coriolis-advection reads
   `barotropic%u_face_x` but does not mutate it.  The split-RK2
   driver (`ocean_dyn`) is the only thing that mutates prognostic
   state.  Kernels stay diag-agnostic for the same reason — the diag
   manager reads state via the registry, not the other way around.

4. **MPI is only ever spoken via `pic_mpi_lib`.** Every
   `rdb_ocean_*` module operates on the local-subdomain
   `ocean_state` and is genuinely MPI-agnostic — the same source
   compiles single-process or multi-rank, and `grep -rE "use
   (mpi|mpi_f08)" src/` returns empty (enforced by the
   `no-mpi-in-rdb` pre-commit hook).  Even the `src/comm/`
   backend goes through the `pic_mpi_lib` wrapper from the pic
   dependency rather than calling MPI directly.  Halo exchanges
   and inter-rank reductions are dispatched via `src/comm/`, which
   has ONE implementation -- no stub twin -- and runs single-rank by
   taking its local paths (`src/comm/README.md`).  Phase 1 will grow
   C-grid halo entry points there (`halo_exchange_face_x_*`,
   `halo_exchange_face_y_*`, `halo_exchange_corner_*` for PV at
   corners); the kernel-side dispatch sites stay identical
   regardless of backend.

These rules are how the handle scales: collaborators can implement a
slot without reading the rest of the tree.  A future portability lint
(`tools/ocean_state_lint.py`) will enforce them statically.

## Slot map

| Slot | Type | File | Phase | Reads from | Writes to |
|---|---|---|---|---|---|
| C-grid barotropic state ✓ | `barotropic_state_t` | `core/rdb_barotropic_state.F90` | 1 | (allocated by) | `h`, `u_face_x`, `v_face_y`, `mass_flux_*`, RK saves |
| C-grid multilayer state ✓ | `multilayer_state_t` | `core/rdb_multilayer_state.F90` | 1 | (allocated by) | per-layer `h_layer`, `u_face_x_layer`, `v_face_y_layer`, vertical `w`, shared `tracers(:)` registry (S, T, optional age/pseudo-salt + any `register_passive_tracer`-appended slot). `register_passive_tracer(grid, name, units, long_name, idx)` grows the registry (6-arg, `idx=0` on refusal); `registry_locked` (set by `enter_data`, cleared by `exit_data`) refuses a post-map registration instead of leaving an unmapped `hTr` on the `mem:separate` GPU build. Budget attribution is `tracer_t%budget_id` (`TRACER_BUDGET_NONE`/`_HEAT`/`_SALT`), not an index comparison — every ocean budget-dispatch site `select case`s on it. Recipe: `docs/howto/add_passive_tracer.md` "Ocean C-grid" section |
| Pseudo-salt verification tracer ✓ (default off) | (free procedures, `rdb_ocean_pseudo_salt` module — no new `_t`; rides the multilayer registry) | `../../tracer/rdb_ocean_pseudo_salt.F90` | — | `multilayer.tracers(idx_salinity)%hTr` (seed), surface salt flux, KPP/EPBL non-local `gamma_s` | `multilayer.tracers(idx_pseudo_salt)%hTr`; `pseudo_salt`/`pseudo_salt_diff` diagnostics. `&ocean_tracers_nml enable_pseudo_salt`; first real consumer of `register_passive_tracer` — seeded to S, given exactly S's surface flux + KPP nonlocal mirror (self-gated inside the kernels that own them, no new dyn-step call site), everything else rides the registry for free. Deviation `pseudo_salt − S` measures the passive-vs-active transport-path error. Fail-loud excluded from SSS restoring + sea-ice (both un-mirrored salinity sources) |
| Split-RK2 driver | `ocean_dyn_t` | `dynamics/split_rk2/rdb_ocean_dyn.F90` | 4 | tendencies from continuity/coriolis/pressure/vmix/lateral | `ubt_sum / vbt_sum / eta_sum` time-mean accumulators, advances state |
| Continuity-PPM ✓ (barotropic + windowed tracer advect) | `continuity_t` | `kernels/continuity_ppm/rdb_continuity.F90` | 2 / P2 | `barotropic.u_face_x`, `multilayer.h_layer` | per-face `mass_flux_x_layer`, `mass_flux_y_layer`; accumulator slots `uhtr`/`vhtr` (face transport m³, ½-weight per RK2 stage) + `t_dyn_rel_adv` (elapsed time since last drain); `continuity_tracer_drain` spends them via swept-average CW-PPM with fixed-budget CFL sub-cycling (MOM6 `DT_TRACER_ADVECT`; `dt_tracer_advect_ratio` knob in `&ocean_vmix_nml`). **Positive-definite continuity** (`&ocean_continuity_nml positive_definite`, default off ⇒ bit-identical): a `2·h_lim` PPM edge floor + a per-donor θ outflux limiter scale the folded `mass_flux_*_layer` so every layer stays `h ≥ h_lim` with **zero mass created** (contrast: MOM6's injecting `max(h,Angstrom)` clamp is NOT ported; the `conservative_floor` borrow stays the backstop). `h_lim = angstrom_h` on VCOORD_LAGRANGIAN else 0; D3 single-source scaling keeps CWC exact; fail-loud vs `&ocean_wetdry_nml enable`; per-call `n_limited_step` + int64 `n_limited_total` drained to the console. See `docs/CLOSURE_MATRIX.md` |
| PV-conserving Coriolis+adv ✓ (Sadourny enstrophy / energy `sadourny_energy` / Arakawa-Hsu `sadourny_hk`) | `coriolis_adv_t` | `kernels/coriolis_adv/rdb_coriolis_adv.F90` | 3 | per-layer u, v, layer thickness | momentum tendency at faces |
| FV pressure force | `ocean_pressure_force_t` | `../../pressure_force/rdb_ocean_pressure_force.F90` | 5d | `multilayer.h_layer`, T, S, `eos` — including `eos%rho0`, which `configure_ocean_pgf` copies into BOTH slot reference densities (see the **PGF reference densities** contract below); `multilayer.p_top` when `&ocean_pgf_nml p_top_in_bc` (FV_MOM6 only, default off ⇒ bit-identical — see the **`p_top` seam contract** below) | momentum tendency at faces; `e_face` for the barotropic `compute_pbce` |
| Surface momentum stress ✓ | `ocean_surface_stress_t` | `../../parameterizations/vertical/rdb_ocean_surface_stress.F90` | 5b | `tau_x`/`tau_y` (wind, or ice-blended via `rdb_ice_ocean_coupler`); `rho0` from `eos%rho0` via `configure_ocean_reference_density` | momentum tendency at `k=nz`; `stress_mag` (cell-centred `\|tau\|`, always allocated, refreshed by EVERY writer of the `tau` pair — the `set_wind_stress_*` setters at configure, the data-forcing seam refresh, and the sea-ice blend on device each outer step, all through `ocean_surface_stress_refresh_mag` — PR-12 dedup, read by both KPP and EPBL instead of each re-deriving it inline, so a `tau` write that skips the refresh freezes both schemes' `u_*`). **Ice-shelf cover:** `ocean_surface_stress_apply_cover(ss, cover_frac)` zeroes the `tau` PAIR on every face touching a covered cell (EITHER-neighbour rule, `1 - max(cover_L, cover_R)`) and refreshes `stress_mag` in the same call -- the mask is applied to the SOURCE, not to each derived view, because `tau` also reaches the implicit vdiff stress fold (`tau_u=ss%tau_x`) and the MLE front sampler raw. Applied once from `configure_ocean_cavity` (after the cover is built, before `enter_data`) and again through `ocean_surface_stress_set_derived`'s optional `cover_frac` at the data-forcing seam; idempotent; cavity off => byte-identical.  **Second stress, Phase 4b:** `stress_shelf` (cell-centred `\|tau_top\|` at an ICE-SHELF base, N/m^2, ALWAYS allocated + mapped + counted, exactly zero without a cavity) is the OTHER half of the upper-boundary momentum flux — see the **`stress_mag` / `stress_shelf` contract** below.
| Ice-shelf TOP drag ✓ (Phase 4a; default off) | `ocean_top_drag_t` | `../../parameterizations/vertical/rdb_ocean_top_drag.F90` | 4a | `multilayer.u_face_x_layer`/`v_face_y_layer`, `h_layer`, `wet_mask`; its own configure-filled FACE cover masks `cover_u`/`cover_v` (the **OR** of the two abutting cells' `metrics.cover_frac`, `max(cover(i-1,j), cover(i,j))`, so the CALVING-FRONT face is dragged) — the SAME face rule `ocean_surface_stress_apply_cover` uses to zero the wind there, enforced face-for-face by `cavity_cover_face_rule_matches_top_drag` rather than left as two copies of one sentence and the cell-centred copy `cover_t`; `rho0` from `eos%rho0` via `configure_ocean_reference_density` (the `stress_top` diagnostic only) | momentum tendency at `k = nz` (`du_drag`/`dv_drag`, layer-only or spread over `htbl` metres), the `k = nz` Rayleigh rate `lambda_top_u/v` for the vdiff fold, and `stress_top` (cell-centred `\|tau_top\|`, N/m², device-resident — copied INLINE by `run_stage`/`run_stage_split` into `surface_stress%stress_shelf`, whence KPP and EPBL take their under-ice `u_*`; see the **`stress_mag` / `stress_shelf` contract** below). Requires `&ocean_cavity_dyn_nml enable` (fail-loud). The tendency is ALSO added into `bt_work%F_slow_u/v` by `add_top_drag_into_F_slow` — a layer tendency left out of that sum is invisible to the barotropic substep AND mis-corrected by `apply_bt_correction`. `&ocean_cavity_melt_nml cdrag_top` must equal `&ocean_tdrag_nml cd`: ONE ice-base drag coefficient, refused at configure if they differ. `enable=.false.` (default) ⇒ placeholder arrays, no kernel, byte-identical |
| Surface heat/salt flux ✓ (PR-12 component-set reshape) | `ocean_surface_flux_t` | `../../parameterizations/vertical/rdb_ocean_surface_flux.F90` | 5b | const scalars (`&ocean_thermo_nml q_heat/q_salt`) +, when `&ocean_forcing_nml enable_components`, the component set (`q_sw/q_lw/q_lat/q_sens/heat_added`, mass fluxes `evap/lprec/fprec/vprec/lrunoff/frunoff/seaice_melt`, their `heat_content_*` enthalpy companions, `salt_flux`, `p_surf_atm`); `rho0` from `eos%rho0` via `configure_ocean_reference_density` | `Q_heat`/`Q_salt` — **derived views**, always the fields every downstream kernel (KPP, EPBL, `apply_tracers`) reads. Components off (default): `Q_heat`/`Q_salt` = the const scalar fill, byte-identical to pre-PR-12. Components on: `ocean_surface_flux_assemble` (the single gate, `vmix_assemble`'s analogue) rebuilds them every thermo step from const + components (`heat_content_massin`/`massout` also assembler-owned outputs) — a filler writes ITS OWN component and MUST NOT write `Q_heat`/`Q_salt` directly, must set `has_heat`/`has_salt` (+ `has_mass_flux`/`has_q_sw` as it fills mass/`q_sw`) host-side, and must register a time-varying component itself in the restart registry (`Q_heat`/`Q_salt` themselves are never registered — derived-field rule). **Component ownership is one filler per component** — the table below the slot map spells it out. The sea-ice coupler (`rdb_ice_ocean_coupler`) writes `salt_flux`/`heat_added` when components are on, the legacy full-overwrite of `Q_salt`/`Q_heat` when off; the ice-shelf cavity (`rdb_ocean_cavity_flux`) writes `heat_cavity`/`salt_cavity` and NOTHING else, precisely because the ice coupler full-overwrites the two it does not touch. `p_surf`/`p_surf_atm` ship zeroed with no consumer yet (PR-17 follow-up). **Ice-shelf cover:** `ocean_surface_flux_assemble`'s optional `cover_frac` multiplies every ATMOSPHERIC contribution (`Q_heat_const`/`Q_salt_const`, `q_sw`/`q_lw`/`q_lat`/`q_sens`/`heat_added`, both `heat_content_mass*`, `salt_flux`) by `1 - cover_frac` and leaves `heat_cavity`/`salt_cavity` alone. The mask lives HERE and not at apply time for two reasons: `Q_heat`/`Q_salt` are what KPP and EPBL read to build `B_0`, so masking later would force both boundary-layer schemes with an atmosphere that is not there; and this is the last point at which the atmospheric bands are still separable from the cavity bands. Components OFF => no assembler, so the static scalar fill is masked once at configure by `ocean_surface_flux_apply_cover_const`. `ocean_surface_flux_apply_sw_penetration` and `ocean_surface_restore_apply_tracers` carry their own optional `cover_frac` (neither routes through `Q_heat`/`Q_salt`). Absent => the original kernel, byte-identical |
| Vertical mixing (KPP) | `ocean_vmix_t` | `../../parameterizations/vertical/rdb_ocean_vmix.F90` | 5b | u, v, T, S, surface forcing; `rho0` from `eos%rho0` via `configure_ocean_reference_density` (the ONE of these copies a kernel reads on-device) | `kv`, `kt`, `ks` on interfaces; non-local `gamma_t/s`. `ks` is DERIVED from `kt` by `vmix_split_kd_heat_salt` (last statement before `vmix_assemble`, PR-20; `ks ≡ kt` until a double-diffusion contributor lands) and consumed by `vdiff_apply_tracers` for salinity + every passive tracer. **α/β source (E4, `&ocean_vmix_nml buoyancy_coeffs`)**: the KPP `B_0` (both passes) and the double-diffusion density ratio historically read the CONSTANT `eos%alpha_T`/`eos%beta_S` whatever the active EOS. `"eos"` swaps in `eos_buoyancy_coeffs` from the ACTIVE EOS — per column for `B_0`, at `ms%p_top` under `&ocean_psurf_nml in_eos` else `eos%p_ref` (a SURFACE flux has no depth of its own); per interface for the density ratio, at the true in-situ pressure seeded from `p_top` and accumulated `g·ρ₀·h` downward, the shape EPBL's stack already has. `"constant"` is the default ⇒ bit-identical, and under `eos="linear"` the two are byte-identical by construction. `p_top_in_eos` mirrors EPBL's `in_eos` gate — NOT redundant with `p_top` being zero, since a cavity fills `p_top` either way. The constant ddiff impl keeps its fully-collapsed `(i,j,k)` launch untouched; the EOS twin is a separate column-structured kernel because the pressure is a running sum down the column |
| EPBL (energetics PBL) ✓ | `ocean_epbl_t` | `../../parameterizations/vertical/rdb_ocean_epbl.F90` | 5b+ | T, S, h, wind stress, surface flux, abs(f) | `kd_int` on interfaces (merged into `vmix%kv/kt` each stage), `mld`, TKE-budget diags |
| kappa-shear (JHL08 interior) ✓ | `ocean_kappa_shear_t` | `../../parameterizations/vertical/rdb_ocean_kappa_shear.F90` | 5b+ | u, v (face→centre), T, S, h, f²; vertex mode (`at_vertex`): native face u/v + metrics wet masks + corner f² | `kd_int`, `tke_int` on interfaces (additive merge into `vmix%kv/kt` each stage; column solve at thermo cadence).  Vertex mode adds `f_corner` + the `kd_corner` carrier (allocated by `init_vertex` at CONFIGURE, mapped by the slot's `enter_data`); corner solve → scatter are two kernels, never fused (`kd_corner` is the race-breaking snapshot); `tke_int` zeroed; Kv goes corner→face — `kd_corner`·`prandtl_turb` feeds `vdiff_apply_momentum kv_corner_source` (cell-centred kv merge suppressed) |
| tidal mixing (St-Laurent/Simmons internal-tide) ✓ | `ocean_tidal_mixing_t` | `../../parameterizations/vertical/rdb_ocean_tidal_mixing.F90` | C (Area-A) | `e_in` (prescribed bottom energy), `h_layer`, T, S, `eos` (N²), `wet_mask` | `kd_int` on interfaces (bottom-intensified exp decay from bed; additive merge into `vmix%kv/kt` each stage; column flux-bookkeeping sweep at thermo cadence). Default off ⇒ bit-identical |
| Wave speed + Rossby radius ✓ (B1; diagnostic, default off) | `ocean_wave_speed_t` | `../../parameterizations/vertical/rdb_ocean_wave_speed.F90` | B1 | `multilayer.rho_layer`, `h_layer`, `wet_mask`, `f_centre`/`beta_centre`, `metrics.dxT` | `cg1`, `rd`, `rd_over_dx` (nx,ny); per-column Sturm–Liouville eigensolve (backtracking merge + fixed-budget bisection); reads `rho_layer` directly (no outer-shim); called once per outer step at thermo cadence (`n_wavespeed`-gated) in `ocean_dyn_step_split`, before both `varmix_compute` sites and `run_meke_step` — feeds B2. `rd_over_dx = rd/metrics.dxT` (metres — NOT `grid%dx`, which is degrees on spherical/supergrid/tripolar). `f_centre`/`beta_centre` (`beta_centre = |∇f_centre|`, the `meke_length_scales` stencil) come from `metrics_fill_coriolis` (planetary dispatch on spherical/tripolar; bit-identical beta-plane elsewhere) via `build_static`, not a hard-coded beta-plane. Unsplit `ocean_dyn_step` does not call it (trap #5, PR-3) |
| Lateral mixing (Leith) ✓ | `ocean_lateral_mix_t` | `../../parameterizations/lateral/rdb_ocean_lateral_mix.F90` | 5a | u, v at faces | `ah_face_*` (Leith populated; Smag/biharmonic deferred) |
| Isopycnal slopes ✓ (diagnostic, default off; mesoscale gate) | `ocean_slopes_t` | `../../parameterizations/lateral/rdb_ocean_isopycnal_slopes.F90` | 5a | `multilayer.h_layer`, T, S, `eos`, metrics (`idxCu`/`idyCv`/`wet_u`/`wet_v`) | `slope_x`/`slope_y` + `n2_u`/`n2_v` at interfaces (neutral slope `S=−∇ρ/∂_zρ`, locally-referenced ρ; `vert_fill_TS` regularizer; bed/surface forced 0). Thermo cadence; `&ocean_slopes_nml enable` default off ⇒ bit-identical. Feeds future GM/Redi/VarMix |
| GM thickness diffusion ✓ ([2] mesoscale; default off) | `ocean_gm_t` | `../../parameterizations/lateral/rdb_ocean_gm.F90` | 5a | `slopes.slope_x`/`slope_y` + `n2_u`/`n2_v`, `h_layer`, metrics | per-face `khth_u`/`khth_v` (constant, CFL-clamped; VarMix/MEKE additive seam) + bolus transports `uhD`/`vhD` folded into `mass_flux_{x,y}_layer` (MOM6 uhtot recurrence + safe-streamfunction + mass-availability limiter, `Σ_k uhD=0` conservative); `gm_src` S²N²κ PE-release (feeds MEKE). Thermo cadence; `&ocean_gm_nml enable` default off ⇒ bit-identical; requires `slopes` on (fail-loud at configure) |
| VarMix coefficients ✓ ([4] mesoscale; default off) | `ocean_varmix_t` | `../../parameterizations/lateral/rdb_ocean_varmix.F90` | 5a | `wavespeed.cg1`, `f_centre`, metrics, `slopes.slope_x/y` + `n2_u/v` | 2D base `khth_u/v` + `khtr_u/v` = `clamp((KHTH + cff·L²·SN)·Res_fn)`: Hallberg-2013 resolution fn (divide-free, `cg1`-based) + Visbeck/Eady `SN` (thickness-weighted, orthogonal slope in S²) + `Res_fn`-before-clamp. `khth_u/v` feed GM's `khth_ext`; `khtr_u/v` feed Redi's per-face KhTr (both consumed). Thermo cadence; `&ocean_varmix_nml enable` default off ⇒ bit-identical; requires `slopes` + `wavespeed` on (fail-loud) |
| Redi neutral diffusion ✓ ([3] mesoscale; default off) | `ocean_redi_t` | `../../parameterizations/lateral/rdb_ocean_redi.F90` | 5b+ | T, S, `eos`, `h_layer`, metrics (rebuilds its OWN interface dRdT/dRdS — does NOT consume the GM slopes) | along-isopycnal (neutral) tracer flux, MOM6 CONTINUOUS variant, two-phase: Phase A `calc_coeffs` (the `2nz+2` neutral-surface sweep → device-resident `uPoL/uKoL/uhEff` + v-mirror; `hEff` converted Pa→m via `H_to_pa=g·ρ₀`) + Phase B per-tracer `neutral_surface_flux` (PPM sublayer Δ, down-gradient sign guard, flux-form scatter with the bottom-up k-flip `knat=nz+1-Ko`). AUGMENTS `hdiff_tracer` (keeps a small along-coordinate Kh). Thermo cadence; `&ocean_redi_nml enable` default off ⇒ bit-identical. Needs `NZ_STACK_MAX≥2·nz` (raised to 256; verified launches at nz=75). R2 flux shipped. VarMix Visbeck KhTr seam WIRED: when VarMix is enabled its per-face `khtr_u/v` override the scalar `khtr` (broadcast to faces otherwise ⇒ scalar path bit-identical). Discontinuous variant + interior_only BL clamp deferred |
| MEKE prognostic eddy energy ✓ ([5] mesoscale; default off) | `ocean_meke_t` | `../../parameterizations/lateral/rdb_ocean_meke.F90` | 5a | `gm.gm_src` (PE release, prev thermo step), `varmix.sn_u/v` (Eady), `wavespeed.rd_over_dx` (deform), `h_layer`/`rho_layer`, metrics | 2D prognostic `meke(nx,ny)` — Strang split (explicit GM/bg source → implicit backward-Euler bottom drag half → harmonic-mass Laplacian + biharmonic → drag half); derived `kh = khcoeff·√(2·γt²·E)·Lmix` (harmonic length-scale sum) fed as `khth_fac·√(kh_i·kh_{i+1})` (geom-mean) into `varmix.khth_u/v` (+ khtr) BEFORE GM's CFL clamp — closes the GM↔MEKE loop (one-step lag). Restart-persistent. Thermo cadence; `&ocean_meke_nml enable` default off ⇒ bit-identical; `khth_fac=khtr_fac=0` ⇒ feedback inert; requires `gm` on (fail-loud). Rhines length live (`alpha_rhines>0`; `beta=|∇f|` from the slot's own `f_centre`, filled at setup from the Coriolis path, scaled by `idxT`/`idyT`); upwind barotropic advection live (`advection_factor>0`; mass-weighted `baroHu` from `mass_flux_*_layer`, conservative — `Σ E·area·mass` telescopes). BBL drag live (`&ocean_meke_nml use_bbl_drag`, default off ⇒ bit-identical): the resolved bed-layer (k=1) eddy speed `|u_bed|²` adds into the bottom-drag rate `ρ₀·i_mass·√(cdrag²·(2·γb²·E + |u_bed|² + uscale²))` (the ρ₀ factor is MOM6's `GV%H_to_RZ`; `ρ₀·i_mass ≈ 1/H` ⇒ rate in 1/s; MOM6 `drag_rate_visc`). Frictional source live (`&ocean_meke_nml frcoeff>=0`): the hvisc slot exposes its KE-dissipation rate `ke_diss = Σ_k ρ_k h_k (u·du_visc + v·dv_visc)` (captured post-compute/pre-apply, gated by `hvisc%compute_ke_diss` so default runs do zero extra work), and MEKE adds `-frcoeff·i_mass·ke_diss` (mean→eddy energy). Default `frcoeff<0` ⇒ inert ⇒ bit-identical. Harmonic backscatter live (`&ocean_meke_nml backscatter` + `backscatter_visc_coeff_ku`, default off ⇒ bit-identical): `meke_step` fills `ku = coeff·√(2·E)·Lmix` (MOM6 `MEKE_VISCOSITY_COEFF_KU`; plain `√(2·E)`, no `γt²` unlike kh; harmonic only, `BS_struct=1`) and `meke_backscatter_apply` subtracts a face-average of `ku` from `lateral_mix%ah_face_x/y` (in `run_stage_split` after `ocean_lateral_mix_compute`, one-step lag) so the NET (resolved − ku) harmonic viscosity goes negative — the momentum energy return; floored at the forward-Euler viscous-CFL lower bound `−0.8·0.5/(dt·(idx²+idy²))` per face (MOM6 `BACKSCATTER_UNDERBOUND`), which bounds the negative mode's growth rate but does not stabilise it alone — a positive biharmonic backstop with a NON-ZERO dissipation coefficient (`nu_4>0` with no flow-aware closure selected, or `smag_ah`+`smag_bi_const>0`, or `leith_biharm`+`c_leith_bi>0`, or a positive `nu_4_bg` floor under either flow-aware closure — selecting a flow-aware closure at its default zero coefficient does NOT satisfy the guard) is mandatory, enforced fail-loud at configure (`has_biharmonic_backstop`); also needs a flow-aware closure active (else `ah_face_*` unused). The biharmonic add-on composes with `stress_tensor` — it does not bypass it. Deferred: biharmonic `Au`, EBT/SQG `BS_struct`. All v1 upstream seams now wired |
| Fox-Kemper ML-eddy restratification ✓ (B5) | `ocean_mle_t` | `../../parameterizations/lateral/rdb_ocean_mle.F90` | 5b+ | `epbl%mld`, `rho_layer`, `h_layer`, metrics (`dy_cu`/`dx_cv`/`f_centre`), surface stress (mixrate form) | `b_ml`/`htot_ml` (2D), per-layer `uhml`/`vhml` folded into `mass_flux_{x,y}_layer` before the continuity divergence (conservative, `Σ_k a(k)=0`; never touches velocities). Thermo cadence; default off ⇒ bit-identical |
| Equation of state (Wright) | `eos_t` | `../../equation_of_state/rdb_eos.F90` | 5d | T, S, p | rho (in-place, by FV-PGF / vmix / diag). Point routines, all `!$acc routine seq` on the flat-POD handle passed BY VALUE: `eos_density_point` (ρ), `eos_specvol_derivs` (dSV/dT, dSV/dS — for consumers that genuinely work in specific volume: EPBL's PE weights, kappa-shear's `g·ρ₀·dSV`), `eos_buoyancy_coeffs` (**E4** — `α = −∂ρ/∂T`, `β = +∂ρ/∂S`, the DIMENSIONAL kg/m³-per-unit convention the `eos_t` members carry, analytic per variant; the LINEAR branch returns those members BIT-FOR-BIT rather than `−ρ²·dSV/dX`, which is what makes `&ocean_vmix_nml buoyancy_coeffs` byte-identical under `eos="linear"`), `eos_density_derivs` (its exact sign twin) and `eos_freezing_point`. The last two are `pure elemental` as well, so a host sweep can evaluate an array in one reference |
| Vertical coordinate + ALE remap ✓ | `ocean_vcoord_t` | `vcoord/rdb_ocean_vcoord.F90` | 5g | `multilayer.h_layer`, `barotropic` H (bathy), `dyn.bt_eta` | `target_h(i,j,k)` per outer step; `rdb_ocean_remap` advances `h_layer` + `tracers(t)%hTr` + face velocities via PPM column remap, wired into `ocean_dyn_step_split_multilayer` as optional `vcoord=`. **`z_top(i,j)` — the rigid-top seam (P6.2).** The slot owns one more static 2-D array: the geopotential depth of the top of the WATER column (m, positive down), `metrics%z_draft` under an ice-shelf cavity and `0` everywhere else. Filled ONCE at configure by `configure_ocean_cavity` (the draft is static), which is why the remap driver's signature is unchanged — it never sees `metrics`. Allocated UNCONDITIONALLY at `(nx_total, ny_total)` with `source = 0.0_wp`, `copyin` in `enter_data`, counted in `bytes()`, freed in `destroy`: unlike `metrics%z_draft` (a `(1,1)` placeholder without a cavity) it is ALWAYS safe to hand to an explicit-shape device dummy. Consumed by `VCOORD_Z_FIXED` alone today (`ocean_vcoord_z_fixed_target`, shared with the IC seed), which measures its nominal interface depths from `z = 0`, vanishes the layers that outcrop into the ice to `zstar_h_min`, and cuts the first live layer at the ice base into a partial top cell (minimum `0.1*h_nominal`, else the sliver merges into the layer below — the mirror of the bed's own sliver rule). `z_top ≡ 0` reproduces the pre-cavity arithmetic BIT-for-bit and that is asserted with `==` in `tests/test_ocean_vcoord_zfixed_cavity.F90`. **v1 envelope: `z_fixed` × cavity is ADIABATIC DYNAMICS ONLY** — on a covered column `k = nz` is an inert filler, so `validate_config` fails loud on melt, top drag, KPP, EPBL, ideal age, GM/Redi/slopes, `kappa_h /= 0` and `regrid_time_scale > 0`, and requires `h_min_cavity >= 2*h_nominal`; atmospheric forcing is NOT refused because the cover mask already zeroes it under ice and an open column has `z_top = 0`. A draft that VARIES is not refused but WARNED (the staircase corrections are a later slice, and those slices need it runnable); only a uniform draft is validated, and it is bit-zero. The shared `k_top(i,j)` that lifts those refusals is the next slice. |
| Tides | `ocean_tides_t` | `forcing/rdb_ocean_tides.F90` | 5e | astronomical clock + bathymetry | `eta_eq`, `eta_sal` |
| Sea ice ✓ (SIS2 port; default off; largest slot in the ocean state — 14 modules, ~6,900 lines) | `ocean_sea_ice_t` | `../ice/state/rdb_ice_state.F90` | (ice ladder PR 0-5) | ocean surface T/S/u/v (`sst_seam`/`ssurf_seam`/`tfw_seam`), wind stress (`tau_a_x/y`), `eos_freezing_point` | 6 category prognostics (`part_size`, `m_ice`, `m_snow`, `enth_ice`, `enth_snow`, `sal_ice`) + Winton column (`rdb_ice_column`) + ITD restore + category transport/`compress_ice` (`transport`) + C-grid EVP (`dynamics`: `u_ice`/`v_ice`/`str_d`/`str_t`/`str_s`/`fxoc`/`fyoc`) + frazil bank (`frazil_heat`); `rdb_ice_ocean_coupler` writes ocean `Q_heat`/`Q_salt`. **A high-quality dynamical core + column model, not yet a sea-ice model** — no ridging/snowfall/flooding/melt-ponds/lateral-melt/IC-path, transmitted shortwave `sw_thru` now coupled to the ocean (PR 31 — `ice_ocean_sw_flux`), momentum not conserved at fractional cover, no freshwater/mass coupling; fail-loud single-rank (`enable` alone, `transport`, `dynamics` each independently guarded, `rdb_config.F90`). Full limits list: `docs/CAPABILITIES_AND_LIMITATIONS.md` "Sea ice"; matrix rows: `docs/CLOSURE_MATRIX.md` "Sea ice" |
| Tides | `ocean_tides_t` | `forcing/rdb_ocean_tides.F90` | 5e | astronomical clock + bathymetry | `eta_eq`, `eta_sal`, `itd_coeff` |
| Surface-pressure loading / inverse barometer ✓ (PR-17, Wunsch & Stammer 1997; `&ocean_psurf_nml enable`, default off ⇒ bit-identical) | `ocean_p_surf_t` | `forcing/rdb_ocean_p_surf.F90` | 5e | `sf%p_surf` (PR-12 component set, needs `enable_components`), `eos%rho0`, `dyn%bt_work%g_bt` | `eta_ib = −p_surf/(ρ₀·g_bt)` + combined `eta_seam`, filled by `p_surf_update_seam` once per outer step and passed as the barotropic `eta_forcing` (see the **`eta_forcing` seam contract** below). Split-solver only; excludes `bt_halo > 0`. With `&ocean_psurf_nml in_eos` (E3, default off ⇒ bit-identical) the SAME `sf%p_surf` is also copied once per outer step into `multilayer_state_t%p_top` (Pa) so the EOS's IN-SITU pressure is measured down from the load rather than from 0 Pa — NOT the potential density `ms%rho_layer`, which stays at the uniform `eos%p_ref`; see the **`p_top` seam contract** below. |
| Porous barriers ✓ (Adcroft 2013; default off) | fields on `ocean_metrics_t` (`use_porous`, `porous_eta_interp`, `porous_mask_depth`, `por_bed`, `por_{dmin,dmax,davg}_{u,v}`, `por_face_area_{u,v}`) | kernels in `state/rdb_ocean_porous.F90`; configure in `state/rdb_ocean_setup.F90::configure_ocean_porous`; per-step refresh `ocean_porous_refresh` in `dynamics/split_rk2/rdb_ocean_dyn.F90` | — | `barotropic.b` (negated once to a topographic HEIGHT), `multilayer.h_layer` | `por_face_area_u/v` (nx+1,ny,nz)/(nx,ny+1,nz) — layer-averaged OPEN-AREA fraction, recomputed once per OUTER step (MOM6 cadence) and MULTIPLIED into `mass_flux_{x,y}_layer` by continuity-PPM (before the BT renormalisation, which also takes the narrowed areas) and into `coriolis_adv.mass_flux_{u,v}` by the TRANSPORT Coriolis forms — PLUS `dy_cu_bt`/`dx_cv_bt` (2D, ALWAYS allocated, byte-equal to `dy_cu`/`dx_cv` when off), the widths the BAROTROPIC substep transports on, scaled by the COLUMN-INTEGRATED open fraction so the barotropic solve is not porous-blind (else the renormalisation to `uhbt` returns the blocked transport). `use_porous=.false.` (default) ⇒ the open-area fields stay at their `(1,1,1)` placeholder, no porous kernel is launched, byte-identical. Fails loud with `&ocean_bt_nml bt_halo > 0` (`bt_wide`'s own `metrics_w` carries no porous stats, so the wide BT loop would silently transport on un-narrowed widths) and with `&ocean_wetdry_nml enable` |
| Ice-shelf cavity statics ✓ (P5.1 geometry + datum, P5.2 load; default off) | fields on `ocean_metrics_t` (`use_cavity`, `z_draft`, `cover_frac`, `p_ice_ref`) | geometry + helpers in `state/rdb_ocean_cavity.F90`; draft fill + grounding in `state/rdb_ocean_state.F90::seed_cavity_draft` (inside the IC seed, before the wet mask); load build + `ms%p_top` assembly + datum assertion in `state/rdb_ocean_setup.F90::configure_ocean_cavity`; per-step re-assembly (psurf seam live only) inline in `dynamics/split_rk2/rdb_ocean_dyn.F90::ocean_dyn_step_split` | — | `&ocean_cavity_dyn_nml`, `barotropic.b`, `pressure_force.rho_ref`, `surface_flux.p_surf` | `z_draft` (m, positive down, ghosts included) — consumed by `configure_ocean_bt_split` as the DATUM `bt_H_ref = b − z_draft` afloat and `0` where GROUNDED (`cavity_datum_impl`), by the layer/`bt_h` seed as the water column, and by `seed_wet_mask_impl` as the grounding decision; `cover_frac` (binary 0/1) — the solve mask the basal-melt slot composes with the wet mask; `p_ice_ref = (rho_ref*GRAVITY)*z_draft` (Pa) — the static half of `multilayer.p_top = p_ice_ref + sf%p_surf`, and thence the FV_MOM6 `pa(nz+1)` top BC (`&ocean_pgf_nml p_top_in_bc`, REQUIRED for a non-uniform draft) and the in-situ EOS (`&ocean_psurf_nml in_eos`). It is deliberately NOT a component of `sf%p_surf`: the datum already carries its barotropic effect and `eta_ib` is built from that total. `use_cavity=.false.` (default) ⇒ all three stay `(1,1)` placeholders and `bt_H_ref = b` byte-identically. See the **cavity datum contract** below |
| Ice-shelf basal melt ✓ (P2b; default off) | `ocean_cavity_flux_t` | `../../parameterizations/vertical/rdb_ocean_cavity_flux.F90` (kernel: `rdb_ocean_cavity_melt.F90`) | 2b | `metrics.cover_frac`, `multilayer.h_layer` + S/T + face velocities + `wet_mask`, `multilayer.p_top` (THE interface pressure), the shared `eos` handle (the liquidus, `&ocean_eos_nml tfreeze_set="isomip"`), and its own configure-filled `f_cor` | the two OWNED surface-flux components `surface_flux.heat_cavity` (= −`q_ocean`, W/m² positive down ⇒ warm water COOLS) and `surface_flux.salt_cavity` (= −`m_mass·(S_far − s_ice)`, a VIRTUAL salt flux ⇒ melting FRESHENS); plus its own 2-D interface state (`t_far`/`s_far`/`u_far`/`v_far`/`ustar`/`t_b`/`s_b`/`melt`/`q_ocean`/`gamma_t`/`gamma_s`/`active`/`status`). `gamma_t`/`gamma_s` are the exchange VELOCITIES the solve converged on (filled by `cavity_melt_point_gamma`, which is the same `cavity_solve_melt` call, not a second solve) — stored rather than re-derived because under `hj99`/`yung25` they are implicit in the interface state, so an `exch_vel_*` diagnostic rebuilt from `u*` would report the NEUTRAL values. Thirteen derived-diag catalog entries read this slot (`melt`, `melt_m_per_yr`, `thermal_driving`, `haline_driving`, `tbdry`, `sbdry`, `tfreeze_ib`, `exch_vel_t`, `exch_vel_s`, `ustar_shelf`, `cavity_melt_status`, plus the geometry pair `z_draft`/`water_column`), NaN outside the cover and fail-loud at configure without their prerequisite knob. Far field sampled over `far_field_depth` METRES below the ice base, thickness-weighted with a partial last layer — never "layer nz". Driven once per thermo step from `engine_step_finalize`, immediately BEFORE `ocean_surface_flux_assemble`. The `do concurrent` over columns lives in the KERNEL module (`cavity_melt_columns_2d`), not here: nvlink cannot resolve an `!$acc routine seq` device symbol out of `librdb_core.so` into a `do concurrent` in another translation unit. Budgets ride the ordinary `heat_budget_surface`/`salt_budget_surface` contributors, so no new accumulator and no stage-weight decision. `enable=.false.` (default) ⇒ fourteen `(1,1)` placeholders, no kernel, byte-identical |
| Barotropic linear wave drag ✓ (Egbert & Ray 2001; Jayne & St Laurent 2001) | fields on `barotropic_workstate_t` (`dyn.bt_work`) | kernels in `kernels/barotropic/rdb_barotropic_coupling.F90`; configure in `state/rdb_ocean_setup.F90::configure_ocean_wave_drag` | — | `lwd_drag_u/v` (static, host-filled at configure from `form="uniform"` or `"roughness_proxy"`; `barotropic.b`, `metrics.wet_T` for the proxy) | MULTIPLIES into `bt_work.bt_rem_u/v` each stage (`compute_bt_rem_wave_drag`); `lwd_enable=.false.` (default) ⇒ arrays unallocated, bit-identical |
| River / discharge | `ocean_river_t` | `forcing/rdb_ocean_river.F90` | 5e | sources NetCDF | `q_mass`, `q_S`, `q_T` distributed fields |
| Open boundary (parent nest) | `ocean_obc_t` | `boundary/rdb_ocean_obc.F90` | 5c | parent NetCDF | ring buffers + FRS blend into edge bands |
| Periodic wrap ✓ | (free procedures, `rdb_ocean_periodic`) | `boundary/rdb_ocean_periodic.F90` | 5c | `bc%periodic_x/y` + prognostics | ghost cells = opposite-interior copies (seam invariant; per-substep inline wraps live in `rdb_barotropic_substep`) |
| Tripolar north fold ✓ | (free procedures, `rdb_ocean_fold` + `rdb_ocean_fold_apply`) | `boundary/rdb_ocean_fold{,_apply}.F90` | M4c | `bc%north_fold` + prognostics + metrics + `f_corner` | reversed-i halo-fill (sign-flip vector normals) + on-row antisymmetric v/corner projection; fires AFTER each periodic wrap (state, BT inline, continuity mid-split, driver init); metric/`f_corner` ghosts folded once at configure |
| Baroclinic OBC ✓ | (free procedures, `rdb_ocean_obc_baroclinic`) | `boundary/rdb_ocean_obc_baroclinic.F90` | 5c | `bc` tags, `bt_ubt_end`, layer state | wall-face per-layer u/v (Flather mean + zero-gradient anomaly), open-edge h/hTr ghost fills, + zero-gradient ghost fills of the per-layer velocities, tracer concentration, and barotropic η **including the ghost×ghost corners**, for open×open AND wall×open corners (only periodic edges keep their wrap). The corner-adjacent Coriolis/KE stencil reads the diagonal corner ghost, so unfilled corners blow up a multilayer domain — at-rest (velocity/η) and under inflow (tracer) |
| Real sponge ✓ (PR-23; default off) | `ocean_sponge_t` | `boundary/rdb_ocean_sponge.F90` | 5c | `bc%<edge>%bc_type==OBC_SPONGE` tags (`damp_source="band"`), seeded IC (`ocean_sponge_snapshot_reference`, `target_source="ic"`) | per-cell `idamp_h/u/v` (1/s) + 3-D `ref_tracer`/`u_ref`/`v_ref`; `ocean_sponge_apply_maps` relaxes momentum toward `u_ref`/`v_ref` (not zero) + every registered tracer toward `ref_tracer` (not a scalar), mirroring S/T into `ms%salt_budget_sponge`/`heat_budget_sponge`. `&ocean_sponge_nml enable=.false.` (default) ⇒ the legacy `bc`-band sponge (`ocean_sponge_apply`/`_tracers`) runs unchanged ⇒ bit-identical. `relax_h` + `damp_source`/`target_source="file"` deferred to PR-23b |
| Diagnostics manager ✓ (registry + remap + time-mean + serial NetCDF) | `ocean_diag_t` | `diag/rdb_ocean_diag.F90` (+ `rdb_ocean_diag_fills.F90`, `rdb_ocean_diag_netcdf.F90`) | 6a–6d | every prognostic + derived field via procedure-pointer fills | per-rank NetCDF file (serial; I/O server hand-off pending MPI) |
| Time-varying NetCDF input reader ✓ | `ocean_data_input_t` | `io/rdb_ocean_data_input.F90` | PR-14 | consumer-registered `(file, var, dest)` triples | opaque `id` from `ocean_data_input_register_2d/_3d/_segment_2d/_segment_3d`; `update_2d/_3d` blends linearly-in-time on-device into the caller's whole, mapped array at a registration-time offset; `fill_static_host{,_3d}` for setup-time `DATA_TIME_STATIC` fills. Target-agnostic — owns no field names/state-slot mapping; `&ocean_data_nml` carries only `max_fields`/`verbose`, each consumer registers from its own namelist group. Pre-regridded files only (no in-core horizontal interp, no vertical remap, no calendar); `t_offset` is the whole time-alignment. Zero registrations ⇒ bit-identical. |
| Restart manager ✓ | `ocean_restart_t` + `restart_registry_t` | `io/rdb_ocean_restart{,_io}.F90` | shipped | every prognostic + per-slot transient (registry walk) | per-rank restart NetCDF (`restart_rank_NNNNNN.nc`) |
| Z-level T/S IC ✓ | (free procedures, `rdb_ocean_z_init` module — no new `_t`; config on `cfg%ocean%zinit`) | `io/rdb_ocean_z_init.F90` | A2 | pre-regridded model-grid T/S NetCDF (`source="file"`) or the analytic affine `lin_*` profile (`source="linear"`), plus seeded `h_layer` + `wet_mask` + (under a cavity) `metrics%z_draft` | overwrites `tracers(idx_S/T)%hTr` at seed time (host-side, before enter_data). Layer-centre depth is GEOPOTENTIAL — measured from `z = 0`, so the column-top offset `z_draft` is added under an ice shelf (`build_z_ctr(h, nz, z_top, z_ctr)`; `z_top = 0` is bit-identical to the pre-cavity arithmetic). File source: linear-in-depth interp + constant tails, interior only, dry columns → namelist land-fill. Linear source: exact, FULL array incl. ghosts + land. NetCDF-gated. |
| 3D scratch buffer ✓ | `scratch_3d_buffer_t` | `../../framework/rdb_scratch_3d.F90` | 1 (then ongoing) | `(n1, n2, n3)` shape from caller | `(stride, stride, nz)` scratch storage with bound init/destroy/enter_data/exit_data; future `ensure_size` |
| Safe-math wrappers | (free procedures, `rdb_safe_math` module) | `../../framework/rdb_safe_math.F90` | library, no production consumer (PR-8 removed `RDB_BITWISE_REPRO`, which was the sole would-be caller) | n/a | plain elemental inlines around the intrinsic (`safe_exp/log/sin/cos/sqrt/pow`); the tested `*_polynomial` implementations are raw material for a future repro PR |

### The `stress_mag` / `stress_shelf` contract — two upper-boundary stresses

`ocean_surface_stress_t` carries TWO cell-centred stress magnitudes and they
are different things.  Both boundary-layer schemes — KPP (`rdb_ocean_vmix`,
two sites) and EPBL (`rdb_ocean_epbl`, one) — read

```
u_*^2 = ( stress_mag + stress_shelf ) / rho_0
```

and NOTHING else supplies their friction velocity.

| | `stress_mag` | `stress_shelf` |
|---|---|---|
| what | `\|tau\|`, and only `\|tau\|` | `\|tau_top\|` at an ICE-SHELF base |
| derived from | the `tau` pair, always | nothing on this slot |
| refreshed by | EVERY writer of `tau` (`set_wind_stress_*`, the data-forcing seam, the sea-ice blend, `apply_cover`), through `ocean_surfstress_refresh_mag` | the RK2 stage drivers (inline), or `engine_step_finalize` |
| zero where | `cover_frac = 1` (the cover mask zeroes `tau` on every face touching a covered cell) | `cover_frac = 0`, and everywhere without a cavity |
| allocation | always, full size | always, full size |

**The sum is not a double count.** The two supports are disjoint under the
binary v1 cover, and each term already carries its own area weight — `tau` is
masked by `1 - max(cover_L, cover_R)` and `stress_top` is multiplied by
`cover_frac` — so the sum is the area-weighted total upper-boundary momentum
flux and generalises unchanged to a fractional cover.

**Why a SECOND field and not a blend into `stress_mag`.** Three reasons, each
of which a folded-in design gets wrong:

1. `stress_mag` is rebuilt from `tau` **from scratch** by every `tau` writer.
   A folded contribution is silently wiped at the next data-forcing bracket or
   sea-ice blend — the exact failure mode the `stress_mag` refresh contract
   above exists to prevent, inverted.
2. The top-drag stress is recomputed EVERY RK2 stage; the `tau` refresh is per
   OUTER step.  A fold would either accumulate across stages or go stale.
3. Keeping `stress_mag` a pure function of `tau` is what lets the existing
   cover gate (`test_ocean_cavity_flux`: "`stress_mag` is EXACTLY zero under
   cover") and the sea-ice gates (`test_ocean_ice_stress_mag`) still mean what
   they say.

**Who fills `stress_shelf`, and the one honest lag.**

```
&ocean_tdrag_nml enable          ->  run_stage / run_stage_split copy
                                     top_drag%stress_top INLINE, in the SAME
                                     stage that computed it and strictly
                                     before vmix_apply_in_stage reads it.
                                     NO LAG.
melt on, &ocean_tdrag_nml off    ->  engine_step_finalize fills it from
                                     rho_0 * cavity_flux%ustar^2 -- the SAME
                                     C_d (the one-drag-coefficient rule), but
                                     at THERMO cadence at the END of the
                                     step, so it reaches the schemes ONE
                                     OUTER STEP LATE.  Documented, not hidden.
neither                          ->  the zero array; `x + 0.0` is `x`
                                     bit-for-bit, so every pre-cavity
                                     configuration is unchanged.
```

Both fills are written as INLINE `do concurrent` loops, never as a call
handing `ss%stress_shelf` to an external subroutine: a host-gated call with a
state array as an actual makes nvfortran treat the array as escaping and
pessimises every `do concurrent` in the calling routine (CLAUDE.md, measured
at +4.8 % for an inert porous pass).  Both are gated on the source slot's
`enable`, not on `present(...)`: a DISABLED slot carries placeholder-sized
arrays.

**A new upper-boundary stress joins here, not in `stress_mag`.**  The test is
the same one the `p_top` seam uses: if the quantity is derived from `tau`, it
belongs in `stress_mag` and its writer owes the refresh; if it is a separate
momentum flux across the ocean's top, it belongs in `stress_shelf` (or a third
always-allocated companion) and its writer owes an inline per-step fill.
Gates: `test_ocean_bl_under_ice`.

### PGF reference densities (`rho0` vs `rho_ref`) — both follow the one configured ρ₀

`ocean_pressure_force_t` carries TWO reference densities, and they are
different things:

- **`rho0`** — the BOUSSINESQ divisor that turns a pressure gradient into an
  acceleration, `du/dt = −(1/ρ₀)·∂p/∂x`. Read by every variant (`inv_rho0` in
  the face passes, `g_over_rho0` in the Montgomery recursion) and by
  `compute_pbce` in the barotropic coupling.
- **`rho_ref`** — the ANOMALY baseline subtracted from layer densities when
  FV_MOM6 builds its `pa` stack (`pa(top) = ρ_ref·g·η`, layer anomaly
  `(ρ_k − ρ_ref)·g·h`), and the surface `g·ρ_ref/ρ₀` in `compute_pbce`.

They stay SEPARATE MEMBERS — one scales, the other shifts, and interchanging
them is a known MOM6 bug class — but roundabout ships no separate
anomaly-reference knob, so `configure_ocean_pgf` sources BOTH from
`ocean_state%eos%rho0`, i.e. from `&ocean_ic_nml rho_0`: the **single ρ₀ of
record** this state already shares with the EOS, EPBL, kappa-shear, tidal
mixing, wave speed, the `eta_ib` surface-pressure seam, GM / MEKE / Redi / MLE
and the isopycnal slopes. Until that wiring landed the slot kept its 1035
type default while the EOS followed the namelist, so a run with
`rho_0 /= 1035` integrated an equation of state and a pressure gradient on two
different reference densities with no warning. Default `rho_0 = 1035`
⇒ bit-identical; the gate is `test_ocean_pgf_rho_ref`.

Both are plain HOST scalars: every consumer reads them host-side (into a local
`inv_rho0`, or by value into a `*_impl`), never through the device-mapped
`pgf` handle, so the configure-time assignment owes no `!$acc update device`
under `mem:separate`. **A new reference density added to this slot must be
wired the same way** — a defaulted-and-never-assigned one is invisible:
it compiles, runs, conserves, and quietly scales the pressure gradient.

### Forcing + KPP reference densities — `configure_ocean_reference_density`

Four more slots kept their own `rho0` copy and were likewise never assigned,
so the same defect reached the forcing terms. They are wired by ONE routine,
`configure_ocean_reference_density` (`state/rdb_ocean_setup.F90`), called from
`engine_setup` after the per-slot configures and before
`ocean_state_enter_data`:

| Slot | What `rho0` divides | Read from |
|---|---|---|
| `ocean_surface_flux_t` | `dt/(ρ₀·cp)` on every surface heat source and `dt/ρ₀` on every surface salt source — **including whatever the sea-ice coupler writes into `Q_heat`/`Q_salt`** | host (folded into `inv_scale`) |
| `ocean_surface_stress_t` | `τ/(ρ₀·h_top)`, plain and DIRECT_STRESS-distributed | host (by value into `surfstress_*_impl`) |
| `ocean_vmix_t` | KPP `N² = −g/ρ₀·∂ρ/∂z`, `u* = √(|τ|/ρ₀)`, the kinematic `q_T = Q_heat/(ρ₀·cp)` / `q_S = Q_salt/ρ₀` behind `B_0`; also the PP81 and convective-adjustment N² | **device** — `this%rho0` inside `do concurrent` bodies |
| `ocean_geothermal_t` | `dt·Q_geo/(ρ₀·cp)` on the lowest massive layer. Lives on the ENGINE, not on `ocean_state_t`, so it is the routine's optional `geo` argument | host (folded into `src_T`) |

`ocean_vmix_t%rho0` is the one of the four that a kernel reads through the
device-mapped handle. It is correct without an `!$acc update device` **only
because the configure pass runs strictly before `ocean_state_enter_data`'s
`copyin`** — exactly the contract the `pp81_*` / `shear2_floor` scalars beside
it already rely on. Anything that writes one of these after `enter_data` owes
the update.

`vdiff_apply_momentum`'s implicit surface-stress fold (`dt·(τ/ρ₀)/h_nz`) takes
ρ₀ as an ARGUMENT — the split driver hands it `surface_stress%rho0`. There is
deliberately no default: omitting it while the fold is active **fails loud**.
It used to fall back to a silent literal `1035`, which is how a run configured
at a different `rho_0` could fold its wind stress in on the wrong density.

Out of scope on purpose, and NOT copies of ρ₀: `RHO_WATER` (console
mass/heat/salt stats, OBC `mass_out`, ice), `&ocean_ice_nml rho_ocean`, the
`ICE_RHO_*` constants, and the dead `&nonhydrostatic_nml rho_0`.

### The `eta_forcing` seam contract (tides + SAL + surface pressure; consumed by future ice loading)

The barotropic substep drives `−G·∇(η − eta_forcing)` where `eta_forcing` is a
cell-centred `(nx, ny)` field in **metres**, valid **including ghosts**,
device-resident, refreshed **once per outer step** and held **static across the
inner substep loop**. It is the sum of every equilibrium-elevation forcing:
`eta_forcing = eta_eq (body tide) + eta_sal (scalar SAL) + eta_ib (surface
pressure)`, with

```
eta_ib = − p_surf / (eos%rho0 · dyn%bt_work%g_bt)
```

**Note the minus sign** (a high depresses SSH), and **note `g_bt`** — the
barotropic-substep gravity the seam actually feeds, NOT `rdb_constants:GRAVITY`
(building the seam with `GRAVITY` mis-scales the response by ~3.4e-4 and would
silently break isostatic cancellation for a future ice load). At most one
producer's combined array is passed per step: `ocean_p_surf_t%eta_seam` folds in
the tide when both are on (`rdb_ocean_dyn.F90`, the 3-way `run_stage_split`
dispatch), else `ocean_tides_t%eta_forcing`. A new surface load adds its own
input component to `ocean_surface_flux_t` and extends the single overwrite that
builds `sf%p_surf` (full overwrite from the pristine `p_surf_atm` base, never
`+=`); `eta_ib` is always built from the assembled total `sf%p_surf`. The seam is
split-solver only and mutually exclusive with `&ocean_bt_nml bt_halo > 0`: an
EXPLICIT width fails loud in `validate_config`, and the `bt_halo` AUTO default
resolves to 0 (psurf is in `bt_halo_auto_exclusion`, alongside the tide that
shares the seam).

**The partition rule (what belongs on this seam and what does not).** The seam
owns the **barotropic** response to a surface load, and owns it *alone*. The
reason is structural, not conventional: `run_stage_split` subtracts the depth
mean of the layer PGF from the barotropic forcing
(`F_bt_u_fast = F_bt_u − ⟨PFu⟩_h`, `rdb_ocean_dyn.F90`) and then folds the
barotropic solution back over the layers, so the column-integrated FV pressure
force is *discarded* and the substep's `−G·∇(η − eta_forcing)` is the only
barotropic term there is. Two consequences you must carry into any load work:

- A **depth-uniform** contribution to the layer PGF is annihilated — it does not
  reach the baroclinic modes (its deviation from the depth mean is zero) and it
  does not reach the barotropic mode (the depth mean is replaced). So putting a
  depth-uniform top load into the PGF boundary condition
  (`&ocean_pgf_nml p_top_in_bc`, see the `p_top` contract below) is **not** a
  double count of this seam: the two are orthogonal by construction, and each
  covers what the other cannot. (Exactly true for the default uniform-Δu BT
  corrector; under `&ocean_bt_nml correction_h_weighted` the redistribution
  leaves a residual `−dt·δPFu·(w_k − 1)` whose depth mean is still zero by
  `⟨w⟩_h = 1`, i.e. a shear-only redistribution of an already-cancelled force.)
- A load a future **datum** absorbs (an ice draft moved into `bt_H_ref`, so that
  `η ≈ 0` at rest) must **NOT** also be added to `sf%p_surf`, because `eta_ib` is
  built from the assembled total and would re-inject it. Only the load
  *anomaly* — what the datum does not already carry — belongs here. The
  configure-time invariant is `ρ₀·g·z_draft + (bt_H_ref − b)·ρ₀·g ≡ 0`.

`compute_pbce` is deliberately load-blind: `pbce` is `∂p_k/∂η`, the response to a
CHANGE in `η`, and a static top load has none, so `p_top` must not enter it (nor
`gtot_*`, nor `e_anom`) — the bc-PGF retro-correction stays exactly what it was.

### The `p_top` seam contract (surface load in the IN-SITU EOS pressure, E3)

**There are three distinct pressures on this path. Conflating any two of them is
a physics bug with no symptom.**

| | field | units | horizontally varying? | what it is |
|---|---|---|---|---|
| gradient of the load | `ocean_p_surf_t%eta_seam` → barotropic `eta_forcing` | m | yes — only `∇` is physical | inverse barometer; a uniform load is provably inert (gauge invariance) |
| magnitude of the load | `multilayer_state_t%p_top` | Pa | **yes** | the top of the column an IN-SITU hydrostatic pressure is measured down from |
| potential-density reference | `eos%p_ref` (`&ocean_eos_nml p_ref`) | Pa | **NO — scalar by design** | the pressure `ms%rho_layer` is referenced to |

The third row is the trap. `ms%rho_layer` is a **potential** density and its
consumers difference it **along a layer** (the Montgomery PGF's
`rho_layer(i) − rho_layer(i−1)`, the FV-lite / FV-MOM6-PCM integrands) and
**vertically** (the vmix N² builders). Give it a reference pressure that varies
with `(i,j)` — a sloping ice draft, say — and two columns holding *identical*
water at the *same* geopotential depth come out with densities differing by
`∂ρ/∂p · Δp_ref` ≈ `4.5e-7 × 5e6` ≈ **2 kg/m³** across a calving front: a large,
entirely fabricated along-layer density gradient, and therefore a fabricated
pressure gradient force — the exact cavity pathology the seam exists to avoid.
So the load goes into **in-situ** pressures only, never into `p_ref`.
`test_ocean_eos_p_top`'s `rho_layer_independent_of_p_top` is the standing guard
(it ramps `p_top` across the domain over uniform water and demands the density
come out exactly uniform); `p_top_gate_off_leaves_p_top_zero` is its end-to-end
twin. Because N² is built by differencing `rho_layer`, N² is *unaffected* by
`in_eos` and stays self-consistent — that is correct, not a gap.

`multilayer_state_t%p_top` is a cell-centred `(nx, ny)` field in **Pa**, `>= 0`,
valid **including ghosts**, device-resident, refreshed **once per outer step**
and held **static across the stages**. It is **always allocated and zero-filled**
(`ms%init`), mapped in `ms%enter_data`, counted in `ms%bytes()` — never
conditionally — so every kernel has ONE code path: no optional dummy, no
assumed-shape dummy, and no host-gated call handing a state array to an external
subroutine. With the knob off `p_top` is the zero array and `p_top(i,j) + p` is
`p` bit-for-bit under IEEE-754.

Rules for a builder that joins this seam:

- **Only an IN-SITU builder may join.** The test is whether the quantity is a
  true per-layer hydrostatic pressure for *that* column. `p_centre_seed` in
  `eos_wright_pgf_column_sweep_impl` is (it joined); `eos%p_ref` is not (it must
  not).
- **The load is the top of the column.** Replace the `p = 0` seed with
  `p = p_top(i,j)` and accumulate downward. The pressure that reaches the EOS
  must be `>= p_top > 0` and increase toward the bed (`k = 1`) —
  `p_top_sign_and_monotone` asserts exactly that by *inverting* the EOS (no
  duplicated coefficient), because MOM6 shipped a NEGATIVE EOS pressure in one
  path for years.
- **The EOS argument and the PGF boundary condition are two separate consumers.**
  `&ocean_psurf_nml in_eos` gives `p_top` to the EOS's IN-SITU pressure
  *arguments*; `&ocean_pgf_nml p_top_in_bc` (P5.0, FV_MOM6 only, default off ⇒
  bit-identical) adds it to the FV_MOM6 pressure-stack *boundary condition*,
  `pa(nz+1) = rho_ref*g*eta_geo + p_top`. Either, both or neither; they are
  independent knobs on one field, and the per-step inline refresh fires on the
  **disjunction** so neither can read a p_top the configure seed left stale.
  Adding the load at `pa(nz+1)` does **not** double-count the `eta_forcing` seam:
  a depth-uniform `p_top` perturbs every layer's `PFu` by the same
  `−(1/ρ₀)∇p_top` (theorem in `compute_fv_mom6_impl`'s docstring) and the split's
  depth-mean replacement annihilates exactly that — see the partition rule in the
  `eta_forcing` contract above. What it buys under a large load is that `pa`, an
  anomaly stack about `rho_ref*g*z`, stays `O(1e4 Pa)` instead of `O(5e6 Pa)`,
  shrinking the `h_neglect` face-divisor leak by the same factor. FV_WRIGHT's
  `p_edge(nz+1) = 0` and FV_LITE's are still untouched — a separate follow-up.
  `mont` hard-zeroes `M(nz)` and has no stack to inject at all, which is why
  `p_top_in_bc` is refused for every form but `fv_mom6`.
- **A COEFFICIENT evaluated at a pressure is not a third consumer of this
  seam — it is a fourth kind of thing, and the rule differs.** `α` and `β`
  (`eos_buoyancy_coeffs`, E4) are read at a pressure but never *become* one,
  and they are never differenced along a layer, so a horizontally-varying
  argument is safe here in a way it is not for `rho_layer`. Two joiners, with
  two different answers, both deliberate:
  - **`&ocean_vmix_nml buoyancy_coeffs="eos"` in the double-diffusion split**
    is an ordinary in-situ joiner: seed from `p_top(i,j)`, accumulate
    `g·ρ₀·h` downward, exactly as above. `p_ref` must NOT be added — the
    interface has a real depth of its own and would be counted twice.
  - **The KPP `B_0`** is the exception the rule needs stating for. It is a
    SURFACE flux with no depth, so there is nothing to accumulate: it takes
    `p_top(i,j)` under `in_eos` and **`eos%p_ref` otherwise**, which is the
    one place `p_ref` legitimately reaches an α — it keeps the coefficient
    referenced to the same pressure as the `ms%rho_layer` the rest of the
    closure differences. Gate on `in_eos`, not on `p_top /= 0`: a cavity
    fills `p_top` with the ice load either way.
- **Fill it from `sf%p_surf`, whole-array.** The refresh is an inline
  `do concurrent` over the WHOLE array, ghosts included, so `p_top` inherits
  exactly the halo validity `p_surf` has and owes no exchange of its own.
- **An unported in-situ builder is refused, never silently mixed.**
  `validate_config` fails loud on `in_eos = .true.` together with any closure
  that still builds a surface-relative in-situ pressure from 0 Pa: today that is
  kappa-shear, tidal mixing, Redi, isopycnal slopes, the PGF in-layer
  reconstruction, and sea ice. Porting one means deleting its line there in the
  same PR. Two pressure conventions inside one time step have no symptom — that
  refusal is the only thing standing between a half-ported seam and a plausible
  wrong answer. A configuration with *no* in-situ consumer at all is **accepted
  with a warning**, not refused: there the knob is honestly inert.
- **Ported so far: FV-Wright, and EPBL (Phase 4b).** EPBL's column stack
  (`epbl_column_kernel`'s `pres`/`p_mid`) now seeds at `ms%p_top(i,j)` under
  `in_eos`, gated by the host scalar `ocean_epbl_t%in_eos` that
  `configure_ocean_epbl` latches before `enter_data`. **The gate is not
  redundant with `p_top` being zero:** a cavity fills `p_top` with the ice load
  whether or not `in_eos` is set, so only the gate keeps an existing
  cavity + EPBL run bit-identical (`test_ocean_bl_under_ice`'s
  `epbl_p_top_off_is_bit_identical` compares a 3 MPa load against a zero one and
  demands byte equality).
- **A stack with two consumers moves BOTH, or neither.** EPBL's `p_mid` is read
  twice — as the in-situ argument of `eos_specvol_derivs`, and as the PE weight
  `dpe_* = dmass·p_mid·dsv_*`. They are the same pressure (the hydrostatic load
  a layer's centre of mass has to lift, ice included, under a shelf that floats
  with the water it sits on), so the seed moves both. Splitting them would put
  two conventions inside one column, which is what this contract exists to
  prevent. The risk that carries — that a UNIFORM load, pure gauge with no
  gradient, would change the mixing energetics — is discharged by test, not by
  assertion: under a LINEAR (pressure-independent) EOS the offset enters only as
  `pec_core → pec_core + P·colht_core`, and `colht_core` (the column-height
  change of a mixing event) is identically zero because mixing at fixed mass
  conserves `Σ mass·T` and `Σ mass·S` and the height is a fixed linear
  functional of them. Measured difference on gfortran 15.1: **exactly zero**
  (`epbl_uniform_p_top_linear_eos_is_gauge_neutral`). Under the nonlinear Wright
  EOS the same load DOES move the answer, through `α(p)`/`β(p)` alone
  (`epbl_p_top_moves_the_nonlinear_eos_only`) — which is both the point of the
  port and the proof the seed reaches the EOS at all.

### The cavity datum contract (`&ocean_cavity_dyn_nml`, P5.1 + P5.2)

**There are two vertical datums on this path and they must not be
conflated.**

| | quantity | zero point | who owns it |
|---|---|---|---|
| free-surface anomaly | `dyn%bt_work%bt_eta = Σ h_layer − bt_H_ref` | the **loaded equilibrium** (the ice base, under a shelf) | the barotropic solve |
| geopotential height | `pgf%e_face`, `eta_geo = −b + Σ h_layer` | `z = 0` (absolute) | the FV pressure force |

Without a cavity the two coincide, which is why the distinction never had
to be drawn before. With one, `bt_eta ≈ 0` while `eta_geo ≈ −z_draft`, and
both are correct.

The static draft is absorbed into the **datum**:

```
(D)  bt_H_ref  = b − z_draft   afloat                (was: bt_H_ref = b)
              = 0              where GROUNDED  (no water column at all)
(P)  p_ice_ref = (rho_ref*GRAVITY) * z_draft
(I)  rho_ref*g*z_draft + (bt_H_ref − b)*rho_ref*g  ==  0     <=>   (D)
     on every WET column
```

(I) is the **counted-once invariant**, asserted at configure
(`configure_ocean_cavity`): whatever the datum absorbs must NOT also be
handed to the `eta_forcing` seam, because `eta_ib` is built from the
assembled total `sf%p_surf` and would re-inject it. Only the load
*anomaly* — what the datum does not already carry — belongs on that seam.
This is the concrete case the partition rule in the `eta_forcing` contract
above was written for.

**The load partition, end to end (P5.2).** `p_ice_ref` has exactly two
destinations, and the split between them is the whole design:

```
BAROTROPIC  ->  the datum (D).  And nothing else.
PRESSURE    ->  ms%p_top = metrics%p_ice_ref + sf%p_surf,
                read by the FV_MOM6 pa(nz+1) top BC
                (&ocean_pgf_nml p_top_in_bc) and by the in-situ EOS
                (&ocean_psurf_nml in_eos).
SEAM        ->  nothing.  eta_ib = -sf%p_surf/(rho0*g_bt) is built from
                sf%p_surf, which the cavity never writes, so the seam
                carries the load ANOMALY only — which, under the
                Boussinesq-isostatic convention, IS sf%p_surf.
```

Consequences worth stating because a test pins each one: an
inverse-barometer run with no cavity is bit-identical (the cavity adds
nothing to `sf%p_surf`); a cavity with no `p_surf` sends the seam
*nothing*, so switching the seam on at `p_surf = 0` under a cavity is
bit-identical to not having it at all
(`test_ocean_cavity_load`, `test_ocean_cavity_equivalence`).

`ms%p_top` is assembled in two places and they must agree: the configure
seed in `configure_ocean_cavity` (which is the FINAL value for a cavity
without the psurf seam — the draft is static, so there is nothing to
refresh) and the once-per-outer-step inline `do concurrent` in
`ocean_dyn_step_split`, which fires only when the psurf seam makes
`sf%p_surf` live. Both build `p_ice_ref + p_surf`, over the WHOLE array
including ghosts.

**`&ocean_pgf_nml p_top_in_bc` is REQUIRED for a cavity whose draft
varies** — refused fail-loud at configure (and mirrored in
`validate_config` on the namelist shape), not auto-enabled: an
answer-changing knob that a second namelist group switches on behind the
user's back is the silent coupling this codebase fails loud on. A draft
that is UNIFORM over the whole array is exempt, and that is a theorem
rather than a courtesy — a load with no gradient is bit-identically inert
in the top BC — which is what keeps the flat-lid datum-equivalence gate
expressible with and without the load.

Rules for anything that joins this seam:

- **The datum is the whole dynamical effect of a STATIC load.** The split
  solver replaces the depth mean of the layer PGF with the barotropic
  solution, so the column-integrated pressure force is discarded and
  `−G·∇(η − η_forcing)` is the only barotropic term there is. A static
  surface load therefore reaches the barotropic mode through the datum (or
  the seam) and through nothing else — putting it only in the PGF would be
  silently inert. The converse is the reason `p_top_in_bc` is a refusal
  and not a blow-up guard: a MISSING load is also annihilated there, so
  the split path loses conditioning (5.4 decades, measured) rather than
  stability. What it protects is the raw `pa` stack, the unsplit driver,
  and `&ocean_bt_nml correction_h_weighted`, where the uniform piece is
  redistributed as a real per-layer shear.
- **The isostatic load is `ρ₀·g·z_draft`, and it is not the true weight.**
  A stratified column's real overburden is `g∫ρ̂`, which differs by
  `−g∫(ρ̂ − ρ₀)`; the gradient of that difference is a residual
  `N²·z_draft·∇z_draft` in the raw PGF (`5.8e-6 m/s²` at `N² = 1e-5`,
  `z_draft = 280 m`, slope `2e-3` — measured, `test_ocean_cavity_load`).
  It is chosen anyway, because `ρ₀·g·z_draft` is the load that makes the
  DISCRETE BAROTROPIC state exactly at rest, and it is what ISOMIP+
  prescribes. The residual is depth-uniform, so the split annihilates it
  too; `draft_source="in_situ"` (the true isostatic solve) is deferred
  and fails loud.
- **`bt_H_ref` is the reference WATER-COLUMN thickness, not the bed.**
  Anything that re-derives it (the API's bathymetry re-injection, a future
  wide-halo BT clone) must re-derive it as `b − z_draft`, or the ice load
  is counted zero times at that site. The wide clone is refused today for
  exactly this reason.
- **The geopotential stack stays absolute.** `configure_ocean_pgf` gets
  the TRUE bed; `e_face(1) = −b` and the column top lands at
  `−z_draft + η` by itself. Do not "helpfully" shift it — a rigid vertical
  shift of the column is inert in exact arithmetic but re-rounds the whole
  `pa` anomaly stack, which is the one thing in the dyn-core that
  legitimately depends on absolute z (see
  `test_ocean_cavity_equivalence`'s module docstring for the measured
  size).
- **Grounding goes through the wet mask, never through a thin film.**
  `b − z_draft < h_min_cavity` ⇒ the column is LAND via
  `seed_wet_mask_impl`, so the existing metric-zeroing land mask and the
  land-state contract below do the rest.
- **A grounded column's datum is `0`, not a negative water column, and
  the counted-once invariant is asserted over the WET columns only.**
  `bt_H_ref` is the reference WATER-COLUMN thickness; a grounded column
  has none, and `b − z_draft` there is negative by hundreds of metres.
  `cavity_datum_impl` writes `0`; `cavity_datum_residual` skips those
  columns. That is not a weakening of (I): (I) is a statement about the
  BAROTROPIC MOMENTUM EQUATION, and a land column has none — every face
  metric on it is zero, `−G·∇(η − η_forcing)` is multiplied by nothing,
  and there is no load on it to count once or twice. Carrying the
  negative value instead bought three things, all unwanted: a phantom
  few-hundred-metre `bt_eta` on every grounded column (masked out of the
  dynamics, visible in `eta` min/max and the `ssh` diagnostic); an ALE
  target built by CANCELLING two draft-sized numbers to recover the land
  column's `nz·H_VANISHED`, so its thickness jittered at `eps·z_draft`
  instead of sitting bit-stably at `H_VANISHED`; and a grounded column
  arithmetically distinguishable from the ordinary land it IS. Gate:
  `test_ocean_cavity_grounded_budget`.

### The land-state contract (`h = H_VANISHED`, `hTr = 0`)

**A LAND T-cell holds `h_layer = H_VANISHED` exactly and `hTr = 0`
exactly — at `t = 0` and at every step after, however it became land.**
Filled by `ocean_state_seed_land_cells` (`state/rdb_ocean_state.F90`).

`h` is pinned AT the D4 vanish marker, so a land layer is on the
VANISHED side of every `h > H_VANISHED` gate in the tree. `hTr = 0` is
the content those gates already give it — in particular the ALE remap's
concentration step (`rdb_ocean_remap::ocean_remap_tracer_field`:
`c = hTr/h` if `h > H_FLOOR` else `c = 0`), which writes `hTr = 0` on
every land column at the first regrid.

**The seed's job is to hand the budget latch the state every later step
will have.** The console `Error` is `[(total(t) − total(0)) + out −
src]/total(0)`, and `total(0)` integrates over land like everything
else — `compute_total_h` / `compute_total_tracer` weight by `areaT`,
which the land mask does NOT zero. So any land content the seed writes
that the first regrid then discards is a step change in the residual
between step 0 and step 1, un-budgeted by construction.

This replaced a "recover `val = hTr/max(h_old, H_VANISHED)` and re-scale
onto the floor" hold that was wrong twice. The `max(...)` divisor makes
the re-scale an exact IDENTITY whenever `h_old ≤ H_VANISHED` —
including `h_old < 0`, which is what an ice-shelf column GROUNDED by
`h_min_cavity` has — so those columns kept a full-column, negative
`hTr` beside a floored `h`, an implied concentration of order `−1e7`
PSU; on `validation_examples/ocean/isomip_plus/ocean0_idealised_draft.nml`
that was **61 % of the initial salt content** and −8 % of the heat, while
mass closed at `−5e-14`. And even where the re-scale DID work (ordinary
land, `h_old > H_VANISHED`), the `val·H_VANISHED` it left was discarded
by the first regrid — a ~1e-7 relative step change in every land-bearing
case, which is the same defect at a magnitude nobody had looked at.

The `0*NaN` hazard the old hold existed to avoid is avoided the same
way: `T = S = hTr/h = 0` is finite.

**Land is provably isolated from wet cells**, and the gate says so
bitwise rather than by argument: `test_ocean_cavity_grounded_budget`'s
`grounded_matches_plain_land_in_every_wet_column` makes the same 64
cells land two ways — grounded under a 600 m draft over a 400 m bed, and
ordinary `b = 0` island bathymetry — and finds every WET column's `h`,
`hS` and `hT` bit-for-bit equal after 20 steps, across land states that
differ by `bt_H_ref` (`−200 m` vs `0`) and `p_top`
(`ρ_ref·g·600` vs `0`).
- **What is left at rest is the sigma-coordinate PGF truncation, and it
  has a formula.** Flat bed + `VCOORD_SIGMA` + a draft slope `s` gives
  `PFu(k) − ⟨PFu⟩_h = N²·D³·(3σ_k²−1)/(12·dx·H̄)` with `D = s·dx` —
  second order in `dx`, CUBIC in the slope, independent of `nz`. Measured
  end to end: `7.08e-8 m/s` after 12000 s at `s = 1e-3`, `N² = 1e-5`
  (the derived `a_peak·t` is `1.04e-7`). With `N² = 0` it vanishes and
  the run sits at `2.6e-13 m/s`. That number is the baseline for the
  sloping-coordinate PGF corrections, which are a later phase.
- **Ordering is load-bearing.** `z_draft` is filled inside
  `ocean_state_seed_from_cfg`, immediately after the bathymetry and BEFORE
  the layer split and the wet-mask seed read `b − z_draft`; it is then
  re-wrapped (periodic/fold) and halo-exchanged alongside `b` and
  `bt_H_ref`, which were latched from the UNWRAPPED arrays. A draft whose
  seam disagreed with the datum's would count the load twice at that face.

### Surface-flux COMPONENT OWNERSHIP (one writer per component)

`ocean_surface_flux_t`'s component set has exactly one legal writer per
field. `Q_heat`/`Q_salt` are DERIVED VIEWS and have exactly one writer
too — `ocean_surface_flux_assemble`. A second writer on any row is a
silent clobber, not a merge.

| Component | Owner | Sign / units | Notes |
|---|---|---|---|
| `q_sw`, `q_lw`, `q_lat`, `q_sens` | a forcing reader (`rdb_ocean_data_forcing`) | W/m², positive DOWN | `q_sw` also gates `has_q_sw` |
| `heat_added` | the sea-ice coupler (`ice_ocean_heat_flux`) | W/m², either sign, positive down | **FULL OVERWRITE.** Do not add a second contributor here |
| `salt_flux` | the sea-ice coupler (`ice_ocean_brine_flux`) | positive SALINIFIES | **FULL OVERWRITE.** Same rule |
| `heat_cavity` | the ice-shelf cavity (`ocean_cavity_flux_step`) | W/m², positive down; `= −q_ocean` | Separate from `heat_added` *because* the ice coupler overwrites that one. Melting cools |
| `salt_cavity` | the ice-shelf cavity (`ocean_cavity_flux_step`) | same as `salt_flux`; `= −m_mass·(S_far − s_ice)` | The fixed-mass dilution equivalent. Under `&ocean_cavity_melt_nml freshwater="virtual"` (default) it IS the meltwater's effect on salinity; under `"mass"` it stays assembled UNCHANGED as the `B_0` buoyancy forcing KPP/EPBL read, and `ocean_cavity_mass_step` removes it again from the tracer. Melting freshens either way |
| mass fluxes + their `heat_content_*` twins | a forcing reader | see the type's docstrings | a mass flux owes its enthalpy companion |
| `heat_content_massin` / `massout` | **the assembler** | W/m² | never a filler |
| `p_surf_atm` | a reader / configure seed | Pa | `p_surf` is the assembled total |
| `Q_heat` / `Q_salt` | **the assembler** | derived | a filler writing these is a bug |

#### Real mass: what IS and what is NOT a column-mass change

The mass-flux components `evap`/`lprec`/`fprec`/`vprec`/`lrunoff`/
`frunoff`/`seaice_melt` and their `heat_content_*` companions remain
**enthalpy + salt bookkeeping only — they do NOT change column mass.**
That has not moved.

What HAS moved is the ice-shelf cavity. With
`&ocean_cavity_melt_nml freshwater="mass"` the basal meltwater is a REAL
Boussinesq VOLUME source on the top layer — `dh = m·dt/ρ₀` added to
`h_layer(:,:,nz)` by `ocean_cavity_mass_step`, in-stage, at thermo
cadence, immediately after the surface-flux apply — with `d(h·S) =
dh·s_ice` and `d(h·T) = dh·T_b`, and the salinity falling by dilution on
its own.  So, precisely:

| path | real column mass? |
|---|---|
| ice-shelf basal melt, `freshwater="mass"` | **YES** — the only one |
| ice-shelf basal melt, `freshwater="virtual"` (default) | no |
| sea ice (`rdb_ice_ocean_coupler`: `salt_flux`, `heat_added`) | no — virtual, and it stays virtual |
| the seven atmospheric/river mass-flux components | no — named follow-up, unchanged |

Consequences worth stating in the contract, because they are seams:

* **The barotropic mode.**  Nothing is added to the barotropic forcing.
  `derive_bt_from_layers` rebuilds `bt_eta = Σ_k h − bt_H_ref` at the top
  of every stage, so the volume is in `bt_eta` one stage later and the
  free surface under the cavity datum simply rises.  `bt_H_ref` stays
  static (it is the ice-base-to-bed datum, not a sea level), and every
  consumer already reads the water column as `bt_H_ref + bt_eta`.
* **`ms%p_top` is prescribed and does not respond.**  The ice draft is
  static, so a rising `η` under a fixed draft does not push the ice up;
  the added volume must leave through an open boundary or raise the OPEN
  ocean's surface.  For a CLOSED domain that is a real sea-level rise —
  ISOMIP+ Ocean0 is ~30 m/yr of melt over ~1e10 m² of shelf into ~4e10 m²
  of open surface, i.e. metres per year — which is why
  `volume_compensation="uniform_open_ocean"` exists.
* **Budgets.**  `multilayer_state_t%mass_src` is the tracked mass SOURCE,
  accumulated with the same per-stage weight as `mass_out`, so the
  console residual `(M − M₀) + mass_out − mass_src` stays at round-off
  while the total legitimately grows.  Salt and heat need no new
  accumulator: every increment rides the existing
  `salt_budget_surface`/`heat_budget_surface` contributors.

Obligations on any filler: write your own component (full overwrite,
never `+=`), push it with `!$acc update device` if you wrote it on the
host, latch `has_heat`/`has_salt` (and `has_mass_flux`/`has_q_sw`)
**host-side** and never from a device reduction, and register the
component in the restart registry if it is time-varying —
`heat_cavity`/`salt_cavity` are (the melt rate is a function of the live
state, and it is written at the end of step N and integrated on step
N+1), so they are registered `optional=.true.`

### `p_top`: the melt liquidus is the THIRD consumer

The `p_top` seam contract above names two producers — the cavity's static
ice load and the `&ocean_psurf_nml` seam, summed as
`ms%p_top = metrics%p_ice_ref + sf%p_surf` — and two consumers, the
FV_MOM6 `pa(nz+1)` top boundary condition and the in-situ EOS pressure.
**The ice-shelf basal-melt liquidus is the third consumer**, and it is
the one that makes the field's name literal: `p_top` IS the pressure at
the ice-ocean interface, which is what `eos_freezing_point(eos, S, p)`
must be evaluated at for the ice pump to exist at all.

`rdb_ocean_cavity_flux` reads it and writes NOTHING to it.  There is one
producer, `configure_ocean_cavity` (re-assembled per outer step in
`ocean_dyn_step_split` when the psurf seam makes `sf%p_surf` live), and a
second writer would be a silent clobber — `configure_ocean_cavity_melt`
runs *after* it, so a "helpful" seed there would drop `sf%p_surf` on
every cavity that had one.  What the melt configure does instead is
ASSERT the producer ran: `cavity_count_unloaded_p_top` counts cells with
`p_top < p_ice_ref` (exact, not a tolerance — `p_surf >= 0` by contract)
and fails loud on any.  That catches a skipped or reordered producer,
which is the failure this seam actually has, rather than trusting the
call order.

Note the consequence for the melt physics: because the liquidus reads the
assembled TOTAL, an atmospheric load under the shelf depresses the
freezing point exactly as the ice load does, with no extra wiring.

## How to pick up a slot

1. **Read the slot's module doc-comment** — it states the algorithm
   variant we plan to use, the hand-off boundaries, and any literature
   citations.
2. **Read this README's row** — what your slot reads from, what it
   writes to.
3. **Read `docs/ROADMAP_OCEAN.md`** for the phase number and any
   constraints from the broader plan.
4. **Stay inside the slot.** If the implementation needs a new
   prognostic field, add it to the matching prognostic state
   (`barotropic_state_t` or `multilayer_state_t`), not to
   your slot type — your slot is for *control* and *workspace*.
5. **Never `use mpi_f08` (or `mpi`) directly — anywhere.** If the
   implementation needs a halo exchange or a reduction, call the
   `src/comm/` layer (`solver_halo_exchange_*`, `halo_allreduce_*`).
   That layer talks to MPI via `pic_mpi_lib` from the pic
   dependency, so even the comm backend stays portable.  The
   `no-mpi-in-rdb` pre-commit hook enforces this repo-wide.
6. **Add a test in `tests/`** that exercises only your slot, using
   `ocean_state%<your_slot>` directly. Don't reach across slots.
7. **If your slot is default-off, gate it for memory** (see Rule 2):
   latch its `enable` in `init_from_config` before `call this%init`,
   gate its `init` + parent `enter_data`/`exit_data` on `enable`, and
   add its arrays to `ocean_state_bytes` so the upfront counted report
   accounts for them (and the reconciliation stays quiet).

## Cross-cutting conventions

- **Vertical indexing**: `k=1` is the bed, `k=nz_ml` is the surface
  layer (ROMS-style bottom-up; see top-level `CLAUDE.md`).
- **Face indexing**: `u_face_x(i, j, k)` is at the interface between
  cell `(i, j, k)` and `(i+1, j, k)`. Allocation runs from `i=1` (west
  outer wall) to `i=nx+1` (east outer wall).
- **GPU**: Phase 1+ will add `!$acc enter data` to each slot's
  `init` after the allocations land. Don't add `!$acc` directives in
  Phase 0e — they go in alongside the real allocations to avoid
  partial-presence errors.
- **Initialisation flag**: every slot type carries `logical ::
  is_init = .false.`.  `init` sets it to `.true.` only after all
  allocations + GPU attachments succeed; `destroy` clears it before
  releasing anything.  **Always check `slot%is_init`, never
  `allocated(slot%some_array)`** — `allocated()` only sees the host
  pointer, not the GPU device-side mapping, so it returns true for
  arrays that are CPU-allocated but not yet on the device.  This is
  also the right guard against double-init and use-after-destroy.
- **Loop pattern (strided hybrid)**: column-local kernels — anything
  with `nz`-sized scratch per `(i, j)` (vmix/KPP, continuity-PPM,
  coriolis_adv with q at corners, lateral-mix vorticity-gradient,
  FV pressure-force integration) — **must be written in the
  hybrid strided form** so the same source delivers cache-friendly
  CPU runs and high-occupancy GPU runs.  The knob is
  `ocean_state%column_stride`; slot workspaces are sized
  `(stride, stride, nz_ml)`.  Skeleton:

  ```fortran
  ! Slot's scratch is a type(scratch_3d_buffer_t), sized at init:
  !   call this%scratch%init(stride, stride, nz_ml, "kpp_ri_scratch")
  ! Phase 1+ scratch_3d_buffer%init also performs !$acc enter data.
  associate (stride => ocean_state%column_stride)
     do i_out = 1, nx, stride
        do j_out = 1, ny, stride
           do concurrent ( &
                  k=1:nz_ml, &
                  j_in=1:min(stride, ny - j_out + 1), &
                  i_in=1:min(stride, nx - i_out + 1))
              i = i_out + i_in - 1
              j = j_out + j_in - 1
              ! kernel body — read state, write this%scratch%data(i_in, j_in, k)
           end do
           ! one or more inner do-concurrent passes reusing this%scratch
        end do
     end do
  end associate
  ```

  At `stride = 1` the inner loop collapses to `(1, 1, nz)` =
  pure-column form, the workspace is effectively 1D, CPU caches are
  happy.  At `stride = nx` the outer loop is one iteration, the
  inner collapses over the whole plane, GPU occupancy is maxed.
  Reference: Kommera & Appelhans (NVIDIA, CESM SEWG 2026), "One
  Codebase with Good Performance on both CPUs and GPUs, for
  Column-based Loop Structures".  Workspaces in the strided
  pattern are declared as `type(scratch_3d_buffer_t)` rather than
  bare allocatables — one allocator, one GPU-mapping site, one
  resize site (`../../framework/rdb_scratch_3d.F90`).
- **VCOORD_LAGRANGIAN grounding stability** (`&ocean_isopycnal_nml`,
  all knobs default OFF ⇒ bit-identical): three opt-in robustness
  controls for the remap-free Lagrangian vertical coordinate:
  `angstrom_h` (Phase 1: min-thickness floor on the h-update,
  MOM6 GV%Angstrom_H analogue), `reset_vanished_u` (Phase 2: zero
  face velocity when both adjacent cells are vanished),
  `cfl_ignore_vanished` (Phase 3: exclude vanished layers from
  MaxCFL / panic / CFL truncation). All gates test
  `coord_type == VCOORD_LAGRANGIAN` specifically — other vcoords
  (sigma, zstar, etc.) remain byte-identical.
- **Initial layer thickness** (`&vcoord_nml thickness_config`, default
  `"sigma"` ⇒ bit-identical): selects what `h_layer` is *seeded* to in
  `ocean_state_seed_from_cfg`, independently of `vcoord_type` (which
  selects the *running* coordinate). `"sigma"` is the historical even
  split of the local depth, `h_layer = b/nz_ml`
  (`seed_h_layer_uniform_impl`). `"uniform_z"` is the MOM6
  `initialize_thickness_uniform` port
  (`seed_h_layer_uniform_z_impl`): uniform z interfaces over the global
  `ocean_max_depth`, clipped bottom-up to the local bathymetry, with
  the remainder collapsed to `max(angstrom_h, 2·H_VANISHED)`. The
  surface layer's target is `z = 0`, so it absorbs the remainder and
  the telescoping column sum is EXACTLY `b(i,j)` — the seed floor
  borrows, it does not inject (unlike the runtime `angstrom_h` clamp).
  A per-k room clamp keeps that exact even when the floor is not small
  compared with `max_depth/nz_ml`; degenerate columns
  (`b <= nz_ml·floor`, including dry beds) fall back to the floored
  even split. **This is the knob that matters for layered/isopycnal
  runs**: with a horizontally-uniform density stack the interfaces ARE
  the isopycnals, so `"sigma"` tilts every isopycnal with the
  bathymetry and releases the resulting APE as a shelf-break gravity
  current. Fail-loud against `&ocean_wetdry_nml enable` (that path pins
  `Σ h_layer = nz·2·H_VANISHED` on emerged columns).
- **D4 thin-layer constant taxonomy (ocean path)**:
  - `H_VANISHED = 1.5e-4 m` — skip/merge marker (layer is so thin
    it should be skipped or merged into a neighbour); do NOT use as
    a positivity floor.
  - `H_DIV_EPS = 1e-20 m` — pure 1/0 armour; only used to avoid
    division by zero in denominators that are never actually tiny.
  - `angstrom_h` (runtime, set from `&ocean_isopycnal_nml`) —
    physical D4 floor for grounding stability; lifts h but NOT hTr
    (non-conservative, bounded injection per floored cell). Distinct
    from `H_VANISHED` (which is a skip/merge marker, not a floor).
    Lives on `continuity_t%angstrom_h`, gated to VCOORD_LAGRANGIAN.
- **Safe math (library, no production consumer)**: `framework/rdb_safe_math.F90`
  carries correct, unit-tested `safe_exp/log/sin/cos/sqrt/pow` wrappers
  + their `*_polynomial` implementations (Cody-Waite reduction + Horner
  Taylor), but the wrappers are today plain elemental inlines around the
  bare intrinsic — no kernel calls them, and no build mode dispatches to
  the polynomial path.  The `RDB_BITWISE_REPRO` CMake option that
  used to select that dispatch was **removed in PR-8**: it validated and
  compiled but changed nothing (zero `safe_*` call sites anywhere in
  `src/`), so keeping it was a documented lie about a capability nothing
  used.  The polynomial implementations remain as tested raw material
  for a future bitwise-reproducibility PR — do not re-add the build
  option without first wiring at least one kernel to call `safe_*`.
- **Units registry — PROPOSED, NOT BUILT.**  The design called for a
  `units_registry_t` on `ocean_state%units` (`framework/rdb_units.F90`)
  carrying units, valid range and CF metadata for every named field,
  with "if you introduce a new field, register it" as a standing
  convention.  None of it exists: there is no `rdb_units.F90`, no
  `units` slot, and nothing to register against — so there is no
  convention to follow here today.  Field metadata currently lives
  where it is used: `tracer_t%units`/`long_name`/`standard_name` for
  the tracer registry, and the diag registry's own per-variable
  attributes for output.  Kept as a design note because the motivation
  still stands (one source of truth for the diag writer, restart
  manager and tracer display, instead of `"m s-1"` literals scattered
  across kernels) — but treat it as unbuilt work, not a rule.
- **`bt_rem_u/v` is a multiplicative accumulator — reset it exactly once
  per stage.** `barotropic_workstate_t%bt_rem_u/v` is allocated
  `source=1.0` once at init; thereafter the ONLY thing that resets it to
  1 is `compute_bt_rem` (`&ocean_bt_nml substep_drag`) or, when that's
  off, `reset_bt_rem` — every other contributor (currently
  `compute_bt_rem_wave_drag`, `&ocean_bt_nml wave_drag`) MULTIPLIES into
  it. `mask_bt_rem` (land faces, idempotent under `*=`) always runs
  last. Skip the reset and a multiplicative contributor compounds
  geometrically across outer steps (`bt_rem = R^n` after `n` stages),
  silently annihilating the barotropic mode — the single most likely way
  to ship a plausible-looking, catastrophically wrong barotropic-drag PR
  (see `test_ocean_wave_drag::wave_drag_no_compounding`). **Any future
  PR adding a third `bt_rem` contributor must extend the reset dispatch**
  in `rdb_ocean_dyn.F90::run_stage_split` (or refactor to an
  unconditional `reset_bt_rem` + N independent contributors — the
  `vmix_assemble` contributor+single-gate idiom is the model to follow).
- **The barotropic substep's live terms need a frozen-reference
  subtraction — never add a fast-loop term without one.** `F_bt_u/v`
  (the depth-mean slow forcing handed to the substep) already contains
  the depth mean of EVERY slow layer tendency; any term the fast loop
  ALSO integrates live is double-counted unless its value at the
  stage-entry BT state is subtracted from the forcing.  Two subtractions
  exist today, both in the `F_slow` assembly in
  `run_stage_split`: the PGF projection (`F_bt_u_fast = F_bt_u −
  depth_mean(PGF)`; the substep's own `−G·∇η` replaces it) and the
  Coriolis/advection reference (`subtract_fast_cor_ref` — MOM6
  `Cor_ref_u/v`; without it the barotropic Coriolis is integrated twice
  and the Δu corrector hands every layer an extra `dt·f·v̄` rotation per
  stage).  Related trap: dissipative slow terms in `F_bt` are FROZEN
  over the outer step, so a grid-scale damping rate `λ·dt ≳ 0.5` is
  applied with reversed phase to fast modes (`ω·dt ≳ 1`) and
  ANTI-damps them — the `&ocean_hvisc_nml bound_kh` clamp exists to
  keep the lateral viscosity out of that regime (600² Lagrangian
  double-gyre blow-up; `test_ocean_bt_cor_ref`,
  `test_ocean_hvisc_kh_bound`).
- **Diagnostics interaction**: kernels are diag-agnostic. Don't add
  `diag%send(...)` calls inside your slot's kernel; the diag manager
  reads from `state` via the registry. If your slot exposes a new
  derived quantity, add it to the diag registry's variable list (in
  Phase 6) — the kernel itself stays untouched.
- **Restart interaction** (shipped): a new prognostic-owning slot joins
  the checkpoint by registering its arrays in
  `ocean_state_build_restart_registry` (in `state/rdb_ocean_state.F90`)
  via `reg%register_2d` / `reg%register_3d` / `reg%register_scalar`
  (or the `register_full_3d` helper). Each entry defaults to **required**
  (a missing field on read is fatal); mark genuinely optional ones with
  `optional=.true.`. Host-only scalars (Chapman corner state) register
  with `register_scalar` and are skipped by the device-pull walk. The
  manager (`ocean_state_restart_write` / `_read`) walks the registry,
  does the `!$acc update self` device pulls on write (device-mapped
  entries only), writes each field as a **full local array** (interior +
  ghosts), and writes a per-rank `restart_rank_NNNNNN.nc` durably (to
  `.tmp` then POSIX-rename) with decomposition + grid/vcoord/tracer +
  schema metadata, all validated on read; you don't write per-slot I/O
  code. Read runs BEFORE `ocean_state_enter_data` so the H→D copy
  carries restored state up.
  v1 ghost-cell policy: structured prognostic + closure arrays are
  saved as full local arrays (physical wall ghosts are genuine owned
  boundary state with no stage-entry rebuild op; periodic/fold/halo
  ghost columns are redundant and refill from the interior on resume).
  Same-decomposition resume only (errors on mismatch); cross-rank
  redistribution is a future offline tool. Diag MEAN/MAX/MIN windows
  reset on restart by design — the bit-exact gate
  (`test_ocean_restart`) covers the prognostic trajectory, not
  mid-window diag aggregates. The ice->ocean EVP momentum mediation is
  restart-carried on the ICE slot (`ice_tau_ocn_x/y` + the
  `ice_tau_ocn_valid` presence scalar, PR 63) and folded into
  `surface_stress%tau_x/y` after the configure-time wind re-seed —
  `ice_ocean_stress_resume_apply` COPIES the exact last-blended value
  instead of reconstructing it, mirroring the `salt_flux_diag -> Q_salt`
  fold; gated by `restart_bit_exact_ice_evp`
  (`tests/test_ocean_ice_restart.F90`).
- **Order-invariant global reductions (PR-32)**: `src/framework/rdb_efp.F90`
  (Extended-Fixed-Point / Hallberg & Adcroft 2014) + the comm-facade
  `halo_allreduce_efp_list` (`src/comm/rdb_halo.F90`)
  are the seam for any future global integral that needs to be
  rank-count- and decomposition-order-reproducible — a plain
  `!$acc parallel loop reduction(+:acc)` + `halo_allreduce_sum` drifts
  at round-off across decompositions.  `rdb_ocean_console_stats.F90`'s
  `&ocean_diag_nml reproducing_sums` knob is the first (and so far
  only) consumer; see `docs/CAPABILITIES_AND_LIMITATIONS.md`'s
  conservation-contract section for the achievable guarantee.  Use
  `efp_real_diff` (not a double subtraction of two EFP-derived reals)
  when the quantity being reported is a DIFFERENCE against a latched
  reference — the exactness of the sum alone does not survive a
  double subtraction. Next contributor needing a global integral:
  reach for this seam before adding a fifteenth naive allreduce.

## What's *not* slotted yet

Open-ended extensions deliberately left without a type until there's a
concrete need:

- **Sea-ice coupling** — the gated `ocean_sea_ice_t` slot now exists
  (`src/core/ice/state/rdb_ice_state.F90`, `&ocean_ice_nml enable`
  default off ⇒ byte-identical) and PR 1 landed the two OCEAN-side
  prerequisites: `eos_freezing_point(S, p)` (linear liquidus
  `T_f = λ1·S + λ2 + λ3·p`, every EOS variant, `rdb_eos`; the triple
  rides the `eos_t` handle and is chosen by `&ocean_eos_nml
  tfreeze_set` — `"seaice"` = SIS2 default ⇒ bit-identical, `"isomip"`
  = the ISOMIP+ cavity set) and the frazil accumulator
  (`rdb_ice_frazil` — once per outer step on the post-RK2-average
  state the driver clamps the surface layer at T_f and banks the
  supercooling deficit `ρ·Cp·h·(T_f−T)⁺` into the restart-carried
  `ice%frazil_heat` (J/m²); the matching hTr increment feeds a
  "frazil_heat" `budgets` contributor + the console heat closure at
  FULL weight, `ocean_frazil_heat_src`).  PR 2 landed the enthalpy
  library (`rdb_ice_enthalpy` — pure elemental closed-form T↔E
  inversions with `Cp_brine = Cp_ice`, exact round-trip; leaf, nothing
  calls it until PR 3).  PR 3a lands the single-column Winton
  thermodynamics: the six category prognostics (`part_size`, `m_ice`,
  `m_snow`, `enth_ice`, `enth_snow`, `sal_ice`, bottom-up k) plus
  `rdb_ice_column` (SEB + conduction solve, ported top-down + flipped
  at the gather/scatter boundary), `rdb_ice_optics` (CSIM4 shortwave),
  and `rdb_ice_mass` (bottom-freeze, melt peel, equal-mass rebalance).
  Exercised only by `tests/test_ocean_ice_column.F90` at PR 3a; PR 3b/
  3c wired the frazil-bank spend + atmospheric-forcing seam + basal
  flux + column driver + ocean-coupling diags into the driver at
  thermo cadence, and PR 4a added the multi-category ITD restore
  (`rdb_ice_itd%ice_adjust_categories`, `ncat>1` mode).  PR 4b lands
  horizontal category ice/snow transport + `compress_ice`
  (`rdb_ice_transport%ice_transport_step`, `&ocean_ice_nml transport`
  default off ⇒ byte-identical) — a category-summed PPM (reusing
  `rdb_continuity`'s five promoted helpers) with a proportionate
  per-category flux split, PCM tracer riding, and the thinnest-first
  `compress_ice` cascade that keeps `Σ_c part_size ≤ 1`.  The SIS2-port
  ladder that fills the slot (enthalpy → Winton thermo → transport →
  EVP) is `PLAN_SEA_ICE.md`; EVP dynamics (PR 5) landed —
  `rdb_ice_evp` writes `ice%u_ice`/`v_ice` from the C-grid EVP momentum
  solve when `&ocean_ice_nml dynamics=.true.`; the ocean-surface-
  velocity sampler remains the `dynamics=.false.` fallback.  PR-58 adds
  `&ocean_ice_nml hlim` — a per-run override of the ITD's category
  lower-thickness edges (`ice_itd_category_bounds`'s optional
  `hlim_vals`, `rdb_ice_state.F90`), for domains (e.g. a thick Antarctic
  pack) where the hardcoded Arctic-tuned SIS2 default table
  (`[1e-10, 0.1, 0.3, 0.7, 1.1, 1.5, 2.0, 2.5]` m) leaves the whole ITD
  sitting in one category.  Unset (`-1.0` sentinel, the default)
  reproduces the SIS2 table bit-identically.  PR 24 lands
  the ANALYTIC initial-condition path (`rdb_ice_init`,
  `&ocean_ice_ic_nml conc_config="uniform"/"latitudes"`, default `"zero"`
  ⇒ byte-identical): a run can now *start* with a live ice pack — mass
  binned host-side against the already-computed `ice%mh_lim` so the
  seeded state is an exact fixed point of `ice_adjust_categories`, and
  `enth_ice`/`enth_snow` set from a namelist temperature through the
  exact `ice_enth_from_ts` inversion — instead of growing one from
  frazil. File-backed ICs (`"file"`) remain deferred to PR-14 (v1.1).
  PR 26 lands the snowfall source term (`&ocean_ice_nml snowfall`,
  default 0 ⇒ byte-identical): `rdb_ice_atm_forcing` spreads a uniform
  frozen-precipitation rate onto the fourth atmospheric-seam field
  `ice%atm_fprec`, and `ice_snow_accumulate` (`rdb_ice_mass`) adds
  `fprec·dt_therm` to `m_snow` inside `ice_column_step`'s resize stage
  (same position as SIS2's `ice_resize_SIS2`, before the melt peels, so
  new snow is meltable the same window).  The share that lands where
  there is no ice — open water at `ncat>1`, ice-free cells at
  `ncat==1`, any category failing the column's own entry gate — is
  NOT orphaned onto a zero-ice category (SIS2's behaviour, which would
  trip `rdb_ice_transport`'s fail-loud `mca_snow>0`-with-`mca_ice<=0`
  reduction here); instead `rdb_ice_snow%ice_snowfall_ocean_share`
  delivers it to the ocean as a latent-heat sink + virtual freshening
  through the existing `heat_flux_diag`/`salt_flux_diag` contributor
  seam — a documented divergence from SIS2, tracked as a PR-16
  (real ice↔ocean mass) cleanup item.
  `&ocean_ice_nml a_face_stress` (default off, requires `dynamics=.true.`)
  weights the EVP momentum balance's wind stress AND ice-ocean drag by
  the face ice concentration `a_u = 0.5*(ci(i-1,j)+ci(i,j))`, closing the
  ice<->ocean momentum budget exactly at every fractional cover (Hibler
  1979 eq. 1); default off leaks `(1-a)*(tau_a-fxoc)` at fractional cover
  (zero at `a in {0,1}`/steady free drift) and drives a spurious Nansen-
  free-drift ghost velocity over ice-free water.
  `&ocean_ice_nml cfl_trunc` (SIS2 `CFL_TRUNCATE`, default 0.5 there, `0`
  here ⇒ byte-identical) clips the FINAL ice velocity to the transport-CFL
  bound `0.95*cfl_trunc*areaT(donor)/(dt_transport*dy_cu)`, demoting
  `ice_transport_step`'s conservation/positivity `error stop` to a
  driver-logged warning + backstop rather than the only defence against a
  runaway EVP drift; `cfl_trunc_dyn_its` (default off, matches SIS2) also
  clips at the bottom of every subcycle. `&ocean_ice_nml project_ci`
  (SIS2 `PROJECT_ICE_CONCENTRATION`, default `.true.` there, `.false.`
  here ⇒ byte-identical) projects `ci` forward along the current
  divergence each subcycle (`ci_proj = ci*exp(-t_cum*sh_Dd)`) and
  recomputes `pres_mice` from it, stiffening the rheology under
  convergence within the call instead of holding the pre-loop strength
  static for all `evp_sub_steps` subcycles — provably inert at a
  saturated (`ci=1`) jam, so it is not a fix for downwind mass pile-up
  (that needs ridging, Tier 3, not yet scheduled).
  PR 27 lands Archimedes snow-ice
  flooding (`&ocean_ice_nml snow_ice`, default off ⇒ byte-identical):
  `ice_snow_ice_flood` (`rdb_ice_mass`) runs in `ice_column_step`'s
  resize stage, between the melt peels and `ice_rebalance_layers` —
  SIS2's exact placement (`ice_resize_SIS2`'s final substantive block,
  `SIS2_ice_thm.F90:1303-1320`). When the snow load submerges the
  snow-ice interface below the waterline, it converts exactly enough
  snow mass into the top ice layer to restore Archimedes equilibrium in
  one non-iterative step, mass-weighting the snow's enthalpy in and
  diluting the layer's bulk salinity — exactly mass/enthalpy/salt-
  conserving within the column.  `snow_to_ice` (SIS2 `SN2IC`) is
  diagnostic-only: SIS2's own formulation exchanges nothing with the
  ocean, so this closure needs no PR-16 seam.  Known limitation
  (inherited from SIS2, not fixed here): the converted ice carries the
  SNOW's enthalpy and ZERO salinity rather than the flooding
  SEAWATER's, so the new ice is too cold and too fresh relative to
  real (seawater-drawing) flooding — that closure is a materially
  larger, separate PR.
- **Ice-shelf cavity basal melt** — would be `ocean_cavity_t`.  The
  PHYSICS has landed as a slot-free leaf module,
  `src/parameterizations/vertical/rdb_ocean_cavity_melt.F90`: scalar,
  `pure`, `!$acc routine seq` three-equation interface thermodynamics
  (Holland & Jenkins 1999; Jenkins et al. 2010; ISOMIP+ /
  Asay-Davis et al. 2016; Yung et al. 2025), with a `Gamma·u*`
  constant-coefficient law, the full H&J99 turbulent+molecular law and
  the Yung et al. (2025) StratFeedback law behind one exchange-law
  enum, insulating / H&J99 advective-diffusive ice conduction, the
  closed-form (cancellation-safe) quadratic in `S_b`, and a
  `CAVITY_MELT_*` status out of every solver instead of an
  `error stop`.  It reads the liquidus off the shared `eos_t` handle
  (`&ocean_eos_nml tfreeze_set`), so it carries no coefficients of its
  own.  **KERNEL ONLY — nothing in the ocean step calls it yet**: no
  namelist group, no state slot, no engine wiring, and therefore no
  bit-identity risk.  Gated by `tests/test_ocean_cavity_melt.F90`
  against a 48-case golden oracle (agreement ~1e-15 relative, asserted
  at 1e-12).  The coupling PR is what needs the slot: per-cell
  `u*` from the surface-layer velocity + `&ocean_ice_nml`-style drag,
  the cavity mask + ice draft (which is also the `p_top` /
  `eta_forcing` seam's customer), the far-field sampling depth in
  METRES rather than layer `nz` (the dominant resolution artefact in
  the literature), and the tracer-budget wiring for the melt heat/salt
  fluxes the kernel already returns.  It should EXTEND the module's own
  `cavity_melt_columns` (a `do concurrent` over columns, living beside
  the scalar solver) rather than write a column loop in the engine: on
  the GPU build nvlink cannot resolve a `!$acc routine seq` device
  symbol out of `librdb_core.so` into a `do concurrent` compiled in a
  different translation unit, so the loop and the routine it calls must
  share an object.
- **Biogeochemistry** — would be `ocean_bgc_t`. The shared
  `tracer_t` registry on `ocean_state%multilayer%tracers(:)` already
  takes BGC tracers as additional entries via
  `multilayer%register_passive_tracer(...)` (see
  `docs/howto/add_passive_tracer.md`); a dedicated slot appears
  only when a real BGC module (NPZD, BLING, COBALT-clone) needs its
  own state beyond the per-tracer fields. Note the ocean diagnostic
  surface is NOT registry-driven (`fill_salinity`/`fill_age`/… each
  name their index) — a namelist-declarable arbitrary-tracer package
  would additionally need a generic registry-driven diag-fill path.
- **Lagrangian particles / floats** — not in the v1 scope.

Adding a new slot follows the same recipe as Phase 0e: declare an
empty `_t`, add it to `ocean_state_t`, wire init/destroy, register in
`src/CMakeLists.txt`, drop a row into the table above.
