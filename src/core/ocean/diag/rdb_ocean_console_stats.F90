!! MOM6-style console-line statistics for the ocean dynamical core.
module rdb_ocean_console_stats
   !! Periodic conservation + stability scalars printed to stdout at the
   !! driver's `status_interval` cadence (separate from the diag manager's
   !! NetCDF). Mirrors MOM6's "MOM Day N:" status line:
   !!
   !!   * `Total mass`    — Σ h·areaT·ρ_0                            (kg)
   !!   * `Total KE`      — Σ 0.5·h·(u_c²+v_c²)·areaT·ρ_0            (J)
   !!   * `Mean salinity` — Σ hS·areaT / Σ h·areaT                  (PSU)
   !!   * `Mean temp`     — Σ hT·areaT / Σ h·areaT                  (°C)
   !!   * `Max CFL`       — max (|u_c|·dt·idxT + |v_c|·dt·idyT)
   !!
   !! Reductions run on device. Tracer sums go through a flat-impl helper
   !! to dodge the registry deep deref (`ms%tracers(idx)%hTr` — NVHPC
   !! can't follow inside an inlined reduction kernel). Initial values
   !! captured on the first call; later calls report relative drift.
   use, intrinsic :: iso_fortran_env, only: int64, real64
   use rdb_constants, only: wp, RHO_WATER
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_console_stats, only: console_stats_t, console_stats_report, &
                                conservation_budget_t
   use rdb_halo, only: halo_allreduce_sum, halo_allreduce_max, halo_allreduce_efp_list
   use rdb_ocean_halo_counters, only: oh_counters_format
   use rdb_ice_column, only: ICE_RHO_ICE
   use rdb_efp, only: efp_t, EFP_DIGITS, EFP_PREC_WIDTH, EFP_MAX_SUMMANDS, &
                      efp_carry, efp_to_real, efp_from_real
   use pic_logger, only: logger => global_logger
   use pic_strings, only: to_string
   implicit none
   private

   public :: console_stats_t
      !! Re-exported from `rdb_console_stats` so the driver keeps a single
      !! import for the ocean path; the reference-snapshot state + the line
      !! format now live in the shared module (coastal uses the same).
   public :: ocean_console_stats_report
   public :: ocean_budget_src, ocean_budget_out, ocean_budget_is_active
   public :: ocean_budget_stage_weight
   public :: ocean_heat_src_sum, ocean_salt_src_sum
   public :: ocean_frazil_heat_src

   real(wp), parameter :: RK2_STAGE_WEIGHT = 0.5_wp
      !! RK2 stage weight: all per-cell budget accumulators (surface,
      !! geothermal, hdiff, horiz_adv) store the raw SUM over both RK2
      !! stages; the outer-step state change from each process is
      !! 0.5·(that sum) because rk2_average halves.  The reporter and
      !! the `ocean_budget_*` helpers (and, via them, the unit tests) all
      !! read this constant, so a mutation here is caught by
      !! `test_ocean_conservation_salt_heat`.  Note: the per-stage
      !! accumulate calls in `rdb_solver.F90` and `rdb_ocean_dyn.F90`
      !! pass the literal `0.5_wp` stage weight directly; those sites are
      !! consistent by convention rather than by reference to this constant.

#ifdef RDB_ENABLE_TESTING
   public :: compute_max_cfl
   public :: compute_ice_totals
      !! Public only for the unit-test suite (`test_ocean_ice_diags`), which
      !! pins this console-side gather copy against the fills' copy.
   public :: efp_decompose_impl
   public :: compute_total_h_efp, compute_total_tracer_efp
   public :: compute_total_ke_efp, compute_ice_totals_efp
   public :: compute_total_h, compute_total_tracer, compute_total_ke
      !! Public only for the unit-test suite (PR-32,
      !! `test_ocean_console_stats_efp`), which pins `efp_decompose_impl`
      !! (the in-module `!$acc routine seq` duplicate of
      !! `rdb_efp::efp_decompose`) against the canonical procedure, and
      !! exercises the EFP reduction kernels on a device-mapped state
      !! against their FP twins (`compute_total_h`/`_tracer`/`_ke`, also
      !! exposed here for that comparison).
#endif

contains

   pure function ocean_budget_src(source_sum, stage_weight) result(src)
      !! Console `src` term for a conserved tracer: the surface (+ any other
      !! source, e.g. geothermal) contribution to the outer-step budget.
      !!
      !!   src = +RK2_STAGE_WEIGHT · source_sum · RHO_WATER
      !!
      !! `source_sum` is the already-area-weighted interior reduction of the
      !! per-cell source budget array(s) — i.e. `compute_total_tracer`
      !! output, PRE-weight and PRE-`RHO_WATER`.  Positive = added to the
      !! ocean by the source.  Sole home of the weight + ρ scaling on the
      !! source side; both the reporter and the unit test call this, so a
      !! mutation of either factor fails the test.
      real(wp), intent(in) :: source_sum
      real(wp), intent(in), optional :: stage_weight
         !! Per-outer-step weight that converts the accumulator's raw sum
         !! into the state change it must match.  Absent ⇒ `RK2_STAGE_WEIGHT`
         !! (0.5), the historical SSP-RK2 value — so every existing caller is
         !! bit-identical.  Supply `ocean_budget_stage_weight(is_pc)` rather
         !! than a literal: that function is the sole home of the mapping.
      real(wp) :: src
      real(wp) :: w
      w = RK2_STAGE_WEIGHT
      if (present(stage_weight)) w = stage_weight
      src = w*source_sum*RHO_WATER
   end function ocean_budget_src

   pure function ocean_budget_out(adv_sum, hdiff_sum, stage_weight) result(out)
      !! Console `out` term for a conserved tracer: the boundary transport
      !! out of the domain, from the horizontal-advection + horizontal-
      !! diffusion budgets.
      !!
      !!   out = −RK2_STAGE_WEIGHT · (adv_sum + hdiff_sum) · RHO_WATER
      !!
      !! `adv_sum` / `hdiff_sum` are the already-area-weighted interior
      !! reductions of the per-cell budget arrays (PRE-weight, PRE-`RHO_WATER`).
      !! The leading MINUS converts an interior loss (negative divergence
      !! integral) into a positive "left the domain" quantity.  Sole home of
      !! the weight + sign + ρ scaling on the transport side.
      real(wp), intent(in) :: adv_sum, hdiff_sum
      real(wp), intent(in), optional :: stage_weight
         !! See `ocean_budget_src`.  Absent ⇒ `RK2_STAGE_WEIGHT` (0.5).
      real(wp) :: out
      real(wp) :: w
      w = RK2_STAGE_WEIGHT
      if (present(stage_weight)) w = stage_weight
      out = -w*(adv_sum + hdiff_sum)*RHO_WATER
   end function ocean_budget_out

   pure function ocean_budget_stage_weight(is_pc) result(w)
      !! The per-outer-step weight that turns the salt/heat budget
      !! ACCUMULATORS into the state change they must account for.
      !!
      !!   ssp_rk2  (is_pc = .false.)  ->  0.5   (= `RK2_STAGE_WEIGHT`)
      !!   pred_corr  (is_pc = .true. )  ->  1.0
      !!
      !! `ms%{salt,heat}_budget_*` are `+=` accumulators filled by the tracer
      !! chain, zeroed once per OUTER step.  Under `split_scheme="ssp_rk2"`
      !! the chain runs on BOTH identical stages, so the accumulator holds
      !! twice the step's contribution and 0.5 converts it.  Under
      !! `split_scheme="pred_corr"` there is ONE prognostic tracer update — the
      !! corrector; the predictor runs `TR_MODE_NONE` and contributes
      !! nothing — so the accumulator already equals the step's contribution
      !! 1:1 and halving it reports exactly half the true source/transport.
      !!
      !! This is the same distinction `ocean_frazil_heat_src` documents for
      !! the frazil accumulator (filled once per outer step ⇒ full weight).
      !!
      !! MEASURED, 2026-09-13: with the 0.5 weight applied to a pred_corr run
      !! the console `src` came out exactly half its ssp_rk2 value on every
      !! surface-flux case (`acc_channel` -1.057E+15 -> -5.285E+14;
      !! `epbl_basin` -4.720E+16 -> -2.360E+16) while the Heat TOTAL matched
      !! to every printed digit, and the closed-budget residual went from
      !! ~1e-13 to ~5e-4. Ten shipped namelists failed the stability suite's
      !! conservation gate on that alone. Sole home of the mapping.
      logical, intent(in) :: is_pc
      real(wp) :: w
      if (is_pc) then
         w = 1.0_wp
      else
         w = RK2_STAGE_WEIGHT
      end if
   end function ocean_budget_stage_weight

   pure function ocean_heat_src_sum(surface_sum, geothermal_sum, sponge_sum) result(s)
      !! Assemble the total HEAT source integral (pre-weight, pre-ρ) that
      !! feeds `ocean_budget_src`: surface flux + geothermal bottom flux +
      !! (PR-23) the map-driven sponge's tracer-relaxation source.
      !!
      !!   s = surface_sum + geothermal_sum + sponge_sum
      !!
      !! Geothermal and the sponge are both unaccounted-for-elsewhere heat
      !! sources when enabled; folding them in keeps the Heat Error a true
      !! numerical-leak residual instead of flagging the source as a
      !! spurious leak (all three terms share the "positive into the ocean"
      !! sign convention). `sponge_sum` is `ms%heat_budget_sponge`'s
      !! area-weighted reduction — zero unless `&ocean_sponge_nml
      !! enable=.true., relax_tracers=.true.` (the legacy band sponge's
      !! tracer sink stays un-instrumented; PR-23 replaced that gate with
      !! this fold, see `ocean_budget_is_active`'s docstring). Sole home of
      !! the geothermal + sponge FOLD — the reporter and the unit tests
      !! both call it, so deleting either term fails
      !! `test_geothermal_folds_into_heat_src` /
      !! `sponge_source_closes_the_heat_budget`.
      real(wp), intent(in) :: surface_sum, geothermal_sum, sponge_sum
      real(wp) :: s
      s = surface_sum + geothermal_sum + sponge_sum
   end function ocean_heat_src_sum

   pure function ocean_salt_src_sum(surface_sum, sponge_sum) result(s)
      !! Assemble the total SALT source integral (pre-weight, pre-ρ) that
      !! feeds `ocean_budget_src`: surface flux + (PR-23) the map-driven
      !! sponge's tracer-relaxation source. Mirror of `ocean_heat_src_sum`
      !! without the geothermal term (salt has no geothermal analogue).
      !! `sponge_sum` is `ms%salt_budget_sponge`'s area-weighted reduction —
      !! zero unless `&ocean_sponge_nml enable=.true., relax_tracers=.true.`.
      real(wp), intent(in) :: surface_sum, sponge_sum
      real(wp) :: s
      s = surface_sum + sponge_sum
   end function ocean_salt_src_sum

   pure function ocean_frazil_heat_src(frazil_sum) result(src)
      !! Console `src` term for the sea-ice frazil clamp (PR 1): the heat
      !! the clamp ADDS to the ocean warming the supercooled surface layer
      !! up to T_f (the matching deficit is banked on `ice%frazil_heat`).
      !!
      !!   src = +frazil_sum · RHO_WATER      (FULL weight — no RK2 half)
      !!
      !! Unlike the per-stage accumulators behind `ocean_budget_src`, the
      !! frazil accumulator is filled ONCE per outer step on the
      !! post-RK2-average state, so its integral already equals the state
      !! change 1:1 and must NOT be halved by `RK2_STAGE_WEIGHT`.  Sole
      !! home of that full-weight scaling — the reporter and the frazil
      !! unit test both call it.
      real(wp), intent(in) :: frazil_sum
      real(wp) :: src
      src = frazil_sum*RHO_WATER
   end function ocean_frazil_heat_src

   pure function ocean_budget_is_active(tracer_idx, horiz_adv_budget_valid, &
                                        redi_with_open_edge) result(active)
      !! Whether a tracer's closed budget (`out`/`src` residual) should be
      !! reported.  `.true.` iff the tracer is registered (`tracer_idx > 0`)
      !! AND the horizontal-advection accumulator is complete AND no
      !! un-instrumented interior source is active.
      !!
      !! `horiz_adv_budget_valid = .false.` declares that some horizontal
      !! tracer transport moved `hTr` without being recorded in
      !! `*_budget_horiz_adv`, so a "closed" residual read off that
      !! accumulator would be WRONG.
      !!
      !! The windowed tracer-advection path (`dt_tracer_advect_ratio > 1`)
      !! USED to be such a producer — the fused `tracer_advect_*` kernels are
      !! skipped in `TR_MODE_ACCUMULATE`, so nothing filled the accumulator.
      !! That fall-back silently reported a live SURFACE HEAT FLUX as a 5e-5
      !! "leak" (the raw-drift column cannot subtract a source), which is what
      !! the four `acc_channel` windowed namelists were failing on.  The
      !! windowed path is now fully instrumented — `continuity_tracer_drain`'s
      !! sub-cycle and BOTH halves of the concentration hold accumulate into
      !! the same arrays (see `DRAIN_BUDGET_POST_AVERAGE_WEIGHT` in
      !! `rdb_continuity`) — so that producer is gone and the driver no longer
      !! passes this argument.  The dummy is kept for the next
      !! un-instrumented transport path (and for the fall-back unit test).
      !!
      !! `redi_with_open_edge = .true.` (optional, default `.false.`): Redi
      !! neutral-diffusion is enabled together with at least one OPEN boundary.
      !! The Redi flux crosses the open edge without being mirrored into the
      !! `out` accumulator.  Falls back to raw drift.
      !!
      !! PR-23 (real sponge) note: this gate used to carry a
      !! `sponge_relax_tracers` fall-back — the legacy band sponge's tracer
      !! sink was an un-instrumented interior source. Both sponge paths are
      !! now instrumented (the legacy band's `hTr` sink was folded into the
      !! same accumulators the map-driven path uses — see
      !! `ms%salt_budget_sponge` / `heat_budget_sponge` and
      !! `ocean_heat_src_sum` / `ocean_salt_src_sum`), so the gate is dead
      !! and has been removed; the closed budget stays active with the
      !! sponge on.
      !!
      !! In the remaining fall-back case the console reverts to raw drift —
      !! `(Q_total − Q0)/Q0` — which is byte-identical to the pre-budget
      !! behaviour (`active = .false.`).  Sole home of the gate logic — the
      !! reporter and the fallback unit test both call it.
      integer, intent(in) :: tracer_idx
      logical, intent(in) :: horiz_adv_budget_valid
      logical, intent(in), optional :: redi_with_open_edge
         !! When `.true.` the Redi flux crosses an open boundary without
         !! accounting; fall back to raw drift (v1 limitation).
      logical :: active
      logical :: redi_gate
      redi_gate = .false.
      if (present(redi_with_open_edge)) redi_gate = redi_with_open_edge
      active = (tracer_idx > 0) .and. horiz_adv_budget_valid .and. .not. redi_gate
   end function ocean_budget_is_active

   subroutine ocean_console_stats_report(this, grid, metrics, ms, t, dt, step, &
                                         horiz_adv_budget_valid, &
                                         redi_with_open_edge, &
                                         cfl_vanish_tol, heat_budget_frazil, &
                                         ice_part_size, ice_m_ice, ice_ncat, &
                                         compute_rank, reproducing_sums, &
                                         budget_stage_weight)
      !! Compute current totals + means + max-CFL and emit a MOM6-style
      !! console block via the shared `console_stats_report` formatter.
      !!
      !! COLLECTIVE: must be called by ALL compute ranks.  Per-rank local
      !! sums are reduced globally (`halo_allreduce_sum` / `_max`) before any
      !! derived quantity is formed; the shared formatter then prints on rank
      !! 0 only (`is_root`) while its NaN / CFL panic `error stop` stays
      !! collective on every rank.  Single-rank allreduce is an identity ⇒
      !! bit-identical to the serial path.
      !!
      !! Reductions run on device; tracer registry indirection is dereferenced
      !! on the host shim before each flat-impl helper.
      !!
      !! Optional `horiz_adv_budget_valid` / `redi_with_open_edge`:
      !! salt/heat closed-budget fall-back gates (see
      !! `ocean_budget_is_active`).  Absent ⇒ budget active / gate off.
      !! (PR-23 removed the `sponge_relax_tracers` fall-back — both sponge
      !! paths' tracer relaxation are now instrumented into
      !! `ms%salt_budget_sponge`/`heat_budget_sponge`, folded by
      !! `ocean_salt_src_sum`/`ocean_heat_src_sum`, so the gate was dead.)
      !!
      !! Phase-3 optional `cfl_vanish_tol`: when present (> 0) the LOCAL MaxCFL
      !! reduction excludes cells with `h_layer <= vanish_tol` before the
      !! global `allreduce_max`.  Absent ⇒ un-gated (bit-identical).
      type(console_stats_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(multilayer_state_t), intent(in) :: ms
      real(wp), intent(in) :: t, dt
      integer, intent(in) :: step
      logical, intent(in), optional :: horiz_adv_budget_valid
         !! When `.false.` the horizontal-advection accumulator is incomplete
         !! (windowed path); revert to raw drift for salt/heat.
      logical, intent(in), optional :: redi_with_open_edge
         !! When `.true.` fall back to raw drift (Redi flux at open edge not
         !! captured in the `out` accumulator).
      real(wp), intent(in), optional :: cfl_vanish_tol
         !! When present: vanish-gated MaxCFL (Phase 3). Absent ⇒ un-gated.
      real(wp), intent(in), optional :: heat_budget_frazil(:, :, :)
         !! Sea-ice frazil heat-budget accumulator (`ocean_sea_ice_t%
         !! heat_budget_frazil`, K·m per cell) — folded into `heat_src`
         !! at FULL weight (`ocean_frazil_heat_src`) so the closed heat
         !! budget still closes with ice on.  Absent (ice off — the
         !! driver passes the slot's unallocated array, which an optional
         !! dummy sees as absent) ⇒ bit-identical.
      real(wp), intent(in), optional :: ice_part_size(:, :, 0:)
         !! Sea-ice category area fractions (`ocean_sea_ice_t%part_size`,
         !! category 0 = open water).  Absent (ice off — the driver passes
         !! the slot's unallocated array, seen as absent) ⇒ no ice line ⇒
         !! bit-identical.  Lower bound 0 declared so category indexing
         !! matches the state array.
      real(wp), intent(in), optional :: ice_m_ice(:, :, :)
         !! Sea-ice mass per category (two-mode convention, `rdb_ice_state`).
      integer, intent(in), optional :: ice_ncat
         !! Category count (selects the lumped vs ITD gather mode).
      integer, intent(in) :: compute_rank
         !! Caller's compute-communicator rank (0 = print, others silent).
      logical, intent(in), optional :: reproducing_sums
         !! PR-32: `.true.` routes the `total_h`/`raw_ke`/`raw_heat`/
         !! `raw_salt`/`raw_age`/ice-area reductions and their cross-rank
         !! combine through the order-invariant EFP path (`rdb_efp` +
         !! `halo_allreduce_efp_list`, ONE collective) instead of the FP
         !! `!$acc parallel loop reduction(+:acc)` kernels + SEVEN separate
         !! `halo_allreduce_sum` calls.  Absent / `.false.` (default) ⇒ the
         !! FP path below runs VERBATIM — byte-identical console output.
         !! Salt/heat closed-budget out/src terms (section (d) below) stay
         !! on the FP path in EITHER case (documented scope reduction —
         !! see the PR-32 report); only the primary totals feeding the
         !! `Error` residual get the EFP + `efp_real_diff` treatment.
      real(wp), intent(in), optional :: budget_stage_weight
         !! Per-outer-step weight for the salt/heat budget accumulators, from
         !! `ocean_budget_stage_weight(is_pc)`.  Absent ⇒ the historical
         !! SSP-RK2 0.5 ⇒ bit-identical.  MUST be 1.0 under
         !! `&ocean_bt_nml split_scheme="pred_corr"`, whose single prognostic
         !! tracer update fills the accumulators once per step rather than
         !! twice; without it the reported `src`/`out` are exactly half and
         !! the closed-budget residual reads ~5e-4 instead of ~1e-13.

      ! (a) per-rank local sums; (b) global after allreduce; (c) derive.
      real(wp) :: total_h, raw_ke, raw_heat, raw_salt, raw_age, max_cfl
      real(wp) :: total_mass, total_ke, total_heat, total_salt
      real(wp) :: mean_S, mean_T, mean_age
      real(wp) :: tmp
      logical :: hav, redi_gate
      type(conservation_budget_t) :: bud
      real(wp) :: b_salt_surf, b_salt_adv, b_salt_hdiff, b_salt_sponge
      real(wp) :: b_heat_surf, b_heat_geo, b_heat_adv, b_heat_hdiff, b_heat_sponge
      real(wp) :: b_heat_frazil
      real(wp) :: bud_w
      logical :: ice_on, use_efp
      real(wp) :: l_wet_area, l_ci_area, l_hi_area
      real(wp) :: g_wet_area, g_ci_area, g_hi_area
      real(wp) :: mean_ci, mean_hi
      character(len=256) :: ice_line
      real(wp) :: g_mass_out, g_mass_src
      ! PR-32 EFP local/global lists — fixed NVAL=9 layout regardless of
      ! gating (ice-off slots stay zero-valued efp_t) so the collective's
      ! size never depends on a per-rank branch.
      integer, parameter :: NVAL_EFP = 10
      integer, parameter :: IX_H = 1, IX_KE = 2, IX_HEAT = 3, IX_SALT = 4, IX_AGE = 5
      integer, parameter :: IX_WET = 6, IX_CI = 7, IX_HI = 8, IX_MOUT = 9
      integer, parameter :: IX_MSRC = 10
      type(efp_t) :: efp_local(NVAL_EFP), efp_global(NVAL_EFP)
      type(efp_t) :: mass_efp_v, salt_efp_v, heat_efp_v

      use_efp = .false.
      if (present(reproducing_sums)) use_efp = reproducing_sums

      ! ---- (a) per-rank local reductions over the subdomain ---------------
      ! Physical integrals are area-weighted: `compute_total_*` fold areaT so
      ! `total_h` is Σ h·areaT [m³] (ghosts excluded ⇒ no double-count on seam).
      ! Phase-3 vanish gate applied to the LOCAL max before allreduce_max —
      ! MaxCFL is untouched by `reproducing_sums` in EITHER branch (max is
      ! already exact/order-invariant in FP; SS2.4 of the plan).
      if (present(cfl_vanish_tol)) then
         max_cfl = compute_max_cfl(ms%u_face_x_layer, ms%v_face_y_layer, &
                                   metrics%idxT, metrics%idyT, dt, grid%nghost, &
                                   ms%h_layer, cfl_vanish_tol)
      else
         max_cfl = compute_max_cfl(ms%u_face_x_layer, ms%v_face_y_layer, &
                                   metrics%idxT, metrics%idyT, dt, grid%nghost)
      end if

      ! Sea-ice conc/thickness area sums (PR ice-diags). `ice_on` is uniform
      ! across ranks (the ice arrays are allocated on every rank iff
      ! ice%enable) so gating the three extra allreduces below on it is
      ! collective-safe.  Pre-zero the global accumulators so the ice-off
      ! path's later `if (ice_on .and. g_wet_area > ...)` read is defined
      ! even if `.and.` does not short-circuit.
      g_wet_area = 0.0_wp
      g_ci_area = 0.0_wp
      g_hi_area = 0.0_wp
      ice_on = present(ice_part_size) .and. present(ice_m_ice) .and. present(ice_ncat)

      if (use_efp) then
         ! ---- EFP path: k-slab-blocked fixed-point reductions, ONE combined
         ! collective (`halo_allreduce_efp_list`) in place of the seven
         ! `halo_allreduce_sum` calls the FP branch below issues for these
         ! same quantities.
         efp_local = efp_t()
         efp_local(IX_H) = compute_total_h_efp(ms%h_layer, metrics%areaT, grid%nghost)
         efp_local(IX_KE) = compute_total_ke_efp(ms%h_layer, ms%u_face_x_layer, &
                                                 ms%v_face_y_layer, metrics%areaT, grid%nghost)
         if (ms%idx_temperature > 0) then
            efp_local(IX_HEAT) = compute_total_tracer_efp(ms%tracers(ms%idx_temperature)%hTr, &
                                                          metrics%areaT, grid%nghost)
         end if
         if (ms%idx_salinity > 0) then
            efp_local(IX_SALT) = compute_total_tracer_efp(ms%tracers(ms%idx_salinity)%hTr, &
                                                          metrics%areaT, grid%nghost)
         end if
         if (ms%idx_age > 0) then
            efp_local(IX_AGE) = compute_total_tracer_efp(ms%tracers(ms%idx_age)%hTr, &
                                                         metrics%areaT, grid%nghost)
         end if
         if (ice_on) then
            call compute_ice_totals_efp(metrics%wet_T, metrics%areaT, ice_part_size, ice_m_ice, &
                                        ice_ncat, grid%nghost, &
                                        efp_local(IX_WET), efp_local(IX_CI), efp_local(IX_HI))
         end if
         efp_local(IX_MOUT) = efp_from_real(real(ms%mass_out, real64))
         efp_local(IX_MSRC) = efp_from_real(real(ms%mass_src, real64))

         call halo_allreduce_efp_list(efp_local, efp_global, NVAL_EFP)

         total_h = real(efp_to_real(efp_global(IX_H)), wp)
         raw_ke = real(efp_to_real(efp_global(IX_KE)), wp)
         raw_heat = 0.0_wp
         if (ms%idx_temperature > 0) raw_heat = real(efp_to_real(efp_global(IX_HEAT)), wp)
         raw_salt = 0.0_wp
         if (ms%idx_salinity > 0) raw_salt = real(efp_to_real(efp_global(IX_SALT)), wp)
         raw_age = 0.0_wp
         if (ms%idx_age > 0) raw_age = real(efp_to_real(efp_global(IX_AGE)), wp)
         if (ice_on) then
            g_wet_area = real(efp_to_real(efp_global(IX_WET)), wp)
            g_ci_area = real(efp_to_real(efp_global(IX_CI)), wp)
            g_hi_area = real(efp_to_real(efp_global(IX_HI)), wp)
         end if
         g_mass_out = real(efp_to_real(efp_global(IX_MOUT)), wp)
         g_mass_src = real(efp_to_real(efp_global(IX_MSRC)), wp)

         tmp = max_cfl
         call halo_allreduce_max(tmp, max_cfl)
      else
         ! ---- FP path — byte-identical to pre-PR-32 -----------------------
         total_h = compute_total_h(ms%h_layer, metrics%areaT, grid%nghost)
         raw_ke = compute_total_ke(ms%h_layer, ms%u_face_x_layer, &
                                   ms%v_face_y_layer, metrics%areaT, grid%nghost)
         if (ms%idx_temperature > 0) then
            raw_heat = compute_total_tracer(ms%tracers(ms%idx_temperature)%hTr, &
                                            metrics%areaT, grid%nghost)
         else
            raw_heat = 0.0_wp
         end if
         if (ms%idx_salinity > 0) then
            raw_salt = compute_total_tracer(ms%tracers(ms%idx_salinity)%hTr, &
                                            metrics%areaT, grid%nghost)
         else
            raw_salt = 0.0_wp
         end if
         ! Ideal-age tracer: passive, reported whenever registered (volume-mean).
         if (ms%idx_age > 0) then
            raw_age = compute_total_tracer(ms%tracers(ms%idx_age)%hTr, &
                                           metrics%areaT, grid%nghost)
         else
            raw_age = 0.0_wp
         end if
         if (ice_on) then
            call compute_ice_totals(metrics%wet_T, metrics%areaT, ice_part_size, ice_m_ice, &
                                    ice_ncat, grid%nghost, l_wet_area, l_ci_area, l_hi_area)
         end if

         ! ---- (b) global allreduce — COLLECTIVE (all compute ranks) -------
         ! 1-rank ⇒ identities. Distinct in/out temps avoid aliasing.
         tmp = total_h
         call halo_allreduce_sum(tmp, total_h)
         tmp = raw_ke
         call halo_allreduce_sum(tmp, raw_ke)
         tmp = raw_heat
         call halo_allreduce_sum(tmp, raw_heat)
         tmp = raw_salt
         call halo_allreduce_sum(tmp, raw_salt)
         tmp = raw_age
         call halo_allreduce_sum(tmp, raw_age)
         tmp = max_cfl
         call halo_allreduce_max(tmp, max_cfl)
         if (ice_on) then
            tmp = l_wet_area
            call halo_allreduce_sum(tmp, g_wet_area)
            tmp = l_ci_area
            call halo_allreduce_sum(tmp, g_ci_area)
            tmp = l_hi_area
            call halo_allreduce_sum(tmp, g_hi_area)
         end if
         ! `g_mass_out` is computed in section (d) below on this branch (as
         ! before PR-32); the EFP branch above computes it in the combined
         ! list instead, so guard against the section-(d) allreduce running
         ! twice.
      end if

      ! ---- (c) derive conserved totals + means from the GLOBAL sums -------
      total_mass = total_h*RHO_WATER
      total_ke = raw_ke*RHO_WATER
      if (ms%idx_temperature > 0) then
         if (total_h > 0.0_wp) then
            mean_T = raw_heat/total_h
         else
            mean_T = 0.0_wp
         end if
         total_heat = raw_heat*RHO_WATER
      else
         total_heat = 0.0_wp
         mean_T = 0.0_wp
      end if
      if (ms%idx_salinity > 0) then
         if (total_h > 0.0_wp) then
            mean_S = raw_salt/total_h
         else
            mean_S = 0.0_wp
         end if
         total_salt = raw_salt*RHO_WATER
      else
         total_salt = 0.0_wp
         mean_S = 0.0_wp
      end if
      if (ms%idx_age > 0 .and. total_h > 0.0_wp) then
         mean_age = raw_age/total_h
      else
         mean_age = 0.0_wp
      end if
      mean_ci = 0.0_wp
      mean_hi = 0.0_wp
      ! `.and.` is not guaranteed to short-circuit in Fortran, so the
      ! ice-off path could read the g_* accumulators (assigned only inside
      ! the `if (ice_on)` allreduce block); they are pre-zeroed at the top
      ! of section (a) to keep this read defined.
      if (ice_on .and. g_wet_area > tiny(0.0_wp)) then
         mean_ci = g_ci_area/g_wet_area
         mean_hi = g_hi_area/g_wet_area
      end if

      ! ---- (d) closed salt/heat/mass budget — allreduced for multi-rank ---
      hav = .true.
      if (present(horiz_adv_budget_valid)) hav = horiz_adv_budget_valid
      redi_gate = .false.
      if (present(redi_with_open_edge)) redi_gate = redi_with_open_edge

      ! Mass out: cumulative open-boundary volume out (scalar accumulator).
      ! `g_mass_out` is already globally combined above on the EFP branch
      ! (packed into the same `halo_allreduce_efp_list` call as the totals);
      ! the FP branch combines it here, as before PR-32.
      if (.not. use_efp) then
         tmp = ms%mass_out
         call halo_allreduce_sum(tmp, g_mass_out)
         tmp = ms%mass_src
         call halo_allreduce_sum(tmp, g_mass_src)
      end if
      bud%mass_out = g_mass_out
      ! Tracked mass SOURCE — the mass twin of `salt_src`/`heat_src`.
      ! Zero on every path but the ice-shelf real-freshwater one
      ! (`&ocean_cavity_melt_nml freshwater="mass"`) and its
      ! `volume_compensation` sink, so the printed budget is unchanged
      ! elsewhere.  Already weighted per RK2 stage at accumulation time,
      ! exactly like `mass_out`, so no `bud_w` appears here.
      bud%mass_src = g_mass_src
      bud%mass_active = ms%mass_out_tracked

      ! Salt / heat closed budget: local per-cell integrals, allreduced, then
      ! the weight+sign+ρ scaling via the pure `ocean_budget_*` helpers.
      ! PR-23: `b_salt_sponge` / `b_heat_sponge` fold in alongside
      ! `b_heat_geo` — zero unless `&ocean_sponge_nml enable=.true.,
      ! relax_tracers=.true.` OR the legacy `sponge_relax_tracers=.true.`.
      bud_w = RK2_STAGE_WEIGHT
      if (present(budget_stage_weight)) bud_w = budget_stage_weight
      if (ocean_budget_is_active(ms%idx_salinity, hav, &
                                 redi_with_open_edge=redi_gate)) then
         b_salt_surf = compute_total_tracer(ms%salt_budget_surface, metrics%areaT, grid%nghost)
         b_salt_sponge = compute_total_tracer(ms%salt_budget_sponge, metrics%areaT, grid%nghost)
         b_salt_adv = compute_total_tracer(ms%salt_budget_horiz_adv, metrics%areaT, grid%nghost)
         b_salt_hdiff = compute_total_tracer(ms%salt_budget_hdiff, metrics%areaT, grid%nghost)
         tmp = b_salt_surf
         call halo_allreduce_sum(tmp, b_salt_surf)
         tmp = b_salt_sponge
         call halo_allreduce_sum(tmp, b_salt_sponge)
         tmp = b_salt_adv
         call halo_allreduce_sum(tmp, b_salt_adv)
         tmp = b_salt_hdiff
         call halo_allreduce_sum(tmp, b_salt_hdiff)
         bud%salt_src = ocean_budget_src(ocean_salt_src_sum(b_salt_surf, b_salt_sponge), &
                                         stage_weight=bud_w)
         bud%salt_out = ocean_budget_out(b_salt_adv, b_salt_hdiff, stage_weight=bud_w)
         bud%salt_active = .true.
      end if
      if (ocean_budget_is_active(ms%idx_temperature, hav, &
                                 redi_with_open_edge=redi_gate)) then
         b_heat_surf = compute_total_tracer(ms%heat_budget_surface, metrics%areaT, grid%nghost)
         b_heat_geo = compute_total_tracer(ms%heat_budget_geothermal, metrics%areaT, grid%nghost)
         b_heat_sponge = compute_total_tracer(ms%heat_budget_sponge, metrics%areaT, grid%nghost)
         b_heat_adv = compute_total_tracer(ms%heat_budget_horiz_adv, metrics%areaT, grid%nghost)
         b_heat_hdiff = compute_total_tracer(ms%heat_budget_hdiff, metrics%areaT, grid%nghost)
         tmp = b_heat_surf
         call halo_allreduce_sum(tmp, b_heat_surf)
         tmp = b_heat_geo
         call halo_allreduce_sum(tmp, b_heat_geo)
         tmp = b_heat_sponge
         call halo_allreduce_sum(tmp, b_heat_sponge)
         tmp = b_heat_adv
         call halo_allreduce_sum(tmp, b_heat_adv)
         tmp = b_heat_hdiff
         call halo_allreduce_sum(tmp, b_heat_hdiff)
         ! Sea-ice frazil source (PR 1) — full weight, see
         ! `ocean_frazil_heat_src`.  Absent / ice-off ⇒ adds 0.
         b_heat_frazil = 0.0_wp
         if (present(heat_budget_frazil)) then
            b_heat_frazil = compute_total_tracer(heat_budget_frazil, &
                                                 metrics%areaT, grid%nghost)
            tmp = b_heat_frazil
            call halo_allreduce_sum(tmp, b_heat_frazil)
         end if
         bud%heat_src = ocean_budget_src(ocean_heat_src_sum(b_heat_surf, b_heat_geo, &
                                                            b_heat_sponge), &
                                         stage_weight=bud_w) &
                        + ocean_frazil_heat_src(b_heat_frazil)
         bud%heat_out = ocean_budget_out(b_heat_adv, b_heat_hdiff, stage_weight=bud_w)
         bud%heat_active = .true.
      end if

      ! ---- (e) hand off to the shared MOM6-style formatter ----------------
      ! COLLECTIVE: all ranks call it (its NaN/CFL panic error-stops on every
      ! rank); `is_root` gates the printed lines to rank 0.
      !
      ! PR-32 SS2.2 fix: the EFP path builds `mass_efp`/`salt_efp`/`heat_efp`
      ! from a SINGLE scalar multiply (`total_mass = total_h*RHO_WATER`,
      ! already formed above) re-decomposed via `efp_from_real` — a lone
      ! multiplication has no accumulation-order issue, so plain FP is fine
      ! for this step; only the SUMMATION that produced `total_h` itself
      ! needed the fixed-point treatment.  `console_stats_report` then
      ! latches these as the EFP reference and forms the `Error` residual
      ! via `efp_real_diff` against it — a difference in FIXED POINT, not a
      ! double subtraction of two already-quantised totals.
      if (use_efp) then
         mass_efp_v = efp_from_real(real(total_mass, real64))
         salt_efp_v = efp_from_real(real(total_salt, real64))
         heat_efp_v = efp_from_real(real(total_heat, real64))
         call console_stats_report(this, t, step, total_mass, total_ke, &
                                   mean_S, mean_T, max_cfl, total_salt, total_heat, &
                                   has_salt=(ms%idx_salinity > 0), &
                                   has_temp=(ms%idx_temperature > 0), &
                                   mean_age=mean_age, has_age=(ms%idx_age > 0), &
                                   budget=bud, is_root=(compute_rank == 0), &
                                   mass_efp=mass_efp_v, salt_efp=salt_efp_v, heat_efp=heat_efp_v)
      else
         call console_stats_report(this, t, step, total_mass, total_ke, &
                                   mean_S, mean_T, max_cfl, total_salt, total_heat, &
                                   has_salt=(ms%idx_salinity > 0), &
                                   has_temp=(ms%idx_temperature > 0), &
                                   mean_age=mean_age, has_age=(ms%idx_age > 0), &
                                   budget=bud, is_root=(compute_rank == 0))
      end if

      ! Sea-ice line — ocean-area-weighted mean conc/thickness (PR ice-diags).
      ! Absent ice arrays (ice off) ⇒ ice_on = .false. ⇒ no line ⇒
      ! bit-identical console output.
      if (ice_on .and. compute_rank == 0) then
         write (ice_line, "('    Ice  : conc ',F7.4,'  thick ',F9.4,' m  (ocean-area mean)')") &
            mean_ci, mean_hi
         call logger%info(trim(ice_line))
      end if

      ! Halo-exchange counter line — cumulative semantic exchange counts,
      ! rank 0 only (the shared formatter owns the physics lines above).
      if (compute_rank == 0) call logger%info("    "//trim(oh_counters_format()))
   end subroutine ocean_console_stats_report

   ! -------- device-side reductions over flat arrays -----------------
   ! Registry indirection is dereferenced once on the host before each call.

   function compute_total_h(h_layer, areaT, nghost) result(total)
      !! `Σ h_layer(i,j,k)·areaT(i,j)` over PHYSICAL cells (ghosts
      !! excluded). Explicit OpenACC reduction — the `sum()` intrinsic on
      !! a present-mapped array silently runs host-side under NVHPC
      !! non-managed mode and returns the stale host shadow.
      real(wp), intent(in) :: h_layer(:, :, :)
      real(wp), intent(in) :: areaT(:, :)
      integer, intent(in) :: nghost
      real(wp) :: total
      real(wp) :: acc
      integer :: i, j, k, nx, ny, nz, i_lo, i_hi, j_lo, j_hi
      nx = min(size(h_layer, 1), size(areaT, 1))
      ny = min(size(h_layer, 2), size(areaT, 2))
      nz = size(h_layer, 3)
      i_lo = nghost + 1
      i_hi = nx - nghost
      j_lo = nghost + 1
      j_hi = ny - nghost
      acc = 0.0_wp
      !$acc parallel loop collapse(3) reduction(+:acc) present(h_layer, areaT)
      do k = 1, nz
         do j = j_lo, j_hi
            do i = i_lo, i_hi
               acc = acc + h_layer(i, j, k)*areaT(i, j)
            end do
         end do
      end do
      total = acc
   end function compute_total_h

   function compute_total_tracer(hTr, areaT, nghost) result(total)
      !! `Σ hTr(i,j,k)·areaT(i,j)` over PHYSICAL cells (ghosts excluded).
      !! Same explicit-reduction pattern as `compute_total_h`.
      real(wp), intent(in) :: hTr(:, :, :)
      real(wp), intent(in) :: areaT(:, :)
      integer, intent(in) :: nghost
      real(wp) :: total
      real(wp) :: acc
      integer :: i, j, k, nx, ny, nz, i_lo, i_hi, j_lo, j_hi
      nx = min(size(hTr, 1), size(areaT, 1))
      ny = min(size(hTr, 2), size(areaT, 2))
      nz = size(hTr, 3)
      i_lo = nghost + 1
      i_hi = nx - nghost
      j_lo = nghost + 1
      j_hi = ny - nghost
      acc = 0.0_wp
      !$acc parallel loop collapse(3) reduction(+:acc) present(hTr, areaT)
      do k = 1, nz
         do j = j_lo, j_hi
            do i = i_lo, i_hi
               acc = acc + hTr(i, j, k)*areaT(i, j)
            end do
         end do
      end do
      total = acc
   end function compute_total_tracer

   subroutine compute_ice_totals(wet_T, areaT, part_size, m_ice, ncat, nghost, &
                                 wet_area, ci_area, hi_area)
      !! Σ wet_T·areaT, Σ ci·areaT and Σ (mice/ICE_RHO_ICE)·areaT over
      !! PHYSICAL cells (ghosts excluded) — ci/mice from the two-mode
      !! per-cell gather (`ice_cell_concentration_impl` convention,
      !! inlined; `test_ocean_ice_diags` pins the fills' copy of the same
      !! math).  One pass, three `reduction(+:)` accumulators.
      real(wp), intent(in) :: wet_T(:, :), areaT(:, :)
      real(wp), intent(in) :: part_size(:, :, 0:)
      real(wp), intent(in) :: m_ice(:, :, :)
      integer, intent(in) :: ncat, nghost
      real(wp), intent(out) :: wet_area, ci_area, hi_area
      real(wp) :: acc_w, acc_c, acc_h, ci, mice
      integer :: i, j, c, nx, ny, i_lo, i_hi, j_lo, j_hi

      nx = min(size(wet_T, 1), size(areaT, 1), size(m_ice, 1))
      ny = min(size(wet_T, 2), size(areaT, 2), size(m_ice, 2))
      i_lo = nghost + 1
      i_hi = nx - nghost
      j_lo = nghost + 1
      j_hi = ny - nghost

      acc_w = 0.0_wp
      acc_c = 0.0_wp
      acc_h = 0.0_wp
      if (ncat == 1) then
         !$acc parallel loop collapse(2) reduction(+:acc_w, acc_c, acc_h) &
         !$acc&         private(ci, mice) present(wet_T, areaT, m_ice)
         do j = j_lo, j_hi
            do i = i_lo, i_hi
               acc_w = acc_w + wet_T(i, j)*areaT(i, j)
               if (wet_T(i, j) > 0.5_wp .and. m_ice(i, j, 1) > 0.0_wp) then
                  mice = m_ice(i, j, 1)
                  ci = 1.0_wp
                  acc_c = acc_c + ci*areaT(i, j)
                  acc_h = acc_h + (mice/ICE_RHO_ICE)*areaT(i, j)
               end if
            end do
         end do
      else
         !$acc parallel loop collapse(2) reduction(+:acc_w, acc_c, acc_h) &
         !$acc&         private(ci, mice, c) present(wet_T, areaT, part_size, m_ice)
         do j = j_lo, j_hi
            do i = i_lo, i_hi
               acc_w = acc_w + wet_T(i, j)*areaT(i, j)
               if (wet_T(i, j) > 0.5_wp) then
                  mice = 0.0_wp
                  ci = 0.0_wp
                  do c = 1, ncat
                     mice = mice + part_size(i, j, c)*m_ice(i, j, c)
                     ci = ci + part_size(i, j, c)
                  end do
                  ci = min(1.0_wp, ci)
                  acc_c = acc_c + ci*areaT(i, j)
                  acc_h = acc_h + (mice/ICE_RHO_ICE)*areaT(i, j)
               end if
            end do
         end do
      end if
      wet_area = acc_w
      ci_area = acc_c
      hi_area = acc_h
   end subroutine compute_ice_totals

   function compute_total_ke(h_layer, u_face, v_face, areaT, nghost) result(total)
      !! Σ 0.5·h·(u_c²+v_c²)·areaT over PHYSICAL cells (ghosts excluded),
      !! using cell-centred face averages. Direct OpenACC reduction.
      real(wp), intent(in) :: h_layer(:, :, :)
      real(wp), intent(in) :: u_face(:, :, :), v_face(:, :, :)
      real(wp), intent(in) :: areaT(:, :)
      integer, intent(in) :: nghost
      real(wp) :: total
      real(wp) :: acc, uc, vc
      integer :: i, j, k, nx, ny, nz, i_lo, i_hi, j_lo, j_hi
      nx = min(size(h_layer, 1), size(u_face, 1) - 1, size(v_face, 1), size(areaT, 1))
      ny = min(size(h_layer, 2), size(u_face, 2), size(v_face, 2) - 1, size(areaT, 2))
      nz = min(size(h_layer, 3), size(u_face, 3), size(v_face, 3))
      i_lo = nghost + 1
      i_hi = nx - nghost
      j_lo = nghost + 1
      j_hi = ny - nghost
      acc = 0.0_wp
      !$acc parallel loop collapse(3) reduction(+:acc) &
      !$acc&         private(uc, vc) present(h_layer, u_face, v_face, areaT)
      do k = 1, nz
         do j = j_lo, j_hi
            do i = i_lo, i_hi
               uc = 0.5_wp*(u_face(i, j, k) + u_face(i + 1, j, k))
               vc = 0.5_wp*(v_face(i, j, k) + v_face(i, j + 1, k))
               acc = acc + 0.5_wp*h_layer(i, j, k)*(uc*uc + vc*vc)*areaT(i, j)
            end do
         end do
      end do
      total = acc
   end function compute_total_ke

   ! -------- PR-32: EFP (order-invariant) reduction twins ------------
   ! Behind `&ocean_diag_nml reproducing_sums` (default .false. ⇒ the FP
   ! kernels above run verbatim, byte-identical).  See `rdb_efp` for the
   ! fixed-point decomposition + the cross-rank combine
   ! (`halo_allreduce_efp_list`), and `docs/CAPABILITIES_AND_LIMITATIONS.md`
   ! for the achievable guarantee.

   pure subroutine efp_decompose_impl(r, e1, e2, e3, e4, e5, e6)
      !! In-module `!$acc routine seq` duplicate of `rdb_efp::efp_decompose`
      !! -- six SCALAR outputs (not an `int64(6)` array) so the call sites
      !! below can accumulate directly into six `reduction(+:e1..e6)`
      !! clauses (OpenACC has no portable array reduction).  Duplicated
      !! rather than called from `rdb_efp` because NVHPC's device codegen
      !! does not inline a `pure !$acc routine seq` helper across a module
      !! boundary (CLAUDE.md Gotchas); `test_efp_impl_matches_canonical`
      !! (`RDB_ENABLE_TESTING`-gated) pins this copy bin-for-bit against
      !! the canonical procedure over the same magnitude table --
      !! `compute_ice_totals`'s docstring documents the identical pattern
      !! for `ice_cell_concentration_impl`.
      !!
      !! No NaN / overflow flags (unlike the canonical `efp_decompose`):
      !! console summands here are physical products (`h*areaT`,
      !! `hTr*areaT`, kinetic energy) already covered by the console's own
      !! NaN panic (`rdb_console_stats.F90`'s `panic_on_nan`); this
      !! duplicate exists purely for the well-behaved-value bin arithmetic
      !! the pinning test exercises.
      !$acc routine seq
      real(real64), intent(in) :: r
      integer(int64), intent(out) :: e1, e2, e3, e4, e5, e6
      real(real64) :: rs, s
      real(real64), parameter :: PR1 = 2.0_real64**(2*EFP_PREC_WIDTH)
      real(real64), parameter :: PR2 = 2.0_real64**(1*EFP_PREC_WIDTH)
      real(real64), parameter :: PR3 = 1.0_real64
      real(real64), parameter :: PR4 = 2.0_real64**(-1*EFP_PREC_WIDTH)
      real(real64), parameter :: PR5 = 2.0_real64**(-2*EFP_PREC_WIDTH)
      real(real64), parameter :: PR6 = 2.0_real64**(-3*EFP_PREC_WIDTH)
      real(real64), parameter :: IPR1 = 1.0_real64/PR1
      real(real64), parameter :: IPR2 = 1.0_real64/PR2
      real(real64), parameter :: IPR3 = 1.0_real64/PR3
      real(real64), parameter :: IPR4 = 1.0_real64/PR4
      real(real64), parameter :: IPR5 = 1.0_real64/PR5
      real(real64), parameter :: IPR6 = 1.0_real64/PR6

      s = 1.0_real64
      rs = r
      if (rs < 0.0_real64) then
         s = -1.0_real64
         rs = -rs
      end if

      e1 = int(s*aint(rs*IPR1), int64)
      rs = rs - real(abs(e1), real64)*PR1
      e2 = int(s*aint(rs*IPR2), int64)
      rs = rs - real(abs(e2), real64)*PR2
      e3 = int(s*aint(rs*IPR3), int64)
      rs = rs - real(abs(e3), real64)*PR3
      e4 = int(s*aint(rs*IPR4), int64)
      rs = rs - real(abs(e4), real64)*PR4
      e5 = int(s*aint(rs*IPR5), int64)
      rs = rs - real(abs(e5), real64)*PR5
      e6 = int(s*aint(rs*IPR6), int64)
   end subroutine efp_decompose_impl

   subroutine efp_summands_guard(nx, ny, name)
      !! Fail-loud (never silent) guard: a k-slab's physical cell count
      !! must not exceed `EFP_MAX_SUMMANDS`, else a bin could overflow
      !! `int64` before the next `efp_carry`.  1.34e8 is an 11500^2
      !! single-rank layer -- unreachable today, but the check costs one
      !! comparison at status cadence (CLAUDE.md: an unchecked bound "is a
      !! silent-corruption path exactly like the NZ_STACK_MAX one").
      integer, intent(in) :: nx, ny
      character(len=*), intent(in) :: name
      if (int(nx, int64)*int(ny, int64) > EFP_MAX_SUMMANDS) then
         call logger%error("============================================")
         call logger%error("[panic] "//name//": k-slab cell count exceeds EFP_MAX_SUMMANDS")
         call logger%error("============================================")
         error stop "EFP k-slab reduction: EFP_MAX_SUMMANDS exceeded"
      end if
   end subroutine efp_summands_guard

   function compute_total_h_efp(h_layer, areaT, nghost) result(total)
      !! EFP twin of `compute_total_h`: order-invariant fixed-point
      !! Sigma h_layer*areaT over PHYSICAL cells.  K-SLAB BLOCKED: a host
      !! loop over `k`, one device `reduction(+:e1..e6)` per slab, then a
      !! host-side `efp_carry` combining the slab into the running total
      !! -- keeps each device reduction block within `EFP_MAX_SUMMANDS`
      !! (SS3.3/SS6.3 of the plan; MOM6's i/j block-partition arithmetic is
      !! NOT ported -- the k-slab is simpler and sufficient at
      !! `EFP_PREC_WIDTH = 36`).  Ghost exclusion + extent clamping copied
      !! VERBATIM from `compute_total_h` -- a divergence here would
      !! silently change what is summed between the FP and EFP paths.
      real(wp), intent(in) :: h_layer(:, :, :)
      real(wp), intent(in) :: areaT(:, :)
      integer, intent(in) :: nghost
      type(efp_t) :: total
      integer :: i, j, k, nx, ny, nz, i_lo, i_hi, j_lo, j_hi
      integer(int64) :: e1, e2, e3, e4, e5, e6
      integer(int64) :: d1, d2, d3, d4, d5, d6
      integer(int64) :: slab_e(EFP_DIGITS)
      real(real64) :: val

      nx = min(size(h_layer, 1), size(areaT, 1))
      ny = min(size(h_layer, 2), size(areaT, 2))
      nz = size(h_layer, 3)
      i_lo = nghost + 1
      i_hi = nx - nghost
      j_lo = nghost + 1
      j_hi = ny - nghost
      call efp_summands_guard(i_hi - i_lo + 1, j_hi - j_lo + 1, "compute_total_h_efp")

      total%v = 0_int64
      do k = 1, nz
         e1 = 0_int64
         e2 = 0_int64
         e3 = 0_int64
         e4 = 0_int64
         e5 = 0_int64
         e6 = 0_int64
         !$acc parallel loop collapse(2) reduction(+:e1,e2,e3,e4,e5,e6) &
         !$acc&         private(val, d1, d2, d3, d4, d5, d6) present(h_layer, areaT)
         do j = j_lo, j_hi
            do i = i_lo, i_hi
               val = real(h_layer(i, j, k), real64)*real(areaT(i, j), real64)
               call efp_decompose_impl(val, d1, d2, d3, d4, d5, d6)
               e1 = e1 + d1
               e2 = e2 + d2
               e3 = e3 + d3
               e4 = e4 + d4
               e5 = e5 + d5
               e6 = e6 + d6
            end do
         end do
         slab_e = [e1, e2, e3, e4, e5, e6]
         call efp_carry(slab_e)
         total%v = total%v + slab_e
         call efp_carry(total%v)
      end do
   end function compute_total_h_efp

   function compute_total_tracer_efp(hTr, areaT, nghost) result(total)
      !! EFP twin of `compute_total_tracer`.  See `compute_total_h_efp`
      !! for the k-slab blocking design; ghost exclusion + extent clamping
      !! copied verbatim from `compute_total_tracer`.
      real(wp), intent(in) :: hTr(:, :, :)
      real(wp), intent(in) :: areaT(:, :)
      integer, intent(in) :: nghost
      type(efp_t) :: total
      integer :: i, j, k, nx, ny, nz, i_lo, i_hi, j_lo, j_hi
      integer(int64) :: e1, e2, e3, e4, e5, e6
      integer(int64) :: d1, d2, d3, d4, d5, d6
      integer(int64) :: slab_e(EFP_DIGITS)
      real(real64) :: val

      nx = min(size(hTr, 1), size(areaT, 1))
      ny = min(size(hTr, 2), size(areaT, 2))
      nz = size(hTr, 3)
      i_lo = nghost + 1
      i_hi = nx - nghost
      j_lo = nghost + 1
      j_hi = ny - nghost
      call efp_summands_guard(i_hi - i_lo + 1, j_hi - j_lo + 1, "compute_total_tracer_efp")

      total%v = 0_int64
      do k = 1, nz
         e1 = 0_int64
         e2 = 0_int64
         e3 = 0_int64
         e4 = 0_int64
         e5 = 0_int64
         e6 = 0_int64
         !$acc parallel loop collapse(2) reduction(+:e1,e2,e3,e4,e5,e6) &
         !$acc&         private(val, d1, d2, d3, d4, d5, d6) present(hTr, areaT)
         do j = j_lo, j_hi
            do i = i_lo, i_hi
               val = real(hTr(i, j, k), real64)*real(areaT(i, j), real64)
               call efp_decompose_impl(val, d1, d2, d3, d4, d5, d6)
               e1 = e1 + d1
               e2 = e2 + d2
               e3 = e3 + d3
               e4 = e4 + d4
               e5 = e5 + d5
               e6 = e6 + d6
            end do
         end do
         slab_e = [e1, e2, e3, e4, e5, e6]
         call efp_carry(slab_e)
         total%v = total%v + slab_e
         call efp_carry(total%v)
      end do
   end function compute_total_tracer_efp

   function compute_total_ke_efp(h_layer, u_face, v_face, areaT, nghost) result(total)
      !! EFP twin of `compute_total_ke`.  See `compute_total_h_efp` for the
      !! k-slab blocking design; face-averaging + extent clamping copied
      !! verbatim from `compute_total_ke`.
      real(wp), intent(in) :: h_layer(:, :, :)
      real(wp), intent(in) :: u_face(:, :, :), v_face(:, :, :)
      real(wp), intent(in) :: areaT(:, :)
      integer, intent(in) :: nghost
      type(efp_t) :: total
      integer :: i, j, k, nx, ny, nz, i_lo, i_hi, j_lo, j_hi
      integer(int64) :: e1, e2, e3, e4, e5, e6
      integer(int64) :: d1, d2, d3, d4, d5, d6
      integer(int64) :: slab_e(EFP_DIGITS)
      real(real64) :: val
      real(wp) :: uc, vc

      nx = min(size(h_layer, 1), size(u_face, 1) - 1, size(v_face, 1), size(areaT, 1))
      ny = min(size(h_layer, 2), size(u_face, 2), size(v_face, 2) - 1, size(areaT, 2))
      nz = min(size(h_layer, 3), size(u_face, 3), size(v_face, 3))
      i_lo = nghost + 1
      i_hi = nx - nghost
      j_lo = nghost + 1
      j_hi = ny - nghost
      call efp_summands_guard(i_hi - i_lo + 1, j_hi - j_lo + 1, "compute_total_ke_efp")

      total%v = 0_int64
      do k = 1, nz
         e1 = 0_int64
         e2 = 0_int64
         e3 = 0_int64
         e4 = 0_int64
         e5 = 0_int64
         e6 = 0_int64
         !$acc parallel loop collapse(2) reduction(+:e1,e2,e3,e4,e5,e6) &
         !$acc&         private(val, uc, vc, d1, d2, d3, d4, d5, d6) &
         !$acc&         present(h_layer, u_face, v_face, areaT)
         do j = j_lo, j_hi
            do i = i_lo, i_hi
               uc = 0.5_wp*(u_face(i, j, k) + u_face(i + 1, j, k))
               vc = 0.5_wp*(v_face(i, j, k) + v_face(i, j + 1, k))
               val = real(0.5_wp*h_layer(i, j, k)*(uc*uc + vc*vc)*areaT(i, j), real64)
               call efp_decompose_impl(val, d1, d2, d3, d4, d5, d6)
               e1 = e1 + d1
               e2 = e2 + d2
               e3 = e3 + d3
               e4 = e4 + d4
               e5 = e5 + d5
               e6 = e6 + d6
            end do
         end do
         slab_e = [e1, e2, e3, e4, e5, e6]
         call efp_carry(slab_e)
         total%v = total%v + slab_e
         call efp_carry(total%v)
      end do
   end function compute_total_ke_efp

   subroutine compute_ice_totals_efp(wet_T, areaT, part_size, m_ice, ncat, nghost, &
                                     wet_area_efp, ci_area_efp, hi_area_efp)
      !! EFP twin of `compute_ice_totals`.  2D-only (no k-slab blocking
      !! needed -- one "slab"), so a single `efp_summands_guard` call
      !! suffices.  Three SEPARATE single-pass reductions (one per
      !! accumulator) rather than one 18-scalar combined reduction clause
      !! -- ice diagnostics are a status-cadence cold path (SS11.11: "the
      !! EFP path costs ~nz times more kernel launches... unmeasurable at
      !! status cadence"), and three simple reductions are easier to keep
      !! correct than one with 18 live accumulators.  Gather logic +
      !! extent clamping copied verbatim from `compute_ice_totals`.
      real(wp), intent(in) :: wet_T(:, :), areaT(:, :)
      real(wp), intent(in) :: part_size(:, :, 0:)
      real(wp), intent(in) :: m_ice(:, :, :)
      integer, intent(in) :: ncat, nghost
      type(efp_t), intent(out) :: wet_area_efp, ci_area_efp, hi_area_efp
      integer :: i, j, c, nx, ny, i_lo, i_hi, j_lo, j_hi
      integer(int64) :: ew1, ew2, ew3, ew4, ew5, ew6
      integer(int64) :: ec1, ec2, ec3, ec4, ec5, ec6
      integer(int64) :: eh1, eh2, eh3, eh4, eh5, eh6
      integer(int64) :: dw1, dw2, dw3, dw4, dw5, dw6
      integer(int64) :: dc1, dc2, dc3, dc4, dc5, dc6
      integer(int64) :: dh1, dh2, dh3, dh4, dh5, dh6
      real(real64) :: val_w, val_c, val_h
      real(wp) :: ci, mice

      nx = min(size(wet_T, 1), size(areaT, 1), size(m_ice, 1))
      ny = min(size(wet_T, 2), size(areaT, 2), size(m_ice, 2))
      i_lo = nghost + 1
      i_hi = nx - nghost
      j_lo = nghost + 1
      j_hi = ny - nghost
      call efp_summands_guard(i_hi - i_lo + 1, j_hi - j_lo + 1, "compute_ice_totals_efp")

      ew1 = 0_int64
      ew2 = 0_int64
      ew3 = 0_int64
      ew4 = 0_int64
      ew5 = 0_int64
      ew6 = 0_int64
      ec1 = 0_int64
      ec2 = 0_int64
      ec3 = 0_int64
      ec4 = 0_int64
      ec5 = 0_int64
      ec6 = 0_int64
      eh1 = 0_int64
      eh2 = 0_int64
      eh3 = 0_int64
      eh4 = 0_int64
      eh5 = 0_int64
      eh6 = 0_int64

      if (ncat == 1) then
         !$acc parallel loop collapse(2) &
         !$acc&    reduction(+:ew1,ew2,ew3,ew4,ew5,ew6,ec1,ec2,ec3,ec4,ec5,ec6, &
         !$acc&              eh1,eh2,eh3,eh4,eh5,eh6) &
         !$acc&    private(ci, mice, val_w, val_c, val_h, &
         !$acc&            dw1, dw2, dw3, dw4, dw5, dw6, &
         !$acc&            dc1, dc2, dc3, dc4, dc5, dc6, &
         !$acc&            dh1, dh2, dh3, dh4, dh5, dh6) present(wet_T, areaT, m_ice)
         do j = j_lo, j_hi
            do i = i_lo, i_hi
               val_w = real(wet_T(i, j)*areaT(i, j), real64)
               call efp_decompose_impl(val_w, dw1, dw2, dw3, dw4, dw5, dw6)
               ew1 = ew1 + dw1
               ew2 = ew2 + dw2
               ew3 = ew3 + dw3
               ew4 = ew4 + dw4
               ew5 = ew5 + dw5
               ew6 = ew6 + dw6
               if (wet_T(i, j) > 0.5_wp .and. m_ice(i, j, 1) > 0.0_wp) then
                  mice = m_ice(i, j, 1)
                  ci = 1.0_wp
                  val_c = real(ci*areaT(i, j), real64)
                  val_h = real((mice/ICE_RHO_ICE)*areaT(i, j), real64)
                  call efp_decompose_impl(val_c, dc1, dc2, dc3, dc4, dc5, dc6)
                  ec1 = ec1 + dc1
                  ec2 = ec2 + dc2
                  ec3 = ec3 + dc3
                  ec4 = ec4 + dc4
                  ec5 = ec5 + dc5
                  ec6 = ec6 + dc6
                  call efp_decompose_impl(val_h, dh1, dh2, dh3, dh4, dh5, dh6)
                  eh1 = eh1 + dh1
                  eh2 = eh2 + dh2
                  eh3 = eh3 + dh3
                  eh4 = eh4 + dh4
                  eh5 = eh5 + dh5
                  eh6 = eh6 + dh6
               end if
            end do
         end do
      else
         !$acc parallel loop collapse(2) &
         !$acc&    reduction(+:ew1,ew2,ew3,ew4,ew5,ew6,ec1,ec2,ec3,ec4,ec5,ec6, &
         !$acc&              eh1,eh2,eh3,eh4,eh5,eh6) &
         !$acc&    private(ci, mice, c, val_w, val_c, val_h, &
         !$acc&            dw1, dw2, dw3, dw4, dw5, dw6, &
         !$acc&            dc1, dc2, dc3, dc4, dc5, dc6, &
         !$acc&            dh1, dh2, dh3, dh4, dh5, dh6) &
         !$acc&    present(wet_T, areaT, part_size, m_ice)
         do j = j_lo, j_hi
            do i = i_lo, i_hi
               val_w = real(wet_T(i, j)*areaT(i, j), real64)
               call efp_decompose_impl(val_w, dw1, dw2, dw3, dw4, dw5, dw6)
               ew1 = ew1 + dw1
               ew2 = ew2 + dw2
               ew3 = ew3 + dw3
               ew4 = ew4 + dw4
               ew5 = ew5 + dw5
               ew6 = ew6 + dw6
               if (wet_T(i, j) > 0.5_wp) then
                  mice = 0.0_wp
                  ci = 0.0_wp
                  do c = 1, ncat
                     mice = mice + part_size(i, j, c)*m_ice(i, j, c)
                     ci = ci + part_size(i, j, c)
                  end do
                  ci = min(1.0_wp, ci)
                  val_c = real(ci*areaT(i, j), real64)
                  val_h = real((mice/ICE_RHO_ICE)*areaT(i, j), real64)
                  call efp_decompose_impl(val_c, dc1, dc2, dc3, dc4, dc5, dc6)
                  ec1 = ec1 + dc1
                  ec2 = ec2 + dc2
                  ec3 = ec3 + dc3
                  ec4 = ec4 + dc4
                  ec5 = ec5 + dc5
                  ec6 = ec6 + dc6
                  call efp_decompose_impl(val_h, dh1, dh2, dh3, dh4, dh5, dh6)
                  eh1 = eh1 + dh1
                  eh2 = eh2 + dh2
                  eh3 = eh3 + dh3
                  eh4 = eh4 + dh4
                  eh5 = eh5 + dh5
                  eh6 = eh6 + dh6
               end if
            end do
         end do
      end if

      wet_area_efp%v = [ew1, ew2, ew3, ew4, ew5, ew6]
      call efp_carry(wet_area_efp%v)
      ci_area_efp%v = [ec1, ec2, ec3, ec4, ec5, ec6]
      call efp_carry(ci_area_efp%v)
      hi_area_efp%v = [eh1, eh2, eh3, eh4, eh5, eh6]
      call efp_carry(hi_area_efp%v)
   end subroutine compute_ice_totals_efp

   pure function cfl_cell_value(u_l, u_r, v_l, v_r, idx, idy, dt) result(cfl)
      !! Per-cell advective CFL from the C-grid face-velocity pairs and metric
      !! inverses: `cfl = (|u_c|·idx + |v_c|·idy)·dt`, with `u_c`/`v_c` the
      !! face averages.  `!$acc routine seq` so it inlines into the reduction
      !! loops below — same module ⇒ NVHPC keeps it inlined (and it's a
      !! status-cadence cold path regardless).  Sole home of the CFL formula,
      !! shared by the gated + un-gated `compute_max_cfl` loops.
      !$acc routine seq
      real(wp), intent(in) :: u_l, u_r, v_l, v_r, idx, idy, dt
      real(wp) :: cfl, uc, vc
      uc = 0.5_wp*(u_l + u_r)
      vc = 0.5_wp*(v_l + v_r)
      cfl = (abs(uc)*idx + abs(vc)*idy)*dt
   end function cfl_cell_value

   function compute_max_cfl(u_face, v_face, idxT, idyT, dt, nghost, h_layer, vanish_tol) &
      result(max_cfl)
      !! max (|u_c|·dt·idxT + |v_c|·dt·idyT) over PHYSICAL cells (ghosts
      !! excluded). idxT/idyT are metric inverses (= 1/dx,1/dy on uniform).
      !!
      !! Phase-3 optional gate: when `h_layer` + `vanish_tol` are BOTH
      !! present, cells where `h_layer(i,j,k) <= vanish_tol` are skipped.
      !! A vanished cell carries no real momentum; its face-averaged
      !! velocity spike should not trigger a CFL panic.  When absent
      !! (default) ⇒ un-gated path runs verbatim ⇒ bit-identical.
      !!
      !! Two-loop form (un-gated / gated dispatched externally before
      !! calling) preferred over referencing an absent optional inside the
      !! `!$acc parallel loop` region — keeps the parallel body clean.
      real(wp), intent(in) :: u_face(:, :, :), v_face(:, :, :)
      real(wp), intent(in) :: idxT(:, :), idyT(:, :)
      real(wp), intent(in) :: dt
      integer, intent(in) :: nghost
      real(wp), intent(in), optional :: h_layer(:, :, :)
         !! Centre-cell thickness (m). Required together with `vanish_tol`.
      real(wp), intent(in), optional :: vanish_tol
         !! Cells with h_layer <= vanish_tol are excluded. Required
         !! together with `h_layer`.
      real(wp) :: max_cfl
      real(wp) :: acc
      integer :: i, j, k, nx, ny, nz, i_lo, i_hi, j_lo, j_hi
      logical :: gate

      nx = min(size(u_face, 1) - 1, size(v_face, 1), size(idxT, 1))
      ny = min(size(u_face, 2), size(v_face, 2) - 1, size(idyT, 2))
      nz = min(size(u_face, 3), size(v_face, 3))
      i_lo = nghost + 1
      i_hi = nx - nghost
      j_lo = nghost + 1
      j_hi = ny - nghost
      gate = present(h_layer) .and. present(vanish_tol)

      if (gate) then
         ! Gated path: skip cells where h_layer <= vanish_tol.
         ! Separate loop so the compiler keeps the !$acc parallel region clean
         ! (no absent-optional reference inside the parallel body).
         acc = 0.0_wp
         !$acc parallel loop collapse(3) reduction(max:acc) &
         !$acc&         present(u_face, v_face, h_layer, idxT, idyT)
         do k = 1, nz
            do j = j_lo, j_hi
               do i = i_lo, i_hi
                  if (h_layer(i, j, k) > vanish_tol) then
                     acc = max(acc, cfl_cell_value( &
                               u_face(i, j, k), u_face(i + 1, j, k), &
                               v_face(i, j, k), v_face(i, j + 1, k), &
                               idxT(i, j), idyT(i, j), dt))
                  end if
               end do
            end do
         end do
      else
         ! Un-gated path: all layers contribute.  Byte-identical to the
         ! pre-Phase-3 implementation.
         acc = 0.0_wp
         !$acc parallel loop collapse(3) reduction(max:acc) &
         !$acc&         present(u_face, v_face, idxT, idyT)
         do k = 1, nz
            do j = j_lo, j_hi
               do i = i_lo, i_hi
                  acc = max(acc, cfl_cell_value( &
                            u_face(i, j, k), u_face(i + 1, j, k), &
                            v_face(i, j, k), v_face(i, j + 1, k), &
                            idxT(i, j), idyT(i, j), dt))
               end do
            end do
         end do
      end if
      max_cfl = acc
   end function compute_max_cfl

end module rdb_ocean_console_stats
