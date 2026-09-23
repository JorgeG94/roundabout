!! Continuity-PPM kernel state + barotropic step.
module rdb_continuity
   !! Holds scheme-variant flags and reusable workspace for the
   !! C-grid continuity-PPM (Lin & Rood) thickness-flux kernel, plus
   !! the public step routines that the split-explicit driver calls
   !! per barotropic substep.  The kernel is the primary mass-flux
   !! producer for the ocean path: it consumes face velocities +
   !! cell-centred thickness and emits per-face mass fluxes that the
   !! tracer-advection kernels then reuse under CWC for free.
   !!
   !! Phase 2b: method-of-lines PPM face reconstruction (Colella &
   !! Woodward 1984, eq 1.6 face value + eq 1.10 monotonic limiter)
   !! on the barotropic C-grid state.  Closed-wall BC.  Cells too
   !! close to a wall (< 2 cells from the boundary) fall back to
   !! first-order (h_L = h_R = h_centre); the resulting kernel
   !! preserves uniform fields bit-for-bit (the constancy-preservation
   !! property the lake-at-rest tests guard) and propagates Gaussian
   !! humps with < 5% peak diffusion over their own width.
   use rdb_constants, only: wp, NZ_STACK_MAX, H_DIV_EPS
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_barotropic_state, only: barotropic_state_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_tracer, only: TRACER_BUDGET_HEAT, TRACER_BUDGET_SALT
   use rdb_scratch_3d, only: scratch_3d_buffer_t, &
                             scratch_3d_buffer_enter_data_impl, &
                             scratch_3d_buffer_exit_data_impl
   use rdb_ocean_boundary_types, only: ocean_bc_state_t, OBC_WALL, OBC_CLAMPED
   use rdb_ocean_periodic, only: ocean_periodic_wrap_centre_3d, &
                                 ocean_periodic_wrap_face_x_3d, &
                                 ocean_periodic_wrap_face_y_3d
   use rdb_ocean_fold_apply, only: ocean_fold_wrap_centre_3d_state
   use rdb_ocean_fold, only: fold_north_centre, fold_north_u_face, fold_north_v_face
   use rdb_ocean_mle, only: ocean_mle_t, mle_fold_x, mle_fold_y
   use rdb_ocean_gm, only: ocean_gm_t, gm_fold_x, gm_fold_y
   use rdb_ocean_halo, only: ocean_halo_centre, &
                             ocean_halo_is_decomposed_x, ocean_halo_is_decomposed_y
   use rdb_profiler, only: profiler_start, profiler_stop
   ! WENO tracer-face reconstruction ladder — reused verbatim from the
   ! coastal multilayer path (Q6).  Only the pure `!$acc routine seq`
   ! swept-average face helpers + the rung-degradation predicate are
   ! imported; the coastal driver kernel is not.  Called cross-module
   ! from `do concurrent` (the helpers carry `!$acc routine seq`).
   use rdb_recon_weno, only: plm_face_swept, weno5_face_swept, &
                             weno7_face_swept, weno9_face_swept, &
                             recon_rung_for_face, &
                             TRACER_RECON_PPM
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_mem_report, only: arr_bytes
   implicit none
   private

   ! Phase 2 (6b) windowed-drain floors.  DRAIN_MIN_H guards Tr = hTr/h
   ! at vanishing thickness (concentration space); DRAIN_MIN_VOL guards
   ! the per-pass Courant denominator areaT·hprev.  Same role as MOM6's
   ! Angstrom_H / a small volume floor.
   real(wp), parameter :: DRAIN_MIN_H = 1.0e-10_wp
   real(wp), parameter :: DRAIN_MIN_VOL = 1.0e-10_wp
   ! Fixed-budget pass count for the upfront availability limiter
   ! (drain_avail_limit).  Each pass is one FCT inflow-scaling sweep; the
   ! count bounds how far the per-cell positivity constraint diffuses
   ! across the C-grid stencil.  8 passes hold vol > 0 on the kitchensink.
   integer, parameter :: AVAIL_LIMIT_PASS = 8
   ! Windowed-advect budget bookkeeping weight.  The shared per-cell tracer
   ! budget accumulators (`ms%heat_budget_horiz_adv` / `salt_budget_horiz_adv`)
   ! carry the RAW SUM over both RK2 stages, because `rk2_average` halves the
   ! state change and the console reporter re-applies that same
   ! `RK2_STAGE_WEIGHT = 0.5` (see `ocean_budget_out`).  Every windowed-drain
   ! write to `hTr` that happens OUTSIDE a stage -- the end-of-window
   ! concentration un-hold and the drain sub-cycle itself, both fired once per
   ! outer step on the already-averaged state -- therefore has to be
   ! pre-multiplied by 2 so the reporter's 0.5 recovers it 1:1.  In-stage
   ! writes (the per-stage concentration hold) use weight 1, like every other
   ! in-stage accumulator.  Mutating this breaks
   ! `windowed_drain_with_surface_flux` (test_ocean_conservation_salt_heat).
   real(wp), parameter :: DRAIN_BUDGET_POST_AVERAGE_WEIGHT = 2.0_wp
   real(wp), parameter :: DRAIN_BUDGET_IN_STAGE_WEIGHT = 1.0_wp

   public :: continuity_t

   ! Production API — imported by rdb_ocean_dyn / rdb_ocean_state.
   ! Barotropic continuity flux compute/apply (fast barotropic loop).
   public :: continuity_compute_fluxes_barotropic
   public :: continuity_apply_fluxes_barotropic
   ! Production entry point: interleaved continuity+tracer directional (Lie)
   ! split, CWC-consistent.  Driven by rdb_ocean_dyn (split RK2).
   public :: continuity_tracer_step_split

   ! Sea-ice PR 4b: the five SIS2-equivalent PPM stencil helpers
   ! (mirror-at-land, swept-volume face flux, van-Leer slope, CW84 cell
   ! limiter, positivity limiter) promoted from test-only/private to
   ! production public — `rdb_ice_transport` reuses them verbatim
   ! (do NOT duplicate their bodies; SPEC_ice-pr4b-transport.md §2) to
   ! reconstruct the category-SUMMED ice/snow mass PPM parabola exactly
   ! as `continuity_zonal_flux`/`continuity_meridional_flux` do for the
   ! ocean layers.  All five are `!$acc routine seq`, safe to call from a
   ! `do concurrent` in any module.
   public :: ppm_mirror_h
   public :: volcfl_face
   public :: ppm_limited_slope
   public :: ppm_cell_limiter
   public :: ppm_limit_pos

   ! Phase 2 (6b) windowed horizontal tracer advection: accumulate the
   ! per-stage mass fluxes (TR_MODE_ACCUMULATE) over the DT_TRACER_ADVECT
   ! window, then drain them once at the window boundary.
   public :: continuity_tracer_drain
   integer, parameter, public :: TR_MODE_ADVECT = 0
      !! Every-step fused mode: advance h AND advect tracers each call
      !! (the historical default; ratio = 1 bypass).
   integer, parameter, public :: TR_MODE_ACCUMULATE = 1
   integer, parameter, public :: TR_MODE_NONE = 2
      !! h-only continuity: skip tracer advection AND window accumulation
      !! entirely.  Used by the pred_corr predictor (SPEC §2 P9 — MOM6's
      !! predictor continuity advances hp but never touches tracers; the
      !! predictor state is discarded except for u_av / h_av).
      !! Windowed mode (ratio > 1): advance h, accumulate
      !! 0.5·mass_flux·dt into uhtr/vhtr, SKIP the per-step tracer advect
      !! (hTr held frozen until the boundary drain).

#ifdef RDB_ENABLE_TESTING
   ! ----------------------------------------------------------------------
   ! Exposed to the unit-test suite ONLY (built with -DRDB_ENABLE_TESTING;
   ! private in production builds).  None of these are imported by production
   ! code — they are reached intra-module on the hot path, or exist purely as
   ! test oracles — so guarding the publicity keeps the production API minimal
   ! while leaving the implementations untouched.
   ! ----------------------------------------------------------------------
   ! Directionally-split (MOM6-style) building blocks, called INTERNALLY by
   ! continuity_tracer_step_split; exposed for isolated unit testing.
   ! (ppm_mirror_h is now unconditionally public above — sea-ice PR 4b.)
   public :: continuity_zonal_flux
   public :: continuity_meridional_flux
   public :: continuity_apply_zonal
   public :: continuity_apply_meridional
   ! Positive-definite per-donor outflux limiters — exposed so the u_cor/v_cor
   ! re-matching test can snapshot pre/post flux + velocity across the limiter.
   public :: pd_limit_zonal_impl
   public :: pd_limit_meridional_impl
   public :: tracer_advect_zonal
   public :: tracer_advect_meridional
   ! Per-direction tracer-advect IMPLs — exposed so the conservation test can
   ! drive the budget-fill WITH vs WITHOUT the optional budget_adv arg on
   ! identical inputs (proving the fill loop never perturbs hTr — the
   ! present(budget_adv)-ABSENT path, which the wrapper-level S/T tests can't
   ! contrast because they always pass budget_adv for S and T).
   public :: tracer_advect_zonal_one_impl
   public :: tracer_advect_meridional_one_impl
   ! (ppm_limit_pos is now unconditionally public above — sea-ice PR 4b —
   ! load-bearing PPM positivity limiter used internally on the hot path.)
   ! No production caller — the unsplit (non-directionally-split) compute/
   ! apply/advect reference path, the oracle the split production path is
   ! validated against; continuity_step_split is the continuity-only wrapper
   ! (production runs the tracer-interleaved continuity_tracer_step_split).
   public :: continuity_compute_fluxes
   public :: continuity_apply_fluxes
   public :: tracer_advect
   public :: continuity_step_split
   ! Phase 2 (6b) windowed-drain internals — exposed for the swept-flux
   ! oracle + reconstruction unit tests.
   public :: drain_swept_flux_x
   public :: drain_swept_flux_y
   public :: drain_parabola_x
   public :: drain_parabola_y
   public :: drain_reconstruct_hprev
   ! Q6 WENO drain internals — exposed for the reconstruction unit tests.
   public :: drain_fill_conc
   public :: drain_swept_flux_x_weno
   public :: drain_swept_flux_y_weno
#endif

   type :: continuity_t
      logical :: is_init = .false.
         !! True between `init` and `destroy`.  Prefer this to
         !! `allocated(...)` — tracks GPU device attachment too.
      real(wp) :: angstrom_h = 0.0_wp
         !! Phase-1 Lagrangian minimum-thickness floor (m). Set from
         !! cfg%ocean%isopycnal%angstrom_h at setup, but only PASSED to the
         !! h-update kernels when the active vcoord is VCOORD_LAGRANGIAN
         !! (gated at the dyn call site). 0 ⇒ off ⇒ bit-identical.
         !! R7: floor lifts h but NOT hTr, so Tr=hTr/h shifts on a floored
         !! layer. Harmless for the adiabatic isopycnal config; thermo-on
         !! isopycnal correctness is OUT OF SCOPE for v1.
      logical :: conservative_floor = .false.
         !! When `.true.` (set from cfg%ocean%isopycnal%conservative_floor at
         !! setup; requires angstrom_h > 0 + VCOORD_LAGRANGIAN) the injecting
         !! `max(h, angstrom_h)` floor is skipped in the h-update and replaced
         !! by the conservative per-column borrow in `rdb_ocean_min_thickness`,
         !! invoked from the dyn continuity site. Uses `mt_h_new` as scratch.
         !! Default `.false.` ⇒ legacy injecting floor ⇒ bit-identical.
      real(wp) :: h_min = 1.0e-6_wp
         !! Lower clip on cell-centred thickness during the update.
         !! Also used as the floor for `ppm_limit_pos` when
         !! `use_ppm_limit_pos = .true.` — same semantic role as
         !! MOM6's `GV%Angstrom_H` for the vanishing-layer limiter.
      real(wp) :: cfl_max = 0.5_wp
         !! Soft cap on per-face CFL before falling back to upwind.
      logical :: renorm_legacy_single_step = .false.
         !! Use the pre-Newton single-linear-step uhbt renormalisation
         !! (donors picked at the UNCORRECTED velocity, no CFL bracket,
         !! no iteration) — wet/dry composition, set from
         !! `cfg%ocean%wetdry%enable` at setup.  The Newton form's donor
         !! RE-PICK + CFL bracket are what stabilise Lagrangian grounding
         !! (LAGRANGIAN_PGF_BUG.md §0), but under wet/dry they interact
         !! with the drying-front face gating and drive a drying column's
         !! `h_layer` negative (`test_ocean_wetdry_driver`); wet/dry owns
         !! its positivity via `ppm_limit_pos` + the BT limiter and ran
         !! validated on the single-step form, so it keeps it
         !! (bit-identical there).
      logical :: renorm_consistent_flux = .true.
         !! `&ocean_continuity_nml renorm_consistent_flux`.  When `.true.`
         !! the `uhbt`/`vhbt` renormalisation evaluates a layer whose upwind
         !! donor FLIPS under the correction as `(u0 + du)·h_face(new
         !! donor)`, i.e. as the flux of its corrected velocity, and
         !! brackets the Newton solve with bisection (MOM6
         !! `zonal_flux_adjust`).  The historical model
         !! `flux0 + du·h_face(new donor)` keeps the OLD donor's `u0·h_old`
         !! and is DISCONTINUOUS at the flip, by `u0·(h_new − h_old)·w`:
         !! whenever `uhbt` falls in that gap (a face where the corrected
         !! velocity must change sign across a thickness jump — a sigma
         !! layer over a bathymetric step, where `h_old ≠ h_new` after the
         !! PPM limiter flattens both edges), Newton has NO root, cycles
         !! for `RENORM_MAXIT` iterations and hands continuity a layer
         !! transport of the wrong SIGN.  The layer `η` then departs from
         !! the barotropic `η_end` by O(η) at the step every such step —
         !! a spurious η dipole that pumps the (undamped) barotropic
         !! grid-scale mode.  Default `.true.` (MOM6 behaviour); faces
         !! where no donor flips are bit-identical to the historical
         !! model, which `.false.` restores.
      logical :: use_ppm_limit_pos = .false.
         !! MOM6 `PPM_limit_pos` analogue.  When `.true.`, the PPM
         !! face-thickness reconstruction in continuity adds a
         !! positivity-preserving limiter that shrinks h_left /
         !! h_right toward h_centre whenever the parabolic fit would
         !! produce an interior minimum below `h_min`.  At the
         !! `h_centre ≤ h_min` limit the reconstruction collapses
         !! to a constant (= upwind for that cell), bounding the
         !! mass flux by the actual layer thickness.  Default off
         !! keeps the pre-knob behaviour bit-identical.  Driver writes
         !! from `cfg%ocean%continuity%ppm_limit_pos`.
      logical :: vol_cfl = .false.
         !! MOM6 `vol_CFL` analogue.  When `.false.` (default) the
         !! continuity flux uses the PPM downwind EDGE value (the
         !! CFL→0 limit), bit-identical to the pre-knob behaviour.
         !! When `.true.` the donor-side face thickness is the
         !! swept-volume integral of the reconstructed parabola
         !! (`volcfl_face`), adding the missing O(CFL) term.  Fixes
         !! the dt-sensitive near-bed residual at steep shelf breaks.
         !! Driver writes from `cfg%ocean%continuity%vol_cfl`.
      logical :: positive_definite = .false.
         !! Positive-definite split continuity master switch.  When `.true.`
         !! the split layer continuity (`continuity_tracer_step_split` via
         !! `continuity_zonal_flux`/`continuity_meridional_flux`) keeps every
         !! layer `>= h_lim`: P1 floors the PPM reconstruction edges at
         !! `2·h_lim` (MOM6's positive-definite reconstruction); P2 (later)
         !! scales down per-donor outfluxes so no mass is created.  Set from
         !! `cfg%ocean%continuity%positive_definite` at setup.  Host-side
         !! control scalar — read into a local before each DC loop (never
         !! dereferenced inside a device loop).  Default `.false.` ⇒ untaken
         !! branches only ⇒ bit-identical.
      real(wp) :: h_lim = 0.0_wp
         !! Positive-definite thickness floor (m).  Derived at setup:
         !! `cfg%ocean%isopycnal%angstrom_h` on VCOORD_LAGRANGIAN, else 0.
         !! Only consumed under `positive_definite = .true.`; at `h_lim = 0`
         !! (every non-Lagrangian coord) the P1 `max(edge, 2·h_lim)` floor is
         !! inert (edges are already `>= 0`), but the P2 outflux limiter still
         !! engages to keep every layer `>= 0` (avail = max(h, 0)).
      logical :: hTr_holds_conc = .false.
         !! True while the frozen tracer content `hTr` has been re-weighted
         !! onto the CURRENT `h_layer` so that the concentration
         !! `T = hTr/h_layer` is the (invariant) window-start value.  Set by
         !! `continuity_tracer_step_split` in TR_MODE_ACCUMULATE, cleared by
         !! `continuity_tracer_drain` once it has converted `hTr` back onto
         !! the reconstructed window-start thickness `hprev`.
         !!
         !! It exists so the drain's input contract stays EXPLICIT: `.false.`
         !! means "hTr is content on the window-start thickness" (the
         !! white-box unit tests that seed the drain directly, and the
         !! historical behaviour); `.true.` means "hTr is content on the
         !! current thickness".  Host-only — never read inside a kernel.
      logical :: windowed_advection = .false.
         !! Gates the ALLOCATION of the Phase-2 windowed-advection state (13
         !! 3D arrays: `uhtr`/`vhtr` accumulators + the (6b) drain
         !! workspace, ~4.2 GB at 1000x800x50) — consumed only when
         !! `dt_tracer_advect_ratio > 1`.  Latched from cfg BEFORE
         !! `init(grid)` by `ocean_state_init_from_config` (the same
         !! conditional-allocation contract as the default-off closures);
         !! the `enter_data`/`exit_data`/drain paths already guard on
         !! `allocated()`.
         !!
         !! Default `.false.`: the workspace is opt-in, so a `ct%init(...)`
         !! on any path that does NOT run through
         !! `ocean_state_init_from_config` costs nothing.  (It used to
         !! default `.true.` for the convenience of direct test call sites —
         !! which meant every such path silently paid the multi-GB
         !! allocation whether or not the drain could ever run.  Those call
         !! sites now set the flag explicitly before `init`.)
      integer :: n_limited_step = 0
         !! P2 diagnostic counter: number of interior faces whose mass flux
         !! the positive-definite limiter scaled (θ_donor < 1) this outer-step
         !! call, summed over the zonal + meridional passes.  Host-side scalar
         !! (the `!$acc parallel loop reduction` returns to the host); zeroed
         !! at the top of every `continuity_tracer_step_split` call.  P3 drains
         !! it to the console stats line.  0 ⇒ no limiting fired (the healthy
         !! case; MOM6-style graceful degradation must be loud).
      integer(int64) :: n_limited_total = 0_int64
         !! P3 running total of `n_limited_step` across every
         !! `continuity_tracer_step_split` call (accumulated once per call, at
         !! the end — one add per split call, mirroring `dyn%ntrunc_total`).
         !! The driver drains its per-report DELTA to the console next to the
         !! CFL-truncation line.  Grows ONLY when `positive_definite = .true.`
         !! (else the passes are skipped), so a nonzero total is itself the
         !! signal the limiter is active.  Loud-by-design caveat: a
         !! `uniform_z` isopycnal (VCOORD_LAGRANGIAN) stack shows a permanently
         !! LARGE, steadily-growing count because ~most layers sit AT the floor
         !! and θ=0 correctly freezes their (massless) outflow every call —
         !! that is not a pathology.  A healthy zstar/sigma run sits at 0.
         !! int64 (not int32 like `dyn%ntrunc_total`) precisely because that
         !! floored-stack case can accumulate > 2·10⁹ over a long run.
      integer :: tracer_recon = TRACER_RECON_PPM
         !! Face-reconstruction scheme for the WINDOWED tracer-advection
         !! drain (Q6; `dt_tracer_advect_ratio > 1` path only).  0 = CW-PPM
         !! (default, bit-identical), 1/2/3 = WENO5/7/9-Z swept-average
         !! (rung-adaptive, degrades near land/walls).  Set at configure
         !! from `&ocean_vmix_nml tracer_recon`.  Host-side control knob:
         !! the drain reads it to pick a face kernel — never dereferenced
         !! inside a device loop, so it rides the struct copyin.
         !! The every-step (ratio = 1) advect path is unaffected (stays
         !! CW-PPM) — WENO is a drain-only scheme in this wave.

      ! ---- Face-reconstruction workspace ----
      ! Phase 2b sizes these for the barotropic case (n3 = 1).  When
      ! Phase 5 adds layers, the multilayer kernel re-`init`s them
      ! at n3 = nz_ml.  Naming convention:
      !
      !   h_face_left_x(i, j, k)  — value AT east face i extrapolated
      !                             from the LEFT-side cell (i-1, j, k),
      !                             i.e. h_R of cell i-1.
      !   h_face_right_x(i, j, k) — value AT east face i extrapolated
      !                             from the RIGHT-side cell (i, j, k),
      !                             i.e. h_L of cell i.
      !
      ! Upwind: the kernel picks left if u >= 0 (left cell donates),
      ! right if u < 0.
      type(scratch_3d_buffer_t) :: h_face_left_x
         !! Left-state thickness at east faces.
      type(scratch_3d_buffer_t) :: h_face_right_x
         !! Right-state thickness at east faces.
      type(scratch_3d_buffer_t) :: h_face_left_y
         !! Left-state thickness at north faces.
      type(scratch_3d_buffer_t) :: h_face_right_y
         !! Right-state thickness at north faces.
      type(scratch_3d_buffer_t) :: mt_h_new
         !! Cell-centred `(nx,ny,nz)` scratch for the conservative
         !! minimum-thickness borrow (`conservative_floor`).  Holds the
         !! floor-only target thickness field between the h-update and the
         !! h_layer overwrite.  Unused (but allocated + mapped) when the knob
         !! is off — the DC kernels never touch it in that case.
      type(scratch_3d_buffer_t) :: mt_grounded
         !! `(nx,ny,1)` grounded-column mask (1.0 = any layer below floor)
         !! for the borrow's early-exit restructure: built in one coalesced
         !! pass, consumed by every borrow kernel in place of per-face
         !! column re-scans.  Same lifetime/mapping as `mt_h_new`.
      type(scratch_3d_buffer_t) :: pd_theta
         !! `(nx,ny,nz)` cell-centred per-donor availability factor θ(i,j,k)
         !! for the P2 positive-definite outflux limiter.  Built once per
         !! direction pass (over the FULL range incl. ghosts, so any cell that
         !! can donate to a swept face has a current θ), then each interior
         !! face is scaled by its upwind donor's θ.  Allocated + device-mapped
         !! UNCONDITIONALLY (like `mt_h_new`/`mt_grounded`) — the enter_data
         !! contract has no config visibility, and at `nx·ny·nz` reals the
         !! footprint matches `mt_h_new` already carried; unused (but present)
         !! when `positive_definite = .false.`.

      ! ---- Phase 2 flux-accumulator slots (DT_TRACER_ADVECT) ----
      ! Persistent device-resident accumulators of the per-stage C-grid
      ! face transports (`mass_flux_x/y_layer`).  When
      ! `dt_tracer_advect_ratio > 1` the windowed horizontal tracer
      ! advect (Phase 6b) drains these over the accumulation window;
      ! reconstruction `hprev = areaT·h_end + div(uhtr)` closes
      ! continuity by construction.  Allocated-but-unused at the default
      ! ratio = 1 (the every-step path bypasses them) — they are still
      ! mapped in enter_data so a ratio>1 run finds them present.
      ! Shapes mirror mass_flux_x/y_layer: east-face (nx+1,ny,nz),
      ! north-face (nx,ny+1,nz).
      real(wp), allocatable :: uhtr(:, :, :)
         !! Accumulated zonal face transport (m^3, area-weighted ·dt).
      real(wp), allocatable :: vhtr(:, :, :)
         !! Accumulated meridional face transport (m^3, area-weighted ·dt).
      real(wp) :: t_dyn_rel_adv = 0.0_wp
         !! Elapsed dynamics time (s) accumulated since the last tracer
         !! advect / accumulator reset.  Adds `dt` once per outer step.

      ! ---- Phase 2 (6b) windowed-drain workspace ----
      ! Persistent device-resident scratch for the fixed-budget PPM drain
      ! that spends `uhtr/vhtr` at the DT_TRACER_ADVECT boundary.  All
      ! allocated-but-unused at ratio = 1 (the every-step bypass never
      ! touches them); mapped in enter_data so a ratio > 1 run finds them
      ! present.  Sized at the same C-grid face / centre shapes as the
      ! prognostics.  Per-pass re-reconstruction (V2) writes Tr_work /
      ! pal / par / pa6 each sub-cycle pass (Reichl & Hallberg / MOM6
      ! ADVECT_PPM pattern); hprev_work carries the evolving thickness.
      real(wp), allocatable :: hprev_work(:, :, :)
      real(wp), allocatable :: h_win_start(:, :, :)
         !! Layer thickness captured when an accumulation window OPENS.
         !! Paired with `hTr_holds_conc`: the per-stage concentration hold
         !! re-weights `hTr` onto the evolving `h_layer`, and the drain undoes
         !! it with THIS array, which returns `hTr` bit-for-bit to the frozen
         !! window-start content the drain has always consumed.  (Undoing via
         !! the reconstructed `hprev` instead would be equivalent only where
         !! `hprev == h_win_start` exactly, and would silently convert any
         !! reconstruction residual — e.g. from a thin-layer `h_min` clip —
         !! into a tracer-mass drift.)
         !! Evolving (drained) layer thickness during the sub-cycle (m).
      real(wp), allocatable :: uhr_x(:, :, :)
         !! Remaining unspent zonal transport this window (m^3).
      real(wp), allocatable :: uhr_y(:, :, :)
         !! Remaining unspent meridional transport this window (m^3).
      real(wp), allocatable :: uhh_x(:, :, :)
         !! Per-pass limited zonal transport portion (m^3).
      real(wp), allocatable :: uhh_y(:, :, :)
         !! Per-pass limited meridional transport portion (m^3).
      real(wp), allocatable :: tr_flux_x(:, :, :)
         !! Per-pass zonal tracer flux F (m^3 · concentration).
      real(wp), allocatable :: tr_flux_y(:, :, :)
         !! Per-pass meridional tracer flux F (m^3 · concentration).
      real(wp), allocatable :: tr_work(:, :, :)
         !! Current concentration Tr = hTr/hprev_work (rebuilt per pass).
      real(wp), allocatable :: pal(:, :, :)
         !! CW parabola left-edge value per cell (rebuilt per pass).
      real(wp), allocatable :: par(:, :, :)
         !! CW parabola right-edge value per cell (rebuilt per pass).
      real(wp), allocatable :: pa6(:, :, :)
         !! CW parabola curvature a6 = 6·Tr − 3·(aL+aR) (rebuilt per pass).
   contains
      procedure, non_overridable :: init => continuity_init
      procedure, non_overridable :: destroy => continuity_destroy
      procedure, non_overridable :: enter_data => continuity_enter_data
      procedure, non_overridable :: exit_data => continuity_exit_data
      procedure, non_overridable :: bytes => continuity_bytes
   end type continuity_t

   integer, parameter :: RENORM_MAXIT = 8
      !! Newton iterations for the `Sum_k uh_k == uhbt` reconciliation.  The
      !! flux is a NONLINEAR function of the correction `du` because the donor
      !! cell is picked by `sign(u)`: adding `du` can flip a layer's velocity
      !! sign, changing its `h_face`, so a single linear step does not actually
      !! land on `uhbt`.  MOM6 solves the same equation with Newton + bisection
      !! (`zonal_flux_adjust`, up to 20 its).
      !! When no donor flips -- the overwhelmingly common case -- iteration 2
      !! sees zero residual and exits, so the result is BIT-IDENTICAL to the
      !! previous single-step form.
   real(wp), parameter :: RENORM_VR_MIN = 1.0e-12_wp
      !! Floor below which a layer's viscous remnant γ_k is treated as zero
      !! when bracketing `du` — such a layer receives no barotropic increment,
      !! so it constrains nothing and must not divide the CFL bound.
   real(wp), parameter :: RENORM_CFL = 0.25_wp
      !! CFL cap on the reconciliation correction (MOM6 CONTINUITY_CFL_LIMIT,
      !! default 0.5).  MOM6 brackets its Newton solve by this and will ACCEPT
      !! a residual `uhbt` mismatch rather than hand any layer a super-CFL
      !! velocity.  Without the
      !! bracket the solve can satisfy the transport constraint by assigning a
      !! huge `du` to a near-massless layer -- which is exactly the grounded
      !! sliver case here.
   real(wp), parameter :: RENORM_TOL = 1.0e-12_wp
      !! Relative convergence tolerance on `|uhbt - Sum_k uh_k|`.
   integer, parameter :: RENORM_MAXIT_CONSISTENT = 20
      !! Iteration cap when `renorm_consistent_flux` is on (MOM6
      !! `zonal_flux_adjust` also allows 20).  Newton lands on a
      !! single-kink root in <= 3 iterations; the head-room is for the
      !! bisection fallback on a many-layer face with several donor flips.

contains

   subroutine continuity_init(this, grid, nz_ml)
      !! Allocate the 4 face-reconstruction scratch buffers sized at
      !! (nx_face, ny_face, nz).  Default nz=1 covers the barotropic
      !! kernel; passing `nz_ml` sizes them for the multilayer
      !! kernel without forcing a separate init routine.  Ocean
      !! init passes `state%multilayer%nz_ml` when the multilayer
      !! state is in play.
      class(continuity_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      integer, intent(in), optional :: nz_ml
      integer :: nx, ny, nz

      nx = grid%nx_total
      ny = grid%ny_total
      nz = 1
      if (present(nz_ml)) nz = nz_ml

      ! East-face shapes: (nx+1, ny, nz)
      call this%h_face_left_x%init(nx + 1, ny, nz, "continuity_h_face_left_x")
      call this%h_face_right_x%init(nx + 1, ny, nz, "continuity_h_face_right_x")
      ! North-face shapes: (nx, ny+1, nz)
      call this%h_face_left_y%init(nx, ny + 1, nz, "continuity_h_face_left_y")
      call this%h_face_right_y%init(nx, ny + 1, nz, "continuity_h_face_right_y")
      ! Cell-centred conservative min-thickness target scratch: (nx, ny, nz)
      call this%mt_h_new%init(nx, ny, nz, "continuity_mt_h_new")
      call this%mt_grounded%init(nx, ny, 1, "continuity_mt_grounded")
      ! P2 positive-definite outflux-limiter θ scratch (unconditional).
      call this%pd_theta%init(nx, ny, nz, "continuity_pd_theta")

      this%t_dyn_rel_adv = 0.0_wp
      ! Phase 2 flux accumulators + (6b) windowed-drain workspace — 13
      ! 3D arrays (~3.7 GiB at 600²x100) consumed ONLY by the
      ! TR_MODE_ACCUMULATE path (`dt_tracer_advect_ratio > 1`), so their
      ! allocation is gated on `windowed_advection` (latched from cfg
      ! BEFORE init by `ocean_state_init_from_config`, the same
      ! conditional-allocation contract as the default-off closures).
      ! Default .false. — direct `ct%init(...)` call sites that DO drive
      ! the drain opt in explicitly.  Plain host allocation + zero (no `do
      ! concurrent` before enter_data — that would force per-loop H<->D
      ! round-trips; setup code is host-side); mapped onto the device in
      ! continuity_enter_data_impl behind its existing `allocated()`
      ! guards.
      if (this%windowed_advection) then
         allocate (this%uhtr(nx + 1, ny, nz))
         allocate (this%vhtr(nx, ny + 1, nz))
         this%uhtr = 0.0_wp
         this%vhtr = 0.0_wp
         allocate (this%hprev_work(nx, ny, nz), source=0.0_wp)
         allocate (this%h_win_start(nx, ny, nz), source=0.0_wp)
         allocate (this%uhr_x(nx + 1, ny, nz), source=0.0_wp)
         allocate (this%uhr_y(nx, ny + 1, nz), source=0.0_wp)
         allocate (this%uhh_x(nx + 1, ny, nz), source=0.0_wp)
         allocate (this%uhh_y(nx, ny + 1, nz), source=0.0_wp)
         allocate (this%tr_flux_x(nx + 1, ny, nz), source=0.0_wp)
         allocate (this%tr_flux_y(nx, ny + 1, nz), source=0.0_wp)
         allocate (this%tr_work(nx, ny, nz), source=0.0_wp)
         allocate (this%pal(nx, ny, nz), source=0.0_wp)
         allocate (this%par(nx, ny, nz), source=0.0_wp)
         allocate (this%pa6(nx, ny, nz), source=0.0_wp)
      end if

      this%is_init = .true.
   end subroutine continuity_init

   subroutine continuity_destroy(this)
      class(continuity_t), intent(inout) :: this
      this%is_init = .false.
      call this%h_face_left_x%destroy()
      call this%h_face_right_x%destroy()
      call this%h_face_left_y%destroy()
      call this%h_face_right_y%destroy()
      call this%mt_h_new%destroy()
      call this%mt_grounded%destroy()
      call this%pd_theta%destroy()
      if (allocated(this%uhtr)) deallocate (this%uhtr)
      if (allocated(this%vhtr)) deallocate (this%vhtr)
      if (allocated(this%hprev_work)) deallocate (this%hprev_work)
      if (allocated(this%h_win_start)) deallocate (this%h_win_start)
      if (allocated(this%uhr_x)) deallocate (this%uhr_x)
      if (allocated(this%uhr_y)) deallocate (this%uhr_y)
      if (allocated(this%uhh_x)) deallocate (this%uhh_x)
      if (allocated(this%uhh_y)) deallocate (this%uhh_y)
      if (allocated(this%tr_flux_x)) deallocate (this%tr_flux_x)
      if (allocated(this%tr_flux_y)) deallocate (this%tr_flux_y)
      if (allocated(this%tr_work)) deallocate (this%tr_work)
      if (allocated(this%pal)) deallocate (this%pal)
      if (allocated(this%par)) deallocate (this%par)
      if (allocated(this%pa6)) deallocate (this%pa6)
   end subroutine continuity_destroy

   subroutine continuity_enter_data(this)
      !! Bare `copyin(this)` removed (stack-descriptor map → AMD cross-slot
      !! overlap; see ocean_surfstress_enter_data).  The face buffers attach
      !! below; ct-descriptor presence (so DCs touching ct%h_face_left_x%data
      !! don't per-launch memcpy) comes from the root copyin(state) in
      !! ocean_state_enter_data.  A V100 A/B with copyin(this) gone is
      !! bit-identical and faster overall, so the root copy fully covers it.
      class(continuity_t), intent(inout) :: this
      select type (this)
      type is (continuity_t)
         call continuity_enter_data_impl(this)
      end select
   end subroutine continuity_enter_data

   subroutine continuity_enter_data_impl(this)
      type(continuity_t), intent(inout) :: this
      call scratch_3d_buffer_enter_data_impl(this%h_face_left_x)
      call scratch_3d_buffer_enter_data_impl(this%h_face_right_x)
      call scratch_3d_buffer_enter_data_impl(this%h_face_left_y)
      call scratch_3d_buffer_enter_data_impl(this%h_face_right_y)
      call scratch_3d_buffer_enter_data_impl(this%mt_h_new)
      call scratch_3d_buffer_enter_data_impl(this%mt_grounded)
      call scratch_3d_buffer_enter_data_impl(this%pd_theta)
      ! Phase 2 flux accumulators (allocated-but-unused at ratio=1; the
      ! windowed drain in Phase 6b reads/writes them on the device).
      if (allocated(this%uhtr)) then
         !$acc enter data copyin(this%uhtr, this%vhtr)
      end if
      ! Phase 2 (6b) windowed-drain workspace.
      if (allocated(this%hprev_work)) then
         !$acc enter data copyin(this%h_win_start)
         !$acc enter data copyin(this%hprev_work, this%uhr_x, this%uhr_y, &
         !$acc                   this%uhh_x, this%uhh_y, this%tr_flux_x, &
         !$acc                   this%tr_flux_y, this%tr_work, this%pal, &
         !$acc                   this%par, this%pa6)
      end if
   end subroutine continuity_enter_data_impl

   subroutine continuity_exit_data(this)
      class(continuity_t), intent(inout) :: this
      select type (this)
      type is (continuity_t)
         call continuity_exit_data_impl(this)
      end select
   end subroutine continuity_exit_data

   subroutine continuity_exit_data_impl(this)
      type(continuity_t), intent(inout) :: this
      if (allocated(this%hprev_work)) then
         !$acc exit data delete(this%pa6, this%par, this%pal, this%tr_work, &
         !$acc                  this%tr_flux_y, this%tr_flux_x, this%uhh_y, &
         !$acc                  this%uhh_x, this%uhr_y, this%uhr_x, this%hprev_work)
         !$acc exit data delete(this%h_win_start)
      end if
      if (allocated(this%uhtr)) then
         !$acc exit data delete(this%vhtr, this%uhtr)
      end if
      call scratch_3d_buffer_exit_data_impl(this%mt_h_new)
      call scratch_3d_buffer_exit_data_impl(this%mt_grounded)
      call scratch_3d_buffer_exit_data_impl(this%pd_theta)
      call scratch_3d_buffer_exit_data_impl(this%h_face_left_x)
      call scratch_3d_buffer_exit_data_impl(this%h_face_right_x)
      call scratch_3d_buffer_exit_data_impl(this%h_face_left_y)
      call scratch_3d_buffer_exit_data_impl(this%h_face_right_y)
   end subroutine continuity_exit_data_impl

   pure subroutine continuity_compute_fluxes_barotropic(grid, metrics, this, bs)
      !! PPM face reconstruction + per-face mass flux + cell-centred
      !! flux divergence for the barotropic C-grid state.
      !!
      !! Curvilinear (design §2): the per-face TRANSPORT is
      !! `uh = u·h_face·dy_cu` [m³/s] (east) / `vh = v·h_face·dx_cv`
      !! [m³/s] (north), and the divergence is `Δuh·iareaT`.  On
      !! uniform Cartesian `dy_cu = dy`, `dx_cv = dx`, `iareaT =
      !! 1/(dx·dy)`, so `Δ(u·h·dy)/(dx·dy) = Δ(u·h)/dx` bitwise — the
      !! old `inv_dx`/`inv_dy` form.  The PPM reconstruction is a
      !! dimensionless h-difference (no dx enters the slope), so it is
      !! untouched.  `mass_flux_x/y` now carry the m³/s transport, and
      !! every downstream consumer (tracer advect, uhbt renormalise)
      !! reads the same width-weighted flux for tracer consistency.
      !!
      !! Algorithm (per direction, x shown; y mirrors):
      !!
      !!   1. Pass 1: for each cell i with valid 5-point stencil
      !!      (i.e. cells [i-2, i+2] all exist), compute the limited
      !!      slopes δh_{i-1}, δh_i, δh_{i+1} and from them the
      !!      face-left value h_L(i) (= h at i-1/2) and face-right
      !!      value h_R(i) (= h at i+1/2) via CW eq 1.6, then apply
      !!      eq 1.10 monotonic limiter.  Store h_L at
      !!      h_face_right_x(i) (the right state at face i = left
      !!      edge of cell i) and h_R at h_face_left_x(i+1) (left
      !!      state at face i+1 = right edge of cell i).
      !!
      !!   2. Pass 2: for each east face i, pick upwind based on
      !!      u_face_x sign and emit mass_flux_x = u * h_face.
      !!      Domain-wall faces (i=1, nx+1) force mass_flux_x = 0
      !!      (closed-wall BC).
      !!
      !!   3. Pass 3: cell-centred flux divergence flux_h.
      !!
      !! Cells closer than 2 to the boundary use first-order
      !! (h_L = h_R = h_centre) — the stencil is short and the
      !! flux there gets gated by the closed-wall BC anyway.
      !!
      !! Loop order: j-then-i for NVHPC GPU coalescing (CLAUDE.md
      !! memory feedback_do_concurrent_order).
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(continuity_t), intent(inout) :: this
      type(barotropic_state_t), intent(inout) :: bs

      integer :: i, j, nx, ny
      real(wp) :: dh_m1, dh_0, dh_p1, h_left, h_right, u, v, h_face
      real(wp) :: hm2, hm1, h0, hp1, hp2
      logical :: do_pos
      real(wp) :: h_min_pos

      nx = grid%nx_total
      ny = grid%ny_total
      do_pos = this%use_ppm_limit_pos
      h_min_pos = this%h_min

      ! ============================================================
      ! X-DIRECTION reconstruction (5-point stencil per cell)
      ! ============================================================
      ! Interior cells i = 3..nx-2 — full 5-point PPM
      do concurrent(j=1:ny, i=3:nx - 2) &
         local(dh_m1, dh_0, dh_p1, h_left, h_right, &
               hm2, hm1, h0, hp1, hp2)
         ! Mirror-h: replace a LAND neighbour's held floor-h with the
         ! local cell's so the parabola sees a reflected coast (C2).
         h0 = bs%h(i, j)
         hm1 = ppm_mirror_h(bs%h(i - 1, j), h0, metrics%wet_T(i - 1, j))
         hp1 = ppm_mirror_h(bs%h(i + 1, j), h0, metrics%wet_T(i + 1, j))
         hm2 = ppm_mirror_h(bs%h(i - 2, j), hm1, metrics%wet_T(i - 2, j))
         hp2 = ppm_mirror_h(bs%h(i + 2, j), hp1, metrics%wet_T(i + 2, j))
         call ppm_limited_slope(hm2, hm1, h0, dh_m1)
         call ppm_limited_slope(hm1, h0, hp1, dh_0)
         call ppm_limited_slope(h0, hp1, hp2, dh_p1)
         ! Slope-flatten: zero the centre slope if the 3-cell stencil
         ! touches land (as MOM6 does); bit-identical for all-wet (×1).
         dh_0 = dh_0*metrics%wet_T(i - 1, j)*metrics%wet_T(i, j)*metrics%wet_T(i + 1, j)
         h_left = 0.5_wp*(hm1 + h0) - (dh_0 - dh_m1)/6.0_wp
         h_right = 0.5_wp*(h0 + hp1) - (dh_p1 - dh_0)/6.0_wp
         call ppm_cell_limiter(h0, h_left, h_right)
         if (do_pos) call ppm_limit_pos(h0, h_left, h_right, h_min_pos)
         this%h_face_right_x%data(i, j, 1) = h_left
         this%h_face_left_x%data(i + 1, j, 1) = h_right
      end do
      ! Boundary cells (i = 1, 2, nx-1, nx): 1st-order — h_face_*
      ! at the four border faces (1, 2, nx, nx+1) just take the
      ! abutting cell's centre value.  Face 1 and face nx+1 are
      ! walls (mass_flux forced to 0); face 2 and face nx use a
      ! one-sided downwind reconstruction equivalent to 1st-order.
      do concurrent(j=1:ny)
         this%h_face_left_x%data(1, j, 1) = bs%h(1, j)
         this%h_face_right_x%data(1, j, 1) = bs%h(1, j)
         this%h_face_left_x%data(2, j, 1) = bs%h(1, j)
         this%h_face_right_x%data(2, j, 1) = bs%h(2, j)
         ! Face 3: cell 2's right edge (h_face_left at face 3) falls
         ! back to 1st order — its 5-point PPM stencil needs cell 0.
         ! Same for cell 2's left-edge contribution at face 3
         ! (h_face_right at face 3).
         this%h_face_left_x%data(3, j, 1) = bs%h(2, j)
         this%h_face_right_x%data(3, j, 1) = bs%h(2, j)
         ! Face nx-1: mirror of face 3.
         this%h_face_right_x%data(nx - 1, j, 1) = bs%h(nx - 1, j)
         this%h_face_left_x%data(nx, j, 1) = bs%h(nx - 1, j)
         this%h_face_right_x%data(nx, j, 1) = bs%h(nx, j)
         this%h_face_left_x%data(nx + 1, j, 1) = bs%h(nx, j)
         this%h_face_right_x%data(nx + 1, j, 1) = bs%h(nx, j)
      end do

      ! X-direction face transport (upwind pick from h_face_left/right_x).
      ! uh = u·h_face·dy_cu(i,j) [m³/s] — face-width-weighted (design §2).
      do concurrent(j=1:ny, i=2:nx) local(u, h_face)
         u = bs%u_face_x(i, j)
         if (u >= 0.0_wp) then
            h_face = this%h_face_left_x%data(i, j, 1)
         else
            h_face = this%h_face_right_x%data(i, j, 1)
         end if
         bs%mass_flux_x(i, j) = u*h_face*metrics%dy_cu(i, j)
      end do
      do concurrent(j=1:ny)
         bs%mass_flux_x(1, j) = 0.0_wp
         bs%mass_flux_x(nx + 1, j) = 0.0_wp
      end do

      ! ============================================================
      ! Y-DIRECTION reconstruction (mirror of X)
      ! ============================================================
      do concurrent(j=3:ny - 2, i=1:nx) &
         local(dh_m1, dh_0, dh_p1, h_left, h_right, &
               hm2, hm1, h0, hp1, hp2)
         h0 = bs%h(i, j)
         hm1 = ppm_mirror_h(bs%h(i, j - 1), h0, metrics%wet_T(i, j - 1))
         hp1 = ppm_mirror_h(bs%h(i, j + 1), h0, metrics%wet_T(i, j + 1))
         hm2 = ppm_mirror_h(bs%h(i, j - 2), hm1, metrics%wet_T(i, j - 2))
         hp2 = ppm_mirror_h(bs%h(i, j + 2), hp1, metrics%wet_T(i, j + 2))
         call ppm_limited_slope(hm2, hm1, h0, dh_m1)
         call ppm_limited_slope(hm1, h0, hp1, dh_0)
         call ppm_limited_slope(h0, hp1, hp2, dh_p1)
         dh_0 = dh_0*metrics%wet_T(i, j - 1)*metrics%wet_T(i, j)*metrics%wet_T(i, j + 1)
         h_left = 0.5_wp*(hm1 + h0) - (dh_0 - dh_m1)/6.0_wp
         h_right = 0.5_wp*(h0 + hp1) - (dh_p1 - dh_0)/6.0_wp
         call ppm_cell_limiter(h0, h_left, h_right)
         if (do_pos) call ppm_limit_pos(h0, h_left, h_right, h_min_pos)
         this%h_face_right_y%data(i, j, 1) = h_left
         this%h_face_left_y%data(i, j + 1, 1) = h_right
      end do
      do concurrent(i=1:nx)
         this%h_face_left_y%data(i, 1, 1) = bs%h(i, 1)
         this%h_face_right_y%data(i, 1, 1) = bs%h(i, 1)
         this%h_face_left_y%data(i, 2, 1) = bs%h(i, 1)
         this%h_face_right_y%data(i, 2, 1) = bs%h(i, 2)
         this%h_face_left_y%data(i, 3, 1) = bs%h(i, 2)
         this%h_face_right_y%data(i, 3, 1) = bs%h(i, 2)
         this%h_face_right_y%data(i, ny - 1, 1) = bs%h(i, ny - 1)
         this%h_face_left_y%data(i, ny, 1) = bs%h(i, ny - 1)
         this%h_face_right_y%data(i, ny, 1) = bs%h(i, ny)
         this%h_face_left_y%data(i, ny + 1, 1) = bs%h(i, ny)
         this%h_face_right_y%data(i, ny + 1, 1) = bs%h(i, ny)
      end do

      do concurrent(j=2:ny, i=1:nx) local(v, h_face)
         v = bs%v_face_y(i, j)
         if (v >= 0.0_wp) then
            h_face = this%h_face_left_y%data(i, j, 1)
         else
            h_face = this%h_face_right_y%data(i, j, 1)
         end if
         bs%mass_flux_y(i, j) = v*h_face*metrics%dx_cv(i, j)
      end do
      do concurrent(i=1:nx)
         bs%mass_flux_y(i, 1) = 0.0_wp
         bs%mass_flux_y(i, ny + 1) = 0.0_wp
      end do

      ! ============================================================
      ! Cell-centred flux divergence (transport divergence · iareaT)
      ! ============================================================
      do concurrent(j=1:ny, i=1:nx)
         bs%flux_h(i, j) = &
            ((bs%mass_flux_x(i + 1, j) - bs%mass_flux_x(i, j)) + &
             (bs%mass_flux_y(i, j + 1) - bs%mass_flux_y(i, j)))*metrics%iareaT(i, j)
      end do
   end subroutine continuity_compute_fluxes_barotropic

   pure subroutine continuity_apply_fluxes_barotropic(bs, dt)
      !! Forward-Euler step: h <- h - dt * flux_h.  The full
      !! split-explicit RK2 scheme (Phase 4) wraps two of these calls
      !! around an RK2 averaging pass; for Phase 2 this single-stage
      !! step is enough to exercise the kernel under the lake-at-rest,
      !! Gaussian-hump, and mass-conservation tests.
      type(barotropic_state_t), intent(inout) :: bs
      real(wp), intent(in) :: dt
      integer :: i, j, nx, ny

      nx = size(bs%h, 1)
      ny = size(bs%h, 2)

      do concurrent(j=1:ny, i=1:nx)
         bs%h(i, j) = bs%h(i, j) - dt*bs%flux_h(i, j)
      end do
   end subroutine continuity_apply_fluxes_barotropic

   pure subroutine continuity_compute_fluxes(grid, metrics, this, ms)
      !! **Test-only** (no production caller): the unsplit reference path,
      !! kept as the oracle the split production path is checked against.
      !! Multilayer counterpart to
      !! `continuity_compute_fluxes_barotropic`: identical PPM
      !! reconstruction + upwind face pick + flux divergence, lifted
      !! per-layer.  Each k-slice is independent (the PPM stencil
      !! reads only the same k), so the do-concurrent kernels
      !! parallelize over (k, j, i) simultaneously for GPU
      !! occupancy.
      !!
      !! Workspaces (`this%h_face_*_x/y`) must have been initialised
      !! with `nz_ml` matching `ms%nz_ml` — handled by passing the
      !! optional `nz_ml` to `continuity_init`.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(continuity_t), intent(inout) :: this
      type(multilayer_state_t), intent(inout) :: ms

      integer :: i, j, k, nx, ny, nz
      real(wp) :: dh_m1, dh_0, dh_p1, h_left, h_right, u, v, h_face
      real(wp) :: hm2, hm1, h0, hp1, hp2
      logical :: do_pos
      real(wp) :: h_min_pos

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml
      do_pos = this%use_ppm_limit_pos
      h_min_pos = this%h_min

      ! ============================================================
      ! X-DIRECTION reconstruction (5-point stencil per cell)
      ! ============================================================
      do concurrent(k=1:nz, j=1:ny, i=3:nx - 2) &
         local(dh_m1, dh_0, dh_p1, h_left, h_right, hm2, hm1, h0, hp1, hp2)
         ! Mirror-h at land neighbours (C2); bit-identical for all-wet.
         h0 = ms%h_layer(i, j, k)
         hm1 = ppm_mirror_h(ms%h_layer(i - 1, j, k), h0, metrics%wet_T(i - 1, j))
         hp1 = ppm_mirror_h(ms%h_layer(i + 1, j, k), h0, metrics%wet_T(i + 1, j))
         hm2 = ppm_mirror_h(ms%h_layer(i - 2, j, k), hm1, metrics%wet_T(i - 2, j))
         hp2 = ppm_mirror_h(ms%h_layer(i + 2, j, k), hp1, metrics%wet_T(i + 2, j))
         call ppm_limited_slope(hm2, hm1, h0, dh_m1)
         call ppm_limited_slope(hm1, h0, hp1, dh_0)
         call ppm_limited_slope(h0, hp1, hp2, dh_p1)
         dh_0 = dh_0*metrics%wet_T(i - 1, j)*metrics%wet_T(i, j)*metrics%wet_T(i + 1, j)
         h_left = 0.5_wp*(hm1 + h0) - (dh_0 - dh_m1)/6.0_wp
         h_right = 0.5_wp*(h0 + hp1) - (dh_p1 - dh_0)/6.0_wp
         call ppm_cell_limiter(h0, h_left, h_right)
         if (do_pos) call ppm_limit_pos(h0, h_left, h_right, h_min_pos)
         this%h_face_right_x%data(i, j, k) = h_left
         this%h_face_left_x%data(i + 1, j, k) = h_right
      end do
      ! Boundary cells: 1st-order fallback
      do concurrent(k=1:nz, j=1:ny)
         this%h_face_left_x%data(1, j, k) = ms%h_layer(1, j, k)
         this%h_face_right_x%data(1, j, k) = ms%h_layer(1, j, k)
         this%h_face_left_x%data(2, j, k) = ms%h_layer(1, j, k)
         this%h_face_right_x%data(2, j, k) = ms%h_layer(2, j, k)
         ! Face 3: cell 2's right edge falls back to 1st order (its
         ! 5-point stencil needs cell 0); cell 3's left edge is set
         ! by the interior loop above.
         this%h_face_left_x%data(3, j, k) = ms%h_layer(2, j, k)
         this%h_face_right_x%data(3, j, k) = ms%h_layer(2, j, k)
         ! Face nx-1: mirror of face 3.  Cell nx-1's left edge falls
         ! back to 1st order; cell nx-2's right edge came from the
         ! interior loop.
         this%h_face_right_x%data(nx - 1, j, k) = ms%h_layer(nx - 1, j, k)
         this%h_face_left_x%data(nx, j, k) = ms%h_layer(nx - 1, j, k)
         this%h_face_right_x%data(nx, j, k) = ms%h_layer(nx, j, k)
         this%h_face_left_x%data(nx + 1, j, k) = ms%h_layer(nx, j, k)
         this%h_face_right_x%data(nx + 1, j, k) = ms%h_layer(nx, j, k)
      end do

      ! X-direction face transport (upwind pick) — width-weighted dy_cu
      do concurrent(k=1:nz, j=1:ny, i=2:nx) local(u, h_face)
         u = ms%u_face_x_layer(i, j, k)
         if (u >= 0.0_wp) then
            h_face = this%h_face_left_x%data(i, j, k)
         else
            h_face = this%h_face_right_x%data(i, j, k)
         end if
         ms%mass_flux_x_layer(i, j, k) = u*h_face*metrics%dy_cu(i, j)
      end do
      do concurrent(k=1:nz, j=1:ny)
         ms%mass_flux_x_layer(1, j, k) = 0.0_wp
         ms%mass_flux_x_layer(nx + 1, j, k) = 0.0_wp
      end do
      ! Physical wall zeroing — see `continuity_zonal_flux` for the
      ! full rationale.  Must mirror the split form's wall closure or
      ! `test_split_zonal_only_matches_unsplit` breaks.
      do concurrent(k=1:nz, j=1:ny)
         ms%mass_flux_x_layer(grid%nghost + 1, j, k) = 0.0_wp
         ms%mass_flux_x_layer(grid%nghost + grid%nx_phys + 1, j, k) = 0.0_wp
      end do

      ! ---- Porous barriers (Adcroft 2013) ----
      ! Narrow the layer transport by the OPEN-AREA fraction of the face.
      ! Host-side gate: with the knob off there is no kernel launch and
      ! the loops above are textually unchanged, so the whole path is
      ! byte-identical.  Applied AFTER the wall zeroing (0 stays 0) and
      ! BEFORE the barotropic renormalisation, which must see the narrowed
      ! transports it is constraining.
      !
      ! WRITTEN INLINE, not as a call to `porous_narrow_3d`.  Handing
      ! `ms%mass_flux_*_layer` to an external subroutine as an
      ! `intent(inout)` actual makes nvfortran treat the array as ESCAPING,
      ! which pessimises every `do concurrent` in this routine — even
      ! though the branch never runs with the knob off.  Measured on a
      ! 600x600x50 default-path (porous OFF) double-gyre, V100: the call
      ! form costs `ocean_continuity` 10.58 s vs 9.11 s inline (+4.8% on
      ! total solver time vs origin/main; the inline form is +0.5%, i.e.
      ! noise).  `rdb_coriolis_adv` keeps the shared `porous_narrow_3d`
      ! helper — measured there at no cost.
      if (metrics%use_porous) then
         do concurrent(k=1:nz, j=1:ny, i=1:nx + 1)
            ms%mass_flux_x_layer(i, j, k) = ms%mass_flux_x_layer(i, j, k)* &
                                            metrics%por_face_area_u(i, j, k)
         end do
      end if

      ! ---- z-level closed faces: see the composition rule on
      ! `ocean_metrics_t%open_v`.  A SEPARATE pass, not composed into
      ! `por_face_area_u`: the two gates are independent and the porous
      ! fraction is refreshed per outer step while this mask is static.
      ! Inline for the same escaping-array reason as the porous pass.
      if (metrics%use_closed_faces) then
         do concurrent(k=1:nz, j=1:ny, i=1:nx + 1)
            ms%mass_flux_x_layer(i, j, k) = ms%mass_flux_x_layer(i, j, k)* &
                                            metrics%open_u(i, j, k)
         end do
      end if

      ! ============================================================
      ! Y-DIRECTION reconstruction
      ! ============================================================
      do concurrent(k=1:nz, j=3:ny - 2, i=1:nx) &
         local(dh_m1, dh_0, dh_p1, h_left, h_right, hm2, hm1, h0, hp1, hp2)
         h0 = ms%h_layer(i, j, k)
         hm1 = ppm_mirror_h(ms%h_layer(i, j - 1, k), h0, metrics%wet_T(i, j - 1))
         hp1 = ppm_mirror_h(ms%h_layer(i, j + 1, k), h0, metrics%wet_T(i, j + 1))
         hm2 = ppm_mirror_h(ms%h_layer(i, j - 2, k), hm1, metrics%wet_T(i, j - 2))
         hp2 = ppm_mirror_h(ms%h_layer(i, j + 2, k), hp1, metrics%wet_T(i, j + 2))
         call ppm_limited_slope(hm2, hm1, h0, dh_m1)
         call ppm_limited_slope(hm1, h0, hp1, dh_0)
         call ppm_limited_slope(h0, hp1, hp2, dh_p1)
         dh_0 = dh_0*metrics%wet_T(i, j - 1)*metrics%wet_T(i, j)*metrics%wet_T(i, j + 1)
         h_left = 0.5_wp*(hm1 + h0) - (dh_0 - dh_m1)/6.0_wp
         h_right = 0.5_wp*(h0 + hp1) - (dh_p1 - dh_0)/6.0_wp
         call ppm_cell_limiter(h0, h_left, h_right)
         if (do_pos) call ppm_limit_pos(h0, h_left, h_right, h_min_pos)
         this%h_face_right_y%data(i, j, k) = h_left
         this%h_face_left_y%data(i, j + 1, k) = h_right
      end do
      do concurrent(k=1:nz, i=1:nx)
         this%h_face_left_y%data(i, 1, k) = ms%h_layer(i, 1, k)
         this%h_face_right_y%data(i, 1, k) = ms%h_layer(i, 1, k)
         this%h_face_left_y%data(i, 2, k) = ms%h_layer(i, 1, k)
         this%h_face_right_y%data(i, 2, k) = ms%h_layer(i, 2, k)
         this%h_face_left_y%data(i, 3, k) = ms%h_layer(i, 2, k)
         this%h_face_right_y%data(i, 3, k) = ms%h_layer(i, 2, k)
         this%h_face_right_y%data(i, ny - 1, k) = ms%h_layer(i, ny - 1, k)
         this%h_face_left_y%data(i, ny, k) = ms%h_layer(i, ny - 1, k)
         this%h_face_right_y%data(i, ny, k) = ms%h_layer(i, ny, k)
         this%h_face_left_y%data(i, ny + 1, k) = ms%h_layer(i, ny, k)
         this%h_face_right_y%data(i, ny + 1, k) = ms%h_layer(i, ny, k)
      end do

      do concurrent(k=1:nz, j=2:ny, i=1:nx) local(v, h_face)
         v = ms%v_face_y_layer(i, j, k)
         if (v >= 0.0_wp) then
            h_face = this%h_face_left_y%data(i, j, k)
         else
            h_face = this%h_face_right_y%data(i, j, k)
         end if
         ms%mass_flux_y_layer(i, j, k) = v*h_face*metrics%dx_cv(i, j)
      end do
      do concurrent(k=1:nz, i=1:nx)
         ms%mass_flux_y_layer(i, 1, k) = 0.0_wp
         ms%mass_flux_y_layer(i, ny + 1, k) = 0.0_wp
      end do
      ! Physical wall zeroing — mirrors split form's
      ! `continuity_meridional_flux`.
      do concurrent(k=1:nz, i=1:nx)
         ms%mass_flux_y_layer(i, grid%nghost + 1, k) = 0.0_wp
         ms%mass_flux_y_layer(i, grid%nghost + grid%ny_phys + 1, k) = 0.0_wp
      end do

      ! ---- Porous barriers: see the zonal twin (incl. why it is inline) ----
      if (metrics%use_porous) then
         do concurrent(k=1:nz, j=1:ny + 1, i=1:nx)
            ms%mass_flux_y_layer(i, j, k) = ms%mass_flux_y_layer(i, j, k)* &
                                            metrics%por_face_area_v(i, j, k)
         end do
      end if

      ! ---- z-level closed faces: see the zonal twin ----
      if (metrics%use_closed_faces) then
         do concurrent(k=1:nz, j=1:ny + 1, i=1:nx)
            ms%mass_flux_y_layer(i, j, k) = ms%mass_flux_y_layer(i, j, k)* &
                                            metrics%open_v(i, j, k)
         end do
      end if

      ! ============================================================
      ! Per-layer flux divergence (transport divergence · iareaT)
      ! ============================================================
      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         ms%flux_h_layer(i, j, k) = &
            ((ms%mass_flux_x_layer(i + 1, j, k) - &
              ms%mass_flux_x_layer(i, j, k)) + &
             (ms%mass_flux_y_layer(i, j + 1, k) - &
              ms%mass_flux_y_layer(i, j, k)))*metrics%iareaT(i, j)
      end do
   end subroutine continuity_compute_fluxes

   pure subroutine continuity_apply_fluxes(ms, dt, h_min)
      !! **Test-only** (no production caller): unsplit apply, paired with
      !! `continuity_compute_fluxes` as the split path's reference oracle.
      !! Per-layer forward-Euler thickness update.
      !!
      !! Optional `h_min` (m): when > 0 and the active vcoord is
      !! VCOORD_LAGRANGIAN, applies max(h_new, h_min) on the h-update so
      !! grounding layers cannot go below the floor.  0.0 (default absent) ⇒
      !! original update verbatim (bit-identical).  The `mass_budget_continuity`
      !! accumulator records the RAW divergence tendency regardless — the floor
      !! injection shows up in the mass Error diagnostic (R2).
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: dt
      real(wp), intent(in), optional :: h_min
         !! Minimum-thickness floor (m). 0 or absent ⇒ off ⇒ bit-identical.
      integer :: i, j, k, nx, ny, nz
      real(wp) :: h_min_use

      h_min_use = 0.0_wp
      if (present(h_min)) h_min_use = h_min

      nx = size(ms%h_layer, 1)
      ny = size(ms%h_layer, 2)
      nz = size(ms%h_layer, 3)

      if (h_min_use > 0.0_wp) then
         do concurrent(k=1:nz, j=1:ny, i=1:nx)
            ms%h_layer(i, j, k) = max(ms%h_layer(i, j, k) - dt*ms%flux_h_layer(i, j, k), h_min_use)
            ms%mass_budget_continuity(i, j, k) = &
               ms%mass_budget_continuity(i, j, k) - dt*ms%flux_h_layer(i, j, k)
         end do
      else
         do concurrent(k=1:nz, j=1:ny, i=1:nx)
            ms%h_layer(i, j, k) = ms%h_layer(i, j, k) - dt*ms%flux_h_layer(i, j, k)
            ms%mass_budget_continuity(i, j, k) = &
               ms%mass_budget_continuity(i, j, k) - dt*ms%flux_h_layer(i, j, k)
         end do
      end if
   end subroutine continuity_apply_fluxes

   pure subroutine continuity_zonal_flux(grid, metrics, this, ms, dt, uhbt, bc, &
                                         visc_rem, u_cor)
      !! Zonal (x-only) PPM reconstruction + per-face mass flux on
      !! the multilayer C-grid.  Companion to
      !! `continuity_meridional_flux` for the
      !! directionally-split (Lie) continuity step.  Writes
      !! `ms%mass_flux_x_layer` and leaves `mass_flux_y_layer` /
      !! `flux_h_layer` untouched.  Wall faces at i=1 and i=nx+1
      !! are zeroed (closed-wall BC).
      !!
      !! Stencil + boundary treatment identical to the X half of
      !! `continuity_compute_fluxes` — the body was
      !! lifted verbatim and the Y block dropped.
      !!
      !! Optional `uhbt(i, j)` is the time-mean barotropic-substep transport
      !! at the east face.  When present, the per-layer mass fluxes
      !! are renormalised by a uniform velocity correction
      !! `du(i, j) = (uhbt - Σ_k mass_flux) / Σ_k h_face` so the
      !! vertical sum matches `uhbt` exactly.  This is the MOM6
      !! continuity-with-uhbt pattern.
      !! With this constraint, the slow continuity advances h_layer
      !! consistently with the barotropic-substep's η end-state — no
      !! post-rescale needed.  Absent ⇒ unconstrained PPM transport.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(continuity_t), intent(inout) :: this
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: dt
         !! Time increment (s).  Used only for the swept-volume CFL when
         !! `this%vol_cfl = .true.`; ignored (bit-identical) otherwise.
      real(wp), intent(in), optional :: uhbt(:, :)
         !! Time-mean east-face transport from barotropic substep (m³/s),
         !! shape `(nx+1, ny)`.  Width-weighted (carries `dy_cu`) so it
         !! constrains the same transport `mass_flux_x_layer` now holds.
      type(ocean_bc_state_t), intent(in), optional :: bc
         !! Per-edge OBC tags.  Default (absent) -> OBC_WALL on both
         !! ends.  Non-WALL tags skip the wall-zero step at that edge.
      ! assumed-shape-ok: pure passthrough to the renormaliser, which carries
      ! the same waiver; never indexed here.
      real(wp), intent(in), optional :: visc_rem(:, :, :)
         !! Per-layer viscous remnant γ_k, forwarded to
         !! `renormalise_zonal_flux_to_uhbt`.  Absent ⇒ γ ≡ 1, bit-identical.
      ! assumed-shape-ok: pure passthrough to the renormaliser.
      real(wp), intent(inout), optional :: u_cor(:, :, :)
         !! Forwarded MOM6 `u_cor` destination — a SEPARATE time-mean field,
         !! never the prognostic.  Absent ⇒ flux-only, bit-identical.

      integer :: i, j, k, nx, ny, nz
      integer :: bc_w_tag, bc_e_tag
      real(wp) :: dh_m1, dh_0, dh_p1, h_left, h_right, u, h_face
      real(wp) :: hm2, hm1, h0, hp1, hp2
      real(wp) :: cfl_d, curv3_d, dh_d
      logical :: do_pos, do_volcfl, do_pd
      logical :: has_w_flux, has_e_flux, renorm_skip_walls
      real(wp) :: h_min_pos, two_h_lim

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml
      do_pos = this%use_ppm_limit_pos
      do_volcfl = this%vol_cfl
      h_min_pos = this%h_min
      ! P1 positive-definite reconstruction floor (MOM6's positive-definite
      ! PPM).  Scalars copied to locals so the DC loop never walks the
      ! continuity_t descriptor.  Off ⇒ untaken branch ⇒ bit-identical.
      do_pd = this%positive_definite
      two_h_lim = 2.0_wp*this%h_lim

      do concurrent(k=1:nz, j=1:ny, i=3:nx - 2) &
         local(dh_m1, dh_0, dh_p1, h_left, h_right, hm2, hm1, h0, hp1, hp2)
         h0 = ms%h_layer(i, j, k)
         hm1 = ppm_mirror_h(ms%h_layer(i - 1, j, k), h0, metrics%wet_T(i - 1, j))
         hp1 = ppm_mirror_h(ms%h_layer(i + 1, j, k), h0, metrics%wet_T(i + 1, j))
         hm2 = ppm_mirror_h(ms%h_layer(i - 2, j, k), hm1, metrics%wet_T(i - 2, j))
         hp2 = ppm_mirror_h(ms%h_layer(i + 2, j, k), hp1, metrics%wet_T(i + 2, j))
         call ppm_limited_slope(hm2, hm1, h0, dh_m1)
         call ppm_limited_slope(hm1, h0, hp1, dh_0)
         call ppm_limited_slope(h0, hp1, hp2, dh_p1)
         dh_0 = dh_0*metrics%wet_T(i - 1, j)*metrics%wet_T(i, j)*metrics%wet_T(i + 1, j)
         h_left = 0.5_wp*(hm1 + h0) - (dh_0 - dh_m1)/6.0_wp
         h_right = 0.5_wp*(h0 + hp1) - (dh_p1 - dh_0)/6.0_wp
         call ppm_cell_limiter(h0, h_left, h_right)
         if (do_pos) call ppm_limit_pos(h0, h_left, h_right, h_min_pos)
         if (do_pd) then
            h_left = max(h_left, two_h_lim)
            h_right = max(h_right, two_h_lim)
         end if
         this%h_face_right_x%data(i, j, k) = h_left
         this%h_face_left_x%data(i + 1, j, k) = h_right
      end do
      do concurrent(k=1:nz, j=1:ny)
         this%h_face_left_x%data(1, j, k) = ms%h_layer(1, j, k)
         this%h_face_right_x%data(1, j, k) = ms%h_layer(1, j, k)
         this%h_face_left_x%data(2, j, k) = ms%h_layer(1, j, k)
         this%h_face_right_x%data(2, j, k) = ms%h_layer(2, j, k)
         ! Face 3: cell 2's right edge falls back to 1st order (its
         ! 5-point stencil needs cell 0); cell 3's left edge is set
         ! by the interior loop above.
         this%h_face_left_x%data(3, j, k) = ms%h_layer(2, j, k)
         this%h_face_right_x%data(3, j, k) = ms%h_layer(2, j, k)
         ! Face nx-1: mirror of face 3.  Cell nx-1's left edge falls
         ! back to 1st order; cell nx-2's right edge came from the
         ! interior loop.
         this%h_face_right_x%data(nx - 1, j, k) = ms%h_layer(nx - 1, j, k)
         this%h_face_left_x%data(nx, j, k) = ms%h_layer(nx - 1, j, k)
         this%h_face_right_x%data(nx, j, k) = ms%h_layer(nx, j, k)
         this%h_face_left_x%data(nx + 1, j, k) = ms%h_layer(nx, j, k)
         this%h_face_right_x%data(nx + 1, j, k) = ms%h_layer(nx, j, k)
      end do

      do concurrent(k=1:nz, j=1:ny, i=2:nx) &
         local(u, h_face, cfl_d, curv3_d, dh_d)
         u = ms%u_face_x_layer(i, j, k)
         if (u >= 0.0_wp) then
            ! Donor = cell i-1; downwind (east) edge = h_face_left_x(i).
            h_face = this%h_face_left_x%data(i, j, k)
            if (do_volcfl) then
               ! Swept-oriented donor edge diff dh = h_L - h_R (west-east).
               dh_d = this%h_face_right_x%data(i - 1, j, k) - h_face
               curv3_d = this%h_face_right_x%data(i - 1, j, k) + h_face &
                         - 2.0_wp*ms%h_layer(i - 1, j, k)
               cfl_d = u*dt*metrics%dy_cu(i, j)*metrics%iareaT(i - 1, j)
               h_face = volcfl_face(h_face, dh_d, curv3_d, cfl_d)
            end if
         else
            ! Donor = cell i; downwind (west) edge = h_face_right_x(i).
            h_face = this%h_face_right_x%data(i, j, k)
            if (do_volcfl) then
               ! Swept-oriented donor edge diff dh = h_R - h_L (east-west).
               dh_d = this%h_face_left_x%data(i + 1, j, k) - h_face
               curv3_d = h_face + this%h_face_left_x%data(i + 1, j, k) &
                         - 2.0_wp*ms%h_layer(i, j, k)
               cfl_d = (-u)*dt*metrics%dy_cu(i, j)*metrics%iareaT(i, j)
               h_face = volcfl_face(h_face, dh_d, curv3_d, cfl_d)
            end if
         end if
         ms%mass_flux_x_layer(i, j, k) = u*h_face*metrics%dy_cu(i, j)
      end do
      do concurrent(k=1:nz, j=1:ny)
         ms%mass_flux_x_layer(1, j, k) = 0.0_wp
         ms%mass_flux_x_layer(nx + 1, j, k) = 0.0_wp
      end do
      ! Zero the PHYSICAL wall faces (not just the array edges).  The
      ! ocean slow-path apply kernels (surface stress / Coriolis / PGF
      ! / drag) don't enforce u_face_x = 0 at the physical interior
      ! boundary, so even with `u^n_wall = 0` the velocity applies
      ! between stages can drive a small `u^{stage1}_wall` that the
      ! next stage's `mass_flux = u·h_face` would treat as a wall
      ! flux.  Without this zeroing, the PPM scheme transports
      ! tracer mass across the physical wall into ghost cells, and
      ! `sum(hTr)` drifts by ~1e-6 over hundreds of outer steps.
      ! Mirrors MOM6's `mask2dCu` face mask convention.  Under MPI
      ! decomposition the gate now reads bc%has_west / bc%has_east
      ! (the flag mirrors decomp%has_west / has_east set at driver
      ! setup) so a subdomain seam is never hard-zeroed here — the
      ! halo exchange corrects it.
      !
      ! OBC dispatch: non-WALL edges keep the computed mass flux so
      ! the downstream transport sees the barotropic-substep's open-boundary
      ! velocity.  Caller must supply matching `uhbt` at those faces
      ! (the driver no longer zeros it for non-WALL edges).
      bc_w_tag = OBC_WALL
      bc_e_tag = OBC_WALL
      ! Physical-edge flags: default .true. => single-rank bit-identical
      ! behaviour; .false. at an MPI seam skips the hard zero (O0, D0).
      has_w_flux = .true.
      has_e_flux = .true.
      if (present(bc)) then
         bc_w_tag = bc%west%bc_type
         bc_e_tag = bc%east%bc_type
         has_w_flux = bc%has_west
         has_e_flux = bc%has_east
      end if
      do concurrent(k=1:nz, j=1:ny)
         if (bc_w_tag == OBC_WALL .and. has_w_flux) ms%mass_flux_x_layer(grid%nghost + 1, j, k) = 0.0_wp
         if (bc_e_tag == OBC_WALL .and. has_e_flux) ms%mass_flux_x_layer(grid%nghost + grid%nx_phys + 1, j, k) = 0.0_wp
      end do

      ! ---- Porous barriers (Adcroft 2013) ----
      ! Narrow the layer transport by the OPEN-AREA fraction of the face.
      ! Host-side gate: with the knob off there is no kernel launch and
      ! the loops above are textually unchanged, so the whole path is
      ! byte-identical.  Applied AFTER the wall zeroing (0 stays 0) and
      ! BEFORE the barotropic renormalisation, which must see the narrowed
      ! transports it is constraining.
      !
      ! WRITTEN INLINE, not as a call to `porous_narrow_3d`.  Handing
      ! `ms%mass_flux_*_layer` to an external subroutine as an
      ! `intent(inout)` actual makes nvfortran treat the array as ESCAPING,
      ! which pessimises every `do concurrent` in this routine — even
      ! though the branch never runs with the knob off.  Measured on a
      ! 600x600x50 default-path (porous OFF) double-gyre, V100: the call
      ! form costs `ocean_continuity` 10.58 s vs 9.11 s inline (+4.8% on
      ! total solver time vs origin/main; the inline form is +0.5%, i.e.
      ! noise).  `rdb_coriolis_adv` keeps the shared `porous_narrow_3d`
      ! helper — measured there at no cost.
      if (metrics%use_porous) then
         do concurrent(k=1:nz, j=1:ny, i=1:nx + 1)
            ms%mass_flux_x_layer(i, j, k) = ms%mass_flux_x_layer(i, j, k)* &
                                            metrics%por_face_area_u(i, j, k)
         end do
      end if

      ! ---- z-level closed faces: see the composition rule on
      ! `ocean_metrics_t%open_v`.  Also before the renormalisation — the
      ! renormaliser must distribute `uhbt` over the OPEN layers ONLY, so
      ! it takes the SAME mask as a weight (`use_open`/`open` below).
      if (metrics%use_closed_faces) then
         do concurrent(k=1:nz, j=1:ny, i=1:nx + 1)
            ms%mass_flux_x_layer(i, j, k) = ms%mass_flux_x_layer(i, j, k)* &
                                            metrics%open_u(i, j, k)
         end do
      end if

      ! ---- MOM6-style transport constraint ----
      ! Renormalise the per-layer mass flux so its vertical sum
      ! matches the barotropic-substep's time-mean transport.  Single Newton
      ! pass: assuming the upwind donor doesn't flip under the
      ! correction (true for small `du`), `Σ_k (u_k + du) · h_face_k =
      ! uhbt` solves to `du = (uhbt - Σ_k u_k·h_face_k) / Σ_k h_face_k`.
      ! Apply `du` per layer; the sign of u_face decides the donor
      ! (left vs right h_face), unchanged from above.
      !
      ! Wall faces (i=1, i=nx+1, plus the physical-wall zeroing
      ! above) keep mass_flux at zero and are skipped via uhbt(I,j)=0
      ! at those indices (caller's barotropic substep initialises and only
      ! accumulates on owned interior faces).
      if (present(uhbt)) then
         ! Any non-WALL zonal edge (periodic OR open-class Flather/tidal/
         ! Chapman/clamped) carries genuine transport at the physical-wall
         ! face: its mass flux was NOT zeroed above, and the barotropic
         ! substep supplies a nonzero uhbt there.  Those faces must be
         ! renormalised like interior faces so the per-layer fluxes sum to
         ! the BT transport — otherwise the boundary leaks mode mismatch.
         ! For a genuine closed WALL (uhbt = 0, mass_flux = 0) the wall face
         ! is skipped (default skip_walls=.true.): nothing to renormalise.
         ! All-WALL → skip_walls stays .true. → bit-identical to pre-OBC.
         ! `visc_rem` / `u_cor` forward straight through: an ABSENT
         ! optional passed to an optional dummy stays absent, so the no-knob
         ! path reaches the renormaliser exactly as before (bit-identical).
         ! ONE call site per `por` actual, not one per BC branch.  `skip_w`
         ! and the `has_*` flags are now computed here and always passed:
         ! when `skip_walls` is .false. the callee never consults `has_*`
         ! (the test is `skip_w .and. (...)`), so folding the two BC
         ! branches into one call is exact, not merely equivalent.  Keeping
         ! them separate cost a second inlined copy of the renormaliser's
         ! `do concurrent` per branch, which is where the measured
         ! default-path regression lived.
         !
         ! The `por` actual still differs per branch because a knob-off run
         ! must NOT hand over the `(1,1,1)` placeholder — see the `por`
         ! dummy's docstring.  `h_face_left_x` is the inert stand-in: right
         ! shape, already mapped, read-only in the callee, never indexed
         ! when `use_por` is .false.
         renorm_skip_walls = (bc_w_tag == OBC_WALL .and. bc_e_tag == OBC_WALL)
         ! FOUR branches, one per (porous, closed-faces) combination: each
         ! 3-D weight has to reach the callee as a full-size, device-present
         ! explicit-shape actual, and the knob-off ones hand over a read-only
         ! PPM edge buffer as the inert stand-in (see the `por` / `open_f`
         ! docstrings).  The DEFAULT path is the last branch and is textually
         ! the call this routine has always made, so it stays bit-identical.
         if (metrics%use_porous) then
            if (metrics%use_closed_faces) then
               call renormalise_zonal_flux_to_uhbt(grid, metrics, this, ms, uhbt, dt, &
                                                   skip_walls=renorm_skip_walls, &
                                                   has_west=has_w_flux, has_east=has_e_flux, &
                                                   visc_rem=visc_rem, u_cor=u_cor, &
                                                   use_por=.true., por=metrics%por_face_area_u, &
                                                   use_open=.true., open_f=metrics%open_u)
            else
               call renormalise_zonal_flux_to_uhbt(grid, metrics, this, ms, uhbt, dt, &
                                                   skip_walls=renorm_skip_walls, &
                                                   has_west=has_w_flux, has_east=has_e_flux, &
                                                   visc_rem=visc_rem, u_cor=u_cor, &
                                                   use_por=.true., por=metrics%por_face_area_u, &
                                                   use_open=.false., open_f=this%h_face_right_x%data)
            end if
         else if (metrics%use_closed_faces) then
            call renormalise_zonal_flux_to_uhbt(grid, metrics, this, ms, uhbt, dt, &
                                                skip_walls=renorm_skip_walls, &
                                                has_west=has_w_flux, has_east=has_e_flux, &
                                                visc_rem=visc_rem, u_cor=u_cor, &
                                                use_por=.false., por=this%h_face_left_x%data, &
                                                use_open=.true., open_f=metrics%open_u)
         else
            call renormalise_zonal_flux_to_uhbt(grid, metrics, this, ms, uhbt, dt, &
                                                skip_walls=renorm_skip_walls, &
                                                has_west=has_w_flux, has_east=has_e_flux, &
                                                visc_rem=visc_rem, u_cor=u_cor, &
                                                use_por=.false., por=this%h_face_left_x%data, &
                                                use_open=.false., open_f=this%h_face_right_x%data)
         end if
      end if
   end subroutine continuity_zonal_flux

   pure subroutine renormalise_zonal_flux_to_uhbt(grid, metrics, this, ms, uhbt, dt, skip_walls, &
                                                  has_west, has_east, visc_rem, u_cor, &
                                                  use_por, por, use_open, open_f)
      !! Apply a uniform per-face velocity correction so
      !! `Σ_k mass_flux_x_layer(i, j, k) = uhbt(i, j)` at every face.
      !! Helper for `continuity_zonal_flux`.
      !! `skip_walls` (default true) bypasses the physical-wall faces
      !! where mass_flux was zeroed; pass `.false.` for periodic axes
      !! so the wall faces (which carry real transport) are renormalised.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(continuity_t), intent(in) :: this
         !! Read-only here (only the PPM edge buffers are consulted).
         !! `intent(in)` is load-bearing, not tidiness: the knob-off call
         !! sites pass one of THESE buffers as the inert `por` stand-in,
         !! and `intent(in)` turns "the callee never defines it" from a
         !! comment into a compiler-enforced invariant, so the argument
         !! association can never become aliasing.
      type(multilayer_state_t), intent(inout) :: ms
      ! assumed-shape-ok: uhbt is a face-sized array (nx+1, ny); explicit-shape
      ! would require a separate nu=nx+1 argument that the caller doesn't pass.
      real(wp), intent(in) :: uhbt(:, :)
      real(wp), intent(in) :: dt
         !! Outer-step dt (s), for the CFL bracket on `du`.
      logical, intent(in), optional :: skip_walls
         !! When .true. (default), cycle the physical-wall faces.
         !! When .false. (periodic axis), include them in the renorm.
      logical, intent(in), optional :: has_west, has_east
         !! Physical-edge flags (default .true. = single-rank behaviour,
         !! bit-identical).  .false. at an MPI seam: the local "wall
         !! position" face `i = nghost+1` (west) / `i = nghost+nx_phys+1`
         !! (east) is a REAL interior face carrying transport, so it MUST
         !! be renormalised like any interior face (O4 seam fix — the
         !! un-renormalised seam face left per-layer fluxes inconsistent
         !! with `uhbt`; the Eulerian-z h-rescale hid it in h while
         !! tracers rode the raw fluxes, breaking hTr/h at the seam).
      ! assumed-shape-ok: mirrors the `uhbt` waiver above — visc_rem is a
      ! (nx+1, ny, nz) face array the caller holds on `bt_work`, and the
      ! routine is cadence-bounded (once per stage, not per substep).
      real(wp), intent(in), optional :: visc_rem(:, :, :)
         !! Per-layer viscous remnant γ_k weighting the barotropic increment
         !! (MOM6 `u_cor = u + du·visc_rem`; Jacobian
         !! `duhdu = dy·h_marg·visc_rem`).  ABSENT ⇒ γ ≡ 1 ⇒ the historical
         !! uniform-`du` form, bit-identical.
      ! assumed-shape-ok: face array, same waiver as `uhbt`.
      real(wp), intent(inout), optional :: u_cor(:, :, :)
         !! MOM6's `u_cor` return: the
         !! transport-consistent velocity `u0 + du*gamma_k`, i.e. the velocity
         !! that yields `uhbt` as the depth-integrated transport.
         !! **This is a SEPARATE time-mean field — it must NEVER be the
         !! prognostic velocity.** MOM6 writes it to `u_av`
         !! and evaluates the slow
         !! tendencies on it; writing it back into the prognostic kills the
         !! run by step 5.  See docs/MOM6_SPLIT_RK2_SPEC.md §5 trap 1.
         !! Absent => flux-only, bit-identical.

      logical, intent(in) :: use_por
         !! Porous barriers active.  `.false.` => `por` is never indexed and
         !! the arithmetic below stays byte-identical to the un-narrowed form.
      real(wp), intent(in) :: por(grid%nx_total + 1, grid%ny_total, ms%nz_ml)
         !! Layer-averaged open-area fraction at this stagger (nondim),
         !! read ONLY when `use_por`.
         !!
         !! `use_por = .false.` callers must still pass a face-sized array
         !! that is genuinely device-present, and NOT the knob-off
         !! `(1,1,1)` placeholder: nvfortran builds the `do concurrent`
         !! data clause from the LOOP BOUNDS, not from the descriptor, so a
         !! placeholder is reported "partially present" and aborts under
         !! `mem:separate` even though the branch that indexes it is never
         !! taken.  The call sites therefore hand over one of this
         !! routine's own read-only PPM edge buffers as an inert stand-in
         !! (right shape, already mapped, never DEFINED here, so no
         !! argument aliasing) -- which keeps the two full-size open-area
         !! fields off the allocation list entirely for a default run.

      logical, intent(in) :: use_open
         !! z-level closed faces active (`&vcoord_nml zfixed_closed_faces`).
         !! `.false.` => `open_f` is never indexed and the arithmetic below
         !! stays byte-identical to the un-masked form.
      real(wp), intent(in) :: open_f(grid%nx_total + 1, grid%ny_total, ms%nz_ml)
         !! Per-layer 0/1 face-open mask at this stagger, read ONLY when
         !! `use_open`.  It enters the SAME weight `wk` the porous fraction
         !! does -- `wk = dy_cu * por * open` -- which is the whole reason
         !! the barotropic transport is distributed over the OPEN layers
         !! only: a closed layer gets `wk = 0`, so it receives no `du` and
         !! contributes nothing to `sum_h`, and `sum_k mass_flux = uhbt`
         !! stays the exact fixed point the Newton solve iterates to.
         !!
         !! `use_open = .false.` callers must still pass a face-sized,
         !! genuinely device-present array and NOT the `(1,1,1)`
         !! placeholder -- see the `por` dummy's docstring for why; the
         !! call sites hand over `h_face_right_x` as the inert stand-in.

      integer :: i, j, k, nx, ny, nz, iter
      real(wp) :: u, h_face, sum_flux, sum_h, du, target, w, wk
      real(wp) :: flux0(NZ_STACK_MAX), u0(NZ_STACK_MAX)
      real(wp) :: u_lim, du_hi, du_lo, vr_k, du_k
      real(wp) :: b_lo, b_hi, du_new
      integer :: maxit
      logical :: consistent
      logical :: skip_w, has_w, has_e, use_vr, upd_u

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml
      skip_w = .true.
      if (present(skip_walls)) skip_w = skip_walls
      has_w = .true.
      if (present(has_west)) has_w = has_west
      has_e = .true.
      if (present(has_east)) has_e = has_east
      upd_u = present(u_cor)
      ! gamma-weighting and the u_cor write-back are ONE behaviour: together they
      ! make this renormalisation the barotropic correction (MOM6 does both in a
      ! single `continuity(..., uhbt, visc_rem, u_cor)` call).  In the historical
      ! mode `apply_bt_correction` owns the gamma-weighted Delta-u and this routine
      ! must stay the unweighted flux-only renormaliser — so both are off.
      use_vr = present(visc_rem) .and. upd_u
      consistent = this%renorm_consistent_flux
      maxit = RENORM_MAXIT
      if (consistent) maxit = RENORM_MAXIT_CONSISTENT

      ! Skip the array-edge faces (i=1, i=nx+1) and, unless skip_walls
      ! is false (periodic) or the edge is an MPI seam (has_* false),
      ! the two physical walls.
      do concurrent(j=1:ny, i=2:nx) &
         local(k, iter, u, h_face, sum_flux, sum_h, du, target, w, wk, flux0, u0, &
               u_lim, du_hi, du_lo, vr_k, du_k, b_lo, b_hi, du_new)
         ! Bypass physical walls when requested (mass_flux already 0 there
         ! for wall BCs; for periodic, real transport is present).
         if (skip_w .and. &
             ((i == grid%nghost + 1 .and. has_w) .or. &
              (i == grid%nghost + grid%nx_phys + 1 .and. has_e))) cycle
         w = metrics%dy_cu(i, j)
         target = uhbt(i, j)
         do k = 1, nz
            u0(k) = ms%u_face_x_layer(i, j, k)
            flux0(k) = ms%mass_flux_x_layer(i, j, k)
         end do
         ! Legacy single-step form (wet/dry composition — see the
         ! `renorm_legacy_single_step` docstring): one linear correction
         ! with donors picked at the uncorrected velocity, bit-identical
         ! to the pre-Newton renormalisation.
         if (this%renorm_legacy_single_step) then
            sum_flux = 0.0_wp
            sum_h = 0.0_wp
            do k = 1, nz
               wk = w
               if (use_por) wk = w*por(i, j, k)
               if (use_open) wk = wk*open_f(i, j, k)
               if (u0(k) >= 0.0_wp) then
                  h_face = this%h_face_left_x%data(i, j, k)
               else
                  h_face = this%h_face_right_x%data(i, j, k)
               end if
               sum_flux = sum_flux + flux0(k)
               sum_h = sum_h + h_face*wk
            end do
            if (sum_h > 0.0_wp) then
               du = (target - sum_flux)/sum_h
               do k = 1, nz
                  wk = w
                  if (use_por) wk = w*por(i, j, k)
                  if (use_open) wk = wk*open_f(i, j, k)
                  if (u0(k) >= 0.0_wp) then
                     h_face = this%h_face_left_x%data(i, j, k)
                  else
                     h_face = this%h_face_right_x%data(i, j, k)
                  end if
                  ms%mass_flux_x_layer(i, j, k) = flux0(k) + du*h_face*wk
               end do
            end if
            cycle
         end if
         ! CFL bracket: the corrected velocity of EVERY layer must stay within
         ! RENORM_CFL.  A residual `uhbt` mismatch is preferable to handing a
         ! near-massless layer a super-CFL velocity (MOM6 does the same).
         !
         ! `max(..., H_DIV_EPS)` is a 1/0 guard and nothing else.  `idxCu` is
         ! `1/dxCu` MULTIPLIED BY `wet_u` in `metrics_apply_land_mask`, so it
         ! is EXACTLY zero on every land face — and this loop visits land
         ! faces: it skips only the array edges and the two physical walls,
         ! never the interior coastline.  Unguarded that is `0.25/0`, which
         ! traps under `-ffpe-trap=zero` and otherwise makes `u_lim = +Inf`,
         ! poisoning `du_hi`/`du_lo` with infinities on a face where the
         ! answer is not used at all (`sum_h = 0` there, the Newton residual
         ! is identically zero, and the flux stays zero because `h_face·w` is
         ! zero).  On a WET face `dt·idxCu = dt/dxCu` is a physical rate,
         ! decades above `H_DIV_EPS = 1e-20`, so the `max` selects the true
         ! operand and the result is bit-identical.  `max` rather than `+`
         ! deliberately: an additive guard perturbs the last bit once
         ! `dt/dxCu` drops near `1e-16/1e-20`, which a long-dx, short-dt
         ! configuration can reach.
         u_lim = RENORM_CFL/max(dt*metrics%idxCu(i, j), H_DIV_EPS)
         du_hi = huge(1.0_wp)
         du_lo = -huge(1.0_wp)
         do k = 1, nz
            if (use_vr) then
               ! Layer increment is `du·γ_k`, so the CFL bound on `du` is
               ! divided by γ_k.  A layer the implicit friction has fully
               ! damped (γ→0) receives no increment and imposes no bound.
               vr_k = visc_rem(i, j, k)
               if (vr_k > RENORM_VR_MIN) then
                  du_hi = min(du_hi, (u_lim - u0(k))/vr_k)
                  du_lo = max(du_lo, (-u_lim - u0(k))/vr_k)
               end if
            else
               du_hi = min(du_hi, u_lim - u0(k))
               du_lo = max(du_lo, -u_lim - u0(k))
            end if
         end do
         du_hi = max(du_hi, 0.0_wp)
         du_lo = min(du_lo, 0.0_wp)
         ! Newton on `du` with the DONOR RE-PICKED each iteration.  With
         ! `visc_rem` present the per-layer increment is `du·γ_k` and the
         ! Jacobian carries the same weight — MOM6 `duhdu = dy·h_marg·visc_rem`.
         ! `du_k` is assigned (not multiplied
         ! by 1) on the γ-absent path so the expression tree — and therefore
         ! FP contraction — is unchanged: bit-identical by construction.
         du = 0.0_wp
         sum_h = 0.0_wp
         b_lo = -huge(1.0_wp)
         b_hi = huge(1.0_wp)
         do iter = 1, maxit
            sum_flux = 0.0_wp
            sum_h = 0.0_wp
            do k = 1, nz
               wk = w
               if (use_por) wk = w*por(i, j, k)
               if (use_open) wk = wk*open_f(i, j, k)
               du_k = du
               if (use_vr) du_k = du*visc_rem(i, j, k)
               if (u0(k) + du_k >= 0.0_wp) then
                  h_face = this%h_face_left_x%data(i, j, k)
               else
                  h_face = this%h_face_right_x%data(i, j, k)
               end if
               ! `consistent`: a layer whose donor FLIPPED under the
               ! correction carries the flux of its corrected velocity
               ! through its NEW donor, `(u0+du_k)·h_face·wk` — continuous
               ! (→ 0 from both sides) where the historical
               ! `flux0 + du_k·h_face` jumps by `u0·(h_new − h_old)·wk`.
               ! Unflipped layers keep the historical expression.
               if (consistent .and. ((u0(k) + du_k >= 0.0_wp) .neqv. (u0(k) >= 0.0_wp))) then
                  sum_flux = sum_flux + (u0(k) + du_k)*h_face*wk
               else
                  sum_flux = sum_flux + (flux0(k) + du_k*h_face*wk)
               end if
               if (use_vr) then
                  sum_h = sum_h + visc_rem(i, j, k)*h_face*wk
               else
                  sum_h = sum_h + h_face*wk
               end if
            end do
            if (sum_h <= 0.0_wp) exit
            if (abs(target - sum_flux) <= RENORM_TOL*max(1.0_wp, abs(target))) exit
            if (consistent) then
               ! Monotone, continuous F(du): keep a bracket around the root
               ! and bisect whenever the Newton step leaves it (MOM6
               ! `zonal_flux_adjust`: Newton + bisection).
               if (sum_flux < target) then
                  b_lo = max(b_lo, du)
               else
                  b_hi = min(b_hi, du)
               end if
               du_new = du + (target - sum_flux)/sum_h
               if (du_new <= b_lo .or. du_new >= b_hi) then
                  if (b_lo > -huge(1.0_wp) .and. b_hi < huge(1.0_wp)) du_new = 0.5_wp*(b_lo + b_hi)
               end if
               du = min(max(du_new, du_lo), du_hi)
            else
               du = min(max(du + (target - sum_flux)/sum_h, du_lo), du_hi)
            end if
         end do
         if (sum_h > 0.0_wp) then
            do k = 1, nz
               wk = w
               if (use_por) wk = w*por(i, j, k)
               if (use_open) wk = wk*open_f(i, j, k)
               du_k = du
               if (use_vr) du_k = du*visc_rem(i, j, k)
               if (u0(k) + du_k >= 0.0_wp) then
                  h_face = this%h_face_left_x%data(i, j, k)
               else
                  h_face = this%h_face_right_x%data(i, j, k)
               end if
               if (consistent .and. ((u0(k) + du_k >= 0.0_wp) .neqv. (u0(k) >= 0.0_wp))) then
                  ms%mass_flux_x_layer(i, j, k) = (u0(k) + du_k)*h_face*wk
               else
                  ms%mass_flux_x_layer(i, j, k) = flux0(k) + du_k*h_face*wk
               end if
               ! MOM6 `u_cor(I,j,k) = u(I,j,k) + du(I)*visc_rem(I,k)`.
               ! Writing this makes the renormalisation the SOLE barotropic
               ! correction; without it the velocity and the flux carry
               ! differently-weighted corrections that can disagree.
               if (upd_u) u_cor(i, j, k) = u0(k) + du_k
            end do
         end if
      end do
   end subroutine renormalise_zonal_flux_to_uhbt

   pure subroutine continuity_meridional_flux(grid, metrics, this, ms, dt, vhbt, bc, &
                                              visc_rem, v_cor)
      !! Meridional (y-only) PPM reconstruction + per-face mass
      !! flux.  Mirror of `continuity_zonal_flux`, with the same
      !! optional `vhbt` transport-constraint renormalisation.
      !! Writes `ms%mass_flux_y_layer`.  Walls at j=1 and j=ny+1
      !! zeroed.  In the Lie split this runs *after* the zonal
      !! apply, so it reconstructs against the already-updated
      !! `ms%h_layer`.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(continuity_t), intent(inout) :: this
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: dt
         !! Time increment (s).  Used only for the swept-volume CFL when
         !! `this%vol_cfl = .true.`; ignored (bit-identical) otherwise.
      real(wp), intent(in), optional :: vhbt(:, :)
         !! Time-mean north-face transport from barotropic substep (m³/s),
         !! shape `(nx, ny+1)`.  Width-weighted (carries `dx_cv`).
      type(ocean_bc_state_t), intent(in), optional :: bc
      ! assumed-shape-ok: pure passthrough to the renormaliser.
      real(wp), intent(in), optional :: visc_rem(:, :, :)
         !! Per-layer viscous remnant gamma_k. Absent => 1, bit-identical.
      ! assumed-shape-ok: pure passthrough to the renormaliser.
      real(wp), intent(inout), optional :: v_cor(:, :, :)
         !! Forwarded MOM6 `v_cor` destination — separate time-mean field.

      integer :: i, j, k, nx, ny, nz
      integer :: bc_s_tag, bc_n_tag
      real(wp) :: dh_m1, dh_0, dh_p1, h_left, h_right, v, h_face
      real(wp) :: hm2, hm1, h0, hp1, hp2
      real(wp) :: cfl_d, curv3_d, dh_d
      logical :: do_pos, do_volcfl, do_pd
      logical :: has_s_flux, has_n_flux, renorm_skip_walls
      real(wp) :: h_min_pos, two_h_lim

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml
      do_pos = this%use_ppm_limit_pos
      do_volcfl = this%vol_cfl
      h_min_pos = this%h_min
      ! P1 positive-definite reconstruction floor (MOM6's positive-definite
      ! PPM).  Scalars copied to locals so the DC loop never walks the
      ! continuity_t descriptor.  Off ⇒ untaken branch ⇒ bit-identical.
      do_pd = this%positive_definite
      two_h_lim = 2.0_wp*this%h_lim

      do concurrent(k=1:nz, j=3:ny - 2, i=1:nx) &
         local(dh_m1, dh_0, dh_p1, h_left, h_right, hm2, hm1, h0, hp1, hp2)
         h0 = ms%h_layer(i, j, k)
         hm1 = ppm_mirror_h(ms%h_layer(i, j - 1, k), h0, metrics%wet_T(i, j - 1))
         hp1 = ppm_mirror_h(ms%h_layer(i, j + 1, k), h0, metrics%wet_T(i, j + 1))
         hm2 = ppm_mirror_h(ms%h_layer(i, j - 2, k), hm1, metrics%wet_T(i, j - 2))
         hp2 = ppm_mirror_h(ms%h_layer(i, j + 2, k), hp1, metrics%wet_T(i, j + 2))
         call ppm_limited_slope(hm2, hm1, h0, dh_m1)
         call ppm_limited_slope(hm1, h0, hp1, dh_0)
         call ppm_limited_slope(h0, hp1, hp2, dh_p1)
         dh_0 = dh_0*metrics%wet_T(i, j - 1)*metrics%wet_T(i, j)*metrics%wet_T(i, j + 1)
         h_left = 0.5_wp*(hm1 + h0) - (dh_0 - dh_m1)/6.0_wp
         h_right = 0.5_wp*(h0 + hp1) - (dh_p1 - dh_0)/6.0_wp
         call ppm_cell_limiter(h0, h_left, h_right)
         if (do_pos) call ppm_limit_pos(h0, h_left, h_right, h_min_pos)
         if (do_pd) then
            h_left = max(h_left, two_h_lim)
            h_right = max(h_right, two_h_lim)
         end if
         this%h_face_right_y%data(i, j, k) = h_left
         this%h_face_left_y%data(i, j + 1, k) = h_right
      end do
      do concurrent(k=1:nz, i=1:nx)
         this%h_face_left_y%data(i, 1, k) = ms%h_layer(i, 1, k)
         this%h_face_right_y%data(i, 1, k) = ms%h_layer(i, 1, k)
         this%h_face_left_y%data(i, 2, k) = ms%h_layer(i, 1, k)
         this%h_face_right_y%data(i, 2, k) = ms%h_layer(i, 2, k)
         this%h_face_left_y%data(i, 3, k) = ms%h_layer(i, 2, k)
         this%h_face_right_y%data(i, 3, k) = ms%h_layer(i, 2, k)
         this%h_face_right_y%data(i, ny - 1, k) = ms%h_layer(i, ny - 1, k)
         this%h_face_left_y%data(i, ny, k) = ms%h_layer(i, ny - 1, k)
         this%h_face_right_y%data(i, ny, k) = ms%h_layer(i, ny, k)
         this%h_face_left_y%data(i, ny + 1, k) = ms%h_layer(i, ny, k)
         this%h_face_right_y%data(i, ny + 1, k) = ms%h_layer(i, ny, k)
      end do

      do concurrent(k=1:nz, j=2:ny, i=1:nx) &
         local(v, h_face, cfl_d, curv3_d, dh_d)
         v = ms%v_face_y_layer(i, j, k)
         if (v >= 0.0_wp) then
            ! Donor = cell j-1; downwind (north) edge = h_face_left_y(j).
            h_face = this%h_face_left_y%data(i, j, k)
            if (do_volcfl) then
               ! Swept-oriented donor edge diff dh = h_S - h_N (south-north).
               dh_d = this%h_face_right_y%data(i, j - 1, k) - h_face
               curv3_d = this%h_face_right_y%data(i, j - 1, k) + h_face &
                         - 2.0_wp*ms%h_layer(i, j - 1, k)
               cfl_d = v*dt*metrics%dx_cv(i, j)*metrics%iareaT(i, j - 1)
               h_face = volcfl_face(h_face, dh_d, curv3_d, cfl_d)
            end if
         else
            ! Donor = cell j; downwind (south) edge = h_face_right_y(j).
            h_face = this%h_face_right_y%data(i, j, k)
            if (do_volcfl) then
               ! Swept-oriented donor edge diff dh = h_N - h_S (north-south).
               dh_d = this%h_face_left_y%data(i, j + 1, k) - h_face
               curv3_d = h_face + this%h_face_left_y%data(i, j + 1, k) &
                         - 2.0_wp*ms%h_layer(i, j, k)
               cfl_d = (-v)*dt*metrics%dx_cv(i, j)*metrics%iareaT(i, j)
               h_face = volcfl_face(h_face, dh_d, curv3_d, cfl_d)
            end if
         end if
         ms%mass_flux_y_layer(i, j, k) = v*h_face*metrics%dx_cv(i, j)
      end do
      do concurrent(k=1:nz, i=1:nx)
         ms%mass_flux_y_layer(i, 1, k) = 0.0_wp
         ms%mass_flux_y_layer(i, ny + 1, k) = 0.0_wp
      end do
      ! Physical wall zeroing — see `continuity_zonal_flux` for the
      ! detailed rationale.  Without this, v_face_y at the physical
      ! south/north walls picks up Coriolis / wind contributions and
      ! the PPM transport leaks tracer through to ghost cells.
      ! OBC dispatch matches the zonal helper — see header.
      ! Under MPI decomposition bc%has_south / bc%has_north gates the
      ! zeroing so a subdomain seam face is left for the halo exchange.
      bc_s_tag = OBC_WALL
      bc_n_tag = OBC_WALL
      ! Physical-edge flags: default .true. => single-rank bit-identical
      ! behaviour; .false. at an MPI seam skips the hard zero (O0, D0).
      has_s_flux = .true.
      has_n_flux = .true.
      if (present(bc)) then
         bc_s_tag = bc%south%bc_type
         bc_n_tag = bc%north%bc_type
         has_s_flux = bc%has_south
         has_n_flux = bc%has_north
      end if
      do concurrent(k=1:nz, i=1:nx)
         if (bc_s_tag == OBC_WALL .and. has_s_flux) ms%mass_flux_y_layer(i, grid%nghost + 1, k) = 0.0_wp
         if (bc_n_tag == OBC_WALL .and. has_n_flux) ms%mass_flux_y_layer(i, grid%nghost + grid%ny_phys + 1, k) = 0.0_wp
      end do

      ! ---- Porous barriers: see the zonal twin (incl. why it is inline) ----
      if (metrics%use_porous) then
         do concurrent(k=1:nz, j=1:ny + 1, i=1:nx)
            ms%mass_flux_y_layer(i, j, k) = ms%mass_flux_y_layer(i, j, k)* &
                                            metrics%por_face_area_v(i, j, k)
         end do
      end if

      ! ---- z-level closed faces: see the zonal twin ----
      if (metrics%use_closed_faces) then
         do concurrent(k=1:nz, j=1:ny + 1, i=1:nx)
            ms%mass_flux_y_layer(i, j, k) = ms%mass_flux_y_layer(i, j, k)* &
                                            metrics%open_v(i, j, k)
         end do
      end if

      ! MOM6-style transport constraint — see continuity_zonal_flux
      ! for the rationale.
      if (present(vhbt)) then
         ! Any non-WALL meridional edge (periodic OR open-class) carries
         ! genuine transport at the physical-wall face and must be
         ! renormalised — see continuity_zonal_flux for the full rationale.
         ! Genuine closed WALL is skipped (default); all-WALL stays
         ! bit-identical to pre-OBC behaviour.
         ! One call site per `por` actual — see the zonal twin for why
         ! folding the two BC branches together is exact and why it
         ! matters for the default-path cost.
         renorm_skip_walls = (bc_s_tag == OBC_WALL .and. bc_n_tag == OBC_WALL)
         ! Four branches — see the zonal twin for why.
         if (metrics%use_porous) then
            if (metrics%use_closed_faces) then
               call renormalise_meridional_flux_to_vhbt(grid, metrics, this, ms, vhbt, dt, &
                                                        skip_walls=renorm_skip_walls, &
                                                        has_south=has_s_flux, has_north=has_n_flux, &
                                                        visc_rem=visc_rem, v_cor=v_cor, &
                                                        use_por=.true., por=metrics%por_face_area_v, &
                                                        use_open=.true., open_f=metrics%open_v)
            else
               call renormalise_meridional_flux_to_vhbt(grid, metrics, this, ms, vhbt, dt, &
                                                        skip_walls=renorm_skip_walls, &
                                                        has_south=has_s_flux, has_north=has_n_flux, &
                                                        visc_rem=visc_rem, v_cor=v_cor, &
                                                        use_por=.true., por=metrics%por_face_area_v, &
                                                        use_open=.false., open_f=this%h_face_right_y%data)
            end if
         else if (metrics%use_closed_faces) then
            call renormalise_meridional_flux_to_vhbt(grid, metrics, this, ms, vhbt, dt, &
                                                     skip_walls=renorm_skip_walls, &
                                                     has_south=has_s_flux, has_north=has_n_flux, &
                                                     visc_rem=visc_rem, v_cor=v_cor, &
                                                     use_por=.false., por=this%h_face_left_y%data, &
                                                     use_open=.true., open_f=metrics%open_v)
         else
            call renormalise_meridional_flux_to_vhbt(grid, metrics, this, ms, vhbt, dt, &
                                                     skip_walls=renorm_skip_walls, &
                                                     has_south=has_s_flux, has_north=has_n_flux, &
                                                     visc_rem=visc_rem, v_cor=v_cor, &
                                                     use_por=.false., por=this%h_face_left_y%data, &
                                                     use_open=.false., open_f=this%h_face_right_y%data)
         end if
      end if
   end subroutine continuity_meridional_flux

   pure subroutine renormalise_meridional_flux_to_vhbt(grid, metrics, this, ms, vhbt, dt, skip_walls, &
                                                       has_south, has_north, visc_rem, v_cor, &
                                                       use_por, por, use_open, open_f)
      !! Apply a uniform per-face velocity correction so
      !! `Σ_k mass_flux_y_layer(i, j, k) = vhbt(i, j)` at every face.
      !! Mirror of `renormalise_zonal_flux_to_uhbt`.  See that routine
      !! for the `skip_walls` and `has_*` (MPI seam, O4 fix) semantics.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(continuity_t), intent(in) :: this
         !! Read-only here (only the PPM edge buffers are consulted).
         !! `intent(in)` is load-bearing, not tidiness: the knob-off call
         !! sites pass one of THESE buffers as the inert `por` stand-in,
         !! and `intent(in)` turns "the callee never defines it" from a
         !! comment into a compiler-enforced invariant, so the argument
         !! association can never become aliasing.
      type(multilayer_state_t), intent(inout) :: ms
      ! assumed-shape-ok: vhbt is a face-sized array (nx, ny+1); explicit-shape
      ! would require a separate nv=ny+1 argument that the caller doesn't pass.
      real(wp), intent(in) :: vhbt(:, :)
      real(wp), intent(in) :: dt
         !! Outer-step dt (s), for the CFL bracket on `dv`.
      logical, intent(in), optional :: skip_walls
      logical, intent(in), optional :: has_south, has_north
         !! Physical-edge flags (default .true. = single-rank behaviour,
         !! bit-identical); .false. at a y-decomposition seam forces the
         !! local wall-position face to be renormalised as interior.
      ! assumed-shape-ok: mirrors the `vhbt` waiver above.
      real(wp), intent(in), optional :: visc_rem(:, :, :)
         !! Per-layer viscous remnant γ_k. Absent ⇒ γ ≡ 1, bit-identical.
         !! See `renormalise_zonal_flux_to_uhbt` for the full rationale.
      ! assumed-shape-ok: face array, same waiver as `vhbt`.
      real(wp), intent(inout), optional :: v_cor(:, :, :)
         !! MOM6 `v_cor`.  Separate time-mean field, NEVER the prognostic —
         !! see the zonal twin and docs/MOM6_SPLIT_RK2_SPEC.md §5 trap 1.

      logical, intent(in) :: use_por
         !! Porous barriers active.  `.false.` => `por` is never indexed and
         !! the arithmetic below stays byte-identical to the un-narrowed form.
      real(wp), intent(in) :: por(grid%nx_total, grid%ny_total + 1, ms%nz_ml)
         !! Layer-averaged open-area fraction at this stagger (nondim),
         !! read ONLY when `use_por`.
         !!
         !! `use_por = .false.` callers must still pass a face-sized array
         !! that is genuinely device-present, and NOT the knob-off
         !! `(1,1,1)` placeholder: nvfortran builds the `do concurrent`
         !! data clause from the LOOP BOUNDS, not from the descriptor, so a
         !! placeholder is reported "partially present" and aborts under
         !! `mem:separate` even though the branch that indexes it is never
         !! taken.  The call sites therefore hand over one of this
         !! routine's own read-only PPM edge buffers as an inert stand-in
         !! (right shape, already mapped, never DEFINED here, so no
         !! argument aliasing) -- which keeps the two full-size open-area
         !! fields off the allocation list entirely for a default run.

      logical, intent(in) :: use_open
         !! z-level closed faces active.  `.false.` => `open_f` is never
         !! indexed; byte-identical to the un-masked form.  See the zonal
         !! twin for the full rationale.
      real(wp), intent(in) :: open_f(grid%nx_total, grid%ny_total + 1, ms%nz_ml)
         !! Per-layer 0/1 face-open mask at this stagger, read ONLY when
         !! `use_open`; enters the SAME weight `wk = dx_cv * por * open`.
         !! Knob-off callers pass `h_face_right_y` as the inert stand-in
         !! (never the `(1,1,1)` placeholder).

      integer :: i, j, k, nx, ny, nz, iter
      real(wp) :: v, h_face, sum_flux, sum_h, dv, target, w, wk
      real(wp) :: flux0(NZ_STACK_MAX), v0(NZ_STACK_MAX)
      real(wp) :: v_lim, dv_hi, dv_lo, vr_k, dv_k
      real(wp) :: b_lo, b_hi, dv_new
      integer :: maxit
      logical :: consistent
      logical :: skip_w, has_s, has_n, use_vr, upd_v

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml
      skip_w = .true.
      if (present(skip_walls)) skip_w = skip_walls
      has_s = .true.
      if (present(has_south)) has_s = has_south
      has_n = .true.
      if (present(has_north)) has_n = has_north
      upd_v = present(v_cor)
      ! See the zonal twin: gamma-weighting and v_cor act together.
      use_vr = present(visc_rem) .and. upd_v
      consistent = this%renorm_consistent_flux
      maxit = RENORM_MAXIT
      if (consistent) maxit = RENORM_MAXIT_CONSISTENT

      do concurrent(j=2:ny, i=1:nx) &
         local(k, iter, v, h_face, sum_flux, sum_h, dv, target, w, wk, flux0, v0, v_lim, dv_hi, dv_lo, &
               vr_k, dv_k, b_lo, b_hi, dv_new)
         if (skip_w .and. &
             ((j == grid%nghost + 1 .and. has_s) .or. &
              (j == grid%nghost + grid%ny_phys + 1 .and. has_n))) cycle
         w = metrics%dx_cv(i, j)
         target = vhbt(i, j)
         do k = 1, nz
            v0(k) = ms%v_face_y_layer(i, j, k)
            flux0(k) = ms%mass_flux_y_layer(i, j, k)
         end do
         ! Legacy single-step form (wet/dry composition — mirror of the
         ! zonal branch; see `renorm_legacy_single_step`).
         if (this%renorm_legacy_single_step) then
            sum_flux = 0.0_wp
            sum_h = 0.0_wp
            do k = 1, nz
               wk = w
               if (use_por) wk = w*por(i, j, k)
               if (use_open) wk = wk*open_f(i, j, k)
               if (v0(k) >= 0.0_wp) then
                  h_face = this%h_face_left_y%data(i, j, k)
               else
                  h_face = this%h_face_right_y%data(i, j, k)
               end if
               sum_flux = sum_flux + flux0(k)
               sum_h = sum_h + h_face*wk
            end do
            if (sum_h > 0.0_wp) then
               dv = (target - sum_flux)/sum_h
               do k = 1, nz
                  wk = w
                  if (use_por) wk = w*por(i, j, k)
                  if (use_open) wk = wk*open_f(i, j, k)
                  if (v0(k) >= 0.0_wp) then
                     h_face = this%h_face_left_y%data(i, j, k)
                  else
                     h_face = this%h_face_right_y%data(i, j, k)
                  end if
                  ms%mass_flux_y_layer(i, j, k) = flux0(k) + dv*h_face*wk
               end do
            end if
            cycle
         end if
         ! CFL bracket, mirror of the zonal routine — including its
         ! land-face 1/0 guard: `idyCv` is zeroed at `wet_v == 0` by
         ! `metrics_apply_land_mask`, exactly as `idxCu` is at `wet_u == 0`.
         v_lim = RENORM_CFL/max(dt*metrics%idyCv(i, j), H_DIV_EPS)
         dv_hi = huge(1.0_wp)
         dv_lo = -huge(1.0_wp)
         do k = 1, nz
            if (use_vr) then
               vr_k = visc_rem(i, j, k)
               if (vr_k > RENORM_VR_MIN) then
                  dv_hi = min(dv_hi, (v_lim - v0(k))/vr_k)
                  dv_lo = max(dv_lo, (-v_lim - v0(k))/vr_k)
               end if
            else
               dv_hi = min(dv_hi, v_lim - v0(k))
               dv_lo = max(dv_lo, -v_lim - v0(k))
            end if
         end do
         dv_hi = max(dv_hi, 0.0_wp)
         dv_lo = min(dv_lo, 0.0_wp)
         ! Newton on `dv` with the DONOR RE-PICKED each iteration (see the
         ! RENORM_MAXIT docstring; mirror of the zonal routine).
         dv = 0.0_wp
         sum_h = 0.0_wp
         b_lo = -huge(1.0_wp)
         b_hi = huge(1.0_wp)
         do iter = 1, maxit
            sum_flux = 0.0_wp
            sum_h = 0.0_wp
            do k = 1, nz
               wk = w
               if (use_por) wk = w*por(i, j, k)
               if (use_open) wk = wk*open_f(i, j, k)
               dv_k = dv
               if (use_vr) dv_k = dv*visc_rem(i, j, k)
               if (v0(k) + dv_k >= 0.0_wp) then
                  h_face = this%h_face_left_y%data(i, j, k)
               else
                  h_face = this%h_face_right_y%data(i, j, k)
               end if
               ! `consistent`: a layer whose donor FLIPPED under the
               ! correction carries the flux of its corrected velocity
               ! through its NEW donor, `(v0+dv_k)·h_face·wk` — continuous
               ! (→ 0 from both sides) where the historical
               ! `flux0 + dv_k·h_face` jumps by `v0·(h_new − h_old)·wk`.
               ! Unflipped layers keep the historical expression.
               if (consistent .and. ((v0(k) + dv_k >= 0.0_wp) .neqv. (v0(k) >= 0.0_wp))) then
                  sum_flux = sum_flux + (v0(k) + dv_k)*h_face*wk
               else
                  sum_flux = sum_flux + (flux0(k) + dv_k*h_face*wk)
               end if
               if (use_vr) then
                  sum_h = sum_h + visc_rem(i, j, k)*h_face*wk
               else
                  sum_h = sum_h + h_face*wk
               end if
            end do
            if (sum_h <= 0.0_wp) exit
            if (abs(target - sum_flux) <= RENORM_TOL*max(1.0_wp, abs(target))) exit
            if (consistent) then
               ! Monotone, continuous F(dv): keep a bracket around the root
               ! and bisect whenever the Newton step leaves it (MOM6
               ! `zonal_flux_adjust`: Newton + bisection).
               if (sum_flux < target) then
                  b_lo = max(b_lo, dv)
               else
                  b_hi = min(b_hi, dv)
               end if
               dv_new = dv + (target - sum_flux)/sum_h
               if (dv_new <= b_lo .or. dv_new >= b_hi) then
                  if (b_lo > -huge(1.0_wp) .and. b_hi < huge(1.0_wp)) dv_new = 0.5_wp*(b_lo + b_hi)
               end if
               dv = min(max(dv_new, dv_lo), dv_hi)
            else
               dv = min(max(dv + (target - sum_flux)/sum_h, dv_lo), dv_hi)
            end if
         end do
         if (sum_h > 0.0_wp) then
            do k = 1, nz
               wk = w
               if (use_por) wk = w*por(i, j, k)
               if (use_open) wk = wk*open_f(i, j, k)
               dv_k = dv
               if (use_vr) dv_k = dv*visc_rem(i, j, k)
               if (v0(k) + dv_k >= 0.0_wp) then
                  h_face = this%h_face_left_y%data(i, j, k)
               else
                  h_face = this%h_face_right_y%data(i, j, k)
               end if
               if (consistent .and. ((v0(k) + dv_k >= 0.0_wp) .neqv. (v0(k) >= 0.0_wp))) then
                  ms%mass_flux_y_layer(i, j, k) = (v0(k) + dv_k)*h_face*wk
               else
                  ms%mass_flux_y_layer(i, j, k) = flux0(k) + dv_k*h_face*wk
               end if
               ! MOM6 `v_cor = v + dv·visc_rem` — see the zonal twin.
               if (upd_v) v_cor(i, j, k) = v0(k) + dv_k
            end do
         end if
      end do
   end subroutine renormalise_meridional_flux_to_vhbt

   pure subroutine continuity_apply_zonal(grid, metrics, ms, dt, h_min)
      !! Apply the zonal (x-flux) thickness update:
      !!   h(i,j,k) ← h(i,j,k) - dt · (Φx(i+1,j,k) - Φx(i,j,k)) · iareaT
      !! Overwrites `flux_h_layer` with the x-divergence so the
      !! companion meridional apply can accumulate the total.  Φx is
      !! the width-weighted transport (m³/s); `iareaT` closes the
      !! divergence to a per-area rate (= `inv_dx` on uniform metrics).
      !!
      !! Optional `h_min` (m): when > 0, applies max(h_new, h_min) on the
      !! h-update (Phase-1 Lagrangian floor). `mass_budget_continuity` records
      !! the RAW divergence regardless — the floor injection shows up in the
      !! mass Error diagnostic (R2). 0 or absent ⇒ bit-identical.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: dt
      real(wp), intent(in), optional :: h_min
         !! Minimum-thickness floor (m). 0 or absent ⇒ off ⇒ bit-identical.

      integer :: i, j, k, nx, ny, nz
      real(wp) :: div_x, h_min_use

      h_min_use = 0.0_wp
      if (present(h_min)) h_min_use = h_min

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml

      if (h_min_use > 0.0_wp) then
         do concurrent(k=1:nz, j=1:ny, i=1:nx)
            div_x = (ms%mass_flux_x_layer(i + 1, j, k) - &
                     ms%mass_flux_x_layer(i, j, k))*metrics%iareaT(i, j)
            ms%flux_h_layer(i, j, k) = div_x
            ms%h_layer(i, j, k) = max(ms%h_layer(i, j, k) - dt*div_x, h_min_use)
            ms%mass_budget_continuity(i, j, k) = &
               ms%mass_budget_continuity(i, j, k) - dt*div_x
         end do
      else
         do concurrent(k=1:nz, j=1:ny, i=1:nx)
            div_x = (ms%mass_flux_x_layer(i + 1, j, k) - &
                     ms%mass_flux_x_layer(i, j, k))*metrics%iareaT(i, j)
            ms%flux_h_layer(i, j, k) = div_x
            ms%h_layer(i, j, k) = ms%h_layer(i, j, k) - dt*div_x
            ms%mass_budget_continuity(i, j, k) = &
               ms%mass_budget_continuity(i, j, k) - dt*div_x
         end do
      end if
   end subroutine continuity_apply_zonal

   pure subroutine continuity_apply_meridional(grid, metrics, ms, dt, h_min)
      !! Apply the meridional (y-flux) thickness update on top of
      !! the zonally-updated state:
      !!   h(i,j,k) ← h(i,j,k) - dt · (Φy(i,j+1,k) - Φy(i,j,k)) · iareaT
      !! Adds the y-divergence to `flux_h_layer` so the field ends
      !! the split step holding the *total* horizontal divergence
      !! that the vertical-advection kernel consumes
      !! (`w_interface(k+1) = w(k) - flux_h_layer(k)`).  Φy carries
      !! `dx_cv`; `iareaT` = `inv_dy` on uniform metrics.
      !!
      !! Optional `h_min` (m): when > 0, applies max(h_new, h_min) on the
      !! h-update (Phase-1 Lagrangian floor). `mass_budget_continuity` records
      !! the RAW divergence regardless — the floor injection shows up in the
      !! mass Error diagnostic (R2). 0 or absent ⇒ bit-identical.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: dt
      real(wp), intent(in), optional :: h_min
         !! Minimum-thickness floor (m). 0 or absent ⇒ off ⇒ bit-identical.

      integer :: i, j, k, nx, ny, nz
      real(wp) :: div_y, h_min_use

      h_min_use = 0.0_wp
      if (present(h_min)) h_min_use = h_min

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml

      if (h_min_use > 0.0_wp) then
         do concurrent(k=1:nz, j=1:ny, i=1:nx)
            div_y = (ms%mass_flux_y_layer(i, j + 1, k) - &
                     ms%mass_flux_y_layer(i, j, k))*metrics%iareaT(i, j)
            ms%flux_h_layer(i, j, k) = ms%flux_h_layer(i, j, k) + div_y
            ms%h_layer(i, j, k) = max(ms%h_layer(i, j, k) - dt*div_y, h_min_use)
            ms%mass_budget_continuity(i, j, k) = &
               ms%mass_budget_continuity(i, j, k) - dt*div_y
         end do
      else
         do concurrent(k=1:nz, j=1:ny, i=1:nx)
            div_y = (ms%mass_flux_y_layer(i, j + 1, k) - &
                     ms%mass_flux_y_layer(i, j, k))*metrics%iareaT(i, j)
            ms%flux_h_layer(i, j, k) = ms%flux_h_layer(i, j, k) + div_y
            ms%h_layer(i, j, k) = ms%h_layer(i, j, k) - dt*div_y
            ms%mass_budget_continuity(i, j, k) = &
               ms%mass_budget_continuity(i, j, k) - dt*div_y
         end do
      end if
   end subroutine continuity_apply_meridional

   subroutine continuity_step_split(grid, metrics, this, ms, dt)
      !! **Test-only** (no production caller): continuity-only split wrapper;
      !! production runs the tracer-interleaved `continuity_tracer_step_split`.
      !! Directionally-split (Lie) PPM continuity step over `dt`:
      !!
      !!   1. zonal flux Φx from current h
      !!   2. apply x-divergence:  h ← h - dt·∂Φx/∂x
      !!   3. meridional flux Φy from the *updated* h
      !!   4. apply y-divergence:  h ← h - dt·∂Φy/∂y
      !!
      !! Standard Lie-split PPM continuity operator order.
      !! Compared to the unsplit
      !! `continuity_compute_fluxes` +
      !! `continuity_apply_fluxes` pair, this form
      !! relaxes the 2D-combined CFL constraint to the per-
      !! direction CFL ≤ 1.  For uniform `h_layer` and divergence-
      !! free flow the two forms are bit-identical (no advection of
      !! variation between substeps); they diverge under non-trivial
      !! advection at O(dt²·||∇h||) — same order as a vanilla Strang
      !! vs Lie split.
      !!
      !! `flux_h_layer` ends the step holding the *total* horizontal
      !! divergence summed over both substeps, ready for the
      !! vertical-advection consumer.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(continuity_t), intent(inout) :: this
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: dt

      call continuity_zonal_flux(grid, metrics, this, ms, dt)
      call continuity_apply_zonal(grid, metrics, ms, dt)
      call continuity_meridional_flux(grid, metrics, this, ms, dt)
      call continuity_apply_meridional(grid, metrics, ms, dt)
   end subroutine continuity_step_split

   pure subroutine tracer_advect_zonal(grid, metrics, this, ms, dt, bc)
      !! Zonal half of the direction-split tracer advection.
      !! Mirrors `tracer_advect` but updates `hTr` using
      !! only the x-direction tracer mass flux.  Companion to
      !! `tracer_advect_meridional`.  Both are called
      !! interleaved with the continuity substeps by
      !! `continuity_tracer_step_split` to preserve CWC.
      !!
      !! OBC inflow handling: when `bc%<edge>%bc_type == OBC_CLAMPED`
      !! and `bc%<edge>%clamped_tracer(it)` carries a prescribed value,
      !! the ghost-cell `hTr` is set so the downstream upwind pick
      !! reads the prescribed boundary value on inflow.  Override is
      !! a no-op for OBC_WALL (no flux crosses the face anyway).
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(continuity_t), intent(inout) :: this
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: dt
      type(ocean_bc_state_t), intent(in), optional :: bc

      integer :: it, i_w_ghost, i_e_ghost, i, j, k
      real(wp) :: clamped_tr

      if (.not. allocated(ms%tracers)) return

      ! Set ALL ghost columns (i = 1..nghost) on the west side to the
      ! clamped value so the impl's hard-coded boundary-halo writes
      ! (Tr_face_*_x(1..3, ...) all read hTr(1) or hTr(2)) all see the
      ! prescribed value, regardless of which is upwind.
      i_w_ghost = grid%nghost
      i_e_ghost = grid%nghost + grid%nx_phys + 1

      do it = 1, size(ms%tracers)
         if (.not. ms%tracers(it)%do_horizontal_advection) cycle
         if (present(bc)) then
            if (bc%west%bc_type == OBC_CLAMPED .and. &
                allocated(bc%west%clamped_tracer) .and. &
                it <= size(bc%west%clamped_tracer)) then
               clamped_tr = bc%west%clamped_tracer(it)
               do concurrent(k=1:ms%nz_ml, j=1:grid%ny_total, i=1:grid%nghost)
                  ms%tracers(it)%hTr(i, j, k) = clamped_tr*ms%h_layer(i, j, k)
               end do
            end if
            if (bc%east%bc_type == OBC_CLAMPED .and. &
                allocated(bc%east%clamped_tracer) .and. &
                it <= size(bc%east%clamped_tracer)) then
               clamped_tr = bc%east%clamped_tracer(it)
               do concurrent(k=1:ms%nz_ml, j=1:grid%ny_total, &
                             i=grid%nghost + grid%nx_phys + 1:grid%nx_total)
                  ms%tracers(it)%hTr(i, j, k) = clamped_tr*ms%h_layer(i, j, k)
               end do
            end if
         end if
         ! Budget dispatch: pass the matching budget accumulator for whichever
         ! tracer opted into a budget slot (heat/salt); age/passive tracers
         ! (budget_id = NONE) get no budget arg (absent ⇒ inert).
         select case (ms%tracers(it)%budget_id)
         case (TRACER_BUDGET_HEAT)
            call tracer_advect_zonal_one_impl( &
               grid%nx_total, grid%ny_total, ms%nz_ml, &
               dt, metrics%iareaT, metrics%wet_T, &
               ms%h_layer, &
               ms%tracers(it)%hTr, &
               ms%mass_flux_x_layer, &
               this%h_face_left_x%data, &
               this%h_face_right_x%data, &
               budget_adv=ms%heat_budget_horiz_adv)
         case (TRACER_BUDGET_SALT)
            call tracer_advect_zonal_one_impl( &
               grid%nx_total, grid%ny_total, ms%nz_ml, &
               dt, metrics%iareaT, metrics%wet_T, &
               ms%h_layer, &
               ms%tracers(it)%hTr, &
               ms%mass_flux_x_layer, &
               this%h_face_left_x%data, &
               this%h_face_right_x%data, &
               budget_adv=ms%salt_budget_horiz_adv)
         case default
            call tracer_advect_zonal_one_impl( &
               grid%nx_total, grid%ny_total, ms%nz_ml, &
               dt, metrics%iareaT, metrics%wet_T, &
               ms%h_layer, &
               ms%tracers(it)%hTr, &
               ms%mass_flux_x_layer, &
               this%h_face_left_x%data, &
               this%h_face_right_x%data)
         end select
      end do
   end subroutine tracer_advect_zonal

   pure subroutine tracer_advect_meridional(grid, metrics, this, ms, dt, bc)
      !! Meridional half of the direction-split tracer advection.
      !! Reads the post-zonal `h_layer` (since
      !! `continuity_apply_zonal` has already updated h
      !! in the split flow) and the just-computed
      !! `mass_flux_y_layer` from
      !! `continuity_meridional_flux`.
      !!
      !! Same OBC_CLAMPED ghost-cell override as the zonal helper —
      !! see `tracer_advect_zonal` for the rationale.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(continuity_t), intent(inout) :: this
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: dt
      type(ocean_bc_state_t), intent(in), optional :: bc

      integer :: it, j_s_ghost, j_n_ghost, i, j, k
      real(wp) :: clamped_tr

      if (.not. allocated(ms%tracers)) return

      j_s_ghost = grid%nghost
      j_n_ghost = grid%nghost + grid%ny_phys + 1

      do it = 1, size(ms%tracers)
         if (.not. ms%tracers(it)%do_horizontal_advection) cycle
         if (present(bc)) then
            if (bc%south%bc_type == OBC_CLAMPED .and. &
                allocated(bc%south%clamped_tracer) .and. &
                it <= size(bc%south%clamped_tracer)) then
               clamped_tr = bc%south%clamped_tracer(it)
               do concurrent(k=1:ms%nz_ml, i=1:grid%nx_total, j=1:grid%nghost)
                  ms%tracers(it)%hTr(i, j, k) = clamped_tr*ms%h_layer(i, j, k)
               end do
            end if
            if (bc%north%bc_type == OBC_CLAMPED .and. &
                allocated(bc%north%clamped_tracer) .and. &
                it <= size(bc%north%clamped_tracer)) then
               clamped_tr = bc%north%clamped_tracer(it)
               do concurrent(k=1:ms%nz_ml, i=1:grid%nx_total, &
                             j=grid%nghost + grid%ny_phys + 1:grid%ny_total)
                  ms%tracers(it)%hTr(i, j, k) = clamped_tr*ms%h_layer(i, j, k)
               end do
            end if
         end if
         ! Budget dispatch: pass the matching budget accumulator for whichever
         ! tracer opted into a budget slot (heat/salt); age/passive tracers
         ! (budget_id = NONE) get no budget arg (absent ⇒ inert).
         select case (ms%tracers(it)%budget_id)
         case (TRACER_BUDGET_HEAT)
            call tracer_advect_meridional_one_impl( &
               grid%nx_total, grid%ny_total, ms%nz_ml, &
               dt, metrics%iareaT, metrics%wet_T, &
               ms%h_layer, &
               ms%tracers(it)%hTr, &
               ms%mass_flux_y_layer, &
               this%h_face_left_y%data, &
               this%h_face_right_y%data, &
               budget_adv=ms%heat_budget_horiz_adv)
         case (TRACER_BUDGET_SALT)
            call tracer_advect_meridional_one_impl( &
               grid%nx_total, grid%ny_total, ms%nz_ml, &
               dt, metrics%iareaT, metrics%wet_T, &
               ms%h_layer, &
               ms%tracers(it)%hTr, &
               ms%mass_flux_y_layer, &
               this%h_face_left_y%data, &
               this%h_face_right_y%data, &
               budget_adv=ms%salt_budget_horiz_adv)
         case default
            call tracer_advect_meridional_one_impl( &
               grid%nx_total, grid%ny_total, ms%nz_ml, &
               dt, metrics%iareaT, metrics%wet_T, &
               ms%h_layer, &
               ms%tracers(it)%hTr, &
               ms%mass_flux_y_layer, &
               this%h_face_left_y%data, &
               this%h_face_right_y%data)
         end select
      end do
   end subroutine tracer_advect_meridional

   subroutine continuity_tracer_step_split(grid, metrics, this, ms, dt, uhbt, vhbt, bc, mle, mle_fold_active, &
                                           tracer_mode, gm, h_min, visc_rem_u, visc_rem_v, u_cor, v_cor)
      !! Production entry point for the directionally-split
      !! continuity + tracer step.  Interleaves the two so the
      !! CWC discrete theorem holds in the split form:
      !!
      !!   1. zonal_flux        — Φx from h^n
      !!   2. tracer_advect_zonal — hTr ← hTr - dt·∂(Φx·T)/∂x at h^n
      !!   3. apply_zonal       — h ← h^n - dt·∂Φx/∂x  (= h^*)
      !!   4. meridional_flux   — Φy from h^*
      !!   5. tracer_advect_meridional — hTr ← hTr - dt·∂(Φy·T)/∂y at h^*
      !!   6. apply_meridional  — h ← h^* - dt·∂Φy/∂y  (= h^{n+1})
      !!
      !! Uniform T preserved: after step 2, hTr = (h - dt·div_x)·T;
      !! after step 3, h = h - dt·div_x, so hTr/h = T still.  After
      !! step 5, hTr = (h^* - dt·div_y)·T = h^{n+1}·T.  After step 6,
      !! hTr/h = T.  Same CWC theorem as the unsplit form, lifted
      !! per direction.
      !!
      !! `flux_h_layer` ends the step holding the total horizontal
      !! divergence (sum of x and y substeps) — that's what the
      !! vertical-advection kernel consumes for w_interface.
      !!
      !! Optional `uhbt, vhbt`: time-mean barotropic-substep transports.  When
      !! supplied, the per-layer mass fluxes are renormalised so
      !! `Σ_k Φx_k = uhbt` and `Σ_k Φy_k = vhbt`, making the slow
      !! continuity advance `h_layer` consistently with the fast
      !! loop's `η_end` — MOM6's split-explicit pattern.  The same
      !! constrained fluxes feed tracer advection, so per-column
      !! `T = hTr/h` stays uniform under the constraint.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(continuity_t), intent(inout) :: this
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: dt
      real(wp), intent(in), optional :: uhbt(:, :)
      real(wp), intent(in), optional :: vhbt(:, :)
      type(ocean_bc_state_t), intent(in), optional :: bc
         !! When present, per-edge OBC tags gate the wall-zero step
         !! inside the flux kernels.  OBC_WALL keeps the Phase 3
         !! closure; OBC_OPEN (and other non-wall tags) leaves the
         !! computed mass flux at the wall face for the downstream
         !! transport.  Absent ⇒ closed-wall everywhere.
      type(ocean_mle_t), intent(in), optional :: mle
         !! Fox-Kemper mixed-layer-eddy transports (B5).  When present
         !! and enabled, `mle%uhml`/`vhml` are folded into the per-layer
         !! mass fluxes AFTER each direction's flux fill and BEFORE the
         !! matching tracer advect + divergence — so the augmented flux
         !! transports both h and tracers (conservative; velocity
         !! untouched).  Absent / disabled ⇒ bit-identical no-op.
      logical, intent(in), optional :: mle_fold_active
         !! Gates the Fox-Kemper fold to the THERMO cadence.  Absent or
         !! `.true.` ⇒ the fold applies (bit-identical default — the case
         !! at `dt_therm_ratio = 1`, where every step is a thermo step).
         !! `.false.` skips the fold so the stale FK transports (computed
         !! once per thermo interval) are NOT re-applied on the
         !! intervening non-thermo outer steps when `dt_therm_ratio > 1`.
      ! assumed-shape-ok: pure passthroughs to the renormaliser.
      real(wp), intent(in), optional :: visc_rem_u(:, :, :), visc_rem_v(:, :, :)
         !! Per-layer viscous remnant gamma_k on east / north faces.  Forwarded
         !! to the flux renormalisers, where it weights the barotropic
         !! increment (MOM6 `u_cor = u + du*visc_rem`).  Absent => gamma == 1,
         !! bit-identical.
      ! assumed-shape-ok: pure passthroughs.
      real(wp), intent(inout), optional :: u_cor(:, :, :), v_cor(:, :, :)
         !! MOM6 `u_cor`/`v_cor` destinations — the step TIME-MEAN velocity
         !! (`u_av`/`v_av`), never the prognostic.  Absent => flux-only.
      integer, intent(in), optional :: tracer_mode
         !! Phase 2 (6b) windowed-advection mode.  `TR_MODE_ADVECT`
         !! (default, absent) ⇒ the historical fused path: advance h AND
         !! advect tracers each call (bit-identical to pre-6b).
         !! `TR_MODE_ACCUMULATE` ⇒ advance h, accumulate
         !! `0.5·mass_flux·dt` into `this%uhtr/vhtr` (one += per RK2
         !! stage, weight 0.5 baked in — closes the reconstruction against
         !! the RK2-averaged h), and SKIP the per-step tracer advect so
         !! `hTr` stays frozen until the boundary drain.
      type(ocean_gm_t), intent(in), optional :: gm
         !! Gent-McWilliams thickness-diffusion transports (capability [2]).
         !! Folded into the per-layer mass fluxes exactly like `mle` — after
         !! each direction's flux fill and before the matching tracer advect
         !! + divergence (conservative; `Sum_k uhD = 0`).  Gated by the same
         !! `mle_fold_active` thermo-cadence flag.  Absent / disabled ⇒
         !! bit-identical no-op.
      real(wp), intent(in), optional :: h_min
         !! Phase-1 Lagrangian minimum-thickness floor (m). When > 0, passed
         !! to `continuity_apply_zonal`/`_meridional` to clamp h_new >= h_min.
         !! Absent or 0 ⇒ off ⇒ bit-identical.

      integer :: it
      logical :: per_x, per_y, do_mle_fold, fold_wall
      integer :: nx, ny, nz, nx_phys, ny_phys, nghost
      integer :: mode
      integer :: ii, jj, kk
      integer :: it_cw
      integer :: bc_w_tag, bc_e_tag, bc_s_tag, bc_n_tag
      real(wp) :: h_min_use

      mode = TR_MODE_ADVECT
      if (present(tracer_mode)) mode = tracer_mode

      h_min_use = 0.0_wp
      if (present(h_min)) h_min_use = h_min

      ! P2 positive-definite limiter: reset the per-call limited-face counter
      ! (accumulated across the zonal + meridional passes below).  Host scalar.
      this%n_limited_step = 0

      per_x = .false.
      per_y = .false.
      if (present(bc)) then
         per_x = bc%periodic_x .and. .not. ocean_halo_is_decomposed_x()
         per_y = bc%periodic_y .and. .not. ocean_halo_is_decomposed_y()
      end if
      ! Fox-Kemper fold defaults ON (bit-identical for callers that do not
      ! pass the gate); the dyn step passes `is_thermo_step()` to suppress
      ! the fold on non-thermo steps when dt_therm_ratio > 1.
      do_mle_fold = .true.
      if (present(mle_fold_active)) do_mle_fold = mle_fold_active
      ! Whether a GM/MLE bolus fold actually contributes this call.  When it
      ! does, the augmented flux must be re-closed at no-normal-flow WALL
      ! faces (the fold adds bolus transport at every face, including the
      ! physical wall the resolved flux already zeroed — otherwise the bolus
      ! bleeds tracer mass into the ghost halo across the wall).
      fold_wall = .false.
      if (present(gm)) fold_wall = fold_wall .or. gm%enable
      if (present(mle)) fold_wall = fold_wall .or. mle%enable
      fold_wall = fold_wall .and. do_mle_fold
      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml
      nx_phys = grid%nx_phys
      ny_phys = grid%ny_phys
      nghost = grid%nghost

      ! Windowed-advect concentration hold (step 1 of 2).  In
      ! TR_MODE_ACCUMULATE the horizontal tracer advect is deferred to the
      ! end-of-window drain, so `hTr` must not move — but `h_layer` does,
      ! every stage.  Snapshot the pre-continuity thickness so the paired
      ! rescale at the bottom of this routine can hold `T = hTr/h_layer`
      ! fixed instead of holding `hTr` fixed.  `hprev_work` is idle here:
      ! the drain is the only other consumer and it runs at outer-step end.
      if (mode == TR_MODE_ACCUMULATE) then
         ! First accumulate stage of a window: latch the window-start
         ! thickness the drain will use to undo the hold exactly.
         if (.not. this%hTr_holds_conc) then
            call drain_copy_3d(nx, ny, nz, ms%h_layer, this%h_win_start)
         end if
         call drain_copy_3d(nx, ny, nz, ms%h_layer, this%hprev_work)
      end if

      if (present(uhbt) .and. present(bc)) then
         call continuity_zonal_flux(grid, metrics, this, ms, dt, uhbt=uhbt, bc=bc, &
                                    visc_rem=visc_rem_u, u_cor=u_cor)
      else if (present(uhbt)) then
         call continuity_zonal_flux(grid, metrics, this, ms, dt, uhbt=uhbt, &
                                    visc_rem=visc_rem_u, u_cor=u_cor)
      else if (present(bc)) then
         call continuity_zonal_flux(grid, metrics, this, ms, dt, bc=bc)
      else
         call continuity_zonal_flux(grid, metrics, this, ms, dt)
      end if
      ! FK MLE fold (B5): add uhml into the zonal mass flux so it advects
      ! both h and tracers and enters the divergence.  No-op if absent.
      ! NOTE: this fold runs every dynamics call; `mle_compute_transports`
      ! runs only at thermo cadence.  When dt_therm_ratio > 1 the stale
      ! uhml/vhml are re-folded on non-thermo steps (over-applies FK at
      ! dynamics cadence — known limitation, safe/conservative via
      ! sum_k a(k) = 0; to be gated when sub-thermo cadence is exercised).
      if (present(mle) .and. do_mle_fold) then
         call mle_fold_x(mle, ms%mass_flux_x_layer, nx + 1, ny, nz)
      end if
      ! GM thickness-diffusion fold (capability [2]): same thermo-cadence
      ! gate as the FK fold; conservative because Sum_k uhD = 0.
      if (present(gm) .and. do_mle_fold) then
         call gm_fold_x(gm, ms%mass_flux_x_layer, nx + 1, ny, nz)
      end if
      ! No-normal-flow wall closure for the bolus folds (mirrors the
      ! resolved-flux wall zeroing in `continuity_zonal_flux`): the GM/MLE
      ! folds above add `uhD`/`uhml` at the physical WALL faces, which the
      ! resolved flux had zeroed.  Without re-zeroing, the bolus transports
      ! tracer mass across the wall into the ghost halo and the
      ! physical-domain `sum(hTr)` drifts.  BC-aware: periodic / open edges
      ! keep the folded transport.
      if (fold_wall) then
         bc_w_tag = OBC_WALL
         bc_e_tag = OBC_WALL
         if (present(bc)) then
            bc_w_tag = bc%west%bc_type
            bc_e_tag = bc%east%bc_type
         end if
         do concurrent(kk=1:nz, jj=1:ny)
            if (bc_w_tag == OBC_WALL) ms%mass_flux_x_layer(nghost + 1, jj, kk) = 0.0_wp
            if (bc_e_tag == OBC_WALL) ms%mass_flux_x_layer(nghost + nx_phys + 1, jj, kk) = 0.0_wp
         end do
      end if
      ! P2 positive-definite outflux limiter (zonal): scale the OUTGOING
      ! east-face mass fluxes so no donor drains below h_lim.  Applied to the
      ! FOLDED total (after the GM/MLE bolus folds + wall closure) so the
      ! h-apply, tracer advect, uhtr accumulation, and the corrector's
      ! `use_state_fluxes` reads all consume the SAME limited flux (D3).
      ! Off ⇒ skipped ⇒ bit-identical.
      if (this%positive_definite) then
         ! v1.1: forward u_cor so the limiter re-scales the captured
         ! transport-matched velocity (the u_av family) by the same θ as
         ! the flux — optional-forwarding propagates absence.  Keeps u_av
         ! consistent with the LIMITED fluxes the pred_corr corrector's
         ! use_state_fluxes CorAdCalc transports with (the flux↔velocity
         ! match is load-bearing).
         call pd_limit_zonal_impl(nx, ny, nz, dt, this%h_lim, metrics%iareaT, &
                                  ms%h_layer, ms%mass_flux_x_layer, &
                                  this%pd_theta%data, this%n_limited_step, &
                                  u_cor=u_cor)
      end if
      if (mode == TR_MODE_ACCUMULATE) then
         ! Windowed mode: accumulate this stage's zonal mass flux (with
         ! the RK2 0.5 weight baked in) and SKIP the tracer advect so
         ! hTr stays frozen.  hprev = h^{n+1} + div(uhtr) then closes
         ! the reconstruction against the RK2-averaged h.
         call accumulate_flux_x(nx + 1, ny, nz, dt, ms%mass_flux_x_layer, this%uhtr)
      else if (mode /= TR_MODE_NONE .and. present(bc)) then
         call tracer_advect_zonal(grid, metrics, this, ms, dt, bc=bc)
      else if (mode /= TR_MODE_NONE) then
         call tracer_advect_zonal(grid, metrics, this, ms, dt)
      end if
      if (h_min_use > 0.0_wp) then
         call continuity_apply_zonal(grid, metrics, ms, dt, h_min=h_min_use)
      else
         call continuity_apply_zonal(grid, metrics, ms, dt)
      end if

      ! Mid-Lie-split MPI seam exchange (O3 correctness fix): after the zonal
      ! apply updates h_layer (and hTr in ADVECT mode) the meridional flux
      ! reconstruction reads h ghost columns that the NEIGHBOUR rank's zonal
      ! apply has updated — but those ghosts were never re-exchanged.  Without
      ! this exchange an x-rank seam leaks ~3e-9/day global mass.
      ! D0-unconditional: single-rank non-periodic => no-op, single-rank
      ! periodic => local wrap, multi-rank => messages.  Must run BEFORE the
      ! periodic wrap below so both exchanges see the same post-zonal state.
      !
      ! NOTE: this exchange is also inside the "ocean_continuity" compute region
      ! opened by the dyn caller (rdb_ocean_dyn.F90); ocean_comms_ml here isolates
      ! the comm share, so sums of compute+comms slightly over-close by this term.
      ! This is a known, accepted double-attribution — no stop/restart of the
      ! outer region from inside this module.
      call profiler_start("ocean_comms_ml")
      call ocean_halo_centre(ms%h_layer, nz)
      if (allocated(ms%tracers)) then
         do it = 1, size(ms%tracers)
            if (.not. allocated(ms%tracers(it)%hTr)) cycle
            call ocean_halo_centre(ms%tracers(it)%hTr, nz)
         end do
      end if
      call profiler_stop("ocean_comms_ml")

      ! Mid-Lie-split ghost wrap (design §1.5): re-wrap h_layer (and
      ! tracers) after the zonal apply so the meridional reconstruction
      ! reads current ghost values.  Both per_x and per_y checked: even
      ! for periodic-x only, per_y ghosts can inherit stale values
      ! accumulated under the zonal update.  Two cheap DC kernels, no
      ! correctness traps.
      if (per_x .or. per_y) then
         ! Batched async wrap (queue 1): h_layer + every tracer issued without
         ! per-call sync, then synced ONCE below — pipelines the tiny ghost-slab
         ! launches (otherwise launch-latency-bound).  Wait before the fold,
         ! which reads these wrapped ghosts.
         call ocean_periodic_wrap_centre_3d(ms%h_layer, nx, ny, nz, &
                                            nx_phys, ny_phys, nghost, per_x, per_y, no_wait=.true.)
         if (allocated(ms%tracers)) then
            do it = 1, size(ms%tracers)
               if (.not. allocated(ms%tracers(it)%hTr)) cycle
               call ocean_periodic_wrap_centre_3d(ms%tracers(it)%hTr, nx, ny, nz, &
                                                  nx_phys, ny_phys, nghost, per_x, per_y, no_wait=.true.)
            end do
         end if
         !$acc wait(1)
      end if
      ! Mid-Lie-split north fold (Appendix A): re-fold the centre fields
      ! (h_layer + tracers) AFTER the periodic wrap so the meridional
      ! reconstruction reads fold-consistent north ghosts.  No-op when not
      ! folding.
      if (present(bc)) call ocean_fold_wrap_centre_3d_state(grid, bc, ms)

      if (present(vhbt) .and. present(bc)) then
         call continuity_meridional_flux(grid, metrics, this, ms, dt, vhbt=vhbt, bc=bc, &
                                         visc_rem=visc_rem_v, v_cor=v_cor)
      else if (present(vhbt)) then
         call continuity_meridional_flux(grid, metrics, this, ms, dt, vhbt=vhbt, &
                                         visc_rem=visc_rem_v, v_cor=v_cor)
      else if (present(bc)) then
         call continuity_meridional_flux(grid, metrics, this, ms, dt, bc=bc)
      else
         call continuity_meridional_flux(grid, metrics, this, ms, dt)
      end if
      ! FK MLE fold (B5): add vhml into the meridional mass flux.
      ! Same thermo-cadence / dynamics-fold limitation as the zonal fold above.
      if (present(mle) .and. do_mle_fold) then
         call mle_fold_y(mle, ms%mass_flux_y_layer, nx, ny + 1, nz)
      end if
      if (present(gm) .and. do_mle_fold) then
         call gm_fold_y(gm, ms%mass_flux_y_layer, nx, ny + 1, nz)
      end if
      ! No-normal-flow wall closure for the meridional bolus fold (see the
      ! zonal block above for the rationale).
      if (fold_wall) then
         bc_s_tag = OBC_WALL
         bc_n_tag = OBC_WALL
         if (present(bc)) then
            bc_s_tag = bc%south%bc_type
            bc_n_tag = bc%north%bc_type
         end if
         do concurrent(kk=1:nz, ii=1:nx)
            if (bc_s_tag == OBC_WALL) ms%mass_flux_y_layer(ii, nghost + 1, kk) = 0.0_wp
            if (bc_n_tag == OBC_WALL) ms%mass_flux_y_layer(ii, nghost + ny_phys + 1, kk) = 0.0_wp
         end do
      end if
      ! P2 positive-definite outflux limiter (meridional): mirror of the zonal
      ! pass, on the post-zonal-apply h* availability.  Same D3 single-source
      ! scaling of the folded total.  Off ⇒ skipped ⇒ bit-identical.
      if (this%positive_definite) then
         ! v1.1: forward v_cor — see the zonal twin.
         call pd_limit_meridional_impl(nx, ny, nz, dt, this%h_lim, metrics%iareaT, &
                                       ms%h_layer, ms%mass_flux_y_layer, &
                                       this%pd_theta%data, this%n_limited_step, &
                                       v_cor=v_cor)
      end if
      if (mode == TR_MODE_ACCUMULATE) then
         call accumulate_flux_y(nx, ny + 1, nz, dt, ms%mass_flux_y_layer, this%vhtr)
      else if (mode /= TR_MODE_NONE .and. present(bc)) then
         call tracer_advect_meridional(grid, metrics, this, ms, dt, bc=bc)
      else if (mode /= TR_MODE_NONE) then
         call tracer_advect_meridional(grid, metrics, this, ms, dt)
      end if
      if (h_min_use > 0.0_wp) then
         call continuity_apply_meridional(grid, metrics, ms, dt, h_min=h_min_use)
      else
         call continuity_apply_meridional(grid, metrics, ms, dt)
      end if
      ! Windowed-advect concentration hold (step 2 of 2).  Both applies have
      ! advanced `h_layer`; re-weight the frozen tracer content onto it so the
      ! concentration every downstream consumer reads (`ocean_eos_compute` →
      ! ρ → PGF, vdiff, hdiff, vertical advect, diagnostics) is EXACTLY the
      ! window-start value, as it is in MOM6 (whose prognostic `Tr%t` is a
      ! concentration and is therefore thickness-invariant for free).
      ! Without this, `T` drifts by the full window thickness divergence and
      ! the resulting grid-scale buoyancy error closes an exponentially
      ! growing EOS→PGF→divergence loop.  Composes with `rk2_average`:
      ! `hTr0 = T·h^n` and `hTr = T·h^(2)` average to `T·h^(n+1)`, so `T` is
      ! still exactly `T`.  `continuity_tracer_drain` converts back before it
      ! spends the accumulated transports.
      if (mode == TR_MODE_ACCUMULATE .and. allocated(ms%tracers)) then
         do it_cw = 1, size(ms%tracers)
            if (.not. allocated(ms%tracers(it_cw)%hTr)) cycle
            if (.not. ms%tracers(it_cw)%do_horizontal_advection) cycle
            ! Budget dispatch (mirrors `tracer_advect_zonal`): heat/salt get
            ! the hold's own content change recorded so the console budget
            ! closes at EVERY report, not only on window boundaries.  Tracers
            ! without a budget slot take the budget-free twin.
            select case (ms%tracers(it_cw)%budget_id)
            case (TRACER_BUDGET_HEAT)
               call drain_rescale_hTr_budget(nx, ny, nz, ms%h_layer, this%hprev_work, &
                                             DRAIN_BUDGET_IN_STAGE_WEIGHT, &
                                             ms%tracers(it_cw)%hTr, ms%heat_budget_horiz_adv)
            case (TRACER_BUDGET_SALT)
               call drain_rescale_hTr_budget(nx, ny, nz, ms%h_layer, this%hprev_work, &
                                             DRAIN_BUDGET_IN_STAGE_WEIGHT, &
                                             ms%tracers(it_cw)%hTr, ms%salt_budget_horiz_adv)
            case default
               call drain_rescale_hTr(nx, ny, nz, ms%h_layer, this%hprev_work, &
                                      ms%tracers(it_cw)%hTr)
            end select
         end do
         this%hTr_holds_conc = .true.
      end if

      ! P3: fold this call's limited-face count into the running total (one
      ! add per split call, mirroring dyn%ntrunc_total).  n_limited_step is 0
      ! when positive_definite is off, so this is a no-op there.
      this%n_limited_total = this%n_limited_total + this%n_limited_step
   end subroutine continuity_tracer_step_split

   subroutine pd_limit_zonal_impl(nx, ny, nz, dt, h_lim, iareaT, h_layer, &
                                  mass_flux_x, theta, n_limited, u_cor)
      !! Positive-definite per-donor outflux limiter — zonal (x) pass (P2,
      !! plan §P2 / decision D2).  Scales the per-layer east-face mass fluxes
      !! DOWN so no donor cell loses more thickness than it holds above the
      !! floor `h_lim`: guarantees `h_layer >= h_lim` after
      !! `continuity_apply_zonal`, with ZERO mass created — outfluxes shrink,
      !! thickness is never inflated (the deliberate contrast with MOM6's
      !! `max(h, Angstrom)` injection; the conservative borrow stays the
      !! backstop).  Two device passes:
      !!
      !!   1. per-cell `θ(i,j,k) = min(1, avail/demand)` over the FULL index
      !!      range (incl. ghosts), with `demand = dt·outflow·iareaT` the
      !!      thickness (m) the two x-faces would drain and
      !!      `avail = max(h − h_lim, 0)`.  `outflow` counts only the OUTGOING
      !!      part of each face — east face `i+1` when its flux is positive,
      !!      west face `i` when its flux is negative — because in a Lie pass
      !!      only outgoing faces drain the cell.  θ construction (the
      !!      `merge(avail/max(demand,H_DIV_EPS), 1, demand>avail)` form with
      !!      pure 1/0 armour) mirrors the barotropic wet/dry sweep.
      !!   2. each interior face `i ∈ 2..nx` is multiplied by its UPWIND
      !!      donor's θ (west cell `i−1` when the face flux ≥ 0, else east
      !!      cell `i`) — donor-side only, NO two-cell min: a face's mass
      !!      leaves exactly one donor per direction pass, and that donor's
      !!      OTHER outgoing face is scaled by the SAME θ, so realised total
      !!      outflow ≤ `avail` is already guaranteed by θ_donor alone.
      !!
      !! The θ range spans ghosts so a seam/periodic donor in a ghost column
      !! carries a current factor.  MPI-seam determinism (deferred multi-rank):
      !! two ranks sharing a seam face compute the same θ from the same
      !! exchanged donor `h_layer` and the same face flux, so they scale the
      !! shared face identically.  `n_limited` accumulates the count of scaled
      !! faces (θ_donor < 1, flux ≠ 0) via an `!$acc parallel loop reduction`
      !! (a `do concurrent` + `sum()` would silently read the stale host shadow
      !! on the mem:separate GPU build).
      !!
      !! v1.1 — u_cor re-matching (optional `u_cor`).  `split_scheme="pred_corr"`
      !! captures the MOM6 `u_cor` (the transport-matched velocity the fluxes
      !! correspond to) inside `continuity_zonal_flux` BEFORE this limiter runs;
      !! the `use_state_fluxes` corrector then evaluates at `u_cor` but transports
      !! with the LIMITED `mass_flux_x`.  After a face is scaled by θ its captured
      !! `u_cor` no longer corresponds — `flux = h_face·u·metric` with `h_face`
      !! fixed ⇒ `u` scales linearly with θ, so `u_cor *= th_face` restores the
      !! flux↔velocity match.  Left unfixed this is a flux/velocity
      !! inconsistency (anti-damping class) that surfaces as a late-onset
      !! CFL blow-up on the split-scheme path.  The re-scale runs in
      !! a `do concurrent` guarded by a host `present` flag (the renorm's proven
      !! optional-device-array pattern) BEFORE the mass-flux sweep, so both read
      !! the same UNSCALED-flux donor sign.  Absent ⇒ bit-identical to v1.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: dt, h_lim
      real(wp), intent(in) :: iareaT(nx, ny)
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      real(wp), intent(inout) :: mass_flux_x(nx + 1, ny, nz)
      real(wp), intent(inout) :: theta(nx, ny, nz)
      integer, intent(inout) :: n_limited
      real(wp), intent(inout), optional :: u_cor(nx + 1, ny, nz)
         !! MOM6 `u_cor` capture (transport-matched velocity); when present,
         !! re-scaled by the SAME per-face θ as `mass_flux_x` so it stays
         !! consistent with the limited flux the mom6-scheme corrector reads.

      integer :: i, j, k, n_lim
      real(wp) :: out_e, out_w, outflow, demand, avail, th_face
      logical :: do_ucor

      do_ucor = present(u_cor)

      do concurrent(k=1:nz, j=1:ny, i=1:nx) &
         local(out_e, out_w, outflow, demand, avail)
         out_e = max(mass_flux_x(i + 1, j, k), 0.0_wp)
         out_w = max(-mass_flux_x(i, j, k), 0.0_wp)
         outflow = out_e + out_w
         demand = dt*outflow*iareaT(i, j)
         avail = max(h_layer(i, j, k) - h_lim, 0.0_wp)
         theta(i, j, k) = merge(avail/max(demand, H_DIV_EPS), 1.0_wp, demand > avail)
      end do

      ! v1.1: re-match u_cor to the limited flux, reading the UNSCALED-flux
      ! donor sign (runs before the mass-flux sweep below).  Host-flag-guarded
      ! do concurrent = the renorm's proven optional-device-array pattern.
      if (do_ucor) then
         do concurrent(k=1:nz, j=1:ny, i=2:nx) local(th_face)
            if (mass_flux_x(i, j, k) >= 0.0_wp) then
               th_face = theta(i - 1, j, k)
            else
               th_face = theta(i, j, k)
            end if
            u_cor(i, j, k) = u_cor(i, j, k)*th_face
         end do
      end if

      n_lim = 0
      do concurrent(k=1:nz, j=1:ny, i=2:nx) local(th_face) reduce(+:n_lim)
         if (mass_flux_x(i, j, k) >= 0.0_wp) then
            th_face = theta(i - 1, j, k)
         else
            th_face = theta(i, j, k)
         end if
         if (th_face < 1.0_wp .and. mass_flux_x(i, j, k) /= 0.0_wp) n_lim = n_lim + 1
         mass_flux_x(i, j, k) = mass_flux_x(i, j, k)*th_face
      end do
      n_limited = n_limited + n_lim
   end subroutine pd_limit_zonal_impl

   subroutine pd_limit_meridional_impl(nx, ny, nz, dt, h_lim, iareaT, h_layer, &
                                       mass_flux_y, theta, n_limited, v_cor)
      !! Positive-definite per-donor outflux limiter — meridional (y) pass
      !! (P2).  Mirror of `pd_limit_zonal_impl`; reads the post-zonal-apply
      !! `h_layer` (= h*), which is exactly the availability the second Lie
      !! pass must respect, and scales `mass_flux_y` so `h_layer >= h_lim`
      !! holds after `continuity_apply_meridional`.  Cell `(i,j,k)` outflow =
      !! (north face `j+1` when positive) + (south face `j` when negative);
      !! interior face `j ∈ 2..ny` scaled by its upwind donor's θ (south cell
      !! `j−1` when the face flux ≥ 0, else north cell `j`).  See the zonal
      !! twin for the θ construction, ghost range, MPI-seam determinism, and
      !! the v1.1 `v_cor` re-matching rationale (MOM6 `v_cor` scaled by the same
      !! per-face θ so it stays consistent with the limited flux).
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: dt, h_lim
      real(wp), intent(in) :: iareaT(nx, ny)
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      real(wp), intent(inout) :: mass_flux_y(nx, ny + 1, nz)
      real(wp), intent(inout) :: theta(nx, ny, nz)
      integer, intent(inout) :: n_limited
      real(wp), intent(inout), optional :: v_cor(nx, ny + 1, nz)
         !! MOM6 `v_cor` capture; when present, re-scaled by the SAME per-face θ
         !! as `mass_flux_y` so it stays consistent with the limited flux the
         !! mom6-scheme corrector reads.

      integer :: i, j, k, n_lim
      real(wp) :: out_n, out_s, outflow, demand, avail, th_face
      logical :: do_vcor

      do_vcor = present(v_cor)

      do concurrent(k=1:nz, j=1:ny, i=1:nx) &
         local(out_n, out_s, outflow, demand, avail)
         out_n = max(mass_flux_y(i, j + 1, k), 0.0_wp)
         out_s = max(-mass_flux_y(i, j, k), 0.0_wp)
         outflow = out_n + out_s
         demand = dt*outflow*iareaT(i, j)
         avail = max(h_layer(i, j, k) - h_lim, 0.0_wp)
         theta(i, j, k) = merge(avail/max(demand, H_DIV_EPS), 1.0_wp, demand > avail)
      end do

      ! v1.1: re-match v_cor to the limited flux (UNSCALED-flux donor sign),
      ! before the mass-flux sweep.  See the zonal twin.
      if (do_vcor) then
         do concurrent(k=1:nz, j=2:ny, i=1:nx) local(th_face)
            if (mass_flux_y(i, j, k) >= 0.0_wp) then
               th_face = theta(i, j - 1, k)
            else
               th_face = theta(i, j, k)
            end if
            v_cor(i, j, k) = v_cor(i, j, k)*th_face
         end do
      end if

      n_lim = 0
      do concurrent(k=1:nz, j=2:ny, i=1:nx) local(th_face) reduce(+:n_lim)
         if (mass_flux_y(i, j, k) >= 0.0_wp) then
            th_face = theta(i, j - 1, k)
         else
            th_face = theta(i, j, k)
         end if
         if (th_face < 1.0_wp .and. mass_flux_y(i, j, k) /= 0.0_wp) n_lim = n_lim + 1
         mass_flux_y(i, j, k) = mass_flux_y(i, j, k)*th_face
      end do
      n_limited = n_limited + n_lim
   end subroutine pd_limit_meridional_impl

   pure subroutine accumulate_flux_x(nx_face, ny, nz, dt, mass_flux_x, uhtr)
      !! Accumulate one RK2 stage's zonal mass flux into the windowed
      !! accumulator with the ½ RK2 weight baked in:
      !! `uhtr += 0.5·mass_flux_x·dt`.  `mass_flux_x` is already
      !! area-weighted (m³/s = u·h_face·dy_cu); the product is m³.
      integer, intent(in) :: nx_face, ny, nz
      real(wp), intent(in) :: dt
      real(wp), intent(in) :: mass_flux_x(nx_face, ny, nz)
      real(wp), intent(inout) :: uhtr(nx_face, ny, nz)
      integer :: i, j, k
      do concurrent(k=1:nz, j=1:ny, i=1:nx_face)
         uhtr(i, j, k) = uhtr(i, j, k) + 0.5_wp*mass_flux_x(i, j, k)*dt
      end do
   end subroutine accumulate_flux_x

   pure subroutine accumulate_flux_y(nx, ny_face, nz, dt, mass_flux_y, vhtr)
      !! Meridional analogue of `accumulate_flux_x`.
      integer, intent(in) :: nx, ny_face, nz
      real(wp), intent(in) :: dt
      real(wp), intent(in) :: mass_flux_y(nx, ny_face, nz)
      real(wp), intent(inout) :: vhtr(nx, ny_face, nz)
      integer :: i, j, k
      do concurrent(k=1:nz, j=1:ny_face, i=1:nx)
         vhtr(i, j, k) = vhtr(i, j, k) + 0.5_wp*mass_flux_y(i, j, k)*dt
      end do
   end subroutine accumulate_flux_y

   pure subroutine tracer_advect_zonal_one_impl(nx, ny, nz, dt, iareaT, wet_T, h, hTr, &
                                                mass_flux_x, &
                                                Tr_face_left_x, Tr_face_right_x, &
                                                budget_adv)
      !! Zonal half of `tracer_advect_one_impl`.  Same three-pass
      !! pattern (PPM reconstruction → upwind pick into tracer mass
      !! flux → forward-Euler update) but only the x-direction half.
      !! Reads the input `h` for the Tr = hTr/h reconstruction.  The
      !! tracer transport `mass_flux_x·Tr_face` inherits `dy_cu`, and
      !! the divergence closes with `iareaT` (= `inv_dx` on uniform).
      !! Mirror-T at land neighbours (C2): a held land column's tracer
      !! is reflected to the local cell so the wet-side face value is
      !! unbiased; bit-identical for all-wet (`wet_T≡1`).
      !!
      !! Optional `budget_adv`: when present, a separate `do concurrent`
      !! loop accumulates `−dt·(Tr_face_left_x(i+1)−Tr_face_left_x(i))·iareaT`
      !! into `budget_adv` (same divergence written to `hTr`).  The
      !! prognostic hTr update loop is UNCHANGED — byte-identical when
      !! `budget_adv` is absent.  Used by the console salt/heat closure.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: dt
      real(wp), intent(in) :: iareaT(nx, ny)
      real(wp), intent(in) :: wet_T(nx, ny)
      real(wp), intent(in) :: h(nx, ny, nz)
      real(wp), intent(inout) :: hTr(nx, ny, nz)
      real(wp), intent(in) :: mass_flux_x(nx + 1, ny, nz)
      real(wp), intent(inout) :: Tr_face_left_x(nx + 1, ny, nz)
      real(wp), intent(inout) :: Tr_face_right_x(nx + 1, ny, nz)
      real(wp), intent(inout), optional :: budget_adv(nx, ny, nz)
         !! Per-cell accumulator for the horizontal-advection budget
         !! (same sign/units as hTr).  When present, the zonal-flux
         !! divergence is added (+=) here after the prognostic update.
         !! Absent ⇒ inert (byte-identical to the pre-feature build).

      integer :: i, j, k
      real(wp) :: Tr_m2, Tr_m1, Tr_0, Tr_p1, Tr_p2
      real(wp) :: dh_m1, dh_0, dh_p1, Tr_left, Tr_right
      real(wp) :: mass_x, Tr_face

      do concurrent(k=1:nz, j=1:ny, i=3:nx - 2) &
         local(Tr_m2, Tr_m1, Tr_0, Tr_p1, Tr_p2, &
               dh_m1, dh_0, dh_p1, Tr_left, Tr_right)
         Tr_0 = hTr(i, j, k)/h(i, j, k)
         Tr_m1 = ppm_mirror_h(hTr(i - 1, j, k)/h(i - 1, j, k), Tr_0, wet_T(i - 1, j))
         Tr_p1 = ppm_mirror_h(hTr(i + 1, j, k)/h(i + 1, j, k), Tr_0, wet_T(i + 1, j))
         Tr_m2 = ppm_mirror_h(hTr(i - 2, j, k)/h(i - 2, j, k), Tr_m1, wet_T(i - 2, j))
         Tr_p2 = ppm_mirror_h(hTr(i + 2, j, k)/h(i + 2, j, k), Tr_p1, wet_T(i + 2, j))
         call ppm_limited_slope(Tr_m2, Tr_m1, Tr_0, dh_m1)
         call ppm_limited_slope(Tr_m1, Tr_0, Tr_p1, dh_0)
         call ppm_limited_slope(Tr_0, Tr_p1, Tr_p2, dh_p1)
         dh_0 = dh_0*wet_T(i - 1, j)*wet_T(i, j)*wet_T(i + 1, j)
         Tr_left = 0.5_wp*(Tr_m1 + Tr_0) - (dh_0 - dh_m1)/6.0_wp
         Tr_right = 0.5_wp*(Tr_0 + Tr_p1) - (dh_p1 - dh_0)/6.0_wp
         call ppm_cell_limiter(Tr_0, Tr_left, Tr_right)
         Tr_face_right_x(i, j, k) = Tr_left
         Tr_face_left_x(i + 1, j, k) = Tr_right
      end do
      do concurrent(k=1:nz, j=1:ny)
         Tr_face_left_x(1, j, k) = hTr(1, j, k)/h(1, j, k)
         Tr_face_right_x(1, j, k) = hTr(1, j, k)/h(1, j, k)
         Tr_face_left_x(2, j, k) = hTr(1, j, k)/h(1, j, k)
         Tr_face_right_x(2, j, k) = hTr(2, j, k)/h(2, j, k)
         Tr_face_left_x(3, j, k) = hTr(2, j, k)/h(2, j, k)
         Tr_face_right_x(3, j, k) = hTr(2, j, k)/h(2, j, k)
         Tr_face_right_x(nx - 1, j, k) = hTr(nx - 1, j, k)/h(nx - 1, j, k)
         Tr_face_left_x(nx, j, k) = hTr(nx - 1, j, k)/h(nx - 1, j, k)
         Tr_face_right_x(nx, j, k) = hTr(nx, j, k)/h(nx, j, k)
         Tr_face_left_x(nx + 1, j, k) = hTr(nx, j, k)/h(nx, j, k)
         Tr_face_right_x(nx + 1, j, k) = hTr(nx, j, k)/h(nx, j, k)
      end do

      do concurrent(k=1:nz, j=1:ny, i=1:nx + 1) local(mass_x, Tr_face)
         mass_x = mass_flux_x(i, j, k)
         if (mass_x >= 0.0_wp) then
            Tr_face = Tr_face_left_x(i, j, k)
         else
            Tr_face = Tr_face_right_x(i, j, k)
         end if
         Tr_face_left_x(i, j, k) = mass_x*Tr_face
      end do

      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         hTr(i, j, k) = hTr(i, j, k) - dt* &
                        (Tr_face_left_x(i + 1, j, k) - Tr_face_left_x(i, j, k))*iareaT(i, j)
      end do

      ! Budget accumulator: separate guarded loop so the prognostic update above
      ! is UNCHANGED (byte-identical when budget_adv is absent).  Mirrors the
      ! same divergence that was just applied to hTr — the console uses it to
      ! close the salt/heat conservation residual.
      if (present(budget_adv)) then
         do concurrent(k=1:nz, j=1:ny, i=1:nx)
            budget_adv(i, j, k) = budget_adv(i, j, k) - dt* &
                                  (Tr_face_left_x(i + 1, j, k) - Tr_face_left_x(i, j, k))*iareaT(i, j)
         end do
      end if
   end subroutine tracer_advect_zonal_one_impl

   pure subroutine tracer_advect_meridional_one_impl(nx, ny, nz, dt, iareaT, wet_T, h, hTr, &
                                                     mass_flux_y, &
                                                     Tr_face_left_y, Tr_face_right_y, &
                                                     budget_adv)
      !! Meridional half of `tracer_advect_one_impl`.  Same shape
      !! as the zonal impl, applied to y.  In the split flow, `h`
      !! here is the post-zonal-apply thickness so the Tr = hTr/h
      !! reconstruction stays consistent with what continuity used
      !! in `continuity_meridional_flux`.  `iareaT` = `inv_dy` on
      !! uniform metrics; `mass_flux_y` carries `dx_cv`.  Mirror-T at
      !! land neighbours (C2); bit-identical for all-wet.
      !!
      !! Optional `budget_adv`: accumulates `−dt·(Tr_face_left_y(i,j+1)
      !! −Tr_face_left_y(i,j))·iareaT` into `budget_adv`.  Absent ⇒
      !! inert.  Called in sequence after the zonal impl so the two
      !! directions compose via `+=` into the same accumulator array.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: dt
      real(wp), intent(in) :: iareaT(nx, ny)
      real(wp), intent(in) :: wet_T(nx, ny)
      real(wp), intent(in) :: h(nx, ny, nz)
      real(wp), intent(inout) :: hTr(nx, ny, nz)
      real(wp), intent(in) :: mass_flux_y(nx, ny + 1, nz)
      real(wp), intent(inout) :: Tr_face_left_y(nx, ny + 1, nz)
      real(wp), intent(inout) :: Tr_face_right_y(nx, ny + 1, nz)
      real(wp), intent(inout), optional :: budget_adv(nx, ny, nz)
         !! Per-cell meridional-advection budget accumulator (PSU·m or
         !! °C·m per cell).  Added to the same array as the zonal half
         !! so the net entry covers both directions.  Absent ⇒ inert.

      integer :: i, j, k
      real(wp) :: Tr_m2, Tr_m1, Tr_0, Tr_p1, Tr_p2
      real(wp) :: dh_m1, dh_0, dh_p1, Tr_left, Tr_right
      real(wp) :: mass_y, Tr_face

      do concurrent(k=1:nz, j=3:ny - 2, i=1:nx) &
         local(Tr_m2, Tr_m1, Tr_0, Tr_p1, Tr_p2, &
               dh_m1, dh_0, dh_p1, Tr_left, Tr_right)
         Tr_0 = hTr(i, j, k)/h(i, j, k)
         Tr_m1 = ppm_mirror_h(hTr(i, j - 1, k)/h(i, j - 1, k), Tr_0, wet_T(i, j - 1))
         Tr_p1 = ppm_mirror_h(hTr(i, j + 1, k)/h(i, j + 1, k), Tr_0, wet_T(i, j + 1))
         Tr_m2 = ppm_mirror_h(hTr(i, j - 2, k)/h(i, j - 2, k), Tr_m1, wet_T(i, j - 2))
         Tr_p2 = ppm_mirror_h(hTr(i, j + 2, k)/h(i, j + 2, k), Tr_p1, wet_T(i, j + 2))
         call ppm_limited_slope(Tr_m2, Tr_m1, Tr_0, dh_m1)
         call ppm_limited_slope(Tr_m1, Tr_0, Tr_p1, dh_0)
         call ppm_limited_slope(Tr_0, Tr_p1, Tr_p2, dh_p1)
         dh_0 = dh_0*wet_T(i, j - 1)*wet_T(i, j)*wet_T(i, j + 1)
         Tr_left = 0.5_wp*(Tr_m1 + Tr_0) - (dh_0 - dh_m1)/6.0_wp
         Tr_right = 0.5_wp*(Tr_0 + Tr_p1) - (dh_p1 - dh_0)/6.0_wp
         call ppm_cell_limiter(Tr_0, Tr_left, Tr_right)
         Tr_face_right_y(i, j, k) = Tr_left
         Tr_face_left_y(i, j + 1, k) = Tr_right
      end do
      do concurrent(k=1:nz, i=1:nx)
         Tr_face_left_y(i, 1, k) = hTr(i, 1, k)/h(i, 1, k)
         Tr_face_right_y(i, 1, k) = hTr(i, 1, k)/h(i, 1, k)
         Tr_face_left_y(i, 2, k) = hTr(i, 1, k)/h(i, 1, k)
         Tr_face_right_y(i, 2, k) = hTr(i, 2, k)/h(i, 2, k)
         Tr_face_left_y(i, 3, k) = hTr(i, 2, k)/h(i, 2, k)
         Tr_face_right_y(i, 3, k) = hTr(i, 2, k)/h(i, 2, k)
         Tr_face_right_y(i, ny - 1, k) = hTr(i, ny - 1, k)/h(i, ny - 1, k)
         Tr_face_left_y(i, ny, k) = hTr(i, ny - 1, k)/h(i, ny - 1, k)
         Tr_face_right_y(i, ny, k) = hTr(i, ny, k)/h(i, ny, k)
         Tr_face_left_y(i, ny + 1, k) = hTr(i, ny, k)/h(i, ny, k)
         Tr_face_right_y(i, ny + 1, k) = hTr(i, ny, k)/h(i, ny, k)
      end do

      do concurrent(k=1:nz, j=1:ny + 1, i=1:nx) local(mass_y, Tr_face)
         mass_y = mass_flux_y(i, j, k)
         if (mass_y >= 0.0_wp) then
            Tr_face = Tr_face_left_y(i, j, k)
         else
            Tr_face = Tr_face_right_y(i, j, k)
         end if
         Tr_face_left_y(i, j, k) = mass_y*Tr_face
      end do

      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         hTr(i, j, k) = hTr(i, j, k) - dt* &
                        (Tr_face_left_y(i, j + 1, k) - Tr_face_left_y(i, j, k))*iareaT(i, j)
      end do

      ! Budget accumulator: same guarded pattern as the zonal half.
      if (present(budget_adv)) then
         do concurrent(k=1:nz, j=1:ny, i=1:nx)
            budget_adv(i, j, k) = budget_adv(i, j, k) - dt* &
                                  (Tr_face_left_y(i, j + 1, k) - Tr_face_left_y(i, j, k))*iareaT(i, j)
         end do
      end if
   end subroutine tracer_advect_meridional_one_impl

   subroutine tracer_advect(grid, metrics, this, ms, dt)
      !! **Test-only** (no production caller): unsplit 2-D tracer advection,
      !! the reference oracle for split `tracer_advect_zonal`/`_meridional`.
      !! Per-layer PPM tracer advection — iterates over the tracer
      !! registry on the multilayer C-grid state and forwards each
      !! tracer's hTr array to the flat-impl below.  Extendable by
      !! construction: appending a new entry to `ms%tracers(:)` (BGC,
      !! sediment, passive scalar) drops it in without touching this
      !! routine.  Per-tracer behaviour gates on
      !! `tracer_t%do_horizontal_advection` — set to .false. for
      !! tracers that should be diagnostic / forced externally.
      !!
      !! Consistency-with-continuity (CWC): the kernel consumes the
      !! same `mass_flux_*_layer` that continuity-PPM produced, so a
      !! uniform tracer remains uniform under non-zero flow (the
      !! standard discrete CWC theorem).  Verified by
      !! `tracer_cwc_uniform` in `test_ocean_tracer_adv`.
      !!
      !! Scratch is borrowed from `continuity_t` — those four
      !! face-reconstruction buffers are unused at this point in the
      !! step (continuity_compute_fluxes already wrote and consumed
      !! them in the upwind pick), so we re-use the same allocations
      !! for tracer face values + tracer mass fluxes without
      !! introducing a new slot.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(continuity_t), intent(inout) :: this
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: dt

      integer :: it

      if (.not. allocated(ms%tracers)) return

      do it = 1, size(ms%tracers)
         if (.not. ms%tracers(it)%do_horizontal_advection) cycle
         call tracer_advect_one_impl( &
            grid%nx_total, grid%ny_total, ms%nz_ml, &
            dt, metrics%iareaT, metrics%wet_T, &
            ms%h_layer, &
            ms%tracers(it)%hTr, &
            ms%mass_flux_x_layer, &
            ms%mass_flux_y_layer, &
            this%h_face_left_x%data, &
            this%h_face_right_x%data, &
            this%h_face_left_y%data, &
            this%h_face_right_y%data)
      end do
   end subroutine tracer_advect

   pure subroutine tracer_advect_one_impl(nx, ny, nz, dt, iareaT, wet_T, h, hTr, &
                                          mass_flux_x, mass_flux_y, &
                                          Tr_face_left_x, Tr_face_right_x, &
                                          Tr_face_left_y, Tr_face_right_y)
      !! One-tracer PPM advection.  Flat-impl: takes bare 3D arrays
      !! (no derived-type deref inside do-concurrent), so NVHPC
      !! stdpar handles it cleanly even for tracers stored in an
      !! array-of-derived-types registry.
      !!
      !! Three passes per direction:
      !!   1. PPM reconstruction on Tr = hTr/h (Colella-Woodward
      !!      monotonic limiter) — same helpers as continuity-PPM.
      !!      Stores left/right face states in the scratch buffers.
      !!   2. Upwind pick using mass_flux sign — overwrites the same
      !!      scratch buffer with mass_flux * Tr_face (the tracer
      !!      mass flux at each face).
      !!   3. Flux divergence + forward-Euler update of hTr.
      !!
      !! Wall faces are gated by mass_flux (= 0 at walls from
      !! continuity), so tracer flux through walls is identically
      !! zero — no special wall handling needed in this kernel.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: dt
      real(wp), intent(in) :: iareaT(nx, ny)
      real(wp), intent(in) :: wet_T(nx, ny)
      real(wp), intent(in) :: h(nx, ny, nz)
      real(wp), intent(inout) :: hTr(nx, ny, nz)
      real(wp), intent(in) :: mass_flux_x(nx + 1, ny, nz)
      real(wp), intent(in) :: mass_flux_y(nx, ny + 1, nz)
      real(wp), intent(inout) :: Tr_face_left_x(nx + 1, ny, nz)
      real(wp), intent(inout) :: Tr_face_right_x(nx + 1, ny, nz)
      real(wp), intent(inout) :: Tr_face_left_y(nx, ny + 1, nz)
      real(wp), intent(inout) :: Tr_face_right_y(nx, ny + 1, nz)

      integer :: i, j, k
      real(wp) :: Tr_m2, Tr_m1, Tr_0, Tr_p1, Tr_p2
      real(wp) :: dh_m1, dh_0, dh_p1, Tr_left, Tr_right
      real(wp) :: mass_x, mass_y, Tr_face

      ! ===== X-direction PPM reconstruction =====
      ! Interior cells with full 5-point stencil.  Mirror-T at land
      ! neighbours (C2); bit-identical for all-wet.
      do concurrent(k=1:nz, j=1:ny, i=3:nx - 2) &
         local(Tr_m2, Tr_m1, Tr_0, Tr_p1, Tr_p2, &
               dh_m1, dh_0, dh_p1, Tr_left, Tr_right)
         Tr_0 = hTr(i, j, k)/h(i, j, k)
         Tr_m1 = ppm_mirror_h(hTr(i - 1, j, k)/h(i - 1, j, k), Tr_0, wet_T(i - 1, j))
         Tr_p1 = ppm_mirror_h(hTr(i + 1, j, k)/h(i + 1, j, k), Tr_0, wet_T(i + 1, j))
         Tr_m2 = ppm_mirror_h(hTr(i - 2, j, k)/h(i - 2, j, k), Tr_m1, wet_T(i - 2, j))
         Tr_p2 = ppm_mirror_h(hTr(i + 2, j, k)/h(i + 2, j, k), Tr_p1, wet_T(i + 2, j))
         call ppm_limited_slope(Tr_m2, Tr_m1, Tr_0, dh_m1)
         call ppm_limited_slope(Tr_m1, Tr_0, Tr_p1, dh_0)
         call ppm_limited_slope(Tr_0, Tr_p1, Tr_p2, dh_p1)
         dh_0 = dh_0*wet_T(i - 1, j)*wet_T(i, j)*wet_T(i + 1, j)
         Tr_left = 0.5_wp*(Tr_m1 + Tr_0) - (dh_0 - dh_m1)/6.0_wp
         Tr_right = 0.5_wp*(Tr_0 + Tr_p1) - (dh_p1 - dh_0)/6.0_wp
         call ppm_cell_limiter(Tr_0, Tr_left, Tr_right)
         Tr_face_right_x(i, j, k) = Tr_left
         Tr_face_left_x(i + 1, j, k) = Tr_right
      end do
      ! Boundary cells (i=1, 2, nx-1, nx): 1st-order fallback (face
      ! value = abutting cell's concentration).  Mass flux is zero
      ! at i=1 and i=nx+1 walls anyway, so the wall faces' Tr value
      ! doesn't enter the divergence.
      do concurrent(k=1:nz, j=1:ny)
         Tr_face_left_x(1, j, k) = hTr(1, j, k)/h(1, j, k)
         Tr_face_right_x(1, j, k) = hTr(1, j, k)/h(1, j, k)
         Tr_face_left_x(2, j, k) = hTr(1, j, k)/h(1, j, k)
         Tr_face_right_x(2, j, k) = hTr(2, j, k)/h(2, j, k)
         ! Face 3: cell 2's 5-point stencil is short, so its right
         ! edge falls back to 1st order — match continuity's fix.
         Tr_face_left_x(3, j, k) = hTr(2, j, k)/h(2, j, k)
         Tr_face_right_x(3, j, k) = hTr(2, j, k)/h(2, j, k)
         ! Face nx-1: mirror of face 3.
         Tr_face_right_x(nx - 1, j, k) = hTr(nx - 1, j, k)/h(nx - 1, j, k)
         Tr_face_left_x(nx, j, k) = hTr(nx - 1, j, k)/h(nx - 1, j, k)
         Tr_face_right_x(nx, j, k) = hTr(nx, j, k)/h(nx, j, k)
         Tr_face_left_x(nx + 1, j, k) = hTr(nx, j, k)/h(nx, j, k)
         Tr_face_right_x(nx + 1, j, k) = hTr(nx, j, k)/h(nx, j, k)
      end do

      ! Upwind pick at east faces -> overwrite Tr_face_left_x with
      ! the per-face tracer mass flux (mass_flux_x * Tr_face_upwind).
      do concurrent(k=1:nz, j=1:ny, i=1:nx + 1) local(mass_x, Tr_face)
         mass_x = mass_flux_x(i, j, k)
         if (mass_x >= 0.0_wp) then
            Tr_face = Tr_face_left_x(i, j, k)
         else
            Tr_face = Tr_face_right_x(i, j, k)
         end if
         Tr_face_left_x(i, j, k) = mass_x*Tr_face
      end do

      ! ===== Y-direction PPM reconstruction =====
      do concurrent(k=1:nz, j=3:ny - 2, i=1:nx) &
         local(Tr_m2, Tr_m1, Tr_0, Tr_p1, Tr_p2, &
               dh_m1, dh_0, dh_p1, Tr_left, Tr_right)
         Tr_0 = hTr(i, j, k)/h(i, j, k)
         Tr_m1 = ppm_mirror_h(hTr(i, j - 1, k)/h(i, j - 1, k), Tr_0, wet_T(i, j - 1))
         Tr_p1 = ppm_mirror_h(hTr(i, j + 1, k)/h(i, j + 1, k), Tr_0, wet_T(i, j + 1))
         Tr_m2 = ppm_mirror_h(hTr(i, j - 2, k)/h(i, j - 2, k), Tr_m1, wet_T(i, j - 2))
         Tr_p2 = ppm_mirror_h(hTr(i, j + 2, k)/h(i, j + 2, k), Tr_p1, wet_T(i, j + 2))
         call ppm_limited_slope(Tr_m2, Tr_m1, Tr_0, dh_m1)
         call ppm_limited_slope(Tr_m1, Tr_0, Tr_p1, dh_0)
         call ppm_limited_slope(Tr_0, Tr_p1, Tr_p2, dh_p1)
         dh_0 = dh_0*wet_T(i, j - 1)*wet_T(i, j)*wet_T(i, j + 1)
         Tr_left = 0.5_wp*(Tr_m1 + Tr_0) - (dh_0 - dh_m1)/6.0_wp
         Tr_right = 0.5_wp*(Tr_0 + Tr_p1) - (dh_p1 - dh_0)/6.0_wp
         call ppm_cell_limiter(Tr_0, Tr_left, Tr_right)
         Tr_face_right_y(i, j, k) = Tr_left
         Tr_face_left_y(i, j + 1, k) = Tr_right
      end do
      do concurrent(k=1:nz, i=1:nx)
         Tr_face_left_y(i, 1, k) = hTr(i, 1, k)/h(i, 1, k)
         Tr_face_right_y(i, 1, k) = hTr(i, 1, k)/h(i, 1, k)
         Tr_face_left_y(i, 2, k) = hTr(i, 1, k)/h(i, 1, k)
         Tr_face_right_y(i, 2, k) = hTr(i, 2, k)/h(i, 2, k)
         Tr_face_left_y(i, 3, k) = hTr(i, 2, k)/h(i, 2, k)
         Tr_face_right_y(i, 3, k) = hTr(i, 2, k)/h(i, 2, k)
         Tr_face_right_y(i, ny - 1, k) = hTr(i, ny - 1, k)/h(i, ny - 1, k)
         Tr_face_left_y(i, ny, k) = hTr(i, ny - 1, k)/h(i, ny - 1, k)
         Tr_face_right_y(i, ny, k) = hTr(i, ny, k)/h(i, ny, k)
         Tr_face_left_y(i, ny + 1, k) = hTr(i, ny, k)/h(i, ny, k)
         Tr_face_right_y(i, ny + 1, k) = hTr(i, ny, k)/h(i, ny, k)
      end do

      ! Upwind pick at north faces
      do concurrent(k=1:nz, j=1:ny + 1, i=1:nx) local(mass_y, Tr_face)
         mass_y = mass_flux_y(i, j, k)
         if (mass_y >= 0.0_wp) then
            Tr_face = Tr_face_left_y(i, j, k)
         else
            Tr_face = Tr_face_right_y(i, j, k)
         end if
         Tr_face_left_y(i, j, k) = mass_y*Tr_face
      end do

      ! ===== Apply: forward-Euler hTr update (transport div · iareaT) =====
      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         hTr(i, j, k) = hTr(i, j, k) - dt*( &
                        (Tr_face_left_x(i + 1, j, k) - Tr_face_left_x(i, j, k)) + &
                        (Tr_face_left_y(i, j + 1, k) - Tr_face_left_y(i, j, k)))*iareaT(i, j)
      end do
   end subroutine tracer_advect_one_impl

   pure elemental function ppm_mirror_h(h_nbr, h_loc, w_nbr) result(h_out)
      !$acc routine seq
      !! Mirror-h at a land neighbour (spec §14 C2 / MOM6's
      !! reflected-coast PPM): substitute the LOCAL cell's
      !! thickness (or tracer) for a LAND neighbour's held floor value
      !! so the PPM parabola sees a flat, reflected coast and the
      !! wet-side face value is not biased by the dry column.  Branchless:
      !!   `h_out = w_nbr·h_nbr + (1-w_nbr)·h_loc`.
      !! Wet neighbour (`w_nbr=1`) ⇒ `h_out = h_nbr` (literal no-op);
      !! land neighbour (`w_nbr=0`) ⇒ `h_out = h_loc`.
      !!
      !! Consumer: `rdb_ice_transport` (sea-ice PR 4b) reuses this
      !! verbatim for the category-summed ice/snow PPM parabola —
      !! SIS2's mask2dT substitution (`SIS_continuity.F90:1478-1479`).
      real(wp), intent(in) :: h_nbr  !! neighbour value
      real(wp), intent(in) :: h_loc  !! local-cell value (the mirror target)
      real(wp), intent(in) :: w_nbr  !! neighbour wet mask (0/1)
      real(wp) :: h_out
      h_out = w_nbr*h_nbr + (1.0_wp - w_nbr)*h_loc
   end function ppm_mirror_h

   pure elemental function volcfl_face(h_edge, dh, curv3, cfl) result(h_face)
      !$acc routine seq
      !! MOM6 swept-volume continuity-PPM face thickness (Lin & Rood
      !! / MOM_continuity_PPM `flux_elem`).  Integrates the donor
      !! cell's reconstructed parabola over the swept volume rather
      !! than sampling the edge value, adding the O(CFL) correction:
      !!
      !!   `h_face = h_edge + CFL·(0.5·dh + curv3·(CFL − 1.5))`
      !!
      !! `h_edge` is the donor's downwind-facing PPM edge value (the
      !! CFL→0 limit), `dh` is the donor's left-minus-right edge
      !! difference *in the swept (upwind→downwind) orientation*
      !! (MOM6 `h_L − h_R` for u>0, `h_R_p1 − h_L_p1` for u<0), and
      !! `curv3 = (h_L + h_R) − 2·h_centre` is the parabola curvature
      !! measure.  `CFL = |u|·dt·dy_Cu·IareaT` of the DONOR cell.
      !!
      !! At `CFL = 0` this returns exactly `h_edge` — i.e. the current
      !! CFL-free edge-value pick — so `vol_cfl = .false.` (which never
      !! calls this) and the `CFL = 0` limit are bit-identical.
      !!
      !! Consumer: `rdb_ice_transport` (sea-ice PR 4b) reuses this for the
      !! category-summed total-mass face transport — bit-for-bit SIS2's
      !! `zonal_mass_flux`/`meridional_mass_flux` face integral
      !! (`SIS_continuity.F90:1168`).
      real(wp), intent(in) :: h_edge  !! downwind-facing donor PPM edge value
      real(wp), intent(in) :: dh      !! swept-oriented donor edge difference
      real(wp), intent(in) :: curv3   !! donor parabola curvature (h_L+h_R−2h)
      real(wp), intent(in) :: cfl     !! donor-cell Courant number (>= 0)
      real(wp) :: h_face
      h_face = h_edge + cfl*(0.5_wp*dh + curv3*(cfl - 1.5_wp))
   end function volcfl_face

   pure subroutine ppm_limited_slope(h_im1, h_i, h_ip1, dh)
      !$acc routine seq
      !! Van Leer monotonized centred slope for cell i.  Returns 0
      !! at local extrema (sign change between left and right
      !! differences) and the slope-limited centred derivative
      !! otherwise.  Standard PPM convention; see Colella-Woodward
      !! 1984.
      !!
      !! Consumer: `rdb_ice_transport` (sea-ice PR 4b), same role as here
      !! (SIS2 Lin-94 slope + limit, `SIS_continuity.F90:1464-1468`).
      real(wp), intent(in) :: h_im1, h_i, h_ip1
      real(wp), intent(out) :: dh
      real(wp) :: dh_centered, dh_left, dh_right

      dh_left = h_i - h_im1
      dh_right = h_ip1 - h_i
      dh_centered = 0.5_wp*(dh_left + dh_right)
      if (dh_left*dh_right > 0.0_wp) then
         dh = sign(min(abs(dh_centered), &
                       2.0_wp*abs(dh_left), &
                       2.0_wp*abs(dh_right)), &
                   dh_centered)
      else
         dh = 0.0_wp
      end if
   end subroutine ppm_limited_slope

   pure subroutine ppm_cell_limiter(h_centre, h_left, h_right)
      !$acc routine seq
      !! Colella-Woodward 1984 eq 1.10 monotonic limiter on the
      !! parabolic profile in a single cell.  Three branches:
      !!
      !!   1. Local extremum in cell (h_centre lies outside
      !!      [min(h_left,h_right), max(h_left,h_right)]): flatten
      !!      the parabola — h_left = h_right = h_centre.
      !!   2. "Overshoot" at the left edge — reset h_left so the
      !!      parabola's minimum/maximum lies at the right edge.
      !!   3. "Overshoot" at the right edge — symmetric.
      !!
      !! Consumer: `rdb_ice_transport` (sea-ice PR 4b) — SIS2's
      !! `PPM_limit_CW84`.
      real(wp), intent(in) :: h_centre
      real(wp), intent(inout) :: h_left, h_right
      real(wp) :: dh_lr, h_six

      dh_lr = h_right - h_left
      h_six = 6.0_wp*(h_centre - 0.5_wp*(h_left + h_right))
      if ((h_right - h_centre)*(h_centre - h_left) <= 0.0_wp) then
         h_left = h_centre
         h_right = h_centre
      else if (dh_lr*h_six > dh_lr*dh_lr) then
         h_left = 3.0_wp*h_centre - 2.0_wp*h_right
      else if (dh_lr*h_six < -dh_lr*dh_lr) then
         h_right = 3.0_wp*h_centre - 2.0_wp*h_left
      end if
   end subroutine ppm_cell_limiter

   pure subroutine ppm_limit_pos(h_centre, h_left, h_right, h_min)
      !$acc routine seq
      !! Positivity-preserving limiter on the PPM reconstruction.
      !! Mirrors MOM6's `PPM_limit_pos`:
      !! when the parabolic fit predicts a minimum interior to the
      !! cell that dips below `h_min`, shrink h_left / h_right toward
      !! h_centre so the minimum sits at exactly `h_min`.  Pure
      !! scalar form per cell; runs after `ppm_cell_limiter` so the
      !! monotonic-limited reconstruction is the input.
      !!
      !! Algorithm:
      !!   curv = 3·(h_L + h_R − 2·h_in)         ! +ve ⇒ interior min
      !!   if curv > 0 and |dh| < curv:           ! min inside cell
      !!     if h_in ≤ h_min: flatten (h_L = h_R = h_in)
      !!     elif 12·curv·(h_in − h_min) < curv² + 3·dh²:
      !!        loc_scale = 12·curv·(h_in − h_min) / (curv² + 3·dh²) ∈ (0,1)
      !!        h_L = h_in + loc_scale·(h_L − h_in)
      !!        h_R = h_in + loc_scale·(h_R − h_in)
      !!
      !! `h_min = 0` ⇒ pure positivity (parabola can't go negative
      !! inside the cell).  Larger `h_min` ⇒ harder floor; matches
      !! MOM6's `GV%Angstrom_H` for the reduced-gravity setup.
      !!
      !! No-op when curv ≤ 0 (maximum interior, or linear / monotone
      !! profile) or |dh| ≥ curv (minimum outside the cell, edges
      !! already control).
      !!
      !! Consumer: `rdb_ice_transport` (sea-ice PR 4b) — SIS2 runs this
      !! unconditionally on the PD continuity scheme
      !! (`SIS_continuity.F90:1585`, `PPM_limit_pos`), so the ice-mass
      !! reconstruction always applies it (not optional there, unlike
      !! the ocean's `use_ppm_limit_pos` knob).
      real(wp), intent(in) :: h_centre, h_min
      real(wp), intent(inout) :: h_left, h_right
      real(wp) :: curv, dh, loc_scale

      curv = 3.0_wp*((h_left + h_right) - 2.0_wp*h_centre)
      if (curv > 0.0_wp) then
         dh = h_right - h_left
         if (abs(dh) < curv) then
            if (h_centre <= h_min) then
               h_left = h_centre
               h_right = h_centre
            else if (12.0_wp*curv*(h_centre - h_min) < (curv*curv + 3.0_wp*dh*dh)) then
               loc_scale = 12.0_wp*curv*(h_centre - h_min)/(curv*curv + 3.0_wp*dh*dh)
               h_left = h_centre + loc_scale*(h_left - h_centre)
               h_right = h_centre + loc_scale*(h_right - h_centre)
            end if
         end if
      end if
   end subroutine ppm_limit_pos

   subroutine continuity_tracer_drain(grid, metrics, this, ms, ratio, bc)
      !! Phase-2 (6b) windowed horizontal tracer-advection drain.
      !!
      !! Spends the accumulated face transports `this%uhtr/vhtr` (built
      !! over `ratio` outer steps in `TR_MODE_ACCUMULATE`) onto the
      !! frozen tracer mass `ms%tracers(:)%hTr` using a fixed-budget
      !! Colella-Woodward swept-average PPM sub-cycle (Adcroft &
      !! Hallberg 2006; MOM6 ADVECT_PPM).  Conservation is structural:
      !!
      !!   hprev = areaT·h_end + div(uhtr)        (≈ window-start thickness)
      !!   Tr_start = hTr_frozen / hprev          (window-start concentration)
      !!
      !! the drain evolves hprev → h_end while moving tracer with each
      !! limited transport portion; at the end hprev == h_end so
      !! Σ(hTr) = Σ(Tr_start·hprev) = Σ(hTr_frozen) to round-off.  Do NOT
      !! seed the concentration from hTr/h_end (that is the drifted value).
      !!
      !! Per-pass re-reconstruction (V2): the PPM parabola is rebuilt from
      !! the CURRENT Tr at the start of every sub-cycle pass.  The
      !! `hup/hlos/min_h` two-test limiter sets the drained transport
      !! `uhh` per pass; the swept-average parabola sets the concentration
      !! multiplying it — the two are orthogonal.
      !!
      !! Single-rank scope: periodic-x / north-fold seams are handled via
      !! the existing ocean periodic/fold ghost wraps.  Multi-rank C-grid
      !! halos ride E1; gated behind the TODO(E1) note below.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(continuity_t), intent(inout) :: this
      type(multilayer_state_t), intent(inout) :: ms
      integer, intent(in) :: ratio
         !! `dt_tracer_advect_ratio` for this window; sets the
         !! fixed-budget pass count `max_iter = 2·ratio+1` (ratio is an
         !! integer ⇒ ceil(ratio) = ratio).
      type(ocean_bc_state_t), intent(in), optional :: bc

      integer :: nx, ny, nz, nx_phys, ny_phys, nghost
      integer :: max_iter, ipass, it
      integer :: i, j, k
      real(wp) :: cfl_max, denom, fmax
      logical :: per_x, per_y, fold_n

      if (.not. allocated(ms%tracers)) return

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml
      nx_phys = grid%nx_phys
      ny_phys = grid%ny_phys
      nghost = grid%nghost

      per_x = .false.
      per_y = .false.
      fold_n = .false.
      if (present(bc)) then
         per_x = bc%periodic_x .and. .not. ocean_halo_is_decomposed_x()
         per_y = bc%periodic_y .and. .not. ocean_halo_is_decomposed_y()
         fold_n = bc%north_fold
      end if

      ! ---- Halo of the accumulators (single-rank seams) ----
      ! hprev reads neighbour uhtr/vhtr across the ghost, so wrap them
      ! first.  TODO(E1): the multi-rank C-grid face halo of uhtr/vhtr
      ! rides the MPI halo work; single-rank periodic/fold is in scope.
      call drain_wrap_face_x(this%uhtr, nx, ny, nz, nx_phys, ny_phys, nghost, &
                             per_x, per_y, fold_n)
      call drain_wrap_face_y(this%vhtr, nx, ny, nz, nx_phys, ny_phys, nghost, &
                             per_x, per_y, fold_n)

      ! ---- Conservative availability limiter (Fox-Kemper × windowed-advect) ----
      ! Tighten the accumulated window transports so every cell's
      ! reconstructed window-start volume stays ≥ areaT·h_min — the
      ! combined FK + resolved transport can otherwise overdraw a thin
      ! z* surface layer (vol < 0), tripping the non-conservative
      ! clamp/hatch below and leaking ~5%/day of tracer.  Limited uhtr/vhtr
      ! feed BOTH the reconstruction AND the sub-cycle seed (drain_copy_3d
      ! below), so the drain still telescopes exactly onto h_end and
      ! Σ(areaT·hTr) is conserved.  Inert (bit-identical) when no cell
      ! would otherwise go vol ≤ areaT·h_min (ratio=1 / FK off / benign
      ! windows).  tr_work is reused as the per-cell scale scratch — it is
      ! overwritten by drain_parabola_* in the sub-cycle before any read.
      ! Fixed-budget FCT tighten (GPU-uniform, no data-dependent exit);
      ! AVAIL_LIMIT_PASS bounds the constraint diffusion across the stencil.
      do ipass = 1, AVAIL_LIMIT_PASS
         call drain_avail_limit(nx, ny, nz, metrics%areaT, this%h_min, &
                                ms%h_layer, this%uhtr, this%vhtr, this%tr_work)
         call drain_wrap_centre(this%tr_work, nx, ny, nz, nx_phys, ny_phys, &
                                nghost, per_x, per_y, fold_n)
         call drain_avail_scale_x(nx, ny, nz, this%tr_work, this%uhtr)
         call drain_avail_scale_y(nx, ny, nz, this%tr_work, this%vhtr)
         call drain_wrap_face_x(this%uhtr, nx, ny, nz, nx_phys, ny_phys, nghost, &
                                per_x, per_y, fold_n)
         call drain_wrap_face_y(this%vhtr, nx, ny, nz, nx_phys, ny_phys, nghost, &
                                per_x, per_y, fold_n)
      end do

      ! ---- Reconstruct hprev = areaT·h_end + div(uhtr) (volume → thickness) ----
      ! + vanishing-layer hatch.  Done on the FROZEN hTr grid: the cell's
      ! window-start thickness, against which the frozen tracer mass is a
      ! consistent concentration.
      call drain_reconstruct_hprev(nx, ny, nz, metrics%areaT, metrics%iareaT, &
                                   ms%h_layer, this%uhtr, this%vhtr, this%hprev_work)
      call drain_wrap_centre(this%hprev_work, nx, ny, nz, nx_phys, ny_phys, nghost, &
                             per_x, per_y, fold_n)

      ! Undo the per-stage concentration hold applied by
      ! `continuity_tracer_step_split`: `hTr` currently carries the frozen
      ! window-start concentration weighted by the CURRENT (window-end)
      ! thickness, and the sub-cycle below needs it weighted by the
      ! window-START thickness `hprev` it just reconstructed.  This restores
      ! exactly the content the drain consumed before the concentration hold
      ! existed, so the drain — and everything it produces — is unchanged.
      if (this%hTr_holds_conc .and. allocated(ms%tracers)) then
         do it = 1, size(ms%tracers)
            if (.not. allocated(ms%tracers(it)%hTr)) cycle
            if (.not. ms%tracers(it)%do_horizontal_advection) cycle
            ! Post-`rk2_average` ⇒ POST_AVERAGE weight (see the parameter's
            ! comment).  Together with the in-stage holds recorded above, the
            ! hold/un-hold pair cancels exactly in the accumulator.
            select case (ms%tracers(it)%budget_id)
            case (TRACER_BUDGET_HEAT)
               call drain_rescale_hTr_budget(nx, ny, nz, this%h_win_start, ms%h_layer, &
                                             DRAIN_BUDGET_POST_AVERAGE_WEIGHT, &
                                             ms%tracers(it)%hTr, ms%heat_budget_horiz_adv)
            case (TRACER_BUDGET_SALT)
               call drain_rescale_hTr_budget(nx, ny, nz, this%h_win_start, ms%h_layer, &
                                             DRAIN_BUDGET_POST_AVERAGE_WEIGHT, &
                                             ms%tracers(it)%hTr, ms%salt_budget_horiz_adv)
            case default
               call drain_rescale_hTr(nx, ny, nz, this%h_win_start, ms%h_layer, &
                                      ms%tracers(it)%hTr)
            end select
         end do
         this%hTr_holds_conc = .false.
      end if

      ! ---- Seam-wrap the frozen tracer mass before the first reconstruction ----
      ! The PPM parabola (drain_parabola_*) reads the seam-adjacent cells'
      ! tracer concentration Tr = hTr/hprev at the ghost band (a ±2 stencil).
      ! hprev is already periodic/fold-consistent (wrapped above), but the
      ! frozen hTr the drain inherits is NOT guaranteed periodic in its ghost
      ! halo (the dynamics wrap h_layer, not necessarily the windowed hTr), so
      ! a seam-physical cell's reconstruction read a different Tr than its
      ! translated interior image — breaking bit-exact translation invariance
      ! (the per-pass updates re-wrap hTr at the end, but the FIRST pass's
      ! parabola already consumed the stale ghosts).  Wrap every advected
      ! tracer's hTr once here so pass 1 sees the same periodic/fold ghosts
      ! the later passes do.  No-op on non-periodic walls (bit-identical).
      do it = 1, size(ms%tracers)
         if (.not. allocated(ms%tracers(it)%hTr)) cycle
         if (.not. ms%tracers(it)%do_horizontal_advection) cycle
         call drain_wrap_centre(ms%tracers(it)%hTr, nx, ny, nz, nx_phys, ny_phys, &
                                nghost, per_x, per_y, fold_n)
      end do

      ! ---- Seed remaining transport from the accumulators ----
      call drain_copy_3d(nx + 1, ny, nz, this%uhtr, this%uhr_x)
      call drain_copy_3d(nx, ny + 1, nz, this%vhtr, this%uhr_y)

      ! ---- Fixed-budget sub-cycle, budget sized to the ACTUAL courant ----
      ! The worst case is `2*ratio+1` (per-step CFL ~ 1), but a real run's
      ! accumulated tracer courant is much smaller, so most of those passes are
      ! exact no-ops (once a face's transport is drained the hup/hlos limiter
      ! yields uhh=0 and the pass changes nothing).  Reduce the max per-cell
      ! accumulated courant ONCE (grid-uniform — every column runs the same
      ! `max_iter`, so no warp divergence; no per-pass early-exit) and size the
      ! fixed loop with the SAME `2*ceil(C)+1` form, capped at the worst case.
      ! Dropping the no-op passes is BIT-IDENTICAL (verified vs `2*ratio+1`).
      cfl_max = 0.0_wp
      do concurrent(k=1:nz, j=1:ny, i=1:nx) reduce(max:cfl_max)
         denom = metrics%areaT(i, j)*max(this%hprev_work(i, j, k), this%h_min)
         fmax = max(abs(this%uhr_x(i, j, k)), abs(this%uhr_x(i + 1, j, k)), &
                    abs(this%uhr_y(i, j, k)), abs(this%uhr_y(i, j + 1, k)))
         cfl_max = max(cfl_max, fmax/denom)
      end do
      ! Passes needed = ceiling(cfl_max / per-pass capacity).  The drain_limit_x
      ! two-test limiter guarantees each pass moves at least 0.5·hup (the
      ! `max(0.5*hup, hup-hlos)` floor), i.e. worst-case capacity 0.5 CFL/pass,
      ! so ceiling(2·cfl_max) passes are provably sufficient (and exact in the
      ! divergent worst case).  This is much tighter than the old
      ! `2*ceiling(cfl_max)+1`: at the default ratio = 1 a stable step has
      ! cfl < 0.5 ⇒ ONE pass (was 3 — two full-grid no-op passes).  Still
      ! grid-uniform (one reduction, no per-pass divergence) and capped at the
      ! worst case `2*ratio+1`.
      max_iter = ceiling(2.0_wp*cfl_max)
      if (max_iter > 2*ratio + 1) max_iter = 2*ratio + 1
      if (max_iter < 1) max_iter = 1
      do ipass = 1, max_iter
         ! Zonal sub-pass --------------------------------------------------
         call drain_limit_x(nx, ny, nz, metrics%areaT, this%h_min, &
                            this%uhr_x, this%hprev_work, this%uhh_x)
         ! Wrap the limited per-pass transport so the swept flux at the
         ! periodic seam / north fold reads the SAME uhh from both sides of
         ! the wrapped face (drain_limit_x writes only interior faces 2..nx
         ! + zeros the array edges; without the wrap the two images of a
         ! seam face disagree ⇒ the flux divergence does not telescope and
         ! Σ(areaT·hTr) leaks at the seam).  No-op on non-periodic walls.
         call drain_wrap_face_x(this%uhh_x, nx, ny, nz, nx_phys, ny_phys, nghost, &
                                per_x, per_y, fold_n)
         do it = 1, size(ms%tracers)
            if (.not. allocated(ms%tracers(it)%hTr)) cycle
            if (.not. ms%tracers(it)%do_horizontal_advection) cycle
            if (this%tracer_recon == TRACER_RECON_PPM) then
               call drain_parabola_x(nx, ny, nz, metrics%wet_T, ms%tracers(it)%hTr, this%hprev_work, &
                                     this%tr_work, this%pal, this%par, this%pa6)
               ! Wrap the parabola coefficients so a seam-face donor that lands
               ! in a ghost column carries the SAME (full-PPM) reconstruction as
               ! its physical image — the ghost band is otherwise PCM, which
               ! makes the two seam-face flux evaluations disagree (seam leak).
               ! Batch the 3 independent parabola-coeff wraps async on queue 1,
               ! sync once before the swept flux reads them.
               call drain_wrap_centre(this%pal, nx, ny, nz, nx_phys, ny_phys, nghost, &
                                      per_x, per_y, fold_n, no_wait=.true.)
               call drain_wrap_centre(this%par, nx, ny, nz, nx_phys, ny_phys, nghost, &
                                      per_x, per_y, fold_n, no_wait=.true.)
               call drain_wrap_centre(this%pa6, nx, ny, nz, nx_phys, ny_phys, nghost, &
                                      per_x, per_y, fold_n, no_wait=.true.)
               !$acc wait(1)
               call drain_swept_flux_x(nx, ny, nz, metrics%areaT, this%uhh_x, &
                                       this%hprev_work, this%pal, this%par, this%pa6, &
                                       this%tr_flux_x)
            else
               ! Q6 WENO drain: reconstruct the swept-average donor-edge
               ! concentration with the WENO ladder instead of the CW
               ! parabola.  Build Tr = hTr/hprev, wrap it so seam-adjacent
               ! stencils read consistent ghosts (the WENO reach is wider
               ! than PPM's — the wrap must cover the widest rung), then
               ! evaluate the swept face flux directly from Tr.
               call drain_fill_conc(nx, ny, nz, ms%tracers(it)%hTr, this%hprev_work, this%tr_work)
               call drain_wrap_centre(this%tr_work, nx, ny, nz, nx_phys, ny_phys, nghost, &
                                      per_x, per_y, fold_n)
               call drain_swept_flux_x_weno(nx, ny, nz, nghost, per_x, metrics%areaT, &
                                            this%uhh_x, this%hprev_work, this%tr_work, &
                                            metrics%wet_T, this%tracer_recon, this%tr_flux_x)
            end if
            select case (ms%tracers(it)%budget_id)
            case (TRACER_BUDGET_HEAT)
               call drain_update_tracer_x_budget(nx, ny, nz, metrics%iareaT, this%tr_flux_x, &
                                                 DRAIN_BUDGET_POST_AVERAGE_WEIGHT, &
                                                 ms%tracers(it)%hTr, ms%heat_budget_horiz_adv)
            case (TRACER_BUDGET_SALT)
               call drain_update_tracer_x_budget(nx, ny, nz, metrics%iareaT, this%tr_flux_x, &
                                                 DRAIN_BUDGET_POST_AVERAGE_WEIGHT, &
                                                 ms%tracers(it)%hTr, ms%salt_budget_horiz_adv)
            case default
               call drain_update_tracer_x(nx, ny, nz, metrics%iareaT, this%tr_flux_x, &
                                          ms%tracers(it)%hTr)
            end select
            call drain_wrap_centre(ms%tracers(it)%hTr, nx, ny, nz, nx_phys, ny_phys, &
                                   nghost, per_x, per_y, fold_n)
         end do
         ! Advance the (shared) thickness + remaining transport once.
         call drain_update_h_x(nx, ny, nz, metrics%iareaT, this%uhh_x, this%hprev_work)
         call drain_subtract_3d(nx + 1, ny, nz, this%uhh_x, this%uhr_x)
         call drain_wrap_centre(this%hprev_work, nx, ny, nz, nx_phys, ny_phys, nghost, &
                                per_x, per_y, fold_n)
         call drain_wrap_face_x(this%uhr_x, nx, ny, nz, nx_phys, ny_phys, nghost, &
                                per_x, per_y, fold_n)

         ! Meridional sub-pass ---------------------------------------------
         call drain_limit_y(nx, ny, nz, metrics%areaT, this%h_min, &
                            this%uhr_y, this%hprev_work, this%uhh_y)
         ! Wrap the limited meridional transport (see the zonal note above).
         call drain_wrap_face_y(this%uhh_y, nx, ny, nz, nx_phys, ny_phys, nghost, &
                                per_x, per_y, fold_n)
         do it = 1, size(ms%tracers)
            if (.not. allocated(ms%tracers(it)%hTr)) cycle
            if (.not. ms%tracers(it)%do_horizontal_advection) cycle
            if (this%tracer_recon == TRACER_RECON_PPM) then
               call drain_parabola_y(nx, ny, nz, metrics%wet_T, ms%tracers(it)%hTr, this%hprev_work, &
                                     this%tr_work, this%pal, this%par, this%pa6)
               ! Wrap the meridional parabola coefficients (see the zonal note).
               ! Batch the 3 independent parabola-coeff wraps async on queue 1,
               ! sync once before the swept flux reads them.
               call drain_wrap_centre(this%pal, nx, ny, nz, nx_phys, ny_phys, nghost, &
                                      per_x, per_y, fold_n, no_wait=.true.)
               call drain_wrap_centre(this%par, nx, ny, nz, nx_phys, ny_phys, nghost, &
                                      per_x, per_y, fold_n, no_wait=.true.)
               call drain_wrap_centre(this%pa6, nx, ny, nz, nx_phys, ny_phys, nghost, &
                                      per_x, per_y, fold_n, no_wait=.true.)
               !$acc wait(1)
               call drain_swept_flux_y(nx, ny, nz, metrics%areaT, this%uhh_y, &
                                       this%hprev_work, this%pal, this%par, this%pa6, &
                                       this%tr_flux_y)
            else
               ! Q6 WENO drain (meridional; see the zonal branch note).
               call drain_fill_conc(nx, ny, nz, ms%tracers(it)%hTr, this%hprev_work, this%tr_work)
               call drain_wrap_centre(this%tr_work, nx, ny, nz, nx_phys, ny_phys, nghost, &
                                      per_x, per_y, fold_n)
               call drain_swept_flux_y_weno(nx, ny, nz, nghost, per_y, metrics%areaT, &
                                            this%uhh_y, this%hprev_work, this%tr_work, &
                                            metrics%wet_T, this%tracer_recon, this%tr_flux_y)
            end if
            select case (ms%tracers(it)%budget_id)
            case (TRACER_BUDGET_HEAT)
               call drain_update_tracer_y_budget(nx, ny, nz, metrics%iareaT, this%tr_flux_y, &
                                                 DRAIN_BUDGET_POST_AVERAGE_WEIGHT, &
                                                 ms%tracers(it)%hTr, ms%heat_budget_horiz_adv)
            case (TRACER_BUDGET_SALT)
               call drain_update_tracer_y_budget(nx, ny, nz, metrics%iareaT, this%tr_flux_y, &
                                                 DRAIN_BUDGET_POST_AVERAGE_WEIGHT, &
                                                 ms%tracers(it)%hTr, ms%salt_budget_horiz_adv)
            case default
               call drain_update_tracer_y(nx, ny, nz, metrics%iareaT, this%tr_flux_y, &
                                          ms%tracers(it)%hTr)
            end select
            call drain_wrap_centre(ms%tracers(it)%hTr, nx, ny, nz, nx_phys, ny_phys, &
                                   nghost, per_x, per_y, fold_n)
         end do
         call drain_update_h_y(nx, ny, nz, metrics%iareaT, this%uhh_y, this%hprev_work)
         call drain_subtract_3d(nx, ny + 1, nz, this%uhh_y, this%uhr_y)
         call drain_wrap_centre(this%hprev_work, nx, ny, nz, nx_phys, ny_phys, nghost, &
                                per_x, per_y, fold_n)
         call drain_wrap_face_y(this%uhr_y, nx, ny, nz, nx_phys, ny_phys, nghost, &
                                per_x, per_y, fold_n)
      end do

      ! ---- Reset accumulators + clock for the next window ----
      call drain_zero_3d(nx + 1, ny, nz, this%uhtr)
      call drain_zero_3d(nx, ny + 1, nz, this%vhtr)
      this%t_dyn_rel_adv = 0.0_wp
   end subroutine continuity_tracer_drain

   ! ====================================================================
   ! Drain seam helpers (single-rank periodic + north-fold wraps).
   ! ====================================================================

   subroutine drain_wrap_centre(fld, nx, ny, nz, nx_phys, ny_phys, nghost, &
                                per_x, per_y, fold_n, no_wait)
      !! Periodic wrap (+ north fold) of a cell-centred drain field.
      !! `no_wait` (optional): when .true. AND not folding, the periodic wrap
      !! is issued async on queue 1 without syncing, so a caller can batch
      !! several independent wraps (e.g. the pal/par/pa6 parabola triple) and
      !! `!$acc wait(1)` once.  Ignored when fold_n (the fold reads the wrapped
      !! field, so the periodic wrap must complete first).
      integer, intent(in) :: nx, ny, nz, nx_phys, ny_phys, nghost
      real(wp), intent(inout) :: fld(nx, ny, nz)
      logical, intent(in) :: per_x, per_y, fold_n
      logical, intent(in), optional :: no_wait
      logical :: nw
      nw = .false.
      if (present(no_wait)) nw = no_wait .and. .not. fold_n
      if (per_x .or. per_y) then
         call ocean_periodic_wrap_centre_3d(fld, nx, ny, nz, &
                                            nx_phys, ny_phys, nghost, per_x, per_y, no_wait=nw)
      end if
      if (fold_n) call fold_north_centre(fld, nx, ny, nz, nx_phys, ny_phys, nghost)
   end subroutine drain_wrap_centre

   subroutine drain_wrap_face_x(fld, nx, ny, nz, nx_phys, ny_phys, nghost, &
                                per_x, per_y, fold_n)
      !! Periodic wrap (+ north fold) of an x-face drain field (nx+1,ny,nz).
      integer, intent(in) :: nx, ny, nz, nx_phys, ny_phys, nghost
      real(wp), intent(inout) :: fld(nx + 1, ny, nz)
      logical, intent(in) :: per_x, per_y, fold_n
      if (per_x .or. per_y) then
         call ocean_periodic_wrap_face_x_3d(fld, nx + 1, ny, nz, &
                                            nx_phys, ny_phys, nghost, per_x, per_y)
      end if
      if (fold_n) call fold_north_u_face(fld, nx + 1, ny, nz, nx_phys, ny_phys, nghost)
   end subroutine drain_wrap_face_x

   subroutine drain_wrap_face_y(fld, nx, ny, nz, nx_phys, ny_phys, nghost, &
                                per_x, per_y, fold_n)
      !! Periodic wrap (+ north fold) of a y-face drain field (nx,ny+1,nz).
      integer, intent(in) :: nx, ny, nz, nx_phys, ny_phys, nghost
      real(wp), intent(inout) :: fld(nx, ny + 1, nz)
      logical, intent(in) :: per_x, per_y, fold_n
      if (per_x .or. per_y) then
         call ocean_periodic_wrap_face_y_3d(fld, nx, ny + 1, nz, &
                                            nx_phys, ny_phys, nghost, per_x, per_y)
      end if
      if (fold_n) call fold_north_v_face(fld, nx, ny + 1, nz, nx_phys, ny_phys, nghost)
   end subroutine drain_wrap_face_y

   ! ====================================================================
   ! Drain compute kernels (flat-impl, explicit-shape, do concurrent).
   ! Face convention: x-face i lies between cell (i-1) and cell (i);
   ! positive uhtr(i) ⇒ rightward ⇒ donor = cell (i-1).  Mirrors the
   ! every-step advect (divergence hTr(i) -= (F(i+1)-F(i))·iareaT).
   ! Transport (uhr/uhh) and the limiter work in area-weighted VOLUME
   ! (m^3); thickness `hprev_work` is m; `Vprev = areaT·hprev_work`.
   ! ====================================================================

   pure subroutine drain_reconstruct_hprev(nx, ny, nz, areaT, iareaT, h_end, &
                                           uhtr, vhtr, hprev)
      !! hprev = max(0, areaT·h_end + div(uhtr,vhtr)) · iareaT, then the
      !! vanishing-layer hatch `hprev += max(0, 1e-13·hprev − h_end)`
      !! (Adcroft & Hallberg 2006; reuse of VANISHING_LAYER_TOL thinking).
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: areaT(nx, ny), iareaT(nx, ny)
      real(wp), intent(in) :: h_end(nx, ny, nz)
      real(wp), intent(in) :: uhtr(nx + 1, ny, nz)
      real(wp), intent(in) :: vhtr(nx, ny + 1, nz)
      real(wp), intent(inout) :: hprev(nx, ny, nz)
      integer :: i, j, k
      real(wp) :: vol, eps_h
      do concurrent(k=1:nz, j=1:ny, i=1:nx) local(vol, eps_h)
         vol = areaT(i, j)*h_end(i, j, k) &
               + (uhtr(i + 1, j, k) - uhtr(i, j, k)) &
               + (vhtr(i, j + 1, k) - vhtr(i, j, k))
         hprev(i, j, k) = max(0.0_wp, vol)*iareaT(i, j)
         eps_h = max(0.0_wp, 1.0e-13_wp*hprev(i, j, k) - h_end(i, j, k))
         hprev(i, j, k) = hprev(i, j, k) + eps_h
      end do
   end subroutine drain_reconstruct_hprev

   pure subroutine drain_avail_limit(nx, ny, nz, areaT, h_min, h_end, &
                                     uhtr, vhtr, scratch)
      !! Conservative upfront availability limiter on the accumulated window
      !! transports `uhtr/vhtr`, applied BEFORE drain_reconstruct_hprev.
      !!
      !! Guarantees, for every cell, that the reconstructed window-start
      !! volume stays positive (≥ areaT·h_min):
      !!
      !!   vol(i,j,k) = areaT·h_end + (uhtr(i+1)-uhtr(i)) + (vhtr(j+1)-vhtr(j))
      !!              ≥ areaT·h_min
      !!
      !! so the non-conservative max(0,·) clamp + vanishing-layer hatch in
      !! drain_reconstruct_hprev never fire.  Without this, the combined
      !! Fox-Kemper (2008,2011) overturning + resolved/barotropic transport
      !! can overdraw a thin z* surface layer (vol < 0), the clamp destroys
      !! volume and the hatch re-inflates hprev WITHOUT matching tracer mass
      !! → ~5%/day tracer leak.  MOM6 avoids this by sizing its availability
      !! cap against the combined transport; this is the windowed-drain
      !! analogue (general, protects against any overdraw source).
      !!
      !! Mechanism (Zalesak/FCT-style inflow scaling, GPU-uniform fixed
      !! budget): the cell deficit is covered by shrinking the cell's
      !! INFLOW-face transports.  Reducing an inflow magnitude raises the
      !! receiving cell's vol AND the donor neighbour's vol (it is that
      !! neighbour's outflow), so the iteration is monotone in every cell's
      !! vol and converges.  Scaling is applied multiplicatively to the
      !! transport, so the limited uhtr/vhtr are used CONSISTENTLY for both
      !! the reconstruction and the sub-cycle ⇒ the drain still telescopes
      !! exactly onto h_end and Σ(areaT·hTr) is conserved to round-off.
      !!
      !! Inert when every vol > areaT·h_min already (scale = 1 everywhere),
      !! so ratio=1 / FK-off / non-overdrawing windows are bit-identical.
      !!
      !! ONE FCT pass — the caller loops it `AVAIL_MAX_PASS` times with a
      !! seam re-wrap of `scratch` and the faces between passes so the
      !! constraint diffuses consistently across periodic/fold ghosts.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: areaT(nx, ny)
      real(wp), intent(in) :: h_min
      real(wp), intent(in) :: h_end(nx, ny, nz)
      real(wp), intent(inout) :: uhtr(nx + 1, ny, nz)
      real(wp), intent(inout) :: vhtr(nx, ny + 1, nz)
      real(wp), intent(inout) :: scratch(nx, ny, nz)
         !! Per-cell inflow scale factor in [0,1] (centre-shaped slot).
      integer :: i, j, k
      real(wp) :: vol, vmin, inflow, deficit, sfac
      ! ---- Per-cell inflow scale factor ----
      do concurrent(k=1:nz, j=1:ny, i=1:nx) &
         local(vol, vmin, inflow, deficit)
         vmin = areaT(i, j)*h_min
         vol = areaT(i, j)*h_end(i, j, k) &
               + (uhtr(i + 1, j, k) - uhtr(i, j, k)) &
               + (vhtr(i, j + 1, k) - vhtr(i, j, k))
         ! Inflow into this cell across its four faces (volume, ≥ 0):
         !   west  face uhtr(i)   inflow if > 0
         !   east  face uhtr(i+1) inflow if < 0
         !   south face vhtr(j)   inflow if > 0
         !   north face vhtr(j+1) inflow if < 0
         inflow = max(0.0_wp, uhtr(i, j, k)) &
                  + max(0.0_wp, -uhtr(i + 1, j, k)) &
                  + max(0.0_wp, vhtr(i, j, k)) &
                  + max(0.0_wp, -vhtr(i, j + 1, k))
         if (vol < vmin .and. inflow > 0.0_wp) then
            deficit = vmin - vol
            scratch(i, j, k) = max(0.0_wp, (inflow - deficit)/inflow)
         else
            scratch(i, j, k) = 1.0_wp
         end if
      end do
   end subroutine drain_avail_limit

   pure subroutine drain_avail_scale_x(nx, ny, nz, scratch, uhtr)
      !! Scale each interior x-face transport by its INFLOW-receiving cell's
      !! factor (drain_avail_limit step 2).  Face i between cell (i-1) and
      !! cell (i): uhtr(i)>0 ⇒ receiver i, uhtr(i)<0 ⇒ receiver i-1.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: scratch(nx, ny, nz)
      real(wp), intent(inout) :: uhtr(nx + 1, ny, nz)
      integer :: i, j, k
      real(wp) :: sfac
      do concurrent(k=1:nz, j=1:ny, i=2:nx) local(sfac)
         if (uhtr(i, j, k) > 0.0_wp) then
            sfac = scratch(i, j, k)
         else if (uhtr(i, j, k) < 0.0_wp) then
            sfac = scratch(i - 1, j, k)
         else
            sfac = 1.0_wp
         end if
         uhtr(i, j, k) = uhtr(i, j, k)*sfac
      end do
   end subroutine drain_avail_scale_x

   pure subroutine drain_avail_scale_y(nx, ny, nz, scratch, vhtr)
      !! Meridional analogue of drain_avail_scale_x.  Face j between cell
      !! (j-1) and cell (j): vhtr(j)>0 ⇒ receiver j, vhtr(j)<0 ⇒ receiver j-1.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: scratch(nx, ny, nz)
      real(wp), intent(inout) :: vhtr(nx, ny + 1, nz)
      integer :: i, j, k
      real(wp) :: sfac
      do concurrent(k=1:nz, j=2:ny, i=1:nx) local(sfac)
         if (vhtr(i, j, k) > 0.0_wp) then
            sfac = scratch(i, j, k)
         else if (vhtr(i, j, k) < 0.0_wp) then
            sfac = scratch(i, j - 1, k)
         else
            sfac = 1.0_wp
         end if
         vhtr(i, j, k) = vhtr(i, j, k)*sfac
      end do
   end subroutine drain_avail_scale_y

   pure subroutine drain_rescale_hTr(nx, ny, nz, h_new, h_old, hTr)
      !! Re-weight a tracer's thickness-weighted content onto a new layer
      !! thickness, holding the CONCENTRATION fixed:
      !!
      !!   hTr := hTr · h_new / h_old      (⇒ hTr/h_new ≡ hTr/h_old)
      !!
      !! Used by the windowed (`dt_tracer_advect_ratio > 1`) path, whose
      !! prognostic is the CONTENT `hTr` while MOM6's is the CONCENTRATION
      !! `Tr`.  Freezing a content across a window over which continuity
      !! keeps advancing `h` silently corrupts every consumer that derives
      !! `T = hTr/h_layer` — the EOS above all — by `δT/T = −δh/h`.  Freezing
      !! a concentration (what MOM6 does) does not.  See
      !! `continuity_tracer_step_split` / `continuity_tracer_drain`.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: h_new(nx, ny, nz), h_old(nx, ny, nz)
      real(wp), intent(inout) :: hTr(nx, ny, nz)
      integer :: i, j, k
      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         hTr(i, j, k) = hTr(i, j, k)*h_new(i, j, k)/max(h_old(i, j, k), DRAIN_MIN_H)
      end do
   end subroutine drain_rescale_hTr

   pure subroutine drain_rescale_hTr_budget(nx, ny, nz, h_new, h_old, w, hTr, budget_adv)
      !! `drain_rescale_hTr` + the closed-budget fill.  The concentration
      !! hold / un-hold pair is NOT content-conserving cell by cell --
      !! `hTr := Tr·h` moves `Σ areaT·hTr` by `Σ areaT·Tr·δh`, which is only
      !! zero when `Tr` is uniform -- so both halves have to be recorded or
      !! the closed budget is only valid on window boundaries and wobbles at
      !! every mid-window report.  Recording both makes them cancel exactly,
      !! since the un-hold is the arithmetic inverse of the accumulated hold.
      !!
      !! `w` is `DRAIN_BUDGET_IN_STAGE_WEIGHT` for the per-stage hold inside
      !! `continuity_tracer_step_split` and
      !! `DRAIN_BUDGET_POST_AVERAGE_WEIGHT` for the un-hold in
      !! `continuity_tracer_drain` (which runs after `rk2_average`).
      !!
      !! The `hTr` expression is character-for-character the budget-free
      !! twin's, so the prognostic path stays bit-identical.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: h_new(nx, ny, nz), h_old(nx, ny, nz)
      real(wp), intent(in) :: w
      real(wp), intent(inout) :: hTr(nx, ny, nz)
      real(wp), intent(inout) :: budget_adv(nx, ny, nz)
      integer :: i, j, k
      real(wp) :: hnew_tr
      do concurrent(k=1:nz, j=1:ny, i=1:nx) local(hnew_tr)
         hnew_tr = hTr(i, j, k)*h_new(i, j, k)/max(h_old(i, j, k), DRAIN_MIN_H)
         budget_adv(i, j, k) = budget_adv(i, j, k) + w*(hnew_tr - hTr(i, j, k))
         hTr(i, j, k) = hnew_tr
      end do
   end subroutine drain_rescale_hTr_budget

   pure subroutine drain_copy_3d(n1, n2, n3, src, dst)
      !! dst = src (explicit-shape device copy).
      integer, intent(in) :: n1, n2, n3
      real(wp), intent(in) :: src(n1, n2, n3)
      real(wp), intent(inout) :: dst(n1, n2, n3)
      integer :: i, j, k
      do concurrent(k=1:n3, j=1:n2, i=1:n1)
         dst(i, j, k) = src(i, j, k)
      end do
   end subroutine drain_copy_3d

   pure subroutine drain_zero_3d(n1, n2, n3, fld)
      !! fld = 0.
      integer, intent(in) :: n1, n2, n3
      real(wp), intent(inout) :: fld(n1, n2, n3)
      integer :: i, j, k
      do concurrent(k=1:n3, j=1:n2, i=1:n1)
         fld(i, j, k) = 0.0_wp
      end do
   end subroutine drain_zero_3d

   pure subroutine drain_subtract_3d(n1, n2, n3, sub, fld)
      !! fld = fld - sub  (uhr -= uhh).
      integer, intent(in) :: n1, n2, n3
      real(wp), intent(in) :: sub(n1, n2, n3)
      real(wp), intent(inout) :: fld(n1, n2, n3)
      integer :: i, j, k
      do concurrent(k=1:n3, j=1:n2, i=1:n1)
         fld(i, j, k) = fld(i, j, k) - sub(i, j, k)
      end do
   end subroutine drain_subtract_3d

   pure subroutine drain_limit_x(nx, ny, nz, areaT, h_min, uhr_x, hprev, uhh_x)
      !! MOM6 hup/hlos/min_h two-test limiter on the zonal face transport
      !! (volume units).  Face i between cell (i-1) and cell (i).
      !! Positive flow (uhr_x(i) > 0), donor = cell (i-1):
      !!   hup  = areaT(i-1)·hprev(i-1) − areaT(i-1)·min_h
      !!   hlos = max(0, −uhr_x(i-1))   (already-committed outflow via the
      !!          donor's OTHER (west) face)
      !!   cap when (hup−hlos)−uhr < 0 AND 0.5·hup−uhr < 0.
      !! Negative flow mirror, donor = cell (i).
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: areaT(nx, ny)
      real(wp), intent(in) :: h_min
      real(wp), intent(in) :: uhr_x(nx + 1, ny, nz)
      real(wp), intent(in) :: hprev(nx, ny, nz)
      real(wp), intent(inout) :: uhh_x(nx + 1, ny, nz)
      integer :: i, j, k
      real(wp) :: uhr, hup, hlos
      ! Interior faces 2..nx (each has a left cell i-1 and right cell i).
      do concurrent(k=1:nz, j=1:ny, i=2:nx) local(uhr, hup, hlos)
         uhr = uhr_x(i, j, k)
         if (uhr > 0.0_wp) then
            hup = areaT(i - 1, j)*hprev(i - 1, j, k) - areaT(i - 1, j)*h_min
            hlos = max(0.0_wp, -uhr_x(i - 1, j, k))
            if (((hup - hlos) - uhr < 0.0_wp) .and. (0.5_wp*hup - uhr < 0.0_wp)) then
               uhh_x(i, j, k) = max(0.0_wp, max(0.5_wp*hup, hup - hlos))
            else
               uhh_x(i, j, k) = uhr
            end if
         else if (uhr < 0.0_wp) then
            hup = areaT(i, j)*hprev(i, j, k) - areaT(i, j)*h_min
            hlos = max(0.0_wp, uhr_x(i + 1, j, k))
            if (((hup - hlos) + uhr < 0.0_wp) .and. (0.5_wp*hup + uhr < 0.0_wp)) then
               uhh_x(i, j, k) = -max(0.0_wp, max(0.5_wp*hup, hup - hlos))
            else
               uhh_x(i, j, k) = uhr
            end if
         else
            uhh_x(i, j, k) = 0.0_wp
         end if
      end do
      ! Boundary faces 1 and nx+1: no interior donor on one side ⇒ no flux.
      do concurrent(k=1:nz, j=1:ny)
         uhh_x(1, j, k) = 0.0_wp
         uhh_x(nx + 1, j, k) = 0.0_wp
      end do
   end subroutine drain_limit_x

   pure subroutine drain_limit_y(nx, ny, nz, areaT, h_min, uhr_y, hprev, uhh_y)
      !! Meridional analogue of drain_limit_x.  Face j between cell (i,j-1)
      !! and cell (i,j); positive donor = cell (i,j-1).
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: areaT(nx, ny)
      real(wp), intent(in) :: h_min
      real(wp), intent(in) :: uhr_y(nx, ny + 1, nz)
      real(wp), intent(in) :: hprev(nx, ny, nz)
      real(wp), intent(inout) :: uhh_y(nx, ny + 1, nz)
      integer :: i, j, k
      real(wp) :: uhr, hup, hlos
      do concurrent(k=1:nz, j=2:ny, i=1:nx) local(uhr, hup, hlos)
         uhr = uhr_y(i, j, k)
         if (uhr > 0.0_wp) then
            hup = areaT(i, j - 1)*hprev(i, j - 1, k) - areaT(i, j - 1)*h_min
            hlos = max(0.0_wp, -uhr_y(i, j - 1, k))
            if (((hup - hlos) - uhr < 0.0_wp) .and. (0.5_wp*hup - uhr < 0.0_wp)) then
               uhh_y(i, j, k) = max(0.0_wp, max(0.5_wp*hup, hup - hlos))
            else
               uhh_y(i, j, k) = uhr
            end if
         else if (uhr < 0.0_wp) then
            hup = areaT(i, j)*hprev(i, j, k) - areaT(i, j)*h_min
            hlos = max(0.0_wp, uhr_y(i, j + 1, k))
            if (((hup - hlos) + uhr < 0.0_wp) .and. (0.5_wp*hup + uhr < 0.0_wp)) then
               uhh_y(i, j, k) = -max(0.0_wp, max(0.5_wp*hup, hup - hlos))
            else
               uhh_y(i, j, k) = uhr
            end if
         else
            uhh_y(i, j, k) = 0.0_wp
         end if
      end do
      do concurrent(k=1:nz, i=1:nx)
         uhh_y(i, 1, k) = 0.0_wp
         uhh_y(i, ny + 1, k) = 0.0_wp
      end do
   end subroutine drain_limit_y

   pure subroutine drain_parabola_x(nx, ny, nz, wet_T, hTr, hprev, tr, aL, aR, a6)
      !! Rebuild the per-cell zonal CW PPM parabola from the CURRENT Tr =
      !! hTr/hprev (V2 — per pass).  Interior cells (3..nx-2) use the
      !! limited PPM edges; the 2-cell boundary band falls back to PCM
      !! (aL=aR=Tr ⇒ swept reduces to the donor value, 1st order),
      !! matching tracer_advect_zonal_one_impl's near-wall band.
      !! TODO(MOM6-fidelity): MOM6 advect_tracer keeps full PPM up to the wall
      !! (dropping to PCM only at genuine local extrema / zero `mask2dCu`
      !! faces), so we are 1st-order in the 2 cells nearest a true WALL where
      !! MOM6 is PPM-with-mask (periodic seams are fine — the wrap restores
      !! full PPM).  Tracked divergence; revisit if near-wall tracer
      !! diffusion matters.
      !! a6 = 6·Tr − 3·(aL+aR).  Mirror-T at land neighbours (C2);
      !! bit-identical for all-wet (`wet_T≡1`).
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: wet_T(nx, ny)
      real(wp), intent(in) :: hTr(nx, ny, nz), hprev(nx, ny, nz)
      real(wp), intent(inout) :: tr(nx, ny, nz), aL(nx, ny, nz), aR(nx, ny, nz)
      real(wp), intent(inout) :: a6(nx, ny, nz)
      integer :: i, j, k
      real(wp) :: Tr_m2, Tr_m1, Tr_0, Tr_p1, Tr_p2
      real(wp) :: dh_m1, dh_0, dh_p1, Tr_left, Tr_right

      ! Concentration field (guard zero thickness).
      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         tr(i, j, k) = hTr(i, j, k)/max(hprev(i, j, k), DRAIN_MIN_H)
      end do
      ! Interior PPM parabola.
      do concurrent(k=1:nz, j=1:ny, i=3:nx - 2) &
         local(Tr_m2, Tr_m1, Tr_0, Tr_p1, Tr_p2, dh_m1, dh_0, dh_p1, Tr_left, Tr_right)
         Tr_0 = tr(i, j, k)
         Tr_m1 = ppm_mirror_h(tr(i - 1, j, k), Tr_0, wet_T(i - 1, j))
         Tr_p1 = ppm_mirror_h(tr(i + 1, j, k), Tr_0, wet_T(i + 1, j))
         Tr_m2 = ppm_mirror_h(tr(i - 2, j, k), Tr_m1, wet_T(i - 2, j))
         Tr_p2 = ppm_mirror_h(tr(i + 2, j, k), Tr_p1, wet_T(i + 2, j))
         call ppm_limited_slope(Tr_m2, Tr_m1, Tr_0, dh_m1)
         call ppm_limited_slope(Tr_m1, Tr_0, Tr_p1, dh_0)
         call ppm_limited_slope(Tr_0, Tr_p1, Tr_p2, dh_p1)
         dh_0 = dh_0*wet_T(i - 1, j)*wet_T(i, j)*wet_T(i + 1, j)
         Tr_left = 0.5_wp*(Tr_m1 + Tr_0) - (dh_0 - dh_m1)/6.0_wp
         Tr_right = 0.5_wp*(Tr_0 + Tr_p1) - (dh_p1 - dh_0)/6.0_wp
         call ppm_cell_limiter(Tr_0, Tr_left, Tr_right)
         aL(i, j, k) = Tr_left
         aR(i, j, k) = Tr_right
         a6(i, j, k) = 6.0_wp*Tr_0 - 3.0_wp*(Tr_left + Tr_right)
      end do
      ! Boundary band (i=1,2 and nx-1,nx): PCM donor (aL=aR=Tr, a6=0).
      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         if (i <= 2 .or. i >= nx - 1) then
            aL(i, j, k) = tr(i, j, k)
            aR(i, j, k) = tr(i, j, k)
            a6(i, j, k) = 0.0_wp
         end if
      end do
   end subroutine drain_parabola_x

   pure subroutine drain_parabola_y(nx, ny, nz, wet_T, hTr, hprev, tr, aL, aR, a6)
      !! Meridional analogue of drain_parabola_x.  aL = south-edge,
      !! aR = north-edge value of each cell.  Mirror-T at land
      !! neighbours (C2); bit-identical for all-wet.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: wet_T(nx, ny)
      real(wp), intent(in) :: hTr(nx, ny, nz), hprev(nx, ny, nz)
      real(wp), intent(inout) :: tr(nx, ny, nz), aL(nx, ny, nz), aR(nx, ny, nz)
      real(wp), intent(inout) :: a6(nx, ny, nz)
      integer :: i, j, k
      real(wp) :: Tr_m2, Tr_m1, Tr_0, Tr_p1, Tr_p2
      real(wp) :: dh_m1, dh_0, dh_p1, Tr_left, Tr_right

      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         tr(i, j, k) = hTr(i, j, k)/max(hprev(i, j, k), DRAIN_MIN_H)
      end do
      do concurrent(k=1:nz, j=3:ny - 2, i=1:nx) &
         local(Tr_m2, Tr_m1, Tr_0, Tr_p1, Tr_p2, dh_m1, dh_0, dh_p1, Tr_left, Tr_right)
         Tr_0 = tr(i, j, k)
         Tr_m1 = ppm_mirror_h(tr(i, j - 1, k), Tr_0, wet_T(i, j - 1))
         Tr_p1 = ppm_mirror_h(tr(i, j + 1, k), Tr_0, wet_T(i, j + 1))
         Tr_m2 = ppm_mirror_h(tr(i, j - 2, k), Tr_m1, wet_T(i, j - 2))
         Tr_p2 = ppm_mirror_h(tr(i, j + 2, k), Tr_p1, wet_T(i, j + 2))
         call ppm_limited_slope(Tr_m2, Tr_m1, Tr_0, dh_m1)
         call ppm_limited_slope(Tr_m1, Tr_0, Tr_p1, dh_0)
         call ppm_limited_slope(Tr_0, Tr_p1, Tr_p2, dh_p1)
         dh_0 = dh_0*wet_T(i, j - 1)*wet_T(i, j)*wet_T(i, j + 1)
         Tr_left = 0.5_wp*(Tr_m1 + Tr_0) - (dh_0 - dh_m1)/6.0_wp
         Tr_right = 0.5_wp*(Tr_0 + Tr_p1) - (dh_p1 - dh_0)/6.0_wp
         call ppm_cell_limiter(Tr_0, Tr_left, Tr_right)
         aL(i, j, k) = Tr_left
         aR(i, j, k) = Tr_right
         a6(i, j, k) = 6.0_wp*Tr_0 - 3.0_wp*(Tr_left + Tr_right)
      end do
      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         if (j <= 2 .or. j >= ny - 1) then
            aL(i, j, k) = tr(i, j, k)
            aR(i, j, k) = tr(i, j, k)
            a6(i, j, k) = 0.0_wp
         end if
      end do
   end subroutine drain_parabola_y

   pure subroutine drain_swept_flux_x(nx, ny, nz, areaT, uhh, hprev, aL, aR, a6, F)
      !! MOM6 swept-average CW parabola flux for the zonal faces.
      !! Face i, donor = cell (i-1) if uhh>0 else cell (i).  Per-pass
      !! Courant CFL = |uhh| / (areaT·hprev) on the donor, clamped [0,1].
      !!   uhh >= 0: F = uhh·( aR − 0.5·CFL·((aR−aL) − a6·(1 − ⅔·CFL)) )
      !!   uhh <  0: F = uhh·( aL + 0.5·CFL·((aR−aL) + a6·(1 − ⅔·CFL)) )
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: areaT(nx, ny)
      real(wp), intent(in) :: uhh(nx + 1, ny, nz)
      real(wp), intent(in) :: hprev(nx, ny, nz)
      real(wp), intent(in) :: aL(nx, ny, nz), aR(nx, ny, nz), a6(nx, ny, nz)
      real(wp), intent(inout) :: F(nx + 1, ny, nz)
      integer :: i, j, k
      real(wp) :: u, cfl, vol, conc
      ! Array-edge faces (1, nx+1) have no donor cell -- the u>0 donor at
      ! face 1 is cell 0 and the u<0 donor at face nx+1 is cell nx+1, both
      ! off the [1,nx] cell arrays.  They feed only ghost cells the caller
      ! re-wraps, so set them to zero.  Explicit do concurrent (not the
      ! F(1,:,:)=0 array-section assignment, which does not reliably
      ! offload under stdpar).
      do concurrent(k=1:nz, j=1:ny)
         F(1, j, k) = 0.0_wp
         F(nx + 1, j, k) = 0.0_wp
      end do
      ! Interior faces 2..nx: the u>0 donor i-1 >= 1 and the u<0 donor i
      ! <= nx are both in range.  (Was i=1:nx+1, which read cell 0 / nx+1
      ! at the edge faces -- a latent OOB that only faults once the array
      ! is page-aligned at large grids.)
      do concurrent(k=1:nz, j=1:ny, i=2:nx) local(u, cfl, vol, conc)
         u = uhh(i, j, k)
         if (u > 0.0_wp) then
            vol = max(areaT(i - 1, j)*hprev(i - 1, j, k), DRAIN_MIN_VOL)
            cfl = min(u/vol, 1.0_wp)
            conc = aR(i - 1, j, k) - 0.5_wp*cfl* &
                   ((aR(i - 1, j, k) - aL(i - 1, j, k)) &
                    - a6(i - 1, j, k)*(1.0_wp - (2.0_wp/3.0_wp)*cfl))
            F(i, j, k) = u*conc
         else if (u < 0.0_wp) then
            vol = max(areaT(i, j)*hprev(i, j, k), DRAIN_MIN_VOL)
            cfl = min(-u/vol, 1.0_wp)
            conc = aL(i, j, k) + 0.5_wp*cfl* &
                   ((aR(i, j, k) - aL(i, j, k)) &
                    + a6(i, j, k)*(1.0_wp - (2.0_wp/3.0_wp)*cfl))
            F(i, j, k) = u*conc
         else
            F(i, j, k) = 0.0_wp
         end if
      end do
   end subroutine drain_swept_flux_x

   pure subroutine drain_swept_flux_y(nx, ny, nz, areaT, uhh, hprev, aL, aR, a6, F)
      !! Meridional analogue of drain_swept_flux_x.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: areaT(nx, ny)
      real(wp), intent(in) :: uhh(nx, ny + 1, nz)
      real(wp), intent(in) :: hprev(nx, ny, nz)
      real(wp), intent(in) :: aL(nx, ny, nz), aR(nx, ny, nz), a6(nx, ny, nz)
      real(wp), intent(inout) :: F(nx, ny + 1, nz)
      integer :: i, j, k
      real(wp) :: u, cfl, vol, conc
      ! Array-edge faces (1, ny+1) have no donor cell (cell 0 / ny+1 are
      ! off-array); they feed only re-wrapped ghosts, so zero them.
      ! Explicit do concurrent (not F(:,1,:)=0 array syntax -- offload).
      do concurrent(k=1:nz, i=1:nx)
         F(i, 1, k) = 0.0_wp
         F(i, ny + 1, k) = 0.0_wp
      end do
      ! Interior faces 2..ny: donors j-1 >= 1 (u>0) and j <= ny (u<0) in
      ! range.  (Was j=1:ny+1, reading cell 0 / ny+1 -- latent OOB.)
      do concurrent(k=1:nz, j=2:ny, i=1:nx) local(u, cfl, vol, conc)
         u = uhh(i, j, k)
         if (u > 0.0_wp) then
            vol = max(areaT(i, j - 1)*hprev(i, j - 1, k), DRAIN_MIN_VOL)
            cfl = min(u/vol, 1.0_wp)
            conc = aR(i, j - 1, k) - 0.5_wp*cfl* &
                   ((aR(i, j - 1, k) - aL(i, j - 1, k)) &
                    - a6(i, j - 1, k)*(1.0_wp - (2.0_wp/3.0_wp)*cfl))
            F(i, j, k) = u*conc
         else if (u < 0.0_wp) then
            vol = max(areaT(i, j)*hprev(i, j, k), DRAIN_MIN_VOL)
            cfl = min(-u/vol, 1.0_wp)
            conc = aL(i, j, k) + 0.5_wp*cfl* &
                   ((aR(i, j, k) - aL(i, j, k)) &
                    + a6(i, j, k)*(1.0_wp - (2.0_wp/3.0_wp)*cfl))
            F(i, j, k) = u*conc
         else
            F(i, j, k) = 0.0_wp
         end if
      end do
   end subroutine drain_swept_flux_y

   ! ====================================================================
   ! Q6 WENO drain helpers (reconstruction-ladder reuse).
   !
   ! The WENO path replaces the CW-PPM parabola (drain_parabola_* +
   ! drain_swept_flux_*) with a per-face swept-average WENO reconstruction
   ! evaluated straight from the concentration field Tr = hTr/hprev.  The
   ! swept-average convention is IDENTICAL to the CW-PPM path:
   !
   !   * the donor cell is (f-1) for u>0 and (f) for u<0 (f = the face
   !     index, uhh(f) sits between cells f-1 and f);
   !   * `cfl = min(|uhh|/(areaT·hprev_donor), 1)` is the donor's per-pass
   !     Courant number — this is the `sigma` the coastal helpers integrate
   !     the reconstruction over [1/2 - sigma, 1/2] of the donor's DOWNWIND
   !     edge;
   !   * the u<0 branch is mirror-exact — the stencil is gathered in
   !     downwind-positive order (step d = -1) so the same helper computes
   !     the donor's (physically LEFT) downwind edge.
   !
   ! Land / walls: interior land carries zero face transport (metric mask)
   ! so a land face contributes F=0 regardless of the reconstruction; land
   ! CELLS inside a live stencil are mirror-reflected (ppm_mirror_h, as the
   ! CW parabola does) so they never inject a spurious extremum.  Domain
   ! walls degrade through the rung ladder (recon_rung_for_face) exactly as
   ! the coastal kernel does — the availability counts are the coastal
   ! formulas, so the deepest stencil read stays in [1, nx]/[1, ny] for any
   ! nghost >= the configure-enforced per-rung minimum (weno5→3, 7→4, 9→5).
   ! ====================================================================

   pure subroutine drain_fill_conc(nx, ny, nz, hTr, hprev, tr)
      !! Concentration field Tr = hTr / max(hprev, DRAIN_MIN_H) for the WENO
      !! drain (the CW path builds the same field as the first loop of
      !! drain_parabola_*; factored out so the WENO path can reuse it
      !! without the parabola coefficients).
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: hTr(nx, ny, nz), hprev(nx, ny, nz)
      real(wp), intent(inout) :: tr(nx, ny, nz)
      integer :: i, j, k
      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         tr(i, j, k) = hTr(i, j, k)/max(hprev(i, j, k), DRAIN_MIN_H)
      end do
   end subroutine drain_fill_conc

   pure function weno_face_conc_x(nx, ny, nz, tr, wet_T, cc, jj, kk, d, &
                                  avail_up, avail_down, rung_max, cfl) result(conc)
      !! Swept-average WENO donor concentration at a zonal face.  `cc` is
      !! the donor cell (i-index), `d = +1` (u>0, downwind toward +i) or
      !! `d = -1` (u<0, downwind toward -i).  Gathers the mirrored/clamped
      !! stencil in downwind-positive order and dispatches to the coastal
      !! swept-average face helper at the highest feasible rung.
      !$acc routine seq
      integer, intent(in) :: nx, ny, nz, cc, jj, kk, d, avail_up, avail_down, rung_max
      real(wp), intent(in) :: tr(nx, ny, nz), wet_T(nx, ny)
      real(wp), intent(in) :: cfl
      real(wp) :: conc
      integer :: rung, ip, im
      real(wp) :: q0, qp1, qp2, qp3, qp4, qm1, qm2, qm3, qm4
      rung = recon_rung_for_face(avail_up, avail_down, rung_max)
      q0 = tr(cc, jj, kk)
      if (rung <= 0) then
         conc = q0                     ! donor cell (0th order)
         return
      end if
      ip = min(max(cc + d, 1), nx)
      im = min(max(cc - d, 1), nx)
      qp1 = ppm_mirror_h(tr(ip, jj, kk), q0, wet_T(ip, jj))
      qm1 = ppm_mirror_h(tr(im, jj, kk), q0, wet_T(im, jj))
      if (rung == 1) then
         conc = plm_face_swept(qm1, q0, qp1, cfl)
         return
      end if
      ip = min(max(cc + 2*d, 1), nx)
      im = min(max(cc - 2*d, 1), nx)
      qp2 = ppm_mirror_h(tr(ip, jj, kk), qp1, wet_T(ip, jj))
      qm2 = ppm_mirror_h(tr(im, jj, kk), qm1, wet_T(im, jj))
      if (rung == 2) then
         conc = weno5_face_swept(qm2, qm1, q0, qp1, qp2, cfl)
         return
      end if
      ip = min(max(cc + 3*d, 1), nx)
      im = min(max(cc - 3*d, 1), nx)
      qp3 = ppm_mirror_h(tr(ip, jj, kk), qp2, wet_T(ip, jj))
      qm3 = ppm_mirror_h(tr(im, jj, kk), qm2, wet_T(im, jj))
      if (rung == 3) then
         conc = weno7_face_swept(qm3, qm2, qm1, q0, qp1, qp2, qp3, cfl)
         return
      end if
      ip = min(max(cc + 4*d, 1), nx)
      im = min(max(cc - 4*d, 1), nx)
      qp4 = ppm_mirror_h(tr(ip, jj, kk), qp3, wet_T(ip, jj))
      qm4 = ppm_mirror_h(tr(im, jj, kk), qm3, wet_T(im, jj))
      conc = weno9_face_swept(qm4, qm3, qm2, qm1, q0, qp1, qp2, qp3, qp4, cfl)
   end function weno_face_conc_x

   pure function weno_face_conc_y(nx, ny, nz, tr, wet_T, ii, cc, kk, d, &
                                  avail_up, avail_down, rung_max, cfl) result(conc)
      !! Meridional analogue of weno_face_conc_x.  `cc` is the donor cell
      !! (j-index); the stencil steps along j with `d = +1` (u>0) / `-1` (u<0).
      !$acc routine seq
      integer, intent(in) :: nx, ny, nz, ii, cc, kk, d, avail_up, avail_down, rung_max
      real(wp), intent(in) :: tr(nx, ny, nz), wet_T(nx, ny)
      real(wp), intent(in) :: cfl
      real(wp) :: conc
      integer :: rung, jp, jm
      real(wp) :: q0, qp1, qp2, qp3, qp4, qm1, qm2, qm3, qm4
      rung = recon_rung_for_face(avail_up, avail_down, rung_max)
      q0 = tr(ii, cc, kk)
      if (rung <= 0) then
         conc = q0
         return
      end if
      jp = min(max(cc + d, 1), ny)
      jm = min(max(cc - d, 1), ny)
      qp1 = ppm_mirror_h(tr(ii, jp, kk), q0, wet_T(ii, jp))
      qm1 = ppm_mirror_h(tr(ii, jm, kk), q0, wet_T(ii, jm))
      if (rung == 1) then
         conc = plm_face_swept(qm1, q0, qp1, cfl)
         return
      end if
      jp = min(max(cc + 2*d, 1), ny)
      jm = min(max(cc - 2*d, 1), ny)
      qp2 = ppm_mirror_h(tr(ii, jp, kk), qp1, wet_T(ii, jp))
      qm2 = ppm_mirror_h(tr(ii, jm, kk), qm1, wet_T(ii, jm))
      if (rung == 2) then
         conc = weno5_face_swept(qm2, qm1, q0, qp1, qp2, cfl)
         return
      end if
      jp = min(max(cc + 3*d, 1), ny)
      jm = min(max(cc - 3*d, 1), ny)
      qp3 = ppm_mirror_h(tr(ii, jp, kk), qp2, wet_T(ii, jp))
      qm3 = ppm_mirror_h(tr(ii, jm, kk), qm2, wet_T(ii, jm))
      if (rung == 3) then
         conc = weno7_face_swept(qm3, qm2, qm1, q0, qp1, qp2, qp3, cfl)
         return
      end if
      jp = min(max(cc + 4*d, 1), ny)
      jm = min(max(cc - 4*d, 1), ny)
      qp4 = ppm_mirror_h(tr(ii, jp, kk), qp3, wet_T(ii, jp))
      qm4 = ppm_mirror_h(tr(ii, jm, kk), qm3, wet_T(ii, jm))
      conc = weno9_face_swept(qm4, qm3, qm2, qm1, q0, qp1, qp2, qp3, qp4, cfl)
   end function weno_face_conc_y

   pure subroutine drain_swept_flux_x_weno(nx, ny, nz, nghost, periodic, areaT, uhh, &
                                           hprev, tr, wet_T, recon, F)
      !! WENO analogue of drain_swept_flux_x: swept-average zonal face flux
      !! from the concentration field `tr` using the WENO rung ladder.
      !! `recon` is TRACER_RECON_WENO5/7/9 (1/2/3); the internal rung_max is
      !! recon + 1 (WENO5→2, WENO7→3, WENO9→4).
      !!
      !! `periodic`: on a periodic axis the ghost band is a wrapped copy of
      !! the interior, so the full-rung stencil is valid everywhere — force
      !! full rung (no position-based degradation).  This keeps the two
      !! images of a wrapped seam face bit-identical, so the interior flux
      !! telescopes and Σ(areaT·hTr) is conserved.  The per-rung nghost
      !! minimum (fail-loud at configure) guarantees the full-rung reads at
      !! the seam faces stay in [1, nx].  On a non-periodic (wall) axis the
      !! position ladder degrades toward the array edge; wall faces carry no
      !! transport, so a reduced-rung reconstruction there costs no accuracy.
      integer, intent(in) :: nx, ny, nz, nghost, recon
      logical, intent(in) :: periodic
      real(wp), intent(in) :: areaT(nx, ny)
      real(wp), intent(in) :: uhh(nx + 1, ny, nz)
      real(wp), intent(in) :: hprev(nx, ny, nz)
      real(wp), intent(in) :: tr(nx, ny, nz)
      real(wp), intent(in) :: wet_T(nx, ny)
      real(wp), intent(inout) :: F(nx + 1, ny, nz)
      integer :: i, j, k, rung_max, avail_up, avail_down
      real(wp) :: u, cfl, vol, conc
      rung_max = recon + 1
      do concurrent(k=1:nz, j=1:ny)
         F(1, j, k) = 0.0_wp
         F(nx + 1, j, k) = 0.0_wp
      end do
      do concurrent(k=1:nz, j=1:ny, i=2:nx) local(u, cfl, vol, conc, avail_up, avail_down)
         u = uhh(i, j, k)
         if (u > 0.0_wp) then
            ! donor = cell i-1, downwind = +i
            vol = max(areaT(i - 1, j)*hprev(i - 1, j, k), DRAIN_MIN_VOL)
            cfl = min(u/vol, 1.0_wp)
            if (periodic) then
               avail_up = nx
               avail_down = nx
            else
               avail_up = i - nghost
               avail_down = nx - nghost - i + 2
            end if
            conc = weno_face_conc_x(nx, ny, nz, tr, wet_T, i - 1, j, k, 1, &
                                    avail_up, avail_down, rung_max, cfl)
            F(i, j, k) = u*conc
         else if (u < 0.0_wp) then
            ! donor = cell i, downwind = -i
            vol = max(areaT(i, j)*hprev(i, j, k), DRAIN_MIN_VOL)
            cfl = min(-u/vol, 1.0_wp)
            if (periodic) then
               avail_up = nx
               avail_down = nx
            else
               avail_up = nx - nghost - i + 2
               avail_down = i - nghost
            end if
            conc = weno_face_conc_x(nx, ny, nz, tr, wet_T, i, j, k, -1, &
                                    avail_up, avail_down, rung_max, cfl)
            F(i, j, k) = u*conc
         else
            F(i, j, k) = 0.0_wp
         end if
      end do
   end subroutine drain_swept_flux_x_weno

   pure subroutine drain_swept_flux_y_weno(nx, ny, nz, nghost, periodic, areaT, uhh, &
                                           hprev, tr, wet_T, recon, F)
      !! Meridional analogue of drain_swept_flux_x_weno.
      integer, intent(in) :: nx, ny, nz, nghost, recon
      logical, intent(in) :: periodic
      real(wp), intent(in) :: areaT(nx, ny)
      real(wp), intent(in) :: uhh(nx, ny + 1, nz)
      real(wp), intent(in) :: hprev(nx, ny, nz)
      real(wp), intent(in) :: tr(nx, ny, nz)
      real(wp), intent(in) :: wet_T(nx, ny)
      real(wp), intent(inout) :: F(nx, ny + 1, nz)
      integer :: i, j, k, rung_max, avail_up, avail_down
      real(wp) :: u, cfl, vol, conc
      rung_max = recon + 1
      do concurrent(k=1:nz, i=1:nx)
         F(i, 1, k) = 0.0_wp
         F(i, ny + 1, k) = 0.0_wp
      end do
      do concurrent(k=1:nz, j=2:ny, i=1:nx) local(u, cfl, vol, conc, avail_up, avail_down)
         u = uhh(i, j, k)
         if (u > 0.0_wp) then
            vol = max(areaT(i, j - 1)*hprev(i, j - 1, k), DRAIN_MIN_VOL)
            cfl = min(u/vol, 1.0_wp)
            if (periodic) then
               avail_up = ny
               avail_down = ny
            else
               avail_up = j - nghost
               avail_down = ny - nghost - j + 2
            end if
            conc = weno_face_conc_y(nx, ny, nz, tr, wet_T, i, j - 1, k, 1, &
                                    avail_up, avail_down, rung_max, cfl)
            F(i, j, k) = u*conc
         else if (u < 0.0_wp) then
            vol = max(areaT(i, j)*hprev(i, j, k), DRAIN_MIN_VOL)
            cfl = min(-u/vol, 1.0_wp)
            if (periodic) then
               avail_up = ny
               avail_down = ny
            else
               avail_up = ny - nghost - j + 2
               avail_down = j - nghost
            end if
            conc = weno_face_conc_y(nx, ny, nz, tr, wet_T, i, j, k, -1, &
                                    avail_up, avail_down, rung_max, cfl)
            F(i, j, k) = u*conc
         else
            F(i, j, k) = 0.0_wp
         end if
      end do
   end subroutine drain_swept_flux_y_weno

   pure subroutine drain_update_tracer_x(nx, ny, nz, iareaT, F, hTr)
      !! hTr(i) -= (F(i+1) - F(i))·iareaT  (zonal flux divergence; dt is
      !! already baked into uhh⊂uhtr, so no dt here).
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: iareaT(nx, ny)
      real(wp), intent(in) :: F(nx + 1, ny, nz)
      real(wp), intent(inout) :: hTr(nx, ny, nz)
      integer :: i, j, k
      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         hTr(i, j, k) = hTr(i, j, k) - (F(i + 1, j, k) - F(i, j, k))*iareaT(i, j)
      end do
   end subroutine drain_update_tracer_x

   pure subroutine drain_update_tracer_y(nx, ny, nz, iareaT, F, hTr)
      !! Meridional analogue of drain_update_tracer_x.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: iareaT(nx, ny)
      real(wp), intent(in) :: F(nx, ny + 1, nz)
      real(wp), intent(inout) :: hTr(nx, ny, nz)
      integer :: i, j, k
      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         hTr(i, j, k) = hTr(i, j, k) - (F(i, j + 1, k) - F(i, j, k))*iareaT(i, j)
      end do
   end subroutine drain_update_tracer_y

   pure subroutine drain_update_tracer_x_budget(nx, ny, nz, iareaT, F, w, hTr, budget_adv)
      !! `drain_update_tracer_x` + the closed-budget fill: the SAME increment
      !! written to `hTr` is accumulated (times `w`) into `budget_adv`, so the
      !! console `out` term sees the horizontal tracer transport the windowed
      !! drain performs.  Without this the drain moves tracer that
      !! `ms%*_budget_horiz_adv` never records, and the Heat/Salt `Error`
      !! column has to fall back to raw drift (which then reports a live
      !! surface flux as a "leak").  `w` is
      !! `DRAIN_BUDGET_POST_AVERAGE_WEIGHT` for every drain call site.
      !!
      !! The `hTr` arithmetic is written EXACTLY as in the budget-free twin,
      !! so a tracer with a budget slot and one without still integrate
      !! bit-identically.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: iareaT(nx, ny)
      real(wp), intent(in) :: F(nx + 1, ny, nz)
      real(wp), intent(in) :: w
      real(wp), intent(inout) :: hTr(nx, ny, nz)
      real(wp), intent(inout) :: budget_adv(nx, ny, nz)
      integer :: i, j, k
      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         hTr(i, j, k) = hTr(i, j, k) - (F(i + 1, j, k) - F(i, j, k))*iareaT(i, j)
         budget_adv(i, j, k) = budget_adv(i, j, k) &
                               - w*(F(i + 1, j, k) - F(i, j, k))*iareaT(i, j)
      end do
   end subroutine drain_update_tracer_x_budget

   pure subroutine drain_update_tracer_y_budget(nx, ny, nz, iareaT, F, w, hTr, budget_adv)
      !! Meridional analogue of `drain_update_tracer_x_budget`.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: iareaT(nx, ny)
      real(wp), intent(in) :: F(nx, ny + 1, nz)
      real(wp), intent(in) :: w
      real(wp), intent(inout) :: hTr(nx, ny, nz)
      real(wp), intent(inout) :: budget_adv(nx, ny, nz)
      integer :: i, j, k
      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         hTr(i, j, k) = hTr(i, j, k) - (F(i, j + 1, k) - F(i, j, k))*iareaT(i, j)
         budget_adv(i, j, k) = budget_adv(i, j, k) &
                               - w*(F(i, j + 1, k) - F(i, j, k))*iareaT(i, j)
      end do
   end subroutine drain_update_tracer_y_budget

   pure subroutine drain_update_h_x(nx, ny, nz, iareaT, uhh, hprev)
      !! hprev(i) -= (uhh(i+1) - uhh(i))·iareaT  (volume div → thickness).
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: iareaT(nx, ny)
      real(wp), intent(in) :: uhh(nx + 1, ny, nz)
      real(wp), intent(inout) :: hprev(nx, ny, nz)
      integer :: i, j, k
      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         hprev(i, j, k) = hprev(i, j, k) - (uhh(i + 1, j, k) - uhh(i, j, k))*iareaT(i, j)
      end do
   end subroutine drain_update_h_x

   pure subroutine drain_update_h_y(nx, ny, nz, iareaT, uhh, hprev)
      !! Meridional analogue of drain_update_h_x.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: iareaT(nx, ny)
      real(wp), intent(in) :: uhh(nx, ny + 1, nz)
      real(wp), intent(inout) :: hprev(nx, ny, nz)
      integer :: i, j, k
      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         hprev(i, j, k) = hprev(i, j, k) - (uhh(i, j + 1, k) - uhh(i, j, k))*iareaT(i, j)
      end do
   end subroutine drain_update_h_y

   pure function continuity_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the continuity-PPM slot
      !! (0 when unallocated). One arr_bytes term per array — add a
      !! term here when a new allocatable joins the type.
      class(continuity_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = this%h_face_left_x%bytes() &
               + this%h_face_right_x%bytes() &
               + this%h_face_left_y%bytes() &
               + this%h_face_right_y%bytes() &
               + this%mt_h_new%bytes() &
               + this%mt_grounded%bytes() &
               + this%pd_theta%bytes() &
               + arr_bytes(this%uhtr) &
               + arr_bytes(this%vhtr) &
               + arr_bytes(this%hprev_work) &
               + arr_bytes(this%h_win_start) &
               + arr_bytes(this%uhr_x) &
               + arr_bytes(this%uhr_y) &
               + arr_bytes(this%uhh_x) &
               + arr_bytes(this%uhh_y) &
               + arr_bytes(this%tr_flux_x) &
               + arr_bytes(this%tr_flux_y) &
               + arr_bytes(this%tr_work) &
               + arr_bytes(this%pal) &
               + arr_bytes(this%par) &
               + arr_bytes(this%pa6)
   end function continuity_bytes

end module rdb_continuity
