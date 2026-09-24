!! Barotropic working-state slot for the split-explicit ocean driver.
module rdb_barotropic_workstate
   !! Transient barotropic state the split driver reads/writes between
   !! the slow baroclinic step and the fast barotropic substep loop.
   !! Transient (not restartable), unlike the prognostic
   !! `barotropic_state_t`.  C-grid stagger: scalars at centres,
   !! u at east faces, v at north faces, ζ at corners.  `F_slow_*`,
   !! `F_bt_*`, `ubt_at_n`, `vbt_at_n`, etc. are allocated only when
   !! `init` is passed `nz_ml` (barotropic-only unit tests skip them).
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_mem_report, only: arr_bytes
   implicit none
   private

   public :: barotropic_workstate_t
   public :: local_BT_cont_u_type, local_BT_cont_v_type
   ! Non-polymorphic device-attach entry points: owning slots call these
   ! directly with a `type(...)` actual so the OpenMP map base is the heap
   ! object, not a polymorphic stack box (AMD libomptarget cross-slot-overlap
   ! fix).  Type-bound `enter_data`/`exit_data` are select-type wrappers.
   public :: barotropic_workstate_enter_data_impl
   public :: barotropic_workstate_exit_data_impl

   type :: local_BT_cont_u_type
      !! Per-u-face coefficient pack for the piecewise-cubic
      !! barotropic-continuity flux closure (`find_uhbt` four-branch
      !! transport function consumes these).
      real(wp) :: FA_u_EE = 0.0_wp
         !! Marginal face area in the saturated far-east-draw
         !! regime `u < uBT_EE` [m].
      real(wp) :: FA_u_E0 = 0.0_wp
         !! Effective face area for near-zero u with u < 0 [m].
      real(wp) :: FA_u_W0 = 0.0_wp
         !! Effective face area for near-zero u with u > 0 [m].
      real(wp) :: FA_u_WW = 0.0_wp
         !! Marginal face area in the saturated far-west-draw
         !! regime `u > uBT_WW` [m].
      real(wp) :: uBT_WW = 0.0_wp
         !! Positive velocity threshold beyond which `find_uhbt`
         !! switches to the saturated linear branch [m/s].
         !! Must be ≥ 0.
      real(wp) :: uBT_EE = 0.0_wp
         !! Negative velocity threshold beyond which `find_uhbt`
         !! switches to the saturated linear branch [m/s].
         !! Must be ≤ 0.
      real(wp) :: uh_crvW = 0.0_wp
         !! Cubic correction in the near-zero positive branch
         !! providing C¹ continuity at `uBT_WW` [s²/m].
      real(wp) :: uh_crvE = 0.0_wp
         !! Cubic correction in the near-zero negative branch
         !! providing C¹ continuity at `uBT_EE` [s²/m].
      real(wp) :: uh_WW = 0.0_wp
         !! Mass transport at `u = uBT_WW` so the saturated
         !! branch matches the cubic [m²/s].
      real(wp) :: uh_EE = 0.0_wp
         !! Mass transport at `u = uBT_EE` so the saturated
         !! branch matches the cubic [m²/s].
   end type local_BT_cont_u_type

   type :: local_BT_cont_v_type
      !! Meridional mirror of `local_BT_cont_u_type` — north-draw
      !! and south-draw branches around v=0.
      real(wp) :: FA_v_NN = 0.0_wp
      real(wp) :: FA_v_N0 = 0.0_wp
      real(wp) :: FA_v_S0 = 0.0_wp
      real(wp) :: FA_v_SS = 0.0_wp
      real(wp) :: vBT_SS = 0.0_wp
         !! Positive (southward-draw) threshold [m/s]; must be ≥ 0.
      real(wp) :: vBT_NN = 0.0_wp
         !! Negative (northward-draw) threshold [m/s]; must be ≤ 0.
      real(wp) :: vh_crvS = 0.0_wp
      real(wp) :: vh_crvN = 0.0_wp
      real(wp) :: vh_SS = 0.0_wp
      real(wp) :: vh_NN = 0.0_wp
   end type local_BT_cont_v_type

   type :: barotropic_workstate_t
      logical :: is_init = .false.
         !! True between `init` and `destroy`.  Tracks GPU device
         !! attachment too — prefer to `allocated(...)` which only
         !! sees the host pointer.

      logical :: bt_substep_drag = .false.
         !! When `.true.`, driver fills `bt_rem_u/v` with the
         !! multiplicative BT-substep drag factor that damps `bt_ubt`
         !! / `bt_vbt` each inner step.  Default off ⇒ `bt_rem` ≡ 1,
         !! multiplication is a no-op (bit-identical).

      logical :: bt_correction_visc_rem = .false.
      logical :: bt_forcing_visc_rem = .false.
         !! MOM6 `wt_u` parity for the BT forcing assembly: weight the
         !! `F_bt_u/v` depth-mean (and the PGF-projection subtraction) by
         !! `h_face·visc_rem(k)` (`&ocean_bt_nml forcing_visc_rem`).
      logical :: bt_renorm_visc_rem = .false.
         !! MOM6 continuity-inversion parity (SPEC S2b): γ-weighted
         !! transport-matching renormaliser (`&ocean_bt_nml
         !! renorm_visc_rem`) — `visc_rem_u/v` forwarded into the slow
         !! continuity so `u_cor = u + du·γ_k` and the fluxes carry the
         !! same weights.
         !! When `.true.`, the h-weighted BT corrector multiplies its
         !! per-layer weight by `visc_rem_*(k)` — biasing the Δu
         !! distribution toward layers LESS damped by vertical
         !! viscosity.  Requires `bt_correction_h_weighted = .true.`.

      logical :: bt_correction_bc_pgf = .false.
         !! When `.true.`, `apply_bt_correction` adds the per-layer
         !! baroclinic-PGF retro-correction on top of the uniform /
         !! h-weighted Δu.  Requires `ocean_pgf_form = "fv_mom6"` and
         !! the `pbce` / `gtot_*` / `e_anom` / `eta_PF` fields filled.

      logical :: bt_correction_h_weighted = .false.
         !! When `.true.`, `apply_bt_correction` distributes the
         !! per-face Δu across layers proportional to
         !! `h_face(k) / ⟨h⟩_h` instead of uniformly.  Default off
         !! bit-identical.

      real(wp) :: g_bt = 9.81_wp
         !! Acceleration in the barotropic-substep η-gradient PGF
         !! (`a = -g_bt · ∇η`), m/s².  Defaults to full gravity
         !! (bit-identical for Mont/FV/coastal).  For
         !! `OPGF_VARIANT_GPRIME` the driver overrides with
         !! `pgf%gprime_gfs` so the BT mode runs at the reduced-gravity
         !! speed `sqrt(g_FS · H)`; else the gprime knobs net out to
         !! full `g`.

      real(wp) :: bebt = 0.1_wp
         !! MOM6 `BEBT` (default 0.1, as MOM6).  Continuity flux uses
         !! `ubt_trans = (1+bebt)·ubt^n − bebt·ubt^{n-1}` — the
         !! `BT_PROJECT_VELOCITY` spelling.  Because η is updated BEFORE
         !! the velocity in this loop, it is the same scheme as MOM6's
         !! default (`BT_PROJECT_VELOCITY = .false.`: predictor η, then
         !! `(1−bebt)·ubt^n + bebt·ubt^{n+1}` transport): this loop's η is
         !! MOM6's `eta_pred` and the velocity sequence is identical, so
         !! the per-substep damping `|λ|² = 1 − bebt·a²` and the stability
         !! limit `a ≤ 2/√(1+2·bebt)` are MOM6's.  `bebt = 0` ⇒
         !! `ubt_trans = ubt^n` (neutral forward-backward Euler).  Set via
         !! `&ocean_bt_nml bebt`.

      logical :: substep_zeta_ke = .true.
         !! Live `(ζ_bt+f)·v − ∇KE` in the fast loop (default,
         !! bit-identical).  `.false.` = MOM6-parity planetary-only
         !! substeps; the `subtract_fast_cor_ref` reference reduces to
         !! `f·v̄` to match.  See `&ocean_bt_nml substep_zeta_ke`.

      ! ---- Primary barotropic-substep state ----
      real(wp), allocatable :: bt_eta(:, :)
         !! Barotropic SSH at cell centres, shape (nx, ny).  Defined
         !! as `sum_k(h_layer) - bt_H_ref`.
      real(wp), allocatable :: bt_ubt(:, :)
         !! Depth-mean u at east faces, shape (nx+1, ny).
      real(wp), allocatable :: bt_vbt(:, :)
         !! Depth-mean v at north faces, shape (nx, ny+1).
      real(wp), allocatable :: bt_H_ref(:, :)
         !! Reference column thickness (m) at cell centres, shape
         !! (nx, ny).  Constant for Eulerian-z; set by the driver
         !! from the initial multilayer state.  Used as `H` in the
         !! linearized continuity `∂η/∂t = -div(H · u_bt)`.

      ! ---- Fast-mode time-averaging accumulators ----
      real(wp), allocatable :: ubt_sum(:, :)
      real(wp), allocatable :: vbt_sum(:, :)
      real(wp), allocatable :: eta_sum(:, :)

      ! ---- Time-mean depth-integrated transport ----
      ! Accumulated each fast substep as `(H_ref + η_inst) · u_bt_inst`
      ! at faces, time-averaged at loop end; consumed by slow continuity
      ! as a barotropic constraint so per-layer mass fluxes vertically
      ! sum to the substep transport (makes `sum_k(h_layer) = H_ref +
      ! η_end`, so `apply_bt_correction`'s h-rescale a no-op).
      real(wp), allocatable :: uhbt_sum(:, :)
         !! East-face transport accumulator, shape `(nx+1, ny)`.
      real(wp), allocatable :: vhbt_sum(:, :)
         !! North-face transport accumulator, shape `(nx, ny+1)`.
      real(wp), allocatable :: bt_uhbt(:, :)
         !! Time-mean east-face transport (m²/s), shape `(nx+1, ny)`.
         !! `uhbt_sum / n_inner`.
      real(wp), allocatable :: bt_vhbt(:, :)
         !! Time-mean north-face transport (m²/s), shape `(nx, ny+1)`.

      ! ---- End-of-barotropic-substep snapshot ----
      ! Saved BEFORE the time-mean overwrite at fast-loop end.  Used by
      ! `apply_bt_correction` so the recombined per-layer momentum AND
      ! the h_layer rescale both see the end-of-step barotropic mode
      ! (Hallberg 2009 split-explicit convention); mixing end-step
      ! velocity with time-mean SSH corrupts gravity-wave dispersion.
      real(wp), allocatable :: bt_ubt_end(:, :)
      real(wp), allocatable :: bt_vbt_end(:, :)
      real(wp), allocatable :: bt_eta_end(:, :)

      ! ---- Nonlinear barotropic-substep scratch ----
      real(wp), allocatable :: bt_zeta_corner(:, :)
         !! Barotropic relative vorticity at corners, shape (nx+1, ny+1).
      real(wp), allocatable :: bt_ke_centre(:, :)
         !! Barotropic kinetic energy at cell centres, shape (nx, ny).
      real(wp), allocatable :: bt_eta_new(:, :)
         !! Per-substep η^{n+1} scratch, shape (nx, ny).  Jacobi
         !! buffer so the face-thickness reads in Pass 1 stay race-
         !! free under the Fortran 2018 `do concurrent` semantics.

      ! ---- Coriolis/advection reference velocity (MOM6 `ubt_Cor`) ----
      ! The barotropic velocity at which `subtract_fast_cor_ref`
      ! evaluates the reference `(ζ+f)·v − ∇KE` it removes from the
      ! substep forcing.  It MUST be the depth mean of the SAME layer
      ! velocity the slow `cor%pv_flux_*` in `F_bt` was evaluated on,
      ! or the uncancelled residual `f × (v̄_ref − v̄_slow)` is injected
      ! into every barotropic substep as a near-constant forcing and
      ! pumps the basin's gravest Poincaré seiche exponentially.
      !   * `ssp_rk2` — slow tendencies are evaluated on the prognostic
      !     `u^n`, so this is a copy of the stage-entry `bt_ubt/bt_vbt`
      !     (bit-identical to reading `bt_ubt/bt_vbt` directly).
      !   * `pred_corr` — slow tendencies are evaluated on `u_av/v_av`,
      !     so this is the depth mean of `u_av/v_av` under the SAME
      !     weights the forcing depth-mean used (h, or h·visc_rem when
      !     `&ocean_bt_nml forcing_visc_rem`).  MOM6 `MOM_barotropic`
      !     builds `Cor_ref_u` from `ubt_Cor = Σ_k wt_u·U_Cor` with
      !     `U_Cor = u_av`, i.e. the same velocity `CorAdCalc` used.
      real(wp), allocatable :: cor_ref_u(:, :)
         !! u-face Coriolis/advection reference velocity, shape (nx+1, ny).
      real(wp), allocatable :: cor_ref_v(:, :)
         !! v-face counterpart, shape (nx, ny+1).

      ! ---- BEBT projection: previous-substep velocity snapshots ----
      ! Hold u^{n-1} / v^{n-1} between substeps for the η-update
      ! extrapolation `(1+bebt)·u^n − bebt·u^{n-1}`.  Initialised to
      ! `bt_ubt / bt_vbt` at the top of each outer step so the first
      ! substep has zero extrapolation (bit-identical when `bebt = 0`).
      real(wp), allocatable :: bt_ubt_prev(:, :)
         !! u-face velocity from the previous substep, shape (nx+1, ny).
      real(wp), allocatable :: bt_vbt_prev(:, :)
         !! v-face counterpart, shape (nx, ny+1).

      ! ---- Split-driver slow-tendency accumulators ----
      ! Only allocated when `init` is called with the optional `nz_ml`.
      real(wp), allocatable :: F_slow_u(:, :, :)
         !! Per-face slow u-acceleration sum, shape (nx+1, ny, nz_ml).
      real(wp), allocatable :: F_slow_v(:, :, :)
         !! Per-face slow v-acceleration sum, shape (nx, ny+1, nz_ml).
      real(wp), allocatable :: F_bt_u(:, :)
         !! Depth-mean of `F_slow_u`, shape (nx+1, ny).  Used in the
         !! `apply_bt_correction` subtraction to clean the slow-apply's
         !! bt projection out of every layer before installing the
         !! barotropic-substep bt mode.
      real(wp), allocatable :: F_bt_v(:, :)
         !! Depth-mean of `F_slow_v`, shape (nx, ny+1).
      real(wp), allocatable :: F_bt_u_fast(:, :)
         !! `F_bt_u` minus the bt projection of the PGF, shape
         !! (nx+1, ny).  Passed as the substep's `force_u` so the
         !! substep's own `-G·∂η/∂x` is the only bt PGF on the bt mode
         !! (else the slow + internal PGF stack to `-2G·∂η/∂x`,
         !! doubling the effective gravity-wave speed).
      real(wp), allocatable :: F_bt_v_fast(:, :)
         !! v counterpart, shape (nx, ny+1).
      real(wp), allocatable :: ubt_at_n(:, :)
         !! Depth-mean u at the start of the outer step, shape
         !! (nx+1, ny).  Used in the recombine step:
         !!   u^{n+1}(k) = u^*(k) + (⟨u_bt⟩ - ubt_at_n - dt·F_bt_u)
      real(wp), allocatable :: vbt_at_n(:, :)
         !! v counterpart, shape (nx, ny+1).

      ! ---- bc-PGF per-layer correction ----
      ! Per-layer baroclinic-PGF retro-correction for the η change
      ! during the BT substep: the slow PGF used `eta_PF`, but the BT
      ! substep evolves η to `bt_eta_end`, so each layer would feel a
      ! PGF based on stale `eta_PF`.  Adds back the layer-dependent
      ! response, scaled by `pbce(k)` minus column-mean `gtot_face`.
      ! Fills only when `ocean_bt_correction_bc_pgf = .true.`; else
      ! holds init zeros and the corrector skips it (bit-identical).
      real(wp), allocatable :: pbce(:, :, :)
         !! Per-layer pressure-anomaly gravity coefficient (m/s²),
         !! shape (nx, ny, nz_ml).  Column-mean equals `gtot_face`.
         !! Built by `compute_pbce`; requires the FV_MOM6 PGF.
      real(wp), allocatable :: gtot_E(:, :)
         !! Depth-weighted column average of `pbce` evaluated at the
         !! east face of cell `(i, j)`, shape (nx, ny).  Built by
         !! `compute_gtot_faces`.
      real(wp), allocatable :: gtot_W(:, :)
         !! West-face counterpart, shape (nx, ny).
      real(wp), allocatable :: gtot_N(:, :)
         !! North-face counterpart, shape (nx, ny).
      real(wp), allocatable :: gtot_S(:, :)
         !! South-face counterpart, shape (nx, ny).
      real(wp), allocatable :: e_anom(:, :)
         !! SSH anomaly relative to `eta_PF`, shape (nx, ny).
         !! `e_anom = 0.5·(bt_eta_end + bt_eta) − eta_PF`.
      real(wp), allocatable :: eta_PF(:, :)
         !! Snapshot of `bt_eta` taken just before the slow PGF is
         !! computed (= η the PGF "saw").  Updated each RK2 stage.
         !! Shape (nx, ny).

      ! ---- BT corrector visc_rem weights ----
      ! Per-face per-layer fraction of velocity remaining after viscous
      ! damping over one outer step; biases the h-weighted corrector's
      ! per-layer Δu toward LESS-damped layers.  Default 1.0 ⇒ h-only
      ! path bit-identically.
      real(wp), allocatable :: visc_rem_u(:, :, :)
         !! u-face per-layer visc_rem, shape (nx+1, ny, nz_ml).
      real(wp), allocatable :: visc_rem_v(:, :, :)
         !! v-face per-layer visc_rem, shape (nx, ny+1, nz_ml).

      ! ---- bt-substep multiplicative drag damping ----
      ! Per-face damping factor in [0,1], applied inside the BT substep
      ! loop after each u/v update.  Computed once per stage from the
      ! linear-drag coefficient + face Htot:
      !   bt_rem_face = Htot / (Htot + r · HBBL · dt_inner)
      ! Cumulative over substeps ⇒ exp(-r·HBBL·t/Htot)-style damping.
      ! `bt_substep_drag` off ⇒ arrays hold 1, no-op (bit-identical).
      real(wp), allocatable :: bt_rem_u(:, :)
         !! u-face damping factor, shape (nx+1, ny).
      real(wp), allocatable :: bt_rem_v(:, :)
         !! v-face counterpart, shape (nx, ny+1).

      ! ---- Barotropic linear (Rayleigh) wave drag (Egbert & Ray 2001;
      ! Jayne & St Laurent 2001) ----
      ! Static per-face piston velocity `r_H` [m/s] representing the
      ! barotropic-to-internal-tide energy sink.  MULTIPLIED into
      ! `bt_rem_u/v` by `compute_bt_rem_wave_drag` exactly as MOM6
      ! composes `lin_drag_u` with the viscous remnant.
      ! `lwd_enable = .false.` (default)
      ! ⇒ arrays stay unallocated and every BT path is bit-identical —
      ! same "arrays stay unallocated" contract as `use_bt_cont_type`
      ! below.  `lwd_` (not `wd_`) because `wd_` is taken by wet/dry.
      logical :: lwd_enable = .false.
         !! Driver writes from `&ocean_bt_nml wave_drag`.
      real(wp), allocatable :: lwd_drag_u(:, :)
         !! u-face piston velocity `r_H` [m/s], shape (nx+1, ny).
         !! Static after `configure_ocean_wave_drag`; allocated only
         !! when `lwd_enable = .true.`.
      real(wp), allocatable :: lwd_drag_v(:, :)
         !! v-face counterpart, shape (nx, ny+1).

      ! ---- BT_cont_type flux-bounded continuity ----
      ! Per-face piecewise-cubic flux closure replacing the naive
      ! `uh = u · h_face` (consumed by `find_uhbt`).  When
      ! `use_bt_cont_type = .false.` the arrays stay unallocated and
      ! every BT path is bit-identical.
      logical :: use_bt_cont_type = .false.
         !! Driver writes from the `ocean_use_bt_cont_type` namelist.
      type(local_BT_cont_u_type), allocatable :: BTCL_u(:, :)
         !! Per-u-face flux-closure coefficients, shape (nx+1, ny).
         !! Allocated only when `use_bt_cont_type = .true.`.
      type(local_BT_cont_v_type), allocatable :: BTCL_v(:, :)
         !! Per-v-face counterpart, shape (nx, ny+1).

      ! ---- Upstream-PPM face thickness for the BT chain ----
      ! Per-face column-sum of `h_layer` at the upstream face side,
      ! built once per outer step from `u/v_face_*_layer` signs and
      ! consumed across the BT chain (`derive_bt_from_layers`, substep,
      ! `apply_bt_correction`) so it shares the per-layer PPM
      ! face-thickness convention — kills phantom bed-layer velocity at
      ! slopes.  `use_upstream_h_face = .false.` ⇒ unallocated, legacy
      ! centred-h (bit-identical).
      logical :: use_upstream_h_face = .false.
         !! Driver writes from the `ocean_bt_upstream_h_face` namelist.
      real(wp), allocatable :: h_face_up_x(:, :)
         !! East-face upstream-h column sum (m), shape (nx+1, ny).
         !! `Σ_k h_layer(upstream_cell, j, k)` where the per-layer
         !! upstream selection follows `u_face_x_layer` sign.
         !! Allocated only when `use_upstream_h_face = .true.`.
      real(wp), allocatable :: h_face_up_y(:, :)
         !! North-face counterpart, shape (nx, ny+1).  Built from
         !! `v_face_y_layer` sign.

      ! ---- Dynamic wetting/drying (docs/ocean_wetdry_plan.md) ----
      ! Knobs + workspaces for the positive-definite wet/dry BT substep
      ! branch.  `wetdry_enable = .false.` (default) leaves every wd_*
      ! array unallocated and the substep on its unmodified centred-face
      ! path — byte-identical.  Allocated by `configure_ocean_wetdry`
      ! (rdb_ocean_setup) when `&ocean_wetdry_nml enable` is on.
      logical :: wetdry_enable = .false.
         !! Driver writes from `&ocean_wetdry_nml enable`.
      real(wp) :: wd_dry_depth = 0.05_wp
         !! Total-depth dry threshold (m); cells with `D < wd_dry_depth`
         !! are dynamically dry.  From `&ocean_wetdry_nml dry_depth`.
      real(wp) :: wd_rewet_depth = 0.10_wp
         !! Hysteresis re-wet threshold (m), > `wd_dry_depth`.  From
         !! `&ocean_wetdry_nml rewet_depth`.
      real(wp), allocatable :: wd_wet_dyn(:, :)
         !! Dynamic cell wet mask (1 = wet, 0 = dry), shape (nx, ny).
         !! PERSISTENT hysteresis state: updated every BT substep from
         !! `D = bt_H_ref + bt_eta` (wet above `wd_rewet_depth`, dry
         !! below `wd_dry_depth`, held in between).  Seeded from the
         !! initial D at configure.  Composes multiplicatively ON TOP of
         !! the static land masks (`metrics%wet_u/v` + zeroed metrics) —
         !! a static-land face can never be dynamically opened.
      real(wp), allocatable :: wd_theta(:, :)
         !! Per-cell positive-definite outflow limiter factor in [0, 1],
         !! shape (nx, ny).  `theta = min(1, available_volume /
         !! substep_outflow_volume)`; each face flux is scaled by
         !! `min(theta_L, theta_R)`, which guarantees `D >= 0` every
         !! substep (a cell drains at most what it holds).  Deep water
         !! ⇒ theta ≡ 1 ⇒ the limiter is exactly inert.
      real(wp), allocatable :: wd_flux_x(:, :)
         !! East-face provisional-then-limited volume flux (m³/s), shape
         !! (nx+1, ny).  Pass A fills `h_up · ubt_trans · dy_cu` with
         !! UPWIND face thickness + the FROUDE_CAP thin-face velocity
         !! guard; Pass C scales it by `min(theta_L, theta_R)` in place.
      real(wp), allocatable :: wd_flux_y(:, :)
         !! North-face counterpart, shape (nx, ny+1).
      real(wp), allocatable :: wd_open_u(:, :)
         !! Dynamic u-face open mask (1 = open, 0 = blocked), shape
         !! (nx+1, ny).  Bed-blocking (C-grid analogue of the coastal
         !! hydrostatic reconstruction): a face into a dry cell is a
         !! wall unless the wet side's surface stands above the dry
         !! side's bed elevation + dry_depth.  Rewetting starts from
         !! rest at the face.  Consumed by the substep Pass 2 gate and
         !! by the driver's layer-velocity masking.
      real(wp), allocatable :: wd_open_v(:, :)
         !! v-face counterpart, shape (nx, ny+1).
   contains
      procedure, non_overridable :: init => barotropic_workstate_init
      procedure, non_overridable :: destroy => barotropic_workstate_destroy
      procedure, non_overridable :: enter_data => barotropic_workstate_enter_data
      procedure, non_overridable :: exit_data => barotropic_workstate_exit_data
      procedure, non_overridable :: bytes => barotropic_workstate_bytes
   end type barotropic_workstate_t

contains

   subroutine barotropic_workstate_init(this, grid, nz_ml)
      !! Allocate the 2D barotropic-substep arrays.  Pass `nz_ml` to also
      !! allocate the split-driver slow-tendency accumulators; omit
      !! when only the 2D barotropic substep is needed (unit tests).
      class(barotropic_workstate_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      integer, intent(in), optional :: nz_ml
      integer :: nx, ny, nz

      nx = grid%nx_total
      ny = grid%ny_total

      allocate (this%bt_eta(nx, ny), source=0.0_wp)
      allocate (this%bt_H_ref(nx, ny), source=0.0_wp)
      allocate (this%bt_ubt(nx + 1, ny), source=0.0_wp)
      allocate (this%bt_vbt(nx, ny + 1), source=0.0_wp)
      allocate (this%eta_sum(nx, ny), source=0.0_wp)
      allocate (this%ubt_sum(nx + 1, ny), source=0.0_wp)
      allocate (this%vbt_sum(nx, ny + 1), source=0.0_wp)
      allocate (this%uhbt_sum(nx + 1, ny), source=0.0_wp)
      allocate (this%vhbt_sum(nx, ny + 1), source=0.0_wp)
      allocate (this%bt_uhbt(nx + 1, ny), source=0.0_wp)
      allocate (this%bt_vhbt(nx, ny + 1), source=0.0_wp)
      allocate (this%bt_ubt_end(nx + 1, ny), source=0.0_wp)
      allocate (this%bt_vbt_end(nx, ny + 1), source=0.0_wp)
      allocate (this%bt_eta_end(nx, ny), source=0.0_wp)
      allocate (this%bt_zeta_corner(nx + 1, ny + 1), source=0.0_wp)
      allocate (this%bt_ke_centre(nx, ny), source=0.0_wp)
      allocate (this%bt_eta_new(nx, ny), source=0.0_wp)
      allocate (this%cor_ref_u(nx + 1, ny), source=0.0_wp)
      allocate (this%cor_ref_v(nx, ny + 1), source=0.0_wp)
      allocate (this%bt_ubt_prev(nx + 1, ny), source=0.0_wp)
      allocate (this%bt_vbt_prev(nx, ny + 1), source=0.0_wp)

      if (present(nz_ml)) then
         nz = nz_ml
         allocate (this%F_slow_u(nx + 1, ny, nz), source=0.0_wp)
         allocate (this%F_slow_v(nx, ny + 1, nz), source=0.0_wp)
         allocate (this%F_bt_u(nx + 1, ny), source=0.0_wp)
         allocate (this%F_bt_v(nx, ny + 1), source=0.0_wp)
         allocate (this%F_bt_u_fast(nx + 1, ny), source=0.0_wp)
         allocate (this%F_bt_v_fast(nx, ny + 1), source=0.0_wp)
         allocate (this%ubt_at_n(nx + 1, ny), source=0.0_wp)
         allocate (this%vbt_at_n(nx, ny + 1), source=0.0_wp)
         allocate (this%pbce(nx, ny, nz), source=0.0_wp)
         allocate (this%gtot_E(nx, ny), source=0.0_wp)
         allocate (this%gtot_W(nx, ny), source=0.0_wp)
         allocate (this%gtot_N(nx, ny), source=0.0_wp)
         allocate (this%gtot_S(nx, ny), source=0.0_wp)
         allocate (this%e_anom(nx, ny), source=0.0_wp)
         allocate (this%eta_PF(nx, ny), source=0.0_wp)
         allocate (this%visc_rem_u(nx + 1, ny, nz), source=1.0_wp)
         allocate (this%visc_rem_v(nx, ny + 1, nz), source=1.0_wp)
      end if

      ! bt_rem_u/v are referenced unconditionally (no-op when ≡ 1);
      ! allocate outside the nz_ml block so barotropic-only unit tests
      ! still get them.
      allocate (this%bt_rem_u(nx + 1, ny), source=1.0_wp)
      allocate (this%bt_rem_v(nx, ny + 1), source=1.0_wp)

      ! BT_cont_type coefficient packs are allocated lazily by the
      ! driver after the namelist toggle is read.

      this%is_init = .true.
   end subroutine barotropic_workstate_init

   subroutine barotropic_workstate_destroy(this)
      class(barotropic_workstate_t), intent(inout) :: this
      this%is_init = .false.
      if (allocated(this%bt_eta)) deallocate (this%bt_eta)
      if (allocated(this%bt_H_ref)) deallocate (this%bt_H_ref)
      if (allocated(this%bt_ubt)) deallocate (this%bt_ubt)
      if (allocated(this%bt_vbt)) deallocate (this%bt_vbt)
      if (allocated(this%ubt_sum)) deallocate (this%ubt_sum)
      if (allocated(this%vbt_sum)) deallocate (this%vbt_sum)
      if (allocated(this%eta_sum)) deallocate (this%eta_sum)
      if (allocated(this%uhbt_sum)) deallocate (this%uhbt_sum)
      if (allocated(this%vhbt_sum)) deallocate (this%vhbt_sum)
      if (allocated(this%bt_uhbt)) deallocate (this%bt_uhbt)
      if (allocated(this%bt_vhbt)) deallocate (this%bt_vhbt)
      if (allocated(this%bt_ubt_end)) deallocate (this%bt_ubt_end)
      if (allocated(this%bt_vbt_end)) deallocate (this%bt_vbt_end)
      if (allocated(this%bt_eta_end)) deallocate (this%bt_eta_end)
      if (allocated(this%bt_zeta_corner)) deallocate (this%bt_zeta_corner)
      if (allocated(this%bt_ke_centre)) deallocate (this%bt_ke_centre)
      if (allocated(this%bt_eta_new)) deallocate (this%bt_eta_new)
      if (allocated(this%cor_ref_u)) deallocate (this%cor_ref_u)
      if (allocated(this%cor_ref_v)) deallocate (this%cor_ref_v)
      if (allocated(this%bt_ubt_prev)) deallocate (this%bt_ubt_prev)
      if (allocated(this%bt_vbt_prev)) deallocate (this%bt_vbt_prev)
      if (allocated(this%F_slow_u)) deallocate (this%F_slow_u)
      if (allocated(this%F_slow_v)) deallocate (this%F_slow_v)
      if (allocated(this%F_bt_u)) deallocate (this%F_bt_u)
      if (allocated(this%F_bt_v)) deallocate (this%F_bt_v)
      if (allocated(this%F_bt_u_fast)) deallocate (this%F_bt_u_fast)
      if (allocated(this%F_bt_v_fast)) deallocate (this%F_bt_v_fast)
      if (allocated(this%ubt_at_n)) deallocate (this%ubt_at_n)
      if (allocated(this%vbt_at_n)) deallocate (this%vbt_at_n)
      if (allocated(this%pbce)) deallocate (this%pbce)
      if (allocated(this%gtot_E)) deallocate (this%gtot_E)
      if (allocated(this%gtot_W)) deallocate (this%gtot_W)
      if (allocated(this%gtot_N)) deallocate (this%gtot_N)
      if (allocated(this%gtot_S)) deallocate (this%gtot_S)
      if (allocated(this%e_anom)) deallocate (this%e_anom)
      if (allocated(this%eta_PF)) deallocate (this%eta_PF)
      if (allocated(this%visc_rem_u)) deallocate (this%visc_rem_u)
      if (allocated(this%visc_rem_v)) deallocate (this%visc_rem_v)
      if (allocated(this%bt_rem_u)) deallocate (this%bt_rem_u)
      if (allocated(this%bt_rem_v)) deallocate (this%bt_rem_v)
      if (allocated(this%lwd_drag_u)) deallocate (this%lwd_drag_u)
      if (allocated(this%lwd_drag_v)) deallocate (this%lwd_drag_v)
      if (allocated(this%BTCL_u)) deallocate (this%BTCL_u)
      if (allocated(this%BTCL_v)) deallocate (this%BTCL_v)
      if (allocated(this%h_face_up_x)) deallocate (this%h_face_up_x)
      if (allocated(this%h_face_up_y)) deallocate (this%h_face_up_y)
      if (allocated(this%wd_wet_dyn)) deallocate (this%wd_wet_dyn)
      if (allocated(this%wd_theta)) deallocate (this%wd_theta)
      if (allocated(this%wd_flux_x)) deallocate (this%wd_flux_x)
      if (allocated(this%wd_flux_y)) deallocate (this%wd_flux_y)
      if (allocated(this%wd_open_u)) deallocate (this%wd_open_u)
      if (allocated(this%wd_open_v)) deallocate (this%wd_open_v)
   end subroutine barotropic_workstate_destroy

   subroutine barotropic_workstate_enter_data(this)
      !! Attaches component arrays only — no bare `copyin(this)` (the
      !! polymorphic stack descriptor caused AMD libomptarget cross-slot
      !! overlap).  The workstate descriptor reaches the device via the
      !! root `copyin(state%ocean)` (bt_work is inline all the way up);
      !! the per-slot copyin was a competing second mapping that cost
      !! ~36% of the fast loop.
      class(barotropic_workstate_t), intent(inout) :: this
      select type (this)
      type is (barotropic_workstate_t)
         call barotropic_workstate_enter_data_impl(this)
      end select
   end subroutine barotropic_workstate_enter_data

   subroutine barotropic_workstate_enter_data_impl(this)
      type(barotropic_workstate_t), intent(inout) :: this
      !$acc enter data copyin(this%bt_eta, this%bt_H_ref)
      !$acc enter data copyin(this%bt_ubt, this%bt_vbt)
      !$acc enter data copyin(this%eta_sum, this%ubt_sum, this%vbt_sum)
      !$acc enter data copyin(this%uhbt_sum, this%vhbt_sum, this%bt_uhbt, this%bt_vhbt)
      !$acc enter data copyin(this%bt_ubt_end, this%bt_vbt_end, this%bt_eta_end)
      !$acc enter data copyin(this%bt_zeta_corner, this%bt_ke_centre, this%bt_eta_new)
      !$acc enter data copyin(this%cor_ref_u, this%cor_ref_v)
      !$acc enter data copyin(this%bt_ubt_prev, this%bt_vbt_prev)
      ! bt_rem_u/v: always present (barotropic-only path uses them too)
      !$acc enter data copyin(this%bt_rem_u, this%bt_rem_v)
      ! Wave-drag piston-velocity maps: filled on the host at configure
      ! time and never written on the device, so this MUST be `copyin`
      ! (not `create`) — see CLAUDE.md gotcha (2).  Lazy: allocated only
      ! when `lwd_enable`.
      if (allocated(this%lwd_drag_u)) then
         !$acc enter data copyin(this%lwd_drag_u, this%lwd_drag_v)
      end if
      if (allocated(this%F_slow_u)) then
         !$acc enter data copyin(this%F_slow_u, this%F_slow_v)
         !$acc enter data copyin(this%F_bt_u, this%F_bt_v)
         !$acc enter data copyin(this%F_bt_u_fast, this%F_bt_v_fast)
         !$acc enter data copyin(this%ubt_at_n, this%vbt_at_n)
         !$acc enter data copyin(this%pbce)
         !$acc enter data copyin(this%gtot_E, this%gtot_W, this%gtot_N, this%gtot_S)
         !$acc enter data copyin(this%e_anom, this%eta_PF)
         !$acc enter data copyin(this%visc_rem_u, this%visc_rem_v)
      end if
      ! BTCL_u/v are arrays of derived type with POD scalar components,
      ! so copying the array body is enough (parent before components
      ! on enter, reverse on exit).
      if (allocated(this%BTCL_u)) then
         !$acc enter data copyin(this%BTCL_u)
      end if
      if (allocated(this%BTCL_v)) then
         !$acc enter data copyin(this%BTCL_v)
      end if
      ! Upstream-h-face slots — same lazy pattern as BTCL_u/v.
      if (allocated(this%h_face_up_x)) then
         !$acc enter data copyin(this%h_face_up_x)
      end if
      if (allocated(this%h_face_up_y)) then
         !$acc enter data copyin(this%h_face_up_y)
      end if
      ! Wet/dry workspaces — same lazy pattern (allocated only when
      ! `wetdry_enable`); missing this attach = the per-launch memcpy
      ! explosion foot-gun, so every wd_* array is listed.
      if (allocated(this%wd_wet_dyn)) then
         !$acc enter data copyin(this%wd_wet_dyn, this%wd_theta)
         !$acc enter data copyin(this%wd_flux_x, this%wd_flux_y)
         !$acc enter data copyin(this%wd_open_u, this%wd_open_v)
      end if
   end subroutine barotropic_workstate_enter_data_impl

   subroutine barotropic_workstate_exit_data(this)
      class(barotropic_workstate_t), intent(inout) :: this
      select type (this)
      type is (barotropic_workstate_t)
         call barotropic_workstate_exit_data_impl(this)
      end select
   end subroutine barotropic_workstate_exit_data

   subroutine barotropic_workstate_exit_data_impl(this)
      type(barotropic_workstate_t), intent(inout) :: this
      !$acc exit data delete(this%bt_eta, this%bt_H_ref)
      !$acc exit data delete(this%bt_ubt, this%bt_vbt)
      !$acc exit data delete(this%eta_sum, this%ubt_sum, this%vbt_sum)
      !$acc exit data delete(this%uhbt_sum, this%vhbt_sum, this%bt_uhbt, this%bt_vhbt)
      !$acc exit data delete(this%bt_ubt_end, this%bt_vbt_end, this%bt_eta_end)
      !$acc exit data delete(this%cor_ref_u, this%cor_ref_v)
      !$acc exit data delete(this%bt_zeta_corner, this%bt_ke_centre, this%bt_eta_new)
      !$acc exit data delete(this%bt_ubt_prev, this%bt_vbt_prev)
      !$acc exit data delete(this%bt_rem_u, this%bt_rem_v)
      if (allocated(this%lwd_drag_u)) then
         !$acc exit data delete(this%lwd_drag_u, this%lwd_drag_v)
      end if
      if (allocated(this%F_slow_u)) then
         !$acc exit data delete(this%F_slow_u, this%F_slow_v)
         !$acc exit data delete(this%F_bt_u, this%F_bt_v)
         !$acc exit data delete(this%F_bt_u_fast, this%F_bt_v_fast)
         !$acc exit data delete(this%ubt_at_n, this%vbt_at_n)
         !$acc exit data delete(this%pbce)
         !$acc exit data delete(this%gtot_E, this%gtot_W, this%gtot_N, this%gtot_S)
         !$acc exit data delete(this%e_anom, this%eta_PF)
         !$acc exit data delete(this%visc_rem_u, this%visc_rem_v)
      end if
      if (allocated(this%BTCL_u)) then
         !$acc exit data delete(this%BTCL_u)
      end if
      if (allocated(this%BTCL_v)) then
         !$acc exit data delete(this%BTCL_v)
      end if
      if (allocated(this%h_face_up_x)) then
         !$acc exit data delete(this%h_face_up_x)
      end if
      if (allocated(this%h_face_up_y)) then
         !$acc exit data delete(this%h_face_up_y)
      end if
      if (allocated(this%wd_wet_dyn)) then
         !$acc exit data delete(this%wd_open_u, this%wd_open_v)
         !$acc exit data delete(this%wd_flux_x, this%wd_flux_y)
         !$acc exit data delete(this%wd_wet_dyn, this%wd_theta)
      end if
      ! Parent detach last — reverses the enter-data ordering.
   end subroutine barotropic_workstate_exit_data_impl

   pure function barotropic_workstate_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the barotropic fast-loop work state (BTCL_u/v derived-type coeffs excluded) slot (0 when
      !! unallocated).
      class(barotropic_workstate_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = arr_bytes(this%bt_eta) &
               + arr_bytes(this%bt_ubt) &
               + arr_bytes(this%cor_ref_u) &
               + arr_bytes(this%cor_ref_v) &
               + arr_bytes(this%bt_vbt) &
               + arr_bytes(this%bt_H_ref) &
               + arr_bytes(this%ubt_sum) &
               + arr_bytes(this%vbt_sum) &
               + arr_bytes(this%eta_sum) &
               + arr_bytes(this%uhbt_sum) &
               + arr_bytes(this%vhbt_sum) &
               + arr_bytes(this%bt_uhbt) &
               + arr_bytes(this%bt_vhbt) &
               + arr_bytes(this%bt_ubt_end) &
               + arr_bytes(this%bt_vbt_end) &
               + arr_bytes(this%bt_eta_end) &
               + arr_bytes(this%bt_zeta_corner) &
               + arr_bytes(this%bt_ke_centre) &
               + arr_bytes(this%bt_eta_new) &
               + arr_bytes(this%bt_ubt_prev) &
               + arr_bytes(this%bt_vbt_prev) &
               + arr_bytes(this%F_slow_u) &
               + arr_bytes(this%F_slow_v) &
               + arr_bytes(this%F_bt_u) &
               + arr_bytes(this%F_bt_v) &
               + arr_bytes(this%F_bt_u_fast) &
               + arr_bytes(this%F_bt_v_fast) &
               + arr_bytes(this%ubt_at_n) &
               + arr_bytes(this%vbt_at_n) &
               + arr_bytes(this%pbce) &
               + arr_bytes(this%gtot_E) &
               + arr_bytes(this%gtot_W) &
               + arr_bytes(this%gtot_N) &
               + arr_bytes(this%gtot_S) &
               + arr_bytes(this%e_anom) &
               + arr_bytes(this%eta_PF) &
               + arr_bytes(this%visc_rem_u) &
               + arr_bytes(this%visc_rem_v) &
               + arr_bytes(this%bt_rem_u) &
               + arr_bytes(this%bt_rem_v) &
               + arr_bytes(this%lwd_drag_u) &
               + arr_bytes(this%lwd_drag_v) &
               + arr_bytes(this%h_face_up_x) &
               + arr_bytes(this%h_face_up_y) &
               + arr_bytes(this%wd_wet_dyn) &
               + arr_bytes(this%wd_theta) &
               + arr_bytes(this%wd_flux_x) &
               + arr_bytes(this%wd_flux_y) &
               + arr_bytes(this%wd_open_u) &
               + arr_bytes(this%wd_open_v)
   end function barotropic_workstate_bytes

end module rdb_barotropic_workstate
