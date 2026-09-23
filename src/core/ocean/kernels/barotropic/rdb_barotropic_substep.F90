!! Barotropic-substep kernels for the split-explicit ocean driver.
module rdb_barotropic_substep
   !! Two kernels:
   !!
   !!   `barotropic_substep_linear`  — linearized shallow-water with closed
   !!     walls, no Coriolis, no advection.  Used by the
   !!     `test_ocean_barotropic_substep` analytic checks (rest state, constant
   !!     forcing, gravity-wave standing mode) and as the unit-test
   !!     anchor for the dynamics.  Not called by the production
   !!     split driver.
   !!
   !!   `barotropic_substep_nonlinear`   — production substep with
   !!     free-surface continuity (`(H_ref + η)` face thickness),
   !!     Coriolis (`(ζ + f)·v_perp` Sadourny enstrophy-conserving),
   !!     and vector-invariant advection.  Called by the split RK2
   !!     driver inside each outer baroclinic step.
   !!
   !! Both kernels read `force_u`, `force_v` as constant slow forcing
   !! for the inner substeps and write the time-mean (sum / n_steps)
   !! back into `bt_eta`, `bt_work%bt_ubt`, `bt_work%bt_vbt` on exit.
   !! Per-substep running sums land in `bt_work%eta_sum`, `bt_work%ubt_sum`,
   !! `bt_work%vbt_sum`.  The nonlinear path additionally uses
   !! `bt_work%bt_zeta_corner`, `bt_work%bt_ke_centre`, `bt_eta_new` as
   !! per-substep scratch.
   !!
   !! Closed-wall convention on the outer boundary:
   !!   * u_bt at i=1 and i=nx+1 (array edges) stay at zero
   !!   * v_bt at j=1 and j=ny+1 (array edges) stay at zero
   !!   * u_bt at i=nghost+1 and i=nghost+nx_phys+1 (physical walls)
   !!     stay at zero — slow continuity closes here, so the
   !!     barotropic substep must too or `sum_k(h_layer)` drifts from `H+bt_eta`
   !!   * v_bt at j=nghost+1 and j=nghost+ny_phys+1 (physical walls)
   !!     stay at zero — same reason
   !!   * ζ at the outer corner ring is clamped to zero so (ζ+f)·v
   !!     reduces to the f-plane Coriolis term there
   !!
   !! Stability: Forward-Backward Euler ordering — η update first
   !! (consumes u^n, v^n), then u/v updates (consume the just-updated
   !! η).  Stable on the gravity-wave eigenmode for
   !! `dt_inner · √(g·H) / dx ≤ 1`.
   !!
   !! Performance footgun (still relevant): the η-update do-concurrent
   !! uses four distinct `local()` face-thickness variables
   !! (h_face_E/W/N/S) rather than reusing a single `h_face`.  gfortran
   !! 15.1 + `-O3 -funroll-loops` was observed to corrupt the flux
   !! divergence when h_face was reused across `if/else` branches.
   !! See `feedback_gfortran_local_reassign.md` in claude memory.
   use rdb_constants, only: wp, FROUDE_CAP, H_DIV_EPS
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_barotropic_workstate, only: barotropic_workstate_t
   use rdb_bt_cont_type, only: find_uhbt, find_vhbt
   ! rdb_coriolis_adv removed: cor was demoted to a plain f_corner(:,:) dummy
   use rdb_ocean_boundary_types, only: ocean_bc_state_t, OBC_WALL, OBC_OPEN, OBC_CLAMPED, &
                                       OBC_TIDAL, OBC_CHAPMAN, OBC_PERIODIC
   use rdb_ocean_halo, only: ocean_halo_bt_group_2d, ocean_halo_face_x, &
                             ocean_halo_bt_group_2d_wide, &
                             ocean_halo_is_decomposed, &
                             ocean_halo_is_decomposed_x, ocean_halo_is_decomposed_y
   use rdb_ocean_halo_counters, only: oh_count_bt_u_mid, &
                                      oh_count_suppress_on, oh_count_suppress_off
   use rdb_profiler, only: profiler_start, profiler_stop
   implicit none
   private

   public :: barotropic_substep_linear
   public :: barotropic_substep_nonlinear
   public :: barotropic_substep_nonlinear_interior

contains

   pure subroutine barotropic_substep_linear(grid, metrics, bt_work, force_u, force_v, n_steps, dt_inner, &
                                             eta_forcing)
      !! Public only for the unit-test suite (no production module imports it);
      !! ignore when developing production code in other modules.
      !! Linearized barotropic-substep kernel.  Forward-Euler steps
      !! the barotropic state (eta, ubt, vbt) at dt_inner for n_steps
      !! against the constant slow forcing `force_u`, `force_v` (each
      !! at the same C-grid location as ubt, vbt).  Accumulates the
      !! per-step running sums into `ubt_sum`, `vbt_sum`, `eta_sum`;
      !! at the end divides by `n_steps` and stores the time-mean
      !! back into `bt_eta`, `bt_ubt`, `bt_vbt` for the caller.
      !!
      !! Dynamics (closed walls, no Coriolis — linearized gravity-
      !! wave only for unit tests):
      !!
      !!   ∂η/∂t   = -∂(H_ref · u_bt)/∂x - ∂(H_ref · v_bt)/∂y
      !!   ∂u_bt/∂t = -g · ∂η/∂x + force_u
      !!   ∂v_bt/∂t = -g · ∂η/∂y + force_v
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(barotropic_workstate_t), intent(inout) :: bt_work
      ! Explicit-shape (not assumed-shape): a `(:, :)` dummy carries an array
      ! descriptor that must be device-resident inside the `do concurrent`
      ! kernel.  stdpar manages that, but the OpenMP-target variant reads a
      ! non-mapped descriptor → CUDA illegal/misaligned access in the fast
      ! loop.  Explicit shape passes base + dims (no descriptor) and also
      ! avoids the per-launch descriptor-walk memcpys.
      real(wp), intent(in) :: force_u(grid%nx_total + 1, grid%ny_total)
      real(wp), intent(in) :: force_v(grid%nx_total, grid%ny_total + 1)
      integer, intent(in) :: n_steps
      real(wp), intent(in) :: dt_inner
      real(wp), intent(in), optional :: eta_forcing(grid%nx_total, grid%ny_total)
         !! Optional equilibrium-tide elevation (m); when present the PGF
         !! drives grad(eta - eta_forcing).  Absent => bit-identical.

      integer :: i, j, n, nx, ny
      logical :: tide_on
      real(wp) :: inv_n
      real(wp) :: d_eta
      real(wp) :: h_face_E, h_face_W, h_face_N, h_face_S
      real(wp) :: flux_x_R, flux_x_L, flux_y_N, flux_y_S, div_h_u
      real(wp) :: G   !! η-gradient PGF coefficient, sourced from
                      !! `bt_work%g_bt`.  Defaults to 9.81 (full gravity)
                      !! but the gprime path can lower it to g_FS.

      nx = grid%nx_total
      ny = grid%ny_total
      G = bt_work%g_bt
      tide_on = present(eta_forcing)

      ! Reset running sums.
      do concurrent(j=1:ny, i=1:nx)
         bt_work%eta_sum(i, j) = 0.0_wp
      end do
      do concurrent(j=1:ny, i=1:nx + 1)
         bt_work%ubt_sum(i, j) = 0.0_wp
      end do
      do concurrent(j=1:ny + 1, i=1:nx)
         bt_work%vbt_sum(i, j) = 0.0_wp
      end do

      do n = 1, n_steps
         ! ---- Pass 1: eta update at every cell ----
         do concurrent(j=1:ny, i=1:nx) &
            local(h_face_E, h_face_W, h_face_N, h_face_S, &
                  flux_x_R, flux_x_L, flux_y_N, flux_y_S, div_h_u)
            if (i < nx) then
               h_face_E = 0.5_wp*(bt_work%bt_H_ref(i, j) + bt_work%bt_H_ref(i + 1, j))
            else
               h_face_E = bt_work%bt_H_ref(i, j)
            end if
            if (i > 1) then
               h_face_W = 0.5_wp*(bt_work%bt_H_ref(i - 1, j) + bt_work%bt_H_ref(i, j))
            else
               h_face_W = bt_work%bt_H_ref(i, j)
            end if
            if (j < ny) then
               h_face_N = 0.5_wp*(bt_work%bt_H_ref(i, j) + bt_work%bt_H_ref(i, j + 1))
            else
               h_face_N = bt_work%bt_H_ref(i, j)
            end if
            if (j > 1) then
               h_face_S = 0.5_wp*(bt_work%bt_H_ref(i, j - 1) + bt_work%bt_H_ref(i, j))
            else
               h_face_S = bt_work%bt_H_ref(i, j)
            end if
            if (bt_work%use_bt_cont_type) then
               ! MOM6 piecewise-cubic flux closure: transport is capped
               ! by the upstream column's per-layer h sum so the BT mode
               ! can't pump mass through a face where the upstream
               ! column lacks the height to supply it.  BTCL_u/v were
               ! built from the slow ML snapshot in
               ! `set_local_BT_cont_types`.
               flux_x_R = find_uhbt(bt_work%bt_ubt(i + 1, j), bt_work%BTCL_u(i + 1, j))*metrics%dy_cu_bt(i + 1, j)
               flux_x_L = find_uhbt(bt_work%bt_ubt(i, j), bt_work%BTCL_u(i, j))*metrics%dy_cu_bt(i, j)
               flux_y_N = find_vhbt(bt_work%bt_vbt(i, j + 1), bt_work%BTCL_v(i, j + 1))*metrics%dx_cv_bt(i, j + 1)
               flux_y_S = find_vhbt(bt_work%bt_vbt(i, j), bt_work%BTCL_v(i, j))*metrics%dx_cv_bt(i, j)
            else if (bt_work%use_upstream_h_face) then
               ! Upstream-PPM h_face from the slow ML snapshot.  Same
               ! convention slow continuity uses, so the BT mode's
               ! mass flux and the per-layer mass flux carry the
               ! same face thickness — eliminates the centred-vs-
               ! upstream mismatch that leaves phantom bed-layer
               ! velocity at slopes.  Built once per outer step in
               ! `compute_h_face_upstream`.  bt_H_ref unused on
               ! this branch — the upstream column sum carries the
               ! total thickness (incl. η at top of stage).
               flux_x_R = bt_work%h_face_up_x(i + 1, j)*bt_work%bt_ubt(i + 1, j)*metrics%dy_cu_bt(i + 1, j)
               flux_x_L = bt_work%h_face_up_x(i, j)*bt_work%bt_ubt(i, j)*metrics%dy_cu_bt(i, j)
               flux_y_N = bt_work%h_face_up_y(i, j + 1)*bt_work%bt_vbt(i, j + 1)*metrics%dx_cv_bt(i, j + 1)
               flux_y_S = bt_work%h_face_up_y(i, j)*bt_work%bt_vbt(i, j)*metrics%dx_cv_bt(i, j)
            else
               flux_x_R = h_face_E*bt_work%bt_ubt(i + 1, j)*metrics%dy_cu_bt(i + 1, j)
               flux_x_L = h_face_W*bt_work%bt_ubt(i, j)*metrics%dy_cu_bt(i, j)
               flux_y_N = h_face_N*bt_work%bt_vbt(i, j + 1)*metrics%dx_cv_bt(i, j + 1)
               flux_y_S = h_face_S*bt_work%bt_vbt(i, j)*metrics%dx_cv_bt(i, j)
            end if
            ! Conservative transport divergence · iareaT (= inv_dx/inv_dy on uniform).
            div_h_u = ((flux_x_R - flux_x_L) + (flux_y_N - flux_y_S))*metrics%iareaT(i, j)
            bt_work%bt_eta(i, j) = bt_work%bt_eta(i, j) - dt_inner*div_h_u
         end do

         ! ---- Pass 2: ubt update at interior east faces ----
         ! Single-rank reference kernel; physical-edge gating lives in
         ! the nonlinear production variant (barotropic_substep_nonlinear).
         do concurrent(j=1:ny, i=2:nx) local(d_eta)
            d_eta = bt_work%bt_eta(i, j) - bt_work%bt_eta(i - 1, j)
            if (tide_on) d_eta = d_eta - (eta_forcing(i, j) - eta_forcing(i - 1, j))
            bt_work%bt_ubt(i, j) = bt_work%bt_ubt(i, j) + dt_inner*( &
                                   -G*d_eta*metrics%idxCu(i, j) + &
                                   force_u(i, j))
         end do
         do concurrent(j=1:ny)
            bt_work%bt_ubt(1, j) = 0.0_wp
            bt_work%bt_ubt(nx + 1, j) = 0.0_wp
            ! Physical-wall closure — slow continuity blocks flow at
            ! these faces, so the substep loop must too.  See header.
            bt_work%bt_ubt(grid%nghost + 1, j) = 0.0_wp
            bt_work%bt_ubt(grid%nghost + grid%nx_phys + 1, j) = 0.0_wp
         end do

         ! ---- Pass 3: vbt update at interior north faces ----
         do concurrent(j=2:ny, i=1:nx) local(d_eta)
            d_eta = bt_work%bt_eta(i, j) - bt_work%bt_eta(i, j - 1)
            if (tide_on) d_eta = d_eta - (eta_forcing(i, j) - eta_forcing(i, j - 1))
            bt_work%bt_vbt(i, j) = bt_work%bt_vbt(i, j) + dt_inner*( &
                                   -G*d_eta*metrics%idyCv(i, j) + &
                                   force_v(i, j))
         end do
         do concurrent(i=1:nx)
            bt_work%bt_vbt(i, 1) = 0.0_wp
            bt_work%bt_vbt(i, ny + 1) = 0.0_wp
            bt_work%bt_vbt(i, grid%nghost + 1) = 0.0_wp
            bt_work%bt_vbt(i, grid%nghost + grid%ny_phys + 1) = 0.0_wp
         end do

         ! ---- Accumulate running sums ----
         do concurrent(j=1:ny, i=1:nx)
            bt_work%eta_sum(i, j) = bt_work%eta_sum(i, j) + bt_work%bt_eta(i, j)
         end do
         do concurrent(j=1:ny, i=1:nx + 1)
            bt_work%ubt_sum(i, j) = bt_work%ubt_sum(i, j) + bt_work%bt_ubt(i, j)
         end do
         do concurrent(j=1:ny + 1, i=1:nx)
            bt_work%vbt_sum(i, j) = bt_work%vbt_sum(i, j) + bt_work%bt_vbt(i, j)
         end do
      end do

      ! Snapshot end-of-loop η/u/v BEFORE the time-mean overwrite.
      ! `apply_bt_correction` uses all three so the recombined
      ! per-layer momentum and the h_layer rescale both live at
      ! t + dt_outer (Hallberg 2009).  Mixing end-step velocity with
      ! time-mean SSH puts the two legs half a substep out of phase
      ! and corrupts the gravity-wave dispersion.
      do concurrent(j=1:ny, i=1:nx)
         bt_work%bt_eta_end(i, j) = bt_work%bt_eta(i, j)
      end do
      do concurrent(j=1:ny, i=1:nx + 1)
         bt_work%bt_ubt_end(i, j) = bt_work%bt_ubt(i, j)
      end do
      do concurrent(j=1:ny + 1, i=1:nx)
         bt_work%bt_vbt_end(i, j) = bt_work%bt_vbt(i, j)
      end do

      ! Time-mean: divide running sum by n_steps, store back into bt_*.
      inv_n = 1.0_wp/real(n_steps, wp)
      do concurrent(j=1:ny, i=1:nx)
         bt_work%bt_eta(i, j) = bt_work%eta_sum(i, j)*inv_n
      end do
      do concurrent(j=1:ny, i=1:nx + 1)
         bt_work%bt_ubt(i, j) = bt_work%ubt_sum(i, j)*inv_n
      end do
      do concurrent(j=1:ny + 1, i=1:nx)
         bt_work%bt_vbt(i, j) = bt_work%vbt_sum(i, j)*inv_n
      end do
   end subroutine barotropic_substep_linear

   subroutine barotropic_substep_nonlinear(grid, bt_work, force_u, force_v, n_steps, dt_inner, &
                                           bt_eta, bt_H_ref, bt_eta_new, bt_ke_centre, eta_sum, bt_eta_end, &
                                           bt_ubt, bt_ubt_prev, bt_rem_u, ubt_sum, uhbt_sum, bt_uhbt, bt_ubt_end, &
                                           bt_vbt, bt_vbt_prev, bt_rem_v, vbt_sum, vhbt_sum, bt_vhbt, bt_vbt_end, &
                                           bt_zeta_corner, f_corner, &
                                           area_cu, area_cv, dx_cu, dx_cv, dy_cu, dy_cv, &
                                           iarea_bu, iarea_t, idx_cu, idy_cv, &
                                           bc, t, eta_forcing, bt_halo)
      !! Nonlinear barotropic substep.  Same time-mean accumulator
      !! pattern as `barotropic_substep_linear` but the per-substep
      !! dynamics include the three barotropic nonlinearities that
      !! matter for MOM6-grade physics:
      !!
      !!   1. Free-surface continuity — fluxes use `(H_ref + η)`
      !!      averaged to faces, not just `H_ref`.
      !!
      !!   2. Coriolis on the barotropic mode — plain 4-point
      !!      velocity average at the perpendicular face times
      !!      `f_at_face` from `coriolis_adv_t%f_corner`.
      !!
      !!   3. Vector-invariant momentum advection — relative
      !!      vorticity ζ at corners and KE at cell centres
      !!      recomputed each substep.  Sadourny enstrophy-conserving
      !!      `(ζ + f) · v_perp - ∂(g·η + KE)/∂x` on the u-face.
      !!
      !! OBC dispatch: when `bc` is absent (existing callers), every
      !! physical wall face is hard-zeroed (Phase 3 closure).  When
      !! `bc` is present, each edge dispatches on its tag:
      !!   OBC_WALL / default  -> hard zero (closed)
      !!   OBC_OPEN            -> Flather radiation with η_ext = 0
      !!                          `u_face = ±√(g/H) · η_interior`
      !!                          (sign chosen so the OUTWARD-NORMAL
      !!                          component carries energy out).
      !!
      !! The 21 fast-loop 2D array arguments are passed explicitly so a
      !! march-in caller (Phase 3c) can supply wider shadow arrays while
      !! keeping the same arithmetic body unchanged.  At `bt_halo=0` every
      !! caller passes `bt_work%<field>` — same memory, bit-identical.
      type(hgrid_t), intent(in) :: grid
      ! Promoted metric arrays (explicit-shape, device-resident, added to the
      ! present() clauses below) — replaces the per-substep `metrics%<field>`
      ! derived-type-dummy derefs that forced NVHPC into per-launch
      ! present_or_copyin (descriptor copies dominated the barotropic profile).
      ! Shapes mirror ocean_metrics_t: Cu=(nx+1,ny), Cv=(nx,ny+1), Bu=(nx+1,ny+1).
      real(wp), intent(in) :: area_cu(grid%nx_total + 1, grid%ny_total)
      real(wp), intent(in) :: area_cv(grid%nx_total, grid%ny_total + 1)
      real(wp), intent(in) :: dx_cu(grid%nx_total + 1, grid%ny_total)
      ! Porous barriers: the caller passes `metrics%dx_cv_bt` /
      ! `dy_cu_bt` here, NOT the slow-path widths.  They are a byte copy
      ! of `dx_cv` / `dy_cu` unless porous barriers are on, in which case
      ! they carry the column-integrated open fraction so the barotropic
      ! transport sees the barrier too (otherwise the per-layer
      ! renormalisation to `uhbt` would hand the blocked transport back).
      real(wp), intent(in) :: dx_cv(grid%nx_total, grid%ny_total + 1)
      real(wp), intent(in) :: dy_cu(grid%nx_total + 1, grid%ny_total)
      real(wp), intent(in) :: dy_cv(grid%nx_total, grid%ny_total + 1)
      real(wp), intent(in) :: iarea_bu(grid%nx_total + 1, grid%ny_total + 1)
      real(wp), intent(in) :: iarea_t(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: idx_cu(grid%nx_total + 1, grid%ny_total)
      real(wp), intent(in) :: idy_cv(grid%nx_total, grid%ny_total + 1)
      type(barotropic_workstate_t), intent(inout) :: bt_work
         !! Carries scalars (g_bt, bebt, wetdry_enable, use_bt_cont_type,
         !! use_upstream_h_face, bt_substep_drag, wd_*) and the BTCL_u/v,
         !! h_face_up_x/y, wd_* arrays that are NOT promoted.  The
         !! 22 fast-loop 2D arrays below replace the bt_work%<field> and
         !! cor%<field> derefs; bt_work itself is still needed here.
      ! Explicit-shape (not assumed-shape): a `(:, :)` dummy carries an array
      ! descriptor that must be device-resident inside the `do concurrent`
      ! kernel.  stdpar manages that, but the OpenMP-target variant reads a
      ! non-mapped descriptor → CUDA illegal/misaligned access in the fast
      ! loop.  Explicit shape passes base + dims (no descriptor) and also
      ! avoids the per-launch descriptor-walk memcpys.
      real(wp), intent(in) :: force_u(grid%nx_total + 1, grid%ny_total)
      real(wp), intent(in) :: force_v(grid%nx_total, grid%ny_total + 1)
      integer, intent(in) :: n_steps
      real(wp), intent(in) :: dt_inner
      ! ---- Promoted fast-loop 2D arrays (explicit-shape, same pattern as force_u/v) ----
      ! Centre scalars (nx_total, ny_total):
      real(wp), intent(inout) :: bt_eta(grid%nx_total, grid%ny_total)
         !! Barotropic SSH at cell centres; read+written every substep.
      real(wp), intent(in) :: bt_H_ref(grid%nx_total, grid%ny_total)
         !! Reference column thickness (m); read-only inside the fast loop.
      real(wp), intent(inout) :: bt_eta_new(grid%nx_total, grid%ny_total)
         !! Per-substep η^{n+1} Jacobi scratch; written in Pass 1, read in η-swap.
      real(wp), intent(inout) :: bt_ke_centre(grid%nx_total, grid%ny_total)
         !! Barotropic KE at cell centres; written in Pass 1, read in Pass 2b/2c.
      real(wp), intent(inout) :: eta_sum(grid%nx_total, grid%ny_total)
         !! η time-mean accumulator; zeroed at entry, accumulated each substep.
      real(wp), intent(inout) :: bt_eta_end(grid%nx_total, grid%ny_total)
         !! End-of-loop η snapshot (before time-mean overwrite).
      ! East-face u (nx_total+1, ny_total):
      real(wp), intent(inout) :: bt_ubt(grid%nx_total + 1, grid%ny_total)
         !! Depth-mean u at east faces; read+written every substep.
      real(wp), intent(inout) :: bt_ubt_prev(grid%nx_total + 1, grid%ny_total)
         !! u^{n-1} BEBT projection snapshot; read+written when bebt > 0.
      real(wp), intent(in) :: bt_rem_u(grid%nx_total + 1, grid%ny_total)
         !! BT-substep multiplicative drag factor at u-faces; read-only.
      real(wp), intent(inout) :: ubt_sum(grid%nx_total + 1, grid%ny_total)
         !! u time-mean accumulator; zeroed at entry, accumulated each substep.
      real(wp), intent(inout) :: uhbt_sum(grid%nx_total + 1, grid%ny_total)
         !! Depth-integrated u transport accumulator; zeroed and accumulated.
      real(wp), intent(inout) :: bt_uhbt(grid%nx_total + 1, grid%ny_total)
         !! Time-mean east-face transport (m²/s); written at loop end.
      real(wp), intent(inout) :: bt_ubt_end(grid%nx_total + 1, grid%ny_total)
         !! End-of-loop u snapshot (before time-mean overwrite).
      ! North-face v (nx_total, ny_total+1):
      real(wp), intent(inout) :: bt_vbt(grid%nx_total, grid%ny_total + 1)
         !! Depth-mean v at north faces; read+written every substep.
      real(wp), intent(inout) :: bt_vbt_prev(grid%nx_total, grid%ny_total + 1)
         !! v^{n-1} BEBT projection snapshot; read+written when bebt > 0.
      real(wp), intent(in) :: bt_rem_v(grid%nx_total, grid%ny_total + 1)
         !! BT-substep multiplicative drag factor at v-faces; read-only.
      real(wp), intent(inout) :: vbt_sum(grid%nx_total, grid%ny_total + 1)
         !! v time-mean accumulator; zeroed at entry, accumulated each substep.
      real(wp), intent(inout) :: vhbt_sum(grid%nx_total, grid%ny_total + 1)
         !! Depth-integrated v transport accumulator; zeroed and accumulated.
      real(wp), intent(inout) :: bt_vhbt(grid%nx_total, grid%ny_total + 1)
         !! Time-mean north-face transport (m²/s); written at loop end.
      real(wp), intent(inout) :: bt_vbt_end(grid%nx_total, grid%ny_total + 1)
         !! End-of-loop v snapshot (before time-mean overwrite).
      ! Corner ζ (nx_total+1, ny_total+1):
      real(wp), intent(inout) :: bt_zeta_corner(grid%nx_total + 1, grid%ny_total + 1)
         !! Barotropic relative vorticity at corners; written in Pass 1, read in Pass 2b/2c.
      real(wp), intent(in) :: f_corner(grid%nx_total + 1, grid%ny_total + 1)
         !! Coriolis parameter at corners (from coriolis_adv_t%f_corner); read-only.
         !! Passed explicitly so a wide-halo caller (Phase 3c) can supply a
         !! wide f_corner without touching the arithmetic body.
      type(ocean_bc_state_t), intent(inout), optional :: bc
         !! `intent(inout)` so the routine can update the persistent
         !! `eta_old_chapman_*` scalars used by the OBC_CHAPMAN
         !! radiation BC.  Other BC types (WALL/OPEN/CLAMPED/TIDAL)
         !! treat `bc` as read-only — the inout intent is for the
         !! Chapman state lifecycle.
      real(wp), intent(in), optional :: t
         !! Outer-step wall time (s) used to evaluate the per-edge
         !! tidal constituent table for OBC_TIDAL.  Treated as a
         !! constant across the barotropic substeps (tidal periods are
         !! orders of magnitude longer than the inner dt).  Absent
         !! ⇒ t = 0; the OBC_TIDAL formula degenerates to OBC_OPEN.
      real(wp), intent(in), optional :: eta_forcing(grid%nx_total, grid%ny_total)
         !! Optional equilibrium-tide elevation (m); when present the PGF
         !! drives grad(eta - eta_forcing).  Absent => bit-identical.
      integer, intent(in), optional :: bt_halo
         !! Wide-halo march-in width (Phase 3c).  When > 0 the caller has
         !! supplied wide shadow arrays and `grid` has nghost = nghost + bt_halo.
         !! The mid-substep u exchange is absorbed and the per-substep group
         !! exchange is replaced by one wide exchange every `bt_halo/2` substeps.
         !! Absent or 0 => v1 per-substep exchange path (bit-identical).

      integer :: i, j, n, nx, ny, gw
      logical :: tide_on
      real(wp) :: inv_n
      real(wp) :: d_eta
      real(wp) :: h_face_E, h_face_W, h_face_N, h_face_S
      real(wp) :: flux_x_R, flux_x_L, flux_y_N, flux_y_S, div_h_u
      real(wp) :: ubt_R, ubt_L, vbt_N, vbt_S
         !! BEBT-projected face velocities used in the η update flux.
         !! `(1 + bebt)·ubt^n − bebt·ubt^{n-1}`; collapses to `ubt^n`
         !! when `bebt = 0`.
      real(wp) :: zeta_at_u, zeta_at_v, f_at_u, f_at_v
      real(wp) :: v_at_u, u_at_v, ke_grad_x, ke_grad_y
      real(wp) :: G   !! η-gradient PGF coefficient, sourced from
                      !! `bt_work%g_bt`.  See the linear substep above
                      !! for rationale.
      real(wp) :: bebt
         !! BEBT velocity-projection coefficient, sourced from
         !! `bt_work%bebt`.  `bebt = 0` ⇒ pure forward-backward Euler
         !! transport (bit-identical to pre-knob behaviour); MOM6's
         !! `BT_PROJECT_VELOCITY=True` default is 0.2.
      real(wp) :: w_nl
         !! Live-nonlinear-term weight from `bt_work%substep_zeta_ke`:
         !! 1 (default) ⇒ Pass 2b/2c integrate live `ζ_bt`/∇KE
         !! (`1.0*x = x` exactly — bit-identical); 0 ⇒ MOM6-parity
         !! planetary-only substeps (ζ/KE stay frozen in `force_u/v`;
         !! the Cor_ref subtraction is reduced to match in
         !! `subtract_fast_cor_ref`).
      integer :: bc_w, bc_e, bc_s, bc_n
      logical :: has_w, has_e, has_s, has_n
         !! Physical-domain-edge flags cached from bc%has_west/east/south/north.
         !! .true. (default) => this edge is a physical domain boundary and
         !! wall/BC closures apply.  .false. => this edge is an MPI seam;
         !! the halo fill corrects it, so closures are skipped here.
         !! Default .true. preserves single-rank bit-identity.
      logical :: per_x, per_y
         !! Periodic-axis flags cached from bc%periodic_x/y before the
         !! `do concurrent` loops so they are loop-invariant scalars.
      logical :: do_fold
         !! Tripolar north-fold flag cached from bc%north_fold.  Gates the
         !! inline fold DC loops in the fast loop; .false. ⇒ bit-identical.
      integer :: nf_isum, nf_jsum_c, nf_jsum_v, nf_jfold, nf_jlo_c, nf_imid
         !! Cached fold index constants (centre/v-face/u-face maps).
      integer :: i_w_face, i_e_face, j_s_face, j_n_face
      integer :: i_w_int, i_e_int, j_s_int, j_n_int
      integer :: eta_gx_lo, eta_gx_hi
      real(wp) :: clamped_u_w, clamped_u_e, clamped_v_s, clamped_v_n
      real(wp) :: eta_target_w, eta_target_e, eta_target_s, eta_target_n
      real(wp) :: t_now
      integer :: nc
      real(wp) :: eta_int_mean_w, eta_int_mean_e, eta_int_mean_s, eta_int_mean_n
      real(wp) :: chapman_alpha
      ! Full-Flather (§4, v2): pre-cached scalars for the interior DC.
      logical  :: use_ff       !! .true. = full-Flather form, .false. = legacy.
      real(wp) :: ext_u_w, ext_u_e, ext_v_s, ext_v_n  !! Exterior velocities.
      logical  :: use_nodal_bc
         !! OBC tidal nodal/astronomical correction (C3), cached from
         !! bc%tidal_nodal.  .false. ⇒ legacy static-phase OBC tidal sum
         !! (bit-identical).  .true. ⇒ per-constituent f_c amplitude factor +
         !! (V_c+u_c) equilibrium/nodal phase, with tidal_phase read as a LAG.
      ! Per-iteration locals for full-Flather computation (used in local() clauses).
      real(wp) :: D_int_w, Cg_w, cfl_w, u_inlet_w, ssh_in_w
      real(wp) :: D_int_e, Cg_e, cfl_e, u_inlet_e, ssh_in_e
      real(wp) :: D_int_s, Cg_s, cfl_s, v_inlet_s, ssh_in_s
      real(wp) :: D_int_n, Cg_n, cfl_n, v_inlet_n, ssh_in_n
      ! Wet/dry branch (docs/ocean_wetdry_plan.md).  All locals are
      ! single-assignment per iteration (feedback_gfortran_local_reassign:
      ! gfortran 15.1 corrupts do-concurrent locals reassigned across
      ! if/else branches — hence merge() single-write forms below).
      logical :: wd_on
         !! Loop-invariant cache of `bt_work%wetdry_enable`.
      real(wp) :: h_up_w, h_up_e, h_up_f, cmax_f, u_cap_f
         !! Pass 1w-a upwind face thickness + FROUDE_CAP guard locals.
      real(wp) :: out_rate, out_vol, avail_v
         !! Pass 1w-b per-cell outflow accounting (m³/s, m-of-depth, m).
      real(wp) :: th_face
         !! Pass 1w-c per-face limiter factor `min(theta_L, theta_R)`.
      real(wp) :: d_new
         !! Pass 1w-d post-update total depth for the hysteresis mask.
      real(wp) :: wet_l, wet_r, zb_l, zb_r, open_f
         !! Pass 2 bed-blocking gate locals.
      integer :: iw_f, ie_f, js_f, jn_f
         !! Clamped neighbour-cell indices at array-edge faces.
      logical :: marchin
         !! .true. when the wide-halo march-in path is active (bt_halo > 0).
      integer :: num_cycles
         !! Substeps between wide grouped exchanges (bt_halo/2); loop-invariant.
      integer :: ins_w, ins_e, ins_s, ins_n
         !! Wide-halo PHYSICAL-EDGE insets (cells).  bt_halo at a physical
         !! (non-seam, non-periodic, non-fold) edge under march-in, else 0.
         !! Why: v1's array-edge closures (forced-zero outer faces, one-sided
         !! h_face / v_at_u / u_at_v stencils, ζ ring clamp) sit at the NORMAL
         !! array edge, and the wall-adjacent interior weakly consumes the
         !! ghost-band evolution they shape (ζ/KE chain creeps 2 cells/substep
         !! through a closed wall — only u faces are zeroed there).  On the
         !! wide grid those closures would land bt_halo cells further out,
         !! the ghost-band dynamics change, and the difference is dynamically
         !! UNSTABLE (measured: KE reldiff 4e-11 at 1 step -> 8.4 at 100
         !! steps on the MPI dyn gate).  The insets re-anchor every
         !! array-edge closure at the EFFECTIVE (normal-array) edge, making
         !! the band evolution bit-identical to v1 at physical edges while
         !! seam / periodic / fold edges keep the full wide march.
      integer :: ilo_c, ihi_c, jlo_c, jhi_c
         !! Effective array-edge CELL indices: first/last cell of the
         !! v1-equivalent band (`1 + ins_w`, `nx - ins_e`, ...).  At
         !! bt_halo = 0 these are exactly 1 / nx / 1 / ny — every use below
         !! reduces to the v1 constant, bit-identical.
      logical :: in_band
         !! Pass-1 per-iteration local: cell lies inside the effective band.

      nx = grid%nx_total
      ny = grid%ny_total
      G = bt_work%g_bt
      bebt = bt_work%bebt
      w_nl = merge(1.0_wp, 0.0_wp, bt_work%substep_zeta_ke)
      tide_on = present(eta_forcing)
      wd_on = bt_work%wetdry_enable
      ! NOTE: no single-expression `present(x) .and. x > 0` — Fortran does
      ! not short-circuit .and., and nvfortran evaluates the absent optional
      ! (nil deref, segfault); gfortran happened to tolerate it.
      marchin = .false.
      if (present(bt_halo)) marchin = bt_halo > 0
      num_cycles = 0
      if (marchin) num_cycles = bt_halo/2

      ! Seed the BEBT velocity-projection snapshots from the incoming
      ! state so the first substep's projection is a no-op
      ! (`ubt_trans = ubt^0`).  After the first substep we refresh
      ! these from `bt_ubt / bt_vbt` between Passes 1 and 2b/2c so
      ! the next substep sees the correct `ubt^{n-1}`.  Skip when
      ! `bebt = 0` to keep the per-call cost of the no-knob path at zero.
      if (bebt > 0.0_wp) then
         do concurrent(j=1:ny, i=1:nx + 1)
            bt_ubt_prev(i, j) = bt_ubt(i, j)
         end do
         do concurrent(j=1:ny + 1, i=1:nx)
            bt_vbt_prev(i, j) = bt_vbt(i, j)
         end do
      end if

      ! Physical-wall geometry: face indices + the interior column
      ! adjacent to each wall.  Computed up front so the Chapman BC
      ! branch below (which reads `bt_eta(i_w_int, …)` etc.) doesn't
      ! pick up uninitialised values — used to live below the BC
      ! block and segfaulted on multicore stdpar with even one OBC
      ! edge set to OBC_CHAPMAN.
      i_w_face = grid%nghost + 1
      i_e_face = grid%nghost + grid%nx_phys + 1
      j_s_face = grid%nghost + 1
      j_n_face = grid%nghost + grid%ny_phys + 1
      i_w_int = i_w_face            ! cell column just east of west wall
      i_e_int = i_e_face - 1        ! cell column just west of east wall
      j_s_int = j_s_face            ! cell row    just north of south wall
      j_n_int = j_n_face - 1        ! cell row    just south of north wall

      ! Cache per-edge tags as scalars so the do-concurrent wall
      ! closures don't deref the optional struct on every iteration.
      bc_w = OBC_WALL
      bc_e = OBC_WALL
      bc_s = OBC_WALL
      bc_n = OBC_WALL
      ! Physical-edge defaults: .true. => single-rank bit-identical behaviour.
      has_w = .true.
      has_e = .true.
      has_s = .true.
      has_n = .true.
      per_x = .false.
      per_y = .false.
      do_fold = .false.
      ! Fold index constants (storage maps, Appendix A).
      nf_isum = 2*grid%nghost + grid%nx_phys + 1      ! centre/v i-map sum
      nf_jsum_c = 2*grid%nghost + 2*grid%ny_phys + 1   ! T/u j-halo sum
      nf_jsum_v = 2*grid%nghost + 2*grid%ny_phys       ! v j-map sum
      nf_jfold = grid%nghost + grid%ny_phys            ! v self-conjugate row
      nf_jlo_c = grid%nghost + grid%ny_phys + 1        ! first T/u north halo row
      nf_imid = grid%nghost + (grid%nx_phys + 1)/2     ! v on-row west-half end
      clamped_u_w = 0.0_wp
      clamped_u_e = 0.0_wp
      clamped_v_s = 0.0_wp
      clamped_v_n = 0.0_wp
      eta_target_w = 0.0_wp
      eta_target_e = 0.0_wp
      eta_target_s = 0.0_wp
      eta_target_n = 0.0_wp
      ! Full-Flather: default to legacy (bit-identical when bc absent).
      use_ff = .false.
      ext_u_w = 0.0_wp
      ext_u_e = 0.0_wp
      ext_v_s = 0.0_wp
      ext_v_n = 0.0_wp
      use_nodal_bc = .false.
      t_now = 0.0_wp
      if (present(t)) t_now = t
      if (present(bc)) then
         bc_w = bc%west%bc_type
         bc_e = bc%east%bc_type
         bc_s = bc%south%bc_type
         bc_n = bc%north%bc_type
         per_x = bc%periodic_x .and. .not. ocean_halo_is_decomposed_x()
         per_y = bc%periodic_y .and. .not. ocean_halo_is_decomposed_y()
         do_fold = bc%north_fold
         ! Physical-edge flags: .false. at an MPI seam (halo fills it);
         ! .true. at a physical domain edge (wall/BC closure applies).
         has_w = bc%has_west
         has_e = bc%has_east
         has_s = bc%has_south
         has_n = bc%has_north
         clamped_u_w = bc%west%clamped_u
         clamped_u_e = bc%east%clamped_u
         clamped_v_s = bc%south%clamped_v
         clamped_v_n = bc%north%clamped_v
         ! Full-Flather exterior velocities (§4, v2).
         use_ff = bc%use_full_flather
         ext_u_w = bc%ext_u_west
         ext_u_e = bc%ext_u_east
         ext_v_s = bc%ext_v_south
         ext_v_n = bc%ext_v_north
         ! OBC_TIDAL eta target — edge-uniform sum of constituents.
         ! Frozen during the barotropic substeps (tidal periods are
         ! O(hr), substep-block duration is O(s), error is negligible).
         !
         ! Nodal/astronomical correction (C3): when use_nodal_bc, apply the
         ! per-constituent 18.6-yr amplitude factor tidal_fnodal (f_c) and the
         ! equilibrium+nodal phase tidal_arg (V_c+u_c), with tidal_phase read as
         ! a Greenwich phase LAG (subtracted, MOM6 OBC_TIDE convention).  When
         ! off, tidal_fnodal≡1 / tidal_arg≡0 and the legacy additive-phase sum
         ! is reproduced exactly (bit-identical).  The branch is loop-invariant
         ! and uniform across the edge (~free on GPU); kept inside the loop.
         use_nodal_bc = bc%tidal_nodal
         do nc = 1, bc%west%n_tidal_constituents
            if (use_nodal_bc) then
               eta_target_w = eta_target_w + bc%west%tidal_fnodal(nc)*bc%west%tidal_amp(nc)* &
                              cos(bc%west%tidal_omega(nc)*t_now + bc%west%tidal_arg(nc) &
                                  - bc%west%tidal_phase(nc))
            else
               eta_target_w = eta_target_w + bc%west%tidal_amp(nc)* &
                              cos(bc%west%tidal_omega(nc)*t_now + bc%west%tidal_phase(nc))
            end if
         end do
         do nc = 1, bc%east%n_tidal_constituents
            if (use_nodal_bc) then
               eta_target_e = eta_target_e + bc%east%tidal_fnodal(nc)*bc%east%tidal_amp(nc)* &
                              cos(bc%east%tidal_omega(nc)*t_now + bc%east%tidal_arg(nc) &
                                  - bc%east%tidal_phase(nc))
            else
               eta_target_e = eta_target_e + bc%east%tidal_amp(nc)* &
                              cos(bc%east%tidal_omega(nc)*t_now + bc%east%tidal_phase(nc))
            end if
         end do
         do nc = 1, bc%south%n_tidal_constituents
            if (use_nodal_bc) then
               eta_target_s = eta_target_s + bc%south%tidal_fnodal(nc)*bc%south%tidal_amp(nc)* &
                              cos(bc%south%tidal_omega(nc)*t_now + bc%south%tidal_arg(nc) &
                                  - bc%south%tidal_phase(nc))
            else
               eta_target_s = eta_target_s + bc%south%tidal_amp(nc)* &
                              cos(bc%south%tidal_omega(nc)*t_now + bc%south%tidal_phase(nc))
            end if
         end do
         do nc = 1, bc%north%n_tidal_constituents
            if (use_nodal_bc) then
               eta_target_n = eta_target_n + bc%north%tidal_fnodal(nc)*bc%north%tidal_amp(nc)* &
                              cos(bc%north%tidal_omega(nc)*t_now + bc%north%tidal_arg(nc) &
                                  - bc%north%tidal_phase(nc))
            else
               eta_target_n = eta_target_n + bc%north%tidal_amp(nc)* &
                              cos(bc%north%tidal_omega(nc)*t_now + bc%north%tidal_phase(nc))
            end if
         end do

         ! Chapman target — Sommerfeld-style radiation update from
         ! the persistent edge-mean η to the current interior η, at
         ! the gravity-wave phase speed.  Edge-uniform scalar form
         ! for v1 (per-cell adaptive Orlanski is a follow-on).
         !
         !     eta_target = eta_old + α · (eta_int_mean − eta_old)
         !     α = c · n_steps · dt_inner / dx  where c = √(g·H_ref)
         !
         ! α is clamped to [0, 1]; α = 1 ⇒ wall fully tracks interior
         ! (perfect outflow), α = 0 ⇒ wall holds eta_old (closed-ish).
         if (bc_w == OBC_CHAPMAN .or. bc_e == OBC_CHAPMAN .or. &
             bc_s == OBC_CHAPMAN .or. bc_n == OBC_CHAPMAN) then
            chapman_alpha = min(1.0_wp, &
                                sqrt(G*max(bt_H_ref(grid%nghost + 1, grid%nghost + 1), 1.0e-6_wp))* &
                                real(n_steps, wp)*dt_inner*idx_cu(grid%nghost + 1, grid%nghost + 1))
            if (bc_w == OBC_CHAPMAN) then
               eta_int_mean_w = sum(bt_eta(i_w_int, &
                                           grid%nghost + 1:grid%nghost + grid%ny_phys))/ &
                                real(grid%ny_phys, wp)
               eta_target_w = bc%eta_old_chapman_w + &
                              chapman_alpha*(eta_int_mean_w - bc%eta_old_chapman_w)
            end if
            if (bc_e == OBC_CHAPMAN) then
               eta_int_mean_e = sum(bt_eta(i_e_int, &
                                           grid%nghost + 1:grid%nghost + grid%ny_phys))/ &
                                real(grid%ny_phys, wp)
               eta_target_e = bc%eta_old_chapman_e + &
                              chapman_alpha*(eta_int_mean_e - bc%eta_old_chapman_e)
            end if
            if (bc_s == OBC_CHAPMAN) then
               eta_int_mean_s = sum(bt_eta(grid%nghost + 1:grid%nghost + grid%nx_phys, &
                                           j_s_int))/ &
                                real(grid%nx_phys, wp)
               eta_target_s = bc%eta_old_chapman_s + &
                              chapman_alpha*(eta_int_mean_s - bc%eta_old_chapman_s)
            end if
            if (bc_n == OBC_CHAPMAN) then
               eta_int_mean_n = sum(bt_eta(grid%nghost + 1:grid%nghost + grid%nx_phys, &
                                           j_n_int))/ &
                                real(grid%nx_phys, wp)
               eta_target_n = bc%eta_old_chapman_n + &
                              chapman_alpha*(eta_int_mean_n - bc%eta_old_chapman_n)
            end if
         end if
      end if

      ! Physical-edge insets (see the ins_* declaration comment).  Needs the
      ! cached has_* / per_* / do_fold flags, so computed after the bc block.
      ! Seams (has_* = .false.), periodic edges (local wrap owns the band),
      ! and the tripolar fold keep inset 0 — full wide march.
      ins_w = 0
      ins_e = 0
      ins_s = 0
      ins_n = 0
      if (marchin) then
         if (has_w .and. .not. per_x) ins_w = bt_halo
         if (has_e .and. .not. per_x) ins_e = bt_halo
         if (has_s .and. .not. per_y) ins_s = bt_halo
         if (has_n .and. .not. (per_y .or. do_fold)) ins_n = bt_halo
      end if
      ilo_c = 1 + ins_w
      ihi_c = nx - ins_e
      jlo_c = 1 + ins_s
      jhi_c = ny - ins_n

      ! One structured data region enclosing the setup kernels, the whole
      ! `do n = 1, n_steps` substep loop, and the end-of-loop time-mean pass so
      ! the ~24 promoted explicit-shape dummy arrays resolve PRESENT once per
      ! substep-call instead of triggering a present_or_copyin descriptor upload
      ! at every per-substep `!$acc kernels` launch (the dominant np2 scaling
      ! brake: two async regions x ~65 substeps x 2 stages x 576 steps, one
      ! upload per promoted array dummy each).  Only the always-mapped promoted
      ! arrays are listed; the lazily-mapped wetdry / BTCL / upstream-h
      ! `bt_work` components and the optional `eta_forcing` keep their implicit
      ! present_or_copyin handling — they are genuinely unmapped when their
      ! feature is off, so a blanket `default(present)` would fault the default
      ! (wetdry-off, tides-off) configuration at region entry.
      !$acc data present(force_u, force_v, bt_eta, bt_H_ref, bt_eta_new, &
      !$acc              bt_ke_centre, eta_sum, bt_eta_end, bt_ubt, bt_ubt_prev, &
      !$acc              bt_rem_u, ubt_sum, uhbt_sum, bt_uhbt, bt_ubt_end, &
      !$acc              bt_vbt, bt_vbt_prev, bt_rem_v, vbt_sum, vhbt_sum, &
      !$acc              bt_vhbt, bt_vbt_end, bt_zeta_corner, f_corner, &
      !$acc              area_cu, area_cv, dx_cu, dx_cv, dy_cu, dy_cv, &
      !$acc              iarea_bu, iarea_t, idx_cu, idy_cv)
      !$acc kernels async(1)
      do concurrent(j=1:ny, i=1:nx)
         eta_sum(i, j) = 0.0_wp
      end do
      do concurrent(j=1:ny, i=1:nx + 1)
         ubt_sum(i, j) = 0.0_wp
         uhbt_sum(i, j) = 0.0_wp
      end do
      do concurrent(j=1:ny + 1, i=1:nx)
         vbt_sum(i, j) = 0.0_wp
         vhbt_sum(i, j) = 0.0_wp
      end do
      !$acc end kernels

      do n = 1, n_steps
         ! Drained by the !$acc wait(1) after the n_steps loop.
         ! Explicit present() on the always-mapped promoted arrays drops the
         ! per-launch copy-fallback pre-check (the lazily-mapped bt_work
         ! components / optional eta_forcing keep implicit present_or_copyin).
         !$acc kernels async(1) &
         !$acc   present(force_u, force_v, bt_eta, bt_H_ref, bt_eta_new, &
         !$acc           bt_ke_centre, eta_sum, bt_eta_end, bt_ubt, bt_ubt_prev, &
         !$acc           bt_rem_u, ubt_sum, uhbt_sum, bt_uhbt, bt_ubt_end, &
         !$acc           bt_vbt, bt_vbt_prev, bt_rem_v, vbt_sum, vhbt_sum, &
         !$acc           bt_vhbt, bt_vbt_end, bt_zeta_corner, f_corner, &
         !$acc           area_cu, area_cv, dx_cu, dx_cv, dy_cu, dy_cv, &
         !$acc           iarea_bu, iarea_t, idx_cu, idy_cv)
         if (wd_on) then
            ! ---- Pass 1w (wet/dry): positive-definite upwind flux form ----
            ! Replaces the centred Pass 1 below when `&ocean_wetdry_nml
            ! enable` (docs/ocean_wetdry_plan.md §3).  Four sweeps on the
            ! same async queue:
            !   1w-a  provisional face fluxes — UPWIND total depth (mass is
            !         drawn from the cell that has it) + the FROUDE_CAP
            !         thin-face runaway-velocity guard,
            !   1w-b  per-cell outflow limiter theta = min(1, avail/outflow),
            !   1w-c  scale each face by min(theta_L, theta_R) — guarantees
            !         D >= 0 every substep (one sweep suffices: a cell's
            !         realised outflow <= theta_c · outflow_c <= avail_c;
            !         inflow only adds) — and accumulate the LIMITED
            !         transport into uhbt_sum/vhbt_sum so the downstream BT
            !         correction distributes what actually flowed,
            !   1w-d  divergence -> eta_new (+ round-off armour at -H_ref),
            !         hysteresis wet-mask update, and the same KE / zeta
            !         fills the centred path computes (Pass 2 reads both).
            ! Pure-inflow (rewetting) cells keep theta = 1 by construction —
            ! the limiter can never block a dry cell from refilling.
            do concurrent(j=1:ny, i=1:nx + 1) &
               local(ubt_R, iw_f, ie_f, h_up_w, h_up_e, h_up_f, cmax_f, u_cap_f)
               ubt_R = (1.0_wp + bebt)*bt_ubt(i, j) - bebt*bt_ubt_prev(i, j)
               iw_f = max(i - 1, 1)
               ie_f = min(i, nx)
               h_up_w = max(bt_H_ref(iw_f, j) + bt_eta(iw_f, j), 0.0_wp)
               h_up_e = max(bt_H_ref(ie_f, j) + bt_eta(ie_f, j), 0.0_wp)
               h_up_f = merge(h_up_w, h_up_e, ubt_R >= 0.0_wp)
               ! FROUDE_CAP is a runaway guard, NOT a physics limiter — the
               ! Fr=1 variant measurably destroyed supercritical runup
               ! (prototype: u error 102%); the sqrt floor only keeps the
               ! argument positive at zero-depth faces.
               cmax_f = FROUDE_CAP*sqrt(G*max(h_up_f, 0.5_wp*bt_work%wd_dry_depth))
               u_cap_f = sign(min(abs(ubt_R), cmax_f), ubt_R)
               bt_work%wd_flux_x(i, j) = h_up_f*u_cap_f*dy_cu(i, j)
            end do
            do concurrent(j=1:ny + 1, i=1:nx) &
               local(vbt_N, js_f, jn_f, h_up_w, h_up_e, h_up_f, cmax_f, u_cap_f)
               vbt_N = (1.0_wp + bebt)*bt_vbt(i, j) - bebt*bt_vbt_prev(i, j)
               js_f = max(j - 1, 1)
               jn_f = min(j, ny)
               h_up_w = max(bt_H_ref(i, js_f) + bt_eta(i, js_f), 0.0_wp)
               h_up_e = max(bt_H_ref(i, jn_f) + bt_eta(i, jn_f), 0.0_wp)
               h_up_f = merge(h_up_w, h_up_e, vbt_N >= 0.0_wp)
               cmax_f = FROUDE_CAP*sqrt(G*max(h_up_f, 0.5_wp*bt_work%wd_dry_depth))
               u_cap_f = sign(min(abs(vbt_N), cmax_f), vbt_N)
               bt_work%wd_flux_y(i, j) = h_up_f*u_cap_f*dx_cv(i, j)
            end do
            do concurrent(j=1:ny, i=1:nx) &
               local(out_rate, out_vol, avail_v)
               out_rate = max(bt_work%wd_flux_x(i + 1, j), 0.0_wp) + &
                          max(-bt_work%wd_flux_x(i, j), 0.0_wp) + &
                          max(bt_work%wd_flux_y(i, j + 1), 0.0_wp) + &
                          max(-bt_work%wd_flux_y(i, j), 0.0_wp)
               out_vol = dt_inner*out_rate*iarea_t(i, j)
               avail_v = max(bt_H_ref(i, j) + bt_eta(i, j), 0.0_wp)
               ! theta = 1 unless this substep would drain more depth than
               ! the cell holds; H_DIV_EPS is pure 1/0 armour (out_vol >
               ! avail_v >= 0 implies out_vol > 0 whenever the ratio is used).
               bt_work%wd_theta(i, j) = merge(avail_v/max(out_vol, H_DIV_EPS), 1.0_wp, &
                                              out_vol > avail_v)
            end do
            do concurrent(j=1:ny, i=1:nx + 1) local(th_face)
               th_face = min(bt_work%wd_theta(max(i - 1, 1), j), &
                             bt_work%wd_theta(min(i, nx), j))
               bt_work%wd_flux_x(i, j) = bt_work%wd_flux_x(i, j)*th_face
               uhbt_sum(i, j) = uhbt_sum(i, j) + bt_work%wd_flux_x(i, j)
            end do
            do concurrent(j=1:ny + 1, i=1:nx) local(th_face)
               th_face = min(bt_work%wd_theta(i, max(j - 1, 1)), &
                             bt_work%wd_theta(i, min(j, ny)))
               bt_work%wd_flux_y(i, j) = bt_work%wd_flux_y(i, j)*th_face
               vhbt_sum(i, j) = vhbt_sum(i, j) + bt_work%wd_flux_y(i, j)
            end do
            do concurrent(j=1:ny, i=1:nx) local(div_h_u, d_new)
               div_h_u = ((bt_work%wd_flux_x(i + 1, j) - bt_work%wd_flux_x(i, j)) + &
                          (bt_work%wd_flux_y(i, j + 1) - bt_work%wd_flux_y(i, j)))* &
                         iarea_t(i, j)
               ! Round-off armour only: the limiter keeps D >= 0 analytically;
               ! the max() clamps float round-off at exactly-drained cells.
               ! It must NEVER inject finite water (the thin-film floor broke
               ! conservation at 2e-2 in the prototype; no-floor conserves to
               ! round-off).
               bt_eta_new(i, j) = max(bt_eta(i, j) - dt_inner*div_h_u, &
                                      -bt_H_ref(i, j))
               d_new = bt_H_ref(i, j) + bt_eta_new(i, j)
               ! Hysteresis wet mask: wet above rewet_depth, dry below
               ! dry_depth, HOLD previous state in the band (kills wet/dry
               ! front chatter).  Own-cell read+write only — race-free.
               if (d_new > bt_work%wd_rewet_depth) then
                  bt_work%wd_wet_dyn(i, j) = 1.0_wp
               else if (d_new < bt_work%wd_dry_depth) then
                  bt_work%wd_wet_dyn(i, j) = 0.0_wp
               end if
               ! KE + interior zeta: identical to the centred path (Pass 2b/2c
               ! consume both regardless of the Pass-1 branch).
               bt_ke_centre(i, j) = 0.25_wp*iarea_t(i, j)*( &
                                    area_cu(i, j)*bt_ubt(i, j)**2 + &
                                    area_cu(i + 1, j)*bt_ubt(i + 1, j)**2 + &
                                    area_cv(i, j)*bt_vbt(i, j)**2 + &
                                    area_cv(i, j + 1)*bt_vbt(i, j + 1)**2)
               if (i >= 2 .and. j >= 2) then
                  bt_zeta_corner(i, j) = &
                     ((bt_vbt(i, j)*dy_cv(i, j) - bt_vbt(i - 1, j)*dy_cv(i - 1, j)) - &
                      (bt_ubt(i, j)*dx_cu(i, j) - bt_ubt(i, j - 1)*dx_cu(i, j - 1)))* &
                     iarea_bu(i, j)
               end if
            end do
         else
            ! ---- Pass 1: η update with (H_ref + η) face thickness ----
            ! Each cell owns its east (i+1) and north (j+1) faces for
            ! the transport accumulator, so the per-face sums are
            ! race-free.  West/south walls (face indices 1) stay at
            ! zero from init — the closed-wall BC makes them.
            !
            ! Face velocity uses the MOM6 BEBT projection:
            !   `ubt_trans = (1 + bebt)·ubt^n − bebt·ubt^{n-1}`
            ! — a forward-time extrapolation that lets the η evolution
            ! anticipate the velocity update later in the substep.  At
            ! `bebt = 0` (pure FB; the default is MOM6's 0.1) the formula collapses to
            ! `ubt^n` and Pass 1 is bit-identical to the pre-knob FBE
            ! scheme.  See `bt_work%bebt` doc for rationale.
            do concurrent(j=1:ny, i=1:nx) &
               local(h_face_E, h_face_W, h_face_N, h_face_S, &
                     ubt_R, ubt_L, vbt_N, vbt_S, &
                     flux_x_R, flux_x_L, flux_y_N, flux_y_S, div_h_u, in_band)
               ! March-in physical-edge emulation: cells beyond the effective
               ! band do not exist in v1 — skip them entirely (their η/KE/ζ
               ! are never consumed by band cells once the closures below are
               ! anchored at ilo_c/ihi_c/jlo_c/jhi_c, and skipping keeps their
               ! scattered uhbt_sum/vhbt_sum writes from contaminating the
               ! effective-edge face sums).  At bt_halo = 0 the band is the
               ! whole array — the branch is uniformly true, bit-identical.
               in_band = i >= ilo_c .and. i <= ihi_c .and. &
                         j >= jlo_c .and. j <= jhi_c
               if (in_band) then
               if (i < ihi_c) then
                  h_face_E = 0.5_wp*((bt_H_ref(i, j) + bt_eta(i, j)) + &
                                     (bt_H_ref(i + 1, j) + bt_eta(i + 1, j)))
               else
                  h_face_E = bt_H_ref(i, j) + bt_eta(i, j)
               end if
               if (i > ilo_c) then
                  h_face_W = 0.5_wp*((bt_H_ref(i - 1, j) + bt_eta(i - 1, j)) + &
                                     (bt_H_ref(i, j) + bt_eta(i, j)))
               else
                  h_face_W = bt_H_ref(i, j) + bt_eta(i, j)
               end if
               if (j < jhi_c) then
                  h_face_N = 0.5_wp*((bt_H_ref(i, j) + bt_eta(i, j)) + &
                                     (bt_H_ref(i, j + 1) + bt_eta(i, j + 1)))
               else
                  h_face_N = bt_H_ref(i, j) + bt_eta(i, j)
               end if
               if (j > jlo_c) then
                  h_face_S = 0.5_wp*((bt_H_ref(i, j - 1) + bt_eta(i, j - 1)) + &
                                     (bt_H_ref(i, j) + bt_eta(i, j)))
               else
                  h_face_S = bt_H_ref(i, j) + bt_eta(i, j)
               end if
               ubt_R = (1.0_wp + bebt)*bt_ubt(i + 1, j) - bebt*bt_ubt_prev(i + 1, j)
               ubt_L = (1.0_wp + bebt)*bt_ubt(i, j) - bebt*bt_ubt_prev(i, j)
               vbt_N = (1.0_wp + bebt)*bt_vbt(i, j + 1) - bebt*bt_vbt_prev(i, j + 1)
               vbt_S = (1.0_wp + bebt)*bt_vbt(i, j) - bebt*bt_vbt_prev(i, j)
               if (bt_work%use_bt_cont_type) then
                  ! Flux-bounded BT continuity — see the linear-substep
                  ! comment for rationale.  Fed the BEBT-projected face
                  ! velocities so BT_cont and BT_PROJECT_VELOCITY compose
                  ! (bebt=0 ⇒ projection is a no-op, bit-identical).
                  flux_x_R = find_uhbt(ubt_R, bt_work%BTCL_u(i + 1, j))*dy_cu(i + 1, j)
                  flux_x_L = find_uhbt(ubt_L, bt_work%BTCL_u(i, j))*dy_cu(i, j)
                  flux_y_N = find_vhbt(vbt_N, bt_work%BTCL_v(i, j + 1))*dx_cv(i, j + 1)
                  flux_y_S = find_vhbt(vbt_S, bt_work%BTCL_v(i, j))*dx_cv(i, j)
               else if (bt_work%use_upstream_h_face) then
                  ! Upstream-PPM h_face from the slow ML snapshot.  Held
                  ! constant across the substep — `h_face_up_x` was filled
                  ! from `Σ_k h_layer(upstream)` at the top of the stage,
                  ! which equals `(H_ref + η)_upstream_at_top`.  We deliberately
                  ! drop the centred-η anomaly term the plan's pseudocode
                  ! shows because adding `0.5·(η_W + η_E)` here would
                  ! double-count `η_at_top` that's already inside
                  ! `h_face_up_x`.  The η evolution during the substep
                  ! shifts the face thickness by O(δη) ~ sub-cm vs H ~
                  ! 100s of metres — well below the upstream-vs-centred
                  ! correction we're after.  Uses the BEBT-
                  ! projected velocities for consistency with the other
                  ! branches (bebt=0 ⇒ bit-identical).
                  flux_x_R = bt_work%h_face_up_x(i + 1, j)*ubt_R*dy_cu(i + 1, j)
                  flux_x_L = bt_work%h_face_up_x(i, j)*ubt_L*dy_cu(i, j)
                  flux_y_N = bt_work%h_face_up_y(i, j + 1)*vbt_N*dx_cv(i, j + 1)
                  flux_y_S = bt_work%h_face_up_y(i, j)*vbt_S*dx_cv(i, j)
               else
                  flux_x_R = h_face_E*ubt_R*dy_cu(i + 1, j)
                  flux_x_L = h_face_W*ubt_L*dy_cu(i, j)
                  flux_y_N = h_face_N*vbt_N*dx_cv(i, j + 1)
                  flux_y_S = h_face_S*vbt_S*dx_cv(i, j)
               end if
               ! Conservative transport divergence · iareaT (= inv_dx/inv_dy on uniform).
               div_h_u = ((flux_x_R - flux_x_L) + (flux_y_N - flux_y_S))*iarea_t(i, j)
               bt_eta_new(i, j) = bt_eta(i, j) - dt_inner*div_h_u
               ! MOM6-style transport accumulator.  Cell (i,j) owns its
               ! east face (i+1, j) and north face (i, j+1).  Uses the
               ! BEBT-projected face velocity for consistency with the
               ! η flux divergence above.
               uhbt_sum(i + 1, j) = uhbt_sum(i + 1, j) + flux_x_R
               vhbt_sum(i, j + 1) = vhbt_sum(i, j + 1) + flux_y_N
               ! KE (centre) and interior ζ (NE corner) read only the
               ! pre-update ubt/vbt and write disjoint arrays, so they ride
               ! along Pass 1's sweep — bit-identical, one fewer launch each.
               ! Wall-ζ closure stays below (it overwrites the corners).
               bt_ke_centre(i, j) = 0.25_wp*iarea_t(i, j)*( &
                                    area_cu(i, j)*bt_ubt(i, j)**2 + &
                                    area_cu(i + 1, j)*bt_ubt(i + 1, j)**2 + &
                                    area_cv(i, j)*bt_vbt(i, j)**2 + &
                                    area_cv(i, j + 1)*bt_vbt(i, j + 1)**2)
               if (i >= ilo_c + 1 .and. j >= jlo_c + 1) then
                  bt_zeta_corner(i, j) = &
                     ((bt_vbt(i, j)*dy_cv(i, j) - bt_vbt(i - 1, j)*dy_cv(i - 1, j)) - &
                      (bt_ubt(i, j)*dx_cu(i, j) - bt_ubt(i, j - 1)*dx_cu(i, j - 1)))* &
                     iarea_bu(i, j)
               end if
               end if   ! in_band
            end do
         end if
         ! η-swap + ζ wall closure — one barrier sweep over the corner
         ! grid. Both must finish before Pass 2b reads η/ζ neighbours, and
         ! they write disjoint arrays (η centres vs ζ corner ring), so the
         ! merge is bit-identical. Guards reproduce the old per-loop index
         ! sets; the interior ζ written in Pass 1 is untouched here.
         do concurrent(j=1:ny + 1, i=1:nx + 1)
            if (i <= nx .and. j <= ny) bt_eta(i, j) = bt_eta_new(i, j)
            ! ζ outer ring → 0 (free-slip); physical walls → 0 so (ζ+f)·v
            ! reduces to f·v there.
            ! §3 (v2): zero ζ at boundary corner line for ALL non-PERIODIC tags
            ! (OPEN, WALL, TIDAL, CLAMPED, CHAPMAN, …).  At OPEN/TIDAL/CHAPMAN
            ! edges the ghost velocities evolve under unphysical array-edge
            ! dynamics and inject noise into the first-interior (ζ+f)·v.
            ! Zeroing matches the "zero relative vorticity at the open boundary"
            ! closure — the Coriolis term reduces to f·v exactly as at WALL.
            ! PERIODIC retains the computed value so Coriolis advects across
            ! the seam with full vorticity.
            ! KE at the first ghost cell feeds ke_grad at non-WALL faces, but
            ! for every non-WALL tag the wall-face velocity is overwritten by
            ! the Flather dispatch after Pass 2b/2c, so the Pass-2b KE
            ! contribution is discarded — KE needs no change here (design §3).
            ! Array outer ring: outside the ghost band, always land,
            ! decomposition-invariant.  Under march-in the EFFECTIVE ring
            ! (normal-array edge) is clamped at physical edges — ilo_c/ihi_c
            ! reduce to 1/nx at bt_halo = 0 (bit-identical).
            if (i == ilo_c .or. i == ihi_c + 1) bt_zeta_corner(i, j) = 0.0_wp
            if (j == jlo_c .or. j == jhi_c + 1) bt_zeta_corner(i, j) = 0.0_wp
            ! Physical-wall corner-zeta closure: gate on has_* so an MPI
            ! seam face is left for the halo to correct (bit-identical at
            ! has_*=.true. / single-rank; no-op at seam with has_*=.false.).
            if (bc_w /= OBC_PERIODIC .and. has_w .and. i == grid%nghost + 1) bt_zeta_corner(i, j) = 0.0_wp
            if (bc_e /= OBC_PERIODIC .and. has_e .and. i == grid%nghost + grid%nx_phys + 1) bt_zeta_corner(i, j) = 0.0_wp
            if (bc_s /= OBC_PERIODIC .and. has_s .and. j == grid%nghost + 1) bt_zeta_corner(i, j) = 0.0_wp
            if (bc_n /= OBC_PERIODIC .and. has_n .and. j == grid%nghost + grid%ny_phys + 1) bt_zeta_corner(i, j) = 0.0_wp
         end do
         ! Periodic η ghost-wrap (design §1.5 step 3).  Required before
         ! Pass 2b because the u-update at the west wall face reads
         ! bt_eta(nghost, j) — a ghost whose Pass-1 update used the
         ! zeroed array-edge fluxes and thus holds garbage without this wrap.
         ! Inline DC (not a subroutine call) to stay inside async(1).
         ! Two sequential passes (x first, then y): a fused loop races at
         ! the doubly-periodic corner ghosts — the y-branch reads x-ghost
         ! columns that other iterations' x-branch is writing.  The y-pass
         ! reading the already-x-wrapped columns is also what makes the
         ! corner ghosts correct (same convention as rdb_ocean_periodic).
         if (per_x) then
            do concurrent(j=1:ny, i=1:nx)
               if (i <= grid%nghost) then
                  bt_eta(i, j) = bt_eta(i + grid%nx_phys, j)
               end if
               if (i > grid%nx_phys + grid%nghost) then
                  bt_eta(i, j) = bt_eta(i - grid%nx_phys, j)
               end if
            end do
         end if
         if (per_y) then
            do concurrent(j=1:ny, i=1:nx)
               if (j <= grid%nghost) then
                  bt_eta(i, j) = bt_eta(i, j + grid%ny_phys)
               end if
               if (j > grid%ny_phys + grid%nghost) then
                  bt_eta(i, j) = bt_eta(i, j - grid%ny_phys)
               end if
            end do
         end if
         ! Tripolar north-fold of η (centre, copy): runs AFTER the periodic-x
         ! η wrap (Appendix A: fold reads cyclically-wrapped corner columns).
         ! Fills the north halo rows from the reflected interior; the on-line
         ! T row is strictly below the seam (pure halo image), so no on-row op.
         if (do_fold) then
            do concurrent(j=nf_jlo_c:ny, i=1:nx)
               bt_eta(i, j) = bt_eta(nf_isum - i, nf_jsum_c - j)
            end do
         end if

         ! Open-edge η ghost fill — zero-gradient extrapolation from the first
         ! interior cell.  The missing analogue of the periodic wrap above: at
         ! an OPEN/TIDAL/CHAPMAN/CLAMPED edge the η ghost is otherwise left to
         ! evolve under Pass-1's zeroed array-edge fluxes (see the "garbage"
         ! note before the periodic wrap).  Left unfilled, the η CORNER ghost
         ! feeds the downstream barotropic correction (target_h = H_ref +
         ! bt_eta_end) at the corner column, which blows up h_layer there.
         ! x-pass runs over the FULL j extent (covers corner rows); y-pass over
         ! the cross-extent CLIPPED to the physical span unless the adjacent
         ! x-edge is also open — so an open×open corner is filled by the y-pass
         ! reading the x-filled column, while a periodic/wall x-edge ghost is
         ! never clobbered (keeps periodic+open mixed configs, e.g. the Eady
         ! channel, bit-identical).  Must run before Pass 2b reads bt_eta(nghost).
         ! WALL/PERIODIC tags are not open-ish ⇒ no DC launched ⇒ no-op.
         ! x-pass: west / east ghost columns over the full j extent.
         if (bc_w == OBC_OPEN .or. bc_w == OBC_TIDAL .or. &
             bc_w == OBC_CHAPMAN .or. bc_w == OBC_CLAMPED) then
            do concurrent(j=1:ny, i=1:grid%nghost)
               bt_eta(i, j) = bt_eta(i_w_int, j)
            end do
         end if
         if (bc_e == OBC_OPEN .or. bc_e == OBC_TIDAL .or. &
             bc_e == OBC_CHAPMAN .or. bc_e == OBC_CLAMPED) then
            do concurrent(j=1:ny, i=i_e_int + 1:nx)
               bt_eta(i, j) = bt_eta(i_e_int, j)
            end do
         end if
         ! y-pass: south / north ghost rows over the corner-safe i extent
         ! Extend into any NON-PERIODIC x-edge ghost (wall or open) so an open
         ! south/north edge fills its η ghost CORNER even against a wall side;
         ! only PERIODIC is excluded (owns its ghost via the wrap).
         eta_gx_lo = i_w_int
         eta_gx_hi = i_e_int
         if (bc_w /= OBC_PERIODIC) eta_gx_lo = 1
         if (bc_e /= OBC_PERIODIC) eta_gx_hi = nx
         if (bc_s == OBC_OPEN .or. bc_s == OBC_TIDAL .or. &
             bc_s == OBC_CHAPMAN .or. bc_s == OBC_CLAMPED) then
            do concurrent(j=1:grid%nghost, i=eta_gx_lo:eta_gx_hi)
               bt_eta(i, j) = bt_eta(i, j_s_int)
            end do
         end if
         if (bc_n == OBC_OPEN .or. bc_n == OBC_TIDAL .or. &
             bc_n == OBC_CHAPMAN .or. bc_n == OBC_CLAMPED) then
            do concurrent(j=j_n_int + 1:ny, i=eta_gx_lo:eta_gx_hi)
               bt_eta(i, j) = bt_eta(i, j_n_int)
            end do
         end if

         ! BEBT projection prep: stash the current (`ubt^n`) before
         ! Pass 2b/2c overwrites it with `ubt^{n+1}`.  Next substep's
         ! Pass 1 will then read `ubt_prev = ubt^n` for the
         ! `(1+bebt)·ubt^{n+1} − bebt·ubt^n` extrapolation.  Skip
         ! when `bebt = 0` to keep the no-knob path at zero extra
         ! per-step memory bandwidth.
         if (bebt > 0.0_wp) then
            do concurrent(j=1:ny, i=1:nx + 1)
               bt_ubt_prev(i, j) = bt_ubt(i, j)
            end do
            do concurrent(j=1:ny + 1, i=1:nx)
               bt_vbt_prev(i, j) = bt_vbt(i, j)
            end do
         end if

         ! ---- Pass 2b: u_bt update at interior east faces ----
         do concurrent(j=1:ny, i=2:nx) &
            local(zeta_at_u, f_at_u, v_at_u, ke_grad_x, d_eta, &
                  wet_l, wet_r, zb_l, zb_r, open_f)
            zeta_at_u = w_nl*0.5_wp*(bt_zeta_corner(i, j) + bt_zeta_corner(i, j + 1))
            f_at_u = 0.5_wp*(f_corner(i, j) + f_corner(i, j + 1))
            ! One-sided rows anchored at the EFFECTIVE array edge (jlo_c/jhi_c
            ! reduce to 1/ny at bt_halo = 0 — bit-identical).
            if (j > jlo_c .and. j < jhi_c) then
               v_at_u = 0.25_wp*(bt_vbt(i - 1, j) + bt_vbt(i - 1, j + 1) + &
                                 bt_vbt(i, j) + bt_vbt(i, j + 1))
            else if (j == jlo_c) then
               v_at_u = 0.5_wp*(bt_vbt(i - 1, j + 1) + bt_vbt(i, j + 1))
            else
               v_at_u = 0.5_wp*(bt_vbt(i - 1, j) + bt_vbt(i, j))
            end if
            ke_grad_x = w_nl*(bt_ke_centre(i, j) - bt_ke_centre(i - 1, j))*idx_cu(i, j)
            d_eta = bt_eta(i, j) - bt_eta(i - 1, j)
            if (tide_on) d_eta = d_eta - (eta_forcing(i, j) - eta_forcing(i - 1, j))
            ! Multiplicative drag damping (MOM6 bt_rem_u):
            ! ubt_new = bt_rem · (ubt_old + dt·forces).  When the
            ! `bt_substep_drag` knob is off, `bt_rem_u` stays at 1.0
            ! (init value) and this collapses to the standard FBE update.
            bt_ubt(i, j) = bt_rem_u(i, j)*( &
                           bt_ubt(i, j) + dt_inner*( &
                           (zeta_at_u + f_at_u)*v_at_u &
                           - G*d_eta*idx_cu(i, j) &
                           - ke_grad_x &
                           + force_u(i, j)))
            if (wd_on) then
               ! Bed-blocking gate (wet/dry, plan §3.4): a face into a dry
               ! cell is a WALL unless the wet side's surface stands above
               ! the dry side's bed elevation (+ dry_depth headroom) —
               ! without this the drying-bank face feels a spurious PGF from
               ! the ghost surface (eta_dry ≈ bed) and the front lags the
               ! analytic shoreline.  Uniform branch inside the one DC loop
               ! (splitting per-case doubles GPU launches).  open_f is a
               ! single-assignment merge — no local reassignment.
               wet_l = bt_work%wd_wet_dyn(i - 1, j)
               wet_r = bt_work%wd_wet_dyn(i, j)
               zb_l = -bt_H_ref(i - 1, j)
               zb_r = -bt_H_ref(i, j)
               open_f = merge(1.0_wp, 0.0_wp, &
                              (wet_l > 0.5_wp .and. wet_r > 0.5_wp) .or. &
                              (wet_l > 0.5_wp .and. &
                               bt_eta(i - 1, j) > zb_r + bt_work%wd_dry_depth) .or. &
                              (wet_r > 0.5_wp .and. &
                               bt_eta(i, j) > zb_l + bt_work%wd_dry_depth))
               bt_work%wd_open_u(i, j) = open_f
               bt_ubt(i, j) = bt_ubt(i, j)*open_f
            end if
         end do
         ! BC dispatch for west and east u-faces.  has_w/has_e wrap the
         ! ENTIRE select-case: at an MPI seam (has_*=.false.) no physical-BC
         ! treatment of any kind fires — not even OPEN/TIDAL/CLAMPED/CHAPMAN
         ! branches, which carry GLOBAL-domain tags and would impose physical
         ! values on an interior seam face.  The halo exchange corrects the
         ! face after the substep.  At a physical edge (has_*=.true., the
         ! default) behaviour is unchanged (bit-identical).  The guards remain
         ! INSIDE the do concurrent so the loop stays on the async(1) ACC
         ! kernels region; has_* are loop-invariant scalars cached above.
         ! Array outer-ring zeros (i=1, i=nx+1) are decomposition-invariant
         ! and must stay UNGATED — they are written unconditionally above.
         do concurrent(j=1:ny) &
            local(D_int_w, Cg_w, cfl_w, u_inlet_w, ssh_in_w, &
                  D_int_e, Cg_e, cfl_e, u_inlet_e, ssh_in_e)
            ! Effective array-edge face zeros (ilo_c = 1, ihi_c + 1 = nx + 1
            ! at bt_halo = 0 — bit-identical).  Under march-in these land on
            ! the NORMAL array-edge faces at physical edges, exactly where v1
            ! zeroes them.
            bt_ubt(ilo_c, j) = 0.0_wp
            bt_ubt(ihi_c + 1, j) = 0.0_wp
            ! Physical-wall closure — dispatch on bc tag.  Default
            ! (OBC_WALL) keeps the Phase 3 hard zero; OBC_OPEN applies
            ! Flather radiation with η_ext = 0.  Outward-normal sign
            ! convention: west wall outward = -x (so u_wall < 0 for
            ! a west-going wave with η > 0); east wall outward = +x.
            ! OBC_PERIODIC: leave the computed face values from Pass 2b
            ! interior update and wrap ghost faces inline below.
            ! Full-Flather form (§4, v2; Flather 1976, half-characteristic):
            !   cfl     = dt_inner · Cg / dx,   Cg = sqrt(G · D_int)
            !   u_inlet = cfl·ubt(I-1,j) + (1-cfl)·ubt(I,j)
            !   ssh_in  = η(i_int,j) + (0.5-cfl)·(η(i_int,j) - η(i_int-1,j))
            !   u_b = 0.5·[(u_inlet + u_ext) + (Cg/D_int)·(ssh_in - η_target)]
            !
            ! Sign convention (outward-normal):
            !   West: outward = -x.  Legacy form: u = -sqrt(G/D)·(η-η_tgt).
            !         Full form: for u_ext=0, cfl→0, ssh_in→η_int:
            !           u_b = 0.5·[u(I,j) + (-Cg/D)·(η_int - η_tgt)]
            !         Factor-½ vs legacy is expected (the half-characteristic form;
            !         the full form is the physically correct one, legacy is an
            !         approximation that doubles the response — this is why the
            !         knob defaults to legacy: changing it alters existing results).
            !   East: outward = +x.  Legacy: u = +sqrt(G/D)·(η-η_tgt).
            !         Full form gives u_b = 0.5·[u(I,j) + (+Cg/D)·(η_int - η_tgt)].
            !   South/North: analogous with v and ±y.
            if (has_w) then
               select case (bc_w)
               case (OBC_OPEN)
                  if (use_ff) then
                     D_int_w = max(bt_H_ref(i_w_int, j) + bt_eta(i_w_int, j), 1.0e-6_wp)
                     Cg_w = sqrt(G*D_int_w)
                     cfl_w = min(dt_inner*Cg_w*idx_cu(i_w_face, j), 1.0_wp)
                     u_inlet_w = cfl_w*bt_ubt(i_w_int + 1, j) + &
                                 (1.0_wp - cfl_w)*bt_ubt(i_w_face, j)
                     ssh_in_w = bt_eta(i_w_int, j) + &
                                (0.5_wp - cfl_w)*(bt_eta(i_w_int, j) - &
                                                  bt_eta(i_w_int + 1, j))
                     bt_ubt(i_w_face, j) = 0.5_wp*( &
                                           (u_inlet_w + ext_u_w) - (Cg_w/D_int_w)*ssh_in_w)
                  else
                     bt_ubt(i_w_face, j) = &
                        -sqrt(G/max(bt_H_ref(i_w_int, j) + bt_eta(i_w_int, j), 1.0e-6_wp))* &
                        bt_eta(i_w_int, j)
                  end if
               case (OBC_TIDAL, OBC_CHAPMAN)
                  if (use_ff) then
                     D_int_w = max(bt_H_ref(i_w_int, j) + bt_eta(i_w_int, j), 1.0e-6_wp)
                     Cg_w = sqrt(G*D_int_w)
                     cfl_w = min(dt_inner*Cg_w*idx_cu(i_w_face, j), 1.0_wp)
                     u_inlet_w = cfl_w*bt_ubt(i_w_int + 1, j) + &
                                 (1.0_wp - cfl_w)*bt_ubt(i_w_face, j)
                     ssh_in_w = bt_eta(i_w_int, j) + &
                                (0.5_wp - cfl_w)*(bt_eta(i_w_int, j) - &
                                                  bt_eta(i_w_int + 1, j))
                     bt_ubt(i_w_face, j) = 0.5_wp*( &
                                           (u_inlet_w + ext_u_w) - (Cg_w/D_int_w)*(ssh_in_w - eta_target_w))
                  else
                     bt_ubt(i_w_face, j) = &
                        -sqrt(G/max(bt_H_ref(i_w_int, j) + bt_eta(i_w_int, j), 1.0e-6_wp))* &
                        (bt_eta(i_w_int, j) - eta_target_w)
                  end if
               case (OBC_CLAMPED)
                  bt_ubt(i_w_face, j) = clamped_u_w
               case (OBC_PERIODIC)
                  ! Leave the computed value — it will be overwritten by the
                  ! ghost wrap below; kept here to not write zero over it.
                  continue
               case default
                  bt_ubt(i_w_face, j) = 0.0_wp
               end select
            end if
            if (has_e) then
               select case (bc_e)
               case (OBC_OPEN)
                  if (use_ff) then
                     D_int_e = max(bt_H_ref(i_e_int, j) + bt_eta(i_e_int, j), 1.0e-6_wp)
                     Cg_e = sqrt(G*D_int_e)
                     cfl_e = min(dt_inner*Cg_e*idx_cu(i_e_face, j), 1.0_wp)
                     u_inlet_e = cfl_e*bt_ubt(i_e_int, j) + &
                                 (1.0_wp - cfl_e)*bt_ubt(i_e_face, j)
                     ssh_in_e = bt_eta(i_e_int, j) + &
                                (0.5_wp - cfl_e)*(bt_eta(i_e_int, j) - &
                                                  bt_eta(i_e_int - 1, j))
                     bt_ubt(i_e_face, j) = 0.5_wp*( &
                                           (u_inlet_e + ext_u_e) + (Cg_e/D_int_e)*ssh_in_e)
                  else
                     bt_ubt(i_e_face, j) = &
                        +sqrt(G/max(bt_H_ref(i_e_int, j) + bt_eta(i_e_int, j), 1.0e-6_wp))* &
                        bt_eta(i_e_int, j)
                  end if
               case (OBC_TIDAL, OBC_CHAPMAN)
                  if (use_ff) then
                     D_int_e = max(bt_H_ref(i_e_int, j) + bt_eta(i_e_int, j), 1.0e-6_wp)
                     Cg_e = sqrt(G*D_int_e)
                     cfl_e = min(dt_inner*Cg_e*idx_cu(i_e_face, j), 1.0_wp)
                     u_inlet_e = cfl_e*bt_ubt(i_e_int, j) + &
                                 (1.0_wp - cfl_e)*bt_ubt(i_e_face, j)
                     ssh_in_e = bt_eta(i_e_int, j) + &
                                (0.5_wp - cfl_e)*(bt_eta(i_e_int, j) - &
                                                  bt_eta(i_e_int - 1, j))
                     bt_ubt(i_e_face, j) = 0.5_wp*( &
                                           (u_inlet_e + ext_u_e) + (Cg_e/D_int_e)*(ssh_in_e - eta_target_e))
                  else
                     bt_ubt(i_e_face, j) = &
                        +sqrt(G/max(bt_H_ref(i_e_int, j) + bt_eta(i_e_int, j), 1.0e-6_wp))* &
                        (bt_eta(i_e_int, j) - eta_target_e)
                  end if
               case (OBC_CLAMPED)
                  bt_ubt(i_e_face, j) = clamped_u_e
               case (OBC_PERIODIC)
                  ! Leave the computed value.
                  continue
               case default
                  bt_ubt(i_e_face, j) = 0.0_wp
               end select
            end if
         end do
         ! Periodic u ghost-wrap + belt-and-braces seam copy (design §1.5 step 4).
         ! Inline DC inside async(1) — wraps ghost faces and ensures
         ! u(i_w_face) == u(i_e_face) bit-for-bit.  Generic in nghost: the
         ! g-loop covers all nghost ghost faces on each side (nghost is a
         ! loop-invariant scalar, legal inside the async region).
         if (per_x) then
            do concurrent(j=1:ny, gw=1:grid%nghost)
               ! West ghost faces i=1..nghost ← interior face i+nx_phys
               bt_ubt(gw, j) = bt_ubt(gw + grid%nx_phys, j)
               ! East ghost faces i=nx_phys+nghost+2..nx+1 ← interior face i-nx_phys.
               ! i_e_face = nghost+nx_phys+1, so the ghost faces beyond it are
               ! i_e_face+1 .. i_e_face+nghost == nx+1.
               bt_ubt(i_e_face + gw, j) = &
                  bt_ubt(i_e_face + gw - grid%nx_phys, j)
            end do
            ! Belt-and-braces: enforce i_w == i_e seam identity.
            do concurrent(j=1:ny)
               bt_ubt(i_e_face, j) = bt_ubt(i_w_face, j)
            end do
         end if
         ! Periodic-y wrap of the u GHOST ROWS.  ubt is centre-type in y
         ! (rows 1..nghost and ny_phys+nghost+1..ny are ghosts).  Without
         ! this, next substep's Pass-1 ζ at the y-seam corners reads ubt
         ! ghost rows that evolved under the zeroed array-edge ζ ring —
         ! the corruption creeps inward one row per substep and reaches
         ! the seam vorticity from substep nghost+2 onward.  Runs after
         ! the x-wrap so doubly-periodic corner faces read the already
         ! x-wrapped columns.
         if (per_y) then
            do concurrent(gw=1:grid%nghost, i=1:nx + 1)
               bt_ubt(i, gw) = bt_ubt(i, gw + grid%ny_phys)
               bt_ubt(i, grid%ny_phys + grid%nghost + gw) = &
                  bt_ubt(i, grid%nghost + gw)
            end do
         end if
         ! Tripolar north-fold of ubt (u-face Cu, NEGATE — true vector).
         ! u i-map is the symmetric f' = ni+2-f; reads the already periodic-x
         ! wrapped faces.  Halo rows only (u points lie strictly below the seam).
         if (do_fold) then
            do concurrent(j=nf_jlo_c:ny, i=1:nx + 1)
               bt_ubt(i, j) = -bt_ubt(nf_isum + 1 - i, nf_jsum_c - j)
            end do
         end if

         ! Mid-substep u seam exchange (D4 mirror of the periodic u ghost-wrap
         ! above, design 1.5 step 4).  The substep is ALTERNATING forward-
         ! backward: Pass 2c's Coriolis u_at_v consumes the u that Pass 2b
         ! JUST WROTE.  Under decomposition the seam-adjacent u (owned seam
         ! face + ghost faces/rows) must carry the OWNER's post-2b values
         ! before 2c reads them -- the end-of-substep grouped exchange is too
         ! late.  Skipping this exchange is an O4-measured seam instability:
         ! uniform-T seamount at rest (f-plane, auto n_inner ~65, px=2) grows
         ! spurious KE at ~e^1.5/step from round-off and NaNs by step ~24;
         ! with the exchange it holds machine-zero KE.  Flat-bottom, f=0, or
         ! small-n_inner configs are blind to it (why the O2/O3 gates missed
         ! it).  Gated on a REAL decomposition: single-rank (incl. periodic,
         ! whose local wrap ran above) keeps the async(1) queue un-drained --
         ! bit-identical, no extra sync.
         !$acc end kernels
         ! Mid-substep u seam exchange (D4 mirror of the periodic u ghost-wrap
         ! above).  Absorbed at bt_halo>0: the 2-cell/substep stencil budget
         ! already covers Pass 2c's u_at_v consumption inside the wide band.
         if (.not. marchin .and. ocean_halo_is_decomposed()) then
            call profiler_start("ocean_comms_bt")
            !$acc wait(1)
            call oh_count_bt_u_mid()
            call oh_count_suppress_on()
            call ocean_halo_face_x(bt_ubt)
            call oh_count_suppress_off()
            call profiler_stop("ocean_comms_bt")
         end if
         !$acc kernels async(1) &
         !$acc   present(force_u, force_v, bt_eta, bt_H_ref, bt_eta_new, &
         !$acc           bt_ke_centre, eta_sum, bt_eta_end, bt_ubt, bt_ubt_prev, &
         !$acc           bt_rem_u, ubt_sum, uhbt_sum, bt_uhbt, bt_ubt_end, &
         !$acc           bt_vbt, bt_vbt_prev, bt_rem_v, vbt_sum, vhbt_sum, &
         !$acc           bt_vhbt, bt_vbt_end, bt_zeta_corner, f_corner, &
         !$acc           area_cu, area_cv, dx_cu, dx_cv, dy_cu, dy_cv, &
         !$acc           iarea_bu, iarea_t, idx_cu, idy_cv)
         ! ---- Pass 2c: v_bt update at interior north faces ----
         do concurrent(j=2:ny, i=1:nx) &
            local(zeta_at_v, f_at_v, u_at_v, ke_grad_y, d_eta, &
                  wet_l, wet_r, zb_l, zb_r, open_f)
            zeta_at_v = w_nl*0.5_wp*(bt_zeta_corner(i, j) + bt_zeta_corner(i + 1, j))
            f_at_v = 0.5_wp*(f_corner(i, j) + f_corner(i + 1, j))
            ! One-sided columns anchored at the EFFECTIVE array edge
            ! (ilo_c/ihi_c reduce to 1/nx at bt_halo = 0 — bit-identical).
            if (i > ilo_c .and. i < ihi_c) then
               u_at_v = 0.25_wp*(bt_ubt(i, j - 1) + bt_ubt(i + 1, j - 1) + &
                                 bt_ubt(i, j) + bt_ubt(i + 1, j))
            else if (i == ilo_c) then
               u_at_v = 0.5_wp*(bt_ubt(i + 1, j - 1) + bt_ubt(i + 1, j))
            else
               u_at_v = 0.5_wp*(bt_ubt(i, j - 1) + bt_ubt(i, j))
            end if
            ke_grad_y = w_nl*(bt_ke_centre(i, j) - bt_ke_centre(i, j - 1))*idy_cv(i, j)
            d_eta = bt_eta(i, j) - bt_eta(i, j - 1)
            if (tide_on) d_eta = d_eta - (eta_forcing(i, j) - eta_forcing(i, j - 1))
            ! Multiplicative drag damping (MOM6 bt_rem_v).  See u-side
            ! comment above; same no-op-when-knob-off semantics.
            bt_vbt(i, j) = bt_rem_v(i, j)*( &
                           bt_vbt(i, j) + dt_inner*( &
                           -(zeta_at_v + f_at_v)*u_at_v &
                           - G*d_eta*idy_cv(i, j) &
                           - ke_grad_y &
                           + force_v(i, j)))
            if (wd_on) then
               ! Bed-blocking gate — v-face mirror of the Pass-2b comment.
               wet_l = bt_work%wd_wet_dyn(i, j - 1)
               wet_r = bt_work%wd_wet_dyn(i, j)
               zb_l = -bt_H_ref(i, j - 1)
               zb_r = -bt_H_ref(i, j)
               open_f = merge(1.0_wp, 0.0_wp, &
                              (wet_l > 0.5_wp .and. wet_r > 0.5_wp) .or. &
                              (wet_l > 0.5_wp .and. &
                               bt_eta(i, j - 1) > zb_r + bt_work%wd_dry_depth) .or. &
                              (wet_r > 0.5_wp .and. &
                               bt_eta(i, j) > zb_l + bt_work%wd_dry_depth))
               bt_work%wd_open_v(i, j) = open_f
               bt_vbt(i, j) = bt_vbt(i, j)*open_f
            end if
         end do
         ! BC dispatch for south and north v-faces.  has_s/has_n wrap the
         ! ENTIRE select-case (same rationale as the west/east comment
         ! above): a seam face receives NO physical-BC treatment of any kind;
         ! the halo corrects it after the substep.  Guards remain INSIDE the
         ! do concurrent to stay on the async(1) ACC kernels region; has_*
         ! are loop-invariant scalars cached above.  Array outer-ring zeros
         ! (j=1, j=ny+1) are decomposition-invariant and stay UNGATED.
         do concurrent(i=1:nx) &
            local(D_int_s, Cg_s, cfl_s, v_inlet_s, ssh_in_s, &
                  D_int_n, Cg_n, cfl_n, v_inlet_n, ssh_in_n)
            ! Effective array-edge face zeros (jlo_c = 1, jhi_c + 1 = ny + 1
            ! at bt_halo = 0 — bit-identical; see the u-side comment).
            bt_vbt(i, jlo_c) = 0.0_wp
            bt_vbt(i, jhi_c + 1) = 0.0_wp
            if (has_s) then
               select case (bc_s)
               case (OBC_OPEN)
                  if (use_ff) then
                     D_int_s = max(bt_H_ref(i, j_s_int) + bt_eta(i, j_s_int), 1.0e-6_wp)
                     Cg_s = sqrt(G*D_int_s)
                     cfl_s = min(dt_inner*Cg_s*idy_cv(i, j_s_face), 1.0_wp)
                     v_inlet_s = cfl_s*bt_vbt(i, j_s_int + 1) + &
                                 (1.0_wp - cfl_s)*bt_vbt(i, j_s_face)
                     ssh_in_s = bt_eta(i, j_s_int) + &
                                (0.5_wp - cfl_s)*(bt_eta(i, j_s_int) - &
                                                  bt_eta(i, j_s_int + 1))
                     bt_vbt(i, j_s_face) = 0.5_wp*( &
                                           (v_inlet_s + ext_v_s) - (Cg_s/D_int_s)*ssh_in_s)
                  else
                     bt_vbt(i, j_s_face) = &
                        -sqrt(G/max(bt_H_ref(i, j_s_int) + bt_eta(i, j_s_int), 1.0e-6_wp))* &
                        bt_eta(i, j_s_int)
                  end if
               case (OBC_TIDAL, OBC_CHAPMAN)
                  if (use_ff) then
                     D_int_s = max(bt_H_ref(i, j_s_int) + bt_eta(i, j_s_int), 1.0e-6_wp)
                     Cg_s = sqrt(G*D_int_s)
                     cfl_s = min(dt_inner*Cg_s*idy_cv(i, j_s_face), 1.0_wp)
                     v_inlet_s = cfl_s*bt_vbt(i, j_s_int + 1) + &
                                 (1.0_wp - cfl_s)*bt_vbt(i, j_s_face)
                     ssh_in_s = bt_eta(i, j_s_int) + &
                                (0.5_wp - cfl_s)*(bt_eta(i, j_s_int) - &
                                                  bt_eta(i, j_s_int + 1))
                     bt_vbt(i, j_s_face) = 0.5_wp*( &
                                           (v_inlet_s + ext_v_s) - (Cg_s/D_int_s)*(ssh_in_s - eta_target_s))
                  else
                     bt_vbt(i, j_s_face) = &
                        -sqrt(G/max(bt_H_ref(i, j_s_int) + bt_eta(i, j_s_int), 1.0e-6_wp))* &
                        (bt_eta(i, j_s_int) - eta_target_s)
                  end if
               case (OBC_CLAMPED)
                  bt_vbt(i, j_s_face) = clamped_v_s
               case (OBC_PERIODIC)
                  continue
               case default
                  bt_vbt(i, j_s_face) = 0.0_wp
               end select
            end if
            if (has_n) then
               select case (bc_n)
               case (OBC_OPEN)
                  if (use_ff) then
                     D_int_n = max(bt_H_ref(i, j_n_int) + bt_eta(i, j_n_int), 1.0e-6_wp)
                     Cg_n = sqrt(G*D_int_n)
                     cfl_n = min(dt_inner*Cg_n*idy_cv(i, j_n_face), 1.0_wp)
                     v_inlet_n = cfl_n*bt_vbt(i, j_n_int) + &
                                 (1.0_wp - cfl_n)*bt_vbt(i, j_n_face)
                     ssh_in_n = bt_eta(i, j_n_int) + &
                                (0.5_wp - cfl_n)*(bt_eta(i, j_n_int) - &
                                                  bt_eta(i, j_n_int - 1))
                     bt_vbt(i, j_n_face) = 0.5_wp*( &
                                           (v_inlet_n + ext_v_n) + (Cg_n/D_int_n)*ssh_in_n)
                  else
                     bt_vbt(i, j_n_face) = &
                        +sqrt(G/max(bt_H_ref(i, j_n_int) + bt_eta(i, j_n_int), 1.0e-6_wp))* &
                        bt_eta(i, j_n_int)
                  end if
               case (OBC_TIDAL, OBC_CHAPMAN)
                  if (use_ff) then
                     D_int_n = max(bt_H_ref(i, j_n_int) + bt_eta(i, j_n_int), 1.0e-6_wp)
                     Cg_n = sqrt(G*D_int_n)
                     cfl_n = min(dt_inner*Cg_n*idy_cv(i, j_n_face), 1.0_wp)
                     v_inlet_n = cfl_n*bt_vbt(i, j_n_int) + &
                                 (1.0_wp - cfl_n)*bt_vbt(i, j_n_face)
                     ssh_in_n = bt_eta(i, j_n_int) + &
                                (0.5_wp - cfl_n)*(bt_eta(i, j_n_int) - &
                                                  bt_eta(i, j_n_int - 1))
                     bt_vbt(i, j_n_face) = 0.5_wp*( &
                                           (v_inlet_n + ext_v_n) + (Cg_n/D_int_n)*(ssh_in_n - eta_target_n))
                  else
                     bt_vbt(i, j_n_face) = &
                        +sqrt(G/max(bt_H_ref(i, j_n_int) + bt_eta(i, j_n_int), 1.0e-6_wp))* &
                        (bt_eta(i, j_n_int) - eta_target_n)
                  end if
               case (OBC_CLAMPED)
                  bt_vbt(i, j_n_face) = clamped_v_n
               case (OBC_PERIODIC)
                  continue
               case default
                  bt_vbt(i, j_n_face) = 0.0_wp
               end select
            end if
         end do
         ! Periodic-x wrap of the v GHOST COLUMNS.  vbt is centre-type in
         ! x (columns 1..nghost and nx_phys+nghost+1..nx are ghosts).
         ! Mirror of the ubt ghost-row wrap above: without it, next
         ! substep's ζ at the x-seam corners reads vbt ghost columns that
         ! diverged from their interior partners (zeroed array-edge ζ
         ! ring), and the corruption creeps into the seam momentum from
         ! substep nghost+2 onward — confirmed by the shifted-domain
         ! bit-identity test at N_INNER > nghost+1.  Runs before the
         ! y-wrap so doubly-periodic corner faces wrap correctly.
         if (per_x) then
            do concurrent(j=1:ny + 1, gw=1:grid%nghost)
               bt_vbt(gw, j) = bt_vbt(gw + grid%nx_phys, j)
               bt_vbt(grid%nx_phys + grid%nghost + gw, j) = &
                  bt_vbt(grid%nghost + gw, j)
            end do
         end if
         ! Periodic v ghost-wrap + belt-and-braces seam copy (design §1.5 step 5).
         ! Generic in nghost (g-loop covers all ghost faces on each side).
         if (per_y) then
            do concurrent(i=1:nx, gw=1:grid%nghost)
               ! South ghost faces j=1..nghost ← interior face j+ny_phys
               bt_vbt(i, gw) = bt_vbt(i, gw + grid%ny_phys)
               ! North ghost faces j=ny_phys+nghost+2..ny+1 ← interior face j-ny_phys
               bt_vbt(i, j_n_face + gw) = &
                  bt_vbt(i, j_n_face + gw - grid%ny_phys)
            end do
            ! Belt-and-braces: enforce j_s == j_n seam identity.
            do concurrent(i=1:nx)
               bt_vbt(i, j_n_face) = bt_vbt(i, j_s_face)
            end do
         end if
         ! Tripolar north-fold of vbt (v-face Cv, NEGATE + on-row projection).
         ! Two ops (Appendix A): (1) halo rows above the fold line filled from
         ! the reflected+negated interior; (2) the self-conjugate fold row
         ! j=nf_jfold antisymmetrised — west half overwritten from the negated
         ! east mirror, self-fixed column (odd ni) → 0.  Runs after periodic-x.
         if (do_fold) then
            do concurrent(j=nf_jfold + 1:ny + 1, i=1:nx)
               bt_vbt(i, j) = -bt_vbt(nf_isum - i, nf_jsum_v - j)
            end do
            do concurrent(i=grid%nghost + 1:nf_imid)
               if (nf_isum - i == i) then
                  bt_vbt(i, nf_jfold) = 0.0_wp
               else
                  bt_vbt(i, nf_jfold) = -bt_vbt(nf_isum - i, nf_jfold)
               end if
            end do
         end if

         ! Time-mean accumulators — three disjoint, independent writes
         ! fused over the staggered union (η centres, u/v faces) with
         ! per-array index guards; one launch instead of three.
         do concurrent(j=1:ny + 1, i=1:nx + 1)
            if (i <= nx .and. j <= ny) then
               eta_sum(i, j) = eta_sum(i, j) + bt_eta(i, j)
            end if
            if (j <= ny) then
               ubt_sum(i, j) = ubt_sum(i, j) + bt_ubt(i, j)
            end if
            if (i <= nx) then
               vbt_sum(i, j) = vbt_sum(i, j) + bt_vbt(i, j)
            end if
         end do
         !$acc end kernels
         ! End-of-substep ghost refresh.
         ! v1 path (marchin=.false.): drain async(1) then one normal grouped
         ! exchange per substep (bit-identical; D0).
         ! March-in path (marchin=.true.): fire a WIDE grouped exchange every
         ! num_cycles substeps (bt_halo/2); skip the last substep (n==n_steps)
         ! since copy_out immediately follows.  The wide band absorbs the
         ! mid-substep u exchange.
         if (marchin) then
            if (mod(n, num_cycles) == 0 .and. n < n_steps) then
               call profiler_start("ocean_comms_bt")
               !$acc wait(1)
               call ocean_halo_bt_group_2d_wide(bt_eta, bt_ubt, bt_vbt, grid%nghost)
               call profiler_stop("ocean_comms_bt")
            end if
         else
            call profiler_start("ocean_comms_bt")
            !$acc wait(1)
            call ocean_halo_bt_group_2d(bt_eta, bt_ubt, bt_vbt)
            call profiler_stop("ocean_comms_bt")
         end if
      end do

      ! Snapshot end-of-loop η/u/v BEFORE the time-mean overwrite.
      ! See barotropic_substep_linear for the Hallberg 2009 rationale.
      ! Stays on async(1) from the fast loop; one drain after the time-mean.
      inv_n = 1.0_wp/real(n_steps, wp)
      !$acc kernels async(1)
      do concurrent(j=1:ny, i=1:nx)
         bt_eta_end(i, j) = bt_eta(i, j)
      end do
      do concurrent(j=1:ny, i=1:nx + 1)
         bt_ubt_end(i, j) = bt_ubt(i, j)
      end do
      do concurrent(j=1:ny + 1, i=1:nx)
         bt_vbt_end(i, j) = bt_vbt(i, j)
      end do

      do concurrent(j=1:ny, i=1:nx)
         bt_eta(i, j) = eta_sum(i, j)*inv_n
      end do
      do concurrent(j=1:ny, i=1:nx + 1)
         bt_ubt(i, j) = ubt_sum(i, j)*inv_n
         bt_uhbt(i, j) = uhbt_sum(i, j)*inv_n
      end do
      do concurrent(j=1:ny + 1, i=1:nx)
         bt_vbt(i, j) = vbt_sum(i, j)*inv_n
         bt_vhbt(i, j) = vhbt_sum(i, j)*inv_n
      end do
      !$acc end kernels
      !$acc wait(1)   ! single drain: the whole substep ran on async(1)
      !$acc end data

      ! Refresh the persistent Chapman state from the end-of-step η.
      ! Next outer-step's substep block reads this as `eta_old` when
      ! computing the radiation target.  All host-side: the scalar
      ! `bc%eta_old_chapman_*` is consumed by host code at the top of
      ! the next substep-block call and never read from inside a `do
      ! concurrent` (where it would need a device copy).  Pull
      ! `bt_eta_end` to host so the reduction reads fresh values.
      if (present(bc)) then
         if (bc_w == OBC_CHAPMAN .or. bc_e == OBC_CHAPMAN .or. &
             bc_s == OBC_CHAPMAN .or. bc_n == OBC_CHAPMAN) then
            !$acc update self(bt_eta_end)
            if (bc_w == OBC_CHAPMAN) then
               bc%eta_old_chapman_w = sum(bt_eta_end(i_w_int, &
                                                     grid%nghost + 1:grid%nghost + grid%ny_phys))/ &
                                      real(grid%ny_phys, wp)
            end if
            if (bc_e == OBC_CHAPMAN) then
               bc%eta_old_chapman_e = sum(bt_eta_end(i_e_int, &
                                                     grid%nghost + 1:grid%nghost + grid%ny_phys))/ &
                                      real(grid%ny_phys, wp)
            end if
            if (bc_s == OBC_CHAPMAN) then
               bc%eta_old_chapman_s = sum(bt_eta_end(grid%nghost + 1:grid%nghost + grid%nx_phys, &
                                                     j_s_int))/ &
                                      real(grid%nx_phys, wp)
            end if
            if (bc_n == OBC_CHAPMAN) then
               bc%eta_old_chapman_n = sum(bt_eta_end(grid%nghost + 1:grid%nghost + grid%nx_phys, &
                                                     j_n_int))/ &
                                      real(grid%nx_phys, wp)
            end if
         end if
      end if
   end subroutine barotropic_substep_nonlinear

   subroutine barotropic_substep_nonlinear_interior(grid, metrics, bt_work, f_corner, &
                                                    n_steps, dt_inner, bc, t, eta_forcing)
      !! Interior (normal-width) entry point for the nonlinear barotropic
      !! fast loop.  Unpacks the `bt_work` fast-loop arrays and forwards them
      !! to `barotropic_substep_nonlinear` (bt_halo=0), so the ~20-array
      !! plumbing lives here once instead of at every call site.  The optional
      !! `bc`/`t`/`eta_forcing` propagate by absence (F2018 15.5.2.13), so this
      !! single entry reproduces the former present()-branch call variants
      !! bit-for-bit.  See `bt_wide_substep` for the wide-halo march-in twin.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(barotropic_workstate_t), intent(inout) :: bt_work
      real(wp), intent(in) :: f_corner(grid%nx_total + 1, grid%ny_total + 1)
         !! Coriolis at corners (`coriolis_adv_t%f_corner`); forwarded read-only.
      integer, intent(in) :: n_steps
      real(wp), intent(in) :: dt_inner
      type(ocean_bc_state_t), intent(inout), optional :: bc
      real(wp), intent(in), optional :: t
      real(wp), intent(in), optional :: eta_forcing(grid%nx_total, grid%ny_total)

      call barotropic_substep_nonlinear(grid, bt_work, &
                                        bt_work%F_bt_u_fast, bt_work%F_bt_v_fast, &
                                        n_steps, dt_inner, &
                                        bt_eta=bt_work%bt_eta, bt_H_ref=bt_work%bt_H_ref, &
                                        bt_eta_new=bt_work%bt_eta_new, bt_ke_centre=bt_work%bt_ke_centre, &
                                        eta_sum=bt_work%eta_sum, bt_eta_end=bt_work%bt_eta_end, &
                                        bt_ubt=bt_work%bt_ubt, bt_ubt_prev=bt_work%bt_ubt_prev, &
                                        bt_rem_u=bt_work%bt_rem_u, ubt_sum=bt_work%ubt_sum, &
                                        uhbt_sum=bt_work%uhbt_sum, bt_uhbt=bt_work%bt_uhbt, &
                                        bt_ubt_end=bt_work%bt_ubt_end, &
                                        bt_vbt=bt_work%bt_vbt, bt_vbt_prev=bt_work%bt_vbt_prev, &
                                        bt_rem_v=bt_work%bt_rem_v, vbt_sum=bt_work%vbt_sum, &
                                        vhbt_sum=bt_work%vhbt_sum, bt_vhbt=bt_work%bt_vhbt, &
                                        bt_vbt_end=bt_work%bt_vbt_end, &
                                        bt_zeta_corner=bt_work%bt_zeta_corner, &
                                        f_corner=f_corner, &
                                        area_cu=metrics%areaCu, area_cv=metrics%areaCv, &
                                        dx_cu=metrics%dxCu, dx_cv=metrics%dx_cv_bt, &
                                        dy_cu=metrics%dy_cu_bt, dy_cv=metrics%dyCv, &
                                        iarea_bu=metrics%iareaBu, iarea_t=metrics%iareaT, &
                                        idx_cu=metrics%idxCu, idy_cv=metrics%idyCv, &
                                        bc=bc, t=t, eta_forcing=eta_forcing)
   end subroutine barotropic_substep_nonlinear_interior

end module rdb_barotropic_substep
