!! ALE remap orchestrator for the ocean dynamical core (centre-cell pass).
module rdb_ocean_remap
   !! Drives the conservative vertical remap of `multilayer.h_layer` + every
   !! registered `multilayer.tracers(t)%c` from the current (Lagrangian) grid
   !! to `vcoord%target_h`. The kernel is the per-column `remap_column` from
   !! `rdb_remap_column` (shared with coastal); this module wires it across the
   !! registered tracer slot list plus the face-velocity pass.
   !!
   !! Per-column conservation: sum_k(h*c) is preserved to round-off (the
   !! content the bounding moves and the few-ulp sliver of a column-total
   !! mismatch land in the remap budget slot, which is how they are measured).
   !!
   !! **The tracer remap computes on the stored CONCENTRATION** (White &
   !! Adcroft 2008; Griffies et al. 2020).  The PPM reconstruction (Colella &
   !! Woodward 1984) and the overlap integral run on `c`, giving `c_new`
   !! directly; a vanished source layer is just a cell value in the stencil
   !! (it keeps its last `c`); the target is bounded by the source cells that
   !! contributed, and an empty target takes the mean of the source cell
   !! containing its bottom interface.  There are no read/write pool merges
   !! and no defect fold around the kernel.
#ifdef LFORTRAN_PASSING
   use rdb_constants, only: wp, REMAP_PPM, H_VANISHED, H_DIV_EPS
#else
   use rdb_constants, only: NZ_STACK_MAX, wp, REMAP_PPM, H_VANISHED, H_DIV_EPS
#endif
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_grid, only: hgrid_t
   use rdb_remap_column, only: remap_column
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_tracer, only: TRACER_BUDGET_HEAT, TRACER_BUDGET_SALT
   use rdb_ocean_vcoord, only: ocean_vcoord_t, VCOORD_EULERIAN_Z, VCOORD_LAGRANGIAN, &
                               VCOORD_RHO, VCOORD_HYCOM
   use rdb_eos, only: eos_t
   implicit none
   private

   ! Degenerate-COLUMN guard for `rescale_anomaly_ke` — a face column whose
   ! whole depth is at or below the vanish marker has no barotropic mean to
   ! rescale about.  Aliased to `H_VANISHED` (same value) so there is ONE
   ! definition of "vanished" in the tree; it used to be a bare `1.5e-4_wp`
   ! literal here, which is a third definition waiting to drift from the
   ! constant of record.
   !
   ! Every test of the marker is a STRICT `>`: a layer sitting exactly ON it
   ! reads as vanished, which is what the geometric vcoord families rely on
   ! (see `rdb_vcoord :: vcoord_h_min_role`).
   real(wp), parameter :: H_FLOOR = H_VANISHED

   real(wp), parameter, public :: OCEAN_REMAP_PRECOND_RTOL = 1.0e-9_wp
      !! Relative tolerance on `sum(h_old) == sum(target_h)` for the
      !! precondition assertion.  The target builders reach the sum by a
      !! different arithmetic route than the continuity update does, so a
      !! few ulp of drift per layer is expected and is not the failure mode
      !! being hunted: the violations that matter (a degenerate column that
      !! manufactures `nz*h_min` of water, a short target that deletes its
      !! tail) are percent-level, decades above this.

   public :: ocean_apply_ale_remap_centres
   public :: ocean_remap_scan_preconditions
   public :: ocean_apply_ale_remap_faces
   public :: ocean_apply_ale_remap_step
   public :: ocean_remap_tracer_column   ! exposed for unit tests
   public :: ocean_remap_tracer_field   ! reference kernel, exposed for the stream test
   public :: ocean_remap_tracer_field_stream
   public :: remap_use_stream

#ifdef LFORTRAN_PASSING
   integer, parameter :: NZ_STACK_MAX = 64
      !! LFortran 0.64 workaround: module-local copy of the rdb_constants value
      !! (an imported parameter used as an explicit-shape dummy bound inside a
      !! PURE call becomes an impure getter under LFortran). Keep in sync (=64).
#endif

contains

   pure subroutine ocean_apply_ale_remap_centres(grid, vcoord, ms, bt_eta, bt_H_ref, method, eos, dt)
      !! Orchestrate the centre-cell pass of the ALE remap step (h_layer +
      !! tracers). Public only for the unit-test suite.
      !! Sequence: skip if EULERIAN_Z/LAGRANGIAN; snapshot column total + h_old;
      !! build `vcoord%target_h`; remap each tracer h_old→target_h via the PPM
      !! column kernel; set h_layer = target_h; recompute
      !! bt_eta = sum_k(h_layer) - bt_H_ref. Face velocities remapped separately.
      !! Takes bt_eta/bt_H_ref directly (not ocean_dyn_t) so this module sits
      !! below the split driver in the dependency tree.
      type(hgrid_t), intent(in) :: grid
      type(ocean_vcoord_t), intent(inout) :: vcoord
      type(multilayer_state_t), intent(inout) :: ms
      ! assumed-shape-ok: outer-driver allocatable; grid%nx_total available; deferred.
      real(wp), intent(inout) :: bt_eta(:, :)
         !! Free-surface anomaly η at cell centres (m).  Updated to
         !! `sum_k(h_layer) - bt_H_ref` on exit (step 7).
      real(wp), intent(in) :: bt_H_ref(:, :)  ! assumed-shape-ok: same as bt_eta above
         !! Reference column depth H (m, positive-down).  Constant for
         !! the ocean path; passed in so the routine doesn't need a
         !! handle on the split-driver state.
      type(eos_t), intent(in), optional :: eos
         !! Device-resident EOS handle — required only for `VCOORD_RHO`
         !! (isopycnal density inversion); ignored by the geometric coords.
      real(wp), intent(in), optional :: dt
         !! Outer/thermo timestep (s) for the grid time-filter.  Absent or
         !! `regrid_time_scale = 0` (default) ⇒ filter skipped ⇒
         !! bit-identical.  See `ocean_apply_ale_remap_step`.
      integer, intent(in), optional :: method
         !! REMAP_PCM / REMAP_PLM / REMAP_PPM. Defaults to PPM (parabolic stencil
         !! cuts the spurious vertical mixing a 1st-order limiter introduces).

      integer :: t, m, nx, ny, nz, i, j, k
      real(wp) :: wtd
      logical :: do_tfilter, use_stream

      ! Eulerian-z and Lagrangian/isopycnal both skip remap: the former
      ! holds h at H·dsig via vert-advection cancellation, the latter
      ! lets h evolve freely (target = current h).
      if (vcoord%coord_type == VCOORD_EULERIAN_Z .or. &
          vcoord%coord_type == VCOORD_LAGRANGIAN) return
      if (.not. vcoord%is_init) return
      if (.not. allocated(ms%h_layer)) return

      m = REMAP_PPM
      if (present(method)) m = method
      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml

      ! --- 2. Current column total (persistent device-mapped scratch on vcoord)
      do concurrent(j=1:ny, i=1:nx) local(k)
         vcoord%remap_total_h(i, j) = 0.0_wp
         do k = 1, nz
            vcoord%remap_total_h(i, j) = vcoord%remap_total_h(i, j) &
                                         + ms%h_layer(i, j, k)
         end do
      end do

      ! --- 3. Snapshot h_old on-device (host source= reads stale memory after a
      ! dynamics-side OpenACC kernel). Before target_h so VCOORD_RHO can read it.
      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         vcoord%remap_h_old(i, j, k) = ms%h_layer(i, j, k)
      end do

      ! --- 4. Populate target_h (geometric vs isopycnal dispatch) ---
      do concurrent(j=1:ny, i=1:nx)
         vcoord%remap_h_ref(i, j) = vcoord%remap_total_h(i, j) - bt_eta(i, j)
      end do
      if ((vcoord%coord_type == VCOORD_RHO .or. vcoord%coord_type == VCOORD_HYCOM) &
          .and. present(eos) .and. allocated(ms%tracers) &
          .and. ms%idx_temperature > 0 .and. ms%idx_salinity > 0) then
         call build_ts_concentration(nx, ny, nz, vcoord%remap_h_old, &
                                     ms%tracers(ms%idx_temperature)%c, &
                                     ms%tracers(ms%idx_salinity)%c, &
                                     vcoord%remap_conc_t, vcoord%remap_conc_s)
         ! HYCOM = RHO inversion + z*-floor/monotonize deltas (hybrid=.true.);
         ! pure RHO passes hybrid=.false. (bit-identical to RHO-only kernel).
         call vcoord%compute_target_h_rho(vcoord%remap_h_ref, bt_eta, &
                                          vcoord%remap_conc_t, vcoord%remap_conc_s, eos, &
                                          hybrid=(vcoord%coord_type == VCOORD_HYCOM))
      else
         call vcoord%compute_target_h(vcoord%remap_h_ref, bt_eta)
      end if

      ! --- 4b. Grid time-filter (White & Adcroft 2008): relax target toward it
      ! from the old grid by wtd = dt/(τ+dt). Convex blend ⇒ column total
      ! conserved. τ=0 (default) or dt absent ⇒ skipped, bit-identical.
      do_tfilter = vcoord%regrid_time_scale > 0.0_wp .and. present(dt)
      if (do_tfilter) then
         wtd = dt/(vcoord%regrid_time_scale + dt)
         do concurrent(k=1:nz, j=1:ny, i=1:nx)
            vcoord%target_h(i, j, k) = vcoord%remap_h_old(i, j, k) &
                                       + wtd*(vcoord%target_h(i, j, k) - vcoord%remap_h_old(i, j, k))
         end do
      end if

      ! --- 5. Remap every tracer ---
      use_stream = remap_use_stream(vcoord, nz, m, vcoord%remap_boundary_extrap, &
                                    vcoord%remap_nonuniform_weights)
      if (allocated(ms%tracers)) then
         do t = 1, size(ms%tracers)
            if (.not. allocated(ms%tracers(t)%c)) cycle
            select case (ms%tracers(t)%budget_id)
            case (TRACER_BUDGET_HEAT)
               if (use_stream) then
                  call ocean_remap_tracer_field_stream( &
                     nx, ny, nz, vcoord%remap_h_old, vcoord%target_h, ms%tracers(t)%c, &
                     vcoord%remap_c_old, vcoord%remap_z_old, vcoord%remap_z_new, &
                     vcoord%remap_q_l, vcoord%remap_q_r, vcoord%remap_q_6, &
                     budget=ms%heat_budget_remap)
               else
                  call ocean_remap_tracer_field( &
                     nx, ny, nz, vcoord%remap_h_old, vcoord%target_h, ms%tracers(t)%c, m, &
                     vcoord%remap_boundary_extrap, vcoord%remap_nonuniform_weights, &
                     budget=ms%heat_budget_remap)
               end if
            case (TRACER_BUDGET_SALT)
               if (use_stream) then
                  call ocean_remap_tracer_field_stream( &
                     nx, ny, nz, vcoord%remap_h_old, vcoord%target_h, ms%tracers(t)%c, &
                     vcoord%remap_c_old, vcoord%remap_z_old, vcoord%remap_z_new, &
                     vcoord%remap_q_l, vcoord%remap_q_r, vcoord%remap_q_6, &
                     budget=ms%salt_budget_remap)
               else
                  call ocean_remap_tracer_field( &
                     nx, ny, nz, vcoord%remap_h_old, vcoord%target_h, ms%tracers(t)%c, m, &
                     vcoord%remap_boundary_extrap, vcoord%remap_nonuniform_weights, &
                     budget=ms%salt_budget_remap)
               end if
            case default
               if (use_stream) then
                  call ocean_remap_tracer_field_stream( &
                     nx, ny, nz, vcoord%remap_h_old, vcoord%target_h, ms%tracers(t)%c, &
                     vcoord%remap_c_old, vcoord%remap_z_old, vcoord%remap_z_new, &
                     vcoord%remap_q_l, vcoord%remap_q_r, vcoord%remap_q_6)
               else
                  call ocean_remap_tracer_field( &
                     nx, ny, nz, vcoord%remap_h_old, vcoord%target_h, ms%tracers(t)%c, m, &
                     vcoord%remap_boundary_extrap, vcoord%remap_nonuniform_weights)
               end if
            end select
         end do
      end if

      ! --- 6. h_layer = target_h; capture mass-budget delta ---
      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         ms%mass_budget_remap(i, j, k) = ms%mass_budget_remap(i, j, k) &
                                         + (vcoord%target_h(i, j, k) - vcoord%remap_h_old(i, j, k))
         ms%h_layer(i, j, k) = vcoord%target_h(i, j, k)
      end do

      ! --- 7. Recompute bt_eta from new sum (round-off-tight) ---
      do concurrent(j=1:ny, i=1:nx) local(k)
         bt_eta(i, j) = -bt_H_ref(i, j)
         do k = 1, nz
            bt_eta(i, j) = bt_eta(i, j) + ms%h_layer(i, j, k)
         end do
      end do
   end subroutine ocean_apply_ale_remap_centres

   pure subroutine ocean_remap_tracer_field(nx, ny, nz, h_old, h_new, c, method, &
                                            bnd_extrap, nonunif, budget, c_new_out, clamp_out)
      !! Flat-impl REFERENCE tracer remap on the stored concentration. Per
      !! (i,j) column: donor-substituted `c_don` (`rdb_vl_column_conc`) →
      !! per-column remap kernel on `c_don` → each filler's excess
      !! `c_old − c_don` added back over its overlap
      !! (`remap_add_filler_excess`) → bounding + empty-target fallback
      !! (`remap_bound_column`) → `c = c_new`.  A filler's stale value never
      !! enters a live cell's reconstruction stencil, and the column content
      !! `sum(h_old·c_old)` is conserved exactly.
      !! `budget` (optional): the per-cell `h_new·c_new − h_old·c_old`
      !! increment (two standalone products) is accumulated into the slot
      !! (heat/salt remap deltas) before overwriting.
      !! Flat-arg so GPU codegen doesn't chase the array-of-derived-types pointer.
      !!
      !! **Concentration form** (White & Adcroft 2008; Griffies et al. 2020).
      !! The reconstruction is built on `c_old` and integrated over the
      !! overlap of old and new layers; a source filler reads its donor in the
      !! stencil and carries its own excess as a uniform offset, so there is NO
      !! read-side pool merge, NO write-side pool merge and NO defect fold
      !! around the kernel.  The target concentration is bounded by the range
      !! of the source cells that contributed to it, and a target with
      !! `h_new <= 0` takes the mean of the source cell containing its bottom
      !! interface; both are in `remap_bound_column`.  The content the bounding
      !! moves, `h_new·(c_bounded − c_raw)`, is part of `budget`:
      !! conservation is MEASURED by the slot, not assumed.
      !!
      !! **The `(h_old, h_new)` pair is the call boundary**: every caller
      !! passes `vcoord%remap_h_old` and `vcoord%target_h` (the thickness the
      !! regrid wrote into `ms%h_layer` right after).
      integer, intent(in) :: nx, ny, nz, method
      real(wp), intent(in) :: h_old(nx, ny, nz)
      real(wp), intent(in) :: h_new(nx, ny, nz)
      real(wp), intent(inout) :: c(nx, ny, nz)
         !! Tracer concentration: read as `c_old`, overwritten with `c_new`.
      logical, intent(in) :: bnd_extrap
         !! Linear-exact boundary-cell reconstruction in the column kernel
         !! (`&vcoord_nml remap_boundary_extrap`); `.false.` ⇒ the PCM
         !! flatten.
      logical, intent(in) :: nonunif
         !! Non-uniform-grid PLM/PPM reconstruction weights in the column
         !! kernel (`&vcoord_nml remap_nonuniform_weights`); `.false.` ⇒ the
         !! equal-thickness specialisations.
      real(wp), intent(inout), optional :: budget(nx, ny, nz)
      real(wp), intent(inout), optional :: c_new_out(nx, ny, nz)
         !! Diagnostic / test hook: the bounded target concentration.
      real(wp), intent(inout), optional :: clamp_out(nx, ny, nz)
         !! Diagnostic / test hook: the content the bounding moved,
         !! `h_new·(c_bounded − c_raw)`.
      integer :: i, j, k
      real(wp) :: h_old_col(NZ_STACK_MAX), h_new_col(NZ_STACK_MAX)
      real(wp) :: c_old_col(NZ_STACK_MAX), c_new_col(NZ_STACK_MAX), c_don_col(NZ_STACK_MAX)
      real(wp) :: clamp_col(NZ_STACK_MAX)
      real(wp) :: p_new, p_old

      ! Gate the budget write INSIDE the one loop (as vdiff does); splitting
      ! present(budget) into two loops makes NVHPC compile the no-budget branch
      ! ~25x slower.
      do concurrent(j=1:ny, i=1:nx) &
         local(k, h_old_col, h_new_col, c_old_col, c_new_col, c_don_col, clamp_col, p_new, p_old)
         do k = 1, nz
            h_old_col(k) = h_old(i, j, k)
            h_new_col(k) = h_new(i, j, k)
            c_old_col(k) = c(i, j, k)
         end do
         call rdb_vl_column_conc(nz, h_old_col, c_old_col, c_don_col)
         call remap_column(method, nz, &
                           h_old_col(1:nz), h_new_col(1:nz), &
                           c_don_col(1:nz), c_new_col(1:nz), bnd_extrap, nonunif)
         call remap_add_filler_excess(nz, h_old_col, h_new_col, c_don_col, c_old_col, c_new_col)
         call remap_bound_column(nz, h_old_col, h_new_col, c_don_col, c_old_col, c_new_col, &
                                 clamp_col, bnd_extrap)
         do k = 1, nz
            if (present(budget)) then
               p_new = h_new_col(k)*c_new_col(k)
               p_old = h_old_col(k)*c_old_col(k)
               budget(i, j, k) = budget(i, j, k) + (p_new - p_old)
            end if
            if (present(c_new_out)) c_new_out(i, j, k) = c_new_col(k)
            if (present(clamp_out)) clamp_out(i, j, k) = clamp_col(k)
            c(i, j, k) = c_new_col(k)
         end do
      end do
   end subroutine ocean_remap_tracer_field

   pure subroutine remap_bound_column(nz, h_old_col, h_new_col, c_old, c_own, c_new, clamp_amt, &
                                      bnd_extrap)
      !$acc routine seq
      !! Concentration-form closure of one remapped column (reference path;
      !! the stream kernel 3 does the same per cell).  For each target layer
      !! `k`, with the source cells that CONTRIBUTE to it (`overlap > 0` and
      !! `h_old > 0`, exactly the cells the overlap integral adds up):
      !!  * `h_new(k) > 0` and at least one contributor: clamp `c_new(k)` to
      !!    `[min, max]` of `c_old` over the contributors' reconstruction
      !!    stencils (each contributor `ko` with `ko-1` and `ko+1`) — the
      !!    bounding of White & Adcroft (2008) / Griffies et al. (2020).
      !!    Analytically a no-op for a Colella & Woodward (1984) limited
      !!    parabola, which cannot leave its stencil's range (bounding by the
      !!    contributor MEANS alone would flatten every target inside one
      !!    source cell to PCM); it removes the
      !!    `O(eps·|c|·h_src/h_new)` overshoot a few-ulp-thin target can pick
      !!    up in the quotient.  A non-finite raw value is left alone (an
      !!    `if/else` clamp would launder it into a plausible `lo`/`hi`).
      !!  * `h_new(k) <= 0`, or no contributor (a sliver above the source
      !!    column): `c_new(k)` is the mean of the source cell containing the
      !!    target's bottom interface (the topmost cell if it lies above the
      !!    column) — the water that interface sits in, never a made-up `0`.
      !! `clamp_amt(k) = h_new(k)·(c_new_final − c_new_raw)` is the content
      !! the closure moves; it is part of `h_new·c_new − h_old·c_old`, hence of the
      !! remap budget slot.
      !!
      !! **Boundary cells under `bnd_extrap`.**  The linear-exact boundary
      !! closure (`boundary_half_jump` in `rdb_remap_column`) gives cells `1`
      !! and `nz` the straight line `c ± d`, `|d| <= |c_nbr − c|`, whose OUTER
      !! half lies outside `[c_nbr, c]` by construction — that is the point of
      !! the closure (a profile linear in z is reproduced exactly).  Bounding
      !! such a contributor by its neighbour means alone clamped every target
      !! inside that outer half (a 20 m surface target under a 67 m source
      !! cell: the whole linear-S column gained 1e-5 of its salt on the first
      !! remap of the flat at-rest control, and the clamped bed/surface
      !! density set rx0 slopes in motion).  Its range is therefore widened by
      !! the mirror point `2c − c_nbr`, which contains the line for every
      !! method that uses the closure (PLM/PPM/PPM_H4/PQM), so the bounding is
      !! again a no-op on an exact reconstruction.
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_old_col(NZ_STACK_MAX)
      real(wp), intent(in) :: h_new_col(NZ_STACK_MAX)
      real(wp), intent(in) :: c_old(NZ_STACK_MAX)
         !! Donor-substituted source concentration (`rdb_vl_column_conc`):
         !! the column the reconstruction was built on.
      real(wp), intent(in) :: c_own(NZ_STACK_MAX)
         !! STORED source concentration.  A filler contributor's profile is its
         !! donor parabola shifted by `c_own − c_old`, so its range is shifted
         !! the same; the empty-target fallback reads the stored value.
      real(wp), intent(inout) :: c_new(NZ_STACK_MAX)
         !! Raw target concentration in, bounded concentration out.
      real(wp), intent(out) :: clamp_amt(NZ_STACK_MAX)
      logical, intent(in) :: bnd_extrap
         !! The column kernel ran with the linear-exact boundary closure
         !! (`&vcoord_nml remap_boundary_extrap`): widen the range of the two
         !! boundary cells to the line they carry (see above).
      real(wp) :: z_old(NZ_STACK_MAX), z_new(NZ_STACK_MAX)
      real(wp) :: zacc_o, zacc_n, zn_lo, zn_hi, zo_lo, zo_hi, overlap
      real(wp) :: c_lo, c_hi, c_raw, cs_lo, cs_hi
      integer :: k, ko, ko0, n_hit

      zacc_o = 0.0_wp
      zacc_n = 0.0_wp
      do k = 1, nz
         zacc_o = zacc_o + h_old_col(k)
         zacc_n = zacc_n + h_new_col(k)
         z_old(k) = zacc_o
         z_new(k) = zacc_n
      end do

      ko0 = 1
      do k = 1, nz
         zn_hi = z_new(k)
         zn_lo = 0.0_wp
         if (k > 1) zn_lo = z_new(k - 1)
         ! First source cell whose top lies above the target's bottom (z_old
         ! is non-decreasing, so the pointer only moves up with k).
         do while (ko0 <= nz)
            if (z_old(ko0) > zn_lo) exit
            ko0 = ko0 + 1
         end do
         c_raw = c_new(k)
         n_hit = 0
         c_lo = 0.0_wp
         c_hi = 0.0_wp
         if (h_new_col(k) > 0.0_wp) then
            ko = ko0
            zo_lo = 0.0_wp
            if (ko > 1 .and. ko <= nz) zo_lo = z_old(ko - 1)
            do while (ko <= nz)
               if (zo_lo >= zn_hi) exit
               zo_hi = z_old(ko)
               overlap = min(zn_hi, zo_hi) - max(zn_lo, zo_lo)
               if (overlap > 0.0_wp .and. h_old_col(ko) > 0.0_wp) then
                  ! Range of the contributor's reconstruction stencil
                  ! (ko-1, ko, ko+1): the limited parabola cannot leave it.
                  cs_lo = min(c_old(max(ko - 1, 1)), min(c_old(ko), c_old(min(ko + 1, nz))))
                  cs_hi = max(c_old(max(ko - 1, 1)), max(c_old(ko), c_old(min(ko + 1, nz))))
                  if (bnd_extrap .and. nz >= 2) then
                     if (ko == 1) then
                        cs_lo = min(cs_lo, 2.0_wp*c_old(1) - c_old(2))
                        cs_hi = max(cs_hi, 2.0_wp*c_old(1) - c_old(2))
                     else if (ko == nz) then
                        cs_lo = min(cs_lo, 2.0_wp*c_old(nz) - c_old(nz - 1))
                        cs_hi = max(cs_hi, 2.0_wp*c_old(nz) - c_old(nz - 1))
                     end if
                  end if
                  if (h_old_col(ko) <= H_VANISHED) then
                     cs_lo = cs_lo + (c_own(ko) - c_old(ko))
                     cs_hi = cs_hi + (c_own(ko) - c_old(ko))
                  end if
                  if (n_hit == 0) then
                     c_lo = cs_lo
                     c_hi = cs_hi
                  else
                     c_lo = min(c_lo, cs_lo)
                     c_hi = max(c_hi, cs_hi)
                  end if
                  n_hit = n_hit + 1
               end if
               zo_lo = zo_hi
               ko = ko + 1
            end do
         end if
         if (n_hit > 0) then
            if (ieee_is_finite(c_raw)) c_new(k) = max(c_lo, min(c_hi, c_raw))
         else
            c_new(k) = c_own(min(ko0, nz))
         end if
         clamp_amt(k) = h_new_col(k)*(c_new(k) - c_raw)
      end do
   end subroutine remap_bound_column

   pure subroutine remap_add_filler_excess(nz, h_old_col, h_new_col, c_don, c_own, c_new)
      !$acc routine seq
      !! Reference-path half of the filler rule (the stream kernel 2 shifts the
      !! filler's parabola instead).  The column kernel ran on the
      !! donor-substituted `c_don`; each source filler `m` (`0 < h_old <=
      !! H_VANISHED`) owns `h_old(m)·(c_own(m) − c_don(m))` more, carried as a
      !! uniform offset over its interval: every target `k` with `h_new > 0`
      !! gains `s·overlap/h_new(k)`.  Exactly restores `sum(h_old·c_own)`.
      !! A column without fillers is untouched (bit-identical).
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_old_col(NZ_STACK_MAX), h_new_col(NZ_STACK_MAX)
      real(wp), intent(in) :: c_don(NZ_STACK_MAX), c_own(NZ_STACK_MAX)
      real(wp), intent(inout) :: c_new(NZ_STACK_MAX)
      real(wp) :: zo_lo, zo_hi, zn_lo, zn_hi, overlap, s
      integer :: k, m

      zo_lo = 0.0_wp
      do m = 1, nz
         zo_hi = zo_lo + h_old_col(m)
         if (h_old_col(m) > 0.0_wp .and. h_old_col(m) <= H_VANISHED) then
            s = c_own(m) - c_don(m)
            if (s /= 0.0_wp) then
               zn_lo = 0.0_wp
               do k = 1, nz
                  zn_hi = zn_lo + h_new_col(k)
                  overlap = min(zn_hi, zo_hi) - max(zn_lo, zo_lo)
                  if (overlap > 0.0_wp .and. h_new_col(k) > 0.0_wp) then
                     c_new(k) = c_new(k) + s*overlap/h_new_col(k)
                  end if
                  zn_lo = zn_hi
               end do
            end if
         end if
         zo_lo = zo_hi
      end do
   end subroutine remap_add_filler_excess

   pure function remap_use_stream(vcoord, nz, method, bnd_extrap, nonunif) result(use_stream)
      !! Host-side dispatch of the tracer remap: `.true.` when the streaming
      !! kernel (`ocean_remap_tracer_field_stream`) can and should run.  The
      !! stream implements exactly the default PPM column (`REMAP_PPM`,
      !! `nz >= 3`, uniform-weight edges, PCM boundary cells) and needs its
      !! global scratch buffers; everything else (other methods,
      !! `nz <= 2`, `remap_boundary_extrap`, `remap_nonuniform_weights`,
      !! `&vcoord_nml remap_kernel = "reference"`) takes the stack-based
      !! reference kernel.  The stream kernel itself contains no call into
      !! stack-based code, hence the decision is taken here, not per column.
      type(ocean_vcoord_t), intent(in) :: vcoord
      integer, intent(in) :: nz, method
      logical, intent(in) :: bnd_extrap, nonunif
      logical :: use_stream
      use_stream = vcoord%remap_stream .and. method == REMAP_PPM .and. nz >= 3 &
                   .and. .not. bnd_extrap .and. .not. nonunif &
                   .and. allocated(vcoord%remap_c_old) &
                   .and. allocated(vcoord%remap_z_old) &
                   .and. allocated(vcoord%remap_z_new) &
                   .and. allocated(vcoord%remap_q_l) &
                   .and. allocated(vcoord%remap_q_r) &
                   .and. allocated(vcoord%remap_q_6)
   end function remap_use_stream

   pure function remap_ppm_edge(ko, nz, c_m1, c_0, c_p1, c_p2) result(e_hi)
      !$acc routine seq
      !! Raw (unlimited) PPM edge between source cells `ko` and `ko+1` from
      !! the four-cell window `c(ko-1..ko+2)`: Colella & Woodward (1984) eq.
      !! 1.9 for uniform cells, the two-point mean next to either boundary.
      integer, intent(in) :: ko, nz
      real(wp), intent(in) :: c_m1, c_0, c_p1, c_p2
      real(wp) :: e_hi
      if (ko == nz) then
         e_hi = c_0
      else if (ko == 1 .or. ko == nz - 1) then
         e_hi = 0.5_wp*(c_0 + c_p1)
      else
         e_hi = (7.0_wp/12.0_wp)*(c_0 + c_p1) &
                - (1.0_wp/12.0_wp)*(c_m1 + c_p2)
      end if
   end function remap_ppm_edge

   pure function remap_c_at(nx, ny, nz, i, j, kk, c_old) result(c)
      !$acc routine seq
      !! Source concentration of layer `kk` of column `(i,j)`, `0` outside
      !! `1..nz`.
      integer, intent(in) :: nx, ny, nz, i, j, kk
      real(wp), intent(in) :: c_old(nx, ny, nz)
      real(wp) :: c
      if (kk < 1 .or. kk > nz) then
         c = 0.0_wp
      else
         c = c_old(i, j, kk)
      end if
   end function remap_c_at

   pure function remap_c_don_at(nx, ny, nz, i, j, kk, c_old, h_old) result(c)
      !$acc routine seq
      !! DONOR-substituted source concentration of layer `kk` of column
      !! `(i,j)` (`0` outside `1..nz`): the stored value on a live layer; on a
      !! filler (`h_old <= H_VANISHED`) the nearest live layer above it, or the
      !! topmost live layer for a run reaching the top; `0` in a column with
      !! no live layer.  The per-layer form of `rdb_vl_column_conc`, so the
      !! stream kernels need no extra scratch.  A filler's stale value must not
      !! enter a LIVE neighbour's PPM edge stencil.
      integer, intent(in) :: nx, ny, nz, i, j, kk
      real(wp), intent(in) :: c_old(nx, ny, nz), h_old(nx, ny, nz)
      real(wp) :: c
      integer :: m, m_don
      ! Kernel 2 only.  Called from kernel 3, nvfortran's device code was
      ! wrong by O(1) even on filler-free columns (host builds correct), so
      ! kernel 3 carries its own inline walk.
      c = 0.0_wp
      if (kk >= 1 .and. kk <= nz) then
         if (h_old(i, j, kk) > H_VANISHED) then
            c = c_old(i, j, kk)
         else
            ! Nearest live layer above: the last live one seen walking down to kk+1.
            m_don = 0
            do m = nz, kk + 1, -1
               if (h_old(i, j, m) > H_VANISHED) m_don = m
            end do
            ! None above: the topmost live layer below (last seen walking up).
            if (m_don == 0) then
               do m = 1, kk - 1
                  if (h_old(i, j, m) > H_VANISHED) m_don = m
               end do
            end if
            if (m_don > 0) c = c_old(i, j, m_don)
         end if
      end if
   end function remap_c_don_at

   pure subroutine remap_ppm_cell(ko, nz, c_m1, c_0, c_p1, c_p2, e_lo, e_hi, q_l, q_r, q_6)
      !$acc routine seq
      !! One source cell of the uniform-weight PPM column, from its window:
      !! the RAW right edge `e_hi` shared with cell `ko+1`, and the Colella &
      !! Woodward (1984) limited parabola `(q_l, q_r, q_6)` of cell `ko`.
      !! Statement-for-statement the edge estimate + limiter of
      !! `remap_column_ppm` (steps 1-2) for one `k`; `e_lo` is the RAW left
      !! edge (`remap_ppm_edge` of cell `ko-1`; ignored for cell 1, whose left
      !! edge is `c_0`), NOT the previous cell's limited right edge.
      integer, intent(in) :: ko, nz
      real(wp), intent(in) :: c_m1, c_0, c_p1, c_p2, e_lo
      real(wp), intent(out) :: e_hi, q_l, q_r, q_6
      real(wp) :: q_min, q_max, dq, dq_l, dq_r

      e_hi = remap_ppm_edge(ko, nz, c_m1, c_0, c_p1, c_p2)

      q_min = c_0
      q_max = c_0
      if (ko > 1) then
         q_min = min(q_min, c_m1)
         q_max = max(q_max, c_m1)
      end if
      if (ko < nz) then
         q_min = min(q_min, c_p1)
         q_max = max(q_max, c_p1)
      end if
      if (ko == 1) then
         q_l = max(q_min, min(q_max, c_0))
      else
         q_l = max(q_min, min(q_max, e_lo))
      end if
      q_r = max(q_min, min(q_max, e_hi))
      dq = q_r - q_l
      dq_l = c_0 - q_l
      dq_r = q_r - c_0
      if (dq_l*dq_r <= 0.0_wp) then
         q_l = c_0
         q_r = c_0
      else
         q_6 = 6.0_wp*c_0 - 3.0_wp*(q_l + q_r)
         if (abs(q_6) > abs(dq)) then
            if (q_6*dq > 0.0_wp) then
               q_l = 3.0_wp*c_0 - 2.0_wp*q_r
            else
               q_r = 3.0_wp*c_0 - 2.0_wp*q_l
            end if
         end if
      end if
      q_6 = 6.0_wp*c_0 - 3.0_wp*(q_l + q_r)
   end subroutine remap_ppm_cell

   pure subroutine ocean_remap_tracer_field_stream(nx, ny, nz, h_old, h_new, c, &
                                                   c_old, z_old, z_new, &
                                                   q_lft, q_rgt, q_six, budget, &
                                                   c_new_out, clamp_out)
      !! Parallel-in-`k` twin of `ocean_remap_tracer_field` for the default
      !! PPM column (`REMAP_PPM`, `nz >= 3`, `bnd_extrap = nonunif = .false.`),
      !! in concentration form: same bits as the reference on a filler-free
      !! column, no `nz`-sized local array, and almost all the work spread
      !! over one thread per CELL instead of one thread per column.  The
      !! caller (`remap_use_stream`) guarantees the preconditions; nothing here
      !! calls stack-based code.
      !!
      !! **Kernels** (all inside this routine; global scratch in between):
      !!  1. `do concurrent(j,i)`, serial `k`: interface prefix sums
      !!     `z_old(k) = z_old(k-1) + h_old(k)`, `z_new` likewise (from 0,
      !!     increasing `k`, as `remap_column_ppm`) and the source
      !!     concentration scratch `c_old = c` per layer.  No pooling, no
      !!     filler pre-scan: a filler is just a cell value.
      !!  2. `do concurrent(k,j,i)`, one thread per SOURCE cell: raw edges
      !!     from the window `c_old(k-2..k+2)`, Colella & Woodward (1984)
      !!     limiter, the limited parabola stored in `q_lft`/`q_rgt`/`q_six`.
      !!  3. `do concurrent(k,j,i)`, one thread per TARGET cell: the overlap
      !!     integral (White & Adcroft 2008) over the contributing source
      !!     cells, `c_new = integral/h_new`, bounded to the range of those
      !!     cells' reconstruction stencils (an empty target, `h_new <= 0`, or one with no
      !!     contributor takes the mean of the source cell containing its
      !!     bottom interface), then the commit: `budget += c_new·h_new − c·h_old`
      !!     and `c = c_new`.  Kernel 3 reads `c` of its OWN cell only (all
      !!     neighbours went through `c_old` in kernel 1), so the commit needs
      !!     no separate pass.
      !!
      !! **Why kernel 3 equals the reference walk.**  The reference computes
      !! each target's `integral` as `0.0` plus, in increasing source index
      !! `ko`, the contribution of every source cell it visits, where a visit
      !! adds nothing unless `overlap = min(z_new(k), z_old(ko)) -
      !! max(z_new(k-1), z_old(ko-1)) > 0` (and `dz_old(ko) > 0`), and its
      !! `exit` / `ko_start` logic only skips cells that cannot overlap.  For
      !! non-negative thicknesses (the remap's precondition, see
      !! `ocean_remap_scan_preconditions`) `z_old` is non-decreasing, so the
      !! cells with `overlap > 0` are a CONTIGUOUS range of `ko`: those with
      !! `z_old(ko) > z_new(k-1)` (a suffix, whose first element
      !! the gallop-and-bisect at the top of kernel 3 finds exactly; exact ties give `overlap <= 0`
      !! and are skipped) and `z_old(ko-1) < z_new(k)` (a prefix, the loop
      !! leaves at the first `z_old(ko-1) >= z_new(k)`).  Kernel 3 therefore
      !! adds the same terms, from `0.0`, in the same increasing order, with
      !! the same contribution expression, and the same contributor set
      !! feeds the bounding as in `remap_bound_column`.  Only inputs
      !! violating the preconditions (negative thickness) can differ.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: h_old(nx, ny, nz)
      real(wp), intent(in) :: h_new(nx, ny, nz)
      real(wp), intent(inout) :: c(nx, ny, nz)
         !! Tracer concentration: read as the source `c_old`, overwritten with `c_new`.
      real(wp), intent(inout) :: c_old(nx, ny, nz)
         !! Scratch: snapshot of the stored source concentration (kernels 2-3
         !! read a filler through its donor, `remap_c_don_at`).
      real(wp), intent(inout) :: z_old(nx, ny, nz)
         !! Scratch: source interface `z_old(k)` (top of layer `k`; the bed is 0).
      real(wp), intent(inout) :: z_new(nx, ny, nz)
         !! Scratch: target interface `z_new(k)`.
      real(wp), intent(inout) :: q_lft(nx, ny, nz)
         !! Scratch: limited left edge of each source cell.
      real(wp), intent(inout) :: q_rgt(nx, ny, nz)
         !! Scratch: limited right edge of each source cell.
      real(wp), intent(inout) :: q_six(nx, ny, nz)
         !! Scratch: parabola curvature `q6` of each source cell.
      real(wp), intent(inout), optional :: budget(nx, ny, nz)
      real(wp), intent(inout), optional :: c_new_out(nx, ny, nz)
         !! Diagnostic / test hook: the bounded target concentration.
      real(wp), intent(inout), optional :: clamp_out(nx, ny, nz)
         !! Diagnostic / test hook: the content the bounding moved,
         !! `h_new·(c_bounded − c_raw)`.
      integer :: i, j, k, ko, lo, hi, mid, stp, k_fb, n_hit, km, kp, m, m_don
      real(wp) :: c_m2, c_m1, c_0, c_p1, c_p2, e_lo, e_hi, q_l, q_r, q_6
      real(wp) :: zo_lo, zo_hi, dzo, zn_lo, zn_hi, dzn
      real(wp) :: z_lo, z_hi, overlap, xi_lo, xi_hi, integral, c_raw, c_new, p_new, p_old
      real(wp) :: c_lo, c_hi, cs_lo, cs_hi, zacc_o, zacc_n, s_fill

      ! ---- kernel 1: interfaces + a snapshot of the source concentration ----
      ! (`c` is overwritten in kernel 3 while neighbouring targets still read
      ! the source stencil, so the source is copied to the scratch `c_old`.)
      do concurrent(j=1:ny, i=1:nx) local(k, zacc_o, zacc_n)
         zacc_o = 0.0_wp
         zacc_n = 0.0_wp
         do k = 1, nz
            zacc_o = zacc_o + h_old(i, j, k)
            zacc_n = zacc_n + h_new(i, j, k)
            z_old(i, j, k) = zacc_o
            z_new(i, j, k) = zacc_n
            c_old(i, j, k) = c(i, j, k)
         end do
      end do

      ! ---- kernel 2: limited parabola of every source cell ----
      ! Every stencil read is DONOR-substituted (`remap_c_don_at`), so a
      ! filler's stale value never tilts a live cell's parabola.  A filler's
      ! own profile is the donor parabola shifted by `s = c_fill − c_donor`
      ! (q_l, q_r += s; q6 unchanged), which integrates to the filler's OWN
      ! content: the remap stays exactly conservative per source cell, and a
      ! filler that maps onto itself keeps its value.
      do concurrent(j=1:ny, k=1:nz, i=1:nx) &
         local(c_m2, c_m1, c_0, c_p1, c_p2, e_lo, e_hi, q_l, q_r, q_6, s_fill)
         c_m2 = remap_c_don_at(nx, ny, nz, i, j, k - 2, c_old, h_old)
         c_m1 = remap_c_don_at(nx, ny, nz, i, j, k - 1, c_old, h_old)
         c_0 = remap_c_don_at(nx, ny, nz, i, j, k, c_old, h_old)
         c_p1 = remap_c_don_at(nx, ny, nz, i, j, k + 1, c_old, h_old)
         c_p2 = remap_c_don_at(nx, ny, nz, i, j, k + 2, c_old, h_old)
         e_lo = 0.0_wp
         if (k > 1) e_lo = remap_ppm_edge(k - 1, nz, c_m2, c_m1, c_0, c_p1)
         call remap_ppm_cell(k, nz, c_m1, c_0, c_p1, c_p2, e_lo, e_hi, q_l, q_r, q_6)
         if (h_old(i, j, k) <= H_VANISHED) then
            s_fill = c_old(i, j, k) - c_0
            q_l = q_l + s_fill
            q_r = q_r + s_fill
         end if
         q_lft(i, j, k) = q_l
         q_rgt(i, j, k) = q_r
         q_six(i, j, k) = q_6
      end do

      ! ---- kernel 3: overlap integral + bounding + commit, one thread per target cell ----
      ! One loop for budget / no budget (splitting present(budget) makes NVHPC
      ! compile the no-budget branch ~25x slower).
      do concurrent(j=1:ny, k=1:nz, i=1:nx) &
         local(ko, lo, hi, mid, stp, k_fb, n_hit, q_l, q_r, q_6, zo_lo, zo_hi, dzo, zn_lo, zn_hi, dzn, &
               z_lo, z_hi, overlap, xi_lo, xi_hi, integral, c_raw, c_new, p_new, p_old, c_lo, c_hi, cs_lo, cs_hi, &
               c_m1, c_0, c_p1, s_fill, km, kp, m, m_don)
         dzn = h_new(i, j, k)
         zn_hi = z_new(i, j, k)
         zn_lo = 0.0_wp
         if (k > 1) zn_lo = z_new(i, j, k - 1)
         ! First source cell whose top lies above the target's bottom: the
         ! smallest m with z_old(m) > zn_lo (z_old is non-decreasing;
         ! m = nz+1 if none).  Gallop outward from the same index k --
         ! neighbouring grids are close, so the lanes of a warp read the
         ! same layer and stay coalesced -- then bisect.  The answer is
         ! independent of the starting guess.
         ko = min(k, nz)
         if (z_old(i, j, ko) > zn_lo) then
            hi = ko
            lo = 0
            stp = 1
            do while (hi > 1)
               mid = hi - stp
               if (mid < 1) exit
               if (z_old(i, j, mid) > zn_lo) then
                  hi = mid
                  stp = 2*stp
               else
                  lo = mid
                  exit
               end if
            end do
         else
            lo = ko
            hi = nz + 1
            stp = 1
            do while (lo < nz)
               mid = lo + stp
               if (mid > nz) exit
               if (z_old(i, j, mid) > zn_lo) then
                  hi = mid
                  exit
               end if
               lo = mid
               stp = 2*stp
            end do
         end if
         do while (hi - lo > 1)
            mid = (lo + hi)/2
            if (z_old(i, j, mid) > zn_lo) then
               hi = mid
            else
               lo = mid
            end if
         end do
         ko = hi
         k_fb = min(ko, nz)
         c_raw = 0.0_wp
         n_hit = 0
         c_lo = 0.0_wp
         c_hi = 0.0_wp
         if (dzn > 0.0_wp) then
            zo_lo = 0.0_wp
            if (ko > 1 .and. ko <= nz) zo_lo = z_old(i, j, ko - 1)
            integral = 0.0_wp
            do while (ko <= nz)
               if (zo_lo >= zn_hi) exit
               zo_hi = z_old(i, j, ko)
               z_lo = max(zn_lo, zo_lo)
               z_hi = min(zn_hi, zo_hi)
               overlap = z_hi - z_lo
               dzo = h_old(i, j, ko)
               if (overlap > 0.0_wp .and. dzo > 0.0_wp) then
                  q_l = q_lft(i, j, ko)
                  q_r = q_rgt(i, j, ko)
                  q_6 = q_six(i, j, ko)
                  xi_lo = (z_lo - zo_lo)/dzo
                  xi_hi = (z_hi - zo_lo)/dzo
                  integral = integral + dzo*( &
                             (xi_hi - xi_lo)*q_l &
                             + 0.5_wp*(xi_hi*xi_hi - xi_lo*xi_lo)*(q_r - q_l + q_6) &
                             - (xi_hi*xi_hi*xi_hi - xi_lo*xi_lo*xi_lo)*q_6/3.0_wp)
                  ! Range of the contributor's (donor-substituted) stencil;
                  ! a filler's profile is shifted by its `s_fill`, so is its
                  ! range.  The donor walk is INLINE: nvfortran miscompiles a
                  ! `remap_c_don_at` call in this kernel (kernel 2 is fine).
                  ! One walk suffices.  Live `ko`: a filler below reads `ko`;
                  ! a filler above reads the first live layer above it, else
                  ! `ko`.  Filler `ko`: its donor `c_0`, which its filler
                  ! neighbours share.
                  km = max(ko - 1, 1)
                  kp = min(ko + 1, nz)
                  if (dzo > H_VANISHED) then
                     c_0 = c_old(i, j, ko)
                     m_don = ko
                     if (h_old(i, j, kp) <= H_VANISHED) then
                        do m = nz, kp + 1, -1
                           if (h_old(i, j, m) > H_VANISHED) m_don = m
                        end do
                     end if
                     c_p1 = c_old(i, j, kp)
                     if (h_old(i, j, kp) <= H_VANISHED) c_p1 = c_old(i, j, m_don)
                     c_m1 = c_old(i, j, km)
                     if (h_old(i, j, km) <= H_VANISHED) c_m1 = c_0
                  else
                     m_don = 0
                     do m = nz, ko + 1, -1
                        if (h_old(i, j, m) > H_VANISHED) m_don = m
                     end do
                     if (m_don == 0) then
                        do m = 1, ko - 1
                           if (h_old(i, j, m) > H_VANISHED) m_don = m
                        end do
                     end if
                     c_0 = 0.0_wp
                     if (m_don > 0) c_0 = c_old(i, j, m_don)
                     c_m1 = c_0
                     if (h_old(i, j, km) > H_VANISHED) c_m1 = c_old(i, j, km)
                     c_p1 = c_0
                     if (h_old(i, j, kp) > H_VANISHED) c_p1 = c_old(i, j, kp)
                  end if
                  cs_lo = min(c_m1, min(c_0, c_p1))
                  cs_hi = max(c_m1, max(c_0, c_p1))
                  if (dzo <= H_VANISHED) then
                     s_fill = c_old(i, j, ko) - c_0
                     cs_lo = cs_lo + s_fill
                     cs_hi = cs_hi + s_fill
                  end if
                  if (n_hit == 0) then
                     c_lo = cs_lo
                     c_hi = cs_hi
                  else
                     c_lo = min(c_lo, cs_lo)
                     c_hi = max(c_hi, cs_hi)
                  end if
                  n_hit = n_hit + 1
               end if
               zo_lo = zo_hi
               ko = ko + 1
            end do
            c_raw = integral/dzn
         end if
         c_new = c_raw
         if (n_hit > 0) then
            if (ieee_is_finite(c_raw)) c_new = max(c_lo, min(c_hi, c_raw))
         else
            c_new = c_old(i, j, k_fb)
         end if
         if (present(budget)) then
            p_new = h_new(i, j, k)*c_new
            p_old = h_old(i, j, k)*c_old(i, j, k)
            budget(i, j, k) = budget(i, j, k) + (p_new - p_old)
         end if
         if (present(c_new_out)) c_new_out(i, j, k) = c_new
         if (present(clamp_out)) clamp_out(i, j, k) = dzn*(c_new - c_raw)
         c(i, j, k) = c_new
      end do
   end subroutine ocean_remap_tracer_field_stream

   subroutine ocean_apply_ale_remap_faces(grid, h_old, h_new, u_face_x, v_face_y, method, &
                                          conserve_ke, zlevel_faces, bnd_extrap, nonunif)
      !! Face-velocity pass of the ALE remap. Remaps u_face_x_layer and
      !! v_face_y_layer h_old→h_new using arithmetic-mean face thicknesses and
      !! the per-column kernel. Public only for the unit-test suite.
      !! Velocity treated as the face "concentration" (analogous to the tracer c at
      !! centres); per-face conservation sum_k(h_face·u_face) preserved
      !! (momentum-conserving). Outer-wall faces take the adjacent cell's
      !! thickness verbatim (no across-cell to average).
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: h_old(:, :, :)
      real(wp), intent(in) :: h_new(:, :, :)
      real(wp), intent(inout) :: u_face_x(:, :, :)
         !! Eastward face velocity, shape `(nx+1, ny, nz)`.
      real(wp), intent(inout) :: v_face_y(:, :, :)
         !! Northward face velocity, shape `(nx, ny+1, nz)`.
      integer, intent(in), optional :: method
      logical, intent(in), optional :: conserve_ke
         !! Enable the KE-conserving baroclinic-anomaly rescale (default
         !! `.false.` ⇒ momentum-only remap, bit-identical).
      logical, intent(in), optional :: zlevel_faces
         !! `&vcoord_nml zfixed_closed_faces` — build the face columns as
         !! `min(h_L, h_R)` and drop the closed layers, so the remap
         !! neither fills nor drains an inert filler.  Default `.false.`
         !! ⇒ bit-identical.  See `remap_x_face_velocity`.
      logical, intent(in), optional :: bnd_extrap
         !! Linear-exact boundary-cell reconstruction (default `.false.`
         !! ⇒ the PCM flatten, bit-identical).
      logical, intent(in), optional :: nonunif
         !! Non-uniform-grid PLM/PPM weights (default `.false.` ⇒ the
         !! equal-thickness specialisations, bit-identical).

      integer :: m, nx, ny, nz
      logical :: ke, zf, be, nu
      m = REMAP_PPM
      if (present(method)) m = method
      ke = .false.
      if (present(conserve_ke)) ke = conserve_ke
      zf = .false.
      if (present(zlevel_faces)) zf = zlevel_faces
      be = .false.
      if (present(bnd_extrap)) be = bnd_extrap
      nu = .false.
      if (present(nonunif)) nu = nonunif
      nx = grid%nx_total
      ny = grid%ny_total
      nz = size(u_face_x, 3)

      call remap_x_face_velocity(nx, ny, nz, h_old, h_new, u_face_x, m, ke, zf, be, nu)
      call remap_y_face_velocity(nx, ny, nz, h_old, h_new, v_face_y, m, ke, zf, be, nu)
   end subroutine ocean_apply_ale_remap_faces

   pure subroutine remap_x_face_velocity(nx, ny, nz, h_old, h_new, u_face_x, method, &
                                         conserve_ke, zlevel_faces, bnd_extrap, nonunif)
      !! Flat-impl x-face remap. East faces at i+1/2; u_face_x(1..nx+1) covers
      !! west wall (I=1), interior (I=2..nx), east wall (I=nx+1).
      !! `conserve_ke` (default .false.): rescale the baroclinic anomaly so column
      !! anomaly KE matches pre-remap, barotropic mean preserved. See
      !! `rescale_anomaly_ke`.
      integer, intent(in) :: nx, ny, nz, method
      real(wp), intent(in) :: h_old(nx, ny, nz)
      real(wp), intent(in) :: h_new(nx, ny, nz)
      real(wp), intent(inout) :: u_face_x(nx + 1, ny, nz)
      logical, intent(in) :: conserve_ke
      logical, intent(in) :: zlevel_faces
         !! `&vcoord_nml zfixed_closed_faces` — z-level partial steps.
         !! Builds the face column as the OVERLAP of the two cell
         !! columns, `min(h_L, h_R)`, instead of their arithmetic mean,
         !! and then drops any layer at or below `H_VANISHED` from BOTH
         !! the source and the target column.  A layer that is an inert
         !! filler on one side then has an exactly-zero target thickness,
         !! which is the one case every `remap_column` variant already
         !! short-circuits (`dz_new <= 0 ⇒ q_new = 0`), so no momentum is
         !! poured INTO a closed layer and none is carried OUT of one.
         !!
         !! With the arithmetic mean a one-sided filler gives
         !! `0.5·h_live > 0` — an ordinary half-thickness live layer as
         !! far as the column kernel is concerned — and the remap happily
         !! deposits momentum in a cell that geometrically holds no
         !! water.  `min` is also the standard partial-cell face
         !! thickness, so this is not a special case bolted on: it is the
         !! face thickness a z-level model has all along.
         !!
         !! `.false.` (default) ⇒ every expression is textually the
         !! arithmetic-mean form ⇒ bit-identical.
      logical, intent(in) :: bnd_extrap
         !! Linear-exact boundary-cell reconstruction in the column kernel.
      logical, intent(in) :: nonunif
         !! Non-uniform-grid PLM/PPM weights in the column kernel.
      integer :: I, j, k
      real(wp) :: h_old_face(NZ_STACK_MAX), h_new_face(NZ_STACK_MAX)
      real(wp) :: u_old_col(NZ_STACK_MAX), u_new_col(NZ_STACK_MAX)

      do concurrent(j=1:ny, I=1:nx + 1) &
         local(k, h_old_face, h_new_face, u_old_col, u_new_col)
         do k = 1, nz
            if (I == 1) then
               h_old_face(k) = h_old(1, j, k)
               h_new_face(k) = h_new(1, j, k)
            else if (I == nx + 1) then
               h_old_face(k) = h_old(nx, j, k)
               h_new_face(k) = h_new(nx, j, k)
            else
               h_old_face(k) = 0.5_wp*(h_old(I - 1, j, k) + h_old(I, j, k))
               h_new_face(k) = 0.5_wp*(h_new(I - 1, j, k) + h_new(I, j, k))
            end if
            if (zlevel_faces) then
               if (I == 1) then
                  h_old_face(k) = h_old(1, j, k)
                  h_new_face(k) = h_new(1, j, k)
               else if (I == nx + 1) then
                  h_old_face(k) = h_old(nx, j, k)
                  h_new_face(k) = h_new(nx, j, k)
               else
                  h_old_face(k) = min(h_old(I - 1, j, k), h_old(I, j, k))
                  h_new_face(k) = min(h_new(I - 1, j, k), h_new(I, j, k))
               end if
               if (h_old_face(k) <= H_VANISHED .or. h_new_face(k) <= H_VANISHED) then
                  h_old_face(k) = 0.0_wp
                  h_new_face(k) = 0.0_wp
               end if
            end if
            u_old_col(k) = u_face_x(I, j, k)
         end do
         call remap_column(method, nz, &
                           h_old_face(1:nz), h_new_face(1:nz), &
                           u_old_col(1:nz), u_new_col(1:nz), bnd_extrap, nonunif)
         if (conserve_ke) then
            call rescale_anomaly_ke(nz, h_old_face, h_new_face, u_old_col, u_new_col)
         end if
         do k = 1, nz
            u_face_x(I, j, k) = u_new_col(k)
         end do
      end do
   end subroutine remap_x_face_velocity

   pure subroutine remap_y_face_velocity(nx, ny, nz, h_old, h_new, v_face_y, method, &
                                         conserve_ke, zlevel_faces, bnd_extrap, nonunif)
      !! Flat-impl y-face remap, mirror of `remap_x_face_velocity`.  See
      !! that routine for the `conserve_ke` / `bnd_extrap` / `nonunif`
      !! semantics.
      integer, intent(in) :: nx, ny, nz, method
      real(wp), intent(in) :: h_old(nx, ny, nz)
      real(wp), intent(in) :: h_new(nx, ny, nz)
      real(wp), intent(inout) :: v_face_y(nx, ny + 1, nz)
      logical, intent(in) :: conserve_ke
      logical, intent(in) :: zlevel_faces
         !! `&vcoord_nml zfixed_closed_faces` — z-level partial steps.
         !! Builds the face column as the OVERLAP of the two cell
         !! columns, `min(h_L, h_R)`, instead of their arithmetic mean,
         !! and then drops any layer at or below `H_VANISHED` from BOTH
         !! the source and the target column.  A layer that is an inert
         !! filler on one side then has an exactly-zero target thickness,
         !! which is the one case every `remap_column` variant already
         !! short-circuits (`dz_new <= 0 ⇒ q_new = 0`), so no momentum is
         !! poured INTO a closed layer and none is carried OUT of one.
         !!
         !! With the arithmetic mean a one-sided filler gives
         !! `0.5·h_live > 0` — an ordinary half-thickness live layer as
         !! far as the column kernel is concerned — and the remap happily
         !! deposits momentum in a cell that geometrically holds no
         !! water.  `min` is also the standard partial-cell face
         !! thickness, so this is not a special case bolted on: it is the
         !! face thickness a z-level model has all along.
         !!
         !! `.false.` (default) ⇒ every expression is textually the
         !! arithmetic-mean form ⇒ bit-identical.
      logical, intent(in) :: bnd_extrap
      logical, intent(in) :: nonunif
      integer :: i, J, k
      real(wp) :: h_old_face(NZ_STACK_MAX), h_new_face(NZ_STACK_MAX)
      real(wp) :: v_old_col(NZ_STACK_MAX), v_new_col(NZ_STACK_MAX)

      do concurrent(J=1:ny + 1, i=1:nx) &
         local(k, h_old_face, h_new_face, v_old_col, v_new_col)
         do k = 1, nz
            if (J == 1) then
               h_old_face(k) = h_old(i, 1, k)
               h_new_face(k) = h_new(i, 1, k)
            else if (J == ny + 1) then
               h_old_face(k) = h_old(i, ny, k)
               h_new_face(k) = h_new(i, ny, k)
            else
               h_old_face(k) = 0.5_wp*(h_old(i, J - 1, k) + h_old(i, J, k))
               h_new_face(k) = 0.5_wp*(h_new(i, J - 1, k) + h_new(i, J, k))
            end if
            if (zlevel_faces) then
               if (J == 1) then
                  h_old_face(k) = h_old(i, 1, k)
                  h_new_face(k) = h_new(i, 1, k)
               else if (J == ny + 1) then
                  h_old_face(k) = h_old(i, ny, k)
                  h_new_face(k) = h_new(i, ny, k)
               else
                  h_old_face(k) = min(h_old(i, J - 1, k), h_old(i, J, k))
                  h_new_face(k) = min(h_new(i, J - 1, k), h_new(i, J, k))
               end if
               if (h_old_face(k) <= H_VANISHED .or. h_new_face(k) <= H_VANISHED) then
                  h_old_face(k) = 0.0_wp
                  h_new_face(k) = 0.0_wp
               end if
            end if
            v_old_col(k) = v_face_y(i, J, k)
         end do
         call remap_column(method, nz, &
                           h_old_face(1:nz), h_new_face(1:nz), &
                           v_old_col(1:nz), v_new_col(1:nz), bnd_extrap, nonunif)
         if (conserve_ke) then
            call rescale_anomaly_ke(nz, h_old_face, h_new_face, v_old_col, v_new_col)
         end if
         do k = 1, nz
            v_face_y(i, J, k) = v_new_col(k)
         end do
      end do
   end subroutine remap_y_face_velocity

   pure subroutine rescale_anomaly_ke(nz, h_old_face, h_new_face, u_old_col, u_new_col)
      !$acc routine seq
      !$omp declare target
      !! KE-conserving rescale of a remapped face-velocity column.
      !! The column remap conserves momentum (Σ h·u) but not KE; restore it by
      !! rescaling ONLY the baroclinic anomaly (Adcroft & Hallberg 2006):
      !!   scale = sqrt(KE_old_anom/KE_new_anom), clamped to [0, 1.25]
      !!   u_new(k) = u_bar_new + scale·(u_new(k) - u_bar_new)
      !! Barotropic mean u_bar preserved verbatim; degenerate columns untouched.
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_old_face(NZ_STACK_MAX)
      real(wp), intent(in) :: h_new_face(NZ_STACK_MAX)
      real(wp), intent(in) :: u_old_col(NZ_STACK_MAX)
      real(wp), intent(inout) :: u_new_col(NZ_STACK_MAX)
      integer :: k
      real(wp) :: h_old_sum, h_new_sum, mom_old, mom_new
      real(wp) :: u_bar_old, u_bar_new, ke_old, ke_new, anom, scale_fac

      h_old_sum = 0.0_wp
      h_new_sum = 0.0_wp
      mom_old = 0.0_wp
      mom_new = 0.0_wp
      do k = 1, nz
         h_old_sum = h_old_sum + h_old_face(k)
         h_new_sum = h_new_sum + h_new_face(k)
         mom_old = mom_old + h_old_face(k)*u_old_col(k)
         mom_new = mom_new + h_new_face(k)*u_new_col(k)
      end do
      ! Dry / vanishing column: nothing to rescale.
      if (h_old_sum <= H_FLOOR .or. h_new_sum <= H_FLOOR) return
      u_bar_old = mom_old/h_old_sum
      u_bar_new = mom_new/h_new_sum

      ke_old = 0.0_wp
      ke_new = 0.0_wp
      do k = 1, nz
         anom = u_old_col(k) - u_bar_old
         ke_old = ke_old + 0.5_wp*h_old_face(k)*anom*anom
         anom = u_new_col(k) - u_bar_new
         ke_new = ke_new + 0.5_wp*h_new_face(k)*anom*anom
      end do
      ! No anomaly KE on either side ⇒ pure barotropic column, leave as-is.
      if (ke_new <= 0.0_wp .or. ke_old <= 0.0_wp) return
      scale_fac = sqrt(ke_old/ke_new)
      if (scale_fac > 1.25_wp) scale_fac = 1.25_wp
      if (scale_fac < 0.0_wp) scale_fac = 0.0_wp
      do k = 1, nz
         u_new_col(k) = u_bar_new + scale_fac*(u_new_col(k) - u_bar_new)
      end do
   end subroutine rescale_anomaly_ke

   pure subroutine ocean_apply_ale_remap_step(grid, vcoord, ms, bt_eta, bt_H_ref, method, eos, dt)
      !! Top-level entry the driver calls between outer steps: snapshot h_old once,
      !! remap centres (h_layer + tracers) AND faces using the same snapshot, then
      !! re-derive bt_eta. Returns early for EULERIAN_Z/LAGRANGIAN.
      !! `eos` (optional): required ONLY for VCOORD_RHO (isopycnal density inversion).
      !! `dt` (optional, s): only for the grid time-filter (regrid_time_scale > 0);
      !! absent or τ=0 (default) ⇒ filter skipped, bit-identical.
      type(hgrid_t), intent(in) :: grid
      type(ocean_vcoord_t), intent(inout) :: vcoord
      type(multilayer_state_t), intent(inout) :: ms
      ! assumed-shape-ok: outer-driver allocatable; grid%nx_total available; deferred.
      real(wp), intent(inout) :: bt_eta(:, :)
      real(wp), intent(in) :: bt_H_ref(:, :)  ! assumed-shape-ok: outer-driver allocatable; grid%nx_total available; deferred
      integer, intent(in), optional :: method
      type(eos_t), intent(in), optional :: eos
      real(wp), intent(in), optional :: dt

      integer :: t, m, nx, ny, nz, i, j, k
      real(wp) :: wtd
      logical :: do_tfilter, conserve_ke, bnd_extrap, nonunif, use_stream

      ! Eulerian-z and Lagrangian/isopycnal both skip remap: the former
      ! holds h at H·dsig via vert-advection cancellation, the latter
      ! lets h evolve freely (target = current h).
      if (vcoord%coord_type == VCOORD_EULERIAN_Z .or. &
          vcoord%coord_type == VCOORD_LAGRANGIAN) return
      if (.not. vcoord%is_init) return
      if (.not. allocated(ms%h_layer)) return

      m = REMAP_PPM
      if (present(method)) m = method
      bnd_extrap = vcoord%remap_boundary_extrap
      nonunif = vcoord%remap_nonuniform_weights
      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml

      ! 1. Column total + target_h (persistent scratch on vcoord — per-call
      ! allocates here generate H↔D transfers, not being in the present table).
      do concurrent(j=1:ny, i=1:nx) local(k)
         vcoord%remap_total_h(i, j) = 0.0_wp
         do k = 1, nz
            vcoord%remap_total_h(i, j) = vcoord%remap_total_h(i, j) &
                                         + ms%h_layer(i, j, k)
         end do
      end do
      do concurrent(j=1:ny, i=1:nx)
         vcoord%remap_h_ref(i, j) = vcoord%remap_total_h(i, j) - bt_eta(i, j)
      end do

      ! 2. Snapshot h_old on-device (host source= reads stale memory after a
      ! dynamics-side OpenACC kernel). Before target_h for VCOORD_RHO's inversion.
      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         vcoord%remap_h_old(i, j, k) = ms%h_layer(i, j, k)
      end do

      ! 3. Populate target_h. Geometric coords use compute_target_h; isopycnal
      ! needs per-layer T/S concentrations + EOS via compute_target_h_rho.
      ! Concentrations copied into persistent scratch.
      if ((vcoord%coord_type == VCOORD_RHO .or. vcoord%coord_type == VCOORD_HYCOM) &
          .and. present(eos) .and. allocated(ms%tracers) &
          .and. ms%idx_temperature > 0 .and. ms%idx_salinity > 0) then
         call build_ts_concentration(nx, ny, nz, vcoord%remap_h_old, &
                                     ms%tracers(ms%idx_temperature)%c, &
                                     ms%tracers(ms%idx_salinity)%c, &
                                     vcoord%remap_conc_t, vcoord%remap_conc_s)
         ! HYCOM = RHO inversion + z*-floor/monotonize deltas (hybrid=.true.);
         ! pure RHO passes hybrid=.false. (bit-identical to RHO-only kernel).
         call vcoord%compute_target_h_rho(vcoord%remap_h_ref, bt_eta, &
                                          vcoord%remap_conc_t, vcoord%remap_conc_s, eos, &
                                          hybrid=(vcoord%coord_type == VCOORD_HYCOM))
      else
         call vcoord%compute_target_h(vcoord%remap_h_ref, bt_eta)
      end if

      ! 3b. Grid time-filter (White & Adcroft 2008): relax target toward it from
      ! the old grid by wtd = dt/(τ+dt). Convex blend ⇒ column total conserved.
      ! τ=0 (default) or dt absent ⇒ skipped, bit-identical.
      do_tfilter = vcoord%regrid_time_scale > 0.0_wp .and. present(dt)
      if (do_tfilter) then
         wtd = dt/(vcoord%regrid_time_scale + dt)
         do concurrent(k=1:nz, j=1:ny, i=1:nx)
            vcoord%target_h(i, j, k) = vcoord%remap_h_old(i, j, k) &
                                       + wtd*(vcoord%target_h(i, j, k) - vcoord%remap_h_old(i, j, k))
         end do
      end if

      ! `vcoord%remap_h_old` and `vcoord%target_h` are the exact pair every
      ! column kernel below consumes, and neither is overwritten by steps 4-7
      ! — so `ocean_remap_scan_preconditions` can assert them from the
      ! (impure) driver AFTER this call.  See `&vcoord_nml
      ! remap_check_preconditions` and `rdb_ocean_dyn :: ocean_dyn_step_split`.

      ! 4. Tracer remap (centre cells)
      use_stream = remap_use_stream(vcoord, nz, m, bnd_extrap, nonunif)
      if (allocated(ms%tracers)) then
         do t = 1, size(ms%tracers)
            if (.not. allocated(ms%tracers(t)%c)) cycle
            select case (ms%tracers(t)%budget_id)
            case (TRACER_BUDGET_HEAT)
               if (use_stream) then
                  call ocean_remap_tracer_field_stream( &
                     nx, ny, nz, vcoord%remap_h_old, vcoord%target_h, ms%tracers(t)%c, &
                     vcoord%remap_c_old, vcoord%remap_z_old, vcoord%remap_z_new, &
                     vcoord%remap_q_l, vcoord%remap_q_r, vcoord%remap_q_6, &
                     budget=ms%heat_budget_remap)
               else
                  call ocean_remap_tracer_field( &
                     nx, ny, nz, vcoord%remap_h_old, vcoord%target_h, ms%tracers(t)%c, m, &
                     bnd_extrap, nonunif, budget=ms%heat_budget_remap)
               end if
            case (TRACER_BUDGET_SALT)
               if (use_stream) then
                  call ocean_remap_tracer_field_stream( &
                     nx, ny, nz, vcoord%remap_h_old, vcoord%target_h, ms%tracers(t)%c, &
                     vcoord%remap_c_old, vcoord%remap_z_old, vcoord%remap_z_new, &
                     vcoord%remap_q_l, vcoord%remap_q_r, vcoord%remap_q_6, &
                     budget=ms%salt_budget_remap)
               else
                  call ocean_remap_tracer_field( &
                     nx, ny, nz, vcoord%remap_h_old, vcoord%target_h, ms%tracers(t)%c, m, &
                     bnd_extrap, nonunif, budget=ms%salt_budget_remap)
               end if
            case default
               if (use_stream) then
                  call ocean_remap_tracer_field_stream( &
                     nx, ny, nz, vcoord%remap_h_old, vcoord%target_h, ms%tracers(t)%c, &
                     vcoord%remap_c_old, vcoord%remap_z_old, vcoord%remap_z_new, &
                     vcoord%remap_q_l, vcoord%remap_q_r, vcoord%remap_q_6)
               else
                  call ocean_remap_tracer_field( &
                     nx, ny, nz, vcoord%remap_h_old, vcoord%target_h, ms%tracers(t)%c, m, &
                     bnd_extrap, nonunif)
               end if
            end select
         end do
      end if

      ! 5. Face-velocity remap (h_old → target_h). KE-conserving anomaly rescale
      ! gated on the vcoord knob (default off ⇒ momentum-only, bit-identical).
      conserve_ke = vcoord%remap_vel_conserve_ke
      if (allocated(ms%u_face_x_layer) .and. allocated(ms%v_face_y_layer)) then
         call remap_x_face_velocity(nx, ny, nz, vcoord%remap_h_old, vcoord%target_h, &
                                    ms%u_face_x_layer, m, conserve_ke, &
                                    vcoord%zfixed_closed_faces, bnd_extrap, nonunif)
         call remap_y_face_velocity(nx, ny, nz, vcoord%remap_h_old, vcoord%target_h, &
                                    ms%v_face_y_layer, m, conserve_ke, &
                                    vcoord%zfixed_closed_faces, bnd_extrap, nonunif)
      end if

      ! 6. h_layer = target_h; capture mass-budget delta
      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         ms%mass_budget_remap(i, j, k) = ms%mass_budget_remap(i, j, k) &
                                         + (vcoord%target_h(i, j, k) - vcoord%remap_h_old(i, j, k))
         ms%h_layer(i, j, k) = vcoord%target_h(i, j, k)
      end do

      ! 7. Re-derive bt_eta
      do concurrent(j=1:ny, i=1:nx) local(k)
         bt_eta(i, j) = -bt_H_ref(i, j)
         do k = 1, nz
            bt_eta(i, j) = bt_eta(i, j) + ms%h_layer(i, j, k)
         end do
      end do
   end subroutine ocean_apply_ale_remap_step

   pure subroutine ocean_remap_scan_preconditions(nx, ny, nz, h_old, h_new, rel_tol, &
                                                  n_bad, worst_rel, worst_neg)
      !! Scan every column for the two remap preconditions and report how
      !! badly they are missed — the domain-wide counterpart of
      !! `rdb_remap_column :: remap_column_preconditions_ok`.
      !!
      !! `pure` and flat-arg so the caller owns the fail-loud decision and
      !! the reduction can run on-device; three scalars come back rather than
      !! a per-column field, so a thermo-cadence call costs one pass plus a
      !! tiny D→H.  Public for the unit-test suite.
      !!
      !! CALLER CONTRACT under `-gpu=...,mem:separate`: `h_old` and `h_new`
      !! must be device-present — the reduction declares them `present` rather
      !! than letting them copy in, so a caller that forgot the map ABORTS
      !! instead of silently scanning stale host values and reporting a clean
      !! bill of health.  Production satisfies it through
      !! `ocean_vcoord_enter_data_impl`, which maps `remap_h_old` and
      !! `target_h`; a host-side caller must add its own
      !! `!$acc enter data copyin(...)` plus an `update device` after every
      !! host edit.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: h_old(nx, ny, nz)
         !! Source thicknesses (the pre-remap `h_layer` snapshot).
      real(wp), intent(in) :: h_new(nx, ny, nz)
         !! Target thicknesses (`vcoord%target_h`, after any time filter).
      real(wp), intent(in) :: rel_tol
         !! Relative tolerance on the column-total match.
      integer, intent(out) :: n_bad
         !! Number of columns violating either precondition.
      real(wp), intent(out) :: worst_rel
         !! Largest relative column-total mismatch over the domain.
      real(wp), intent(out) :: worst_neg
         !! Most negative thickness found (0 when there is none).

      integer :: i, j, k
      real(wp) :: s_old, s_new, rel, hmin

      n_bad = 0
      worst_rel = 0.0_wp
      worst_neg = 0.0_wp
      do concurrent(j=1:ny, i=1:nx) local(k, s_old, s_new, rel, hmin) &
         reduce(+:n_bad) reduce(max:worst_rel) reduce(min:worst_neg)
         s_old = 0.0_wp
         s_new = 0.0_wp
         hmin = 0.0_wp
         do k = 1, nz
            s_old = s_old + h_old(i, j, k)
            s_new = s_new + h_new(i, j, k)
            hmin = min(hmin, h_old(i, j, k), h_new(i, j, k))
         end do
         ! Relative to the LARGER total, so a land column (both zero)
         ! scores 0 rather than tripping on a 0/0.
         rel = abs(s_new - s_old)/max(abs(s_old), abs(s_new), H_DIV_EPS)
         worst_rel = max(worst_rel, rel)
         worst_neg = min(worst_neg, hmin)
         if (hmin < 0.0_wp .or. rel > rel_tol) n_bad = n_bad + 1
      end do
   end subroutine ocean_remap_scan_preconditions

   pure subroutine build_ts_concentration(nx, ny, nz, h_old, hTr_T, hTr_S, conc_t, conc_s)
      !! Copy the stored layer-mean T/S concentrations into the VCOORD_RHO /
      !! HYCOM density-inversion scratch (the tracers ARE concentrations, so
      !! this is a plain copy; `h_old` is kept for the call signature).
      !! Flat-impl, explicit-shape; one cadence-bounded launch per remap.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: h_old(nx, ny, nz)
      real(wp), intent(in) :: hTr_T(nx, ny, nz)
      real(wp), intent(in) :: hTr_S(nx, ny, nz)
      real(wp), intent(out) :: conc_t(nx, ny, nz)
      real(wp), intent(out) :: conc_s(nx, ny, nz)
      integer :: i, j, k
      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         conc_t(i, j, k) = hTr_T(i, j, k)
         conc_s(i, j, k) = hTr_S(i, j, k)
      end do
   end subroutine build_ts_concentration

   subroutine ocean_remap_tracer_column(nz, h_old, h_new, c_inout, method, c_new_out, clamp_out)
      !! Single-column unit-test entry: wraps `remap_column` with the SAME
      !! pipeline the production reference kernel runs (donor column ->
      !! `remap_column` -> `remap_add_filler_excess` -> `remap_bound_column`
      !! -> `c = c_new`), so this entry cannot drift from
      !! `ocean_remap_tracer_field`. Production callers go through that one.
      integer, intent(in) :: nz, method
      real(wp), intent(in) :: h_old(nz), h_new(nz)
      real(wp), intent(inout) :: c_inout(nz)
         !! Column concentration: read as `c_old`, overwritten with `c_new`.
      real(wp), intent(out), optional :: c_new_out(nz)
         !! The bounded target concentration.
      real(wp), intent(out), optional :: clamp_out(nz)
         !! The content the bounding moved, `h_new·(c_bounded − c_raw)`.
      real(wp) :: c_old(NZ_STACK_MAX), c_new(NZ_STACK_MAX), clamp_col(NZ_STACK_MAX)
      real(wp) :: c_don(NZ_STACK_MAX)
      real(wp) :: h_old_col(NZ_STACK_MAX), h_new_col(NZ_STACK_MAX)
      integer :: k
      do k = 1, nz
         h_old_col(k) = h_old(k)
         h_new_col(k) = h_new(k)
         c_old(k) = c_inout(k)
      end do
      call rdb_vl_column_conc(nz, h_old_col, c_old, c_don)
      call remap_column(method, nz, h_old, h_new, c_don(1:nz), c_new(1:nz))
      call remap_add_filler_excess(nz, h_old_col, h_new_col, c_don, c_old, c_new)
      call remap_bound_column(nz, h_old_col, h_new_col, c_don, c_old, c_new, clamp_col, .false.)
      do k = 1, nz
         c_inout(k) = c_new(k)
      end do
      if (present(c_new_out)) c_new_out(1:nz) = c_new(1:nz)
      if (present(clamp_out)) clamp_out(1:nz) = clamp_col(1:nz)
   end subroutine ocean_remap_tracer_column

#include "rdb_vanished_layer.inc"

end module rdb_ocean_remap
