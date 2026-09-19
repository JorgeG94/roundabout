!! ALE remap orchestrator for the ocean dynamical core (centre-cell pass).
module rdb_ocean_remap
   !! Drives the conservative vertical remap of `multilayer.h_layer` + every
   !! registered `multilayer.tracers(t)%hTr` from the current (Lagrangian) grid
   !! to `vcoord%target_h`. The kernel is the per-column `remap_column` from
   !! `rdb_remap_column` (shared with coastal); this module wires it across the
   !! registered tracer slot list plus the face-velocity pass.
   !!
   !! Per-column conservation: sum_k(c_old·h_old) = sum_k(c_new·h_new) to machine
   !! precision (modulo the c = hTr/h step, which a vanishing-layer guard protects).
#ifdef LFORTRAN_PASSING
   use rdb_constants, only: wp, REMAP_PPM, H_VANISHED
#else
   use rdb_constants, only: NZ_STACK_MAX, wp, REMAP_PPM, H_VANISHED
#endif
   use rdb_grid, only: hgrid_t
   use rdb_remap_column, only: remap_column
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_tracer, only: TRACER_BUDGET_HEAT, TRACER_BUDGET_SALT
   use rdb_ocean_vcoord, only: ocean_vcoord_t, VCOORD_EULERIAN_Z, VCOORD_LAGRANGIAN, &
                               VCOORD_RHO, VCOORD_HYCOM
   use rdb_eos, only: eos_t
   implicit none
   private

   ! Vanishing-layer guard for the `c = hTr / h` step — the D4 skip/merge
   ! marker, NOT a positivity floor: below it the layer's concentration is
   ! taken as 0 rather than recovered from a near-zero divisor.  Aliased to
   ! `H_VANISHED` (same value) so there is ONE definition of "vanished" in
   ! the tree; it used to be a bare `1.5e-4_wp` literal here, which is a
   ! third definition waiting to drift from the constant of record.  Every
   ! test of it is a STRICT `>`: a layer sitting exactly ON the marker reads
   ! as vanished, which is what the geometric vcoord families rely on (see
   ! `rdb_vcoord :: vcoord_h_min_role`).
   real(wp), parameter :: H_FLOOR = H_VANISHED

   public :: ocean_apply_ale_remap_centres
   public :: ocean_apply_ale_remap_faces
   public :: ocean_apply_ale_remap_step
   public :: ocean_remap_tracer_column   ! exposed for unit tests

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
      logical :: do_tfilter

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
                                     ms%tracers(ms%idx_temperature)%hTr, &
                                     ms%tracers(ms%idx_salinity)%hTr, &
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
      if (allocated(ms%tracers)) then
         do t = 1, size(ms%tracers)
            if (.not. allocated(ms%tracers(t)%hTr)) cycle
            select case (ms%tracers(t)%budget_id)
            case (TRACER_BUDGET_HEAT)
               call ocean_remap_tracer_field( &
                  nx, ny, nz, vcoord%remap_h_old, vcoord%target_h, ms%tracers(t)%hTr, m, &
                  budget=ms%heat_budget_remap)
            case (TRACER_BUDGET_SALT)
               call ocean_remap_tracer_field( &
                  nx, ny, nz, vcoord%remap_h_old, vcoord%target_h, ms%tracers(t)%hTr, m, &
                  budget=ms%salt_budget_remap)
            case default
               call ocean_remap_tracer_field( &
                  nx, ny, nz, vcoord%remap_h_old, vcoord%target_h, ms%tracers(t)%hTr, m)
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

   pure subroutine ocean_remap_tracer_field(nx, ny, nz, h_old, h_new, hTr, method, budget)
      !! Flat-impl tracer remap. Per (i,j) column: c = hTr/h (vanishing-layer-
      !! guarded) → per-column remap kernel → hTr_new = c_new·h_new. Conservative.
      !! `budget` (optional): when present, the per-cell hTr_new−hTr_old increment
      !! is accumulated into the slot (heat/salt remap deltas) before overwriting.
      !! Flat-arg so GPU codegen doesn't chase the array-of-derived-types pointer.
      integer, intent(in) :: nx, ny, nz, method
      real(wp), intent(in) :: h_old(nx, ny, nz)
      real(wp), intent(in) :: h_new(nx, ny, nz)
      real(wp), intent(inout) :: hTr(nx, ny, nz)
      real(wp), intent(inout), optional :: budget(nx, ny, nz)
      integer :: i, j, k
      real(wp) :: h_old_col(NZ_STACK_MAX), h_new_col(NZ_STACK_MAX)
      real(wp) :: c_old_col(NZ_STACK_MAX), c_new_col(NZ_STACK_MAX)
      real(wp) :: hTr_col(NZ_STACK_MAX)
      real(wp) :: hTr_new

      ! Gate the budget write INSIDE the one loop (as vdiff does); splitting
      ! present(budget) into two loops makes NVHPC compile the no-budget branch
      ! ~25x slower. Bit-identical to the split form.
      do concurrent(j=1:ny, i=1:nx) &
         local(k, h_old_col, h_new_col, c_old_col, c_new_col, hTr_col, hTr_new)
         do k = 1, nz
            h_old_col(k) = h_old(i, j, k)
            h_new_col(k) = h_new(i, j, k)
            hTr_col(k) = hTr(i, j, k)
            if (h_old_col(k) > H_FLOOR) then
               c_old_col(k) = hTr_col(k)/h_old_col(k)
            else
               c_old_col(k) = 0.0_wp
            end if
         end do
         call remap_column(method, nz, &
                           h_old_col(1:nz), h_new_col(1:nz), &
                           c_old_col(1:nz), c_new_col(1:nz))
         if (present(budget)) then
            do k = 1, nz
               hTr_new = c_new_col(k)*h_new_col(k)
               budget(i, j, k) = budget(i, j, k) + (hTr_new - hTr(i, j, k))
               hTr(i, j, k) = hTr_new
            end do
         else
            do k = 1, nz
               hTr(i, j, k) = c_new_col(k)*h_new_col(k)
            end do
         end if
      end do
   end subroutine ocean_remap_tracer_field

   subroutine ocean_apply_ale_remap_faces(grid, h_old, h_new, u_face_x, v_face_y, method, conserve_ke)
      !! Face-velocity pass of the ALE remap. Remaps u_face_x_layer and
      !! v_face_y_layer h_old→h_new using arithmetic-mean face thicknesses and
      !! the per-column kernel. Public only for the unit-test suite.
      !! Velocity treated as the face "concentration" (analogous to T=hTr/h at
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

      integer :: m, nx, ny, nz
      logical :: ke
      m = REMAP_PPM
      if (present(method)) m = method
      ke = .false.
      if (present(conserve_ke)) ke = conserve_ke
      nx = grid%nx_total
      ny = grid%ny_total
      nz = size(u_face_x, 3)

      call remap_x_face_velocity(nx, ny, nz, h_old, h_new, u_face_x, m, ke)
      call remap_y_face_velocity(nx, ny, nz, h_old, h_new, v_face_y, m, ke)
   end subroutine ocean_apply_ale_remap_faces

   pure subroutine remap_x_face_velocity(nx, ny, nz, h_old, h_new, u_face_x, method, conserve_ke)
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
            u_old_col(k) = u_face_x(I, j, k)
         end do
         call remap_column(method, nz, &
                           h_old_face(1:nz), h_new_face(1:nz), &
                           u_old_col(1:nz), u_new_col(1:nz))
         if (conserve_ke) then
            call rescale_anomaly_ke(nz, h_old_face, h_new_face, u_old_col, u_new_col)
         end if
         do k = 1, nz
            u_face_x(I, j, k) = u_new_col(k)
         end do
      end do
   end subroutine remap_x_face_velocity

   pure subroutine remap_y_face_velocity(nx, ny, nz, h_old, h_new, v_face_y, method, conserve_ke)
      !! Flat-impl y-face remap, mirror of `remap_x_face_velocity`.  See
      !! that routine for the `conserve_ke` semantics.
      integer, intent(in) :: nx, ny, nz, method
      real(wp), intent(in) :: h_old(nx, ny, nz)
      real(wp), intent(in) :: h_new(nx, ny, nz)
      real(wp), intent(inout) :: v_face_y(nx, ny + 1, nz)
      logical, intent(in) :: conserve_ke
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
            v_old_col(k) = v_face_y(i, J, k)
         end do
         call remap_column(method, nz, &
                           h_old_face(1:nz), h_new_face(1:nz), &
                           v_old_col(1:nz), v_new_col(1:nz))
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
      logical :: do_tfilter, conserve_ke

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
      ! Concentrations built into persistent scratch (guarded c = hTr/h).
      if ((vcoord%coord_type == VCOORD_RHO .or. vcoord%coord_type == VCOORD_HYCOM) &
          .and. present(eos) .and. allocated(ms%tracers) &
          .and. ms%idx_temperature > 0 .and. ms%idx_salinity > 0) then
         call build_ts_concentration(nx, ny, nz, vcoord%remap_h_old, &
                                     ms%tracers(ms%idx_temperature)%hTr, &
                                     ms%tracers(ms%idx_salinity)%hTr, &
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

      ! 4. Tracer remap (centre cells)
      if (allocated(ms%tracers)) then
         do t = 1, size(ms%tracers)
            if (.not. allocated(ms%tracers(t)%hTr)) cycle
            select case (ms%tracers(t)%budget_id)
            case (TRACER_BUDGET_HEAT)
               call ocean_remap_tracer_field( &
                  nx, ny, nz, vcoord%remap_h_old, vcoord%target_h, ms%tracers(t)%hTr, m, &
                  budget=ms%heat_budget_remap)
            case (TRACER_BUDGET_SALT)
               call ocean_remap_tracer_field( &
                  nx, ny, nz, vcoord%remap_h_old, vcoord%target_h, ms%tracers(t)%hTr, m, &
                  budget=ms%salt_budget_remap)
            case default
               call ocean_remap_tracer_field( &
                  nx, ny, nz, vcoord%remap_h_old, vcoord%target_h, ms%tracers(t)%hTr, m)
            end select
         end do
      end if

      ! 5. Face-velocity remap (h_old → target_h). KE-conserving anomaly rescale
      ! gated on the vcoord knob (default off ⇒ momentum-only, bit-identical).
      conserve_ke = vcoord%remap_vel_conserve_ke
      if (allocated(ms%u_face_x_layer) .and. allocated(ms%v_face_y_layer)) then
         call remap_x_face_velocity(nx, ny, nz, vcoord%remap_h_old, vcoord%target_h, &
                                    ms%u_face_x_layer, m, conserve_ke)
         call remap_y_face_velocity(nx, ny, nz, vcoord%remap_h_old, vcoord%target_h, &
                                    ms%v_face_y_layer, m, conserve_ke)
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

   pure subroutine build_ts_concentration(nx, ny, nz, h_old, hTr_T, hTr_S, conc_t, conc_s)
      !! Build layer-mean T/S concentrations (c = hTr/h) from extensive tracer
      !! content + pre-remap thicknesses, for the VCOORD_RHO density inversion.
      !! Guarded against H_VANISHED (sub-floor layer ⇒ c = 0). Flat-impl,
      !! explicit-shape; one cadence-bounded launch per remap.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: h_old(nx, ny, nz)
      real(wp), intent(in) :: hTr_T(nx, ny, nz)
      real(wp), intent(in) :: hTr_S(nx, ny, nz)
      real(wp), intent(out) :: conc_t(nx, ny, nz)
      real(wp), intent(out) :: conc_s(nx, ny, nz)
      integer :: i, j, k
      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         if (h_old(i, j, k) > H_VANISHED) then
            conc_t(i, j, k) = hTr_T(i, j, k)/h_old(i, j, k)
            conc_s(i, j, k) = hTr_S(i, j, k)/h_old(i, j, k)
         else
            conc_t(i, j, k) = 0.0_wp
            conc_s(i, j, k) = 0.0_wp
         end if
      end do
   end subroutine build_ts_concentration

   subroutine ocean_remap_tracer_column(nz, h_old, h_new, hTr_inout, method)
      !! Single-column unit-test entry: wraps `remap_column` with the
      !! c = hTr/h ↔ hTr_new = c_new·h_new pattern. Production callers go through
      !! `ocean_remap_tracer_field`.
      integer, intent(in) :: nz, method
      real(wp), intent(in) :: h_old(nz), h_new(nz)
      real(wp), intent(inout) :: hTr_inout(nz)
      real(wp) :: c_old(NZ_STACK_MAX), c_new(NZ_STACK_MAX)
      integer :: k
      do k = 1, nz
         if (h_old(k) > H_FLOOR) then
            c_old(k) = hTr_inout(k)/h_old(k)
         else
            c_old(k) = 0.0_wp
         end if
      end do
      call remap_column(method, nz, h_old, h_new, c_old(1:nz), c_new(1:nz))
      do k = 1, nz
         hTr_inout(k) = c_new(k)*h_new(k)
      end do
   end subroutine ocean_remap_tracer_column

end module rdb_ocean_remap
