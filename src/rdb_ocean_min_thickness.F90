module rdb_ocean_min_thickness
   !! Conservative minimum-layer-thickness adjustment for the isopycnal
   !! (`VCOORD_LAGRANGIAN`) ocean path.
   !!
   !! Motivation: in the remap-free Lagrangian coordinate an eddy can displace
   !! an interface so far that a layer thins to ~0 ("outcropping"), and then
   !! `u = hu/h` blows up.  The existing `angstrom_h` floor lifts a sub-floor
   !! layer with `max(h_new, angstrom_h)`, which INJECTS mass (non-conservative:
   !! total column thickness, SSH, and tracer mass all drift up on every floored
   !! step) and the injected volume shapes the eddy field.
   !!
   !! This module instead performs a CONSERVATIVE adjustment: when a layer
   !! thins below the floor, the deficit is borrowed from the surplus layers of
   !! the SAME water column — the interface is moved, no mass is created.  Per
   !! water column the following are conserved to round-off:
   !!   * total thickness   Sigma_k h_k   (SSH unchanged)
   !!   * momentum          Sigma_k h_face_k . u_face_k   (per C-grid face)
   !!   * every tracer mass Sigma_k h_k . Tr_k
   !!
   !! The adjustment is a STRICT no-op on any column where all layers already
   !! meet the floor: those columns (and their faces) are left byte-unchanged so
   !! the isopycnal interface structure — and hence the baroclinic mode we are
   !! modelling — is never pinned.  This is the load-bearing constraint that
   !! distinguishes this from a z*/sigma remap.
   !!
   !! Implementation: reuse the conservative per-column `remap_column` primitive
   !! (the same engine the ALE remap uses) with a "floor-only" target column —
   !! each sub-floor layer inflated to the floor, the excess drawn conservatively
   !! from the surplus layers so the column total is preserved.
   !!
   !! COST STRUCTURE (the 31%-of-runtime restructure, 2026-07-24): the borrow
   !! is a no-op on every non-grounded column, but the original implementation
   !! still paid ~8 full-field passes per call on a healthy domain
   !! (unconditional target build + copy-back, per-tracer concentration
   !! divisions before the gate, and per-face NEIGHBOUR COLUMN re-scans with
   !! k-strided access).  Now: ONE coalesced pass builds a 2D grounded mask +
   !! a global count; zero grounded columns ⇒ the whole call returns after
   !! that single read; otherwise every kernel consults the mask (2 loads) and
   !! only grounded columns / active faces do work — and reads fall back to
   !! `h_old` on non-grounded neighbours, where `h_new ≡ h_old` by
   !! construction.  Byte-identical results in all cases.
   use rdb_constants, only: wp, NZ_STACK_MAX, REMAP_PPM
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_remap_column, only: remap_column

   implicit none
   private

   public :: ocean_apply_conservative_min_thickness
   public :: min_thickness_target_column   ! exposed for unit tests

   real(wp), parameter :: H_CONC_FLOOR = 1.0e-20_wp
      !! Pure 1/0 armour when forming a layer concentration c = q/h_old; only
      !! guards genuinely-zero source layers (never a physical thickness).

contains

   pure subroutine min_thickness_target_column(nz, h_old, h_floor, h_new, grounded)
      !$acc routine seq
      !$omp declare target
      !! Build the floor-only conservative target thickness column.
      !!
      !! `grounded` (out): `.true.` iff any layer is strictly below `floor`.
      !! When `.false.`, `h_new` is an exact copy of `h_old` (no-op) and the
      !! caller must leave the column untouched.
      !!
      !! When grounded, each sub-floor layer is inflated to `floor` and the
      !! required volume is removed from the surplus layers (`h_old > floor`)
      !! in proportion to their surplus, so `Sigma h_new == Sigma h_old` and
      !! every returned `h_new(k) >= floor` (feasible whenever the column can
      !! hold `nz*floor`).  The largest-surplus layer absorbs the summation
      !! round-off so the column total is exact to a single ULP.
      !!
      !! Degenerate fallback: if the whole column cannot hold `nz*floor`
      !! (never reached in a deep isopycnal ocean), the mass is spread
      !! uniformly — still conservative, floor not guaranteed.
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_old(nz)
      real(wp), intent(in) :: h_floor
      real(wp), intent(out) :: h_new(nz)
      logical, intent(out) :: grounded

      integer :: k, k_donor
      real(wp) :: total, surplus_sum, deficit_sum, surplus_k
      real(wp) :: give, sum_new, best_surplus

      grounded = .false.
      total = 0.0_wp
      do k = 1, nz
         h_new(k) = h_old(k)
         total = total + h_old(k)
         if (h_old(k) < h_floor) grounded = .true.
      end do
      if (.not. grounded) return

      ! Degenerate: not enough water to floor every layer -> spread uniformly.
      if (total <= real(nz, wp)*h_floor) then
         do k = 1, nz
            h_new(k) = total/real(nz, wp)
         end do
         return
      end if

      ! Surplus available above the floor and total deficit below it.
      surplus_sum = 0.0_wp
      deficit_sum = 0.0_wp
      do k = 1, nz
         surplus_k = h_old(k) - h_floor
         if (surplus_k > 0.0_wp) then
            surplus_sum = surplus_sum + surplus_k
         else
            deficit_sum = deficit_sum - surplus_k
         end if
      end do
      ! total > nz*floor guarantees surplus_sum >= deficit_sum (feasible).

      ! Inflate sub-floor layers to the floor; draw the deficit from surplus
      ! layers proportionally to their surplus.
      do k = 1, nz
         surplus_k = h_old(k) - h_floor
         if (surplus_k > 0.0_wp) then
            give = (surplus_k/surplus_sum)*deficit_sum
            h_new(k) = h_old(k) - give
         else
            h_new(k) = h_floor
         end if
      end do

      ! Absorb summation round-off into the largest-surplus layer so the
      ! column total is preserved to a single ULP (exact conservation of h).
      k_donor = 1
      best_surplus = -1.0_wp
      sum_new = 0.0_wp
      do k = 1, nz
         sum_new = sum_new + h_new(k)
         surplus_k = h_old(k) - h_floor
         if (surplus_k > best_surplus) then
            best_surplus = surplus_k
            k_donor = k
         end if
      end do
      h_new(k_donor) = h_new(k_donor) + (total - sum_new)
   end subroutine min_thickness_target_column

   subroutine ocean_apply_conservative_min_thickness(grid, ms, h_new, grounded_mask, &
                                                     h_floor)
      !! Apply the conservative minimum-thickness adjustment in place on `ms`.
      !!
      !! `h_new` is caller-owned device-resident scratch shaped
      !! `(nx_total, ny_total, nz_ml)` (a `scratch_3d_buffer_t%data` slot);
      !! `grounded_mask` is `(nx_total, ny_total, 1)` scratch (`ct%mt_grounded`).
      !! On grounded columns `h_new` receives the floor-only target thickness;
      !! non-grounded columns are never written (their target ≡ h_old).
      !!
      !! Order (each step a separate device launch so the h_old reads all
      !! complete before h_layer is overwritten):
      !!   0. build the 2D grounded mask + global count (ONE coalesced pass);
      !!      count == 0 ⇒ return — the entire adjustment is a byte-no-op.
      !!   1. build the target field h_new on grounded columns,
      !!   2. remap every tracer   (h_old -> h_new) on grounded columns,
      !!   3. remap the face velocities (h_old -> h_new) on active faces,
      !!   4. h_layer := h_new on grounded columns.
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(inout) :: h_new(grid%nx_total, grid%ny_total, ms%nz_ml)
      real(wp), intent(inout) :: grounded_mask(grid%nx_total, grid%ny_total, 1)
      real(wp), intent(in) :: h_floor
         !! Minimum layer thickness (m); the isopycnal `angstrom_h`.

      integer :: nx, ny, nz, t, n_grounded

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml

      call build_grounded_mask(nx, ny, nz, ms%h_layer, h_floor, grounded_mask, &
                               n_grounded)
      if (n_grounded == 0) return

      call build_target_field(nx, ny, nz, ms%h_layer, h_floor, grounded_mask, h_new)

      if (allocated(ms%tracers)) then
         do t = 1, size(ms%tracers)
            if (.not. allocated(ms%tracers(t)%hTr)) cycle
            call remap_tracer_grounded(nx, ny, nz, ms%h_layer, h_new, grounded_mask, &
                                       ms%tracers(t)%hTr)
         end do
      end if

      call remap_x_face_grounded(nx, ny, nz, ms%h_layer, h_new, grounded_mask, &
                                 ms%u_face_x_layer)
      call remap_y_face_grounded(nx, ny, nz, ms%h_layer, h_new, grounded_mask, &
                                 ms%v_face_y_layer)

      call assign_h_layer(nx, ny, nz, h_new, grounded_mask, ms%h_layer)
   end subroutine ocean_apply_conservative_min_thickness

   subroutine build_grounded_mask(nx, ny, nz, h_old, h_floor, mask, n_grounded)
      !! ONE coalesced pass: mask(i,j,1) = 1.0 iff any layer of column (i,j)
      !! is strictly below the floor; `n_grounded` counts them (explicit
      !! OpenACC reduction — a `sum()` on a present-mapped array would run
      !! host-side under NVHPC non-managed mode and read the stale shadow).
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: h_old(nx, ny, nz)
      real(wp), intent(in) :: h_floor
      real(wp), intent(out) :: mask(nx, ny, 1)
      integer, intent(out) :: n_grounded

      integer :: i, j, k, n
      logical :: g

      n = 0
      do concurrent(j=1:ny, i=1:nx) local(g, k) reduce(+:n)
         g = .false.
         do k = 1, nz
            if (h_old(i, j, k) < h_floor) g = .true.
         end do
         if (g) then
            mask(i, j, 1) = 1.0_wp
            n = n + 1
         else
            mask(i, j, 1) = 0.0_wp
         end if
      end do
      n_grounded = n
   end subroutine build_grounded_mask

   pure subroutine build_target_field(nx, ny, nz, h_old, h_floor, mask, h_new)
      !! Flat-impl: per grounded (i,j) column build the floor-only conservative
      !! target.  Non-grounded columns are SKIPPED — their h_new is never read
      !! downstream (every consumer falls back to h_old via the mask).
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: h_old(nx, ny, nz)
      real(wp), intent(in) :: h_floor
      real(wp), intent(in) :: mask(nx, ny, 1)
      real(wp), intent(out) :: h_new(nx, ny, nz)
      integer :: i, j, k
      real(wp) :: h_old_col(NZ_STACK_MAX), h_new_col(NZ_STACK_MAX)
      logical :: grounded

      do concurrent(j=1:ny, i=1:nx) local(k, h_old_col, h_new_col, grounded)
         if (mask(i, j, 1) > 0.5_wp) then
            do k = 1, nz
               h_old_col(k) = h_old(i, j, k)
            end do
            call min_thickness_target_column(nz, h_old_col(1:nz), h_floor, &
                                             h_new_col(1:nz), grounded)
            do k = 1, nz
               h_new(i, j, k) = h_new_col(k)
            end do
         end if
      end do
   end subroutine build_target_field

   pure subroutine remap_tracer_grounded(nx, ny, nz, h_old, h_new, mask, hTr)
      !! Conservative tracer remap gated on the grounded mask.  Non-grounded
      !! columns are skipped -> hTr byte-unchanged; the concentration
      !! divisions only happen on grounded columns.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: h_old(nx, ny, nz)
      real(wp), intent(in) :: h_new(nx, ny, nz)
      real(wp), intent(in) :: mask(nx, ny, 1)
      real(wp), intent(inout) :: hTr(nx, ny, nz)
      integer :: i, j, k
      real(wp) :: h_old_col(NZ_STACK_MAX), h_new_col(NZ_STACK_MAX)
      real(wp) :: c_old(NZ_STACK_MAX), c_new(NZ_STACK_MAX)

      do concurrent(j=1:ny, i=1:nx) &
         local(k, h_old_col, h_new_col, c_old, c_new)
         if (mask(i, j, 1) > 0.5_wp) then
            do k = 1, nz
               h_old_col(k) = h_old(i, j, k)
               h_new_col(k) = h_new(i, j, k)
               if (h_old_col(k) > H_CONC_FLOOR) then
                  c_old(k) = hTr(i, j, k)/h_old_col(k)
               else
                  c_old(k) = 0.0_wp
               end if
            end do
            call remap_column(REMAP_PPM, nz, h_old_col(1:nz), h_new_col(1:nz), &
                              c_old(1:nz), c_new(1:nz))
            do k = 1, nz
               hTr(i, j, k) = c_new(k)*h_new_col(k)
            end do
         end if
      end do
   end subroutine remap_tracer_grounded

   pure subroutine remap_x_face_grounded(nx, ny, nz, h_old, h_new, mask, u_face_x)
      !! Conservative east-face velocity remap gated on the grounded mask.  A
      !! face is active iff either adjacent cell is grounded (2 mask loads —
      !! no neighbour-column re-scan); inactive faces are left byte-unchanged.
      !! Face thickness is the arithmetic mean of the two adjacent cells
      !! (outer walls take the single interior cell); the target side reads
      !! h_new only where the mask is set — elsewhere h_new ≡ h_old by
      !! construction, so h_old is read directly.  `remap_column` preserves
      !! the per-face momentum Sigma h_face . u.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: h_old(nx, ny, nz)
      real(wp), intent(in) :: h_new(nx, ny, nz)
      real(wp), intent(in) :: mask(nx, ny, 1)
      real(wp), intent(inout) :: u_face_x(nx + 1, ny, nz)
      integer :: I, j, k
      real(wp) :: h_old_face(NZ_STACK_MAX), h_new_face(NZ_STACK_MAX)
      real(wp) :: u_old(NZ_STACK_MAX), u_new(NZ_STACK_MAX)
      real(wp) :: hnw, hne
      logical :: active, gw, ge

      do concurrent(j=1:ny, I=1:nx + 1) &
         local(k, h_old_face, h_new_face, u_old, u_new, active, gw, ge, hnw, hne)
         gw = .false.
         ge = .false.
         if (I >= 2) gw = mask(I - 1, j, 1) > 0.5_wp
         if (I <= nx) ge = mask(I, j, 1) > 0.5_wp
         active = gw .or. ge
         if (active) then
            do k = 1, nz
               if (I == 1) then
                  h_old_face(k) = h_old(1, j, k)
                  h_new_face(k) = merge(h_new(1, j, k), h_old(1, j, k), ge)
               else if (I == nx + 1) then
                  h_old_face(k) = h_old(nx, j, k)
                  h_new_face(k) = merge(h_new(nx, j, k), h_old(nx, j, k), gw)
               else
                  h_old_face(k) = 0.5_wp*(h_old(I - 1, j, k) + h_old(I, j, k))
                  hnw = merge(h_new(I - 1, j, k), h_old(I - 1, j, k), gw)
                  hne = merge(h_new(I, j, k), h_old(I, j, k), ge)
                  h_new_face(k) = 0.5_wp*(hnw + hne)
               end if
               u_old(k) = u_face_x(I, j, k)
            end do
            call remap_column(REMAP_PPM, nz, h_old_face(1:nz), h_new_face(1:nz), &
                              u_old(1:nz), u_new(1:nz))
            do k = 1, nz
               u_face_x(I, j, k) = u_new(k)
            end do
         end if
      end do
   end subroutine remap_x_face_grounded

   pure subroutine remap_y_face_grounded(nx, ny, nz, h_old, h_new, mask, v_face_y)
      !! Conservative north-face velocity remap, mirror of the x routine.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: h_old(nx, ny, nz)
      real(wp), intent(in) :: h_new(nx, ny, nz)
      real(wp), intent(in) :: mask(nx, ny, 1)
      real(wp), intent(inout) :: v_face_y(nx, ny + 1, nz)
      integer :: i, J, k
      real(wp) :: h_old_face(NZ_STACK_MAX), h_new_face(NZ_STACK_MAX)
      real(wp) :: v_old(NZ_STACK_MAX), v_new(NZ_STACK_MAX)
      real(wp) :: hns, hnn
      logical :: active, gs, gn

      do concurrent(J=1:ny + 1, i=1:nx) &
         local(k, h_old_face, h_new_face, v_old, v_new, active, gs, gn, hns, hnn)
         gs = .false.
         gn = .false.
         if (J >= 2) gs = mask(i, J - 1, 1) > 0.5_wp
         if (J <= ny) gn = mask(i, J, 1) > 0.5_wp
         active = gs .or. gn
         if (active) then
            do k = 1, nz
               if (J == 1) then
                  h_old_face(k) = h_old(i, 1, k)
                  h_new_face(k) = merge(h_new(i, 1, k), h_old(i, 1, k), gn)
               else if (J == ny + 1) then
                  h_old_face(k) = h_old(i, ny, k)
                  h_new_face(k) = merge(h_new(i, ny, k), h_old(i, ny, k), gs)
               else
                  h_old_face(k) = 0.5_wp*(h_old(i, J - 1, k) + h_old(i, J, k))
                  hns = merge(h_new(i, J - 1, k), h_old(i, J - 1, k), gs)
                  hnn = merge(h_new(i, J, k), h_old(i, J, k), gn)
                  h_new_face(k) = 0.5_wp*(hns + hnn)
               end if
               v_old(k) = v_face_y(i, J, k)
            end do
            call remap_column(REMAP_PPM, nz, h_old_face(1:nz), h_new_face(1:nz), &
                              v_old(1:nz), v_new(1:nz))
            do k = 1, nz
               v_face_y(i, J, k) = v_new(k)
            end do
         end if
      end do
   end subroutine remap_y_face_grounded

   pure subroutine assign_h_layer(nx, ny, nz, h_new, mask, h_layer)
      !! Copy the target field into h_layer on grounded columns only —
      !! non-grounded columns' h_new was never written, and their old
      !! h_layer is already the exact target.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: h_new(nx, ny, nz)
      real(wp), intent(in) :: mask(nx, ny, 1)
      real(wp), intent(inout) :: h_layer(nx, ny, nz)
      integer :: i, j, k
      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         if (mask(i, j, 1) > 0.5_wp) then
            h_layer(i, j, k) = h_new(i, j, k)
         end if
      end do
   end subroutine assign_h_layer

end module rdb_ocean_min_thickness
