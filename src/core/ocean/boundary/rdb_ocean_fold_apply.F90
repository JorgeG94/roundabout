!! State-level tripolar north-fold seam application for the ocean dyn-core.
module rdb_ocean_fold_apply
   !! Orchestration over the pure seam operators in `rdb_ocean_fold`:
   !! dereferences the state slots (outer-shim per-tracer loop) and calls
   !! the explicit-shape fold kernels.
   !!
   !! Ordering contract: periodic-x is wrapped FIRST, fold SECOND — every
   !! fold routine is called AFTER the matching periodic wrap so it reads
   !! the already cyclically-wrapped corner columns.
   !!
   !! Stagger map:
   !!   * h_layer, η, tracer hTr  → centre fold (copy, no sign flip)
   !!   * u_face_x_layer, bt_ubt  → u-face fold (negate — true vector)
   !!   * v_face_y_layer, bt_vbt  → v-face fold (negate + on-row
   !!     antisymmetric projection of the fold line, storage row
   !!     nghost+ny_phys+1 — see the `rdb_ocean_fold` header)
   !!
   !! Every routine no-ops when `bc%north_fold` is .false. ⇒ non-tripolar
   !! runs stay bit-identical.
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_boundary_types, only: ocean_bc_state_t
   use rdb_ocean_fold, only: fold_north_centre, fold_north_u_face, &
                             fold_north_v_face
   implicit none
   private

   public :: ocean_fold_wrap_state
   public :: ocean_fold_wrap_centre_3d_state
   public :: ocean_fold_wrap_eta_2d

contains

   subroutine ocean_fold_wrap_state(grid, bc, ms)
      !! Fold the north seam of h_layer, u/v layer faces, and every
      !! registered tracer.  Call AFTER `ocean_periodic_wrap_state`.
      !! No-op when `bc%north_fold` is .false.
      type(hgrid_t), intent(in) :: grid
      type(ocean_bc_state_t), intent(in) :: bc
      type(multilayer_state_t), intent(inout) :: ms

      integer :: it
      integer :: nx, ny, nz, nx_phys, ny_phys, nghost

      if (.not. bc%north_fold) return

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml
      nx_phys = grid%nx_phys
      ny_phys = grid%ny_phys
      nghost = grid%nghost

      ! Centre (T): h_layer.
      call fold_north_centre(ms%h_layer, nx, ny, nz, &
                             nx_phys, ny_phys, nghost)
      ! u-face (Cu): negate.
      call fold_north_u_face(ms%u_face_x_layer, nx + 1, ny, nz, &
                             nx_phys, ny_phys, nghost)
      ! v-face (Cv): negate + on-row antisymmetric projection.
      call fold_north_v_face(ms%v_face_y_layer, nx, ny + 1, nz, &
                             nx_phys, ny_phys, nghost)

      ! Per-tracer loop OUTSIDE the DC kernels (outer-shim for array-of-DTs).
      if (allocated(ms%tracers)) then
         do it = 1, size(ms%tracers)
            if (.not. allocated(ms%tracers(it)%hTr)) cycle
            call fold_north_centre(ms%tracers(it)%hTr, nx, ny, nz, &
                                   nx_phys, ny_phys, nghost)
         end do
      end if
   end subroutine ocean_fold_wrap_state

   subroutine ocean_fold_wrap_centre_3d_state(grid, bc, ms)
      !! Fold ONLY h_layer + tracers (centre fields) — the continuity
      !! mid-split site, which re-wraps the centre fields between the
      !! zonal and meridional Lie-split halves.  No-op when not folding.
      type(hgrid_t), intent(in) :: grid
      type(ocean_bc_state_t), intent(in) :: bc
      type(multilayer_state_t), intent(inout) :: ms

      integer :: it
      integer :: nx, ny, nz, nx_phys, ny_phys, nghost

      if (.not. bc%north_fold) return

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml
      nx_phys = grid%nx_phys
      ny_phys = grid%ny_phys
      nghost = grid%nghost

      call fold_north_centre(ms%h_layer, nx, ny, nz, &
                             nx_phys, ny_phys, nghost)
      if (allocated(ms%tracers)) then
         do it = 1, size(ms%tracers)
            if (.not. allocated(ms%tracers(it)%hTr)) cycle
            call fold_north_centre(ms%tracers(it)%hTr, nx, ny, nz, &
                                   nx_phys, ny_phys, nghost)
         end do
      end if
   end subroutine ocean_fold_wrap_centre_3d_state

   subroutine ocean_fold_wrap_eta_2d(grid, bc, eta)
      !! Fold a 2D cell-centred η field (driver-level SSH wrap site).
      !! No-op when not folding.
      type(hgrid_t), intent(in) :: grid
      type(ocean_bc_state_t), intent(in) :: bc
      real(wp), intent(inout) :: eta(:, :)

      if (.not. bc%north_fold) return

      call fold_north_centre(eta, grid%nx_total, grid%ny_total, &
                             grid%nx_phys, grid%ny_phys, grid%nghost)
   end subroutine ocean_fold_wrap_eta_2d

end module rdb_ocean_fold_apply
