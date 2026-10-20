!! Periodic-boundary ghost-wrap helpers for the ocean dyn-core.
module rdb_ocean_periodic
   !! GPU-resident periodic ghost-wrap helpers. Seam invariant after every
   !! wrap call: (a) every ghost cell equals its interior partner, and
   !! (b) u(i_w, ·) == u(i_e, ·) bit-for-bit (i_w = nghost+1,
   !! i_e = nghost+nx_phys+1).
   !!
   !! Explicit-shape dummies (never assumed-shape), j-outer / i-inner
   !! `do concurrent` for NVHPC GPU coalescing. The fast barotropic substep
   !! does NOT call these — it wraps inline in its `!$acc kernels async(1)`
   !! region.
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_boundary_types, only: ocean_bc_state_t
   ! NOTE (D4): this kernel module must NOT import rdb_ocean_halo — that would
   ! create a USE cycle (rdb_ocean_halo uses these wrap kernels).  The
   ! decomposed-axis skip is passed IN by the caller via the optional
   ! skip_x/skip_y arguments below, computed at the (higher-level) call site
   ! from ocean_halo_is_decomposed_x/y().
   implicit none
   private

   public :: ocean_periodic_wrap_centre_2d
   public :: ocean_periodic_wrap_centre_3d
   public :: ocean_periodic_wrap_face_x_2d
   public :: ocean_periodic_wrap_face_x_3d
   public :: ocean_periodic_wrap_face_y_2d
   public :: ocean_periodic_wrap_face_y_3d
   public :: ocean_periodic_wrap_state

contains

   pure subroutine ocean_periodic_wrap_centre_2d(fld, nx_total, ny_total, &
                                                 nx_phys, ny_phys, nghost, &
                                                 wrap_x, wrap_y)
      !! Fill ghost cells of a cell-centred 2D field (e.g. η, bt_H_ref)
      !! with the periodically-matching interior values.
      !! Explicit-shape dummies avoid per-launch descriptor-walk memcpys.
      integer, intent(in) :: nx_total, ny_total, nx_phys, ny_phys, nghost
      real(wp), intent(inout) :: fld(nx_total, ny_total)
         !! Cell-centred field, shape (nx_total, ny_total).
      logical, intent(in) :: wrap_x
         !! Wrap ghost columns (west ↔ east).
      logical, intent(in) :: wrap_y
         !! Wrap ghost rows (south ↔ north).

      integer :: i, j

      ! X-wrap first, then Y-wrap in a separate loop.  The two passes must
      ! not share a single do-concurrent because the Y-wrap reads ghost
      ! columns that the X-wrap writes (data dependency across iterations).
      if (wrap_x) then
         do concurrent(j=1:ny_total, i=1:nx_total)
            ! West ghosts: fld(1..nghost, j) := fld(nx_phys+1..nx_phys+nghost, j)
            if (i <= nghost) then
               fld(i, j) = fld(i + nx_phys, j)
            end if
            ! East ghosts: fld(nx_phys+nghost+1..nx_total, j) := fld(nghost+1..2*nghost, j)
            if (i > nx_phys + nghost) then
               fld(i, j) = fld(i - nx_phys, j)
            end if
         end do
      end if
      if (wrap_y) then
         do concurrent(j=1:ny_total, i=1:nx_total)
            ! South ghosts: fld(i, 1..nghost) := fld(i, ny_phys+1..ny_phys+nghost)
            if (j <= nghost) then
               fld(i, j) = fld(i, j + ny_phys)
            end if
            ! North ghosts: fld(i, ny_phys+nghost+1..ny_total) := fld(i, nghost+1..2*nghost)
            if (j > ny_phys + nghost) then
               fld(i, j) = fld(i, j - ny_phys)
            end if
         end do
      end if
   end subroutine ocean_periodic_wrap_centre_2d

   subroutine ocean_periodic_wrap_centre_3d(fld, nx_total, ny_total, nz, &
                                            nx_phys, ny_phys, nghost, &
                                            wrap_x, wrap_y, no_wait)
      !! Fill ghost cells of a cell-centred 3D field (e.g. h_layer, hTr, T, S).
      !! `no_wait` (optional, default .false.): when .true., loops issue on
      !! OpenACC queue 1 and the routine returns WITHOUT syncing, so a batched
      !! caller can pipeline many tiny ghost-slab wraps and `!$acc wait(1)`
      !! once. Default ⇒ self-contained blocking wrap. Not `pure` (directives).
      integer, intent(in) :: nx_total, ny_total, nz, nx_phys, ny_phys, nghost
      real(wp), intent(inout) :: fld(nx_total, ny_total, nz)
         !! Cell-centred field, shape (nx_total, ny_total, nz).
      logical, intent(in) :: wrap_x
      logical, intent(in) :: wrap_y
      logical, intent(in), optional :: no_wait

      integer :: i, j, k
      logical :: lwait

      lwait = .true.   ! default: blocking (wait) — safe for non-batched callers
      if (present(no_wait)) lwait = .not. no_wait

      ! X-wrap first, then Y-wrap in a separate loop.  Two passes needed
      ! because the Y-wrap reads the x-ghost columns the X-pass just wrote.
      ! Same queue ⇒ ordered ⇒ the X→Y dependency holds.
      !$acc kernels async(1)
      if (wrap_x) then
         do concurrent(k=1:nz, j=1:ny_total, i=1:nx_total)
            if (i <= nghost) then
               fld(i, j, k) = fld(i + nx_phys, j, k)
            end if
            if (i > nx_phys + nghost) then
               fld(i, j, k) = fld(i - nx_phys, j, k)
            end if
         end do
      end if
      if (wrap_y) then
         do concurrent(k=1:nz, j=1:ny_total, i=1:nx_total)
            if (j <= nghost) then
               fld(i, j, k) = fld(i, j + ny_phys, k)
            end if
            if (j > ny_phys + nghost) then
               fld(i, j, k) = fld(i, j - ny_phys, k)
            end if
         end do
      end if
      !$acc end kernels
      if (lwait) then
         !$acc wait(1)
      end if
   end subroutine ocean_periodic_wrap_centre_3d

   pure subroutine ocean_periodic_wrap_face_x_2d(fld, nx_face, ny_total, &
                                                 nx_phys, ny_phys, nghost, &
                                                 wrap_x, wrap_y)
      !! Fill ghost faces of a 2D x-face field (e.g. bt_ubt, shape nx_total+1).
      !! Also copies the west physical-wall face value onto the east physical-wall
      !! face (belt-and-braces seam invariant, §1.1 property b).
      !! nx_face = nx_total + 1.
      integer, intent(in) :: nx_face, ny_total, nx_phys, ny_phys, nghost
      real(wp), intent(inout) :: fld(nx_face, ny_total)
         !! x-face field, shape (nx_total+1, ny_total).
      logical, intent(in) :: wrap_x
      logical, intent(in) :: wrap_y

      integer :: i, j
      integer :: i_w, i_e

      i_w = nghost + 1          ! west physical-wall face index
      i_e = nghost + nx_phys + 1  ! east physical-wall face index

      ! X-wrap (including belt-and-braces seam copy) in its own DC loop,
      ! then Y-wrap in a separate loop to avoid cross-iteration data dependencies.
      if (wrap_x) then
         do concurrent(j=1:ny_total, i=1:nx_face)
            ! West ghost faces i=1..nghost  ← interior face i+nx_phys
            if (i <= nghost) then
               fld(i, j) = fld(i + nx_phys, j)
            end if
            ! East ghost faces i=nx_phys+nghost+2..nx_face ← interior face i-nx_phys
            if (i > nx_phys + nghost + 1) then
               fld(i, j) = fld(i - nx_phys, j)
            end if
            ! Belt-and-braces: copy west wall face onto east wall face.
            ! Ensures property (b) even if a kernel breaks expression symmetry.
            if (i == i_e) fld(i, j) = fld(i_w, j)
         end do
      end if
      if (wrap_y) then
         do concurrent(j=1:ny_total, i=1:nx_face)
            if (j <= nghost) then
               fld(i, j) = fld(i, j + ny_phys)
            end if
            if (j > ny_phys + nghost) then
               fld(i, j) = fld(i, j - ny_phys)
            end if
         end do
      end if
   end subroutine ocean_periodic_wrap_face_x_2d

   subroutine ocean_periodic_wrap_face_x_3d(fld, nx_face, ny_total, nz, &
                                            nx_phys, ny_phys, nghost, &
                                            wrap_x, wrap_y, no_wait)
      !! Fill ghost faces of a 3D x-face field (e.g. u_face_x_layer).
      !! nx_face = nx_total + 1.  `no_wait` (optional): see
      !! `ocean_periodic_wrap_centre_3d` — batched async(1), sync once at the
      !! caller.  Not `pure` (async/wait directives).
      integer, intent(in) :: nx_face, ny_total, nz, nx_phys, ny_phys, nghost
      real(wp), intent(inout) :: fld(nx_face, ny_total, nz)
         !! x-face field, shape (nx_total+1, ny_total, nz).
      logical, intent(in) :: wrap_x
      logical, intent(in) :: wrap_y
      logical, intent(in), optional :: no_wait

      integer :: i, j, k
      integer :: i_w, i_e
      logical :: lwait

      lwait = .true.   ! default: blocking (wait) — safe for non-batched callers
      if (present(no_wait)) lwait = .not. no_wait
      i_w = nghost + 1
      i_e = nghost + nx_phys + 1

      ! Two-pass: X-wrap (+ belt-and-braces) first, Y-wrap second.
      !$acc kernels async(1)
      if (wrap_x) then
         do concurrent(k=1:nz, j=1:ny_total, i=1:nx_face)
            if (i <= nghost) then
               fld(i, j, k) = fld(i + nx_phys, j, k)
            end if
            if (i > nx_phys + nghost + 1) then
               fld(i, j, k) = fld(i - nx_phys, j, k)
            end if
            if (i == i_e) fld(i, j, k) = fld(i_w, j, k)
         end do
      end if
      if (wrap_y) then
         do concurrent(k=1:nz, j=1:ny_total, i=1:nx_face)
            if (j <= nghost) then
               fld(i, j, k) = fld(i, j + ny_phys, k)
            end if
            if (j > ny_phys + nghost) then
               fld(i, j, k) = fld(i, j - ny_phys, k)
            end if
         end do
      end if
      !$acc end kernels
      if (lwait) then
         !$acc wait(1)
      end if
   end subroutine ocean_periodic_wrap_face_x_3d

   pure subroutine ocean_periodic_wrap_face_y_2d(fld, nx_total, ny_face, &
                                                 nx_phys, ny_phys, nghost, &
                                                 wrap_x, wrap_y)
      !! Fill ghost faces of a 2D y-face field (e.g. bt_vbt, shape ny_total+1).
      !! Also copies south physical-wall face onto north physical-wall face.
      !! ny_face = ny_total + 1.
      integer, intent(in) :: nx_total, ny_face, nx_phys, ny_phys, nghost
      real(wp), intent(inout) :: fld(nx_total, ny_face)
         !! y-face field, shape (nx_total, ny_total+1).
      logical, intent(in) :: wrap_x
      logical, intent(in) :: wrap_y

      integer :: i, j
      integer :: j_s, j_n

      j_s = nghost + 1
      j_n = nghost + ny_phys + 1

      ! Two-pass: X-wrap first, Y-wrap (+ belt-and-braces) second.
      if (wrap_x) then
         do concurrent(j=1:ny_face, i=1:nx_total)
            if (i <= nghost) then
               fld(i, j) = fld(i + nx_phys, j)
            end if
            if (i > nx_phys + nghost) then
               fld(i, j) = fld(i - nx_phys, j)
            end if
         end do
      end if
      if (wrap_y) then
         do concurrent(j=1:ny_face, i=1:nx_total)
            ! South ghost faces j=1..nghost ← interior face j+ny_phys
            if (j <= nghost) then
               fld(i, j) = fld(i, j + ny_phys)
            end if
            ! North ghost faces j=ny_phys+nghost+2..ny_face ← interior face j-ny_phys
            if (j > ny_phys + nghost + 1) then
               fld(i, j) = fld(i, j - ny_phys)
            end if
            ! Belt-and-braces: copy south wall face onto north wall face.
            if (j == j_n) fld(i, j) = fld(i, j_s)
         end do
      end if
   end subroutine ocean_periodic_wrap_face_y_2d

   subroutine ocean_periodic_wrap_face_y_3d(fld, nx_total, ny_face, nz, &
                                            nx_phys, ny_phys, nghost, &
                                            wrap_x, wrap_y, no_wait)
      !! Fill ghost faces of a 3D y-face field (e.g. v_face_y_layer).
      !! ny_face = ny_total + 1.  `no_wait` (optional): see
      !! `ocean_periodic_wrap_centre_3d`.  Not `pure` (async/wait directives).
      integer, intent(in) :: nx_total, ny_face, nz, nx_phys, ny_phys, nghost
      real(wp), intent(inout) :: fld(nx_total, ny_face, nz)
         !! y-face field, shape (nx_total, ny_total+1, nz).
      logical, intent(in) :: wrap_x
      logical, intent(in) :: wrap_y
      logical, intent(in), optional :: no_wait

      integer :: i, j, k
      integer :: j_s, j_n
      logical :: lwait

      lwait = .true.
      if (present(no_wait)) lwait = .not. no_wait
      j_s = nghost + 1
      j_n = nghost + ny_phys + 1

      ! Two-pass: X-wrap first, Y-wrap (+ belt-and-braces) second.
      if (wrap_x) then
         !$acc kernels async(1)
         do concurrent(k=1:nz, j=1:ny_face, i=1:nx_total)
            if (i <= nghost) then
               fld(i, j, k) = fld(i + nx_phys, j, k)
            end if
            if (i > nx_phys + nghost) then
               fld(i, j, k) = fld(i - nx_phys, j, k)
            end if
         end do
         !$acc end kernels
      end if
      if (wrap_y) then
         !$acc kernels async(1)
         do concurrent(k=1:nz, j=1:ny_face, i=1:nx_total)
            if (j <= nghost) then
               fld(i, j, k) = fld(i, j + ny_phys, k)
            end if
            if (j > ny_phys + nghost + 1) then
               fld(i, j, k) = fld(i, j - ny_phys, k)
            end if
            if (j == j_n) fld(i, j, k) = fld(i, j_s, k)
         end do
         !$acc end kernels
      end if
      if (lwait) then
         !$acc wait(1)
      end if
   end subroutine ocean_periodic_wrap_face_y_3d

   subroutine ocean_periodic_wrap_state(grid, bc, ms, skip_x, skip_y)
      !! Convenience wrapper: wrap h_layer, u/v layer faces, and every
      !! registered tracer.  Called at stage entry (before
      !! `derive_bt_from_layers`) and after continuity (before hdiff).
      !!
      !! Per-tracer loop is OUTSIDE the DC kernels (outer-shim pattern:
      !! array-of-derived-types cannot be dereferenced on device).
      !!
      !! No-op when neither `bc%periodic_x` nor `bc%periodic_y` is set.
      !!
      !! D4: `skip_x` / `skip_y` (optional, default .false.): when .true.,
      !! the local wrap on that axis is suppressed.  The caller passes
      !! `ocean_halo_is_decomposed_x/y()` so that a multi-rank periodic-x
      !! decomposition does not double-wrap the seam (the MPI halo already
      !! owns those ghost columns).  Single-rank ⇒ caller passes .false.
      !! (or omits) ⇒ bit-identical.
      !! NOTE: rdb_ocean_periodic must NOT import rdb_ocean_halo (USE cycle —
      !! the halo backends import this module's wrap kernels).  The skip is
      !! passed in by the caller.
      type(hgrid_t), intent(in) :: grid
      type(ocean_bc_state_t), intent(in) :: bc
      type(multilayer_state_t), intent(inout) :: ms
      logical, intent(in), optional :: skip_x
         !! When .true., suppress the local x-axis wrap (D4 multi-rank).
      logical, intent(in), optional :: skip_y
         !! When .true., suppress the local y-axis wrap (D4 multi-rank).

      integer :: it
      integer :: nx, ny, nz, nx_phys, ny_phys, nghost
      logical :: per_x, per_y

      per_x = bc%periodic_x
      per_y = bc%periodic_y
      if (.not. per_x .and. .not. per_y) return

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml
      nx_phys = grid%nx_phys
      ny_phys = grid%ny_phys
      nghost = grid%nghost

      ! D4: skip the LOCAL periodic wrap on any axis split across MPI ranks.
      ! On a decomposed axis the ghost columns are owned by a NEIGHBOUR rank —
      ! the local wrap kernel copies the wrong (local-subdomain) columns into
      ! those ghosts, overwriting what ocean_halo_* just filled correctly.
      ! The MPI halo already ran before this call, so the ghosts are correct;
      ! just leave them alone.  Single-rank ⇒ skip_x absent/.false.
      ! ⇒ local wrap governs as before (bit-identical).
      if (present(skip_x)) per_x = per_x .and. .not. skip_x
      if (present(skip_y)) per_y = per_y .and. .not. skip_y
      if (.not. per_x .and. .not. per_y) return

      ! Batched async wrap: every per-field wrap issues on queue 1
      ! (no_wait=.true.); the whole batch is synced ONCE below. Pipelines the
      ! tiny launch-latency-bound ghost-slab kernels; same queue ⇒ ordered.
      call ocean_periodic_wrap_centre_3d(ms%h_layer, nx, ny, nz, &
                                         nx_phys, ny_phys, nghost, per_x, per_y, no_wait=.true.)
      call ocean_periodic_wrap_face_x_3d(ms%u_face_x_layer, nx + 1, ny, nz, &
                                         nx_phys, ny_phys, nghost, per_x, per_y, no_wait=.true.)
      call ocean_periodic_wrap_face_y_3d(ms%v_face_y_layer, nx, ny + 1, nz, &
                                         nx_phys, ny_phys, nghost, per_x, per_y, no_wait=.true.)

      ! Per-tracer loop outside DCs — outer-shim pattern for array-of-DTs.
      if (allocated(ms%tracers)) then
         do it = 1, size(ms%tracers)
            if (.not. allocated(ms%tracers(it)%hTr)) cycle
            call ocean_periodic_wrap_centre_3d(ms%tracers(it)%hTr, nx, ny, nz, &
                                               nx_phys, ny_phys, nghost, per_x, per_y, no_wait=.true.)
         end do
      end if

      ! Single sync for the whole batch (no-op on non-OpenACC builds).
      !$acc wait(1)
   end subroutine ocean_periodic_wrap_state

end module rdb_ocean_periodic
