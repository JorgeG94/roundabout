!! Tripolar north-fold seam-exchange helpers for the ocean dyn-core.
module rdb_ocean_fold
   !! Single-rank discrete tripolar north-fold exchange (Murray 1996). Pure
   !! seam operators only (state orchestration lives in `rdb_ocean_fold_apply`).
   !! Free procedures, explicit-shape dummies, j-outer / i-inner `do concurrent`.
   !!
   !! Seam geometry: the fold line is the Cv/Bu line at the top of T-row
   !! `nj` (= ny_phys). T and u(Cu) images are pure halo (rows `j > nj`);
   !! v(Cv) and corner(Bu) at `j = nj` lie ON the self-conjugate line. So an
   !! exchange is two operations: (1) halo-fill rows `j > nj` for T (copy)
   !! and u (copy + sign flip); (2) on-row antisymmetric projection for v
   !! and corners at `j = nj`.
   !!
   !! Index maps — physical `i ∈ 1..ni`, last T-row `j = nj`:
   !! | Stagger | i-map (phys)  | j-map (phys) | on-line |
   !! |---------|---------------|--------------|---------|
   !! | T       | i' = ni+1-i   | j' = 2nj-j+1 | no      |
   !! | u (Cu)  | i' = ni+2-f   | j' = 2nj-j+1 | no      |
   !! | v (Cv)  | i' = ni+1-i   | j' = 2nj-j   | YES j=nj|
   !! | corner  | i' = ni+2-c   | j' = 2nj-j   | YES j=nj|
   !! u/corner use `ni+2-f` (not `ni-i`) because face storage is symmetric:
   !! extent `nx_total+1`, `u_face_x(f)` is the WEST face of T-cell `f`.
   !!
   !! On-line projection (v, corners): the two i-halves at `j = nj` duplicate
   !! the same physical points; enforce `v(i,nj) = -v(i',nj)`, i' = ni+1-i,
   !! by overwriting the west half from the negated east mirror. The
   !! self-fixed column (odd ni only) is set to 0.
   !!
   !! Caller ordering (MANDATORY): periodic-x wrap FIRST, then the fold, so
   !! the fold reads already cyclically-wrapped ghost columns at the corners.
   !!
   !! All helpers stay `pure` + `do concurrent` so they run both on
   !! device-mapped arrays (stdpar) and during host setup of metric ghosts.
   use rdb_constants, only: wp
   implicit none
   private

   ! Rank-generic public API: each generic resolves at compile time to the
   ! 2D or 3D specific by array rank (static dispatch, GPU-safe, no vtable).
   public :: fold_north_centre
   public :: fold_north_u_face
   public :: fold_north_v_face
   public :: fold_north_corner

   interface fold_north_centre
      module procedure fold_north_centre_2d, fold_north_centre_3d
   end interface fold_north_centre

   interface fold_north_u_face
      module procedure fold_north_u_face_2d, fold_north_u_face_3d
   end interface fold_north_u_face

   interface fold_north_v_face
      module procedure fold_north_v_face_2d, fold_north_v_face_3d
   end interface fold_north_v_face

   interface fold_north_corner
      module procedure fold_north_corner_2d
   end interface fold_north_corner

contains

   ! ================================================================
   ! T-stagger (cell centre): pure halo-fill, rows j > nj.
   !   i' = ni+1-i   (storage: i' = 2*nghost+ni+1 - i)
   !   j' = 2nj-j+1  (storage: j' = 2*nghost+2*nj+1 - j)
   ! Scalars copy unchanged (no sign flip).
   ! ================================================================

   pure subroutine fold_north_centre_2d(fld, nx_total, ny_total, &
                                        nx_phys, ny_phys, nghost)
      !! Fill the north halo of a 2D cell-centred field by the T-fold.
      integer, intent(in) :: nx_total, ny_total, nx_phys, ny_phys, nghost
      real(wp), intent(inout) :: fld(nx_total, ny_total)
         !! Cell-centred field, shape (nx_total, ny_total).

      integer :: i, j, isum, jsum, j_lo

      isum = 2*nghost + nx_phys + 1
      jsum = 2*nghost + 2*ny_phys + 1
      j_lo = nghost + ny_phys + 1   ! first north halo row (storage)

      do concurrent(j=j_lo:ny_total, i=1:nx_total)
         fld(i, j) = fld(isum - i, jsum - j)
      end do
   end subroutine fold_north_centre_2d

   pure subroutine fold_north_centre_3d(fld, nx_total, ny_total, nz, &
                                        nx_phys, ny_phys, nghost)
      !! 3D T-fold halo-fill — identical per level.
      integer, intent(in) :: nx_total, ny_total, nz, nx_phys, ny_phys, nghost
      real(wp), intent(inout) :: fld(nx_total, ny_total, nz)
         !! Cell-centred 3D field, shape (nx_total, ny_total, nz).

      integer :: i, j, k, isum, jsum, j_lo

      isum = 2*nghost + nx_phys + 1
      jsum = 2*nghost + 2*ny_phys + 1
      j_lo = nghost + ny_phys + 1

      do concurrent(k=1:nz, j=j_lo:ny_total, i=1:nx_total)
         fld(i, j, k) = fld(isum - i, jsum - j, k)
      end do
   end subroutine fold_north_centre_3d

   ! ================================================================
   ! u-stagger (Cu, x-face, extent nx_total+1): halo-fill rows j > nj.
   !   f' = ni+2-f  (storage: f' = 2*nghost+ni+2 - f)  [sym storage]
   !   j' = 2nj-j+1
   ! True vector component → NEGATE across the fold.
   ! ================================================================

   pure subroutine fold_north_u_face_2d(u, nx_face, ny_total, &
                                        nx_phys, ny_phys, nghost)
      !! Fill the north halo of a 2D x-face (Cu) field; sign-flipped.
      integer, intent(in) :: nx_face, ny_total, nx_phys, ny_phys, nghost
      real(wp), intent(inout) :: u(nx_face, ny_total)
         !! x-face field, shape (nx_total+1, ny_total).

      integer :: i, j, fsum, jsum, j_lo

      fsum = 2*nghost + nx_phys + 2
      jsum = 2*nghost + 2*ny_phys + 1
      j_lo = nghost + ny_phys + 1

      do concurrent(j=j_lo:ny_total, i=1:nx_face)
         u(i, j) = -u(fsum - i, jsum - j)
      end do
   end subroutine fold_north_u_face_2d

   pure subroutine fold_north_u_face_3d(u, nx_face, ny_total, nz, &
                                        nx_phys, ny_phys, nghost)
      !! 3D x-face (Cu) north-halo fill; sign-flipped, per-level identical.
      integer, intent(in) :: nx_face, ny_total, nz, nx_phys, ny_phys, nghost
      real(wp), intent(inout) :: u(nx_face, ny_total, nz)
         !! x-face 3D field, shape (nx_total+1, ny_total, nz).

      integer :: i, j, k, fsum, jsum, j_lo

      fsum = 2*nghost + nx_phys + 2
      jsum = 2*nghost + 2*ny_phys + 1
      j_lo = nghost + ny_phys + 1

      do concurrent(k=1:nz, j=j_lo:ny_total, i=1:nx_face)
         u(i, j, k) = -u(fsum - i, jsum - j, k)
      end do
   end subroutine fold_north_u_face_3d

   ! ================================================================
   ! v-stagger (Cv, y-face): TWO operations.
   !   (1) halo-fill rows j > j_fold (reflected + negated):
   !         i' = ni+1-i, j' = 2nj-j  (storage j' = 2*nghost+2*nj - j)
   !   (2) ON-LINE projection at j = j_fold (= nghost+nj):
   !         v(i,j_fold) = -v(i', j_fold), i' = ni+1-i
   !       overwrite the WEST half from the negated east-mirror; self-fixed
   !       column (odd ni) → 0.
   ! ================================================================

   pure subroutine fold_north_v_face_2d(v, nx_total, ny_face, &
                                        nx_phys, ny_phys, nghost)
      !! 2D y-face (Cv) fold: north-halo fill (negated) + on-line
      !! antisymmetric projection at the fold row.  Used for the
      !! barotropic `bt_vbt` field in the BT fast loop.
      integer, intent(in) :: nx_total, ny_face, nx_phys, ny_phys, nghost
      real(wp), intent(inout) :: v(nx_total, ny_face)
         !! y-face 2D field, shape (nx_total, ny_total+1).

      integer :: i, j, isum, jsum, j_fold, i_lo, i_mid, ip

      isum = 2*nghost + nx_phys + 1
      jsum = 2*nghost + 2*ny_phys
      j_fold = nghost + ny_phys
      i_lo = nghost + 1

      ! (1) Halo rows strictly beyond the fold row.
      do concurrent(j=j_fold + 1:ny_face, i=1:nx_total)
         v(i, j) = -v(isum - i, jsum - j)
      end do

      ! (2) On-line antisymmetric projection at j = j_fold.
      i_mid = nghost + (nx_phys + 1)/2
      do concurrent(i=i_lo:i_mid)
         ip = isum - i
         if (ip == i) then
            v(i, j_fold) = 0.0_wp
         else
            v(i, j_fold) = -v(ip, j_fold)
         end if
      end do
   end subroutine fold_north_v_face_2d

   pure subroutine fold_north_v_face_3d(v, nx_total, ny_face, nz, &
                                        nx_phys, ny_phys, nghost)
      !! 3D y-face (Cv) fold: north-halo fill (negated) + on-line
      !! antisymmetric projection at the fold row.
      integer, intent(in) :: nx_total, ny_face, nz, nx_phys, ny_phys, nghost
      real(wp), intent(inout) :: v(nx_total, ny_face, nz)
         !! y-face 3D field, shape (nx_total, ny_total+1, nz).

      integer :: i, j, k, isum, jsum, j_fold, i_lo, i_mid, ip

      isum = 2*nghost + nx_phys + 1
      jsum = 2*nghost + 2*ny_phys
      j_fold = nghost + ny_phys          ! the self-conjugate fold row
      i_lo = nghost + 1                  ! first physical column

      ! (1) Halo rows strictly beyond the fold row.
      do concurrent(k=1:nz, j=j_fold + 1:ny_face, i=1:nx_total)
         v(i, j, k) = -v(isum - i, jsum - j, k)
      end do

      ! (2) On-line antisymmetric projection at j = j_fold.
      !     West half (i in i_lo .. i_mid) overwritten from negated mirror.
      !     i_mid is the middle physical column; for even ni it is the last
      !     west-half column (no fixed point), for odd ni it is the fixed
      !     column which is forced to 0.
      i_mid = nghost + (nx_phys + 1)/2
      do concurrent(k=1:nz, i=i_lo:i_mid)
         ip = isum - i
         if (ip == i) then
            v(i, j_fold, k) = 0.0_wp
         else
            v(i, j_fold, k) = -v(ip, j_fold, k)
         end if
      end do
   end subroutine fold_north_v_face_3d

   ! ================================================================
   ! corner-stagger (Bu, vorticity / PV diag): on-line self-conjugate row
   ! + north halos.  Symmetric storage (extent nx_total+1):
   !   c' = ni+2-c, j' = 2nj-j (halo) / fold row at j = nghost+nj.
   ! negate=.true. for true-vector corner components; negate=.false. for
   ! scalars (vorticity is a pseudoscalar — invariant here, so it copies).
   ! ================================================================

   pure subroutine fold_north_corner_2d(fld, nx_face, ny_face, &
                                        nx_phys, ny_phys, nghost, negate)
      !! 2D Bu-corner fold: north-halo fill + on-line projection.
      integer, intent(in) :: nx_face, ny_face, nx_phys, ny_phys, nghost
      real(wp), intent(inout) :: fld(nx_face, ny_face)
         !! Corner field, shape (nx_total+1, ny_total+1).
      logical, intent(in) :: negate
         !! .true. → negate (true-vector component); .false. → copy (scalar
         !! / pseudoscalar vorticity).

      integer :: i, j, fsum, jsum, j_fold, i_lo, i_mid, ip
      real(wp) :: sgn

      fsum = 2*nghost + nx_phys + 2
      jsum = 2*nghost + 2*ny_phys
      j_fold = nghost + ny_phys
      i_lo = nghost + 1
      sgn = merge(-1.0_wp, 1.0_wp, negate)

      ! Halo rows strictly beyond the fold row.
      do concurrent(j=j_fold + 1:ny_face, i=1:nx_face)
         fld(i, j) = sgn*fld(fsum - i, jsum - j)
      end do

      ! On-line projection at j = j_fold (west half from mirror).
      i_mid = nghost + (nx_phys + 2)/2
      do concurrent(i=i_lo:i_mid)
         ip = fsum - i
         if (negate .and. ip == i) then
            fld(i, j_fold) = 0.0_wp
         else
            fld(i, j_fold) = sgn*fld(ip, j_fold)
         end if
      end do
   end subroutine fold_north_corner_2d

end module rdb_ocean_fold
