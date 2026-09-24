!! Tripolar north-fold seam-exchange helpers for the ocean dyn-core.
module rdb_ocean_fold
   !! Single-rank discrete tripolar north-fold exchange (Murray 1996). Pure
   !! seam operators only (state orchestration lives in `rdb_ocean_fold_apply`).
   !! Free procedures, explicit-shape dummies, j-outer / i-inner `do concurrent`.
   !!
   !! ## Roundabout staggering (the load-bearing input to every map below)
   !!
   !! Continuous grid coordinates: T-cell (i,j), physical i ∈ 1..ni,
   !! j ∈ 1..nj, occupies [i-1,i]×[j-1,j]; storage index = nghost + physical.
   !!   * T  `h(i,j)`         centre      (i-1/2, j-1/2)
   !!   * Cu `u(i,j)`         WEST face   (i-1,   j-1/2)   extent nx+1
   !!   * Cv `v(i,j)`         SOUTH face  (i-1/2, j-1)     extent ny+1
   !!   * Bu `q(i,j)`         SW corner   (i-1,   j-1)     extent (nx+1,ny+1)
   !! (`rdb_multilayer_state` `u_face_x_layer`/`v_face_y_layer` docstrings;
   !! `coriolis_adv` "zeta_corner(i,j) sits at (i-1/2,j-1/2)" relative to
   !! T(i,j); `metrics%wet_q` "SW corner of T-cell (i,j)".)
   !!
   !! ## The fold
   !!
   !! The fold line is y = nj, the NORTH edge of T-row nj.  A point (x, y)
   !! north of it is the point (ni - x, 2nj - y) (i-periodic, period ni),
   !! reached through a 180° rotation of the local frame: both unit vectors
   !! reverse, so a true-vector component (u, v, a face flux) NEGATES and a
   !! scalar — and the pseudoscalar vorticity / PV (rotation preserves
   !! orientation) — COPIES.
   !!
   !! Solving x' = ni - x, y' = 2nj - y for each stagger's storage index:
   !! | Stagger | x(i)  | y(j)  | i-map (storage)       | j-map (storage)       | fold-line row |
   !! |---------|-------|-------|-----------------------|-----------------------|---------------|
   !! | T       | i-½   | j-½   | i' = 2ng+ni+1 - i     | j' = 2ng+2nj+1 - j    | none          |
   !! | u (Cu)  | i-1   | j-½   | i' = 2ng+ni+2 - i     | j' = 2ng+2nj+1 - j    | none          |
   !! | v (Cv)  | i-½   | j-1   | i' = 2ng+ni+1 - i     | j' = 2ng+2nj+2 - j    | ng+nj+1       |
   !! | corner  | i-1   | j-1   | i' = 2ng+ni+2 - i     | j' = 2ng+2nj+2 - j    | ng+nj+1       |
   !!
   !! T and u points never lie on y = nj, so their exchange is a pure
   !! halo fill of rows j > ng+nj.  v and corners have a row ON the fold
   !! line: storage row `ng+nj+1` — the south face of the first ghost row,
   !! i.e. the NORTH face of the last physical T-row.  (The pre-2026-09
   !! code used `jsum = 2ng+2nj`, `j_fold = ng+nj`: MOM6's NORTH-face /
   !! NE-corner rule applied to roundabout's SOUTH-face / SW-corner
   !! storage.  It antisymmetrised the south face of the last T-row — an
   !! ordinary interior face — and filled the true fold-line face from
   !! `-v(i', ng+nj-1)`, so every cell of the last row took an unrelated
   !! flux through its north face: the global-tripolar B1 mass leak.)
   !!
   !! Cross-check against MOM6 (symmetric memory, `pass_vector` /
   !! FMS `mpp_update_domains` with a folded north edge, CGRID_NE):
   !! MOM6 `v(i,J)` is the north face of cell j = roundabout `v(i,J+1)`,
   !! MOM6 `q(I,J)` NE corner = roundabout `q(I+1,J+1)`, MOM6 `u(I,j)` east
   !! face = roundabout `u(I+1,j)`.  MOM6's fold maps v(i, nj+d) ←
   !! -v(ni+1-i, nj-d), q(I, nj+d) ← q(ni-I, nj-d), u(I, nj+d) ←
   !! -u(ni-I, nj+1-d); shifting by the index offsets gives exactly the
   !! table above (MOM6's fold-line row J = nj ↔ roundabout row nj+1).
   !!
   !! ## The duplicated-DOF fold-line row (v and corners)
   !!
   !! On row ng+nj+1 the storage slots i and i' (v: i' = 2ng+ni+1-i;
   !! corner: i' = 2ng+ni+2-i) are the SAME physical face / vertex seen from
   !! the two sides of the fold, with opposite orientation.  So a single
   !! DOF is stored twice and must satisfy v(i) = -v(i') (a true normal
   !! velocity / normal flux through one edge, which leaves cell (i,nj)
   !! northward and ENTERS cell (i',nj) from its north) and q(i) = q(i')
   !! for scalars/vorticity.  The dynamics updates both slots
   !! independently; the projection overwrites the WEST half from the
   !! (negated) east mirror so the row is exactly (anti)symmetric.  This is
   !! what makes the cross-fold mass flux telescope: Σ_i F(i, ng+nj+1)
   !! pairs off to zero.  Self-conjugate slots (i = i'):
   !!   * v: only when ni is odd (column (ni+1)/2) → 0 (a normal velocity
   !!     equal to minus itself).  ni even (every real tripolar grid) has
   !!     none.
   !!   * corner: c = ni/2+1 and c = 1 ≡ ni+1 (periodic images) — the two
   !!     bipoles.  Vector corner components → 0; scalars copy.
   !! Rows above the fold line (j > ng+nj+1) are pure mirrored images.
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
   ! v-stagger (Cv, SOUTH y-face): TWO operations.
   !   (1) halo-fill rows j > j_fold (reflected + negated):
   !         i' = ni+1-i, j' = 2nj+2-j  (storage j' = 2*nghost+2*nj+2 - j)
   !   (2) ON-LINE projection at j = j_fold (= nghost+nj+1, the north face
   !       of the last physical T-row):
   !         v(i,j_fold) = -v(i', j_fold), i' = ni+1-i
   !       overwrite the WEST half (and its periodic ghost images) from the
   !       negated east-mirror; self-fixed column (odd ni) → 0.
   ! ================================================================

   pure subroutine fold_north_v_face_2d(v, nx_total, ny_face, &
                                        nx_phys, ny_phys, nghost)
      !! 2D y-face (Cv) fold: north-halo fill (negated) + on-line
      !! antisymmetric projection at the fold row.  Used for the
      !! barotropic `bt_vbt` field in the BT fast loop.
      integer, intent(in) :: nx_total, ny_face, nx_phys, ny_phys, nghost
      real(wp), intent(inout) :: v(nx_total, ny_face)
         !! y-face 2D field, shape (nx_total, ny_total+1).

      integer :: i, j, isum, jsum, j_fold, i_lo, p, pm

      isum = 2*nghost + nx_phys + 1
      jsum = 2*nghost + 2*ny_phys + 2   ! v is SOUTH-face: y = j-1
      j_fold = nghost + ny_phys + 1     ! north face of the last T-row
      i_lo = nghost + 1

      ! (1) Halo rows strictly beyond the fold row.
      do concurrent(j=j_fold + 1:ny_face, i=1:nx_total)
         v(i, j) = -v(isum - i, jsum - j)
      end do

      ! (2) On-line antisymmetric projection at j = j_fold, over EVERY
      !     storage column (periodic ghosts included, so the row stays
      !     periodic-consistent without a second wrap): a column whose
      !     physical index p lies in the west half takes -v of the east
      !     mirror p' = ni+1-p; the east half is the (read-only) source.
      do concurrent(i=1:nx_total) local(p, pm)
         p = modulo(i - i_lo, nx_phys) + 1
         pm = nx_phys + 1 - p
         if (p < pm) then
            v(i, j_fold) = -v(nghost + pm, j_fold)
         else if (p == pm) then
            v(i, j_fold) = 0.0_wp
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

      integer :: i, j, k, isum, jsum, j_fold, i_lo, p, pm

      isum = 2*nghost + nx_phys + 1
      jsum = 2*nghost + 2*ny_phys + 2    ! v is SOUTH-face: y = j-1
      j_fold = nghost + ny_phys + 1      ! the self-conjugate fold row
      i_lo = nghost + 1                  ! first physical column

      ! (1) Halo rows strictly beyond the fold row.
      do concurrent(k=1:nz, j=j_fold + 1:ny_face, i=1:nx_total)
         v(i, j, k) = -v(isum - i, jsum - j, k)
      end do

      ! (2) On-line antisymmetric projection at j = j_fold (see the 2D
      !     twin): every storage column whose physical index p is in the
      !     west half takes -v of its east mirror p' = ni+1-p; the
      !     self-conjugate column (odd ni only) is zeroed; the east half is
      !     the read-only source, so the kernel is race-free.
      do concurrent(k=1:nz, i=1:nx_total) local(p, pm)
         p = modulo(i - i_lo, nx_phys) + 1
         pm = nx_phys + 1 - p
         if (p < pm) then
            v(i, j_fold, k) = -v(nghost + pm, j_fold, k)
         else if (p == pm) then
            v(i, j_fold, k) = 0.0_wp
         end if
      end do
   end subroutine fold_north_v_face_3d

   ! ================================================================
   ! corner-stagger (Bu, SW corner; vorticity / PV / f): on-line
   ! self-conjugate row + north halos.  Symmetric storage (extent nx_total+1):
   !   c' = ni+2-c, j' = 2nj+2-j (halo) / fold row at j = nghost+nj+1.
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

      integer :: i, j, fsum, jsum, j_fold, i_lo, p, pm
      real(wp) :: sgn

      fsum = 2*nghost + nx_phys + 2
      jsum = 2*nghost + 2*ny_phys + 2   ! corner is SW: y = j-1
      j_fold = nghost + ny_phys + 1
      i_lo = nghost + 1
      sgn = merge(-1.0_wp, 1.0_wp, negate)

      ! Halo rows strictly beyond the fold row.
      do concurrent(j=j_fold + 1:ny_face, i=1:nx_face)
         fld(i, j) = sgn*fld(fsum - i, jsum - j)
      end do

      ! On-line projection at j = j_fold over EVERY storage column
      ! (periodic ghosts included).  Physical corner index p ∈ 1..ni
      ! (c = ni+1 is the periodic image of c = 1); mirror p' = ni+2-p
      ! modulo ni.  The two self-conjugate corners p = 1 and p = ni/2+1 are
      ! the bipoles: a vector component is zeroed there, a scalar is left
      ! as is.  West-half columns take sgn*mirror; the east half is the
      ! read-only source.
      do concurrent(i=1:nx_face) local(p, pm)
         p = modulo(i - i_lo, nx_phys) + 1
         pm = modulo(nx_phys + 1 - p, nx_phys) + 1
         if (p < pm) then
            fld(i, j_fold) = sgn*fld(nghost + pm, j_fold)
         else if (p == pm .and. negate) then
            fld(i, j_fold) = 0.0_wp
         end if
      end do
   end subroutine fold_north_corner_2d

end module rdb_ocean_fold
