!! Tests for the tripolar north-fold seam exchange (design Appendix A).
!! Hard-codes the ni=8 acceptance vectors (translated to physical+ghost
!! storage indexing).  Mirrors test_ocean_periodic.F90 conventions:
!! test-drive, small grids, explicit-shape helper calls, host asserts
!! after `acc update self` for the GPU-resident case.
!!
!! Index conventions (ng=nghost): storage index = ng + physical index.
!!   T:      i'=ni+1-i  -> storage isum=2ng+ni+1 ; j halo jsum=2ng+2nj+1
!!   u (Cu): f'=ni+2-f  -> storage fsum=2ng+ni+2 ; (sym storage, nx+1)
!!   v (Cv): i'=ni+1-i  ; fold row j_fold=ng+nj ; halo jsum=2ng+2nj
!!   corner: c'=ni+2-c  ; fold row j_fold=ng+nj
module test_ocean_fold
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_fold, only: fold_north_centre, fold_north_u_face, &
                             fold_north_v_face, fold_north_corner
   implicit none
   private

   public :: collect_ocean_fold_tests

   integer, parameter :: NGHOST = 3
   integer, parameter :: NI = 8, NJ = 6

contains

   subroutine collect_ocean_fold_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("fold_centre_acceptance", test_centre_acceptance), &
                  new_unittest("fold_vector_halo", test_vector_halo), &
                  new_unittest("fold_v_online_antisym", test_v_online), &
                  new_unittest("fold_involution", test_involution), &
                  new_unittest("fold_u_map_pairs", test_u_map_pairs), &
                  new_unittest("fold_centre_3d_levels", test_centre_3d_levels), &
                  new_unittest("fold_cyclic_corner", test_cyclic_corner), &
                  new_unittest("fold_corner_scalar_vs_vector", test_corner_modes), &
                  new_unittest("fold_centre_gpu", test_centre_gpu) &
                  ]
   end subroutine collect_ocean_fold_tests

   subroutine make_grid(grid)
      type(hgrid_t), intent(out) :: grid
      call grid%init(NI, NJ, NGHOST, 1.0_wp, 1.0_wp)
   end subroutine make_grid

   ! Periodic-x wrap of a 2D centre array (callers apply this BEFORE the fold).
   subroutine wrapx_centre(a, nxt, nyt)
      integer, intent(in) :: nxt, nyt
      real(wp), intent(inout) :: a(nxt, nyt)
      integer :: i, j
      do j = 1, nyt
         do i = 1, NGHOST
            a(i, j) = a(i + NI, j)
            a(NI + NGHOST + i, j) = a(NGHOST + i, j)
         end do
      end do
   end subroutine wrapx_centre

   ! Periodic-x wrap of a 2D x-face array (extent nxt+1).
   subroutine wrapx_face(a, nxf, nyt)
      integer, intent(in) :: nxf, nyt
      real(wp), intent(inout) :: a(nxf, nyt)
      integer :: i, j
      do j = 1, nyt
         do i = 1, NGHOST
            a(i, j) = a(i + NI, j)                 ! west ghost faces
         end do
         do i = NGHOST + NI + 2, nxf               ! east ghost faces
            a(i, j) = a(i - NI, j)
         end do
      end do
   end subroutine wrapx_face

   ! -----------------------------------------------------------------
   ! Test 1: T-centre acceptance — f(i,j)=phys_i gives halo = 9-i.
   ! -----------------------------------------------------------------
   subroutine test_centre_acceptance(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      integer :: nxt, nyt, i, j, ip, jlo
      real(wp), allocatable :: fld(:, :)

      checks: block
         call make_grid(grid)
         nxt = grid%nx_total; nyt = grid%ny_total
         allocate (fld(nxt, nyt), source=0.0_wp)
         ! Seed physical cells with their physical i index.
         do j = 1, nyt
            do ip = 1, NI
               fld(NGHOST + ip, j) = real(ip, wp)
            end do
         end do
         call wrapx_centre(fld, nxt, nyt)
         call fold_north_centre(fld, nxt, nyt, NI, NJ, NGHOST)

         ! First north halo row = ng+nj+1; Appendix vector: halo(phys i)=9-i.
         jlo = NGHOST + NJ + 1
         do ip = 1, NI
            call check(error, fld(NGHOST + ip, jlo) == real(9 - ip, wp), &
                       "T-fold halo != 9-i acceptance vector")
            if (allocated(error)) exit checks
         end do

         ! Constant scalar is unchanged (copy, no sign flip).
         fld = 7.5_wp
         call fold_north_centre(fld, nxt, nyt, NI, NJ, NGHOST)
         do j = jlo, nyt
            do i = 1, nxt
               call check(error, fld(i, j) == 7.5_wp, "T-fold: constant scalar changed")
               if (allocated(error)) exit checks
            end do
         end do
         deallocate (fld)
      end block checks
   end subroutine test_centre_acceptance

   ! -----------------------------------------------------------------
   ! Test 2: vector halo (u,v)=(1,1) -> (-1,-1) in vector halos.
   ! -----------------------------------------------------------------
   subroutine test_vector_halo(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      integer :: nxt, nyt, nxf, nyf, i, j, jlo
      real(wp), allocatable :: u(:, :), v(:, :, :)

      checks: block
         call make_grid(grid)
         nxt = grid%nx_total; nyt = grid%ny_total
         nxf = nxt + 1; nyf = nyt + 1

         allocate (u(nxf, nyt), source=1.0_wp)
         call fold_north_u_face(u, nxf, nyt, NI, NJ, NGHOST)
         jlo = NGHOST + NJ + 1
         do j = jlo, nyt
            do i = 1, nxf
               call check(error, u(i, j) == -1.0_wp, "u-fold halo != -1")
               if (allocated(error)) exit checks
            end do
         end do
         deallocate (u)

         ! v halo rows strictly beyond the fold row -> -1.
         allocate (v(nxt, nyf, 1), source=1.0_wp)
         call fold_north_v_face(v, nxt, nyf, 1, NI, NJ, NGHOST)
         do j = NGHOST + NJ + 1, nyf
            do i = 1, nxt
               call check(error, v(i, j, 1) == -1.0_wp, "v-fold halo != -1")
               if (allocated(error)) exit checks
            end do
         end do
         deallocate (v)
      end block checks
   end subroutine test_vector_halo

   ! -----------------------------------------------------------------
   ! Test 3: v on-line antisymmetry v(i,nj) + v(9-i,nj) = 0; fixed pt = 0.
   ! -----------------------------------------------------------------
   subroutine test_v_online(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      integer :: nxt, nyt, nyf, j, ip, jf
      real(wp), allocatable :: v(:, :, :)
      real(wp) :: s

      checks: block
         call make_grid(grid)
         nxt = grid%nx_total; nyt = grid%ny_total; nyf = nyt + 1
         allocate (v(nxt, nyf, 1), source=0.0_wp)
         do j = 1, nyf
            do ip = 1, NI
               v(NGHOST + ip, j, 1) = real(ip, wp)
            end do
         end do
         call wrapx_centre(v(:, :, 1), nxt, nyf)
         call fold_north_v_face(v, nxt, nyf, 1, NI, NJ, NGHOST)

         jf = NGHOST + NJ
         do ip = 1, NI
            s = v(NGHOST + ip, jf, 1) + v(NGHOST + (NI + 1 - ip), jf, 1)
            call check(error, abs(s) < 1.0e-13_wp, &
                       "v on-line not antisymmetric: v(i)+v(9-i) /= 0")
            if (allocated(error)) exit checks
         end do
         deallocate (v)
      end block checks
   end subroutine test_v_online

   ! -----------------------------------------------------------------
   ! Test 4: involution — fold twice == once (per stagger), byte-identical.
   ! -----------------------------------------------------------------
   subroutine test_involution(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      integer :: nxt, nyt, nxf, nyf, i, j, jlo
      real(wp), allocatable :: c(:, :), c2(:, :), u(:, :), u2(:, :)

      checks: block
         call make_grid(grid)
         nxt = grid%nx_total; nyt = grid%ny_total
         nxf = nxt + 1; nyf = nyt + 1

         ! Centre: deterministic pattern, periodic-wrapped.
         allocate (c(nxt, nyt), c2(nxt, nyt))
         do j = 1, nyt
            do i = 1, nxt
               c(i, j) = sin(real(i, wp))*cos(real(j, wp)) + real(i*j, wp)*0.01_wp
            end do
         end do
         call wrapx_centre(c, nxt, nyt)
         c2 = c
         call fold_north_centre(c, nxt, nyt, NI, NJ, NGHOST)
         ! Apply twice: re-wrap (composition rule) then fold.
         call fold_north_centre(c2, nxt, nyt, NI, NJ, NGHOST)
         call wrapx_centre(c2, nxt, nyt)
         call fold_north_centre(c2, nxt, nyt, NI, NJ, NGHOST)
         jlo = NGHOST + NJ + 1
         do j = jlo, nyt
            do i = 1, nxt
               call check(error, c(i, j) == c2(i, j), "centre fold not involutive on halo")
               if (allocated(error)) exit checks
            end do
         end do
         deallocate (c, c2)

         ! u-face involution (vector, sign flip): twice == once on halo.
         allocate (u(nxf, nyt), u2(nxf, nyt))
         do j = 1, nyt
            do i = 1, nxf
               u(i, j) = cos(real(i + j, wp))
            end do
         end do
         call wrapx_face(u, nxf, nyt)
         u2 = u
         call fold_north_u_face(u, nxf, nyt, NI, NJ, NGHOST)
         call fold_north_u_face(u2, nxf, nyt, NI, NJ, NGHOST)
         call wrapx_face(u2, nxf, nyt)
         call fold_north_u_face(u2, nxf, nyt, NI, NJ, NGHOST)
         do j = jlo, nyt
            do i = 1, nxf
               call check(error, u(i, j) == u2(i, j), "u fold not involutive on halo")
               if (allocated(error)) exit checks
            end do
         end do
         deallocate (u, u2)
      end block checks
   end subroutine test_involution

   ! -----------------------------------------------------------------
   ! Test 5: u-map pairs for OUR symmetric storage: f' = ni+2-f.
   !   Pairs (phys faces 1..9): 1<->9, 2<->8, 3<->7, 4<->6, 5<->5.
   !   Derivation: u_face_x(f) = WEST face of T-cell f; face f borders
   !   cells {f-1,f}; fold -> {ni+2-f, ni+1-f}; shared boundary = west
   !   face of cell ni+2-f.  Sym storage (nx+1) keeps the east-boundary
   !   face so the map is the palindromic ni+2-f, NOT the nonsym ni-i.
   ! -----------------------------------------------------------------
   subroutine test_u_map_pairs(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      integer :: nxt, nyt, nxf, f, fp, jlo, j_src
      real(wp), allocatable :: u(:, :)

      checks: block
         call make_grid(grid)
         nxt = grid%nx_total; nyt = grid%ny_total; nxf = nxt + 1
         allocate (u(nxf, nyt), source=0.0_wp)
         ! Seed each physical face f with value f, on the source row that the
         ! first halo row reflects from: jsum-(ng+nj+1) = ng+nj (the fold-
         ! adjacent owned row).  Seed all rows = f so the source is unambiguous.
         do f = 1, NI + 1
            u(NGHOST + f, :) = real(f, wp)
         end do
         call wrapx_face(u, nxf, nyt)
         call fold_north_u_face(u, nxf, nyt, NI, NJ, NGHOST)
         jlo = NGHOST + NJ + 1
         j_src = NGHOST + NJ            ! row the halo reflects from (value = f)
         ! halo(phys f) = -u_src(phys f') with f' = ni+2-f.
         do f = 1, NI + 1
            fp = NI + 2 - f
            call check(error, u(NGHOST + f, jlo) == -real(fp, wp), &
                       "u-map pair mismatch: halo(f) != -(ni+2-f)")
            if (allocated(error)) exit checks
         end do
         ! Explicit named pairs.
         call check(error, abs(u(NGHOST + 1, jlo)) == 9.0_wp, "u pair 1<->9")
         if (allocated(error)) exit checks
         call check(error, abs(u(NGHOST + 5, jlo)) == 5.0_wp, "u fixed pt 5<->5")
         deallocate (u)
      end block checks
   end subroutine test_u_map_pairs

   ! -----------------------------------------------------------------
   ! Test 6: 3D per-level independence — k-dependent field folds per level.
   ! -----------------------------------------------------------------
   subroutine test_centre_3d_levels(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      integer :: nxt, nyt, nz, i, j, k, jlo, isum, jsum
      real(wp), allocatable :: fld(:, :, :)

      checks: block
         call make_grid(grid)
         nxt = grid%nx_total; nyt = grid%ny_total; nz = 4
         allocate (fld(nxt, nyt, nz))
         do k = 1, nz
            do j = 1, nyt
               do i = 1, nxt
                  fld(i, j, k) = real(i, wp) + 100.0_wp*real(j, wp) + 10000.0_wp*real(k, wp)
               end do
            end do
         end do
         call fold_north_centre(fld, nxt, nyt, nz, NI, NJ, NGHOST)
         isum = 2*NGHOST + NI + 1
         jsum = 2*NGHOST + 2*NJ + 1
         jlo = NGHOST + NJ + 1
         do k = 1, nz
            do j = jlo, nyt
               do i = 1, nxt
                  call check(error, fld(i, j, k) == &
                             real(isum - i, wp) + 100.0_wp*real(jsum - j, wp) + 10000.0_wp*real(k, wp), &
                             "3D fold: level not independent / wrong source")
                  if (allocated(error)) exit checks
               end do
            end do
         end do
         deallocate (fld)
      end block checks
   end subroutine test_centre_3d_levels

   ! -----------------------------------------------------------------
   ! Test 7: fold∘cyclic at the two top corners (u-face).
   !   Seed each face with its STORAGE index, periodic-wrap, fold.  The
   !   halo at phys face 1 (storage ng+1) reflects from storage fsum-(ng+1)
   !   = ng+ni+1 (the east seam face), value = ng+ni+1, negated.
   ! -----------------------------------------------------------------
   subroutine test_cyclic_corner(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      integer :: nxt, nyt, nxf, i, j, jlo
      real(wp), allocatable :: u(:, :)

      checks: block
         call make_grid(grid)
         nxt = grid%nx_total; nyt = grid%ny_total; nxf = nxt + 1
         allocate (u(nxf, nyt))
         do j = 1, nyt
            do i = 1, nxf
               u(i, j) = real(i, wp)
            end do
         end do
         call wrapx_face(u, nxf, nyt)
         call fold_north_u_face(u, nxf, nyt, NI, NJ, NGHOST)
         jlo = NGHOST + NJ + 1
         ! West top corner: phys face 1 -> -(ng+ni+1) = -(3+8+1) = -12.
         call check(error, u(NGHOST + 1, jlo) == -real(NGHOST + NI + 1, wp), &
                    "fold-cyclic west corner wrong")
         if (allocated(error)) exit checks
         ! East top corner: phys face ni+1=9 -> -(ng+1) = -4.
         call check(error, u(NGHOST + NI + 1, jlo) == -real(NGHOST + 1, wp), &
                    "fold-cyclic east corner wrong")
         deallocate (u)
      end block checks
   end subroutine test_cyclic_corner

   ! -----------------------------------------------------------------
   ! Test 8: corner scalar (negate=.false. copies) vs vector (negate=.true.
   !   flips + zeroes self-fixed column).
   ! -----------------------------------------------------------------
   subroutine test_corner_modes(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      integer :: nxt, nyt, nxf, nyf, jf, c, cp
      real(wp), allocatable :: s(:, :), w(:, :)

      checks: block
         call make_grid(grid)
         nxt = grid%nx_total; nyt = grid%ny_total
         nxf = nxt + 1; nyf = nyt + 1
         jf = NGHOST + NJ

         ! Scalar (vorticity diag) — on-line copy, no sign flip.
         allocate (s(nxf, nyf), source=0.0_wp)
         do c = 1, NI + 1
            s(NGHOST + c, jf) = real(c, wp)
         end do
         call wrapx_face(s(:, jf:jf), nxf, 1)
         call fold_north_corner(s, nxf, nyf, NI, NJ, NGHOST, negate=.false.)
         ! c' = ni+2-c; scalar copies (sign +): s(c)=s(c') after projection.
         do c = 1, (NI + 2)/2
            cp = NI + 2 - c
            call check(error, s(NGHOST + c, jf) == s(NGHOST + cp, jf), &
                       "corner scalar: on-line copy not symmetric")
            if (allocated(error)) exit checks
         end do
         deallocate (s)

         ! Vector corner — negate + zero self-fixed column (c=5 for ni=8).
         allocate (w(nxf, nyf), source=0.0_wp)
         do c = 1, NI + 1
            w(NGHOST + c, jf) = real(c, wp)
         end do
         call wrapx_face(w(:, jf:jf), nxf, 1)
         call fold_north_corner(w, nxf, nyf, NI, NJ, NGHOST, negate=.true.)
         call check(error, w(NGHOST + 5, jf) == 0.0_wp, &
                    "corner vector self-fixed column (c=5) not zeroed")
         if (allocated(error)) exit checks
         ! Antisymmetry on-line: w(c) + w(ni+2-c) = 0.
         do c = 1, (NI + 2)/2
            cp = NI + 2 - c
            call check(error, abs(w(NGHOST + c, jf) + w(NGHOST + cp, jf)) < 1.0e-13_wp, &
                       "corner vector not antisymmetric on-line")
            if (allocated(error)) exit checks
         end do
         deallocate (w)
      end block checks
   end subroutine test_corner_modes

   ! -----------------------------------------------------------------
   ! Test 9: GPU-resident centre fold — map, fold on device, update self.
   ! -----------------------------------------------------------------
   subroutine test_centre_gpu(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      integer :: nxt, nyt, nz, i, j, k, ip, jlo
      real(wp), allocatable :: fld(:, :, :)

      checks: block
         call make_grid(grid)
         nxt = grid%nx_total; nyt = grid%ny_total; nz = 3
         allocate (fld(nxt, nyt, nz), source=0.0_wp)
         do k = 1, nz
            do j = 1, nyt
               do ip = 1, NI
                  fld(NGHOST + ip, j, k) = real(ip, wp) + 1000.0_wp*real(k, wp)
               end do
            end do
         end do
         call wrapx_centre(fld(:, :, 1), nxt, nyt)
         do k = 2, nz
            call wrapx_centre(fld(:, :, k), nxt, nyt)
         end do

         !$acc enter data copyin(fld)
         call fold_north_centre(fld, nxt, nyt, nz, NI, NJ, NGHOST)
         !$acc update self(fld)
         !$acc exit data delete(fld)

         jlo = NGHOST + NJ + 1
         do k = 1, nz
            do ip = 1, NI
               call check(error, fld(NGHOST + ip, jlo, k) == &
                          real(9 - ip, wp) + 1000.0_wp*real(k, wp), &
                          "GPU centre fold: wrong halo value")
               if (allocated(error)) exit checks
            end do
         end do
         deallocate (fld)
      end block checks
   end subroutine test_centre_gpu

end module test_ocean_fold
