!! MPI integration test for the ocean staggered halo exchange.
!!
!! Branches on the launched rank count:
!!   nprocs == 4 : 2x2 non-periodic grid.  GLOBAL-INDEX-SIGNATURE fills
!!     (value = gi*10000 + gj*100 + k) so every ghost's expected value is
!!     computable independently of which rank supplied it — any index-map
!!     error in a pack/unpack fails loudly.  Covers centre_2d, centre_3d,
!!     face_x_2d/3d, face_y_2d/3d, bt_group_2d; asserts side ghost bands,
!!     CORNER ghost blocks (the D2 two-pass invariant), and the duplicated
!!     seam faces.
!!   nprocs == 2 : px=2, periodic_x.  Strongest check: the same signature
!!     field is built once as a full-domain single-rank array wrapped with
!!     the ocean_periodic_wrap_* kernels (the 1-rank reference), and once
!!     as rank-local arrays exchanged via ocean_halo_*; the rank's FULL
!!     local window (ghosts included) must match the reference window
!!     BIT-EXACTLY.  Uses a periodic-consistent fill (signature of
!!     mod(gi-1, NXG)+1) so the two identified copies of the wrap face
!!     hold equal values regardless of ownership direction.
!!   nprocs == 1 : D0 gate.  Periodic: each ocean_halo_* exchange must be
!!     bit-identical to calling the matching ocean_periodic_wrap_* kernel
!!     directly.  Non-periodic: every exchange must leave the array
!!     bit-unchanged (no-op).
!!   other       : skip + pass.
!!
!! All exchanges here pass `device_resident=.false.` (host mode): the test
!! arrays are plain host locals with no `!$acc enter data` mapping, so on a
!! GPU+MPI build the default device path would abort on `present(fld)`.
!! The device-resident path is covered end-to-end on the GPU+MPI build by
!! test_ocean_dyn_mpi's per-step exchanges (state device-mapped there).
#ifdef RDB_ENABLE_MPI
program test_halo_ocean_mpi
   use rdb_constants, only: wp
   use rdb_decomp, only: decomp_t, decomp_init
   use rdb_ocean_halo, only: ocean_halo_init, ocean_halo_destroy, &
                             ocean_halo_reserve, &
                             ocean_halo_centre, &
                             ocean_halo_face_x, &
                             ocean_halo_face_y, &
                             ocean_halo_bt_group_2d, &
                             ocean_halo_centre_2d_wide, &
                             ocean_halo_face_x_2d_wide, &
                             ocean_halo_face_y_2d_wide, &
                             ocean_halo_bt_group_2d_wide
   use rdb_ocean_periodic, only: ocean_periodic_wrap_centre_2d, &
                                 ocean_periodic_wrap_face_x_2d, &
                                 ocean_periodic_wrap_face_y_2d
   use rdb_comm_env, only: comm_env_init, comm_env_setup_roles, &
                           comm_env_finalize, comm_env_rank, comm_env_size, &
                           comm_env_compute_comm
   use pic_mpi_lib, only: comm_t, allreduce, MPI_SUM
   implicit none

   integer, parameter :: NXG = 12
      !! Global physical cells in x
   integer, parameter :: NYG = 10
      !! Global physical cells in y
   integer, parameter :: NG = 3
      !! Ghost width
   integer, parameter :: NZ3 = 3
      !! Layers for the 3D legs
   real(wp), parameter :: SENTINEL = -777.0_wp
      !! Ghost-init value; must survive wherever no neighbour exists
   ! Wide-halo test parameters.
   ! np2 (px=2): nxl=6, nyl=10.  ng_wide=NG+2=5 <= nxl=6 (OK).
   ! np4 (2x2): nxl=6, nyl=5.   ng_wide=NG+1=4 <= nxl=6, nyl=5 (OK).
   integer, parameter :: NG_WIDE_NP2 = NG + 2
      !! Wide ghost width for np2 legs (=5; nxl=6 >= 5)
   integer, parameter :: NG_WIDE_NP4 = NG + 1
      !! Wide ghost width for np4 legs (=4; nyl=5 >= 4)

   integer :: rank, nprocs, n_fail, total_fail
   type(decomp_t) :: decomp
   type(comm_t) :: comm
   integer :: nxl, nyl, nxt, nyt
   logical :: fold_x = .false.
      !! Periodic-consistent signature: fold gi by mod(gi-1, NXG)+1

   call comm_env_init()
   call comm_env_setup_roles(.false.)
   rank = comm_env_rank()
   nprocs = comm_env_size()
   n_fail = 0

   select case (nprocs)
   case (4)
      call run_4rank()
      call run_4rank_wide()
   case (2)
      call run_2rank_periodic()
      call run_2rank_periodic_wide()
   case (1)
      call run_1rank()
      call run_1rank_wide()
   case default
      if (rank == 0) write (*, *) "SKIP: test supports 1/2/4 ranks, got", nprocs
   end select

   comm = comm_env_compute_comm()
   call comm%barrier()
   if (n_fail > 0) then
      write (*, *) "Rank", rank, ":", n_fail, "check(s) FAILED (nprocs=", nprocs, ")"
   else
      write (*, *) "Rank", rank, ": all checks PASSED (nprocs=", nprocs, ")"
   end if
   total_fail = n_fail
   call allreduce(comm, total_fail, MPI_SUM)
   call comm_env_finalize()
   if (total_fail > 0) error stop 1

contains

   ! =================================================================
   ! Signature + comparison helpers
   ! =================================================================

   pure function sig(gi, gj, k) result(v)
      !! Global-index signature.  When fold_x is set (periodic-x legs)
      !! gi is folded onto 1..NXG so the two identified copies of a
      !! wrap face/cell carry equal values.
      integer, intent(in) :: gi, gj, k
      real(wp) :: v
      integer :: gii
      gii = gi
      if (fold_x) gii = mod(gi - 1 + 4*NXG, NXG) + 1
      v = real(gii, wp)*10000.0_wp + real(gj, wp)*100.0_wp + real(k, wp)
   end function sig

   pure function zone(idx, n_phys) result(z)
      !! -1 = low ghost band, 0 = physical span, +1 = high ghost band.
      !! For face-staggered axes pass n_phys = n_local+1 (physical faces
      !! span NG+1 .. NG+n_local+1, including both seam copies).
      integer, intent(in) :: idx, n_phys
      integer :: z
      if (idx <= NG) then
         z = -1
      else if (idx > NG + n_phys) then
         z = 1
      else
         z = 0
      end if
   end function zone

   pure function dir_ok_x(z) result(ok)
      !! True if the x ghost band z is reachable (neighbour exists).
      integer, intent(in) :: z
      logical :: ok
      select case (z)
      case (-1)
         ok = .not. decomp%has_west
      case (1)
         ok = .not. decomp%has_east
      case default
         ok = .true.
      end select
   end function dir_ok_x

   pure function dir_ok_y(z) result(ok)
      !! True if the y ghost band z is reachable (neighbour exists).
      integer, intent(in) :: z
      logical :: ok
      select case (z)
      case (-1)
         ok = .not. decomp%has_south
      case (1)
         ok = .not. decomp%has_north
      case default
         ok = .true.
      end select
   end function dir_ok_y

   subroutine expect_eq(got, expect, label, i, j, k)
      !! Bit-exact comparison (all values are exactly representable).
      real(wp), intent(in) :: got, expect
      character(len=*), intent(in) :: label
      integer, intent(in) :: i, j, k
      if (got /= expect) then
         n_fail = n_fail + 1
         if (n_fail <= 8) write (*, *) "FAIL ", label, " rank=", rank, &
            " i=", i, " j=", j, " k=", k, " expect=", expect, " got=", got
      end if
   end subroutine expect_eq

   ! =================================================================
   ! Fills (physical cells/faces only; everything else = SENTINEL)
   ! =================================================================

   subroutine fill_centre_2d(fld, k)
      real(wp), intent(inout) :: fld(:, :)
      integer, intent(in) :: k
      integer :: i, j, gi, gj
      fld = SENTINEL
      do j = NG + 1, NG + nyl
         do i = NG + 1, NG + nxl
            gi = i - NG + decomp%i_start - 1
            gj = j - NG + decomp%j_start - 1
            fld(i, j) = sig(gi, gj, k)
         end do
      end do
   end subroutine fill_centre_2d

   subroutine fill_face_x_2d(fld, k)
      !! Physical x-faces span i = NG+1 .. NG+nxl+1 (both seam copies
      !! filled with the same global-face signature).
      real(wp), intent(inout) :: fld(:, :)
      integer, intent(in) :: k
      integer :: i, j, gf, gj
      fld = SENTINEL
      do j = NG + 1, NG + nyl
         do i = NG + 1, NG + nxl + 1
            gf = i - NG - 1 + decomp%i_start
            gj = j - NG + decomp%j_start - 1
            fld(i, j) = sig(gf, gj, k)
         end do
      end do
   end subroutine fill_face_x_2d

   subroutine fill_face_y_2d(fld, k)
      real(wp), intent(inout) :: fld(:, :)
      integer, intent(in) :: k
      integer :: i, j, gi, gfj
      fld = SENTINEL
      do j = NG + 1, NG + nyl + 1
         do i = NG + 1, NG + nxl
            gi = i - NG + decomp%i_start - 1
            gfj = j - NG - 1 + decomp%j_start
            fld(i, j) = sig(gi, gfj, k)
         end do
      end do
   end subroutine fill_face_y_2d

   ! =================================================================
   ! Checkers (non-periodic legs): full-window sweep.  A position is
   ! REACHABLE iff each ghost band it sits in has a live neighbour;
   ! reachable => global signature, unreachable => SENTINEL untouched.
   ! Corner blocks need BOTH neighbours (D2 two-pass transitivity).
   ! Duplicated seam faces are physical positions in BOTH owners'
   ! windows, so the sweep asserts both copies == the owner signature.
   ! =================================================================

   subroutine check_centre_2d(fld, k, label)
      real(wp), intent(in) :: fld(:, :)
      integer, intent(in) :: k
      character(len=*), intent(in) :: label
      integer :: i, j, gi, gj
      real(wp) :: expect
      do j = 1, nyt
         do i = 1, nxt
            gi = i - NG + decomp%i_start - 1
            gj = j - NG + decomp%j_start - 1
            if (dir_ok_x(zone(i, nxl)) .and. dir_ok_y(zone(j, nyl))) then
               expect = sig(gi, gj, k)
            else
               expect = SENTINEL
            end if
            call expect_eq(fld(i, j), expect, label, i, j, k)
         end do
      end do
   end subroutine check_centre_2d

   subroutine check_face_x_2d(fld, k, label)
      real(wp), intent(in) :: fld(:, :)
      integer, intent(in) :: k
      character(len=*), intent(in) :: label
      integer :: i, j, gf, gj
      real(wp) :: expect
      do j = 1, nyt
         do i = 1, nxt + 1
            gf = i - NG - 1 + decomp%i_start
            gj = j - NG + decomp%j_start - 1
            if (dir_ok_x(zone(i, nxl + 1)) .and. dir_ok_y(zone(j, nyl))) then
               expect = sig(gf, gj, k)
            else
               expect = SENTINEL
            end if
            call expect_eq(fld(i, j), expect, label, i, j, k)
         end do
      end do
   end subroutine check_face_x_2d

   subroutine check_face_y_2d(fld, k, label)
      real(wp), intent(in) :: fld(:, :)
      integer, intent(in) :: k
      character(len=*), intent(in) :: label
      integer :: i, j, gi, gfj
      real(wp) :: expect
      do j = 1, nyt + 1
         do i = 1, nxt
            gi = i - NG + decomp%i_start - 1
            gfj = j - NG - 1 + decomp%j_start
            if (dir_ok_x(zone(i, nxl)) .and. dir_ok_y(zone(j, nyl + 1))) then
               expect = sig(gi, gfj, k)
            else
               expect = SENTINEL
            end if
            call expect_eq(fld(i, j), expect, label, i, j, k)
         end do
      end do
   end subroutine check_face_y_2d

   ! =================================================================
   ! nprocs == 4 : 2x2 non-periodic, signature fills
   ! =================================================================

   subroutine run_4rank()
      real(wp), allocatable :: eta(:, :), ubt(:, :), vbt(:, :)
      real(wp), allocatable :: c3(:, :, :), fx3(:, :, :), fy3(:, :, :)
      integer :: k

      call decomp_init(decomp, NXG, NYG, 2, 2, rank)
      nxl = decomp%nx_local
      nyl = decomp%ny_local
      nxt = nxl + 2*NG
      nyt = nyl + 2*NG
      fold_x = .false.
      call ocean_halo_init(decomp, NG, .false., .false.)
      call ocean_halo_reserve(NZ3, 0)

      allocate (eta(nxt, nyt), ubt(nxt + 1, nyt), vbt(nxt, nyt + 1))
      allocate (c3(nxt, nyt, NZ3), fx3(nxt + 1, nyt, NZ3), fy3(nxt, nyt + 1, NZ3))

      ! --- centre 2D ---
      call fill_centre_2d(eta, 0)
      call ocean_halo_centre(eta, device_resident=.false.)
      call check_centre_2d(eta, 0, "centre_2d")

      ! --- centre 3D ---
      do k = 1, NZ3
         call fill_centre_2d(c3(:, :, k), k)
      end do
      call ocean_halo_centre(c3, NZ3, device_resident=.false.)
      do k = 1, NZ3
         call check_centre_2d(c3(:, :, k), k, "centre_3d")
      end do

      ! --- face_x 2D ---
      call fill_face_x_2d(ubt, 0)
      call ocean_halo_face_x(ubt, device_resident=.false.)
      call check_face_x_2d(ubt, 0, "face_x_2d")

      ! --- face_x 3D ---
      do k = 1, NZ3
         call fill_face_x_2d(fx3(:, :, k), k)
      end do
      call ocean_halo_face_x(fx3, NZ3, device_resident=.false.)
      do k = 1, NZ3
         call check_face_x_2d(fx3(:, :, k), k, "face_x_3d")
      end do

      ! --- face_y 2D ---
      call fill_face_y_2d(vbt, 0)
      call ocean_halo_face_y(vbt, device_resident=.false.)
      call check_face_y_2d(vbt, 0, "face_y_2d")

      ! --- face_y 3D ---
      do k = 1, NZ3
         call fill_face_y_2d(fy3(:, :, k), k)
      end do
      call ocean_halo_face_y(fy3, NZ3, device_resident=.false.)
      do k = 1, NZ3
         call check_face_y_2d(fy3(:, :, k), k, "face_y_3d")
      end do

      ! --- bt_group (all three in one call) ---
      call fill_centre_2d(eta, 0)
      call fill_face_x_2d(ubt, 0)
      call fill_face_y_2d(vbt, 0)
      call ocean_halo_bt_group_2d(eta, ubt, vbt, device_resident=.false.)
      call check_centre_2d(eta, 0, "bt_group_eta")
      call check_face_x_2d(ubt, 0, "bt_group_ubt")
      call check_face_y_2d(vbt, 0, "bt_group_vbt")

      call ocean_halo_destroy()
      deallocate (eta, ubt, vbt, c3, fx3, fy3)
   end subroutine run_4rank

   ! =================================================================
   ! nprocs == 2 : px=2, periodic-x — bit-exact against the wrapped
   ! single-rank reference field (periodic-consistent signature fill).
   ! =================================================================

   subroutine run_2rank_periodic()
      real(wp), allocatable :: fld(:, :), ref(:, :)
      integer :: i, j, off

      call decomp_init(decomp, NXG, NYG, 2, 1, rank)
      nxl = decomp%nx_local
      nyl = decomp%ny_local
      nxt = nxl + 2*NG
      nyt = nyl + 2*NG
      fold_x = .true.
      call ocean_halo_init(decomp, NG, .true., .false.)
      call ocean_halo_reserve(NZ3, 0)
      off = decomp%i_start - 1   ! local i  ->  reference index i + off

      ! --- centre ---
      allocate (ref(NXG + 2*NG, NYG + 2*NG))
      ref = SENTINEL
      do j = NG + 1, NG + NYG
         do i = NG + 1, NG + NXG
            ref(i, j) = sig(i - NG, j - NG, 0)
         end do
      end do
      call ocean_periodic_wrap_centre_2d(ref, NXG + 2*NG, NYG + 2*NG, &
                                         NXG, NYG, NG, .true., .false.)
      allocate (fld(nxt, nyt))
      call fill_centre_2d(fld, 0)
      call ocean_halo_centre(fld, device_resident=.false.)
      do j = 1, nyt
         do i = 1, nxt
            call expect_eq(fld(i, j), ref(i + off, j), "p2_centre", i, j, 0)
         end do
      end do
      deallocate (fld, ref)

      ! --- face_x (the wrap-seam duplicated face lives on this axis) ---
      allocate (ref(NXG + 2*NG + 1, NYG + 2*NG))
      ref = SENTINEL
      do j = NG + 1, NG + NYG
         do i = NG + 1, NG + NXG + 1
            ref(i, j) = sig(i - NG, j - NG, 0)   ! sig folds NXG+1 -> 1
         end do
      end do
      call ocean_periodic_wrap_face_x_2d(ref, NXG + 2*NG + 1, NYG + 2*NG, &
                                         NXG, NYG, NG, .true., .false.)
      allocate (fld(nxt + 1, nyt))
      call fill_face_x_2d(fld, 0)
      call ocean_halo_face_x(fld, device_resident=.false.)
      do j = 1, nyt
         do i = 1, nxt + 1
            call expect_eq(fld(i, j), ref(i + off, j), "p2_face_x", i, j, 0)
         end do
      end do
      deallocate (fld, ref)

      ! --- face_y (y faces, x-periodic exchange) ---
      allocate (ref(NXG + 2*NG, NYG + 2*NG + 1))
      ref = SENTINEL
      do j = NG + 1, NG + NYG + 1
         do i = NG + 1, NG + NXG
            ref(i, j) = sig(i - NG, j - NG, 0)
         end do
      end do
      call ocean_periodic_wrap_face_y_2d(ref, NXG + 2*NG, NYG + 2*NG + 1, &
                                         NXG, NYG, NG, .true., .false.)
      allocate (fld(nxt, nyt + 1))
      call fill_face_y_2d(fld, 0)
      call ocean_halo_face_y(fld, device_resident=.false.)
      do j = 1, nyt + 1
         do i = 1, nxt
            call expect_eq(fld(i, j), ref(i + off, j), "p2_face_y", i, j, 0)
         end do
      end do
      deallocate (fld, ref)

      call ocean_halo_destroy()
   end subroutine run_2rank_periodic

   ! =================================================================
   ! nprocs == 1 : D0 gate — halo == direct wrap kernel (periodic),
   ! then halo == no-op (non-periodic)
   ! =================================================================

   subroutine run_1rank()
      real(wp), allocatable :: a(:, :), b(:, :)
      real(wp), allocatable :: ea(:, :), ua(:, :), va(:, :)
      real(wp), allocatable :: eb(:, :), ub(:, :), vb(:, :)
      integer :: i, j

      call decomp_init(decomp, NXG, NYG, 1, 1, rank)
      nxl = decomp%nx_local
      nyl = decomp%ny_local
      nxt = nxl + 2*NG
      nyt = nyl + 2*NG
      fold_x = .false.

      ! ---- Periodic both axes: halo must equal the wrap kernel ----
      call ocean_halo_init(decomp, NG, .true., .true.)
      call ocean_halo_reserve(1, 0)

      ! centre
      allocate (a(nxt, nyt), b(nxt, nyt))
      call fill_centre_2d(a, 0)
      b = a
      call ocean_halo_centre(a, device_resident=.false.)
      call ocean_periodic_wrap_centre_2d(b, nxt, nyt, nxl, nyl, NG, .true., .true.)
      do j = 1, nyt
         do i = 1, nxt
            call expect_eq(a(i, j), b(i, j), "p1_centre", i, j, 0)
         end do
      end do
      deallocate (a, b)

      ! face_x
      allocate (a(nxt + 1, nyt), b(nxt + 1, nyt))
      call fill_face_x_2d(a, 0)
      b = a
      call ocean_halo_face_x(a, device_resident=.false.)
      call ocean_periodic_wrap_face_x_2d(b, nxt + 1, nyt, nxl, nyl, NG, .true., .true.)
      do j = 1, nyt
         do i = 1, nxt + 1
            call expect_eq(a(i, j), b(i, j), "p1_face_x", i, j, 0)
         end do
      end do
      deallocate (a, b)

      ! face_y
      allocate (a(nxt, nyt + 1), b(nxt, nyt + 1))
      call fill_face_y_2d(a, 0)
      b = a
      call ocean_halo_face_y(a, device_resident=.false.)
      call ocean_periodic_wrap_face_y_2d(b, nxt, nyt + 1, nxl, nyl, NG, .true., .true.)
      do j = 1, nyt + 1
         do i = 1, nxt
            call expect_eq(a(i, j), b(i, j), "p1_face_y", i, j, 0)
         end do
      end do
      deallocate (a, b)

      ! bt_group: exchange vs direct wraps of all three arrays
      allocate (ea(nxt, nyt), ua(nxt + 1, nyt), va(nxt, nyt + 1))
      allocate (eb(nxt, nyt), ub(nxt + 1, nyt), vb(nxt, nyt + 1))
      call fill_centre_2d(ea, 0)
      call fill_face_x_2d(ua, 0)
      call fill_face_y_2d(va, 0)
      eb = ea
      ub = ua
      vb = va
      call ocean_halo_bt_group_2d(ea, ua, va, device_resident=.false.)
      call ocean_periodic_wrap_centre_2d(eb, nxt, nyt, nxl, nyl, NG, .true., .true.)
      call ocean_periodic_wrap_face_x_2d(ub, nxt + 1, nyt, nxl, nyl, NG, .true., .true.)
      call ocean_periodic_wrap_face_y_2d(vb, nxt, nyt + 1, nxl, nyl, NG, .true., .true.)
      do j = 1, nyt
         do i = 1, nxt
            call expect_eq(ea(i, j), eb(i, j), "p1_bt_eta", i, j, 0)
         end do
      end do
      do j = 1, nyt
         do i = 1, nxt + 1
            call expect_eq(ua(i, j), ub(i, j), "p1_bt_ubt", i, j, 0)
         end do
      end do
      do j = 1, nyt + 1
         do i = 1, nxt
            call expect_eq(va(i, j), vb(i, j), "p1_bt_vbt", i, j, 0)
         end do
      end do
      deallocate (ea, ua, va, eb, ub, vb)

      call ocean_halo_destroy()

      ! ---- Non-periodic: every exchange is a bit-exact no-op ----
      call ocean_halo_init(decomp, NG, .false., .false.)
      call ocean_halo_reserve(1, 0)

      allocate (a(nxt, nyt), b(nxt, nyt))
      call fill_centre_2d(a, 0)
      b = a
      call ocean_halo_centre(a, device_resident=.false.)
      do j = 1, nyt
         do i = 1, nxt
            call expect_eq(a(i, j), b(i, j), "np1_centre", i, j, 0)
         end do
      end do
      deallocate (a, b)

      allocate (a(nxt + 1, nyt), b(nxt + 1, nyt))
      call fill_face_x_2d(a, 0)
      b = a
      call ocean_halo_face_x(a, device_resident=.false.)
      do j = 1, nyt
         do i = 1, nxt + 1
            call expect_eq(a(i, j), b(i, j), "np1_face_x", i, j, 0)
         end do
      end do
      deallocate (a, b)

      allocate (a(nxt, nyt + 1), b(nxt, nyt + 1))
      call fill_face_y_2d(a, 0)
      b = a
      call ocean_halo_face_y(a, device_resident=.false.)
      do j = 1, nyt + 1
         do i = 1, nxt
            call expect_eq(a(i, j), b(i, j), "np1_face_y", i, j, 0)
         end do
      end do
      deallocate (a, b)

      call ocean_halo_destroy()
   end subroutine run_1rank

   ! =================================================================
   ! WIDE LEG: nprocs == 4, 2x2 non-periodic.
   ! Uses ng_wide = NG+1 = 4.  nxl=6, nyl=5 for 2x2 of NXG=12, NYG=10.
   ! Fill + exchange at the wide band; check the FULL wide window including
   ! ghost bands out to ng_wide (side bands and corners).
   ! =================================================================

   subroutine run_4rank_wide()
      real(wp), allocatable :: eta_w(:, :), ubt_w(:, :), vbt_w(:, :)
      integer :: nxt_w, nyt_w, ngw

      ngw = NG_WIDE_NP4   ! 4
      call decomp_init(decomp, NXG, NYG, 2, 2, rank)
      nxl = decomp%nx_local
      nyl = decomp%ny_local
      nxt = nxl + 2*NG
      nyt = nyl + 2*NG
      nxt_w = nxl + 2*ngw
      nyt_w = nyl + 2*ngw
      fold_x = .false.
      call ocean_halo_init(decomp, NG, .false., .false.)
      call ocean_halo_reserve(1, ngw)

      allocate (eta_w(nxt_w, nyt_w), ubt_w(nxt_w + 1, nyt_w), vbt_w(nxt_w, nyt_w + 1))

      ! --- wide centre 2D ---
      call fill_centre_2d_wide(eta_w, 0, ngw)
      call ocean_halo_centre_2d_wide(eta_w, ngw, device_resident=.false.)
      call check_centre_2d_wide(eta_w, 0, "w4_centre", ngw)

      ! --- wide face_x 2D ---
      call fill_face_x_2d_wide(ubt_w, 0, ngw)
      call ocean_halo_face_x_2d_wide(ubt_w, ngw, device_resident=.false.)
      call check_face_x_2d_wide(ubt_w, 0, "w4_face_x", ngw)

      ! --- wide face_y 2D ---
      call fill_face_y_2d_wide(vbt_w, 0, ngw)
      call ocean_halo_face_y_2d_wide(vbt_w, ngw, device_resident=.false.)
      call check_face_y_2d_wide(vbt_w, 0, "w4_face_y", ngw)

      ! --- wide bt_group ---
      call fill_centre_2d_wide(eta_w, 0, ngw)
      call fill_face_x_2d_wide(ubt_w, 0, ngw)
      call fill_face_y_2d_wide(vbt_w, 0, ngw)
      call ocean_halo_bt_group_2d_wide(eta_w, ubt_w, vbt_w, ngw, device_resident=.false.)
      call check_centre_2d_wide(eta_w, 0, "w4_bt_eta", ngw)
      call check_face_x_2d_wide(ubt_w, 0, "w4_bt_ubt", ngw)
      call check_face_y_2d_wide(vbt_w, 0, "w4_bt_vbt", ngw)

      call ocean_halo_destroy()
      deallocate (eta_w, ubt_w, vbt_w)
   end subroutine run_4rank_wide

   ! =================================================================
   ! WIDE LEG: nprocs == 2, px=2, periodic-x.
   ! Uses ng_wide = NG+2 = 5.  nxl=6, nyl=10 for px=2 of NXG=12, NYG=10.
   ! Check: the rank's FULL wide local window matches the single-rank
   ! reference BIT-EXACTLY (same as run_2rank_periodic but at ng_wide).
   ! =================================================================

   subroutine run_2rank_periodic_wide()
      real(wp), allocatable :: fld(:, :), ref(:, :)
      integer :: i, j, off, ngw, nxt_w, nyt_w

      ngw = NG_WIDE_NP2   ! 5; nxl=6 >= 5 OK
      call decomp_init(decomp, NXG, NYG, 2, 1, rank)
      nxl = decomp%nx_local
      nyl = decomp%ny_local
      nxt = nxl + 2*NG
      nyt = nyl + 2*NG
      nxt_w = nxl + 2*ngw
      nyt_w = nyl + 2*ngw
      fold_x = .true.
      call ocean_halo_init(decomp, NG, .true., .false.)
      call ocean_halo_reserve(1, ngw)
      off = decomp%i_start - 1

      ! --- wide centre (periodic-x check) ---
      allocate (ref(NXG + 2*ngw, NYG + 2*ngw))
      ref = SENTINEL
      do j = ngw + 1, ngw + NYG
         do i = ngw + 1, ngw + NXG
            ref(i, j) = sig(i - ngw, j - ngw, 0)
         end do
      end do
      call ocean_periodic_wrap_centre_2d(ref, NXG + 2*ngw, NYG + 2*ngw, &
                                         NXG, NYG, ngw, .true., .false.)
      allocate (fld(nxt_w, nyt_w))
      call fill_centre_2d_wide(fld, 0, ngw)
      call ocean_halo_centre_2d_wide(fld, ngw, device_resident=.false.)
      do j = 1, nyt_w
         do i = 1, nxt_w
            call expect_eq(fld(i, j), ref(i + off, j), "wp2_centre", i, j, 0)
         end do
      end do
      deallocate (fld, ref)

      ! --- wide face_x ---
      allocate (ref(NXG + 2*ngw + 1, NYG + 2*ngw))
      ref = SENTINEL
      do j = ngw + 1, ngw + NYG
         do i = ngw + 1, ngw + NXG + 1
            ref(i, j) = sig(i - ngw, j - ngw, 0)
         end do
      end do
      call ocean_periodic_wrap_face_x_2d(ref, NXG + 2*ngw + 1, NYG + 2*ngw, &
                                         NXG, NYG, ngw, .true., .false.)
      allocate (fld(nxt_w + 1, nyt_w))
      call fill_face_x_2d_wide(fld, 0, ngw)
      call ocean_halo_face_x_2d_wide(fld, ngw, device_resident=.false.)
      do j = 1, nyt_w
         do i = 1, nxt_w + 1
            call expect_eq(fld(i, j), ref(i + off, j), "wp2_face_x", i, j, 0)
         end do
      end do
      deallocate (fld, ref)

      ! --- wide face_y ---
      allocate (ref(NXG + 2*ngw, NYG + 2*ngw + 1))
      ref = SENTINEL
      do j = ngw + 1, ngw + NYG + 1
         do i = ngw + 1, ngw + NXG
            ref(i, j) = sig(i - ngw, j - ngw, 0)
         end do
      end do
      call ocean_periodic_wrap_face_y_2d(ref, NXG + 2*ngw, NYG + 2*ngw + 1, &
                                         NXG, NYG, ngw, .true., .false.)
      allocate (fld(nxt_w, nyt_w + 1))
      call fill_face_y_2d_wide(fld, 0, ngw)
      call ocean_halo_face_y_2d_wide(fld, ngw, device_resident=.false.)
      do j = 1, nyt_w + 1
         do i = 1, nxt_w
            call expect_eq(fld(i, j), ref(i + off, j), "wp2_face_y", i, j, 0)
         end do
      end do
      deallocate (fld, ref)

      call ocean_halo_destroy()
   end subroutine run_2rank_periodic_wide

   ! =================================================================
   ! WIDE LEG: nprocs == 1.
   ! Periodic: wide exchange == direct wrap kernel at ng_wide.
   ! Non-periodic: wide exchange is a no-op (bit-unchanged).
   ! ng_wide = NG+2 = 5 for periodic (NXG=NYG_phys, nxl=12, ok).
   ! ng_wide = NG+1 = 4 for non-periodic (no wrap guards needed).
   ! =================================================================

   subroutine run_1rank_wide()
      real(wp), allocatable :: a(:, :), b(:, :)
      integer :: i, j, ngw_p, ngw_np
      integer :: nxt_p, nyt_p, nxt_np, nyt_np

      ! Periodic wide (ng_wide = NG+2 = 5; nxl=12 >= 5 OK for single-rank)
      ngw_p = NG + 2
      nxt_p = NXG + 2*ngw_p
      nyt_p = NYG + 2*ngw_p

      ! Non-periodic wide (ng_wide = NG+1 = 4)
      ngw_np = NG + 1
      nxt_np = NXG + 2*ngw_np
      nyt_np = NYG + 2*ngw_np

      call decomp_init(decomp, NXG, NYG, 1, 1, rank)
      nxl = decomp%nx_local
      nyl = decomp%ny_local
      fold_x = .false.

      ! ---- Periodic (both axes): wide exchange == wrap kernel ----
      call ocean_halo_init(decomp, NG, .true., .true.)
      call ocean_halo_reserve(1, ngw_p)

      ! wide centre
      allocate (a(nxt_p, nyt_p), b(nxt_p, nyt_p))
      call fill_centre_2d_wide(a, 0, ngw_p)
      b = a
      call ocean_halo_centre_2d_wide(a, ngw_p, device_resident=.false.)
      call ocean_periodic_wrap_centre_2d(b, nxt_p, nyt_p, NXG, NYG, ngw_p, .true., .true.)
      do j = 1, nyt_p
         do i = 1, nxt_p
            call expect_eq(a(i, j), b(i, j), "wp1_centre", i, j, 0)
         end do
      end do
      deallocate (a, b)

      ! wide face_x
      allocate (a(nxt_p + 1, nyt_p), b(nxt_p + 1, nyt_p))
      call fill_face_x_2d_wide(a, 0, ngw_p)
      b = a
      call ocean_halo_face_x_2d_wide(a, ngw_p, device_resident=.false.)
      call ocean_periodic_wrap_face_x_2d(b, nxt_p + 1, nyt_p, NXG, NYG, ngw_p, .true., .true.)
      do j = 1, nyt_p
         do i = 1, nxt_p + 1
            call expect_eq(a(i, j), b(i, j), "wp1_face_x", i, j, 0)
         end do
      end do
      deallocate (a, b)

      ! wide face_y
      allocate (a(nxt_p, nyt_p + 1), b(nxt_p, nyt_p + 1))
      call fill_face_y_2d_wide(a, 0, ngw_p)
      b = a
      call ocean_halo_face_y_2d_wide(a, ngw_p, device_resident=.false.)
      call ocean_periodic_wrap_face_y_2d(b, nxt_p, nyt_p + 1, NXG, NYG, ngw_p, .true., .true.)
      do j = 1, nyt_p + 1
         do i = 1, nxt_p
            call expect_eq(a(i, j), b(i, j), "wp1_face_y", i, j, 0)
         end do
      end do
      deallocate (a, b)

      call ocean_halo_destroy()

      ! ---- Non-periodic: wide exchange is a no-op ----
      call ocean_halo_init(decomp, NG, .false., .false.)
      call ocean_halo_reserve(1, ngw_np)

      ! wide centre
      allocate (a(nxt_np, nyt_np), b(nxt_np, nyt_np))
      call fill_centre_2d_wide(a, 0, ngw_np)
      b = a
      call ocean_halo_centre_2d_wide(a, ngw_np, device_resident=.false.)
      do j = 1, nyt_np
         do i = 1, nxt_np
            call expect_eq(a(i, j), b(i, j), "wnp1_centre", i, j, 0)
         end do
      end do
      deallocate (a, b)

      ! wide face_x
      allocate (a(nxt_np + 1, nyt_np), b(nxt_np + 1, nyt_np))
      call fill_face_x_2d_wide(a, 0, ngw_np)
      b = a
      call ocean_halo_face_x_2d_wide(a, ngw_np, device_resident=.false.)
      do j = 1, nyt_np
         do i = 1, nxt_np + 1
            call expect_eq(a(i, j), b(i, j), "wnp1_face_x", i, j, 0)
         end do
      end do
      deallocate (a, b)

      ! wide face_y
      allocate (a(nxt_np, nyt_np + 1), b(nxt_np, nyt_np + 1))
      call fill_face_y_2d_wide(a, 0, ngw_np)
      b = a
      call ocean_halo_face_y_2d_wide(a, ngw_np, device_resident=.false.)
      do j = 1, nyt_np + 1
         do i = 1, nxt_np
            call expect_eq(a(i, j), b(i, j), "wnp1_face_y", i, j, 0)
         end do
      end do
      deallocate (a, b)

      call ocean_halo_destroy()
   end subroutine run_1rank_wide

   ! =================================================================
   ! Wide fills: like fill_*_2d but with an explicit ng_wide argument
   ! =================================================================

   subroutine fill_centre_2d_wide(fld, k, ngw)
      real(wp), intent(inout) :: fld(:, :)
      integer, intent(in) :: k, ngw
      integer :: i, j, gi, gj
      fld = SENTINEL
      do j = ngw + 1, ngw + nyl
         do i = ngw + 1, ngw + nxl
            gi = i - ngw + decomp%i_start - 1
            gj = j - ngw + decomp%j_start - 1
            fld(i, j) = sig(gi, gj, k)
         end do
      end do
   end subroutine fill_centre_2d_wide

   subroutine fill_face_x_2d_wide(fld, k, ngw)
      real(wp), intent(inout) :: fld(:, :)
      integer, intent(in) :: k, ngw
      integer :: i, j, gf, gj
      fld = SENTINEL
      do j = ngw + 1, ngw + nyl
         do i = ngw + 1, ngw + nxl + 1
            gf = i - ngw - 1 + decomp%i_start
            gj = j - ngw + decomp%j_start - 1
            fld(i, j) = sig(gf, gj, k)
         end do
      end do
   end subroutine fill_face_x_2d_wide

   subroutine fill_face_y_2d_wide(fld, k, ngw)
      real(wp), intent(inout) :: fld(:, :)
      integer, intent(in) :: k, ngw
      integer :: i, j, gi, gfj
      fld = SENTINEL
      do j = ngw + 1, ngw + nyl + 1
         do i = ngw + 1, ngw + nxl
            gi = i - ngw + decomp%i_start - 1
            gfj = j - ngw - 1 + decomp%j_start
            fld(i, j) = sig(gi, gfj, k)
         end do
      end do
   end subroutine fill_face_y_2d_wide

   ! =================================================================
   ! Wide checkers (non-periodic np4): full wide-window sweep.
   ! Reachability uses the same dir_ok_* helpers; ghost depth is ngw.
   ! A position is reachable iff each ghost band it sits in has a
   ! live neighbour; corners need BOTH neighbours (D2 transitivity).
   ! =================================================================

   subroutine check_centre_2d_wide(fld, k, label, ngw)
      real(wp), intent(in) :: fld(:, :)
      integer, intent(in) :: k, ngw
      character(len=*), intent(in) :: label
      integer :: i, j, gi, gj, nxt_w, nyt_w
      real(wp) :: expect
      nxt_w = nxl + 2*ngw
      nyt_w = nyl + 2*ngw
      do j = 1, nyt_w
         do i = 1, nxt_w
            gi = i - ngw + decomp%i_start - 1
            gj = j - ngw + decomp%j_start - 1
            if (dir_ok_x(zone_w(i, nxl, ngw)) .and. dir_ok_y(zone_w(j, nyl, ngw))) then
               expect = sig(gi, gj, k)
            else
               expect = SENTINEL
            end if
            call expect_eq(fld(i, j), expect, label, i, j, k)
         end do
      end do
   end subroutine check_centre_2d_wide

   subroutine check_face_x_2d_wide(fld, k, label, ngw)
      real(wp), intent(in) :: fld(:, :)
      integer, intent(in) :: k, ngw
      character(len=*), intent(in) :: label
      integer :: i, j, gf, gj, nxt_w, nyt_w
      real(wp) :: expect
      nxt_w = nxl + 2*ngw
      nyt_w = nyl + 2*ngw
      do j = 1, nyt_w
         do i = 1, nxt_w + 1
            gf = i - ngw - 1 + decomp%i_start
            gj = j - ngw + decomp%j_start - 1
            if (dir_ok_x(zone_w(i, nxl + 1, ngw)) .and. dir_ok_y(zone_w(j, nyl, ngw))) then
               expect = sig(gf, gj, k)
            else
               expect = SENTINEL
            end if
            call expect_eq(fld(i, j), expect, label, i, j, k)
         end do
      end do
   end subroutine check_face_x_2d_wide

   subroutine check_face_y_2d_wide(fld, k, label, ngw)
      real(wp), intent(in) :: fld(:, :)
      integer, intent(in) :: k, ngw
      character(len=*), intent(in) :: label
      integer :: i, j, gi, gfj, nxt_w, nyt_w
      real(wp) :: expect
      nxt_w = nxl + 2*ngw
      nyt_w = nyl + 2*ngw
      do j = 1, nyt_w + 1
         do i = 1, nxt_w
            gi = i - ngw + decomp%i_start - 1
            gfj = j - ngw - 1 + decomp%j_start
            if (dir_ok_x(zone_w(i, nxl, ngw)) .and. dir_ok_y(zone_w(j, nyl + 1, ngw))) then
               expect = sig(gi, gfj, k)
            else
               expect = SENTINEL
            end if
            call expect_eq(fld(i, j), expect, label, i, j, k)
         end do
      end do
   end subroutine check_face_y_2d_wide

   pure function zone_w(idx, n_phys, ngw) result(z)
      !! Zone function generalised for arbitrary ghost width ngw.
      integer, intent(in) :: idx, n_phys, ngw
      integer :: z
      if (idx <= ngw) then
         z = -1
      else if (idx > ngw + n_phys) then
         z = 1
      else
         z = 0
      end if
   end function zone_w

end program test_halo_ocean_mpi
#endif
