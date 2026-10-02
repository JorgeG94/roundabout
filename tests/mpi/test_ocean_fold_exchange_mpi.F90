!! MPI test of the distributed tripolar north-fold exchange
!! (`rdb_ocean_fold_exchange`) against the serial fold kernels.
!!
!! Every rank builds the SAME global, index-encoded, x-periodic field
!! (value = fid*1e6 + c*1e3 + J*10 + k + 1/4, exact in wp; c the reduced
!! physical column, J the global storage row), folds it with the serial
!! `rdb_ocean_fold` kernel (the 1-rank reference), and folds its own tile
!! window through the px-dispatching entry points.  Two legs per field:
!!   (a) OWNER-ONLY: the window holds only the tile's OWNED points (T/v:
!!       cells ng+1..ng+w; u/corner: faces ng+2..ng+w+1, never the D1
!!       west-seam copy), everything else NaN — every slot the fold writes
!!       (all north ghost rows; the fold-line row's west half and, for a
!!       vector, its self-conjugate column) must equal the reference
!!       bitwise.  Proves the owner routing needs no prior x exchange.
!!   (b) FULL WINDOW: x ghosts as the x exchange leaves them, north ghosts
!!       NaN — the whole storage window must equal the reference bitwise
!!       (nothing else is touched; east-half fold-row slots stay).
!! Fields: T, u, v (2D + 3D), corner (copy + negate), and one GROUP (T3d,
!! u3d, v3d, corner2d negated, T2d in one message per peer).
!!
!! Rank-count legs (x-only splits unless noted; ng = 3):
!!   np 1: px = 1 — the dispatcher must hit the local kernels (leg (b)).
!!   np 2: ni 30, 31 (odd: self-conjugate v column), 8.
!!   np 3: ni 30 (self-mirror middle tile), 32 (straddling self-mirror),
!!         31, 10.
!!   np 4: 4x1 with ni 30 (8/8/7/7), 31, 13 (tiles narrower than ng+1),
!!         and 2x2 (only the north rank row folds).
!! Host mode (`device_resident=.false.`): the arrays are host locals.  The
!! device path is exercised end-to-end by the engine tests on the GPU build.
#ifdef RDB_ENABLE_MPI
program test_ocean_fold_exchange_mpi
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use rdb_constants, only: wp
   use rdb_decomp, only: decomp_t, decomp_init
   use rdb_ocean_fold, only: fold_north_centre, fold_north_u_face, &
                             fold_north_v_face, fold_north_corner
   use rdb_ocean_fold_exchange, only: ocean_fold_exchange_init, &
                                      ocean_fold_exchange_destroy, &
                                      ocean_fold_exchange_reserve, &
                                      ocean_fold_is_distributed, &
                                      ocean_fold_begin, ocean_fold_pack, &
                                      ocean_fold_exchange, ocean_fold_unpack, &
                                      ocean_fold_end, ocean_fold_north_centre, &
                                      ocean_fold_north_u_face, ocean_fold_north_v_face, &
                                      ocean_fold_north_corner, FOLD_STAG_T, FOLD_STAG_U, &
                                      FOLD_STAG_V, FOLD_STAG_CORNER
   use rdb_comm_env, only: comm_env_init, comm_env_setup_roles, &
                           comm_env_finalize, comm_env_rank, comm_env_size, &
                           comm_env_compute_comm
   use pic_mpi_lib, only: comm_t, allreduce, MPI_SUM
   implicit none

   integer, parameter :: NG = 3
      !! Ghost width (the production minimum).
   integer, parameter :: NJ_G = 9
      !! Global rows: a 2-row split leaves a north tile of 4 = ng+1 rows.
   integer, parameter :: NZ = 3
      !! Layers of the 3D legs.

   integer :: rank, nprocs, n_fail, total_fail, n_checked
   type(comm_t) :: comm
   real(wp) :: nan

   call comm_env_init()
   call comm_env_setup_roles(.false.)
   rank = comm_env_rank()
   nprocs = comm_env_size()
   n_fail = 0
   n_checked = 0
   nan = ieee_value(1.0_wp, ieee_quiet_nan)

   select case (nprocs)
   case (1)
      call run_case(30, 1, 1)
      call run_case(31, 1, 1)
   case (2)
      call run_case(30, 2, 1)
      call run_case(31, 2, 1)
      call run_case(8, 2, 1)
   case (3)
      call run_case(30, 3, 1)
      call run_case(32, 3, 1)
      call run_case(31, 3, 1)
      call run_case(10, 3, 1)
   case (4)
      call run_case(30, 4, 1)
      call run_case(31, 4, 1)
      call run_case(13, 4, 1)
      call run_case(30, 2, 2)
   case default
      if (rank == 0) write (*, *) "SKIP: test supports 1/2/3/4 ranks, got", nprocs
   end select

   comm = comm_env_compute_comm()
   call comm%barrier()
   if (n_fail > 0) then
      write (*, *) "Rank", rank, ":", n_fail, "check(s) FAILED (nprocs=", nprocs, ")"
   else
      write (*, *) "Rank", rank, ": all checks PASSED (nprocs=", nprocs, ", slots=", &
         n_checked, ")"
   end if
   total_fail = n_fail
   call allreduce(comm, total_fail, MPI_SUM)
   call comm_env_finalize()
   if (total_fail > 0) error stop 1

contains

   subroutine run_case(ni, px, py)
      integer, intent(in) :: ni, px, py
      type(decomp_t) :: d

      call decomp_init(d, ni, NJ_G, px, py, rank)
      call ocean_fold_exchange_init(d, NG, d%has_north)
      ! Largest group below: T3d + u3d + v3d + corner2d + T2d.
      call ocean_fold_exchange_reserve((NG + 1)*(3*NZ + 2))
      if (d%has_north .and. (px > 1 .neqv. ocean_fold_is_distributed())) then
         call fail_msg(ni, px, py, 0, "ocean_fold_is_distributed disagrees with px")
      end if

      if (d%has_north) then
         call one_field(d, ni, px, py, FOLD_STAG_T, 1, .false., 1)
         call one_field(d, ni, px, py, FOLD_STAG_T, NZ, .false., 2)
         call one_field(d, ni, px, py, FOLD_STAG_U, 1, .true., 3)
         call one_field(d, ni, px, py, FOLD_STAG_U, NZ, .true., 4)
         call one_field(d, ni, px, py, FOLD_STAG_V, 1, .true., 5)
         call one_field(d, ni, px, py, FOLD_STAG_V, NZ, .true., 6)
         call one_field(d, ni, px, py, FOLD_STAG_CORNER, 1, .false., 7)
         call one_field(d, ni, px, py, FOLD_STAG_CORNER, 1, .true., 8)
         if (px > 1) call group_leg(d, ni, px, py)
      end if
      call ocean_fold_exchange_destroy()
   end subroutine run_case

   ! ----------------------------------------------------------------
   ! Geometry helpers
   ! ----------------------------------------------------------------

   pure subroutine stag_ext(stag, xext, yext)
      integer, intent(in) :: stag
      integer, intent(out) :: xext, yext
      xext = merge(1, 0, stag == FOLD_STAG_U .or. stag == FOLD_STAG_CORNER)
      yext = merge(1, 0, stag == FOLD_STAG_V .or. stag == FOLD_STAG_CORNER)
   end subroutine stag_ext

   pure real(wp) function encode(fid, c, jg, k) result(v)
      integer, intent(in) :: fid, c, jg, k
      v = real(fid, wp)*1.0e6_wp + real(c, wp)*1.0e3_wp + real(jg, wp)*10.0_wp + &
          real(k, wp) + 0.25_wp
   end function encode

   subroutine make_global(stag, nz, negate, fid, ni, g, ref)
      !! Pre-fold global field (north ghosts NaN) and its serial fold.
      integer, intent(in) :: stag, nz, fid, ni
      logical, intent(in) :: negate
      real(wp), allocatable, intent(out) :: g(:, :, :), ref(:, :, :)
      integer :: xext, yext, nxg, nyg, i, j, k

      call stag_ext(stag, xext, yext)
      nxg = ni + 2*NG + xext
      nyg = NJ_G + 2*NG + yext
      allocate (g(nxg, nyg, nz), ref(nxg, nyg, nz))
      do k = 1, nz
         do j = 1, nyg
            do i = 1, nxg
               if (j > NG + NJ_G + yext) then
                  g(i, j, k) = nan
               else
                  g(i, j, k) = encode(fid, modulo(i - NG - 1, ni) + 1, j, k)
               end if
            end do
         end do
      end do
      ref = g
      select case (stag)
      case (FOLD_STAG_T)
         call fold_north_centre(ref, nxg, nyg, nz, ni, NJ_G, NG)
      case (FOLD_STAG_U)
         call fold_north_u_face(ref, nxg, nyg, nz, ni, NJ_G, NG)
      case (FOLD_STAG_V)
         call fold_north_v_face(ref, nxg, nyg, nz, ni, NJ_G, NG)
      case (FOLD_STAG_CORNER)
         call fold_north_corner(ref(:, :, 1), nxg, nyg, ni, NJ_G, NG, negate)
      end select
   end subroutine make_global

   subroutine make_window(d, stag, g, owner_only, loc)
      !! This tile's storage window of the pre-fold global field.
      type(decomp_t), intent(in) :: d
      integer, intent(in) :: stag
      real(wp), intent(in) :: g(:, :, :)
      logical, intent(in) :: owner_only
      real(wp), allocatable, intent(out) :: loc(:, :, :)
      integer :: xext, yext, nxa, nya, i, j, k, lo, hi

      call stag_ext(stag, xext, yext)
      nxa = d%nx_local + 2*NG + xext
      nya = d%ny_local + 2*NG + yext
      allocate (loc(nxa, nya, size(g, 3)))
      do k = 1, size(g, 3)
         do j = 1, nya
            do i = 1, nxa
               loc(i, j, k) = g(d%i_start - 1 + i, d%j_start - 1 + j, k)
            end do
         end do
      end do
      if (owner_only) then
         ! Owned points only; rows above the last pre-fold row are already NaN.
         lo = NG + 1 + xext
         hi = NG + d%nx_local + xext
         do k = 1, size(loc, 3)
            do j = 1, nya
               do i = 1, nxa
                  if (i < lo .or. i > hi .or. j <= NG) loc(i, j, k) = nan
               end do
            end do
         end do
      end if
   end subroutine make_window

   logical function fold_writes(d, stag, ni, negate, i, j) result(w)
      !! Does the fold write storage slot (i, j) of this tile's window?
      type(decomp_t), intent(in) :: d
      integer, intent(in) :: stag, ni, i, j
      logical, intent(in) :: negate
      integer :: xext, yext, p, pm

      call stag_ext(stag, xext, yext)
      if (j > NG + d%ny_local + yext) then
         w = .true.
      else if (yext == 1 .and. j == NG + d%ny_local + 1) then
         p = modulo(d%i_start - 1 + i - NG - 1, ni) + 1
         if (xext == 1) then
            pm = modulo(ni + 1 - p, ni) + 1
         else
            pm = ni + 1 - p
         end if
         w = (p < pm) .or. (p == pm .and. negate)
      else
         w = .false.
      end if
   end function fold_writes

   subroutine fold_dispatch(stag, nz, negate, loc, d, ni)
      integer, intent(in) :: stag, nz, ni
      logical, intent(in) :: negate
      real(wp), intent(inout) :: loc(:, :, :)
      type(decomp_t), intent(in) :: d
      integer :: nxa, nya

      nxa = size(loc, 1)
      nya = size(loc, 2)
      select case (stag)
      case (FOLD_STAG_T)
         if (nz == 1) then
            call ocean_fold_north_centre(loc(:, :, 1), nxa, nya, ni, d%ny_local, NG, &
                                         device_resident=.false.)
         else
            call ocean_fold_north_centre(loc, nxa, nya, nz, ni, d%ny_local, NG, &
                                         device_resident=.false.)
         end if
      case (FOLD_STAG_U)
         if (nz == 1) then
            call ocean_fold_north_u_face(loc(:, :, 1), nxa, nya, ni, d%ny_local, NG, &
                                         device_resident=.false.)
         else
            call ocean_fold_north_u_face(loc, nxa, nya, nz, ni, d%ny_local, NG, &
                                         device_resident=.false.)
         end if
      case (FOLD_STAG_V)
         if (nz == 1) then
            call ocean_fold_north_v_face(loc(:, :, 1), nxa, nya, ni, d%ny_local, NG, &
                                         device_resident=.false.)
         else
            call ocean_fold_north_v_face(loc, nxa, nya, nz, ni, d%ny_local, NG, &
                                         device_resident=.false.)
         end if
      case (FOLD_STAG_CORNER)
         call ocean_fold_north_corner(loc(:, :, 1), nxa, nya, ni, d%ny_local, NG, negate, &
                                      device_resident=.false.)
      end select
   end subroutine fold_dispatch

   subroutine compare(d, stag, ni, negate, loc, ref, owner_only, px, py, fid)
      type(decomp_t), intent(in) :: d
      integer, intent(in) :: stag, ni, px, py, fid
      logical, intent(in) :: negate, owner_only
      real(wp), intent(in) :: loc(:, :, :), ref(:, :, :)
      integer :: i, j, k
      character(len=64) :: where

      do k = 1, size(loc, 3)
         do j = 1, size(loc, 2)
            do i = 1, size(loc, 1)
               if (owner_only) then
                  if (.not. fold_writes(d, stag, ni, negate, i, j)) cycle
               end if
               n_checked = n_checked + 1
               if (.not. (loc(i, j, k) == ref(d%i_start - 1 + i, d%j_start - 1 + j, k))) then
                  write (where, '(a,l1,a,i0,a,i0,a,i0)') "owner_only=", owner_only, &
                     " i=", i, " j=", j, " k=", k
                  call fail_msg(ni, px, py, fid, trim(where))
                  return
               end if
            end do
         end do
      end do
   end subroutine compare

   subroutine one_field(d, ni, px, py, stag, nz, negate, fid)
      type(decomp_t), intent(in) :: d
      integer, intent(in) :: ni, px, py, stag, nz, fid
      logical, intent(in) :: negate
      real(wp), allocatable :: g(:, :, :), ref(:, :, :), loc(:, :, :)

      call make_global(stag, nz, negate, fid, ni, g, ref)
      ! (b) full window.
      call make_window(d, stag, g, .false., loc)
      call fold_dispatch(stag, nz, negate, loc, d, ni)
      call compare(d, stag, ni, negate, loc, ref, .false., px, py, fid)
      ! (a) owner-only window: the exchange path must not read anything
      ! but owned points (px = 1 keeps the local kernels, which read the
      ! periodic ghost columns by design — leg (b) covers it).
      if (px > 1) then
         call make_window(d, stag, g, .true., loc)
         call fold_dispatch(stag, nz, negate, loc, d, ni)
         call compare(d, stag, ni, negate, loc, ref, .true., px, py, fid)
      end if
   end subroutine one_field

   subroutine group_leg(d, ni, px, py)
      !! Five fields of mixed stagger and family in ONE message per peer.
      type(decomp_t), intent(in) :: d
      integer, intent(in) :: ni, px, py
      real(wp), allocatable :: gt3(:, :, :), rt3(:, :, :), lt3(:, :, :)
      real(wp), allocatable :: gu3(:, :, :), ru3(:, :, :), lu3(:, :, :)
      real(wp), allocatable :: gv3(:, :, :), rv3(:, :, :), lv3(:, :, :)
      real(wp), allocatable :: gc2(:, :, :), rc2(:, :, :), lc2(:, :, :)
      real(wp), allocatable :: gt2(:, :, :), rt2(:, :, :), lt2(:, :, :)
      integer :: leg
      logical :: own

      call make_global(FOLD_STAG_T, NZ, .false., 11, ni, gt3, rt3)
      call make_global(FOLD_STAG_U, NZ, .true., 12, ni, gu3, ru3)
      call make_global(FOLD_STAG_V, NZ, .true., 13, ni, gv3, rv3)
      call make_global(FOLD_STAG_CORNER, 1, .true., 14, ni, gc2, rc2)
      call make_global(FOLD_STAG_T, 1, .false., 15, ni, gt2, rt2)
      do leg = 1, 2
         own = (leg == 2)
         call make_window(d, FOLD_STAG_T, gt3, own, lt3)
         call make_window(d, FOLD_STAG_U, gu3, own, lu3)
         call make_window(d, FOLD_STAG_V, gv3, own, lv3)
         call make_window(d, FOLD_STAG_CORNER, gc2, own, lc2)
         call make_window(d, FOLD_STAG_T, gt2, own, lt2)
         call ocean_fold_begin((NG + 1)*(3*NZ + 2))
         call ocean_fold_pack(lt3, size(lt3, 1), size(lt3, 2), NZ, FOLD_STAG_T, .false.)
         call ocean_fold_pack(lu3, size(lu3, 1), size(lu3, 2), NZ, FOLD_STAG_U, .false.)
         call ocean_fold_pack(lv3, size(lv3, 1), size(lv3, 2), NZ, FOLD_STAG_V, .false.)
         call ocean_fold_pack(lc2(:, :, 1), size(lc2, 1), size(lc2, 2), FOLD_STAG_CORNER, .false.)
         call ocean_fold_pack(lt2(:, :, 1), size(lt2, 1), size(lt2, 2), FOLD_STAG_T, .false.)
         call ocean_fold_exchange(.false.)
         call ocean_fold_unpack(lt3, size(lt3, 1), size(lt3, 2), NZ, FOLD_STAG_T, .false., .false.)
         call ocean_fold_unpack(lu3, size(lu3, 1), size(lu3, 2), NZ, FOLD_STAG_U, .true., .false.)
         call ocean_fold_unpack(lv3, size(lv3, 1), size(lv3, 2), NZ, FOLD_STAG_V, .true., .false.)
         call ocean_fold_unpack(lc2(:, :, 1), size(lc2, 1), size(lc2, 2), FOLD_STAG_CORNER, &
                                .true., .false.)
         call ocean_fold_unpack(lt2(:, :, 1), size(lt2, 1), size(lt2, 2), FOLD_STAG_T, &
                                .false., .false.)
         call ocean_fold_end()
         call compare(d, FOLD_STAG_T, ni, .false., lt3, rt3, own, px, py, 11)
         call compare(d, FOLD_STAG_U, ni, .true., lu3, ru3, own, px, py, 12)
         call compare(d, FOLD_STAG_V, ni, .true., lv3, rv3, own, px, py, 13)
         call compare(d, FOLD_STAG_CORNER, ni, .true., lc2, rc2, own, px, py, 14)
         call compare(d, FOLD_STAG_T, ni, .false., lt2, rt2, own, px, py, 15)
      end do
   end subroutine group_leg

   subroutine fail_msg(ni, px, py, fid, what)
      integer, intent(in) :: ni, px, py, fid
      character(len=*), intent(in) :: what
      n_fail = n_fail + 1
      write (*, '(a,i0,a,i0,a,i0,a,i0,a,i0,2a)') "FAIL rank ", rank, ": ni=", ni, " px=", px, &
         " py=", py, " fid=", fid, " ", what
   end subroutine fail_msg

end program test_ocean_fold_exchange_mpi
#endif
