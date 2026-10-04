!! Serial tests of the distributed-fold routing plan (`rdb_ocean_fold_plan`).
!!
!! No MPI: every tile's plan is a pure function of (ni, px, ng, rx), so the
!! whole north rank row is emulated in one process.  Three layers:
!!   1. the plan's tile split / owner closed forms agree with `decomp_init`
!!      and a brute-force owner search;
!!   2. every (sender, receiver) list pair agrees entry-for-entry, covers
!!      each destination column exactly once, and its source is the owned
!!      mirror point;
!!   3. an emulated exchange over index-encoded fields — senders expose
!!      ONLY their owned points (everything else NaN-poisoned) — reproduces
!!      the serial fold kernels of `rdb_ocean_fold` bit for bit on every
!!      storage column of every tile window, for all four staggers, odd and
!!      even ni, odd px and uneven tile widths.
module test_ocean_fold_plan
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_decomp, only: decomp_t, decomp_init
   use rdb_ocean_fold, only: fold_north_centre, fold_north_u_face, &
                             fold_north_v_face, fold_north_corner
   use rdb_ocean_fold_plan, only: fold_plan_t, fold_plan_build, fold_tile_extent, &
                                  fold_tile_owner, &
                                  fold_stagger_family, fold_stagger_nrows, fold_row_map, &
                                  FOLD_STAG_T, FOLD_STAG_U, FOLD_STAG_V, FOLD_STAG_CORNER, &
                                  FOLD_FAM_T, FOLD_FAM_U, FOLD_NFAM, FOLD_ROW_WEST, &
                                  FOLD_ROW_SELF, FOLD_PLAN_OK, &
                                  FOLD_PLAN_ERR_ARGS
   implicit none
   private

   public :: collect_ocean_fold_plan_tests

   ! Sweep: even and odd ni (odd ni has a self-conjugate v column), the
   ! plan-§2.4 cases (30/4 uneven 8/8/7/7, 30/3 self-mirror middle tile,
   ! 32/3 straddling self-mirror, 10/3) and widths below ng+1.
   integer, parameter :: NI_LIST(*) = [8, 9, 10, 11, 30, 31, 32]
   integer, parameter :: PX_MAX = 6
   integer, parameter :: NG_LIST(*) = [1, 2, 3]

contains

   subroutine collect_ocean_fold_plan_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("fold_plan_tile_extent_matches_decomp", test_tile_extent), &
                  new_unittest("fold_plan_owner_brute_force", test_owner), &
                  new_unittest("fold_plan_rejects_bad_args", test_bad_args), &
                  new_unittest("fold_plan_pairs_consistent", test_pairs), &
                  new_unittest("fold_plan_section_2_4_cases", test_cases_2_4), &
                  new_unittest("fold_plan_emulated_exchange_bitwise", test_emulated) &
                  ]
   end subroutine collect_ocean_fold_plan_tests

   subroutine test_tile_extent(error)
      type(error_type), allocatable, intent(out) :: error
      type(decomp_t) :: d
      integer :: ni, px, rx, a, w

      do ni = 1, 40
         do px = 1, min(ni, 9)
            do rx = 0, px - 1
               call decomp_init(d, ni, 4, px, 1, rx)
               call fold_tile_extent(ni, px, rx, a, w)
               call check(error, a == d%i_start .and. w == d%nx_local, &
                          "fold_tile_extent disagrees with decomp_init")
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_tile_extent

   subroutine test_owner(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: ni, px, c, rx, a, w, brute

      do ni = 1, 40
         do px = 1, min(ni, 9)
            do c = 1, ni
               brute = -1
               do rx = 0, px - 1
                  call fold_tile_extent(ni, px, rx, a, w)
                  if (c >= a .and. c <= a + w - 1) brute = rx
               end do
               call check(error, fold_tile_owner(ni, px, c) == brute, &
                          "fold_tile_owner disagrees with the brute-force search")
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_owner

   subroutine test_bad_args(error)
      type(error_type), allocatable, intent(out) :: error
      type(fold_plan_t) :: plan
      integer :: status
      logical :: freed

      call fold_plan_build(plan, 3, 4, 3, 0, status)       ! px > ni
      call check(error, status == FOLD_PLAN_ERR_ARGS, "px > ni must be refused")
      if (allocated(error)) return
      call fold_plan_build(plan, 30, 4, 3, 4, status)      ! rx out of range
      call check(error, status == FOLD_PLAN_ERR_ARGS, "rx >= px must be refused")
      if (allocated(error)) return
      call fold_plan_build(plan, 30, 4, 0, 0, status)      ! ng = 0
      call check(error, status == FOLD_PLAN_ERR_ARGS, "ng = 0 must be refused")
      if (allocated(error)) return
      call fold_plan_build(plan, 30, 4, 3, 3, status)
      call check(error, status == FOLD_PLAN_OK, "a valid plan must build")
      if (allocated(error)) return
      call plan%destroy()
      freed = .not. (allocated(plan%send_col) .or. allocated(plan%peer_rx))
      call check(error, freed .and. plan%npeer == 0, "destroy must free the plan's lists")
      if (allocated(error)) return
      call plan%destroy()   ! idempotent on an empty plan
   end subroutine test_bad_args

   ! Every (sender, receiver) list pair agrees, every destination column is
   ! covered exactly once, and each source is the owned mirror point.
   subroutine test_pairs(error)
      type(error_type), allocatable, intent(out) :: error
      type(fold_plan_t), allocatable :: plans(:)
      integer :: ii, ig, ni, px, ng, r, s, fam, status, e, kr, ks, pr, ps
      integer :: a_r, w_r, a_s, w_s, ncol, gdst, gsrc, c, m, nent
      integer, allocatable :: hits(:)

      do ii = 1, size(NI_LIST)
         ni = NI_LIST(ii)
         do px = 1, min(PX_MAX, ni)
            do ig = 1, size(NG_LIST)
               ng = NG_LIST(ig)
               if (allocated(plans)) deallocate (plans)
               allocate (plans(0:px - 1))
               do r = 0, px - 1
                  call fold_plan_build(plans(r), ni, px, ng, r, status)
                  call check(error, status == FOLD_PLAN_OK, "plan build failed")
                  if (allocated(error)) return
               end do
               do r = 0, px - 1
                  call fold_tile_extent(ni, px, r, a_r, w_r)
                  do fam = 1, FOLD_NFAM
                     ncol = w_r + 2*ng
                     if (fam == FOLD_FAM_U) ncol = ncol + 1
                     call check(error, plans(r)%nrecv(fam) == ncol, &
                                "receive lists must cover the whole window")
                     if (allocated(error)) return
                     if (allocated(hits)) deallocate (hits)
                     allocate (hits(ncol))
                     hits = 0
                     do pr = 1, plans(r)%npeer
                        s = plans(r)%peer_rx(pr)
                        call fold_tile_extent(ni, px, s, a_s, w_s)
                        ps = peer_index(plans(s), r)
                        nent = plans(r)%recv_n(pr, fam)
                        if (nent == 0) cycle
                        call check(error, ps > 0, "sender does not list the receiver")
                        if (allocated(error)) return
                        call check(error, plans(s)%send_n(ps, fam) == nent, &
                                   "send/receive counts disagree")
                        if (allocated(error)) return
                        do e = 1, nent
                           kr = plans(r)%recv_start(pr, fam) + e
                           ks = plans(s)%send_start(ps, fam) + e
                           call check(error, plans(r)%recv_e(kr) == e .and. &
                                      plans(s)%send_e(ks) == e .and. &
                                      plans(r)%recv_peer(kr) == pr .and. &
                                      plans(s)%send_peer(ks) == ps, "entry bookkeeping")
                           if (allocated(error)) return
                           hits(plans(r)%recv_col(kr)) = hits(plans(r)%recv_col(kr)) + 1
                           ! Destination global column (reduced) and its mirror.
                           gdst = a_r + plans(r)%recv_col(kr) - ng - 1
                           c = modulo(gdst - 1, ni) + 1
                           gsrc = a_s + plans(s)%send_col(ks) - ng - 1
                           if (fam == FOLD_FAM_U) then
                              m = ni + 2 - c
                              ! Owned faces only: ng+2 .. ng+w+1 (D1).
                              call check(error, plans(s)%send_col(ks) >= ng + 2 .and. &
                                         plans(s)%send_col(ks) <= ng + w_s + 1, &
                                         "u/corner source is not an owned face")
                           else
                              m = ni + 1 - c
                              call check(error, plans(s)%send_col(ks) >= ng + 1 .and. &
                                         plans(s)%send_col(ks) <= ng + w_s, &
                                         "T/v source is not an owned cell")
                           end if
                           if (allocated(error)) return
                           call check(error, gsrc == m, "source is not the mirror point")
                           if (allocated(error)) return
                        end do
                     end do
                     call check(error, all(hits == 1), &
                                "a destination column is not covered exactly once")
                     if (allocated(error)) return
                  end do
               end do
            end do
         end do
      end do
   end subroutine test_pairs

   ! The worked cases of plan §2.4 (ng = 3): peer sets as stated there.
   subroutine test_cases_2_4(error)
      type(error_type), allocatable, intent(out) :: error
      type(fold_plan_t) :: plan
      integer :: status

      ! ni = 30, px = 4 (8/8/7/7): tile 0's window draws from itself
      ! (its ghost band wraps), tile 2 and tile 3.
      call fold_plan_build(plan, 30, 4, 3, 0, status)
      call check(error, same_list(recv_peers(plan, FOLD_FAM_T), [0, 2, 3]), &
                 "30/4 tile 0 T peers")
      if (allocated(error)) return
      ! ni = 30, px = 3: the middle tile mirrors onto itself; its ghosts
      ! reach both neighbours.
      call fold_plan_build(plan, 30, 3, 3, 1, status)
      call check(error, same_list(recv_peers(plan, FOLD_FAM_T), [0, 1, 2]), &
                 "30/3 middle tile T peers")
      if (allocated(error)) return
      ! Physical mirror of the middle tile is itself only: count the
      ! self-routed entries among its physical columns.
      call check(error, plan%recv_n(plan%self_peer, FOLD_FAM_T) >= 10, &
                 "30/3 middle tile must self-route its physical columns")
      if (allocated(error)) return
      ! ni = 10, px = 3 (4/3/3): tile 0's physical cells 1..4 mirror to
      ! 7..10, held by tiles 1 and 2.
      call fold_plan_build(plan, 10, 3, 3, 0, status)
      call check(error, any(recv_peers(plan, FOLD_FAM_T) == 1) .and. &
                 any(recv_peers(plan, FOLD_FAM_T) == 2), "10/3 tile 0 straddles tiles 1 and 2")
   end subroutine test_cases_2_4

   ! Emulated exchange vs the serial fold kernels, bitwise.
   subroutine test_emulated(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: ii, ig, ni, px, ng, stag
      logical :: neg

      do ii = 1, size(NI_LIST)
         ni = NI_LIST(ii)
         do px = 1, min(PX_MAX, ni)
            do ig = 1, size(NG_LIST)
               ng = NG_LIST(ig)
               do stag = FOLD_STAG_T, FOLD_STAG_CORNER
                  ! T copies and u/v negate in the serial kernels; the
                  ! corner kernel takes either.
                  neg = (stag == FOLD_STAG_U .or. stag == FOLD_STAG_V)
                  call emulate(error, stag, neg, ni, px, ng)
                  if (allocated(error)) return
                  if (stag == FOLD_STAG_CORNER) then
                     call emulate(error, stag, .true., ni, px, ng)
                     if (allocated(error)) return
                  end if
               end do
            end do
         end do
      end do
   end subroutine test_emulated

   subroutine emulate(error, stag, neg, ni, px, ng)
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: stag, ni, px, ng
      logical, intent(in) :: neg

      type(fold_plan_t), allocatable :: plans(:)
      real(wp), allocatable :: g(:, :), ref(:, :), loc(:, :), snd(:, :)
      real(wp) :: nan, val
      integer :: nyl, nxg, nyg, xext, yext, ilast_row, status, fam, nrow
      integer :: r, s, a, w, a_s, w_s, i, j, pr, ps, e, kr, ks, rr, srow, drow
      integer :: own_lo, own_hi
      character(len=96) :: what

      nan = ieee_value(1.0_wp, ieee_quiet_nan)
      ! Shortest legal north tile: ng+1 rows (the fold reads ng rows below
      ! the fold line, plus the v fold row).
      nyl = ng + 1
      fam = fold_stagger_family(stag)
      nrow = fold_stagger_nrows(stag, ng)
      xext = 0
      if (fam == FOLD_FAM_U) xext = 1
      yext = 0
      if (stag == FOLD_STAG_V .or. stag == FOLD_STAG_CORNER) yext = 1
      nxg = ni + 2*ng + xext
      nyg = nyl + 2*ng + yext
      ! Last pre-fold storage row: the last T row, or the fold-line row.
      ilast_row = ng + nyl + yext
      write (what, '(a,i0,a,i0,a,i0,a,i0,a,l1)') "stag=", stag, " ni=", ni, " px=", px, &
         " ng=", ng, " neg=", neg

      ! Global array: index-encoded, periodic in x by construction
      ! (every storage column holds its reduced physical column's value,
      ! so face ni+1 equals face 1 bitwise), north ghosts poisoned.
      allocate (g(nxg, nyg), ref(nxg, nyg))
      do j = 1, nyg
         do i = 1, nxg
            if (j > ilast_row) then
               g(i, j) = nan
            else
               g(i, j) = encode(modulo(i - ng - 1, ni) + 1, j)
            end if
         end do
      end do
      ref = g
      select case (stag)
      case (FOLD_STAG_T)
         call fold_north_centre(ref, nxg, nyg, ni, nyl, ng)
      case (FOLD_STAG_U)
         call fold_north_u_face(ref, nxg, nyg, ni, nyl, ng)
      case (FOLD_STAG_V)
         call fold_north_v_face(ref, nxg, nyg, ni, nyl, ng)
      case (FOLD_STAG_CORNER)
         call fold_north_corner(ref, nxg, nyg, ni, nyl, ng, neg)
      end select

      allocate (plans(0:px - 1))
      do r = 0, px - 1
         call fold_plan_build(plans(r), ni, px, ng, r, status)
      end do

      do r = 0, px - 1
         call fold_tile_extent(ni, px, r, a, w)
         ! Receiver window: what the x exchange leaves (every column, rows
         ! up to the fold row), north ghosts poisoned.
         allocate (loc(w + 2*ng + xext, nyg))
         do j = 1, nyg
            do i = 1, size(loc, 1)
               loc(i, j) = g(a - 1 + i, j)
            end do
         end do
         do pr = 1, plans(r)%npeer
            s = plans(r)%peer_rx(pr)
            if (plans(r)%recv_n(pr, fam) == 0) cycle
            ps = peer_index(plans(s), r)
            ! Sender window exposing ONLY its owned points.
            call fold_tile_extent(ni, px, s, a_s, w_s)
            allocate (snd(w_s + 2*ng + xext, nyg))
            snd = nan
            own_lo = ng + 1 + xext
            own_hi = ng + w_s + xext
            do j = ng + 1, ilast_row
               do i = own_lo, own_hi
                  snd(i, j) = g(a_s - 1 + i, j)
               end do
            end do
            do e = 1, plans(r)%recv_n(pr, fam)
               kr = plans(r)%recv_start(pr, fam) + e
               ks = plans(s)%send_start(ps, fam) + e
               do rr = 1, nrow
                  call fold_row_map(stag, ng, nyl, rr, srow, drow)
                  val = snd(plans(s)%send_col(ks), srow)
                  if (neg) val = -val
                  if (srow == drow) then
                     ! Fold-line row: projection rule.
                     if (plans(r)%recv_cls(kr) == FOLD_ROW_WEST) then
                        loc(plans(r)%recv_col(kr), drow) = val
                     else if (plans(r)%recv_cls(kr) == FOLD_ROW_SELF .and. neg) then
                        loc(plans(r)%recv_col(kr), drow) = 0.0_wp
                     end if
                  else
                     loc(plans(r)%recv_col(kr), drow) = val
                  end if
               end do
            end do
            deallocate (snd)
         end do
         ! Bitwise vs the serial fold over the WHOLE window (NaN fails).
         do j = 1, nyg
            do i = 1, size(loc, 1)
               if (.not. (loc(i, j) == ref(a - 1 + i, j))) then
                  call check(error, .false., "emulated fold /= serial fold: "// &
                             trim(what))
                  return
               end if
            end do
         end do
         deallocate (loc)
      end do
   end subroutine emulate

   pure real(wp) function encode(c, j) result(v)
      !! Index-encoded value, exact in wp, never zero.
      integer, intent(in) :: c, j
      v = real(1000*c + j, wp) + 0.25_wp
   end function encode

   pure integer function peer_index(plan, rx) result(p)
      type(fold_plan_t), intent(in) :: plan
      integer, intent(in) :: rx
      integer :: k
      p = 0
      do k = 1, plan%npeer
         if (plan%peer_rx(k) == rx) p = k
      end do
   end function peer_index

   pure logical function same_list(a, b) result(same)
      integer, intent(in) :: a(:), b(:)
      same = size(a) == size(b)
      if (same) same = all(a == b)
   end function same_list

   pure function recv_peers(plan, fam) result(peers)
      type(fold_plan_t), intent(in) :: plan
      integer, intent(in) :: fam
      integer, allocatable :: peers(:)
      integer :: k
      allocate (peers(0))
      do k = 1, plan%npeer
         if (plan%recv_n(k, fam) > 0) peers = [peers, plan%peer_rx(k)]
      end do
   end function recv_peers

end module test_ocean_fold_plan
