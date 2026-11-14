!! Tests of the tracer accessor pair `ice_trcr_get/set` and of the in-place
!! aggregate `itd_aggregate_cell` against the packed `itd_aggregate`.
!!
!! The in-place form runs the same additions in the same order as the packed
!! form, so on CPU builds the two agree BITWISE (checked with `==` on finite
!! results and with `ieee_is_nan` on the poisoned cell).  On the offload build
!! FMA contraction applies to host code too, so the comparison is a scale-aware
!! 1e-12 there (measured worst difference is quoted in the failure message).
!! A `do concurrent` run of the in-place kernel on mapped arrays is compared
!! to the host packed aggregate with the same tolerance rule.
module test_cice_agg_cell
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
   use, intrinsic :: iso_fortran_env, only: int32, int64
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use rdb_constants, only: wp
   use rdb_ice_dims, only: ICE_NTRCR_MAX
   use rdb_ice_errors, only: ICE_ERR_NONFINITE
   use rdb_ice_params, only: ice_ntrcr, ice_nt_Tsfc, ice_nt_qice, ice_nt_qsno, ice_nt_sice, &
                             ice_trcr_get, ice_trcr_set, ice_trcr_set_mean, ICE_QNAN
   use rdb_ice_itd_shift, only: itd_aggregate, itd_aggregate_cell
   use lcg_deterministic, only: lcg_next, lcg_unit
   implicit none
   private

   public :: collect_cice_agg_cell_tests

   real(wp), parameter :: PUNY = 1.0e-11_wp
   real(wp), parameter :: TF = -1.8_wp

contains

   subroutine collect_cice_agg_cell_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("qnan_constant_is_nan", test_qnan), &
                  new_unittest("trcr_accessor_roundtrip", test_roundtrip), &
                  new_unittest("trcr_set_mean_matches_set", test_set_mean), &
                  new_unittest("aggregate_cell_vs_packed_nilyr7", test_nilyr7), &
                  new_unittest("aggregate_cell_vs_packed_nilyr4", test_nilyr4), &
                  new_unittest("aggregate_cell_device", test_device) &
                  ]
   end subroutine collect_cice_agg_cell_tests

   subroutine fill_state(nx, ny, ncat, nilyr, nslyr, seed, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon)
      !! Random state: ~30 % empty categories, some snow-free, tracers in plausible ranges.
      integer, intent(in) :: nx, ny, ncat, nilyr, nslyr
      integer(int64), intent(in) :: seed
      real(wp), intent(out) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat)
      real(wp), intent(out) :: Tsfcn(nx, ny, ncat), qicen(nx, ny, ncat, nilyr)
      real(wp), intent(out) :: sicen(nx, ny, ncat, nilyr), qsnon(nx, ny, ncat, nslyr)
      integer :: i, j, n, k
      integer(int64) :: st
      real(wp) :: u

      st = lcg_next(seed)
      do n = 1, ncat
         do j = 1, ny
            do i = 1, nx
               st = lcg_next(st); u = lcg_unit(st)
               if (u < 0.3_wp) then
                  aicen(i, j, n) = 0.0_wp
                  vicen(i, j, n) = 0.0_wp
                  vsnon(i, j, n) = 0.0_wp
               else
                  aicen(i, j, n) = 0.05_wp + 0.15_wp*lcg_unit(lcg_next(st))
                  st = lcg_next(st)
                  vicen(i, j, n) = aicen(i, j, n)*(0.1_wp + 3.0_wp*lcg_unit(lcg_next(st)))
                  st = lcg_next(st)
                  vsnon(i, j, n) = aicen(i, j, n)*0.4_wp*lcg_unit(lcg_next(st))
                  st = lcg_next(st)
               end if
               st = lcg_next(st); Tsfcn(i, j, n) = -30.0_wp*lcg_unit(st)
               do k = 1, nilyr
                  st = lcg_next(st); qicen(i, j, n, k) = -3.0e8_wp + 1.0e8_wp*lcg_unit(st)
                  st = lcg_next(st); sicen(i, j, n, k) = 0.1_wp + 5.9_wp*lcg_unit(st)
               end do
               do k = 1, nslyr
                  st = lcg_next(st); qsnon(i, j, n, k) = -1.5e8_wp + 0.3e8_wp*lcg_unit(st)
               end do
            end do
         end do
      end do
   end subroutine fill_state

   subroutine test_roundtrip(error)
      !! `ice_trcr_get` returns exactly what `ice_trcr_set` stored, at the index the
      !! packed model assigns (Tsfc, qice, qsno, sice), and NaN / no-op outside it.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nx = 3, ny = 2, ncat = 2, nilyr = 7, nslyr = 1
      real(wp) :: T(nx, ny, ncat), qi(nx, ny, ncat, nilyr), si(nx, ny, ncat, nilyr), qs(nx, ny, ncat, nslyr)
      real(wp) :: v
      integer :: it, ntr

      T = 0.0_wp; qi = 0.0_wp; si = 0.0_wp; qs = 0.0_wp
      ntr = ice_ntrcr(nilyr, nslyr)
      do it = 1, ntr
         call ice_trcr_set(it, 2, 1, 2, nx, ny, ncat, nilyr, nslyr, 100.0_wp + real(it, wp), T, qi, si, qs)
      end do
      call check(error, T(2, 1, 2) == 100.0_wp + real(ice_nt_Tsfc(), wp), "Tsfc slot")
      if (allocated(error)) return
      call check(error, qi(2, 1, 2, nilyr) == 100.0_wp + real(ice_nt_qice() + nilyr - 1, wp), "last qice slot")
      if (allocated(error)) return
      call check(error, qs(2, 1, 2, 1) == 100.0_wp + real(ice_nt_qsno(nilyr), wp), "qsno slot")
      if (allocated(error)) return
      call check(error, si(2, 1, 2, 1) == 100.0_wp + real(ice_nt_sice(nilyr, nslyr), wp), "first sice slot")
      if (allocated(error)) return
      call check(error, si(2, 1, 2, nilyr) == 100.0_wp + real(ntr, wp), "last sice slot")
      if (allocated(error)) return
      ! nothing else was touched
      call check(error, count(T /= 0.0_wp) + count(qi /= 0.0_wp) + count(si /= 0.0_wp) + count(qs /= 0.0_wp) == ntr, &
                 "exactly ntrcr entries written")
      if (allocated(error)) return
      do it = 1, ntr
         v = ice_trcr_get(it, 2, 1, 2, nx, ny, ncat, nilyr, nslyr, T, qi, si, qs)
         if (.not. (v == 100.0_wp + real(it, wp))) then
            call test_failed(error, "get after set")
            return
         end if
      end do
      v = ice_trcr_get(0, 2, 1, 2, nx, ny, ncat, nilyr, nslyr, T, qi, si, qs)
      call check(error, ieee_is_nan(v), "index 0 gives NaN")
      if (allocated(error)) return
      v = ice_trcr_get(ntr + 1, 2, 1, 2, nx, ny, ncat, nilyr, nslyr, T, qi, si, qs)
      call check(error, ieee_is_nan(v), "index ntrcr+1 gives NaN")
      if (allocated(error)) return
      call ice_trcr_set(ntr + 1, 1, 1, 1, nx, ny, ncat, nilyr, nslyr, 5.0_wp, T, qi, si, qs)
      call check(error, T(1, 1, 1) == 0.0_wp .and. si(1, 1, 1, 1) == 0.0_wp, "out-of-range set is a no-op")
   end subroutine test_roundtrip

   subroutine test_qnan(error)
      !! The device-safe NaN constant really is a NaN in the working precision (a
      !! 64-bit pattern `transfer`red into a 32-bit real reads as 0.0, not NaN).
      type(error_type), allocatable, intent(out) :: error

      call check(error, ieee_is_nan(ICE_QNAN), "ICE_QNAN is NaN")
   end subroutine test_qnan

   subroutine test_set_mean(error)
      !! `ice_trcr_set_mean` writes the same slots as `ice_trcr_set` with ncat = 1.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nx = 3, ny = 2, nilyr = 7, nslyr = 1
      real(wp) :: T1(nx, ny, 1), q1(nx, ny, 1, nilyr), s1(nx, ny, 1, nilyr), qs1(nx, ny, 1, nslyr)
      real(wp) :: T2(nx, ny), q2(nx, ny, nilyr), s2(nx, ny, nilyr), qs2(nx, ny, nslyr)
      integer :: it, ntr

      T1 = 0.0_wp; q1 = 0.0_wp; s1 = 0.0_wp; qs1 = 0.0_wp
      T2 = 0.0_wp; q2 = 0.0_wp; s2 = 0.0_wp; qs2 = 0.0_wp
      ntr = ice_ntrcr(nilyr, nslyr)
      do it = 0, ntr + 1                                   ! 0 and ntr+1 are ignored by design
         call ice_trcr_set(it, 2, 1, 1, nx, ny, 1, nilyr, nslyr, 10.0_wp + real(it, wp), T1, q1, s1, qs1)
         call ice_trcr_set_mean(it, 2, 1, nx, ny, nilyr, nslyr, 10.0_wp + real(it, wp), T2, q2, qs2, s2)
      end do
      call check(error, all(T1(:, :, 1) == T2) .and. all(q1(:, :, 1, :) == q2) .and. all(s1(:, :, 1, :) == s2) &
                 .and. all(qs1(:, :, 1, :) == qs2), "set_mean fills the same slots as set")
      if (allocated(error)) return
      call check(error, count(T2 /= 0.0_wp) + count(q2 /= 0.0_wp) + count(s2 /= 0.0_wp) + count(qs2 /= 0.0_wp) == ntr, &
                 "exactly ntrcr entries written")
   end subroutine test_set_mean

   pure function mean_get(it, i, j, nx, ny, nilyr, nslyr, Tc, qc, sc, qsc) result(v)
      !! Cell-mean tracer `it` (packed order) read back from the cell-mean arrays.
      integer, intent(in) :: it, i, j, nx, ny, nilyr, nslyr
      real(wp), intent(in) :: Tc(nx, ny), qc(nx, ny, nilyr), sc(nx, ny, nilyr), qsc(nx, ny, nslyr)
      real(wp) :: v

      if (it == ice_nt_Tsfc()) then
         v = Tc(i, j)
      else if (it < ice_nt_qsno(nilyr)) then
         v = qc(i, j, it - ice_nt_qice() + 1)
      else if (it < ice_nt_sice(nilyr, nslyr)) then
         v = qsc(i, j, it - ice_nt_qsno(nilyr) + 1)
      else
         v = sc(i, j, it - ice_nt_sice(nilyr, nslyr) + 1)
      end if
   end function mean_get

   subroutine test_nilyr7(error)
      type(error_type), allocatable, intent(out) :: error

      call run_compare(error, 7, 1)
   end subroutine test_nilyr7

   subroutine test_nilyr4(error)
      type(error_type), allocatable, intent(out) :: error

      call run_compare(error, 4, 1)
   end subroutine test_nilyr4

   subroutine host_packed(i, j, nx, ny, ncat, nilyr, nslyr, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, &
                          aice, vice, vsno, aice0, trcr, ierr)
      !! Reference: gather the column packed and call the packed aggregate.
      integer, intent(in) :: i, j, nx, ny, ncat, nilyr, nslyr
      real(wp), intent(in) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat)
      real(wp), intent(in) :: Tsfcn(nx, ny, ncat), qicen(nx, ny, ncat, nilyr)
      real(wp), intent(in) :: sicen(nx, ny, ncat, nilyr), qsnon(nx, ny, ncat, nslyr)
      real(wp), intent(out) :: aice, vice, vsno, aice0, trcr(ICE_NTRCR_MAX)
      integer(int32), intent(inout) :: ierr
      real(wp) :: a(ncat), v(ncat), s(ncat), tr(ICE_NTRCR_MAX, ncat)
      integer :: n, it

      tr = 0.0_wp
      do n = 1, ncat
         a(n) = aicen(i, j, n); v(n) = vicen(i, j, n); s(n) = vsnon(i, j, n)
         do it = 1, ice_ntrcr(nilyr, nslyr)
            tr(it, n) = ice_trcr_get(it, i, j, n, nx, ny, ncat, nilyr, nslyr, Tsfcn, qicen, sicen, qsnon)
         end do
      end do
      call itd_aggregate(ncat, ICE_NTRCR_MAX, nilyr, nslyr, a, v, s, tr, TF, PUNY, aice, vice, vsno, aice0, trcr, ierr)
   end subroutine host_packed

   subroutine unpack_cell(i, j, nx, ny, nilyr, nslyr, Tc, qc, sc, qsc, trcr)
      !! Cell-mean outputs of the in-place form in packed order.
      integer, intent(in) :: i, j, nx, ny, nilyr, nslyr
      real(wp), intent(in) :: Tc(nx, ny), qc(nx, ny, nilyr), sc(nx, ny, nilyr), qsc(nx, ny, nslyr)
      real(wp), intent(out) :: trcr(ICE_NTRCR_MAX)
      integer :: it

      trcr = 0.0_wp
      do it = 1, ice_ntrcr(nilyr, nslyr)
         trcr(it) = mean_get(it, i, j, nx, ny, nilyr, nslyr, Tc, qc, sc, qsc)
      end do
   end subroutine unpack_cell

   logical function same(a, b, tol)
      !! Bitwise (tol = 0) or scale-aware agreement; NaN matches NaN only.
      real(wp), intent(in) :: a, b, tol

      if (ieee_is_nan(a) .or. ieee_is_nan(b)) then
         same = ieee_is_nan(a) .and. ieee_is_nan(b)
      else if (tol == 0.0_wp) then
         same = (a == b)
      else
         same = abs(a - b) <= tol*max(abs(a), abs(b), 1.0_wp)
      end if
   end function same

   subroutine run_compare(error, nilyr, nslyr)
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: nilyr, nslyr
      integer, parameter :: nx = 24, ny = 20, ncat = 5
      real(wp) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat), Tsfcn(nx, ny, ncat)
      real(wp) :: qicen(nx, ny, ncat, nilyr), sicen(nx, ny, ncat, nilyr), qsnon(nx, ny, ncat, nslyr)
      real(wp) :: Tc(nx, ny), qc(nx, ny, nilyr), sc(nx, ny, nilyr), qsc(nx, ny, nslyr)
      real(wp) :: aice, vice, vsno, aice0, ref(4), got(4), trp(ICE_NTRCR_MAX), trc(ICE_NTRCR_MAX), tol
      real(wp) :: worst
      integer(int32) :: ierr_p, ierr_c(nx, ny)
      integer :: i, j, it, nbad_err

#ifdef RDB_GPU_OFFLOAD
      tol = 1.0e-12_wp
#else
      tol = 0.0_wp
#endif
      call fill_state(nx, ny, ncat, nilyr, nslyr, 987654321_int64 + int(nilyr, int64), &
                      aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon)
      ! one cell with a non-finite area: both forms must poison and flag
      aicen(3, 4, 2) = transfer(int(z'7FF8000000000000', int64), 1.0_wp)
      ! one all-empty cell
      aicen(5, 5, :) = 0.0_wp; vicen(5, 5, :) = 0.0_wp; vsnon(5, 5, :) = 0.0_wp
      ierr_c = 0_int32
      Tc = 0.0_wp; qc = 0.0_wp; sc = 0.0_wp; qsc = 0.0_wp
      worst = 0.0_wp
      nbad_err = 0
      do j = 1, ny
         do i = 1, nx
            ierr_p = 0_int32
            call host_packed(i, j, nx, ny, ncat, nilyr, nslyr, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, &
                             ref(1), ref(2), ref(3), ref(4), trp, ierr_p)
            call itd_aggregate_cell(i, j, nx, ny, ncat, nilyr, nslyr, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, &
                                    TF, PUNY, got(1), got(2), got(3), got(4), Tc, qc, sc, qsc, ierr_c)
            call unpack_cell(i, j, nx, ny, nilyr, nslyr, Tc, qc, sc, qsc, trc)
            do it = 1, 4
               if (.not. same(ref(it), got(it), tol)) then
                  call test_failed(error, "cell scalar aggregate differs from the packed form")
                  return
               end if
            end do
            do it = 1, ice_ntrcr(nilyr, nslyr)
               if (.not. same(trp(it), trc(it), tol)) then
                  call test_failed(error, "cell-mean tracer differs from the packed form")
                  return
               end if
               if (ieee_is_finite(trp(it)) .and. ieee_is_finite(trc(it))) &
                  worst = max(worst, abs(trp(it) - trc(it))/max(abs(trp(it)), 1.0_wp))
            end do
            if ((ierr_p /= 0_int32) .neqv. (ierr_c(i, j) /= 0_int32)) nbad_err = nbad_err + 1
            if (ierr_p /= ierr_c(i, j)) nbad_err = nbad_err + 1
         end do
      end do
      call check(error, nbad_err, 0, "error words agree with the packed form")
      if (allocated(error)) return
      call check(error, iand(ierr_c(3, 4), ICE_ERR_NONFINITE) /= 0_int32, "NaN cell flagged")
      if (allocated(error)) return
      call check(error, count(ierr_c /= 0_int32), 1, "only the NaN cell is flagged")
      if (allocated(error)) return
      call check(error, worst <= tol, "worst relative difference within tolerance")
   end subroutine run_compare

   subroutine test_device(error)
      !! The in-place kernel inside one `do concurrent` over mapped arrays
      !! against the host packed aggregate.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nx = 24, ny = 20, ncat = 5, nilyr = 7, nslyr = 1
      real(wp) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat), Tsfcn(nx, ny, ncat)
      real(wp) :: qicen(nx, ny, ncat, nilyr), sicen(nx, ny, ncat, nilyr), qsnon(nx, ny, ncat, nslyr)
      real(wp) :: Tc(nx, ny), qc(nx, ny, nilyr), sc(nx, ny, nilyr), qsc(nx, ny, nslyr)
      real(wp) :: ai(nx, ny), vi(nx, ny), vs(nx, ny), a0(nx, ny)
      real(wp) :: ref(4), trp(ICE_NTRCR_MAX), trc(ICE_NTRCR_MAX), tol
      integer(int32) :: ierr_p, ierr_d(nx, ny)
      integer :: i, j, it

#ifdef RDB_GPU_OFFLOAD
      tol = 1.0e-12_wp
#else
      tol = 0.0_wp
#endif
      call fill_state(nx, ny, ncat, nilyr, nslyr, 13579_int64, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon)
      ierr_d = 0_int32
      ai = 0.0_wp; vi = 0.0_wp; vs = 0.0_wp; a0 = 0.0_wp
      Tc = 0.0_wp; qc = 0.0_wp; sc = 0.0_wp; qsc = 0.0_wp
      !$acc enter data copyin(aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, ierr_d) create(ai, vi, vs, a0, Tc, qc, sc, qsc)
      do concurrent(j=1:ny, i=1:nx)
         call itd_aggregate_cell(i, j, nx, ny, ncat, nilyr, nslyr, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, &
                                 TF, PUNY, ai(i, j), vi(i, j), vs(i, j), a0(i, j), Tc, qc, sc, qsc, ierr_d)
      end do
      !$acc update self(ai, vi, vs, a0, Tc, qc, sc, qsc, ierr_d)
      !$acc exit data delete(aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, ierr_d, ai, vi, vs, a0, Tc, qc, sc, qsc)
      do j = 1, ny
         do i = 1, nx
            ierr_p = 0_int32
            call host_packed(i, j, nx, ny, ncat, nilyr, nslyr, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, &
                             ref(1), ref(2), ref(3), ref(4), trp, ierr_p)
            if (.not. (same(ref(1), ai(i, j), tol) .and. same(ref(2), vi(i, j), tol) .and. &
                       same(ref(3), vs(i, j), tol) .and. same(ref(4), a0(i, j), tol))) then
               call test_failed(error, "device scalar aggregate differs from host packed")
               return
            end if
            call unpack_cell(i, j, nx, ny, nilyr, nslyr, Tc, qc, sc, qsc, trc)
            do it = 1, ice_ntrcr(nilyr, nslyr)
               if (.not. same(trp(it), trc(it), tol)) then
                  call test_failed(error, "device cell-mean tracer differs from host packed")
                  return
               end if
            end do
         end do
      end do
      call check(error, count(ierr_d /= 0_int32), 0, "device error words clean")
   end subroutine test_device

end module test_cice_agg_cell
