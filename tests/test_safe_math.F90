!! Unit tests for `rdb_safe_math`: the position-independent transcendental library.
!!
!! The module `#include`s the bodies itself (the form a production kernel module will
!! use) in addition to using the `rdb_safe_math` facade, so it can compare the two.
!!
!! What each group proves:
!!   * golden bits (contraction-OFF CPU builds only) -- a fixed table
!!     (`test_safe_math_golden.inc`) of arguments and the exact IEEE bit patterns of
!!     the results, asserted with NO tolerance, the same table on every toolchain.
!!     Produced by gfortran `-ffp-contract=off`, reproduced bit for bit by nvfortran
!!     `-Mnofma` and ifx `-no-fma`.  The GPU build contracts a*b+c and differs from
!!     the table by 1-2 ulp, so the table is NOT asserted there.
!!   * accuracy -- max ulp error over a sweep against 128-bit reals where the
!!     compiler has them (gfortran, ifx); on nvfortran / flang against the double
!!     intrinsic with extra slack.
!!   * special values -- NaN / Inf / zero / domain-error propagation.
!!   * position independence -- the included bodies evaluated in a `do concurrent`
!!     over contiguous arrays of several lengths and start offsets (elements fall in
!!     the vector body and in the remainder) are bit-identical to one-element
!!     launches of the same kernel (device scalar calls on the GPU build); on CPU
!!     builds also bit-identical to the facade's scalar calls.
!!
!! Kernels use explicit-shape dummies and `!$acc enter data` / `update`, inert on host.
module test_safe_math
   use, intrinsic :: iso_fortran_env, only: int64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite, ieee_value, &
                                                                               ieee_quiet_nan, ieee_positive_inf
   use testdrive, only: new_unittest, unittest_type, error_type, check, skip_test
   use rdb_constants, only: wp
   use rdb_safe_math, only: safe_exp, safe_log, safe_sin, safe_cos, safe_tan, safe_atan, &
                            safe_atan2, safe_tanh, safe_cosh, safe_acosh, safe_pow, &
                            safe_cbrt, safe_erfc
   implicit none
   private

   public :: collect_safe_math_tests

   integer, parameter :: NFUN = 13
   integer, parameter :: NSPEC = 8
   integer, parameter :: F_EXP = 1, F_LOG = 2, F_SIN = 3, F_COS = 4, F_TAN = 5, F_ATAN = 6, &
                         F_ATAN2 = 7, F_TANH = 8, F_COSH = 9, F_ACOSH = 10, F_POW = 11, &
                         F_CBRT = 12, F_ERFC = 13
   character(len=6), parameter :: FNAME(NFUN) = [character(len=6) :: 'exp', 'log', 'sin', 'cos', 'tan', &
                                                                    'atan', 'atan2', 'tanh', 'cosh', 'acosh', 'pow', 'cbrt', 'erfc']
   !! documented accuracy bound per function, in ulp
   real(wp), parameter :: ULP_BOUND(NFUN) = [1.3_wp, 3.0_wp, 3.0_wp, 3.0_wp, 5.0_wp, 2.5_wp, 3.0_wp, &
                                             4.0_wp, 2.5_wp, 5.0_wp, 2.0_wp, 1.5_wp, 7.0_wp]

#include "test_safe_math_golden.inc"

   ! the library bodies, included module-locally exactly as a kernel module would
#include "rdb_safe_math_consts.inc"

contains

#include "rdb_safe_math_core.inc"
#include "rdb_safe_math_pow.inc"
#include "rdb_safe_math_trig.inc"
#include "rdb_safe_math_atan.inc"
#include "rdb_safe_math_hyp.inc"
#include "rdb_safe_math_erfc.inc"

   subroutine collect_safe_math_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("safe_math_golden_exp", test_golden_exp), &
                  new_unittest("safe_math_golden_log", test_golden_log), &
                  new_unittest("safe_math_golden_sin", test_golden_sin), &
                  new_unittest("safe_math_golden_cos", test_golden_cos), &
                  new_unittest("safe_math_golden_tan", test_golden_tan), &
                  new_unittest("safe_math_golden_atan", test_golden_atan), &
                  new_unittest("safe_math_golden_atan2", test_golden_atan2), &
                  new_unittest("safe_math_golden_tanh", test_golden_tanh), &
                  new_unittest("safe_math_golden_cosh", test_golden_cosh), &
                  new_unittest("safe_math_golden_acosh", test_golden_acosh), &
                  new_unittest("safe_math_golden_pow", test_golden_pow), &
                  new_unittest("safe_math_golden_pow_large_exponent", test_golden_powbig), &
                  new_unittest("safe_math_golden_cbrt", test_golden_cbrt), &
                  new_unittest("safe_math_golden_erfc", test_golden_erfc), &
                  new_unittest("safe_math_golden_close", test_golden_close), &
                  new_unittest("safe_math_accuracy_ulp", test_accuracy), &
                  new_unittest("safe_math_exact_landmarks", test_landmarks), &
                  new_unittest("safe_math_special_values", test_special_values), &
                  new_unittest("safe_math_position_independent", test_position_independence) &
                  ]
   end subroutine collect_safe_math_tests

   ! -----------------------------------------------------------------
   ! scalar dispatcher
   ! -----------------------------------------------------------------

   function eval1(fid, x, y) result(r)
      !! Scalar call of function `fid` (`y` is the second argument of atan2 / pow).
      integer, intent(in) :: fid
      real(wp), intent(in) :: x, y
      real(wp) :: r
      select case (fid)
      case (F_EXP); r = safe_exp(x)
      case (F_LOG); r = safe_log(x)
      case (F_SIN); r = safe_sin(x)
      case (F_COS); r = safe_cos(x)
      case (F_TAN); r = safe_tan(x)
      case (F_ATAN); r = safe_atan(x)
      case (F_ATAN2); r = safe_atan2(y, x)
      case (F_TANH); r = safe_tanh(x)
      case (F_COSH); r = safe_cosh(x)
      case (F_ACOSH); r = safe_acosh(x)
      case (F_POW); r = safe_pow(x, y)
      case (F_CBRT); r = safe_cbrt(x)
      case default; r = safe_erfc(x)
      end select
   end function eval1

   function bits(x) result(b)
      real(wp), intent(in) :: x
      integer(int64) :: b
      b = transfer(x, 1_int64)
   end function bits

   function from_hex(h) result(x)
      !! Real64 whose IEEE bit pattern is the 16 hex digits `h`.
      character(len=16), intent(in) :: h
      real(wp) :: x
      integer(int64) :: b
      read (h, '(Z16)') b
      x = transfer(b, 1.0_wp)
   end function from_hex

   ! -----------------------------------------------------------------
   ! golden bits
   ! -----------------------------------------------------------------

   subroutine golden_check(fid, tab, error)
      !! Every row of `tab` (argument(s), result) reproduces to the last bit.
      integer, intent(in) :: fid
      character(len=16), intent(in) :: tab(:, :)
      type(error_type), allocatable, intent(out) :: error
      integer :: i, nbad, ncol, firstbad
      real(wp) :: x, y, want, got
      character(len=16) :: hw, hg
#ifdef RDB_GPU_OFFLOAD
      ! The GPU build contracts a*b+c (host and device): the table, produced by the
      ! contraction-off CPU builds, is not a contract here.  It is reported as SKIPPED
      ! (not silently passed); `safe_math_golden_close` still checks every row to 16 ulp,
      ! and position independence, accuracy and special values are asserted.
      if (size(tab) > 0) then
         call skip_test(error, "golden bits are a contraction-off CPU contract; GPU build contracts")
         return
      end if
#endif
      ncol = size(tab, 1)
      nbad = 0
      firstbad = 0
      hw = ''
      hg = ''
      do i = 1, size(tab, 2)
         x = from_hex(tab(1, i))
         y = 0.0_wp
         if (ncol == 3) y = from_hex(tab(2, i))
         want = from_hex(tab(ncol, i))
         got = eval1(fid, x, y)
         if (bits(got) /= bits(want)) then
            nbad = nbad + 1
            if (firstbad == 0) then
               firstbad = i
               write (hw, '(Z16.16)') bits(want)
               write (hg, '(Z16.16)') bits(got)
            end if
         end if
      end do
      call check(error, nbad, 0, trim(FNAME(fid))//": golden bit patterns differ; first row "// &
                 itoa(firstbad)//" want "//hw//" got "//hg)
   end subroutine golden_check

   function itoa(i) result(s)
      integer, intent(in) :: i
      character(len=12) :: s
      write (s, '(i0)') i
   end function itoa

   subroutine test_golden_exp(error)
      type(error_type), allocatable, intent(out) :: error
      call golden_check(F_EXP, G_EXP, error)
   end subroutine test_golden_exp

   subroutine test_golden_log(error)
      type(error_type), allocatable, intent(out) :: error
      call golden_check(F_LOG, G_LOG, error)
   end subroutine test_golden_log

   subroutine test_golden_sin(error)
      type(error_type), allocatable, intent(out) :: error
      call golden_check(F_SIN, G_SIN, error)
   end subroutine test_golden_sin

   subroutine test_golden_cos(error)
      type(error_type), allocatable, intent(out) :: error
      call golden_check(F_COS, G_COS, error)
   end subroutine test_golden_cos

   subroutine test_golden_tan(error)
      type(error_type), allocatable, intent(out) :: error
      call golden_check(F_TAN, G_TAN, error)
   end subroutine test_golden_tan

   subroutine test_golden_atan(error)
      type(error_type), allocatable, intent(out) :: error
      call golden_check(F_ATAN, G_ATAN, error)
   end subroutine test_golden_atan

   subroutine test_golden_atan2(error)
      type(error_type), allocatable, intent(out) :: error
      call golden_check(F_ATAN2, G_ATAN2, error)
   end subroutine test_golden_atan2

   subroutine test_golden_tanh(error)
      type(error_type), allocatable, intent(out) :: error
      call golden_check(F_TANH, G_TANH, error)
   end subroutine test_golden_tanh

   subroutine test_golden_cosh(error)
      type(error_type), allocatable, intent(out) :: error
      call golden_check(F_COSH, G_COSH, error)
   end subroutine test_golden_cosh

   subroutine test_golden_acosh(error)
      type(error_type), allocatable, intent(out) :: error
      call golden_check(F_ACOSH, G_ACOSH, error)
   end subroutine test_golden_acosh

   subroutine test_golden_pow(error)
      type(error_type), allocatable, intent(out) :: error
      call golden_check(F_POW, G_POW, error)
   end subroutine test_golden_pow

   subroutine test_golden_powbig(error)
      type(error_type), allocatable, intent(out) :: error
      call golden_check(F_POW, G_POWBIG, error)
   end subroutine test_golden_powbig

   subroutine test_golden_cbrt(error)
      type(error_type), allocatable, intent(out) :: error
      call golden_check(F_CBRT, G_CBRT, error)
   end subroutine test_golden_cbrt

   subroutine test_golden_erfc(error)
      type(error_type), allocatable, intent(out) :: error
      call golden_check(F_ERFC, G_ERFC, error)
   end subroutine test_golden_erfc

   function golden_far(fid, tab) result(nfar)
      !! Number of rows of `tab` whose result is more than 16 ulp from the golden value
      !! (any contraction setting stays within a few ulp; a zero golden value must be zero).
      integer, intent(in) :: fid
      character(len=16), intent(in) :: tab(:, :)
      integer :: nfar
      integer :: i, ncol
      real(wp) :: x, y, want, got
      ncol = size(tab, 1)
      nfar = 0
      do i = 1, size(tab, 2)
         x = from_hex(tab(1, i))
         y = 0.0_wp
         if (ncol == 3) y = from_hex(tab(2, i))
         want = from_hex(tab(ncol, i))
         got = eval1(fid, x, y)
         if (want == 0.0_wp) then
            if (got /= 0.0_wp) nfar = nfar + 1
         else if (.not. (abs(got - want) <= 16.0_wp*spacing(want))) then
            nfar = nfar + 1
         end if
      end do
   end function golden_far

   subroutine test_golden_close(error)
      !! Every golden row is reproduced to 16 ulp under ANY contraction setting (this is
      !! the golden-table check that stays on for the GPU build).
      type(error_type), allocatable, intent(out) :: error
      integer :: nfar
      nfar = golden_far(F_EXP, G_EXP) + golden_far(F_LOG, G_LOG) + golden_far(F_SIN, G_SIN) + &
             golden_far(F_COS, G_COS) + golden_far(F_TAN, G_TAN) + golden_far(F_ATAN, G_ATAN) + &
             golden_far(F_ATAN2, G_ATAN2) + golden_far(F_TANH, G_TANH) + golden_far(F_COSH, G_COSH) + &
             golden_far(F_ACOSH, G_ACOSH) + golden_far(F_POW, G_POW) + golden_far(F_POW, G_POWBIG) + &
             golden_far(F_CBRT, G_CBRT) + golden_far(F_ERFC, G_ERFC)
      call check(error, nfar, 0, "golden rows more than 16 ulp from the table")
   end subroutine test_golden_close

   ! -----------------------------------------------------------------
   ! accuracy
   ! -----------------------------------------------------------------

   function lcg(s) result(u)
      !! Deterministic uniform [0,1) from a 32-bit LCG state (pure integer arithmetic).
      integer(int64), intent(inout) :: s
      real(wp) :: u
      s = iand(s*1664525_int64 + 1013904223_int64, 4294967295_int64)
      u = real(s, wp)/4294967296.0_wp
   end function lcg

   subroutine sample_args(fid, k, s, x, y)
      !! The k-th sweep argument of function `fid` (valid range of the function).
      integer, intent(in) :: fid, k
      integer(int64), intent(inout) :: s
      real(wp), intent(out) :: x, y
      real(wp) :: t
      t = lcg(s)
      y = 0.0_wp
      select case (fid)
      case (F_EXP); x = -700.0_wp + 1400.0_wp*t
         if (mod(k, 3) == 0) x = -40.0_wp + 80.0_wp*t
      case (F_LOG); x = scale(1.0_wp + t, nint(60.0_wp*lcg(s)) - 30)
         if (mod(k, 3) == 0) x = 0.5_wp + 1.5_wp*t
      case (F_SIN, F_COS, F_TAN); x = -10.0_wp + 20.0_wp*t
         if (mod(k, 3) == 1) x = 1.6e6_wp*(2.0_wp*t - 1.0_wp)
         if (mod(k, 3) == 2) x = 4.5e5_wp*t
      case (F_ATAN); x = -10.0_wp + 20.0_wp*t
         if (mod(k, 2) == 0) then
            x = scale(1.0_wp + t, nint(60.0_wp*lcg(s)) - 30)
            if (lcg(s) < 0.5_wp) x = -x
         end if
      case (F_ATAN2); x = -10.0_wp + 20.0_wp*t
         y = -10.0_wp + 20.0_wp*lcg(s)
      case (F_TANH); x = -25.0_wp + 50.0_wp*t
         if (mod(k, 2) == 0) x = -0.5_wp + t
      case (F_COSH); x = -5.0_wp + 10.0_wp*t
         if (mod(k, 3) == 0) x = -710.0_wp + 1420.0_wp*t
      case (F_ACOSH); x = 1.0_wp + 3.0_wp*t
         if (mod(k, 2) == 0) x = 1.0_wp + scale(1.0_wp + t, nint(40.0_wp*lcg(s)) - 38)
      case (F_POW); x = scale(1.0_wp + t, nint(20.0_wp*lcg(s)) - 10)
         y = -4.0_wp + 8.0_wp*lcg(s)
      case (F_CBRT); x = scale(1.0_wp + t, nint(100.0_wp*lcg(s)) - 50)
         if (lcg(s) < 0.5_wp) x = -x
      case default; x = -6.0_wp + 36.0_wp*t
         if (mod(k, 3) == 0) x = -1.0_wp + 2.0_wp*t
      end select
   end subroutine sample_args

   function ulp_error(fid, x, y) result(e)
      !! |safe - exact| in units of the exact result's ulp.  The reference is a
      !! 128-bit evaluation where the compiler has one, else the double intrinsic
      !! (then the bound the caller applies is widened).
      integer, intent(in) :: fid
      real(wp), intent(in) :: x, y
      real(wp) :: e
      real(wp) :: got, rd, u
#if !defined(__NVCOMPILER) && !defined(__flang__)
      integer, parameter :: qp = selected_real_kind(33, 4931)
      real(qp) :: q, r
      q = real(x, qp)
      select case (fid)
      case (F_EXP); r = exp(q)
      case (F_LOG); r = log(q)
      case (F_SIN); r = sin(q)
      case (F_COS); r = cos(q)
      case (F_TAN); r = tan(q)
      case (F_ATAN); r = atan(q)
      case (F_ATAN2); r = atan2(real(y, qp), real(x, qp))
      case (F_TANH); r = tanh(q)
      case (F_COSH); r = cosh(q)
      case (F_ACOSH); r = acosh(q)
      case (F_POW); r = q**real(y, qp)
      case (F_CBRT); r = sign(abs(q)**(1.0_qp/3.0_qp), q)
      case default; r = erfc(q)
      end select
      got = eval1(fid, x, y)
      rd = real(r, wp)
      u = spacing(rd)
      e = real(abs((real(got, qp) - r)/real(u, qp)), wp)
#else
      select case (fid)
      case (F_EXP); rd = exp(x)
      case (F_LOG); rd = log(x)
      case (F_SIN); rd = sin(x)
      case (F_COS); rd = cos(x)
      case (F_TAN); rd = tan(x)
      case (F_ATAN); rd = atan(x)
      case (F_ATAN2); rd = atan2(y, x)
      case (F_TANH); rd = tanh(x)
      case (F_COSH); rd = cosh(x)
      case (F_ACOSH); rd = acosh(x)
      case (F_POW); rd = x**y
      case (F_CBRT); rd = x
      case default; rd = erfc(x)
      end select
      got = eval1(fid, x, y)
      if (fid == F_CBRT) then
         ! x**(1/3) in double is itself several ulp off for large |ln x|: judge the root by
         ! its cube instead (an error of e ulp in the root is 3e ulp in the cube; the
         ! cube's own rounding adds up to ~1.5 ulp, hence the slack the caller applies)
         e = abs(got*got*got - x)/(3.0_wp*spacing(x))
      else
         u = spacing(rd)
         e = abs(got - rd)/u
      end if
#endif
   end function ulp_error

   subroutine test_accuracy(error)
      !! Max ulp error of every function over 3000 sweep points stays within its
      !! documented bound (plus 3 ulp of slack when the reference is the double
      !! intrinsic).
      type(error_type), allocatable, intent(out) :: error
      integer :: fid, k
      integer(int64) :: s
      real(wp) :: x, y, e, emax, slack
#if defined(__NVCOMPILER) || defined(__flang__)
      ! No 128-bit real on this compiler (this includes the nvfortran GPU build): the
      ! reference is the double intrinsic, itself ~1 ulp, so the documented bound cannot
      ! be asserted tightly.  3 ulp of slack makes this a SANITY bound, labelled as such.
      slack = 3.0_wp
#else
      slack = 0.0_wp
#endif
      do fid = 1, NFUN
         s = 123456789_int64 + 7_int64*fid
         emax = 0.0_wp
         do k = 1, 3000
            call sample_args(fid, k, s, x, y)
            e = ulp_error(fid, x, y)
            if (ieee_is_nan(e)) then
               call check(error, .false., trim(FNAME(fid))//": NaN in the ulp error (NaN result for a valid argument)")
               return
            end if
            emax = max(emax, e)
         end do
         if (slack > 0.0_wp) then
            call check(error, emax <= ULP_BOUND(fid) + slack, trim(FNAME(fid))// &
                       ": max ulp error exceeds its bound (sanity bound, no 128-bit reference on this compiler)")
         else
            call check(error, emax <= ULP_BOUND(fid), trim(FNAME(fid))//": max ulp error exceeds its documented bound")
         end if
         if (allocated(error)) return
      end do
   end subroutine test_accuracy

   ! -----------------------------------------------------------------
   ! exact landmarks and special values
   ! -----------------------------------------------------------------

   subroutine test_landmarks(error)
      !! Results that are exactly representable come out exact.
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: Z = 0.0_wp, ONE = 1.0_wp
      call check(error, safe_exp(Z), ONE, "exp(0)")
      if (allocated(error)) return
      call check(error, safe_log(ONE), Z, "log(1)")
      if (allocated(error)) return
      call check(error, safe_sin(Z), Z, "sin(0)")
      if (allocated(error)) return
      call check(error, safe_cos(Z), ONE, "cos(0)")
      if (allocated(error)) return
      call check(error, safe_atan(Z), Z, "atan(0)")
      if (allocated(error)) return
      call check(error, safe_tanh(Z), Z, "tanh(0)")
      if (allocated(error)) return
      call check(error, safe_cosh(Z), ONE, "cosh(0)")
      if (allocated(error)) return
      call check(error, safe_acosh(ONE), Z, "acosh(1)")
      if (allocated(error)) return
      call check(error, safe_cbrt(27.0_wp), 3.0_wp, "cbrt(27)")
      if (allocated(error)) return
      call check(error, safe_cbrt(-8.0_wp), -2.0_wp, "cbrt(-8)")
      if (allocated(error)) return
      call check(error, safe_cbrt(0.125_wp), 0.5_wp, "cbrt(1/8)")
      if (allocated(error)) return
      call check(error, safe_pow(2.0_wp, 10.0_wp), 1024.0_wp, "2**10")
      if (allocated(error)) return
      call check(error, safe_pow(2.0_wp, -3.0_wp), 0.125_wp, "2**-3")
      if (allocated(error)) return
      call check(error, safe_pow(10.0_wp, 2.0_wp), 100.0_wp, "10**2")
      if (allocated(error)) return
      call check(error, safe_pow(4.0_wp, 0.5_wp), 2.0_wp, "4**0.5")
      if (allocated(error)) return
      call check(error, safe_pow(-2.0_wp, 3.0_wp), -8.0_wp, "(-2)**3")
      if (allocated(error)) return
      call check(error, safe_erfc(Z), ONE, "erfc(0)")
   end subroutine test_landmarks

   subroutine test_special_values(error)
      !! NaN / Inf / zero / domain errors propagate; nothing is laundered into a
      !! finite number.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: nan, inf, z, nz, one
      one = 1.0_wp
      z = 0.0_wp
      nz = -z
      nan = ieee_value(one, ieee_quiet_nan)
      inf = ieee_value(one, ieee_positive_inf)

      ! NaN in -> NaN out, every one-argument function
      call check(error, all([ieee_is_nan(safe_exp(nan)), ieee_is_nan(safe_log(nan)), ieee_is_nan(safe_sin(nan)), &
                             ieee_is_nan(safe_cos(nan)), ieee_is_nan(safe_tan(nan)), ieee_is_nan(safe_atan(nan)), &
                             ieee_is_nan(safe_tanh(nan)), ieee_is_nan(safe_cosh(nan)), ieee_is_nan(safe_acosh(nan)), &
                             ieee_is_nan(safe_cbrt(nan)), ieee_is_nan(safe_erfc(nan))]), "NaN in must give NaN out")
      if (allocated(error)) return
      call check(error, all([ieee_is_nan(safe_atan2(nan, one)), ieee_is_nan(safe_atan2(one, nan)), &
                             ieee_is_nan(safe_pow(nan, one)), ieee_is_nan(safe_pow(2.0_wp, nan))]), &
                 "NaN in atan2 / pow must give NaN out")
      if (allocated(error)) return

      ! exp / log
      call check(error, safe_exp(inf) > huge(one) .and. .not. ieee_is_finite(safe_exp(inf)), "exp(+Inf) = +Inf")
      if (allocated(error)) return
      call check(error, safe_exp(-inf), z, "exp(-Inf) = 0")
      if (allocated(error)) return
     call check(error,.not. ieee_is_finite(safe_exp(710.0_wp)) .and. safe_exp(710.0_wp) > z, "exp(710) overflows to +Inf, not huge")
      if (allocated(error)) return
      call check(error, safe_exp(-746.0_wp), z, "exp(-746) underflows to 0")
      if (allocated(error)) return
      call check(error,.not. ieee_is_finite(safe_log(z)) .and. safe_log(z) < z, "log(0) = -Inf")
      if (allocated(error)) return
      call check(error, ieee_is_nan(safe_log(-one)), "log(-1) = NaN")
      if (allocated(error)) return
      call check(error,.not. ieee_is_finite(safe_log(inf)) .and. safe_log(inf) > z, "log(+Inf) = +Inf")
      if (allocated(error)) return

      ! trig: beyond the valid range is NaN, never a wrong finite number
      call check(error, ieee_is_nan(safe_sin(2.0e6_wp)) .and. ieee_is_nan(safe_cos(-2.0e6_wp)) .and. &
                 ieee_is_nan(safe_tan(1.0e15_wp)), "trig beyond 1.6e6 rad is NaN")
      if (allocated(error)) return
      call check(error, ieee_is_nan(safe_sin(inf)) .and. ieee_is_nan(safe_cos(-inf)), "sin/cos(+-Inf) = NaN")
      if (allocated(error)) return
      call check(error, sign(1.0_wp, safe_sin(nz)) < z .and. sign(1.0_wp, safe_tan(nz)) < z .and. &
                 sign(1.0_wp, safe_atan(nz)) < z, "sin/tan/atan(-0) keep the sign of zero")
      if (allocated(error)) return

      ! atan / atan2
      call check(error, abs(safe_atan(inf) - 1.5707963267948966_wp) < 1.0e-15_wp .and. &
                 abs(safe_atan(-inf) + 1.5707963267948966_wp) < 1.0e-15_wp, "atan(+-Inf) = +-pi/2")
      if (allocated(error)) return
      call check(error, abs(safe_atan2(z, -one) - 3.141592653589793_wp) < 1.0e-15_wp .and. &
                 abs(safe_atan2(nz, -one) + 3.141592653589793_wp) < 1.0e-15_wp .and. safe_atan2(z, one) == z .and. &
                 sign(1.0_wp, safe_atan2(nz, one)) < z, "atan2 signed-zero table")
      if (allocated(error)) return
      call check(error, abs(safe_atan2(inf, inf) - 0.7853981633974483_wp) < 1.0e-15_wp .and. &
                 abs(safe_atan2(inf, -inf) - 2.356194490192345_wp) < 1.0e-15_wp .and. &
                 abs(safe_atan2(one, z) - 1.5707963267948966_wp) < 1.0e-15_wp, "atan2 infinities / x = 0")
      if (allocated(error)) return

      ! tanh / cosh / acosh / cbrt / erfc
      call check(error, safe_tanh(inf) == one .and. safe_tanh(-inf) == -one, "tanh(+-Inf) = +-1")
      if (allocated(error)) return
      call check(error,.not. ieee_is_finite(safe_cosh(inf)) .and. .not. ieee_is_finite(safe_cosh(711.0_wp)) .and. &
                 safe_cosh(-inf) > z, "cosh overflows to +Inf")
      if (allocated(error)) return
      call check(error, ieee_is_nan(safe_acosh(0.5_wp)) .and. ieee_is_nan(safe_acosh(-inf)) .and. &
                 .not. ieee_is_finite(safe_acosh(inf)), "acosh domain / +Inf")
      if (allocated(error)) return
      call check(error,.not. ieee_is_finite(safe_cbrt(inf)) .and. safe_cbrt(-inf) < z .and. safe_cbrt(z) == z .and. &
                 sign(1.0_wp, safe_cbrt(nz)) < z, "cbrt(+-Inf), cbrt(+-0)")
      if (allocated(error)) return
      call check(error, safe_erfc(inf) == z .and. safe_erfc(-inf) == 2.0_wp .and. safe_erfc(30.0_wp) == z .and. &
                 safe_erfc(-30.0_wp) == 2.0_wp, "erfc limits")
      if (allocated(error)) return

      ! pow: the IEEE-754 / C99 table
      call check(error, ieee_is_nan(safe_pow(nan, z)) .and. ieee_is_nan(safe_pow(one, nan)) .and. &
                 safe_pow(2.0_wp, z) == one .and. safe_pow(one, 5.0_wp) == one, &
                 "pow: NaN in any argument is NaN (unlike C99); pow(x,0) = pow(1,y) = 1 otherwise")
      if (allocated(error)) return
      call check(error, safe_pow(z, 2.0_wp) == z .and. .not. ieee_is_finite(safe_pow(z, -2.0_wp)) .and. &
                 sign(1.0_wp, safe_pow(nz, 3.0_wp)) < z .and. sign(1.0_wp, safe_pow(nz, 2.0_wp)) > z, "pow(+-0, y)")
      if (allocated(error)) return
      call check(error, ieee_is_nan(safe_pow(-8.0_wp, 1.0_wp/3.0_wp)) .and. safe_pow(-2.0_wp, 2.0_wp) == 4.0_wp .and. &
                 safe_pow(-2.0_wp, 3.0_wp) == -8.0_wp, "pow(negative base)")
      if (allocated(error)) return
      call check(error,.not. ieee_is_finite(safe_pow(2.0_wp, inf)) .and. safe_pow(0.5_wp, inf) == z .and. &
                 safe_pow(2.0_wp, -inf) == z .and. safe_pow(-one, inf) == one, "pow(x, +-Inf)")
      if (allocated(error)) return
      call check(error,.not. ieee_is_finite(safe_pow(inf, 2.0_wp)) .and. safe_pow(inf, -2.0_wp) == z .and. &
                 safe_pow(-inf, 3.0_wp) < z .and. .not. ieee_is_finite(safe_pow(-inf, 3.0_wp)), "pow(+-Inf, y)")
      if (allocated(error)) return
      call check(error,.not. ieee_is_finite(safe_pow(2.0_wp, 2000.0_wp)) .and. safe_pow(2.0_wp, -2000.0_wp) == z, &
                 "pow overflow -> Inf, underflow -> 0")
   end subroutine test_special_values

   ! -----------------------------------------------------------------
   ! position independence
   ! -----------------------------------------------------------------

   function special_args() result(v)
      !! The special arguments planted in elements 1..NSPEC of the kernel arrays.
      real(wp) :: v(NSPEC)
      real(wp) :: zero
      zero = 0.0_wp
      v = [ieee_value(zero, ieee_quiet_nan), ieee_value(zero, ieee_positive_inf), &
           -ieee_value(zero, ieee_positive_inf), zero, -zero, 1.0e30_wp, transfer(2251799813685248_int64, zero), -1.0_wp]
   end function special_args

   function special_violations(fid, r) result(nbad)
      !! Number of elements 1..NSPEC of the kernel result `r` that break the expected
      !! class for function `fid`.  Codes: N NaN, P +Inf, M -Inf, Z zero (either sign),
      !! 1 / 2 / m the values 1 / 2 / -1, F finite, - not checked.  Columns: x = NaN, +Inf,
      !! -Inf, +0, -0, 1e30, subnormal (2**-1023), -1 (atan2 / pow: only the NaN column is exact).
      integer, intent(in) :: fid
      real(wp), intent(in) :: r(*)
      integer :: nbad
      character(len=NSPEC), parameter :: EXPECT(NFUN) = [ &
                                         "NPZ11P1F", "NPNMMF-N", "NNNZZNZF", "NNN11N1F", "NNNZZNZF", "NFFZZFZF", &
                                         "N-------", "N1mZZ1ZF", "NPP11P1F", "NPNNNF-N", "N-------", "NPMZZFZm", &
                                         "NZ211Z1F"]
      integer :: j
      character :: c
      logical :: ok
      nbad = 0
      do j = 1, NSPEC
         c = EXPECT(fid) (j:j)
         select case (c)
         case ("N"); ok = ieee_is_nan(r(j))
         case ("P"); ok = .not. ieee_is_finite(r(j)) .and. .not. ieee_is_nan(r(j)) .and. r(j) > 0.0_wp
         case ("M"); ok = .not. ieee_is_finite(r(j)) .and. .not. ieee_is_nan(r(j)) .and. r(j) < 0.0_wp
         case ("Z"); ok = r(j) == 0.0_wp
         case ("1"); ok = r(j) == 1.0_wp
         case ("2"); ok = r(j) == 2.0_wp
         case ("m"); ok = r(j) == -1.0_wp
         case ("F"); ok = ieee_is_finite(r(j))
         case default; ok = .true.
         end select
         if (.not. ok) nbad = nbad + 1
      end do
   end function special_violations

   subroutine kern(fid, n, off, x, y2, r)
      !! Apply function `fid` to elements off+1 .. off+n inside `do concurrent`.
      !! One small loop per function, as a production kernel would be.
      integer, intent(in) :: fid, n, off
      real(wp), intent(in) :: x(off + n), y2(off + n)
      real(wp), intent(inout) :: r(off + n)
      integer :: i
      select case (fid)
      case (F_EXP)
         do concurrent(i=off + 1:off + n)
            r(i) = rdb_sm_exp(x(i))
         end do
      case (F_LOG)
         do concurrent(i=off + 1:off + n)
            r(i) = rdb_sm_log(x(i))
         end do
      case (F_SIN)
         do concurrent(i=off + 1:off + n)
            r(i) = rdb_sm_sin(x(i))
         end do
      case (F_COS)
         do concurrent(i=off + 1:off + n)
            r(i) = rdb_sm_cos(x(i))
         end do
      case (F_TAN)
         do concurrent(i=off + 1:off + n)
            r(i) = rdb_sm_tan(x(i))
         end do
      case (F_ATAN)
         do concurrent(i=off + 1:off + n)
            r(i) = rdb_sm_atan(x(i))
         end do
      case (F_ATAN2)
         do concurrent(i=off + 1:off + n)
            r(i) = rdb_sm_atan2(y2(i), x(i))
         end do
      case (F_TANH)
         do concurrent(i=off + 1:off + n)
            r(i) = rdb_sm_tanh(x(i))
         end do
      case (F_COSH)
         do concurrent(i=off + 1:off + n)
            r(i) = rdb_sm_cosh(x(i))
         end do
      case (F_ACOSH)
         do concurrent(i=off + 1:off + n)
            r(i) = rdb_sm_acosh(x(i))
         end do
      case (F_POW)
         do concurrent(i=off + 1:off + n)
            r(i) = rdb_sm_pow(x(i), y2(i))
         end do
      case (F_CBRT)
         do concurrent(i=off + 1:off + n)
            r(i) = rdb_sm_cbrt(x(i))
         end do
      case default
         do concurrent(i=off + 1:off + n)
            r(i) = rdb_sm_erfc(x(i))
         end do
      end select
   end subroutine kern

   subroutine test_position_independence(error)
      !! Vector-body / remainder / start-offset independence of the INCLUDED bodies: the
      !! array result of every function equals one-element launches of the same kernel
      !! (device scalar calls on the GPU build), for several lengths and offsets; on CPU
      !! builds it also equals the facade's scalar calls (separate compilation unit).
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NMAX = 200, NOFF = 8
      integer, parameter :: LENS(10) = [1, 2, 3, 4, 5, 7, 8, 9, 31, 128]
      real(wp) :: x(NMAX + NOFF), y2(NMAX + NOFF), r(NMAX + NOFF), r1(NMAX + NOFF)
      integer :: fid, il, off, n, i, nbad, nbad_fac, k, nspec
      integer(int64) :: s
      nspec = 0
      nbad = 0
      nbad_fac = 0
      do fid = 1, NFUN
         s = 987654321_int64 + 13_int64*fid
         do k = 1, NMAX + NOFF
            call sample_args(fid, k, s, x(k), y2(k))
         end do
         ! elements 1..8: NaN, +Inf, -Inf, +0, -0, out of range, subnormal, -1
         x(1:NSPEC) = special_args()
         r = 0.0_wp
         r1 = 0.0_wp
         !$acc enter data copyin(x, y2, r, r1)
         ! one-element launches: the scalar reference
         do i = 1, NMAX + NOFF
            call kern(fid, 1, i - 1, x, y2, r1)
         end do
         !$acc update self(r1)
         nspec = nspec + special_violations(fid, r1)
         do il = 1, size(LENS)
            n = LENS(il)
            do off = 0, NOFF
               call kern(fid, n, off, x, y2, r)
               !$acc update self(r)
               do i = off + 1, off + n
                  if (bits(r(i)) /= bits(r1(i))) nbad = nbad + 1
#ifndef RDB_GPU_OFFLOAD
                  if (bits(r(i)) /= bits(eval1(fid, x(i), y2(i)))) nbad_fac = nbad_fac + 1
#endif
               end do
            end do
         end do
         !$acc exit data delete(x, y2, r, r1)
      end do
      call check(error, nbad, 0, "array (do concurrent) results differ from one-element launches")
      if (allocated(error)) return
      call check(error, nspec, 0, "special values (NaN, Inf, 0, out of range, subnormal) wrong in the kernel results")
      if (allocated(error)) return
      call check(error, nbad_fac, 0, "included copy differs from the rdb_safe_math facade")
   end subroutine test_position_independence

end module test_safe_math
