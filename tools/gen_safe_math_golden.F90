!! Prints tests/test_safe_math_golden.inc: the golden bit-pattern table of
!! tests/test_safe_math.F90.
!!
!! Build and run it on a CONTRACTION-OFF toolchain (it is a stand-alone driver, not part
!! of the CMake build), with the math library object built the same way:
!!
!!   gfortran -cpp -O3 -ffp-contract=off -ffree-line-length-none -std=gnu -Isrc/shared_module_utilities  &
!!       <stub rdb_constants: integer, parameter :: wp = selected_real_kind(15,307)>  &
!!       src/framework/rdb_safe_math.F90 tools/gen_safe_math_golden.F90 && ./a.out > tests/test_safe_math_golden.inc
!!
!! Regenerate ONLY for a deliberate change of an algorithm or table (then say so in the
!! commit).  The arguments come from an integer xorshift generator and exact power-of-two
!! scalings -- no library transcendental -- so every compiler draws the same arguments.
!! Rows avoid subnormal arguments; atan2 and pow rows with subnormal results cover the flush rule; exp rows in [-745,-708] cover the flush-to-zero region.
program gen_safe_math_golden
   use, intrinsic :: iso_fortran_env, only: real64, int64
   use rdb_safe_math, only: safe_exp, safe_log, safe_sin, safe_cos, safe_tan, safe_atan, safe_atan2, &
                            safe_tanh, safe_cosh, safe_acosh, safe_pow, safe_cbrt, safe_erfc
   implicit none
   integer, parameter :: wp = real64
   integer, parameter :: NF = 14
   integer, parameter :: NG = 48
   character(len=6), parameter :: FN(NF) = [character(len=6) :: "EXP", "LOG", "SIN", "COS", "TAN", "ATAN", "ATAN2", &
                                                                           "TANH", "COSH", "ACOSH", "POW", "CBRT", "ERFC", "POWBIG"]
   integer(int64) :: s
   integer :: f, i, ncol
   write (*, "(a)") "! test_safe_math_golden.inc -- GENERATED golden bit patterns for tests/test_safe_math.F90."
   write (*, "(a)") "! Each row is the argument(s) and the result as 16 hex digits (IEEE binary64 bits); atan2 rows"
   write (*, "(a)") "! are (x, y, atan2(y, x)), pow rows (x, y, x**y).  Produced by the contraction-off builds"
   write (*, "(a)") "! (gfortran -ffp-contract=off, nvfortran -Mnofma), which reproduce it bit for bit.  Rows avoid"
   write (*, "(a)") "! subnormal arguments; the library flushes results below tiny to exactly 0 on every compiler."
   write (*, "(a)") "! Regenerate only for a deliberate change: tools/gen_safe_math_golden.F90."
   do f = 1, NF
      s = 88172645463325252_int64 + 7919_int64*f
      ncol = 2
      if (f == 7 .or. f == 11 .or. f == 14) ncol = 3
      write (*, "(a,a,a,i0)") "   integer, parameter :: NG_", trim(FN(f)), " = ", NG
      write (*, "(a,a,a,i0,a,a,a)") "   character(len=16), parameter :: G_", trim(FN(f)), "(", ncol, ", NG_", trim(FN(f)), &
         ") = reshape([character(len=16) :: &"
      do i = 1, NG
         call row(f, i, ncol, i == NG)
      end do
   end do
contains
   function rnd() result(r)
      real(wp) :: r
      s = ieor(s, ishft(s, 13))
      s = ieor(s, ishft(s, -7))
      s = ieor(s, ishft(s, 17))
      r = real(ishft(s, -11), wp)/9007199254740992.0_wp
   end function rnd

   function uni(a, b) result(u)
      real(wp), intent(in) :: a, b
      real(wp) :: u
      u = a + (b - a)*rnd()
   end function uni

   function logu(a, b) result(lu)
      !! log-uniform in [10**a, 10**b] to within a factor 2, from exact power-of-two scalings
      real(wp), intent(in) :: a, b
      real(wp) :: lu
      real(wp) :: m, e
      m = 1.0_wp + rnd()
      e = uni(a*3.3219280948873623_wp, b*3.3219280948873623_wp)
      lu = scale(m, nint(e))
   end function logu

   subroutine pick(f, i, x, y)
      integer, intent(in) :: f, i
      real(wp), intent(out) :: x, y
      integer, parameter :: NLM = 12, NAR = 9
      real(wp), parameter :: PI = 3.14159265358979323846_wp
      real(wp), parameter :: LY(NLM) = [0.5_wp, 10.0_wp, -2.0_wp, 3.0_wp, 0.3333333333333333_wp, 0.5_wp, &
                                        -1000.0_wp, 0.5_wp, 2.0_wp, 0.5_wp, 1.5_wp, 1.1_wp]
      real(wp), parameter :: AX2(NAR) = [1.0_wp, 1.0_wp, -1.0_wp, -1.0_wp, 0.0_wp, 2.0_wp, -5.0_wp, 1.0e-200_wp, 1.0e10_wp]
      real(wp), parameter :: AY2(NAR) = [1.0_wp, -1.0_wp, 1.0_wp, -1.0_wp, 3.0_wp, 0.0_wp, 0.0_wp, 1.0e-200_wp, 1.0e-300_wp]
      real(wp) :: lm(NLM)
      y = 0.0_wp
      select case (f)
      case (1)
         lm = [0.0_wp, 1.0_wp, -1.0_wp, 0.5_wp, 0.693147180559945309_wp, 709.78_wp, -707.0_wp, 700.0_wp, -700.0_wp, &
               1.0e-10_wp, -1.0e-300_wp, 1.0e-17_wp]
         if (i <= 12) then
            x = lm(i)
         else if (i <= 32) then
            x = uni(-40.0_wp, 40.0_wp)
         else if (i <= 40) then
            x = uni(-745.0_wp, -708.0_wp)
         else
            x = uni(-708.0_wp, 709.0_wp)
         end if
      case (2)
         lm = [1.0_wp, 2.0_wp, 0.5_wp, 2.718281828459045_wp, 1.0e-300_wp, 1.0e300_wp, tiny(1.0_wp), &
               1.0000000000000002_wp, 0.9999999999999999_wp, 3.0e-308_wp, 10.0_wp, 0.1_wp]
         if (i <= 12) then
            x = lm(i)
         else if (i <= 34) then
            x = logu(-30.0_wp, 30.0_wp)
         else
            x = uni(0.5_wp, 2.0_wp)
         end if
      case (3, 4, 5)
         lm = [0.0_wp, 1.0_wp, -1.0_wp, 1.5707963267948966_wp, PI, 1.0e3_wp, 4.5e5_wp, 1.5e6_wp, -4.5e5_wp, 1.0e-5_wp, &
               0.7853981633974483_wp, 6.283185307179586_wp]
         if (i <= 12) then
            x = lm(i)
         else if (i <= 28) then
            x = uni(-10.0_wp, 10.0_wp)
         else if (i <= 40) then
            x = uni(0.0_wp, 4.5e5_wp)
         else
            x = uni(-1.6e6_wp, 1.6e6_wp)
         end if
      case (6)
         lm = [0.0_wp, 1.0_wp, -1.0_wp, 1.0e-10_wp, 1.0e10_wp, 0.2679491924311227_wp, 0.5_wp, 2.0_wp, 1.0e300_wp, &
               -1.0e-300_wp, 3.0_wp, 0.3_wp]
         if (i <= 12) then
            x = lm(i)
         else if (i <= 30) then
            x = uni(-10.0_wp, 10.0_wp)
         else
            x = logu(-20.0_wp, 20.0_wp)
            y = uni(-1.0_wp, 1.0_wp)
            x = sign(x, y)
            y = 0.0_wp
         end if
      case (7)
         if (i <= NAR) then
            x = AX2(i)
            y = AY2(i)
         else
            x = uni(-10.0_wp, 10.0_wp)
            y = uni(-10.0_wp, 10.0_wp)
         end if
      case (8)
         lm = [0.0_wp, 1.0_wp, -1.0_wp, 1.0e-10_wp, 0.1_wp, 0.18_wp, 0.5_wp, 20.0_wp, -22.0_wp, 5.0e-9_wp, 3.0_wp, -0.01_wp]
         if (i <= 12) then
            x = lm(i)
         else if (i <= 30) then
            x = uni(-3.0_wp, 3.0_wp)
         else
            x = uni(-25.0_wp, 25.0_wp)
         end if
      case (9)
         lm = [0.0_wp, 1.0_wp, -1.0_wp, 1.0e-10_wp, 700.0_wp, 705.0_wp, -709.0_wp, 710.0_wp, 20.0_wp, 0.5_wp, 3.0_wp, -0.01_wp]
         if (i <= 12) then
            x = lm(i)
         else if (i <= 30) then
            x = uni(-5.0_wp, 5.0_wp)
         else
            x = uni(-710.0_wp, 710.0_wp)
         end if
      case (10)
         lm = [1.0_wp, 2.0_wp, 1.0000000000000002_wp, 1.5_wp, 10.0_wp, 1.0e8_wp, 1.0e10_wp, 3.0_wp, 1.001_wp, 1.0e300_wp, &
               268435456.0_wp, 1.9999999_wp]
         if (i <= 12) then
            x = lm(i)
         else if (i <= 28) then
            x = 1.0_wp + uni(0.0_wp, 1.0_wp)
         else
            x = 1.0_wp + logu(-12.0_wp, 12.0_wp)
         end if
      case (11)
         lm = [2.0_wp, 2.0_wp, 10.0_wp, 0.5_wp, 3.0_wp, 1.0e-3_wp, 2.0_wp, 9.0_wp, 7.0_wp, 100.0_wp, 0.25_wp, 1.1_wp]
         if (i <= 12) then
            x = lm(i)
            y = LY(i)
         else if (i == 13) then
            x = 1.0e-155_wp                  ! x**2 = 1e-310: flushed to exactly 0
            y = 2.0_wp
         else if (i == 14) then
            x = 2.0_wp                       ! 2**-1060 is subnormal: flushed to exactly 0
            y = -1060.0_wp
         else if (i <= 30) then
            x = logu(-3.0_wp, 3.0_wp)
            y = uni(-4.0_wp, 4.0_wp)
         else
            x = uni(0.01_wp, 1.0_wp)
            y = 1.0_wp/real(1 + int(uni(0.0_wp, 8.0_wp)), wp)
         end if
      case (12)
         lm = [1.0_wp, 8.0_wp, -8.0_wp, 27.0_wp, 2.0_wp, 1.0e-300_wp, 1.0e300_wp, tiny(1.0_wp), 0.125_wp, 1.0e-5_wp, &
               1000.0_wp, -0.001_wp]
         if (i <= 12) then
            x = lm(i)
         else if (i <= 30) then
            x = uni(-10.0_wp, 10.0_wp)
         else
            x = logu(-30.0_wp, 30.0_wp)
            y = uni(-1.0_wp, 1.0_wp)
            x = sign(x, y)
            y = 0.0_wp
         end if
      case (13)
         lm = [0.0_wp, 0.5_wp, 1.0_wp, -1.0_wp, 2.0_wp, 4.0_wp, 10.0_wp, 26.0_wp, -3.0_wp, 1.0e-10_wp, 0.25_wp, 26.5_wp]
         if (i <= 12) then
            x = lm(i)
         else if (i <= 24) then
            x = uni(-1.0_wp, 1.0_wp)
         else if (i <= 40) then
            x = uni(0.5_wp, 6.0_wp)
         else
            x = uni(4.0_wp, 27.0_wp)
         end if
      case default
         x = logu(-3.0_wp, 3.0_wp)
         y = uni(-60.0_wp, 60.0_wp)
      end select
   end subroutine pick

   function ev(f, x, y) result(r)
      integer, intent(in) :: f
      real(wp), intent(in) :: x, y
      real(wp) :: r
      select case (f)
      case (1)
         r = safe_exp(x)
      case (2)
         r = safe_log(x)
      case (3)
         r = safe_sin(x)
      case (4)
         r = safe_cos(x)
      case (5)
         r = safe_tan(x)
      case (6)
         r = safe_atan(x)
      case (7)
         r = safe_atan2(y, x)
      case (8)
         r = safe_tanh(x)
      case (9)
         r = safe_cosh(x)
      case (10)
         r = safe_acosh(x)
      case (11, 14)
         r = safe_pow(x, y)
      case (12)
         r = safe_cbrt(x)
      case default
         r = safe_erfc(x)
      end select
   end function ev

   subroutine row(f, i, ncol, last)
      integer, intent(in) :: f, i, ncol
      logical, intent(in) :: last
      real(wp) :: x, y, r
      call pick(f, i, x, y)
      r = ev(f, x, y)
      if (ncol == 3) then
         write (*, "(a)", advance="no") '      "'//hex(x)//'", "'//hex(y)//'", "'//hex(r)//'"'
      else
         write (*, "(a)", advance="no") '      "'//hex(x)//'", "'//hex(r)//'"'
      end if
      if (last) then
         write (*, "(a,i0,a,a,a)") "], [", ncol, ", NG_", trim(FN(f)), "])"
      else
         write (*, "(a)") ", &"
      end if
   end subroutine row

   function hex(x) result(h)
      real(wp), intent(in) :: x
      character(len=16) :: h
      write (h, "(Z16.16)") transfer(x, 1_int64)
   end function hex
end program gen_safe_math_golden
