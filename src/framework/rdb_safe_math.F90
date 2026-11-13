!! Position-independent transcendental functions: the host / setup / test face of
!! `src/shared_module_utilities/rdb_safe_math_*.inc`.
module rdb_safe_math
   !! **What this is.**  `exp`, `log`, `sin`, `cos`, `tan`, `atan`, `atan2`, `tanh`,
   !! `cosh`, `acosh`, `pow` (real `**` real), `cbrt` and `erfc`, written ONLY in
   !! IEEE `+ - * /`, `sqrt`, `anint`/`nint`, `scale`/`exponent`/`fraction` on
   !! finite operands and compile-time tables.  The value of `safe_exp(x)` is a
   !! function of `x` alone: it does not change with the element's position in a
   !! vector body or remainder loop, with the tile or rank the element lives on,
   !! with the compiler's vector math library (`libmvec`, `libpgmath`), or with the
   !! optimisation level.  The intrinsic `exp(x)` has none of those guarantees --
   !! a decomposed run and a one-rank run can round the same cell differently --
   !! which is why a model kernel that must reproduce bitwise across decompositions
   !! calls these instead.
   !!
   !! **Always on.**  There is no build option and no namelist knob: the `safe_*`
   !! names ARE the reproducible implementations.
   !!
   !! **One source text.**  The bodies live in `rdb_safe_math_{core,pow,trig,atan,hyp,erfc}.inc`
   !! (constants in `rdb_safe_math_consts.inc`, included once) so a kernel module can `#include` its own
   !! module-local copy (the form every toolchain inlines into a `do concurrent`
   !! body); this module includes the same text and exports elemental wrappers for
   !! host code, setup code and tests.  Both carry `!$acc routine seq`.
   !!
   !! **Special values** propagate (NaN in, NaN out; `exp(+Inf) = +Inf`,
   !! `log(0) = -Inf`, `log(-1) = NaN`, `pow` per IEEE-754/C99 except that a NaN
   !! argument always gives NaN) -- no `huge` sentinels.  No function returns a
   !! subnormal: every result with |r| < `tiny` is exactly zero (the sign of a flushed zero is not guaranteed: ifx returns +0 where gfortran keeps the sign) (the flush
   !! rule, `rdb_sm_flush`), so compilers that flush denormals and those that do not
   !! agree.  Subnormal ARGUMENTS: `log` of one is finite unless the compiler is in
   !! denormals-are-zero mode (then `-Inf`); `cbrt` reads it as 0; the others see
   !! whatever the compiler gives them and flush their result.  **Valid ranges** are in each function's docstring; `sin`, `cos`,
   !! `tan` return NaN beyond 1.6e6 rad (a phase that grows with time, e.g. tidal
   !! `omega*t`, must be reduced on the host with `modulo(phase, 2*pi)` first).
   !!
   !! **Contraction policy (what is and is not bit-reproducible).**  The library adds NO
   !! per-file compiler flag: it follows each build's global fused-multiply-add setting.
   !! The CPU builds of this repository compile with contraction OFF (gfortran
   !! `-ffp-contract=off`, nvfortran `-Mnofma`; ifx needs `-no-fma`); there the values
   !! are bitwise the same on every compiler and at every vector position, and the
   !! golden-bit test `tests/test_safe_math.F90` pins them.  The GPU build (and any
   !! contracting build) differs from that table by 1-2 ulp on some arguments, so the
   !! test does NOT assert the table there; it asserts that device array results equal
   !! device scalar calls of the same included copy, the ulp bounds, and the special
   !! values.  Working precision must be real64; the constants are real64 splits and
   !! the build stops otherwise (`SM_WP_IS_REAL64`).
   !!
   !! **Not provided**: `sqrt` is IEEE-exact on every target and needs no wrapper
   !! beyond `safe_sqrt`; integer powers `x**n` are exact products.
   use, intrinsic :: iso_fortran_env, only: int64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
   use rdb_constants, only: wp
   implicit none
   private

   public :: safe_exp, safe_log, safe_sin, safe_cos, safe_tan
   public :: safe_atan, safe_atan2, safe_tanh, safe_cosh, safe_acosh
   public :: safe_pow, safe_cbrt, safe_erfc, safe_sqrt

#include "rdb_safe_math_consts.inc"

contains

   elemental function safe_exp(x) result(y)
      !! `exp(x)`; see `rdb_sm_exp`.
      !$acc routine seq
      real(wp), intent(in) :: x
      real(wp) :: y
      y = rdb_sm_exp(x)
   end function safe_exp

   elemental function safe_log(x) result(y)
      !! `log(x)`; see `rdb_sm_log`.
      !$acc routine seq
      real(wp), intent(in) :: x
      real(wp) :: y
      y = rdb_sm_log(x)
   end function safe_log

   elemental function safe_sin(x) result(y)
      !! `sin(x)`, |x| <= 1.6e6; see `rdb_sm_sin`.
      !$acc routine seq
      real(wp), intent(in) :: x
      real(wp) :: y
      y = rdb_sm_sin(x)
   end function safe_sin

   elemental function safe_cos(x) result(y)
      !! `cos(x)`, |x| <= 1.6e6; see `rdb_sm_cos`.
      !$acc routine seq
      real(wp), intent(in) :: x
      real(wp) :: y
      y = rdb_sm_cos(x)
   end function safe_cos

   elemental function safe_tan(x) result(y)
      !! `tan(x)`, |x| <= 1.6e6; see `rdb_sm_tan`.
      !$acc routine seq
      real(wp), intent(in) :: x
      real(wp) :: y
      y = rdb_sm_tan(x)
   end function safe_tan

   elemental function safe_atan(x) result(y)
      !! `atan(x)`; see `rdb_sm_atan`.
      !$acc routine seq
      real(wp), intent(in) :: x
      real(wp) :: y
      y = rdb_sm_atan(x)
   end function safe_atan

   elemental function safe_atan2(y, x) result(r)
      !! `atan2(y, x)`; see `rdb_sm_atan2`.
      !$acc routine seq
      real(wp), intent(in) :: y, x
      real(wp) :: r
      r = rdb_sm_atan2(y, x)
   end function safe_atan2

   elemental function safe_tanh(x) result(y)
      !! `tanh(x)`; see `rdb_sm_tanh`.
      !$acc routine seq
      real(wp), intent(in) :: x
      real(wp) :: y
      y = rdb_sm_tanh(x)
   end function safe_tanh

   elemental function safe_cosh(x) result(y)
      !! `cosh(x)`; see `rdb_sm_cosh`.
      !$acc routine seq
      real(wp), intent(in) :: x
      real(wp) :: y
      y = rdb_sm_cosh(x)
   end function safe_cosh

   elemental function safe_acosh(x) result(y)
      !! `acosh(x)`, x >= 1; see `rdb_sm_acosh`.
      !$acc routine seq
      real(wp), intent(in) :: x
      real(wp) :: y
      y = rdb_sm_acosh(x)
   end function safe_acosh

   elemental function safe_pow(base, expo) result(y)
      !! `base**expo` for real arguments; see `rdb_sm_pow`.
      !$acc routine seq
      real(wp), intent(in) :: base, expo
      real(wp) :: y
      y = rdb_sm_pow(base, expo)
   end function safe_pow

   elemental function safe_cbrt(x) result(y)
      !! Real cube root; see `rdb_sm_cbrt`.
      !$acc routine seq
      real(wp), intent(in) :: x
      real(wp) :: y
      y = rdb_sm_cbrt(x)
   end function safe_cbrt

   elemental function safe_erfc(x) result(y)
      !! `erfc(x)`; see `rdb_sm_erfc`.
      !$acc routine seq
      real(wp), intent(in) :: x
      real(wp) :: y
      y = rdb_sm_erfc(x)
   end function safe_erfc

   elemental function safe_sqrt(x) result(y)
      !! `sqrt(x)`: IEEE-754 requires a correctly rounded result, so the intrinsic
      !! is already reproducible on every conforming target.
      !$acc routine seq
      real(wp), intent(in) :: x
      real(wp) :: y
      y = sqrt(x)
   end function safe_sqrt

#include "rdb_safe_math_core.inc"
#include "rdb_safe_math_pow.inc"
#include "rdb_safe_math_trig.inc"
#include "rdb_safe_math_atan.inc"
#include "rdb_safe_math_hyp.inc"
#include "rdb_safe_math_erfc.inc"

end module rdb_safe_math
