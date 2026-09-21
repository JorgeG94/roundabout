!! Deterministic pseudo-random values for tests, without signed overflow.
module lcg_deterministic
   !! A linear congruential generator whose arithmetic is DEFINED.
   !!
   !! Tests here need value streams that are reproducible bit-for-bit
   !! across compilers -- `random_number()` is not, so several suites grew
   !! their own copy of Knuth's MMIX LCG:
   !!
   !!     s = 6364136223846793005*seed + 1442695040888963407
   !!
   !! That product overflows `integer(int64)` on essentially every call.
   !! Wraparound is exactly what an LCG wants, but Fortran does not define
   !! signed integer overflow, and GCC treats it as undefined behaviour --
   !! `-Waggressive-loop-optimizations` reports it, having PROVEN the
   !! overflow and reserved the right to optimise as though it cannot
   !! happen. A generator chosen for cross-compiler reproducibility was
   !! therefore the one construct in the suite that could not guarantee it.
   !!
   !! The state here is the low 48 bits, which is all the callers ever
   !! consumed, and the multiply is done in 24-bit halves so that no
   !! intermediate can exceed 2**50 (measured worst case 5.8e14, a factor
   !! ~1.6e4 below `huge(0_int64)`). Every operation is defined, and the
   !! result is bit-identical to what the wrapping form produced:
   !! `(A*seed + C) mod 2**48` depends only on `seed mod 2**48`.
   !!
   !! NOT yet used by `test_ocean_ppm_h4_remap` / `test_remap_pqm`, which
   !! carry the same UB in a different shape: they take
   !! `mod(state*A + C, 9223372036854775783_8)` -- a modulus of an ALREADY
   !! WRAPPED product, which cannot be reproduced without emulating the
   !! wraparound in full 64-bit. Substituting this generator there is not
   !! value-preserving, and `test_order_of_accuracy` asserts a fitted
   !! convergence slope with a hand-tuned 0.1 margin between two schemes,
   !! so it is sensitive to the particular grid REALISATION rather than to
   !! grid quality (measured: old and new grids agree to within 1% on both
   !! max/min and worst adjacent jump, and the assertion still failed).
   !! Migrating them needs either a defined 64-bit wraparound helper or a
   !! more robust assertion, plus a local test run -- tracked separately.
   use, intrinsic :: iso_fortran_env, only: int64, real64
   implicit none
   private

   public :: lcg_next
   public :: lcg_unit

   integer(int64), parameter :: M48 = 281474976710655_int64  !! 2**48 - 1
   integer(int64), parameter :: M24 = 16777215_int64         !! 2**24 - 1
   integer(int64), parameter :: A_HI = 16002380_int64        !! (A mod 2**48) >> 24
   integer(int64), parameter :: A_LO = 9797421_int64         !! (A mod 2**48) & M24
   integer(int64), parameter :: C48 = 135785246851407_int64  !! C mod 2**48
   real(real64), parameter :: TWO48 = 281474976710656.0_real64  !! 2**48

contains

   pure function lcg_next(state) result(next)
      !! Advance the 48-bit LCG state. Any `state` is accepted; only its
      !! low 48 bits affect the result, so a caller may seed with a loop
      !! index (or a negative value) without masking first.
      integer(int64), intent(in) :: state
      integer(int64) :: next
      integer(int64) :: s, s_hi, s_lo, mid

      s = iand(state, M48)
      s_hi = ishft(s, -24)
      s_lo = iand(s, M24)
      ! `mid` carries only the 24 bits that survive into the high half;
      ! masking it here is what keeps the final sum below 2**50.
      mid = iand(A_LO*s_hi + A_HI*s_lo, M24)
      next = iand(A_LO*s_lo + ishft(mid, 24) + C48, M48)
   end function lcg_next

   pure function lcg_unit(state) result(u)
      !! Map a 48-bit state to [0, 1). Uses the full 48 bits, not a
      !! low-order remainder -- the low bits of an LCG modulo a power of
      !! two have short periods and make poor value streams.
      integer(int64), intent(in) :: state
      real(real64) :: u

      u = real(iand(state, M48), real64)/TWO48
   end function lcg_unit

end module lcg_deterministic
