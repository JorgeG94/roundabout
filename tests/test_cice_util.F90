!! Shared helpers for the `rdb_test_cice_*` sea-ice tests (first slice).
!!
!! Linked from the small `rdb_cice_test_util` library so ice tests do not need
!! the ocean core.  Provides deterministic random columns, a per-line
!! conservation ledger and a golden-vector comparison.
module test_cice_util
   use, intrinsic :: iso_fortran_env, only: int64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_constants, only: wp
   use lcg_deterministic, only: lcg_next, lcg_unit
   use testdrive, only: error_type, test_failed
   implicit none
   private

   public :: ice_random_columns, ice_ledger_t, ice_golden_check, ice_golden_first_bad, ICE_LEDGER_MAX

   integer, parameter :: ICE_LEDGER_MAX = 16
      !! Maximum number of ledger lines.

   type :: ice_ledger_t
      !! Per-line closure bookkeeping: each line records the largest
      !! `|residual| / max(scale, tiny)` seen over all columns.
      integer :: n = 0
      character(len=24) :: name(ICE_LEDGER_MAX) = ''
      real(wp) :: worst(ICE_LEDGER_MAX) = 0.0_wp
      logical :: overflow = .false.   !! more than ICE_LEDGER_MAX lines were requested
   contains
      procedure :: add => ledger_add
      procedure :: worst_of => ledger_worst_of
      procedure :: closed => ledger_closed
   end type ice_ledger_t

contains

   subroutine ice_random_columns(n, nilyr, seed, hi, hs, Tsfc, S, T)
      !! `n` deterministic pseudo-random columns, bit-identical across
      !! compilers (48-bit LCG).  Thickness `hi` in [0.05, 5] m (denser at the
      !! thin end), snow `hs` in [0, 0.5] m (30 % bare), surface temperature
      !! `Tsfc` in [-40, -0.1] degC, layer salinity `S(k,c)` in [0.1, 6] ppt
      !! and layer temperature `T(k,c)` strictly between `Tsfc` and the layer's
      !! melting point (`-0.054*S`), top-down.
      integer, intent(in) :: n, nilyr
      integer(int64), intent(in) :: seed
      real(wp), intent(out) :: hi(n), hs(n), Tsfc(n), S(nilyr, n), T(nilyr, n)
      integer :: c, k
      integer(int64) :: st
      real(wp) :: u, tm

      st = lcg_next(seed)
      do c = 1, n
         st = lcg_next(st); u = lcg_unit(st)
         hi(c) = 0.05_wp + 4.95_wp*u*u
         st = lcg_next(st); u = lcg_unit(st)
         hs(c) = 0.0_wp
         if (u > 0.3_wp) then
            st = lcg_next(st); hs(c) = 0.5_wp*lcg_unit(st)
         end if
         st = lcg_next(st); Tsfc(c) = -0.1_wp - 39.9_wp*lcg_unit(st)
         do k = 1, nilyr
            st = lcg_next(st); S(k, c) = 0.1_wp + 5.9_wp*lcg_unit(st)
            tm = -0.054_wp*S(k, c)
            st = lcg_next(st); u = 0.02_wp + 0.96_wp*lcg_unit(st)
            T(k, c) = tm + u*(min(Tsfc(c), tm) - tm)
            if (T(k, c) >= tm) T(k, c) = tm - 0.01_wp
         end do
      end do
   end subroutine ice_random_columns

   subroutine ledger_add(self, line, residual, scale)
      !! Record `|residual|/max(|scale|, tiny)` on ledger line `line` (created
      !! on first use).
      class(ice_ledger_t), intent(inout) :: self
      character(len=*), intent(in) :: line
      real(wp), intent(in) :: residual, scale
      integer :: i
      real(wp) :: r

      i = find_line(self, line)
      if (i == 0) then
         if (self%n >= ICE_LEDGER_MAX) then
            self%overflow = .true.      ! seen by worst_of / closed: never a silent pass
            return
         end if
         self%n = self%n + 1
         i = self%n
         self%name(i) = line
      end if
      r = abs(residual)/max(abs(scale), tiny(1.0_wp))
      ! Non-finite residuals are recorded as `huge` explicitly: `max(worst, NaN)`
      ! hides a NaN, and a bare comparison is NaN-blind under relaxed FP.
      if (.not. ieee_is_finite(r)) r = huge(1.0_wp)
      if (r > self%worst(i)) self%worst(i) = r
   end subroutine ledger_add

   function ledger_worst_of(self, line) result(w)
      !! Worst relative residual of `line` (0 if the line was never added;
      !! `huge` if any residual was non-finite or the ledger overflowed).
      class(ice_ledger_t), intent(in) :: self
      character(len=*), intent(in) :: line
      real(wp) :: w
      integer :: i

      i = find_line(self, line)
      w = 0.0_wp
      if (i > 0) w = self%worst(i)
      if (self%overflow) w = huge(1.0_wp)
   end function ledger_worst_of

   function ledger_closed(self, tol) result(ok)
      !! True when every line closed to within `tol` (relative).
      class(ice_ledger_t), intent(in) :: self
      real(wp), intent(in) :: tol
      logical :: ok

      ok = .not. self%overflow .and. all(self%worst(1:self%n) <= tol)
   end function ledger_closed

   function find_line(self, line) result(i)
      type(ice_ledger_t), intent(in) :: self
      character(len=*), intent(in) :: line
      integer :: i, j

      i = 0
      do j = 1, self%n
         if (trim(self%name(j)) == trim(line)) then
            i = j
            exit
         end if
      end do
   end function find_line

   subroutine ice_golden_check(error, label, got, expected, rtol, atol)
      !! Compare `got` against a vector of closed-form / reference values:
      !! passes where `|got-expected| <= atol + rtol*|expected|`; on the first
      !! failure sets `error` with the label, index and both values.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: label
      real(wp), intent(in) :: got(:), expected(:), rtol, atol
      integer :: i
      character(len=160) :: msg

      if (size(got) /= size(expected)) then
         call test_failed(error, trim(label)//": length mismatch")
         return
      end if
      i = ice_golden_first_bad(got, expected, rtol, atol)
      if (i > 0) then
         write (msg, '(a,a,i0,2(a,es23.15))') trim(label), ": index ", i, &
            " got ", got(i), " expected ", expected(i)
         call test_failed(error, trim(msg))
      end if
   end subroutine ice_golden_check

   pure function ice_golden_first_bad(got, expected, rtol, atol) result(i)
      !! Index of the first element with `|got-expected| > atol + rtol*|expected|`
      !! (a NaN counts as bad), or 0 when all agree.
      real(wp), intent(in) :: got(:), expected(:), rtol, atol
      integer :: i, j

      i = 0
      do j = 1, size(got)
         if (.not. ieee_is_finite(got(j)) .or. &
             .not. (abs(got(j) - expected(j)) <= atol + rtol*abs(expected(j)))) then
            i = j
            return
         end if
      end do
   end function ice_golden_first_bad

end module test_cice_util
