! rdb_safe_math_core.inc -- error-free transformations, exp and log (bodies).
!
! PART OF the position-independent math include set (see rdb_safe_math_consts.inc and
! src/shared_module_utilities/README.md).  Provides: rdb_sm_two_sum, rdb_sm_two_prod,
! rdb_sm_exp_core, rdb_sm_exp, rdb_sm_reduce_log, rdb_sm_log, rdb_sm_log_dd.
! REQUIRES: rdb_safe_math_consts.inc in the specification part (include it ONCE).
! Required BY: _pow, _hyp, _erfc (include this file first).
!
! WHAT THE CONTRACT IS.  Every procedure in the set is a pure scalar function of its
! argument built ONLY from IEEE +,-,*,/ , sqrt, anint/nint, scale/exponent/fraction
! on finite normal operands, comparisons and compile-time tables.  No libm, no
! vector-library entry, no `transfer` in an executable statement (SM_INF and
! SM_NAN are constant expressions), no `shifta`, no FMA intrinsic, no `ieee_value`.
! So the value does not depend on where an element sits in a vector body / remainder
! loop, on the tile it belongs to, or on which math library the compiler links.  It
! depends on the compiler only through FMA contraction of a*b+c: the golden bits are a
! contract for contraction-OFF CPU builds (the project's global CPU default); the
! library adds no per-file flag and follows whatever the build does.
!
! SPECIAL VALUES propagate: NaN in gives NaN out, never a finite sentinel.  Each
! procedure tests ieee_is_nan / ieee_is_finite BEFORE any exponent extraction or
! integer conversion (exponent() of a non-finite value is garbage on a device).
! THE FLUSH RULE.  No procedure returns a subnormal result: every result with
! |r| < tiny(1.0_wp) is returned as exactly zero (the sign of a flushed zero is not guaranteed: ifx returns +0 where gfortran keeps the sign), through rdb_sm_flush below
! (exp, pow and erfc flush their own underflow the same way).  Compilers that flush
! denormals (nvfortran -fast) and compilers that do not therefore agree bitwise on
! results.  ARGUMENTS: a subnormal argument is not specially treated except by log
! (see there) and cbrt (read as zero); a compiler in denormals-are-zero mode reads it
! as 0 itself, so results for subnormal ARGUMENTS other than flushed ones are only
! compiler-independent where the function is continuous through 0 (all but log, and
! pow with a subnormal base).
!
! Method references: Cody & Waite (1980), Software Manual for the Elementary
! Functions (argument reduction with split constants; exp, log, sin, cos, atan,
! tanh); Muller (2016), Elementary Functions, 3rd ed. (ch. 4, 5, 11, 14);
! Dekker (1971) and Veltkamp (1968) (exact sum / product); Kahan (1991),
! Computing a Real Cube Root; Kahan (1983), log1p; Cody (1969), Math. Comp. 23
! (the erfc range split).  Fitted tables: tools/gen_safe_math_coeffs.py.

   pure function rdb_sm_flush(r_sm) result(y_sm)
      !! The flush rule: a result with |r| < `tiny` becomes zero (sign of a flushed zero not guaranteed); NaN and
      !! non-finite values pass through.
      !$acc routine seq
      real(wp), intent(in) :: r_sm
      real(wp) :: y_sm
      if (abs(r_sm) < tiny(r_sm)) then
         y_sm = sign(0.0_wp, r_sm)
      else
         y_sm = r_sm
      end if
   end function rdb_sm_flush

   pure subroutine rdb_sm_two_sum(a_sm, b_sm, s_sm, e_sm)
      !! Exact sum: `s + e == a + b`, `s = fl(a+b)` (Knuth).
      !$acc routine seq
      real(wp), intent(in) :: a_sm, b_sm
      real(wp), intent(out) :: s_sm, e_sm
      real(wp) :: bb_sm
      s_sm = a_sm + b_sm
      bb_sm = s_sm - a_sm
      e_sm = (a_sm - (s_sm - bb_sm)) + (b_sm - bb_sm)
   end subroutine rdb_sm_two_sum

   pure subroutine rdb_sm_two_prod(a_sm, b_sm, p_sm, e_sm)
      !! Exact product without FMA: `p + e == a * b`, `p = fl(a*b)` (Dekker/Veltkamp
      !! split; valid while |a|,|b| < 1e290 and the product does not under/overflow).
      !$acc routine seq
      real(wp), intent(in) :: a_sm, b_sm
      real(wp), intent(out) :: p_sm, e_sm
      real(wp) :: c_sm, ah_sm, al_sm, bh_sm, bl_sm
      c_sm = SM_SPLIT*a_sm
      ah_sm = c_sm - (c_sm - a_sm)
      al_sm = a_sm - ah_sm
      c_sm = SM_SPLIT*b_sm
      bh_sm = c_sm - (c_sm - b_sm)
      bl_sm = b_sm - bh_sm
      p_sm = a_sm*b_sm
      e_sm = (((ah_sm*bh_sm - p_sm) + ah_sm*bl_sm) + al_sm*bh_sm) + al_sm*bl_sm
   end subroutine rdb_sm_two_prod

   pure function rdb_sm_exp_core(x_sm, c_sm) result(y_sm)
      !! `exp(x + c)` for finite `x` in [-745.2, 709.8] and a small correction `c`
      !! (|c| <= ~2e-5: erfc hands in up to 1.1e-5, pow below one ulp of the
      !! head, so |r| stays within ln2/2 + |c| where the polynomial is accurate; pass
      !! 0 for plain exp).  A result below `tiny` is returned as exactly 0.
      !! `x` is reduced by k*ln2 (k = nearest integer) with the split ln2, `c` is added to the REDUCED argument, so the
      !! pair (x, c) carries a double-double exponent through (see `rdb_sm_pow`).
      !$acc routine seq
      real(wp), intent(in) :: x_sm, c_sm
      real(wp) :: y_sm
      real(wp) :: kr_sm, r_sm, p_sm
      integer :: k_sm, i_sm
      kr_sm = anint(x_sm*SM_INV_LN2)
      k_sm = nint(kr_sm)
      ! one term per statement: see the note in rdb_sm_log
      r_sm = x_sm - kr_sm*SM_LN2_HI                  ! exact
      r_sm = r_sm - kr_sm*SM_LN2_LO
      r_sm = r_sm + c_sm
      p_sm = SM_EXP_C(13)
      do i_sm = 12, 0, -1
         p_sm = p_sm*r_sm + SM_EXP_C(i_sm)
      end do
      ! 2**k in two steps near the top of the range: scale() overflowing is
      ! processor-dependent, a multiplication overflowing is IEEE.
      if (k_sm > 1000) then
         y_sm = 2.0_wp*scale(p_sm, k_sm - 1)
      else
         y_sm = scale(p_sm, k_sm)
      end if
      ! no subnormal results: below tiny the value is exactly 0, so compilers that flush
      ! denormals (nvfortran -fast) and compilers that do not agree bitwise, and no
      ! consumer ever meets a denormal operand
      if (y_sm < tiny(y_sm)) y_sm = 0.0_wp
   end function rdb_sm_exp_core

   pure function rdb_sm_exp(x_sm) result(y_sm)
      !! Exponential.  Defined for every real: NaN -> NaN, `+Inf` -> `+Inf`,
      !! `-Inf` -> 0, overflow -> `+Inf` (x > 709.78271289338), underflow -> 0
      !! (x < -745.13321910194).  Results below `tiny` (2.2e-308, i.e. x below about
      !! -708.4) are flushed to exactly 0 on every compiler: the library never produces
      !! a subnormal, so gfortran and nvfortran `-fast` (which flushes) agree bitwise
      !! and no consumer stalls on denormals.  Accuracy <= 1.2 ulp over
      !! the whole range (Cody & Waite 1980 ch. 6: k = nint(x/ln2), r = x - k ln2
      !! with the split ln2, degree-13 Taylor polynomial in |r| <= ln2/2).
      !$acc routine seq
      real(wp), intent(in) :: x_sm
      real(wp) :: y_sm
      if (ieee_is_nan(x_sm)) then
         y_sm = x_sm
      else if (x_sm > SM_EXP_MAX) then
         y_sm = SM_INF
      else if (x_sm < SM_EXP_MIN) then
         y_sm = 0.0_wp
      else
         y_sm = rdb_sm_exp_core(x_sm, 0.0_wp)
      end if
   end function rdb_sm_exp

   pure subroutine rdb_sm_reduce_log(x_sm, m_sm, k_sm)
      !! `x = m * 2**k` with m in [sqrt(1/2), sqrt(2)); `x` finite and > 0 and, if it is
      !! subnormal, not read as zero by the compiler (callers guard that: log and
      !! log_dd return `-Inf` for a subnormal that the compiler's denormals-are-zero
      !! mode turns into 0).  Subnormals are scaled by 2**54 first.
      !$acc routine seq
      real(wp), intent(in) :: x_sm
      real(wp), intent(out) :: m_sm
      integer, intent(out) :: k_sm
      real(wp) :: xs_sm
      integer :: kadj_sm
      if (x_sm < tiny(x_sm)) then
         xs_sm = x_sm*SM_TWO54
         kadj_sm = -54
      else
         xs_sm = x_sm
         kadj_sm = 0
      end if
      k_sm = exponent(xs_sm) + kadj_sm
      m_sm = fraction(xs_sm)
      if (m_sm < SM_SQRT_HALF) then
         m_sm = 2.0_wp*m_sm
         k_sm = k_sm - 1
      end if
   end subroutine rdb_sm_reduce_log

   pure function rdb_sm_log(x_sm) result(y_sm)
      !! Natural logarithm.  NaN -> NaN, x < 0 -> NaN, 0 -> `-Inf`, `+Inf` ->
      !! `+Inf`.  A subnormal x is evaluated exactly (scaled by 2**54) unless the
      !! compiler is in denormals-are-zero mode, where it is read as 0 and the result
      !! is `-Inf`.  Accuracy <= 3 ulp for every positive finite x (Cody & Waite 1980
      !! ch. 5 / Muller 2016 ch. 4: x = m 2**k, m in [sqrt(1/2), sqrt(2)),
      !! log m = 2 atanh f, f = (m-1)/(m+1), odd series in f, k ln2 with the split).
      !$acc routine seq
      real(wp), intent(in) :: x_sm
      real(wp) :: y_sm
      real(wp) :: m_sm, f_sm, u_sm, p_sm, kr_sm, s_sm, w_sm
      integer :: k_sm, i_sm
      if (ieee_is_nan(x_sm)) then
         y_sm = x_sm
      else if (x_sm < 0.0_wp) then
         y_sm = SM_NAN
      else if (x_sm == 0.0_wp) then
         y_sm = -SM_INF
      else if (x_sm < tiny(x_sm) .and. x_sm*SM_TWO54 == 0.0_wp) then
         y_sm = -SM_INF                    ! subnormal read as zero (denormals-are-zero)
      else if (.not. ieee_is_finite(x_sm)) then
         y_sm = x_sm
      else
         call rdb_sm_reduce_log(x_sm, m_sm, k_sm)
         f_sm = (m_sm - 1.0_wp)/(m_sm + 1.0_wp)
         u_sm = f_sm*f_sm
         p_sm = SM_LOG_C(9)
         do i_sm = 8, 1, -1
            p_sm = p_sm*u_sm + SM_LOG_C(i_sm)
         end do
         kr_sm = real(k_sm, wp)
         ! log m = 2f + 2 f u (1/3 + u/5 + ...): the head 2f is exact, the 1% tail
         ! carries the polynomial rounding; k ln2 joins as head + tail.  Three
         ! statements on purpose: nvfortran -O3 re-associated the same four terms
         ! written as ONE expression (-Mnofprelaxed restored it), which moved the
         ! result by 1 ulp; split like this, the golden table holds at -O3 -fast
         ! with no extra flag.
         s_sm = kr_sm*SM_LN2_HI + 2.0_wp*f_sm
         w_sm = (2.0_wp*f_sm)*u_sm*p_sm + kr_sm*SM_LN2_LO
         y_sm = s_sm + w_sm
      end if
   end function rdb_sm_log

   pure subroutine rdb_sm_log_dd(x_sm, lh_sm, ll_sm)
      !! `log(x) ~ lh + ll` for finite x > 0: the head `lh` is the leading part (not
      !! normalised against `ll`); `lh + ll` is correct to ~1e-18 relative, i.e.
      !! about 8 more bits than a plain double, which is what `rdb_sm_pow` needs so
      !! that y*log(x) does not amplify the logarithm's rounding error.
      !! Only the series tail (about 1% of the result) is in plain double.
      !$acc routine seq
      real(wp), intent(in) :: x_sm
      real(wp), intent(out) :: lh_sm, ll_sm
      real(wp) :: m_sm, a_sm, bh_sm, bl_sm, q1_sm, q2_sm, p1_sm, p2_sm, rr_sm
      real(wp) :: fh_sm, fl_sm, u_sm, p_sm, tail_sm, kr_sm, s_sm, e_sm
      integer :: k_sm, i_sm
      if (x_sm < tiny(x_sm) .and. x_sm*SM_TWO54 == 0.0_wp) then
         lh_sm = -SM_INF                   ! subnormal read as zero (denormals-are-zero)
         ll_sm = 0.0_wp
         return
      end if
      call rdb_sm_reduce_log(x_sm, m_sm, k_sm)
      a_sm = m_sm - 1.0_wp                      ! exact (Sterbenz)
      call rdb_sm_two_sum(m_sm, 1.0_wp, bh_sm, bl_sm)
      q1_sm = a_sm/bh_sm                        ! f = (m-1)/(m+1) = q1 + q2 + O(2**-106)
      call rdb_sm_two_prod(q1_sm, bh_sm, p1_sm, p2_sm)
      rr_sm = ((a_sm - p1_sm) - p2_sm) - q1_sm*bl_sm
      q2_sm = rr_sm/bh_sm
      fh_sm = q1_sm + q2_sm
      fl_sm = q2_sm - (fh_sm - q1_sm)
      u_sm = fh_sm*fh_sm
      p_sm = SM_LOG_C(9)
      do i_sm = 8, 1, -1
         p_sm = p_sm*u_sm + SM_LOG_C(i_sm)
      end do
      tail_sm = (2.0_wp*fh_sm)*u_sm*p_sm
      kr_sm = real(k_sm, wp)
      call rdb_sm_two_sum(kr_sm*SM_LN2_HI, 2.0_wp*fh_sm, s_sm, e_sm)
      lh_sm = s_sm
      ll_sm = e_sm + ((2.0_wp*fl_sm + tail_sm) + kr_sm*SM_LN2_LO)
   end subroutine rdb_sm_log_dd
