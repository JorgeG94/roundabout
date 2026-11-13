! rdb_safe_math_pow.inc -- pow (real ** real) and cbrt (bodies).
! Provides: rdb_sm_pow, rdb_sm_cbrt.
! REQUIRES: rdb_safe_math_consts.inc (spec part) and rdb_safe_math_core.inc, included
! BEFORE this file in the same `contains`.  Contract: see rdb_safe_math_core.inc.

   pure function rdb_sm_pow(x_sm, y_sm) result(r_sm)
      !! `x**y` for real x, y with the IEEE-754 / C99 `pow` special-value table EXCEPT
      !! that NaN in ANY argument gives NaN: C99's `pow(NaN, 0) = pow(1, NaN) = 1` is
      !! deliberately NOT followed, because the project rule is that non-finite input
      !! must propagate, never be laundered into a finite result.  `pow(x, +-0) = 1`
      !! and `pow(1, y) = 1` for every non-NaN argument; a negative base with a non-integer exponent is NaN; a negative base with
      !! an integer exponent takes the sign (-1)**y; `pow(+-0, y)` and the infinite
      !! bases and exponents as in C99 F.9.4.4.  For finite x > 0 the result is
      !! exp(y log x) with the logarithm carried as a head + tail pair and the
      !! product y*log x formed exactly (Dekker 1971; Muller 2016 ch. 14), so the
      !! logarithm's rounding is not multiplied by y: <= ~2 ulp for |y log x| <= 32,
      !! growing by about 1 ulp per 100 beyond (the series tail is plain double;
      !! 4.2 ulp measured at |y log x| = 345).
      !! Overflow -> `+Inf`; a result below `tiny` (2.2e-308) is exactly 0 (see
      !! `rdb_sm_exp`).
      !$acc routine seq
      real(wp), intent(in) :: x_sm, y_sm
      real(wp) :: r_sm
      real(wp) :: ax_sm, ay_sm, lh_sm, ll_sm, t_sm, ph_sm, pl_sm, c_sm, c2_sm, e_sm, h_sm
      logical :: yint_sm, yodd_sm, xneg_sm
      if (ieee_is_nan(x_sm) .or. ieee_is_nan(y_sm)) then
         r_sm = SM_NAN
         return
      end if
      if (y_sm == 0.0_wp .or. x_sm == 1.0_wp) then
         r_sm = 1.0_wp
         return
      end if
      ax_sm = abs(x_sm)
      ay_sm = abs(y_sm)
      ! integer-ness and parity of y (every |y| >= 2**53 is an even integer)
      if (ieee_is_finite(y_sm)) then
         if (ay_sm >= 9007199254740992.0_wp) then
            yint_sm = .true.
            yodd_sm = .false.
         else
            yint_sm = (ay_sm == anint(ay_sm))
            h_sm = 0.5_wp*ay_sm
            yodd_sm = yint_sm .and. (h_sm /= anint(h_sm))
         end if
      else
         yint_sm = .false.
         yodd_sm = .false.
      end if
      xneg_sm = (x_sm < 0.0_wp) .or. (x_sm == 0.0_wp .and. sign(1.0_wp, x_sm) < 0.0_wp)
      if (.not. ieee_is_finite(y_sm)) then
         if (ax_sm == 1.0_wp) then
            r_sm = 1.0_wp
         else if ((ax_sm < 1.0_wp) .eqv. (y_sm > 0.0_wp)) then
            r_sm = 0.0_wp
         else
            r_sm = SM_INF
         end if
         return
      end if
      if (ax_sm == 0.0_wp) then
         if (y_sm > 0.0_wp) then
            r_sm = 0.0_wp
         else
            r_sm = SM_INF
         end if
         if (xneg_sm .and. yodd_sm) r_sm = -r_sm
         return
      end if
      if (.not. ieee_is_finite(x_sm)) then
         if (y_sm > 0.0_wp) then
            r_sm = SM_INF
         else
            r_sm = 0.0_wp
         end if
         if (xneg_sm .and. yodd_sm) r_sm = -r_sm
         return
      end if
      if (xneg_sm .and. .not. yint_sm) then
         r_sm = SM_NAN
         return
      end if
      ! finite ax > 0, finite y: exp(y log ax)
      call rdb_sm_log_dd(ax_sm, lh_sm, ll_sm)
      t_sm = y_sm*lh_sm
      if (t_sm > 800.0_wp) then
         r_sm = SM_INF
      else if (t_sm < -800.0_wp) then
         r_sm = 0.0_wp
      else
         call rdb_sm_two_prod(y_sm, lh_sm, ph_sm, pl_sm)
         c_sm = pl_sm + y_sm*ll_sm
         ! fold the (possibly sizeable) tail into the head so the correction handed to
         ! exp_core is below one ulp of the head
         call rdb_sm_two_sum(ph_sm, c_sm, e_sm, c2_sm)
         if (e_sm > SM_EXP_MAX) then
            r_sm = SM_INF
         else if (e_sm < SM_EXP_MIN) then
            r_sm = 0.0_wp
         else
            r_sm = rdb_sm_exp_core(e_sm, c2_sm)
         end if
      end if
      if (xneg_sm .and. yodd_sm) r_sm = -r_sm
   end function rdb_sm_pow

   pure function rdb_sm_cbrt(x_sm) result(y_sm)
      !! Real cube root, odd in x.  NaN -> NaN, +-Inf -> +-Inf, +-0 -> +-0.  A subnormal
      !! argument is treated as zero (+-0 returned) on every compiler, so the result
      !! does not depend on the compiler's denormals-are-zero mode.
      !! <= 1.5 ulp (Kahan 1991; Muller 2016 ch. 5): the binary exponent is split
      !! modulo 3, the mantissa part in [1/2, 4) gets a degree-4 seed (0.5% error)
      !! and two Halley steps `y <- y + y (m - y**3) / (2 y**3 + m)`.
      !$acc routine seq
      real(wp), intent(in) :: x_sm
      real(wp) :: y_sm
      real(wp) :: ax_sm, m_sm, y0_sm, y3_sm
      integer :: e_sm, q_sm, rem_sm, i_sm
      if (ieee_is_nan(x_sm) .or. .not. ieee_is_finite(x_sm)) then
         y_sm = x_sm
         return
      end if
      ax_sm = abs(x_sm)
      if (ax_sm < tiny(ax_sm)) then
         y_sm = sign(0.0_wp, x_sm)
         return
      end if
      e_sm = exponent(ax_sm)
      rem_sm = modulo(e_sm, 3)
      q_sm = (e_sm - rem_sm)/3
      m_sm = scale(fraction(ax_sm), rem_sm)     ! in [1/2, 4)
      y0_sm = SM_CBRT_C(4)
      do i_sm = 3, 0, -1
         y0_sm = y0_sm*m_sm + SM_CBRT_C(i_sm)
      end do
      do i_sm = 1, 2
         y3_sm = (y0_sm*y0_sm)*y0_sm
         y0_sm = y0_sm + y0_sm*((m_sm - y3_sm)/((y3_sm + y3_sm) + m_sm))
      end do
      y_sm = scale(y0_sm, q_sm)
      if (x_sm < 0.0_wp) y_sm = -y_sm
   end function rdb_sm_cbrt
