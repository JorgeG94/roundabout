! rdb_safe_math_hyp.inc -- tanh, cosh, acosh (bodies).
! Provides: rdb_sm_expm1_small, rdb_sm_tanh, rdb_sm_cosh, rdb_sm_log1p, rdb_sm_acosh.
! REQUIRES: rdb_safe_math_consts.inc (spec part) and rdb_safe_math_core.inc, included
! BEFORE this file in the same `contains`.  Contract: see rdb_safe_math_core.inc.

   pure function rdb_sm_expm1_small(x_sm) result(y_sm)
      !! `exp(x) - 1` for |x| <= 0.36 without cancellation: `x + x**2 * E(x)`.
      !$acc routine seq
      real(wp), intent(in) :: x_sm
      real(wp) :: y_sm
      real(wp) :: p_sm
      integer :: i_sm
      p_sm = SM_EXP_C(13)
      do i_sm = 12, 2, -1
         p_sm = p_sm*x_sm + SM_EXP_C(i_sm)
      end do
      y_sm = x_sm + (x_sm*x_sm)*p_sm
   end function rdb_sm_expm1_small

   pure function rdb_sm_tanh(x_sm) result(y_sm)
      !! Hyperbolic tangent.  NaN -> NaN, +-Inf -> +-1.  <= 4 ulp.  For |x| >= 0.18
      !! (via t = expm1(2|x|)): tanh = t/(t+2), with expm1 from the exp
      !! reduction below 0.36 and `exp(.) - 1` above (Cody & Waite 1980 ch. 12).
      !$acc routine seq
      real(wp), intent(in) :: x_sm
      real(wp) :: y_sm
      real(wp) :: ax_sm, t_sm
      if (ieee_is_nan(x_sm)) then
         y_sm = x_sm
         return
      end if
      ax_sm = abs(x_sm)
      if (ax_sm < 7.0e-9_wp) then
         y_sm = rdb_sm_flush(x_sm)
      else if (ax_sm > 22.0_wp) then
         y_sm = sign(1.0_wp, x_sm)
      else
         if (ax_sm < 0.18_wp) then
            t_sm = rdb_sm_expm1_small(2.0_wp*ax_sm)
         else
            t_sm = rdb_sm_exp_core(2.0_wp*ax_sm, 0.0_wp) - 1.0_wp
         end if
         y_sm = sign(t_sm/(t_sm + 2.0_wp), x_sm)
      end if
   end function rdb_sm_tanh

   pure function rdb_sm_cosh(x_sm) result(y_sm)
      !! Hyperbolic cosine.  NaN -> NaN, +-Inf -> +Inf, overflow (|x| > 710.4758)
      !! -> +Inf.  <= 2.5 ulp.
      !$acc routine seq
      real(wp), intent(in) :: x_sm
      real(wp) :: y_sm
      real(wp) :: ax_sm, e_sm
      if (ieee_is_nan(x_sm)) then
         y_sm = x_sm
         return
      end if
      ax_sm = abs(x_sm)
      if (ax_sm > 710.5_wp) then
         y_sm = SM_INF
      else if (ax_sm > 700.0_wp) then
         ! 0.5 exp(a) = exp(a - 1) * (e/2): a-1 is exact here
         y_sm = rdb_sm_exp_core(ax_sm - 1.0_wp, 0.0_wp)*SM_HALF_E
      else
         e_sm = rdb_sm_exp_core(ax_sm, 0.0_wp)
         y_sm = 0.5_wp*(e_sm + 1.0_wp/e_sm)
      end if
   end function rdb_sm_cosh

   pure function rdb_sm_log1p(t_sm) result(y_sm)
      !! log(1+t) for t >= 0 without cancellation (Kahan 1983): the logarithm of the
      !! ROUNDED 1+t times t over the exact increment (u-1).
      !$acc routine seq
      real(wp), intent(in) :: t_sm
      real(wp) :: y_sm
      real(wp) :: u_sm
      u_sm = 1.0_wp + t_sm
      if (u_sm == 1.0_wp) then
         y_sm = t_sm
      else
         y_sm = rdb_sm_log(u_sm)*(t_sm/(u_sm - 1.0_wp))
      end if
   end function rdb_sm_log1p

   pure function rdb_sm_acosh(x_sm) result(y_sm)
      !! Inverse hyperbolic cosine for x >= 1.  NaN -> NaN, x < 1 -> NaN,
      !! `+Inf` -> `+Inf`.  <= 5 ulp.  With t = x-1: log1p(t + sqrt(2t + t*t)) up
      !! to x = 2, log(2x - 1/(x + sqrt(x*x-1))) beyond, log(x) + ln2 above 2**28
      !! (Cody & Waite 1980 ch. 13).
      !$acc routine seq
      real(wp), intent(in) :: x_sm
      real(wp) :: y_sm
      real(wp) :: t_sm
      if (ieee_is_nan(x_sm)) then
         y_sm = x_sm
      else if (x_sm < 1.0_wp) then
         y_sm = SM_NAN
      else if (.not. ieee_is_finite(x_sm)) then
         y_sm = x_sm
      else if (x_sm <= 2.0_wp) then
         t_sm = x_sm - 1.0_wp
         y_sm = rdb_sm_log1p(t_sm + sqrt((t_sm + t_sm) + t_sm*t_sm))
      else if (x_sm <= 268435456.0_wp) then
         y_sm = rdb_sm_log((x_sm + x_sm) - 1.0_wp/(x_sm + sqrt(x_sm*x_sm - 1.0_wp)))
      else
         y_sm = rdb_sm_log(x_sm) + (SM_LN2_HI + SM_LN2_LO)
      end if
   end function rdb_sm_acosh
