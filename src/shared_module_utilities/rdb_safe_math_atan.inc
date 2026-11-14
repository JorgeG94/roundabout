! rdb_safe_math_atan.inc -- atan and atan2 (bodies).
! Provides: rdb_sm_atan, rdb_sm_atan2.
! REQUIRES: rdb_safe_math_consts.inc (spec part) and rdb_safe_math_core.inc (for
! rdb_sm_flush), included BEFORE this file in the same `contains`.
! Contract: see rdb_safe_math_core.inc.

   pure function rdb_sm_atan(x_sm) result(y_sm)
      !! Arctangent in (-pi/2, pi/2).  NaN -> NaN, +-Inf -> +-pi/2, atan(-0) = -0.
      !! <= 2.5 ulp.  Method (Cody & Waite 1980 ch. 11): |x| > 1 -> pi/2 - atan(1/|x|);
      !! t > 2-sqrt(3) -> pi/6 + atan((sqrt(3) t - 1)/(sqrt(3) + t)); then
      !! a*P(a**2), a fitted degree-11 polynomial on |a| <= 2-sqrt(3).
      !$acc routine seq
      real(wp), intent(in) :: x_sm
      real(wp) :: y_sm
      real(wp) :: ax_sm, t_sm, a_sm, u_sm, p_sm
      integer :: i_sm
      logical :: inv_sm, mid_sm
      if (ieee_is_nan(x_sm)) then
         y_sm = x_sm
         return
      end if
      ax_sm = abs(x_sm)
      inv_sm = ax_sm > 1.0_wp
      if (inv_sm) then
         t_sm = 1.0_wp/ax_sm
      else
         t_sm = ax_sm
      end if
      mid_sm = t_sm > SM_TAN_PIO12
      if (mid_sm) then
         a_sm = (SM_SQRT3*t_sm - 1.0_wp)/(SM_SQRT3 + t_sm)
      else
         a_sm = t_sm
      end if
      u_sm = a_sm*a_sm
      p_sm = SM_ATAN_C(11)
      do i_sm = 10, 0, -1
         p_sm = p_sm*u_sm + SM_ATAN_C(i_sm)
      end do
      y_sm = a_sm*p_sm
      if (mid_sm) y_sm = SM_PIO6_HI + (y_sm + SM_PIO6_LO)
      if (inv_sm) y_sm = (SM_PIO2_HI - y_sm) + SM_PIO2_LO
      y_sm = rdb_sm_flush(sign(y_sm, x_sm))
   end function rdb_sm_atan

   pure function rdb_sm_atan2(y_sm, x_sm) result(r_sm)
      !! Four-quadrant arctangent in [-pi, pi] with the C99 special-value table
      !! (signed zeros, infinities); NaN if either argument is NaN.  <= 3 ulp.
      !$acc routine seq
      real(wp), intent(in) :: y_sm, x_sm
      real(wp) :: r_sm
      real(wp) :: ay_sm, ax_sm
      logical :: xneg_sm
      if (ieee_is_nan(x_sm) .or. ieee_is_nan(y_sm)) then
         r_sm = SM_NAN
         return
      end if
      ay_sm = abs(y_sm)
      ax_sm = abs(x_sm)
      xneg_sm = sign(1.0_wp, x_sm) < 0.0_wp
      if (ay_sm == 0.0_wp) then
         if (xneg_sm) then
            r_sm = SM_PI_HI + SM_PI_LO
         else
            r_sm = 0.0_wp
         end if
      else if (ax_sm == 0.0_wp) then
         r_sm = SM_PIO2_HI + SM_PIO2_LO
      else if (.not. ieee_is_finite(x_sm)) then
         if (.not. ieee_is_finite(y_sm)) then
            if (xneg_sm) then
               r_sm = 3.0_wp*SM_PIO4_HI + 3.0_wp*SM_PIO4_LO     ! 3 pi/4
            else
               r_sm = SM_PIO4_HI + SM_PIO4_LO
            end if
         else if (xneg_sm) then
            r_sm = SM_PI_HI + SM_PI_LO
         else
            r_sm = 0.0_wp
         end if
      else if (.not. ieee_is_finite(y_sm)) then
         r_sm = SM_PIO2_HI + SM_PIO2_LO
      else
         r_sm = rdb_sm_atan(ay_sm/ax_sm)
         if (xneg_sm) r_sm = (SM_PI_HI - r_sm) + SM_PI_LO
      end if
      r_sm = rdb_sm_flush(sign(r_sm, y_sm))
   end function rdb_sm_atan2
