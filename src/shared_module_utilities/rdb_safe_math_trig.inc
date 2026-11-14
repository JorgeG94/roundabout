! rdb_safe_math_trig.inc -- sin, cos, tan (bodies).
! Provides: rdb_sm_reduce_pio2, rdb_sm_sin_poly, rdb_sm_cos_poly, rdb_sm_trig,
! rdb_sm_sin, rdb_sm_cos, rdb_sm_tan.
! REQUIRES: rdb_safe_math_consts.inc (spec part) and rdb_safe_math_core.inc (for
! rdb_sm_flush), included BEFORE this file in the same `contains`.
! Contract: see rdb_safe_math_core.inc.

   pure subroutine rdb_sm_reduce_pio2(x_sm, r_sm, q_sm)
      !! `x = k pi/2 + r`, |r| <= pi/4, `q = k mod 4`, for |x| <= SM_TRIG_MAX.
      !! Three-part Cody & Waite (1980, ch. 8) reduction: k*P1 and k*P2 are exact
      !! (31- and 32-bit constants, |k| < 2**20) so only the last, rounded, subtraction
      !! contributes error.
      !$acc routine seq
      real(wp), intent(in) :: x_sm
      real(wp), intent(out) :: r_sm
      integer, intent(out) :: q_sm
      real(wp) :: kr_sm
      kr_sm = anint(x_sm*SM_INV_PIO2)
      q_sm = iand(nint(kr_sm), 3)
      if (kr_sm == 0.0_wp) then
         r_sm = x_sm                                 ! keeps the sign of a zero argument
      else
         r_sm = x_sm - kr_sm*SM_PIO2_1               ! k*P1 exact, so is the difference
         r_sm = r_sm - kr_sm*SM_PIO2_2               ! likewise
         r_sm = r_sm - kr_sm*SM_PIO2_3
      end if
   end subroutine rdb_sm_reduce_pio2

   pure function rdb_sm_sin_poly(r_sm) result(y_sm)
      !! sin(r), |r| <= pi/4: degree-17 odd Taylor polynomial.
      !$acc routine seq
      real(wp), intent(in) :: r_sm
      real(wp) :: y_sm
      real(wp) :: u_sm, p_sm
      integer :: i_sm
      u_sm = r_sm*r_sm
      p_sm = SM_SIN_C(8)
      do i_sm = 7, 0, -1
         p_sm = p_sm*u_sm + SM_SIN_C(i_sm)
      end do
      y_sm = r_sm*p_sm
   end function rdb_sm_sin_poly

   pure function rdb_sm_cos_poly(r_sm) result(y_sm)
      !! cos(r), |r| <= pi/4: degree-16 even Taylor polynomial.
      !$acc routine seq
      real(wp), intent(in) :: r_sm
      real(wp) :: y_sm
      real(wp) :: u_sm, p_sm
      integer :: i_sm
      u_sm = r_sm*r_sm
      p_sm = SM_COS_C(8)
      do i_sm = 7, 0, -1
         p_sm = p_sm*u_sm + SM_COS_C(i_sm)
      end do
      y_sm = p_sm
   end function rdb_sm_cos_poly

   pure function rdb_sm_trig(x_sm, shift_sm) result(y_sm)
      !! sin (shift 0) or cos (shift 1) of x; shared body of the two entry points.
      !$acc routine seq
      real(wp), intent(in) :: x_sm
      integer, intent(in) :: shift_sm
      real(wp) :: y_sm
      real(wp) :: r_sm
      integer :: q_sm
      if (ieee_is_nan(x_sm)) then
         y_sm = x_sm
      else if (.not. (abs(x_sm) <= SM_TRIG_MAX)) then
         y_sm = SM_NAN
      else
         call rdb_sm_reduce_pio2(x_sm, r_sm, q_sm)
         q_sm = iand(q_sm + shift_sm, 3)
         if (q_sm == 0) then
            y_sm = rdb_sm_sin_poly(r_sm)
         else if (q_sm == 1) then
            y_sm = rdb_sm_cos_poly(r_sm)
         else if (q_sm == 2) then
            y_sm = -rdb_sm_sin_poly(r_sm)
         else
            y_sm = -rdb_sm_cos_poly(r_sm)
         end if
         y_sm = rdb_sm_flush(y_sm)
      end if
   end function rdb_sm_trig

   pure function rdb_sm_sin(x_sm) result(y_sm)
      !! Sine.  Valid for |x| <= 1.6e6 rad (k = nint(2x/pi) < 2**20, where the
      !! reduction products are exact); <= 3 ulp there (relative, including near
      !! the zeros).  NaN -> NaN; +-Inf and |x| > 1.6e6 -> NaN: the caller removes
      !! whole periods first (a tidal phase `omega*t` exceeds the range after ~370
      !! years of M2).  A caller whose phase grows with model time (tides: omega*t) must
      !! reduce it on the host, `modulo(phase, 2*pi)`, before the kernel.
      !! Method: Cody & Waite (1980) ch. 8.
      !$acc routine seq
      real(wp), intent(in) :: x_sm
      real(wp) :: y_sm
      y_sm = rdb_sm_trig(x_sm, 0)
   end function rdb_sm_sin

   pure function rdb_sm_cos(x_sm) result(y_sm)
      !! Cosine.  Same range, accuracy and special values as `rdb_sm_sin`.
      !$acc routine seq
      real(wp), intent(in) :: x_sm
      real(wp) :: y_sm
      y_sm = rdb_sm_trig(x_sm, 1)
   end function rdb_sm_cos

   pure function rdb_sm_tan(x_sm) result(y_sm)
      !! Tangent, sin/cos from the same reduced argument.  Same range and special
      !! values as `rdb_sm_sin`; <= 5 ulp away from the poles (Cody & Waite 1980
      !! ch. 9).
      !$acc routine seq
      real(wp), intent(in) :: x_sm
      real(wp) :: y_sm
      real(wp) :: r_sm
      integer :: q_sm
      if (ieee_is_nan(x_sm)) then
         y_sm = x_sm
      else if (.not. (abs(x_sm) <= SM_TRIG_MAX)) then
         y_sm = SM_NAN
      else
         call rdb_sm_reduce_pio2(x_sm, r_sm, q_sm)
         if (iand(q_sm, 1) == 0) then
            y_sm = rdb_sm_sin_poly(r_sm)/rdb_sm_cos_poly(r_sm)
         else
            y_sm = -rdb_sm_cos_poly(r_sm)/rdb_sm_sin_poly(r_sm)
         end if
         y_sm = rdb_sm_flush(y_sm)
      end if
   end function rdb_sm_tan
