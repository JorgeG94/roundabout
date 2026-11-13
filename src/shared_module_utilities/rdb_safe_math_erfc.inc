! rdb_safe_math_erfc.inc -- erfc (body).
! Provides: rdb_sm_erfc.
! REQUIRES: rdb_safe_math_consts.inc (spec part) and rdb_safe_math_core.inc, included
! BEFORE this file in the same `contains`.  Contract: see rdb_safe_math_core.inc.

   pure function rdb_sm_erfc(x_sm) result(y_sm)
      !! Complementary error function.  NaN -> NaN, `+Inf` -> 0, `-Inf` -> 2,
      !! x > 27.3 -> 0 (the true value is below 1e-324); a result below `tiny` is
      !! exactly 0.  <= 7 ulp (relative
      !! error <= 1e-15 measured over the range; two Horner chains and a quotient).
      !! |x| < 0.5: 1 - x*E(x**2), E a degree-11 fit of erf(sqrt(v))/sqrt(v).
      !! 0.5 <= |x| < 4: exp(-x**2) R(x), R = exp(x**2) erfc(x) a (6,6) rational
      !! fit on [0.5,2] and on [2,4]; |x| >= 4: exp(-x**2)/(|x| sqrt(pi)) (1 + w F(w)),
      !! w = 1/x**2, F a (5,5) rational fit.  exp(-x**2) is formed as one
      !! exponential of the exactly split argument (Cody 1969: the same range split;
      !! tables of this library are fitted anew, tools/gen_safe_math_coeffs.py).
      !$acc routine seq
      real(wp), intent(in) :: x_sm
      real(wp) :: y_sm
      real(wp) :: ax_sm, v_sm, p_sm, q_sm, c_sm, xh_sm, xl_sm, ex_sm, w_sm, rax_sm
      integer :: i_sm
      if (ieee_is_nan(x_sm)) then
         y_sm = x_sm
         return
      end if
      ax_sm = abs(x_sm)
      if (ax_sm < 0.5_wp) then
         v_sm = x_sm*x_sm
         p_sm = SM_ERF_C(11)
         do i_sm = 10, 0, -1
            p_sm = p_sm*v_sm + SM_ERF_C(i_sm)
         end do
         y_sm = 1.0_wp - x_sm*p_sm
         return
      end if
      if (ax_sm > 27.3_wp) then
         ex_sm = 0.0_wp
      else
         ! exp(-x**2): xh has 26 significant bits so xh*xh is exact
         c_sm = SM_SPLIT*ax_sm
         xh_sm = c_sm - (c_sm - ax_sm)
         xl_sm = ax_sm - xh_sm
         ex_sm = rdb_sm_exp_core(-(xh_sm*xh_sm), -(xl_sm*(ax_sm + xh_sm)))
      end if
      if (ax_sm < 2.0_wp) then
         p_sm = SM_ERFC_A_P(6)
         q_sm = SM_ERFC_A_Q(6)
         do i_sm = 5, 0, -1
            p_sm = p_sm*ax_sm + SM_ERFC_A_P(i_sm)
            q_sm = q_sm*ax_sm + SM_ERFC_A_Q(i_sm)
         end do
         ex_sm = ex_sm*(p_sm/q_sm)
      else if (ax_sm < 4.0_wp) then
         p_sm = SM_ERFC_B_P(6)
         q_sm = SM_ERFC_B_Q(6)
         do i_sm = 5, 0, -1
            p_sm = p_sm*ax_sm + SM_ERFC_B_P(i_sm)
            q_sm = q_sm*ax_sm + SM_ERFC_B_Q(i_sm)
         end do
         ex_sm = ex_sm*(p_sm/q_sm)
      else
         rax_sm = 1.0_wp/ax_sm
         w_sm = rax_sm*rax_sm
         p_sm = SM_ERFC_C_P(5)
         q_sm = SM_ERFC_C_Q(5)
         do i_sm = 4, 0, -1
            p_sm = p_sm*w_sm + SM_ERFC_C_P(i_sm)
            q_sm = q_sm*w_sm + SM_ERFC_C_Q(i_sm)
         end do
         ex_sm = ((ex_sm*SM_INV_SQRTPI)*rax_sm)*(1.0_wp + w_sm*(p_sm/q_sm))
      end if
      ! no subnormal results (see rdb_sm_exp)
      if (ex_sm < tiny(ex_sm)) ex_sm = 0.0_wp
      if (x_sm < 0.0_wp) then
         y_sm = 2.0_wp - ex_sm
      else
         y_sm = ex_sm
      end if
   end function rdb_sm_erfc
