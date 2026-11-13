!! Flux-bounded barotropic continuity helpers: piecewise-cubic mass-transport
!! closures that bound face transport as velocity grows, plus their derivatives
!! and inverses. Operate on the per-face coefficient packs on
!! `barotropic_workstate_t`. Hallberg & Adcroft (2009).
module rdb_bt_cont_type
   !! `find_uhbt`/`find_duhbt_du` are `pure` + `!$acc routine seq` (callable from
   !! the BT substep `do concurrent`); `uhbt_to_ubt` keeps a Newton iteration and
   !! is host/single-thread only (OBC setup).
   use rdb_constants, only: wp
   use rdb_barotropic_workstate, only: local_BT_cont_u_type, local_BT_cont_v_type
   implicit none
   private

   public :: find_uhbt, find_duhbt_du, uhbt_to_ubt
   public :: find_vhbt, find_dvhbt_dv, vhbt_to_vbt

contains

   pure function find_uhbt(u, BTC) result(uhbt)
      !! Zonal mass transport through a u-face given face velocity `u`.
      !! C¹ continuous in `u` (cubic near zero, linear saturation beyond).
      !$acc routine seq
      real(wp), intent(in) :: u
      type(local_BT_cont_u_type), intent(in) :: BTC
      real(wp) :: uhbt

      if (u == 0.0_wp) then
         uhbt = 0.0_wp
      else if (u < BTC%uBT_EE) then
         uhbt = (u - BTC%uBT_EE)*BTC%FA_u_EE + BTC%uh_EE
      else if (u < 0.0_wp) then
         uhbt = u*(BTC%FA_u_E0 + BTC%uh_crvE*u**2)
      else if (u <= BTC%uBT_WW) then
         uhbt = u*(BTC%FA_u_W0 + BTC%uh_crvW*u**2)
      else
         uhbt = (u - BTC%uBT_WW)*BTC%FA_u_WW + BTC%uh_WW
      end if
   end function find_uhbt

   pure function find_duhbt_du(u, BTC) result(duhbt_du)
      !! Marginal zonal face area `d(uhbt)/du`. At `u = 0` returns the average
      !! of the two cubic-branch slopes (discontinuity harmless — only consumed
      !! via `max(…, h_neglect)`).
      !$acc routine seq
      real(wp), intent(in) :: u
      type(local_BT_cont_u_type), intent(in) :: BTC
      real(wp) :: duhbt_du

      if (u == 0.0_wp) then
         duhbt_du = 0.5_wp*(BTC%FA_u_E0 + BTC%FA_u_W0)
      else if (u < BTC%uBT_EE) then
         duhbt_du = BTC%FA_u_EE
      else if (u < 0.0_wp) then
         duhbt_du = BTC%FA_u_E0 + 3.0_wp*BTC%uh_crvE*u**2
      else if (u <= BTC%uBT_WW) then
         duhbt_du = BTC%FA_u_W0 + 3.0_wp*BTC%uh_crvW*u**2
      else
         duhbt_du = BTC%FA_u_WW
      end if
   end function find_duhbt_du

   pure function uhbt_to_ubt(uhbt, BTC) result(ubt)
      !! Invert `find_uhbt`: recover `u` from a target transport `uhbt`.
      !! Saturated branches close in one line; cubic branches use Newton +
      !! bisection fallback to tol·|uhbt|. Hallberg & Adcroft (2009).
      real(wp), intent(in) :: uhbt
      type(local_BT_cont_u_type), intent(in) :: BTC
      real(wp) :: ubt

      real(wp) :: ubt_min, ubt_max
      real(wp) :: uhbt_err, derr_du
      real(wp) :: uherr_min, uherr_max
      real(wp), parameter :: tol = 1.0e-10_wp
      integer, parameter :: max_itt = 20
      integer :: itt

      if (uhbt == 0.0_wp) then
         ubt = 0.0_wp
      else if (uhbt < BTC%uh_EE) then
         ubt = BTC%uBT_EE + (uhbt - BTC%uh_EE)/BTC%FA_u_EE
      else if (uhbt < 0.0_wp) then
         ubt_min = BTC%uBT_EE
         uherr_min = BTC%uh_EE - uhbt
         ubt_max = 0.0_wp
         uherr_max = -uhbt
         ubt = BTC%uBT_EE*(uhbt/BTC%uh_EE)
         do itt = 1, max_itt
            uhbt_err = ubt*(BTC%FA_u_E0 + BTC%uh_crvE*ubt**2) - uhbt
            if (abs(uhbt_err) < tol*abs(uhbt)) exit
            if (uhbt_err > 0.0_wp) then
               ubt_max = ubt
               uherr_max = uhbt_err
            end if
            if (uhbt_err < 0.0_wp) then
               ubt_min = ubt
               uherr_min = uhbt_err
            end if
            derr_du = BTC%FA_u_E0 + 3.0_wp*BTC%uh_crvE*ubt**2
            if ((uhbt_err >= derr_du*(ubt - ubt_min)) .or. &
                (-uhbt_err >= derr_du*(ubt_max - ubt)) .or. (derr_du <= 0.0_wp)) then
               ubt = ubt_max + (ubt_min - ubt_max)*(uherr_max/(uherr_max - uherr_min))
            else
               ubt = ubt - uhbt_err/derr_du
               if (abs(uhbt_err) < (0.01_wp*tol)*abs(ubt_min*derr_du)) exit
            end if
         end do
      else if (uhbt <= BTC%uh_WW) then
         ubt_min = 0.0_wp
         uherr_min = -uhbt
         ubt_max = BTC%uBT_WW
         uherr_max = BTC%uh_WW - uhbt
         ubt = BTC%uBT_WW*(uhbt/BTC%uh_WW)
         do itt = 1, max_itt
            uhbt_err = ubt*(BTC%FA_u_W0 + BTC%uh_crvW*ubt**2) - uhbt
            if (abs(uhbt_err) < tol*abs(uhbt)) exit
            if (uhbt_err > 0.0_wp) then
               ubt_max = ubt
               uherr_max = uhbt_err
            end if
            if (uhbt_err < 0.0_wp) then
               ubt_min = ubt
               uherr_min = uhbt_err
            end if
            derr_du = BTC%FA_u_W0 + 3.0_wp*BTC%uh_crvW*ubt**2
            if ((uhbt_err >= derr_du*(ubt - ubt_min)) .or. &
                (-uhbt_err >= derr_du*(ubt_max - ubt)) .or. (derr_du <= 0.0_wp)) then
               ubt = ubt_min + (ubt_max - ubt_min)*(-uherr_min/(uherr_max - uherr_min))
            else
               ubt = ubt - uhbt_err/derr_du
               if (abs(uhbt_err) < (0.01_wp*tol)*(ubt_max*derr_du)) exit
            end if
         end do
      else
         ubt = BTC%uBT_WW + (uhbt - BTC%uh_WW)/BTC%FA_u_WW
      end if
   end function uhbt_to_ubt

   pure function find_vhbt(v, BTC) result(vhbt)
      !! Meridional mirror of `find_uhbt`.
      !$acc routine seq
      real(wp), intent(in) :: v
      type(local_BT_cont_v_type), intent(in) :: BTC
      real(wp) :: vhbt

      if (v == 0.0_wp) then
         vhbt = 0.0_wp
      else if (v < BTC%vBT_NN) then
         vhbt = (v - BTC%vBT_NN)*BTC%FA_v_NN + BTC%vh_NN
      else if (v < 0.0_wp) then
         vhbt = v*(BTC%FA_v_N0 + BTC%vh_crvN*v**2)
      else if (v <= BTC%vBT_SS) then
         vhbt = v*(BTC%FA_v_S0 + BTC%vh_crvS*v**2)
      else
         vhbt = (v - BTC%vBT_SS)*BTC%FA_v_SS + BTC%vh_SS
      end if
   end function find_vhbt

   pure function find_dvhbt_dv(v, BTC) result(dvhbt_dv)
      !! Meridional mirror of `find_duhbt_du`.
      !$acc routine seq
      real(wp), intent(in) :: v
      type(local_BT_cont_v_type), intent(in) :: BTC
      real(wp) :: dvhbt_dv

      if (v == 0.0_wp) then
         dvhbt_dv = 0.5_wp*(BTC%FA_v_N0 + BTC%FA_v_S0)
      else if (v < BTC%vBT_NN) then
         dvhbt_dv = BTC%FA_v_NN
      else if (v < 0.0_wp) then
         dvhbt_dv = BTC%FA_v_N0 + 3.0_wp*BTC%vh_crvN*v**2
      else if (v <= BTC%vBT_SS) then
         dvhbt_dv = BTC%FA_v_S0 + 3.0_wp*BTC%vh_crvS*v**2
      else
         dvhbt_dv = BTC%FA_v_SS
      end if
   end function find_dvhbt_dv

   pure function vhbt_to_vbt(vhbt, BTC) result(vbt)
      !! Meridional mirror of `uhbt_to_ubt`.
      real(wp), intent(in) :: vhbt
      type(local_BT_cont_v_type), intent(in) :: BTC
      real(wp) :: vbt

      real(wp) :: vbt_min, vbt_max
      real(wp) :: vhbt_err, derr_dv
      real(wp) :: vherr_min, vherr_max
      real(wp), parameter :: tol = 1.0e-10_wp
      integer, parameter :: max_itt = 20
      integer :: itt

      if (vhbt == 0.0_wp) then
         vbt = 0.0_wp
      else if (vhbt < BTC%vh_NN) then
         vbt = BTC%vBT_NN + (vhbt - BTC%vh_NN)/BTC%FA_v_NN
      else if (vhbt < 0.0_wp) then
         vbt_min = BTC%vBT_NN
         vherr_min = BTC%vh_NN - vhbt
         vbt_max = 0.0_wp
         vherr_max = -vhbt
         vbt = BTC%vBT_NN*(vhbt/BTC%vh_NN)
         do itt = 1, max_itt
            vhbt_err = vbt*(BTC%FA_v_N0 + BTC%vh_crvN*vbt**2) - vhbt
            if (abs(vhbt_err) < tol*abs(vhbt)) exit
            if (vhbt_err > 0.0_wp) then
               vbt_max = vbt
               vherr_max = vhbt_err
            end if
            if (vhbt_err < 0.0_wp) then
               vbt_min = vbt
               vherr_min = vhbt_err
            end if
            derr_dv = BTC%FA_v_N0 + 3.0_wp*BTC%vh_crvN*vbt**2
            if ((vhbt_err >= derr_dv*(vbt - vbt_min)) .or. &
                (-vhbt_err >= derr_dv*(vbt_max - vbt)) .or. (derr_dv <= 0.0_wp)) then
               vbt = vbt_max + (vbt_min - vbt_max)*(vherr_max/(vherr_max - vherr_min))
            else
               vbt = vbt - vhbt_err/derr_dv
               if (abs(vhbt_err) < (0.01_wp*tol)*abs(vbt_min*derr_dv)) exit
            end if
         end do
      else if (vhbt <= BTC%vh_SS) then
         vbt_min = 0.0_wp
         vherr_min = -vhbt
         vbt_max = BTC%vBT_SS
         vherr_max = BTC%vh_SS - vhbt
         vbt = BTC%vBT_SS*(vhbt/BTC%vh_SS)
         do itt = 1, max_itt
            vhbt_err = vbt*(BTC%FA_v_S0 + BTC%vh_crvS*vbt**2) - vhbt
            if (abs(vhbt_err) < tol*abs(vhbt)) exit
            if (vhbt_err > 0.0_wp) then
               vbt_max = vbt
               vherr_max = vhbt_err
            end if
            if (vhbt_err < 0.0_wp) then
               vbt_min = vbt
               vherr_min = vhbt_err
            end if
            derr_dv = BTC%FA_v_S0 + 3.0_wp*BTC%vh_crvS*vbt**2
            if ((vhbt_err >= derr_dv*(vbt - vbt_min)) .or. &
                (-vhbt_err >= derr_dv*(vbt_max - vbt)) .or. (derr_dv <= 0.0_wp)) then
               vbt = vbt_min + (vbt_max - vbt_min)*(-vherr_min/(vherr_max - vherr_min))
            else
               vbt = vbt - vhbt_err/derr_dv
               if (abs(vhbt_err) < (0.01_wp*tol)*(vbt_max*derr_dv)) exit
            end if
         end do
      else
         vbt = BTC%vBT_SS + (vhbt - BTC%vh_SS)/BTC%FA_v_SS
      end if
   end function vhbt_to_vbt

end module rdb_bt_cont_type
