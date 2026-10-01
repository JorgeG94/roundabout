! -*- f90 -*-  (vim: set ft=fortran: )
!
! rdb_rel_vort_corner.inc — the C-grid corner relative-vorticity RULE, as an
! inlinable body: the circulation-form ζ at one Bu-point (the SW corner of
! T-cell (i,j), position (i-1/2,j-1/2)):
!
!   zeta = ( v_e*dyCv_e - v_w*dyCv_w - (u_n*dxCu_n - u_s*dxCu_s) ) * iareaBu
!
! multiplied by the C1 slip factor `(1-2*ns)*wet_q + 2*ns` (`ns=0` free-slip
! default -> factor=`wet_q`, land corner contributes zero rel-vort; `ns=1`
! no-slip -> factor=`2-wet_q`, MOM6's image vorticity).  On uniform square
! metrics this reduces to `(Δv)/dx - (Δu)/dy` bitwise.
!
! `#include` this file in the `contains` section of any module that computes
! this corner ζ inside a `do concurrent` kernel.  It is the SAME formula
! `rdb_coriolis_adv` evaluates inline for its `q_corner` scratch (Pass 1 of
! the barotropic, Sadourny, HK and Sadourny-energy kernels): a diagnostic
! that reports ζ (the `vorticity_z` derived diag, `rdb_ocean_diag_derived`)
! must agree bit-for-bit with the dynamics' own PV term, so this is a rule,
! not a routine.
!
! Today this file has ONE consumer (the diag).  `rdb_coriolis_adv` is NOT
! yet migrated: its four inline copies are marked "inline twin" and must be
! kept in step by hand until the migration (a separate refactor that needs
! a golden check and GPU timing) lands.
!
! The consuming module must already `use rdb_constants, only: wp`.
! Everything here is `pure` and device-callable (`!$acc routine seq`).
!
! Do NOT add state, logging, or anything impure to this file.

   pure function rdb_rvc_zeta_corner(v_e, v_w, dyCv_e, dyCv_w, &
                                      u_n, u_s, dxCu_n, dxCu_s, &
                                      iareaBu, wet_q, ns) result(zeta)
      !$acc routine seq
      !! Circulation-form relative vorticity ζ = dv/dx - du/dy at one
      !! C-grid (Bu) corner, MOM6's `RV` / roundabout `coriolis_adv`'s
      !! `q_corner` — see this file's header for the formula and the
      !! consumers that must agree on it.
      real(wp), intent(in) :: v_e, v_w, dyCv_e, dyCv_w
         !! v at the east/west v-faces abutting this corner, and their
         !! `dyCv` face lengths (m).
      real(wp), intent(in) :: u_n, u_s, dxCu_n, dxCu_s
         !! u at the north/south u-faces abutting this corner, and their
         !! `dxCu` face lengths (m).
      real(wp), intent(in) :: iareaBu
         !! 1/(corner cell area), m^-2.
      real(wp), intent(in) :: wet_q
         !! Corner (Bu) open mask, 0 (land corner) or 1 (all 4
         !! surrounding T-cells wet).
      real(wp), intent(in) :: ns
         !! Slip selector: 0.0 = free-slip (default), 1.0 = no-slip.
      real(wp) :: zeta

      zeta = ((1.0_wp - 2.0_wp*ns)*wet_q + 2.0_wp*ns)* &
             ((v_e*dyCv_e - v_w*dyCv_w) - (u_n*dxCu_n - u_s*dxCu_s))* &
             iareaBu
   end function rdb_rvc_zeta_corner
