! -*- f90 -*-  (vim: set ft=fortran: )
!
! rdb_roquet_spv.inc — the VALUE of the Roquet et al. (2015) specific-volume
! polynomial, split into its (T, S) part and its pressure part, as inlinable
! bodies.
!
! `#include` this file in the `contains` section of any module that evaluates
! Roquet SpV inside a `do concurrent` kernel.  `rdb_eos` includes it (its
! `roquet_spv_ts_coeffs` / `roquet_spv_value` are these bodies), and so do the
! FV-MOM6 in-situ PGF kernels (`rdb_ocean_pressure_force`), whose vertical Boole
! rules evaluate one parcel at five pressures.  Measured (nvfortran 26.5, cc70):
! the PGF kernels' call to `rdb_eos`'s out-of-line `roquet_spv_ts_coeffs` was
! NOT inlined across the module boundary (a real `call` in the PTX, its four
! results through the stack); the module-local copy this file gives is.
!
!   SV(T, S, p) = (sv0 + p*(sv1 + p*(sv2 + p*sv3))) + SV00p(p),
!   SV00p(p)    = p*(ROQ_V00 + p*(ROQ_V01 + ... + p*ROQ_V05)),
!
! `sv0..sv3` from `rdb_roq_ts_coeffs` (model variables: potential temperature
! degC, practical salinity PSU -> SR -> CT inside), `SV(p)` from
! `rdb_roq_spv_p`.  Term for term the value half of `rdb_eos`'s
! `roquet_spv_point` (same SR / CT conversion, same Horner nesting), so the
! pair reproduces its `sv` exactly.
!
! The consuming module must already `use rdb_constants, only: wp` and have the
! Roquet coefficient table in scope: `rdb_eos` declares it; every other
! consumer imports from `rdb_eos` the `SPV*` coefficients below plus
! `ROQ_CP0, ROQ_CT_SFAC, ROQ_R1_S0, ROQ_RDELTAS, ROQ_SR_FACTOR` and
! `ROQ_V00 .. ROQ_V05`.  Everything here is `pure` and device-callable.
!
! Do NOT add state, logging, or anything impure to this file.

   pure subroutine rdb_roq_ts_coeffs(t_roq, s_roq, sv0_roq, sv1_roq, sv2_roq, sv3_roq)
      !$acc routine seq
      !! The (T, S)-dependent coefficients of the Roquet SpV polynomial viewed
      !! as a polynomial in PRESSURE.  Everything expensive in the EOS -- two
      !! sqrt, the degree-7 PT->CT polynomial and the ~50-term (zs, zt) sums
      !! -- is here; a caller evaluating ONE parcel at several pressures pays
      !! it once and `rdb_roq_spv_p` per pressure.
      real(wp), intent(in)  :: t_roq
         !! Potential temperature (degC).
      real(wp), intent(in)  :: s_roq
         !! Practical salinity (PSU).
      real(wp), intent(out) :: sv0_roq
         !! Pressure-independent part, `sv_ts0 + sv_0s0` (m^3/kg).
      real(wp), intent(out) :: sv1_roq, sv2_roq, sv3_roq
         !! Coefficients of `p`, `p^2`, `p^3` (Pa-powers folded in).

      real(wp) :: zt_roq, zs_roq, x2_roq, xx_roq, yy_roq, hh_roq
      real(wp) :: c0_roq, c1_roq, c2_roq, c3_roq, c4_roq, c5_roq, c6_roq, c7_roq
      real(wp) :: svts0_roq, sv0s0_roq

      ! CT = ct_from_pt(SR, PT) -- identical to `roquet_spv_point`.
      x2_roq = max(ROQ_CT_SFAC*(s_roq*ROQ_SR_FACTOR), 1.0e-20_wp)
      xx_roq = sqrt(x2_roq)
      yy_roq = t_roq*0.025_wp
      c0_roq = 61.01362420681071_wp &
               + x2_roq*(268.5520265845071_wp &
                         + xx_roq*(937.2099110620707_wp &
                                   + xx_roq*(-1687.914374187449_wp + xx_roq*246.9598888781377_wp)))
      c1_roq = 168776.46138048015_wp &
               + x2_roq*(-12019.028203559312_wp &
                         + xx_roq*(588.1802812170108_wp &
                                   + xx_roq*(936.3206544460336_wp + xx_roq*123.59576582457964_wp)))
      c2_roq = -2735.2785605119625_wp &
               + x2_roq*(3734.858026725145_wp &
                         + xx_roq*(248.39476522971285_wp &
                                   + xx_roq*(-942.7827304544439_wp + xx_roq*(-48.5891069025409_wp))))
      c3_roq = 2574.2164453821433_wp &
               + x2_roq*(-2046.7671145057618_wp &
                         + xx_roq*(-3.871557904936333_wp + xx_roq*369.4389437509002_wp))
      c4_roq = -1536.6644434977543_wp &
               + x2_roq*(465.28655623126450_wp + xx_roq*(-2.6268019854268356_wp + xx_roq*(-33.83664947895248_wp)))
      c5_roq = 545.7340497931629_wp &
               + x2_roq*(-0.6370820302831379_wp + xx_roq*(-9.987880382780322_wp))
      c6_roq = -50.91091728474331_wp + x2_roq*(-10.650848542359153_wp)
      c7_roq = -18.30489878927802_wp
      hh_roq = c0_roq + yy_roq*(c1_roq + yy_roq*(c2_roq + yy_roq*(c3_roq + yy_roq*(c4_roq &
                                                                                   + yy_roq*(c5_roq + yy_roq*(c6_roq + yy_roq*c7_roq))))))
      zt_roq = hh_roq/ROQ_CP0
      zs_roq = sqrt(abs(s_roq*ROQ_SR_FACTOR + ROQ_RDELTAS)*ROQ_R1_S0)

      sv3_roq = SPV003 + (zs_roq*SPV103 + zt_roq*SPV013)
      sv2_roq = SPV002 + (zs_roq*(SPV102 + zs_roq*SPV202) &
                          + zt_roq*(SPV012 + (zs_roq*SPV112 + zt_roq*SPV022)))
      sv1_roq = SPV001 + (zs_roq*(SPV101 + zs_roq*(SPV201 + zs_roq*(SPV301 + zs_roq*SPV401))) &
                          + zt_roq*(SPV011 + (zs_roq*(SPV111 + zs_roq*(SPV211 + zs_roq*SPV311)) &
                                              + zt_roq*(SPV021 + (zs_roq*(SPV121 + zs_roq*SPV221) &
                                                                  + zt_roq*(SPV031 + (zs_roq*SPV131 &
                                                                                      + zt_roq*SPV041)))))))
      svts0_roq = zt_roq*(SPV010 &
                          + (zs_roq*(SPV110 + zs_roq*(SPV210 + zs_roq*(SPV310 + zs_roq*(SPV410 + zs_roq*SPV510)))) &
                             + zt_roq*(SPV020 + (zs_roq*(SPV120 + zs_roq*(SPV220 + zs_roq*(SPV320 + zs_roq*SPV420))) &
                                                 + zt_roq*(SPV030 + (zs_roq*(SPV130 + zs_roq*(SPV230 + zs_roq*SPV330)) &
                                                                     + zt_roq*(SPV040 + (zs_roq*(SPV140 + zs_roq*SPV240) &
                                                                                         + zt_roq*(SPV050 &
                                                                                                   + (zs_roq*SPV150 &
                                                                                                      + zt_roq*SPV060))))))))))
      sv0s0_roq = SPV000 + zs_roq*(SPV100 + zs_roq*(SPV200 + zs_roq*(SPV300 + zs_roq*(SPV400 &
                                                                                      + zs_roq*(SPV500 + zs_roq*SPV600)))))
      sv0_roq = svts0_roq + sv0s0_roq
   end subroutine rdb_roq_ts_coeffs

   pure function rdb_roq_spv_p(sv0_roq, sv1_roq, sv2_roq, sv3_roq, p_roq) result(sv_roq)
      !$acc routine seq
      !! Roquet specific volume (m^3/kg) at pressure `p_roq` (Pa) of the parcel
      !! whose (T, S) coefficients `rdb_roq_ts_coeffs` returned: the degree-3
      !! (T, S)-dependent pressure polynomial plus the (T, S)-independent
      !! degree-6 reference profile `SV00p`.
      real(wp), intent(in) :: sv0_roq, sv1_roq, sv2_roq, sv3_roq
         !! `rdb_roq_ts_coeffs` results.
      real(wp), intent(in) :: p_roq
         !! Pressure (Pa).
      real(wp) :: sv_roq

      sv_roq = (sv0_roq + p_roq*(sv1_roq + p_roq*(sv2_roq + p_roq*sv3_roq))) &
               + p_roq*(ROQ_V00 + p_roq*(ROQ_V01 + p_roq*(ROQ_V02 + p_roq*(ROQ_V03 &
                                                                        + p_roq*(ROQ_V04 + p_roq*ROQ_V05)))))
   end function rdb_roq_spv_p
