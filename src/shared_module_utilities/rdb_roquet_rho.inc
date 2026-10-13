! -*- f90 -*-  (vim: set ft=fortran: )
!
! rdb_roquet_rho.inc — the VALUE of the Roquet et al. (2015) in-situ-DENSITY
! polynomial (the "polyTEOS10_55t" fit, ΔS = 32 g/kg — independent of, and
! NOT bit-identical to, the specific-volume "polyTEOS10_75t" fit ΔS = 24 g/kg
! that `rdb_roquet_spv.inc` carries), split into its (T, S) part and its
! pressure part, as inlinable bodies. Same structure as `rdb_roquet_spv.inc`
! (mirror it when editing either).
!
! `#include` this file in the `contains` section of any module that
! evaluates Roquet in-situ density inside a `do concurrent` kernel, AFTER
! `#include "rdb_roquet_spv.inc"` in the SAME module -- `rdb_roq_rho_ts_coeffs`
! below calls that file's `rdb_roq_ct_of` for the shared PT->CT conversion
! (see that file's header: it is one of three copies of the TEOS-10
! polynomial, and does not need a fourth). `rdb_eos` includes both (its
! `roquet_rho_ts_coeffs` / `roquet_rho_value` are these bodies), and so do
! the FV-MOM6 in-situ PGF kernels (`rdb_ocean_pressure_force`), whose
! vertical Boole rules evaluate one parcel at five pressures.
!
! Citation: Roquet, F., Madec, G., McDougall, T. J., Barker, P. M. (2015):
! "Accurate polynomial expressions for the density and specific volume of
! seawater using the TEOS-10 standard." Ocean Modelling 90:29-43.
!
!   rho(T, S, p) = (r0 + p*(r1 + p*(r2 + p*r3))) + R00p(p),
!   R00p(p)      = p*(ROQ_RHO_R00 + p*(ROQ_RHO_R01 + ... + p*ROQ_RHO_R05)),
!
! `r0..r3` from `rdb_roq_rho_ts_coeffs` (model variables: potential
! temperature degC, practical salinity PSU -> SR -> CT inside), `rho(p)`
! from `rdb_roq_rho_p`. Term for term the value half of `rdb_eos`'s
! `roquet_rho_point` (same SR / CT conversion, same Horner nesting), so the
! pair reproduces its `rho` exactly. UNLIKE the SpV pair, this fit returns
! density directly -- no `1/SV` reciprocal per pressure point.
!
! The consuming module must already `use rdb_constants, only: wp` and have
! the Roquet coefficient tables in scope: `rdb_eos` declares them; every
! other consumer imports from `rdb_eos` the `RHO*` coefficients below plus
! `ROQ_RHO_R00 .. ROQ_RHO_R05`, `ROQ_RHO_RDELTAS`, and (shared with the SpV
! inc) `ROQ_CP0, ROQ_CT_SFAC, ROQ_R1_S0, ROQ_SR_FACTOR`. Everything here is
! `pure` and device-callable.
!
! Do NOT add state, logging, or anything impure to this file.

   pure subroutine rdb_roq_rho_ts_coeffs(t_roq, s_roq, r0_roq, r1_roq, r2_roq, r3_roq)
      !$acc routine seq
      !! The (T, S)-dependent coefficients of the Roquet in-situ-density
      !! polynomial viewed as a polynomial in PRESSURE. Everything expensive
      !! in the EOS -- the shared PT->CT conversion (`rdb_roq_ct_of`, one
      !! sqrt + the degree-7 PT->CT polynomial) and the ~49-term (zs, zt)
      !! sums -- is here; a caller evaluating ONE parcel at several
      !! pressures pays it once and `rdb_roq_rho_p` per pressure (a degree-3
      !! Horner plus the shared degree-6 reference profile, no division).
      real(wp), intent(in)  :: t_roq
         !! Potential temperature (degC).
      real(wp), intent(in)  :: s_roq
         !! Practical salinity (PSU).
      real(wp), intent(out) :: r0_roq
         !! Pressure-independent part, `rho_ts0 + rho_0s0` (kg/m^3).
      real(wp), intent(out) :: r1_roq, r2_roq, r3_roq
         !! Coefficients of `p`, `p^2`, `p^3` (Pa-powers folded in).

      real(wp) :: zt_roq, zs_roq
      real(wp) :: rts0_roq, r0s0_roq

      ! CT, through the one shared helper (`rdb_roq_ct_of`, `rdb_roquet_spv.inc`,
      ! included before this file in the consuming module). This fit's ΔS
      ! offset (32, vs the SpV fit's 24) is distinct and folded into zs here.
      zt_roq = rdb_roq_ct_of(t_roq, s_roq)
      zs_roq = sqrt(abs(s_roq*ROQ_SR_FACTOR + ROQ_RHO_RDELTAS)*ROQ_R1_S0)

      r3_roq = RHO003 + (zs_roq*RHO103 + zt_roq*RHO013)
      r2_roq = RHO002 + (zs_roq*(RHO102 + zs_roq*RHO202) &
                         + zt_roq*(RHO012 + (zs_roq*RHO112 + zt_roq*RHO022)))
      r1_roq = RHO001 + (zs_roq*(RHO101 + zs_roq*(RHO201 + zs_roq*(RHO301 + zs_roq*RHO401))) &
                         + zt_roq*(RHO011 + (zs_roq*(RHO111 + zs_roq*(RHO211 + zs_roq*RHO311)) &
                                             + zt_roq*(RHO021 + (zs_roq*(RHO121 + zs_roq*RHO221) &
                                                                 + zt_roq*(RHO031 + (zs_roq*RHO131 &
                                                                                     + zt_roq*RHO041)))))))
      rts0_roq = zt_roq*(RHO010 &
                         + (zs_roq*(RHO110 + zs_roq*(RHO210 + zs_roq*(RHO310 + zs_roq*(RHO410 + zs_roq*RHO510)))) &
                            + zt_roq*(RHO020 + (zs_roq*(RHO120 + zs_roq*(RHO220 + zs_roq*(RHO320 + zs_roq*RHO420))) &
                                                + zt_roq*(RHO030 + (zs_roq*(RHO130 + zs_roq*(RHO230 + zs_roq*RHO330)) &
                                                                    + zt_roq*(RHO040 + (zs_roq*(RHO140 + zs_roq*RHO240) &
                                                                                        + zt_roq*(RHO050 &
                                                                                                  + (zs_roq*RHO150 &
                                                                                                     + zt_roq*RHO060))))))))))
      r0s0_roq = RHO000 + zs_roq*(RHO100 + zs_roq*(RHO200 + zs_roq*(RHO300 + zs_roq*(RHO400 &
                                                                                     + zs_roq*(RHO500 + zs_roq*RHO600)))))
      r0_roq = rts0_roq + r0s0_roq
   end subroutine rdb_roq_rho_ts_coeffs

   pure function rdb_roq_rho_p(r0_roq, r1_roq, r2_roq, r3_roq, p_roq) result(rho_roq)
      !$acc routine seq
      !! Roquet in-situ density (kg/m^3) at pressure `p_roq` (Pa) of the
      !! parcel whose (T, S) coefficients `rdb_roq_rho_ts_coeffs` returned:
      !! the degree-3 (T, S)-dependent pressure polynomial plus the (T,
      !! S)-independent degree-6 reference profile `R00p`. No reciprocal --
      !! this fit returns density directly, unlike `rdb_roq_spv_p`'s
      !! specific volume.
      real(wp), intent(in) :: r0_roq, r1_roq, r2_roq, r3_roq
         !! `rdb_roq_rho_ts_coeffs` results.
      real(wp), intent(in) :: p_roq
         !! Pressure (Pa).
      real(wp) :: rho_roq

      rho_roq = (r0_roq + p_roq*(r1_roq + p_roq*(r2_roq + p_roq*r3_roq))) &
                + p_roq*(ROQ_RHO_R00 + p_roq*(ROQ_RHO_R01 + p_roq*(ROQ_RHO_R02 + p_roq*(ROQ_RHO_R03 &
                                                                        + p_roq*(ROQ_RHO_R04 + p_roq*ROQ_RHO_R05)))))
   end function rdb_roq_rho_p
