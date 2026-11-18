!! Atmosphere-ice turbulent boundary layer: Monin-Obukhov similarity over sea
!! ice, the ice branch of the Icepack atmospheric boundary-layer routine.
module rdb_ice_atmo
   ! SPDX-License-Identifier: BSD-3-Clause
   !
   ! Ported from CICE/Icepack (CICE Consortium, https://github.com/CICE-Consortium).
   ! Upstream: Icepack 9a1f3e4 (Icepack 1.5.3+14), module icepack_atmo
   ! (atmo_boundary_layer, ice branch; compute_stability_parameter,
   ! compute_stability_function and the unstable profile functions).
   ! Translated to roundabout's do-concurrent/GPU form and modified by the
   ! roundabout developers; see NOTICE for the full notice.
   !
   ! Copyright 1998-2026, Triad National Security, LLC. All rights reserved.
   ! Redistribution and use in source and binary forms, with or without
   ! modification, are permitted under the BSD-3-Clause conditions reproduced in
   ! NOTICE. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS
   ! AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES ARE
   ! DISCLAIMED (full disclaimer in NOTICE).
   !
   ! Original authors (from the upstream headers): Elizabeth C. Hunke (LANL),
   ! William H. Lipscomb (LANL), Bruce P. Briegleb (NCAR), Michael Winton (GFDL),
   ! Nicole Jeffery (LANL).
   ! Physics: Jordan, Andreas and Makshtas (1999, JGR 104 7785) for the ice
   ! surface layer and the windless sensible-heat term; stable profile
   ! functions Holtslag and de Bruin (1988, J. Appl. Meteorol. 27 689); unstable
   ! profile functions Paulson (1970, J. Appl. Meteorol. 9 857); constant-
   ! roughness neutral drag as in Large and Yeager (2004, NCAR/TN-460+STR);
   ! high-frequency momentum coupling Roberts et al. (2014, Ocean Modelling
   ! 79 28).
   !
   !! **What it computes.**  For one ice surface (one cell and category) at the
   !! start-of-step surface temperature `Tsf`: the bulk heat and moisture
   !! transfer coefficients `shcoef` [W/m2/K] and `lhcoef` [W/m2 per kg/kg], the
   !! wind stress on the ice `strx, stry` [N/m2, on the model axes] and the
   !! potential-temperature and humidity differences `delt` [K], `delq`
   !! [kg/kg].  All fluxes are positive downward.  The thermodynamics freezes
   !! the coefficients and evaluates the surface fluxes itself
   !! (`rdb_ice_thermo_shared`).
   !!
   !! **Algorithm** (neutral first guess, then `natmiter` Monin-Obukhov passes
   !! in upstream's order: stability parameter at the momentum and scalar
   !! heights, clamped to `|z/L| <= 10`; profile functions; shifted exchange
   !! coefficients `rd, rh, re`; updated `ustar, tstar, qstar`).  The
   !! iteration stops early when `|ustar - ustar_prev|/ustar <= atmiter_conv`.
   !! The neutral drag is the constant-roughness `rdn = kappa/ln(zref/iceruf)`,
   !! derived once per configuration by `ice_atm_params_finalize` into `atm%rdn`
   !! (the kernel does no `iceruf` range check).
   !! With `highfreq` the wind-speed floor is 0.5 m/s and the stress is
   !! `rhoa*rd^2*|Uatm - Uice|*(Uatm - Uice)`; otherwise the floor is 1 m/s and
   !! the stress is `rhoa*ustar*rd*Uatm`.  The sensible coefficient carries the
   !! windless `+1 W/m2/K` of Jordan et al. (1999).
   !!
   !! **Kept upstream quirks.**  (1) The stable branch is used at exactly
   !! `z/L = 0` and `psi_stable(0) = +0.025` (not 0), while the unstable
   !! branch tends to 0 as `z/L -> 0-`, so the "neutral" drag is 0.5 percent
   !! above `rdn^2` and the coefficient jumps across `delt = 0`.  (2) The stress
   !! is not normalised by the floored wind speed (`strx = rhoa*ustar*rd*uatm`).
   !! (3) The surface saturation humidity is `qqqice*exp(-TTTice/TsfK)` with a
   !! divide (the thermodynamics routine multiplies by a reciprocal).
   !!
   !! **Differences from upstream.**  Ice branch only (the open-ocean branch is
   !! the bulk-forcing module's); no optional arguments, no module globals;
   !! `atmbndy = 'similarity'` only; no form drag; no isotopes; the 2 m
   !! reference diagnostics `Tref, Qref, Uref` are dropped.  Every
   !! transcendental is the reproducible `rdb_sm_*` body, so results do not
   !! depend on the compiler's libm.  The ice and the open-ocean schemes differ
   !! on purpose (constant roughness and Holtslag-de Bruin stable functions here;
   !! polynomial drag and Businger-Dyer there); do not share code.
   !!
   !! **Caveats for callers and later tests.**
   !! (a) Near neutral the stable/unstable branch switch is a jump: `psi` goes
   !! from 0 to +0.025 at `z/L = 0`.  `-0` is unreachable (`z/L` is a clamped
   !! product), but FMA contraction on the GPU build against none on the CPU build
   !! can put `|z/L| ~ 1e-20` on different sides of zero, a ~0.5 percent stress
   !! difference; device-parity tolerances of later PRs must allow for it.
   !! (b) Five passes do not converge stable cases: the stress at `delt = +3 K` is
   !! 4.4e-4 relative from the 50-pass answer (as upstream).
   !! (c) Below `zlvl ~ 7e-3 m` the drag denominator reaches zero (as upstream).
   !! (d) The input-prep PR must supply: `potT` and `rhoa` by the CICE fallback
   !! formulas (`potT = Tair*(1e5/p)^0.286`, `rhoa = p/(287.058 (1 + 0.608 Qa) Tair)`);
   !! `zlvl = z_wind`, `zlvs = 10` fixed; winds rotated to the model axes with
   !! `wind = hypot(uatm, vatm)`; under `highfreq` (OM3 uses `.true.`) the
   !! T-averaged start-of-step ice velocity; per-category start-of-step `Tsfc`; the
   !! `aicen_init > puny` gate; and the stress merged with weight `aicen_init`.
   !! The config PR must refuse `atmbndy /= similarity`, `formdrag`,
   !! `calc_strair = .false.`, `natmiter < 0`, `iceruf` outside `(0, zref)` and
   !! `zlvl`/`zlvs <= iceruf`, and call `ice_atm_params_finalize` after reading
   !! `iceruf`.
   !!
   !! **Errors.**  A non-finite input, a non-positive `TsfK, potT, rhoa, zlvl,
   !! zlvs`, an unset `atm%rdn`, or a non-finite output sets `ICE_ERR_NONFINITE` in `ierr` (OR-ed,
   !! never cleared) and every output is returned as zero (defined, never a
   !! laundered clamp).  `natmiter <= 0` leaves the neutral coefficients.
   !!
   !! Every procedure is `pure` with `!$acc routine seq`; scalars in, scalars out.
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
   use, intrinsic :: iso_fortran_env, only: int32, int64
   use rdb_constants, only: wp
   use rdb_ice_errors, only: ICE_ERR_NONFINITE
   use rdb_ice_params, only: ice_constants_t, ice_atm_params_t
   use rdb_ice_thermo_shared, only: ICE_QQQICE, ICE_TTTICE
   implicit none
   private

   public :: ice_psi_m, ice_psi_h, ice_psi_stable, ice_stability_parameter
   public :: ice_atm_boundary_ice
   public :: ICE_ATMO_HOL_MAX

   real(wp), parameter :: ICE_ATMO_HOL_MAX = 10.0_wp
      !! Clamp on the stability parameter `|z/L|`.
   real(wp), parameter :: PIH = 3.14159265358979323846_wp*0.5_wp
      !! pi/2 (upstream `pih = pi*p5`).

#include "rdb_safe_math_consts.inc"

contains

#include "rdb_safe_math_core.inc"
#include "rdb_safe_math_atan.inc"

   pure function ice_psi_stable(hol) result(psi)
      !! Stable profile function (momentum and scalar alike), Holtslag and de
      !! Bruin (1988): `-(0.7 h + 0.75 (h - 14.3) exp(-0.35 h) + 10.7)`.
      !! Note `psi_stable(0) = +0.025`.
      !$acc routine seq
      real(wp), intent(in) :: hol
      real(wp) :: psi

      psi = -(0.7_wp*hol + 0.75_wp*(hol - 14.3_wp)*rdb_sm_exp(-0.35_wp*hol) + 10.7_wp)
   end function ice_psi_stable

   pure function ice_psi_m(hol) result(psi)
      !! Momentum stability function at `hol = z/L`: the stable form for
      !! `hol >= +0` (sign test as upstream, so `-0` is unstable) and Paulson's
      !! `ln((1+X(2+X))(1+X^2)/8) - 2 atan X + pi/2` otherwise, with
      !! `X = sqrt(max(sqrt|1 - 16 hol|, 1))`.
      !$acc routine seq
      real(wp), intent(in) :: hol
      real(wp) :: psi
      real(wp) :: stable, xd

      stable = 0.5_wp + sign(0.5_wp, hol)
      if (stable > 0.5_wp) then
         psi = ice_psi_stable(hol)
      else
         xd = sqrt(max(sqrt(abs(1.0_wp - 16.0_wp*hol)), 1.0_wp))
         psi = rdb_sm_log((1.0_wp + xd*(2.0_wp + xd))*(1.0_wp + xd*xd)/8.0_wp) &
               - 2.0_wp*rdb_sm_atan(xd) + PIH
      end if
   end function ice_psi_m

   pure function ice_psi_h(hol) result(psi)
      !! Scalar (heat and moisture) stability function: the stable form for
      !! `hol >= +0`, `2 ln((1 + X^2)/2)` otherwise.
      !$acc routine seq
      real(wp), intent(in) :: hol
      real(wp) :: psi
      real(wp) :: stable, xd

      stable = 0.5_wp + sign(0.5_wp, hol)
      if (stable > 0.5_wp) then
         psi = ice_psi_stable(hol)
      else
         xd = sqrt(max(sqrt(abs(1.0_wp - 16.0_wp*hol)), 1.0_wp))
         psi = 2.0_wp*rdb_sm_log((1.0_wp + xd*xd)/2.0_wp)
      end if
   end function ice_psi_h

   pure function ice_stability_parameter(zlvl, thva, ustar, tstar, qstar, qden, c) result(hol)
      !! `z/L = kappa g z (tstar/thva + qstar/qden)/ustar^2` with
      !! `qden = 1/zvir + Qa` (formed once per call by the caller), clamped
      !! to `[-10, 10]`.  The caller guarantees finite inputs and `ustar > 0`
      !! (the clamp is a min/max select that would launder a NaN).
      !$acc routine seq
      real(wp), intent(in) :: zlvl, thva, ustar, tstar, qstar, qden
      type(ice_constants_t), intent(in) :: c
      real(wp) :: hol

      hol = c%vonkar*c%gravit*zlvl*(tstar/thva + qstar/qden)/ustar**2
      hol = sign(min(abs(hol), ICE_ATMO_HOL_MAX), hol)
   end function ice_stability_parameter

   pure subroutine ice_atm_boundary_ice(Tsf, potT, uatm, vatm, uvel, vvel, wind, zlvl, zlvs, Qa, rhoa, &
                                        c, atm, shcoef, lhcoef, strx, stry, delt, delq, ierr)
      !! Exchange coefficients, wind stress and the air-surface differences for
      !! one ice surface.  `Tsf` [degC] start-of-step surface temperature,
      !! `potT` [K], `uatm, vatm` [m/s] wind on the model axes, `uvel, vvel`
      !! [m/s] ice velocity at the cell (read only with `atm%highfreq`),
      !! `wind` [m/s] wind speed (read only without `highfreq`), `zlvl, zlvs`
      !! [m] momentum and scalar reference heights, `Qa` [kg/kg], `rhoa`
      !! [kg/m3].  Outputs as in the module header.
      !$acc routine seq
      real(wp), intent(in) :: Tsf, potT, uatm, vatm, uvel, vvel, wind, zlvl, zlvs, Qa, rhoa
      type(ice_constants_t), intent(in) :: c
      type(ice_atm_params_t), intent(in) :: atm
      real(wp), intent(out) :: shcoef, lhcoef, strx, stry, delt, delq
      integer(int32), intent(inout) :: ierr
      integer :: k
      real(wp) :: qden, TsfK, qsat, ssq, thva, alzm, alzs, cpvir, cp, umin, vmag, du, dv, dumag
      real(wp) :: rdn, rhn, ren, rd, rh, re, ustar, ustar_prev, tstar, qstar, holm, hols, tau
      logical :: ok

      shcoef = 0.0_wp
      lhcoef = 0.0_wp
      strx = 0.0_wp
      stry = 0.0_wp
      delt = 0.0_wp
      delq = 0.0_wp

      ok = ieee_is_finite(Tsf) .and. ieee_is_finite(potT) .and. ieee_is_finite(uatm) &
           .and. ieee_is_finite(vatm) .and. ieee_is_finite(Qa) .and. ieee_is_finite(rhoa) &
           .and. ieee_is_finite(zlvl) .and. ieee_is_finite(zlvs)
      if (atm%highfreq) then
         ok = ok .and. ieee_is_finite(uvel) .and. ieee_is_finite(vvel)
      else
         ok = ok .and. ieee_is_finite(wind)
      end if
      if (ok) then
         TsfK = Tsf + c%Tffresh
         ok = TsfK > 0.0_wp .and. potT > 0.0_wp .and. rhoa > 0.0_wp .and. zlvl > 0.0_wp &
              .and. zlvs > 0.0_wp .and. atm%rdn > 0.0_wp
      end if
      if (.not. ok) then
         ierr = ior(ierr, ICE_ERR_NONFINITE)
         return
      end if

      cpvir = c%cp_wv/c%cp_air - 1.0_wp
      if (atm%highfreq) then
         umin = 0.5_wp
         du = uatm - uvel
         dv = vatm - vvel
         dumag = sqrt(du**2 + dv**2)
         vmag = max(umin, dumag)
      else
         umin = 1.0_wp
         du = uatm
         dv = vatm
         dumag = 0.0_wp
         vmag = max(umin, wind)
      end if
      rdn = atm%rdn   ! neutral coefficient, once per configuration (ice_atm_params_finalize)

      delt = potT - TsfK
      qsat = ICE_QQQICE*rdb_sm_exp(-ICE_TTTICE/TsfK)
      ssq = qsat/rhoa
      thva = potT*(1.0_wp + c%zvir*Qa)
      delq = Qa - ssq
      alzm = rdb_sm_log(zlvl/c%zref)
      alzs = rdb_sm_log(zlvs/c%zref)
      cp = c%cp_air*(1.0_wp + cpvir*ssq)
      qden = 1.0_wp/c%zvir + Qa

      ! Neutral first estimate.
      rhn = rdn
      ren = rdn
      rd = rdn
      rh = rhn
      re = ren
      ustar = rdn*vmag
      tstar = rhn*delt
      qstar = ren*delq
      ustar_prev = 2.0_wp*ustar

      do k = 1, atm%natmiter
         if (.not. (abs(ustar - ustar_prev)/ustar > atm%atmiter_conv)) exit
         ustar_prev = ustar
         holm = ice_stability_parameter(zlvl, thva, ustar, tstar, qstar, qden, c)
         hols = ice_stability_parameter(zlvs, thva, ustar, tstar, qstar, qden, c)
         rd = rdn/(1.0_wp + rdn/c%vonkar*(alzm - ice_psi_m(holm)))
         rh = rhn/(1.0_wp + rhn/c%vonkar*(alzs - ice_psi_h(hols)))
         re = ren/(1.0_wp + ren/c%vonkar*(alzs - ice_psi_h(hols)))
         ustar = rd*vmag
         tstar = rh*delt
         qstar = re*delq
      end do

      if (atm%highfreq) then
         tau = rhoa*rd*rd
         strx = tau*dumag*du
         stry = tau*dumag*dv
      else
         tau = rhoa*ustar*rd
         strx = tau*uatm
         stry = tau*vatm
      end if

      shcoef = rhoa*ustar*cp*rh + 1.0_wp
      lhcoef = rhoa*ustar*c%Lsub*re

      if (.not. (ieee_is_finite(shcoef) .and. ieee_is_finite(lhcoef) .and. ieee_is_finite(strx) &
                 .and. ieee_is_finite(stry) .and. ieee_is_finite(delt) .and. ieee_is_finite(delq))) then
         ierr = ior(ierr, ICE_ERR_NONFINITE)
         shcoef = 0.0_wp
         lhcoef = 0.0_wp
         strx = 0.0_wp
         stry = 0.0_wp
         delt = 0.0_wp
         delq = 0.0_wp
      end if
   end subroutine ice_atm_boundary_ice

end module rdb_ice_atmo
