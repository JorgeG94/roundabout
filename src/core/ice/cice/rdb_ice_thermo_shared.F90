!! Shared thermodynamic kernels of the Bitz-Lipscomb column: the bulk surface
!! heat flux and its Tsf-derivative, the saturation humidity over ice, and the
!! conservative relayering of an enthalpy-like column quantity.
module rdb_ice_thermo_shared
   ! SPDX-License-Identifier: BSD-3-Clause
   !
   ! Ported from CICE/Icepack (CICE Consortium, https://github.com/CICE-Consortium).
   ! Upstream: Icepack 9a1f3e4 (Icepack 1.5.3+14), module icepack_therm_shared
   ! (surface_heat_flux, dsurface_heat_flux_dTsf, adjust_enthalpy).
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
   ! William H. Lipscomb (LANL), C. M. Bitz (UW).
   ! Physics: Bitz and Lipscomb (1999, JGR 104 15669); saturation vapour
   ! pressure over ice as the single-exponential fit used by CICE (Hunke et al.
   ! 2015).
   !
   !! **Surface heat flux** (W/m2, positive INTO the ice):
   !! `fsurf = fsw + eps*flw - eps*sigma*TsfK^4 + shcoef*(potT - TsfK)
   !! + lhcoef*(Qa - Qsfc)`, with `Qsfc = qsat(TsfK)/rhoa` and
   !! `qsat = qqqice*exp(-TTTice/TsfK)`.  `dsurface_heat_flux_dTsf` returns
   !! the analytic derivative of every term.  The exponential is the
   !! reproducible `rdb_sm_exp` (included body), so the value does not depend on
   !! the compiler's libm.
   !!
   !! **`sfc_kind` test hook.**  `ICE_SFC_BULK` (0) is the production bulk
   !! formula.  `ICE_SFC_LINEAR` (1) replaces the whole surface balance by the
   !! linearised upward-positive form `fsurf = -(sf_0 + dsf_dT*Tsf)`, which
   !! makes the steady conduction problem exactly solvable in tests; the
   !! component fluxes and their derivatives are returned as zero.
   !!
   !! **`adjust_enthalpy`.**  Conserving relayering of a per-layer quantity
   !! (enthalpy, salinity) from unequal old layers to equal new layers by
   !! overlap integration; the column total of `h*q` and the thickness are
   !! preserved.  Locals are `ICE_NILYR_MAX`-sized (no allocation).
   !!
   !! Every procedure is `pure` with `!$acc routine seq`; the physical
   !! constants arrive in the `ice_constants_t` bundle.  The error word
   !! (`ierr`, `rdb_ice_errors` bits) is only ever OR-ed into: a non-finite
   !! input or output sets `ICE_ERR_NONFINITE`; so does a surface temperature at
   !! or below absolute zero and a layer count outside `1..ICE_NILYR_MAX` (the
   !! registry has no better-fitting bit and this PR adds none).  A non-finite
   !! `z1`, `z2`, `qn`, `hn` or `hlyr` leaves the column untouched.
   !!
   !! **Constants.**  `ICE_QQQICE`/`ICE_TTTICE` live here, not in
   !! `ice_constants_t`, because no constant set overrides them.  Upstream's
   !! atmosphere routine evaluates `exp(-TTT/TsfK)` with a divide and also needs
   !! the OCEAN pair, so the atmosphere-boundary-layer PR (D8) that reuses
   !! `ice_qsat_surface` differs from upstream by about 4e-15 relative.
   !! `sfc_kind` is a TEST HOOK only and must never become a namelist knob.
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
   use, intrinsic :: iso_fortran_env, only: int32, int64
   use rdb_constants, only: wp
   use rdb_ice_dims, only: ICE_NILYR_MAX
   use rdb_ice_errors, only: ICE_ERR_NONFINITE
   use rdb_ice_params, only: ice_constants_t
   implicit none
   private

   public :: ICE_SFC_BULK, ICE_SFC_LINEAR, ICE_QQQICE, ICE_TTTICE
   public :: surface_heat_flux, dsurface_heat_flux_dTsf
   public :: ice_qsat_surface, ice_dqsat_dTsf
   public :: adjust_enthalpy

   integer, parameter :: ICE_SFC_BULK = 0
      !! `sfc_kind`: Icepack bulk surface balance.
   integer, parameter :: ICE_SFC_LINEAR = 1
      !! `sfc_kind`: linearised test surface `fsurf = -(sf_0 + dsf_dT*Tsf)`.

   real(wp), parameter :: ICE_QQQICE = 11637800.0_wp
      !! Saturation-humidity prefactor over ice [kg/m3] (Icepack `qqqice`).
   real(wp), parameter :: ICE_TTTICE = 5897.8_wp
      !! Saturation-humidity exponent over ice [K] (Icepack `TTTice`).

#include "rdb_safe_math_consts.inc"

contains

#include "rdb_safe_math_core.inc"

   pure function ice_qsat_surface(TsfK, rhoa) result(Qsfc)
      !! Saturated specific humidity [kg/kg] at the ice surface at `TsfK` [K]:
      !! `qqqice*exp(-TTTice/TsfK)/rhoa`.  The caller guarantees `TsfK > 0`.
      !$acc routine seq
      real(wp), intent(in) :: TsfK, rhoa
      real(wp) :: Qsfc
      real(wp) :: qsat, tmpvar

      tmpvar = 1.0_wp/TsfK
      qsat = ICE_QQQICE*rdb_sm_exp(-ICE_TTTICE*tmpvar)
      Qsfc = qsat/rhoa
   end function ice_qsat_surface

   pure function ice_dqsat_dTsf(TsfK, rhoa) result(dQsfc)
      !! `dQsfc/dTsf` [kg/kg/K]: `TTTice/TsfK^2 * qsat/rhoa`.
      !$acc routine seq
      real(wp), intent(in) :: TsfK, rhoa
      real(wp) :: dQsfc
      real(wp) :: qsat, tmpvar

      tmpvar = 1.0_wp/TsfK
      qsat = ICE_QQQICE*rdb_sm_exp(-ICE_TTTICE*tmpvar)
      dQsfc = ICE_TTTICE*tmpvar*tmpvar*(qsat/rhoa)
   end function ice_dqsat_dTsf

   pure subroutine surface_heat_flux(Tsf, fswsfc, rhoa, flw, potT, Qa, shcoef, lhcoef, &
                                     emissivity, c, sfc_kind, sf_0, dsf_dT, &
                                     flwoutn, fsensn, flatn, fsurfn, ierr)
      !! Net surface heat flux into the ice at surface temperature `Tsf` [degC]
      !! and its parts.  `fswsfc` absorbed shortwave, `flw` incoming longwave,
      !! `potT` air potential temperature [K], `Qa` air specific humidity,
      !! `shcoef/lhcoef` bulk transfer coefficients [W/m2/K, W/m2 per kg/kg].
      !! `sf_0, dsf_dT` are read only for `sfc_kind = ICE_SFC_LINEAR`.
      !$acc routine seq
      real(wp), intent(in) :: Tsf, fswsfc, rhoa, flw, potT, Qa, shcoef, lhcoef
      real(wp), intent(in) :: emissivity
      type(ice_constants_t), intent(in) :: c
      integer, intent(in) :: sfc_kind
      real(wp), intent(in) :: sf_0, dsf_dT
      real(wp), intent(out) :: flwoutn, fsensn, flatn, fsurfn
      integer(int32), intent(inout) :: ierr
      real(wp) :: TsfK, Qsfc, flwdabs

      if (sfc_kind == ICE_SFC_LINEAR) then
         flwoutn = 0.0_wp
         fsensn = 0.0_wp
         flatn = 0.0_wp
         fsurfn = -(sf_0 + dsf_dT*Tsf)
      else
         TsfK = Tsf + c%Tffresh
         if (.not. (TsfK > 0.0_wp)) ierr = ior(ierr, ICE_ERR_NONFINITE)
         Qsfc = ice_qsat_surface(TsfK, rhoa)

         flwdabs = emissivity*flw
         flwoutn = -emissivity*c%stefan_boltzmann*TsfK**4

         fsensn = shcoef*(potT - TsfK)
         flatn = lhcoef*(Qa - Qsfc)

         fsurfn = fswsfc + flwdabs + flwoutn + fsensn + flatn
      end if
      if (.not. ieee_is_finite(fsurfn)) ierr = ior(ierr, ICE_ERR_NONFINITE)
   end subroutine surface_heat_flux

   pure subroutine dsurface_heat_flux_dTsf(Tsf, rhoa, shcoef, lhcoef, emissivity, c, &
                                           sfc_kind, dsf_dT, &
                                           dfsurfn_dTsf, dflwoutn_dTsf, dfsensn_dTsf, dflatn_dTsf, ierr)
      !! Derivatives with respect to `Tsf` of the terms of `surface_heat_flux`
      !! (`dfsurfn_dTsf` is their sum).  Under `ICE_SFC_LINEAR` the total is
      !! `-dsf_dT` and the parts are zero.
      !$acc routine seq
      real(wp), intent(in) :: Tsf, rhoa, shcoef, lhcoef, emissivity
      type(ice_constants_t), intent(in) :: c
      integer, intent(in) :: sfc_kind
      real(wp), intent(in) :: dsf_dT
      real(wp), intent(out) :: dfsurfn_dTsf, dflwoutn_dTsf, dfsensn_dTsf, dflatn_dTsf
      integer(int32), intent(inout) :: ierr
      real(wp) :: TsfK, dQsfc_dTsf

      if (sfc_kind == ICE_SFC_LINEAR) then
         dflwoutn_dTsf = 0.0_wp
         dfsensn_dTsf = 0.0_wp
         dflatn_dTsf = 0.0_wp
         dfsurfn_dTsf = -dsf_dT
      else
         TsfK = Tsf + c%Tffresh
         if (.not. (TsfK > 0.0_wp)) ierr = ior(ierr, ICE_ERR_NONFINITE)
         dQsfc_dTsf = ice_dqsat_dTsf(TsfK, rhoa)

         dflwoutn_dTsf = -emissivity*c%stefan_boltzmann*4.0_wp*TsfK**3
         dfsensn_dTsf = -shcoef
         dflatn_dTsf = -lhcoef*dQsfc_dTsf

         dfsurfn_dTsf = dflwoutn_dTsf + dfsensn_dTsf + dflatn_dTsf
      end if
      if (.not. ieee_is_finite(dfsurfn_dTsf)) ierr = ior(ierr, ICE_ERR_NONFINITE)
   end subroutine dsurface_heat_flux_dTsf

   pure subroutine adjust_enthalpy(nlyr, z1, z2, hlyr, hn, qn, c, ierr)
      !! Conserving relayering of the layer quantity `qn(1:nlyr)` (top-down):
      !! `z1(1:nlyr+1)` are the interface depths of the old unequal layers,
      !! `z2(1:nlyr+1)` those of the new equal layers of thickness `hlyr`,
      !! `hn` the total thickness.  For `hn <= puny` the column is
      !! homogenised to its mean (zero for `hn <= 0`).
      !$acc routine seq
      integer, intent(in) :: nlyr
      real(wp), intent(in) :: z1(nlyr + 1), z2(nlyr + 1)
      real(wp), intent(in) :: hlyr, hn
      real(wp), intent(inout) :: qn(nlyr)
      type(ice_constants_t), intent(in) :: c
      integer(int32), intent(inout) :: ierr
      integer :: k, k1, k2
      real(wp) :: hovlp, rhlyr, qtot
      real(wp) :: hq(ICE_NILYR_MAX)

      ! No distinct registry bit fits a layer-count range error: NONFINITE.
      if (nlyr < 1 .or. nlyr > ICE_NILYR_MAX) then
         ierr = ior(ierr, ICE_ERR_NONFINITE)
         return
      end if
      if (.not. (ieee_is_finite(hn) .and. ieee_is_finite(hlyr))) then
         ierr = ior(ierr, ICE_ERR_NONFINITE)
         return
      end if
      do k = 1, nlyr
         if (.not. (ieee_is_finite(qn(k)) .and. ieee_is_finite(z1(k)) .and. ieee_is_finite(z2(k)))) then
            ierr = ior(ierr, ICE_ERR_NONFINITE)
            return
         end if
      end do
      if (.not. (ieee_is_finite(z1(nlyr + 1)) .and. ieee_is_finite(z2(nlyr + 1)))) then
         ierr = ior(ierr, ICE_ERR_NONFINITE)
         return
      end if

      if (hn > c%puny) then
         rhlyr = 1.0_wp/hlyr
         do k2 = 1, nlyr
            hq(k2) = 0.0_wp
         end do
         k1 = 1
         k2 = 1
         do while (k1 <= nlyr .and. k2 <= nlyr)
            hovlp = min(z1(k1 + 1), z2(k2 + 1)) - max(z1(k1), z2(k2))
            hovlp = max(hovlp, 0.0_wp)
            hq(k2) = hq(k2) + hovlp*qn(k1)
            if (z1(k1 + 1) > z2(k2 + 1)) then
               k2 = k2 + 1
            else
               k1 = k1 + 1
            end if
         end do
         do k = 1, nlyr
            qn(k) = hq(k)*rhlyr
         end do
      else
         qtot = 0.0_wp
         do k = 1, nlyr
            qtot = qtot + qn(k)*(z1(k + 1) - z1(k))
         end do
         if (hn > 0.0_wp) then
            do k = 1, nlyr
               qn(k) = qtot/hn
            end do
         else
            do k = 1, nlyr
               qn(k) = 0.0_wp
            end do
         end if
      end if
      do k = 1, nlyr
         if (.not. ieee_is_finite(qn(k))) then
            ierr = ior(ierr, ICE_ERR_NONFINITE)
            exit
         end if
      end do
   end subroutine adjust_enthalpy

end module rdb_ice_thermo_shared
