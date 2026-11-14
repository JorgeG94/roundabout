!! CCSM3 sea-ice shortwave: category albedos and absorbed solar radiation (scalar kernels).
!!
!! **What this is.**  The two column-local arithmetic kernels of the CCSM3 shortwave scheme
!! (Briegleb et al. 2004, NCAR/TN-463+STR): `ice_albedos_ccsm3` (thickness- and
!! melt-dependent ice and snow albedos in the visible and near-IR bands) and
!! `ice_absorbed_ccsm3` (the absorbed / penetrating / transmitted split and the Beer-law
!! layer absorption).  They are `pure`, carry `!$acc routine seq`, hold no state and are
!! callable from inside a `do concurrent` kernel.  The drivers, state and aggregation are
!! a later PR.
!!
! SPDX-License-Identifier: BSD-3-Clause
!
! Ported from CICE/Icepack (CICE Consortium, https://github.com/CICE-Consortium).
! Upstream: Icepack 9a1f3e4 (Icepack 1.5.3+14), module icepack_shortwave
! (routines compute_albedos, absorbed_solar; the ccsm3 branch of
! icepack_step_radiation).
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
! Original authors (from the upstream module): Bruce P. Briegleb (NCAR),
! Elizabeth C. Hunke (LANL), William H. Lipscomb (LANL).
! Physics: Briegleb et al. (2004, NCAR/TN-463+STR); Maykut and Untersteiner
! (1971, JGR 76, 1550); Bitz and Lipscomb (1999, JGR 104, 15669).
module rdb_ice_shortwave
   !! Scalar CCSM3 albedo and absorbed-solar kernels.
   !!
   !! **Physics.**  Per ice category, with `hi` the mean ice thickness and `hs` the mean
   !! snow thickness over the ice-covered area:
   !!
   !! * bare ice: `fh = min(atan(4 hi)/atan(4 ahmax), 1)`, `albedo = albice*fh + albocn*(1-fh)`,
   !!   reduced by `dalb_mlt*fT` with `fT = min(-Tsfc/dT_melt - 1, 0)` (zero below -1 degC,
   !!   linear to -1 at 0 degC) and floored at the ocean albedo;
   !! * snow (only if `hs > puny`): cold snow albedo less `dalb_mltv*fT` (visible) or
   !!   `dalb_mlti*fT` (near-IR); snow cover fraction `asnow = hs/(hs + snowpatch)`;
   !! * direct and diffuse albedos are identical in this scheme, so only visible and
   !!   near-IR albedos exist;
   !! * all near-IR and the snow-covered / reflected-out part of the visible is absorbed
   !!   at the surface; a fraction `i0vis` of the bare-ice visible absorbed flux penetrates
   !!   and is attenuated by Beer's law with extinction `kappav` (Maykut and Untersteiner
   !!   1971); nothing is absorbed inside the snow.
   !!
   !! **Energy identity.**  `SW_in = reflected + fswsfc + fswint + fswthru` and
   !! `sum_k Iswabs(k) = fswint`, by construction (telescoping layer transmittances).
   !!
   !! **Conventions.**  Temperatures in degC (`Timelt = 0`), fluxes in W/m2 over the
   !! ice-covered area of the category, layer index `k = 1` is the TOP ice layer (Icepack
   !! order, the documented exception to the ocean's bottom-up rule).  The caller guards
   !! `aicen > puny` (an empty category has albedo `ICE_ALBOCN` and zero fluxes) and
   !! supplies finite `Tsfc`.
   !!
   !! **Transcendentals.**  `exp` and `atan` come from the safe-math included bodies
   !! (`rdb_sm_exp`, `rdb_sm_atan`).  `atan(4 ahmax)` is per configuration, not per cell:
   !! compute it ONCE with `ice_albedo_fhtan` and pass it down as `fhtan`.  `ahmax = 0` divides
   !! by zero here exactly as upstream; refusing it is owed by the config PR (D1.3, `validate_config`).
   !!
   !! **Band split of the transmitted light.**  `fswthru` is all visible; its direct/diffuse
   !! parts are recoverable as `fswthru*swvdr/SW_v` and `fswthru*swvdf/SW_v` (identical albedos).
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_constants, only: wp
   implicit none
   private

   public :: ice_albedo_fhtan, ice_albedos_ccsm3, ice_absorbed_ccsm3
   public :: ICE_ALBOCN, ICE_KAPPAV, ICE_I0VIS, ICE_SNOWPATCH_DEFAULT, ICE_DT_MELT
   public :: ICE_DALB_MLT, ICE_DALB_MLTV, ICE_DALB_MLTI, ICE_SW_PUNY, ICE_SW_TIMELT
   public :: ICE_ALBICEV_DEFAULT, ICE_ALBICEI_DEFAULT, ICE_ALBSNOWV_DEFAULT
   public :: ICE_ALBSNOWI_DEFAULT, ICE_AHMAX_DEFAULT

   ! NOTE: the constants below are local to this module until the shared ice parameter
   ! module (D1.1, `rdb_ice_params`) lands; it then supersedes them and these become
   ! `use` associations.  Values: Icepack parameter defaults / CCSM3 (Briegleb et al. 2004).
   real(wp), parameter :: ICE_ALBOCN = 0.06_wp
      !! Open-ocean albedo, and the floor of every ice albedo (fraction).
   real(wp), parameter :: ICE_KAPPAV = 1.4_wp
      !! Visible extinction coefficient in ice, wavelength < 700 nm (1/m).
   real(wp), parameter :: ICE_I0VIS = 0.70_wp
      !! Fraction of the bare-ice absorbed visible flux that penetrates the surface.
   real(wp), parameter :: ICE_SNOWPATCH_DEFAULT = 0.02_wp
      !! Default snow-cover fraction length scale, `asnow = hs/(hs + snowpatch)` (m).  The
      !! kernel takes it as an argument: 0.02 is the upstream default, 0.005 the CESM-coupled set.
   real(wp), parameter :: ICE_DT_MELT = 1.0_wp
      !! Temperature range below the melt point over which the melt albedo reduction ramps (degC).
   real(wp), parameter :: ICE_DALB_MLT = -0.075_wp
      !! Ice albedo change per `ICE_DT_MELT` of warming, both bands (fraction).
   real(wp), parameter :: ICE_DALB_MLTV = -0.1_wp
      !! Snow visible albedo change per `ICE_DT_MELT` of warming (fraction).
   real(wp), parameter :: ICE_DALB_MLTI = -0.15_wp
      !! Snow near-IR albedo change per `ICE_DT_MELT` of warming (fraction).
   real(wp), parameter :: ICE_SW_PUNY = 1.0e-11_wp
      !! Threshold on snow thickness (and area) below which snow is absent.
   real(wp), parameter :: ICE_SW_TIMELT = 0.0_wp
      !! Melting temperature of the ice top surface (degC).

   ! Upstream code defaults for the namelist knobs (OM3 uses 0.86/0.44/0.98/0.70/0.1).
   real(wp), parameter :: ICE_ALBICEV_DEFAULT = 0.78_wp
      !! Cold bare thick ice albedo, visible.
   real(wp), parameter :: ICE_ALBICEI_DEFAULT = 0.36_wp
      !! Cold bare thick ice albedo, near-IR.
   real(wp), parameter :: ICE_ALBSNOWV_DEFAULT = 0.98_wp
      !! Cold snow albedo, visible.
   real(wp), parameter :: ICE_ALBSNOWI_DEFAULT = 0.70_wp
      !! Cold snow albedo, near-IR.
   real(wp), parameter :: ICE_AHMAX_DEFAULT = 0.3_wp
      !! Ice thickness above which the ice albedo is constant (m).

#include "rdb_safe_math_consts.inc"

contains

#include "rdb_safe_math_core.inc"
#include "rdb_safe_math_atan.inc"

   pure function ice_albedo_fhtan(ahmax) result(fhtan)
      !! `atan(4 ahmax)`, the normaliser of the thin-ice albedo ramp.  Constant per
      !! configuration: evaluate once on the host and pass it to `ice_albedos_ccsm3`.
      !$acc routine seq
      real(wp), intent(in) :: ahmax
         !! Thickness above which the ice albedo is constant (m).
      real(wp) :: fhtan
      fhtan = rdb_sm_atan(4.0_wp*ahmax)
   end function ice_albedo_fhtan

   pure subroutine ice_albedos_ccsm3(hi, hs, tsfc, albicev, albicei, albsnowv, albsnowi, fhtan, snowpatch, &
                                     alv_i, ali_i, alv_s, ali_s, asnow, alv, ali)
      !! Visible / near-IR albedos of the bare-ice part, the snow part and the category
      !! mean (direct == diffuse).  Briegleb et al. (2004).
      !$acc routine seq
      real(wp), intent(in) :: hi
         !! Mean ice thickness over the ice area (m).
      real(wp), intent(in) :: hs
         !! Mean snow thickness over the ice area (m).
      real(wp), intent(in) :: tsfc
         !! Surface temperature (degC).
      real(wp), intent(in) :: albicev, albicei
         !! Cold bare thick ice albedo, visible / near-IR.
      real(wp), intent(in) :: albsnowv, albsnowi
         !! Cold snow albedo, visible / near-IR.
      real(wp), intent(in) :: fhtan
         !! `ice_albedo_fhtan(ahmax)`.
      real(wp), intent(in) :: snowpatch
         !! Snow-cover length scale (m); `ICE_SNOWPATCH_DEFAULT` = 0.02, CESM-coupled 0.005.
      real(wp), intent(out) :: alv_i, ali_i
         !! Bare-ice albedo, visible / near-IR.
      real(wp), intent(out) :: alv_s, ali_s
         !! Snow albedo, visible / near-IR (`ICE_ALBOCN` when there is no snow).
      real(wp), intent(out) :: asnow
         !! Fractional snow cover.
      real(wp), intent(out) :: alv, ali
         !! Category-mean albedo, visible / near-IR.
      real(wp) :: fh, albo, ft
      fh = min(rdb_sm_atan(4.0_wp*hi)/fhtan, 1.0_wp)
      albo = ICE_ALBOCN*(1.0_wp - fh)
      alv_i = albicev*fh + albo
      ali_i = albicei*fh + albo
      ! temperature dependence: 0 at or below -dT_melt, -1 at the melt point
      ft = min((ICE_SW_TIMELT - tsfc)/ICE_DT_MELT - 1.0_wp, 0.0_wp)
      alv_i = alv_i - ICE_DALB_MLT*ft
      ali_i = ali_i - ICE_DALB_MLT*ft
      ! no sub-ocean albedo for thin, bare, melting ice
      alv_i = max(alv_i, ICE_ALBOCN)
      ali_i = max(ali_i, ICE_ALBOCN)
      if (hs > ICE_SW_PUNY) then
         alv_s = albsnowv - ICE_DALB_MLTV*ft
         ali_s = albsnowi - ICE_DALB_MLTI*ft
         asnow = hs/(hs + snowpatch)
      else
         alv_s = ICE_ALBOCN
         ali_s = ICE_ALBOCN
         asnow = 0.0_wp
      end if
      alv = alv_i*(1.0_wp - asnow) + alv_s*asnow
      ali = ali_i*(1.0_wp - asnow) + ali_s*asnow
   end subroutine ice_albedos_ccsm3

   pure subroutine ice_absorbed_ccsm3(nilyr, hi, asnow, sw_v, sw_n, alv_i, ali_i, alv_s, ali_s, &
                                      fswsfc, fswint, fswthru, iswabs)
      !! Shortwave absorbed at the surface, in the ice interior (per layer) and
      !! transmitted to the ocean.  Only visible light penetrates; the snow layer absorbs
      !! nothing internally (`Sswabs = 0` for this scheme, not returned).
      !$acc routine seq
      integer, intent(in) :: nilyr
         !! Number of ice layers (layer 1 is the top).
      real(wp), intent(in) :: hi
         !! Mean ice thickness over the ice area (m).
      real(wp), intent(in) :: asnow
         !! Fractional snow cover from `ice_albedos_ccsm3`.
      real(wp), intent(in) :: sw_v, sw_n
         !! Incident shortwave, visible / near-IR (direct + diffuse), W/m2.
      real(wp), intent(in) :: alv_i, ali_i, alv_s, ali_s
         !! Albedos from `ice_albedos_ccsm3`.
      real(wp), intent(out) :: fswsfc
         !! Absorbed at the ice/snow surface (W/m2).
      real(wp), intent(out) :: fswint
         !! Absorbed in the ice interior (W/m2).
      real(wp), intent(out) :: fswthru
         !! Transmitted through the ice into the ocean (W/m2, all visible).
      real(wp), intent(out) :: iswabs(nilyr)
         !! Absorbed in ice layer k, k = 1 at the top (W/m2).
      real(wp) :: swabs, fswpen, hilyr, trantop, tranbot
      integer :: k
      swabs = sw_v*((1.0_wp - alv_i)*(1.0_wp - asnow) + (1.0_wp - alv_s)*asnow) &
              + sw_n*((1.0_wp - ali_i)*(1.0_wp - asnow) + (1.0_wp - ali_s)*asnow)
      fswpen = sw_v*(1.0_wp - alv_i)*(1.0_wp - asnow)*ICE_I0VIS
      fswsfc = swabs - fswpen
      hilyr = hi/real(nilyr, wp)
      trantop = 1.0_wp
      tranbot = 1.0_wp
      do k = 1, nilyr
         tranbot = rdb_sm_exp(-ICE_KAPPAV*hilyr*real(k, wp))
         iswabs(k) = fswpen*(trantop - tranbot)
         trantop = tranbot
      end do
      fswthru = fswpen*tranbot
      fswint = fswpen - fswthru
   end subroutine ice_absorbed_ccsm3

end module rdb_ice_shortwave
