!! Column thickness bookkeeping of the Bitz-Lipscomb thermodynamics: snow-ice
!! formation at the waterline (`ice_freeboard`), the conservative return to
!! equal layer thicknesses (`ice_relayer`) and the conversion of the column
!! thicknesses back to the category volumes (`ice_update_state_vthermo`).
!! The thickness-change driver of the next PR joins this module.
module rdb_ice_thickness
   ! SPDX-License-Identifier: BSD-3-Clause
   !
   ! Ported from CICE/Icepack (CICE Consortium, https://github.com/CICE-Consortium).
   ! Upstream: Icepack 9a1f3e4 (Icepack 1.5.3+14), module icepack_therm_vertical
   ! (freeboard, update_state_vthermo, and the relayering block of
   ! thickness_changes).
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
   ! Original authors (from the upstream headers): William H. Lipscomb (LANL),
   ! Elizabeth C. Hunke (LANL), C. M. Bitz (UW).
   ! Physics: Bitz and Lipscomb (1999, JGR 104 15669); Maykut and McPhee (1995,
   ! JGR 100 24691).
   !
   !! **Layer convention.**  Ice and snow layers are indexed TOP-DOWN (k = 1 is
   !! the top), the Icepack order; `nilyr <= ICE_NILYR_MAX`, `nslyr <=
   !! ICE_NSLYR_MAX`.  All arguments are explicit-shape; locals are sized by the
   !! `ICE_*_MAX` bounds (no allocation).  Every procedure is `pure` with
   !! `!$acc routine seq`, takes the constants in the `intent(in)` bundle and
   !! ORs its error bits into `ierr`.  A non-finite input sets
   !! `ICE_ERR_NONFINITE` and leaves every state argument untouched (results
   !! are formed in locals and committed only when finite); the pure outputs
   !! `hilyr`/`hslyr` of `ice_relayer` are then set to 0.
   !!
   !! **Freeboard (snow-ice formation).**  With ice thickness `hin`, snow `hsn`
   !! and densities `rhoi, rhos, rhow`, the snow sits below the waterline when
   !! `wk1 = hsn - hin*(rhow-rhoi)/rhos > puny` (and `hsn > puny`).  Then
   !! `dhsn = min(wk1*rhoi/rhow, hsn)` of snow is removed, layer by layer from
   !! the LAST snow layer up (upstream order; immaterial for `nslyr = 1`), and
   !! `dhin = dhsn*rhos/rhoi` of ice is added to the top ice layer with the
   !! snow's own enthalpy (`dzi(1)*q1 + sum(dhs*qsn)` over the new thickness).
   !! Snow mass lost equals ice mass gained, `rhos*dhsn = rhoi*dhin`, and the
   !! interface lands on the waterline in one non-iterative step (`wk1 = 0`
   !! afterwards unless the snow ran out).  `snoice` accumulates `dhin` (the
   !! history diagnostic); `dsnow` accumulates the snow change (negative here).
   !!
   !! **Salt and fresh water (decision M-D3-3) - hand-off to the thermodynamic
   !! driver.**  `ice_freeboard` books nothing to the ocean, as upstream.  The
   !! driver books, AFTER `thickness_changes` returns, from
   !! `dhi = hin - worki` and `dhs = hsn - works - hsn_new` (`worki`, `works`
   !! the thicknesses before the thickness change, `hsn_new` the new snowfall
   !! thickness, which is EXCLUDED so snowfall is never an ocean freshwater
   !! sink): `fresh = fresh + evapn - (rhoi*dhi + rhos*dhs)/dt` and
   !! `fsalt = fsalt - rhoi*dhi*Sref*1e-3/dt`, both ACCUMULATED (`+=`), with
   !! `dhi` the TOTAL ice change (growth, melt AND the snow-ice `dhin`).  The
   !! prognostic-salinity salt-flux option is a separate branch and must fail
   !! loud.  Snow-ice formation then gives zero fresh water and salt
   !! `-rhoi*dhin*Sref*1e-3/dt` leaving the ocean (`rhoi*dhin = rhos*dhsn`).
   !! The driver also passes `Tbot` (set to `Tf`) as the `Tf` of
   !! `ice_update_state_vthermo`.  Ported as is; this module adds no flux.
   !!
   !! **Dropped.**  The snow-grain tracers (`massice`, `massliq`, `smice`,
   !! `smliq`, `rsnw`) and the mushy-physics salinity relayering are not ported
   !! (`snwgrain = .false.`, `ktherm = 1`), and the dead `ice age` lines.
   !!
   !! **Relayering.**  `ice_relayer` is the block of upstream's thickness
   !! change that follows the growth/melt: `hilyr = hin/nilyr` (0 if
   !! `hin <= 0`, and then `hin` is set to exactly 0), the same for snow, and
   !! `adjust_enthalpy` from the unequal old layers to the equal new ones.  It
   !! conserves `sum(h*q)` and `hin`; on return `dzi = hilyr` and `dzs = hslyr`
   !! (upstream leaves them stale; nothing reads them afterwards).  For
   !! `nslyr = 1` the snow enthalpy is not touched.
   !!
   !! The `nslyr > 1` snow relayering branch is unreachable while
   !! `ICE_NSLYR_MAX = 1` and has no test.
   !!
   !! **State update.**  `ice_update_state_vthermo` zeroes the category when
   !! `hin <= 0` exactly (area, volumes, all layer enthalpies, and `Tsf = Tf`)
   !! and otherwise sets `vicen = aicen*hin`, `vsnon = aicen*hsn`.  Quirk kept:
   !! when only the SNOW vanishes (`hsn = 0`, `hin > 0`), `vsnon` becomes 0 but
   !! the snow enthalpy and `Tsf` are left alone.  Downstream aggregates must
   !! therefore weight the snow enthalpy by `vsnon`; the vertical-profile
   !! initialisation of the next step resets `zqsn = -rhos*Lfresh` when the
   !! snow layer is below its threshold, which is how upstream discards the
   !! stale value.
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use, intrinsic :: iso_fortran_env, only: int32
   use rdb_constants, only: wp
   use rdb_ice_dims, only: ICE_NILYR_MAX, ICE_NSLYR_MAX
   use rdb_ice_errors, only: ICE_ERR_NONFINITE
   use rdb_ice_params, only: ice_constants_t
   use rdb_ice_thermo_shared, only: adjust_enthalpy
   implicit none
   private

   public :: ice_freeboard, ice_relayer, ice_update_state_vthermo

contains

   pure subroutine ice_freeboard(nilyr, nslyr, c, hin, hsn, zqin, zqsn, dzi, dzs, snoice, dsnow, ierr)
      !! Snow-ice formation when the snow load pushes the ice surface below
      !! sea level.  `hin, hsn` [m] total thicknesses, `zqin, zqsn` layer
      !! enthalpies [J/m3], `dzi, dzs` layer thicknesses [m] (all updated);
      !! `snoice` [m] and `dsnow` [m] are accumulated.  Bitwise no-op when the
      !! freeboard is not negative.
      !$acc routine seq
      integer, intent(in) :: nilyr, nslyr
      type(ice_constants_t), intent(in) :: c
      real(wp), intent(inout) :: hin, hsn
      real(wp), intent(inout) :: zqin(nilyr)
      real(wp), intent(in) :: zqsn(nslyr)
      real(wp), intent(inout) :: dzi(nilyr), dzs(nslyr)
      real(wp), intent(inout) :: snoice, dsnow
      integer(int32), intent(inout) :: ierr
      integer :: k
      real(wp) :: dhin, dhsn, hqs, wk1, dhs
      real(wp) :: hin_l, hsn_l, snoice_l, dsnow_l, zq1, dz1
      real(wp) :: dzs_l(ICE_NSLYR_MAX)

      if (nilyr < 1 .or. nilyr > ICE_NILYR_MAX .or. nslyr < 1 .or. nslyr > ICE_NSLYR_MAX) then
         ierr = ior(ierr, ICE_ERR_NONFINITE)
         return
      end if
      if (.not. (ieee_is_finite(hin) .and. ieee_is_finite(hsn) .and. ieee_is_finite(snoice) &
                 .and. ieee_is_finite(dsnow) .and. ieee_is_finite(zqin(1)) .and. ieee_is_finite(dzi(1)))) then
         ierr = ior(ierr, ICE_ERR_NONFINITE)
         return
      end if
      do k = 1, nslyr
         if (.not. (ieee_is_finite(zqsn(k)) .and. ieee_is_finite(dzs(k)))) then
            ierr = ior(ierr, ICE_ERR_NONFINITE)
            return
         end if
      end do

      dhin = 0.0_wp
      dhsn = 0.0_wp
      hqs = 0.0_wp
      wk1 = hsn - hin*(c%rhow - c%rhoi)/c%rhos
      if (wk1 > c%puny .and. hsn > c%puny) then
         dhsn = min(wk1*c%rhoi/c%rhow, hsn)
         dhin = dhsn*c%rhos/c%rhoi
      end if
      if (.not. (dhin > c%puny)) return

      hin_l = hin
      hsn_l = hsn
      snoice_l = snoice
      dsnow_l = dsnow
      do k = 1, nslyr
         dzs_l(k) = dzs(k)
      end do
      do k = nslyr, 1, -1
         dhs = min(dhsn, dzs_l(k))
         hsn_l = hsn_l - dhs
         dsnow_l = dsnow_l - dhs
         dzs_l(k) = dzs_l(k) - dhs
         dhsn = dhsn - dhs
         dhsn = max(dhsn, 0.0_wp)
         hqs = hqs + dhs*zqsn(k)
      end do
      wk1 = dzi(1) + dhin
      hin_l = hin_l + dhin
      zq1 = (dzi(1)*zqin(1) + hqs)/wk1
      dz1 = wk1
      snoice_l = snoice_l + dhin

      if (.not. (ieee_is_finite(zq1) .and. ieee_is_finite(hin_l) .and. ieee_is_finite(hsn_l))) then
         ierr = ior(ierr, ICE_ERR_NONFINITE)
         return
      end if
      hin = hin_l
      hsn = hsn_l
      snoice = snoice_l
      dsnow = dsnow_l
      zqin(1) = zq1
      dzi(1) = dz1
      do k = 1, nslyr
         dzs(k) = dzs_l(k)
      end do
   end subroutine ice_freeboard

   pure subroutine ice_relayer(nilyr, nslyr, c, dzi, dzs, hin, hsn, zqin, zqsn, hilyr, hslyr, ierr)
      !! Return the ice and snow columns to equal-thickness layers, conserving
      !! the thickness-weighted enthalpy.  `dzi, dzs` are the unequal layer
      !! thicknesses on entry and the new equal ones on exit; `hin, hsn` are
      !! set to exactly 0 when not positive; `hilyr, hslyr` are the new layer
      !! thicknesses.
      !$acc routine seq
      integer, intent(in) :: nilyr, nslyr
      type(ice_constants_t), intent(in) :: c
      real(wp), intent(inout) :: dzi(nilyr), dzs(nslyr)
      real(wp), intent(inout) :: hin, hsn
      real(wp), intent(inout) :: zqin(nilyr), zqsn(nslyr)
      real(wp), intent(out) :: hilyr, hslyr
      integer(int32), intent(inout) :: ierr
      integer :: k
      integer(int32) :: ierr_l
      real(wp) :: hin_l, hsn_l, hil, hsl
      real(wp) :: zi1(ICE_NILYR_MAX + 1), zi2(ICE_NILYR_MAX + 1)
      real(wp) :: zs1(ICE_NSLYR_MAX + 1), zs2(ICE_NSLYR_MAX + 1)
      real(wp) :: zqi_l(ICE_NILYR_MAX), zqs_l(ICE_NSLYR_MAX)

      hilyr = 0.0_wp
      hslyr = 0.0_wp
      if (nilyr < 1 .or. nilyr > ICE_NILYR_MAX .or. nslyr < 1 .or. nslyr > ICE_NSLYR_MAX) then
         ierr = ior(ierr, ICE_ERR_NONFINITE)
         return
      end if
      if (.not. (ieee_is_finite(hin) .and. ieee_is_finite(hsn))) then
         ierr = ior(ierr, ICE_ERR_NONFINITE)
         return
      end if
      do k = 1, nilyr
         if (.not. (ieee_is_finite(zqin(k)) .and. ieee_is_finite(dzi(k)))) then
            ierr = ior(ierr, ICE_ERR_NONFINITE)
            return
         end if
      end do
      do k = 1, nslyr
         if (.not. (ieee_is_finite(zqsn(k)) .and. ieee_is_finite(dzs(k)))) then
            ierr = ior(ierr, ICE_ERR_NONFINITE)
            return
         end if
      end do

      hin_l = hin
      hsn_l = hsn
      if (hin_l > 0.0_wp) then
         hil = hin_l/real(nilyr, wp)
      else
         hin_l = 0.0_wp
         hil = 0.0_wp
      end if
      if (hsn_l > 0.0_wp) then
         hsl = hsn_l/real(nslyr, wp)
      else
         hsn_l = 0.0_wp
         hsl = 0.0_wp
      end if

      zi1(1) = 0.0_wp
      zi1(1 + nilyr) = hin_l
      zi2(1) = 0.0_wp
      zi2(1 + nilyr) = hin_l
      do k = 1, nilyr - 1
         zi1(k + 1) = zi1(k) + dzi(k)
         zi2(k + 1) = zi2(k) + hil
      end do
      do k = 1, nilyr
         zqi_l(k) = zqin(k)
      end do
      ierr_l = 0_int32
      call adjust_enthalpy(nilyr, zi1, zi2, hil, hin_l, zqi_l, c, ierr_l)

      do k = 1, nslyr
         zqs_l(k) = zqsn(k)
      end do
      if (nslyr > 1) then
         zs1(1) = 0.0_wp
         zs1(1 + nslyr) = hsn_l
         zs2(1) = 0.0_wp
         zs2(1 + nslyr) = hsn_l
         do k = 1, nslyr - 1
            zs1(k + 1) = zs1(k) + dzs(k)
            zs2(k + 1) = zs2(k) + hsl
         end do
         call adjust_enthalpy(nslyr, zs1, zs2, hsl, hsn_l, zqs_l, c, ierr_l)
      end if

      if (ierr_l /= 0_int32) then
         ierr = ior(ierr, ierr_l)
         return
      end if
      hin = hin_l
      hsn = hsn_l
      hilyr = hil
      hslyr = hsl
      do k = 1, nilyr
         zqin(k) = zqi_l(k)
         dzi(k) = hil
      end do
      do k = 1, nslyr
         zqsn(k) = zqs_l(k)
         dzs(k) = hsl
      end do
   end subroutine ice_relayer

   pure subroutine ice_update_state_vthermo(nilyr, nslyr, Tf, hin, hsn, Tsf, zqin, zqsn, aicen, vicen, vsnon, ierr)
      !! Category volumes from the column thicknesses.  `Tf` freezing
      !! temperature [degC] (assigned to `Tsf` when the ice is gone), `hin,
      !! hsn` [m], `aicen` area fraction (unchanged except when the ice melted
      !! out, then 0), `vicen, vsnon` [m] outputs.
      !$acc routine seq
      integer, intent(in) :: nilyr, nslyr
      real(wp), intent(in) :: Tf, hin, hsn
      real(wp), intent(inout) :: Tsf
      real(wp), intent(inout) :: zqin(nilyr), zqsn(nslyr)
      real(wp), intent(inout) :: aicen, vicen, vsnon
      integer(int32), intent(inout) :: ierr
      integer :: k

      if (.not. (ieee_is_finite(Tf) .and. ieee_is_finite(hin) .and. ieee_is_finite(hsn) &
                 .and. ieee_is_finite(Tsf) .and. ieee_is_finite(aicen))) then
         ierr = ior(ierr, ICE_ERR_NONFINITE)
         return
      end if
      do k = 1, nilyr
         if (.not. ieee_is_finite(zqin(k))) then
            ierr = ior(ierr, ICE_ERR_NONFINITE)
            return
         end if
      end do
      do k = 1, nslyr
         if (.not. ieee_is_finite(zqsn(k))) then
            ierr = ior(ierr, ICE_ERR_NONFINITE)
            return
         end if
      end do

      if (hin <= 0.0_wp) then
         aicen = 0.0_wp
         vicen = 0.0_wp
         vsnon = 0.0_wp
         Tsf = Tf
         do k = 1, nilyr
            zqin(k) = 0.0_wp
         end do
         do k = 1, nslyr
            zqsn(k) = 0.0_wp
         end do
      else
         vicen = aicen*hin
         vsnon = aicen*hsn
      end if
   end subroutine ice_update_state_vthermo

end module rdb_ice_thickness
