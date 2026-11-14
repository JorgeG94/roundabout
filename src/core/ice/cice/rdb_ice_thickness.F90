!! Column thickness bookkeeping of the Bitz-Lipscomb thermodynamics: the
!! growth and melt of the column from the flux imbalances
!! (`ice_thickness_changes`), its energy check (`ice_conservation_check_vthermo`),
!! snow-ice formation at the waterline (`ice_freeboard`), the conservative return
!! to equal layer thicknesses (`ice_relayer`) and the conversion of the column
!! thicknesses back to the category volumes (`ice_update_state_vthermo`).
module rdb_ice_thickness
   ! SPDX-License-Identifier: BSD-3-Clause
   !
   ! Ported from CICE/Icepack (CICE Consortium, https://github.com/CICE-Consortium).
   ! Upstream: Icepack 9a1f3e4 (Icepack 1.5.3+14), module icepack_therm_vertical
   ! (thickness_changes, conservation_check_vthermo, freeboard,
   ! update_state_vthermo).
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
   !! **Thickness changes (`ice_thickness_changes`, ktherm = 1).**  After the
   !! temperature solve the column holds equal layers (`hilyr`, `hslyr`) with
   !! enthalpies `zqin`, `zqsn` [J/m3] and the surface and conductive fluxes
   !! (W/m2, positive downward into the ice) `fsurfn`, `fcondtopn`, `fcondbotn`,
   !! the basal ocean flux `fbot` (<= 0) and the latent flux `flatn`.  The
   !! energies available over the step are `esub = max(-flatn dt, 0)`
   !! (sublimation), `econ = min(-flatn dt, 0)` (condensation),
   !! `etop_mlt = max((fsurfn - fcondtopn) dt, 0)`, `ebot_mlt = max((fcondbotn -
   !! fbot) dt, 0)` and `ebot_gro = min((fcondbotn - fbot) dt, 0)`.  They are
   !! spent in this order, each consumption feeding the next (the order is part
   !! of the answer): condensation onto the snow (or, with no snow, onto the top
   !! ice layer); congelation at the base, the new ice entering at the BL99
   !! enthalpy `qbot = -rhoi (cp_ice (Tmlts - Tbot) + Lfresh (1 - Tmlts/Tbot) -
   !! cp_ocn Tmlts)`, `Tmlts = -S_base depressT`, capped at `-rhoi Lfresh/2`, and
   !! mixed thickness-weighted into the bottom layer only if that layer ends
   !! thicker than `puny` (otherwise its enthalpy is kept, as upstream); per
   !! snow layer, sublimation then top melt; per ice layer top-down,
   !! sublimation then top melt; per ice layer bottom-up, bottom melt; the snow
   !! melted from below once the ice is gone; the energy left over in
   !! `esub + etop_mlt + ebot_mlt` returns to the ocean in `fhocnn = fbot +
   !! (...)/dt`; the new snowfall (`hsn_new = fsnow dt/rhos` at `-rhos
   !! Lfresh`); then `ice_freeboard` and `ice_relayer`.  A layer melts by
   !! `max(-dz, E/q)` with `q < 0` the layer enthalpy, so a layer never melts
   !! below zero thickness and the energy it cannot take stays in the budget.
   !! `efinal = -evap Lvap + sum(hslyr zqsn) + sum(hilyr zqin)` is formed with
   !! `evap` still in kg/m2 and the sublimation/condensation fluxes
   !! `evapn, evapsn, evapin` returned per second afterwards.  The history
   !! accumulators `meltt` (top ice melt), `melts` (snow melt), `meltb` (basal
   !! melt), `congel`, `snoice`, `dsnow` [m] are ACCUMULATED (`+=`), so the
   !! caller zeroes them once per step.  `l_brine = .false.` keeps upstream's
   !! fresh-ice branch (snow or ice layers warmer than 0 degC are melted back to
   !! 0 degC and the base ice is formed at `-rhoi (Lfresh - cp_ice Tbot)`); the
   !! melt cascade then still reads the pre-clamp layer enthalpy, as upstream.
   !!
   !! **Errors.**  A non-finite input, `dt <= 0`, a non-negative snow
   !! enthalpy on a layer that carries snow (`hslyr > hs_min/nslyr`, the
   !! threshold of the vertical profile; the snow melt divides by it) or a layer
   !! count outside the bounds sets `ICE_ERR_NONFINITE` before any output is
   !! written (`dt <= 0` and the layer count have no bit of their own and are
   !! reported under it).  A snow layer at or below that threshold whose
   !! enthalpy is not negative (a raw snow-free category, `qsno = 0`) has it
   !! WRITTEN BACK as `-rhos*Lfresh` (the melting-point value the profile
   !! writes; idempotent after it) and is a valid input; a negative value is
   !! kept as is.  Under `l_brine`, `Tbot >
   !! -puny` (a base at or above the freezing point of fresh water, where the
   !! growth enthalpy divides by `Tbot`; slightly stricter than upstream in
   !! (-puny, 0)) sets `ICE_ERR_TBOT`; a layer
   !! thickness below `-puny` on entry or after the melt sets
   !! `ICE_ERR_NEG_THICK`; a non-finite result sets `ICE_ERR_NONFINITE`.  The
   !! outputs are then unspecified and the driver does not commit the column.
   !! Nothing is buffered: outputs are written in place.  The latent-heat
   !! flux correction, the mushy-layer physics, the snow-grain tracers, the
   !! melt-onset and freeze-onset days (history diagnostics, owned by a later
   !! PR) and the excess-energy arguments of calc_Tsfc = false
   !! are not ported (all inert for ktherm = 1, calc_Tsfc = true).
   !!
   !! **Conservation check.**  `ice_conservation_check_vthermo` evaluates
   !! `ferr = |efinal - einit - einp|/dt` with the energy input over the step
   !! `einp = (fsurfn - flatn + fswint - fhocnn - fsnow Lfresh) dt` and sets
   !! `ICE_ERR_THERMO_CONSERV` when `ferr > 1.1 ferrmax` (`ICE_FERRMAX` = 1e-3
   !! W/m2, upstream).  It never aborts: the driver reduces the bit.
   !!
   !! **Cell entry.**  `ice_thickness_changes_cell(i, j, n, ...)` gathers one
   !! `(i, j, n)` column from the `(nx, ny, ncat[, nlyr])` state arrays, calls
   !! the column routine and scatters once; the indices are passed BY VALUE.
   !! Every output is per category `(nx, ny, ncat)` (no cell accumulator), so
   !! the driver runs it from `do concurrent (n, j, i)`; the merge into the
   !! cell fluxes happens in the next per-cell kernel.
   !!
   !! **Layer convention.**  Ice and snow layers are indexed TOP-DOWN (k = 1 is
   !! the top), the Icepack order; `nilyr <= ICE_NILYR_MAX`, `nslyr <=
   !! ICE_NSLYR_MAX`.  All arguments are explicit-shape; locals are sized by the
   !! `ICE_*_MAX` bounds (no allocation).  Every procedure is `pure` with
   !! `!$acc routine seq`, takes the constants in the `intent(in)` bundle and
   !! ORs its error bits into `ierr`.  A non-finite input sets
   !! `ICE_ERR_NONFINITE` and leaves every state argument untouched (results
   !! are formed in locals and committed only when finite); the pure outputs
   !! `hilyr`/`hslyr` of `ice_relayer` are then set to 0.  (The thickness-change
   !! routine above, written later, does not buffer.)
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
   use rdb_ice_errors, only: ICE_ERR_NONFINITE, ICE_ERR_THERMO_CONSERV, ICE_ERR_TBOT, ICE_ERR_NEG_THICK
   use rdb_ice_params, only: ice_constants_t
   use rdb_ice_thermo_shared, only: adjust_enthalpy
   implicit none
   private

   public :: ICE_FERRMAX
   public :: ice_thickness_changes, ice_thickness_changes_cell
   public :: ice_conservation_check_vthermo, ice_conservation_check_vthermo_cell
   public :: ice_freeboard, ice_relayer, ice_update_state_vthermo

   real(wp), parameter :: ICE_FERRMAX = 1.0e-3_wp
      !! Tolerated energy-conservation error [W/m2] of the thermodynamic step (upstream value).

contains

   pure subroutine ice_thickness_changes(nilyr, nslyr, c, l_brine, dt, Tbot, fbot, flatn, fsurfn, fcondtopn, fcondbotn, &
                                         fsnow, zSin_bot, hin, hilyr, hsn, hslyr, zqin, zqsn, &
                                         meltt, melts, meltb, congel, snoice, dsnow, &
                                         hsn_new, fhocnn, evapn, evapsn, evapin, efinal, ierr)
      !! Growth, melt, sublimation/condensation and snowfall of one column over
      !! `dt` [s] from the flux imbalances, then freeboard and relayering.
      !! In: `Tbot` [degC] base temperature (= the freezing point), `fbot` [W/m2]
      !! ocean-to-ice flux (<= 0), `flatn, fsurfn, fcondtopn, fcondbotn` [W/m2]
      !! (positive downward), `fsnow` [kg/m2/s] snowfall, `zSin_bot` [ppt]
      !! salinity of the bottom layer.  State, updated in place: `hin, hsn` [m]
      !! (total thicknesses), `hilyr, hslyr` [m] (layer thicknesses; equal
      !! layers on entry and exit), `zqin, zqsn` [J/m3].  Accumulated [m]:
      !! `meltt, melts, meltb, congel, snoice, dsnow`.  Out: `hsn_new` [m]
      !! (new snow), `fhocnn` [W/m2] (+ into the ocean), `evapn, evapsn, evapin`
      !! [kg/m2/s] (total, over snow, over ice; < 0 sublimation), `efinal`
      !! [J/m2] (the column energy after the step, for the conservation check).
      !$acc routine seq
      integer, intent(in) :: nilyr, nslyr
      type(ice_constants_t), intent(in) :: c
      logical, intent(in) :: l_brine
      real(wp), intent(in) :: dt, Tbot, fbot, flatn, fsurfn, fcondtopn, fcondbotn, fsnow, zSin_bot
      real(wp), intent(inout) :: hin, hilyr, hsn, hslyr
      real(wp), intent(inout) :: zqin(nilyr), zqsn(nslyr)
      real(wp), intent(inout) :: meltt, melts, meltb, congel, snoice, dsnow
      real(wp), intent(out) :: hsn_new, fhocnn, evapn, evapsn, evapin, efinal
      integer(int32), intent(inout) :: ierr
      integer :: k
      real(wp) :: esub, econ, etop_mlt, ebot_mlt, ebot_gro
      real(wp) :: dhi, dhs, wk1, qbot, qbotmax, qsub, hqtot, hstot, zqsnew, Tmlts, Ts, Ti
      real(wp) :: dzi(ICE_NILYR_MAX), dzs(ICE_NSLYR_MAX), qm(ICE_NILYR_MAX)

      hsn_new = 0.0_wp
      fhocnn = 0.0_wp
      evapn = 0.0_wp
      evapsn = 0.0_wp
      evapin = 0.0_wp
      efinal = 0.0_wp
      if (nilyr < 1 .or. nilyr > ICE_NILYR_MAX .or. nslyr < 1 .or. nslyr > ICE_NSLYR_MAX) then
         ierr = ior(ierr, ICE_ERR_NONFINITE)
         return
      end if
      if (.not. (ieee_is_finite(dt) .and. ieee_is_finite(Tbot) .and. ieee_is_finite(fbot) .and. ieee_is_finite(flatn) &
                 .and. ieee_is_finite(fsurfn) .and. ieee_is_finite(fcondtopn) .and. ieee_is_finite(fcondbotn) &
                 .and. ieee_is_finite(fsnow) .and. ieee_is_finite(zSin_bot) .and. ieee_is_finite(hin) &
                 .and. ieee_is_finite(hsn) .and. ieee_is_finite(hilyr) .and. ieee_is_finite(hslyr) &
                 .and. ieee_is_finite(meltt) .and. ieee_is_finite(melts) .and. ieee_is_finite(meltb) &
                 .and. ieee_is_finite(congel) .and. ieee_is_finite(snoice) .and. ieee_is_finite(dsnow))) then
         ierr = ior(ierr, ICE_ERR_NONFINITE)
         return
      end if
      if (.not. (dt > 0.0_wp)) then
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
         ! The snow melt divides by the snow enthalpy: for a layer that carries snow a non-negative value
         ! is not a snow state.  A layer at or below the profile's snow threshold is reset below.
         if (.not. ieee_is_finite(zqsn(k))) then
            ierr = ior(ierr, ICE_ERR_NONFINITE)
            return
         end if
         if (hslyr > c%hs_min/real(nslyr, wp) .and. zqsn(k) >= 0.0_wp) then
            ierr = ior(ierr, ICE_ERR_NONFINITE)
            return
         end if
      end do
      ! (dt <= 0 and a layer count out of range have no bit of their own: reported as ICE_ERR_NONFINITE.)
      if (l_brine .and. Tbot > -c%puny) then
         ierr = ior(ierr, ICE_ERR_TBOT)
         return
      end if
      if (hilyr < -c%puny .or. hslyr < -c%puny) then
         ierr = ior(ierr, ICE_ERR_NEG_THICK)
         return
      end if

      ! A snow layer too thin to carry a temperature takes the melting-point enthalpy (what the vertical
      ! profile writes; idempotent after it), so a raw snow-free category (qsno = 0) is a valid input.
      if (hslyr <= c%hs_min/real(nslyr, wp)) then
         do k = 1, nslyr
            if (zqsn(k) >= 0.0_wp) zqsn(k) = -c%rhos*c%Lfresh
         end do
      end if

      dhi = 0.0_wp
      dhs = 0.0_wp
      do k = 1, nilyr
         dzi(k) = hilyr
         qm(k) = zqin(k)
      end do
      do k = 1, nslyr
         dzs(k) = hslyr
      end do

      ! Fresh-ice option: layers warmer than 0 degC are melted back to 0 degC (the melt
      ! cascade below keeps the pre-clamp enthalpy `qm`, as upstream).
      if (.not. l_brine) then
         do k = 1, nslyr
            Ts = (c%Lfresh + zqsn(k)/c%rhos)/c%cp_ice
            if (Ts > 0.0_wp) then
               dhs = c%cp_ice*Ts*dzs(k)/c%Lfresh
               dzs(k) = dzs(k) - dhs
               melts = melts + dhs
               zqsn(k) = -c%rhos*c%Lfresh
            end if
         end do
         do k = 1, nilyr
            Ti = (c%Lfresh + zqin(k)/c%rhoi)/c%cp_ice
            if (Ti > 0.0_wp) then
               dhi = c%cp_ice*Ti*dzi(k)/c%Lfresh
               dzi(k) = dzi(k) - dhi
               zqin(k) = -c%rhoi*c%Lfresh
            end if
         end do
      end if

      ! Energy available at the surfaces [J/m2].
      wk1 = -flatn*dt
      esub = max(wk1, 0.0_wp)
      econ = min(wk1, 0.0_wp)
      wk1 = (fsurfn - fcondtopn)*dt
      etop_mlt = max(wk1, 0.0_wp)
      wk1 = (fcondbotn - fbot)*dt
      ebot_mlt = max(wk1, 0.0_wp)
      ebot_gro = min(wk1, 0.0_wp)

      ! Condensation: onto the snow, or onto the ice when there is none.
      if (hsn > c%puny) then
         dhs = econ/(zqsn(1) - c%rhos*c%Lvap)
         dzs(1) = dzs(1) + dhs
         evapn = evapn + dhs*c%rhos
         evapsn = evapsn + dhs*c%rhos
      else
         dhi = econ/(qm(1) - c%rhoi*c%Lvap)
         dzi(1) = dzi(1) + dhi
         evapn = evapn + dhi*c%rhoi
         evapin = evapin + dhi*c%rhoi
      end if

      ! Congelation at the base.
      if (l_brine) then
         Tmlts = -zSin_bot*c%depressT
         qbotmax = -0.5_wp*c%rhoi*c%Lfresh
         qbot = -c%rhoi*(c%cp_ice*(Tmlts - Tbot) + c%Lfresh*(1.0_wp - Tmlts/Tbot) - c%cp_ocn*Tmlts)
         qbot = min(qbot, qbotmax)
      else
         qbot = -c%rhoi*(-c%cp_ice*Tbot + c%Lfresh)
      end if
      dhi = ebot_gro/qbot
      hqtot = dzi(nilyr)*zqin(nilyr) + dhi*qbot
      dzi(nilyr) = dzi(nilyr) + dhi
      if (dzi(nilyr) > c%puny) then
         zqin(nilyr) = hqtot/dzi(nilyr)
         qm(nilyr) = zqin(nilyr)
      end if
      congel = congel + dhi

      ! Snow: sublimation, then top melt.
      do k = 1, nslyr
         qsub = zqsn(k) - c%rhos*c%Lvap
         dhs = max(-dzs(k), esub/qsub)
         dzs(k) = dzs(k) + dhs
         esub = esub - dhs*qsub
         esub = max(esub, 0.0_wp)
         evapn = evapn + dhs*c%rhos
         evapsn = evapsn + dhs*c%rhos
         dhs = max(-dzs(k), etop_mlt/zqsn(k))
         dzs(k) = dzs(k) + dhs
         etop_mlt = etop_mlt - dhs*zqsn(k)
         etop_mlt = max(etop_mlt, 0.0_wp)
         melts = melts - dhs
      end do

      ! Ice, top-down: sublimation, then top melt.
      do k = 1, nilyr
         qsub = qm(k) - c%rhoi*c%Lvap
         dhi = max(-dzi(k), esub/qsub)
         dzi(k) = dzi(k) + dhi
         esub = esub - dhi*qsub
         esub = max(esub, 0.0_wp)
         evapn = evapn + dhi*c%rhoi
         evapin = evapin + dhi*c%rhoi
         if (qm(k) < 0.0_wp) then
            dhi = max(-dzi(k), etop_mlt/qm(k))
         else
            qm(k) = 0.0_wp
            dhi = -dzi(k)
         end if
         dzi(k) = dzi(k) + dhi
         etop_mlt = max(etop_mlt - dhi*qm(k), 0.0_wp)
         meltt = meltt - dhi
      end do

      ! Ice, bottom-up: basal melt.
      do k = nilyr, 1, -1
         if (qm(k) < 0.0_wp) then
            dhi = max(-dzi(k), ebot_mlt/qm(k))
         else
            qm(k) = 0.0_wp
            dhi = -dzi(k)
         end if
         dzi(k) = dzi(k) + dhi
         ebot_mlt = max(ebot_mlt - dhi*qm(k), 0.0_wp)
         meltb = meltb - dhi
      end do

      ! Snow melted from below (reached only once the ice is gone).
      do k = nslyr, 1, -1
         dhs = max(-dzs(k), ebot_mlt/zqsn(k))
         dzs(k) = dzs(k) + dhs
         ebot_mlt = ebot_mlt - dhs*zqsn(k)
         ebot_mlt = max(ebot_mlt, 0.0_wp)
         melts = melts - dhs
      end do

      do k = 1, nilyr
         if (dzi(k) < -c%puny) then
            ierr = ior(ierr, ICE_ERR_NEG_THICK)
            return
         end if
      end do
      do k = 1, nslyr
         if (dzs(k) < -c%puny) then
            ierr = ior(ierr, ICE_ERR_NEG_THICK)
            return
         end if
      end do

      ! Energy not used by the melt returns to the ocean.
      fhocnn = fbot + (esub + etop_mlt + ebot_mlt)/dt

      ! New snowfall at 0 degC.
      if (fsnow > 0.0_wp) then
         hsn_new = fsnow/c%rhos*dt
         zqsnew = -c%rhos*c%Lfresh
         hstot = dzs(1) + hsn_new
         if (hstot > 0.0_wp) then
            zqsn(1) = (dzs(1)*zqsn(1) + hsn_new*zqsnew)/hstot
            zqsn(1) = min(zqsn(1), zqsnew)
            dzs(1) = hstot
         end if
      end if

      hin = 0.0_wp
      hsn = 0.0_wp
      do k = 1, nilyr
         hin = hin + dzi(k)
      end do
      do k = 1, nslyr
         hsn = hsn + dzs(k)
         dsnow = dsnow + dzs(k) - hslyr
      end do

      call ice_freeboard(nilyr, nslyr, c, hin, hsn, zqin, zqsn, dzi, dzs, snoice, dsnow, ierr)
      call ice_relayer(nilyr, nslyr, c, dzi, dzs, hin, hsn, zqin, zqsn, hilyr, hslyr, ierr)

      ! Column energy after the step (evap still a mass), then the fluxes per second.
      efinal = -evapn*c%Lvap
      evapn = evapn/dt
      evapsn = evapsn/dt
      evapin = evapin/dt
      do k = 1, nslyr
         efinal = efinal + hslyr*zqsn(k)
      end do
      do k = 1, nilyr
         efinal = efinal + hilyr*zqin(k)
      end do

      if (.not. (ieee_is_finite(hin) .and. ieee_is_finite(hsn) .and. ieee_is_finite(hilyr) .and. ieee_is_finite(hslyr) &
                 .and. ieee_is_finite(fhocnn) .and. ieee_is_finite(efinal) .and. ieee_is_finite(evapn) &
                 .and. ieee_is_finite(evapsn) .and. ieee_is_finite(evapin) &
                 .and. ieee_is_finite(dsnow) .and. ieee_is_finite(meltt) .and. ieee_is_finite(meltb) &
                 .and. ieee_is_finite(melts) .and. ieee_is_finite(congel) .and. ieee_is_finite(snoice))) then
         ierr = ior(ierr, ICE_ERR_NONFINITE)
      end if
      do k = 1, nilyr
         if (.not. ieee_is_finite(zqin(k))) ierr = ior(ierr, ICE_ERR_NONFINITE)
      end do
      do k = 1, nslyr
         if (.not. ieee_is_finite(zqsn(k))) ierr = ior(ierr, ICE_ERR_NONFINITE)
      end do
   end subroutine ice_thickness_changes

   pure subroutine ice_conservation_check_vthermo(c, dt, fsurfn, flatn, fhocnn, fswint, fsnow, einit, efinal, ferrmax, ferr, ierr)
      !! Energy conservation of the thermodynamic step of one column: the
      !! energy change `efinal - einit` [J/m2] must equal the net input `einp =
      !! (fsurfn - flatn + fswint - fhocnn - fsnow Lfresh) dt` (`fswint` the
      !! shortwave absorbed in the ice interior, `fhocnn` the heat returned to
      !! the ocean).  Sets `ICE_ERR_THERMO_CONSERV` when `ferr = |efinal - einit -
      !! einp|/dt` [W/m2] exceeds `1.1 ferrmax`; a non-finite input sets
      !! `ICE_ERR_NONFINITE` and `ferr = 0`.
      !$acc routine seq
      type(ice_constants_t), intent(in) :: c
      real(wp), intent(in) :: dt, fsurfn, flatn, fhocnn, fswint, fsnow, einit, efinal, ferrmax
      real(wp), intent(out) :: ferr
      integer(int32), intent(inout) :: ierr
      real(wp) :: einp

      ferr = 0.0_wp
      if (.not. (ieee_is_finite(dt) .and. ieee_is_finite(fsurfn) .and. ieee_is_finite(flatn) .and. ieee_is_finite(fhocnn) &
                 .and. ieee_is_finite(fswint) .and. ieee_is_finite(fsnow) .and. ieee_is_finite(einit) &
                 .and. ieee_is_finite(efinal) .and. dt > 0.0_wp)) then
         ierr = ior(ierr, ICE_ERR_NONFINITE)
         return
      end if
      einp = (fsurfn - flatn + fswint - fhocnn - fsnow*c%Lfresh)*dt
      ferr = abs(efinal - einit - einp)/dt
      if (.not. ieee_is_finite(ferr)) then
         ierr = ior(ierr, ICE_ERR_NONFINITE)
      else if (ferr > 1.1_wp*ferrmax) then
         ierr = ior(ierr, ICE_ERR_THERMO_CONSERV)
      end if
   end subroutine ice_conservation_check_vthermo

   pure subroutine ice_thickness_changes_cell(i, j, n, nx, ny, ncat, nilyr, nslyr, c, l_brine, dt, &
                                              Tbot, fbot, fsnow, flatn, fsurfn, fcondtopn, fcondbotn, sicen, &
                                              hin, hsn, hilyr, hslyr, qicen, qsnon, &
                                              meltt, melts, meltb, congel, snoice, dsnow, &
                                              hsn_new, fhocnn, evapn, evapsn, evapin, efinal, ierr)
      !! One `(i, j, n)` column of `ice_thickness_changes` on the state arrays:
      !! the column is gathered once, the column routine run, and the results
      !! scattered once.  `i, j, n` are passed by value (a compiler fault
      !! passes a bare non-innermost `do concurrent` index by reference as its
      !! lower bound).  `Tbot, fbot, fsnow` are per cell `(nx, ny)`; `flatn,
      !! fsurfn, fcondtopn, fcondbotn, hin, hsn, hilyr, hslyr` and every
      !! output are per category `(nx, ny, ncat)`; `sicen`, `qicen` are `(nx, ny,
      !! ncat, nilyr)` (only the bottom layer of `sicen` is read), `qsnon` is
      !! `(nx, ny, ncat, nslyr)`.  Updated in place: `hin, hsn, hilyr, hslyr,
      !! qicen, qsnon` and the accumulators `meltt, melts, meltb, congel, snoice,
      !! dsnow`; written: `hsn_new, fhocnn, evapn, evapsn, evapin, efinal`; the
      !! error word of the category is OR-ed into `ierr(i, j, n)`.  When a bit is
      !! set the scattered state is UNSPECIFIED (mid-way returns happen after some
      !! writes); that is acceptable only because every bit this routine sets
      !! (`ICE_ERR_NONFINITE`, `ICE_ERR_TBOT`, `ICE_ERR_NEG_THICK`) is FATAL in
      !! `rdb_ice_errors` and the driver does not commit the column.  None of
      !! them may become a counted (non-fatal) bit without revisiting this routine.
      !$acc routine seq
      integer, value :: i, j, n
      integer, intent(in) :: nx, ny, ncat, nilyr, nslyr
      type(ice_constants_t), intent(in) :: c
      logical, intent(in) :: l_brine
      real(wp), intent(in) :: dt
      real(wp), intent(in) :: Tbot(nx, ny), fbot(nx, ny), fsnow(nx, ny)
      real(wp), intent(in) :: flatn(nx, ny, ncat), fsurfn(nx, ny, ncat), fcondtopn(nx, ny, ncat), fcondbotn(nx, ny, ncat)
      real(wp), intent(in) :: sicen(nx, ny, ncat, nilyr)
      real(wp), intent(inout) :: hin(nx, ny, ncat), hsn(nx, ny, ncat), hilyr(nx, ny, ncat), hslyr(nx, ny, ncat)
      real(wp), intent(inout) :: qicen(nx, ny, ncat, nilyr), qsnon(nx, ny, ncat, nslyr)
      real(wp), intent(inout) :: meltt(nx, ny, ncat), melts(nx, ny, ncat), meltb(nx, ny, ncat)
      real(wp), intent(inout) :: congel(nx, ny, ncat), snoice(nx, ny, ncat), dsnow(nx, ny, ncat)
      real(wp), intent(out) :: hsn_new(nx, ny, ncat), fhocnn(nx, ny, ncat), evapn(nx, ny, ncat)
      real(wp), intent(out) :: evapsn(nx, ny, ncat), evapin(nx, ny, ncat), efinal(nx, ny, ncat)
      integer(int32), intent(inout) :: ierr(nx, ny, ncat)
      integer :: k
      integer(int32) :: ie
      real(wp) :: zqi(ICE_NILYR_MAX), zqs(ICE_NSLYR_MAX)
      real(wp) :: h_in, h_sn, h_il, h_sl, mt, ms, mb, cg, si, ds, hn, fh, ev, evs, evi, ef

      ie = 0_int32
      do k = 1, nilyr
         zqi(k) = qicen(i, j, n, k)
      end do
      do k = 1, nslyr
         zqs(k) = qsnon(i, j, n, k)
      end do
      h_in = hin(i, j, n)
      h_sn = hsn(i, j, n)
      h_il = hilyr(i, j, n)
      h_sl = hslyr(i, j, n)
      mt = meltt(i, j, n)
      ms = melts(i, j, n)
      mb = meltb(i, j, n)
      cg = congel(i, j, n)
      si = snoice(i, j, n)
      ds = dsnow(i, j, n)
      call ice_thickness_changes(nilyr, nslyr, c, l_brine, dt, Tbot(i, j), fbot(i, j), flatn(i, j, n), fsurfn(i, j, n), &
                                 fcondtopn(i, j, n), fcondbotn(i, j, n), fsnow(i, j), sicen(i, j, n, nilyr), &
                                 h_in, h_il, h_sn, h_sl, zqi, zqs, mt, ms, mb, cg, si, ds, &
                                 hn, fh, ev, evs, evi, ef, ie)
      hin(i, j, n) = h_in
      hsn(i, j, n) = h_sn
      hilyr(i, j, n) = h_il
      hslyr(i, j, n) = h_sl
      do k = 1, nilyr
         qicen(i, j, n, k) = zqi(k)
      end do
      do k = 1, nslyr
         qsnon(i, j, n, k) = zqs(k)
      end do
      meltt(i, j, n) = mt
      melts(i, j, n) = ms
      meltb(i, j, n) = mb
      congel(i, j, n) = cg
      snoice(i, j, n) = si
      dsnow(i, j, n) = ds
      hsn_new(i, j, n) = hn
      fhocnn(i, j, n) = fh
      evapn(i, j, n) = ev
      evapsn(i, j, n) = evs
      evapin(i, j, n) = evi
      efinal(i, j, n) = ef
      ierr(i, j, n) = ior(ierr(i, j, n), ie)
   end subroutine ice_thickness_changes_cell

   pure subroutine ice_conservation_check_vthermo_cell(i, j, n, nx, ny, ncat, c, dt, fsurfn, flatn, fhocnn, fswint, fsnow, &
                                                       einit, efinal, ferrmax, ferr, ierr)
      !! The conservation check of one `(i, j, n)` column on the state arrays
      !! (`fsnow` per cell; the rest per category); `ferr(i, j, n)` [W/m2] is
      !! written and the error word OR-ed into `ierr(i, j, n)`.
      !$acc routine seq
      integer, value :: i, j, n
      integer, intent(in) :: nx, ny, ncat
      type(ice_constants_t), intent(in) :: c
      real(wp), intent(in) :: dt, ferrmax
      real(wp), intent(in) :: fsurfn(nx, ny, ncat), flatn(nx, ny, ncat), fhocnn(nx, ny, ncat), fswint(nx, ny, ncat)
      real(wp), intent(in) :: fsnow(nx, ny), einit(nx, ny, ncat), efinal(nx, ny, ncat)
      real(wp), intent(out) :: ferr(nx, ny, ncat)
      integer(int32), intent(inout) :: ierr(nx, ny, ncat)
      integer(int32) :: ie
      real(wp) :: fe

      ie = 0_int32
      call ice_conservation_check_vthermo(c, dt, fsurfn(i, j, n), flatn(i, j, n), fhocnn(i, j, n), fswint(i, j, n), &
                                          fsnow(i, j), einit(i, j, n), efinal(i, j, n), ferrmax, fe, ie)
      ferr(i, j, n) = fe
      ierr(i, j, n) = ior(ierr(i, j, n), ie)
   end subroutine ice_conservation_check_vthermo_cell

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
