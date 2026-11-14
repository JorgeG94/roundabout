!! Post-thermodynamics / post-ridging clean-up of the ice-thickness distribution: rebin, zap of residual ice and snow.
module rdb_ice_itd_cleanup
   ! SPDX-License-Identifier: BSD-3-Clause
   !
   ! Ported from CICE/Icepack (CICE Consortium, https://github.com/CICE-Consortium).
   ! Upstream: Icepack 9a1f3e4 (Icepack 1.5.3+14), module icepack_itd
   ! (cleanup_itd, zap_small_areas, zap_snow, zap_snow_temperature).
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
   ! Original authors (from the upstream headers): William H. Lipscomb and
   ! Elizabeth C. Hunke (LANL), C. M. Bitz (UW).
   ! Physics: Lipscomb (2001, JGR 106 13989); Bitz et al. (2001, JGR 106 2441).
   !
   !! **What the clean-up does, in order** (per cell; the order is the upstream
   !! one, so the answers match bitwise on CPU builds):
   !!
   !! 1. the cell ice area `aice = sum_n aicen` is formed; outside
   !!    `[-puny, 1 + puny]` the cell is refused (`ICE_ERR_AICE_RANGE`);
   !! 2. if `dorebin` and `aice > puny`, `itd_rebin_cell` moves every category
   !!    whose mean thickness left its bin (the partner of the remapping: when
   !!    the remapping is abandoned for a cell it leaves the state untouched and
   !!    this step rebins);
   !! 3. `zap_small_areas`: a category with `|aicen| <= puny` that still carries
   !!    area, ice volume or snow volume is removed, and ALL categories are
   !!    removed when the cell is almost ice-free (`aice < max(itd_area_min,
   !!    puny)` or `aice*rhoi < max(itd_mass_min, puny)`).  A category with
   !!    `aicen < -puny` is refused (`ICE_ERR_AICEN_NEG`).  A round-off excess
   !!    `1 < aice < 1 + puny` scales every category's area, ice volume and snow
   !!    volume by `1/aice` (tracers unchanged), removing the excess to the ocean;
   !! 4. `zap_snow_temperature`: the snow of a category with `aicen > puny` is
   !!    removed when any snow layer is colder than `-100 degC` or warmer than
   !!    the melting point allowed by round-off (`Tmax`, a function of the layer
   !!    enthalpy; a snow cover thinner than `hs_min` is judged on the reference
   !!    enthalpy `-rhos*Lfresh` against `Tmax = puny`);
   !! 5. the zap fluxes are ADDED to the cell flux accumulators;
   !! 6. the cell aggregates are refreshed (`itd_aggregate_cell`).
   !!
   !! **What reaches the ocean.**  A removed category hands its mass, salt and
   !! heat to the ocean, per second and per unit CELL area (a zap is a per-cell
   !! term: not weighted by `aicen_init`), with the Icepack signs of the flux
   !! slot: `fresh += rhoi*vicen/dt` (+ snow `rhos*vsnon/dt`), positive into the
   !! ocean; `fsalt += rhoi*vicen*ice_ref_salinity*0.001/dt`, positive
   !! salinifying (`saltflux_option = "constant"`, the only form the
   !! configuration accepts); `fhocn += sum_k qice(k)*vicen/nilyr/dt` (+ snow
   !! `qsno*vsnon/nslyr/dt`), NEGATIVE because the enthalpies are negative (the
   !! ocean spends heat melting what was removed).  Hence the closure the tests
   !! pin: mass, salt and energy removed from the state equal `dt` times the
   !! increments of `fresh`, `fsalt`, `fhocn`.  The remapping and the rebin move
   !! quantities between categories and emit nothing.
   !!
   !! **Quirks kept verbatim.**  (a) Upstream's `zap_small_areas` updates the open
   !! water `aice0` for a removed category but NOT `aice`; its round-off-excess test
   !! reads the `aice` of the entry (before the rebin).  The port keeps that
   !! decision on the entry `aice` and then recomputes the TRUE aggregates from the
   !! final state (the intermediate `aice`, `aice0` of the packed form are
   !! upstream's).  (b) An empty cell (`aice < itd_area_min`) is "zapped" in all
   !! categories: no flux (the removed quantities are zero) but the category
   !! tracers are normalised (`Tsfc = Tf`, `qice = qsno = sice = 0`).  (c) The
   !! round-off scaling passes `aice - 1` and the entry `aice` of the cell to
   !! every category.
   !!
   !! **Errors.**  Non-finite state or tracer, a non-finite `Tf`, `dt` not finite
   !! and positive, or `ncat > ICE_NCAT_MAX`: `ICE_ERR_NONFINITE`, nothing is
   !! touched, nothing is zapped (a comparison with NaN is false, so a NaN
   !! category would otherwise be neither "small" nor out of range: it is caught
   !! up front instead).  `ICE_ERR_AICE_RANGE` and `ICE_ERR_AICEN_NEG` leave the
   !! state as the rebin left it and add no flux.  All three are FATAL and
   !! UNRECOVERABLE: the kernel works in place, so it cannot "not commit" -- after
   !! a fatal bit the rebinned state is already written, the aggregates are stale
   !! (not refreshed) and the fluxes are not incremented.  The run must stop; the
   !! cell is not to be used.  `aicen < -puny` is checked for every category
   !! BEFORE any zap (upstream interleaves check and zap, and aborts anyway).
   !! The re-check of `aice > 1 + puny` inside upstream's `zap_small_areas` is
   !! unreachable (the entry check used the same `aice`) and is not ported.
   !!
   !! **Forms.**  `itd_cleanup` works on one packed column (`aicen(ncat)`,
   !! `trcrn(ldt, ncat)`; host reference and tests) and returns upstream's
   !! intermediate `aice`, `aice0`.  `ice_cleanup_itd_cell` is the GPU form: cell
   !! `(i, j)` of the SoA state, IN PLACE, no column gathered (a category is
   !! zeroed where it lies), the fluxes added into `fresh(nx,ny)`, `fsalt`,
   !! `fhocn`, the true aggregates refreshed through `itd_aggregate_cell` when
   !! `do_aggregate` is true.  Upstream's clean-up does NOT aggregate (the step
   !! driver aggregates once after the thermodynamics and once after the ridging):
   !! a driver that does so passes `do_aggregate = .false.` and the aggregate
   !! outputs are then left untouched.
   !!
   !! **Two exact shortcuts of the cell form** (bitwise equal to the general path;
   !! the packed form has none, it is the reference).  (1) A cell whose every
   !! category area, ice volume and snow volume is exactly 0 (with `limit_aice`) is
   !! "zapped" in all categories by the general path with zero quantities: no
   !! flux except `x + 0` (a `-0` flux becomes `+0`), tracers normalised.  The
   !! shortcut writes each normalised value only where its bit pattern differs,
   !! and the known aggregates directly.  (2) `itd_rebin_cell` is not called when
   !! no category is out of its bin by the very predicates the rebin uses (its
   !! gather and write-back would be an identity).
   !! `ice_cleanup_itd_step` is the host shim of the one `do concurrent` over the
   !! owned cells.  Packed and cell forms run the same arithmetic in the same
   !! order (bitwise equal on CPU builds).
   !!
   !! Upstream branches dropped as inert for BL99 (`ktherm = 1`), no ponds, no
   !! aerosols, no isotopes, no floe-size distribution, no BGC: the pond water
   !! flux (`tr_pond_topo`, `fpond`), aerosol and isotope fluxes and zeroing
   !! (`tr_aero`, `tr_iso`), the skeleton and layer BGC zaps (`skl_bgc`,
   !! `z_tracers`), `first_ice` (read only by BGC and mushy thermodynamics), the
   !! brine-height reset (`tr_brine`), the snow-density and snow-grain tracer
   !! resets (`snwredist`, `snwgrain`), and `saltflux_option = "prognostic"` (the
   !! configuration refuses it; the salt flux uses the reference salinity).
   !! Deviation: the optional-argument defaults (`limit_aice`, `dorebin`) are
   !! explicit logical arguments.
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use, intrinsic :: iso_fortran_env, only: int32
   use rdb_constants, only: wp
   use rdb_ice_dims, only: ICE_NCAT_MAX
   use rdb_ice_errors, only: ICE_ERR_NONFINITE, ICE_ERR_AICE_RANGE, ICE_ERR_AICEN_NEG
   use rdb_ice_params, only: ice_params_t, ice_ntrcr, ice_nt_qice, ice_nt_qsno, ice_nt_sice, ICE_QNAN
   use rdb_ice_itd_shift, only: itd_rebin, itd_rebin_cell, itd_aggregate_cell
   implicit none
   private

   public :: ITD_ZAP_TMIN
   public :: itd_zap_snow, itd_zap_snow_cell, itd_cleanup, ice_cleanup_itd_cell, ice_cleanup_itd_step

   real(wp), parameter :: ITD_ZAP_TMIN = -100.0_wp
      !! Coldest snow temperature [degC] kept by the snow-temperature zap.
   real(wp), parameter :: P001 = 0.001_wp
      !! Salinity [ppt] to mass fraction, as upstream writes it.

contains

   pure function bits_differ(x, ref) result(differ)
      !! `x` is not the bit pattern of `ref` (finite `x`, `ref`): different value or different sign of zero.
      !$acc routine seq
      real(wp), intent(in) :: x, ref
      logical :: differ

      differ = x /= ref .or. sign(1.0_wp, x) /= sign(1.0_wp, ref)
   end function bits_differ

   pure subroutine itd_zap_snow(nslyr, ldt, iqs, dt, rhos, trcr, vsnon, dfresh, dfhocn)
      !! Remove the snow of one category (packed column `trcr(ldt)`, first snow
      !! enthalpy at index `iqs`): `dfhocn += qsno*vsnon/nslyr/dt` per layer
      !! (negative), `dfresh += rhos*vsnon/dt`; the enthalpies and `vsnon` are zeroed.
      !$acc routine seq
      integer, intent(in) :: nslyr, ldt, iqs
      real(wp), intent(in) :: dt, rhos
      real(wp), intent(inout) :: trcr(ldt)
      real(wp), intent(inout) :: vsnon, dfresh, dfhocn

      integer :: k
      real(wp) :: xtmp

      do k = 1, nslyr
         xtmp = trcr(iqs + k - 1)/dt*vsnon/real(nslyr, wp)
         dfhocn = dfhocn + xtmp
         trcr(iqs + k - 1) = 0.0_wp
      end do
      xtmp = (rhos*vsnon)/dt
      dfresh = dfresh + xtmp
      vsnon = 0.0_wp
   end subroutine itd_zap_snow

   pure subroutine itd_zap_snow_cell(i, j, n, nx, ny, ncat, nslyr, dt, rhos, vsnon, qsnon, dfresh, dfhocn)
      !! In-place twin of `itd_zap_snow` for category `n` of cell `(i, j)` of the
      !! SoA state `vsnon(nx,ny,ncat)`, `qsnon(nx,ny,ncat,nslyr)`.
      !$acc routine seq
      integer, value :: i, j, n
      integer, intent(in) :: nx, ny, ncat, nslyr
      real(wp), intent(in) :: dt, rhos
      real(wp), intent(inout) :: vsnon(nx, ny, ncat), qsnon(nx, ny, ncat, nslyr)
      real(wp), intent(inout) :: dfresh, dfhocn

      integer :: k
      real(wp) :: xtmp

      do k = 1, nslyr
         xtmp = qsnon(i, j, n, k)/dt*vsnon(i, j, n)/real(nslyr, wp)
         dfhocn = dfhocn + xtmp
         qsnon(i, j, n, k) = 0.0_wp
      end do
      xtmp = (rhos*vsnon(i, j, n))/dt
      dfresh = dfresh + xtmp
      vsnon(i, j, n) = 0.0_wp
   end subroutine itd_zap_snow_cell

   pure subroutine itd_cleanup(ncat, ldt, nilyr, nslyr, hin_max, dt, Tf, rhoi, rhos, Lfresh, cp_ice, hs_min, puny, &
                               area_min, mass_min, sref, dorebin, limit_aice, aicen, vicen, vsnon, trcrn, &
                               aice, aice0, fresh, fsalt, fhocn, ierr)
      !! Packed clean-up of one column (see the module note).  `aice`, `aice0`
      !! return upstream's values (`aice0` raised by every removed area, `aice` the
      !! entry sum or 1 after a round-off scaling); `fresh`, `fsalt`, `fhocn` are
      !! incremented by the zap fluxes.  `ldt` is the leading extent of `trcrn`.
      !$acc routine seq
      integer, intent(in) :: ncat, ldt, nilyr, nslyr
      real(wp), intent(in) :: hin_max(0:ncat)
      real(wp), intent(in) :: dt, Tf, rhoi, rhos, Lfresh, cp_ice, hs_min, puny
      real(wp), intent(in) :: area_min, mass_min, sref
      logical, intent(in) :: dorebin, limit_aice
      real(wp), intent(inout) :: aicen(ncat), vicen(ncat), vsnon(ncat), trcrn(ldt, ncat)
      real(wp), intent(out) :: aice, aice0
      real(wp), intent(inout) :: fresh, fsalt, fhocn
      integer(int32), intent(inout) :: ierr

      integer :: n, k, it, ntr, iq, iqs, isi
      logical :: bad, zap_residual, zap_cat, l_zap
      integer(int32) :: ierr0, ierr_r
      real(wp) :: dfresh, dfsalt, dfhocn, xtmp, hsn, zqsn, zTsn, Tmax, rnslyr

      ntr = ice_ntrcr(nilyr, nslyr)
      iq = ice_nt_qice()
      iqs = ice_nt_qsno(nilyr)
      isi = ice_nt_sice(nilyr, nslyr)
      bad = ncat < 1 .or. ncat > ICE_NCAT_MAX .or. nilyr < 1 .or. nslyr < 1 .or. ntr > ldt &
            .or. .not. ieee_is_finite(Tf) .or. .not. ieee_is_finite(dt)
      if (.not. bad) bad = .not. (dt > 0.0_wp)
      if (.not. bad) then
         do n = 1, ncat
            if (.not. (ieee_is_finite(aicen(n)) .and. ieee_is_finite(vicen(n)) .and. ieee_is_finite(vsnon(n)))) bad = .true.
            do it = 1, ntr
               if (.not. ieee_is_finite(trcrn(it, n))) bad = .true.
            end do
         end do
      end if
      if (bad) then
         ierr = ior(ierr, ICE_ERR_NONFINITE)
         aice = ICE_QNAN
         aice0 = ICE_QNAN
         return
      end if

      dfresh = 0.0_wp
      dfsalt = 0.0_wp
      dfhocn = 0.0_wp
      aice = 0.0_wp
      do n = 1, ncat
         aice = aice + aicen(n)
      end do
      aice0 = max(1.0_wp - aice, 0.0_wp)

      if (limit_aice) then
         if (aice > 1.0_wp + puny .or. aice < -puny) then
            ierr = ior(ierr, ICE_ERR_AICE_RANGE)
            return
         end if
      end if

      if (dorebin .and. aice > puny) then
         ierr0 = ierr
         ierr_r = 0_int32
         call itd_rebin(ncat, ldt, nilyr, nslyr, hin_max, puny, Tf, aicen, vicen, vsnon, trcrn, ierr_r)
         ierr = ior(ierr0, ierr_r)
         if (ierr_r /= 0_int32) return
      end if

      if (limit_aice) then
         do n = 1, ncat
            if (aicen(n) < -puny) then
               ierr = ior(ierr, ICE_ERR_AICEN_NEG)
               return
            end if
         end do
         zap_residual = .false.
         if (aice < max(area_min, puny) .or. aice*rhoi < max(mass_min, puny)) zap_residual = .true.
         do n = 1, ncat
            zap_cat = .false.
            if (abs(aicen(n)) <= puny .and. &
                (abs(aicen(n)) /= 0.0_wp .or. abs(vicen(n)) /= 0.0_wp .or. abs(vsnon(n)) /= 0.0_wp)) zap_cat = .true.
            if (zap_cat .or. zap_residual) then
               do k = 1, nilyr
                  xtmp = trcrn(iq + k - 1, n)/dt*vicen(n)/real(nilyr, wp)
                  dfhocn = dfhocn + xtmp
                  trcrn(iq + k - 1, n) = 0.0_wp
               end do
               xtmp = (rhoi*vicen(n))/dt
               dfresh = dfresh + xtmp
               xtmp = rhoi*vicen(n)*sref*P001/dt
               dfsalt = dfsalt + xtmp
               aice0 = aice0 + aicen(n)
               aicen(n) = 0.0_wp
               vicen(n) = 0.0_wp
               trcrn(1, n) = Tf
               call itd_zap_snow(nslyr, ldt, iqs, dt, rhos, trcrn(1:ldt, n), vsnon(n), dfresh, dfhocn)
               do k = 1, nilyr
                  trcrn(isi + k - 1, n) = 0.0_wp
               end do
            end if
         end do

         ! round-off excess of the cell area: every category loses the share (aice - 1)/aice
         ! (upstream re-checks aice > 1 + puny here; unreachable, the entry check used the same aice)
         if (aice > 1.0_wp .and. aice < 1.0_wp + puny) then
            do n = 1, ncat
               do k = 1, nilyr
                  xtmp = trcrn(iq + k - 1, n)*vicen(n)/real(nilyr, wp)*(aice - 1.0_wp)/aice/dt
                  dfhocn = dfhocn + xtmp
               end do
               do k = 1, nslyr
                  xtmp = trcrn(iqs + k - 1, n)*vsnon(n)/real(nslyr, wp)*(aice - 1.0_wp)/aice/dt
                  dfhocn = dfhocn + xtmp
               end do
               xtmp = (rhoi*vicen(n) + rhos*vsnon(n))*(aice - 1.0_wp)/aice/dt
               dfresh = dfresh + xtmp
               xtmp = rhoi*vicen(n)*sref*P001*(aice - 1.0_wp)/aice/dt
               dfsalt = dfsalt + xtmp
               aicen(n) = aicen(n)*(1.0_wp/aice)
               vicen(n) = vicen(n)*(1.0_wp/aice)
               vsnon(n) = vsnon(n)*(1.0_wp/aice)
            end do
            aice = 1.0_wp
            aice0 = 0.0_wp
         end if
      end if

      rnslyr = real(nslyr, wp)
      do n = 1, ncat
         l_zap = .false.
         if (aicen(n) > puny) then
            hsn = vsnon(n)/aicen(n)
            do k = 1, nslyr
               if (hsn > hs_min) then
                  zqsn = trcrn(iqs + k - 1, n)
                  Tmax = -zqsn*puny*rnslyr/(rhos*cp_ice*vsnon(n))
               else
                  zqsn = -rhos*Lfresh
                  Tmax = puny
               end if
               zTsn = (Lfresh + zqsn/rhos)/cp_ice
               if (zTsn < ITD_ZAP_TMIN .or. zTsn > Tmax) l_zap = .true.
            end do
         end if
         if (l_zap) call itd_zap_snow(nslyr, ldt, iqs, dt, rhos, trcrn(1:ldt, n), vsnon(n), dfresh, dfhocn)
      end do

      fresh = fresh + dfresh
      fsalt = fsalt + dfsalt
      fhocn = fhocn + dfhocn
   end subroutine itd_cleanup

   pure subroutine ice_cleanup_itd_cell(i, j, nx, ny, ncat, nilyr, nslyr, hin_max, dt, Tf, rhoi, rhos, Lfresh, cp_ice, &
                                        hs_min, puny, area_min, mass_min, sref, dorebin, limit_aice, do_aggregate, &
                                        aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, &
                                        aice, vice, vsno, aice0, Tsfc_c, qice_c, sice_c, qsno_c, &
                                        fresh, fsalt, fhocn, ierr)
      !! GPU cell form of `itd_cleanup` for cell `(i, j)` of the SoA state
      !! (`aicen, vicen, vsnon (nx,ny,ncat)`, `Tsfcn`, `qicen, sicen
      !! (nx,ny,ncat,nilyr)`, `qsnon (nx,ny,ncat,nslyr)`), IN PLACE.  On success
      !! the zap fluxes are added into `fresh(i,j)`, `fsalt(i,j)`, `fhocn(i,j)`
      !! and the true aggregates are written to `aice, vice, vsno, aice0 (i,j)`
      !! and the cell-mean tracers `Tsfc_c, qice_c, sice_c, qsno_c` (the packed
      !! order of `itd_aggregate_cell`).  See the module note for the errors.
      !$acc routine seq
      integer, value :: i, j
      logical, value :: dorebin, limit_aice, do_aggregate
      integer, intent(in) :: nx, ny, ncat, nilyr, nslyr
      real(wp), intent(in) :: hin_max(0:ncat)
      real(wp), intent(in) :: dt, Tf, rhoi, rhos, Lfresh, cp_ice, hs_min, puny
      real(wp), intent(in) :: area_min, mass_min, sref
      real(wp), intent(inout) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat)
      real(wp), intent(inout) :: Tsfcn(nx, ny, ncat), qicen(nx, ny, ncat, nilyr)
      real(wp), intent(inout) :: sicen(nx, ny, ncat, nilyr), qsnon(nx, ny, ncat, nslyr)
      real(wp), intent(inout) :: aice(nx, ny), vice(nx, ny), vsno(nx, ny), aice0(nx, ny)
      real(wp), intent(inout) :: Tsfc_c(nx, ny), qice_c(nx, ny, nilyr), sice_c(nx, ny, nilyr), qsno_c(nx, ny, nslyr)
      real(wp), intent(inout) :: fresh(nx, ny), fsalt(nx, ny), fhocn(nx, ny)
      integer(int32), intent(inout) :: ierr(nx, ny)

      integer :: n, k
      logical :: bad, zap_residual, zap_cat, l_zap, empty, need_rebin
      integer(int32) :: ierr0, ierr_r
      real(wp) :: dfresh, dfsalt, dfhocn, xtmp, hsn, zqsn, zTsn, Tmax, rnslyr
      real(wp) :: aic, aic0, aic_l, vic_l, vso_l, aic0_l

      bad = ncat < 1 .or. ncat > ICE_NCAT_MAX .or. nilyr < 1 .or. nslyr < 1 &
            .or. .not. ieee_is_finite(Tf) .or. .not. ieee_is_finite(dt)
      if (.not. bad) bad = .not. (dt > 0.0_wp)
      if (.not. bad) then
         do n = 1, ncat
            if (.not. (ieee_is_finite(aicen(i, j, n)) .and. ieee_is_finite(vicen(i, j, n)) .and. &
                       ieee_is_finite(vsnon(i, j, n)) .and. ieee_is_finite(Tsfcn(i, j, n)))) bad = .true.
            do k = 1, nilyr
               if (.not. (ieee_is_finite(qicen(i, j, n, k)) .and. ieee_is_finite(sicen(i, j, n, k)))) bad = .true.
            end do
            do k = 1, nslyr
               if (.not. ieee_is_finite(qsnon(i, j, n, k))) bad = .true.
            end do
         end do
      end if
      if (bad) then
         ierr(i, j) = ior(ierr(i, j), ICE_ERR_NONFINITE)
         aice(i, j) = ICE_QNAN
         vice(i, j) = ICE_QNAN
         vsno(i, j) = ICE_QNAN
         aice0(i, j) = ICE_QNAN
         return
      end if

      if (limit_aice) then
         empty = .true.
         do n = 1, ncat
            if (aicen(i, j, n) /= 0.0_wp .or. vicen(i, j, n) /= 0.0_wp .or. vsnon(i, j, n) /= 0.0_wp) empty = .false.
         end do
         if (empty) then
            ! all-empty cell: the general path zaps every category with zero quantities (see the module note)
            do n = 1, ncat
               if (bits_differ(aicen(i, j, n), 0.0_wp)) aicen(i, j, n) = 0.0_wp
               if (bits_differ(vicen(i, j, n), 0.0_wp)) vicen(i, j, n) = 0.0_wp
               if (bits_differ(vsnon(i, j, n), 0.0_wp)) vsnon(i, j, n) = 0.0_wp
               if (bits_differ(Tsfcn(i, j, n), Tf)) Tsfcn(i, j, n) = Tf
               do k = 1, nilyr
                  if (bits_differ(qicen(i, j, n, k), 0.0_wp)) qicen(i, j, n, k) = 0.0_wp
                  if (bits_differ(sicen(i, j, n, k), 0.0_wp)) sicen(i, j, n, k) = 0.0_wp
               end do
               do k = 1, nslyr
                  if (bits_differ(qsnon(i, j, n, k), 0.0_wp)) qsnon(i, j, n, k) = 0.0_wp
               end do
            end do
            ! x + (+0): only a negative zero changes (written with abs: `x + 0.0` may be folded away by -fast)
            if (fresh(i, j) == 0.0_wp) fresh(i, j) = abs(fresh(i, j))
            if (fsalt(i, j) == 0.0_wp) fsalt(i, j) = abs(fsalt(i, j))
            if (fhocn(i, j) == 0.0_wp) fhocn(i, j) = abs(fhocn(i, j))
            if (do_aggregate) then
               aice(i, j) = 0.0_wp
               vice(i, j) = 0.0_wp
               vsno(i, j) = 0.0_wp
               aice0(i, j) = 1.0_wp
               Tsfc_c(i, j) = Tf
               do k = 1, nilyr
                  qice_c(i, j, k) = 0.0_wp
                  sice_c(i, j, k) = 0.0_wp
               end do
               do k = 1, nslyr
                  qsno_c(i, j, k) = 0.0_wp
               end do
            end if
            return
         end if
      end if

      dfresh = 0.0_wp
      dfsalt = 0.0_wp
      dfhocn = 0.0_wp
      aic = 0.0_wp
      do n = 1, ncat
         aic = aic + aicen(i, j, n)
      end do
      aic0 = max(1.0_wp - aic, 0.0_wp)

      if (limit_aice) then
         if (aic > 1.0_wp + puny .or. aic < -puny) then
            ierr(i, j) = ior(ierr(i, j), ICE_ERR_AICE_RANGE)
            return
         end if
      end if

      need_rebin = .false.
      if (aicen(i, j, 1) > puny) then
         if (vicen(i, j, 1)/aicen(i, j, 1) <= hin_max(0) .and. hin_max(0) > 0.0_wp) need_rebin = .true.
      end if
      do n = 1, ncat - 1
         if (aicen(i, j, n) > puny) then
            if (vicen(i, j, n)/aicen(i, j, n) > hin_max(n)) need_rebin = .true.
         end if
         if (aicen(i, j, n + 1) > puny) then
            if (vicen(i, j, n + 1)/aicen(i, j, n + 1) <= hin_max(n)) need_rebin = .true.
         end if
      end do
      if (dorebin .and. aic > puny .and. need_rebin) then
         ierr0 = ierr(i, j)
         ierr(i, j) = 0_int32
         call itd_rebin_cell(i, j, nx, ny, ncat, nilyr, nslyr, hin_max, puny, Tf, &
                             aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, ierr)
         ierr_r = ierr(i, j)
         ierr(i, j) = ior(ierr0, ierr_r)
         if (ierr_r /= 0_int32) return
      end if

      if (limit_aice) then
         do n = 1, ncat
            if (aicen(i, j, n) < -puny) then
               ierr(i, j) = ior(ierr(i, j), ICE_ERR_AICEN_NEG)
               return
            end if
         end do
         zap_residual = .false.
         if (aic < max(area_min, puny) .or. aic*rhoi < max(mass_min, puny)) zap_residual = .true.
         do n = 1, ncat
            zap_cat = .false.
            if (abs(aicen(i, j, n)) <= puny .and. &
                (abs(aicen(i, j, n)) /= 0.0_wp .or. abs(vicen(i, j, n)) /= 0.0_wp .or. &
                 abs(vsnon(i, j, n)) /= 0.0_wp)) zap_cat = .true.
            if (zap_cat .or. zap_residual) then
               do k = 1, nilyr
                  xtmp = qicen(i, j, n, k)/dt*vicen(i, j, n)/real(nilyr, wp)
                  dfhocn = dfhocn + xtmp
                  qicen(i, j, n, k) = 0.0_wp
               end do
               xtmp = (rhoi*vicen(i, j, n))/dt
               dfresh = dfresh + xtmp
               xtmp = rhoi*vicen(i, j, n)*sref*P001/dt
               dfsalt = dfsalt + xtmp
               aic0 = aic0 + aicen(i, j, n)
               aicen(i, j, n) = 0.0_wp
               vicen(i, j, n) = 0.0_wp
               Tsfcn(i, j, n) = Tf
               call itd_zap_snow_cell(i, j, n, nx, ny, ncat, nslyr, dt, rhos, vsnon, qsnon, dfresh, dfhocn)
               do k = 1, nilyr
                  sicen(i, j, n, k) = 0.0_wp
               end do
            end if
         end do

         ! (upstream re-checks aice > 1 + puny here; unreachable, the entry check used the same aice)
         if (aic > 1.0_wp .and. aic < 1.0_wp + puny) then
            do n = 1, ncat
               do k = 1, nilyr
                  xtmp = qicen(i, j, n, k)*vicen(i, j, n)/real(nilyr, wp)*(aic - 1.0_wp)/aic/dt
                  dfhocn = dfhocn + xtmp
               end do
               do k = 1, nslyr
                  xtmp = qsnon(i, j, n, k)*vsnon(i, j, n)/real(nslyr, wp)*(aic - 1.0_wp)/aic/dt
                  dfhocn = dfhocn + xtmp
               end do
               xtmp = (rhoi*vicen(i, j, n) + rhos*vsnon(i, j, n))*(aic - 1.0_wp)/aic/dt
               dfresh = dfresh + xtmp
               xtmp = rhoi*vicen(i, j, n)*sref*P001*(aic - 1.0_wp)/aic/dt
               dfsalt = dfsalt + xtmp
               aicen(i, j, n) = aicen(i, j, n)*(1.0_wp/aic)
               vicen(i, j, n) = vicen(i, j, n)*(1.0_wp/aic)
               vsnon(i, j, n) = vsnon(i, j, n)*(1.0_wp/aic)
            end do
         end if
      end if

      rnslyr = real(nslyr, wp)
      do n = 1, ncat
         l_zap = .false.
         if (aicen(i, j, n) > puny) then
            hsn = vsnon(i, j, n)/aicen(i, j, n)
            do k = 1, nslyr
               if (hsn > hs_min) then
                  zqsn = qsnon(i, j, n, k)
                  Tmax = -zqsn*puny*rnslyr/(rhos*cp_ice*vsnon(i, j, n))
               else
                  zqsn = -rhos*Lfresh
                  Tmax = puny
               end if
               zTsn = (Lfresh + zqsn/rhos)/cp_ice
               if (zTsn < ITD_ZAP_TMIN .or. zTsn > Tmax) l_zap = .true.
            end do
         end if
         if (l_zap) call itd_zap_snow_cell(i, j, n, nx, ny, ncat, nslyr, dt, rhos, vsnon, qsnon, dfresh, dfhocn)
      end do

      fresh(i, j) = fresh(i, j) + dfresh
      fsalt(i, j) = fsalt(i, j) + dfsalt
      fhocn(i, j) = fhocn(i, j) + dfhocn

      if (do_aggregate) then
         call itd_aggregate_cell(i, j, nx, ny, ncat, nilyr, nslyr, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, Tf, puny, &
                                 aic_l, vic_l, vso_l, aic0_l, Tsfc_c, qice_c, sice_c, qsno_c, ierr)
         aice(i, j) = aic_l
         vice(i, j) = vic_l
         vsno(i, j) = vso_l
         aice0(i, j) = aic0_l
      end if
   end subroutine ice_cleanup_itd_cell

   subroutine ice_cleanup_itd_step(nx, ny, ncat, nilyr, nslyr, i_lo, i_hi, j_lo, j_hi, p, hin_max, dt, Tf, &
                                   dorebin, limit_aice, do_aggregate, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, &
                                   aice, vice, vsno, aice0, Tsfc_c, qice_c, sice_c, qsno_c, &
                                   fresh, fsalt, fhocn, ierr)
      !! Host shim of the clean-up sweep over the owned cells `i_lo:i_hi, j_lo:j_hi`
      !! (every cell, ice or not: an empty cell only has its tracers normalised).
      !! `Tf(nx, ny)` is the freezing temperature per cell, `dt` the thermodynamic
      !! step [s].  Every array must be device-present on a GPU build; the parameter
      !! bundle is read on the host (scalars are passed into the kernel).  Cells
      !! outside the range are not touched.
      integer, intent(in) :: nx, ny, ncat, nilyr, nslyr, i_lo, i_hi, j_lo, j_hi
      type(ice_params_t), intent(in) :: p
      real(wp), intent(in) :: hin_max(0:ncat), Tf(nx, ny)
      real(wp), intent(in) :: dt
      logical, intent(in) :: dorebin, limit_aice, do_aggregate
      real(wp), intent(inout) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat)
      real(wp), intent(inout) :: Tsfcn(nx, ny, ncat), qicen(nx, ny, ncat, nilyr)
      real(wp), intent(inout) :: sicen(nx, ny, ncat, nilyr), qsnon(nx, ny, ncat, nslyr)
      real(wp), intent(inout) :: aice(nx, ny), vice(nx, ny), vsno(nx, ny), aice0(nx, ny)
      real(wp), intent(inout) :: Tsfc_c(nx, ny), qice_c(nx, ny, nilyr), sice_c(nx, ny, nilyr), qsno_c(nx, ny, nslyr)
      real(wp), intent(inout) :: fresh(nx, ny), fsalt(nx, ny), fhocn(nx, ny)
      integer(int32), intent(inout) :: ierr(nx, ny)

      call cleanup_itd_impl(nx, ny, ncat, nilyr, nslyr, i_lo, i_hi, j_lo, j_hi, dt, p%c%rhoi, p%c%rhos, p%c%Lfresh, &
                            p%c%cp_ice, p%c%hs_min, p%c%puny, p%c%itd_area_min, p%c%itd_mass_min, &
                            p%thermo%ice_ref_salinity, dorebin, limit_aice, do_aggregate, hin_max, Tf, &
                            aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, &
                            aice, vice, vsno, aice0, Tsfc_c, qice_c, sice_c, qsno_c, fresh, fsalt, fhocn, ierr)
   end subroutine ice_cleanup_itd_step

   subroutine cleanup_itd_impl(nx, ny, ncat, nilyr, nslyr, i_lo, i_hi, j_lo, j_hi, dt, rhoi, rhos, Lfresh, cp_ice, &
                               hs_min, puny, area_min, mass_min, sref, dorebin, limit_aice, do_aggregate, hin_max, Tf, &
                               aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, &
                               aice, vice, vsno, aice0, Tsfc_c, qice_c, sice_c, qsno_c, fresh, fsalt, fhocn, ierr)
      integer, intent(in) :: nx, ny, ncat, nilyr, nslyr, i_lo, i_hi, j_lo, j_hi
      real(wp), intent(in) :: dt, rhoi, rhos, Lfresh, cp_ice, hs_min, puny, area_min, mass_min, sref
      logical, intent(in) :: dorebin, limit_aice, do_aggregate
      real(wp), intent(in) :: hin_max(0:ncat), Tf(nx, ny)
      real(wp), intent(inout) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat)
      real(wp), intent(inout) :: Tsfcn(nx, ny, ncat), qicen(nx, ny, ncat, nilyr)
      real(wp), intent(inout) :: sicen(nx, ny, ncat, nilyr), qsnon(nx, ny, ncat, nslyr)
      real(wp), intent(inout) :: aice(nx, ny), vice(nx, ny), vsno(nx, ny), aice0(nx, ny)
      real(wp), intent(inout) :: Tsfc_c(nx, ny), qice_c(nx, ny, nilyr), sice_c(nx, ny, nilyr), qsno_c(nx, ny, nslyr)
      real(wp), intent(inout) :: fresh(nx, ny), fsalt(nx, ny), fhocn(nx, ny)
      integer(int32), intent(inout) :: ierr(nx, ny)
      integer :: i, j

      !$acc data present(hin_max, Tf, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, aice, vice, vsno, aice0) &
      !$acc& present(Tsfc_c, qice_c, sice_c, qsno_c, fresh, fsalt, fhocn, ierr)
      do concurrent(j=j_lo:j_hi, i=i_lo:i_hi)
         call ice_cleanup_itd_cell(i, j, nx, ny, ncat, nilyr, nslyr, hin_max, dt, Tf(i, j), rhoi, rhos, Lfresh, cp_ice, &
                                   hs_min, puny, area_min, mass_min, sref, dorebin, limit_aice, do_aggregate, &
                                   aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, &
                                   aice, vice, vsno, aice0, Tsfc_c, qice_c, sice_c, qsno_c, &
                                   fresh, fsalt, fhocn, ierr)
      end do
      !$acc end data
   end subroutine cleanup_itd_impl

end module rdb_ice_itd_cleanup
