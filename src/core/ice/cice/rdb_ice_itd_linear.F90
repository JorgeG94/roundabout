!! Linear remapping of the ice-thickness distribution after thermodynamics.
module rdb_ice_itd_linear
   ! SPDX-License-Identifier: BSD-3-Clause
   !
   ! Ported from CICE/Icepack (CICE Consortium, https://github.com/CICE-Consortium).
   ! Upstream: Icepack 9a1f3e4 (Icepack 1.5.3+14), module icepack_therm_itd
   ! (linear_itd).
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
   !! **Remapping.**  After the vertical thermodynamics every category `n` has
   !! changed its mean thickness from `hicen_init(n) = vicen_init/aicen_init`
   !! to `hicen(n)` by `dhicen(n)`.  `itd_linear_plan` moves the category
   !! boundaries `hbnew(n) = hin_max(n) + (growth rate at the boundary)` by
   !! interpolating `dhicen` linearly in the initial thickness between the two
   !! neighbouring categories (Lipscomb 2001), fits `g(h)` (area per unit
   !! thickness; `itd_fit_line`) in every category on the NEW boundaries, and
   !! integrates `g` and `h*g` over the part of each new category that now lies
   !! across a FIXED boundary `hin_max(n)`: that is the area `daice(n)` and
   !! volume `dvice(n)` handed from `donor(n)` to the other side.  The
   !! category-1 treatment is the one of upstream: a melting category 1
   !! (`dhicen(1) < 0`) loses the area of the sliver of `g(h)` below `|dh0|`
   !! (volume conserved, so the remaining ice thickens; no flux), a growing one
   !! lifts the lower edge `hbnew(0) = min(dh0, hin_max(1))`.  The transfers are
   !! trimmed (nothing for increments below `puny` of the donor, the whole
   !! category above `1 - puny`) and applied by `itd_shift_ice` with the snow
   !! enthalpy offset `rhos*Lfresh`.  Last, `hi_min`: a category 1 thinner than
   !! `hi_min` gives up area (volume conserved) until its thickness is `hi_min`.
   !!
   !! **Fallback.**  The remap is abandoned for the cell (state untouched,
   !! `ICE_ERR_REMAP_FALLBACK` counted) when a new boundary no longer lies between the
   !! mean thicknesses of its two categories or outside `[hin_max(n-1),
   !! hin_max(n+1)]` (very fast growth or melt), or when the initial mean
   !! thicknesses of two populated neighbours are not increasing
   !! (`ICE_ERR_ITD_HBNEW_ORDER`, which upstream aborts on; counted unless the
   !! fatal knob folds it in).  The caller then rebins (`itd_rebin`) -- the
   !! post-thermo cleanup does that for every cell -- or asks for it here with
   !! `rebin_fb = .true.`.
   !!
   !! **Validity of the output.**  The remapped state is valid only AFTER the post-thermo
   !! cleanup.  A category drained by a transfer carries an ill-conditioned intensive tracer
   !! (a small difference of two large extensive products divided by a small area) until the
   !! small-area zap removes it, and a fallback cell may have categories outside their bins
   !! until it is rebinned.  The therm2 driver PR must therefore test that a fallback cell ends
   !! inside its bins after the cleanup.  An emptied category's snow enthalpy is `-rhos*Lfresh`
   !! after a remap (upstream's offset) but 0 after a rebin.  An out-of-range `ncat` (`< 1` or
   !! `> ICE_NCAT_MAX`) is reported under `ICE_ERR_NONFINITE`, like every other invalid input.
   !!
   !! **Errors and state.**  A non-finite input (area, volume, tracer, `Tf`)
   !! ORs `ICE_ERR_NONFINITE` and leaves the state untouched.  A refused shift
   !! (`ICE_ERR_ITD_SHIFT_RANGE`, fatal) leaves the cell UNSPECIFIED (the cell
   !! form happens to leave it untouched, the packed form may already carry the
   !! category-1 melt loss; callers must not rely on either) and `aice`, `aice0`
   !! NaN.  No output is double-buffered.
   !!
   !! **Forms.**  `itd_linear_itd` works on one packed column
   !! (`aicen(ncat)`, `trcrn(ldt, ncat)`; host reference and tests).
   !! `ice_linear_itd_cell` is the GPU form: cell `(i, j)` of the SoA state
   !! arrays, in place, tracers through `ice_trcr_get/set`, three `ncat`
   !! vectors (`aicen, vicen, vsnon`) of per-thread scratch plus the plan
   !! (`ICE_NCAT_MAX`-sized locals, no automatic arrays); both are
   !! `pure` and `!$acc routine seq`.  `ice_linear_itd_step` is the host shim
   !! of the flat `do concurrent` over the owned cells.  The two forms run the
   !! same operations in the same order (bitwise equal on CPU builds).
   !!
   !! Upstream branches dropped as inert for BL99 (`ktherm = 1`), no ponds
   !! (`fpond`, `dpnd_melt`), no snow-density tracer (`snwredist`, `rhosmin`),
   !! no brine height, no strata tracers; and the optional `conserv_check`
   !! (reproduced as test ledgers).  Deviation: upstream leaves `hbnew(n)`
   !! undefined when the order check aborts; here it is `hin_max(n)` and the
   !! remap is abandoned.
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use, intrinsic :: iso_fortran_env, only: int32
   use rdb_constants, only: wp
   use rdb_ice_dims, only: ICE_NCAT_MAX
   use rdb_ice_errors, only: ICE_ERR_NONFINITE, ICE_ERR_ITD_HBNEW_ORDER, ICE_ERR_ITD_SHIFT_RANGE, &
                             ICE_ERR_REMAP_FALLBACK
   use rdb_ice_params, only: ice_params_t, ice_ntrcr, ice_trcr_get, ICE_QNAN
   use rdb_ice_itd_shift, only: itd_fit_line, itd_aggregate_area, itd_shift_ice, itd_shift_ice_cell, &
                                itd_rebin, itd_rebin_cell
   implicit none
   private

   public :: itd_linear_plan, itd_linear_itd, ice_linear_itd_cell, ice_linear_itd_step

   real(wp), parameter :: P333 = 1.0_wp/3.0_wp
      !! One third as `icepack_parameters` defines `p333` (`c1/c3`).

contains

   pure subroutine itd_linear_plan(ncat, hin_max, puny, aicen_init, vicen_init, aicen, vicen, &
                                   hicen, remap, donor, daice, dvice, ierr)
      !! Boundary motion and transfers of the linear remapping on packed
      !! `ncat`-vectors (see the module note).  `aicen(1)` is reduced by the
      !! melt-area loss of category 1 (the only state it modifies); `hicen` is
      !! the post-thermo thickness (category 1 after that loss).  `remap =
      !! .false.`: nothing is to be shifted (`donor`, `daice`, `dvice` are 0).
      !! `ncat` outside `1..ICE_NCAT_MAX` ORs `ICE_ERR_NONFINITE` and sets
      !! `remap = .false.`.  `ICE_ERR_ITD_HBNEW_ORDER` is OR-ed for a violated
      !! initial-thickness order.
      !$acc routine seq
      integer, intent(in) :: ncat
      real(wp), intent(in) :: puny
      real(wp), intent(in) :: hin_max(0:ncat)
      real(wp), intent(in) :: aicen_init(ncat), vicen_init(ncat), vicen(ncat)
      real(wp), intent(inout) :: aicen(ncat)
      real(wp), intent(out) :: hicen(ncat)
      logical, intent(out) :: remap
      integer, intent(out) :: donor(ncat)
      real(wp), intent(out) :: daice(ncat), dvice(ncat)
      integer(int32), intent(inout) :: ierr

      integer :: n, nd
      real(wp) :: slope, dh0, da0, damax, etamin, etamax, x1, x2, x3, wk1, wk2
      real(wp) :: hbnew(0:ICE_NCAT_MAX)
      real(wp) :: g0(ICE_NCAT_MAX), g1(ICE_NCAT_MAX), hL(ICE_NCAT_MAX), hR(ICE_NCAT_MAX)
      real(wp) :: hicen_init(ICE_NCAT_MAX), dhicen(ICE_NCAT_MAX)

      remap = .false.
      do n = 1, ncat
         donor(n) = 0
         daice(n) = 0.0_wp
         dvice(n) = 0.0_wp
         hicen(n) = 0.0_wp
      end do
      if (ncat < 1 .or. ncat > ICE_NCAT_MAX) then
         ierr = ior(ierr, ICE_ERR_NONFINITE)
         return
      end if
      remap = .true.

      ! thickness change in each category
      do n = 1, ncat
         if (aicen_init(n) > puny) then
            hicen_init(n) = vicen_init(n)/aicen_init(n)
         else
            hicen_init(n) = 0.0_wp
         end if
         if (aicen(n) > puny) then
            hicen(n) = vicen(n)/aicen(n)
            dhicen(n) = hicen(n) - hicen_init(n)
         else
            hicen(n) = 0.0_wp
            dhicen(n) = 0.0_wp
         end if
      end do

      ! new category boundaries from the growth rates
      hbnew(0) = hin_max(0)
      do n = 1, ncat - 1
         if (hicen_init(n) > puny .and. hicen_init(n + 1) > puny) then
            if ((hicen_init(n + 1) - hicen_init(n)) > 0.0_wp) then
               ! interpolate between adjacent category growth rates
               slope = (dhicen(n + 1) - dhicen(n))/(hicen_init(n + 1) - hicen_init(n))
               hbnew(n) = hin_max(n) + dhicen(n) + slope*(hin_max(n) - hicen_init(n))
            else
               ierr = ior(ierr, ICE_ERR_ITD_HBNEW_ORDER)
               remap = .false.
               hbnew(n) = hin_max(n)
            end if
         else if (hicen_init(n) > puny) then
            hbnew(n) = hin_max(n) + dhicen(n)
         else if (hicen_init(n + 1) > puny) then
            hbnew(n) = hin_max(n) + dhicen(n + 1)
         else
            hbnew(n) = hin_max(n)
         end if

         ! each boundary must lie between the mean thicknesses of its categories
         if (aicen(n) > puny .and. hicen(n) >= hbnew(n)) then
            remap = .false.
         else if (aicen(n + 1) > puny .and. hicen(n + 1) <= hbnew(n)) then
            remap = .false.
         end if

         ! and between hin_max(n-1) and hin_max(n+1)
         if (hbnew(n) > hin_max(n + 1)) remap = .false.
         if (hbnew(n) < hin_max(n - 1)) remap = .false.
      end do

      if (aicen(ncat) > puny) then
         hbnew(ncat) = 3.0_wp*hicen(ncat) - 2.0_wp*hbnew(ncat - 1)
      else
         hbnew(ncat) = hin_max(ncat)
      end if
      hbnew(ncat) = max(hbnew(ncat), hin_max(ncat - 1))

      if (.not. remap) return

      ! g(h) of category 1 at the start of the step (hicen = hicen_init)
      call itd_fit_line(aicen(1), hicen_init(1), hbnew(0), hin_max(1), puny, g0(1), g1(1), hL(1), hR(1))

      ! area lost by melting of thin (category 1) ice
      if (aicen(1) > puny) then
         dh0 = dhicen(1)
         if (dh0 < 0.0_wp) then   ! remove area from category 1
            dh0 = min(-dh0, hin_max(1))   ! dh0 --> |dh0|
            ! integrate g(1) from 0 to dh0 (right limit; the left limit is 0)
            etamax = min(dh0, hR(1)) - hL(1)
            if (etamax > 0.0_wp) then
               x1 = etamax
               x2 = 0.5_wp*etamax*etamax
               da0 = g1(1)*x2 + g0(1)*x1   ! ice area removed
               ! constrain the new thickness <= hicen_init
               damax = aicen(1)*(1.0_wp - hicen(1)/hicen_init(1))   ! damax > 0
               da0 = min(da0, damax)
               ! remove area, conserving volume
               hicen(1) = hicen(1)*aicen(1)/(aicen(1) - da0)
               aicen(1) = aicen(1) - da0
            end if
         else                  ! dh0 >= 0
            hbnew(0) = min(dh0, hin_max(1))   ! shift hbnew(0) to the right
         end if
      end if

      ! g(h) of every category on its new boundaries
      do n = 1, ncat
         call itd_fit_line(aicen(n), hicen(n), hbnew(n - 1), hbnew(n), puny, g0(n), g1(n), hL(n), hR(n))
         donor(n) = 0
         daice(n) = 0.0_wp
         dvice(n) = 0.0_wp
      end do

      ! area and volume to be shifted across each boundary
      do n = 1, ncat - 1
         if (hbnew(n) > hin_max(n)) then   ! transfer from n to n+1
            etamin = max(hin_max(n), hL(n)) - hL(n)
            etamax = min(hbnew(n), hR(n)) - hL(n)
            donor(n) = n
         else                              ! transfer from n+1 to n
            etamin = 0.0_wp
            etamax = min(hin_max(n), hR(n + 1)) - hL(n + 1)
            donor(n) = n + 1
         end if

         if (etamax > etamin) then
            x1 = etamax - etamin
            wk1 = etamin*etamin
            wk2 = etamax*etamax
            x2 = 0.5_wp*(wk2 - wk1)
            wk1 = wk1*etamin
            wk2 = wk2*etamax
            x3 = P333*(wk2 - wk1)
            nd = donor(n)
            daice(n) = g1(nd)*x2 + g0(nd)*x1
            dvice(n) = g1(nd)*x3 + g0(nd)*x2 + daice(n)*hL(nd)
         end if

         ! if daice or dvice is very small, shift no ice
         nd = donor(n)
         if (daice(n) < aicen(nd)*puny) then
            daice(n) = 0.0_wp
            dvice(n) = 0.0_wp
            donor(n) = 0
         end if
         if (dvice(n) < vicen(nd)*puny) then
            daice(n) = 0.0_wp
            dvice(n) = 0.0_wp
            donor(n) = 0
         end if

         ! if daice is close to aicen or dvice is close to vicen, shift the entire category
         if (daice(n) > aicen(nd)*(1.0_wp - puny)) then
            daice(n) = aicen(nd)
            dvice(n) = vicen(nd)
         end if
         if (dvice(n) > vicen(nd)*(1.0_wp - puny)) then
            daice(n) = aicen(nd)
            dvice(n) = vicen(nd)
         end if
      end do
   end subroutine itd_linear_plan

   pure subroutine itd_linear_itd(ncat, ldt, nilyr, nslyr, hin_max, aicen_init, vicen_init, &
                                  aicen, vicen, vsnon, trcrn, Tf, rhos, Lfresh, hi_min, puny, rebin_fb, &
                                  aice, aice0, ierr)
      !! Packed column `linear_itd`: `aicen, vicen, vsnon` and `trcrn(ldt, ncat)`
      !! (packed `Tsfc | qice | qsno | sice`) in place after the thermodynamics,
      !! `aicen_init, vicen_init` the start-of-step values.  `aice`, `aice0` are
      !! the cell aggregates after the remap.  See the module note for the
      !! fallback and the error semantics.
      !$acc routine seq
      integer, intent(in) :: ncat, ldt, nilyr, nslyr
      real(wp), intent(in) :: hin_max(0:ncat), aicen_init(ncat), vicen_init(ncat)
      real(wp), intent(inout) :: aicen(ncat), vicen(ncat), vsnon(ncat)
      real(wp), intent(inout) :: trcrn(ldt, ncat)
      real(wp), intent(in) :: Tf, rhos, Lfresh, hi_min, puny
      logical, intent(in) :: rebin_fb
      real(wp), intent(out) :: aice, aice0
      integer(int32), intent(inout) :: ierr

      integer :: n, it, ntr
      integer :: donor(ICE_NCAT_MAX)
      logical :: remap, bad
      real(wp) :: da0
      real(wp) :: daice(ICE_NCAT_MAX), dvice(ICE_NCAT_MAX), hicen(ICE_NCAT_MAX)
      integer(int32) :: ierr_l

      ntr = ice_ntrcr(nilyr, nslyr)
      bad = ntr > ldt .or. ncat < 1 .or. ncat > ICE_NCAT_MAX .or. nilyr < 1 .or. nslyr < 1 &
            .or. .not. ieee_is_finite(Tf)
      if (.not. bad) then
         do n = 1, ncat
            if (.not. (ieee_is_finite(aicen_init(n)) .and. ieee_is_finite(vicen_init(n)) .and. &
                       ieee_is_finite(aicen(n)) .and. ieee_is_finite(vicen(n)) .and. &
                       ieee_is_finite(vsnon(n)))) bad = .true.
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

      call itd_linear_plan(ncat, hin_max, puny, aicen_init, vicen_init, aicen, vicen, &
                           hicen, remap, donor, daice, dvice, ierr)
      if (remap) then
         ierr_l = 0_int32
         call itd_shift_ice(ncat, ldt, nilyr, nslyr, puny, Tf, rhos*Lfresh, donor, daice, dvice, &
                            aicen, vicen, vsnon, hicen, trcrn, ierr_l)
         if (ierr_l /= 0_int32) then
            ierr = ior(ierr, ierr_l)
            aice = ICE_QNAN
            aice0 = ICE_QNAN
            return
         end if
         ! make sure hice(1) >= the minimum ice thickness hi_min
         if (hi_min > 0.0_wp .and. aicen(1) > puny .and. hicen(1) < hi_min) then
            da0 = aicen(1)*(1.0_wp - hicen(1)/hi_min)
            aicen(1) = aicen(1) - da0
            hicen(1) = hi_min
         end if
      else
         ierr = ior(ierr, ICE_ERR_REMAP_FALLBACK)
         if (rebin_fb) call itd_rebin(ncat, ldt, nilyr, nslyr, hin_max, puny, Tf, aicen, vicen, vsnon, trcrn, ierr)
      end if
      call itd_aggregate_area(ncat, aicen, aice, aice0)
   end subroutine itd_linear_itd

   pure subroutine ice_linear_itd_cell(i, j, nx, ny, ncat, nilyr, nslyr, hin_max, aicen_init, vicen_init, &
                                       aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, &
                                       Tf, rhos, Lfresh, hi_min, puny, rebin_fb, aice, aice0, ierr)
      !! GPU cell form of `itd_linear_itd` for cell `(i, j)` of the SoA state
      !! (`aicen, vicen, vsnon (nx,ny,ncat)`, `Tsfcn`, `qicen, sicen
      !! (nx,ny,ncat,nilyr)`, `qsnon (nx,ny,ncat,nslyr)`), updated in place;
      !! `aicen_init`, `vicen_init` the start-of-step snapshots.  The tracers
      !! are read once (finiteness guard) and then updated through the
      !! accessors; the area and volumes are gathered into three `ncat` vectors
      !! and written back once (only when the remap ran).  `ierr(nx, ny)` is the
      !! per-cell error word; `aice`, `aice0` the aggregates.
      !$acc routine seq
      integer, value :: i, j
      logical, value :: rebin_fb
      integer, intent(in) :: nx, ny, ncat, nilyr, nslyr
      real(wp), intent(in) :: hin_max(0:ncat)
      real(wp), intent(in) :: aicen_init(nx, ny, ncat), vicen_init(nx, ny, ncat)
      real(wp), intent(inout) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat)
      real(wp), intent(inout) :: Tsfcn(nx, ny, ncat), qicen(nx, ny, ncat, nilyr)
      real(wp), intent(inout) :: sicen(nx, ny, ncat, nilyr), qsnon(nx, ny, ncat, nslyr)
      real(wp), intent(in) :: Tf, rhos, Lfresh, hi_min, puny
      real(wp), intent(out) :: aice, aice0
      integer(int32), intent(inout) :: ierr(nx, ny)

      integer :: n, it, ntr
      integer :: donor(ICE_NCAT_MAX)
      logical :: remap, bad
      real(wp) :: a(ICE_NCAT_MAX), v(ICE_NCAT_MAX), s(ICE_NCAT_MAX)
      real(wp) :: ai(ICE_NCAT_MAX), vi(ICE_NCAT_MAX)
      real(wp) :: da0
      real(wp) :: daice(ICE_NCAT_MAX), dvice(ICE_NCAT_MAX), hicen(ICE_NCAT_MAX)
      integer(int32) :: ierr_l

      ntr = ice_ntrcr(nilyr, nslyr)
      bad = ncat < 1 .or. ncat > ICE_NCAT_MAX .or. nilyr < 1 .or. nslyr < 1 .or. .not. ieee_is_finite(Tf)
      if (.not. bad) then
         do n = 1, ncat
            if (.not. (ieee_is_finite(aicen_init(i, j, n)) .and. ieee_is_finite(vicen_init(i, j, n)) .and. &
                       ieee_is_finite(aicen(i, j, n)) .and. ieee_is_finite(vicen(i, j, n)) .and. &
                       ieee_is_finite(vsnon(i, j, n)))) bad = .true.
            do it = 1, ntr
               if (.not. ieee_is_finite(ice_trcr_get(it, i, j, n, nx, ny, ncat, nilyr, nslyr, &
                                                     Tsfcn, qicen, sicen, qsnon))) bad = .true.
            end do
         end do
      end if
      if (bad) then
         ierr(i, j) = ior(ierr(i, j), ICE_ERR_NONFINITE)
         aice = ICE_QNAN
         aice0 = ICE_QNAN
         return
      end if

      do n = 1, ncat
         a(n) = aicen(i, j, n)
         v(n) = vicen(i, j, n)
         s(n) = vsnon(i, j, n)
         ai(n) = aicen_init(i, j, n)
         vi(n) = vicen_init(i, j, n)
      end do

      call itd_linear_plan(ncat, hin_max, puny, ai, vi, a, v, hicen, remap, donor, daice, dvice, ierr(i, j))
      if (remap) then
         ierr_l = 0_int32
         call itd_shift_ice_cell(i, j, nx, ny, ncat, nilyr, nslyr, puny, Tf, rhos*Lfresh, donor, daice, dvice, &
                                 a, v, s, hicen, Tsfcn, qicen, sicen, qsnon, ierr_l)
         if (ierr_l /= 0_int32) then
            ierr(i, j) = ior(ierr(i, j), ierr_l)
            aice = ICE_QNAN
            aice0 = ICE_QNAN
            return
         end if
         if (hi_min > 0.0_wp .and. a(1) > puny .and. hicen(1) < hi_min) then
            da0 = a(1)*(1.0_wp - hicen(1)/hi_min)
            a(1) = a(1) - da0
            hicen(1) = hi_min
         end if
         do n = 1, ncat
            aicen(i, j, n) = a(n)
            vicen(i, j, n) = v(n)
            vsnon(i, j, n) = s(n)
         end do
      else
         ierr(i, j) = ior(ierr(i, j), ICE_ERR_REMAP_FALLBACK)
         if (rebin_fb) then
            call itd_rebin_cell(i, j, nx, ny, ncat, nilyr, nslyr, hin_max, puny, Tf, &
                                aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, ierr)
            do n = 1, ncat
               a(n) = aicen(i, j, n)
            end do
         end if
      end if
      call itd_aggregate_area(ncat, a, aice, aice0)
   end subroutine ice_linear_itd_cell

   subroutine ice_linear_itd_step(nx, ny, ncat, nilyr, nslyr, i_lo, i_hi, j_lo, j_hi, p, hin_max, Tf, &
                                  aicen_init, vicen_init, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, &
                                  rebin_fb, aice, aice0, ierr)
      !! Host shim of the remap sweep over the owned cells `i_lo:i_hi, j_lo:j_hi`:
      !! cells with `sum(aicen) > puny` run `ice_linear_itd_cell` (upstream's
      !! gate); the others only get `aice`, `aice0`; a non-finite total flags
      !! `ICE_ERR_NONFINITE`.  `Tf(nx, ny)` is the freezing temperature per cell.
      !! Every array must be device-present on a GPU build; the parameter
      !! bundle is read on the host (scalars are passed into the kernel).
      integer, intent(in) :: nx, ny, ncat, nilyr, nslyr, i_lo, i_hi, j_lo, j_hi
      type(ice_params_t), intent(in) :: p
      real(wp), intent(in) :: hin_max(0:ncat), Tf(nx, ny)
      real(wp), intent(in) :: aicen_init(nx, ny, ncat), vicen_init(nx, ny, ncat)
      real(wp), intent(inout) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat)
      real(wp), intent(inout) :: Tsfcn(nx, ny, ncat), qicen(nx, ny, ncat, nilyr)
      real(wp), intent(inout) :: sicen(nx, ny, ncat, nilyr), qsnon(nx, ny, ncat, nslyr)
      logical, intent(in) :: rebin_fb
      real(wp), intent(inout) :: aice(nx, ny), aice0(nx, ny)
      integer(int32), intent(inout) :: ierr(nx, ny)

      call linear_itd_impl(nx, ny, ncat, nilyr, nslyr, i_lo, i_hi, j_lo, j_hi, &
                           p%c%rhos, p%c%Lfresh, p%thermo%hi_min, p%c%puny, rebin_fb, hin_max, Tf, &
                           aicen_init, vicen_init, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, &
                           aice, aice0, ierr)
   end subroutine ice_linear_itd_step

   subroutine linear_itd_impl(nx, ny, ncat, nilyr, nslyr, i_lo, i_hi, j_lo, j_hi, &
                              rhos, Lfresh, hi_min, puny, rebin_fb, hin_max, Tf, &
                              aicen_init, vicen_init, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, &
                              aice, aice0, ierr)
      integer, intent(in) :: nx, ny, ncat, nilyr, nslyr, i_lo, i_hi, j_lo, j_hi
      real(wp), intent(in) :: rhos, Lfresh, hi_min, puny
      logical, intent(in) :: rebin_fb
      real(wp), intent(in) :: hin_max(0:ncat), Tf(nx, ny)
      real(wp), intent(in) :: aicen_init(nx, ny, ncat), vicen_init(nx, ny, ncat)
      real(wp), intent(inout) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat)
      real(wp), intent(inout) :: Tsfcn(nx, ny, ncat), qicen(nx, ny, ncat, nilyr)
      real(wp), intent(inout) :: sicen(nx, ny, ncat, nilyr), qsnon(nx, ny, ncat, nslyr)
      real(wp), intent(inout) :: aice(nx, ny), aice0(nx, ny)
      integer(int32), intent(inout) :: ierr(nx, ny)
      integer :: i, j, n
      real(wp) :: tot

      do concurrent(j=j_lo:j_hi, i=i_lo:i_hi) local(n, tot)
         tot = 0.0_wp
         do n = 1, ncat
            tot = tot + aicen(i, j, n)
         end do
         if (.not. ieee_is_finite(tot)) then
            ierr(i, j) = ior(ierr(i, j), ICE_ERR_NONFINITE)
            aice(i, j) = tot
            aice0(i, j) = tot
         else if (tot > puny) then
            call ice_linear_itd_cell(i, j, nx, ny, ncat, nilyr, nslyr, hin_max, aicen_init, vicen_init, &
                                     aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, &
                                     Tf(i, j), rhos, Lfresh, hi_min, puny, rebin_fb, aice(i, j), aice0(i, j), ierr)
         else
            aice(i, j) = tot
            aice0(i, j) = max(1.0_wp - tot, 0.0_wp)
         end if
      end do
   end subroutine linear_itd_impl

end module rdb_ice_itd_linear
