!! Thickness-distribution column kernels: line fit (part 1) and cell aggregates.
module rdb_ice_itd_shift
   ! SPDX-License-Identifier: BSD-3-Clause
   !
   ! Ported from CICE/Icepack (CICE Consortium, https://github.com/CICE-Consortium).
   ! Upstream: Icepack 9a1f3e4 (Icepack 1.5.3+14), modules icepack_therm_itd
   ! (line fit) and icepack_itd (aggregate, aggregate_area).
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
   !! **Line fit.**  `itd_fit_line` represents the thickness distribution of one
   !! category on `[hbL, hbR]` as `g(eta) = g0 + g1*eta`, `eta = h - hL`, over a
   !! support `[hL, hR]` inside the category bounds.  The support is the whole
   !! category when its mean thickness lies in the central third, and is cut to
   !! `[hbL, 3*hice - 2*hbL]` (mean in the left third) or
   !! `[3*hice - 2*hbR, hbR]` (right third) otherwise, so that `g >= 0` and the
   !! mean is reproduced: `int g = aicen`, `int h*g = aicen*hice`.
   !!
   !! *Deliberate, answers-neutral deviation from upstream:* a mean exactly at a
   !! bound makes the cut support empty (`hR <= hL`); upstream then divides by
   !! zero and returns NaN `g0`, `g1`.  Here `g0 = g1 = 0` (with `hL`, `hR` as
   !! cut), so every output of a finite input is finite.  Upstream only reads
   !! `g0`, `g1` behind `etamax > etamin`, which cannot hold when `hR <= hL`.
   !! A support collapsed by one ulp (`hR - hL` a few 1e-16) still gives the
   !! huge upstream values, kept verbatim, as the caller's transfer sanity tests
   !! absorb them.  A non-finite INPUT is not absorbed: all four outputs become
   !! non-finite so the error cannot be mistaken for "empty".
   !!
   !! **Aggregates.**  `itd_aggregate_area` and `itd_aggregate` fold the `ncat`
   !! categories of one column into the cell values.  Extensive tracer
   !! `trcr(it)` is `sum_n trcrn(it,n)*w_n / W` where `w_n` is `aicen`, `vicen`
   !! or `vsnon` and `W` its sum, by the tracer's dependence (`ice_trcr_depend`;
   !! one-hot, which is what upstream's `trcr_base` table encodes).  The
   !! accumulation runs tracer-outer, category-inner into ONE scalar, so no
   !! work array exists (upstream allocates one per call).
   !!
   !! All routines are `pure`, `!$acc routine seq`, explicit-shape, with the
   !! category axis the FIRST (contiguous) one: `aicen(ncat)`,
   !! `trcrn(ldt, ncat)`, where `ldt` is the LEADING EXTENT OF THE CALLER'S
   !! ARRAY (for the fixed-size per-thread local, `ICE_NTRCR_MAX`) and the active
   !! tracer count `ice_ntrcr(nilyr, nslyr) <= ldt` is derived inside.  Callers
   !! on `(nx,ny,ncat)` state gather the column into per-thread locals; no array
   !! section is passed.
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
   use rdb_constants, only: wp
   use, intrinsic :: iso_fortran_env, only: int32
   use rdb_ice_errors, only: ICE_ERR_NONFINITE
   use rdb_ice_params, only: ice_trcr_depend, ice_intensive, ice_ntrcr
   implicit none
   private

   public :: itd_fit_line, itd_aggregate_area, itd_aggregate

   real(wp), parameter :: THIRD = 1.0_wp/3.0_wp
      !! One third, multiplied exactly as upstream does (never divided by 3).

contains

   pure subroutine itd_fit_line(aicen, hice, hbL, hbR, puny, g0, g1, hL, hR)
      !! Line fit `g(eta) = g0 + g1*eta` of one category (see the module note).
      !! `aicen` area, `hice` mean thickness, `hbL`/`hbR` category bounds, all
      !! [m] or 1.  Empty (`aicen <= puny`) or too narrow (`hbR - hbL <= puny`)
      !! categories give all four outputs 0.
      !$acc routine seq
      real(wp), intent(in) :: aicen, hice, hbL, hbR, puny
      real(wp), intent(out) :: g0, g1, hL, hR

      real(wp) :: h13, h23, dhr, wk1, wk2, poison

      poison = aicen + hice + hbL + hbR
      if (.not. ieee_is_finite(poison)) then
         g0 = poison
         g1 = poison
         hL = poison
         hR = poison
         return
      end if

      if (aicen > puny .and. hbR - hbL > puny) then
         hL = hbL
         hR = hbR
         ! cut the support when the mean is outside the central third
         h13 = THIRD*(2.0_wp*hL + hR)
         h23 = THIRD*(hL + 2.0_wp*hR)
         if (hice < h13) then
            hR = 3.0_wp*hice - 2.0_wp*hL
         else if (hice > h23) then
            hL = 3.0_wp*hice - 2.0_wp*hR
         end if
         if (hR > hL) then
            dhr = 1.0_wp/(hR - hL)
            wk1 = 6.0_wp*aicen*dhr
            wk2 = (hice - hL)*dhr
            g0 = wk1*(2.0_wp/3.0_wp - wk2)
            g1 = 2.0_wp*dhr*wk1*(wk2 - 0.5_wp)
         else
            g0 = 0.0_wp
            g1 = 0.0_wp
         end if
      else
         g0 = 0.0_wp
         g1 = 0.0_wp
         hL = 0.0_wp
         hR = 0.0_wp
      end if
   end subroutine itd_fit_line

   pure subroutine itd_aggregate_area(ncat, aicen, aice, aice0)
      !! Cell ice area `aice = sum aicen` and open water `aice0 = max(1-aice, 0)`.
      !! A non-finite `aicen` gives a non-finite `aice` and `aice0` (no clamp
      !! launders it).
      !$acc routine seq
      integer, intent(in) :: ncat
      real(wp), intent(in) :: aicen(ncat)
      real(wp), intent(out) :: aice, aice0
      integer :: n

      aice = 0.0_wp
      do n = 1, ncat
         aice = aice + aicen(n)
      end do
      if (ieee_is_finite(aice)) then
         aice0 = max(1.0_wp - aice, 0.0_wp)
      else
         aice0 = aice
      end if
   end subroutine itd_aggregate_area

   pure subroutine itd_aggregate(ncat, ldt, nilyr, nslyr, aicen, vicen, vsnon, trcrn, Tf, puny, &
                                 aice, vice, vsno, aice0, trcr, ierr)
      !! Cell aggregates of one column: `aice`, `vice`, `vsno`, `aice0` and the
      !! intensive tracer column `trcr(1:ntr)`, `ntr = ice_ntrcr(nilyr, nslyr)`
      !! (packed `Tsfc | qice | qsno | sice`).  `ldt` is the leading extent of the
      !! caller's `trcrn` and `trcr` arrays (`ICE_NTRCR_MAX` for a fixed-size
      !! per-thread local); entries `ntr+1:ldt` are not touched.  An area tracer on
      !! a cell with `aice <= puny` takes `Tf` for `Tsfc` and 0 for any other; a
      !! volume tracer on an empty cell is 0.  Any non-finite input (state or
      !! tracer) makes the affected outputs non-finite and ORs `ICE_ERR_NONFINITE`
      !! into `ierr`.  `ntr > ldt` (an array too short for the layer counts) sets
      !! `ICE_ERR_NONFINITE`, makes every output non-finite and returns; no
      !! tracer is ever given a snow-volume weight by default.
      !$acc routine seq
      integer, intent(in) :: ncat, ldt, nilyr, nslyr
      real(wp), intent(in) :: aicen(ncat), vicen(ncat), vsnon(ncat)
      real(wp), intent(in) :: trcrn(ldt, ncat)
      real(wp), intent(in) :: Tf, puny
      real(wp), intent(out) :: aice, vice, vsno, aice0
      real(wp), intent(out) :: trcr(ldt)
      integer(int32), intent(inout) :: ierr

      integer :: n, it, d, ntr
      real(wp) :: atr, w, empty, poison

      ntr = ice_ntrcr(nilyr, nslyr)
      if (ntr > ldt .or. ntr < 1 .or. ncat < 1) then
         poison = ieee_value(1.0_wp, ieee_quiet_nan)
         aice = poison
         vice = poison
         vsno = poison
         aice0 = poison
         do it = 1, ldt
            trcr(it) = poison
         end do
         ierr = ior(ierr, ICE_ERR_NONFINITE)
         return
      end if

      call itd_aggregate_area(ncat, aicen, aice, aice0)
      vice = 0.0_wp
      vsno = 0.0_wp
      do n = 1, ncat
         vice = vice + vicen(n)
         vsno = vsno + vsnon(n)
      end do

      poison = aice + vice + vsno
      do it = 1, ntr
         d = ice_trcr_depend(it, nilyr, nslyr)
         atr = 0.0_wp
         do n = 1, ncat
            if (d == 0) then
               w = aicen(n)
            else if (d == 1) then
               w = vicen(n)
            else if (d == 2) then
               w = vsnon(n)
            else
               w = ieee_value(1.0_wp, ieee_quiet_nan)   ! unreachable for it <= ice_ntrcr
            end if
            atr = atr + trcrn(it, n)*w
         end do
         empty = 0.0_wp
         if (it == 1) empty = Tf
         if (.not. ieee_is_finite(poison)) then
            trcr(it) = poison
         else if (.not. ieee_is_finite(atr)) then
            trcr(it) = atr
         else
            trcr(it) = ice_intensive(atr, d, aice, vice, vsno, puny, empty)
         end if
         if (.not. ieee_is_finite(trcr(it))) ierr = ior(ierr, ICE_ERR_NONFINITE)
      end do
      if (.not. ieee_is_finite(aice0)) ierr = ior(ierr, ICE_ERR_NONFINITE)
   end subroutine itd_aggregate

end module rdb_ice_itd_shift
