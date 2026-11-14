!! Thickness-distribution column kernels: line fit, cell aggregates, category shift and rebin.
module rdb_ice_itd_shift
   ! SPDX-License-Identifier: BSD-3-Clause
   !
   ! Ported from CICE/Icepack (CICE Consortium, https://github.com/CICE-Consortium).
   ! Upstream: Icepack 9a1f3e4 (Icepack 1.5.3+14), modules icepack_therm_itd
   ! (line fit) and icepack_itd (aggregate, aggregate_area, shift_ice, rebin).
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
   !!
   !! **In-place form.**  `itd_aggregate_cell` is the same aggregate for cell
   !! `(i, j)` of the full SoA state arrays (`aicen(nx,ny,ncat)`, `Tsfcn`,
   !! `qicen(nx,ny,ncat,nilyr)`, `sicen`, `qsnon(nx,ny,ncat,nslyr)`), reading
   !! tracers through `ice_trcr_get` so no `(ntrcr, ncat)` block is gathered
   !! (the packed form costs ~1.5 KB of per-thread stack at `ICE_NTRCR_MAX`).
   !! It runs the same additions in the same order, so its results are
   !! BITWISE equal to the packed form.  The cell-mean tracers are written to
   !! `(nx,ny[,nlyr])` arrays through `ice_trcr_set_mean`.  Cell-coupled passes (ITD shift,
   !! ridging) should use this form.
   !!
   !! **Category shift and rebin.**  `itd_shift_ice` moves the area `daice(n)`
   !! and ice volume `dvice(n)` across category boundary `n` (`n = 1 .. ncat-1`,
   !! donor `donor(n)` is `n` or `n+1`, 0 for none) together with the snow
   !! volume in proportion to the AREA moved (`dvsnow = vsnon_entry(nd)*daice/
   !! aicen_entry(nd)`) and every tracer: tracer `it` carries the weight
   !! `trcr_base` (area, ice volume or snow volume, one-hot per
   !! `ice_trcr_depend`), the DONOR's entry value moves, and the new intensive
   !! value is recovered by division (`ice_intensive`), so the extensive totals
   !! are conserved to round-off.  Every category's tracers are recomputed by
   !! that division, empty or not (upstream does; a tracer of an untouched
   !! category therefore differs from its entry value by round-off).  Increments
   !! slightly outside `[0, aicen]` (within `puny`) are trimmed; beyond that
   !! the shift is refused with `ICE_ERR_ITD_SHIFT_RANGE` and nothing is
   !! modified.  `snow_offset` is added to the snow enthalpy before and
   !! subtracted after (the remapping passes `rhos*Lfresh` so that the carried
   !! quantity is non-negative; the rebin passes 0).
   !!
   !! `itd_rebin` shifts whole categories whose mean thickness left their bin:
   !! an upward sweep over `n = 1 .. ncat-1`, then a downward one.  It is the
   !! fallback of the remapping and is bitwise the upstream algorithm
   !! (including the `hin_max(0)` floor of category 1, inert for `kitd = 1`).
   !!
   !! Each routine comes in a PACKED form (`aicen(ncat)`, `trcrn(ldt, ncat)`,
   !! a host-side and test reference) and a CELL form that works on the SoA
   !! state through the tracer accessors with only `O(ncat)` per-thread scratch
   !! (`itd_shift_ice_cell`, `itd_rebin_cell`).  The two run the same additions
   !! in the same order and agree bitwise on CPU builds.
   !!
   !! Upstream branches dropped as inert for BL99 (`ktherm = 1`), no ponds, no
   !! floe-size distribution, no bio: the `n_trcr_strata` multiplications (OM3
   !! has no strata tracer), the brine-height, snow-density, snow-grain and
   !! aerosol/isotope tracer branches of `icepack_compute_tracers`, and the
   !! `warnings` text.
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_constants, only: wp
   use rdb_ice_dims, only: ICE_NCAT_MAX
   use, intrinsic :: iso_fortran_env, only: int32
   use rdb_ice_errors, only: ICE_ERR_NONFINITE, ICE_ERR_ITD_SHIFT_RANGE
   use rdb_ice_params, only: ice_trcr_depend, ice_intensive, ice_ntrcr, ice_trcr_get, ice_trcr_set, ice_trcr_set_mean, ICE_QNAN
   implicit none
   private

   public :: itd_fit_line, itd_aggregate_area, itd_aggregate, itd_aggregate_cell
   public :: itd_shift_amounts, itd_shift_tracers, itd_shift_tracers_cell, itd_shift_ice, itd_shift_ice_cell
   public :: itd_rebin, itd_rebin_cell

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
      if (ntr > ldt .or. ntr < 1 .or. ncat < 1 .or. nilyr < 1 .or. nslyr < 1) then
         poison = ICE_QNAN
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
               w = ICE_QNAN   ! unreachable for it <= ice_ntrcr
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

   pure subroutine itd_aggregate_cell(i, j, nx, ny, ncat, nilyr, nslyr, aicen, vicen, vsnon, &
                                      Tsfcn, qicen, sicen, qsnon, Tf, puny, &
                                      aice, vice, vsno, aice0, Tsfc_c, qice_c, sice_c, qsno_c, ierr)
      !! In-place twin of `itd_aggregate` for cell `(i, j)` of full SoA state
      !! arrays (see the module note).  Scalars `aice`, `vice`, `vsno`, `aice0`
      !! are returned; the intensive cell-mean tracers are written to
      !! `Tsfc_c(i,j)`, `qice_c(i,j,1:nilyr)`, `sice_c(i,j,1:nilyr)` and
      !! `qsno_c(i,j,1:nslyr)` (the packed order `Tsfc | qice | qsno | sice`).
      !! Same error semantics as the packed form: any non-finite input or output
      !! ORs `ICE_ERR_NONFINITE` into `ierr(i,j)`.  Invalid counts (`ncat`, `nilyr` or
      !! `nslyr` < 1) poison only the four scalars `aice, vice, vsno, aice0` and set the
      !! error bit: the tracer arrays are NOT written, because with `nilyr < 1` or
      !! `nslyr < 1` the packed tracer indices (and the array extents) do not exist.
      !! Outputs are written directly (no
      !! double buffering): after `ierr /= 0` they are unspecified.
      !$acc routine seq
      integer, value :: i, j
      integer, intent(in) :: nx, ny, ncat, nilyr, nslyr
      real(wp), intent(in) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat)
      real(wp), intent(in) :: Tsfcn(nx, ny, ncat), qicen(nx, ny, ncat, nilyr)
      real(wp), intent(in) :: sicen(nx, ny, ncat, nilyr), qsnon(nx, ny, ncat, nslyr)
      real(wp), intent(in) :: Tf, puny
      real(wp), intent(out) :: aice, vice, vsno, aice0
      real(wp), intent(inout) :: Tsfc_c(nx, ny), qice_c(nx, ny, nilyr)
      real(wp), intent(inout) :: sice_c(nx, ny, nilyr), qsno_c(nx, ny, nslyr)
      integer(int32), intent(inout) :: ierr(nx, ny)

      integer :: n, it, d, ntr
      real(wp) :: atr, w, empty, poison, tv

      ntr = ice_ntrcr(nilyr, nslyr)
      if (ntr < 1 .or. ncat < 1 .or. nilyr < 1 .or. nslyr < 1) then
         poison = ICE_QNAN
         aice = poison
         vice = poison
         vsno = poison
         aice0 = poison
         ierr(i, j) = ior(ierr(i, j), ICE_ERR_NONFINITE)
         return
      end if

      aice = 0.0_wp
      vice = 0.0_wp
      vsno = 0.0_wp
      do n = 1, ncat
         aice = aice + aicen(i, j, n)
      end do
      if (ieee_is_finite(aice)) then
         aice0 = max(1.0_wp - aice, 0.0_wp)
      else
         aice0 = aice
      end if
      do n = 1, ncat
         vice = vice + vicen(i, j, n)
         vsno = vsno + vsnon(i, j, n)
      end do

      poison = aice + vice + vsno
      do it = 1, ntr
         d = ice_trcr_depend(it, nilyr, nslyr)
         atr = 0.0_wp
         do n = 1, ncat
            if (d == 0) then
               w = aicen(i, j, n)
            else if (d == 1) then
               w = vicen(i, j, n)
            else if (d == 2) then
               w = vsnon(i, j, n)
            else
               w = ICE_QNAN   ! unreachable for it <= ice_ntrcr
            end if
            atr = atr + ice_trcr_get(it, i, j, n, nx, ny, ncat, nilyr, nslyr, Tsfcn, qicen, sicen, qsnon)*w
         end do
         empty = 0.0_wp
         if (it == 1) empty = Tf
         if (.not. ieee_is_finite(poison)) then
            tv = poison
         else if (.not. ieee_is_finite(atr)) then
            tv = atr
         else
            tv = ice_intensive(atr, d, aice, vice, vsno, puny, empty)
         end if
         if (.not. ieee_is_finite(tv)) ierr(i, j) = ior(ierr(i, j), ICE_ERR_NONFINITE)
         call ice_trcr_set_mean(it, i, j, nx, ny, nilyr, nslyr, tv, Tsfc_c, qice_c, qsno_c, sice_c)
      end do
      if (.not. ieee_is_finite(aice0)) ierr(i, j) = ior(ierr(i, j), ICE_ERR_NONFINITE)
   end subroutine itd_aggregate_cell

   pure subroutine itd_shift_amounts(ncat, puny, donor, daice, dvice, aicen, vicen, vsnon, hicen, &
                                     a0, v0, s0, dvsn, ok)
      !! Area, ice-volume and snow-volume part of `itd_shift_ice` on packed
      !! `ncat`-vectors.  Trims or refuses the increments (see the module note),
      !! then moves them; returns the entry values `a0, v0, s0` and the snow
      !! volume moved across each boundary `dvsn`, which the tracer pass needs,
      !! and the refreshed thickness `hicen`.  `ok = .false.` (increment outside
      !! its tolerance band): `aicen, vicen, vsnon` are untouched, `daice`, `dvice`
      !! may have been trimmed, and the other outputs are unspecified.
      !$acc routine seq
      integer, intent(in) :: ncat
      real(wp), intent(in) :: puny
      integer, intent(in) :: donor(ncat)
      real(wp), intent(inout) :: daice(ncat), dvice(ncat)
      real(wp), intent(inout) :: aicen(ncat), vicen(ncat), vsnon(ncat)
      real(wp), intent(out) :: hicen(ncat), a0(ncat), v0(ncat), s0(ncat), dvsn(ncat)
      logical, intent(out) :: ok

      integer :: n, nd, nr
      real(wp) :: worka

      ok = .true.
      do n = 1, ncat
         a0(n) = aicen(n)
         v0(n) = vicen(n)
         s0(n) = vsnon(n)
         dvsn(n) = 0.0_wp
         hicen(n) = 0.0_wp
      end do

      ! daice or dvice out of range, allowing for round-off
      do n = 1, ncat - 1
         if (donor(n) > 0) then
            nd = donor(n)
            if (daice(n) < 0.0_wp) then
               if (daice(n) > -puny*aicen(nd)) then
                  daice(n) = 0.0_wp   ! shift no ice
                  dvice(n) = 0.0_wp
               else
                  ok = .false.
               end if
            end if
            if (dvice(n) < 0.0_wp) then
               if (dvice(n) > -puny*vicen(nd)) then
                  daice(n) = 0.0_wp
                  dvice(n) = 0.0_wp
               else
                  ok = .false.
               end if
            end if
            if (daice(n) > aicen(nd)*(1.0_wp - puny)) then
               if (daice(n) < aicen(nd)*(1.0_wp + puny)) then
                  daice(n) = aicen(nd)
                  dvice(n) = vicen(nd)
               else
                  ok = .false.
               end if
            end if
            if (dvice(n) > vicen(nd)*(1.0_wp - puny)) then
               if (dvice(n) < vicen(nd)*(1.0_wp + puny)) then
                  daice(n) = aicen(nd)
                  dvice(n) = vicen(nd)
               else
                  ok = .false.
               end if
            end if
         end if
         if (.not. ok) return
      end do

      ! transfer area and volumes between categories
      do n = 1, ncat - 1
         if (daice(n) > 0.0_wp .and. donor(n) > 0) then
            nd = donor(n)
            if (nd == n) then
               nr = nd + 1
            else
               nr = n
            end if
            aicen(nd) = aicen(nd) - daice(n)
            aicen(nr) = aicen(nr) + daice(n)
            vicen(nd) = vicen(nd) - dvice(n)
            vicen(nr) = vicen(nr) + dvice(n)
            worka = daice(n)/a0(nd)
            dvsn(n) = s0(nd)*worka        ! snow moves in proportion to the area moved
            vsnon(nd) = vsnon(nd) - dvsn(n)
            vsnon(nr) = vsnon(nr) + dvsn(n)
         end if
      end do

      do n = 1, ncat
         if (aicen(n) > puny) then
            hicen(n) = vicen(n)/aicen(n)
         else
            hicen(n) = 0.0_wp
         end if
      end do
   end subroutine itd_shift_amounts

   pure subroutine itd_shift_tracers(ncat, ldt, nilyr, nslyr, puny, Tf, snow_offset, donor, daice, dvice, dvsn, &
                                     a0, v0, s0, aicen, vicen, vsnon, trcrn)
      !! Tracer part of `itd_shift_ice` on a packed column `trcrn(ldt, ncat)`
      !! (entries `ice_ntrcr(nilyr, nslyr)+1:ldt` are not touched).  `daice`,
      !! `dvice`, `dvsn`, `a0, v0, s0` are the outputs of `itd_shift_amounts`;
      !! `aicen, vicen, vsnon` the categories AFTER the transfer.
      !$acc routine seq
      integer, intent(in) :: ncat, ldt, nilyr, nslyr
      real(wp), intent(in) :: puny, Tf, snow_offset
      integer, intent(in) :: donor(ncat)
      real(wp), intent(in) :: daice(ncat), dvice(ncat), dvsn(ncat), a0(ncat), v0(ncat), s0(ncat)
      real(wp), intent(in) :: aicen(ncat), vicen(ncat), vsnon(ncat)
      real(wp), intent(inout) :: trcrn(ldt, ncat)

      integer :: it, n, nd, nr, d, ntr
      real(wp) :: off, empty, w, amt, datr
      real(wp) :: atr(ICE_NCAT_MAX), t0(ICE_NCAT_MAX)

      ntr = ice_ntrcr(nilyr, nslyr)
      do it = 1, ntr
         d = ice_trcr_depend(it, nilyr, nslyr)
         off = 0.0_wp
         if (d == 2) off = snow_offset
         empty = 0.0_wp
         if (it == 1) empty = Tf
         do n = 1, ncat
            t0(n) = trcrn(it, n) + off
            if (d == 0) then
               w = a0(n)
            else if (d == 1) then
               w = v0(n)
            else
               w = s0(n)
            end if
            atr(n) = t0(n)*w
         end do
         do n = 1, ncat - 1
            if (daice(n) > 0.0_wp .and. donor(n) > 0) then
               nd = donor(n)
               if (nd == n) then
                  nr = nd + 1
               else
                  nr = n
               end if
               if (d == 0) then
                  amt = daice(n)
               else if (d == 1) then
                  amt = dvice(n)
               else
                  amt = dvsn(n)
               end if
               datr = t0(nd)*amt
               atr(nd) = atr(nd) - datr
               atr(nr) = atr(nr) + datr
            end if
         end do
         do n = 1, ncat
            trcrn(it, n) = ice_intensive(atr(n), d, aicen(n), vicen(n), vsnon(n), puny, empty)
            if (d == 2) trcrn(it, n) = trcrn(it, n) - off
         end do
      end do
   end subroutine itd_shift_tracers

   pure subroutine itd_shift_tracers_cell(i, j, nx, ny, ncat, nilyr, nslyr, puny, Tf, snow_offset, donor, daice, dvice, &
                                          dvsn, a0, v0, s0, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon)
      !! In-place twin of `itd_shift_tracers` for cell `(i, j)` of the SoA tracer
      !! arrays (read and written through `ice_trcr_get/set`; per-thread scratch
      !! is `2*ICE_NCAT_MAX` reals).  `aicen, vicen, vsnon` are the packed
      !! post-transfer category vectors of this cell.  Bitwise equal to the
      !! packed form on CPU builds.
      !$acc routine seq
      integer, value :: i, j
      integer, intent(in) :: nx, ny, ncat, nilyr, nslyr
      real(wp), intent(in) :: puny, Tf, snow_offset
      integer, intent(in) :: donor(ncat)
      real(wp), intent(in) :: daice(ncat), dvice(ncat), dvsn(ncat), a0(ncat), v0(ncat), s0(ncat)
      real(wp), intent(in) :: aicen(ncat), vicen(ncat), vsnon(ncat)
      real(wp), intent(inout) :: Tsfcn(nx, ny, ncat), qicen(nx, ny, ncat, nilyr)
      real(wp), intent(inout) :: sicen(nx, ny, ncat, nilyr), qsnon(nx, ny, ncat, nslyr)

      integer :: it, n, nd, nr, d, ntr
      real(wp) :: off, empty, w, amt, datr, tv
      real(wp) :: atr(ICE_NCAT_MAX), t0(ICE_NCAT_MAX)

      ntr = ice_ntrcr(nilyr, nslyr)
      do it = 1, ntr
         d = ice_trcr_depend(it, nilyr, nslyr)
         off = 0.0_wp
         if (d == 2) off = snow_offset
         empty = 0.0_wp
         if (it == 1) empty = Tf
         do n = 1, ncat
            t0(n) = ice_trcr_get(it, i, j, n, nx, ny, ncat, nilyr, nslyr, Tsfcn, qicen, sicen, qsnon) + off
            if (d == 0) then
               w = a0(n)
            else if (d == 1) then
               w = v0(n)
            else
               w = s0(n)
            end if
            atr(n) = t0(n)*w
         end do
         do n = 1, ncat - 1
            if (daice(n) > 0.0_wp .and. donor(n) > 0) then
               nd = donor(n)
               if (nd == n) then
                  nr = nd + 1
               else
                  nr = n
               end if
               if (d == 0) then
                  amt = daice(n)
               else if (d == 1) then
                  amt = dvice(n)
               else
                  amt = dvsn(n)
               end if
               datr = t0(nd)*amt
               atr(nd) = atr(nd) - datr
               atr(nr) = atr(nr) + datr
            end if
         end do
         do n = 1, ncat
            tv = ice_intensive(atr(n), d, aicen(n), vicen(n), vsnon(n), puny, empty)
            if (d == 2) tv = tv - off
            call ice_trcr_set(it, i, j, n, nx, ny, ncat, nilyr, nslyr, tv, Tsfcn, qicen, sicen, qsnon)
         end do
      end do
   end subroutine itd_shift_tracers_cell

   pure subroutine itd_shift_ice(ncat, ldt, nilyr, nslyr, puny, Tf, snow_offset, donor, daice, dvice, &
                                 aicen, vicen, vsnon, hicen, trcrn, ierr)
      !! Packed `shift_ice`: `itd_shift_amounts` then `itd_shift_tracers`.  A
      !! refused increment ORs `ICE_ERR_ITD_SHIFT_RANGE` into `ierr` and leaves
      !! the state untouched.  `hicen` is the refreshed thickness on success.
      !$acc routine seq
      integer, intent(in) :: ncat, ldt, nilyr, nslyr
      real(wp), intent(in) :: puny, Tf, snow_offset
      integer, intent(in) :: donor(ncat)
      real(wp), intent(inout) :: daice(ncat), dvice(ncat)
      real(wp), intent(inout) :: aicen(ncat), vicen(ncat), vsnon(ncat)
      real(wp), intent(out) :: hicen(ncat)
      real(wp), intent(inout) :: trcrn(ldt, ncat)
      integer(int32), intent(inout) :: ierr

      real(wp) :: a0(ICE_NCAT_MAX), v0(ICE_NCAT_MAX), s0(ICE_NCAT_MAX), dvsn(ICE_NCAT_MAX)
      logical :: ok

      call itd_shift_amounts(ncat, puny, donor, daice, dvice, aicen, vicen, vsnon, hicen, a0, v0, s0, dvsn, ok)
      if (.not. ok) then
         ierr = ior(ierr, ICE_ERR_ITD_SHIFT_RANGE)
         return
      end if
      call itd_shift_tracers(ncat, ldt, nilyr, nslyr, puny, Tf, snow_offset, donor, daice, dvice, dvsn, &
                             a0, v0, s0, aicen, vicen, vsnon, trcrn)
   end subroutine itd_shift_ice

   pure subroutine itd_shift_ice_cell(i, j, nx, ny, ncat, nilyr, nslyr, puny, Tf, snow_offset, donor, daice, dvice, &
                                      aicen, vicen, vsnon, hicen, Tsfcn, qicen, sicen, qsnon, ierr)
      !! Cell form of `itd_shift_ice`: area, ice and snow volume of the cell
      !! are the packed vectors `aicen, vicen, vsnon` (the caller's per-cell
      !! gather of three `ncat` vectors), the tracers are the SoA arrays at
      !! `(i, j)`, updated in place.  On a refused increment `ierr` (the cell's word) gets
      !! `ICE_ERR_ITD_SHIFT_RANGE` and nothing is modified.
      !$acc routine seq
      integer, value :: i, j
      integer, intent(in) :: nx, ny, ncat, nilyr, nslyr
      real(wp), intent(in) :: puny, Tf, snow_offset
      integer, intent(in) :: donor(ncat)
      real(wp), intent(inout) :: daice(ncat), dvice(ncat)
      real(wp), intent(inout) :: aicen(ncat), vicen(ncat), vsnon(ncat)
      real(wp), intent(out) :: hicen(ncat)
      real(wp), intent(inout) :: Tsfcn(nx, ny, ncat), qicen(nx, ny, ncat, nilyr)
      real(wp), intent(inout) :: sicen(nx, ny, ncat, nilyr), qsnon(nx, ny, ncat, nslyr)
      integer(int32), intent(inout) :: ierr

      real(wp) :: a0(ICE_NCAT_MAX), v0(ICE_NCAT_MAX), s0(ICE_NCAT_MAX), dvsn(ICE_NCAT_MAX)
      logical :: ok

      call itd_shift_amounts(ncat, puny, donor, daice, dvice, aicen, vicen, vsnon, hicen, a0, v0, s0, dvsn, ok)
      if (.not. ok) then
         ierr = ior(ierr, ICE_ERR_ITD_SHIFT_RANGE)
         return
      end if
      call itd_shift_tracers_cell(i, j, nx, ny, ncat, nilyr, nslyr, puny, Tf, snow_offset, donor, daice, dvice, &
                                  dvsn, a0, v0, s0, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon)
   end subroutine itd_shift_ice_cell

   pure subroutine itd_rebin(ncat, ldt, nilyr, nslyr, hin_max, puny, Tf, aicen, vicen, vsnon, trcrn, ierr)
      !! Packed `rebin`: shift whole categories whose mean thickness left their
      !! bin `(hin_max(n-1), hin_max(n)]`, up then down (see the module note).
      !! Each shift is a full `itd_shift_ice` with `snow_offset = 0`.  A refused
      !! shift (cannot happen for whole-category moves of finite states) ORs
      !! `ICE_ERR_ITD_SHIFT_RANGE` and returns with the earlier shifts applied (the state is
      !! then UNSPECIFIED).  The snow enthalpy of an emptied category is 0 after a rebin (no
      !! `rhos*Lfresh` offset, unlike the remapping).
      !$acc routine seq
      integer, intent(in) :: ncat, ldt, nilyr, nslyr
      real(wp), intent(in) :: puny, Tf
      real(wp), intent(in) :: hin_max(0:ncat)
      real(wp), intent(inout) :: aicen(ncat), vicen(ncat), vsnon(ncat)
      real(wp), intent(inout) :: trcrn(ldt, ncat)
      integer(int32), intent(inout) :: ierr

      integer :: n
      integer :: donor(ICE_NCAT_MAX)
      integer(int32) :: ierr_l
      real(wp) :: daice(ICE_NCAT_MAX), dvice(ICE_NCAT_MAX), hicen(ICE_NCAT_MAX)

      do n = 1, ncat
         donor(n) = 0
         daice(n) = 0.0_wp
         dvice(n) = 0.0_wp
         if (aicen(n) > puny) then
            hicen(n) = vicen(n)/aicen(n)
         else
            hicen(n) = 0.0_wp
         end if
      end do

      ! make sure the thickness of category 1 is at least hin_max(0)
      if (aicen(1) > puny) then
         if (hicen(1) <= hin_max(0) .and. hin_max(0) > 0.0_wp) then
            aicen(1) = vicen(1)/hin_max(0)
            hicen(1) = hin_max(0)
         end if
      end if

      ! move thin categories up
      do n = 1, ncat - 1
         if (aicen(n) > puny .and. hicen(n) > hin_max(n)) then
            donor(n) = n
            daice(n) = aicen(n)
            dvice(n) = vicen(n)
            ierr_l = 0_int32
            call itd_shift_ice(ncat, ldt, nilyr, nslyr, puny, Tf, 0.0_wp, donor, daice, dvice, &
                               aicen, vicen, vsnon, hicen, trcrn, ierr_l)
            if (ierr_l /= 0_int32) then
               ierr = ior(ierr, ierr_l)
               return
            end if
            donor(n) = 0
            daice(n) = 0.0_wp
            dvice(n) = 0.0_wp
         end if
      end do

      ! move thick categories down
      do n = ncat - 1, 1, -1
         if (aicen(n + 1) > puny .and. hicen(n + 1) <= hin_max(n)) then
            donor(n) = n + 1
            daice(n) = aicen(n + 1)
            dvice(n) = vicen(n + 1)
            ierr_l = 0_int32
            call itd_shift_ice(ncat, ldt, nilyr, nslyr, puny, Tf, 0.0_wp, donor, daice, dvice, &
                               aicen, vicen, vsnon, hicen, trcrn, ierr_l)
            if (ierr_l /= 0_int32) then
               ierr = ior(ierr, ierr_l)
               return
            end if
            donor(n) = 0
            daice(n) = 0.0_wp
            dvice(n) = 0.0_wp
         end if
      end do
   end subroutine itd_rebin

   pure subroutine itd_rebin_cell(i, j, nx, ny, ncat, nilyr, nslyr, hin_max, puny, Tf, &
                                  aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, ierr)
      !! Cell form of `itd_rebin` on the SoA state, in place: area and volumes
      !! are gathered into three `ncat` vectors and written back, tracers move
      !! through `ice_trcr_get/set`.  Bitwise equal to the packed form on CPU
      !! builds.  On a refused shift `ierr(i, j)` gets `ICE_ERR_ITD_SHIFT_RANGE`
      !! and the cell state is UNSPECIFIED (the same as the packed form).
      !$acc routine seq
      integer, value :: i, j
      integer, intent(in) :: nx, ny, ncat, nilyr, nslyr
      real(wp), intent(in) :: puny, Tf
      real(wp), intent(in) :: hin_max(0:ncat)
      real(wp), intent(inout) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat)
      real(wp), intent(inout) :: Tsfcn(nx, ny, ncat), qicen(nx, ny, ncat, nilyr)
      real(wp), intent(inout) :: sicen(nx, ny, ncat, nilyr), qsnon(nx, ny, ncat, nslyr)
      integer(int32), intent(inout) :: ierr(nx, ny)

      integer :: n
      integer :: donor(ICE_NCAT_MAX)
      integer(int32) :: ierr_l
      real(wp) :: a(ICE_NCAT_MAX), v(ICE_NCAT_MAX), s(ICE_NCAT_MAX)
      real(wp) :: daice(ICE_NCAT_MAX), dvice(ICE_NCAT_MAX), hicen(ICE_NCAT_MAX)

      do n = 1, ncat
         a(n) = aicen(i, j, n)
         v(n) = vicen(i, j, n)
         s(n) = vsnon(i, j, n)
         donor(n) = 0
         daice(n) = 0.0_wp
         dvice(n) = 0.0_wp
         if (a(n) > puny) then
            hicen(n) = v(n)/a(n)
         else
            hicen(n) = 0.0_wp
         end if
      end do

      if (a(1) > puny) then
         if (hicen(1) <= hin_max(0) .and. hin_max(0) > 0.0_wp) then
            a(1) = v(1)/hin_max(0)
            hicen(1) = hin_max(0)
         end if
      end if

      ierr_l = 0_int32
      do n = 1, ncat - 1
         if (a(n) > puny .and. hicen(n) > hin_max(n)) then
            donor(n) = n
            daice(n) = a(n)
            dvice(n) = v(n)
            call itd_shift_ice_cell(i, j, nx, ny, ncat, nilyr, nslyr, puny, Tf, 0.0_wp, donor, daice, dvice, &
                                    a, v, s, hicen, Tsfcn, qicen, sicen, qsnon, ierr_l)
            if (ierr_l /= 0_int32) then
               ierr(i, j) = ior(ierr(i, j), ierr_l)
               return
            end if
            donor(n) = 0
            daice(n) = 0.0_wp
            dvice(n) = 0.0_wp
         end if
      end do
      do n = ncat - 1, 1, -1
         if (a(n + 1) > puny .and. hicen(n + 1) <= hin_max(n)) then
            donor(n) = n + 1
            daice(n) = a(n + 1)
            dvice(n) = v(n + 1)
            call itd_shift_ice_cell(i, j, nx, ny, ncat, nilyr, nslyr, puny, Tf, 0.0_wp, donor, daice, dvice, &
                                    a, v, s, hicen, Tsfcn, qicen, sicen, qsnon, ierr_l)
            if (ierr_l /= 0_int32) then
               ierr(i, j) = ior(ierr(i, j), ierr_l)
               return
            end if
            donor(n) = 0
            daice(n) = 0.0_wp
            dvice(n) = 0.0_wp
         end if
      end do

      do n = 1, ncat
         aicen(i, j, n) = a(n)
         vicen(i, j, n) = v(n)
         vsnon(i, j, n) = s(n)
      end do
   end subroutine itd_rebin_cell

end module rdb_ice_itd_shift
