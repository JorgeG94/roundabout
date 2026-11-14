!! Lateral (floe-edge) melt of sea ice: the melted fraction `rside` per step
!! and its application to every thickness category with the matching
!! ocean fluxes.
module rdb_ice_lateral_melt
   ! SPDX-License-Identifier: BSD-3-Clause
   !
   ! Ported from CICE/Icepack (CICE Consortium, https://github.com/CICE-Consortium).
   ! Upstream: Icepack 9a1f3e4 (Icepack 1.5.3+14), module icepack_therm_itd
   ! (the lateral-melt application) and the `rside` formula that sits in
   ! icepack_therm_vertical (frzmlt_bottom_lateral).
   ! Translated to roundabout's do-concurrent/GPU form (pure per-cell routines,
   ! integer error flags, no optional arguments or allocatables) and modified by
   ! the roundabout developers; see NOTICE for the full notice.
   !
   ! Copyright 1998-2026, Triad National Security, LLC. All rights reserved.
   ! Redistribution and use in source and binary forms, with or without
   ! modification, are permitted under the BSD-3-Clause conditions reproduced in
   ! NOTICE. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS
   ! AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES ARE
   ! DISCLAIMED (full disclaimer in NOTICE).
   !
   ! Original authors (from the upstream headers): William H. Lipscomb,
   ! Elizabeth C. Hunke (LANL), C. M. Bitz (UW).
   ! Physics: Maykut and Perovich (1987, JGR 92 7032); Steele (1992, JGR 97
   ! 17729).
   !
   !! **Physics.**  Floes melt at their edges at the speed
   !! `wlat = m1 * dT**m2` (`dT = max(Tw - Tf, 0)`, `m1 = 1.6e-6 m/s/K**m2`,
   !! `m2 = 1.36`; Maykut and Perovich 1987).  For floes of effective diameter
   !! `D` and shape factor `alpha` the fraction of ice area lost in one step
   !! is `rside = pi * wlat * dt / (alpha * D)` (Steele 1992), clipped to
   !! [0, 1].  `rside` multiplies the area and the ice and snow volumes of
   !! every category alike (no floe-size distribution); the intensive tracers
   !! (surface temperature, enthalpies, salinity) are unchanged, i.e. the
   !! energy per unit volume is constant while the energy per unit area falls.
   !!
   !! **Split with the basal-heat PR (decision M-D4-a).**  This module owns the
   !! pure `rside` function and the application.  It does NOT own the gates
   !! (`aice > puny`, melt potential `frzmlt < 0`), the basal flux `fbot`, or the
   !! limiter that scales `fbot` and `rside` together so the melt energy never
   !! exceeds the ocean's melt potential: those stay with the vertical-thermo
   !! flux-balance routine, which calls `ice_lateral_rside` and then applies
   !! its limiter factor to the result before calling the application here.
   !! Hand-off to that routine: in gated-off cells (`aice <= puny`, `frzmlt >= 0`)
   !! it must write `rside = 0` AND `wlat = 0` (upstream returns before setting
   !! `wlat` when `rside` clips to 0); it must compute the lateral heat flux
   !! `fside` from the state at the TOP of the vertical-thermo step, before any
   !! melt of this step.
   !!
   !! **Flux contract (per unit CELL area, per second, every producer ADDS).**
   !! With `v`, `vs` the ice and snow volumes before the melt, summed over the
   !! categories:
   !! `fresh += (rhoi*v + rhos*vs) * rside / dt` (> 0, into the ocean),
   !! `fsalt += rhoi * v * ice_ref_salinity * 1e-3 * rside / dt` (> 0),
   !! `fhocn += [sum_k qice_k*v/nilyr + sum_k qsno_k*vs/nslyr] * rside / dt`
   !! (< 0: the melt energy is taken from the ocean), `meltl += v * rside` [m].
   !! `aice0` is not updated here; the caller re-aggregates after the whole
   !! thermodynamic tail.  `rside == 0` is a bitwise no-op.
   !!
   !! **Dropped.**  Floe-size-distribution lateral melt, ponds, aerosols,
   !! isotopes, biogeochemistry and the prognostic-salinity salt flux
   !! (`saltflux_option = "prognostic"`): none is carried.
   !!
   !! **Errors (`ierr` bits, OR-ed in, never an abort).**  `ICE_ERR_NONFINITE`:
   !! a non-finite input (the cell is left untouched); `ICE_ERR_AICE_RANGE`:
   !! `rside` outside [0, 1]; `ICE_ERR_AICEN_NEG`: a category area below
   !! `-puny`.  In the last two the cell is also left untouched.
   use, intrinsic :: iso_fortran_env, only: int32, int64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
   use rdb_constants, only: wp
   use rdb_ice_errors, only: ICE_ERR_NONFINITE, ICE_ERR_AICE_RANGE, ICE_ERR_AICEN_NEG
   use rdb_ice_params, only: ice_params_t
   implicit none
   private

   public :: ice_lateral_rside, ice_lateral_melt_cell, ice_lateral_melt_step
   public :: LATERAL_MELT_M1, LATERAL_MELT_M2

   real(wp), parameter :: LATERAL_MELT_M1 = 1.6e-6_wp
      !! Melt-rate coefficient of Maykut and Perovich (1987) [m/s/K**m2].
   real(wp), parameter :: LATERAL_MELT_M2 = 1.36_wp
      !! Melt-rate exponent of Maykut and Perovich (1987).

#include "rdb_safe_math_consts.inc"

contains

#include "rdb_safe_math_core.inc"
#include "rdb_safe_math_pow.inc"

   pure subroutine ice_lateral_rside(dt, sst, Tf, floediam, floeshape, pi, rside, wlat, ierr)
      !! Fraction of the ice melted laterally in one step and the melt speed.
      !! `dT = sst - Tf`; with `dT <= 0` both outputs are exactly 0 and no
      !! power is evaluated.  The real-exponent power is taken once per call
      !! (once per cell, never per category).  A non-finite `dt`, `sst` or
      !! `Tf`, or a NaN quotient, sets `ICE_ERR_NONFINITE` and returns
      !! `rside = wlat = 0`.  Upstream behaviour kept for `floediam`/`floeshape`
      !! at or below zero: only 0/0-type inputs give NaN; a zero product gives
      !! +Inf, clipped to `rside = 1`, and a negative one clips to 0.  The
      !! configuration must refuse `floediam <= 0` and `floeshape <= 0`.
      !$acc routine seq
      real(wp), intent(in) :: dt, sst, Tf, floediam, floeshape, pi
      real(wp), intent(out) :: rside, wlat
      integer(int32), intent(inout) :: ierr
      real(wp) :: dtemp

      rside = 0.0_wp
      wlat = 0.0_wp
      if (.not. (ieee_is_finite(dt) .and. ieee_is_finite(sst) .and. ieee_is_finite(Tf))) then
         ierr = ior(ierr, ICE_ERR_NONFINITE)
         return
      end if
      dtemp = sst - Tf
      if (dtemp <= 0.0_wp) return
      wlat = LATERAL_MELT_M1*rdb_sm_pow(dtemp, LATERAL_MELT_M2)
      rside = wlat*dt*pi/(floeshape*floediam)
      if (ieee_is_nan(rside)) then
         rside = 0.0_wp
         wlat = 0.0_wp
         ierr = ior(ierr, ICE_ERR_NONFINITE)
         return
      end if
      rside = max(0.0_wp, min(rside, 1.0_wp))
   end subroutine ice_lateral_rside

   pure subroutine ice_lateral_melt_cell(i, j, nx, ny, ncat, nilyr, nslyr, rside, dt, &
                                         rhoi, rhos, ice_ref_salinity, puny, &
                                         aicen, vicen, vsnon, qicen, qsnon, &
                                         fresh, fsalt, fhocn, meltl, ierr)
      !! Melt fraction `rside` of every category of cell `(i, j)` and ADD the
      !! resulting fluxes to the cell accumulators (see the module header for
      !! signs).  Area and the ice and snow volumes are scaled by `1 - rside`;
      !! the enthalpy arrays are read only (intensive, unchanged).  Category
      !! `n` is visited in order and the fluxes are added in upstream's order.
      !$acc routine seq
      integer, intent(in) :: i, j, nx, ny, ncat, nilyr, nslyr
      real(wp), intent(in) :: rside, dt, rhoi, rhos, ice_ref_salinity, puny
      real(wp), intent(inout) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat)
      real(wp), intent(in) :: qicen(nx, ny, ncat, nilyr), qsnon(nx, ny, ncat, nslyr)
      real(wp), intent(inout) :: fresh(nx, ny), fsalt(nx, ny), fhocn(nx, ny), meltl(nx, ny)
      integer(int32), intent(inout) :: ierr(nx, ny)
      integer :: n, k
      logical :: bad
      real(wp) :: fr, fs, fh, ml, vi, vs, f1

      if (rside == 0.0_wp) return          ! bitwise no-op
      if (.not. ieee_is_finite(rside)) then
         ierr(i, j) = ior(ierr(i, j), ICE_ERR_NONFINITE)
         return
      end if
      ! ICE_ERR_AICE_RANGE is reused for an out-of-range rside (only a caller bug can trigger it).
      if (rside < 0.0_wp .or. rside > 1.0_wp) then
         ierr(i, j) = ior(ierr(i, j), ICE_ERR_AICE_RANGE)
         return
      end if
      ! Entry guard: every input is checked before anything is written.
      bad = .false.
      do n = 1, ncat
         if (.not. (ieee_is_finite(aicen(i, j, n)) .and. ieee_is_finite(vicen(i, j, n)) &
                    .and. ieee_is_finite(vsnon(i, j, n)))) bad = .true.
         do k = 1, nilyr
            if (.not. ieee_is_finite(qicen(i, j, n, k))) bad = .true.
         end do
         do k = 1, nslyr
            if (.not. ieee_is_finite(qsnon(i, j, n, k))) bad = .true.
         end do
      end do
      if (bad) then
         ierr(i, j) = ior(ierr(i, j), ICE_ERR_NONFINITE)
         return
      end if
      do n = 1, ncat
         if (aicen(i, j, n) < -puny) bad = .true.
      end do
      if (bad) then
         ierr(i, j) = ior(ierr(i, j), ICE_ERR_AICEN_NEG)
         return
      end if

      fr = fresh(i, j)
      fs = fsalt(i, j)
      fh = fhocn(i, j)
      ml = meltl(i, j)
      f1 = 1.0_wp - rside
      do n = 1, ncat
         vi = vicen(i, j, n)               ! volumes before the melt
         vs = vsnon(i, j, n)
         fr = fr + (rhoi*vi + rhos*vs)*rside/dt
         fs = fs + rhoi*vi*ice_ref_salinity*0.001_wp*rside/dt
         ml = ml + vi*rside
         aicen(i, j, n) = aicen(i, j, n)*f1
         vicen(i, j, n) = vi*f1
         vsnon(i, j, n) = vs*f1
         do k = 1, nilyr
            fh = fh + qicen(i, j, n, k)*rside/dt*vi/real(nilyr, wp)
         end do
         do k = 1, nslyr
            fh = fh + qsnon(i, j, n, k)*rside/dt*vs/real(nslyr, wp)
         end do
      end do
      fresh(i, j) = fr
      fsalt(i, j) = fs
      fhocn(i, j) = fh
      meltl(i, j) = ml
   end subroutine ice_lateral_melt_cell

   subroutine ice_lateral_melt_step(nx, ny, ncat, nilyr, nslyr, i_lo, i_hi, j_lo, j_hi, dt, p, rside, &
                                    aicen, vicen, vsnon, qicen, qsnon, fresh, fsalt, fhocn, meltl, ierr)
      !! Column driver: lateral melt on the owned cells `i_lo:i_hi, j_lo:j_hi`
      !! with the per-cell fraction `rside(nx, ny)`.  Every array must be
      !! device-present on a GPU build; the parameter bundle is read on the
      !! host (scalars are passed into the kernel), so it needs no mapping.
      integer, intent(in) :: nx, ny, ncat, nilyr, nslyr, i_lo, i_hi, j_lo, j_hi
      real(wp), intent(in) :: dt
      type(ice_params_t), intent(in) :: p
      real(wp), intent(in) :: rside(nx, ny)
      real(wp), intent(inout) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat)
      real(wp), intent(in) :: qicen(nx, ny, ncat, nilyr), qsnon(nx, ny, ncat, nslyr)
      real(wp), intent(inout) :: fresh(nx, ny), fsalt(nx, ny), fhocn(nx, ny), meltl(nx, ny)
      integer(int32), intent(inout) :: ierr(nx, ny)

      call lateral_melt_impl(nx, ny, ncat, nilyr, nslyr, i_lo, i_hi, j_lo, j_hi, dt, &
                             p%c%rhoi, p%c%rhos, p%thermo%ice_ref_salinity, p%c%puny, rside, &
                             aicen, vicen, vsnon, qicen, qsnon, fresh, fsalt, fhocn, meltl, ierr)
   end subroutine ice_lateral_melt_step

   subroutine lateral_melt_impl(nx, ny, ncat, nilyr, nslyr, i_lo, i_hi, j_lo, j_hi, dt, &
                                rhoi, rhos, sref, puny, rside, &
                                aicen, vicen, vsnon, qicen, qsnon, fresh, fsalt, fhocn, meltl, ierr)
      integer, intent(in) :: nx, ny, ncat, nilyr, nslyr, i_lo, i_hi, j_lo, j_hi
      real(wp), intent(in) :: dt, rhoi, rhos, sref, puny
      real(wp), intent(in) :: rside(nx, ny)
      real(wp), intent(inout) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat)
      real(wp), intent(in) :: qicen(nx, ny, ncat, nilyr), qsnon(nx, ny, ncat, nslyr)
      real(wp), intent(inout) :: fresh(nx, ny), fsalt(nx, ny), fhocn(nx, ny), meltl(nx, ny)
      integer(int32), intent(inout) :: ierr(nx, ny)
      integer :: i, j

      do concurrent(j=j_lo:j_hi, i=i_lo:i_hi)
         call ice_lateral_melt_cell(i, j, nx, ny, ncat, nilyr, nslyr, rside(i, j), dt, &
                                    rhoi, rhos, sref, puny, aicen, vicen, vsnon, qicen, qsnon, &
                                    fresh, fsalt, fhocn, meltl, ierr)
      end do
   end subroutine lateral_melt_impl

end module rdb_ice_lateral_melt
