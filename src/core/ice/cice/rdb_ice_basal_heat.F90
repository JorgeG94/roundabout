!! Heat flux from the ocean to the ice base (`fbot`) and the lateral-melt
!! fraction (`rside`, `wlat`), limited together by the ocean melt potential.
module rdb_ice_basal_heat
   ! SPDX-License-Identifier: BSD-3-Clause
   !
   ! Ported from CICE/Icepack (CICE Consortium, https://github.com/CICE-Consortium).
   ! Upstream: Icepack 9a1f3e4 (Icepack 1.5.3+14), module icepack_therm_vertical
   ! (frzmlt_bottom_lateral).
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
   ! Original authors (from the upstream headers): C. M. Bitz (UW),
   ! William H. Lipscomb, Elizabeth C. Hunke (LANL).
   ! Physics: Maykut and McPhee (1995, JGR 100 24691); McPhee (1992, J. Geophys. Res.
   ! 97(C4), 5365-5379);
   ! Maykut and Perovich (1987, JGR 92 7032); Steele (1992, JGR 97 17729).
   !
   !! **Physics.**  In a cell with ice (`aice > puny`) whose ocean has a melt
   !! potential (`frzmlt < 0`, W/m2 per CELL area) the ocean heats the ice base
   !! through the turbulent boundary layer under the ice (Maykut and McPhee
   !! 1995): `fbot = -cp_ocn*rhow*Ch*u*(sst - Tf)`, bounded below by `frzmlt`,
   !! with `Ch = 0.006` the constant Stanton number and `u*` the friction
   !! velocity of the ice-ocean stress, `u* = max(sqrt(|tau|/rhow), ustar_min)`
   !! (`tau` per ICE area at the tracer point, see **Stress contract**).  The floes also melt
   !! at their edges, `rside` from `ice_lateral_rside` (that module owns the
   !! formula; this one calls it).  The two sinks are then limited together: the
   !! lateral heat flux is `fside = sum_n rside*(Ei_n + Es_n)/dt` (<= 0, `E` the
   !! enthalpy per unit area of category `n` before the melt) and
   !! `xtmp = clamp(frzmlt/(fbot + fside - puny), 0, 1)` multiplies both
   !! `fbot` and `rside`, so that the melt never takes more energy than the
   !! ocean offers.
   !!
   !! **Signs.**  `fbot <= 0` and `fbot >= frzmlt`: heat taken FROM the ocean
   !! and spent on basal melt, per unit ICE area (the vertical thermodynamics
   !! reads it as `ebot_mlt = (fcondbot - fbot)*dt`).  `frzmlt` and `fside` are
   !! per CELL area while `fbot` is per ice area, so for `aice < 1` the limiter
   !! is conservative (kept verbatim).  `Tbot = Tf` always.
   !!
   !! **Gates and exits.**  With `aice <= puny` or `frzmlt >= 0` the outputs are
   !! `fbot = rside = wlat = 0` (and `xtmp = 1`).  When `rside` evaluates to
   !! exactly 0 the limiter is skipped (upstream returns there), `fbot` keeps
   !! its `frzmlt` bound only, and `wlat = 0` (not the unclipped speed).
   !!
   !! **Heat-transfer option.**  Only `fbot_xfer_type = 'constant'` (OM3) is
   !! ported; `'Cdn_ocn'` needs the form-drag coefficient and is not offered
   !! (no knob exists, so no refusal either).
   !!
   !! **Stress contract (for the glue).**  `tau_x/tau_y` is the ice-ocean stress
   !! PER ICE AREA at the T point, built the way the dynamics does it: the
   !! per-cell-area stress on each face divided by the FACE ice area where that
   !! exceeds `puny` (0 otherwise), halo-updated, and only THEN averaged to T.
   !! It is the PREVIOUS step's field (one-step lag).  Dividing a T-averaged
   !! stress by the T-cell ice area is not equivalent and must not be used.
   !!
   !! **Errors (`ierr` bits, OR-ed in, never an abort).**  `ICE_ERR_NONFINITE`: a
   !! non-finite scalar input, friction velocity (e.g. `tau**2` overflow), melt
   !! energy or lateral flux; `ICE_ERR_AICEN_NEG`: a category area below `-puny`
   !! (cell entry only).  Both bits are fatal.  On an error bit the outputs are
   !! UNSPECIFIED (the implementation happens to zero them) and the driver must
   !! not commit the cell.  Preconditions, owned by the caller (the routine does
   !! not test them): `floediam > 0`, `floeshape > 0`, `ustar_min >= 0`, and the
   !! constants (`rhow`, `cp_ocn`, `pi`, `puny`) taken from a validated constants
   !! set (`p%c`).  `floeshape <= 0` is refused by the configuration on branches
   !! that carry that refusal; `ustar_min` has no knob and comes from
   !! `p%thermo%ustar_min`.
   use, intrinsic :: iso_fortran_env, only: int32
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_constants, only: wp
   use rdb_ice_dims, only: ICE_NCAT_MAX
   use rdb_ice_errors, only: ICE_ERR_NONFINITE, ICE_ERR_AICEN_NEG
   use rdb_ice_lateral_melt, only: ice_lateral_rside
   use rdb_ice_params, only: ice_params_t
   implicit none
   private

   public :: ice_frzmlt_bottom_lateral, ice_frzmlt_bottom_lateral_cell, ice_basal_heat_step
   public :: ICE_FBOT_CH

   real(wp), parameter :: ICE_FBOT_CH = 0.006_wp
      !! Constant Stanton number of the basal heat transfer
      !! (`fbot_xfer_type = 'constant'`; McPhee 1992).

contains

   pure subroutine ice_frzmlt_bottom_lateral(ncat, dt, aice, frzmlt, sst, Tf, tau_x, tau_y, &
                                             ustar_min, floediam, floeshape, rhow, cp_ocn, pi, puny, &
                                             etotn, Tbot, fbot, rside, wlat, xtmp, ierr)
      !! Column form: one cell, categories reduced to the per-category melt
      !! energy `etotn(n) = sum_k qsno*vsno/nslyr + sum_k qice*vice/nilyr`
      !! [J/m2, <= 0], summed snow first then ice in layer order (the caller owns
      !! that order; it is what makes the result match upstream bit for bit).
      !! `etotn` is read only inside the melt gate.  `rside` returned is the
      !! limited fraction applied to every category alike.
      !$acc routine seq
      integer, intent(in) :: ncat
      real(wp), intent(in) :: dt, aice, frzmlt, sst, Tf, tau_x, tau_y
      real(wp), intent(in) :: ustar_min, floediam, floeshape, rhow, cp_ocn, pi, puny
      real(wp), intent(in) :: etotn(ncat)
      real(wp), intent(out) :: Tbot, fbot, rside, wlat, xtmp
      integer(int32), intent(inout) :: ierr
      integer :: n
      real(wp) :: deltaT, ustar, cpchr, fside
      integer(int32) :: ie

      Tbot = Tf
      fbot = 0.0_wp
      rside = 0.0_wp
      wlat = 0.0_wp
      xtmp = 1.0_wp
      if (.not. (ieee_is_finite(dt) .and. ieee_is_finite(aice) .and. ieee_is_finite(frzmlt) .and. &
                 ieee_is_finite(sst) .and. ieee_is_finite(Tf) .and. ieee_is_finite(tau_x) .and. &
                 ieee_is_finite(tau_y))) then
         ierr = ior(ierr, ICE_ERR_NONFINITE)
         return
      end if
      if (.not. (aice > puny .and. frzmlt < 0.0_wp)) return      ! nothing can melt

      deltaT = max(sst - Tbot, 0.0_wp)
      ustar = sqrt(sqrt(tau_x**2 + tau_y**2)/rhow)
      ustar = max(ustar, ustar_min)
      if (.not. ieee_is_finite(ustar)) then      ! tau**2 overflow: 0*Inf would be laundered by max/min
         ierr = ior(ierr, ICE_ERR_NONFINITE)
         fbot = 0.0_wp
         return
      end if
      cpchr = -cp_ocn*rhow*ICE_FBOT_CH
      fbot = cpchr*deltaT*ustar                  ! < 0
      fbot = max(fbot, frzmlt)                   ! frzmlt <= fbot <= 0

      ie = 0_int32
      call ice_lateral_rside(dt, sst, Tf, floediam, floeshape, pi, rside, wlat, ie)
      if (rside == 0.0_wp) then                  ! upstream returns here, wlat stays 0
         wlat = 0.0_wp
         if (ie /= 0_int32) then                 ! non-finite quotient: refuse the cell
            ierr = ior(ierr, ie)
            fbot = 0.0_wp
         end if
         return
      end if

      fside = 0.0_wp
      do n = 1, ncat
         fside = fside + rside*etotn(n)/dt
      end do
      if (.not. ieee_is_finite(fside)) then
         ierr = ior(ierr, ICE_ERR_NONFINITE)
         fbot = 0.0_wp
         rside = 0.0_wp
         wlat = 0.0_wp
         return
      end if
      xtmp = frzmlt/(fbot + fside - puny)
      xtmp = min(xtmp, 1.0_wp)
      xtmp = max(xtmp, 0.0_wp)
      fbot = fbot*xtmp
      rside = rside*xtmp
   end subroutine ice_frzmlt_bottom_lateral

   pure subroutine ice_frzmlt_bottom_lateral_cell(i, j, nx, ny, ncat, nilyr, nslyr, dt, &
                                                  ustar_min, floediam, floeshape, rhow, cp_ocn, pi, puny, &
                                                  aicen, vicen, vsnon, qicen, qsnon, &
                                                  frzmlt, sst, Tf, tau_x, tau_y, fbot, rside, wlat, ierr)
      !! Cell entry: reads the category state of cell `(i, j)` in place and
      !! writes `fbot(i,j)`, `rside(i,j)`, `wlat(i,j)` (`Tbot = Tf` is not
      !! stored).  `tau_x/tau_y` are the ice-ocean stress at the TRACER point per
      !! ICE area, previous step (see the module header, **Stress contract**).  Call from one `do concurrent (j, i)` over owned cells.
      !$acc routine seq
      integer, value :: i, j
      integer, intent(in) :: nx, ny, ncat, nilyr, nslyr
      real(wp), intent(in) :: dt, ustar_min, floediam, floeshape, rhow, cp_ocn, pi, puny
      real(wp), intent(in) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat)
      real(wp), intent(in) :: qicen(nx, ny, ncat, nilyr), qsnon(nx, ny, ncat, nslyr)
      real(wp), intent(in) :: frzmlt(nx, ny), sst(nx, ny), Tf(nx, ny), tau_x(nx, ny), tau_y(nx, ny)
      real(wp), intent(inout) :: fbot(nx, ny), rside(nx, ny), wlat(nx, ny)
      integer(int32), intent(inout) :: ierr(nx, ny)
      integer :: n, k
      logical :: neg
      real(wp) :: aice, etot, Tbot, xtmp, fb, rs, wl
      real(wp) :: etotn(ICE_NCAT_MAX)
      integer(int32) :: ie

      aice = 0.0_wp
      neg = .false.
      do n = 1, ncat
         aice = aice + aicen(i, j, n)
         if (aicen(i, j, n) < -puny) neg = .true.
      end do
      if (neg) then
         ierr(i, j) = ior(ierr(i, j), ICE_ERR_AICEN_NEG)
         fbot(i, j) = 0.0_wp
         rside(i, j) = 0.0_wp
         wlat(i, j) = 0.0_wp
         return
      end if
      do n = 1, ncat
         etotn(n) = 0.0_wp
      end do
      if (aice > puny .and. frzmlt(i, j) < 0.0_wp) then
         do n = 1, ncat
            etot = 0.0_wp
            do k = 1, nslyr
               etot = etot + qsnon(i, j, n, k)*vsnon(i, j, n)/real(nslyr, wp)
            end do
            do k = 1, nilyr
               etot = etot + qicen(i, j, n, k)*vicen(i, j, n)/real(nilyr, wp)
            end do
            etotn(n) = etot
         end do
      end if
      ie = 0_int32
      call ice_frzmlt_bottom_lateral(ncat, dt, aice, frzmlt(i, j), sst(i, j), Tf(i, j), tau_x(i, j), tau_y(i, j), &
                                     ustar_min, floediam, floeshape, rhow, cp_ocn, pi, puny, &
                                     etotn, Tbot, fb, rs, wl, xtmp, ie)
      fbot(i, j) = fb
      rside(i, j) = rs
      wlat(i, j) = wl
      if (ie /= 0_int32) ierr(i, j) = ior(ierr(i, j), ie)
   end subroutine ice_frzmlt_bottom_lateral_cell

   subroutine ice_basal_heat_step(nx, ny, ncat, nilyr, nslyr, i_lo, i_hi, j_lo, j_hi, dt, p, &
                                  aicen, vicen, vsnon, qicen, qsnon, frzmlt, sst, Tf, tau_x, tau_y, &
                                  fbot, rside, wlat, ierr)
      !! Driver over the owned cells `i_lo:i_hi, j_lo:j_hi`: one `do concurrent`
      !! calling the cell entry.  Every array must be device-present on a GPU
      !! build; the parameter bundle is read on the host (scalars go into the
      !! kernel), so it needs no mapping.
      integer, intent(in) :: nx, ny, ncat, nilyr, nslyr, i_lo, i_hi, j_lo, j_hi
      real(wp), intent(in) :: dt
      type(ice_params_t), intent(in) :: p
      real(wp), intent(in) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat)
      real(wp), intent(in) :: qicen(nx, ny, ncat, nilyr), qsnon(nx, ny, ncat, nslyr)
      real(wp), intent(in) :: frzmlt(nx, ny), sst(nx, ny), Tf(nx, ny), tau_x(nx, ny), tau_y(nx, ny)
      real(wp), intent(inout) :: fbot(nx, ny), rside(nx, ny), wlat(nx, ny)
      integer(int32), intent(inout) :: ierr(nx, ny)

      call basal_heat_impl(nx, ny, ncat, nilyr, nslyr, i_lo, i_hi, j_lo, j_hi, dt, &
                           p%thermo%ustar_min, p%thermo%floediam, p%thermo%floeshape, &
                           p%c%rhow, p%c%cp_ocn, p%c%pi, p%c%puny, &
                           aicen, vicen, vsnon, qicen, qsnon, frzmlt, sst, Tf, tau_x, tau_y, &
                           fbot, rside, wlat, ierr)
   end subroutine ice_basal_heat_step

   subroutine basal_heat_impl(nx, ny, ncat, nilyr, nslyr, i_lo, i_hi, j_lo, j_hi, dt, &
                              ustar_min, floediam, floeshape, rhow, cp_ocn, pi, puny, &
                              aicen, vicen, vsnon, qicen, qsnon, frzmlt, sst, Tf, tau_x, tau_y, &
                              fbot, rside, wlat, ierr)
      integer, intent(in) :: nx, ny, ncat, nilyr, nslyr, i_lo, i_hi, j_lo, j_hi
      real(wp), intent(in) :: dt, ustar_min, floediam, floeshape, rhow, cp_ocn, pi, puny
      real(wp), intent(in) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat)
      real(wp), intent(in) :: qicen(nx, ny, ncat, nilyr), qsnon(nx, ny, ncat, nslyr)
      real(wp), intent(in) :: frzmlt(nx, ny), sst(nx, ny), Tf(nx, ny), tau_x(nx, ny), tau_y(nx, ny)
      real(wp), intent(inout) :: fbot(nx, ny), rside(nx, ny), wlat(nx, ny)
      integer(int32), intent(inout) :: ierr(nx, ny)
      integer :: i, j

      do concurrent(j=j_lo:j_hi, i=i_lo:i_hi)
         call ice_frzmlt_bottom_lateral_cell(i, j, nx, ny, ncat, nilyr, nslyr, dt, &
                                             ustar_min, floediam, floeshape, rhow, cp_ocn, pi, puny, &
                                             aicen, vicen, vsnon, qicen, qsnon, frzmlt, sst, Tf, tau_x, tau_y, &
                                             fbot, rside, wlat, ierr)
      end do
   end subroutine basal_heat_impl

end module rdb_ice_basal_heat
