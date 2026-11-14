!! Ridging step: redistribution of ice area, volume, snow and tracers between thickness categories, the closing and
!! opening limiters and the convergence iteration (Thorndike et al. 1975; Lipscomb et al. 2007).
module rdb_ice_ridge_shift
   ! SPDX-License-Identifier: BSD-3-Clause
   !
   ! Ported from CICE/Icepack (CICE Consortium, https://github.com/CICE-Consortium).
   ! Upstream: Icepack 9a1f3e4 (Icepack 1.5.3+14), module icepack_mechred
   ! (ridge_ice, ridge_prep, ridge_shift).
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
   ! Elizabeth C. Hunke (LANL).
   ! Physics: Thorndike et al. (1975, JGR 80 4501); Hibler (1980, MWR 108 1943);
   ! Flato and Hibler (1995, JGR 100 18611); Lipscomb et al. (2007, JGR 112 C03S91).
   !
   !! **What one call does** (`ice_ridge_ice_cell`, cell `(i, j)` of the SoA state, IN PLACE).
   !!
   !! 1. *Closing and opening rates.*  `asum = aice0 + sum_n aicen(n)` (it differs from 1 after
   !!    transport).  The net closing rate is `closing_net = Cs rdg_shear + rdg_conv` (Flato and
   !!    Hibler 1995; `Cs = 0.25`); the advective divergence `divu_adv = (1 - asum)/dt`; for
   !!    `divu_adv < 0` the closing is raised to `-divu_adv`; the opening `opning = closing_net +
   !!    divu_adv` is then non-negative by construction and `(opning - closing_net) dt = 1 - asum`.
   !! 2. *One pass* (repeated up to `itermax`, production 20): `ridge_itd_cell` (participation
   !!    `apartic`, ridge geometry, `aksum`), the gross closing `closing_gross = closing_net/aksum`,
   !!    cut (together with the opening) so that no participating category, nor the open water, loses
   !!    more area than it has; the open water becomes `aice0 - apartic(0) closing_gross dt +
   !!    opning dt` (`< -puny`: `ICE_ERR_RIDGE_AICE0`, a smaller negative is clipped to 0).  Each
   !!    ridging category `n` (`aicen > puny`, `apartic(n) > 0`, `closing_gross > 0`) loses the area
   !!    `ardg1n = apartic(n) closing_gross dt` (more than `aicen + puny`: `ICE_ERR_RIDGE_OVERSHOOT`),
   !!    the proportional ice volume `virdgn` and snow volume `vsrdgn`; the ridge area is
   !!    `ardg2n = ardg1n/krdg(n)`.  Receiving category `m` gets the fractions `farea(n,m)`,
   !!    `fvol(n,m)` of the new ridge area and volume (`krdg_redist = 0` Hibler 1980: uniform
   !!    between `hrmin` and `hrmax`; `= 1` Lipscomb et al. 2007: negative exponential beyond
   !!    `hrmin`, exponent clamped at 100, the top category taking the tail), and the fraction
   !!    `fsnowrdg = 0.5` of the ridged snow volume.  The other half of the snow is thrown into the
   !!    ocean (`msnow_mlt = sum rhos vsrdgn (1 - fsnowrdg)`, `esnow_mlt = sum (1 - fsnowrdg)
   !!    vsrdgn qsno/nslyr`).  Tracers move with their extensive carriers (area tracer `Tsfc`:
   !!    `ardg1n` out, `ardg2n farea` in, so it is not conserved; ice-volume tracers `qice`, `sice`:
   !!    `virdgn` out, `virdgn fvol` in, conserved exactly; snow tracer `qsno`: `vsrdgn` out,
   !!    `vsrdgn fvol fsnowrdg` in), and EVERY category's intensive tracers are then rederived
   !!    from the new extensive content (empty area tracer: `Tf`; empty volume tracer: 0), exactly as
   !!    upstream does, so a category the pass did not touch still moves by round-off (`x a/a`).
   !! 3. *Convergence.*  `asum` is recomputed; within `puny` of 1 the cell is done, else the closing
   !!    and opening are reset from the remaining misfit (`max(0, -divu_adv)`, `max(0, divu_adv)`) and
   !!    the next pass runs on the UPDATED state.  Still not converged after `itermax` passes:
   !!    `ICE_ERR_RIDGE_NOCONV`.  Typical: one pass; two when a limiter fired.
   !! 4. *Booking* (on success only): `fresh += msnow_mlt/(ndtd dt)` (positive into the ocean, kg
   !!    m-2 s-1, per unit CELL area), `fhocn += esnow_mlt/(ndtd dt)` (negative: the ocean cools to
   !!    melt the snow), and the diagnostics `dardg1dt += ardg1/dt`, `dardg2dt += ardg2/dt`,
   !!    `dvirdgdt += virdg/dt` (m/s), `opening += aopen/dt` where `aopen` is the opening of the
   !!    LAST pass and `ardg1, ardg2, virdg, msnow_mlt, esnow_mlt` accumulate over passes.  Upstream
   !!    ASSIGNS the four diagnostics; the flux slot is zeroed once per thermodynamic step and every
   !!    producer adds, which is the same for `ndtd = 1` (the only supported value).  `fsalt` is not
   !!    touched (snow carries no salt, ice salt moves inside the pack).  Everything else reaching the
   !!    ocean after ridging comes from the clean-up that `ice_step_ridge` runs next.
   !!
   !! **Exact early-outs** (the only places that differ from upstream in the state).  (a) A QUIESCENT
   !! cell -- `closing_net == 0` and `opning == 0` after step 1, i.e. `asum == 1` exactly and no
   !! deformation -- returns with the state bit for bit unchanged, `niter = 0`; upstream runs one
   !! empty pass there, which only rederives the tracers (a round-off `x a/a` on every category and the
   !! normalisation of empty categories).  (b) An ICE-FREE cell -- every `aicen`, `vicen`, `vsnon` exactly 0
   !! and `asum == 1` exactly -- returns unchanged, `niter = 0`; upstream would move `aice0` by
   !! `aice0 - c dt + c dt`, one ulp.  Only the OPENING diagnostic is booked, formed as upstream does
   !! (net to gross closing through the open-water participation, the open-water limiter, the opening cut in
   !! step).  In a quiescent cell the tracers of an exactly empty category keep their stale values (upstream
   !! normalises them): every consumer weights a tracer by its zero area or volume, and the clean-up zaps
   !! anything in `(0, puny]`.  The supported `ndtd` is 1: the diagnostics are added once per call, so the
   !! driver refuses `ndtd /= 1`.  Both save the divergent 20-pass loop on the open
   !! ocean.  Neither reads the tracers.
   !!
   !! **Reordering of the redistribution** (`ridge_shift`).  Upstream loops source `n` outer and
   !! destination `nr` inner into a `(ntrcr, ncat)` extensive array.  Here, per destination `m`, the
   !! extensive value starts from its own old content and visits the sources `n = 1..ncat` in the same
   !! order (`n == m`: subtract the removal first, then add the receipt), so every destination sees the
   !! very same sequence of floating-point operations: the result is bitwise the upstream one
   !! (arithmetic only; the exponentials come from `rdb_sm_exp` and agree with libm to a few ulp, see the
   !! test).  The new tracer of destination `m` is written straight into the state, which would destroy
   !! the donor value a later destination needs; so the tracers are processed ONE AT A TIME: all donor
   !! values of that tracer are gathered first (`xd(1:ncat)`), then every destination written.  The
   !! `farea(n,m)`, `fvol(n,m)` fractions (one exponential pair each) are evaluated once into two
   !! `(ICE_NCAT_MAX, ICE_NCAT_MAX)` tables and reused by every tracer.  `hin_max` is `intent(in)`:
   !! upstream overwrites the top bound with `1e8`; here the top category's bound is a local `RIDGE_BIG`
   !! (the redistribution `krdg_redist = 1` never reads it; the clean-up's rebin never reads it either).
   !!
   !! **Errors.**  Non-finite `aice0`, `Tf`, `dt`, `rdg_conv`, `rdg_shear`, any `aicen`, `vicen`, `vsnon`
   !! (checked first, state untouched), or any tracer of a cell that is not an early-out (checked next,
   !! state untouched): `ICE_ERR_NONFINITE`.  A comparison with NaN is false, so these would otherwise
   !! be skipped silently.  A bad dimension (`ncat` outside `1..ICE_NCAT_MAX`, `nilyr`, `nslyr < 1`,
   !! `itermax < 1`, `dt <= 0`) is also `ICE_ERR_NONFINITE`.  `ICE_ERR_RIDGE_AICE0`,
   !! `ICE_ERR_RIDGE_OVERSHOOT`, `ICE_ERR_RIDGE_DEGENERATE` (new: `aksum <= 0` with `closing_net > 0`, an
   !! inconsistent empty cell where upstream divides by zero) and `ICE_ERR_RIDGE_NOCONV` are FATAL and
   !! UNRECOVERABLE: the fluxes and diagnostics are not booked and the state is that of the last
   !! COMPLETED pass (the first three are detected before the pass writes anything, so the state is not
   !! half updated; NOCONV leaves the state of the last pass), but no consumer may use the cell.  With
   !! `closing_net <= 0` and `aksum <= 0` (opening only, an empty cell) the closing is taken as 0.
   !!
   !! **`ICE_ERR_RIDGE_OVERSHOOT` is unreachable through this interface** (kept as upstream's guard):
   !! the limiter loop cuts `closing_gross` until `apartic(n) closing_gross dt <= aicen(n)` for every
   !! participating category (and only decreases it afterwards), so `ardg1n` exceeds `aicen + puny`
   !! only by a few ulp of `aicen`, never by `puny`; the ridging test (`aicen > puny`) is the limiter's
   !! own.  Non-finite input never gets that far.  There is no input that sets it, hence no test.
   !!
   !! **Per-thread memory** (fixed `ICE_NCAT_MAX = 8` locals, no automatic arrays): 6 working vectors
   !! of 8, the ridge geometry (`apartic` 9, `hrmin`, `hrmax`, `hrexp`, `krdg` 8 each), 4 donor-amount
   !! vectors, the two 8 x 8 fraction tables, one tracer-donor vector, one `ridge_itd_cell` work vector
   !! (10) and a flag vector: about 270 doubles.  See the test for the measured stack.
   !!
   !! **Upstream branches dropped, and why they are inert.**  The conservation-check block
   !! (`conserv_check`, `column_sum`, `column_conservation_check`; OM3 sets `conserv_check = .false.`;
   !! the tests carry the same closure instead); the `opening`/`closing` argument branch (EAP dynamics
   !! only); the optional per-category diagnostics (`aparticn`, `krdgn`, `aredistn`, `vredistn`,
   !! `dardg*ndt`, `dvirdgndt`, `araftn`, `vraftn`, `mraftn`) which the flux slot does not carry; level-ice
   !! tracers (`alvl`, `vlvl`: `tr_lvl` is off, so the `nt_alvl`/`nt_vlvl` exclusion is vacuous) and
   !! melt ponds (`tr_pond_*`, `mpond`, `fpond`, `dpnd_ridge`: ponds are not set in the OM3 configuration,
   !! default off), the brine-height tracer (`tr_brine`) and the extra tracer strata
   !! (`n_trcr_strata`, none without `tr_brine`/ponds), aerosols (`tr_aero`, `maero`), isotopes
   !! (`tr_iso`, `miso`), biogeochemistry (`z_tracers`, `nbtrcr`, `mbio`), the snow-density and grain
   !! tracers (`snwredist`, `snwgrain`), and the final `abs(asum - 1) > puny` abort (unreachable: the loop
   !! only leaves through the test that makes it false).  The `fiso_ocn`, `faero_ocn`, `flux_bio`
   !! arguments go with them.  Deviation: `aksum <= 0` is an error bit (upstream divides by zero).
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
   use, intrinsic :: iso_fortran_env, only: int32, int64
   use rdb_constants, only: wp
   use rdb_ice_dims, only: ICE_NCAT_MAX
   use rdb_ice_errors, only: ICE_ERR_NONFINITE, ICE_ERR_RIDGE_NOCONV, ICE_ERR_RIDGE_AICE0, ICE_ERR_RIDGE_OVERSHOOT, &
                             ICE_ERR_RIDGE_DEGENERATE
   use rdb_ice_params, only: ice_ridge_params_t, ice_ntrcr, ice_trcr_get, ice_trcr_set, ice_trcr_depend, ice_intensive
   use rdb_ice_ridge_itd, only: ridge_itd_cell, RIDGE_CS, RIDGE_FSNOWRDG, RIDGE_EXP_ARGMAX, RIDGE_BIG
   implicit none
   private

   public :: ice_ridge_ice_cell, ridge_redist_fractions

#include "rdb_safe_math_consts.inc"

contains

#include "rdb_safe_math_core.inc"

   pure subroutine ridge_redist_fractions(kredist, nr, ncat, hrmin, hrmax, hrexp, hlo, hhi, farea, fvol)
      !! Fractions `farea`, `fvol` of the new ridge area and volume (of a ridging category with `hrmin`, `hrmax`
      !! [`kredist = 0`] or `hrmin`, `hrexp` [`kredist = 1`]) that land in category `nr` of `ncat`, whose
      !! bounds are `hlo = hin_max(nr-1)` and `hhi = hin_max(nr)` (the top category's `hhi` is replaced by
      !! `RIDGE_BIG` here).  Over `nr = 1..ncat` each sums to 1.
      !$acc routine seq
      integer, value :: kredist, nr, ncat
      real(wp), intent(in) :: hrmin, hrmax, hrexp, hlo, hhi
      real(wp), intent(out) :: farea, fvol
      real(wp) :: hL, hR, dhr, dhr2, hi1, hexp, expL, expR, expL_arg, expR_arg, hb

      if (kredist == 0) then
         hb = hhi
         if (nr == ncat) hb = RIDGE_BIG
         dhr = hrmax - hrmin
         dhr2 = hrmax*hrmax - hrmin*hrmin
         if (hrmin >= hb .or. hrmax <= hlo) then
            hL = 0.0_wp
            hR = 0.0_wp
         else
            hL = max(hrmin, hlo)
            hR = min(hrmax, hb)
         end if
         farea = (hR - hL)/dhr
         fvol = (hR*hR - hL*hL)/dhr2
      else
         hi1 = hrmin
         hexp = hrexp
         if (nr < ncat) then
            if (hi1 >= hhi) then
               farea = 0.0_wp
               fvol = 0.0_wp
            else
               hL = max(hi1, hlo)
               hR = hhi
               expL_arg = min(((hL - hi1)/hexp), RIDGE_EXP_ARGMAX)
               expR_arg = min(((hR - hi1)/hexp), RIDGE_EXP_ARGMAX)
               expL = rdb_sm_exp(-(expL_arg))
               expR = rdb_sm_exp(-(expR_arg))
               farea = expL - expR
               fvol = ((hL + hexp)*expL - (hR + hexp)*expR)/(hi1 + hexp)
            end if
         else
            hL = max(hi1, hlo)
            expL_arg = min(((hL - hi1)/hexp), RIDGE_EXP_ARGMAX)
            expL = rdb_sm_exp(-(expL_arg))
            farea = expL
            fvol = (hL + hexp)*expL/(hi1 + hexp)
         end if
      end if
   end subroutine ridge_redist_fractions

   pure subroutine ice_ridge_ice_cell(i, j, nx, ny, ncat, nilyr, nslyr, rp, puny, rhos, dt, ndtd, itermax, hin_max, Tf, &
                                      rdg_conv, rdg_shear, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, aice0, &
                                      fresh, fhocn, dardg1dt, dardg2dt, dvirdgdt, opening, ierr, niter)
      !! Ridging of cell `(i, j)` of the SoA state (`aicen, vicen, vsnon (nx,ny,ncat)`, `Tsfcn`, `qicen, sicen
      !! (nx,ny,ncat,nilyr)`, `qsnon (nx,ny,ncat,nslyr)`, `aice0 (nx,ny)`) IN PLACE; see the module note for the
      !! algorithm, the early-outs and the errors.  `dt` is the dynamics step, `ndtd` the number of dynamics
      !! subcycles per thermodynamic step (the fluxes are per second of `ndtd dt`), `itermax` the pass cap
      !! (production `RIDGE_NITERMAX = 20`; a smaller value is the non-convergence test hook), `hin_max(0:ncat)` the
      !! category bounds, `Tf` the cell freezing point (the value of an empty area tracer), `rdg_conv`, `rdg_shear`
      !! [1/s] the deformation rates.  The booked outputs are `fresh`, `fhocn` and the four diagnostics `dardg1dt`,
      !! `dardg2dt`, `dvirdgdt`, `opening` (all added to the `(nx,ny)` cell arrays) and `niter(i,j)`-style scalar
      !! `niter` (0 = early-out, else the number of passes run).  `ierr(i,j)` receives the bits.
      !$acc routine seq
      integer, value :: i, j, ndtd, itermax
      integer, intent(in) :: nx, ny, ncat, nilyr, nslyr
      type(ice_ridge_params_t), intent(in) :: rp
      real(wp), intent(in) :: puny, rhos, dt
      real(wp), intent(in) :: hin_max(0:ncat)
      real(wp), intent(in) :: Tf, rdg_conv, rdg_shear
      real(wp), intent(inout) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat)
      real(wp), intent(inout) :: Tsfcn(nx, ny, ncat), qicen(nx, ny, ncat, nilyr)
      real(wp), intent(inout) :: sicen(nx, ny, ncat, nilyr), qsnon(nx, ny, ncat, nslyr)
      real(wp), intent(inout) :: aice0(nx, ny)
      real(wp), intent(inout) :: fresh(nx, ny), fhocn(nx, ny), dardg1dt(nx, ny), dardg2dt(nx, ny), dvirdgdt(nx, ny)
      real(wp), intent(inout) :: opening(nx, ny)
      integer(int32), intent(inout) :: ierr(nx, ny)
      integer, intent(out) :: niter

      integer :: n, m, k, it, kdep, kredist, pass, ntr, nrdg
      logical :: iterate_ridging, empty
      logical :: rdg(ICE_NCAT_MAX)
      real(wp) :: an(ICE_NCAT_MAX), vn(ICE_NCAT_MAX), sn(ICE_NCAT_MAX)
      real(wp) :: an2(ICE_NCAT_MAX), vn2(ICE_NCAT_MAX), sn2(ICE_NCAT_MAX)
      real(wp) :: apartic(0:ICE_NCAT_MAX), hrmin(ICE_NCAT_MAX), hrmax(ICE_NCAT_MAX), hrexp(ICE_NCAT_MAX), krdg(ICE_NCAT_MAX)
      real(wp) :: ardg1n(ICE_NCAT_MAX), ardg2n(ICE_NCAT_MAX), virdgn(ICE_NCAT_MAX), vsrdgn(ICE_NCAT_MAX)
      real(wp) :: fa(ICE_NCAT_MAX, ICE_NCAT_MAX), fv(ICE_NCAT_MAX, ICE_NCAT_MAX), xd(ICE_NCAT_MAX)
      real(wp) :: a0, asum, aksum, closing_net, divu_adv, opning, closing_gross, tmpfac, wk1, aopen
      real(wp) :: ardg1, ardg2, virdg, msnow_mlt, esnow_mlt, afrac, esrdgn, poison, dti
      real(wp) :: x_a, x_v, x_s, ext, work, base, rem, coef, empty_value, newv, farea_, fvol_, hlo, hhi
      real(wp), parameter :: ONE = 1.0_wp, ZERO = 0.0_wp

      niter = 0
      if (ncat < 1 .or. ncat > ICE_NCAT_MAX .or. nilyr < 1 .or. nslyr < 1 .or. itermax < 1 .or. ndtd < 1) then
         ierr(i, j) = ior(ierr(i, j), ICE_ERR_NONFINITE)
         return
      end if

      ! state and rates, finiteness first (a NaN would make every later comparison false)
      a0 = aice0(i, j)
      poison = a0 + rdg_conv + rdg_shear + Tf + dt
      do n = 1, ncat
         an(n) = aicen(i, j, n)
         vn(n) = vicen(i, j, n)
         sn(n) = vsnon(i, j, n)
         poison = poison + an(n) + vn(n) + sn(n)
      end do
      if (.not. ieee_is_finite(poison) .or. .not. (dt > ZERO)) then
         ierr(i, j) = ior(ierr(i, j), ICE_ERR_NONFINITE)
         return
      end if

      ! opening and closing rates (ridge_prep)
      asum = a0
      do n = 1, ncat
         asum = asum + an(n)
      end do
      closing_net = RIDGE_CS*rdg_shear + rdg_conv
      divu_adv = (ONE - asum)/dt
      if (divu_adv < ZERO) closing_net = max(closing_net, -divu_adv)
      opning = closing_net + divu_adv

      ! exact early-outs: nothing closes or opens; no ice and the area is exactly 1
      if (closing_net == ZERO .and. opning == ZERO) return
      if (asum == ONE) then
         empty = .true.
         do n = 1, ncat
            if (an(n) /= ZERO .or. vn(n) /= ZERO .or. sn(n) /= ZERO) empty = .false.
         end do
         if (empty) then
            ! Ice-free cell under deformation: the state is left alone (upstream would move `aice0` by an ulp), but
            ! the opening diagnostic is booked exactly as upstream forms it: the net closing becomes the gross
            ! closing through the open-water participation, is cut so that no more than the open water is removed,
            ! and the opening is cut in step (`opening = min(c, aice0/dt)` for a convergence `c`).
            call ridge_itd_cell(ncat, rp, puny, a0, an, vn, aksum, apartic, hrmin, hrmax, hrexp, krdg)
            if (aksum > ZERO) then
               closing_gross = closing_net/aksum
               if (apartic(0) > ZERO) then
                  wk1 = apartic(0)*closing_gross*dt
                  if (wk1 > a0) then
                     tmpfac = a0/wk1
                     opning = opning*tmpfac
                  end if
               end if
               opening(i, j) = opening(i, j) + (opning*dt)*(ONE/dt)
            end if
            return
         end if
      end if

      ! tracers of a ridging cell must be finite too
      poison = ZERO
      do n = 1, ncat
         poison = poison + Tsfcn(i, j, n)
         do k = 1, nilyr
            poison = poison + qicen(i, j, n, k) + sicen(i, j, n, k)
         end do
         do k = 1, nslyr
            poison = poison + qsnon(i, j, n, k)
         end do
      end do
      if (.not. ieee_is_finite(poison)) then
         ierr(i, j) = ior(ierr(i, j), ICE_ERR_NONFINITE)
         return
      end if

      kredist = rp%krdg_redist
      ntr = ice_ntrcr(nilyr, nslyr)
      ardg1 = ZERO
      ardg2 = ZERO
      virdg = ZERO
      aopen = ZERO
      msnow_mlt = ZERO
      esnow_mlt = ZERO

      do pass = 1, itermax
         niter = pass

         ! thickness distribution of the ridging ice and of the new ridges
         call ridge_itd_cell(ncat, rp, puny, a0, an, vn, aksum, apartic, hrmin, hrmax, hrexp, krdg)

         ! gross closing rate; the ridging step proper (ridge_shift)
         closing_gross = ZERO
         if (aksum > ZERO) closing_gross = closing_net/aksum
         if (.not. (aksum > ZERO) .and. closing_net > ZERO) then
            ierr(i, j) = ior(ierr(i, j), ICE_ERR_RIDGE_DEGENERATE)
            return
         end if

         ! no more than all of the open water or of a category is removed; the opening is cut in step
         if (apartic(0) > ZERO) then
            wk1 = apartic(0)*closing_gross*dt
            if (wk1 > a0) then
               tmpfac = a0/wk1
               closing_gross = closing_gross*tmpfac
               opning = opning*tmpfac
            end if
         end if
         do n = 1, ncat
            if (an(n) > puny .and. apartic(n) > ZERO) then
               wk1 = apartic(n)*closing_gross*dt
               if (wk1 > an(n)) then
                  tmpfac = an(n)/wk1
                  closing_gross = closing_gross*tmpfac
                  opning = opning*tmpfac
               end if
            end if
         end do

         ! open water
         a0 = a0 - apartic(0)*closing_gross*dt + opning*dt
         if (a0 < -puny) then
            ierr(i, j) = ior(ierr(i, j), ICE_ERR_RIDGE_AICE0)
            return
         end if
         if (a0 < ZERO) a0 = ZERO
         aopen = opning*dt

         ! area, volume and snow of the ridging ice of each category
         nrdg = 0
         do n = 1, ncat
            rdg(n) = .false.
            ardg1n(n) = ZERO
            ardg2n(n) = ZERO
            virdgn(n) = ZERO
            vsrdgn(n) = ZERO
            if (an(n) > puny .and. apartic(n) > ZERO .and. closing_gross > ZERO) then
               ardg1n(n) = apartic(n)*closing_gross*dt
               if (ardg1n(n) > an(n) + puny) then
                  ierr(i, j) = ior(ierr(i, j), ICE_ERR_RIDGE_OVERSHOOT)
                  return
               end if
               ardg1n(n) = min(an(n), ardg1n(n))
               ardg2n(n) = ardg1n(n)/krdg(n)
               afrac = ardg1n(n)/an(n)
               virdgn(n) = vn(n)*afrac
               vsrdgn(n) = sn(n)*afrac
               rdg(n) = .true.
               nrdg = nrdg + 1
               ardg1 = ardg1 + ardg1n(n)
               ardg2 = ardg2 + ardg2n(n)
               virdg = virdg + virdgn(n)
               ! the part of the snow that does not stay on the ridge goes into the ocean
               msnow_mlt = msnow_mlt + rhos*vsrdgn(n)*(ONE - RIDGE_FSNOWRDG)
               do k = 1, nslyr
                  esrdgn = vsrdgn(n)*qsnon(i, j, n, k)/real(nslyr, wp)
                  esnow_mlt = esnow_mlt + esrdgn*(ONE - RIDGE_FSNOWRDG)
               end do
            end if
         end do

         ! new area, ice volume and snow volume of every category (destination order)
         do m = 1, ncat
            x_a = an(m)
            x_v = vn(m)
            x_s = sn(m)
            do n = 1, ncat
               if (rdg(n)) then
                  if (n == m) then
                     x_a = x_a - ardg1n(n)
                     x_v = x_v - virdgn(n)
                     x_s = x_s - vsrdgn(n)
                  end if
                  hlo = hin_max(m - 1)
                  hhi = hin_max(m)
                  call ridge_redist_fractions(kredist, m, ncat, hrmin(n), hrmax(n), hrexp(n), hlo, hhi, farea_, fvol_)
                  fa(n, m) = farea_
                  fv(n, m) = fvol_
                  x_a = x_a + farea_*ardg2n(n)
                  x_v = x_v + fvol_*virdgn(n)
                  x_s = x_s + fvol_*vsrdgn(n)*RIDGE_FSNOWRDG
               end if
            end do
            an2(m) = x_a
            vn2(m) = x_v
            sn2(m) = x_s
         end do

         ! tracers, one at a time: gather every donor value, then write every destination
         do it = 1, ntr
            kdep = ice_trcr_depend(it, nilyr, nslyr)
            empty_value = ZERO
            if (it == 1) empty_value = Tf
            do n = 1, ncat
               xd(n) = ice_trcr_get(it, i, j, n, nx, ny, ncat, nilyr, nslyr, Tsfcn, qicen, sicen, qsnon)
            end do
            do m = 1, ncat
               if (kdep == 0) then
                  base = an(m)
               else if (kdep == 1) then
                  base = vn(m)
               else
                  base = sn(m)
               end if
               ext = xd(m)*base
               do n = 1, ncat
                  if (rdg(n)) then
                     if (n == m) then
                        if (kdep == 0) then
                           rem = ardg1n(n)
                        else if (kdep == 1) then
                           rem = virdgn(n)
                        else
                           rem = vsrdgn(n)
                        end if
                        work = rem*xd(n)
                        ext = ext - work
                     end if
                     if (kdep == 0) then
                        coef = ardg2n(n)*fa(n, m)
                     else if (kdep == 1) then
                        coef = virdgn(n)*fv(n, m)
                     else
                        coef = vsrdgn(n)*fv(n, m)*RIDGE_FSNOWRDG
                     end if
                     work = coef*xd(n)
                     ext = ext + work
                  end if
               end do
               newv = ice_intensive(ext, kdep, an2(m), vn2(m), sn2(m), puny, empty_value)
               call ice_trcr_set(it, i, j, m, nx, ny, ncat, nilyr, nslyr, newv, Tsfcn, qicen, sicen, qsnon)
            end do
         end do

         ! commit the pass
         do n = 1, ncat
            an(n) = an2(n)
            vn(n) = vn2(n)
            sn(n) = sn2(n)
            aicen(i, j, n) = an(n)
            vicen(i, j, n) = vn(n)
            vsnon(i, j, n) = sn(n)
         end do
         aice0(i, j) = a0

         ! total area 1?  If not (a limiter cut the rates), ridge again with the rates of the misfit
         asum = a0
         do n = 1, ncat
            asum = asum + an(n)
         end do
         if (abs(asum - ONE) < puny) then
            iterate_ridging = .false.
            closing_net = ZERO
            opning = ZERO
         else
            iterate_ridging = .true.
            divu_adv = (ONE - asum)/dt
            closing_net = max(ZERO, -divu_adv)
            opning = max(ZERO, divu_adv)
         end if
         if (.not. iterate_ridging) exit
         if (pass == itermax) then
            ierr(i, j) = ior(ierr(i, j), ICE_ERR_RIDGE_NOCONV)
            return
         end if
      end do

      ! diagnostics, then the snow thrown into the ocean (per second of the thermodynamic step)
      dti = ONE/dt
      dardg1dt(i, j) = dardg1dt(i, j) + ardg1*dti
      dardg2dt(i, j) = dardg2dt(i, j) + ardg2*dti
      dvirdgdt(i, j) = dvirdgdt(i, j) + virdg*dti
      opening(i, j) = opening(i, j) + aopen*dti
      dti = ONE/(real(ndtd, wp)*dt)
      fresh(i, j) = fresh(i, j) + msnow_mlt*dti
      fhocn(i, j) = fhocn(i, j) + esnow_mlt*dti
   end subroutine ice_ridge_ice_cell

end module rdb_ice_ridge_shift
