!! Ridging preparation (participation and redistribution functions) and ice strength.
module rdb_ice_ridge_itd
   ! SPDX-License-Identifier: BSD-3-Clause
   !
   ! Ported from CICE/Icepack (CICE Consortium, https://github.com/CICE-Consortium).
   ! Upstream: Icepack 9a1f3e4 (Icepack 1.5.3+14), module icepack_mechred
   ! (asum_ridging, ridge_itd, icepack_ice_strength).
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
   ! Physics: Thorndike et al. (1975, JGR 80 4501); Rothrock (1975, JGR 80 4514);
   ! Hibler (1979, JPO 9 817; 1980, MWR 108 1943); Flato and Hibler (1995, JGR
   ! 100 18611); Lipscomb et al. (2007, JGR 112 C03S91).
   !
   !! **Participation function.**  `apartic(n)` is the fraction of the area
   !! lost to ridging and closing that comes from category `n` (`n = 0` is open
   !! water, which takes part only when `aice0 > puny`; categories with
   !! `aicen <= puny` are skipped).  With `G(n)` the normalised cumulative area
   !! over categories `0..n` (`G(ncat) = 1`):
   !!   * `krdg_partic = 0` (Thorndike et al. 1975): `b(G) = (2/G*)(1 - G/G*)`
   !!     for `G < G*` (`G* = 0.15`), zero above;
   !!   * `krdg_partic = 1` (Lipscomb et al. 2007): `b(G) = exp(-G/a*)/(a*(1 -
   !!     exp(-1/a*)))`, `a* = 0.05`; `apartic(n) = (e(n-1) - e(n))/(1 -
   !!     exp(-1/a*))`, which telescopes to a sum of exactly 1.
   !!
   !! **Redistribution.**  A ridging category of mean thickness `hi` raises
   !! ridges no thinner than `hrmin = min(2 hi, hi + maxraft)` (rafting for thin
   !! ice).  `krdg_redist = 0` (Hibler 1980, modified): uniform on
   !! `[hrmin, hrmax]`, `hrmax = max(2 sqrt(Hstar hi), hrmin + puny)`;
   !! `krdg_redist = 1` (Lipscomb et al. 2007): negative exponential beyond
   !! `hrmin` with e-folding thickness `hrexp = mu_rdg sqrt(hi)` (`hi` floored
   !! at `puny`).  `krdg = (mean ridge thickness)/hi`, and `aksum = apartic(0) +
   !! sum_n apartic(n)(1 - 1/krdg(n))` is the net area removed per unit area
   !! participating.  Unlike upstream, `hin_max` is not touched: the thickness
   !! categories enter only the later `ridge_shift` step, not this preparation.
   !!
   !! **Strength.**  `kstrength = 1` (Rothrock 1975): the potential-energy
   !! change integrated over those functions, `P = Cf Cp sum_n apartic(n)
   !! (h2rdg(n)/krdg(n) - hi(n)^2) / aksum`, `h2rdg` the mean `h^2` of the new
   !! ridge.  `kstrength = 0` (Hibler 1979): `P = Pstar vice exp(-Cstar (1 -
   !! aice))`.  An empty or inconsistent cell with `aksum <= 0` (upstream divides
   !! by zero) gets `strength = 0` and `ICE_ERR_RIDGE_DEGENERATE` (set by the
   !! strength kernel with no closing involved, unlike the ridging step's use of
   !! the bit).  A strength of exactly zero means free drift: no consumer may
   !! divide by it.
   !!
   !! All procedures are `pure`, `!$acc routine seq`, explicit-shape, with
   !! fixed-size (`ICE_NCAT_MAX`) locals and no allocation.  `ridge_itd_cell`
   !! needs finite inputs (the strength cell guards them).  `exp` is the
   !! reproducible `rdb_sm_exp`, so values agree with a libm evaluation to a
   !! few ulp, not bitwise.  The constants upstream fixes are module parameters
   !! (`RIDGE_*`); nothing derived is stored in `ice_ridge_params_t`, so nothing
   !! can go stale: `1/G*`, `1/a*` are parameters, the exponential normaliser is one
   !! `exp` per `ridge_itd_cell` call outside the category loop, and `Cp` is three
   !! flops from the constants bundle in `ice_strength_cell`.
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
   use, intrinsic :: iso_fortran_env, only: int32, int64
   use rdb_constants, only: wp
   use rdb_ice_dims, only: ICE_NCAT_MAX
   use rdb_ice_errors, only: ICE_ERR_NONFINITE, ICE_ERR_RIDGE_DEGENERATE
   use rdb_ice_params, only: ice_ridge_params_t, ice_constants_t
   implicit none
   private

   public :: ice_asum_cell, ridge_itd_cell, ice_strength_cell
   public :: RIDGE_CS, RIDGE_FSNOWRDG, RIDGE_GSTAR, RIDGE_ASTAR, RIDGE_MAXRAFT, RIDGE_HSTAR
   public :: RIDGE_EXP_ARGMAX, RIDGE_NITERMAX, RIDGE_BIG

#include "rdb_safe_math_consts.inc"

   real(wp), parameter :: RIDGE_CS = 0.25_wp
      !! Fraction of shear energy contributing to ridging (fixed upstream, not a knob).
   real(wp), parameter :: RIDGE_FSNOWRDG = 0.5_wp
      !! Snow fraction that survives ridging.
   real(wp), parameter :: RIDGE_GSTAR = 0.15_wp
      !! Maximum `G(h)` that participates, `krdg_partic = 0`.
   real(wp), parameter :: RIDGE_ASTAR = 0.05_wp
      !! E-folding scale of the participation, `krdg_partic = 1`.
   real(wp), parameter :: RIDGE_MAXRAFT = 1.0_wp
      !! Maximum thickness of ice that rafts [m].
   real(wp), parameter :: RIDGE_HSTAR = 25.0_wp
      !! Mean thickness of ridged ice, `krdg_redist = 0` [m].
   real(wp), parameter :: RIDGE_EXP_ARGMAX = 100.0_wp
      !! Exponent clamp of the redistribution functions (used by the shift step).
   integer, parameter :: RIDGE_NITERMAX = 20
      !! Ridging iteration limit (used by the ridging step).
   real(wp), parameter :: RIDGE_BIG = 1.0e8_wp
      !! Stand-in for the top category bound (upstream overwrites `hin_max(ncat)`).
   real(wp), parameter :: GSTARI = 1.0_wp/RIDGE_GSTAR
   real(wp), parameter :: ASTARI = 1.0_wp/RIDGE_ASTAR

   real(wp), parameter :: P333 = 1.0_wp/3.0_wp
      !! One third, multiplied exactly as upstream does (never divided by 3).

contains

#include "rdb_safe_math_core.inc"

   pure function ice_asum_cell(ncat, aicen, aice0) result(asum)
      !! Ice plus open-water area `aice0 + sum aicen` (may exceed 1 after
      !! transport), summed in the upstream order.
      !$acc routine seq
      integer, intent(in) :: ncat
      real(wp), intent(in) :: aice0
      real(wp), intent(in) :: aicen(ncat)
      real(wp) :: asum
      integer :: n

      asum = aice0
      do n = 1, ncat
         asum = asum + aicen(n)
      end do
   end function ice_asum_cell

   pure subroutine ridge_itd_cell(ncat, rp, puny, aice0, aicen, vicen, aksum, apartic, hrmin, hrmax, hrexp, krdg)
      !! Participation function `apartic(0:ncat)`, ridge geometry `hrmin`,
      !! `hrmax` (`krdg_redist = 0`), `hrexp` (`krdg_redist = 1`), `krdg` of each
      !! category and `aksum` (see the module note).  `ncat <= ICE_NCAT_MAX`;
      !! inputs finite.
      !$acc routine seq
      integer, intent(in) :: ncat
      type(ice_ridge_params_t), intent(in) :: rp
      real(wp), intent(in) :: puny, aice0
      real(wp), intent(in) :: aicen(ncat), vicen(ncat)
      real(wp), intent(out) :: aksum
      real(wp), intent(out) :: apartic(0:ncat), hrmin(ncat), hrmax(ncat), hrexp(ncat), krdg(ncat)

      integer :: n, kpartic, kredist
      real(wp) :: Gsum(-1:ICE_NCAT_MAX)
      real(wp) :: work, hi, hrmean, xtmp, mu

      kpartic = rp%krdg_partic
      kredist = rp%krdg_redist
      mu = rp%mu_rdg

      ! cumulative area G(n) over categories 0..n, small categories ignored
      Gsum(-1) = 0.0_wp
      if (aice0 > puny) then
         Gsum(0) = aice0
      else
         Gsum(0) = Gsum(-1)
      end if
      apartic(0) = 0.0_wp
      do n = 1, ncat
         apartic(n) = 0.0_wp
         hrmin(n) = 0.0_wp
         hrmax(n) = 0.0_wp
         hrexp(n) = 0.0_wp
         krdg(n) = 1.0_wp
         if (aicen(n) > puny) then
            Gsum(n) = Gsum(n - 1) + aicen(n)
         else
            Gsum(n) = Gsum(n - 1)
         end if
      end do

      if (Gsum(ncat) > 0.0_wp) then
         work = 1.0_wp/Gsum(ncat)
         do n = 0, ncat
            Gsum(n) = Gsum(n)*work
         end do
      end if

      if (kpartic == 0) then
         ! Thorndike et al. (1975): b(G) = (2/G*)(1 - G/G*), integrated over each category
         do n = 0, ncat
            if (Gsum(n) < RIDGE_GSTAR) then
               apartic(n) = GSTARI*(Gsum(n) - Gsum(n - 1))*(2.0_wp - GSTARI*(Gsum(n - 1) + Gsum(n)))
            else if (Gsum(n - 1) < RIDGE_GSTAR) then
               apartic(n) = GSTARI*(RIDGE_GSTAR - Gsum(n - 1))*(2.0_wp - GSTARI*(RIDGE_GSTAR + Gsum(n - 1)))
            end if
         end do
      else if (kpartic == 1) then
         ! Lipscomb et al. (2007): b(G) = exp(-G/a*); Gsum is reused for the exponentials.
         ! exp(-0) = 1 exactly, so the n = -1 term is the normalisation itself.
         ! normaliser 1/(1 - exp(-1/a*)): one exp per call, outside the category loop
         xtmp = 1.0_wp/(1.0_wp - rdb_sm_exp(-ASTARI))
         Gsum(-1) = xtmp
         do n = 0, ncat
            Gsum(n) = rdb_sm_exp(-(Gsum(n)*ASTARI))*xtmp
            apartic(n) = Gsum(n - 1) - Gsum(n)
         end do
      end if

      if (kredist == 0) then
         ! Hibler (1980), modified: uniform on [hrmin, hrmax]
         do n = 1, ncat
            if (aicen(n) > puny) then
               hi = vicen(n)/aicen(n)
               hrmin(n) = min(2.0_wp*hi, hi + RIDGE_MAXRAFT)
               hrmax(n) = 2.0_wp*sqrt(RIDGE_HSTAR*hi)
               hrmax(n) = max(hrmax(n), hrmin(n) + puny)
               hrmean = 0.5_wp*(hrmin(n) + hrmax(n))
               krdg(n) = hrmean/hi
            end if
         end do
      else
         ! Lipscomb et al. (2007): negative exponential beyond hrmin
         do n = 1, ncat
            if (aicen(n) > puny) then
               hi = vicen(n)/aicen(n)
               hi = max(hi, puny)
               hrmin(n) = min(2.0_wp*hi, hi + RIDGE_MAXRAFT)
               hrexp(n) = mu*sqrt(hi)
               krdg(n) = (hrmin(n) + hrexp(n))/hi
            end if
         end do
      end if

      ! net ice area removed per unit area participating
      aksum = apartic(0)
      do n = 1, ncat
         aksum = aksum + apartic(n)*(1.0_wp - 1.0_wp/krdg(n))
      end do
   end subroutine ridge_itd_cell

   pure subroutine ice_strength_cell(ncat, rp, c, aice, vice, aice0, aicen, vicen, strength, ierr)
      !! Ice strength `strength` [N/m] of one column (`kstrength` 0 Hibler,
      !! otherwise 1 Rothrock; see the module note).  `aice`, `vice` are the
      !! cell totals, `aice0` the open-water fraction, `aicen`/`vicen` the
      !! category area and volume; `c` supplies `puny` and the constants of `Cp`.  Non-finite input: `strength` is non-finite
      !! and `ICE_ERR_NONFINITE` is ORed into `ierr`.
      !$acc routine seq
      integer, intent(in) :: ncat
      type(ice_ridge_params_t), intent(in) :: rp
      type(ice_constants_t), intent(in) :: c
      real(wp), intent(in) :: aice, vice, aice0
      real(wp), intent(in) :: aicen(ncat), vicen(ncat)
      real(wp), intent(out) :: strength
      integer(int32), intent(inout) :: ierr

      integer :: n
      real(wp) :: poison, aksum, hi, h2rdg, dh2rdg, ssum, puny, Cp
      real(wp) :: apartic(0:ICE_NCAT_MAX), hrmin(ICE_NCAT_MAX), hrmax(ICE_NCAT_MAX), hrexp(ICE_NCAT_MAX), krdg(ICE_NCAT_MAX)

      puny = c%puny
      poison = aice + vice + aice0
      do n = 1, ncat
         poison = poison + aicen(n) + vicen(n)
      end do
      if (.not. ieee_is_finite(poison)) then
         strength = poison
         ierr = ior(ierr, ICE_ERR_NONFINITE)
         return
      end if

      if (rp%kstrength == 1) then
         ! Rothrock (1975): potential-energy change over the participation/redistribution functions
         call ridge_itd_cell(ncat, rp, puny, aice0, aicen, vicen, aksum, apartic, hrmin, hrmax, hrexp, krdg)
         ! potential-energy constant (g/2)(rhow-rhoi) rhoi/rhow, upstream operation order
         Cp = 0.5_wp*c%gravit*(c%rhow - c%rhoi)*c%rhoi/c%rhow
         ssum = 0.0_wp
         if (rp%krdg_redist == 0) then
            do n = 1, ncat
               if (aicen(n) > puny .and. apartic(n) > 0.0_wp) then
                  hi = vicen(n)/aicen(n)
                  h2rdg = P333*((hrmax(n)*hrmax(n))*hrmax(n) - (hrmin(n)*hrmin(n))*hrmin(n))/(hrmax(n) - hrmin(n))
                  dh2rdg = -hi*hi + h2rdg/krdg(n)
                  ssum = ssum + apartic(n)*dh2rdg
               end if
            end do
         else
            do n = 1, ncat
               if (aicen(n) > puny .and. apartic(n) > 0.0_wp) then
                  hi = vicen(n)/aicen(n)
                  h2rdg = hrmin(n)*hrmin(n) + 2.0_wp*hrmin(n)*hrexp(n) + 2.0_wp*hrexp(n)*hrexp(n)
                  dh2rdg = -hi*hi + h2rdg/krdg(n)
                  ssum = ssum + apartic(n)*dh2rdg
               end if
            end do
         end if
         if (aksum > 0.0_wp) then
            strength = rp%Cf*Cp*ssum/aksum
         else
            strength = 0.0_wp
            ierr = ior(ierr, ICE_ERR_RIDGE_DEGENERATE)
         end if
      else
         ! Hibler (1979)
         strength = rp%Pstar*vice*rdb_sm_exp(-(rp%Cstar*(1.0_wp - aice)))
      end if
   end subroutine ice_strength_cell

end module rdb_ice_ridge_itd
