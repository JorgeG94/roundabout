!! The Bitz-Lipscomb (1999) vertical heat solve: the layer conductivities, the
!! tridiagonal system rows for the snow/ice column, the Thomas solver and the
!! nonlinear iteration on the surface temperature that drives them
!! (`ice_temperature_changes`).
module rdb_ice_bl99
   ! SPDX-License-Identifier: BSD-3-Clause
   !
   ! Ported from CICE/Icepack (CICE Consortium, https://github.com/CICE-Consortium).
   ! Upstream: Icepack 9a1f3e4 (Icepack 1.5.3+14), module icepack_therm_bl99
   ! (conductivity, get_matrix_elements_calc_Tsfc, tridiag_solver).
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
   ! C. M. Bitz (UW), Elizabeth C. Hunke (LANL).
   ! Physics: Bitz and Lipscomb (1999, JGR 104 15669); Maykut and Untersteiner
   ! (1971, JGR 76 1550); Pringle et al. (2007, JGR 112 C04017).
   !
   !! **Column and rows.**  Rows are top-down: row 1 is the surface/snow dummy,
   !! rows `2..nslyr+1` the snow layers, rows `nslyr+2..nslyr+nilyr+1` the ice
   !! layers (`k = 1` the TOP ice layer, the documented exception to the
   !! ocean's bottom-up rule).  `kh(1:nmat)`, `nmat = nilyr + nslyr + 1`, are the
   !! interface conductances [W/m2/K]: `kh(1)` snow top, `kh(1+nslyr)` the
   !! snow/ice interface (or the ice top without snow), `kh(nmat)` the ice
   !! base.  They are evaluated once per step from the start-of-step
   !! temperatures and held fixed over the iteration.
   !!
   !! **Conductivity.**  `ICE_CONDUCT_MU71`:
   !! `kice + betak*S/min(-puny, T)` (Maykut and Untersteiner 1971 with the
   !! Wettlaufer 1991 constants, `betak = 0.13`); `ICE_CONDUCT_BUBBLY`
   !! (Pringle et al. 2007): `(2.11 - 0.011*T + 0.09*S/min(-puny, T))*rhoi/917`;
   !! both floored at `ICE_KIMIN = 0.10`.  Interface values are harmonic means
   !! `2*k1*k2/((k1+k2)*h)` of equal-thickness layers; the ends are `2*k/h`.
   !!
   !! **Matrix rows** (backward Euler, `etai = dt/(rhoi*ci*hilyr)`):
   !! `-etai*kh(k)*T(k-1) + (1 + etai*(kh(k)+kh(k+1)))*T(k) - etai*kh(k+1)*T(k+1)
   !! = Tinit + etai*Iswabs`, the bottom row adding `etai*kh*Tbot`.  A COLD
   !! surface adds the Newton-linearised surface balance as row 1 (or as the
   !! ice-top row without snow): `(dfsurf_dT - kh)*Tsf + kh*T1 = dfsurf_dT*Tsf_old
   !! - fsurfn`.  A MELTING surface (`l_cold = .false.`) pins `Tsf` and moves the
   !! `kh*Tsf` term to the right-hand side of the top layer, zeroing its
   !! sub-diagonal.  Properties the scheme guarantees (tested): the sum of a
   !! row's three coefficients is `1` for every layer row whose two interfaces
   !! both couple to unknown neighbours (a uniform temperature is a fixed point
   !! of the interior), `1 + etai*kh(nmat)` for the bottom ice row and
   !! `1 + etai*kh(top)` for a layer whose top neighbour is the pinned `Tsf`
   !! (the excess is exactly the coupling moved to the right-hand side), and
   !! `dfsurf_dT` for the cold surface row; the layer rows are strictly
   !! diagonally dominant, so the Thomas solve needs no pivoting.
   !! The `nslyr > 1` branches (interior snow interfaces and rows) are kept as
   !! upstream wrote them but cannot run while `ICE_NSLYR_MAX = 1`, so they are
   !! untested.
   !!
   !! The `calc_Tsfc = .false.` (surface temperature known) builder is not
   !! ported (the configuration refuses it).
   !!
   !! **The iteration** (`ice_temperature_changes`, `calc_Tsfc = .true.`).
   !! Conductances are evaluated once from the start-of-step temperatures.
   !! Each pass (at most `nitermax = 100`) updates the brine heat capacity
   !! `ci = cp_ice - Lfresh*Tmlt/(T*T_init)` from the latest iterate and the
   !! start value (the energy-conserving trick of BL99, their eq. 6-7), evaluates
   !! the bulk surface flux and its derivative at `Tsf`, the conductive flux
   !! `fcondtopn = kh_top*(Tsf - T_top)`, resets `Tsf = -puny` when a
   !! non-negative `Tsf` would need `fsurfn < fcondtopn` (the surface wants to
   !! cool, so the cold branch is re-entered), builds the rows of the COLD
   !! (`Tsf < 0`, an unknown) or MELTING (`Tsf = 0`, a pinned boundary value)
   !! branch and solves.  The result is then limited and tested against five
   !! conditions: (1) `Tsf > puny` is reset to 0 (and the interior averaged
   !! with its start values), (2) an oscillating `Tsf` (`-dTsf/dTsf_prev > 0.5`)
   !! halves `dTsf` and averages, (3) `|dTsf| <= 5e-4`, (4) `Tsf >= 0` needs
   !! `fsurfn >= fcondtopn`, (5) the energy error `ferr <= 0.9*ferrmax`
   !! (`ferrmax = 1e-3 W/m2`; with brine the scheme conserves only to this bound,
   !! not to round-off).  With brine, layer and snow temperatures are limited at
   !! their melting temperatures (snow at 0), the excess enthalpy `dqmat`
   !! leaves through the bottom flux (`fcondbot += einex_sfc_calc/dt`) and the
   !! conductance of limited layers is reduced when condition (5) fails.
   !! After convergence the fluxes are corrected for the last `Tsf` update
   !! (`flwoutn += dTsf_prev*dflwout_dT`, likewise sensible and latent).
   !! `sw_redist` first moves shortwave absorbed in a layer beyond `sw_frac`
   !! of what warms it to melting up into `fswsfc` (bounded by `fswint`).
   !! Not converged after `nitermax` passes sets `ICE_ERR_TSFC_NOCONV` and
   !! returns the last iterate (it still applies the post-loop flux correction
   !! that upstream skips by aborting first; NOCONV is fatal, so that is not
   !! observable).  `zqin`, `zqsn` and the flux outputs are written, never read.
   !! Non-finite INPUT sets `ICE_ERR_NONFINITE`
   !! before anything is written; a non-finite intermediate sets it and
   !! returns with the outputs unspecified (the caller does not commit a
   !! column with a fatal bit).  The solver works on its in/out arguments in
   !! place (no buffered copy of the column); constants are validated once at
   !! configure, not per call.
   !! Ported scope: `calc_Tsfc = .true.`, `semi_implicit_Tsfc = .false.`;
   !! the known-`Tsf` builder and its flux-reduction bookkeeping are dropped.
   !!
   !! **Kernel contract.**  Every procedure is `pure` with `!$acc routine seq`,
   !! explicit-shape dummies sized by the layer counts (callers pass arrays of
   !! at least that size, typically `ICE_NMAT_MAX`), no allocation.  Constants
   !! arrive in the `ice_constants_t` bundle.  `ierr` is only OR-ed into: a
   !! non-finite input, a layer count outside `1..ICE_NILYR_MAX`
   !! (`1..ICE_NSLYR_MAX` for snow), a non-positive layer thickness where it is
   !! divided by, or a singular pivot sets `ICE_ERR_NONFINITE` (the registry has
   !! no better-fitting bit) and the output arrays are left untouched.  The
   !! entry checks are what stop NaN from being laundered by `min`/comparisons.
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use, intrinsic :: iso_fortran_env, only: int32
   use rdb_constants, only: wp
   use rdb_ice_dims, only: ICE_NILYR_MAX, ICE_NSLYR_MAX
   use rdb_ice_errors, only: ICE_ERR_NONFINITE, ICE_ERR_TSFC_NOCONV
   use rdb_ice_params, only: ice_constants_t, ICE_CONDUCT_MU71
   use rdb_ice_thermo_shared, only: surface_heat_flux, dsurface_heat_flux_dTsf
   implicit none
   private

   public :: ICE_NMAT_MAX, ICE_KIMIN, ICE_BETAK
   public :: ice_layer_conductivity, ice_conductivity
   public :: ice_matrix_calc_Tsfc, ice_tridiag_solver
   public :: ice_conductivity_core, ice_matrix_calc_Tsfc_core, ice_tridiag_solver_core
   public :: ice_temperature_changes, ice_temperature_changes_cell
   public :: ICE_NITERMAX, ICE_TSF_ERRMAX, ICE_FERRMAX, ICE_SW_FRAC, ICE_SW_DTEMP

   integer, parameter :: ICE_NMAT_MAX = ICE_NILYR_MAX + ICE_NSLYR_MAX + 1
      !! Largest matrix dimension (`nilyr + nslyr + 1`).
   real(wp), parameter :: ICE_KIMIN = 0.10_wp
      !! Minimum conductivity of saline ice [W/m/K].
   real(wp), parameter :: ICE_BETAK = 0.13_wp
      !! Salinity coefficient of the MU71 conductivity [W/m/K/ppt].
   real(wp), parameter :: ICE_K_BUBBLY_RHOREF = 917.0_wp
      !! Reference density of the bubbly-brine fit [kg/m3].
   integer, parameter :: ICE_NITERMAX = 100
      !! Iteration cap of the temperature solve (upstream `nitermax`).
   real(wp), parameter :: ICE_TSF_ERRMAX = 5.0e-4_wp
      !! Largest accepted `|dTsf|` of the last pass [K] (upstream `Tsf_errmax`).
   real(wp), parameter :: ICE_FERRMAX = 1.0e-3_wp
      !! Energy-error bound [W/m2]; passes are accepted at `0.9*ferrmax`.
   real(wp), parameter :: ICE_SW_FRAC = 0.9_wp
      !! Default fraction of the warm-to-melting energy a layer may absorb
      !! before `sw_redist` moves the excess shortwave to the surface.
   real(wp), parameter :: ICE_SW_DTEMP = 0.02_wp
      !! Default `sw_redist` temperature offset from melting [K].

contains

   pure function ice_layer_conductivity(T, S, c, conduct) result(k)
      !! Conductivity [W/m/K] of an ice layer at temperature `T` [degC] and
      !! salinity `S` [ppt] for the `conduct` option, floored at `ICE_KIMIN`.
      !! `T` must be finite (the guard `min(-puny, T)` is NaN-blind); a value
      !! other than `ICE_CONDUCT_MU71` selects the bubbly form, as upstream.
      !$acc routine seq
      real(wp), intent(in) :: T, S
      type(ice_constants_t), intent(in) :: c
      integer, intent(in) :: conduct
      real(wp) :: k

      if (conduct == ICE_CONDUCT_MU71) then
         k = c%kice + ICE_BETAK*S/min(-c%puny, T)
      else
         k = (2.11_wp - 0.011_wp*T + 0.09_wp*S/min(-c%puny, T))*c%rhoi/ICE_K_BUBBLY_RHOREF
      end if
      k = max(k, ICE_KIMIN)
   end function ice_layer_conductivity

   pure subroutine ice_conductivity(l_snow, nilyr, nslyr, hilyr, hslyr, zTin, zSin, c, conduct, kh, ierr)
      !! Interface conductances `kh(1:nilyr+nslyr+1)` of the snow/ice column
      !! (see the module header).  Without snow (`l_snow = .false.`) the snow
      !! interfaces are zero and `kh(1+nslyr)` is `2*k1/hilyr`.  Requires
      !! `hilyr > 0` and, with snow, `hslyr > 0`.  An entry-check failure
      !! leaves `kh` untouched; a non-finite result is only found after `kh`
      !! was written, so on any later error `kh` is unspecified.
      !$acc routine seq
      logical, intent(in) :: l_snow
      integer, intent(in) :: nilyr, nslyr, conduct
      real(wp), intent(in) :: hilyr, hslyr
      real(wp), intent(in) :: zTin(nilyr), zSin(nilyr)
      type(ice_constants_t), intent(in) :: c
      real(wp), intent(inout) :: kh(nilyr + nslyr + 1)
      integer(int32), intent(inout) :: ierr
      integer :: k

      if (nilyr < 1 .or. nilyr > ICE_NILYR_MAX .or. nslyr < 1 .or. nslyr > ICE_NSLYR_MAX) then
         ierr = ior(ierr, ICE_ERR_NONFINITE)
         return
      end if
      if (.not. (hilyr > 0.0_wp .and. ieee_is_finite(hilyr) .and. ieee_is_finite(hslyr) &
                 .and. ieee_is_finite(c%puny) .and. ieee_is_finite(c%kice) .and. ieee_is_finite(c%ksno) &
                 .and. ieee_is_finite(c%rhoi))) then
         ierr = ior(ierr, ICE_ERR_NONFINITE)
         return
      end if
      if (l_snow .and. .not. (hslyr > 0.0_wp)) then
         ierr = ior(ierr, ICE_ERR_NONFINITE)
         return
      end if
      do k = 1, nilyr
         if (.not. (ieee_is_finite(zTin(k)) .and. ieee_is_finite(zSin(k)))) then
            ierr = ior(ierr, ICE_ERR_NONFINITE)
            return
         end if
      end do

      call ice_conductivity_core(l_snow, nilyr, nslyr, hilyr, hslyr, zTin, zSin, c, conduct, kh)
      do k = 1, nilyr + nslyr + 1
         if (.not. ieee_is_finite(kh(k))) then
            ierr = ior(ierr, ICE_ERR_NONFINITE)
            return
         end if
      end do
   end subroutine ice_conductivity

   pure subroutine ice_conductivity_core(l_snow, nilyr, nslyr, hilyr, hslyr, zTin, zSin, c, conduct, kh)
      !! `ice_conductivity` without any check: the caller has validated the
      !! layer counts, `hilyr > 0`, `hslyr` and the temperatures / salinities
      !! (the temperature solve does it once at entry).  A non-finite result
      !! is caught downstream, by the pivot and the solution test of
      !! `ice_tridiag_solver_core`.
      !$acc routine seq
      logical, intent(in) :: l_snow
      integer, intent(in) :: nilyr, nslyr, conduct
      real(wp), intent(in) :: hilyr, hslyr
      real(wp), intent(in) :: zTin(nilyr), zSin(nilyr)
      type(ice_constants_t), intent(in) :: c
      real(wp), intent(inout) :: kh(nilyr + nslyr + 1)
      integer :: k
      real(wp) :: kslyr
      real(wp) :: kilyr(ICE_NILYR_MAX)

      kslyr = c%ksno
      do k = 1, nilyr
         kilyr(k) = ice_layer_conductivity(zTin(k), zSin(k), c, conduct)
      end do

      if (l_snow) then
         kh(1) = 2.0_wp*kslyr/hslyr
         kh(1 + nslyr) = 2.0_wp*kslyr*kilyr(1)/(kslyr*hilyr + kilyr(1)*hslyr)
      else
         kh(1) = 0.0_wp
         kh(1 + nslyr) = 2.0_wp*kilyr(1)/hilyr
      end if
      kh(1 + nslyr + nilyr) = 2.0_wp*kilyr(nilyr)/hilyr

      if (nslyr > 1) then
         do k = 2, nslyr
            if (l_snow) then
               kh(k) = 2.0_wp*kslyr*kslyr/((kslyr + kslyr)*hslyr)
            else
               kh(k) = 0.0_wp
            end if
         end do
      end if

      do k = 2, nilyr
         kh(k + nslyr) = 2.0_wp*kilyr(k - 1)*kilyr(k)/((kilyr(k - 1) + kilyr(k))*hilyr)
      end do
   end subroutine ice_conductivity_core

   pure subroutine ice_matrix_calc_Tsfc(l_snow, l_cold, nilyr, nslyr, Tsf, Tbot, fsurfn, dfsurf_dT, &
                                        Tin_init, Tsn_init, kh, Sswabs, Iswabs, etai, etas, &
                                        sbdiag, diag, spdiag, rhs, ierr)
      !! Rows of the tridiagonal system `sbdiag*T(k-1) + diag*T(k) +
      !! spdiag*T(k+1) = rhs` for the column when the surface temperature is
      !! solved for.  Four cases: cold/melting surface, with/without snow.
      !! Rows that carry no unknown (the snow/surface dummies) solve `T = 0`.
      !! `Tsf` is the previous iterate (cold branch) or the pinned melting
      !! value (melting branch); `fsurfn`, `dfsurf_dT` the net surface flux and
      !! its derivative; `Sswabs`/`Iswabs` absorbed shortwave [W/m2];
      !! `etas`/`etai` = `dt/(rho*cp*h)` of the snow/ice layers.
      !$acc routine seq
      logical, intent(in) :: l_snow, l_cold
      integer, intent(in) :: nilyr, nslyr
      real(wp), intent(in) :: Tsf, Tbot, fsurfn, dfsurf_dT
      real(wp), intent(in) :: Tin_init(nilyr), Tsn_init(nslyr), Sswabs(nslyr), Iswabs(nilyr)
      real(wp), intent(in) :: etai(nilyr), etas(nslyr), kh(nilyr + nslyr + 1)
      real(wp), intent(inout) :: sbdiag(nilyr + nslyr + 1), diag(nilyr + nslyr + 1)
      real(wp), intent(inout) :: spdiag(nilyr + nslyr + 1), rhs(nilyr + nslyr + 1)
      integer(int32), intent(inout) :: ierr
      integer :: k, nmat

      if (nilyr < 1 .or. nilyr > ICE_NILYR_MAX .or. nslyr < 1 .or. nslyr > ICE_NSLYR_MAX) then
         ierr = ior(ierr, ICE_ERR_NONFINITE)
         return
      end if
      nmat = nilyr + nslyr + 1
      if (.not. (ieee_is_finite(Tsf) .and. ieee_is_finite(Tbot) .and. ieee_is_finite(fsurfn) &
                 .and. ieee_is_finite(dfsurf_dT))) then
         ierr = ior(ierr, ICE_ERR_NONFINITE)
         return
      end if
      do k = 1, nilyr
         if (.not. (ieee_is_finite(Tin_init(k)) .and. ieee_is_finite(Iswabs(k)) .and. ieee_is_finite(etai(k)))) then
            ierr = ior(ierr, ICE_ERR_NONFINITE)
            return
         end if
      end do
      do k = 1, nslyr
         if (.not. (ieee_is_finite(Tsn_init(k)) .and. ieee_is_finite(Sswabs(k)) .and. ieee_is_finite(etas(k)))) then
            ierr = ior(ierr, ICE_ERR_NONFINITE)
            return
         end if
      end do
      do k = 1, nmat
         if (.not. ieee_is_finite(kh(k))) then
            ierr = ior(ierr, ICE_ERR_NONFINITE)
            return
         end if
      end do

      call ice_matrix_calc_Tsfc_core(l_snow, l_cold, nilyr, nslyr, Tsf, Tbot, fsurfn, dfsurf_dT, &
                                     Tin_init, Tsn_init, kh, Sswabs, Iswabs, etai, etas, &
                                     sbdiag, diag, spdiag, rhs)
   end subroutine ice_matrix_calc_Tsfc

   pure subroutine ice_matrix_calc_Tsfc_core(l_snow, l_cold, nilyr, nslyr, Tsf, Tbot, fsurfn, dfsurf_dT, &
                                             Tin_init, Tsn_init, kh, Sswabs, Iswabs, etai, etas, &
                                             sbdiag, diag, spdiag, rhs)
      !! `ice_matrix_calc_Tsfc` without any check: the caller has validated the
      !! layer counts and the finiteness of the inputs (or relies on the
      !! pivot / solution test of `ice_tridiag_solver_core`, which a
      !! non-finite coefficient cannot pass).
      !$acc routine seq
      logical, intent(in) :: l_snow, l_cold
      integer, intent(in) :: nilyr, nslyr
      real(wp), intent(in) :: Tsf, Tbot, fsurfn, dfsurf_dT
      real(wp), intent(in) :: Tin_init(nilyr), Tsn_init(nslyr), Sswabs(nslyr), Iswabs(nilyr)
      real(wp), intent(in) :: etai(nilyr), etas(nslyr), kh(nilyr + nslyr + 1)
      real(wp), intent(inout) :: sbdiag(nilyr + nslyr + 1), diag(nilyr + nslyr + 1)
      real(wp), intent(inout) :: spdiag(nilyr + nslyr + 1), rhs(nilyr + nslyr + 1)
      integer :: k, ki, kr

      ! Dummy rows (solution T = 0) for the surface and snow rows.
      do k = 1, nslyr + 1
         sbdiag(k) = 0.0_wp
         diag(k) = 1.0_wp
         spdiag(k) = 0.0_wp
         rhs(k) = 0.0_wp
      end do

      ! Surface equation, cold surface (with or without snow).
      if (l_cold) then
         if (l_snow) then
            k = 1
         else
            k = 1 + nslyr
         end if
         kr = k
         sbdiag(kr) = 0.0_wp
         diag(kr) = dfsurf_dT - kh(k)
         spdiag(kr) = kh(k)
         rhs(kr) = dfsurf_dT*Tsf - fsurfn
      end if

      ! Top snow layer.
      if (l_snow) then
         if (l_cold) then
            sbdiag(2) = -etas(1)*kh(1)
            spdiag(2) = -etas(1)*kh(2)
            diag(2) = 1.0_wp + etas(1)*(kh(1) + kh(2))
            rhs(2) = Tsn_init(1) + etas(1)*Sswabs(1)
         else
            sbdiag(2) = 0.0_wp
            spdiag(2) = -etas(1)*kh(2)
            diag(2) = 1.0_wp + etas(1)*(kh(1) + kh(2))
            rhs(2) = Tsn_init(1) + etas(1)*kh(1)*Tsf + etas(1)*Sswabs(1)
         end if
      end if

      ! Remaining snow layers.
      if (nslyr > 1) then
         do k = 2, nslyr
            kr = k + 1
            if (l_snow) then
               sbdiag(kr) = -etas(k)*kh(k)
               spdiag(kr) = -etas(k)*kh(k + 1)
               diag(kr) = 1.0_wp + etas(k)*(kh(k) + kh(k + 1))
               rhs(kr) = Tsn_init(k) + etas(k)*Sswabs(k)
            end if
         end do
      end if

      if (nilyr > 1) then
         ! Top ice layer.
         ki = 1
         k = ki + nslyr
         kr = k + 1
         if (l_snow .or. l_cold) then
            sbdiag(kr) = -etai(ki)*kh(k)
            spdiag(kr) = -etai(ki)*kh(k + 1)
            diag(kr) = 1.0_wp + etai(ki)*(kh(k) + kh(k + 1))
            rhs(kr) = Tin_init(ki) + etai(ki)*Iswabs(ki)
         else
            sbdiag(kr) = 0.0_wp
            spdiag(kr) = -etai(ki)*kh(k + 1)
            diag(kr) = 1.0_wp + etai(ki)*(kh(k) + kh(k + 1))
            rhs(kr) = Tin_init(ki) + etai(ki)*Iswabs(ki) + etai(ki)*kh(k)*Tsf
         end if

         ! Bottom ice layer.
         ki = nilyr
         k = ki + nslyr
         kr = k + 1
         sbdiag(kr) = -etai(ki)*kh(k)
         spdiag(kr) = 0.0_wp
         diag(kr) = 1.0_wp + etai(ki)*(kh(k) + kh(k + 1))
         rhs(kr) = Tin_init(ki) + etai(ki)*Iswabs(ki) + etai(ki)*kh(k + 1)*Tbot
      else
         ! Single ice layer.
         ki = 1
         k = ki + nslyr
         kr = k + 1
         if (l_snow .or. l_cold) then
            sbdiag(kr) = -etai(ki)*kh(k)
            spdiag(kr) = 0.0_wp
            diag(kr) = 1.0_wp + etai(ki)*(kh(k) + kh(k + 1))
            rhs(kr) = Tin_init(ki) + etai(ki)*Iswabs(ki) + etai(ki)*kh(k + 1)*Tbot
         else
            sbdiag(kr) = 0.0_wp
            spdiag(kr) = 0.0_wp
            diag(kr) = 1.0_wp + etai(ki)*(kh(k) + kh(k + 1))
            rhs(kr) = Tin_init(ki) + etai(ki)*Iswabs(ki) + etai(ki)*kh(k)*Tsf + etai(ki)*kh(k + 1)*Tbot
         end if
      end if

      ! Interior ice layers.
      do ki = 2, nilyr - 1
         k = ki + nslyr
         kr = k + 1
         sbdiag(kr) = -etai(ki)*kh(k)
         spdiag(kr) = -etai(ki)*kh(k + 1)
         diag(kr) = 1.0_wp + etai(ki)*(kh(k) + kh(k + 1))
         rhs(kr) = Tin_init(ki) + etai(ki)*Iswabs(ki)
      end do

   end subroutine ice_matrix_calc_Tsfc_core

   pure subroutine ice_tridiag_solver(nmat, sbdiag, diag, spdiag, rhs, xout, ierr)
      !! Thomas algorithm for `sbdiag(k)*xout(k-1) + diag(k)*xout(k) +
      !! spdiag(k)*xout(k+1) = rhs(k)`, `k = 1..nmat` (`sbdiag(1)` and
      !! `spdiag(nmat)` unused), no pivoting: the BL99 rows are diagonally
      !! dominant.  A non-finite coefficient or a zero/non-finite pivot sets
      !! `ICE_ERR_NONFINITE` and is returned with `xout` unspecified (non-finite
      !! coefficients are caught before `xout` is written).
      !$acc routine seq
      integer, intent(in) :: nmat
      real(wp), intent(in) :: sbdiag(nmat), diag(nmat), spdiag(nmat), rhs(nmat)
      real(wp), intent(inout) :: xout(nmat)
      integer(int32), intent(inout) :: ierr
      integer :: k

      if (nmat < 1 .or. nmat > ICE_NMAT_MAX) then
         ierr = ior(ierr, ICE_ERR_NONFINITE)
         return
      end if
      do k = 1, nmat
         if (.not. (ieee_is_finite(diag(k)) .and. ieee_is_finite(rhs(k)) .and. ieee_is_finite(sbdiag(k)) &
                    .and. ieee_is_finite(spdiag(k)))) then
            ierr = ior(ierr, ICE_ERR_NONFINITE)
            return
         end if
      end do

      call ice_tridiag_solver_core(nmat, sbdiag, diag, spdiag, rhs, xout, ierr)
   end subroutine ice_tridiag_solver

   pure subroutine ice_tridiag_solver_core(nmat, sbdiag, diag, spdiag, rhs, xout, ierr)
      !! `ice_tridiag_solver` without the entry scan of the coefficients: the
      !! caller has validated `1 <= nmat <= ICE_NMAT_MAX`; `sbdiag(1)` and
      !! `spdiag(nmat)` are never read.  A non-finite
      !! coefficient cannot pass unnoticed: it makes a pivot non-finite, or
      !! (`rhs`, `spdiag`) a solution entry, and both are tested, with
      !! `ieee_is_finite`, so the same `ICE_ERR_NONFINITE` is raised.  `xout`
      !! is written as the sweep goes, so it is unspecified on error.
      !$acc routine seq
      integer, intent(in) :: nmat
      real(wp), intent(in) :: sbdiag(nmat), diag(nmat), spdiag(nmat), rhs(nmat)
      real(wp), intent(inout) :: xout(nmat)
      integer(int32), intent(inout) :: ierr
      integer :: k
      real(wp) :: wbeta
      real(wp) :: wgamma(ICE_NMAT_MAX)

      wbeta = diag(1)
      if (.not. (ieee_is_finite(wbeta) .and. abs(wbeta) > tiny(1.0_wp))) then
         ierr = ior(ierr, ICE_ERR_NONFINITE)
         return
      end if
      xout(1) = rhs(1)/wbeta
      wgamma(1) = 0.0_wp
      do k = 2, nmat
         wgamma(k) = spdiag(k - 1)/wbeta
         wbeta = diag(k) - sbdiag(k)*wgamma(k)
         if (.not. (ieee_is_finite(wbeta) .and. abs(wbeta) > tiny(1.0_wp))) then
            ierr = ior(ierr, ICE_ERR_NONFINITE)
            return
         end if
         xout(k) = (rhs(k) - sbdiag(k)*xout(k - 1))/wbeta
      end do
      do k = nmat - 1, 1, -1
         xout(k) = xout(k) - wgamma(k + 1)*xout(k + 1)
      end do

      do k = 1, nmat
         if (.not. ieee_is_finite(xout(k))) then
            ierr = ior(ierr, ICE_ERR_NONFINITE)
            return
         end if
      end do
   end subroutine ice_tridiag_solver_core

   pure subroutine ice_temperature_changes(nilyr, nslyr, dt, rhoa, flw, potT, Qa, shcoef, lhcoef, &
                                           emissivity, sfc_kind, sf_0, dsf_dT, c, conduct, l_brine, &
                                           sw_redist, sw_frac, sw_dtemp, nitermax, &
                                           hilyr, hslyr, zSin, Tbot, einit, &
                                           fswsfc, fswint, Sswabs, Iswabs, zqin, zTin, zqsn, zTsn, Tsf, &
                                           fsensn, flatn, flwoutn, fsurfn, fcondtopn, fcondbot, &
                                           einex_sfc_flux, diag, ierr)
      !! New snow/ice temperatures, surface temperature and the conductive
      !! fluxes of one column over `dt` (see the module header).  Layers are
      !! top-down.  Inputs: atmosphere `rhoa, flw, potT, Qa, shcoef, lhcoef`,
      !! `emissivity`, the surface balance selector `sfc_kind` (`sf_0, dsf_dT`
      !! only for the linear test surface), `conduct`, `l_brine`, the
      !! `sw_redist` switch and its `sw_frac, sw_dtemp`, `nitermax`, the layer
      !! thicknesses, salinities, `Tbot` (ice base temperature) and `einit`
      !! (initial column enthalpy, J/m2).  In/out: absorbed shortwave `fswsfc,
      !! fswint, Sswabs, Iswabs` (changed only by `sw_redist`), enthalpies
      !! `zqin, zqsn` [J/m3], temperatures `zTin, zTsn` and `Tsf` (the start
      !! value on entry).  Outputs: `fsensn, flatn, flwoutn` (downward
      !! positive, as the surface flux), `fsurfn` (net surface flux excluding
      !! the conduction), `fcondtopn, fcondbot` [W/m2], `einex_sfc_flux`
      !! (zero for `calc_Tsfc`).  `diag = (iterations, ferr, dTsf_prev)`.
      !$acc routine seq
      integer, intent(in) :: nilyr, nslyr
      real(wp), intent(in) :: dt, rhoa, flw, potT, Qa, shcoef, lhcoef
      real(wp), intent(in) :: emissivity, sf_0, dsf_dT
      integer, intent(in) :: sfc_kind, conduct, nitermax
      type(ice_constants_t), intent(in) :: c
      logical, intent(in) :: l_brine, sw_redist
      real(wp), intent(in) :: sw_frac, sw_dtemp
      real(wp), intent(in) :: hilyr, hslyr, Tbot, einit
      real(wp), intent(in) :: zSin(nilyr)
      real(wp), intent(inout) :: fswsfc, fswint
      real(wp), intent(inout) :: Sswabs(nslyr), Iswabs(nilyr)
      real(wp), intent(inout) :: zqin(nilyr), zTin(nilyr), zqsn(nslyr), zTsn(nslyr)
      real(wp), intent(inout) :: Tsf
      real(wp), intent(inout) :: fsensn, flatn, flwoutn, fsurfn, fcondtopn, fcondbot, einex_sfc_flux
      real(wp), intent(inout) :: diag(3)
      integer(int32), intent(inout) :: ierr

      real(wp), parameter :: P5 = 0.5_wp
      integer :: k, niter, nmat, niter_done
      logical :: l_snow, l_cold, converged
      real(wp) :: Tsf_start, dTsf, dTsf_prev, dfsurf_dT, avg_Tsi, avg_Tsf, enew
      real(wp) :: dfsens_dT, dflat_dT, dflwout_dT, dt_rhoi_hlyr, einex_sfc_calc, ferr
      real(wp) :: ci, Iswabs_tmp, Sswabs_tmp, dswabs, frac, dTmat
      real(wp) :: Tin_init(ICE_NILYR_MAX), Tin_start(ICE_NILYR_MAX), dqmat(ICE_NILYR_MAX)
      real(wp) :: Tmlts(ICE_NILYR_MAX), etai(ICE_NILYR_MAX), LfTm(ICE_NILYR_MAX)
      real(wp) :: Tsn_init(ICE_NSLYR_MAX), Tsn_start(ICE_NSLYR_MAX), etas(ICE_NSLYR_MAX)
      real(wp) :: sbdiag(ICE_NMAT_MAX), diagm(ICE_NMAT_MAX), spdiag(ICE_NMAT_MAX)
      real(wp) :: rhs(ICE_NMAT_MAX), Tmat(ICE_NMAT_MAX), kh(ICE_NMAT_MAX)
      logical :: reduce_kh(ICE_NILYR_MAX)
      integer(int32) :: ierl

      ! ---- entry checks: a non-finite or out-of-range input flags and returns
      if (nilyr < 1 .or. nilyr > ICE_NILYR_MAX .or. nslyr < 1 .or. nslyr > ICE_NSLYR_MAX) then
         ierr = ior(ierr, ICE_ERR_NONFINITE)
         return
      end if
      if (.not. (ieee_is_finite(dt) .and. ieee_is_finite(rhoa) .and. ieee_is_finite(flw) &
                 .and. ieee_is_finite(potT) .and. ieee_is_finite(Qa) .and. ieee_is_finite(shcoef) &
                 .and. ieee_is_finite(lhcoef) .and. ieee_is_finite(emissivity) .and. ieee_is_finite(sf_0) &
                 .and. ieee_is_finite(dsf_dT) .and. ieee_is_finite(sw_frac) .and. ieee_is_finite(sw_dtemp) &
                 .and. ieee_is_finite(hilyr) .and. ieee_is_finite(hslyr) .and. ieee_is_finite(Tbot) &
                 .and. ieee_is_finite(einit) .and. ieee_is_finite(fswsfc) .and. ieee_is_finite(fswint) &
                 .and. ieee_is_finite(Tsf))) then
         ierr = ior(ierr, ICE_ERR_NONFINITE)
         return
      end if
      if (.not. (dt > 0.0_wp .and. hilyr > 0.0_wp .and. hslyr >= 0.0_wp .and. nitermax >= 1)) then
         ierr = ior(ierr, ICE_ERR_NONFINITE)
         return
      end if
      do k = 1, nilyr
         if (.not. (ieee_is_finite(zSin(k)) .and. ieee_is_finite(zTin(k)) &
                    .and. ieee_is_finite(Iswabs(k)))) then
            ierr = ior(ierr, ICE_ERR_NONFINITE)
            return
         end if
      end do
      do k = 1, nslyr
         if (.not. (ieee_is_finite(zTsn(k)) .and. ieee_is_finite(Sswabs(k)))) then
            ierr = ior(ierr, ICE_ERR_NONFINITE)
            return
         end if
      end do

      ! ---- initialise.  The temperature, enthalpy and shortwave arguments are
      !      overwritten in place at every pass; after an error return they hold
      !      the last iterate and must not be used (the cell entry gathers the
      !      column into local buffers and commits them only on success).
      ierl = 0_int32
      nmat = nilyr + nslyr + 1
      einex_sfc_calc = 0.0_wp
      dTsf = 0.0_wp
      dTsf_prev = 0.0_wp
      dfsurf_dT = 0.0_wp
      dfsens_dT = 0.0_wp
      dflat_dT = 0.0_wp
      dflwout_dT = 0.0_wp
      fcondbot = 0.0_wp
      fcondtopn = 0.0_wp
      fsurfn = 0.0_wp
      flwoutn = 0.0_wp
      fsensn = 0.0_wp
      flatn = 0.0_wp
      ferr = 0.0_wp
      enew = 0.0_wp
      l_cold = .true.
      l_snow = hslyr > c%hs_min/real(nslyr, wp)
      dt_rhoi_hlyr = dt/(c%rhoi*hilyr)

      do k = 1, nslyr
         Tsn_init(k) = zTsn(k)
         Tsn_start(k) = zTsn(k)
         if (l_snow) then
            etas(k) = dt/(c%rhos*c%cp_ice*hslyr)
         else
            etas(k) = 0.0_wp
         end if
      end do
      do k = 1, nilyr
         Tin_init(k) = zTin(k)
         Tin_start(k) = zTin(k)
         Tmlts(k) = -zSin(k)*c%depressT
         LfTm(k) = c%Lfresh*Tmlts(k)
      end do

      call ice_conductivity_core(l_snow, nilyr, nslyr, hilyr, hslyr, zTin, zSin, c, conduct, kh)

      ! ---- shortwave redistribution: energy beyond what warms a layer to
      !      near melting goes to the surface instead
      if (sw_redist) then
         do k = 1, nilyr
            Iswabs_tmp = 0.0_wp
            if (Tin_init(k) <= Tmlts(k) - sw_dtemp) then
               if (l_brine) then
                  ci = c%cp_ice - c%Lfresh*Tmlts(k)/(Tin_init(k)**2)
                  Iswabs_tmp = min(Iswabs(k), sw_frac*(Tmlts(k) - Tin_init(k))*ci/dt_rhoi_hlyr)
               else
                  ci = c%cp_ice
                  Iswabs_tmp = min(Iswabs(k), sw_frac*(-Tin_init(k))*ci/dt_rhoi_hlyr)
               end if
            end if
            if (Iswabs_tmp < c%puny) Iswabs_tmp = 0.0_wp

            dswabs = min(Iswabs(k) - Iswabs_tmp, fswint)
            fswsfc = fswsfc + dswabs
            fswint = fswint - dswabs
            Iswabs(k) = Iswabs_tmp
         end do

         do k = 1, nslyr
            if (l_snow) then
               Sswabs_tmp = 0.0_wp
               if (Tsn_init(k) <= -sw_dtemp) then
                  Sswabs_tmp = min(Sswabs(k), -sw_frac*Tsn_init(k)/etas(k))
               end if
               if (Sswabs_tmp < c%puny) Sswabs_tmp = 0.0_wp

               dswabs = min(Sswabs(k) - Sswabs_tmp, fswint)
               fswsfc = fswsfc + dswabs
               fswint = fswint - dswabs
               Sswabs(k) = Sswabs_tmp
            end if
         end do
      end if

      ! ---- iterate on the temperatures
      converged = .false.
      niter_done = 0
      do niter = 1, nitermax
         niter_done = niter
         converged = .true.
         dfsurf_dT = 0.0_wp
         avg_Tsi = 0.0_wp
         enew = 0.0_wp
         einex_sfc_calc = 0.0_wp

         ! heat capacity from the latest iterate and the start value
         do k = 1, nilyr
            if (l_brine) then
               ci = c%cp_ice - LfTm(k)/(zTin(k)*Tin_init(k))
            else
               ci = c%cp_ice
            end if
            etai(k) = dt_rhoi_hlyr/ci
         end do

         ! surface flux and derivative at the current Tsf
         call surface_heat_flux(Tsf, fswsfc, rhoa, flw, potT, Qa, shcoef, lhcoef, emissivity, c, &
                                sfc_kind, sf_0, dsf_dT, flwoutn, fsensn, flatn, fsurfn, ierl)
         call dsurface_heat_flux_dTsf(Tsf, rhoa, shcoef, lhcoef, emissivity, c, sfc_kind, dsf_dT, &
                                      dfsurf_dT, dflwout_dT, dfsens_dT, dflat_dT, ierl)
         if (ierl /= 0_int32) then
            ierr = ior(ierr, ierl)
            return
         end if

         ! conductive flux at the top; a surface that wants to cool re-enters
         ! the cold branch
         if (l_snow) then
            fcondtopn = kh(1)*(Tsf - zTsn(1))
         else
            fcondtopn = kh(1 + nslyr)*(Tsf - zTin(1))
         end if
         if (Tsf >= 0.0_wp .and. fsurfn < fcondtopn) Tsf = -c%puny

         Tsf_start = Tsf
         l_cold = Tsf < 0.0_wp

         call ice_matrix_calc_Tsfc_core(l_snow, l_cold, nilyr, nslyr, Tsf, Tbot, fsurfn, dfsurf_dT, &
                                        Tin_init, Tsn_init, kh, Sswabs, Iswabs, etai, etas, &
                                        sbdiag, diagm, spdiag, rhs)
         call ice_tridiag_solver_core(nmat, sbdiag, diagm, spdiag, rhs, Tmat, ierl)
         if (ierl /= 0_int32) then
            ierr = ior(ierr, ierl)
            return
         end if

         ! reload Tsf from the solution
         if (l_cold) then
            if (l_snow) then
               Tsf = Tmat(1)
            else
               Tsf = Tmat(1 + nslyr)
            end if
         else
            Tsf = 0.0_wp
         end if

         dTsf = Tsf - Tsf_start
         avg_Tsf = 0.0_wp

         ! condition 1: Tsf > 0; condition 2: oscillating Tsf
         if (Tsf > c%puny) then
            Tsf = 0.0_wp
            dTsf = -Tsf_start
            if (l_brine) avg_Tsi = 1.0_wp
            converged = .false.
         else if (niter > 1 .and. Tsf_start <= -c%puny .and. abs(dTsf) > c%puny &
                  .and. abs(dTsf_prev) > c%puny .and. -dTsf/(dTsf_prev + c%puny*c%puny) > P5) then
            if (l_brine) then
               avg_Tsf = 1.0_wp
               avg_Tsi = 1.0_wp
            end if
            dTsf = P5*dTsf
            converged = .false.
         end if

         Tsf = Tsf + avg_Tsf*P5*(Tsf_start - Tsf)

         do k = 1, nslyr
            if (l_snow) then
               zTsn(k) = Tmat(k + 1)
            else
               zTsn(k) = 0.0_wp
            end if
            if (l_brine .and. zTsn(k) > 0.0_wp) zTsn(k) = min(zTsn(k), 0.0_wp)
            zTsn(k) = zTsn(k) + avg_Tsi*P5*(Tsn_start(k) - zTsn(k))
            zqsn(k) = -c%rhos*(c%Lfresh - c%cp_ice*zTsn(k))
            enew = enew + hslyr*zqsn(k)
            Tsn_start(k) = zTsn(k)
         end do

         do k = 1, nilyr
            dqmat(k) = 0.0_wp
            reduce_kh(k) = .false.
         end do
         do k = 1, nilyr
            zTin(k) = Tmat(k + 1 + nslyr)
            if (l_brine .and. zTin(k) > Tmlts(k) - c%puny) then
               dTmat = zTin(k) - Tmlts(k)
               dqmat(k) = c%rhoi*dTmat*(c%cp_ice - c%Lfresh*Tmlts(k)/zTin(k)**2)
               zTin(k) = Tmlts(k)
               reduce_kh(k) = .true.
            end if
            zTin(k) = zTin(k) + avg_Tsi*P5*(Tin_start(k) - zTin(k))
            if (l_brine) then
               zqin(k) = -c%rhoi*(c%cp_ice*(Tmlts(k) - zTin(k)) &
                                  + c%Lfresh*(1.0_wp - Tmlts(k)/zTin(k)) - c%cp_ocn*Tmlts(k))
            else
               zqin(k) = -c%rhoi*(-c%cp_ice*zTin(k) + c%Lfresh)
            end if
            enew = enew + hilyr*zqin(k)
            einex_sfc_calc = einex_sfc_calc + hilyr*dqmat(k)
            Tin_start(k) = zTin(k)
         end do

         ! condition 3: large change in Tsf
         if (abs(dTsf) > ICE_TSF_ERRMAX) converged = .false.

         ! condition 4: Tsf >= 0 with fsurfn < fcondtopn
         fsurfn = fsurfn + dTsf*dfsurf_dT
         if (l_snow) then
            fcondtopn = kh(1)*(Tsf - zTsn(1))
         else
            fcondtopn = kh(1 + nslyr)*(Tsf - zTin(1))
         end if
         if (Tsf >= 0.0_wp .and. fsurfn < fcondtopn) converged = .false.
         dTsf_prev = dTsf

         ! condition 5: energy conservation
         fcondbot = kh(1 + nslyr + nilyr)*(zTin(nilyr) - Tbot)
         fcondbot = fcondbot + einex_sfc_calc/dt
         ferr = abs((enew - einit)/dt - (fcondtopn - fcondbot + fswint))

         ! a NaN/Inf anywhere here would compare false everywhere and
         ! masquerade as convergence
         if (.not. (ieee_is_finite(ferr) .and. ieee_is_finite(Tsf) .and. ieee_is_finite(dTsf) &
                    .and. ieee_is_finite(fsurfn) .and. ieee_is_finite(enew))) then
            ierr = ior(ierr, ICE_ERR_NONFINITE)
            return
         end if

         if (ferr > 0.9_wp*ICE_FERRMAX) then
            converged = .false.
            do k = 1, nilyr
               if (reduce_kh(k) .and. dqmat(k) > 0.0_wp) then
                  frac = max(0.5_wp*(1.0_wp - ferr/abs(fcondtopn - fcondbot)), 0.1_wp)
                  kh(k + nslyr + 1) = kh(k + nslyr + 1)*frac
                  kh(k + nslyr) = kh(k + nslyr + 1)
               end if
            end do
         end if

         if (converged) exit
      end do

      ! fluxes were evaluated before the last Tsf update
      flwoutn = flwoutn + dTsf_prev*dflwout_dT
      fsensn = fsensn + dTsf_prev*dfsens_dT
      flatn = flatn + dTsf_prev*dflat_dT

      if (.not. (ieee_is_finite(flwoutn) .and. ieee_is_finite(fsensn) .and. ieee_is_finite(flatn) &
                 .and. ieee_is_finite(fcondtopn) .and. ieee_is_finite(fcondbot))) then
         ierr = ior(ierr, ICE_ERR_NONFINITE)
         return
      end if

      einex_sfc_flux = 0.0_wp
      diag(1) = real(niter_done, wp)
      diag(2) = ferr
      diag(3) = dTsf_prev
      if (.not. converged) ierr = ior(ierr, ICE_ERR_TSFC_NOCONV)
   end subroutine ice_temperature_changes

   pure subroutine ice_temperature_changes_cell(i, j, n, nx, ny, ncat, nilyr, nslyr, &
                                                dt, rhoa, flw, potT, Qa, shcoef, lhcoef, &
                                                emissivity, sfc_kind, sf_0, dsf_dT, c, conduct, l_brine, &
                                                sw_redist, sw_frac, sw_dtemp, nitermax, Tbot, &
                                                hilyr, hslyr, zSin, einit, &
                                                fswsfc, fswint, Sswabs, Iswabs, zqin, zTin, zqsn, zTsn, Tsf, &
                                                fsensn, flatn, flwoutn, fsurfn, fcondtopn, fcondbot, &
                                                einex_sfc_flux, niter, ierr)
      !! The category `(i, j, n)` of grid-shaped state arrays through
      !! `ice_temperature_changes`: the driver shape
      !! `do concurrent (n = 1:ncat, j = js:je, i = is:ie)`.  The column is
      !! gathered ONCE into column-sized locals (the layer arrays are strided
      !! `(nx, ny, ncat, nlyr)`, the core wants contiguous layers), solved with
      !! the whole iteration inside, and scattered once.  State arrays are
      !! updated in place; the atmosphere and the controls are scalars (the
      !! caller reads its per-cell fields at `(i, j)`).  `niter` (integer) receives the
      !! pass count; `ierr(i, j, n)` is OR-ed into.  On an error the state
      !! arrays are untouched (non-finite input) or hold the last iterate
      !! (`ICE_ERR_TSFC_NOCONV`), as for the core routine.
      !$acc routine seq
      integer, value :: i, j, n
         !! Loop indices BY VALUE: by reference, nvfortran -stdpar=gpu reads the wrong cells.
      integer, intent(in) :: nx, ny, ncat, nilyr, nslyr
      real(wp), intent(in) :: dt, rhoa, flw, potT, Qa, shcoef, lhcoef, emissivity, sf_0, dsf_dT, Tbot
      integer, intent(in) :: sfc_kind, conduct, nitermax
      type(ice_constants_t), intent(in) :: c
      logical, intent(in) :: l_brine, sw_redist
      real(wp), intent(in) :: sw_frac, sw_dtemp
      real(wp), intent(in) :: hilyr(nx, ny, ncat), hslyr(nx, ny, ncat), einit(nx, ny, ncat)
      real(wp), intent(in) :: zSin(nx, ny, ncat, nilyr)
      real(wp), intent(inout) :: fswsfc(nx, ny, ncat), fswint(nx, ny, ncat), Tsf(nx, ny, ncat)
      real(wp), intent(inout) :: Sswabs(nx, ny, ncat, nslyr), zqsn(nx, ny, ncat, nslyr), zTsn(nx, ny, ncat, nslyr)
      real(wp), intent(inout) :: Iswabs(nx, ny, ncat, nilyr), zqin(nx, ny, ncat, nilyr), zTin(nx, ny, ncat, nilyr)
      real(wp), intent(inout) :: fsensn(nx, ny, ncat), flatn(nx, ny, ncat), flwoutn(nx, ny, ncat)
      real(wp), intent(inout) :: fsurfn(nx, ny, ncat), fcondtopn(nx, ny, ncat), fcondbot(nx, ny, ncat)
      real(wp), intent(inout) :: einex_sfc_flux(nx, ny, ncat)
      integer(int32), intent(inout) :: niter(nx, ny, ncat), ierr(nx, ny, ncat)
      integer :: k
      real(wp) :: S_c(ICE_NILYR_MAX), T_c(ICE_NILYR_MAX), q_c(ICE_NILYR_MAX), I_c(ICE_NILYR_MAX)
      real(wp) :: Ts_c(ICE_NSLYR_MAX), qs_c(ICE_NSLYR_MAX), Ss_c(ICE_NSLYR_MAX)
      real(wp) :: Tsf_c, fsw_c, fsi_c, fsens_c, flat_c, flwout_c, fsurf_c, fct_c, fcb_c, einex_c
      real(wp) :: diag_c(3)
      integer(int32) :: ierr_c

      if (nilyr < 1 .or. nilyr > ICE_NILYR_MAX .or. nslyr < 1 .or. nslyr > ICE_NSLYR_MAX) then
         ierr(i, j, n) = ior(ierr(i, j, n), ICE_ERR_NONFINITE)
         return
      end if
      ! Gather what the solver READS; the enthalpies and the flux outputs are only written.
      do k = 1, nilyr
         S_c(k) = zSin(i, j, n, k)
         T_c(k) = zTin(i, j, n, k)
         I_c(k) = Iswabs(i, j, n, k)
      end do
      do k = 1, nslyr
         Ts_c(k) = zTsn(i, j, n, k)
         Ss_c(k) = Sswabs(i, j, n, k)
      end do
      Tsf_c = Tsf(i, j, n)
      fsw_c = fswsfc(i, j, n)
      fsi_c = fswint(i, j, n)
      ierr_c = 0_int32

      call ice_temperature_changes(nilyr, nslyr, dt, rhoa, flw, potT, Qa, shcoef, lhcoef, &
                                   emissivity, sfc_kind, sf_0, dsf_dT, c, conduct, l_brine, &
                                   sw_redist, sw_frac, sw_dtemp, nitermax, &
                                   hilyr(i, j, n), hslyr(i, j, n), S_c, Tbot, einit(i, j, n), &
                                   fsw_c, fsi_c, Ss_c, I_c, q_c, T_c, qs_c, Ts_c, Tsf_c, &
                                   fsens_c, flat_c, flwout_c, fsurf_c, fct_c, fcb_c, einex_c, diag_c, ierr_c)

      ! Commit unless THIS call raised the non-finite bit (bits set before the call do not block it).
      if (iand(ierr_c, ICE_ERR_NONFINITE) == 0_int32) then
         do k = 1, nilyr
            zTin(i, j, n, k) = T_c(k)
            zqin(i, j, n, k) = q_c(k)
            Iswabs(i, j, n, k) = I_c(k)
         end do
         do k = 1, nslyr
            zTsn(i, j, n, k) = Ts_c(k)
            zqsn(i, j, n, k) = qs_c(k)
            Sswabs(i, j, n, k) = Ss_c(k)
         end do
         Tsf(i, j, n) = Tsf_c
         fswsfc(i, j, n) = fsw_c
         fswint(i, j, n) = fsi_c
         fsensn(i, j, n) = fsens_c
         flatn(i, j, n) = flat_c
         flwoutn(i, j, n) = flwout_c
         fsurfn(i, j, n) = fsurf_c
         fcondtopn(i, j, n) = fct_c
         fcondbot(i, j, n) = fcb_c
         einex_sfc_flux(i, j, n) = einex_c
         niter(i, j, n) = nint(diag_c(1), int32)
      end if
      ierr(i, j, n) = ior(ierr(i, j, n), ierr_c)
   end subroutine ice_temperature_changes_cell

end module rdb_ice_bl99
