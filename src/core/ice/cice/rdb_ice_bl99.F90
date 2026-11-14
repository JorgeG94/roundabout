!! Building blocks of the Bitz-Lipscomb (1999) vertical heat solve: the layer
!! conductivities, the tridiagonal system rows for the snow/ice column and the
!! Thomas solver.  The iteration that drives them lives in the surface-balance
!! PR that extends this module.
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
   use rdb_ice_errors, only: ICE_ERR_NONFINITE
   use rdb_ice_params, only: ice_constants_t, ICE_CONDUCT_MU71
   implicit none
   private

   public :: ICE_NMAT_MAX, ICE_KIMIN, ICE_BETAK
   public :: ice_layer_conductivity, ice_conductivity
   public :: ice_matrix_calc_Tsfc, ice_tridiag_solver

   integer, parameter :: ICE_NMAT_MAX = ICE_NILYR_MAX + ICE_NSLYR_MAX + 1
      !! Largest matrix dimension (`nilyr + nslyr + 1`).
   real(wp), parameter :: ICE_KIMIN = 0.10_wp
      !! Minimum conductivity of saline ice [W/m/K].
   real(wp), parameter :: ICE_BETAK = 0.13_wp
      !! Salinity coefficient of the MU71 conductivity [W/m/K/ppt].
   real(wp), parameter :: ICE_K_BUBBLY_RHOREF = 917.0_wp
      !! Reference density of the bubbly-brine fit [kg/m3].

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
      !! `hilyr > 0` and, with snow, `hslyr > 0`.
      !$acc routine seq
      logical, intent(in) :: l_snow
      integer, intent(in) :: nilyr, nslyr, conduct
      real(wp), intent(in) :: hilyr, hslyr
      real(wp), intent(in) :: zTin(nilyr), zSin(nilyr)
      type(ice_constants_t), intent(in) :: c
      real(wp), intent(inout) :: kh(nilyr + nslyr + 1)
      integer(int32), intent(inout) :: ierr
      integer :: k
      real(wp) :: kslyr
      real(wp) :: kilyr(ICE_NILYR_MAX), khl(ICE_NMAT_MAX)

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

      kslyr = c%ksno
      do k = 1, nilyr
         kilyr(k) = ice_layer_conductivity(zTin(k), zSin(k), c, conduct)
      end do

      if (l_snow) then
         khl(1) = 2.0_wp*kslyr/hslyr
         khl(1 + nslyr) = 2.0_wp*kslyr*kilyr(1)/(kslyr*hilyr + kilyr(1)*hslyr)
      else
         khl(1) = 0.0_wp
         khl(1 + nslyr) = 2.0_wp*kilyr(1)/hilyr
      end if
      khl(1 + nslyr + nilyr) = 2.0_wp*kilyr(nilyr)/hilyr

      if (nslyr > 1) then
         do k = 2, nslyr
            if (l_snow) then
               khl(k) = 2.0_wp*kslyr*kslyr/((kslyr + kslyr)*hslyr)
            else
               khl(k) = 0.0_wp
            end if
         end do
      end if

      do k = 2, nilyr
         khl(k + nslyr) = 2.0_wp*kilyr(k - 1)*kilyr(k)/((kilyr(k - 1) + kilyr(k))*hilyr)
      end do

      do k = 1, nilyr + nslyr + 1
         if (.not. ieee_is_finite(khl(k))) then
            ierr = ior(ierr, ICE_ERR_NONFINITE)
            return
         end if
      end do
      do k = 1, nilyr + nslyr + 1
         kh(k) = khl(k)
      end do
   end subroutine ice_conductivity

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
      integer :: k, ki, kr, nmat
      real(wp) :: sb(ICE_NMAT_MAX), dg(ICE_NMAT_MAX), sp(ICE_NMAT_MAX), rh(ICE_NMAT_MAX)

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

      ! Dummy rows (solution T = 0) for the surface and snow rows.
      do k = 1, nslyr + 1
         sb(k) = 0.0_wp
         dg(k) = 1.0_wp
         sp(k) = 0.0_wp
         rh(k) = 0.0_wp
      end do

      ! Surface equation, cold surface (with or without snow).
      if (l_cold) then
         if (l_snow) then
            k = 1
         else
            k = 1 + nslyr
         end if
         kr = k
         sb(kr) = 0.0_wp
         dg(kr) = dfsurf_dT - kh(k)
         sp(kr) = kh(k)
         rh(kr) = dfsurf_dT*Tsf - fsurfn
      end if

      ! Top snow layer.
      if (l_snow) then
         if (l_cold) then
            sb(2) = -etas(1)*kh(1)
            sp(2) = -etas(1)*kh(2)
            dg(2) = 1.0_wp + etas(1)*(kh(1) + kh(2))
            rh(2) = Tsn_init(1) + etas(1)*Sswabs(1)
         else
            sb(2) = 0.0_wp
            sp(2) = -etas(1)*kh(2)
            dg(2) = 1.0_wp + etas(1)*(kh(1) + kh(2))
            rh(2) = Tsn_init(1) + etas(1)*kh(1)*Tsf + etas(1)*Sswabs(1)
         end if
      end if

      ! Remaining snow layers.
      if (nslyr > 1) then
         do k = 2, nslyr
            kr = k + 1
            if (l_snow) then
               sb(kr) = -etas(k)*kh(k)
               sp(kr) = -etas(k)*kh(k + 1)
               dg(kr) = 1.0_wp + etas(k)*(kh(k) + kh(k + 1))
               rh(kr) = Tsn_init(k) + etas(k)*Sswabs(k)
            end if
         end do
      end if

      if (nilyr > 1) then
         ! Top ice layer.
         ki = 1
         k = ki + nslyr
         kr = k + 1
         if (l_snow .or. l_cold) then
            sb(kr) = -etai(ki)*kh(k)
            sp(kr) = -etai(ki)*kh(k + 1)
            dg(kr) = 1.0_wp + etai(ki)*(kh(k) + kh(k + 1))
            rh(kr) = Tin_init(ki) + etai(ki)*Iswabs(ki)
         else
            sb(kr) = 0.0_wp
            sp(kr) = -etai(ki)*kh(k + 1)
            dg(kr) = 1.0_wp + etai(ki)*(kh(k) + kh(k + 1))
            rh(kr) = Tin_init(ki) + etai(ki)*Iswabs(ki) + etai(ki)*kh(k)*Tsf
         end if

         ! Bottom ice layer.
         ki = nilyr
         k = ki + nslyr
         kr = k + 1
         sb(kr) = -etai(ki)*kh(k)
         sp(kr) = 0.0_wp
         dg(kr) = 1.0_wp + etai(ki)*(kh(k) + kh(k + 1))
         rh(kr) = Tin_init(ki) + etai(ki)*Iswabs(ki) + etai(ki)*kh(k + 1)*Tbot
      else
         ! Single ice layer.
         ki = 1
         k = ki + nslyr
         kr = k + 1
         if (l_snow .or. l_cold) then
            sb(kr) = -etai(ki)*kh(k)
            sp(kr) = 0.0_wp
            dg(kr) = 1.0_wp + etai(ki)*(kh(k) + kh(k + 1))
            rh(kr) = Tin_init(ki) + etai(ki)*Iswabs(ki) + etai(ki)*kh(k + 1)*Tbot
         else
            sb(kr) = 0.0_wp
            sp(kr) = 0.0_wp
            dg(kr) = 1.0_wp + etai(ki)*(kh(k) + kh(k + 1))
            rh(kr) = Tin_init(ki) + etai(ki)*Iswabs(ki) + etai(ki)*kh(k)*Tsf + etai(ki)*kh(k + 1)*Tbot
         end if
      end if

      ! Interior ice layers.
      do ki = 2, nilyr - 1
         k = ki + nslyr
         kr = k + 1
         sb(kr) = -etai(ki)*kh(k)
         sp(kr) = -etai(ki)*kh(k + 1)
         dg(kr) = 1.0_wp + etai(ki)*(kh(k) + kh(k + 1))
         rh(kr) = Tin_init(ki) + etai(ki)*Iswabs(ki)
      end do

      do k = 1, nmat
         sbdiag(k) = sb(k)
         diag(k) = dg(k)
         spdiag(k) = sp(k)
         rhs(k) = rh(k)
      end do
   end subroutine ice_matrix_calc_Tsfc

   pure subroutine ice_tridiag_solver(nmat, sbdiag, diag, spdiag, rhs, xout, ierr)
      !! Thomas algorithm for `sbdiag(k)*x(k-1) + diag(k)*x(k) +
      !! spdiag(k)*x(k+1) = rhs(k)`, `k = 1..nmat` (`sbdiag(1)` and
      !! `spdiag(nmat)` unused), no pivoting: the BL99 rows are diagonally
      !! dominant.  A non-finite coefficient or a zero/non-finite pivot sets
      !! `ICE_ERR_NONFINITE` and leaves `xout` untouched.
      !$acc routine seq
      integer, intent(in) :: nmat
      real(wp), intent(in) :: sbdiag(nmat), diag(nmat), spdiag(nmat), rhs(nmat)
      real(wp), intent(inout) :: xout(nmat)
      integer(int32), intent(inout) :: ierr
      integer :: k
      real(wp) :: wbeta
      real(wp) :: wgamma(ICE_NMAT_MAX), x(ICE_NMAT_MAX)

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

      wbeta = diag(1)
      if (.not. (abs(wbeta) > tiny(1.0_wp))) then
         ierr = ior(ierr, ICE_ERR_NONFINITE)
         return
      end if
      x(1) = rhs(1)/wbeta
      wgamma(1) = 0.0_wp
      do k = 2, nmat
         wgamma(k) = spdiag(k - 1)/wbeta
         wbeta = diag(k) - sbdiag(k)*wgamma(k)
         if (.not. (abs(wbeta) > tiny(1.0_wp))) then
            ierr = ior(ierr, ICE_ERR_NONFINITE)
            return
         end if
         x(k) = (rhs(k) - sbdiag(k)*x(k - 1))/wbeta
      end do
      do k = nmat - 1, 1, -1
         x(k) = x(k) - wgamma(k + 1)*x(k + 1)
      end do

      do k = 1, nmat
         if (.not. ieee_is_finite(x(k))) then
            ierr = ior(ierr, ICE_ERR_NONFINITE)
            return
         end if
      end do
      do k = 1, nmat
         xout(k) = x(k)
      end do
   end subroutine ice_tridiag_solver

end module rdb_ice_bl99
