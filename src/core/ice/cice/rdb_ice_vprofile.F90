!! Vertical profile of a snow/ice category at the start of the thermodynamic
!! step: layer thicknesses, layer temperatures from the stored enthalpies, the
!! brine-pocket melting temperatures, the range checks and the initial energy.
module rdb_ice_vprofile
   ! SPDX-License-Identifier: BSD-3-Clause
   !
   ! Ported from CICE/Icepack (CICE Consortium, https://github.com/CICE-Consortium).
   ! Upstream: Icepack 9a1f3e4 (Icepack 1.5.3+14), module icepack_therm_vertical
   ! (init_vertical_profile).
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
   ! C. M. Bitz (UW).
   ! Physics: Bitz and Lipscomb (1999, JGR 104 15669).
   !
   !! **What it does.**  From the category state (`aicen`, `vicen`, `vsnon`) and
   !! the layer enthalpies it returns the ice and snow thicknesses and layer
   !! thicknesses (`hin = vicen/aicen`, `hilyr = hin/nilyr`, `hslyr =
   !! hsn/nslyr`), the layer temperatures from the BL99 quadratic (ice) and
   !! `(Lfresh + q/rhos)/cp_ice` (snow), the melting temperatures `Tmlts =
   !! -depressT*S`, and `einit = sum(h*q)` relative to 0 degC.  Snow too thin to
   !! carry a temperature (`hslyr <= hs_min/nslyr`) is reset to the melting
   !! reference enthalpy `-rhos*Lfresh`.  Round-off excursions are repaired as
   !! upstream: a snow or fresh-ice temperature above 0 degC becomes 0 with the
   !! reference enthalpy; with brine, a layer at or above its melting point is
   !! put `puny` below it with `q = -S*depressT*rhoi*cp_ocn`.
   !!
   !! **Range checks (fatal bits).**  Snow `T > Tmax` (the round-off allowance
   !! `-q*puny*nslyr/(rhos*cp_ice*vsnon)`, `puny` for thin snow) or `T <
   !! ICE_TMIN` sets `ICE_ERR_TSN_RANGE`; `S < min_salin - puny` sets
   !! `ICE_ERR_SIN_LOW`; ice `T > Tmax` (`Tmlts` with brine, else the upstream
   !! allowance `-q*puny/(rhos*cp_ice*vicen)`, kept as is, snow density and
   !! all) or `T < ICE_TMIN` sets `ICE_ERR_TIN_RANGE`; a negative discriminant
   !! of the enthalpy inversion sets `ICE_ERR_QIN_SQRT`.  The first failing
   !! stage wins and returns (snow, then ice layers top-down).  The mushy
   !! thermodynamics option is not ported.
   !!
   !! **Kernel contract.**  `pure`, `!$acc routine seq`, explicit shapes,
   !! `ICE_NILYR_MAX`/`ICE_NSLYR_MAX` locals, constants in the `intent(in)`
   !! bundle, `ierr` only OR-ed into.  ERROR-RETURN CONTRACT: on ANY error
   !! (entry check or range check) the layer arrays `zqin, zqsn, zTin, zTsn,
   !! Tmlts` are left untouched and ALL scalar outputs (`hin, hilyr, hsn,
   !! hslyr, einit`) are zero.  Non-finite input, or a layer count out of
   !! range, sets `ICE_ERR_NONFINITE`; `aicen <= 0` or `vicen <= 0` (the caller
   !! gates on ice-covered categories; both divide) sets `ICE_ERR_NEG_THICK`.
   !! A finite `vsnon <= 0` (round-off negative snow volume) is "no snow", as
   !! upstream: the thin-snow reset applies.  Every comparison is guarded by an
   !! entry finiteness check because comparisons with NaN are false.
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use, intrinsic :: iso_fortran_env, only: int32
   use rdb_constants, only: wp
   use rdb_ice_dims, only: ICE_NILYR_MAX, ICE_NSLYR_MAX
   use rdb_ice_errors, only: ICE_ERR_NONFINITE, ICE_ERR_NEG_THICK, ICE_ERR_TSN_RANGE, ICE_ERR_SIN_LOW, &
                             ICE_ERR_TIN_RANGE, ICE_ERR_QIN_SQRT
   use rdb_ice_params, only: ice_constants_t, ice_T_from_q, ice_Tsno_from_q
   implicit none
   private

   public :: ICE_TMIN, ice_init_vertical_profile

   real(wp), parameter :: ICE_TMIN = -100.0_wp
      !! Lowest admissible internal temperature [degC].

contains

   pure subroutine ice_init_vertical_profile(aicen, vicen, vsnon, nilyr, nslyr, c, l_brine, zSin, &
                                             zqin, zqsn, zTin, zTsn, Tmlts, &
                                             hin, hilyr, hsn, hslyr, einit, ierr)
      !! See the module header.  `zqin`/`zqsn` [J/m3] are updated in place only
      !! where the profile repairs them; `zTin`, `zTsn`, `Tmlts` are outputs;
      !! `zSin(1:nilyr)` the prescribed layer salinities [ppt].
      !$acc routine seq
      real(wp), intent(in) :: aicen, vicen, vsnon
      integer, intent(in) :: nilyr, nslyr
      type(ice_constants_t), intent(in) :: c
      logical, intent(in) :: l_brine
      real(wp), intent(in) :: zSin(nilyr)
      real(wp), intent(inout) :: zqin(nilyr), zqsn(nslyr)
      real(wp), intent(inout) :: zTin(nilyr), zTsn(nslyr), Tmlts(nilyr)
      real(wp), intent(out) :: hin, hilyr, hsn, hslyr, einit
      integer(int32), intent(inout) :: ierr
      integer :: k
      real(wp) :: rnslyr, Tmax, ein
      real(wp) :: qn(ICE_NILYR_MAX), Tn(ICE_NILYR_MAX), Tm(ICE_NILYR_MAX)
      real(wp) :: qs(ICE_NSLYR_MAX), Ts(ICE_NSLYR_MAX)
      logical :: bad

      hin = 0.0_wp
      hilyr = 0.0_wp
      hsn = 0.0_wp
      hslyr = 0.0_wp
      einit = 0.0_wp

      if (nilyr < 1 .or. nilyr > ICE_NILYR_MAX .or. nslyr < 1 .or. nslyr > ICE_NSLYR_MAX) then
         ierr = ior(ierr, ICE_ERR_NONFINITE)
         return
      end if
      if (.not. (ieee_is_finite(aicen) .and. ieee_is_finite(vicen) .and. ieee_is_finite(vsnon))) then
         ierr = ior(ierr, ICE_ERR_NONFINITE)
         return
      end if
      if (.not. (aicen > 0.0_wp .and. vicen > 0.0_wp)) then
         ierr = ior(ierr, ICE_ERR_NEG_THICK)
         return
      end if
      do k = 1, nilyr
         if (.not. (ieee_is_finite(zqin(k)) .and. ieee_is_finite(zSin(k)))) then
            ierr = ior(ierr, ICE_ERR_NONFINITE)
            return
         end if
      end do
      do k = 1, nslyr
         if (.not. ieee_is_finite(zqsn(k))) then
            ierr = ior(ierr, ICE_ERR_NONFINITE)
            return
         end if
      end do

      rnslyr = real(nslyr, wp)
      ein = 0.0_wp

      ! Thicknesses.
      hin = vicen/aicen
      hsn = vsnon/aicen
      hilyr = hin/real(nilyr, wp)
      hslyr = hsn/rnslyr
      if (.not. (ieee_is_finite(hin) .and. ieee_is_finite(hsn))) then
         hin = 0.0_wp
         hilyr = 0.0_wp
         hsn = 0.0_wp
         hslyr = 0.0_wp
         einit = 0.0_wp
         ierr = ior(ierr, ICE_ERR_NONFINITE)
         return
      end if

      ! Snow: enthalpy of thin snow reset, temperature, bounds, round-off repair.
      do k = 1, nslyr
         qs(k) = zqsn(k)
         if (hslyr > c%hs_min/rnslyr) then
            Tmax = -qs(k)*c%puny*rnslyr/(c%rhos*c%cp_ice*vsnon)
         else
            qs(k) = -c%rhos*c%Lfresh
            Tmax = c%puny
         end if
         Ts(k) = ice_Tsno_from_q(qs(k), c)
         if (Ts(k) > Tmax .or. Ts(k) < ICE_TMIN .or. .not. ieee_is_finite(Ts(k))) then
            hin = 0.0_wp
            hilyr = 0.0_wp
            hsn = 0.0_wp
            hslyr = 0.0_wp
            einit = 0.0_wp
            ierr = ior(ierr, ICE_ERR_TSN_RANGE)
            return
         end if
         if (Ts(k) > 0.0_wp) then
            Ts(k) = 0.0_wp
            qs(k) = -c%rhos*c%Lfresh
         end if
         ein = ein + hslyr*qs(k)
      end do

      ! Ice layers, top-down.
      do k = 1, nilyr
         if (zSin(k) < c%min_salin - c%puny) then
            hin = 0.0_wp
            hilyr = 0.0_wp
            hsn = 0.0_wp
            hslyr = 0.0_wp
            einit = 0.0_wp
            ierr = ior(ierr, ICE_ERR_SIN_LOW)
            return
         end if
         Tm(k) = -zSin(k)*c%depressT
         qn(k) = zqin(k)
         call ice_T_from_q(qn(k), zSin(k), c, l_brine, Tn(k), bad)
         if (bad) then
            hin = 0.0_wp
            hilyr = 0.0_wp
            hsn = 0.0_wp
            hslyr = 0.0_wp
            einit = 0.0_wp
            ierr = ior(ierr, ICE_ERR_QIN_SQRT)
            return
         end if
         if (l_brine) then
            Tmax = Tm(k)
         else
            Tmax = -qn(k)*c%puny/(c%rhos*c%cp_ice*vicen)
         end if
         if (Tn(k) > Tmax .or. Tn(k) < ICE_TMIN .or. .not. ieee_is_finite(Tn(k))) then
            hin = 0.0_wp
            hilyr = 0.0_wp
            hsn = 0.0_wp
            hslyr = 0.0_wp
            einit = 0.0_wp
            ierr = ior(ierr, ICE_ERR_TIN_RANGE)
            return
         end if
         if (Tn(k) > 0.0_wp) then
            Tn(k) = 0.0_wp
            qn(k) = -c%rhoi*c%Lfresh
         end if
         if (Tn(k) >= -zSin(k)*c%depressT) then
            Tn(k) = -zSin(k)*c%depressT - c%puny
            qn(k) = -zSin(k)*c%depressT*c%rhoi*c%cp_ocn
         end if
         ein = ein + hilyr*qn(k)
      end do

      if (.not. ieee_is_finite(ein)) then
         hin = 0.0_wp
         hilyr = 0.0_wp
         hsn = 0.0_wp
         hslyr = 0.0_wp
         einit = 0.0_wp
         ierr = ior(ierr, ICE_ERR_NONFINITE)
         return
      end if

      do k = 1, nslyr
         zqsn(k) = qs(k)
         zTsn(k) = Ts(k)
      end do
      do k = 1, nilyr
         zqin(k) = qn(k)
         zTin(k) = Tn(k)
         Tmlts(k) = Tm(k)
      end do
      einit = ein

   end subroutine ice_init_vertical_profile

end module rdb_ice_vprofile
