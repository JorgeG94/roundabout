!! Thickness-category bounds of the linear ice-thickness distribution.
module rdb_ice_itd_params
   ! SPDX-License-Identifier: BSD-3-Clause
   !
   ! Ported from CICE/Icepack (CICE Consortium, https://github.com/CICE-Consortium).
   ! Upstream: Icepack 9a1f3e4 (Icepack 1.5.3+14), module icepack_itd
   ! (category-limit initialisation).
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
   ! Physics: Lipscomb (2001, JGR 106 13989); Bitz et al. (2001, JGR 106 2441);
   ! Thorndike et al. (1975, JGR 80 4501).
   !
   !! **Category limits.**  `itd_category_bounds` fills `hin_max(0:ncat)`, the
   !! UPPER edges of the thickness categories (`hin_max(0)` the lower edge of
   !! category 1), by one of five formulas selected by `kcatbound`:
   !!
   !! * `-1`  single category (`ncat = 1` only), edges `0` and `100` m (the top
   !!         edge is `999.9` m once `kitd = 1` applies the rule below);
   !! * ` 0`  the original Lipscomb (2001) tanh spacing,
   !!         `H(n) = H(n-1) + cc1 + cc2*(1 + tanh(cc3*((n-1)/ncat - 1)))`,
   !!         with `cc1 = 3/ncat, cc2 = 15*cc1, cc3 = 3, H(0) = 0` for the
   !!         linear remapping (`kitd = 1`) and `cc1 = max(1.1/ncat, hi_min),
   !!         cc2 = 25*cc1, cc3 = 2.25, H(0) = hi_min` for the delta-function
   !!         distribution (`kitd = 0`);
   !! * ` 1`  `H(n) = n*(d1 + (n-1)*d2)`, `d1 = 3/ncat`, `d2 = 0.5/ncat`
   !!         (round numbers; the upstream default);
   !! * ` 2`  World Meteorological Organization classes (`ncat` = 5, 6 or 7);
   !! * ` 3`  asymptotic widths `b1*(n + b3*n*(n+1)/(2*ncat*(ncat-1)))`,
   !!         `b1 = 0.1`, `b2 = 3`, `b3 = max(ncat*(ncat-1), 2*b2/b1)`
   !!         (`ncat >= 2`: the formula divides by `ncat*(ncat-1)`).
   !!
   !! With `kitd = 1` the top edge is then set to the arbitrary large `999.9` m.
   !! The routine is host/setup code (it calls `tanh`, through the safe-math
   !! facade).  Inconsistent selections return a status instead of aborting;
   !! the configuration validator turns a non-zero status into a refusal.
   use rdb_constants, only: wp
   use rdb_safe_math, only: safe_tanh
   use rdb_ice_dims, only: ICE_NCAT_MAX
   use rdb_ice_params, only: ice_params_t
   implicit none
   private

   public :: ITD_BOUNDS_OK, ITD_BOUNDS_BAD_KCATBOUND, ITD_BOUNDS_BAD_NCAT, ITD_BOUNDS_BAD_KITD
   public :: ITD_HIN_TOP
   public :: itd_category_bounds, ice_itd_bounds

   integer, parameter :: ITD_BOUNDS_OK = 0
      !! Bounds were built.
   integer, parameter :: ITD_BOUNDS_BAD_KCATBOUND = 1
      !! `kcatbound` is not one of -1..3.
   integer, parameter :: ITD_BOUNDS_BAD_NCAT = 2
      !! `ncat` is outside `1..ICE_NCAT_MAX`, or incompatible with `kcatbound`
      !! (`-1` needs `ncat = 1`, `2` needs `ncat` in 5, 6, 7).  `hin_max(0:1)` is
      !! still set for `kcatbound = -1` (the two upstream values), all else zero.
   integer, parameter :: ITD_BOUNDS_BAD_KITD = 3
      !! `kitd` is not 0 or 1.

   integer, parameter :: N_WMO5 = 5
   integer, parameter :: N_WMO6 = 6
   integer, parameter :: N_WMO7 = 7

   real(wp), parameter :: ITD_HIN_TOP = 999.9_wp
      !! Top edge of the last category under the linear remapping (`kitd = 1`).

contains

   subroutine itd_category_bounds(kcatbound, kitd, ncat, hi_min, hin_max, status)
      !! Category limits `hin_max(0:ncat)` [m] for `kcatbound` (-1..3) and `kitd`
      !! (0 delta function, 1 linear remapping); `hi_min` [m] enters only the
      !! `kitd = 0` branch of `kcatbound = 0`.  `status` is `ITD_BOUNDS_*`.
      integer, intent(in) :: kcatbound, kitd, ncat
      real(wp), intent(in) :: hi_min
      real(wp), intent(out) :: hin_max(0:ncat)
      integer, intent(out) :: status

      integer :: n
      real(wp) :: cc1, cc2, cc3, x1, rn, rncat, d1, d2, b1, b2, b3
      real(wp), parameter :: wmo5(N_WMO5) = [0.30_wp, 0.70_wp, 1.20_wp, 2.00_wp, 999.0_wp]
         !! thinnest 3 classes combined
      real(wp), parameter :: wmo6(N_WMO6) = [0.15_wp, 0.30_wp, 0.70_wp, 1.20_wp, 2.00_wp, 999.0_wp]
         !! thinnest 2 classes combined
      real(wp), parameter :: wmo7(N_WMO7) = [0.10_wp, 0.15_wp, 0.30_wp, 0.70_wp, 1.20_wp, 2.00_wp, 999.0_wp]
         !! all classes

      status = ITD_BOUNDS_OK
      hin_max = 0.0_wp
      if (kitd /= 0 .and. kitd /= 1) then
         status = ITD_BOUNDS_BAD_KITD
         return
      end if
      if (kcatbound < -1 .or. kcatbound > 3) then
         status = ITD_BOUNDS_BAD_KCATBOUND
         return
      end if
      if (ncat < 1 .or. ncat > ICE_NCAT_MAX) then
         status = ITD_BOUNDS_BAD_NCAT
         return
      end if

      rncat = real(ncat, wp)
      d1 = 3.0_wp/rncat
      d2 = 0.5_wp/rncat
      b1 = 0.1_wp          ! asymptotic category width [m]
      b2 = 3.0_wp          ! thickness at which the participation function is small [m]
      b3 = max(rncat*(rncat - 1.0_wp), 2.0_wp*b2/b1)

      select case (kcatbound)
      case (-1)
         hin_max(0) = 0.0_wp
         hin_max(1) = 100.0_wp
         if (ncat /= 1) status = ITD_BOUNDS_BAD_NCAT
      case (0)
         if (kitd == 1) then
            cc1 = 3.0_wp/rncat
            cc2 = 15.0_wp*cc1
            cc3 = 3.0_wp
            hin_max(0) = 0.0_wp
         else
            cc1 = max(1.1_wp/rncat, hi_min)
            cc2 = 25.0_wp*cc1
            cc3 = 2.25_wp
            hin_max(0) = hi_min
         end if
         do n = 1, ncat
            x1 = real(n - 1, wp)/rncat
            hin_max(n) = hin_max(n - 1) + cc1 + cc2*(1.0_wp + safe_tanh(cc3*(x1 - 1.0_wp)))
         end do
      case (1)
         hin_max(0) = 0.0_wp
         do n = 1, ncat
            rn = real(n, wp)
            hin_max(n) = rn*(d1 + (rn - 1.0_wp)*d2)
         end do
      case (2)
         hin_max(0) = 0.0_wp
         select case (ncat)
         case (5)
            hin_max(1:N_WMO5) = wmo5
         case (6)
            hin_max(1:N_WMO6) = wmo6
         case (7)
            hin_max(1:N_WMO7) = wmo7
         case default
            status = ITD_BOUNDS_BAD_NCAT
         end select
      case (3)
         hin_max(0) = 0.0_wp
         if (ncat == 1) status = ITD_BOUNDS_BAD_NCAT   ! the formula divides by ncat*(ncat-1)
         do n = 1, ncat
            rn = real(n, wp)
            if (ncat > 1) hin_max(n) = b1*(rn + b3*rn*(rn + 1.0_wp)/(2.0_wp*rncat*(rncat - 1.0_wp)))
         end do
      case default
         status = ITD_BOUNDS_BAD_KCATBOUND   ! unreachable: validated above
      end select

      if (kitd == 1 .and. status == ITD_BOUNDS_OK) hin_max(ncat) = ITD_HIN_TOP
   end subroutine itd_category_bounds

   subroutine ice_itd_bounds(p, hin_max, status)
      !! `itd_category_bounds` driven by the parameter bundle (`p%ncat`,
      !! `p%itd%kitd`, `p%itd%kcatbound`, `p%thermo%hi_min`).  `hin_max` must
      !! have room for `0:p%ncat`.
      type(ice_params_t), intent(in) :: p
      real(wp), intent(out) :: hin_max(0:)
      integer, intent(out) :: status

      if (size(hin_max) < p%ncat + 1) then
         hin_max = 0.0_wp
         status = ITD_BOUNDS_BAD_NCAT
         return
      end if
      call itd_category_bounds(p%itd%kcatbound, p%itd%kitd, p%ncat, p%thermo%hi_min, &
                               hin_max(0:p%ncat), status)
   end subroutine ice_itd_bounds

end module rdb_ice_itd_params
