!! Ridging drivers: the ice-strength field.
module rdb_ice_ridge
   ! SPDX-License-Identifier: BSD-3-Clause
   !
   ! Ported from CICE/Icepack (CICE Consortium, https://github.com/CICE-Consortium).
   ! Upstream: Icepack 9a1f3e4 (Icepack 1.5.3+14), module icepack_mechred
   ! (icepack_ice_strength).
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
   ! Physics: Rothrock (1975, JGR 80 4514); Hibler (1979, JPO 9 817).
   !
   !! **`ice_compute_strength`.**  Fills `strength(i,j)` [N/m] over the index
   !! window `(is:ie, js:je)` (halo cells included when the caller widens the
   !! window: each cell is computed locally from the state it is given, no
   !! exchange) from the category state, one `do concurrent` work item per cell
   !! calling `ice_strength_cell`.  A cell outside `wet_mask`, or empty
   !! (`aice <= puny`: open ocean, skipped), gets exactly zero.  The test is
   !! explicit (`ieee_is_nan(aice) .or. aice > puny`) so a NaN `aice` is NOT skipped and reaches the cell kernel's
   !! non-finite guard.  Per-cell error bits are ORed into `ice_err`.  The
   !! dynamics' own ice mask (`dyn_area_min`, mass) is applied by the caller.
   !! All arrays must be device-resident on the offload build.
   use, intrinsic :: iso_fortran_env, only: int32
   use rdb_constants, only: wp
   use rdb_ice_dims, only: ICE_NCAT_MAX
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
   use rdb_ice_params, only: ice_ridge_params_t, ice_constants_t, ice_ridge_params_valid
   use rdb_ice_ridge_itd, only: ice_strength_cell
   implicit none
   private

   public :: ice_compute_strength

contains

   pure subroutine strength_at(nx, ny, ncat, i, j, rp, c, aice, vice, aice0, aicen, vicen, strength, ierr)
      !! One cell: gather the category column from the `(nx,ny,ncat)` arrays
      !! into fixed-size locals and call the cell kernel.
      !$acc routine seq
      integer, intent(in) :: nx, ny, ncat, i, j
      type(ice_ridge_params_t), intent(in) :: rp
      type(ice_constants_t), intent(in) :: c
      real(wp), intent(in) :: aice(nx, ny), vice(nx, ny), aice0(nx, ny)
      real(wp), intent(in) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat)
      real(wp), intent(out) :: strength
      integer(int32), intent(inout) :: ierr
      real(wp) :: an(ICE_NCAT_MAX), vn(ICE_NCAT_MAX)
      integer :: n

      do n = 1, ncat
         an(n) = aicen(i, j, n)
         vn(n) = vicen(i, j, n)
      end do
      call ice_strength_cell(ncat, rp, c, aice(i, j), vice(i, j), aice0(i, j), an, vn, strength, ierr)
   end subroutine strength_at

   subroutine ice_compute_strength(nx, ny, ncat, is, ie, js, je, rp, c, wet_mask, &
                                   aice, vice, aice0, aicen, vicen, strength, ice_err)
      !! Ice strength field over `(is:ie, js:je)` (see the module note).
      !! `ncat` must satisfy `1 <= ncat <= ICE_NCAT_MAX` and the options of `rp` be 0 or 1 (validated at
      !! configure); a violation stops the run, it is never reached by a
      !! validated configuration.
      integer, intent(in) :: nx, ny, ncat, is, ie, js, je
      type(ice_ridge_params_t), intent(in) :: rp
      type(ice_constants_t), intent(in) :: c
      logical, intent(in) :: wet_mask(nx, ny)
      real(wp), intent(in) :: aice(nx, ny), vice(nx, ny), aice0(nx, ny)
      real(wp), intent(in) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat)
      real(wp), intent(inout) :: strength(nx, ny)
      integer(int32), intent(inout) :: ice_err(nx, ny)
      integer :: i, j
      integer(int32) :: ierr_c
      real(wp) :: s_c, puny

      puny = c%puny
      if (.not. ice_ridge_params_valid(rp)) error stop "ice_compute_strength: kstrength/krdg_partic/krdg_redist must be 0 or 1"
      if (ncat < 1 .or. ncat > ICE_NCAT_MAX) error stop "ice_compute_strength: ncat outside 1..ICE_NCAT_MAX"

      do concurrent(j=js:je, i=is:ie) local(ierr_c, s_c)
         if (wet_mask(i, j) .and. (ieee_is_nan(aice(i, j)) .or. aice(i, j) > puny)) then
            ierr_c = ice_err(i, j)
            call strength_at(nx, ny, ncat, i, j, rp, c, aice, vice, aice0, aicen, vicen, s_c, ierr_c)
            strength(i, j) = s_c
            ice_err(i, j) = ierr_c
         else
            strength(i, j) = 0.0_wp
         end if
      end do
   end subroutine ice_compute_strength

end module rdb_ice_ridge
