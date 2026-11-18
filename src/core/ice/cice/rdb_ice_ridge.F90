!! Ridging drivers: the ice-strength field and the ridging step (ridge, then clean-up).
module rdb_ice_ridge
   ! SPDX-License-Identifier: BSD-3-Clause
   !
   ! Ported from CICE/Icepack (CICE Consortium, https://github.com/CICE-Consortium).
   ! Upstream: Icepack 9a1f3e4 (Icepack 1.5.3+14), module icepack_mechred
   ! (icepack_ice_strength, icepack_step_ridge).
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
   !!
   !! **`ice_step_ridge`.**  The ridging step of one dynamics step, run right after the transport
   !! and before the next aggregation or thermodynamics: ONE `do concurrent` over the owned WET
   !! cells, whose body is the ridging cell kernel (`ice_ridge_ice_cell`, `rdb_ice_ridge_shift`)
   !! followed, in the same thread, by the clean-up cell kernel (`ice_cleanup_itd_cell`; thermodynamic
   !! step `ndtd dt`, rebin on, area limit on, the cell aggregates refreshed once).  Cells are
   !! independent, so this is bitwise the two-launch sequence ridge sweep then clean-up sweep, with one
   !! launch and one read of the state.  A land cell is not touched at all (its `Tf` may be anything,
   !! its aggregates are not refreshed, `niter = 0`).  A cell that ALREADY carries a fatal error bit on entry
   !! is skipped (ridging and clean-up), and so is the clean-up of a cell whose ridging raised a fatal bit:
   !! the caller must stop the run on any fatal bit instead of using the cell.
   !! Inputs: the deformation rates `rdg_conv`, `rdg_shear` [1/s] of the same step (two `(nx,ny)` fields the
   !! dynamics must keep alive across the transport), `dt` the dynamics step, `ndtd` the subcycle count:
   !! it MUST be 1 (`ok = .false.` and nothing touched otherwise, as for an invalid ridging parameter set,
   !! `ncat` outside `1..ICE_NCAT_MAX`, or `dt <= 0`; checked once per call at entry, never in the launch).
   !! It updates the category state in place and ADDS to `fresh` (snow thrown into the ocean, plus the
   !! clean-up zaps), `fsalt` (the zaps only), `fhocn`, and to the diagnostics `dardg1dt`, `dardg2dt`,
   !! `dvirdgdt`, `opening`; `niter(i,j)` is the number of ridging passes (0: early-out, land or skipped).
   !! The pass-count histogram bit (`ICE_ERR_RIDGE_PASSES`) is set by the wiring PR from `niter`.
   use, intrinsic :: iso_fortran_env, only: int32
   use rdb_constants, only: wp
   use rdb_ice_dims, only: ICE_NCAT_MAX
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
   use rdb_ice_errors, only: ICE_ERR_FATAL_MASK
   use rdb_ice_params, only: ice_params_t, ice_ridge_params_t, ice_constants_t, ice_ridge_params_valid
   use rdb_ice_ridge_itd, only: ice_strength_cell, RIDGE_NITERMAX
   use rdb_ice_ridge_shift, only: ice_ridge_ice_cell
   use rdb_ice_itd_cleanup, only: ice_cleanup_itd_cell
   implicit none
   private

   public :: ice_compute_strength, ice_step_ridge

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

   subroutine ice_step_ridge(nx, ny, ncat, nilyr, nslyr, i_lo, i_hi, j_lo, j_hi, p, hin_max, dt, ndtd, wet_mask, Tf, &
                             rdg_conv, rdg_shear, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, &
                             aice, vice, vsno, aice0, Tsfc_c, qice_c, sice_c, qsno_c, &
                             fresh, fsalt, fhocn, dardg1dt, dardg2dt, dvirdgdt, opening, niter, ice_err, ok)
      !! Ridging then clean-up of every wet cell of `i_lo:i_hi, j_lo:j_hi` in one launch (see the module note).
      !! `hin_max(0:ncat)` is the category bound table, `Tf(nx,ny)` the freezing point.  Every array must be
      !! device-present on a GPU build.  `ok = .false.` (nothing touched) when `p%ridge` is invalid, `ncat` is
      !! outside `1..ICE_NCAT_MAX`, `ndtd /= 1` or `dt` is not positive and finite.
      integer, intent(in) :: nx, ny, ncat, nilyr, nslyr, i_lo, i_hi, j_lo, j_hi, ndtd
      type(ice_params_t), intent(in) :: p
      real(wp), intent(in) :: dt
      real(wp), intent(in) :: hin_max(0:ncat)
      logical, intent(in) :: wet_mask(nx, ny)
      real(wp), intent(in) :: Tf(nx, ny), rdg_conv(nx, ny), rdg_shear(nx, ny)
      real(wp), intent(inout) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat)
      real(wp), intent(inout) :: Tsfcn(nx, ny, ncat), qicen(nx, ny, ncat, nilyr)
      real(wp), intent(inout) :: sicen(nx, ny, ncat, nilyr), qsnon(nx, ny, ncat, nslyr)
      real(wp), intent(inout) :: aice(nx, ny), vice(nx, ny), vsno(nx, ny), aice0(nx, ny)
      real(wp), intent(inout) :: Tsfc_c(nx, ny), qice_c(nx, ny, nilyr), sice_c(nx, ny, nilyr), qsno_c(nx, ny, nslyr)
      real(wp), intent(inout) :: fresh(nx, ny), fsalt(nx, ny), fhocn(nx, ny)
      real(wp), intent(inout) :: dardg1dt(nx, ny), dardg2dt(nx, ny), dvirdgdt(nx, ny), opening(nx, ny)
      integer(int32), intent(inout) :: niter(nx, ny)
      integer(int32), intent(inout) :: ice_err(nx, ny)
      logical, intent(out) :: ok

      ok = ice_ridge_params_valid(p%ridge) .and. ncat >= 1 .and. ncat <= ICE_NCAT_MAX .and. ndtd == 1 .and. dt > 0.0_wp &
           .and. dt <= huge(dt)
      if (.not. ok) return

      call ridge_step_impl(nx, ny, ncat, nilyr, nslyr, i_lo, i_hi, j_lo, j_hi, p%ridge, p%c%puny, p%c%rhoi, p%c%rhos, &
                           p%c%Lfresh, p%c%cp_ice, p%c%hs_min, p%c%itd_area_min, p%c%itd_mass_min, &
                           p%thermo%ice_ref_salinity, dt, ndtd, RIDGE_NITERMAX, hin_max, wet_mask, Tf, rdg_conv, rdg_shear, &
                           aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, aice, vice, vsno, aice0, Tsfc_c, qice_c, sice_c, &
                           qsno_c, fresh, fsalt, fhocn, dardg1dt, dardg2dt, dvirdgdt, opening, niter, ice_err)
   end subroutine ice_step_ridge

   subroutine ridge_step_impl(nx, ny, ncat, nilyr, nslyr, i_lo, i_hi, j_lo, j_hi, rp, puny, rhoi, rhos, Lfresh, cp_ice, &
                              hs_min, area_min, mass_min, sref, dt, ndtd, itermax, hin_max, wet_mask, Tf, rdg_conv, &
                              rdg_shear, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, aice, vice, vsno, aice0, &
                              Tsfc_c, qice_c, sice_c, qsno_c, fresh, fsalt, fhocn, dardg1dt, dardg2dt, dvirdgdt, &
                              opening, niter, ierr)
      !! The one `do concurrent` launch of the ridging step: ridge cell, then clean-up cell, per wet cell.
      integer, intent(in) :: nx, ny, ncat, nilyr, nslyr, i_lo, i_hi, j_lo, j_hi, ndtd, itermax
      type(ice_ridge_params_t), intent(in) :: rp
      real(wp), intent(in) :: puny, rhoi, rhos, Lfresh, cp_ice, hs_min, area_min, mass_min, sref, dt
      real(wp), intent(in) :: hin_max(0:ncat)
      logical, intent(in) :: wet_mask(nx, ny)
      real(wp), intent(in) :: Tf(nx, ny), rdg_conv(nx, ny), rdg_shear(nx, ny)
      real(wp), intent(inout) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat)
      real(wp), intent(inout) :: Tsfcn(nx, ny, ncat), qicen(nx, ny, ncat, nilyr)
      real(wp), intent(inout) :: sicen(nx, ny, ncat, nilyr), qsnon(nx, ny, ncat, nslyr)
      real(wp), intent(inout) :: aice(nx, ny), vice(nx, ny), vsno(nx, ny), aice0(nx, ny)
      real(wp), intent(inout) :: Tsfc_c(nx, ny), qice_c(nx, ny, nilyr), sice_c(nx, ny, nilyr), qsno_c(nx, ny, nslyr)
      real(wp), intent(inout) :: fresh(nx, ny), fsalt(nx, ny), fhocn(nx, ny)
      real(wp), intent(inout) :: dardg1dt(nx, ny), dardg2dt(nx, ny), dvirdgdt(nx, ny), opening(nx, ny)
      integer(int32), intent(inout) :: niter(nx, ny)
      integer(int32), intent(inout) :: ierr(nx, ny)
      integer :: i, j, nit
      real(wp) :: dtt

      dtt = real(ndtd, wp)*dt
      !$acc data present(hin_max, wet_mask, Tf, rdg_conv, rdg_shear, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon) &
      !$acc& present(aice, vice, vsno, aice0, Tsfc_c, qice_c, sice_c, qsno_c, fresh, fsalt, fhocn) &
      !$acc& present(dardg1dt, dardg2dt, dvirdgdt, opening, niter, ierr)
      do concurrent(j=j_lo:j_hi, i=i_lo:i_hi) local(nit)
         niter(i, j) = 0
         if (wet_mask(i, j)) then
            if (iand(ierr(i, j), ICE_ERR_FATAL_MASK) == 0_int32) then
               call ice_ridge_ice_cell(i, j, nx, ny, ncat, nilyr, nslyr, rp, puny, rhos, dt, ndtd, itermax, hin_max, &
                                       Tf(i, j), rdg_conv(i, j), rdg_shear(i, j), aicen, vicen, vsnon, Tsfcn, qicen, &
                                       sicen, qsnon, aice0, fresh, fhocn, dardg1dt, dardg2dt, dvirdgdt, opening, ierr, nit)
               niter(i, j) = nit
               if (iand(ierr(i, j), ICE_ERR_FATAL_MASK) == 0_int32) then
                  call ice_cleanup_itd_cell(i, j, nx, ny, ncat, nilyr, nslyr, hin_max, dtt, Tf(i, j), rhoi, rhos, Lfresh, &
                                            cp_ice, hs_min, puny, area_min, mass_min, sref, .true., .true., .true., &
                                            aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, aice, vice, vsno, aice0, &
                                            Tsfc_c, qice_c, sice_c, qsno_c, fresh, fsalt, fhocn, ierr)
               end if
            end if
         end if
      end do
      !$acc end data
   end subroutine ridge_step_impl

end module rdb_ice_ridge
