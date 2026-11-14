!! Thermodynamic tail of an ice step (`therm2`): rain, thickness-distribution
!! remap, new ice, lateral melt and clean-up, in staged launches over the SoA state.
module rdb_ice_therm2
   ! SPDX-License-Identifier: BSD-3-Clause
   !
   ! Ported from CICE/Icepack (CICE Consortium, https://github.com/CICE-Consortium).
   ! Upstream: Icepack 9a1f3e4 (Icepack 1.5.3+14), module icepack_therm_itd
   ! (icepack_step_therm2: the ORDER of the stages and their conditions; the stages
   ! themselves live in the modules that own them).
   ! Translated to roundabout's do-concurrent/GPU form (one pure per-cell routine
   ! over the SoA state, integer error flags, no optional arguments or
   ! allocatables) and modified by the roundabout developers; see NOTICE for the
   ! full notice.
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
   ! Physics: Lipscomb (2001, JGR 106 13989); Bitz and Lipscomb (1999, JGR 104
   ! 15669); Maykut and Perovich (1987, JGR 92 7032); Steele (1992, JGR 97 17729).
   !
   !! **What it does, per owned wet cell, in upstream's order.**
   !! 1. *Rain.*  `fresh += frain * aice`, where `aice` is the cell aggregate the
   !!    state carried INTO the step (the start-of-step ice concentration, read
   !!    before anything below changes it): rain falls through the ice cover and
   !!    reaches the ocean at once, it is never stored in the ice.  `frain` is per
   !!    unit cell area of the atmosphere flux, so the product is per cell area.
   !! 2. *Aggregate and remap* (`ice_linear_itd_cell`, Lipscomb 2001).  The cell
   !!    area `sum(aicen)` is formed; where it exceeds `puny` (and the
   !!    thickness-distribution scheme is the linear remap, `kitd = 1`) the
   !!    categories are remapped from the start-of-step snapshots
   !!    `aicen_init`/`vicen_init` (which `ice_save_init` took BEFORE therm1)
   !!    to the post-thermodynamics thicknesses.  The in-remap rebin is OFF
   !!    (`rebin_fb = .false.`): the clean-up (6) rebins every cell, including the
   !!    ones whose remap fell back.
   !! 3. *New ice* (`ice_add_new_ice_cell`, Bitz and Lipscomb 1999) where the
   !!    freezing potential `frzmlt > 0`.
   !! 4. *Lateral melt* (`ice_lateral_melt_cell`, Steele 1992) with the per-cell
   !!    melted fraction `rside` that therm1's basal-heat routine produced.  It
   !!    runs AFTER the new ice, so it melts a share of the new ice too.
   !! 5. *Clean-up* (`ice_cleanup_itd_cell`): rebin of every category that left
   !!    its bin, zap of vanishing categories and of snow with a nonsensical
   !!    temperature (fluxes ADDED to `fresh`, `fsalt`, `fhocn`), and, as the LAST
   !!    action, the cell aggregates `aice, vice, vsno, aice0` and the cell-mean
   !!    tracers (the aggregate upstream's step driver runs once after therm2).
   !!    `dorebin = .true.` and `limit_aice = .true.`, the values upstream's
   !!    therm2 gets from the optional-argument defaults.
   !!
   !! **The `frzmlt` sign split.**  `frzmlt(nx, ny)` [W m-2] is ONE signed field
   !! (MASTER_PLAN C.4/8): positive where the ocean-side frazil clamp banked heat
   !! (`bank/dt`: the freezing potential), negative where the ocean can melt ice
   !! (the melt potential), never both.  This driver reads ONLY its positive part.
   !! The negative part was consumed upstream of it, in therm1, where the
   !! basal-heat routine caps `fbot + fside` at `-frzmlt` (the limiter `xtmp`) and
   !! produced `rside`; this driver gets the result of that as `rside` and never
   !! looks at a negative `frzmlt` (a negative or zero value is an exact no-op of
   !! the new-ice stage).  The bank is spent where `frzmlt > 0` (reset to 0: all
   !! of it became ice) and kept elsewhere (a cell with `frzmlt <= 0` but a
   !! positive bank is a caller inconsistency; the bank is conserved, not lost).
   !! New ice never writes `fhocn`: the ocean already paid the latent heat when
   !! it warmed the supercooled water back to freezing (no double count; the
   !! test pins it).
   !!
   !! **Launches.**  Four `do concurrent (j, i)` sweeps over the owned wet cells, each calling an
   !! existing `pure`, `!$acc routine seq` cell kernel in place on the SoA state (A: rain, finite check
   !! and remap; B: new ice and the bank; C: lateral melt; D: clean-up), or three when B and C share one
   !! sweep (`fuse_bc`, the default, see the shim).  Every stage cell starts with the same READ-ONLY gate: a cell
   !! that carries a FATAL error bit is skipped, so a cell that turns fatal in one launch is skipped by
   !! the following ones, and a cell fatal on entry is untouched bitwise (state, fluxes, bank).  A
   !! COUNTED bit (`ICE_ERR_REMAP_FALLBACK`) does not stop the cell: the clean-up rebin is exactly the
   !! repair for it.  Stage A flags `ICE_ERR_NONFINITE` and leaves the cell untouched when `frain`,
   !! `aice`, `rside`, `frzmlt` or `Tf` is not finite, or the area sum is not; a non-finite or
   !! non-positive `dt` is refused ONCE in the host shim (every owned cell is flagged, nothing else is
   !! done).  The remap works on cell-local aggregates: the stored `aice`, `aice0` are not written
   !! until the clean-up of stage D (the new-ice kernel does write `aice0(i, j)`, which D overwrites).
   !!
   !! **Exact early-out.**  Stage C skips a cell whose every category area, ice volume and snow volume
   !! is exactly 0 (lateral melt would add `+0` to the fluxes and scale zeros: bitwise the same); the
   !! clean-up has its own exact empty-cell shortcut.
   !!
   !! **Why staged, not one fused kernel.**  A fused per-cell kernel (all four stages in one thread)
   !! was built first and measured on a V100: 228 registers and 3144 B of stack per thread (the remap
   !! alone 193 and 2456 B), and 12-20 % SLOWER than the staged sweeps (32.3 vs 26.8 ms on a random
   !! 512 x 512 x 5 x 7 state, 4.75 vs 4.25 ms on a coherent ~85 % ice-free one).  The staged form
   !! keeps each kernel at its own footprint.  Launches-per-step timings are in the test.
   !!
   !! **Dropped upstream branches** (each inert for this port): the single-category `reduce_area`
   !! (configure refuses `ncat < 2`); the floe-size-distribution welding and lateral-growth tracers;
   !! the melt-pond water flux (`fpond`, `dpnd_melt`); aerosol, isotope and biogeochemistry fluxes;
   !! `frz_onset`; `first_ice`.  `wlat` is not an input: upstream reads it only for the floe-size
   !! distribution; the melted fraction `rside` carries everything the lateral melt needs.
   !!
   !! **Caller contract** (the step-sequencing PR).  `ice_flux_zero` and `ice_save_init` first; therm1
   !! (vertical thermodynamics, its merge into `fresh/fsalt/fhocn`, and the basal-heat routine that
   !! fills `rside`); then this.  It supplies, per cell: `frzmlt` (signed, above), `bank` (J m-2 of
   !! cell), `Tf` (the ICE-side linear freezing temperature `-depressT*sss`, not the ocean
   !! equation-of-state freezing point), `rside` (0 where therm1 gated the lateral melt off), `frain`
   !! (kg m-2 s-1, per cell area), `wet`, `dt`, `hin_max(0:ncat)` and `salinz(nilyr)`.  New ice spends
   !! `frzmlt*dt`: the caller forms `frzmlt = bank/dt` with the SAME `dt` (the volume then equals the
   !! bank's to about 1 ulp, not bitwise).  The ocean must be debited ONLY through `fhocn` (never by
   !! `-frzmlt`: the bank already paid the new ice).  A positive bank with `frzmlt <= 0` is KEPT.
   !! A later PR (ridging, transport) aggregates again from the state it changes; the `aice, vice,
   !! vsno, aice0` left here are exact for the state left here.
   use, intrinsic :: iso_fortran_env, only: int32
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_constants, only: wp
   use rdb_ice_errors, only: ICE_ERR_NONFINITE, ICE_ERR_FATAL_MASK
   use rdb_ice_params, only: ice_params_t
   use rdb_ice_itd_linear, only: ice_linear_itd_cell
   use rdb_ice_new_ice, only: ice_add_new_ice_cell
   use rdb_ice_lateral_melt, only: ice_lateral_melt_cell
   use rdb_ice_itd_cleanup, only: ice_cleanup_itd_cell
   implicit none
   private

   public :: ice_step_therm2
   public :: ice_therm2_remap_cell, ice_therm2_newice_cell, ice_therm2_lateral_cell, ice_therm2_cleanup_cell

contains

   pure subroutine ice_therm2_remap_cell(i, j, nx, ny, ncat, nilyr, nslyr, do_linear_itd, rhos, Lfresh, hi_min, puny, &
                                         hin_max, Tf, frzmlt, frain, rside, aice, aicen_init, vicen_init, &
                                         aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, fresh, ierr)
      !! Stage A for cell `(i, j)`: finite check of `aice, frain, rside, frzmlt, Tf` and of the area sum
      !! (`ICE_ERR_NONFINITE`, cell untouched), the rain term `fresh += frain*aice` with the ENTRY
      !! aggregate `aice(i, j)` (read only), then the linear remap where `do_linear_itd` and the area
      !! sum exceeds `puny`.  Skipped for a cell with a fatal bit.
      !$acc routine seq
      integer, value :: i, j
      integer, intent(in) :: nx, ny, ncat, nilyr, nslyr
      logical, intent(in) :: do_linear_itd
      real(wp), intent(in) :: rhos, Lfresh, hi_min, puny
      real(wp), intent(in) :: hin_max(0:ncat)
      real(wp), intent(in) :: Tf, frzmlt, frain, rside
      real(wp), intent(in) :: aice(nx, ny), aicen_init(nx, ny, ncat), vicen_init(nx, ny, ncat)
      real(wp), intent(inout) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat)
      real(wp), intent(inout) :: Tsfcn(nx, ny, ncat), qicen(nx, ny, ncat, nilyr)
      real(wp), intent(inout) :: sicen(nx, ny, ncat, nilyr), qsnon(nx, ny, ncat, nslyr)
      real(wp), intent(inout) :: fresh(nx, ny)
      integer(int32), intent(inout) :: ierr(nx, ny)

      integer :: n
      real(wp) :: a_in, tot, aice_l, aice0_l

      if (iand(ierr(i, j), ICE_ERR_FATAL_MASK) /= 0_int32) return
      a_in = aice(i, j)
      if (.not. (ieee_is_finite(a_in) .and. ieee_is_finite(frain) .and. ieee_is_finite(rside) .and. &
                 ieee_is_finite(frzmlt) .and. ieee_is_finite(Tf))) then
         ierr(i, j) = ior(ierr(i, j), ICE_ERR_NONFINITE)
         return
      end if
      fresh(i, j) = fresh(i, j) + frain*a_in
      tot = 0.0_wp
      do n = 1, ncat
         tot = tot + aicen(i, j, n)
      end do
      if (.not. ieee_is_finite(tot)) then
         ierr(i, j) = ior(ierr(i, j), ICE_ERR_NONFINITE)
         return
      end if
      if (do_linear_itd .and. tot > puny) then
         call ice_linear_itd_cell(i, j, nx, ny, ncat, nilyr, nslyr, hin_max, aicen_init, vicen_init, &
                                  aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, &
                                  Tf, rhos, Lfresh, hi_min, puny, .false., aice_l, aice0_l, ierr)
      end if
   end subroutine ice_therm2_remap_cell

   pure subroutine ice_therm2_newice_cell(i, j, nx, ny, ncat, nilyr, dt, rhoi, Lfresh, puny, sref, hfrazilmin, &
                                          update_ocn_f, hin_max, salinz, Tf, frzmlt, bank, aicen, vicen, Tsfcn, qicen, &
                                          sicen, aice0, fresh, fsalt, frazil_vol, ierr)
      !! Stage B: new ice where `frzmlt > 0` (the bank is reset only when no fatal bit is left), skipped
      !! for a cell with a fatal bit.
      !$acc routine seq
      integer, value :: i, j
      integer, intent(in) :: nx, ny, ncat, nilyr
      real(wp), intent(in) :: dt, rhoi, Lfresh, puny, sref, hfrazilmin
      logical, intent(in) :: update_ocn_f
      real(wp), intent(in) :: hin_max(0:ncat), salinz(nilyr)
      real(wp), intent(in) :: Tf, frzmlt
      real(wp), intent(inout) :: bank(nx, ny)
      real(wp), intent(inout) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), Tsfcn(nx, ny, ncat)
      real(wp), intent(inout) :: qicen(nx, ny, ncat, nilyr), sicen(nx, ny, ncat, nilyr)
      real(wp), intent(inout) :: aice0(nx, ny), fresh(nx, ny), fsalt(nx, ny), frazil_vol(nx, ny)
      integer(int32), intent(inout) :: ierr(nx, ny)

      if (iand(ierr(i, j), ICE_ERR_FATAL_MASK) /= 0_int32) return
      if (frzmlt > 0.0_wp) then
         call ice_add_new_ice_cell(i, j, nx, ny, ncat, nilyr, dt, frzmlt, Tf, hin_max(1), rhoi, Lfresh, puny, &
                                   hfrazilmin, sref, update_ocn_f, salinz, aicen, vicen, Tsfcn, qicen, sicen, &
                                   aice0, fresh, fsalt, frazil_vol, ierr)
         if (iand(ierr(i, j), ICE_ERR_FATAL_MASK) == 0_int32) bank(i, j) = 0.0_wp
      end if
   end subroutine ice_therm2_newice_cell

   pure subroutine ice_therm2_lateral_cell(i, j, nx, ny, ncat, nilyr, nslyr, dt, rhoi, rhos, sref, puny, rside, &
                                           aicen, vicen, vsnon, qicen, qsnon, fresh, fsalt, fhocn, meltl, ierr)
      !! Stage C: lateral melt where `rside /= 0` and the cell holds anything, skipped for a fatal cell.
      !$acc routine seq
      integer, value :: i, j
      integer, intent(in) :: nx, ny, ncat, nilyr, nslyr
      real(wp), intent(in) :: dt, rhoi, rhos, sref, puny, rside
      real(wp), intent(inout) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat)
      real(wp), intent(in) :: qicen(nx, ny, ncat, nilyr), qsnon(nx, ny, ncat, nslyr)
      real(wp), intent(inout) :: fresh(nx, ny), fsalt(nx, ny), fhocn(nx, ny), meltl(nx, ny)
      integer(int32), intent(inout) :: ierr(nx, ny)

      integer :: n
      logical :: empty

      if (iand(ierr(i, j), ICE_ERR_FATAL_MASK) /= 0_int32) return
      if (rside == 0.0_wp) return
      empty = .true.
      do n = 1, ncat
         if (.not. (aicen(i, j, n) == 0.0_wp .and. vicen(i, j, n) == 0.0_wp .and. vsnon(i, j, n) == 0.0_wp)) empty = .false.
      end do
      if (empty .and. rside >= 0.0_wp .and. rside <= 1.0_wp) return      ! bitwise a no-op; a bad rside still flags
      call ice_lateral_melt_cell(i, j, nx, ny, ncat, nilyr, nslyr, rside, dt, rhoi, rhos, sref, puny, &
                                 aicen, vicen, vsnon, qicen, qsnon, fresh, fsalt, fhocn, meltl, ierr)
   end subroutine ice_therm2_lateral_cell

   pure subroutine ice_therm2_cleanup_cell(i, j, nx, ny, ncat, nilyr, nslyr, hin_max, dt, Tf, rhoi, rhos, Lfresh, cp_ice, &
                                           hs_min, puny, area_min, mass_min, sref, aicen, vicen, vsnon, Tsfcn, qicen, &
                                           sicen, qsnon, aice, vice, vsno, aice0, Tsfc_c, qice_c, sice_c, qsno_c, &
                                           fresh, fsalt, fhocn, ierr)
      !! Stage D: the clean-up (rebin on, `limit_aice` on) and the one aggregate of the step, skipped for a
      !! fatal cell.
      !$acc routine seq
      integer, value :: i, j
      integer, intent(in) :: nx, ny, ncat, nilyr, nslyr
      real(wp), intent(in) :: hin_max(0:ncat)
      real(wp), intent(in) :: dt, Tf, rhoi, rhos, Lfresh, cp_ice, hs_min, puny, area_min, mass_min, sref
      real(wp), intent(inout) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat)
      real(wp), intent(inout) :: Tsfcn(nx, ny, ncat), qicen(nx, ny, ncat, nilyr)
      real(wp), intent(inout) :: sicen(nx, ny, ncat, nilyr), qsnon(nx, ny, ncat, nslyr)
      real(wp), intent(inout) :: aice(nx, ny), vice(nx, ny), vsno(nx, ny), aice0(nx, ny)
      real(wp), intent(inout) :: Tsfc_c(nx, ny), qice_c(nx, ny, nilyr), sice_c(nx, ny, nilyr), qsno_c(nx, ny, nslyr)
      real(wp), intent(inout) :: fresh(nx, ny), fsalt(nx, ny), fhocn(nx, ny)
      integer(int32), intent(inout) :: ierr(nx, ny)

      if (iand(ierr(i, j), ICE_ERR_FATAL_MASK) /= 0_int32) return
      call ice_cleanup_itd_cell(i, j, nx, ny, ncat, nilyr, nslyr, hin_max, dt, Tf, rhoi, rhos, Lfresh, cp_ice, &
                                hs_min, puny, area_min, mass_min, sref, .true., .true., .true., &
                                aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, &
                                aice, vice, vsno, aice0, Tsfc_c, qice_c, sice_c, qsno_c, &
                                fresh, fsalt, fhocn, ierr)
   end subroutine ice_therm2_cleanup_cell

   subroutine ice_step_therm2(nx, ny, ncat, nilyr, nslyr, i_lo, i_hi, j_lo, j_hi, dt, p, hin_max, salinz, &
                              wet, Tf, frzmlt, frain, rside, bank, &
                              aicen_init, vicen_init, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, &
                              aice, vice, vsno, aice0, Tsfc_c, qice_c, sice_c, qsno_c, &
                              fresh, fsalt, fhocn, meltl, frazil_vol, ierr, fuse_bc)
      !! Host shim of the `therm2` sweeps over the owned cells `i_lo:i_hi, j_lo:j_hi` (see the module
      !! header).  Every array, `hin_max`, `salinz` and `wet` included, must be device-present on an
      !! offload build; the parameter bundle is read on the host.  The remap runs when `p%itd%kitd == 1`.
      !! `fuse_bc` (default `.true.`, measured marginally faster): run the new ice and the lateral melt in ONE
      !! sweep (three launches; `.false.` gives four).  A non-finite or non-positive `dt` flags every owned cell and does nothing else.
      integer, intent(in) :: nx, ny, ncat, nilyr, nslyr, i_lo, i_hi, j_lo, j_hi
      real(wp), intent(in) :: dt
      type(ice_params_t), intent(in) :: p
      real(wp), intent(in) :: hin_max(0:ncat), salinz(nilyr)
      logical, intent(in) :: wet(nx, ny)
      real(wp), intent(in) :: Tf(nx, ny), frzmlt(nx, ny), frain(nx, ny), rside(nx, ny)
      real(wp), intent(inout) :: bank(nx, ny)
      real(wp), intent(in) :: aicen_init(nx, ny, ncat), vicen_init(nx, ny, ncat)
      real(wp), intent(inout) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat)
      real(wp), intent(inout) :: Tsfcn(nx, ny, ncat), qicen(nx, ny, ncat, nilyr)
      real(wp), intent(inout) :: sicen(nx, ny, ncat, nilyr), qsnon(nx, ny, ncat, nslyr)
      real(wp), intent(inout) :: aice(nx, ny), vice(nx, ny), vsno(nx, ny), aice0(nx, ny)
      real(wp), intent(inout) :: Tsfc_c(nx, ny), qice_c(nx, ny, nilyr), sice_c(nx, ny, nilyr), qsno_c(nx, ny, nslyr)
      real(wp), intent(inout) :: fresh(nx, ny), fsalt(nx, ny), fhocn(nx, ny), meltl(nx, ny), frazil_vol(nx, ny)
      integer(int32), intent(inout) :: ierr(nx, ny)
      logical, intent(in), optional :: fuse_bc
      logical :: fuse

      fuse = .true.
      if (present(fuse_bc)) fuse = fuse_bc
      if (.not. (ieee_is_finite(dt) .and. dt > 0.0_wp)) then
         call flag_all(nx, ny, i_lo, i_hi, j_lo, j_hi, ierr)
         return
      end if
      call therm2_impl(nx, ny, ncat, nilyr, nslyr, i_lo, i_hi, j_lo, j_hi, dt, p%itd%kitd == 1, fuse, &
                       p%c%rhoi, p%c%rhos, p%c%Lfresh, p%c%cp_ice, p%c%hs_min, p%c%puny, p%c%itd_area_min, &
                       p%c%itd_mass_min, p%thermo%ice_ref_salinity, p%thermo%hfrazilmin, p%thermo%hi_min, &
                       p%thermo%update_ocn_f, hin_max, salinz, wet, Tf, frzmlt, frain, rside, bank, &
                       aicen_init, vicen_init, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, &
                       aice, vice, vsno, aice0, Tsfc_c, qice_c, sice_c, qsno_c, &
                       fresh, fsalt, fhocn, meltl, frazil_vol, ierr)
   end subroutine ice_step_therm2

   subroutine flag_all(nx, ny, i_lo, i_hi, j_lo, j_hi, ierr)
      integer, intent(in) :: nx, ny, i_lo, i_hi, j_lo, j_hi
      integer(int32), intent(inout) :: ierr(nx, ny)
      integer :: i, j

      !$acc data present(ierr)
      do concurrent(j=j_lo:j_hi, i=i_lo:i_hi)
         ierr(i, j) = ior(ierr(i, j), ICE_ERR_NONFINITE)
      end do
      !$acc end data
   end subroutine flag_all

   subroutine therm2_impl(nx, ny, ncat, nilyr, nslyr, i_lo, i_hi, j_lo, j_hi, dt, do_linear_itd, fuse_bc, &
                          rhoi, rhos, Lfresh, cp_ice, hs_min, puny, area_min, mass_min, sref, &
                          hfrazilmin, hi_min, update_ocn_f, hin_max, salinz, wet, Tf, frzmlt, frain, rside, bank, &
                          aicen_init, vicen_init, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, &
                          aice, vice, vsno, aice0, Tsfc_c, qice_c, sice_c, qsno_c, &
                          fresh, fsalt, fhocn, meltl, frazil_vol, ierr)
      integer, intent(in) :: nx, ny, ncat, nilyr, nslyr, i_lo, i_hi, j_lo, j_hi
      real(wp), intent(in) :: dt
      logical, intent(in) :: do_linear_itd, fuse_bc, update_ocn_f
      real(wp), intent(in) :: rhoi, rhos, Lfresh, cp_ice, hs_min, puny, area_min, mass_min, sref, hfrazilmin, hi_min
      real(wp), intent(in) :: hin_max(0:ncat), salinz(nilyr)
      logical, intent(in) :: wet(nx, ny)
      real(wp), intent(in) :: Tf(nx, ny), frzmlt(nx, ny), frain(nx, ny), rside(nx, ny)
      real(wp), intent(inout) :: bank(nx, ny)
      real(wp), intent(in) :: aicen_init(nx, ny, ncat), vicen_init(nx, ny, ncat)
      real(wp), intent(inout) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat)
      real(wp), intent(inout) :: Tsfcn(nx, ny, ncat), qicen(nx, ny, ncat, nilyr)
      real(wp), intent(inout) :: sicen(nx, ny, ncat, nilyr), qsnon(nx, ny, ncat, nslyr)
      real(wp), intent(inout) :: aice(nx, ny), vice(nx, ny), vsno(nx, ny), aice0(nx, ny)
      real(wp), intent(inout) :: Tsfc_c(nx, ny), qice_c(nx, ny, nilyr), sice_c(nx, ny, nilyr), qsno_c(nx, ny, nslyr)
      real(wp), intent(inout) :: fresh(nx, ny), fsalt(nx, ny), fhocn(nx, ny), meltl(nx, ny), frazil_vol(nx, ny)
      integer(int32), intent(inout) :: ierr(nx, ny)
      integer :: i, j

      !$acc data present(hin_max, salinz, wet, Tf, frzmlt, frain, rside, bank, aicen_init, vicen_init, aicen, vicen, vsnon) &
      !$acc& present(Tsfcn, qicen, sicen, qsnon, aice, vice, vsno, aice0, Tsfc_c, qice_c, sice_c, qsno_c) &
      !$acc& present(fresh, fsalt, fhocn, meltl, frazil_vol, ierr)
      ! A: finite check, rain, remap
      do concurrent(j=j_lo:j_hi, i=i_lo:i_hi)
         if (wet(i, j)) then
            call ice_therm2_remap_cell(i, j, nx, ny, ncat, nilyr, nslyr, do_linear_itd, rhos, Lfresh, hi_min, puny, &
                                       hin_max, Tf(i, j), frzmlt(i, j), frain(i, j), rside(i, j), aice, aicen_init, &
                                       vicen_init, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, fresh, ierr)
         end if
      end do
      if (fuse_bc) then
         ! B + C in one sweep
         do concurrent(j=j_lo:j_hi, i=i_lo:i_hi)
            if (wet(i, j)) then
               call ice_therm2_newice_cell(i, j, nx, ny, ncat, nilyr, dt, rhoi, Lfresh, puny, sref, hfrazilmin, &
                                           update_ocn_f, hin_max, salinz, Tf(i, j), frzmlt(i, j), bank, aicen, vicen, &
                                           Tsfcn, qicen, sicen, aice0, fresh, fsalt, frazil_vol, ierr)
               call ice_therm2_lateral_cell(i, j, nx, ny, ncat, nilyr, nslyr, dt, rhoi, rhos, sref, puny, rside(i, j), &
                                            aicen, vicen, vsnon, qicen, qsnon, fresh, fsalt, fhocn, meltl, ierr)
            end if
         end do
      else
         ! B: new ice
         do concurrent(j=j_lo:j_hi, i=i_lo:i_hi)
            if (wet(i, j)) then
               call ice_therm2_newice_cell(i, j, nx, ny, ncat, nilyr, dt, rhoi, Lfresh, puny, sref, hfrazilmin, &
                                           update_ocn_f, hin_max, salinz, Tf(i, j), frzmlt(i, j), bank, aicen, vicen, &
                                           Tsfcn, qicen, sicen, aice0, fresh, fsalt, frazil_vol, ierr)
            end if
         end do
         ! C: lateral melt
         do concurrent(j=j_lo:j_hi, i=i_lo:i_hi)
            if (wet(i, j)) then
               call ice_therm2_lateral_cell(i, j, nx, ny, ncat, nilyr, nslyr, dt, rhoi, rhos, sref, puny, rside(i, j), &
                                            aicen, vicen, vsnon, qicen, qsnon, fresh, fsalt, fhocn, meltl, ierr)
            end if
         end do
      end if
      ! D: clean-up and the aggregates
      do concurrent(j=j_lo:j_hi, i=i_lo:i_hi)
         if (wet(i, j)) then
            call ice_therm2_cleanup_cell(i, j, nx, ny, ncat, nilyr, nslyr, hin_max, dt, Tf(i, j), rhoi, rhos, Lfresh, &
                                         cp_ice, hs_min, puny, area_min, mass_min, sref, aicen, vicen, vsnon, Tsfcn, &
                                         qicen, sicen, qsnon, aice, vice, vsno, aice0, Tsfc_c, qice_c, sice_c, qsno_c, &
                                         fresh, fsalt, fhocn, ierr)
         end if
      end do
      !$acc end data
   end subroutine therm2_impl

end module rdb_ice_therm2
