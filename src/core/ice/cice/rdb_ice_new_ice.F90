!! New (frazil) sea ice formed in open water from the ocean's freezing
!! potential, partitioned between category 1 and the existing ice.
module rdb_ice_new_ice
   ! SPDX-License-Identifier: BSD-3-Clause
   !
   ! Ported from CICE/Icepack (CICE Consortium, https://github.com/CICE-Consortium).
   ! Upstream: Icepack 9a1f3e4 (Icepack 1.5.3+14), module icepack_therm_itd
   ! (add_new_ice, the Bitz-Lipscomb branch).
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
   ! Original authors (from the upstream headers): William H. Lipscomb,
   ! Elizabeth C. Hunke (LANL), C. M. Bitz (UW).
   ! Physics: Bitz and Lipscomb (1999, JGR 104 15669); Lipscomb (2001, JGR 106).
   !
   !! **Physics.**  The ocean has already booked the latent heat of the frazil
   !! it formed (the supercooled water was warmed back to freezing); what it
   !! returns to the ice is the freezing potential `frzmlt` [W m-2, > 0].  The
   !! ice turns it into a volume of new ice per cell area
   !! `V = frzmlt*dt/(rhoi*Lfresh)` (the new-ice enthalpy is the pure latent
   !! content `qi0new = -rhoi*Lfresh`, Bitz and Lipscomb 1999), and NEVER
   !! writes `fhocn`: the heat is not charged to the ocean a second time.
   !!
   !! **Partition** (`hi0max = 0.9*hin_max(1)`, kept off the category bound;
   !! `aice0` open water, `aice = sum(aicen)`):
   !! - `aice0 > puny`: `hi0new = max(V/aice0, hfrazilmin)`.  If `hi0new > hi0max`
   !!   and `aice0 + puny < 1` (some ice is present) the new ice fills ALL the
   !!   open water at `hi0new = hi0max` and the surplus volume is spread over
   !!   the existing categories as a uniform thickness `hsurp = vsurp/aice`
   !!   (volume `hsurp*aicen(n)` into category `n`, area unchanged).  Otherwise
   !!   the new ice covers `V/hi0new` of the cell (a thin frazil covers less
   !!   than the open water: `hfrazilmin` bounds the new-ice AREA from above,
   !!   never the volume; with no ice present the new ice may exceed `hi0max`
   !!   and sits in category 1 above its bound until the next remap).
   !! - `aice0 <= puny` (closed ice): the whole volume is spread, `hsurp = V/aice`.
   !! New area and volume go into category 1 only; surface temperature becomes
   !! the area-weighted mean with the freezing temperature `Tf`, capped at 0;
   !! layer enthalpy and salinity are volume-weighted means of the old layer
   !! and the new ice (`qi0new`, the BL99 salinity profile `salinz(k)`).  The
   !! spread step runs FIRST, then the category-1 merge sees the spread volume.
   !!
   !! **Fluxes** (per unit CELL area and second, every producer ADDS; signs are
   !! Icepack's).  With `update_ocn_f`: `fresh += -rhoi*V/dt` (< 0),
   !! `fsalt += ice_ref_salinity*1e-3*(-rhoi*V/dt)` (< 0): the ocean is charged
   !! the mass and the reference-salinity salt of the ice it lost.  Without it
   !! both increments are exactly zero.  NOTHING in this model adds the water
   !! and salt back: with `.false.` the ice gains mass and salt the ocean never
   !! loses.  That is legal but NOT conservative until the coupling PR handles
   !! it.  The salt is a flux convention, not a closure: the ice stores
   !! `salinz`, the ocean is charged `ice_ref_salinity`.  `frazil_vol += V`
   !! (the volume before the partition; closes to the bank to about 1 ulp,
   !! not bitwise, since the caller forms `frzmlt = bank/dt`; upstream assigns it, equal to `+=` after
   !! the per-step zero of the flux slot, and ONE call per step is the contract).
   !!
   !! **Dropped upstream branches** (each inert for this port): the mushy
   !! thermodynamics (`ktherm = 2`: frazil salinity from `sss`, the
   !! `fresh_ice_correction` of `cpl_frazil`, `update_vertical_tracers`); the
   !! `cpl_frazil = "external"` bypass and `saltflux_option = "prognostic"`
   !! (configure refuses both); the floe-size distribution (`fsd_lateral_growth`,
   !! `fsd_add_new_ice`, new ice into every category); the tracers iage, level
   !! ice, first-year, aerosols, isotopes and melt ponds; the biogeochemistry
   !! call; `frz_onset`; the optional `conserv_check` (it becomes the closure
   !! tests of the suite).  The `ncat = 1` branch (`hi0max = bignum`) is kept but dead: configure
   !! refuses `ncat < 2`.
   !!
   !! **Bank (roundabout addition, not upstream).**  The ocean-side frazil
   !! clamp accumulates its heat in a bank [J m-2 of cell].  `ice_new_ice_step`
   !! spends it: where the caller's `frzmlt > 0` and the cell is wet and the
   !! kernel raises no fatal bit, the bank is reset to 0 (all of it became ice,
   !! `hfrazilmin` never discards volume); a dry cell, or a cell with a fatal
   !! bit, keeps its bank (conserved), and a `frzmlt <= 0` cell keeps it too.
   !! The caller forms `frzmlt` (`bank/dt` where banked, else the melt side).
   !!
   !! **GPU shape.**  The kernel couples the categories of ONE cell: it works
   !! IN PLACE on the SoA state through `(i, j, n, k)` indexing (no gathered
   !! column, no local array at all; `i`, `j` by value).  Category 1 is always
   !! touched; the other categories only when surplus volume is spread.  The
   !! driver is ONE `do concurrent (j, i)`.  Layer index `k = 1` is the ice TOP.
   !!
   !! **Errors (`ierr` bits, OR-ed in, never an abort).**  For `frzmlt <= 0`
   !! nothing is read or written (bitwise no-op; a NaN `frzmlt` sets
   !! `ICE_ERR_NONFINITE`).  Otherwise: a non-finite `Tf`, `aicen`, `vicen`,
   !! `Tsfcn(1)`, `qicen` or `sicen` sets `ICE_ERR_NONFINITE`; an
   !! `aicen < -puny` at entry, or an exit open water `< -puny`, sets
   !! `ICE_ERR_NEW_ICE_AREA`.  Entry failures leave the cell untouched;
   !! the exit failure leaves the outputs UNSPECIFIED (the cell is fatal).
   !! `aice0` is recomputed from `aicen`, never read.
   use, intrinsic :: iso_fortran_env, only: int32
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_constants, only: wp
   use rdb_ice_errors, only: ICE_ERR_NONFINITE, ICE_ERR_NEW_ICE_AREA, ICE_ERR_FATAL_MASK
   use rdb_ice_params, only: ice_params_t
   implicit none
   private

   public :: ice_add_new_ice_cell, ice_new_ice_step
   public :: NEW_ICE_HI0MAX_FRAC

   real(wp), parameter :: NEW_ICE_HI0MAX_FRAC = 0.9_wp
      !! `hi0max = NEW_ICE_HI0MAX_FRAC*hin_max(1)` (not too close to the bound).
   real(wp), parameter :: NEW_ICE_BIGNUM = 1.0e30_wp
      !! `hi0max` for a single category (upstream `bignum`).

contains

   pure subroutine ice_add_new_ice_cell(i, j, nx, ny, ncat, nilyr, dt, frzmlt, Tf, hin_max1, &
                                        rhoi, Lfresh, puny, hfrazilmin, ice_ref_salinity, update_ocn_f, &
                                        salinz, aicen, vicen, Tsfcn, qicen, sicen, aice0, &
                                        fresh, fsalt, frazil_vol, ierr)
      !! Add the new ice of freezing potential `frzmlt` [W m-2] to cell `(i, j)`
      !! (see the module header).  `hin_max1 = hin_max(1)`; `salinz(nilyr)` is
      !! the BL99 salinity profile of new ice.  `Tf` is the ICE-side freezing
      !! temperature (the linear `-depressT*sss` form), never the ocean
      !! equation-of-state freezing point.  `aice = sum aicen` (upstream's order)
      !! and the open water `max(1 - aice, 0)` are RECOMPUTED here: the incoming
      !! `aice0(i, j)` is never read (a stale value cannot overfill the cell);
      !! it is an OUTPUT, written (the post-fill open water) whenever
      !! `frzmlt > 0` and no entry error occurs, and untouched otherwise.
      !! `aicen, vicen, Tsfcn, qicen, sicen` are updated in place,
      !! `fresh, fsalt, frazil_vol` are ADDED to.
      !$acc routine seq
      integer, value :: i, j
      integer, intent(in) :: nx, ny, ncat, nilyr
      real(wp), intent(in) :: dt, frzmlt, Tf, hin_max1, rhoi, Lfresh, puny, hfrazilmin, ice_ref_salinity
      logical, intent(in) :: update_ocn_f
      real(wp), intent(in) :: salinz(nilyr)
      real(wp), intent(inout) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), Tsfcn(nx, ny, ncat)
      real(wp), intent(inout) :: qicen(nx, ny, ncat, nilyr), sicen(nx, ny, ncat, nilyr)
      real(wp), intent(inout) :: aice0(nx, ny), fresh(nx, ny), fsalt(nx, ny), frazil_vol(nx, ny)
      integer(int32), intent(inout) :: ierr(nx, ny)
      integer :: n, k
      logical :: bad
      real(wp) :: aice, hi0max, qi0new, vi0new, ai0mod, hi0new, ai0new, hsurp, vsurp, vtmp
      real(wp) :: dfresh, dfsalt, area1, vice1, dan, vin1, a0, vfrz

      if (.not. ieee_is_finite(frzmlt)) then
         ierr(i, j) = ior(ierr(i, j), ICE_ERR_NONFINITE)
         return
      end if
      if (.not. (frzmlt > 0.0_wp)) return          ! fnew = 0: bitwise no-op

      ! Entry guard: every input read is checked before anything is written.
      bad = .not. (ieee_is_finite(Tf) .and. ieee_is_finite(Tsfcn(i, j, 1)))
      do n = 1, ncat
         if (.not. (ieee_is_finite(aicen(i, j, n)) .and. ieee_is_finite(vicen(i, j, n)))) bad = .true.
         do k = 1, nilyr
            if (.not. (ieee_is_finite(qicen(i, j, n, k)) .and. ieee_is_finite(sicen(i, j, n, k)))) bad = .true.
         end do
      end do
      if (bad) then
         ierr(i, j) = ior(ierr(i, j), ICE_ERR_NONFINITE)
         return
      end if
      bad = .false.
      aice = 0.0_wp
      do n = 1, ncat
         if (aicen(i, j, n) < -puny) bad = .true.
         aice = aice + aicen(i, j, n)
      end do
      if (bad) then
         ierr(i, j) = ior(ierr(i, j), ICE_ERR_NEW_ICE_AREA)
         return
      end if
      a0 = max(1.0_wp - aice, 0.0_wp)

      if (ncat > 1) then
         hi0max = hin_max1*NEW_ICE_HI0MAX_FRAC
      else
         hi0max = NEW_ICE_BIGNUM
      end if
      qi0new = -rhoi*Lfresh
      vi0new = -frzmlt*dt/qi0new                    ! qi0new < 0
      vfrz = vi0new                                 ! volume before the partition

      ! Partition of the volume (scalars only; nothing is written yet).
      hsurp = 0.0_wp
      hi0new = 0.0_wp
      ai0new = 0.0_wp
      vin1 = 0.0_wp
      dan = 0.0_wp
      if (vi0new > 0.0_wp) then
         ai0mod = a0
         if (ai0mod > puny) then
            hi0new = max(vi0new/ai0mod, hfrazilmin)
            if (hi0new > hi0max .and. ai0mod + puny < 1.0_wp) then
               hi0new = hi0max
               ai0new = ai0mod
               vsurp = vi0new - ai0new*hi0new
               hsurp = vsurp/aice
               vi0new = ai0new*hi0new
            else
               ai0new = vi0new/hi0new
            end if
         else                                        ! closed ice: spread everything
            hsurp = vi0new/aice
            vi0new = 0.0_wp
         end if
         dan = ai0new
         vin1 = vi0new
      end if

      ! Fluxes to the ocean: mass and salt only, never heat.
      frazil_vol(i, j) = frazil_vol(i, j) + vfrz
      dfresh = 0.0_wp
      if (update_ocn_f) dfresh = -rhoi*vfrz/dt
      dfsalt = ice_ref_salinity*0.001_wp*dfresh
      fresh(i, j) = fresh(i, j) + dfresh
      fsalt(i, j) = fsalt(i, j) + dfsalt

      if (hsurp > 0.0_wp) then                       ! spread over every category
         do n = 1, ncat
            vsurp = hsurp*aicen(i, j, n)
            vtmp = vicen(i, j, n) + vsurp
            vicen(i, j, n) = vtmp
            do k = 1, nilyr
               vsurp = hsurp*aicen(i, j, n)
               vtmp = vicen(i, j, n) - vsurp
               if (vicen(i, j, n) > 0.0_wp) then
                  qicen(i, j, n, k) = (qicen(i, j, n, k)*vtmp + qi0new*vsurp)/vicen(i, j, n)
                  sicen(i, j, n, k) = (sicen(i, j, n, k)*vtmp + salinz(k)*vsurp)/vicen(i, j, n)
               end if
            end do
         end do
      end if

      if (dan > 0.0_wp .and. vin1 > 0.0_wp) then      ! new ice into category 1
         area1 = aicen(i, j, 1)
         vice1 = vicen(i, j, 1)
         aicen(i, j, 1) = aicen(i, j, 1) + dan
         a0 = a0 - dan
         vicen(i, j, 1) = vicen(i, j, 1) + vin1
         Tsfcn(i, j, 1) = (Tsfcn(i, j, 1)*area1 + Tf*dan)/aicen(i, j, 1)
         Tsfcn(i, j, 1) = min(Tsfcn(i, j, 1), 0.0_wp)
         do k = 1, nilyr
            if (vicen(i, j, 1) > 0.0_wp) then
               qicen(i, j, 1, k) = (qicen(i, j, 1, k)*vice1 + qi0new*vin1)/vicen(i, j, 1)
               sicen(i, j, 1, k) = (sicen(i, j, 1, k)*vice1 + salinz(k)*vin1)/vicen(i, j, 1)
            end if
         end do
      end if
      aice0(i, j) = a0
      if (a0 < -puny) ierr(i, j) = ior(ierr(i, j), ICE_ERR_NEW_ICE_AREA)
   end subroutine ice_add_new_ice_cell

   subroutine ice_new_ice_step(nx, ny, ncat, nilyr, i_lo, i_hi, j_lo, j_hi, dt, p, hin_max, salinz, &
                               wet, frzmlt, bank, Tf, aicen, vicen, Tsfcn, qicen, sicen, aice0, &
                               fresh, fsalt, frazil_vol, ierr)
      !! Step driver over the owned cells `i_lo:i_hi, j_lo:j_hi`: new ice where
      !! `wet` and `frzmlt > 0`, then the bank reset described in the module
      !! header (a cell already carrying a fatal bit is not processed and keeps its bank;
      !! the reset is gated on the WHOLE error word, so a fatal bit raised now keeps it too).
      !! `hin_max` is read on the HOST (static after setup); every other
      !! array, `salinz` and `wet` included, must be device-present on an offload
      !! build.  `p` is read on the host (scalars go into the kernel).
      integer, intent(in) :: nx, ny, ncat, nilyr, i_lo, i_hi, j_lo, j_hi
      real(wp), intent(in) :: dt
      type(ice_params_t), intent(in) :: p
      real(wp), intent(in) :: hin_max(0:ncat), salinz(nilyr)
      logical, intent(in) :: wet(nx, ny)
      real(wp), intent(in) :: frzmlt(nx, ny), Tf(nx, ny)
      real(wp), intent(inout) :: bank(nx, ny)
      real(wp), intent(inout) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), Tsfcn(nx, ny, ncat)
      real(wp), intent(inout) :: qicen(nx, ny, ncat, nilyr), sicen(nx, ny, ncat, nilyr)
      real(wp), intent(inout) :: aice0(nx, ny), fresh(nx, ny), fsalt(nx, ny), frazil_vol(nx, ny)
      integer(int32), intent(inout) :: ierr(nx, ny)

      call new_ice_impl(nx, ny, ncat, nilyr, i_lo, i_hi, j_lo, j_hi, dt, p%c%rhoi, p%c%Lfresh, p%c%puny, &
                        p%thermo%hfrazilmin, p%thermo%ice_ref_salinity, p%thermo%update_ocn_f, hin_max(1), &
                        salinz, wet, frzmlt, bank, Tf, aicen, vicen, Tsfcn, qicen, sicen, aice0, &
                        fresh, fsalt, frazil_vol, ierr)
   end subroutine ice_new_ice_step

   subroutine new_ice_impl(nx, ny, ncat, nilyr, i_lo, i_hi, j_lo, j_hi, dt, rhoi, Lfresh, puny, &
                           hfrazilmin, sref, update_ocn_f, hin_max1, salinz, wet, frzmlt, bank, Tf, &
                           aicen, vicen, Tsfcn, qicen, sicen, aice0, fresh, fsalt, frazil_vol, ierr)
      integer, intent(in) :: nx, ny, ncat, nilyr, i_lo, i_hi, j_lo, j_hi
      real(wp), intent(in) :: dt, rhoi, Lfresh, puny, hfrazilmin, sref, hin_max1
      logical, intent(in) :: update_ocn_f
      real(wp), intent(in) :: salinz(nilyr)
      logical, intent(in) :: wet(nx, ny)
      real(wp), intent(in) :: frzmlt(nx, ny), Tf(nx, ny)
      real(wp), intent(inout) :: bank(nx, ny)
      real(wp), intent(inout) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), Tsfcn(nx, ny, ncat)
      real(wp), intent(inout) :: qicen(nx, ny, ncat, nilyr), sicen(nx, ny, ncat, nilyr)
      real(wp), intent(inout) :: aice0(nx, ny), fresh(nx, ny), fsalt(nx, ny), frazil_vol(nx, ny)
      integer(int32), intent(inout) :: ierr(nx, ny)
      integer :: i, j

      !$acc data present(salinz, wet, frzmlt, bank, Tf, aicen, vicen, Tsfcn, qicen, sicen, aice0, fresh, fsalt, frazil_vol, ierr)
      do concurrent(j=j_lo:j_hi, i=i_lo:i_hi)
         if (wet(i, j) .and. iand(ierr(i, j), ICE_ERR_FATAL_MASK) == 0_int32) then
            call ice_add_new_ice_cell(i, j, nx, ny, ncat, nilyr, dt, frzmlt(i, j), Tf(i, j), hin_max1, &
                                      rhoi, Lfresh, puny, hfrazilmin, sref, update_ocn_f, salinz, &
                                      aicen, vicen, Tsfcn, qicen, sicen, aice0, fresh, fsalt, frazil_vol, ierr)
            if (frzmlt(i, j) > 0.0_wp .and. iand(ierr(i, j), ICE_ERR_FATAL_MASK) == 0_int32) then
               bank(i, j) = 0.0_wp
            end if
         end if
      end do
      !$acc end data
   end subroutine new_ice_impl

end module rdb_ice_new_ice
