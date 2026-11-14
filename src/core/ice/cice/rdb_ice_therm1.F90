!! The vertical thermodynamics step of the CICE-port sea ice: per category the
!! atmosphere boundary layer, the Bitz-Lipscomb temperature solve, the thickness
!! changes and their ocean fluxes, then the merge into the cell fluxes.
module rdb_ice_therm1
   ! SPDX-License-Identifier: BSD-3-Clause
   !
   ! Ported from CICE/Icepack (CICE Consortium, https://github.com/CICE-Consortium).
   ! Upstream: Icepack 9a1f3e4 (Icepack 1.5.3+14), module icepack_therm_vertical
   ! (icepack_step_therm1 and thermo_vertical: the order of the stages and the
   ! per-category flux bookkeeping) and icepack_flux (merge_fluxes).
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
   ! Original authors (from the upstream headers): William H. Lipscomb,
   ! Elizabeth C. Hunke (LANL), C. M. Bitz (UW).
   ! Physics: Bitz and Lipscomb (1999, JGR 104 15669); Maykut and McPhee (1995,
   ! JGR 100 24691); Jordan, Andreas and Makshtas (1999, JGR 104 7785).
   !
   !! **What one step is.**  `ice_step_therm1` runs, in this order and each as
   !! ONE kernel launch over the owned cells: (1) the basal heat flux and the
   !! lateral-melt fraction per cell (`ice_basal_heat_step`: `fbot, rside,
   !! wlat`); (2) the category-parallel column thermodynamics
   !! (`ice_therm1_impl`, `do concurrent (n, j, i)`); (3) the merge of the
   !! per-category fluxes into the cell accumulators with the start-of-step area
   !! `aicen_init` as weight (the existing `ice_flux_merge_impl`); (4) the stress
   !! merge and the OR-reduction of the per-category error words
   !! (`ice_therm1_merge_impl`).  Nothing is gathered to the host, nothing is
   !! launched per category.
   !!
   !! **The per-category column** (`ice_therm1_cell`, in the upstream order).  A
   !! category with `aicen_init <= puny` is skipped: its STATE is not touched
   !! (bitwise), and its per-category outputs are set to zero so the merge adds
   !! nothing for it.  Otherwise, on the start-of-step column: the atmosphere
   !! boundary layer at the start-of-step surface temperature (`ice_atm_boundary_ice`:
   !! the coefficients `shcoef, lhcoef` and the wind stress); the vertical profile
   !! (`ice_init_vertical_profile`: thicknesses, layer temperatures, `einit`);
   !! the temperature solve (`ice_temperature_changes`, which also redistributes
   !! the absorbed shortwave when `sw_redist`); the thickness changes
   !! (`ice_thickness_changes`); the energy check
   !! (`ice_conservation_check_vthermo`); the ocean-flux bookkeeping below; and
   !! the category state update (`ice_update_state_vthermo`).  The column is
   !! gathered ONCE from the `(nx, ny, ncat[, nlyr])` state arrays into fixed
   !! `ICE_*_MAX` locals, every stage works on those locals, and the state is
   !! scattered ONCE at the end.  Only per-category arrays are written, never a
   !! cell `(nx, ny)` accumulator (different categories would race on it).
   !!
   !! **Ocean-flux bookkeeping** (upstream `thermo_vertical`, after the thickness
   !! change).  With `dhi = hin - worki` the TOTAL ice thickness change (growth,
   !! melt and snow-ice) and `dhs = hsn - works - hsn_new` the snow change with
   !! the new snowfall excluded: `freshn = evapn - (rhoi*dhi + rhos*dhs)/dt`,
   !! `fsaltn = -rhoi*dhi*Sref*1e-3/dt`, `fhocnn` as returned by the thickness
   !! change (`fbot` plus the leftover energy), `fswabsn = fswsfcn + fswintn +
   !! fswthrun` (after the shortwave redistribution).  The history accumulators
   !! `melttn, meltbn, congeln, snoicen` [m] are SET per category each step.
   !! Signs are Icepack's (see `rdb_ice_flux_state`): positive into the ocean for
   !! `freshn, fsaltn, fhocnn`; `fsensn, flatn` positive downward; `flwoutn`
   !! negative; `evapn > 0` condensation.
   !!
   !! **The merge.**  `x += sum_n aicen_init(n)*xn(n)` for `fresh, fsalt, fhocn,
   !! evap, fsens, flat, fswabs, meltt, meltb, congel, snoice` and `flwout +=
   !! sum_n aicen_init*(flwoutn - (1 - emis)*flw)` (the existing kernel; the
   !! weight is the PRE-thermodynamics area, so a category that melts out is
   !! still weighted).  The wind stress on the ice is merged the same way and
   !! SET (`strair_x_T = sum_n aicen_init*strairxn`, likewise y): nothing else
   !! produces it, so there is no accumulation to zero first.  Valid on the owned
   !! cells only.  `fswthru`/`fswthru_vis` are the shortwave aggregate's
   !! (`sum aicen*fswthrun`), not merged here.
   !!
   !! **What the caller owes.**  Per cell `(nx, ny)`: the atmosphere (`uatm, vatm`
   !! [m/s] on the model axes, `wind` [m/s] (read without `highfreq`), `potT` [K],
   !! `Qa` [kg/kg], `rhoa` [kg/m3], `flw` [W/m2]), the snowfall `fsnow`
   !! [kg/m2/s] (it lands on the ice categories at the cell value: the area
   !! weighting is the merge's), the T-averaged start-of-step ice velocity `uvel,
   !! vvel` [m/s] (read with `highfreq`), the freezing point `Tf` [degC], and for
   !! the basal heat `frzmlt` [W/m2], `sst` [degC] and the previous step's
   !! ice-ocean stress per ICE area `tau_x, tau_y` at the T point (see
   !! `rdb_ice_basal_heat`).  Scalars: the step `dt`, the heights `zlvl, zlvs` [m].
   !! Per category: the start-of-step area `aicen_init` (`ice_save_init`), the
   !! state, and the shortwave absorbed in the category from the radiation
   !! driver: `fswsfcn, fswintn, fswthrun` [W/m2], `Iswabsn` per ice layer and
   !! `Sswabsn` per snow layer [W/m2] (all ZERO when the radiation driver is not
   !! wired).  They are updated in place only by `sw_redist` (upstream does the
   !! same).  `sw_redist` is a shim argument (the namelist knob does not exist
   !! yet); the iteration cap and thresholds are the BL99 constants.
   !!
   !! **Errors.**  Per-category word `ierr_n(i, j, n)` is WRITTEN (not OR-ed) each
   !! call; the merge kernel ORs the categories of a cell into `ice_err(i, j)`
   !! (so no kernel writes a shared word from the category loop).  Every bit
   !! this step can set is FATAL: a flagged category is NOT committed (state and
   !! shortwave arrays untouched), its per-category flux outputs are zero, and
   !! other cells are unaffected.  The cell word `ice_err` is read-only in the
   !! category launch: a cell that already carries a fatal bit after the basal-heat
   !! launch (a non-finite stress or enthalpy, ...) is SKIPPED whole (no category
   !! commits with `fbot = 0`).  The stages of a category chain on separate
   !! per-stage words (a compiler fault with one word threaded through the calls,
   !! see the cell routine).  A non-finite `aicen_init` is
   !! flagged (it is never silently treated as "no ice"); any non-finite input
   !! of a category, or a non-finite result, flags it with `ICE_ERR_NONFINITE`.
   !! A non-finite cell forcing flags every category of that cell with ice.
   !!
   !! **Not ported.**  Melt ponds, aerosols, isotopes, the first-year area and age
   !! tracers, snow grains, the form-drag atmosphere, `calc_Tsfc = .false.` and
   !! `calc_strair = .false.` (the configuration refuses them), the reference-height
   !! diagnostics `Tref, Qref, Uref`, and the `fsurf, fcondtop, fcondbot` cell
   !! diagnostics.  `sfc_kind, sf_0, dsf_dT` of the cell and flat routines are the
   !! BL99 test hook (`ICE_SFC_LINEAR` fixes the surface balance for closed-form
   !! tests); `ice_step_therm1` always uses the production bulk formula.
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use, intrinsic :: iso_fortran_env, only: int32
   use rdb_constants, only: wp
   use rdb_ice_atmo, only: ice_atm_boundary_ice
   use rdb_ice_basal_heat, only: ice_basal_heat_step
   use rdb_ice_bl99, only: ice_temperature_changes, ICE_NITERMAX, ICE_SW_FRAC, ICE_SW_DTEMP
   use rdb_ice_dims, only: ICE_NCAT_MAX, ICE_NILYR_MAX, ICE_NSLYR_MAX
   use rdb_ice_errors, only: ICE_ERR_NONFINITE, ICE_ERR_FATAL_MASK
   use rdb_ice_flux, only: ice_flux_merge_impl
   use rdb_ice_params, only: ice_constants_t, ice_atm_params_t, ice_params_t
   use rdb_ice_thermo_shared, only: ICE_SFC_BULK
   use rdb_ice_thickness, only: ice_thickness_changes, ice_conservation_check_vthermo, ice_update_state_vthermo, &
                                ICE_FERRMAX
   use rdb_ice_vprofile, only: ice_init_vertical_profile
   implicit none
   private

   public :: ice_step_therm1, ice_therm1_impl, ice_therm1_cell, ice_therm1_merge_impl

contains

   pure subroutine ice_therm1_cell(i, j, n, nx, ny, ncat, nilyr, nslyr, dt, c, atm, &
                                   conduct, l_brine, emissivity, ref_salinity, sw_redist, sw_frac, sw_dtemp, nitermax, &
                                   sfc_kind, sf_0, dsf_dT, zlvl, zlvs, &
                                   uatm, vatm, wind, potT, Qa, rhoa, flw, fsnow, uvel, vvel, Tf, fbot, &
                                   aicen_init, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, &
                                   fswsfcn, fswintn, fswthrun, Iswabsn, Sswabsn, &
                                   freshn, fsaltn, fhocnn, evapn, fsensn, flatn, flwoutn, fswabsn, &
                                   melttn, meltbn, congeln, snoicen, strairxn, strairyn, ice_err, ierr_n)
      !! One `(i, j, n)` category: the whole upstream per-category sequence (see
      !! the module header).  `i, j, n` are passed BY VALUE (a compiler fault
      !! passes a bare non-innermost `do concurrent` index by reference as its
      !! lower bound).  The stages chain on SEPARATE per-stage error words (`ie_a ... ie_u`) merged into `ie`
      !! after each call: with one word threaded through all six `intent(inout)` calls, nvfortran 26.5
      !! (`-stdpar=gpu -fast`) let later stages run after the atmosphere layer had flagged the cell.
      !! Cell inputs are `(nx, ny)`, category arrays `(nx, ny, ncat)`,
      !! layer arrays `(nx, ny, ncat, nilyr|nslyr)`.  Updated in place on success:
      !! `aicen, vicen, vsnon, Tsfcn, qicen, qsnon` and, only with `sw_redist`,
      !! `fswsfcn, fswintn, Iswabsn, Sswabsn`.  Written: the 12 flux outputs
      !! `freshn ... snoicen`, the stress `strairxn, strairyn` [N/m2, on the model
      !! axes, on the ice] and the category error word `ierr_n`.  A skipped or
      !! flagged category writes zeros to the flux outputs and leaves the state alone.
      !$acc routine seq
      integer, value :: i, j, n
      integer, intent(in) :: nx, ny, ncat, nilyr, nslyr, conduct, nitermax, sfc_kind
      real(wp), intent(in) :: dt
      type(ice_constants_t), intent(in) :: c
      type(ice_atm_params_t), intent(in) :: atm
      logical, intent(in) :: l_brine, sw_redist
      real(wp), intent(in) :: emissivity, ref_salinity, sw_frac, sw_dtemp, sf_0, dsf_dT, zlvl, zlvs
      real(wp), intent(in) :: uatm(nx, ny), vatm(nx, ny), wind(nx, ny), potT(nx, ny), Qa(nx, ny), rhoa(nx, ny)
      real(wp), intent(in) :: flw(nx, ny), fsnow(nx, ny), uvel(nx, ny), vvel(nx, ny), Tf(nx, ny), fbot(nx, ny)
      real(wp), intent(in) :: aicen_init(nx, ny, ncat), fswthrun(nx, ny, ncat)
      real(wp), intent(inout) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat), Tsfcn(nx, ny, ncat)
      real(wp), intent(inout) :: qicen(nx, ny, ncat, nilyr), qsnon(nx, ny, ncat, nslyr)
      real(wp), intent(in) :: sicen(nx, ny, ncat, nilyr)
      real(wp), intent(inout) :: fswsfcn(nx, ny, ncat), fswintn(nx, ny, ncat)
      real(wp), intent(inout) :: Iswabsn(nx, ny, ncat, nilyr), Sswabsn(nx, ny, ncat, nslyr)
      real(wp), intent(inout) :: freshn(nx, ny, ncat), fsaltn(nx, ny, ncat), fhocnn(nx, ny, ncat), evapn(nx, ny, ncat)
      real(wp), intent(inout) :: fsensn(nx, ny, ncat), flatn(nx, ny, ncat), flwoutn(nx, ny, ncat), fswabsn(nx, ny, ncat)
      real(wp), intent(inout) :: melttn(nx, ny, ncat), meltbn(nx, ny, ncat), congeln(nx, ny, ncat), snoicen(nx, ny, ncat)
      real(wp), intent(inout) :: strairxn(nx, ny, ncat), strairyn(nx, ny, ncat)
      integer(int32), intent(in) :: ice_err(nx, ny)
      integer(int32), intent(inout) :: ierr_n(nx, ny, ncat)

      integer :: k
      integer(int32) :: ie, ie_a, ie_p, ie_t, ie_h, ie_c, ie_u
      logical :: go, active, commit
      real(wp) :: a_init, a, vi, vs, Tsf, Tbot, hin, hilyr, hsn, hslyr, einit, worki, works, dhi, dhs
      real(wp) :: shcoef, lhcoef, strx, stry, delt, delq, fsw, fsi, fsens_c, flat_c, flwout_c, fsurf_c, fct_c, fcb_c
      real(wp) :: einex, meltt, melts, meltb, congel, snoice, dsnow, hsn_new, fhoc, evap, evaps, evapi, efinal, ferr
      real(wp) :: fresh_c, fsalt_c, fswabs_c, chk
      real(wp) :: zqi(ICE_NILYR_MAX), zSi(ICE_NILYR_MAX), zTi(ICE_NILYR_MAX), Tml(ICE_NILYR_MAX), Isw(ICE_NILYR_MAX)
      real(wp) :: zqs(ICE_NSLYR_MAX), zTs(ICE_NSLYR_MAX), Ssw(ICE_NSLYR_MAX)
      real(wp) :: diag(3)

      ie = 0_int32
      commit = .false.
      a_init = aicen_init(i, j, n)
      if (.not. ieee_is_finite(a_init)) then      ! never read as "no ice"
         ie = ICE_ERR_NONFINITE
      end if
      ! no ice at the start of the step (or an invalid word above): the state is not touched
      active = ie == 0_int32 .and. a_init > c%puny .and. iand(ice_err(i, j), ICE_ERR_FATAL_MASK) == 0_int32

      if (active) then
         ! ---- gather the column once
         a = aicen(i, j, n)
         vi = vicen(i, j, n)
         vs = vsnon(i, j, n)
         Tsf = Tsfcn(i, j, n)
         Tbot = Tf(i, j)
         fsw = fswsfcn(i, j, n)
         fsi = fswintn(i, j, n)
         do k = 1, nilyr
            zqi(k) = qicen(i, j, n, k)
            zSi(k) = sicen(i, j, n, k)
            Isw(k) = Iswabsn(i, j, n, k)
            zTi(k) = 0.0_wp
            Tml(k) = 0.0_wp
         end do
         do k = 1, nslyr
            zqs(k) = qsnon(i, j, n, k)
            Ssw(k) = Sswabsn(i, j, n, k)
            zTs(k) = 0.0_wp
         end do
         diag(1) = 0.0_wp
         diag(2) = 0.0_wp
         diag(3) = 0.0_wp
         meltt = 0.0_wp
         melts = 0.0_wp
         meltb = 0.0_wp
         congel = 0.0_wp
         snoice = 0.0_wp
         dsnow = 0.0_wp
         fsens_c = 0.0_wp
         flat_c = 0.0_wp
         flwout_c = 0.0_wp
         fsurf_c = 0.0_wp
         fct_c = 0.0_wp
         fcb_c = 0.0_wp
         einex = 0.0_wp
         hsn_new = 0.0_wp
         fhoc = 0.0_wp
         evap = 0.0_wp
         evaps = 0.0_wp
         evapi = 0.0_wp
         efinal = 0.0_wp
         ferr = 0.0_wp
         einit = 0.0_wp
         worki = 0.0_wp
         works = 0.0_wp

         ! The forcing is read here, so a non-finite value of the cell flags this category; the
         ! callee checks of the atmosphere and the solver cover the rest.  `fbot, fsnow, Tbot,
         ! fswthrun` are checked explicitly because the later stages read them after comparisons
         ! that would let a NaN pass.
         go = ieee_is_finite(Tbot) .and. ieee_is_finite(fbot(i, j)) .and. ieee_is_finite(fsnow(i, j)) &
              .and. ieee_is_finite(fswthrun(i, j, n)) .and. ieee_is_finite(flw(i, j)) &
              .and. ieee_is_finite(potT(i, j)) .and. ieee_is_finite(Qa(i, j)) .and. ieee_is_finite(rhoa(i, j)) &
              .and. ieee_is_finite(uatm(i, j)) .and. ieee_is_finite(vatm(i, j)) .and. ieee_is_finite(wind(i, j))
         if (atm%highfreq) go = go .and. ieee_is_finite(uvel(i, j)) .and. ieee_is_finite(vvel(i, j))
         if (.not. go) ie = ior(ie, ICE_ERR_NONFINITE)

         ! ---- 1. atmosphere boundary layer at the start-of-step surface temperature
         shcoef = 0.0_wp
         lhcoef = 0.0_wp
         strx = 0.0_wp
         stry = 0.0_wp
         delt = 0.0_wp
         delq = 0.0_wp
         if (ie == 0_int32) then
            ie_a = 0_int32
            call ice_atm_boundary_ice(Tsf, potT(i, j), uatm(i, j), vatm(i, j), uvel(i, j), vvel(i, j), wind(i, j), &
                                      zlvl, zlvs, Qa(i, j), rhoa(i, j), c, atm, shcoef, lhcoef, strx, stry, delt, delq, ie_a)
            ie = ior(ie, ie_a)
         end if

         ! ---- 2. vertical profile
         if (ie == 0_int32) then
            ie_p = 0_int32
            call ice_init_vertical_profile(a, vi, vs, nilyr, nslyr, c, l_brine, zSi, zqi, zqs, zTi, zTs, Tml, &
                                           hin, hilyr, hsn, hslyr, einit, ie_p)
            ie = ior(ie, ie_p)
         end if

         ! ---- 3. temperature solve
         if (ie == 0_int32) then
            ie_t = 0_int32
            worki = hin
            works = hsn
            call ice_temperature_changes(nilyr, nslyr, dt, rhoa(i, j), flw(i, j), potT(i, j), Qa(i, j), shcoef, lhcoef, &
                                         emissivity, sfc_kind, sf_0, dsf_dT, c, conduct, l_brine, &
                                         sw_redist, sw_frac, sw_dtemp, nitermax, &
                                         hilyr, hslyr, zSi, Tbot, einit, &
                                         fsw, fsi, Ssw, Isw, zqi, zTi, zqs, zTs, Tsf, &
                                         fsens_c, flat_c, flwout_c, fsurf_c, fct_c, fcb_c, einex, diag, ie_t)
            ie = ior(ie, ie_t)
         end if

         ! ---- 4. thickness changes and the energy check
         if (ie == 0_int32) then
            ie_h = 0_int32
            call ice_thickness_changes(nilyr, nslyr, c, l_brine, dt, Tbot, fbot(i, j), flat_c, fsurf_c, fct_c, fcb_c, &
                                       fsnow(i, j), zSi(nilyr), hin, hilyr, hsn, hslyr, zqi, zqs, &
                                       meltt, melts, meltb, congel, snoice, dsnow, &
                                       hsn_new, fhoc, evap, evaps, evapi, efinal, ie_h)
            ie = ior(ie, ie_h)
         end if
         if (ie == 0_int32) then
            ie_c = 0_int32
            call ice_conservation_check_vthermo(c, dt, fsurf_c, flat_c, fhoc, fsi, fsnow(i, j), einit, efinal, ICE_FERRMAX, &
                                                ferr, ie_c)
            ie = ior(ie, ie_c)
         end if

         ! ---- 5. ocean-flux bookkeeping and the category state
         fresh_c = 0.0_wp
         fsalt_c = 0.0_wp
         fswabs_c = 0.0_wp
         if (ie == 0_int32) then
            ie_u = 0_int32
            dhi = hin - worki
            dhs = hsn - works - hsn_new
            fresh_c = evap - (c%rhoi*dhi + c%rhos*dhs)/dt
            fsalt_c = -c%rhoi*dhi*ref_salinity*0.001_wp/dt
            fswabs_c = fsw + fsi + fswthrun(i, j, n)
            call ice_update_state_vthermo(nilyr, nslyr, Tbot, hin, hsn, Tsf, zqi, zqs, a, vi, vs, ie_u)
            ie = ior(ie, ie_u)
         end if
         if (ie == 0_int32) then
            chk = ((fresh_c + fsalt_c) + (fhoc + evap)) + ((fsens_c + flat_c) + (flwout_c + fswabs_c)) &
                  + ((meltt + meltb) + (congel + snoice)) + (strx + stry)
            if (.not. ieee_is_finite(chk)) ie = ior(ie, ICE_ERR_NONFINITE)
         end if
         commit = iand(ie, ICE_ERR_FATAL_MASK) == 0_int32
      end if

      ! ---- scatter once; a flagged or skipped category is not committed (its flux outputs are zero)
      if (commit) then
         aicen(i, j, n) = a
         vicen(i, j, n) = vi
         vsnon(i, j, n) = vs
         Tsfcn(i, j, n) = Tsf
         do k = 1, nilyr
            qicen(i, j, n, k) = zqi(k)
         end do
         do k = 1, nslyr
            qsnon(i, j, n, k) = zqs(k)
         end do
         if (sw_redist) then
            fswsfcn(i, j, n) = fsw
            fswintn(i, j, n) = fsi
            do k = 1, nilyr
               Iswabsn(i, j, n, k) = Isw(k)
            end do
            do k = 1, nslyr
               Sswabsn(i, j, n, k) = Ssw(k)
            end do
         end if
         freshn(i, j, n) = fresh_c
         fsaltn(i, j, n) = fsalt_c
         fhocnn(i, j, n) = fhoc
         evapn(i, j, n) = evap
         fsensn(i, j, n) = fsens_c
         flatn(i, j, n) = flat_c
         flwoutn(i, j, n) = flwout_c
         fswabsn(i, j, n) = fswabs_c
         melttn(i, j, n) = meltt
         meltbn(i, j, n) = meltb
         congeln(i, j, n) = congel
         snoicen(i, j, n) = snoice
         strairxn(i, j, n) = strx
         strairyn(i, j, n) = stry
      else
         freshn(i, j, n) = 0.0_wp
         fsaltn(i, j, n) = 0.0_wp
         fhocnn(i, j, n) = 0.0_wp
         evapn(i, j, n) = 0.0_wp
         fsensn(i, j, n) = 0.0_wp
         flatn(i, j, n) = 0.0_wp
         flwoutn(i, j, n) = 0.0_wp
         fswabsn(i, j, n) = 0.0_wp
         melttn(i, j, n) = 0.0_wp
         meltbn(i, j, n) = 0.0_wp
         congeln(i, j, n) = 0.0_wp
         snoicen(i, j, n) = 0.0_wp
         strairxn(i, j, n) = 0.0_wp
         strairyn(i, j, n) = 0.0_wp
      end if
      ierr_n(i, j, n) = ie

   end subroutine ice_therm1_cell

   subroutine ice_therm1_impl(nx, ny, ncat, nilyr, nslyr, is, ie, js, je, dt, c, atm, &
                              conduct, l_brine, emissivity, ref_salinity, sw_redist, sw_frac, sw_dtemp, nitermax, &
                              sfc_kind, sf_0, dsf_dT, zlvl, zlvs, &
                              uatm, vatm, wind, potT, Qa, rhoa, flw, fsnow, uvel, vvel, Tf, fbot, &
                              aicen_init, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, &
                              fswsfcn, fswintn, fswthrun, Iswabsn, Sswabsn, &
                              freshn, fsaltn, fhocnn, evapn, fsensn, flatn, flwoutn, fswabsn, &
                              melttn, meltbn, congeln, snoicen, strairxn, strairyn, ice_err, ierr_n)
      !! The category-parallel launch over the owned cells `is:ie, js:je`: ONE
      !! `do concurrent (n, j, i)` calling `ice_therm1_cell` (arguments as there).
      !! Every array must be device-present on an offload build (`present`
      !! below fails loudly otherwise); `c` and `atm` are the all-scalar bundles
      !! of the (once-mapped) parameter set.
      integer, intent(in) :: nx, ny, ncat, nilyr, nslyr, is, ie, js, je, conduct, nitermax, sfc_kind
      real(wp), intent(in) :: dt
      type(ice_constants_t), intent(in) :: c
      type(ice_atm_params_t), intent(in) :: atm
      logical, intent(in) :: l_brine, sw_redist
      real(wp), intent(in) :: emissivity, ref_salinity, sw_frac, sw_dtemp, sf_0, dsf_dT, zlvl, zlvs
      real(wp), intent(in) :: uatm(nx, ny), vatm(nx, ny), wind(nx, ny), potT(nx, ny), Qa(nx, ny), rhoa(nx, ny)
      real(wp), intent(in) :: flw(nx, ny), fsnow(nx, ny), uvel(nx, ny), vvel(nx, ny), Tf(nx, ny), fbot(nx, ny)
      real(wp), intent(in) :: aicen_init(nx, ny, ncat), fswthrun(nx, ny, ncat)
      real(wp), intent(inout) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat), Tsfcn(nx, ny, ncat)
      real(wp), intent(inout) :: qicen(nx, ny, ncat, nilyr), qsnon(nx, ny, ncat, nslyr)
      real(wp), intent(in) :: sicen(nx, ny, ncat, nilyr)
      real(wp), intent(inout) :: fswsfcn(nx, ny, ncat), fswintn(nx, ny, ncat)
      real(wp), intent(inout) :: Iswabsn(nx, ny, ncat, nilyr), Sswabsn(nx, ny, ncat, nslyr)
      real(wp), intent(inout) :: freshn(nx, ny, ncat), fsaltn(nx, ny, ncat), fhocnn(nx, ny, ncat), evapn(nx, ny, ncat)
      real(wp), intent(inout) :: fsensn(nx, ny, ncat), flatn(nx, ny, ncat), flwoutn(nx, ny, ncat), fswabsn(nx, ny, ncat)
      real(wp), intent(inout) :: melttn(nx, ny, ncat), meltbn(nx, ny, ncat), congeln(nx, ny, ncat), snoicen(nx, ny, ncat)
      real(wp), intent(inout) :: strairxn(nx, ny, ncat), strairyn(nx, ny, ncat)
      integer(int32), intent(in) :: ice_err(nx, ny)
      integer(int32), intent(inout) :: ierr_n(nx, ny, ncat)

      integer :: i, j, n

      !$acc data present(uatm, vatm, wind, potT, Qa, rhoa, flw, fsnow, uvel, vvel, Tf, fbot) &
      !$acc& present(aicen_init, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon) &
      !$acc& present(fswsfcn, fswintn, fswthrun, Iswabsn, Sswabsn) &
      !$acc& present(freshn, fsaltn, fhocnn, evapn, fsensn, flatn, flwoutn, fswabsn) &
      !$acc& present(melttn, meltbn, congeln, snoicen, strairxn, strairyn, ice_err, ierr_n)
      do concurrent(n=1:ncat, j=js:je, i=is:ie)
         call ice_therm1_cell(i, j, n, nx, ny, ncat, nilyr, nslyr, dt, c, atm, &
                              conduct, l_brine, emissivity, ref_salinity, sw_redist, sw_frac, sw_dtemp, nitermax, &
                              sfc_kind, sf_0, dsf_dT, zlvl, zlvs, &
                              uatm, vatm, wind, potT, Qa, rhoa, flw, fsnow, uvel, vvel, Tf, fbot, &
                              aicen_init, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, &
                              fswsfcn, fswintn, fswthrun, Iswabsn, Sswabsn, &
                              freshn, fsaltn, fhocnn, evapn, fsensn, flatn, flwoutn, fswabsn, &
                              melttn, meltbn, congeln, snoicen, strairxn, strairyn, ice_err, ierr_n)
      end do
      !$acc end data
   end subroutine ice_therm1_impl

   pure subroutine ice_therm1_merge_impl(nx, ny, ncat, is, ie, js, je, aicen_init, strairxn, strairyn, ierr_n, &
                                         strair_x_T, strair_y_T, ice_err)
      !! The per-cell part of the merge the flux kernel does not do: the wind stress
      !! on the ice `strair_{x,y}_T = sum_{n=1..ncat} aicen_init*strair{x,y}n`
      !! (SET, fixed category order, owned cells only; N/m2 per cell area) and the
      !! OR of the per-category error words into `ice_err(i, j)`.  A non-finite
      !! stress sum is propagated into the cell field and flagged
      !! `ICE_ERR_NONFINITE`.  One `do concurrent (j, i)` with the category loop inside.
      integer, intent(in) :: nx, ny, ncat, is, ie, js, je
      real(wp), intent(in) :: aicen_init(nx, ny, ncat), strairxn(nx, ny, ncat), strairyn(nx, ny, ncat)
      integer(int32), intent(in) :: ierr_n(nx, ny, ncat)
      real(wp), intent(inout) :: strair_x_T(nx, ny), strair_y_T(nx, ny)
      integer(int32), intent(inout) :: ice_err(nx, ny)

      integer :: i, j, n
      integer(int32) :: w
      real(wp) :: a, sx, sy

      !$acc data present(aicen_init, strairxn, strairyn, ierr_n, strair_x_T, strair_y_T, ice_err)
      do concurrent(j=js:je, i=is:ie) local(n, w, a, sx, sy)
         sx = 0.0_wp
         sy = 0.0_wp
         w = 0_int32
         !$acc loop seq
         do n = 1, ncat
            a = aicen_init(i, j, n)
            sx = sx + a*strairxn(i, j, n)
            sy = sy + a*strairyn(i, j, n)
            w = ior(w, ierr_n(i, j, n))
         end do
         strair_x_T(i, j) = sx
         strair_y_T(i, j) = sy
         if (.not. ieee_is_finite(sx + sy)) w = ior(w, ICE_ERR_NONFINITE)
         ice_err(i, j) = ior(ice_err(i, j), w)
      end do
      !$acc end data
   end subroutine ice_therm1_merge_impl

   subroutine ice_step_therm1(nx, ny, ncat, nilyr, nslyr, is, ie, js, je, dt, p, sw_redist, zlvl, zlvs, &
                              uatm, vatm, wind, potT, Qa, rhoa, flw, fsnow, uvel, vvel, Tf, sst, frzmlt, tau_x, tau_y, &
                              aicen_init, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, &
                              fswsfcn, fswintn, fswthrun, Iswabsn, Sswabsn, &
                              fresh, fsalt, fhocn, evap, fsens, flat, flwout, fswabs, meltt, meltb, congel, snoice, &
                              freshn, fsaltn, fhocnn, evapn, fsensn, flatn, flwoutn, fswabsn, &
                              melttn, meltbn, congeln, snoicen, strairxn, strairyn, strair_x_T, strair_y_T, &
                              fbot, rside, wlat, ierr_n, ice_err, istat)
      !! The host shim of the whole step (see the module header): basal heat
      !! launch, category launch, flux merge launch, stress/error merge launch.
      !! The parameter bundle `p` is unpacked into scalars and its two all-scalar
      !! sub-bundles (`p%c`, `p%atm`) here; on an offload build `p` must have been
      !! mapped ONCE at setup (`enter data copyin(p)`) so the launch does not copy it.
      !! `is:ie, js:je` are the owned cells (`ice_therm1_cell` requires the
      !! ghost-extent arrays `(nx, ny, ...)`).  The cell accumulators `fresh ...
      !! snoice` are ADDED to (zeroed once per step by `ice_flux_zero`), the
      !! per-category outputs, `strair_*_T`, `fbot, rside, wlat` and `ierr_n` are
      !! overwritten, and `ice_err` is OR-ed into.  Cell fluxes are valid on the
      !! owned cells only.  `istat` is 1 for a layer or category count outside the `ICE_*_MAX`
      !! bounds and 2 for a non-finite or non-positive `dt, zlvl, zlvs`; nothing runs then.
      integer, intent(in) :: nx, ny, ncat, nilyr, nslyr, is, ie, js, je
      real(wp), intent(in) :: dt, zlvl, zlvs
      type(ice_params_t), intent(in) :: p
      logical, intent(in) :: sw_redist
      real(wp), intent(in) :: uatm(nx, ny), vatm(nx, ny), wind(nx, ny), potT(nx, ny), Qa(nx, ny), rhoa(nx, ny)
      real(wp), intent(in) :: flw(nx, ny), fsnow(nx, ny), uvel(nx, ny), vvel(nx, ny), Tf(nx, ny), sst(nx, ny)
      real(wp), intent(in) :: frzmlt(nx, ny), tau_x(nx, ny), tau_y(nx, ny)
      real(wp), intent(in) :: aicen_init(nx, ny, ncat), fswthrun(nx, ny, ncat)
      real(wp), intent(inout) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat), Tsfcn(nx, ny, ncat)
      real(wp), intent(inout) :: qicen(nx, ny, ncat, nilyr), qsnon(nx, ny, ncat, nslyr)
      real(wp), intent(in) :: sicen(nx, ny, ncat, nilyr)
      real(wp), intent(inout) :: fswsfcn(nx, ny, ncat), fswintn(nx, ny, ncat)
      real(wp), intent(inout) :: Iswabsn(nx, ny, ncat, nilyr), Sswabsn(nx, ny, ncat, nslyr)
      real(wp), intent(inout) :: fresh(nx, ny), fsalt(nx, ny), fhocn(nx, ny), evap(nx, ny), fsens(nx, ny), flat(nx, ny)
      real(wp), intent(inout) :: flwout(nx, ny), fswabs(nx, ny), meltt(nx, ny), meltb(nx, ny), congel(nx, ny), snoice(nx, ny)
      real(wp), intent(inout) :: freshn(nx, ny, ncat), fsaltn(nx, ny, ncat), fhocnn(nx, ny, ncat), evapn(nx, ny, ncat)
      real(wp), intent(inout) :: fsensn(nx, ny, ncat), flatn(nx, ny, ncat), flwoutn(nx, ny, ncat), fswabsn(nx, ny, ncat)
      real(wp), intent(inout) :: melttn(nx, ny, ncat), meltbn(nx, ny, ncat), congeln(nx, ny, ncat), snoicen(nx, ny, ncat)
      real(wp), intent(inout) :: strairxn(nx, ny, ncat), strairyn(nx, ny, ncat), strair_x_T(nx, ny), strair_y_T(nx, ny)
      real(wp), intent(inout) :: fbot(nx, ny), rside(nx, ny), wlat(nx, ny)
      integer(int32), intent(inout) :: ierr_n(nx, ny, ncat), ice_err(nx, ny)
      integer, intent(out) :: istat

      ! configuration validated ONCE here, never per category: nothing runs (and `istat /= 0`) on a bad one
      istat = 0
      if (ncat < 1 .or. ncat > ICE_NCAT_MAX .or. nilyr < 1 .or. nilyr > ICE_NILYR_MAX &
          .or. nslyr < 1 .or. nslyr > ICE_NSLYR_MAX) istat = 1
      if (.not. (ieee_is_finite(dt) .and. ieee_is_finite(zlvl) .and. ieee_is_finite(zlvs))) then
         istat = 2
      else if (.not. (dt > 0.0_wp .and. zlvl > 0.0_wp .and. zlvs > 0.0_wp)) then
         istat = 2
      end if
      if (istat /= 0) return
      call ice_basal_heat_step(nx, ny, ncat, nilyr, nslyr, is, ie, js, je, dt, p, &
                               aicen, vicen, vsnon, qicen, qsnon, frzmlt, sst, Tf, tau_x, tau_y, &
                               fbot, rside, wlat, ice_err)
      call ice_therm1_impl(nx, ny, ncat, nilyr, nslyr, is, ie, js, je, dt, p%c, p%atm, &
                           p%thermo%conduct, p%thermo%l_brine, p%thermo%emissivity, p%thermo%ice_ref_salinity, &
                           sw_redist, ICE_SW_FRAC, ICE_SW_DTEMP, ICE_NITERMAX, ICE_SFC_BULK, 0.0_wp, 0.0_wp, zlvl, zlvs, &
                           uatm, vatm, wind, potT, Qa, rhoa, flw, fsnow, uvel, vvel, Tf, fbot, &
                           aicen_init, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, &
                           fswsfcn, fswintn, fswthrun, Iswabsn, Sswabsn, &
                           freshn, fsaltn, fhocnn, evapn, fsensn, flatn, flwoutn, fswabsn, &
                           melttn, meltbn, congeln, snoicen, strairxn, strairyn, ice_err, ierr_n)
      call ice_flux_merge_impl(nx, ny, ncat, is, ie, js, je, p%thermo%emissivity, p%c%puny, flw, aicen_init, &
                               freshn, fsaltn, fhocnn, evapn, fsensn, flatn, flwoutn, fswabsn, &
                               melttn, meltbn, congeln, snoicen, &
                               fresh, fsalt, fhocn, evap, fsens, flat, flwout, fswabs, &
                               meltt, meltb, congel, snoice, ice_err)
      call ice_therm1_merge_impl(nx, ny, ncat, is, ie, js, je, aicen_init, strairxn, strairyn, ierr_n, &
                                 strair_x_T, strair_y_T, ice_err)
   end subroutine ice_step_therm1

end module rdb_ice_therm1
