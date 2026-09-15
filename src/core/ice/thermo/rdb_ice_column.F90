!! Single-column Winton NkIce-layer ice thermodynamics (SIS2 port, PR 3a).
module rdb_ice_column
   !! Port of `ice_temp_SIS2` + `laytemp_SIS2` + `update_lay_enth`
   !! (SIS2_ice_thm.F90:169-945, Apache-2.0) under the `ICE_CP_BRINE ==
   !! ICE_CP_ICE` simplification (`rdb_ice_enthalpy` module docstring):
   !! every per-layer implicit solve and T<->E inversion is a
   !! closed-form quadratic — no Newton / false-position iteration
   !! anywhere. Ported line-by-line from the validated stdlib-only
   !! Python prototype `tmp_local_artifacts/ice_pr3a_prototype/
   !! sis2_column.py`; when any formula here and SIS2 itself seem to
   !! disagree, the prototype (which reproduces SIS2's own commented-out
   !! `col_check` energy-closure diagnostic to ~1e-14 fractional) is
   !! the tiebreaker.
   !!
   !! **Index convention — TRAP #2.** Internal columns are TOP-DOWN,
   !! identical to SIS2 and the prototype: index `0` = snow, `1..nk` =
   !! ice, top to bottom. Roundabout's state (`ocean_sea_ice_t`) is
   !! BOTTOM-UP: `enth_ice(..., 1)` = ice bottom (ocean side),
   !! `enth_ice(..., nk_ice)` = ice top (atm/snow side). The flip
   !! happens ONLY at the gather/scatter boundary in `ice_column_step`
   !! (`local(k) = state(nk+1-k)`, the `kg = nz+1-k` idiom of
   !! `kappa_shear_column_driver`, rdb_ocean_kappa_shear.F90:366) —
   !! nothing inside `ice_temp_sis2` / `laytemp_sis2` / `update_lay_enth`
   !! ever sees the bottom-up convention.
   !!
   !! **Ocean-freeze enthalpy — TRAP #1.** The bottom-freeze ocean-side
   !! enthalpy is the LIQUID formula `ice_enthalpy_liquid(sst, s_surf)`
   !! (SIS_slow_thermo.F90:981), never the frozen/mushy
   !! `ice_enth_from_ts(tfw, sice)` — the latter is ~8-9x more negative
   !! at typical sea-ice salinities, inflating freeze mass per Joule and
   !! making Stefan growth ~3x too fast. Lands in `ice_column_step`
   !! step 4 (§3.4).
   !!
   !! **`fb` is a post-hoc residual, not a matrix BC — TRAP #3.** The
   !! conduction matrix's bottom row always couples to the freezing
   !! temperature `tfw` (`cc(nk+1) = 2*kk*dtt`); `fb` (ocean->ice heat
   !! flux) enters exactly once, AFTER the conservative enthalpy update,
   !! as `bmelt = bmelt + (dtt*fb - tflux_bot)`. See `ice_temp_sis2`.
   !!
   !! **`bb(k)` two-branch live formula — TRAP #4.** Not dead code: see
   !! `ice_temp_sis2`'s `bb` computation, the `(Cp_brine - Cp_ice)` term
   !! kept explicit even though it is 0 under the simplification
   !! (matches the prototype's chosen style, sis2_column.py:214).
   !!
   !! **`nk_ice == 2` quasi-conservative double-pass — TRAP #5.** After
   !! the up/down tridiagonal estimate, every layer temperature is
   !! re-solved via `laytemp_sis2` (SIS2_ice_thm.F90:372-388) — this is
   !! what pulls the Stefan-problem error under 1%; see `ice_temp_sis2`.
   !! `ICE_CP_BRINE == ICE_CP_ICE` is asserted at slot init
   !! (`ocean_sea_ice_init`, `rdb_ice_state`) so the unported
   !! Newton/false-position branches (SIS2_ice_thm.F90:625-692,
   !! 837-870, 1876-1933) are provably unreachable.
   !!
   !! **Snow lands AFTER optics + conduction — TRAP #6 (PR 26).**
   !! `ice_column_step`'s `snow` argument is added in step 4 (resize,
   !! via `ice_snow_accumulate`), which runs strictly after step 2
   !! (optics, `ice_optics_csim4`) and step 3 (conduction,
   !! `ice_temp_sis2`) already used the PRE-snowfall `m_snow`. This is
   !! SIS2's fast/slow split (`ice_resize_SIS2` runs after the
   !! conduction solve, SIS2_ice_thm.F90:1122), NOT an oversight: new
   !! snow IS meltable in the same window (`ice_top_melt_peel` starts
   !! at k=0), but its albedo and conduction effect are felt only on
   !! the NEXT window. Do not move the add earlier chasing "why doesn't
   !! the albedo respond immediately".
   !!
   !! **Flood after melt, before rebalance — TRAP #7 (PR 27).** The
   !! Archimedes freeboard snow-ice flood (`ice_snow_ice_flood`) runs
   !! AFTER `ice_bottom_melt_peel` (so freshly-converted mass is not
   !! re-melted this step and `m_i` reflects the post-melt column) and
   !! BEFORE `ice_rebalance_layers` (else the new mass in layer 1 is
   !! never redistributed across the `nk` layers) — exactly SIS2's order
   !! (`ice_resize_SIS2` -> `rebalance_ice_layers`,
   !! `SIS_slow_thermo.F90:998,1008`). It writes local index 1 (the
   !! TOP ice layer, TRAP #2) — do not "helpfully" index `nk`.
   !!
   !! Everything pure; every per-column routine `!$acc routine seq`;
   !! explicit-shape dummies with integer dims declared before the
   !! arrays that use them (decl-order); fixed-size locals capped by
   !! `ICE_NK_MAX` (kappa-shear `NZ_STACK_MAX` precedent) in the driver.
   use rdb_constants, only: wp, H_VANISHED
   use rdb_ice_enthalpy, only: ICE_LAT_FUS, ICE_CP_ICE, ICE_CP_WATER, ICE_CP_BRINE, &
                               ICE_DTF_DS, ICE_ENTH_LIQ_0, ICE_LIQ_LIM, ICE_NK_MAX, &
                               ice_t_freeze, ice_enth_from_ts, ice_temp_from_en_s, &
                               ice_enthalpy_liquid_freeze, ice_enthalpy_liquid
   use rdb_ice_optics, only: ice_optics_csim4
   use rdb_ice_mass, only: ice_snow_accumulate, ice_bottom_freeze, ice_top_melt_peel, &
                           ice_bottom_melt_peel, ice_rebalance_layers, ice_snow_ice_flood
   implicit none
   private

   public :: laytemp_sis2, update_lay_enth, ice_temp_sis2, ice_column_step, &
             ice_thermo_columns
   public :: ICE_K_ICE, ICE_K_SNOW, ICE_RHO_ICE, ICE_RHO_SNOW, ICE_RHO_OCEAN, &
             ICE_H_LO_LIM, ICE_TEMP_RANGE_EST, ICE_BULK_SALINITY, ICE_NK_MAX

   ! ---- Constants (get_param defaults, SIS2_ice_thm.F90:1534-1617 unless noted) ----
   real(wp), parameter :: ICE_K_ICE = 2.03_wp
      !! Bulk ice thermal conductivity (W/m/K) — SIS2 `ICE_CONDUCTIVITY`.
   real(wp), parameter :: ICE_K_SNOW = 0.31_wp
      !! Bulk snow thermal conductivity (W/m/K) — SIS2 `SNOW_CONDUCTIVITY`.
   real(wp), parameter :: ICE_RHO_ICE = 905.0_wp
      !! Nominal sea-ice density (kg/m³).
   real(wp), parameter :: ICE_RHO_SNOW = 330.0_wp
      !! Nominal snow density (kg/m³).
   real(wp), parameter :: ICE_RHO_OCEAN = 1030.0_wp
      !! Nominal seawater reference density (kg/m³) — SIS2 `RHO_OCEAN`.
      !! Consumed by the PR-27 Archimedes freeboard flood
      !! (`ice_column_step` passes `ICE_RHO_ICE/ICE_RHO_OCEAN` into
      !! `ice_snow_ice_flood`). Deliberately independent of
      !! `&ocean_ice_nml rho_ocean` (the EVP ice-drag reference density,
      !! `rdb_ice_evp`) — do not unify them, that would silently couple
      !! the flood threshold to an EVP tuning knob.
   real(wp), parameter :: ICE_H_LO_LIM = 0.0_wp
      !! `MIN_H_FOR_TEMP_CALC` (m) — floor applied in the effective
      !! layer-thickness expressions of `ice_temp_sis2`. Kept in the
      !! algebra at 0 per the prototype (sis2_thermo.py:41).
   real(wp), parameter :: ICE_TEMP_RANGE_EST = 40.0_wp
      !! `temp_range_est` default (K) — feeds `heat_flux_err_rat`
      !! (SIS2_ice_thm.F90:397ff; prototype sis2_column.py:291).
   real(wp), parameter :: ICE_BULK_SALINITY = 4.0_wp
      !! `ICE_BULK_SALINITY` (SIS_slow_thermo.F90:1604) — new-ice
      !! salinity used by `ice_bottom_freeze`'s `salin_freeze`, AND the
      !! `sal_ice` state-init value (`rdb_ice_state`). NOT 5 — resolves
      !! to the prototype/SIS2 default 4.0 (run_validation.py:21
      !! `sice_val=4.0`).
   ! ICE_NK_MAX moved to the `rdb_ice_enthalpy` leaf (imported above and
   ! re-exported via the `public` list) so `rdb_ice_mass` can size its
   ! fixed device-stack locals by it without a module cycle
   ! (`rdb_ice_column` USES `rdb_ice_mass`). See the docstring there.

contains

   pure function laytemp_sis2(m, t_fr, qf, bf, tp, dtt) result(new_temp)
      !! Per-layer implicit heat-budget solve for the new layer
      !! temperature — SIS2 `laytemp_SIS2` (SIS2_ice_thm.F90:544-700),
      !! `ICE_CP_BRINE == ICE_CP_ICE` closed-form branches only (the
      !! Newton/false-position refinement at :625-692 is dead code
      !! under the simplification and is deliberately NOT ported).
      !! Port of prototype `sis2_column.py:21-55`.
      !$acc routine seq
      real(wp), intent(in) :: m
         !! Layer mass (kg/m²).
      real(wp), intent(in) :: t_fr
         !! Layer freezing temperature (degC); 0 for snow/fresh water.
      real(wp), intent(in) :: qf
         !! Forcing heat flux into the layer (W/m²).
      real(wp), intent(in) :: bf
         !! Implicit coupling coefficient to the neighbour temperature
         !! (W/m²/K).
      real(wp), intent(in) :: tp
         !! Previous-step layer temperature (degC).
      real(wp), intent(in) :: dtt
         !! Timestep (s).
      real(wp) :: new_temp

      real(wp) :: e0, aa, bb, cc, disc

      if (t_fr == 0.0_wp) then
         ! Fresh water / snow linear branch (SIS2:585-592).
         new_temp = (m*ICE_CP_ICE*tp + qf*dtt)/(m*ICE_CP_ICE + bf*dtt)
      else
         if (tp >= t_fr) then
            e0 = ICE_CP_WATER*(tp - t_fr)
         else
            ! (Cp_brine - Cp_ice) term vanishes under the simplification.
            e0 = ICE_CP_ICE*(tp - t_fr) - ICE_LAT_FUS*(1.0_wp - t_fr/tp)
         end if

         if (m*e0 + dtt*(qf - bf*t_fr) >= 0.0_wp) then
            ! Layer would be fully melted -> pin to freezing (SIS2:606).
            new_temp = t_fr
         else
            aa = m*ICE_CP_ICE + bf*dtt
            bb = -(m*((e0 + ICE_LAT_FUS) + ICE_CP_ICE*t_fr) + qf*dtt)
            cc = m*ICE_LAT_FUS*t_fr
            disc = max(bb*bb - 4.0_wp*aa*cc, 0.0_wp)
            if (bb >= 0.0_wp) then
               new_temp = -(bb + sqrt(disc))/(2.0_wp*aa)
            else
               new_temp = (2.0_wp*cc)/(-bb + sqrt(disc))
            end if
            ! Cp_ice == Cp_brine -> the quadratic root IS the final
            ! answer; the Newton/false-position loop is not ported.
         end if
      end if

      new_temp = min(new_temp, t_fr)
   end function laytemp_sis2

   pure subroutine update_lay_enth(m_lay, sice, enth, ftop, ht_body, fbot, &
                                   dftop_dt, dfbot_dt, dtt, hf_err_rat, &
                                   extra_heat, new_temp, has_temp_max, temp_max)
      !! Conservative per-layer implicit enthalpy update — SIS2
      !! `update_lay_enth` (SIS2_ice_thm.F90:704-945), closed-form
      !! branches only. Port of prototype `sis2_column.py:58-135`.
      !! Four solution branches (massless layer; pin-to-max with
      !! banked `extra_enth`; fresh `sice==0` linear; salty quadratic),
      !! then the three-way explicit-vs-conservation-inverted flux
      !! bookkeeping (prototype :117-133, incl. the `denom > 0` guard).
      !!
      !! `temp_max` is optional in the SIS2 signature; here it is a
      !! `has_temp_max` logical + `temp_max` value pair (device-routine
      !! `optional` dummies are avoided — same-module call sites only).
      !$acc routine seq
      real(wp), intent(in) :: m_lay
         !! Layer mass (kg/m²).
      real(wp), intent(in) :: sice
         !! Layer bulk salinity (PSU); 0 for snow/fresh.
      real(wp), intent(inout) :: enth
         !! Layer specific enthalpy (J/kg); in = prior step, out = new.
      real(wp), intent(inout) :: ftop
         !! Heat flux at the layer's top interface (W/m²); in = prior
         !! estimate, out = updated (explicit or conservation-inverted).
      real(wp), intent(in) :: ht_body
         !! In-layer heating (solar absorption) (W/m²).
      real(wp), intent(inout) :: fbot
         !! Heat flux at the layer's bottom interface (W/m²); in/out as
         !! `ftop`.
      real(wp), intent(in) :: dftop_dt
         !! d(ftop)/d(new_temp) (W/m²/K).
      real(wp), intent(in) :: dfbot_dt
         !! d(fbot)/d(new_temp) (W/m²/K).
      real(wp), intent(in) :: dtt
         !! Timestep (s).
      real(wp), intent(in) :: hf_err_rat
         !! Precomputed `heat_flux_err_rat` (degC*s/J) deciding explicit
         !! vs conservation-inverted flux bookkeeping.
      real(wp), intent(out) :: extra_heat
         !! Banked excess heat when pinned to `temp_max` (J/m²).
      real(wp), intent(out) :: new_temp
         !! Resulting layer temperature (degC).
      logical, intent(in) :: has_temp_max
         !! True when an explicit `temp_max` clamp applies (snow-branch
         !! call site); false uses the freezing point as the max.
      real(wp), intent(in) :: temp_max
         !! Explicit temperature ceiling (degC), used only when
         !! `has_temp_max`.

      real(wp) :: ftop_in, fbot_in, htg, fb, t_fr, enth_fp
      real(wp) :: max_temp, max_enth, enth_in, extra_enth
      real(wp) :: en_j, aa, bb, cc, disc, dt_denth
      real(wp) :: denom, dflux_dtot_dt

      ftop_in = ftop
      fbot_in = fbot
      htg = (ht_body + ftop_in) - fbot_in
      fb = -(dftop_dt - dfbot_dt)

      extra_heat = 0.0_wp
      extra_enth = 0.0_wp
      if (sice > 0.0_wp) then
         t_fr = ice_t_freeze(sice)
         enth_fp = ice_enthalpy_liquid_freeze(sice)
      else
         t_fr = 0.0_wp
         enth_fp = ice_enth_from_ts(0.0_wp, 0.0_wp)
      end if

      max_temp = t_fr
      max_enth = enth_fp
      if (has_temp_max) then
         if (temp_max < t_fr) then
            max_temp = temp_max
            max_enth = ice_enth_from_ts(temp_max, sice)
         end if
      end if

      enth_in = enth

      if (m_lay == 0.0_wp) then
         new_temp = min(htg/fb, max_temp)
         enth = ice_enth_from_ts(new_temp, sice)
      else if (dtt*(htg - fb*max_temp) >= m_lay*(max_enth - enth_in)) then
         ! Heat applied would push the layer above max_temp -> pin and
         ! bank the excess heat.
         extra_enth = m_lay*(enth_in - max_enth) + dtt*(htg - fb*max_temp)
         extra_heat = extra_enth
         new_temp = max_temp
         enth = max_enth
      else if (sice == 0.0_wp) then
         dt_denth = 1.0_wp/ICE_CP_ICE
         enth = enth_fp + (dtt*htg + m_lay*(enth_in - enth_fp))/ &
                (m_lay + dtt*(fb*dt_denth))
         new_temp = dt_denth*((dtt*htg + m_lay*(enth_in - enth_fp))/ &
                              (m_lay + dtt*(fb*dt_denth)))
      else
         en_j = enth_in - ice_enthalpy_liquid(0.0_wp, 0.0_wp)
         aa = m_lay*ICE_CP_ICE + fb*dtt
         bb = -(m_lay*((en_j - (ICE_CP_WATER - ICE_CP_ICE)*t_fr) + ICE_LAT_FUS) + htg*dtt)
         cc = m_lay*ICE_LAT_FUS*t_fr
         disc = max(bb*bb - 4.0_wp*aa*cc, 0.0_wp)
         if (bb >= 0.0_wp) then
            new_temp = -(bb + sqrt(disc))/(2.0_wp*aa)
         else
            new_temp = (2.0_wp*cc)/(-bb + sqrt(disc))
         end if
         ! Cp_ice == Cp_brine -> "keep this solution" (SIS2:837).
         enth = ice_enth_from_ts(new_temp, sice)
      end if

      ! Decide explicit vs. conservation-inverted flux bookkeeping
      ! (SIS2:913-941; prototype :117-133).
      if (abs(hf_err_rat*dftop_dt) <= m_lay) then
         ftop = ftop_in + dftop_dt*new_temp
         if (hf_err_rat*dfbot_dt <= m_lay) then
            fbot = fbot_in + dfbot_dt*new_temp
         else
            fbot = (ht_body + ftop) - (m_lay*(enth - enth_in) + extra_enth)/dtt
         end if
      else if (hf_err_rat*dfbot_dt <= m_lay) then
         fbot = fbot_in + dfbot_dt*new_temp
         ftop = (fbot - ht_body) + (m_lay*(enth - enth_in) + extra_enth)/dtt
      else
         denom = dfbot_dt - dftop_dt
         if (denom > 0.0_wp) then
            dflux_dtot_dt = (htg - (m_lay*(enth - enth_in) + extra_enth)/dtt)/denom
         else
            dflux_dtot_dt = 0.0_wp
         end if
         ftop = ftop_in + dftop_dt*dflux_dtot_dt
         fbot = fbot_in + dfbot_dt*dflux_dtot_dt
      end if
   end subroutine update_lay_enth

   pure subroutine ice_temp_sis2(nk, m_snow, m_ice_tot, sice, enthalpy, &
                                 sf_0, dsf_dt, sol, tfw, fb, dtt, &
                                 tsurf, tmelt, bmelt, &
                                 col_enth_in, col_enth_out, sum_sol, &
                                 tflux_sfc, tflux_bot)
      !! SEB + vertical-conduction column solve — SIS2 `ice_temp_SIS2`
      !! (SIS2_ice_thm.F90:169-540). Port of prototype
      !! `sis2_column.py:151-412`. TOP-DOWN column (index 0 = snow,
      !! 1..nk = ice top->bottom) — see module docstring TRAP #2.
      !!
      !! The five diag outputs (`col_enth_in`, `col_enth_out`,
      !! `sum_sol`, `tflux_sfc`, `tflux_bot`) are ALWAYS computed (a
      !! handful of flops on a small column) and feed the
      !! `column_energy_closure` test's identity: `col_enth_out -
      !! col_enth_in == sum_sol + tflux_sfc + tflux_bot` (`tflux_bot`
      !! ADDED — already the signed contribution, prototype :198-227).
      !! `col_enth_out` is measured AFTER the conservative update but
      !! BEFORE the liq-lim clamp (prototype `col_enth2b`) — the clamp
      !! moves energy into tmelt/bmelt, outside this identity.
      !$acc routine seq
      integer, intent(in) :: nk
         !! Number of ice layers (declared first — decl-order).
      real(wp), intent(in) :: m_snow
         !! Snow mass per unit area (kg/m²).
      real(wp), intent(in) :: m_ice_tot
         !! Total ice mass per unit area (kg/m²).
      real(wp), intent(in) :: sice(nk)
         !! TOP-DOWN ice bulk salinities (PSU).
      real(wp), intent(inout) :: enthalpy(0:nk)
         !! TOP-DOWN specific enthalpies (J/kg): 0 = snow, 1..nk = ice.
      real(wp), intent(in) :: sf_0
         !! Linearized SEB intercept (W/m²), upward-positive: `SF(T) =
         !! sf_0 + dsf_dt*T`.
      real(wp), intent(in) :: dsf_dt
         !! Linearized SEB slope (W/m²/K), upward-positive.
      real(wp), intent(in) :: sol(0:nk)
         !! Absorbed solar per layer (W/m²), TOP-DOWN.
      real(wp), intent(in) :: tfw
         !! Seawater freezing temperature at the ice base (degC).
      real(wp), intent(in) :: fb
         !! Ocean -> ice-base heat flux (W/m²); post-hoc residual only
         !! — TRAP #3.
      real(wp), intent(in) :: dtt
         !! Timestep (s).
      real(wp), intent(out) :: tsurf
         !! Surface skin temperature (degC).
      real(wp), intent(inout) :: tmelt
         !! Accumulated top melting energy (J/m²); caller zeroes per step.
      real(wp), intent(inout) :: bmelt
         !! Accumulated bottom melting/freezing energy (J/m²); caller
         !! zeroes per step.
      real(wp), intent(out) :: col_enth_in
         !! Column enthalpy Σ m_lay*enth BEFORE anything (diag).
      real(wp), intent(out) :: col_enth_out
         !! Column enthalpy Σ m_lay*enth AFTER the conservative update,
         !! BEFORE the liq-lim clamp (diag).
      real(wp), intent(out) :: sum_sol
         !! Σ sol*dtt over the column (diag, J/m²).
      real(wp), intent(out) :: tflux_sfc
         !! Time-integrated surface heat flux into the column (diag,
         !! J/m²).
      real(wp), intent(out) :: tflux_bot
         !! Time-integrated basal heat flux into the column (diag,
         !! J/m²).

      real(wp) :: temp_ic(0:ICE_NK_MAX), tfi(ICE_NK_MAX)
      real(wp) :: ml_ice, ml_snow, hl_ice_eff, hsnow_eff, tsf
      real(wp) :: kk, k10, k0a, k0skin, k0a_x_ta
      real(wp) :: m_lay(0:ICE_NK_MAX)
      real(wp) :: bb(0:ICE_NK_MAX), cc(0:ICE_NK_MAX + 1), cc_bb(0:ICE_NK_MAX)
      real(wp) :: temp_est(0:ICE_NK_MAX)
      real(wp) :: heat_flux_int(-1:ICE_NK_MAX)
      real(wp) :: b_denom_1, i_bb, comp_rat, tsurf_est, m_pond
      real(wp) :: heat_flux_err_rat
      real(wp) :: e_extra, e_extra_sum, ftop_new, fbot_new, snow_temp_max, snow_temp_new
      real(wp) :: enth_liq_lim, i_liq_lim
      integer :: k

      ! ---- T<->E inversion of the incoming state (top-down) ----
      temp_ic(0) = ice_temp_from_en_s(enthalpy(0), 0.0_wp)
      do k = 1, nk
         temp_ic(k) = ice_temp_from_en_s(enthalpy(k), sice(k))
      end do

      ml_ice = m_ice_tot/real(nk, wp)
      ml_snow = m_snow
      do k = 1, nk
         tfi(k) = ice_t_freeze(sice(k))
      end do

      hl_ice_eff = max(ml_ice/ICE_RHO_ICE, ICE_H_LO_LIM)
      hsnow_eff = ml_snow/ICE_RHO_SNOW + max(1.0e-35_wp, 1.0e-20_wp*ICE_H_LO_LIM)

      tsf = tfi(1)
      if (ml_snow > 0.0_wp) tsf = 0.0_wp

      kk = ICE_K_ICE/hl_ice_eff
      k10 = 2.0_wp*(ICE_K_SNOW*ICE_K_ICE)/(hl_ice_eff*ICE_K_SNOW + hsnow_eff*ICE_K_ICE)
      k0a = (ICE_K_SNOW*dsf_dt)/(0.5_wp*dsf_dt*hsnow_eff + ICE_K_SNOW)
      k0skin = 2.0_wp*ICE_K_SNOW/hsnow_eff
      k0a_x_ta = (ICE_K_SNOW*sf_0)/(0.5_wp*dsf_dt*hsnow_eff + ICE_K_SNOW)

      m_lay(0) = ml_snow
      do k = 1, nk
         m_lay(k) = ml_ice
      end do

      col_enth_in = 0.0_wp
      do k = 0, nk
         col_enth_in = col_enth_in + m_lay(k)*enthalpy(k)
      end do

      ! ---- Effective layer heat capacities bb(k) — TRAP #4 ----
      bb(0) = ml_snow*ICE_CP_ICE
      do k = 1, nk
         if (tfi(k) >= 0.0_wp) then
            bb(k) = ml_ice*ICE_CP_ICE
         else if (temp_ic(k) < tfi(k)) then
            bb(k) = ml_ice*(ICE_CP_ICE - (tfi(k)/temp_ic(k)**2)* &
                            (ICE_LAT_FUS - (ICE_CP_BRINE - ICE_CP_ICE)*temp_ic(k)))
         else
            bb(k) = ml_ice*(ICE_CP_BRINE - ICE_LAT_FUS/tfi(k))
         end if
      end do

      ! ---- Coupling coefficients cc — TRAP #3 (bottom couples to tfw) ----
      cc(0) = k0a*dtt
      cc(1) = k10*dtt
      do k = 2, nk
         cc(k) = kk*dtt
      end do
      cc(nk + 1) = 2.0_wp*kk*dtt

      ! ---- UP sweep ----
      b_denom_1 = bb(nk) + cc(nk + 1)
      i_bb = 1.0_wp/(b_denom_1 + cc(nk))
      temp_est(nk) = ((sol(nk)*dtt + bb(nk)*temp_ic(nk)) + cc(nk + 1)*tfw)*i_bb
      comp_rat = b_denom_1*i_bb
      cc_bb(nk) = cc(nk)*i_bb

      do k = nk - 1, 1, -1
         b_denom_1 = bb(k) + comp_rat*cc(k + 1)
         i_bb = 1.0_wp/(b_denom_1 + cc(k))
         temp_est(k) = ((sol(k)*dtt + bb(k)*temp_ic(k)) + cc(k + 1)*temp_est(k + 1))*i_bb
         comp_rat = b_denom_1*i_bb
         cc_bb(k) = cc(k)*i_bb
      end do

      b_denom_1 = bb(0) + comp_rat*cc(1)
      i_bb = 1.0_wp/(b_denom_1 + cc(0))
      temp_est(0) = (((sol(0)*dtt + bb(0)*temp_ic(0)) - k0a_x_ta*dtt) + cc(1)*temp_est(1))*i_bb

      tsurf_est = (k0skin*temp_est(0) - sf_0)/(dsf_dt + k0skin)

      m_pond = 0.0_wp
      if (tsurf_est > tsf .or. m_pond > 0.0_wp) then
         tsurf_est = tsf
         i_bb = 1.0_wp/(b_denom_1 + k0skin*dtt)
         temp_est(0) = min(tsf, &
                           (((sol(0)*dtt + bb(0)*temp_ic(0)) + k0skin*dtt*tsf) + &
                            cc(1)*temp_est(1))*i_bb)
      end if

      ! ---- DOWN sweep ----
      do k = 1, nk
         temp_est(k) = min(temp_est(k) + cc_bb(k)*temp_est(k - 1), tfi(k))
      end do

      ! ---- Quasi-conservative re-solve via laytemp_sis2 — TRAP #5 ----
      if (nk == 1) then
         temp_est(1) = laytemp_sis2(ml_ice, tfi(1), &
                                    sol(1) + (2.0_wp*kk*tfw + k10*temp_est(0)), &
                                    2.0_wp*kk + k10, temp_ic(1), dtt)
      else
         temp_est(nk) = laytemp_sis2(ml_ice, tfi(nk), &
                                     sol(nk) + kk*(2.0_wp*tfw + temp_est(nk - 1)), &
                                     3.0_wp*kk, temp_ic(nk), dtt)
         do k = nk - 1, 2, -1
            temp_est(k) = laytemp_sis2(ml_ice, tfi(k), &
                                       sol(k) + kk*(temp_est(k - 1) + temp_est(k + 1)), &
                                       2.0_wp*kk, temp_ic(k), dtt)
         end do
         temp_est(1) = laytemp_sis2(ml_ice, tfi(1), &
                                    sol(1) + (kk*temp_est(2) + k10*temp_est(0)), &
                                    kk + k10, temp_ic(1), dtt)
      end if

      temp_est(0) = laytemp_sis2(ml_snow, 0.0_wp, &
                                 sol(0) + (k10*temp_est(1) - k0a_x_ta), &
                                 k10 + k0a, temp_ic(0), dtt)
      tsurf = (k0skin*temp_est(0) - sf_0)/(dsf_dt + k0skin)

      ! ---- Conservative DOWN pass: actually update enthalpies ----
      heat_flux_err_rat = 0.7071_wp*dtt*ICE_TEMP_RANGE_EST/ &
                          (ICE_TEMP_RANGE_EST*ICE_CP_ICE + ICE_LAT_FUS)

      e_extra_sum = 0.0_wp
      sum_sol = 0.0_wp
      do k = 0, nk
         sum_sol = sum_sol + sol(k)
      end do
      sum_sol = sum_sol*dtt

      if (tsurf > tsf .or. m_pond > 0.0_wp) then
         tsurf = tsf
         if (ml_snow > 0.0_wp) then
            heat_flux_int(-1) = k0skin*tsf
            heat_flux_int(0) = -k10*temp_est(1)
            call update_lay_enth(ml_snow, 0.0_wp, enthalpy(0), heat_flux_int(-1), &
                                 sol(0), heat_flux_int(0), -k0skin, k10, dtt, &
                                 heat_flux_err_rat, e_extra, snow_temp_new, &
                                 .false., 0.0_wp)
            tmelt = tmelt + e_extra - dtt*((sf_0 + dsf_dt*tsf) + heat_flux_int(-1))
            e_extra_sum = e_extra_sum + e_extra
            tflux_sfc = dtt*heat_flux_int(-1)
         else
            enthalpy(0) = ice_enth_from_ts(tsf, 0.0_wp)
            heat_flux_int(0) = k10*(tsf - temp_est(1))
            heat_flux_int(-1) = heat_flux_int(0)
            tmelt = tmelt + dtt*((sol(0) - (sf_0 + dsf_dt*tsf)) - heat_flux_int(0))
            tflux_sfc = dtt*heat_flux_int(0)
         end if
      else
         heat_flux_int(-1) = -k0a_x_ta
         heat_flux_int(0) = -k10*temp_est(1)
         snow_temp_max = (tsf*(dsf_dt + k0skin) + sf_0)/k0skin
         call update_lay_enth(ml_snow, 0.0_wp, enthalpy(0), heat_flux_int(-1), &
                              sol(0), heat_flux_int(0), -k0a, k10, dtt, &
                              heat_flux_err_rat, e_extra, snow_temp_new, &
                              .true., snow_temp_max)
         tsurf = (k0skin*snow_temp_new - sf_0)/(dsf_dt + k0skin)
         e_extra_sum = e_extra_sum + e_extra
         tmelt = tmelt + e_extra
         tflux_sfc = dtt*heat_flux_int(-1)
      end if

      do k = 1, nk - 1
         heat_flux_int(k) = -kk*temp_est(k + 1)
         ftop_new = heat_flux_int(k - 1)
         fbot_new = heat_flux_int(k)
         call update_lay_enth(ml_ice, sice(k), enthalpy(k), ftop_new, &
                              sol(k), fbot_new, 0.0_wp, kk, dtt, &
                              heat_flux_err_rat, e_extra, snow_temp_new, &
                              .false., 0.0_wp)
         heat_flux_int(k - 1) = ftop_new
         heat_flux_int(k) = fbot_new
         e_extra_sum = e_extra_sum + e_extra
         if (k <= nk/2) then
            tmelt = tmelt + e_extra
         else
            bmelt = bmelt + e_extra
         end if
      end do

      heat_flux_int(nk) = -2.0_wp*kk*tfw
      ftop_new = heat_flux_int(nk - 1)
      fbot_new = heat_flux_int(nk)
      call update_lay_enth(ml_ice, sice(nk), enthalpy(nk), ftop_new, &
                           sol(nk), fbot_new, 0.0_wp, 2.0_wp*kk, dtt, &
                           heat_flux_err_rat, e_extra, snow_temp_new, &
                           .false., 0.0_wp)
      heat_flux_int(nk - 1) = ftop_new
      heat_flux_int(nk) = fbot_new
      e_extra_sum = e_extra_sum + e_extra
      bmelt = bmelt + e_extra
      ! ---- END conservative update of enthalpy ----

      col_enth_out = 0.0_wp
      do k = 0, nk
         col_enth_out = col_enth_out + m_lay(k)*enthalpy(k)
      end do

      tflux_bot = -heat_flux_int(nk)*dtt

      ! TRAP #3: fb enters ONLY here, as a post-hoc bmelt residual.
      bmelt = bmelt + (dtt*fb - tflux_bot)

      ! ---- Excess-heat clamp to liq_lim (SIS2:500-524) ----
      enth_liq_lim = ice_enth_from_ts(0.0_wp, 0.0_wp)
      if (enthalpy(0) > enth_liq_lim) then
         e_extra = (enthalpy(0) - enth_liq_lim)*ml_snow
         tmelt = tmelt + e_extra
         enthalpy(0) = enth_liq_lim
      end if

      i_liq_lim = 1.0_wp/ICE_LIQ_LIM
      do k = 1, nk
         enth_liq_lim = ice_enth_from_ts(tfi(k)*i_liq_lim, sice(k))
         if (enthalpy(k) > enth_liq_lim) then
            e_extra = (enthalpy(k) - enth_liq_lim)*ml_ice
            enthalpy(k) = enth_liq_lim
            if (k <= nk/2) then
               tmelt = tmelt + e_extra
            else
               bmelt = bmelt + e_extra
            end if
         end if
      end do
   end subroutine ice_temp_sis2

   pure subroutine ice_column_step(nk, m_snow, m_ice_tot, enth_snow_pt, &
                                   enth_ice_bu, sal_ice_bu, &
                                   sf_0, dsf_dt, sw_dn, tfw, fb, sst, s_surf, dtt, &
                                   do_snow_ice, &
                                   snow, tsurf, h2o_ocn_to_ice, h2o_ice_to_ocn, &
                                   heat_to_ocn, sw_thru, snow_to_ice)
      !! Per-(cell,category) orchestrator: gather (bottom-up -> top-down
      !! flip, TRAP #2), optics, conduction (`ice_temp_sis2`), resize
      !! (snow add, bottom-freeze, top/bottom melt peel, rebalance),
      !! scatter (flip back).
      !!
      !! PR 26 TRAP: `snow` is added inside step 4 (resize), i.e. AFTER
      !! the optics (step 2) and conduction (step 3) already ran on the
      !! PRE-snowfall `m_snow`. That is SIS2's fast/slow split
      !! (`ice_resize_SIS2` runs after the conduction solve,
      !! SIS2_ice_thm.F90:1122) — new snow IS meltable in this same
      !! window (`ice_top_melt_peel` starts at k=0), but its albedo and
      !! conduction effect are felt only on the NEXT window. Do not
      !! "helpfully" move the add earlier.
      !$acc routine seq
      integer, intent(in) :: nk
         !! Number of ice layers (declared first — decl-order).
      real(wp), intent(inout) :: m_snow
         !! Snow mass per unit area (kg/m²).
      real(wp), intent(inout) :: m_ice_tot
         !! Total ice mass per unit area (kg/m²).
      real(wp), intent(inout) :: enth_snow_pt
         !! Snow specific enthalpy (J/kg).
      real(wp), intent(inout) :: enth_ice_bu(nk)
         !! BOTTOM-UP ice specific enthalpies (J/kg) — state order,
         !! `enth_ice_bu(1)` = ice bottom.
      real(wp), intent(inout) :: sal_ice_bu(nk)
         !! BOTTOM-UP ice bulk salinities (PSU) — state order.
      real(wp), intent(in) :: sf_0
         !! Linearized SEB intercept (W/m²), upward-positive.
      real(wp), intent(in) :: dsf_dt
         !! Linearized SEB slope (W/m²/K), upward-positive.
      real(wp), intent(in) :: sw_dn
         !! Downwelling shortwave at the surface (W/m²).
      real(wp), intent(in) :: tfw
         !! Seawater freezing temperature at the ice base (degC).
      real(wp), intent(in) :: fb
         !! Ocean -> ice-base heat flux (W/m²).
      real(wp), intent(in) :: sst
         !! Sea-surface temperature (degC) — feeds the TRAP-#1 liquid
         !! ocean enthalpy.
      real(wp), intent(in) :: s_surf
         !! Sea-surface salinity (PSU) — feeds the TRAP-#1 liquid
         !! ocean enthalpy (unused by the linear formula but kept for
         !! call-site parity, `ice_enthalpy_liquid`).
      real(wp), intent(in) :: dtt
         !! Timestep (s).
      logical, intent(in) :: do_snow_ice
         !! Archimedes freeboard flood gate (`&ocean_ice_nml snow_ice`,
         !! PR 27) — `.false.` (the default) is a bit-identical no-op:
         !! `snow_to_ice` stays 0 and `ice_snow_ice_flood` is never
         !! called.
      real(wp), intent(in) :: snow
         !! New snow mass this window (kg/m²), `fprec*dtt` — PR 26 source
         !! term, `ice_snow_accumulate`'s `snow` argument. 0 (the
         !! `&ocean_ice_nml snowfall=0` default) is a bit-identical no-op.
      real(wp), intent(out) :: tsurf
         !! Surface skin temperature (degC).
      real(wp), intent(out) :: h2o_ocn_to_ice
         !! Mass flux frozen from the ocean onto the ice base (kg/m²).
      real(wp), intent(out) :: h2o_ice_to_ocn
         !! Meltwater mass flux to the ocean (kg/m²), top + bottom peel.
      real(wp), intent(out) :: heat_to_ocn
         !! Leftover melt heat dumped to the ocean (J/m²), top + bottom.
      real(wp), intent(out) :: sw_thru
         !! Shortwave transmitted through the ice to the ocean (W/m²).
      real(wp), intent(out) :: snow_to_ice
         !! Mass converted from snow to the top ice layer this call
         !! (kg/m²), >= 0 — SIS2 `SN2IC` (PR 27). 0 when `do_snow_ice`
         !! is `.false.` or the column is not flooded.

      real(wp) :: enth_loc(0:ICE_NK_MAX), sal_loc(ICE_NK_MAX)
      real(wp) :: albedo, abs_sfc, abs_snow, abs_ocn, abs_int, pen
      real(wp) :: abs_ice_lay(ICE_NK_MAX), sol(0:ICE_NK_MAX)
      real(wp) :: ts_opt, sw_tot, sf_0_eff
      real(wp) :: tmelt, bmelt
      real(wp) :: col_enth_in, col_enth_out, sum_sol, tflux_sfc, tflux_bot
      real(wp) :: m_lay(0:ICE_NK_MAX), enthalpy(0:ICE_NK_MAX + 1), salin(0:ICE_NK_MAX)
      real(wp) :: enth_ocean, salin_freeze, mtot_ice
      integer :: k

      ! ---- 1. Gather + flip (TRAP #2) ----
      do k = 1, nk
         enth_loc(k) = enth_ice_bu(nk + 1 - k)
         sal_loc(k) = sal_ice_bu(nk + 1 - k)
      end do
      enth_loc(0) = enth_snow_pt
      if (m_snow == 0.0_wp) then
         ! Massless snow slot: re-seed from the top ice layer's
         ! temperature at 0 salinity (SIS_slow_thermo.F90:977).
         enth_loc(0) = ice_enth_from_ts(ice_temp_from_en_s(enth_loc(1), sal_loc(1)), 0.0_wp)
      end if

      ! ---- 2. Optics -> sol ----
      if (m_snow > 0.0_wp) then
         ts_opt = ice_temp_from_en_s(enth_loc(0), 0.0_wp)
      else
         ts_opt = ice_temp_from_en_s(enth_loc(1), sal_loc(1))
         ! TODO(PR-3b): carry a true prognostic Tskin; this reuses the
         ! top-ice-layer temperature as a skin-temp stand-in.
      end if
      call ice_optics_csim4(nk, m_snow/ICE_RHO_SNOW, m_ice_tot/ICE_RHO_ICE, &
                            ts_opt, sal_loc(1), albedo, abs_sfc, abs_snow, &
                            abs_ice_lay(1:nk), abs_ocn, abs_int, pen)

      sw_tot = (1.0_wp - albedo)*sw_dn
      sf_0_eff = sf_0 - abs_sfc*sw_tot
      sol(0) = abs_snow*sw_tot
      do k = 1, nk
         sol(k) = abs_ice_lay(k)*sw_tot
      end do
      sw_thru = abs_ocn*sw_tot

      ! ---- 3. Conduction ----
      tmelt = 0.0_wp
      bmelt = 0.0_wp
      enthalpy(0:nk) = enth_loc(0:nk)
      call ice_temp_sis2(nk, m_snow, m_ice_tot, sal_loc(1:nk), enthalpy(0:nk), &
                         sf_0_eff, dsf_dt, sol(0:nk), tfw, fb, dtt, &
                         tsurf, tmelt, bmelt, &
                         col_enth_in, col_enth_out, sum_sol, tflux_sfc, tflux_bot)
      enth_loc(0:nk) = enthalpy(0:nk)

      ! ---- 4. Resize ----
      m_lay(0) = m_snow
      do k = 1, nk
         m_lay(k) = m_ice_tot/real(nk, wp)
      end do
      enth_ocean = ice_enthalpy_liquid(sst, s_surf)  ! TRAP #1
      salin_freeze = ICE_BULK_SALINITY

      salin(0) = 0.0_wp
      salin(1:nk) = sal_loc(1:nk)
      enthalpy(0:nk) = enth_loc(0:nk)

      ! PR 26: snow source term — SIS2's ice_resize_SIS2 runs this FIRST,
      ! before the melt peels (SIS2_ice_thm.F90:1122), so new snow is
      ! meltable in this same window. `enthalpy(0)` is unchanged by the
      ! add (see ice_snow_accumulate's docstring).
      call ice_snow_accumulate(nk, m_lay(0:nk), snow)

      ! Negative-top-melt fold (do_pond=false path, SIS2:1159-1163) —
      ! unreachable in the PR-3a gates but required for SIS2 parity.
      if (tmelt < 0.0_wp) then
         bmelt = bmelt + tmelt
         tmelt = 0.0_wp
      end if

      call ice_bottom_freeze(nk, m_lay(0:nk), enthalpy(0:nk), salin(0:nk), &
                             bmelt, enth_ocean, salin_freeze, h2o_ocn_to_ice)

      heat_to_ocn = 0.0_wp
      h2o_ice_to_ocn = 0.0_wp
      call ice_top_melt_peel(nk, m_lay(0:nk), enthalpy(0:nk), salin(0:nk), tmelt, &
                             heat_to_ocn, h2o_ice_to_ocn)
      call ice_bottom_melt_peel(nk, m_lay(0:nk), enthalpy(0:nk), salin(0:nk), bmelt, &
                                heat_to_ocn, h2o_ice_to_ocn)

      ! PR 27: Archimedes freeboard flood — AFTER the melt peels (so the
      ! freshly-converted mass is not re-melted this step and m_i
      ! reflects the post-melt column) and BEFORE ice_rebalance_layers
      ! (so the new mass in layer 1 IS redistributed across the nk
      ! layers) — SIS2's exact order (TRAP #7). `do_snow_ice=.false.`
      ! (the default) is a bit-identical no-op by inspection.
      snow_to_ice = 0.0_wp
      if (do_snow_ice) then
         call ice_snow_ice_flood(nk, m_lay(0:nk), enthalpy(0:nk), salin(0:nk), &
                                 ICE_RHO_ICE/ICE_RHO_OCEAN, snow_to_ice)
      end if

      call ice_rebalance_layers(nk, m_lay(0:nk), enthalpy(0:nk), salin(0:nk), mtot_ice)

      ! ---- 5. Scatter + flip back ----
      m_ice_tot = mtot_ice
      m_snow = m_lay(0)
      enth_snow_pt = enthalpy(0)
      do k = 1, nk
         enth_ice_bu(nk + 1 - k) = enthalpy(k)
         sal_ice_bu(nk + 1 - k) = salin(k)
      end do
   end subroutine ice_column_step

   pure subroutine ice_thermo_columns(nghost, nx, ny, ncat, nk, dtt, do_snow_ice, wet_mask, &
                                      m_ice, m_snow, enth_ice, enth_snow, sal_ice, &
                                      sf_0, dsf_dt, sw_dn, fprec, tfw, fb, sst, s_surf, &
                                      tsurf_out, h2o_ocn_to_ice, h2o_ice_to_ocn, &
                                      heat_to_ocn, sw_thru, snow_to_ice)
      !! `do concurrent` cell driver: PHYSICAL cells only, inner `if`
      !! gate (never a masked DC header), serial `do cat` loop inside.
      !! Per (i,j,cat): outputs zeroed unconditionally, then gated on
      !! `wet_mask > 0.5 .and. m_ice > ICE_RHO_ICE*H_VANISHED` (dynamic-
      !! vanish taxonomy: skip intact, never clamp/divide a vanished
      !! column). `part_size` is deliberately NOT an argument — thermo
      !! is per unit ice area.
      integer, intent(in) :: nghost, nx, ny, ncat, nk
         !! Grid + category + layer extents (declared first — decl-order).
      real(wp), intent(in) :: dtt
         !! Timestep (s).
      logical, intent(in) :: do_snow_ice
         !! Archimedes freeboard flood gate (`&ocean_ice_nml snow_ice`,
         !! PR 27), captured by value into the DC loop (same treatment
         !! as `dtt`). `.false.` (the default) is a bit-identical no-op.
      real(wp), intent(in) :: wet_mask(nx, ny)
         !! Ocean wet mask (>0.5 = wet).
      real(wp), intent(inout) :: m_ice(nx, ny, ncat)
         !! Total ice mass per unit area per category (kg/m²).
      real(wp), intent(inout) :: m_snow(nx, ny, ncat)
         !! Snow mass per unit area per category (kg/m²).
      real(wp), intent(inout) :: enth_ice(nx, ny, ncat, nk)
         !! Ice specific enthalpy (J/kg), BOTTOM-UP (k=1 = ice bottom).
      real(wp), intent(inout) :: enth_snow(nx, ny, ncat, 1)
         !! Snow specific enthalpy (J/kg).
      real(wp), intent(inout) :: sal_ice(nx, ny, ncat, nk)
         !! Ice bulk salinity (PSU), BOTTOM-UP.
      real(wp), intent(in) :: sf_0(nx, ny)
         !! Linearized SEB intercept (W/m²), upward-positive.
      real(wp), intent(in) :: dsf_dt(nx, ny)
         !! Linearized SEB slope (W/m²/K), upward-positive.
      real(wp), intent(in) :: sw_dn(nx, ny)
         !! Downwelling shortwave at the surface (W/m²).
      real(wp), intent(in) :: fprec(nx, ny)
         !! Frozen-precipitation rate onto the ice top (kg/m²/s), >= 0 —
         !! PR 26 snowfall seam (`ice%atm_fprec`). Passed to
         !! `ice_column_step` as `fprec(i,j)*dtt`; 0 (the
         !! `&ocean_ice_nml snowfall=0` default) is a bit-identical no-op.
      real(wp), intent(in) :: tfw(nx, ny)
         !! Seawater freezing temperature at the ice base (degC).
      real(wp), intent(in) :: fb(nx, ny)
         !! Ocean -> ice-base heat flux (W/m²).
      real(wp), intent(in) :: sst(nx, ny)
         !! Sea-surface temperature (degC).
      real(wp), intent(in) :: s_surf(nx, ny)
         !! Sea-surface salinity (PSU).
      real(wp), intent(inout) :: tsurf_out(nx, ny, ncat)
         !! Surface skin temperature (degC).
      real(wp), intent(inout) :: h2o_ocn_to_ice(nx, ny, ncat)
         !! Mass flux frozen from the ocean onto the ice base (kg/m²).
      real(wp), intent(inout) :: h2o_ice_to_ocn(nx, ny, ncat)
         !! Meltwater mass flux to the ocean (kg/m²).
      real(wp), intent(inout) :: heat_to_ocn(nx, ny, ncat)
         !! Leftover melt heat dumped to the ocean (J/m²).
      real(wp), intent(inout) :: sw_thru(nx, ny, ncat)
         !! Shortwave transmitted through the ice to the ocean (W/m²).
      real(wp), intent(inout) :: snow_to_ice(nx, ny, ncat)
         !! Mass converted from snow to the top ice layer this call
         !! (kg/m²), >= 0 — SIS2 `SN2IC` (PR 27). Zeroed unconditionally
         !! every thermo step, same lifecycle as `sw_thru`.

      integer :: i, j, cat, i_lo, i_hi, j_lo, j_hi
      real(wp) :: m_snow_pt, m_ice_pt, enth_snow_pt
      real(wp) :: enth_ice_col(ICE_NK_MAX), sal_ice_col(ICE_NK_MAX)
      real(wp) :: tsurf_pt, h2o_ocn_to_ice_pt, h2o_ice_to_ocn_pt, heat_to_ocn_pt, sw_thru_pt
      real(wp) :: snow_to_ice_pt

      i_lo = nghost + 1
      i_hi = nx - nghost
      j_lo = nghost + 1
      j_hi = ny - nghost

      do concurrent(j=j_lo:j_hi, i=i_lo:i_hi) &
         local(cat, m_snow_pt, m_ice_pt, enth_snow_pt, enth_ice_col, sal_ice_col, &
               tsurf_pt, h2o_ocn_to_ice_pt, h2o_ice_to_ocn_pt, heat_to_ocn_pt, sw_thru_pt, &
               snow_to_ice_pt)
         do cat = 1, ncat
            tsurf_out(i, j, cat) = 0.0_wp
            h2o_ocn_to_ice(i, j, cat) = 0.0_wp
            h2o_ice_to_ocn(i, j, cat) = 0.0_wp
            heat_to_ocn(i, j, cat) = 0.0_wp
            sw_thru(i, j, cat) = 0.0_wp
            snow_to_ice(i, j, cat) = 0.0_wp

            if (wet_mask(i, j) > 0.5_wp .and. m_ice(i, j, cat) > ICE_RHO_ICE*H_VANISHED) then
               m_snow_pt = m_snow(i, j, cat)
               m_ice_pt = m_ice(i, j, cat)
               enth_snow_pt = enth_snow(i, j, cat, 1)
               enth_ice_col(1:nk) = enth_ice(i, j, cat, 1:nk)
               sal_ice_col(1:nk) = sal_ice(i, j, cat, 1:nk)

               call ice_column_step(nk, m_snow_pt, m_ice_pt, enth_snow_pt, &
                                    enth_ice_col(1:nk), sal_ice_col(1:nk), &
                                    sf_0(i, j), dsf_dt(i, j), sw_dn(i, j), &
                                    tfw(i, j), fb(i, j), sst(i, j), s_surf(i, j), dtt, &
                                    do_snow_ice, &
                                    fprec(i, j)*dtt, tsurf_pt, h2o_ocn_to_ice_pt, &
                                    h2o_ice_to_ocn_pt, heat_to_ocn_pt, sw_thru_pt, &
                                    snow_to_ice_pt)

               m_snow(i, j, cat) = m_snow_pt
               m_ice(i, j, cat) = m_ice_pt
               enth_snow(i, j, cat, 1) = enth_snow_pt
               enth_ice(i, j, cat, 1:nk) = enth_ice_col(1:nk)
               sal_ice(i, j, cat, 1:nk) = sal_ice_col(1:nk)
               tsurf_out(i, j, cat) = tsurf_pt
               h2o_ocn_to_ice(i, j, cat) = h2o_ocn_to_ice_pt
               h2o_ice_to_ocn(i, j, cat) = h2o_ice_to_ocn_pt
               heat_to_ocn(i, j, cat) = heat_to_ocn_pt
               sw_thru(i, j, cat) = sw_thru_pt
               snow_to_ice(i, j, cat) = snow_to_ice_pt
            end if
         end do
      end do
   end subroutine ice_thermo_columns

end module rdb_ice_column
