!! Sea-ice column mass bookkeeping (SIS2 port, PR 3a).
module rdb_ice_mass
   !! Ports of the subset of `ice_resize_SIS2` (SIS2_ice_thm.F90:1010-1338)
   !! needed by the Winton column step: bottom freezing, top/bottom melt
   !! peel, and equal-mass layer rebalancing
   !! (`rebalance_ice_layers`, SIS2_ice_thm.F90:1448-1512). Ported from
   !! the validated prototype `sis2_resize.py` (Apache-2.0 source
   !! attribution as above).
   !!
   !! PR 26 ports the snow source term (`ice_snow_accumulate`, below) —
   !! the first branch SIS2 runs in `ice_resize_SIS2`
   !! (SIS2_ice_thm.F90:1122), so it sits first in this module too.
   !! PR 27 ports the Archimedes freeboard snow-ice flooding conversion
   !! (`ice_snow_ice_flood`, below) — the FINAL substantive block of
   !! `ice_resize_SIS2` (SIS2_ice_thm.F90:1303-1320), so it sits last in
   !! this module too. Pond/evap/rain mass paths remain deliberately NOT
   !! ported. TODO(PR-35): pond mass paths. Evap/rain need atmospheric
   !! fields Roundabout does not carry (PR-55 at the earliest).
   !!
   !! All procedures pure, `!$acc routine seq`, explicit-shape,
   !! TOP-DOWN columns (index 0 = snow, 1..nk = ice top->bottom) —
   !! identical convention to `rdb_ice_column`.
   use rdb_constants, only: wp
   use rdb_ice_enthalpy, only: ICE_LAT_FUS, ICE_LIQ_LIM, ICE_NK_MAX, &
                               ice_enthalpy_liquid_freeze
   implicit none
   private

   public :: ice_snow_accumulate
   public :: ice_bottom_freeze, ice_top_melt_peel, ice_bottom_melt_peel, &
             ice_rebalance_layers, ice_snow_ice_flood

contains

   pure subroutine ice_snow_accumulate(nk, m_lay, snow)
      !! Snow source-term branch of `ice_resize_SIS2`
      !! (SIS2_ice_thm.F90:1122): `m_lay(0) = m_lay(0) + snow`. The
      !! snow layer's specific enthalpy `enthalpy(0)` is UNCHANGED by
      !! snowfall -- SIS2 is explicit that this "should do nothing"
      !! (SIS_slow_thermo.F90:1032-1033) and books the implied energy
      !! against the atmosphere (`enth_snowfall = snow*enthalpy(0)` ->
      !! `Enth_Mass_in_atm`). Roundabout's atmosphere is a prescribed slab
      !! with no energy budget to charge, so that bookkeeping is not
      !! ported -- a documented v1 divergence (PLAN_PR26_snowfall.md, S3).
      !!
      !! The SIS2 guard `if (mtot_ice == 0.0) m_lay(0) = 0.0`
      !! (SIS2_ice_thm.F90:1121, needed because SIS2 gates its column on
      !! `part_size > 0` and so can reach a category with zero ice mass)
      !! is deliberately NOT ported: `ice_column_step`'s caller only
      !! ever reaches this routine on a category that already passed
      !! the column's own entry gate (`m_ice > ICE_RHO_ICE*H_VANISHED`,
      !! `rdb_ice_column.F90`), so `mtot_ice > 0` always holds here. The
      !! ice-free share of the snowfall is routed to the ocean by
      !! `rdb_ice_snow%ice_snowfall_ocean_share` instead of being added
      !! as orphan snow (which would trip `rdb_ice_transport`'s
      !! fail-loud `mca_snow > 0` where `mca_ice <= 0` reduction).
      !$acc routine seq
      integer, intent(in) :: nk
         !! Number of ice layers (declared first -- decl-order).
      real(wp), intent(inout) :: m_lay(0:nk)
         !! Layer masses (kg/m^2), 0 = snow, 1..nk = ice.
      real(wp), intent(in) :: snow
         !! New snow mass this thermo window (kg/m^2), `fprec*dt_therm`.
         !! `snowfall=0` => `snow==0.0_wp` => the guard below makes this
         !! call bit-identical to a no-op by inspection, not IEEE luck.

      if (snow /= 0.0_wp) then
         m_lay(0) = m_lay(0) + snow
      end if
   end subroutine ice_snow_accumulate

   pure subroutine ice_bottom_freeze(nk, m_lay, enthalpy, salin, bmelt, &
                                     enth_ocean, salin_freeze, h2o_ocn_to_ice)
      !! Bottom-freezing branch of `ice_resize_SIS2`
      !! (SIS2_ice_thm.F90:1164-1188; prototype sis2_resize.py:20-45).
      !! When `bmelt < 0` (net upward heat deficit at the base), freeze
      !! ocean water onto the bottom ice layer: `enth_freeze =
      !! min(enthalpy(nk), enth_ocean - min_denth_freeze)`,
      !! `m_freeze = -bmelt/(enth_ocean - enth_freeze)`, mass-weighted
      !! mix of the bottom layer's enthalpy + salinity, `bmelt` reset
      !! to 0.
      !$acc routine seq
      integer, intent(in) :: nk
         !! Number of ice layers (declared first — decl-order).
      real(wp), intent(inout) :: m_lay(0:nk)
         !! Layer masses (kg/m²), 0 = snow, 1..nk = ice.
      real(wp), intent(inout) :: enthalpy(0:nk)
         !! Layer specific enthalpies (J/kg).
      real(wp), intent(inout) :: salin(0:nk)
         !! Layer bulk salinities (PSU).
      real(wp), intent(inout) :: bmelt
         !! Accumulated bottom melting/freezing energy (J/m²); reset
         !! to 0 on exit when freezing occurred.
      real(wp), intent(in) :: enth_ocean
         !! Ocean-water specific enthalpy at the ice base (J/kg) —
         !! TRAP #1: the LIQUID formula `ice_enthalpy_liquid(sst,
         !! s_surf)`, never the frozen/mushy `ice_enth_from_ts`.
      real(wp), intent(in) :: salin_freeze
         !! Salinity of newly frozen ice (PSU) — bulk-salinity mode:
         !! `ICE_BULK_SALINITY`.
      real(wp), intent(out) :: h2o_ocn_to_ice
         !! Mass flux frozen from the ocean onto the ice base (kg/m²).

      real(wp) :: min_denth_freeze, enth_freeze, m_freeze

      h2o_ocn_to_ice = 0.0_wp
      if (bmelt < 0.0_wp) then
         min_denth_freeze = ICE_LAT_FUS*(1.0_wp - ICE_LIQ_LIM)
         enth_freeze = min(enthalpy(nk), enth_ocean - min_denth_freeze)
         m_freeze = -bmelt/(enth_ocean - enth_freeze)

         enthalpy(nk) = (m_lay(nk)*enthalpy(nk) + m_freeze*enth_freeze)/ &
                        (m_lay(nk) + m_freeze)
         salin(nk) = (m_lay(nk)*salin(nk) + m_freeze*salin_freeze)/ &
                     (m_lay(nk) + m_freeze)
         m_lay(nk) = m_lay(nk) + m_freeze
         h2o_ocn_to_ice = m_freeze
         bmelt = 0.0_wp
      end if
   end subroutine ice_bottom_freeze

   pure subroutine ice_top_melt_peel(nk, m_lay, enthalpy, salin, tmelt, &
                                     heat_to_ocn, h2o_ice_to_ocn)
      !! Top melt peel (SIS2_ice_thm.F90:1217-1242; prototype
      !! sis2_resize.py:48-80). Peels mass from k=0 (snow) upward
      !! through k=nk until `tmelt` is spent; a massless layer is
      !! skipped; the partial layer takes `m_melt = melt_left/(enth_fr
      !! - enthalpy)`; any leftover melt energy (all layers exhausted)
      !! drains to `heat_to_ocn`. Snow/pond-free path — `tmelt` is
      !! assumed already >= 0 on entry (the caller folds any negative
      !! top-melt into `bmelt` upstream, SIS2:1159-1163).
      !$acc routine seq
      integer, intent(in) :: nk
         !! Number of ice layers (declared first — decl-order).
      real(wp), intent(inout) :: m_lay(0:nk)
         !! Layer masses (kg/m²), 0 = snow, 1..nk = ice.
      real(wp), intent(inout) :: enthalpy(0:nk)
         !! Layer specific enthalpies (J/kg).
      real(wp), intent(in) :: salin(0:nk)
         !! Layer bulk salinities (PSU) — index 0 (snow) unused.
      real(wp), intent(in) :: tmelt
         !! Accumulated top melting energy (J/m²), assumed >= 0.
      real(wp), intent(inout) :: heat_to_ocn
         !! Leftover melt energy after all layers exhausted (J/m²) —
         !! accumulator, caller zeroes once per step.
      real(wp), intent(inout) :: h2o_ice_to_ocn
         !! Meltwater mass flux to the ocean (kg/m²) — accumulator.

      real(wp) :: enth_fr(0:ICE_NK_MAX)
      real(wp) :: melt_left, avail, m_melt
      integer :: k

      enth_fr(0) = ice_enthalpy_liquid_freeze(0.0_wp)
      do k = 1, nk
         enth_fr(k) = ice_enthalpy_liquid_freeze(salin(k))
      end do

      melt_left = tmelt
      if (melt_left > 0.0_wp) then
         do k = 0, nk
            if (m_lay(k) <= 0.0_wp) cycle
            avail = m_lay(k)*(enth_fr(k) - enthalpy(k))
            if (melt_left < avail) then
               m_melt = melt_left/(enth_fr(k) - enthalpy(k))
               melt_left = 0.0_wp
            else
               m_melt = m_lay(k)
               melt_left = melt_left - avail
            end if
            m_lay(k) = m_lay(k) - m_melt
            h2o_ice_to_ocn = h2o_ice_to_ocn + m_melt
            if (melt_left <= 0.0_wp) exit
         end do
         heat_to_ocn = heat_to_ocn + melt_left
      end if
   end subroutine ice_top_melt_peel

   pure subroutine ice_bottom_melt_peel(nk, m_lay, enthalpy, salin, bmelt, &
                                        heat_to_ocn, h2o_ice_to_ocn)
      !! Bottom melt peel (SIS2_ice_thm.F90:1246-1271; prototype
      !! sis2_resize.py:83-111). Same peel as `ice_top_melt_peel` but
      !! from k=nk down to k=0. The prototype's separate `ablation`
      !! return is dropped (it is `h2o_ice_to_ocn` restricted to this
      !! call; PR 3b can re-derive it if needed).
      !$acc routine seq
      integer, intent(in) :: nk
         !! Number of ice layers (declared first — decl-order).
      real(wp), intent(inout) :: m_lay(0:nk)
         !! Layer masses (kg/m²), 0 = snow, 1..nk = ice.
      real(wp), intent(inout) :: enthalpy(0:nk)
         !! Layer specific enthalpies (J/kg).
      real(wp), intent(in) :: salin(0:nk)
         !! Layer bulk salinities (PSU) — index 0 (snow) unused.
      real(wp), intent(inout) :: bmelt
         !! Accumulated bottom melting energy (J/m²), consumed by the
         !! peel; NOTE unlike `ice_bottom_freeze` this is intent
         !! inout only for interface symmetry — the peel does not
         !! reset it (caller passes the post-freeze residual).
      real(wp), intent(inout) :: heat_to_ocn
         !! Leftover melt energy after all layers exhausted (J/m²) —
         !! accumulator, caller zeroes once per step.
      real(wp), intent(inout) :: h2o_ice_to_ocn
         !! Meltwater mass flux to the ocean (kg/m²) — accumulator.

      real(wp) :: enth_fr(0:ICE_NK_MAX)
      real(wp) :: melt_left, avail, m_melt
      integer :: k

      enth_fr(0) = ice_enthalpy_liquid_freeze(0.0_wp)
      do k = 1, nk
         enth_fr(k) = ice_enthalpy_liquid_freeze(salin(k))
      end do

      melt_left = bmelt
      if (melt_left > 0.0_wp) then
         do k = nk, 0, -1
            if (m_lay(k) <= 0.0_wp) cycle
            avail = m_lay(k)*(enth_fr(k) - enthalpy(k))
            if (melt_left < avail) then
               m_melt = melt_left/(enth_fr(k) - enthalpy(k))
               melt_left = 0.0_wp
            else
               m_melt = m_lay(k)
               melt_left = melt_left - avail
            end if
            m_lay(k) = m_lay(k) - m_melt
            h2o_ice_to_ocn = h2o_ice_to_ocn + m_melt
            if (melt_left <= 0.0_wp) exit
         end do
         heat_to_ocn = heat_to_ocn + melt_left
      end if
   end subroutine ice_bottom_melt_peel

   pure subroutine ice_rebalance_layers(nk, m_lay, enthalpy, salin, mtot_ice)
      !! Equal-mass repartition of ice layers 1..nk, mass-weighting
      !! enthalpy and salinity (SIS2_ice_thm.F90:1448-1512; prototype
      !! sis2_resize.py:114-161). The snow slot (`m_lay(0)`,
      !! `enthalpy(0)`) is untouched. `mtot_ice = sum(m_lay(1:nk))` is
      !! returned; `mtot_ice == 0` is an early exit leaving the ice
      !! layers untouched (there is no ice to rebalance). The k1/k2
      !! two-pointer drain loop follows the SIS2/prototype branch
      !! order exactly, including the
      !! `(m_ice_avg - mlay_new(k2) > src_m(k1)) .or. (k2 == nk)` test.
      !$acc routine seq
      integer, intent(in) :: nk
         !! Number of ice layers (declared first — decl-order).
      real(wp), intent(inout) :: m_lay(0:nk)
         !! Layer masses (kg/m²), 0 = snow (untouched), 1..nk = ice.
      real(wp), intent(inout) :: enthalpy(0:nk)
         !! Layer specific enthalpies (J/kg).
      real(wp), intent(inout) :: salin(0:nk)
         !! Layer bulk salinities (PSU).
      real(wp), intent(out) :: mtot_ice
         !! Summed ice mass (kg/m²), `sum(m_lay(1:nk))` on entry.

      real(wp) :: mlay_new(ICE_NK_MAX), enth_ice_new(ICE_NK_MAX), sal_ice_new(ICE_NK_MAX)
      real(wp) :: src_m(ICE_NK_MAX)
      real(wp) :: m_ice_avg, m_transfer
      integer :: k1, k2, k

      mtot_ice = 0.0_wp
      do k = 1, nk
         mtot_ice = mtot_ice + m_lay(k)
      end do
      if (mtot_ice == 0.0_wp) return

      do k = 1, nk
         mlay_new(k) = 0.0_wp
         enth_ice_new(k) = 0.0_wp
         sal_ice_new(k) = 0.0_wp
         src_m(k) = m_lay(k)
      end do
      m_ice_avg = mtot_ice/real(nk, wp)

      k1 = 1
      k2 = 1
      do
         if (mlay_new(k2) >= m_ice_avg .and. k2 < nk) then
            k2 = k2 + 1
         else if (src_m(k1) <= 0.0_wp) then
            k1 = k1 + 1
         else if ((m_ice_avg - mlay_new(k2) > src_m(k1)) .or. (k2 == nk)) then
            m_transfer = src_m(k1)
            enth_ice_new(k2) = enth_ice_new(k2) + m_transfer*enthalpy(k1)
            sal_ice_new(k2) = sal_ice_new(k2) + m_transfer*salin(k1)
            mlay_new(k2) = mlay_new(k2) + m_transfer
            src_m(k1) = 0.0_wp
            k1 = k1 + 1
         else
            m_transfer = m_ice_avg - mlay_new(k2)
            enth_ice_new(k2) = enth_ice_new(k2) + m_transfer*enthalpy(k1)
            sal_ice_new(k2) = sal_ice_new(k2) + m_transfer*salin(k1)
            mlay_new(k2) = m_ice_avg
            src_m(k1) = src_m(k1) - m_transfer
            k2 = k2 + 1
         end if
         if (k1 > nk) exit
      end do

      do k = 1, nk
         if (mlay_new(k) > 0.0_wp) then
            enthalpy(k) = enth_ice_new(k)/mlay_new(k)
            salin(k) = sal_ice_new(k)/mlay_new(k)
         end if
         m_lay(k) = mlay_new(k)
      end do
   end subroutine ice_rebalance_layers

   pure subroutine ice_snow_ice_flood(nk, m_lay, enthalpy, salin, rho_ratio, snow_to_ice)
      !! Archimedes freeboard snow-ice flooding — the FINAL substantive
      !! block of `ice_resize_SIS2` (SIS2_ice_thm.F90:1303-1320; PR 27).
      !! Standard closure: Leppäranta (1983), *A growth model for black
      !! ice, snow ice and snow thickness in subarctic basins*, Nordic
      !! Hydrology 14, 59-70; Fichefet & Morales Maqueda (1997), JGR 102,
      !! 12609-12646 §2.3.
      !!
      !! The column floats, displacing its own mass of seawater. If the
      !! ice alone cannot support the snow load (`m_submerged =
      !! (m_i+m_s)*rho_ratio > m_i`, `rho_ratio = ICE_RHO_ICE/
      !! ICE_RHO_OCEAN`), the snow-ice interface is below the waterline
      !! ("flooded"): convert `snow_to_ice = min(m_submerged - m_i,
      !! m_lay(0))` kg/m^2 of snow into the TOP ice layer (local index 1
      !! — TOP-DOWN column, TRAP #2 in `rdb_ice_column`'s module
      !! docstring; the caller flips this to the state's bottom-up
      !! `enth_ice(...,nk)` at the scatter boundary). One non-iterative
      !! step suffices: the conversion is 1:1 in mass, so `m_i+m_lay(0)`
      !! is invariant and `m_submerged` does not move, so the interface
      !! lands EXACTLY on the waterline (see PLAN_PR27_snow_ice_
      !! flooding.md §3 for the algebraic proof + a worked golden).
      !!
      !! Roundabout has no ponds (`m_pond` is a hardwired dead local in
      !! `rdb_ice_column`), so the `min(...)` clamp is PROVABLY DEAD
      !! CODE here: it binds iff `m_lay(0) < -m_i`, impossible for
      !! non-negative masses (see the plan §9.2 sweep test). Kept for
      !! SIS2 parity — ponds would resurrect it.
      !!
      !! Mass-, enthalpy- and salt-conserving WITHIN THE COLUMN, exactly
      !! (`(m_0-s)*E_0 + (m_1+s)*[(m_1*E_1+s*E_0)/(m_1+s)] = m_0*E_0 +
      !! m_1*E_1`, identity; the salinity dilution `S_1*m_1/(m_1+s)` IS
      !! the mass-weighted mix since `S_snow == 0`, `salin(0)`).
      !! `enthalpy(1)` and `salin(1)` both DECREASE — snow's cold, fresh
      !! mass lands in the ice.
      !!
      !! **Known physical incompleteness, inherited from SIS2 (not fixed
      !! here — true seawater flooding is a materially larger, separate
      !! closure, see the plan's §12/§14 Q1).** Real flooding draws
      !! SEAWATER into the pore space, which then refreezes, releasing
      !! latent heat and rejecting brine. SIS2 instead moves snow mass
      !! with the SNOW's enthalpy and ZERO salinity — snow enthalpy
      !! (<= -L_f ~ -3.34e5 J/kg) is far more negative than near-freezing
      !! seawater, so the new ice is too COLD; `S_snow == 0` means it is
      !! too FRESH. `snow_to_ice` is diagnostic-only (SIS2's SN2IC):
      !! flooding exchanges NOTHING with the ocean (SIS2_ice_thm.F90:1273,
      !! "There are no further heat or mass losses or gains by the
      !! ice+snow").
      !$acc routine seq
      integer, intent(in) :: nk
         !! Number of ice layers (declared first — decl-order).
      real(wp), intent(inout) :: m_lay(0:nk)
         !! Layer masses (kg/m²), 0 = snow, 1..nk = ice.
      real(wp), intent(inout) :: enthalpy(0:nk)
         !! Layer specific enthalpies (J/kg).
      real(wp), intent(inout) :: salin(0:nk)
         !! Layer bulk salinities (PSU).
      real(wp), intent(in) :: rho_ratio
         !! `ICE_RHO_ICE/ICE_RHO_OCEAN` (nondim, < 1) — computed by the
         !! caller (`rdb_ice_column`), which already has both density
         !! parameters in scope; keeps this module free of a circular
         !! `use rdb_ice_column` (that module `use`s this one).
      real(wp), intent(out) :: snow_to_ice
         !! Mass converted from snow to the top ice layer this call
         !! (kg/m²), >= 0 — SIS2 `SN2IC`.

      real(wp) :: m_i, m_submerged

      m_i = sum(m_lay(1:nk))
      m_submerged = (m_i + m_lay(0))*rho_ratio

      snow_to_ice = 0.0_wp
      if (m_submerged > m_i) then
         snow_to_ice = min(m_submerged - m_i, m_lay(0))
      end if

      ! Denominator `m_lay(1) + snow_to_ice` is provably > 0 whenever
      ! this branch is taken: entering it forces `m_lay(0) > 0` (since
      ! rho_ratio < 1, m_submerged > m_i needs m_lay(0) > 0), hence
      ! snow_to_ice = min(m_lay(0), ...) > 0 — even if m_lay(1) == 0
      ! (top layer melted out this step). The inner guard below is free
      ! and makes the no-snow no-op provable by inspection; it is NOT
      ! `H_DIV_EPS` armour (CLAUDE.md taxonomy — the branch guard above
      ! already rules out the zero case).
      if (snow_to_ice > 0.0_wp) then
         m_lay(0) = m_lay(0) - snow_to_ice
         enthalpy(1) = (m_lay(1)*enthalpy(1) + snow_to_ice*enthalpy(0))/ &
                       (m_lay(1) + snow_to_ice)
         salin(1) = salin(1)*m_lay(1)/(m_lay(1) + snow_to_ice)
         m_lay(1) = m_lay(1) + snow_to_ice
      end if
   end subroutine ice_snow_ice_flood

end module rdb_ice_mass
