!! Frazil-heat-bank spend: new-ice formation from the ocean-side
!! supercooling bank (SIS2 port, PR 3b).
module rdb_ice_frazil_uptake
   !! Port of `add_frazil_SIS2` (SIS2_ice_thm.F90:1342-1444, Apache-2.0),
   !! bulk-salinity mode only. Spends the WHOLE per-cell frazil bank
   !! (`ocean_sea_ice_t%frazil_heat`, banked by `rdb_ice_frazil`'s
   !! surface-freezing clamp) once per thermo step by depositing new ice
   !! mass, evenly split across the `nk_ice` layers of category 1
   !! (`cat = 1` — the only category this v1 model fills; `part_size` is
   !! untouched, same convention as PR 3a).
   !!
   !! **TRAP #1** (module docstring precedent: `rdb_ice_column`,
   !! `rdb_ice_mass%ice_bottom_freeze`): the ocean-water specific
   !! enthalpy consumed by the freeze is the LIQUID form
   !! `enth_ocean = ice_enthalpy_liquid(sst, s_surf)`, never the
   !! frozen/mushy `ice_enth_from_ts` — see `ice_frazil_uptake_column`.
   !!
   !! **Index-flip discipline** (PR 3a TRAP #2 precedent): the state
   !! arrays `enth_ice`/`sal_ice` are BOTTOM-UP (`k=1` = ice bottom); the
   !! SIS2 algorithm spends the bank TOP-DOWN. The flip happens ONLY at
   !! the gather/scatter boundary inside `ice_frazil_uptake_column` —
   !! nothing in the per-layer spend loop ever sees the bottom-up
   !! convention.
   !!
   !! **Two-mode convention (PR 4a, `rdb_ice_state` module docstring)**:
   !! `ncat == 1` — legacy lumped mode, dispatched to the EXISTING
   !! `ice_frazil_uptake_impl` byte-for-byte UNCHANGED: ALL frazil ice
   !! deposits into category 1, treated as per unit CELL area (kg per m²
   !! of cell) — self-consistent with the frazil bank (J per m² of cell)
   !! and the brine-rejection flux. `ncat > 1` — SIS2 ITD mode, dispatched
   !! to `ice_frazil_uptake_multicat_impl`: the bank ANNEXES open water
   !! into the thinnest occupied-or-empty category (`k_merge`, SIS2
   !! `SIS_slow_thermo.F90:1121-1146`) via an area-weighted dilution at
   !! constant mass (`part(k_merge) += part(0)`, thickness drops, mass
   !! doesn't — the SIS2 area-creation move), then spends the bank on
   !! that category's column PER UNIT ICE AREA
   !! (`frazil_col = frazil_heat/part(k_merge)`, SIS_slow_thermo.F90:
   !! 1181-1186) via the SAME UNCHANGED `ice_frazil_uptake_column`. SIS2's
   !! default `SIS2_FILLING_FRAZIL=.true.` thin-category-fill mode (which
   !! would spread new ice across MULTIPLE thin categories instead of one
   !! `k_merge`) is NOT ported — v1 always merges into a single category,
   !! documented divergence.
   !!
   !! **Salt bookkeeping** (Boussinesq virtual-flux convention): freezing
   !! `m_frozen` kg/m² of seawater at salinity `s_surf` into ice that
   !! keeps `salt_to_ice = m_frozen * ICE_BULK_SALINITY` yields a diag
   !! rate `salt_flux_diag = (m_frozen*s_surf - salt_to_ice)/dt_therm`
   !! [PSU·kg/m²/s]. The ocean's water MASS is not reduced (volume-
   !! conserving virtual-salt-flux ocean — MOM6/SIS2 default); the
   !! closed-budget identity is exact by construction:
   !! `rho0*Delta(sum_k hTr_S) + salt_into_ice == m_frozen*s_surf` per
   !! cell (`rho0 = sf%rho0`; the surface-flux apply gives
   !! `Delta(hS) = Q_salt*dt_therm/rho0`). Freshwater/mass coupling is
   !! PR-3c+ territory. NOTE: when `s_surf < ICE_BULK_SALINITY` the
   !! brine flux goes NEGATIVE (freezing freshens the ocean locally) —
   !! this is SIS2-faithful and deliberately not clamped.
   !!
   !! **Energy accounting**: Q_heat is NOT written by this module. The
   !! frazil latent heat was already credited to the ocean by the PR-1
   !! surface clamp (which warmed the surface to T_f when it banked the
   !! deficit); spending the bank here as ice latent heat is the closing
   !! half of that exchange — writing a Q_heat here would double-count.
   !! Melt-side Q_heat/Q_salt arrive with PR 3c.
   !!
   !! Mirrors the outer-shim + flat-impl + `!$acc routine seq`
   !! column-worker structure of `rdb_ice_frazil` /
   !! `rdb_ice_column%ice_thermo_columns`.
   use rdb_constants, only: wp, H_VANISHED
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_eos, only: eos_t, eos_freezing_point
   use rdb_ice_enthalpy, only: ICE_LAT_FUS, ICE_LIQ_LIM, ice_enthalpy_liquid, &
                               ice_enth_from_ts, ice_t_freeze
   use rdb_ice_column, only: ICE_BULK_SALINITY, ICE_NK_MAX
   use rdb_ice_mass, only: ice_rebalance_layers
   use rdb_ice_state, only: ocean_sea_ice_t
   implicit none
   private

   public :: ice_frazil_uptake

   real(wp), parameter, public :: ICE_FRAZIL_T_OFFSET = 0.5_wp
      !! SIS2 `frazil_temp_offset` default (degC) — SIS2_ice_thm.F90:130.
      !! The per-layer frazil crystal forms `ICE_FRAZIL_T_OFFSET` degC
      !! BELOW the local layer freezing point, matching the observed
      !! slight supercooling of newly nucleated frazil ice.

contains

   pure subroutine ice_frazil_uptake(grid, eos, ms, ice, dt_therm)
      !! Outer shim (outer-shim + flat-impl pattern): dereference the
      !! tracer registry (`ms%tracers(idx)%hTr`) on the HOST and forward
      !! bare arrays to the device kernel — NVHPC stdpar cannot follow
      !! the array-of-derived-types indirection inside a do-concurrent
      !! body (same rule as `ice_frazil_accumulate`). No-op when either
      !! S or T is unregistered.
      !!
      !! PR 4a dispatch (module docstring two-mode convention): `ncat==1`
      !! calls the EXISTING `ice_frazil_uptake_impl` byte-for-byte
      !! UNCHANGED (bit-identity with PR 3c); `ncat>1` calls the new
      !! `ice_frazil_uptake_multicat_impl`. The branch lives ONLY here —
      !! never inside a kernel.
      type(hgrid_t), intent(in) :: grid
      type(eos_t), intent(in) :: eos
      type(multilayer_state_t), intent(in) :: ms
         !! READ-ONLY: SST/SSS are SAMPLED here, never written — the
         !! ocean-side effect of the uptake (brine rejection) arrives
         !! only via `Q_salt` on the NEXT thermo window, through
         !! `ice_ocean_brine_flux` + `ocean_surface_flux_apply_tracers`.
      type(ocean_sea_ice_t), intent(inout) :: ice
      real(wp), intent(in) :: dt_therm
         !! Effective thermo timestep (s) — `ocean_dyn_t%therm_dt(dt)`.

      integer :: idx_T, idx_S

      idx_T = ms%idx_temperature
      idx_S = ms%idx_salinity
      if (idx_T <= 0 .or. idx_S <= 0) return

      if (ice%ncat == 1) then
         call ice_frazil_uptake_impl(ms%tracers(idx_T)%hTr, ms%tracers(idx_S)%hTr, &
                                     ms%h_layer, ms%wet_mask, eos, &
                                     ice%frazil_heat, ice%m_ice, ice%m_snow, &
                                     ice%enth_ice, ice%sal_ice, &
                                     ice%m_frozen_diag, ice%salt_flux_diag, &
                                     dt_therm, grid%nghost, ice%ncat, ice%nk_ice, &
                                     ms%nz_ml, grid%nx_total, grid%ny_total)
      else
         call ice_frazil_uptake_multicat_impl(ms%tracers(idx_T)%hTr, ms%tracers(idx_S)%hTr, &
                                              ms%h_layer, ms%wet_mask, eos, &
                                              ice%frazil_heat, ice%part_size, ice%m_ice, &
                                              ice%m_snow, ice%enth_ice, ice%sal_ice, &
                                              ice%m_frozen_diag, ice%salt_flux_diag, &
                                              dt_therm, grid%nghost, ice%ncat, ice%nk_ice, &
                                              ms%nz_ml, grid%nx_total, grid%ny_total)
      end if
   end subroutine ice_frazil_uptake

   pure subroutine ice_frazil_uptake_impl(hTr_T, hTr_S, h_layer, wet_mask, eos, &
                                          frazil_heat, m_ice, m_snow, enth_ice, sal_ice, &
                                          m_frozen_diag, salt_flux_diag, &
                                          dt_therm, nghost, ncat, nk, nz, nx, ny)
      !! Device kernel over PHYSICAL cells (ghosts excluded — same
      !! physical-cells-only contract as `ice_frazil_accumulate_impl` and
      !! `ice_thermo_columns`). Per wet, non-vanished, banked cell:
      !! sample SST/SSS at `k = nz`, compute the seawater freezing point,
      !! spend the WHOLE bank on category-1's column via
      !! `ice_frazil_uptake_column`, and reset `frazil_heat` to 0 (fully
      !! spent). The per-cell diags (`m_frozen_diag`, `salt_flux_diag`)
      !! are zeroed UNCONDITIONALLY first, then overwritten under the
      !! gate — an unbanked/dry/land/vanished cell reports zero, not a
      !! stale value from a prior window.
      !!
      !! Decl-order: all integer dims declared before the explicit-shape
      !! arrays that use them. Inner `if` gate only (never a masked
      !! `do concurrent` header, per the dc-to-omp constraint).
      !!
      !! Bank on dry/vanished columns is NOT zeroed here — it stays
      !! banked (conserved) until the column is wet+intact enough to
      !! spend it, matching the "never discard a bank" contract of
      !! `rdb_ice_frazil`.
      integer, intent(in) :: nghost, ncat, nk, nz, nx, ny
      real(wp), intent(in) :: hTr_T(nx, ny, nz)
      real(wp), intent(in) :: hTr_S(nx, ny, nz)
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      real(wp), intent(in) :: wet_mask(nx, ny)
      type(eos_t), intent(in) :: eos
      real(wp), intent(inout) :: frazil_heat(nx, ny)
      real(wp), intent(inout) :: m_ice(nx, ny, ncat)
      real(wp), intent(inout) :: m_snow(nx, ny, ncat)
      real(wp), intent(inout) :: enth_ice(nx, ny, ncat, nk)
      real(wp), intent(inout) :: sal_ice(nx, ny, ncat, nk)
      real(wp), intent(inout) :: m_frozen_diag(nx, ny)
      real(wp), intent(inout) :: salt_flux_diag(nx, ny)
      real(wp), intent(in) :: dt_therm

      integer :: i, j, i_lo, i_hi, j_lo, j_hi
      real(wp) :: h, sst, s_surf, tfw
      real(wp) :: m_snow_pt, m_ice_pt, m_frozen_pt, salt_to_ice_pt
      real(wp) :: enth_ice_col(ICE_NK_MAX), sal_ice_col(ICE_NK_MAX)

      i_lo = nghost + 1
      i_hi = nx - nghost
      j_lo = nghost + 1
      j_hi = ny - nghost

      do concurrent(j=j_lo:j_hi, i=i_lo:i_hi) local(h, sst, s_surf, tfw, &
                                                    m_snow_pt, m_ice_pt, m_frozen_pt, &
                                                    salt_to_ice_pt, enth_ice_col, sal_ice_col)
         m_frozen_diag(i, j) = 0.0_wp
         salt_flux_diag(i, j) = 0.0_wp
         h = h_layer(i, j, nz)
         if (wet_mask(i, j) > 0.5_wp .and. h > H_VANISHED .and. frazil_heat(i, j) > 0.0_wp) then
            sst = hTr_T(i, j, nz)/h
            s_surf = hTr_S(i, j, nz)/h
            tfw = eos_freezing_point(eos, s_surf, 0.0_wp)

            m_snow_pt = m_snow(i, j, 1)
            m_ice_pt = m_ice(i, j, 1)
            enth_ice_col(1:nk) = enth_ice(i, j, 1, 1:nk)
            sal_ice_col(1:nk) = sal_ice(i, j, 1, 1:nk)

            call ice_frazil_uptake_column(nk, frazil_heat(i, j), tfw, sst, s_surf, &
                                          m_snow_pt, m_ice_pt, enth_ice_col(1:nk), &
                                          sal_ice_col(1:nk), m_frozen_pt, salt_to_ice_pt)

            m_ice(i, j, 1) = m_ice_pt
            enth_ice(i, j, 1, 1:nk) = enth_ice_col(1:nk)
            sal_ice(i, j, 1, 1:nk) = sal_ice_col(1:nk)
            m_frozen_diag(i, j) = m_frozen_pt
            salt_flux_diag(i, j) = (m_frozen_pt*s_surf - salt_to_ice_pt)/dt_therm
            frazil_heat(i, j) = 0.0_wp
         end if
      end do
   end subroutine ice_frazil_uptake_impl

   pure subroutine ice_frazil_uptake_multicat_impl(hTr_T, hTr_S, h_layer, wet_mask, eos, &
                                                   frazil_heat, part_size, m_ice, m_snow, &
                                                   enth_ice, sal_ice, &
                                                   m_frozen_diag, salt_flux_diag, &
                                                   dt_therm, nghost, ncat, nk, nz, nx, ny)
      !! ncat>1 SIS2 ITD-mode frazil spend. Port of SIS2
      !! `SIS_slow_thermo.F90:1121-1146 + 1181-1186` (non-filling mode —
      !! SIS2's default `SIS2_FILLING_FRAZIL=.true.` thin-category fill is
      !! DEFERRED, see module docstring). Per banked cell (same
      !! wet/non-vanished/bank>0 gate as `ice_frazil_uptake_impl`):
      !!   1. k_merge scan (SIS2:1124-1129): first category `c` with
      !!      `part(0) + part(c) > 0.01`; falls back to `k_merge = 1` if
      !!      no category qualifies (SIS2's `k_merge` default).
      !!   2. Open-water annexation (SIS2:1131-1145): if `part(0) > 0`,
      !!      dilute category `k_merge`'s thickness at CONSTANT MASS —
      !!      `m_ice(k_merge)`/`m_snow(k_merge)` scale by
      !!      `part(k_merge)/(part(k_merge)+part(0))`, `part(k_merge)`
      !!      absorbs all of `part(0)`, `part(0)` -> 0. `enth`/`sal` are
      !!      per-MASS intensive — untouched by an area-only dilution.
      !!   3. Per-ice-area spend (SIS2:1181-1186):
      !!      `frazil_col = frazil_heat/part(k_merge)` (J per m² of
      !!      category area; the denominator is > 0 by construction —
      !!      either an occupied category was found, or step 2 just grew
      !!      `part(k_merge)` from the Σpart=1 invariant), then the SAME
      !!      UNCHANGED `ice_frazil_uptake_column` spends it on category
      !!      `k_merge`'s column.
      !!   4. Diags (per cell, part-weighted back to CELL-area units to
      !!      match the ncat==1 diag convention that the couplers and
      !!      `rdb_ice_thermo_driver` consume): `m_frozen_diag =
      !!      part(k_merge)*m_frozen_pt`; `salt_flux_diag =
      !!      part(k_merge)*(m_frozen_pt*s_surf - salt_to_ice_pt)/dt_therm`;
      !!      `frazil_heat` reset to 0 (fully spent). Both diags are
      !!      zeroed UNCONDITIONALLY at loop top, same contract as
      !!      `ice_frazil_uptake_impl`.
      !!
      !! Decl-order: all integer dims declared before the explicit-shape
      !! arrays that use them. Inner `if` gate only, serial `do cat` scan
      !! for k_merge — never a masked `do concurrent` header.
      integer, intent(in) :: nghost, ncat, nk, nz, nx, ny
      real(wp), intent(in) :: hTr_T(nx, ny, nz)
      real(wp), intent(in) :: hTr_S(nx, ny, nz)
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      real(wp), intent(in) :: wet_mask(nx, ny)
      type(eos_t), intent(in) :: eos
      real(wp), intent(inout) :: frazil_heat(nx, ny)
      real(wp), intent(inout) :: part_size(nx, ny, 0:ncat)
      real(wp), intent(inout) :: m_ice(nx, ny, ncat)
      real(wp), intent(inout) :: m_snow(nx, ny, ncat)
      real(wp), intent(inout) :: enth_ice(nx, ny, ncat, nk)
      real(wp), intent(inout) :: sal_ice(nx, ny, ncat, nk)
      real(wp), intent(inout) :: m_frozen_diag(nx, ny)
      real(wp), intent(inout) :: salt_flux_diag(nx, ny)
      real(wp), intent(in) :: dt_therm

      integer :: i, j, i_lo, i_hi, j_lo, j_hi
      integer :: cat, k_merge
      real(wp) :: h, sst, s_surf, tfw, i_part
      real(wp) :: m_snow_pt, m_ice_pt, m_frozen_pt, salt_to_ice_pt, frazil_col
      real(wp) :: enth_ice_col(ICE_NK_MAX), sal_ice_col(ICE_NK_MAX)

      i_lo = nghost + 1
      i_hi = nx - nghost
      j_lo = nghost + 1
      j_hi = ny - nghost

      do concurrent(j=j_lo:j_hi, i=i_lo:i_hi) local(h, sst, s_surf, tfw, i_part, &
                                                    cat, k_merge, &
                                                    m_snow_pt, m_ice_pt, m_frozen_pt, &
                                                    salt_to_ice_pt, frazil_col, &
                                                    enth_ice_col, sal_ice_col)
         m_frozen_diag(i, j) = 0.0_wp
         salt_flux_diag(i, j) = 0.0_wp
         h = h_layer(i, j, nz)
         if (wet_mask(i, j) > 0.5_wp .and. h > H_VANISHED .and. frazil_heat(i, j) > 0.0_wp) then
            sst = hTr_T(i, j, nz)/h
            s_surf = hTr_S(i, j, nz)/h
            tfw = eos_freezing_point(eos, s_surf, 0.0_wp)

            ! ---- 1. k_merge scan (SIS2:1124-1129) ----
            k_merge = 1
            do cat = 1, ncat
               if (part_size(i, j, 0) + part_size(i, j, cat) > 0.01_wp) then
                  k_merge = cat
                  exit
               end if
            end do

            ! ---- 2. Open-water annexation (SIS2:1131-1145) ----
            if (part_size(i, j, 0) > 0.0_wp) then
               i_part = 1.0_wp/(part_size(i, j, k_merge) + part_size(i, j, 0))
               m_snow(i, j, k_merge) = (m_snow(i, j, k_merge)*part_size(i, j, k_merge))*i_part
               m_ice(i, j, k_merge) = (m_ice(i, j, k_merge)*part_size(i, j, k_merge))*i_part
               part_size(i, j, k_merge) = part_size(i, j, k_merge) + part_size(i, j, 0)
               part_size(i, j, 0) = 0.0_wp
            end if

            ! ---- 3. Per-ice-area spend (SIS2:1181-1186) ----
            frazil_col = frazil_heat(i, j)/part_size(i, j, k_merge)

            m_snow_pt = m_snow(i, j, k_merge)
            m_ice_pt = m_ice(i, j, k_merge)
            enth_ice_col(1:nk) = enth_ice(i, j, k_merge, 1:nk)
            sal_ice_col(1:nk) = sal_ice(i, j, k_merge, 1:nk)

            call ice_frazil_uptake_column(nk, frazil_col, tfw, sst, s_surf, &
                                          m_snow_pt, m_ice_pt, enth_ice_col(1:nk), &
                                          sal_ice_col(1:nk), m_frozen_pt, salt_to_ice_pt)

            m_ice(i, j, k_merge) = m_ice_pt
            enth_ice(i, j, k_merge, 1:nk) = enth_ice_col(1:nk)
            sal_ice(i, j, k_merge, 1:nk) = sal_ice_col(1:nk)

            ! ---- 4. Diags (per cell) ----
            m_frozen_diag(i, j) = part_size(i, j, k_merge)*m_frozen_pt
            salt_flux_diag(i, j) = part_size(i, j, k_merge) &
                                   *(m_frozen_pt*s_surf - salt_to_ice_pt)/dt_therm
            frazil_heat(i, j) = 0.0_wp
         end if
      end do
   end subroutine ice_frazil_uptake_multicat_impl

   pure subroutine ice_frazil_uptake_column(nk, frazil, tfw, sst, s_surf, m_snow, &
                                            m_ice_tot, enth_ice_bu, sal_ice_bu, &
                                            m_frozen, salt_to_ice)
      !! Per-(cell,category-1) column worker: gather + flip
      !! (bottom-up -> top-down, TRAP #2 discipline), spend the whole
      !! `frazil` bank evenly over the `nk` layers (SIS2:1402), rebalance,
      !! scatter + flip back. Port of `add_frazil_SIS2`
      !! (SIS2_ice_thm.F90:1342-1444), bulk-salinity mode
      !! (`salin_freeze = ICE_BULK_SALINITY`, SIS_slow_thermo.F90:
      !! 1207-1209) — the `ice_rel_salin` mode is NOT ported.
      !!
      !! Snow slot (`m_snow`, index 0 in the local top-down column) is
      !! carried through untouched — the frazil spend loop only ever
      !! touches `k = 1..nk` (ice layers); `enthalpy(0)`/`salin(0)` are
      !! placeholder zeros that `ice_rebalance_layers` neither reads nor
      !! writes.
      !!
      !! Empty-column degeneracy: `m_lay(k) = 0` on entry makes the
      !! mass-weighted mix reduce exactly to
      !! `enthalpy(k) = enth_frazil`, `salin(k) = ICE_BULK_SALINITY` —
      !! SIS2's post-rebalance `mH_ice == 0 => enthalpy_liquid_freeze`
      !! reset (SIS2_ice_thm.F90:1420-1424) is unreachable here (a spend
      !! only runs when `frazil > 0`, which always deposits mass) and is
      !! deliberately NOT ported.
      !$acc routine seq
      integer, intent(in) :: nk
         !! Number of ice layers (declared first — decl-order).
      real(wp), intent(in) :: frazil
         !! Whole per-cell frazil bank to spend this call (J/m² of cell).
      real(wp), intent(in) :: tfw
         !! Seawater freezing temperature at the surface (degC).
      real(wp), intent(in) :: sst
         !! Sea-surface temperature (degC) — feeds the TRAP-#1 liquid
         !! ocean enthalpy.
      real(wp), intent(in) :: s_surf
         !! Sea-surface salinity (PSU) — feeds the TRAP-#1 liquid ocean
         !! enthalpy (unused by the linear formula, kept for call-site
         !! parity with `ice_enthalpy_liquid`).
      real(wp), intent(in) :: m_snow
         !! Snow mass per unit CELL area (kg/m²) — untouched, carried
         !! through only to seed the local column's slot 0.
      real(wp), intent(inout) :: m_ice_tot
         !! Total category-1 ice mass per unit CELL area (kg/m²); in =
         !! prior step, out = post-freeze.
      real(wp), intent(inout) :: enth_ice_bu(nk)
         !! BOTTOM-UP ice specific enthalpies (J/kg) — state order,
         !! `enth_ice_bu(1)` = ice bottom.
      real(wp), intent(inout) :: sal_ice_bu(nk)
         !! BOTTOM-UP ice bulk salinities (PSU) — state order.
      real(wp), intent(out) :: m_frozen
         !! Total new-ice mass formed this call (kg/m² of cell).
      real(wp), intent(out) :: salt_to_ice
         !! Salt content retained by the new ice (kg/m² of cell) —
         !! `m_frozen * ICE_BULK_SALINITY`.

      real(wp) :: m_lay(0:ICE_NK_MAX), enthalpy(0:ICE_NK_MAX), salin(0:ICE_NK_MAX)
      real(wp) :: frazil_per_layer, enth_ocean, min_denth, t_frazil, enth_frazil, m_frazil
      real(wp) :: mtot_ice
      integer :: k

      ! ---- Gather + flip (TRAP #2): bottom-up state -> top-down local ----
      do k = 1, nk
         enthalpy(k) = enth_ice_bu(nk + 1 - k)
         salin(k) = sal_ice_bu(nk + 1 - k)
      end do
      m_lay(0) = m_snow
      enthalpy(0) = 0.0_wp
      salin(0) = 0.0_wp
      do k = 1, nk
         m_lay(k) = m_ice_tot/real(nk, wp)
      end do

      ! ---- Spend the whole bank, evenly split, top-down (SIS2:1402) ----
      frazil_per_layer = frazil/real(nk, wp)
      enth_ocean = ice_enthalpy_liquid(sst, s_surf)   ! TRAP #1
      min_denth = ICE_LAT_FUS*(1.0_wp - ICE_LIQ_LIM)
      m_frozen = 0.0_wp
      salt_to_ice = 0.0_wp
      do k = 1, nk
         t_frazil = min(tfw, ice_t_freeze(salin(k)) - ICE_FRAZIL_T_OFFSET)
         enth_frazil = min(ice_enth_from_ts(t_frazil, salin(k)), enth_ocean - min_denth)
         m_frazil = frazil_per_layer/(enth_ocean - enth_frazil)
         enthalpy(k) = (m_lay(k)*enthalpy(k) + m_frazil*enth_frazil)/(m_lay(k) + m_frazil)
         salin(k) = (m_lay(k)*salin(k) + m_frazil*ICE_BULK_SALINITY)/(m_lay(k) + m_frazil)
         m_lay(k) = m_lay(k) + m_frazil
         m_frozen = m_frozen + m_frazil
         salt_to_ice = salt_to_ice + m_frazil*ICE_BULK_SALINITY
      end do

      call ice_rebalance_layers(nk, m_lay(0:nk), enthalpy(0:nk), salin(0:nk), mtot_ice)

      ! ---- Scatter + flip back: top-down local -> bottom-up state ----
      m_ice_tot = mtot_ice
      do k = 1, nk
         enth_ice_bu(nk + 1 - k) = enthalpy(k)
         sal_ice_bu(nk + 1 - k) = salin(k)
      end do
   end subroutine ice_frazil_uptake_column

end module rdb_ice_frazil_uptake
