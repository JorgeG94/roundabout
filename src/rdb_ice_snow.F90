!! Ice-free share of the snowfall source term, delivered to the ocean
!! (SIS2 port, PR 26).
module rdb_ice_snow
   !! Complements `rdb_ice_mass%ice_snow_accumulate` (the on-ice snow add):
   !! the share of `&ocean_ice_nml snowfall` that lands where there is no
   !! ice (open water at `ncat>1`, ice-free cells at `ncat==1`, and any
   !! category that fails the column's own entry gate,
   !! `m_ice(i,j,c) > ICE_RHO_ICE*H_VANISHED`) does not accumulate — it
   !! melts on contact with seawater. This module delivers that mass's
   !! latent heat + virtual freshening to the ocean through the EXISTING
   !! `heat_flux_diag`/`salt_flux_diag` contributor seam
   !! (`rdb_ice_thermo_driver.F90`'s ordering contract), the same one
   !! `ice_frazil_uptake` and the melt-side reduce kernels already use.
   !! Without this, `snowfall > 0` over open water would silently
   !! annihilate mass and energy (the `sw_thru` failure mode,
   !! PLAN_PR26_snowfall.md §2 / roadmap §5 trap 15).
   !!
   !! **Physics** (PLAN_PR26_snowfall.md §3). Per unit CELL area, the mass
   !! delivered to the ocean over the thermo window is
   !!   m_ocn = snow_part_ocn * atm_fprec * dt_therm         [kg/m^2]
   !! It arrives as fresh solid at 0 degC, specific enthalpy
   !!   enth_fall = ice_enth_from_ts(0, 0) = ICE_ENTH_LIQ_0 - ICE_LAT_FUS
   !!             = -3.34e5 J/kg (SIS2's "-LI" convention for frozen
   !!               precip entering the ocean, SIS_sum_output.F90:783),
   !! module constant `ENTH_SNOWFALL` below. Bringing it to the ocean's
   !! own liquid enthalpy `enth_ocean = ice_enthalpy_liquid(sst, s_surf)`
   !! costs the ocean `dE = m_ocn*(enth_ocean - ENTH_SNOWFALL)` [J/m^2],
   !! so in the `heat_flux_diag` sign convention (positive DOWN into the
   !! ocean, `rdb_ice_thermo_driver.F90:14`):
   !!   heat_flux_diag -= m_ocn*(enth_ocean - ENTH_SNOWFALL)/dt_therm
   !! and, in the ice model's existing virtual-salt convention
   !! (`salt_flux_diag += m_net*(s_surf - S_ice)/dt_therm`, `m_net`
   !! positive = freeze), fresh water added to the ocean is
   !! `m_net = -m_ocn` with `S_snow = 0`:
   !!   salt_flux_diag -= m_ocn*ssurf_seam/dt_therm
   !! The ocean's water MASS is not increased -- Roundabout's ocean is a
   !! volume-conserving virtual-salt-flux model today, the same
   !! convention `rdb_ice_frazil_uptake`'s docstring states for the
   !! freeze side. PR-16 converts this (and the melt/frazil paths) to
   !! real mass; see the `TODO(PR-16)` marker below and this module's
   !! binding seam spec.
   !!
   !! **Contributor ordering (mandatory, driver-enforced)**: this module
   !! MUST run AFTER `ice_thermo_driver_step` (which unconditionally
   !! zeroes `heat_flux_diag` and writes the melt-side `salt_flux_diag`
   !! contribution) and BEFORE `ice_ocean_brine_flux`/`ice_ocean_heat_flux`
   !! (which overwrite `Q_salt`/`Q_heat` from `salt_flux_diag`/
   !! `heat_flux_diag`) -- see `rdb_driver.F90`'s mandated-order comment
   !! block. It ADDS to both diags, never zeroes them (third contributor
   !! on `heat_flux_diag`, after the column; second on `salt_flux_diag`,
   !! after the melt-side reduce).
   !!
   !! **Binding seam spec (PLAN_PR26_snowfall.md §13, owner: PR-16).**
   !! `ocean_sea_ice_t%fprec_ocn_diag(nx_total, ny_total)`: real(wp),
   !! kg/m^2/s, per unit CELL area, positive = frozen fresh water
   !! delivered to the ocean surface, zero on land and fully ice-covered
   !! cells. SCRATCH: zeroed + rewritten every thermo window by this
   !! module; not restart-carried; device-mapped `copyin` by
   !! `ocean_sea_ice_enter_data_impl`. When PR-16 makes ice<->ocean mass
   !! real, it consumes `fprec_ocn_diag` as a `net_massin` source at
   !! `h_layer(:,:,nz)` and MUST remove this module's virtual-salt term
   !! (marked `TODO(PR-16)` below) or the dilution double-counts; the
   !! heat term stays (becomes the `heat_content_fprec` companion, not a
   !! duplicate).
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ice_enthalpy, only: ICE_ENTH_LIQ_0, ICE_LAT_FUS, ice_enthalpy_liquid
   use rdb_ice_state, only: ocean_sea_ice_t
   implicit none
   private

   public :: ice_snowfall_ocean_share

   real(wp), parameter :: ENTH_SNOWFALL = ICE_ENTH_LIQ_0 - ICE_LAT_FUS
      !! Specific enthalpy of frozen precipitation entering the ocean
      !! (J/kg) -- `ice_enth_from_ts(0, 0)` evaluated in closed form
      !! (fresh water, T<=0: `(ICE_ENTH_LIQ_0 - ICE_LAT_FUS) + ICE_CP_ICE*0`).
      !! SIS2's "-LI" convention (SIS_sum_output.F90:783).

contains

   pure subroutine ice_snowfall_ocean_share(grid, ice, dt_therm)
      !! Outer shim (outer-shim + flat-impl pattern): forward the ice
      !! slot's arrays to the device kernel. Called only when
      !! `ice%has_snowfall` (driver gate, `rdb_driver.F90`) -- default
      !! `snowfall=0` never reaches this module.
      type(hgrid_t), intent(in) :: grid
      type(ocean_sea_ice_t), intent(inout) :: ice
      real(wp), intent(in) :: dt_therm
         !! Effective thermo timestep (s) -- `ocean_dyn_t%therm_dt(dt)`.

      call ice_snowfall_ocean_share_impl(ice%snow_part_ocn, ice%atm_fprec, &
                                         ice%sst_seam, ice%ssurf_seam, &
                                         ice%fprec_ocn_diag, ice%heat_flux_diag, &
                                         ice%salt_flux_diag, dt_therm, grid%nghost, &
                                         ice%nx_total, ice%ny_total)
   end subroutine ice_snowfall_ocean_share

   pure subroutine ice_snowfall_ocean_share_impl(snow_part_ocn, atm_fprec, sst_seam, &
                                                 ssurf_seam, fprec_ocn_diag, &
                                                 heat_flux_diag, salt_flux_diag, &
                                                 dt_therm, nghost, nx, ny)
      !! Device kernel over PHYSICAL cells (ghosts excluded). Gated on
      !! `snow_part_ocn(i,j) > 0` rather than a separate `wet_mask` arg:
      !! `ice_snow_part_ocn_fill_impl` (`rdb_ice_thermo_driver`) already
      !! zeroes `snow_part_ocn` on land/dry/fully-ice-covered cells, so
      !! that field IS the wet-and-ice-free gate this kernel needs --
      !! functionally identical to (and cheaper than) re-deriving a
      !! `wet_mask > 0.5` test here, and keeps this module's signature to
      !! the `(grid, ice, dt_therm)` seam (no `ms` dependency).
      !!
      !! `fprec_ocn_diag` is zeroed unconditionally at loop top (this
      !! module owns it outright, same contract as `heat_flux_diag`/
      !! `m_melt_diag` in `ice_thermo_driver_reduce_impl`).
      !! `heat_flux_diag`/`salt_flux_diag` are ADDED to, never zeroed --
      !! the ordering contract this module's docstring states, now the
      !! second (`salt_flux_diag`) / third (`heat_flux_diag`) contributor.
      !!
      !! Reuses the `sst_seam`/`ssurf_seam` sample `ice_compute_basal_
      !! flux` filled earlier in the same thermo window (module
      !! docstring) -- does not re-derive from `ms`/`eos`.
      !!
      !! Decl-order: all integer dims declared before the explicit-shape
      !! arrays that use them. Inner `if` gate only (never a masked
      !! `do concurrent` header).
      integer, intent(in) :: nghost, nx, ny
      real(wp), intent(in) :: snow_part_ocn(nx, ny)
      real(wp), intent(in) :: atm_fprec(nx, ny)
      real(wp), intent(in) :: sst_seam(nx, ny)
      real(wp), intent(in) :: ssurf_seam(nx, ny)
      real(wp), intent(inout) :: fprec_ocn_diag(nx, ny)
      real(wp), intent(inout) :: heat_flux_diag(nx, ny)
      real(wp), intent(inout) :: salt_flux_diag(nx, ny)
      real(wp), intent(in) :: dt_therm

      integer :: i, j, i_lo, i_hi, j_lo, j_hi
      real(wp) :: m_ocn, enth_ocean

      i_lo = nghost + 1
      i_hi = nx - nghost
      j_lo = nghost + 1
      j_hi = ny - nghost

      do concurrent(j=j_lo:j_hi, i=i_lo:i_hi) local(m_ocn, enth_ocean)
         fprec_ocn_diag(i, j) = 0.0_wp
         if (snow_part_ocn(i, j) > 0.0_wp) then
            fprec_ocn_diag(i, j) = snow_part_ocn(i, j)*atm_fprec(i, j)
            m_ocn = fprec_ocn_diag(i, j)*dt_therm
            enth_ocean = ice_enthalpy_liquid(sst_seam(i, j), ssurf_seam(i, j))
            heat_flux_diag(i, j) = heat_flux_diag(i, j) &
                                   - m_ocn*(enth_ocean - ENTH_SNOWFALL)/dt_therm
            ! TODO(PR-16): virtual-salt term; delete when net_massin lands
            ! (fprec_ocn_diag becomes a real freshwater source instead).
            salt_flux_diag(i, j) = salt_flux_diag(i, j) - m_ocn*ssurf_seam(i, j)/dt_therm
         end if
      end do
   end subroutine ice_snowfall_ocean_share_impl

end module rdb_ice_snow
