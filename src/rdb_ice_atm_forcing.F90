!! v1 slab-atmosphere restoring filler for the ice atmospheric-forcing seam
!! (SIS2 port, PR 3c; PR 26 adds the fourth field, `atm_fprec`) + the
!! JRA55-do-driven `ice_atm_forcing_from_atm` seam (OM3 forcing wave-3,
!! C2c-3b).
module rdb_ice_atm_forcing
   !! Fills ocean_sea_ice_t's coupleable seam
   !! (atm_sf0/atm_dsfdt/atm_sw_dn/atm_fprec) either from the
   !! &ocean_ice_nml scalar restoring knobs (`ice_atm_forcing_restoring`)
   !! or from the prescribed atmospheric state (`ice_atm_forcing_from_atm`,
   !! C2c-3b) -- the engine dispatches on `atm_state_t%enable`. SF(T) is
   !! always upward-positive, linearised at a surface temperature `T`:
   !! SF(T) = sf_0 + dsf_dt*T. The seam is identical whether this filler,
   !! the atm-driven one, or a live coupler fills it -- that is what makes
   !! the model coupleable (PLAN_ICE_PR3c S"coupleable seam";
   !! PLAN_PR26_snowfall.md S13 fixes atm_fprec's shape/units/lifecycle
   !! for PR-55 to fill from data; PLAN_C2c.md S4 "Ice seam").
   !!
   !! **`ice_atm_forcing_restoring` (v1, unchanged).** SCALAR/UNIFORM:
   !! `SF(T) = lambda*(T - T_air) => sf_0 = -lambda*T_air`, `dsf_dt =
   !! lambda`, `sw_dn = sw_down` (uniform), `atm_fprec = snowfall`
   !! (uniform). Ghost rows are filled too (harmless -- ice_thermo_columns
   !! only steps physical cells). Default-off byte-identity: never called
   !! unless ice%enable.
   !!
   !! **`ice_atm_forcing_from_atm` (C2c-3b).** A linearisation of the
   !! surface energy balance at `T_s = 0 degC` (`BULK_T0_KELVIN`) using
   !! Large & Yeager's (2004, 2009) over-ice constant transfer
   !! coefficients (`C_D = C_H = C_E = 1.63e-3` -- no Monin-Obukhov
   !! stability iteration, unlike the ocean's `bulk_flux_column`; this is
   !! the simpler CORE/OMIP "sea-ice" bulk formula, not the open-water
   !! one): with `rho` the ambient air density (ideal-gas law off `psl`/
   !! `tas`/`huss`, same form `bulk_flux_column` uses) and `wv =
   !! max(|(uas,vas)|, BULK_WIND_MIN)` the atmospheric wind speed (no ice
   !! velocity subtracted -- same "no ice-specific bulk drag law" posture
   !! as `tau_a_x/y`, D7 in `rdb_ice_state`):
   !!
   !!   sf_0   = eps*sigma*T0^4 - rlds
   !!            + rho*cp*C_H*wv*(T0 - tas)
   !!            + rho*L_s*C_E*wv*(qsat_ice(T0) - huss)
   !!   dsf_dt = 4*eps*sigma*T0^3
   !!            + rho*cp*C_H*wv
   !!            + rho*L_s*C_E*wv*dqsat_ice_dT(T0)
   !!
   !! (`eps` is the atm-forcing bulk-flux emissivity knob,
   !! `&ocean_bulk_flux_nml emissivity` -- LY09's over-ice formula does
   !! not distinguish an ocean- from an ice-emissivity, so this reuses the
   !! one knob rather than adding a second; `qsat_ice` is
   !! `bulk_flux_qsat_ice`, `L_s` is `BULK_LATENT_SUBLIMATION`).
   !! `atm_sw_dn = rsds`, `atm_fprec = prsn` -- same fields the ocean
   !! consumer reads, same "ghosts untouched" contract as `atm_state_t`
   !! itself (its own docstring: no consumer reads a ghost cell; a ghost
   !! `tas = 0` would divide `rho` by zero). Unlike the restoring filler,
   !! this kernel therefore loops the PHYSICAL interior only (the
   !! `atm%nghost` offset idiom `atm_state_scatter_impl`/
   !! `bulk_forcing_core_impl` already use) -- ghost cells of the four
   !! seam arrays simply keep whatever the restoring filler (if ever
   !! called earlier) or the `source=0.0_wp` allocation left there, which
   !! is harmless by the same "ice_thermo_columns only steps physical
   !! cells" contract.
   use rdb_constants, only: wp
   use rdb_ice_state, only: ocean_sea_ice_t
#ifndef RDB_NO_NETCDF
   use rdb_atm_state, only: atm_state_t
#endif
   use rdb_ocean_bulk_flux, only: bulk_flux_config_t, bulk_flux_qsat_ice, &
                                  BULK_WIND_MIN, BULK_T0_KELVIN, BULK_STEFAN_BOLTZMANN, &
                                  BULK_GAS_CONST_AIR, BULK_TVQ_AIR, BULK_CP_AIR, &
                                  BULK_LATENT_SUBLIMATION, BULK_QSAT_ICE_COEF_B
   implicit none
   private
   public :: ice_atm_forcing_restoring
#ifndef RDB_NO_NETCDF
   public :: ice_atm_forcing_from_atm
#endif
contains
   pure subroutine ice_atm_forcing_restoring(ice, air_temp, restore_lambda, sw_down, snowfall)
      !! Outer shim: forward the four seam arrays + four scalars to the
      !! device kernel. No registry indirection here, but keep the shim+_impl
      !! split for the explicit-shape device-kernel discipline.
      type(ocean_sea_ice_t), intent(inout) :: ice
      real(wp), intent(in) :: air_temp, restore_lambda, sw_down, snowfall
      call ice_atm_forcing_restoring_impl(ice%atm_sf0, ice%atm_dsfdt, ice%atm_sw_dn, &
                                          ice%atm_fprec, air_temp, restore_lambda, &
                                          sw_down, snowfall, ice%nx_total, ice%ny_total)
   end subroutine ice_atm_forcing_restoring

   pure subroutine ice_atm_forcing_restoring_impl(sf0, dsfdt, sw, fpr, air_temp, &
                                                  restore_lambda, sw_down, snowfall, nx, ny)
      !! Full-array fill (incl. ghosts). Explicit-shape + decl-order (dims first).
      !! Runs on device-resident seam arrays (mapped by ice%enter_data).
      integer, intent(in) :: nx, ny
      real(wp), intent(inout) :: sf0(nx, ny), dsfdt(nx, ny), sw(nx, ny), fpr(nx, ny)
      real(wp), intent(in) :: air_temp, restore_lambda, sw_down, snowfall
      integer :: i, j
      do concurrent(j=1:ny, i=1:nx)
         sf0(i, j) = -restore_lambda*air_temp
         dsfdt(i, j) = restore_lambda
         sw(i, j) = sw_down
         fpr(i, j) = snowfall
      end do
   end subroutine ice_atm_forcing_restoring_impl

#ifndef RDB_NO_NETCDF
   pure subroutine ice_atm_forcing_from_atm(ice, atm, cfg)
      !! Outer shim: forward the four seam arrays, the nine-field atm
      !! state's own cached geometry (`atm%nx_phys/ny_phys/nghost/
      !! nx_total/ny_total`) and the bulk-flux emissivity knob to the
      !! device kernel. A no-op unless `atm%enable` (mirrors every other
      !! atm-provider consumer's own guard, even though the engine
      !! dispatcher only calls this when it already holds).
      type(ocean_sea_ice_t), intent(inout) :: ice
      type(atm_state_t), intent(in) :: atm
      type(bulk_flux_config_t), intent(in) :: cfg
      if (.not. atm%enable) return
      call ice_atm_forcing_from_atm_impl(ice%atm_sf0, ice%atm_dsfdt, ice%atm_sw_dn, &
                                         ice%atm_fprec, &
                                         atm%uas, atm%vas, atm%tas, atm%huss, atm%psl, &
                                         atm%rlds, atm%rsds, atm%prsn, cfg%emissivity, &
                                         atm%nx_total, atm%ny_total, atm%nghost, &
                                         atm%nx_phys, atm%ny_phys)
   end subroutine ice_atm_forcing_from_atm
#endif

   pure subroutine ice_atm_forcing_from_atm_impl(sf0, dsfdt, sw, fpr, uas, vas, tas, huss, &
                                                 psl, rlds, rsds, prsn, emissivity, &
                                                 nx, ny, ng, nxp, nyp)
      !! Physical-interior fill (module docstring: ghosts untouched here,
      !! unlike the restoring filler -- the atm fields are zero in the
      !! ghost ring and would divide `rho_air` by zero). Explicit-shape +
      !! decl-order (dims first).
      integer, intent(in) :: nx, ny, ng, nxp, nyp
      real(wp), intent(inout) :: sf0(nx, ny), dsfdt(nx, ny), sw(nx, ny), fpr(nx, ny)
      real(wp), intent(in) :: uas(nx, ny), vas(nx, ny), tas(nx, ny), huss(nx, ny), psl(nx, ny)
      real(wp), intent(in) :: rlds(nx, ny), rsds(nx, ny), prsn(nx, ny)
      real(wp), intent(in) :: emissivity
      integer :: i, j, ii, jj
      real(wp) :: rho_air, wv, qsat0, dqsat0_dt, lw_t3, sens_coef, lat_coef
      real(wp), parameter :: ICE_BULK_CH = 1.63e-3_wp
         !! Large & Yeager (2004, 2009) constant over-ice Stanton number
         !! -- no stability correction (contrast the ocean path's
         !! Monin-Obukhov-iterated `ch`).
      real(wp), parameter :: ICE_BULK_CE = 1.63e-3_wp
         !! Ditto, the constant over-ice Dalton number. LY09 sets
         !! `C_D = C_H = C_E = 1.63e-3` over ice; `C_D` itself is not
         !! needed here (no atm-ice stress term in this seam -- D7 in
         !! `rdb_ice_state`, a separate deferred item).

      do concurrent(j=1:nyp, i=1:nxp) &
         local(ii, jj, rho_air, wv, qsat0, dqsat0_dt, lw_t3, sens_coef, lat_coef)
         ii = ng + i
         jj = ng + j
         rho_air = psl(ii, jj)/(BULK_GAS_CONST_AIR*tas(ii, jj)*(1.0_wp + BULK_TVQ_AIR*huss(ii, jj)))
         wv = max(sqrt(uas(ii, jj)**2 + vas(ii, jj)**2), BULK_WIND_MIN)
         qsat0 = bulk_flux_qsat_ice(BULK_T0_KELVIN, rho_air)
         dqsat0_dt = qsat0*BULK_QSAT_ICE_COEF_B/(BULK_T0_KELVIN*BULK_T0_KELVIN)
         ! lw_t3 = eps*sigma*T0**3, reused for BOTH the sf_0 emission term
         ! (lw_t3*T0 = eps*sigma*T0**4) and the dsf_dt slope (4*lw_t3 =
         ! 4*eps*sigma*T0**3) -- one evaluation of the cubic, not two.
         lw_t3 = emissivity*BULK_STEFAN_BOLTZMANN*BULK_T0_KELVIN**3
         sens_coef = rho_air*BULK_CP_AIR*ICE_BULK_CH*wv
         lat_coef = rho_air*BULK_LATENT_SUBLIMATION*ICE_BULK_CE*wv

         sf0(ii, jj) = lw_t3*BULK_T0_KELVIN - rlds(ii, jj) + &
                       sens_coef*(BULK_T0_KELVIN - tas(ii, jj)) + &
                       lat_coef*(qsat0 - huss(ii, jj))
         dsfdt(ii, jj) = 4.0_wp*lw_t3 + sens_coef + lat_coef*dqsat0_dt
         sw(ii, jj) = rsds(ii, jj)
         fpr(ii, jj) = prsn(ii, jj)
      end do
   end subroutine ice_atm_forcing_from_atm_impl
end module rdb_ice_atm_forcing
