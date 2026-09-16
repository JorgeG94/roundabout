!! v1 slab-atmosphere restoring filler for the ice atmospheric-forcing seam
!! (SIS2 port, PR 3c; PR 26 adds the fourth field, `atm_fprec`).
module rdb_ice_atm_forcing
   !! Fills ocean_sea_ice_t's coupleable seam
   !! (atm_sf0/atm_dsfdt/atm_sw_dn/atm_fprec) from the &ocean_ice_nml scalar
   !! restoring knobs. SF(T) = lambda*(T - T_air) => sf_0 = -lambda*T_air
   !! (upward-positive SEB intercept), dsf_dt = lambda (thermostat slope),
   !! sw_dn = sw_down (uniform), atm_fprec = snowfall (uniform, PR 26). The
   !! seam is identical whether this filler, a bulk formula, or a live
   !! coupler fills it -- that is what makes the model coupleable
   !! (PLAN_ICE_PR3c S"coupleable seam"; PLAN_PR26_snowfall.md S13 fixes
   !! atm_fprec's shape/units/lifecycle for PR-55 to fill from data).
   !!
   !! v1 SCALAR/UNIFORM: every physical cell gets the same
   !! (sf_0, dsf_dt, sw_dn, fprec). 2D-field / file inputs are the bulk-flux
   !! subsystem's job (deferred). Ghost rows are filled too (harmless --
   !! ice_thermo_columns only steps physical cells). Default-off byte-
   !! identity: never called unless ice%enable.
   use rdb_constants, only: wp
   use rdb_ice_state, only: ocean_sea_ice_t
   implicit none
   private
   public :: ice_atm_forcing_restoring
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
end module rdb_ice_atm_forcing
