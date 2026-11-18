!! Device-side global sums behind the conservation ledger (roundabout-original).
module rdb_ice_ledger_sums
   !! Area-weighted sums over the OWNED cells `is:ie, js:je` of the quantities of
   !! MASTER_PLAN C.5, each kernel a `do concurrent ... reduce(+: ...)` with at
   !! most three accumulators (the proven footprint; wide reductions miscompile
   !! on the device).  Arrays are explicit-shape and device-present on an offload
   !! build; the result is a handful of scalars, the ONLY host<-device traffic of
   !! the ledger.
   !!
   !! **Reproducibility.**  A `reduce(+)` has no fixed summation order, so these
   !! totals are accurate to round-off but not bitwise reproducible across
   !! backends or rank counts.  The ledger reports residuals relative to the
   !! flux scale, never exact zeros; a multi-rank PR replaces them with the EFP
   !! reproducing sums behind the same interface (`rdb_ice_ledger`'s
   !! collective hook).
   !!
   !! Non-finite input propagates through the plain sums; the caller tests the
   !! results with `ieee_is_finite`.
   use rdb_constants, only: wp
   implicit none
   private

   public :: ice_ledger_mass_salt_sums_impl, ice_ledger_energy_sum_impl
   public :: ice_ledger_ocean_flux_sums_impl, ice_ledger_atm_flux_sums_impl

contains

   pure subroutine ice_ledger_mass_salt_sums_impl(nx, ny, ncat, is, ie, js, je, area, vicen, vsnon, &
                                                  rhoi, rhos, sref, mass, salt)
      !! `mass = sum A*(rhoi*vice + rhos*vsno)` [kg], `salt = sum A*rhoi*vice*sref*1e-3`
      !! [kg]; `vice`, `vsno` summed over the categories from the category arrays
      !! (not the cached aggregates, which may be stale).
      integer, intent(in) :: nx, ny, ncat, is, ie, js, je
      real(wp), intent(in) :: area(nx, ny), vicen(nx, ny, ncat), vsnon(nx, ny, ncat)
      real(wp), intent(in) :: rhoi, rhos, sref
      real(wp), intent(out) :: mass, salt

      integer :: i, j, n
      real(wp) :: vi, vs

      mass = 0.0_wp
      salt = 0.0_wp
      !$acc data present(area, vicen, vsnon)
      do concurrent(j=js:je, i=is:ie) local(n, vi, vs) reduce(+:mass, salt)
         vi = 0.0_wp
         vs = 0.0_wp
         !$acc loop seq
         do n = 1, ncat
            vi = vi + vicen(i, j, n)
            vs = vs + vsnon(i, j, n)
         end do
         mass = mass + area(i, j)*(rhoi*vi + rhos*vs)
         salt = salt + area(i, j)*(rhoi*vi*sref*1.0e-3_wp)
      end do
      !$acc end data
   end subroutine ice_ledger_mass_salt_sums_impl

   pure subroutine ice_ledger_energy_sum_impl(nx, ny, ncat, nilyr, nslyr, is, ie, js, je, area, &
                                              vicen, vsnon, qicen, qsnon, energy)
      !! `energy = sum A * sum_n( sum_k qicen*vicen/nilyr + sum_k qsnon*vsnon/nslyr )` [J].
      integer, intent(in) :: nx, ny, ncat, nilyr, nslyr, is, ie, js, je
      real(wp), intent(in) :: area(nx, ny), vicen(nx, ny, ncat), vsnon(nx, ny, ncat)
      real(wp), intent(in) :: qicen(nx, ny, ncat, nilyr), qsnon(nx, ny, ncat, nslyr)
      real(wp), intent(out) :: energy

      integer :: i, j, n, k
      real(wp) :: e

      energy = 0.0_wp
      !$acc data present(area, vicen, vsnon, qicen, qsnon)
      do concurrent(j=js:je, i=is:ie) local(n, k, e) reduce(+:energy)
         e = 0.0_wp
         !$acc loop seq
         do n = 1, ncat
            do k = 1, nilyr
               e = e + qicen(i, j, n, k)*vicen(i, j, n)/real(nilyr, wp)
            end do
            do k = 1, nslyr
               e = e + qsnon(i, j, n, k)*vsnon(i, j, n)/real(nslyr, wp)
            end do
         end do
         energy = energy + area(i, j)*e
      end do
      !$acc end data
   end subroutine ice_ledger_energy_sum_impl

   pure subroutine ice_ledger_ocean_flux_sums_impl(nx, ny, is, ie, js, je, area, fresh, fsalt, fhocn, &
                                                   s_fresh, s_fsalt, s_fhocn)
      !! Area integrals `sum A*fresh`, `sum A*fsalt`, `sum A*fhocn` (flux units times m2).
      integer, intent(in) :: nx, ny, is, ie, js, je
      real(wp), intent(in) :: area(nx, ny), fresh(nx, ny), fsalt(nx, ny), fhocn(nx, ny)
      real(wp), intent(out) :: s_fresh, s_fsalt, s_fhocn

      integer :: i, j

      s_fresh = 0.0_wp
      s_fsalt = 0.0_wp
      s_fhocn = 0.0_wp
      !$acc data present(area, fresh, fsalt, fhocn)
      do concurrent(j=js:je, i=is:ie) reduce(+:s_fresh, s_fsalt, s_fhocn)
         s_fresh = s_fresh + area(i, j)*fresh(i, j)
         s_fsalt = s_fsalt + area(i, j)*fsalt(i, j)
         s_fhocn = s_fhocn + area(i, j)*fhocn(i, j)
      end do
      !$acc end data
   end subroutine ice_ledger_ocean_flux_sums_impl

   pure subroutine ice_ledger_atm_flux_sums_impl(nx, ny, is, ie, js, je, Lvap, area, fswthru, evap, fswabs, &
                                                 flwout, fsens, s_fswthru, s_evap, s_atm)
      !! Icepack signs.  `sum A*fswthru`, `sum A*evap` (+ = condensation, mass
      !! GAIN) and the net atmosphere heat into the ice
      !! `sum A*(fswabs + flwout + fsens + Lvap*evap)`: `flwout` and `fsens` are
      !! positive DOWNWARD and the merged `flwout` is already net of the reflected
      !! incoming longwave.  `flat` is deliberately absent: it cancels in
      !! upstream's `fsurf - flat`.  The `+ Lvap*evap` term is there because
      !! upstream's final column energy carries the vapour at `-Lvap*evap*dt`
      !! and the energy check closes ice plus vapour, so the stored-state
      !! enthalpy changes by `einp + Lvap*evap*dt`.  The incoming-longwave and
      !! `fswthru` terms are added by the caller (`H_ext`) and the ledger.
      integer, intent(in) :: nx, ny, is, ie, js, je
      real(wp), intent(in) :: Lvap
      real(wp), intent(in) :: area(nx, ny), fswthru(nx, ny), evap(nx, ny), fswabs(nx, ny)
      real(wp), intent(in) :: flwout(nx, ny), fsens(nx, ny)
      real(wp), intent(out) :: s_fswthru, s_evap, s_atm

      integer :: i, j

      s_fswthru = 0.0_wp
      s_evap = 0.0_wp
      s_atm = 0.0_wp
      !$acc data present(area, fswthru, evap, fswabs, flwout, fsens)
      do concurrent(j=js:je, i=is:ie) reduce(+:s_fswthru, s_evap, s_atm)
         s_fswthru = s_fswthru + area(i, j)*fswthru(i, j)
         s_evap = s_evap + area(i, j)*evap(i, j)
         s_atm = s_atm + area(i, j)*(((fswabs(i, j) + flwout(i, j)) + fsens(i, j)) + Lvap*evap(i, j))
      end do
      !$acc end data
   end subroutine ice_ledger_atm_flux_sums_impl

end module rdb_ice_ledger_sums
