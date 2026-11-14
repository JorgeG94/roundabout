!! Flat kernels of the ocean-ice flux accumulators (roundabout-original).
module rdb_ice_flux
   !! Explicit-shape `do concurrent` bodies behind `ice_flux_zero`,
   !! `ice_flux_merge_categories` and `ice_save_init` (the type-aware shims live in
   !! `rdb_ice_flux_state`).  Nothing here needs the ocean grid: every array
   !! arrives as an explicit-shape dummy, device-present on an offload build.
   !!
   !! **Why a per-category set exists.**  The vertical thermodynamics, shortwave
   !! and atmosphere kernels run category-parallel, `do concurrent (n, j, i)`;
   !! different `n` would add into one `(nx, ny)` cell (a race, and atomics would
   !! make the sum order non-reproducible).  They therefore write the
   !! per-category outputs `(nx, ny, ncat)` and `ice_flux_merge_impl` folds them
   !! into the cell accumulators in ONE per-cell kernel with the category loop
   !! INLINE in the fixed order `n = 1..ncat` (a reproducible sum), each
   !! category weighted by its start-of-step area `aicen_init` (the weight is the
   !! PRE-thermodynamics area: a category that melts out is still weighted).
   !!
   !! **Kernels.**
   !! - `ice_flux_zero_impl`: every accumulator and every per-category output,
   !!   ghosts included, in ONE launch (a single `do concurrent (j, i)` with the
   !!   category loop inside).
   !! - `ice_flux_merge_impl`: ADDS `sum_n aicen_init*xn` into the twelve cell
   !!   accumulators over the cells `is:ie, js:je`; pre-existing cell values are
   !!   kept.  A non-finite contribution propagates into the cell fields (no
   !!   clamp can launder it) and ORs `ICE_ERR_NONFINITE` into `ice_err`.
   !! - `ice_save_init_impl`: the start-of-step snapshots
   !!   `aicen/vicen/vsnon/Tsfcn/aice -> *_init`, ghosts included, ONE launch.
   !!
   !! No kernel here calls a procedure with a loop index, so the by-value index
   !! rule (nvfortran 26.5 `do concurrent` bug) is satisfied by construction.
   use, intrinsic :: iso_fortran_env, only: int32
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_constants, only: wp
   use rdb_ice_errors, only: ICE_ERR_NONFINITE
   implicit none
   private

   public :: ice_flux_zero_impl, ice_flux_merge_impl, ice_save_init_impl

contains

   pure subroutine ice_flux_zero_impl(nx, ny, ncat, &
                                      fresh, fsalt, fhocn, fswthru, fswthru_vis, evap, fsens, flat, flwout, fswabs, &
                                      frazil_vol, meltt, meltb, meltl, congel, snoice, &
                                      dardg1dt, dardg2dt, dvirdgdt, opening, &
                                      freshn, fsaltn, fhocnn, evapn, fsensn, flatn, flwoutn, fswabsn, &
                                      melttn, meltbn, congeln, snoicen)
      !! Zero the 20 cell accumulators and the 12 per-category outputs over the
      !! whole ghosted extent, in one launch.
      integer, intent(in) :: nx, ny, ncat
      real(wp), intent(inout) :: fresh(nx, ny), fsalt(nx, ny), fhocn(nx, ny), fswthru(nx, ny)
      real(wp), intent(inout) :: fswthru_vis(nx, ny), evap(nx, ny), fsens(nx, ny), flat(nx, ny)
      real(wp), intent(inout) :: flwout(nx, ny), fswabs(nx, ny), frazil_vol(nx, ny), meltt(nx, ny)
      real(wp), intent(inout) :: meltb(nx, ny), meltl(nx, ny), congel(nx, ny), snoice(nx, ny)
      real(wp), intent(inout) :: dardg1dt(nx, ny), dardg2dt(nx, ny), dvirdgdt(nx, ny), opening(nx, ny)
      real(wp), intent(inout) :: freshn(nx, ny, ncat), fsaltn(nx, ny, ncat), fhocnn(nx, ny, ncat)
      real(wp), intent(inout) :: evapn(nx, ny, ncat), fsensn(nx, ny, ncat), flatn(nx, ny, ncat)
      real(wp), intent(inout) :: flwoutn(nx, ny, ncat), fswabsn(nx, ny, ncat), melttn(nx, ny, ncat)
      real(wp), intent(inout) :: meltbn(nx, ny, ncat), congeln(nx, ny, ncat), snoicen(nx, ny, ncat)

      integer :: i, j, n

      !$acc data present(fresh, fsalt, fhocn, fswthru, fswthru_vis, evap, fsens, flat, flwout, fswabs) &
      !$acc& present(frazil_vol, meltt, meltb, meltl, congel, snoice, dardg1dt, dardg2dt, dvirdgdt, opening) &
      !$acc& present(freshn, fsaltn, fhocnn, evapn, fsensn, flatn, flwoutn, fswabsn) &
      !$acc& present(melttn, meltbn, congeln, snoicen)
      do concurrent(j=1:ny, i=1:nx) local(n)
         fresh(i, j) = 0.0_wp
         fsalt(i, j) = 0.0_wp
         fhocn(i, j) = 0.0_wp
         fswthru(i, j) = 0.0_wp
         fswthru_vis(i, j) = 0.0_wp
         evap(i, j) = 0.0_wp
         fsens(i, j) = 0.0_wp
         flat(i, j) = 0.0_wp
         flwout(i, j) = 0.0_wp
         fswabs(i, j) = 0.0_wp
         frazil_vol(i, j) = 0.0_wp
         meltt(i, j) = 0.0_wp
         meltb(i, j) = 0.0_wp
         meltl(i, j) = 0.0_wp
         congel(i, j) = 0.0_wp
         snoice(i, j) = 0.0_wp
         dardg1dt(i, j) = 0.0_wp
         dardg2dt(i, j) = 0.0_wp
         dvirdgdt(i, j) = 0.0_wp
         opening(i, j) = 0.0_wp
         !$acc loop seq
         do n = 1, ncat
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
         end do
      end do
      !$acc end data
   end subroutine ice_flux_zero_impl

   pure subroutine ice_flux_merge_impl(nx, ny, ncat, is, ie, js, je, emis, puny, flw, aicen_init, &
                                       freshn, fsaltn, fhocnn, evapn, fsensn, flatn, flwoutn, fswabsn, &
                                       melttn, meltbn, congeln, snoicen, &
                                       fresh, fsalt, fhocn, evap, fsens, flat, flwout, fswabs, &
                                       meltt, meltb, congel, snoice, ice_err)
      !! `x(i,j) = x(i,j) + sum_{n=1..ncat} aicen_init(i,j,n) * xn(i,j,n)` for the
      !! twelve merged fields, over the cells `is:ie, js:je`, with ONE exception:
      !! the outgoing longwave merges the portion that is not reflected incoming
      !! longwave, `flwout += sum_n aicen_init*(flwoutn - (1 - emis)*flw)` (the
      !! upstream merge, which only merges categories with `aicen_init > puny`: the reflected
      !! term is added for those only; `flw` is the cell's incoming longwave, `emis` the
      !! surface emissivity).  The sum is formed first (fixed category order), then added
      !! to the cell value.  The result is valid on the cells `is:ie, js:je` ONLY.
      integer, intent(in) :: nx, ny, ncat, is, ie, js, je
      real(wp), intent(in) :: emis, puny
      real(wp), intent(in) :: flw(nx, ny), aicen_init(nx, ny, ncat)
      real(wp), intent(in) :: freshn(nx, ny, ncat), fsaltn(nx, ny, ncat), fhocnn(nx, ny, ncat)
      real(wp), intent(in) :: evapn(nx, ny, ncat), fsensn(nx, ny, ncat), flatn(nx, ny, ncat)
      real(wp), intent(in) :: flwoutn(nx, ny, ncat), fswabsn(nx, ny, ncat), melttn(nx, ny, ncat)
      real(wp), intent(in) :: meltbn(nx, ny, ncat), congeln(nx, ny, ncat), snoicen(nx, ny, ncat)
      real(wp), intent(inout) :: fresh(nx, ny), fsalt(nx, ny), fhocn(nx, ny), evap(nx, ny)
      real(wp), intent(inout) :: fsens(nx, ny), flat(nx, ny), flwout(nx, ny), fswabs(nx, ny)
      real(wp), intent(inout) :: meltt(nx, ny), meltb(nx, ny), congel(nx, ny), snoice(nx, ny)
      integer(int32), intent(inout) :: ice_err(nx, ny)

      integer :: i, j, n
      real(wp) :: a, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, chk

      !$acc data present(flw, aicen_init, freshn, fsaltn, fhocnn, evapn, fsensn, flatn, flwoutn, fswabsn) &
      !$acc& present(melttn, meltbn, congeln, snoicen) &
      !$acc& present(fresh, fsalt, fhocn, evap, fsens, flat, flwout, fswabs) &
      !$acc& present(meltt, meltb, congel, snoice, ice_err)
      do concurrent(j=js:je, i=is:ie) local(n, a, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, chk)
         s1 = 0.0_wp
         s2 = 0.0_wp
         s3 = 0.0_wp
         s4 = 0.0_wp
         s5 = 0.0_wp
         s6 = 0.0_wp
         s7 = 0.0_wp
         s8 = 0.0_wp
         s9 = 0.0_wp
         s10 = 0.0_wp
         s11 = 0.0_wp
         s12 = 0.0_wp
         !$acc loop seq
         do n = 1, ncat
            a = aicen_init(i, j, n)
            s1 = s1 + a*freshn(i, j, n)
            s2 = s2 + a*fsaltn(i, j, n)
            s3 = s3 + a*fhocnn(i, j, n)
            s4 = s4 + a*evapn(i, j, n)
            s5 = s5 + a*fsensn(i, j, n)
            s6 = s6 + a*flatn(i, j, n)
            if (a > puny) then   ! the reflected incoming longwave belongs to categories the thermodynamics ran
               s7 = s7 + a*(flwoutn(i, j, n) - (1.0_wp - emis)*flw(i, j))
            end if
            s8 = s8 + a*fswabsn(i, j, n)
            s9 = s9 + a*melttn(i, j, n)
            s10 = s10 + a*meltbn(i, j, n)
            s11 = s11 + a*congeln(i, j, n)
            s12 = s12 + a*snoicen(i, j, n)
         end do
         fresh(i, j) = fresh(i, j) + s1
         fsalt(i, j) = fsalt(i, j) + s2
         fhocn(i, j) = fhocn(i, j) + s3
         evap(i, j) = evap(i, j) + s4
         fsens(i, j) = fsens(i, j) + s5
         flat(i, j) = flat(i, j) + s6
         flwout(i, j) = flwout(i, j) + s7
         fswabs(i, j) = fswabs(i, j) + s8
         meltt(i, j) = meltt(i, j) + s9
         meltb(i, j) = meltb(i, j) + s10
         congel(i, j) = congel(i, j) + s11
         snoice(i, j) = snoice(i, j) + s12
         chk = ((s1 + s2) + (s3 + s4)) + ((s5 + s6) + (s7 + s8)) + ((s9 + s10) + (s11 + s12))
         if (.not. ieee_is_finite(chk)) ice_err(i, j) = ior(ice_err(i, j), ICE_ERR_NONFINITE)
      end do
      !$acc end data
   end subroutine ice_flux_merge_impl

   pure subroutine ice_save_init_impl(nx, ny, ncat, aicen, vicen, vsnon, Tsfcn, aice, &
                                      aicen_init, vicen_init, vsnon_init, Tsfcn_init, aice_init)
      !! Copy the category state and the cell aggregate into their start-of-step
      !! twins over the whole ghosted extent, in one launch.
      integer, intent(in) :: nx, ny, ncat
      real(wp), intent(in) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat)
      real(wp), intent(in) :: Tsfcn(nx, ny, ncat), aice(nx, ny)
      real(wp), intent(inout) :: aicen_init(nx, ny, ncat), vicen_init(nx, ny, ncat)
      real(wp), intent(inout) :: vsnon_init(nx, ny, ncat), Tsfcn_init(nx, ny, ncat), aice_init(nx, ny)

      integer :: i, j, n

      !$acc data present(aicen, vicen, vsnon, Tsfcn, aice, aicen_init, vicen_init, vsnon_init, Tsfcn_init, aice_init)
      do concurrent(j=1:ny, i=1:nx) local(n)
         aice_init(i, j) = aice(i, j)
         !$acc loop seq
         do n = 1, ncat
            aicen_init(i, j, n) = aicen(i, j, n)
            vicen_init(i, j, n) = vicen(i, j, n)
            vsnon_init(i, j, n) = vsnon(i, j, n)
            Tsfcn_init(i, j, n) = Tsfcn(i, j, n)
         end do
      end do
      !$acc end data
   end subroutine ice_save_init_impl

end module rdb_ice_flux
