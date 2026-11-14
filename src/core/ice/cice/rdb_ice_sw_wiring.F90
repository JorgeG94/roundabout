!! Shortwave wiring of the CICE-port sea ice: radiation -> vertical thermodynamics, and the
!! transmitted-shortwave flux producer (roundabout-original sequencing; no column physics).
module rdb_ice_sw_wiring
   !! One ice step's shortwave sequence, in the order the radiation driver documents:
   !!
   !! 1. `ice_step_radiation` (after the forcing import and `ice_save_init`, so `aicen ==
   !!    aicen_init`): fills the slot `rad` (`fswsfcn, fswintn, fswthrun, Iswabsn`).
   !! 2. `ice_step_therm1` on the slot's per-category arrays plus a snow absorption array that
   !!    is ZERO: under the CCSM3 scheme (Briegleb et al. 2004) snow absorbs no shortwave
   !!    itself (the visible light that penetrates a snow cover is dropped, the snow albedo
   !!    carries the whole surface part), so `Sswabs = 0`.  The array is owned ONCE here
   !!    (`ice_sw_wiring_workspace_ensure`, device-mapped, never allocated per call); the
   !!    `sw_redist` pass can only write zero back into it.
   !! 3. The transmitted-shortwave producer `ice_sw_thru_merge_impl`: `fswthru += sum_n
   !!    aicen_init(n)*fswthrun(n)`, `fswthru_vis` likewise (one per-cell launch, fixed
   !!    category order, the weights are the same start-of-step areas as therm1's flux merge).
   !!
   !! **`fswthrun` is not changed by the thermodynamics.**  Both here (`intent(in)` all the way
   !! down) and upstream (`thermo_vertical` takes it `intent(in)`; the redistribution moves
   !! absorption between the surface and the layers, never to or from the through-ice flux; a
   !! melted-through column's leftover energy goes to `fhocn`, not into `fswthru`), so the
   !! producer reads the radiation value and the ordering against therm1 does not matter.  It
   !! runs after therm1 only to see which categories therm1 committed.
   !!
   !! **Which categories contribute.**  A category contributes iff therm1 committed it, read
   !! off its own output: `fswabsn > 0`.  A skipped category (ice-free, a cell carrying a fatal
   !! error bit before the step, a flagged column) has `fswabsn = 0` and adds nothing to the
   !! `fswabs` merge, so it must add nothing here either or the heat ledger (`fswabs -
   !! fswthru`) would not close.  A committed column with `fswabsn = 0` has `fswthrun = 0`
   !! too (`fswabsn = fswsfcn + fswintn + fswthrun`, all non-negative), so the gate loses
   !! nothing.  A non-finite `fswthrun` of a committed column cannot occur (therm1 flags it).
   !!
   !! **Bands.**  Under this scheme all transmitted light is visible (near-IR is absorbed in
   !! the top ice millimetres), so `fswthru_vis == fswthru` bitwise; the second accumulator is
   !! kept separate so a band-resolved scheme replaces one line.
   !!
   !! **Errors.**  `ierr` is the radiation code (`ICE_RAD_*`; nothing ran, therm1 included,
   !! `istat = -1`); a corrupt cell is flagged by the radiation aggregate in `ice_err` and
   !! therm1 then skips that cell whole.  `istat` is therm1's.  Neither may be read as "no
   !! shortwave".
   use, intrinsic :: iso_fortran_env, only: int32
   use rdb_constants, only: wp
   use rdb_ice_params, only: ice_params_t, ice_sw_params_t
   use rdb_ice_rad_state, only: ice_rad_state_t
   use rdb_ice_radiation, only: ice_step_radiation, ICE_RAD_OK
   use rdb_ice_therm1, only: ice_step_therm1
   implicit none
   private

   public :: ice_sw_thru_merge_impl, ice_step_sw_therm1
   public :: ice_sw_wiring_workspace_ensure, ice_sw_wiring_cleanup, ice_sw_wiring_zero_snow_ok

   real(wp), allocatable, save :: sswabs_zero(:, :, :, :)
      !! The persistent zero snow absorption `(nx, ny, ncat, nslyr)`, device-mapped.

contains

   subroutine ice_sw_wiring_workspace_ensure(nx, ny, ncat, nslyr)
      !! Allocate (zero) and device-map the snow-absorption array for these extents; a no-op
      !! when it already has them.  Host-side; call at setup or lazily from the shim.
      integer, intent(in) :: nx, ny, ncat, nslyr

      if (allocated(sswabs_zero)) then
         if (all(shape(sswabs_zero) == [nx, ny, ncat, nslyr])) return
         call ice_sw_wiring_cleanup()
      end if
      allocate (sswabs_zero(nx, ny, ncat, nslyr), source=0.0_wp)
      !$acc enter data copyin(sswabs_zero)
   end subroutine ice_sw_wiring_workspace_ensure

   subroutine ice_sw_wiring_cleanup()
      !! Release the device copy and the array.
      if (allocated(sswabs_zero)) then
         !$acc exit data delete(sswabs_zero)
         deallocate (sswabs_zero)
      end if
   end subroutine ice_sw_wiring_cleanup

   function ice_sw_wiring_zero_snow_ok() result(ok)
      !! Host check (test hook): the DEVICE copy of the snow absorption is still exactly zero.
      logical :: ok

      ok = allocated(sswabs_zero)
      if (.not. ok) return
      !$acc update self(sswabs_zero)
      ok = all(sswabs_zero == 0.0_wp)
   end function ice_sw_wiring_zero_snow_ok

   pure subroutine ice_sw_thru_merge_impl(nx, ny, ncat, is, ie, js, je, puny, aicen_init, fswthrun, fswabsn, &
                                          fswthru, fswthru_vis)
      !! ONE per-cell launch over `is:ie, js:je`: `fswthru += sum_n aicen_init*fswthrun` over the
      !! categories with `aicen_init > puny` that therm1 committed (`fswabsn > 0`, module header),
      !! in the fixed order `n = 1..ncat`; `fswthru_vis` the same sum (all transmitted light is
      !! visible).  ADDITIVE; each work item writes only its own cell.  Arrays device-present.
      integer, intent(in) :: nx, ny, ncat, is, ie, js, je
      real(wp), intent(in) :: puny
      real(wp), intent(in) :: aicen_init(nx, ny, ncat), fswthrun(nx, ny, ncat), fswabsn(nx, ny, ncat)
      real(wp), intent(inout) :: fswthru(nx, ny), fswthru_vis(nx, ny)
      integer :: i, j, n
      real(wp) :: a, s

      !$acc data present(aicen_init, fswthrun, fswabsn, fswthru, fswthru_vis)
      do concurrent(j=js:je, i=is:ie) local(n, a, s)
         s = 0.0_wp
         do n = 1, ncat
            a = aicen_init(i, j, n)
            if (a > puny .and. fswabsn(i, j, n) > 0.0_wp) s = s + a*fswthrun(i, j, n)
         end do
         fswthru(i, j) = fswthru(i, j) + s
         fswthru_vis(i, j) = fswthru_vis(i, j) + s
      end do
      !$acc end data
   end subroutine ice_sw_thru_merge_impl

   subroutine ice_step_sw_therm1(rad, swp, nx, ny, ncat, nilyr, nslyr, is, ie, js, je, dt, p, sw_redist, zlvl, zlvs, &
                                 sw_down, uatm, vatm, wind, potT, Qa, rhoa, flw, fsnow, uvel, vvel, Tf, sst, frzmlt, &
                                 tau_x, tau_y, aicen_init, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, &
                                 fresh, fsalt, fhocn, evap, fsens, flat, flwout, fswabs, meltt, meltb, congel, snoice, &
                                 freshn, fsaltn, fhocnn, evapn, fsensn, flatn, flwoutn, fswabsn, &
                                 melttn, meltbn, congeln, snoicen, strairxn, strairyn, strair_x_T, strair_y_T, &
                                 fbot, rside, wlat, ierr_n, ice_err, fswthru, fswthru_vis, ierr, istat)
      !! Radiation -> therm1 -> transmitted-shortwave producer for the cells `is:ie, js:je` (module
      !! header).  The arguments are `ice_step_therm1`'s with the five shortwave inputs replaced
      !! by the radiation slot `rad` and the broadband `sw_down` [W/m2], plus the two cell
      !! accumulators `fswthru, fswthru_vis` (ADDED to; zeroed once per step by `ice_flux_zero`).
      !! Every array must be device-mapped; nothing is gathered to the host.
      !! `ierr /= ICE_RAD_OK` (radiation refused; therm1 NOT run, `istat = -1`) or `istat /= 0`
      !! (therm1 refused its configuration) is fatal for the step.
      type(ice_rad_state_t), intent(inout) :: rad
      type(ice_sw_params_t), intent(in) :: swp
      integer, intent(in) :: nx, ny, ncat, nilyr, nslyr, is, ie, js, je
      real(wp), intent(in) :: dt, zlvl, zlvs
      type(ice_params_t), intent(in) :: p
      logical, intent(in) :: sw_redist
      real(wp), intent(in) :: sw_down(nx, ny)
      real(wp), intent(in) :: uatm(nx, ny), vatm(nx, ny), wind(nx, ny), potT(nx, ny), Qa(nx, ny), rhoa(nx, ny)
      real(wp), intent(in) :: flw(nx, ny), fsnow(nx, ny), uvel(nx, ny), vvel(nx, ny), Tf(nx, ny), sst(nx, ny)
      real(wp), intent(in) :: frzmlt(nx, ny), tau_x(nx, ny), tau_y(nx, ny)
      real(wp), intent(in) :: aicen_init(nx, ny, ncat)
      real(wp), intent(inout) :: aicen(nx, ny, ncat), vicen(nx, ny, ncat), vsnon(nx, ny, ncat), Tsfcn(nx, ny, ncat)
      real(wp), intent(inout) :: qicen(nx, ny, ncat, nilyr), qsnon(nx, ny, ncat, nslyr)
      real(wp), intent(in) :: sicen(nx, ny, ncat, nilyr)
      real(wp), intent(inout) :: fresh(nx, ny), fsalt(nx, ny), fhocn(nx, ny), evap(nx, ny), fsens(nx, ny), flat(nx, ny)
      real(wp), intent(inout) :: flwout(nx, ny), fswabs(nx, ny), meltt(nx, ny), meltb(nx, ny), congel(nx, ny), snoice(nx, ny)
      real(wp), intent(inout) :: freshn(nx, ny, ncat), fsaltn(nx, ny, ncat), fhocnn(nx, ny, ncat), evapn(nx, ny, ncat)
      real(wp), intent(inout) :: fsensn(nx, ny, ncat), flatn(nx, ny, ncat), flwoutn(nx, ny, ncat), fswabsn(nx, ny, ncat)
      real(wp), intent(inout) :: melttn(nx, ny, ncat), meltbn(nx, ny, ncat), congeln(nx, ny, ncat), snoicen(nx, ny, ncat)
      real(wp), intent(inout) :: strairxn(nx, ny, ncat), strairyn(nx, ny, ncat), strair_x_T(nx, ny), strair_y_T(nx, ny)
      real(wp), intent(inout) :: fbot(nx, ny), rside(nx, ny), wlat(nx, ny)
      integer(int32), intent(inout) :: ierr_n(nx, ny, ncat), ice_err(nx, ny)
      real(wp), intent(inout) :: fswthru(nx, ny), fswthru_vis(nx, ny)
      integer, intent(out) :: ierr, istat

      istat = -1
      call ice_step_radiation(rad, swp, is, ie, js, je, nx, ny, ncat, nilyr, aicen, vicen, vsnon, Tsfcn, sw_down, &
                              ice_err, ierr)
      if (ierr /= ICE_RAD_OK) return
      call ice_sw_wiring_workspace_ensure(nx, ny, ncat, nslyr)
      call ice_step_therm1(nx, ny, ncat, nilyr, nslyr, is, ie, js, je, dt, p, sw_redist, zlvl, zlvs, &
                           uatm, vatm, wind, potT, Qa, rhoa, flw, fsnow, uvel, vvel, Tf, sst, frzmlt, tau_x, tau_y, &
                           aicen_init, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, &
                           rad%fswsfcn, rad%fswintn, rad%fswthrun, rad%Iswabsn, sswabs_zero, &
                           fresh, fsalt, fhocn, evap, fsens, flat, flwout, fswabs, meltt, meltb, congel, snoice, &
                           freshn, fsaltn, fhocnn, evapn, fsensn, flatn, flwoutn, fswabsn, &
                           melttn, meltbn, congeln, snoicen, strairxn, strairyn, strair_x_T, strair_y_T, &
                           fbot, rside, wlat, ierr_n, ice_err, istat)
      if (istat /= 0) return
      call ice_sw_thru_merge_impl(nx, ny, ncat, is, ie, js, je, p%c%puny, aicen_init, rad%fswthrun, fswabsn, &
                                  fswthru, fswthru_vis)
   end subroutine ice_step_sw_therm1

end module rdb_ice_sw_wiring
