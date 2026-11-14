!! The SHIPPED conservation ledger (`ice_ledger_mirror_t`) closing a 1000-step seasonal run of the vertical
!! thermodynamics step on the real state and flux slots (D3.4).
!!
!! Five categories, a cold start, a sunny melt season with ocean heat, snowfall and `sw_redist`, driven through
!! `ice_step_therm1`; every step is bracketed by `ice_flux_zero`, `ice_save_init`, `led%open` and `led%close`.
!! Water and salt close to round-off.  The heat line closes to the BL99 solver tolerance: a pass is accepted at
!! `ferr <= 0.9 ferrmax` and the energy check at `1.1 ferrmax` (`ferrmax = 1e-3 W/m2` per ice area), so the cell
!! residual of a step is bounded by `1.1 ferrmax dt A sum(aicen_init)`; the worst ratio to that bound is reported.
!! The heat line carries the vapour term `+Lvap*evap` (upstream's final column energy holds the vapour at
!! `-Lvap*evap*dt`); without it the same run leaves a residual of about 1e-4 of the throughput.
module test_cice_therm1_ledger
   use, intrinsic :: iso_fortran_env, only: int32, output_unit
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ice_params, only: ice_params_t, ice_params_default, ICE_CONSTANTS_ICEPACK, ice_q_from_T, ice_qsno_from_T, &
                             ice_salinity_profile
   use rdb_ice_cice_state, only: ice_cice_state_t
   use rdb_ice_flux_state, only: ice_flux_t, ice_flux_zero, ice_save_init
   use rdb_ice_ledger, only: ice_ledger_mirror_t, ICE_LEDGER_WATER, ICE_LEDGER_SALT, ICE_LEDGER_HEAT
   use rdb_ice_therm1, only: ice_step_therm1
   use rdb_ice_thickness, only: ICE_FERRMAX
   implicit none
   private

   public :: collect_cice_therm1_ledger_tests

   integer, parameter :: NG = 1, NCAT = 5, NIL = 4, NX = 3, NY = 3, IC = 2, JC = 2
   real(wp), parameter :: TWOPI = 6.283185307179586_wp, AREA = 1.0e6_wp, DT = 3600.0_wp

contains

   subroutine collect_cice_therm1_ledger_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("shipped_ledger_closes_1000_step_season", test_season)]
   end subroutine collect_cice_therm1_ledger_tests

   subroutine test_season(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: g
      type(ice_params_t) :: p
      type(ice_cice_state_t) :: st
      type(ice_flux_t) :: fx
      type(ice_ledger_mirror_t) :: led
      real(wp) :: uatm(NX, NY), vatm(NX, NY), wind(NX, NY), potT(NX, NY), Qa(NX, NY), rhoa(NX, NY), flw(NX, NY)
      real(wp) :: fsnow(NX, NY), uvel(NX, NY), vvel(NX, NY), Tf(NX, NY), sst(NX, NY), frzmlt(NX, NY), tau_x(NX, NY)
      real(wp) :: tau_y(NX, NY), sx(NX, NY), sy(NX, NY), fbot(NX, NY), rside(NX, NY), wlat(NX, NY), area_a(NX, NY)
      real(wp) :: fswsfcn(NX, NY, NCAT), fswintn(NX, NY, NCAT), fswthrun(NX, NY, NCAT), sxn(NX, NY, NCAT), syn(NX, NY, NCAT)
      real(wp) :: Isw(NX, NY, NCAT, NIL), Ssw(NX, NY, NCAT, 1)
      integer(int32) :: ierr_n(NX, NY, NCAT)
      real(wp), parameter :: A(NCAT) = [0.05_wp, 0.15_wp, 0.30_wp, 0.25_wp, 0.10_wp]
      real(wp), parameter :: H(NCAT) = [0.12_wp, 0.45_wp, 1.0_wp, 1.9_wp, 3.1_wp]
      real(wp), parameter :: HS(NCAT) = [0.0_wp, 0.03_wp, 0.08_wp, 0.15_wp, 0.2_wp]
      integer, parameter :: NSTEP = 1000
      integer :: ierr, istat, s, n, k, nmelt, nsnow
      real(wp) :: Sp(NIL + 1), Tz, ph, Tair, SWr, qs, wk(NIL), wsum, a0, fthr, w_ext, h_ext, bound
      real(wp) :: worst_w, worst_s, worst_h, worst_ratio, rw, rs, rh

      call g%init(1, 1, NG, 1000.0_wp, 1000.0_wp)
      p = ice_params_default(ICE_CONSTANTS_ICEPACK, NCAT, NIL, 1)
      call st%init(g, p, ierr)
      call check(error, ierr, 0, "state init")
      if (allocated(error)) return
      call fx%init(g, NCAT, ierr)
      call check(error, ierr, 0, "flux init")
      if (allocated(error)) return
      call ice_salinity_profile(NIL, p%c%saltmax, p%c%min_salin, Sp)
      do n = 1, NCAT
         st%aicen(IC, JC, n) = A(n)
         st%vicen(IC, JC, n) = A(n)*H(n)
         st%vsnon(IC, JC, n) = A(n)*HS(n)
         st%Tsfcn(IC, JC, n) = -9.0_wp
         do k = 1, NIL
            st%sicen(IC, JC, n, k) = Sp(k)
            Tz = -7.0_wp + 5.0_wp*real(k - 1, wp)/real(NIL - 1, wp)
            st%qicen(IC, JC, n, k) = ice_q_from_T(Tz, Sp(k), p%c, p%thermo%l_brine)
         end do
         st%qsnon(IC, JC, n, 1) = ice_qsno_from_T(-4.5_wp, p%c)
      end do
      area_a = AREA
      uatm = 6.0_wp; vatm = 0.0_wp; wind = 6.0_wp; rhoa = 1.3_wp; uvel = 0.0_wp; vvel = 0.0_wp
      Tf = -1.8_wp; tau_x = 0.03_wp; tau_y = 0.0_wp
      fswsfcn = 0.0_wp; fswintn = 0.0_wp; fswthrun = 0.0_wp; Isw = 0.0_wp; Ssw = 0.0_wp
      sxn = 0.0_wp; syn = 0.0_wp; sx = 0.0_wp; sy = 0.0_wp; fbot = 0.0_wp; rside = 0.0_wp; wlat = 0.0_wp
      ierr_n = 0_int32
      potT = 263.0_wp; Qa = 6.0e-4_wp; flw = 200.0_wp; fsnow = 0.0_wp; sst = -1.8_wp; frzmlt = 0.0_wp

      call st%enter_data()
      call fx%enter_data()
      !$acc enter data copyin(area_a, uatm, vatm, wind, potT, Qa, rhoa, flw, fsnow, uvel, vvel, Tf, sst, frzmlt, tau_x, tau_y) &
      !$acc& copyin(sx, sy, fbot, rside, wlat, fswsfcn, fswintn, fswthrun, sxn, syn, Isw, Ssw, ierr_n)
      worst_w = 0.0_wp
      worst_s = 0.0_wp
      worst_h = 0.0_wp
      worst_ratio = 0.0_wp
      nmelt = 0
      nsnow = 0
      do s = 1, NSTEP
         ph = TWOPI*real(s - 1, wp)/real(NSTEP, wp)
         Tair = 266.0_wp - 18.0_wp*cos(ph)
         SWr = 260.0_wp*(0.5_wp - 0.5_wp*cos(ph))**2
         potT = Tair
         qs = 11637800.0_wp*exp(-5897.8_wp/Tair)/1.3_wp
         Qa = 0.8_wp*qs
         flw = 0.8_wp*5.67e-8_wp*Tair**4
         fsnow = merge(1.5e-4_wp, 0.0_wp, Tair < 270.0_wp .and. mod(s, 50) < 5)
         frzmlt = -30.0_wp*max(0.0_wp, -cos(ph))
         sst = Tf + merge(0.2_wp, 0.0_wp, SWr > 100.0_wp)
         wsum = 0.0_wp
         do k = 1, NIL
            wk(k) = 0.5_wp**k
            wsum = wsum + wk(k)
         end do
         do n = 1, NCAT
            fswsfcn(:, :, n) = 0.25_wp*SWr
            fswthrun(:, :, n) = 0.03_wp*SWr*real(n, wp)/5.0_wp
            fswintn(:, :, n) = 0.08_wp*SWr
            do k = 1, NIL
               Isw(:, :, n, k) = 0.08_wp*SWr*wk(k)/wsum
            end do
         end do
         !$acc update device(potT, Qa, flw, fsnow, frzmlt, sst, fswsfcn, fswthrun, fswintn, Isw)

         call ice_flux_zero(fx)
         call ice_save_init(st)
         call led%open(st, area_a, DT, 1)
         !$acc update self(st%aicen_init)
         call ice_step_therm1(st%nx, st%ny, NCAT, NIL, 1, IC, IC, JC, JC, DT, st%p, .true., 10.0_wp, 10.0_wp, &
                              uatm, vatm, wind, potT, Qa, rhoa, flw, fsnow, uvel, vvel, Tf, sst, frzmlt, tau_x, tau_y, &
                              st%aicen_init, st%aicen, st%vicen, st%vsnon, st%Tsfcn, st%qicen, st%sicen, st%qsnon, &
                              fswsfcn, fswintn, fswthrun, Isw, Ssw, &
                              fx%fresh, fx%fsalt, fx%fhocn, fx%evap, fx%fsens, fx%flat, fx%flwout, fx%fswabs, fx%meltt, &
                              fx%meltb, fx%congel, fx%snoice, &
                              fx%freshn, fx%fsaltn, fx%fhocnn, fx%evapn, fx%fsensn, fx%flatn, fx%flwoutn, fx%fswabsn, &
                              fx%melttn, fx%meltbn, fx%congeln, fx%snoicen, sxn, syn, sx, sy, fbot, rside, wlat, &
                              ierr_n, st%ice_err, istat)
         call check(error, istat, 0, "step accepted")
         if (allocated(error)) exit
         ! the radiation driver's aggregate: fswthru = sum aicen_init * fswthrun
         a0 = 0.0_wp
         fthr = 0.0_wp
         do n = 1, NCAT
            if (st%aicen_init(IC, JC, n) > p%c%puny) then
               a0 = a0 + st%aicen_init(IC, JC, n)
               fthr = fthr + st%aicen_init(IC, JC, n)*fswthrun(IC, JC, n)
            end if
         end do
         fx%fswthru(IC, JC) = fthr
         !$acc update device(fx%fswthru)
         w_ext = DT*AREA*fsnow(IC, JC)*a0
         h_ext = DT*AREA*(a0*flw(IC, JC) - p%c%Lfresh*fsnow(IC, JC)*a0)
         call led%close(st, fx, area_a, w_ext, h_ext)
         !$acc update self(st%ice_err)
         if (st%ice_err(IC, JC) /= 0_int32) then
            call test_failed(error, "error word raised during the season")
            exit
         end if
         if (.not. led%valid) then
            call test_failed(error, "the ledger did not close")
            exit
         end if
         rw = led%rel_residual(ICE_LEDGER_WATER)
         rs = led%rel_residual(ICE_LEDGER_SALT)
         rh = led%rel_residual(ICE_LEDGER_HEAT)
         if (.not. (ieee_is_finite(rw) .and. ieee_is_finite(rs) .and. ieee_is_finite(rh))) then
            call test_failed(error, "non-finite ledger residual")
            exit
         end if
         worst_w = max(worst_w, rw)
         worst_s = max(worst_s, rs)
         worst_h = max(worst_h, rh)
         bound = 1.1_wp*ICE_FERRMAX*DT*AREA*a0 + 1.0e-15_wp*led%scale(ICE_LEDGER_HEAT)
         worst_ratio = max(worst_ratio, abs(led%res(ICE_LEDGER_HEAT))/bound)
         !$acc update self(fx%meltt)
         if (fx%meltt(IC, JC) > 1.0e-9_wp) nmelt = nmelt + 1
         !$acc update self(st%vsnon)
         if (sum(st%vsnon(IC, JC, :)) > 0.0_wp) nsnow = nsnow + 1
      end do
      !$acc exit data delete(area_a, uatm, vatm, wind, potT, Qa, rhoa, flw, fsnow, uvel, vvel, Tf, sst, frzmlt, tau_x, tau_y) &
      !$acc& delete(sx, sy, fbot, rside, wlat, fswsfcn, fswintn, fswthrun, sxn, syn, Isw, Ssw, ierr_n)
      call fx%exit_data()
      call st%exit_data()
      if (allocated(error)) return
    write (output_unit, '(a,3es11.3,a,es10.3,a,i0,a,i0)') "shipped ledger, 1000 steps, worst relative residual water/salt/heat: ", &
         worst_w, worst_s, worst_h, "; worst heat residual / BL99 bound ", worst_ratio, "; melting steps ", nmelt, &
         ", steps with snow ", nsnow
      call check(error, nmelt > 50 .and. nsnow > 50, "the season exercises melt and snow")
      if (allocated(error)) return
      call check(error, worst_w < 1.0e-9_wp, "water closes (bounded by the round-off of the state total M, ~1e-16 M over the step throughput)")
      if (allocated(error)) return
      call check(error, worst_s < 1.0e-9_wp, "salt closes (same round-off of the state total)")
      if (allocated(error)) return
      call check(error, worst_ratio <= 1.0_wp, "heat closes within the BL99 solver bound")
   end subroutine test_season

end module test_cice_therm1_ledger
