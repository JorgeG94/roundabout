module test_cice_atmo
   !! Tests of `rdb_ice_atmo` (A15): profile-function closed forms and the
   !! documented neutral discontinuity `psi_stable(0) = +0.025`; the neutral
   !! known answer of the whole kernel (exactly zero air-surface differences,
   !! reached with a 0.05 K surface where the saturation humidity underflows to
   !! exactly zero) in closed form for both constant sets and for unequal
   !! momentum and scalar heights; wind floors and the stress frames; Galilean
   !! invariance of the high-frequency stress; clamp saturation at `|z/L| = 10`
   !! (closed form, identical results for twice the temperature difference);
   !! monotone drag through neutral; iteration-count truncation; humidity sign;
   !! NaN-proof error flagging; and a `do concurrent` device test over mapped
   !! arrays.  Expected values are closed forms written here, never the
   !! kernel's own output.
   use, intrinsic :: iso_fortran_env, only: int32, int64, output_unit
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan, ieee_positive_inf
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use rdb_constants, only: wp
   use rdb_ice_errors, only: ICE_ERR_NONFINITE
   use rdb_ice_params, only: ice_constants_t, ice_constants, ice_atm_params_t, ice_params_t, ice_params_default, &
                             ICE_CONSTANTS_ICEPACK, ICE_CONSTANTS_CESM
   use rdb_ice_thermo_shared, only: ice_qsat_surface
   use rdb_ice_atmo, only: ice_psi_m, ice_psi_h, ice_psi_stable, ice_atm_boundary_ice
   implicit none
   private

   public :: collect_cice_atmo_tests

   real(wp), parameter :: LN_20000 = 9.903487552536127_wp   ! ln(zref/iceruf), zref = 10, iceruf = 5e-4
   real(wp), parameter :: LN_0P2 = -1.6094379124341003_wp   ! ln(2/10)
   real(wp), parameter :: LN_1P5 = 0.4054651081081644_wp    ! ln(15/10)
   real(wp), parameter :: T_SINGLE = 0.05_wp - 273.15_wp    ! 0.05 K surface: qsat underflows to exactly 0

contains

   subroutine collect_cice_atmo_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("profile_functions_closed_form", test_psi), &
                  new_unittest("neutral_discontinuity_pinned", test_neutral_jump), &
                  new_unittest("neutral_known_answer_both_sets", test_neutral), &
                  new_unittest("neutral_unequal_heights", test_neutral_heights), &
                  new_unittest("wind_floor_and_stress_frames", test_frames), &
                  new_unittest("highfreq_galilean_invariance", test_galilean), &
                  new_unittest("stability_clamps", test_clamps), &
                  new_unittest("drag_monotone_through_neutral", test_monotone), &
                  new_unittest("iteration_truncation", test_iteration), &
                  new_unittest("humidity_sign", test_humidity), &
                  new_unittest("nan_and_range_flags", test_nan), &
                  new_unittest("device_kernel_matches_host", test_device), &
                  new_unittest("driver_sweep_over_cell_categories", test_driver) &
                  ]
   end subroutine collect_cice_atmo_tests

   subroutine run(Tsf, potT, uatm, vatm, uvel, vvel, zlvl, zlvs, Qa, rhoa, c, atm, out, ierr)
      !! One call; `out = [shcoef, lhcoef, strx, stry, delt, delq]`.
      real(wp), intent(in) :: Tsf, potT, uatm, vatm, uvel, vvel, zlvl, zlvs, Qa, rhoa
      type(ice_constants_t), intent(in) :: c
      type(ice_atm_params_t), intent(in) :: atm
      real(wp), intent(out) :: out(6)
      integer(int32), intent(inout) :: ierr

      call ice_atm_boundary_ice(Tsf, potT, uatm, vatm, uvel, vvel, hypot(uatm, vatm), zlvl, zlvs, Qa, rhoa, &
                                c, atm, out(1), out(2), out(3), out(4), out(5), out(6), ierr)
   end subroutine run

   logical function close_to(x, y, rtol, atol)
      !! NaN-proof comparison: a non-finite `x` is never close.
      real(wp), intent(in) :: x, y, rtol, atol

      close_to = .false.
      if (.not. ieee_is_finite(x)) return
      close_to = abs(x - y) <= atol + rtol*abs(y)
   end function close_to

   function rdn_of(c, atm) result(r)
      type(ice_constants_t), intent(in) :: c
      type(ice_atm_params_t), intent(in) :: atm
      real(wp) :: r
      r = c%vonkar/LN_20000   ! default iceruf = 5e-4 m (atm unused: closed forms are for the default)
   end function rdn_of

   subroutine test_psi(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: H(4) = [-0.1_wp, -1.0_wp, -5.0_wp, -10.0_wp]
      real(wp), parameter :: PM(4) = [0.283613711212781_wp, 1.1162322497683264_wp, 2.0684370555523786_wp, &
                                      2.549267894070169_wp]
      real(wp), parameter :: PH(4) = [0.5342837819484254_wp, 1.881227284214417_wp, 3.2188758248682006_wp, &
                                      3.8468290966691465_wp]
      real(wp), parameter :: HS(3) = [0.5_wp, 2.0_wp, 10.0_wp]
      real(wp), parameter :: PS(3) = [-2.3616198350387023_wp, -7.519000572524245_wp, -17.602613438463024_wp]
      integer :: i

      do i = 1, 4
         call check(error, close_to(ice_psi_m(H(i)), PM(i), 1.0e-13_wp, 0.0_wp), "psi_m unstable closed form")
         if (allocated(error)) return
         call check(error, close_to(ice_psi_h(H(i)), PH(i), 1.0e-13_wp, 0.0_wp), "psi_h unstable closed form")
         if (allocated(error)) return
      end do
      do i = 1, 3
         call check(error, close_to(ice_psi_m(HS(i)), PS(i), 1.0e-13_wp, 0.0_wp), "psi_m stable closed form")
         if (allocated(error)) return
         call check(error, close_to(ice_psi_h(HS(i)), PS(i), 1.0e-13_wp, 0.0_wp), "psi_h stable closed form")
         if (allocated(error)) return
         call check(error, close_to(ice_psi_stable(HS(i)), PS(i), 1.0e-13_wp, 0.0_wp), "psi_stable closed form")
         if (allocated(error)) return
      end do
      ! Stable functions are the same for momentum and scalar; unstable ones are not.
      call check(error, ice_psi_m(-1.0_wp) /= ice_psi_h(-1.0_wp), "unstable functions differ")
   end subroutine test_psi

   subroutine test_neutral_jump(error)
      !! psi_stable(0) = +0.025 (not 0) and the stable branch owns exactly +0, while the
      !! unstable branch tends to 0: the upstream discontinuity, pinned.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: negzero

      negzero = sign(0.0_wp, -1.0_wp)
      call check(error, close_to(ice_psi_stable(0.0_wp), 0.025_wp, 0.0_wp, 1.0e-12_wp), "psi_stable(0) = 0.025")
      if (allocated(error)) return
      call check(error, close_to(ice_psi_m(0.0_wp), 0.025_wp, 0.0_wp, 1.0e-12_wp), "psi_m(+0) = 0.025")
      if (allocated(error)) return
      call check(error, close_to(ice_psi_h(0.0_wp), 0.025_wp, 0.0_wp, 1.0e-12_wp), "psi_h(+0) = 0.025")
      if (allocated(error)) return
      call check(error, close_to(ice_psi_m(1.0e-9_wp), 0.025_wp, 0.0_wp, 1.0e-8_wp), "psi_m(0+) -> 0.025")
      if (allocated(error)) return
      call check(error, close_to(ice_psi_m(-1.0e-9_wp), 0.0_wp, 0.0_wp, 1.0e-8_wp), "psi_m(0-) -> 0")
      if (allocated(error)) return
      call check(error, close_to(ice_psi_h(-1.0e-9_wp), 0.0_wp, 0.0_wp, 1.0e-8_wp), "psi_h(0-) -> 0")
      if (allocated(error)) return
      call check(error, close_to(ice_psi_m(negzero), 0.0_wp, 0.0_wp, 1.0e-15_wp), "psi_m(-0) is the unstable branch")
   end subroutine test_neutral_jump

   subroutine test_neutral(error)
      !! delt = delq = 0 exactly (TsfK = potT, qsat underflows to 0 at 0.05 K, Qa = 0), so z/L = +0 and the
      !! stable branch gives psi = 0.025 at every pass:
      !! rd = rdn/(1 + rdn/kappa (alzm - 0.025)), ustar = rd vmag,
      !! shcoef = rhoa ustar cp_air rd + 1, lhcoef = rhoa ustar Lsub rd, strx = rhoa ustar rd uatm.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      type(ice_constants_t) :: c
      integer :: iset, ih
      real(wp) :: out(6), rdn, rd, ustar, vmag, potT, rho, uatm, vatm, ref(4)
      integer(int32) :: ierr

      rho = 1.3_wp
      do iset = ICE_CONSTANTS_ICEPACK, ICE_CONSTANTS_CESM
         p = ice_params_default(iset, 5, 4, 1)
         c = p%c
         rdn = rdn_of(c, p%atm)
         call check(error, close_to(rdn*rdn, 1.631e-3_wp, 1.0e-3_wp, 0.0_wp), "rdn^2 = 1.631e-3 to 4 digits")
         if (allocated(error)) return
         rd = rdn/(1.0_wp + rdn/c%vonkar*(0.0_wp - 0.025_wp))
         call check(error, close_to(rd/rdn, 1.0025_wp, 5.0e-5_wp, 0.0_wp), "neutral drag 0.5 percent above rdn^2")
         if (allocated(error)) return
         potT = T_SINGLE + c%Tffresh
         do ih = 0, 1
            p%atm%highfreq = (ih == 1)
            uatm = 4.0_wp
            vatm = -3.0_wp
            ierr = 0_int32
            if (ih == 0) then
               vmag = 5.0_wp
               ustar = rd*vmag
               ref = [rho*ustar*c%cp_air*rd + 1.0_wp, rho*ustar*c%Lsub*rd, rho*ustar*rd*uatm, rho*ustar*rd*vatm]
               call run(T_SINGLE, potT, uatm, vatm, 0.0_wp, 0.0_wp, 10.0_wp, 10.0_wp, 0.0_wp, rho, c, p%atm, out, ierr)
            else
               ! ice at rest: relative wind = wind
               vmag = 5.0_wp
               ustar = rd*vmag
               ref = [rho*ustar*c%cp_air*rd + 1.0_wp, rho*ustar*c%Lsub*rd, rho*rd*rd*5.0_wp*uatm, rho*rd*rd*5.0_wp*vatm]
               call run(T_SINGLE, potT, uatm, vatm, 0.0_wp, 0.0_wp, 10.0_wp, 10.0_wp, 0.0_wp, rho, c, p%atm, out, ierr)
            end if
            call check(error, ierr == 0_int32, "no error bits")
            if (allocated(error)) return
            call check(error, out(5) == 0.0_wp .and. out(6) == 0.0_wp, "delt and delq exactly zero")
            if (allocated(error)) return
            call check(error, close_to(out(1), ref(1), 1.0e-13_wp, 0.0_wp), "shcoef neutral closed form")
            if (allocated(error)) return
            call check(error, close_to(out(2), ref(2), 1.0e-13_wp, 0.0_wp), "lhcoef neutral closed form")
            if (allocated(error)) return
            call check(error, close_to(out(3), ref(3), 1.0e-13_wp, 0.0_wp), "strx neutral closed form")
            if (allocated(error)) return
            call check(error, close_to(out(4), ref(4), 1.0e-13_wp, 0.0_wp), "stry neutral closed form")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_neutral

   subroutine test_neutral_heights(error)
      !! zlvl = 2 shifts alzm to ln(0.2); zlvs = 15 with zlvl = 10 gives rh /= rd.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      type(ice_constants_t) :: c
      real(wp) :: out(6), rdn, rd, rh, re, ustar, potT, rho
      integer(int32) :: ierr

      p = ice_params_default(ICE_CONSTANTS_ICEPACK, 5, 4, 1)
      c = p%c
      rdn = rdn_of(c, p%atm)
      rho = 1.3_wp
      potT = T_SINGLE + c%Tffresh
      ierr = 0_int32
      call run(T_SINGLE, potT, 4.0_wp, -3.0_wp, 0.0_wp, 0.0_wp, 2.0_wp, 2.0_wp, 0.0_wp, rho, c, p%atm, out, ierr)
      rd = rdn/(1.0_wp + rdn/c%vonkar*(LN_0P2 - 0.025_wp))
      ustar = rd*5.0_wp
      call check(error, close_to(out(1), rho*ustar*c%cp_air*rd + 1.0_wp, 1.0e-13_wp, 0.0_wp), "shcoef at z = 2 m")
      if (allocated(error)) return
      call check(error, close_to(out(3), rho*ustar*rd*4.0_wp, 1.0e-13_wp, 0.0_wp), "strx at z = 2 m")
      if (allocated(error)) return

      call run(T_SINGLE, potT, 4.0_wp, -3.0_wp, 0.0_wp, 0.0_wp, 10.0_wp, 15.0_wp, 0.0_wp, rho, c, p%atm, out, ierr)
      rd = rdn/(1.0_wp + rdn/c%vonkar*(0.0_wp - 0.025_wp))
      rh = rdn/(1.0_wp + rdn/c%vonkar*(LN_1P5 - 0.025_wp))
      re = rh
      ustar = rd*5.0_wp
      call check(error, close_to(out(1), rho*ustar*c%cp_air*rh + 1.0_wp, 1.0e-13_wp, 0.0_wp), "shcoef with zlvs /= zlvl")
      if (allocated(error)) return
      call check(error, close_to(out(2), rho*ustar*c%Lsub*re, 1.0e-13_wp, 0.0_wp), "lhcoef with zlvs /= zlvl")
      if (allocated(error)) return
      call check(error, close_to(out(3), rho*ustar*rd*4.0_wp, 1.0e-13_wp, 0.0_wp), "stress uses the momentum height only")
   end subroutine test_neutral_heights

   subroutine test_frames(error)
      !! Wind floor (1 m/s; 0.5 m/s with highfreq acting on the relative speed only), stress
      !! directions, zero relative velocity, ice moving with the wind.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      type(ice_constants_t) :: c
      real(wp) :: out(6), o2(6), rdn, rd, potT, rho
      integer(int32) :: ierr

      p = ice_params_default(ICE_CONSTANTS_ICEPACK, 5, 4, 1)
      c = p%c
      rdn = rdn_of(c, p%atm)
      rd = rdn/(1.0_wp + rdn/c%vonkar*(0.0_wp - 0.025_wp))
      rho = 1.3_wp
      potT = T_SINGLE + c%Tffresh
      ierr = 0_int32

      ! weak wind 0.3 m/s without highfreq: vmag = 1; stress = rho ustar rd uatm with ustar = rd * 1.
      call run(T_SINGLE, potT, 0.3_wp, 0.0_wp, 0.0_wp, 0.0_wp, 10.0_wp, 10.0_wp, 0.0_wp, rho, c, p%atm, out, ierr)
      call check(error, close_to(out(3), rho*rd*rd*1.0_wp*0.3_wp, 1.0e-13_wp, 0.0_wp), "floored vmag in the stress")
      if (allocated(error)) return
      call check(error, close_to(out(1), rho*(rd*1.0_wp)*c%cp_air*rd + 1.0_wp, 1.0e-13_wp, 0.0_wp), "floored vmag in shcoef")
      if (allocated(error)) return
      call check(error, out(4) == 0.0_wp, "stry zero for vatm = 0")
      if (allocated(error)) return

      ! highfreq: stress is parallel to the relative wind; coefficients use the 0.5 m/s floor.
      p%atm%highfreq = .true.
      call run(T_SINGLE, potT, 0.3_wp, 0.0_wp, 0.0_wp, 0.0_wp, 10.0_wp, 10.0_wp, 0.0_wp, rho, c, p%atm, out, ierr)
      call check(error, close_to(out(3), rho*rd*rd*0.3_wp*0.3_wp, 1.0e-13_wp, 0.0_wp), "stress from the unfloored |dU|")
      if (allocated(error)) return
      call check(error, close_to(out(1), rho*(rd*0.5_wp)*c%cp_air*rd + 1.0_wp, 1.0e-13_wp, 0.0_wp), "0.5 m/s floor in shcoef")
      if (allocated(error)) return
      call run(T_SINGLE, potT, 4.0_wp, -3.0_wp, 1.0_wp, 2.0_wp, 10.0_wp, 10.0_wp, 0.0_wp, rho, c, p%atm, out, ierr)
      call check(error, close_to(out(3)*(-5.0_wp), out(4)*3.0_wp, 1.0e-13_wp, 1.0e-30_wp), "stress parallel to Uatm - Uice")
      if (allocated(error)) return
      call check(error, out(3) > 0.0_wp .and. out(4) < 0.0_wp, "signs of stress follow dU = (3, -5)")
      if (allocated(error)) return

      ! zero relative velocity: no stress, coefficients finite and set by the 0.5 m/s floor.
      call run(T_SINGLE, potT, 4.0_wp, -3.0_wp, 4.0_wp, -3.0_wp, 10.0_wp, 10.0_wp, 0.0_wp, rho, c, p%atm, out, ierr)
      call check(error, out(3) == 0.0_wp .and. out(4) == 0.0_wp, "zero relative velocity: zero stress")
      if (allocated(error)) return
      call check(error, ieee_is_finite(out(1)) .and. ieee_is_finite(out(2)) .and. out(1) > 1.0_wp, "finite coefficients")
      if (allocated(error)) return
      call check(error, close_to(out(1), rho*(rd*0.5_wp)*c%cp_air*rd + 1.0_wp, 1.0e-13_wp, 0.0_wp), "floor sets shcoef")
      if (allocated(error)) return

      ! ice moving with the wind reduces the stress against ice at rest.
      call run(T_SINGLE, potT, 4.0_wp, -3.0_wp, 3.0_wp, -2.0_wp, 10.0_wp, 10.0_wp, 0.0_wp, rho, c, p%atm, out, ierr)
      call run(T_SINGLE, potT, 4.0_wp, -3.0_wp, 0.0_wp, 0.0_wp, 10.0_wp, 10.0_wp, 0.0_wp, rho, c, p%atm, o2, ierr)
      call check(error, hypot(out(3), out(4)) < 0.2_wp*hypot(o2(3), o2(4)), "moving ice feels less stress")
   end subroutine test_frames

   subroutine test_galilean(error)
      !! Adding one constant to the wind and the ice velocity leaves the stress and the coefficients
      !! unchanged (highfreq), over stable, unstable and near-neutral surfaces.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      real(wp) :: o1(6), o2(6), s, mag
      integer :: i, j
      integer(int32) :: ierr
      real(wp), parameter :: TS(5) = [-20.0_wp, -2.0_wp, -10.0_wp, 0.0_wp, -35.0_wp]
      real(wp), parameter :: DT(5) = [12.0_wp, -30.0_wp, 0.3_wp, 6.0_wp, 25.0_wp]
      real(wp), parameter :: SH(3) = [7.3_wp, -2.1_wp, 21.0_wp]
      real(wp), parameter :: SV(3) = [-2.1_wp, 11.0_wp, 5.5_wp]

      p = ice_params_default(ICE_CONSTANTS_ICEPACK, 5, 4, 1)
      p%atm%highfreq = .true.
      ierr = 0_int32
      do i = 1, 5
         do j = 1, 3
            s = TS(i) + p%c%Tffresh
            call run(TS(i), s + DT(i), 4.0_wp, -3.0_wp, 0.4_wp, 0.2_wp, 10.0_wp, 10.0_wp, 5.0e-4_wp, 1.3_wp, &
                     p%c, p%atm, o1, ierr)
            call run(TS(i), s + DT(i), 4.0_wp + SH(j), -3.0_wp + SV(j), 0.4_wp + SH(j), 0.2_wp + SV(j), 10.0_wp, &
                     10.0_wp, 5.0e-4_wp, 1.3_wp, p%c, p%atm, o2, ierr)
            mag = maxval(abs(o1(1:4)))
            call check(error, close_to(o2(1), o1(1), 1.0e-12_wp, 1.0e-13_wp), "shcoef invariant")
            if (allocated(error)) return
            call check(error, close_to(o2(2), o1(2), 1.0e-12_wp, 1.0e-13_wp), "lhcoef invariant")
            if (allocated(error)) return
            call check(error, close_to(o2(3), o1(3), 1.0e-12_wp, 1.0e-13_wp*mag), "strx invariant")
            if (allocated(error)) return
            call check(error, close_to(o2(4), o1(4), 1.0e-12_wp, 1.0e-13_wp*mag), "stry invariant")
            if (allocated(error)) return
         end do
      end do
      ! The non-highfreq stress is NOT invariant (it uses the absolute wind): pin that.
      p%atm%highfreq = .false.
      s = -10.0_wp + p%c%Tffresh
      call run(-10.0_wp, s + 0.3_wp, 4.0_wp, -3.0_wp, 0.0_wp, 0.0_wp, 10.0_wp, 10.0_wp, 5.0e-4_wp, 1.3_wp, p%c, p%atm, o1, ierr)
      call run(-10.0_wp, s + 0.3_wp, 11.0_wp, -3.0_wp, 7.0_wp, 0.0_wp, 10.0_wp, 10.0_wp, 5.0e-4_wp, 1.3_wp, p%c, p%atm, o2, ierr)
      call check(error, abs(o2(3)) > 2.0_wp*abs(o1(3)), "no Galilean invariance without highfreq")
   end subroutine test_galilean

   subroutine test_clamps(error)
      !! |z/L| saturates at 10: the stress is the closed form with psi(+-10) and is identical for
      !! twice the temperature difference.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      type(ice_constants_t) :: c
      real(wp) :: o30(6), o60(6), rdn, rd, ts, rho
      real(wp), parameter :: PSI_S10 = -17.602613438463024_wp
      real(wp), parameter :: PSI_M10 = 2.549267894070169_wp
      integer(int32) :: ierr

      p = ice_params_default(ICE_CONSTANTS_ICEPACK, 5, 4, 1)
      c = p%c
      rdn = rdn_of(c, p%atm)
      rho = 1.3_wp
      ierr = 0_int32
      ts = -20.0_wp + c%Tffresh

      call run(-20.0_wp, ts + 30.0_wp, 2.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 10.0_wp, 10.0_wp, 5.0e-4_wp, rho, c, p%atm, o30, ierr)
      call run(-20.0_wp, ts + 60.0_wp, 2.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 10.0_wp, 10.0_wp, 5.0e-4_wp, rho, c, p%atm, o60, ierr)
      rd = rdn/(1.0_wp + rdn/c%vonkar*(0.0_wp - PSI_S10))
      call check(error, close_to(o30(3), rho*rd*rd*2.0_wp*2.0_wp, 1.0e-12_wp, 0.0_wp), "stable clamp closed form")
      if (allocated(error)) return
      call check(error, close_to(rd*rd/(rdn*rdn), 0.13_wp, 5.0e-2_wp, 0.0_wp), "Cd/Cdn about 0.13 at the stable clamp")
      if (allocated(error)) return
      call check(error, close_to(o60(3), o30(3), 1.0e-13_wp, 0.0_wp), "stable stress identical for delt 30 and 60")
      if (allocated(error)) return
      call check(error, close_to(o60(1), o30(1), 1.0e-13_wp, 0.0_wp), "stable shcoef identical for delt 30 and 60")
      if (allocated(error)) return
      call check(error, o60(5) == 60.0_wp .or. close_to(o60(5), 60.0_wp, 1.0e-12_wp, 0.0_wp), "delt reports the true difference")
      if (allocated(error)) return

      call run(-20.0_wp, ts - 30.0_wp, 2.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 10.0_wp, 10.0_wp, 5.0e-4_wp, rho, c, p%atm, o30, ierr)
      call run(-20.0_wp, ts - 60.0_wp, 2.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 10.0_wp, 10.0_wp, 5.0e-4_wp, rho, c, p%atm, o60, ierr)
      rd = rdn/(1.0_wp + rdn/c%vonkar*(0.0_wp - PSI_M10))
      call check(error, close_to(o30(3), rho*rd*rd*2.0_wp*2.0_wp, 1.0e-12_wp, 0.0_wp), "unstable clamp closed form")
      if (allocated(error)) return
      call check(error, rd*rd/(rdn*rdn) > 1.0_wp, "Cd/Cdn > 1 when unstable")
      if (allocated(error)) return
      call check(error, close_to(o60(3), o30(3), 1.0e-13_wp, 0.0_wp), "unstable stress identical for delt -30 and -60")
      if (allocated(error)) return
      call check(error, close_to(o60(1), o30(1), 1.0e-13_wp, 0.0_wp), "unstable shcoef identical for delt -30 and -60")
      if (allocated(error)) return
      call check(error, ierr == 0_int32, "no error bits")
   end subroutine test_clamps

   subroutine test_monotone(error)
      !! The drag coefficient falls with delt from strongly unstable to strongly stable (dry
      !! surface so delq is negligible); the exact-zero 0.5 percent jump is the stable/unstable boundary.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      real(wp), parameter :: DT(8) = [-30.0_wp, -10.0_wp, -2.0_wp, -0.5_wp, 0.5_wp, 1.0_wp, 2.0_wp, 4.0_wp]
      real(wp) :: out(6), cd(8), ts
      integer :: i
      integer(int32) :: ierr

      p = ice_params_default(ICE_CONSTANTS_ICEPACK, 5, 4, 1)
      ierr = 0_int32
      ts = -80.0_wp + p%c%Tffresh
      do i = 1, 8
         call run(-80.0_wp, ts + DT(i), 3.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 10.0_wp, 10.0_wp, 0.0_wp, 1.3_wp, &
                  p%c, p%atm, out, ierr)
         cd(i) = out(3)/(1.3_wp*3.0_wp*3.0_wp)
         call check(error, ieee_is_finite(cd(i)), "finite drag")
         if (allocated(error)) return
      end do
      do i = 1, 7
         call check(error, cd(i + 1) < cd(i), "drag decreases with delt")
         if (allocated(error)) return
      end do
      call check(error, ierr == 0_int32, "no error bits")
   end subroutine test_monotone

   subroutine test_iteration(error)
      !! 5 passes against 50: unstable and weakly stable cases agree to 1e-6 relative; the
      !! moderately stable case (delt = +3 K) is NOT converged after 5 passes (measured 4e-4
      !! relative in the stress, an upstream property) and is bounded between 1e-5 and 1e-2; a single
      !! pass differs from the converged answer; natmiter = 0 keeps the neutral coefficients.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p5, p50, p1
      real(wp), parameter :: DT(4) = [3.0_wp, -3.0_wp, 0.5_wp, -0.5_wp]
      real(wp) :: rel
      real(wp) :: o5(6), o50(6), o1(6), ts
      integer :: i
      integer(int32) :: ierr

      p5 = ice_params_default(ICE_CONSTANTS_ICEPACK, 5, 4, 1)
      p50 = p5
      p50%atm%natmiter = 50
      p1 = p5
      p1%atm%natmiter = 1
      ierr = 0_int32
      ts = -10.0_wp + p5%c%Tffresh
      do i = 1, 4
         call run(-10.0_wp, ts + DT(i), 4.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 10.0_wp, 10.0_wp, 5.0e-4_wp, 1.3_wp, &
                  p5%c, p5%atm, o5, ierr)
         call run(-10.0_wp, ts + DT(i), 4.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 10.0_wp, 10.0_wp, 5.0e-4_wp, 1.3_wp, &
                  p50%c, p50%atm, o50, ierr)
         call run(-10.0_wp, ts + DT(i), 4.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 10.0_wp, 10.0_wp, 5.0e-4_wp, 1.3_wp, &
                  p1%c, p1%atm, o1, ierr)
         if (i == 1) then
            rel = abs(o5(3) - o50(3))/abs(o50(3))
            call check(error, ieee_is_finite(rel) .and. rel > 1.0e-5_wp .and. rel < 1.0e-2_wp, &
                       "delt = +3: 5 passes not converged, bounded")
            if (allocated(error)) return
         else
            call check(error, close_to(o5(3), o50(3), 1.0e-6_wp, 0.0_wp), "stress: 5 passes = 50 passes to 1e-6")
            if (allocated(error)) return
            call check(error, close_to(o5(1), o50(1), 1.0e-6_wp, 0.0_wp), "shcoef: 5 passes = 50 passes to 1e-6")
            if (allocated(error)) return
         end if
         call check(error, ieee_is_finite(o1(3)), "one pass is finite")
         if (allocated(error)) return
      end do
      call run(-10.0_wp, ts + 3.0_wp, 4.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 10.0_wp, 10.0_wp, 5.0e-4_wp, 1.3_wp, &
               p1%c, p1%atm, o1, ierr)
      call run(-10.0_wp, ts + 3.0_wp, 4.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 10.0_wp, 10.0_wp, 5.0e-4_wp, 1.3_wp, &
               p50%c, p50%atm, o50, ierr)
      call check(error, abs(o1(3) - o50(3)) > 1.0e-4_wp*abs(o50(3)), "single pass is not converged")
      if (allocated(error)) return
      ! natmiter = 0 keeps the neutral coefficients (and is finite).
      p1%atm%natmiter = 0
      call run(-10.0_wp, ts + 3.0_wp, 4.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 10.0_wp, 10.0_wp, 5.0e-4_wp, 1.3_wp, &
               p1%c, p1%atm, o1, ierr)
      call check(error, ieee_is_finite(o1(1)) .and. ieee_is_finite(o1(3)), "natmiter = 0 finite")
      if (allocated(error)) return
      call check(error, ierr == 0_int32, "no error bits")
   end subroutine test_iteration

   subroutine test_humidity(error)
      !! Warm moist air over cold ice: delq > 0 and lhcoef > 0 (condensation, downward positive);
      !! cold dry air: delq < 0.  delq = Qa - qsat(TsfK)/rhoa to the reciprocal-versus-divide ulps.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      real(wp) :: out(6), ts, qs
      integer(int32) :: ierr

      p = ice_params_default(ICE_CONSTANTS_ICEPACK, 5, 4, 1)
      ierr = 0_int32
      ts = -10.0_wp + p%c%Tffresh
      qs = ice_qsat_surface(ts, 1.3_wp)
      call run(-10.0_wp, ts + 5.0_wp, 4.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 10.0_wp, 10.0_wp, 1.5_wp*qs, 1.3_wp, &
               p%c, p%atm, out, ierr)
      call check(error, out(6) > 0.0_wp .and. out(2) > 0.0_wp, "moist air: condensation")
      if (allocated(error)) return
      call check(error, close_to(out(6), 0.5_wp*qs, 1.0e-12_wp, 0.0_wp), "delq = Qa - Qsfc")
      if (allocated(error)) return
      call run(-10.0_wp, ts + 5.0_wp, 4.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 10.0_wp, 10.0_wp, 0.2_wp*qs, 1.3_wp, &
               p%c, p%atm, out, ierr)
      call check(error, out(6) < 0.0_wp, "dry air: sublimation, delq < 0")
      if (allocated(error)) return
      call check(error, close_to(out(5), 5.0_wp, 1.0e-12_wp, 0.0_wp), "delt = potT - TsfK")
   end subroutine test_humidity

   subroutine test_nan(error)
      !! A non-finite or non-physical input flags ICE_ERR_NONFINITE, ORs into the word, and returns
      !! every output as exactly zero; an input the branch does not read is not an error.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      real(wp) :: out(6), nan, inf, ts
      integer(int32) :: ierr
      integer :: i
      logical :: flagged(12)

      nan = ieee_value(1.0_wp, ieee_quiet_nan)
      inf = ieee_value(1.0_wp, ieee_positive_inf)
      p = ice_params_default(ICE_CONSTANTS_ICEPACK, 5, 4, 1)
      ts = -10.0_wp + p%c%Tffresh
      flagged = .false.
      do i = 1, 12
         out = 7.0_wp
         ierr = 0_int32
         select case (i)
         case (1); call run(nan, ts, 4.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 10.0_wp, 10.0_wp, 1.0e-3_wp, 1.3_wp, p%c, p%atm, out, ierr)
         case (2); call run(-10.0_wp, nan, 4.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 10.0_wp, 10.0_wp, 1.0e-3_wp, 1.3_wp, p%c, p%atm, out, ierr)
         case (3); call run(-10.0_wp, ts, nan, 1.0_wp, 0.0_wp, 0.0_wp, 10.0_wp, 10.0_wp, 1.0e-3_wp, 1.3_wp, p%c, p%atm, out, ierr)
         case (4); call run(-10.0_wp, ts, 4.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 10.0_wp, 10.0_wp, nan, 1.3_wp, p%c, p%atm, out, ierr)
         case (5); call run(-10.0_wp, ts, 4.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 10.0_wp, 10.0_wp, 1.0e-3_wp, inf, p%c, p%atm, out, ierr)
         case (6); call run(-10.0_wp, ts, 4.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, nan, 10.0_wp, 1.0e-3_wp, 1.3_wp, p%c, p%atm, out, ierr)
         case (7); call run(-300.0_wp, ts, 4.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 10.0_wp, 10.0_wp, 1.0e-3_wp, 1.3_wp, p%c, p%atm, out, ierr)
         case (8); call run(-10.0_wp, ts, 4.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 10.0_wp, 10.0_wp, 1.0e-3_wp, 0.0_wp, p%c, p%atm, out, ierr)
         case (9); call run(-10.0_wp, ts, 4.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 10.0_wp, -1.0_wp, 1.0e-3_wp, 1.3_wp, p%c, p%atm, out, ierr)
         case (10)
            ! the wind speed enters only without highfreq: a NaN ice velocity is fine there
            call run(-10.0_wp, ts, 4.0_wp, 1.0_wp, nan, nan, 10.0_wp, 10.0_wp, 1.0e-3_wp, 1.3_wp, p%c, p%atm, out, ierr)
         case (11)
            p%atm%highfreq = .true.
            call run(-10.0_wp, ts, 4.0_wp, 1.0_wp, nan, 0.0_wp, 10.0_wp, 10.0_wp, 1.0e-3_wp, 1.3_wp, p%c, p%atm, out, ierr)
            p%atm%highfreq = .false.
         case (12)
            ierr = 4_int32
            call run(nan, ts, 4.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 10.0_wp, 10.0_wp, 1.0e-3_wp, 1.3_wp, p%c, p%atm, out, ierr)
            call check(error, iand(ierr, 4_int32) == 4_int32, "existing error bits preserved")
            if (allocated(error)) return
         end select
         flagged(i) = iand(ierr, ICE_ERR_NONFINITE) /= 0_int32
         if (i == 10) then
            call check(error,.not. flagged(i), "unused ice velocity not flagged")
            if (allocated(error)) return
            call check(error, ieee_is_finite(out(3)) .and. out(3) /= 7.0_wp, "result produced")
            if (allocated(error)) return
         else
            call check(error, flagged(i), "bad input flagged")
            if (allocated(error)) return
            call check(error, all(out == 0.0_wp), "outputs defined as zero")
            if (allocated(error)) return
         end if
      end do
   end subroutine test_nan

   subroutine test_device(error)
      !! One `do concurrent` over mapped arrays (constants and options copied in as bundles) matches the
      !! same calls on the host to a tolerance (FMA on the device build).
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      type(ice_constants_t) :: c
      type(ice_atm_params_t) :: atm
      integer, parameter :: N = 2048
      real(wp) :: Tsf(N), potT(N), ua(N), va(N), ui(N), vi(N), wnd(N), Qa(N), rho(N)
      real(wp) :: sh(N), lh(N), sx(N), sy(N), dt(N), dq(N)
      real(wp) :: sh_h(N), lh_h(N), sx_h(N), sy_h(N), dt_h(N), dq_h(N)
      integer(int32) :: ierr(N), ierr_l, ierr_h(N)
      integer :: i
      real(wp) :: r

      p = ice_params_default(ICE_CONSTANTS_CESM, 5, 4, 1)
      p%atm%highfreq = .true.
      c = p%c
      atm = p%atm
      do i = 1, N
         r = real(mod(i*7919, 1000), wp)/1000.0_wp
         Tsf(i) = -35.0_wp + 35.0_wp*r
         potT(i) = Tsf(i) + c%Tffresh + (-20.0_wp + 40.0_wp*real(mod(i*104729, 1000), wp)/1000.0_wp)
         ua(i) = -10.0_wp + 20.0_wp*real(mod(i*31, 997), wp)/997.0_wp
         va(i) = -10.0_wp + 20.0_wp*real(mod(i*17, 991), wp)/991.0_wp
         ui(i) = 0.5_wp*ua(i)/10.0_wp
         vi(i) = -0.3_wp*va(i)/10.0_wp
         wnd(i) = hypot(ua(i), va(i))
         Qa(i) = 2.0e-3_wp*real(mod(i*13, 101), wp)/101.0_wp
         rho(i) = 1.2_wp + 0.2_wp*r
      end do
      ierr = 0_int32
      ierr_h = 0_int32
      do i = 1, N
         ierr_l = 0_int32
         call ice_atm_boundary_ice(Tsf(i), potT(i), ua(i), va(i), ui(i), vi(i), wnd(i), 10.0_wp, 10.0_wp, Qa(i), &
                                   rho(i), c, atm, sh_h(i), lh_h(i), sx_h(i), sy_h(i), dt_h(i), dq_h(i), ierr_l)
         ierr_h(i) = ierr_l
      end do

      !$acc enter data copyin(c, atm, Tsf, potT, ua, va, ui, vi, wnd, Qa, rho, ierr) create(sh, lh, sx, sy, dt, dq)
      do concurrent(i=1:N) local(ierr_l)
         ierr_l = 0_int32
         call ice_atm_boundary_ice(Tsf(i), potT(i), ua(i), va(i), ui(i), vi(i), wnd(i), 10.0_wp, 10.0_wp, Qa(i), &
                                   rho(i), c, atm, sh(i), lh(i), sx(i), sy(i), dt(i), dq(i), ierr_l)
         ierr(i) = ierr_l
      end do
      !$acc update self(sh, lh, sx, sy, dt, dq, ierr)
      !$acc exit data delete(c, atm, Tsf, potT, ua, va, ui, vi, wnd, Qa, rho, ierr, sh, lh, sx, sy, dt, dq)

      call check(error, all(ierr == 0_int32) .and. all(ierr_h == 0_int32), "no error bits")
      if (allocated(error)) return
      do i = 1, N
         if (.not. (ieee_is_finite(sh(i)) .and. ieee_is_finite(lh(i)) .and. ieee_is_finite(sx(i)) &
                    .and. ieee_is_finite(sy(i)) .and. ieee_is_finite(dt(i)) .and. ieee_is_finite(dq(i)))) then
            call test_failed(error, "non-finite device output")
            return
         end if
         if (.not. (close_to(sh(i), sh_h(i), 1.0e-11_wp, 1.0e-10_wp) .and. close_to(lh(i), lh_h(i), 1.0e-11_wp, 1.0e-8_wp) &
                    .and. close_to(sx(i), sx_h(i), 1.0e-11_wp, 1.0e-12_wp) &
                    .and. close_to(sy(i), sy_h(i), 1.0e-11_wp, 1.0e-12_wp) &
                    .and. close_to(dt(i), dt_h(i), 1.0e-12_wp, 1.0e-10_wp) &
                    .and. close_to(dq(i), dq_h(i), 1.0e-11_wp, 1.0e-15_wp))) then
            call test_failed(error, "device result differs from host")
            return
         end if
      end do
   end subroutine test_device

   subroutine sweep_impl(nx, ny, ncat, is, ie, js, je, c, atm, Tsf, potT, ua, va, ui, vi, wnd, Qa, rho, &
                         sh, lh, sx, sy, dt, dq, ierr)
      !! Driver loop: one `do concurrent` over (category, j, i); explicit-shape arrays, scalars into the
      !! per-surface kernel.  Every array is device-resident under `mem:separate`.
      integer, intent(in) :: nx, ny, ncat, is, ie, js, je
      type(ice_constants_t), intent(in) :: c
      type(ice_atm_params_t), intent(in) :: atm
      real(wp), intent(in) :: Tsf(nx, ny, ncat), potT(nx, ny, ncat), ua(nx, ny), va(nx, ny), ui(nx, ny), vi(nx, ny)
      real(wp), intent(in) :: wnd(nx, ny), Qa(nx, ny), rho(nx, ny)
      real(wp), intent(inout) :: sh(nx, ny, ncat), lh(nx, ny, ncat), sx(nx, ny, ncat), sy(nx, ny, ncat)
      real(wp), intent(inout) :: dt(nx, ny, ncat), dq(nx, ny, ncat)
      integer(int32), intent(inout) :: ierr(nx, ny, ncat)
      integer :: i, j, n
      integer(int32) :: ierr_l

      do concurrent(n=1:ncat, j=js:je, i=is:ie) local(ierr_l)
         ierr_l = 0_int32
         call ice_atm_boundary_ice(Tsf(i, j, n), potT(i, j, n), ua(i, j), va(i, j), ui(i, j), vi(i, j), wnd(i, j), &
                                   10.0_wp, 10.0_wp, Qa(i, j), rho(i, j), c, atm, sh(i, j, n), lh(i, j, n), &
                                   sx(i, j, n), sy(i, j, n), dt(i, j, n), dq(i, j, n), ierr_l)
         ierr(i, j, n) = ierr_l
      end do
   end subroutine sweep_impl

   subroutine test_driver(error)
      !! 512 x 512 x 5 cell-categories through `sweep_impl`; sampled host comparison; reports throughput.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 512, NY = 512, NC = 5, NREP = 5
      type(ice_params_t) :: p
      type(ice_constants_t) :: c
      type(ice_atm_params_t) :: atm
      real(wp), allocatable :: Tsf(:, :, :), potT(:, :, :), sh(:, :, :), lh(:, :, :), sx(:, :, :), sy(:, :, :)
      real(wp), allocatable :: dt(:, :, :), dq(:, :, :), ua(:, :), va(:, :), ui(:, :), vi(:, :), wnd(:, :), Qa(:, :), rho(:, :)
      integer(int32), allocatable :: ierr(:, :, :)
      real(wp) :: o(6), r, secs
      integer(int32) :: ierr_l
      integer(int64) :: t0, t1, rate
      integer :: i, j, n, k

      p = ice_params_default(ICE_CONSTANTS_ICEPACK, NC, 7, 1)
      p%atm%highfreq = .true.
      c = p%c
      atm = p%atm
      allocate (Tsf(NX, NY, NC), potT(NX, NY, NC), sh(NX, NY, NC), lh(NX, NY, NC), sx(NX, NY, NC), sy(NX, NY, NC), &
                dt(NX, NY, NC), dq(NX, NY, NC), ua(NX, NY), va(NX, NY), ui(NX, NY), vi(NX, NY), wnd(NX, NY), &
                Qa(NX, NY), rho(NX, NY), ierr(NX, NY, NC))
      do n = 1, NC
         do j = 1, NY
            do i = 1, NX
               r = real(mod(i*7919 + j*104729 + n*15485863, 1000), wp)/1000.0_wp
               Tsf(i, j, n) = -35.0_wp + 35.0_wp*r
               potT(i, j, n) = Tsf(i, j, n) + c%Tffresh - 20.0_wp + 40.0_wp*real(mod(i*31 + j*17 + n, 997), wp)/997.0_wp
            end do
         end do
      end do
      do j = 1, NY
         do i = 1, NX
            ua(i, j) = -10.0_wp + 20.0_wp*real(mod(i*31 + j, 997), wp)/997.0_wp
            va(i, j) = -10.0_wp + 20.0_wp*real(mod(i + j*17, 991), wp)/991.0_wp
            ui(i, j) = 0.05_wp*ua(i, j)
            vi(i, j) = -0.03_wp*va(i, j)
            wnd(i, j) = hypot(ua(i, j), va(i, j))
            Qa(i, j) = 2.0e-3_wp*real(mod(i*13 + j, 101), wp)/101.0_wp
            rho(i, j) = 1.2_wp + 0.2_wp*real(mod(i + j, 10), wp)/10.0_wp
         end do
      end do
      sh = 0.0_wp
      lh = 0.0_wp
      sx = 0.0_wp
      sy = 0.0_wp
      dt = 0.0_wp
      dq = 0.0_wp
      ierr = 0_int32

      !$acc enter data copyin(c, atm, Tsf, potT, ua, va, ui, vi, wnd, Qa, rho, ierr) create(sh, lh, sx, sy, dt, dq)
      call sweep_impl(NX, NY, NC, 1, NX, 1, NY, c, atm, Tsf, potT, ua, va, ui, vi, wnd, Qa, rho, sh, lh, sx, sy, dt, dq, ierr)
      call system_clock(t0, rate)
      do k = 1, NREP
         call sweep_impl(NX, NY, NC, 1, NX, 1, NY, c, atm, Tsf, potT, ua, va, ui, vi, wnd, Qa, rho, sh, lh, sx, sy, dt, dq, &
                         ierr)
      end do
      call system_clock(t1)
      !$acc update self(sh, lh, sx, sy, dt, dq, ierr)
      !$acc exit data delete(c, atm, Tsf, potT, ua, va, ui, vi, wnd, Qa, rho, ierr, sh, lh, sx, sy, dt, dq)
      secs = real(t1 - t0, wp)/real(rate, wp)
      write (output_unit, '(a,f12.4,a,es11.3,a)') "atmo driver sweep: ", secs/real(NREP, wp)*1.0e3_wp, " ms/sweep, ", &
         real(NX*NY*NC, wp)*real(NREP, wp)/secs, " cell-categories/s"

      call check(error, all(ierr == 0_int32), "no error bits")
      if (allocated(error)) return
      do k = 1, NX*NY*NC, 101
         n = (k - 1)/(NX*NY) + 1
         j = (k - 1 - (n - 1)*NX*NY)/NX + 1
         i = k - (n - 1)*NX*NY - (j - 1)*NX
         ierr_l = 0_int32
         call ice_atm_boundary_ice(Tsf(i, j, n), potT(i, j, n), ua(i, j), va(i, j), ui(i, j), vi(i, j), wnd(i, j), &
                                   10.0_wp, 10.0_wp, Qa(i, j), rho(i, j), c, atm, o(1), o(2), o(3), o(4), o(5), o(6), ierr_l)
         if (.not. (ieee_is_finite(sh(i, j, n)) .and. ieee_is_finite(sx(i, j, n)) .and. ieee_is_finite(dq(i, j, n)))) then
            call test_failed(error, "non-finite driver output")
            return
         end if
         ! Branch switch at |z/L| ~ 1e-20 can differ under FMA; sampled inputs are not that close to neutral.
         if (.not. (close_to(sh(i, j, n), o(1), 1.0e-11_wp, 1.0e-10_wp) .and. close_to(sx(i, j, n), o(3), 1.0e-11_wp, 1.0e-12_wp) &
                    .and. close_to(sy(i, j, n), o(4), 1.0e-11_wp, 1.0e-12_wp) &
                    .and. close_to(lh(i, j, n), o(2), 1.0e-11_wp, 1.0e-8_wp))) then
            call test_failed(error, "driver result differs from host")
            return
         end if
      end do
   end subroutine test_driver

end module test_cice_atmo
