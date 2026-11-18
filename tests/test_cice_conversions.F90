!! Enthalpy <-> temperature conversions of the Bitz and Lipscomb (1999)
!! thermodynamics (analytical test A2) and the shared test helpers.
module test_cice_conversions
   use, intrinsic :: iso_fortran_env, only: int64
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_ice_params, only: ice_constants_t, ice_constants, ICE_CONSTANTS_ICEPACK, ICE_CONSTANTS_CESM, &
                             ice_q_from_T, ice_T_from_q, ice_qsno_from_T, ice_qsno_init_from_T, ice_Tsno_from_q
   use test_cice_util, only: ice_random_columns, ice_ledger_t, ice_golden_first_bad, ICE_LEDGER_MAX
   implicit none
   private

   public :: collect_cice_conversions_tests

contains

   subroutine collect_cice_conversions_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("fresh_ice_closed_form", test_fresh_closed_form), &
                  new_unittest("brine_closed_form", test_brine_closed_form), &
                  new_unittest("melting_point_enthalpy", test_melting_point), &
                  new_unittest("heat_capacity_derivative", test_derivative), &
                  new_unittest("round_trip_random_columns", test_round_trip), &
                  new_unittest("inverse_clips_at_melting_point", test_clip), &
                  new_unittest("inverse_flags_bad_input", test_bad_input), &
                  new_unittest("snow_enthalpy", test_snow), &
                  new_unittest("test_util_ledger_and_golden", test_util) &
                  ]
   end subroutine collect_cice_conversions_tests

   subroutine test_fresh_closed_form(error)
      !! No brine: q = rhoi*(cp_ice*T - Lfresh), T = (Lfresh + q/rhoi)/cp_ice.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      real(wp) :: q, T
      logical :: bad

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      ! Upstream defaults: rhoi 917, cp_ice 2106, Lfresh 3.34e5.
      q = ice_q_from_T(-10.0_wp, 0.0_wp, c, .false.)
      call check(error, abs(q - 917.0_wp*(2106.0_wp*(-10.0_wp) - 3.34e5_wp)) < 1.0e-3_wp, "fresh q(-10)")
      if (allocated(error)) return
      call ice_T_from_q(q, 0.0_wp, c, .false., T, bad)
      call check(error,.not. bad .and. abs(T + 10.0_wp) < 1.0e-12_wp, "fresh T(q)")
      if (allocated(error)) return
      ! Brine on with S = 0 collapses to the same fresh form.
      q = ice_q_from_T(-10.0_wp, 0.0_wp, c, .true.)
      call check(error, abs(q - 917.0_wp*(2106.0_wp*(-10.0_wp) - 3.34e5_wp)) < 1.0e-3_wp, "brine, S=0")
   end subroutine test_fresh_closed_form

   subroutine test_brine_closed_form(error)
      !! S = 5, T = -2: Tmlt = -0.27, then the BL99 expression written out.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      real(wp) :: q, ref, tm

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      tm = -0.27_wp
      ref = -917.0_wp*(2106.0_wp*(tm + 2.0_wp) + 3.34e5_wp*(1.0_wp - tm/(-2.0_wp)) - 4218.0_wp*tm)
      q = ice_q_from_T(-2.0_wp, 5.0_wp, c, .true.)
      call check(error, abs(q - ref) <= 1.0e-12_wp*abs(ref), "brine q(-2, 5 ppt)")
      if (allocated(error)) return
      ! CESM constants change the result through cp_ocn, cp_ice, Lfresh.
      c = ice_constants(ICE_CONSTANTS_CESM)
      ref = -917.0_wp*(2117.27_wp*(tm + 2.0_wp) + 3.337e5_wp*(1.0_wp - tm/(-2.0_wp)) - 3996.0_wp*tm)
      q = ice_q_from_T(-2.0_wp, 5.0_wp, c, .true.)
      call check(error, abs(q - ref) <= 1.0e-12_wp*abs(ref), "brine q, cesm set")
   end subroutine test_brine_closed_form

   subroutine test_melting_point(error)
      !! At T = Tmlt the latent term vanishes: q = rhoi*cp_ocn*Tmlt (the
      !! cp_ocn*Tmlt reference-state offset).
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      real(wp) :: q, tm

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      tm = -0.054_wp*3.2_wp
      q = ice_q_from_T(tm, 3.2_wp, c, .true.)
      call check(error, abs(q - 917.0_wp*4218.0_wp*tm) <= 1.0e-12_wp*abs(q), "q(Tmlt) = rhoi cp_ocn Tmlt")
   end subroutine test_melting_point

   subroutine test_derivative(error)
      !! dq/dT = rhoi*(cp_ice - Lfresh*Tmlt/T**2) (BL99 effective heat capacity),
      !! checked by a central difference.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      real(wp) :: h, fd, exact, tm, T0
      integer :: i
      real(wp), parameter :: TT(3) = [-30.0_wp, -5.0_wp, -1.0_wp]

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      tm = -0.054_wp*5.0_wp
      h = 1.0e-4_wp
      do i = 1, 3
         T0 = TT(i)
         fd = (ice_q_from_T(T0 + h, 5.0_wp, c, .true.) - ice_q_from_T(T0 - h, 5.0_wp, c, .true.))/(2.0_wp*h)
         exact = 917.0_wp*(2106.0_wp - 3.34e5_wp*tm/(T0*T0))
         call check(error, abs(fd - exact) <= 1.0e-6_wp*abs(exact), "dq/dT")
         if (allocated(error)) return
      end do
   end subroutine test_derivative

   subroutine test_round_trip(error)
      !! 1e4 random columns, both constant sets: T -> q -> T and q -> T -> q.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: N = 10000, NL = 4
      real(wp) :: hi(N), hs(N), Tsfc(N), S(NL, N), T(NL, N)
      type(ice_ledger_t) :: ledger
      type(ice_constants_t) :: c
      real(wp) :: q, T2, q2
      logical :: bad
      integer :: set, col, k
      integer, parameter :: SETS(2) = [ICE_CONSTANTS_ICEPACK, ICE_CONSTANTS_CESM]

      call ice_random_columns(N, NL, 20261113_int64, hi, hs, Tsfc, S, T)
      do set = 1, 2
         c = ice_constants(SETS(set))
         do col = 1, N
            do k = 1, NL
               q = ice_q_from_T(T(k, col), S(k, col), c, .true.)
               call ice_T_from_q(q, S(k, col), c, .true., T2, bad)
               if (bad) call ledger%add("bad flag", 1.0_wp, 1.0_wp)
               call ledger%add("T roundtrip", T2 - T(k, col), 1.0_wp)
               q2 = ice_q_from_T(T2, S(k, col), c, .true.)
               call ledger%add("q roundtrip", q2 - q, q)
            end do
         end do
      end do
      call check(error, ledger%worst_of("bad flag") == 0.0_wp, "no discriminant failures")
      if (allocated(error)) return
      ! Measured worst |dT| is 8.5e-14 over 2e6 columns.  The root cancels near
      ! Tmlt (condition number ~400 at S = 12), and the GPU build contracts
      ! multiply-adds, so allow ~10x headroom: 1e-12.
      call check(error, ledger%worst_of("T roundtrip") < 1.0e-12_wp, "T -> q -> T")
      if (allocated(error)) return
      ! q -> T -> q is looser: just below Tmlt (small S) dq/dT = rhoi*(cp_ice -
      ! Lfresh*Tmlt/T**2) reaches ~5e10 J/m3/K, so a 1e-13 degC root error is
      ! ~1e-10 relative in q.
      call check(error, ledger%worst_of("q roundtrip") < 1.0e-9_wp, "q -> T -> q")
   end subroutine test_round_trip

   subroutine test_clip(error)
      !! Enthalpy above the melting enthalpy gives T = Tmlt (the min() clip).
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      real(wp) :: qm, T
      logical :: bad

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      qm = ice_q_from_T(-0.054_wp*4.0_wp, 4.0_wp, c, .true.)
      call ice_T_from_q(qm + 1.0e7_wp, 4.0_wp, c, .true., T, bad)
      call check(error,.not. bad .and. T == -0.054_wp*4.0_wp, "clipped at Tmlt")
   end subroutine test_clip

   subroutine test_bad_input(error)
      !! A negative discriminant (S < 0 is unphysical) and a NaN enthalpy set
      !! the flag instead of aborting.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      real(wp) :: T, q, nan
      logical :: bad

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      ! S = -10: Tmlt = 0.54; choose q so b = 0 and the discriminant is -4 cp Lfresh Tmlt < 0.
      q = 917.0_wp*((c%cp_ocn - c%cp_ice)*0.54_wp - c%Lfresh)
      call ice_T_from_q(q, -10.0_wp, c, .true., T, bad)
      call check(error, bad, "negative discriminant flagged")
      if (allocated(error)) return
      nan = ieee_value(nan, ieee_quiet_nan)
      call ice_T_from_q(nan, 5.0_wp, c, .true., T, bad)
      call check(error, bad, "NaN enthalpy flagged")
   end subroutine test_bad_input

   subroutine test_snow(error)
      !! qsno = -rhos*(Lfresh - cp_ice*min(0,T)); T = (Lfresh + q/rhos)/cp_ice.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      real(wp) :: q

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      q = ice_qsno_from_T(-15.0_wp, c)
      call check(error, abs(q + 330.0_wp*(3.34e5_wp + 2106.0_wp*15.0_wp)) < 1.0e-3_wp, "snow q(-15)")
      if (allocated(error)) return
      call check(error, abs(ice_Tsno_from_q(q, c) + 15.0_wp) < 1.0e-12_wp, "snow T(q)")
      if (allocated(error)) return
      ! Above freezing the enthalpy saturates at the melting value.
      call check(error, ice_qsno_init_from_T(3.0_wp, c) == ice_qsno_init_from_T(0.0_wp, c), &
                 "init form caps at 0 degC")
      if (allocated(error)) return
      call check(error, ice_qsno_init_from_T(-15.0_wp, c) == q, "init form == solver form below 0 degC")
      if (allocated(error)) return
      call check(error, ice_qsno_from_T(3.0_wp, c) > ice_qsno_from_T(0.0_wp, c), "solver form is not clamped")
   end subroutine test_snow

   subroutine test_util(error)
      !! The helpers themselves: determinism, ledger bookkeeping, golden compare.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: N = 50, NL = 3
      real(wp) :: hi(N), hs(N), Tsfc(N), S(NL, N), T(NL, N)
      real(wp) :: hi2(N), hs2(N), Tsfc2(N), S2(NL, N), T2(NL, N)
      type(ice_ledger_t) :: l, l2, l3
      character(len=16) :: nm
      integer :: c, k

      call ice_random_columns(N, NL, 7_int64, hi, hs, Tsfc, S, T)
      call ice_random_columns(N, NL, 7_int64, hi2, hs2, Tsfc2, S2, T2)
      call check(error, all(hi == hi2) .and. all(T == T2) .and. all(S == S2), "deterministic")
      if (allocated(error)) return
      call check(error, all(hi >= 0.05_wp .and. hi <= 5.0_wp) .and. all(hs >= 0.0_wp .and. hs <= 0.5_wp), &
                 "thickness ranges")
      if (allocated(error)) return
      do c = 1, N
         do k = 1, NL
            call check(error, T(k, c) < -0.054_wp*S(k, c), "layer below its melting point")
            if (allocated(error)) return
         end do
      end do
      call l%add("mass", 1.0e-16_wp, 1.0_wp)
      call l%add("heat", 5.0e-14_wp, 1.0e3_wp)
      call l%add("mass", 2.0e-16_wp, 1.0_wp)
      call check(error, l%n == 2 .and. abs(l%worst_of("mass") - 2.0e-16_wp) < 1.0e-30_wp, "ledger lines")
      if (allocated(error)) return
      call check(error, l%closed(1.0e-15_wp) .and. .not. l%closed(1.0e-17_wp), "ledger closure")
      if (allocated(error)) return
      ! A NaN residual is never hidden, and overflowing the ledger is never silent.
      call l2%add("nan", ieee_value(1.0_wp, ieee_quiet_nan), 1.0_wp)
      call check(error,.not. l2%closed(1.0_wp) .and. l2%worst_of("nan") > 1.0e300_wp, "NaN residual fails")
      if (allocated(error)) return
      do c = 1, ICE_LEDGER_MAX + 1
         write (nm, '(a,i0)') "line", c
         call l3%add(trim(nm), 0.0_wp, 1.0_wp)
      end do
      call check(error,.not. l3%closed(1.0_wp) .and. l3%worst_of("line1") > 1.0e300_wp, "ledger overflow fails")
      if (allocated(error)) return
      call check(error, ice_golden_first_bad([1.0_wp, 2.0_wp], [1.0_wp, 2.0_wp + 1.0e-3_wp], 1.0e-6_wp, 0.0_wp) == 2, &
                 "golden check detects a mismatch")
      call check(error, ice_golden_first_bad([1.0_wp, 2.0_wp], [1.0_wp, 2.0_wp], 0.0_wp, 0.0_wp) == 0, &
                 "golden check accepts equality")
   end subroutine test_util

end module test_cice_conversions
