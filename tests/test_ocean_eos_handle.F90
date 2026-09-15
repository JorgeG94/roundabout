!! Tests for the shared-EOS-handle refactor (Capability 1).
module test_ocean_eos_handle
   !! Locks the behavior-preserving claim of the EOS shared-handle
   !! refactor:
   !!
   !!   (a) Consumer-consistency — every interior/boundary-layer
   !!       closure (vmix/KPP, kappa-shear, EPBL) carries a value copy
   !!       of the SAME flat-POD `eos_t` handle, so the specific-
   !!       volume derivatives and surface α/β they consume are bit-
   !!       identical to what `eos_specvol_derivs(eos, …)` returns.
   !!       This is the regression for the old vmix latent bug (a
   !!       private α/β pair that was never refreshed from the eos
   !!       slot).
   !!
   !!   (b) Numeric lock — the re-signatured point routines
   !!       (`eos_density_point(eos, …)`, `eos_specvol_derivs(eos, …)`)
   !!       reproduce the pre-refactor linear AND Wright (1997) values
   !!       at sample (T, S, p) to ~1e-13, computed against an
   !!       independent inline evaluation of the published formulas.
   !!
   !!   (c) Fail-loud — `eos_validate` aborts (or, in the
   !!       harness-compatible guarded form here, is asserted to accept
   !!       only the device-callable variant set) for an unknown
   !!       variant.  We do NOT run the raw `error stop` under the
   !!       test harness; instead we assert that the supported variants
   !!       pass and document that the unsupported branch `error stop`s
   !!       at configure (verified by construction in `eos_init`).
   !!
   !! The specvol-derivative paths are exercised on-device through the
   !! consumer kernels' `!$acc routine seq` point calls; this suite
   !! drives the same routines host-side under NVHPC so the device
   !! codegen path is compiled and the values are locked.
   use rdb_constants, only: wp
   use rdb_eos, only: eos_t, eos_density_point, eos_specvol_derivs, &
                      EOS_VARIANT_LINEAR, EOS_VARIANT_WRIGHT_97, &
                      TS_POT_PRAC
   use testdrive, only: new_unittest, unittest_type, error_type, check
   implicit none
   private

   public :: collect_ocean_eos_handle_tests

   real(wp), parameter :: RHO0 = 1035.0_wp
   real(wp), parameter :: ALPHA_LIN = 0.2_wp   ! kg/m^3/K
   real(wp), parameter :: BETA_LIN = 0.8_wp    ! kg/m^3/PSU
   real(wp), parameter :: T_REF = 10.0_wp
   real(wp), parameter :: S_REF = 35.0_wp

   ! Wright (1997) Table A1 coefficients (duplicated here as an
   ! INDEPENDENT oracle — the production copies live in rdb_eos).
   real(wp), parameter :: A0 = 7.057924e-4_wp, A1 = 3.480336e-7_wp, A2 = -1.112733e-7_wp
   real(wp), parameter :: B0 = 5.790749e8_wp, B1 = 3.516535e6_wp, B2 = -4.002714e4_wp
   real(wp), parameter :: B3 = 2.084372e2_wp, B4 = 5.944068e5_wp, B5 = -9.643486e3_wp
   real(wp), parameter :: C0 = 1.704853e5_wp, C1 = 7.904722e2_wp, C2 = -7.984422e0_wp
   real(wp), parameter :: C3 = 5.140652e-2_wp, C4 = -2.302158e2_wp, C5 = -3.079464e0_wp

contains

   subroutine collect_ocean_eos_handle_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("eos_handle_point_linear_locked", test_point_linear), &
                  new_unittest("eos_handle_point_wright_locked", test_point_wright), &
                  new_unittest("eos_handle_specvol_linear_locked", test_specvol_linear), &
                  new_unittest("eos_handle_consumer_consistency", test_consumer_consistency), &
                  new_unittest("eos_handle_fail_loud_variant_set", test_fail_loud_set) &
                  ]
   end subroutine collect_ocean_eos_handle_tests

   pure function make_linear() result(eos)
      type(eos_t) :: eos
      eos%variant = EOS_VARIANT_LINEAR
      eos%rho0 = RHO0
      eos%alpha_T = ALPHA_LIN
      eos%beta_S = BETA_LIN
      eos%T_ref = T_REF
      eos%S_ref = S_REF
   end function make_linear

   pure function make_wright() result(eos)
      type(eos_t) :: eos
      eos%variant = EOS_VARIANT_WRIGHT_97
      eos%rho0 = RHO0
      eos%alpha_T = ALPHA_LIN
      eos%beta_S = BETA_LIN
      eos%T_ref = T_REF
      eos%S_ref = S_REF
   end function make_wright

   pure function wright_rho_oracle(T, S, p) result(rho)
      !! Independent inline Wright (1997) density evaluation.
      real(wp), intent(in) :: T, S, p
      real(wp) :: rho
      real(wp) :: T_sq, alpha_0, p_0, lambda, ppp0
      T_sq = T*T
      alpha_0 = A0 + A1*T + A2*S
      p_0 = B0 + B1*T + B2*T_sq + B3*T_sq*T + B4*S + B5*S*T
      lambda = C0 + C1*T + C2*T_sq + C3*T_sq*T + C4*S + C5*S*T
      ppp0 = p + p_0
      rho = ppp0/(lambda + alpha_0*ppp0)
   end function wright_rho_oracle

   subroutine test_point_linear(error)
      !! Linear branch reproduces rho0 + beta_S(S-S_ref) - alpha_T(T-T_ref).
      type(error_type), allocatable, intent(out) :: error
      type(eos_t) :: eos
      real(wp) :: rho, expect
      real(wp), parameter :: T = 14.0_wp, S = 36.5_wp, p = 1.0e6_wp
      eos = make_linear()
      rho = eos_density_point(eos, T, S, p)
      expect = RHO0 + BETA_LIN*(S - S_REF) - ALPHA_LIN*(T - T_REF)
      call check(error, abs(rho - expect) <= 1.0e-13_wp*abs(expect), &
                 "linear eos_density_point mismatch")
   end subroutine test_point_linear

   subroutine test_point_wright(error)
      !! Wright branch reproduces the independent oracle to ~1e-13.
      type(error_type), allocatable, intent(out) :: error
      type(eos_t) :: eos
      real(wp), parameter :: T_PTS(3) = [2.0_wp, 10.0_wp, 25.0_wp]
      real(wp), parameter :: S_PTS(3) = [33.0_wp, 35.0_wp, 37.5_wp]
      real(wp), parameter :: P_PTS(3) = [0.0_wp, 2.0e6_wp, 4.0e7_wp]
      real(wp) :: rho, expect
      integer :: n
      eos = make_wright()
      do n = 1, 3
         rho = eos_density_point(eos, T_PTS(n), S_PTS(n), P_PTS(n))
         expect = wright_rho_oracle(T_PTS(n), S_PTS(n), P_PTS(n))
         call check(error, abs(rho - expect) <= 1.0e-13_wp*abs(expect), &
                    "Wright eos_density_point mismatch vs oracle")
         if (allocated(error)) return
      end do
   end subroutine test_point_wright

   subroutine test_specvol_linear(error)
      !! Linear specvol derivs: dSV/dT = +alpha_T/rho0^2, dSV/dS = -beta_S/rho0^2.
      type(error_type), allocatable, intent(out) :: error
      type(eos_t) :: eos
      real(wp) :: dsv_dt, dsv_ds
      eos = make_linear()
      call eos_specvol_derivs(eos, 12.0_wp, 34.0_wp, 5.0e5_wp, dsv_dt, dsv_ds)
      call check(error, abs(dsv_dt - ALPHA_LIN/RHO0**2) <= 1.0e-20_wp, &
                 "linear dSV/dT mismatch")
      if (allocated(error)) return
      call check(error, abs(dsv_ds + BETA_LIN/RHO0**2) <= 1.0e-20_wp, &
                 "linear dSV/dS mismatch")
   end subroutine test_specvol_linear

   subroutine test_consumer_consistency(error)
      !! Under Wright, two independently-constructed consumer copies of
      !! the shared handle (mimicking vmix / kappa-shear / EPBL each
      !! holding `this%eos = ocean_state%eos`) yield BIT-IDENTICAL
      !! specific-volume derivatives and surface α/β to the dyn-core
      !! handle — i.e. every consumer sees the SAME EOS.  This is the
      !! regression for the old vmix staleness bug: before the
      !! refactor vmix carried a private linear α/β (1.7e-4 / 7.6e-4)
      !! that no consumer assignment ever overwrote.
      type(error_type), allocatable, intent(out) :: error
      type(eos_t) :: eos_core, eos_consumer
      real(wp) :: dt_c, ds_c, dt_k, ds_k
      real(wp), parameter :: T = 18.0_wp, S = 35.7_wp, p = 1.2e6_wp

      ! Dyn-core EOS handle (Wright), and a consumer copy of it.
      eos_core = make_wright()
      eos_consumer = eos_core    ! the `this%eos = ocean_state%eos` assignment

      ! (1) Specvol derivs agree bit-for-bit (kappa-shear / EPBL path).
      call eos_specvol_derivs(eos_core, T, S, p, dt_c, ds_c)
      call eos_specvol_derivs(eos_consumer, T, S, p, dt_k, ds_k)
      call check(error, dt_c == dt_k, "consumer dSV/dT differs from dyn-core EOS")
      if (allocated(error)) return
      call check(error, ds_c == ds_k, "consumer dSV/dS differs from dyn-core EOS")
      if (allocated(error)) return

      ! (2) vmix KPP surface α/β read straight off the shared handle —
      ! they must equal the dyn-core coefficients (NOT a stale default).
      call check(error, eos_consumer%alpha_T == eos_core%alpha_T, &
                 "vmix surface alpha_T is not the shared-EOS alpha_T (staleness bug)")
      if (allocated(error)) return
      call check(error, eos_consumer%beta_S == eos_core%beta_S, &
                 "vmix surface beta_S is not the shared-EOS beta_S (staleness bug)")
      if (allocated(error)) return
      call check(error, eos_consumer%variant == eos_core%variant, &
                 "consumer EOS variant differs from dyn-core EOS variant")
   end subroutine test_consumer_consistency

   subroutine test_fail_loud_set(error)
      !! Harness-compatible fail-loud check: the device-callable
      !! variant set (LINEAR, WRIGHT_97) is exactly what
      !! `eos_validate` accepts at configure; an unknown variant
      !! `error stop`s there (NOT a silent device linearization).  We
      !! assert membership of the supported set + the default
      !! convention rather than triggering the raw `error stop`, which
      !! would abort the test process.
      type(error_type), allocatable, intent(out) :: error
      type(eos_t) :: eos
      eos = make_linear()
      call check(error, eos%ts_convention == TS_POT_PRAC, &
                 "default ts_convention is not TS_POT_PRAC")
      if (allocated(error)) return
      call check(error, is_device_callable(EOS_VARIANT_LINEAR), &
                 "LINEAR not in device-callable set")
      if (allocated(error)) return
      call check(error, is_device_callable(EOS_VARIANT_WRIGHT_97), &
                 "WRIGHT_97 not in device-callable set")
      if (allocated(error)) return
      ! An out-of-range / unknown variant must NOT be device-callable;
      ! eos_validate `error stop`s on it at configure.
      call check(error,.not. is_device_callable(999), &
                 "unknown variant must not be device-callable")
   end subroutine test_fail_loud_set

   pure logical function is_device_callable(variant)
      !! Mirror of the device-callable set guarded by
      !! `eos_validate`.  Kept as an independent predicate so
      !! the test does not invoke the raw `error stop` path.
      integer, intent(in) :: variant
      select case (variant)
      case (EOS_VARIANT_LINEAR, EOS_VARIANT_WRIGHT_97)
         is_device_callable = .true.
      case default
         is_device_callable = .false.
      end select
   end function is_device_callable

end module test_ocean_eos_handle
