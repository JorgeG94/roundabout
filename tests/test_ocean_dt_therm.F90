!! Unit tests for the DT_THERM ratio helpers on `ocean_dyn_t`:
!! `is_thermo_step()` and `therm_dt()`.  Phase 7 of the dycore-
!! stabilizer plan (MOM6 DT_THERM analogue).
module test_ocean_dt_therm
!! COVERAGE GAP (2026-09-11) — READ BEFORE TRUSTING THIS FILE.
!!
!! Every test below checks BOOKKEEPING: that the thermo cadence fires on the
!! right steps and reports the right `dt`. None checks that the SOLUTION is
!! correct. All five passed continuously while `dt_therm_ratio > 1` was
!! destabilising the model: gating the EOS on `is_thermo_step()`
!! zero-order-held the baroclinic PGF for tau = (ratio-1)*dt, giving a
!! delayed-restoring-force instability of the grid-scale internal gravity
!! wave (sigma ~ omega^2*tau/2, unstable for ANY tau > 0). ~14 shipped
!! namelists were exposed, surviving only by viscosity margin. Fixed by
!! recomputing the EOS every dynamics step.
!!
!! The gate that WOULD have caught it: an unforced adiabatic run's kinetic
!! energy must not grow. On the eady front the control holds
!! En 2.676E-04 -> 2.454E-04 over 25 days, while every broken variant grew
!! 30-60x. That needs a driver-level integration, so it is NOT implemented
!! here. Reproducer and measurements:
!! `tmp_local_artifacts/eady_hunt/FINDINGS_windowed_advect.md`.
!!
!! If you add cadence knobs, add a PHYSICS assertion, not another parity
!! check.

   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_ocean_dyn, only: ocean_dyn_t
   implicit none
   private

   public :: collect_ocean_dt_therm_tests

contains

   subroutine collect_ocean_dt_therm_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("ratio_one_always_thermo_step", test_ratio_one), &
                  new_unittest("ratio_one_dt_unchanged", test_ratio_one_dt), &
                  new_unittest("ratio_two_alternates", test_ratio_two), &
                  new_unittest("ratio_two_scales_dt", test_ratio_two_dt), &
                  new_unittest("ratio_three_aligns_with_zero", test_ratio_three) &
                  ]
   end subroutine collect_ocean_dt_therm_tests

   subroutine test_ratio_one(error)
      !! `dt_therm_ratio = 1` (default) must report every step as a
      !! thermo step regardless of step count — bit-identical to
      !! pre-Phase-7 behaviour.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_dyn_t) :: dyn
      integer :: s
      logical :: all_true

      dyn%dt_therm_ratio = 1
      all_true = .true.
      do s = 0, 10
         dyn%outer_step_count = s
         if (.not. dyn%is_thermo_step()) all_true = .false.
      end do
      call check(error, all_true, &
                 "ratio = 1 must report every step as thermo")
   end subroutine test_ratio_one

   subroutine test_ratio_one_dt(error)
      !! With ratio = 1, `therm_dt(dt) = dt` exactly.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_dyn_t) :: dyn
      real(wp) :: dt_obs

      dyn%dt_therm_ratio = 1
      dt_obs = dyn%therm_dt(1200.0_wp)
      call check(error, abs(dt_obs - 1200.0_wp) < 1.0e-14_wp, &
                 "ratio = 1: therm_dt should equal dt")
   end subroutine test_ratio_one_dt

   subroutine test_ratio_two(error)
      !! `dt_therm_ratio = 2`: thermo fires on even step counts
      !! (0, 2, 4, ...) and skips odd ones.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_dyn_t) :: dyn
      logical :: ok
      integer :: s

      dyn%dt_therm_ratio = 2
      ok = .true.
      do s = 0, 9
         dyn%outer_step_count = s
         if (mod(s, 2) == 0 .and. .not. dyn%is_thermo_step()) ok = .false.
         if (mod(s, 2) /= 0 .and. dyn%is_thermo_step()) ok = .false.
      end do
      call check(error, ok, "ratio = 2 must alternate thermo on/off by step parity")
   end subroutine test_ratio_two

   subroutine test_ratio_two_dt(error)
      !! With ratio = 2, `therm_dt(dt) = 2·dt` so the gated kernels
      !! advance by the cumulative interval since their last fire.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_dyn_t) :: dyn
      real(wp) :: dt_obs

      dyn%dt_therm_ratio = 2
      dt_obs = dyn%therm_dt(1200.0_wp)
      call check(error, abs(dt_obs - 2400.0_wp) < 1.0e-14_wp, &
                 "ratio = 2: therm_dt should equal 2·dt")
   end subroutine test_ratio_two_dt

   subroutine test_ratio_three(error)
      !! `dt_therm_ratio = 3`: thermo fires on steps 0, 3, 6, 9, ...
      !! and `therm_dt = 3·dt`.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_dyn_t) :: dyn
      integer :: s, count_thermo
      real(wp) :: dt_obs

      dyn%dt_therm_ratio = 3
      count_thermo = 0
      do s = 0, 11
         dyn%outer_step_count = s
         if (dyn%is_thermo_step()) count_thermo = count_thermo + 1
      end do
      ! Steps 0..11 = 12 steps; thermo fires at 0, 3, 6, 9 = 4 times.
      call check(error, count_thermo == 4, &
                 "ratio = 3 over 12 steps should fire exactly 4 times")
      dt_obs = dyn%therm_dt(600.0_wp)
      if (.not. allocated(error)) call check(error, abs(dt_obs - 1800.0_wp) < 1.0e-14_wp, &
                                             "ratio = 3: therm_dt should equal 3·dt")
   end subroutine test_ratio_three

end module test_ocean_dt_therm
