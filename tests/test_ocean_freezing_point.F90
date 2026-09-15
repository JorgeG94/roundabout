!! Unit tests for `eos_freezing_point` (`rdb_eos`) — the
!! ocean-side sea-ice prerequisite (PLAN_SEA_ICE.md, PR 1).
!!
!! v1 is the SIS2/MOM6 linear liquidus for every EOS variant:
!!   T_f = TFR_S_COEFF·S + TFR_P_COEFF·p
!! with TFR_S_COEFF = −0.054 °C/(g/kg) and TFR_P_COEFF = −7.53e-8 °C/Pa.
!!
!! Cases:
!!   * Reference point — T_f(S=35, p=0) = −1.89 °C exactly.
!!   * Fresh water — T_f(0, 0) = 0.
!!   * Monotonic decreasing in S over [0, 40].
!!   * Pressure depression — T_f(S, p>0) < T_f(S, 0) by TFR_P_COEFF·p.
!!   * Elemental broadcast — one call over a salinity ARRAY matches the
!!     scalar loop element-by-element.
!!   * Variant-invariant (v1) — linear / Wright / Roquet handles all
!!     return the same liquidus.
module test_ocean_freezing_point
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_eos, only: eos_t, eos_freezing_point, &
                      TFR_S_COEFF, TFR_P_COEFF, &
                      EOS_VARIANT_LINEAR, EOS_VARIANT_WRIGHT_97, &
                      EOS_VARIANT_ROQUET_SPV
   implicit none
   private

   public :: collect_ocean_freezing_point_tests

contains

   subroutine collect_ocean_freezing_point_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("freezing_point_reference_s35", test_reference_s35), &
                  new_unittest("freezing_point_fresh_water_zero", test_fresh_zero), &
                  new_unittest("freezing_point_monotonic_in_s", test_monotonic_s), &
                  new_unittest("freezing_point_pressure_depression", test_pressure), &
                  new_unittest("freezing_point_elemental_array", test_elemental_array), &
                  new_unittest("freezing_point_variant_invariant_v1", test_variants) &
                  ]
   end subroutine collect_ocean_freezing_point_tests

   subroutine make_eos(eos)
      !! Minimal EOS handle — the point function only reads the flat POD,
      !! so a default-init handle (linear variant) suffices.
      type(eos_t), intent(out) :: eos
      type(hgrid_t) :: grid
      call grid%init(4, 4, 2, 1.0_wp, 1.0_wp)
      call eos%init(grid)
   end subroutine make_eos

   subroutine test_reference_s35(error)
      !! SIS2 reference point: T_f(35 PSU, surface) = −0.054·35 = −1.89 °C.
      type(error_type), allocatable, intent(out) :: error
      type(eos_t) :: eos
      real(wp) :: t_f

      call make_eos(eos)
      t_f = eos_freezing_point(eos, 35.0_wp, 0.0_wp)
      call check(error, abs(t_f - (-1.89_wp)) < 1.0e-12_wp, &
                 "T_f(S=35, p=0) must be -1.89 degC")
      call eos%destroy()
   end subroutine test_reference_s35

   subroutine test_fresh_zero(error)
      !! Fresh water at the surface freezes at exactly 0 °C.
      type(error_type), allocatable, intent(out) :: error
      type(eos_t) :: eos
      real(wp) :: t_f

      call make_eos(eos)
      t_f = eos_freezing_point(eos, 0.0_wp, 0.0_wp)
      call check(error, abs(t_f) < 1.0e-15_wp, "T_f(0, 0) must be 0")
      call eos%destroy()
   end subroutine test_fresh_zero

   subroutine test_monotonic_s(error)
      !! Strictly decreasing in S: saltier water freezes colder.
      type(error_type), allocatable, intent(out) :: error
      type(eos_t) :: eos
      real(wp) :: t_prev, t_cur, s
      integer :: is
      logical :: monotone

      call make_eos(eos)
      monotone = .true.
      t_prev = eos_freezing_point(eos, 0.0_wp, 0.0_wp)
      do is = 1, 40
         s = real(is, wp)
         t_cur = eos_freezing_point(eos, s, 0.0_wp)
         if (t_cur >= t_prev) monotone = .false.
         t_prev = t_cur
      end do
      call check(error, monotone, "T_f must be strictly decreasing in S")
      call eos%destroy()
   end subroutine test_monotonic_s

   subroutine test_pressure(error)
      !! Pressure depression: at p = 1e7 Pa (~1000 m) T_f drops by
      !! TFR_P_COEFF·p = −0.753 °C relative to the surface value.
      type(error_type), allocatable, intent(out) :: error
      type(eos_t) :: eos
      real(wp), parameter :: P_DEEP = 1.0e7_wp
      real(wp) :: t_sfc, t_deep

      call make_eos(eos)
      t_sfc = eos_freezing_point(eos, 35.0_wp, 0.0_wp)
      t_deep = eos_freezing_point(eos, 35.0_wp, P_DEEP)
      call check(error, abs((t_deep - t_sfc) - TFR_P_COEFF*P_DEEP) < 1.0e-12_wp, &
                 "pressure depression must be TFR_P_COEFF*p")
      if (allocated(error)) return
      call check(error, t_deep < t_sfc, "deep freezing point must be colder")
      call eos%destroy()
   end subroutine test_pressure

   subroutine test_elemental_array(error)
      !! `elemental` contract: one call over an S ARRAY (scalar eos + p
      !! broadcast) returns the per-element scalar results.
      type(error_type), allocatable, intent(out) :: error
      type(eos_t) :: eos
      real(wp) :: s_arr(5), t_arr(5)
      real(wp) :: max_err
      integer :: k

      call make_eos(eos)
      s_arr = [0.0_wp, 5.0_wp, 20.0_wp, 35.0_wp, 40.0_wp]
      t_arr = eos_freezing_point(eos, s_arr, 0.0_wp)
      max_err = 0.0_wp
      do k = 1, size(s_arr)
         max_err = max(max_err, &
                       abs(t_arr(k) - eos_freezing_point(eos, s_arr(k), 0.0_wp)))
      end do
      call check(error, max_err < 1.0e-15_wp, &
                 "elemental array call must match scalar calls exactly")
      if (allocated(error)) return
      call check(error, abs(t_arr(4) - TFR_S_COEFF*35.0_wp) < 1.0e-12_wp, &
                 "array element S=35 must give -1.89 degC")
      call eos%destroy()
   end subroutine test_elemental_array

   subroutine test_variants(error)
      !! v1 contract: the linear liquidus applies to ALL EOS variants
      !! (matching MOM6's TFREEZE_FORM="LINEAR" default under any density
      !! branch).  When a TEOS-10 branch lands, this test pins the
      !! variants it must NOT change.
      type(error_type), allocatable, intent(out) :: error
      type(eos_t) :: eos
      real(wp) :: t_lin, t_wright, t_roquet

      call make_eos(eos)
      eos%variant = EOS_VARIANT_LINEAR
      t_lin = eos_freezing_point(eos, 35.0_wp, 0.0_wp)
      eos%variant = EOS_VARIANT_WRIGHT_97
      t_wright = eos_freezing_point(eos, 35.0_wp, 0.0_wp)
      eos%variant = EOS_VARIANT_ROQUET_SPV
      t_roquet = eos_freezing_point(eos, 35.0_wp, 0.0_wp)
      call check(error, abs(t_lin - t_wright) < 1.0e-15_wp .and. &
                 abs(t_lin - t_roquet) < 1.0e-15_wp, &
                 "v1: every EOS variant must return the linear liquidus")
      call eos%destroy()
   end subroutine test_variants

end module test_ocean_freezing_point
