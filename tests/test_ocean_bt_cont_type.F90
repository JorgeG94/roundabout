!! Unit tests for the MOM6 BT_cont_type flux-bounded continuity
!! helpers in `rdb_bt_cont_type`.  Covers tests #1-4 and #8:
!!
!!   1. `find_uhbt(0, BTC) = 0` for any BTC.
!!   2. `find_uhbt` is monotone on [uBT_EE, uBT_WW] for a
!!      realistic synthetic BTC.
!!   3. The derivative `find_duhbt_du` is continuous at
!!      `uBT_EE` and `uBT_WW` (to roundoff).
!!   4. Linear limit: with `uh_crv=0` and saturation thresholds
!!      pushed to the bounds, `find_uhbt(u, BTC) = u·FA_u_W0`
!!      for u > 0 (and the symmetric statement for u < 0).
!!   8. `uhbt_to_ubt(find_uhbt(u, BTC), BTC) ≈ u` to Newton
!!      tolerance over a sweep of u inside the cubic interval.
!!
!! Each test uses a synthetic `local_BT_cont_u_type` built in
!! `make_realistic_btc` — a 1000 m × 4000 m face with bed-supplied
!! saturation thresholds asymmetric in the east/west direction.
!! The C¹-matching coefficients (`uh_crvE/W`, `uh_EE/WW`) are
!! computed analytically so the BTC is internally consistent.
module test_ocean_bt_cont_type
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_barotropic_workstate, only: barotropic_workstate_t, &
                                       local_BT_cont_u_type
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_bt_cont_type, only: find_uhbt, find_duhbt_du, uhbt_to_ubt
   use rdb_barotropic_coupling, only: set_local_BT_cont_types
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   implicit none
   private

   public :: collect_ocean_bt_cont_type_tests

contains

   subroutine collect_ocean_bt_cont_type_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("bt_cont_zero_velocity", test_zero_velocity), &
                  new_unittest("bt_cont_monotone", test_monotone), &
                  new_unittest("bt_cont_c1_continuity", test_c1_continuity), &
                  new_unittest("bt_cont_linear_limit", test_linear_limit), &
                  new_unittest("bt_cont_roundtrip", test_roundtrip), &
                  new_unittest("bt_cont_saturation_caps_transport", test_saturation_cap), &
                  new_unittest("bt_cont_set_flat_bath", test_set_flat_bath), &
                  new_unittest("bt_cont_set_vanished_bed", test_set_vanished_bed), &
                  new_unittest("bt_cont_set_off_is_noop", test_set_off_is_noop) &
                  ]
   end subroutine collect_ocean_bt_cont_type_tests

   subroutine make_realistic_btc(BTC)
      !! Build a BTC for a "shelf-break" face: deep west side,
      !! shallow east side.  The cubic branch curvatures and
      !! intercepts are computed so the function is C¹ at the
      !! threshold velocities — i.e. `uhbt` and `duhbt/du` match
      !! across the boundary.
      type(local_BT_cont_u_type), intent(out) :: BTC

      real(wp) :: u_th_W, u_th_E

      BTC%FA_u_W0 = 1000.0_wp   ! near-zero face area, u > 0 (m)
      BTC%FA_u_WW = 600.0_wp    ! saturated face area, u > uBT_WW
      BTC%FA_u_E0 = 800.0_wp    ! near-zero face area, u < 0
      BTC%FA_u_EE = 300.0_wp    ! saturated face area, u < uBT_EE

      BTC%uBT_WW = 1.5_wp       ! threshold (m/s, positive)
      BTC%uBT_EE = -2.0_wp      ! threshold (m/s, negative)

      u_th_W = BTC%uBT_WW
      u_th_E = BTC%uBT_EE
      ! C¹ matching at the positive threshold:
      !   value: u_th·(FA_W0 + crvW·u_th²) = (u_th - u_th)·FA_WW + uh_WW
      !          ⇒  uh_WW = u_th·(FA_W0 + crvW·u_th²)
      !   slope: FA_W0 + 3·crvW·u_th² = FA_WW
      !          ⇒  crvW = (FA_WW - FA_W0)/(3·u_th²)
      BTC%uh_crvW = (BTC%FA_u_WW - BTC%FA_u_W0)/(3.0_wp*u_th_W*u_th_W)
      BTC%uh_WW = u_th_W*(BTC%FA_u_W0 + BTC%uh_crvW*u_th_W*u_th_W)
      BTC%uh_crvE = (BTC%FA_u_EE - BTC%FA_u_E0)/(3.0_wp*u_th_E*u_th_E)
      BTC%uh_EE = u_th_E*(BTC%FA_u_E0 + BTC%uh_crvE*u_th_E*u_th_E)
   end subroutine make_realistic_btc

   subroutine test_zero_velocity(error)
      type(error_type), allocatable, intent(out) :: error
      type(local_BT_cont_u_type) :: BTC
      real(wp) :: uhbt

      call make_realistic_btc(BTC)
      uhbt = find_uhbt(0.0_wp, BTC)
      call check(error, abs(uhbt) < 1.0e-30_wp, "find_uhbt(0) is not zero")
      if (allocated(error)) return

      ! Also assert for an all-zero (default-initialised) BTC.
      block
         type(local_BT_cont_u_type) :: BTC_zero
         uhbt = find_uhbt(0.0_wp, BTC_zero)
         call check(error, abs(uhbt) < 1.0e-30_wp, &
                    "find_uhbt(0) on default BTC is not zero")
      end block
   end subroutine test_zero_velocity

   subroutine test_monotone(error)
      !! Sweep u over [-2.5, +2.0] and verify `find_uhbt` is
      !! monotonically non-decreasing.  Hits all four branches:
      !!   u < uBT_EE         (saturated negative)
      !!   uBT_EE ≤ u < 0     (cubic negative)
      !!   0 ≤ u ≤ uBT_WW     (cubic positive)
      !!   u > uBT_WW         (saturated positive)
      type(error_type), allocatable, intent(out) :: error
      type(local_BT_cont_u_type) :: BTC
      integer, parameter :: N = 201
      real(wp) :: u, du, prev_uhbt, uhbt
      integer :: i

      call make_realistic_btc(BTC)
      du = 4.5_wp/real(N - 1, wp)
      prev_uhbt = find_uhbt(-2.5_wp, BTC)
      do i = 2, N
         u = -2.5_wp + real(i - 1, wp)*du
         uhbt = find_uhbt(u, BTC)
         call check(error, uhbt >= prev_uhbt - 1.0e-12_wp, &
                    "find_uhbt not monotone")
         if (allocated(error)) return
         prev_uhbt = uhbt
      end do
   end subroutine test_monotone

   subroutine test_c1_continuity(error)
      !! Check that `find_duhbt_du` matches across each threshold
      !! to within roundoff.  By construction `crvW/E` are picked
      !! so the cubic-side slope equals the saturated slope at
      !! `u = ±uBT_WW/EE`.
      type(error_type), allocatable, intent(out) :: error
      type(local_BT_cont_u_type) :: BTC
      real(wp), parameter :: eps_th = 1.0e-8_wp
      real(wp) :: d_cubic, d_satur

      call make_realistic_btc(BTC)

      ! Positive threshold: cubic side at u = uBT_WW − ε, saturated at uBT_WW + ε.
      d_cubic = find_duhbt_du(BTC%uBT_WW - eps_th, BTC)
      d_satur = find_duhbt_du(BTC%uBT_WW + eps_th, BTC)
      call check(error, abs(d_cubic - d_satur) < 1.0e-6_wp*max(1.0_wp, abs(BTC%FA_u_WW)), &
                 "C1 mismatch at uBT_WW")
      if (allocated(error)) return

      ! Negative threshold.
      d_cubic = find_duhbt_du(BTC%uBT_EE + eps_th, BTC)
      d_satur = find_duhbt_du(BTC%uBT_EE - eps_th, BTC)
      call check(error, abs(d_cubic - d_satur) < 1.0e-6_wp*max(1.0_wp, abs(BTC%FA_u_EE)), &
                 "C1 mismatch at uBT_EE")
   end subroutine test_c1_continuity

   subroutine test_linear_limit(error)
      !! With `uh_crvW = 0` AND `uBT_WW` pushed beyond the test
      !! velocity, the positive cubic branch collapses to
      !!   uhbt = u·FA_u_W0
      !! — the naive `uh = u·h_face` form we want the off-switch
      !! to reproduce bit-identically.
      type(error_type), allocatable, intent(out) :: error
      type(local_BT_cont_u_type) :: BTC
      real(wp) :: u, uhbt, expected

      BTC%FA_u_W0 = 1234.5_wp
      BTC%FA_u_E0 = 1234.5_wp
      BTC%FA_u_WW = 1234.5_wp
      BTC%FA_u_EE = 1234.5_wp
      BTC%uBT_WW = 1.0e6_wp
      BTC%uBT_EE = -1.0e6_wp
      BTC%uh_crvW = 0.0_wp
      BTC%uh_crvE = 0.0_wp
      BTC%uh_WW = BTC%uBT_WW*BTC%FA_u_W0
      BTC%uh_EE = BTC%uBT_EE*BTC%FA_u_E0

      u = 0.7_wp
      uhbt = find_uhbt(u, BTC)
      expected = u*BTC%FA_u_W0
      call check(error, abs(uhbt - expected) < 1.0e-12_wp*abs(expected), &
                 "linear-limit mismatch on positive side")
      if (allocated(error)) return

      u = -0.5_wp
      uhbt = find_uhbt(u, BTC)
      expected = u*BTC%FA_u_E0
      call check(error, abs(uhbt - expected) < 1.0e-12_wp*abs(expected), &
                 "linear-limit mismatch on negative side")
   end subroutine test_linear_limit

   subroutine test_roundtrip(error)
      !! For a sweep of test velocities in [uBT_EE, uBT_WW] (the
      !! interesting cubic range), confirm that
      !!   uhbt_to_ubt(find_uhbt(u, BTC), BTC) ≈ u
      !! within Newton-iteration tolerance.  Also tests the
      !! saturated branches at u = 2·uBT_WW and u = 2·uBT_EE.
      type(error_type), allocatable, intent(out) :: error
      type(local_BT_cont_u_type) :: BTC
      real(wp), parameter :: u_samples(*) = [ &
                             real(wp) :: &
                             -3.5_wp, -2.5_wp, -1.5_wp, -0.7_wp, -0.05_wp, &
                             0.05_wp, 0.4_wp, 1.0_wp, 1.5_wp, 2.5_wp]
      real(wp) :: u, uhbt, ubt_back
      integer :: i

      call make_realistic_btc(BTC)
      do i = 1, size(u_samples)
         u = u_samples(i)
         uhbt = find_uhbt(u, BTC)
         ubt_back = uhbt_to_ubt(uhbt, BTC)
         call check(error, abs(ubt_back - u) < 1.0e-8_wp*max(1.0_wp, abs(u)), &
                    "roundtrip residual exceeds Newton tol")
         if (allocated(error)) return
      end do
   end subroutine test_roundtrip

   subroutine make_btc_setup(grid, metrics, ms, bt_work, nz, h_west_per_k, h_east_per_k)
      !! Spin up a 3×1×nz multilayer state with two cells (west at
      !! `i=1`, east at `i=2`) and a single u-face at `i=2`.  Caller
      !! supplies the per-layer h on each side; the helper builds the
      !! grid, allocates the multilayer state, fills h_layer, allocates
      !! BTCL_u/v with `use_bt_cont_type=.true.`, and is ready for a
      !! call to `set_local_BT_cont_types`.
      type(hgrid_t), intent(out) :: grid
      type(ocean_metrics_t), intent(out) :: metrics
      type(multilayer_state_t), intent(out) :: ms
      type(barotropic_workstate_t), intent(out) :: bt_work
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_west_per_k(nz)
      real(wp), intent(in) :: h_east_per_k(nz)

      integer :: nx_phys, ny_phys, nghost, k

      nx_phys = 2
      ny_phys = 1
      nghost = 0
      call grid%init(nx_phys, ny_phys, nghost, 1000.0_wp, 1000.0_wp)

      ms%nz_ml = nz
      call ms%init(grid)

      do k = 1, nz
         ms%h_layer(1, 1, k) = h_west_per_k(k)
         ms%h_layer(2, 1, k) = h_east_per_k(k)
      end do

      call bt_work%init(grid, nz_ml=nz)
      bt_work%use_bt_cont_type = .true.
      allocate (bt_work%BTCL_u(grid%nx_total + 1, grid%ny_total))
      allocate (bt_work%BTCL_v(grid%nx_total, grid%ny_total + 1))
      call make_cartesian_metrics(metrics, grid)
   end subroutine make_btc_setup

   subroutine test_set_flat_bath(error)
      !! Two identical columns (each NZ=3 layers × 1000 m) with no
      !! vanishing bed.  Expect at the shared u-face:
      !!   FA_u_W0 = FA_u_E0 = FA_u_WW = FA_u_EE = 3000 m
      !!   uh_crvW = uh_crvE = 0
      !!   uh_WW = uBT_WW · 3000, uh_EE = uBT_EE · 3000.
      !! Hits plan test #5.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(barotropic_workstate_t) :: bt_work
      type(ocean_metrics_t) :: metrics
      real(wp) :: H_total, dt
      type(local_BT_cont_u_type) :: BTC

      H_total = 3000.0_wp
      dt = 1200.0_wp
      call make_btc_setup(grid, metrics, ms, bt_work, 3, &
                          [1000.0_wp, 1000.0_wp, 1000.0_wp], &
                          [1000.0_wp, 1000.0_wp, 1000.0_wp])

      call set_local_BT_cont_types(grid, metrics, bt_work, ms, dt)
      BTC = bt_work%BTCL_u(2, 1)

      call check(error, abs(BTC%FA_u_W0 - H_total) < 1.0e-10_wp, "FA_u_W0 ≠ H_total")
      if (allocated(error)) return
      call check(error, abs(BTC%FA_u_E0 - H_total) < 1.0e-10_wp, "FA_u_E0 ≠ H_total")
      if (allocated(error)) return
      call check(error, abs(BTC%FA_u_WW - H_total) < 1.0e-10_wp, "FA_u_WW ≠ H_total")
      if (allocated(error)) return
      call check(error, abs(BTC%FA_u_EE - H_total) < 1.0e-10_wp, "FA_u_EE ≠ H_total")
      if (allocated(error)) return
      call check(error, abs(BTC%uh_crvW) < 1.0e-12_wp, "uh_crvW ≠ 0 on flat bath")
      if (allocated(error)) return
      call check(error, abs(BTC%uh_crvE) < 1.0e-12_wp, "uh_crvE ≠ 0 on flat bath")
      if (allocated(error)) return
      ! Closed-form for uh_WW when FA_W0=FA_WW=H: uh_WW = uBT_WW · H.
      call check(error, abs(BTC%uh_WW - BTC%uBT_WW*H_total) < 1.0e-8_wp, &
                 "uh_WW ≠ uBT_WW·H on flat bath")
      if (allocated(error)) return
      call check(error, abs(BTC%uh_EE - BTC%uBT_EE*H_total) < 1.0e-8_wp, &
                 "uh_EE ≠ uBT_EE·H on flat bath")

      call ms%destroy()
      call bt_work%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_set_flat_bath

   subroutine test_set_vanished_bed(error)
      !! Spoon-margin east-west margin.  West column has the bed layer
      !! vanished (h_bed=0, h_surf=1000); east column has both layers
      !! present (h_bed=30, h_surf=1000).  Expected at the u-face:
      !!   FA_u_W0 = FA_u_E0 = h_face_bed + h_face_surf
      !!                     = 0.5·(0+30) + 0.5·(1000+1000)
      !!                     = 15 + 1000 = 1015
      !!   FA_u_WW = sum_k h_west(k) = 0 + 1000 = 1000   ← upstream west
      !!   FA_u_EE = sum_k h_east(k) = 30 + 1000 = 1030  ← upstream east
      !! Crucially `FA_u_WW < FA_u_W0` — the bed layer can't supply
      !! strong eastward flow because the west cell's bed layer is gone.
      !! Hits plan test #6.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(barotropic_workstate_t) :: bt_work
      type(ocean_metrics_t) :: metrics
      type(local_BT_cont_u_type) :: BTC

      call make_btc_setup(grid, metrics, ms, bt_work, 2, &
                          [0.0_wp, 1000.0_wp], &       ! west: vanished bed
                          [30.0_wp, 1000.0_wp])         ! east: full

      call set_local_BT_cont_types(grid, metrics, bt_work, ms, 1200.0_wp)
      BTC = bt_work%BTCL_u(2, 1)

      call check(error, abs(BTC%FA_u_W0 - 1015.0_wp) < 1.0e-10_wp, &
                 "FA_u_W0 ≠ centered-face sum")
      if (allocated(error)) return
      call check(error, abs(BTC%FA_u_WW - 1000.0_wp) < 1.0e-10_wp, &
                 "FA_u_WW ≠ west-column h sum")
      if (allocated(error)) return
      call check(error, abs(BTC%FA_u_EE - 1030.0_wp) < 1.0e-10_wp, &
                 "FA_u_EE ≠ east-column h sum")
      if (allocated(error)) return
      ! THE point of the closure: saturated face area on the west-draw
      ! side is BELOW the centered-face value.
      call check(error, BTC%FA_u_WW < BTC%FA_u_W0, &
                 "FA_u_WW must drop below FA_u_W0 at a vanished-bed margin")
      if (allocated(error)) return
      ! Curvature sign follows: westward-draw cubic must bend down.
      call check(error, BTC%uh_crvW < 0.0_wp, "uh_crvW should be negative")
      if (allocated(error)) return
      ! East-draw side has more head room → positive curvature.
      call check(error, BTC%uh_crvE > 0.0_wp, "uh_crvE should be positive")

      call ms%destroy()
      call bt_work%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_set_vanished_bed

   subroutine test_set_off_is_noop(error)
      !! Verify the producer is a true no-op when `use_bt_cont_type =
      !! .false.` — BTCL_u/v aren't allocated and the routine returns
      !! without touching anything.  Calling on an unallocated workstate
      !! field must not crash.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(barotropic_workstate_t) :: bt_work
      type(ocean_metrics_t) :: metrics
      integer :: nx_phys, ny_phys, nghost

      nx_phys = 2
      ny_phys = 1
      nghost = 0
      call grid%init(nx_phys, ny_phys, nghost, 1000.0_wp, 1000.0_wp)
      ms%nz_ml = 2
      call ms%init(grid)
      ms%h_layer(:, :, :) = 100.0_wp
      call bt_work%init(grid, nz_ml=2)
      ! Knob OFF — BTCL_u/v not allocated.
      bt_work%use_bt_cont_type = .false.

      call set_local_BT_cont_types(grid, metrics, bt_work, ms, 1200.0_wp)
      call check(error,.not. allocated(bt_work%BTCL_u), &
                 "BTCL_u allocated when knob is off")
      if (allocated(error)) return
      call check(error,.not. allocated(bt_work%BTCL_v), &
                 "BTCL_v allocated when knob is off")

      call ms%destroy()
      call bt_work%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_set_off_is_noop

   subroutine test_saturation_cap(error)
      !! Once `|u| > uBT_*` the marginal transport per unit velocity
      !! collapses to `FA_u_*` — i.e. doubling u outside the cap
      !! does NOT double the transport.  This is the actual
      !! property the BT_cont closure is designed to deliver.
      type(error_type), allocatable, intent(out) :: error
      type(local_BT_cont_u_type) :: BTC
      real(wp) :: uhbt_lo, uhbt_hi, slope_naive, slope_actual

      call make_realistic_btc(BTC)
      ! Sample two velocities deep in the saturated positive branch.
      uhbt_lo = find_uhbt(BTC%uBT_WW + 0.5_wp, BTC)
      uhbt_hi = find_uhbt(BTC%uBT_WW + 1.5_wp, BTC)
      slope_naive = BTC%FA_u_W0
      slope_actual = (uhbt_hi - uhbt_lo)/1.0_wp
      call check(error, slope_actual < slope_naive, &
                 "saturated transport slope not below naive slope")
      if (allocated(error)) return
      call check(error, abs(slope_actual - BTC%FA_u_WW) < 1.0e-10_wp, &
                 "saturated transport slope ≠ FA_u_WW")
   end subroutine test_saturation_cap

end module test_ocean_bt_cont_type
