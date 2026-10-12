module test_ocean_teos10_convert
   !! Tests for the pure TEOS-10 conversion module `rdb_teos10_convert`.
   !!
   !! Check values extracted ONCE (python3.10 + netCDF4, never pip) from
   !! the third-party TEOS-10/GSW toolbox's own reference-check dataset
   !! `gsw_data_v3_0.nc` (`SA/SP/CT/p_chck_cast` + the `*_from_*`/
   !! `*_freezing_poly` result arrays it ships, sat_frac=0.5 as its own
   !! check driver uses) into the `parameter` tables below -- 11 cast
   !! points spanning low/high salinity, warm/cold, shallow/deep. The
   !! values are the dataset's AS SHIPPED (`gsw_data_v3_0.nc`, unmodified,
   !! not recomputed here), so they are an independent reference.
   !!
   !! Tolerances: each is a fixed bound set ABOVE the GSW dataset's own
   !! `computation_accuracy` attribute for that function (the toolbox's
   !! declared reproducibility across platforms), with headroom for
   !! FMA/reassociation differences between gfortran, nvfortran and ifx.
   !! Measured errors on gfortran are ~1e-15 (2e-14 for SR).
   !! Refs: IOC/SCOR/IAPSO (2010) TEOS-10 manual; McDougall & Barker (2011).
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_constants, only: wp
   use rdb_teos10_convert, only: sr_from_sp, sp_from_sr, ct_from_pt, pt_from_ct, &
                                 ct_freezing_poly, pt_freezing_poly, ct_from_pt_impl, &
                                 pt_from_ct_impl, ct_freezing_poly_impl
   use testdrive, only: new_unittest, unittest_type, error_type, check
   implicit none
   private

   public :: collect_ocean_teos10_convert_tests

   integer, parameter :: N = 11
   real(wp), parameter :: SAT_FRAC = 0.5_wp  !! GSW check driver's own value.

   real(wp), parameter :: SA(N) = [ &
                          34.50749946569206_wp, &
                          34.6075463959262_wp, &
                          34.751452325212476_wp, &
                          34.86361922550161_wp, &
                          34.60569590953979_wp, &
                          34.800921106392714_wp, &
                          34.79356770686205_wp, &
                          34.879244950054655_wp, &
                          10.389468455026284_wp, &
                          34.60668764696172_wp, &
                          6.669904340924572_wp]
   real(wp), parameter :: SP(N) = [ &
                          34.345367_wp, &
                          34.442738000000006_wp, &
                          34.574997_wp, &
                          34.683853000000006_wp, &
                          34.443054000000004_wp, &
                          34.633102_wp, &
                          34.614994_wp, &
                          34.700863000000005_wp, &
                          10.279548_wp, &
                          34.444057_wp, &
                          6.568259000000001_wp]
   real(wp), parameter :: CT(N) = [ &
                          27.994827331978655_wp, &
                          12.005103039522865_wp, &
                          3.2883333327081914_wp, &
                          1.1343476986712904_wp, &
                          27.279898123599327_wp, &
                          8.968233107420147_wp, &
                          2.4970408562072413_wp, &
                          0.8870403114493046_wp, &
                          4.61139799279878_wp, &
                          27.320893180033327_wp, &
                          10.50276773530223_wp]
   real(wp), parameter :: P(N) = [ &
                          0.0_wp, &
                          252.0_wp, &
                          1314.0_wp, &
                          4069.0_wp, &
                          30.0_wp, &
                          404.0_wp, &
                          1771.0_wp, &
                          4840.0_wp, &
                          101.0_wp, &
                          0.0_wp, &
                          0.0_wp]
   real(wp), parameter :: SR_FROM_SP_REF(N) = [ &
                          34.50732012484802_wp, &
                          34.60515027084345_wp, &
                          34.73803292871088_wp, &
                          34.84740223140345_wp, &
                          34.60546776091888_wp, &
                          34.796411918688015_wp, &
                          34.778218531707445_wp, &
                          34.86449244084345_wp, &
                          10.328020474340576_wp, &
                          34.60647549049373_wp, &
                          6.599231156153147_wp]
   real(wp), parameter :: PT_FROM_CT_REF(N) = [ &
                          27.962_wp, &
                          12.004161576491267_wp, &
                          3.2888557950200528_wp, &
                          1.134061398239781_wp, &
                          27.253077797851773_wp, &
                          8.97167513742059_wp, &
                          2.4973809758874443_wp, &
                          0.886620221881316_wp, &
                          4.409035820643079_wp, &
                          27.294_wp, &
                          10.045999999999994_wp]
   real(wp), parameter :: CT_FROM_PT_REF(N) = [ &
                          27.994827331978655_wp, &
                          12.005103039522865_wp, &
                          3.2883333327081914_wp, &
                          1.1343476986712904_wp, &
                          27.279898123599327_wp, &
                          8.968233107420147_wp, &
                          2.4970408562072413_wp, &
                          0.8870403114493046_wp, &
                          4.61139799279878_wp, &
                          27.320893180033327_wp, &
                          10.50276773530223_wp]
   real(wp), parameter :: CT_FREEZING_REF(N) = [ &
                          -1.8791652991889622_wp, &
                          -2.079029915212826_wp, &
                          -2.9314243188431504_wp, &
                          -5.318308325801521_wp, &
                          -1.9078508208601117_wp, &
                          -2.2085298703693867_wp, &
                          -3.309927008868304_wp, &
                          -6.032822204520544_wp, &
                          -0.6176071419076397_wp, &
                          -1.8849339490842347_wp, &
                          -0.34436067697948836_wp]
   real(wp), parameter :: T_FREEZING_REF(N) = [ &
                          -1.8826599878941543_wp, &
                          -2.077968097708414_wp, &
                          -2.907274910075863_wp, &
                          -5.207621634554549_wp, &
                          -1.9107060315809388_wp, &
                          -2.204319649514458_wp, &
                          -3.2738997935678396_wp, &
                          -5.892403110925523_wp, &
                          -0.6333932622374288_wp, &
                          -1.8882928863900617_wp, &
                          -0.35965171025197_wp]

contains

   subroutine collect_ocean_teos10_convert_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("teos10_sr_from_sp", test_sr_from_sp), &
                  new_unittest("teos10_zero_salinity_finite", test_zero_salinity_finite), &
                  new_unittest("teos10_sp_from_sr_inverse", test_sp_from_sr_inverse), &
                  new_unittest("teos10_pt_from_ct", test_pt_from_ct), &
                  new_unittest("teos10_ct_from_pt", test_ct_from_pt), &
                  new_unittest("teos10_round_trip", test_round_trip), &
                  new_unittest("teos10_ct_freezing_poly", test_ct_freezing_poly), &
                  new_unittest("teos10_pt_freezing_poly_consistency", test_pt_freezing_poly_consistency), &
                  new_unittest("teos10_gpu_residency", test_gpu_residency) &
                  ]
   end subroutine collect_ocean_teos10_convert_tests

   subroutine test_sr_from_sp(error)
      !! |SR_from_SP - GSW| <= 1e-9 g/kg: GSW declares 1.3e-10.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: mx
      mx = maxval(abs(sr_from_sp(SP) - SR_FROM_SP_REF))
      call check(error, mx <= 1.0e-9_wp, "sr_from_sp vs GSW check cast")
   end subroutine test_sr_from_sp

   subroutine test_sp_from_sr_inverse(error)
      !! sp_from_sr is the exact algebraic inverse of sr_from_sp.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: mx
      mx = maxval(abs(sp_from_sr(sr_from_sp(SP)) - SP))
      call check(error, mx <= 1.0e-12_wp, "sp_from_sr(sr_from_sp(sp)) /= sp")
   end subroutine test_sp_from_sr_inverse

   subroutine test_pt_from_ct(error)
      !! |pt_from_ct(SA, CT) - GSW| <= 1e-8 degC: GSW declares 6.1e-10.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: mx
      mx = maxval(abs(pt_from_ct(SA, CT) - PT_FROM_CT_REF))
      call check(error, mx <= 1.0e-8_wp, "pt_from_ct vs GSW check cast")
   end subroutine test_pt_from_ct

   subroutine test_ct_from_pt(error)
      !! |ct_from_pt(SA, PT_FROM_CT_REF) - GSW| <= 1e-8 degC (GSW declares
      !! 6.3e-10) --
      !! PT_FROM_CT_REF is itself the GSW-reference pt_from_CT(SA,CT), so
      !! this exercises ct_from_pt directly against the matching
      !! reference CT_from_pt output, exactly as the GSW check driver
      !! wires the two functions together.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: mx
      mx = maxval(abs(ct_from_pt(SA, PT_FROM_CT_REF) - CT_FROM_PT_REF))
      call check(error, mx <= 1.0e-8_wp, "ct_from_pt vs GSW check cast")
   end subroutine test_ct_from_pt

   subroutine test_round_trip(error)
      !! ct_from_pt(SA, pt_from_ct(SA, CT)) == CT to <= 1e-10 degC: the
      !! 1.5-step Newton-McDougall inverse converges to ~1e-12 degC
      !! (McDougall & Barker 2011), so 1e-10 checks inversion, not GSW.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: mx
      mx = maxval(abs(ct_from_pt(SA, pt_from_ct(SA, CT)) - CT))
      call check(error, mx <= 1.0e-10_wp, "ct_from_pt(pt_from_ct(CT)) round trip")
   end subroutine test_round_trip

   subroutine test_ct_freezing_poly(error)
      !! |ct_freezing_poly(SA, p[dbar], 0.5) - GSW| <= 1e-8 degC: GSW
      !! declares 2.3e-11.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: mx
      mx = maxval(abs(ct_freezing_poly(SA, P, SAT_FRAC) - CT_FREEZING_REF))
      call check(error, mx <= 1.0e-8_wp, "ct_freezing_poly vs GSW check cast")
   end subroutine test_ct_freezing_poly

   subroutine test_pt_freezing_poly_consistency(error)
      !! At p=0, in-situ temperature == potential temperature, so GSW's
      !! own t_freezing_poly(SA,0,sat_frac) check value doubles as the
      !! independent reference for pt_freezing_poly(SA,0,sat_frac) here
      !! (<= 1e-8 degC; GSW declares 2.2e-11).  Three of the eleven cast points have
      !! p=0 (indices 1, 10, 11).
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: mx
      real(wp) :: pt_f(3)
      integer, parameter :: IDX(3) = [1, 10, 11]
      pt_f = pt_freezing_poly(SA(IDX), 0.0_wp, SAT_FRAC)
      mx = maxval(abs(pt_f - T_FREEZING_REF(IDX)))
      call check(error, mx <= 1.0e-8_wp, "pt_freezing_poly vs t_freezing_poly at p=0")
   end subroutine test_pt_freezing_poly_consistency

   subroutine test_zero_salinity_finite(error)
      !! Fresh water (SA = SP = 0) is a legal input (rivers, ice melt): every
      !! primitive returns a finite value there, and CT/pt round-trip.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: ct, pt, ctf, ptf
      ct = ct_from_pt(0.0_wp, 10.0_wp)
      pt = pt_from_ct(0.0_wp, ct)
      ctf = ct_freezing_poly(0.0_wp, 0.0_wp, SAT_FRAC)
      ptf = pt_freezing_poly(0.0_wp, 0.0_wp, SAT_FRAC)
      call check(error, ieee_is_finite(ct) .and. ieee_is_finite(pt) .and. &
                 ieee_is_finite(ctf) .and. ieee_is_finite(ptf) .and. &
                 ieee_is_finite(sr_from_sp(0.0_wp)), "non-finite result at zero salinity")
      if (allocated(error)) return
      call check(error, abs(pt - 10.0_wp) <= 1.0e-10_wp, "zero-salinity CT/pt round trip")
   end subroutine test_zero_salinity_finite

   subroutine test_gpu_residency(error)
      !! `mem:separate` contract (see test_open_boundary_out_closes in
      !! test_ocean_conservation_salt_heat.F90): every array a `do
      !! concurrent` kernel touches must be explicitly mapped -- inputs
      !! `copyin`, pure kernel OUTPUTS `create` (never read before the
      !! kernel writes them, so no host-value loss), pulled back with
      !! `update self` before the host comparison.  Inert directives on
      !! host/multicore builds; exercises the actual device code path
      !! on the GPU build.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: M = 6
      real(wp) :: sa_l(M), pt_l(M), p_l(M), sat_l(M)
      real(wp) :: ct_dev(M), pt_dev(M), ctf_dev(M)
      real(wp) :: ct_host(M), pt_host(M), ctf_host(M)
      integer :: i

      sa_l = SA(1:M)
      pt_l = PT_FROM_CT_REF(1:M)
      p_l = P(1:M)
      sat_l = SAT_FRAC

      !$acc enter data copyin(sa_l, pt_l, p_l, sat_l) create(ct_dev, pt_dev, ctf_dev)
      ! Library kernels: a test's own kernel cannot device-link a library
      ! `routine seq` function.
      call ct_from_pt_impl(M, sa_l, pt_l, ct_dev)
      call pt_from_ct_impl(M, sa_l, ct_dev, pt_dev)
      call ct_freezing_poly_impl(M, sa_l, p_l, sat_l, ctf_dev)
      !$acc update self(ct_dev, pt_dev, ctf_dev)
      !$acc exit data delete(sa_l, pt_l, p_l, sat_l, ct_dev, pt_dev, ctf_dev)

      do i = 1, M
         ct_host(i) = ct_from_pt(sa_l(i), pt_l(i))
         pt_host(i) = pt_from_ct(sa_l(i), ct_host(i))
         ctf_host(i) = ct_freezing_poly(sa_l(i), p_l(i), sat_l(i))
      end do

      call check(error, maxval(abs(ct_dev - ct_host)) < 1.0e-12_wp, &
                 "ct_from_pt device result != host result")
      if (allocated(error)) return
      call check(error, maxval(abs(pt_dev - pt_host)) < 1.0e-12_wp, &
                 "pt_from_ct device result != host result")
      if (allocated(error)) return
      call check(error, maxval(abs(ctf_dev - ctf_host)) < 1.0e-12_wp, &
                 "ct_freezing_poly device result != host result")
   end subroutine test_gpu_residency

end module test_ocean_teos10_convert
