!! Golden comparisons of the sea-ice constants, salinity profile, freezing
!! point and enthalpy conversions against reference tables from the Icepack
!! column physics (both constant sets).
!!
!! The tables are embedded as 17-significant-digit literals; they were produced
!! by the Icepack oracle build (upstream routines called directly, except the
!! forward enthalpy, whose formula is the one in the upstream initialisation).
!! The CESM-set values rest on SHR_CONST numbers from an old CAM constants
!! file and are unverified against the CESM version ACCESS-NRI builds
!! (interim choice).
!!
!! Tolerances.  Constants: relative 1e-15 (derived Lfresh = Lsub - Lvap and
!! Tocnfrz = -34*0.054 may differ by one rounding).  Freezing point: exact
!! multiplication, relative 1e-15.  Forward enthalpy: same operations in the
!! same order, relative 1e-14 (headroom for a fused multiply-add on the GPU
!! build).  Inverse temperature: the quadratic root cancels (relative
!! conditioning ~ |q/(rhoi cp_ice T)|), so 1e-12 relative.  Salinity profile:
!! the reproducible library cos/pow differ from the system libm by ulps,
!! relative 1e-13.
module test_cice_golden
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_ice_params, only: ice_constants_t, ice_constants, ICE_CONSTANTS_ICEPACK, ICE_CONSTANTS_CESM, &
                             ICE_TFRZ_LINEAR_SALT, ice_params_default, ice_q_from_T, ice_T_from_q, &
                             ice_qsno_init_from_T, ice_Tsno_from_q, ice_salinity_profile, ice_salinity_at, &
                             ice_sea_freezing_temperature, ice_params_t
   use test_cice_util, only: ice_golden_check
   implicit none
   private

   public :: collect_cice_golden_tests

   real(wp), parameter :: K_ICEPACK(18) = [ &
                          9.1700000000000000E+002_wp, &
                          3.3000000000000000E+002_wp, &
                          1.0260000000000000E+003_wp, &
                          2.1060000000000000E+003_wp, &
                          4.2180000000000000E+003_wp, &
                          3.3400000000000000E+005_wp, &
                          2.8350000000000000E+006_wp, &
                          2.5010000000000000E+006_wp, &
                          5.3999999999999999E-002_wp, &
                          2.7314999999999998E+002_wp, &
                          9.8061600000000002E+000_wp, &
                          5.6699999999999998E-008_wp, &
                          2.0000000000000000E-002_wp, &
                          9.9999999999999994E-012_wp, &
                          4.0000000000000000E+000_wp, &
                          -1.8000000000000000E+000_wp, &
                          3.2000000000000002E+000_wp, &
                          1.0000000000000001E-001_wp &
                          ]
   real(wp), parameter :: EI_T_ICEPACK(15) = [ &
                          -4.0000000000000000E+001_wp, &
                          -1.0000000000000000E+001_wp, &
                          -2.0000000000000000E+000_wp, &
                          -5.0000000000000000E-001_wp, &
                          -4.0000000000000000E+001_wp, &
                          -1.0000000000000000E+001_wp, &
                          -2.0000000000000000E+000_wp, &
                          -5.0000000000000000E-001_wp, &
                          -4.0000000000000000E+001_wp, &
                          -1.0000000000000000E+001_wp, &
                          -2.0000000000000000E+000_wp, &
                          -5.0000000000000000E-001_wp, &
                          -4.0000000000000000E+001_wp, &
                          -1.0000000000000000E+001_wp, &
                          -2.0000000000000000E+000_wp &
                          ]
   real(wp), parameter :: EI_S_ICEPACK(15) = [ &
                          1.0000000000000001E-001_wp, &
                          1.0000000000000001E-001_wp, &
                          1.0000000000000001E-001_wp, &
                          1.0000000000000001E-001_wp, &
                          3.2000000000000002E+000_wp, &
                          3.2000000000000002E+000_wp, &
                          3.2000000000000002E+000_wp, &
                          3.2000000000000002E+000_wp, &
                          5.0000000000000000E+000_wp, &
                          5.0000000000000000E+000_wp, &
                          5.0000000000000000E+000_wp, &
                          5.0000000000000000E+000_wp, &
                          1.2000000000000000E+001_wp, &
                          1.2000000000000000E+001_wp, &
                          1.2000000000000000E+001_wp &
                          ]
   real(wp), parameter :: EI_Q_ICEPACK(15) = [ &
                          -3.8349519067159998E+008_wp, &
                          -3.2543508808160001E+008_wp, &
                          -3.0932391160160005E+008_wp, &
                          -3.0394625680159998E+008_wp, &
                          -3.8253762149120003E+008_wp, &
                          -3.2063219861120003E+008_wp, &
                          -2.8401264725119996E+008_wp, &
                          -2.0172858665120003E+008_wp, &
                          -3.8198161357999998E+008_wp, &
                          -3.1784342407999998E+008_wp, &
                          -2.6931578407999998E+008_wp, &
                          -1.4237639107999998E+008_wp, &
                          -3.7981936059200001E+008_wp, &
                          -3.0699818979200000E+008_wp, &
                          -2.1216131619199997E+008_wp &
                          ]
   real(wp), parameter :: EI_TB_ICEPACK(15) = [ &
                          -3.9999999999999986E+001_wp, &
                          -1.0000000000000007E+001_wp, &
                          -2.0000000000000044E+000_wp, &
                          -4.9999999999999245E-001_wp, &
                          -4.0000000000000028E+001_wp, &
                          -1.0000000000000014E+001_wp, &
                          -1.9999999999999982E+000_wp, &
                          -5.0000000000000000E-001_wp, &
                          -4.0000000000000000E+001_wp, &
                          -1.0000000000000000E+001_wp, &
                          -2.0000000000000000E+000_wp, &
                          -5.0000000000000000E-001_wp, &
                          -4.0000000000000007E+001_wp, &
                          -9.9999999999999964E+000_wp, &
                          -2.0000000000000000E+000_wp &
                          ]
   real(wp), parameter :: ES_T_ICEPACK(5) = [ &
                          -4.0000000000000000E+001_wp, &
                          -1.0000000000000000E+001_wp, &
                          -2.0000000000000000E+000_wp, &
                          -5.0000000000000000E-001_wp, &
                          -1.0000000000000001E-001_wp &
                          ]
   real(wp), parameter :: ES_Q_ICEPACK(5) = [ &
                          -1.3801920000000000E+008_wp, &
                          -1.1716980000000000E+008_wp, &
                          -1.1160996000000000E+008_wp, &
                          -1.1056749000000000E+008_wp, &
                          -1.1028949799999999E+008_wp &
                          ]
   real(wp), parameter :: ES_TB_ICEPACK(5) = [ &
                          -4.0000000000000000E+001_wp, &
                          -1.0000000000000000E+001_wp, &
                          -2.0000000000000000E+000_wp, &
                          -5.0000000000000000E-001_wp, &
                          -9.9999999999988945E-002_wp &
                          ]
   real(wp), parameter :: FZ_S_ICEPACK(15) = [ &
                          0.0000000000000000E+000_wp, &
                          1.0000000000000001E-001_wp, &
                          1.0000000000000000E+000_wp, &
                          3.2000000000000002E+000_wp, &
                          4.0000000000000000E+000_wp, &
                          5.0000000000000000E+000_wp, &
                          1.0000000000000000E+001_wp, &
                          2.0000000000000000E+001_wp, &
                          3.0000000000000000E+001_wp, &
                          3.2000000000000000E+001_wp, &
                          3.4000000000000000E+001_wp, &
                          3.4700000000000003E+001_wp, &
                          3.5000000000000000E+001_wp, &
                          3.6000000000000000E+001_wp, &
                          4.0000000000000000E+001_wp &
                          ]
   real(wp), parameter :: FZ_TF_ICEPACK(15) = [ &
                          -0.0000000000000000E+000_wp, &
                          -5.4000000000000003E-003_wp, &
                          -5.3999999999999999E-002_wp, &
                          -1.7280000000000001E-001_wp, &
                          -2.1600000000000000E-001_wp, &
                          -2.7000000000000002E-001_wp, &
                          -5.4000000000000004E-001_wp, &
                          -1.0800000000000001E+000_wp, &
                          -1.6199999999999999E+000_wp, &
                          -1.7280000000000000E+000_wp, &
                          -1.8360000000000001E+000_wp, &
                          -1.8738000000000001E+000_wp, &
                          -1.8899999999999999E+000_wp, &
                          -1.9440000000000000E+000_wp, &
                          -2.1600000000000001E+000_wp &
                          ]
   real(wp), parameter :: SP_ICEPACK(5) = [ &
                          6.4920185640031347E-001_wp, &
                          2.3545811119274354E+000_wp, &
                          3.0310920697984636E+000_wp, &
                          3.1892977722169995E+000_wp, &
                          3.2000000000000002E+000_wp &
                          ]
   real(wp), parameter :: SF_ZN(20) = [ &
                          5.0000000000000003E-002_wp, &
                          1.0000000000000001E-001_wp, &
                          1.4999999999999999E-001_wp, &
                          2.0000000000000001E-001_wp, &
                          2.5000000000000000E-001_wp, &
                          2.9999999999999999E-001_wp, &
                          3.4999999999999998E-001_wp, &
                          4.0000000000000002E-001_wp, &
                          4.5000000000000001E-001_wp, &
                          5.0000000000000000E-001_wp, &
                          5.5000000000000004E-001_wp, &
                          5.9999999999999998E-001_wp, &
                          6.5000000000000002E-001_wp, &
                          6.9999999999999996E-001_wp, &
                          7.5000000000000000E-001_wp, &
                          8.0000000000000004E-001_wp, &
                          8.4999999999999998E-001_wp, &
                          9.0000000000000002E-001_wp, &
                          9.4999999999999996E-001_wp, &
                          1.0000000000000000E+000_wp &
                          ]
   real(wp), parameter :: SF_S(20) = [ &
                          1.5500592074108291E-001_wp, &
                          4.6314589457992561E-001_wp, &
                          8.4559860595624137E-001_wp, &
                          1.2437492048558576E+000_wp, &
                          1.6191224998882505E+000_wp, &
                          1.9513130405273573E+000_wp, &
                          2.2328583740964834E+000_wp, &
                          2.4641989425431818E+000_wp, &
                          2.6498978336212478E+000_wp, &
                          2.7961969999707410E+000_wp, &
                          2.9096063116676256E+000_wp, &
                          2.9961808128568856E+000_wp, &
                          3.0612141515423730E+000_wp, &
                          3.1091636022895912E+000_wp, &
                          3.1436921978536883E+000_wp, &
                          3.1677615982897374E+000_wp, &
                          3.1837396235436115E+000_wp, &
                          3.1935043544288795E+000_wp, &
                          3.1985368546216209E+000_wp, &
                          3.2000000000000002E+000_wp &
                          ]
   real(wp), parameter :: K_CESM(18) = [ &
                          9.1700000000000000E+002_wp, &
                          3.3000000000000000E+002_wp, &
                          1.0260000000000000E+003_wp, &
                          2.1172700000000000E+003_wp, &
                          3.9960000000000000E+003_wp, &
                          3.3370000000000000E+005_wp, &
                          2.8347000000000000E+006_wp, &
                          2.5010000000000000E+006_wp, &
                          5.3999999999999999E-002_wp, &
                          2.7314999999999998E+002_wp, &
                          9.8061600000000002E+000_wp, &
                          5.6699999999999998E-008_wp, &
                          5.0000000000000001E-003_wp, &
                          9.9999999999999994E-012_wp, &
                          4.0000000000000000E+000_wp, &
                          -1.8360000000000001E+000_wp, &
                          3.2000000000000002E+000_wp, &
                          1.0000000000000001E-001_wp &
                          ]
   real(wp), parameter :: EI_T_CESM(15) = [ &
                          -4.0000000000000000E+001_wp, &
                          -1.0000000000000000E+001_wp, &
                          -2.0000000000000000E+000_wp, &
                          -5.0000000000000000E-001_wp, &
                          -4.0000000000000000E+001_wp, &
                          -1.0000000000000000E+001_wp, &
                          -2.0000000000000000E+000_wp, &
                          -5.0000000000000000E-001_wp, &
                          -4.0000000000000000E+001_wp, &
                          -1.0000000000000000E+001_wp, &
                          -2.0000000000000000E+000_wp, &
                          -5.0000000000000000E-001_wp, &
                          -4.0000000000000000E+001_wp, &
                          -1.0000000000000000E+001_wp, &
                          -2.0000000000000000E+000_wp &
                          ]
   real(wp), parameter :: EI_S_CESM(15) = [ &
                          1.0000000000000001E-001_wp, &
                          1.0000000000000001E-001_wp, &
                          1.0000000000000001E-001_wp, &
                          1.0000000000000001E-001_wp, &
                          3.2000000000000002E+000_wp, &
                          3.2000000000000002E+000_wp, &
                          3.2000000000000002E+000_wp, &
                          3.2000000000000002E+000_wp, &
                          5.0000000000000000E+000_wp, &
                          5.0000000000000000E+000_wp, &
                          5.0000000000000000E+000_wp, &
                          5.0000000000000000E+000_wp, &
                          1.2000000000000000E+001_wp, &
                          1.2000000000000000E+001_wp, &
                          1.2000000000000000E+001_wp &
                          ]
   real(wp), parameter :: EI_Q_CESM(15) = [ &
                          -3.8363235630371398E+008_wp, &
                          -3.2526232742921400E+008_wp, &
                          -3.0906906844521397E+008_wp, &
                          -3.0367814007021397E+008_wp, &
                          -3.8264013011884797E+008_wp, &
                          -3.2042823483484799E+008_wp, &
                          -2.8374502166684800E+008_wp, &
                          -2.0151676510184801E+008_wp, &
                          -3.8206399878569996E+008_wp, &
                          -3.1762134236069995E+008_wp, &
                          -2.6904073644069999E+008_wp, &
                          -1.4219725705570003E+008_wp, &
                          -3.7982348804567999E+008_wp, &
                          -3.0670564940568000E+008_wp, &
                          -2.1185740500567999E+008_wp &
                          ]
   real(wp), parameter :: EI_TB_CESM(15) = [ &
                          -3.9999999999999979E+001_wp, &
                          -9.9999999999999893E+000_wp, &
                          -1.9999999999999953E+000_wp, &
                          -4.9999999999999883E-001_wp, &
                          -3.9999999999999979E+001_wp, &
                          -1.0000000000000000E+001_wp, &
                          -2.0000000000000040E+000_wp, &
                          -5.0000000000000100E-001_wp, &
                          -3.9999999999999979E+001_wp, &
                          -9.9999999999999947E+000_wp, &
                          -1.9999999999999987E+000_wp, &
                          -4.9999999999999756E-001_wp, &
                          -3.9999999999999993E+001_wp, &
                          -1.0000000000000000E+001_wp, &
                          -2.0000000000000004E+000_wp &
                          ]
   real(wp), parameter :: ES_T_CESM(5) = [ &
                          -4.0000000000000000E+001_wp, &
                          -1.0000000000000000E+001_wp, &
                          -2.0000000000000000E+000_wp, &
                          -5.0000000000000000E-001_wp, &
                          -1.0000000000000001E-001_wp &
                          ]
   real(wp), parameter :: ES_Q_CESM(5) = [ &
                          -1.3806896400000000E+008_wp, &
                          -1.1710799100000000E+008_wp, &
                          -1.1151839819999999E+008_wp, &
                          -1.1047034955000000E+008_wp, &
                          -1.1019086991000001E+008_wp &
                          ]
   real(wp), parameter :: ES_TB_CESM(5) = [ &
                          -3.9999999999999993E+001_wp, &
                          -1.0000000000000005E+001_wp, &
                          -1.9999999999999902E+000_wp, &
                          -5.0000000000000444E-001_wp, &
                          -1.0000000000000638E-001_wp &
                          ]
   real(wp), parameter :: FZ_S_CESM(15) = [ &
                          0.0000000000000000E+000_wp, &
                          1.0000000000000001E-001_wp, &
                          1.0000000000000000E+000_wp, &
                          3.2000000000000002E+000_wp, &
                          4.0000000000000000E+000_wp, &
                          5.0000000000000000E+000_wp, &
                          1.0000000000000000E+001_wp, &
                          2.0000000000000000E+001_wp, &
                          3.0000000000000000E+001_wp, &
                          3.2000000000000000E+001_wp, &
                          3.4000000000000000E+001_wp, &
                          3.4700000000000003E+001_wp, &
                          3.5000000000000000E+001_wp, &
                          3.6000000000000000E+001_wp, &
                          4.0000000000000000E+001_wp &
                          ]
   real(wp), parameter :: FZ_TF_CESM(15) = [ &
                          -0.0000000000000000E+000_wp, &
                          -5.4000000000000003E-003_wp, &
                          -5.3999999999999999E-002_wp, &
                          -1.7280000000000001E-001_wp, &
                          -2.1600000000000000E-001_wp, &
                          -2.7000000000000002E-001_wp, &
                          -5.4000000000000004E-001_wp, &
                          -1.0800000000000001E+000_wp, &
                          -1.6199999999999999E+000_wp, &
                          -1.7280000000000000E+000_wp, &
                          -1.8360000000000001E+000_wp, &
                          -1.8738000000000001E+000_wp, &
                          -1.8899999999999999E+000_wp, &
                          -1.9440000000000000E+000_wp, &
                          -2.1600000000000001E+000_wp &
                          ]
   real(wp), parameter :: SP_CESM(5) = [ &
                          6.4920185640031347E-001_wp, &
                          2.3545811119274354E+000_wp, &
                          3.0310920697984636E+000_wp, &
                          3.1892977722169995E+000_wp, &
                          3.2000000000000002E+000_wp &
                          ]

contains

   subroutine collect_cice_golden_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("golden_constants", test_constants), &
                  new_unittest("golden_freezing_point", test_freezing), &
                  new_unittest("golden_salinity_profile", test_salinity), &
                  new_unittest("golden_enthalpy_ice", test_enthalpy_ice), &
                  new_unittest("golden_enthalpy_snow", test_enthalpy_snow) &
                  ]
   end subroutine collect_cice_golden_tests

   function const_vector(p) result(v)
      type(ice_params_t), intent(in) :: p
      type(ice_constants_t) :: c
      real(wp) :: v(18)

      c = p%c

      v = [c%rhoi, c%rhos, c%rhow, c%cp_ice, c%cp_ocn, c%Lfresh, c%Lsub, c%Lvap, c%depressT, &
           c%Tffresh, c%gravit, c%stefan_boltzmann, c%snowpatch, c%puny, p%thermo%ice_ref_salinity, c%Tocnfrz, &
           c%saltmax, c%min_salin]
   end function const_vector

   subroutine test_constants(error)
      type(error_type), allocatable, intent(out) :: error

      call ice_golden_check(error, "constants, upstream set", &
                            const_vector(ice_params_default(ICE_CONSTANTS_ICEPACK, 5, 4, 1)), K_ICEPACK, 1.0e-15_wp, 0.0_wp)
      if (allocated(error)) return
      call ice_golden_check(error, "constants, cesm set", &
                            const_vector(ice_params_default(ICE_CONSTANTS_CESM, 5, 4, 1)), K_CESM, 1.0e-15_wp, 0.0_wp)
   end subroutine test_constants

   subroutine test_freezing(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: tf(15)
      type(ice_constants_t) :: c
      logical :: bad
      integer :: i

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      do i = 1, 15
         call ice_sea_freezing_temperature(FZ_S_ICEPACK(i), ICE_TFRZ_LINEAR_SALT, c%depressT, c%Tocnfrz, tf(i), bad)
      end do
      call ice_golden_check(error, "linear_salt Tf, upstream set", tf, FZ_TF_ICEPACK, 1.0e-15_wp, 0.0_wp)
      if (allocated(error)) return
      c = ice_constants(ICE_CONSTANTS_CESM)
      do i = 1, 15
         call ice_sea_freezing_temperature(FZ_S_CESM(i), ICE_TFRZ_LINEAR_SALT, c%depressT, c%Tocnfrz, tf(i), bad)
      end do
      call ice_golden_check(error, "linear_salt Tf, cesm set", tf, FZ_TF_CESM, 1.0e-15_wp, 0.0_wp)
   end subroutine test_freezing

   subroutine test_salinity(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: s(5), sf(size(SF_ZN))

      call ice_salinity_profile(4, 3.2_wp, 0.1_wp, s)
      call ice_golden_check(error, "4-layer profile (upstream set)", s, SP_ICEPACK, 1.0e-13_wp, 0.0_wp)
      if (allocated(error)) return
      call ice_salinity_profile(4, 3.2_wp, 0.1_wp, s)
      call ice_golden_check(error, "4-layer profile (cesm set)", s, SP_CESM, 1.0e-13_wp, 0.0_wp)
      if (allocated(error)) return
      sf = ice_salinity_at(SF_ZN, 3.2_wp)
      call ice_golden_check(error, "raw profile function", sf, SF_S, 1.0e-13_wp, 1.0e-15_wp)
   end subroutine test_salinity

   subroutine test_enthalpy_ice(error)
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      real(wp) :: q(15), tb(15)
      logical :: bad
      integer :: i

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      do i = 1, 15
         q(i) = ice_q_from_T(EI_T_ICEPACK(i), EI_S_ICEPACK(i), c, .true.)
      end do
      call ice_golden_check(error, "q(T,S) upstream set", q, EI_Q_ICEPACK, 1.0e-14_wp, 0.0_wp)
      if (allocated(error)) return
      do i = 1, 15
         call ice_T_from_q(EI_Q_ICEPACK(i), EI_S_ICEPACK(i), c, .true., tb(i), bad)
         call check(error,.not. bad, "no flag")
         if (allocated(error)) return
      end do
      call ice_golden_check(error, "T(q,S) upstream set", tb, EI_TB_ICEPACK, 1.0e-12_wp, 0.0_wp)
      if (allocated(error)) return

      c = ice_constants(ICE_CONSTANTS_CESM)
      do i = 1, 15
         q(i) = ice_q_from_T(EI_T_CESM(i), EI_S_CESM(i), c, .true.)
      end do
      call ice_golden_check(error, "q(T,S) cesm set", q, EI_Q_CESM, 1.0e-14_wp, 0.0_wp)
      if (allocated(error)) return
      do i = 1, 15
         call ice_T_from_q(EI_Q_CESM(i), EI_S_CESM(i), c, .true., tb(i), bad)
         call check(error,.not. bad, "no flag")
         if (allocated(error)) return
      end do
      call ice_golden_check(error, "T(q,S) cesm set", tb, EI_TB_CESM, 1.0e-12_wp, 0.0_wp)
   end subroutine test_enthalpy_ice

   subroutine test_enthalpy_snow(error)
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      real(wp) :: q(5), tb(5)
      integer :: i

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      do i = 1, 5
         q(i) = ice_qsno_init_from_T(ES_T_ICEPACK(i), c)
         tb(i) = ice_Tsno_from_q(ES_Q_ICEPACK(i), c)
      end do
      call ice_golden_check(error, "qsno(T) upstream set", q, ES_Q_ICEPACK, 1.0e-14_wp, 0.0_wp)
      if (allocated(error)) return
      call ice_golden_check(error, "Tsno(q) upstream set", tb, ES_TB_ICEPACK, 1.0e-12_wp, 1.0e-12_wp)
      if (allocated(error)) return
      c = ice_constants(ICE_CONSTANTS_CESM)
      do i = 1, 5
         q(i) = ice_qsno_init_from_T(ES_T_CESM(i), c)
         tb(i) = ice_Tsno_from_q(ES_Q_CESM(i), c)
      end do
      call ice_golden_check(error, "qsno(T) cesm set", q, ES_Q_CESM, 1.0e-14_wp, 0.0_wp)
      if (allocated(error)) return
      call ice_golden_check(error, "Tsno(q) cesm set", tb, ES_TB_CESM, 1.0e-12_wp, 1.0e-12_wp)
   end subroutine test_enthalpy_snow

end module test_cice_golden
