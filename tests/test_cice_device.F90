!! Device exercise of the `!$acc routine seq` helpers that take the constants
!! bundle `type(ice_constants_t), intent(in)`: enthalpy <-> temperature (ice
!! and snow) and the freezing point, called from `do concurrent` over mapped
!! arrays and compared with the host results.
!!
!! On the GPU build (`mem:separate`) every array the kernel touches is mapped
!! explicitly.  Tolerance, not bitwise: the device build contracts
!! multiply-adds.  On host builds the directives are inert and the test still
!! checks that the kernel form equals the scalar form.
module test_cice_device
   use, intrinsic :: iso_fortran_env, only: int64
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_ice_params, only: ice_constants_t, ice_constants, ICE_CONSTANTS_ICEPACK, ICE_CONSTANTS_CESM, &
                             ICE_TFRZ_LINEAR_SALT, ice_q_from_T, ice_T_from_q, ice_qsno_from_T, &
                             ice_Tsno_from_q, ice_sea_freezing_temperature
   use test_cice_util, only: ice_random_columns, ice_golden_first_bad
   implicit none
   private

   public :: collect_cice_device_tests

   integer, parameter :: N = 4096, NL = 1

contains

   subroutine collect_cice_device_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("constants_bundle_in_kernel_icepack", test_icepack), &
                  new_unittest("constants_bundle_in_kernel_cesm", test_cesm) &
                  ]
   end subroutine collect_cice_device_tests

   subroutine test_icepack(error)
      type(error_type), allocatable, intent(out) :: error

      call run_set(error, ICE_CONSTANTS_ICEPACK)
   end subroutine test_icepack

   subroutine test_cesm(error)
      type(error_type), allocatable, intent(out) :: error

      call run_set(error, ICE_CONSTANTS_CESM)
   end subroutine test_cesm

   subroutine run_set(error, set)
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: set
      type(ice_constants_t) :: c
      real(wp) :: hi(N), hs(N), Tsfc(N), S(NL, N), T(NL, N)
      real(wp) :: Ti(N), Si(N), Tsn(N)
      real(wp) :: q_d(N), Tb_d(N), qs_d(N), Ts_d(N), tf_d(N)
      real(wp) :: q_h(N), Tb_h(N), qs_h(N), Ts_h(N), tf_h(N)
      integer :: bad_d(N), bad_h(N), i, nbad
      logical :: b, b2

      c = ice_constants(set)
      call ice_random_columns(N, NL, 424242_int64, hi, hs, Tsfc, S, T)
      Ti = T(1, :)
      Si = S(1, :)
      Tsn = Tsfc
      ! Host reference (scalar calls).
      do i = 1, N
         q_h(i) = ice_q_from_T(Ti(i), Si(i), c, .true.)
         call ice_T_from_q(q_h(i), Si(i), c, .true., Tb_h(i), b)
         bad_h(i) = merge(1, 0, b)
         qs_h(i) = ice_qsno_from_T(Tsn(i), c)
         Ts_h(i) = ice_Tsno_from_q(qs_h(i), c)
         call ice_sea_freezing_temperature(Si(i), ICE_TFRZ_LINEAR_SALT, c%depressT, c%Tocnfrz, tf_h(i), b2)
      end do

      q_d = 0.0_wp
      Tb_d = 0.0_wp
      qs_d = 0.0_wp
      Ts_d = 0.0_wp
      tf_d = 0.0_wp
      bad_d = -1
      !$acc enter data copyin(Ti, Si, Tsn) create(q_d, Tb_d, qs_d, Ts_d, tf_d, bad_d)
      do concurrent(i=1:N) local(b, b2)
         q_d(i) = ice_q_from_T(Ti(i), Si(i), c, .true.)
         call ice_T_from_q(q_d(i), Si(i), c, .true., Tb_d(i), b)
         bad_d(i) = merge(1, 0, b)
         qs_d(i) = ice_qsno_from_T(Tsn(i), c)
         Ts_d(i) = ice_Tsno_from_q(qs_d(i), c)
         call ice_sea_freezing_temperature(Si(i), ICE_TFRZ_LINEAR_SALT, c%depressT, c%Tocnfrz, tf_d(i), b2)
      end do
      !$acc update self(q_d, Tb_d, qs_d, Ts_d, tf_d, bad_d)
      !$acc exit data delete(Ti, Si, Tsn, q_d, Tb_d, qs_d, Ts_d, tf_d, bad_d)

      nbad = sum(bad_d)
      call check(error, nbad == 0 .and. sum(bad_h) == 0, "no flags raised")
      if (allocated(error)) return
      call check(error, ice_golden_first_bad(q_d, q_h, 1.0e-13_wp, 0.0_wp) == 0, "q(T,S) device == host")
      if (allocated(error)) return
      call check(error, ice_golden_first_bad(Tb_d, Tb_h, 0.0_wp, 1.0e-12_wp) == 0, "T(q,S) device == host")
      if (allocated(error)) return
      call check(error, ice_golden_first_bad(qs_d, qs_h, 1.0e-13_wp, 0.0_wp) == 0, "snow q device == host")
      if (allocated(error)) return
      call check(error, ice_golden_first_bad(Ts_d, Ts_h, 0.0_wp, 1.0e-12_wp) == 0, "snow T device == host")
      if (allocated(error)) return
      call check(error, ice_golden_first_bad(tf_d, tf_h, 1.0e-14_wp, 0.0_wp) == 0, "Tf device == host")
   end subroutine run_set

end module test_cice_device
