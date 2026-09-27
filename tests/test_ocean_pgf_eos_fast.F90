!! Tests for the fast EOS entry points.
!!
!! value_and_fused_entry_points -- `roquet_spv_value` is the `sv` of
!! `roquet_spv_point`, and `eos_density_specvol_derivs` returns what
!! `eos_density_point` + `eos_specvol_derivs` return, for every variant.
module test_ocean_pgf_eos_fast
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_eos, only: eos_t, eos_density_point, eos_specvol_derivs, &
                      eos_density_specvol_derivs, roquet_spv_value, &
                      EOS_VARIANT_LINEAR, EOS_VARIANT_WRIGHT_97, EOS_VARIANT_ROQUET_SPV
   use pic_logger, only: global_logger
   implicit none
   private

   public :: collect_ocean_pgf_eos_fast_tests

   real(wp), parameter :: RHO0 = 1035.0_wp
   integer, parameter :: N_T = 5, N_S = 3
   real(wp), parameter :: T_SET(N_T) = [-1.9_wp, 2.0_wp, 10.0_wp, 20.0_wp, 30.0_wp]
      !! Potential temperatures (degC): freezing to tropical surface.
   real(wp), parameter :: S_SET(N_S) = [32.0_wp, 34.7_wp, 37.0_wp]
      !! Practical salinities (PSU).
   real(wp), parameter :: TOL_ULP = 4.0_wp*epsilon(1.0_wp)
      !! Relative bound for "the same expression, compiled twice" -- equal
      !! but for an FMA contraction choice.

contains

   subroutine collect_ocean_pgf_eos_fast_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("value_and_fused_entry_points", test_entry_points) &
                  ]
   end subroutine collect_ocean_pgf_eos_fast_tests

   subroutine make_eos(eos, variant)
      type(eos_t), intent(out) :: eos
      integer, intent(in) :: variant
      eos%variant = variant
      eos%rho0 = RHO0
      eos%alpha_T = 0.17_wp
      eos%beta_S = 0.78_wp
      eos%is_init = .true.
   end subroutine make_eos

   subroutine test_entry_points(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: VARS(3) = [EOS_VARIANT_LINEAR, EOS_VARIANT_WRIGHT_97, &
                                       EOS_VARIANT_ROQUET_SPV]
      real(wp), parameter :: P_SET(4) = [0.0_wp, 1.0e6_wp, 2.0e7_wp, 6.0e7_wp]
      type(eos_t) :: eos
      real(wp) :: sv_v, sv_p, d1, d2, rho_a, rho_b, dt_a, ds_a, dt_b, ds_b, err, max_err
      integer :: it, is, ip, iv
      character(len=200) :: msg

      max_err = 0.0_wp
      do ip = 1, size(P_SET)
         do is = 1, N_S
            do it = 1, N_T
               sv_v = roquet_spv_value(T_SET(it), S_SET(is), P_SET(ip))
               call make_eos(eos, EOS_VARIANT_ROQUET_SPV)
               rho_b = eos_density_point(eos, T_SET(it), S_SET(is), P_SET(ip))
               call eos_specvol_derivs(eos, T_SET(it), S_SET(is), P_SET(ip), d1, d2)
               sv_p = 1.0_wp/rho_b
               max_err = max(max_err, abs(sv_v - sv_p)/sv_p)
               do iv = 1, size(VARS)
                  call make_eos(eos, VARS(iv))
                  rho_a = eos_density_point(eos, T_SET(it), S_SET(is), P_SET(ip))
                  call eos_specvol_derivs(eos, T_SET(it), S_SET(is), P_SET(ip), dt_a, ds_a)
                  call eos_density_specvol_derivs(eos, T_SET(it), S_SET(is), P_SET(ip), &
                                                  rho_b, dt_b, ds_b)
                  err = max(abs(rho_a - rho_b)/abs(rho_a), abs(dt_a - dt_b)/abs(dt_a), &
                            abs(ds_a - ds_b)/abs(ds_a))
                  max_err = max(max_err, err)
               end do
            end do
         end do
      end do
      write (msg, '(a,es10.3)') "value-only / fused EOS vs separate calls: max rel diff ", &
         max_err
      call global_logger%info(trim(msg))
      call check(error, max_err <= TOL_ULP, trim(msg))
   end subroutine test_entry_points

end module test_ocean_pgf_eos_fast
