!! Unit tests for the linear ocean EOS kernel
!! (rdb_eos Phase 5c).  Converts the registered S, T tracers
!! into a cell-centred density field via
!!
!!   rho = rho_0 + beta_S * (S - S_ref) - alpha_T * (T - T_ref)
!!
!! where S = hS / h, T = hT / h.  Outer-shim + flat-impl pattern
!! avoids the NVHPC stdpar deep-deref on the tracer registry.
!!
!! Cases:
!!   * Reference-only state — set S = S_ref, T = T_ref everywhere;
!!     rho must equal rho_0 to round-off.  Sanity check that the
!!     EOS produces no spurious anomaly at the reference point.
!!   * Pure salinity anomaly — set S = S_ref + dS, T = T_ref;
!!     rho - rho_0 must equal beta_S * dS to round-off.  Confirms
!!     the salinity branch.
!!   * Pure temperature anomaly — set T = T_ref + dT, S = S_ref;
!!     rho - rho_0 must equal -alpha_T * dT.  Confirms the
!!     temperature branch (with the correct sign).
!!   * Combined anomaly with layer-dependent values — different
!!     dS_k + dT_k per layer; rho_layer must match the linear
!!     formula in every layer.  Catches cross-layer index bugs.
module test_ocean_eos
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_eos, only: eos_t
   use rdb_ocean_eos_compute, only: ocean_eos_compute
   implicit none
   private

   public :: collect_ocean_eos_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 3

contains

   subroutine collect_ocean_eos_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("eos_at_reference", test_at_reference), &
                  new_unittest("eos_pure_salinity_anomaly", test_pure_salt), &
                  new_unittest("eos_pure_temperature_anomaly", test_pure_temp), &
                  new_unittest("eos_combined_per_layer", test_combined_per_layer) &
                  ]
   end subroutine collect_ocean_eos_tests

   subroutine make_grid(grid)
      type(hgrid_t), intent(out) :: grid
      call grid%init(8, 6, NGHOST, 1.0_wp, 1.0_wp)
   end subroutine make_grid

   subroutine run_eos(ms, eos)
      type(multilayer_state_t), intent(inout) :: ms
      type(eos_t), intent(inout) :: eos
      !$acc enter data copyin(ms)
      call ms%enter_data()
      call ocean_eos_compute(eos, ms)
      call ms%exit_data()
      !$acc exit data delete(ms)
   end subroutine run_eos

   ! -----------------------------------------------------------------
   ! Cases
   ! -----------------------------------------------------------------

   subroutine test_at_reference(error)
      !! S = S_ref, T = T_ref -> rho = rho_0 to round-off.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      real(wp), parameter :: H0 = 10.0_wp
      real(wp) :: max_diff

      call make_grid(grid)
      ms%nz_ml = NZ
      call ms%init(grid)
      call eos%init(grid)
      ms%h_layer = H0
      ms%tracers(ms%idx_salinity)%hTr = eos%S_ref*H0
      ms%tracers(ms%idx_temperature)%hTr = eos%T_ref*H0

      call run_eos(ms, eos)

      max_diff = maxval(abs(ms%rho_layer - eos%rho0))
      call check(error, max_diff < 1.0e-12_wp, &
                 "EOS at reference state: rho deviates from rho_0")

      call eos%destroy()
      call ms%destroy()
   end subroutine test_at_reference

   subroutine test_pure_salt(error)
      !! S = S_ref + dS, T = T_ref -> rho - rho_0 = beta_S * dS.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      real(wp), parameter :: H0 = 5.0_wp
      real(wp), parameter :: DS = 2.5_wp
      real(wp) :: max_diff, expected_anom

      call make_grid(grid)
      ms%nz_ml = NZ
      call ms%init(grid)
      call eos%init(grid)
      ms%h_layer = H0
      ms%tracers(ms%idx_salinity)%hTr = (eos%S_ref + DS)*H0
      ms%tracers(ms%idx_temperature)%hTr = eos%T_ref*H0

      call run_eos(ms, eos)

      expected_anom = eos%beta_S*DS
      max_diff = maxval(abs(ms%rho_layer - (eos%rho0 + expected_anom)))
      call check(error, max_diff < 1.0e-12_wp, &
                 "pure salinity anomaly: rho-rho_0 != beta_S*dS")

      call eos%destroy()
      call ms%destroy()
   end subroutine test_pure_salt

   subroutine test_pure_temp(error)
      !! S = S_ref, T = T_ref + dT -> rho - rho_0 = -alpha_T * dT.
      !! Note the sign — temperature anomaly *reduces* density.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      real(wp), parameter :: H0 = 5.0_wp
      real(wp), parameter :: DT = 4.0_wp
      real(wp) :: max_diff, expected_anom

      call make_grid(grid)
      ms%nz_ml = NZ
      call ms%init(grid)
      call eos%init(grid)
      ms%h_layer = H0
      ms%tracers(ms%idx_salinity)%hTr = eos%S_ref*H0
      ms%tracers(ms%idx_temperature)%hTr = (eos%T_ref + DT)*H0

      call run_eos(ms, eos)

      expected_anom = -eos%alpha_T*DT
      max_diff = maxval(abs(ms%rho_layer - (eos%rho0 + expected_anom)))
      call check(error, max_diff < 1.0e-12_wp, &
                 "pure temperature anomaly: rho-rho_0 != -alpha_T*dT")

      call eos%destroy()
      call ms%destroy()
   end subroutine test_pure_temp

   subroutine test_combined_per_layer(error)
      !! Per-layer dS_k and dT_k anomalies + non-uniform h per layer.
      !! Asserts the full linear formula matches in every cell.
      !! Different (dS_k, dT_k) per layer specifically guards against
      !! cross-layer index bugs in the kernel.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      real(wp) :: h_k(NZ), dS_k(NZ), dT_k(NZ), expected_rho(NZ)
      real(wp) :: max_diff
      integer :: k
      checks: block

         call make_grid(grid)
         ms%nz_ml = NZ
         call ms%init(grid)
         call eos%init(grid)

         h_k = [3.0_wp, 5.0_wp, 7.0_wp]
         dS_k = [+1.0_wp, -0.5_wp, +2.0_wp]
         dT_k = [-2.0_wp, +1.5_wp, -3.0_wp]
         do k = 1, NZ
            ms%h_layer(:, :, k) = h_k(k)
            ms%tracers(ms%idx_salinity)%hTr(:, :, k) = (eos%S_ref + dS_k(k))*h_k(k)
            ms%tracers(ms%idx_temperature)%hTr(:, :, k) = (eos%T_ref + dT_k(k))*h_k(k)
            expected_rho(k) = eos%rho0 + eos%beta_S*dS_k(k) - eos%alpha_T*dT_k(k)
         end do

         call run_eos(ms, eos)

         do k = 1, NZ
            max_diff = maxval(abs(ms%rho_layer(:, :, k) - expected_rho(k)))
            call check(error, max_diff < 1.0e-12_wp, &
                       "per-layer EOS: rho deviates from analytic")
            if (allocated(error)) exit checks
         end do

      end block checks
      call eos%destroy()
      call ms%destroy()
   end subroutine test_combined_per_layer

end module test_ocean_eos
