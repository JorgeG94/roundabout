!! Unit tests for `seed_eady_ic` — the Eady-front IC overlay
!! invoked when `ic_config = "eady"`.
!!
!! Verifies:
!!   * dT/dz (vertical stratification) matches `cfg%ocean%ic%eady_dT_dz` to
!!     round-off (averaged across cells to wash out perturbation noise).
!!   * dT/dy (meridional gradient) matches `cfg%ocean%ic%eady_dT_dy`.
!!   * u_face_x shear matches the thermal-wind relation
!!     dU/dz = g · alpha_T / (rho_0 · f) · dT/dy.
!!   * Surface u(z = -dz/2) = +dU/dz · H/2 and bed u(z = -H + dz/2)
!!     = -dU/dz · H/2 — i.e. zero-mean shear about z_mid = -H/2.
!!   * Perturbation amplitude stays within ±eady_pert_amp/2.
module test_ocean_eady_ic
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY
   use rdb_config, only: config_t
   use rdb_grid, only: hgrid_t
   use rdb_ocean_state, only: ocean_state_t, ocean_state_seed_from_cfg
   implicit none
   private

   public :: collect_ocean_eady_ic_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NX_PHYS = 20
   integer, parameter :: NY_PHYS = 20
   integer, parameter :: NZ = 10

contains

   subroutine collect_ocean_eady_ic_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("eady_ic_dT_dz_matches_cfg", test_dT_dz), &
                  new_unittest("eady_ic_dT_dy_matches_cfg", test_dT_dy), &
                  new_unittest("eady_ic_thermal_wind_balance", test_thermal_wind), &
                  new_unittest("eady_ic_zero_mean_shear", test_zero_mean_shear), &
                  new_unittest("eady_ic_perturbation_bounded", test_pert_bounded) &
                  ]
   end subroutine collect_ocean_eady_ic_tests

   subroutine setup_eady(cfg, grid, state)
      type(config_t), intent(out) :: cfg
      type(hgrid_t), intent(out) :: grid
      type(ocean_state_t), intent(out) :: state

      cfg%sim_type = "ocean"
      cfg%ocean%topo%topo_config = "flat"
      cfg%ocean%ic%ic_config = "eady"
      cfg%ocean%topo%max_depth = 1000.0_wp
      cfg%coriolis_f = 1.0e-4_wp
      cfg%ocean%ic%alpha_T = 0.17_wp
      cfg%ocean%ic%rho_0 = 1035.0_wp
      cfg%ocean%ic%eady_dT_dy = -2.0e-5_wp
      cfg%ocean%ic%eady_dT_dz = 0.01_wp
      cfg%ocean%ic%eady_T_ref = 10.0_wp
      cfg%ocean%ic%eady_pert_amp = 1.0e-3_wp
      cfg%ocean%ic%eady_pert_seed = 42
      cfg%initial_salinity = 35.0_wp
      cfg%initial_temperature = 10.0_wp
      cfg%T_init_surface = 0.0_wp
      cfg%T_init_bottom = 0.0_wp

      call grid%init(NX_PHYS, NY_PHYS, NGHOST, 5000.0_wp, 5000.0_wp)
      state%multilayer%nz_ml = NZ
      call state%init(grid)
      ! Match the wiring done in state_init_from_config so the Eady IC
      ! sees the configured EOS params.
      state%eos%alpha_T = cfg%ocean%ic%alpha_T
      state%eos%rho0 = cfg%ocean%ic%rho_0

      call ocean_state_seed_from_cfg(state, grid, cfg)
   end subroutine setup_eady

   subroutine test_dT_dz(error)
      !! Verify vertical T gradient matches eady_dT_dz.  Average over a
      !! cell strip to wash out the small random perturbation.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      real(wp) :: T_k(NZ), dz_layer, computed
      integer :: k, idx_T, i_lo, i_hi, j_lo, j_hi

      call setup_eady(cfg, grid, state)
      idx_T = state%multilayer%idx_temperature

      i_lo = NGHOST + 1
      i_hi = NGHOST + NX_PHYS
      j_lo = NGHOST + 1
      j_hi = NGHOST + NY_PHYS

      do k = 1, NZ
         ! Concentration = hTr / h_layer, averaged over the physical interior.
         T_k(k) = sum(state%multilayer%tracers(idx_T)%hTr(i_lo:i_hi, j_lo:j_hi, k) &
                      /state%multilayer%h_layer(i_lo:i_hi, j_lo:j_hi, k)) &
                  /real(NX_PHYS*NY_PHYS, wp)
      end do
      dz_layer = cfg%ocean%topo%max_depth/real(NZ, wp)
      computed = (T_k(NZ) - T_k(1))/(real(NZ - 1, wp)*dz_layer)
      ! Tolerance set 6 orders below the value being verified — covers
      ! the residual O(1e-8) finite-sample mean of the noise without
      ! masking real algorithm errors.
      call check(error, abs(computed - cfg%ocean%ic%eady_dT_dz) < 1.0e-7_wp, &
                 "eady IC dT/dz must match cfg%ocean%ic%eady_dT_dz")
      call state%destroy()
   end subroutine test_dT_dz

   subroutine test_dT_dy(error)
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      real(wp) :: T_lo, T_hi, computed
      integer :: idx_T, j_lo, j_hi, k_mid

      call setup_eady(cfg, grid, state)
      idx_T = state%multilayer%idx_temperature

      ! Pick a layer in the middle so dT/dz contribution is roughly
      ! symmetric, then average over x to wash out the perturbation.
      k_mid = NZ/2
      j_lo = NGHOST + 2
      j_hi = NGHOST + NY_PHYS - 1
      T_lo = sum(state%multilayer%tracers(idx_T)%hTr(NGHOST + 1:NGHOST + NX_PHYS, j_lo, k_mid) &
                 /state%multilayer%h_layer(NGHOST + 1:NGHOST + NX_PHYS, j_lo, k_mid)) &
             /real(NX_PHYS, wp)
      T_hi = sum(state%multilayer%tracers(idx_T)%hTr(NGHOST + 1:NGHOST + NX_PHYS, j_hi, k_mid) &
                 /state%multilayer%h_layer(NGHOST + 1:NGHOST + NX_PHYS, j_hi, k_mid)) &
             /real(NX_PHYS, wp)
      computed = (T_hi - T_lo)/(real(j_hi - j_lo, wp)*grid%dy)
      call check(error, abs(computed - cfg%ocean%ic%eady_dT_dy) < 1.0e-8_wp, &
                 "eady IC dT/dy must match cfg%ocean%ic%eady_dT_dy")
      call state%destroy()
   end subroutine test_dT_dy

   subroutine test_thermal_wind(error)
      !! Verify the u(z) shear matches thermal-wind:
      !!   dU/dz = g · alpha_T / (rho_0 · f) · dT/dy
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      real(wp) :: u_top, u_bot, computed_dUdz, expected_dUdz, dz_layer

      call setup_eady(cfg, grid, state)

      ! u depends only on z, so any (i, j) gives the same answer.
      u_bot = state%multilayer%u_face_x_layer(NGHOST + 1, NGHOST + 1, 1)
      u_top = state%multilayer%u_face_x_layer(NGHOST + 1, NGHOST + 1, NZ)
      dz_layer = cfg%ocean%topo%max_depth/real(NZ, wp)
      computed_dUdz = (u_top - u_bot)/(real(NZ - 1, wp)*dz_layer)
      expected_dUdz = -GRAVITY*cfg%ocean%ic%alpha_T/(cfg%ocean%ic%rho_0*cfg%coriolis_f) &
                      *cfg%ocean%ic%eady_dT_dy
      call check(error, abs(computed_dUdz - expected_dUdz) < 1.0e-12_wp, &
                 "eady IC u_face shear must match thermal-wind dU/dz")
      call state%destroy()
   end subroutine test_thermal_wind

   subroutine test_zero_mean_shear(error)
      !! u(z=-dz/2) and u(z=-H+dz/2) should be equal in magnitude with
      !! opposite signs (jet centred at z_mid = -H/2).
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      real(wp) :: u_top, u_bot

      call setup_eady(cfg, grid, state)
      u_bot = state%multilayer%u_face_x_layer(NGHOST + 1, NGHOST + 1, 1)
      u_top = state%multilayer%u_face_x_layer(NGHOST + 1, NGHOST + 1, NZ)
      call check(error, abs(u_bot + u_top) < 1.0e-12_wp, &
                 "eady IC must produce zero-mean shear about z = -H/2")
      call state%destroy()
   end subroutine test_zero_mean_shear

   subroutine test_pert_bounded(error)
      !! After subtracting the analytical T(y, z), residual must stay
      !! within ±eady_pert_amp/2 inside the physical interior.  At the
      !! perturbation-free outer rows the residual is exactly zero.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      real(wp) :: T_analytic, residual, H, dz_layer, y_mid, y_phys, z_k
      real(wp) :: max_res
      integer :: i, j, k, idx_T, j_phys, j_lo, j_hi

      call setup_eady(cfg, grid, state)
      idx_T = state%multilayer%idx_temperature
      H = cfg%ocean%topo%max_depth
      dz_layer = H/real(NZ, wp)
      y_mid = 0.5_wp*real(NY_PHYS, wp)*grid%dy
      j_lo = NGHOST + 2
      j_hi = NGHOST + NY_PHYS - 1
      max_res = 0.0_wp
      do k = 1, NZ
         z_k = -H + (real(k, wp) - 0.5_wp)*dz_layer
         do j = j_lo, j_hi
            j_phys = j - NGHOST
            y_phys = (real(j_phys, wp) - 0.5_wp)*grid%dy
            T_analytic = cfg%ocean%ic%eady_T_ref + cfg%ocean%ic%eady_dT_dz*z_k &
                         + cfg%ocean%ic%eady_dT_dy*(y_phys - y_mid)
            do i = NGHOST + 1, NGHOST + NX_PHYS
               residual = state%multilayer%tracers(idx_T)%hTr(i, j, k) &
                          /state%multilayer%h_layer(i, j, k) - T_analytic
               max_res = max(max_res, abs(residual))
            end do
         end do
      end do
      call check(error, max_res <= 0.5_wp*cfg%ocean%ic%eady_pert_amp + 1.0e-12_wp, &
                 "eady IC perturbation must stay within ±amp/2")
      call state%destroy()
   end subroutine test_pert_bounded

end module test_ocean_eady_ic
