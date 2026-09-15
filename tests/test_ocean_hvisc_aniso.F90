!! Gap 5 — anisotropic viscosity + live velocity-scale viscosity for
!! the ocean horizontal-viscosity closure.
!!
!! Two sub-features under one capability:
!!
!! (a) LIVE velocity-scale viscosity (`kh_vel_scale_live`, MOM6
!!     `KH_VEL_SCALE` Kh = U·Δ made flow-aware via the live face
!!     speed).  An `A_vel = vel_scale·dx·|u|` contribution
!!     max-combined into the per-face harmonic viscosity.  Tested
!!     against the analytic `vel_scale·dx·|u|` on a uniform grid with
!!     closure=NONE and ah_bg=0.
!!
!! (b) ANISOTROPIC viscosity (`kh_aniso`, Smith & McWilliams 2003) on
!!     the stress-tensor path.  A direction tensor splits the
!!     viscosity into tension/shear coefficients + cross terms.
!!     Tested by direction-dependent dissipation: a pure i-tension
!!     flow is damped MORE with the grid-i anisotropy direction
!!     `(1,0)` (tension coefficient gains kh_aniso) than with the
!!     diagonal direction `(1,1)` (only the shear coefficient gains
!!     kh_aniso, but the shear strain is zero for this flow).
!!
!! Plus the mandatory both-off bit-identity gate.
module test_ocean_hvisc_aniso
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_lateral_mix, only: ocean_lateral_mix_t, &
                                    ocean_lateral_mix_compute, &
                                    LMIX_NONE
   use rdb_ocean_horizontal_viscosity, only: ocean_horizontal_viscosity_t, &
                                             ocean_horizontal_viscosity_compute_tendencies, &
                                             ocean_hvisc_set_aniso_direction, &
                                             aniso_mode_is_implemented
   use testdrive, only: error_type, check, new_unittest, unittest_type
   implicit none
   private

   public :: collect_ocean_hvisc_aniso_tests

contains

   subroutine collect_ocean_hvisc_aniso_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("vel_scale_live_matches_analytic", test_vel_scale_analytic), &
                  new_unittest("vel_scale_live_off_is_bit_identical", test_vel_scale_off_identical), &
                  new_unittest("aniso_direction_dependent_dissipation", test_aniso_direction), &
                  new_unittest("aniso_off_is_bit_identical", test_aniso_off_identical), &
                  new_unittest("aniso_mode_fail_loud", test_aniso_mode_fail_loud) &
                  ]
   end subroutine collect_ocean_hvisc_aniso_tests

   subroutine setup_state(grid, ms, nx, ny, nz, dx)
      type(hgrid_t), intent(inout) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: dx
      call grid%init(nx, ny, 1, dx, dx)
      ms%nz_ml = nz
      call ms%init(grid)
      ms%h_layer = 10.0_wp
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp
   end subroutine setup_state

   ! ----------------------------------------------------------------
   ! (a) Live velocity-scale viscosity.
   ! ----------------------------------------------------------------
   subroutine test_vel_scale_analytic(error)
      !! Uniform u-flow + closure=NONE + ah_bg=0.  Then every interior
      !! u-face viscosity equals the analytic `vel_scale·dx·|u|`
      !! (L_grid = dx on the uniform square grid).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_lateral_mix_t) :: lmix
      real(wp), parameter :: DX = 50000.0_wp
      real(wp), parameter :: U = 0.4_wp
      real(wp), parameter :: VEL_SCALE = 0.01_wp
      integer, parameter :: NX = 10, NY = 10, NZ = 1
      real(wp) :: expected, observed

      call setup_state(grid, ms, NX, NY, NZ, DX)
      ms%u_face_x_layer = U     ! uniform |u| everywhere
      ms%v_face_y_layer = 0.0_wp
      call lmix%init(grid, nz_ml=NZ)
      lmix%closure = LMIX_NONE
      lmix%ah_bg = 0.0_wp
      lmix%ah_max = 1.0e15_wp
      lmix%kh_vel_scale_live = VEL_SCALE

      call make_cartesian_metrics(metrics, grid)
      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(lmix)
      call lmix%enter_data()

      call ocean_lateral_mix_compute(grid, metrics, lmix, ms)
      !$acc update self(lmix%ah_face_x)

      call lmix%exit_data()
      !$acc exit data delete(lmix)
      call ms%exit_data()
      !$acc exit data delete(ms)
      call destroy_cartesian_metrics(metrics)

      expected = VEL_SCALE*DX*U
      ! Interior u-face (away from the i=1 wall column).
      observed = lmix%ah_face_x(NX/2, NY/2, 1)

      call check(error, abs(observed - expected) <= 1.0e-6_wp*expected, &
                 "live vel-scale ah_face_x must equal vel_scale*dx*|u|")

      call lmix%destroy()
      call ms%destroy()
   end subroutine test_vel_scale_analytic

   subroutine test_vel_scale_off_identical(error)
      !! kh_vel_scale_live = 0 with closure=NONE leaves ah_face_x at its
      !! initialised background — i.e. the live path is never entered.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_lateral_mix_t) :: lmix
      real(wp), parameter :: DX = 50000.0_wp, AH_BG = 123.0_wp
      integer, parameter :: NX = 8, NY = 8, NZ = 1
      real(wp) :: max_diff

      call setup_state(grid, ms, NX, NY, NZ, DX)
      ms%u_face_x_layer = 0.4_wp
      call lmix%init(grid, nz_ml=NZ)
      lmix%closure = LMIX_NONE
      lmix%ah_bg = AH_BG
      lmix%ah_face_x = AH_BG    ! mimic the configure-time background fill
      lmix%ah_face_y = AH_BG
      lmix%kh_vel_scale_live = 0.0_wp   ! OFF

      call make_cartesian_metrics(metrics, grid)
      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(lmix)
      call lmix%enter_data()

      call ocean_lateral_mix_compute(grid, metrics, lmix, ms)
      !$acc update self(lmix%ah_face_x, lmix%ah_face_y)

      call lmix%exit_data()
      !$acc exit data delete(lmix)
      call ms%exit_data()
      !$acc exit data delete(ms)
      call destroy_cartesian_metrics(metrics)

      max_diff = max(maxval(abs(lmix%ah_face_x - AH_BG)), &
                     maxval(abs(lmix%ah_face_y - AH_BG)))

      call check(error, max_diff == 0.0_wp, &
                 "kh_vel_scale_live=0 must leave the face viscosity untouched")

      call lmix%destroy()
      call ms%destroy()
   end subroutine test_vel_scale_off_identical

   ! ----------------------------------------------------------------
   ! (b) Anisotropic viscosity.
   ! ----------------------------------------------------------------
   subroutine run_aniso_case(grid, ms, metrics, nx, ny, nz, dx, nu_h, &
                             kh_aniso, n1, n2, du_visc)
      !! Drive the stress-tensor hvisc kernel on a fixed pure-i-tension
      !! flow with the given anisotropy magnitude + direction, and
      !! return the resulting east-face viscous tendency.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: dx, nu_h, kh_aniso, n1, n2
      type(hgrid_t), intent(inout) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_metrics_t), intent(inout) :: metrics
      real(wp), allocatable, intent(out) :: du_visc(:, :, :)
      type(ocean_horizontal_viscosity_t) :: hv
      real(wp), parameter :: DT = 1200.0_wp
      real(wp), parameter :: SHEAR = 0.02_wp
      integer :: i, j

      call setup_state(grid, ms, nx, ny, nz, dx)
      ! Pure i-tension flow: u_face_x(i,j) = SHEAR * i  (divergent in i),
      ! v = 0 ⇒ sh_xx = du/dx /= 0, sh_xy = 0.
      do j = 1, grid%ny_total
         do i = 1, grid%nx_total + 1
            ms%u_face_x_layer(i, j, 1) = SHEAR*real(i, wp)
         end do
      end do
      ms%v_face_y_layer = 0.0_wp

      call hv%init(grid, nz_ml=nz)
      hv%nu_h = nu_h
      hv%nu_4 = 0.0_wp
      hv%stress_tensor = .true.
      hv%kh_aniso = kh_aniso
      call ocean_hvisc_set_aniso_direction(hv, n1, n2)

      call make_cartesian_metrics(metrics, grid)
      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(hv)
      call hv%enter_data()

      call ocean_horizontal_viscosity_compute_tendencies(grid, metrics, hv, ms, dt=DT)
      !$acc update self(hv%du_visc%data)

      allocate (du_visc, source=hv%du_visc%data)

      call hv%exit_data()
      !$acc exit data delete(hv)
      call ms%exit_data()
      !$acc exit data delete(ms)
      call destroy_cartesian_metrics(metrics)

      call hv%destroy()
      call ms%destroy()
   end subroutine run_aniso_case

   subroutine test_aniso_direction(error)
      !! Direction-dependent dissipation.  For a pure i-tension flow
      !! (sh_xx /= 0, sh_xy = 0):
      !!   - grid-i `(1,0)`: tension coefficient gains kh_aniso ⇒ the
      !!     viscous tendency magnitude RISES above the isotropic
      !!     baseline.
      !!   - diagonal `(1,1)`: n1n2 = 1, n1²−n2² = 0 ⇒ tension gains 0,
      !!     only the shear coefficient gains kh_aniso, and the cross
      !!     term is zero — but the shear strain is zero for this flow,
      !!     so the tendency stays at the isotropic baseline.
      !! Hence dissipation along i > dissipation along the diagonal:
      !! the closure is genuinely anisotropic.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: DX = 50000.0_wp, NU_H = 1000.0_wp, KH_A = 5000.0_wp
      integer, parameter :: NX = 12, NY = 12, NZ = 1
      real(wp), allocatable :: du_iso(:, :, :), du_i(:, :, :), du_diag(:, :, :)
      real(wp) :: mag_iso, mag_i, mag_diag

      ! Isotropic baseline (kh_aniso = 0).
      call run_aniso_case(grid, ms, metrics, NX, NY, NZ, DX, NU_H, 0.0_wp, 1.0_wp, 0.0_wp, du_iso)
      ! Anisotropy along grid-i.
      call run_aniso_case(grid, ms, metrics, NX, NY, NZ, DX, NU_H, KH_A, 1.0_wp, 0.0_wp, du_i)
      ! Anisotropy along the diagonal.
      call run_aniso_case(grid, ms, metrics, NX, NY, NZ, DX, NU_H, KH_A, 1.0_wp, 1.0_wp, du_diag)

      mag_iso = maxval(abs(du_iso))
      mag_i = maxval(abs(du_i))
      mag_diag = maxval(abs(du_diag))

      ! (1) grid-i anisotropy strictly increases the tension dissipation.
      call check(error, mag_i > mag_iso*(1.0_wp + 1.0e-6_wp), &
                 "grid-i anisotropy must raise the i-tension viscous tendency above isotropic")
      if (allocated(error)) return

      ! (2) diagonal anisotropy leaves this flow's tendency at the
      !     isotropic baseline (tension add = 0, shear strain = 0).
      call check(error, abs(mag_diag - mag_iso) <= 1.0e-8_wp*mag_iso + 1.0e-12_wp, &
                 "diagonal anisotropy must NOT change the pure-i-tension tendency")
      if (allocated(error)) return

      ! (3) Direction matters: i-tendency differs from diagonal-tendency.
      call check(error, mag_i > mag_diag*(1.0_wp + 1.0e-6_wp), &
                 "dissipation must be direction-dependent (i > diagonal for i-tension flow)")
   end subroutine test_aniso_direction

   subroutine test_aniso_off_identical(error)
      !! kh_aniso = 0 reproduces the isotropic stress-tensor tendency
      !! bit-for-bit regardless of the (unused) direction vector.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: DX = 50000.0_wp, NU_H = 1000.0_wp
      integer, parameter :: NX = 12, NY = 12, NZ = 1
      real(wp), allocatable :: du_a(:, :, :), du_b(:, :, :)
      real(wp) :: max_diff

      ! Two runs, kh_aniso=0 but different direction vectors — must agree.
      call run_aniso_case(grid, ms, metrics, NX, NY, NZ, DX, NU_H, 0.0_wp, 1.0_wp, 0.0_wp, du_a)
      call run_aniso_case(grid, ms, metrics, NX, NY, NZ, DX, NU_H, 0.0_wp, 1.0_wp, 1.0_wp, du_b)

      max_diff = maxval(abs(du_a - du_b))
      call check(error, max_diff == 0.0_wp, &
                 "kh_aniso=0 must be bit-identical regardless of direction")
   end subroutine test_aniso_off_identical

   subroutine test_aniso_mode_fail_loud(error)
      !! Only anisotropy mode 0 (grid-relative `aniso_dir`) has an
      !! implemented direction tensor.  This asserts the
      !! `aniso_mode_is_implemented` predicate that `validate_config`
      !! aborts on — exercising the fail-loud guard's code path (the
      !! `error stop` itself can't run in-process).
      type(error_type), allocatable, intent(out) :: error

      call check(error, aniso_mode_is_implemented(0), &
                 "aniso_mode=0 (grid-relative) must be implemented")
      if (allocated(error)) return
      call check(error,.not. aniso_mode_is_implemented(1), &
                 "aniso_mode=1 must report unimplemented (drives fail-loud abort)")
      if (allocated(error)) return
      call check(error,.not. aniso_mode_is_implemented(-1), &
                 "negative aniso_mode must report unimplemented")
   end subroutine test_aniso_mode_fail_loud

end module test_ocean_hvisc_aniso
