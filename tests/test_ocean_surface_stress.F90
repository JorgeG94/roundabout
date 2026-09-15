!! Unit tests for the ocean surface-stress kernel
!! (rdb_ocean_surface_stress).  Surface wind stress acts only on
!! the surface-most layer (k = nz under the ROMS-style convention).
!!
!! Cases:
!!   * Linear momentum input — uniform tau_x, uniform h_top.
!!     After N steps the surface-layer velocity at an interior
!!     face must equal `N * dt * tau_x / (rho_0 * h_top)` exactly.
!!   * Surface-only support — non-trivial multilayer u, v.  After
!!     one stress step every layer k < nz must remain at its IC.
!!     Mirror of the bottom-drag bed-only test.
!!   * Zero wind no-op — tau = 0 with non-trivial IC must leave
!!     u, v untouched.
module test_ocean_surface_stress
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t, &
                                       ocean_surface_stress_compute_tendencies, &
                                       ocean_surface_stress_apply_tendencies
   implicit none
   private

   public :: collect_ocean_surface_stress_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 3

contains

   subroutine collect_ocean_surface_stress_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("linear_momentum_input", test_linear_input), &
                  new_unittest("surface_only_support", test_surface_only), &
                  new_unittest("zero_wind_no_op", test_zero_wind), &
                  new_unittest("scalar_setter_recovers_uniform", &
                               test_scalar_setter_uniform), &
                  new_unittest("spatially_varying_wind", &
                               test_spatially_varying_wind), &
                  new_unittest("direct_stress_thin_sbl_matches_bed_only", &
                               test_direct_stress_thin), &
                  new_unittest("direct_stress_spans_multi_layer", &
                               test_direct_stress_multi) &
                  ]
   end subroutine collect_ocean_surface_stress_tests

   subroutine make_grid(grid, nx_phys, ny_phys, dx, dy)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      real(wp), intent(in) :: dx, dy
      call grid%init(nx_phys, ny_phys, NGHOST, dx, dy)
   end subroutine make_grid

   subroutine map_in(ms, ss)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_surface_stress_t), intent(inout) :: ss
      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(ss)
      call ss%enter_data()
   end subroutine map_in

   subroutine map_out(ms, ss)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_surface_stress_t), intent(inout) :: ss
      call ss%exit_data()
      !$acc exit data delete(ss)
      call ms%exit_data()
      !$acc exit data delete(ms)
   end subroutine map_out

   subroutine test_linear_input(error)
      !! tau_x = TAU0, uniform h_top = H0, zero initial u, v.
      !! After N steps the surface-layer u at an interior face
      !! must equal N * dt * TAU0 / (rho0 * H0) — pure first-
      !! integral of a constant acceleration.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_surface_stress_t) :: ss
      real(wp), parameter :: TAU0 = 0.1_wp
      real(wp), parameter :: H0 = 10.0_wp
      real(wp), parameter :: RHO0 = 1035.0_wp
      real(wp), parameter :: DT = 0.5_wp
      integer, parameter :: N_STEPS = 6
      real(wp) :: u_expected, u_obs, max_lower
      integer :: step, nx, ny
      checks: block

         call make_grid(grid, 10, 8, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ss%init(grid, nz_ml=NZ)
         ss%rho0 = RHO0
         call ss%set_wind_stress_const(TAU0, 0.0_wp)
         nx = grid%nx_total
         ny = grid%ny_total

         ms%h_layer = H0
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp

         call map_in(ms, ss)
         do step = 1, N_STEPS
            call ocean_surface_stress_compute_tendencies(grid, ss, ms)
            call ocean_surface_stress_apply_tendencies(ss, ms, DT)
         end do
         call map_out(ms, ss)

         u_expected = real(N_STEPS, wp)*DT*TAU0/(RHO0*H0)
         u_obs = ms%u_face_x_layer(nx/2, ny/2, NZ)
         max_lower = maxval(abs(ms%u_face_x_layer(:, :, 1:NZ - 1))) + &
                     maxval(abs(ms%v_face_y_layer(:, :, 1:NZ - 1)))

         call check(error, abs(u_obs - u_expected) < 1.0e-14_wp, &
                    "linear momentum input: surface u off analytic")
         if (allocated(error)) exit checks
         call check(error, max_lower < 1.0e-12_wp, &
                    "linear input: stress leaked into layers k<nz")

      end block checks
      call ss%destroy(); call ms%destroy()
   end subroutine test_linear_input

   subroutine test_surface_only(error)
      !! Non-trivial multilayer u, v.  One stress step: every layer
      !! k < nz must be byte-identical to the IC.  Catches whole-
      !! column writes (mirror of the bottom-drag bed-only test).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_surface_stress_t) :: ss
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp), parameter :: DT = 0.05_wp
      real(wp), allocatable :: u_ic(:, :, :), v_ic(:, :, :)
      real(wp) :: max_du_lower, max_dv_lower
      integer :: i, j, k, nx, ny
      checks: block

         call make_grid(grid, 12, 10, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ss%init(grid, nz_ml=NZ)
         call ss%set_wind_stress_const(0.2_wp, -0.1_wp)
         nx = grid%nx_total
         ny = grid%ny_total

         ms%h_layer = 10.0_wp
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx + 1
                  ms%u_face_x_layer(i, j, k) = sin(2.0_wp*PI*real(i + k, wp)/real(nx, wp))
               end do
            end do
            do j = 1, ny + 1
               do i = 1, nx
                  ms%v_face_y_layer(i, j, k) = cos(2.0_wp*PI*real(j + k, wp)/real(ny, wp))
               end do
            end do
         end do
         allocate (u_ic, source=ms%u_face_x_layer)
         allocate (v_ic, source=ms%v_face_y_layer)

         call map_in(ms, ss)
         call ocean_surface_stress_compute_tendencies(grid, ss, ms)
         call ocean_surface_stress_apply_tendencies(ss, ms, DT)
         call map_out(ms, ss)

         max_du_lower = maxval(abs(ms%u_face_x_layer(:, :, 1:NZ - 1) - u_ic(:, :, 1:NZ - 1)))
         max_dv_lower = maxval(abs(ms%v_face_y_layer(:, :, 1:NZ - 1) - v_ic(:, :, 1:NZ - 1)))

         call check(error, max_du_lower < 1.0e-12_wp, &
                    "surface-only: u in layers k<nz was touched")
         if (allocated(error)) exit checks
         call check(error, max_dv_lower < 1.0e-12_wp, &
                    "surface-only: v in layers k<nz was touched")

      end block checks
      deallocate (u_ic, v_ic)
      call ss%destroy(); call ms%destroy()
   end subroutine test_surface_only

   subroutine test_zero_wind(error)
      !! tau_x = tau_y = 0 → tendency stays zero; apply no-op.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_surface_stress_t) :: ss
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp), parameter :: DT = 0.05_wp
      real(wp), allocatable :: u_ic(:, :, :), v_ic(:, :, :)
      integer :: i, j, k, nx, ny
      checks: block

         call make_grid(grid, 10, 8, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ss%init(grid, nz_ml=NZ)
         call ss%set_wind_stress_const(0.0_wp, 0.0_wp)
         nx = grid%nx_total
         ny = grid%ny_total

         ms%h_layer = 10.0_wp
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx + 1
                  ms%u_face_x_layer(i, j, k) = 0.2_wp*sin(PI*real(i, wp)/real(nx, wp))
               end do
            end do
            do j = 1, ny + 1
               do i = 1, nx
                  ms%v_face_y_layer(i, j, k) = 0.1_wp*cos(PI*real(j, wp)/real(ny, wp))
               end do
            end do
         end do
         allocate (u_ic, source=ms%u_face_x_layer)
         allocate (v_ic, source=ms%v_face_y_layer)

         call map_in(ms, ss)
         call ocean_surface_stress_compute_tendencies(grid, ss, ms)
         call ocean_surface_stress_apply_tendencies(ss, ms, DT)
         call map_out(ms, ss)

         call check(error, maxval(abs(ms%u_face_x_layer - u_ic)) < 1.0e-12_wp, &
                    "zero wind: u changed")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(ms%v_face_y_layer - v_ic)) < 1.0e-12_wp, &
                    "zero wind: v changed")

      end block checks
      deallocate (u_ic, v_ic)
      call ss%destroy(); call ms%destroy()
   end subroutine test_zero_wind

   subroutine test_scalar_setter_uniform(error)
      !! `set_wind_stress_const(tau_x, tau_y)` must fill the 2D
      !! tau_x / tau_y arrays uniformly, producing the same
      !! tendency at every interior surface face.  Verifies the
      !! backward-compat scalar pathway still gives a spatially-
      !! homogeneous response.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_surface_stress_t) :: ss
      real(wp), parameter :: TAU = 0.15_wp
      real(wp), parameter :: H0 = 8.0_wp
      real(wp), parameter :: RHO0 = 1035.0_wp
      real(wp) :: du_expected, max_dev_u, max_dev_v
      integer :: nx, ny
      checks: block

         call make_grid(grid, 10, 8, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ss%init(grid, nz_ml=NZ)
         ss%rho0 = RHO0
         call ss%set_wind_stress_const(TAU, -TAU)
         ms%h_layer = H0
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         nx = grid%nx_total
         ny = grid%ny_total

         call map_in(ms, ss)
         call ocean_surface_stress_compute_tendencies(grid, ss, ms)
         !$acc update self(ss%du_stress%data, ss%dv_stress%data)
         call map_out(ms, ss)

         du_expected = TAU/(RHO0*H0)
         ! Interior u-faces (i = 2..nx) at k = nz must all equal du_expected.
         max_dev_u = maxval(abs(ss%du_stress%data(2:nx, :, NZ) - du_expected))
         ! Interior v-faces (j = 2..ny) at k = nz must all equal -du_expected
         ! (tau_y = -TAU, same h_top, same rho0).
         max_dev_v = maxval(abs(ss%dv_stress%data(:, 2:ny, NZ) + du_expected))

         call check(error, max_dev_u < 1.0e-14_wp, &
                    "scalar setter did not yield uniform du_stress")
         if (allocated(error)) exit checks
         call check(error, max_dev_v < 1.0e-14_wp, &
                    "scalar setter did not yield uniform dv_stress")

      end block checks
      call ss%destroy(); call ms%destroy()
   end subroutine test_scalar_setter_uniform

   subroutine test_spatially_varying_wind(error)
      !! Spatially varying `tau_x(i, j)` — sinusoidal in j — must
      !! produce a `du_stress` field that mirrors the input pattern
      !! (after the `1 / (rho_0 * h_top)` normalization).  Probes
      !! the new 2D-aware pathway end-to-end: array assignment in
      !! host memory → enter_data → kernel reads via 2D index.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_surface_stress_t) :: ss
      real(wp), parameter :: TAU0 = 0.2_wp
      real(wp), parameter :: H0 = 10.0_wp
      real(wp), parameter :: RHO0 = 1035.0_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp) :: tau_at_j, expected, max_err
      integer :: i, j, nx, ny
      checks: block

         call make_grid(grid, 10, 8, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call ss%init(grid, nz_ml=NZ)
         ss%rho0 = RHO0
         nx = grid%nx_total
         ny = grid%ny_total

         ! Sinusoidal tau_x(i, j) = TAU0 * sin(π·j/ny), independent of i.
         do j = 1, ny
            tau_at_j = TAU0*sin(PI*real(j, wp)/real(ny, wp))
            do i = 1, nx + 1
               ss%tau_x(i, j) = tau_at_j
            end do
         end do
         ss%tau_y = 0.0_wp

         ms%h_layer = H0
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp

         call map_in(ms, ss)
         call ocean_surface_stress_compute_tendencies(grid, ss, ms)
         !$acc update self(ss%du_stress%data, ss%dv_stress%data)
         call map_out(ms, ss)

         ! Each interior u-face (i = 2..nx, j) at k = nz must equal
         ! tau_x(i, j) / (rho_0 * H0).
         max_err = 0.0_wp
         do j = 1, ny
            do i = 2, nx
               expected = ss%tau_x(i, j)/(RHO0*H0)
               max_err = max(max_err, abs(ss%du_stress%data(i, j, NZ) - expected))
            end do
         end do

         call check(error, max_err < 1.0e-14_wp, &
                    "spatial wind: du_stress did not mirror tau_x(i,j)")
         if (allocated(error)) exit checks

         ! tau_y = 0 → dv_stress(:, 2:ny, nz) = 0.
         call check(error, &
                    maxval(abs(ss%dv_stress%data(:, 2:ny, NZ))) < 1.0e-14_wp, &
                    "spatial wind: dv_stress non-zero with tau_y = 0")

      end block checks
      call ss%destroy(); call ms%destroy()
   end subroutine test_spatially_varying_wind

   subroutine test_direct_stress_thin(error)
      !! DIRECT_STRESS with `hmix_stress` smaller than the surface
      !! layer thickness — the SBL is fully contained in layer nz.
      !! The distributed kernel must collapse to the bed-only result
      !! at the surface layer, with zero acceleration in layers k<nz.
      !!
      !! For TAU0=0.1, RHO0=1035, H0=10, hmix_stress=2:
      !!   bed-only:        a_nz = TAU0/(RHO0*H0) = 9.66e-6 m/s²
      !!   distributed:     a_nz = TAU0/(RHO0*hmix_stress) * (hmix_stress/H0)
      !!                          = TAU0/(RHO0*H0)
      !!   → bit-identical.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_surface_stress_t) :: ss
      real(wp), parameter :: TAU0 = 0.1_wp, H0 = 10.0_wp, RHO0 = 1035.0_wp
      real(wp), parameter :: HMIX = 2.0_wp, DT = 0.5_wp
      real(wp) :: u_expected, u_obs, max_lower
      integer :: nx, ny

      call make_grid(grid, 8, 6, 1.0_wp, 1.0_wp)
      ms%nz_ml = NZ
      call ms%init(grid)
      call ss%init(grid, nz_ml=NZ)
      ss%rho0 = RHO0
      ss%direct_stress = .true.
      ss%hmix_stress = HMIX
      call ss%set_wind_stress_const(TAU0, 0.0_wp)
      nx = grid%nx_total
      ny = grid%ny_total

      ms%h_layer = H0
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp

      call map_in(ms, ss)
      call ocean_surface_stress_compute_tendencies(grid, ss, ms)
      call ocean_surface_stress_apply_tendencies(ss, ms, DT)
      call map_out(ms, ss)

      u_expected = DT*TAU0/(RHO0*H0)
      u_obs = ms%u_face_x_layer(nx/2, ny/2, NZ)
      max_lower = maxval(abs(ms%u_face_x_layer(:, :, 1:NZ - 1))) + &
                  maxval(abs(ms%v_face_y_layer(:, :, 1:NZ - 1)))

      call check(error, abs(u_obs - u_expected) < 1.0e-14_wp, &
                 "DIRECT_STRESS thin SBL: surface u not bit-equiv to bed-only")
      if (.not. allocated(error)) call check(error, max_lower < 1.0e-12_wp, &
                                             "DIRECT_STRESS thin SBL: leaked into layers k<nz")
      call ss%destroy(); call ms%destroy()
   end subroutine test_direct_stress_thin

   subroutine test_direct_stress_multi(error)
      !! DIRECT_STRESS with `hmix_stress` spanning two layers.
      !! H_k = 1 m per layer (NZ = 3), hmix_stress = 1.5 m.  The
      !! surface layer (k = NZ) takes the top 1 m; layer k = NZ-1
      !! takes the remaining 0.5 m.  Total impulse summed over k
      !! must equal tau · dt / rho0 (per unit area), partitioned
      !! according to the fractional SBL coverage.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_surface_stress_t) :: ss
      real(wp), parameter :: TAU0 = 0.5_wp, HK = 1.0_wp, RHO0 = 1000.0_wp
      real(wp), parameter :: HMIX = 1.5_wp, DT = 1.0_wp
      real(wp) :: u_nz, u_nzm1, sum_h_u, expected_impulse
      integer :: nx, ny, k_top

      call make_grid(grid, 8, 6, 1.0_wp, 1.0_wp)
      ms%nz_ml = NZ
      call ms%init(grid)
      call ss%init(grid, nz_ml=NZ)
      ss%rho0 = RHO0
      ss%direct_stress = .true.
      ss%hmix_stress = HMIX
      call ss%set_wind_stress_const(TAU0, 0.0_wp)
      nx = grid%nx_total
      ny = grid%ny_total
      k_top = NZ

      ms%h_layer = HK
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp

      call map_in(ms, ss)
      call ocean_surface_stress_compute_tendencies(grid, ss, ms)
      call ocean_surface_stress_apply_tendencies(ss, ms, DT)
      call map_out(ms, ss)

      u_nz = ms%u_face_x_layer(nx/2, ny/2, k_top)
      u_nzm1 = ms%u_face_x_layer(nx/2, ny/2, k_top - 1)
      ! Per-layer formula: u_k = tau · dt / (rho · hmix) · (h_in_sbl_k / h_k)
      !   k=NZ:   h_in_sbl = 1.0, so u = 0.5·1.0/(1000·1.5)·(1.0/1.0) = 3.33e-4
      !   k=NZ-1: h_in_sbl = 0.5, so u = 0.5·1.0/(1000·1.5)·(0.5/1.0) = 1.67e-4
      call check(error, abs(u_nz - 0.5_wp*1.0_wp/(RHO0*HMIX)) < 1.0e-14_wp, &
                 "DIRECT_STRESS multi-SBL: surface layer u wrong")
      if (.not. allocated(error)) call check(error, &
                                             abs(u_nzm1 - 0.5_wp*1.0_wp/(RHO0*HMIX)*0.5_wp) < 1.0e-14_wp, &
                                             "DIRECT_STRESS multi-SBL: sub-surface layer u wrong")
      ! Sum of (h_k · u_k) over SBL = tau · dt / rho0 (total impulse per area).
      sum_h_u = HK*u_nz + HK*u_nzm1
      expected_impulse = TAU0*DT/RHO0
      if (.not. allocated(error)) call check(error, abs(sum_h_u - expected_impulse) < 1.0e-14_wp, &
                                             "DIRECT_STRESS multi-SBL: total impulse not conserved")
      call ss%destroy(); call ms%destroy()
   end subroutine test_direct_stress_multi

end module test_ocean_surface_stress
