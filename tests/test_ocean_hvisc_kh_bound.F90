!! Per-face harmonic-viscosity CFL clamp on the velocity-Laplacian
!! paths (`&ocean_hvisc_nml bound_kh`, MOM6 `BOUND_KH` analogue).
!!
!! When `bound_kh = .true.` the effective Laplacian viscosity at each
!! face is clamped to
!!     `ν_max = bound_coef · 0.125 / (dt · (1/dx² + 1/dy²))`
!! — one quarter of the forward-Euler stability limit.  Beyond simple
!! FE stability this keeps the FROZEN depth-mean viscous forcing that
!! the split-RK2 barotropic mode receives (via `F_bt`) out of the
!! phase-reversed anti-damping regime for grid-scale gravity modes
!! (the 600² Lagrangian double-gyre rim blow-up).
!!
!! Covers (scalar `nu_h` velocity-Laplacian path, uniform grid):
!!   1. `nu_h` far above the bound + `bound_kh` → tendency equals a run
!!      with `nu_h` set exactly to the analytic bound.
!!   2. `nu_h` below the bound + `bound_kh` → byte-identical to the
!!      unbounded run (the default-off bit-identity gate).
module test_ocean_hvisc_kh_bound
   use rdb_constants, only: wp, PI
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_horizontal_viscosity, only: ocean_horizontal_viscosity_t, &
                                             ocean_horizontal_viscosity_compute_tendencies
   use testdrive, only: error_type, check, new_unittest, unittest_type
   implicit none
   private

   public :: collect_ocean_hvisc_kh_bound_tests

   real(wp), parameter :: DX = 1000.0_wp       !! fine grid spacing (m)
   real(wp), parameter :: DT = 600.0_wp        !! outer step (s)
   real(wp), parameter :: BOUND_COEF = 0.8_wp  !! default CFL safety margin
   integer, parameter :: NX = 12, NY = 12, NZ = 2

contains

   subroutine collect_ocean_hvisc_kh_bound_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("kh_bound_large_nu_clamped", test_large_nu_clamped), &
                  new_unittest("kh_bound_small_nu_untouched", test_small_nu_untouched) &
                  ]
   end subroutine collect_ocean_hvisc_kh_bound_tests

   pure function kh_bound() result(kh_max)
      !! Analytic per-face harmonic ceiling for the uniform grid
      !! (`idx = idy = 1/DX`), matching `hvisc_kh_cfl_bound`.
      real(wp) :: kh_max
      real(wp) :: idx, k2
      idx = 1.0_wp/DX
      k2 = idx*idx + idx*idx
      kh_max = BOUND_COEF*0.125_wp/(DT*k2)
   end function kh_bound

   subroutine run_scalar_laplacian(nu_h, do_bound, du_out)
      !! Run the scalar-`nu_h` velocity-Laplacian compute on a fresh
      !! state and return the east-face tendency.
      real(wp), intent(in) :: nu_h
      logical, intent(in) :: do_bound
      real(wp), allocatable, intent(out) :: du_out(:, :, :)
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_horizontal_viscosity_t) :: hv
      integer :: i, j, k

      call grid%init(NX, NY, 1, DX, DX)
      ms%nz_ml = NZ
      call ms%init(grid)
      ms%h_layer = 10.0_wp
      ms%v_face_y_layer = 0.0_wp
      ! Short-wavelength sinusoid in i: non-zero Laplacian.
      do k = 1, NZ
         do j = 1, grid%ny_total
            do i = 1, grid%nx_total + 1
               ms%u_face_x_layer(i, j, k) = sin(2.0_wp*PI*real(i - 1, wp)/real(NX, wp))
            end do
         end do
      end do
      call hv%init(grid, nz_ml=NZ)
      hv%nu_h = nu_h
      hv%nu_4 = 0.0_wp
      hv%bound_coef = BOUND_COEF
      hv%bound_kh = do_bound
      call make_cartesian_metrics(metrics, grid)
      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(hv)
      call hv%enter_data()

      call ocean_horizontal_viscosity_compute_tendencies(grid, metrics, hv, ms, dt=DT)
      !$acc update self(hv%du_visc%data)
      du_out = hv%du_visc%data

      call hv%exit_data()
      !$acc exit data delete(hv)
      call ms%exit_data()
      !$acc exit data delete(ms)
      call destroy_cartesian_metrics(metrics)
      call hv%destroy()
      call ms%destroy()
   end subroutine run_scalar_laplacian

   subroutine test_large_nu_clamped(error)
      !! `nu_h` 100× above the bound with `bound_kh` must produce the
      !! SAME tendency as `nu_h` set exactly to the bound (uniform grid
      !! ⇒ every interior face clamps to the same value).
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: du_big(:, :, :), du_lim(:, :, :)
      real(wp) :: bound, max_diff, ref
      bound = kh_bound()
      call run_scalar_laplacian(100.0_wp*bound, .true., du_big)
      call run_scalar_laplacian(bound, .false., du_lim)
      max_diff = maxval(abs(du_big - du_lim))
      ref = maxval(abs(du_lim))
      call check(error, ref > 0.0_wp, "Laplacian tendency must be non-trivial")
      if (allocated(error)) return
      call check(error, max_diff <= 1.0e-9_wp*ref, &
                 "nu_h above the bound must clamp to the analytic kh ceiling")
   end subroutine test_large_nu_clamped

   subroutine test_small_nu_untouched(error)
      !! `nu_h` below the bound: `bound_kh = .true.` must be
      !! byte-identical to `bound_kh = .false.` (bit-identity gate).
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: du_on(:, :, :), du_off(:, :, :)
      real(wp) :: bound
      bound = kh_bound()
      call run_scalar_laplacian(0.5_wp*bound, .true., du_on)
      call run_scalar_laplacian(0.5_wp*bound, .false., du_off)
      call check(error, maxval(abs(du_on)) > 0.0_wp, &
                 "Laplacian tendency must be non-trivial")
      if (allocated(error)) return
      call check(error, maxval(abs(du_on - du_off)) == 0.0_wp, &
                 "nu_h below the bound must be untouched by bound_kh")
   end subroutine test_small_nu_untouched

end module test_ocean_hvisc_kh_bound
