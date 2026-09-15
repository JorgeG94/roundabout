!! Per-cell biharmonic CFL clamp on the ocean horizontal-viscosity
!! kernel (Gap 4 — doc↔code reconcile).
!!
!! The scalar `nu_4` and the flow-aware Smag_AH `nu4_face_*` biharmonic
!! coefficients are clamped per face to the explicit-biharmonic CFL
!! ceiling
!!     `ν₄ · dt · ((π/dx)² + (π/dy)²)² ≤ 2`   (project constant),
!! scaled by the MOM6 `bound_coef` safety margin, on top of the static
!! `nu4_max` ceiling.
!!
!! Covers:
!!   1. Fine grid + large `nu_4` (above the CFL limit) → the effective
!!      coefficient is clamped to the CFL bound (the tendency matches a
!!      run whose `nu_4` is set exactly to the analytic bound, on a
!!      uniform grid where every interior face shares the same bound).
!!   2. `nu_4` below the CFL limit → untouched: the tendency is linear
!!      in `nu_4` (halving `nu_4` halves the tendency to round-off),
!!      proving the clamp does not bite (bit-identical merge gate).
module test_ocean_hvisc_bih_cfl
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

   public :: collect_ocean_hvisc_bih_cfl_tests

   real(wp), parameter :: DX = 1000.0_wp     !! fine grid spacing (m)
   real(wp), parameter :: DT = 600.0_wp      !! outer step (s)
   real(wp), parameter :: BOUND_COEF = 0.8_wp  !! default CFL safety margin
   integer, parameter :: NX = 12, NY = 12, NZ = 2

contains

   subroutine collect_ocean_hvisc_bih_cfl_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("bih_cfl_large_nu4_clamped_to_bound", test_large_nu4_clamped), &
                  new_unittest("bih_cfl_small_nu4_untouched", test_small_nu4_untouched) &
                  ]
   end subroutine collect_ocean_hvisc_bih_cfl_tests

   pure function cfl_bound() result(nu4_max_cfl)
      !! Analytic per-face CFL ceiling on `nu_4` for the uniform fine
      !! grid (`idx = idy = 1/DX`), matching `hvisc_nu4_cfl_bound`.
      real(wp) :: nu4_max_cfl
      real(wp) :: idx, k2
      idx = 1.0_wp/DX
      k2 = (PI*idx)*(PI*idx) + (PI*idx)*(PI*idx)
      nu4_max_cfl = BOUND_COEF*2.0_wp/(DT*k2*k2)
   end function cfl_bound

   subroutine setup(grid, ms, metrics, hv)
      !! Build a fine uniform-cartesian state with a non-harmonic
      !! u-field (so ∇⁴u ≠ 0) and device-map everything.
      type(hgrid_t), intent(inout) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_metrics_t), intent(inout) :: metrics
      type(ocean_horizontal_viscosity_t), intent(inout) :: hv
      integer :: i, j, k
      call grid%init(NX, NY, 1, DX, DX)
      ms%nz_ml = NZ
      call ms%init(grid)
      ms%h_layer = 10.0_wp
      ms%v_face_y_layer = 0.0_wp
      ! Non-harmonic u-field: a short-wavelength sinusoid in i has a
      ! non-zero biharmonic, so du_visc is non-trivial.
      do k = 1, NZ
         do j = 1, grid%ny_total
            do i = 1, grid%nx_total + 1
               ms%u_face_x_layer(i, j, k) = sin(2.0_wp*PI*real(i - 1, wp)/real(NX, wp))
            end do
         end do
      end do
      call hv%init(grid, nz_ml=NZ)
      hv%nu_h = 0.0_wp
      hv%bound_coef = BOUND_COEF
      call make_cartesian_metrics(metrics, grid)
      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(hv)
      call hv%enter_data()
   end subroutine setup

   subroutine teardown(grid, ms, metrics, hv)
      type(hgrid_t), intent(inout) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_metrics_t), intent(inout) :: metrics
      type(ocean_horizontal_viscosity_t), intent(inout) :: hv
      call hv%exit_data()
      !$acc exit data delete(hv)
      call ms%exit_data()
      !$acc exit data delete(ms)
      call destroy_cartesian_metrics(metrics)
      call hv%destroy()
      call ms%destroy()
   end subroutine teardown

   subroutine run_scalar_biharmonic(nu_4, du_out)
      !! Run the scalar-`nu_4` biharmonic compute on a fresh state and
      !! return the resulting east-face tendency.
      real(wp), intent(in) :: nu_4
      real(wp), allocatable, intent(out) :: du_out(:, :, :)
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_horizontal_viscosity_t) :: hv
      call setup(grid, ms, metrics, hv)
      hv%nu_4 = nu_4
      call ocean_horizontal_viscosity_compute_tendencies(grid, metrics, hv, ms, dt=DT)
      !$acc update self(hv%du_visc%data)
      du_out = hv%du_visc%data
      call teardown(grid, ms, metrics, hv)
   end subroutine run_scalar_biharmonic

   subroutine test_large_nu4_clamped(error)
      !! A `nu_4` set 100× above the CFL bound must produce the SAME
      !! tendency as `nu_4` set exactly to the bound — on a uniform
      !! grid every interior face clamps to the same value, so the two
      !! tendencies coincide to round-off.
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: du_big(:, :, :), du_lim(:, :, :)
      real(wp) :: bound, max_diff, ref
      bound = cfl_bound()
      call run_scalar_biharmonic(100.0_wp*bound, du_big)
      call run_scalar_biharmonic(bound, du_lim)
      max_diff = maxval(abs(du_big - du_lim))
      ref = maxval(abs(du_lim))
      ! The clamped run must be a strictly non-trivial tendency (so the
      ! test actually exercises the biharmonic) and must match the
      ! at-the-limit run to round-off.
      call check(error, ref > 0.0_wp, "biharmonic tendency must be non-trivial")
      if (allocated(error)) return
      call check(error, max_diff <= 1.0e-9_wp*ref, &
                 "nu_4 above the CFL limit must clamp to the CFL bound")
   end subroutine test_large_nu4_clamped

   subroutine test_small_nu4_untouched(error)
      !! A `nu_4` well below the CFL bound must NOT be clamped: the
      !! biharmonic tendency is linear in `nu_4`, so halving `nu_4`
      !! must halve the tendency to round-off (bit-identical gate).
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: du_a(:, :, :), du_b(:, :, :)
      real(wp) :: bound, nu_a, nu_b, max_diff, ref
      bound = cfl_bound()
      nu_a = 0.25_wp*bound   ! safely below the limit
      nu_b = 0.50_wp*bound   ! still below, exactly 2× nu_a
      call run_scalar_biharmonic(nu_a, du_a)
      call run_scalar_biharmonic(nu_b, du_b)
      ! Unclamped ⇒ du_b == 2 * du_a to round-off.
      max_diff = maxval(abs(du_b - 2.0_wp*du_a))
      ref = maxval(abs(du_b))
      call check(error, ref > 0.0_wp, "biharmonic tendency must be non-trivial")
      if (allocated(error)) return
      call check(error, max_diff <= 1.0e-12_wp*ref, &
                 "nu_4 below the CFL limit must be left untouched (linear in nu_4)")
   end subroutine test_small_nu4_untouched

end module test_ocean_hvisc_bih_cfl
