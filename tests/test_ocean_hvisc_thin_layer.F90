!! Thin-layer treatment of the stress-form viscosity (`bound_thin_layer`).
!!
!! (a) A 2-layer column with one vanished face thickness (two adjacent
!!     zero-thickness cells): Kh at the vanished T is exactly 0, at a T
!!     whose smallest face thickness is half its own it is exactly
!!     0.5 x the unaffected Kh_max, and the corner stress on the vanished
!!     face is exactly 0 (the arithmetic mean would leave it nonzero).
!! (b) `hvisc_hq_harmonic` returns the face thickness for uniform h and
!!     exactly 0 when any one face thickness is 0.
!! (c) Uniform h: knob on == knob off (r = 1, harmonic = arithmetic).
!! (d) Miniature of the OM3 blow-up: a vanished layer next to live ones with
!!     a large Kh stays finite and bounded with the knob on.
module test_ocean_hvisc_thin_layer
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_horizontal_viscosity, only: ocean_horizontal_viscosity_t, &
                                             ocean_horizontal_viscosity_compute_tendencies, &
                                             hvisc_hq_harmonic
   use testdrive, only: error_type, check, new_unittest, unittest_type
   implicit none
   private

   public :: collect_ocean_hvisc_thin_layer_tests

   integer, parameter :: NX = 14, NY = 14
   real(wp), parameter :: DX = 50000.0_wp, DT = 1200.0_wp

contains

   subroutine collect_ocean_hvisc_thin_layer_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("hq_harmonic_values", test_hq_values), &
                  new_unittest("vanished_face_kh_and_hq", test_vanished_face), &
                  new_unittest("uniform_h_knob_on_equals_off", test_uniform), &
                  new_unittest("vanished_patch_stays_finite", test_blowup_miniature) &
                  ]
   end subroutine collect_ocean_hvisc_thin_layer_tests

   subroutine run_case(patch, thin, nu, ah_t, str_xy, du, dv)
      !! 2-layer column, layer 1 optionally carrying a vanished patch
      !! (cells i=6:7, j=7), a swirling velocity field, stress path.
      logical, intent(in) :: patch, thin
      real(wp), intent(in) :: nu
      real(wp), allocatable, intent(out) :: ah_t(:, :, :), str_xy(:, :, :), du(:, :, :), dv(:, :, :)
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_horizontal_viscosity_t) :: hv
      integer :: i, j, k, nxt, nyt

      call grid%init(NX, NY, 1, DX, DX)
      nxt = grid%nx_total
      nyt = grid%ny_total
      ms%nz_ml = 2
      call ms%init(grid)
      ms%h_layer = 10.0_wp
      if (patch) then
         ms%h_layer(6, 7, 1) = 0.0_wp
         ms%h_layer(7, 7, 1) = 0.0_wp
      end if
      do k = 1, 2
         do j = 1, nyt
            do i = 1, nxt + 1
               ms%u_face_x_layer(i, j, k) = 0.3_wp*sin(0.9_wp*real(i, wp) + 1.3_wp*real(j, wp))
            end do
         end do
         do j = 1, nyt + 1
            do i = 1, nxt
               ms%v_face_y_layer(i, j, k) = 0.3_wp*cos(1.1_wp*real(i, wp) - 0.7_wp*real(j, wp))
            end do
         end do
      end do
      call hv%init(grid, nz_ml=2)
      hv%nu_h = nu
      hv%stress_tensor = .true.
      hv%bound_thin_layer = thin
      call make_cartesian_metrics(metrics, grid)
      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(hv)
      call hv%enter_data()
      call ocean_horizontal_viscosity_compute_tendencies(grid, metrics, hv, ms, dt=DT)
      !$acc update self(hv%ah_t%data, hv%str_xy%data, hv%du_visc%data, hv%dv_visc%data)
      allocate (ah_t, source=hv%ah_t%data)
      allocate (str_xy, source=hv%str_xy%data)
      allocate (du, source=hv%du_visc%data)
      allocate (dv, source=hv%dv_visc%data)
      call hv%exit_data()
      !$acc exit data delete(hv)
      call ms%exit_data()
      !$acc exit data delete(ms)
      call destroy_cartesian_metrics(metrics)
      call hv%destroy()
      call ms%destroy()
   end subroutine run_case

   subroutine test_hq_values(error)
      type(error_type), allocatable, intent(out) :: error
      call check(error, abs(hvisc_hq_harmonic(7.0_wp, 7.0_wp, 7.0_wp, 7.0_wp) - 7.0_wp) <= 1.0e-14_wp, &
                 "uniform face thickness must give h_q = h")
      if (allocated(error)) return
      call check(error, hvisc_hq_harmonic(0.0_wp, 7.0_wp, 7.0_wp, 7.0_wp) == 0.0_wp, "zero huS => h_q = 0")
      if (allocated(error)) return
      call check(error, hvisc_hq_harmonic(7.0_wp, 7.0_wp, 7.0_wp, 0.0_wp) == 0.0_wp, "zero hvE => h_q = 0")
      if (allocated(error)) return
      call check(error, hvisc_hq_harmonic(7.0_wp, 0.0_wp, 3.0_wp, 5.0_wp) == 0.0_wp, "zero huN => h_q = 0")
   end subroutine test_hq_values

   subroutine test_vanished_face(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: at_on(:, :, :), sq_on(:, :, :), du(:, :, :), dv(:, :, :)
      real(wp), allocatable :: at_off(:, :, :), sq_off(:, :, :)
      real(wp) :: kmax

      call run_case(.true., .true., 1.0e9_wp, at_on, sq_on, du, dv)
      call run_case(.true., .false., 1.0e9_wp, at_off, sq_off, du, dv)
      kmax = at_off(3, 3, 1)       ! far from the patch: the plain CFL bound
      call check(error, kmax > 0.0_wp .and. at_on(3, 3, 1) == kmax, "far T must keep the full Kh_max")
      if (allocated(error)) return
      ! T (7,7) is vanished, T (8,7) has face thicknesses (5, 10, 10, 10).
      call check(error, at_on(7, 7, 1) == 0.0_wp, "vanished T: Kh = r*Kh_max = 0 exactly")
      if (allocated(error)) return
      call check(error, abs(at_on(8, 7, 1) - 0.5_wp*kmax) <= 1.0e-14_wp*kmax, &
                 "T with half-thickness face: Kh = 0.5*Kh_max")
      if (allocated(error)) return
      call check(error, at_off(7, 7, 1) == kmax, "knob off: no thin-layer ratio")
      if (allocated(error)) return
      ! Corner (7,8) has huS = h_u(7,7) = 0 (cells 6,7 of row 7 are vanished).
      call check(error, sq_on(7, 8, 1) == 0.0_wp, "corner stress on the vanished face is exactly 0")
      if (allocated(error)) return
      call check(error, abs(sq_off(7, 8, 1)) > 0.0_wp, "knob off: arithmetic h_q leaves a nonzero corner stress")
   end subroutine test_vanished_face

   subroutine test_uniform(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: at_on(:, :, :), sq_on(:, :, :), du_on(:, :, :), dv_on(:, :, :)
      real(wp), allocatable :: at_off(:, :, :), sq_off(:, :, :), du_off(:, :, :), dv_off(:, :, :)
      call run_case(.false., .true., 1.0e9_wp, at_on, sq_on, du_on, dv_on)
      call run_case(.false., .false., 1.0e9_wp, at_off, sq_off, du_off, dv_off)
      call check(error, maxval(abs(at_on - at_off)) == 0.0_wp, "uniform h: r = 1, Kh identical")
      if (allocated(error)) return
      call check(error, maxval(abs(du_on - du_off)) <= 1.0e-12_wp*maxval(abs(du_off)) .and. &
                 maxval(abs(dv_on - dv_off)) <= 1.0e-12_wp*maxval(abs(dv_off)), &
                 "uniform h: harmonic h_q reproduces the arithmetic tendency")
   end subroutine test_uniform

   subroutine test_blowup_miniature(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: at(:, :, :), sq(:, :, :), du_on(:, :, :), dv_on(:, :, :)
      real(wp), allocatable :: du_off(:, :, :), dv_off(:, :, :)
      real(wp) :: worst_on, worst_off
      call run_case(.true., .true., 1.0e9_wp, at, sq, du_on, dv_on)
      call run_case(.true., .false., 1.0e9_wp, at, sq, du_off, dv_off)
      worst_on = DT*max(maxval(abs(du_on)), maxval(abs(dv_on)))
      worst_off = DT*max(maxval(abs(du_off)), maxval(abs(dv_off)))
      call check(error, worst_on == worst_on .and. worst_on < 1.0_wp, &
                 "knob on: velocity change per step stays bounded and finite near a vanished patch")
      if (allocated(error)) return
      call check(error, worst_off > 10.0_wp*worst_on, &
                 "positive control: knob off, the same patch gives a much larger (unstable) step")
   end subroutine test_blowup_miniature

end module test_ocean_hvisc_thin_layer
