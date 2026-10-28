!! Static Laplacian floors AT the T and corner points (stress-tensor path,
!! `&ocean_hvisc_nml kh_bg_at_points`).
!!
!! (a) Point values: on a stretched (dx /= dy) lon-lat sector the built
!!     `kh_bg_t`/`kh_bg_q` equal `max(nu_h, c_K*Delta, K_s*|sin phi|^p)`
!!     recomputed here from the point's own latitude and spacings, with an
!!     independent form of the harmonic length, to 1e-14 relative — at a
!!     T row and at a corner row sitting on the equator (sin phi = 0).
!!     The face-mean-then-max form differs from the point value, which is
!!     asserted too (the reason the knob exists).
!! (b) Seed: a floor below `nu_h` is a no-op (bit-identical to knob off); a
!!     floor above `nu_h` raises the tendency.
!! (c) Solid-body rotation u = -Omega*y, v = Omega*x with a strongly varying
!!     Kh field gives zero viscous tendency: the stress form never
!!     differentiates Kh.  A deformation flow on the same field does not.
module test_ocean_hvisc_kh_points
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, make_spherical_metrics, &
                                 destroy_cartesian_metrics
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_horizontal_viscosity, only: ocean_horizontal_viscosity_t, &
                                             ocean_horizontal_viscosity_compute_tendencies
   use testdrive, only: error_type, check, new_unittest, unittest_type
   implicit none
   private

   public :: collect_ocean_hvisc_kh_points_tests

   real(wp), parameter :: PI_L = 3.14159265358979323846264338327950288_wp

contains

   subroutine collect_ocean_hvisc_kh_points_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("point_values_q_equator", test_points_q_equator), &
                  new_unittest("point_values_t_equator", test_points_t_equator), &
                  new_unittest("floor_below_nu_h_is_noop", test_floor_noop), &
                  new_unittest("solid_body_rotation_zero_tendency", test_solid_body) &
                  ]
   end subroutine collect_ocean_hvisc_kh_points_tests

   subroutine point_case(error, lat_south_in_dlat)
      !! Build the point statics on a stretched lon-lat sector and compare.
      type(error_type), allocatable, intent(out) :: error
      real(wp), intent(in) :: lat_south_in_dlat
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ocean_horizontal_viscosity_t) :: hv
      real(wp), parameter :: DLON = 0.5_wp, DLAT = 2.0_wp, CK = 0.01_wp, KS = 2000.0_wp
      real(wp), parameter :: PW = 4.0_wp, NU = 3.0_wp
      integer, parameter :: NX = 8, NY = 12
      integer :: i, j, iq, jq, it, jt, nxt, nyt
      real(wp) :: exp_v, dlt, worst, face_form
      logical :: found_eq

      call grid%init(NX, NY, 1, 1.0_wp, 1.0_wp)
      nxt = grid%nx_total
      nyt = grid%ny_total
      call make_spherical_metrics(metrics, grid, 10.0_wp, lat_south_in_dlat*DLAT, DLON, DLAT, 6.371e6_wp)
      call hv%init(grid, nz_ml=1)
      hv%nu_h = NU
      call hv%build_kh_bg_points(metrics%geolatT, metrics%geolatBu, metrics%dx2h, metrics%dy2h, &
                                 metrics%dx2q, metrics%dy2q, CK, KS, PW)

      worst = 0.0_wp
      do j = 1, nyt
         do i = 1, nxt
            dlt = 1.0_wp/sqrt(0.5_wp*(1.0_wp/metrics%dxT(i, j)**2 + 1.0_wp/metrics%dyT(i, j)**2))
            exp_v = max(NU, CK*dlt, KS*abs(sin(metrics%geolatT(i, j)*PI_L/180.0_wp))**PW)
            worst = max(worst, abs(hv%kh_bg_t(i, j) - exp_v)/exp_v)
         end do
      end do
      do j = 1, nyt + 1
         do i = 1, nxt + 1
            dlt = 1.0_wp/sqrt(0.5_wp*(1.0_wp/metrics%dxBu(i, j)**2 + 1.0_wp/metrics%dyBu(i, j)**2))
            exp_v = max(NU, CK*dlt, KS*abs(sin(metrics%geolatBu(i, j)*PI_L/180.0_wp))**PW)
            worst = max(worst, abs(hv%kh_bg_q(i, j) - exp_v)/exp_v)
         end do
      end do
      call check(error, worst <= 1.0e-14_wp, "kh_bg_t/q must equal the point-wise recomputation to 1e-14")
      if (allocated(error)) return

      ! Equator row: the sin floor is zero, so the value is the other floors.
      found_eq = .false.
      do j = 1, nyt + 1
         if (abs(metrics%geolatBu(nxt/2, j)) < 1.0e-10_wp) then
            found_eq = .true.
            jq = j
            iq = nxt/2
            dlt = 1.0_wp/sqrt(0.5_wp*(1.0_wp/metrics%dxBu(iq, jq)**2 + 1.0_wp/metrics%dyBu(iq, jq)**2))
            call check(error, abs(hv%kh_bg_q(iq, jq) - max(NU, CK*dlt)) <= 1.0e-14_wp*max(NU, CK*dlt), &
                       "equator corner: sin floor must vanish")
            if (allocated(error)) return
         end if
      end do
      do j = 1, nyt
         if (abs(metrics%geolatT(nxt/2, j)) < 1.0e-10_wp) then
            found_eq = .true.
            jt = j
            it = nxt/2
            dlt = 1.0_wp/sqrt(0.5_wp*(1.0_wp/metrics%dxT(it, jt)**2 + 1.0_wp/metrics%dyT(it, jt)**2))
            call check(error, abs(hv%kh_bg_t(it, jt) - max(NU, CK*dlt)) <= 1.0e-14_wp*max(NU, CK*dlt), &
                       "equator T cell: sin floor must vanish")
            if (allocated(error)) return
         end if
      end do
      call check(error, found_eq, "test grid must contain an equator row")
      if (allocated(error)) return

      ! Convexity: the face-averaged form (two u-faces at the T latitude, two
      ! v-faces at the bounding corner latitudes) differs from the point value.
      j = nyt - 2
      face_form = 0.25_wp*(2.0_wp*KS*abs(sin(metrics%geolatT(nxt/2, j)*PI_L/180.0_wp))**PW &
                           + KS*abs(sin(metrics%geolatBu(nxt/2, j)*PI_L/180.0_wp))**PW &
                           + KS*abs(sin(metrics%geolatBu(nxt/2, j + 1)*PI_L/180.0_wp))**PW)
      exp_v = KS*abs(sin(metrics%geolatT(nxt/2, j)*PI_L/180.0_wp))**PW
      call check(error, abs(face_form - exp_v) > 1.0e-3_wp*exp_v .and. exp_v > 0.0_wp, &
                 "face-averaged form must differ from the point value (convex floor)")

      call hv%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine point_case

   subroutine test_points_q_equator(error)
      type(error_type), allocatable, intent(out) :: error
      call point_case(error, -5.0_wp)    ! corner row j=6 on the equator
   end subroutine test_points_q_equator

   subroutine test_points_t_equator(error)
      type(error_type), allocatable, intent(out) :: error
      call point_case(error, -4.5_wp)    ! T row on the equator
   end subroutine test_points_t_equator

   subroutine run_stress(du, mode, kh_lo, kh_hi, use_knob)
      !! Cartesian stress-tensor tendency for `mode` ('solid' | 'deform')
      !! with point statics varying from kh_lo to kh_hi across the domain.
      real(wp), allocatable, intent(out) :: du(:, :, :)
      character(len=*), intent(in) :: mode
      real(wp), intent(in) :: kh_lo, kh_hi
      logical, intent(in) :: use_knob
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_horizontal_viscosity_t) :: hv
      real(wp), parameter :: DX = 50000.0_wp, OM = 1.0e-5_wp, DT = 1200.0_wp, NU = 100.0_wp
      integer, parameter :: NX = 14, NY = 14
      integer :: i, j, nxt, nyt
      real(wp), allocatable :: geo_t(:, :), geo_q(:, :)

      call grid%init(NX, NY, 1, DX, DX)
      nxt = grid%nx_total
      nyt = grid%ny_total
      ms%nz_ml = 1
      call ms%init(grid)
      ms%h_layer = 10.0_wp
      do j = 1, nyt
         do i = 1, nxt + 1
            if (mode == 'solid') then
               ms%u_face_x_layer(i, j, 1) = -OM*(real(j, wp) - 0.5_wp)*DX
            else
               ms%u_face_x_layer(i, j, 1) = 0.02_wp*real(i, wp)**2
            end if
         end do
      end do
      do j = 1, nyt + 1
         do i = 1, nxt
            if (mode == 'solid') then
               ms%v_face_y_layer(i, j, 1) = OM*(real(i, wp) - 0.5_wp)*DX
            else
               ms%v_face_y_layer(i, j, 1) = 0.0_wp
            end if
         end do
      end do
      call hv%init(grid, nz_ml=1)
      hv%nu_h = NU
      hv%stress_tensor = .true.
      call make_cartesian_metrics(metrics, grid)
      if (use_knob) then
         ! Synthetic strongly varying floors (kh_lo..kh_hi), set directly:
         ! latitude fields are not needed for the stencil check.
         allocate (geo_t(nxt, nyt), geo_q(nxt + 1, nyt + 1))
         geo_t = 0.0_wp
         geo_q = 0.0_wp
         call hv%build_kh_bg_points(geo_t, geo_q, metrics%dx2h, metrics%dy2h, metrics%dx2q, metrics%dy2q, &
                                    0.0_wp, 0.0_wp, 4.0_wp)
         do j = 1, nyt
            do i = 1, nxt
               hv%kh_bg_t(i, j) = kh_lo + (kh_hi - kh_lo)*real(i*j, wp)/real(nxt*nyt, wp)
            end do
         end do
         do j = 1, nyt + 1
            do i = 1, nxt + 1
               hv%kh_bg_q(i, j) = kh_hi - (kh_hi - kh_lo)*real(i*j, wp)/real((nxt + 1)*(nyt + 1), wp)
            end do
         end do
      end if
      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(hv)
      call hv%enter_data()
      call ocean_horizontal_viscosity_compute_tendencies(grid, metrics, hv, ms, dt=DT)
      !$acc update self(hv%du_visc%data)
      allocate (du, source=hv%du_visc%data)
      call hv%exit_data()
      !$acc exit data delete(hv)
      call ms%exit_data()
      !$acc exit data delete(ms)
      call destroy_cartesian_metrics(metrics)
      call hv%destroy()
      call ms%destroy()
   end subroutine run_stress

   subroutine test_floor_noop(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: d_off(:, :, :), d_lo(:, :, :), d_hi(:, :, :)
      call run_stress(d_off, 'deform', 0.0_wp, 0.0_wp, .false.)
      call run_stress(d_lo, 'deform', 1.0_wp, 50.0_wp, .true.)    ! all below nu_h = 100
      call run_stress(d_hi, 'deform', 1000.0_wp, 5000.0_wp, .true.)
      call check(error, maxval(abs(d_lo - d_off)) == 0.0_wp, &
                 "a point floor below nu_h must leave the tendency bit-identical")
      if (allocated(error)) return
      call check(error, maxval(abs(d_hi)) > 2.0_wp*maxval(abs(d_off)), &
                 "a point floor above nu_h must raise the viscous tendency")
   end subroutine test_floor_noop

   subroutine test_solid_body(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: d_solid(:, :, :), d_def(:, :, :)
      call run_stress(d_solid, 'solid', 100.0_wp, 1.0e5_wp, .true.)
      call run_stress(d_def, 'deform', 100.0_wp, 1.0e5_wp, .true.)
      call check(error, maxval(abs(d_solid(5:size(d_solid, 1) - 4, 5:size(d_solid, 2) - 4, :))) <= 1.0e-15_wp, &
                 "solid-body rotation must give zero tendency under a strongly varying Kh")
      if (allocated(error)) return
      call check(error, maxval(abs(d_def(5:size(d_def, 1) - 4, 5:size(d_def, 2) - 4, :))) > 1.0e-7_wp, &
                 "positive control: a deformation flow must give a nonzero tendency")
   end subroutine test_solid_body

end module test_ocean_hvisc_kh_points
