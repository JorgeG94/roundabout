!! OM3 lateral-viscosity floors and the Coriolis bound on the biharmonic
!! Smagorinsky (`kh_sin_lat`, `ah_vel_scale`, `bound_coriolis_biharm`).
!!
!! Covers (each against an independent closed form):
!!   1. `kh_sin_lat`: face floor = 2000·|sin φ|^4 at 0, 30, 60, -60 deg;
!!      a larger competing floor (`ah_bg`, or a prefilled field) wins.
!!   2. `ah_vel_scale`: uniform dx=dy=10 km gives Ah_bg = 0.01·(1e4)^3 = 1e10;
!!      an anisotropic cell matches 0.01·(Δ²)^(3/2), Δ² the harmonic mean.
!!   3. `bound_coriolis_biharm`: uniform shear matches
!!      |S|·(C·Δ⁴ + Δ⁶·f_max·|S|/(5V²)); f = 0 is bit-identical to off.
!!   4. Configure refusals.
module test_ocean_hvisc_om3_floors
   use rdb_constants, only: wp
   use rdb_config, only: config_t, read_config_from_string, validate_config
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_lateral_mix, only: ocean_lateral_mix_t, &
                                    ocean_lateral_mix_compute_smag_ah, &
                                    ocean_lateral_mix_compute_lat_floor
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_CONFIG_VALIDATE
   use testdrive, only: error_type, check, new_unittest, unittest_type
   implicit none
   private

   public :: collect_ocean_hvisc_om3_floors_tests

contains

   subroutine collect_ocean_hvisc_om3_floors_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("kh_sin_lat_floor_matches_sine_power", test_kh_sin_lat), &
                  new_unittest("ah_vel_scale_uniform_1e10", test_ah_vel_uniform), &
                  new_unittest("ah_vel_scale_anisotropic_harmonic_mean", test_ah_vel_aniso), &
                  new_unittest("bound_coriolis_biharm_matches_formula", test_cor_formula), &
                  new_unittest("bound_coriolis_biharm_f0_bit_identical", test_cor_f0), &
                  new_unittest("om3_floor_knobs_fail_loud", test_refusals) &
                  ]
   end subroutine collect_ocean_hvisc_om3_floors_tests

   subroutine setup_grid_ms(grid, ms, nx, ny, nz, dx, dy)
      type(hgrid_t), intent(inout) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: dx, dy
      call grid%init(nx, ny, 1, dx, dy)
      ms%nz_ml = nz
      call ms%init(grid)
      ms%h_layer = 10.0_wp
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp
   end subroutine setup_grid_ms

   subroutine map_in(ms, metrics, grid, lmix)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_metrics_t), intent(inout) :: metrics
      type(hgrid_t), intent(in) :: grid
      type(ocean_lateral_mix_t), intent(inout) :: lmix
      call make_cartesian_metrics(metrics, grid)
      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(lmix)
      call lmix%enter_data()
   end subroutine map_in

   subroutine map_out(ms, metrics, lmix)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_metrics_t), intent(inout) :: metrics
      type(ocean_lateral_mix_t), intent(inout) :: lmix
      call lmix%exit_data()
      !$acc exit data delete(lmix)
      call ms%exit_data()
      !$acc exit data delete(ms)
      call destroy_cartesian_metrics(metrics)
   end subroutine map_out

   subroutine test_kh_sin_lat(error)
      !! Rows of corners at 0, 30, 60, -60 deg; u-face rows j=1,3,5,7 sit at
      !! exactly those face-centre latitudes.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_lateral_mix_t) :: lmix
      real(wp), allocatable :: lat_bu(:, :), f_bu(:, :)
      real(wp), parameter :: LATS(4) = [0.0_wp, 30.0_wp, 60.0_wp, -60.0_wp]
      real(wp), parameter :: EXPECT(4) = [0.0_wp, 125.0_wp, 1125.0_wp, 1125.0_wp]
      integer, parameter :: NX = 4, NY = 8, NZ = 1
      integer :: n, jrow
      real(wp) :: worst_a, worst_b, worst_c, got

      call setup_grid_ms(grid, ms, NX, NY, NZ, 1.0e4_wp, 1.0e4_wp)
      call lmix%init(grid, nz_ml=NZ)
      allocate (lat_bu(grid%nx_total + 1, grid%ny_total + 1), f_bu(grid%nx_total + 1, grid%ny_total + 1))
      f_bu = 0.0_wp
      do jrow = 1, grid%ny_total + 1
         lat_bu(:, jrow) = LATS(min(4, (jrow + 1)/2))   ! rows (1,2)->0, (3,4)->30, (5,6)->60, rest -60
      end do
      lmix%kh_sin_lat = 2000.0_wp
      lmix%kh_pwr_of_sine = 4.0_wp
      lmix%ah_bg = 0.0_wp
      call lmix%build_statics(lat_bu, f_bu)
      call map_in(ms, metrics, grid, lmix)

      ! (a) nothing else fills the field: the floor itself (ah_bg = 0).
      call ocean_lateral_mix_compute_lat_floor(lmix, seed_bg=.true.)
      !$acc update self(lmix%ah_face_x, lmix%ah_face_y)
      worst_a = 0.0_wp
      do n = 2, 4
         jrow = 2*n - 1
         got = lmix%ah_face_x(2, jrow, 1)
         worst_a = max(worst_a, abs(got - EXPECT(n))/EXPECT(n))
      end do
      call check(error, lmix%ah_face_x(2, 1, 1) == 0.0_wp, "equator floor must be exactly 0")
      if (allocated(error)) return
      call check(error, worst_a < 1.0e-14_wp, "floor must equal 2000*|sin(lat)|^4 (30, 60, -60 deg)")
      if (allocated(error)) return
      ! v-face at the 60 deg row (corner rows 5,6 -> v-face j=5 averages rows 5 and 5).
      got = lmix%ah_face_y(2, 5, 1)
      call check(error, abs(got - 1125.0_wp)/1125.0_wp < 1.0e-14_wp, "v-face floor at 60 deg")
      if (allocated(error)) return

      ! (b) a larger competing floor (ah_bg = 500) wins where it is larger.
      lmix%ah_bg = 500.0_wp
      !$acc update device(lmix%ah_bg)
      call ocean_lateral_mix_compute_lat_floor(lmix, seed_bg=.true.)
      !$acc update self(lmix%ah_face_x)
      worst_b = abs(lmix%ah_face_x(2, 1, 1) - 500.0_wp) + abs(lmix%ah_face_x(2, 3, 1) - 500.0_wp)
      call check(error, worst_b == 0.0_wp, "ah_bg=500 must beat the 0 and 125 floors")
      if (allocated(error)) return
      call check(error, abs(lmix%ah_face_x(2, 5, 1) - 1125.0_wp)/1125.0_wp < 1.0e-14_wp, &
                 "the 1125 floor must beat ah_bg=500")
      if (allocated(error)) return

      ! (c) non-seeded: max against a prefilled 700 field.
      lmix%ah_face_x = 700.0_wp
      lmix%ah_face_y = 700.0_wp
      !$acc update device(lmix%ah_face_x, lmix%ah_face_y)
      call ocean_lateral_mix_compute_lat_floor(lmix, seed_bg=.false.)
      !$acc update self(lmix%ah_face_x)
      worst_c = abs(lmix%ah_face_x(2, 1, 1) - 700.0_wp) + abs(lmix%ah_face_x(2, 3, 1) - 700.0_wp)
      call check(error, worst_c == 0.0_wp .and. &
                 abs(lmix%ah_face_x(2, 7, 1) - 1125.0_wp)/1125.0_wp < 1.0e-14_wp, &
                 "floor must max against the existing field")
      call map_out(ms, metrics, lmix)
      call lmix%destroy()
      call ms%destroy()
   end subroutine test_kh_sin_lat

   subroutine run_ah_vel(nx, ny, dx, dy, ah_vel, nu4_bg, lmix, grid, ms, metrics)
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: dx, dy, ah_vel, nu4_bg
      type(ocean_lateral_mix_t), intent(inout) :: lmix
      type(hgrid_t), intent(inout) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_metrics_t), intent(inout) :: metrics
      call setup_grid_ms(grid, ms, nx, ny, 1, dx, dy)
      call lmix%init(grid, nz_ml=1)
      lmix%smag_ah_active = .true.
      lmix%smag_bi_const = 0.06_wp
      lmix%nu4_bg = nu4_bg
      lmix%nu4_max = 1.0e20_wp
      lmix%ah_vel_scale = ah_vel
      call map_in(ms, metrics, grid, lmix)
      call ocean_lateral_mix_compute_smag_ah(grid, metrics, lmix, ms)
      !$acc update self(lmix%nu4_face_x, lmix%nu4_face_y)
      call map_out(ms, metrics, lmix)
   end subroutine run_ah_vel

   subroutine test_ah_vel_uniform(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_lateral_mix_t) :: lmix
      integer, parameter :: NX = 8, NY = 8
      real(wp) :: rel

      call run_ah_vel(NX, NY, 1.0e4_wp, 1.0e4_wp, 0.01_wp, 0.0_wp, lmix, grid, ms, metrics)
      rel = max(maxval(abs(lmix%nu4_face_x(2:NX, 2:NY - 1, :) - 1.0e10_wp)), &
                maxval(abs(lmix%nu4_face_y(2:NX - 1, 2:NY, :) - 1.0e10_wp)))/1.0e10_wp
      call check(error, rel < 1.0e-14_wp, "0.01*(1e4)^3 must be 1e10 on a 10 km grid")
      if (allocated(error)) return
      call lmix%destroy()
      call ms%destroy()

      ! A larger nu_4_bg wins over the velocity-scale floor.
      call run_ah_vel(NX, NY, 1.0e4_wp, 1.0e4_wp, 0.01_wp, 3.0e10_wp, lmix, grid, ms, metrics)
      rel = maxval(abs(lmix%nu4_face_x(2:NX, 2:NY - 1, :) - 3.0e10_wp))/3.0e10_wp
      call check(error, rel < 1.0e-14_wp, "nu_4_bg=3e10 must beat the 1e10 velocity floor")
      call lmix%destroy()
      call ms%destroy()
   end subroutine test_ah_vel_uniform

   subroutine test_ah_vel_aniso(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_lateral_mix_t) :: lmix
      integer, parameter :: NX = 8, NY = 8
      real(wp), parameter :: DX = 1.0e4_wp, DY = 2.0e4_wp
      real(wp) :: d2, expect, rel

      d2 = 2.0_wp*DX**2*DY**2/(DX**2 + DY**2)   ! 1.6e8 m^2
      expect = 0.01_wp*d2*sqrt(d2)
      call run_ah_vel(NX, NY, DX, DY, 0.01_wp, 0.0_wp, lmix, grid, ms, metrics)
      rel = max(maxval(abs(lmix%nu4_face_x(2:NX, 2:NY - 1, :) - expect)), &
                maxval(abs(lmix%nu4_face_y(2:NX - 1, 2:NY, :) - expect)))/expect
      call check(error, rel < 1.0e-14_wp, "anisotropic cell must follow the harmonic-mean Delta^2")
      call lmix%destroy()
      call ms%destroy()
   end subroutine test_ah_vel_aniso

   subroutine run_cor(f_scale, bound_on, lmix, ms, grid, metrics, nx, ny, dx, shear, vcor)
      real(wp), intent(in) :: f_scale, dx, shear, vcor
      logical, intent(in) :: bound_on
      type(ocean_lateral_mix_t), intent(inout) :: lmix
      type(multilayer_state_t), intent(inout) :: ms
      type(hgrid_t), intent(inout) :: grid
      type(ocean_metrics_t), intent(inout) :: metrics
      integer, intent(in) :: nx, ny
      real(wp), allocatable :: f_bu(:, :), lat_bu(:, :)
      integer :: i, j

      call setup_grid_ms(grid, ms, nx, ny, 1, dx, dx)
      do j = 1, ny
         do i = 1, nx + 1
            ms%u_face_x_layer(i, j, 1) = shear*real(j, wp)
         end do
      end do
      call lmix%init(grid, nz_ml=1)
      lmix%smag_ah_active = .true.
      lmix%smag_bi_const = 0.06_wp
      lmix%nu4_bg = 0.0_wp
      lmix%nu4_max = 1.0e30_wp
      lmix%bound_coriolis_biharm = bound_on
      lmix%bound_coriolis_vel = vcor
      allocate (f_bu(grid%nx_total + 1, grid%ny_total + 1), lat_bu(grid%nx_total + 1, grid%ny_total + 1))
      lat_bu = 0.0_wp
      do j = 1, grid%ny_total + 1
         f_bu(:, j) = f_scale*(1.0_wp + 0.1_wp*real(j, wp))   ! f grows with j: max at the north corner
      end do
      call lmix%build_statics(lat_bu, f_bu)
      call map_in(ms, metrics, grid, lmix)
      call ocean_lateral_mix_compute_smag_ah(grid, metrics, lmix, ms)
      !$acc update self(lmix%nu4_face_x)
      call map_out(ms, metrics, lmix)
   end subroutine run_cor

   subroutine test_cor_formula(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_lateral_mix_t) :: lmix
      real(wp), allocatable :: a_off(:, :), a_on(:, :)
      integer, parameter :: NX = 12, NY = 12
      real(wp), parameter :: DX = 5.0e4_wp, SHEAR = 0.01_wp, VCOR = 6.0_wp, F0 = 1.0e-4_wp
      real(wp), parameter :: CB = 0.06_wp
      real(wp) :: s_mag, d2, expect_on, expect_extra, fmax, worst_on, worst_extra
      integer :: i, j

      s_mag = SHEAR/DX          ! uniform u = SHEAR*j: |D_S| = du/dy = SHEAR/dy
      d2 = DX*DX
      call run_cor(F0, .false., lmix, ms, grid, metrics, NX, NY, DX, SHEAR, VCOR)
      a_off = lmix%nu4_face_x(:, :, 1)
      call lmix%destroy()
      call ms%destroy()
      call run_cor(F0, .true., lmix, ms, grid, metrics, NX, NY, DX, SHEAR, VCOR)
      a_on = lmix%nu4_face_x(:, :, 1)
      call lmix%destroy()
      call ms%destroy()

      worst_on = 0.0_wp
      worst_extra = 0.0_wp
      do j = 3, NY - 2
         do i = 3, NX - 1
            fmax = F0*(1.0_wp + 0.1_wp*real(j + 1, wp))      ! largest of the cell's 4 corners
            expect_on = s_mag*(CB*d2*d2 + d2*d2*d2*fmax*s_mag/(5.0_wp*VCOR*VCOR))
            expect_extra = (d2*d2*d2)*fmax*s_mag*s_mag/(5.0_wp*VCOR*VCOR)
            worst_on = max(worst_on, abs(a_on(i, j) - expect_on)/expect_on)
            worst_extra = max(worst_extra, abs((a_on(i, j) - a_off(i, j)) - expect_extra)/expect_extra)
         end do
      end do
      call check(error, worst_on < 1.0e-13_wp, "Ah_Sm must match |S|(C D^4 + D^6 fmax |S|/(5V^2))")
      if (allocated(error)) return
      call check(error, worst_extra < 1.0e-13_wp, "bounded - unbounded must be exactly the Coriolis term")
      if (allocated(error)) return
      call check(error, minval(a_on(3:NX - 1, 3:NY - 2) - a_off(3:NX - 1, 3:NY - 2)) > 0.0_wp, &
                 "the bound must raise the viscosity")
   end subroutine test_cor_formula

   subroutine test_cor_f0(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_lateral_mix_t) :: lmix
      real(wp), allocatable :: a_off(:, :), a_on(:, :)
      integer, parameter :: NX = 12, NY = 12

      call run_cor(0.0_wp, .false., lmix, ms, grid, metrics, NX, NY, 5.0e4_wp, 0.01_wp, 6.0_wp)
      a_off = lmix%nu4_face_x(:, :, 1)
      call lmix%destroy()
      call ms%destroy()
      call run_cor(0.0_wp, .true., lmix, ms, grid, metrics, NX, NY, 5.0e4_wp, 0.01_wp, 6.0_wp)
      a_on = lmix%nu4_face_x(:, :, 1)
      call lmix%destroy()
      call ms%destroy()
      call check(error, all(a_on == a_off) .and. maxval(a_off) > 0.0_wp, &
                 "f=0 must reproduce the unbounded Smagorinsky biharmonic bit for bit")
   end subroutine test_cor_f0

   subroutine test_refusals(error)
      type(error_type), allocatable, intent(out) :: error
      call expect(error, "&ocean_hvisc_nml smag_ah=.true., bound_coriolis_biharm=.true. /", .true., &
                  "bound_coriolis_biharm with smag_ah")
      if (allocated(error)) return
      call expect(error, "&ocean_hvisc_nml bound_coriolis_biharm=.true. /", .false., &
                  "bound_coriolis_biharm without smag_ah")
      if (allocated(error)) return
      call expect(error, "&ocean_hvisc_nml ah_vel_scale=0.01 /", .false., "ah_vel_scale without smag_ah")
      if (allocated(error)) return
      call expect(error, "&ocean_hvisc_nml smag_ah=.true., bound_coriolis_biharm=.true., "// &
                  "bound_coriolis_vel=0.0 /", .false., "bound_coriolis_vel = 0")
      if (allocated(error)) return
      call expect(error, "&ocean_hvisc_nml kh_sin_lat=2000.0 /", .false., "kh_sin_lat on cartesian")
   end subroutine test_refusals

   subroutine expect(error, nml, valid, what)
      type(error_type), allocatable, intent(inout) :: error
      character(len=*), intent(in) :: nml, what
      logical, intent(in) :: valid
      type(config_t) :: cfg
      integer :: ierr
      call read_config_from_string(nml, cfg, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) then
         call check(error, .false., what//": namelist must parse")
         return
      end if
      ierr = -999
      call validate_config(cfg, ierr=ierr)
      if (valid) then
         call check(error, ierr == OCEAN_STATUS_OK, what//": must be accepted")
      else
         call check(error, ierr == OCEAN_STATUS_ERR_CONFIG_VALIDATE, what//": must be refused")
      end if
   end subroutine expect

end module test_ocean_hvisc_om3_floors
