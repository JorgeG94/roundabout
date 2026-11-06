!! Unit tests for the reduced-gravity / gprime PGF kernel
!! (`OPGF_VARIANT_GPRIME`).  MOM6 COORD_CONFIG="gprime" analogue,
!! restricted to NK = 2 in Tier-1.
!!
!! Convention reminder: k = 1 is the bottom layer (heavier, ρ₂),
!! k = NZ = 2 is the surface layer (lighter, ρ₁).
!!
!!     a_top = -g_FS · ∂η/∂x,                  η = h_1 + h_2 - b
!!     a_bot = -g_FS · ∂η/∂x - g_int · ∂e/∂x,    e = h_1 - b (interface height)
!!
!! Cases:
!!   * Uniform h_layer → zero PGF in both layers.
!!   * Linear surface tilt with uniform bottom layer thickness →
!!     constant `-g_FS · slope` in BOTH layers, no internal mode.
!!   * Linear interface tilt with uniform total depth → BC mode only,
!!     top layer untouched, bottom layer gets `-g_int · slope`.
!!   * Sloping bed, flat surface AND flat interface (the layered resting
!!     state over topography) → zero PGF in both layers.  Every case above
!!     has a flat bed, so a kernel that differenced `h_1` instead of the
!!     interface height `h_1 - b` passed them all while pushing the
!!     double-gyre spoon's bed layer with a spurious `g_int · ∇b`.
!!   * Sloping bed under a tilted interface → `-g_int · (interface slope)`,
!!     not `-g_int · ∂h_1/∂x`.
module test_ocean_pgf_gprime
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_ocean_pressure_force, only: ocean_pressure_force_t, &
                                       ocean_pressure_force_compute, &
                                       OPGF_VARIANT_GPRIME
   implicit none
   private

   public :: collect_ocean_pgf_gprime_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 2

contains

   subroutine collect_ocean_pgf_gprime_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("gprime_uniform_rest_zero_pgf", test_uniform_rest), &
                  new_unittest("gprime_surface_tilt_pure_bt", test_surface_tilt), &
                  new_unittest("gprime_interface_tilt_pure_bc", test_interface_tilt), &
                  new_unittest("gprime_sloping_bed_flat_interface_rest", test_sloping_bed_rest), &
                  new_unittest("gprime_sloping_bed_interface_tilt", test_sloping_bed_interface_tilt) &
                  ]
   end subroutine collect_ocean_pgf_gprime_tests

   subroutine make_grid(grid, nx_phys, ny_phys, dx)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      real(wp), intent(in) :: dx
      call grid%init(nx_phys, ny_phys, NGHOST, dx, dx)
   end subroutine make_grid

   subroutine setup_pgf(grid, ms, pgf, gfs, gint)
      type(hgrid_t), intent(inout) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_pressure_force_t), intent(inout) :: pgf
      real(wp), intent(in) :: gfs, gint

      ms%nz_ml = NZ
      call ms%init(grid)
      call pgf%init(grid, nz_ml=NZ)
      pgf%variant = OPGF_VARIANT_GPRIME
      pgf%gprime_gfs = gfs
      pgf%gprime_gint = gint
   end subroutine setup_pgf

   subroutine map_in(grid, ms, pgf, metrics)
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_pressure_force_t), intent(inout) :: pgf
      type(ocean_metrics_t), intent(inout) :: metrics
      call make_cartesian_metrics(metrics, grid)
      !$omp target enter data map(to: ms)
      call ms%enter_data()
      !$omp target enter data map(to: pgf)
      call pgf%enter_data()
   end subroutine map_in

   subroutine map_out(ms, pgf, metrics)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_pressure_force_t), intent(inout) :: pgf
      type(ocean_metrics_t), intent(inout) :: metrics
      call pgf%exit_data()
      !$omp target exit data map(delete: pgf)
      call ms%exit_data()
      !$omp target exit data map(delete: ms)
      call destroy_cartesian_metrics(metrics)
   end subroutine map_out

   subroutine test_uniform_rest(error)
      !! Uniform h_layer ⇒ zero PGF everywhere.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_pressure_force_t) :: pgf
      type(ocean_metrics_t) :: metrics
      real(wp) :: max_dpdx, max_dpdy

      call make_grid(grid, 6, 6, 1000.0_wp)
      call setup_pgf(grid, ms, pgf, gfs=9.81_wp, gint=0.01_wp)
      ms%h_layer = 100.0_wp

      call map_in(grid, ms, pgf, metrics)
      call ocean_pressure_force_compute(grid, metrics, pgf, ms)
      !$omp target update from(pgf%dpdx_face%data, pgf%dpdy_face%data)
      call map_out(ms, pgf, metrics)

      max_dpdx = maxval(abs(pgf%dpdx_face%data))
      max_dpdy = maxval(abs(pgf%dpdy_face%data))

      call check(error, max_dpdx < 1.0e-14_wp .and. max_dpdy < 1.0e-14_wp, &
                 "Uniform h_layer must produce zero gprime PGF")

      call pgf%destroy(); call ms%destroy()
   end subroutine test_uniform_rest

   subroutine test_surface_tilt(error)
      !! Surface tilts uniformly across the basin; h_1 (bottom) is
      !! constant; h_2 (top) varies linearly with x.  Both layers
      !! feel `-g_FS · slope`; internal interface is flat ⇒ no BC
      !! contribution to the bottom layer.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_pressure_force_t) :: pgf
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: DX = 1000.0_wp, GFS = 9.81_wp, GINT = 0.01_wp
      real(wp), parameter :: SLOPE = 1.0e-5_wp  ! dh/dx
      real(wp) :: a_expected, a_top_obs, a_bot_obs
      integer :: i, j, nx, ny, i_int, j_int

      call make_grid(grid, 6, 6, DX)
      call setup_pgf(grid, ms, pgf, gfs=GFS, gint=GINT)
      nx = grid%nx_total
      ny = grid%ny_total

      ! h_1 (bottom) uniform; h_2 (top) varies linearly in x.
      ms%h_layer(:, :, 1) = 100.0_wp
      do j = 1, ny
         do i = 1, nx
            ms%h_layer(i, j, 2) = 100.0_wp + SLOPE*real(i, wp)*DX
         end do
      end do

      call map_in(grid, ms, pgf, metrics)
      call ocean_pressure_force_compute(grid, metrics, pgf, ms)
      !$omp target update from(pgf%dpdx_face%data, pgf%dpdy_face%data)
      call map_out(ms, pgf, metrics)

      a_expected = -GFS*SLOPE
      ! Pick an interior u-face well away from walls.
      i_int = nx/2
      j_int = ny/2
      a_top_obs = pgf%dpdx_face%data(i_int, j_int, 2)
      a_bot_obs = pgf%dpdx_face%data(i_int, j_int, 1)

      call check(error, abs(a_top_obs - a_expected) < 1.0e-12_wp, &
                 "Surface tilt: top layer u-face PGF off analytic")
      if (.not. allocated(error)) call check(error, abs(a_bot_obs - a_expected) < 1.0e-12_wp, &
                                             "Surface tilt: bottom layer must also see -g_FS · slope")

      call pgf%destroy(); call ms%destroy()
   end subroutine test_surface_tilt

   subroutine test_interface_tilt(error)
      !! Interface tilts; total depth h_1 + h_2 stays constant.  Top
      !! layer feels zero PGF (∇SSH = 0); bottom layer feels
      !! `-g_int · slope` from the interface gradient.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_pressure_force_t) :: pgf
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: DX = 1000.0_wp, GFS = 9.81_wp, GINT = 0.01_wp
      real(wp), parameter :: SLOPE = 1.0e-5_wp
      real(wp) :: a_top_obs, a_bot_obs, h_bot_x
      integer :: i, j, nx, ny, i_int, j_int

      call make_grid(grid, 6, 6, DX)
      call setup_pgf(grid, ms, pgf, gfs=GFS, gint=GINT)
      nx = grid%nx_total
      ny = grid%ny_total

      ! h_1 + h_2 = 200 (constant).  h_1 = 100 + slope·x, so h_2 =
      ! 100 - slope·x.  Top layer's SSH gradient cancels; bottom
      ! layer feels just the interface tilt.
      do j = 1, ny
         do i = 1, nx
            h_bot_x = 100.0_wp + SLOPE*real(i, wp)*DX
            ms%h_layer(i, j, 1) = h_bot_x
            ms%h_layer(i, j, 2) = 200.0_wp - h_bot_x
         end do
      end do

      call map_in(grid, ms, pgf, metrics)
      call ocean_pressure_force_compute(grid, metrics, pgf, ms)
      !$omp target update from(pgf%dpdx_face%data, pgf%dpdy_face%data)
      call map_out(ms, pgf, metrics)

      i_int = nx/2
      j_int = ny/2
      a_top_obs = pgf%dpdx_face%data(i_int, j_int, 2)
      a_bot_obs = pgf%dpdx_face%data(i_int, j_int, 1)

      call check(error, abs(a_top_obs) < 1.0e-12_wp, &
                 "Interface tilt: top layer should see no PGF")
      if (.not. allocated(error)) call check(error, abs(a_bot_obs - (-GINT*SLOPE)) < 1.0e-12_wp, &
                                             "Interface tilt: bottom layer should feel -g_int · slope")

      call pgf%destroy(); call ms%destroy()
   end subroutine test_interface_tilt

   subroutine test_sloping_bed_rest(error)
      !! Bed slopes in x AND y (`b = 500 + SX·x + SY·y`); the surface
      !! (`η = 0`) and the interface (`e = -300`) are flat, so
      !! `h_2 = 300` and `h_1 = b - 300`.  At rest over topography the
      !! layered PGF is exactly zero in both layers, on every face.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_pressure_force_t) :: pgf
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: DX = 1000.0_wp, GFS = 9.81_wp, GINT = 0.01_wp
      real(wp), parameter :: SX = 4.0e-3_wp, SY = -2.0e-3_wp  ! db/dx, db/dy
      real(wp), allocatable :: b(:, :)
      real(wp) :: max_dpdx, max_dpdy
      integer :: i, j, nx, ny

      call make_grid(grid, 6, 6, DX)
      call setup_pgf(grid, ms, pgf, gfs=GFS, gint=GINT)
      nx = grid%nx_total
      ny = grid%ny_total
      allocate (b(nx, ny))
      do j = 1, ny
         do i = 1, nx
            b(i, j) = 500.0_wp + SX*real(i, wp)*DX + SY*real(j, wp)*DX
            ms%h_layer(i, j, 2) = 300.0_wp
            ms%h_layer(i, j, 1) = b(i, j) - 300.0_wp
         end do
      end do
      call pgf%set_bathymetry(b)

      call map_in(grid, ms, pgf, metrics)
      call ocean_pressure_force_compute(grid, metrics, pgf, ms)
      !$omp target update from(pgf%dpdx_face%data, pgf%dpdy_face%data)
      call map_out(ms, pgf, metrics)

      ! A spurious `g_int · ∇b` would be ~4e-5 here; round-off on O(1e3)
      ! thicknesses differenced over 1 km is ~1e-16.
      max_dpdx = maxval(abs(pgf%dpdx_face%data))
      max_dpdy = maxval(abs(pgf%dpdy_face%data))
      call check(error, max_dpdx < 1.0e-14_wp .and. max_dpdy < 1.0e-14_wp, &
                 "Sloping bed, flat surface + interface: gprime PGF must be zero")

      call pgf%destroy()
      call ms%destroy()
   end subroutine test_sloping_bed_rest

   subroutine test_sloping_bed_interface_tilt(error)
      !! Bed slopes (`db/dx = SB`), surface flat, interface tilted
      !! (`e = -300 + SE·x`, i.e. `de/dx = SE`).  Top layer: zero.
      !! Bottom layer: `-g_int · SE` — the INTERFACE slope, not
      !! `-g_int · ∂h_1/∂x = -g_int · (SE + SB)`.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_pressure_force_t) :: pgf
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: DX = 1000.0_wp, GFS = 9.81_wp, GINT = 0.01_wp
      real(wp), parameter :: SB = 4.0e-3_wp, SE = 1.0e-5_wp
      real(wp), allocatable :: b(:, :)
      real(wp) :: a_top_obs, a_bot_obs, x
      integer :: i, j, nx, ny, i_int, j_int

      call make_grid(grid, 6, 6, DX)
      call setup_pgf(grid, ms, pgf, gfs=GFS, gint=GINT)
      nx = grid%nx_total
      ny = grid%ny_total
      allocate (b(nx, ny))
      do j = 1, ny
         do i = 1, nx
            x = real(i, wp)*DX
            b(i, j) = 500.0_wp + SB*x
            ms%h_layer(i, j, 2) = 300.0_wp - SE*x        ! surface at 0, interface at -300 + SE·x
            ms%h_layer(i, j, 1) = b(i, j) - ms%h_layer(i, j, 2)
         end do
      end do
      call pgf%set_bathymetry(b)

      call map_in(grid, ms, pgf, metrics)
      call ocean_pressure_force_compute(grid, metrics, pgf, ms)
      !$omp target update from(pgf%dpdx_face%data, pgf%dpdy_face%data)
      call map_out(ms, pgf, metrics)

      i_int = nx/2
      j_int = ny/2
      a_top_obs = pgf%dpdx_face%data(i_int, j_int, 2)
      a_bot_obs = pgf%dpdx_face%data(i_int, j_int, 1)
      call check(error, abs(a_top_obs) < 1.0e-12_wp, &
                 "Sloping bed + interface tilt: top layer should see no PGF")
      if (.not. allocated(error)) call check(error, abs(a_bot_obs - (-GINT*SE)) < 1.0e-12_wp, &
                                             "Sloping bed + interface tilt: bottom layer should feel "// &
                                             "-g_int · (interface slope), not -g_int · dh_1/dx")

      call pgf%destroy()
      call ms%destroy()
   end subroutine test_sloping_bed_interface_tilt

end module test_ocean_pgf_gprime
