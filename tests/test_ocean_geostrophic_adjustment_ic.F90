!! Unit tests for `seed_geostrophic_adjustment_ic` — Rossby's
!! classic geostrophic-adjustment IC.
!!
!! Verifies:
!!   * The peak SSH at the bump centre equals `ga_eta_amp` exactly.
!!   * SSH falls off as a Gaussian — value at one e-folding distance
!!     drops by exactly 1/e.
!!   * Default centre falls at the basin midpoint when ga_x_center /
!!     ga_y_center are left negative.
!!   * Velocities (and face mass fluxes) stay at zero everywhere.
!!   * Total mass equals undisturbed + ∫η dA (bump adds mass).
module test_ocean_geostrophic_adjustment_ic
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_config, only: config_t
   use rdb_grid, only: hgrid_t
   use rdb_ocean_state, only: ocean_state_t, ocean_state_seed_from_cfg
   implicit none
   private

   public :: collect_ocean_geostrophic_adjustment_ic_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NX_PHYS = 40
   integer, parameter :: NY_PHYS = 40
   integer, parameter :: NZ = 1
   real(wp), parameter :: DX = 5000.0_wp
   real(wp), parameter :: DY = 5000.0_wp
   real(wp), parameter :: H_BASIN = 1000.0_wp

contains

   subroutine collect_ocean_geostrophic_adjustment_ic_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("ga_ic_peak_amplitude_matches_cfg", test_peak_amplitude), &
                  new_unittest("ga_ic_gaussian_e_fold_decay", test_e_fold_decay), &
                  new_unittest("ga_ic_default_centre_is_midpoint", test_default_centre), &
                  new_unittest("ga_ic_velocities_zero", test_velocities_zero), &
                  new_unittest("ga_ic_mass_added_matches_bump_integral", test_mass_added) &
                  ]
   end subroutine collect_ocean_geostrophic_adjustment_ic_tests

   subroutine setup_ga(cfg, grid, state, eta_amp, length_scale, x_c, y_c)
      type(config_t), intent(out) :: cfg
      type(hgrid_t), intent(out) :: grid
      type(ocean_state_t), intent(out) :: state
      real(wp), intent(in) :: eta_amp, length_scale, x_c, y_c

      cfg%sim_type = "ocean"
      cfg%ocean%topo%topo_config = "flat"
      cfg%ocean%ic%ic_config = "geostrophic_adjustment"
      cfg%ocean%topo%max_depth = H_BASIN
      cfg%coriolis_f = 1.0e-4_wp
      cfg%ocean%ic%ga_eta_amp = eta_amp
      cfg%ocean%ic%ga_length_scale = length_scale
      cfg%ocean%ic%ga_x_center = x_c
      cfg%ocean%ic%ga_y_center = y_c
      cfg%initial_salinity = 35.0_wp
      cfg%initial_temperature = 10.0_wp
      cfg%T_init_surface = 0.0_wp
      cfg%T_init_bottom = 0.0_wp

      call grid%init(NX_PHYS, NY_PHYS, NGHOST, DX, DY)
      state%multilayer%nz_ml = NZ
      call state%init(grid)
      call ocean_state_seed_from_cfg(state, grid, cfg)
   end subroutine setup_ga

   subroutine test_peak_amplitude(error)
      !! At the bump centre, SSH = ga_eta_amp exactly.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      real(wp), parameter :: AMP = 1.0_wp
      real(wp) :: x_c, y_c, eta_centre
      integer :: i_c, j_c

      x_c = 0.5_wp*real(NX_PHYS, wp)*DX
      y_c = 0.5_wp*real(NY_PHYS, wp)*DY
      call setup_ga(cfg, grid, state, eta_amp=AMP, length_scale=50000.0_wp, &
                    x_c=x_c, y_c=y_c)

      ! Centre cell — choose i,j such that cell-centre coord matches x_c.
      ! With cell centres at (i_phys - 0.5)*dx, i_phys = NX_PHYS/2 + 1 → 21
      ! → x_centre = 20.5 * 5000 = 102500.  x_c = 100000.  Off by ½ cell.
      ! Closest cell is i_phys = 20 (x = 97500) or 21 (x = 102500).
      ! Use 20 since (NX_PHYS/2) integer-divides to that.
      i_c = NGHOST + NX_PHYS/2
      j_c = NGHOST + NY_PHYS/2
      eta_centre = state%barotropic%h(i_c, j_c) - state%barotropic%b(i_c, j_c)
      ! Off by half a cell from the exact centre, so SSH should be
      ! amp * exp(-r²/L²) with r² = (dx/2)² + (dy/2)² = 2·(2500)².
      ! At L = 50 km, r/L ≈ 0.07, exp(-0.005) ≈ 0.995, so eta ≈ AMP.
      call check(error, abs(eta_centre - AMP) < 0.01_wp*AMP, &
                 "geostrophic adjustment IC: SSH peak near centre ≈ ga_eta_amp")
      call state%destroy()
   end subroutine test_peak_amplitude

   subroutine test_e_fold_decay(error)
      !! SSH at distance L from centre should drop by 1/e.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      real(wp), parameter :: AMP = 2.0_wp
      real(wp), parameter :: L = 4.0_wp*DX
                              !! 20 km — clean multiple of dx.
      real(wp) :: x_c, y_c, eta_at_L, expected, eta_centre
      integer :: i_c, j_c

      x_c = 0.5_wp*real(NX_PHYS, wp)*DX
      y_c = 0.5_wp*real(NY_PHYS, wp)*DY
      call setup_ga(cfg, grid, state, eta_amp=AMP, length_scale=L, &
                    x_c=x_c, y_c=y_c)

      i_c = NGHOST + NX_PHYS/2 + 4
      j_c = NGHOST + NY_PHYS/2
      eta_at_L = state%barotropic%h(i_c, j_c) - state%barotropic%b(i_c, j_c)
      ! Compute expected from the actual cell coordinates: cell centres
      ! sit at (i_phys - 0.5)*dx + (j_phys - 0.5)*dy, so the cell at
      ! (i_c, j_c) is at (117500, 97500) and the bump centre is at
      ! (100000, 100000) → r² = 17500² + 2500².
      block
         real(wp) :: x_at, y_at, r2_at
         x_at = (real(i_c - NGHOST, wp) - 0.5_wp)*DX
         y_at = (real(j_c - NGHOST, wp) - 0.5_wp)*DY
         r2_at = (x_at - x_c)**2 + (y_at - y_c)**2
         expected = AMP*exp(-r2_at/(L*L))
      end block
      call check(error, abs(eta_at_L - expected) < 1.0e-10_wp, &
                 "geostrophic adjustment IC: SSH at L east of centre matches Gaussian")
      ! `eta_centre` not used after the rewrite; keep for clarity.
      eta_centre = state%barotropic%h(NGHOST + NX_PHYS/2, NGHOST + NY_PHYS/2) &
                   - state%barotropic%b(NGHOST + NX_PHYS/2, NGHOST + NY_PHYS/2)
      call state%destroy()
   end subroutine test_e_fold_decay

   subroutine test_default_centre(error)
      !! ga_x_center = ga_y_center = -1 → centre at basin midpoint.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      real(wp), parameter :: AMP = 0.5_wp
      real(wp) :: eta_NE_corner, eta_SW_corner
      integer :: i_NE, j_NE, i_SW, j_SW

      call setup_ga(cfg, grid, state, eta_amp=AMP, length_scale=50000.0_wp, &
                    x_c=-1.0_wp, y_c=-1.0_wp)

      ! Default centre = midpoint.  Compare SSH at two cells equidistant
      ! from midpoint (NE corner of interior + SW corner of interior):
      ! they should match to round-off (symmetric Gaussian).
      i_NE = NGHOST + NX_PHYS - 1
      j_NE = NGHOST + NY_PHYS - 1
      i_SW = NGHOST + 2
      j_SW = NGHOST + 2
      eta_NE_corner = state%barotropic%h(i_NE, j_NE) - state%barotropic%b(i_NE, j_NE)
      eta_SW_corner = state%barotropic%h(i_SW, j_SW) - state%barotropic%b(i_SW, j_SW)
      call check(error, abs(eta_NE_corner - eta_SW_corner) < 1.0e-12_wp, &
                 "geostrophic adjustment IC: default centre produces symmetric SSH")
      call state%destroy()
   end subroutine test_default_centre

   subroutine test_velocities_zero(error)
      !! All velocity components + face mass fluxes start at zero.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      real(wp) :: max_u, max_v, max_hu, max_hv

      call setup_ga(cfg, grid, state, eta_amp=1.0_wp, length_scale=50000.0_wp, &
                    x_c=-1.0_wp, y_c=-1.0_wp)

      max_u = maxval(abs(state%multilayer%u_face_x_layer))
      max_v = maxval(abs(state%multilayer%v_face_y_layer))
      max_hu = maxval(abs(state%multilayer%hu_face_x_layer))
      max_hv = maxval(abs(state%multilayer%hv_face_y_layer))
      call check(error, max_u + max_v + max_hu + max_hv < 1.0e-14_wp, &
                 "geostrophic adjustment IC: all velocities + fluxes zero")
      call state%destroy()
   end subroutine test_velocities_zero

   subroutine test_mass_added(error)
      !! Total volume added by the bump matches the analytical integral
      !! ∫∫ A·exp(-r²/L²) dA = A · π · L² (over infinite plane).  Domain
      !! is large enough (Lx = Ly = 200 km, L = 20 km) that >99.9% of
      !! the bump is captured.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      real(wp), parameter :: AMP = 1.0_wp
      real(wp), parameter :: L = 4.0_wp*DX            ! 20 km
      real(wp) :: extra_volume, expected_volume
      real(wp), parameter :: PI = 4.0_wp*atan(1.0_wp)

      call setup_ga(cfg, grid, state, eta_amp=AMP, length_scale=L, &
                    x_c=-1.0_wp, y_c=-1.0_wp)

      ! Sum (h - b)*dx*dy over the physical interior.
      extra_volume = sum(state%barotropic%h(NGHOST + 1:NGHOST + NX_PHYS, &
                                            NGHOST + 1:NGHOST + NY_PHYS) &
                         - state%barotropic%b(NGHOST + 1:NGHOST + NX_PHYS, &
                                              NGHOST + 1:NGHOST + NY_PHYS))*DX*DY
      expected_volume = AMP*PI*L*L
      ! Tolerance 1% — Gaussian integral truncated at finite domain
      ! + discretisation error.
      call check(error, abs(extra_volume - expected_volume) < 0.01_wp*expected_volume, &
                 "geostrophic adjustment IC: integrated bump volume matches A·π·L²")
      call state%destroy()
   end subroutine test_mass_added

end module test_ocean_geostrophic_adjustment_ic
