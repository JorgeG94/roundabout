!! Unit tests for the ocean curvilinear metrics slot (`rdb_ocean_metrics`).
!! Design §5 subset for M1 (no kernel consumers):
!!   * T2a cartesian: every metric uniform, area = dx*dy, Idx*dx = 1
!!     exactly, ratio bundle === 1.
!!   * T2b spherical spot checks vs hand-computed values at several
!!     staggers incl. a ghost row.
!!   * T2c sum(areaT) over the physical domain === the analytic
!!     sum_j R^2 cos(phi_j) dlambda dphi (the analytic-derivative form's
!!     EXACT identity — NOT the spherical-cap integral).
!!   * T2d Adcroft reciprocal: zero-width -> zero inverse, no NaN.
!!   * T4 planetary f at corner + centre vs 2*Omega*sin(lat) spot
!!     values, and beta_plane fill BIT-IDENTICAL to the existing
!!     coriolis_adv / EPBL fills.
!!   * lifecycle: init/destroy twice safely.
module test_ocean_metrics
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t, metrics_finalize, &
                                metrics_fill_cartesian, metrics_fill_spherical, &
                                metrics_fill_coriolis, adcroft_recip, &
                                CORIOLIS_SCHEME_BETA_PLANE, CORIOLIS_SCHEME_PLANETARY
   implicit none
   private

   public :: collect_ocean_metrics_tests

   integer, parameter :: NGHOST = 2
   real(wp), parameter :: DEG2RAD = 3.14159265358979323846_wp/180.0_wp

contains

   subroutine collect_ocean_metrics_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("metrics_cartesian_uniform", test_cartesian), &
                  new_unittest("metrics_spherical_spot_values", test_spherical_spot), &
                  new_unittest("metrics_spherical_area_identity", test_area_identity), &
                  new_unittest("metrics_adcroft_zero_width", test_adcroft), &
                  new_unittest("metrics_coriolis_planetary", test_coriolis_planetary), &
                  new_unittest("metrics_coriolis_beta_bit_identical", test_coriolis_beta), &
                  new_unittest("metrics_coriolis_fill_offset", test_coriolis_fill_offset), &
                  new_unittest("metrics_lifecycle_double", test_lifecycle) &
                  ]
   end subroutine collect_ocean_metrics_tests

   function make_grid(nx_phys, ny_phys, dx, dy) result(g)
      integer, intent(in) :: nx_phys, ny_phys
      real(wp), intent(in) :: dx, dy
      type(hgrid_t) :: g
      call g%init(nx_phys, ny_phys, NGHOST, dx, dy)
   end function make_grid

   ! ---- T2a: cartesian — uniform, area=dx*dy, Idx*dx=1, ratios===1 ----
   subroutine test_cartesian(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_metrics_t) :: m
      type(hgrid_t) :: g
      real(wp), parameter :: dx = 250.0_wp, dy = 400.0_wp
      integer :: i, j

      g = make_grid(8, 6, dx, dy)
      call m%init(g)
      call metrics_fill_cartesian(m, g, dx, dy)
      call metrics_finalize(m)

      ! Every length constant (incl. ghosts — check full arrays).
      call check(error, all(m%dxT == dx), "dxT not uniform dx"); if (allocated(error)) return
      call check(error, all(m%dyT == dy), "dyT not uniform dy"); if (allocated(error)) return
      call check(error, all(m%dxCu == dx), "dxCu not dx"); if (allocated(error)) return
      call check(error, all(m%dyCv == dy), "dyCv not dy"); if (allocated(error)) return
      call check(error, all(m%dxBu == dx), "dxBu not dx"); if (allocated(error)) return
      call check(error, all(m%dy_cu == dy), "dy_cu not dy"); if (allocated(error)) return
      call check(error, all(m%dx_cv == dx), "dx_cv not dx"); if (allocated(error)) return

      ! area = dx*dy everywhere.
      call check(error, all(abs(m%areaT - dx*dy) < 1.0e-9_wp), "areaT /= dx*dy")
      if (allocated(error)) return
      call check(error, all(abs(m%areaBu - dx*dy) < 1.0e-9_wp), "areaBu /= dx*dy")
      if (allocated(error)) return

      ! Idx*dx = 1 EXACTLY (single source).
      do j = 1, size(m%dxT, 2)
         do i = 1, size(m%dxT, 1)
            call check(error, m%idxT(i, j)*m%dxT(i, j) == 1.0_wp, "idxT*dxT /= 1")
            if (allocated(error)) return
            call check(error, m%idyT(i, j)*m%dyT(i, j) == 1.0_wp, "idyT*dyT /= 1")
            if (allocated(error)) return
         end do
      end do
      call check(error, all(m%iareaT*m%areaT == 1.0_wp), "iareaT*areaT /= 1")
      if (allocated(error)) return

      ! Ratio bundle = the geometric length ratio (dy/dx, dx/dy) — it
      ! reduces to 1 only when cells are SQUARE (dx==dy).  Here dx/=dy,
      ! so assert the exact rectangular values; dx2h == dx^2 etc.
      call check(error, all(m%dy_dxT == dy/dx), "dy_dxT /= dy/dx")
      if (allocated(error)) return
      call check(error, all(m%dx_dyT == dx/dy), "dx_dyT /= dx/dy")
      if (allocated(error)) return
      call check(error, all(m%dy_dxBu == dy/dx), "dy_dxBu /= dy/dx")
      if (allocated(error)) return
      call check(error, all(m%dx_dyBu == dx/dy), "dx_dyBu /= dx/dy")
      if (allocated(error)) return
      call check(error, all(m%dx2h == dx*dx), "dx2h /= dx^2"); if (allocated(error)) return
      call check(error, all(m%dy2q == dy*dy), "dy2q /= dy^2"); if (allocated(error)) return
      call m%destroy()

      ! On a SQUARE uniform grid the whole ratio bundle === 1 exactly
      ! (the hvisc tensor reduces to the scalar form — the M2c gate).
      g = make_grid(5, 5, dx, dx)
      call m%init(g)
      call metrics_fill_cartesian(m, g, dx, dx)
      call metrics_finalize(m)
      call check(error, all(m%dy_dxT == 1.0_wp), "square dy_dxT /= 1")
      if (allocated(error)) return
      call check(error, all(m%dx_dyT == 1.0_wp), "square dx_dyT /= 1")
      if (allocated(error)) return
      call check(error, all(m%dy_dxBu == 1.0_wp), "square dy_dxBu /= 1")
      if (allocated(error)) return
      call check(error, all(m%dx_dyBu == 1.0_wp), "square dx_dyBu /= 1")
      call m%destroy()
   end subroutine test_cartesian

   ! ---- T2b: spherical spot checks at known staggers + a ghost row ----
   subroutine test_spherical_spot(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_metrics_t) :: m
      type(hgrid_t) :: g
      real(wp), parameter :: R = 6.378e6_wp
      real(wp), parameter :: lon_w = 0.0_wp, lat_s = 0.0_wp
      real(wp), parameter :: dlon = 1.0_wp, dlat = 1.0_wp  ! degrees
      real(wp) :: expect, latv
      integer :: i, j, ng

      ! Sized so the lat-30 row falls in the interior.  T-row j has
      ! lat = lat_s + (j - ng - 0.5)*dlat.  For lat = 29.5 -> j-ng = 30.
      g = make_grid(4, 64, dlon, dlat)
      ng = NGHOST
      call m%init(g)
      call metrics_fill_spherical(m, g, lon_w, lat_s, dlon, dlat, R)
      call metrics_finalize(m)

      ! dy is latitude-independent: R * dlat_rad.
      expect = R*dlat*DEG2RAD
      call check(error, abs(m%dyT(1 + ng, 1 + ng) - expect) < 1.0e-6_wp*expect, &
                 "dyT /= R*dlat_rad"); if (allocated(error)) return

      ! dxT at the T-row whose centre latitude is 29.5 deg.
      j = ng + 30   ! lat = (30 - 0.5)*1 = 29.5
      latv = lat_s + (real(j - ng, wp) - 0.5_wp)*dlat
      expect = R*cos(latv*DEG2RAD)*dlon*DEG2RAD
      call check(error, abs(m%dxT(1 + ng, j) - expect) < 1.0e-6_wp*expect, &
                 "dxT spot (lat 29.5) wrong"); if (allocated(error)) return
      call check(error, abs(m%geolatT(1 + ng, j) - 29.5_wp) < 1.0e-9_wp, &
                 "geolatT spot wrong"); if (allocated(error)) return

      ! Corner (Bu) length uses the CORNER latitude: lat_b = (j-ng-1)*dlat.
      latv = lat_s + real(j - ng - 1, wp)*dlat   ! = 29.0
      expect = R*cos(latv*DEG2RAD)*dlon*DEG2RAD
      call check(error, abs(m%dxBu(1 + ng, j) - expect) < 1.0e-6_wp*expect, &
                 "dxBu spot (corner lat 29.0) wrong"); if (allocated(error)) return
      call check(error, abs(m%geolatBu(1 + ng, j) - 29.0_wp) < 1.0e-9_wp, &
                 "geolatBu spot wrong"); if (allocated(error)) return

      ! GHOST row j=1 (below the south edge): lat_centre = (1-ng-0.5)*dlat
      ! = (1-2-0.5) = -1.5 deg.  The formula must extend there.
      latv = lat_s + (real(1 - ng, wp) - 0.5_wp)*dlat
      expect = R*cos(latv*DEG2RAD)*dlon*DEG2RAD
      call check(error, abs(m%dxT(1, 1) - expect) < 1.0e-6_wp*expect, &
                 "dxT ghost row not filled by formula"); if (allocated(error)) return
      call check(error, m%dxT(1, 1) > 0.0_wp, "dxT ghost row is zero")
      if (allocated(error)) return

      call m%destroy()
   end subroutine test_spherical_spot

   ! ---- T2c: sum(areaT) over the physical domain === analytic sum ----
   subroutine test_area_identity(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_metrics_t) :: m
      type(hgrid_t) :: g
      real(wp), parameter :: R = 6.378e6_wp
      real(wp), parameter :: lon_w = 10.0_wp, lat_s = -20.0_wp
      real(wp), parameter :: dlon = 0.5_wp, dlat = 0.5_wp
      integer, parameter :: NXP = 12, NYP = 20
      real(wp) :: sum_area, sum_ana, latv
      integer :: i, j, ng

      ! NOTE: for the ANALYTIC-DERIVATIVE metric form the exact identity
      ! is sum dxT(j)*dyT == sum R^2 cos(phi_j) dlambda dphi (a per-row
      ! rectangle sum at the cell-centre latitude) — NOT the
      ! spherical-cap integral R^2 dlambda (sin phi_N - sin phi_S),
      ! which the analytic form only approximates.  We assert against
      ! the rectangle sum.
      g = make_grid(NXP, NYP, dlon, dlat)
      ng = NGHOST
      call m%init(g)
      call metrics_fill_spherical(m, g, lon_w, lat_s, dlon, dlat, R)
      call metrics_finalize(m)

      sum_area = 0.0_wp
      sum_ana = 0.0_wp
      do j = ng + 1, ng + NYP
         latv = lat_s + (real(j - ng, wp) - 0.5_wp)*dlat
         do i = ng + 1, ng + NXP
            sum_area = sum_area + m%areaT(i, j)
            sum_ana = sum_ana + R*R*cos(latv*DEG2RAD)*(dlon*DEG2RAD)*(dlat*DEG2RAD)
         end do
      end do

      call check(error, abs(sum_area - sum_ana) <= 1.0e-6_wp*sum_ana, &
                 "sum(areaT) /= analytic rectangle sum")
      call m%destroy()
   end subroutine test_area_identity

   ! ---- T2d: Adcroft reciprocal — zero-width -> zero inverse ----
   subroutine test_adcroft(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_metrics_t) :: m
      type(hgrid_t) :: g
      real(wp), parameter :: dx = 100.0_wp, dy = 100.0_wp

      ! Bare function: 0 -> 0, finite -> 1/x.
      call check(error, adcroft_recip(0.0_wp) == 0.0_wp, "recip(0) /= 0")
      if (allocated(error)) return
      call check(error, adcroft_recip(4.0_wp) == 0.25_wp, "recip(4) /= 0.25")
      if (allocated(error)) return

      ! Inject a zero-width face and finalize: inverse must be 0, no NaN.
      g = make_grid(4, 4, dx, dy)
      call m%init(g)
      call metrics_fill_cartesian(m, g, dx, dy)
      m%dxCu(2, 2) = 0.0_wp
      m%areaT(3, 3) = 0.0_wp
      call metrics_finalize(m)

      call check(error, m%idxCu(2, 2) == 0.0_wp, "zero-width dxCu gave nonzero idxCu")
      if (allocated(error)) return
      call check(error, m%iareaT(3, 3) == 0.0_wp, "zero areaT gave nonzero iareaT")
      if (allocated(error)) return
      ! NaN guard (x /= x is true only for NaN).
      call check(error,.not. any(m%idxCu /= m%idxCu), "NaN in idxCu")
      if (allocated(error)) return
      call check(error,.not. any(m%iareaT /= m%iareaT), "NaN in iareaT")

      call m%destroy()
   end subroutine test_adcroft

   ! ---- T4: planetary f at corner + centre vs 2*Omega*sin(lat) ----
   subroutine test_coriolis_planetary(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_metrics_t) :: m
      type(hgrid_t) :: g
      real(wp), parameter :: R = 6.378e6_wp, omega = 7.2921e-5_wp
      real(wp), parameter :: lon_w = 0.0_wp, lat_s = 0.0_wp, dlon = 1.0_wp, dlat = 1.0_wp
      real(wp), allocatable :: f_corner(:, :), f_centre(:, :)
      real(wp) :: expect, latv
      integer :: nx, ny, j, ng

      g = make_grid(4, 60, dlon, dlat)
      nx = g%nx_total; ny = g%ny_total; ng = NGHOST
      call m%init(g)
      call metrics_fill_spherical(m, g, lon_w, lat_s, dlon, dlat, R)
      call metrics_finalize(m)

      allocate (f_corner(nx + 1, ny + 1), source=0.0_wp)
      allocate (f_centre(nx, ny), source=0.0_wp)
      call metrics_fill_coriolis(m, CORIOLIS_SCHEME_PLANETARY, 0.0_wp, 0.0_wp, &
                                 0.0_wp, omega, g, f_corner, f_centre)

      ! Centre at the T-row of lat 29.5.
      j = ng + 30
      latv = lat_s + (real(j - ng, wp) - 0.5_wp)*dlat
      expect = abs(2.0_wp*omega*sin(latv*DEG2RAD))
      call check(error, abs(f_centre(1 + ng, j) - expect) < 1.0e-12_wp, &
                 "planetary f_centre /= 2*Om*sin(lat)"); if (allocated(error)) return

      ! Corner at lat 29.0.
      latv = lat_s + real(j - ng - 1, wp)*dlat
      expect = 2.0_wp*omega*sin(latv*DEG2RAD)
      call check(error, abs(f_corner(1 + ng, j) - expect) < 1.0e-12_wp, &
                 "planetary f_corner /= 2*Om*sin(lat)")

      call m%destroy()
   end subroutine test_coriolis_planetary

   ! ---- T4: beta_plane fill BIT-IDENTICAL to the existing fills ----
   subroutine test_coriolis_beta(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_metrics_t) :: m
      type(hgrid_t) :: g
      real(wp), parameter :: dx = 1000.0_wp, dy = 2000.0_wp
      real(wp), parameter :: f_0 = 7.0e-5_wp, beta = 2.0e-11_wp, y_ref = 5.0e4_wp
      real(wp), allocatable :: f_corner(:, :), f_centre(:, :)
      real(wp) :: y, expect
      integer :: nx, ny, i, j, ng

      g = make_grid(6, 8, dx, dy)
      nx = g%nx_total; ny = g%ny_total; ng = NGHOST
      call m%init(g)
      call metrics_fill_cartesian(m, g, dx, dy)
      call metrics_finalize(m)

      allocate (f_corner(nx + 1, ny + 1), source=0.0_wp)
      allocate (f_centre(nx, ny), source=0.0_wp)
      call metrics_fill_coriolis(m, CORIOLIS_SCHEME_BETA_PLANE, f_0, beta, &
                                 y_ref, 0.0_wp, g, f_corner, f_centre)

      ! Corner: reproduce coriolis_adv_set_beta_plane EXACTLY.
      !   y = (j - 1 - nghost)*dy ; f = f_0 + beta*(y - y_ref).
      do j = 1, ny + 1
         y = real(j - 1 - ng, wp)*dy
         do i = 1, nx + 1
            expect = f_0 + beta*(y - y_ref)
            call check(error, f_corner(i, j) == expect, &
                       "beta f_corner not bit-identical to coriolis_adv")
            if (allocated(error)) return
         end do
      end do

      ! Centre: reproduce EPBL / kappa-shear set_f_centre EXACTLY.
      !   y = (j - nghost - 0.5)*dy ; f = abs(f_0 + beta*(y - y_ref)).
      do j = 1, ny
         y = (real(j - ng, wp) - 0.5_wp)*dy
         do i = 1, nx
            expect = abs(f_0 + beta*(y - y_ref))
            call check(error, f_centre(i, j) == expect, &
                       "beta f_centre not bit-identical to set_f_centre")
            if (allocated(error)) return
         end do
      end do

      call m%destroy()
   end subroutine test_coriolis_beta

   ! ---- O4 offset: subdomain fill == matching window of the global fill ----
   subroutine test_coriolis_fill_offset(error)
      !! Verify that `metrics_fill_coriolis` (beta_plane) and
      !! `metrics_fill_spherical` produce, for a half-domain grid with a
      !! non-zero `j_offset_global`, EXACTLY the same values as the
      !! corresponding physical rows of the full global grid.
      !! This is the regression gate for the O4 py-split Coriolis bug:
      !! before the fix every rank restarted y from 0.
      type(error_type), allocatable, intent(out) :: error
      ! ---- shared grid params ----
      integer, parameter :: NXP_FULL = 8, NYP_FULL = 12
      integer, parameter :: NYP_HALF = 6        ! northern half: rows 7..12
      integer, parameter :: NY_OFFSET = 6       ! j_offset_global for the half grid
      real(wp), parameter :: DX = 1000.0_wp, DY = 1000.0_wp
      real(wp), parameter :: F_0 = 1.0e-4_wp, BETA = 1.0e-11_wp, Y_REF = 3000.0_wp
      real(wp), parameter :: R = 6.378e6_wp
      real(wp), parameter :: LON_W = 10.0_wp, LAT_S = -20.0_wp
      real(wp), parameter :: DLON = 1.0_wp, DLAT = 1.0_wp
      ! beta-plane leg objects
      type(ocean_metrics_t) :: mg, mh
      type(hgrid_t) :: gg, gh
      real(wp), allocatable :: fc_g(:, :), fctr_g(:, :)
      real(wp), allocatable :: fc_h(:, :), fctr_h(:, :)
      ! spherical leg objects
      type(ocean_metrics_t) :: msg, msh
      type(hgrid_t) :: gsg, gsh
      ! scalars
      integer :: nx_g, ny_g, nx_h, ny_h, ng, i, j

      ng = NGHOST

      ! ---- beta-plane leg ----
      ! Full global grid (offset stays 0 — default).
      call gg%init(NXP_FULL, NYP_FULL, ng, DX, DY)
      nx_g = gg%nx_total; ny_g = gg%ny_total
      call mg%init(gg)
      call metrics_fill_cartesian(mg, gg, DX, DY)
      call metrics_finalize(mg)
      allocate (fc_g(nx_g + 1, ny_g + 1), source=0.0_wp)
      allocate (fctr_g(nx_g, ny_g), source=0.0_wp)
      call metrics_fill_coriolis(mg, CORIOLIS_SCHEME_BETA_PLANE, F_0, BETA, &
                                 Y_REF, 0.0_wp, gg, fc_g, fctr_g)

      ! Northern half: local ny = NYP_HALF, j_offset = NY_OFFSET.
      call gh%init(NXP_FULL, NYP_HALF, ng, DX, DY)
      gh%j_offset_global = NY_OFFSET
      nx_h = gh%nx_total; ny_h = gh%ny_total
      call mh%init(gh)
      call metrics_fill_cartesian(mh, gh, DX, DY)
      call metrics_finalize(mh)
      allocate (fc_h(nx_h + 1, ny_h + 1), source=0.0_wp)
      allocate (fctr_h(nx_h, ny_h), source=0.0_wp)
      call metrics_fill_coriolis(mh, CORIOLIS_SCHEME_BETA_PLANE, F_0, BETA, &
                                 Y_REF, 0.0_wp, gh, fc_h, fctr_h)

      ! Assert corner physical rows of the half-grid == global rows shifted
      ! by NY_OFFSET.  Physical corner j range: 1..ny_h+1 maps globally to
      ! 1+NY_OFFSET..ny_h+1+NY_OFFSET.
      do j = 1, ny_h + 1
         do i = 1, nx_h + 1
            call check(error, fc_h(i, j) == fc_g(i, j + NY_OFFSET), &
                       "beta f_corner offset mismatch at j="//achar(48 + j))
            if (allocated(error)) return
         end do
      end do

      ! Assert centre physical rows.
      do j = 1, ny_h
         do i = 1, nx_h
            call check(error, fctr_h(i, j) == fctr_g(i, j + NY_OFFSET), &
                       "beta f_centre offset mismatch at j="//achar(48 + j))
            if (allocated(error)) return
         end do
      end do

      call mg%destroy()
      call mh%destroy()

      ! ---- spherical leg ----
      ! Reuse NXP_FULL x NYP_FULL / NYP_HALF layout for geolatT + geolonT.
      call gsg%init(NXP_FULL, NYP_FULL, ng, DLON, DLAT)
      call msg%init(gsg)
      call metrics_fill_spherical(msg, gsg, LON_W, LAT_S, DLON, DLAT, R)
      call metrics_finalize(msg)

      call gsh%init(NXP_FULL, NYP_HALF, ng, DLON, DLAT)
      gsh%j_offset_global = NY_OFFSET
      gsh%i_offset_global = 0
      call msh%init(gsh)
      call metrics_fill_spherical(msh, gsh, LON_W, LAT_S, DLON, DLAT, R)
      call metrics_finalize(msh)

      ! Physical T cells: local j maps to global j + NY_OFFSET.
      do j = ng + 1, ng + NYP_HALF
         do i = ng + 1, ng + NXP_FULL
            call check(error, msh%geolatT(i, j) == msg%geolatT(i, j + NY_OFFSET), &
                       "spherical geolatT offset mismatch"); if (allocated(error)) return
            call check(error, msh%geolonT(i, j) == msg%geolonT(i, j + NY_OFFSET), &
                       "spherical geolonT offset mismatch"); if (allocated(error)) return
         end do
      end do

      call msg%destroy()
      call msh%destroy()
   end subroutine test_coriolis_fill_offset

   ! ---- lifecycle: init/destroy twice safely ----
   subroutine test_lifecycle(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_metrics_t) :: m
      type(hgrid_t) :: g
      real(wp), parameter :: dx = 50.0_wp, dy = 60.0_wp

      g = make_grid(4, 4, dx, dy)
      call m%init(g)
      call check(error, m%is_init, "is_init false after first init")
      if (allocated(error)) return
      call m%destroy()
      call check(error,.not. m%is_init, "is_init true after destroy")
      if (allocated(error)) return

      call m%init(g)
      call metrics_fill_cartesian(m, g, dx, dy)
      call metrics_finalize(m)
      call check(error, m%is_init, "is_init false after second init")
      if (allocated(error)) return
      call check(error, all(m%dxT == dx), "second init not filled")
      call m%destroy()
   end subroutine test_lifecycle

end module test_ocean_metrics
