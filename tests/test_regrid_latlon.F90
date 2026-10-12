!! Unit tests for `rdb_regrid_latlon` (OM3 wave-2 PR-C2b): in-model
!! bilinear and first-order conservative regridding weights from a
!! regular lat-lon source grid onto the model T-grid, and the single ELL
!! apply kernel. Synthetic sources only, no NetCDF.
!!
!! Analytic source function: f(lon,lat) = sin(lat) cos(2 lon) +
!! 0.3 cos(5 lat) (lon, lat in radians inside the trig calls).
!!
!! Coverage:
!!   regrid_bilinear_constant / regrid_conservative_constant  -- test 1
!!   regrid_bilinear_linear_exact                             -- test 2
!!   regrid_bilinear_convergence                              -- test 3
!!   regrid_bilinear_pole_clamp                                -- test 4
!!   regrid_conservative_global_integral                       -- test 5
!!   regrid_conservative_nmax_cap                               -- test 6
!!   regrid_tile_bitid_serial                                   -- test 7
!!   regrid_apply_device                                        -- test 8
module test_regrid_latlon
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, DEG2RAD
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_spherical_metrics, make_tripolar_metrics, &
                                 destroy_cartesian_metrics
   use rdb_regrid_latlon, only: regrid_weights_t, &
                                regrid_bnds_to_edges, regrid_weights_bilinear, &
                                regrid_weights_conservative, regrid_apply_2d_impl
   implicit none
   private

   public :: collect_regrid_latlon_tests

contains

   subroutine collect_regrid_latlon_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("regrid_bilinear_constant", test_bilinear_constant), &
                  new_unittest("regrid_conservative_constant", test_conservative_constant), &
                  new_unittest("regrid_bilinear_linear_exact", test_bilinear_linear_exact), &
                  new_unittest("regrid_bilinear_convergence", test_bilinear_convergence), &
                  new_unittest("regrid_bilinear_pole_clamp", test_bilinear_pole_clamp), &
                  new_unittest("regrid_conservative_global_integral", &
                               test_conservative_global_integral), &
                  new_unittest("regrid_conservative_nmax_cap", test_conservative_nmax_cap), &
                  new_unittest("regrid_tile_bitid_serial", test_tile_bitid_serial), &
                  new_unittest("regrid_apply_device", test_apply_device) &
                  ]
   end subroutine collect_regrid_latlon_tests

   ! =====================================================================
   ! Shared helpers
   ! =====================================================================

   pure function f_analytic(lon_deg, lat_deg) result(val)
      real(wp), intent(in) :: lon_deg, lat_deg
      real(wp) :: val
      val = sin(lat_deg*DEG2RAD)*cos(2.0_wp*lon_deg*DEG2RAD) + &
            0.3_wp*cos(5.0_wp*lat_deg*DEG2RAD)
   end function f_analytic

   pure subroutine centers_from_span(n, lo, hi, c)
      integer, intent(in) :: n
      real(wp), intent(in) :: lo, hi
      real(wp), intent(out) :: c(n)
      real(wp) :: dx
      integer :: k
      dx = (hi - lo)/real(n, wp)
      do k = 1, n
         c(k) = lo + (real(k, wp) - 0.5_wp)*dx
      end do
   end subroutine centers_from_span

   pure subroutine edges_from_span(n, lo, hi, e)
      integer, intent(in) :: n
      real(wp), intent(in) :: lo, hi
      real(wp), intent(out) :: e(n + 1)
      real(wp) :: dx
      integer :: k
      dx = (hi - lo)/real(n, wp)
      do k = 1, n + 1
         e(k) = lo + real(k - 1, wp)*dx
      end do
   end subroutine edges_from_span

   ! =====================================================================
   ! Test 1: constant source -> every physical target == 1, padding OK.
   ! =====================================================================

   subroutine test_bilinear_constant(error)
      type(error_type), allocatable, intent(out) :: error
      call check_constant_one(error, .false., -10.0_wp, -20.0_wp, 1.0_wp, 1.0_wp, 40, 30)
      if (allocated(error)) return
      call check_constant_one(error, .true., 0.0_wp, -90.0_wp, 5.0_wp, 4.0_wp, 72, 40)
   end subroutine test_bilinear_constant

   subroutine test_conservative_constant(error)
      type(error_type), allocatable, intent(out) :: error
      call check_constant_one(error, .false., -10.0_wp, -20.0_wp, 1.0_wp, 1.0_wp, 40, 30, &
                              conservative=.true.)
      if (allocated(error)) return
      call check_constant_one(error, .true., 0.0_wp, -90.0_wp, 5.0_wp, 4.0_wp, 72, 40, &
                              conservative=.true.)
   end subroutine test_conservative_constant

   subroutine check_constant_one(error, tripolar, lon_w, lat_s, dlon, dlat, nxp, nyp, &
                                 conservative)
      type(error_type), allocatable, intent(out) :: error
      logical, intent(in) :: tripolar
      real(wp), intent(in) :: lon_w, lat_s, dlon, dlat
      integer, intent(in) :: nxp, nyp
      logical, intent(in), optional :: conservative
      logical :: do_cons

      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(regrid_weights_t) :: wts
      integer, parameter :: NLON = 20, NLAT = 15
      real(wp) :: lon_c(NLON), lat_c(NLAT), lon_e(NLON + 1), lat_e(NLAT + 1)
      real(wp), allocatable :: src(:, :), dst(:, :)
      integer :: ierr, i, j, k

      do_cons = .false.
      if (present(conservative)) do_cons = conservative

      call grid%init(nxp, nyp, 2, 1000.0_wp, 1000.0_wp)
      if (tripolar) then
         call make_tripolar_metrics(metrics, grid, lon_w, lat_s, dlon, dlat, 6.378e6_wp, &
                                    65.0_wp, 0.0_wp)
      else
         call make_spherical_metrics(metrics, grid, lon_w, lat_s, dlon, dlat, 6.378e6_wp)
      end if

      ! Gaussian-like source axis (uniform in sin(lat)) -- exercises the
      ! non-uniform-latitude bisection path even on a trivial constant
      ! field.
      call centers_from_span(NLON, 0.0_wp, 360.0_wp, lon_c)
      call gaussian_like_centers(NLAT, lat_c)

      if (do_cons) then
         call edges_from_span(NLON, 0.0_wp, 360.0_wp, lon_e)
         call gaussian_like_edges(NLAT, lat_e)
         call regrid_weights_conservative(wts, NLON, NLAT, lon_e, lat_e, &
                                          grid%nx_phys, grid%ny_phys, grid%nghost, &
                                          grid%nx_total, grid%ny_total, &
                                          metrics%geolonBu, metrics%geolatBu, &
                                          metrics%geolonT, metrics%geolatT, ierr)
      else
         call regrid_weights_bilinear(wts, NLON, NLAT, lon_c, lat_c, &
                                      grid%nx_phys, grid%ny_phys, grid%nghost, &
                                      grid%nx_total, grid%ny_total, &
                                      metrics%geolonT, metrics%geolatT, ierr)
      end if
      call check(error, ierr, 0)
      if (allocated(error)) return

      allocate (src(wts%nsx, wts%nsy), source=1.0_wp)
      allocate (dst(wts%nx, wts%ny), source=0.0_wp)
      call regrid_apply_2d_impl(wts%nx, wts%ny, wts%nmax, wts%nsx, wts%nsy, &
                                wts%idx_i, wts%idx_j, wts%w, src, dst)

      do j = 1, wts%ny
         do i = 1, wts%nx
            call check(error, abs(dst(i, j) - 1.0_wp) < 1.0e-14_wp)
            if (allocated(error)) return
            do k = 1, wts%nmax
               call check(error, wts%idx_i(i, j, k) >= 1 .and. wts%idx_i(i, j, k) <= wts%nsx)
               if (allocated(error)) return
               call check(error, wts%idx_j(i, j, k) >= 1 .and. wts%idx_j(i, j, k) <= wts%nsy)
               if (allocated(error)) return
            end do
         end do
      end do

      call destroy_cartesian_metrics(metrics)
   end subroutine check_constant_one

   pure subroutine gaussian_like_centers(n, c)
      !! Uniform-in-sin(lat) centres -- a coarse Gaussian-grid proxy.
      integer, intent(in) :: n
      real(wp), intent(out) :: c(n)
      real(wp) :: mu_c
      real(wp) :: mu_e(n + 1)
      integer :: k
      call edges_from_span(n, -1.0_wp, 1.0_wp, mu_e)
      do k = 1, n
         mu_c = 0.5_wp*(mu_e(k) + mu_e(k + 1))
         c(k) = asin(mu_c)/DEG2RAD
      end do
   end subroutine gaussian_like_centers

   pure subroutine gaussian_like_edges(n, e)
      integer, intent(in) :: n
      real(wp), intent(out) :: e(n + 1)
      real(wp) :: mu_e(n + 1)
      integer :: k
      call edges_from_span(n, -1.0_wp, 1.0_wp, mu_e)
      do k = 1, n + 1
         e(k) = asin(mu_e(k))/DEG2RAD
      end do
   end subroutine gaussian_like_edges

   ! =====================================================================
   ! Test 2: bilinear reproduces a linear function exactly.
   ! =====================================================================

   subroutine test_bilinear_linear_exact(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(regrid_weights_t) :: wts
      integer, parameter :: NLON = 40, NLAT = 30, NXP = 20, NYP = 15
      real(wp), parameter :: A = 1.0_wp, B = 0.02_wp, C = 0.05_wp
      real(wp) :: lon_c(NLON), lat_c(NLAT)
      real(wp), allocatable :: src(:, :), dst(:, :)
      integer :: ierr, i, j, jg
      real(wp) :: lon_t, lat_t, expect

      call grid%init(NXP, NYP, 2, 1000.0_wp, 1000.0_wp)
      call make_spherical_metrics(metrics, grid, 100.0_wp, 10.0_wp, 1.0_wp, 1.0_wp, 6.378e6_wp)

      call centers_from_span(NLON, 60.0_wp, 260.0_wp, lon_c)
      call centers_from_span(NLAT, -10.0_wp, 50.0_wp, lat_c)

      call regrid_weights_bilinear(wts, NLON, NLAT, lon_c, lat_c, &
                                   grid%nx_phys, grid%ny_phys, grid%nghost, &
                                   grid%nx_total, grid%ny_total, &
                                   metrics%geolonT, metrics%geolatT, ierr)
      call check(error, ierr, 0)
      if (allocated(error)) return

      allocate (src(wts%nsx, wts%nsy))
      do j = 1, wts%nsy
         jg = j + wts%j_lo - 1
         do i = 1, wts%nsx
            src(i, j) = A + B*lon_c(i) + C*lat_c(jg)
         end do
      end do

      allocate (dst(wts%nx, wts%ny), source=0.0_wp)
      call regrid_apply_2d_impl(wts%nx, wts%ny, wts%nmax, wts%nsx, wts%nsy, &
                                wts%idx_i, wts%idx_j, wts%w, src, dst)

      do j = 1, wts%ny
         lat_t = metrics%geolatT(grid%nghost + 1, grid%nghost + j)
         do i = 1, wts%nx
            lon_t = metrics%geolonT(grid%nghost + i, grid%nghost + j)
            expect = A + B*lon_t + C*lat_t
            call check(error, abs(dst(i, j) - expect) < 1.0e-12_wp)
            if (allocated(error)) return
         end do
      end do

      call destroy_cartesian_metrics(metrics)
   end subroutine test_bilinear_linear_exact

   ! =====================================================================
   ! Test 3: O(h^2) convergence, dateline-crossing target.
   ! =====================================================================

   subroutine test_bilinear_convergence(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      integer, parameter :: NXP = 50, NYP = 30
      integer, parameter :: NRES = 3
      integer :: nlon_r(NRES), nlat_r(NRES)
      real(wp) :: errs(NRES)
      integer :: r

      nlon_r = [36, 72, 144]
      nlat_r = nlon_r/2

      call grid%init(NXP, NYP, 2, 1000.0_wp, 1000.0_wp)
      call make_spherical_metrics(metrics, grid, -5.0_wp, -20.0_wp, 1.0_wp, 1.0_wp, 6.378e6_wp)

      do r = 1, NRES
         call conv_one_resolution(grid, metrics, nlon_r(r), nlat_r(r), errs(r))
      end do

      call check(error, errs(1)/errs(2) >= 3.6_wp)
      if (allocated(error)) then
         call destroy_cartesian_metrics(metrics)
         return
      end if
      call check(error, errs(2)/errs(3) >= 3.6_wp)

      call destroy_cartesian_metrics(metrics)
   end subroutine test_bilinear_convergence

   subroutine conv_one_resolution(grid, metrics, nlon, nlat, err_max)
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      integer, intent(in) :: nlon, nlat
      real(wp), intent(out) :: err_max
      type(regrid_weights_t) :: wts
      real(wp), allocatable :: lon_c(:), lat_c(:), src(:, :), dst(:, :)
      integer :: ierr, i, j, jg
      real(wp) :: lon_t, lat_t

      allocate (lon_c(nlon), lat_c(nlat))
      call centers_from_span(nlon, 0.0_wp, 360.0_wp, lon_c)
      call centers_from_span(nlat, -90.0_wp, 90.0_wp, lat_c)

      call regrid_weights_bilinear(wts, nlon, nlat, lon_c, lat_c, &
                                   grid%nx_phys, grid%ny_phys, grid%nghost, &
                                   grid%nx_total, grid%ny_total, &
                                   metrics%geolonT, metrics%geolatT, ierr)

      allocate (src(wts%nsx, wts%nsy))
      do j = 1, wts%nsy
         jg = j + wts%j_lo - 1
         do i = 1, wts%nsx
            src(i, j) = f_analytic(lon_c(i), lat_c(jg))
         end do
      end do

      allocate (dst(wts%nx, wts%ny), source=0.0_wp)
      call regrid_apply_2d_impl(wts%nx, wts%ny, wts%nmax, wts%nsx, wts%nsy, &
                                wts%idx_i, wts%idx_j, wts%w, src, dst)

      err_max = 0.0_wp
      do j = 1, wts%ny
         lat_t = metrics%geolatT(grid%nghost + 1, grid%nghost + j)
         do i = 1, wts%nx
            lon_t = metrics%geolonT(grid%nghost + i, grid%nghost + j)
            err_max = max(err_max, abs(dst(i, j) - f_analytic(lon_t, lat_t)))
         end do
      end do
   end subroutine conv_one_resolution

   ! =====================================================================
   ! Test 4: poleward of the last source centre clamps to the edge row.
   ! =====================================================================

   subroutine test_bilinear_pole_clamp(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(regrid_weights_t) :: wts
      integer, parameter :: NLON = 36, NLAT = 20, NXP = 4, NYP = 1
      real(wp) :: lon_c(NLON), lat_c(NLAT)
      real(wp), allocatable :: src(:, :), dst(:, :)
      integer :: ierr, i, j, k
      real(wp) :: lon_t

      call grid%init(NXP, NYP, 2, 1000.0_wp, 1000.0_wp)
      ! One physical row centred exactly at 89.9 degrees.
      call make_spherical_metrics(metrics, grid, 0.0_wp, 89.4_wp, 10.0_wp, 1.0_wp, 6.378e6_wp)

      call centers_from_span(NLON, 0.0_wp, 360.0_wp, lon_c)
      do k = 1, NLAT
         lat_c(k) = -89.57_wp + real(k - 1, wp)*(2.0_wp*89.57_wp)/real(NLAT - 1, wp)
      end do

      call regrid_weights_bilinear(wts, NLON, NLAT, lon_c, lat_c, &
                                   grid%nx_phys, grid%ny_phys, grid%nghost, &
                                   grid%nx_total, grid%ny_total, &
                                   metrics%geolonT, metrics%geolatT, ierr)
      call check(error, ierr, 0)
      if (allocated(error)) return

      do j = 1, wts%ny
         do i = 1, wts%nx
            do k = 1, wts%nmax
               call check(error, wts%idx_j(i, j, k) >= 1 .and. wts%idx_j(i, j, k) <= wts%nsy)
               if (allocated(error)) return
            end do
         end do
      end do

      ! Source field is pure longitude (linear in lon, lat-independent):
      ! clamp-to-edge-row must reproduce it exactly.
      allocate (src(wts%nsx, wts%nsy))
      do j = 1, wts%nsy
         do i = 1, wts%nsx
            src(i, j) = lon_c(i)
         end do
      end do
      allocate (dst(wts%nx, wts%ny), source=0.0_wp)
      call regrid_apply_2d_impl(wts%nx, wts%ny, wts%nmax, wts%nsx, wts%nsy, &
                                wts%idx_i, wts%idx_j, wts%w, src, dst)

      do j = 1, wts%ny
         do i = 1, wts%nx
            lon_t = metrics%geolonT(grid%nghost + i, grid%nghost + j)
            call check(error, abs(dst(i, j) - lon_t) < 1.0e-9_wp)
            if (allocated(error)) return
         end do
      end do

      call destroy_cartesian_metrics(metrics)
   end subroutine test_bilinear_pole_clamp

   ! =====================================================================
   ! Test 5: conservative global integral + coverage (spherical and
   ! tripolar global targets).
   ! =====================================================================

   subroutine test_conservative_global_integral(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: err36, err90

      call check_global_case(error, .false., 36, 18, err36)
      if (allocated(error)) return
      call check_global_case(error, .true., 36, 18, err90)
      if (allocated(error)) return
      ! O(h) sanity on the spherical target: a finer source halves the
      ! clipped-fragment discretisation error.
      call check_global_case(error, .false., 90, 45, err90)
      if (allocated(error)) return
      call check(error, err36/max(err90, 1.0e-30_wp) >= 1.8_wp)
   end subroutine test_conservative_global_integral

   subroutine check_global_case(error, tripolar, nlon, nlat, err_out)
      type(error_type), allocatable, intent(out) :: error
      logical, intent(in) :: tripolar
      integer, intent(in) :: nlon, nlat
      real(wp), intent(out) :: err_out

      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(regrid_weights_t) :: wts
      integer, parameter :: NXP_S = 72, NYP_S = 36
      integer, parameter :: NXP_T = 72, NYP_T = 40
      real(wp), allocatable :: lon_c(:), lat_c(:), lon_e(:), lat_e(:)
      real(wp), allocatable :: lon_bnds(:, :), lat_bnds(:, :)
      real(wp), allocatable :: src(:, :), dst(:, :)
      real(wp) :: dlon, dlat, sum_src_area, sum_tgt_area, sum_src_f, sum_tgt_f
      real(wp) :: total_area_t, lon_t, lat_t
      integer :: ierr, i, j, jg, nxp, nyp

      allocate (lon_c(nlon), lat_c(nlat), lon_e(nlon + 1), lat_e(nlat + 1))
      allocate (lon_bnds(2, nlon), lat_bnds(2, nlat))
      call centers_from_span(nlon, 0.0_wp, 360.0_wp, lon_c)
      call centers_from_span(nlat, -90.0_wp, 90.0_wp, lat_c)
      dlon = 360.0_wp/real(nlon, wp)
      dlat = 180.0_wp/real(nlat, wp)
      do i = 1, nlon
         lon_bnds(1, i) = real(i - 1, wp)*dlon
         lon_bnds(2, i) = real(i, wp)*dlon
      end do
      do j = 1, nlat
         lat_bnds(1, j) = -90.0_wp + real(j - 1, wp)*dlat
         lat_bnds(2, j) = -90.0_wp + real(j, wp)*dlat
      end do
      call regrid_bnds_to_edges(nlon, lon_bnds, lon_e)
      call regrid_bnds_to_edges(nlat, lat_bnds, lat_e)

      if (tripolar) then
         nxp = NXP_T
         nyp = NYP_T
      else
         nxp = NXP_S
         nyp = NYP_S
      end if
      call grid%init(nxp, nyp, 2, 1000.0_wp, 1000.0_wp)
      if (tripolar) then
         call make_tripolar_metrics(metrics, grid, 0.0_wp, -90.0_wp, 5.0_wp, 4.0_wp, &
                                    6.378e6_wp, 65.0_wp, 0.0_wp)
      else
         call make_spherical_metrics(metrics, grid, -180.0_wp, -90.0_wp, 5.0_wp, 5.0_wp, &
                                     6.378e6_wp)
      end if

      call regrid_weights_conservative(wts, nlon, nlat, lon_e, lat_e, &
                                       grid%nx_phys, grid%ny_phys, grid%nghost, &
                                       grid%nx_total, grid%ny_total, &
                                       metrics%geolonBu, metrics%geolatBu, &
                                       metrics%geolonT, metrics%geolatT, ierr)
      call check(error, ierr, 0)
      if (allocated(error)) then
         call destroy_cartesian_metrics(metrics)
         return
      end if
      call check(error, wts%n_fallback, 0)
      if (allocated(error)) then
         call destroy_cartesian_metrics(metrics)
         return
      end if

      ! Coverage: sum of target areas == sum of source areas == 4*pi
      ! (both in (rad, sin(phi)) units, i.e. without the R^2 factor).
      total_area_t = 0.0_wp
      do j = 1, wts%ny
         do i = 1, wts%nx
            total_area_t = total_area_t + wts%area_t(i, j)
         end do
      end do
      sum_src_area = 0.0_wp
      do j = 1, nlat
         sum_src_area = sum_src_area + (lon_e(nlon + 1) - lon_e(1))*DEG2RAD* &
                        (sin(lat_e(j + 1)*DEG2RAD) - sin(lat_e(j)*DEG2RAD))
      end do
      call check(error, abs(total_area_t - sum_src_area) <= 1.0e-12_wp*sum_src_area)
      if (allocated(error)) then
         call destroy_cartesian_metrics(metrics)
         return
      end if
      call check(error, abs(sum_src_area - 4.0_wp*acos(-1.0_wp)) <= 1.0e-12_wp)
      if (allocated(error)) then
         call destroy_cartesian_metrics(metrics)
         return
      end if

      ! Global integral: sum(area_t * dst) == sum(A_s * f_s) to a tight
      ! relative tolerance.
      allocate (src(wts%nsx, wts%nsy))
      sum_src_f = 0.0_wp
      do j = 1, wts%nsy
         jg = j + wts%j_lo - 1
         do i = 1, wts%nsx
            src(i, j) = f_analytic(lon_c(i), lat_c(jg))
            sum_src_f = sum_src_f + src(i, j)*dlon*DEG2RAD* &
                        (sin(lat_e(jg + 1)*DEG2RAD) - sin(lat_e(jg)*DEG2RAD))
         end do
      end do

      allocate (dst(wts%nx, wts%ny), source=0.0_wp)
      call regrid_apply_2d_impl(wts%nx, wts%ny, wts%nmax, wts%nsx, wts%nsy, &
                                wts%idx_i, wts%idx_j, wts%w, src, dst)

      sum_tgt_f = 0.0_wp
      err_out = 0.0_wp
      do j = 1, wts%ny
         lat_t = metrics%geolatT(grid%nghost + 1, grid%nghost + j)
         do i = 1, wts%nx
            sum_tgt_f = sum_tgt_f + wts%area_t(i, j)*dst(i, j)
            lon_t = metrics%geolonT(grid%nghost + i, grid%nghost + j)
            err_out = max(err_out, abs(dst(i, j) - f_analytic(lon_t, lat_t)))
         end do
      end do
      call check(error, abs(sum_tgt_f - sum_src_f) <= 1.0e-12_wp*max(1.0_wp, abs(sum_src_f)))

      call destroy_cartesian_metrics(metrics)
   end subroutine check_global_case

   ! =====================================================================
   ! Test 6: nmax cap fail-loud path.
   ! =====================================================================

   subroutine test_conservative_nmax_cap(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(regrid_weights_t) :: wts
      integer, parameter :: NLON = 64, NLAT = 32, NXP = 2, NYP = 2
      real(wp) :: lon_e(NLON + 1), lat_e(NLAT + 1)
      integer :: ierr

      call grid%init(NXP, NYP, 2, 1000.0_wp, 1000.0_wp)
      call make_spherical_metrics(metrics, grid, 0.0_wp, 0.0_wp, 4.0_wp, 4.0_wp, 6.378e6_wp)

      call edges_from_span(NLON, 0.0_wp, 36.0_wp, lon_e)
      call edges_from_span(NLAT, 0.0_wp, 18.0_wp, lat_e)

      call regrid_weights_conservative(wts, NLON, NLAT, lon_e, lat_e, &
                                       grid%nx_phys, grid%ny_phys, grid%nghost, &
                                       grid%nx_total, grid%ny_total, &
                                       metrics%geolonBu, metrics%geolatBu, &
                                       metrics%geolonT, metrics%geolatT, ierr, &
                                       nmax_cap_override=10)
      call check(error, ierr /= 0)

      call destroy_cartesian_metrics(metrics)
   end subroutine test_conservative_nmax_cap

   ! =====================================================================
   ! Test 7: one tile vs 2x2 tiles, bit-identical, on the tripolar grid.
   ! =====================================================================

   subroutine test_tile_bitid_serial(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NLON = 40, NLAT = 20
      integer, parameter :: NXG = 72, NYG = 40
      real(wp) :: lon_c(NLON), lat_c(NLAT), lon_e(NLON + 1), lat_e(NLAT + 1)
      real(wp), allocatable :: src(:, :)
      type(hgrid_t) :: grid_g
      type(ocean_metrics_t) :: metrics_g
      type(regrid_weights_t) :: wts_g
      real(wp), allocatable :: dst_g_bl(:, :), dst_g_cs(:, :)
      integer :: ierr, i, j, jg

      call centers_from_span(NLON, 0.0_wp, 360.0_wp, lon_c)
      call centers_from_span(NLAT, -90.0_wp, 90.0_wp, lat_c)
      call edges_from_span(NLON, 0.0_wp, 360.0_wp, lon_e)
      call edges_from_span(NLAT, -90.0_wp, 90.0_wp, lat_e)

      call grid_g%init(NXG, NYG, 2, 1000.0_wp, 1000.0_wp)
      call make_tripolar_metrics(metrics_g, grid_g, 0.0_wp, -90.0_wp, 5.0_wp, 4.0_wp, &
                                 6.378e6_wp, 65.0_wp, 0.0_wp)

      call regrid_weights_bilinear(wts_g, NLON, NLAT, lon_c, lat_c, &
                                   grid_g%nx_phys, grid_g%ny_phys, grid_g%nghost, &
                                   grid_g%nx_total, grid_g%ny_total, &
                                   metrics_g%geolonT, metrics_g%geolatT, ierr)
      call check(error, ierr, 0)
      if (allocated(error)) then
         call destroy_cartesian_metrics(metrics_g)
         return
      end if

      allocate (src(NLON, NLAT))
      do j = 1, NLAT
         do i = 1, NLON
            src(i, j) = f_analytic(lon_c(i), lat_c(j))
         end do
      end do

      allocate (dst_g_bl(wts_g%nx, wts_g%ny), source=0.0_wp)
      call apply_global_window(wts_g, src, NLON, NLAT, dst_g_bl)

      call regrid_weights_conservative(wts_g, NLON, NLAT, lon_e, lat_e, &
                                       grid_g%nx_phys, grid_g%ny_phys, grid_g%nghost, &
                                       grid_g%nx_total, grid_g%ny_total, &
                                       metrics_g%geolonBu, metrics_g%geolatBu, &
                                       metrics_g%geolonT, metrics_g%geolatT, ierr)
      call check(error, ierr, 0)
      if (allocated(error)) then
         call destroy_cartesian_metrics(metrics_g)
         return
      end if
      allocate (dst_g_cs(wts_g%nx, wts_g%ny), source=0.0_wp)
      call apply_global_window(wts_g, src, NLON, NLAT, dst_g_cs)

      call destroy_cartesian_metrics(metrics_g)

      call check_tiles_match(error, NLON, NLAT, lon_c, lat_c, lon_e, lat_e, src, &
                             dst_g_bl, dst_g_cs)
   end subroutine test_tile_bitid_serial

   subroutine apply_global_window(wts, src_full, nlon, nlat, dst)
      type(regrid_weights_t), intent(in) :: wts
      integer, intent(in) :: nlon, nlat
      real(wp), intent(in) :: src_full(nlon, nlat)
      real(wp), intent(inout) :: dst(wts%nx, wts%ny)
      real(wp) :: src_win(wts%nsx, wts%nsy)
      integer :: j
      do j = 1, wts%nsy
         src_win(:, j) = src_full(:, j + wts%j_lo - 1)
      end do
      call regrid_apply_2d_impl(wts%nx, wts%ny, wts%nmax, wts%nsx, wts%nsy, &
                                wts%idx_i, wts%idx_j, wts%w, src_win, dst)
   end subroutine apply_global_window

   subroutine check_tiles_match(error, nlon, nlat, lon_c, lat_c, lon_e, lat_e, src, &
                                dst_g_bl, dst_g_cs)
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: nlon, nlat
      real(wp), intent(in) :: lon_c(nlon), lat_c(nlat), lon_e(nlon + 1), lat_e(nlat + 1)
      real(wp), intent(in) :: src(nlon, nlat)
      real(wp), intent(in) :: dst_g_bl(:, :), dst_g_cs(:, :)

      integer, parameter :: NXG = 72, NYG = 40
      integer, parameter :: NXT = 36, NYT = 20
      type(hgrid_t) :: grid_t
      type(ocean_metrics_t) :: metrics_t
      type(regrid_weights_t) :: wts_t
      real(wp), allocatable :: dst_t(:, :)
      integer :: ti, tj, i_off, j_off, ierr, i, j

      do tj = 0, 1
         do ti = 0, 1
            i_off = ti*NXT
            j_off = tj*NYT
            call grid_t%init(NXT, NYT, 2, 1000.0_wp, 1000.0_wp)
            grid_t%i_offset_global = i_off
            grid_t%j_offset_global = j_off
            grid_t%nx_global = NXG
            grid_t%ny_global = NYG
            call make_tripolar_metrics(metrics_t, grid_t, 0.0_wp, -90.0_wp, 5.0_wp, 4.0_wp, &
                                       6.378e6_wp, 65.0_wp, 0.0_wp)

            call regrid_weights_bilinear(wts_t, nlon, nlat, lon_c, lat_c, &
                                         grid_t%nx_phys, grid_t%ny_phys, grid_t%nghost, &
                                         grid_t%nx_total, grid_t%ny_total, &
                                         metrics_t%geolonT, metrics_t%geolatT, ierr)
            allocate (dst_t(wts_t%nx, wts_t%ny), source=0.0_wp)
            call apply_global_window(wts_t, src, nlon, nlat, dst_t)
            do j = 1, NYT
               do i = 1, NXT
                  call check(error, dst_t(i, j) == dst_g_bl(i + i_off, j + j_off))
                  if (allocated(error)) then
                     call destroy_cartesian_metrics(metrics_t)
                     return
                  end if
               end do
            end do
            deallocate (dst_t)

            call regrid_weights_conservative(wts_t, nlon, nlat, lon_e, lat_e, &
                                             grid_t%nx_phys, grid_t%ny_phys, grid_t%nghost, &
                                             grid_t%nx_total, grid_t%ny_total, &
                                             metrics_t%geolonBu, metrics_t%geolatBu, &
                                             metrics_t%geolonT, metrics_t%geolatT, ierr)
            allocate (dst_t(wts_t%nx, wts_t%ny), source=0.0_wp)
            call apply_global_window(wts_t, src, nlon, nlat, dst_t)
            do j = 1, NYT
               do i = 1, NXT
                  call check(error, dst_t(i, j) == dst_g_cs(i + i_off, j + j_off))
                  if (allocated(error)) then
                     call destroy_cartesian_metrics(metrics_t)
                     return
                  end if
               end do
            end do
            deallocate (dst_t)

            call destroy_cartesian_metrics(metrics_t)
         end do
      end do
   end subroutine check_tiles_match

   ! =====================================================================
   ! Test 8: the apply kernel on the GPU build matches the host result.
   ! =====================================================================

   subroutine test_apply_device(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(regrid_weights_t) :: wts
      integer, parameter :: NLON = 30, NLAT = 20, NXP = 15, NYP = 10
      real(wp) :: lon_c(NLON), lat_c(NLAT)
      real(wp), allocatable :: src(:, :), dst_host(:, :), dst_dev(:, :)
      integer :: ierr, i, j, jg

      call grid%init(NXP, NYP, 2, 1000.0_wp, 1000.0_wp)
      call make_spherical_metrics(metrics, grid, 10.0_wp, -20.0_wp, 1.0_wp, 1.0_wp, 6.378e6_wp)

      call centers_from_span(NLON, 0.0_wp, 360.0_wp, lon_c)
      call centers_from_span(NLAT, -90.0_wp, 90.0_wp, lat_c)

      call regrid_weights_bilinear(wts, NLON, NLAT, lon_c, lat_c, &
                                   grid%nx_phys, grid%ny_phys, grid%nghost, &
                                   grid%nx_total, grid%ny_total, &
                                   metrics%geolonT, metrics%geolatT, ierr)
      call check(error, ierr, 0)
      if (allocated(error)) then
         call destroy_cartesian_metrics(metrics)
         return
      end if

      allocate (src(wts%nsx, wts%nsy))
      do j = 1, wts%nsy
         jg = j + wts%j_lo - 1
         do i = 1, wts%nsx
            src(i, j) = f_analytic(lon_c(i), lat_c(jg))
         end do
      end do

      allocate (dst_host(wts%nx, wts%ny), source=0.0_wp)
      call regrid_apply_2d_impl(wts%nx, wts%ny, wts%nmax, wts%nsx, wts%nsy, &
                                wts%idx_i, wts%idx_j, wts%w, src, dst_host)

      allocate (dst_dev(wts%nx, wts%ny), source=0.0_wp)
      call wts%enter_data()
      !$acc enter data copyin(src) create(dst_dev)
      call regrid_apply_2d_impl(wts%nx, wts%ny, wts%nmax, wts%nsx, wts%nsy, &
                                wts%idx_i, wts%idx_j, wts%w, src, dst_dev)
      !$acc update self(dst_dev)
      !$acc exit data delete(src, dst_dev)
      call wts%exit_data()

      do j = 1, wts%ny
         do i = 1, wts%nx
            call check(error, abs(dst_dev(i, j) - dst_host(i, j)) < 1.0e-14_wp)
            if (allocated(error)) then
               call destroy_cartesian_metrics(metrics)
               return
            end if
         end do
      end do

      call destroy_cartesian_metrics(metrics)
   end subroutine test_apply_device

end module test_regrid_latlon
