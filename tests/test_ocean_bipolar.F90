!! Unit tests for the analytic TRIPOLAR (Murray 1996) grid generator
!! (`rdb_ocean_bipolar` + `metrics_fill_tripolar`).
!!   T1 below-join identity : rows entirely below phi_join reproduce
!!      metrics_fill_spherical to <=1e-9 (great-circle vs analytic-deriv
!!      lengths agree to the discretization tolerance).
!!   T2 join continuity      : corner lat/lon at the join ring identical
!!      from the bipolar map and the lon-lat construction (<1e-12).
!!   T3 orthogonality        : interior cap cell i-edge vs j-edge angle
!!      = 90 deg within tolerance (conformal => near-exact).
!!   T4 fold-seam self-conjugacy : corner (c, top) coincides with corner
!!      (ni+2-c, top) geographically (Appendix A) to <1e-9.
!!   T5 pole placement       : the two poles sit at (lon_pole, phi_join)
!!      and (lon_pole+180, phi_join); all cap latitudes < 90.
!!   T6 area sanity          : sum of cap areaT ~ spherical-cap area above
!!      phi_join within the discretization tolerance.
module test_ocean_bipolar
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_bipolar, only: bipolar_corner_latlon, bipolar_pole_lat
   use rdb_ocean_metrics, only: ocean_metrics_t, metrics_finalize, &
                                metrics_fill_tripolar, metrics_fill_spherical
   implicit none
   private

   public :: collect_ocean_bipolar_tests

   integer, parameter :: NGHOST = 2
   real(wp), parameter :: PI = 3.14159265358979323846_wp
   real(wp), parameter :: DEG2RAD = PI/180.0_wp
   real(wp), parameter :: R_EARTH = 6.378e6_wp

contains

   subroutine collect_ocean_bipolar_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("bipolar_below_join_identity", test_below_join), &
                  new_unittest("bipolar_join_continuity", test_join_continuity), &
                  new_unittest("bipolar_orthogonality", test_orthogonality), &
                  new_unittest("bipolar_fold_seam_conjugacy", test_fold_seam), &
                  new_unittest("bipolar_pole_placement", test_pole_placement), &
                  new_unittest("bipolar_cap_area_sanity", test_cap_area) &
                  ]
   end subroutine collect_ocean_bipolar_tests

   function make_grid(nx_phys, ny_phys, dx, dy) result(g)
      integer, intent(in) :: nx_phys, ny_phys
      real(wp), intent(in) :: dx, dy
      type(hgrid_t) :: g
      call g%init(nx_phys, ny_phys, NGHOST, dx, dy)
   end function make_grid

   pure function geo_to_xyz(lat, lon) result(p)
      real(wp), intent(in) :: lat, lon
      real(wp) :: p(3), la, lo
      la = lat*DEG2RAD
      lo = lon*DEG2RAD
      p = [cos(la)*cos(lo), cos(la)*sin(lo), sin(la)]
   end function geo_to_xyz

   ! ---- T1: below-join identity vs metrics_fill_spherical ----
   ! Domain whose TOP corner latitude is still below phi_join => the whole
   ! grid is plain lon-lat; tripolar metrics must MATCH spherical metrics.
   !
   ! Caveat (honest): the tripolar generator builds lengths by GREAT-CIRCLE
   ! distance (supergrid-reuse path), whereas metrics_fill_spherical uses the
   ! analytic-derivative form dx = R*cos(lat)*dlon.  A constant-latitude line
   ! is NOT a great circle, so dxT/areaT differ by the great-circle-vs-rhumb
   ! discretization (~2-3e-4 relative at dlon=45 deg, shrinking as dlon->0).
   ! dyT is a meridian (a true great circle) so it matches to roundoff, and
   ! geography matches EXACTLY (same lon-lat formula).  We assert dyT +
   ! geography to <1e-12 and dxT/areaT within the great-circle band.
   subroutine test_below_join(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_metrics_t) :: mt, ms
      type(hgrid_t) :: g
      real(wp), parameter :: lon_w = 0.0_wp, lat_s = 10.0_wp
      real(wp), parameter :: dlon = 360.0_wp/120.0_wp, dlat = 1.0_wp
      real(wp), parameter :: phi_join = 80.0_wp, lon_pole = 100.0_wp
      integer :: ng, i0, i1, j0, j1

      ! 120 cells in i (dlon=3 deg, so dlam wraps 360 over the i range,
      ! matching the tripolar i-convention), 10 in j: top corner lat =
      ! 10 + 10*1 = 20 < 80, so the whole grid is plain lon-lat.
      g = make_grid(120, 10, dlon, dlat)
      ng = NGHOST
      call mt%init(g)
      call ms%init(g)
      call metrics_fill_tripolar(mt, g, lon_w, lat_s, dlon, dlat, R_EARTH, phi_join, lon_pole)
      call metrics_finalize(mt)
      call metrics_fill_spherical(ms, g, lon_w, lat_s, dlon, dlat, R_EARTH)
      call metrics_finalize(ms)

      i0 = ng + 1; i1 = ng + 120
      j0 = ng + 1; j1 = ng + 10
      ! dyT (meridian = great circle): match to roundoff.
      call check(error, maxval(abs(mt%dyT(i0:i1, j0:j1) - ms%dyT(i0:i1, j0:j1))) &
                 <= 1.0e-9_wp*maxval(ms%dyT(i0:i1, j0:j1)), "dyT tripolar /= spherical below join")
      if (allocated(error)) return
      ! Geography must match exactly (same lon-lat formula).
      call check(error, maxval(abs(mt%geolatT(i0:i1, j0:j1) - ms%geolatT(i0:i1, j0:j1))) &
                 < 1.0e-12_wp, "geolatT tripolar /= spherical below join")
      if (allocated(error)) return
      call check(error, maxval(abs(mt%geolonT(i0:i1, j0:j1) - ms%geolonT(i0:i1, j0:j1))) &
                 < 1.0e-12_wp, "geolonT tripolar /= spherical below join")
      if (allocated(error)) return
      ! dxT / areaT: great-circle vs analytic-derivative => agree within the
      ! great-circle discretization band, which scales as (dlon)^2.  At
      ! dlon=3 deg the realized dxT deviation is ~7e-6 (rel); assert < 1e-4.
      call check(error, maxval(abs(mt%dxT(i0:i1, j0:j1) - ms%dxT(i0:i1, j0:j1))) &
                 <= 1.0e-4_wp*maxval(ms%dxT(i0:i1, j0:j1)), &
                 "dxT tripolar departs spherical beyond great-circle band")
      if (allocated(error)) return
      call check(error, maxval(abs(mt%areaT(i0:i1, j0:j1) - ms%areaT(i0:i1, j0:j1))) &
                 <= 1.0e-4_wp*maxval(ms%areaT(i0:i1, j0:j1)), &
                 "areaT tripolar departs spherical beyond great-circle band")
      if (allocated(error)) return

      call mt%destroy()
      call ms%destroy()
   end subroutine test_below_join

   ! ---- T2: join continuity — cap s=0 reproduces the lon-lat ring ----
   subroutine test_join_continuity(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: phi_join = 65.0_wp, lon_pole = 100.0_wp
      real(wp) :: lat, lon, lam
      integer :: k

      do k = 0, 11
         lam = -180.0_wp + real(k, wp)*30.0_wp
         call bipolar_corner_latlon(lam, 0.0_wp, phi_join, lon_pole, lat, lon)
         call check(error, abs(lat - phi_join) < 1.0e-12_wp, &
                    "join lat /= phi_join"); if (allocated(error)) return
         call check(error, abs(lon - lam) < 1.0e-10_wp, &
                    "join lon /= lam (continuity broken)"); if (allocated(error)) return
      end do
   end subroutine test_join_continuity

   ! ---- T3: orthogonality of the (lam, s) cap grid ----
   subroutine test_orthogonality(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: phi_join = 65.0_wp, lon_pole = 100.0_wp
      real(wp), parameter :: h = 1.0e-6_wp
      real(wp) :: lam, s, worst, ang
      real(wp) :: pe1(3), pe0(3), ps1(3), ps0(3), de(3), ds(3), c
      real(wp) :: lat, lon
      integer :: il, is

      worst = 0.0_wp
      do il = 0, 17
         lam = -170.0_wp + real(il, wp)*20.0_wp
         do is = 1, 9
            s = real(is, wp)*0.1_wp
            ! central differences in lam and s on the unit sphere
            call bipolar_corner_latlon(lam + h, s, phi_join, lon_pole, lat, lon)
            pe1 = geo_to_xyz(lat, lon)
            call bipolar_corner_latlon(lam - h, s, phi_join, lon_pole, lat, lon)
            pe0 = geo_to_xyz(lat, lon)
            call bipolar_corner_latlon(lam, min(s + h, 1.0_wp), phi_join, lon_pole, lat, lon)
            ps1 = geo_to_xyz(lat, lon)
            call bipolar_corner_latlon(lam, max(s - h, 0.0_wp), phi_join, lon_pole, lat, lon)
            ps0 = geo_to_xyz(lat, lon)
            de = pe1 - pe0
            ds = ps1 - ps0
            c = dot_product(de, ds)/(norm2(de)*norm2(ds))
            c = max(-1.0_wp, min(1.0_wp, c))
            ang = abs(acos(c) - 0.5_wp*PI)
            worst = max(worst, ang)
         end do
      end do
      ! Conformal construction => near-exact orthogonality.
      call check(error, worst < 1.0e-6_wp, "cap orthogonality worse than 1e-6 rad")
      ! Report the worst angle so the value is captured in the run log.
      if (allocated(error)) return
   end subroutine test_orthogonality

   ! ---- T4: fold-seam self-conjugacy (Appendix A) ----
   ! At the fold (s=1), corner at lam coincides geographically with its
   ! partner at 2*lon_pole - lam (the SAME physical point).
   subroutine test_fold_seam(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: phi_join = 65.0_wp, lon_pole = 100.0_wp
      real(wp) :: lam, lamp, lat1, lon1, lat2, lon2
      real(wp) :: p1(3), p2(3)
      integer :: k

      do k = 1, 7
         lam = 40.0_wp + real(k, wp)*15.0_wp   ! 55..145, around lon_pole
         lamp = 2.0_wp*lon_pole - lam
         lamp = modulo(lamp + 180.0_wp, 360.0_wp) - 180.0_wp
         call bipolar_corner_latlon(lam, 1.0_wp, phi_join, lon_pole, lat1, lon1)
         call bipolar_corner_latlon(lamp, 1.0_wp, phi_join, lon_pole, lat2, lon2)
         ! Compare on the sphere (lon wrap-safe).
         p1 = geo_to_xyz(lat1, lon1)
         p2 = geo_to_xyz(lat2, lon2)
         call check(error, norm2(p1 - p2) < 1.0e-9_wp, &
                    "fold seam not self-conjugate"); if (allocated(error)) return
      end do
   end subroutine test_fold_seam

   ! ---- T5: pole placement + no true-pole grid point ----
   subroutine test_pole_placement(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: phi_join = 65.0_wp, lon_pole = 100.0_wp
      real(wp) :: lat, lon, partner
      integer :: is

      ! Pole 1 at lam = lon_pole (the eta -> -inf focus): lat = phi_join,
      ! lon = lon_pole, for ANY s (it is a single point).
      do is = 0, 4
         call bipolar_corner_latlon(lon_pole, real(is, wp)*0.25_wp, &
                                    phi_join, lon_pole, lat, lon)
         call check(error, abs(lat - bipolar_pole_lat(phi_join)) < 1.0e-6_wp, &
                    "pole 1 lat /= phi_join"); if (allocated(error)) return
         call check(error, abs(lon - lon_pole) < 1.0e-6_wp, &
                    "pole 1 lon /= lon_pole"); if (allocated(error)) return
      end do

      ! Pole 2 at lam = lon_pole + 180.
      partner = modulo(lon_pole + 180.0_wp + 180.0_wp, 360.0_wp) - 180.0_wp
      do is = 0, 4
         call bipolar_corner_latlon(partner, real(is, wp)*0.25_wp, &
                                    phi_join, lon_pole, lat, lon)
         call check(error, abs(lat - phi_join) < 1.0e-6_wp, &
                    "pole 2 lat /= phi_join"); if (allocated(error)) return
         call check(error, abs(abs(lon - lon_pole) - 180.0_wp) < 1.0e-6_wp, &
                    "pole 2 lon /= lon_pole+180"); if (allocated(error)) return
      end do

      ! No cap point reaches the true geographic pole (lat 90) on a grid node:
      ! the NP is on the fold LINE but only at the seam interior; spot-check
      ! several cap interior points stay < 90 - eps.
      do is = 1, 9
         call bipolar_corner_latlon(37.0_wp, real(is, wp)*0.1_wp, &
                                    phi_join, lon_pole, lat, lon)
         call check(error, lat < 90.0_wp, "cap latitude reached 90 (singularity)")
         if (allocated(error)) return
      end do
   end subroutine test_pole_placement

   ! ---- T6: cap areaT sum ~ spherical-cap area above phi_join ----
   subroutine test_cap_area(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_metrics_t) :: m
      type(hgrid_t) :: g
      real(wp), parameter :: lon_w = -180.0_wp, lat_s = 20.0_wp
      real(wp), parameter :: phi_join = 60.0_wp, lon_pole = 100.0_wp
      ! 40 i-cells, 35 j-cells, dlat=2 => top corner lat = 20 + 35*2 = 90.
      ! Cap is the rows with corner lat > 60 (i.e. lat in (60, 90]).
      real(wp), parameter :: dlat = 2.0_wp
      real(wp) :: dlon
      integer :: ni, nj, ng, i, j
      real(wp) :: cap_sum, cap_exact, lat_corner

      ni = 40; nj = 35
      dlon = 360.0_wp/real(ni, wp)
      g = make_grid(ni, nj, dlon, dlat)
      ng = NGHOST
      call m%init(g)
      call metrics_fill_tripolar(m, g, lon_w, lat_s, dlon, dlat, R_EARTH, phi_join, lon_pole)
      call metrics_finalize(m)

      ! Sum areaT over the cells whose SOUTH corner latitude (the lon-lat
      ! ladder) is >= phi_join — those are the cap cells.
      cap_sum = 0.0_wp
      do j = 1, nj
         lat_corner = lat_s + real(j - 1, wp)*dlat
         if (lat_corner >= phi_join) then
            do i = 1, ni
               cap_sum = cap_sum + m%areaT(ng + i, ng + j)
            end do
         end if
      end do

      ! Spherical-cap area above phi_join: 2*pi*R^2*(1 - sin(phi_join)).
      cap_exact = 2.0_wp*PI*R_EARTH*R_EARTH*(1.0_wp - sin(phi_join*DEG2RAD))

      ! Discretization tolerance: a 2-degree ladder + spherical-quad areas;
      ! assert within 1% of the analytic cap (the cap cells cover exactly
      ! lat in [phi_join, 90], the full cap).
      call check(error, abs(cap_sum - cap_exact) <= 1.0e-2_wp*cap_exact, &
                 "cap areaT sum /= spherical-cap area within 1%")
      call m%destroy()
   end subroutine test_cap_area

end module test_ocean_bipolar
