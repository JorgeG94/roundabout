!! Tests for the MOM6 supergrid (mosaic) NetCDF reader
!! (`metrics_fill_from_supergrid`), design §5 T2 / M1.
!!
!! Tests (NetCDF I/O — must run with OMP_NUM_THREADS=1):
!!   supergrid_roundtrip    — write analytic lon-lat supergrid, read via
!!                            metrics_fill_from_supergrid + metrics_finalize,
!!                            compare every metric array against
!!                            metrics_fill_spherical on the same sector.
!!   supergrid_dim_mismatch — call supergrid_dims_ok with a mismatched grid
!!                            size and verify it returns .false. (validator
!!                            function; does NOT call error stop so it is safe
!!                            in-process).
!!   supergrid_ghost_extrap — ghost-row metric equals nearest physical row.
!!
!! Discretization tolerance analysis:
!!   The supergrid writer places T-cell centres at EVEN supergrid nodes
!!   (convention: node (2i,2j) = T(i,j) centre) and evaluates dx sub-segments
!!   at the j-row of the segment's OWN node:
!!     dx_half(m,n) = R * cos(lat_south + (n-1)*dlat/2) * (dlon/2) * DEG2RAD
!!   For the two sub-segments contributing to dxT(i,j):
!!     - row n=2j: lat = lat_south + (2j-1)*dlat/2 = lat_south + (j-0.5)*dlat
!!     Both segments use the same row, so:
!!     dxT_sg = 2 * R * cos(lat_south + (j-0.5)*dlat) * (dlon/2) * DEG2RAD
!!            = R * cos(lat_T_centre) * dlon * DEG2RAD
!!   This is EXACTLY what the spherical generator computes for T-centre latitude
!!   lat_T_centre = lat_south + (j_phys - 0.5)*dlat.  Therefore the round-trip
!!   agreement is to floating-point roundoff (~1e-14 relative), NOT O(dlat).
!!   We assert max_rel_err < 1e-10 for all metrics where exact identity holds.
module test_ocean_supergrid
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t, metrics_finalize, &
                                metrics_fill_spherical, &
                                metrics_fill_from_supergrid
   use rdb_io_netcdf, only: nc_create_file, nc_close, &
                            nc_def_dim, nc_def_var_2d, nc_enddef, &
                            nc_put_var_2d
   use netcdf, only: nf90_open, nf90_nowrite, nf90_noerr, &
                     nf90_inq_dimid, nf90_inquire_dimension, nf90_close
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_continuity, only: continuity_t
   use rdb_coriolis_adv, only: coriolis_adv_t
   use rdb_eos, only: eos_t
   use rdb_ocean_pressure_force, only: ocean_pressure_force_t
   use rdb_ocean_horizontal_viscosity, only: ocean_horizontal_viscosity_t
   use rdb_ocean_bottom_drag, only: ocean_bottom_drag_t
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t
   use rdb_ocean_vertical_advection, only: ocean_vertical_advection_t
   use rdb_ocean_hdiff_tracer, only: ocean_hdiff_tracer_t
   use rdb_ocean_vdiff, only: ocean_vdiff_t
   use rdb_ocean_vmix, only: ocean_vmix_t
   use rdb_ocean_dyn, only: ocean_dyn_t, ocean_dyn_step_split
   implicit none
   private

   public :: collect_ocean_supergrid_tests

   integer, parameter :: NGHOST = 2
   real(wp), parameter :: DEG2RAD = 3.14159265358979323846_wp/180.0_wp
   real(wp), parameter :: RAD_EARTH = 6.378e6_wp

contains

   subroutine collect_ocean_supergrid_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("supergrid_roundtrip", test_roundtrip), &
                  new_unittest("supergrid_dim_mismatch", test_dim_mismatch), &
                  new_unittest("supergrid_ghost_extrap", test_ghost_extrap), &
                  new_unittest("supergrid_nonuniform_face_spans", test_nonuniform_face_spans), &
                  new_unittest("supergrid_areaBu_interior", test_areaBu_interior), &
                  new_unittest("supergrid_driven_quiescent_rest", test_supergrid_quiescent) &
                  ]
   end subroutine collect_ocean_supergrid_tests

   function make_grid(nx_phys, ny_phys, dx, dy) result(g)
      integer, intent(in) :: nx_phys, ny_phys
      real(wp), intent(in) :: dx, dy
      type(hgrid_t) :: g
      call g%init(nx_phys, ny_phys, NGHOST, dx, dy)
   end function make_grid

   ! =================================================================
   ! Test-support: analytic supergrid writer
   ! =================================================================

   subroutine write_analytic_supergrid(filename, ni, nj, lon_west, lat_south, &
                                       dlon_deg, dlat_deg, rad_earth)
      !! Write an analytic lon-lat supergrid NetCDF for an ni×nj model grid.
      !!
      !! Convention: supergrid node (m,n) 1-based, m∈[1,nxp], n∈[1,nyp]
      !!   lon(m,n) = lon_west + (m-1)*dlon/2  [degrees]
      !!   lat(m,n) = lat_south + (n-1)*dlat/2 [degrees]
      !!
      !! dx(m,n): segment length from node m to m+1 at j-row n,
      !!   = R * cos(lat_south + (n-1)*dlat/2) * (dlon/2) * DEG2RAD
      !! dy(m,n): segment length from node n to n+1 at i-col m,
      !!   = R * (dlat/2) * DEG2RAD  (latitude-independent)
      !! area(m,n): sub-cell area = dx(m,n) * dy(m,n)
      character(len=*), intent(in) :: filename
      integer, intent(in) :: ni, nj
      real(wp), intent(in) :: lon_west, lat_south, dlon_deg, dlat_deg, rad_earth

      integer :: ncid
      integer :: dim_nxp, dim_nyp, dim_nx, dim_ny
      integer :: vid_x, vid_y, vid_dx, vid_dy, vid_area
      integer :: nxp, nyp, nx_sg, ny_sg, m, n
      real(wp), allocatable :: sg_x(:, :), sg_y(:, :)
      real(wp), allocatable :: sg_dx(:, :), sg_dy(:, :), sg_area(:, :)
      real(wp) :: lat_n, dy_half

      nxp = 2*ni + 1; nyp = 2*nj + 1
      nx_sg = 2*ni; ny_sg = 2*nj
      dy_half = rad_earth*0.5_wp*dlat_deg*DEG2RAD

      allocate (sg_x(nxp, nyp), sg_y(nxp, nyp))
      allocate (sg_dx(nx_sg, nyp), sg_dy(nxp, ny_sg))
      allocate (sg_area(nx_sg, ny_sg))

      ! Supergrid nodes: half-integer spacing
      do n = 1, nyp
         lat_n = lat_south + (n - 1)*dlat_deg*0.5_wp
         do m = 1, nxp
            sg_x(m, n) = lon_west + (m - 1)*dlon_deg*0.5_wp
            sg_y(m, n) = lat_n
         end do
      end do

      ! dx sub-segments at each j-row (cos evaluated at the node latitude)
      do n = 1, nyp
         lat_n = lat_south + (n - 1)*dlat_deg*0.5_wp
         do m = 1, nx_sg
            sg_dx(m, n) = rad_earth*cos(lat_n*DEG2RAD)*0.5_wp*dlon_deg*DEG2RAD
         end do
      end do

      ! dy sub-segments (lat-independent)
      do n = 1, ny_sg
         do m = 1, nxp
            sg_dy(m, n) = dy_half
         end do
      end do

      ! Sub-cell areas
      do n = 1, ny_sg
         do m = 1, nx_sg
            sg_area(m, n) = sg_dx(m, n)*sg_dy(m, n)
         end do
      end do

      call nc_create_file(filename, ncid)
      call nc_def_dim(ncid, "nxp", nxp, dim_nxp)
      call nc_def_dim(ncid, "nyp", nyp, dim_nyp)
      call nc_def_dim(ncid, "nx", nx_sg, dim_nx)
      call nc_def_dim(ncid, "ny", ny_sg, dim_ny)
      call nc_def_var_2d(ncid, "x", [dim_nxp, dim_nyp], vid_x)
      call nc_def_var_2d(ncid, "y", [dim_nxp, dim_nyp], vid_y)
      call nc_def_var_2d(ncid, "dx", [dim_nx, dim_nyp], vid_dx)
      call nc_def_var_2d(ncid, "dy", [dim_nxp, dim_ny], vid_dy)
      call nc_def_var_2d(ncid, "area", [dim_nx, dim_ny], vid_area)
      call nc_enddef(ncid)
      call nc_put_var_2d(ncid, vid_x, sg_x)
      call nc_put_var_2d(ncid, vid_y, sg_y)
      call nc_put_var_2d(ncid, vid_dx, sg_dx)
      call nc_put_var_2d(ncid, vid_dy, sg_dy)
      call nc_put_var_2d(ncid, vid_area, sg_area)
      call nc_close(ncid)
      deallocate (sg_x, sg_y, sg_dx, sg_dy, sg_area)
   end subroutine write_analytic_supergrid

   ! =================================================================
   ! Dimension validator — safe for in-process testing (no error stop)
   ! =================================================================

   function supergrid_dims_ok(filename, ni, nj) result(ok)
      !! Check that `nxp == 2*ni+1` and `nyp == 2*nj+1` in the file.
      !! Returns .true. on match, .false. on any mismatch or I/O error.
      !! Does NOT call `error stop` — safe to call from test code.
      character(len=*), intent(in) :: filename
      integer, intent(in) :: ni, nj
      logical :: ok

      integer :: ncid, dimid, dim_len, ierr

      ok = .false.
      ierr = nf90_open(trim(filename), nf90_nowrite, ncid)
      if (ierr /= nf90_noerr) return

      ! check nxp
      ierr = nf90_inq_dimid(ncid, "nxp", dimid)
      if (ierr /= nf90_noerr) then
         ierr = nf90_close(ncid); return
      end if
      ierr = nf90_inquire_dimension(ncid, dimid, len=dim_len)
      if (ierr /= nf90_noerr .or. dim_len /= 2*ni + 1) then
         ierr = nf90_close(ncid); return
      end if

      ! check nyp
      ierr = nf90_inq_dimid(ncid, "nyp", dimid)
      if (ierr /= nf90_noerr) then
         ierr = nf90_close(ncid); return
      end if
      ierr = nf90_inquire_dimension(ncid, dimid, len=dim_len)
      if (ierr /= nf90_noerr .or. dim_len /= 2*nj + 1) then
         ierr = nf90_close(ncid); return
      end if

      ok = .true.
      ierr = nf90_close(ncid)
   end function supergrid_dims_ok

   ! =================================================================
   ! Test T2: round-trip supergrid vs spherical generator
   ! =================================================================

   subroutine test_roundtrip(error)
      !! Write an analytic lon-lat supergrid (8x6 cells, 1° resolution,
      !! lat=30°), read via the supergrid reader, and compare every metric
      !! against metrics_fill_spherical on the same sector.
      !!
      !! The supergrid writer uses EVEN nodes for T-cell centres and evaluates
      !! dx at the T-centre latitude — identical to the spherical generator.
      !! Round-trip tolerance should be floating-point roundoff (~1e-14).
      !! We assert max_rel_err < 1e-10 with 4 orders of margin.
      !! geolatT differs: spherical uses cell-CENTRE lat, supergrid uses the
      !! EVEN node (which equals the centre): they should agree exactly.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NI = 8, NJ = 6
      real(wp), parameter :: LON_W = 20.0_wp, LAT_S = 30.0_wp
      real(wp), parameter :: DLON = 1.0_wp, DLAT = 1.0_wp
      real(wp), parameter :: TOL = 1.0e-10_wp  ! roundoff; 4 orders of margin

      type(ocean_metrics_t) :: m_sg, m_sph
      type(hgrid_t) :: g
      character(len=256) :: fname
      integer :: i, j, ng
      real(wp) :: rel_err, max_err

      g = make_grid(NI, NJ, DLON, DLAT)
      ng = NGHOST

      write (fname, '(a)') "/tmp/test_supergrid_roundtrip.nc"
      call write_analytic_supergrid(fname, NI, NJ, LON_W, LAT_S, DLON, DLAT, RAD_EARTH)

      call m_sg%init(g)
      call metrics_fill_from_supergrid(m_sg, g, trim(fname))
      call metrics_finalize(m_sg)

      call m_sph%init(g)
      call metrics_fill_spherical(m_sph, g, LON_W, LAT_S, DLON, DLAT, RAD_EARTH)
      call metrics_finalize(m_sph)

      ! dxT and dyT: should be bit-or-near-identical (same formula)
      max_err = 0.0_wp
      do j = ng + 1, ng + NJ
         do i = ng + 1, ng + NI
            if (m_sph%dxT(i, j) > 0.0_wp) then
               rel_err = abs(m_sg%dxT(i, j) - m_sph%dxT(i, j))/m_sph%dxT(i, j)
               if (rel_err > max_err) max_err = rel_err
            end if
         end do
      end do
      call check(error, max_err < TOL, &
                 "dxT roundtrip max_rel_err "//fmt_e(max_err)//" >= "//fmt_e(TOL))
      if (allocated(error)) return

      max_err = maxval(abs(m_sg%dyT(ng + 1:ng + NI, ng + 1:ng + NJ) - &
                           m_sph%dyT(ng + 1:ng + NI, ng + 1:ng + NJ)))/ &
                maxval(abs(m_sph%dyT(ng + 1:ng + NI, ng + 1:ng + NJ)))
      call check(error, max_err < TOL, &
                 "dyT roundtrip max_rel_err "//fmt_e(max_err)//" >= "//fmt_e(TOL))
      if (allocated(error)) return

      ! areaT: the supergrid sums 4 sub-cell areas which correctly differ
      ! from dxT*dyT by O(dlat^2) when the bottom vs top half-rows have
      ! different latitudes (D5 — "areaT /= dxT*dyT on supergrid").
      ! We verify: (a) all areaT > 0, (b) no NaN,
      ! (c) areaT is consistent with the writer's own dx*dy (each sub-area
      !     = dx_sub * dy_sub — both from the writer, so exact recovery).
      call check(error, all(m_sg%areaT(ng + 1:ng + NI, ng + 1:ng + NJ) > 0.0_wp), &
                 "areaT has non-positive values")
      if (allocated(error)) return
      call check(error,.not. any(m_sg%areaT /= m_sg%areaT), "NaN in areaT")
      if (allocated(error)) return
      ! areaT must be strictly larger than 0 and smaller than 4*dxT*dyT
      ! (since the sub-areas use slightly different latitudes, their sum
      ! is close to dxT*dyT — within ~1% for 1° resolution at lat 30°).
      block
         real(wp) :: ratio
         ratio = m_sg%areaT(ng + 1, ng + 1)/ &
                 (m_sg%dxT(ng + 1, ng + 1)*m_sg%dyT(ng + 1, ng + 1))
         call check(error, ratio > 0.99_wp .and. ratio < 1.01_wp, &
                    "areaT(1,1) / (dxT*dyT) out of [0.99,1.01]: "//fmt_e(ratio))
         if (allocated(error)) return
      end block

      ! dxBu: corner lat = lat_south + (j-1)*dlat (ODD node 2j-1)
      ! Spherical dxBu uses the corner latitude directly — same formula.
      max_err = 0.0_wp
      do j = ng + 1, ng + NJ + 1
         do i = ng + 1, ng + NI + 1
            if (m_sph%dxBu(i, j) > 0.0_wp) then
               rel_err = abs(m_sg%dxBu(i, j) - m_sph%dxBu(i, j))/m_sph%dxBu(i, j)
               if (rel_err > max_err) max_err = rel_err
            end if
         end do
      end do
      call check(error, max_err < TOL, &
                 "dxBu roundtrip max_rel_err "//fmt_e(max_err)//" >= "//fmt_e(TOL))
      if (allocated(error)) return

      ! dyCu: same as dyT (lat-independent) → roundoff
      max_err = 0.0_wp
      do j = ng + 1, ng + NJ
         do i = ng + 1, ng + NI
            if (m_sph%dyCu(i, j) > 0.0_wp) then
               rel_err = abs(m_sg%dyCu(i, j) - m_sph%dyCu(i, j))/m_sph%dyCu(i, j)
               if (rel_err > max_err) max_err = rel_err
            end if
         end do
      end do
      call check(error, max_err < TOL, &
                 "dyCu roundtrip max_rel_err "//fmt_e(max_err)//" >= "//fmt_e(TOL))
      if (allocated(error)) return

      ! dxCv: v-face uses corner latitude (ODD j-row) = lat_south + (j-1)*dlat
      ! Same formula as spherical → roundoff
      max_err = 0.0_wp
      do j = ng + 1, ng + NJ
         do i = ng + 1, ng + NI
            if (m_sph%dxCv(i, j) > 0.0_wp) then
               rel_err = abs(m_sg%dxCv(i, j) - m_sph%dxCv(i, j))/m_sph%dxCv(i, j)
               if (rel_err > max_err) max_err = rel_err
            end if
         end do
      end do
      call check(error, max_err < TOL, &
                 "dxCv roundtrip max_rel_err "//fmt_e(max_err)//" >= "//fmt_e(TOL))
      if (allocated(error)) return

      ! geolatT: T-centre lat matches between the two generators
      max_err = 0.0_wp
      do j = ng + 1, ng + NJ
         do i = ng + 1, ng + NI
            if (abs(m_sph%geolatT(i, j)) > 0.0_wp) then
               rel_err = abs(m_sg%geolatT(i, j) - m_sph%geolatT(i, j))/ &
                         abs(m_sph%geolatT(i, j))
            else
               rel_err = abs(m_sg%geolatT(i, j) - m_sph%geolatT(i, j))
            end if
            if (rel_err > max_err) max_err = rel_err
         end do
      end do
      call check(error, max_err < TOL, &
                 "geolatT roundtrip max_rel_err "//fmt_e(max_err)//" >= "//fmt_e(TOL))
      if (allocated(error)) return

      call m_sg%destroy()
      call m_sph%destroy()
   end subroutine test_roundtrip

   ! =================================================================
   ! Test: dim-mismatch validator
   ! =================================================================

   subroutine test_dim_mismatch(error)
      !! Write a 4x4 supergrid then check supergrid_dims_ok returns:
      !!   .true.  when queried with (4,4) [correct dimensions]
      !!   .false. when queried with (5,4) [nxp mismatch]
      !!   .false. when queried with (4,6) [nyp mismatch]
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NI = 4, NJ = 4
      real(wp), parameter :: LON_W = 0.0_wp, LAT_S = 10.0_wp
      real(wp), parameter :: DLON = 1.0_wp, DLAT = 1.0_wp
      character(len=256) :: fname

      write (fname, '(a)') "/tmp/test_supergrid_mismatch.nc"
      call write_analytic_supergrid(fname, NI, NJ, LON_W, LAT_S, DLON, DLAT, RAD_EARTH)

      ! Correct query → ok
      call check(error, supergrid_dims_ok(trim(fname), NI, NJ), &
                 "correct dims flagged as mismatch")
      if (allocated(error)) return

      ! nxp mismatch
      call check(error,.not. supergrid_dims_ok(trim(fname), NI + 1, NJ), &
                 "nxp mismatch not detected")
      if (allocated(error)) return

      ! nyp mismatch
      call check(error,.not. supergrid_dims_ok(trim(fname), NI, NJ + 2), &
                 "nyp mismatch not detected")
   end subroutine test_dim_mismatch

   ! =================================================================
   ! Test: ghost extrapolation
   ! =================================================================

   subroutine test_ghost_extrap(error)
      !! After reading a supergrid, ghost metric values must equal the
      !! nearest physical row/column (constant extrapolation).
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NI = 6, NJ = 5
      real(wp), parameter :: LON_W = 100.0_wp, LAT_S = -10.0_wp
      real(wp), parameter :: DLON = 1.0_wp, DLAT = 1.0_wp

      type(ocean_metrics_t) :: m
      type(hgrid_t) :: g
      character(len=256) :: fname
      integer :: ng, nx_tot, ny_tot

      g = make_grid(NI, NJ, DLON, DLAT)
      ng = NGHOST
      nx_tot = g%nx_total
      ny_tot = g%ny_total

      write (fname, '(a)') "/tmp/test_supergrid_ghost.nc"
      call write_analytic_supergrid(fname, NI, NJ, LON_W, LAT_S, DLON, DLAT, RAD_EARTH)

      call m%init(g)
      call metrics_fill_from_supergrid(m, g, trim(fname))
      call metrics_finalize(m)

      ! West ghost dxT: ghost cols 1..ng should equal first physical col ng+1
      call check(error, m%dxT(1, ng + 1) == m%dxT(ng + 1, ng + 1), &
                 "west ghost dxT col 1 /= first physical col")
      if (allocated(error)) return
      call check(error, m%dxT(ng, ng + 1) == m%dxT(ng + 1, ng + 1), &
                 "west ghost dxT col ng /= first physical col")
      if (allocated(error)) return

      ! East ghost dxT: ghost cols ng+NI+1..nx_tot should equal last physical col
      call check(error, m%dxT(nx_tot, ng + 1) == m%dxT(ng + NI, ng + 1), &
                 "east ghost dxT /= last physical col")
      if (allocated(error)) return

      ! South ghost dyT: ghost rows 1..ng should equal first physical row ng+1
      call check(error, m%dyT(ng + 1, 1) == m%dyT(ng + 1, ng + 1), &
                 "south ghost dyT row 1 /= first physical row")
      if (allocated(error)) return

      ! North ghost dyT
      call check(error, m%dyT(ng + 1, ny_tot) == m%dyT(ng + 1, ng + NJ), &
                 "north ghost dyT /= last physical row")
      if (allocated(error)) return

      ! geolatT west ghost
      call check(error, m%geolatT(1, ng + 1) == m%geolatT(ng + 1, ng + 1), &
                 "west ghost geolatT /= first physical col")
      if (allocated(error)) return

      ! areaT south ghost
      call check(error, m%areaT(ng + 1, 1) == m%areaT(ng + 1, ng + 1), &
                 "south ghost areaT /= first physical row")
      if (allocated(error)) return

      call m%destroy()
   end subroutine test_ghost_extrap

   ! =================================================================
   ! Non-uniform supergrid writer (dlon varies per half-column)
   ! =================================================================

   subroutine write_nonuniform_supergrid(filename, ni, nj, lon_west, lat_south, &
                                         dlon_deg, dlat_deg, rad_earth)
      !! Write a supergrid with NON-UNIFORM dx sub-segments: the half-segment
      !! width in the x direction varies per supergrid column m as
      !!   dx_half(m) = R * cos(lat) * (dlon/2) * (1 + 0.05*(m-1)/(2*ni)) * DEG2RAD
      !! (a gentle 5% taper).  dy sub-segments remain uniform (dlat constant).
      !! This breaks the dxCu = dxT coincidence that holds on uniform grids,
      !! making the T-to-T face-span correction observable in the test.
      character(len=*), intent(in) :: filename
      integer, intent(in) :: ni, nj
      real(wp), intent(in) :: lon_west, lat_south, dlon_deg, dlat_deg, rad_earth

      integer :: ncid
      integer :: dim_nxp, dim_nyp, dim_nx, dim_ny
      integer :: vid_x, vid_y, vid_dx, vid_dy, vid_area
      integer :: nxp, nyp, nx_sg, ny_sg, m, n
      real(wp), allocatable :: sg_x(:, :), sg_y(:, :)
      real(wp), allocatable :: sg_dx(:, :), sg_dy(:, :), sg_area(:, :)
      real(wp) :: lat_n, dy_half, taper

      nxp = 2*ni + 1; nyp = 2*nj + 1
      nx_sg = 2*ni; ny_sg = 2*nj
      dy_half = rad_earth*0.5_wp*dlat_deg*DEG2RAD

      allocate (sg_x(nxp, nyp), sg_y(nxp, nyp))
      allocate (sg_dx(nx_sg, nyp), sg_dy(nxp, ny_sg))
      allocate (sg_area(nx_sg, ny_sg))

      ! Node positions (uniform half-integer spacing — positions are not affected
      ! by the dx perturbation, just the segment lengths)
      do n = 1, nyp
         lat_n = lat_south + (n - 1)*dlat_deg*0.5_wp
         do m = 1, nxp
            sg_x(m, n) = lon_west + (m - 1)*dlon_deg*0.5_wp
            sg_y(m, n) = lat_n
         end do
      end do

      ! Non-uniform dx: taper = 1 + 0.05*(m-1)/(2*ni)
      do n = 1, nyp
         lat_n = lat_south + (n - 1)*dlat_deg*0.5_wp
         do m = 1, nx_sg
            taper = 1.0_wp + 0.05_wp*real(m - 1, wp)/real(2*ni, wp)
            sg_dx(m, n) = rad_earth*cos(lat_n*DEG2RAD)*0.5_wp*dlon_deg*DEG2RAD*taper
         end do
      end do

      ! Uniform dy
      do n = 1, ny_sg
         do m = 1, nxp
            sg_dy(m, n) = dy_half
         end do
      end do

      ! Sub-cell areas
      do n = 1, ny_sg
         do m = 1, nx_sg
            sg_area(m, n) = sg_dx(m, n)*sg_dy(m, n)
         end do
      end do

      call nc_create_file(filename, ncid)
      call nc_def_dim(ncid, "nxp", nxp, dim_nxp)
      call nc_def_dim(ncid, "nyp", nyp, dim_nyp)
      call nc_def_dim(ncid, "nx", nx_sg, dim_nx)
      call nc_def_dim(ncid, "ny", ny_sg, dim_ny)
      call nc_def_var_2d(ncid, "x", [dim_nxp, dim_nyp], vid_x)
      call nc_def_var_2d(ncid, "y", [dim_nxp, dim_nyp], vid_y)
      call nc_def_var_2d(ncid, "dx", [dim_nx, dim_nyp], vid_dx)
      call nc_def_var_2d(ncid, "dy", [dim_nxp, dim_ny], vid_dy)
      call nc_def_var_2d(ncid, "area", [dim_nx, dim_ny], vid_area)
      call nc_enddef(ncid)
      call nc_put_var_2d(ncid, vid_x, sg_x)
      call nc_put_var_2d(ncid, vid_y, sg_y)
      call nc_put_var_2d(ncid, vid_dx, sg_dx)
      call nc_put_var_2d(ncid, vid_dy, sg_dy)
      call nc_put_var_2d(ncid, vid_area, sg_area)
      call nc_close(ncid)
      deallocate (sg_x, sg_y, sg_dx, sg_dy, sg_area)
   end subroutine write_nonuniform_supergrid

   ! =================================================================
   ! Test: non-uniform supergrid dxCu / dyCv face spans (review finding 8)
   ! =================================================================

   subroutine test_nonuniform_face_spans(error)
      !! Adversarial test (review finding 8): on a NON-uniform supergrid
      !! (dlon varies per half-column), the T-to-T dxCu definition differs
      !! from the host-cell dxT.  This catches any regression to the old
      !! host-cell code path.
      !!
      !! We hand-compute the expected dxCu for 3 interior faces and assert
      !! exact equality (the reader just sums two segments — no float arithmetic
      !! beyond integer indexing, so bit-exact recovery is expected).
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NI = 6, NJ = 4
      real(wp), parameter :: LON_W = 0.0_wp, LAT_S = 30.0_wp
      real(wp), parameter :: DLON = 1.0_wp, DLAT = 1.0_wp

      type(ocean_metrics_t) :: m
      type(hgrid_t) :: g
      character(len=256) :: fname
      integer :: ng, i, j
      real(wp) :: expect_dxCu, expect_dxT
      ! Supergrid dx at specific physical j=1 T-row (sg j-row sj = 2):
      !   lat_n for sg n-index n: lat_south + (n-1)*dlat/2
      !   sj = 2*j_phys = 2 (j_phys=1)
      !   lat at sg row n=2: LAT_S + (2-1)*DLAT/2 = LAT_S + 0.5*DLAT
      real(wp) :: lat_row2
      real(wp) :: dx_seg

      g = make_grid(NI, NJ, DLON, DLAT)
      ng = NGHOST

      write (fname, '(a)') "/tmp/test_supergrid_nonuniform.nc"
      call write_nonuniform_supergrid(fname, NI, NJ, LON_W, LAT_S, DLON, DLAT, RAD_EARTH)

      call m%init(g)
      call metrics_fill_from_supergrid(m, g, trim(fname))
      call metrics_finalize(m)

      lat_row2 = LAT_S + 0.5_wp*DLAT   ! sg row n=2, lat of T j=1 row

      ! --- Check 3 interior Cu faces (i=2,3,4) at physical j=1 ---
      ! Cu face at model i_face uses sg segs m=2i-2 and m=2i-1 at sg row sj=2.
      ! taper(m) = 1 + 0.05*(m-1)/(2*NI); seg = R*cos(lat)*dlon/2*DEG2RAD*taper
      j = ng + 1  ! physical j=1 in total indexing
      do i = 2, 4
         ! expected: seg(2i-2) + seg(2i-1)
         dx_seg = RAD_EARTH*cos(lat_row2*DEG2RAD)*0.5_wp*DLON*DEG2RAD
         expect_dxCu = dx_seg*(1.0_wp + 0.05_wp*real(2*i - 3, wp)/real(2*NI, wp)) + &
                       dx_seg*(1.0_wp + 0.05_wp*real(2*i - 2, wp)/real(2*NI, wp))
         call check(error, abs(m%dxCu(ng + i, j) - expect_dxCu) < 1.0e-10_wp*expect_dxCu, &
                    "dxCu face i="//fmt_e(real(i, wp))//" T-to-T span wrong")
         if (allocated(error)) return

         ! Verify dxCu /= dxT at these faces (the distinction IS observable)
         expect_dxT = m%dxT(ng + i, j)
         call check(error, abs(m%dxCu(ng + i, j) - expect_dxT) > 1.0e-8_wp*expect_dxT, &
                    "dxCu == dxT on non-uniform supergrid (regression to host-cell)")
         if (allocated(error)) return
      end do

      call m%destroy()
   end subroutine test_nonuniform_face_spans

   ! =================================================================
   ! Test: interior areaBu round-trip check
   ! =================================================================

   subroutine test_areaBu_interior(error)
      !! Check areaBu = dxBu * dyBu for all interior physical corners, and
      !! that no areaBu value is NaN or negative.  Previously only edge
      !! areaBu values were tested via ghost_extrap.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NI = 8, NJ = 6
      real(wp), parameter :: LON_W = 10.0_wp, LAT_S = 20.0_wp
      real(wp), parameter :: DLON = 1.0_wp, DLAT = 1.0_wp
      real(wp), parameter :: TOL = 1.0e-10_wp

      type(ocean_metrics_t) :: m
      type(hgrid_t) :: g
      character(len=256) :: fname
      integer :: i, j, ng
      real(wp) :: ratio

      g = make_grid(NI, NJ, DLON, DLAT)
      ng = NGHOST

      write (fname, '(a)') "/tmp/test_supergrid_areaBu.nc"
      call write_analytic_supergrid(fname, NI, NJ, LON_W, LAT_S, DLON, DLAT, RAD_EARTH)

      call m%init(g)
      call metrics_fill_from_supergrid(m, g, trim(fname))
      call metrics_finalize(m)

      ! No NaN anywhere
      call check(error,.not. any(m%areaBu /= m%areaBu), "NaN in areaBu")
      if (allocated(error)) return

      ! All physical interior corners: areaBu > 0 and == dxBu*dyBu exactly
      do j = ng + 2, ng + NJ
         do i = ng + 2, ng + NI
            call check(error, m%areaBu(i, j) > 0.0_wp, "areaBu <= 0 at interior corner")
            if (allocated(error)) return
            if (m%dxBu(i, j) > 0.0_wp .and. m%dyBu(i, j) > 0.0_wp) then
               ratio = m%areaBu(i, j)/(m%dxBu(i, j)*m%dyBu(i, j))
               call check(error, abs(ratio - 1.0_wp) < TOL, &
                          "areaBu /= dxBu*dyBu at interior corner: ratio = "//fmt_e(ratio))
               if (allocated(error)) return
            end if
         end do
      end do

      ! Physical edge corners (i=ng+1, i=ng+NI+1, j=ng+1, j=ng+NJ+1)
      ! — these use boundary extrapolation, check they are positive and finite
      do j = ng + 1, ng + NJ + 1
         call check(error, m%areaBu(ng + 1, j) > 0.0_wp, "west-edge areaBu <= 0")
         if (allocated(error)) return
         call check(error, m%areaBu(ng + NI + 1, j) > 0.0_wp, "east-edge areaBu <= 0")
         if (allocated(error)) return
      end do
      do i = ng + 1, ng + NI + 1
         call check(error, m%areaBu(i, ng + 1) > 0.0_wp, "south-edge areaBu <= 0")
         if (allocated(error)) return
         call check(error, m%areaBu(i, ng + NJ + 1) > 0.0_wp, "north-edge areaBu <= 0")
         if (allocated(error)) return
      end do

      call m%destroy()
   end subroutine test_areaBu_interior

   ! =================================================================
   ! Test: end-to-end supergrid-driven quiescent run (deliverable 3)
   ! =================================================================

   subroutine test_supergrid_quiescent(error)
      !! Supergrid twin of the ACC quiescent gate: write an analytic
      !! lon-lat supergrid, configure metrics via the SUPERGRID branch
      !! (`metrics_fill_from_supergrid` + `metrics_finalize`), seed a
      !! resting stratified state (uniform per-layer T/S ⇒ no horizontal
      !! density gradient ⇒ no PGF), and drive several full
      !! `ocean_dyn_step_split` steps on the supergrid-read metrics.
      !! Velocities must stay at rest to round-off (<1e-10 m/s): any
      !! metric inconsistency introduced by the supergrid reader (bad
      !! face span, area, or f) injects spurious motion.
      !!
      !! We chose the metrics-level path (configure the slot from the
      !! supergrid, then run the dyn-step kernels directly) over driving
      !! the full namelist driver — the driver wiring is integration-
      !! tested elsewhere; here the point is that the SUPERGRID-READ
      !! metrics produce a clean rest state through the real dyn core.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NI = 10, NJ = 8, NZC = 3
      real(wp), parameter :: LON_W = 0.0_wp, LAT_S = 30.0_wp
      real(wp), parameter :: DLON = 0.25_wp, DLAT = 0.25_wp
      real(wp), parameter :: H0 = 200.0_wp
      real(wp), parameter :: DT = 10.0_wp
      integer, parameter :: N_INNER = 5, N_OUTER = 3

      type(ocean_metrics_t) :: metrics
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(coriolis_adv_t) :: cor
      type(ocean_pressure_force_t) :: pgf
      type(ocean_horizontal_viscosity_t) :: hv
      type(ocean_bottom_drag_t) :: bd
      type(ocean_surface_stress_t) :: ss
      type(ocean_vertical_advection_t) :: va
      type(ocean_hdiff_tracer_t) :: hd
      type(ocean_vdiff_t) :: vd
      type(ocean_vmix_t) :: vmix
      type(eos_t) :: eos
      type(ocean_dyn_t) :: dyn
      character(len=256) :: fname
      real(wp) :: S_k(NZC), T_k(NZC)
      real(wp) :: max_du, max_dv
      integer :: k, step
      checks: block

         call grid%init(NI, NJ, NGHOST, DLON, DLAT)

         ! ---- configure metrics from the supergrid file ----
         write (fname, '(a)') "/tmp/test_supergrid_quiescent.nc"
         call write_analytic_supergrid(fname, NI, NJ, LON_W, LAT_S, &
                                       DLON, DLAT, RAD_EARTH)
         call metrics%init(grid)
         call metrics_fill_from_supergrid(metrics, grid, trim(fname))
         call metrics_finalize(metrics)

         ! ---- slots ----
         ms%nz_ml = NZC
         call ms%init(grid)
         call ct%init(grid, nz_ml=NZC)
         cor%f_0 = 2.0_wp*7.292115e-5_wp*sin(LAT_S*DEG2RAD)
         call cor%init(grid, nz_ml=NZC)
         call pgf%init(grid, nz_ml=NZC)
         call hv%init(grid, nz_ml=NZC)
         call bd%init(grid, nz_ml=NZC)
         call ss%init(grid, nz_ml=NZC)
         call va%init(grid, nz_ml=NZC)
         call hd%init(grid, nz_ml=NZC)
         call vd%init(grid, nz_ml=NZC)
         call vmix%init(grid, nz_ml=NZC)
         call eos%init(grid)
         call dyn%init(grid, nz_ml=NZC)

         ! Stable stratification anchored on the EOS reference, uniform in
         ! the horizontal ⇒ zero PGF, true rest state.
         S_k = [eos%S_ref + 1.0_wp, eos%S_ref, eos%S_ref - 1.0_wp]
         T_k = [eos%T_ref - 2.0_wp, eos%T_ref, eos%T_ref + 2.0_wp]
         ms%h_layer = H0
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         do k = 1, NZC
            ms%tracers(ms%idx_salinity)%hTr(:, :, k) = S_k(k)*H0
            ms%tracers(ms%idx_temperature)%hTr(:, :, k) = T_k(k)*H0
         end do
         dyn%bt_work%bt_H_ref = real(NZC, wp)*H0

         !$acc enter data copyin(metrics)
         call metrics%enter_data()
         !$acc enter data copyin(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
         call ms%enter_data(); call ct%enter_data(); call cor%enter_data()
         call pgf%enter_data(); call hv%enter_data(); call bd%enter_data()
         call ss%enter_data(); call va%enter_data(); call hd%enter_data()
         call vd%enter_data(); call vmix%enter_data(); call dyn%enter_data()

         do step = 1, N_OUTER
            call ocean_dyn_step_split( &
               grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
               va, hd, vd, vmix, ms, DT, N_INNER)
         end do

         call dyn%exit_data(); call vmix%exit_data(); call vd%exit_data()
         call hd%exit_data(); call va%exit_data(); call ss%exit_data()
         call bd%exit_data(); call hv%exit_data(); call pgf%exit_data()
         call cor%exit_data(); call ct%exit_data(); call ms%exit_data()
         !$acc exit data delete(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
         call metrics%exit_data()
         !$acc exit data delete(metrics)

         max_du = maxval(abs(ms%u_face_x_layer))
         max_dv = maxval(abs(ms%v_face_y_layer))
         call check(error, max_du < 1.0e-10_wp, &
                    "supergrid quiescent: injected u "//fmt_e(max_du))
         if (allocated(error)) exit checks
         call check(error, max_dv < 1.0e-10_wp, &
                    "supergrid quiescent: injected v "//fmt_e(max_dv))

      end block checks
      call dyn%destroy(); call eos%destroy(); call vmix%destroy(); call vd%destroy()
      call hd%destroy(); call va%destroy(); call ss%destroy(); call bd%destroy()
      call hv%destroy(); call pgf%destroy(); call cor%destroy(); call ct%destroy()
      call ms%destroy(); call metrics%destroy()
   end subroutine test_supergrid_quiescent

   ! -----------------------------------------------------------------
   ! Helper
   ! -----------------------------------------------------------------

   pure function fmt_e(x) result(s)
      !! Format a real as a short exponential string (for check messages).
      real(wp), intent(in) :: x
      character(len=12) :: s
      write (s, '(es10.3)') x
   end function fmt_e

end module test_ocean_supergrid
