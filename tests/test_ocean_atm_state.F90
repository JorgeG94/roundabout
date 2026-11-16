!! Tests for the JRA55-do-style prescribed atmospheric-state provider
!! (`rdb_atm_state`, OM3 wave-3 C2c-2a/C2c-2b).  Synthetic files written
!! IN-TEST (NetCDF I/O -- must run with OMP_NUM_THREADS=1).
!!
!! Every synthetic field is f(lon, lat, t) = A + B*lon + C*lat + D*t
!! (degrees for lon/lat) -- affine, so BOTH the bilinear spatial regrid
!! and the linear time-blend reproduce it exactly (to round-off)
!! regardless of where the query point/time falls inside the source
!! window, which is what makes an END-TO-END analytic check possible
!! without re-deriving `rdb_regrid_latlon`'s own accuracy (covered by
!! `test_regrid_latlon`).
!!
!! Tests:
!!   atm_state_time_interp_exactness  -- §5 test 1a: single file, every
!!                                       one of the nine fields recovers
!!                                       A+B*lon+C*lat+D*t at an off-node
!!                                       query time, <= 1e-10 absolute.
!!   atm_state_file_boundary_bracket  -- §5 test 1b: two files, query
!!                                       straddles the file boundary.
!!   atm_state_ryf_wrap_bitid         -- §5 test 1c: time_mode=
!!                                       "repeat_year"; querying one full
!!                                       period later reproduces the same
!!                                       wrapped value (mirrors
!!                                       `test_ocean_data_input`'s
!!                                       `test_filelist_ryf_wrap`).
!!   atm_state_cadence_asymmetry      -- §5 test 2: the five states (one
!!                                       time axis) and four fluxes
!!                                       (a DIFFERENT time axis) blended
!!                                       at the same query time each
!!                                       recover their OWN analytic
!!                                       value -- proof every field uses
!!                                       its own reader id/bracket, not
!!                                       a shared one.
!!   atm_state_gpu_resident           -- §5 test 7: enter_data + update
!!                                       entirely under `!$acc`, `bytes()`
!!                                       reconciled against the mapped
!!                                       arrays' own sizes (mem:separate
!!                                       contract template).
!!   atm_state_disabled_bit_identity  -- §5 test 8: `enable = .false.`
!!                                       leaves every slot untouched
!!                                       (`this%enable` stays `.false.`,
!!                                       `atm_state_update` a no-op).
!!   atm_sw_held_record_edge          -- §5 test 3a: `sw_interp="coszen"`,
!!                                       `sw_offset_s=+5400` -- every
!!                                       query in bucket k reads
!!                                       `rsds_rec` as EXACTLY record k
!!                                       (the record stamped at that
!!                                       bucket's centre; held, not
!!                                       blended toward the next one).
!!   atm_sw_energy_identity           -- §5 test 3b: a day at a sunset +
!!                                       sunrise cell, dt = 600/1800/3600
!!                                       -- `Sum rsds*dt` over EACH 3 h
!!                                       bucket equals `rsds_rec*3h` to
!!                                       <= 1e-9 relative, and `rsds <=
!!                                       rsds_rec*3h/dt` everywhere.
!!   atm_sw_failing_geometry          -- regression: 34.7S 74.5W, day of
!!                                       year 3, the bucket holding the
!!                                       local sunset (the cell that blew
!!                                       the 1-degree Southern Ocean up
!!                                       at t = 48 h with `rsds` ~ 1e4 x
!!                                       the record).
!!   atm_sw_polar                     -- polar night: `rsds` exactly 0;
!!                                       polar day: energy-exact, bounded
!!                                       and smooth step to step.
!!   atm_sw_record_alignment          -- a record that is non-zero in
!!                                       exactly ONE bucket deposits its
!!                                       energy in that bucket only.
!!   atm_sw_axis_checks               -- configure refuses `sw_offset_s = -5400`
!!                                       and a spacing/offset mismatch;
!!                                       6-hourly records run exactly.
!!   atm_sw_shape                     -- step-mean `cosz` against an
!!                                       independent 1e4-node quadrature
!!                                       (sunset step, lon -280, a short
!!                                       night, polar night/day).
!!   atm_sw_coszen_analytic           -- §5 test 3d: equinox, lat=lon=0
!!                                       -- local noon gives `cosz`
!!                                       within 1e-3 of 1.0, midnight
!!                                       gives `cosz = 0` (clamped).
!!
!!   §5 test 3c ("sw_interp='linear' reproduces test 1") needs no
!!   separate test: `test_time_interp_exactness` above already runs with
!!   `base_cfg`'s `sw_interp = "linear"` and asserts `rsds` (not
!!   `rsds_rec`) against the plain affine blend -- unchanged by C2c-2b.
!!
!! Scope note: only the fp64 exactness tolerance (plan §5 test 1: "<=
!! 1e-12 relative on the fp64 variant") is exercised here -- the fp32
!! ("<= 2e-6") file variant is not written by this chunk's test file
!! (budget; every synthetic file here is `NC_WP` = working precision).
!! `fluxes_regrid = "bilinear"` throughout (the conservative weight
!! MATH is `test_regrid_latlon`'s job; these tests exercise the
!! reader-blend -> regrid -> scatter WIRING only).
module test_ocean_atm_state
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_config, only: config_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_spherical_metrics, destroy_cartesian_metrics
   use rdb_ocean_data_input, only: ocean_data_input_t, ocean_data_input_update_all
   use rdb_atm_state, only: atm_state_t, atm_state_configure, atm_state_update, atm_sw_window
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   use rdb_io_netcdf, only: nc_create_file, nc_close, nc_def_dim, nc_def_var_3d, &
                            nc_def_var_2d, rdb_def_var_1d, nc_enddef, rdb_put_var_1d, &
                            nc_put_att
   use netcdf, only: nf90_put_var
   implicit none
   private

   public :: collect_ocean_atm_state_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NX_PHYS = 3, NY_PHYS = 2
   real(wp), parameter :: LON_WEST = 20.0_wp, LAT_SOUTH = -5.0_wp
   real(wp), parameter :: DLON_MODEL = 10.0_wp, DLAT_MODEL = 10.0_wp
   real(wp), parameter :: RAD_EARTH = 6.378e6_wp
   integer, parameter :: NLON_SRC = 8, NLAT_SRC = 6
   real(wp), parameter :: TOL_EXACT = 1.0e-10_wp
   real(wp), parameter :: TEST_DT = 1800.0_wp
      !! Model step handed to `atm_state_update` where the test does not
      !! care about the coszen window (`sw_interp = "linear"` ignores it).
   real(wp), parameter :: BUCKET_S = 10800.0_wp

contains

   subroutine collect_ocean_atm_state_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("atm_state_time_interp_exactness", &
                               test_time_interp_exactness), &
                  new_unittest("atm_state_file_boundary_bracket", &
                               test_file_boundary_bracket), &
                  new_unittest("atm_state_ryf_wrap_bitid", test_ryf_wrap_bitid), &
                  new_unittest("atm_state_cadence_asymmetry", test_cadence_asymmetry), &
                  new_unittest("atm_state_gpu_resident", test_gpu_resident), &
                  new_unittest("atm_state_disabled_bit_identity", &
                               test_disabled_bit_identity), &
                  new_unittest("atm_sw_held_record_edge", test_sw_held_record_edge), &
                  new_unittest("atm_sw_energy_identity", test_sw_energy_identity), &
                  new_unittest("atm_sw_failing_geometry", test_sw_failing_geometry), &
                  new_unittest("atm_sw_polar", test_sw_polar), &
                  new_unittest("atm_sw_record_alignment", test_sw_record_alignment), &
                  new_unittest("atm_sw_axis_checks", test_sw_axis_checks), &
                  new_unittest("atm_sw_shape", test_sw_shape), &
                  new_unittest("atm_sw_coszen_analytic", test_sw_coszen_analytic) &
                  ]
   end subroutine collect_ocean_atm_state_tests

   ! =====================================================================
   ! Test-support
   ! =====================================================================

   pure subroutine uniform_axis(n, lo, hi, c, bnds)
      !! Uniform centers + CF bnds(2,n) over [lo, hi].
      integer, intent(in) :: n
      real(wp), intent(in) :: lo, hi
      real(wp), intent(out) :: c(n), bnds(2, n)
      real(wp) :: dx
      integer :: k
      dx = (hi - lo)/real(n, wp)
      do k = 1, n
         bnds(1, k) = lo + real(k - 1, wp)*dx
         bnds(2, k) = lo + real(k, wp)*dx
         c(k) = 0.5_wp*(bnds(1, k) + bnds(2, k))
      end do
   end subroutine uniform_axis

   subroutine write_atm_grid_file(filename, nlon, nlat, lon_c, lat_c, lon_bnds, lat_bnds, &
                                  nt, t_axis, tunits, varnames, coefs, t_data)
      !! One NetCDF file carrying lon/lat/lon_bnds/lat_bnds/time plus
      !! `size(varnames)` field variables, each
      !! `f(lon,lat,t) = coefs(1,v) + coefs(2,v)*lon + coefs(3,v)*lat +
      !! coefs(4,v)*t`, where `t` is `t_data` when present (the ABSOLUTE
      !! time the analytic check uses) or `t_axis` itself otherwise --
      !! distinct for a multi-file axis, where `t_axis` is this file's
      !! own LOCAL record values (relative to `tunits`'s epoch) but the
      !! data must be a function of the GLOBAL absolute time for the
      !! cross-file check to recover a single continuous `f`.
      character(len=*), intent(in) :: filename
      integer, intent(in) :: nlon, nlat, nt
      real(wp), intent(in) :: lon_c(nlon), lat_c(nlat)
      real(wp), intent(in) :: lon_bnds(2, nlon), lat_bnds(2, nlat)
      real(wp), intent(in) :: t_axis(nt)
      character(len=*), intent(in) :: tunits
      character(len=*), intent(in) :: varnames(:)
      real(wp), intent(in) :: coefs(4, size(varnames))
      real(wp), intent(in), optional :: t_data(nt)

      integer :: ncid, dim_lon, dim_lat, dim_time, dim_bnds
      integer :: vid_lon, vid_lat, vid_lonb, vid_latb, vid_time
      integer, allocatable :: vid_f(:)
      real(wp), allocatable :: f(:, :, :)
      real(wp) :: t_for_data(nt)
      integer :: nv, v, i, j, k, ierr

      nv = size(varnames)
      allocate (vid_f(nv))
      allocate (f(nlon, nlat, nt))
      if (present(t_data)) then
         t_for_data = t_data
      else
         t_for_data = t_axis
      end if

      call nc_create_file(filename, ncid)
      call nc_def_dim(ncid, "lon", nlon, dim_lon)
      call nc_def_dim(ncid, "lat", nlat, dim_lat)
      call nc_def_dim(ncid, "bnds", 2, dim_bnds)
      call nc_def_dim(ncid, "time", nt, dim_time)
      call rdb_def_var_1d(ncid, "lon", dim_lon, vid_lon)
      call rdb_def_var_1d(ncid, "lat", dim_lat, vid_lat)
      call nc_def_var_2d(ncid, "lon_bnds", [dim_bnds, dim_lon], vid_lonb)
      call nc_def_var_2d(ncid, "lat_bnds", [dim_bnds, dim_lat], vid_latb)
      call rdb_def_var_1d(ncid, "time", dim_time, vid_time)
      call nc_put_att(ncid, vid_time, "units", tunits)
      do v = 1, nv
         call nc_def_var_3d(ncid, trim(varnames(v)), [dim_lon, dim_lat, dim_time], vid_f(v))
      end do
      call nc_enddef(ncid)

      call rdb_put_var_1d(ncid, vid_lon, lon_c)
      call rdb_put_var_1d(ncid, vid_lat, lat_c)
      ierr = nf90_put_var(ncid, vid_lonb, lon_bnds)
      ierr = nf90_put_var(ncid, vid_latb, lat_bnds)
      call rdb_put_var_1d(ncid, vid_time, t_axis)

      do v = 1, nv
         do k = 1, nt
            do j = 1, nlat
               do i = 1, nlon
                  f(i, j, k) = coefs(1, v) + coefs(2, v)*lon_c(i) + coefs(3, v)*lat_c(j) + &
                               coefs(4, v)*t_for_data(k)
               end do
            end do
         end do
         ierr = nf90_put_var(ncid, vid_f(v), f)
      end do
      call nc_close(ncid)
   end subroutine write_atm_grid_file

   subroutine write_filelist_txt(path, nc_paths)
      character(len=*), intent(in) :: path
      character(len=*), intent(in) :: nc_paths(:)
      integer :: u, k
      open (newunit=u, file=trim(path), status='replace', action='write')
      do k = 1, size(nc_paths)
         write (u, '(A)') trim(nc_paths(k))
      end do
      close (u)
   end subroutine write_filelist_txt

   subroutine make_model(grid, metrics)
      !! Small spherical T-grid whose physical-cell lon/lat all land
      !! strictly inside the source window's interior brackets (never at
      !! a node, never clamped at a source edge) -- see the module
      !! docstring for why that makes the affine field exactly
      !! recoverable regardless of the regrid weights' actual values.
      type(hgrid_t), intent(out) :: grid
      type(ocean_metrics_t), intent(inout) :: metrics
      call grid%init(NX_PHYS, NY_PHYS, NGHOST, 1000.0_wp, 1000.0_wp)
      call make_spherical_metrics(metrics, grid, LON_WEST, LAT_SOUTH, &
                                  DLON_MODEL, DLAT_MODEL, RAD_EARTH)
   end subroutine make_model

   subroutine base_cfg(cfg, start_date, filelist_path)
      !! Every `&ocean_atm_forcing_nml` knob the provider reads, pointed
      !! at ONE shared filelist for all nine fields (the synthetic files
      !! below always carry all nine variable names). `var_<field>`
      !! keeps its default (== the JRA55-do tag), matching the variable
      !! names `write_atm_grid_file` writes.
      type(config_t), intent(inout) :: cfg
      character(len=*), intent(in) :: start_date, filelist_path
      cfg%start_date = start_date
      cfg%ocean%atm_forcing%enable = .true.
      cfg%ocean%atm_forcing%sw_interp = "linear"
      cfg%ocean%atm_forcing%fluxes_regrid = "bilinear"
      cfg%ocean%atm_forcing%time_mode = "interannual"
      cfg%ocean%atm_forcing%filelist_tas = filelist_path
      cfg%ocean%atm_forcing%filelist_huss = filelist_path
      cfg%ocean%atm_forcing%filelist_uas = filelist_path
      cfg%ocean%atm_forcing%filelist_vas = filelist_path
      cfg%ocean%atm_forcing%filelist_psl = filelist_path
      cfg%ocean%atm_forcing%filelist_rlds = filelist_path
      cfg%ocean%atm_forcing%filelist_rsds = filelist_path
      cfg%ocean%atm_forcing%filelist_prra = filelist_path
      cfg%ocean%atm_forcing%filelist_prsn = filelist_path
   end subroutine base_cfg

   pure function query_lon_lat(ip, jp, metrics) result(ll)
      !! Physical T-cell (ip,jp) (1-based, PHYSICAL index) lon/lat.
      integer, intent(in) :: ip, jp
      type(ocean_metrics_t), intent(in) :: metrics
      real(wp) :: ll(2)
      ll(1) = metrics%geolonT(ip + NGHOST, jp + NGHOST)
      ll(2) = metrics%geolatT(ip + NGHOST, jp + NGHOST)
   end function query_lon_lat

   ! =====================================================================
   ! §5 test 1a: time-interp exactness, every one of the nine fields.
   ! =====================================================================

   subroutine test_time_interp_exactness(error)
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: VARNAMES(9) = &
                                     [character(len=4) :: "tas", "huss", "uas", "vas", "psl", &
                                                           "rlds", "rsds", "prra", "prsn"]
      real(wp), parameter :: A = 1.0_wp, B = 0.02_wp, C = 0.05_wp, D = 0.4_wp

      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(config_t) :: cfg
      type(ocean_data_input_t) :: reader
      type(atm_state_t) :: atm
      real(wp) :: lon_c(NLON_SRC), lat_c(NLAT_SRC)
      real(wp) :: lon_bnds(2, NLON_SRC), lat_bnds(2, NLAT_SRC)
      real(wp) :: t_axis(3), coefs(4, 9)
      real(wp) :: tq, ll(2), expect
      integer :: v, i, j, ierr

      call uniform_axis(NLON_SRC, 0.0_wp, 360.0_wp, lon_c, lon_bnds)
      call uniform_axis(NLAT_SRC, -90.0_wp, 90.0_wp, lat_c, lat_bnds)
      t_axis = [0.0_wp, 1.0_wp, 2.0_wp]
      do v = 1, 9
         coefs(:, v) = [A + real(v, wp), B, C, D]
      end do

      call write_atm_grid_file("tmp_local_artifacts/atm_exact.nc", NLON_SRC, NLAT_SRC, &
                               lon_c, lat_c, lon_bnds, lat_bnds, 3, t_axis, &
                               "days since 2000-01-01 00:00:00", VARNAMES, coefs)
      call write_filelist_txt("tmp_local_artifacts/atm_exact_filelist.txt", &
                              ["tmp_local_artifacts/atm_exact.nc"])

      call make_model(grid, metrics)
      call base_cfg(cfg, "2000-01-01", "tmp_local_artifacts/atm_exact_filelist.txt")

      call reader%init()
      call atm_state_configure(cfg, reader, grid, metrics, atm, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "atm_state_configure failed")
      if (allocated(error)) return
      call map_atm(reader, atm)

      tq = 1.37_wp*86400.0_wp
      call ocean_data_input_update_all(reader, tq)
      call atm_state_update(atm, reader, metrics, tq, TEST_DT)
      call sync_atm(atm)

      do j = 1, NY_PHYS
         do i = 1, NX_PHYS
            ll = query_lon_lat(i, j, metrics)
            do v = 1, 9
               expect = coefs(1, v) + coefs(2, v)*ll(1) + coefs(3, v)*ll(2) + &
                        coefs(4, v)*(tq/86400.0_wp)
               call check(error, abs(field_value(atm, v, i + NGHOST, j + NGHOST) - expect) &
                          < TOL_EXACT, "exactness: field "//trim(VARNAMES(v))//" wrong")
               if (allocated(error)) return
            end do
         end do
      end do

      call unmap_atm(reader, atm)
      call destroy_cartesian_metrics(metrics)
      call atm%destroy()
      call reader%destroy()
   end subroutine test_time_interp_exactness

   subroutine map_atm(reader, atm)
      !! Map the reader and the atm state for a device-side update, as the
      !! engine's `enter_data` does (mem:separate: nothing reaches the
      !! device implicitly).
      type(ocean_data_input_t), intent(inout) :: reader
      type(atm_state_t), intent(inout) :: atm
      !$acc enter data copyin(reader)
      call reader%enter_data()
      !$acc enter data copyin(atm)
      call atm%enter_data()
   end subroutine map_atm

   subroutine sync_atm(atm)
      !! Bring the nine regridded fields back to the host for the asserts
      !! (component arrays only, never the aggregate).
      type(atm_state_t), intent(inout) :: atm
      !$acc update self(atm%tas, atm%huss, atm%uas, atm%vas, atm%psl)
      !$acc update self(atm%rlds, atm%rsds, atm%prra, atm%prsn)
      !$acc update self(atm%rsds_rec, atm%cosz, atm%cosz_mean)
   end subroutine sync_atm

   subroutine unmap_atm(reader, atm)
      !! Inverse of `map_atm`.
      type(ocean_data_input_t), intent(inout) :: reader
      type(atm_state_t), intent(inout) :: atm
      call atm%exit_data()
      !$acc exit data delete(atm)
      call reader%exit_data()
      !$acc exit data delete(reader)
   end subroutine unmap_atm

   pure function field_value(atm, v, i, j) result(val)
      !! Index into the nine named field arrays by VARNAMES order --
      !! avoids a pointer array over non-contiguous allocatable
      !! components (array-of-derived-type indirection rule, CLAUDE.md).
      type(atm_state_t), intent(in) :: atm
      integer, intent(in) :: v, i, j
      real(wp) :: val
      select case (v)
      case (1); val = atm%tas(i, j)
      case (2); val = atm%huss(i, j)
      case (3); val = atm%uas(i, j)
      case (4); val = atm%vas(i, j)
      case (5); val = atm%psl(i, j)
      case (6); val = atm%rlds(i, j)
      case (7); val = atm%rsds(i, j)
      case (8); val = atm%prra(i, j)
      case (9); val = atm%prsn(i, j)
      case default; val = huge(1.0_wp)
      end select
   end function field_value

   ! =====================================================================
   ! §5 test 1b: file-boundary bracket.
   ! =====================================================================

   subroutine test_file_boundary_bracket(error)
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: VARNAMES(9) = &
                                     [character(len=4) :: "tas", "huss", "uas", "vas", "psl", &
                                                           "rlds", "rsds", "prra", "prsn"]
      real(wp), parameter :: A = 2.0_wp, B = 0.01_wp, C = 0.03_wp, D = 0.7_wp

      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(config_t) :: cfg
      type(ocean_data_input_t) :: reader
      type(atm_state_t) :: atm
      real(wp) :: lon_c(NLON_SRC), lat_c(NLAT_SRC)
      real(wp) :: lon_bnds(2, NLON_SRC), lat_bnds(2, NLAT_SRC)
      real(wp) :: t_axis1(2), t_axis2(2), coefs(4, 9)
      real(wp) :: tq, ll(2), expect
      integer :: v, ierr

      call uniform_axis(NLON_SRC, 0.0_wp, 360.0_wp, lon_c, lon_bnds)
      call uniform_axis(NLAT_SRC, -90.0_wp, 90.0_wp, lat_c, lat_bnds)
      ! File 1: days 0,1 since 2000-01-01. File 2: days 0,1 since
      ! 2000-01-03 == absolute days 2,3 -- a continuous axis straddled
      ! by a query at day 1.5 (bracket n0 in file 1, n1 in file 2).
      t_axis1 = [0.0_wp, 1.0_wp]
      t_axis2 = [0.0_wp, 1.0_wp]
      do v = 1, 9
         coefs(:, v) = [A + real(v, wp), B, C, D]
      end do

      call write_atm_grid_file("tmp_local_artifacts/atm_bound_1.nc", NLON_SRC, NLAT_SRC, &
                               lon_c, lat_c, lon_bnds, lat_bnds, 2, t_axis1, &
                               "days since 2000-01-01 00:00:00", VARNAMES, coefs)
      call write_atm_grid_file("tmp_local_artifacts/atm_bound_2.nc", NLON_SRC, NLAT_SRC, &
                               lon_c, lat_c, lon_bnds, lat_bnds, 2, t_axis2, &
                               "days since 2000-01-03 00:00:00", VARNAMES, coefs, &
                               t_data=[2.0_wp, 3.0_wp])
      call write_filelist_txt("tmp_local_artifacts/atm_bound_filelist.txt", &
                              ["tmp_local_artifacts/atm_bound_1.nc", &
                               "tmp_local_artifacts/atm_bound_2.nc"])

      call make_model(grid, metrics)
      call base_cfg(cfg, "2000-01-01", "tmp_local_artifacts/atm_bound_filelist.txt")

      call reader%init()
      call atm_state_configure(cfg, reader, grid, metrics, atm, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "atm_state_configure failed")
      if (allocated(error)) return
      call map_atm(reader, atm)

      tq = 1.5_wp*86400.0_wp
      call ocean_data_input_update_all(reader, tq)
      call atm_state_update(atm, reader, metrics, tq, TEST_DT)
      call sync_atm(atm)

      ll = query_lon_lat(1, 1, metrics)
      expect = coefs(1, 1) + coefs(2, 1)*ll(1) + coefs(3, 1)*ll(2) + coefs(4, 1)*1.5_wp
      call check(error, abs(atm%tas(1 + NGHOST, 1 + NGHOST) - expect) < TOL_EXACT, &
                 "file-boundary bracket: tas wrong at the cross-file query")
      if (allocated(error)) return

      call unmap_atm(reader, atm)
      call destroy_cartesian_metrics(metrics)
      call atm%destroy()
      call reader%destroy()
   end subroutine test_file_boundary_bracket

   ! =====================================================================
   ! §5 test 1c: RYF wrap bit-identity.
   ! =====================================================================

   subroutine test_ryf_wrap_bitid(error)
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: VARNAMES(9) = &
                                     [character(len=4) :: "tas", "huss", "uas", "vas", "psl", &
                                                           "rlds", "rsds", "prra", "prsn"]
      real(wp), parameter :: A = 3.0_wp, D = 1.3_wp
         !! B = C = 0 -- spatially uniform, isolating the RYF date-wrap
         !! arithmetic (mirrors `test_filelist_ryf_wrap`'s own style).

      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(config_t) :: cfg
      type(ocean_data_input_t) :: reader
      type(atm_state_t) :: atm
      real(wp) :: lon_c(NLON_SRC), lat_c(NLAT_SRC)
      real(wp) :: lon_bnds(2, NLON_SRC), lat_bnds(2, NLAT_SRC)
      real(wp) :: t_axis(4), coefs(4, 9)
      real(wp) :: tq, v1, v2
      integer :: v, ierr

      call uniform_axis(NLON_SRC, 0.0_wp, 360.0_wp, lon_c, lon_bnds)
      call uniform_axis(NLAT_SRC, -90.0_wp, 90.0_wp, lat_c, lat_bnds)
      t_axis = [0.0_wp, 100.0_wp, 200.0_wp, 300.0_wp]
      do v = 1, 9
         coefs(:, v) = [A, 0.0_wp, 0.0_wp, D]
      end do

      call write_atm_grid_file("tmp_local_artifacts/atm_ryf.nc", NLON_SRC, NLAT_SRC, &
                               lon_c, lat_c, lon_bnds, lat_bnds, 4, t_axis, &
                               "days since 1958-01-01 00:00:00", VARNAMES, coefs)
      call write_filelist_txt("tmp_local_artifacts/atm_ryf_filelist.txt", &
                              ["tmp_local_artifacts/atm_ryf.nc"])

      call make_model(grid, metrics)
      call base_cfg(cfg, "1958-01-01", "tmp_local_artifacts/atm_ryf_filelist.txt")
      cfg%ocean%atm_forcing%time_mode = "repeat_year"
      cfg%ocean%atm_forcing%ryf_start_date = "1958-01-01"
      cfg%ocean%atm_forcing%ryf_period_days = 400.0_wp

      call reader%init()
      call atm_state_configure(cfg, reader, grid, metrics, atm, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "atm_state_configure failed")
      if (allocated(error)) return
      call map_atm(reader, atm)

      tq = 350.0_wp*86400.0_wp
      call ocean_data_input_update_all(reader, tq)
      call atm_state_update(atm, reader, metrics, tq, TEST_DT)
      call sync_atm(atm)
      v1 = atm%tas(1 + NGHOST, 1 + NGHOST)

      tq = (350.0_wp + 400.0_wp)*86400.0_wp
      call ocean_data_input_update_all(reader, tq)
      call atm_state_update(atm, reader, metrics, tq, TEST_DT)
      call sync_atm(atm)
      v2 = atm%tas(1 + NGHOST, 1 + NGHOST)

      call check(error, v1 == v2, &
                 "ryf wrap: one period later must reproduce the same wrapped value")
      if (allocated(error)) return

      call unmap_atm(reader, atm)
      call destroy_cartesian_metrics(metrics)
      call atm%destroy()
      call reader%destroy()
   end subroutine test_ryf_wrap_bitid

   ! =====================================================================
   ! §5 test 2: cadence asymmetry -- states vs fluxes on DIFFERENT time
   ! axes, each must recover its own analytic value at a shared query.
   ! =====================================================================

   subroutine test_cadence_asymmetry(error)
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: STATE_VARS(5) = &
                                     [character(len=4) :: "tas", "huss", "uas", "vas", "psl"]
      character(len=*), parameter :: FLUX_VARS(4) = &
                                     [character(len=4) :: "rlds", "rsds", "prra", "prsn"]
      real(wp), parameter :: A_S = 10.0_wp, D_S = 2.0_wp
         !! states: 3hrPt-like stamps at 0, 3, 6, 9, 12 h.
      real(wp), parameter :: A_F = 50.0_wp, D_F = 5.0_wp
         !! fluxes: 3hr-mean-like stamps at 1.5, 4.5, 7.5, 10.5, 13.5 h
         !! -- offset by half the cadence, so the SAME query time falls
         !! in a DIFFERENT bracket with a different blend weight.

      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(config_t) :: cfg
      type(ocean_data_input_t) :: reader
      type(atm_state_t) :: atm
      real(wp) :: lon_c(NLON_SRC), lat_c(NLAT_SRC)
      real(wp) :: lon_bnds(2, NLON_SRC), lat_bnds(2, NLAT_SRC)
      real(wp) :: t_state(5), t_flux(5), coefs_s(4, 5), coefs_f(4, 4)
      real(wp) :: tq_hours, tq, expect
      integer :: v, ierr

      call uniform_axis(NLON_SRC, 0.0_wp, 360.0_wp, lon_c, lon_bnds)
      call uniform_axis(NLAT_SRC, -90.0_wp, 90.0_wp, lat_c, lat_bnds)
      t_state = [0.0_wp, 3.0_wp, 6.0_wp, 9.0_wp, 12.0_wp]
      t_flux = [1.5_wp, 4.5_wp, 7.5_wp, 10.5_wp, 13.5_wp]
      do v = 1, 5
         coefs_s(:, v) = [A_S + real(v, wp), 0.0_wp, 0.0_wp, D_S]
      end do
      do v = 1, 4
         coefs_f(:, v) = [A_F + real(v, wp), 0.0_wp, 0.0_wp, D_F]
      end do

      call write_atm_grid_file("tmp_local_artifacts/atm_cadence_state.nc", NLON_SRC, NLAT_SRC, &
                               lon_c, lat_c, lon_bnds, lat_bnds, 5, t_state, &
                               "hours since 2000-01-01 00:00:00", STATE_VARS, coefs_s)
      call write_atm_grid_file("tmp_local_artifacts/atm_cadence_flux.nc", NLON_SRC, NLAT_SRC, &
                               lon_c, lat_c, lon_bnds, lat_bnds, 5, t_flux, &
                               "hours since 2000-01-01 00:00:00", FLUX_VARS, coefs_f)
      call write_filelist_txt("tmp_local_artifacts/atm_cadence_state_filelist.txt", &
                              ["tmp_local_artifacts/atm_cadence_state.nc"])
      call write_filelist_txt("tmp_local_artifacts/atm_cadence_flux_filelist.txt", &
                              ["tmp_local_artifacts/atm_cadence_flux.nc"])

      call make_model(grid, metrics)
      cfg%start_date = "2000-01-01"
      cfg%ocean%atm_forcing%enable = .true.
      cfg%ocean%atm_forcing%sw_interp = "linear"
      cfg%ocean%atm_forcing%fluxes_regrid = "bilinear"
      cfg%ocean%atm_forcing%time_mode = "interannual"
      cfg%ocean%atm_forcing%filelist_tas = "tmp_local_artifacts/atm_cadence_state_filelist.txt"
      cfg%ocean%atm_forcing%filelist_huss = "tmp_local_artifacts/atm_cadence_state_filelist.txt"
      cfg%ocean%atm_forcing%filelist_uas = "tmp_local_artifacts/atm_cadence_state_filelist.txt"
      cfg%ocean%atm_forcing%filelist_vas = "tmp_local_artifacts/atm_cadence_state_filelist.txt"
      cfg%ocean%atm_forcing%filelist_psl = "tmp_local_artifacts/atm_cadence_state_filelist.txt"
      cfg%ocean%atm_forcing%filelist_rlds = "tmp_local_artifacts/atm_cadence_flux_filelist.txt"
      cfg%ocean%atm_forcing%filelist_rsds = "tmp_local_artifacts/atm_cadence_flux_filelist.txt"
      cfg%ocean%atm_forcing%filelist_prra = "tmp_local_artifacts/atm_cadence_flux_filelist.txt"
      cfg%ocean%atm_forcing%filelist_prsn = "tmp_local_artifacts/atm_cadence_flux_filelist.txt"

      call reader%init()
      call atm_state_configure(cfg, reader, grid, metrics, atm, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "atm_state_configure failed")
      if (allocated(error)) return
      call map_atm(reader, atm)

      tq_hours = 5.0_wp
      tq = tq_hours*3600.0_wp
      call ocean_data_input_update_all(reader, tq)
      call atm_state_update(atm, reader, metrics, tq, TEST_DT)
      call sync_atm(atm)

      ! tas: bracket [3,6]h, weight (5-3)/3 = 0.667 -- exact recovery
      ! proves the STATE reader id/bracket was used.
      expect = coefs_s(1, 1) + coefs_s(4, 1)*tq_hours
      call check(error, abs(atm%tas(1 + NGHOST, 1 + NGHOST) - expect) < TOL_EXACT, &
                 "cadence asymmetry: tas (3hrPt-like axis) wrong")
      if (allocated(error)) return

      ! rlds: bracket [4.5,7.5]h, weight (5-4.5)/3 = 0.167 -- a DIFFERENT
      ! bracket/weight than tas's at the SAME query time; exact recovery
      ! proves rlds used its OWN (flux) reader id, not tas's.
      expect = coefs_f(1, 1) + coefs_f(4, 1)*tq_hours
      call check(error, abs(atm%rlds(1 + NGHOST, 1 + NGHOST) - expect) < TOL_EXACT, &
                 "cadence asymmetry: rlds (3hr-mean-like axis) wrong")
      if (allocated(error)) return

      call unmap_atm(reader, atm)
      call destroy_cartesian_metrics(metrics)
      call atm%destroy()
      call reader%destroy()
   end subroutine test_cadence_asymmetry

   ! =====================================================================
   ! §5 test 7: GPU residency + bytes().
   ! =====================================================================

   subroutine test_gpu_resident(error)
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: VARNAMES(9) = &
                                     [character(len=4) :: "tas", "huss", "uas", "vas", "psl", &
                                                           "rlds", "rsds", "prra", "prsn"]
      real(wp), parameter :: A = 4.0_wp, B = 0.01_wp, C = 0.02_wp, D = 0.6_wp

      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(config_t) :: cfg
      type(ocean_data_input_t) :: reader
      type(atm_state_t) :: atm
      real(wp) :: lon_c(NLON_SRC), lat_c(NLAT_SRC)
      real(wp) :: lon_bnds(2, NLON_SRC), lat_bnds(2, NLAT_SRC)
      real(wp) :: t_axis(3), coefs(4, 9)
      real(wp) :: tq, ll(2), expect
      integer(kind=8) :: measured, reported
      integer :: v, ierr

      call uniform_axis(NLON_SRC, 0.0_wp, 360.0_wp, lon_c, lon_bnds)
      call uniform_axis(NLAT_SRC, -90.0_wp, 90.0_wp, lat_c, lat_bnds)
      t_axis = [0.0_wp, 1.0_wp, 2.0_wp]
      do v = 1, 9
         coefs(:, v) = [A + real(v, wp), B, C, D]
      end do

      call write_atm_grid_file("tmp_local_artifacts/atm_gpu.nc", NLON_SRC, NLAT_SRC, &
                               lon_c, lat_c, lon_bnds, lat_bnds, 3, t_axis, &
                               "days since 2000-01-01 00:00:00", VARNAMES, coefs)
      call write_filelist_txt("tmp_local_artifacts/atm_gpu_filelist.txt", &
                              ["tmp_local_artifacts/atm_gpu.nc"])

      call make_model(grid, metrics)
      call base_cfg(cfg, "2000-01-01", "tmp_local_artifacts/atm_gpu_filelist.txt")

      call reader%init()
      call atm_state_configure(cfg, reader, grid, metrics, atm, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "atm_state_configure failed")
      if (allocated(error)) return

      ! `bytes()` reconciled against the SAME arrays `enter_data` maps --
      ! measured independently here, not just re-derived from `atm`'s own
      ! stored sizes, so a term silently dropped from either side shows up.
      measured = 8_8*(int(size(atm%tas), 8) + int(size(atm%huss), 8) + &
                      int(size(atm%uas), 8) + int(size(atm%vas), 8) + &
                      int(size(atm%psl), 8) + int(size(atm%rlds), 8) + &
                      int(size(atm%prra), 8) + int(size(atm%prsn), 8) + &
                      int(size(atm%rsds), 8) + int(size(atm%rsds_rec), 8) + &
                      int(size(atm%cosz), 8) + int(size(atm%cosz_mean), 8) + &
                      int(size(atm%win_bl), 8) + &
                      int(size(atm%win_cv), 8) + int(size(atm%tgt_bl), 8) + &
                      int(size(atm%tgt_cv), 8)) + atm%wb%bytes() + atm%wc%bytes()
      reported = atm%bytes()
      call check(error, reported == measured, "atm_state%bytes() does not match the mapped arrays")
      if (allocated(error)) return

      !$acc enter data copyin(reader)
      call reader%enter_data()
      !$acc enter data copyin(atm)
      call atm%enter_data()

      tq = 1.37_wp*86400.0_wp
      call ocean_data_input_update_all(reader, tq)
      call atm_state_update(atm, reader, metrics, tq, TEST_DT)

      ! The test's OWN `update self` -- atm_state_update itself issues no
      ! host-visible update, by the mem:separate contract (CLAUDE.md /
      ! `test_open_boundary_out_closes` template).
      !$acc update self(atm%tas, atm%huss, atm%uas, atm%vas, atm%psl)
      !$acc update self(atm%rlds, atm%rsds, atm%prra, atm%prsn)

      ll = query_lon_lat(1, 1, metrics)
      expect = coefs(1, 1) + coefs(2, 1)*ll(1) + coefs(3, 1)*ll(2) + coefs(4, 1)*(tq/86400.0_wp)
      call check(error, abs(atm%tas(1 + NGHOST, 1 + NGHOST) - expect) < TOL_EXACT, &
                 "gpu_resident: tas wrong after the mapped update")
      if (allocated(error)) return

      call atm%exit_data()
      !$acc exit data delete(atm)
      call reader%exit_data()
      !$acc exit data delete(reader)

      call destroy_cartesian_metrics(metrics)
      call atm%destroy()
      call reader%destroy()
   end subroutine test_gpu_resident

   ! =====================================================================
   ! §5 test 8: bit-identity gate -- enable=.false. is a total no-op.
   ! =====================================================================

   subroutine test_disabled_bit_identity(error)
      type(error_type), allocatable, intent(out) :: error

      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(config_t) :: cfg
      type(ocean_data_input_t) :: reader
      type(atm_state_t) :: atm
      integer :: ierr

      call make_model(grid, metrics)
      ! Default `config_t`: `&ocean_atm_forcing_nml enable = .false.`.

      call reader%init()
      call atm_state_configure(cfg, reader, grid, metrics, atm, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "disabled configure must not fail")
      if (allocated(error)) return
      call check(error,.not. atm%enable, "atm%enable must stay .false. when the knob is off")
      if (allocated(error)) return
      call check(error,.not. allocated(atm%tas), &
                 "disabled atm_state must allocate NOTHING")
      if (allocated(error)) return
      call check(error, atm%bytes() == 0_8, "disabled atm_state must report 0 bytes")
      if (allocated(error)) return

      ! A no-op update must not touch (or crash on) the unallocated state.
      call ocean_data_input_update_all(reader, 123.0_wp)
      call atm_state_update(atm, reader, metrics, 123.0_wp, TEST_DT)

      call destroy_cartesian_metrics(metrics)
      call atm%destroy()
      call reader%destroy()
   end subroutine test_disabled_bit_identity

   ! =====================================================================
   ! §5 test 3a: held-record edge -- [00:00,03:00) reads record 1 exactly.
   ! =====================================================================

   subroutine test_sw_held_record_edge(error)
      !! `sw_interp="coszen"` registers `rsds` with `hold=.true.` and
      !! `t_offset=sw_offset_s=+5400`.  The synthetic SW file's records
      !! are spatially uniform (`B=C=0`) and carry the record's OWN raw
      !! stamp time as the value (`D=1`), so a blended (non-held) read
      !! would vary with the query time inside the bracket while a HELD
      !! read stays exactly one record's value for every query in a
      !! bucket.  The records are stamped at the bucket CENTRES (01:30,
      !! 04:30, 07:30 -- JRA55-do), so every query in [00:00,03:00) reads
      !! record 1 (5400.0) and every query in [03:00,06:00) reads record
      !! 2 (16200.0): the record of the bucket that CONTAINS the query
      !! time, switching exactly at the 03:00 boundary.
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: OTHER8(8) = &
                                     [character(len=4) :: "tas", "huss", "uas", "vas", "psl", &
                                                           "rlds", "prra", "prsn"]
      character(len=*), parameter :: RSDS1(1) = [character(len=4) :: "rsds"]
      real(wp), parameter :: D = 1.0_wp

      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(config_t) :: cfg
      type(ocean_data_input_t) :: reader
      type(atm_state_t) :: atm
      real(wp) :: lon_c(NLON_SRC), lat_c(NLAT_SRC)
      real(wp) :: lon_bnds(2, NLON_SRC), lat_bnds(2, NLAT_SRC)
      real(wp) :: t_axis_other(4), t_axis_rsds(3), coefs_other(4, 8), coefs_rsds(4, 1)
      real(wp) :: tq, rsds_rec_val
      integer :: v, ierr, k
      real(wp), parameter :: TQ_SAMPLES(6) = [0.0_wp, 3600.0_wp, 10799.0_wp, &
                                              10800.0_wp, 14400.0_wp, 21599.0_wp]
      real(wp), parameter :: EXPECTED(6) = [5400.0_wp, 5400.0_wp, 5400.0_wp, &
                                            16200.0_wp, 16200.0_wp, 16200.0_wp]

      call uniform_axis(NLON_SRC, 0.0_wp, 360.0_wp, lon_c, lon_bnds)
      call uniform_axis(NLAT_SRC, -90.0_wp, 90.0_wp, lat_c, lat_bnds)

      ! The OTHER eight fields (oor_mode="limit") need a time axis that
      ! comfortably covers every `TQ_SAMPLES` query (0, 3600, 10799) --
      ! unrelated to the `rsds` hold-edge mechanism this test targets.
      t_axis_other = [0.0_wp, 10800.0_wp, 21600.0_wp, 32400.0_wp]
      do v = 1, 8
         coefs_other(:, v) = [0.0_wp, 0.0_wp, 0.0_wp, D]
      end do
      call write_atm_grid_file("tmp_local_artifacts/atm_hold_other.nc", NLON_SRC, NLAT_SRC, &
                               lon_c, lat_c, lon_bnds, lat_bnds, 4, t_axis_other, &
                               "seconds since 2000-01-01 00:00:00", OTHER8, coefs_other)
      call write_filelist_txt("tmp_local_artifacts/atm_hold_other_filelist.txt", &
                              ["tmp_local_artifacts/atm_hold_other.nc"])

      ! `rsds` alone, on the file-native 3-hour-mean axis (01:30, 04:30,
      ! 07:30) -- `hold=.true.` + `t_offset=+5400` is what makes a query
      ! anywhere in a bucket read that bucket's centre-stamped record.
      t_axis_rsds = [5400.0_wp, 16200.0_wp, 27000.0_wp]
      coefs_rsds(:, 1) = [0.0_wp, 0.0_wp, 0.0_wp, D]
      call write_atm_grid_file("tmp_local_artifacts/atm_hold_rsds.nc", NLON_SRC, NLAT_SRC, &
                               lon_c, lat_c, lon_bnds, lat_bnds, 3, t_axis_rsds, &
                               "seconds since 2000-01-01 00:00:00", RSDS1, coefs_rsds)
      call write_filelist_txt("tmp_local_artifacts/atm_hold_rsds_filelist.txt", &
                              ["tmp_local_artifacts/atm_hold_rsds.nc"])

      call make_model(grid, metrics)
      call base_cfg(cfg, "2000-01-01", "tmp_local_artifacts/atm_hold_other_filelist.txt")
      cfg%ocean%atm_forcing%filelist_rsds = "tmp_local_artifacts/atm_hold_rsds_filelist.txt"
      cfg%ocean%atm_forcing%sw_interp = "coszen"
      cfg%ocean%atm_forcing%sw_offset_s = 5400.0_wp

      call reader%init()
      call atm_state_configure(cfg, reader, grid, metrics, atm, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "atm_state_configure failed")
      if (allocated(error)) return
      call map_atm(reader, atm)

      do k = 1, size(TQ_SAMPLES)
         tq = TQ_SAMPLES(k)
         call ocean_data_input_update_all(reader, tq)
         call atm_state_update(atm, reader, metrics, tq, 1.0_wp)   ! 1 s window: every sample time is on its grid
         call sync_atm(atm)
         rsds_rec_val = atm%rsds_rec(1 + NGHOST, 1 + NGHOST)
         call check(error, abs(rsds_rec_val - EXPECTED(k)) < TOL_EXACT, &
                    "held record: rsds_rec must equal the raw value of the record "// &
                    "stamped at the centre of the bucket containing the query time")
         if (allocated(error)) then
            call unmap_atm(reader, atm)
            call destroy_cartesian_metrics(metrics)
            call atm%destroy()
            call reader%destroy()
            return
         end if
      end do

      call unmap_atm(reader, atm)
      call destroy_cartesian_metrics(metrics)
      call atm%destroy()
      call reader%destroy()
   end subroutine test_sw_held_record_edge

   ! =====================================================================
   ! §5 test 3b..: bucket-energy, bound, polar and alignment gates for the
   ! `coszen` diurnal redistribution.  All run the full provider (reader,
   ! hold, regrid, kernel) on ONE physical cell at a chosen lon/lat.
   ! =====================================================================

   subroutine sw_run(lon, lat, start_date, rec, dt, nsteps, rsds_t, cosz_t, mean_t, wlen_t, &
                     used, ierr, bucket, ratio, dfirst, offset)
      !! Run `nsteps` model steps of length `dt` from t = 0 (UTC 00:00 of
      !! `start_date`) on a 1x1 grid centred at (`lon`, `lat`) whose `rsds`
      !! file holds `size(rec)` records stamped at the averaging-interval
      !! CENTRES (`bucket`, default 3 h), record `k` = `rec(k)` uniformly in
      !! space.  `ratio`/`dfirst` emulate `dt_therm_ratio`/`diabatic_first`:
      !! each step's window comes from `atm_sw_window` exactly as the engine
      !! calls it.  Returns the per-step `rsds`, window-mean `cosz`,
      !! interval-mean `cosz_mean`, window length and whether the step is a
      !! thermo step (the only ones whose flux is applied).  Maps every
      !! device array the kernel touches (`map_atm`).
      real(wp), intent(in) :: lon, lat, rec(:), dt
      character(len=*), intent(in) :: start_date
      integer, intent(in) :: nsteps
      real(wp), intent(out) :: rsds_t(nsteps), cosz_t(nsteps), mean_t(nsteps), wlen_t(nsteps)
      logical, intent(out) :: used(nsteps)
      integer, intent(out) :: ierr
      real(wp), intent(in), optional :: bucket, offset
      integer, intent(in), optional :: ratio
      logical, intent(in), optional :: dfirst

      character(len=*), parameter :: OTHER8(8) = &
                                     [character(len=4) :: "tas", "huss", "uas", "vas", "psl", &
                                                           "rlds", "prra", "prsn"]
      character(len=*), parameter :: RSDS1(1) = [character(len=4) :: "rsds"]
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(config_t) :: cfg
      type(ocean_data_input_t) :: reader
      type(atm_state_t) :: atm
      real(wp) :: lon_c(NLON_SRC), lat_c(NLAT_SRC)
      real(wp) :: lon_bnds(2, NLON_SRC), lat_bnds(2, NLAT_SRC)
      real(wp) :: t_axis_other(2), coefs_other(4, 8), coefs_rsds(4, 1), tq
      real(wp) :: bk, shift, wl
      real(wp), allocatable :: t_axis_rsds(:)
      character(len=64) :: units
      integer :: v, k, nrec, nr
      logical :: df, thermo

      bk = BUCKET_S
      if (present(bucket)) bk = bucket
      nr = 1
      if (present(ratio)) nr = ratio
      df = .true.
      if (present(dfirst)) df = dfirst
      nrec = size(rec)
      allocate (t_axis_rsds(nrec))
      call uniform_axis(NLON_SRC, 0.0_wp, 360.0_wp, lon_c, lon_bnds)
      call uniform_axis(NLAT_SRC, -90.0_wp, 90.0_wp, lat_c, lat_bnds)
      units = "seconds since "//start_date//" 00:00:00"
      do k = 1, nrec
         t_axis_rsds(k) = (real(k, wp) - 0.5_wp)*bk
      end do
      do v = 1, 8
         coefs_other(:, v) = [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
      end do
      t_axis_other = [-100000.0_wp, 400000.0_wp]
      call write_atm_grid_file("tmp_local_artifacts/atm_sw_other.nc", NLON_SRC, NLAT_SRC, &
                               lon_c, lat_c, lon_bnds, lat_bnds, 2, t_axis_other, &
                               trim(units), OTHER8, coefs_other)
      call write_filelist_txt("tmp_local_artifacts/atm_sw_other_filelist.txt", &
                              ["tmp_local_artifacts/atm_sw_other.nc"])
      ! value = 1 * t_data(k) = rec(k): piecewise, uniform in space.
      coefs_rsds(:, 1) = [0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp]
      call write_atm_grid_file("tmp_local_artifacts/atm_sw_rsds.nc", NLON_SRC, NLAT_SRC, &
                               lon_c, lat_c, lon_bnds, lat_bnds, nrec, t_axis_rsds, &
                               trim(units), RSDS1, coefs_rsds, t_data=rec)
      call write_filelist_txt("tmp_local_artifacts/atm_sw_rsds_filelist.txt", &
                              ["tmp_local_artifacts/atm_sw_rsds.nc"])

      ! One physical T-cell centred exactly at (lon, lat).
      call grid%init(1, 1, NGHOST, 1000.0_wp, 1000.0_wp)
      call make_spherical_metrics(metrics, grid, lon - 1.0_wp, lat - 1.0_wp, 2.0_wp, 2.0_wp, &
                                  RAD_EARTH)
      call base_cfg(cfg, start_date, "tmp_local_artifacts/atm_sw_other_filelist.txt")
      cfg%ocean%atm_forcing%filelist_rsds = "tmp_local_artifacts/atm_sw_rsds_filelist.txt"
      cfg%ocean%atm_forcing%sw_interp = "coszen"
      cfg%ocean%atm_forcing%sw_offset_s = 0.5_wp*bk
      if (present(offset)) cfg%ocean%atm_forcing%sw_offset_s = offset

      call reader%init()
      call atm_state_configure(cfg, reader, grid, metrics, atm, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) then
         call destroy_cartesian_metrics(metrics)
         call reader%destroy()
         return
      end if
      call map_atm(reader, atm)
      do k = 1, nsteps
         tq = real(k - 1, wp)*dt
         if (nr <= 1) then
            thermo = .true.
         else if (df) then
            thermo = (mod(k - 1, nr) == 0)
         else
            thermo = (mod(k, nr) == 0)
         end if
         call atm_sw_window(dt, nr, df, thermo, shift, wl)
         call ocean_data_input_update_all(reader, tq)
         call atm_state_update(atm, reader, metrics, tq, dt, shift, wl)
         call sync_atm(atm)
         rsds_t(k) = atm%rsds(1 + NGHOST, 1 + NGHOST)
         cosz_t(k) = atm%cosz(1 + NGHOST, 1 + NGHOST)
         mean_t(k) = atm%cosz_mean(1 + NGHOST, 1 + NGHOST)
         wlen_t(k) = wl
         used(k) = thermo
      end do
      call unmap_atm(reader, atm)
      call destroy_cartesian_metrics(metrics)
      call atm%destroy()
      call reader%destroy()
   end subroutine sw_run

   subroutine bucket_audit(rsds_t, mean_t, wlen_t, used, rec, bucket, dt, nsteps, &
                           energy_err, max_ratio, max_frac, n_partial)
      !! Per-bucket energy error `|Sum rsds*wlen - e_b| / (rec_b*bucket)` over
      !! the steps whose flux is applied (`used`; max over buckets;
      !! `e_b = rec_b*bucket` where the bucket's mean cosine is >= 1e-9 and 0
      !! where the Sun never rises in it), the largest `rsds/rec`, the
      !! largest single-window energy fraction `rsds*wlen/(rec*bucket)` (the
      !! bound: <= 1), and the number of PARTLY lit buckets (a
      !! sunrise/sunset bucket: 0 < mean < 0.3 x the day's largest).
      integer, intent(in) :: nsteps
      real(wp), intent(in) :: rsds_t(nsteps), mean_t(nsteps), wlen_t(nsteps), rec(:)
      real(wp), intent(in) :: bucket, dt
      logical, intent(in) :: used(nsteps)
      real(wp), intent(out) :: energy_err, max_ratio, max_frac
      integer, intent(out) :: n_partial

      integer :: k, b, nb
      real(wp) :: expect, mmax
      real(wp) :: e(size(rec)), bmean(size(rec))

      nb = size(rec)
      e = 0.0_wp
      bmean = 0.0_wp
      energy_err = 0.0_wp
      max_ratio = 0.0_wp
      max_frac = 0.0_wp
      n_partial = 0
      do k = 1, nsteps
         b = int(real(k - 1, wp)*dt/bucket) + 1
         if (b > nb) exit
         bmean(b) = max(bmean(b), mean_t(k))
         if (.not. used(k)) cycle
         e(b) = e(b) + rsds_t(k)*wlen_t(k)
         if (rec(b) > 0.0_wp) then
            max_ratio = max(max_ratio, rsds_t(k)/rec(b))
            max_frac = max(max_frac, rsds_t(k)*wlen_t(k)/(rec(b)*bucket))
         end if
      end do
      mmax = maxval(bmean)
      do b = 1, nb
         expect = 0.0_wp
         if (bmean(b) >= 1.0e-9_wp) expect = rec(b)*bucket
         if (rec(b) > 0.0_wp) energy_err = max(energy_err, abs(e(b) - expect)/(rec(b)*bucket))
         if (bmean(b) > 1.0e-9_wp .and. bmean(b) < 0.3_wp*mmax) n_partial = n_partial + 1
      end do
   end subroutine bucket_audit

   subroutine test_sw_energy_identity(error)
      !! One UTC day at a mid-latitude cell (34.7 S, 74.5 W, Jan 3: the
      !! day contains a sunset bucket and a sunrise bucket), uniform
      !! record 250 W/m2: for EVERY 3-hour bucket `Sum rsds*window = rsds_rec*3h`
      !! (to round-off -- the window and bucket means are differences of
      !! one closed-form antiderivative) and no window carries more than the
      !! bucket's energy.  Swept over dt = 600/1800/3600 s, and over
      !! `dt_therm_ratio` = 2 and 3 for BOTH `diabatic_first` settings (the
      !! flux sampled on a thermo step is applied over `N*dt`; only thermo
      !! steps count).
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: REC_VAL = 250.0_wp
      real(wp), parameter :: DTS(6) = [600.0_wp, 1800.0_wp, 3600.0_wp, 900.0_wp, 900.0_wp, 900.0_wp]
      integer, parameter :: RATIOS(6) = [1, 1, 1, 2, 3, 3]
      logical, parameter :: DFIRST(6) = [.true., .true., .true., .true., .true., .false.]
      real(wp) :: rec(8), eerr, mratio, mfrac, dt
      real(wp), allocatable :: r(:), c(:), m(:), w(:)
      logical, allocatable :: u(:)
      integer :: id, ns, ierr, npart, pass, nr
      logical :: df

      rec = REC_VAL
      do pass = 1, 8
         ! passes 1-6: the table; 7: ratio 2 not diabatic_first; 8: ratio 2 N=2 dt 1800
         if (pass <= 6) then
            id = pass
            dt = DTS(id)
            nr = RATIOS(id)
            df = DFIRST(id)
         else if (pass == 7) then
            dt = 900.0_wp
            nr = 2
            df = .false.
         else
            dt = 1800.0_wp
            nr = 2
            df = .true.
         end if
         ns = nint(86400.0_wp/dt)
         allocate (r(ns), c(ns), m(ns), w(ns), u(ns))
         call sw_run(285.5_wp, -34.7_wp, "2000-01-03", rec, dt, ns, r, c, m, w, u, ierr, &
                     ratio=nr, dfirst=df)
         call check(error, ierr == OCEAN_STATUS_OK, "atm_state_configure failed")
         if (allocated(error)) return
         call bucket_audit(r, m, w, u, rec, BUCKET_S, dt, ns, eerr, mratio, mfrac, npart)
         call check(error, eerr <= 1.0e-9_wp, &
                    "energy identity: per-bucket Sum rsds*window must equal rsds_rec*3h "// &
                    "(relative error <= 1e-9) at every dt / dt_therm_ratio / diabatic_first")
         if (allocated(error)) return
         call check(error, mfrac <= 1.0_wp + 1.0e-9_wp, &
                    "bound: one window may never carry more than the bucket's energy")
         if (allocated(error)) return
         call check(error, npart >= 2, &
                    "test geometry must contain both a sunrise and a sunset bucket")
         if (allocated(error)) return
         call check(error, minval(r) >= 0.0_wp, "rsds must be non-negative")
         if (allocated(error)) return
         deallocate (r, c, m, w, u)
      end do
   end subroutine test_sw_energy_identity

   subroutine test_sw_failing_geometry(error)
      !! The cell that blew the 1-degree Southern Ocean up (74.5 W, 34.7 S,
      !! day of year 3), bucket [00:00,03:00) UTC holding its local sunset.
      !! The old kernel divided the INSTANTANEOUS cosine at t = 0 (0.024,
      !! the Sun still up) by a 12-node bucket mean of 1.6e-6 -- `rsds`
      !! = 1.5e4 x the record -- and read the PREVIOUS bucket's record.
      !! Now: first step <= 6 x record, the bucket's sunlit sliver carries
      !! exactly record x 3 h, and nothing after sunset.
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: DT = 1800.0_wp
      real(wp) :: rec(8), r(6), c(6), m(6), w(6), e
      logical :: u(6)
      integer :: ierr

      rec = 0.0_wp
      rec(1) = 249.28_wp
      call sw_run(285.5_wp, -34.7_wp, "1958-01-03", rec, DT, 6, r, c, m, w, u, ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "atm_state_configure failed")
      if (allocated(error)) return
      call check(error, m(1) > 1.0e-9_wp .and. m(1) < 1.0e-2_wp, &
                 "failing geometry: bucket 1 must be a thin sunlit sliver (sunset)")
      if (allocated(error)) return
      call check(error, maxval(r) <= 6.0_wp*rec(1)*(1.0_wp + 1.0e-9_wp), &
                 "failing geometry: rsds must stay <= 6 x the record (old kernel: ~1.5e4 x)")
      if (allocated(error)) return
      e = sum(r)*DT
      call check(error, abs(e - rec(1)*BUCKET_S) <= 1.0e-9_wp*rec(1)*BUCKET_S, &
                 "failing geometry: the bucket's energy must equal record x 3 h")
      if (allocated(error)) return
      call check(error, all(r(4:6) == 0.0_wp), &
                 "failing geometry: exactly zero once the Sun has set for the whole step")
   end subroutine test_sw_failing_geometry

   subroutine test_sw_polar(error)
      !! 80 N: polar night (Jan 3) -> `rsds` and `cosz_mean` EXACTLY 0 even
      !! with a non-zero record; polar day (Jul 1) -> the Sun never sets,
      !! the day's energy is exact per bucket, `rsds` is bounded and
      !! changes smoothly from step to step inside a bucket (no clamp kink).
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: DT = 1800.0_wp
      real(wp) :: rec(8), r(48), c(48), m(48), w(48), eerr, mratio, mfrac, jump
      logical :: u(48)
      integer :: ierr, npart, k

      rec = 300.0_wp
      call sw_run(10.0_wp, 80.0_wp, "2000-01-03", rec, DT, 48, r, c, m, w, u, ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "atm_state_configure failed")
      if (allocated(error)) return
      call check(error, all(r == 0.0_wp) .and. all(m == 0.0_wp) .and. all(c == 0.0_wp), &
                 "polar night: rsds, cosz and the bucket mean cosine must be exactly 0")
      if (allocated(error)) return

      call sw_run(10.0_wp, 80.0_wp, "2000-07-01", rec, DT, 48, r, c, m, w, u, ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "atm_state_configure failed")
      if (allocated(error)) return
      call bucket_audit(r, m, w, u, rec, BUCKET_S, DT, 48, eerr, mratio, mfrac, npart)
      call check(error, eerr <= 1.0e-9_wp, "polar day: per-bucket energy must be exact")
      if (allocated(error)) return
      call check(error, minval(m) > 0.0_wp .and. minval(r) > 0.0_wp, &
                 "polar day: the Sun never sets, rsds must stay positive")
      if (allocated(error)) return
      call check(error, mratio <= 2.0_wp, "polar day: rsds must stay within 2 x the record")
      if (allocated(error)) return
      jump = 0.0_wp
      ! WITHIN a bucket (rsds_rec*A/B_b: the bucket-to-bucket step is the
      ! ratio of two buckets' mean cosines, inherent to a held record).
      do k = 2, 48
         if (mod(k - 1, 6) == 0) cycle
         jump = max(jump, abs(r(k) - r(k - 1))/rec(1))
      end do
      call check(error, jump < 0.15_wp, "polar day: rsds must vary smoothly step to step")
   end subroutine test_sw_polar

   subroutine test_sw_record_alignment(error)
      !! Record 3 (the [06:00,09:00) UTC mean) is non-zero, every other
      !! record 0, at a cell (60 E, 20 N, Jan 3) whose Sun is up through
      !! the whole 03-12 UTC window: ALL of the energy `rec*3h` must land in
      !! bucket 3 and none in the others.  (A one-bucket-late record -- the
      !! old `sw_offset_s = -5400` -- put it in bucket 4.)
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: DT = 1800.0_wp
      real(wp), parameter :: X = 400.0_wp
      real(wp) :: rec(8), r(48), c(48), m(48), w(48), eb(8)
      logical :: u(48)
      integer :: ierr, b, k

      rec = 0.0_wp
      rec(3) = X
      call sw_run(60.0_wp, 20.0_wp, "2000-01-03", rec, DT, 48, r, c, m, w, u, ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "atm_state_configure failed")
      if (allocated(error)) return
      call check(error, m(13) > 0.1_wp .and. m(19) > 0.1_wp .and. m(7) > 0.1_wp, &
                 "alignment test geometry: the Sun must be up in buckets 2, 3 and 4")
      if (allocated(error)) return
      eb = 0.0_wp
      do b = 1, 8
         do k = (b - 1)*6 + 1, b*6
            eb(b) = eb(b) + r(k)*DT
         end do
      end do
      call check(error, abs(eb(3) - X*BUCKET_S) <= 1.0e-9_wp*X*BUCKET_S, &
                 "alignment: the whole record energy must land in its own bucket")
      if (allocated(error)) return
      call check(error, maxval(abs(eb([1, 2, 4, 5, 6, 7, 8]))) == 0.0_wp, &
                 "alignment: no energy may land in any other bucket")
   end subroutine test_sw_record_alignment

   subroutine test_sw_axis_checks(error)
      !! Configure-time checks against the `rsds` file: (a) a namelist still
      !! carrying `sw_offset_s = -5400` is REFUSED (it would silently
      !! reproduce the one-bucket-late record); (b) the averaging interval
      !! is derived from the file -- 6-hourly records with `sw_offset_s =
      !! 10800` run and conserve energy per 6 h interval; (c) 6-hourly
      !! records with the 3-hourly offset 5400 are refused.
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: DT = 1800.0_wp
      real(wp) :: rec(8), r(48), c(48), m(48), w(48), eerr, mratio, mfrac
      logical :: u(48)
      integer :: ierr, npart

      rec = 250.0_wp
      call sw_run(285.5_wp, -34.7_wp, "2000-01-03", rec, DT, 48, r, c, m, w, u, ierr, &
                  offset=-5400.0_wp)
      call check(error, ierr /= OCEAN_STATUS_OK, &
                 "sw_offset_s = -5400 with centre-stamped 3-hourly records must be refused")
      if (allocated(error)) return

      call sw_run(285.5_wp, -34.7_wp, "2000-01-03", rec, DT, 48, r, c, m, w, u, ierr, &
                  bucket=21600.0_wp, offset=5400.0_wp)
      call check(error, ierr /= OCEAN_STATUS_OK, &
                 "6-hourly records with sw_offset_s = 5400 must be refused")
      if (allocated(error)) return

      rec = 250.0_wp
      call sw_run(285.5_wp, -34.7_wp, "2000-01-03", rec, DT, 48, r, c, m, w, u, ierr, &
                  bucket=21600.0_wp)
      call check(error, ierr == OCEAN_STATUS_OK, "6-hourly records with offset 10800 must run")
      if (allocated(error)) return
      call bucket_audit(r, m, w, u, rec(1:4), 21600.0_wp, DT, 48, eerr, mratio, mfrac, npart)
      call check(error, eerr <= 1.0e-9_wp .and. mfrac <= 1.0_wp + 1.0e-9_wp, &
                 "6-hourly records: per-interval energy must be exact and bounded")
   end subroutine test_sw_axis_checks

   pure function ref_step_mean(lon, lat, doy, t0, t1, bucket) result(cz)
      !! Independent reference: mean of max(0, sin(phi) sin(dec) + cos(phi)
      !! cos(dec) cos(omega)) over [t0, t1) seconds after UTC midnight by
      !! 10^4-node midpoint quadrature.  Spencer (1971) declination and
      !! equation of time at the midpoint of the `bucket` holding `t0`
      !! (written out here, not taken from the kernel); omega from the UTC
      !! hour, the equation of time and the longitude.
      real(wp), intent(in) :: lon, lat, t0, t1, bucket
      integer, intent(in) :: doy
      real(wp) :: cz

      real(wp), parameter :: PI = 3.14159265358979323846_wp
      integer, parameter :: NQ = 10000
      real(wp) :: tmid, hf, g, dec, eq, tau, om, c, acc
      integer :: i

      tmid = (floor(t0/bucket) + 0.5_wp)*bucket
      hf = tmid/3600.0_wp
      g = 2.0_wp*PI*(real(doy - 1, wp) + (hf - 12.0_wp)/24.0_wp)/365.0_wp
      dec = 0.006918_wp - 0.399912_wp*cos(g) + 0.070257_wp*sin(g) - 0.006758_wp*cos(2*g) &
            + 0.000907_wp*sin(2*g) - 0.002697_wp*cos(3*g) + 0.001480_wp*sin(3*g)
      eq = 229.18_wp*(0.000075_wp + 0.001868_wp*cos(g) - 0.032077_wp*sin(g) &
                      - 0.014615_wp*cos(2*g) - 0.040849_wp*sin(2*g))
      acc = 0.0_wp
      do i = 1, NQ
         tau = t0 + (real(i, wp) - 0.5_wp)/real(NQ, wp)*(t1 - t0)
         om = (15.0_wp*tau/3600.0_wp + eq/4.0_wp + lon - 180.0_wp)*PI/180.0_wp
         c = sin(lat*PI/180.0_wp)*sin(dec) + cos(lat*PI/180.0_wp)*cos(dec)*cos(om)
         acc = acc + max(0.0_wp, c)
      end do
      cz = acc/real(NQ, wp)
   end function ref_step_mean

   subroutine test_sw_shape(error)
      !! The energy identity telescopes for ANY monotone antiderivative, so
      !! it cannot see the SHAPE.  Compare the kernel's step-mean `cosz`
      !! (48 half-hour steps of a day) with `ref_step_mean` at: a mid-latitude
      !! cell with sunset + sunrise steps, longitude -280 (an unwrapped
      !! longitude), 65 N in late June (a short night: a bucket holding both
      !! a sunset and the next sunrise), polar night, polar day.  Every
      !! step the reference calls dark must be EXACTLY 0 (covers a night
      !! step across local midnight).
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: DT = 1800.0_wp
      real(wp), parameter :: LONS(5) = [285.5_wp, -280.0_wp, 337.5_wp, 10.0_wp, 10.0_wp]
      real(wp), parameter :: LATS(5) = [-34.7_wp, 40.0_wp, 65.0_wp, 80.0_wp, 80.0_wp]
      integer, parameter :: DOYS(5) = [3, 3, 177, 3, 183]
      character(len=10), parameter :: DATES(5) = [character(len=10) :: "2000-01-03", &
                                                                                         "2000-01-03", "2000-06-25", "2000-01-03", &
                                                                                                                   "2000-07-01"]
      real(wp) :: rec(8), r(48), c(48), m(48), w(48), ref, worst, t0
      logical :: u(48)
      integer :: ic, k, ierr, ndark_ok

      rec = 250.0_wp
      do ic = 1, 5
         call sw_run(LONS(ic), LATS(ic), DATES(ic), rec, DT, 48, r, c, m, w, u, ierr)
         call check(error, ierr == OCEAN_STATUS_OK, "atm_state_configure failed")
         if (allocated(error)) return
         worst = 0.0_wp
         ndark_ok = 0
         do k = 1, 48
            t0 = real(k - 1, wp)*DT
            ref = ref_step_mean(LONS(ic), LATS(ic), DOYS(ic), t0, t0 + DT, BUCKET_S)
            worst = max(worst, abs(c(k) - ref))
            if (ref == 0.0_wp) then
               call check(error, c(k) == 0.0_wp, &
                          "shape: a step the reference calls dark must be exactly 0")
               if (allocated(error)) return
               ndark_ok = ndark_ok + 1
            end if
         end do
         call check(error, worst <= 2.0e-6_wp, &
                    "shape: step-mean cosz must match the independent 1e4-node quadrature")
         if (allocated(error)) return
         if (ic == 4) then
            call check(error, ndark_ok == 48, "polar night: every step dark")
            if (allocated(error)) return
         end if
         if (ic == 3) then
            call check(error, ndark_ok >= 1 .and. ndark_ok < 12, &
                       "65 N late June: a short night must be present")
            if (allocated(error)) return
         end if
      end do
   end subroutine test_sw_shape

   ! =====================================================================
   ! §5 test 3d: coszen vs the analytic equinox/noon/midnight case.
   ! =====================================================================

   subroutine test_sw_coszen_analytic(error)
      !! Equinox (2000-03-20, day-of-year 80), one physical T-cell
      !! centred exactly at (lon, lat) = (0, 0): local noon (UTC 12:00,
      !! lon=0) must give `cosz` within the Spencer fit's own 1e-3 of the
      !! ideal overhead-sun value 1.0 (`cosz` is the mean over the 60 s step
      !! window `[t, t + 60 s)`); local midnight (UTC 00:00) must
      !! give `cosz = 0` (clamped -- the raw cosine is negative).
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: VARNAMES(9) = &
                                     [character(len=4) :: "tas", "huss", "uas", "vas", "psl", &
                                                           "rlds", "rsds", "prra", "prsn"]
      integer, parameter :: NXP1 = 1, NYP1 = 1

      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(config_t) :: cfg
      type(ocean_data_input_t) :: reader
      type(atm_state_t) :: atm
      real(wp) :: lon_c(NLON_SRC), lat_c(NLAT_SRC)
      real(wp) :: lon_bnds(2, NLON_SRC), lat_bnds(2, NLAT_SRC)
      real(wp) :: t_axis(3), coefs(4, 9)
      real(wp) :: tq, cosz_val
      integer :: v, ierr

      call uniform_axis(NLON_SRC, 0.0_wp, 360.0_wp, lon_c, lon_bnds)
      call uniform_axis(NLAT_SRC, -90.0_wp, 90.0_wp, lat_c, lat_bnds)
      t_axis = [-0.5_wp, 0.5_wp, 1.5_wp]   ! daily means stamped at the interval centre
      do v = 1, 9
         coefs(:, v) = [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
      end do

      call write_atm_grid_file("tmp_local_artifacts/atm_coszen.nc", NLON_SRC, NLAT_SRC, &
                               lon_c, lat_c, lon_bnds, lat_bnds, 3, t_axis, &
                               "days since 2000-03-20 00:00:00", VARNAMES, coefs)
      call write_filelist_txt("tmp_local_artifacts/atm_coszen_filelist.txt", &
                              ["tmp_local_artifacts/atm_coszen.nc"])

      ! One physical T-cell centred exactly at (0 deg lon, 0 deg lat):
      ! lon_west + 0.5*dlon = 0, lat_south + 0.5*dlat = 0.
      call grid%init(NXP1, NYP1, NGHOST, 1000.0_wp, 1000.0_wp)
      call make_spherical_metrics(metrics, grid, -1.0_wp, -1.0_wp, 2.0_wp, 2.0_wp, RAD_EARTH)
      call base_cfg(cfg, "2000-03-20", "tmp_local_artifacts/atm_coszen_filelist.txt")
      cfg%ocean%atm_forcing%sw_interp = "coszen"
      cfg%ocean%atm_forcing%sw_offset_s = 43200.0_wp   ! half the daily interval

      call reader%init()
      call atm_state_configure(cfg, reader, grid, metrics, atm, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "atm_state_configure failed")
      if (allocated(error)) return
      call map_atm(reader, atm)

      ! Local noon (UTC 12:00, lon=0).
      tq = 12.0_wp*3600.0_wp
      call ocean_data_input_update_all(reader, tq)
      call atm_state_update(atm, reader, metrics, tq, 60.0_wp)
      call sync_atm(atm)
      cosz_val = atm%cosz(1 + NGHOST, 1 + NGHOST)
      call check(error, abs(cosz_val - 1.0_wp) < 1.0e-3_wp, &
                 "coszen analytic: equinox local noon at the equator must give "// &
                 "cosz within 1e-3 of 1.0")
      if (allocated(error)) then
         call unmap_atm(reader, atm)
         call destroy_cartesian_metrics(metrics)
         call atm%destroy()
         call reader%destroy()
         return
      end if

      ! Local midnight (UTC 00:00, lon=0).
      tq = 0.0_wp
      call ocean_data_input_update_all(reader, tq)
      call atm_state_update(atm, reader, metrics, tq, 60.0_wp)
      call sync_atm(atm)
      cosz_val = atm%cosz(1 + NGHOST, 1 + NGHOST)
      call check(error, abs(cosz_val) < 1.0e-12_wp, &
                 "coszen analytic: equinox local midnight at the equator must give "// &
                 "cosz = 0 (clamped)")
      if (allocated(error)) then
         call unmap_atm(reader, atm)
         call destroy_cartesian_metrics(metrics)
         call atm%destroy()
         call reader%destroy()
         return
      end if

      call unmap_atm(reader, atm)
      call destroy_cartesian_metrics(metrics)
      call atm%destroy()
      call reader%destroy()
   end subroutine test_sw_coszen_analytic

end module test_ocean_atm_state
