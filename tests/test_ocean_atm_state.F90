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
!!                                       `sw_offset_s=-5400` -- every
!!                                       query in [00:00,03:00) reads
!!                                       `rsds_rec` as EXACTLY record 1
!!                                       (held, not blended toward
!!                                       record 2).
!!   atm_sw_energy_identity           -- §5 test 3b: one day of 15-min
!!                                       steps (the kernel's own 12-point
!!                                       quadrature nodes) -- `Sum
!!                                       rsds*dt` over the day equals
!!                                       `Sum rsds_rec*3h` over the 8
!!                                       buckets to <= 1e-10 relative.
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
   use rdb_atm_state, only: atm_state_t, atm_state_configure, atm_state_update
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
      call atm_state_update(atm, reader, metrics, tq)
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
      !$omp target enter data map(to: reader)
      call reader%enter_data()
      !$omp target enter data map(to: atm)
      call atm%enter_data()
   end subroutine map_atm

   subroutine sync_atm(atm)
      !! Bring the nine regridded fields back to the host for the asserts
      !! (component arrays only, never the aggregate).
      type(atm_state_t), intent(inout) :: atm
      !$omp target update from(atm%tas, atm%huss, atm%uas, atm%vas, atm%psl)
      !$omp target update from(atm%rlds, atm%rsds, atm%prra, atm%prsn)
      !$omp target update from(atm%rsds_rec, atm%cosz, atm%cosz_mean)
   end subroutine sync_atm

   subroutine unmap_atm(reader, atm)
      !! Inverse of `map_atm`.
      type(ocean_data_input_t), intent(inout) :: reader
      type(atm_state_t), intent(inout) :: atm
      call atm%exit_data()
      !$omp target exit data map(delete: atm)
      call reader%exit_data()
      !$omp target exit data map(delete: reader)
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
      call atm_state_update(atm, reader, metrics, tq)
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
      call atm_state_update(atm, reader, metrics, tq)
      call sync_atm(atm)
      v1 = atm%tas(1 + NGHOST, 1 + NGHOST)

      tq = (350.0_wp + 400.0_wp)*86400.0_wp
      call ocean_data_input_update_all(reader, tq)
      call atm_state_update(atm, reader, metrics, tq)
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
      call atm_state_update(atm, reader, metrics, tq)
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

      !$omp target enter data map(to: reader)
      call reader%enter_data()
      !$omp target enter data map(to: atm)
      call atm%enter_data()

      tq = 1.37_wp*86400.0_wp
      call ocean_data_input_update_all(reader, tq)
      call atm_state_update(atm, reader, metrics, tq)

      ! The test's OWN `update self` -- atm_state_update itself issues no
      ! host-visible update, by the mem:separate contract (CLAUDE.md /
      ! `test_open_boundary_out_closes` template).
      !$omp target update from(atm%tas, atm%huss, atm%uas, atm%vas, atm%psl)
      !$omp target update from(atm%rlds, atm%rsds, atm%prra, atm%prsn)

      ll = query_lon_lat(1, 1, metrics)
      expect = coefs(1, 1) + coefs(2, 1)*ll(1) + coefs(3, 1)*ll(2) + coefs(4, 1)*(tq/86400.0_wp)
      call check(error, abs(atm%tas(1 + NGHOST, 1 + NGHOST) - expect) < TOL_EXACT, &
                 "gpu_resident: tas wrong after the mapped update")
      if (allocated(error)) return

      call atm%exit_data()
      !$omp target exit data map(delete: atm)
      call reader%exit_data()
      !$omp target exit data map(delete: reader)

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
      call atm_state_update(atm, reader, metrics, 123.0_wp)

      call destroy_cartesian_metrics(metrics)
      call atm%destroy()
      call reader%destroy()
   end subroutine test_disabled_bit_identity

   ! =====================================================================
   ! §5 test 3a: held-record edge -- [00:00,03:00) reads record 1 exactly.
   ! =====================================================================

   subroutine test_sw_held_record_edge(error)
      !! `sw_interp="coszen"` registers `rsds` with `hold=.true.` and
      !! `t_offset=sw_offset_s=-5400`.  The synthetic SW file's records
      !! are spatially uniform (`B=C=0`) and carry the record's OWN raw
      !! stamp time as the value (`D=1`), so a blended (non-held) read
      !! would vary with the query time inside the bracket while a HELD
      !! read stays exactly the record-1 value (5400.0) for every query
      !! across [00:00,03:00) -- including right at the 03:00 boundary
      !! (the half-open bracket convention; see the module docstring's
      !! `rsds` paragraph for why record 1 -- not a "bucket 1" value --
      !! is what the oor="extend" edge clamps to at run start).
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: OTHER8(8) = &
                                     [character(len=4) :: "tas", "huss", "uas", "vas", "psl", &
                                                           "rlds", "prra", "prsn"]
      character(len=*), parameter :: RSDS1(1) = [character(len=4) :: "rsds"]
      real(wp), parameter :: D = 1.0_wp
      real(wp), parameter :: RECORD1_VALUE = 5400.0_wp

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
      real(wp), parameter :: TQ_SAMPLES(3) = [0.0_wp, 3600.0_wp, 10799.0_wp]

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
      ! 07:30) -- `hold=.true.` + `t_offset=-5400` is what makes a query
      ! anywhere in [00:00,03:00) read record 1's raw value, via the
      ! oor_mode="extend" run-start edge (see the module docstring).
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
      cfg%ocean%atm_forcing%sw_offset_s = -5400.0_wp

      call reader%init()
      call atm_state_configure(cfg, reader, grid, metrics, atm, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "atm_state_configure failed")
      if (allocated(error)) return
      call map_atm(reader, atm)

      do k = 1, 3
         tq = TQ_SAMPLES(k)
         call ocean_data_input_update_all(reader, tq)
         call atm_state_update(atm, reader, metrics, tq)
         call sync_atm(atm)
         rsds_rec_val = atm%rsds_rec(1 + NGHOST, 1 + NGHOST)
         call check(error, abs(rsds_rec_val - RECORD1_VALUE) < TOL_EXACT, &
                    "held record: rsds_rec must equal record 1's raw value across "// &
                    "[00:00,03:00)")
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
   ! §5 test 3b: energy identity over one day of 15-min (quadrature-node)
   ! steps.
   ! =====================================================================

   subroutine test_sw_energy_identity(error)
      !! One day (8 three-hour buckets), stepped at the kernel's own
      !! 12-point midpoint-quadrature cadence (15 min) so each step's
      !! instantaneous `cosz` is bit-identical to the matching term of
      !! that same step's `cosz_mean` quadrature -- making `Sum rsds*dt`
      !! over the day equal `Sum rsds_rec*3h` over the 8 buckets by
      !! construction (the module docstring's `rsds` paragraph), to
      !! round-off rather than by any physical argument. Model epoch
      !! chosen off an equinox so neither side is trivially zero.
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: OTHER8(8) = &
                                     [character(len=4) :: "tas", "huss", "uas", "vas", "psl", &
                                                           "rlds", "prra", "prsn"]
      character(len=*), parameter :: RSDS1(1) = [character(len=4) :: "rsds"]
      real(wp), parameter :: A = 200.0_wp, D = 0.001_wp
      real(wp), parameter :: BUCKET_S = 10800.0_wp, DT_STEP = 900.0_wp
      integer, parameter :: N_BUCKETS = 8, N_STEPS = 96

      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(config_t) :: cfg
      type(ocean_data_input_t) :: reader
      type(atm_state_t) :: atm
      real(wp) :: lon_c(NLON_SRC), lat_c(NLAT_SRC)
      real(wp) :: lon_bnds(2, NLON_SRC), lat_bnds(2, NLAT_SRC)
      real(wp) :: t_axis_other(2), t_axis_rsds(N_BUCKETS), coefs_other(4, 8), coefs_rsds(4, 1)
      real(wp) :: tq, lhs, rhs, rsds_val, rsds_rec_val
      integer :: v, ierr, k, b, nsub

      call uniform_axis(NLON_SRC, 0.0_wp, 360.0_wp, lon_c, lon_bnds)
      call uniform_axis(NLAT_SRC, -90.0_wp, 90.0_wp, lat_c, lat_bnds)
      do b = 1, N_BUCKETS
         t_axis_rsds(b) = (real(b, wp) - 0.5_wp)*BUCKET_S
      end do
      do v = 1, 8
         coefs_other(:, v) = [A + real(v, wp), 0.0_wp, 0.0_wp, D]
      end do
      coefs_rsds(:, 1) = [A + 9.0_wp, 0.0_wp, 0.0_wp, D]

      ! OTHER eight fields (oor_mode="limit"): a 2-point axis comfortably
      ! wide enough for every query this test makes -- their actual
      ! values are irrelevant here, only `rsds`/`rsds_rec` are checked.
      t_axis_other = [-100000.0_wp, 200000.0_wp]
      call write_atm_grid_file("tmp_local_artifacts/atm_energy_other.nc", NLON_SRC, NLAT_SRC, &
                               lon_c, lat_c, lon_bnds, lat_bnds, 2, t_axis_other, &
                               "seconds since 2000-06-21 00:00:00", OTHER8, coefs_other)
      call write_filelist_txt("tmp_local_artifacts/atm_energy_other_filelist.txt", &
                              ["tmp_local_artifacts/atm_energy_other.nc"])

      ! `rsds` alone, on the file-native 3-hour-mean axis spanning the day.
      call write_atm_grid_file("tmp_local_artifacts/atm_energy_rsds.nc", NLON_SRC, NLAT_SRC, &
                               lon_c, lat_c, lon_bnds, lat_bnds, N_BUCKETS, t_axis_rsds, &
                               "seconds since 2000-06-21 00:00:00", RSDS1, coefs_rsds)
      call write_filelist_txt("tmp_local_artifacts/atm_energy_rsds_filelist.txt", &
                              ["tmp_local_artifacts/atm_energy_rsds.nc"])

      call make_model(grid, metrics)
      call base_cfg(cfg, "2000-06-21", "tmp_local_artifacts/atm_energy_other_filelist.txt")
      cfg%ocean%atm_forcing%filelist_rsds = "tmp_local_artifacts/atm_energy_rsds_filelist.txt"
      cfg%ocean%atm_forcing%sw_interp = "coszen"
      cfg%ocean%atm_forcing%sw_offset_s = -5400.0_wp

      call reader%init()
      call atm_state_configure(cfg, reader, grid, metrics, atm, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "atm_state_configure failed")
      if (allocated(error)) return
      call map_atm(reader, atm)

      ! `nsub` = fine steps per bucket (12, the kernel's own quadrature
      ! node count) -- a single pass accumulates LHS every step and RHS
      ! once per bucket (its FIRST sub-step's `rsds_rec`/`cosz_mean`,
      ! both constant across the whole bucket under `hold`).  This
      ! sidesteps having to independently re-derive which file RECORD a
      ! given bucket's `rsds_rec` reads (see
      ! `test_sw_held_record_edge`'s docstring for why that mapping is
      ! offset by one bucket after the very first) -- the identity only
      ! needs `rsds_rec`/`cosz_mean` constant per bucket, guaranteed by
      ! construction regardless of the record index. The RHS mirrors
      ! the production "0 where cosz_mean < 1e-6" rule (a bucket that is
      ! entirely night contributes 0 to BOTH sides, not `rsds_rec*3h` --
      ! this synthetic `rsds_rec` is an arbitrary affine-in-time number,
      ! not a real dataset's already-near-zero nighttime mean).
      nsub = N_STEPS/N_BUCKETS
      lhs = 0.0_wp
      rhs = 0.0_wp
      do k = 1, N_STEPS
         tq = (real(k, wp) - 0.5_wp)*DT_STEP
         call ocean_data_input_update_all(reader, tq)
         call atm_state_update(atm, reader, metrics, tq)
         call sync_atm(atm)
         rsds_val = atm%rsds(1 + NGHOST, 1 + NGHOST)
         lhs = lhs + rsds_val*DT_STEP
         if (mod(k - 1, nsub) == 0) then
            if (atm%cosz_mean(1 + NGHOST, 1 + NGHOST) >= 1.0e-6_wp) then
               rsds_rec_val = atm%rsds_rec(1 + NGHOST, 1 + NGHOST)
               rhs = rhs + rsds_rec_val*BUCKET_S
            end if
         end if
      end do

      call check(error, abs(lhs - rhs) <= 1.0e-10_wp*abs(rhs), &
                 "energy identity: Sum rsds*dt over the day must equal "// &
                 "Sum rsds_rec*3h over the buckets to <= 1e-10 relative")
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
   end subroutine test_sw_energy_identity

   ! =====================================================================
   ! §5 test 3d: coszen vs the analytic equinox/noon/midnight case.
   ! =====================================================================

   subroutine test_sw_coszen_analytic(error)
      !! Equinox (2000-03-20, day-of-year 80), one physical T-cell
      !! centred exactly at (lon, lat) = (0, 0): local noon (UTC 12:00,
      !! lon=0) must give `cosz` within the Spencer fit's own 1e-3 of the
      !! ideal overhead-sun value 1.0; local midnight (UTC 00:00) must
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
      real(wp) :: t_axis(2), coefs(4, 9)
      real(wp) :: tq, cosz_val
      integer :: v, ierr

      call uniform_axis(NLON_SRC, 0.0_wp, 360.0_wp, lon_c, lon_bnds)
      call uniform_axis(NLAT_SRC, -90.0_wp, 90.0_wp, lat_c, lat_bnds)
      t_axis = [0.0_wp, 1.0_wp]
      do v = 1, 9
         coefs(:, v) = [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
      end do

      call write_atm_grid_file("tmp_local_artifacts/atm_coszen.nc", NLON_SRC, NLAT_SRC, &
                               lon_c, lat_c, lon_bnds, lat_bnds, 2, t_axis, &
                               "days since 2000-03-20 00:00:00", VARNAMES, coefs)
      call write_filelist_txt("tmp_local_artifacts/atm_coszen_filelist.txt", &
                              ["tmp_local_artifacts/atm_coszen.nc"])

      ! One physical T-cell centred exactly at (0 deg lon, 0 deg lat):
      ! lon_west + 0.5*dlon = 0, lat_south + 0.5*dlat = 0.
      call grid%init(NXP1, NYP1, NGHOST, 1000.0_wp, 1000.0_wp)
      call make_spherical_metrics(metrics, grid, -1.0_wp, -1.0_wp, 2.0_wp, 2.0_wp, RAD_EARTH)
      call base_cfg(cfg, "2000-03-20", "tmp_local_artifacts/atm_coszen_filelist.txt")
      cfg%ocean%atm_forcing%sw_interp = "coszen"

      call reader%init()
      call atm_state_configure(cfg, reader, grid, metrics, atm, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "atm_state_configure failed")
      if (allocated(error)) return
      call map_atm(reader, atm)

      ! Local noon (UTC 12:00, lon=0).
      tq = 12.0_wp*3600.0_wp
      call ocean_data_input_update_all(reader, tq)
      call atm_state_update(atm, reader, metrics, tq)
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
      call atm_state_update(atm, reader, metrics, tq)
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
