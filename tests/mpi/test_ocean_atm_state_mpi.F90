!! MPI bit-identity gate for the JRA55-do-style atmospheric-state provider
!! (`rdb_atm_state`) + the ocean bulk-flux consumer (`rdb_ocean_bulk_forcing`),
!! OM3 forcing wave-3 C2c-3b -- PLAN_C2c.md §5 test 6.
!!
!! Same bit-identity contract as `test_ocean_decomp_bitid_mpi`: a decomposed
!! run IS the serial run, cell-for-cell, bitwise -- but scoped to the three
!! fields this chunk's brief names (`atm%tas`, `ss%tau_x/tau_y`,
!! `sf%q_sens`) through the LOWER-LEVEL building blocks
!! (`atm_state_configure`/`atm_state_update` + hand-built `bulk_forcing_t`/
!! `multilayer_state_t`/`ocean_surface_stress_t`/`ocean_surface_flux_t`,
!! exactly as `tests/test_ocean_atm_state.F90` and
!! `tests/test_ocean_bulk_forcing.F90` already exercise them at single rank)
!! rather than the full production engine -- no bathymetry/EOS/time
!! integration is needed to reach any of the three decomposition-sensitive
!! code paths this test exists for:
!!   1. `atm_state_configure`'s per-rank regrid-weight build
!!      (`regrid_weights_bilinear`/`_conservative`), which depends on THIS
!!      rank's own `metrics%geolonT/geolatT` range;
!!   2. `atm_state_update`'s per-rank windowed file read (the reader's
!!      `win_i0/win_j0/win_nx/win_ny` override, C2c-1) + blend + regrid +
!!      scatter into the physical interior;
!!   3. `bulk_forcing_apply`'s cell-centred stress -> scalar halo -> face
!!      average -> REAL vector seam refresh (`ocean_seam_refresh_surface_
!!      stress`) -- the one path here that genuinely crosses an MPI rank
!!      seam (interior tile-tile faces), as opposed to 1-2 which are
!!      per-rank-local once each rank's own window is correct.
!!
!! On EACH launch (ctest runs np = 1, 4):
!!   1. every physical rank independently builds the SERIAL reference
!!      (px = py = 1, the whole 12x12 global sector on one logical tile --
!!      no cross-rank communication, same "run compute_size=1 on every
!!      rank" trick `test_ocean_decomp_bitid_mpi::run_one` uses);
!!   2. at nprocs = 4 only, every physical rank ALSO builds the real
!!      2x2-decomposed run (px = py = 2, crank = this rank) and the four
!!      ranks genuinely exchange halos;
!!   3. the decomposed run's owned window (ghosts excluded; the two face
!!      fields' owned window includes both of the tile's edge faces, same
!!      convention as `test_ocean_decomp_bitid_mpi::compare`) is compared,
!!      BITWISE, against the matching slice of the reference.
!! At nprocs = 1 there is nothing to decompose -- the run is a smoke test
!! (serial reference built and sane) matching the `_1rank` leg's role for
!! every other MPI test in this suite that needs >= 2 ranks for real teeth.
!!
!! Spatially-varying synthetic JRA-like fields (affine in lon/lat/time,
!! `tests/test_ocean_atm_state.F90`'s own convention) are written to ONE
!! shared NetCDF file (all nine CMOR variable names) by rank 0 before either
!! run, so `uas`/`vas`/`tas`/... -- and hence the bulk-flux outputs -- vary
!! across the sector and across the rank seam, giving the comparison real
!! teeth (a spatially UNIFORM field, as the single-rank unit tests use,
!! would pass even with a broken per-rank window).
#ifdef RDB_ENABLE_MPI
program test_ocean_atm_state_mpi
   use, intrinsic :: iso_fortran_env, only: int64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_decomp, only: decomp_t, decomp_init
   use rdb_comm_env, only: comm_env_init, comm_env_setup_roles, comm_env_finalize, &
                           comm_env_rank, comm_env_size, comm_env_compute_comm
   use pic_mpi_lib, only: comm_t, allreduce, MPI_SUM
   use rdb_ocean_halo, only: ocean_halo_init, ocean_halo_destroy
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_spherical_metrics, destroy_cartesian_metrics
   use rdb_config, only: config_t
   use rdb_ocean_data_input, only: ocean_data_input_t, ocean_data_input_update_all
   use rdb_atm_state, only: atm_state_t, atm_state_configure, atm_state_update
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   use rdb_ocean_boundary_types, only: ocean_bc_state_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ice_state, only: ocean_sea_ice_t
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t
   use rdb_ocean_bulk_forcing, only: bulk_forcing_t, bulk_forcing_apply
   use rdb_io_netcdf, only: nc_create_file, nc_close, nc_def_dim, nc_def_var_2d, &
                            nc_def_var_3d, rdb_def_var_1d, nc_enddef, rdb_put_var_1d, &
                            nc_put_att
   use netcdf, only: nf90_put_var
   implicit none

   integer, parameter :: NX_G = 12, NY_G = 12, NGHOST = 3
   real(wp), parameter :: DLON_MODEL = 3.0_wp, DLAT_MODEL = 3.0_wp
   real(wp), parameter :: LON_WEST = 100.0_wp, LAT_SOUTH = -15.0_wp
   real(wp), parameter :: RAD_EARTH = 6.378e6_wp
   integer, parameter :: NLON_SRC = 8, NLAT_SRC = 6
   real(wp), parameter :: SST_DEGC = 15.0_wp, H0 = 10.0_wp
   character(len=*), parameter :: ATM_FILE = "atm_state_mpi_test.nc"
   character(len=*), parameter :: FILELIST = "atm_state_mpi_test_filelist.txt"
   character(len=*), parameter :: VARNAMES(9) = &
                                  [character(len=4) :: "tas", "huss", "uas", "vas", "psl", &
                                                        "rlds", "rsds", "prra", "prsn"]
   ! f(lon,lat,t) = A + B*lon + C*lat + D*t, one row per field (K, kg/kg,
   ! m/s, m/s, Pa, W/m^2, W/m^2, kg/m^2/s, kg/m^2/s) -- real spatial
   ! variation across the 12x12 sector AND across the 2x2 rank seam.
   real(wp), parameter :: COEFS(4, 9) = reshape([ &
                                                290.0_wp, 0.05_wp, 0.08_wp, 0.02_wp, &
                                                8.0e-3_wp, 1.0e-5_wp, 2.0e-5_wp, 1.0e-6_wp, &
                                                5.0_wp, 0.02_wp, 0.03_wp, 0.01_wp, &
                                                3.0_wp, 0.01_wp, 0.02_wp, 0.01_wp, &
                                                1.01e5_wp, 2.0_wp, 3.0_wp, 0.5_wp, &
                                                300.0_wp, 0.3_wp, 0.2_wp, 0.05_wp, &
                                                200.0_wp, 0.2_wp, 0.3_wp, 0.05_wp, &
                                                2.0e-5_wp, 2.0e-8_wp, 3.0e-8_wp, 1.0e-8_wp, &
                                                1.0e-6_wp, 1.0e-9_wp, 2.0e-9_wp, 1.0e-10_wp], &
                                                [4, 9])
   real(wp), parameter :: TQ_DAYS = 5.0_wp

   type :: snap_t
      real(wp), allocatable :: tas(:, :), tau_x(:, :), tau_y(:, :), q_sens(:, :)
      integer :: io = 0, jo = 0, nxl = 0, nyl = 0, ng = 0
   end type snap_t

   integer :: rank, nprocs, n_fail, total_fail, nbad
   type(comm_t) :: comm
   type(snap_t) :: ref, dec

   call comm_env_init()
   call comm_env_setup_roles(.false.)
   rank = comm_env_rank()
   nprocs = comm_env_size()
   comm = comm_env_compute_comm()
   n_fail = 0

   if (nprocs /= 1 .and. nprocs /= 4) then
      error stop "test_ocean_atm_state_mpi: must run at 1 or 4 ranks"
   end if

   if (rank == 0) call write_atm_file()
   call comm%barrier()

   ! Serial reference: every physical rank independently builds the WHOLE
   ! 12x12 sector on one logical tile (px=py=1, crank=0) -- no cross-rank
   ! communication, same trick `test_ocean_decomp_bitid_mpi::run_one` uses
   ! for its csize=1 call.
   call run_one(1, 1, 0, ref)

   ! Teeth: a reference that came out all-zero (e.g. kernels run on unmapped
   ! state under mem:separate) would compare equal to an equally broken
   ! decomposed run.  Require a live, finite stress and sensible-heat field.
   if (.not. all(ieee_is_finite(ref%tau_x)) .or. .not. all(ieee_is_finite(ref%q_sens)) .or. &
       maxval(abs(ref%tau_x)) <= 0.0_wp .or. maxval(abs(ref%tau_y)) <= 0.0_wp .or. &
       maxval(abs(ref%q_sens)) <= 0.0_wp) then
      if (rank == 0) write (*, '(a)') "test_ocean_atm_state_mpi: serial reference "// &
         "tau_x/tau_y/q_sens is zero or non-finite -- comparison would be vacuous"
      n_fail = n_fail + 1
   end if

   if (nprocs == 4) then
      call run_one(2, 2, rank, dec)
      call compare2d(dec%tas, dec%io, dec%jo, dec%nxl, dec%nyl, dec%ng, ref%tas, &
                     "atm%tas", nbad)
      n_fail = n_fail + nbad
      call compare2d(dec%tau_x, dec%io, dec%jo, dec%nxl, dec%nyl, dec%ng, ref%tau_x, &
                     "ss%tau_x", nbad)
      n_fail = n_fail + nbad
      call compare2d(dec%tau_y, dec%io, dec%jo, dec%nxl, dec%nyl, dec%ng, ref%tau_y, &
                     "ss%tau_y", nbad)
      n_fail = n_fail + nbad
      call compare2d(dec%q_sens, dec%io, dec%jo, dec%nxl, dec%nyl, dec%ng, ref%q_sens, &
                     "sf%q_sens", nbad)
      n_fail = n_fail + nbad
   end if

   call comm%barrier()
   total_fail = n_fail
   call allreduce(comm, total_fail, MPI_SUM)
   if (rank == 0) then
      if (total_fail > 0) then
         write (*, '(a,i0,a,i0,a)') "test_ocean_atm_state_mpi: ", total_fail, &
            " check(s) FAILED (nprocs=", nprocs, ")"
      else
         write (*, '(a,i0,a)') "test_ocean_atm_state_mpi: all checks PASSED (nprocs=", &
            nprocs, ")"
      end if
   end if
   call comm_env_finalize()
   if (total_fail > 0) error stop 1

contains

   subroutine write_atm_file()
      !! One NetCDF file, all nine CMOR variable names, two time records
      !! (`t_axis` days since the epoch below) -- same affine-field
      !! convention as `tests/test_ocean_atm_state.F90::write_atm_grid_file`.
      real(wp) :: lon_c(NLON_SRC), lat_c(NLAT_SRC)
      real(wp) :: lon_bnds(2, NLON_SRC), lat_bnds(2, NLAT_SRC)
      real(wp) :: t_axis(2), f(NLON_SRC, NLAT_SRC, 2)
      integer :: ncid, dim_lon, dim_lat, dim_time, dim_bnds
      integer :: vid_lon, vid_lat, vid_lonb, vid_latb, vid_time, vid_f(9)
      integer :: v, i, j, k, ierr

      call uniform_axis(NLON_SRC, 0.0_wp, 360.0_wp, lon_c, lon_bnds)
      call uniform_axis(NLAT_SRC, -90.0_wp, 90.0_wp, lat_c, lat_bnds)
      t_axis = [0.0_wp, 10.0_wp]

      call nc_create_file(ATM_FILE, ncid)
      call nc_def_dim(ncid, "lon", NLON_SRC, dim_lon)
      call nc_def_dim(ncid, "lat", NLAT_SRC, dim_lat)
      call nc_def_dim(ncid, "bnds", 2, dim_bnds)
      call nc_def_dim(ncid, "time", 2, dim_time)
      call rdb_def_var_1d(ncid, "lon", dim_lon, vid_lon)
      call rdb_def_var_1d(ncid, "lat", dim_lat, vid_lat)
      call nc_def_var_2d(ncid, "lon_bnds", [dim_bnds, dim_lon], vid_lonb)
      call nc_def_var_2d(ncid, "lat_bnds", [dim_bnds, dim_lat], vid_latb)
      call rdb_def_var_1d(ncid, "time", dim_time, vid_time)
      call nc_put_att(ncid, vid_time, "units", "days since 2000-01-01 00:00:00")
      do v = 1, 9
         call nc_def_var_3d(ncid, trim(VARNAMES(v)), [dim_lon, dim_lat, dim_time], vid_f(v))
      end do
      call nc_enddef(ncid)

      call rdb_put_var_1d(ncid, vid_lon, lon_c)
      call rdb_put_var_1d(ncid, vid_lat, lat_c)
      ierr = nf90_put_var(ncid, vid_lonb, lon_bnds)
      ierr = nf90_put_var(ncid, vid_latb, lat_bnds)
      call rdb_put_var_1d(ncid, vid_time, t_axis)

      do v = 1, 9
         do k = 1, 2
            do j = 1, NLAT_SRC
               do i = 1, NLON_SRC
                  f(i, j, k) = COEFS(1, v) + COEFS(2, v)*lon_c(i) + COEFS(3, v)*lat_c(j) + &
                               COEFS(4, v)*t_axis(k)
               end do
            end do
         end do
         ierr = nf90_put_var(ncid, vid_f(v), f)
      end do
      call nc_close(ncid)

      call write_filelist_txt(FILELIST, [ATM_FILE])
   end subroutine write_atm_file

   pure subroutine uniform_axis(n, lo, hi, c, bnds)
      integer, intent(in) :: n
      real(wp), intent(in) :: lo, hi
      real(wp), intent(out) :: c(n), bnds(2, n)
      real(wp) :: dxv
      integer :: k
      dxv = (hi - lo)/real(n, wp)
      do k = 1, n
         bnds(1, k) = lo + real(k - 1, wp)*dxv
         bnds(2, k) = lo + real(k, wp)*dxv
         c(k) = 0.5_wp*(bnds(1, k) + bnds(2, k))
      end do
   end subroutine uniform_axis

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

   subroutine run_one(px, py, crank, snap)
      !! Build this (px,py,crank) tile's grid/metrics/halo, configure +
      !! update the atm provider, hand-build the bulk-flux consumer's
      !! inputs (uniform SST/current -- not the point of this test) and
      !! call `bulk_forcing_apply`, then snapshot the four fields' FULL
      !! local arrays (ghosts included) to host.
      integer, intent(in) :: px, py, crank
      type(snap_t), intent(out) :: snap

      type(decomp_t) :: decomp
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(config_t) :: cfg
      type(ocean_data_input_t) :: reader
      type(atm_state_t) :: atm
      type(ocean_bc_state_t) :: bc
      type(multilayer_state_t) :: ms
      type(ocean_sea_ice_t) :: ice
      type(bulk_forcing_t) :: bf
      type(ocean_surface_stress_t) :: ss
      type(ocean_surface_flux_t) :: sf
      real(wp) :: tq
      integer :: ierr

      call decomp_init(decomp, NX_G, NY_G, px, py, crank)
      call grid%init(decomp%nx_local, decomp%ny_local, NGHOST, DLON_MODEL, DLAT_MODEL)
      grid%i_offset_global = decomp%i_start - 1
      grid%j_offset_global = decomp%j_start - 1
      grid%nx_global = decomp%nx_global
      grid%ny_global = decomp%ny_global

      call ocean_halo_destroy()
      call ocean_halo_init(decomp, NGHOST, .false., .false.)

      call make_spherical_metrics(metrics, grid, LON_WEST, LAT_SOUTH, &
                                  DLON_MODEL, DLAT_MODEL, RAD_EARTH)

      cfg%start_date = "2000-01-01"
      cfg%ocean%atm_forcing%enable = .true.
      cfg%ocean%atm_forcing%sw_interp = "linear"
      cfg%ocean%atm_forcing%fluxes_regrid = "bilinear"
      cfg%ocean%atm_forcing%time_mode = "interannual"
      cfg%ocean%atm_forcing%filelist_tas = FILELIST
      cfg%ocean%atm_forcing%filelist_huss = FILELIST
      cfg%ocean%atm_forcing%filelist_uas = FILELIST
      cfg%ocean%atm_forcing%filelist_vas = FILELIST
      cfg%ocean%atm_forcing%filelist_psl = FILELIST
      cfg%ocean%atm_forcing%filelist_rlds = FILELIST
      cfg%ocean%atm_forcing%filelist_rsds = FILELIST
      cfg%ocean%atm_forcing%filelist_prra = FILELIST
      cfg%ocean%atm_forcing%filelist_prsn = FILELIST

      call reader%init()
      call atm_state_configure(cfg, reader, grid, metrics, atm, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) then
         error stop "test_ocean_atm_state_mpi: atm_state_configure failed"
      end if

      ms%nz_ml = 1
      call ms%init(grid)
      ms%h_layer = H0
      ms%tracers(ms%idx_temperature)%hTr = SST_DEGC*H0
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp

      bf%enable = .true.
      bf%nx_total = grid%nx_total
      bf%ny_total = grid%ny_total
      bf%cfg_atm%wind_relative = .true.
      bf%cfg_atm%ocean_albedo = 0.066_wp
      bf%cfg_bulk%n_iter = 5
      bf%cfg_bulk%z_wind = 10.0_wp
      bf%cfg_bulk%z_ta = 2.0_wp
      bf%cfg_bulk%emissivity = 0.97_wp
      allocate (bf%taux_c(grid%nx_total, grid%ny_total), source=0.0_wp)
      allocate (bf%tauy_c(grid%nx_total, grid%ny_total), source=0.0_wp)
      allocate (bf%open_water(grid%nx_total, grid%ny_total), source=1.0_wp)

      call ss%init(grid)
      call sf%init(grid)
      call sf%set_components(grid, .true.)

      !$acc enter data copyin(reader)
      call reader%enter_data()
      !$acc enter data copyin(atm)
      call atm%enter_data()

      tq = TQ_DAYS*86400.0_wp
      call ocean_data_input_update_all(reader, tq)
      call atm_state_update(atm, reader, metrics, tq)
      !$acc update self(atm%tas)

      ! Production device residency (as `test_ocean_bulk_forcing::run_apply`):
      ! the GPU build is mem:separate, so every object `bulk_forcing_apply`
      ! touches must be mapped -- the stress seam refresh in particular.
      ! `ice` stays disabled here, so it is not mapped.
      !$acc enter data copyin(ms, bf, ss, sf)
      call ms%enter_data()
      call bf%enter_data()
      call ss%enter_data()
      call sf%enter_data()

      call bulk_forcing_apply(grid, metrics, bc, atm, ms, ice, bf, ss, sf)
      !$acc update self(ss%tau_x, ss%tau_y)
      !$acc update self(sf%q_sens)

      call sf%exit_data()
      call ss%exit_data()
      call bf%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, bf, ss, sf)

      allocate (snap%tas(size(atm%tas, 1), size(atm%tas, 2)))
      snap%tas = atm%tas
      allocate (snap%tau_x(size(ss%tau_x, 1), size(ss%tau_x, 2)))
      snap%tau_x = ss%tau_x
      allocate (snap%tau_y(size(ss%tau_y, 1), size(ss%tau_y, 2)))
      snap%tau_y = ss%tau_y
      allocate (snap%q_sens(size(sf%q_sens, 1), size(sf%q_sens, 2)))
      snap%q_sens = sf%q_sens
      snap%io = grid%i_offset_global
      snap%jo = grid%j_offset_global
      snap%nxl = grid%nx_phys
      snap%nyl = grid%ny_phys
      snap%ng = grid%nghost

      call atm%exit_data()
      !$acc exit data delete(atm)
      call reader%exit_data()
      !$acc exit data delete(reader)
      call atm%destroy()
      call reader%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine run_one

   subroutine compare2d(dec_a, io, jo, nxl, nyl, ng, ref_a, tag, nbad_out)
      !! Bitwise comparison of `dec_a`'s owned window against the matching
      !! slice of the (global-sized) reference array `ref_a`.  `dec_a` may
      !! be a cell array (same extent as `(nxl+2*ng, nyl+2*ng)`) or a face
      !! array one wider along its stagger axis -- `ex`/`ey` pick that up
      !! automatically, same convention as
      !! `test_ocean_decomp_bitid_mpi::compare`.
      real(wp), intent(in) :: dec_a(:, :), ref_a(:, :)
      integer, intent(in) :: io, jo, nxl, nyl, ng
      character(len=*), intent(in) :: tag
      integer, intent(out) :: nbad_out
      integer :: i, j, ex, ey, nb
      character(len=160) :: first

      ex = size(dec_a, 1) - (nxl + 2*ng)
      ey = size(dec_a, 2) - (nyl + 2*ng)
      nb = 0
      first = ""
      do j = ng + 1, ng + nyl + ey
         do i = ng + 1, ng + nxl + ex
            if (transfer(dec_a(i, j), 0_int64) /= transfer(ref_a(i + io, j + jo), 0_int64)) then
               if (nb == 0) write (first, '(a,2(i0,1x),a,es24.16,a,es24.16)') &
                  "first at local (i,j)=", i, j, ": ", dec_a(i, j), " vs ", ref_a(i + io, j + jo)
               nb = nb + 1
            end if
         end do
      end do
      if (nb > 0) then
         write (*, '(a,i0,3a,i0,2a)') "  rank ", rank, " ", tag, ": mismatches=", nb, "  ", &
            trim(first)
      end if
      nbad_out = nb
   end subroutine compare2d

end program test_ocean_atm_state_mpi
#endif
