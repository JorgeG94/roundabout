!! The CICE-port incremental-remapping transport on several ranks IS the serial run, bitwise
!! (D10.6; I3 part (c): a blob through the tripolar fold seam, decomposed == serial).
!!
!! On EACH launch (ctest runs np = 1, 2, 4), for every case below the program
!!   1. builds the geometry through the production engine (`engine_setup`: grid, metrics, land
!!      mask, the C-grid ice metrics, the halo, the fold exchange), once as the SERIAL reference
!!      (`compute_size = 1` on every rank, a private self-communicator for the collectives of
!!      the test itself) and once per px x py factorisation of the launched rank count
!!      (np = 2: 2x1, 1x2; np = 4: 4x1, 1x4, 2x2);
!!   2. seeds a category state and a face velocity from FUNCTIONS OF THE GLOBAL INDEX (so every
!!      layout holds the same physical state), then runs `NSTEP` transport steps through
!!      `ice_step_transport` exactly as the engine's hooks sequence them: the face velocity pair
!!      through `ice_halo_faces`, the category state through `ice_halo_category_state` before every
!!      step, the transport seam bound by `ice_transport_seam_bind`, `transport_check` on (EFP
!!      conservation sums, monotonicity bounds);
!!   3. compares, bitwise, every owned cell of the final category state, the step result (`ierr`,
!!      `worst_rel`: the conservation EFP sums reduced over the ranks), and, on the fold cases, that
!!      the north fluxes of the fold row are exactly antisymmetric under the fold mirror (the paired
!!      face: what leaves one cell enters the other in the same bits).
!!
!! Cases: `channel` (periodic x, walls south and north, a seamount: the blob crosses the wrap seam
!! and every rank seam); `island` (a closed basin with an interior land block: land across seams);
!! `fold` (the tripolar cap, a blob in the east half driven north through the fold row into the west
!! half; `px = 1` runs the local pairing, `px > 1` the distributed fold exchange); `fold_upwind`
!! (the same with `advection = "upwind"`: its fold-row flux is not paired, so the transport REFUSES the step with
!! `TRANSPORT_CFG_FOLD` on the serial run and on every layout, with the same word).
!!
!! Negative controls (each must FAIL its comparison, so the comparison can fail): the category state
!! fold skipped, the face velocity fold skipped, the seam's reductions off (the EFP conservation
!! numbers of a tile are then not the global ones), and, on the fold cases, the flux pairing off (the
!! fold row is then NOT exactly antisymmetric).  The `fold` case is also checked non-vacuous: the
!! reference has ice in the west half of the top rows that started empty.
#ifdef RDB_ENABLE_MPI
program test_cice_remap_decomp_mpi
   use, intrinsic :: iso_fortran_env, only: int32, int64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
   use rdb_constants, only: wp
   use rdb_config, only: config_t, read_config_from_toml_string, validate_config
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   use rdb_ocean_engine, only: ocean_engine_t, engine_setup, engine_enter_data, engine_exit_data, engine_teardown
   use rdb_comm_env, only: comm_env_init, comm_env_setup_roles, comm_env_finalize, comm_env_rank, comm_env_size, &
                           comm_env_compute_comm, comm_env_push_compute_comm, comm_env_pop_compute_comm
   use pic_mpi_lib, only: comm_t, allreduce, MPI_SUM
   use rdb_ice_params, only: ice_params_t, ice_params_default, ICE_CONSTANTS_ICEPACK, ICE_ADVECT_REMAP, ICE_ADVECT_UPWIND
   use rdb_ice_transport_cice, only: transport_cice_ws_t, transport_result_t, ice_step_transport, ICE_TRANSPORT_POST_NONE, &
                                     TRANSPORT_CFG_OK, TRANSPORT_CFG_FOLD
   use rdb_ice_cgrid_metrics, only: ice_halo_faces
   use rdb_ice_transport_seam_comm, only: ice_transport_seam_bind, ice_halo_category_state
   use rdb_ocean_boundary_types, only: ocean_bc_state_t
   implicit none

   integer, parameter :: NCAT = 5, NIL = 3, NSTEP = 8
   real(wp), parameter :: DT = 1000.0_wp
   character(len=*), parameter :: NL = new_line("a")

   type :: snap_t
      !! Owned-cell window of one run (host copies), its offsets and the step results.
      integer :: nxl = 0, nyl = 0, ng = 0, io = 0, jo = 0, nxg = 0, nyg = 0
      real(wp), allocatable :: aicen(:, :, :), vicen(:, :, :), vsnon(:, :, :), Tsfcn(:, :, :), aice0(:, :)
      real(wp), allocatable :: qicen(:, :, :, :), sicen(:, :, :, :), qsnon(:, :, :, :)
      integer(int32) :: ierr = 0
      integer :: cfg = 0, nbad_steps = 0, n_err_cells = 0
      real(wp) :: worst_rel = 0.0_wp
      integer :: worst_q = 0
      real(wp) :: mass0 = 0.0_wp, mass1 = 0.0_wp
         !! Global sum of `aicen*tarea` before and after the run (plain sum, for the round-off check).
      real(wp) :: anti_max = -1.0_wp, flux_max = -1.0_wp, anti_t = -1.0_wp, flux_t = -1.0_wp
         !! Fold row of the last step's north area fluxes: `max |F(i) + F(mirror)|` and `max |F|` (-1: no fold).
      logical :: ok = .false.
   end type snap_t

   integer :: rank, nprocs, n_fail, total_fail
   type(comm_t) :: comm
   character(len=64) :: case_filter = ""
      !! `RDB_REMAP_CASES=<substring>`: run only the cases whose label contains it (development filter).

   call comm_env_init()
   call comm_env_setup_roles(.false.)
   rank = comm_env_rank()
   nprocs = comm_env_size()
   comm = comm_env_compute_comm()
   n_fail = 0
   call get_environment_variable("RDB_REMAP_CASES", case_filter)

   call run_case("channel", .false.)
   call run_case("island", .false.)
   call run_case("fold", .false.)
   call run_case("fold_upwind", .false.)
   call run_case("fold_smallbuf", .false.)
   call run_case("fold_odd", .false.)
   call run_case("fold_odd_smallbuf", .false.)
   call run_case("channel_nan", .false.)
   call run_case("channel_cfl", .false.)
   call run_case("fold_nan", .false.)
   call run_case("fold", .true.)
   call run_case("fold_odd", .true.)

   call comm%barrier()
   total_fail = n_fail
   call allreduce(comm, total_fail, MPI_SUM)
   if (rank == 0) then
      if (total_fail > 0) then
         write (*, '(a,i0,a,i0,a)') "test_cice_remap_decomp_mpi: ", total_fail, " check(s) FAILED (nprocs=", nprocs, ")"
      else
         write (*, '(a,i0,a)') "test_cice_remap_decomp_mpi: all checks PASSED (nprocs=", nprocs, ")"
      end if
   end if
   call comm_env_finalize()
   if (total_fail > 0) error stop 1

contains

   ! ------------------------------------------------------------------
   ! namelists
   ! ------------------------------------------------------------------

   function case_nml(label, px, py) result(nml)
      !! A sea-ice case with the C-grid EVP slot (it owns the C-grid ice metrics) on a px x py process grid.
      character(len=*), intent(in) :: label
      integer, intent(in) :: px, py
      character(len=:), allocatable :: nml
      character(len=16) :: spx, spy
      character(len=:), allocatable :: grid, topo, bc

      write (spx, '(i0)') px
      write (spy, '(i0)') py
      select case (label)
      case ("fold", "fold_upwind", "fold_smallbuf", "fold_odd", "fold_odd_smallbuf", "fold_nan")
         grid = "[grid]"//NL// &
                "nx = "//merge("29", "30", index(label, "odd") > 0)//NL// &
                "ny = 24"//NL// &
                "nghost = 3"//NL// &
                "dx = "//merge("12.413793", "12.0     ", index(label, "odd") > 0)//NL// &
                "dy = 1.0"//NL// &
                "[ocean.grid]"//NL// &
                'grid_config = "tripolar"'//NL// &
                "lon_west = 0.0"//NL// &
                "lat_south = 59.0"//NL// &
                "phi_join = 74.0"//NL// &
                "lon_pole = 0.0"//NL// &
                "rad_earth = 6.378e6"//NL// &
                'coriolis_scheme = "planetary"'//NL
         topo = "[ocean.topo]"//NL// &
                'topo_config = "double_drake"'//NL// &
                "max_depth = 1000.0"//NL// &
                "slope_scale = 0.2"//NL
         bc = "[ocean.bc]"//NL// &
              'west = "periodic"'//NL// &
              'east = "periodic"'//NL// &
              'south = "wall"'//NL// &
              'north = "tripolar_fold"'//NL
      case ("island")
         grid = "[grid]"//NL// &
                "nx = 26"//NL// &
                "ny = 18"//NL// &
                "nghost = 3"//NL// &
                "dx = 20000.0"//NL// &
                "dy = 20000.0"//NL
         topo = "[ocean.topo]"//NL// &
                'topo_config = "island"'//NL// &
                "max_depth = 1000.0"//NL// &
                "slope_scale = 0.15"//NL
         bc = "[ocean.bc]"//NL// &
              'west = "wall"'//NL// &
              'east = "wall"'//NL// &
              'south = "wall"'//NL// &
              'north = "wall"'//NL
      case default   ! channel, channel_nan, channel_cfl
         grid = "[grid]"//NL// &
                "nx = 26"//NL// &
                "ny = 18"//NL// &
                "nghost = 3"//NL// &
                "dx = 20000.0"//NL// &
                "dy = 20000.0"//NL
         topo = "[ocean.topo]"//NL// &
                'topo_config = "seamount"'//NL// &
                "max_depth = 1000.0"//NL// &
                "edge_depth = 800.0"//NL// &
                "slope_scale = 60000.0"//NL
         bc = "[ocean.bc]"//NL// &
              'west = "periodic"'//NL// &
              'east = "periodic"'//NL// &
              'south = "wall"'//NL// &
              'north = "wall"'//NL
      end select
      nml = "[sim]"//NL// &
            'sim_type = "ocean"'//NL// &
            "[mpi]"//NL// &
            "px = "//trim(spx)//NL// &
            "py = "//trim(spy)//NL// &
            "[time]"//NL// &
            "t_end = 86400.0"//NL// &
            "dt_fixed = 900.0"//NL// &
            "[nonhydrostatic]"//NL// &
            "nz_layers = 4"//NL// &
            "[tracer]"//NL// &
            "initial_temperature = -1.83"//NL// &
            "initial_salinity = 34.0"//NL// &
            "T_init_surface = -1.82"//NL// &
            "T_init_bottom = -1.835"//NL// &
            "[ocean.bt]"//NL// &
            "auto_n_inner = true"//NL// &
            'split_scheme = "pred_corr"'//NL// &
            "[ocean.hvisc]"//NL// &
            "nu_h = 200.0"//NL// &
            'lateral_closure = "smagorinsky"'//NL// &
            "smag_ah = true"//NL// &
            "[ocean.diag]"//NL// &
            "enabled = false"//NL// &
            grid// &
            "[physics]"//NL// &
            "coriolis_f = 1.4e-4"//NL// &
            "wind_stress_x = 0.1"//NL// &
            "wind_stress_y = 0.04"//NL// &
            "[vcoord]"//NL// &
            'vcoord_type = "sigma"'//NL// &
            topo// &
            "[ocean.thermo]"//NL// &
            "enable_thermodynamics = true"//NL// &
            "q_heat = -100.0"//NL// &
            "[ocean.vmix]"//NL// &
            "dt_therm_ratio = 2"//NL// &
            "[ocean.ice]"//NL// &
            "enable = true"//NL// &
            "ncat = 5"//NL// &
            "nk_ice = 2"//NL// &
            "air_temp = -20.0"//NL// &
            "restore_lambda = 20.0"//NL// &
            "sw_down = 0.0"//NL// &
            "dynamics = true"//NL// &
            "evp_sub_steps = 30"//NL// &
            "cfl_trunc = 0.01"//NL// &
            "project_ci = true"//NL// &
            "transport = true"//NL// &
            "adv_substeps = 2"//NL// &
            "[ocean.ice_ic]"//NL// &
            'conc_config = "uniform"'//NL// &
            "h_ice = 1.0"//NL// &
            "conc = 0.7"//NL// &
            "[ocean.ice_ridge]"//NL// &
            "kstrength = 0"//NL// &
            "[ocean.ice_dyn]"//NL// &
            'dyn_scheme = "cice_c"'//NL// &
            "ndte = 10"//NL// &
            bc
   end function case_nml

   ! ------------------------------------------------------------------
   ! seeding functions of the GLOBAL index
   ! ------------------------------------------------------------------

   pure function hash01(a, b, c) result(r)
      !! Deterministic pseudo-random number in [0, 1) from three integers.
      integer, intent(in) :: a, b, c
      real(wp) :: r
      integer(int64) :: h

      h = iand(int(a, int64)*73856093_int64 + int(b, int64)*19349663_int64 + int(c, int64)*83492791_int64 + 12345_int64, &
               1099511627775_int64)
      h = ieor(h, ishft(h, -13))
      h = iand(h*1274126177_int64, 1099511627775_int64)
      h = ieor(h, ishft(h, -16))
      r = real(iand(h, 1048575_int64), wp)/1048576.0_wp
   end function hash01

   pure function blob(gi, gj, nxg, nyg, fold_case) result(c)
      !! Ice concentration pattern of the global cell (gi, gj): three bumps: one in the east half near the top
      !! rows (the fold case drives it through the fold), one over the periodic edge at mid height, one on the
      !! middle meridian (it sits on the 2x1 / 2x2 seams).
      integer, intent(in) :: gi, gj, nxg, nyg
      logical, intent(in) :: fold_case
      real(wp) :: c
      real(wp) :: r1, r2, r3, x1, y1, x2, y2, x3, y3

      x1 = 0.78_wp*real(nxg, wp)
      y1 = real(nyg, wp) - merge(2.0_wp, 4.0_wp, fold_case)
      x2 = 0.08_wp*real(nxg, wp)
      y2 = 0.45_wp*real(nyg, wp)
      x3 = 0.5_wp*real(nxg, wp)
      y3 = 0.55_wp*real(nyg, wp)
      r1 = sqrt((real(gi, wp) - x1)**2 + (real(gj, wp) - y1)**2)
      r2 = sqrt((real(gi, wp) - x2)**2 + (real(gj, wp) - y2)**2)
      r3 = sqrt((real(gi, wp) - x3)**2 + (real(gj, wp) - y3)**2)
      c = 0.0_wp
      if (r1 < 6.5_wp) c = max(c, 0.9_wp*sqrt(1.0_wp - r1/6.5_wp))
      if (r2 < 6.5_wp) c = max(c, 0.8_wp*sqrt(1.0_wp - r2/6.5_wp))
      if (r3 < 5.5_wp) c = max(c, 0.7_wp*sqrt(1.0_wp - r3/5.5_wp))
   end function blob

   pure function wave(x) result(w)
      !! A smooth period-1 wave of amplitude about 1 from polynomials only (no libm: a vectorised libm rounds by loop
      !! position, which would make the SEED depend on the tile width).
      real(wp), intent(in) :: x
      real(wp) :: w, t
      t = x - floor(x)
      w = 6.4_wp*t*(1.0_wp - t)*(1.0_wp - 2.0_wp*t)
   end function wave

   pure function ufrac(gif, gj, nxg, nyg) result(f)
      !! Eastward displacement of an E face as a fraction of the cell width per step (global face index `gif`
      !! = the global index of the cell whose WEST face it is).
      integer, intent(in) :: gif, gj, nxg, nyg
      real(wp) :: f
      real(wp) :: x, y
      x = real(gif - 1, wp)/real(nxg, wp)
      y = real(gj, wp)/real(nyg, wp)
      f = 0.07_wp*wave(x + 0.06_wp) + 0.04_wp*wave(y + 0.25_wp)
   end function ufrac

   pure function vfrac(gi, gjf, nxg, nyg) result(f)
      !! Northward displacement of an N face as a fraction of the cell height per step (`gjf` = the global index
      !! of the cell whose SOUTH face it is).  Strongly northward near the top rows (into the fold).
      integer, intent(in) :: gi, gjf, nxg, nyg
      real(wp) :: f
      real(wp) :: x
      x = real(gi, wp)/real(nxg, wp)
      f = 0.03_wp*wave(x + 0.42_wp) &
          + 0.10_wp*max(0.0_wp, real(gjf - (nyg - 6), wp))/6.0_wp + 0.03_wp
   end function vfrac

   ! ------------------------------------------------------------------
   ! one run
   ! ------------------------------------------------------------------

   subroutine run_one(label, px, py, csize, crank, ctrl, snap)
      !! Configure, seed, step `NSTEP` times, snapshot the owned window.  `ctrl`: 0 normal; 1 category state
      !! fold skipped; 2 face velocity fold skipped; 3 the seam's reductions off; 4 the flux pairing off.
      character(len=*), intent(in) :: label
      integer, intent(in) :: px, py, csize, crank, ctrl
      type(snap_t), intent(out) :: snap
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      type(ice_params_t) :: p
      type(transport_cice_ws_t) :: w
      type(transport_result_t) :: res
      type(ocean_bc_state_t) :: bc_state, bc_vel
      type(comm_t) :: real_comm, self_comm
      real(wp), allocatable :: u(:, :), v(:, :), Tf(:, :), hin_max(:), aice0(:, :), aicen(:, :, :), vicen(:, :, :)
      real(wp), allocatable :: vsnon(:, :, :), Tsfcn(:, :, :), qicen(:, :, :, :), sicen(:, :, :, :), qsnon(:, :, :, :)
      real(wp), allocatable :: aice(:, :), vice(:, :), vsno(:, :), Tsfc_c(:, :), qice_c(:, :, :), sice_c(:, :, :)
      real(wp), allocatable :: qsno_c(:, :, :), fresh(:, :), fsalt(:, :), fhocn(:, :)
      integer(int32), allocatable :: ice_err(:, :)
      integer :: ierr, nx, ny, ng, is, ie, js, je, i, j, n, k, step, gi, gj, nxg, nyg
      real(wp) :: s, tot
      logical :: fold_case, upwind
      integer :: inj, nt_probe
      real(wp) :: vscale
      integer(int64) :: bmax
      type(transport_cice_ws_t) :: probe

      snap%ok = .false.
      fold_case = index(label, "fold") == 1
      upwind = label == "fold_upwind"
      inj = 0
      vscale = 1.0_wp
      if (index(label, "_nan") > 0) inj = 1
      if (index(label, "_cfl") > 0) inj = 2
      call read_config_from_toml_string(case_nml(label, px, py), cfg, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call validate_config(cfg, ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call engine_setup(engine, cfg, ierr, compute_rank=crank, compute_size=csize)
      if (ierr /= OCEAN_STATUS_OK) return
      call engine_enter_data(engine, cfg)

      associate (g => engine%grid, cm => engine%state%ice_cgrid)
         nx = g%nx_total
         ny = g%ny_total
         ng = g%nghost
         is = ng + 1
         ie = ng + g%nx_phys
         js = ng + 1
         je = ng + g%ny_phys
         nxg = g%nx_global
         nyg = g%ny_global
         snap%nxl = g%nx_phys
         snap%nyl = g%ny_phys
         snap%ng = ng
         snap%io = g%i_offset_global
         snap%jo = g%j_offset_global
         snap%nxg = nxg
         snap%nyg = nyg

         p = ice_params_default(ICE_CONSTANTS_ICEPACK, NCAT, NIL, 1)
         p%transport%advection = merge(ICE_ADVECT_UPWIND, ICE_ADVECT_REMAP, upwind)
         p%transport%check = .true.
         if (index(label, "smallbuf") > 0) then
            ! a flux-buffer budget that gives tiles of different width DIFFERENT local batch sizes (x tile widths of
            ! 40 ghosted columns: 13-wide tiles take 3 category indices, 14-wide ones 2): the batch must still be uniform
            call probe%init(nx, ny, NCAT, NIL, 1, .true., ierr, remap=.true.)
            nt_probe = probe%ntrace
            call probe%destroy()
            bmax = int(4 + 5*nt_probe, int64)*40_int64*int(ny, int64)*8_int64
            call w%init(nx, ny, NCAT, NIL, 1, .true., ierr, remap=.true., remap_flux_bytes_max=bmax)
         else
            call w%init(nx, ny, NCAT, NIL, 1, .true., ierr, remap=.not. upwind)
         end if
         if (ierr /= 0) then
            write (*, '(a,i0)') "FAIL: transport workspace init ierr ", ierr
            call engine_exit_data(engine)
            call engine_teardown(engine)
            return
         end if
         bc_state = engine%state%bc
         bc_vel = engine%state%bc
         if (ctrl == 1) bc_state%north_fold = .false.
         if (ctrl == 2) bc_vel%north_fold = .false.
         call ice_transport_seam_bind(w%seam, g, engine%state%bc, w%fw%nbatch, w%ntrace)
         if (ctrl == 3) w%seam%multirank = .false.
         if (ctrl == 4) w%seam%pair_fold = .false.
         if (ctrl == 5) w%seam%pair_tracers = .false.

         allocate (u(nx + 1, ny), v(nx, ny + 1), Tf(nx, ny), hin_max(0:NCAT))
         allocate (aice0(nx, ny), aicen(nx, ny, NCAT), vicen(nx, ny, NCAT), vsnon(nx, ny, NCAT), Tsfcn(nx, ny, NCAT))
         allocate (qicen(nx, ny, NCAT, NIL), sicen(nx, ny, NCAT, NIL), qsnon(nx, ny, NCAT, 1))
         allocate (aice(nx, ny), vice(nx, ny), vsno(nx, ny), Tsfc_c(nx, ny), qice_c(nx, ny, NIL), sice_c(nx, ny, NIL))
         allocate (qsno_c(nx, ny, 1), fresh(nx, ny), fsalt(nx, ny), fhocn(nx, ny), ice_err(nx, ny))
         Tf = -1.8_wp
         hin_max = 0.0_wp
         aice = 0.0_wp; vice = 0.0_wp; vsno = 0.0_wp; Tsfc_c = 0.0_wp; qice_c = 0.0_wp; sice_c = 0.0_wp; qsno_c = 0.0_wp
         fresh = 0.0_wp; fsalt = 0.0_wp; fhocn = 0.0_wp; ice_err = 0_int32
         aicen = 0.0_wp; vicen = 0.0_wp; vsnon = 0.0_wp; Tsfcn = -1.8_wp; qicen = 0.0_wp; sicen = 0.0_wp; qsnon = 0.0_wp

         ! the state, a function of the global index (owned cells of wet T points)
         do j = js, je
            do i = is, ie
               gi = g%i_offset_global + (i - ng)
               gj = g%j_offset_global + (j - ng)
               if (cm%hm(i, j) < 0.5_wp) cycle
               s = blob(gi, gj, nxg, nyg, fold_case)
               if (s <= 0.0_wp) cycle
               tot = 0.0_wp
               do n = 1, NCAT
                  aicen(i, j, n) = 0.2_wp + hash01(gi, gj, n)
                  tot = tot + aicen(i, j, n)
               end do
               do n = 1, NCAT
                  aicen(i, j, n) = s*aicen(i, j, n)/tot
                  vicen(i, j, n) = aicen(i, j, n)*(0.3_wp + 0.35_wp*real(n - 1, wp) + 0.1_wp*hash01(gi, gj, 10 + n))
                  vsnon(i, j, n) = aicen(i, j, n)*(0.04_wp + 0.06_wp*hash01(gi, gj, 20 + n))
                  Tsfcn(i, j, n) = -3.0_wp + 2.0_wp*hash01(gi, gj, 30 + n)
                  do k = 1, NIL
                     qicen(i, j, n, k) = -(2.0e8_wp + 6.0e7_wp*hash01(gi, gj, 40 + 7*n + k))
                     sicen(i, j, n, k) = 3.0_wp + 4.0_wp*hash01(gi, gj, 60 + 7*n + k)
                  end do
                  qsnon(i, j, n, 1) = -(1.0e8_wp + 2.0e7_wp*hash01(gi, gj, 80 + n))
               end do
            end do
         end do
         do j = 1, ny
            do i = 1, nx
               aice0(i, j) = merge(max(1.0_wp - sum(aicen(i, j, :)), 0.0_wp), 0.0_wp, cm%hm(i, j) > 0.5_wp)
            end do
         end do

         ! the face velocities: the owned faces from the global face index, scaled by the narrower of the two cells
         ! the face joins, so the displacement is a fixed fraction of it (the CFL check's own measure); closed faces zero; ghosts by the production halo
         u = 0.0_wp
         v = 0.0_wp
         do j = js, je
            do i = is, ie + 1
               gi = g%i_offset_global + (i - ng)
               gj = g%j_offset_global + (j - ng)
               u(i, j) = vscale*ufrac(gi, gj, nxg, nyg)*min(cm%dxT(i - 1, j), cm%dxT(i, j))/DT*cm%epm(i - 1, j)
            end do
         end do
         do j = js, je + 1
            do i = is, ie
               gi = g%i_offset_global + (i - ng)
               gj = g%j_offset_global + (j - ng)
               v(i, j) = vscale*vfrac(gi, gj, nxg, nyg)*min(cm%dyT(i, j - 1), cm%dyT(i, j))/DT*cm%npm(i, j - 1)
            end do
         end do
         if (inj > 0) then
            ! one rank only: a face near the east end of the middle rows
            do j = js, je
               do i = is, ie + 1
                  if (g%i_offset_global + (i - ng) == nxg - 2 .and. g%j_offset_global + (j - ng) == 4) then
                     if (inj == 1) u(i, j) = ieee_value(1.0_wp, ieee_quiet_nan)
                     if (inj == 2) u(i, j) = 3.0_wp*max(cm%dxT(i - 1, j), cm%dxT(i, j))/DT
                  end if
               end do
            end do
         end if
         ! mem:separate: every array the transport kernels touch must be device-present (the engine maps its own
         ! state, not this test's local fields or workspace).  Host edits are all above, so the copyin carries them.
         !$acc enter data copyin(u, v, Tf, hin_max, aice0, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon)
         !$acc enter data copyin(aice, vice, vsno, Tsfc_c, qice_c, sice_c, qsno_c, fresh, fsalt, fhocn, ice_err)
         call w%enter_data()
         call ice_halo_faces(g, bc_vel, u, v, .true.)

         tot = 0.0_wp
         do j = js, je
            do i = is, ie
               tot = tot + sum(aicen(i, j, :))*cm%tarea(i, j)
            end do
         end do
         snap%mass0 = tot

         if (csize == 1) then
            real_comm = comm_env_compute_comm()
            self_comm = real_comm%split_by(real_comm%rank())
            call comm_env_push_compute_comm(self_comm)
         end if

         snap%nbad_steps = 0
         do step = 1, NSTEP
            call ice_halo_category_state(g, bc_state, nx, ny, NCAT, NIL, 1, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, &
                                         aice0, .true.)
            call ice_step_transport(w, p, DT, ICE_TRANSPORT_POST_NONE, is, ie, js, je, u, v, cm%HTE, cm%HTN, cm%tarea, &
                                    cm%dxT, cm%dyT, cm%epm, cm%npm, cm%hm, Tf, hin_max, aice0, aicen, vicen, vsnon, Tsfcn, &
                                    qicen, sicen, qsnon, aice, vice, vsno, Tsfc_c, qice_c, sice_c, qsno_c, fresh, fsalt, &
                                    fhocn, ice_err, res, dxU=cm%dxU, dyU=cm%dyU, earea=cm%earea, narea=cm%narea, uvm=cm%uvm)
            snap%ierr = ior(snap%ierr, res%ierr)
            snap%cfg = max(snap%cfg, res%cfg)
            snap%worst_rel = res%worst_rel
            snap%worst_q = res%worst_q
            if (.not. res%stepped) then
               snap%nbad_steps = snap%nbad_steps + 1
               !$acc update self(ice_err)
               write (*, '(a,i0,a,i0,a,z0,5(a,i0))') "  [transport refused] step ", step, " cfg ", res%cfg, " ierr 0x", res%ierr, &
                  ": n_cfl ", res%n_cfl, " n_geom ", res%n_geom, " n_nonfinite ", res%n_nonfinite, " n_neg ", res%n_neg, &
                  " n_mono ", res%n_mono
               do j = 1, ny
                  do i = 1, nx
                   if (ice_err(i, j) /= 0_int32) write (*, '(a,4(i0,1x),4es11.3)') "    cell (i j gi gj) dxT dyT dxE dyN: ", i, j, &
                      g%i_offset_global + i - ng, g%j_offset_global + j - ng, cm%dxT(i, j), cm%dyT(i, j), cm%dxE(i, j), cm%dyN(i, j)
                  end do
               end do
            end if
            if (.not. res%stepped) exit
         end do
         !$acc update self(aice0, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon)
         !$acc update self(aice, vice, vsno, Tsfc_c, qice_c, sice_c, qsno_c, fresh, fsalt, fhocn, ice_err)
         snap%n_err_cells = count(ice_err(is:ie, js:je) /= 0_int32)
         call ice_halo_category_state(g, bc_state, nx, ny, NCAT, NIL, 1, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon, &
                                      aice0, .true.)

         tot = 0.0_wp
         do j = js, je
            do i = is, ie
               tot = tot + sum(aicen(i, j, :))*cm%tarea(i, j)
            end do
         end do
         snap%mass1 = tot
         ! global sums for the round-off conservation check (collective: every rank, any layout)
         call global_sum(snap%mass0)
         call global_sum(snap%mass1)

         if (fold_case .and. .not. upwind) then
            !$acc update self(w%fw%mflx_n, w%fw%mtflx_n)
         end if
         if (fold_case) call fold_row_stats(w, engine%state%bc%north_fold, g, ctrl, snap%anti_max, snap%flux_max, snap%anti_t, snap%flux_t, upwind)
         if (csize == 1) then
            call comm_env_pop_compute_comm()
            call self_comm%finalize()
         end if

         call w%exit_data()
         !$acc exit data delete(aice, vice, vsno, Tsfc_c, qice_c, sice_c, qsno_c, fresh, fsalt, fhocn, ice_err)
         !$acc exit data delete(u, v, Tf, hin_max, aice0, aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon)
         allocate (snap%aicen(snap%nxl, snap%nyl, NCAT), snap%vicen(snap%nxl, snap%nyl, NCAT))
         allocate (snap%vsnon(snap%nxl, snap%nyl, NCAT), snap%Tsfcn(snap%nxl, snap%nyl, NCAT), snap%aice0(snap%nxl, snap%nyl))
         allocate (snap%qicen(snap%nxl, snap%nyl, NCAT, NIL), snap%sicen(snap%nxl, snap%nyl, NCAT, NIL))
         allocate (snap%qsnon(snap%nxl, snap%nyl, NCAT, 1))
         snap%aicen = aicen(is:ie, js:je, :)
         snap%vicen = vicen(is:ie, js:je, :)
         snap%vsnon = vsnon(is:ie, js:je, :)
         snap%Tsfcn = Tsfcn(is:ie, js:je, :)
         snap%aice0 = aice0(is:ie, js:je)
         snap%qicen = qicen(is:ie, js:je, :, :)
         snap%sicen = sicen(is:ie, js:je, :, :)
         snap%qsnon = qsnon(is:ie, js:je, :, :)
         snap%ok = .true.
      end associate
      call w%destroy()
      call engine_exit_data(engine)
      call engine_teardown(engine)
   end subroutine run_one

   subroutine global_sum(x)
      !! In-place sum over the compute communicator (the pushed self-communicator in the serial reference).
      real(wp), intent(inout) :: x
      type(comm_t) :: c

      c = comm_env_compute_comm()
      if (c%size() == 1) return
      call allreduce(c, x, MPI_SUM)
   end subroutine global_sum

   subroutine fold_row_stats(w, has_fold, g, ctrl, anti_max, flux_max, anti_t, flux_t, upwind)
      !! Gather the fold row of the last step's north AREA fluxes (all category indices of the one batch)
      !! over the ranks and measure `max |F(i) + F(ni+1-i)|` (zero when the face is paired) and `max |F|`.
      use rdb_grid, only: hgrid_t
      type(transport_cice_ws_t), intent(in) :: w
      logical, intent(in) :: has_fold, upwind
      type(hgrid_t), intent(in) :: g
      integer, intent(in) :: ctrl
      real(wp), intent(out) :: anti_max, flux_max, anti_t, flux_t
      real(wp), allocatable :: loc(:), glob(:), loct(:), globt(:)
      integer :: t, nt
      type(comm_t) :: c
      integer :: ni, nb, b, i, im, je, is

      anti_max = -1.0_wp
      flux_max = -1.0_wp
      anti_t = -1.0_wp
      flux_t = -1.0_wp
      if (upwind) return
      ni = g%nx_global
      nb = NCAT + 1
      if (w%fw%nbatch < nb) return
      je = g%nghost + g%ny_phys
      is = g%nghost + 1
      nt = w%ntrace
      allocate (loc(ni*nb), glob(ni*nb), loct(ni*nb*nt), globt(ni*nb*nt))
      loc = 0.0_wp
      loct = 0.0_wp
      if (has_fold) then
         do b = 1, nb
            do i = 1, g%nx_phys
               loc((b - 1)*ni + g%i_offset_global + i) = w%fw%mflx_n(is - 1 + i, je, b)
            end do
         end do
         do b = 1, nb
            do t = 1, nt
               do i = 1, g%nx_phys
                  loct(((b - 1)*nt + t - 1)*ni + g%i_offset_global + i) = w%fw%mtflx_n(is - 1 + i, je, t, b)
               end do
            end do
         end do
      end if
      c = comm_env_compute_comm()
      glob = loc
      globt = loct
      if (c%size() > 1) then
         call allreduce(c, loc, glob, op=MPI_SUM)
         call allreduce(c, loct, globt, op=MPI_SUM)
      end if
      anti_max = 0.0_wp
      flux_max = 0.0_wp
      do b = 1, nb
         do i = 1, ni
            im = ni + 1 - i
            anti_max = max(anti_max, abs(glob((b - 1)*ni + i) + glob((b - 1)*ni + im)))
            flux_max = max(flux_max, abs(glob((b - 1)*ni + i)))
         end do
      end do
      anti_t = 0.0_wp
      flux_t = 0.0_wp
      do b = 1, nb*nt
         do i = 1, ni
            im = ni + 1 - i
            anti_t = max(anti_t, abs(globt((b - 1)*ni + i) + globt((b - 1)*ni + im)))
            flux_t = max(flux_t, abs(globt((b - 1)*ni + i)))
         end do
      end do
      if (ctrl < 0) anti_max = -1.0_wp
   end subroutine fold_row_stats

   ! ------------------------------------------------------------------
   ! comparison
   ! ------------------------------------------------------------------

   function same(a, b) result(ok)
      real(wp), intent(in) :: a, b
      logical :: ok
      ok = ieee_is_finite(a) .and. ieee_is_finite(b) .and. a == b
   end function same

   subroutine compare(dec, ref, nbad, first)
      !! Bitwise comparison of every owned cell (the reference window at the tile's global offsets) and of
      !! the step results.  `first` describes the first mismatch.
      type(snap_t), intent(in) :: dec, ref
      integer, intent(out) :: nbad
      character(len=*), intent(out) :: first
      integer :: i, j, n, k, gi, gj
      logical :: have

      nbad = 0
      first = ""
      have = .false.
      do n = 1, NCAT
         do j = 1, dec%nyl
            do i = 1, dec%nxl
               gi = dec%io + i
               gj = dec%jo + j
               call chk("aicen", dec%aicen(i, j, n), ref%aicen(gi, gj, n), i, j, nbad, have, first)
               call chk("vicen", dec%vicen(i, j, n), ref%vicen(gi, gj, n), i, j, nbad, have, first)
               call chk("vsnon", dec%vsnon(i, j, n), ref%vsnon(gi, gj, n), i, j, nbad, have, first)
               call chk("Tsfcn", dec%Tsfcn(i, j, n), ref%Tsfcn(gi, gj, n), i, j, nbad, have, first)
               call chk("qsnon", dec%qsnon(i, j, n, 1), ref%qsnon(gi, gj, n, 1), i, j, nbad, have, first)
               do k = 1, NIL
                  call chk("qicen", dec%qicen(i, j, n, k), ref%qicen(gi, gj, n, k), i, j, nbad, have, first)
                  call chk("sicen", dec%sicen(i, j, n, k), ref%sicen(gi, gj, n, k), i, j, nbad, have, first)
               end do
            end do
         end do
      end do
      do j = 1, dec%nyl
         do i = 1, dec%nxl
            call chk("aice0", dec%aice0(i, j), ref%aice0(dec%io + i, dec%jo + j), i, j, nbad, have, first)
         end do
      end do
      if (dec%ierr /= ref%ierr .or. dec%cfg /= ref%cfg .or. dec%nbad_steps /= ref%nbad_steps) then
         nbad = nbad + 1
         if (.not. have) write (first, '(a)') "step result (ierr/cfg/stepped) differs"
         have = .true.
      end if
      if (.not. same(dec%worst_rel, ref%worst_rel) .or. dec%worst_q /= ref%worst_q) then
         nbad = nbad + 1
         if (.not. have) write (first, '(a,es12.4,a,es12.4)') "conservation EFP result differs: ", dec%worst_rel, " vs ", &
            ref%worst_rel
         have = .true.
      end if
      if (dec%anti_t /= ref%anti_t .and. dec%anti_t >= 0.0_wp .and. ref%anti_t >= 0.0_wp) then
         nbad = nbad + 1
         if (.not. have) write (first, '(a,es12.4,a,es12.4)') "fold-row tracer flux antisymmetry differs: ", dec%anti_t, " vs ", &
            ref%anti_t
         have = .true.
      end if
      if (dec%anti_max /= ref%anti_max .and. dec%anti_max >= 0.0_wp .and. ref%anti_max >= 0.0_wp) then
         nbad = nbad + 1
         if (.not. have) write (first, '(a,es12.4,a,es12.4)') "fold-row flux antisymmetry differs: ", dec%anti_max, " vs ", &
            ref%anti_max
         have = .true.
      end if

   end subroutine compare

   subroutine chk(tag, a, b, i, j, nbad, have, first)
      !! One owned-value comparison; the first mismatch is described in `first`.
      character(len=*), intent(in) :: tag
      real(wp), intent(in) :: a, b
      integer, intent(in) :: i, j
      integer, intent(inout) :: nbad
      logical, intent(inout) :: have
      character(len=*), intent(inout) :: first

      if (.not. same(a, b)) then
         nbad = nbad + 1
         if (.not. have) then
            write (first, '(a,a,i0,a,i0,a,es22.15,a,es22.15)') trim(tag), " owned (", i, ",", j, ") ", a, " vs ", b
            have = .true.
         end if
      end if
   end subroutine chk

   ! ------------------------------------------------------------------
   ! cases
   ! ------------------------------------------------------------------

   subroutine run_case(label, with_controls)
      !! Reference + every factorisation; `with_controls` instead runs the negative controls.
      character(len=*), intent(in) :: label
      logical, intent(in) :: with_controls
      type(snap_t) :: ref, dec
      integer :: px, nbad, ctrl, nglob, ncells_w, nseam, nfold
      character(len=160) :: tag, first
      logical :: fold_case, upwind, fault
      real(wp) :: rel
      integer :: i, j

      if (len_trim(case_filter) > 0) then
         if (index(label, trim(case_filter)) == 0) return
      end if
      fold_case = index(label, "fold") == 1
      upwind = label == "fold_upwind"
      fault = index(label, "_nan") > 0 .or. index(label, "_cfl") > 0
      call run_one(label, 1, 1, 1, 0, 0, ref)
      if (.not. ref%ok) then
         write (*, '(3a,i0)') "FAIL ", label, ": serial reference failed on rank ", rank
         n_fail = n_fail + 1
         return
      end if
      if (fault) then
         ! fault injection on one rank: the reference (one rank) and every layout must refuse the step with the same word
         if (ref%nbad_steps /= 1 .or. ref%ierr == 0) then
            write (*, '(3a)') "FAIL ", label, ": the injected fault did not refuse the reference step"
            n_fail = n_fail + 1
            return
         end if
         do px = nprocs, 1, -1
            if (mod(nprocs, px) /= 0) cycle
            write (tag, '(a,1x,i0,a,i0)') trim(label), px, "x", nprocs/px
            call run_one(label, px, nprocs/px, nprocs, rank, 0, dec)
            call compare(dec, ref, nbad, first)
            nseam = nbad
            call global_isum(nseam)
            if (.not. dec%ok .or. nseam > 0 .or. dec%nbad_steps /= 1) then
         write (*, '(5a,i0,a,z0,a,z0)') "FAIL ", trim(tag), " rank ", "", ": ", nseam, " ierr dec 0x", dec%ierr, " ref 0x", ref%ierr
               n_fail = n_fail + 1
            else if (rank == 0) then
               write (*, '(3a,z0)') "ok   ", trim(tag), ": every rank refused the step with the serial word 0x", ref%ierr
            end if
         end do
         return
      end if
      if (label == "fold_upwind") then
         ! upwind across the fold is refused (its fold-row face flux is not paired): serial and every layout alike
         if (ref%cfg /= TRANSPORT_CFG_FOLD .or. ref%nbad_steps /= 1) then
            write (*, '(3a,i0)') "FAIL ", label, ": the serial run was not refused with TRANSPORT_CFG_FOLD, cfg ", ref%cfg
            n_fail = n_fail + 1
            return
         end if
         do px = nprocs, 1, -1
            if (mod(nprocs, px) /= 0) cycle
            write (tag, '(a,1x,i0,a,i0)') trim(label), px, "x", nprocs/px
            call run_one(label, px, nprocs/px, nprocs, rank, 0, dec)
            ! every rank refuses (the verdict is reduced); the ranks holding the fold row name TRANSPORT_CFG_FOLD, the
            ! others a generic configuration error
            nseam = merge(0, 1, dec%cfg /= TRANSPORT_CFG_OK .and. dec%nbad_steps == 1)
            nfold = merge(1, 0, dec%cfg == TRANSPORT_CFG_FOLD)
            call global_isum(nseam)
            call global_isum(nfold)
            if (nseam > 0 .or. nfold == 0) then
               write (*, '(3a)') "FAIL ", trim(tag), ": not every rank refused, or no rank named TRANSPORT_CFG_FOLD"
               n_fail = n_fail + 1
            else if (rank == 0) then
               write (*, '(3a)') "ok   ", trim(tag), ": every rank refused upwind across the fold (TRANSPORT_CFG_FOLD)"
            end if
         end do
         return
      end if
      ! the reference itself: accepted, finite, conserving to round-off, and the ice really moved
      if (ref%ierr /= 0 .or. ref%cfg /= TRANSPORT_CFG_OK .or. ref%nbad_steps /= 0 .or. ref%n_err_cells /= 0) then
         write (*, '(3a,z0,a,i0,a,i0,a,i0)') "FAIL ", label, ": reference step refused: ierr 0x", ref%ierr, " cfg ", ref%cfg, &
            " unstepped ", ref%nbad_steps, " cells ", ref%n_err_cells
         n_fail = n_fail + 1
         return
      end if
      rel = abs(ref%mass1 - ref%mass0)/max(abs(ref%mass0), tiny(1.0_wp))
      if (.not. (rel < 1.0e-12_wp .and. ref%mass0 > 0.0_wp)) then
         write (*, '(3a,es10.3)') "FAIL ", label, ": reference ice area not conserved, relative change ", rel
         n_fail = n_fail + 1
      end if
      if (fold_case .and. .not. with_controls) then
         ! non-vacuity: the west half of the top rows started empty and the blob came through the fold row
         ncells_w = 0
         do j = ref%nyl - 3, ref%nyl
            do i = 1, ref%nxl/2
               if (ref%aicen(i, j, 1) > 1.0e-6_wp) ncells_w = ncells_w + 1
            end do
         end do
         if (ncells_w == 0) then
            write (*, '(3a)') "FAIL ", label, ": no ice crossed the fold into the west half (vacuous test)"
            n_fail = n_fail + 1
         end if
         if (.not. upwind .and. ref%flux_max >= 0.0_wp) then
            if (.not. (ref%flux_max > 0.0_wp)) then
               write (*, '(3a)') "FAIL ", label, ": the reference fold-row flux is zero (vacuous test)"
               n_fail = n_fail + 1
            end if
            if (.not. (ref%flux_t > 0.0_wp)) then
               write (*, '(3a)') "FAIL ", label, ": the reference fold-row tracer flux is zero (vacuous test)"
               n_fail = n_fail + 1
            end if
            if (ref%anti_max /= 0.0_wp .or. ref%anti_t /= 0.0_wp) then
               write (*, '(3a,es10.3)') "FAIL ", label, ": reference fold row not exactly antisymmetric: ", ref%anti_max
               n_fail = n_fail + 1
            end if
         end if
      end if

      if (with_controls) then
         call run_controls(label, ref)
         return
      end if

      nglob = 0
      do px = nprocs, 1, -1
         if (mod(nprocs, px) /= 0) cycle
         write (tag, '(a,1x,i0,a,i0)') trim(label), px, "x", nprocs/px
         call run_one(label, px, nprocs/px, nprocs, rank, 0, dec)
         if (.not. dec%ok) then
            write (*, '(3a,i0)') "FAIL ", trim(tag), ": decomposed run failed on rank ", rank
            n_fail = n_fail + 1
            cycle
         end if
         call compare(dec, ref, nbad, first)
         if (abs(dec%mass1 - dec%mass0)/max(abs(dec%mass0), tiny(1.0_wp)) > 1.0e-12_wp) then
            write (*, '(3a,i0)') "FAIL ", trim(tag), ": decomposed ice area not conserved, rank ", rank
            nbad = nbad + 1
         end if
         nseam = 0
         if (dec%io > 0) then
            if (.not. any(ref%aicen(dec%io:dec%io + 1, :, :) > 1.0e-6_wp)) nseam = nseam + 1
         end if
         if (dec%jo > 0) then
            if (.not. any(ref%aicen(:, dec%jo:dec%jo + 1, :) > 1.0e-6_wp)) nseam = nseam + 1
         end if
         call global_isum(nseam)
         if (nseam > 0) then
            write (*, '(3a,i0,a)') "FAIL ", trim(tag), ": ", nseam, " rank seam(s) with no ice in the reference (vacuous seam)"
            nbad = nbad + 1
         end if
         if (nbad > 0) then
            write (*, '(5a,i0,2a)') "FAIL ", trim(tag), " rank ", "", ": ", nbad, " mismatch(es); first: ", trim(first)
            n_fail = n_fail + 1
         else if (rank == 0) then
            write (*, '(3a)') "ok   ", trim(tag), ": owned cells bitwise equal to the serial reference"
         end if
      end do
   end subroutine run_case

   subroutine run_controls(label, ref)
      !! Each deliberately broken run must differ from the reference (fold case, decomposed only).
      character(len=*), intent(in) :: label
      type(snap_t), intent(in) :: ref
      type(snap_t) :: dec
      integer :: ctrl, nbad, px, nbad_all
      character(len=160) :: first
      character(len=40) :: what
      logical :: want

      do ctrl = 1, 5
         select case (ctrl)
         case (1)
            what = "category state fold skipped"
         case (2)
            what = "face velocity fold skipped"
         case (3)
            what = "seam reductions off"
         case (4)
            what = "flux pairing off"
         case default
            what = "tracer flux pairing off"
         end select
         ! the layout that exercises the control: any decomposed one (px = nprocs, py = 1 and, for the
         ! pairing, the serial run itself: it is a property of the fold row, not of the decomposition)
         if (ctrl >= 4) then
            call run_one(label, 1, 1, 1, 0, ctrl, dec)
            if (.not. dec%ok) then
               write (*, '(3a)') "FAIL control ", trim(what), ": run failed"
               n_fail = n_fail + 1
               cycle
            end if
            if (.not. (ctrl_defect(dec, ctrl) > 0.0_wp)) then
               write (*, '(3a,es10.3)') "FAIL control ", trim(what), ": the unpaired fold row is already exactly antisymmetric: ", &
                  ctrl_defect(dec, ctrl)
               n_fail = n_fail + 1
            else if (rank == 0) then
            write (*, '(3a,es10.3)') "ok   control ", trim(what), ": unpaired fold row antisymmetry defect ", ctrl_defect(dec, ctrl)
            end if
            if (nprocs >= 2) then
               ! the distributed pairing too: with it off, the fold row of the px = nprocs run is not antisymmetric
               call run_one(label, nprocs, 1, nprocs, rank, ctrl, dec)
               nbad_all = 0
               if (dec%ok) then
                  if (.not. (ctrl_defect(dec, ctrl) > 0.0_wp)) nbad_all = 1
               else
                  nbad_all = 1
               end if
               call global_isum(nbad_all)
               if (nbad_all > 0) then
                  write (*, '(3a)') "FAIL control ", trim(what), ": the split fold row is already antisymmetric without pairing"
                  n_fail = n_fail + 1
               else if (rank == 0) then
                  write (*, '(3a,es10.3)') "ok   control ", trim(what), " (split fold row): defect ", ctrl_defect(dec, ctrl)
               end if
            end if
            cycle
         end if
         if (nprocs < 2) cycle
         px = nprocs
         call run_one(label, px, 1, nprocs, rank, ctrl, dec)
         if (.not. dec%ok) then
            write (*, '(3a)') "FAIL control ", trim(what), ": run failed"
            n_fail = n_fail + 1
            cycle
         end if
         call compare(dec, ref, nbad, first)
         nbad_all = nbad
         call global_isum(nbad_all)
         want = nbad_all > 0
         if (.not. want) then
            write (*, '(3a)') "FAIL control ", trim(what), ": the comparison did NOT detect the broken seam"
            n_fail = n_fail + 1
         else if (rank == 0) then
            write (*, '(3a,i0,a)') "ok   control ", trim(what), ": detected (", nbad_all, " mismatches over all ranks)"
         end if
      end do
   end subroutine run_controls

   function ctrl_defect(dec, ctrl) result(d)
      !! The antisymmetry defect a pairing control must produce: the tracer fluxes for control 5 (the area flux stays
      !! paired there and is required exactly antisymmetric), area and tracers otherwise.
      type(snap_t), intent(in) :: dec
      integer, intent(in) :: ctrl
      real(wp) :: d
      if (ctrl == 5) then
         d = dec%anti_t
         if (dec%anti_max /= 0.0_wp) d = -1.0_wp
      else
         d = min(dec%anti_max, dec%anti_t)
      end if
   end function ctrl_defect

   subroutine global_isum(n)
      integer, intent(inout) :: n
      type(comm_t) :: c

      c = comm_env_compute_comm()
      if (c%size() == 1) return
      call allreduce(c, n, MPI_SUM)
   end subroutine global_isum

end program test_cice_remap_decomp_mpi
#else
program test_cice_remap_decomp_mpi
   print *, "test_cice_remap_decomp_mpi: MPI build only"
end program test_cice_remap_decomp_mpi
#endif
