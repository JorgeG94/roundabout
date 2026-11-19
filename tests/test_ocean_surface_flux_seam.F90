!! Surface-flux ghost columns on a periodic seam, through the real engine.
!!
!! The boundary-layer schemes (KPP, EPBL) loop the FULL array and read the
!! assembled `Q_heat` / `Q_salt` (and `sw_net` under `sw_source = 'q_sw'`)
!! in the GHOST columns; the momentum vertical-viscosity face average then
!! reads the ghost-column `kv` at the first owned face.  The assembler sums
!! components, and the atmosphere-driven bulk consumer writes its components
!! (`q_sens`, `q_lat`, `q_sw`, `q_lw`, `evap`, ...) on PHYSICAL cells only, so
!! without a refresh the ghost `Q_heat` lacks the whole bulk part and the
!! seam column mixes as if it were unforced.  An MPI seam shows the same
!! thing as a decomposed run that differs from the serial one
!! (`bulk_flux_bl*` in `test_ocean_decomp_bitid_mpi`); on one rank the
!! periodic wrap is the seam, which is what this test exercises.
!!
!! Cases: KPP, EPBL, KPP with `sw_source = 'q_sw'`, EPBL with `sw_source = 'q_sw'`
!! and `epbl_sw_ctke` (the penetrating-shortwave TKE ledger reads `sw_net`),
!! KPP with a data-override SALT file on top of the bulk forcing (a non-zero
!! `Q_salt` exchange), and a components-OFF data-override HEAT file (the
!! reader exchanges the `Q_heat` it writes, no bulk consumer).  Plus
!! `api_setters`: `rdb_ocean_set_heat_flux`/`set_salt_flux`/`set_wind` write
!! physical cells only and must leave valid ghosts.
!!
!! Two gates per case:
!!   * `ghost_flux_is_the_wrapped_interior` -- after the steps, the west
!!     ghost columns of `Q_heat` and `Q_salt` equal the east physical
!!     columns and vice versa, over the physical rows.
!!     Under `sw_source = 'q_sw'` the same holds for `sw_net`.
!!   * `shift_invariant` -- the same ocean on a domain rolled by half a period
!!     (the atmosphere and the data files have zonal wavenumber 2, so the roll
!!     is exact) gives the same `kv`, `u`, `v`, `h`, `T`, `S`, `eta` cell for
!!     cell: the answer cannot depend on where the seam falls.
!!
!! NON-VACUITY (the way this bug hid): the cold dry wind over 20 C water
!! cools the surface so hard that the boundary layer is deeper than the top
!! interface, so the overlay's `kv` at the first interface below the surface
!! is orders above the interior value in every column, and the ghost-column
!! flux error (~1e2 W/m2) is far above its tolerance.  Both are asserted.
!!
!! Own ctest binary: `rdb_ocean_api` keeps one live handle per process.
module test_ocean_surface_flux_seam
   use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double, c_null_ptr, c_f_pointer
   use rdb_constants, only: wp
   use rdb_ocean_api, only: rdb_ocean_create_from_string, rdb_ocean_step, &
                            rdb_ocean_destroy, rdb_ocean_refresh_host, &
                            rdb_ocean_get_grid_info, rdb_ocean_get_bt_eta_ptr, &
                            rdb_ocean_set_heat_flux, rdb_ocean_set_salt_flux, &
                            rdb_ocean_set_wind, rdb_ocean_get_tau_x_ptr, &
                            rdb_ocean_get_tau_y_ptr, &
                            rdb_ocean_get_u_face_x_layer_ptr, &
                            rdb_ocean_get_v_face_y_layer_ptr, rdb_ocean_get_h_layer_ptr, &
                            rdb_ocean_get_q_heat_ptr, rdb_ocean_get_q_salt_ptr, &
                            rdb_ocean_get_kv_ptr
   use rdb_handle, only: ocean_handle_t, handle_check, HANDLE_OK
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   use rdb_io_netcdf, only: nc_create_file, nc_close, nc_def_dim, nc_def_var_2d, &
                            nc_def_var_3d, nc_enddef, rdb_def_var_1d, rdb_put_var_1d, nc_put_att
   use netcdf, only: nf90_put_var
   use testdrive, only: error_type, check, new_unittest, unittest_type
   implicit none
   private

   public :: collect_ocean_surface_flux_seam_tests

   integer, parameter :: NX = 16
      !! Physical cells in x (periodic).  Even, so the half-period roll is exact.
   integer, parameter :: NY = 8
   integer, parameter :: NZ = 10
   integer, parameter :: SHIFT = NX/2
   integer, parameter :: N_STEPS = 6
   real(wp), parameter :: REL_TOL = 1.0e-10_wp
   character(len=*), parameter :: ATM_FILE = "tmp_local_artifacts/seam_flux_atm.nc"
   character(len=*), parameter :: ATM_LIST = "tmp_local_artifacts/seam_flux_atm_files.txt"
   character(len=*), parameter :: HEAT_FILE = "tmp_local_artifacts/seam_flux_heat.nc"
   character(len=*), parameter :: SALT_FILE = "tmp_local_artifacts/seam_flux_salt.nc"
   character(len=1), parameter :: NL = new_line("a")

   type :: run_out_t
      real(wp), allocatable :: qh(:, :), qs(:, :)
         !! `Q_heat` / `Q_salt`, FULL extent (ghosts included).
      real(wp), allocatable :: sw(:, :)
         !! `sw_net`, FULL extent (zero-size when the case has no `q_sw` source).
      real(wp), allocatable :: kv(:, :, :), u(:, :, :), v(:, :, :), h(:, :, :)
      real(wp), allocatable :: tt(:, :, :), ss(:, :, :), eta(:, :)
         !! Physical interior (u: faces 1..NX; v: all NY+1 rows; kv: NZ+1 interfaces).
      integer :: ng = 0
   end type run_out_t

contains

   subroutine collect_ocean_surface_flux_seam_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("seam_flux_ghost_kpp", test_ghost_kpp), &
                  new_unittest("seam_flux_ghost_epbl", test_ghost_epbl), &
                  new_unittest("seam_flux_ghost_kpp_sw", test_ghost_kpp_sw), &
                  new_unittest("seam_flux_shift_invariant_kpp", test_shift_kpp), &
                  new_unittest("seam_flux_shift_invariant_epbl", test_shift_epbl), &
                  new_unittest("seam_flux_shift_invariant_kpp_sw", test_shift_kpp_sw), &
                  new_unittest("seam_flux_ghost_epbl_sw", test_ghost_epbl_sw), &
                  new_unittest("seam_flux_shift_invariant_epbl_sw", test_shift_epbl_sw), &
                  new_unittest("seam_flux_ghost_kpp_salt", test_ghost_kpp_salt), &
                  new_unittest("seam_flux_shift_invariant_kpp_salt", test_shift_kpp_salt), &
                  new_unittest("seam_flux_ghost_dataovr_heat", test_ghost_dataovr_heat), &
                  new_unittest("seam_flux_shift_invariant_dataovr_heat", test_shift_dataovr_heat), &
                  new_unittest("seam_flux_api_setters", test_api_setters), &
                  new_unittest("seam_flux_ghost_cavity_periodic", test_ghost_cavity) &
                  ]
   end subroutine collect_ocean_surface_flux_seam_tests

   function case_nml(kind, lon_west) result(nml)
      !! Re-entrant 100 m basin on ten 10 m sigma layers.  `kind`:
      !!   "kpp" / "epbl"      bulk forcing (cold dry atmosphere file) + that scheme;
      !!   "kpp_sw"/"epbl_sw"  same with `sw_source = 'q_sw'`, `sw_pen_frac = 0.6`
      !!                       (EPBL also `epbl_sw_ctke`);
      !!   "kpp_salt"          "kpp" plus a data-override SALT file;
      !!   "dataovr_heat"      NO bulk forcing, components OFF, KPP, a data-override
      !!                       HEAT file written straight into `Q_heat`;
      !!   "api"               components OFF, nothing driving (API-setter test).
      character(len=*), intent(in) :: kind
      real(wp), intent(in) :: lon_west
      character(len=:), allocatable :: nml
      character(len=32) :: slon
      character(len=:), allocatable :: thermo, bl, forcing, ovr, wind
      logical :: bulk, is_epbl, is_sw

      write (slon, "(f8.2)") lon_west
      bulk = (kind /= "dataovr_heat" .and. kind /= "api")
      is_epbl = (index(kind, "epbl") == 1)
      is_sw = (index(kind, "_sw") > 0)
      thermo = "&ocean_thermo_nml enable_thermodynamics = .true. /"//NL
      if (is_sw) thermo = "&ocean_thermo_nml enable_thermodynamics = .true., "// &
                          "sw_pen_frac = 0.6, sw_source = 'q_sw' /"//NL
      if (is_sw .and. is_epbl) thermo = "&ocean_thermo_nml enable_thermodynamics = .true., "// &
                                        "sw_pen_frac = 0.6, sw_source = 'q_sw', "// &
                                        "epbl_sw_ctke = .true. /"//NL
      bl = ""
      if (is_epbl) bl = "&ocean_vmix_nml use_kpp = .false. /"//NL// &
                        "&ocean_epbl_nml enable = .true. /"//NL
      forcing = ""
      if (bulk) forcing = "&ocean_forcing_nml enable_components = .true. /"//NL// &
                          "&ocean_bulk_flux_nml enable = .true., n_iter = 5, z_wind = 10.0, "// &
                          "z_ta = 10.0 /"//NL// &
                          "&ocean_atm_forcing_nml enable = .true., time_mode = 'interannual', "// &
                          "sw_interp = 'linear', "// &
                          "filelist_tas = '"//ATM_LIST//"', filelist_huss = '"//ATM_LIST// &
                          "', filelist_uas = '"//ATM_LIST//"', filelist_vas = '"//ATM_LIST// &
                          "', filelist_psl = '"//ATM_LIST//"', filelist_rlds = '"//ATM_LIST// &
                          "', filelist_rsds = '"//ATM_LIST//"', filelist_prra = '"//ATM_LIST// &
                          "', filelist_prsn = '"//ATM_LIST//"' /"//NL
      wind = ""
      if (kind == "dataovr_heat") wind = "&physics_nml wind_stress_x = 0.2, wind_stress_y = 0.05 /"//NL
      ovr = ""
      if (kind == "kpp_salt") ovr = "&ocean_dataovr_nml enable = .true., time_mode = 'cyclic', "// &
                                    "cycle_period = 1.0e6, salt_file = '"//SALT_FILE// &
                                    "', salt_var = 'field' /"//NL
      if (kind == "dataovr_heat") ovr = "&ocean_dataovr_nml enable = .true., time_mode = 'cyclic', "// &
                                        "cycle_period = 1.0e6, heat_file = '"//HEAT_FILE// &
                                        "', heat_var = 'field' /"//NL
      nml = "&sim_nml sim_type = 'ocean' /"//NL// &
            "&time_nml t_end = 1.0e9, dt_fixed = 900.0, start_date = '2000-06-21' /"//NL// &
            "&nonhydrostatic_nml nz_layers = 10 /"//NL// &
            "&ocean_bt_nml auto_n_inner = .true. /"//NL// &
            "&ocean_hvisc_nml nu_h = 200.0, lateral_closure = 'smagorinsky', "// &
            "smag_ah = .true. /"//NL// &
            "&ocean_diag_nml enabled = .false. /"//NL// &
            "&ocean_eos_nml eos = 'wright' /"//NL// &
            thermo//forcing//ovr//bl//wind// &
            "&tracer_nml initial_temperature = 12.0, initial_salinity = 35.0, "// &
            "T_init_surface = 20.0, T_init_bottom = 19.0 /"//NL// &
            "&grid_nml nx = 16, ny = 8, nghost = 3, dx = 22.5, dy = 1.0 /"//NL// &
            "&ocean_grid_nml grid_config = 'spherical', lon_west = "//trim(adjustl(slon))// &
            ", lat_south = 30.0, rad_earth = 6.371e6, coriolis_scheme = 'planetary' /"//NL// &
            "&vcoord_nml vcoord_type = 'sigma' /"//NL// &
            "&ocean_topo_nml max_depth = 100.0 /"//NL// &
            "&ocean_bc_nml west = 'periodic', east = 'periodic', south = 'wall', "// &
            "north = 'wall' /"//NL// &
            "&output_nml output_to_file = .false. /"//NL
   end function case_nml

   subroutine run_case(error, kind, lon_west, out)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: kind
      real(wp), intent(in) :: lon_west
      type(run_out_t), intent(out) :: out
      type(c_ptr) :: handle, ptr
      integer(c_int) :: status, nx_p, ny_p, nz_p, ng, nxt, nyt, nzt, gen
      character(len=:), allocatable :: nml
      real(wp), pointer :: p2(:, :), p3(:, :, :)
      type(ocean_handle_t), pointer :: h
      integer :: it, nzm

      call write_flux_files()
      nml = case_nml(kind, lon_west)
      handle = c_null_ptr
      status = rdb_ocean_create_from_string(nml, len(nml, kind=c_int), handle)
      call check(error, status == OCEAN_STATUS_OK, kind//": the namelist must build")
      if (allocated(error)) return
      status = rdb_ocean_get_grid_info(handle, nx_p, ny_p, nz_p, ng)
      status = rdb_ocean_step(handle, int(N_STEPS, c_int))
      call check(error, status == OCEAN_STATUS_OK, kind//": the run must step cleanly")
      if (allocated(error)) return
      status = rdb_ocean_refresh_host(handle)
      out%ng = ng
      status = rdb_ocean_get_q_heat_ptr(handle, ptr, nxt, nyt, gen)
      call c_f_pointer(ptr, p2, [nxt, nyt])
      out%qh = p2
      status = rdb_ocean_get_q_salt_ptr(handle, ptr, nxt, nyt, gen)
      call c_f_pointer(ptr, p2, [nxt, nyt])
      out%qs = p2
      status = rdb_ocean_get_kv_ptr(handle, ptr, nxt, nyt, nzt, gen)
      call c_f_pointer(ptr, p3, [nxt, nyt, nzt])
      out%kv = p3(ng + 1:ng + NX, ng + 1:ng + NY, 1:NZ + 1)
      status = rdb_ocean_get_u_face_x_layer_ptr(handle, ptr, nxt, nyt, nzt, gen)
      call c_f_pointer(ptr, p3, [nxt, nyt, nzt])
      out%u = p3(ng + 1:ng + NX, ng + 1:ng + NY, 1:NZ)
      status = rdb_ocean_get_v_face_y_layer_ptr(handle, ptr, nxt, nyt, nzt, gen)
      call c_f_pointer(ptr, p3, [nxt, nyt, nzt])
      out%v = p3(ng + 1:ng + NX, ng + 1:ng + NY + 1, 1:NZ)
      status = rdb_ocean_get_h_layer_ptr(handle, ptr, nxt, nyt, nzt, gen)
      call c_f_pointer(ptr, p3, [nxt, nyt, nzt])
      out%h = p3(ng + 1:ng + NX, ng + 1:ng + NY, 1:NZ)
      status = rdb_ocean_get_bt_eta_ptr(handle, ptr, nxt, nyt, gen)
      call c_f_pointer(ptr, p2, [nxt, nyt])
      out%eta = p2(ng + 1:ng + NX, ng + 1:ng + NY)
      ! Tracers and `sw_net` are not API getters: read them through the handle.
      if (handle_check(handle, h) /= HANDLE_OK) then
         call check(error, .false., kind//": handle lost")
         return
      end if
      nzm = h%state%multilayer%nz_ml
      it = h%state%multilayer%idx_temperature
      !$omp target update from(h%state%multilayer%tracers(it)%c)
      out%tt = h%state%multilayer%tracers(it)%c(ng + 1:ng + NX, ng + 1:ng + NY, 1:nzm)
      it = h%state%multilayer%idx_salinity
      !$omp target update from(h%state%multilayer%tracers(it)%c)
      out%ss = h%state%multilayer%tracers(it)%c(ng + 1:ng + NX, ng + 1:ng + NY, 1:nzm)
      allocate (out%sw(0, 0))
      if (h%state%surface_flux%sw_from_qsw) then
         !$omp target update from(h%state%surface_flux%sw_net)
         out%sw = h%state%surface_flux%sw_net
      end if
      status = rdb_ocean_destroy(handle)
   end subroutine run_case

   subroutine check_overlay_active(error, o, what)
      !! Non-vacuity: the boundary layer must reach below the first interface
      !! below the surface in the seam columns, and the surface must be
      !! forced hard.
      type(error_type), allocatable, intent(out) :: error
      type(run_out_t), intent(in) :: o
      character(len=*), intent(in) :: what
      real(wp) :: kv_seam, kv_mid, qmax
      integer :: ng
      character(len=200) :: msg

      ng = o%ng
      qmax = maxval(abs(o%qh(ng + 1:ng + NX, ng + 1:ng + NY)))
      call check(error, qmax > 50.0_wp, what//": the bulk heat flux is not strong "// &
                 "(max|Q_heat| < 50 W/m2) -- the case is vacuous")
      if (allocated(error)) return
      ! kv interface NZ = the first interface below the surface (NZ+1).
      kv_seam = minval(o%kv(1, 2:NY - 1, NZ))
      kv_mid = maxval(o%kv(NX/2, 2:NY - 1, 1:NZ - 3))
      write (msg, "(a,es10.3,a,es10.3,a,es10.3)") " kv_seam ", kv_seam, " kv_mid ", kv_mid, &
         " maxQ ", qmax
      call check(error, kv_seam > 1.0e-3_wp .and. kv_seam > 10.0_wp*kv_mid, &
                 what//trim(msg)//": the boundary-layer overlay is dormant at the seam column "// &
                 "(kv at the first interface too small) -- the case is vacuous")
   end subroutine check_overlay_active

   subroutine check_ghosts(error, o, what, want_salt)
      type(error_type), allocatable, intent(out) :: error
      type(run_out_t), intent(in) :: o
      character(len=*), intent(in) :: what
      logical, intent(in) :: want_salt
      real(wp) :: dh, ds, dw, qmax, smax
      integer :: ng, g, j
      character(len=240) :: msg

      call check_overlay_active(error, o, what)
      if (allocated(error)) return
      ng = o%ng
      dh = 0.0_wp
      ds = 0.0_wp
      do j = ng + 1, ng + NY
         do g = 1, ng
            dh = max(dh, abs(o%qh(g, j) - o%qh(g + NX, j)), &
                     abs(o%qh(ng + NX + g, j) - o%qh(ng + g, j)))
            ds = max(ds, abs(o%qs(g, j) - o%qs(g + NX, j)), &
                     abs(o%qs(ng + NX + g, j) - o%qs(ng + g, j)))
         end do
      end do
      qmax = maxval(abs(o%qh(ng + 1:ng + NX, ng + 1:ng + NY)))
      smax = maxval(abs(o%qs(ng + 1:ng + NX, ng + 1:ng + NY)))
      dw = 0.0_wp
      if (size(o%sw, 1) > 0) then
         do j = ng + 1, ng + NY
            do g = 1, ng
               dw = max(dw, abs(o%sw(g, j) - o%sw(g + NX, j)), &
                        abs(o%sw(ng + NX + g, j) - o%sw(ng + g, j)))
            end do
         end do
         call check(error, maxval(abs(o%sw(ng + 1:ng + NX, ng + 1:ng + NY))) > 10.0_wp, &
                    what//": sw_net is not forced (max < 10 W/m2) -- vacuous")
         if (allocated(error)) return
      end if
      if (want_salt) then
         call check(error, smax > 1.0e-6_wp, what//": Q_salt is ~0 -- the salt gate is vacuous")
         if (allocated(error)) return
      end if
      write (msg, "(a,es10.3,a,es10.3,a,es10.3,a,es10.3,a,es10.3)") ": |dQ_heat| = ", dh, &
         " (max|Q_heat| ", qmax, ")  |dQ_salt| = ", ds, " (max|Q_salt| ", smax, ")  |dsw_net| = ", dw
      call check(error, dh <= 1.0e-12_wp*qmax .and. ds <= 1.0e-12_wp*max(smax, tiny(1.0_wp)) .and. &
                 dw <= 1.0e-12_wp*max(1.0_wp, maxval(abs(o%sw))), &
                 what//trim(msg)//" -- the ghost surface flux is not the wrapped "// &
                 "interior: the boundary layer reads a stale ghost flux at the seam")
   end subroutine check_ghosts

   subroutine check_shift(error, a, b, what)
      !! Run B is rolled by `SHIFT`: B(i + SHIFT) must equal A(i).
      type(error_type), allocatable, intent(out) :: error
      type(run_out_t), intent(in) :: a, b
      character(len=*), intent(in) :: what
      real(wp) :: umax, kmax, du, dk, dh, dv, dt, ds, de
      integer :: i, ib
      character(len=300) :: msg

      call check_overlay_active(error, a, what)
      if (allocated(error)) return
      umax = max(maxval(abs(a%u)), maxval(abs(a%v)))
      kmax = maxval(a%kv)
      call check(error, umax > 1.0e-6_wp, what//": the ocean did not move -- vacuous")
      if (allocated(error)) return
      du = 0.0_wp
      dv = 0.0_wp
      dk = 0.0_wp
      dh = 0.0_wp
      dt = 0.0_wp
      ds = 0.0_wp
      de = 0.0_wp
      do i = 1, NX
         ib = modulo(i - 1 + SHIFT, NX) + 1
         dt = max(dt, maxval(abs(a%tt(i, :, :) - b%tt(ib, :, :))))
         ds = max(ds, maxval(abs(a%ss(i, :, :) - b%ss(ib, :, :))))
         de = max(de, maxval(abs(a%eta(i, :) - b%eta(ib, :))))
         du = max(du, maxval(abs(a%u(i, :, :) - b%u(ib, :, :))))
         dv = max(dv, maxval(abs(a%v(i, :, :) - b%v(ib, :, :))))
         dk = max(dk, maxval(abs(a%kv(i, :, :) - b%kv(ib, :, :))))
         dh = max(dh, maxval(abs(a%h(i, :, :) - b%h(ib, :, :))))
      end do
      write (msg, "(a,es10.3,a,es10.3,a,es10.3,a,es10.3,a,es10.3,a,3es10.3)") ": max|u| = ", umax, &
         "  |du| = ", du, "  |dv| = ", dv, "  |dkv| = ", dk, "  |dh| = ", dh, &
         "  |dT,dS,deta| = ", dt, ds, de
      call check(error, max(du, dv) <= REL_TOL*umax .and. dk <= REL_TOL*kmax .and. &
                 dh <= REL_TOL*100.0_wp .and. dt <= REL_TOL*20.0_wp .and. &
                 ds <= REL_TOL*35.0_wp .and. de <= REL_TOL*max(1.0_wp, maxval(abs(a%eta))), &
                 what//trim(msg)//" -- the run is not "// &
                 "shift-invariant; the periodic seam saw a stale ghost flux")
   end subroutine check_shift

   subroutine ghost_case(error, kind)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: kind
      type(run_out_t) :: a

      call run_case(error, kind, 0.0_wp, a)
      if (allocated(error)) return
      call check_ghosts(error, a, kind, kind == "kpp_salt")
   end subroutine ghost_case

   subroutine shift_case(error, kind)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: kind
      type(run_out_t) :: a, b

      call run_case(error, kind, 0.0_wp, a)
      if (allocated(error)) return
      call run_case(error, kind, 180.0_wp, b)
      if (allocated(error)) return
      call check_shift(error, a, b, kind)
   end subroutine shift_case

   subroutine test_ghost_kpp(error)
      type(error_type), allocatable, intent(out) :: error
      call ghost_case(error, "kpp")
   end subroutine test_ghost_kpp

   subroutine test_ghost_epbl(error)
      type(error_type), allocatable, intent(out) :: error
      call ghost_case(error, "epbl")
   end subroutine test_ghost_epbl

   subroutine test_ghost_kpp_sw(error)
      type(error_type), allocatable, intent(out) :: error
      call ghost_case(error, "kpp_sw")
   end subroutine test_ghost_kpp_sw

   subroutine test_shift_kpp(error)
      type(error_type), allocatable, intent(out) :: error
      call shift_case(error, "kpp")
   end subroutine test_shift_kpp

   subroutine test_shift_epbl(error)
      type(error_type), allocatable, intent(out) :: error
      call shift_case(error, "epbl")
   end subroutine test_shift_epbl

   subroutine test_shift_kpp_sw(error)
      type(error_type), allocatable, intent(out) :: error
      call shift_case(error, "kpp_sw")
   end subroutine test_shift_kpp_sw

   subroutine test_ghost_epbl_sw(error)
      type(error_type), allocatable, intent(out) :: error
      call ghost_case(error, "epbl_sw")
   end subroutine test_ghost_epbl_sw

   subroutine test_shift_epbl_sw(error)
      type(error_type), allocatable, intent(out) :: error
      call shift_case(error, "epbl_sw")
   end subroutine test_shift_epbl_sw

   subroutine test_ghost_kpp_salt(error)
      type(error_type), allocatable, intent(out) :: error
      call ghost_case(error, "kpp_salt")
   end subroutine test_ghost_kpp_salt

   subroutine test_shift_kpp_salt(error)
      type(error_type), allocatable, intent(out) :: error
      call shift_case(error, "kpp_salt")
   end subroutine test_shift_kpp_salt

   subroutine test_ghost_dataovr_heat(error)
      type(error_type), allocatable, intent(out) :: error
      call ghost_case(error, "dataovr_heat")
   end subroutine test_ghost_dataovr_heat

   subroutine test_shift_dataovr_heat(error)
      type(error_type), allocatable, intent(out) :: error
      call shift_case(error, "dataovr_heat")
   end subroutine test_shift_dataovr_heat

   subroutine test_api_setters(error)
      !! `rdb_ocean_set_heat_flux` / `set_salt_flux` / `set_wind` write the
      !! physical interior only; the ghost band must then be the wrapped
      !! interior (periodic x), and the wind's `stress_mag` re-derived.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle, ptr
      integer(c_int) :: status, nx_p, ny_p, nz_p, ng, nxt, nyt, gen
      character(len=:), allocatable :: nml
      real(c_double) :: q(NX, NY), tx(NX + 1, NY), ty(NX, NY + 1)
      real(wp), pointer :: p2(:, :)
      real(wp) :: dq, ds, dtx, dty
      integer :: i, j, g

      nml = case_nml("api", 0.0_wp)
      handle = c_null_ptr
      status = rdb_ocean_create_from_string(nml, len(nml, kind=c_int), handle)
      call check(error, status == OCEAN_STATUS_OK, "api: the namelist must build")
      if (allocated(error)) return
      status = rdb_ocean_get_grid_info(handle, nx_p, ny_p, nz_p, ng)
      do j = 1, NY
         do i = 1, NX
            q(i, j) = 100.0_wp + 10.0_wp*real(i, wp) + real(j, wp)
         end do
      end do
      do j = 1, NY
         do i = 1, NX + 1
            tx(i, j) = 0.01_wp*real(i, wp) + 0.001_wp*real(j, wp)
         end do
      end do
      do j = 1, NY + 1
         do i = 1, NX
            ty(i, j) = 0.02_wp*real(i, wp) - 0.001_wp*real(j, wp)
         end do
      end do
      status = rdb_ocean_set_heat_flux(handle, q, int(NX, c_int), int(NY, c_int))
      call check(error, status == OCEAN_STATUS_OK, "api: set_heat_flux")
      status = rdb_ocean_set_salt_flux(handle, 1.0e-4_wp*q, int(NX, c_int), int(NY, c_int))
      status = rdb_ocean_set_wind(handle, tx, ty, int(NX, c_int), int(NY, c_int))
      call check(error, status == OCEAN_STATUS_OK, "api: set_wind")
      if (allocated(error)) return
      status = rdb_ocean_refresh_host(handle)
      dq = 0.0_wp
      ds = 0.0_wp
      dtx = 0.0_wp
      dty = 0.0_wp
      status = rdb_ocean_get_q_heat_ptr(handle, ptr, nxt, nyt, gen)
      call c_f_pointer(ptr, p2, [nxt, nyt])
      do j = ng + 1, ng + NY
         do g = 1, ng
            dq = max(dq, abs(p2(g, j) - p2(g + NX, j)), abs(p2(ng + NX + g, j) - p2(ng + g, j)))
         end do
      end do
      call check(error, dq == 0.0_wp .and. maxval(abs(p2(ng + 1:ng + NX, ng + 1:ng + NY))) > 100.0_wp, &
                 "api: set_heat_flux left stale ghost Q_heat")
      if (allocated(error)) return
      status = rdb_ocean_get_q_salt_ptr(handle, ptr, nxt, nyt, gen)
      call c_f_pointer(ptr, p2, [nxt, nyt])
      do j = ng + 1, ng + NY
         do g = 1, ng
            ds = max(ds, abs(p2(g, j) - p2(g + NX, j)), abs(p2(ng + NX + g, j) - p2(ng + g, j)))
         end do
      end do
      call check(error, ds == 0.0_wp .and. maxval(abs(p2(ng + 1:ng + NX, ng + 1:ng + NY))) > 1.0e-3_wp, &
                 "api: set_salt_flux left stale ghost Q_salt")
      if (allocated(error)) return
      status = rdb_ocean_get_tau_x_ptr(handle, ptr, nxt, nyt, gen)
      call c_f_pointer(ptr, p2, [nxt, nyt])
      do j = ng + 1, ng + NY
         do g = 1, ng
            dtx = max(dtx, abs(p2(g, j) - p2(g + NX, j)), &
                      abs(p2(ng + NX + 1 + g, j) - p2(ng + 1 + g, j)))
         end do
      end do
      call check(error, dtx == 0.0_wp .and. maxval(abs(p2)) > 0.01_wp, &
                 "api: set_wind left stale ghost tau_x faces")
      if (allocated(error)) return
      status = rdb_ocean_get_tau_y_ptr(handle, ptr, nxt, nyt, gen)
      call c_f_pointer(ptr, p2, [nxt, nyt])
      do j = ng + 1, ng + NY + 1
         do g = 1, ng
            dty = max(dty, abs(p2(g, j) - p2(g + NX, j)), abs(p2(ng + NX + g, j) - p2(ng + g, j)))
         end do
      end do
      call check(error, dty == 0.0_wp .and. maxval(abs(p2)) > 0.01_wp, &
                 "api: set_wind left stale ghost tau_y faces")
      status = rdb_ocean_destroy(handle)
   end subroutine test_api_setters

   function cavity_nml() result(nml)
      !! The decomposition harness's `cavity_periodic` (under `ssp_rk2`, where the
      !! unexchanged ghost lags the owner by ~8e-2 W/m2; ~1e-7 under pred_corr): a 26 x 18, 20 km sigma
      !! basin, re-entrant in x, under a sloping ice-shelf draft with virtual
      !! melt, KPP live under the ice.  The shelf box sits against the wrap, so
      !! the cavity bands (`heat_cavity`/`salt_cavity`, written on physical
      !! cells only) differ across the seam unless the assembled flux is exchanged.
      character(len=:), allocatable :: nml
      nml = "&sim_nml sim_type = 'ocean' /"//NL// &
            "&time_nml t_end = 1.0e9, dt_fixed = 900.0 /"//NL// &
            "&nonhydrostatic_nml nz_layers = 4 /"//NL// &
            "&tracer_nml initial_salinity = 34.5, T_init_surface = 3.0, T_init_bottom = 1.0 /"//NL// &
            "&ocean_bt_nml auto_n_inner = .true., split_scheme = 'ssp_rk2' /"//NL// &
            "&ocean_hvisc_nml nu_h = 200.0, lateral_closure = 'smagorinsky', smag_ah = .true. /"//NL// &
            "&ocean_diag_nml enabled = .false. /"//NL// &
            "&grid_nml nx = 26, ny = 18, nghost = 3, dx = 20000.0, dy = 20000.0 /"//NL// &
            "&physics_nml coriolis_f = -1.409e-4 /"//NL// &
            "&vcoord_nml vcoord_type = 'sigma' /"//NL// &
            "&ocean_topo_nml topo_config = 'seamount', max_depth = 1000.0, edge_depth = 400.0, "// &
            "slope_scale = 80000.0, wind_config = 'constant', taux_magnitude = 0.05 /"//NL// &
            "&ocean_cavity_dyn_nml enable = .true., draft_config = 'linear', "// &
            "draft_depth = 500.0, draft_slope = -4.0e-4, draft_x0 = 0.0, "// &
            "draft_x1 = 380000.0, draft_y1 = 260000.0, h_min_cavity = 100.0 /"//NL// &
            "&ocean_forcing_nml enable_components = .true. /"//NL// &
            "&ocean_cavity_melt_nml enable = .true., exchange_law = 'const_gamma', "// &
            "gamma_t = 2.2e-2, gamma_s = -1.0, cdrag_top = 2.5e-3, "// &
            "ice_conduction = 'insulating', s_ice = 0.0, far_field_depth = 10.0 /"//NL// &
            "&ocean_tdrag_nml enable = .true., form = 'quadratic', cd = 2.5e-3 /"//NL// &
            "&ocean_pgf_nml form = 'fv_mom6', p_top_in_bc = .true. /"//NL// &
            "&ocean_eos_nml eos = 'linear', tfreeze_set = 'isomip' /"//NL// &
            "&ocean_bc_nml west = 'periodic', east = 'periodic', south = 'wall', "// &
            "north = 'wall' /"//NL// &
            "&output_nml output_to_file = .false. /"//NL
   end function cavity_nml

   subroutine test_ghost_cavity(error)
      !! Ice-shelf cavity, re-entrant in x: after the steps the ghost columns of
      !! `Q_heat`/`Q_salt` equal the wrapped interior.  Before the assembled-flux
      !! exchange the ghost differed by ~6e-2 W/m2 (3e-4 relative): the cavity
      !! bands are physical-span-only components and KPP reads the ghost flux.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NXC = 26, NYC = 18
      type(c_ptr) :: handle, ptr
      integer(c_int) :: status, nx_p, ny_p, nz_p, ng, nxt, nyt, gen
      character(len=:), allocatable :: nml
      real(wp), pointer :: p2(:, :)
      real(wp) :: dh, ds, qmax, smax, hcmax
      type(ocean_handle_t), pointer :: h
      integer :: j, g, istep
      character(len=200) :: msg

      nml = cavity_nml()
      handle = c_null_ptr
      status = rdb_ocean_create_from_string(nml, len(nml, kind=c_int), handle)
      call check(error, status == OCEAN_STATUS_OK, "cavity: the namelist must build")
      if (allocated(error)) return
      status = rdb_ocean_get_grid_info(handle, nx_p, ny_p, nz_p, ng)
      call check(error, handle_check(handle, h) == HANDLE_OK, "cavity: handle")
      if (allocated(error)) return
      ! The ghost lags the owner by one thermo step while the melt spins up, so
      ! the worst mismatch is early: check after EVERY one of the first steps.
      dh = 0.0_wp
      qmax = 0.0_wp
      hcmax = 0.0_wp
      do istep = 1, 24
         status = rdb_ocean_step(handle, 1_c_int)
         call check(error, status == OCEAN_STATUS_OK, "cavity: the run must step cleanly")
         if (allocated(error)) return
         status = rdb_ocean_refresh_host(handle)
         !$omp target update from(h%state%surface_flux%heat_cavity)
         hcmax = max(hcmax, maxval(abs(h%state%surface_flux%heat_cavity(ng + 1:ng + NXC, ng + 1:ng + NYC))))
         status = rdb_ocean_get_q_heat_ptr(handle, ptr, nxt, nyt, gen)
         call c_f_pointer(ptr, p2, [nxt, nyt])
         qmax = max(qmax, maxval(abs(p2(ng + 1:ng + NXC, ng + 1:ng + NYC))))
         do j = ng + 1, ng + NYC
            do g = 1, ng
               dh = max(dh, abs(p2(g, j) - p2(g + NXC, j)), abs(p2(ng + NXC + g, j) - p2(ng + g, j)))
            end do
         end do
      end do
      call check(error, hcmax > 1.0_wp .and. qmax > 10.0_wp, &
                 "cavity: the melt heat flux is not live -- the case is vacuous")
      if (allocated(error)) return
      write (msg, "(a,es10.3,a,es10.3,a,es10.3)") "cavity: |dQ_heat| = ", dh, "  max|Q_heat| ", qmax, &
         "  max|heat_cavity| ", hcmax
      call check(error, dh <= 1.0e-12_wp*qmax, trim(msg)//" -- the ghost surface flux is not the "// &
                 "wrapped interior")
      if (allocated(error)) return
      status = rdb_ocean_get_q_salt_ptr(handle, ptr, nxt, nyt, gen)
      call c_f_pointer(ptr, p2, [nxt, nyt])
      smax = maxval(abs(p2(ng + 1:ng + NXC, ng + 1:ng + NYC)))
      ds = 0.0_wp
      do j = ng + 1, ng + NYC
         do g = 1, ng
            ds = max(ds, abs(p2(g, j) - p2(g + NXC, j)), abs(p2(ng + NXC + g, j) - p2(ng + g, j)))
         end do
      end do
      call check(error, smax > 0.0_wp .and. ds <= 1.0e-12_wp*smax, &
                 "cavity: the ghost Q_salt is not the wrapped interior (or the melt salt is zero)")
      status = rdb_ocean_destroy(handle)
   end subroutine test_ghost_cavity

   subroutine write_flux_files()
      !! The atmosphere plus the data-override heat / salt files (wavenumber 2
      !! in x: period 8 cells, so a half-period roll leaves them unchanged).
      call write_atm_file(ATM_FILE, ATM_LIST)
      call write_field_file(HEAT_FILE, -800.0_wp, 0.3_wp)
      call write_field_file(SALT_FILE, 2.0e-4_wp, 0.5_wp)
   end subroutine write_flux_files

   subroutine write_field_file(path, mean, rel)
      !! `field(x, y, time) = mean*(1 + rel*sin(2*pi*2*(i - 1/2)/NX))*(1 + 0.02*j)`,
      !! two records (identical), time axis 0 and 1e6 s.
      character(len=*), intent(in) :: path
      real(wp), intent(in) :: mean, rel
      real(wp), parameter :: PI = 3.14159265358979323846_wp
      real(wp) :: f(NX, NY, 2)
      integer :: ncid, dx_, dy_, dt_, vf, vt, i, j, k, ierr

      do k = 1, 2
         do j = 1, NY
            do i = 1, NX
               f(i, j, k) = mean*(1.0_wp + rel*sin(2.0_wp*PI*2.0_wp*(real(i, wp) - 0.5_wp)/real(NX, wp)))* &
                            (1.0_wp + 0.02_wp*real(j, wp))
            end do
         end do
      end do
      call nc_create_file(path, ncid)
      call nc_def_dim(ncid, "x", NX, dx_)
      call nc_def_dim(ncid, "y", NY, dy_)
      call nc_def_dim(ncid, "time", 2, dt_)
      call nc_def_var_3d(ncid, "field", [dx_, dy_, dt_], vf)
      call rdb_def_var_1d(ncid, "time", dt_, vt)
      call nc_enddef(ncid)
      ierr = nf90_put_var(ncid, vf, f)
      call rdb_put_var_1d(ncid, vt, [0.0_wp, 1.0e6_wp])
      call nc_close(ncid)
   end subroutine write_field_file

   subroutine write_atm_file(path, listname)
      !! One global 10-degree lon-lat file, two time records, the nine CMOR
      !! names: a COLD (255 K) DRY atmosphere with a zonal wavenumber-2
      !! structure (period 180 degrees, so a half-period roll of the model
      !! domain sees the identical atmosphere) and a meridional dependence.
      character(len=*), intent(in) :: path, listname
      integer, parameter :: NLON = 36, NLAT = 18
      real(wp), parameter :: PI = 3.14159265358979323846_wp
      character(len=4), parameter :: VN(9) = [character(len=4) :: "tas", "huss", "uas", &
                                                                                       "vas", "psl", "rlds", "rsds", "prra", "prsn"]
      real(wp) :: lon(NLON), lat(NLAT), lonb(2, NLON), latb(2, NLAT), tax(2), f(NLON, NLAT, 2)
      real(wp) :: sx, sy, dlo, dla
      integer :: ncid, dlon_, dlat_, dbnd, dtim, vlon, vlat, vlonb, vlatb, vtim, vf(9)
      integer :: v, i, j, k, ierr, u

      dlo = 360.0_wp/real(NLON, wp)
      dla = 180.0_wp/real(NLAT, wp)
      do i = 1, NLON
         lonb(1, i) = dlo*real(i - 1, wp)
         lonb(2, i) = dlo*real(i, wp)
         lon(i) = 0.5_wp*(lonb(1, i) + lonb(2, i))
      end do
      do j = 1, NLAT
         latb(1, j) = -90.0_wp + dla*real(j - 1, wp)
         latb(2, j) = -90.0_wp + dla*real(j, wp)
         lat(j) = 0.5_wp*(latb(1, j) + latb(2, j))
      end do
      tax = [0.0_wp, 10.0_wp]
      call nc_create_file(path, ncid)
      call nc_def_dim(ncid, "lon", NLON, dlon_)
      call nc_def_dim(ncid, "lat", NLAT, dlat_)
      call nc_def_dim(ncid, "bnds", 2, dbnd)
      call nc_def_dim(ncid, "time", 2, dtim)
      call rdb_def_var_1d(ncid, "lon", dlon_, vlon)
      call rdb_def_var_1d(ncid, "lat", dlat_, vlat)
      call nc_def_var_2d(ncid, "lon_bnds", [dbnd, dlon_], vlonb)
      call nc_def_var_2d(ncid, "lat_bnds", [dbnd, dlat_], vlatb)
      call rdb_def_var_1d(ncid, "time", dtim, vtim)
      call nc_put_att(ncid, vtim, "units", "days since 2000-06-20 00:00:00")
      do v = 1, 9
         call nc_def_var_3d(ncid, trim(VN(v)), [dlon_, dlat_, dtim], vf(v))
      end do
      call nc_enddef(ncid)
      call rdb_put_var_1d(ncid, vlon, lon)
      call rdb_put_var_1d(ncid, vlat, lat)
      ierr = nf90_put_var(ncid, vlonb, lonb)
      ierr = nf90_put_var(ncid, vlatb, latb)
      call rdb_put_var_1d(ncid, vtim, tax)
      do v = 1, 9
         do k = 1, 2
            do j = 1, NLAT
               do i = 1, NLON
                  sx = sin(2.0_wp*PI*lon(i)/180.0_wp + 0.3_wp)
                  sy = cos(2.0_wp*PI*lat(j)/29.0_wp - 0.5_wp)
                  select case (v)
                  case (1)
                     f(i, j, k) = 255.0_wp + 4.0_wp*sx + 0.15_wp*(lat(j) - 20.0_wp)
                  case (2)
                     f(i, j, k) = 6.0e-4_wp*(1.0_wp + 0.2_wp*sy + 0.1_wp*sx)
                  case (3)
                     f(i, j, k) = 6.0_wp + 4.0_wp*sx + 0.05_wp*(lat(j) - 20.0_wp)
                  case (4)
                     f(i, j, k) = -3.0_wp + 3.0_wp*sy + 2.0_wp*sx
                  case (5)
                     f(i, j, k) = 1.01e5_wp - 300.0_wp*sx + 5.0_wp*lat(j)
                  case (6)
                     f(i, j, k) = 200.0_wp + 15.0_wp*sy + 20.0_wp*sx
                  case (7)
                     f(i, j, k) = 180.0_wp + 40.0_wp*sx*sy + 0.5_wp*lat(j)
                  case (8)
                     f(i, j, k) = 2.0e-5_wp*(1.0_wp + 0.5_wp*sx)
                  case default
                     f(i, j, k) = 1.0e-6_wp*(1.0_wp + 0.5_wp*sy)
                  end select
               end do
            end do
         end do
         ierr = nf90_put_var(ncid, vf(v), f)
      end do
      call nc_close(ncid)
      open (newunit=u, file=listname, status='replace', action='write')
      write (u, '(a)') trim(path)
      close (u)
   end subroutine write_atm_file

end module test_ocean_surface_flux_seam
