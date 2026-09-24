!! Periodic-seam consistency of LOADED bathymetry, through the real setup.
!!
!! A bathymetry that does not come from a formula — the NetCDF file loader
!! (`&ocean_topo_nml topo_config = "file"`) or the API's staged array
!! (`rdb_ocean_stage_bathymetry`) — has only its physical interior defined;
!! the ghosts are filled by CONSTANT EXTRAPOLATION of the nearest edge
!! column.  On a periodic-x domain that is the wrong value at the seam: the
!! west ghosts must hold the EAST edge columns and vice versa.  The engine
!! re-wraps `b` before the first step, but a setup-time consumer that ran
!! before that wrap kept a copy of the extrapolated ghosts — the
!! pressure-gradient force's own bathymetry (`ocean_pressure_force_t%b`,
!! which FV-MOM6 and gprime read at every face to place the bottom
!! interface).  The symptom is a spurious jet on the seam face: 1.9 m/s
!! within three hours on the 1-degree global probe, following the seam when
!! the grid is shifted.
!!
!! The gate is SHIFT INVARIANCE, the defining property of a periodic
!! direction: run A over a bathymetry with a 400 m cliff at the seam, run B
!! over the same bathymetry rolled by `nx/2` (the cliff now in the interior,
!! the seam on a smooth ramp); after `N_STEPS` the rolled-back B must equal
!! A.  The ocean is linearly stratified on sigma layers over the ramp, so it
!! moves (sloping isopycnals — that motion is the same physics in both runs
!! and is what makes the comparison non-vacuous), and the rolled comparison
!! does not care that the motion is not zero.  A stale seam ghost shows up
!! as a difference concentrated on the seam faces.
!!
!! A uniform-density rest test would NOT catch it: a depth-uniform spurious
!! PGF is exactly what the split's barotropic correction removes, so only
!! the baroclinic (stratified) part of the seam error survives.
!!
!! The file cases write their NetCDF under `tmp_local_artifacts/` (ctest's
!! working directory is the build's `tests/`), per the repo artifact rule.
!!
!! Own ctest binary: `rdb_ocean_api` keeps one live handle per process.
module test_ocean_periodic_seam_file
   use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double, c_null_ptr, c_f_pointer
   use rdb_constants, only: wp
   use rdb_ocean_api, only: rdb_ocean_create_from_string, rdb_ocean_create_pending, &
                            rdb_ocean_create_finalize, rdb_ocean_stage_bathymetry, &
                            rdb_ocean_stage_topology, rdb_ocean_step, rdb_ocean_destroy, &
                            rdb_ocean_refresh_host, rdb_ocean_get_grid_info, &
                            rdb_ocean_get_u_face_x_layer_ptr, &
                            rdb_ocean_get_v_face_y_layer_ptr, rdb_ocean_get_h_layer_ptr
   use rdb_ocean_bathymetry_inject, only: BATHY_CONVENTION_DEPTH_POSITIVE_DOWN
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   use rdb_io_netcdf, only: nc_create_file, nc_close, nc_def_dim, nc_def_var_2d, &
                            nc_enddef, nc_put_var_2d
   use testdrive, only: error_type, check, new_unittest, unittest_type
   implicit none
   private

   public :: collect_ocean_periodic_seam_file_tests

   integer, parameter :: NX = 16
      !! Physical cells in x (the periodic direction).  Even, so the roll
      !! `NX/2` is exact.
   integer, parameter :: NY = 8
      !! Physical cells in y.
   integer, parameter :: NZ = 4
      !! Layers.
   integer, parameter :: SHIFT = NX/2
      !! Roll applied to run B's bathymetry.
   integer, parameter :: N_STEPS = 24
      !! Outer steps (2 h at dt = 300 s) — the probe's jet was O(1 m/s)
      !! within 6 steps.
   real(wp), parameter :: B_WEST = 1000.0_wp
      !! Depth of the westernmost physical column (m).
   real(wp), parameter :: B_EAST = 1400.0_wp
      !! Depth of the easternmost physical column (m): a 400 m cliff at the
      !! seam.
   real(wp), parameter :: REL_TOL = 1.0e-10_wp
      !! Shift-invariance tolerance, relative to max |u|.  The two runs do
      !! the same arithmetic on the same numbers, cell for cell.

   type :: run_out_t
      !! Physical-interior copy of one run's end state.
      real(wp), allocatable :: u(:, :, :)
         !! u faces 1..NX (face NX+1 is face 1 on a periodic axis).
      real(wp), allocatable :: v(:, :, :)
         !! v faces, all NY+1 rows.
      real(wp), allocatable :: h(:, :, :)
         !! Layer thickness.
   end type run_out_t

contains

   subroutine collect_ocean_periodic_seam_file_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("file_bathy_seam_shift_invariant_cartesian", test_file_cartesian), &
                  new_unittest("file_bathy_seam_shift_invariant_spherical", test_file_spherical), &
                  new_unittest("staged_bathy_topology_seam_shift_invariant", test_staged_topology) &
                  ]
   end subroutine collect_ocean_periodic_seam_file_tests

   pure function seam_ramp(roll) result(b)
      !! Zonal ramp from `B_WEST` to `B_EAST` (a cliff at the seam) plus a
      !! gentle meridional tilt, circularly rolled by `roll` columns:
      !! `b_rolled(i) = b(i - roll)`.  Positive-down depth, interior only.
      integer, intent(in) :: roll
      real(wp) :: b(NX, NY)
      integer :: i, j, i_src
      do j = 1, NY
         do i = 1, NX
            i_src = modulo(i - 1 - roll, NX) + 1
            b(i, j) = B_WEST + (B_EAST - B_WEST)*real(i_src - 1, wp)/real(NX - 1, wp) &
                      + 20.0_wp*real(j - 1, wp)
         end do
      end do
   end function seam_ramp

   subroutine write_bathy_nc(filename, b)
      !! Minimal bathymetry file for `load_bathymetry_into_array`:
      !! dimensions `x`, `y`, variable `depth(x, y)` positive-down.
      character(len=*), intent(in) :: filename
      real(wp), intent(in) :: b(:, :)
      integer :: ncid, dim_x, dim_y, vid
      call nc_create_file(filename, ncid)
      call nc_def_dim(ncid, "x", size(b, 1), dim_x)
      call nc_def_dim(ncid, "y", size(b, 2), dim_y)
      call nc_def_var_2d(ncid, "depth", [dim_x, dim_y], vid)
      call nc_enddef(ncid)
      call nc_put_var_2d(ncid, vid, b)
      call nc_close(ncid)
   end subroutine write_bathy_nc

   function common_nml() result(nml)
      !! Everything but the grid, the edges and the bathymetry source:
      !! linearly stratified T on sigma layers, FV-MOM6 PGF (reads the PGF's
      !! own bathymetry copy), sigma coordinate.
      character(len=:), allocatable :: nml
      character(len=1), parameter :: nl = new_line("a")
      nml = "&sim_nml sim_type = 'ocean' /"//nl// &
            "&time_nml t_end = 1.0e9, dt_fixed = 300.0, cfl_interval = 1 /"//nl// &
            "&nonhydrostatic_nml nz_layers = 4 /"//nl// &
            "&vcoord_nml vcoord_type = 'sigma' /"//nl// &
            "&tracer_nml initial_temperature = 10.0, initial_salinity = 35.0, "// &
            "T_init_surface = 20.0, T_init_bottom = 4.0 /"//nl// &
            "&ocean_pgf_nml form = 'fv_mom6' /"//nl// &
            "&ocean_coriolis_nml form = 'sadourny_energy' /"//nl// &
            "&ocean_hvisc_nml nu_h = 100.0 /"//nl// &
            "&ocean_bdrag_nml form = 'quadratic', cd = 3.0e-3 /"//nl// &
            "&ocean_bt_nml auto_n_inner = .false., n_inner = 20 /"//nl// &
            "&ocean_diag_nml enabled = .false. /"//nl
   end function common_nml

   function case_nml(grid_kind, bathy_file, walls_x) result(nml)
      !! `grid_kind` "cartesian" (20 km f-plane channel) or "spherical"
      !! (a 360-degree zonal band, 22.5 x 4 degree cells).  `bathy_file`
      !! blank ⇒ no `topo_config = "file"` (the staged-bathymetry case).
      !! `walls_x` tags the east/west edges WALL in the namelist, so
      !! periodicity can only come from the staged topology.
      character(len=*), intent(in) :: grid_kind, bathy_file
      logical, intent(in) :: walls_x
      character(len=:), allocatable :: nml
      character(len=1), parameter :: nl = new_line("a")
      character(len=:), allocatable :: ew

      ew = "periodic"
      if (walls_x) ew = "wall"
      nml = common_nml()//"&ocean_bc_nml west = '"//ew//"', east = '"//ew//"', "// &
            "south = 'wall', north = 'wall' /"//nl
      if (grid_kind == "spherical") then
         nml = nml//"&grid_nml nx = 16, ny = 8, nghost = 3, dx = 22.5, dy = 4.0 /"//nl// &
               "&ocean_grid_nml grid_config = 'spherical', lon_west = -300.0, "// &
               "lat_south = 20.0, rad_earth = 6.371e6, coriolis_scheme = 'planetary' /"//nl
      else
         nml = nml//"&grid_nml nx = 16, ny = 8, nghost = 3, dx = 20000.0, dy = 20000.0 /"//nl// &
               "&physics_nml coriolis_f = 1.0e-4 /"//nl
      end if
      if (len_trim(bathy_file) > 0) then
         nml = nml//"&ocean_topo_nml topo_config = 'file', max_depth = 2000.0 /"//nl// &
               "&output_nml bathymetry_file = '"//bathy_file//"', output_to_file = .false. /"//nl
      else
         nml = nml//"&ocean_topo_nml max_depth = 2000.0 /"//nl// &
               "&output_nml output_to_file = .false. /"//nl
      end if
   end function case_nml

   subroutine step_and_capture(error, handle, what, out)
      !! Step `N_STEPS` and copy the physical interior of u, v, h out.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr), intent(in) :: handle
      character(len=*), intent(in) :: what
      type(run_out_t), intent(out) :: out
      type(c_ptr) :: ptr
      integer(c_int) :: status, nx_p, ny_p, nz_p, ng, nxt, nyt, nzt, gen
         !! `nxt`/`nyt`/`nzt`: the pointer extents (Fortran is case-blind —
         !! `nx` here would shadow the module's `NX`).
      real(wp), pointer :: p3(:, :, :)

      status = rdb_ocean_get_grid_info(handle, nx_p, ny_p, nz_p, ng)
      status = rdb_ocean_step(handle, int(N_STEPS, c_int))
      call check(error, status == OCEAN_STATUS_OK, what//": the run must step cleanly")
      if (allocated(error)) return
      status = rdb_ocean_refresh_host(handle)
      status = rdb_ocean_get_u_face_x_layer_ptr(handle, ptr, nxt, nyt, nzt, gen)
      call c_f_pointer(ptr, p3, [nxt, nyt, nzt])
      out%u = p3(ng + 1:ng + NX, ng + 1:ng + NY, 1:NZ)
      status = rdb_ocean_get_v_face_y_layer_ptr(handle, ptr, nxt, nyt, nzt, gen)
      call c_f_pointer(ptr, p3, [nxt, nyt, nzt])
      out%v = p3(ng + 1:ng + NX, ng + 1:ng + NY + 1, 1:NZ)
      status = rdb_ocean_get_h_layer_ptr(handle, ptr, nxt, nyt, nzt, gen)
      call c_f_pointer(ptr, p3, [nxt, nyt, nzt])
      out%h = p3(ng + 1:ng + NX, ng + 1:ng + NY, 1:NZ)
   end subroutine step_and_capture

   subroutine check_shift_invariant(error, a, b, what)
      !! Run B was rolled by `SHIFT`: B(i + SHIFT) must equal A(i).  u
      !! faces roll like centres (face i is the west face of cell i).
      type(error_type), allocatable, intent(out) :: error
      type(run_out_t), intent(in) :: a, b
      character(len=*), intent(in) :: what
      real(wp) :: umax, du, dv, dh, d
      integer :: i, ib, i_worst
      character(len=200) :: msg

      umax = max(maxval(abs(a%u)), maxval(abs(a%v)))
      call check(error, umax > 1.0e-6_wp, what//": the stratified ocean did not move — "// &
                 "the comparison is vacuous")
      if (allocated(error)) return
      du = 0.0_wp
      dv = 0.0_wp
      dh = 0.0_wp
      i_worst = 0
      do i = 1, NX
         ib = modulo(i - 1 + SHIFT, NX) + 1
         d = maxval(abs(a%u(i, :, :) - b%u(ib, :, :)))
         if (d > du) then
            du = d
            i_worst = i
         end if
         dv = max(dv, maxval(abs(a%v(i, :, :) - b%v(ib, :, :))))
         dh = max(dh, maxval(abs(a%h(i, :, :) - b%h(ib, :, :))))
      end do
      write (msg, "(a,es10.3,a,es10.3,a,es10.3,a,es10.3,a,i0)") ": max|u| = ", umax, &
         "  |du| = ", du, "  |dv| = ", dv, "  |dh| = ", dh, "  worst u face i = ", i_worst
      call check(error, max(du, dv) <= REL_TOL*umax .and. dh <= REL_TOL*B_EAST, &
                 what//trim(msg)//" — the run is not shift-invariant; the periodic "// &
                 "seam saw a stale ghost")
   end subroutine check_shift_invariant

   subroutine run_file(error, grid_kind, fname, roll, out)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: grid_kind, fname
      integer, intent(in) :: roll
      type(run_out_t), intent(out) :: out
      type(c_ptr) :: handle
      integer(c_int) :: status
      character(len=:), allocatable :: nml

      call write_bathy_nc(fname, seam_ramp(roll))
      nml = case_nml(grid_kind, fname, walls_x=.false.)
      handle = c_null_ptr
      status = rdb_ocean_create_from_string(nml, len(nml, kind=c_int), handle)
      call check(error, status == OCEAN_STATUS_OK, grid_kind//": the namelist must build")
      if (allocated(error)) return
      call step_and_capture(error, handle, grid_kind//" file bathymetry", out)
      status = rdb_ocean_destroy(handle)
   end subroutine run_file

   subroutine run_file_pair(error, grid_kind)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: grid_kind
      type(run_out_t) :: a, b
      call run_file(error, grid_kind, "tmp_local_artifacts/seam_bathy_"//grid_kind//"_a.nc", 0, a)
      if (allocated(error)) return
      call run_file(error, grid_kind, "tmp_local_artifacts/seam_bathy_"//grid_kind//"_b.nc", &
                    SHIFT, b)
      if (allocated(error)) return
      call check_shift_invariant(error, a, b, grid_kind//" file bathymetry, periodic-x")
   end subroutine run_file_pair

   subroutine test_file_cartesian(error)
      type(error_type), allocatable, intent(out) :: error
      call run_file_pair(error, "cartesian")
   end subroutine test_file_cartesian

   subroutine test_file_spherical(error)
      !! The probe reproduced the jet on `spherical` (no fold) as well.
      type(error_type), allocatable, intent(out) :: error
      call run_file_pair(error, "spherical")
   end subroutine test_file_spherical

   subroutine run_staged(error, roll, out)
      !! The API path: staged bathymetry (ghosts constant-extrapolated by
      !! `bathymetry_fill_ghosts_array`) + staged periodic-x topology over a
      !! namelist that tags the east/west edges WALL — so the seed has to
      !! learn the periodicity from the stage, not from `&ocean_bc_nml`.
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: roll
      type(run_out_t), intent(out) :: out
      type(c_ptr) :: handle
      integer(c_int) :: status
      character(len=:), allocatable :: nml
      real(c_double) :: b(NX, NY)

      b = real(seam_ramp(roll), c_double)
      nml = case_nml("cartesian", "", walls_x=.true.)
      handle = c_null_ptr
      status = rdb_ocean_create_pending(nml, len(nml, kind=c_int), handle)
      call check(error, status == OCEAN_STATUS_OK, "staged: create_pending")
      if (allocated(error)) return
      body: block
         status = rdb_ocean_stage_bathymetry(handle, b, int(NX, c_int), int(NY, c_int), &
                                             int(BATHY_CONVENTION_DEPTH_POSITIVE_DOWN, c_int))
         call check(error, status == OCEAN_STATUS_OK, "staged: stage_bathymetry")
         if (allocated(error)) exit body
         status = rdb_ocean_stage_topology(handle, 1_c_int, 0_c_int)
         call check(error, status == OCEAN_STATUS_OK, "staged: stage_topology")
         if (allocated(error)) exit body
         status = rdb_ocean_create_finalize(handle)
         call check(error, status == OCEAN_STATUS_OK, "staged: create_finalize")
         if (allocated(error)) exit body
         call step_and_capture(error, handle, "staged bathymetry", out)
      end block body
      status = rdb_ocean_destroy(handle)
   end subroutine run_staged

   subroutine test_staged_topology(error)
      type(error_type), allocatable, intent(out) :: error
      type(run_out_t) :: a, b
      call run_staged(error, 0, a)
      if (allocated(error)) return
      call run_staged(error, SHIFT, b)
      if (allocated(error)) return
      call check_shift_invariant(error, a, b, "staged bathymetry + staged periodic-x topology")
   end subroutine test_staged_topology

end module test_ocean_periodic_seam_file
