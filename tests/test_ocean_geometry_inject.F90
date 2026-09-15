!! Unit tests for P2.5 pre-create geometry injection: the two-phase
!! create (`rdb_ocean_create_pending` -> `rdb_ocean_stage_*` ->
!! `rdb_ocean_create_finalize`), bathymetry sign validation + ghost
!! fill (`rdb_ocean_bathymetry_inject`), grid-owns-periodicity topology
!! injection (`ocean_bc_state_set_topology`), and the folded
!! `required_halo` query. See `docs/ocean_python_api_plan.md` S5b/S6 and
!! `tmp_local_artifacts/python_ffi_scope/06_python_surface_design.md` D6.2.
!!
!! Mixed style, deliberately: the create()-path tests drive the real
!! `bind(c)` ABI (the same entry points a Python/C caller uses), while the
!! topology + required_halo tests call the underlying Fortran directly
!! (`ocean_bc_state_set_topology`, `required_halo`) — both are pure
!! functions of their inputs with no handle/process-global state, so a
!! direct unit test is the more precise tool (same precedent as
!! `test_ocean_boundary.F90`).
!!
!! Own ctest binary: the single-live-handle module state
!! (`rdb_ocean_api`'s `g_handle_live`) must start from a fresh process,
!! same reasoning as `test_rdb_ocean_api` / `test_ocean_api_p2`.
module test_ocean_geometry_inject
   use, intrinsic :: iso_c_binding, only: c_associated, c_char, c_double, c_f_pointer, &
                                                                             c_int, c_null_ptr, c_ptr
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_state, only: set_bathymetry_seamount
   use rdb_ocean_boundary_types, only: ocean_bc_state_t, ocean_bc_state_init, &
                                       ocean_bc_state_destroy, ocean_bc_state_set_topology, &
                                       OBC_WALL, OBC_PERIODIC
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_SETUP
   use rdb_ocean_halo_width, only: required_halo
   implicit none
   private

   public :: collect_ocean_geometry_inject_tests

   interface
      function rdb_ocean_create_from_string(nml_text, nml_len, handle_out) &
         result(status) bind(c, name="rdb_ocean_create_from_string")
         import :: c_char, c_int, c_ptr
         implicit none
         integer(c_int), intent(in), value :: nml_len
         character(kind=c_char), intent(in) :: nml_text(nml_len)
         type(c_ptr), intent(out) :: handle_out
         integer(c_int) :: status
      end function rdb_ocean_create_from_string

      function rdb_ocean_create_pending(nml_text, nml_len, handle_out) &
         result(status) bind(c, name="rdb_ocean_create_pending")
         import :: c_char, c_int, c_ptr
         implicit none
         integer(c_int), intent(in), value :: nml_len
         character(kind=c_char), intent(in) :: nml_text(nml_len)
         type(c_ptr), intent(out) :: handle_out
         integer(c_int) :: status
      end function rdb_ocean_create_pending

      function rdb_ocean_create_finalize(c_handle) result(status) &
         bind(c, name="rdb_ocean_create_finalize")
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), intent(inout) :: c_handle
         integer(c_int) :: status
      end function rdb_ocean_create_finalize

      function rdb_ocean_stage_bathymetry(c_handle, b_data, nx_p, ny_p, convention) &
         result(status) bind(c, name="rdb_ocean_stage_bathymetry")
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), intent(in), value :: c_handle
         integer(c_int), intent(in), value :: nx_p, ny_p
         real(c_double), intent(in) :: b_data(nx_p, ny_p)
         integer(c_int), intent(in), value :: convention
         integer(c_int) :: status
      end function rdb_ocean_stage_bathymetry

      function rdb_ocean_destroy(c_handle) result(status) &
         bind(c, name="rdb_ocean_destroy")
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), intent(inout) :: c_handle
         integer(c_int) :: status
      end function rdb_ocean_destroy

      function rdb_ocean_get_b_ptr(c_handle, ptr, nx, ny, gen) result(status) &
         bind(c, name="rdb_ocean_get_b_ptr")
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), intent(in), value :: c_handle
         type(c_ptr), intent(out) :: ptr
         integer(c_int), intent(out) :: nx, ny, gen
         integer(c_int) :: status
      end function rdb_ocean_get_b_ptr

      function rdb_ocean_get_grid_info(c_handle, nx, ny, nz, nghost) result(status) &
         bind(c, name="rdb_ocean_get_grid_info")
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), intent(in), value :: c_handle
         integer(c_int), intent(out) :: nx, ny, nz, nghost
         integer(c_int) :: status
      end function rdb_ocean_get_grid_info

      function rdb_ocean_required_halo(pv_adv_scheme, pv_adv_scheme_len, &
                                       tracer_recon, tracer_recon_len, &
                                       periodic, tripolar_fold, decomposed, &
                                       kappa_shear_at_vertex) result(ng) &
         bind(c, name="rdb_ocean_required_halo")
         import :: c_char, c_int
         implicit none
         integer(c_int), intent(in), value :: pv_adv_scheme_len
         character(kind=c_char), intent(in) :: pv_adv_scheme(pv_adv_scheme_len)
         integer(c_int), intent(in), value :: tracer_recon_len
         character(kind=c_char), intent(in) :: tracer_recon(tracer_recon_len)
         integer(c_int), intent(in), value :: periodic, tripolar_fold, decomposed
         integer(c_int), intent(in), value :: kappa_shear_at_vertex
         integer(c_int) :: ng
      end function rdb_ocean_required_halo
   end interface

   integer(c_int), parameter :: OK = 0
   integer(c_int), parameter :: ERR_ALREADY_EXISTS = 11
   integer(c_int), parameter :: ERR_NOT_INITIALISED = 12
   integer(c_int), parameter :: ERR_BAD_SHAPE = 14
   integer(c_int), parameter :: ERR_BATHYMETRY_SIGN = 15
   integer(c_int), parameter :: ERR_NOT_PENDING = 16

   integer(c_int), parameter :: BATHY_DEPTH_POSITIVE_DOWN = 1
   integer(c_int), parameter :: BATHY_HEIGHT_POSITIVE_UP = 2

   ! Shared tiny grid for every test below.
   integer, parameter :: ONX = 8, ONY = 6, ONZ = 2, ONGHOST = 3
   real(wp), parameter :: ODX = 1000.0_wp, ODY = 1000.0_wp
   real(wp), parameter :: SEAMOUNT_MAX_DEPTH = 500.0_wp
   real(wp), parameter :: SEAMOUNT_PEAK_DEPTH = 50.0_wp
   real(wp), parameter :: SEAMOUNT_HALF_WIDTH = 3000.0_wp

contains

   subroutine collect_ocean_geometry_inject_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("bathymetry_inject_matches_formula", &
                               test_bathymetry_inject_matches_formula), &
                  new_unittest("bathymetry_sign_reject_all_negative", &
                               test_bathymetry_sign_reject), &
                  new_unittest("bathymetry_ghost_rows_filled", &
                               test_ghost_rows_filled), &
                  new_unittest("bathymetry_stage_requires_pending_handle", &
                               test_stage_requires_pending), &
                  new_unittest("bathymetry_bad_shape_deferred_to_finalize", &
                               test_bad_shape_deferred), &
                  new_unittest("topology_injection_sets_periodic_and_edges", &
                               test_topology_injection), &
                  new_unittest("topology_injection_rejects_thin_halo", &
                               test_topology_injection_thin_halo), &
                  new_unittest("required_halo_folds_scattered_minimums", &
                               test_required_halo_minimums), &
                  new_unittest("required_halo_abi_matches_fortran", &
                               test_required_halo_abi_parity) &
                  ]
   end subroutine collect_ocean_geometry_inject_tests

   function seamount_nml() result(txt)
      !! Namelist form of the SAME seamount this file builds directly via
      !! `set_bathymetry_seamount` below — the "two routes, same field"
      !! reference trajectory.
      character(len=:), allocatable :: txt
      txt = '&sim_nml sim_type = "ocean" /'//new_line("a")// &
            "&grid_nml nx = "//itoa(ONX)//", ny = "//itoa(ONY)// &
            ", dx = "//rtoa(ODX)//", dy = "//rtoa(ODY)// &
            ", nghost = "//itoa(ONGHOST)//" /"//new_line("a")// &
            "&nonhydrostatic_nml nz_layers = "//itoa(ONZ)//" /"//new_line("a")// &
            "&time_nml t_end = 86400.0, dt_fixed = 300.0 /"//new_line("a")// &
            '&ocean_topo_nml topo_config = "seamount", max_depth = '//rtoa(SEAMOUNT_MAX_DEPTH)// &
            ", edge_depth = "//rtoa(SEAMOUNT_PEAK_DEPTH)// &
            ", slope_scale = "//rtoa(SEAMOUNT_HALF_WIDTH)//" /"//new_line("a")// &
            "&ocean_bt_nml auto_n_inner = .true. /"//new_line("a")// &
            "&tracer_nml initial_salinity = 35.0, initial_temperature = 12.0 /"//new_line("a")// &
            "&ocean_diag_nml enabled = .false. /"//new_line("a")// &
            "&output_nml output_to_file = .false. /"//new_line("a")
   end function seamount_nml

   function inject_nml(topo_max_depth) result(txt)
      !! Base namelist for the injection-path handles: `topo_config="flat"`
      !! at a DIFFERENT depth than anything injected, so a silent fallback
      !! to the namelist path (injection not actually taking effect) would
      !! be caught by comparison against the seamount reference.
      real(wp), intent(in) :: topo_max_depth
      character(len=:), allocatable :: txt
      txt = '&sim_nml sim_type = "ocean" /'//new_line("a")// &
            "&grid_nml nx = "//itoa(ONX)//", ny = "//itoa(ONY)// &
            ", dx = "//rtoa(ODX)//", dy = "//rtoa(ODY)// &
            ", nghost = "//itoa(ONGHOST)//" /"//new_line("a")// &
            "&nonhydrostatic_nml nz_layers = "//itoa(ONZ)//" /"//new_line("a")// &
            "&time_nml t_end = 86400.0, dt_fixed = 300.0 /"//new_line("a")// &
            '&ocean_topo_nml topo_config = "flat", max_depth = '//rtoa(topo_max_depth)//" /"// &
            new_line("a")// &
            "&ocean_bt_nml auto_n_inner = .true. /"//new_line("a")// &
            "&tracer_nml initial_salinity = 35.0, initial_temperature = 12.0 /"//new_line("a")// &
            "&ocean_diag_nml enabled = .false. /"//new_line("a")// &
            "&output_nml output_to_file = .false. /"//new_line("a")
   end function inject_nml

   function itoa(n) result(s)
      integer, intent(in) :: n
      character(len=:), allocatable :: s
      character(len=32) :: buf
      write (buf, "(I0)") n
      s = trim(buf)
   end function itoa

   function rtoa(x) result(s)
      real(wp), intent(in) :: x
      character(len=:), allocatable :: s
      character(len=32) :: buf
      write (buf, "(F16.4)") x
      s = trim(adjustl(buf))
   end function rtoa

   subroutine seamount_interior(b_interior)
      !! Build the reference seamount interior array directly via the
      !! SAME production formula the namelist path dispatches to
      !! (`ocean_state_seed_from_cfg`'s "seamount" case) — the "two
      !! routes" both tests in this file compare against.
      real(wp), intent(out) :: b_interior(ONX, ONY)
      type(hgrid_t) :: grid
      real(wp), allocatable :: b_full(:, :)

      call grid%init(ONX, ONY, ONGHOST, ODX, ODY)
      allocate (b_full(grid%nx_total, grid%ny_total))
      call set_bathymetry_seamount(b_full, grid, SEAMOUNT_MAX_DEPTH, SEAMOUNT_PEAK_DEPTH, &
                                   SEAMOUNT_HALF_WIDTH)
      b_interior = b_full(ONGHOST + 1:ONGHOST + ONX, ONGHOST + 1:ONGHOST + ONY)
   end subroutine seamount_interior

   subroutine test_bathymetry_inject_matches_formula(error)
      !! Inject the seamount's interior array through the ABI
      !! (create_pending -> stage_bathymetry -> create_finalize) and
      !! assert the resulting `%barotropic%b` interior matches, cell for
      !! cell, what `topo_config="seamount"` produces through the ordinary
      !! namelist path — same field, two routes (the headline P2.5 test).
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle, ptr
      character(len=:), allocatable :: nml
      integer(c_int) :: status, nx, ny, gen, nghost_out, nz_out
      real(wp) :: b_interior(ONX, ONY)
      real(c_double) :: b_data(ONX, ONY)
      real(wp), pointer :: bp(:, :)
      real(wp), allocatable :: b_formula_snapshot(:, :)
      integer :: ng_formula

      call seamount_interior(b_interior)
      b_data = real(b_interior, c_double)

      ! ---- Route 1: the ordinary namelist formula path. Only ONE ocean
      ! handle may be live at a time (g_handle_live) -- snapshot the
      ! interior into an OWNED array and destroy this handle before
      ! creating route 2's, exactly as test_ocean_api_p2's
      ! test_setter_no_clobber does. ----
      nml = seamount_nml()
      handle = c_null_ptr
      status = rdb_ocean_create_from_string(nml, len(nml, kind=c_int), handle)
      route1: block
         call check(error, status == OK, "formula-path create succeeds")
         if (allocated(error)) exit route1
         status = rdb_ocean_get_b_ptr(handle, ptr, nx, ny, gen)
         call check(error, status == OK, "formula-path get_b_ptr succeeds")
         if (allocated(error)) exit route1
         call c_f_pointer(ptr, bp, [int(nx), int(ny)])
         status = rdb_ocean_get_grid_info(handle, nx, ny, nz_out, nghost_out)
         call check(error, status == OK, "formula-path get_grid_info succeeds")
         if (allocated(error)) exit route1
         ng_formula = int(nghost_out)
         b_formula_snapshot = bp(ng_formula + 1:ng_formula + ONX, ng_formula + 1:ng_formula + ONY)
      end block route1
      status = rdb_ocean_destroy(handle)
      if (allocated(error)) then
         call check(error, status == OK, "cleanup destroy (formula)")
         return
      end if

      ! ---- Route 2: pre-create geometry injection. ----
      nml = inject_nml(1.0_wp)   ! deliberately different from the seamount
      handle = c_null_ptr
      status = rdb_ocean_create_pending(nml, len(nml, kind=c_int), handle)
      route2: block
         call check(error, status == OK, "injection-path create_pending succeeds")
         if (allocated(error)) exit route2
         status = rdb_ocean_stage_bathymetry(handle, b_data, int(ONX, c_int), &
                                             int(ONY, c_int), BATHY_DEPTH_POSITIVE_DOWN)
         call check(error, status == OK, "stage_bathymetry succeeds")
         if (allocated(error)) exit route2
         status = rdb_ocean_create_finalize(handle)
         call check(error, status == OK, "create_finalize succeeds")
         if (allocated(error)) exit route2
         status = rdb_ocean_get_b_ptr(handle, ptr, nx, ny, gen)
         call check(error, status == OK, "injection-path get_b_ptr succeeds")
         if (allocated(error)) exit route2
         call c_f_pointer(ptr, bp, [int(nx), int(ny)])

         status = rdb_ocean_get_grid_info(handle, nx, ny, nz_out, nghost_out)
         call check(error, status == OK, "get_grid_info succeeds")
         if (allocated(error)) exit route2

         ! Compare the physical INTERIOR only -- ghosts are filled by
         ! DIFFERENT mechanisms on the two routes (the seamount formula
         ! extends its own smooth function into ghost cells; injection
         ! uses constant extrapolation), so they are not expected to
         ! agree there (see bathymetry_ghost_rows_filled for that check).
         associate (ng => int(nghost_out))
            call check(error, &
                       maxval(abs(b_formula_snapshot - &
                                  bp(ng + 1:ng + ONX, ng + 1:ng + ONY))) < 1.0e-8_wp, &
                       "injected interior matches the namelist seamount formula, cell for cell")
            if (allocated(error)) exit route2
            call check(error, &
                       maxval(abs(bp(ng + 1:ng + ONX, ng + 1:ng + ONY) - b_interior)) < 1.0e-8_wp, &
                       "injected interior matches the direct set_bathymetry_seamount reference")
         end associate
      end block route2

      status = rdb_ocean_destroy(handle)
      if (.not. allocated(error)) call check(error, status == OK, "cleanup destroy (inject)")
   end subroutine test_bathymetry_inject_matches_formula

   subroutine test_bathymetry_sign_reject(error)
      !! D6.2's headline safety property: a GEBCO-style array (negative
      !! height) staged under the WRONG convention (claiming it is
      !! already positive-down depth) normalises to zero wet cells and
      !! MUST be rejected with RDB_OCEAN_ERR_BATHYMETRY_SIGN at
      !! create_finalize -- not silently accepted as an all-land basin.
      !! The complementary CORRECT-convention call on the identical array
      !! must succeed, to prove the rejection is about the sign
      !! mismatch, not the array itself.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle
      character(len=:), allocatable :: nml
      integer(c_int) :: status
      integer :: i
      real(c_double) :: b_gebco(ONX, ONY)

      ! GEBCO-style: negative heights, e.g. -50 m to -500 m.
      b_gebco = -reshape([(50.0_c_double + 10.0_c_double*real(mod(i, 46), c_double), &
                           i=1, ONX*ONY)], [ONX, ONY])

      ! ---- Wrong convention: claim positive-down on a negative array. ----
      nml = inject_nml(1.0_wp)
      handle = c_null_ptr
      status = rdb_ocean_create_pending(nml, len(nml, kind=c_int), handle)
      wrong: block
         call check(error, status == OK, "create_pending succeeds")
         if (allocated(error)) exit wrong
         status = rdb_ocean_stage_bathymetry(handle, b_gebco, int(ONX, c_int), int(ONY, c_int), &
                                             BATHY_DEPTH_POSITIVE_DOWN)
         call check(error, status == OK, "stage_bathymetry (wrong convention) itself still succeeds "// &
                    "-- validation happens at finalize, not staging")
         if (allocated(error)) exit wrong
         status = rdb_ocean_create_finalize(handle)
         call check(error, status == ERR_BATHYMETRY_SIGN, &
                    "wrong convention on an all-negative array is rejected with "// &
                    "RDB_OCEAN_ERR_BATHYMETRY_SIGN, not silently accepted")
         if (allocated(error)) exit wrong
         call check(error,.not. c_associated(handle), &
                    "the whole handle was destroyed on this failure (F9)")
      end block wrong
      if (c_associated(handle)) status = rdb_ocean_destroy(handle)
      if (allocated(error)) return

      ! ---- Correct convention: the SAME array, flipped, must succeed. ----
      handle = c_null_ptr
      status = rdb_ocean_create_pending(nml, len(nml, kind=c_int), handle)
      right: block
         call check(error, status == OK, "create_pending succeeds (correct-convention run)")
         if (allocated(error)) exit right
         status = rdb_ocean_stage_bathymetry(handle, b_gebco, int(ONX, c_int), int(ONY, c_int), &
                                             BATHY_HEIGHT_POSITIVE_UP)
         call check(error, status == OK, "stage_bathymetry (correct convention) succeeds")
         if (allocated(error)) exit right
         status = rdb_ocean_create_finalize(handle)
         call check(error, status == OK, &
                    "correct convention (height_positive_up) on the SAME array succeeds -- "// &
                    "the rejection above was about the sign mismatch, not the data")
      end block right
      status = rdb_ocean_destroy(handle)
      if (.not. allocated(error)) call check(error, status == OK, "cleanup destroy")
   end subroutine test_bathymetry_sign_reject

   subroutine test_ghost_rows_filled(error)
      !! CLAUDE.md's formula-bathymetry ghost-fill gotcha, for the
      !! injection path: an unfilled ghost row leaves b=0 there, the EOS
      !! falls back to rho=rho_0, and the spurious density jump at the
      !! wall-adjacent face e-folds the basin in ~12 h. A UNIFORM
      !! interior makes the expected ghost value unambiguous (constant
      !! extrapolation of a constant is that same constant), so this
      !! checks both "non-zero" AND "correct", not just "non-zero".
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle, ptr
      character(len=:), allocatable :: nml
      integer(c_int) :: status, nx, ny, gen, nghost_out, nz_out
      real(c_double), parameter :: DEPTH = 300.0_c_double
      real(c_double) :: b_data(ONX, ONY)
      real(wp), pointer :: b(:, :)

      b_data = DEPTH
      nml = inject_nml(1.0_wp)
      handle = c_null_ptr
      status = rdb_ocean_create_pending(nml, len(nml, kind=c_int), handle)
      checks: block
         call check(error, status == OK, "create_pending succeeds")
         if (allocated(error)) exit checks
         status = rdb_ocean_stage_bathymetry(handle, b_data, int(ONX, c_int), int(ONY, c_int), &
                                             BATHY_DEPTH_POSITIVE_DOWN)
         call check(error, status == OK, "stage_bathymetry succeeds")
         if (allocated(error)) exit checks
         status = rdb_ocean_create_finalize(handle)
         call check(error, status == OK, "create_finalize succeeds")
         if (allocated(error)) exit checks
         status = rdb_ocean_get_grid_info(handle, nx, ny, nz_out, nghost_out)
         call check(error, status == OK, "get_grid_info succeeds")
         if (allocated(error)) exit checks
         status = rdb_ocean_get_b_ptr(handle, ptr, nx, ny, gen)
         call check(error, status == OK, "get_b_ptr succeeds")
         if (allocated(error)) exit checks
         call c_f_pointer(ptr, b, [int(nx), int(ny)])

         associate (ng => int(nghost_out))
            ! West ghost column, south ghost row, and the SW ghost corner.
            call check(error, abs(b(1, ng + 1) - real(DEPTH, wp)) < 1.0e-9_wp, &
                       "west ghost column filled (non-zero, == interior depth)")
            if (allocated(error)) exit checks
            call check(error, abs(b(ng + 1, 1) - real(DEPTH, wp)) < 1.0e-9_wp, &
                       "south ghost row filled (non-zero, == interior depth)")
            if (allocated(error)) exit checks
            call check(error, abs(b(1, 1) - real(DEPTH, wp)) < 1.0e-9_wp, &
                       "SW ghost corner filled (non-zero, == interior depth)")
            if (allocated(error)) exit checks
            ! East ghost column, north ghost row.
            call check(error, abs(b(ng + ONX + ng, ng + 1) - real(DEPTH, wp)) < 1.0e-9_wp, &
                       "east ghost column filled (non-zero, == interior depth)")
            if (allocated(error)) exit checks
            call check(error, abs(b(ng + 1, ng + ONY + ng) - real(DEPTH, wp)) < 1.0e-9_wp, &
                       "north ghost row filled (non-zero, == interior depth)")
         end associate
      end block checks
      status = rdb_ocean_destroy(handle)
      if (.not. allocated(error)) call check(error, status == OK, "cleanup destroy")
   end subroutine test_ghost_rows_filled

   subroutine test_stage_requires_pending(error)
      !! A `rdb_ocean_stage_*` call on a FULLY-CREATED (not pending)
      !! handle is rejected with RDB_OCEAN_ERR_NOT_PENDING, not
      !! silently ignored or applied mid-run.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle
      character(len=:), allocatable :: nml
      integer(c_int) :: status
      real(c_double) :: b_data(ONX, ONY)

      b_data = 200.0_c_double
      nml = inject_nml(200.0_wp)
      handle = c_null_ptr
      status = rdb_ocean_create_from_string(nml, len(nml, kind=c_int), handle)
      checks: block
         call check(error, status == OK, "ordinary create succeeds")
         if (allocated(error)) exit checks
         status = rdb_ocean_stage_bathymetry(handle, b_data, int(ONX, c_int), int(ONY, c_int), &
                                             BATHY_DEPTH_POSITIVE_DOWN)
         call check(error, status == ERR_NOT_PENDING, &
                    "stage_bathymetry on a fully-created (non-pending) handle is rejected")
      end block checks
      status = rdb_ocean_destroy(handle)
      if (.not. allocated(error)) call check(error, status == OK, "cleanup destroy")
   end subroutine test_stage_requires_pending

   subroutine test_bad_shape_deferred(error)
      !! Shape cannot be validated at stage_bathymetry time (the grid
      !! does not exist yet) -- it surfaces as RDB_OCEAN_ERR_BAD_SHAPE
      !! at create_finalize instead, once the grid is built.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle
      character(len=:), allocatable :: nml
      integer(c_int) :: status
      real(c_double) :: b_data(ONX + 1, ONY)   ! deliberately wrong nx_p

      b_data = 200.0_c_double
      nml = inject_nml(200.0_wp)
      handle = c_null_ptr
      status = rdb_ocean_create_pending(nml, len(nml, kind=c_int), handle)
      checks: block
         call check(error, status == OK, "create_pending succeeds")
         if (allocated(error)) exit checks
         status = rdb_ocean_stage_bathymetry(handle, b_data, int(ONX + 1, c_int), int(ONY, c_int), &
                                             BATHY_DEPTH_POSITIVE_DOWN)
         call check(error, status == OK, "stage_bathymetry with a wrong-but-self-consistent "// &
                    "shape still succeeds -- shape is checked against the grid at finalize")
         if (allocated(error)) exit checks
         status = rdb_ocean_create_finalize(handle)
         call check(error, status == ERR_BAD_SHAPE, &
                    "shape mismatch against the grid surfaces at create_finalize as BAD_SHAPE")
      end block checks
      if (c_associated(handle)) status = rdb_ocean_destroy(handle)
   end subroutine test_bad_shape_deferred

   subroutine test_topology_injection(error)
      !! `ocean_bc_state_set_topology` (P2.5): forces per-dimension
      !! periodicity and back-fills the edge tags on the periodic axis
      !! ONLY -- the non-periodic axis keeps whatever it already had.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_bc_state_t) :: bc
      integer :: ierr

      call grid%init(ONX, ONY, ONGHOST, ODX, ODY)
      call ocean_bc_state_init(bc, grid, nz_ml=ONZ)

      checks: block
         call check(error, bc%west%bc_type == OBC_WALL .and. bc%south%bc_type == OBC_WALL, &
                    "fresh bc state defaults to WALL on every edge")
         if (allocated(error)) exit checks

         call ocean_bc_state_set_topology(bc, periodic_x=.true., periodic_y=.false., ierr=ierr)
         call check(error, ierr == OCEAN_STATUS_OK, "set_topology(periodic_x=T, periodic_y=F) succeeds")
         if (allocated(error)) exit checks

         call check(error, bc%periodic_x, "periodic_x set")
         if (allocated(error)) exit checks
         call check(error,.not. bc%periodic_y, "periodic_y left false")
         if (allocated(error)) exit checks
         call check(error, bc%west%bc_type == OBC_PERIODIC, "west edge back-filled to OBC_PERIODIC")
         if (allocated(error)) exit checks
         call check(error, bc%east%bc_type == OBC_PERIODIC, "east edge back-filled to OBC_PERIODIC")
         if (allocated(error)) exit checks
         call check(error, bc%south%bc_type == OBC_WALL, &
                    "south edge UNTOUCHED (non-periodic axis keeps its existing physical BC)")
         if (allocated(error)) exit checks
         call check(error, bc%north%bc_type == OBC_WALL, &
                    "north edge UNTOUCHED (non-periodic axis keeps its existing physical BC)")
      end block checks
      call ocean_bc_state_destroy(bc)
   end subroutine test_topology_injection

   subroutine test_topology_injection_thin_halo(error)
      !! Periodic topology on a halo too thin for the PPM/biharmonic
      !! stencil (nghost < 3) is rejected fail-loud, mirroring
      !! `ocean_bc_validate_periodic`'s rule (c) -- this entry point
      !! bypasses that routine so must enforce the rule itself.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_bc_state_t) :: bc
      integer :: ierr

      call grid%init(ONX, ONY, 2, ODX, ODY)   ! nghost = 2 < 3
      call ocean_bc_state_init(bc, grid, nz_ml=ONZ)

      call ocean_bc_state_set_topology(bc, periodic_x=.true., periodic_y=.false., ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_ERR_SETUP, &
                 "periodic topology with nghost < 3 is rejected fail-loud")
      call ocean_bc_state_destroy(bc)
   end subroutine test_topology_injection_thin_halo

   subroutine test_required_halo_minimums(error)
      !! The folded query, checked against every one of the six scattered
      !! rules it collects.
      type(error_type), allocatable, intent(out) :: error

      call check(error, required_halo() == 2, "baseline (nothing selected) is 2")
      if (allocated(error)) return
      call check(error, required_halo(pv_adv_scheme="weno5") >= 3, "pv_adv weno5 >= 3")
      if (allocated(error)) return
      call check(error, required_halo(pv_adv_scheme="weno7") >= 4, "pv_adv weno7 >= 4")
      if (allocated(error)) return
      call check(error, required_halo(tracer_recon="weno5") >= 3, "tracer weno5 >= 3")
      if (allocated(error)) return
      call check(error, required_halo(tracer_recon="weno7") >= 4, "tracer weno7 >= 4")
      if (allocated(error)) return
      call check(error, required_halo(tracer_recon="weno9") >= 5, "tracer weno9 >= 5")
      if (allocated(error)) return
      call check(error, required_halo(periodic=.true.) >= 3, "periodic topology >= 3")
      if (allocated(error)) return
      call check(error, required_halo(tripolar_fold=.true.) >= 3, "tripolar fold >= 3")
      if (allocated(error)) return
      call check(error, required_halo(decomposed=.true.) >= 3, "MPI-decomposed >= 3")
      if (allocated(error)) return
      call check(error, required_halo(kappa_shear_at_vertex=.true.) >= 2, "kappa-shear vertex >= 2")
      if (allocated(error)) return
      call check(error, required_halo(pv_adv_scheme="weno7", tracer_recon="weno9") == 5, &
                 "max() over multiple simultaneous rules")
   end subroutine test_required_halo_minimums

   subroutine test_required_halo_abi_parity(error)
      !! The stateless C ABI entry point must agree with the Fortran
      !! function it wraps.
      type(error_type), allocatable, intent(out) :: error
      integer(c_int) :: ng5, ng7, ng_none

      ! NOTE: the placeholder character(len=1) "x" below (rather than a
      ! literal "") is deliberate -- nvfortran rejects a genuinely
      ! zero-length character LITERAL sequence-associated onto an
      ! explicit-shape character(kind=c_char) array dummy
      ! (NVFORTRAN-S-0188 "type mismatch"), even though the length
      ! ACTUALLY read is governed entirely by the separate `_len`
      ! argument below (0 here) -- c_to_f_string only ever touches
      ! `_len` elements of the array, so the placeholder's content is
      ! never read.
      ng5 = rdb_ocean_required_halo("weno5", 5_c_int, "x", 0_c_int, 0_c_int, 0_c_int, 0_c_int, 0_c_int)
      ng7 = rdb_ocean_required_halo("weno7", 5_c_int, "x", 0_c_int, 0_c_int, 0_c_int, 0_c_int, 0_c_int)
      ng_none = rdb_ocean_required_halo("x", 0_c_int, "x", 0_c_int, 0_c_int, 0_c_int, 0_c_int, 0_c_int)

      call check(error, int(ng5) == required_halo(pv_adv_scheme="weno5"), &
                 "ABI weno5 matches Fortran required_halo")
      if (allocated(error)) return
      call check(error, int(ng7) == required_halo(pv_adv_scheme="weno7"), &
                 "ABI weno7 matches Fortran required_halo")
      if (allocated(error)) return
      call check(error, int(ng_none) == required_halo(), "ABI baseline matches Fortran required_halo")
   end subroutine test_required_halo_abi_parity

end module test_ocean_geometry_inject
