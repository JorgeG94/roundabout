!! Unit tests for the ocean solver-creation path's optional `ierr` contract
!! (`rdb_ocean_status.F90`, P0 of the Python runtime API plan).  Closes the
!! P0.1 adversarial-review finding F3: before this file, none of the ~99
!! `ierr`-present branches P0 added were exercised by the test suite — the
!! 177/177-green result P0 shipped with only ever ran the ABSENT half of
!! every new branch.
!!
!! Also serves as the F1 regression gate: `metrics_fill_from_supergrid` on a
!! missing NetCDF file, called with `ierr` PRESENT, must RETURN
!! `OCEAN_STATUS_ERR_IO` — not `error stop` the process.  Before the F1 fix,
!! this test would abort the whole ctest binary instead of failing cleanly.
!!
!! NetCDF I/O (the F1 gate opens/reads a file) — must run with
!! OMP_NUM_THREADS=1, same as every other NetCDF-touching test.
module test_ocean_status_returns
   use, intrinsic :: iso_c_binding, only: c_associated, c_char, c_int, c_null_ptr, c_ptr
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_config, only: config_t, read_config_from_string, validate_config
   use rdb_eos, only: eos_t, eos_validate
   use rdb_ocean_metrics, only: ocean_metrics_t, metrics_fill_from_supergrid
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_CONFIG_PARSE, &
                               OCEAN_STATUS_ERR_CONFIG_VALIDATE, OCEAN_STATUS_ERR_SETUP, &
                               OCEAN_STATUS_ERR_IO
   use rdb_error_ring, only: error_ring_clear, error_ring_get, error_ring_count, &
                             ERROR_RING_MSG_LEN
   implicit none
   private

   public :: collect_ocean_status_returns_tests

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

      function rdb_ocean_destroy(c_handle) result(status) &
         bind(c, name="rdb_ocean_destroy")
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), intent(inout) :: c_handle
         integer(c_int) :: status
      end function rdb_ocean_destroy

      function rdb_ocean_last_error(idx, buf, cap) result(len_out) &
         bind(c, name="rdb_ocean_last_error")
         import :: c_char, c_int
         implicit none
         integer(c_int), intent(in), value :: idx
         integer(c_int), intent(in), value :: cap
         character(kind=c_char), intent(out) :: buf(cap)
         integer(c_int) :: len_out
      end function rdb_ocean_last_error
   end interface

   integer(c_int), parameter :: ERR_IO = 5

contains

   subroutine collect_ocean_status_returns_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("status_bad_nml_returns_config_parse", &
                               test_bad_nml_returns_config_parse), &
                  new_unittest("status_bad_crossknob_returns_config_validate", &
                               test_bad_crossknob_returns_config_validate), &
                  new_unittest("status_bad_eos_variant_returns_setup", &
                               test_bad_eos_variant_returns_setup), &
                  new_unittest("status_missing_supergrid_file_returns_io_no_abort", &
                               test_missing_supergrid_file_returns_io), &
                  new_unittest("status_good_config_returns_ok", &
                               test_good_config_returns_ok), &
                  new_unittest("status_ring_carries_specific_message", &
                               test_ring_carries_specific_message), &
                  new_unittest("status_ring_survives_threading_wrapper", &
                               test_ring_survives_threading_wrapper) &
                  ]
   end subroutine collect_ocean_status_returns_tests

   subroutine test_bad_nml_returns_config_parse(error)
      !! `read_config_from_string` with `ierr` PRESENT on an unknown
      !! namelist key must return `ERR_CONFIG_PARSE`, not `error stop`.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      character(len=:), allocatable :: nml
      integer :: ierr

      nml = "&grid_nml nx = 8, ny = 6, totally_bogus_key = 1 /"//new_line("a")
      ierr = -999
      call read_config_from_string(nml, cfg, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_ERR_CONFIG_PARSE, &
                 "unknown namelist key returns ERR_CONFIG_PARSE")
   end subroutine test_bad_nml_returns_config_parse

   subroutine test_bad_crossknob_returns_config_validate(error)
      !! `validate_config` with `ierr` PRESENT on a semantically-invalid
      !! cross-knob (nx < 1) must return `ERR_CONFIG_VALIDATE`, not
      !! `error stop`.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      integer :: ierr

      cfg%nx = 0
      ierr = -999
      call validate_config(cfg, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_ERR_CONFIG_VALIDATE, &
                 "nx < 1 returns ERR_CONFIG_VALIDATE")
   end subroutine test_bad_crossknob_returns_config_validate

   subroutine test_bad_eos_variant_returns_setup(error)
      !! `eos_validate` with `ierr` PRESENT on an unrecognised
      !! `eos%variant` must return `ERR_SETUP`, not `error stop` — the
      !! `configure_ocean_*` setup-chain stage of the status taxonomy.
      type(error_type), allocatable, intent(out) :: error
      type(eos_t) :: eos
      integer :: ierr

      eos%variant = -999
      ierr = -999
      call eos_validate(eos, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_ERR_SETUP, &
                 "unknown eos%variant returns ERR_SETUP")
   end subroutine test_bad_eos_variant_returns_setup

   subroutine test_missing_supergrid_file_returns_io(error)
      !! F1 regression gate: `metrics_fill_from_supergrid` on a missing
      !! file, with `ierr` PRESENT, must RETURN `OCEAN_STATUS_ERR_IO`.
      !! Before the F1 fix, `nc_check` (`rdb_io_netcdf.F90`) sat under
      !! this call with no status path and `error stop`-ed unconditionally
      !! — this test would have aborted the whole ctest binary rather than
      !! failing cleanly.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      integer :: ierr

      call grid%init(8, 6, 3, 1.0_wp, 1.0_wp)
      call metrics%init(grid)
      ierr = -999
      call metrics_fill_from_supergrid(metrics, grid, &
                                       "tests/does_not_exist_p0_1_f1_gate.nc", ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_ERR_IO, &
                 "missing supergrid file returns ERR_IO instead of aborting")
   end subroutine test_missing_supergrid_file_returns_io

   subroutine test_good_config_returns_ok(error)
      !! The success path: a well-formed namelist string + a valid config
      !! must both report `OCEAN_STATUS_OK` through `ierr`, not merely
      !! "no abort" — the ABSENT-`ierr` legacy tests never assert the
      !! numeric value on success either.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      character(len=:), allocatable :: nml
      integer :: ierr

      nml = '&sim_nml sim_type = "ocean" /'//new_line("a")// &
            "&grid_nml nx = 8, ny = 6, dx = 2000.0, dy = 2000.0 /"//new_line("a")// &
            "&nonhydrostatic_nml nz_layers = 3 /"//new_line("a")// &
            "&time_nml t_end = 86400.0, dt_fixed = 300.0 /"//new_line("a")
      ierr = -999
      call read_config_from_string(nml, cfg, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "well-formed namelist parses to OCEAN_STATUS_OK")
      if (allocated(error)) return

      ierr = -999
      call validate_config(cfg, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "well-formed config validates to OCEAN_STATUS_OK")
   end subroutine test_good_config_returns_ok

   subroutine test_ring_carries_specific_message(error)
      !! P2 error ring: `fail()` (called internally by `eos_validate` on an
      !! unknown `eos%variant`, migrated in this same phase) must push the
      !! SPECIFIC reason onto the ring, retrievable at index 0 — not just
      !! set the coarse `ierr` code. `RuntimeError: error 3` is exactly
      !! what this ring exists to prevent (D4.1,
      !! `06_python_surface_design.md`).
      type(error_type), allocatable, intent(out) :: error
      type(eos_t) :: eos
      integer :: ierr
      character(len=ERROR_RING_MSG_LEN) :: msg

      call error_ring_clear()
      call check(error, error_ring_count() == 0, "ring starts empty after clear")
      if (allocated(error)) return

      eos%variant = -999
      ierr = -999
      call eos_validate(eos, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_ERR_SETUP, "unknown eos%variant returns ERR_SETUP")
      if (allocated(error)) return

      call check(error, error_ring_count() >= 1, "the failure pushed at least one ring entry")
      if (allocated(error)) return
      msg = error_ring_get(0)
      call check(error, index(msg, "eos%variant") > 0, &
                 "ring[0] carries the SPECIFIC reason ('unknown eos%variant'), "// &
                 "not just the coarse ERR_SETUP code")
   end subroutine test_ring_carries_specific_message

   subroutine test_ring_survives_threading_wrapper(error)
      !! P2 error ring, the 23-sites case (P0 review F2): the missing-
      !! supergrid-file path threads `ierr` from `metrics_fill_from_supergrid`
      !! up through the `configure_ocean_metrics` wrapper — which, per the
      !! ALREADY-FIXED F2 pattern, does NOT push its own generic message
      !! when `ierr` is present (it only re-raises the callee's code). So
      !! the ring's most recent entry must be the SPECIFIC NetCDF reason
      !! ("NetCDF error in opening ...: No such file or directory"), not a
      !! generic "configure_ocean_metrics failed" wrapper sentence — this
      !! is exactly the case a single last-error SLOT would have gotten
      !! wrong (an outer wrapper overwriting the specific inner one) and a
      !! RING gets right.
      !!
      !! Drives the real C-ABI (`rdb_ocean_create_from_string` +
      !! `rdb_ocean_last_error`), not the bare Fortran routine, so this
      !! is the actual Python-facing contract, not just the ring module in
      !! isolation.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle
      character(len=:), allocatable :: nml
      character(kind=c_char) :: buf(256)
      integer(c_int) :: status, msg_len
      character(len=256) :: msg
      integer :: i

      nml = '&sim_nml sim_type = "ocean" /'//new_line("a")// &
            "&grid_nml nx = 8, ny = 6, dx = 2000.0, dy = 2000.0 /"//new_line("a")// &
            "&nonhydrostatic_nml nz_layers = 3 /"//new_line("a")// &
            "&time_nml t_end = 86400.0, dt_fixed = 300.0 /"//new_line("a")// &
            "&ocean_topo_nml max_depth = 200.0 /"//new_line("a")// &
            "&ocean_bt_nml auto_n_inner = .true. /"//new_line("a")// &
            "&tracer_nml initial_salinity = 35.0, initial_temperature = 12.0 /"//new_line("a")// &
            "&ocean_diag_nml enabled = .false. /"//new_line("a")// &
            "&output_nml output_to_file = .false. /"//new_line("a")// &
            '&ocean_grid_nml grid_config = "supergrid", '// &
            'supergrid_file = "tests/does_not_exist_p0_1_f1_gate.nc" /'//new_line("a")

      handle = c_null_ptr
      status = rdb_ocean_create_from_string(nml, len(nml, kind=c_int), handle)
      call check(error, status == ERR_IO, "create() with a missing supergrid_file returns ERR_IO")
      if (allocated(error)) return
      call check(error,.not. c_associated(handle), "rejected create leaves handle_out null")
      if (allocated(error)) return

      msg_len = rdb_ocean_last_error(0_c_int, buf, int(len(msg), c_int))
      call check(error, msg_len > 0, "the ring has a message at index 0")
      if (allocated(error)) return

      msg = ""
      do i = 1, min(int(msg_len), len(msg))
         msg(i:i) = buf(i)
      end do

      call check(error, index(msg, "NetCDF") > 0, &
                 "ring[0] is the SPECIFIC NetCDF reason (mentions 'NetCDF'), "// &
                 "surviving the configure_ocean_metrics threading wrapper")
      if (allocated(error)) return
      call check(error, index(msg, "configure_ocean_metrics") == 0, &
                 "ring[0] is NOT the generic outer-wrapper text — "// &
                 "the F2 regression this ring design exists to close")
      if (allocated(error)) return

      ! An out-of-range index is "nothing there", not an error/crash.
      msg_len = rdb_ocean_last_error(999_c_int, buf, int(len(msg), c_int))
      call check(error, msg_len == 0_c_int, "an out-of-range ring index returns length 0")
   end subroutine test_ring_survives_threading_wrapper

end module test_ocean_status_returns
