!! Unit tests for the ocean C ABI lifecycle skeleton (Python runtime API
!! plan, P1: `src/api/rdb_handle.F90` + `src/api/rdb_ocean_api.F90`).
!!
!! Drives the SAME `bind(c)` entry points a Python/C caller uses: in-memory
!! namelist create, fixed-dt stepping, the scalar introspection getters, and
!! destroy — plus the single-live-handle guard and the "no `error stop` on a
!! bad config" contract. Modelled on the pre-carve-out precedent
!! (`tmp_local_artifacts/python_ffi_scope/recovered/tests/test_rdb_api_ocean.F90`),
!! trimmed to this phase's lifecycle-only surface (no state-array accessors —
!! those are P2).
!!
!! Own ctest binary: handle lifecycle + the single-live-handle module state
!! (`rdb_ocean_api`'s `g_handle_live`) must start from a fresh process, same
!! reasoning as the other `rdb_api*` suites in the recovered precedent.
module test_rdb_ocean_api
   use, intrinsic :: iso_c_binding, only: c_associated, c_char, c_double, c_f_pointer, c_int, &
                                                                             c_null_ptr, c_ptr
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   implicit none
   private

   public :: collect_rdb_ocean_api_tests

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

      function rdb_ocean_step(c_handle, n_steps) result(status) &
         bind(c, name="rdb_ocean_step")
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), intent(in), value :: c_handle
         integer(c_int), intent(in), value :: n_steps
         integer(c_int) :: status
      end function rdb_ocean_step

      function rdb_ocean_destroy(c_handle) result(status) &
         bind(c, name="rdb_ocean_destroy")
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), intent(inout) :: c_handle
         integer(c_int) :: status
      end function rdb_ocean_destroy

      function rdb_ocean_get_time(c_handle, t_out) result(status) &
         bind(c, name="rdb_ocean_get_time")
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), intent(in), value :: c_handle
         real(c_double), intent(out) :: t_out
         integer(c_int) :: status
      end function rdb_ocean_get_time

      function rdb_ocean_get_step_count(c_handle, step_out) result(status) &
         bind(c, name="rdb_ocean_get_step_count")
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), intent(in), value :: c_handle
         integer(c_int), intent(out) :: step_out
         integer(c_int) :: status
      end function rdb_ocean_get_step_count

      function rdb_ocean_get_grid_info(c_handle, nx, ny, nz, nghost) &
         result(status) bind(c, name="rdb_ocean_get_grid_info")
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), intent(in), value :: c_handle
         integer(c_int), intent(out) :: nx, ny, nz, nghost
         integer(c_int) :: status
      end function rdb_ocean_get_grid_info

      function rdb_ocean_get_total_mass(c_handle, m_out) result(status) &
         bind(c, name="rdb_ocean_get_total_mass")
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), intent(in), value :: c_handle
         real(c_double), intent(out) :: m_out
         integer(c_int) :: status
      end function rdb_ocean_get_total_mass

      function rdb_working_precision() result(bytes) &
         bind(c, name="rdb_working_precision")
         import :: c_int
         implicit none
         integer(c_int) :: bytes
      end function rdb_working_precision

      function rdb_ocean_refresh_host(c_handle) result(status) &
         bind(c, name="rdb_ocean_refresh_host")
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), intent(in), value :: c_handle
         integer(c_int) :: status
      end function rdb_ocean_refresh_host

      function rdb_ocean_get_q_heat_ptr(c_handle, ptr, nx, ny, gen) result(status) &
         bind(c, name="rdb_ocean_get_q_heat_ptr")
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), intent(in), value :: c_handle
         type(c_ptr), intent(out) :: ptr
         integer(c_int), intent(out) :: nx, ny, gen
         integer(c_int) :: status
      end function rdb_ocean_get_q_heat_ptr

      function rdb_ocean_get_q_salt_ptr(c_handle, ptr, nx, ny, gen) result(status) &
         bind(c, name="rdb_ocean_get_q_salt_ptr")
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), intent(in), value :: c_handle
         type(c_ptr), intent(out) :: ptr
         integer(c_int), intent(out) :: nx, ny, gen
         integer(c_int) :: status
      end function rdb_ocean_get_q_salt_ptr
   end interface

   ! Status codes (mirrors include/rdb_ocean.h / rdb_ocean_status.F90 —
   ! kept as local literals so this test exercises the ABI's numeric
   ! contract, not just the Fortran module's named parameters).
   integer(c_int), parameter :: OK = 0
   integer(c_int), parameter :: ERR_CONFIG_PARSE = 1
   integer(c_int), parameter :: ERR_CONFIG_VALIDATE = 2
   integer(c_int), parameter :: ERR_IO = 5
   integer(c_int), parameter :: ERR_ALREADY_EXISTS = 11

   ! Tiny quiescent ocean column: flat topo, all-WALL boundaries (default),
   ! zero wind/forcing -> mass must be exactly conserved step to step.
   integer, parameter :: ONX = 8, ONY = 6, ONZ = 3
   real(wp), parameter :: ODT = 300.0_wp

contains

   subroutine collect_rdb_ocean_api_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("ocean_create_step_destroy", test_create_step_destroy), &
                  new_unittest("ocean_destroy_idempotent", test_destroy_idempotent), &
                  new_unittest("ocean_second_create_rejected", test_second_create_rejected), &
                  new_unittest("ocean_create_bad_config", test_create_bad_config), &
                  new_unittest("ocean_create_missing_supergrid_file", &
                               test_create_missing_supergrid_file), &
                  new_unittest("ocean_ice_advances_through_abi", &
                               test_ice_advances_through_abi) &
                  ]
   end subroutine collect_rdb_ocean_api_tests

   function ocean_nml() result(txt)
      !! In-memory namelist for the quiescent ocean column. auto_n_inner is
      !! mandatory: the default n_inner=0 / auto_n_inner=.false. leaves the
      !! barotropic substep count unresolved, which create() rejects
      !! (OCEAN_STATUS_ERR_SETUP) rather than silently running unsplit.
      character(len=:), allocatable :: txt
      txt = '&sim_nml sim_type = "ocean" /'//new_line("a")// &
            "&grid_nml nx = 8, ny = 6, dx = 2000.0, dy = 2000.0 /"//new_line("a")// &
            "&nonhydrostatic_nml nz_layers = 3 /"//new_line("a")// &
            "&time_nml t_end = 86400.0, dt_fixed = 300.0 /"//new_line("a")// &
            "&ocean_topo_nml max_depth = 200.0 /"//new_line("a")// &
            "&ocean_bt_nml auto_n_inner = .true. /"//new_line("a")// &
            "&tracer_nml initial_salinity = 35.0, initial_temperature = 12.0 /"//new_line("a")// &
            "&ocean_diag_nml enabled = .false. /"//new_line("a")// &
            "&output_nml output_to_file = .false. /"//new_line("a")
   end function ocean_nml

   subroutine test_create_step_destroy(error)
      !! One-shot create -> step -> introspect -> destroy. Asserts time and
      !! step count advance by exactly the expected amount, total mass stays
      !! finite and (quiescent closed basin) is conserved to machine
      !! precision, and grid/precision queries report the configured shape.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle
      character(len=:), allocatable :: nml
      integer(c_int) :: status, nx, ny, nz, nghost, step_out, wp_bytes
      real(c_double) :: t_out, mass0, mass1

      handle = c_null_ptr
      nml = ocean_nml()
      status = rdb_ocean_create_from_string(nml, len(nml, kind=c_int), handle)
      checks: block
         call check(error, status == OK, "create_from_string succeeds")
         if (allocated(error)) exit checks

         wp_bytes = rdb_working_precision()
         call check(error, wp_bytes == 4_c_int .or. wp_bytes == 8_c_int, &
                    "working precision is 4 or 8 bytes")
         if (allocated(error)) exit checks

         status = rdb_ocean_get_grid_info(handle, nx, ny, nz, nghost)
         call check(error, status == OK, "get_grid_info succeeds")
         if (allocated(error)) exit checks
         call check(error, nx == int(ONX, c_int) .and. ny == int(ONY, c_int) &
                    .and. nz == int(ONZ, c_int), "grid shape matches the namelist")
         if (allocated(error)) exit checks

         status = rdb_ocean_get_time(handle, t_out)
         call check(error, status == OK .and. abs(t_out) < 1.0e-12_wp, &
                    "time starts at 0")
         if (allocated(error)) exit checks
         status = rdb_ocean_get_step_count(handle, step_out)
         call check(error, status == OK .and. step_out == 0_c_int, &
                    "step_count starts at 0")
         if (allocated(error)) exit checks

         status = rdb_ocean_get_total_mass(handle, mass0)
         call check(error, status == OK, "get_total_mass succeeds (pre-step)")
         if (allocated(error)) exit checks
         call check(error, ieee_is_finite(real(mass0, wp)) .and. mass0 > 0.0_c_double, &
                    "initial total mass is finite and positive")
         if (allocated(error)) exit checks

         status = rdb_ocean_step(handle, 4_c_int)
         call check(error, status == OK, "step(4) succeeds")
         if (allocated(error)) exit checks

         status = rdb_ocean_get_time(handle, t_out)
         call check(error, status == OK .and. abs(t_out - 4.0_wp*ODT) < 1.0e-9_wp, &
                    "time == 4*dt after step(4)")
         if (allocated(error)) exit checks
         status = rdb_ocean_get_step_count(handle, step_out)
         call check(error, status == OK .and. step_out == 4_c_int, &
                    "step_count == 4 after step(4)")
         if (allocated(error)) exit checks

         status = rdb_ocean_get_total_mass(handle, mass1)
         call check(error, status == OK, "get_total_mass succeeds (post-step)")
         if (allocated(error)) exit checks
         call check(error, ieee_is_finite(real(mass1, wp)), "post-step mass is finite")
         if (allocated(error)) exit checks
         ! Quiescent closed (all-WALL) basin, zero wind/forcing: mass is
         ! conserved to within floating-point roundoff, not merely "close".
         call check(error, abs(mass1 - mass0) < 1.0e-6_wp*abs(mass0), &
                    "total mass conserved across steps in a quiescent closed basin")
         if (allocated(error)) exit checks

         ! step(0) / step(negative) are no-ops, not errors.
         status = rdb_ocean_step(handle, 0_c_int)
         call check(error, status == OK, "step(0) is a no-op success")
         if (allocated(error)) exit checks
         status = rdb_ocean_get_step_count(handle, step_out)
         call check(error, step_out == 4_c_int, "step(0) does not advance step_count")
      end block checks

      status = rdb_ocean_destroy(handle)
      call check(error, status == OK, "destroy succeeds")
      call check(error,.not. c_associated(handle), "destroy nulls the handle")
   end subroutine test_create_step_destroy

   subroutine test_destroy_idempotent(error)
      !! destroy() on an already-destroyed handle is a no-op success, and on
      !! a handle that was never created (null c_ptr) likewise — the
      !! contract a Python __del__ relies on.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle
      character(len=:), allocatable :: nml
      integer(c_int) :: status

      handle = c_null_ptr
      nml = ocean_nml()
      status = rdb_ocean_create_from_string(nml, len(nml, kind=c_int), handle)
      call check(error, status == OK, "create succeeds")
      if (allocated(error)) return

      status = rdb_ocean_destroy(handle)
      call check(error, status == OK, "first destroy succeeds")
      if (allocated(error)) return

      status = rdb_ocean_destroy(handle)
      call check(error, status == OK, "second destroy on the same (now-null) handle is a no-op")
      if (allocated(error)) return

      handle = c_null_ptr
      status = rdb_ocean_destroy(handle)
      call check(error, status == OK, "destroy on a handle that was never created is a no-op")
   end subroutine test_destroy_idempotent

   subroutine test_second_create_rejected(error)
      !! Multi-instance is impossible for the ocean path (module-scope
      !! topology/device state in rdb_ocean_halo et al. — see
      !! rdb_ocean_api's header comment). A second create while one handle
      !! is live must fail with the distinct ALREADY_EXISTS code, leave the
      !! first handle untouched, and leave *handle_out null — never silently
      !! corrupt the first handle. After destroying the first, a second
      !! create must succeed (the guard lifts).
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: h1, h2
      character(len=:), allocatable :: nml
      integer(c_int) :: status, step_out

      h1 = c_null_ptr
      h2 = c_null_ptr
      nml = ocean_nml()
      status = rdb_ocean_create_from_string(nml, len(nml, kind=c_int), h1)
      checks: block
         call check(error, status == OK, "first create succeeds")
         if (allocated(error)) exit checks

         status = rdb_ocean_create_from_string(nml, len(nml, kind=c_int), h2)
         call check(error, status == ERR_ALREADY_EXISTS, &
                    "second create while one is live returns ALREADY_EXISTS")
         if (allocated(error)) exit checks
         call check(error,.not. c_associated(h2), "rejected create leaves handle_out null")
         if (allocated(error)) exit checks

         ! First handle still fully usable.
         status = rdb_ocean_step(h1, 1_c_int)
         call check(error, status == OK, "first handle still steppable after a rejected create")
         if (allocated(error)) exit checks
         status = rdb_ocean_get_step_count(h1, step_out)
         call check(error, status == OK .and. step_out == 1_c_int, &
                    "first handle's state is untouched by the rejected create")
      end block checks

      status = rdb_ocean_destroy(h1)
      call check(error, status == OK, "destroy of the first handle succeeds")
      if (allocated(error)) return

      ! Guard lifts once the live handle is gone.
      status = rdb_ocean_create_from_string(nml, len(nml, kind=c_int), h2)
      call check(error, status == OK, "create succeeds again after the live handle is destroyed")
      if (allocated(error)) return
      status = rdb_ocean_destroy(h2)
      call check(error, status == OK, "cleanup: destroy the second handle")
   end subroutine test_second_create_rejected

   subroutine test_create_bad_config(error)
      !! A malformed namelist (unknown key -> strict schema parse failure)
      !! or a semantically-invalid one (dt_fixed <= 0, which create() checks
      !! itself, matching `bench_ocean`'s own guard) must return a non-zero
      !! status — NOT abort the process. This is the whole point of the P0
      !! `ierr` plumbing: create() is a library entry point, not a driver,
      !! and cannot call `error stop` on bad user input.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle
      character(len=:), allocatable :: nml
      integer(c_int) :: status

      ! Unknown key -> strict schema parse failure.
      handle = c_null_ptr
      nml = "&grid_nml nx = 8, ny = 6, totally_bogus_key = 1 /"//new_line("a")
      status = rdb_ocean_create_from_string(nml, len(nml, kind=c_int), handle)
      call check(error, status == ERR_CONFIG_PARSE, &
                 "unknown namelist key returns ERR_CONFIG_PARSE, not an abort")
      if (allocated(error)) return
      call check(error,.not. c_associated(handle), "bad-config create leaves handle_out null")
      if (allocated(error)) return

      ! Semantically invalid: dt_fixed defaults to 0.0 and is never set here
      ! -> create()'s own guard (mirrors bench_ocean's `dt_fixed <= 0` check).
      handle = c_null_ptr
      nml = "&grid_nml nx = 8, ny = 6, dx = 2000.0, dy = 2000.0 /"//new_line("a")
      status = rdb_ocean_create_from_string(nml, len(nml, kind=c_int), handle)
      call check(error, status == ERR_CONFIG_VALIDATE, &
                 "dt_fixed <= 0 returns ERR_CONFIG_VALIDATE, not an abort")
      if (allocated(error)) return
      call check(error,.not. c_associated(handle), "bad-config create leaves handle_out null")
      if (allocated(error)) return

      ! A good create still works after two rejected ones -- the module-
      ! level single-handle guard was never engaged by a failed create.
      handle = c_null_ptr
      nml = ocean_nml()
      status = rdb_ocean_create_from_string(nml, len(nml, kind=c_int), handle)
      call check(error, status == OK, "a good create still succeeds after prior rejections")
      if (allocated(error)) return
      status = rdb_ocean_destroy(handle)
      call check(error, status == OK, "cleanup destroy")
   end subroutine test_create_bad_config

   subroutine test_create_missing_supergrid_file(error)
      !! End-to-end proof for P0.1 review F1: `create()` with
      !! `grid_config="supergrid"` pointing at a NetCDF file that does not
      !! exist must return `ERR_IO` — NOT abort the host process. Before
      !! the F1 fix, `nc_check` (`rdb_io_netcdf.F90`) sat under
      !! `metrics_fill_from_supergrid` with no status path and
      !! unconditionally `error stop`-ped, so this exact scenario (a wrong
      !! file path — by far the most likely real Python-caller mistake)
      !! would have killed the whole process instead of returning a status
      !! this test could even observe.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle
      character(len=:), allocatable :: nml
      integer(c_int) :: status

      handle = c_null_ptr
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
      status = rdb_ocean_create_from_string(nml, len(nml, kind=c_int), handle)
      call check(error, status == ERR_IO, &
                 "create() with a nonexistent supergrid_file returns ERR_IO, not an abort")
      if (allocated(error)) return
      call check(error,.not. c_associated(handle), "rejected create leaves handle_out null")
      if (allocated(error)) return

      ! The single-live-handle guard must not have latched from the failed
      ! create -- a good create right after must still succeed.
      handle = c_null_ptr
      nml = ocean_nml()
      status = rdb_ocean_create_from_string(nml, len(nml, kind=c_int), handle)
      call check(error, status == OK, "a good create still succeeds after the ERR_IO rejection")
      if (allocated(error)) return
      status = rdb_ocean_destroy(handle)
      call check(error, status == OK, "cleanup destroy")
   end subroutine test_create_missing_supergrid_file

   subroutine test_ice_advances_through_abi(error)
      !! P2.4b regression gate: before this phase, `engine_step`/
      !! `engine_step_finalize` never ran the sea-ice per-step block, so a
      !! caller that enabled ice via `engine_setup` got it initialised but
      !! never advanced ("configured but not executed" — the exact bug
      !! class `engine_step_ice` closes, see
      !! `src/core/ocean/rdb_ocean_engine.F90`'s module docstring).
      !! `rdb_ocean_step` is the ONLY C-ABI entry point that advances
      !! time, so this drives the fix end-to-end through the ABI, not by
      !! calling into any Fortran-internal engine routine directly.
      !!
      !! `initial_temperature` is seeded BELOW the salinity-dependent
      !! freezing point, so `ice_frazil_accumulate` deterministically banks
      !! heat on the very first outer step (no multi-day cooling spin-up
      !! needed — see `rdb_ice_frazil.F90`) and, since `dt_therm_ratio`
      !! defaults to 1 (thermo cadence fires every step), the thermo chain
      !! inside `engine_step_ice` spends that bank on the SAME step,
      !! feeding real ice-derived salt/heat back into the ocean via
      !! `ice_ocean_brine_flux`/`ice_ocean_heat_flux` — both of which
      !! OVERWRITE `Q_salt`/`Q_heat` (`rdb_ice_ocean_coupler.F90`).
      !!
      !! Q_heat/Q_salt (not ice conc/thickness directly — those have no
      !! C-ABI accessor yet) are the ABI-visible witness: with ice enabled
      !! and `&ocean_forcing_nml enable_components` at its default off,
      !! NOTHING else in the step path ever touches these two fields after
      !! `engine_setup` seeds them from the constant `&ocean_thermo_nml`
      !! knobs — `ocean_surface_flux_assemble` is a no-op without
      !! components. So pre-fix, Q_heat/Q_salt are bit-identical to their
      !! setup value after any number of `rdb_ocean_step` calls;
      !! post-fix they change on step 1. Each read goes through
      !! `rdb_ocean_refresh_host` first — the D<->H contract
      !! (CLAUDE.md, `docs/ocean_python_api_plan.md` S3): a getter alone
      !! never triggers a device->host copy, so skipping the refresh would
      !! read stale host memory on the GPU build and could pass/fail for
      !! the wrong reason.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle
      character(len=:), allocatable :: nml
      integer(c_int) :: status
      type(c_ptr) :: qptr
      integer(c_int) :: qnx, qny, qgen
      real(wp), pointer :: q2(:, :)
      real(wp) :: q_heat_0, q_salt_0, q_heat_1, q_salt_1

      handle = c_null_ptr
      nml = '&sim_nml sim_type = "ocean" /'//new_line("a")// &
            "&grid_nml nx = 8, ny = 6, dx = 2000.0, dy = 2000.0 /"//new_line("a")// &
            "&nonhydrostatic_nml nz_layers = 3 /"//new_line("a")// &
            "&time_nml t_end = 86400.0, dt_fixed = 300.0 /"//new_line("a")// &
            "&ocean_topo_nml max_depth = 200.0 /"//new_line("a")// &
            "&ocean_bt_nml auto_n_inner = .true. /"//new_line("a")// &
            ! S=34 psu -> T_f ~ -1.9 degC; -3.0 is comfortably below it so
            ! frazil banks on step 1 regardless of the exact EOS used.
            "&tracer_nml initial_salinity = 34.0, initial_temperature = -3.0 /"//new_line("a")// &
            "&ocean_ice_nml enable = .true. /"//new_line("a")// &
            "&ocean_diag_nml enabled = .false. /"//new_line("a")// &
            "&output_nml output_to_file = .false. /"//new_line("a")
      status = rdb_ocean_create_from_string(nml, len(nml, kind=c_int), handle)
      checks: block
         call check(error, status == OK, "ice-enabled create succeeds")
         if (allocated(error)) exit checks

         status = rdb_ocean_refresh_host(handle)
         call check(error, status == OK, "refresh_host succeeds (pre-step)")
         if (allocated(error)) exit checks

         status = rdb_ocean_get_q_heat_ptr(handle, qptr, qnx, qny, qgen)
         call check(error, status == OK, "get_q_heat_ptr succeeds (pre-step)")
         if (allocated(error)) exit checks
         call c_f_pointer(qptr, q2, [int(qnx), int(qny)])
         q_heat_0 = sum(q2)/real(size(q2), wp)

         status = rdb_ocean_get_q_salt_ptr(handle, qptr, qnx, qny, qgen)
         call check(error, status == OK, "get_q_salt_ptr succeeds (pre-step)")
         if (allocated(error)) exit checks
         call c_f_pointer(qptr, q2, [int(qnx), int(qny)])
         q_salt_0 = sum(q2)/real(size(q2), wp)

         status = rdb_ocean_step(handle, 2_c_int)
         call check(error, status == OK, "step(2) succeeds with ice enabled")
         if (allocated(error)) exit checks

         status = rdb_ocean_refresh_host(handle)
         call check(error, status == OK, "refresh_host succeeds (post-step)")
         if (allocated(error)) exit checks

         status = rdb_ocean_get_q_heat_ptr(handle, qptr, qnx, qny, qgen)
         call check(error, status == OK, "get_q_heat_ptr succeeds (post-step)")
         if (allocated(error)) exit checks
         call c_f_pointer(qptr, q2, [int(qnx), int(qny)])
         q_heat_1 = sum(q2)/real(size(q2), wp)
         call check(error, ieee_is_finite(q_heat_1), "post-step Q_heat is finite")
         if (allocated(error)) exit checks

         status = rdb_ocean_get_q_salt_ptr(handle, qptr, qnx, qny, qgen)
         call check(error, status == OK, "get_q_salt_ptr succeeds (post-step)")
         if (allocated(error)) exit checks
         call c_f_pointer(qptr, q2, [int(qnx), int(qny)])
         q_salt_1 = sum(q2)/real(size(q2), wp)
         call check(error, ieee_is_finite(q_salt_1), "post-step Q_salt is finite")
         if (allocated(error)) exit checks

         ! The regression gate: pre-P2.4b this never changes (engine_step_ice
         ! did not exist, so nothing ever wrote these fields after setup).
         call check(error, abs(q_heat_1 - q_heat_0) > 1.0e-10_wp .or. &
                    abs(q_salt_1 - q_salt_0) > 1.0e-10_wp, &
                    "sea-ice per-step physics mediates a change into Q_heat/Q_salt "// &
                    "through rdb_ocean_step (engine_step_ice ran)")
      end block checks

      status = rdb_ocean_destroy(handle)
      call check(error, status == OK, "cleanup destroy")
   end subroutine test_ice_advances_through_abi

end module test_rdb_ocean_api
