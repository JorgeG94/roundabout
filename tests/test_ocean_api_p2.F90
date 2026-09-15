!! Unit tests for the P2 accessor surface (`src/api/rdb_ocean_api.F90`):
!! the D<->H contract's regression gate + the lifetime fact the whole
!! silent-re-sync design rests on. See `docs/ocean_python_api_plan.md` S3
!! and `tmp_local_artifacts/python_ffi_scope/06_python_surface_design.md`
!! D3.4/D3.5.
!!
!! CAVEAT, stated where it matters (CLAUDE.md, docs/ocean_python_api_plan.md
!! S6): on a host-only build every `!$acc` directive is an inert comment,
!! so BOTH tests here pass on this build by construction and prove nothing
!! about the GPU `mem:separate` path — see each test's docstring for what
!! it would actually be catching on a GPU build.
!!
!! Own ctest binary: the single-live-handle module state
!! (`rdb_ocean_api`'s `g_handle_live`) must start from a fresh process,
!! same reasoning as `test_rdb_ocean_api`.
module test_ocean_api_p2
   use, intrinsic :: iso_c_binding, only: c_associated, c_char, c_double, c_f_pointer, &
                                                                             c_int, c_null_ptr, c_ptr
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   implicit none
   private

   public :: collect_ocean_api_p2_tests

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

      function rdb_ocean_refresh_host(c_handle) result(status) &
         bind(c, name="rdb_ocean_refresh_host")
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), intent(in), value :: c_handle
         integer(c_int) :: status
      end function rdb_ocean_refresh_host

      function rdb_ocean_get_h_layer_ptr(c_handle, ptr, nx, ny, nz, gen) &
         result(status) bind(c, name="rdb_ocean_get_h_layer_ptr")
         import :: c_int, c_ptr
         implicit none
         type(c_ptr), intent(in), value :: c_handle
         type(c_ptr), intent(out) :: ptr
         integer(c_int), intent(out) :: nx, ny, nz, gen
         integer(c_int) :: status
      end function rdb_ocean_get_h_layer_ptr

      function rdb_ocean_set_wind(c_handle, taux_data, tauy_data, nx_p, ny_p) &
         result(status) bind(c, name="rdb_ocean_set_wind")
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), intent(in), value :: c_handle
         integer(c_int), intent(in), value :: nx_p, ny_p
         real(c_double), intent(in) :: taux_data(nx_p + 1, ny_p)
         real(c_double), intent(in) :: tauy_data(nx_p, ny_p + 1)
         integer(c_int) :: status
      end function rdb_ocean_set_wind

      function rdb_ocean_set_bathymetry(c_handle, b_data, nx_p, ny_p) &
         result(status) bind(c, name="rdb_ocean_set_bathymetry")
         import :: c_double, c_int, c_ptr
         implicit none
         type(c_ptr), intent(in), value :: c_handle
         integer(c_int), intent(in), value :: nx_p, ny_p
         real(c_double), intent(in) :: b_data(nx_p, ny_p)
         integer(c_int) :: status
      end function rdb_ocean_set_bathymetry
   end interface

   integer(c_int), parameter :: OK = 0

   ! Same tiny quiescent ocean column as test_rdb_ocean_api.
   integer, parameter :: ONX = 8, ONY = 6, ONZ = 3
   real(wp), parameter :: ODEPTH = 200.0_wp

contains

   subroutine collect_ocean_api_p2_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("ocean_setter_no_clobber", test_setter_no_clobber), &
                  new_unittest("ocean_host_ptr_stable_across_step", &
                               test_host_ptr_stable_across_step) &
                  ]
   end subroutine collect_ocean_api_p2_tests

   function ocean_nml() result(txt)
      !! Same in-memory namelist as test_rdb_ocean_api's — auto_n_inner is
      !! mandatory (the default n_inner=0 leaves the barotropic substep
      !! count unresolved and create() rejects it).
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

   subroutine test_setter_no_clobber(error)
      !! Ported from the recovered `feat/api-ocean-write` precedent
      !! (`tmp_local_artifacts/python_ffi_scope/recovered/tests/
      !! test_rdb_api_ocean_write.F90:378`), the Q9-blocker regression: a
      !! forcing setter called mid-run with NO intervening field read must
      !! NOT rewind live device prognostics to a stale host snapshot. A
      !! BROAD "push everything host->device" would (the old
      !! `ocean_engine_update_device` trap the recovered code documents);
      !! the narrow `rdb_ocean_set_wind` / `rdb_ocean_set_bathymetry`
      !! (pushing ONLY `tau_x/tau_y` / `b`+`bt_H_ref`) do not.
      !!
      !! Adapted for the single-live-handle constraint this tree enforces
      !! (the recovered version held TWO live handles at once, which is
      !! impossible here — `rdb_ocean_api`'s module-scope
      !! `g_handle_live` guard): run the reference trajectory to
      !! completion and SNAPSHOT `h_layer` into an OWNED array before
      !! destroying that handle, then run the no-clobber trajectory in a
      !! fresh handle and compare against the saved snapshot.
      !!
      !! On a GPU `mem:separate` build this is exactly what would fail if
      !! a future "fix" made either setter push more than its own array:
      !! the mid-run re-set would silently re-copy the STALE host
      !! prognostics (never refreshed since the first 10 steps) onto the
      !! device, erasing the second 10 steps' worth of state advance. On
      !! this host build `!$acc update device` is a no-op comment, so the
      !! two trajectories are bit-exact regardless — this test guards a
      !! FUTURE regression in the narrow-push discipline, not a bug being
      !! demonstrated live here.
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle, ptr
      character(len=:), allocatable :: nml
      integer(c_int) :: status, nx, ny, nz, gen
      real(c_double), allocatable :: taux(:, :), tauy(:, :), bflat(:, :)
      real(wp), pointer :: v(:, :, :)
      real(wp), allocatable :: h_ref(:, :, :)

      nml = ocean_nml()
      allocate (taux(ONX + 1, ONY), tauy(ONX, ONY + 1), bflat(ONX, ONY))
      taux = 0.15_c_double
      tauy = 0.0_c_double
      bflat = real(ODEPTH, c_double)

      ! ---- Reference: one clean wind set, then 20 uninterrupted steps. ----
      handle = c_null_ptr
      status = rdb_ocean_create_from_string(nml, len(nml, kind=c_int), handle)
      refblock: block
         call check(error, status == OK, "reference create succeeds")
         if (allocated(error)) exit refblock
         status = rdb_ocean_set_wind(handle, taux, tauy, int(ONX, c_int), int(ONY, c_int))
         call check(error, status == OK, "reference set_wind succeeds")
         if (allocated(error)) exit refblock
         status = rdb_ocean_step(handle, 20_c_int)
         call check(error, status == OK, "reference 20 steps succeed")
         if (allocated(error)) exit refblock
         status = rdb_ocean_refresh_host(handle)
         call check(error, status == OK, "reference refresh_host succeeds")
         if (allocated(error)) exit refblock
         status = rdb_ocean_get_h_layer_ptr(handle, ptr, nx, ny, nz, gen)
         call check(error, status == OK, "reference get_h_layer_ptr succeeds")
         if (allocated(error)) exit refblock
         call c_f_pointer(ptr, v, [int(nx), int(ny), int(nz)])
         h_ref = v   ! owned snapshot -- v itself dangles once handle is destroyed
      end block refblock
      status = rdb_ocean_destroy(handle)
      if (allocated(error)) then
         call check(error, status == OK, "cleanup destroy (reference)")
         return
      end if

      ! ---- No-clobber run: wind, step 10, then RE-SET wind + bathymetry
      ! with NO field read in between (host prognostics are stale at that
      ! point -- host_is_current went false at the first step()), then
      ! step 10 more. ----
      handle = c_null_ptr
      status = rdb_ocean_create_from_string(nml, len(nml, kind=c_int), handle)
      runblock: block
         call check(error, status == OK, "run create succeeds")
         if (allocated(error)) exit runblock
         status = rdb_ocean_set_wind(handle, taux, tauy, int(ONX, c_int), int(ONY, c_int))
         call check(error, status == OK, "run set_wind succeeds")
         if (allocated(error)) exit runblock
         status = rdb_ocean_step(handle, 10_c_int)
         call check(error, status == OK, "run first 10 steps succeed")
         if (allocated(error)) exit runblock

         status = rdb_ocean_set_wind(handle, taux, tauy, int(ONX, c_int), int(ONY, c_int))
         call check(error, status == OK, "mid-run re-set wind (no read) succeeds")
         if (allocated(error)) exit runblock
         status = rdb_ocean_set_bathymetry(handle, bflat, int(ONX, c_int), int(ONY, c_int))
         call check(error, status == OK, "mid-run re-set bathymetry (no read) succeeds")
         if (allocated(error)) exit runblock

         status = rdb_ocean_step(handle, 10_c_int)
         call check(error, status == OK, "run last 10 steps succeed")
         if (allocated(error)) exit runblock

         status = rdb_ocean_refresh_host(handle)
         call check(error, status == OK, "run refresh_host succeeds")
         if (allocated(error)) exit runblock
         status = rdb_ocean_get_h_layer_ptr(handle, ptr, nx, ny, nz, gen)
         call check(error, status == OK, "run get_h_layer_ptr succeeds")
         if (allocated(error)) exit runblock
         call c_f_pointer(ptr, v, [int(nx), int(ny), int(nz)])

         call check(error, maxval(abs(v - h_ref)) < 1.0e-12_wp, &
                    "h_layer bit-exact vs the uninterrupted reference -- "// &
                    "no mid-run clobber from the narrow wind/bathymetry setters")
      end block runblock
      status = rdb_ocean_destroy(handle)
      call check(error, status == OK, "cleanup destroy (run)")
   end subroutine test_setter_no_clobber

   subroutine test_host_ptr_stable_across_step(error)
      !! The fact D3.4's silent-re-sync design rests on ("This must be
      !! asserted in a P2 test"): the Fortran host allocation for a
      !! P2-exposed array is made ONCE (at `ocean_state_init`/
      !! `enter_data`) and `!$acc update self` writes back INTO THE SAME
      !! BUFFER -- so a pointer/numpy view taken before a step stays valid
      !! across a later refresh, and the getter never needs re-issuing a
      !! new address. If this test fails, the design is wrong: a Python
      !! `Field`'s cached buffer would need to re-fetch the pointer on
      !! every refresh, not just re-check `generation`.
      !!
      !! CAVEAT: on this host-only build `!$acc update self` is an inert
      !! comment, so this only proves Fortran's own `allocatable` does not
      !! reallocate `h_layer` across `step()`+`refresh_host()` calls --
      !! it says NOTHING about the GPU `mem:separate` path, which is
      !! where a stray re-`allocate` in a kernel workspace helper would
      !! actually break the contract. Must also be run on the NVHPC GPU
      !! build before this design is trusted end-to-end (CLAUDE.md).
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle, ptr1, ptr2
      character(len=:), allocatable :: nml
      integer(c_int) :: status, nx1, ny1, nz1, gen1, nx2, ny2, nz2, gen2

      nml = ocean_nml()
      handle = c_null_ptr
      status = rdb_ocean_create_from_string(nml, len(nml, kind=c_int), handle)
      checks: block
         call check(error, status == OK, "create succeeds")
         if (allocated(error)) exit checks

         status = rdb_ocean_get_h_layer_ptr(handle, ptr1, nx1, ny1, nz1, gen1)
         call check(error, status == OK, "get_h_layer_ptr (pre-step) succeeds")
         if (allocated(error)) exit checks
         call check(error, c_associated(ptr1), "pre-step pointer is non-null")
         if (allocated(error)) exit checks

         status = rdb_ocean_step(handle, 3_c_int)
         call check(error, status == OK, "step(3) succeeds")
         if (allocated(error)) exit checks
         status = rdb_ocean_refresh_host(handle)
         call check(error, status == OK, "refresh_host succeeds")
         if (allocated(error)) exit checks

         status = rdb_ocean_get_h_layer_ptr(handle, ptr2, nx2, ny2, nz2, gen2)
         call check(error, status == OK, "get_h_layer_ptr (post-step) succeeds")
         if (allocated(error)) exit checks

         call check(error, gen2 > gen1, "generation (outer_step_count) advances after step+refresh")
         if (allocated(error)) exit checks
         call check(error, nx1 == nx2 .and. ny1 == ny2 .and. nz1 == nz2, &
                    "extents unchanged across step+refresh")
         if (allocated(error)) exit checks
         call check(error, c_associated(ptr1, ptr2), &
                    "h_layer's host allocation address is STABLE across "// &
                    "step()+refresh_host() -- the fact the D3.4 lazy-sync "// &
                    "design rests on")
      end block checks
      status = rdb_ocean_destroy(handle)
      call check(error, status == OK, "cleanup destroy")
   end subroutine test_host_ptr_stable_across_step

end module test_ocean_api_p2
