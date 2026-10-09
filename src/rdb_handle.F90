!! Opaque handle type for the ocean C ABI (Python runtime API plan, P1).
!!
!! Transplanted from the pre-coastal-split `feat/api-ocean-write` branch
!! (`tmp_local_artifacts/python_ffi_scope/recovered/fortran_api_ocean_write/rdb_handle.F90`),
!! stripped of the unstructured/coastal twin — this tree is ocean-only, so
!! there is exactly one handle kind and one magic sentinel.
!!
!! This module is deliberately mechanical: heap-allocate, sentinel-tag,
!! resolve, idempotent-destroy. It knows nothing about ocean policy (the
!! single-live-handle rule lives in `rdb_ocean_api`, which is the only
!! caller) — "the handle table stays generic" (docs/ocean_python_api_plan.md).
module rdb_handle
   use, intrinsic :: iso_c_binding, only: c_associated, c_f_pointer, c_loc, &
                                                                             c_null_ptr, c_ptr
   use rdb_constants, only: wp
   use rdb_config, only: config_t
   use rdb_grid, only: hgrid_t
   use rdb_ocean_state, only: ocean_state_t
   use rdb_ocean_engine, only: ocean_engine_t
   implicit none
   private

   public :: ocean_handle_t
   public :: handle_create, handle_destroy, handle_check
   public :: HANDLE_OK, HANDLE_ERR_NULL, HANDLE_ERR_INVALID

   !! Magic-number sentinel. Arbitrary but unlikely-in-uninitialised-memory,
   !! and distinct from MAGIC_DESTROYED so a stale/garbage pointer and an
   !! already-destroyed handle are both rejected by `handle_check`.
   integer, parameter :: MAGIC_OCEAN = int(z'72616B4F')      !! 'rakO'
   integer, parameter :: MAGIC_DESTROYED = int(z'DEADBEEF')

   !! handle_check return codes.
   integer, parameter :: HANDLE_OK = 0
   integer, parameter :: HANDLE_ERR_NULL = 1
   integer, parameter :: HANDLE_ERR_INVALID = 2

   type :: ocean_handle_t
      !! One ocean simulation: config + the P2.4 setup/step/teardown
      !! engine (god-state + grid + geothermal + bc_source + decomp +
      !! sea-ice params — `rdb_ocean_engine`). `state`/`grid` below are
      !! thin views onto `engine%state`/`engine%grid` kept ONLY so the
      !! P2 accessor bodies (`h%state%...`, `h%grid%...`) did not all
      !! need touching when this handle was rewired onto the shared
      !! engine (P2.4): `state` is POINTER-associated to `engine%state`
      !! right after `engine_enter_data` (stable for the handle's
      !! lifetime — `engine` is a component of this heap-allocated,
      !! pointer-accessed handle, never copied/reallocated); `grid` is a
      !! ONE-TIME VALUE COPY of `engine%grid` right after `engine_setup`
      !! (safe because grid metadata is read-only after setup — no
      !! kernel ever mutates `nx_phys`/`dx`/... mid-run). The standalone
      !! `sf` slot from P1 is gone: `engine_setup` now fully configures
      !! `engine%state%surface_flux` from the namelist (P2.4 closes the
      !! "API surface flux is a separate, minimally-seeded object" gap),
      !! so `rdb_ocean_get/set_*_flux` read/write `h%state%surface_flux`
      !! directly, same as the driver.
      integer :: magic = 0
      type(config_t) :: cfg
      type(ocean_engine_t) :: engine
         !! Not declared `target` — component attributes may not include
         !! `target` (illegal syntax) — but this needs none: `h` is only
         !! ever reached through a `pointer` variable (`c_f_pointer` off
         !! the opaque `c_ptr` handle), and a pointer's entire pointee,
         !! subobjects included, is a valid pointer-association target
         !! for as long as that association persists. `h%state =>
         !! h%engine%state` below relies on exactly this.
      type(hgrid_t) :: grid
      type(ocean_state_t), pointer :: state => null()
      integer :: n_inner = 0
         !! Barotropic fast-loop substep count, mirrored from
         !! `engine%n_inner` once at create — kept as a separate field
         !! only because the C ABI's OWN stricter "n_inner must resolve
         !! to >= 1" contract (P1) is enforced here, not inside the
         !! shared `engine_setup` (which tolerates n_inner < 1 the same
         !! way `driver_run_ocean` does — see `ocean_engine_t`'s
         !! docstring).
      real(wp) :: t_current = 0.0_wp
      logical :: is_initialised = .false.
         !! True once create() has completed device mapping. `step`/
         !! `get_*` reject a handle that never got this far.
      logical :: is_pending = .false.
         !! P2.5: true from `rdb_ocean_create_pending` until
         !! `rdb_ocean_create_finalize` completes (success or failure —
         !! failure destroys the whole handle, per the F9 rollback
         !! contract). The `rdb_ocean_stage_*` geometry-injection calls
         !! require this; `is_initialised` (above) still gates every
         !! step/getter/setter exactly as before P2.5.
      logical :: device_mapped = .false.
         !! True between `ocean_state_enter_data` and `ocean_state_exit_data`
         !! — tells `destroy` whether an exit_data pairing is owed.
      logical :: host_is_current = .true.
         !! P2 lazy-sync gate (`docs/ocean_python_api_plan.md`,
         !! `06_python_surface_design.md` D3.4): true iff the HOST copies of
         !! the API-exposed prognostic/forcing arrays match the (possibly
         !! more advanced) device copies. `rdb_ocean_step` clears it
         !! unconditionally — it does not sync; `rdb_ocean_refresh_host`
         !! is the only thing that sets it back to true (via a leaf-array
         !! `!$acc update self`, never the aggregate `state`). Starts `true`:
         !! right after `create()` the host seed IS what `enter_data` copied
         !! to the device, so the two agree with nothing to refresh yet.
   end type ocean_handle_t

contains

   function handle_create() result(c_handle)
      !! Allocate a new ocean handle and return its opaque C pointer. Caller
      !! (`rdb_ocean_create_from_string`) is responsible for eventually
      !! passing it to `handle_destroy` — including on every early-return
      !! failure path, so a handle that never reached `is_initialised` does
      !! not leak.
      type(c_ptr) :: c_handle
      type(ocean_handle_t), pointer :: h

      allocate (h)
      h%magic = MAGIC_OCEAN
      c_handle = c_loc(h)
   end function handle_create

   subroutine handle_destroy(c_handle)
      !! Free an ocean handle. Idempotent on null/already-destroyed/garbage
      !! handles — silently no-ops rather than aborting, so a Python
      !! `__del__` can call this blind. Caller is responsible for having
      !! already unwound device residency and called `engine_teardown`
      !! (`engine%state%destroy()` releases the god-state's Fortran
      !! allocations) — `h%state` is only a pointer VIEW onto
      !! `h%engine%state` (see `ocean_handle_t`'s docstring), never
      !! separately allocated, so there is nothing to deallocate through
      !! it here; this only releases the handle's own heap block
      !! (which takes `engine` — a value component — down with it).
      type(c_ptr), intent(inout) :: c_handle
      type(ocean_handle_t), pointer :: h

      if (.not. c_associated(c_handle)) return
      call c_f_pointer(c_handle, h)
      if (h%magic /= MAGIC_OCEAN) return   !! garbage or already destroyed
      h%magic = MAGIC_DESTROYED
      h%state => null()
      deallocate (h)
      c_handle = c_null_ptr
   end subroutine handle_destroy

   function handle_check(c_handle, h) result(status)
      !! Resolve an opaque `c_ptr` to a Fortran pointer and verify it.
      !!
      !! On success: `h` is associated and `status == HANDLE_OK`.
      !! On failure: `h` is null and `status` is one of `HANDLE_ERR_*`.
      !! Callers in `rdb_ocean_api` translate this into the ocean status
      !! code table (`OCEAN_STATUS_ERR_BAD_HANDLE`).
      type(c_ptr), intent(in) :: c_handle
      type(ocean_handle_t), pointer, intent(out) :: h
      integer :: status

      h => null()
      if (.not. c_associated(c_handle)) then
         status = HANDLE_ERR_NULL
         return
      end if
      call c_f_pointer(c_handle, h)
      if (h%magic /= MAGIC_OCEAN) then
         h => null()
         status = HANDLE_ERR_INVALID
         return
      end if
      status = HANDLE_OK
   end function handle_check

end module rdb_handle
