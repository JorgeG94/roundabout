!! Ocean C ABI — lifecycle (P1) + accessors and the error ring (P2).
!!
!! `bind(c)` entry points for create / step / destroy, introspection (time,
!! step count, grid shape, scalar diagnostics), the error ring, and the P2
!! state-array getters/setters. See the "P2 — accessors and the error ring"
!! section below the P1 lifecycle code for the D<->H contract these follow
!! (`docs/ocean_python_api_plan.md` S3, `06_python_surface_design.md` D3/D4).
!!
!! Every entry point returns `integer(c_int)` status — never `error stop`
!! (a library cannot abort its host process). Status codes reuse
!! `rdb_ocean_status` (config/setup/IC-seed failures reported by
!! `read_config_from_string` / `validate_config` / the `configure_ocean_*`
!! chain / `ocean_state_seed_from_cfg` via their optional `ierr`, P0) plus
!! three handle-specific codes added here.
!!
!! Setup mirrors `benchmarks/bench_ocean.F90` exactly — grid init ->
!! `ocean_state_t%init_from_config` -> vcoord params -> seed IC -> register
!! tracers -> surface-flux slot -> the `configure_ocean_*` chain -> resolve
!! `n_inner` -> `ocean_state_enter_data` + `sf` device mapping. That is the
!! proven init ordering (also used by `driver_run_ocean`); this module does
!! not re-derive it.
!!
!! Single live ocean handle, enforced. Several modules hold process-global
!! topology/device state keyed to "the one running ocean simulation"
!! (`rdb_ocean_halo`, `rdb_ice_ocean_coupler`, `rdb_ocean_data_input`,
!! `rdb_profiler`, `rdb_ocean_diag_fills`) — a second concurrent create would
!! silently corrupt the first rather than merely waste memory. `g_handle_live`
!! below is the guard; it is ocean-specific POLICY, which is why it lives
!! here and not in the generic `rdb_handle`.
!!
!! Never exposed: `comm_env_finalize` (irreversible, calls `MPI_Finalize`).
!! `comm_env_init` IS called, on every create -- see `build_pending_handle`.
!!
!! Failure contract (P0.1 review F9): `rdb_ocean_create_from_string`
!! calls `handle_destroy` on EVERY non-`OCEAN_STATUS_OK` branch below,
!! before returning the status code — never partially and never retries.
!! `handle_destroy` deallocates the WHOLE `ocean_handle_t` (`cfg` +
!! `state` pointer target together, `rdb_handle.F90`), so a
!! `configure_ocean_*` step that writes a knob (e.g. flips `enable`)
!! before validating it and failing partway through leaves nothing
!! observable: the entire handle carrying that half-configured `cfg`/
!! `state` is thrown away, not reused. This is why individual
!! `configure_ocean_*` procedures are not each required to be
!! transactional/rollback-safe on their own — the ONE rollback boundary
!! is "destroy the whole handle" at the `create()` level, and it must
!! stay that way: adding a retry-with-the-same-handle path would revive
!! the half-configured-state hazard the review flagged.
module rdb_ocean_api
   use, intrinsic :: iso_c_binding, only: c_char, c_double, c_int, c_loc, c_null_ptr, c_ptr
   use, intrinsic :: iso_fortran_env, only: output_unit
   use rdb_constants, only: wp
   use rdb_handle, only: ocean_handle_t, handle_create, handle_destroy, handle_check, &
                         HANDLE_OK
   use rdb_config, only: read_config_from_string, validate_config
   use rdb_ocean_engine, only: engine_setup, engine_enter_data, engine_step, &
                               engine_step_ice, engine_step_finalize, engine_exit_data, &
                               engine_teardown
   use rdb_ocean_dyn, only: ocean_dyn_flush_tracer_window
   use rdb_comm_env, only: comm_env_init
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_CONFIG_VALIDATE, &
                               OCEAN_STATUS_ERR_SETUP, OCEAN_STATUS_ERR_BAD_HANDLE, &
                               OCEAN_STATUS_ERR_ALREADY_EXISTS, &
                               OCEAN_STATUS_ERR_NOT_INITIALISED, &
                               OCEAN_STATUS_ERR_NOT_FOUND, OCEAN_STATUS_ERR_BAD_SHAPE, &
                               OCEAN_STATUS_ERR_NOT_PENDING
   use rdb_error_ring, only: fail, error_ring_push, error_ring_get, error_ring_count, &
                             error_ring_clear, ERROR_RING_MSG_LEN
   use rdb_ocean_eos_compute, only: ocean_eos_compute
   ! NetCDF-free (deliberately NOT `rdb_bathymetry`, which requires
   ! RDB_ENABLE_NETCDF=ON purely because it also contains the NetCDF
   ! reader — see rdb_ocean_bathymetry_inject's module docstring). Renamed
   ! on import so the existing rdb_ocean_set_bathymetry call sites below
   ! (P2) are untouched.
   use rdb_ocean_bathymetry_inject, only: fill_bathymetry_ghosts_array => bathymetry_fill_ghosts_array
   use rdb_ocean_periodic, only: ocean_periodic_wrap_centre_2d
   use rdb_ocean_fold_apply, only: ocean_fold_wrap_eta_2d
   ! ---- P2.5: pre-create geometry injection ----
   use rdb_ocean_halo_width, only: required_halo
   ! ---- P7: diagnostics discoverability + in-memory access ----
   use rdb_ocean_diag_derived, only: derived_catalog_size, derived_catalog_name
   use rdb_ocean_diag_fills, only: canonical_diag_catalog_size, canonical_diag_catalog_name
   implicit none
   private

   public :: rdb_ocean_create_from_string
   public :: rdb_ocean_step
   public :: rdb_ocean_destroy
   ! ---- P2.5: pre-create geometry injection ----
   public :: rdb_ocean_create_pending
   public :: rdb_ocean_create_finalize
   public :: rdb_ocean_stage_bathymetry
   public :: rdb_ocean_stage_metrics
   public :: rdb_ocean_stage_topology
   public :: rdb_ocean_required_halo
   public :: rdb_ocean_get_time
   public :: rdb_ocean_get_step_count
   public :: rdb_ocean_get_grid_info
   public :: rdb_ocean_get_total_mass
   public :: rdb_working_precision
   ! ---- P2: error ring ----
   public :: rdb_ocean_last_error
   public :: rdb_flush_logs
   ! ---- P2: D<->H refresh ----
   public :: rdb_ocean_refresh_host
   ! ---- P2: getters ----
   public :: rdb_ocean_get_h_layer_ptr
   public :: rdb_ocean_get_u_face_x_layer_ptr
   public :: rdb_ocean_get_v_face_y_layer_ptr
   public :: rdb_ocean_get_hu_ptr
   public :: rdb_ocean_get_hv_ptr
   public :: rdb_ocean_get_w_interface_ptr
   public :: rdb_ocean_get_rho_layer_ptr
   public :: rdb_ocean_get_b_ptr
   public :: rdb_ocean_get_bt_eta_ptr
   public :: rdb_ocean_get_tau_x_ptr
   public :: rdb_ocean_get_tau_y_ptr
   public :: rdb_ocean_get_q_heat_ptr
   public :: rdb_ocean_get_q_salt_ptr
   public :: rdb_ocean_get_kv_ptr
   public :: rdb_ocean_get_kt_ptr
   public :: rdb_ocean_get_ks_ptr
   public :: rdb_ocean_get_wet_t_ptr
   ! ---- P2: tracers by name ----
   public :: rdb_ocean_get_tracer_count
   public :: rdb_ocean_list_tracers
   public :: rdb_ocean_get_tracer_ptr
   ! ---- P2: narrow setters ----
   public :: rdb_ocean_set_tracer
   public :: rdb_ocean_set_h
   public :: rdb_ocean_set_u
   public :: rdb_ocean_set_v
   public :: rdb_ocean_set_bathymetry
   public :: rdb_ocean_set_wind
   public :: rdb_ocean_set_heat_flux
   public :: rdb_ocean_set_salt_flux
   ! ---- P2: scalar diagnostics ----
   public :: rdb_ocean_get_kinetic_energy
   ! ---- P7: diagnostics discoverability ----
   public :: rdb_ocean_derived_catalog_size
   public :: rdb_ocean_derived_catalog_name
   public :: rdb_ocean_canonical_catalog_size
   public :: rdb_ocean_canonical_catalog_name
   public :: rdb_ocean_get_diag_count
   public :: rdb_ocean_list_diags
   ! ---- P7: in-memory diagnostic access ----
   public :: rdb_ocean_get_diagnostic_ptr

   logical, save :: g_handle_live = .false.
      !! True from a successful create() until destroy(). Guards the single-
      !! live-ocean-handle invariant (see module header). Module-`save`
      !! state is the correct tool here — it mirrors the process-global
      !! singletons it is protecting against, not a design smell.

   real(wp), parameter :: RHO0_DIAG = 1025.0_wp
      !! Reference density used ONLY to give the total-mass diagnostic
      !! (`rdb_ocean_get_total_mass`) units of kg. The solver itself is
      !! Boussinesq (continuity conserves volume/thickness, not mass), so
      !! this constant does not feed back into the dynamics anywhere — it
      !! exists purely so the diagnostic is dimensionally a mass rather than
      !! a bare volume, and cancels out of any before/after conservation
      !! check a caller does with it.

contains

   ! ================================================================
   ! Lifecycle
   ! ================================================================

   function rdb_ocean_create_from_string(nml_text, nml_len, handle_out) &
      result(status) bind(c, name="rdb_ocean_create_from_string")
      !! Build a config from an in-memory namelist string (no filesystem
      !! touch), set up the ocean dyn-core exactly as `bench_ocean` /
      !! `driver_run_ocean` do, map it onto the device, and hand back an
      !! opaque handle. On ANY failure, `handle_out` is `c_null_ptr` and the
      !! partially-built handle (if one was allocated) is freed — never a
      !! half-initialised handle escaping to the caller.
      integer(c_int), intent(in), value :: nml_len
      character(kind=c_char), intent(in) :: nml_text(nml_len)
      type(c_ptr), intent(out) :: handle_out
      integer(c_int) :: status

      type(c_ptr) :: c_handle
      type(ocean_handle_t), pointer :: h
      character(len=:), allocatable :: text

      handle_out = c_null_ptr
      call error_ring_clear()

      if (g_handle_live) then
         call error_ring_push("rdb_ocean_create_from_string: a live ocean handle "// &
                              "already exists in this process (multi-instance is not "// &
                              "supported) — destroy it before creating another")
         status = int(OCEAN_STATUS_ERR_ALREADY_EXISTS, c_int)
         return
      end if

      call c_to_f_string(nml_text, nml_len, text)

      call build_pending_handle(text, c_handle, h, status)
      if (status /= int(OCEAN_STATUS_OK, c_int)) return   ! already destroyed
      g_handle_live = .true.

      ! ---- Setup (P2.4) — driver_run_ocean's exact 21-stage sequence,
      ! via the shared engine. Single-rank: compute_rank=0, compute_size=1
      ! (defaults); no restart on the C ABI yet (a future phase's concern,
      ! not silently dropped — `engine_setup`'s `restart_file` is simply
      ! not passed here). This is where the API used to run its OWN
      ! independent 9-stage subset (configure_ocean_metrics/_forcing/
      ! _drag/_vmix/_lateral/_pgf/_bt/_bt_split/_land_mask only) — every
      ! other stage (bc/diag/hdiff/porous/p_surf/sponge/tides/tracers/
      ! wave_drag/wetdry/dataovr/halo_init) is now reachable through here
      ! too. No geometry staged on this single-call path (P2.5's
      ! rdb_ocean_create_pending + rdb_ocean_stage_* do that) — so
      ! engine_setup falls through to its ordinary namelist-driven path,
      ! byte-identical to before P2.5. ----
      call complete_ocean_create(c_handle, h, status)
      handle_out = c_handle   ! c_null_ptr iff complete_ocean_create failed
   end function rdb_ocean_create_from_string

   subroutine build_pending_handle(text, c_handle, h, status)
      !! Shared prefix of `rdb_ocean_create_from_string` AND
      !! `rdb_ocean_create_pending` (P2.5): allocate a handle, parse +
      !! validate the config, check `dt_fixed`. Does NOT touch
      !! `g_handle_live` — callers set it themselves once they know which
      !! of the two flows they are in (immediate `complete_ocean_create` vs
      !! staying pending for geometry injection).
      !!
      !! On success: `c_handle`/`h` are a valid, freshly allocated handle
      !! (`h%is_pending`/`h%is_initialised` both still `.false.`).
      !! On failure: the handle has ALREADY been destroyed (matches the F9
      !! "destroy the whole handle on any failure" contract) — `c_handle`
      !! is left as returned by `handle_create`/`handle_destroy` (the
      !! caller must not touch it further) and `h` is null.
      character(len=*), intent(in) :: text
      type(c_ptr), intent(out) :: c_handle
      type(ocean_handle_t), pointer, intent(out) :: h
      integer(c_int), intent(out) :: status

      integer :: ierr, hstat

      ! Bring the comm env up before ANYTHING touches it.  The API is
      ! single-rank by design (compute_rank=0, compute_size=1), but
      ! "single-rank" is not "no MPI": in an RDB_ENABLE_MPI=ON build the
      ! setup path reaches `comm_world()`, and reaching it before MPI is
      ! initialised is not a soft failure --
      !
      !   *** The MPI_Comm_f2c() function was called before MPI_INIT was
      !   *** invoked.  This is disallowed by the MPI standard.
      !
      ! -- the process aborts.  Both create entry points funnel through
      ! here, so this is the one place that needs it, and `comm_env_init`
      ! is idempotent, so paying it per create costs nothing.
      !
      ! NOT a fix for being embedded in a host that has ALREADY called
      ! MPI_Init (mpi4py, a C driver of its own).  `pic_mpi_init` calls
      ! `MPI_Init_thread` unconditionally and pic-mpi exposes no
      ! `MPI_Initialized` wrapper to guard on, so that case would still
      ! double-initialise.  Closing it needs an upstream wrapper: the
      ! `no-mpi-in-rdb` rule means this file cannot ask MPI directly.
      call comm_env_init()

      c_handle = handle_create()
      hstat = handle_check(c_handle, h)
      if (hstat /= HANDLE_OK) then
         ! Unreachable in practice (we just created it), but never trust a
         ! resolve blindly — fail loud rather than dereference a null h.
         call handle_destroy(c_handle)
         h => null()
         status = int(OCEAN_STATUS_ERR_BAD_HANDLE, c_int)
         return
      end if

      call read_config_from_string(text, h%cfg, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) then
         status = int(ierr, c_int)
         call handle_destroy(c_handle)
         h => null()
         return
      end if

      call validate_config(h%cfg, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) then
         status = int(ierr, c_int)
         call handle_destroy(c_handle)
         h => null()
         return
      end if

      if (h%cfg%dt_fixed <= 0.0_wp) then
         status = int(OCEAN_STATUS_ERR_CONFIG_VALIDATE, c_int)
         call handle_destroy(c_handle)
         h => null()
         return
      end if

      status = int(OCEAN_STATUS_OK, c_int)
   end subroutine build_pending_handle

   subroutine complete_ocean_create(c_handle, h, status)
      !! Shared tail of `rdb_ocean_create_from_string` AND
      !! `rdb_ocean_create_finalize` (P2.5): `engine_setup` ->
      !! resolved-`n_inner` check -> `engine_enter_data` -> finalize handle
      !! bookkeeping. Consumes whatever geometry `h%engine%staged_*` fields
      !! carry (empty/unset for the single-call path — byte-identical to
      !! before P2.5). On any failure, destroys the WHOLE handle (F9) and
      !! clears the single-live-handle guard; `c_handle` is `c_null_ptr` and
      !! `h` is null on return in that case, matching `rdb_ocean_destroy`'s
      !! out-null convention.
      type(c_ptr), intent(inout) :: c_handle
      type(ocean_handle_t), pointer, intent(inout) :: h
      integer(c_int), intent(out) :: status

      integer :: ierr

      call engine_setup(h%engine, h%cfg, ierr)
      if (ierr /= OCEAN_STATUS_OK) then
         status = int(ierr, c_int)
         call handle_destroy(c_handle)
         g_handle_live = .false.
         h => null()
         return
      end if

      if (h%engine%n_inner < 1) then
         call error_ring_push("rdb_ocean_create: barotropic substep count "// &
                              "n_inner is unresolved (< 1) — set &ocean_bt_nml auto_n_inner "// &
                              "= .true. or n_inner explicitly")
         status = int(OCEAN_STATUS_ERR_SETUP, c_int)
         call handle_destroy(c_handle)
         g_handle_live = .false.
         h => null()
         return
      end if
      h%n_inner = h%engine%n_inner

      ! ---- Device placement ----
      call engine_enter_data(h%engine, h%cfg)
      ! `h%grid`/`h%state` are views onto the engine (see `ocean_handle_t`'s
      ! docstring) so every existing P2 accessor body below (`h%grid%...`,
      ! `h%state%...`) keeps working unchanged.
      h%grid = h%engine%grid
      h%state => h%engine%state
      h%device_mapped = .true.

      h%is_initialised = .true.
      h%is_pending = .false.
      status = int(OCEAN_STATUS_OK, c_int)
   end subroutine complete_ocean_create

   function rdb_ocean_create_pending(nml_text, nml_len, handle_out) &
      result(status) bind(c, name="rdb_ocean_create_pending")
      !! P2.5 phase 1 of 2: build + validate a config and allocate a
      !! handle EXACTLY like `rdb_ocean_create_from_string`, but stop
      !! there — does NOT run `engine_setup` or map the device. Claims the
      !! single-live-handle guard immediately (a second concurrent
      !! create/create_pending is refused from this point on, even though
      !! `engine_setup` has not yet touched any process-global state — the
      !! invariant is "one handle mid-create at a time", not merely "one
      !! finished one").
      !!
      !! Between this call and `rdb_ocean_create_finalize`, the
      !! `rdb_ocean_stage_*` entry points below inject geometry
      !! (bathymetry, grid metrics, topology) that `engine_setup` consumes
      !! — exactly the window passive-tracer registration already needed
      !! (`registry_locked` closes at `enter_data`, `rdb_multilayer_state.F90`).
      !! A caller with no geometry to inject has no reason to use this
      !! entry point — `rdb_ocean_create_from_string` remains the
      !! normal one-call path and is completely unaffected by this phase.
      integer(c_int), intent(in), value :: nml_len
      character(kind=c_char), intent(in) :: nml_text(nml_len)
      type(c_ptr), intent(out) :: handle_out
      integer(c_int) :: status

      type(c_ptr) :: c_handle
      type(ocean_handle_t), pointer :: h
      character(len=:), allocatable :: text

      handle_out = c_null_ptr
      call error_ring_clear()

      if (g_handle_live) then
         call error_ring_push("rdb_ocean_create_pending: a live ocean handle "// &
                              "already exists in this process (multi-instance is not "// &
                              "supported) — destroy it before creating another")
         status = int(OCEAN_STATUS_ERR_ALREADY_EXISTS, c_int)
         return
      end if

      call c_to_f_string(nml_text, nml_len, text)

      call build_pending_handle(text, c_handle, h, status)
      if (status /= int(OCEAN_STATUS_OK, c_int)) return   ! already destroyed

      h%is_pending = .true.
      g_handle_live = .true.
      handle_out = c_handle
   end function rdb_ocean_create_pending

   function rdb_ocean_create_finalize(c_handle) result(status) &
      bind(c, name="rdb_ocean_create_finalize")
      !! P2.5 phase 2 of 2: complete a handle started by
      !! `rdb_ocean_create_pending` (optionally staged with
      !! `rdb_ocean_stage_*` geometry in between) — runs `engine_setup`
      !! (consuming any staged geometry) through device mapping, exactly
      !! like the tail of `rdb_ocean_create_from_string`.
      !!
      !! On success, `c_handle` is unchanged (same pointer value, now fully
      !! initialised) — call the ordinary P2 accessors/`rdb_ocean_step`
      !! on it from here. On failure, the WHOLE handle is destroyed (F9)
      !! and `c_handle` is set to `c_null_ptr` — same out-null convention
      !! as `rdb_ocean_destroy`, since (unlike `create_from_string`,
      !! whose caller only ever sees a handle on success) a caller here
      !! already holds a handle value going in.
      type(c_ptr), intent(inout) :: c_handle
      integer(c_int) :: status

      type(ocean_handle_t), pointer :: h
      integer :: hstat

      call error_ring_clear()
      hstat = handle_check(c_handle, h)
      if (hstat /= HANDLE_OK) then
         status = int(OCEAN_STATUS_ERR_BAD_HANDLE, c_int)
         return
      end if
      if (.not. h%is_pending) then
         call error_ring_push("rdb_ocean_create_finalize: handle is not in the "// &
                              "pending (staging) window — call rdb_ocean_create_pending "// &
                              "first (this handle may already have been finalised)")
         status = int(OCEAN_STATUS_ERR_NOT_PENDING, c_int)
         return
      end if

      call complete_ocean_create(c_handle, h, status)
   end function rdb_ocean_create_finalize

   function rdb_ocean_stage_bathymetry(c_handle, b_data, nx_p, ny_p, convention) &
      result(status) bind(c, name="rdb_ocean_stage_bathymetry")
      !! P2.5: stage an interior-sized `(nx_p, ny_p)` bathymetry array on
      !! the PENDING handle `c_handle` (`rdb_ocean_create_pending`).
      !! Consumed by `engine_setup` inside `rdb_ocean_create_finalize` —
      !! see `rdb_ocean_bathymetry_inject` for the sign-normalisation +
      !! wet-fraction validation the array goes through THERE (not here:
      !! shape can't be checked against `nx_phys`/`ny_phys` until the grid
      !! exists, which `engine_setup` builds).
      type(c_ptr), intent(in), value :: c_handle
      integer(c_int), intent(in), value :: nx_p, ny_p
      real(c_double), intent(in) :: b_data(nx_p, ny_p)
      integer(c_int), intent(in), value :: convention
         !! `BATHY_CONVENTION_DEPTH_POSITIVE_DOWN` (1) or
         !! `_HEIGHT_POSITIVE_UP` (2) — REQUIRED, no default (D6.2: sign is
         !! the single most dangerous argument in this API). Any other
         !! value is rejected at `create_finalize` time (not here).
      integer(c_int) :: status

      type(ocean_handle_t), pointer :: h

      status = resolve_ocean_pending(c_handle, h)
      if (status /= OCEAN_STATUS_OK) return

      if (allocated(h%engine%staged_bathymetry)) deallocate (h%engine%staged_bathymetry)
      allocate (h%engine%staged_bathymetry(nx_p, ny_p), source=real(b_data, wp))
      h%engine%staged_bathymetry_convention = int(convention)
      status = int(OCEAN_STATUS_OK, c_int)
   end function rdb_ocean_stage_bathymetry

   function rdb_ocean_stage_metrics(c_handle, x, y, dx, dy, area, nxp, nyp, nx, ny) &
      result(status) bind(c, name="rdb_ocean_stage_metrics")
      !! P2.5: stage in-memory MOM6-style supergrid arrays on the PENDING
      !! handle `c_handle` — `metrics_assemble_from_supergrid_arrays`'s
      !! exact layout (`x`/`y`: `(nxp,nyp)` degrees; `dx`: `(nx,nyp)` m;
      !! `dy`: `(nxp,ny)` m; `area`: `(nx,ny)` m^2, where
      !! `nxp=2*nx_phys+1`, `nyp=2*ny_phys+1`, `nx=2*nx_phys`,
      !! `ny=2*ny_phys`). Consumed inside `rdb_ocean_create_finalize`,
      !! bypassing `cfg%ocean%grid%grid_config` entirely — exports the
      !! SAME assembler the NetCDF supergrid reader and the analytic
      !! tripolar generator already use, so this and a mosaic-file grid
      !! produce identical metrics for the identical arrays. Shape is
      !! validated at `create_finalize` time (against the grid, which does
      !! not exist yet here).
      type(c_ptr), intent(in), value :: c_handle
      integer(c_int), intent(in), value :: nxp, nyp, nx, ny
      real(c_double), intent(in) :: x(nxp, nyp), y(nxp, nyp)
      real(c_double), intent(in) :: dx(nx, nyp), dy(nxp, ny)
      real(c_double), intent(in) :: area(nx, ny)
      integer(c_int) :: status

      type(ocean_handle_t), pointer :: h

      status = resolve_ocean_pending(c_handle, h)
      if (status /= OCEAN_STATUS_OK) return

      if (allocated(h%engine%staged_metrics_x)) then
         deallocate (h%engine%staged_metrics_x, h%engine%staged_metrics_y, &
                     h%engine%staged_metrics_dx, h%engine%staged_metrics_dy, &
                     h%engine%staged_metrics_area)
      end if
      allocate (h%engine%staged_metrics_x(nxp, nyp), source=real(x, wp))
      allocate (h%engine%staged_metrics_y(nxp, nyp), source=real(y, wp))
      allocate (h%engine%staged_metrics_dx(nx, nyp), source=real(dx, wp))
      allocate (h%engine%staged_metrics_dy(nxp, ny), source=real(dy, wp))
      allocate (h%engine%staged_metrics_area(nx, ny), source=real(area, wp))
      status = int(OCEAN_STATUS_OK, c_int)
   end function rdb_ocean_stage_metrics

   function rdb_ocean_stage_topology(c_handle, periodic_x, periodic_y) &
      result(status) bind(c, name="rdb_ocean_stage_topology")
      !! P2.5: stage Oceananigans-style grid topology (per-dimension
      !! periodicity — "the grid owns periodicity", not the per-edge
      !! `&ocean_bc_nml` tags) on the PENDING handle `c_handle`. Consumed
      !! inside `rdb_ocean_create_finalize` via
      !! `ocean_bc_state_set_topology`, run AFTER the namelist edge tags
      !! are parsed and OVERRIDING whatever they derived for
      !! `periodic_x`/`periodic_y` (an axis NOT marked periodic here keeps
      !! whatever physical BC the namelist gave it).
      type(c_ptr), intent(in), value :: c_handle
      integer(c_int), intent(in), value :: periodic_x, periodic_y
         !! Nonzero = periodic on that axis.
      integer(c_int) :: status

      type(ocean_handle_t), pointer :: h

      status = resolve_ocean_pending(c_handle, h)
      if (status /= OCEAN_STATUS_OK) return

      h%engine%staged_periodic_x = (periodic_x /= 0_c_int)
      h%engine%staged_periodic_y = (periodic_y /= 0_c_int)
      h%engine%has_staged_topology = .true.
      status = int(OCEAN_STATUS_OK, c_int)
   end function rdb_ocean_stage_topology

   function rdb_ocean_required_halo(pv_adv_scheme, pv_adv_scheme_len, &
                                    tracer_recon, tracer_recon_len, &
                                    periodic, tripolar_fold, decomposed, &
                                    kappa_shear_at_vertex) result(ng) &
      bind(c, name="rdb_ocean_required_halo")
      !! Stateless (no handle required — like `rdb_working_precision`):
      !! the minimum `nghost` for a given scheme/topology selection,
      !! folding the six scattered per-scheme/per-topology minimums behind
      !! `rdb_ocean_halo_width::required_halo` so a Python caller sizing
      !! a grid's halo before `create()` agrees with what `engine_setup`
      !! will itself enforce. A zero-length scheme string means "unset"
      !! (baseline 2), matching the Fortran function's
      !! absent-optional-argument behaviour.
      integer(c_int), intent(in), value :: pv_adv_scheme_len
      character(kind=c_char), intent(in) :: pv_adv_scheme(pv_adv_scheme_len)
      integer(c_int), intent(in), value :: tracer_recon_len
      character(kind=c_char), intent(in) :: tracer_recon(tracer_recon_len)
      integer(c_int), intent(in), value :: periodic, tripolar_fold, decomposed
      integer(c_int), intent(in), value :: kappa_shear_at_vertex
      integer(c_int) :: ng

      character(len=:), allocatable :: pv_s, tr_s

      call c_to_f_string(pv_adv_scheme, pv_adv_scheme_len, pv_s)
      call c_to_f_string(tracer_recon, tracer_recon_len, tr_s)
      ng = int(required_halo(pv_adv_scheme=pv_s, tracer_recon=tr_s, &
                             periodic=(periodic /= 0_c_int), &
                             tripolar_fold=(tripolar_fold /= 0_c_int), &
                             decomposed=(decomposed /= 0_c_int), &
                             kappa_shear_at_vertex=(kappa_shear_at_vertex /= 0_c_int)), c_int)
   end function rdb_ocean_required_halo

   function rdb_ocean_step(c_handle, n_steps) result(status) &
      bind(c, name="rdb_ocean_step")
      !! Advance `n_steps` fixed-`dt` (`cfg%dt_fixed`) outer steps via the
      !! shared `engine_step` / `engine_step_ice` / `engine_step_finalize`
      !! sequence (P2.4 + P2.4b) — the SAME calls `driver_run_ocean`'s time
      !! loop makes: the dyn-core advance, then sea-ice per-step physics
      !! (`engine_step_ice` — a no-op when `&ocean_ice_nml enable = .false.`,
      !! so this is bit-identical to before P2.4b for every non-ice config),
      !! then surface-flux-component assembly / the diag step. `n_steps <=
      !! 0` is a successful no-op (mirrors an empty range, not an error).
      type(c_ptr), intent(in), value :: c_handle
      integer(c_int), intent(in), value :: n_steps
      integer(c_int) :: status

      type(ocean_handle_t), pointer :: h
      integer :: i
      real(wp) :: dt
      integer :: step_ierr

      status = resolve_ocean(c_handle, h)
      if (status /= OCEAN_STATUS_OK) return
      if (n_steps <= 0_c_int) return

      dt = h%cfg%dt_fixed
      do i = 1, int(n_steps)
         call engine_step(h%engine, dt, h%t_current, ierr=step_ierr)
         if (step_ierr /= OCEAN_STATUS_OK) then
            status = int(step_ierr, c_int)
            return
         end if
         call engine_step_ice(h%engine, h%cfg, dt, h%t_current, ierr=step_ierr)
         if (step_ierr /= OCEAN_STATUS_OK) then
            status = int(step_ierr, c_int)
            return
         end if
         call engine_step_finalize(h%engine, dt, h%t_current, ierr=step_ierr)
         if (step_ierr /= OCEAN_STATUS_OK) then
            status = int(step_ierr, c_int)
            return
         end if
         h%t_current = h%t_current + dt
      end do
      ! P2 D<->H contract (docs/ocean_python_api_plan.md S3): a step only
      ! CLEARS the lazy-sync flag, never syncs — the next accessor that
      ! needs a live host view calls rdb_ocean_refresh_host itself.
      h%host_is_current = .false.
      status = int(OCEAN_STATUS_OK, c_int)
   end function rdb_ocean_step

   function rdb_ocean_destroy(c_handle) result(status) &
      bind(c, name="rdb_ocean_destroy")
      !! Idempotent: a null/already-destroyed/garbage handle is a no-op
      !! success (so a Python `__del__` can call this blind). On a live
      !! handle: unwind device residency (`engine_exit_data`) then release
      !! the god-state's host-side allocations + the process-global ocean-
      !! halo module state (`engine_teardown`), then free the handle itself
      !! and clear the single-live-handle guard.
      type(c_ptr), intent(inout) :: c_handle
      integer(c_int) :: status

      type(ocean_handle_t), pointer :: h
      integer :: hstat

      hstat = handle_check(c_handle, h)
      if (hstat /= HANDLE_OK) then
         ! Already destroyed, never created, or garbage: idempotent no-op.
         c_handle = c_null_ptr
         status = int(OCEAN_STATUS_OK, c_int)
         return
      end if

      if (h%device_mapped) then
         call engine_exit_data(h%engine)
         h%device_mapped = .false.
      end if
      ! P2.5: a handle that never got past rdb_ocean_create_pending (no
      ! finalize call, abandoned mid-staging) never ran engine_setup —
      ! engine%state was never init_from_config'd, ocean_halo_init/
      ! engine%geo%init never ran, so there is nothing for engine_teardown
      ! to tear down. Guard on engine%is_setup (set only at the tail of
      ! engine_setup) rather than calling it unconditionally.
      if (h%engine%is_setup) call engine_teardown(h%engine)
      h%state => null()

      call handle_destroy(c_handle)
      g_handle_live = .false.
      status = int(OCEAN_STATUS_OK, c_int)
   end function rdb_ocean_destroy

   ! ================================================================
   ! Queries
   ! ================================================================

   function rdb_ocean_get_time(c_handle, t_out) result(status) &
      bind(c, name="rdb_ocean_get_time")
      !! Current simulation time (seconds).
      type(c_ptr), intent(in), value :: c_handle
      real(c_double), intent(out) :: t_out
      integer(c_int) :: status

      type(ocean_handle_t), pointer :: h

      t_out = 0.0_c_double
      status = resolve_ocean(c_handle, h)
      if (status /= OCEAN_STATUS_OK) return
      t_out = real(h%t_current, c_double)
   end function rdb_ocean_get_time

   function rdb_ocean_get_step_count(c_handle, step_out) result(status) &
      bind(c, name="rdb_ocean_get_step_count")
      !! Outer-step counter. Reads `dyn%outer_step_count` directly (the
      !! kernel's own bookkeeping) rather than keeping a second counter on
      !! the handle, so there is exactly one source of truth.
      type(c_ptr), intent(in), value :: c_handle
      integer(c_int), intent(out) :: step_out
      integer(c_int) :: status

      type(ocean_handle_t), pointer :: h

      step_out = 0_c_int
      status = resolve_ocean(c_handle, h)
      if (status /= OCEAN_STATUS_OK) return
      step_out = int(h%state%dyn%outer_step_count, c_int)
   end function rdb_ocean_get_step_count

   function rdb_ocean_get_grid_info(c_handle, nx, ny, nz, nghost) result(status) &
      bind(c, name="rdb_ocean_get_grid_info")
      !! Physical (interior, ghost-excluded) grid shape + ghost width.
      type(c_ptr), intent(in), value :: c_handle
      integer(c_int), intent(out) :: nx, ny, nz, nghost
      integer(c_int) :: status

      type(ocean_handle_t), pointer :: h

      nx = 0
      ny = 0
      nz = 0
      nghost = 0
      status = resolve_ocean(c_handle, h)
      if (status /= OCEAN_STATUS_OK) return
      nx = int(h%grid%nx_phys, c_int)
      ny = int(h%grid%ny_phys, c_int)
      nz = int(h%state%multilayer%nz_ml, c_int)
      nghost = int(h%grid%nghost, c_int)
   end function rdb_ocean_get_grid_info

   function rdb_ocean_get_total_mass(c_handle, m_out) result(status) &
      bind(c, name="rdb_ocean_get_total_mass")
      !! One scalar diagnostic — total water mass over the physical domain
      !! (`sum(h_layer) * dx * dy * RHO0_DIAG`) — so a caller can prove the
      !! solver actually advanced (and, in a closed quiescent/wall basin,
      !! that it is conserving mass) without any state-array accessor (P2).
      !! `!$acc update self` on the leaf array via `associate` (never the
      !! aggregate `ocean_state_t`/`multilayer_state_t`) before summing, per
      !! the D<->H contract — inert on a host build, load-bearing on GPU.
      type(c_ptr), intent(in), value :: c_handle
      real(c_double), intent(out) :: m_out
      integer(c_int) :: status

      type(ocean_handle_t), pointer :: h
      real(wp) :: total

      m_out = 0.0_c_double
      status = resolve_ocean(c_handle, h)
      if (status /= OCEAN_STATUS_OK) return

      associate (hlayer => h%state%multilayer%h_layer, ng => h%grid%nghost, &
                 nxp => h%grid%nx_phys, nyp => h%grid%ny_phys)
         !$acc update self(hlayer)
         total = sum(hlayer(ng + 1:ng + nxp, ng + 1:ng + nyp, :))
      end associate

      m_out = real(total*h%grid%dx*h%grid%dy*RHO0_DIAG, c_double)
   end function rdb_ocean_get_total_mass

   function rdb_working_precision() result(bytes) &
      bind(c, name="rdb_working_precision")
      !! Bytes per `wp` (4 or 8) — so a Python caller resolves float32 vs
      !! float64 at load time instead of guessing. No handle needed: this
      !! is a build-time constant.
      integer(c_int) :: bytes

      bytes = int(storage_size(1.0_wp)/8, c_int)
   end function rdb_working_precision

   ! ================================================================
   ! P2 — the error ring
   ! ================================================================

   function rdb_ocean_last_error(idx, buf, cap) result(len_out) &
      bind(c, name="rdb_ocean_last_error")
      !! Read the error ring at ring-relative index `idx` (0 = most
      !! recent push). Copies up to `cap` bytes of the trimmed message into
      !! `buf`, NUL-terminating if room remains, and returns the FULL
      !! trimmed message length — a `snprintf`-style contract: a returned
      !! length >= `cap` means the copy was truncated. Returns 0 (buf
      !! untouched) if `idx` is out of `[0, count)` or `cap <= 0`.
      !!
      !! The ring is the ONLY channel a caller should trust for the
      !! SPECIFIC reason behind a coarse status code (D4.3,
      !! `06_python_surface_design.md`) — never scrape stdout.
      integer(c_int), intent(in), value :: idx
      integer(c_int), intent(in), value :: cap
      character(kind=c_char), intent(out) :: buf(cap)
      integer(c_int) :: len_out

      character(len=ERROR_RING_MSG_LEN) :: msg
      integer :: n, ncopy, i

      len_out = 0_c_int
      if (idx < 0_c_int .or. idx >= int(error_ring_count(), c_int) .or. cap <= 0_c_int) return

      msg = error_ring_get(int(idx))
      n = len_trim(msg)
      len_out = int(n, c_int)
      ncopy = min(n, int(cap))
      do i = 1, ncopy
         buf(i) = msg(i:i)
      end do
      if (ncopy < cap) buf(ncopy + 1) = achar(0)
   end function rdb_ocean_last_error

   subroutine rdb_flush_logs() bind(c, name="rdb_flush_logs")
      !! Flush the buffered log stream (unit 6 / stdout — what
      !! `pic_logger`'s `global_logger` writes to; it exposes no flush of
      !! its own). Call this in a `finally` around a C entry point if the
      !! human-readable log needs to be ordered relative to Python's own
      !! stdout (D4.4) — the ring (above) is still the only channel to
      !! trust for the SPECIFIC failure reason.
      flush (output_unit)
   end subroutine rdb_flush_logs

   ! ================================================================
   ! P2 — D<->H refresh
   ! ================================================================

   function rdb_ocean_refresh_host(c_handle) result(status) &
      bind(c, name="rdb_ocean_refresh_host")
      !! Lazy device->host refresh of every P2-exposed leaf array, gated on
      !! `host_is_current` (a second call with nothing new on device is a
      !! cheap flag check, not a re-copy). `!$acc update self` runs on
      !! LEAF component names only, via `associate` — never the aggregate
      !! `ocean_state_t`/sub-state derived type (the documented
      !! `rdb_ocean_dyn.F90:2949` segfault: an aggregate D->H copy
      !! overwrites the host allocatable descriptors with DEVICE
      !! addresses). Inert on a host-only build (no `!$acc` support) —
      !! see CLAUDE.md's GPU-verification caveat.
      type(c_ptr), intent(in), value :: c_handle
      integer(c_int) :: status

      type(ocean_handle_t), pointer :: h

      status = resolve_ocean(c_handle, h)
      if (status /= OCEAN_STATUS_OK) return
      call ocean_handle_refresh_host(h)
      status = int(OCEAN_STATUS_OK, c_int)
   end function rdb_ocean_refresh_host

   ! ================================================================
   ! P2 — getters: (ptr, extents..., generation)
   !
   ! None of these refresh the host copy themselves — per the D3.4 lazy-
   ! sync design, the pointer + `generation` (== `dyn%outer_step_count` at
   ! issue) are cheap to hand back every call; the CALLER compares
   ! `generation` against what it last saw and, on a mismatch, calls
   ! `rdb_ocean_refresh_host` itself before dereferencing. This is what
   ! keeps "pay only when you touch it" true: a getter call alone never
   ! triggers a D->H copy.
   !
   ! `mass_flux_x/y_layer` and `flux_h_layer` are deliberately NOT exposed
   ! — they are `!$acc enter data create(...)`-mapped (rdb_multilayer_state
   ! .F90), so there is no valid host copy at any point in the run
   ! (docs/ocean_python_api_plan.md S3).
   ! ================================================================

   function rdb_ocean_get_h_layer_ptr(c_handle, ptr, nx, ny, nz, gen) result(status) &
      bind(c, name="rdb_ocean_get_h_layer_ptr")
      !! Layer thickness (m), cell-centred, FULL extent (ghosts included):
      !! `multilayer%h_layer`, shape `(nx_total, ny_total, nz_ml)`.
      type(c_ptr), intent(in), value :: c_handle
      type(c_ptr), intent(out) :: ptr
      integer(c_int), intent(out) :: nx, ny, nz, gen
      integer(c_int) :: status

      type(ocean_handle_t), pointer :: h
      real(wp), pointer :: tmp(:, :, :)

      ptr = c_null_ptr
      nx = 0
      ny = 0
      nz = 0
      gen = 0
      status = resolve_ocean(c_handle, h)
      if (status /= OCEAN_STATUS_OK) return
      tmp => h%state%multilayer%h_layer
      call fill_getter_3d(ptr, nx, ny, nz, gen, tmp, h%state%dyn%outer_step_count)
   end function rdb_ocean_get_h_layer_ptr

   function rdb_ocean_get_u_face_x_layer_ptr(c_handle, ptr, nx, ny, nz, gen) result(status) &
      bind(c, name="rdb_ocean_get_u_face_x_layer_ptr")
      !! West-face x-velocity (m/s): `multilayer%u_face_x_layer`, shape
      !! `(nx_total+1, ny_total, nz_ml)`.
      type(c_ptr), intent(in), value :: c_handle
      type(c_ptr), intent(out) :: ptr
      integer(c_int), intent(out) :: nx, ny, nz, gen
      integer(c_int) :: status

      type(ocean_handle_t), pointer :: h
      real(wp), pointer :: tmp(:, :, :)

      ptr = c_null_ptr
      nx = 0
      ny = 0
      nz = 0
      gen = 0
      status = resolve_ocean(c_handle, h)
      if (status /= OCEAN_STATUS_OK) return
      tmp => h%state%multilayer%u_face_x_layer
      call fill_getter_3d(ptr, nx, ny, nz, gen, tmp, h%state%dyn%outer_step_count)
   end function rdb_ocean_get_u_face_x_layer_ptr

   function rdb_ocean_get_v_face_y_layer_ptr(c_handle, ptr, nx, ny, nz, gen) result(status) &
      bind(c, name="rdb_ocean_get_v_face_y_layer_ptr")
      !! South-face y-velocity (m/s): `multilayer%v_face_y_layer`, shape
      !! `(nx_total, ny_total+1, nz_ml)`.
      type(c_ptr), intent(in), value :: c_handle
      type(c_ptr), intent(out) :: ptr
      integer(c_int), intent(out) :: nx, ny, nz, gen
      integer(c_int) :: status

      type(ocean_handle_t), pointer :: h
      real(wp), pointer :: tmp(:, :, :)

      ptr = c_null_ptr
      nx = 0
      ny = 0
      nz = 0
      gen = 0
      status = resolve_ocean(c_handle, h)
      if (status /= OCEAN_STATUS_OK) return
      tmp => h%state%multilayer%v_face_y_layer
      call fill_getter_3d(ptr, nx, ny, nz, gen, tmp, h%state%dyn%outer_step_count)
   end function rdb_ocean_get_v_face_y_layer_ptr

   function rdb_ocean_get_hu_ptr(c_handle, ptr, nx, ny, nz, gen) result(status) &
      bind(c, name="rdb_ocean_get_hu_ptr")
      !! West-face x transport (h*u, m^2/s): `multilayer%hu_face_x_layer`,
      !! same shape/stagger as `u_face_x_layer`.
      type(c_ptr), intent(in), value :: c_handle
      type(c_ptr), intent(out) :: ptr
      integer(c_int), intent(out) :: nx, ny, nz, gen
      integer(c_int) :: status

      type(ocean_handle_t), pointer :: h
      real(wp), pointer :: tmp(:, :, :)

      ptr = c_null_ptr
      nx = 0
      ny = 0
      nz = 0
      gen = 0
      status = resolve_ocean(c_handle, h)
      if (status /= OCEAN_STATUS_OK) return
      tmp => h%state%multilayer%hu_face_x_layer
      call fill_getter_3d(ptr, nx, ny, nz, gen, tmp, h%state%dyn%outer_step_count)
   end function rdb_ocean_get_hu_ptr

   function rdb_ocean_get_hv_ptr(c_handle, ptr, nx, ny, nz, gen) result(status) &
      bind(c, name="rdb_ocean_get_hv_ptr")
      !! South-face y transport (h*v, m^2/s): `multilayer%hv_face_y_layer`,
      !! same shape/stagger as `v_face_y_layer`.
      type(c_ptr), intent(in), value :: c_handle
      type(c_ptr), intent(out) :: ptr
      integer(c_int), intent(out) :: nx, ny, nz, gen
      integer(c_int) :: status

      type(ocean_handle_t), pointer :: h
      real(wp), pointer :: tmp(:, :, :)

      ptr = c_null_ptr
      nx = 0
      ny = 0
      nz = 0
      gen = 0
      status = resolve_ocean(c_handle, h)
      if (status /= OCEAN_STATUS_OK) return
      tmp => h%state%multilayer%hv_face_y_layer
      call fill_getter_3d(ptr, nx, ny, nz, gen, tmp, h%state%dyn%outer_step_count)
   end function rdb_ocean_get_hv_ptr

   function rdb_ocean_get_w_interface_ptr(c_handle, ptr, nx, ny, nz, gen) result(status) &
      bind(c, name="rdb_ocean_get_w_interface_ptr")
      !! Vertical velocity at layer interfaces (m/s):
      !! `multilayer%w_interface`, shape `(nx_total, ny_total, nz_ml+1)`,
      !! k=1 bed .. k=nz_ml+1 surface.
      type(c_ptr), intent(in), value :: c_handle
      type(c_ptr), intent(out) :: ptr
      integer(c_int), intent(out) :: nx, ny, nz, gen
      integer(c_int) :: status

      type(ocean_handle_t), pointer :: h
      real(wp), pointer :: tmp(:, :, :)

      ptr = c_null_ptr
      nx = 0
      ny = 0
      nz = 0
      gen = 0
      status = resolve_ocean(c_handle, h)
      if (status /= OCEAN_STATUS_OK) return
      tmp => h%state%multilayer%w_interface
      call fill_getter_3d(ptr, nx, ny, nz, gen, tmp, h%state%dyn%outer_step_count)
   end function rdb_ocean_get_w_interface_ptr

   function rdb_ocean_get_rho_layer_ptr(c_handle, ptr, nx, ny, nz, gen) result(status) &
      bind(c, name="rdb_ocean_get_rho_layer_ptr")
      !! In-situ density (kg/m^3): `multilayer%rho_layer`, same shape as
      !! `h_layer`. Filled by the EOS each thermo step (and by every P2
      !! setter that touches T/S/h).
      type(c_ptr), intent(in), value :: c_handle
      type(c_ptr), intent(out) :: ptr
      integer(c_int), intent(out) :: nx, ny, nz, gen
      integer(c_int) :: status

      type(ocean_handle_t), pointer :: h
      real(wp), pointer :: tmp(:, :, :)

      ptr = c_null_ptr
      nx = 0
      ny = 0
      nz = 0
      gen = 0
      status = resolve_ocean(c_handle, h)
      if (status /= OCEAN_STATUS_OK) return
      tmp => h%state%multilayer%rho_layer
      call fill_getter_3d(ptr, nx, ny, nz, gen, tmp, h%state%dyn%outer_step_count)
   end function rdb_ocean_get_rho_layer_ptr

   function rdb_ocean_get_b_ptr(c_handle, ptr, nx, ny, gen) result(status) &
      bind(c, name="rdb_ocean_get_b_ptr")
      !! Bathymetry (m, POSITIVE DOWN — `eta = sum(h) - b`):
      !! `barotropic%b`, shape `(nx_total, ny_total)`.
      type(c_ptr), intent(in), value :: c_handle
      type(c_ptr), intent(out) :: ptr
      integer(c_int), intent(out) :: nx, ny, gen
      integer(c_int) :: status

      type(ocean_handle_t), pointer :: h
      real(wp), pointer :: tmp(:, :)

      ptr = c_null_ptr
      nx = 0
      ny = 0
      gen = 0
      status = resolve_ocean(c_handle, h)
      if (status /= OCEAN_STATUS_OK) return
      tmp => h%state%barotropic%b
      call fill_getter_2d(ptr, nx, ny, gen, tmp, h%state%dyn%outer_step_count)
   end function rdb_ocean_get_b_ptr

   function rdb_ocean_get_bt_eta_ptr(c_handle, ptr, nx, ny, gen) result(status) &
      bind(c, name="rdb_ocean_get_bt_eta_ptr")
      !! Sea-surface height (m): `dyn%bt_work%bt_eta`, shape
      !! `(nx_total, ny_total)`. Diagnostic (`sum_k(h_layer) - bt_H_ref`),
      !! not independently settable — write `h` or `b` instead.
      type(c_ptr), intent(in), value :: c_handle
      type(c_ptr), intent(out) :: ptr
      integer(c_int), intent(out) :: nx, ny, gen
      integer(c_int) :: status

      type(ocean_handle_t), pointer :: h
      real(wp), pointer :: tmp(:, :)

      ptr = c_null_ptr
      nx = 0
      ny = 0
      gen = 0
      status = resolve_ocean(c_handle, h)
      if (status /= OCEAN_STATUS_OK) return
      tmp => h%state%dyn%bt_work%bt_eta
      call fill_getter_2d(ptr, nx, ny, gen, tmp, h%state%dyn%outer_step_count)
   end function rdb_ocean_get_bt_eta_ptr

   function rdb_ocean_get_tau_x_ptr(c_handle, ptr, nx, ny, gen) result(status) &
      bind(c, name="rdb_ocean_get_tau_x_ptr")
      !! East-face wind stress (N/m^2): `surface_stress%tau_x`, shape
      !! `(nx_total+1, ny_total)`.
      type(c_ptr), intent(in), value :: c_handle
      type(c_ptr), intent(out) :: ptr
      integer(c_int), intent(out) :: nx, ny, gen
      integer(c_int) :: status

      type(ocean_handle_t), pointer :: h
      real(wp), pointer :: tmp(:, :)

      ptr = c_null_ptr
      nx = 0
      ny = 0
      gen = 0
      status = resolve_ocean(c_handle, h)
      if (status /= OCEAN_STATUS_OK) return
      tmp => h%state%surface_stress%tau_x
      call fill_getter_2d(ptr, nx, ny, gen, tmp, h%state%dyn%outer_step_count)
   end function rdb_ocean_get_tau_x_ptr

   function rdb_ocean_get_tau_y_ptr(c_handle, ptr, nx, ny, gen) result(status) &
      bind(c, name="rdb_ocean_get_tau_y_ptr")
      !! North-face wind stress (N/m^2): `surface_stress%tau_y`, shape
      !! `(nx_total, ny_total+1)`.
      type(c_ptr), intent(in), value :: c_handle
      type(c_ptr), intent(out) :: ptr
      integer(c_int), intent(out) :: nx, ny, gen
      integer(c_int) :: status

      type(ocean_handle_t), pointer :: h
      real(wp), pointer :: tmp(:, :)

      ptr = c_null_ptr
      nx = 0
      ny = 0
      gen = 0
      status = resolve_ocean(c_handle, h)
      if (status /= OCEAN_STATUS_OK) return
      tmp => h%state%surface_stress%tau_y
      call fill_getter_2d(ptr, nx, ny, gen, tmp, h%state%dyn%outer_step_count)
   end function rdb_ocean_get_tau_y_ptr

   function rdb_ocean_get_q_heat_ptr(c_handle, ptr, nx, ny, gen) result(status) &
      bind(c, name="rdb_ocean_get_q_heat_ptr")
      !! Net surface heat flux (W/m^2): `h%state%surface_flux%Q_heat`. P2
      !! called this the "STANDALONE `sf` slot" (a separate, minimally-
      !! seeded object from `ocean_state%surface_flux`) because the P1
      !! step call never ran `ocean_surface_flux_assemble`. P2.4 unifies
      !! setup+step across all three callers via the shared engine
      !! (`rdb_ocean_engine`), so `rdb_ocean_step` now DOES run the
      !! assembler (via `engine_step_finalize`) against this SAME field —
      !! a write here is live-consumed the same way the driver's is.
      !! Shape `(nx_total, ny_total)`.
      type(c_ptr), intent(in), value :: c_handle
      type(c_ptr), intent(out) :: ptr
      integer(c_int), intent(out) :: nx, ny, gen
      integer(c_int) :: status

      type(ocean_handle_t), pointer :: h
      real(wp), pointer :: tmp(:, :)

      ptr = c_null_ptr
      nx = 0
      ny = 0
      gen = 0
      status = resolve_ocean(c_handle, h)
      if (status /= OCEAN_STATUS_OK) return
      tmp => h%state%surface_flux%Q_heat
      call fill_getter_2d(ptr, nx, ny, gen, tmp, h%state%dyn%outer_step_count)
   end function rdb_ocean_get_q_heat_ptr

   function rdb_ocean_get_q_salt_ptr(c_handle, ptr, nx, ny, gen) result(status) &
      bind(c, name="rdb_ocean_get_q_salt_ptr")
      !! Net surface salt flux: `h%state%surface_flux%Q_salt` — see
      !! `rdb_ocean_get_q_heat_ptr` for the P2.4 unification note.
      !! Shape `(nx_total, ny_total)`.
      type(c_ptr), intent(in), value :: c_handle
      type(c_ptr), intent(out) :: ptr
      integer(c_int), intent(out) :: nx, ny, gen
      integer(c_int) :: status

      type(ocean_handle_t), pointer :: h
      real(wp), pointer :: tmp(:, :)

      ptr = c_null_ptr
      nx = 0
      ny = 0
      gen = 0
      status = resolve_ocean(c_handle, h)
      if (status /= OCEAN_STATUS_OK) return
      tmp => h%state%surface_flux%Q_salt
      call fill_getter_2d(ptr, nx, ny, gen, tmp, h%state%dyn%outer_step_count)
   end function rdb_ocean_get_q_salt_ptr

   function rdb_ocean_get_kv_ptr(c_handle, ptr, nx, ny, nz, gen) result(status) &
      bind(c, name="rdb_ocean_get_kv_ptr")
      !! Vertical viscosity (m^2/s): `vmix%kv`, shape
      !! `(nx_total, ny_total, nz_ml+1)` (layer INTERFACES).
      type(c_ptr), intent(in), value :: c_handle
      type(c_ptr), intent(out) :: ptr
      integer(c_int), intent(out) :: nx, ny, nz, gen
      integer(c_int) :: status

      type(ocean_handle_t), pointer :: h
      real(wp), pointer :: tmp(:, :, :)

      ptr = c_null_ptr
      nx = 0
      ny = 0
      nz = 0
      gen = 0
      status = resolve_ocean(c_handle, h)
      if (status /= OCEAN_STATUS_OK) return
      tmp => h%state%vmix%kv
      call fill_getter_3d(ptr, nx, ny, nz, gen, tmp, h%state%dyn%outer_step_count)
   end function rdb_ocean_get_kv_ptr

   function rdb_ocean_get_kt_ptr(c_handle, ptr, nx, ny, nz, gen) result(status) &
      bind(c, name="rdb_ocean_get_kt_ptr")
      !! Vertical heat diffusivity (m^2/s): `vmix%kt`, shape
      !! `(nx_total, ny_total, nz_ml+1)` (layer INTERFACES).
      type(c_ptr), intent(in), value :: c_handle
      type(c_ptr), intent(out) :: ptr
      integer(c_int), intent(out) :: nx, ny, nz, gen
      integer(c_int) :: status

      type(ocean_handle_t), pointer :: h
      real(wp), pointer :: tmp(:, :, :)

      ptr = c_null_ptr
      nx = 0
      ny = 0
      nz = 0
      gen = 0
      status = resolve_ocean(c_handle, h)
      if (status /= OCEAN_STATUS_OK) return
      tmp => h%state%vmix%kt
      call fill_getter_3d(ptr, nx, ny, nz, gen, tmp, h%state%dyn%outer_step_count)
   end function rdb_ocean_get_kt_ptr

   function rdb_ocean_get_ks_ptr(c_handle, ptr, nx, ny, nz, gen) result(status) &
      bind(c, name="rdb_ocean_get_ks_ptr")
      !! Vertical salt (+ every passive tracer) diffusivity (m^2/s):
      !! `vmix%ks`, shape `(nx_total, ny_total, nz_ml+1)` (layer
      !! INTERFACES). `ks ≡ kt` unless `&ocean_ddiff_nml` double diffusion
      !! is enabled.
      type(c_ptr), intent(in), value :: c_handle
      type(c_ptr), intent(out) :: ptr
      integer(c_int), intent(out) :: nx, ny, nz, gen
      integer(c_int) :: status

      type(ocean_handle_t), pointer :: h
      real(wp), pointer :: tmp(:, :, :)

      ptr = c_null_ptr
      nx = 0
      ny = 0
      nz = 0
      gen = 0
      status = resolve_ocean(c_handle, h)
      if (status /= OCEAN_STATUS_OK) return
      tmp => h%state%vmix%ks
      call fill_getter_3d(ptr, nx, ny, nz, gen, tmp, h%state%dyn%outer_step_count)
   end function rdb_ocean_get_ks_ptr

   function rdb_ocean_get_wet_t_ptr(c_handle, ptr, nx, ny, gen) result(status) &
      bind(c, name="rdb_ocean_get_wet_t_ptr")
      !! Wet mask at T points (1 = wet, 0 = land): `metrics%wet_T`, shape
      !! `(nx_total, ny_total)`. Read-only (no setter — land masking is a
      !! configure-time / bathymetry concern).
      type(c_ptr), intent(in), value :: c_handle
      type(c_ptr), intent(out) :: ptr
      integer(c_int), intent(out) :: nx, ny, gen
      integer(c_int) :: status

      type(ocean_handle_t), pointer :: h
      real(wp), pointer :: tmp(:, :)

      ptr = c_null_ptr
      nx = 0
      ny = 0
      gen = 0
      status = resolve_ocean(c_handle, h)
      if (status /= OCEAN_STATUS_OK) return
      tmp => h%state%metrics%wet_T
      call fill_getter_2d(ptr, nx, ny, gen, tmp, h%state%dyn%outer_step_count)
   end function rdb_ocean_get_wet_t_ptr

   ! ================================================================
   ! P2 — tracers by name (not index; budget attribution is
   ! `tracer_t%budget_id`, not an index comparison, and the recovered
   ! code had none)
   ! ================================================================

   function rdb_ocean_get_tracer_count(c_handle, count_out) result(status) &
      bind(c, name="rdb_ocean_get_tracer_count")
      !! Number of registered tracers (S, T, + any passive tracers) —
      !! bounds for `rdb_ocean_list_tracers`'s index argument.
      type(c_ptr), intent(in), value :: c_handle
      integer(c_int), intent(out) :: count_out
      integer(c_int) :: status

      type(ocean_handle_t), pointer :: h

      count_out = 0
      status = resolve_ocean(c_handle, h)
      if (status /= OCEAN_STATUS_OK) return
      if (allocated(h%state%multilayer%tracers)) then
         count_out = int(size(h%state%multilayer%tracers), c_int)
      end if
   end function rdb_ocean_get_tracer_count

   function rdb_ocean_list_tracers(c_handle, idx, buf, cap) result(len_out) &
      bind(c, name="rdb_ocean_list_tracers")
      !! Name of the tracer at registry index `idx` (0-based). Same
      !! `snprintf`-style contract as `rdb_ocean_last_error`: copies up
      !! to `cap` bytes, NUL-terminates if room remains, returns the FULL
      !! trimmed name length. Returns 0 (buf untouched) on a bad handle,
      !! an out-of-range `idx`, or `cap <= 0`.
      type(c_ptr), intent(in), value :: c_handle
      integer(c_int), intent(in), value :: idx
      integer(c_int), intent(in), value :: cap
      character(kind=c_char), intent(out) :: buf(cap)
      integer(c_int) :: len_out

      type(ocean_handle_t), pointer :: h
      integer :: hstat, n, ncopy, i
      character(len=32) :: name

      len_out = 0_c_int
      hstat = resolve_ocean(c_handle, h)
      if (hstat /= OCEAN_STATUS_OK) return
      if (.not. allocated(h%state%multilayer%tracers)) return
      if (idx < 0_c_int .or. idx >= int(size(h%state%multilayer%tracers), c_int)) return
      if (cap <= 0_c_int) return

      name = h%state%multilayer%tracers(int(idx) + 1)%name
      n = len_trim(name)
      len_out = int(n, c_int)
      ncopy = min(n, int(cap))
      do i = 1, ncopy
         buf(i) = name(i:i)
      end do
      if (ncopy < cap) buf(ncopy + 1) = achar(0)
   end function rdb_ocean_list_tracers

   function rdb_ocean_get_tracer_ptr(c_handle, name, name_len, ptr, nx, ny, nz, gen) &
      result(status) bind(c, name="rdb_ocean_get_tracer_ptr")
      !! Raw `h*Tr` store (NOT concentration — D3.2) for the tracer named
      !! `name`, by NAME (never index — the recovered code had no tracer
      !! accessors, and index comparison is already the wrong idiom in the
      !! Fortran itself). `OCEAN_STATUS_ERR_NOT_FOUND` if no registered
      !! tracer matches. Shape `(nx_total, ny_total, nz_ml)`.
      type(c_ptr), intent(in), value :: c_handle
      integer(c_int), intent(in), value :: name_len
      character(kind=c_char), intent(in) :: name(name_len)
      type(c_ptr), intent(out) :: ptr
      integer(c_int), intent(out) :: nx, ny, nz, gen
      integer(c_int) :: status

      type(ocean_handle_t), pointer :: h
      real(wp), pointer :: tmp(:, :, :)
      character(len=:), allocatable :: fname
      integer :: it, found, ierr_local

      ptr = c_null_ptr
      nx = 0
      ny = 0
      nz = 0
      gen = 0
      status = resolve_ocean(c_handle, h)
      if (status /= OCEAN_STATUS_OK) return

      call c_to_f_string(name, name_len, fname)
      found = tracer_index_by_name(h, fname)
      if (found == 0) then
         call fail("rdb_ocean_get_tracer_ptr: no tracer named '"//fname//"'", &
                   ierr_local, OCEAN_STATUS_ERR_NOT_FOUND)
         status = int(ierr_local, c_int)
         return
      end if
      tmp => h%state%multilayer%tracers(found)%hTr
      call fill_getter_3d(ptr, nx, ny, nz, gen, tmp, h%state%dyn%outer_step_count)
   end function rdb_ocean_get_tracer_ptr

   ! ================================================================
   ! P2 — narrow setters. NARROW, NEVER BROAD (D3.5): each pushes only the
   ! device array(s) it just mutated, so calling one mid-run with no
   ! intervening read/refresh cannot rewind live device state to a stale
   ! host snapshot — the exact regression `ocean_setter_no_clobber` guards.
   ! No hTr -> Tr division happens here (D3.2) — the tracer setter takes
   ! concentration and forms `hTr = h_layer * value` on the physical
   ! interior using the CURRENT thickness (refreshed first); every other
   ! setter is a plain overwrite of its own array.
   ! ================================================================

   function rdb_ocean_set_h(c_handle, h_data, nx_p, ny_p, nz_p) result(status) &
      bind(c, name="rdb_ocean_set_h")
      !! Overwrite `h_layer` on the physical interior (k=1 bed .. k=nz
      !! surface), narrow-push it, then recompute `rho_layer` from the new
      !! thickness against whatever T/S currently sit on device (pulling
      !! the recomputed density back to host so the host stays
      !! authoritative). No other prognostic slot is touched.
      type(c_ptr), intent(in), value :: c_handle
      integer(c_int), intent(in), value :: nx_p, ny_p, nz_p
      real(c_double), intent(in) :: h_data(nx_p, ny_p, nz_p)
      integer(c_int) :: status

      type(ocean_handle_t), pointer :: h
      integer :: ng, i, j, k, ierr_local

      status = resolve_ocean(c_handle, h)
      if (status /= OCEAN_STATUS_OK) return
      if (.not. shape_matches_interior(h, nx_p, ny_p, nz_p)) then
         call fail("rdb_ocean_set_h: shape mismatch against the physical interior", &
                   ierr_local, OCEAN_STATUS_ERR_BAD_SHAPE)
         status = int(ierr_local, c_int)
         return
      end if

      ng = h%grid%nghost
      associate (hl => h%state%multilayer%h_layer)
         do k = 1, nz_p
            do j = 1, ny_p
               do i = 1, nx_p
                  hl(ng + i, ng + j, k) = real(h_data(i, j, k), wp)
               end do
            end do
         end do
         !$acc update device(hl)
      end associate
      call ocean_eos_compute(h%state%eos, h%state%multilayer)
      associate (rl => h%state%multilayer%rho_layer)
         !$acc update self(rl)
      end associate
      status = int(OCEAN_STATUS_OK, c_int)
   end function rdb_ocean_set_h

   function rdb_ocean_set_u(c_handle, u_data, nx_p, ny_p, nz_p) result(status) &
      bind(c, name="rdb_ocean_set_u")
      !! Overwrite `u_face_x_layer` on the physical interior faces
      !! (`nx_p+1` west faces bracketing `nx_p` physical columns) and
      !! narrow-push it. Does NOT re-derive `hu_face_x_layer` (a separate
      !! prognostic transport slot) — that stays whatever it was until the
      !! next step recomputes it.
      type(c_ptr), intent(in), value :: c_handle
      integer(c_int), intent(in), value :: nx_p, ny_p, nz_p
      real(c_double), intent(in) :: u_data(nx_p + 1, ny_p, nz_p)
      integer(c_int) :: status

      type(ocean_handle_t), pointer :: h
      integer :: ng, ierr_local

      status = resolve_ocean(c_handle, h)
      if (status /= OCEAN_STATUS_OK) return
      if (.not. shape_matches_interior(h, nx_p, ny_p, nz_p)) then
         call fail("rdb_ocean_set_u: shape mismatch against the physical interior", &
                   ierr_local, OCEAN_STATUS_ERR_BAD_SHAPE)
         status = int(ierr_local, c_int)
         return
      end if

      ng = h%grid%nghost
      associate (uf => h%state%multilayer%u_face_x_layer)
         uf(ng + 1:ng + nx_p + 1, ng + 1:ng + ny_p, :) = real(u_data, wp)
         !$acc update device(uf)
      end associate
      status = int(OCEAN_STATUS_OK, c_int)
   end function rdb_ocean_set_u

   function rdb_ocean_set_v(c_handle, v_data, nx_p, ny_p, nz_p) result(status) &
      bind(c, name="rdb_ocean_set_v")
      !! Overwrite `v_face_y_layer` on the physical interior faces
      !! (`ny_p+1` south faces bracketing `ny_p` physical rows) and
      !! narrow-push it. See `rdb_ocean_set_u` for the
      !! `hv_face_y_layer` caveat (not re-derived).
      type(c_ptr), intent(in), value :: c_handle
      integer(c_int), intent(in), value :: nx_p, ny_p, nz_p
      real(c_double), intent(in) :: v_data(nx_p, ny_p + 1, nz_p)
      integer(c_int) :: status

      type(ocean_handle_t), pointer :: h
      integer :: ng, ierr_local

      status = resolve_ocean(c_handle, h)
      if (status /= OCEAN_STATUS_OK) return
      if (.not. shape_matches_interior(h, nx_p, ny_p, nz_p)) then
         call fail("rdb_ocean_set_v: shape mismatch against the physical interior", &
                   ierr_local, OCEAN_STATUS_ERR_BAD_SHAPE)
         status = int(ierr_local, c_int)
         return
      end if

      ng = h%grid%nghost
      associate (vf => h%state%multilayer%v_face_y_layer)
         vf(ng + 1:ng + nx_p, ng + 1:ng + ny_p + 1, :) = real(v_data, wp)
         !$acc update device(vf)
      end associate
      status = int(OCEAN_STATUS_OK, c_int)
   end function rdb_ocean_set_v

   function rdb_ocean_set_bathymetry(c_handle, b_data, nx_p, ny_p) result(status) &
      bind(c, name="rdb_ocean_set_bathymetry")
      !! Overwrite bathymetry `b` (m, POSITIVE DOWN) on the physical
      !! interior, fill ghosts, re-wrap the periodic/fold seam, re-derive
      !! `bt_H_ref = b` (the mode-split contract) and re-wrap IT too, then
      !! narrow-push `b` (+ `bt_H_ref`) — mirrors the init-time sequence at
      !! `rdb_driver.F90` exactly, so a mid-run perturbation sees the same
      !! seam/BT-reference treatment the setup path does. No prognostic
      !! slot is touched.
      !!
      !! CAVEAT (mid-run, not a re-seed): does NOT re-derive the ALE
      !! `z_ref` table (ZSTAR_FULL) or the static land mask + C-grid face
      !! metrics — for those, drive the case from a namelist/restart
      !! instead (matches the recovered precedent's documented limit).
      type(c_ptr), intent(in), value :: c_handle
      integer(c_int), intent(in), value :: nx_p, ny_p
      real(c_double), intent(in) :: b_data(nx_p, ny_p)
      integer(c_int) :: status

      type(ocean_handle_t), pointer :: h
      integer :: ng, ierr_local
      logical :: wrap_seam

      status = resolve_ocean(c_handle, h)
      if (status /= OCEAN_STATUS_OK) return
      if (nx_p /= h%grid%nx_phys .or. ny_p /= h%grid%ny_phys) then
         call fail("rdb_ocean_set_bathymetry: shape mismatch against the physical interior", &
                   ierr_local, OCEAN_STATUS_ERR_BAD_SHAPE)
         status = int(ierr_local, c_int)
         return
      end if

      ng = h%grid%nghost
      wrap_seam = h%state%bc%periodic_x .or. h%state%bc%periodic_y .or. h%state%bc%north_fold

      associate (b => h%state%barotropic%b)
         b(ng + 1:ng + nx_p, ng + 1:ng + ny_p) = real(b_data, wp)
         call fill_bathymetry_ghosts_array(b, h%grid)
         if (wrap_seam) then
            call ocean_periodic_wrap_centre_2d(b, h%grid%nx_total, h%grid%ny_total, &
                                               h%grid%nx_phys, h%grid%ny_phys, h%grid%nghost, &
                                               h%state%bc%periodic_x, h%state%bc%periodic_y)
            call ocean_fold_wrap_eta_2d(h%grid, h%state%bc, b)
         end if
         !$acc update device(b)
      end associate

      if (h%state%dyn%n_inner >= 1) then
         associate (bref => h%state%dyn%bt_work%bt_H_ref)
            bref = h%state%barotropic%b
            if (wrap_seam) then
               call ocean_periodic_wrap_centre_2d(bref, h%grid%nx_total, h%grid%ny_total, &
                                                  h%grid%nx_phys, h%grid%ny_phys, h%grid%nghost, &
                                                  h%state%bc%periodic_x, h%state%bc%periodic_y)
               call ocean_fold_wrap_eta_2d(h%grid, h%state%bc, bref)
            end if
            !$acc update device(bref)
         end associate
      end if
      status = int(OCEAN_STATUS_OK, c_int)
   end function rdb_ocean_set_bathymetry

   function rdb_ocean_set_wind(c_handle, taux_data, tauy_data, nx_p, ny_p) result(status) &
      bind(c, name="rdb_ocean_set_wind")
      !! Overwrite the C-grid surface wind stress (N/m^2) mid-run:
      !! `taux_data` is `(nx_p+1, ny_p)` (east faces), `tauy_data` is
      !! `(nx_p, ny_p+1)` (north faces) over the physical interior; ghost
      !! faces are left untouched (0 from init — wall faces see no
      !! spurious stress). Narrow-pushes `tau_x`/`tau_y` only. Settable
      !! repeatedly for a time-varying wind schedule.
      type(c_ptr), intent(in), value :: c_handle
      integer(c_int), intent(in), value :: nx_p, ny_p
      real(c_double), intent(in) :: taux_data(nx_p + 1, ny_p)
      real(c_double), intent(in) :: tauy_data(nx_p, ny_p + 1)
      integer(c_int) :: status

      type(ocean_handle_t), pointer :: h
      integer :: ng, ierr_local

      status = resolve_ocean(c_handle, h)
      if (status /= OCEAN_STATUS_OK) return
      if (nx_p /= h%grid%nx_phys .or. ny_p /= h%grid%ny_phys) then
         call fail("rdb_ocean_set_wind: shape mismatch against the physical interior", &
                   ierr_local, OCEAN_STATUS_ERR_BAD_SHAPE)
         status = int(ierr_local, c_int)
         return
      end if

      ng = h%grid%nghost
      associate (tx => h%state%surface_stress%tau_x, ty => h%state%surface_stress%tau_y)
         tx(ng + 1:ng + nx_p + 1, ng + 1:ng + ny_p) = real(taux_data, wp)
         ty(ng + 1:ng + nx_p, ng + 1:ng + ny_p + 1) = real(tauy_data, wp)
         !$acc update device(tx, ty)
      end associate
      status = int(OCEAN_STATUS_OK, c_int)
   end function rdb_ocean_set_wind

   function rdb_ocean_set_heat_flux(c_handle, q_data, nx_p, ny_p) result(status) &
      bind(c, name="rdb_ocean_set_heat_flux")
      !! Overwrite net surface heat flux `h%state%surface_flux%Q_heat`
      !! (W/m^2) on the physical interior, narrow-push. See
      !! `rdb_ocean_get_q_heat_ptr` for the P2.4 unification note.
      type(c_ptr), intent(in), value :: c_handle
      integer(c_int), intent(in), value :: nx_p, ny_p
      real(c_double), intent(in) :: q_data(nx_p, ny_p)
      integer(c_int) :: status

      type(ocean_handle_t), pointer :: h
      integer :: ng, ierr_local

      status = resolve_ocean(c_handle, h)
      if (status /= OCEAN_STATUS_OK) return
      if (nx_p /= h%grid%nx_phys .or. ny_p /= h%grid%ny_phys) then
         call fail("rdb_ocean_set_heat_flux: shape mismatch against the physical interior", &
                   ierr_local, OCEAN_STATUS_ERR_BAD_SHAPE)
         status = int(ierr_local, c_int)
         return
      end if

      ng = h%grid%nghost
      associate (qh => h%state%surface_flux%Q_heat)
         qh(ng + 1:ng + nx_p, ng + 1:ng + ny_p) = real(q_data, wp)
         !$acc update device(qh)
      end associate
      status = int(OCEAN_STATUS_OK, c_int)
   end function rdb_ocean_set_heat_flux

   function rdb_ocean_set_salt_flux(c_handle, q_data, nx_p, ny_p) result(status) &
      bind(c, name="rdb_ocean_set_salt_flux")
      !! Overwrite net surface salt flux `h%state%surface_flux%Q_salt` on
      !! the physical interior, narrow-push. See
      !! `rdb_ocean_get_q_heat_ptr` for the P2.4 unification note.
      type(c_ptr), intent(in), value :: c_handle
      integer(c_int), intent(in), value :: nx_p, ny_p
      real(c_double), intent(in) :: q_data(nx_p, ny_p)
      integer(c_int) :: status

      type(ocean_handle_t), pointer :: h
      integer :: ng, ierr_local

      status = resolve_ocean(c_handle, h)
      if (status /= OCEAN_STATUS_OK) return
      if (nx_p /= h%grid%nx_phys .or. ny_p /= h%grid%ny_phys) then
         call fail("rdb_ocean_set_salt_flux: shape mismatch against the physical interior", &
                   ierr_local, OCEAN_STATUS_ERR_BAD_SHAPE)
         status = int(ierr_local, c_int)
         return
      end if

      ng = h%grid%nghost
      associate (qs => h%state%surface_flux%Q_salt)
         qs(ng + 1:ng + nx_p, ng + 1:ng + ny_p) = real(q_data, wp)
         !$acc update device(qs)
      end associate
      status = int(OCEAN_STATUS_OK, c_int)
   end function rdb_ocean_set_salt_flux

   function rdb_ocean_set_tracer(c_handle, name, name_len, data, nx_p, ny_p, nz_p) &
      result(status) bind(c, name="rdb_ocean_set_tracer")
      !! Set the tracer named `name` to `data` (its own units — degC for
      !! temperature, PSU for salinity) on the physical interior. Verbatim
      !! the recovered `ocean_set_tracer_impl` sequence, generalised from
      !! hardcoded S/T to any registered tracer BY NAME: flush the
      !! windowed-advection accumulator -> refresh host (need the LIVE
      !! `h_layer` to form `hTr = h_layer*value`) -> mutate `hTr` on the
      !! physical interior -> narrow-push `hTr` -> recompute `rho_layer`
      !! (pulled back to host so it stays authoritative).
      !! `OCEAN_STATUS_ERR_NOT_FOUND` if no registered tracer matches.
      type(c_ptr), intent(in), value :: c_handle
      integer(c_int), intent(in), value :: name_len
      character(kind=c_char), intent(in) :: name(name_len)
      integer(c_int), intent(in), value :: nx_p, ny_p, nz_p
      real(c_double), intent(in) :: data(nx_p, ny_p, nz_p)
      integer(c_int) :: status

      type(ocean_handle_t), pointer :: h
      character(len=:), allocatable :: fname
      integer :: ng, i, j, k, it, ierr_local

      status = resolve_ocean(c_handle, h)
      if (status /= OCEAN_STATUS_OK) return
      if (.not. shape_matches_interior(h, nx_p, ny_p, nz_p)) then
         call fail("rdb_ocean_set_tracer: shape mismatch against the physical interior", &
                   ierr_local, OCEAN_STATUS_ERR_BAD_SHAPE)
         status = int(ierr_local, c_int)
         return
      end if

      call c_to_f_string(name, name_len, fname)
      it = tracer_index_by_name(h, fname)
      if (it == 0) then
         call fail("rdb_ocean_set_tracer: no tracer named '"//fname//"'", &
                   ierr_local, OCEAN_STATUS_ERR_NOT_FOUND)
         status = int(ierr_local, c_int)
         return
      end if

      ! Drain any open windowed tracer-advection accumulation BEFORE
      ! reading/overwriting the tracer, so pending transport (accumulated
      ! against the OLD field) is never later applied to the NEW one.
      ! Hard no-op at dt_tracer_advect_ratio<=1 (default; bit-identical).
      call ocean_dyn_flush_tracer_window(h%grid, h%state%metrics, h%state%dyn, &
                                         h%state%continuity, h%state%multilayer, &
                                         bc=h%state%bc)
      call ocean_handle_refresh_host(h)

      ng = h%grid%nghost
      associate (hl => h%state%multilayer%h_layer, htr => h%state%multilayer%tracers(it)%hTr)
         do k = 1, nz_p
            do j = 1, ny_p
               do i = 1, nx_p
                  htr(ng + i, ng + j, k) = hl(ng + i, ng + j, k)*real(data(i, j, k), wp)
               end do
            end do
         end do
         !$acc update device(htr)
      end associate
      call ocean_eos_compute(h%state%eos, h%state%multilayer)
      associate (rl => h%state%multilayer%rho_layer)
         !$acc update self(rl)
      end associate
      h%host_is_current = .true.
      status = int(OCEAN_STATUS_OK, c_int)
   end function rdb_ocean_set_tracer

   ! ================================================================
   ! P2 — scalar diagnostics
   ! ================================================================

   function rdb_ocean_get_kinetic_energy(c_handle, ke_out) result(status) &
      bind(c, name="rdb_ocean_get_kinetic_energy")
      !! Total kinetic energy (J, up to the Boussinesq reference-density
      !! factor — matches `rdb_ocean_get_total_mass`'s convention of
      !! leaving `rho0` out) over the physical interior:
      !! `sum(0.5 * h_layer * (u_centre^2 + v_centre^2)) * dx * dy`, faces
      !! averaged to centres — same formula as
      !! `rdb_ocean_budgets::budget_total_ke`, computed inline (cartesian
      !! `dx`/`dy` only; P2.5 generalises to curvilinear `areaT`). Unlike
      !! the raw-pointer getters above, this refreshes the host itself —
      !! it hands back a NUMBER, not a pointer a caller could otherwise
      !! defer syncing for.
      type(c_ptr), intent(in), value :: c_handle
      real(c_double), intent(out) :: ke_out
      integer(c_int) :: status

      type(ocean_handle_t), pointer :: h
      real(wp) :: total, uc, vc
      integer :: i, j, k, ng, nxp, nyp

      ke_out = 0.0_c_double
      status = resolve_ocean(c_handle, h)
      if (status /= OCEAN_STATUS_OK) return

      call ocean_handle_refresh_host(h)
      ng = h%grid%nghost
      nxp = h%grid%nx_phys
      nyp = h%grid%ny_phys
      total = 0.0_wp
      associate (ms => h%state%multilayer)
         do k = 1, ms%nz_ml
            do j = 1, nyp
               do i = 1, nxp
                  uc = 0.5_wp*(ms%u_face_x_layer(ng + i, ng + j, k) + &
                               ms%u_face_x_layer(ng + i + 1, ng + j, k))
                  vc = 0.5_wp*(ms%v_face_y_layer(ng + i, ng + j, k) + &
                               ms%v_face_y_layer(ng + i, ng + j + 1, k))
                  total = total + ms%h_layer(ng + i, ng + j, k)*0.5_wp*(uc*uc + vc*vc)
               end do
            end do
         end do
      end associate
      ke_out = real(total*h%grid%dx*h%grid%dy, c_double)
      status = int(OCEAN_STATUS_OK, c_int)
   end function rdb_ocean_get_kinetic_energy

   ! ================================================================
   ! P7 — diagnostics discoverability. Two STATIC catalogs (no handle
   ! needed — a Python user should be able to ask "what CAN I request"
   ! before ever calling create()) plus the LIVE registered set on one
   ! handle ("what IS this instance actually emitting"). Mirrors the
   ! `rdb_ocean_list_tracers` name-by-index / snprintf-style contract
   ! throughout: copies up to `cap` bytes, NUL-terminates if room
   ! remains, returns the FULL trimmed name length; 0 (buf untouched) on
   ! an out-of-range index or `cap <= 0`.
   ! ================================================================

   function rdb_ocean_derived_catalog_size() result(n) &
      bind(c, name="rdb_ocean_derived_catalog_size")
      !! Number of names in the DERIVED diagnostic catalog (vorticity_z,
      !! ke_total, mld_density, ... — `rdb_ocean_diag_derived`'s static
      !! table, opt-in via `&ocean_diag_nml diags`). No handle required:
      !! this is build-time information, reachable before `create()`.
      integer(c_int) :: n
      n = int(derived_catalog_size(), c_int)
   end function rdb_ocean_derived_catalog_size

   function rdb_ocean_derived_catalog_name(idx, buf, cap) result(len_out) &
      bind(c, name="rdb_ocean_derived_catalog_name")
      !! Name of derived-catalog entry `idx` (0-based, < the size above).
      integer(c_int), intent(in), value :: idx
      integer(c_int), intent(in), value :: cap
      character(kind=c_char), intent(out) :: buf(cap)
      integer(c_int) :: len_out

      character(len=64) :: name
      integer :: n, ncopy, i

      len_out = 0_c_int
      if (idx < 0_c_int .or. idx >= int(derived_catalog_size(), c_int)) return
      if (cap <= 0_c_int) return

      name = derived_catalog_name(int(idx) + 1)
      n = len_trim(name)
      len_out = int(n, c_int)
      ncopy = min(n, int(cap))
      do i = 1, ncopy
         buf(i) = name(i:i)
      end do
      if (ncopy < cap) buf(ncopy + 1) = achar(0)
   end function rdb_ocean_derived_catalog_name

   function rdb_ocean_canonical_catalog_size() result(n) &
      bind(c, name="rdb_ocean_canonical_catalog_size")
      !! Number of names in the CANONICAL diagnostic catalog (SSH,
      !! temperature, salinity, u, v, KE, ... — the set
      !! `register_default_diags` MAY register at setup; some entries
      !! are gated by another namelist group, e.g. `temperature` needs
      !! `&ocean_thermo_nml enable_thermodynamics`). No handle required.
      !! Use `rdb_ocean_get_diag_count`/`rdb_ocean_list_diags` on a
      !! live handle to see what actually registered.
      integer(c_int) :: n
      n = int(canonical_diag_catalog_size(), c_int)
   end function rdb_ocean_canonical_catalog_size

   function rdb_ocean_canonical_catalog_name(idx, buf, cap) result(len_out) &
      bind(c, name="rdb_ocean_canonical_catalog_name")
      !! Name of canonical-catalog entry `idx` (0-based, < the size above).
      integer(c_int), intent(in), value :: idx
      integer(c_int), intent(in), value :: cap
      character(kind=c_char), intent(out) :: buf(cap)
      integer(c_int) :: len_out

      character(len=16) :: name
      integer :: n, ncopy, i

      len_out = 0_c_int
      if (idx < 0_c_int .or. idx >= int(canonical_diag_catalog_size(), c_int)) return
      if (cap <= 0_c_int) return

      name = canonical_diag_catalog_name(int(idx) + 1)
      n = len_trim(name)
      len_out = int(n, c_int)
      ncopy = min(n, int(cap))
      do i = 1, ncopy
         buf(i) = name(i:i)
      end do
      if (ncopy < cap) buf(ncopy + 1) = achar(0)
   end function rdb_ocean_canonical_catalog_name

   function rdb_ocean_get_diag_count(c_handle, count_out) result(status) &
      bind(c, name="rdb_ocean_get_diag_count")
      !! Number of diagnostics REGISTERED on this live instance right
      !! now (canonical + derived + anything the `&ocean_diag_nml diags`
      !! token list added) — bounds for `rdb_ocean_list_diags`'s
      !! index argument. This is the "selected" set, as opposed to the
      !! two static "available" catalogs above.
      type(c_ptr), intent(in), value :: c_handle
      integer(c_int), intent(out) :: count_out
      integer(c_int) :: status

      type(ocean_handle_t), pointer :: h

      count_out = 0
      status = resolve_ocean(c_handle, h)
      if (status /= OCEAN_STATUS_OK) return
      count_out = int(h%state%diag%nvars, c_int)
   end function rdb_ocean_get_diag_count

   function rdb_ocean_list_diags(c_handle, idx, buf, cap) result(len_out) &
      bind(c, name="rdb_ocean_list_diags")
      !! Name of the registered diagnostic at index `idx` (0-based).
      !! Same snprintf-style contract as `rdb_ocean_list_tracers`.
      type(c_ptr), intent(in), value :: c_handle
      integer(c_int), intent(in), value :: idx
      integer(c_int), intent(in), value :: cap
      character(kind=c_char), intent(out) :: buf(cap)
      integer(c_int) :: len_out

      type(ocean_handle_t), pointer :: h
      integer :: hstat, n, ncopy, i
      character(len=64) :: name

      len_out = 0_c_int
      hstat = resolve_ocean(c_handle, h)
      if (hstat /= OCEAN_STATUS_OK) return
      if (idx < 0_c_int .or. idx >= int(h%state%diag%nvars, c_int)) return
      if (cap <= 0_c_int) return

      name = h%state%diag%vars(int(idx) + 1)%name
      n = len_trim(name)
      len_out = int(n, c_int)
      ncopy = min(n, int(cap))
      do i = 1, ncopy
         buf(i) = name(i:i)
      end do
      if (ncopy < cap) buf(ncopy + 1) = achar(0)
   end function rdb_ocean_list_diags

   ! ================================================================
   ! P7 — in-memory diagnostic access. Returns the diag manager's OWN
   ! `output_buffer` for a registered diagnostic, by name — no NetCDF
   ! round-trip. `output_buffer` is pulled host-ward unconditionally at
   ! every cadence fire (`ocean_diag_step`, the `!$acc update self
   ! ... if_present` right before `v%fire_count` increments), so unlike
   ! the P2 state getters this needs no separate `refresh_host` call:
   ! by the time `rdb_ocean_step` returns, every diagnostic that
   ! fired this call is already host-current. `fire_count` (not
   ! `dyn%outer_step_count`) is the generation — a diagnostic on a
   ! multi-hour cadence should not look "stale" to Python between fires.
   ! ================================================================

   function rdb_ocean_get_diagnostic_ptr(c_handle, name, name_len, ptr, &
                                         nx, ny, nz, gen) result(status) &
      bind(c, name="rdb_ocean_get_diagnostic_ptr")
      !! Raw `output_buffer` for the diagnostic named `name` (LAYER vgrid:
      !! `(nx_total, ny_total, nz_ml)` for a layered var, `(nx_total,
      !! ny_total, 1)` for a 2D var; a non-LAYER `output_vgrid` reports
      !! the remapped shape, e.g. `nz` z-levels). `OCEAN_STATUS_ERR_NOT_FOUND`
      !! if `name` is not currently REGISTERED on this instance (it may
      !! still be a legal name on one of the two static catalogs, just
      !! gated off or not selected — see `rdb_ocean_get_diag_count`/
      !! `rdb_ocean_list_diags` to discover what IS registered).
      type(c_ptr), intent(in), value :: c_handle
      integer(c_int), intent(in), value :: name_len
      character(kind=c_char), intent(in) :: name(name_len)
      type(c_ptr), intent(out) :: ptr
      integer(c_int), intent(out) :: nx, ny, nz, gen
      integer(c_int) :: status

      type(ocean_handle_t), pointer :: h
      real(wp), pointer :: tmp(:, :, :)
      character(len=:), allocatable :: fname
      integer :: iv, found, ierr_local

      ptr = c_null_ptr
      nx = 0
      ny = 0
      nz = 0
      gen = 0
      status = resolve_ocean(c_handle, h)
      if (status /= OCEAN_STATUS_OK) return

      call c_to_f_string(name, name_len, fname)
      found = 0
      do iv = 1, h%state%diag%nvars
         if (trim(h%state%diag%vars(iv)%name) == fname) then
            found = iv
            exit
         end if
      end do
      if (found == 0) then
         call fail("rdb_ocean_get_diagnostic_ptr: no diagnostic named '"//fname// &
                   "' is currently registered (see rdb_ocean_list_diags)", &
                   ierr_local, OCEAN_STATUS_ERR_NOT_FOUND)
         status = int(ierr_local, c_int)
         return
      end if
      if (.not. allocated(h%state%diag%vars(found)%output_buffer)) then
         call fail("rdb_ocean_get_diagnostic_ptr: '"//fname// &
                   "' has no output_buffer allocated", ierr_local, OCEAN_STATUS_ERR_NOT_FOUND)
         status = int(ierr_local, c_int)
         return
      end if
      tmp => h%state%diag%vars(found)%output_buffer
      ptr = c_loc(tmp(1, 1, 1))
      nx = int(size(tmp, 1), c_int)
      ny = int(size(tmp, 2), c_int)
      nz = int(size(tmp, 3), c_int)
      gen = int(h%state%diag%vars(found)%fire_count, c_int)
   end function rdb_ocean_get_diagnostic_ptr

   ! ================================================================
   ! Internal helpers
   ! ================================================================

   function resolve_ocean(c_handle, h) result(status)
      !! Resolve + verify a handle is a live, fully-initialised ocean
      !! simulation. Every query/step entry point funnels through this so
      !! the bad-handle and not-yet-initialised cases are reported with one
      !! consistent status code each, in one place.
      type(c_ptr), intent(in) :: c_handle
      type(ocean_handle_t), pointer, intent(out) :: h
      integer(c_int) :: status

      integer :: hstat

      hstat = handle_check(c_handle, h)
      if (hstat /= HANDLE_OK) then
         h => null()
         status = int(OCEAN_STATUS_ERR_BAD_HANDLE, c_int)
         return
      end if
      if (.not. h%is_initialised) then
         h => null()
         status = int(OCEAN_STATUS_ERR_NOT_INITIALISED, c_int)
         return
      end if
      status = int(OCEAN_STATUS_OK, c_int)
   end function resolve_ocean

   function resolve_ocean_pending(c_handle, h) result(status)
      !! P2.5: resolve + verify a handle is a live, PENDING (not yet
      !! finalised) ocean simulation — the window every
      !! `rdb_ocean_stage_*` geometry-injection call requires. Mirrors
      !! `resolve_ocean`'s shape, checking `is_pending` instead of
      !! `is_initialised`.
      type(c_ptr), intent(in) :: c_handle
      type(ocean_handle_t), pointer, intent(out) :: h
      integer(c_int) :: status

      integer :: hstat

      hstat = handle_check(c_handle, h)
      if (hstat /= HANDLE_OK) then
         h => null()
         status = int(OCEAN_STATUS_ERR_BAD_HANDLE, c_int)
         return
      end if
      if (.not. h%is_pending) then
         h => null()
         status = int(OCEAN_STATUS_ERR_NOT_PENDING, c_int)
         return
      end if
      status = int(OCEAN_STATUS_OK, c_int)
   end function resolve_ocean_pending

   subroutine c_to_f_string(c_str, c_len, f_str)
      !! Convert a `bind(c)` `character(kind=c_char)` buffer + explicit
      !! length into a Fortran allocatable string. No null-termination
      !! assumption — `c_len` is authoritative (matches the recovered
      !! precedent's convention).
      integer(c_int), intent(in) :: c_len
      character(kind=c_char), intent(in) :: c_str(c_len)
      character(len=:), allocatable, intent(out) :: f_str

      integer :: i

      allocate (character(len=c_len) :: f_str)
      do i = 1, c_len
         f_str(i:i) = c_str(i)
      end do
   end subroutine c_to_f_string

   subroutine fill_getter_3d(ptr, nx, ny, nz, gen, arr, step)
      !! Shared tail of every 3D P2 getter: base-address `c_loc`, actual
      !! array extents (never assumed from grid metadata — read straight
      !! off the pointer), and the generation stamp. `arr` must already be
      !! pointer-associated with a live, non-empty state array.
      type(c_ptr), intent(out) :: ptr
      integer(c_int), intent(out) :: nx, ny, nz, gen
      real(wp), pointer, intent(in) :: arr(:, :, :)
      integer, intent(in) :: step

      ptr = c_loc(arr(1, 1, 1))
      nx = int(size(arr, 1), c_int)
      ny = int(size(arr, 2), c_int)
      nz = int(size(arr, 3), c_int)
      gen = int(step, c_int)
   end subroutine fill_getter_3d

   subroutine fill_getter_2d(ptr, nx, ny, gen, arr, step)
      !! 2D counterpart of `fill_getter_3d`.
      type(c_ptr), intent(out) :: ptr
      integer(c_int), intent(out) :: nx, ny, gen
      real(wp), pointer, intent(in) :: arr(:, :)
      integer, intent(in) :: step

      ptr = c_loc(arr(1, 1))
      nx = int(size(arr, 1), c_int)
      ny = int(size(arr, 2), c_int)
      gen = int(step, c_int)
   end subroutine fill_getter_2d

   subroutine ocean_handle_refresh_host(h)
      !! Lazy D->H refresh of every P2-exposed leaf array, gated on
      !! `h%host_is_current`. `!$acc update self` on LEAF component names
      !! only, via `associate` — never the aggregate `ocean_state_t`/
      !! sub-state derived type (`rdb_ocean_dyn.F90:2949`). Inert on a
      !! host-only build.
      type(ocean_handle_t), pointer, intent(inout) :: h

      integer :: it

      if (h%host_is_current) return

      associate (ms => h%state%multilayer, bt => h%state%barotropic, &
                 bw => h%state%dyn%bt_work, ss => h%state%surface_stress, &
                 vm => h%state%vmix, mt => h%state%metrics, sf => h%state%surface_flux)
         !$acc update self(ms%h_layer, ms%u_face_x_layer, ms%v_face_y_layer, &
         !$acc&            ms%hu_face_x_layer, ms%hv_face_y_layer, &
         !$acc&            ms%w_interface, ms%rho_layer)
         if (allocated(ms%tracers)) then
            do it = 1, size(ms%tracers)
               !$acc update self(ms%tracers(it)%hTr)
            end do
         end if
         !$acc update self(bt%b, bw%bt_eta)
         !$acc update self(ss%tau_x, ss%tau_y)
         !$acc update self(sf%Q_heat, sf%Q_salt)
         !$acc update self(vm%kv, vm%kt, vm%ks)
         !$acc update self(mt%wet_T)
      end associate
      h%host_is_current = .true.
   end subroutine ocean_handle_refresh_host

   pure function tracer_index_by_name(h, name) result(idx)
      !! Registry index (1-based) of the tracer named `name`, or 0 if none
      !! matches. Name comparison is `trim`-both-sides (registry names are
      !! fixed-length `character(len=32)`).
      type(ocean_handle_t), pointer, intent(in) :: h
      character(len=*), intent(in) :: name
      integer :: idx

      integer :: it

      idx = 0
      if (.not. allocated(h%state%multilayer%tracers)) return
      do it = 1, size(h%state%multilayer%tracers)
         if (trim(h%state%multilayer%tracers(it)%name) == trim(name)) then
            idx = it
            return
         end if
      end do
   end function tracer_index_by_name

   pure function shape_matches_interior(h, nx_p, ny_p, nz_p) result(ok)
      !! True iff `(nx_p, ny_p, nz_p)` matches the live handle's physical
      !! interior shape (`grid%nx_phys`, `grid%ny_phys`,
      !! `multilayer%nz_ml`) — the shape every P2 3D setter's caller-
      !! supplied array must have.
      type(ocean_handle_t), pointer, intent(in) :: h
      integer(c_int), intent(in) :: nx_p, ny_p, nz_p
      logical :: ok

      ok = (nx_p == h%grid%nx_phys) .and. (ny_p == h%grid%ny_phys) .and. &
           (nz_p == h%state%multilayer%nz_ml)
   end function shape_matches_interior

end module rdb_ocean_api
