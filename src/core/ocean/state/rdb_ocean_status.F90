!! Named status codes for the ocean solver-creation path's optional `ierr`
!! out-arguments (P0 of the Python runtime API plan,
!! `docs/ocean_python_api_plan.md`).
module rdb_ocean_status
   !! A library cannot abort its host process. Every procedure on the
   !! solver-creation path (config parse/validate, the `configure_ocean_*`
   !! setup chain, IC seeding) that used to `error stop` now takes an
   !! `optional, intent(out) :: ierr` — present ⇒ return one of these codes
   !! instead of aborting; absent ⇒ unchanged legacy `error stop` behaviour
   !! (every existing caller keeps aborting exactly as it does today).
   !!
   !! Coarse, stage-level codes: the specific reason is always logged via
   !! `global_logger%error` immediately before `ierr` is set (unchanged
   !! message text), so the code only needs to say which STAGE failed —
   !! mirroring the existing `restart_mismatch_kind` precedent in
   !! `rdb_ocean_state.F90`.
   !!
   !! `OCEAN_STATUS_ERR_IO`'s "a dimension mismatch" clause (below) is
   !! deliberate, not the "taxonomy drift" an early adversarial review
   !! (P0.1, F4) flagged: a mismatch between a FILE's declared dimensions
   !! and the model's config is a file-content problem, distinct from the
   !! namelist-only checks (`rad_earth`, domain size, tripolar knobs) that
   !! `configure_ocean_metrics` itself resolves as `OCEAN_STATUS_ERR_SETUP`
   !! with no file involved at all.
   !!
   !! Layering note (P0.1 review F8): this module is registered in the
   !! CMake "core infrastructure" block (alongside `rdb_config`/`rdb_eos`,
   !! both of which `use` it) even though it lives under
   !! `src/core/ocean/state/` and is named `rdb_ocean_*`. This is a
   !! DELIBERATE choice, not an oversight: this tree is ocean-only
   !! (`sim_type='ocean'` is the only accepted regime — see the repo
   !! CLAUDE.md), so "the one status vocabulary" and "the one regime's
   !! status vocabulary" are the same thing today, and a rename to a
   !! regime-neutral `rdb_status` would be pure churn until a second
   !! regime returns to this tree. Revisit if/when that happens.
   implicit none
   private

   public :: OCEAN_STATUS_OK
   public :: OCEAN_STATUS_ERR_CONFIG_PARSE
   public :: OCEAN_STATUS_ERR_CONFIG_VALIDATE
   public :: OCEAN_STATUS_ERR_SETUP
   public :: OCEAN_STATUS_ERR_IC_SEED
   public :: OCEAN_STATUS_ERR_IO
   public :: OCEAN_STATUS_ERR_BAD_HANDLE
   public :: OCEAN_STATUS_ERR_ALREADY_EXISTS
   public :: OCEAN_STATUS_ERR_NOT_INITIALISED
   public :: OCEAN_STATUS_ERR_NOT_FOUND
   public :: OCEAN_STATUS_ERR_BAD_SHAPE
   public :: OCEAN_STATUS_ERR_BATHYMETRY_SIGN
   public :: OCEAN_STATUS_ERR_NOT_PENDING
   public :: OCEAN_STATUS_ERR_RESTART_SCHEMA
   public :: OCEAN_STATUS_ERR_RESTART_DECOMP
   public :: OCEAN_STATUS_ERR_RESTART_GRID

   integer, parameter :: OCEAN_STATUS_OK = 0
      !! No error.
   integer, parameter :: OCEAN_STATUS_ERR_CONFIG_PARSE = 1
      !! `read_config`/`read_config_from_string`: strict namelist/schema
      !! parse failure (unknown group/key, type/range/enum violation).
   integer, parameter :: OCEAN_STATUS_ERR_CONFIG_VALIDATE = 2
      !! `validate_config`: a cross-knob semantic check failed.
   integer, parameter :: OCEAN_STATUS_ERR_SETUP = 3
      !! The `configure_ocean_*` setup chain (`rdb_ocean_setup.F90`) or one
      !! of its direct callees (metrics/EOS/BC validation) rejected the
      !! resolved configuration.
   integer, parameter :: OCEAN_STATUS_ERR_IC_SEED = 4
      !! Initial-condition seeding (`ocean_state_seed_from_cfg` and its
      !! analytical-IC helpers) rejected the resolved configuration.
   integer, parameter :: OCEAN_STATUS_ERR_IO = 5
      !! A setup-time file read failed (missing NetCDF at build time, a
      !! malformed supergrid/z-level-IC file, a dimension mismatch).

   ! ---- Handle-specific codes (Python runtime API P1, `rdb_ocean_api`) ----
   ! Numbered from 10 so the config/setup range above (1-5) has room to grow
   ! without colliding with the handle-lifecycle range.
   integer, parameter :: OCEAN_STATUS_ERR_BAD_HANDLE = 10
      !! The `c_ptr` handle is null, stale (already destroyed), or does not
      !! carry the ocean-handle magic sentinel (`rdb_handle`).
   integer, parameter :: OCEAN_STATUS_ERR_ALREADY_EXISTS = 11
      !! `rdb_ocean_create_from_string` was called while a live ocean
      !! handle already exists. Multi-instance is impossible in this tree
      !! (several modules hold process-global topology/device state keyed to
      !! "the one running ocean simulation" — see `rdb_ocean_api`'s header)
      !! so a second create fails loud instead of silently corrupting the
      !! first.
   integer, parameter :: OCEAN_STATUS_ERR_NOT_INITIALISED = 12
      !! The handle resolved (valid magic) but never completed `create()`
      !! (device mapping not done) — a step/query call landed on it anyway.
   integer, parameter :: OCEAN_STATUS_ERR_NOT_FOUND = 13
      !! P2 accessors: `rdb_ocean_get_tracer_ptr`/`rdb_ocean_set_tracer`
      !! were given a name that is not in the tracer registry.
   integer, parameter :: OCEAN_STATUS_ERR_BAD_SHAPE = 14
      !! P2 setters: the caller-supplied array extents do not match the
      !! live state's physical interior shape.
   integer, parameter :: OCEAN_STATUS_ERR_BATHYMETRY_SIGN = 15
      !! P2.5 geometry injection (`rdb_ocean_bathymetry_inject`): a staged
      !! bathymetry array normalised to zero wet cells under the given sign
      !! convention — almost certainly the OTHER convention was intended
      !! (D6.2: the single most dangerous argument in the geometry API).
   integer, parameter :: OCEAN_STATUS_ERR_NOT_PENDING = 16
      !! P2.5 geometry injection: a `rdb_ocean_stage_*` call (or
      !! `rdb_ocean_create_finalize`) landed on a handle that is not
      !! currently in the pending (staging) window — either it was never
      !! created via `rdb_ocean_create_pending`, or `create_finalize`
      !! already completed it.

   ! ---- Restart-resume codes (`ocean_state_restart_read`) ----
   ! Numbered from 20, disjoint from BOTH the config/setup range (1-5) and
   ! the handle-lifecycle range (10-12).  Named codes for
   ! `restart_mismatch_kind`'s internal 1/2/3 (P0.1 review F4: the
   ! PUBLIC `ierr` this routine hands back to `create()` must not reuse
   ! the SAME integers as `OCEAN_STATUS_ERR_CONFIG_PARSE` /
   ! `_CONFIG_VALIDATE` / `_SETUP` — the INTERNAL 1/2/3 convention shared
   ! by `ocean_restart_check_decomp`/`ocean_restart_read_local`
   ! (`rdb_ocean_restart_io.F90`) stays as-is; only the value surfaced
   ! through THIS module's public `ierr` is translated).
   integer, parameter :: OCEAN_STATUS_ERR_RESTART_SCHEMA = 20
      !! Restart file's `schema_version` attribute does not match the
      !! build's `RESTART_SCHEMA_VERSION`.
   integer, parameter :: OCEAN_STATUS_ERR_RESTART_DECOMP = 21
      !! Restart file's decomposition metadata (`px`/`py`/`i_start`/
      !! `j_start`/`nx_local`/`ny_local`) does not match the live `decomp`.
   integer, parameter :: OCEAN_STATUS_ERR_RESTART_GRID = 22
      !! Restart file's grid/vcoord/tracer metadata (`nz_ml`/`nghost`/
      !! `vcoord_type`/`n_tracers`/tracer indices) does not match the
      !! live state.

end module rdb_ocean_status
