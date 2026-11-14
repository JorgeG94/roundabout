!! The ice step sequencer of the CICE-port sea ice (`scheme = "cice"`; roundabout-original).
module rdb_ice_step_cice
   !! `ice_step_cice` owns the ORDER of one ice step, the cadence switches and the
   !! error enforcement between stage groups.  It computes nothing: every stage is a
   !! call to an existing driver, or, for a stage whose implementation another PR
   !! delivers, a PROCEDURE HOOK that defaults to "stage missing".  It lives in the
   !! MAIN library (not `rdb_ice_core_objs`) because it takes the state slots
   !! (`ice_cice_state_t` composes the ocean grid and restart registry).  The engine
   !! calls it once per outer step (`ice_cice_step`, `rdb_ice_cice_hooks`).
   !!
   !! **The stage MANIFEST** (`ice_step_plan_t` lists it as data; the test asserts it
   !! against a literal list and checks `res%executed` against it; the step itself is
   !! straight-line code in this order, not a loop over the manifest).  One ice step,
   !! `dt` = the outer ocean step, `ndtd = 1`, every outer step:
   !!
   !! | # | stage | kind | group |
   !! |---|-------|------|-------|
   !! | 1 | `ocean_import` | hook | prep |
   !! | 2 | `velocity_to_T` | hook, only if `highfreq` | prep |
   !! | 3 | `atm_prepare` | hook | prep |
   !! | 4 | `flux_zero` | `ice_flux_zero` | prep |
   !! | 5 | `save_init` | `ice_save_init` | prep |
   !! | 6 | `radiation` | `ice_step_sw_therm1`, first half | thermo |
   !! | 7 | `therm1` | `ice_step_sw_therm1`, second half | thermo |
   !! | 8 | `therm2` | `ice_step_therm2` | thermo |
   !! | 9 | `halo_state` | hook | dynamics |
   !! | 10 | `aggregate_ghost` | `ice_update_state(.true.)` | dynamics |
   !! | 11 | `stress_to_faces` | hook | dynamics |
   !! | 12 | `dynamics` | hook (strength + EVP) | dynamics |
   !! | 13 | `transport` | hook, post mode NONE | transport |
   !! | 14 | `ridging` | hook | close |
   !! | 15 | `halo_state_post` | hook | close |
   !! | 16 | `aggregate` | `ice_update_state(.true.)` | close |
   !! | 17 | `ocean_handoff` | hook | out |
   !!
   !! Stages 6 and 7 are ONE call (`ice_step_sw_therm1`, the single owner of radiation,
   !! vertical thermodynamics and the through-ice shortwave merge); they stay two manifest
   !! entries so the radiation placement (on the step-start state, after the snapshot) is
   !! visible, and its two status values tell which half failed.
   !!
   !! **Why this order.**
   !! - *Thermodynamics before dynamics* (the reference model's order within a step).
   !!   The ridging consumes the deformation rates of THIS step's dynamics and the
   !!   transported state, and the next step's basal heat reads this step's ice-ocean
   !!   stress, so the dynamics comes last among the state changes.
   !! - *Inputs first, then zero, then save.*  The three input hooks write only the
   !!   work slot, so every refusal that depends on them (the forcing heights) happens
   !!   before the first state or flux write; the flux slot is then zeroed once and
   !!   every producer adds; the snapshot (`aicen_init`) is the weight of every merge
   !!   and the input of the linear ITD remap; the shortwave runs on the step-start
   !!   state (`aicen == aicen_init`, `sw_lag = .false.`).
   !! - *therm1 before therm2.*  therm1 produces `rside`, which the lateral melt of
   !!   therm2 consumes; therm2 (rain, remap, new ice, lateral melt, cleanup) ends with
   !!   ONE aggregate (the cleanup's).
   !! - *Halo, ghost aggregate, strength, EVP.*  The dynamics reads the category state
   !!   and the aggregates in the ghost ring (the strength at the face neighbours, the
   !!   face ice areas), so after therm2 the state ring is refreshed (`halo_state`:
   !!   exchange, periodic wrap and fold on one rank too) and the aggregates are rebuilt
   !!   over the ghosted extent (`aggregate_ghost`).  The `dynamics` hook then computes
   !!   `ice_compute_strength` from that post-thermodynamics, haloed state and runs the
   !!   EVP.  `stress_to_faces` (therm1's `strair_x_T/y_T` onto the faces) precedes it.
   !! - *Transport, ridging, halo, aggregate.*  Transport runs with post mode
   !!   `ICE_TRANSPORT_POST_NONE`; the error word is enforced right after it so the
   !!   ridging never runs on a state transport flagged; nothing aggregates between
   !!   transport and ridging (the ridging closes the ADVECTED open water `aice0`; an
   !!   aggregate before it would replace it by `1 - aice`, and its own cleanup removes
   !!   the slivers transport creates).  Transport updates the owned cells only, so the
   !!   ring is refreshed (`halo_state_post`) before the final ghosted aggregate.
   !! - *`ocean_handoff` last.*  The flux averager accumulates EVERY step and applies
   !!   into the ocean's surface fluxes only on the ocean thermo step
   !!   (`ctx%thermo_step`); the ocean is debited only through `fhocn`.
   !!
   !! **Slot for the tendency diagnostics.**  A snapshot of the state before the
   !! dynamics (for the dynamics/thermodynamics tendency diagnostics) belongs between
   !! `aggregate_ghost` and `stress_to_faces`; it is not a stage yet.  A PR adding it
   !! inserts one manifest entry there (and its id), nothing else moves.
   !!
   !! **The freezing potential and the frazil bank (contract).**  `frzmlt` is ONE signed
   !! field.  The `ocean_import` hook forms it EVERY ice step from the CURRENT ocean-side
   !! bank held in `work%bank_ocean`: copy it into `work%bank`, set the positive part
   !! `frzmlt = work%bank/dt` (zero once spent) and the negative part per the coupling
   !! layer's melt rule.  The new-ice stage (in therm2) spends `work%bank` (zeroes it where
   !! it made ice).  The SINGLE write-back owner is the `ocean_handoff` hook: it copies the
   !! spent `work%bank` back into the ocean-side bank (`work%bank_ocean`, which the ocean
   !! frazil provider aliases) every step, so a bank is frozen exactly once however many
   !! ice steps a thermo window spans.  The ocean is debited only through `fhocn`, never
   !! through `frzmlt`.
   !!
   !! **The ice-ocean stress (contract).**  `work%tau_x, work%tau_y` (the stress of the
   !! ice on the ocean per ICE area at T, last step's dynamics output) are PERSISTENT
   !! state, not scratch: written by the `dynamics` hook, read by the basal heat of the
   !! next step, never zeroed by the sequencer, registered for restart by
   !! `work%register_restart` (tags `ice_strocnTx`, `ice_strocnTy`).  The dynamics
   !! coupling PR owns the field's physics; it may move the arrays into its own slot but
   !! must keep the tags.
   !!
   !! **The wet mask.**  `work%wet` is produced by the setup from the land mask through
   !! `work%set_wet`; the step is refused if it was never set.
   !!
   !! **Cadence.**  The whole ice step runs every outer step.  `thermo_step`
   !! (`dyn%is_thermo_step()` in the engine) is handed to the hooks as `ctx%thermo_step`;
   !! only `ocean_handoff` acts on it (apply the averaged flux), `ocean_import` forms
   !! `frzmlt` every step.
   !!
   !! **Errors.**  The cell error word is cleared BEFORE the first stage, so a bit raised
   !! by the import hook survives to the first enforcement point.  The kernels skip a cell
   !! that carries a bit of the COMPILE-TIME fatal mask (`ICE_ERR_FATAL_MASK`); the driver
   !! stops the step on `opts%fatal_mask`, which also holds the bits `fatal_on_error`
   !! promotes (e.g. the ITD boundary-order bit).  Such a promoted bit therefore stops the
   !! step at the end of the stage group that raised it, not inside the group.  The word
   !! is reduced (`ice_err_reduce`, one host sync) after therm2, after the dynamics, after
   !! transport and after the final aggregate.  A fatal bit stops the step before the
   !! next stage (`ICE_STEP_FATAL`; `res%summary` holds the reduction; the CALLER aborts
   !! with `ice_err_enforce`).  A hook or driver that returns a non-zero status stops the
   !! step at that stage (`ICE_STEP_STAGE_FAILED`).
   !!
   !! **Refusals before anything is written.**  `ice_step_cice` checks the step
   !! envelope first and returns without writing a single state or flux value: a
   !! missing REQUIRED stage (named in `res%message`), a test double without
   !! `opts%allow_test_doubles`, `ndtd /= 1`, more than one rank, `nghost < 3`, an ice
   !! step different from the outer step, a non-finite or non-positive `dt`, layer or
   !! category counts outside the compile-time maxima, an unmapped or disagreeing slot,
   !! an unconfigured radiation slot, a wet mask never set.  The forcing heights
   !! `zlvl, zlvs` (filled by `atm_prepare`) are checked right after the input hooks,
   !! still before the first write.  The default hook returns `ICE_HOOK_MISSING`; it is
   !! never a silent skip.
   !!
   !! **Hook contract (what the in-flight PRs must satisfy).**  Every hook is
   !! `subroutine(ctx, status)`; `status = 0` is success.  `ctx` carries the slots by
   !! pointer (`ctx%cice, ctx%flux, ctx%rad, ctx%work`), `ctx%dt`, `ctx%istep`,
   !! `ctx%thermo_step` and `ctx%transport_post` (the post mode the `transport` hook MUST
   !! pass to `ice_step_transport`).  All arrays are device-mapped; a hook launches
   !! kernels and never copies state to the host.
   !! - `ocean_import` (ocean import PR): `work%sst`, `work%Tf` (ice-side linear freezing
   !!   point), `work%bank` and the signed `work%frzmlt` (above).
   !! - `velocity_to_T` (dynamics PR; only when `p%atm%highfreq`): `work%uvel, work%vvel`.
   !! - `atm_prepare` (atmosphere PR): `work%uatm, vatm, wind, potT, Qa, rhoa, flw, fsnow,
   !!   frain, sw_down` and the scalars `work%zlvl, work%zlvs`.
   !! - `halo_state` and `halo_state_post` (halo/multi-rank PR): refresh the ghost ring
   !!   of the category state (`aicen, vicen, vsnon, Tsfcn, qicen, sicen, qsnon`) and
   !!   `aice0` (exchange on split axes, periodic wrap, fold; also on one rank).
   !! - `stress_to_faces` and `dynamics` (dynamics-coupling PR): the first maps
   !!   `work%strair_x_T/y_T` onto the faces; the second computes `ice_compute_strength`
   !!   from the post-thermodynamics, halo-refreshed category state, runs the EVP
   !!   (`ice_evp_c_step`) and fills `work%rdg_conv, work%rdg_shear` and the persistent
   !!   `work%tau_x, work%tau_y`.
   !! - `transport`, `ridging` (the wiring PR: adapters over `ice_step_transport` with
   !!   `post = ctx%transport_post` and `ice_step_ridge`, which also adds the
   !!   `ICE_ERR_RIDGE_PASSES` bit from `niter`).
   !! - `ocean_handoff` (the flux averager PR): accumulate `ctx%flux` every step; apply
   !!   when `ctx%thermo_step`; write the spent `work%bank` back to `work%bank_ocean`.
   !!
   !! Teardown: `ice_sw_wiring_cleanup` (module `rdb_ice_sw_wiring`) releases the zero
   !! snow-absorption array `ice_step_sw_therm1` allocates lazily.
   use, intrinsic :: iso_fortran_env, only: int32, int64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_constants, only: wp
   use rdb_mem_report, only: arr_bytes
   use rdb_ice_params, only: ice_salinity_profile
   use rdb_ice_dims, only: ICE_NCAT_MAX, ICE_NILYR_MAX, ICE_NSLYR_MAX
   use rdb_ocean_restart, only: restart_registry_t
   use rdb_ice_cice_state, only: ice_cice_state_t, ice_update_state
   use rdb_ice_flux_state, only: ice_flux_t, ice_flux_zero, ice_save_init, ice_slots_agree
   use rdb_ice_rad_state, only: ice_rad_state_t
   use rdb_ice_radiation, only: ICE_RAD_OK
   use rdb_ice_sw_wiring, only: ice_step_sw_therm1
   use rdb_ice_therm2, only: ice_step_therm2
   use rdb_ice_transport_cice, only: ICE_TRANSPORT_POST_NONE
   use rdb_ice_errors, only: ice_err_summary_t, ice_err_reduce, ice_err_is_fatal, ice_err_clear, ICE_ERR_FATAL_MASK
   implicit none
   private

   public :: ice_step_cice, ice_step_plan_build
   public :: ice_step_plan_t, ice_step_hooks_t, ice_step_opts_t, ice_step_result_t, ice_step_ctx_t, ice_step_work_t
   public :: ice_stage_missing, ice_stage_name
   public :: ICE_NSTAGE, ICE_STAGE_OCEAN_IMPORT, ICE_STAGE_VELOCITY_TO_T, ICE_STAGE_ATM_PREPARE
   public :: ICE_STAGE_FLUX_ZERO, ICE_STAGE_SAVE_INIT, ICE_STAGE_RADIATION, ICE_STAGE_THERM1, ICE_STAGE_THERM2
   public :: ICE_STAGE_HALO_STATE, ICE_STAGE_AGGREGATE_GHOST, ICE_STAGE_STRESS_TO_FACES, ICE_STAGE_DYNAMICS
   public :: ICE_STAGE_TRANSPORT, ICE_STAGE_RIDGING, ICE_STAGE_HALO_STATE_POST, ICE_STAGE_AGGREGATE
   public :: ICE_STAGE_OCEAN_HANDOFF
   public :: ICE_GROUP_PREP, ICE_GROUP_THERMO, ICE_GROUP_DYNAMICS, ICE_GROUP_TRANSPORT, ICE_GROUP_CLOSE, ICE_GROUP_OUT
   public :: ICE_STEP_OK, ICE_STEP_REFUSED_MISSING, ICE_STEP_REFUSED_DOUBLE, ICE_STEP_REFUSED_ENVELOPE
   public :: ICE_STEP_STAGE_FAILED, ICE_STEP_FATAL, ICE_HOOK_MISSING
   public :: ICE_STEP_INIT_OK, ICE_STEP_INIT_BAD

   integer, parameter :: ICE_NSTAGE = 17
      !! Number of stages of the manifest.
   integer, parameter :: ICE_STAGE_OCEAN_IMPORT = 1
   integer, parameter :: ICE_STAGE_VELOCITY_TO_T = 2
   integer, parameter :: ICE_STAGE_ATM_PREPARE = 3
   integer, parameter :: ICE_STAGE_FLUX_ZERO = 4
   integer, parameter :: ICE_STAGE_SAVE_INIT = 5
   integer, parameter :: ICE_STAGE_RADIATION = 6
   integer, parameter :: ICE_STAGE_THERM1 = 7
   integer, parameter :: ICE_STAGE_THERM2 = 8
   integer, parameter :: ICE_STAGE_HALO_STATE = 9
   integer, parameter :: ICE_STAGE_AGGREGATE_GHOST = 10
   integer, parameter :: ICE_STAGE_STRESS_TO_FACES = 11
   integer, parameter :: ICE_STAGE_DYNAMICS = 12
   integer, parameter :: ICE_STAGE_TRANSPORT = 13
   integer, parameter :: ICE_STAGE_RIDGING = 14
   integer, parameter :: ICE_STAGE_HALO_STATE_POST = 15
   integer, parameter :: ICE_STAGE_AGGREGATE = 16
   integer, parameter :: ICE_STAGE_OCEAN_HANDOFF = 17

   integer, parameter :: ICE_GROUP_PREP = 1
      !! Inputs, zero, snapshot (nothing here raises a cell error bit).
   integer, parameter :: ICE_GROUP_THERMO = 2
      !! Radiation, therm1, therm2: the error word is enforced after it.
   integer, parameter :: ICE_GROUP_DYNAMICS = 3
      !! Halo, ghost aggregate, stress to faces, EVP: enforced after it.
   integer, parameter :: ICE_GROUP_TRANSPORT = 4
      !! Transport alone: enforced after it, before the ridging.
   integer, parameter :: ICE_GROUP_CLOSE = 5
      !! Ridging, halo, final aggregate: enforced after it.
   integer, parameter :: ICE_GROUP_OUT = 6
      !! The hand-off to the ocean.

   integer, parameter :: ICE_STEP_OK = 0
      !! The step ran to the end.
   integer, parameter :: ICE_STEP_REFUSED_MISSING = 1
      !! A required stage has no implementation; nothing was touched.
   integer, parameter :: ICE_STEP_REFUSED_DOUBLE = 2
      !! A test double is installed and `allow_test_doubles` is off; nothing was touched.
   integer, parameter :: ICE_STEP_REFUSED_ENVELOPE = 3
      !! The step envelope is violated (`ndtd`, ranks, `nghost`, `dt`, slots); nothing was touched.
   integer, parameter :: ICE_STEP_STAGE_FAILED = 4
      !! A stage returned a non-zero status (`res%stage`); the step stopped there.
   integer, parameter :: ICE_STEP_FATAL = 5
      !! A fatal bit is in the cell error word after a stage group (`res%summary`).
   integer, parameter :: ICE_HOOK_MISSING = -1
      !! Status the default hook returns.
   integer, parameter :: ICE_STEP_INIT_OK = 0
   integer, parameter :: ICE_STEP_INIT_BAD = 1

   type :: ice_step_work_t
      !! Inputs and scratch of the stages the sequencer calls itself, owned once
      !! (`init`), device-mapped once (`enter_data`).  `(nx, ny)` unless noted;
      !! ghosted extent.  The hooks fill the input rows; the stages read them.
      logical :: is_init = .false.
      logical :: mapped = .false.
      integer :: nx = 0
      integer :: ny = 0
      integer :: ncat = 0
      integer :: nilyr = 0
      integer :: nslyr = 0
      logical :: wet_set = .false.
         !! `set_wet` was called.
      real(wp) :: zlvl = 10.0_wp
         !! Height of the wind/temperature/humidity forcing [m] (set by `atm_prepare`).
      real(wp) :: zlvs = 10.0_wp
         !! Height of the snow-flux reference [m].
      ! ---- inputs ----
      real(wp), allocatable :: sw_down(:, :)
         !! Broadband downward shortwave [W m-2] (atmosphere).
      real(wp), allocatable :: uatm(:, :), vatm(:, :), wind(:, :)
         !! Wind components and speed [m s-1] (atmosphere).
      real(wp), allocatable :: potT(:, :), Qa(:, :), rhoa(:, :)
         !! Potential temperature [K], specific humidity [kg/kg], density [kg m-3].
      real(wp), allocatable :: flw(:, :)
         !! Downward longwave [W m-2].
      real(wp), allocatable :: fsnow(:, :), frain(:, :)
         !! Snowfall and rainfall [kg m-2 s-1] per cell area.
      real(wp), allocatable :: uvel(:, :), vvel(:, :)
         !! T-averaged ice velocity of the last step [m s-1] (`highfreq`).
      real(wp), allocatable :: Tf(:, :), sst(:, :)
         !! Ice-side freezing point and sea-surface temperature [degC] (ocean).
      real(wp), allocatable :: frzmlt(:, :), bank(:, :)
         !! Signed freezing/melt potential [W m-2] and the frazil bank [J m-2] (ocean).
      real(wp), allocatable :: bank_ocean(:, :)
         !! The ocean-side frazil bank [J m-2]: read by `ocean_import`, written back by `ocean_handoff`.
      real(wp), allocatable :: tau_x(:, :), tau_y(:, :)
         !! PERSISTENT: last step's ice-ocean stress per ICE area at T [N m-2]; written by `dynamics`,
         !! restart-registered by `register_restart`, never zeroed by the sequencer.
      logical, allocatable :: wet(:, :)
         !! Ocean (wet) cells; set by `set_wet` (the setup's land mask), never defaulted.
      real(wp), allocatable :: salinz(:)
         !! Salinity profile of new ice `(nilyr)` (filled by `init` from the BL99 profile).
      ! ---- outputs and scratch ----
      real(wp), allocatable :: strair_x_T(:, :), strair_y_T(:, :)
         !! Wind stress on the ice, merged over the categories (therm1).
      real(wp), allocatable :: strairxn(:, :, :), strairyn(:, :, :)
         !! Per-category wind stress `(nx, ny, ncat)`.
      real(wp), allocatable :: fbot(:, :), rside(:, :), wlat(:, :)
         !! Basal heat flux, lateral-melt fraction (therm1 -> therm2) and its rate.
      integer(int32), allocatable :: ierr_n(:, :, :)
         !! Per-category error words of therm1 `(nx, ny, ncat)`.
      real(wp), allocatable :: Tsfc_c(:, :), qice_c(:, :, :), sice_c(:, :, :), qsno_c(:, :, :)
         !! Cell-mean tracer scratch of the cleanup (therm2, ridging).
      real(wp), allocatable :: rdg_conv(:, :), rdg_shear(:, :)
         !! Ridging convergence and shear rates [1/s] of this step's dynamics.
      integer(int32), allocatable :: niter(:, :)
         !! Ridging pass counts (ridging).
   contains
      procedure, non_overridable :: init => work_init
      procedure, non_overridable :: destroy => work_destroy
      procedure, non_overridable :: enter_data => work_enter_data
      procedure, non_overridable :: exit_data => work_exit_data
      procedure, non_overridable :: bytes => step_work_bytes
      procedure, non_overridable :: set_wet => work_set_wet
      procedure, non_overridable :: register_restart => work_register_restart
   end type ice_step_work_t

   type :: ice_step_ctx_t
      !! What a hook sees (module docstring, "Hook contract").
      type(ice_cice_state_t), pointer :: cice => null()
      type(ice_flux_t), pointer :: flux => null()
      type(ice_rad_state_t), pointer :: rad => null()
      type(ice_step_work_t), pointer :: work => null()
      real(wp) :: dt = 0.0_wp
         !! The ice step [s] (== the outer step).
      integer :: istep = 0
         !! The caller's step counter.
      logical :: thermo_step = .false.
         !! The ocean thermodynamics step (`dyn%is_thermo_step()`).
      integer :: transport_post = ICE_TRANSPORT_POST_NONE
         !! Post mode the `transport` hook must pass to `ice_step_transport`.
   end type ice_step_ctx_t

   abstract interface
      subroutine ice_stage_hook_i(ctx, status)
         import :: ice_step_ctx_t
         implicit none
         type(ice_step_ctx_t), intent(inout) :: ctx
         integer, intent(out) :: status
      end subroutine ice_stage_hook_i
   end interface

   type :: ice_step_hooks_t
      !! The stage implementations the sequencer does not own.  Every slot defaults
      !! to `ice_stage_missing`.  Install with `install`; `is_double` marks a test
      !! double (refused unless `opts%allow_test_doubles`).
      procedure(ice_stage_hook_i), pointer, nopass :: ocean_import => ice_stage_missing
      procedure(ice_stage_hook_i), pointer, nopass :: velocity_to_T => ice_stage_missing
      procedure(ice_stage_hook_i), pointer, nopass :: atm_prepare => ice_stage_missing
      procedure(ice_stage_hook_i), pointer, nopass :: halo_state => ice_stage_missing
      procedure(ice_stage_hook_i), pointer, nopass :: stress_to_faces => ice_stage_missing
      procedure(ice_stage_hook_i), pointer, nopass :: dynamics => ice_stage_missing
      procedure(ice_stage_hook_i), pointer, nopass :: transport => ice_stage_missing
      procedure(ice_stage_hook_i), pointer, nopass :: ridging => ice_stage_missing
      procedure(ice_stage_hook_i), pointer, nopass :: halo_state_post => ice_stage_missing
      procedure(ice_stage_hook_i), pointer, nopass :: ocean_handoff => ice_stage_missing
      logical :: is_double(ICE_NSTAGE) = .false.
   contains
      procedure :: install => hooks_install
      procedure :: installed => hooks_installed
      procedure :: reset => hooks_reset
   end type ice_step_hooks_t

   type :: ice_step_opts_t
      !! Step envelope and switches.
      integer :: ndtd = 1
         !! Dynamics substeps per ice step (must be 1).
      integer :: nranks = 1
         !! MPI ranks (must be 1 until the multi-rank PRs).
      real(wp) :: dt_outer = -1.0_wp
         !! The outer ocean step [s]; `> 0` must equal `dt`; `<= 0` skips the comparison.
      logical :: allow_test_doubles = .false.
         !! Accept hooks installed with `is_double`.  Tests only.
      integer(int32) :: fatal_mask = ICE_ERR_FATAL_MASK
         !! Bits that stop the step (`ice_err_fatal_set(fatal_on_error)`).
   end type ice_step_opts_t

   type :: ice_step_plan_t
      !! The stage list in order, as data.
      integer :: n = 0
      integer :: id(ICE_NSTAGE) = 0
      logical :: required(ICE_NSTAGE) = .false.
         !! The stage will run (a not-required stage is skipped by configuration).
      logical :: present(ICE_NSTAGE) = .false.
         !! The stage has an implementation (a driver stage is always present).
      logical :: is_hook(ICE_NSTAGE) = .false.
      logical :: is_double(ICE_NSTAGE) = .false.
      integer :: group(ICE_NSTAGE) = 0
      logical :: enforce_after(ICE_NSTAGE) = .false.
         !! The cell error word is reduced and enforced right after this stage.
      integer :: transport_post = ICE_TRANSPORT_POST_NONE
      logical :: dynamics_computes_strength = .true.
         !! The `dynamics` hook computes the ice strength from the post-thermodynamics, halo-refreshed state.
      character(len=24) :: strength_routine = "ice_compute_strength"
         !! The routine that does it.
      logical :: ridging_uses_advected_aice0 = .true.
         !! Nothing aggregates (rebuilds `aice0`) between transport and ridging: the ridging
         !! closes the open water the transport advected.
   end type ice_step_plan_t

   type :: ice_step_result_t
      !! What a step did.
      integer :: status = ICE_STEP_OK
      integer :: stage = 0
         !! The stage that refused, failed or ended the step (0: none).
      integer :: stage_status = 0
         !! The raw status of that stage (hook status, driver code).
      character(len=160) :: message = ""
      integer :: n_exec = 0
         !! Number of stages executed.
      integer :: executed(ICE_NSTAGE) = 0
         !! Stage ids in the order they ran.
      type(ice_err_summary_t) :: summary
         !! The last error-word reduction.
      integer(int32) :: err_mask = 0_int32
         !! OR of the masks of every reduction of the step.
      type(ice_step_plan_t) :: plan
   end type ice_step_result_t

contains

   ! ------------------------------------------------------------------
   ! hooks
   ! ------------------------------------------------------------------

   subroutine ice_stage_missing(ctx, status)
      !! The default of every hook slot: the stage has no implementation.
      type(ice_step_ctx_t), intent(inout) :: ctx
      integer, intent(out) :: status

      status = ICE_HOOK_MISSING
   end subroutine ice_stage_missing

   subroutine hooks_install(this, stage, proc, is_double)
      !! Install `proc` as the implementation of `stage` (`ICE_STAGE_*`, a hook stage).
      class(ice_step_hooks_t), intent(inout) :: this
      integer, intent(in) :: stage
      procedure(ice_stage_hook_i) :: proc
      logical, intent(in), optional :: is_double

      select case (stage)
      case (ICE_STAGE_OCEAN_IMPORT)
         this%ocean_import => proc
      case (ICE_STAGE_VELOCITY_TO_T)
         this%velocity_to_T => proc
      case (ICE_STAGE_ATM_PREPARE)
         this%atm_prepare => proc
      case (ICE_STAGE_HALO_STATE)
         this%halo_state => proc
      case (ICE_STAGE_HALO_STATE_POST)
         this%halo_state_post => proc
      case (ICE_STAGE_STRESS_TO_FACES)
         this%stress_to_faces => proc
      case (ICE_STAGE_DYNAMICS)
         this%dynamics => proc
      case (ICE_STAGE_TRANSPORT)
         this%transport => proc
      case (ICE_STAGE_RIDGING)
         this%ridging => proc
      case (ICE_STAGE_OCEAN_HANDOFF)
         this%ocean_handoff => proc
      case default
         return
      end select
      this%is_double(stage) = .false.
      if (present(is_double)) this%is_double(stage) = is_double
   end subroutine hooks_install

   subroutine hooks_reset(this)
      !! Every slot back to the `ice_stage_missing` default (a structure constructor with no
      !! arguments would do the same but crashes nvfortran 26.5).
      class(ice_step_hooks_t), intent(inout) :: this

      this%ocean_import => ice_stage_missing
      this%velocity_to_T => ice_stage_missing
      this%atm_prepare => ice_stage_missing
      this%halo_state => ice_stage_missing
      this%halo_state_post => ice_stage_missing
      this%stress_to_faces => ice_stage_missing
      this%dynamics => ice_stage_missing
      this%transport => ice_stage_missing
      this%ridging => ice_stage_missing
      this%ocean_handoff => ice_stage_missing
      this%is_double = .false.
   end subroutine hooks_reset

   function hooks_installed(this, stage) result(yes)
      !! The slot of `stage` holds something other than the default.
      class(ice_step_hooks_t), intent(in) :: this
      integer, intent(in) :: stage
      logical :: yes

      select case (stage)
      case (ICE_STAGE_OCEAN_IMPORT)
         yes = .not. associated(this%ocean_import, ice_stage_missing)
      case (ICE_STAGE_VELOCITY_TO_T)
         yes = .not. associated(this%velocity_to_T, ice_stage_missing)
      case (ICE_STAGE_ATM_PREPARE)
         yes = .not. associated(this%atm_prepare, ice_stage_missing)
      case (ICE_STAGE_HALO_STATE)
         yes = .not. associated(this%halo_state, ice_stage_missing)
      case (ICE_STAGE_HALO_STATE_POST)
         yes = .not. associated(this%halo_state_post, ice_stage_missing)
      case (ICE_STAGE_STRESS_TO_FACES)
         yes = .not. associated(this%stress_to_faces, ice_stage_missing)
      case (ICE_STAGE_DYNAMICS)
         yes = .not. associated(this%dynamics, ice_stage_missing)
      case (ICE_STAGE_TRANSPORT)
         yes = .not. associated(this%transport, ice_stage_missing)
      case (ICE_STAGE_RIDGING)
         yes = .not. associated(this%ridging, ice_stage_missing)
      case (ICE_STAGE_OCEAN_HANDOFF)
         yes = .not. associated(this%ocean_handoff, ice_stage_missing)
      case default
         yes = .true.   ! a driver stage
      end select
   end function hooks_installed

   pure function ice_stage_name(stage) result(name)
      integer, intent(in) :: stage
      character(len=16) :: name

      select case (stage)
      case (ICE_STAGE_OCEAN_IMPORT)
         name = "ocean_import"
      case (ICE_STAGE_FLUX_ZERO)
         name = "flux_zero"
      case (ICE_STAGE_SAVE_INIT)
         name = "save_init"
      case (ICE_STAGE_VELOCITY_TO_T)
         name = "velocity_to_T"
      case (ICE_STAGE_ATM_PREPARE)
         name = "atm_prepare"
      case (ICE_STAGE_RADIATION)
         name = "radiation"
      case (ICE_STAGE_THERM1)
         name = "therm1"
      case (ICE_STAGE_THERM2)
         name = "therm2"
      case (ICE_STAGE_HALO_STATE)
         name = "halo_state"
      case (ICE_STAGE_AGGREGATE_GHOST)
         name = "aggregate_ghost"
      case (ICE_STAGE_HALO_STATE_POST)
         name = "halo_state_post"
      case (ICE_STAGE_STRESS_TO_FACES)
         name = "stress_to_faces"
      case (ICE_STAGE_DYNAMICS)
         name = "dynamics"
      case (ICE_STAGE_TRANSPORT)
         name = "transport"
      case (ICE_STAGE_RIDGING)
         name = "ridging"
      case (ICE_STAGE_AGGREGATE)
         name = "aggregate"
      case (ICE_STAGE_OCEAN_HANDOFF)
         name = "ocean_handoff"
      case default
         name = "?"
      end select
   end function ice_stage_name

   ! ------------------------------------------------------------------
   ! the plan
   ! ------------------------------------------------------------------

   function ice_step_plan_build(hooks, highfreq) result(plan)
      !! The stage list for these hooks.  `highfreq` (`p%atm%highfreq`) makes
      !! `velocity_to_T` required; otherwise it is not run.
      type(ice_step_hooks_t), intent(in) :: hooks
      logical, intent(in) :: highfreq
      type(ice_step_plan_t) :: plan
      integer :: k

      plan%n = ICE_NSTAGE
      plan%transport_post = ICE_TRANSPORT_POST_NONE
      do k = 1, ICE_NSTAGE
         plan%id(k) = k    ! the stage ids ARE the order
         plan%required(k) = .true.
         plan%present(k) = hooks%installed(k)
         plan%is_double(k) = hooks%is_double(k) .and. plan%present(k)
      end do
      plan%required(ICE_STAGE_VELOCITY_TO_T) = highfreq
      plan%is_hook([ICE_STAGE_OCEAN_IMPORT, ICE_STAGE_VELOCITY_TO_T, ICE_STAGE_ATM_PREPARE, ICE_STAGE_HALO_STATE, &
                    ICE_STAGE_STRESS_TO_FACES, ICE_STAGE_DYNAMICS, ICE_STAGE_TRANSPORT, ICE_STAGE_RIDGING, &
                    ICE_STAGE_HALO_STATE_POST, ICE_STAGE_OCEAN_HANDOFF]) = .true.
      plan%group([ICE_STAGE_OCEAN_IMPORT, ICE_STAGE_VELOCITY_TO_T, ICE_STAGE_ATM_PREPARE, ICE_STAGE_FLUX_ZERO, &
                  ICE_STAGE_SAVE_INIT]) = ICE_GROUP_PREP
      plan%group([ICE_STAGE_RADIATION, ICE_STAGE_THERM1, ICE_STAGE_THERM2]) = ICE_GROUP_THERMO
      plan%group([ICE_STAGE_HALO_STATE, ICE_STAGE_AGGREGATE_GHOST, ICE_STAGE_STRESS_TO_FACES, &
                  ICE_STAGE_DYNAMICS]) = ICE_GROUP_DYNAMICS
      plan%group(ICE_STAGE_TRANSPORT) = ICE_GROUP_TRANSPORT
      plan%group([ICE_STAGE_RIDGING, ICE_STAGE_HALO_STATE_POST, ICE_STAGE_AGGREGATE]) = ICE_GROUP_CLOSE
      plan%group(ICE_STAGE_OCEAN_HANDOFF) = ICE_GROUP_OUT
      plan%enforce_after([ICE_STAGE_THERM2, ICE_STAGE_DYNAMICS, ICE_STAGE_TRANSPORT, ICE_STAGE_AGGREGATE]) = .true.
   end function ice_step_plan_build

   ! ------------------------------------------------------------------
   ! work slot
   ! ------------------------------------------------------------------

   subroutine work_init(this, cice, ierr)
      !! Allocate the slot for the extents and the salinity profile of `cice`
      !! (host-zeroed; `wet` is NOT set: call `set_wet`).
      class(ice_step_work_t), intent(inout) :: this
      type(ice_cice_state_t), intent(in) :: cice
      integer, intent(out) :: ierr
      integer :: nx, ny, nc, nl, ns
      real(wp), allocatable :: sp(:)

      ierr = ICE_STEP_INIT_OK
      if (.not. cice%is_init) then
         ierr = ICE_STEP_INIT_BAD
         return
      end if
      if (this%is_init) call this%destroy()
      nx = cice%nx
      ny = cice%ny
      nc = cice%ncat
      nl = cice%nilyr
      ns = cice%nslyr
      this%nx = nx
      this%ny = ny
      this%ncat = nc
      this%nilyr = nl
      this%nslyr = ns
      allocate (this%sw_down(nx, ny), this%uatm(nx, ny), this%vatm(nx, ny), this%wind(nx, ny), source=0.0_wp)
      allocate (this%potT(nx, ny), this%Qa(nx, ny), this%rhoa(nx, ny), this%flw(nx, ny), source=0.0_wp)
      allocate (this%fsnow(nx, ny), this%frain(nx, ny), this%uvel(nx, ny), this%vvel(nx, ny), source=0.0_wp)
      allocate (this%Tf(nx, ny), this%sst(nx, ny), this%frzmlt(nx, ny), this%bank(nx, ny), source=0.0_wp)
      allocate (this%tau_x(nx, ny), this%tau_y(nx, ny), this%bank_ocean(nx, ny), source=0.0_wp)
      allocate (this%wet(nx, ny), source=.false.)
      this%wet_set = .false.
      allocate (this%strair_x_T(nx, ny), this%strair_y_T(nx, ny), source=0.0_wp)
      allocate (this%strairxn(nx, ny, nc), this%strairyn(nx, ny, nc), source=0.0_wp)
      allocate (this%fbot(nx, ny), this%rside(nx, ny), this%wlat(nx, ny), source=0.0_wp)
      allocate (this%ierr_n(nx, ny, nc), source=0_int32)
      allocate (this%Tsfc_c(nx, ny), source=0.0_wp)
      allocate (this%qice_c(nx, ny, nl), this%sice_c(nx, ny, nl), source=0.0_wp)
      allocate (this%qsno_c(nx, ny, ns), source=0.0_wp)
      allocate (this%rdg_conv(nx, ny), this%rdg_shear(nx, ny), source=0.0_wp)
      allocate (this%niter(nx, ny), source=0_int32)
      allocate (this%salinz(nl), source=0.0_wp)
      allocate (sp(nl + 1))
      call ice_salinity_profile(nl, cice%p%c%saltmax, cice%p%c%min_salin, sp)
      this%salinz(1:nl) = sp(1:nl)
      this%is_init = .true.
      this%mapped = .false.
   end subroutine work_init

   subroutine work_destroy(this)
      class(ice_step_work_t), intent(inout) :: this

      if (this%mapped) call this%exit_data()
      this%is_init = .false.
      if (allocated(this%sw_down)) deallocate (this%sw_down, this%uatm, this%vatm, this%wind)
      if (allocated(this%potT)) deallocate (this%potT, this%Qa, this%rhoa, this%flw)
      if (allocated(this%fsnow)) deallocate (this%fsnow, this%frain, this%uvel, this%vvel)
      if (allocated(this%Tf)) deallocate (this%Tf, this%sst, this%frzmlt, this%bank)
      if (allocated(this%tau_x)) deallocate (this%tau_x, this%tau_y, this%bank_ocean)
      if (allocated(this%wet)) deallocate (this%wet)
      if (allocated(this%strair_x_T)) deallocate (this%strair_x_T, this%strair_y_T)
      if (allocated(this%strairxn)) deallocate (this%strairxn, this%strairyn)
      if (allocated(this%fbot)) deallocate (this%fbot, this%rside, this%wlat)
      if (allocated(this%ierr_n)) deallocate (this%ierr_n)
      if (allocated(this%Tsfc_c)) deallocate (this%Tsfc_c)
      if (allocated(this%qice_c)) deallocate (this%qice_c, this%sice_c)
      if (allocated(this%qsno_c)) deallocate (this%qsno_c)
      if (allocated(this%rdg_conv)) deallocate (this%rdg_conv, this%rdg_shear)
      if (allocated(this%niter)) deallocate (this%niter)
      if (allocated(this%salinz)) deallocate (this%salinz)
      this%nx = 0
      this%ny = 0
   end subroutine work_destroy

   pure function step_work_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the slot (one `arr_bytes` term per array).
      class(ice_step_work_t), intent(in) :: this
      integer(int64) :: nbytes

      nbytes = 0_int64
      if (.not. this%is_init) return
      nbytes = arr_bytes(this%sw_down) &
               + arr_bytes(this%uatm) &
               + arr_bytes(this%vatm) &
               + arr_bytes(this%wind) &
               + arr_bytes(this%potT) &
               + arr_bytes(this%Qa) &
               + arr_bytes(this%rhoa) &
               + arr_bytes(this%flw) &
               + arr_bytes(this%fsnow) &
               + arr_bytes(this%frain) &
               + arr_bytes(this%uvel) &
               + arr_bytes(this%vvel) &
               + arr_bytes(this%Tf) &
               + arr_bytes(this%sst) &
               + arr_bytes(this%frzmlt) &
               + arr_bytes(this%bank) &
               + arr_bytes(this%tau_x) &
               + arr_bytes(this%bank_ocean) &
               + arr_bytes(this%tau_y) &
               + arr_bytes(this%wet) &
               + arr_bytes(this%salinz) &
               + arr_bytes(this%strair_x_T) &
               + arr_bytes(this%strair_y_T) &
               + arr_bytes(this%strairxn) &
               + arr_bytes(this%strairyn) &
               + arr_bytes(this%fbot) &
               + arr_bytes(this%rside) &
               + arr_bytes(this%wlat) &
               + arr_bytes(this%ierr_n) &
               + arr_bytes(this%Tsfc_c) &
               + arr_bytes(this%qice_c) &
               + arr_bytes(this%sice_c) &
               + arr_bytes(this%qsno_c) &
               + arr_bytes(this%rdg_conv) &
               + arr_bytes(this%rdg_shear) &
               + arr_bytes(this%niter)
   end function step_work_bytes

   subroutine work_set_wet(this, mask)
      !! The setup's land mask (`.true.` = ocean) into `wet`; pushes it to the device when mapped.
      class(ice_step_work_t), intent(inout) :: this
      logical, intent(in) :: mask(:, :)

      if (.not. this%is_init) return
      if (size(mask, 1) /= this%nx .or. size(mask, 2) /= this%ny) return
      this%wet = mask
      this%wet_set = .true.
      if (this%mapped) then
         !$acc update device(this%wet)
      end if
   end subroutine work_set_wet

   subroutine work_register_restart(this, reg)
      !! Register the persistent ice-ocean stress (`ice_strocnTx`, `ice_strocnTy`), both optional.
      !! Full local arrays: `ng = 0`, extents `(nx, ny)`.  The caller pulls them to the host before a
      !! write and pushes after a read (like every registry entry).
      class(ice_step_work_t), target, intent(in) :: this
      type(restart_registry_t), intent(inout) :: reg

      if (.not. this%is_init) return
      call reg%register_2d("ice_strocnTx", this%tau_x, 0, this%nx, this%ny, optional=.true.)
      call reg%register_2d("ice_strocnTy", this%tau_y, 0, this%nx, this%ny, optional=.true.)
   end subroutine work_register_restart

   subroutine work_enter_data(this)
      !! `copyin` every COMPONENT array (never the aggregate type); the inputs the
      !! hooks fill are set on the host first and pushed by the hooks' own updates.
      class(ice_step_work_t), intent(inout) :: this

      if (.not. this%is_init .or. this%mapped) return
      !$acc enter data copyin(this%sw_down, this%uatm, this%vatm, this%wind, this%potT, this%Qa, this%rhoa, this%flw)
      !$acc enter data copyin(this%fsnow, this%frain, this%uvel, this%vvel, this%Tf, this%sst, this%frzmlt, this%bank)
      !$acc enter data copyin(this%tau_x, this%tau_y, this%bank_ocean, this%wet, this%salinz)
      !$acc enter data copyin(this%strair_x_T, this%strair_y_T, this%strairxn, this%strairyn)
      !$acc enter data copyin(this%fbot, this%rside, this%wlat, this%ierr_n)
      !$acc enter data copyin(this%Tsfc_c, this%qice_c, this%sice_c, this%qsno_c)
      !$acc enter data copyin(this%rdg_conv, this%rdg_shear, this%niter)
      this%mapped = .true.
   end subroutine work_enter_data

   subroutine work_exit_data(this)
      class(ice_step_work_t), intent(inout) :: this

      if (.not. this%is_init .or. .not. this%mapped) return
      !$acc exit data delete(this%sw_down, this%uatm, this%vatm, this%wind, this%potT, this%Qa, this%rhoa, this%flw)
      !$acc exit data delete(this%fsnow, this%frain, this%uvel, this%vvel, this%Tf, this%sst, this%frzmlt, this%bank)
      !$acc exit data delete(this%tau_x, this%tau_y, this%bank_ocean, this%wet, this%salinz)
      !$acc exit data delete(this%strair_x_T, this%strair_y_T, this%strairxn, this%strairyn)
      !$acc exit data delete(this%fbot, this%rside, this%wlat, this%ierr_n)
      !$acc exit data delete(this%Tsfc_c, this%qice_c, this%sice_c, this%qsno_c)
      !$acc exit data delete(this%rdg_conv, this%rdg_shear, this%niter)
      this%mapped = .false.
   end subroutine work_exit_data

   ! ------------------------------------------------------------------
   ! the step
   ! ------------------------------------------------------------------

   subroutine refuse(res, status, stage, msg)
      type(ice_step_result_t), intent(inout) :: res
      integer, intent(in) :: status, stage
      character(len=*), intent(in) :: msg

      res%status = status
      res%stage = stage
      res%message = msg
   end subroutine refuse

   subroutine check_envelope(cice, flux, rad, work, opts, dt, plan, res)
      !! Everything that can be refused before a value is written.
      type(ice_cice_state_t), intent(in) :: cice
      type(ice_flux_t), intent(in) :: flux
      type(ice_rad_state_t), intent(in) :: rad
      type(ice_step_work_t), intent(in) :: work
      type(ice_step_opts_t), intent(in) :: opts
      real(wp), intent(in) :: dt
      type(ice_step_plan_t), intent(in) :: plan
      type(ice_step_result_t), intent(inout) :: res
      integer :: k
      character(len=160) :: list

      ! missing required stages: ALL are named
      list = ""
      do k = 1, ICE_NSTAGE
         if (plan%required(k) .and. .not. plan%present(k)) then
            if (res%stage == 0) res%stage = k
            list = trim(list)//" "//trim(ice_stage_name(k))
         end if
      end do
      if (res%stage /= 0) then
         call refuse(res, ICE_STEP_REFUSED_MISSING, res%stage, "required stage(s) missing:"//trim(list))
         return
      end if
      if (.not. opts%allow_test_doubles) then
         do k = 1, ICE_NSTAGE
            if (plan%required(k) .and. plan%is_double(k)) then
               call refuse(res, ICE_STEP_REFUSED_DOUBLE, k, "stage '"//trim(ice_stage_name(k))// &
                           "' is a test double (not accepted outside tests)")
               return
            end if
         end do
      end if
      if (opts%ndtd /= 1) then
         call refuse(res, ICE_STEP_REFUSED_ENVELOPE, 0, "ndtd /= 1 is not supported")
         return
      end if
      if (opts%nranks /= 1) then
         call refuse(res, ICE_STEP_REFUSED_ENVELOPE, 0, "scheme='cice' is single-rank until the multi-rank PRs")
         return
      end if
      if (.not. (ieee_is_finite(dt) .and. dt > 0.0_wp)) then
         call refuse(res, ICE_STEP_REFUSED_ENVELOPE, 0, "the ice step dt must be finite and positive")
         return
      end if
      if (opts%dt_outer > 0.0_wp .and. .not. (abs(dt - opts%dt_outer) <= 1.0e-9_wp*opts%dt_outer)) then
         call refuse(res, ICE_STEP_REFUSED_ENVELOPE, 0, "the ice step must equal the outer step")
         return
      end if
      if (.not. (cice%is_init .and. flux%is_init .and. rad%is_init .and. work%is_init)) then
         call refuse(res, ICE_STEP_REFUSED_ENVELOPE, 0, "a slot is not initialised")
         return
      end if
      if (.not. ice_slots_agree(flux, cice) .or. cice%nghost < 3) then
         call refuse(res, ICE_STEP_REFUSED_ENVELOPE, 0, "slots disagree, or nghost < 3 (the dynamics needs 3)")
         return
      end if
      if (work%nx /= cice%nx .or. work%ny /= cice%ny .or. work%ncat /= cice%ncat .or. work%nilyr /= cice%nilyr &
          .or. work%nslyr /= cice%nslyr .or. rad%nx /= cice%nx .or. rad%ny /= cice%ny .or. rad%ncat /= cice%ncat &
          .or. rad%nilyr /= cice%nilyr) then
         call refuse(res, ICE_STEP_REFUSED_ENVELOPE, 0, "work or radiation slot extents differ from the state")
         return
      end if
      if (.not. (cice%mapped .and. flux%mapped .and. rad%mapped .and. work%mapped)) then
         call refuse(res, ICE_STEP_REFUSED_ENVELOPE, 0, "a slot is not device-mapped (enter_data)")
         return
      end if
      if (.not. work%wet_set) then
         call refuse(res, ICE_STEP_REFUSED_ENVELOPE, 0, "work%wet was never set (set_wet from the land mask)")
         return
      end if
      if (cice%ncat < 1 .or. cice%ncat > ICE_NCAT_MAX .or. cice%nilyr < 1 .or. cice%nilyr > ICE_NILYR_MAX &
          .or. cice%nslyr < 1 .or. cice%nslyr > ICE_NSLYR_MAX) then
         call refuse(res, ICE_STEP_REFUSED_ENVELOPE, 0, "category or layer count outside the compile-time maxima")
         return
      end if
      if (.not. rad%configured) then
         call refuse(res, ICE_STEP_REFUSED_ENVELOPE, ICE_STAGE_RADIATION, "the radiation slot is not configured")
         return
      end if
   end subroutine check_envelope

   subroutine note_exec(res, stage)
      type(ice_step_result_t), intent(inout) :: res
      integer, intent(in) :: stage

      res%n_exec = res%n_exec + 1
      res%executed(res%n_exec) = stage
   end subroutine note_exec

   subroutine run_hook(proc, ctx, stage, res, ok)
      !! Call one hook; a non-zero status stops the step (`ok = .false.`).
      procedure(ice_stage_hook_i) :: proc
      type(ice_step_ctx_t), intent(inout) :: ctx
      integer, intent(in) :: stage
      type(ice_step_result_t), intent(inout) :: res
      logical, intent(out) :: ok
      integer :: st

      st = 0
      call note_exec(res, stage)
      call proc(ctx, st)
      ok = st == 0
      if (.not. ok) then
         res%stage_status = st
         call refuse(res, ICE_STEP_STAGE_FAILED, stage, "stage '"//trim(ice_stage_name(stage))//"' failed")
      end if
   end subroutine run_hook

   subroutine enforce(cice, opts, res, stage, ok)
      !! Reduce the cell error word (one host sync); a fatal bit ends the step.
      type(ice_cice_state_t), intent(inout) :: cice
      type(ice_step_opts_t), intent(in) :: opts
      type(ice_step_result_t), intent(inout) :: res
      integer, intent(in) :: stage
      logical, intent(out) :: ok
      integer :: is, ie, js, je

      is = cice%nghost + 1
      ie = cice%nghost + cice%nx_phys
      js = cice%nghost + 1
      je = cice%nghost + cice%ny_phys
      call ice_err_reduce(cice%ice_err, cice%nx, cice%ny, is, ie, js, je, cice%i_offset_global, &
                          cice%j_offset_global, cice%nx_global, opts%fatal_mask, res%summary)
      res%err_mask = ior(res%err_mask, res%summary%mask)
      ok = .not. ice_err_is_fatal(res%summary, opts%fatal_mask)
      if (.not. ok) then
         call refuse(res, ICE_STEP_FATAL, stage, "fatal ice error bit(s) after stage '"//trim(ice_stage_name(stage))//"'")
      end if
   end subroutine enforce

   subroutine ice_step_cice(cice, flux, rad, work, hooks, opts, dt, thermo_step, istep, res)
      !! One ice step (module docstring).  `dt` [s] must be the outer step; the
      !! stages run in the order of `res%plan`.  On return `res%status` is
      !! `ICE_STEP_OK` or why not; a refusal wrote nothing.
      type(ice_cice_state_t), intent(inout), target :: cice
      type(ice_flux_t), intent(inout), target :: flux
      type(ice_rad_state_t), intent(inout), target :: rad
      type(ice_step_work_t), intent(inout), target :: work
      type(ice_step_hooks_t), intent(in) :: hooks
      type(ice_step_opts_t), intent(in) :: opts
      real(wp), intent(in) :: dt
      logical, intent(in) :: thermo_step
         !! The ocean thermodynamics step (`dyn%is_thermo_step()`).
      integer, intent(in) :: istep
      type(ice_step_result_t), intent(out) :: res
      type(ice_step_ctx_t) :: ctx
      integer :: is, ie, js, je, nx, ny, nc, nl, ns, ierr, istat
      logical :: ok

      res%plan = ice_step_plan_build(hooks, cice%p%atm%highfreq)
      call check_envelope(cice, flux, rad, work, opts, dt, res%plan, res)
      if (res%status /= ICE_STEP_OK) return

      ctx%cice => cice
      ctx%flux => flux
      ctx%rad => rad
      ctx%work => work
      ctx%dt = dt
      ctx%istep = istep
      ctx%thermo_step = thermo_step
      ctx%transport_post = res%plan%transport_post
      nx = cice%nx
      ny = cice%ny
      nc = cice%ncat
      nl = cice%nilyr
      ns = cice%nslyr
      is = cice%nghost + 1
      ie = cice%nghost + cice%nx_phys
      js = cice%nghost + 1
      je = cice%nghost + cice%ny_phys

      call ice_err_clear(cice%ice_err, nx, ny)

      ! ---- prep group: the inputs (work slot only), then the first writes ----
      call run_hook(hooks%ocean_import, ctx, ICE_STAGE_OCEAN_IMPORT, res, ok)
      if (.not. ok) return
      if (res%plan%required(ICE_STAGE_VELOCITY_TO_T)) then
         call run_hook(hooks%velocity_to_T, ctx, ICE_STAGE_VELOCITY_TO_T, res, ok)
         if (.not. ok) return
      end if
      call run_hook(hooks%atm_prepare, ctx, ICE_STAGE_ATM_PREPARE, res, ok)
      if (.not. ok) return
      if (.not. (ieee_is_finite(work%zlvl) .and. ieee_is_finite(work%zlvs) .and. work%zlvl > 0.0_wp &
                 .and. work%zlvs > 0.0_wp)) then
         call refuse(res, ICE_STEP_REFUSED_ENVELOPE, ICE_STAGE_ATM_PREPARE, &
                     "the forcing heights zlvl, zlvs must be finite and positive")
         return
      end if
      call note_exec(res, ICE_STAGE_FLUX_ZERO)
      call ice_flux_zero(flux)
      call note_exec(res, ICE_STAGE_SAVE_INIT)
      call ice_save_init(cice)

      ! ---- thermodynamics group: radiation + therm1 + through-ice shortwave are ONE call ----
      call note_exec(res, ICE_STAGE_RADIATION)
      call ice_step_sw_therm1(rad, cice%p%sw, nx, ny, nc, nl, ns, is, ie, js, je, dt, cice%p, cice%p%sw%sw_redist, &
                              work%zlvl, work%zlvs, work%sw_down, work%uatm, work%vatm, work%wind, work%potT, &
                              work%Qa, work%rhoa, work%flw, work%fsnow, work%uvel, work%vvel, work%Tf, work%sst, &
                              work%frzmlt, work%tau_x, work%tau_y, &
                              cice%aicen_init, cice%aicen, cice%vicen, cice%vsnon, cice%Tsfcn, cice%qicen, &
                              cice%sicen, cice%qsnon, &
                              flux%fresh, flux%fsalt, flux%fhocn, flux%evap, flux%fsens, flux%flat, flux%flwout, &
                              flux%fswabs, flux%meltt, flux%meltb, flux%congel, flux%snoice, &
                              flux%freshn, flux%fsaltn, flux%fhocnn, flux%evapn, flux%fsensn, flux%flatn, &
                              flux%flwoutn, flux%fswabsn, flux%melttn, flux%meltbn, flux%congeln, flux%snoicen, &
                              work%strairxn, work%strairyn, work%strair_x_T, work%strair_y_T, &
                              work%fbot, work%rside, work%wlat, work%ierr_n, cice%ice_err, &
                              flux%fswthru, flux%fswthru_vis, ierr, istat)
      if (ierr /= ICE_RAD_OK) then
         res%stage_status = ierr
         call refuse(res, ICE_STEP_STAGE_FAILED, ICE_STAGE_RADIATION, "stage 'radiation' failed")
         return
      end if
      call note_exec(res, ICE_STAGE_THERM1)
      if (istat /= 0) then
         res%stage_status = istat
         call refuse(res, ICE_STEP_STAGE_FAILED, ICE_STAGE_THERM1, "stage 'therm1' failed")
         return
      end if
      call note_exec(res, ICE_STAGE_THERM2)
      call ice_step_therm2(nx, ny, nc, nl, ns, is, ie, js, je, dt, cice%p, cice%hin_max, work%salinz, work%wet, &
                           work%Tf, work%frzmlt, work%frain, work%rside, work%bank, &
                           cice%aicen_init, cice%vicen_init, cice%aicen, cice%vicen, cice%vsnon, cice%Tsfcn, &
                           cice%qicen, cice%sicen, cice%qsnon, cice%aice, cice%vice, cice%vsno, cice%aice0, &
                           work%Tsfc_c, work%qice_c, work%sice_c, work%qsno_c, &
                           flux%fresh, flux%fsalt, flux%fhocn, flux%meltl, flux%frazil_vol, cice%ice_err)
      call enforce(cice, opts, res, ICE_STAGE_THERM2, ok)
      if (.not. ok) return

      ! ---- dynamics group: ghost ring and ghost aggregates first ----
      call run_hook(hooks%halo_state, ctx, ICE_STAGE_HALO_STATE, res, ok)
      if (.not. ok) return
      call note_exec(res, ICE_STAGE_AGGREGATE_GHOST)
      call ice_update_state(cice, .true.)
      call run_hook(hooks%stress_to_faces, ctx, ICE_STAGE_STRESS_TO_FACES, res, ok)
      if (.not. ok) return
      call run_hook(hooks%dynamics, ctx, ICE_STAGE_DYNAMICS, res, ok)
      if (.not. ok) return
      call enforce(cice, opts, res, ICE_STAGE_DYNAMICS, ok)
      if (.not. ok) return

      ! ---- transport, enforced alone, then the closing group ----
      call run_hook(hooks%transport, ctx, ICE_STAGE_TRANSPORT, res, ok)
      if (.not. ok) return
      call enforce(cice, opts, res, ICE_STAGE_TRANSPORT, ok)
      if (.not. ok) return
      call run_hook(hooks%ridging, ctx, ICE_STAGE_RIDGING, res, ok)
      if (.not. ok) return
      call run_hook(hooks%halo_state_post, ctx, ICE_STAGE_HALO_STATE_POST, res, ok)
      if (.not. ok) return
      call note_exec(res, ICE_STAGE_AGGREGATE)
      call ice_update_state(cice, .true.)
      call enforce(cice, opts, res, ICE_STAGE_AGGREGATE, ok)
      if (.not. ok) return

      ! ---- hand-off ----
      call run_hook(hooks%ocean_handoff, ctx, ICE_STAGE_OCEAN_HANDOFF, res, ok)
   end subroutine ice_step_cice

end module rdb_ice_step_cice
