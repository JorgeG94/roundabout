!! Top-level ocean-regime state.
module rdb_ocean_state
   !! Composes every type the ocean dynamical core needs into a single
   !! object the driver builds and passes around for `sim_type='ocean'`
   !! runs.  Each component owns its allocations + bound init/destroy.
   !! GPU mapping is dispatched via the free `ocean_state_enter_data` /
   !! `ocean_state_exit_data` orchestrator below — adding a slot
   !! requires wiring it into both.
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_constants, only: wp, LAND_DEPTH_THRESHOLD, GRAVITY, H_VANISHED, H_DIV_EPS, &
                            DEG2RAD, TWO_PI
   use rdb_grid, only: hgrid_t
   use pic_logger, only: logger => global_logger
   use rdb_error_ring, only: fail
   use pic_strings, only: to_string
   use rdb_comm_env, only: comm_env_rank
   use rdb_barotropic_state, only: barotropic_state_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_pseudo_salt, only: ocean_pseudo_salt_register, ocean_pseudo_salt_seed
   use rdb_ocean_dyn, only: ocean_dyn_t
   use rdb_continuity, only: continuity_t
   use rdb_coriolis_adv, only: coriolis_adv_t
   use rdb_ocean_pressure_force, only: ocean_pressure_force_t, parse_opgf_variant, &
                                       OPGF_VARIANT_MONT, OPGF_VARIANT_FV_LITE, &
                                       OPGF_VARIANT_FV_WRIGHT, OPGF_VARIANT_FV_MOM6
   use rdb_ocean_vmix, only: ocean_vmix_t
   use rdb_ocean_epbl, only: ocean_epbl_t
   use rdb_ocean_wave_speed, only: ocean_wave_speed_t
   use rdb_ocean_kappa_shear, only: ocean_kappa_shear_t
   use rdb_ocean_tidal_mixing, only: ocean_tidal_mixing_t
   use rdb_ocean_lateral_mix, only: ocean_lateral_mix_t
   use rdb_ocean_isopycnal_slopes, only: ocean_slopes_t
   use rdb_ocean_mle, only: ocean_mle_t
   use rdb_ocean_gm, only: ocean_gm_t
   use rdb_ocean_redi, only: ocean_redi_t
   use rdb_ocean_varmix, only: ocean_varmix_t
   use rdb_ocean_meke, only: ocean_meke_t
   use rdb_ocean_horizontal_viscosity, only: ocean_horizontal_viscosity_t
   use rdb_ocean_bottom_drag, only: ocean_bottom_drag_t
   use rdb_ocean_top_drag, only: ocean_top_drag_t
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t
   use rdb_ocean_cavity_flux, only: ocean_cavity_flux_t
   use rdb_ocean_vertical_advection, only: ocean_vertical_advection_t
   use rdb_ocean_hdiff_tracer, only: ocean_hdiff_tracer_t
   use rdb_ocean_vdiff, only: ocean_vdiff_t
   use rdb_eos, only: eos_t
   use rdb_ocean_tides, only: ocean_tides_t
   use rdb_ocean_p_surf, only: ocean_p_surf_t
   use rdb_ice_state, only: ocean_sea_ice_t
   use rdb_ocean_obc, only: ocean_obc_t
   use rdb_ocean_boundary_types, only: ocean_bc_state_t, ocean_bc_state_init, &
                                       ocean_bc_state_destroy, &
                                       ocean_bc_state_enter_data, &
                                       ocean_bc_state_exit_data
   use rdb_ocean_sponge, only: ocean_sponge_t
   use rdb_ocean_diag, only: ocean_diag_t
   use rdb_ocean_restart, only: ocean_restart_t, restart_registry_t
#ifndef RDB_NO_NETCDF
   use rdb_ocean_restart_io, only: ocean_restart_write_local, ocean_restart_read_local, &
                                   ocean_restart_metadata_t
   use rdb_ocean_data_input, only: ocean_data_input_t
   use rdb_ocean_data_forcing, only: ocean_data_forcing_t
#endif
   use rdb_decomp, only: decomp_t
   use rdb_ocean_vcoord, only: ocean_vcoord_t, parse_ocean_vcoord_type
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_ocean_cavity, only: parse_cavity_draft_config, parse_cavity_draft_source, &
                               CAVITY_DRAFT_NONE, CAVITY_DRAFT_FLAT, CAVITY_DRAFT_LINEAR, &
                               CAVITY_SOURCE_DRAFT, CAVITY_SOURCE_THICKNESS, &
                               set_draft_flat, set_draft_linear, &
                               cavity_water_column_impl, cavity_apply_land_exclusion, &
                               cavity_count_grounded, cavity_fill_cover_frac, &
                               cavity_draft_is_finite_nonneg, CAVITY_BOUND_INF
   use rdb_config, only: config_t, ice_hlim_count
   use rdb_profiler, only: profiler_start, profiler_stop
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_IC_SEED, &
                               OCEAN_STATUS_ERR_IO, OCEAN_STATUS_ERR_RESTART_SCHEMA, &
                               OCEAN_STATUS_ERR_RESTART_DECOMP, OCEAN_STATUS_ERR_RESTART_GRID, &
                               OCEAN_STATUS_ERR_BAD_SHAPE
   use rdb_ocean_bathymetry_inject, only: bathymetry_normalise_sign, bathymetry_fill_ghosts_array
#ifndef RDB_NO_NETCDF
   use rdb_bathymetry, only: load_bathymetry_into_array
   use rdb_ocean_z_init, only: seed_ts_from_zfile, seed_ts_linear_z
#endif
   implicit none
   private

   public :: ocean_state_t
   public :: ocean_state_enter_data
   public :: ocean_state_exit_data
   public :: ocean_state_seed_from_cfg
   public :: ocean_state_seed_land_cells
   public :: ocean_state_build_restart_registry
   public :: ocean_state_restart_write
   public :: ocean_state_restart_read
   public :: seed_eady_ic
   public :: seed_geostrophic_adjustment_ic
   public :: seed_baroclinic_jet_ic
   public :: set_bathymetry_spoon
   public :: set_bathymetry_seamount
   public :: seed_h_layer_uniform_impl
   public :: seed_h_layer_uniform_z_impl
   public :: topo_length_to_grid_units
   public :: set_bathymetry_neverworld2
   public :: set_bathymetry_island
   public :: apply_layer_rho_init
   public :: ocean_linear_layer_density
   public :: pgf_nonoverlap_gate_on

   type :: ocean_state_t
      logical :: is_init = .false.
         !! True once every nested slot has been initialised + GPU
         !! attached.  Set only after the last child `%init` returns.
         !! Prefer this to walking the children with `allocated(...)`.

      ! ---- Prognostic state ----
      type(barotropic_state_t) :: barotropic
         !! Depth-integrated C-grid prognostic state.
      type(multilayer_state_t) :: multilayer
         !! Per-layer C-grid prognostic state.
      type(ocean_metrics_t) :: metrics
         !! Orthogonal curvilinear horizontal metrics (lengths, areas,
         !! inverses, geography, hvisc ratio bundle).  Filled by
         !! `configure_ocean_metrics`.

      ! ---- Kernel + driver state ----
      type(ocean_dyn_t) :: dyn
         !! Split-explicit RK2 driver state.
      type(continuity_t) :: continuity
         !! Continuity-PPM kernel state.
      type(coriolis_adv_t) :: coriolis_adv
         !! PV-conserving Coriolis+advection kernel state.
      type(ocean_pressure_force_t) :: pressure_force
         !! FV pressure-force kernel state.

      ! ---- Physics parameterisations ----
      type(ocean_vmix_t) :: vmix
         !! Vertical mixing (KPP + interior closure).
      type(ocean_epbl_t) :: epbl
      type(ocean_wave_speed_t) :: wavespeed
         !! First-baroclinic wave speed + Rossby deformation radius.
         !! Diagnostic, default off.
      type(ocean_tidal_mixing_t) :: vmix_tidal
         !! St-Laurent/Simmons internal-tide interior diapycnal mixing.
         !! INTERIOR closure: its Kd adds to the surface PBL +
         !! PP81/background/kappa-shear via the additive merge each stage;
         !! `kd_int` refreshes at thermo cadence.  Default off.
      type(ocean_kappa_shear_t) :: kshear
         !! Energetics-based PBL.  Mutually exclusive with the KPP
         !! overlay; folds its interface diffusivity into `vmix%kv` /
         !! `vmix%kt` each stage.
      type(ocean_lateral_mix_t) :: lateral_mix
         !! Lateral mixing closure (Leith / Smagorinsky / biharmonic);
         !! supplies face viscosity coefficients.
      type(ocean_slopes_t) :: slopes
         !! Isopycnal (neutral) slope diagnostics (Griffies 1998) at
         !! C-grid interfaces — gates the mesoscale-eddy params.
         !! Diagnostic, default off; refreshed at thermo cadence.
      type(ocean_mle_t) :: mle
         !! Fox-Kemper mixed-layer-eddy restratification.  Injects
         !! ML-confined overturning transports into the continuity mass
         !! fluxes before the divergence; reads `epbl%mld`.  Default off.
      type(ocean_gm_t) :: gm
         !! Gent-McWilliams thickness diffusion.  Injects eddy bolus
         !! thickness transports into the continuity mass fluxes before
         !! the divergence; reads `slopes%slope_x/y`.  Default off.
      type(ocean_redi_t) :: redi
         !! Redi continuous neutral (along-isopycnal) tracer diffusion.
         !! Two-phase: `redi_calc_coeffs` fills neutral-surface coeffs at
         !! thermo cadence, `redi_apply_flux` adds the rotated tracer flux
         !! after `hdiff_tracer`.  Recomputes its own interface density
         !! derivs (does NOT read `slopes`).  Default off.
      type(ocean_varmix_t) :: varmix
         !! Spatially-varying GM/Redi coefficient fields.  Produces the
         !! pre-CFL base `khth_u/v` (+ `khtr_u/v`) face fields GM consumes
         !! as its optional external base.  Reads slopes + wavespeed.
         !! Default off ⇒ GM uses the scalar khth.
      type(ocean_meke_t) :: meke
         !! Prognostic mesoscale eddy kinetic energy.  Evolves a 2D
         !! eddy-energy field sourced by `gm%gm_src`, damped by implicit
         !! bottom drag, transported laterally, fed back (geom-mean) into
         !! `varmix%khth_u/v` (+ khtr).  Requires GM.  Default off ⇒ no-op;
         !! `khth_fac=khtr_fac=0` ⇒ feedback inert even when enabled.
      type(ocean_horizontal_viscosity_t) :: hvisc
         !! Constant-coefficient Laplacian horizontal-viscosity kernel.
         !! Reads scalar `nu_h`.
      type(ocean_bottom_drag_t) :: bdrag
         !! Bottom-drag kernel (linear Rayleigh or quadratic log-layer).
         !! Acts on the k=1 layer only.
      type(ocean_top_drag_t) :: tdrag
         !! Ice-shelf TOP-drag kernel (`&ocean_tdrag_nml`, default off).
         !! The mirror of `bdrag` at `k = nz`, masked by ice cover on
         !! faces.  Off ⇒ placeholder arrays, no kernel.
      type(ocean_surface_stress_t) :: surface_stress
         !! Surface wind-stress kernel.  Acts on the k=nz layer only.
      type(ocean_cavity_flux_t) :: cavity_flux
         !! Ice-shelf basal-melt slot (`&ocean_cavity_melt_nml`, default
         !! off).  Holds the sampled far field, `u*`, the interface
         !! state, the melt mass flux and the per-column solver status;
         !! its driver fills `surface_flux%heat_cavity`/`salt_cavity`.
         !! Gated: off ⇒ `(1,1)` placeholders, no kernel, bit-identical.
      type(ocean_surface_flux_t) :: surface_flux
         !! Surface heat + salt flux slot.  2D `Q_heat(:,:)` / `Q_salt(:,:)`
         !! (W/m^2 and kg/m^2/s, positive downward / salinifying).  Seeded
         !! at configure from `&ocean_thermo_nml q_heat / q_salt`;
         !! data-override / restoring may overwrite the fields after.
         !! Device-mapped via the orchestrator.
      type(ocean_vertical_advection_t) :: vert_advect
         !! Vertical advection: diagnoses `w_interface` from horizontal
         !! continuity and advects tracers vertically with a matching
         !! `h_layer` update (Eulerian z — horizontal and vertical h
         !! updates cancel).
      type(ocean_hdiff_tracer_t) :: hdiff_tracer
         !! Constant-coefficient Laplacian horizontal-tracer-diffusion
         !! kernel.  Tracer-side counterpart of `hvisc`.
      type(ocean_vdiff_t) :: vdiff
         !! Backward-Euler vertical-diffusion solver (Thomas tridiagonal
         !! per column) for momentum + tracers, using the eddy
         !! coefficients in `vmix%kv` / `vmix%kt` / `vmix%ks`.
      type(eos_t) :: eos
         !! Equation of state (Wright by default).

      ! ---- Forcing ----
      type(ocean_tides_t) :: tides
         !! Equilibrium + SAL + internal-tide drag.
      type(ocean_p_surf_t) :: p_surf
         !! Atmospheric surface-pressure loading / inverse barometer
         !! (PR-17).  Folds `eta_ib = -p_surf/(rho0 g)` into the same
         !! `eta_forcing` barotropic seam the tide composes through.

      ! ---- Sea ice ----
      type(ocean_sea_ice_t) :: ice
         !! Gated sea-ice slot (SIS2 port, `PLAN_SEA_ICE.md`).  PR 0
         !! scaffold: knobs + lifecycle only, no arrays and no physics.
         !! `&ocean_ice_nml enable` default off ⇒ never initialised /
         !! mapped / stepped ⇒ byte-identical.

      ! ---- Boundary ----
      type(ocean_obc_t) :: obc
         !! Open-boundary nest from parent run.
      type(ocean_bc_state_t) :: bc
         !! Per-edge OBC config + per-step boundary data.  Defaults to
         !! all-OBC_WALL (closed wall) until config opens an edge or a
         !! data source populates the `data_*` buffers.
      type(ocean_sponge_t) :: sponge
         !! Map-driven sponge (PR-23): per-cell `idamp_h/u/v` + 3-D
         !! reference state.  `enable` (default `.false.`) gates BOTH the
         !! allocation (gated `init` call below, mirrors epbl/kshear) and
         !! the device mapping (gated `enter_data`/`exit_data` call in the
         !! orchestrator) — a disabled sponge costs nothing.  Default off
         !! ⇒ the legacy `bc%<edge>%sponge_*` band kernels
         !! (`rdb_ocean_sponge::ocean_sponge_apply{,_tracers}`) run
         !! unchanged ⇒ bit-identical.

      ! ---- Vertical coordinate ----
      type(ocean_vcoord_t) :: vcoord
         !! Vertical-coordinate config + per-step ALE remap state.
         !! Defaults to `VCOORD_EULERIAN_Z` (`h_layer` fixed).

      ! ---- I/O ----
      type(ocean_diag_t) :: diag
         !! Diagnostics manager (registry + remap + I/O server hand-off).
#ifndef RDB_NO_NETCDF
      type(ocean_data_input_t) :: data_input
         !! The shared time-varying NetCDF input reader (PR-14).  NetCDF-
         !! only (mirrors `zinit` — the slot compiles out entirely with
         !! `RDB_ENABLE_NETCDF=OFF`, and no consumer path exists that
         !! would need it live in that build).
      type(ocean_data_forcing_t) :: data_forcing
         !! File-backed surface forcing (PR-15) — the `(file, variable)
         !! -> forcing slot` binding over `data_input`.  Holds only
         !! registration ids + logicals, all read host-side, so unlike
         !! every other slot here it deliberately adds NO term to
         !! `ocean_state_enter_data`: there is nothing to map.  The
         !! arrays it writes are mapped by `surface_stress` /
         !! `surface_flux`.
#endif
      type(ocean_restart_t) :: restart
         !! Restart / checkpoint manager.

      ! ---- Top-level scalars ----
      logical :: use_multilayer = .true.
         !! Whether the 3D multilayer state is active (always true on
         !! the ocean path; flag exposed for symmetry with the coastal
         !! state and so 2D-only debug runs can disable it later).
      logical :: use_nonhydrostatic = .false.
         !! NH-on-ocean toggle.  Phase 5f.
      logical :: enable_ideal_age = .false.
         !! MOM6 `USE_IDEAL_AGE_TRACER`.  Set by the driver from cfg
         !! BEFORE calling `init(grid)` so the multilayer state
         !! registers the age tracer at index 3.  Per-step aging +
         !! surface reset is then called by the dyn step driver.
      logical :: enable_pseudo_salt = .false.
         !! `&ocean_tracers_nml enable_pseudo_salt`.  Set by the driver
         !! from cfg BEFORE calling `init(grid)` — pseudo-salt is
         !! registered via `register_passive_tracer` right after
         !! `multilayer%init`, before `ocean_bc_state_init` sizes
         !! `bc%n_tracers` (registration-order contract, see
         !! `rdb_ocean_pseudo_salt`).

      ! ---- Tile policy ----
      ! Single knob controlling the strided hybrid column/plane loop
      ! pattern that every column-local kernel (vmix, continuity,
      ! coriolis_adv, lateral_mix, pressure_force) is written against
      ! — see README "Loop pattern" convention.
      !
      !   column_stride = 1        → CPU-optimal (column-based, 1D
      !                              column workspace, cache-friendly)
      !   column_stride = nx       → GPU-optimal (plane-based, full
      !                              3D parallelism, max occupancy)
      !   1 < column_stride < nx   → tiled hybrid (typically 16–64
      !                              saturates an SM on H100/A100)
      !
      ! Phase 5b will auto-pick based on RDB_PARALLEL_BACKEND:
      ! 1 for `-stdpar=multicore`, 32 for `-stdpar=gpu`.  Reference:
      ! Kommera & Appelhans (NVIDIA, CESM SEWG 2026), "One Codebase
      ! with Good Performance on both CPUs and GPUs, for Column-based
      ! Loop Structures".
      integer :: column_stride = 1
         !! Outer tile size for column-local kernels.  1 = column,
         !! `nx` = plane.  Drives slot workspace allocation shape.
         !! Inform restart-roundtrip + cross-rank assertions whether
         !! bit-for-bit checks are expected to succeed.
   contains
      procedure, non_overridable :: init => ocean_state_init
      procedure, non_overridable :: init_from_config => ocean_state_init_from_config
      procedure, non_overridable :: destroy => ocean_state_destroy
      procedure, non_overridable :: bytes => ocean_state_bytes
   end type ocean_state_t

contains

   subroutine ocean_state_init(this, grid)
      !! Construct the ocean god state.  Each slot's `init` allocates its
      !! own arrays; the DEFAULT-OFF closures (EPBL, kappa-shear, tidal
      !! mixing, GM, Redi, MLE/Fox-Kemper, VarMix, MEKE, isopycnal slopes)
      !! are gated on their `enable` flag so a plain run does not pay their
      !! multi-GB footprint.  This requires the enable flags to be set
      !! BEFORE `init` — `ocean_state_init_from_config` hoists them above
      !! its `call this%init(grid)` for exactly this reason.  The gated
      !! closures' runtime kernels already early-return on `.not. enable`
      !! (and `enter_data`/`exit_data` are gated in the parent walk), so a
      !! gated-off slot is never touched with unallocated arrays.
      class(ocean_state_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid

      logical :: need_slopes

      ! Isopycnal slopes feed GM (and the standalone slope diffusivity);
      ! Redi and MLE compute their own coefficients and do NOT read slopes.
      need_slopes = this%slopes%enable .or. this%gm%enable

      call this%barotropic%init(grid)
      call this%metrics%init(grid)
      if (this%use_multilayer) then
         call this%multilayer%init(grid, with_ideal_age=this%enable_ideal_age)
         ! Pseudo-salt registration MUST happen here — after multilayer%init
         ! (so the registry exists), before enter_data (registry_locked
         ! guard) and, critically, before ocean_bc_state_init below sizes
         ! bc%n_tracers (§2.5 of the plan: registering after that call
         ! silently skips the new tracer at OBC-clamped edges).
         if (this%enable_pseudo_salt) call ocean_pseudo_salt_register(this%multilayer, grid)
         call this%dyn%init(grid, nz_ml=this%multilayer%nz_ml)
      else
         call this%dyn%init(grid)
      end if
      ! the ocean is always multilayer
      if (this%use_multilayer) then
         call this%continuity%init(grid, nz_ml=this%multilayer%nz_ml)
         call this%coriolis_adv%init(grid, nz_ml=this%multilayer%nz_ml)
         call this%pressure_force%init(grid, nz_ml=this%multilayer%nz_ml)
         call this%hvisc%init(grid, nz_ml=this%multilayer%nz_ml)
         call this%bdrag%init(grid, nz_ml=this%multilayer%nz_ml)
         call this%tdrag%init(grid, nz_ml=this%multilayer%nz_ml)
         call this%surface_stress%init(grid, nz_ml=this%multilayer%nz_ml)
         call this%surface_flux%init(grid)
         call this%cavity_flux%init(grid)
         call this%vert_advect%init(grid, nz_ml=this%multilayer%nz_ml)
         call this%hdiff_tracer%init(grid, nz_ml=this%multilayer%nz_ml)
         call this%vdiff%init(grid, nz_ml=this%multilayer%nz_ml)
      else
         call this%continuity%init(grid)
         call this%coriolis_adv%init(grid)
         call this%pressure_force%init(grid)
         call this%hvisc%init(grid)
         call this%bdrag%init(grid)
         call this%tdrag%init(grid)
         call this%surface_stress%init(grid)
         call this%surface_flux%init(grid)
         call this%cavity_flux%init(grid)
         call this%vert_advect%init(grid)
         call this%hdiff_tracer%init(grid)
         call this%vdiff%init(grid)
      end if
      ! Always-on vertical mixing (PP81/KPP base) + horizontal viscosity
      ! (harmonic/biharmonic) back the production envelope, so they are
      ! unconditional.  The eddy/alternative closures below are gated.
      if (this%use_multilayer) then
         call this%vmix%init(grid, nz_ml=this%multilayer%nz_ml)
         call this%lateral_mix%init(grid, nz_ml=this%multilayer%nz_ml)
         if (this%epbl%enable) call this%epbl%init(grid, nz_ml=this%multilayer%nz_ml)
         if (this%kshear%enable) call this%kshear%init(grid, nz_ml=this%multilayer%nz_ml)
         if (this%vmix_tidal%enable) call this%vmix_tidal%init(grid, nz_ml=this%multilayer%nz_ml)
         if (need_slopes) call this%slopes%init(grid, nz_ml=this%multilayer%nz_ml)
         if (this%mle%enable) call this%mle%init(grid, nz_ml=this%multilayer%nz_ml)
         if (this%gm%enable) call this%gm%init(grid, nz_ml=this%multilayer%nz_ml)
         if (this%redi%enable) call this%redi%init(grid, nz_ml=this%multilayer%nz_ml)
         if (this%varmix%enable) call this%varmix%init(grid, nz_ml=this%multilayer%nz_ml)
         if (this%meke%enable) call this%meke%init(grid, nz_ml=this%multilayer%nz_ml)
      else
         call this%vmix%init(grid)
         call this%lateral_mix%init(grid)
         if (this%epbl%enable) call this%epbl%init(grid)
         if (this%kshear%enable) call this%kshear%init(grid)
         if (this%vmix_tidal%enable) call this%vmix_tidal%init(grid)
         if (need_slopes) call this%slopes%init(grid)
         if (this%mle%enable) call this%mle%init(grid)
         if (this%gm%enable) call this%gm%init(grid)
         if (this%redi%enable) call this%redi%init(grid)
         if (this%varmix%enable) call this%varmix%init(grid)
         if (this%meke%enable) call this%meke%init(grid)
      end if
      call this%wavespeed%init(grid)
      call this%eos%init(grid)
      if (this%use_multilayer) then
         call this%vcoord%init(grid, nz_ml=this%multilayer%nz_ml)
      else
         call this%vcoord%init(grid)
      end if
      call this%tides%init(grid)
      call this%p_surf%init(grid)
      if (this%ice%enable) call this%ice%init(grid)
      call this%obc%init(grid)
      if (this%use_multilayer) then
         call ocean_bc_state_init(this%bc, grid, nz_ml=this%multilayer%nz_ml, &
                                  n_tracers=size(this%multilayer%tracers))
         ! PR-23: gated call (mirrors epbl/kshear/...) — `ocean_sponge_init`
         ! only allocates the (potentially large) idamp/ref_tracer arrays
         ! when `enable`, so a disabled sponge (default) never allocates.
         if (this%sponge%enable) then
            call this%sponge%init(grid, nz_ml=this%multilayer%nz_ml, &
                                  n_tracers=size(this%multilayer%tracers))
         end if
      else
         call ocean_bc_state_init(this%bc, grid, nz_ml=1)
      end if
      call this%diag%init(grid)
#ifndef RDB_NO_NETCDF
      call this%data_input%init()
#endif
      call this%restart%init()
      this%is_init = .true.
   end subroutine ocean_state_init

   subroutine ocean_state_init_from_config(this, cfg, grid)
      !! Seed the cfg-derived scalars the slot inits read up front (layer
      !! count, ideal-age toggle), allocate via init(grid), then override the
      !! linear-EOS params (after eos%init has set its defaults).  Carries the
      !! ocean branch that state_init_from_config held before the coastal /
      !! ocean state split — the order (nz_ml + ideal_age before init, eos
      !! after) is load-bearing and matches the pre-split behaviour.
      class(ocean_state_t), intent(inout) :: this
      type(config_t), intent(in) :: cfg
      type(hgrid_t), intent(in) :: grid

      this%multilayer%nz_ml = cfg%nz_layers
      this%enable_ideal_age = cfg%ocean%tracers%enable_ideal_age
      this%enable_pseudo_salt = cfg%ocean%tracers%enable_pseudo_salt

      ! Enable flags govern CONDITIONAL ALLOCATION (ocean_state_init gates
      ! each default-off closure's array allocation on its `enable`), so
      ! they must be latched BEFORE init(grid).  The full per-closure knob
      ! copy still runs after init (below) and re-sets these idempotently;
      ! MEKE and Fox-Kemper/MLE otherwise get their enable only in the
      ! later configure_ocean_* pass, so they are seeded from cfg here too.
      ! Windowed tracer advection (Phase 2/6b): its 13-array workspace
      ! (~3.7 GiB at 600²x100) is only consumed at ratio > 1, so the
      ! allocation gate latches here like the closures below.
      this%continuity%windowed_advection = cfg%ocean%vmix%dt_tracer_advect_ratio > 1
      ! Pressure-force scratch: 11 of the 16 buffers are FV_MOM6- or
      ! in-layer-reconstruction-only (~3.5 GB at 1000x800x50) and are
      ! unreachable on the shipped gprime / montgomery / fv_lite envelope.
      ! `variant` + `reconstruct_for_pressure` therefore latch HERE, ahead
      ! of init(grid), and `scratch_gated` turns the per-variant gate on;
      ! `configure_ocean_pgf` re-sets both idempotently later.  Direct
      ! `pgf%init(...)` call sites leave `scratch_gated = .false.` and keep
      ! the allocate-everything behaviour.
      this%pressure_force%variant = parse_opgf_variant(cfg%ocean%pgf%form)
      this%pressure_force%reconstruct_for_pressure = cfg%ocean%pgf%reconstruct_for_pressure
      ! Grounded-layer PGF gate: for FV_MOM6 the `z_centre` buffer exists only
      ! to feed this gate, so the flag has to be known BEFORE the allocation
      ! gate runs.  It is a pure config-time decision, so latch it here next to
      ! `variant`; `configure_ocean_pgf` recomputes it (and error-stops if the
      ! latch drifted) for the FV_LITE / FV_WRIGHT paths too.
      this%pressure_force%skip_nonoverlap = &
         pgf_nonoverlap_gate_on(cfg, this%pressure_force%variant)
      this%pressure_force%scratch_gated = .true.
      ! Ice-shelf cavity geometry: `use_cavity` gates the allocation of
      ! the three static metrics fields (`z_draft`, `cover_frac`,
      ! `p_ice_ref`), so it latches HERE — and it has to be earlier than
      ! most: the draft is filled inside `ocean_state_seed_from_cfg`,
      ! which is the first thing the engine does after this, long before
      ! any `configure_ocean_*` pass.  Off ⇒ three `(1,1)` placeholders.
      this%metrics%use_cavity = cfg%ocean%cavity_dyn%enable
      ! Ice-shelf basal melt: the melt slot's own gate, latched here for
      ! the same reason (its `init` sizes eleven 2-D arrays off it).
      this%cavity_flux%enable = cfg%ocean%cavity_melt%enable
      ! Ice-shelf top drag: same gate discipline — its `init` sizes the
      ! two face scratch buffers and five 2-D fields off `enable`.
      this%tdrag%enable = cfg%ocean%tdrag%enable
      this%epbl%enable = cfg%ocean%epbl%enable
      this%kshear%enable = cfg%ocean%kshear%enable
      this%vmix_tidal%enable = cfg%ocean%tidal_mixing%enable
      this%slopes%enable = cfg%ocean%slopes%enable
      this%gm%enable = cfg%ocean%gm%enable
      this%redi%enable = cfg%ocean%redi%enable
      this%varmix%enable = cfg%ocean%varmix%enable
      this%meke%enable = cfg%ocean%meke%enable
      this%mle%enable = cfg%ocean%foxkemper%enable
      ! Sea ice (PR 0 scaffold): enable gates the slot's init +
      ! enter_data/exit_data; ncat/nk_ice size the per-category arrays
      ! once PR 3+ allocates them, so all three latch before init like
      ! the closures above.
      this%ice%enable = cfg%ocean%ice%enable
      this%ice%ncat = cfg%ocean%ice%ncat
      this%ice%nk_ice = cfg%ocean%ice%nk_ice
      ! PR 4b: transport latches before init too (gates the workspace
      ! allocation, memory Rule 2 — see rdb_ice_state%init).
      this%ice%transport = cfg%ocean%ice%transport
      ! PR 5: EVP dynamics master flag (gates the driver's ice_evp_step
      ! call + the transport sampler skip; no per-slot allocation gate —
      ! str_d/str_t/str_s/tau_a_*/fxoc/fyoc are unconditional whenever the
      ! ice slot is live, same contract as u_ice/v_ice).
      this%ice%dynamics = cfg%ocean%ice%dynamics
      ! PR-23: enable gates the sponge's (potentially large, 4-D
      ! nx*ny*nz*n_tracers) reference-state allocation, so it latches
      ! before init like the closures above (`ocean_sponge_init` is only
      ! called when `this%sponge%enable`, see below).
      this%sponge%enable = cfg%ocean%sponge%enable
      ! PR-58: hlim latches before init (ice_itd_category_bounds runs
      ! inside ice%init). Unallocated => the hardcoded SIS2 default table.
      block
         integer :: n_hlim
         n_hlim = ice_hlim_count(cfg%ocean%ice%hlim)
         if (n_hlim > 0) this%ice%hlim_cfg = cfg%ocean%ice%hlim(1:n_hlim)  ! allocate-on-assign
      end block
      ! PR 26: snowfall's host-side dispatch gate — latches before init
      ! like transport/dynamics above (no per-slot allocation gate;
      ! atm_fprec/snow_part_ocn/fprec_ocn_diag are unconditional whenever
      ! the ice slot is live, same contract as u_ice/v_ice). This is what
      ! makes `snowfall=0` byte-identical BY CONSTRUCTION.
      this%ice%has_snowfall = (cfg%ocean%ice%snowfall /= 0.0_wp)
      ! PR 27: Archimedes snow-ice flood master flag — latches before
      ! init like has_snowfall/transport/dynamics above (no per-slot
      ! allocation gate; snow_to_ice is unconditional whenever the ice
      ! slot is live, same contract as sw_thru). This is what makes
      ! `snow_ice=.false.` byte-identical BY CONSTRUCTION.
      this%ice%snow_ice = cfg%ocean%ice%snow_ice

      call this%init(grid)
#ifndef RDB_NO_NETCDF
      ! `this%init(grid)` above already built a default-sized (16 slot,
      ! quiet) reader — safe to re-run with the real config now, since
      ! nothing can have been registered between the two calls (grid
      ! init happens before any consumer registration, and registration
      ! itself only happens after this returns — see the seam contract
      ! in rdb_ocean_data_input's module docstring).
      call this%data_input%init(cfg%ocean%data)
#endif
      ! Linear-EOS reference state, all five members off `&ocean_ic_nml`.
      ! This runs before `ocean_state_enter_data` AND before the
      ! `configure_ocean_*` pass that copies the whole flat-POD handle
      ! onto the closure slots that carry their own copy
      ! (`vmix%eos`, `epbl%eos`, `kshear%eos`, `tidal_mixing%eos` —
      ! `rdb_ocean_setup.F90`), so every rider and every device kernel
      ! that takes `eos_t` by value sees the configured values.  Defaults
      ! equal the `eos_t` component defaults ⇒ bit-identical.
      this%eos%alpha_T = cfg%ocean%ic%alpha_T
      this%eos%beta_S = cfg%ocean%ic%beta_S
      this%eos%T_ref = cfg%ocean%ic%T_ref
      this%eos%S_ref = cfg%ocean%ic%S_ref
      this%eos%rho0 = cfg%ocean%ic%rho_0
      ! EPBL master switch is read here (before the configure_ocean_*
      ! pass) because diag registration — which gates the MLD_EPBL /
      ! Kd_EPBL variables on it — runs before configure in the driver.
      this%epbl%enable = cfg%ocean%epbl%enable
      this%wavespeed%enable = cfg%ocean%wavespeed%enable
      this%kshear%enable = cfg%ocean%kshear%enable
      this%slopes%enable = cfg%ocean%slopes%enable
      this%slopes%kd_smooth = cfg%ocean%slopes%kd_smooth
      this%slopes%min_dz_for_n2 = cfg%ocean%slopes%min_dz_for_n2
      this%slopes%rho0 = cfg%ocean%ic%rho_0
      this%gm%enable = cfg%ocean%gm%enable
      this%gm%khth = cfg%ocean%gm%khth
      this%gm%khth_max_cfl = cfg%ocean%gm%khth_max_cfl
      this%gm%khth_slope_max = cfg%ocean%gm%khth_slope_max
      this%gm%rho0 = cfg%ocean%ic%rho_0
      this%redi%enable = cfg%ocean%redi%enable
      this%redi%continuous = cfg%ocean%redi%continuous
      this%redi%khtr = cfg%ocean%redi%khtr
      this%varmix%enable = cfg%ocean%varmix%enable
      this%varmix%use_visbeck = cfg%ocean%varmix%use_visbeck
      this%varmix%resoln_scaled_khth = cfg%ocean%varmix%resoln_scaled_khth
      this%varmix%resoln_scaled_khtr = cfg%ocean%varmix%resoln_scaled_khtr
      this%varmix%gill_equatorial_ld = cfg%ocean%varmix%gill_equatorial_ld
      this%varmix%interpolate_res_fn = cfg%ocean%varmix%interpolate_res_fn
      this%varmix%kh_res_fn_power = cfg%ocean%varmix%kh_res_fn_power
      this%varmix%kh_res_scale_coef = cfg%ocean%varmix%kh_res_scale_coef
      this%varmix%khth = cfg%ocean%varmix%khth
      this%varmix%khtr = cfg%ocean%varmix%khtr
      this%varmix%khth_slope_cff = cfg%ocean%varmix%khth_slope_cff
      this%varmix%khtr_slope_cff = cfg%ocean%varmix%khtr_slope_cff
      this%varmix%khth_min = cfg%ocean%varmix%khth_min
      this%varmix%khth_max = cfg%ocean%varmix%khth_max
      this%varmix%khtr_min = cfg%ocean%varmix%khtr_min
      this%varmix%khtr_max = cfg%ocean%varmix%khtr_max
      this%varmix%visbeck_l_scale = cfg%ocean%varmix%visbeck_l_scale
      this%varmix%visbeck_max_slope = cfg%ocean%varmix%visbeck_max_slope

      ! ---- MEKE (capability [5]) ----
      this%meke%enable = cfg%ocean%meke%enable
      this%meke%gmcoeff = cfg%ocean%meke%gmcoeff
      this%meke%frcoeff = cfg%ocean%meke%frcoeff
      ! Turn the hvisc KE-dissipation diagnostic on ONLY when MEKE's
      ! frictional source is active; otherwise the viscous apply skips it
      ! (zero extra work, bit-identical).
      this%hvisc%compute_ke_diss = cfg%ocean%meke%enable .and. &
                                   (cfg%ocean%meke%frcoeff >= 0.0_wp)
      this%meke%bgsrc = cfg%ocean%meke%bgsrc
      this%meke%damping = cfg%ocean%meke%damping
      this%meke%kh = cfg%ocean%meke%kh
      this%meke%k4 = cfg%ocean%meke%k4
      this%meke%khcoeff = cfg%ocean%meke%khcoeff
      this%meke%cd_scale = cfg%ocean%meke%cd_scale
      this%meke%cb = cfg%ocean%meke%cb
      this%meke%ct = cfg%ocean%meke%ct
      this%meke%min_gamma2 = cfg%ocean%meke%min_gamma2
      this%meke%uscale = cfg%ocean%meke%uscale
      this%meke%use_bbl_drag = cfg%ocean%meke%use_bbl_drag
      this%meke%dtscale = cfg%ocean%meke%dtscale
      this%meke%khth_fac = cfg%ocean%meke%khth_fac
      this%meke%khtr_fac = cfg%ocean%meke%khtr_fac
      this%meke%backscatter = cfg%ocean%meke%backscatter
      this%meke%visc_coeff_ku = cfg%ocean%meke%backscatter_visc_coeff_ku
      this%meke%khmeke_fac = cfg%ocean%meke%khmeke_fac
      this%meke%advection_factor = cfg%ocean%meke%advection_factor
      this%meke%alpha_deform = cfg%ocean%meke%alpha_deform
      this%meke%alpha_rhines = cfg%ocean%meke%alpha_rhines
      this%meke%alpha_eady = cfg%ocean%meke%alpha_eady
      this%meke%alpha_frict = cfg%ocean%meke%alpha_frict
      this%meke%alpha_grid = cfg%ocean%meke%alpha_grid
      ! cdrag: prefer the side-drag coefficient if set (> 0), else the
      ! MEKE-specific knob default.
      if (cfg%ocean%bdrag%cdrag_side > 0.0_wp) then
         this%meke%cdrag = cfg%ocean%bdrag%cdrag_side
      else
         this%meke%cdrag = cfg%ocean%meke%cdrag
      end if
      this%vmix_tidal%enable = cfg%ocean%tidal_mixing%enable
      this%tides%enable = cfg%ocean%tides%enable
      this%p_surf%enable = cfg%ocean%psurf%enable
      this%p_surf%in_eos = cfg%ocean%psurf%in_eos
   end subroutine ocean_state_init_from_config

   pure function ocean_state_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the whole ocean god state, summed
      !! from each slot's own `bytes()` (each term is 0 when that slot's
      !! arrays are unallocated, so the DEFAULT-OFF closures gated in
      !! `ocean_state_init` contribute nothing — the total tracks the
      !! conditional allocation directly).  Reported before `enter_data` and
      !! reconciled against the measured device mapping — a new array added
      !! without a matching `bytes()` term makes the measured map exceed this
      !! count and self-announces the drift (see rdb_mem_report).  Slots with
      !! no allocatables (eos, restart) carry no term.  The `diag%vars`
      !! registry USED to be excluded — it is now counted, per registered
      !! variable, by `diag_var_bytes` (it was the single largest uncounted
      !! device block: ~3.2 GB for the default catalog at 1000x800x50).
      !! `data_input` (PR-14) DOES have
      !! allocatables (`f0`/`f1` per registered field) and IS
      !! device-mapped, so it carries a real term below — zero registered
      !! fields on every shipped namelist keeps that term at 0.
      class(ocean_state_t), intent(in) :: this
      integer(int64) :: nbytes

      nbytes = this%barotropic%bytes() &
               + this%metrics%bytes() &
               + this%dyn%bytes() &
               + this%continuity%bytes() &
               + this%coriolis_adv%bytes() &
               + this%pressure_force%bytes() &
               + this%vmix%bytes() &
               + this%wavespeed%bytes() &
               + this%lateral_mix%bytes() &
               + this%hvisc%bytes() &
               + this%bdrag%bytes() &
               + this%tdrag%bytes() &
               + this%surface_stress%bytes() &
               + this%surface_flux%bytes() &
               + this%cavity_flux%bytes() &
               + this%vert_advect%bytes() &
               + this%hdiff_tracer%bytes() &
               + this%vdiff%bytes() &
               + this%vcoord%bytes() &
               + this%diag%bytes() &
               + this%epbl%bytes() &
               + this%kshear%bytes() &
               + this%vmix_tidal%bytes() &
               + this%slopes%bytes() &
               + this%mle%bytes() &
               + this%gm%bytes() &
               + this%redi%bytes() &
               + this%varmix%bytes() &
               + this%meke%bytes() &
               + this%ice%bytes() &
               + this%tides%bytes() &
               + this%p_surf%bytes() &
               + this%obc%bytes() &
               + this%bc%bytes() &
               + this%sponge%bytes()
#ifndef RDB_NO_NETCDF
      nbytes = nbytes + this%data_input%bytes()
#endif
      if (this%use_multilayer) nbytes = nbytes + this%multilayer%bytes()
   end function ocean_state_bytes

   subroutine ocean_state_enter_data(state)
      !! Attach the ocean god state's allocatables to the device.
      !! Map the parent struct first so the device knows the shape of
      !! `state`, then call each slot's bound `enter_data`.
      !! Only slots that own host allocations need to be visited; the
      !! Phase 0 shells (dyn, continuity, coriolis_adv, …) have no
      !! allocatables yet, so they're skipped until they ship their
      !! own bound methods.
      type(ocean_state_t), intent(inout) :: state

      !$acc enter data copyin(state)
      call profiler_start("ed_barotropic", nvtx_only=.true.)
      call state%barotropic%enter_data()
      call profiler_stop("ed_barotropic")
      call profiler_start("ed_metrics", nvtx_only=.true.)
      call state%metrics%enter_data()
      call profiler_stop("ed_metrics")
      if (state%use_multilayer) then
         call profiler_start("ed_multilayer", nvtx_only=.true.)
         call state%multilayer%enter_data()
         call profiler_stop("ed_multilayer")
      end if
      call profiler_start("ed_dyn", nvtx_only=.true.)
      call state%dyn%enter_data()
      call profiler_stop("ed_dyn")
      call profiler_start("ed_continuity", nvtx_only=.true.)
      call state%continuity%enter_data()
      call profiler_stop("ed_continuity")
      call profiler_start("ed_coriolis_adv", nvtx_only=.true.)
      call state%coriolis_adv%enter_data()
      call profiler_stop("ed_coriolis_adv")
      call profiler_start("ed_pressure_force", nvtx_only=.true.)
      call state%pressure_force%enter_data()
      call profiler_stop("ed_pressure_force")
      call profiler_start("ed_hvisc", nvtx_only=.true.)
      call state%hvisc%enter_data()
      call profiler_stop("ed_hvisc")
      call profiler_start("ed_bdrag", nvtx_only=.true.)
      call state%bdrag%enter_data()
      call state%tdrag%enter_data()
      call profiler_stop("ed_bdrag")
      call profiler_start("ed_surface_stress", nvtx_only=.true.)
      call state%surface_stress%enter_data()
      call profiler_stop("ed_surface_stress")
      call profiler_start("ed_surface_flux", nvtx_only=.true.)
      call state%surface_flux%enter_data()
      call state%cavity_flux%enter_data()
      call profiler_stop("ed_surface_flux")
      call profiler_start("ed_vert_advect", nvtx_only=.true.)
      call state%vert_advect%enter_data()
      call profiler_stop("ed_vert_advect")
      call profiler_start("ed_hdiff_tracer", nvtx_only=.true.)
      call state%hdiff_tracer%enter_data()
      call profiler_stop("ed_hdiff_tracer")
      call profiler_start("ed_vdiff", nvtx_only=.true.)
      call state%vdiff%enter_data()
      call profiler_stop("ed_vdiff")
      call profiler_start("ed_vmix", nvtx_only=.true.)
      call state%vmix%enter_data()
      call profiler_stop("ed_vmix")
      ! Gated closures are only device-mapped when enabled (their init is
      ! likewise gated, so the arrays are unallocated when off; unlike the
      ! kshear/tidal/epbl/meke slots these `enter_data` bodies carry no
      ! internal `allocated` guard, so the gate lives here at the call).
      call profiler_start("ed_epbl", nvtx_only=.true.)
      if (state%epbl%enable) call state%epbl%enter_data()
      call profiler_stop("ed_epbl")
      call profiler_start("ed_wavespeed", nvtx_only=.true.)
      call state%wavespeed%enter_data()
      call profiler_stop("ed_wavespeed")
      call profiler_start("ed_kshear", nvtx_only=.true.)
      if (state%kshear%enable) call state%kshear%enter_data()
      call profiler_stop("ed_kshear")
      call profiler_start("ed_vmix_tidal", nvtx_only=.true.)
      if (state%vmix_tidal%enable) call state%vmix_tidal%enter_data()
      call profiler_stop("ed_vmix_tidal")
      call profiler_start("ed_lateral_mix", nvtx_only=.true.)
      call state%lateral_mix%enter_data()
      call profiler_stop("ed_lateral_mix")
      call profiler_start("ed_slopes", nvtx_only=.true.)
      if (state%slopes%enable .or. state%gm%enable) call state%slopes%enter_data()
      call profiler_stop("ed_slopes")
      call profiler_start("ed_mle", nvtx_only=.true.)
      if (state%mle%enable) call state%mle%enter_data()
      call profiler_stop("ed_mle")
      call profiler_start("ed_gm", nvtx_only=.true.)
      if (state%gm%enable) call state%gm%enter_data()
      call profiler_stop("ed_gm")
      call profiler_start("ed_redi", nvtx_only=.true.)
      if (state%redi%enable) call state%redi%enter_data()
      call profiler_stop("ed_redi")
      call profiler_start("ed_varmix", nvtx_only=.true.)
      if (state%varmix%enable) call state%varmix%enter_data()
      call profiler_stop("ed_varmix")
      call profiler_start("ed_meke", nvtx_only=.true.)
      if (state%meke%enable) call state%meke%enter_data()
      call profiler_stop("ed_meke")
      call profiler_start("ed_vcoord", nvtx_only=.true.)
      call state%vcoord%enter_data()
      call profiler_stop("ed_vcoord")
      call profiler_start("ed_diag", nvtx_only=.true.)
      call state%diag%enter_data()
      call profiler_stop("ed_diag")
      call profiler_start("ed_tides", nvtx_only=.true.)
      call state%tides%enter_data()
      call profiler_stop("ed_tides")
      call profiler_start("ed_psurf", nvtx_only=.true.)
      call state%p_surf%enter_data()
      call profiler_stop("ed_psurf")
      call profiler_start("ed_ice", nvtx_only=.true.)
      if (state%ice%enable) call state%ice%enter_data()
      call profiler_stop("ed_ice")
      call profiler_start("ed_bc", nvtx_only=.true.)
      call ocean_bc_state_enter_data(state%bc)
      call profiler_stop("ed_bc")
      ! PR-23: gated (mirrors epbl/kshear/... — see the "Gated closures"
      ! comment above) — a disabled sponge (default) has unallocated arrays,
      ! so it must not be mapped.
      call profiler_start("ed_sponge", nvtx_only=.true.)
      if (state%sponge%enable) call state%sponge%enter_data()
      call profiler_stop("ed_sponge")
#ifndef RDB_NO_NETCDF
      call profiler_start("ed_data_input", nvtx_only=.true.)
      call state%data_input%enter_data()
      call profiler_stop("ed_data_input")
#endif
   end subroutine ocean_state_enter_data

   subroutine ocean_state_exit_data(state)
      !! Reverse of `ocean_state_enter_data`.  Walk components first,
      !! parent struct last.
      type(ocean_state_t), intent(inout) :: state

#ifndef RDB_NO_NETCDF
      call profiler_start("xd_data_input", nvtx_only=.true.)
      call state%data_input%exit_data()
      call profiler_stop("xd_data_input")
#endif
      call profiler_start("xd_sponge", nvtx_only=.true.)
      if (state%sponge%enable) call state%sponge%exit_data()
      call profiler_stop("xd_sponge")
      call profiler_start("xd_bc", nvtx_only=.true.)
      call ocean_bc_state_exit_data(state%bc)
      call profiler_stop("xd_bc")
      call profiler_start("xd_ice", nvtx_only=.true.)
      if (state%ice%enable) call state%ice%exit_data()
      call profiler_stop("xd_ice")
      call profiler_start("xd_psurf", nvtx_only=.true.)
      call state%p_surf%exit_data()
      call profiler_stop("xd_psurf")
      call profiler_start("xd_tides", nvtx_only=.true.)
      call state%tides%exit_data()
      call profiler_stop("xd_tides")
      call profiler_start("xd_diag", nvtx_only=.true.)
      call state%diag%exit_data()
      call profiler_stop("xd_diag")
      call profiler_start("xd_vcoord", nvtx_only=.true.)
      call state%vcoord%exit_data()
      call profiler_stop("xd_vcoord")
      call profiler_start("xd_meke", nvtx_only=.true.)
      if (state%meke%enable) call state%meke%exit_data()
      call profiler_stop("xd_meke")
      call profiler_start("xd_varmix", nvtx_only=.true.)
      if (state%varmix%enable) call state%varmix%exit_data()
      call profiler_stop("xd_varmix")
      call profiler_start("xd_redi", nvtx_only=.true.)
      if (state%redi%enable) call state%redi%exit_data()
      call profiler_stop("xd_redi")
      call profiler_start("xd_gm", nvtx_only=.true.)
      if (state%gm%enable) call state%gm%exit_data()
      call profiler_stop("xd_gm")
      call profiler_start("xd_mle", nvtx_only=.true.)
      if (state%mle%enable) call state%mle%exit_data()
      call profiler_stop("xd_mle")
      call profiler_start("xd_slopes", nvtx_only=.true.)
      if (state%slopes%enable .or. state%gm%enable) call state%slopes%exit_data()
      call profiler_stop("xd_slopes")
      call profiler_start("xd_lateral_mix", nvtx_only=.true.)
      call state%lateral_mix%exit_data()
      call profiler_stop("xd_lateral_mix")
      call profiler_start("xd_vmix_tidal", nvtx_only=.true.)
      if (state%vmix_tidal%enable) call state%vmix_tidal%exit_data()
      call profiler_stop("xd_vmix_tidal")
      call profiler_start("xd_kshear", nvtx_only=.true.)
      if (state%kshear%enable) call state%kshear%exit_data()
      call profiler_stop("xd_kshear")
      call profiler_start("xd_wavespeed", nvtx_only=.true.)
      call state%wavespeed%exit_data()
      call profiler_stop("xd_wavespeed")
      call profiler_start("xd_epbl", nvtx_only=.true.)
      if (state%epbl%enable) call state%epbl%exit_data()
      call profiler_stop("xd_epbl")
      call profiler_start("xd_vmix", nvtx_only=.true.)
      call state%vmix%exit_data()
      call profiler_stop("xd_vmix")
      call profiler_start("xd_vdiff", nvtx_only=.true.)
      call state%vdiff%exit_data()
      call profiler_stop("xd_vdiff")
      call profiler_start("xd_hdiff_tracer", nvtx_only=.true.)
      call state%hdiff_tracer%exit_data()
      call profiler_stop("xd_hdiff_tracer")
      call profiler_start("xd_vert_advect", nvtx_only=.true.)
      call state%vert_advect%exit_data()
      call profiler_stop("xd_vert_advect")
      call profiler_start("xd_surface_flux", nvtx_only=.true.)
      call state%cavity_flux%exit_data()
      call state%surface_flux%exit_data()
      call profiler_stop("xd_surface_flux")
      call profiler_start("xd_surface_stress", nvtx_only=.true.)
      call state%surface_stress%exit_data()
      call profiler_stop("xd_surface_stress")
      call profiler_start("xd_bdrag", nvtx_only=.true.)
      call state%tdrag%exit_data()
      call state%bdrag%exit_data()
      call profiler_stop("xd_bdrag")
      call profiler_start("xd_hvisc", nvtx_only=.true.)
      call state%hvisc%exit_data()
      call profiler_stop("xd_hvisc")
      call profiler_start("xd_pressure_force", nvtx_only=.true.)
      call state%pressure_force%exit_data()
      call profiler_stop("xd_pressure_force")
      call profiler_start("xd_coriolis_adv", nvtx_only=.true.)
      call state%coriolis_adv%exit_data()
      call profiler_stop("xd_coriolis_adv")
      call profiler_start("xd_continuity", nvtx_only=.true.)
      call state%continuity%exit_data()
      call profiler_stop("xd_continuity")
      call profiler_start("xd_dyn", nvtx_only=.true.)
      call state%dyn%exit_data()
      call profiler_stop("xd_dyn")
      if (state%use_multilayer) then
         call profiler_start("xd_multilayer", nvtx_only=.true.)
         call state%multilayer%exit_data()
         call profiler_stop("xd_multilayer")
      end if
      call profiler_start("xd_metrics", nvtx_only=.true.)
      call state%metrics%exit_data()
      call profiler_stop("xd_metrics")
      call profiler_start("xd_barotropic", nvtx_only=.true.)
      call state%barotropic%exit_data()
      call profiler_stop("xd_barotropic")
      !$acc exit data delete(state)
   end subroutine ocean_state_exit_data

   subroutine ocean_state_destroy(this)
      class(ocean_state_t), intent(inout) :: this
      this%is_init = .false.
      call this%restart%destroy()
#ifndef RDB_NO_NETCDF
      call this%data_input%destroy()
#endif
      call this%diag%destroy()
      call this%sponge%destroy()
      call ocean_bc_state_destroy(this%bc)
      call this%obc%destroy()
      call this%ice%destroy()
      call this%p_surf%destroy()
      call this%tides%destroy()
      call this%vcoord%destroy()
      call this%eos%destroy()
      call this%meke%destroy()
      call this%varmix%destroy()
      call this%redi%destroy()
      call this%gm%destroy()
      call this%mle%destroy()
      call this%slopes%destroy()
      call this%lateral_mix%destroy()
      call this%epbl%destroy()
      call this%wavespeed%destroy()
      call this%kshear%destroy()
      call this%vmix_tidal%destroy()
      call this%hvisc%destroy()
      call this%bdrag%destroy()
      call this%tdrag%destroy()
      call this%cavity_flux%destroy()
      call this%surface_flux%destroy()
      call this%surface_stress%destroy()
      call this%vert_advect%destroy()
      call this%hdiff_tracer%destroy()
      call this%vdiff%destroy()
      call this%vmix%destroy()
      call this%pressure_force%destroy()
      call this%coriolis_adv%destroy()
      call this%continuity%destroy()
      call this%dyn%destroy()
      call this%multilayer%destroy()
      call this%metrics%destroy()
      call this%barotropic%destroy()
   end subroutine ocean_state_destroy

   subroutine ocean_state_seed_from_cfg(state, grid, cfg, ierr, injected_b, injected_b_convention)
      !! Populate the ocean prognostic state with an analytical IC
      !! derived from cfg scalars.  Bathymetry is set per `cfg%ocean%topo%topo_config`
      !! (`"flat"` → uniform `ocean_max_depth`; `"spoon"` → MOM6 spoon
      !! shape) — UNLESS `injected_b` is present, in which case it
      !! overrides `topo_config` entirely (P2.5 pre-create geometry
      !! injection).  Water column thickness `h = b` so the free surface
      !! starts at SSH = 0.  Layers split the local depth evenly
      !! (`h_layer(i,j,k) = b(i,j) / nz_ml`).  Tracers carry per-layer
      !! `h * Tr` at the configured initial T/S.  Velocities zeroed.
      !!
      !! Caller must have run `state%init(grid)` first so every slot's
      !! allocations are in place; this routine only writes into them
      !! on the host.  Run before `ocean_state_enter_data` so the GPU
      !! mapping captures the seeded values.
      !!
      !! Every position-aware formula fill below reads its global offsets
      !! and global extents from `grid` (`i_offset_global` /
      !! `j_offset_global`, `nx_global` / `ny_global`), so each rank seeds
      !! its own window on ONE global analytical IC.  There is deliberately
      !! no `decomp` argument and no optional offset quartet: an optional
      !! that defaults to the LOCAL extents is silently wrong under MPI, and
      !! forgetting to pass it shipped three separate bugs.
      type(ocean_state_t), intent(inout) :: state
      type(hgrid_t), intent(in) :: grid
      type(config_t), intent(in) :: cfg
      integer, intent(out), optional :: ierr
         !! Non-zero on an initial-condition configuration/data failure
         !! when present; absent behaves as today (`error stop`).
      real(wp), intent(in), optional :: injected_b(:, :)
         !! P2.5 pre-create geometry injection: interior-sized
         !! `(grid%nx_phys, grid%ny_phys)` bathymetry, RAW in the caller's
         !! own sign convention (see `injected_b_convention`) — REQUIRES
         !! `injected_b_convention`.  When present, this OVERRIDES
         !! `cfg%ocean%topo%topo_config` entirely: sign-normalised to
         !! positive-down depth + wet-fraction validated
         !! (`bathymetry_normalise_sign`), written into the interior of
         !! `state%barotropic%b`, and ghost-filled by constant
         !! extrapolation (`bathymetry_fill_ghosts_array`) — mirroring
         !! `topo_config='file'` rather than trusting the caller
         !! (CLAUDE.md formula-bathymetry ghost-fill gotcha).  Everything
         !! downstream in this routine (h=b, layer split, wet mask, T/S
         !! seed, IC overlays, the ALE z_ref table) then runs unchanged
         !! against the injected bathymetry.
      integer, intent(in), optional :: injected_b_convention
         !! REQUIRED alongside `injected_b` — no default (D6.2: sign is
         !! the single most dangerous argument in the geometry API). One
         !! of `BATHY_CONVENTION_DEPTH_POSITIVE_DOWN` /
         !! `_HEIGHT_POSITIVE_UP` (`rdb_ocean_bathymetry_inject`).

      integer :: nz_ml, idx_S, idx_T, i, j, k, nx, ny, local_ierr
      real(wp), allocatable :: water(:, :)
         !! Reference water-column thickness the IC seeds work on:
         !! `b − z_draft` under an ice shelf, a byte copy of `b`
         !! otherwise.  Host-only setup scratch, released on return.

      nz_ml = state%multilayer%nz_ml
      idx_S = state%multilayer%idx_salinity
      idx_T = state%multilayer%idx_temperature
      nx = size(state%barotropic%b, 1)
      ny = size(state%barotropic%b, 2)

      ! Bathymetry.  `slope_scale` (the spoon/seamount length scale) is a
      ! metres knob; `set_bathymetry_*` works in GRID coordinate units,
      ! which are metres on Cartesian grids but DEGREES on spherical/
      ! curvilinear ones.  Convert here so the length scale is consistent
      ! with the grid coordinate (without it, a metres scale against a
      ! degrees position collapses the seamount/spoon to a flat basin).
      if (present(injected_b)) then
         ! P2.5 pre-create geometry injection — bypasses topo_config
         ! entirely.  `injected_b_convention` has no default (D6.2).
         if (.not. present(injected_b_convention)) then
            call fail("ocean_state_seed_from_cfg: injected_b requires "// &
                      "injected_b_convention (no default — sign is the single most "// &
                      "dangerous argument; see D6.2 in the Python API design notes)", &
                      ierr, OCEAN_STATUS_ERR_BAD_SHAPE)
            return
         end if
         if (size(injected_b, 1) /= grid%nx_phys .or. size(injected_b, 2) /= grid%ny_phys) then
            call fail("ocean_state_seed_from_cfg: injected_b shape mismatch — got ("// &
                      to_string(size(injected_b, 1))//","//to_string(size(injected_b, 2))// &
                      "), expected the physical interior ("//to_string(grid%nx_phys)//","// &
                      to_string(grid%ny_phys)//")", ierr, OCEAN_STATUS_ERR_BAD_SHAPE)
            return
         end if
         block
            real(wp) :: b_local(grid%nx_phys, grid%ny_phys)
            integer :: ng
            b_local = injected_b
            call bathymetry_normalise_sign(b_local, injected_b_convention, ierr=local_ierr)
            if (local_ierr /= OCEAN_STATUS_OK) then
               if (present(ierr)) then
                  ierr = local_ierr
                  return
               end if
               error stop "ocean_state_seed_from_cfg: injected_b sign normalisation failed"
            end if
            ng = grid%nghost
            state%barotropic%b(ng + 1:ng + grid%nx_phys, ng + 1:ng + grid%ny_phys) = b_local
            call bathymetry_fill_ghosts_array(state%barotropic%b, grid)
         end block
      else
         select case (trim(cfg%ocean%topo%topo_config))
         case ("spoon")
            call set_bathymetry_spoon(state%barotropic%b, grid, &
                                      cfg%ocean%topo%max_depth, cfg%ocean%topo%edge_depth, &
                                      topo_length_to_grid_units(cfg%ocean%topo%slope_scale, &
                                                                cfg%ocean%grid%grid_config, &
                                                                cfg%ocean%grid%rad_earth))
         case ("seamount")
            call set_bathymetry_seamount(state%barotropic%b, grid, &
                                         cfg%ocean%topo%max_depth, cfg%ocean%topo%edge_depth, &
                                         topo_length_to_grid_units(cfg%ocean%topo%slope_scale, &
                                                                   cfg%ocean%grid%grid_config, &
                                                                   cfg%ocean%grid%rad_earth))
         case ("neverworld2")
            call set_bathymetry_neverworld2(state%barotropic%b, grid, &
                                            cfg%ocean%topo%max_depth, &
                                            cfg%ocean%topo%nl_continent_amp, &
                                            cfg%ocean%topo%nl_roughness_amp, &
                                            cfg%ocean%topo%nl_min_depth)
         case ("island")
            ! Flat basin at `max_depth` with a central LAND square (depth
            ! forced below LAND_DEPTH_THRESHOLD).  Produces the interior
            ! land block the static land-mask validation needs.  The square
            ! spans the central `island_frac` fraction of the physical
            ! domain in each direction (default via edge_depth/slope_scale
            ! re-use; see set_bathymetry_island).
            call set_bathymetry_island(state%barotropic%b, grid, &
                                       cfg%ocean%topo%max_depth, &
                                       cfg%ocean%topo%slope_scale)
         case ("double_drake")
            ! Ferreira et al. (2010) two-wall supercontinent with a
            ! reentrant southern channel — the land-mask showcase.
            call set_bathymetry_double_drake(state%barotropic%b, grid, &
                                             cfg%ocean%topo%max_depth, &
                                             cfg%ocean%topo%slope_scale)
         case ("file")
            ! Real bathymetry from NetCDF.  File must be pre-projected onto
            ! the model's Cartesian grid (matching nx_phys × ny_phys); the
            ! loader fills ghost rows by constant extrapolation.  Sign
            ! convention: `b` is bottom depth positive-down (consistent with
            ! the spoon + flat branches).  GEBCO + ETOPO ship elevation
            ! positive-up — the preprocessing script should flip the sign
            ! before writing the NetCDF.
            if (len_trim(cfg%bathymetry_file) == 0) then
               call fail("ocean_state_seed_from_cfg: topo_config='file' requires "// &
                         "bathymetry_file", ierr, OCEAN_STATUS_ERR_IC_SEED)
               return
            end if
#ifndef RDB_NO_NETCDF
            if (present(ierr)) then
               call load_bathymetry_into_array(trim(cfg%bathymetry_file), &
                                               state%barotropic%b, grid, ierr=local_ierr)
               if (local_ierr /= OCEAN_STATUS_OK) then
                  ierr = local_ierr
                  return
               end if
            else
               call load_bathymetry_into_array(trim(cfg%bathymetry_file), &
                                               state%barotropic%b, grid)
            end if
#else
            ! File-loaded bathymetry needs the NetCDF reader in `rdb_bathymetry`,
            ! which isn't compiled in when RDB_ENABLE_NETCDF=OFF.  Fail fast
            ! with a descriptive error so the user knows to rebuild with NetCDF
            ! enabled rather than hitting a less-obvious runtime issue.
            call fail("ocean_state_seed_from_cfg: topo_config='file' requires "// &
                      "RDB_ENABLE_NETCDF=ON at build time (the NetCDF-backed "// &
                      "rdb_bathymetry reader is needed to load the file).", ierr, OCEAN_STATUS_ERR_IO)
            return
#endif
         case default  ! "flat"
            state%barotropic%b = cfg%ocean%topo%max_depth
         end select
      end if

      ! ---- Static ice-shelf cavity geometry (&ocean_cavity_dyn_nml, P5.1) ----
      ! ORDERING IS LOAD-BEARING and this is the only place it can go: the
      ! draft must exist before the water column `b − z_draft` seeds the
      ! layer split and the wet mask (grounding), both of which happen a
      ! few lines below, and long before any `configure_ocean_*` pass.
      ! The formula setters fill the FULL array including ghosts; the
      ! periodic/fold re-wrap + halo exchange that `barotropic%b` gets are
      ! applied to `z_draft` alongside it in `rdb_ocean_engine`.
      !
      ! `water` is the reference water-column thickness every seed below
      ! works on.  With the knob off it is a byte copy of `b`, so there is
      ! ONE code path and the default run is bit-identical.
      allocate (water(nx, ny))
      if (state%metrics%use_cavity) then
         call seed_cavity_draft(state, grid, cfg, ierr=local_ierr)
         if (local_ierr /= OCEAN_STATUS_OK) then
            if (present(ierr)) then
               ierr = local_ierr
               return
            end if
            error stop "ocean_state_seed_from_cfg: ice-shelf cavity geometry failed"
         end if
         call cavity_water_column_impl(water, state%barotropic%b, &
                                       state%metrics%z_draft, nx, ny)
      else
         water = state%barotropic%b
      end if

      ! Water column thickness h = b − z_draft → the free-surface anomaly
      ! `bt_eta = Σ h_layer − bt_H_ref` starts at zero (SSH = 0 without a
      ! cavity; the loaded equilibrium under one).
      state%barotropic%h = water
      state%barotropic%u_face_x = 0.0_wp
      state%barotropic%v_face_y = 0.0_wp
      state%barotropic%hu_face_x = 0.0_wp
      state%barotropic%hv_face_y = 0.0_wp

      ! Per-column layer thickness = depth / nz_ml.  Pass flat
      ! allocatables into a `_impl` helper so the `do concurrent` body
      ! reads/writes plain arrays — registry / array-of-DT indirection
      ! crashes NVHPC's device codegen with CUDA_ERROR_ILLEGAL_ADDRESS.
      ! With wet/dry ON, floor the seed to 2·H_VANISHED so an emerged
      ! intertidal rest-bed (b < 0 ⇒ now a LIVE wet_mask=1 column) never
      ! seeds a negative layer / negative barotropic depth D = Σ h_layer,
      ! and each seeded layer clears the strict `h_old > H_VANISHED`
      ! remap-drain gate (so the seeded S/T survives the first regrid).
      ! Tracers are seeded from THIS floored h_layer just below, so hTr
      ! stays consistent (const·2·H_VANISHED, finite T/S, no negative mass).
      ! `thickness_config = "uniform_z"` swaps the local-depth even split for
      ! MOM6's uniform-z-interface seed (flat resting isopycnals under a
      ! horizontally-uniform density stack).  Validated as mutually exclusive
      ! with wet/dry, so the two branches never both need the emerged-column
      ! floor.  Default "sigma" ⇒ byte-identical to the pre-knob path.
      if (trim(cfg%thickness_config) == "uniform_z") then
         call seed_h_layer_uniform_z_impl(state%multilayer%h_layer, &
                                          water, nz_ml, &
                                          cfg%ocean%topo%max_depth, &
                                          cfg%ocean%isopycnal%angstrom_h)
      else
         call seed_h_layer_uniform_impl(state%multilayer%h_layer, &
                                        water, nz_ml, &
                                        apply_wetdry_floor=cfg%ocean%wetdry%enable)
      end if
      ! Keep the barotropic water-column prognostic non-negative on the
      ! same emerged band (it is inert in the multilayer dyn-core — SSH is
      ! Σ h_layer − b — but is restart-registered as `bt_h`; floor it so a
      ! checkpoint never carries a negative rest depth).  Match the layer
      ! floor exactly: on an emerged column Σ h_layer = nz·2·H_VANISHED, so
      ! floor bt_h to the SAME value (not 0) — otherwise a restart that
      ! re-derives D from bt_h would disagree with Σ h_layer by nz·2·H_VANISHED
      ! on every emerged column.  Knob-off ⇒ h = b (byte-identical).
      if (cfg%ocean%wetdry%enable) then
         state%barotropic%h = max(water, &
                                  real(nz_ml, wp)*2.0_wp*H_VANISHED)
      end if
      state%multilayer%u_face_x_layer = 0.0_wp
      state%multilayer%v_face_y_layer = 0.0_wp
      state%multilayer%hu_face_x_layer = 0.0_wp
      state%multilayer%hv_face_y_layer = 0.0_wp

      ! Wet/dry mask: 1.0 where `b >= cutoff`, 0.0 elsewhere (land).
      ! Default cutoff = LAND_DEPTH_THRESHOLD (byte-identical for all
      ! non-wetdry configs).  When wetdry is enabled, pass `land_cutoff =
      ! -land_margin` so intertidal columns (bed above rest MSL but below the
      ! flood headroom) stay wet_mask=1 and the dynamic wd_wet_dyn gate
      ! handles their wetting/drying instead of the static land mask.
      if (cfg%ocean%wetdry%enable) then
         call seed_wet_mask_impl(state%multilayer%wet_mask, water, &
                                 land_cutoff=-cfg%ocean%wetdry%land_margin)
      else if (state%metrics%use_cavity) then
         ! GROUNDING.  A column with less than `h_min_cavity` of water
         ! under the ice is LAND — routed through the SAME wet-mask seed
         ! the bathymetry uses, so the static metric-zeroing land mask
         ! (`configure_ocean_land_mask`) and the finite land-state hold
         ! (`ocean_state_seed_land_cells`) follow for free.  Never a thin
         ! film of water under grounded ice.  Note the cutoff is applied
         ! to `b − z_draft`, so an ordinary land column (`b` below
         ! LAND_DEPTH_THRESHOLD, draft already zeroed there) is land for
         ! the same reason it always was.
         call seed_wet_mask_impl(state%multilayer%wet_mask, water, &
                                 land_cutoff=cfg%ocean%cavity_dyn%h_min_cavity)
      else
         call seed_wet_mask_impl(state%multilayer%wet_mask, water)
      end if

      ! Tracers carry per-layer h*Tr.  Multiply the (now spatially-
      ! varying) h_layer by the configured uniform scalar.
      if (idx_S > 0) then
         block
            real(wp) :: s_layer(nz_ml)
            real(wp) :: dS_dlayer
            logical :: stratify_s
            ! Linear S(z) when both surface + bottom are specified —
            ! EXACT mirror of the temperature branch below, including
            ! the gate (`both /= 0`), the vertical convention (k=1 is
            ! the bed = S_init_bottom, k=nz_ml the surface =
            ! S_init_surface), the single-layer fall-through, and the
            ! seed helper (so ghosts and land columns are filled the
            ! same way: hTr = S(k)·h_layer everywhere, land included,
            ! since h_layer is already 0/floored there).  The profile is
            ! linear in LAYER INDEX, which under the sigma-style
            ! `h_layer = b/nz_ml` seed is linear in layer-centre depth
            ! on every column — so a sloping bed gets the same endpoint
            ! values with a depth-proportional gradient.
            ! Note the stable polarity is the INVERSE of temperature:
            ! dense/salty water belongs at the bed, so a stable haline
            ! column has `S_init_bottom > S_init_surface`.
            ! Both-zero (the default) ⇒ uniform `initial_salinity`,
            ! byte-identical to the pre-knob path.
            stratify_s = (cfg%S_init_surface /= 0.0_wp) .and. &
                         (cfg%S_init_bottom /= 0.0_wp) .and. &
                         (nz_ml > 1)
            if (stratify_s) then
               dS_dlayer = (cfg%S_init_surface - cfg%S_init_bottom)/ &
                           real(nz_ml - 1, wp)
               do k = 1, nz_ml
                  s_layer(k) = cfg%S_init_bottom + dS_dlayer*real(k - 1, wp)
               end do
               call seed_tracer_stratified_impl( &
                  state%multilayer%tracers(idx_S)%hTr, &
                  state%multilayer%h_layer, s_layer, nz_ml)
            else
               call seed_tracer_uniform_impl( &
                  state%multilayer%tracers(idx_S)%hTr, &
                  state%multilayer%h_layer, &
                  cfg%initial_salinity, nz_ml)
            end if
         end block
      end if
      if (idx_T > 0) then
         block
            real(wp) :: t_layer(nz_ml)
            real(wp) :: dT_dlayer
            logical :: stratify
            ! Linear T(z) when both surface + bottom are specified.
            ! Convention: k=1 is the bed (T = T_init_bottom), k=nz_ml is
            ! the surface (T = T_init_surface).  Single-layer case
            ! falls through to T_init_surface (no profile to build).
            stratify = (cfg%T_init_surface /= 0.0_wp) .and. &
                       (cfg%T_init_bottom /= 0.0_wp) .and. &
                       (nz_ml > 1)
            if (stratify) then
               dT_dlayer = (cfg%T_init_surface - cfg%T_init_bottom)/ &
                           real(nz_ml - 1, wp)
               do k = 1, nz_ml
                  t_layer(k) = cfg%T_init_bottom + dT_dlayer*real(k - 1, wp)
               end do
            else
               t_layer = cfg%initial_temperature
            end if
            call seed_tracer_stratified_impl( &
               state%multilayer%tracers(idx_T)%hTr, &
               state%multilayer%h_layer, t_layer, nz_ml)
         end block
      end if

      ! Optional IC overlay applied after bathymetry + uniform h_layer
      ! + analytical T/S are in place.
      ! ierr threaded down ONLY when THIS routine's own ierr is present:
      ! otherwise each seed_*_ic helper must keep reaching its own
      ! `error stop` (specific text) rather than the generic wrapper
      ! message below (P0.1 review F2).
      if (present(ierr)) then
         select case (trim(cfg%ocean%ic%ic_config))
         case ("eady")
            call seed_eady_ic(state, grid, cfg, ierr=local_ierr)
         case ("geostrophic_adjustment")
            call seed_geostrophic_adjustment_ic(state, grid, cfg, ierr=local_ierr)
         case ("baroclinic_jet")
            call seed_baroclinic_jet_ic(state, grid, cfg, ierr=local_ierr)
         case default
            ! "" — keep the default analytical IC unchanged.
            local_ierr = 0
         end select
         if (local_ierr /= 0) then
            ierr = local_ierr
            return
         end if
      else
         select case (trim(cfg%ocean%ic%ic_config))
         case ("eady")
            call seed_eady_ic(state, grid, cfg)
         case ("geostrophic_adjustment")
            call seed_geostrophic_adjustment_ic(state, grid, cfg)
         case ("baroclinic_jet")
            call seed_baroclinic_jet_ic(state, grid, cfg)
         case default
            ! "" — keep the default analytical IC unchanged.
         end select
      end if

      ! Z-level T/S IC overlay (capability A2).  Runs after bathymetry +
      ! uniform h_layer + wet_mask + analytical T/S are in place (it reads
      ! the seeded column depths + wet_mask) and intentionally OVERWRITES
      ! any analytical T/S.  Default off (`enable = .false.`) preserves
      ! bit-identity.  NetCDF-only: the reader lives in rdb_ocean_z_init,
      ! which only compiles with RDB_ENABLE_NETCDF=ON.
      !
      ! Under a cavity the overlay is handed `metrics%z_draft` so every
      ! layer centre's depth is measured from `z = 0` rather than from the
      ! ice base; see `seed_zinit_overlay`.
      if (cfg%ocean%zinit%enable) then
#ifndef RDB_NO_NETCDF
         if (present(ierr)) then
            call seed_zinit_overlay(state, grid, cfg, ierr=local_ierr)
            if (local_ierr /= 0) then
               ierr = local_ierr
               return
            end if
         else
            call seed_zinit_overlay(state, grid, cfg)
         end if
#else
         call fail("ocean_state_seed_from_cfg: ocean_zinit requires "// &
                   "RDB_ENABLE_NETCDF=ON at build time (the NetCDF-backed "// &
                   "rdb_ocean_z_init reader is needed to load the z-level T/S file).", ierr, OCEAN_STATUS_ERR_IO)
         return
#endif
      end if

      ! Pseudo-salt seed — MUST run after every write to salinity's initial
      ! condition above (analytical IC, the eady/geostrophic/baroclinic-jet
      ! overlays, and critically the z-file overlay, which OVERWRITES any
      ! analytical S) and before ocean_state_enter_data (driver-ordered).
      ! Seeding before the z-file overlay would leave D(t=0) /= 0 in exactly
      ! the configuration (z-file IC) where the diagnostic matters most.
      ! Self-gates on idx_pseudo_salt <= 0 (knob off).
      call ocean_pseudo_salt_seed(state%multilayer)

      ! Direct per-layer density init — MOM6 `COORD_CONFIG="gprime"` analogue.
      ! When `ocean_layer_rho_init` is set (any non-sentinel value), write
      ! `ms%rho_layer(:,:,k) = ocean_layer_rho_init(k)` directly, bypassing
      ! the EOS path entirely.  Intended use case: `ocean_enable_thermodynamics
      ! = .false.` + this knob ⇒ static reduced-gravity stratification, mirroring
      ! MOM6's adiabatic gprime IC.  Without this, `rho_layer` stays at its
      ! alloc-time zero when thermo is off (the dyn step's EOS call is gated
      ! on `therm_active`), and FV-LITE's pressure-stack collapses to zero
      ! ⇒ no PGF, no Sverdrup balance, no WBC dynamics — only the wind +
      ! Coriolis side of the momentum equation produces anything.  Setting
      ! the per-layer ρ here is what turns the PGF back on under the
      ! adiabatic-stratified semantics MOM6 uses for `double_gyre`.
      !
      ! Convention: `k=1` bed → `k=nz_ml` surface (Roundabout bottom-up).  Pass
      ! the heavier value first.  Count of non-sentinel entries must match
      ! `nz_ml`; otherwise we abort with a descriptive error.
      !
      ! `rho_lightest >= 0` selects the linear density-range generator
      ! (MOM6 `COORD_CONFIG="linear"`): build `nz_ml` linearly-spaced
      ! densities and reuse the same write path, so the column scales by
      ! just bumping `nz_layers`.  Mutually exclusive with the explicit
      ! `layer_rho_init` list — set only one.
      if (cfg%ocean%ic%rho_lightest >= 0.0_wp) then
         if (count(cfg%ocean%ic%layer_rho_init >= 0.0_wp) > 0) then
            call fail("ocean_ic: rho_lightest (linear density-range) and "// &
                      "layer_rho_init (explicit list) are mutually exclusive — "// &
                      "set only one.", ierr, OCEAN_STATUS_ERR_IC_SEED)
            return
         end if
         if (present(ierr)) then
            call apply_layer_rho_init(state%multilayer, &
                                      ocean_linear_layer_density(cfg%ocean%ic%rho_lightest, &
                                                                 cfg%ocean%ic%rho_range, nz_ml), &
                                      nz_ml, ierr=local_ierr)
            if (local_ierr /= 0) then
               ierr = local_ierr
               return
            end if
         else
            call apply_layer_rho_init(state%multilayer, &
                                      ocean_linear_layer_density(cfg%ocean%ic%rho_lightest, &
                                                                 cfg%ocean%ic%rho_range, nz_ml), &
                                      nz_ml)
         end if
      else
         if (present(ierr)) then
            call apply_layer_rho_init(state%multilayer, cfg%ocean%ic%layer_rho_init, nz_ml, &
                                      ierr=local_ierr)
            if (local_ierr /= 0) then
               ierr = local_ierr
               return
            end if
         else
            call apply_layer_rho_init(state%multilayer, cfg%ocean%ic%layer_rho_init, nz_ml)
         end if
      end if

      ! Populate the per-column z_ref table for VCOORD_ZSTAR_FULL.  Without
      ! this, `compute_target_h(ZSTAR_FULL)` walks an all-zero table and
      ! the first ALE remap collapses every layer to `zstar_h_min` — bug
      ! caught by the double_gyre adiabatic NK=2 setup (2026-05-26).
      ! Other vcoord types ignore z_ref so the call is safe regardless.
      call state%vcoord%build_zref_full(state%barotropic%b)
      if (present(ierr)) ierr = OCEAN_STATUS_OK
   end subroutine ocean_state_seed_from_cfg

#ifndef RDB_NO_NETCDF
   subroutine seed_zinit_overlay(state, grid, cfg, ierr)
      !! Dispatch the `&ocean_zinit_nml` T/S overlay across its two axes:
      !! the profile SOURCE (`"file"` — the pre-regridded NetCDF reader;
      !! `"linear"` — the analytic affine `lin_*` profile) and whether a
      !! cavity draft is present.
      !!
      !! The draft is the whole reason this is a separate routine.  Both
      !! seeders measure each layer centre's GEOPOTENTIAL depth from
      !! `z = 0`; under an ice shelf the column top is `z_draft` metres
      !! down, so `metrics%z_draft` has to reach them or a `T(z)` profile
      !! lands systematically too shallow (and, under a SLOPING lid,
      !! tilts the isopycnals with the ice base — not a state of rest).
      !! With the cavity off `metrics%z_draft` is the `(1, 1)`
      !! placeholder, so the argument is simply not passed and the
      !! arithmetic is bit-identical to the pre-cavity path.
      !!
      !! NOT `pure`: the seeders it dispatches to read files and log.
      type(ocean_state_t), intent(inout) :: state
      type(hgrid_t), intent(in) :: grid
      type(config_t), intent(in) :: cfg
      integer, intent(out), optional :: ierr
         !! Threaded straight through to the seeder; an ABSENT `ierr`
         !! stays absent there, so each keeps its own `error stop` text.

      logical :: cav

      cav = state%metrics%use_cavity
      if (trim(adjustl(cfg%ocean%zinit%source)) == "linear") then
         if (cav) then
            call seed_ts_linear_z(state%multilayer, cfg%ocean%zinit, ierr, &
                                  state%metrics%z_draft)
         else
            call seed_ts_linear_z(state%multilayer, cfg%ocean%zinit, ierr)
         end if
      else
         if (cav) then
            call seed_ts_from_zfile(state%multilayer, grid, cfg%ocean%zinit, ierr, &
                                    state%metrics%z_draft)
         else
            call seed_ts_from_zfile(state%multilayer, grid, cfg%ocean%zinit, ierr)
         end if
      end if
   end subroutine seed_zinit_overlay
#endif

   subroutine ocean_state_build_restart_registry(state, grid, reg)
      !! Walk the ocean god state and register every field that must
      !! checkpoint for a bit-exact step-(N+1) resume (ROADMAP A1).
      !!
      !! Registered (always):
      !!   barotropic prognostics h, u_face_x, v_face_y;
      !!   multilayer prognostics h_layer, u/v_face_x/y_layer, every
      !!   tracers(:)%hTr (registry-of-registries — iterate the tracer
      !!   registry so passive tracers join automatically).
      !! Registered (conditionally — only when the slot's persistent
      !!   array is allocated):
      !!   vmix%bl_depth — KPP's lagged BL-depth seed; the V_t² wstar3
      !!     term reads the PREVIOUS step's value, so a restart between
      !!     steps reads it stale unless checkpointed.  Allocated
      !!     unconditionally (alloc tracks slot existence, not whether the
      !!     closure is active), so it always registers.
      !!   epbl%mld + epbl%kd_int, kshear%kd_int + kshear%tke_int.
      !!     These refresh only at thermo cadence but are *merged every
      !!     stage*, so a restart between refreshes reads them stale.
      !!     Also allocated unconditionally — register-always; the
      !!     enable-state itself is fixed at config and validated by the
      !!     schema metadata, so a flipped flag is a fresh run, not a
      !!     resume (flags cannot flip mid-run).
      !!   ml_rho_layer — diagnosed each stage from S/T, BUT when
      !!     `dt_therm_ratio > 1` it is NOT recomputed on the slow steps
      !!     that skip the EOS, so the value carried into the next step
      !!     depends on the last thermo step.  Cheapest correct fix:
      !!     checkpoint it (review #2).
      !!   bc%eta_old_chapman_{w,e,s,n} — Chapman radiation persistent
      !!     corner state (host scalars).  bc%tres_* tracer reservoirs —
      !!     registered when allocated (device-mapped, so they join the
      !!     update-self walk).
      !!
      !! Deliberately EXCLUDED:
      !!   sf%Q_heat / sf%Q_salt — DERIVED VIEWS (PR-12), not raw state.
      !!     With the surface-flux component set off (default), they are
      !!     re-seeded from `&ocean_thermo_nml q_heat / q_salt` (or the
      !!     Area-A3/A4 fill path) on every resume, same as before PR-12.
      !!     With the component set on, `ocean_surface_flux_assemble`
      !!     REBUILDS them every thermo step from `Q_heat_const` /
      !!     `Q_salt_const` plus the component fields — they fall under
      !!     the SAME derived-field exclusion rule as `w_interface` /
      !!     `mass_flux_*` below ("diagnosed each stage from the
      !!     prognostics before first use"), just diagnosed from the
      !!     component fields instead.  Registering a derived field
      !!     creates the classic file-vs-namelist ambiguity on resume
      !!     (does a checkpointed `Q_heat` win over an edited namelist
      !!     scalar?) — the correct fix is NOT to register `Q_heat`/
      !!     `Q_salt` but to require that a time-varying COMPONENT
      !!     register ITSELF: any filler (a future reader, the ice
      !!     coupler's `salt_flux`/`heat_added`, ...) that makes a
      !!     component time-varying MUST register it via
      !!     `registry_register_2d(tag, arr, ng, nx, ny, optional=.true.,
      !!     device_mapped=.true.)` from this subroutine.  In v1 no
      !!     component is time-varying on its own (the ice coupler's
      !!     restart-exact refill source is `ice_salt_flux_diag` /
      !!     `ice_heat_flux_diag`, registered below), so nothing
      !!     additional is registered here yet.
      !!   RK saves (h_layer0, hTr0, *_layer0) + bt accumulators
      !!     (bt_eta/ubt/...) — step-internal scratch, re-zeroed at the
      !!     top of every outer step (restart is step-aligned).
      !!   bt_H_ref — recomputed from `b` by configure_ocean_bt_split;
      !!     `b` is bathymetry (static), reconstructed at config time.
      !!   metrics%z_draft / cover_frac / p_ice_ref — the ice-shelf cavity
      !!     statics (`&ocean_cavity_dyn_nml`) fall under the SAME rule,
      !!     and for the same reason: the draft is a PRESCRIBED, static
      !!     geometry, rebuilt at configure by `seed_cavity_draft` from
      !!     the namelist before the restart read runs, exactly as `b` and
      !!     `bt_H_ref` are.  Checkpointing it would create the
      !!     file-vs-namelist ambiguity the derived-field exclusion exists
      !!     to avoid (does a saved draft beat an edited namelist?), and a
      !!     resume whose draft disagreed with its datum would be
      !!     silently wrong — `bt_H_ref = b − z_draft` ties the two
      !!     together, so they must be rebuilt together or not at all.
      !!     A TIME-VARYING draft (a coupled ice sheet) is a different
      !!     field with a different owner and would register itself, the
      !!     same way a time-varying surface-flux component must.
      !!   vcoord target_h / z_ref — recomputed per step (z_ref rebuilt
      !!     from `b` in the seed); never prognostic.
      !!   w_interface, mass_flux_* — diagnosed each stage from the
      !!     prognostics before first use.
      !!   diag accumulators — MEAN/MAX windows RESET on restart by
      !!     design (documented); the bit-exact gate covers the
      !!     prognostic trajectory, not mid-window diag aggregates.
      !!
      !! The dead `ocean_obc_t` scaffold (`state%obc`) is NOT registered:
      !! it is never enabled, never device-mapped, and has no enter_data —
      !! registering it (and issuing `update self` on its host-only
      !! buffers) was a latent GPU crash (review #3/#8).  When OBC lands,
      !! it registers its own live persistent state here.
      type(ocean_state_t), intent(in), target :: state
      type(hgrid_t), intent(in) :: grid
      type(restart_registry_t), intent(inout) :: reg

      integer :: it
      character(len=64) :: tag

      call reg%clear()

      ! Ghost-cell policy (v1): the structured prognostic + closure
      ! arrays are checkpointed as FULL local arrays (interior + ghosts),
      ! registered with ng=0 over the total extent.  Rationale — three
      ! kinds of ghost cell, distinguished for the MPI-readiness contract:
      !   * physical domain-edge (wall) ghosts are GENUINE owned boundary
      !     state — the slow-tendency operators (EOS / face-thickness /
      !     hvisc) read `h_layer` etc. at wall-adjacent ghost faces and
      !     there is NO stage-entry wall-mirror op to re-establish them,
      !     so owned-only loses information and breaks the bit-exact gate.
      !   * periodic / tripolar-fold seam ghosts are redundant — the
      !     stage-entry wrap (`ocean_periodic_wrap_state` / fold) rebuilds
      !     them from the interior every step; saving them is harmless.
      !   * inter-rank halo ghosts (E1, multi-rank) will be refilled by
      !     the halo exchange on resume; under the v1 SAME-decomp rule
      !     those columns round-trip identically anyway.
      ! So a full-array write is bit-exact today and forward-compatible
      ! with E1 (the halo columns simply get overwritten by exchange).
      ! The owned-slice machinery in restart_entry_t is retained for the
      ! ghostless edge buffers (OBC rings) and for a future owned-only
      ! mode once a wall-ghost-rebuild op exists.
      ! --- Barotropic prognostics ---
      call reg%register_2d("bt_h", state%barotropic%h, 0, &
                           size(state%barotropic%h, 1), size(state%barotropic%h, 2))
      call reg%register_2d("bt_u_face_x", state%barotropic%u_face_x, 0, &
                           size(state%barotropic%u_face_x, 1), size(state%barotropic%u_face_x, 2))
      call reg%register_2d("bt_v_face_y", state%barotropic%v_face_y, 0, &
                           size(state%barotropic%v_face_y, 1), size(state%barotropic%v_face_y, 2))

      if (state%use_multilayer) then
         ! --- Multilayer prognostics ---
         call register_full_3d(reg, "ml_h_layer", state%multilayer%h_layer)
         call register_full_3d(reg, "ml_u_face_x_layer", state%multilayer%u_face_x_layer)
         call register_full_3d(reg, "ml_v_face_y_layer", state%multilayer%v_face_y_layer)

         ! --- Tracers (registry-of-registries) ---
         if (allocated(state%multilayer%tracers)) then
            do it = 1, size(state%multilayer%tracers)
               write (tag, "(A,I0)") "tracer_hTr_", it
               call register_full_3d(reg, trim(tag), state%multilayer%tracers(it)%hTr)
            end do
         end if
         ! --- pred_corr step time-means (SPEC S1).  Under
         !     `&ocean_bt_nml split_scheme="pred_corr"` these are CARRIED
         !     STATE, not scratch: the predictor's Coriolis-advection and
         !     horizontal viscosity read u_av/v_av/h_av, and the dyn loop
         !     seeds them from the initial state only at
         !     `outer_step_count == 0`.  A resume restores outer_step_count
         !     from the file, so that seed does NOT re-fire -- without
         !     checkpointing them the first post-restart predictor reads
         !     the allocation zeros (a zero-thickness h_av) and the resume
         !     is not bit-exact.  Caught by `restart_bit_exact_*` the day
         !     pred_corr became the default.  `optional=.true.`: a pre-flip
         !     checkpoint has no such variables and must still resume --
         !     it will be an ssp_rk2 file, where they are unread.
         if (allocated(state%multilayer%u_av_layer)) then
            call reg%register_3d("ml_u_av_layer", state%multilayer%u_av_layer, 0, &
                                 size(state%multilayer%u_av_layer, 1), &
                                 size(state%multilayer%u_av_layer, 2), optional=.true.)
         end if
         if (allocated(state%multilayer%v_av_layer)) then
            call reg%register_3d("ml_v_av_layer", state%multilayer%v_av_layer, 0, &
                                 size(state%multilayer%v_av_layer, 1), &
                                 size(state%multilayer%v_av_layer, 2), optional=.true.)
         end if
         if (allocated(state%multilayer%h_av_layer)) then
            call reg%register_3d("ml_h_av_layer", state%multilayer%h_av_layer, 0, &
                                 size(state%multilayer%h_av_layer, 1), &
                                 size(state%multilayer%h_av_layer, 2), optional=.true.)
         end if
      end if

      ! --- rho_layer (review #2): diagnosed each stage, but NOT
      !     recomputed on dt_therm-skipped slow steps, so the value
      !     carried forward is thermo-step-dependent.  Checkpoint it. ---
      if (state%use_multilayer .and. allocated(state%multilayer%rho_layer)) then
         call register_full_3d(reg, "ml_rho_layer", state%multilayer%rho_layer)
      end if

      ! --- KPP lagged BL-depth seed (review #1): the V_t² wstar3_lagged
      !     term reads the previous step's bl_depth.  Allocated
      !     unconditionally; required for a bit-exact KPP resume. ---
      if (allocated(state%vmix%bl_depth)) then
         call reg%register_2d("vmix_bl_depth", state%vmix%bl_depth, 0, &
                              size(state%vmix%bl_depth, 1), size(state%vmix%bl_depth, 2))
      end if

      ! --- EPBL prev-MLD seed + kd_int (merged every stage) ---
      if (allocated(state%epbl%mld)) then
         call reg%register_2d("epbl_mld", state%epbl%mld, 0, &
                              size(state%epbl%mld, 1), size(state%epbl%mld, 2))
      end if
      if (allocated(state%epbl%kd_int)) then
         call register_full_3d(reg, "epbl_kd_int", state%epbl%kd_int)
      end if

      ! --- MEKE prognostic eddy-energy field (capability [5]): without
      !     this it would cold-reset to 0 each restart, losing the spun-up
      !     eddy field.  Allocated whenever the multilayer path is on. ---
      if (allocated(state%meke%meke)) then
         call reg%register_2d("meke", state%meke%meke, 0, &
                              size(state%meke%meke, 1), size(state%meke%meke, 2))
      end if

      ! --- Wet/dry hysteresis mask (docs/ocean_wetdry_plan.md): persistent
      !     front state (wet/dry/held-in-band).  optional — a restart
      !     written before the knob existed re-seeds from depth at
      !     configure instead of failing. ---
      if (allocated(state%dyn%bt_work%wd_wet_dyn)) then
         call reg%register_2d("wd_wet_dyn", state%dyn%bt_work%wd_wet_dyn, 0, &
                              size(state%dyn%bt_work%wd_wet_dyn, 1), &
                              size(state%dyn%bt_work%wd_wet_dyn, 2), optional=.true.)
      end if

      ! --- Sea-ice frazil bank (PR 1): un-spent supercooling heat the ice
      !     model (PR 3) will consume.  PERSISTENT (accumulates across
      !     steps), so a restart must carry it or banked energy would be
      !     silently discarded.  optional — only ice-on runs write it, and
      !     an ice-on resume from an older checkpoint re-seeds 0. ---
      if (allocated(state%ice%frazil_heat)) then
         call reg%register_2d("ice_frazil_heat", state%ice%frazil_heat, 0, &
                              size(state%ice%frazil_heat, 1), &
                              size(state%ice%frazil_heat, 2), optional=.true.)
      end if

      ! --- Sea-ice PR 3b: brine-rejection flux diag.  `salt_flux_diag` is
      !     REQUIRED for restart-exact Q_salt refill: `rdb_ice_ocean_coupler`
      !     rebuilds Q_salt from this field every thermo step, and the
      !     driver's configure-time resume fold re-applies it so the first
      !     post-resume window sees the same flux the uninterrupted run
      !     would have (Q_salt itself is configure-static and NOT
      !     registered).  `m_frozen_diag` is registered alongside for
      !     post-resume diagnostic continuity only.  optional — an
      !     ice-on resume from a pre-3b checkpoint re-seeds 0 (a one-window
      !     cold-start of the brine flux) rather than failing.
      !
      !     PR 3c extends the same contract to the melt side:
      !     `heat_flux_diag` is REQUIRED for restart-exact Q_heat refill
      !     (mirrors `salt_flux_diag` — `ice_ocean_heat_flux` rebuilds
      !     Q_heat from it every thermo step, and the driver resume fold
      !     re-applies it); `m_melt_diag` is diagnostic-continuity only,
      !     like `m_frozen_diag`.  PR 31 extends the same required-refill
      !     contract to `sw_thru_diag` — `ice_ocean_sw_flux` rebuilds the
      !     ocean shortwave from it every thermo step and the driver resume
      !     fold re-applies it, so it is registered alongside
      !     `heat_flux_diag`.  The coupleable atmospheric-forcing seam
      !     (`atm_sf0`/`atm_dsfdt`/`atm_sw_dn`) and the column-driver
      !     scratch (`fb`/`sst_seam`/`ssurf_seam`/`tfw_seam`/`tsurf_out`/
      !     `h2o_ocn_to_ice`/`h2o_ice_to_ocn`/`heat_to_ocn`/`sw_thru`) are
      !     ALL recomputed every thermo step — deliberately NOT
      !     registered (the per-category `sw_thru` stays unregistered; only
      !     its per-cell reduction `sw_thru_diag` is carried). ---
      if (allocated(state%ice%m_frozen_diag)) then
         call reg%register_2d("ice_m_frozen_diag", state%ice%m_frozen_diag, 0, &
                              size(state%ice%m_frozen_diag, 1), &
                              size(state%ice%m_frozen_diag, 2), optional=.true.)
      end if
      if (allocated(state%ice%salt_flux_diag)) then
         call reg%register_2d("ice_salt_flux_diag", state%ice%salt_flux_diag, 0, &
                              size(state%ice%salt_flux_diag, 1), &
                              size(state%ice%salt_flux_diag, 2), optional=.true.)
      end if
      if (allocated(state%ice%m_melt_diag)) then
         call reg%register_2d("ice_m_melt_diag", state%ice%m_melt_diag, 0, &
                              size(state%ice%m_melt_diag, 1), &
                              size(state%ice%m_melt_diag, 2), optional=.true.)
      end if
      if (allocated(state%ice%heat_flux_diag)) then
         call reg%register_2d("ice_heat_flux_diag", state%ice%heat_flux_diag, 0, &
                              size(state%ice%heat_flux_diag, 1), &
                              size(state%ice%heat_flux_diag, 2), optional=.true.)
      end if
      if (allocated(state%ice%sw_thru_diag)) then
         call reg%register_2d("ice_sw_thru_diag", state%ice%sw_thru_diag, 0, &
                              size(state%ice%sw_thru_diag, 1), &
                              size(state%ice%sw_thru_diag, 2), optional=.true.)
      end if

      ! --- Ice-shelf cavity basal melt (P2b): the TWO OWNED surface-flux
      !     components.  These ARE registered, and the exclusion rule
      !     above says exactly why: `Q_heat`/`Q_salt` are derived views
      !     and stay out, but "any filler that makes a COMPONENT
      !     time-varying MUST register it".  The melt rate is a function
      !     of the live state, so `heat_cavity`/`salt_cavity` are
      !     time-varying — and they are written at the END of outer step
      !     N and integrated on step N+1 (the ice coupler's documented
      !     one-step lag), so without them a warm restart would apply
      !     zero melt for its first step.  `optional=.true.`: an older
      !     checkpoint resumes with the zero seed rather than failing.
      !     The slot's own arrays (`cavity_flux%melt`, `t_b`, ...) are
      !     NOT registered — they are recomputed from state before first
      !     use, the same derived-field rule as `mass_flux_*`.
      if (allocated(state%surface_flux%heat_cavity)) then
         call reg%register_2d("sf_heat_cavity", state%surface_flux%heat_cavity, 0, &
                              size(state%surface_flux%heat_cavity, 1), &
                              size(state%surface_flux%heat_cavity, 2), optional=.true.)
      end if
      if (allocated(state%surface_flux%salt_cavity)) then
         call reg%register_2d("sf_salt_cavity", state%surface_flux%salt_cavity, 0, &
                              size(state%surface_flux%salt_cavity, 1), &
                              size(state%surface_flux%salt_cavity, 2), optional=.true.)
      end if

      ! --- Sea-ice PR 5: C-grid EVP dynamics prognostics.  u_ice/v_ice
      !     become PROGNOSTIC under dynamics=.true. (previously a v1
      !     interim sampler output — cheap to always register).
      !     str_d/str_t/str_s are the H&D stress tensor components;
      !     fxoc/fyoc are carried for DIAGNOSTIC CONTINUITY only (PR 63 —
      !     see below, they are no longer what makes tau_x/tau_y
      !     restart-exact).  All optional: a resume from a pre-5
      !     checkpoint re-seeds 0 (a one-window cold-start of the ice
      !     velocity/stress) rather than failing.  tau_a_x/tau_a_y are
      !     deliberately NOT registered (configure-time snapshot, rebuilt
      !     fresh from the wind-stress config on every run).
      !
      !     PR 63: tau_ocn_x/tau_ocn_y + tau_ocn_valid are what make the
      !     ice->ocean tau MEDIATION restart-exact.  ice_ocean_stress_flux
      !     mirrors its own blended output into tau_ocn_x/y every outer
      !     step; ice_ocean_stress_resume_apply COPIES them into
      !     surface_stress%tau_x/y at configure (after the wind re-seed),
      !     replacing the old deleted resume-fold RECONSTRUCT routine
      !     (which used the checkpoint's POST-thermo ci — wrong whenever a
      !     checkpoint step's thermo/transport changed ci after the blend,
      !     F4).  fxoc/fyoc above are no longer load-bearing for this: they
      !     are write-only-within-a-call accumulators
      !     (`ice_evp_dynamics` zeros them at entry), so the resume fold
      !     was their only cross-step-boundary reader — that reader is
      !     gone now.  tau_ocn_valid is register_scalar (device_mapped
      !     forced .false.) and gated on the SAME allocated(...) check as
      !     the arrays, so an ocean-only restart file never grows this
      !     variable. ---
      if (allocated(state%ice%u_ice)) then
         call reg%register_2d("ice_u_ice", state%ice%u_ice, 0, &
                              size(state%ice%u_ice, 1), &
                              size(state%ice%u_ice, 2), optional=.true.)
      end if
      if (allocated(state%ice%v_ice)) then
         call reg%register_2d("ice_v_ice", state%ice%v_ice, 0, &
                              size(state%ice%v_ice, 1), &
                              size(state%ice%v_ice, 2), optional=.true.)
      end if
      if (allocated(state%ice%str_d)) then
         call reg%register_2d("ice_str_d", state%ice%str_d, 0, &
                              size(state%ice%str_d, 1), &
                              size(state%ice%str_d, 2), optional=.true.)
      end if
      if (allocated(state%ice%str_t)) then
         call reg%register_2d("ice_str_t", state%ice%str_t, 0, &
                              size(state%ice%str_t, 1), &
                              size(state%ice%str_t, 2), optional=.true.)
      end if
      if (allocated(state%ice%str_s)) then
         call reg%register_2d("ice_str_s", state%ice%str_s, 0, &
                              size(state%ice%str_s, 1), &
                              size(state%ice%str_s, 2), optional=.true.)
      end if
      if (allocated(state%ice%fxoc)) then
         call reg%register_2d("ice_fxoc", state%ice%fxoc, 0, &
                              size(state%ice%fxoc, 1), &
                              size(state%ice%fxoc, 2), optional=.true.)
      end if
      if (allocated(state%ice%fyoc)) then
         call reg%register_2d("ice_fyoc", state%ice%fyoc, 0, &
                              size(state%ice%fyoc, 1), &
                              size(state%ice%fyoc, 2), optional=.true.)
      end if
      ! PR 63: the ice->ocean tau mediation seam (see the comment block
      ! above).  The allocated(...) gate is mandatory for all THREE
      ! entries, including the scalar — tau_ocn_valid lives on the type
      ! whether or not the ice slot is live; gating it on the array's
      ! allocation (rather than registering it unconditionally) keeps an
      ! ocean-only restart file free of a stray ice_tau_ocn_valid
      ! variable.
      if (allocated(state%ice%tau_ocn_x)) then
         call reg%register_2d("ice_tau_ocn_x", state%ice%tau_ocn_x, 0, &
                              size(state%ice%tau_ocn_x, 1), &
                              size(state%ice%tau_ocn_x, 2), optional=.true.)
         call reg%register_2d("ice_tau_ocn_y", state%ice%tau_ocn_y, 0, &
                              size(state%ice%tau_ocn_y, 1), &
                              size(state%ice%tau_ocn_y, 2), optional=.true.)
         call reg%register_scalar("ice_tau_ocn_valid", state%ice%tau_ocn_valid, &
                                  optional=.true.)
      end if

      ! --- Sea-ice Winton column prognostics (PR 3a): part_size, m_ice,
      !     m_snow are rank-3 (register directly); enth_ice/sal_ice
      !     (rank-4, per-k bottom-up) and enth_snow (rank-4, nk=1) are
      !     registered as per-k rank-3 slices — the registry is
      !     rank-2/3 only.  All optional: a resume from a pre-3a
      !     checkpoint warns-and-seeds the init values (part_size all
      !     open water, m_ice/m_snow/enth_* = 0, sal_ice =
      !     ICE_BULK_SALINITY) rather than failing. ---
      if (allocated(state%ice%part_size)) then
         call reg%register_3d("ice_part_size", state%ice%part_size, 0, &
                              size(state%ice%part_size, 1), &
                              size(state%ice%part_size, 2), optional=.true.)
      end if
      if (allocated(state%ice%m_ice)) then
         call reg%register_3d("ice_m_ice", state%ice%m_ice, 0, &
                              size(state%ice%m_ice, 1), &
                              size(state%ice%m_ice, 2), optional=.true.)
      end if
      if (allocated(state%ice%m_snow)) then
         call reg%register_3d("ice_m_snow", state%ice%m_snow, 0, &
                              size(state%ice%m_snow, 1), &
                              size(state%ice%m_snow, 2), optional=.true.)
      end if
      if (allocated(state%ice%enth_ice)) then
         do it = 1, size(state%ice%enth_ice, 4)
            call reg%register_3d("ice_enth_ice_k"//to_string(it), &
                                 state%ice%enth_ice(:, :, :, it), 0, &
                                 size(state%ice%enth_ice, 1), &
                                 size(state%ice%enth_ice, 2), optional=.true.)
         end do
      end if
      if (allocated(state%ice%sal_ice)) then
         do it = 1, size(state%ice%sal_ice, 4)
            call reg%register_3d("ice_sal_ice_k"//to_string(it), &
                                 state%ice%sal_ice(:, :, :, it), 0, &
                                 size(state%ice%sal_ice, 1), &
                                 size(state%ice%sal_ice, 2), optional=.true.)
         end do
      end if
      if (allocated(state%ice%enth_snow)) then
         do it = 1, size(state%ice%enth_snow, 4)
            call reg%register_3d("ice_enth_snow_k"//to_string(it), &
                                 state%ice%enth_snow(:, :, :, it), 0, &
                                 size(state%ice%enth_snow, 1), &
                                 size(state%ice%enth_snow, 2), optional=.true.)
         end do
      end if

      ! --- Kappa-shear kd_int / tke_int (merged every stage) ---
      if (allocated(state%kshear%kd_int)) then
         call register_full_3d(reg, "kshear_kd_int", state%kshear%kd_int)
      end if
      if (allocated(state%kshear%tke_int)) then
         call register_full_3d(reg, "kshear_tke_int", state%kshear%tke_int)
      end if

      ! --- Live persistent BC state (review #3) ---
      ! Chapman radiation corner scalars: host-side, always present.
      call reg%register_scalar("bc_eta_old_chapman_w", state%bc%eta_old_chapman_w)
      call reg%register_scalar("bc_eta_old_chapman_e", state%bc%eta_old_chapman_e)
      call reg%register_scalar("bc_eta_old_chapman_s", state%bc%eta_old_chapman_s)
      call reg%register_scalar("bc_eta_old_chapman_n", state%bc%eta_old_chapman_n)
      ! Tracer reservoirs: rank-3 (edge, layer, tracer); device-mapped
      ! when allocated, so they join the update-self walk by default.
      if (allocated(state%bc%tres_west)) then
         call register_full_3d(reg, "bc_tres_west", state%bc%tres_west)
      end if
      if (allocated(state%bc%tres_east)) then
         call register_full_3d(reg, "bc_tres_east", state%bc%tres_east)
      end if
      if (allocated(state%bc%tres_south)) then
         call register_full_3d(reg, "bc_tres_south", state%bc%tres_south)
      end if
      if (allocated(state%bc%tres_north)) then
         call register_full_3d(reg, "bc_tres_north", state%bc%tres_north)
      end if
   end subroutine ocean_state_build_restart_registry

   subroutine register_full_3d(reg, tag, arr)
      !! Register a rank-3 field as a FULL local array (interior +
      !! ghosts) — ng=0 over the total extent.  See the ghost-cell
      !! policy note in `ocean_state_build_restart_registry`.
      type(restart_registry_t), intent(inout) :: reg
      character(len=*), intent(in) :: tag
      real(wp), target, intent(in) :: arr(:, :, :)
      call reg%register_3d(tag, arr, 0, size(arr, 1), size(arr, 2))
   end subroutine register_full_3d

#ifndef RDB_NO_NETCDF
   subroutine ocean_state_fill_restart_metadata(state, grid, meta)
      !! Build the grid/vcoord/tracer fingerprint validated on resume
      !! (review #7).
      type(ocean_state_t), intent(in) :: state
      type(hgrid_t), intent(in) :: grid
      type(ocean_restart_metadata_t), intent(out) :: meta
      meta%nz_ml = state%multilayer%nz_ml
      meta%nghost = grid%nghost
      meta%vcoord_type = state%vcoord%coord_type
      meta%vcoord_name = vcoord_type_name(state%vcoord%coord_type)
      meta%idx_salinity = state%multilayer%idx_salinity
      meta%idx_temperature = state%multilayer%idx_temperature
      if (allocated(state%multilayer%tracers)) then
         meta%n_tracers = size(state%multilayer%tracers)
      else
         meta%n_tracers = 0
      end if
   end subroutine ocean_state_fill_restart_metadata
#endif

   pure function pgf_nonoverlap_gate_on(cfg, variant) result(on)
      !! Single source of truth for "is the grounded-layer PGF gate armed?"
      !! (`&ocean_isopycnal_nml pgf_skip_nonoverlap` →
      !! `ocean_pressure_force_t%skip_nonoverlap`).
      !!
      !! Two call sites must agree BEFORE any array is allocated:
      !! `ocean_state_init_from_config` latches the answer so the PGF
      !! slot's allocation gate can decide whether FV_MOM6 needs `z_centre`
      !! (it is read by nothing else on that path), and `configure_ocean_pgf`
      !! re-evaluates it after the full config pass and fails loud if the two
      !! disagree.  Keeping the predicate in one place is what makes that
      !! guard meaningful.
      !!
      !! ONLY under `VCOORD_LAGRANGIAN`: that is the coordinate where layers
      !! wedge out against the bed onto the `angstrom_h` floor while staying
      !! massive one cell away, which is what makes the face PGF ill-posed and
      !! leaves a spurious gradient AT REST.  `gprime` differences interface
      !! positions directly (no layer-centre Jacobian, no `z_centre` buffer),
      !! so it is the one variant excluded here — `configure_ocean_pgf` warns
      !! for it.  `mont` IS covered: the Montgomery potential is exact at rest
      !! wherever a layer is present on both sides of the face, but a grounded
      !! layer sits at the bed on one side and at its flat-isopycnal height on
      !! the other, so `M` stops being horizontally uniform and the residual
      !! returns.
      use rdb_constants, only: VCOORD_LAGRANGIAN
      type(config_t), intent(in) :: cfg
      integer, intent(in) :: variant
         !! `OPGF_VARIANT_*` code (already parsed from `cfg%ocean%pgf%form`).
      logical :: on
      on = cfg%ocean%isopycnal%pgf_skip_nonoverlap .and. &
           parse_ocean_vcoord_type(cfg%vcoord_type) == VCOORD_LAGRANGIAN .and. &
           (variant == OPGF_VARIANT_MONT .or. &
            variant == OPGF_VARIANT_FV_LITE .or. &
            variant == OPGF_VARIANT_FV_WRIGHT .or. &
            variant == OPGF_VARIANT_FV_MOM6)
   end function pgf_nonoverlap_gate_on

   pure function vcoord_type_name(code) result(name)
      !! Human-readable tag for a VCOORD_* enum (mismatch messages only;
      !! the integer enum is what is validated).
      use rdb_constants, only: VCOORD_SIGMA, VCOORD_ZSIGMA, VCOORD_ZSTAR, &
                               VCOORD_ZSTAR_SIGMA, VCOORD_ZSTAR_FULL, &
                               VCOORD_EULERIAN_Z, VCOORD_LAGRANGIAN
      integer, intent(in) :: code
      character(len=32) :: name
      select case (code)
      case (VCOORD_SIGMA)
         name = "sigma"
      case (VCOORD_ZSIGMA)
         name = "zsigma"
      case (VCOORD_ZSTAR)
         name = "zstar"
      case (VCOORD_ZSTAR_SIGMA)
         name = "zstar_sigma"
      case (VCOORD_ZSTAR_FULL)
         name = "zstar_full"
      case (VCOORD_EULERIAN_Z)
         name = "eulerian_z"
      case (VCOORD_LAGRANGIAN)
         name = "lagrangian"
      case default
         name = "code_"//to_string(code)
      end select
   end function vcoord_type_name

   subroutine ocean_state_restart_write(state, grid, decomp, filename, t, step, ierr)
      !! Build the registry, pull every DEVICE-MAPPED registered array
      !! down with `!$acc update self`, then write the per-rank file
      !! durably (tmp + rename).  Safe to call mid-run (state is
      !! device-resident) — the D->H pulls leave the device copy
      !! authoritative.  `intent(inout)`: the `update self` mutates the
      !! host copy of `state` (review #12), and pointer association into
      !! `state` requires the `target` attribute.
      !!
      !! `ierr` (P7, same treatment as the read side): absent behaves as
      !! today (`error stop` on any NetCDF/rename failure, or on a
      !! no-NetCDF build); present returns `OCEAN_STATUS_ERR_IO` instead.
      type(ocean_state_t), intent(inout), target :: state
      type(hgrid_t), intent(in) :: grid
      type(decomp_t), intent(in) :: decomp
      character(len=*), intent(in) :: filename
      real(wp), intent(in) :: t
      integer, intent(in) :: step
      integer, intent(out), optional :: ierr
#ifndef RDB_NO_NETCDF
      type(restart_registry_t) :: reg
      type(ocean_restart_metadata_t) :: meta
      integer :: e

      if (present(ierr)) ierr = OCEAN_STATUS_OK

      call ocean_state_build_restart_registry(state, grid, reg)
      call ocean_state_fill_restart_metadata(state, grid, meta)

      ! Device sync: pull each DEVICE-MAPPED registered array host-ward.
      ! Host-only entries (rank-0 scalars, never-mapped buffers) are
      ! skipped — `update self` on an unmapped array crashes on GPU
      ! (review #8).
      do e = 1, reg%n
         associate (en => reg%entries(e))
            if (.not. en%device_mapped) cycle
            if (en%rank == 2) then
               !$acc update self(en%p2)
            else if (en%rank == 3) then
               !$acc update self(en%p3)
            end if
         end associate
      end do

      call ocean_restart_write_local(filename, reg, decomp, meta, t, step, &
                                     state%dyn%outer_step_count, ierr=ierr)
#else
      call fail("ocean_state_restart_write: built without NetCDF (RDB_ENABLE_NETCDF=ON required)", &
                ierr, OCEAN_STATUS_ERR_IO)
#endif
   end subroutine ocean_state_restart_write

   subroutine ocean_state_restart_read(state, grid, decomp, filename, t, step, ierr)
      !! Read a per-rank ocean restart into the (host) prognostic arrays.
      !! MUST run BEFORE `ocean_state_enter_data` — the subsequent H->D
      !! copy carries the restored values to the device.  Restores
      !! `dyn%outer_step_count` (dt_therm alignment).
      !!
      !! Error policy (review #4): when the caller PASSES `ierr`, a
      !! decomp/schema/grid mismatch returns it non-zero (the test path).
      !! When the caller OMITS `ierr` (the production driver path), a
      !! mismatch is FATAL — `error stop` naming the kind + the file —
      !! so a mismatched resume dies loudly instead of silently
      !! cold-starting and then clobbering the good checkpoint at the next
      !! cadence.
      type(ocean_state_t), intent(inout), target :: state
      type(hgrid_t), intent(in) :: grid
      type(decomp_t), intent(in) :: decomp
      character(len=*), intent(in) :: filename
      real(wp), intent(out) :: t
      integer, intent(out) :: step
      integer, intent(out), optional :: ierr
#ifndef RDB_NO_NETCDF
      type(restart_registry_t) :: reg
      type(ocean_restart_metadata_t) :: meta
      integer :: osc, local_ierr

      call ocean_state_build_restart_registry(state, grid, reg)
      call ocean_state_fill_restart_metadata(state, grid, meta)
      call ocean_restart_read_local(filename, reg, decomp, meta, t, step, osc, &
                                    ierr=local_ierr)
      if (present(ierr)) then
         ! Translate the INTERNAL 1/2/3 restart_mismatch_kind convention
         ! (private to `ocean_restart_check_decomp`/`ocean_restart_read_local`)
         ! into a named `OCEAN_STATUS_ERR_RESTART_*` code before it reaches
         ! the caller — the raw integers collide with
         ! `OCEAN_STATUS_ERR_CONFIG_PARSE`/`_CONFIG_VALIDATE`/`_SETUP`
         ! (P0.1 review F4).
         if (local_ierr /= 0) then
            ierr = restart_mismatch_status_code(local_ierr)
         else
            ierr = OCEAN_STATUS_OK
         end if
      else if (local_ierr /= 0) then
         call fail("Ocean restart resume aborted: "// &
                   restart_mismatch_kind(local_ierr)//" mismatch reading "// &
                   trim(filename)//".  A mismatched resume would cold-start "// &
                   "and then overwrite the good checkpoint; refusing.", &
                   code=restart_mismatch_status_code(local_ierr))
      end if
      if (local_ierr == 0) state%dyn%outer_step_count = osc
#else
      t = 0.0_wp
      step = 0
      call fail("ocean restart requires RDB_ENABLE_NETCDF=ON", ierr, OCEAN_STATUS_ERR_IO)
      return
#endif
   end subroutine ocean_state_restart_read

   pure function restart_mismatch_kind(code) result(kind)
      !! Map the check-decomp ierr code (INTERNAL 1/2/3 convention —
      !! see `ocean_restart_check_decomp`) to a human-readable kind for
      !! the legacy `error stop` message.
      integer, intent(in) :: code
      character(len=:), allocatable :: kind
      select case (code)
      case (1)
         kind = "schema-version"
      case (2)
         kind = "decomposition"
      case (3)
         kind = "grid/vcoord/tracer"
      case default
         kind = "unknown"
      end select
   end function restart_mismatch_kind

   pure function restart_mismatch_status_code(code) result(status_code)
      !! Map the check-decomp ierr code (INTERNAL 1/2/3 convention) onto
      !! the PUBLIC, collision-free `OCEAN_STATUS_ERR_RESTART_*` codes
      !! returned through `ocean_state_restart_read`'s `ierr` (P0.1
      !! review F4).
      integer, intent(in) :: code
      integer :: status_code
      select case (code)
      case (1)
         status_code = OCEAN_STATUS_ERR_RESTART_SCHEMA
      case (2)
         status_code = OCEAN_STATUS_ERR_RESTART_DECOMP
      case (3)
         status_code = OCEAN_STATUS_ERR_RESTART_GRID
      case default
         status_code = OCEAN_STATUS_ERR_RESTART_GRID
      end select
   end function restart_mismatch_status_code

   subroutine apply_layer_rho_init(ms, layer_rho_init, nz_ml, ierr)
      !! Public only for the unit-test suite (no production module imports it);
      !! ignore when developing production code in other modules.
      !! Helper: write per-layer density from `layer_rho_init(:)` when the
      !! caller has set at least one non-sentinel value.  Sentinel is
      !! `< 0`; we expect exactly `nz_ml` non-sentinel entries listed
      !! in `k=1..nz_ml` order (bed → surface).
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: layer_rho_init(:)
      integer, intent(in) :: nz_ml
      integer, intent(out), optional :: ierr
         !! Non-zero when the non-sentinel count mismatches `nz_ml`, when
         !! present; absent behaves as today (`error stop`).
      integer :: i, j, k, n_set

      n_set = count(layer_rho_init >= 0.0_wp)
      if (n_set == 0) then
         if (present(ierr)) ierr = OCEAN_STATUS_OK
         return  ! sentinel everywhere — leave rho_layer alone
      end if

      if (n_set /= nz_ml) then
         call fail("ocean_layer_rho_init: count of non-sentinel "// &
                   "entries ("//to_string(n_set)//") does not match "// &
                   "nz_layers ("//to_string(nz_ml)//").", ierr, OCEAN_STATUS_ERR_IC_SEED)
         return
      end if

      do k = 1, nz_ml
         do j = 1, size(ms%rho_layer, 2)
            do i = 1, size(ms%rho_layer, 1)
               ms%rho_layer(i, j, k) = layer_rho_init(k)
            end do
         end do
      end do

      ! Root-only: this is a once-per-run init note, not per-rank state — on a
      ! multi-rank run it would otherwise print nranks identical copies.
      if (comm_env_rank() == 0) then
         call logger%info("ocean_layer_rho_init: wrote per-layer ρ directly, "// &
                          "bypassing EOS init.  ρ(k=1, bed) = "// &
                          to_string(layer_rho_init(1))//" kg/m³, "// &
                          "ρ(k="//to_string(nz_ml)//", surf) = "// &
                          to_string(layer_rho_init(nz_ml))//" kg/m³.")
      end if
      if (present(ierr)) ierr = OCEAN_STATUS_OK
   end subroutine apply_layer_rho_init

   pure function ocean_linear_layer_density(rho_lightest, rho_range, nz_ml) result(rho)
      !! Linearly-spaced layer densities for the MOM6 `COORD_CONFIG="linear"`
      !! IC analogue (`set_coord_linear`), re-derived for Roundabout's bottom-up
      !! layer convention.  Returns `nz_ml` densities, bed `k=1` heaviest →
      !! surface `k=nz_ml` lightest, layer-centred:
      !!
      !!     ρ(k) = rho_lightest + rho_range · (nz_ml − k + 0.5)/nz_ml
      !!
      !! so ρ(nz_ml) ≈ rho_lightest (surface, lightest) and ρ(1) ≈
      !! rho_lightest + rho_range (bed, heaviest), with uniform spacing
      !! rho_range/nz_ml between adjacent layers.  MOM6 indexes `k=1` at the
      !! lightest (top) layer; the `(nz_ml − k)` flip here gives the same
      !! physical column under Roundabout's `k=1`-is-bed orientation.
      real(wp), intent(in) :: rho_lightest   !! Surface (lightest) layer density [kg/m³].
      real(wp), intent(in) :: rho_range      !! Top-to-bottom density contrast [kg/m³].
      integer, intent(in) :: nz_ml           !! Number of layers.
      real(wp) :: rho(nz_ml)
      integer :: k

      do k = 1, nz_ml
         rho(k) = rho_lightest + rho_range* &
                  ((real(nz_ml - k, wp) + 0.5_wp)/real(nz_ml, wp))
      end do
   end function ocean_linear_layer_density

   subroutine seed_eady_ic(state, grid, cfg, ierr)
      !! Eady-front overlay IC.  Assumes flat bottom + uniform layer
      !! split already in place from `ocean_state_seed_from_cfg`.
      !!
      !! Sets:
      !!   T(i, j, k) = T_ref + dT_dz · z(k) + dT_dy · (y_phys(j) - y_mid)
      !!                                                       + ε(i, j, k)
      !!   u_face_x(i, j, k) = dU/dz · (z(k) - z_mid)
      !! with `dU/dz` from thermal-wind balance:
      !!     dU/dz = g · alpha_T / (rho_0 · f) · dT_dy
      !!
      !! Convention: k=1 is the bed, k=nz the surface.  z(k=1) =
      !! -H + dz/2, z(k=nz) = -dz/2.  y_phys uses cell centres.
      !!
      !! The perturbation ε is a uniform-random ±eady_pert_amp/2 noise
      !! seeded from `cfg%ocean%ic%eady_pert_seed`, applied to interior cells
      !! (j ∈ [ng+2, ng+ny_phys-1]) so wall-adjacent rows stay clean.
      !!
      !! Requires `coriolis_f /= 0`.  Errors if topo_config is not "flat".
      !!
      !! MPI: `y_mid` (the front centre) uses the GLOBAL meridional extent
      !! `grid%ny_global`, and the local row index is mapped to a global
      !! physical position with `grid%j_offset_global`.  Built from the
      !! LOCAL extents instead, every rank would put the front in the middle
      !! of its own tile.  On a single rank the offset is 0 and global ==
      !! local, so the seed is byte-identical.
      !!
      !! Not decomposition-invariant: the `eady_pert_amp` noise is drawn
      !! per-rank on the LOCAL array, so its realisation depends on the rank
      !! layout (the deterministic Eady profile underneath does not).
      type(ocean_state_t), intent(inout) :: state
      type(hgrid_t), intent(in) :: grid
      type(config_t), intent(in) :: cfg
      integer, intent(out), optional :: ierr
         !! Non-zero on an Eady-IC configuration conflict when present;
         !! absent behaves as today (`error stop`).

      integer :: i, j, k, nx_total, ny_total, nz_ml, ng
      integer :: ny_phys, j_phys, joff
      real(wp) :: H, dz, y_phys, y_mid, z_k, z_mid
      real(wp) :: dUdz, T_local, hk
      real(wp), allocatable :: noise(:, :, :)
      integer :: idx_T, n_seed, seed_val
      integer, allocatable :: seed_buf(:)

      if (trim(cfg%ocean%topo%topo_config) /= "flat") then
         call fail("seed_eady_ic: requires topo_config='flat'", ierr, OCEAN_STATUS_ERR_IC_SEED)
         return
      end if
      if (abs(cfg%coriolis_f) < tiny(1.0_wp)) then
         call fail("seed_eady_ic: requires coriolis_f /= 0", ierr, OCEAN_STATUS_ERR_IC_SEED)
         return
      end if

      nx_total = size(state%barotropic%b, 1)
      ny_total = size(state%barotropic%b, 2)
      nz_ml = state%multilayer%nz_ml
      ng = grid%nghost
      ! LOCAL ny_phys — used only to place the noise-free wall-adjacent
      ! rows, which is a per-rank array-bound question, not a coordinate.
      ny_phys = grid%ny_phys
      joff = grid%j_offset_global
      idx_T = state%multilayer%idx_temperature

      H = cfg%ocean%topo%max_depth
      dz = H/real(nz_ml, wp)
      y_mid = 0.5_wp*real(grid%ny_global, wp)*grid%dy
      z_mid = -0.5_wp*H

      ! Thermal-wind balance: ∂u_g/∂z = (g/(fρ₀))·∂ρ/∂y, with linear
      ! EOS ρ = ρ₀ − α(T−T_ref) giving ∂ρ/∂y = −α·∂T/∂y.  Therefore
      !     dU/dz = -g·α/(f·ρ₀) · dT/dy
      ! (cold-to-north dT/dy < 0 → eastward shear with z, surface jet
      !  westerly — the canonical mid-latitude convention).
      !
      ! Geostrophic balance for the y-independent part of -f·u would
      ! also demand a linearly tilted SSH (∂η/∂y = f·dU/dz·z_mid/g).
      ! Setting that η here AND refreshing h_layer accordingly produces
      ! ghost-row h_layer mismatch and tracer-mass corruption in the
      ! continuity-PPM path — the analytic balance doesn't exactly
      ! match the FV-lite PGF's discrete vertical integral.  Left for
      ! follow-up; for now the IC is hydrostatic-only and IG waves
      ! shed in the first ~24h as the system geostrophically adjusts.
      dUdz = -GRAVITY*cfg%ocean%ic%alpha_T/(cfg%ocean%ic%rho_0*cfg%coriolis_f) &
             *cfg%ocean%ic%eady_dT_dy

      ! Pre-generate noise on host (RNG isn't device-callable).  Fill
      ! the full array including ghosts; we zero ghost + bordering rows
      ! after to keep the perturbation cleanly inside the physical
      ! interior.
      allocate (noise(nx_total, ny_total, nz_ml))
      call random_seed(size=n_seed)
      allocate (seed_buf(n_seed))
      seed_val = cfg%ocean%ic%eady_pert_seed
      seed_buf = [(seed_val + 17*j, j=1, n_seed)]
      call random_seed(put=seed_buf)
      call random_number(noise)
      noise = cfg%ocean%ic%eady_pert_amp*(noise - 0.5_wp)
      ! Zero ghosts + first/last physical row to avoid wall-adjacent
      ! noise that would project onto the gravest mode.
      noise(:, 1:ng + 1, :) = 0.0_wp
      noise(:, ng + ny_phys:, :) = 0.0_wp

      ! Override T(i,j,k) with the Eady analytical profile + noise.
      if (idx_T > 0) then
         do k = 1, nz_ml
            z_k = -H + (real(k, wp) - 0.5_wp)*dz
            do j = 1, ny_total
               ! GLOBAL physical row index of this local row.
               j_phys = j - ng + joff
               y_phys = (real(j_phys, wp) - 0.5_wp)*grid%dy
               T_local = cfg%ocean%ic%eady_T_ref + cfg%ocean%ic%eady_dT_dz*z_k &
                         + cfg%ocean%ic%eady_dT_dy*(y_phys - y_mid)
               do i = 1, nx_total
                  hk = state%multilayer%h_layer(i, j, k)
                  state%multilayer%tracers(idx_T)%hTr(i, j, k) = &
                     (T_local + noise(i, j, k))*hk
               end do
            end do
         end do
      end if

      ! Thermal-wind-balanced u_face_x(z).  U depends on z only (no y, x).
      do k = 1, nz_ml
         z_k = -H + (real(k, wp) - 0.5_wp)*dz
         state%multilayer%u_face_x_layer(:, :, k) = dUdz*(z_k - z_mid)
         ! hu on x-faces: interpolate cell-centred h to the u-face (end faces
         ! take the adjacent cell), mirroring the production convention in
         ! face_depth_mean_u.  u_face/h_layer are staggered (nx+1 faces vs nx
         ! cells), so the old whole-array `u_face * h_layer` was non-conformant
         ! (silently OOB on non-checking compilers; LFortran flags it).
         associate (uf => state%multilayer%u_face_x_layer, &
                    hl => state%multilayer%h_layer, &
                    huf => state%multilayer%hu_face_x_layer)
            huf(1, :, k) = uf(1, :, k)*hl(1, :, k)
            huf(2:nx_total, :, k) = uf(2:nx_total, :, k) &
                                    *0.5_wp*(hl(1:nx_total - 1, :, k) + hl(2:nx_total, :, k))
            huf(nx_total + 1, :, k) = uf(nx_total + 1, :, k)*hl(nx_total, :, k)
         end associate
      end do

      deallocate (noise, seed_buf)
      if (present(ierr)) ierr = OCEAN_STATUS_OK
   end subroutine seed_eady_ic

   subroutine seed_baroclinic_jet_ic(state, grid, cfg, ierr)
      !! Two-layer reduced-gravity baroclinic-instability IC — a
      !! geostrophically-balanced tanh jet in the upper (surface) layer of
      !! a spherical re-entrant channel, plus a front-localised sech²
      !! meander seed that the instability grows.  Reproduces the `bc_inst`
      !! spec (SIM_DETAILS.md §5), mapped to Roundabout's bottom-up layer
      !! convention (k=1 bed, k=nz surface).  The spec's layers FLIP:
      !!
      !!   * spec upper layer 1 (the jet) → Roundabout k=2 (surface)
      !!   * spec lower layer 2 (at rest) → Roundabout k=1 (bed)
      !!
      !! With `y = R·(φ − φ₀)` (m) about the jet centre and
      !! `xm = R·cosφ₀·(λ − λ_west)` (m):
      !!
      !!   ξ(x,y) = Δξ·tanh(y/L) + A_pert·sech²(y/L)·cos(k_x·xm)
      !!   η(y)   = −(g'/g_FS)·ξ                       (free-surface signature)
      !!   h_layer(:,:,1) = H_bed  + ξ                 (lower, k=1)
      !!   h_layer(:,:,2) = H_surf + (η − ξ)           (upper, k=2)
      !!   u_face_x_layer(:,:,2) = (g'/f)·(Δξ/L)·sech²(y/L)  (base tanh only)
      !!   u_face_x_layer(:,:,1) = 0 ;  v_face_y_layer = 0
      !!
      !! `g'` and `g_FS` are read from the SAME reduced-gravity parameters
      !! the gprime PGF consumes (`cfg%ocean%pgf%gprime_gint` /
      !! `gprime_gfs`), so the IC balance and the PGF that maintains it use
      !! one source of truth.  `f = 2Ω·sin(φ)` is the full variable
      !! Coriolis at the u-face (T-row) latitude — NOT the scalar
      !! `coriolis_f`.  The geostrophic jet uses the BASE tanh only; the
      !! perturbation is an unbalanced interface displacement the
      !! instability feeds on.
      !!
      !! Requires `nz_layers == 2`, `&ocean_pgf_nml form="gprime"`, and a
      !! spherical grid.  Fills the FULL arrays incl. ghost rows: φ/λ are
      !! continued analytically past the physical cells, which for a
      !! periodic-x channel wraps exactly because `k_x·L_x = 2π·pert_nx`
      !! (integer `pert_nx`), so the x-ghost columns match the physical
      !! ones.  Wall-y ghosts continue the smooth tanh/sech² profile.
      !!
      !! Host-compute (deterministic, RNG-free): plain host loops writing
      !! `state%multilayer` BEFORE `enter_data`, per the setup-code
      !! convention (a `do concurrent` here would round-trip unmapped
      !! arrays through the device).
      type(ocean_state_t), intent(inout) :: state
      type(hgrid_t), intent(in) :: grid
      type(config_t), intent(in) :: cfg
      integer, intent(out), optional :: ierr
         !! Non-zero on a baroclinic-jet IC configuration conflict when
         !! present; absent behaves as today (`error stop`).
      !!
      !! MPI: the global physical-index offsets and the global domain
      !! extents come off `grid` (`i_offset_global` / `j_offset_global`,
      !! `nx_global` / `ny_global`).  All four are load-bearing: the
      !! jet-centre latitude, the perturbation wavenumber AND both
      !! index→coordinate maps must describe the WHOLE domain, or every
      !! rank reproduces the entire jet inside its own tile (with
      !! `pert_nx = 3`, three wavelengths per TILE instead of three per
      !! DOMAIN).  On a single rank the offsets are 0 and global == local,
      !! so the seed is byte-identical.

      integer :: i, j, k, nx_total, ny_total, nz_ml, ng
      integer :: idx_S, idx_T
      integer :: ioff, joff, nxg, nyg
      real(wp) :: gprime, g_fs, H_total, H_surf, H_bed
      real(wp) :: jet_L, dxi, a_pert, kx, phi0, cos_phi0, Lx
      real(wp) :: rad_earth, omega, lon_west, lat_south, dlon, dlat
      real(wp) :: phi, y_arg, tanh_y, sech2, f_u, u_jet
      real(wp) :: lam_deg, xm, xi_full, eta_local

      nz_ml = state%multilayer%nz_ml
      if (nz_ml /= 2) then
         call fail("seed_baroclinic_jet_ic: requires nz_layers == 2 "// &
                   "(two-layer reduced gravity); got "//to_string(nz_ml), ierr, OCEAN_STATUS_ERR_IC_SEED)
         return
      end if
      if (trim(cfg%ocean%pgf%form) /= "gprime") then
         call fail("seed_baroclinic_jet_ic: requires &ocean_pgf_nml "// &
                   "form='gprime'; got '"//trim(cfg%ocean%pgf%form)//"'", ierr, OCEAN_STATUS_ERR_IC_SEED)
         return
      end if
      if (trim(cfg%ocean%grid%grid_config) /= "spherical") then
         call fail("seed_baroclinic_jet_ic: requires &ocean_grid_nml "// &
                   "grid_config='spherical'; got '"// &
                   trim(cfg%ocean%grid%grid_config)//"'", ierr, OCEAN_STATUS_ERR_IC_SEED)
         return
      end if
      if (cfg%ocean%grid%omega <= 0.0_wp) then
         call fail("seed_baroclinic_jet_ic: requires &ocean_grid_nml "// &
                   "omega > 0 (planetary Coriolis f = 2*Omega*sin(phi))", ierr, OCEAN_STATUS_ERR_IC_SEED)
         return
      end if

      nx_total = grid%nx_total
      ny_total = grid%ny_total
      ng = grid%nghost

      ! Reduced-gravity parameters — read from the gprime PGF's own config
      ! (single source of truth for g' and the free-surface gravity g_FS).
      gprime = cfg%ocean%pgf%gprime_gint
      g_fs = cfg%ocean%pgf%gprime_gfs

      ! Rest layer split: upper (surface, k=2) H_surf; lower (bed, k=1) H_bed.
      H_total = cfg%ocean%topo%max_depth
      H_surf = cfg%ocean%ic%upper_layer_rest
      H_bed = H_total - H_surf
      if (H_surf <= 0.0_wp .or. H_bed <= 0.0_wp) then
         call fail("seed_baroclinic_jet_ic: upper_layer_rest ("// &
                   to_string(H_surf)//" m) must lie in (0, max_depth="// &
                   to_string(H_total)//" m)", ierr, OCEAN_STATUS_ERR_IC_SEED)
         return
      end if

      ioff = grid%i_offset_global
      joff = grid%j_offset_global
      nxg = grid%nx_global
      nyg = grid%ny_global

      jet_L = cfg%ocean%ic%jet_half_width
      dxi = cfg%ocean%ic%interface_amp
      a_pert = cfg%ocean%ic%pert_amp_frac*dxi

      ! Spherical geometry.  dx/dy are dlon/dlat in DEGREES for spherical.
      rad_earth = cfg%ocean%grid%rad_earth
      omega = cfg%ocean%grid%omega
      lon_west = cfg%ocean%grid%lon_west
      lat_south = cfg%ocean%grid%lat_south
      dlon = grid%dx
      dlat = grid%dy
      ! Jet-centre latitude = domain centre; φ = (lat_south + ny/2·dlat).
      phi0 = (lat_south + 0.5_wp*real(nyg, wp)*dlat)*DEG2RAD
      cos_phi0 = cos(phi0)

      ! Zonal perturbation wavenumber.  L_x = R·cosφ₀·(lon span, rad);
      ! integer pert_nx ⇒ k_x·L_x = 2π·pert_nx makes the x-continuation
      ! periodic (ghost columns match the physical wrap).
      Lx = rad_earth*cos_phi0*(real(nxg, wp)*dlon*DEG2RAD)
      kx = TWO_PI*real(cfg%ocean%ic%pert_nx, wp)/Lx

      ! ---- h_layer (cell centres), full array incl. ghosts ----
      do j = 1, ny_total
         phi = (lat_south + (real(j - ng + joff, wp) - 0.5_wp)*dlat)*DEG2RAD
         y_arg = rad_earth*(phi - phi0)/jet_L
         tanh_y = tanh(y_arg)
         sech2 = 1.0_wp/cosh(y_arg)**2
         do i = 1, nx_total
            lam_deg = lon_west + (real(i - ng + ioff, wp) - 0.5_wp)*dlon
            xm = rad_earth*cos_phi0*((lam_deg - lon_west)*DEG2RAD)
            xi_full = dxi*tanh_y + a_pert*sech2*cos(kx*xm)
            eta_local = -(gprime/g_fs)*xi_full
            state%multilayer%h_layer(i, j, 1) = H_bed + xi_full            ! lower (bed)
            state%multilayer%h_layer(i, j, 2) = H_surf + (eta_local - xi_full)  ! upper (surface)
         end do
      end do

      ! ---- balanced zonal jet on the upper-layer u-faces (row-only) ----
      ! u depends on the T-row latitude only (base tanh; NOT the pert).  The
      ! u-face shares its T-row's latitude, so f = 2Ω·sin(φ_row).  Lower
      ! layer + all v-faces stay at rest.
      state%multilayer%u_face_x_layer = 0.0_wp
      state%multilayer%v_face_y_layer = 0.0_wp
      do j = 1, ny_total
         phi = (lat_south + (real(j - ng + joff, wp) - 0.5_wp)*dlat)*DEG2RAD
         y_arg = rad_earth*(phi - phi0)/jet_L
         sech2 = 1.0_wp/cosh(y_arg)**2
         f_u = 2.0_wp*omega*sin(phi)
         u_jet = (gprime/f_u)*(dxi/jet_L)*sech2
         do i = 1, nx_total + 1
            state%multilayer%u_face_x_layer(i, j, 2) = u_jet
         end do
      end do

      ! ---- hu on x-faces: u · face-interpolated h (mirror seed_eady_ic) ----
      ! u_face/h_layer are staggered (nx+1 faces vs nx cells); end faces take
      ! the adjacent cell.  Lower layer u=0 ⇒ hu=0.  hv stays zero (v=0).
      state%multilayer%hv_face_y_layer = 0.0_wp
      do k = 1, nz_ml
         associate (uf => state%multilayer%u_face_x_layer, &
                    hl => state%multilayer%h_layer, &
                    huf => state%multilayer%hu_face_x_layer)
            huf(1, :, k) = uf(1, :, k)*hl(1, :, k)
            huf(2:nx_total, :, k) = uf(2:nx_total, :, k) &
                                    *0.5_wp*(hl(1:nx_total - 1, :, k) + hl(2:nx_total, :, k))
            huf(nx_total + 1, :, k) = uf(nx_total + 1, :, k)*hl(nx_total, :, k)
         end associate
      end do

      ! ---- re-seed S/T uniform against the new h_layer ----
      ! The default seed set hTr = const·h_layer against the OLD uniform
      ! split; refresh so hTr/h_layer stays at the configured scalars.
      idx_S = state%multilayer%idx_salinity
      idx_T = state%multilayer%idx_temperature
      if (idx_S > 0) then
         call seed_tracer_uniform_impl( &
            state%multilayer%tracers(idx_S)%hTr, &
            state%multilayer%h_layer, cfg%initial_salinity, nz_ml)
      end if
      if (idx_T > 0) then
         call seed_tracer_uniform_impl( &
            state%multilayer%tracers(idx_T)%hTr, &
            state%multilayer%h_layer, cfg%initial_temperature, nz_ml)
      end if

      if (comm_env_rank() == 0) then
         call logger%info("seed_baroclinic_jet_ic: 2-layer reduced-gravity jet — "// &
                          "L="//to_string(jet_L)//" m, dxi="//to_string(dxi)// &
                          " m, g'="//to_string(gprime)//" m/s^2, H_surf="// &
                          to_string(H_surf)//" m, H_bed="//to_string(H_bed)// &
                          " m, n_x="//to_string(cfg%ocean%ic%pert_nx))
      end if
      if (present(ierr)) ierr = OCEAN_STATUS_OK
   end subroutine seed_baroclinic_jet_ic

   subroutine seed_geostrophic_adjustment_ic(state, grid, cfg, ierr)
      !! Rossby's classic geostrophic-adjustment problem.  Overlays a
      !! Gaussian SSH bump on a flat-bottom, single-layer (barotropic)
      !! state at rest:
      !!     η(x, y) = A · exp(-r² / L²)
      !!     h(i, j) = b(i, j) + η(i, j)
      !!     u = v = 0
      !! where r is the radial distance from the bump centre.
      !!
      !! Expected evolution (analytical):
      !!   * Length-scale ratio L / Rd with Rd = √(gH)/f governs the
      !!     split between radiated IG-wave energy and retained
      !!     balanced geostrophic ring.
      !!     L >> Rd → most energy retained, surface stays bumped with
      !!                 a balanced anticyclonic flow around it.
      !!     L << Rd → most energy radiates, the bump flattens out.
      !!     L ~ Rd  → partial, ~50/50.
      !!   * Time-series of η at the centre oscillates with period
      !!     2π/f, damping toward the residual value.
      !!   * IG-wave fronts propagate outward at c = √(gH).
      !!
      !! Diagnostic interpretation:
      !!   * Bump amplitude growing in time → numerical instability
      !!     (wrong PGF sign, bad continuity).
      !!   * Oscillation period != 2π/f → wrong Coriolis coupling.
      !!   * Energy not conserved (after damping out the IG transient)
      !!     → spurious sources / sinks in the bt step or Coriolis-adv.
      !!
      !! Requires `topo_config = "flat"` so the bump sits on a
      !! uniform reference depth.  Single-layer is recommended (nz=1)
      !! to isolate the barotropic adjustment; multi-layer works but
      !! the bump distributes uniformly across layers.
      !!
      !! MPI: the default bump centre uses the GLOBAL extents
      !! `grid%nx_global` / `grid%ny_global`, and local cell indices are
      !! mapped to global physical positions with `grid%i_offset_global` /
      !! `grid%j_offset_global`.  Built from the LOCAL extents instead,
      !! every rank would drop a full-amplitude bump in the middle of its
      !! own tile.  On a single rank the offsets are 0 and global == local,
      !! so the seed is byte-identical.
      type(ocean_state_t), intent(inout) :: state
      type(hgrid_t), intent(in) :: grid
      type(config_t), intent(in) :: cfg
      integer, intent(out), optional :: ierr
         !! Non-zero on a geostrophic-adjustment IC configuration conflict
         !! when present; absent behaves as today (`error stop`).

      integer :: i, j, nx_total, ny_total, nz_ml, ng
      integer :: i_phys, j_phys, ioff, joff
      real(wp) :: x_phys, y_phys, x_c, y_c, r2, eta_local, inv_L2

      if (trim(cfg%ocean%topo%topo_config) /= "flat") then
         call fail("seed_geostrophic_adjustment_ic: requires topo_config='flat'", ierr, OCEAN_STATUS_ERR_IC_SEED)
         return
      end if
      if (cfg%ocean%ic%ga_length_scale <= 0.0_wp) then
         call fail("seed_geostrophic_adjustment_ic: ga_length_scale must be > 0", ierr, OCEAN_STATUS_ERR_IC_SEED)
         return
      end if

      nx_total = size(state%barotropic%b, 1)
      ny_total = size(state%barotropic%b, 2)
      nz_ml = state%multilayer%nz_ml
      ng = grid%nghost
      ioff = grid%i_offset_global
      joff = grid%j_offset_global

      ! Default centre = GLOBAL basin midpoint (negative cfg values).
      x_c = cfg%ocean%ic%ga_x_center
      y_c = cfg%ocean%ic%ga_y_center
      if (x_c < 0.0_wp) x_c = 0.5_wp*real(grid%nx_global, wp)*grid%dx
      if (y_c < 0.0_wp) y_c = 0.5_wp*real(grid%ny_global, wp)*grid%dy

      inv_L2 = 1.0_wp/(cfg%ocean%ic%ga_length_scale*cfg%ocean%ic%ga_length_scale)

      ! Override h_total = b + η(x, y) and redistribute the per-layer
      ! thickness evenly.  Velocities + face fluxes left at zero from
      ! the default seed.  Velocities are unchanged at zero, so hu/hv
      ! stay at zero too.
      do j = 1, ny_total
         ! GLOBAL physical indices of this local cell.
         j_phys = j - ng + joff
         y_phys = (real(j_phys, wp) - 0.5_wp)*grid%dy
         do i = 1, nx_total
            i_phys = i - ng + ioff
            x_phys = (real(i_phys, wp) - 0.5_wp)*grid%dx
            r2 = (x_phys - x_c)**2 + (y_phys - y_c)**2
            eta_local = cfg%ocean%ic%ga_eta_amp*exp(-r2*inv_L2)
            state%barotropic%h(i, j) = state%barotropic%b(i, j) + eta_local
         end do
      end do

      call seed_h_layer_uniform_impl(state%multilayer%h_layer, &
                                     state%barotropic%h, nz_ml)
      ! Re-seed S, T so hTr/h_layer stays at the configured uniform
      ! values regardless of the η bump.
      if (state%multilayer%idx_salinity > 0) then
         call seed_tracer_uniform_impl( &
            state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, &
            state%multilayer%h_layer, cfg%initial_salinity, nz_ml)
      end if
      if (state%multilayer%idx_temperature > 0) then
         call seed_tracer_uniform_impl( &
            state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
            state%multilayer%h_layer, cfg%initial_temperature, nz_ml)
      end if
      if (present(ierr)) ierr = OCEAN_STATUS_OK
   end subroutine seed_geostrophic_adjustment_ic

   pure subroutine seed_h_layer_uniform_impl(h_layer, b, nz_ml, apply_wetdry_floor)
      !! Even-split layer thickness per column.  Pulled into a flat-impl
      !! so the `do concurrent` body works on plain allocatables — the
      !! outer shim reaches the multilayer + barotropic slot components
      !! once on the host.
      !!
      !! `apply_wetdry_floor` (wet/dry v2): with dynamic wet/dry ON the
      !! static land cutoff moves to `-land_margin`, so the intertidal band
      !! `b ∈ [-land_margin, LAND_DEPTH_THRESHOLD)` seeds LIVE (wet_mask=1).
      !! An emerged rest-bed (`b < 0`, normal for real coastlines / file
      !! bathymetry with intertidal terrain) would then get NEGATIVE layer
      !! thickness `b/nz_ml`, and since `D = Σ h_layer` (derive_bt_from_layers)
      !! that negative column sum becomes a negative barotropic depth at IC —
      !! the positive-definite limiter PRESERVES D≥0, it does not REPAIR a
      !! negative start.  When the flag is set we floor each layer to
      !! `2·H_VANISHED` (never 0 ⇒ no 1/0), so `Σ h_layer ≥ nz·2·H_VANISHED >
      !! 0` ⇒ D≥0 and `bt_eta = Σ h_layer − bt_H_ref` is consistent, and the
      !! emerged column sits at a vanished ≈0 depth ready to flood (its
      !! `wd_wet_dyn` gate seeds DRY since Σ h_layer < dry_depth).
      !!
      !! Why `2·H_VANISHED` and not `H_VANISHED`: the ALE tracer remap drain
      !! (`ocean_remap_tracer_field`) gates the concentration recovery on
      !! `h_old > H_FLOOR` with H_FLOOR == H_VANISHED and a STRICT `>`, so a
      !! layer sitting at exactly `H_VANISHED` fails the gate ⇒ its
      !! concentration is zeroed ⇒ the seeded S/T content (`hTr = const·h`)
      !! is DESTROYED on the first regrid (every outer step on sigma at the
      !! default `dt_therm_ratio=1`).  Flooring at `2·H_VANISHED` clears the
      !! strict drain gate so the seeded tracer survives — the same
      !! inflation floor `rdb_ocean_vcoord` uses for its z*/z-blend layers.
      !! Deep columns (`b/nz_ml > 2·H_VANISHED`) are untouched — only the
      !! near-zero/emerged band is floored.  Flag absent/false ⇒ literal
      !! `b/nz_ml` (byte-identical to the pre-v2 path).
      ! assumed-shape-ok: init routine called once at startup; size(b,1/2) used
      ! to derive loop bounds (flat-impl over registry-dereferenced allocatables).
      real(wp), intent(inout) :: h_layer(:, :, :)
      real(wp), intent(in)    :: b(:, :)  ! assumed-shape-ok: init routine; size(b,1/2) derives loop bounds
      integer, intent(in)    :: nz_ml
      logical, intent(in), optional :: apply_wetdry_floor
         !! When present and .true., floor each seeded layer to `2·H_VANISHED`
         !! so an emerged intertidal rest-bed never seeds a negative column
         !! (and hence a negative barotropic depth) AND survives the strict
         !! `h_old > H_VANISHED` remap-drain gate.  Absent ⇒ byte-identical.
      integer :: i, j, k, nx, ny
      real(wp) :: inv_nz
      logical :: do_floor
      nx = size(b, 1)
      ny = size(b, 2)
      inv_nz = 1.0_wp/real(nz_ml, wp)
      do_floor = .false.
      if (present(apply_wetdry_floor)) do_floor = apply_wetdry_floor
      ! Plain host loop ON PURPOSE: this runs BEFORE enter_data, so a
      ! do concurrent here makes -stdpar=gpu implicitly round-trip the
      ! (unmapped) arrays through the device per loop — measured 36.9 s
      ! of ic_seed at 5.3M cells.  Setup-only code seeds on the host;
      ! enter_data then maps the seeded values once.
      if (do_floor) then
         do k = 1, nz_ml
            do j = 1, ny
               do i = 1, nx
                  h_layer(i, j, k) = max(b(i, j)*inv_nz, 2.0_wp*H_VANISHED)
               end do
            end do
         end do
      else
         do k = 1, nz_ml
            do j = 1, ny
               do i = 1, nx
                  h_layer(i, j, k) = b(i, j)*inv_nz
               end do
            end do
         end do
      end if
   end subroutine seed_h_layer_uniform_impl

   pure subroutine seed_h_layer_uniform_z_impl(h_layer, b, nz_ml, max_depth, angstrom_h)
      !! `thickness_config = "uniform_z"`: MOM6 `initialize_thickness_uniform`
      !! port.  Lays uniform **z** interfaces
      !! over the GLOBAL `max_depth`, clips them bottom-up against the local
      !! bathymetry, and collapses whatever will not fit to a minimum-thickness
      !! floor:
      !!
      !!     z_top_target(k) = -max_depth * (nz_ml - k) / nz_ml
      !!     h(k) = max(z_top_target(k) - z_bot(k), h_floor)
      !!
      !! walking `k = 1 -> nz_ml` from the bed up (Roundabout's bottom-up
      !! convention; MOM6 walks `k = nz -> 1` from its surface-first index, the
      !! same sweep).  The surface layer's target is `z = 0`, so it absorbs the
      !! remainder and the telescoping sum is EXACTLY `b(i,j)` — the floor
      !! injects nothing, unlike the runtime `angstrom_h` clamp.
      !!
      !! Why this exists.  Under a horizontally-uniform density stack (the
      !! `rho_lightest`/`rho_range` linear-coordinate IC, MOM6
      !! `COORD_CONFIG="linear"`) the layer index IS the density, so the layer
      !! interfaces ARE the isopycnals.  Seeding `b/nz_ml` (`"sigma"`) therefore
      !! makes every isopycnal follow the bathymetry, which stands the whole
      !! `rho_range` contrast up across each shelf break at t=0: on the MOM6
      !! double-gyre spoon that is ~1.9 kg/m3 between the 100 m rim and the
      !! 2000 m interior, i.e. g' ~ 0.018 m/s2 driving a ~1.3 m/s gravity
      !! current over layers only `edge_depth/nz_ml` thick.  `"uniform_z"`
      !! instead gives flat resting isopycnals with the sub-bathymetry layers
      !! collapsed against the bed, which is the layered-model resting state.
      !!
      !! Deep columns are unaffected where `b >= max_depth`: every layer lands
      !! on its z target at `max_depth/nz_ml`, and any excess `b - max_depth`
      !! is absorbed by the bed-most layer (`k = 1`).
      !!
      !! `h_floor = max(angstrom_h, 2*H_VANISHED)`.  The `angstrom_h` term ties
      !! the seed to the same floor the Lagrangian continuity update uses, so a
      !! layer seeded as collapsed reads as vanished to
      !! `isopycnal_vanish_tol()` (`max(angstrom_h, H_VANISHED)`) and is picked
      !! up by `reset_vanished_u` / `cfl_ignore_vanished` from step 1.  The
      !! `2*H_VANISHED` term keeps the floor strictly positive when
      !! `angstrom_h = 0` (knob off) so no seeded layer is ever exactly zero.
      !!
      !! Degenerate columns (`b <= nz_ml * h_floor` — including dry/negative
      !! land bed) cannot hold `nz_ml` floored layers, so they fall back to the
      !! floored even split, matching `seed_h_layer_uniform_impl`'s wet/dry
      !! branch.  Land columns are masked out downstream either way.
      ! assumed-shape-ok: init routine called once at startup; size(b,1/2) used
      ! to derive loop bounds (flat-impl over registry-dereferenced allocatables).
      real(wp), intent(inout) :: h_layer(:, :, :)
      real(wp), intent(in)    :: b(:, :)  ! assumed-shape-ok: init routine; size(b,1/2) derives loop bounds
      integer, intent(in)     :: nz_ml
      real(wp), intent(in)    :: max_depth
         !! Global basin depth the z interfaces are laid over (`ocean_max_depth`).
      real(wp), intent(in)    :: angstrom_h
         !! Isopycnal minimum-thickness floor (`&ocean_isopycnal_nml angstrom_h`).
      integer  :: i, j, k, nx, ny
      real(wp) :: inv_nz, h_floor, depth, z_bot, z_top, z_top_target

      nx = size(b, 1)
      ny = size(b, 2)
      inv_nz = 1.0_wp/real(nz_ml, wp)
      ! Honour `angstrom_h` all the way down.  The previous
      ! `max(angstrom_h, 2*H_VANISHED)` pinned the seeded floor at 3e-4 m,
      ! which made angstrom_h = 1e-6 and 1e-10 indistinguishable.  That matters:
      ! the vanished layers stack along the bed, so each of their interfaces
      ! carries the full topographic slope and contributes g'*db/dx of spurious
      ! PGF at REST.  The acceleration is independent of h, but the transport it
      ! drives scales WITH h -- which is why MOM6 can tolerate the same
      ! acceleration at ANGSTROM = 1e-10.  Fall back to 2*H_VANISHED only when
      ! the knob is off, so no layer is ever seeded at exactly zero.
      if (angstrom_h > 0.0_wp) then
         h_floor = angstrom_h
      else
         h_floor = 2.0_wp*H_VANISHED
      end if

      ! Plain host loop ON PURPOSE — same reasoning as
      ! `seed_h_layer_uniform_impl`: this runs BEFORE enter_data, so a
      ! `do concurrent` here would make -stdpar=gpu round-trip the
      ! (unmapped) arrays through the device once per loop.
      do j = 1, ny
         do i = 1, nx
            depth = b(i, j)
            if (depth <= real(nz_ml, wp)*h_floor) then
               ! Too shallow (or dry) to hold nz_ml floored layers.
               do k = 1, nz_ml
                  h_layer(i, j, k) = max(depth*inv_nz, h_floor)
               end do
            else
               z_bot = -depth
               do k = 1, nz_ml
                  z_top_target = -max_depth*real(nz_ml - k, wp)*inv_nz
                  ! Never thinner than the floor ...
                  z_top = max(z_top_target, z_bot + h_floor)
                  ! ... and never so thick that the (nz_ml - k) layers still
                  ! above it cannot each clear the floor below z = 0.  Without
                  ! this the sweep can run past the surface whenever
                  ! `max_depth/nz_ml` is not comfortably above `h_floor`, and
                  ! the telescoping sum stops equalling `b`.  The column guard
                  ! (`depth > nz_ml*h_floor`) makes the two bounds consistent
                  ! at every k, and pins `z_top = 0` exactly at k = nz_ml.
                  z_top = min(z_top, -real(nz_ml - k, wp)*h_floor)
                  h_layer(i, j, k) = z_top - z_bot
                  z_bot = z_top
               end do
            end if
         end do
      end do
   end subroutine seed_h_layer_uniform_z_impl

   pure subroutine seed_tracer_uniform_impl(hTr, h_layer, tracer_const, nz_ml)
      !! Seed `hTr = const * h_layer` for a uniform-IC tracer (salinity).
      ! assumed-shape-ok: init routine called once at startup; size(hTr,1/2)
      ! used to derive loop bounds (flat-impl over registry-dereferenced allocatables).
      real(wp), intent(inout) :: hTr(:, :, :)
      real(wp), intent(in)    :: h_layer(:, :, :)  ! assumed-shape-ok: init routine; size() derives loop bounds
      real(wp), intent(in)    :: tracer_const
      integer, intent(in)    :: nz_ml
      integer :: i, j, k, nx, ny
      nx = size(hTr, 1)
      ny = size(hTr, 2)
      ! Plain host loop ON PURPOSE: this runs BEFORE enter_data, so a
      ! do concurrent here makes -stdpar=gpu implicitly round-trip the
      ! (unmapped) arrays through the device per loop — measured 36.9 s
      ! of ic_seed at 5.3M cells.  Setup-only code seeds on the host;
      ! enter_data then maps the seeded values once.
      do k = 1, nz_ml
         do j = 1, ny
            do i = 1, nx
               hTr(i, j, k) = tracer_const*h_layer(i, j, k)
            end do
         end do
      end do
   end subroutine seed_tracer_uniform_impl

   pure subroutine seed_tracer_stratified_impl(hTr, h_layer, t_layer, nz_ml)
      !! Seed `hTr(i,j,k) = t_layer(k) * h_layer(i,j,k)` for the
      !! per-layer stratified-IC tracer (temperature).  `t_layer(:)` is
      !! a small 1D array pre-computed on host.
      ! assumed-shape-ok: init routine called once at startup; size(hTr,1/2)
      ! used to derive loop bounds (flat-impl over registry-dereferenced allocatables).
      real(wp), intent(inout) :: hTr(:, :, :)
      real(wp), intent(in)    :: h_layer(:, :, :)  ! assumed-shape-ok: init routine; size() derives loop bounds
      real(wp), intent(in)    :: t_layer(:)
      integer, intent(in)    :: nz_ml
      integer :: i, j, k, nx, ny
      nx = size(hTr, 1)
      ny = size(hTr, 2)
      ! Plain host loop ON PURPOSE: this runs BEFORE enter_data, so a
      ! do concurrent here makes -stdpar=gpu implicitly round-trip the
      ! (unmapped) arrays through the device per loop — measured 36.9 s
      ! of ic_seed at 5.3M cells.  Setup-only code seeds on the host;
      ! enter_data then maps the seeded values once.
      do k = 1, nz_ml
         do j = 1, ny
            do i = 1, nx
               hTr(i, j, k) = t_layer(k)*h_layer(i, j, k)
            end do
         end do
      end do
   end subroutine seed_tracer_stratified_impl

   pure subroutine seed_wet_mask_impl(wet_mask, b, land_cutoff)
      !! 1.0 where `b >= cutoff` (ocean), 0.0 elsewhere (land).
      !! Default (absent `land_cutoff`): `cutoff = LAND_DEPTH_THRESHOLD`
      !! (byte-identical to the pre-v2 path).  When `land_cutoff` is
      !! present, a column is land iff `b < land_cutoff`; the wetdry-aware
      !! seed passes `land_cutoff = -land_margin` so intertidal columns
      !! (bed above rest MSL but below the flood headroom) stay wet_mask=1
      !! and the dynamic wd_wet_dyn gate handles their wetting/drying.
      ! assumed-shape-ok: init routine called once at startup; size(b,1/2) used
      ! to derive loop bounds.
      real(wp), intent(inout) :: wet_mask(:, :)
      real(wp), intent(in)    :: b(:, :)  ! assumed-shape-ok: init routine; size(b,1/2) derives loop bounds
      real(wp), intent(in), optional :: land_cutoff
         !! When present, a column is land iff `b < land_cutoff` (used by
         !! the wetdry-aware seed: land_cutoff = -land_margin lets intertidal
         !! columns near/above rest MSL stay wet_mask=1). Absent ⇒ the default
         !! `b >= LAND_DEPTH_THRESHOLD` ocean test (byte-identical).
      integer :: i, j, nx, ny
      real(wp) :: cutoff
      nx = size(b, 1)
      ny = size(b, 2)
      cutoff = LAND_DEPTH_THRESHOLD
      if (present(land_cutoff)) cutoff = land_cutoff
      ! Plain host loop ON PURPOSE: this runs BEFORE enter_data, so a
      ! do concurrent here makes -stdpar=gpu implicitly round-trip the
      ! (unmapped) arrays through the device per loop — measured 36.9 s
      ! of ic_seed at 5.3M cells.  Setup-only code seeds on the host;
      ! enter_data then maps the seeded values once.
      do j = 1, ny
         do i = 1, nx
            if (b(i, j) >= cutoff) then
               wet_mask(i, j) = 1.0_wp
            else
               wet_mask(i, j) = 0.0_wp
            end if
         end do
      end do
   end subroutine seed_wet_mask_impl

   subroutine ocean_state_seed_land_cells(state, grid)
      !! Hold land T-cells (`wet_mask==0`) at FINITE reference values so
      !! the masked dyn-core never evaluates `0*NaN` (a zeroed face metric
      !! times a NaN land contribution is still NaN).  Spec §13.3 /
      !! `land_mask_final_resolution.md` R-land-state:
      !!   * `h_layer` floored to `H_VANISHED` (never 0 — avoids 1/0);
      !!   * tracers held at their seeded per-layer value (`hTr/h_old`),
      !!     re-scaled onto the floored thickness so `T/S = hTr/h` is
      !!     preserved and finite;
      !!   * layer + barotropic face velocities zeroed at land faces.
      !!
      !! Runs at SETUP, after `metrics_apply_land_mask` has derived
      !! `wet_u/wet_v`, BEFORE `ocean_state_enter_data`.  Plain host loops
      !! (pre-enter_data).  No-op for an all-wet domain (`wet_mask≡1`).
      type(ocean_state_t), intent(inout) :: state
      type(hgrid_t), intent(in) :: grid

      integer :: nz_ml, t

      associate (ms => state%multilayer)
         nz_ml = size(ms%h_layer, 3)

         ! Tracers FIRST: read the old (pre-floor) h to recover T/S, then
         ! re-scale hTr onto the floored thickness.  Must precede the
         ! h_layer floor so the per-layer value is recovered from the
         ! thickness it was seeded against.
         if (allocated(ms%tracers)) then
            do t = 1, size(ms%tracers)
               call seed_land_tracer_hold_impl(ms%tracers(t)%hTr, ms%h_layer, &
                                               ms%wet_mask, nz_ml)
            end do
         end if

         ! h_layer floor on land.
         call seed_land_h_floor_impl(ms%h_layer, ms%wet_mask, nz_ml)

         ! Zero layer face velocities at land faces.
         call seed_land_face_vel_impl(ms%u_face_x_layer, state%metrics%wet_u, nz_ml)
         call seed_land_face_vel_impl(ms%v_face_y_layer, state%metrics%wet_v, nz_ml)
      end associate

      ! Zero the layer mass-transports + barotropic face velocities at land
      ! faces (the transports ride masked metrics anyway, but resetting them
      ! keeps the held land state internally consistent).
      call seed_land_face_vel_impl(state%multilayer%hu_face_x_layer, &
                                   state%metrics%wet_u, nz_ml)
      call seed_land_face_vel_impl(state%multilayer%hv_face_y_layer, &
                                   state%metrics%wet_v, nz_ml)
      call seed_land_face_vel_2d_impl(state%barotropic%u_face_x, state%metrics%wet_u)
      call seed_land_face_vel_2d_impl(state%barotropic%v_face_y, state%metrics%wet_v)
      call seed_land_face_vel_2d_impl(state%barotropic%hu_face_x, state%metrics%wet_u)
      call seed_land_face_vel_2d_impl(state%barotropic%hv_face_y, state%metrics%wet_v)

      ! `grid` is part of the setup contract (shapes come from the state
      ! allocations); not re-read here.
      associate (unused => grid%nx_total)
      end associate
   end subroutine ocean_state_seed_land_cells

   pure subroutine seed_land_tracer_hold_impl(hTr, h_layer, wet_mask, nz)
      !! On land T-cells, recover the seeded per-layer tracer value
      !! `val = hTr/h_old` and re-scale onto the `H_VANISHED` floor so
      !! `hTr = val*H_VANISHED` stays finite and `T/S` is held.  Wet cells
      !! untouched (bit-identical when `wet_mask≡1`).
      integer, intent(in) :: nz
      real(wp), intent(inout) :: hTr(:, :, :)
      real(wp), intent(in) :: h_layer(:, :, :)
      real(wp), intent(in) :: wet_mask(:, :)
      integer :: i, j, k, nx, ny
      real(wp) :: val
      nx = size(hTr, 1)
      ny = size(hTr, 2)
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               if (wet_mask(i, j) == 0.0_wp) then
                  val = hTr(i, j, k)/max(h_layer(i, j, k), H_VANISHED)
                  hTr(i, j, k) = val*H_VANISHED
               end if
            end do
         end do
      end do
   end subroutine seed_land_tracer_hold_impl

   pure subroutine seed_land_h_floor_impl(h_layer, wet_mask, nz)
      !! Floor land-cell layer thickness to `H_VANISHED` (never 0 ⇒ no
      !! 1/0 in any per-layer divide).  Wet cells untouched.
      integer, intent(in) :: nz
      real(wp), intent(inout) :: h_layer(:, :, :)
      real(wp), intent(in) :: wet_mask(:, :)
      integer :: i, j, k, nx, ny
      nx = size(h_layer, 1)
      ny = size(h_layer, 2)
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               if (wet_mask(i, j) == 0.0_wp) h_layer(i, j, k) = H_VANISHED
            end do
         end do
      end do
   end subroutine seed_land_h_floor_impl

   pure subroutine seed_land_face_vel_impl(fld, wet_face, nz)
      !! Zero a per-layer face field (`(nf1,nf2,nz)`) at land faces
      !! (`wet_face==0`).  `wet_face` is the matching `wet_u`/`wet_v`.
      integer, intent(in) :: nz
      real(wp), intent(inout) :: fld(:, :, :)
      real(wp), intent(in) :: wet_face(:, :)
      integer :: i, j, k, nf1, nf2
      nf1 = size(fld, 1)
      nf2 = size(fld, 2)
      do k = 1, nz
         do j = 1, nf2
            do i = 1, nf1
               fld(i, j, k) = fld(i, j, k)*wet_face(i, j)
            end do
         end do
      end do
   end subroutine seed_land_face_vel_impl

   pure subroutine seed_land_face_vel_2d_impl(fld, wet_face)
      !! 2D (barotropic) variant of `seed_land_face_vel_impl`.
      real(wp), intent(inout) :: fld(:, :)
      real(wp), intent(in) :: wet_face(:, :)
      integer :: i, j, nf1, nf2
      nf1 = size(fld, 1)
      nf2 = size(fld, 2)
      do j = 1, nf2
         do i = 1, nf1
            fld(i, j) = fld(i, j)*wet_face(i, j)
         end do
      end do
   end subroutine seed_land_face_vel_2d_impl

   subroutine set_bathymetry_island(b, grid, max_depth, half_frac)
      !! Flat-bottom basin at `max_depth` everywhere, except a central
      !! square of LAND (`b = 0`, well below `LAND_DEPTH_THRESHOLD`).  The
      !! land square is centred on the physical domain and spans the
      !! central `2*half_frac` fraction of each axis (e.g. `half_frac=0.2`
      !! ⇒ the middle 40 % is land).  `half_frac` is taken from
      !! `&ocean_topo_nml slope_scale` at the call site (no new knob).
      !!
      !! Fills the FULL array including ghost rows (the formula-bathy
      !! ghost-fill gotcha): a ghost `b=0` at a wall-adjacent face sends
      !! the EOS into its ρ=ρ₀ fallback → spurious density jump → blowup.
      !! Ghost rows here inherit the flat `max_depth` (the land square is
      !! strictly interior), so every ghost is ocean.
      !!
      !! MPI: the global physical-index offsets and the global domain
      !! extents come off `grid` (`i_offset_global` / `j_offset_global`,
      !! `nx_global` / `ny_global`), so the land square is centred on the
      !! WHOLE domain and each rank carves only the part of it that falls
      !! inside its own tile.  On a single rank the offsets are 0 and
      !! global == local, so the fill is byte-identical.
      real(wp), intent(inout) :: b(:, :)
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: max_depth, half_frac
      integer :: i, j, i_phys, j_phys, ng, ci0, ci1, cj0, cj1, half_i, half_j
      integer :: ioff, joff, nxg, nyg
      real(wp) :: hf

      ioff = grid%i_offset_global
      joff = grid%j_offset_global
      nxg = grid%nx_global
      nyg = grid%ny_global

      ng = grid%nghost
      hf = half_frac
      if (hf <= 0.0_wp) hf = 0.2_wp     ! sensible default if knob left 0
      if (hf > 0.49_wp) hf = 0.49_wp    ! keep at least one wet ring inside

      ! Land-square physical-index bounds (centred on the GLOBAL domain), inclusive.
      half_i = nint(hf*real(nxg, wp))
      half_j = nint(hf*real(nyg, wp))
      ci0 = nxg/2 - half_i + 1
      ci1 = nxg/2 + half_i
      cj0 = nyg/2 - half_j + 1
      cj1 = nyg/2 + half_j

      ! Flat basin everywhere (incl. ghosts), then carve the interior land.
      ! The global physical index of local cell (i,j) is (i_phys + ioff).
      b = max_depth
      do j = 1, size(b, 2)
         j_phys = j - ng
         do i = 1, size(b, 1)
            i_phys = i - ng
            if ((i_phys + ioff) >= ci0 .and. (i_phys + ioff) <= ci1 .and. &
                (j_phys + joff) >= cj0 .and. (j_phys + joff) <= cj1) then
               b(i, j) = 0.0_wp   ! LAND (< LAND_DEPTH_THRESHOLD)
            end if
         end do
      end do
   end subroutine set_bathymetry_island

   subroutine set_bathymetry_double_drake(b, grid, max_depth, half_frac)
      !! "Double Drake" idealised supercontinent (Ferreira, Marshall &
      !! Campin 2010, *J. Climate*): a flat-bottom global ocean at
      !! `max_depth` with TWO thin meridional wall-continents 90° of
      !! longitude apart, each running from the north pole down to a
      !! southern-channel latitude, leaving a reentrant circumpolar
      !! channel (Drake-Passage analogue) to the south.  Used as the
      !! static-land-mask integration showcase on a spherical periodic-x
      !! sector grid.
      !!
      !! Index-space geometry (mapped to the Ferreira config by the
      !! caller's grid lon/lat extent): physical-x is the 360° periodic
      !! longitude; the two walls sit at the western seam (`i_phys=1`)
      !! and a quarter of the way across (`i_phys≈nx_phys/4`, ≡ 90°
      !! apart on a 360° domain).  Each wall is one cell wide and spans
      !! the NORTHERN `(1-2*half_frac)` fraction of the physical-y axis
      !! (north = pole), leaving the southern `2*half_frac` band fully
      !! open as the reentrant channel.  `half_frac` is taken from
      !! `&ocean_topo_nml slope_scale` (no new knob), default 0.2
      !! ⇒ southern 40 % of the basin is the open channel.
      !!
      !! Fills the FULL array incl. ghosts (formula-bathy ghost-fill
      !! gotcha): the walls extend through the y-ghost rows so the
      !! periodic-x seam wall is land in the halo too; the southern
      !! channel + interior basin stay ocean.
      !!
      !! MPI: the global physical-index offsets and the global domain
      !! extents come off `grid` (`i_offset_global` / `j_offset_global`,
      !! `nx_global` / `ny_global`), so the two wall meridians and the
      !! channel latitude are located on the WHOLE domain and each rank
      !! carves only the part inside its own tile.  On a single rank the
      !! offsets are 0 and global == local, so the fill is byte-identical.
      real(wp), intent(inout) :: b(:, :)
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: max_depth, half_frac
      integer :: i, j, i_phys, j_phys, ng, iw1, iw2, j_chan_top
      integer :: ioff, joff, nxg, nyg
      real(wp) :: hf

      ioff = grid%i_offset_global
      joff = grid%j_offset_global
      nxg = grid%nx_global
      nyg = grid%ny_global

      ng = grid%nghost
      hf = half_frac
      if (hf <= 0.0_wp) hf = 0.2_wp
      if (hf > 0.49_wp) hf = 0.49_wp

      ! Two one-cell-wide meridional walls, 90° apart on the GLOBAL 360° axis.
      iw1 = 1                   ! western-seam global wall at global i_phys = 1
      iw2 = max(2, nxg/4)       ! a quarter across ≡ 90° of longitude globally
      ! Channel occupies the southern `2*hf` fraction of the GLOBAL domain;
      ! walls cover the rest (the northern band up to the pole).
      j_chan_top = nint(2.0_wp*hf*real(nyg, wp))

      b = max_depth
      do j = 1, size(b, 2)
         j_phys = j - ng
         do i = 1, size(b, 1)
            i_phys = i - ng
            ! Wall cells: on either meridional line AND north of the
            ! reentrant channel.  Compare GLOBAL physical indices.
            ! j-ghosts north of the channel are wall (walls reach the pole),
            ! south ghosts stay ocean (open channel).
            if (((i_phys + ioff) == iw1 .or. (i_phys + ioff) == iw2) .and. &
                (j_phys + joff) > j_chan_top) then
               b(i, j) = 0.0_wp   ! LAND (< LAND_DEPTH_THRESHOLD)
            end if
         end do
      end do
   end subroutine set_bathymetry_double_drake

   subroutine set_bathymetry_spoon(b, grid, max_depth, edge_depth, slope_scale)
      !! Public only for the unit-test suite (no production module imports it);
      !! ignore when developing production code in other modules.
      !! Fill `b(:,:)` with the MOM6 spoon bathymetry.  In Cartesian
      !! terms (we collapse MOM6's lat/lon factors of `R_earth · π / 180`
      !! into a direct meters scale), the local depth is
      !!
      !!     D(i,j) = D_edge
      !!            + D_0 · sin(π · x_phys / x_len)
      !!                  · (1 − exp((y_phys − y_len) / slope_scale))
      !!
      !! with `D_0 = (max_depth − D_edge) / (1 − exp(−0.5 · y_len /
      !! slope_scale))²`.  Zero at east/west walls, full depth at the
      !! south wall, exponentially decaying to `D_edge` at the north
      !! wall.  Physical interior only; ghost cells retain whatever
      !! value `b` had on entry.
      !!
      !! MPI: the global physical-index offsets and the global domain
      !! extents come off `grid` (`i_offset_global` / `j_offset_global`,
      !! `nx_global` / `ny_global`), so `x_len` / `y_len` and the
      !! index→position map describe the WHOLE domain on every rank.  On a
      !! single rank the offsets are 0 and global == local, so the fill is
      !! byte-identical to the undecomposed formula.
      real(wp), intent(inout) :: b(:, :)
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: max_depth, edge_depth, slope_scale
      real(wp), parameter :: PI = 4.0_wp*atan(1.0_wp)
      real(wp) :: x_len, y_len, x_phys, y_phys, D_0, denom, D_local
      integer :: i, j, i_phys, j_phys, ng
      integer :: ioff, joff

      ioff = grid%i_offset_global
      joff = grid%j_offset_global

      ng = grid%nghost
      x_len = real(grid%nx_global, wp)*grid%dx
      y_len = real(grid%ny_global, wp)*grid%dy
      denom = 1.0_wp - exp(-0.5_wp*y_len/slope_scale)
      D_0 = (max_depth - edge_depth)/(denom*denom)

      ! Fill the full array including ghost rows.  Ghost cells left at
      ! the alloc-time zero cause `h_layer = 0` there, which sends the
      ! EOS into its vanishing-layer fallback (`rho_layer = rho_0`).
      ! That puts a non-physical density jump at every wall-adjacent
      ! face whenever the IC has `T_init /= T_ref` or `S_init /= S_ref`,
      ! and the BPG operator picks it up as a spurious horizontal
      ! pressure gradient.  Confirmed via the seamount test.
      do j = 1, size(b, 2)
         j_phys = j - ng
         ! Global physical y position: local j_phys shifted by joff.
         y_phys = (real(j_phys + joff, wp) - 0.5_wp)*grid%dy
         do i = 1, size(b, 1)
            i_phys = i - ng
            ! Global physical x position: local i_phys shifted by ioff.
            x_phys = (real(i_phys + ioff, wp) - 0.5_wp)*grid%dx
            D_local = edge_depth &
                      + D_0*sin(PI*x_phys/x_len) &
                      *(1.0_wp - exp((y_phys - y_len)/slope_scale))
            ! Safety: the sin term goes negative outside [0, x_len],
            ! and the exp term saturates the depth at the north wall.
            ! Clamp to [edge_depth, max_depth].
            if (D_local < edge_depth) D_local = edge_depth
            if (D_local > max_depth) D_local = max_depth
            b(i, j) = D_local
         end do
      end do
   end subroutine set_bathymetry_spoon

   subroutine seed_cavity_draft(state, grid, cfg, ierr)
      !! Fill `metrics%z_draft` (and its `cover_frac` companion) from
      !! `&ocean_cavity_dyn_nml`, then cross-validate the geometry against
      !! the seeded bathymetry.  Runs from `ocean_state_seed_from_cfg`
      !! IMMEDIATELY after the bathymetry and BEFORE the layer split and
      !! the wet-mask seed, which both read `b − z_draft`.
      !!
      !! Non-`pure` on purpose (the only cavity routine that is): it
      !! reports the grounded / over-land column counts through the
      !! logger and fails loud through the error ring.  The arithmetic it
      !! drives is in `rdb_ocean_cavity`, where every routine IS `pure`.
      !!
      !! UNITS.  The namelist carries metres; the setters work in GRID
      !! coordinate units (metres on Cartesian, DEGREES on
      !! spherical/curvilinear), so the box corners are converted with
      !! `topo_length_to_grid_units` and the dimensionless slope is
      !! converted the inverse way — the same trap `&ocean_topo_nml
      !! slope_scale` documents, where a metres length against a degrees
      !! position collapsed a seamount to a flat basin.
      type(ocean_state_t), intent(inout) :: state
      type(hgrid_t), intent(in) :: grid
      type(config_t), intent(in) :: cfg
      integer, intent(out) :: ierr

      integer :: draft_code, source_code, nx, ny, ng
      integer :: n_over_land, n_grounded, n_interior
      real(wp) :: per_metre, x0_g, x1_g, y0_g, y1_g, slope_g, amp
      real(wp) :: grounded_frac

      ierr = OCEAN_STATUS_OK
      nx = size(state%metrics%z_draft, 1)
      ny = size(state%metrics%z_draft, 2)
      ng = grid%nghost
      if (nx /= size(state%barotropic%b, 1) .or. ny /= size(state%barotropic%b, 2)) then
         call fail("&ocean_cavity_dyn_nml: z_draft is at its placeholder size — "// &
                   "metrics%use_cavity must be latched BEFORE metrics%init (it is "// &
                   "latched in ocean_state_init_from_config)", &
                   ierr, OCEAN_STATUS_ERR_IC_SEED)
         return
      end if

      draft_code = parse_cavity_draft_config(cfg%ocean%cavity_dyn%draft_config)
      source_code = parse_cavity_draft_source(cfg%ocean%cavity_dyn%draft_source)
      ! `validate_config` already refuses every spelling outside the v1
      ! envelope with a full explanation; this is the same gate one level
      ! down, for direct (test / API) callers that bypass it.
      if (draft_code /= CAVITY_DRAFT_NONE .and. draft_code /= CAVITY_DRAFT_FLAT &
          .and. draft_code /= CAVITY_DRAFT_LINEAR) then
         call fail("&ocean_cavity_dyn_nml draft_config='"// &
                   trim(cfg%ocean%cavity_dyn%draft_config)//"' is not available "// &
                   "(none|flat|linear; 'file' needs the static-2-D reader)", &
                   ierr, OCEAN_STATUS_ERR_IC_SEED)
         return
      end if
      if (source_code /= CAVITY_SOURCE_DRAFT .and. source_code /= CAVITY_SOURCE_THICKNESS) then
         call fail("&ocean_cavity_dyn_nml draft_source='"// &
                   trim(cfg%ocean%cavity_dyn%draft_source)//"' is not available "// &
                   "(draft|thickness; 'in_situ' isostasy is deferred)", &
                   ierr, OCEAN_STATUS_ERR_IC_SEED)
         return
      end if

      ! GRID units per metre (1 on Cartesian; degrees-per-metre otherwise).
      ! The "no limit" sentinels pass through UNCONVERTED so they stay
      ! sentinels on every grid.
      per_metre = topo_length_to_grid_units(1.0_wp, cfg%ocean%grid%grid_config, &
                                            cfg%ocean%grid%rad_earth)
      x0_g = cavity_bound_to_grid(cfg%ocean%cavity_dyn%draft_x0, per_metre)
      x1_g = cavity_bound_to_grid(cfg%ocean%cavity_dyn%draft_x1, per_metre)
      y0_g = cavity_bound_to_grid(cfg%ocean%cavity_dyn%draft_y0, per_metre)
      y1_g = cavity_bound_to_grid(cfg%ocean%cavity_dyn%draft_y1, per_metre)
      ! `draft_slope` is d(draft [m]) / d(x [m]); the setter wants
      ! d(draft [m]) / d(x [grid units]) = slope / (grid units per metre).
      slope_g = cfg%ocean%cavity_dyn%draft_slope/per_metre

      ! `draft_source = "thickness"`: the formula amplitude is an ice
      ! THICKNESS, converted by the Boussinesq-isostatic (flotation)
      ! relation `z_draft = rho_ice*h_ice/rho_0`.  Scaling the amplitude
      ! (and the slope with it) is exact because both setters are LINEAR
      ! in the amplitude — and it keeps one draft field, so everything
      ! downstream stays source-agnostic.
      amp = cfg%ocean%cavity_dyn%draft_depth
      if (source_code == CAVITY_SOURCE_THICKNESS) then
         amp = amp*cfg%ocean%cavity_dyn%rho_ice/state%eos%rho0
         slope_g = slope_g*cfg%ocean%cavity_dyn%rho_ice/state%eos%rho0
      end if

      select case (draft_code)
      case (CAVITY_DRAFT_FLAT)
         call set_draft_flat(state%metrics%z_draft, grid, amp, x0_g, x1_g, y0_g, y1_g)
      case (CAVITY_DRAFT_LINEAR)
         call set_draft_linear(state%metrics%z_draft, grid, amp, slope_g, &
                               x0_g, x1_g, y0_g, y1_g)
      case default  ! CAVITY_DRAFT_NONE
         state%metrics%z_draft = 0.0_wp
      end select

      ! Guard the field itself before anything derives geometry from it.
      ! Written as `.not. (z >= 0)` inside the helper so a NaN FAILS
      ! rather than sliding through a `z < 0` test that is false for NaN.
      if (.not. cavity_draft_is_finite_nonneg(state%metrics%z_draft, nx, ny)) then
         call fail("&ocean_cavity_dyn_nml: z_draft must be finite and >= 0 "// &
                   "everywhere (it is a DEPTH below z = 0, positive down)", &
                   ierr, OCEAN_STATUS_ERR_IC_SEED)
         return
      end if

      ! NO ICE OVER LAND: zero the draft wherever the bathymetry already
      ! says land, so land columns keep the datum they always had
      ! (`bt_H_ref = b`) and the counted-once invariant stays exact on
      ! every column.
      call cavity_apply_land_exclusion(state%metrics%z_draft, state%barotropic%b, &
                                       nx, ny, n_over_land)

      call cavity_count_grounded(state%barotropic%b, state%metrics%z_draft, &
                                 cfg%ocean%cavity_dyn%h_min_cavity, ng, &
                                 grid%nx_phys, grid%ny_phys, nx, ny, &
                                 n_grounded, n_interior)
      grounded_frac = real(n_grounded, wp)/real(max(n_interior, 1), wp)
      if (grounded_frac > cfg%ocean%cavity_dyn%grounded_max_frac) then
         call fail("&ocean_cavity_dyn_nml: the prescribed draft grounds "// &
                   to_string(n_grounded)//" of "//to_string(n_interior)// &
                   " interior columns ("//to_string(grounded_frac)//"), above "// &
                   "grounded_max_frac = "// &
                   to_string(cfg%ocean%cavity_dyn%grounded_max_frac)// &
                   ".  Either the draft is too deep for this bathymetry or the "// &
                   "shelf box is in the wrong place (check the UNITS of "// &
                   "draft_x0/x1 — they are metres, converted to grid units).", &
                   ierr, OCEAN_STATUS_ERR_IC_SEED)
         return
      end if

      call cavity_fill_cover_frac(state%metrics%cover_frac, state%metrics%z_draft, nx, ny)

      if (comm_env_rank() == 0) then
         call logger%info("Ice-shelf cavity: ON  draft_config='"// &
                          trim(cfg%ocean%cavity_dyn%draft_config)//"' source='"// &
                          trim(cfg%ocean%cavity_dyn%draft_source)//"' max draft = "// &
                          to_string(maxval(state%metrics%z_draft))//" m, "// &
                          "h_min_cavity = "// &
                          to_string(cfg%ocean%cavity_dyn%h_min_cavity)//" m")
         call logger%info("                  datum bt_H_ref = b - z_draft; "// &
                          to_string(n_grounded)//" of "//to_string(n_interior)// &
                          " interior columns grounded (-> LAND via the wet mask)")
         if (n_over_land > 0) then
            call logger%warning("&ocean_cavity_dyn_nml: draft zeroed on "// &
                                to_string(n_over_land)//" column(s) whose bed is "// &
                                "already land (b < LAND_DEPTH_THRESHOLD) — no ice "// &
                                "over land.  Check the shelf box if that is a surprise.")
         end if
      end if
   end subroutine seed_cavity_draft

   pure function cavity_bound_to_grid(bound_m, per_metre) result(bound_grid)
      !! Convert ONE shelf-box bound from metres to grid coordinate units,
      !! leaving the `CAVITY_BOUND_INF` "no limit" sentinel alone.  Without
      !! the guard the sentinel would be scaled by the degrees-per-metre
      !! factor on a spherical grid and come out as a finite (if absurd)
      !! bound — harmless numerically, but it would stop meaning what it
      !! says, and the next reader would have to re-derive that.
      real(wp), intent(in) :: bound_m
      real(wp), intent(in) :: per_metre
         !! Grid units per metre (1 on Cartesian).
      real(wp) :: bound_grid
      if (abs(bound_m) >= CAVITY_BOUND_INF) then
         bound_grid = bound_m
      else
         bound_grid = bound_m*per_metre
      end if
   end function cavity_bound_to_grid

   pure function topo_length_to_grid_units(length_m, grid_config, rad_earth) result(len_grid)
      !! Convert a metres length scale (`slope_scale` / `half_width`) into the
      !! GRID coordinate units the formula bathymetry setters operate in.
      !! Cartesian grids carry positions in metres, so the scale passes
      !! through unchanged.  Spherical / supergrid / tripolar grids carry
      !! positions in DEGREES, so the metres scale is converted to degrees of
      !! latitude (meridional metres-per-degree = `rad_earth · π/180`).  This
      !! keeps the seamount/spoon length scale commensurate with the grid
      !! coordinate; without it a metres scale divided by a degrees position
      !! underflows the Gaussian/exponential and the basin collapses flat.
      !! Zonal cells are narrower by `cos(lat)`, so a degree-isotropic bump is
      !! mildly elongated zonally in physical space — acceptable for these
      !! idealised topographies; the degeneracy is what this fixes.
      real(wp), intent(in) :: length_m, rad_earth
      character(len=*), intent(in) :: grid_config
      real(wp) :: len_grid
      select case (trim(grid_config))
      case ("spherical", "supergrid", "tripolar")
         len_grid = length_m/(rad_earth*DEG2RAD)
      case default
         len_grid = length_m
      end select
   end function topo_length_to_grid_units

   subroutine set_bathymetry_seamount(b, grid, max_depth, peak_depth, half_width)
      !! Public only for the unit-test suite (no production module imports it);
      !! ignore when developing production code in other modules.
      !! Fill `b(:,:)` with a centred Gaussian seamount bathymetry:
      !!
      !!     D(i,j) = max_depth − (max_depth − peak_depth)
      !!                          · exp(−((x − xc)² + (y − yc)²) / L²)
      !!
      !! with `(xc, yc)` at the centre of the physical domain and
      !! `L = half_width` (gaussian e-folding distance) — expressed in the
      !! GRID coordinate units (metres on Cartesian, degrees on spherical/
      !! curvilinear).  The production dispatch converts the metres
      !! `slope_scale` knob via `topo_length_to_grid_units`; callers that
      !! pass `half_width` directly must match the grid's units.
      !! Reaches `peak_depth` at the bump centre, asymptotes to
      !! `max_depth` far from the bump.
      !!
      !! Designed as a Tier-1.5 bridge test between flat-bottom
      !! analytical setups and full real-bathymetry regional runs.
      !! Tests σ-coord pressure-gradient and Coriolis-advection
      !! kernels under a varying h_layer without the confounds of
      !! a 50× depth ratio + IC stratification + closed-basin
      !! resonance that `tasman_2km.nml` introduces.
      !!
      !! Standard test protocol: uniform T,S (no APE), zero
      !! Coriolis or f-plane only, no wind, walls all sides.  The
      !! expected steady state is **identically zero motion** — any
      !! non-zero u, v, η is a σ-coord-related numerical artefact.
      !!
      !! MPI: the global physical-index offsets and the global domain
      !! extents come off `grid` (`i_offset_global` / `j_offset_global`,
      !! `nx_global` / `ny_global`), so `(xc, yc)` sits at the centre of the
      !! WHOLE domain on every rank rather than the centre of each tile.
      !! On a single rank the offsets are 0 and global == local, so the fill
      !! is byte-identical to the undecomposed formula.
      real(wp), intent(inout) :: b(:, :)
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: max_depth, peak_depth, half_width
      real(wp) :: x_len, y_len, xc, yc, x_phys, y_phys, r2, depression
      integer :: i, j, i_phys, j_phys, ng, nx_total, ny_total
      integer :: ioff, joff

      ioff = grid%i_offset_global
      joff = grid%j_offset_global

      ng = grid%nghost
      x_len = real(grid%nx_global, wp)*grid%dx
      y_len = real(grid%ny_global, wp)*grid%dy
      xc = 0.5_wp*x_len
      yc = 0.5_wp*y_len
      depression = max_depth - peak_depth
      nx_total = size(b, 1)
      ny_total = size(b, 2)

      ! Fill the full array including ghost rows so the EOS, BPG, and
      ! mass-flux operators see a smooth bathymetry on every face they
      ! touch.  Leaving ghosts at zero injects a spurious density jump
      ! at wall-adjacent faces (EOS falls back to `rho_0` where
      ! `h_layer = 0`) — the smoking-gun bug the seamount test was built
      ! to expose.
      do j = 1, ny_total
         j_phys = j - ng
         ! Global physical y position: local j_phys shifted by joff.
         y_phys = (real(j_phys + joff, wp) - 0.5_wp)*grid%dy
         do i = 1, nx_total
            i_phys = i - ng
            ! Global physical x position: local i_phys shifted by ioff.
            x_phys = (real(i_phys + ioff, wp) - 0.5_wp)*grid%dx
            r2 = (x_phys - xc)**2 + (y_phys - yc)**2
            b(i, j) = max_depth - depression*exp(-r2/(half_width*half_width))
         end do
      end do
   end subroutine set_bathymetry_seamount

   pure function nw2_cosbell(x, L) result(c)
      !! Cosine-bell kernel for the Neverworld2 basin: `0.5·(1 + cos(π·min(|x/L|,1)))`.
      !! Peaks at 1 for x=0, decays smoothly to 0 at |x|=L.  Re-derived from
      !! Marques et al. (2022, GMD) "Neverworld2"; MOM6-inspired.
      real(wp), intent(in) :: x, L
      real(wp) :: c
      real(wp), parameter :: PI = 4.0_wp*atan(1.0_wp)
      c = 0.5_wp*(1.0_wp + cos(PI*min(abs(x/L), 1.0_wp)))
   end function nw2_cosbell

   pure function nw2_spike(x, L) result(s)
      !! Sin-spike kernel for the Neverworld2 basin: `1 − sin(π·min(|x/L|,0.5))`.
      !! Equals 1 at x=0 and drops to 0 at |x|=L/2 (the 0.5 cap stops the sine
      !! re-ascending past its first zero).  Re-derived from Marques et al.
      !! (2022, GMD); MOM6-inspired.
      real(wp), intent(in) :: x, L
      real(wp) :: s
      real(wp), parameter :: PI = 4.0_wp*atan(1.0_wp)
      s = 1.0_wp - sin(PI*min(abs(x/L), 0.5_wp))
   end function nw2_spike

   subroutine set_bathymetry_neverworld2(b, grid, max_depth, nl_continent_amp, nl_roughness_amp, min_depth)
      !! Public only for the unit-test suite (no production module imports it);
      !! ignore when developing production code in other modules.
      !!
      !! Fill `b(:,:)` with the **Neverworld2** idealized-basin bathymetry
      !! (Marques et al. 2022, GMD; MOM6-inspired).  A single Pangaea-style
      !! basin spanning a 60°×140° spherical sector with a re-entrant
      !! (periodic-x) southern channel — the Drake-Passage analog.  Depth is
      !! built in **normalized coordinates**
      !!
      !!     x = (i_phys − 0.5)/nx_phys ∈ [0,1],   y = (j_phys − 0.5)/ny_phys ∈ [0,1]
      !!
      !! which equal MOM6's `(lon − west)/len_lon` and `(lat − south)/len_lat`
      !! on a uniform grid, so no metrics access is needed.  The fractional
      !! depth is
      !!
      !!     D_frac(x,y) = 1
      !!        − 1.1·spike(y−1, 0.12)              ! great northern wall
      !!        − 1.1·spike(y,   0.12)              ! Antarctica (south wall)
      !!        − A_c·[ continents + ridges ]        ! A_c = nl_continent_amp
      !!        − A_r·cos(14πx)·sin(14πy)            ! A_r = nl_roughness_amp
      !!        − A_r·cos(20πx)·cos(20πy)
      !!
      !! clamped `D_frac = max(D_frac, 0)`, then `D = D_frac·max_depth`.  The
      !! continent/ridge block carves two meridional barriers (S-America at the
      !! x-edges, Africa mid-basin) that wall the northern gyres, plus the
      !! Drake/Scotia ridge system that sets the channel sill; `A_c=0` gives an
      !! aquaplanet with the southern channel only (useful for bring-up).
      !!
      !! Sign convention: `b` is bottom depth positive-down (matching the spoon
      !! + flat branches).  Land emerges where D_frac→0 (the wet/dry mask is
      !! seeded downstream from `b ≥ LAND_DEPTH_THRESHOLD`).
      !!
      !! GHOST-ROW FOOTGUN (load-bearing): the loop runs over the FULL array
      !! incl. ghost rows.  Leaving ghosts at the alloc-time zero gives
      !! `h_layer = 0` there, sending the EOS into its vanishing-layer fallback
      !! (`rho_layer = rho_0`) → a spurious density jump at wall-adjacent faces
      !! → ~12-h e-fold blowup (the same trap the spoon/seamount setters guard).
      !! Host only — call `enter_data` afterwards.
      !!
      !! MPI: the normalised coordinates divide by the GLOBAL extents
      !! `grid%nx_global` / `grid%ny_global`, and local cell indices are
      !! mapped to global physical indices with `grid%i_offset_global` /
      !! `grid%j_offset_global`, so `x`/`y` sweep [0,1] across the WHOLE
      !! basin.  Normalised by the LOCAL extents instead, every rank would
      !! rebuild the entire Pangaea basin — walls, continents and all —
      !! inside its own tile.  On a single rank the offsets are 0 and
      !! global == local, so the fill is byte-identical.
      real(wp), intent(inout) :: b(:, :)
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: max_depth, nl_continent_amp, nl_roughness_amp
      real(wp), intent(in) :: min_depth
         !! Floor on the resulting depth (m).  The fractional-depth formula
         !! produces b in [0, ~1.1*max_depth]; the `1.1*spike` walls and the
         !! continent terms drive b to 0 (land) at the basin edges.  The ocean
         !! C-grid dyn-core does not yet carry a robust wet/dry path, and true
         !! zero-depth land + the resulting sub-metre surface layers blow up
         !! within ~12 h.  Flooring every cell to `min_depth` (MOM6's
         !! MINIMUM_DEPTH approach) turns the continents into shallow shelves
         !! that steer — rather than hard-block — the flow, keeping the basin
         !! all-wet and stable.  Lower values give stronger topographic
         !! steering at the cost of thinner layers (min_depth/nz); true land
         !! barriers wait on the wet/dry-plumbing follow-up.
      real(wp), parameter :: PI = 4.0_wp*atan(1.0_wp)
      real(wp) :: x, y, d_frac, nxp, nyp
      integer :: i, j, i_phys, j_phys, ng, ioff, joff

      ng = grid%nghost
      ioff = grid%i_offset_global
      joff = grid%j_offset_global
      nxp = real(grid%nx_global, wp)
      nyp = real(grid%ny_global, wp)

      do j = 1, size(b, 2)
         ! GLOBAL physical indices of this local cell.
         j_phys = j - ng + joff
         y = (real(j_phys, wp) - 0.5_wp)/nyp
         do i = 1, size(b, 1)
            i_phys = i - ng + ioff
            x = (real(i_phys, wp) - 0.5_wp)/nxp
            d_frac = 1.0_wp &
                     - 1.1_wp*nw2_spike(y - 1.0_wp, 0.12_wp) &   ! great northern wall
                     - 1.1_wp*nw2_spike(y, 0.12_wp) &            ! Antarctica (south wall)
                     - nl_continent_amp*( &
                     (1.2_wp*nw2_spike(x, 0.2_wp) + 1.2_wp*nw2_spike(x - 1.0_wp, 0.2_wp)) &
                     *nw2_spike(min(0.0_wp, y - 0.3_wp), 0.2_wp) &                          ! South America
                     + 1.2_wp*nw2_spike(x - 0.5_wp, 0.2_wp) &
                     *nw2_spike(min(0.0_wp, y - 0.55_wp), 0.2_wp) &                         ! Africa
                     + 1.2_wp*(nw2_spike(x, 0.12_wp) + nw2_spike(x - 1.0_wp, 0.12_wp)) &
                     *nw2_spike(max(0.0_wp, y - 0.06_wp), 0.12_wp) &                        ! Antarctic Peninsula
                     + 0.1_wp*(nw2_cosbell(x, 0.1_wp) + nw2_cosbell(x - 1.0_wp, 0.1_wp)) &  ! Drake Passage ridge
                     + 0.5_wp*nw2_cosbell(x - 0.16_wp, 0.05_wp)*(nw2_cosbell(y - 0.18_wp, 0.13_wp)**0.4_wp) &  ! Scotia Arc E
                     + 0.4_wp*(nw2_cosbell(x - 0.09_wp, 0.08_wp)**0.4_wp)*nw2_cosbell(y - 0.26_wp, 0.05_wp) &  ! Scotia Arc N
                     + 0.4_wp*(nw2_cosbell(x - 0.08_wp, 0.08_wp)**0.4_wp)*nw2_cosbell(y - 0.1_wp, 0.05_wp)) &  ! Scotia Arc S
                     - nl_roughness_amp*cos(14.0_wp*PI*x)*sin(14.0_wp*PI*y) &  ! roughness
                     - nl_roughness_amp*cos(20.0_wp*PI*x)*cos(20.0_wp*PI*y)
            if (d_frac < 0.0_wp) d_frac = 0.0_wp
            b(i, j) = max(d_frac*max_depth, min_depth)
         end do
      end do
   end subroutine set_bathymetry_neverworld2

end module rdb_ocean_state
