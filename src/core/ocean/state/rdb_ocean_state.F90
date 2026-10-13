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
                            DEG2RAD, TWO_PI, VCOORD_Z_FIXED, VCOORD_ZSTAR_FULL, &
                            VCOORD_ZSTAR
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
                                       ocean_bc_state_exit_data, &
                                       ocean_bc_type_from_string, OBC_PERIODIC, &
                                       OBC_TRIPOLAR_FOLD
   use rdb_ocean_periodic, only: ocean_periodic_wrap_centre_2d
   use rdb_ocean_fold, only: fold_north_centre
   use rdb_ocean_sponge, only: ocean_sponge_t
   use rdb_ocean_diag, only: ocean_diag_t
   use rdb_ocean_restart, only: ocean_restart_t, restart_registry_t
#ifndef RDB_NO_NETCDF
   use rdb_ocean_restart_io, only: ocean_restart_write_local, ocean_restart_read_local, &
                                   ocean_restart_metadata_t
   use rdb_ocean_data_input, only: ocean_data_input_t, ocean_data_input_load_static_2d
   use rdb_ocean_data_forcing, only: ocean_data_forcing_t
   use rdb_atm_state, only: atm_state_t
   use rdb_ocean_bulk_forcing, only: bulk_forcing_t
#endif
   use rdb_decomp, only: decomp_t
   use rdb_ocean_vcoord, only: ocean_vcoord_t, parse_ocean_vcoord_type, &
                               ocean_vcoord_z_fixed_target, &
                               ocean_vcoord_z_fixed_target_uniform, &
                               ocean_vcoord_eta0_target
   use rdb_ocean_metrics, only: ocean_metrics_t
   use rdb_ocean_cavity, only: parse_cavity_draft_config, parse_cavity_draft_source, &
                               CAVITY_DRAFT_NONE, CAVITY_DRAFT_FLAT, CAVITY_DRAFT_LINEAR, &
                               CAVITY_DRAFT_FILE, &
                               CAVITY_SOURCE_DRAFT, CAVITY_SOURCE_THICKNESS, &
                               set_draft_flat, set_draft_linear, &
                               parse_cavity_draft_sign, cavity_draft_apply_sign, &
                               CAVITY_SIGN_INVALID, &
                               cavity_water_column_impl, cavity_apply_land_exclusion, &
                               cavity_count_grounded, cavity_fill_cover_frac, &
                               cavity_draft_is_finite_nonneg, CAVITY_BOUND_INF, &
                               cavity_trim_eta_linear_impl
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
#ifdef RDB_ENABLE_TESTING
   public :: ocean_state_restart_write_drop_field
      !! PR-2 (bt-rem-from-av-rem review): exposed test-only so
      !! `test_ocean_restart_engine` can write a checkpoint with ONE
      !! named registered field held back under a different tag --
      !! simulating an OLDER checkpoint schema that predates that
      !! field (e.g. a pre-PR-2 file with no `vmix_kv`) without any
      !! NetCDF-level file surgery.  No production call site.
#endif
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
   public :: set_bathymetry_isomip_plus
   public :: isomip_plus_bx, isomip_plus_by
   public :: ISOMIP_B0, ISOMIP_B2, ISOMIP_B4, ISOMIP_B6
   public :: ISOMIP_XBAR, ISOMIP_DC, ISOMIP_FC, ISOMIP_WC, ISOMIP_ZB_DEEP
   public :: apply_layer_rho_init
   public :: ocean_linear_layer_density
   public :: pgf_nonoverlap_gate_on

   ! ==================================================================
   ! MISMIP+ / ISOMIP+ bedrock coefficients
   !
   ! Asay-Davis, Cornford, Durand, Galton-Fenzi, Gladstone, Gudmundsson,
   ! Hattermann, Holland, Holland, Holland, Martin, Mathiot, Pattyn and
   ! Seroussi (2016), "Experimental design for three interrelated marine
   ! ice sheet and ocean model intercomparison projects: MISMIP v. 3
   ! (MISMIP+), ISOMIP v. 2 (ISOMIP+) and MISOMIP v. 1 (MISOMIP1)",
   ! Geosci. Model Dev. 9, 2471-2497, doi:10.5194/gmd-9-2471-2016.
   !
   ! Table 1 (MISMIP+ parameters); the bed itself is their Eqs. (1)-(4).
   ! ISOMIP+ reuses it verbatim ("The bathymetry is the same as in
   ! Eq. (1)", their Sect. 3.1.1).  Every value below is METRES and is
   ! an ELEVATION coefficient (z positive UP, sea level at z = 0) —
   ! Roundabout's `b` is the opposite sign convention, which is why the
   ! setter negates at the very end and nowhere else.
   ! ==================================================================
   real(wp), parameter :: ISOMIP_B0 = -150.0_wp
      !! `B0`, bedrock elevation at x = 0 (m).  Table 1.
   real(wp), parameter :: ISOMIP_B2 = -728.8_wp
      !! `B2`, second bedrock coefficient (m).  Table 1.
   real(wp), parameter :: ISOMIP_B4 = 343.91_wp
      !! `B4`, third bedrock coefficient (m).  Table 1.
   real(wp), parameter :: ISOMIP_B6 = -50.57_wp
      !! `B6`, fourth bedrock coefficient (m).  Table 1.
   real(wp), parameter :: ISOMIP_XBAR = 300.0e3_wp
      !! `x_bar`, along-flow length scale of the bedrock (m).  Table 1.
   real(wp), parameter :: ISOMIP_DC = 500.0_wp
      !! `d_c`, depth of the trough relative to the side walls (m).  Table 1.
   real(wp), parameter :: ISOMIP_FC = 4.0e3_wp
      !! `f_c`, characteristic width of the channel side walls (m).  Table 1.
   real(wp), parameter :: ISOMIP_WC = 24.0e3_wp
      !! `w_c`, half-width of the trough (m).  Table 1.
   real(wp), parameter :: ISOMIP_ZB_DEEP = -720.0_wp
      !! `z_b,deep`, maximum depth of the bedrock (m, ELEVATION so
      !! negative).  Table 1.  The dispatch does NOT hard-code this: the
      !! clip comes from `&ocean_topo_nml max_depth`, so the protocol
      !! value is `max_depth = 720.0`.  Exported for the unit test and
      !! for anyone writing an ISOMIP+ namelist.

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
      type(atm_state_t) :: atm
         !! JRA55-do-style prescribed atmospheric state (OM3 C2c-2a) --
         !! `tas, huss, uas, vas, psl, rlds, rsds, prra, prsn` on the
         !! model T-grid. Unlike `data_forcing`, this slot OWNS real
         !! device-mapped arrays (the nine fields + regrid weights +
         !! window/scratch buffers), so it DOES carry terms in
         !! `ocean_state_enter_data`/`exit_data`/`bytes`/`destroy`
         !! below. Default off (`&ocean_atm_forcing_nml enable =
         !! .false.`) => every array stays unallocated => bit-identical.
         !! No consumer yet (the bulk-flux driver + ice seam are C2c-3).
      type(bulk_forcing_t) :: bulk_forcing
         !! Ocean consumer of `atm` (OM3 C2c-3a): persistent cell-centred
         !! stress scratch (`taux_c`/`tauy_c`) + the open-water fraction
         !! used to weight every component when sea ice is enabled.
         !! Default off (`&ocean_bulk_flux_nml enable = .false.`) =>
         !! unallocated => bit-identical. No ice-seam consumer yet
         !! (`ice_atm_forcing_from_atm` is C2c-3b).
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

   interface
      module subroutine ocean_state_seed_from_cfg(state, grid, cfg, ierr, injected_b, injected_b_convention, &
                                                  periodic_x, periodic_y)
            !! See submodule [[rdb_ocean_state_seed]].
         implicit none
         type(ocean_state_t), intent(inout) :: state
         type(hgrid_t), intent(in) :: grid
         type(config_t), intent(in) :: cfg
         integer, intent(out), optional :: ierr
         real(wp), intent(in), optional :: injected_b(:, :)
         integer, intent(in), optional :: injected_b_convention
         logical, intent(in), optional :: periodic_x
         logical, intent(in), optional :: periodic_y
      end subroutine ocean_state_seed_from_cfg

      module subroutine ocean_state_seed_land_cells(state, grid)
            !! See submodule [[rdb_ocean_state_seed]].
         implicit none
         type(ocean_state_t), intent(inout) :: state
         type(hgrid_t), intent(in) :: grid
      end subroutine ocean_state_seed_land_cells

      module subroutine apply_layer_rho_init(ms, layer_rho_init, nz_ml, ierr)
            !! See submodule [[rdb_ocean_state_seed]].
         implicit none
         type(multilayer_state_t), intent(inout) :: ms
         real(wp), intent(in) :: layer_rho_init(:)
         integer, intent(in) :: nz_ml
         integer, intent(out), optional :: ierr
      end subroutine apply_layer_rho_init

      module subroutine seed_eady_ic(state, grid, cfg, ierr)
            !! See submodule [[rdb_ocean_state_seed]].
         implicit none
         type(ocean_state_t), intent(inout) :: state
         type(hgrid_t), intent(in) :: grid
         type(config_t), intent(in) :: cfg
         integer, intent(out), optional :: ierr
      end subroutine seed_eady_ic

      module subroutine seed_baroclinic_jet_ic(state, grid, cfg, ierr)
            !! See submodule [[rdb_ocean_state_seed]].
         implicit none
         type(ocean_state_t), intent(inout) :: state
         type(hgrid_t), intent(in) :: grid
         type(config_t), intent(in) :: cfg
         integer, intent(out), optional :: ierr
      end subroutine seed_baroclinic_jet_ic

      module subroutine seed_geostrophic_adjustment_ic(state, grid, cfg, ierr)
            !! See submodule [[rdb_ocean_state_seed]].
         implicit none
         type(ocean_state_t), intent(inout) :: state
         type(hgrid_t), intent(in) :: grid
         type(config_t), intent(in) :: cfg
         integer, intent(out), optional :: ierr
      end subroutine seed_geostrophic_adjustment_ic

      module subroutine seed_cavity_draft(state, grid, cfg, ierr)
            !! See submodule [[rdb_ocean_state_seed]].
         implicit none
         type(ocean_state_t), intent(inout) :: state
         type(hgrid_t), intent(in) :: grid
         type(config_t), intent(in) :: cfg
         integer, intent(out) :: ierr
      end subroutine seed_cavity_draft

      module subroutine set_bathymetry_island(b, grid, max_depth, half_frac)
            !! See submodule [[rdb_ocean_state_topo]].
         implicit none
         real(wp), intent(inout) :: b(:, :)
         type(hgrid_t), intent(in) :: grid
         real(wp), intent(in) :: max_depth, half_frac
      end subroutine set_bathymetry_island

      module subroutine set_bathymetry_double_drake(b, grid, max_depth, half_frac)
            !! See submodule [[rdb_ocean_state_topo]].
         implicit none
         real(wp), intent(inout) :: b(:, :)
         type(hgrid_t), intent(in) :: grid
         real(wp), intent(in) :: max_depth, half_frac
      end subroutine set_bathymetry_double_drake

      module subroutine set_bathymetry_spoon(b, grid, max_depth, edge_depth, slope_scale)
            !! See submodule [[rdb_ocean_state_topo]].
         implicit none
         real(wp), intent(inout) :: b(:, :)
         type(hgrid_t), intent(in) :: grid
         real(wp), intent(in) :: max_depth, edge_depth, slope_scale
      end subroutine set_bathymetry_spoon

      module subroutine set_bathymetry_seamount(b, grid, max_depth, peak_depth, half_width)
            !! See submodule [[rdb_ocean_state_topo]].
         implicit none
         real(wp), intent(inout) :: b(:, :)
         type(hgrid_t), intent(in) :: grid
         real(wp), intent(in) :: max_depth, peak_depth, half_width
      end subroutine set_bathymetry_seamount

      module subroutine set_bathymetry_isomip_plus(b, grid, max_depth, x_origin, m_per_grid)
            !! See submodule [[rdb_ocean_state_topo]].
         implicit none
         real(wp), intent(inout) :: b(:, :)
         type(hgrid_t), intent(in) :: grid
         real(wp), intent(in) :: max_depth
         real(wp), intent(in) :: x_origin
         real(wp), intent(in) :: m_per_grid
      end subroutine set_bathymetry_isomip_plus

      module subroutine set_bathymetry_neverworld2(b, grid, max_depth, nl_continent_amp, nl_roughness_amp, min_depth)
            !! See submodule [[rdb_ocean_state_topo]].
         implicit none
         real(wp), intent(inout) :: b(:, :)
         type(hgrid_t), intent(in) :: grid
         real(wp), intent(in) :: max_depth, nl_continent_amp, nl_roughness_amp
         real(wp), intent(in) :: min_depth
      end subroutine set_bathymetry_neverworld2

      module subroutine ocean_state_build_restart_registry(state, grid, reg)
            !! See submodule [[rdb_ocean_state_restart]].
         implicit none
         type(ocean_state_t), intent(in), target :: state
         type(hgrid_t), intent(in) :: grid
         type(restart_registry_t), intent(inout) :: reg
      end subroutine ocean_state_build_restart_registry

#ifndef RDB_NO_NETCDF
      module subroutine ocean_state_fill_restart_metadata(state, grid, meta)
            !! See submodule [[rdb_ocean_state_restart]].
         implicit none
         type(ocean_state_t), intent(in) :: state
         type(hgrid_t), intent(in) :: grid
         type(ocean_restart_metadata_t), intent(out) :: meta
      end subroutine ocean_state_fill_restart_metadata
#endif

      module subroutine ocean_state_restart_read(state, grid, decomp, filename, t, step, ierr)
            !! See submodule [[rdb_ocean_state_restart]].
         implicit none
         type(ocean_state_t), intent(inout), target :: state
         type(hgrid_t), intent(in) :: grid
         type(decomp_t), intent(in) :: decomp
         character(len=*), intent(in) :: filename
         real(wp), intent(out) :: t
         integer, intent(out) :: step
         integer, intent(out), optional :: ierr
      end subroutine ocean_state_restart_read
   end interface

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
      nbytes = nbytes + this%data_input%bytes() + this%atm%bytes() + this%bulk_forcing%bytes()
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
      call profiler_start("ed_atm", nvtx_only=.true.)
      call state%atm%enter_data()
      call profiler_stop("ed_atm")
      call profiler_start("ed_bulk_forcing", nvtx_only=.true.)
      call state%bulk_forcing%enter_data()
      call profiler_stop("ed_bulk_forcing")
#endif
   end subroutine ocean_state_enter_data

   subroutine ocean_state_exit_data(state)
      !! Reverse of `ocean_state_enter_data`.  Walk components first,
      !! parent struct last.
      type(ocean_state_t), intent(inout) :: state

#ifndef RDB_NO_NETCDF
      call profiler_start("xd_bulk_forcing", nvtx_only=.true.)
      call state%bulk_forcing%exit_data()
      call profiler_stop("xd_bulk_forcing")
      call profiler_start("xd_atm", nvtx_only=.true.)
      call state%atm%exit_data()
      call profiler_stop("xd_atm")
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
      call this%bulk_forcing%destroy()
      call this%atm%destroy()
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

   subroutine ocean_state_restart_write_drop_field(state, grid, decomp, filename, t, step, &
                                                   drop_tag, ierr)
      !! Test-only sibling of `ocean_state_restart_write`: identical,
      !! except the ONE registered entry whose tag matches `drop_tag` is
      !! written under a mangled tag instead of its real one, so a
      !! subsequent `ocean_state_restart_read`/`ocean_restart_read_local`
      !! -- which looks up variables by their REAL tag -- finds nothing
      !! and takes exactly the "optional field absent" path it takes for
      !! a genuinely older checkpoint.  `error stop`s if `drop_tag` does
      !! not match any registered entry (a typo here must not silently
      !! test nothing).
      type(ocean_state_t), intent(inout), target :: state
      type(hgrid_t), intent(in) :: grid
      type(decomp_t), intent(in) :: decomp
      character(len=*), intent(in) :: filename
      real(wp), intent(in) :: t
      integer, intent(in) :: step
      character(len=*), intent(in) :: drop_tag
      integer, intent(out), optional :: ierr
#ifndef RDB_NO_NETCDF
      type(restart_registry_t) :: reg
      type(ocean_restart_metadata_t) :: meta
      integer :: e
      logical :: found

      if (present(ierr)) ierr = OCEAN_STATUS_OK

      call ocean_state_build_restart_registry(state, grid, reg)
      call ocean_state_fill_restart_metadata(state, grid, meta)

      found = .false.
      do e = 1, reg%n
         associate (en => reg%entries(e))
            ! Device sync (same as ocean_state_restart_write): no `cycle`
            ! here, unlike that routine, because the tag check below must
            ! still run for every entry including host-only ones.
            if (en%device_mapped) then
               if (en%rank == 2) then
                  !$acc update self(en%p2)
               else if (en%rank == 3) then
                  !$acc update self(en%p3)
               end if
            end if
            if (trim(en%tag) == trim(drop_tag)) then
               en%tag = trim(drop_tag)//"_DROPPED_FOR_TEST"
               found = .true.
            end if
         end associate
      end do
      if (.not. found) then
         error stop "ocean_state_restart_write_drop_field: drop_tag not registered"
      end if

      call ocean_restart_write_local(filename, reg, decomp, meta, t, step, &
                                     state%dyn%outer_step_count, ierr=ierr)
#else
      call fail("ocean_state_restart_write_drop_field: built without NetCDF "// &
                "(RDB_ENABLE_NETCDF=ON required)", ierr, OCEAN_STATUS_ERR_IO)
#endif
   end subroutine ocean_state_restart_write_drop_field

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

   pure function isomip_logistic(t) result(r)
      !! `1/(1 + exp(t))`, saturated instead of overflowing.
      !!
      !! The two logistic terms of Eq. (4) reach `|t| ~ 34` over the
      !! ISOMIP+ box including ghost rows, but a caller with a much wider
      !! `y_len` (or a tiny `f_c`) would drive `exp(t)` past the `real64`
      !! overflow at `t ~ 709`.  Saturating at +/-`T_SAT` is exact to the
      !! last bit of `r` on both sides (`1/(1+exp(500))` underflows to 0
      !! and `1/(1+exp(-500))` rounds to 1 anyway), so this costs nothing
      !! and removes an Inf that would propagate as a NaN.
      !!
      !! `t` is built from grid positions and the Table-1 constants, all
      !! finite by construction, so the CLAUDE.md "if/else clamps launder
      !! NaN" trap does not apply: there is no path that feeds this a NaN.
      real(wp), intent(in) :: t
      real(wp) :: r
      real(wp), parameter :: T_SAT = 500.0_wp
      if (t > T_SAT) then
         r = 0.0_wp
      else if (t < -T_SAT) then
         r = 1.0_wp
      else
         r = 1.0_wp/(1.0_wp + exp(t))
      end if
   end function isomip_logistic

   pure function isomip_plus_bx(x) result(bx)
      !! Along-flow bedrock elevation, Asay-Davis et al. (2016) Eq. (2):
      !!
      !!     Bx(x) = B0 + B2*xt**2 + B4*xt**4 + B6*xt**6,   xt = x/x_bar
      !!
      !! `x` is the ABSOLUTE MISMIP+ along-flow coordinate in METRES
      !! (0 at the ice divide), NOT the model's domain-relative x — the
      !! ISOMIP+ ocean box starts at `x = 320 km` (their Table 3 `x0`),
      !! which is what `&ocean_topo_nml x_origin` supplies.  Result is an
      !! ELEVATION (m, positive up), so it is negative everywhere in the
      !! ISOMIP+ box.
      real(wp), intent(in) :: x
      real(wp) :: bx
      real(wp) :: xt, xt2
      xt = x/ISOMIP_XBAR
      xt2 = xt*xt
      bx = ISOMIP_B0 + xt2*(ISOMIP_B2 + xt2*(ISOMIP_B4 + xt2*ISOMIP_B6))
   end function isomip_plus_bx

   pure function isomip_plus_by(y, y_len) result(by)
      !! Across-flow bedrock elevation, Asay-Davis et al. (2016) Eq. (4):
      !!
      !!     By(y) = d_c/(1 + exp(-2*(y - Ly/2 - w_c)/f_c))
      !!           + d_c/(1 + exp( 2*(y - Ly/2 + w_c)/f_c))
      !!
      !! A two-sided logistic trough: ~0 within `|y - Ly/2| < w_c` (the
      !! trough floor) rising to `d_c` on both side walls over the
      !! `f_c` transition.  Note the paper's own caveat (their Fig. 1b):
      !! `By` is an OFFSET relative to the trough, not a transect — the
      !! bed is `Bx + By`, and `Bx` is never zero.
      !!
      !! `y` and `y_len` are METRES.  `y_len` is the model's own domain
      !! width, which for the prescribed ISOMIP+ box (`0 <= y <= 80 km`)
      !! is the paper's `Ly` — taking it from the grid rather than
      !! hard-coding 80 km keeps the trough centred in whatever box the
      !! caller actually built.
      real(wp), intent(in) :: y, y_len
      real(wp) :: by
      real(wp) :: yc
      yc = y - 0.5_wp*y_len
      by = ISOMIP_DC*isomip_logistic(-2.0_wp*(yc - ISOMIP_WC)/ISOMIP_FC) &
           + ISOMIP_DC*isomip_logistic(2.0_wp*(yc + ISOMIP_WC)/ISOMIP_FC)
   end function isomip_plus_by

end module rdb_ocean_state
