!! Ocean slot configuration from a parsed namelist config.
module rdb_ocean_setup
   !! Applies the per-concern `&ocean_*` namelist knobs onto an
   !! `ocean_state_t`'s kernel slots (bottom drag, vmix, lateral
   !! viscosity, vertical coord, PGF, barotropic split).  Pure
   !! cfg → slot wiring + rank-0 logging — no I/O, no NetCDF — so it
   !! is shared by the production driver and the NetCDF-free benchmark.
   !!
   !! Call order after `ocean_state%init_from_config` +
   !! `ocean_state_seed_from_cfg` + `register_default_tracers`:
   !!   configure_ocean_drag / _vmix / _lateral / _pgf / _bt / _bt_split
#ifdef LFORTRAN_PASSING
   use rdb_constants, only: wp, GRAVITY, LAND_DEPTH_THRESHOLD
#else
   use rdb_constants, only: wp, GRAVITY, LAND_DEPTH_THRESHOLD, NZ_STACK_MAX
#endif
   use rdb_constants, only: H_VANISHED, VCOORD_ZSTAR_FULL
   use rdb_config, only: config_t
   use rdb_grid, only: hgrid_t
   use rdb_decomp, only: decomp_t
   use rdb_ocean_state, only: ocean_state_t, ocean_state_seed_land_cells, &
                              pgf_nonoverlap_gate_on
   use rdb_ocean_pressure_force, only: OPGF_VARIANT_GPRIME, OPGF_VARIANT_FV_MOM6, &
                                       OPGF_VARIANT_FV_WRIGHT, parse_opgf_variant
   use rdb_eos, only: parse_eos_variant, eos_validate, &
                      EOS_VARIANT_ROQUET_SPV, &
                      parse_tfreeze_set, eos_apply_tfreeze_set, &
                      TFREEZE_SET_SEAICE
   use rdb_ocean_vcoord, only: parse_ocean_vcoord_type, VCOORD_EULERIAN_Z, &
                               VCOORD_RHO, VCOORD_HYCOM, VCOORD_LAGRANGIAN, &
                               VCOORD_Z_FIXED, ocean_vcoord_z_fixed_target, &
                               ocean_vcoord_closed_face_masks, &
                               ocean_vcoord_k_top_from_target, &
                               ocean_vcoord_set_z_fixed_profile, &
                               ocean_vcoord_count_ledges
   use rdb_vcoord, only: parse_remap_method, parse_z_fixed_profile, z_fixed_nominal_dz, &
                         ZFIXED_PROFILE_UNIFORM, ZFIXED_PROFILE_INVALID, ZFIXED_DZ_OK
   use rdb_ocean_bottom_drag, only: parse_bdrag_variant
   use rdb_ocean_top_drag, only: parse_tdrag_variant, TDRAG_QUADRATIC, &
                                 top_drag_fill_face_cover_impl
   use rdb_ocean_lateral_mix, only: parse_lateral_closure, LMIX_NONE, &
                                    LMIX_LEITH, LMIX_SMAGORINSKY, &
                                    LMIX_BIHARMONIC, LMIX_LEITH_BIHARM
   use rdb_ocean_horizontal_viscosity, only: ocean_hvisc_set_aniso_direction
   use rdb_ocean_surface_stress, only: ocean_surface_stress_set_derived, &
                                       ocean_surface_stress_apply_cover
   use rdb_ocean_surface_flux, only: ocean_surface_flux_apply_cover_const
   use rdb_coriolis_adv, only: parse_pv_variant, PV_VARIANT_SADOURNY_HK, &
                               CORNER_H_CELL_MEAN, CORNER_H_MOM6_AREA, &
                               parse_pv_adv_scheme
   use rdb_ocean_boundary_types, only: ocean_bc_type_from_string, &
                                       ocean_bc_validate_periodic, &
                                       ocean_bc_validate_fold, &
                                       OBC_WALL, OBC_PERIODIC, OBC_OPEN, OBC_TIDAL, &
                                       OBC_CHAPMAN, OBC_NESTED, OBC_CLAMPED, OBC_SPONGE, &
                                       OBC_TRIPOLAR_FOLD, &
                                       OBC_MAX_TIDAL_CONSTITUENTS, &
                                       ocean_bc_face_tag_t, &
                                       obc_tide_nodal_fill, obc_match_constituent
   use rdb_ocean_sponge, only: sponge_band_alpha, SPONGE_RAMP_COSINE, SPONGE_RAMP_LINEAR
   use rdb_ocean_halo, only: ocean_halo_centre, ocean_halo_is_decomposed
   use rdb_ocean_halo_state, only: ocean_seam_refresh_surface_stress
   use rdb_ocean_metrics, only: metrics_finalize, metrics_fill_cartesian, &
                                metrics_fill_spherical, metrics_fill_from_supergrid, &
                                metrics_fill_tripolar, metrics_apply_land_mask, &
                                metrics_fill_coriolis, parse_grid_config, &
                                parse_coriolis_scheme, ocean_metrics_t, &
                                GRID_CONFIG_CARTESIAN, GRID_CONFIG_SPHERICAL, &
                                GRID_CONFIG_SUPERGRID, GRID_CONFIG_TRIPOLAR, &
                                metrics_porous_alloc, metrics_closed_faces_alloc
   use rdb_ocean_cavity, only: cavity_fill_p_ice_ref, &
                               cavity_datum_impl, cavity_datum_residual
   use rdb_ocean_cavity_melt, only: CAVITY_GAMMA_RATIO_ISOMIP
   use rdb_ocean_porous, only: parse_porous_source, parse_porous_eta_interp, &
                               porous_fill_stats_resolved, &
                               porous_stats_are_ordered, &
                               POROUS_SOURCE_RESOLVED, POROUS_SOURCE_FILE
   use rdb_ocean_epbl, only: parse_epbl_mstar_scheme, parse_epbl_vstar_scheme, &
                             parse_epbl_combine, parse_epbl_lt_scheme
   use rdb_ocean_vmix, only: parse_kpp_sw_method, vmix_resolve_kd_min, &
                             bkgnd_henyey_conflicts_profile, &
                             parse_buoyancy_coeffs, BUOY_COEFFS_INVALID
   use rdb_ocean_geothermal, only: ocean_geothermal_t
   use rdb_ocean_fold, only: fold_north_corner
   use rdb_ocean_tides, only: tides_configure_astronomy, tides_build_struct
   use rdb_ocean_p_surf, only: p_surf_configure
   use rdb_ocean_tide_astro, only: tide_name_index, days_since_1900, &
                                   parse_date_string, TIDES_CATALOG_SIZE, &
                                   equilibrium_arguments, nodal_fu, TIDE_NAME
   use rdb_ocean_dyn, only: ocean_dt_tracer_advect_ratios_ok, &
                            SPLIT_SCHEME_SSP_RK2, SPLIT_SCHEME_PRED_CORR
   use rdb_recon_weno, only: parse_tracer_recon, TRACER_RECON_PPM
   use rdb_halo, only: halo_allreduce_min
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_SETUP
   use rdb_error_ring, only: fail
   use pic_logger, only: logger => global_logger
   use pic_strings, only: to_string
   implicit none
   private

#ifdef LFORTRAN_PASSING
   integer, parameter :: NZ_STACK_MAX = 64
      !! LFortran 0.64 workaround: module-local copy of the rdb_constants value
      !! (an imported parameter used as an explicit-shape dummy bound inside a
      !! PURE call becomes an impure getter under LFortran). Keep in sync (=64).
#endif

   public :: configure_ocean_metrics
   public :: configure_ocean_land_mask
   public :: configure_ocean_forcing
   public :: configure_ocean_drag
   public :: configure_ocean_hdiff
   public :: configure_ocean_vmix
   public :: configure_ocean_tracers
   public :: configure_ocean_lateral
   public :: configure_ocean_reference_density
   public :: configure_ocean_pgf
   public :: configure_ocean_bt
   public :: configure_ocean_bt_split
   public :: configure_ocean_tides
   public :: configure_ocean_p_surf
   public :: configure_ocean_wave_drag
   public :: wave_drag_roughness_proxy
   public :: configure_ocean_porous
   public :: configure_ocean_closed_faces
   public :: configure_ocean_k_top
   public :: configure_ocean_z_fixed_profile
   public :: configure_ocean_cavity
   public :: configure_ocean_cavity_melt
   public :: configure_ocean_top_drag
   public :: cavity_resolve_gamma_s
   public :: cavity_count_zero_f
   public :: cavity_count_unloaded_p_top
   public :: configure_ocean_wetdry
   public :: bt_auto_n_inner
   public :: metrics_bt_cfl_length
   public :: bt_auto_n_inner_from_dt
   public :: bt_cfl_dt_wet
   public :: configure_ocean_bc
   public :: configure_ocean_sponge

contains

   subroutine configure_ocean_metrics(cfg, ocean_state, grid, compute_rank, ierr)
      !! Fill the `ocean_metrics_t` slot per `cfg%ocean%grid%grid_config`,
      !! then single-source the inverses + hvisc ratio bundle
      !! (`metrics_finalize`).  Must run BEFORE `ocean_state_enter_data`
      !! (the host fill is what the GPU copyin captures).  For
      !! "spherical", the `&grid_nml dx`/`dy` are reinterpreted as
      !! dlon/dlat in degrees.
      type(config_t), intent(in) :: cfg
      type(ocean_state_t), intent(inout) :: ocean_state
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: compute_rank
      integer, intent(out), optional :: ierr
         !! Non-zero on a grid/geometry configuration failure when present;
         !! absent behaves as today (`error stop`).

      integer :: gc, local_ierr

      gc = parse_grid_config(cfg%ocean%grid%grid_config)
      select case (gc)
      case (GRID_CONFIG_SPHERICAL)
         ! Validation the schema can't express (cross-knob domain check).
         if (cfg%ocean%grid%rad_earth <= 0.0_wp) then
            call fail("configure_ocean_metrics: spherical grid needs rad_earth > 0", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
         ! Include ghost rows: the generator fills ghost rows beyond the physical
         ! edge by formula (see metrics_fill_spherical), so the most extreme
         ! latitude reached is lat_south - (nghost+0.5)*dlat (south ghost centre)
         ! and lat_south + (ny_phys+nghost+0.5)*dlat (north ghost centre).
         ! Near-pole sectors need tripolar treatment (M3 follow-up).
         if (cfg%ocean%grid%lat_south - (real(grid%nghost, wp) + 0.5_wp)*grid%dy < -90.0_wp .or. &
             cfg%ocean%grid%lat_south + real(grid%ny_phys + grid%nghost, wp)*grid%dy + &
             0.5_wp*grid%dy > 90.0_wp) then
            call fail("configure_ocean_metrics: spherical domain (including ghost rows) "// &
                      "exceeds |lat| <= 90. Near-pole sectors need tripolar treatment (M3).", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
         ! grid%dx / grid%dy are dlon / dlat in DEGREES here.
         call metrics_fill_spherical(ocean_state%metrics, grid, &
                                     cfg%ocean%grid%lon_west, cfg%ocean%grid%lat_south, &
                                     grid%dx, grid%dy, cfg%ocean%grid%rad_earth)
         if (compute_rank == 0) then
            call logger%info("Grid config:      spherical lon-lat sector")
         end if
      case (GRID_CONFIG_SUPERGRID)
         if (len_trim(cfg%ocean%grid%supergrid_file) == 0) then
            call fail("configure_ocean_metrics: grid_config='supergrid' "// &
                      "requires supergrid_file to be set in &ocean_grid_nml", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
         ! ierr threaded down ONLY when THIS routine's own ierr is
         ! present: otherwise metrics_fill_from_supergrid must keep
         ! reaching its own `error stop` (specific dimension-mismatch /
         ! NetCDF-failure text) rather than the generic wrapper message
         ! below (P0.1 review F2).
         !
         ! The edge tags decide the ghost-metric topology (periodic-x wrap,
         ! tripolar fold) exactly as they do for the analytic tripolar; the
         ! reader cross-checks the fold against the file's own top row.
         if (present(ierr)) then
            call metrics_fill_from_supergrid(ocean_state%metrics, grid, &
                                             cfg%ocean%grid%supergrid_file, ierr=local_ierr, &
                                             periodic_x=tags_periodic_x(cfg), &
                                             north_fold=tags_north_fold(cfg))
            if (local_ierr /= 0) then
               ierr = local_ierr
               return
            end if
         else
            call metrics_fill_from_supergrid(ocean_state%metrics, grid, &
                                             cfg%ocean%grid%supergrid_file, &
                                             periodic_x=tags_periodic_x(cfg), &
                                             north_fold=tags_north_fold(cfg))
         end if
         if (compute_rank == 0) then
            call logger%info("Grid config:      supergrid (mosaic) from "// &
                             trim(cfg%ocean%grid%supergrid_file))
         end if
      case (GRID_CONFIG_TRIPOLAR)
         ! Cross-knob domain checks the schema can't express.
         if (cfg%ocean%grid%rad_earth <= 0.0_wp) then
            call fail("configure_ocean_metrics: tripolar grid needs rad_earth > 0", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
         if (cfg%ocean%grid%phi_join <= cfg%ocean%grid%lat_south) then
            call fail("configure_ocean_metrics: tripolar phi_join must be "// &
                      "north of lat_south (the cap sits above the lon-lat region)", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
         ! Tripolar requires an east-west PERIODIC grid: the i-direction
         ! wraps the full 360 deg of pseudo-longitude and the north fold
         ! identifies the two halves of the top row.  Without periodic-x
         ! the wrap + fold exchange (M4c) is ill-defined.
         if (trim(cfg%ocean%bc%west) /= "periodic" .or. &
             trim(cfg%ocean%bc%east) /= "periodic") then
            call fail("configure_ocean_metrics: grid_config='tripolar' requires "// &
                      "periodic east-west BCs (&ocean_bc_nml west='periodic' east='periodic')", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
         ! A tripolar grid MUST close the top with the north fold; otherwise
         ! the duplicated top-row DOFs are never reconciled and the cap
         ! diverges.  (The reverse — fold tag requires a tripolar grid — is
         ! enforced by ocean_bc_validate_fold's periodic-w/e requirement plus
         ! this check.)
         if (trim(cfg%ocean%bc%north) /= "tripolar_fold") then
            call fail("configure_ocean_metrics: grid_config='tripolar' requires "// &
                      "north='tripolar_fold' (&ocean_bc_nml north='tripolar_fold')", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
         ! grid%dx / grid%dy are dlon / dlat in DEGREES here (as spherical).
         call metrics_fill_tripolar(ocean_state%metrics, grid, &
                                    cfg%ocean%grid%lon_west, cfg%ocean%grid%lat_south, &
                                    grid%dx, grid%dy, cfg%ocean%grid%rad_earth, &
                                    cfg%ocean%grid%phi_join, cfg%ocean%grid%lon_pole)
         if (compute_rank == 0) then
            call logger%info("Grid config:      tripolar (Murray 1996), phi_join="// &
                             to_string(cfg%ocean%grid%phi_join)//" lon_pole="// &
                             to_string(cfg%ocean%grid%lon_pole))
         end if
      case default   ! GRID_CONFIG_CARTESIAN
         call metrics_fill_cartesian(ocean_state%metrics, grid, grid%dx, grid%dy)
      end select

      call metrics_finalize(ocean_state%metrics)
      if (present(ierr)) ierr = OCEAN_STATUS_OK
   end subroutine configure_ocean_metrics

   pure function tags_periodic_x(cfg) result(per_x)
      !! `.true.` iff `&ocean_bc_nml` tags BOTH west and east `periodic` —
      !! the rule `ocean_bc_state_init` applies, for callers that run before
      !! `configure_ocean_bc` (the grid metrics).
      type(config_t), intent(in) :: cfg
      logical :: per_x
      per_x = ocean_bc_type_from_string(cfg%ocean%bc%west) == OBC_PERIODIC .and. &
              ocean_bc_type_from_string(cfg%ocean%bc%east) == OBC_PERIODIC
   end function tags_periodic_x

   pure function tags_north_fold(cfg) result(fold)
      !! `.true.` iff `&ocean_bc_nml north = "tripolar_fold"`.
      type(config_t), intent(in) :: cfg
      logical :: fold
      fold = ocean_bc_type_from_string(cfg%ocean%bc%north) == OBC_TRIPOLAR_FOLD
   end function tags_north_fold

   subroutine configure_ocean_land_mask(cfg, ocean_state, grid, compute_rank)
      !! Derive the static C-grid land masks from the seeded T-cell
      !! `wet_mask` and zero the 6 face metrics at land faces
      !! (`metrics_apply_land_mask`).  Run AFTER `configure_ocean_metrics`
      !! (the metrics + inverses must exist) AND `configure_ocean_bc` (the
      !! periodic / fold flags drive the halo-aware mask derivation), but
      !! BEFORE `ocean_state_enter_data` (the host edit is what the GPU
      !! copyin captures).
      !!
      !! For a domain with no land (`wet_mask≡1`, the flat-bottom /
      !! analytical path) every mask is 1.0 and the metric multiply is a
      !! literal no-op ⇒ bit-identical to a no-mask build.
      type(config_t), intent(in) :: cfg
      type(ocean_state_t), intent(inout) :: ocean_state
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: compute_rank

      ! Multi-rank seam ghost fill of wet_mask (O3 land x decomp): a subdomain
      ! seam that bisects a continent needs the NEIGHBOUR rank's wet_T in the
      ! seam ghost columns so the C-grid face/corner masks derived inside
      ! metrics_apply_land_mask close the seam faces correctly.  Done HERE (the
      ! caller) rather than inside metrics_apply_land_mask because importing the
      ! comm-layer rdb_ocean_halo into the low-level rdb_ocean_metrics leaf
      ! creates an NVFORTRAN USE cycle (config -> lateral_mix -> metrics).
      ! ocean_halo_init has already run (driver orders it before this call).
      ! Gated on ocean_halo_is_DECOMPOSED (px>1 or py>1) — a genuine
      ! cross-rank seam is the ONLY case wet_mask needs a message: the
      ! periodic ghost wrap is done locally inside metrics_apply_land_mask,
      ! and a single-rank non-periodic run's ghosts are domain-boundary
      ! walls.  Gating on the seam (not periodicity / not merely is_init)
      ! also keeps the MPI runtime untouched in non-MPI ocean unit tests
      ! that init the halo via this setup path but never call MPI_Init
      ! (those are px=py=1 ⇒ skip).  Runs host-side before enter_data; a
      ! GPU+MPI init-time host-halo variant is an O3+ follow-up.
      if (ocean_halo_is_decomposed()) then
         call ocean_halo_centre(ocean_state%multilayer%wet_mask, device_resident=.false.)
      end if

      ! Solid-wall velocity masking (&ocean_bc_nml mask_wall_velocity):
      ! thread the per-edge WALL flags so only genuine solid walls get their
      ! ghost wet_T zeroed (periodic + open/OBC edges keep their ghosts).
      !
      ! The `has_*` conjunct is load-bearing under MPI.  `bc_type` is the
      ! GLOBAL edge tag, identical on every rank; `has_*` is the per-rank
      ! "I actually sit on that domain edge" flag (set from decomp).
      ! Without the AND, every rank zeroed its own west/east/south/north
      ! ghost band, so an interior MPI SEAM was masked as land: wet_u at
      ! the seam face went to 0 and the six masked face metrics with it.
      ! Measured on double_gyre_mom6 (44x40, 10 days, wall basin): px=1
      ! gave En 3.036E-03 / MaxCFL 0.01000, px=2 gave 2.998E-03 / 0.01197
      ! — and MaxCFL is a max reduction, which is exactly order-invariant,
      ! so that spread could not have been reduction roundoff.  With the
      ! conjunct, px=1 and px=2 agree.
      call metrics_apply_land_mask(ocean_state%metrics, &
                                   ocean_state%multilayer%wet_mask, grid, &
                                   ocean_state%bc%periodic_x, &
                                   ocean_state%bc%periodic_y, &
                                   ocean_state%bc%north_fold, &
                                   mask_wall_velocity=cfg%ocean%bc%mask_wall_velocity, &
                                   wall_west=(ocean_state%bc%west%bc_type == OBC_WALL &
                                              .and. ocean_state%bc%has_west), &
                                   wall_east=(ocean_state%bc%east%bc_type == OBC_WALL &
                                              .and. ocean_state%bc%has_east), &
                                   wall_south=(ocean_state%bc%south%bc_type == OBC_WALL &
                                               .and. ocean_state%bc%has_south), &
                                   wall_north=(ocean_state%bc%north%bc_type == OBC_WALL &
                                               .and. ocean_state%bc%has_north))

      ! Hold land T-cells at finite reference values so masked arithmetic
      ! never multiplies 0 by a NaN (0*NaN = NaN).  Land h_layer is floored
      ! to H_VANISHED (never 0), tracers held at the IC reference, layer +
      ! face velocities zeroed on land.
      call ocean_state_seed_land_cells(ocean_state, grid)

      if (compute_rank == 0) then
         associate (unused => cfg%ocean%grid%grid_config)
         end associate
      end if
   end subroutine configure_ocean_land_mask

   pure function metrics_dx_min(metrics, grid) result(dx_min)
      !! Representative minimum grid length over the PHYSICAL region,
      !! taken over both `dxT` and `dyT` (host-side, configure time).
      !! Bit-identical to `min(grid%dx, grid%dy)` on uniform Cartesian.
      type(ocean_metrics_t), intent(in) :: metrics
      type(hgrid_t), intent(in) :: grid
      real(wp) :: dx_min
      integer :: ng, i0, i1, j0, j1
      ng = grid%nghost
      i0 = ng + 1
      i1 = ng + grid%nx_phys
      j0 = ng + 1
      j1 = ng + grid%ny_phys
      dx_min = min(minval(metrics%dxT(i0:i1, j0:j1)), &
                   minval(metrics%dyT(i0:i1, j0:j1)))
   end function metrics_dx_min

   pure function metrics_bt_cfl_length(metrics, grid) result(l_cfl)
      !! 2-D external-gravity-wave CFL length over the PHYSICAL region:
      !!     l_cfl = min_cell  1 / sqrt(1/dxT^2 + 1/dyT^2).
      !! Length scale for the barotropic CFL `c_ext*dt*sqrt(1/dx^2+1/dy^2) <= 1`
      !! (includes the cross-direction term; on uniform Cartesian = dx/sqrt(2)).
      !! Handles anisotropic cells exactly.  Host-side, configure time.
      type(ocean_metrics_t), intent(in) :: metrics
      type(hgrid_t), intent(in) :: grid
      real(wp) :: l_cfl, inv_l2_max
      integer :: ng, i0, i1, j0, j1
      ng = grid%nghost
      i0 = ng + 1
      i1 = ng + grid%nx_phys
      j0 = ng + 1
      j1 = ng + grid%ny_phys
      inv_l2_max = maxval(1.0_wp/metrics%dxT(i0:i1, j0:j1)**2 &
                          + 1.0_wp/metrics%dyT(i0:i1, j0:j1)**2)
      l_cfl = 1.0_wp/sqrt(inv_l2_max)
   end function metrics_bt_cfl_length

   pure function bt_auto_n_inner(dt_outer, cfl_safety, c_ext, l_cfl) result(n_inner)
      !! Smallest barotropic substep count `n_inner` such that the substep
      !! `dt_outer/n_inner` satisfies the 2-D external-gravity-wave CFL:
      !!     dt_bt <= cfl_safety * l_cfl / c_ext,
      !! with `l_cfl` the 2-D CFL length (`metrics_bt_cfl_length`) and `c_ext`
      !! the external wave speed `sqrt(g*H_max)`.  MOM6 set_dtbt analogue, now
      !! with the cross-direction term included (the legacy estimate used a
      !! 1-D length and under-counted n_inner by ~sqrt(2) on square cells,
      !! leaving the effective 2-D CFL at ~0.92 for cfl_bt_safety=0.65 — on the
      !! edge of the forward-backward scheme's stability).
      real(wp), intent(in) :: dt_outer, cfl_safety, c_ext, l_cfl
      integer :: n_inner
      real(wp) :: dt_bt_safe
      dt_bt_safe = cfl_safety*l_cfl/c_ext
      n_inner = bt_auto_n_inner_from_dt(dt_outer, dt_bt_safe)
   end function bt_auto_n_inner

   pure function bt_auto_n_inner_from_dt(dt_outer, dt_bt_safe) result(n_inner)
      !! Smallest `n_inner >= 1` with `dt_outer/n_inner <= dt_bt_safe`,
      !! for an ALREADY-LIMITED safe barotropic substep `dt_bt_safe` (s) —
      !! the per-wet-cell minimum `bt_cfl_dt_wet` returns, reduced across
      !! ranks.  `bt_auto_n_inner` is this with `dt_bt_safe` formed from a
      !! single `(c_ext, l_cfl)` pair.
      real(wp), intent(in) :: dt_outer
         !! Outer (baroclinic) step (s).
      real(wp), intent(in) :: dt_bt_safe
         !! Largest stable barotropic substep, safety factor included (s).
      integer :: n_inner
      n_inner = max(1, ceiling(dt_outer/dt_bt_safe))
   end function bt_auto_n_inner_from_dt

   pure subroutine bt_cfl_dt_wet(nx, ny, i0, i1, j0, j1, b, wet, dxT, dyT, &
                                 cfl_safety, dt_bt, h_at, l_at, n_wet)
      !! Per-WET-CELL external-gravity-wave CFL limit (MOM6 `set_dtbt`):
      !!
      !!     dt_bt = min over wet (i,j) of  cfl_safety * l(i,j) / c(i,j),
      !!     l(i,j) = 1/sqrt(1/dxT(i,j)^2 + 1/dyT(i,j)^2),
      !!     c(i,j) = sqrt(g * max(b(i,j), 1 m)),
      !!
      !! i.e. the LOCAL depth with the LOCAL cell size, over OCEAN only.
      !! MOM6 evaluates the same quantity as `gtot*dt^2*(1/dx^2+1/dy^2)`
      !! per wet point.
      !!
      !! **Why per point.** The former estimate combined the deepest depth
      !! ANYWHERE with the smallest cell ANYWHERE — on the 1° tripolar grid
      !! that was 6000 m of ocean against a 362 m LAND cell at a
      !! land-locked bipole, and bought 1930 substeps where ~32 suffice.
      !! A land cell carries no gravity wave, and a small cell over a
      !! shallow shelf does not see the abyssal wave speed.
      !!
      !! **Bit-identity where the two agree.** Each point's value is
      !! evaluated with EXACTLY the arithmetic the global-extremes estimate
      !! used (`1/sqrt(inv_l2)`, `sqrt(g*max(b,1))`, `safety*l/c`, in that
      !! order).  Every step is monotone under round-to-nearest, so where
      !! the deepest wet column and the smallest wet cell COINCIDE (a
      !! flat-bottomed uniform grid, a wall basin with no land) the minimum
      !! is attained at that point and equals the old number bit-for-bit.
      !!
      !! No wet cell in the window ⇒ `dt_bt = huge`, `n_wet = 0` — the
      !! identity of the cross-rank `min` reduction (a rank that is all
      !! land does not constrain the others).
      !!
      !! Host-side, configure time: plain loops, no `do concurrent`.
      integer, intent(in) :: nx
         !! First extent of the centre arrays (ghosts included).
      integer, intent(in) :: ny
         !! Second extent of the centre arrays (ghosts included).
      integer, intent(in) :: i0
         !! First physical i.
      integer, intent(in) :: i1
         !! Last physical i.
      integer, intent(in) :: j0
         !! First physical j.
      integer, intent(in) :: j1
         !! Last physical j.
      real(wp), intent(in) :: b(nx, ny)
         !! Bed depth (m, positive down).
      real(wp), intent(in) :: wet(nx, ny)
         !! Static wet (1) / land (0) T-cell mask.
      real(wp), intent(in) :: dxT(nx, ny)
         !! T-cell x length (m).
      real(wp), intent(in) :: dyT(nx, ny)
         !! T-cell y length (m).
      real(wp), intent(in) :: cfl_safety
         !! `&ocean_bt_nml cfl_bt_safety`.
      real(wp), intent(out) :: dt_bt
         !! Smallest per-wet-cell safe substep (s); `huge` if no wet cell.
      real(wp), intent(out) :: h_at
         !! Bed depth at the limiting cell (m); 0 if no wet cell.
      real(wp), intent(out) :: l_at
         !! 2-D CFL length at the limiting cell (m); 0 if no wet cell.
      integer, intent(out) :: n_wet
         !! Wet cells scanned.
      integer :: i, j
      real(wp) :: inv_l2, l_ij, c_ij, dt_ij

      dt_bt = huge(1.0_wp)
      h_at = 0.0_wp
      l_at = 0.0_wp
      n_wet = 0
      do j = j0, j1
         do i = i0, i1
            if (wet(i, j) <= 0.5_wp) cycle
            n_wet = n_wet + 1
            inv_l2 = 1.0_wp/dxT(i, j)**2 + 1.0_wp/dyT(i, j)**2
            l_ij = 1.0_wp/sqrt(inv_l2)
            c_ij = sqrt(GRAVITY*max(b(i, j), 1.0_wp))
            dt_ij = cfl_safety*l_ij/c_ij
            if (dt_ij < dt_bt) then
               dt_bt = dt_ij
               h_at = b(i, j)
               l_at = l_ij
            end if
         end do
      end do
   end subroutine bt_cfl_dt_wet

   subroutine fill_coriolis_corner(cfg, metrics, grid, f_corner)
      !! Fill a C-grid corner Coriolis array via the single
      !! generator-driven routine `metrics_fill_coriolis`, honouring
      !! `&ocean_grid_nml coriolis_scheme` (D7).  beta_plane (default)
      !! is BIT-IDENTICAL to the legacy `coriolis_adv_set_beta_plane`;
      !! planetary uses the metrics geography.  The centre output is
      !! discarded here (filled into local scratch).
      type(config_t), intent(in) :: cfg
      type(ocean_metrics_t), intent(in) :: metrics
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(out) :: f_corner(:, :)
      real(wp), allocatable :: f_centre_scratch(:, :)
      allocate (f_centre_scratch(grid%nx_total, grid%ny_total))
      call metrics_fill_coriolis(metrics, &
                                 parse_coriolis_scheme(cfg%ocean%grid%coriolis_scheme), &
                                 cfg%coriolis_f, cfg%ocean%topo%coriolis_beta, &
                                 cfg%ocean%topo%coriolis_y_ref, cfg%ocean%grid%omega, &
                                 grid, f_corner, f_centre_scratch)
      deallocate (f_centre_scratch)
   end subroutine fill_coriolis_corner

   subroutine fill_coriolis_centre(cfg, metrics, grid, f_centre)
      !! Fill a cell-centre |Coriolis| array via `metrics_fill_coriolis`
      !! (D7).  beta_plane (default) is BIT-IDENTICAL to the legacy
      !! EPBL / kappa-shear `set_f_centre`.  The corner output is
      !! discarded here (filled into local scratch).
      type(config_t), intent(in) :: cfg
      type(ocean_metrics_t), intent(in) :: metrics
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(out) :: f_centre(:, :)
      real(wp), allocatable :: f_corner_scratch(:, :)
      allocate (f_corner_scratch(grid%nx_total + 1, grid%ny_total + 1))
      call metrics_fill_coriolis(metrics, &
                                 parse_coriolis_scheme(cfg%ocean%grid%coriolis_scheme), &
                                 cfg%coriolis_f, cfg%ocean%topo%coriolis_beta, &
                                 cfg%ocean%topo%coriolis_y_ref, cfg%ocean%grid%omega, &
                                 grid, f_corner_scratch, f_centre)
      deallocate (f_corner_scratch)
   end subroutine fill_coriolis_centre

   subroutine configure_ocean_forcing(cfg, ocean_state, grid, compute_rank, decomp)
      !! Wire surface wind stress, horizontal-viscosity coefficients, and the
      !! Coriolis beta-plane + PV-scheme variant from cfg into the ocean slots.
      !!
      !! When `decomp` is present the wind-stress setters receive the rank's
      !! global meridional offset and global domain extent so each subdomain
      !! seeds the correct portion of the global wind profile.  Absent ⇒
      !! single-rank byte-identical path.
      type(config_t), intent(in) :: cfg
      type(ocean_state_t), intent(inout) :: ocean_state
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: compute_rank
      type(decomp_t), intent(in), optional :: decomp

      integer :: joff, nyg

      joff = 0
      nyg = grid%ny_phys
      if (present(decomp)) then
         joff = decomp%j_start - 1
         nyg = decomp%ny_global
      end if

      ! Wind stress: dispatch on cfg%ocean%topo%wind_config.
      select case (trim(cfg%ocean%topo%wind_config))
      case ("2gyre")
         call ocean_state%surface_stress%set_wind_stress_2gyre( &
            grid, cfg%ocean%topo%taux_magnitude, &
            j_offset=joff, ny_global=nyg)
      case ("neverworld2")
         call ocean_state%surface_stress%set_wind_stress_neverworld2( &
            grid, cfg%ocean%topo%taux_magnitude, &
            j_offset=joff, ny_global=nyg)
      case default  ! "constant"
         ! Warn loudly when a wind was configured on the knob this path does
         ! NOT read. `taux_magnitude` is only consumed by the "2gyre" and
         ! "neverworld2" profiles above; the constant path reads
         ! `&physics_nml wind_stress_x/y`. Setting the former under
         ! `wind_config="constant"` is silently inert — which is exactly how
         ! `coriolis_coast.nml` shipped a case whose header promises a
         ! "wind-spun-up gyre" that it has never generated: it ran its full
         ! 10 days at En = 0.000E+00.
         if (abs(cfg%ocean%topo%taux_magnitude) > 0.0_wp .and. &
             abs(cfg%wind_stress_x) <= 0.0_wp .and. &
             abs(cfg%wind_stress_y) <= 0.0_wp) then
            call logger%warning( &
               "&ocean_topo_nml taux_magnitude = "// &
               to_string(cfg%ocean%topo%taux_magnitude)//" is INERT under "// &
               'wind_config="constant" (it is only read by "2gyre" / '// &
               '"neverworld2"), and &physics_nml wind_stress_x/y are both '// &
               "zero, so this run has NO wind forcing at all. Set "// &
               "&physics_nml wind_stress_x to apply a constant wind.")
         end if
         call ocean_state%surface_stress%set_wind_stress_const( &
            cfg%wind_stress_x, cfg%wind_stress_y)
      end select

      ! Make the freshly-seeded wind valid in the ghost bands, then derive
      ! stress_mag from it.  Host-side (device_resident=.false.): this runs
      ! in the driver's configure phase, ahead of ocean_state_enter_data.
      !
      ! The exchange matters even though the wind is configure-static.
      ! `set_wind_stress_2gyre` / `_neverworld2` fill PHYSICAL rows only and
      ! leave ghost rows at zero ("so wall faces see no spurious stress"),
      ! which is right at a true wall and wrong at an MPI seam — there the
      ! ghost belongs to the neighbour rank and must carry its wind, not 0.
      ! Nothing else ever repairs it, and both `stress_mag` (KPP/EPBL u_*)
      ! and the MLE corner average read one cell beyond their own.  At a
      ! true wall the exchange is a no-op, so the zero-ghost intent of the
      ! analytic setters is preserved exactly (and `test_wind_2gyre_ghosts`
      ! still holds).
      call ocean_seam_refresh_surface_stress(ocean_state%surface_stress, grid, &
                                             ocean_state%bc, device_resident=.false.)

      ! Horizontal momentum viscosity (default 0 keeps analytical tests bit-
      ! identical; realistic wind-driven runs need a non-zero value).
      ocean_state%hvisc%nu_h = cfg%ocean%hvisc%nu_h
      ocean_state%hvisc%nu_4 = cfg%ocean%hvisc%nu_4
      ocean_state%hvisc%stress_tensor = cfg%ocean%hvisc%stress_tensor
      ocean_state%hvisc%bound_coef = cfg%ocean%hvisc%bound_coef
      ocean_state%hvisc%bound_kh = cfg%ocean%hvisc%bound_kh
      ocean_state%hvisc%no_slip = cfg%ocean%hvisc%no_slip
      ! Fail-loud guard (2026-07-28 dt=800 forensics): bound_kh clamps the
      ! effective per-cell viscosity to
      ! kh_max = bound_coef·0.125/(dt·(1/dx² + 1/dy²)), so a configured
      ! nu_h far above kh_max is silently unreachable and the run is
      ! under-damped relative to its nml intent — the dt=800 blow-up
      ! family was exactly this (bound_coef=0.15 capped the effective nu
      ! at ~177 m²/s against nu_h=10000).  Warning, not error: flow-aware
      ! closures legitimately over-provision constant floors/ceilings.
      ! CARTESIAN only: `grid%dx`/`grid%dy` are the cell size there, but a
      ! placeholder on a curvilinear grid (1 m on a supergrid, degrees on a
      ! spherical sector) — which made the estimate ~1e-5 m2/s and the
      ! warning fire, falsely, on every global run.  Curvilinear grids are
      ! covered per wet cell by the viscous-CFL stability audit.
      if (cfg%ocean%hvisc%bound_kh .and. cfg%dt_fixed > 0.0_wp &
          .and. parse_grid_config(cfg%ocean%grid%grid_config) == GRID_CONFIG_CARTESIAN &
          .and. grid%dx > 0.0_wp .and. grid%dy > 0.0_wp) then
         block
            real(wp) :: kh_max_est
            kh_max_est = cfg%ocean%hvisc%bound_coef*0.125_wp/ &
                         (cfg%dt_fixed*(1.0_wp/grid%dx**2 + 1.0_wp/grid%dy**2))
            if (cfg%ocean%hvisc%nu_h > 10.0_wp*kh_max_est .and. compute_rank == 0) then
               call logger%warning("hvisc: nu_h = "//to_string(cfg%ocean%hvisc%nu_h)// &
                                   " m2/s is unreachable — bound_kh clamps the effective "// &
                                   "viscosity to ~"//to_string(kh_max_est)// &
                                   " m2/s at this dx/dt (bound_coef = "// &
                                   to_string(cfg%ocean%hvisc%bound_coef)// &
                                   ").  The run will be under-damped relative to the nml "// &
                                   "intent; raise bound_coef (MOM6 default 0.8), reduce dt, "// &
                                   "or lower nu_h to what the clamp admits.")
            end if
         end block
      end if
      ! Anisotropic viscosity (Smith & McWilliams 2003) — stress-tensor
      ! path only.  Default kh_aniso=0 ⇒ isotropic, bit-identical.  The
      ! constant direction tensor is precomputed once from aniso_dir
      ! (mode 0 = grid-relative).  `aniso_mode /= 0` is rejected at
      ! configure (validate_config); the else below is a defensive guard.
      ocean_state%hvisc%kh_aniso = cfg%ocean%hvisc%kh_aniso
      if (cfg%ocean%hvisc%aniso_mode == 0) then
         call ocean_hvisc_set_aniso_direction(ocean_state%hvisc, &
                                              cfg%ocean%hvisc%aniso_dir(1), &
                                              cfg%ocean%hvisc%aniso_dir(2))
      else
         call ocean_hvisc_set_aniso_direction(ocean_state%hvisc, 1.0_wp, 0.0_wp)
      end if

      ! Coriolis: beta = 0 → uniform f_0 = cfg%coriolis_f; non-zero beta engages
      ! f(y) = f_0 + beta*(y - y_ref).  Filled by the single generator-driven
      ! routine (D7) — beta_plane is bit-identical to the legacy setter; the
      ! `f_0`/`beta` diagnostics are recorded here (the generator only fills
      ! the array).
      ocean_state%coriolis_adv%f_0 = cfg%coriolis_f
      ocean_state%coriolis_adv%beta = cfg%ocean%topo%coriolis_beta
      call fill_coriolis_corner(cfg, ocean_state%metrics, grid, &
                                ocean_state%coriolis_adv%f_corner)
      ! Tripolar: fold the STATIC f_corner once at configure (it is read
      ! straight from f_corner in the BT substep + Coriolis-adv, never
      ! re-wrapped per step).  f is a SCALAR under the fold — the reflected
      ! corner sits at the SAME latitude, so it COPIES (negate=.false.), no
      ! sign flip.  Periodic-x of the corner columns first (Appendix A
      ! ordering), then the north fold.  No-op for non-tripolar runs.
      if (trim(cfg%ocean%bc%north) == "tripolar_fold") then
         call fill_f_corner_seam_ghosts(grid, ocean_state%coriolis_adv%f_corner)
      end if
      ocean_state%coriolis_adv%pv_variant = parse_pv_variant(cfg%ocean%coriolis%form)
      ocean_state%coriolis_adv%use_hk_correction = &
         (ocean_state%coriolis_adv%pv_variant == PV_VARIANT_SADOURNY_HK)
      ocean_state%coriolis_adv%pv_adv_scheme = &
         parse_pv_adv_scheme(cfg%ocean%coriolis%pv_adv_scheme)
      ! Coastal lateral BC (land mask, C1): free-slip default, shared knob.
      ocean_state%coriolis_adv%no_slip = cfg%ocean%hvisc%no_slip
      ! Mass-consistent CorAdCalc (MOM6 parity); config fail-loud
      ! guarantees form="sadourny_energy" + split_scheme="pred_corr".
      ocean_state%coriolis_adv%state_fluxes = cfg%ocean%coriolis%use_state_fluxes
      ! BOUND_CORIOLIS velocity-form clamp (energy scheme only; config
      ! fail-loud guarantees form="sadourny_energy").
      ocean_state%coriolis_adv%bound_coriolis = cfg%ocean%coriolis%bound_coriolis
      ! PV corner-thickness construction (energy scheme only, config fail-loud).
      if (trim(cfg%ocean%coriolis%corner_h) == "mom6_area") then
         ocean_state%coriolis_adv%corner_h_variant = CORNER_H_MOM6_AREA
      else
         ocean_state%coriolis_adv%corner_h_variant = CORNER_H_CELL_MEAN
      end if
      if (compute_rank == 0) then
         call logger%info("Coriolis form:    "//trim(cfg%ocean%coriolis%form))
         if (cfg%ocean%coriolis%use_state_fluxes) then
            call logger%info("Coriolis CorAdv:  mass-consistent state fluxes "// &
                             "in the mom6 corrector (MOM6 uh/vh parity)")
         end if
         if (cfg%ocean%coriolis%bound_coriolis) then
            call logger%info("Coriolis bound:   BOUND_CORIOLIS ON — energy-scheme accel "// &
                             "clamped to the (f+zeta)*v velocity-form range (MOM6)")
         end if
         if (trim(cfg%ocean%coriolis%corner_h) == "mom6_area") then
            call logger%info("Coriolis corner:  PV corner-h = mom6_area (MOM6 "// &
                             "Area_q/(hArea_q+vol_neglect); round-off-identical to "// &
                             "cell_mean above H_MIN_PV)")
         end if
      end if
   end subroutine configure_ocean_forcing

   subroutine fill_f_corner_seam_ghosts(grid, f_corner)
      !! Periodic-x wrap + north-fold (scalar copy) of the static corner
      !! Coriolis array for a tripolar grid.  f is reflection-invariant
      !! (same latitude at the conjugate corner), so negate=.false.
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(inout) :: f_corner(:, :)
      integer :: ng, ni, i, j
      ng = grid%nghost
      ni = grid%nx_phys
      ! Periodic-x ghost columns (corner/Bu face-type: physical i=ng+1..ng+ni+1).
      do j = 1, size(f_corner, 2)
         do i = 1, ng
            f_corner(i, j) = f_corner(i + ni, j)
            f_corner(ng + ni + 1 + i, j) = f_corner(ng + 1 + i, j)
         end do
      end do
      ! North fold (copy).
      call fold_north_corner(f_corner, grid%nx_total + 1, grid%ny_total + 1, &
                             grid%nx_phys, grid%ny_phys, grid%nghost, negate=.false.)
   end subroutine fill_f_corner_seam_ghosts

   subroutine configure_ocean_drag(cfg, ocean_state, compute_rank, ierr)
      !! Bottom-drag variant + coefficients, the continuity PPM positivity
      !! guard, and the BT-budget diagnostic probe — plus their rank-0 log
      !! lines.  All knobs default off/zero (bit-identical to pre-knob nmls).
      type(config_t), intent(in) :: cfg
      type(ocean_state_t), intent(inout) :: ocean_state
      integer, intent(in) :: compute_rank
      integer, intent(out), optional :: ierr
         !! Non-zero on an isopycnal-floor configuration conflict when
         !! present; absent behaves as today (`error stop`).

      ! EOS variant FIRST: this is the earliest configure_ocean_* call, so
      ! `ocean_state%eos%variant` is set before the EOS-handle is copied into
      ! the vmix / EPBL / kappa-shear slots (in configure_ocean_lateral /
      ! _epbl / _kappa_shear).  The device-callable gate + the FV_WRIGHT
      ! cross-check run later, in configure_ocean_pgf (post pgf-variant parse).
      ocean_state%eos%variant = parse_eos_variant(cfg%ocean%eos%eos)
      ! Potential-density reference pressure.  Assigned HERE, in the
      ! earliest configure_ocean_* call, so the value is already on the
      ! handle when it is copied into the vmix / EPBL / kappa-shear slots
      ! further down — and before `enter_data`, so the mapped copies carry
      ! it.  `eos_t` is a flat POD read by value into the `_impl` calls, so
      ! this host assignment owes no `!$acc update device` under
      ! mem:separate (same contract as `rho0`/`rho_ref` on the PGF slot).
      ! Horizontally uniform by design — see the knob's FORD docstring.
      ocean_state%eos%p_ref = cfg%ocean%eos%p_ref
      ! Freezing-point (liquidus) coefficient set — same site, same
      ! reasons, and one more: the ICE slot reads the god-state handle
      ! directly (`engine%state%eos` is what `ice_frazil_accumulate`,
      ! `ice_frazil_uptake` and `ice_compute_basal_flux` are handed), so
      ! landing the set here puts it on the ONE handle every liquidus
      ! consumer sees, before any copy is taken and before `enter_data`.
      ! An unrecognised string is already fatal in `validate_config`;
      ! `eos_apply_tfreeze_set` leaves the handle untouched for it.
      call eos_apply_tfreeze_set(ocean_state%eos, &
                                 parse_tfreeze_set(cfg%ocean%eos%tfreeze_set))
      if (compute_rank == 0) then
         call logger%info("EOS variant:      "//trim(cfg%ocean%eos%eos))
         if (parse_tfreeze_set(cfg%ocean%eos%tfreeze_set) /= TFREEZE_SET_SEAICE) then
            call logger%info("Liquidus set:     "//trim(cfg%ocean%eos%tfreeze_set)// &
                             " (ISOMIP+ / Asay-Davis et al. 2016 Table 4; the "// &
                             "default 'seaice' SIS2 set is ~0.03 degC warmer at S=34.5)")
         end if
      end if

      ocean_state%bdrag%variant = parse_bdrag_variant(cfg%ocean%bdrag%form)
      ocean_state%bdrag%c_drag = cfg%ocean%bdrag%cd
      ocean_state%bdrag%r_linear = cfg%ocean%bdrag%r
      ocean_state%bdrag%hbbl = cfg%ocean%bdrag%hbbl
      ocean_state%bdrag%drag_bg_vel = cfg%ocean%bdrag%bg_vel
      ocean_state%bdrag%bbl_thick_min = cfg%ocean%bdrag%bbl_thick_min
      ocean_state%bdrag%bed_factor = cfg%ocean%bdrag%bed_factor
      ocean_state%bdrag%channel_drag = cfg%ocean%bdrag%channel_drag
      ocean_state%bdrag%cdrag_side = cfg%ocean%bdrag%cdrag_side
      ocean_state%bdrag%implicit = cfg%ocean%bdrag%implicit

      ! Continuity PPM positivity guard (MOM6 PPM_limit_pos analogue).
      ocean_state%continuity%use_ppm_limit_pos = cfg%ocean%continuity%ppm_limit_pos
      ocean_state%continuity%renorm_consistent_flux = cfg%ocean%continuity%renorm_consistent_flux
      ! Wet/dry composition: the Newton uhbt renormalisation's donor
      ! re-pick + CFL bracket destabilise the drying front (h_layer goes
      ! negative, test_ocean_wetdry_driver) — wet/dry keeps the legacy
      ! single-step form it was validated on.
      ocean_state%continuity%renorm_legacy_single_step = cfg%ocean%wetdry%enable
      ocean_state%continuity%h_min = cfg%ocean%continuity%h_min
      ocean_state%continuity%vol_cfl = cfg%ocean%continuity%vol_cfl
      if (compute_rank == 0 .and. cfg%ocean%continuity%ppm_limit_pos) then
         call logger%info("Continuity PPM:   positivity guard ON (MOM6 PPM_limit_pos), "// &
                          "h_min = "//to_string(cfg%ocean%continuity%h_min)//" m")
      end if
      ! Positive-definite split continuity (MOM6-prevention + Roundabout-conservation).
      ! h_lim = angstrom_h on VCOORD_LAGRANGIAN, else 0 (⇒ P1 floor inert on
      ! non-Lagrangian coords).  Default off ⇒ bit-identical.
      ocean_state%continuity%positive_definite = cfg%ocean%continuity%positive_definite
      if (parse_ocean_vcoord_type(cfg%vcoord_type) == VCOORD_LAGRANGIAN) then
         ocean_state%continuity%h_lim = cfg%ocean%isopycnal%angstrom_h
      else
         ocean_state%continuity%h_lim = 0.0_wp
      end if
      if (compute_rank == 0 .and. cfg%ocean%continuity%positive_definite) then
         call logger%info("Continuity PD:    positive_definite ON (h>=h_lim, zero mass "// &
                          "created), h_lim = "//to_string(ocean_state%continuity%h_lim)//" m")
      end if
      ! Phase-1 Lagrangian minimum-thickness floor (MOM6 Angstrom_H analogue).
      ! Only passed to the h-update kernels when the active vcoord is
      ! VCOORD_LAGRANGIAN (gated in rdb_ocean_dyn at the call site). 0 = off.
      ocean_state%continuity%angstrom_h = cfg%ocean%isopycnal%angstrom_h
      if (compute_rank == 0 .and. cfg%ocean%isopycnal%angstrom_h > 0.0_wp) then
         call logger%info("Isopycnal floor:  angstrom_h = "// &
                          to_string(cfg%ocean%isopycnal%angstrom_h)//" m (Lagrangian only)")
      end if
      ! Conservative minimum-thickness mode: replaces the injecting floor with a
      ! per-column borrow (rdb_ocean_min_thickness), invoked from the dyn
      ! continuity site. Fail-loud: needs angstrom_h > 0 AND VCOORD_LAGRANGIAN.
      ocean_state%continuity%conservative_floor = cfg%ocean%isopycnal%conservative_floor
      if (cfg%ocean%isopycnal%conservative_floor) then
         if (cfg%ocean%isopycnal%angstrom_h <= 0.0_wp) then
            call fail("ocean_isopycnal_nml: conservative_floor=.true. requires "// &
                      "angstrom_h > 0 (nothing to floor to).", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
         if (parse_ocean_vcoord_type(cfg%vcoord_type) /= VCOORD_LAGRANGIAN) then
            call fail("ocean_isopycnal_nml: conservative_floor=.true. is only "// &
                      "valid with vcoord_type='lagrangian' (isopycnal path).", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
         if (compute_rank == 0) then
            call logger%info("Isopycnal floor:  conservative_floor = ON "// &
                             "(per-column borrow, mass/momentum/tracer conserving)")
         end if
      end if
      ! Phase-2/3 Lagrangian grounding-stability knobs on ocean_dyn_t.
      ! Gates tested at the dyn call sites; bit-identical when off.
      ocean_state%dyn%angstrom_h = cfg%ocean%isopycnal%angstrom_h
      ocean_state%dyn%reset_vanished_u = cfg%ocean%isopycnal%reset_vanished_u
      ocean_state%dyn%cfl_ignore_vanished = cfg%ocean%isopycnal%cfl_ignore_vanished
      ocean_state%dyn%check_h_positive = cfg%ocean%isopycnal%check_h_positive
      if (compute_rank == 0) then
         if (cfg%ocean%isopycnal%reset_vanished_u) then
            call logger%info("Isopycnal reset:  reset_vanished_u = ON (Lagrangian only)")
         end if
         if (cfg%ocean%isopycnal%cfl_ignore_vanished) then
            call logger%info("Isopycnal CFL:    cfl_ignore_vanished = ON (Lagrangian only)")
         end if
      end if
      if (compute_rank == 0 .and. cfg%ocean%continuity%vol_cfl) then
         call logger%info("Continuity PPM:   swept-volume face thickness ON (MOM6 vol_CFL)")
      end if

      ! Ghost-band poison debug knob.
      ocean_state%dyn%poison_ghosts = cfg%ocean%mpi%poison_ghosts

      ! BT-budget diagnostic probe — heavy when on (D->H + 16 prints/step).
      ocean_state%dyn%debug_bt_budget = cfg%ocean%debug%budget
      if (compute_rank == 0 .and. cfg%ocean%debug%budget) then
         call logger%info("BT-budget probe:  ENABLED (per-step per-region "// &
                          "power decomposition)")
      end if

      ! KE attribution meter — heavy when on (serialises the apply chain).
      ocean_state%dyn%ke_probe%enable = cfg%ocean%debug%ke_attr
      ocean_state%dyn%ke_probe%start_step = cfg%ocean%debug%ke_attr_start_step
      ocean_state%dyn%ke_probe%end_step = cfg%ocean%debug%ke_attr_end_step
      ocean_state%dyn%chksum_probe%enable = cfg%ocean%debug%chksum
      ocean_state%dyn%chksum_probe%start_step = cfg%ocean%debug%chksum_start_step
      ocean_state%dyn%chksum_probe%end_step = cfg%ocean%debug%chksum_end_step
      ocean_state%dyn%chksum_probe%interior = cfg%ocean%debug%chksum_interior
      if (compute_rank == 0 .and. cfg%ocean%debug%ke_attr) then
         call logger%info("KE-attr meter:    ENABLED (per-segment layer-KE "// &
                          "rows, steps "//to_string(cfg%ocean%debug%ke_attr_start_step)// &
                          " to "//to_string(cfg%ocean%debug%ke_attr_end_step)//")")
      end if
      if (compute_rank == 0) then
         call logger%info("Bottom drag:      "//trim(cfg%ocean%bdrag%form)// &
                          "  C_d="//to_string(cfg%ocean%bdrag%cd)// &
                          "  r="//to_string(cfg%ocean%bdrag%r)//" 1/s")
         if (cfg%ocean%bdrag%hbbl > 0.0_wp) then
            call logger%info("Bottom drag BBL:  HBBL="//to_string(cfg%ocean%bdrag%hbbl)// &
                             " m  bg_vel="//to_string(cfg%ocean%bdrag%bg_vel)// &
                             " m/s  thick_min="//to_string(cfg%ocean%bdrag%bbl_thick_min)//" m")
         else
            call logger%info("Bottom drag BBL:  bed-layer only (HBBL=0)")
            ! Bed-only mode drags layer k = 1.  On a coordinate whose bed-side
            ! layers VANISH (z_fixed, zstar_full) k = 1 is an inert filler in
            ! every column shallower than the deepest nominal interface, so
            ! almost the whole domain runs with NO bottom drag — measured on
            ! the global 1-degree case (validation_examples/ocean/global_1deg).
            ! HBBL mode accumulates thickness from the bed up, skips the
            ! fillers and reaches the live bottom layer.
            if ((cfg%ocean%bdrag%cd > 0.0_wp .or. cfg%ocean%bdrag%r > 0.0_wp) .and. &
                (parse_ocean_vcoord_type(cfg%vcoord_type) == VCOORD_Z_FIXED .or. &
                 parse_ocean_vcoord_type(cfg%vcoord_type) == VCOORD_ZSTAR_FULL)) then
               call logger%warning("&ocean_bdrag_nml hbbl = 0 (bed-layer-only drag) under "// &
                                   "vcoord_type = '"//trim(cfg%vcoord_type)//"': the drag "// &
                                   "acts on layer k = 1, an inert filler in every column "// &
                                   "shallower than the deepest nominal layer, so those "// &
                                   "columns get NO bottom drag.  Set hbbl > 0 (MOM6 "// &
                                   "OM4/OM_1deg: HBBL = 10 m, bg_vel = 0.1 m/s).")
            end if
         end if
         if (cfg%ocean%bdrag%channel_drag) then
            call logger%info("Channel drag:     ON  cdrag_side="// &
                             to_string(cfg%ocean%bdrag%cdrag_side))
         end if
      end if
      if (present(ierr)) ierr = OCEAN_STATUS_OK
   end subroutine configure_ocean_drag

   subroutine configure_ocean_hdiff(cfg, ocean_state, compute_rank)
      !! Along-coordinate tracer Laplacian coefficient (`&ocean_hdiff_nml
      !! kappa_h`).  Scalar copy onto `ocean_state%hdiff_tracer`, mapped
      !! with its parent slot at `ocean_state_enter_data` — no explicit
      !! `!$acc update` needed as long as this runs before that (it does;
      !! see `rdb_driver.F90`).  Default `kappa_h = 0.0` leaves the
      !! kernel's short-circuit intact ⇒ bit-identical.
      type(config_t), intent(in) :: cfg
      type(ocean_state_t), intent(inout) :: ocean_state
      integer, intent(in) :: compute_rank

      ocean_state%hdiff_tracer%kappa_h = cfg%ocean%hdiff%kappa_h
      if (compute_rank == 0 .and. cfg%ocean%hdiff%kappa_h > 0.0_wp) then
         call logger%info("Tracer hdiff:     along-coordinate kappa_h = "// &
                          to_string(cfg%ocean%hdiff%kappa_h)//" m^2/s")
      end if
   end subroutine configure_ocean_hdiff

   subroutine configure_ocean_vmix(cfg, ocean_state, compute_rank, ierr)
      !! Thermodynamics on/off, velocity-truncation clamp (MAXVEL), DIRECT_STRESS
      !! surface-stress distribution, KV_ML_INVZ2 surface-band viscosity,
      !! HARMONIC_VISC face-thickness mean, and the DT_THERM thermo/tracer
      !! cadence — with their rank-0 log lines.
      type(config_t), intent(in) :: cfg
      type(ocean_state_t), intent(inout) :: ocean_state
      integer, intent(in) :: compute_rank
      integer, intent(out), optional :: ierr
         !! Non-zero on a vertical-mixing configuration conflict when
         !! present; absent behaves as today (`error stop`).

      ocean_state%dyn%enable_thermodynamics = cfg%ocean%thermo%enable_thermodynamics

      ! Velocity-truncation clamp (MOM6 MAXVEL); 0 = off.
      ocean_state%dyn%maxvel = cfg%ocean%pgf%maxvel
      if (compute_rank == 0) then
         if (cfg%ocean%pgf%maxvel > 0.0_wp) then
            call logger%info("Velocity clamp:   "//to_string(cfg%ocean%pgf%maxvel)//" m/s")
         end if
      end if

      ! Advective-CFL velocity truncation (E7, MOM6 CFL_trunc); 0 = off.
      ocean_state%dyn%cfl_trunc = cfg%ocean%pgf%cfl_trunc
      if (compute_rank == 0) then
         if (cfg%ocean%pgf%cfl_trunc > 0.0_wp) then
            call logger%info("CFL truncation:   threshold "//to_string(cfg%ocean%pgf%cfl_trunc))
         end if
      end if

      ! DIRECT_STRESS — wind stress spread over the top hmix_stress metres.
      ocean_state%surface_stress%direct_stress = cfg%ocean%vmix%direct_stress
      ocean_state%surface_stress%hmix_stress = cfg%ocean%vmix%hmix_stress
      if (compute_rank == 0) then
         if (cfg%ocean%vmix%direct_stress) then
            call logger%info("Surface stress:   DIRECT (over top "// &
                             to_string(cfg%ocean%vmix%hmix_stress)//" m)")
         end if
      end if

      ! KV_ML_INVZ2 — extra surface-band vertical viscosity (MOM6).
      ocean_state%vmix%kv_ml_invz2 = cfg%ocean%vmix%kv_ml_invz2
      ocean_state%vmix%hmix_fixed = cfg%ocean%vmix%hmix_fixed
      if (compute_rank == 0) then
         if (cfg%ocean%vmix%kv_ml_invz2 > 0.0_wp) then
            call logger%info("KV_ML_INVZ2:      "//to_string(cfg%ocean%vmix%kv_ml_invz2)// &
                             " m²/s in top "//to_string(cfg%ocean%vmix%hmix_fixed)//" m")
         end if
      end if

      ! C11 assembly gate — ceilings / smoothing / guard.  Backgrounds
      ! stay on the slot defaults (= pp81_*_bg, set in ocean_vmix_init)
      ! so the floor is a no-op for the closure path; only the ceilings /
      ! smoothing / guard are config-exposed.
      ocean_state%vmix%kv_max = cfg%ocean%vmix%kv_max
      ocean_state%vmix%kd_max = cfg%ocean%vmix%kd_max
      ocean_state%vmix%kd_smooth_iterations = cfg%ocean%vmix%kd_smooth_iterations
      ocean_state%vmix%vmix_guard = cfg%ocean%vmix%vmix_guard

      ! C7 Bryan-Lewis depth-varying background (default off = bit-identical).
      ocean_state%vmix%bkgnd_profile = cfg%ocean%vmix%bkgnd_profile
      ocean_state%vmix%bkgnd_kd_sfc = cfg%ocean%vmix%bkgnd_kd_sfc
      ocean_state%vmix%bkgnd_kd_deep = cfg%ocean%vmix%bkgnd_kd_deep
      ocean_state%vmix%bkgnd_z0 = cfg%ocean%vmix%bkgnd_z0
      ocean_state%vmix%bkgnd_delta = cfg%ocean%vmix%bkgnd_delta
      ocean_state%vmix%bkgnd_prandtl = cfg%ocean%vmix%bkgnd_prandtl
      ocean_state%vmix%bkgnd_henyey = cfg%ocean%vmix%bkgnd_henyey
      ocean_state%vmix%bkgnd_kd_min = cfg%ocean%vmix%bkgnd_kd_min
      ocean_state%vmix%bkgnd_henyey_n0_2omega = cfg%ocean%vmix%bkgnd_henyey_n0_2omega
      ocean_state%vmix%bkgnd_henyey_max_lat = cfg%ocean%vmix%bkgnd_henyey_max_lat
      ! Backstop for the two `validate_config` guards (Bryan-Lewis and Henyey
      ! mutually exclusive; non-cartesian grid_config required) — reached only
      ! by a caller that bypassed validation, since the namelist path stops
      ! there first.  Repeated here because a Henyey run with either
      ! precondition violated is silent (wrong background field), not noisy.
      if (bkgnd_henyey_conflicts_profile(cfg%ocean%vmix%bkgnd_henyey, &
                                         cfg%ocean%vmix%bkgnd_profile)) then
         call fail("ocean_vmix_nml: bkgnd_henyey and bkgnd_profile are "// &
                   "mutually exclusive background schemes — the Henyey factor "// &
                   "scales the SCALAR background, it does not compose with the "// &
                   "Bryan-Lewis depth profile", ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if
      if (cfg%ocean%vmix%bkgnd_henyey .and. &
          trim(cfg%ocean%grid%grid_config) == "cartesian") then
         call fail("ocean_vmix_nml: bkgnd_henyey requires a non-cartesian "// &
                   "grid_config — geolatT is identically zero on cartesian, so "// &
                   "every column would take the equatorial L(0 deg)=0 and the "// &
                   "background would flatten to a uniform bkgnd_kd_min", ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if
      if (compute_rank == 0) then
         if (cfg%ocean%vmix%bkgnd_profile) then
            call logger%info("Background mixing: Bryan-Lewis profile (Kd_sfc="// &
                             to_string(cfg%ocean%vmix%bkgnd_kd_sfc)//" Kd_deep="// &
                             to_string(cfg%ocean%vmix%bkgnd_kd_deep)//" z0="// &
                             to_string(cfg%ocean%vmix%bkgnd_z0)//" m)")
         end if
         if (cfg%ocean%vmix%bkgnd_henyey) then
            call logger%info("Background mixing: Henyey IGW latitude factor on the "// &
                             "scalar background (N0_2Omega="// &
                             to_string(cfg%ocean%vmix%bkgnd_henyey_n0_2omega)// &
                             " max_lat="//to_string(cfg%ocean%vmix%bkgnd_henyey_max_lat)// &
                             " degN Kd_min="// &
                             to_string(vmix_resolve_kd_min(cfg%ocean%vmix%bkgnd_kd_min, &
                                                           cfg%ocean%vmix%pp81_kappa_bg))// &
                             " m^2/s)")
            ! `pp81_kappa_bg` (not `ocean_state%vmix%kt_bg`) because
            ! `seed_backgrounds` — which re-derives kt_bg from it — runs
            ! later in this same routine; reading the slot here would log
            ! the stale type default whenever the user retuned kappa_bg.
         end if
      end if

      ! smooth_scratch is allocated ONLY when smoothing is requested —
      ! here, after kd_smooth_iterations is known, before enter_data.
      ! The enter_data path guards on allocated() so the not-allocated
      ! case (smoothing off) is handled correctly.
      if (cfg%ocean%vmix%kd_smooth_iterations > 0) then
         if (.not. allocated(ocean_state%vmix%smooth_scratch)) then
            associate (kv => ocean_state%vmix%kv)
               allocate (ocean_state%vmix%smooth_scratch( &
                         size(kv, 1), size(kv, 2), size(kv, 3)), source=0.0_wp)
            end associate
         end if
      end if

      if (compute_rank == 0) then
         if (cfg%ocean%vmix%kd_max < huge(1.0_wp) .or. cfg%ocean%vmix%kv_max < huge(1.0_wp)) then
            call logger%info("Vmix assembly:    kv_max="//to_string(cfg%ocean%vmix%kv_max)// &
                             " kd_max="//to_string(cfg%ocean%vmix%kd_max)//" m²/s")
         end if
         if (cfg%ocean%vmix%kd_smooth_iterations > 0) then
            call logger%info("Vmix assembly:    "// &
                             to_string(cfg%ocean%vmix%kd_smooth_iterations)// &
                             " 1-2-1 smoothing pass(es) on kv/kt (wet-mask aware)")
         end if
         if (cfg%ocean%vmix%vmix_guard) then
            call logger%info("Vmix assembly:    negative/NaN guard ON")
         end if
      end if

      ! HARMONIC_VISC — vdiff face-thickness mean (better at thin layers).
      ocean_state%vdiff%use_harmonic = cfg%ocean%vmix%harmonic_visc
      ! MOM6 HARMONIC_VISC parity for the momentum face thickness (hvel) +
      ! arithmetic h_shear.  Off by default => historical arithmetic h_u.
      ocean_state%vdiff%hvel_mom6 = cfg%ocean%vdiff%hvel_mom6
      ocean_state%vdiff%hbbl_visc = cfg%ocean%vdiff%hbbl_visc
      ! MOM6 bottomdraglaw coupling parity: kv_bbl botfn glue + piston bed
      ! drag (PGF_BUG.md §9).  Preconditions (hvel_mom6, implicit_drag,
      ! linear bdrag form) validated at configure.
      ocean_state%vdiff%bbl_glue = cfg%ocean%vdiff%bbl_glue
      ocean_state%vdiff%bbl_piston = cfg%ocean%vdiff%bbl_piston
      ocean_state%vdiff%hvel_upwind = cfg%ocean%vdiff%hvel_upwind
      if (compute_rank == 0 .and. cfg%ocean%vdiff%hvel_mom6 &
          .and. .not. cfg%ocean%vdiff%hvel_upwind) then
         call logger%info("vdiff: hvel near-bed upwind blend OFF (pure harmonic hvel)")
      end if
      if (compute_rank == 0 .and. cfg%ocean%vdiff%hvel_mom6) then
         call logger%info("vdiff: MOM6 hvel (harmonic + near-bed upwind blend, "// &
                          "arithmetic h_shear) ON")
      end if
      if (compute_rank == 0 .and. cfg%ocean%vdiff%bbl_glue) then
         call logger%info("vdiff: MOM6 bottomdraglaw BBL coupling glue ON "// &
                          "(kv_bbl botfn interfaces + piston bed row)")
      end if
      if (compute_rank == 0) then
         if (cfg%ocean%vmix%harmonic_visc) then
            call logger%info("Vertical visc:    harmonic-mean face thickness")
         end if
      end if

      ! Implicit stress/drag folding into the vdiff tridiagonal
      ! (`&ocean_vdiff_nml`).  When on, the corresponding explicit
      ! pre-solve apply is gated off in run_stage (see rdb_ocean_dyn), and
      ! the bottom-drag slot fills the bed-layer Rayleigh-rate field
      ! `lambda_bot_u/v` that the vdiff diagonal consumes.  Mutual
      ! exclusions validated at configure (validate_config).
      ocean_state%vdiff%implicit_stress = cfg%ocean%vdiff%implicit_stress
      ocean_state%vdiff%implicit_drag = cfg%ocean%vdiff%implicit_drag
      ocean_state%bdrag%implicit_fold = cfg%ocean%vdiff%implicit_drag
      if (compute_rank == 0) then
         if (cfg%ocean%vdiff%implicit_stress) then
            call logger%info("Vertical visc:    wind stress folded into vdiff "// &
                             "surface (k=nz) RHS (implicit_stress)")
         end if
         if (cfg%ocean%vdiff%implicit_drag) then
            call logger%info("Vertical visc:    bottom drag folded into vdiff "// &
                             "bed (k=1) diagonal (implicit_drag)")
         end if
      end if

      ! DT_THERM ratio — thermo/tracer kernels fire every Nth outer step.
      ocean_state%dyn%dt_therm_ratio = max(1, cfg%ocean%vmix%dt_therm_ratio)
      if (compute_rank == 0) then
         if (ocean_state%dyn%dt_therm_ratio > 1) then
            call logger%info("Thermo ratio:     "// &
                             to_string(ocean_state%dyn%dt_therm_ratio)// &
                             "·dt (DT_THERM = ratio · dt_dynamic)")
         end if
      end if

      ! DT_TRACER_ADVECT ratio — horizontal tracer advect fires every
      ! Nth outer step over the accumulated face transports.  Must be
      ! >= 1 and dt_therm_ratio must be an integer multiple of it so the
      ! ALE remap (which fires at DT_THERM cadence) never lands inside an
      ! open accumulation window.  Fail loud, host-side.
      if (cfg%ocean%vmix%dt_tracer_advect_ratio < 1) then
         call fail("ocean_vmix_nml: dt_tracer_advect_ratio must be >= 1 (got "// &
                   to_string(cfg%ocean%vmix%dt_tracer_advect_ratio)//")", ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if
      ! `< 1` already error-stopped above, so the value is >= 1 here.
      ocean_state%dyn%dt_tracer_advect_ratio = cfg%ocean%vmix%dt_tracer_advect_ratio
      if (.not. ocean_dt_tracer_advect_ratios_ok(ocean_state%dyn%dt_therm_ratio, &
                                                 ocean_state%dyn%dt_tracer_advect_ratio)) then
         call fail("ocean_vmix_nml: dt_therm_ratio ("// &
                   to_string(ocean_state%dyn%dt_therm_ratio)// &
                   ") must be an integer multiple of dt_tracer_advect_ratio ("// &
                   to_string(ocean_state%dyn%dt_tracer_advect_ratio)// &
                   ") so the ALE remap never fires mid-accumulation-window", ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if
      if (compute_rank == 0) then
         if (ocean_state%dyn%dt_tracer_advect_ratio > 1) then
            call logger%info("Tracer advect:    "// &
                             to_string(ocean_state%dyn%dt_tracer_advect_ratio)// &
                             "·dt (DT_TRACER_ADVECT = ratio · dt_dynamic)")
         end if
      end if

      ! Q6: windowed-drain tracer face reconstruction (ppm | weno5/7/9).
      ! The per-rung nghost minimum is validated fail-loud in validate_config
      ! (cfg%nghost is the ocean grid nghost); here we only map the parsed
      ! code onto the continuity slot.  ppm (default) ⇒ bit-identical.
      ocean_state%continuity%tracer_recon = &
         parse_tracer_recon(trim(cfg%ocean%vmix%tracer_recon))
      if (compute_rank == 0 .and. &
          ocean_state%continuity%tracer_recon /= TRACER_RECON_PPM) then
         call logger%info("Tracer recon:     "//trim(cfg%ocean%vmix%tracer_recon)// &
                          " (windowed-drain WENO-Z swept-average face)")
      end if
      if (compute_rank == 0) then
         if (cfg%ocean%thermo%enable_thermodynamics) then
            call logger%info("Thermodynamics:   on")
         else
            call logger%info("Thermodynamics:   OFF (adiabatic — EOS/tracer kernels skipped)")
         end if
      end if
      if (present(ierr)) ierr = OCEAN_STATUS_OK
   end subroutine configure_ocean_vmix

   subroutine configure_ocean_tracers(cfg, ocean_state, compute_rank)
      !! `&ocean_tracers_nml` scalar knobs (PR-7): the ideal-age Dirichlet
      !! surface value and its vintage-mode exponential growth rate.
      !! `enable_ideal_age` itself is latched earlier by
      !! `ocean_state_copy_config` (before `init(grid)`, since it gates
      !! tracer-slot allocation) — these two scalars gate nothing, so
      !! they land here in the post-init configure pass.  Both default
      !! to 0 ⇒ bit-identical.
      type(config_t), intent(in) :: cfg
      type(ocean_state_t), intent(inout) :: ocean_state
      integer, intent(in) :: compute_rank

      ocean_state%dyn%ideal_age_young_val = cfg%ocean%tracers%ideal_age_young_val
      ocean_state%dyn%ideal_age_sfc_growth_rate = cfg%ocean%tracers%ideal_age_sfc_growth_rate
      if (compute_rank == 0 .and. cfg%ocean%tracers%enable_ideal_age) then
         call logger%info("Ideal age:        surface reset after rk2_average + remap")
         if (cfg%ocean%tracers%ideal_age_sfc_growth_rate /= 0.0_wp) then
            call logger%info("Ideal age:        vintage mode, young_val="// &
                             to_string(cfg%ocean%tracers%ideal_age_young_val)// &
                             " growth_rate="// &
                             to_string(cfg%ocean%tracers%ideal_age_sfc_growth_rate)//" 1/s")
         end if
      end if
   end subroutine configure_ocean_tracers

   subroutine configure_ocean_lateral(cfg, ocean_state, grid, compute_rank, ierr)
      !! Flow-aware lateral-viscosity closure (Leith / Smagorinsky + biharmonic
      !! Smagorinsky_AH), the vertical coordinate (VCOORD_* code + z_fixed
      !! reference depth), and the PP81/KPP vertical-mixing switches — with
      !! their rank-0 log lines.
      type(config_t), intent(in) :: cfg
      type(ocean_state_t), intent(inout) :: ocean_state
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: compute_rank
      integer, intent(out), optional :: ierr
         !! Non-zero when any of the sub-closures configured here
         !! (EPBL/kappa-shear/tidal-mixing/conv/ddiff/wavespeed/Fox-Kemper)
         !! rejects the resolved configuration, when present; absent
         !! behaves as today (`error stop`).
      integer :: local_ierr

      ! Lateral closure + background/cap coefficients.  ah_bg defaults to the
      ! MOM6 floor max(nu_h, kh_vel_scale·dx) when cfg leaves it negative.
      ocean_state%lateral_mix%closure = parse_lateral_closure(cfg%ocean%hvisc%lateral_closure)
      ocean_state%lateral_mix%no_slip = cfg%ocean%hvisc%no_slip
      ocean_state%lateral_mix%c_smag = cfg%ocean%hvisc%c_smag
      ocean_state%lateral_mix%c_leith = cfg%ocean%hvisc%c_leith
      if (cfg%ocean%hvisc%ah_bg >= 0.0_wp) then
         ocean_state%lateral_mix%ah_bg = cfg%ocean%hvisc%ah_bg
      else
         ocean_state%lateral_mix%ah_bg = &
            max(cfg%ocean%hvisc%nu_h, cfg%ocean%hvisc%kh_vel_scale* &
                metrics_dx_min(ocean_state%metrics, grid))
      end if
      ocean_state%lateral_mix%ah_max = cfg%ocean%hvisc%ah_max
      ocean_state%lateral_mix%kh_vel_scale_live = cfg%ocean%hvisc%kh_vel_scale_live
      ocean_state%lateral_mix%smag_ah_active = cfg%ocean%hvisc%smag_ah
      ocean_state%lateral_mix%smag_bi_const = cfg%ocean%hvisc%smag_bi_const
      ocean_state%lateral_mix%c_leith_bi = cfg%ocean%hvisc%c_leith_bi
      ocean_state%lateral_mix%nu4_bg = cfg%ocean%hvisc%nu_4_bg
      ocean_state%lateral_mix%nu4_max = cfg%ocean%hvisc%nu_4_max
      ocean_state%lateral_mix%resoln_scaled_visc = cfg%ocean%hvisc%resoln_scaled_visc
      if (compute_rank == 0) then
         select case (ocean_state%lateral_mix%closure)
         case (LMIX_NONE)
            call logger%info("Lateral closure:  none (scalar ocean_nu_h)")
         case (LMIX_LEITH)
            call logger%info("Lateral closure:  Leith  C_L="//to_string(cfg%ocean%hvisc%c_leith)// &
                             "  ah_bg="//to_string(ocean_state%lateral_mix%ah_bg)// &
                             "  ah_max="//to_string(cfg%ocean%hvisc%ah_max)//" m²/s")
         case (LMIX_SMAGORINSKY)
            call logger%info("Lateral closure:  Smagorinsky  C_S="//to_string(cfg%ocean%hvisc%c_smag)// &
                             "  ah_bg="//to_string(ocean_state%lateral_mix%ah_bg)// &
                             "  ah_max="//to_string(cfg%ocean%hvisc%ah_max)//" m²/s")
         case (LMIX_LEITH_BIHARM)
            call logger%info("Lateral closure:  Leith-biharmonic  C_lb="//to_string(cfg%ocean%hvisc%c_leith_bi)// &
                             "  nu4_bg="//to_string(ocean_state%lateral_mix%nu4_bg)// &
                             "  nu4_max="//to_string(cfg%ocean%hvisc%nu_4_max)//" m⁴/s")
            ! PR-6: the c_leith_bi<=0 "inert" case is now a configure-time
            ! ABORT in validate_config (leith_biharm_is_inert), so setup is
            ! never reached with an inert leith_biharm — the warning here
            ! would be dead code.
         case (LMIX_BIHARMONIC)
            call logger%info("Lateral closure:  constant biharmonic  nu_4="//to_string(cfg%ocean%hvisc%nu_4)//" m⁴/s")
         case default
            ! LMIX_NONE handled above; any other code is rejected at
            ! configure (validate_config), so this is unreachable.
            call logger%info("Lateral closure:  "//trim(cfg%ocean%hvisc%lateral_closure))
         end select
         if (ocean_state%lateral_mix%smag_ah_active) then
            call logger%info("Biharmonic closure: Smagorinsky_AH  C_b="// &
                             to_string(cfg%ocean%hvisc%smag_bi_const)// &
                             "  nu4_bg="//to_string(cfg%ocean%hvisc%nu_4_bg)// &
                             "  nu4_max="//to_string(cfg%ocean%hvisc%nu_4_max)//" m⁴/s")
         end if
      end if

      ! Vertical coordinate: namelist string → VCOORD_* code; z_fixed needs a
      ! reference total depth.  Without this every vcoord_type silently ran
      ! Eulerian-z (no ALE remap).
      ocean_state%vcoord%coord_type = parse_ocean_vcoord_type(cfg%vcoord_type)
      ocean_state%vcoord%remap_method = parse_remap_method(cfg%remap_method)
      call configure_ocean_z_fixed_profile(cfg, ocean_state, compute_rank, log_it=.false.)
      ! Isopycnal (VCOORD_RHO) target densities: a uniform light->dense
      ! linspace from rho_target_light/dense (MOM6 ALE_COORDINATE_CONFIG=
      ! UNIFORM analogue).  `rho_target(0)` is the lightest (surface)
      ! interface, `rho_target(nz)` the densest (bed).
      ocean_state%vcoord%rho_ref_pressure = cfg%rho_ref_pressure
      if (allocated(ocean_state%vcoord%rho_target)) then
         call configure_rho_target(cfg, ocean_state%vcoord%rho_target, ocean_state%vcoord%nz_ml)
      end if
      ! ALE-regrid refinements (all default-off ⇒ bit-identical).
      ocean_state%vcoord%regrid_time_scale = cfg%regrid_time_scale
      ocean_state%vcoord%remap_vel_conserve_ke = cfg%remap_vel_conserve_ke
      ocean_state%vcoord%remap_boundary_extrap = cfg%remap_boundary_extrap
      ocean_state%vcoord%remap_nonuniform_weights = cfg%remap_nonuniform_weights
      ocean_state%vcoord%remap_check_preconditions = cfg%remap_check_preconditions
      ocean_state%vcoord%check_vanished_content = cfg%check_vanished_content
      if (compute_rank == 0) then
         if (ocean_state%vcoord%coord_type == VCOORD_EULERIAN_Z) then
            call logger%info("Vertical coord:  "//trim(cfg%vcoord_type)// &
                             " → EULERIAN_Z (no ALE remap)")
         else if (ocean_state%vcoord%coord_type == VCOORD_RHO) then
            ! RHO is state-dependent: the regrid silently falls back to a
            ! geometric grid if S/T or the EOS handle are unavailable at
            ! remap time.  Make activation + its requirement visible here
            ! (a hard configure-time S/T check is a deferred follow-up —
            ! S/T are always registered on the default ocean path).
            call logger%info("Vertical coord:  rho (isopycnal, ALE remap; "// &
                             "validation-grade — needs S+T + EOS; rho_ref_p="// &
                             to_string(cfg%rho_ref_pressure)//" Pa, light/dense="// &
                             to_string(cfg%rho_target_light)//"/"// &
                             to_string(cfg%rho_target_dense)//" kg/m^3)")
         else if (ocean_state%vcoord%coord_type == VCOORD_HYCOM) then
            ! HYCOM = the RHO density-space inversion + a z* nominal-floor
            ! sweep (fixed-resolution surface band, isopycnal interior).
            ! Reuses the same rho_target linspace + EOS as RHO; needs S+T.
            call logger%info("Vertical coord:  hycom (hybrid z*/isopycnal, ALE "// &
                             "remap; needs S+T + EOS; rho_ref_p="// &
                             to_string(cfg%rho_ref_pressure)//" Pa, light/dense="// &
                             to_string(cfg%rho_target_light)//"/"// &
                             to_string(cfg%rho_target_dense)//" kg/m^3)")
         else
            call logger%info("Vertical coord:  "//trim(cfg%vcoord_type)// &
                             " (code "//to_string(ocean_state%vcoord%coord_type)// &
                             ", ALE remap enabled)  h_surf_target="// &
                             to_string(cfg%zstar_h_surf_target)//" m")
         end if
         ! Echo the IC thickness seed only when it is NOT the default, so the
         ! banner stays unchanged for every existing config.  Worth surfacing:
         ! "uniform_z" changes the resting isopycnal geometry, not just a
         ! tolerance, and is otherwise invisible after t=0.
         if (trim(cfg%thickness_config) /= "sigma") then
            call logger%info("Initial thickness: "//trim(cfg%thickness_config)// &
                             " (MOM6 uniform-z interfaces over max_depth="// &
                             to_string(cfg%ocean%topo%max_depth)// &
                             " m, clipped to bathymetry; flat resting isopycnals)")
         end if
      end if

      ! Vertical mixing: PP81 interior + KPP boundary-layer overlay.
      ocean_state%vmix%use_closure = cfg%ocean%vmix%use_closure
      ocean_state%vmix%use_kpp = cfg%ocean%vmix%use_kpp
      ! Shared EOS handle: KPP B0 reads the surface α/β from here so
      ! it tracks the dyn-core EOS (was a private never-refreshed copy).
      ocean_state%vmix%eos = ocean_state%eos

      ! PP81 interior closure + KPP BL-depth constants (`&ocean_vmix_nml
      ! pp81_*`/`kpp_*`).  Each default equals the ocean_vmix_t field
      ! default, so an nml that never mentions these keys is bit-identical.
      ocean_state%vmix%pp81_nu0 = cfg%ocean%vmix%pp81_nu0
      ocean_state%vmix%pp81_nu_bg = cfg%ocean%vmix%pp81_nu_bg
      ocean_state%vmix%pp81_kappa_bg = cfg%ocean%vmix%pp81_kappa_bg
      ocean_state%vmix%pp81_alpha = cfg%ocean%vmix%pp81_alpha
      ocean_state%vmix%shear2_floor = cfg%ocean%vmix%shear2_floor
      ocean_state%vmix%ri_crit = cfg%ocean%vmix%kpp_ri_crit
      ocean_state%vmix%cs_nonlocal = cfg%ocean%vmix%kpp_cs_nonlocal
      ocean_state%vmix%c_vt2 = cfg%ocean%vmix%kpp_c_vt2
      ! (PR-21) KPP shortwave-in-BL method (&ocean_thermo_nml
      ! kpp_sw_method).  validate_config already rejected unknown strings.
      ocean_state%vmix%kpp_sw_method = &
         parse_kpp_sw_method(trim(cfg%ocean%thermo%kpp_sw_method))
      ! (E4) Source of the α/β pair behind the KPP `B_0` and the
      ! double-diffusion density ratio (`&ocean_vmix_nml
      ! buoyancy_coeffs`).  `nml_enum allowed=` already rejected an
      ! unknown spelling at parse time; the INVALID arm here is the
      ! belt-and-braces guard that keeps the schema list and
      ! `parse_buoyancy_coeffs` from drifting apart silently (a wrong α is
      ! a physics change with no symptom, so never a fallback).  Both this
      ! and the `p_top` gate below are plain host scalars latched BEFORE
      ! `ocean_state_enter_data`'s `copyin`, which is the same contract
      ! `rho0` and the `pp81_*` scalars beside them rely on.
      ocean_state%vmix%buoyancy_coeffs = &
         parse_buoyancy_coeffs(trim(cfg%ocean%vmix%buoyancy_coeffs))
      if (ocean_state%vmix%buoyancy_coeffs == BUOY_COEFFS_INVALID) then
         call fail("ocean_vmix_nml: buoyancy_coeffs='"// &
                   trim(cfg%ocean%vmix%buoyancy_coeffs)// &
                   "' is not a known source; use 'constant' or 'eos'", &
                   ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if
      ! Mirror of the E3 `&ocean_psurf_nml in_eos` seam gate, taken from
      ! the SAME knob EPBL's `in_eos` takes so the two boundary-layer
      ! schemes measure their in-situ pressures from the same origin.
      ! NOT redundant with `p_top` being zero — a cavity fills `p_top`
      ! with the ice load whether or not `in_eos` is set.
      ocean_state%vmix%p_top_in_eos = cfg%ocean%psurf%in_eos
      ! MANDATORY re-derive: kv_bg/kt_bg/ks_bg and the seeded kv/kt/ks/kd_bg
      ! arrays were set from the TYPE-DEFAULT pp81_* at `init` time (before
      ! this configure call runs); without this, `vmix_assemble`'s
      ! background floor would silently clamp against the old default even
      ! though the user just set a new pp81_nu_bg/pp81_kappa_bg — see
      ! `vmix_seed_backgrounds`'s docstring.  Runs before
      ! `ocean_state_enter_data` (driver), so no `!$acc update` is needed.
      call ocean_state%vmix%seed_backgrounds()

      ! EPBL — energetics-based PBL.  Replaces the KPP overlay when
      ! enabled (mutually exclusive surface schemes); the PP81
      ! interior + background continue underneath.
      if (present(ierr)) then
         call configure_ocean_epbl(cfg, ocean_state, grid, compute_rank, ierr=local_ierr)
         if (local_ierr /= 0) then
            ierr = local_ierr
            return
         end if
      else
         call configure_ocean_epbl(cfg, ocean_state, grid, compute_rank)
      end if

      ! Kappa-shear — JHL08 shear-driven INTERIOR closure.  Coexists
      ! with KPP and EPBL (no mutual exclusion); its kappa is added to
      ! the interior diffusivities.
      if (present(ierr)) then
         call configure_ocean_kappa_shear(cfg, ocean_state, grid, compute_rank, ierr=local_ierr)
         if (local_ierr /= 0) then
            ierr = local_ierr
            return
         end if
      else
         call configure_ocean_kappa_shear(cfg, ocean_state, grid, compute_rank)
      end if

      ! Tidal mixing — St-Laurent/Simmons internal-tide INTERIOR closure.
      ! Coexists with KPP/EPBL and PP81/background/kappa-shear (no mutual
      ! exclusion); its Kd is added to the interior diffusivities.
      if (present(ierr)) then
         call configure_ocean_tidal_mixing(cfg, ocean_state, compute_rank, ierr=local_ierr)
         if (local_ierr /= 0) then
            ierr = local_ierr
            return
         end if
      else
         call configure_ocean_tidal_mixing(cfg, ocean_state, compute_rank)
      end if

      ! Convective adjustment — Brunt-Vaisala-triggered INTERIOR closure
      ! CONTRIBUTOR (max() floor on kv/kt).  Coexists with KPP/EPBL (masks
      ! against whichever BL depth is live); no mutual exclusion.
      if (present(ierr)) then
         call configure_ocean_conv(cfg, ocean_state, compute_rank, ierr=local_ierr)
         if (local_ierr /= 0) then
            ierr = local_ierr
            return
         end if
      else
         call configure_ocean_conv(cfg, ocean_state, compute_rank)
      end if

      ! Double diffusion — salt fingering + diffusive convection, folded
      ! INTO the heat/salt split (asymmetric ks vs kt).  Reads the mirrored
      ! vmix%eos coefficients set above.
      if (present(ierr)) then
         call configure_ocean_ddiff(cfg, ocean_state, compute_rank, ierr=local_ierr)
         if (local_ierr /= 0) then
            ierr = local_ierr
            return
         end if
      else
         call configure_ocean_ddiff(cfg, ocean_state, compute_rank)
      end if

      ! Wave speed — B1 first-baroclinic cg1 + Rd diagnostic.
      if (present(ierr)) then
         call configure_ocean_wavespeed(cfg, ocean_state, grid, ierr=local_ierr)
         if (local_ierr /= 0) then
            ierr = local_ierr
            return
         end if
      else
         call configure_ocean_wavespeed(cfg, ocean_state, grid)
      end if
      ! Fox-Kemper MLE restratification (B5) — reads epbl%mld; configure
      ! after EPBL so the enable check sees the resolved EPBL state.
      if (present(ierr)) then
         call configure_ocean_foxkemper(cfg, ocean_state, compute_rank, ierr=local_ierr)
         if (local_ierr /= 0) then
            ierr = local_ierr
            return
         end if
      else
         call configure_ocean_foxkemper(cfg, ocean_state, compute_rank)
      end if

      ! VarMix (capability [4]) — build the STATIC f2_dx2/beta_dx2/L2 face
      ! terms once from the (filled) metrics + cell-centre Coriolis.  Runs
      ! after wavespeed (cg1) + metrics; the per-step Res_fn/SN/assembly
      ! fire at thermo cadence in the dyn step.  No-op when disabled.
      call configure_ocean_varmix(cfg, ocean_state, grid)

      ! MEKE (capability [5]) — fill the cell-centre |f| from the same
      ! Coriolis path VarMix/EPBL use, so the Rhines length is live when
      ! alpha_rhines > 0.  No-op when MEKE is disabled.  Scalar knobs are
      ! copied earlier by `ocean_state_copy_config`.
      call configure_ocean_meke(cfg, ocean_state, grid)

      if (compute_rank == 0) then
         if (ocean_state%epbl%enable) then
            call logger%info("Vertical mixing: PP81 interior + EPBL boundary layer")
         else if (ocean_state%vmix%use_closure .and. ocean_state%vmix%use_kpp) then
            call logger%info("Vertical mixing: PP81 interior + KPP overlay enabled")
         else if (ocean_state%vmix%use_closure) then
            call logger%info("Vertical mixing: PP81 interior only (KPP off)")
         else
            call logger%info("Vertical mixing: scalar K_v fallback (closure disabled)")
         end if
         if (cfg%ocean%vmix%pp81_nu0 /= 1.0e-2_wp .or. cfg%ocean%vmix%pp81_nu_bg /= 1.0e-4_wp .or. &
             cfg%ocean%vmix%pp81_kappa_bg /= 1.0e-5_wp .or. cfg%ocean%vmix%pp81_alpha /= 5.0_wp .or. &
             cfg%ocean%vmix%shear2_floor /= 1.0e-10_wp) then
            call logger%info("PP81 closure:     nu0="//to_string(cfg%ocean%vmix%pp81_nu0)// &
                             "  nu_bg="//to_string(cfg%ocean%vmix%pp81_nu_bg)// &
                             "  kappa_bg="//to_string(cfg%ocean%vmix%pp81_kappa_bg)// &
                             "  alpha="//to_string(cfg%ocean%vmix%pp81_alpha)// &
                             "  shear2_floor="//to_string(cfg%ocean%vmix%shear2_floor))
         end if
         if (cfg%ocean%vmix%kpp_ri_crit /= 0.3_wp .or. cfg%ocean%vmix%kpp_cs_nonlocal /= 6.3_wp .or. &
             cfg%ocean%vmix%kpp_c_vt2 /= 1.8_wp) then
            call logger%info("KPP constants:    ri_crit="//to_string(cfg%ocean%vmix%kpp_ri_crit)// &
                             "  cs_nonlocal="//to_string(cfg%ocean%vmix%kpp_cs_nonlocal)// &
                             "  c_vt2="//to_string(cfg%ocean%vmix%kpp_c_vt2))
         end if
      end if
      if (present(ierr)) ierr = OCEAN_STATUS_OK
   end subroutine configure_ocean_lateral

   subroutine configure_rho_target(cfg, rho_target, nz_ml)
      !! Populate the isopycnal `rho_target(0:nz_ml)` interface densities as
      !! a uniform light→dense linspace from `rho_target_light` /
      !! `rho_target_dense` (the MOM6 `ALE_COORDINATE_CONFIG=UNIFORM`
      !! analogue for the RHO mode; a stretched RFNC1-style profile is a
      !! P3 follow-up).  `rho_target(0)` is the surface (lightest)
      !! interface, `rho_target(nz_ml)` the bed (densest).  Only consulted
      !! for VCOORD_RHO; harmless otherwise.
      type(config_t), intent(in) :: cfg
      real(wp), intent(inout) :: rho_target(0:)
      integer, intent(in) :: nz_ml
      integer :: k
      real(wp) :: frac

      do k = 0, nz_ml
         if (nz_ml > 0) then
            frac = real(k, wp)/real(nz_ml, wp)
         else
            frac = 0.0_wp
         end if
         rho_target(k) = cfg%rho_target_light &
                         + (cfg%rho_target_dense - cfg%rho_target_light)*frac
      end do
   end subroutine configure_rho_target

   subroutine configure_ocean_epbl(cfg, ocean_state, grid, compute_rank, ierr)
      !! Copy the `&ocean_epbl_nml` knobs onto the EPBL slot — every
      !! field, end to end (don't repeat the KPP ri_crit/c_vt2
      !! dead-config gap).  Also fills `f_centre` from the same
      !! beta-plane parameters the Coriolis slot uses, copies the EOS
      !! hookup, validates the configuration, and resolves the
      !! EPBL-vs-KPP mutual exclusion.
      type(config_t), intent(in) :: cfg
      type(ocean_state_t), intent(inout) :: ocean_state
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: compute_rank
      integer, intent(out), optional :: ierr
         !! Non-zero on an EPBL configuration conflict when present; absent
         !! behaves as today (`error stop`).

      associate (epbl => ocean_state%epbl, ecfg => cfg%ocean%epbl)
         ! `enable` is flipped LAST, only once every validation below has
         ! passed (P0.1 review F9 — `configure_ocean_porous` is the
         ! pattern): setting it up front left a half-configured slot
         ! flagged enabled on any of this routine's many failure returns.
         if (.not. ecfg%enable) then
            epbl%enable = .false.
            if (present(ierr)) ierr = OCEAN_STATUS_OK
            return
         end if

         epbl%mstar_scheme = parse_epbl_mstar_scheme(ecfg%mstar_scheme)
         if (epbl%mstar_scheme < 0) then
            call fail("ocean_epbl_nml: unknown mstar_scheme '"// &
                      trim(ecfg%mstar_scheme)//"' (constant/om4/rh18)", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
         epbl%vstar_scheme = parse_epbl_vstar_scheme(ecfg%vel_scale_scheme)
         if (epbl%vstar_scheme < 0) then
            call fail("ocean_epbl_nml: unknown vel_scale_scheme '"// &
                      trim(ecfg%vel_scale_scheme)//"' (cube_root/rh18)", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
         epbl%combine_mode = parse_epbl_combine(ecfg%combine)
         if (epbl%combine_mode < 0) then
            call fail("ocean_epbl_nml: unknown combine '"// &
                      trim(ecfg%combine)//"' (add/max)", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if

         epbl%mstar_const = ecfg%mstar
         epbl%mstar_cap = ecfg%mstar_cap
         epbl%mstar_coef1 = ecfg%mstar_coef1
         epbl%c_ek = ecfg%c_ek
         epbl%mstar_conv_adj = ecfg%mstar_conv_adj
         epbl%rh18_cn1 = ecfg%rh18_cn1
         epbl%rh18_cn2 = ecfg%rh18_cn2
         epbl%rh18_cn3 = ecfg%rh18_cn3
         epbl%rh18_cs1 = ecfg%rh18_cs1
         epbl%rh18_cs2 = ecfg%rh18_cs2
         epbl%nstar = ecfg%nstar
         epbl%tke_decay = ecfg%tke_decay
         epbl%wstar_ustar_coef = ecfg%wstar_ustar_coef
         epbl%vstar_scale_fac = ecfg%vstar_scale_fac
         epbl%vstar_surf_fac = ecfg%vstar_surf_fac
         epbl%von_karman = ecfg%von_karman
         epbl%ekman_scale_coef = ecfg%ekman_scale_coef
         epbl%min_mix_len = ecfg%min_mix_len
         epbl%mixlen_exponent = ecfg%mixlen_exponent
         epbl%translay_scale = ecfg%translay_scale
         epbl%mld_iteration = ecfg%mld_iteration
         epbl%mld_tol = ecfg%mld_tol
         epbl%mld_max_its = ecfg%mld_max_its
         epbl%mld_bisection = ecfg%mld_bisection
         epbl%mld_use_prev_guess = ecfg%mld_use_prev_guess
         epbl%omega = ecfg%omega
         epbl%omega_frac = ecfg%omega_frac
         epbl%prandtl = ecfg%prandtl
         epbl%tke_diags = ecfg%tke_diags
         ! (PR-21) Penetrating-SW TKE ledger switch lives on
         ! &ocean_thermo_nml (shared with KPP), not &ocean_epbl_nml.
         epbl%epbl_sw_ctke = cfg%ocean%thermo%epbl_sw_ctke

         ! Langmuir (LF17 wind-only path).
         epbl%use_lt = ecfg%use_lt
         if (ecfg%use_lt) then
            epbl%lt_scheme = parse_epbl_lt_scheme(ecfg%lt_scheme)
            if (epbl%lt_scheme < 0) then
               call fail("ocean_epbl_nml: unknown lt_scheme '"// &
                         trim(ecfg%lt_scheme)//"' (rescale/additive)", ierr, OCEAN_STATUS_ERR_SETUP)
               return
            end if
            epbl%lt_enhance_coef = ecfg%lt_enhance_coef
            epbl%lt_enhance_exp = ecfg%lt_enhance_exp
            epbl%lt_max_enhance = ecfg%lt_max_enhance
            epbl%la_frac_hbl = ecfg%la_frac_hbl
            epbl%lt_lac1 = ecfg%lt_lac1
            epbl%lt_lac2 = ecfg%lt_lac2
            epbl%lt_lac3 = ecfg%lt_lac3
            epbl%lt_lac4 = ecfg%lt_lac4
            epbl%lt_lac5 = ecfg%lt_lac5
            if (compute_rank == 0) then
               call logger%info("EPBL Langmuir:    LF17 wind-only, scheme="// &
                                trim(ecfg%lt_scheme)//"  coef="// &
                                to_string(ecfg%lt_enhance_coef)//"  exp="// &
                                to_string(ecfg%lt_enhance_exp))
            end if
         end if

         ! EOS hookup: the EPBL energy weights use the SAME EOS handle
         ! the dyn-core runs (shared flat-POD copy — one source of
         ! truth, no scalar copies that can drift).
         epbl%eos = ocean_state%eos
         epbl%rho0 = ocean_state%eos%rho0

         ! (E3) Top-of-column pressure in the IN-SITU EOS arguments AND
         ! in the PE weight.  `&ocean_psurf_nml in_eos` is the single
         ! gate for the whole `ms%p_top` seam; assigned HERE, at
         ! configure, so it is latched before `ocean_state_enter_data`
         ! and the kernel reads it by value (host scalar -- no
         ! `!$acc update device` owed).  Off (default) ⇒ the stack
         ! starts at 0 Pa, bit-identically.  NOTE this is NOT redundant
         ! with `p_top` being zero: a cavity fills `p_top` with the ice
         ! load whether or not `in_eos` is set.
         epbl%in_eos = cfg%ocean%psurf%in_eos

         ! |f| at cell centres, same scheme the Coriolis slot uses (D7).
         call fill_coriolis_centre(cfg, ocean_state%metrics, grid, epbl%f_centre)

         ! Validation.
         if (.not. cfg%ocean%vmix%use_closure) then
            call fail("ocean_epbl_nml: enable=.true. requires "// &
                      "ocean_vmix_nml use_closure=.true. (EPBL folds "// &
                      "into the interior closure's kv/kt)", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
         if (.not. cfg%ocean%thermo%enable_thermodynamics) then
            call fail("ocean_epbl_nml: enable=.true. requires "// &
                      "active thermodynamics (EPBL reads T, S)", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
         if (epbl%mld_iteration .and. &
             (epbl%translay_scale < 0.0_wp .or. epbl%translay_scale >= 1.0_wp)) then
            call fail("ocean_epbl_nml: translay_scale must be in "// &
                      "[0,1) when mld_iteration is on", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
         ! KPP mutual exclusion: EPBL replaces the overlay.  use_kpp
         ! defaults on, so override with a log line rather than abort.
         if (ocean_state%vmix%use_kpp) then
            ocean_state%vmix%use_kpp = .false.
            if (compute_rank == 0) then
               call logger%info("EPBL enabled: KPP overlay disabled "// &
                                "(mutually exclusive surface schemes)")
            end if
         end if

         if (compute_rank == 0) then
            call logger%info("EPBL:             mstar="//trim(ecfg%mstar_scheme)// &
                             "  nstar="//to_string(ecfg%nstar)// &
                             "  tke_decay="//to_string(ecfg%tke_decay)// &
                             "  vstar="//trim(ecfg%vel_scale_scheme)// &
                             "  combine="//trim(ecfg%combine))
         end if
         epbl%enable = .true.
      end associate
      if (present(ierr)) ierr = OCEAN_STATUS_OK
   end subroutine configure_ocean_epbl

   subroutine configure_ocean_foxkemper(cfg, ocean_state, compute_rank, ierr)
      !! Copy the `&ocean_foxkemper_nml` knobs onto the MLE slot (B5).
      !! Validates: B5 reads `epbl%mld`, so EPBL must be enabled; and the
      !! `resolution_taper` hook is a hard error until B2 lands.
      type(config_t), intent(in) :: cfg
      type(ocean_state_t), intent(inout) :: ocean_state
      integer, intent(in) :: compute_rank
      integer, intent(out), optional :: ierr
         !! Non-zero on a Fox-Kemper MLE configuration conflict when
         !! present; absent behaves as today (`error stop`).

      associate (mle => ocean_state%mle, fcfg => cfg%ocean%foxkemper)
         mle%enable = fcfg%enable
         if (.not. fcfg%enable) then
            if (present(ierr)) ierr = OCEAN_STATUS_OK
            return
         end if

         mle%ce = fcfg%ce
         mle%f_floor = fcfg%f_floor
         mle%mld_decay_time = fcfg%mld_decay_time
         mle%tail_dh = fcfg%tail_dh
         mle%use_mom_mixrate = fcfg%use_mom_mixrate
         mle%resolution_taper = fcfg%resolution_taper
         mle%use_bodner = fcfg%use_bodner
         mle%cr = fcfg%cr
         mle%bodner_mstar = fcfg%bodner_mstar
         mle%bodner_nstar = fcfg%bodner_nstar
         mle%min_wstar2 = fcfg%min_wstar2

         ! B5 takes the mixed-layer depth from epbl%mld — EPBL is the
         ! enabler and must be on.
         if (.not. ocean_state%epbl%enable) then
            call fail("ocean_foxkemper_nml: enable=.true. requires "// &
                      "ocean_epbl_nml enable=.true. (B5 reads epbl%mld)", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
         ! resolution_taper is a B2 hook; the res_fn is not available yet.
         if (mle%resolution_taper) then
            call fail("ocean_foxkemper_nml: resolution_taper=.true. "// &
                      "requires the B2 resolution function (not yet "// &
                      "available) — leave it .false. until B2 lands", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
         if (mle%tail_dh /= 0.0_wp) then
            call fail("ocean_foxkemper_nml: tail_dh /= 0 (mu cubic "// &
                      "tail) is deferred; use the exact mu (tail_dh=0)", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if

         if (compute_rank == 0) then
            if (fcfg%use_bodner) then
               call logger%info("Fox-Kemper MLE:   restratification enabled  "// &
                                "Bodner-2023 frontogenesis arrest  Cr="//to_string(fcfg%cr))
            else
               call logger%info("Fox-Kemper MLE:   restratification enabled  Ce="// &
                                to_string(fcfg%ce)//"  timescale="// &
                                merge("FK11-mixrate", "bare Ce/|f| ", fcfg%use_mom_mixrate))
            end if
         end if
      end associate
      if (present(ierr)) ierr = OCEAN_STATUS_OK
   end subroutine configure_ocean_foxkemper

   subroutine configure_ocean_kappa_shear(cfg, ocean_state, grid, compute_rank, ierr)
      !! Copy the `&ocean_kappa_shear_nml` knobs onto the kappa-shear slot
      !! — every field, end to end (don't repeat the KPP dead-config
      !! gap).  Fills `f_centre` from the same beta-plane parameters the
      !! Coriolis slot uses, copies the EOS hookup, and validates.  No
      !! mutual exclusion: kappa-shear is an interior closure that
      !! coexists with KPP / EPBL and PP81/background.
      type(config_t), intent(in) :: cfg
      type(ocean_state_t), intent(inout) :: ocean_state
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: compute_rank
      integer, intent(out), optional :: ierr
         !! Non-zero on a kappa-shear configuration conflict when present;
         !! absent behaves as today (`error stop`).

      associate (ks => ocean_state%kshear, kcfg => cfg%ocean%kshear)
         ! `enable` is flipped LAST, only once every validation below has
         ! passed (P0.1 review F9 — `configure_ocean_porous` is the
         ! pattern): setting it up front left a half-configured slot
         ! flagged enabled on any of this routine's many failure returns.
         if (.not. kcfg%enable) then
            ks%enable = .false.
            if (present(ierr)) ierr = OCEAN_STATUS_OK
            return
         end if

         ! Fail loud at configure (wavespeed precedent): both kappa-shear
         ! kernels gather into fixed-size NZ_STACK_MAX column arrays, so
         ! nz > NZ_STACK_MAX would silently overrun per-thread stack.
         if (ocean_state%multilayer%nz_ml > NZ_STACK_MAX) then
            call fail("ocean_kappa_shear_nml: nz_layers exceeds "// &
                      "NZ_STACK_MAX (raise NZ_STACK_MAX in rdb_constants)", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if

         ! Every JHL08 knob (no dead-config gaps).
         ks%ri_crit = kcfg%ri_crit
         ks%shearmix_rate = kcfg%shearmix_rate
         ks%fri_curvature = kcfg%fri_curvature
         ks%c_n = kcfg%c_n
         ks%c_s = kcfg%c_s
         ks%lambda = kcfg%lambda
         ks%lz_rescale = kcfg%lz_rescale
         ks%kappa_0 = kcfg%kappa_0
         ks%kappa_seed = kcfg%kappa_seed
         ks%kappa_trunc = kcfg%kappa_trunc
         ks%tke_bg = kcfg%tke_bg
         ks%tol_err = kcfg%tol_err
         ks%max_inner_it = kcfg%max_inner_it
         ks%max_substep_it = kcfg%max_substep_it
         ks%src_max_chg = kcfg%src_max_chg
         ks%prandtl_turb = kcfg%prandtl_turb
         ks%vel_underflow = kcfg%vel_underflow
         ks%massless_merge = kcfg%massless_merge
         ks%vertex_geometric_mean = kcfg%vertex_geometric_mean
         ks%vertex_geomean_kdmin = kcfg%vertex_geomean_kdmin

         ! EOS hookup: the buoyancy derivatives use the SAME EOS handle
         ! the dyn-core runs (shared flat-POD copy — one source of
         ! truth, no scalar copies that can drift).
         ks%eos = ocean_state%eos
         ks%rho0 = ocean_state%eos%rho0

         ! |f| at cell centres, same scheme the Coriolis slot uses (D7).
         call fill_coriolis_centre(cfg, ocean_state%metrics, grid, ks%f_centre)

         ! Vertex form (MOM6 VERTEX_SHEAR): allocate the corner fields
         ! (deliberately not in `init` — the corner carrier is nz+1 full
         ! planes, paid only when selected) and fill the SIGNED corner
         ! f via the same generator-driven scheme (the kernel squares
         ! it).  Runs before `ocean_state_enter_data`, so the new
         ! allocatables map with the slot.
         if (kcfg%at_vertex) then
            ! nghost >= 2 fail-loud (pv_adv_required_nghost precedent):
            ! corners are solved on [2,nx]x[2,ny] only, so the outermost
            ! array ring averages with never-solved zero ring corners.
            ! With nghost >= 2 that ring is entirely ghost cells and the
            ! owned-cell answer is exact; with nghost < 2 the depressed
            ! first-ring kt could leak into owned cells via smoothing.
            if (grid%nghost < 2) then
               call fail("ocean_kappa_shear_nml: at_vertex requires "// &
                         "nghost >= 2 (corner solve covers interior "// &
                         "corners only; got nghost="//to_string(grid%nghost)//")", ierr, OCEAN_STATUS_ERR_SETUP)
               return
            end if
            call ks%init_vertex(grid, ocean_state%multilayer%nz_ml)
            call fill_coriolis_corner(cfg, ocean_state%metrics, grid, ks%f_corner)
         end if

         ! Validation.
         if (.not. cfg%ocean%vmix%use_closure) then
            call fail("ocean_kappa_shear_nml: enable=.true. requires "// &
                      "ocean_vmix_nml use_closure=.true. (kappa-shear "// &
                      "folds into the interior closure's kv/kt)", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
         if (.not. cfg%ocean%thermo%enable_thermodynamics) then
            call fail("ocean_kappa_shear_nml: enable=.true. requires "// &
                      "active thermodynamics (kappa-shear reads T, S)", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
         if (ks%ri_crit <= 0.0_wp) then
            call fail("ocean_kappa_shear_nml: ri_crit must be > 0", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
         if (ks%max_inner_it < 1) then
            call fail("ocean_kappa_shear_nml: max_inner_it must be >= 1", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
         if (ks%max_substep_it < 1) then
            call fail("ocean_kappa_shear_nml: max_substep_it must be >= 1", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
         if (ks%kappa_0 <= 0.0_wp) then
            call fail("ocean_kappa_shear_nml: kappa_0 must be > 0", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
         if (ks%tol_err <= 0.0_wp) then
            call fail("ocean_kappa_shear_nml: tol_err must be > 0", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
         if (ks%prandtl_turb <= 0.0_wp) then
            call fail("ocean_kappa_shear_nml: prandtl_turb must be > 0", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
         if (ks%vertex_geomean_kdmin < 0.0_wp) then
            call fail("ocean_kappa_shear_nml: vertex_geomean_kdmin must be >= 0", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
         if (compute_rank == 0) then
            ! Advisory warnings (MOM6-parity semantics, not errors).
            if (kcfg%vertex_geomean_kdmin > 0.0_wp .and. &
                .not. kcfg%vertex_geometric_mean) then
               call logger%warning("ocean_kappa_shear_nml: vertex_geomean_kdmin "// &
                                   "is inert without vertex_geometric_mean")
            end if
            if (kcfg%at_vertex .and. kcfg%vertex_geometric_mean .and. &
                kcfg%vertex_geomean_kdmin == 0.0_wp) then
               call logger%warning("ocean_kappa_shear_nml: geometric mean with "// &
                                   "kdmin=0 hard-zeros Kd wherever ANY corner is 0 "// &
                                   "(every shear-zone edge); production configs "// &
                                   "use vertex_geomean_kdmin=1e-9")
            end if
         end if

         if (compute_rank == 0) then
            call logger%info("Kappa-shear:      interior closure enabled "// &
                             "(ri_crit="//to_string(kcfg%ri_crit)// &
                             "  prandtl_turb="//to_string(kcfg%prandtl_turb)// &
                             trim(merge("  at_vertex", "           ", kcfg%at_vertex))//")")
         end if
         ks%enable = .true.
      end associate
      if (present(ierr)) ierr = OCEAN_STATUS_OK
   end subroutine configure_ocean_kappa_shear

   subroutine configure_ocean_tidal_mixing(cfg, ocean_state, compute_rank, ierr)
      !! Copy the `&ocean_tidal_mixing_nml` knobs onto the tidal-mixing
      !! slot — every field, end to end (no dead-config gaps).  Seeds the
      !! prescribed bottom energy field `e_in` (v1 uniform path), copies
      !! the shared EOS hookup for the N^2 buoyancy derivatives, and
      !! validates.  No mutual exclusion: tidal mixing is an interior
      !! closure that coexists with KPP / EPBL / PP81 / kappa-shear.
      type(config_t), intent(in) :: cfg
      type(ocean_state_t), intent(inout) :: ocean_state
      integer, intent(in) :: compute_rank
      integer, intent(out), optional :: ierr
         !! Non-zero on a tidal-mixing configuration conflict when present;
         !! absent behaves as today (`error stop`).

      associate (tm => ocean_state%vmix_tidal, tcfg => cfg%ocean%tidal_mixing)
         ! `enable` is flipped LAST, only once every validation below has
         ! passed (P0.1 review F9 — `configure_ocean_porous` is the
         ! pattern): setting it up front left a half-configured slot
         ! flagged enabled on any of this routine's failure returns.
         if (.not. tcfg%enable) then
            tm%enable = .false.
            if (present(ierr)) ierr = OCEAN_STATUS_OK
            return
         end if

         ! Every St-Laurent/Simmons knob (no dead-config gaps).
         tm%gamma = tcfg%gamma
         tm%mu = tcfg%mu
         tm%zeta = tcfg%zeta
         tm%kd_max = tcfg%kd_max
         tm%prandtl_tidal = tcfg%prandtl_tidal
         tm%min_zbot = tcfg%min_zbot
         tm%e_uniform = tcfg%e_uniform
         tm%e_compute = tcfg%e_compute
         tm%kappa_itides = tcfg%kappa_itides
         tm%kappa_h2 = tcfg%kappa_h2
         tm%utide = tcfg%utide
         tm%h2_rough = tcfg%h2_rough
         tm%frac_rough = tcfg%frac_rough
         tm%e_max = tcfg%e_max

         ! EOS hookup: the buoyancy derivatives use the SAME EOS handle
         ! the dyn-core runs (shared flat-POD copy — one source of truth).
         tm%eos = ocean_state%eos
         tm%rho0 = ocean_state%eos%rho0

         ! Seed the prescribed bottom energy field (v1 uniform path).
         call tm%set_e_uniform(tm%e_uniform)

         ! Validation.
         if (.not. cfg%ocean%vmix%use_closure) then
            call fail("ocean_tidal_mixing_nml: enable=.true. requires "// &
                      "ocean_vmix_nml use_closure=.true. (tidal mixing "// &
                      "folds into the interior closure's kv/kt)", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
         if (.not. cfg%ocean%thermo%enable_thermodynamics) then
            call fail("ocean_tidal_mixing_nml: enable=.true. requires "// &
                      "active thermodynamics (tidal mixing reads T, S for N^2)", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
         if (tm%zeta <= 0.0_wp) then
            call fail("ocean_tidal_mixing_nml: zeta must be > 0", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
         if (tm%prandtl_tidal <= 0.0_wp) then
            call fail("ocean_tidal_mixing_nml: prandtl_tidal must be > 0", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if

         if (compute_rank == 0) then
            call logger%info("Tidal mixing:     interior closure enabled "// &
                             "(zeta="//to_string(tcfg%zeta)//" m  "// &
                             "E="//to_string(tcfg%e_uniform)//" W/m^2)")
         end if
         tm%enable = .true.
      end associate
      if (present(ierr)) ierr = OCEAN_STATUS_OK
   end subroutine configure_ocean_tidal_mixing

   subroutine configure_ocean_conv(cfg, ocean_state, compute_rank, ierr)
      !! Copy the `&ocean_conv_nml` knobs onto `ocean_state%vmix` -- every
      !! field, no dead-config gaps.  Convective adjustment adds no new
      !! slot / allocatable (it lives on the already-unconditionally-
      !! allocated `vmix`), so unlike `configure_ocean_tidal_mixing` there
      !! is no separate `enable` latch to set on a distinct sub-object;
      !! `vmix%conv_enable` IS the latch.  No mutual exclusion with
      !! KPP / EPBL: convection masks against whichever BL depth is live
      !! this stage (`rdb_ocean_dyn.F90 vmix_apply_in_stage`).
      type(config_t), intent(in) :: cfg
      type(ocean_state_t), intent(inout) :: ocean_state
      integer, intent(in) :: compute_rank
      integer, intent(out), optional :: ierr
         !! Non-zero on a convective-adjustment configuration conflict when
         !! present; absent behaves as today (`error stop`).

      associate (vmix => ocean_state%vmix, ccfg => cfg%ocean%conv)
         ! `conv_enable` is flipped LAST, only once every validation below
         ! has passed (P0.1 review F9 — `configure_ocean_porous` is the
         ! pattern): setting it up front left a half-configured slot
         ! flagged enabled on any of this routine's failure returns.
         if (.not. ccfg%enable) then
            vmix%conv_enable = .false.
            if (present(ierr)) ierr = OCEAN_STATUS_OK
            return
         end if

         vmix%conv_kd = ccfg%kd_conv
         vmix%conv_prandtl = ccfg%prandtl_conv
         vmix%conv_n2_thresh = ccfg%n2_thresh

         ! Validation.
         if (.not. cfg%ocean%vmix%use_closure) then
            call fail("ocean_conv_nml: enable=.true. requires "// &
                      "ocean_vmix_nml use_closure=.true. (convection "// &
                      "folds into the interior closure's kv/kt)", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
         if (.not. cfg%ocean%thermo%enable_thermodynamics) then
            call fail("ocean_conv_nml: enable=.true. requires "// &
                      "active thermodynamics (convection reads rho_layer, "// &
                      "which is only refreshed by the EOS)", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
         if (ccfg%kd_conv < 0.0_wp) then
            call fail("ocean_conv_nml: kd_conv must be >= 0", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
         if (ccfg%prandtl_conv <= 0.0_wp) then
            call fail("ocean_conv_nml: prandtl_conv must be > 0", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if

         if (compute_rank == 0) then
            call logger%info("Convection:       interior closure enabled "// &
                             "(kd_conv="//to_string(ccfg%kd_conv)//" m^2/s  "// &
                             "n2_thresh="//to_string(ccfg%n2_thresh)//" s^-2)")
         end if
         vmix%conv_enable = .true.
      end associate
      if (present(ierr)) ierr = OCEAN_STATUS_OK
   end subroutine configure_ocean_conv

   subroutine configure_ocean_ddiff(cfg, ocean_state, compute_rank, ierr)
      !! Copy the `&ocean_ddiff_nml` knobs onto `ocean_state%vmix` -- every
      !! field, no dead-config gaps.  Like convection, double diffusion adds
      !! no new slot (it rides the unconditionally-allocated `vmix` and the
      !! already-mirrored `vmix%eos`); `vmix%ddiff_enable` IS the latch.
      !! Folded into `vmix_split_kd_heat_salt`, so it needs the interior
      !! closure pipeline (`use_closure`) and thermodynamics (it reads T/S
      !! and the EOS alpha/beta).  `enable=.false.` => bit-identical.
      type(config_t), intent(in) :: cfg
      type(ocean_state_t), intent(inout) :: ocean_state
      integer, intent(in) :: compute_rank
      integer, intent(out), optional :: ierr
         !! Non-zero on a double-diffusion configuration conflict when
         !! present; absent behaves as today (`error stop`).

      associate (vmix => ocean_state%vmix, dcfg => cfg%ocean%ddiff)
         ! `ddiff_enable` is flipped LAST, only once every validation
         ! below has passed (P0.1 review F9 — `configure_ocean_porous` is
         ! the pattern): setting it up front left a half-configured slot
         ! flagged enabled on any of this routine's failure returns.
         if (.not. dcfg%enable) then
            vmix%ddiff_enable = .false.
            if (present(ierr)) ierr = OCEAN_STATUS_OK
            return
         end if

         vmix%ddiff_strat_param_max = dcfg%strat_param_max
         vmix%ddiff_kappa_s = dcfg%kappa_ddiff_s
         vmix%ddiff_exp1 = dcfg%ddiff_exp1
         vmix%ddiff_exp2 = dcfg%ddiff_exp2
         vmix%ddiff_param1 = dcfg%param1
         vmix%ddiff_param2 = dcfg%param2
         vmix%ddiff_param3 = dcfg%param3
         vmix%ddiff_mol_diff = dcfg%mol_diff
         vmix%ddiff_use_k90 = dcfg%use_k90

         ! Validation.
         if (.not. cfg%ocean%vmix%use_closure) then
            call fail("ocean_ddiff_nml: enable=.true. requires "// &
                      "ocean_vmix_nml use_closure=.true. (double diffusion "// &
                      "folds into the interior closure's heat/salt split)", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
         if (.not. cfg%ocean%thermo%enable_thermodynamics) then
            call fail("ocean_ddiff_nml: enable=.true. requires "// &
                      "active thermodynamics (double diffusion reads T/S "// &
                      "and the EOS alpha/beta)", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
         if (dcfg%strat_param_max <= 1.0_wp) then
            call fail("ocean_ddiff_nml: strat_param_max must be > 1 "// &
                      "(the fingering form divides by strat_param_max - 1)", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
         if (dcfg%kappa_ddiff_s < 0.0_wp .or. dcfg%mol_diff < 0.0_wp) then
            call fail("ocean_ddiff_nml: kappa_ddiff_s and mol_diff must be >= 0", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if

         if (compute_rank == 0) then
            call logger%info("Double diffusion: enabled (Rrho_max="// &
                             to_string(dcfg%strat_param_max)//"  K_f="// &
                             to_string(dcfg%kappa_ddiff_s)//" m^2/s  "// &
                             "conv="//trim(merge("K90 ", "MC76", dcfg%use_k90))//")")
         end if
         vmix%ddiff_enable = .true.
      end associate
      if (present(ierr)) ierr = OCEAN_STATUS_OK
   end subroutine configure_ocean_ddiff

   subroutine configure_ocean_tides(cfg, ocean_state, grid, compute_rank, ierr)
      !! Configure the equilibrium body-force tide slot (C1): parse the
      !! constituent list + reference dates, fill the astronomy catalog
      !! (phase0, nodal f/u), and build the (nx,ny,3) spatial-structure
      !! arrays from `metrics%geolatT/geolonT`.  Host-side setup — runs
      !! AFTER `configure_ocean_metrics` (lat/lon must be filled) and
      !! BEFORE `ocean_state_enter_data`.  `enable=.false.` => no-op,
      !! bit-identical.
      type(config_t), intent(in) :: cfg
      type(ocean_state_t), intent(inout) :: ocean_state
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: compute_rank
      integer, intent(out), optional :: ierr
         !! Non-zero on a tides configuration conflict when present;
         !! absent behaves as today (`error stop`).
      integer :: cat_idx(TIDES_CATALOG_SIZE)
      integer :: nconst, y, m, d, nx, ny, local_ierr
      real(wp) :: dref, dnodal
      logical :: ok
      character(len=16) :: nodal_str

      associate (td => ocean_state%tides, tcfg => cfg%ocean%tides)
         td%enable = tcfg%enable
         if (.not. tcfg%enable) then
            if (present(ierr)) ierr = OCEAN_STATUS_OK
            return
         end if

         td%use_sal = tcfg%use_sal
         td%beta_sal = tcfg%beta_sal
         td%t_epoch = 0.0_wp

         if (td%use_sal .and. td%beta_sal <= 0.0_wp .and. compute_rank == 0) then
            call logger%warning("ocean_tides_nml: use_sal=.true. but "// &
                                "beta_sal<=0 => scalar SAL is inert "// &
                                "(typical beta ~0.085-0.12)")
         end if

         ! Parse the active constituent list -> catalog indices.
         if (present(ierr)) then
            call parse_constituent_list(tcfg%constituents, cat_idx, nconst, ierr=local_ierr)
            if (local_ierr /= 0) then
               ierr = local_ierr
               return
            end if
         else
            call parse_constituent_list(tcfg%constituents, cat_idx, nconst)
         end if
         if (nconst < 1) then
            call fail("ocean_tides_nml: constituents='"// &
                      trim(tcfg%constituents)//"' has no valid entries", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if

         ! Reference date -> day number since 1900-01-01.
         call parse_date_string(tcfg%ref_date, y, m, d, ok)
         if (.not. ok) then
            call fail("ocean_tides_nml: bad ref_date '"// &
                      trim(tcfg%ref_date)//"' (want YYYY-MM-DD)", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
         dref = days_since_1900(y, m, d)

         nodal_str = tcfg%nodal_ref_date
         if (len_trim(nodal_str) == 0) nodal_str = tcfg%ref_date
         call parse_date_string(nodal_str, y, m, d, ok)
         if (.not. ok) then
            call fail("ocean_tides_nml: bad nodal_ref_date '"// &
                      trim(nodal_str)//"' (want YYYY-MM-DD)", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
         dnodal = days_since_1900(y, m, d)

         nx = grid%nx_total
         ny = grid%ny_total
         call tides_configure_astronomy(td, cat_idx(1:nconst), nconst, &
                                        dref, dnodal, tcfg%add_nodal, nx, ny)
         call tides_build_struct(td, ocean_state%metrics%geolatT, &
                                 ocean_state%metrics%geolonT, nx, ny)

         if (compute_rank == 0) then
            call logger%info("Equilibrium tide: C1 body force enabled ("// &
                             to_string(nconst)//" constituents, ref="// &
                             trim(tcfg%ref_date)//")")
         end if
      end associate
      if (present(ierr)) ierr = OCEAN_STATUS_OK
   end subroutine configure_ocean_tides

   subroutine configure_ocean_p_surf(cfg, ocean_state, grid, compute_rank)
      !! Configure the atmospheric surface-pressure loading slot (PR-17):
      !! copy `enable`, take ρ₀ from `ocean_state%eos%rho0` (the single ρ₀
      !! of record — NOT a namelist knob), and allocate the seam fields.
      !! Host-side setup — runs BEFORE `ocean_state_enter_data`.
      !! `enable=.false.` => no-op, bit-identical.  The uniform-`p_surf`
      !! inert warning is emitted in `validate_config`; here we only log the
      !! enable on rank 0.
      type(config_t), intent(in) :: cfg
      type(ocean_state_t), intent(inout) :: ocean_state
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: compute_rank
      integer :: nx, ny

      associate (ps => ocean_state%p_surf, pcfg => cfg%ocean%psurf)
         ps%enable = pcfg%enable
         if (.not. pcfg%enable) return

         ! Single ρ₀ of record: the model's own reference density, so the
         ! inverse barometer scales with whatever `&ocean_ic_nml rho_0`
         ! set.  `eos%rho0` is assigned in `ocean_state_init_from_config`
         ! (runs earlier), so it is finalised here.
         ps%rho0 = ocean_state%eos%rho0

         nx = grid%nx_total
         ny = grid%ny_total
         call p_surf_configure(ps, nx, ny)

         if (compute_rank == 0) then
            call logger%info("Surface pressure: inverse-barometer PGF enabled "// &
                             "(eta_ib = -p_surf/(rho0*g_bt), rho0 = "// &
                             to_string(ps%rho0)//" kg/m^3)")
         end if
      end associate
   end subroutine configure_ocean_p_surf

   subroutine configure_ocean_wave_drag(cfg, ocean_state, grid, compute_rank, ierr)
      !! Configure the barotropic linear (Rayleigh) wave-drag piston-velocity
      !! maps (Egbert & Ray 2001; Jayne & St Laurent 2001) — the bulk energy
      !! sink for the barotropic tide, MOM6 `BT_LINEAR_WAVE_DRAG`. Builds a
      !! host-only h-point `r_h(nx,ny)` map (uniform scalar or a
      !! resolved-bathymetry roughness proxy), scales it, averages h->face
      !! into `bt_work%lwd_drag_u/v`, and leaves the arrays unallocated when
      !! the knob is off (bit-identical). MUST run AFTER bathymetry is set
      !! (`ocean_state%barotropic%b`) and land masking
      !! (`configure_ocean_land_mask`) and BEFORE `ocean_state_enter_data` —
      !! the `!$acc enter data copyin` in `barotropic_workstate_enter_data`
      !! carries these host-filled values to the device (CLAUDE.md gotcha
      !! (2): arrays mapped `create` do not carry pre-map host values, so
      !! this ordering is load-bearing).
      type(config_t), intent(in) :: cfg
      type(ocean_state_t), intent(inout) :: ocean_state
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: compute_rank
      integer, intent(out), optional :: ierr
         !! Non-zero on an unreachable/unimplemented `wave_drag_form` when
         !! present; absent behaves as today (`error stop`).
      real(wp), allocatable :: r_h(:, :)
      integer :: i, j, nx, ny

      associate (bw => ocean_state%dyn%bt_work, bcfg => cfg%ocean%bt)
         bw%lwd_enable = bcfg%wave_drag
         if (.not. bcfg%wave_drag) then
            if (present(ierr)) ierr = OCEAN_STATUS_OK
            return
         end if

         nx = grid%nx_total
         ny = grid%ny_total
         allocate (bw%lwd_drag_u(nx + 1, ny), source=0.0_wp)
         allocate (bw%lwd_drag_v(nx, ny + 1), source=0.0_wp)
         allocate (r_h(nx, ny), source=0.0_wp)

         select case (trim(bcfg%wave_drag_form))
         case ("uniform")
            r_h = bcfg%wave_drag_r_uniform
         case ("roughness_proxy")
            call wave_drag_roughness_proxy(ocean_state%barotropic%b, &
                                           ocean_state%metrics%wet_T, &
                                           nx, ny, grid%nghost, bcfg%wave_drag_kappa, &
                                           bcfg%wave_drag_n_bot, bcfg%wave_drag_h2_max, r_h)
         case default
            ! Unreachable: validate_config fails loud on any other tag
            ! (including "file", registered but not yet implemented).
            call fail("configure_ocean_wave_drag: unreachable wave_drag_form='"// &
                      trim(bcfg%wave_drag_form)//"'", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end select

         r_h = bcfg%wave_drag_scale*r_h

         ! h -> face average (as MOM6 forms its wave-drag face map), built ONCE
         ! here; the array-edge faces (i=1/nx+1, j=1/ny+1) stay at their
         ! source=0.0 init value.
         do i = 2, nx
            do j = 1, ny
               bw%lwd_drag_u(i, j) = 0.5_wp*(r_h(i - 1, j) + r_h(i, j))
            end do
         end do
         do j = 2, ny
            do i = 1, nx
               bw%lwd_drag_v(i, j) = 0.5_wp*(r_h(i, j - 1) + r_h(i, j))
            end do
         end do

         if (compute_rank == 0) then
            call logger%info("Barotropic wave drag: ON  form="//trim(bcfg%wave_drag_form)// &
                             "  scale="//to_string(bcfg%wave_drag_scale)// &
                             "  max(r_H)="//to_string(maxval(bw%lwd_drag_u))// &
                             " m/s  mean(r_H)="// &
                             to_string(sum(bw%lwd_drag_u)/real(size(bw%lwd_drag_u), wp))//" m/s")
         end if
      end associate
      if (present(ierr)) ierr = OCEAN_STATUS_OK
   end subroutine configure_ocean_wave_drag

   subroutine configure_ocean_z_fixed_profile(cfg, ocean_state, compute_rank, log_it)
      !! Resolve the `VCOORD_Z_FIXED` nominal layering onto the vcoord slot:
      !! `z_fixed_h_ref = &ocean_topo_nml max_depth` (the uniform
      !! `max_depth/nz` spacing — the default, byte-identical), or, under
      !! `&vcoord_nml z_fixed_profile = "list" | "tanh"`, the stretched
      !! per-layer tables `z_fixed_zi` / `z_fixed_dz` built by
      !! `rdb_vcoord :: z_fixed_nominal_dz` (`z_fixed_h_ref` then becomes
      !! the profile's total depth).
      !!
      !! Idempotent (it rebuilds from `cfg` each call).  Called twice: by
      !! `engine_setup` BEFORE the IC seed — the cavity `z_fixed` seed lays
      !! `h_layer` from the same target builder and must see the same
      !! profile — and from `configure_ocean_lateral`, which has always
      !! owned `z_fixed_h_ref`.  Both precede `ocean_state_enter_data`, so
      !! the copyin captures the tables.  `validate_config` has already
      !! refused a profile that does not build, or one on another family.
      type(config_t), intent(in) :: cfg
      type(ocean_state_t), intent(inout) :: ocean_state
      integer, intent(in) :: compute_rank
      logical, intent(in) :: log_it
         !! Log the resolved profile (rank 0).
      integer :: code, ierr, nz
      real(wp), allocatable :: dz(:)

      ocean_state%vcoord%z_fixed_h_ref = cfg%ocean%topo%max_depth
      ocean_state%vcoord%z_fixed_use_profile = .false.
      if (.not. ocean_state%vcoord%is_init) return
      if (ocean_state%vcoord%coord_type /= VCOORD_Z_FIXED) return
      code = parse_z_fixed_profile(cfg%z_fixed_profile)
      if (code == ZFIXED_PROFILE_UNIFORM .or. code == ZFIXED_PROFILE_INVALID) return
      nz = ocean_state%vcoord%nz_ml
      allocate (dz(nz))
      call z_fixed_nominal_dz(code, nz, cfg%ocean%topo%max_depth, cfg%z_fixed_dz, &
                              cfg%z_fixed_dz_top, cfg%z_fixed_tanh_center, &
                              cfg%z_fixed_tanh_width, dz, ierr)
      if (ierr /= ZFIXED_DZ_OK) return
      call ocean_vcoord_set_z_fixed_profile(ocean_state%vcoord, dz)
      if (log_it .and. compute_rank == 0) then
         call logger%info("z_fixed profile:  "//trim(cfg%z_fixed_profile)// &
                          " — nominal dz "//to_string(dz(1))//" m (surface) … "// &
                          to_string(dz(nz))//" m (bed), total "// &
                          to_string(ocean_state%vcoord%z_fixed_h_ref)//" m over "// &
                          to_string(nz)//" layers")
      end if
   end subroutine configure_ocean_z_fixed_profile

   subroutine configure_ocean_k_top(cfg, ocean_state, grid, compute_rank)
      !! Fill `ms%k_top` / `k_top_u` / `k_top_v` — the shared index of
      !! the first LIVE layer counting down from the top, and the field
      !! every top-side consumer reads instead of spelling `nz`.
      !!
      !! Under a quasi-geopotential coordinate beneath an ice shelf
      !! (`vcoord_type = "z_fixed"` with `vcoord%z_top > 0`) the layers
      !! whose nominal geopotential range lies INSIDE the ice are inert
      !! fillers at `zstar_h_min`, so on an ice-covered column
      !! `k = nz` is NOT the ice-adjacent layer.  This is the only
      !! producer of that index.
      !!
      !! **Static.** The pattern is read ONCE, from
      !! `ocean_vcoord_z_fixed_target` at `eta = 0` — the same kernel,
      !! the same `bt_H_ref` and the same `z_top` the closed-face mask
      !! is built from, so there is exactly ONE definition of "live" in
      !! the tree and `test_ocean_ktop` can assert the two agree.  Under
      !! `z_fixed` `eta` is absorbed by the first live layer (the partial
      !! top cell) and a filler's target is `zstar_h_min` whatever `eta`
      !! does, so the live/filler pattern does not move and there is
      !! nothing to recompute per stage.
      !!
      !! **Every other coordinate is a literal no-op**: the arrays were
      !! allocated at `source = nz_ml` in `multilayer_state_init`, which
      !! IS the answer wherever nothing vanishes against the top, so the
      !! whole phase is bit-identical off `z_fixed`.
      !!
      !! **Ordering.** Same as `configure_ocean_closed_faces` — AFTER
      !! `configure_ocean_cavity` (`vcoord%z_top`), AFTER
      !! `configure_ocean_bt_split` (`bt_H_ref`), AFTER the periodic-wrap
      !! / halo pass, and BEFORE `ocean_state_enter_data` (the host fill
      !! is what the `copyin` captures).
      type(config_t), intent(in) :: cfg
      type(ocean_state_t), intent(inout) :: ocean_state
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: compute_rank

      integer :: nx, ny, nz, n_filler
      real(wp) :: h_nominal, h_min
      real(wp), allocatable :: tgt(:, :, :), eta0(:, :)

      if (.not. ocean_state%multilayer%is_init) return
      if (ocean_state%vcoord%coord_type /= VCOORD_Z_FIXED) return
      if (ocean_state%vcoord%z_fixed_h_ref <= 0.0_wp) return
      if (.not. ocean_state%dyn%bt_work%is_init) return
      if (.not. cfg%ocean%cavity_dyn%enable) return

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ocean_state%multilayer%nz_ml
      h_nominal = ocean_state%vcoord%z_fixed_h_ref/real(nz, wp)
      h_min = ocean_state%vcoord%zstar_h_min

      allocate (tgt(nx, ny, nz), source=0.0_wp)
      allocate (eta0(nx, ny), source=0.0_wp)
      call ocean_vcoord_z_fixed_target(tgt, ocean_state%dyn%bt_work%bt_H_ref, &
                                       eta0, ocean_state%vcoord%z_top, &
                                       nx, ny, nz, h_nominal, &
                                       ocean_state%vcoord%z_fixed_use_profile, &
                                       ocean_state%vcoord%z_fixed_zi, &
                                       ocean_state%vcoord%z_fixed_dz, h_min)
      call ocean_vcoord_k_top_from_target(ocean_state%multilayer%k_top, &
                                          ocean_state%multilayer%k_top_u, &
                                          ocean_state%multilayer%k_top_v, &
                                          tgt, nx, ny, nz, H_VANISHED)
      n_filler = count(ocean_state%multilayer%k_top < nz)
      deallocate (tgt, eta0)

      if (compute_rank == 0) then
         call logger%info("z_fixed k_top: "//to_string(n_filler)//"/"// &
                          to_string(nx*ny)//" columns carry top-side "// &
                          "fillers (k_top < nz); min k_top = "// &
                          to_string(minval(ocean_state%multilayer%k_top))// &
                          ", nz = "//to_string(nz))
      end if
   end subroutine configure_ocean_k_top

   subroutine configure_ocean_closed_faces(cfg, ocean_state, grid, compute_rank, ierr)
      !! Build the static partial-step z-level FACE-CLOSURE mask
      !! (`&vcoord_nml zfixed_closed_faces`; Adcroft, Hill & Marshall
      !! 1997; Losch 2008 §2.1 for the ice-shelf cavity).
      !!
      !! Under `vcoord_type = "z_fixed"` a layer whose nominal
      !! geopotential range lies inside the bed — or inside the ice draft
      !! — carries an inert FILLER of thickness `zstar_h_min`
      !! (`<= H_VANISHED`).  A velocity face at which layer `k` is a
      !! filler on EITHER side stays OPEN today, and the FV pressure
      !! gradient evaluated across that staircase step drives
      !! `a = |ρ′|·g·Δz_step/(ρ₀·dx)` out of a resting stratified state —
      !! independent of the filler thickness, so no `h`-gate can reach
      !! it.  A z-LEVEL model treats such a face as a WALL for that
      !! layer: no normal velocity, no mass or tracer flux, free-slip.
      !!
      !! This fills `metrics%open_u/open_v` with that wall, ONCE, from
      !! `ocean_vcoord_z_fixed_target` at `η = 0` — the same kernel the
      !! ALE regrid and the IC seed use, so there is no second definition
      !! of "live".  The mask is STATIC: the bed and the draft are
      !! static, and under `z_fixed` `η` is absorbed by the first LIVE
      !! layer (the partial cell), so the live/filler pattern does not
      !! move.
      !!
      !! It also seeds the barotropic face widths `dy_cu_bt`/`dx_cv_bt`
      !! with the OPEN-depth fraction of the face, so the barotropic
      !! solve is not blind to the closed layers.  That seed is refreshed
      !! from the LIVE `h` every outer step by `ocean_porous_refresh`;
      !! this is only the `η = 0` value the first stage reads.
      !!
      !! **Ordering.** MUST run AFTER `configure_ocean_cavity` (which
      !! fills `vcoord%z_top`), after `configure_ocean_bt_split` (which
      !! lays `bt_H_ref`) and after `configure_ocean_land_mask` and the
      !! periodic-wrap / halo pass (the target is built from GHOST-FILLED
      !! `bt_H_ref` and `z_top`, which is what makes the mask correct at
      !! a periodic seam — the physical seam face is an interior index).
      !! And BEFORE `ocean_state_enter_data`: the host fill is what the
      !! `copyin` captures, and `metrics_closed_faces_alloc` reallocs.
      type(config_t), intent(in) :: cfg
      type(ocean_state_t), intent(inout) :: ocean_state
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: compute_rank
      integer, intent(out), optional :: ierr
         !! Non-zero on a configuration conflict when present; absent
         !! behaves as today (`error stop`).

      integer :: nx, ny, nz, i, j, k
      integer :: n_closed_u, n_closed_v, n_open_u, n_open_v, n_ledge
      real(wp) :: h_nominal, h_min, h_face, sum_all, sum_open
      real(wp), allocatable :: tgt(:, :, :), eta0(:, :)

      if (present(ierr)) ierr = OCEAN_STATUS_OK
      if (.not. cfg%zfixed_closed_faces) return

      if (.not. ocean_state%multilayer%is_init) then
         call fail("&vcoord_nml zfixed_closed_faces requires the ocean "// &
                   "multilayer path (the mask is per-layer)", ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if
      if (ocean_state%vcoord%coord_type /= VCOORD_Z_FIXED) then
         call fail("&vcoord_nml zfixed_closed_faces is only defined for "// &
                   "vcoord_type='z_fixed': the live/filler staircase it "// &
                   "closes is a z-level artefact, and on a terrain-following "// &
                   "or z* family every layer is live on every wet face", &
                   ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if
      if (ocean_state%vcoord%z_fixed_h_ref <= 0.0_wp) then
         call fail("&vcoord_nml zfixed_closed_faces needs a resolved "// &
                   "z_fixed_h_ref (set &ocean_topo_nml max_depth): without "// &
                   "it the z_fixed target degenerates to uniform sigma, "// &
                   "there are no fillers, and the mask would close nothing", &
                   ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if
      if (ocean_state%dyn%bt_work%is_init) then
         if (maxval(ocean_state%dyn%bt_work%bt_H_ref) <= 0.0_wp) then
            call fail("&vcoord_nml zfixed_closed_faces: bt_H_ref is empty — "// &
                      "the barotropic datum must be built (&ocean_bt_nml "// &
                      "n_inner >= 1) before the static face mask can be laid", &
                      ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
      end if
      if (cfg%ocean%wetdry%enable) then
         call fail("&vcoord_nml zfixed_closed_faces is incompatible with "// &
                   "&ocean_wetdry_nml enable: wet/dry moves the live/filler "// &
                   "pattern under the running state and the mask is static", &
                   ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if
      if (cfg%ocean%bt%bt_halo > 0) then
         call fail("&vcoord_nml zfixed_closed_faces is incompatible with "// &
                   "&ocean_bt_nml bt_halo > 0: the wide-halo barotropic "// &
                   "clone carries its own metrics and no face mask, so the "// &
                   "wide BT loop would transport through closed faces "// &
                   "(the same argument &ocean_porous_nml already makes)", &
                   ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if
      if (cfg%ocean%gm%enable .or. cfg%ocean%redi%enable .or. &
          cfg%ocean%foxkemper%enable) then
         call fail("&vcoord_nml zfixed_closed_faces does not yet compose "// &
                   "with GM / Redi / MLE: all three form their face fluxes "// &
                   "from a 2-D `wet_u`/`wet_v` gate and fold them into "// &
                   "mass_flux_*_layer AFTER continuity has applied the "// &
                   "per-layer mask, so their transports would leak through "// &
                   "a closed face.  Per-layer seams in those three kernels "// &
                   "are the follow-up slice", ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if
      if (cfg%ocean%hvisc%nu_4 > 0.0_wp .or. cfg%ocean%hvisc%stress_tensor) then
         call fail("&vcoord_nml zfixed_closed_faces does not yet compose "// &
                   "with the BIHARMONIC viscosity (&ocean_hvisc_nml nu_4) or "// &
                   "with stress_tensor: the free-slip closure of a closed "// &
                   "face is implemented for the HARMONIC velocity-Laplacian "// &
                   "kernels only (scalar nu_h and the per-face ah_face_*), "// &
                   "and a mirror-Neumann biharmonic pass or a wet_q-style "// &
                   "corner factor needs its own derivation and its own "// &
                   "test.  Use a harmonic closure, or land that slice first", &
                   ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if
      ! ---- Barotropic paths still on FULL-column weights ----
      !
      ! `derive_bt_from_layers`, `face_depth_mean_*`, `set_cor_ref_velocity`
      ! and the `apply_bt_correction` fold weight by `h_face·open` under
      ! this knob.  Three barotropic paths do NOT: each sums
      ! `0.5·(h_L + h_R)` over EVERY layer.  Their answers are not wrong
      ! by a round-off — the closed layers' thickness is O(h_nominal) at
      ! a staircase face — so they are refused until they are ported,
      ! rather than left to run on the wrong column.
      if (cfg%ocean%bt%correction_bc_pgf) then
         call fail("&vcoord_nml zfixed_closed_faces does not yet compose "// &
                   "with &ocean_bt_nml correction_bc_pgf: compute_pbce, "// &
                   "compute_gtot_faces and the bc-PGF du_bc block in "// &
                   "apply_bt_correction all weight by the FULL column "// &
                   "(0.5*(h_L+h_R) over every layer), so the correction's "// &
                   "depth-mean-zero identity holds on the full column, not "// &
                   "on the OPEN column ubt is the mean of — it would leak a "// &
                   "barotropic increment into the open layers and push one "// &
                   "into the closed ones.  Set correction_bc_pgf=.false.", &
                   ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if
      if (cfg%ocean%bt%substep_drag) then
         call fail("&vcoord_nml zfixed_closed_faces does not yet compose "// &
                   "with &ocean_bt_nml substep_drag: compute_bt_rem builds "// &
                   "the per-face damping Htot/(Htot + r*hbbl*dt_inner) from "// &
                   "the FULL-column face depth, while the barotropic "// &
                   "transport under this knob runs on the OPEN column only, "// &
                   "so a partially closed face is under-damped by the "// &
                   "closed fraction.  Set substep_drag=.false.", &
                   ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if
      if (cfg%ocean%bt%wave_drag) then
         call fail("&vcoord_nml zfixed_closed_faces does not yet compose "// &
                   "with &ocean_bt_nml wave_drag: compute_bt_rem_wave_drag "// &
                   "multiplies Htot/(Htot + r_H*dt_inner) into bt_rem from "// &
                   "the FULL-column face depth, while the barotropic "// &
                   "transport under this knob runs on the OPEN column only, "// &
                   "so a partially closed face is under-damped by the "// &
                   "closed fraction.  Set wave_drag=.false.", &
                   ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ocean_state%multilayer%nz_ml
      h_nominal = ocean_state%vcoord%z_fixed_h_ref/real(nz, wp)
      h_min = ocean_state%vcoord%zstar_h_min

      call metrics_closed_faces_alloc(ocean_state%metrics, grid, nz)

      allocate (tgt(nx, ny, nz), source=0.0_wp)
      allocate (eta0(nx, ny), source=0.0_wp)
      call ocean_vcoord_z_fixed_target(tgt, ocean_state%dyn%bt_work%bt_H_ref, &
                                       eta0, ocean_state%vcoord%z_top, &
                                       nx, ny, nz, h_nominal, &
                                       ocean_state%vcoord%z_fixed_use_profile, &
                                       ocean_state%vcoord%z_fixed_zi, &
                                       ocean_state%vcoord%z_fixed_dz, h_min)
      call ocean_vcoord_closed_face_masks(ocean_state%metrics%open_u, &
                                          ocean_state%metrics%open_v, &
                                          tgt, nx, ny, nz, H_VANISHED)

      ! Seed the BAROTROPIC face widths with the eta = 0 open-depth
      ! fraction.  `ocean_porous_refresh` recomputes them from the live
      ! `h` every outer step; this is only what stage 1 of step 1 reads.
      ! The width and the depth are the same number here: the BT
      ! transport is `ubt * FA * dy_cu_bt` with `FA = sum_k h_face`, so
      ! `dy_cu_bt = dy_cu * (sum_open h_face)/(sum_k h_face)` makes that
      ! product identically `ubt * (sum_open h_face) * dy_cu`.
      do j = 1, ny
         do i = 2, nx
            sum_all = 0.0_wp
            sum_open = 0.0_wp
            do k = 1, nz
               h_face = 0.5_wp*(tgt(i - 1, j, k) + tgt(i, j, k))
               sum_all = sum_all + h_face
               if (ocean_state%metrics%open_u(i, j, k) > 0.0_wp) then
                  sum_open = sum_open + h_face
               end if
            end do
            if (sum_all > 0.0_wp) then
               ocean_state%metrics%dy_cu_bt(i, j) = &
                  ocean_state%metrics%dy_cu(i, j)*(sum_open/sum_all)
            end if
         end do
      end do
      do j = 2, ny
         do i = 1, nx
            sum_all = 0.0_wp
            sum_open = 0.0_wp
            do k = 1, nz
               h_face = 0.5_wp*(tgt(i, j - 1, k) + tgt(i, j, k))
               sum_all = sum_all + h_face
               if (ocean_state%metrics%open_v(i, j, k) > 0.0_wp) then
                  sum_open = sum_open + h_face
               end if
            end do
            if (sum_all > 0.0_wp) then
               ocean_state%metrics%dx_cv_bt(i, j) = &
                  ocean_state%metrics%dx_cv(i, j)*(sum_open/sum_all)
            end if
         end do
      end do

      ! Cheap one-shot census for the configure line: how much of the
      ! array the mask actually closes, and whether it isolated any
      ! water.  Host-side, once, O(nx*ny*nz) — not a per-step diagnostic.
      n_closed_u = 0
      n_open_u = 0
      do k = 1, nz
         do j = 1, ny
            do i = 2, nx
               if (ocean_state%metrics%open_u(i, j, k) > 0.0_wp) then
                  n_open_u = n_open_u + 1
               else
                  n_closed_u = n_closed_u + 1
               end if
            end do
         end do
      end do
      n_closed_v = 0
      n_open_v = 0
      do k = 1, nz
         do j = 2, ny
            do i = 1, nx
               if (ocean_state%metrics%open_v(i, j, k) > 0.0_wp) then
                  n_open_v = n_open_v + 1
               else
                  n_closed_v = n_closed_v + 1
               end if
            end do
         end do
      end do
      n_ledge = ocean_vcoord_count_ledges(ocean_state%metrics%open_u, &
                                          ocean_state%metrics%open_v, &
                                          tgt, nx, ny, nz, H_VANISHED)

      deallocate (tgt, eta0)

      ! Flip the switches LAST — every consumer branches on them, and the
      ! mask has to be in place before any of them can read a closed face.
      ocean_state%metrics%use_closed_faces = .true.
      ocean_state%vcoord%zfixed_closed_faces = .true.
      ocean_state%vdiff%zlevel_faces = .true.

      if (compute_rank == 0) then
         call logger%info("z_fixed closed faces: ON (partial steps, "// &
                          "Adcroft/Hill/Marshall 1997) — u closed "// &
                          to_string(n_closed_u)//"/"// &
                          to_string(n_closed_u + n_open_u)//", v closed "// &
                          to_string(n_closed_v)//"/"// &
                          to_string(n_closed_v + n_open_v)// &
                          ", isolated ledge cells "//to_string(n_ledge))
         if (n_ledge > 0) then
            call logger%warning("z_fixed closed faces: "//to_string(n_ledge)// &
                                " LIVE cells have all four own-layer faces "// &
                                "closed — the mask has isolated water.  They "// &
                                "are inert (no flux in or out, velocity zeroed "// &
                                "every stage) but a one-cell spike in the bed "// &
                                "or the draft is worth looking at.")
         end if
      end if
      if (present(ierr)) ierr = OCEAN_STATUS_OK
   end subroutine configure_ocean_closed_faces

   subroutine configure_ocean_porous(cfg, ocean_state, grid, compute_rank, ierr)
      !! Configure porous barriers (`&ocean_porous_nml`, Adcroft 2013).
      !! Grows the `ocean_metrics_t` porous arrays to full face size and
      !! fills the STATIC along-face `d_min`/`d_max`/`d_avg` statistics;
      !! the layer-averaged open fractions themselves are recomputed on
      !! the device every RK2 stage (they depend on the interface
      !! heights).
      !!
      !! MUST run AFTER bathymetry is set (`ocean_state%barotropic%b`)
      !! and land masking, and BEFORE `ocean_state_enter_data` — the
      !! realloc has to happen before the device map, and the host fill
      !! is what the `copyin` captures.  Leaves the placeholder-sized
      !! arrays untouched when the knob is off (bit-identical).
      type(config_t), intent(in) :: cfg
      type(ocean_state_t), intent(inout) :: ocean_state
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: compute_rank
      integer, intent(out), optional :: ierr
         !! Non-zero on a porous-barrier configuration conflict when
         !! present; absent behaves as today (`error stop`).

      integer :: src, interp, nx, ny, nz

      if (.not. cfg%ocean%porous%enable) then
         if (present(ierr)) ierr = OCEAN_STATUS_OK
         return
      end if

      if (.not. ocean_state%multilayer%is_init) then
         call fail("&ocean_porous_nml enable requires the ocean multilayer "// &
                   "path (the open fractions are per-layer)", ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if

      src = parse_porous_source(cfg%ocean%porous%source)
      select case (src)
      case (POROUS_SOURCE_RESOLVED)
         continue
      case (POROUS_SOURCE_FILE)
         call fail("&ocean_porous_nml source='file' (an offline "// &
                   "subgrid-bathymetry file, MOM6 topog_edge.nc) is not "// &
                   "implemented — it needs the file-forcing backend that "// &
                   "&ocean_bt_nml wave_drag_form='file' also waits on. "// &
                   "Use source='resolved' (a documented proxy).", ierr, OCEAN_STATUS_ERR_SETUP)
         return
      case default
         call fail("&ocean_porous_nml source='"// &
                   trim(cfg%ocean%porous%source)//"' is not recognised", ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end select

      interp = parse_porous_eta_interp(cfg%ocean%porous%eta_interp)
      if (interp < 0) then
         call fail("&ocean_porous_nml eta_interp='"// &
                   trim(cfg%ocean%porous%eta_interp)//"' is not recognised", ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ocean_state%multilayer%nz_ml

      call metrics_porous_alloc(ocean_state%metrics, grid, nz)
      ocean_state%metrics%porous_eta_interp = interp
      ! The namelist knob is a DEPTH (positive below the surface); the
      ! kernel gate compares HEIGHTS, so flip the sign once here.
      ocean_state%metrics%porous_mask_depth = -cfg%ocean%porous%masking_depth

      ! SIGN CONVENTION.  The Adcroft fit works in topographic HEIGHTS
      ! (positive up, negative below the sea surface) so that they compare
      ! directly against interface heights.  On the ocean path
      ! `barotropic%b` holds the reference column DEPTH, positive down
      ! (`SSH = sum(h_layer) - b`), so it is negated once, here, and every
      ! porous array downstream is a height.  Bathymetry is static, so a
      ! snapshot on the metrics slot is exact and lets the per-step
      ! recompute stay a metrics-only kernel.
      ocean_state%metrics%por_bed = -ocean_state%barotropic%b

      call porous_fill_stats_resolved(nx, ny, ocean_state%metrics%por_bed, &
                                      ocean_state%metrics%wet_T, &
                                      ocean_state%metrics%por_dmin_u, &
                                      ocean_state%metrics%por_dmax_u, &
                                      ocean_state%metrics%por_davg_u, &
                                      ocean_state%metrics%por_dmin_v, &
                                      ocean_state%metrics%por_dmax_v, &
                                      ocean_state%metrics%por_davg_v)

      ! READER-BOUNDARY ASSERTION.  `m = (d_avg-d_min)/(d_max-d_min)` is
      ! only a valid Adcroft parameter when `d_min <= d_avg <= d_max`.
      ! The resolved fill cannot break that (d_avg is a convex combination
      ! of the samples d_min/d_max bracket), but a FILE-backed source
      ! could, so the check sits here — at the boundary every source
      ! passes through — rather than inside one filler.
      if (.not. porous_stats_are_ordered(nx + 1, ny, &
                                         ocean_state%metrics%por_dmin_u, &
                                         ocean_state%metrics%por_dmax_u, &
                                         ocean_state%metrics%por_davg_u)) then
         call fail("&ocean_porous_nml: u-face subgrid statistics violate "// &
                   "d_min <= d_avg <= d_max", ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if
      if (.not. porous_stats_are_ordered(nx, ny + 1, &
                                         ocean_state%metrics%por_dmin_v, &
                                         ocean_state%metrics%por_dmax_v, &
                                         ocean_state%metrics%por_davg_v)) then
         call fail("&ocean_porous_nml: v-face subgrid statistics violate "// &
                   "d_min <= d_avg <= d_max", ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if

      ! Flip the master switch LAST: every kernel branches on it, and the
      ! stats must be in place before any of them can read a narrowed width.
      ocean_state%metrics%use_porous = .true.

      if (compute_rank == 0) then
         call logger%info("Porous barriers:  ON (Adcroft 2013), source='"// &
                          trim(cfg%ocean%porous%source)//"' eta_interp='"// &
                          trim(cfg%ocean%porous%eta_interp)//"' masking_depth="// &
                          to_string(cfg%ocean%porous%masking_depth)//" m")
         call logger%info("                  source='resolved' is an along-face "// &
                          "statistic of the RESOLVED bathymetry, a documented "// &
                          "proxy for true subgrid data")
      end if
      if (present(ierr)) ierr = OCEAN_STATUS_OK
   end subroutine configure_ocean_porous

   subroutine wave_drag_roughness_proxy(b, wet_T, nx, ny, nghost, kappa, n_bot, h2_max, r_h)
      !! `form="roughness_proxy"` filler — a DOCUMENTED PLACEHOLDER for
      !! Jayne & St Laurent (2001)'s subgrid `<h^2>`, not a substitute for
      !! it (that needs PR-14's file reader or PR-30's field-valued
      !! roughness). Estimates the subgrid topographic-height variance from
      !! the RESOLVED 2-delta bathymetry increment:
      !!     <h^2>_proxy(i,j) = 1/4*[(b(i+1,j)-b(i-1,j))^2 + (b(i,j+1)-b(i,j-1))^2]
      !! then `r_H = 1/2*kappa*min(<h^2>_proxy, h2_max)*N_bot`. `b` is
      !! bottom elevation, positive UP (`rdb_barotropic_state.F90`) —
      !! differences are sign-independent. Zero on land (`wet_T==0`) and on
      !! the ghost ring (the 2-delta stencil is unavailable there; a
      !! formula-bathymetry path that leaves ghosts unfilled would
      !! otherwise manufacture a spurious cliff at the ghost seam — CLAUDE.md
      !! "Formula bathymetry setters must fill ghost rows").
      integer, intent(in) :: nx, ny, nghost
      real(wp), intent(in) :: b(nx, ny), wet_T(nx, ny)
      real(wp), intent(in) :: kappa, n_bot, h2_max
      real(wp), intent(inout) :: r_h(nx, ny)
      integer :: i, j
      real(wp) :: h2

      r_h = 0.0_wp
      do j = nghost + 1, ny - nghost
         do i = nghost + 1, nx - nghost
            if (wet_T(i, j) <= 0.0_wp) cycle
            h2 = 0.25_wp*((b(i + 1, j) - b(i - 1, j))**2 + (b(i, j + 1) - b(i, j - 1))**2)
            r_h(i, j) = 0.5_wp*kappa*min(h2, h2_max)*n_bot
         end do
      end do
   end subroutine wave_drag_roughness_proxy

   subroutine parse_constituent_list(list, cat_idx, nconst, ierr)
      !! Tokenize a whitespace/comma-separated constituent list into
      !! catalog indices (case-insensitive).  Unknown tokens fail loud;
      !! `nconst` is the count of recognised entries.
      character(len=*), intent(in) :: list
      integer, intent(out) :: cat_idx(:)
      integer, intent(out) :: nconst
      integer, intent(out), optional :: ierr
         !! Non-zero on an unrecognised constituent token when present;
         !! absent behaves as today (`error stop`).
      integer :: i, n, i0, idx
      character(len=len(list)) :: buf
      logical :: in_tok

      buf = list
      ! Normalise commas to spaces.
      do i = 1, len(buf)
         if (buf(i:i) == ",") buf(i:i) = " "
      end do

      nconst = 0
      n = len_trim(buf)
      in_tok = .false.
      i0 = 1
      do i = 1, n + 1
         if (i <= n .and. buf(i:i) /= " ") then
            if (.not. in_tok) then
               in_tok = .true.
               i0 = i
            end if
         else
            if (in_tok) then
               in_tok = .false.
               idx = tide_name_index(buf(i0:i - 1))
               if (idx < 1) then
                  call fail("ocean_tides_nml: unknown constituent '"// &
                            trim(buf(i0:i - 1))//"'", ierr, OCEAN_STATUS_ERR_SETUP)
                  return
               end if
               nconst = nconst + 1
               cat_idx(nconst) = idx
            end if
         end if
      end do
      if (present(ierr)) ierr = OCEAN_STATUS_OK
   end subroutine parse_constituent_list

   subroutine configure_ocean_wavespeed(cfg, ocean_state, grid, ierr)
      !! Copy the `&ocean_wavespeed_nml` knobs onto the wave-speed slot
      !! (B1).  Diagnostic, no mutual exclusion: cg1/Rd read `rho_layer`
      !! directly.  Fills `f_centre` + the static `beta_centre = |grad
      !! f|` field via `fill_coriolis_centre` + `build_static` — the
      !! SAME `metrics_fill_coriolis` path VarMix/MEKE use (planetary on
      !! spherical/tripolar, beta-plane bit-identical elsewhere) —
      !! rather than the legacy hard-coded beta-plane `set_f_centre`.
      !! Copies rho0 from the EOS slot for the Boussinesq `gprime`.  Must
      !! run AFTER `configure_ocean_metrics` (metrics filled + finalized).
      type(config_t), intent(in) :: cfg
      type(ocean_state_t), intent(inout) :: ocean_state
      type(hgrid_t), intent(in) :: grid
      integer, intent(out), optional :: ierr
         !! Non-zero on a wave-speed configuration conflict when present;
         !! absent behaves as today (`error stop`).
      real(wp), allocatable :: f_centre(:, :)

      associate (ws => ocean_state%wavespeed, wcfg => cfg%ocean%wavespeed)
         ws%enable = wcfg%enable
         if (.not. wcfg%enable) then
            if (present(ierr)) ierr = OCEAN_STATUS_OK
            return
         end if

         ! Fail loud at configure: the column kernel indexes fixed-size
         ! NZ_STACK_MAX column arrays; wavespeed_compute never ran before
         ! this PR so nz > NZ_STACK_MAX never overran it.  It can now.
         if (cfg%nz_layers > NZ_STACK_MAX) then
            call fail("ocean_wavespeed_nml: nz_layers exceeds "// &
                      "NZ_STACK_MAX (raise NZ_STACK_MAX in rdb_constants)", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if

         ws%mono_n2_depth = wcfg%mono_n2
         ws%use_ebt = wcfg%use_ebt
         ws%n_wavespeed = wcfg%n_wavespeed
         ws%rho0 = ocean_state%eos%rho0

         ! |f| at centres + the static |grad f| field for Rd.
         allocate (f_centre(grid%nx_total, grid%ny_total))
         call fill_coriolis_centre(cfg, ocean_state%metrics, grid, f_centre)
         call ws%build_static(grid, ocean_state%metrics, f_centre)
         deallocate (f_centre)

         if (wcfg%n_wavespeed < 1) then
            call fail("ocean_wavespeed_nml: n_wavespeed must be >= 1", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
      end associate
      if (present(ierr)) ierr = OCEAN_STATUS_OK
   end subroutine configure_ocean_wavespeed

   subroutine configure_ocean_varmix(cfg, ocean_state, grid)
      !! Build the STATIC VarMix grid terms (`f2_dx2_*`, `beta_dx2_*`,
      !! `l2_*`) on the varmix slot once at configure time, from the filled
      !! curvilinear metrics + the cell-centre Coriolis magnitude.  The slot
      !! scalars are copied earlier by `ocean_state_copy_config`; this only
      !! fills the static arrays (the per-step Res_fn/SN/assembly fire in the
      !! dyn step).  No-op when VarMix is disabled.
      type(config_t), intent(in) :: cfg
      type(ocean_state_t), intent(inout) :: ocean_state
      type(hgrid_t), intent(in) :: grid
      real(wp), allocatable :: f_centre(:, :)

      if (.not. ocean_state%varmix%enable) return
      allocate (f_centre(grid%nx_total, grid%ny_total))
      call fill_coriolis_centre(cfg, ocean_state%metrics, grid, f_centre)
      call ocean_state%varmix%build_static(ocean_state%metrics, f_centre)
      deallocate (f_centre)
   end subroutine configure_ocean_varmix

   subroutine configure_ocean_meke(cfg, ocean_state, grid)
      !! Fill the MEKE slot's cell-centre |Coriolis| from the same
      !! `metrics_fill_coriolis` path VarMix / EPBL use (handles beta-plane
      !! AND spherical), so `beta = |grad f|` for the Rhines length is live
      !! when `alpha_rhines > 0`.  The MEKE scalar knobs are copied earlier
      !! by `ocean_state_copy_config`; this only fills `f_centre`.  Host loop
      !! before `enter_data`.  No-op when MEKE is disabled.
      type(config_t), intent(in) :: cfg
      type(ocean_state_t), intent(inout) :: ocean_state
      type(hgrid_t), intent(in) :: grid
      real(wp), allocatable :: f_centre(:, :)

      if (.not. ocean_state%meke%enable) return
      allocate (f_centre(grid%nx_total, grid%ny_total))
      call fill_coriolis_centre(cfg, ocean_state%metrics, grid, f_centre)
      call ocean_state%meke%set_f_centre(grid, f_centre)
      deallocate (f_centre)
   end subroutine configure_ocean_meke

   subroutine configure_ocean_reference_density(ocean_state, geo)
      !! Fan the ONE configured Boussinesq reference density out to every
      !! remaining slot that carries its own `rho0` copy.
      !!
      !! `&ocean_ic_nml rho_0` lands on `eos%rho0` in
      !! `ocean_state_init_from_config`; that is the ρ₀ of record.  EPBL,
      !! kappa-shear, tidal mixing, wave speed, the `eta_ib` surface-pressure
      !! seam, GM / MEKE / Redi / MLE, the isopycnal slopes and (since the
      !! preceding commit) the PGF all copy it in their own
      !! `configure_ocean_*`.  The four slots wired here were the remainder:
      !! they kept a hard 1035 type default that nothing ever assigned, so a
      !! namelist with `rho_0 /= 1035` ran the EOS on one reference density
      !! and the surface forcing, wind stress and KPP buoyancy on another —
      !! no warning, no fail-loud, just a 1035/ρ₀ scaling on every surface
      !! heat/salt flux, every wind-stress acceleration, and N²/u*/B_0.
      !!
      !!   * `surface_flux%rho0`  — the `dt/(ρ₀·cp)` heat and `dt/ρ₀` salt
      !!     divisors for EVERY surface tracer source, including what the
      !!     sea-ice coupler delivers through `Q_heat`/`Q_salt`.
      !!   * `surface_stress%rho0` — the `τ/(ρ₀·h_top)` acceleration (both
      !!     the top-layer and the DIRECT_STRESS distributed form).
      !!   * `vmix%rho0` — KPP: N² = −g/ρ₀·∂ρ/∂z, u* = √(|τ|/ρ₀), the
      !!     kinematic surface fluxes q_T = Q_heat/(ρ₀·cp), q_S = Q_salt/ρ₀
      !!     feeding B_0, and the PP81 / convective-adjustment N².
      !!   * `geothermal%rho0` — the `dt·Q_geo/(ρ₀·cp)` bed heat source.
      !!     The geothermal slot lives on the engine, not on `ocean_state`,
      !!     so it is passed in (optional: a caller with no geothermal slot
      !!     simply omits it).
      !!
      !! Device contract (`mem:separate`): `surface_flux`, `surface_stress`
      !! and `geothermal` read their `rho0` HOST-side (folded into the
      !! `inv_scale` / `src_T` scalar or passed by value into a `*_impl`), so
      !! those owe nothing.  `vmix%rho0` IS read on-device — `this%rho0`
      !! appears inside the `do concurrent` bodies of `vmix_compute_pp81`,
      !! `vmix_kpp_overlay_impl` and `vmix_convective_impl` — but this
      !! routine runs in the configure phase, strictly BEFORE
      !! `ocean_state_enter_data`'s `copyin`, exactly like the neighbouring
      !! `pp81_*` / `shear2_floor` scalars it sits with.  No
      !! `!$acc update device` is owed.  **A configure step that ever moves
      !! after `enter_data` must add one.**
      type(ocean_state_t), intent(inout) :: ocean_state
      type(ocean_geothermal_t), intent(inout), optional :: geo
         !! Engine-held geothermal slot (the split driver's `geo` argument).

      ocean_state%surface_flux%rho0 = ocean_state%eos%rho0
      ocean_state%surface_stress%rho0 = ocean_state%eos%rho0
      ocean_state%vmix%rho0 = ocean_state%eos%rho0
      if (present(geo)) geo%rho0 = ocean_state%eos%rho0
   end subroutine configure_ocean_reference_density

   subroutine configure_ocean_pgf(cfg, ocean_state, compute_rank, ierr)
      !! Pressure-force variant, the reference densities (`rho0` / `rho_ref`,
      !! both from the single configured ρ₀ — `&ocean_ic_nml rho_0` via
      !! `eos%rho0`), the reduced-gravity (gprime / gfs_scale) knobs, the
      !! bathymetry copy into the PGF slot, and the matching barotropic
      !! fast-loop gravity g_bt for the gprime / FV_MOM6-reduced-GFS paths.
      type(config_t), intent(in) :: cfg
      type(ocean_state_t), intent(inout) :: ocean_state
      integer, intent(in) :: compute_rank
         !! Unused (no rank-0 logging in this helper); kept for a uniform
         !! configure_ocean_* signature.
      integer, intent(out), optional :: ierr
         !! Non-zero on a PGF/EOS configuration conflict when present;
         !! absent behaves as today (`error stop`).
      integer :: local_ierr

      ! ALLOCATION-GATE CONSISTENCY GUARD.  When `scratch_gated` is on, the
      ! PGF slot allocated only the buffers reachable from the `variant` /
      ! `reconstruct_for_pressure` latched by `ocean_state_init_from_config`
      ! BEFORE `init(grid)`.  Both are re-derived here from the same cfg
      ! fields, so they must agree.  If a future edit ever lets them
      ! diverge, the gated-off buffers would be unallocated on a path that
      ! reaches them — which under `-gpu=mem:separate` is a silent wrong
      ! answer, not a crash.  Fail loud instead.
      if (ocean_state%pressure_force%scratch_gated) then
         if (ocean_state%pressure_force%variant /= &
             parse_opgf_variant(cfg%ocean%pgf%form)) then
            call fail("configure_ocean_pgf: pgf variant changed between "// &
                      "the pre-init latch and configure — the scratch "// &
                      "allocation gate was decided on the stale value.", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
         if (ocean_state%pressure_force%reconstruct_for_pressure .neqv. &
             cfg%ocean%pgf%reconstruct_for_pressure) then
            call fail("configure_ocean_pgf: reconstruct_for_pressure "// &
                      "changed between the pre-init latch and configure — "// &
                      "the recon scratch allocation gate is stale.", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
      end if
      ocean_state%pressure_force%variant = parse_opgf_variant(cfg%ocean%pgf%form)

      ! Grounded-layer PGF gate (`&ocean_isopycnal_nml pgf_skip_nonoverlap`).
      ! The predicate lives in `pgf_nonoverlap_gate_on` because the PGF slot's
      ! allocation gate had to evaluate it BEFORE `init` (FV_MOM6 allocates
      ! `z_centre` for this gate and nothing else) — see the pre-init latch in
      ! `ocean_state_init_from_config`.  Re-evaluating it here and
      ! comparing is what keeps that latch honest.
      if (ocean_state%pressure_force%scratch_gated .and. &
          (ocean_state%pressure_force%skip_nonoverlap .neqv. &
           pgf_nonoverlap_gate_on(cfg, ocean_state%pressure_force%variant))) then
         call fail("configure_ocean_pgf: pgf_skip_nonoverlap changed between "// &
                   "the pre-init latch and configure — the z_centre allocation "// &
                   "gate is stale.", ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if
      ocean_state%pressure_force%skip_nonoverlap = &
         pgf_nonoverlap_gate_on(cfg, ocean_state%pressure_force%variant)
      if (ocean_state%pressure_force%skip_nonoverlap) then
         if (compute_rank == 0) then
            call logger%info("Isopycnal PGF:    pgf_skip_nonoverlap = ON "// &
                             "(grounded-layer face PGF zeroed where layers do not overlap)")
         end if
      else if (parse_ocean_vcoord_type(cfg%vcoord_type) == VCOORD_LAGRANGIAN .and. &
               cfg%ocean%isopycnal%pgf_skip_nonoverlap .and. &
               ocean_state%pressure_force%variant /= OPGF_VARIANT_GPRIME) then
         ! Asked for, but not available.  What is left here is `gprime`: the
         ! NK=2 reduced-gravity form differences interface positions directly,
         ! never builds the layer-centre Jacobian the gate corrects and carries
         ! no `z_centre` buffer — enabling the gate would read unallocated
         ! memory, a silent wrong answer under `-gpu=mem:separate` rather than
         ! a crash.  Loud, not fatal: the run is legal, it just keeps the
         ! spurious gradient.
         if (compute_rank == 0) then
            call logger%warning("ocean_isopycnal_nml: pgf_skip_nonoverlap has no "// &
                                "effect with ocean_pgf_nml form='"// &
                                trim(cfg%ocean%pgf%form)//"' (the gate covers mont "// &
                                "/ fv_lite / fv_wright / fv_mom6) — a grounded "// &
                                "isopycnal layer can still drive a spurious "// &
                                "pressure gradient at rest")
         end if
      end if

      ! N2: device-callable-variant gate runs POST-config (eos%variant was
      ! set from the knob in configure_ocean_drag, the first configure call).
      if (present(ierr)) then
         call eos_validate(ocean_state%eos, ierr=local_ierr)
         if (local_ierr /= 0) then
            ierr = local_ierr
            return
         end if
      else
         call eos_validate(ocean_state%eos)
      end if
      ! FV_WRIGHT re-evaluates the Wright (1997) rational EOS inside its
      ! Picard pressure sweep, so it cannot honour a non-Wright rho_layer.
      ! Roquet + FV_WRIGHT is therefore unsupported (FV_LITE / FV_MOM6 /
      ! gprime read rho_layer generically and are fine).  Fail loud — a
      ! follow-up could add a Roquet column sweep.
      if (ocean_state%eos%variant == EOS_VARIANT_ROQUET_SPV .and. &
          ocean_state%pressure_force%variant == OPGF_VARIANT_FV_WRIGHT) then
         call fail("configure_ocean_pgf: eos='roquet_spv' is "// &
                   "incompatible with pgf form='fv_wright' (the FV_WRIGHT "// &
                   "Picard sweep re-evaluates Wright internally). Use "// &
                   "fv_lite, fv_mom6, or gprime with roquet_spv.", ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if

      ! Reference densities.  BOTH PGF reference densities follow the SINGLE
      ! configured ρ₀ (`&ocean_ic_nml rho_0`, landed on `eos%rho0` by
      ! `ocean_state_init_from_config` — the same scalar EPBL, kappa-shear,
      ! tidal mixing, the `eta_ib` surface-pressure seam, MEKE/GM and the
      ! isopycnal slopes all take).  Without this the PGF slot kept its
      ! 1035 type default while the EOS followed the namelist, so a run with
      ! `rho_0 /= 1035` silently integrated an EOS and a pressure gradient on
      ! two different reference densities.
      !
      ! `rho0` (Boussinesq divisor, `du/dt = −(1/ρ₀)∂p/∂x`) and `rho_ref`
      ! (the baseline subtracted from layer densities when building the
      ! FV_MOM6 `pa` anomaly stack, and the `g·ρ_ref/ρ₀` surface value in
      ! `compute_pbce`) stay SEPARATE MEMBERS — the roles differ and
      ! interchanging them is a known MOM6 bug class — but roundabout has no
      ! separate anomaly-reference knob, so both take the one configured ρ₀.
      !
      ! Plain host scalars: every consumer reads them host-side (into a local
      ! `inv_rho0`, or by value into a `*_impl`), so no `!$acc update device`
      ! is owed here — and this runs well before `ocean_state_enter_data` in
      ! any case.
      ocean_state%pressure_force%rho0 = ocean_state%eos%rho0
      ocean_state%pressure_force%rho_ref = ocean_state%eos%rho0

      ocean_state%pressure_force%gprime_gfs = cfg%ocean%pgf%gprime_gfs
      ocean_state%pressure_force%gprime_gint = cfg%ocean%pgf%gprime_gint
      ocean_state%pressure_force%gfs_scale = cfg%ocean%pgf%gfs_scale
      ocean_state%pressure_force%mass_weight = cfg%ocean%pgf%mass_weight
      ocean_state%pressure_force%reconstruct_for_pressure = &
         cfg%ocean%pgf%reconstruct_for_pressure
      ocean_state%pressure_force%recon_scheme = cfg%ocean%pgf%recon_scheme
      ocean_state%pressure_force%p_top_in_bc = cfg%ocean%pgf%p_top_in_bc
      ! The top-of-column load enters the `pa` stack's surface boundary
      ! condition, which only the FV_MOM6 family builds.  Mirrors the
      ! `validate_config` refusal so an in-memory namelist / API caller
      ! that bypasses validation still fails loud instead of running a
      ! silently inert knob.
      if (cfg%ocean%pgf%p_top_in_bc .and. &
          ocean_state%pressure_force%variant /= OPGF_VARIANT_FV_MOM6) then
         call fail("configure_ocean_pgf: p_top_in_bc=.true. requires "// &
                   "form='fv_mom6'. Other PGF variants build no pa(nz+1) "// &
                   "pressure-stack boundary condition for the load to enter.", &
                   ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if
      ! In-layer T/S reconstruction wires into the FV_MOM6 layer-integrated
      ! form only (it carries e_face / pa / intz_dpa).  Fail loud if the
      ! knob is on with any other PGF form.
      if (cfg%ocean%pgf%reconstruct_for_pressure .and. &
          ocean_state%pressure_force%variant /= OPGF_VARIANT_FV_MOM6) then
         call fail("configure_ocean_pgf: reconstruct_for_pressure=.true. "// &
                   "requires form='fv_mom6' (the layer-integrated FV path). "// &
                   "Other PGF variants do not carry the e_face/pa/intz_dpa "// &
                   "stack the reconstruction replaces.", ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if
      ! The PGF's bathymetry copy (gprime recovers ∇η = ∇(Σh) - ∇b, FV-MOM6
      ! places the bottom interface from it) is NOT taken here: `engine_setup`
      ! takes it after the init-time periodic/fold wrap + halo exchange of
      ! `b`, so the seam ghosts it copies are the right ones.
      ! gprime: run the BT fast loop at the reduced free-surface gravity g_FS,
      ! else the BT correction cancels the gprime reduced-gravity tendency.
      if (ocean_state%pressure_force%variant == OPGF_VARIANT_GPRIME) then
         ocean_state%dyn%bt_work%g_bt = ocean_state%pressure_force%gprime_gfs
      end if
      ! BEBT projection (MOM6 BT_PROJECT_VELOCITY); 0 = pure forward-backward.
      ocean_state%dyn%bt_work%bebt = cfg%ocean%bt%bebt
      ! FV_MOM6 with reduced GFS: drop BT gravity to gfs_scale·g (the Pass-5
      ! Montgomery correction supplies the complementary share).  Gate must
      ! NOT condition on nz_layers (NK=1 SSH bug, 2026-05-26).
      if (ocean_state%pressure_force%variant == OPGF_VARIANT_FV_MOM6 .and. &
          ocean_state%pressure_force%gfs_scale < 1.0_wp - 1.0e-12_wp) then
         ocean_state%dyn%bt_work%g_bt = ocean_state%pressure_force%gfs_scale*GRAVITY
      end if
      if (present(ierr)) ierr = OCEAN_STATUS_OK
   end subroutine configure_ocean_pgf

   subroutine configure_ocean_bt(cfg, ocean_state, grid, compute_rank)
      !! Barotropic-substep correction knobs (MOM6 frhatu h-weighting, bc-PGF
      !! retro-correction, bt_rem_u drag damping, visc_rem joint weight), the
      !! BT_cont_type / upstream-PPM h_face workspace allocations, and the
      !! rank-0 PGF/BT configuration log lines.
      type(config_t), intent(in) :: cfg
      type(ocean_state_t), intent(inout) :: ocean_state
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: compute_rank

      ! MOM6 frhatu — h-weighted BT-corrector Δu distribution.
      ocean_state%dyn%bt_work%bt_correction_h_weighted = cfg%ocean%bt%correction_h_weighted
      if (compute_rank == 0 .and. cfg%ocean%bt%correction_h_weighted) then
         call logger%info("BT correction:    h-weighted Δu distribution "// &
                          "(MOM6 frhatu pattern)")
      end if

      ! MOM6 BT_force / eta_PF split: the depth-mean baroclinic PGF forces
      ! the barotropic substep (default on).
      ocean_state%dyn%bt_work%bt_bc_pgf_forcing = cfg%ocean%bt%bc_pgf_forcing
      if (compute_rank == 0 .and. .not. cfg%ocean%bt%bc_pgf_forcing) then
         call logger%warning("&ocean_bt_nml bc_pgf_forcing = .false.: LEGACY split — the "// &
                             "barotropic mode does not feel the depth-mean baroclinic "// &
                             "pressure gradient (no JEBAR / bottom-pressure torque)")
      end if

      ! MOM6 btstep_layer_accel — per-layer bc-PGF retro-correction (needs FV_MOM6).
      ocean_state%dyn%bt_work%bt_correction_bc_pgf = cfg%ocean%bt%correction_bc_pgf
      if (compute_rank == 0 .and. cfg%ocean%bt%correction_bc_pgf) then
         call logger%info("BT correction:    bc-PGF retro-correction ON "// &
                          "(MOM6 btstep_layer_accel; requires fv_mom6 PGF)")
      end if

      ! MOM6 bt_rem_u — multiplicative drag damping inside the BT substep.
      ocean_state%dyn%bt_work%bt_substep_drag = cfg%ocean%bt%substep_drag
      if (compute_rank == 0 .and. cfg%ocean%bt%substep_drag) then
         call logger%info("BT substep:       drag damping ON "// &
                          "(MOM6 bt_rem_u: r·HBBL/(Htot+r·HBBL·dt_bt))")
      end if

      ! MOM6 planetary-only fast loop — drop live ζ_bt/∇KE from the substeps
      ! (they stay frozen inside F_bt_*_fast); Cor_ref reduces to f·v̄.
      ocean_state%dyn%bt_work%substep_zeta_ke = cfg%ocean%bt%substep_zeta_ke
      if (compute_rank == 0 .and. .not. cfg%ocean%bt%substep_zeta_ke) then
         call logger%info("BT substep:       planetary Coriolis only "// &
                          "(substep_zeta_ke=.false.; MOM6 q=f/D parity — "// &
                          "live ζ_bt/∇KE off, frozen copies stay in F_bt)")
      end if

      ! MOM6 frhatu·visc_rem joint weight (no-op unless h-weighted is on).
      ! visc_rem is produced every split-path stage by vdiff_apply_momentum
      ! from the momentum tridiagonal (rdb_ocean_vdiff.F90); it is inert
      ! (≡ 1) unless &ocean_vdiff_nml implicit_drag is also on (validate_config
      ! warns in that case).
      ocean_state%dyn%bt_work%bt_correction_visc_rem = cfg%ocean%bt%correction_visc_rem
      if (compute_rank == 0 .and. cfg%ocean%bt%correction_visc_rem) then
         call logger%info("BT correction:    h·visc_rem joint weight ON "// &
                          "(MOM6 frhatu·visc_rem; visc_rem produced by vdiff)")
      end if
      ! MOM6 wt_u parity for the BT FORCING assembly (PGF_BUG.md §9):
      ! friction-damped layers stop forcing the fast loop.  Requires
      ! correction_visc_rem (validated at configure).
      ocean_state%dyn%bt_work%bt_forcing_visc_rem = cfg%ocean%bt%forcing_visc_rem
      if (compute_rank == 0 .and. cfg%ocean%bt%forcing_visc_rem) then
         call logger%info("BT forcing:       h·visc_rem weight ON (MOM6 wt_u)")
      end if
      ! MOM6 dt·visc_rem·accel parity: the per-layer slow applies are
      ! attenuated by the viscous remnant (requires correction_visc_rem,
      ! validated at configure).
      ocean_state%dyn%accel_visc_rem = cfg%ocean%vdiff%accel_visc_rem
      if (compute_rank == 0 .and. cfg%ocean%vdiff%accel_visc_rem) then
         call logger%info("Velocity applies: visc_rem-attenuated ON "// &
                          "(MOM6 u = u0 + dt*visc_rem*accel)")
      end if
      ! Eagerly allocate the accel_visc_rem stage-entry velocity snapshots
      ! here (ocean-state convention: slots allocated at setup, mapped once
      ! at enter_data — no mid-run lazy alloc), but ONLY when the knob is on
      ! — on a large grid these two face arrays are ~GB, so the default-off
      ! path must not pay for them.  Runs BEFORE ocean_state_enter_data, so
      ! ocean_dyn_enter_data_impl copyin sees them allocated (host source=0).
      if (ocean_state%dyn%accel_visc_rem) then
         allocate (ocean_state%dyn%avr_u0(grid%nx_total + 1, grid%ny_total, &
                                          ocean_state%multilayer%nz_ml), source=0.0_wp)
         allocate (ocean_state%dyn%avr_v0(grid%nx_total, grid%ny_total + 1, &
                                          ocean_state%multilayer%nz_ml), source=0.0_wp)
      end if
      ! SPEC S3: outer split-explicit time-scheme selector.  The enum
      ! registration constrains the string (and rejects both retired
      ! spellings of this scheme, "mom6_pc" and "split_rk2", by name);
      ! the default is the MOM6 predictor-corrector, `pred_corr`.
      if (trim(cfg%ocean%bt%split_scheme) == "pred_corr") then
         ocean_state%dyn%split_scheme = SPLIT_SCHEME_PRED_CORR
      else
         ocean_state%dyn%split_scheme = SPLIT_SCHEME_SSP_RK2
      end if
      ocean_state%dyn%pc_be = cfg%ocean%bt%pc_be
      ! Print the RESOLVED scheme for BOTH branches, not just the non-default
      ! one.  A banner that only speaks up for one value cannot answer "which
      ! scheme did this namelist actually run", which is the only reliable way
      ! to enumerate what a default flip touches -- grepping the namelists
      ! cannot (a key may be absent, commented out, or set in an included
      ! group).
      if (compute_rank == 0) then
         if (ocean_state%dyn%split_scheme == SPLIT_SCHEME_PRED_CORR) then
            call logger%info("Split scheme:     pred_corr (MOM6 predictor-corrector; "// &
                             "tendencies on u_av/h_av; BE = "// &
                             to_string(cfg%ocean%bt%pc_be)//")")
         else
            call logger%info("Split scheme:     ssp_rk2 (EXPERIMENTAL; two "// &
                             "identical stages, SSP average -- grows internal "// &
                             "gravity waves from a stratified rest state)")
         end if
      end if
      ! SPEC S2b: γ-weighted continuity transport-matching inversion.
      ocean_state%dyn%bt_work%bt_renorm_visc_rem = cfg%ocean%bt%renorm_visc_rem
      if (compute_rank == 0 .and. cfg%ocean%bt%renorm_visc_rem) then
         call logger%info("BT renormaliser:  gamma-weighted du + u_cor (MOM6 "// &
                          "continuity inversion parity)")
      end if
      ! PR-8: `correction_visc_rem` couples against `visc_rem_u/v`, which are
      ! only ever filled away from their `source=1.0` default by the implicit
      ! bottom-drag diagonal fold (`&ocean_vdiff_nml implicit_drag`) breaking
      ! vdiff's per-row unit-sum normalisation — without it the knob is
      ! exactly gamma==1 (a x1.0). Warn rather than abort (tidal_mixing
      ! e_uniform=0 precedent): PR-19 (visc_rem) is the named owner that
      ! fills the arrays for other configurations.
      if (compute_rank == 0 .and. cfg%ocean%bt%correction_visc_rem .and. &
          .not. cfg%ocean%vdiff%implicit_drag) then
         call logger%warning("&ocean_bt_nml correction_visc_rem=.true. has no effect "// &
                             "without &ocean_vdiff_nml implicit_drag=.true. "// &
                             "(visc_rem_u/v stay at their source=1.0 default; PR-19 owns filling them)")
      end if

      ! BT_cont_type flux-bounded continuity: allocate the per-face coefficient
      ! packs so enter_data attaches them; default off leaves them unallocated.
      ocean_state%dyn%bt_work%use_bt_cont_type = cfg%ocean%bt%use_cont_type
      if (cfg%ocean%bt%use_cont_type) then
         block
            integer :: nx_w, ny_w
            nx_w = grid%nx_total
            ny_w = grid%ny_total
            allocate (ocean_state%dyn%bt_work%BTCL_u(nx_w + 1, ny_w))
            allocate (ocean_state%dyn%bt_work%BTCL_v(nx_w, ny_w + 1))
         end block
      end if
      ! Upstream-PPM h_face for the BT chain — allocate the per-face slots.
      ocean_state%dyn%bt_work%use_upstream_h_face = cfg%ocean%bt%upstream_h_face
      if (cfg%ocean%bt%upstream_h_face) then
         block
            integer :: nx_w, ny_w
            nx_w = grid%nx_total
            ny_w = grid%ny_total
            allocate (ocean_state%dyn%bt_work%h_face_up_x(nx_w + 1, ny_w), source=0.0_wp)
            allocate (ocean_state%dyn%bt_work%h_face_up_y(nx_w, ny_w + 1), source=0.0_wp)
         end block
      end if
      if (compute_rank == 0) then
         if (ocean_state%pressure_force%variant == OPGF_VARIANT_GPRIME) then
            call logger%info("PGF form:         gprime  gfs="// &
                             to_string(cfg%ocean%pgf%gprime_gfs)//" m/s²  gint="// &
                             to_string(cfg%ocean%pgf%gprime_gint)//" m/s²")
         else
            call logger%info("PGF form:         "//trim(cfg%ocean%pgf%form))
         end if
         if (cfg%ocean%bt%bebt > 0.0_wp) then
            call logger%info("BEBT projection:  bebt="// &
                             to_string(cfg%ocean%bt%bebt)// &
                             " (BT η flux uses (1+bebt)·u^n − bebt·u^{n-1})")
         end if
         if (cfg%ocean%bt%use_cont_type) then
            call logger%info("BT_cont_type:     ON  (flux-bounded BT continuity)")
            if (cfg%ocean%bt%cont_corr_bounds) then
               call logger%info("                  BT_CONT_CORR_BOUNDS=.true. (η-corr bound from BT_cont)")
            end if
         end if
         if (cfg%ocean%bt%upstream_h_face) then
            call logger%info("BT upstream-h:    ON  (BT chain uses upstream-PPM h_face)")
         end if
      end if
   end subroutine configure_ocean_bt

   subroutine configure_ocean_bt_split(cfg, ocean_state, grid, compute_rank)
      !! Auto-derive the barotropic substep count n_inner from the external
      !! gravity-wave CFL (MOM6 set_dtbt) when requested, then latch the
      !! mode-split reference column depth bt_H_ref from the seeded bathymetry.
      !! `cfg` is intent(inout) because auto_n_inner writes cfg%ocean%bt%n_inner.
      type(config_t), intent(inout) :: cfg
      type(ocean_state_t), intent(inout) :: ocean_state
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: compute_rank

      ! Auto-derive n_inner from the external gravity-wave CFL (MOM6 set_dtbt),
      ! evaluated PER WET CELL: the local depth with the local 2-D CFL length
      ! `l = 1/sqrt(1/dx^2+1/dy^2)` (the cross-direction term matters — on
      ! square cells it is the sqrt(2) the legacy 1-D estimate omitted), over
      ! ocean only.  See `bt_cfl_dt_wet` for why the old "deepest anywhere x
      ! smallest anywhere, land included" combination was wrong on a real
      ! global grid, and why this reproduces it bit-for-bit wherever the two
      ! extremes coincide on a wet cell.
      !
      ! `multilayer%wet_mask`, not `metrics%wet_T`: the metric mask is built
      ! later, by `configure_ocean_land_mask`; its interior IS this array.
      if (cfg%ocean%bt%auto_n_inner) then
         block
            integer :: ng_, n_inner_derived, n_wet_local
            real(wp) :: dt_bt_local, dt_bt_safe, h_at, l_at
            ng_ = grid%nghost
            call bt_cfl_dt_wet(size(ocean_state%barotropic%b, 1), &
                               size(ocean_state%barotropic%b, 2), &
                               ng_ + 1, ng_ + grid%nx_phys, &
                               ng_ + 1, ng_ + grid%ny_phys, &
                               ocean_state%barotropic%b, &
                               ocean_state%multilayer%wet_mask, &
                               ocean_state%metrics%dxT, ocean_state%metrics%dyT, &
                               cfg%ocean%bt%cfl_bt_safety, &
                               dt_bt_local, h_at, l_at, n_wet_local)
            ! Reduce to the GLOBAL per-point minimum so every rank derives the
            ! SAME n_inner.  Different n_inner per rank means a different
            ! number of barotropic substeps, which desyncs the per-substep
            ! grouped halo exchanges (Isend/Irecv pairing crosses between
            ! substeps -> MPI_ERR_TRUNCATE at the first N/S exchange; seen on
            ! a spherical py>1 split).  A rank with no wet cell contributes
            ! `huge`, the identity of min.  Single rank: identity stub.
            call halo_allreduce_min(dt_bt_local, dt_bt_safe)
            n_inner_derived = bt_auto_n_inner_from_dt(cfg%dt_fixed, dt_bt_safe)
            if (compute_rank == 0) then
               if (dt_bt_safe >= huge(1.0_wp)) then
                  call logger%warning("Auto n_inner: no wet cell anywhere — "// &
                                      "n_inner = 1 (was "// &
                                      to_string(cfg%ocean%bt%n_inner)//")")
               else if (dt_bt_local == dt_bt_safe) then
                  call logger%info("Auto n_inner (per wet cell): limiting cell H = "// &
                                   to_string(h_at)//" m, l_cfl = "// &
                                   to_string(l_at)//" m, c_ext = "// &
                                   to_string(sqrt(GRAVITY*max(h_at, 1.0_wp)))// &
                                   " m/s, dt_bt = "//to_string(dt_bt_safe)// &
                                   " s → n_inner = "//to_string(n_inner_derived)// &
                                   " (was "//to_string(cfg%ocean%bt%n_inner)//")")
               else
                  call logger%info("Auto n_inner (per wet cell): dt_bt = "// &
                                   to_string(dt_bt_safe)//" s (limited on another "// &
                                   "rank) → n_inner = "//to_string(n_inner_derived)// &
                                   " (was "//to_string(cfg%ocean%bt%n_inner)//")")
               end if
            end if
            cfg%ocean%bt%n_inner = n_inner_derived
         end block
      end if

      ! Mode-split contract: bt_H_ref is the reference WATER-COLUMN
      ! thickness.  Without an ice shelf that is the seeded bathymetry
      ! (Σh_layer = b); under one it is `b − z_draft`, which is what makes
      ! `bt_eta = Σh_layer − bt_H_ref` the deviation from the LOADED
      ! equilibrium (zero at rest) rather than a permanent −z_draft.  That
      ! is Losch (2008) §2.1's own convention, and it is why every
      ! consumer of `D = bt_H_ref + bt_eta` — the BT continuity face
      ! thickness, the Chapman phase speed, the ALE `remap_h_ref` — needs
      ! no cavity branch of its own.  Knob off ⇒ `z_draft ≡ 0` ⇒ the
      ! literal `bt_H_ref = b` this always was.
      !
      ! The snapshot is taken from the UNWRAPPED `b`/`z_draft`; the engine
      ! re-wraps + halo-exchanges `bt_H_ref` right after, alongside both
      ! of its sources.
      !
      ! NOTE on `auto_n_inner` above: it derives each wet cell's external
      ! gravity-wave speed from `b`, the BED depth, not from `b − z_draft`.
      ! Under a shelf that OVERESTIMATES `c_ext` and so buys more
      ! barotropic substeps than the CFL needs — conservative, never
      ! unstable, and at a calving front (where the draft is 0) it is
      ! exactly right.  Left as-is deliberately: tightening it would
      ! change `n_inner` for cavity runs only, which is a separate,
      ! answer-changing decision.
      if (cfg%ocean%bt%n_inner >= 1) then
         if (ocean_state%metrics%use_cavity) then
            ! `cavity_datum_impl`, not the raw `b - z_draft`: a GROUNDED
            ! column has no water column, and its datum is 0 rather than
            ! a negative thickness.  See that routine for why — in short,
            ! it is land, nothing downstream reads its datum, and the
            ! negative value put a phantom few-hundred-metre `bt_eta` on
            ! it and turned the ALE land target into a cancellation.
            call cavity_datum_impl(ocean_state%dyn%bt_work%bt_H_ref, &
                                   ocean_state%barotropic%b, &
                                   ocean_state%metrics%z_draft, &
                                   cfg%ocean%cavity_dyn%h_min_cavity, &
                                   size(ocean_state%barotropic%b, 1), &
                                   size(ocean_state%barotropic%b, 2))
         else
            ocean_state%dyn%bt_work%bt_H_ref = ocean_state%barotropic%b
         end if
      end if
   end subroutine configure_ocean_bt_split

   subroutine configure_ocean_cavity(cfg, ocean_state, grid, compute_rank, ierr)
      !! Build the static ice-shelf cavity LOAD field
      !! `metrics%p_ice_ref = (rho_ref*GRAVITY)*z_draft` (Pa), ASSEMBLE it
      !! into the top-of-column pressure `multilayer_state_t%p_top`, and
      !! assert the counted-once datum invariant.
      !!
      !! The GEOMETRY (`z_draft`, `cover_frac`) is filled much earlier, in
      !! `ocean_state_seed_from_cfg`, because the wet mask and the layer
      !! split are seeded from `b − z_draft`.  What is left for configure
      !! is the part that needs the PGF's reference density, which
      !! `configure_ocean_pgf` / `configure_ocean_reference_density` only
      !! settle later — hence this runs after them and before
      !! `ocean_state_enter_data`, like every other static field the
      !! device map has to capture.
      !!
      !! `rho_ref*GRAVITY` is formed as ONE product, the same one the
      !! FV_MOM6 Pass-1 surface BC forms, so that
      !! `pa(nz+1) = rho_ref*g*(−z_draft) + p_ice_ref` cancels to bit-zero
      !! at rest (exactly when the toolchain rounds the product before the
      !! add; under FMA contraction, to the rounding of `p_ice_ref`).
      !!
      !! ### The partition (P5.2)
      !!
      !! ```
      !! ms%p_top = metrics%p_ice_ref  +  sf%p_surf
      !!            (static ice load)     (atmospheric / anomaly load)
      !! ```
      !!
      !! `p_ice_ref` goes HERE and to the datum (`bt_H_ref = b − z_draft`)
      !! and NOWHERE else — in particular it is never added into
      !! `sf%p_surf`, because `eta_ib = −p_surf/(ρ₀ g_bt)` is built from
      !! the assembled total and the datum already carries exactly this
      !! much.  Only the load ANOMALY reaches the `eta_forcing` seam, and
      !! for the Boussinesq-isostatic default that anomaly IS `sf%p_surf`
      !! (a cavity with no atmospheric load sends the seam nothing at
      !! all).  See `src/core/ocean/README.md`'s `p_top` / `eta_forcing`
      !! seam contracts.
      !!
      !! This seed is the FINAL value for a cavity without the psurf seam
      !! (the draft is static, so there is nothing to refresh); with psurf
      !! enabled, `ocean_dyn_step_split` rebuilds the same sum once per
      !! outer step from the live `sf%p_surf`.  Host-side plain `do`
      !! loops, before `ocean_state_enter_data` maps the result.
      type(config_t), intent(in) :: cfg
      type(ocean_state_t), intent(inout) :: ocean_state
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: compute_rank
      integer, intent(out), optional :: ierr
         !! Non-zero on a cavity configuration conflict when present;
         !! absent behaves as today (`error stop`).

      integer :: nx, ny, i, j
      real(wp) :: resid, datum_tol

      if (.not. ocean_state%metrics%use_cavity) then
         if (present(ierr)) ierr = OCEAN_STATUS_OK
         return
      end if

      nx = size(ocean_state%metrics%z_draft, 1)
      ny = size(ocean_state%metrics%z_draft, 2)

      call cavity_fill_p_ice_ref(ocean_state%metrics%p_ice_ref, &
                                 ocean_state%metrics%z_draft, &
                                 ocean_state%pressure_force%rho_ref*GRAVITY, nx, ny)

      ! P5.2 — THE LOAD MUST HAVE A CONSUMER WHEN IT HAS A GRADIENT.
      !
      ! `&ocean_pgf_nml p_top_in_bc` is the only route by which a cavity
      ! load reaches the pressure stack, and a cavity whose draft VARIES
      ! is not in hydrostatic balance without it: `pa(nz+1)` then carries
      ! the uncancelled `-rho_ref*g*z_draft`, so the stack sits ~5e6 Pa
      ! off its anomaly scale (every `h_neglect` face-divisor leak grows
      ! by the same factor), the UNSPLIT driver — which has no
      ! depth-mean replacement — feels a raw `g*grad(z_draft)` ~ 0.1 m/s^2,
      ! and `&ocean_bt_nml correction_h_weighted` turns the uncancelled
      ! depth-uniform force into a real per-layer shear.  Refused, not
      ! auto-enabled: the namelist should say what the run does.
      !
      ! EXEMPTION, and it is a theorem rather than a courtesy: a draft
      ! that is UNIFORM over the whole array has a load with no gradient,
      ! and a gradient-free `p_top` is bit-identically inert in the top BC
      ! (`test_ocean_pgf_p_top_bc::uniform_p_top_bit_identical`).  Such a
      ! run — the flat-lid datum-equivalence case — is refused by nothing.
      !
      ! Tested on the FILLED array, not on the namelist shape: land
      ! exclusion can zero `z_draft` under a formula that looks uniform.
      ! `validate_config` carries the same refusal on the config-level
      ! predicate so the user meets it before any state is built.
      if (.not. ocean_state%pressure_force%p_top_in_bc) then
         if (maxval(ocean_state%metrics%z_draft) /= &
             minval(ocean_state%metrics%z_draft)) then
            call fail("&ocean_cavity_dyn_nml enable=.true. with a NON-UNIFORM "// &
                      "draft requires &ocean_pgf_nml p_top_in_bc=.true.  The "// &
                      "isostatic load rho_ref*g*z_draft would otherwise never "// &
                      "reach the FV_MOM6 pa(nz+1) surface boundary condition, "// &
                      "leaving the pressure stack ~"// &
                      to_string(maxval(ocean_state%metrics%p_ice_ref))// &
                      " Pa off its anomaly scale and the column out of "// &
                      "hydrostatic balance.  (A UNIFORM draft is exempt: a "// &
                      "load with no gradient is provably inert in the top BC.)", &
                      ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
      end if

      ! P5.2: assemble the top-of-column load.  Rebuilt from scratch (not
      ! `+=` onto whatever the earlier psurf seed left) so the result does
      ! not depend on which of the two seeds ran first, and over the WHOLE
      ! array including ghosts — `p_top` owes no halo exchange of its own
      ! precisely because both of its sources are already ghost-valid.
      if (.not. (size(ocean_state%multilayer%p_top, 1) == nx .and. &
                 size(ocean_state%multilayer%p_top, 2) == ny)) then
         call fail("configure_ocean_cavity: ms%p_top and metrics%p_ice_ref "// &
                   "have different shapes — the top-of-column load cannot be "// &
                   "assembled.", ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if
      do j = 1, ny
         do i = 1, nx
            ocean_state%multilayer%p_top(i, j) = ocean_state%metrics%p_ice_ref(i, j)
         end do
      end do
      if (allocated(ocean_state%surface_flux%p_surf)) then
         do j = 1, ny
            do i = 1, nx
               ocean_state%multilayer%p_top(i, j) = &
                  ocean_state%multilayer%p_top(i, j) + &
                  ocean_state%surface_flux%p_surf(i, j)
            end do
         end do
      end if

      ! ---- The vertical coordinate's rigid-top seam (P6.2) ----
      ! `vcoord%z_top(i,j)` is the geopotential depth of the top of the
      ! WATER column — the ice base under a shelf, `z = 0` elsewhere.
      ! The draft is static, so this is a configure-time copy and the
      ! vcoord slot stays self-contained at run time: the ALE remap
      ! driver never sees `metrics` and the target builder's signature is
      ! unchanged.  Same pattern (and the same reason) as the sponge's
      ! `sp%z_top` in `configure_ocean_sponge`.
      !
      ! Consumed by the `VCOORD_Z_FIXED` branch of `compute_target_h`.
      ! Without a cavity this routine has already returned, so `z_top`
      ! keeps its init-time zero fill and every geometric family
      ! reproduces its pre-cavity arithmetic bit-for-bit.  BEFORE
      ! `ocean_state_enter_data` maps it.
      if (allocated(ocean_state%vcoord%z_top)) then
         if (size(ocean_state%vcoord%z_top, 1) == nx .and. &
             size(ocean_state%vcoord%z_top, 2) == ny) then
            do j = 1, ny
               do i = 1, nx
                  ocean_state%vcoord%z_top(i, j) = ocean_state%metrics%z_draft(i, j)
               end do
            end do
         else
            call fail("configure_ocean_cavity: vcoord%z_top and metrics%z_draft "// &
                      "have different shapes — the vertical coordinate cannot see "// &
                      "the ice base.", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
      end if

      ! The counted-once invariant (I), in metres of reference depth, over
      ! the WET columns:
      !     rho*g*z_draft + (bt_H_ref - b)*rho*g == 0   <=>   bt_H_ref == b - z_draft.
      ! GROUNDED columns are excluded because they are LAND: every face
      ! metric on them is zero, so they carry no barotropic momentum
      ! equation and there is no load on them to count once or twice.
      ! Their datum is deliberately 0 (`cavity_datum_impl`).
      ! Asserted on the common positive factor divided out — scale-free,
      ! and it does not fabricate a product the code never forms.  The
      ! bound is a pure round-off allowance on the ONE subtraction
      ! `b - z_draft`: both operands are O(max depth), so the result
      ! carries at most a few ulp of it.  (It is NOT asserted as bit-zero:
      ! the latch and this check evaluate the same difference in two
      ! places, and an FMA-contracting build is free to round them
      ! differently.)
      if (cfg%ocean%bt%n_inner >= 1) then
         datum_tol = 8.0_wp*epsilon(1.0_wp)* &
                     max(maxval(abs(ocean_state%barotropic%b)), 1.0_wp)
         resid = cavity_datum_residual(ocean_state%dyn%bt_work%bt_H_ref, &
                                       ocean_state%barotropic%b, &
                                       ocean_state%metrics%z_draft, &
                                       cfg%ocean%cavity_dyn%h_min_cavity, nx, ny)
         if (.not. (resid <= datum_tol)) then
            call fail("&ocean_cavity_dyn_nml: the barotropic datum and the ice "// &
                      "draft disagree (max |bt_H_ref - (b - z_draft)| = "// &
                      to_string(resid)//" m > "//to_string(datum_tol)//" m).  The "// &
                      "ice load would then be counted twice, or not at all — check "// &
                      "that z_draft went through the same periodic/fold re-wrap and "// &
                      "halo exchange as b.", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
      end if

      ! ---- Cover mask on the atmospheric forcing (P2c) ----
      ! There is no atmosphere under an ice shelf.  Two static forcing
      ! fields are masked ONCE, here, because this is the first point at
      ! which `cover_frac` exists AND the forcing has been seeded
      ! (`configure_ocean_forcing` runs much earlier, before the cavity
      ! geometry):
      !
      !   * the wind-stress PAIR, masked on every face touching a covered
      !     cell and followed by the `stress_mag` refresh in the same
      !     call — see `ocean_surface_stress_apply_cover`.  Masking the
      !     source rather than the derived views is what also silences
      !     the implicit vdiff stress fold and the MLE front sampler,
      !     which read `ss%tau_x` raw;
      !   * the SCALAR `&ocean_thermo_nml q_heat` / `q_salt` fill, but
      !     ONLY when the component set is off.  With components on those
      !     two are assembler outputs and the mask belongs in
      !     `ocean_surface_flux_assemble` (a second writer here would
      !     break the fill contract); the call below no-ops itself in
      !     that case.
      !
      ! Both are host-side and both run BEFORE `ocean_state_enter_data`,
      ! so the masked values are what the device map captures.  Both are
      ! idempotent.  The time-varying twin of the wind mask lives in
      ! `ocean_seam_refresh_surface_stress` (per data-forcing bracket).
      call ocean_surface_stress_apply_cover(ocean_state%surface_stress, &
                                            ocean_state%metrics%cover_frac)
      call ocean_surface_flux_apply_cover_const(ocean_state%surface_flux, &
                                                ocean_state%metrics%cover_frac)

      if (compute_rank == 0) then
         call logger%info("Ice-shelf cavity: isostatic load p_ice_ref = "// &
                          "rho_ref*g*z_draft, max = "// &
                          to_string(maxval(ocean_state%metrics%p_ice_ref))// &
                          " Pa (rho_ref = "// &
                          to_string(ocean_state%pressure_force%rho_ref)//" kg/m^3)")
         call logger%info("                  assembled into ms%p_top (max = "// &
                          to_string(maxval(ocean_state%multilayer%p_top))// &
                          " Pa = p_ice_ref + sf%p_surf); NOT into sf%p_surf, "// &
                          "which the eta_forcing seam is built from")
         if (ocean_state%pressure_force%p_top_in_bc) then
            call logger%info("                  consumed by the FV_MOM6 pa(nz+1) "// &
                             "surface BC (&ocean_pgf_nml p_top_in_bc)")
         else
            call logger%info("                  the draft is UNIFORM here, so the "// &
                             "load has no gradient and the PGF top BC is provably "// &
                             "inert; &ocean_pgf_nml p_top_in_bc is not required")
         end if
      end if
      if (present(ierr)) ierr = OCEAN_STATUS_OK
   end subroutine configure_ocean_cavity

   subroutine configure_ocean_cavity_melt(cfg, ocean_state, grid, compute_rank, ierr)
      !! Configure the ice-shelf basal-melt slot (`&ocean_cavity_melt_nml`,
      !! P2b): copy the knobs onto the slot's three flat parameter
      !! bundles, resolve the `gamma_s` sentinel, build the per-column
      !! Coriolis array the `hj99` law needs, and CHECK that the
      !! interface pressure the liquidus will read has actually been
      !! loaded.
      !!
      !! Runs immediately after `configure_ocean_cavity` — which is the
      !! SOLE producer of `ms%p_top` — and before
      !! `ocean_state_enter_data`, like every other static fill the
      !! device map has to capture.
      !!
      !! **`ms%p_top`: this routine is a CONSUMER, not a writer.**
      !! `ms%p_top` is THE interface pressure
      !! (`src/core/ocean/README.md`, the `p_top` seam contract), and the
      !! melt liquidus is its THIRD consumer alongside the FV_MOM6
      !! `pa(nz+1)` surface boundary condition and the in-situ EOS.  It
      !! has exactly ONE producer, `configure_ocean_cavity`, which
      !! assembles `ms%p_top = metrics%p_ice_ref + sf%p_surf` (P5.2) —
      !! rebuilt per outer step in `ocean_dyn_step_split` when the psurf
      !! seam makes `sf%p_surf` live.  A second writer here would be a
      !! silent clobber, so there is none: what this routine does instead
      !! is ASSERT that every ice-covered column carries at least its own
      !! isostatic load, which catches a producer that was skipped or
      !! reordered rather than trusting the call order.
      use rdb_ocean_cavity_melt, only: parse_cavity_exchange_law, parse_cavity_ice_mode, &
                                       CAVITY_LAW_HJ99, parse_cavity_freshwater, &
                                       parse_cavity_volume_comp, CAVITY_FW_MASS, &
                                       CAVITY_VC_UNIFORM_OPEN
      type(config_t), intent(in) :: cfg
      type(ocean_state_t), intent(inout) :: ocean_state
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: compute_rank
      integer, intent(out), optional :: ierr
         !! Non-zero on a melt configuration conflict when present;
         !! absent behaves as today (`error stop`).
      integer :: nx, ny, n_zero_f, n_unloaded
      real(wp) :: gamma_s_eff

      if (.not. ocean_state%cavity_flux%enable) then
         if (present(ierr)) ierr = OCEAN_STATUS_OK
         return
      end if

      nx = size(ocean_state%cavity_flux%f_cor, 1)
      ny = size(ocean_state%cavity_flux%f_cor, 2)

      ! ---- knobs -> slot ----
      ocean_state%cavity_flux%far_field_depth = cfg%ocean%cavity_melt%far_field_depth
      ocean_state%cavity_flux%cdrag_top = cfg%ocean%cavity_melt%cdrag_top
      ocean_state%cavity_flux%u_tide = cfg%ocean%cavity_melt%u_tide
      ocean_state%cavity_flux%ustar_min = cfg%ocean%cavity_melt%ustar_min
      ocean_state%cavity_flux%s_ice = cfg%ocean%cavity_melt%s_ice
      ocean_state%cavity_flux%freshwater = &
         parse_cavity_freshwater(cfg%ocean%cavity_melt%freshwater)
      ocean_state%cavity_flux%volume_comp = &
         parse_cavity_volume_comp(cfg%ocean%cavity_melt%volume_compensation)

      ! THE Boussinesq reference density, taken from the surface-flux
      ! slot rather than re-read from `cfg`, because the real-mass path's
      ! salt correction has to be the EXACT negation of the increment
      ! `apply_surface_src_2d_impl` stamped with `dt/sf%rho0`.  One
      ! source, no second literal, and `configure_ocean_reference_density`
      ! has already run (it seeds `sf%rho0` from `&ocean_ic_nml rho_0`).
      ocean_state%cavity_flux%rho0 = ocean_state%surface_flux%rho0
      if (ocean_state%cavity_flux%rho0 <= 0.0_wp) then
         call fail("&ocean_cavity_melt_nml: the Boussinesq reference density "// &
                   "reaching the melt slot is non-positive, so the mass -> volume "// &
                   "conversion m/rho_0 has no value.  configure_ocean_reference_"// &
                   "density must run before configure_ocean_cavity_melt.", ierr, &
                   OCEAN_STATUS_ERR_SETUP)
         return
      end if

      ! `gamma_s` carries the repo's negative "unset" sentinel; resolve
      ! it to the ISOMIP+ `gamma_t/35` exactly once, here, so the value
      ! the kernel is handed is the value the log prints.
      gamma_s_eff = cavity_resolve_gamma_s(cfg%ocean%cavity_melt%gamma_s, &
                                           cfg%ocean%cavity_melt%gamma_t)
      ocean_state%cavity_flux%par%law = &
         parse_cavity_exchange_law(cfg%ocean%cavity_melt%exchange_law)
      ocean_state%cavity_flux%par%gamma_t_coeff = cfg%ocean%cavity_melt%gamma_t
      ocean_state%cavity_flux%par%gamma_s_coeff = gamma_s_eff
      ocean_state%cavity_flux%ice%mode = &
         parse_cavity_ice_mode(cfg%ocean%cavity_melt%ice_conduction)
      ocean_state%cavity_flux%ice%T_ice = cfg%ocean%cavity_melt%t_ice

      ! `cavity_flux%const` is deliberately LEFT at its ISOMIP+ defaults
      ! (Asay-Davis et al. (2016) Table 4): `rho_w`, `c_w`, `alpha_T`,
      ! `beta_S` there are the melt law's own calibrated constants, not
      ! copies of the Boussinesq `rho_0` — the same "out of scope on
      ! purpose" category the README's reference-density table already
      ! lists `RHO_WATER` and `&ocean_ice_nml rho_ocean` under.
      ! Overriding them would silently break ISOMIP+ comparability,
      ! which is the reason this path exists.

      ! ---- per-column Coriolis for the hj99 law ----
      call fill_coriolis_centre(cfg, ocean_state%metrics, grid, ocean_state%cavity_flux%f_cor)
      if (ocean_state%cavity_flux%par%law == CAVITY_LAW_HJ99) then
         n_zero_f = cavity_count_zero_f(ocean_state%cavity_flux%f_cor, &
                                        ocean_state%metrics%cover_frac, &
                                        ocean_state%multilayer%wet_mask, nx, ny)
         if (n_zero_f > 0) then
            call fail("&ocean_cavity_melt_nml exchange_law='hj99': "// &
                      to_string(n_zero_f)//" ice-covered wet column(s) sit at "// &
                      "f = 0.  Holland & Jenkins (1999) eq. (15) takes "// &
                      "ln(.../|f| h_nu) and eq. (18) divides by f*L_O, so the law "// &
                      "has no value there and the kernel would return "// &
                      "CAVITY_MELT_NO_CORIOLIS for every one of them.  Use "// &
                      "exchange_law='const_gamma', or move the cavity off the "// &
                      "f = 0 line.", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
      end if

      ! ---- the interface pressure: ASSERT, never write (see the docstring) ----
      ! `p_surf >= 0` by contract, so a correctly assembled `p_top` is
      ! `>= p_ice_ref` on every cell.  A zero `p_top` under a loaded
      ! draft means the producer never ran (or ran before the load was
      ! built), which would melt against a surface-pressure liquidus and
      ! silently delete the entire ice pump.
      n_unloaded = cavity_count_unloaded_p_top(ocean_state%multilayer%p_top, &
                                               ocean_state%metrics%p_ice_ref, nx, ny)
      if (n_unloaded > 0) then
         call fail("&ocean_cavity_melt_nml: "//to_string(n_unloaded)//" column(s) "// &
                   "have ms%p_top < metrics%p_ice_ref, so the top-of-column load was "// &
                   "never assembled.  The melt liquidus is a CONSUMER of ms%p_top; "// &
                   "its sole producer is configure_ocean_cavity "// &
                   "(p_top = p_ice_ref + sf%p_surf), which must run before this "// &
                   "routine.", ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if

      if (compute_rank == 0) then
         call logger%info("Ice-shelf basal melt: exchange_law='"// &
                          trim(adjustl(cfg%ocean%cavity_melt%exchange_law))// &
                          "', Gamma_T = "//to_string(cfg%ocean%cavity_melt%gamma_t)// &
                          ", Gamma_S = "//to_string(gamma_s_eff)// &
                          ", ice='"// &
                          trim(adjustl(cfg%ocean%cavity_melt%ice_conduction))//"'")
         call logger%info("                      far_field_depth = "// &
                          to_string(cfg%ocean%cavity_melt%far_field_depth)// &
                          " m (METRES below the ice base, not 'layer nz'), "// &
                          "C_d,top = "//to_string(cfg%ocean%cavity_melt%cdrag_top)// &
                          ", u_tide = "//to_string(cfg%ocean%cavity_melt%u_tide)//" m/s")
         call logger%info("                      liquidus evaluated at ms%p_top "// &
                          "(max = "// &
                          to_string(maxval(ocean_state%multilayer%p_top))// &
                          " Pa), assembled by configure_ocean_cavity as "// &
                          "p_ice_ref + sf%p_surf")
         if (ocean_state%cavity_flux%freshwater == CAVITY_FW_MASS) then
            call logger%info("                      freshwater = 'mass': the "// &
                             "meltwater is a REAL Boussinesq volume on the top "// &
                             "layer, dh = m*dt/rho_0 with rho_0 = "// &
                             to_string(ocean_state%cavity_flux%rho0)//" kg/m3.  The "// &
                             "virtual salt flux stays in Q_salt as the KPP/EPBL B_0 "// &
                             "buoyancy forcing and is removed again from the tracer "// &
                             "in the same stage.")
            if (ocean_state%cavity_flux%volume_comp == CAVITY_VC_UNIFORM_OPEN) then
               call logger%info("                      volume_compensation = "// &
                                "'uniform_open_ocean': the domain-integrated melt "// &
                                "volume is removed again each thermo step over the "// &
                                "uncovered wet cells, carrying their own T and S, "// &
                                "and is tracked as a sink in all three budgets.")
            else
               call logger%warning("Ice-shelf basal melt: freshwater='mass' with "// &
                                   "volume_compensation='none' — a CLOSED domain "// &
                                   "gains the melt volume and its sea level rises "// &
                                   "(for an ISOMIP+ Ocean0 box that is metres per "// &
                                   "year).  Correct for a short run or an open "// &
                                   "boundary; otherwise set "// &
                                   "volume_compensation='uniform_open_ocean'.")
            end if
         else
            call logger%warning("Ice-shelf basal melt: the meltwater is a VIRTUAL SALT "// &
                                "FLUX — it carries no mass, so it adds no volume and no "// &
                                "direct buoyancy.  That is a first-order limitation for "// &
                                "cavity circulation; set &ocean_cavity_melt_nml "// &
                                "freshwater='mass' for the real volume source.")
         end if
      end if
      if (present(ierr)) ierr = OCEAN_STATUS_OK
   end subroutine configure_ocean_cavity_melt

   subroutine configure_ocean_top_drag(cfg, ocean_state, grid, compute_rank)
      !! Ice-shelf TOP drag (`&ocean_tdrag_nml`, Phase 4a): copy the
      !! variant + coefficients onto the slot, project the static
      !! cell-centred `metrics%cover_frac` onto the velocity FACES, and
      !! enforce the ONE-`C_d` rule against the melt slot.
      !!
      !! Runs AFTER `configure_ocean_cavity_melt` (which is where the melt
      !! slot's `cdrag_top` is seeded, and which this routine then
      !! overwrites from `&ocean_tdrag_nml cd` — `validate_config` has
      !! already refused a disagreement, so the assignment is structural,
      !! not a silent override) and AFTER the `cover_frac` halo exchange,
      !! and BEFORE `ocean_state_enter_data` so the host-filled face masks
      !! reach the device with the `copyin` map.
      !!
      !! Disabled ⇒ nothing is written, the slot keeps its placeholder
      !! arrays, and the run is bit-identical.
      type(config_t), intent(in) :: cfg
      type(ocean_state_t), intent(inout) :: ocean_state
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: compute_rank

      integer :: nx, ny

      if (.not. cfg%ocean%tdrag%enable) return

      nx = grid%nx_total
      ny = grid%ny_total

      ocean_state%tdrag%variant = parse_tdrag_variant(cfg%ocean%tdrag%form)
      ocean_state%tdrag%c_drag = cfg%ocean%tdrag%cd
      ocean_state%tdrag%r_linear = cfg%ocean%tdrag%r
      ocean_state%tdrag%htbl = cfg%ocean%tdrag%htbl
      ocean_state%tdrag%drag_bg_vel = cfg%ocean%tdrag%bg_vel
      ocean_state%tdrag%tbl_thick_min = cfg%ocean%tdrag%tbl_thick_min
      ocean_state%tdrag%implicit = cfg%ocean%tdrag%implicit
      ! `&ocean_vdiff_nml implicit_top_drag` is the OTHER implicit form:
      ! the fold into the vdiff `k = nz` diagonal.  Latched on both slots
      ! — the top-drag kernel needs it to fill `lambda_top_u/v`, and the
      ! driver reads `td%implicit_fold` to skip the explicit apply.
      ! `validate_config` has already refused the double-count
      ! (`implicit_top_drag` × `&ocean_tdrag_nml implicit`) and the
      ! distributed case (`htbl > 0`).
      ocean_state%tdrag%implicit_fold = cfg%ocean%vdiff%implicit_top_drag
      ocean_state%vdiff%implicit_top_drag = cfg%ocean%vdiff%implicit_top_drag
      ! `rho0` is only ever used to scale the `stress_top` diagnostic into
      ! N/m^2; no dynamics reads it.  The one rho0 of record is
      ! `&ocean_ic_nml rho_0` -> `eos%rho0`, already resolved by
      ! `configure_ocean_reference_density` well before this point.
      ocean_state%tdrag%rho0 = ocean_state%eos%rho0

      ! Static face cover, OR of the two abutting cells (see
      ! `rdb_ocean_top_drag`'s module docstring for why OR at a calving
      ! front), plus the cell-centred copy the `stress_top` diagnostic
      ! reads.  `cover_frac` is `&ocean_cavity_dyn_nml` geometry: filled
      ! in the IC seed, halo-exchanged by the engine, never touched again.
      call top_drag_fill_face_cover_impl(ocean_state%tdrag%cover_u, &
                                         ocean_state%tdrag%cover_v, &
                                         ocean_state%metrics%cover_frac, nx, ny)
      ocean_state%tdrag%cover_t(:, :) = ocean_state%metrics%cover_frac(:, :)

      ! ONE drag coefficient for momentum and melt.  MOM6 carries two
      ! independent top-drag coefficients; we deliberately do not — see
      ! the agreement rule in `validate_config`, which has already failed
      ! loud if the user set two different values.
      if (cfg%ocean%cavity_melt%enable .and. &
          ocean_state%tdrag%variant == TDRAG_QUADRATIC) then
         ocean_state%cavity_flux%cdrag_top = cfg%ocean%tdrag%cd
      end if

      if (compute_rank == 0) then
         call logger%info("Ice-shelf top drag: "//trim(cfg%ocean%tdrag%form)// &
                          "  C_d="//to_string(cfg%ocean%tdrag%cd)// &
                          "  r="//to_string(cfg%ocean%tdrag%r)//" 1/s")
         if (cfg%ocean%tdrag%htbl > 0.0_wp) then
            call logger%info("  top BL:         HTBL="//to_string(cfg%ocean%tdrag%htbl)// &
                             " m  bg_vel="//to_string(cfg%ocean%tdrag%bg_vel)// &
                             " m/s  thick_min="// &
                             to_string(cfg%ocean%tdrag%tbl_thick_min)//" m")
         else
            call logger%info("  top BL:         layer-nz only (HTBL=0)")
         end if
         call logger%info("  covered faces:  u "// &
                          to_string(int(sum(ocean_state%tdrag%cover_u)))//", v "// &
                          to_string(int(sum(ocean_state%tdrag%cover_v)))// &
                          " (OR of the two abutting cells)")
         if (cfg%ocean%cavity_melt%enable) then
            call logger%info("  melt u* C_d:    taken from &ocean_tdrag_nml cd "// &
                             "(one coefficient for momentum and melt)")
         end if
      end if
   end subroutine configure_ocean_top_drag

   pure function cavity_resolve_gamma_s(gamma_s, gamma_t) result(gamma_s_eff)
      !! Resolve the `&ocean_cavity_melt_nml gamma_s` "unset" sentinel to
      !! the ISOMIP+ default `gamma_t/CAVITY_GAMMA_RATIO_ISOMIP`
      !! (Asay-Davis et al. (2016) Table 4 p. 2483; the 35 is Jenkins,
      !! Nicholls & Corr (2010) p. 2309).
      !!
      !! A NEGATIVE `gamma_s` means "not set by the user"; zero and
      !! positive values are taken literally — zero is then refused as a
      !! range error by `validate_config` rather than silently
      !! re-triggering the default, because the three-equation form
      !! divides by `gamma_s`.
      real(wp), intent(in) :: gamma_s
         !! The raw knob; negative = unset sentinel.
      real(wp), intent(in) :: gamma_t
         !! Heat-transfer coefficient the default is a fraction of.
      real(wp) :: gamma_s_eff
      gamma_s_eff = merge(gamma_t/CAVITY_GAMMA_RATIO_ISOMIP, gamma_s, gamma_s < 0.0_wp)
   end function cavity_resolve_gamma_s

   pure function cavity_count_zero_f(f_cor, cover_frac, wet_mask, nx, ny) result(n_zero)
      !! Count ice-covered WET columns sitting at exactly `f = 0` — the
      !! configure-time domain check for `exchange_law = "hj99"`.
      !!
      !! Exactly zero, not "small": the law's failure at `f = 0` is a
      !! division and a logarithm, not a loss of accuracy, and a small
      !! `|f|` is a legitimate (if strongly suppressed) answer.
      integer, intent(in) :: nx
         !! First dimension.
      integer, intent(in) :: ny
         !! Second dimension.
      real(wp), intent(in) :: f_cor(nx, ny)
         !! Cell-centred Coriolis parameter (1/s).
      real(wp), intent(in) :: cover_frac(nx, ny)
         !! Ice-cover fraction.
      real(wp), intent(in) :: wet_mask(nx, ny)
         !! Static wet/land mask.
      integer :: n_zero
      integer :: i, j
      n_zero = 0
      do j = 1, ny
         do i = 1, nx
            if (cover_frac(i, j) > 0.5_wp .and. wet_mask(i, j) > 0.5_wp) then
               if (f_cor(i, j) == 0.0_wp) n_zero = n_zero + 1
            end if
         end do
      end do
   end function cavity_count_zero_f

   pure function cavity_count_unloaded_p_top(p_top, p_ice_ref, nx, ny) result(n_unloaded)
      !! Count cells whose top-of-column pressure does NOT carry the
      !! isostatic ice load — the melt path's guard that
      !! `configure_ocean_cavity` (the sole `ms%p_top` producer) ran, and
      !! ran before this check.
      !!
      !! The test is `p_top < p_ice_ref`, which is exact rather than a
      !! tolerance: `ms%p_top = p_ice_ref + sf%p_surf` with
      !! `sf%p_surf >= 0` by contract, so a loaded column satisfies it
      !! with no rounding argument at all, and an unloaded one misses it
      !! by the whole 5e6 Pa.  A `pure` predicate, so the refusal can be
      !! tested without provoking the `error stop`.
      integer, intent(in) :: nx
         !! First dimension.
      integer, intent(in) :: ny
         !! Second dimension.
      real(wp), intent(in) :: p_top(nx, ny)
         !! `multilayer_state_t%p_top` (Pa).
      real(wp), intent(in) :: p_ice_ref(nx, ny)
         !! `metrics%p_ice_ref` (Pa) = `rho_ref*GRAVITY*z_draft`.
      integer :: n_unloaded
      integer :: i, j
      n_unloaded = 0
      do j = 1, ny
         do i = 1, nx
            if (p_top(i, j) < p_ice_ref(i, j)) n_unloaded = n_unloaded + 1
         end do
      end do
   end function cavity_count_unloaded_p_top

   subroutine configure_ocean_wetdry(cfg, ocean_state, grid, compute_rank)
      !! Dynamic wetting/drying (docs/ocean_wetdry_plan.md): copy the
      !! `&ocean_wetdry_nml` knobs onto the BT workstate, allocate the
      !! wd_* workspaces (lazy — absent when the knob is off, so the
      !! default path carries no new arrays and stays byte-identical),
      !! and seed the hysteresis wet mask from the seeded bathymetry
      !! (`barotropic%b` — the same array `configure_ocean_bt_split`
      !! later latches into `bt_H_ref`; `bt_eta` is still 0 here).
      !! Must run BEFORE the restart read (so `wd_wet_dyn` is allocated
      !! + registered when the registry walk runs and a warm restart
      !! overwrites the seed with the saved front state) and BEFORE
      !! `ocean_state_enter_data` (host seeding; the enter_data walk
      !! attaches whatever is allocated).
      type(config_t), intent(in) :: cfg
      type(ocean_state_t), intent(inout) :: ocean_state
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: compute_rank
      integer :: nx_w, ny_w, i, j, ng, n_intertidal
      real(wp) :: d0

      ocean_state%dyn%bt_work%wetdry_enable = cfg%ocean%wetdry%enable
      if (.not. cfg%ocean%wetdry%enable) return

      ocean_state%dyn%bt_work%wd_dry_depth = cfg%ocean%wetdry%dry_depth
      ocean_state%dyn%bt_work%wd_rewet_depth = cfg%ocean%wetdry%rewet_depth
      nx_w = grid%nx_total
      ny_w = grid%ny_total
      associate (bt_work => ocean_state%dyn%bt_work)
         allocate (bt_work%wd_wet_dyn(nx_w, ny_w), source=1.0_wp)
         allocate (bt_work%wd_theta(nx_w, ny_w), source=1.0_wp)
         allocate (bt_work%wd_flux_x(nx_w + 1, ny_w), source=0.0_wp)
         allocate (bt_work%wd_flux_y(nx_w, ny_w + 1), source=0.0_wp)
         allocate (bt_work%wd_open_u(nx_w + 1, ny_w), source=1.0_wp)
         allocate (bt_work%wd_open_v(nx_w, ny_w + 1), source=1.0_wp)
         ! Seed the hysteresis mask from the seeded bathymetry (host
         ! loop, pre-enter_data; eta = 0 at configure so D0 = b).  Cells
         ! inside the hysteresis band seed WET — the first substep's own
         ! update settles them.  A restart carrying `wd_wet_dyn`
         ! overwrites this seed (this configure runs BEFORE the registry
         ! read — see the driver ordering comment).
         do j = 1, ny_w
            do i = 1, nx_w
               d0 = ocean_state%barotropic%b(i, j)
               ! `d0 < dry_depth` is the live criterion.  The `wet_mask<=0`
               ! term is belt-and-suspenders: with the knob ON a static-land
               ! column has `b < -land_margin < 0 < dry_depth`, so the first
               ! test already fires — it only guards against a caller seeding
               ! wet_mask with a stricter cutoff than the depth test implies.
               if (d0 < cfg%ocean%wetdry%dry_depth .or. &
                   ocean_state%multilayer%wet_mask(i, j) <= 0.0_wp) then
                  bt_work%wd_wet_dyn(i, j) = 0.0_wp
               end if
            end do
         end do
      end associate
      if (compute_rank == 0) then
         call logger%info("Wet/dry:          ON  dry_depth = "// &
                          to_string(cfg%ocean%wetdry%dry_depth)// &
                          " m, rewet_depth = "// &
                          to_string(cfg%ocean%wetdry%rewet_depth)// &
                          " m (upwind BT faces + positive-definite "// &
                          "outflow limiter + bed-blocking gate)")
         ! v2 semantic trap (plan §10.1): with wet/dry ON the static land
         ! seed cutoff moves to `b < -land_margin`, so any column whose bed
         ! is in the band `[-land_margin, LAND_DEPTH_THRESHOLD)` is NOT
         ! static land any more — it becomes a dynamically-dry INTERTIDAL
         ! flat that floods once eta clears dry_depth.  This affects formula
         ! topos that flag land as b ≈ 0 (island / double_drake), file
         ! bathymetry with a b = 1 m land-flag convention, AND generated
         ! topos like neverworld2 that force continents into the band — so
         ! a name-based check is incomplete.  Scan the seeded interior
         ! bathymetry directly and warn (rank 0) with the exact count when
         ! ANY intertidal-band cell exists.  It's a WARNING, not an abort:
         ! intertidal flats can be intentional; the message just makes the
         ! regime change loud and specific.
         ng = grid%nghost
         n_intertidal = 0
         do j = ng + 1, ny_w - ng
            do i = ng + 1, nx_w - ng
               d0 = ocean_state%barotropic%b(i, j)
               if (d0 >= -cfg%ocean%wetdry%land_margin .and. &
                   d0 < LAND_DEPTH_THRESHOLD) then
                  n_intertidal = n_intertidal + 1
               end if
            end do
         end do
         if (n_intertidal > 0) then
            call logger%warning("Wet/dry: "//to_string(n_intertidal)// &
                                " interior cell(s) have bed elevations in the "// &
                                "intertidal band [-land_margin, "// &
                                "LAND_DEPTH_THRESHOLD) — under the wetdry-aware "// &
                                "static seed (land iff b < -land_margin) these "// &
                                "are LIVE intertidal flats that flood when eta "// &
                                "clears dry_depth, NOT static land.  If you "// &
                                "intended static continents, use bed elevations "// &
                                "b < -land_margin (= "// &
                                to_string(-cfg%ocean%wetdry%land_margin)//" m).")
         end if
      end if
   end subroutine configure_ocean_wetdry

   subroutine configure_ocean_bc(cfg, ocean_state, compute_rank, ierr)
      !! Populate `ocean_state%bc` edge tags + per-edge data values from
      !! `cfg%ocean%bc` (read from `&ocean_bc_nml`).  Run this after
      !! `configure_ocean_bt_split` and before `ocean_state_enter_data`
      !! so the BC state is set before the first dyn step.
      !!
      !! All defaults in `ocean_bc_config_t` resolve to OBC_WALL, so
      !! namelist files that omit `&ocean_bc_nml` are bit-identical to
      !! prior behaviour.
      !!
      !! After populating tags, calls `ocean_bc_validate_periodic` when any
      !! axis is periodic so the ghost-width and sponge incompatibility
      !! checks fire at setup time.
      type(config_t), intent(inout) :: cfg
      type(ocean_state_t), intent(inout) :: ocean_state
      integer, intent(in) :: compute_rank
      integer, intent(out), optional :: ierr
         !! Non-zero on a boundary-condition configuration conflict when
         !! present; absent behaves as today (`error stop`).

      integer :: ntidal, n_tracers
      integer :: i, local_ierr

      associate (bc => ocean_state%bc, bc_cfg => cfg%ocean%bc)

         ! ---- Per-edge tag ----
         bc%west%bc_type = ocean_bc_type_from_string(bc_cfg%west)
         bc%east%bc_type = ocean_bc_type_from_string(bc_cfg%east)
         bc%south%bc_type = ocean_bc_type_from_string(bc_cfg%south)
         bc%north%bc_type = ocean_bc_type_from_string(bc_cfg%north)

         ! ---- Clamped Dirichlet values ----
         bc%west%clamped_eta = bc_cfg%west_clamped_eta
         bc%east%clamped_eta = bc_cfg%east_clamped_eta
         bc%south%clamped_eta = bc_cfg%south_clamped_eta
         bc%north%clamped_eta = bc_cfg%north_clamped_eta
         bc%west%clamped_u = bc_cfg%west_clamped_u
         bc%east%clamped_u = bc_cfg%east_clamped_u
         bc%south%clamped_v = bc_cfg%south_clamped_v
         bc%north%clamped_v = bc_cfg%north_clamped_v

         ! ---- Inflow tracer concentrations ----
         ! Populate clamped_tracer(idx_S) and clamped_tracer(idx_T) for
         ! each edge.  The bc_state_init has already allocated these arrays.
         n_tracers = bc%n_tracers
         if (n_tracers > 0 .and. allocated(bc%west%clamped_tracer)) then
            if (ocean_state%multilayer%idx_salinity > 0 .and. &
                ocean_state%multilayer%idx_salinity <= n_tracers) then
               i = ocean_state%multilayer%idx_salinity
               bc%west%clamped_tracer(i) = bc_cfg%west_inflow_S
               bc%east%clamped_tracer(i) = bc_cfg%east_inflow_S
               bc%south%clamped_tracer(i) = bc_cfg%south_inflow_S
               bc%north%clamped_tracer(i) = bc_cfg%north_inflow_S
            end if
            if (ocean_state%multilayer%idx_temperature > 0 .and. &
                ocean_state%multilayer%idx_temperature <= n_tracers) then
               i = ocean_state%multilayer%idx_temperature
               bc%west%clamped_tracer(i) = bc_cfg%west_inflow_T
               bc%east%clamped_tracer(i) = bc_cfg%east_inflow_T
               bc%south%clamped_tracer(i) = bc_cfg%south_inflow_T
               bc%north%clamped_tracer(i) = bc_cfg%north_inflow_T
            end if
            ! Pseudo-salt OBC mirror: same clamped inflow value as
            ! salinity (pseudo-salt receives exactly salinity's boundary
            ! fluxes, §3.2 of the plan).
            if (ocean_state%multilayer%idx_pseudo_salt > 0 .and. &
                ocean_state%multilayer%idx_pseudo_salt <= n_tracers) then
               i = ocean_state%multilayer%idx_pseudo_salt
               bc%west%clamped_tracer(i) = bc_cfg%west_inflow_S
               bc%east%clamped_tracer(i) = bc_cfg%east_inflow_S
               bc%south%clamped_tracer(i) = bc_cfg%south_inflow_S
               bc%north%clamped_tracer(i) = bc_cfg%north_inflow_S
            end if
         end if

         ! ---- Sponge ----
         bc%west%sponge_width = bc_cfg%sponge_width
         bc%east%sponge_width = bc_cfg%sponge_width
         bc%south%sponge_width = bc_cfg%sponge_width
         bc%north%sponge_width = bc_cfg%sponge_width
         bc%west%sponge_strength = bc_cfg%sponge_strength
         bc%east%sponge_strength = bc_cfg%sponge_strength
         bc%south%sponge_strength = bc_cfg%sponge_strength
         bc%north%sponge_strength = bc_cfg%sponge_strength
         bc%west%sponge_relax_tracers = bc_cfg%sponge_relax_tracers
         bc%east%sponge_relax_tracers = bc_cfg%sponge_relax_tracers
         bc%south%sponge_relax_tracers = bc_cfg%sponge_relax_tracers
         bc%north%sponge_relax_tracers = bc_cfg%sponge_relax_tracers

         ! ---- Tidal constituents ----
         ntidal = min(bc_cfg%west_n_tidal, OBC_MAX_TIDAL_CONSTITUENTS)
         bc%west%n_tidal_constituents = ntidal
         bc%west%tidal_amp(1:ntidal) = bc_cfg%west_tidal_amp(1:ntidal)
         bc%west%tidal_phase(1:ntidal) = bc_cfg%west_tidal_phase(1:ntidal)
         bc%west%tidal_omega(1:ntidal) = bc_cfg%west_tidal_omega(1:ntidal)

         ntidal = min(bc_cfg%east_n_tidal, OBC_MAX_TIDAL_CONSTITUENTS)
         bc%east%n_tidal_constituents = ntidal
         bc%east%tidal_amp(1:ntidal) = bc_cfg%east_tidal_amp(1:ntidal)
         bc%east%tidal_phase(1:ntidal) = bc_cfg%east_tidal_phase(1:ntidal)
         bc%east%tidal_omega(1:ntidal) = bc_cfg%east_tidal_omega(1:ntidal)

         ntidal = min(bc_cfg%south_n_tidal, OBC_MAX_TIDAL_CONSTITUENTS)
         bc%south%n_tidal_constituents = ntidal
         bc%south%tidal_amp(1:ntidal) = bc_cfg%south_tidal_amp(1:ntidal)
         bc%south%tidal_phase(1:ntidal) = bc_cfg%south_tidal_phase(1:ntidal)
         bc%south%tidal_omega(1:ntidal) = bc_cfg%south_tidal_omega(1:ntidal)

         ntidal = min(bc_cfg%north_n_tidal, OBC_MAX_TIDAL_CONSTITUENTS)
         bc%north%n_tidal_constituents = ntidal
         bc%north%tidal_amp(1:ntidal) = bc_cfg%north_tidal_amp(1:ntidal)
         bc%north%tidal_phase(1:ntidal) = bc_cfg%north_tidal_phase(1:ntidal)
         bc%north%tidal_omega(1:ntidal) = bc_cfg%north_tidal_omega(1:ntidal)

         ! ---- OBC tidal nodal/astronomical correction (capability C3) ----
         ! When on, bake the 18.6-yr nodal factor f_c and the equilibrium+nodal
         ! phase (V_c + u_c) into each edge's tidal_fnodal/tidal_arg, using the
         ! SHARED &ocean_tides_nml reference epoch so the boundary tide stays
         ! phase-consistent with the interior body tide (C1).  Host setup —
         ! plain loops, no do concurrent.  Default off ⇒ f=1, arg=0 (the
         ! defaults) ⇒ legacy static-phase OBC sum bit-identical.
         bc%tidal_nodal = bc_cfg%obc_tidal_nodal
         if (bc%tidal_nodal) then
            block
               integer :: yr, mo, dy
               logical :: date_ok
               real(wp) :: dref, dnodal
               real(wp) :: v_all(TIDES_CATALOG_SIZE)
               real(wp) :: f_all(TIDES_CATALOG_SIZE), u_all(TIDES_CATALOG_SIZE)
               character(len=16) :: nodal_str

               ! Equilibrium argument V_c at the model-t=0 reference date.
               call parse_date_string(cfg%ocean%tides%ref_date, yr, mo, dy, date_ok)
               if (.not. date_ok) then
                  call fail("OBC tides: unparseable &ocean_tides_nml ref_date '"// &
                            trim(cfg%ocean%tides%ref_date)//"'", ierr, OCEAN_STATUS_ERR_SETUP)
                  return
               end if
               dref = days_since_1900(yr, mo, dy)

               ! Nodal f/u at the nodal reference date ("" ⇒ ref_date).
               nodal_str = cfg%ocean%tides%nodal_ref_date
               if (len_trim(nodal_str) == 0) nodal_str = cfg%ocean%tides%ref_date
               call parse_date_string(nodal_str, yr, mo, dy, date_ok)
               if (.not. date_ok) then
                  call fail("OBC tides: unparseable &ocean_tides_nml nodal_ref_date '"// &
                            trim(nodal_str)//"'", ierr, OCEAN_STATUS_ERR_SETUP)
                  return
               end if
               dnodal = days_since_1900(yr, mo, dy)

               call equilibrium_arguments(dref, v_all)
               call nodal_fu(dnodal, .true., f_all, u_all)

               if (present(ierr)) then
                  call configure_obc_edge_nodal(bc%west, f_all, u_all, v_all, "west", compute_rank, &
                                                ierr=local_ierr)
                  if (local_ierr /= 0) then
                     ierr = local_ierr
                     return
                  end if
               else
                  call configure_obc_edge_nodal(bc%west, f_all, u_all, v_all, "west", compute_rank)
               end if
               if (present(ierr)) then
                  call configure_obc_edge_nodal(bc%east, f_all, u_all, v_all, "east", compute_rank, &
                                                ierr=local_ierr)
                  if (local_ierr /= 0) then
                     ierr = local_ierr
                     return
                  end if
               else
                  call configure_obc_edge_nodal(bc%east, f_all, u_all, v_all, "east", compute_rank)
               end if
               if (present(ierr)) then
                  call configure_obc_edge_nodal(bc%south, f_all, u_all, v_all, "south", compute_rank, &
                                                ierr=local_ierr)
                  if (local_ierr /= 0) then
                     ierr = local_ierr
                     return
                  end if
               else
                  call configure_obc_edge_nodal(bc%south, f_all, u_all, v_all, "south", compute_rank)
               end if
               if (present(ierr)) then
                  call configure_obc_edge_nodal(bc%north, f_all, u_all, v_all, "north", compute_rank, &
                                                ierr=local_ierr)
                  if (local_ierr /= 0) then
                     ierr = local_ierr
                     return
                  end if
               else
                  call configure_obc_edge_nodal(bc%north, f_all, u_all, v_all, "north", compute_rank)
               end if
            end block
         end if

         ! ---- Open-edge tracer reservoirs (§1, v2) ----
         ! Cache length-scale knobs on bc so kernels can read them without
         ! going back to cfg.  Then allocate reservoir arrays when the feature
         ! is enabled and the edge is open-ish.
         bc%res_lscale_out = bc_cfg%res_lscale_out
         bc%res_lscale_in = bc_cfg%res_lscale_in

         if ((bc%res_lscale_out > 0.0_wp .or. bc%res_lscale_in > 0.0_wp) .and. &
             bc%n_tracers > 0) then
            ! West
            if (bc%west%bc_type /= OBC_WALL .and. bc%west%bc_type /= OBC_PERIODIC) then
               allocate (bc%tres_west(bc%ny_total, bc%nz_ml, bc%n_tracers), source=0.0_wp)
            end if
            ! East
            if (bc%east%bc_type /= OBC_WALL .and. bc%east%bc_type /= OBC_PERIODIC) then
               allocate (bc%tres_east(bc%ny_total, bc%nz_ml, bc%n_tracers), source=0.0_wp)
            end if
            ! South
            if (bc%south%bc_type /= OBC_WALL .and. bc%south%bc_type /= OBC_PERIODIC) then
               allocate (bc%tres_south(bc%nx_total, bc%nz_ml, bc%n_tracers), source=0.0_wp)
            end if
            ! North
            if (bc%north%bc_type /= OBC_WALL .and. bc%north%bc_type /= OBC_PERIODIC) then
               allocate (bc%tres_north(bc%nx_total, bc%nz_ml, bc%n_tracers), source=0.0_wp)
            end if

            ! Seed reservoirs from the interior IC concentration.
            ! We seed from clamped_tracer (which was just populated above from
            ! west/east/south/north_inflow_S/T).  This is the best available
            ! "interior" estimate at setup time; the first-step update will
            ! refine it from the actual h_layer state.  Per-layer init is
            ! uniform in k (the IC is typically also uniform in k at setup time).
            if (allocated(bc%tres_west)) then
               do i = 1, bc%n_tracers
                  bc%tres_west(:, :, i) = bc%west%clamped_tracer(i)
               end do
            end if
            if (allocated(bc%tres_east)) then
               do i = 1, bc%n_tracers
                  bc%tres_east(:, :, i) = bc%east%clamped_tracer(i)
               end do
            end if
            if (allocated(bc%tres_south)) then
               do i = 1, bc%n_tracers
                  bc%tres_south(:, :, i) = bc%south%clamped_tracer(i)
               end do
            end if
            if (allocated(bc%tres_north)) then
               do i = 1, bc%n_tracers
                  bc%tres_north(:, :, i) = bc%north%clamped_tracer(i)
               end do
            end if

            if (compute_rank == 0) then
               call logger%info("OBC reservoirs: res_lscale_out=" &
                                //to_string(bc%res_lscale_out)//"m" &
                                //" res_lscale_in="//to_string(bc%res_lscale_in)//"m")
            end if
         end if

         ! ---- Per-layer Orlanski radiation (§2, v2) ----
         ! Parse the radiation_scheme string to an integer tag (0=anomaly, 1=orlanski).
         ! Allocate rx and u_prev arrays on radiating edges when the scheme is orlanski.
         ! Both arrays are initialised to 0; u_prev will be seeded from the IC velocity
         ! on the first ocean_obc_apply_baroclinic call (cold-start: dhdt = 0 ⇒ rx = 0).
         ! Not restart-registered — restarts cold (known gap).
         bc%orlanski_rx_max = bc_cfg%orlanski_rx_max
         bc%orlanski_gamma = bc_cfg%orlanski_gamma
         bc%nudge_tau_in = bc_cfg%nudge_tau_in
         bc%nudge_tau_out = bc_cfg%nudge_tau_out
         select case (trim(bc_cfg%radiation_scheme))
         case ("orlanski")
            bc%radiation_scheme = 1
         case default   ! "anomaly" or anything else: v1 path
            bc%radiation_scheme = 0
         end select

         if (bc%radiation_scheme == 1) then
            ! West
            if (bc%west%bc_type /= OBC_WALL .and. bc%west%bc_type /= OBC_PERIODIC .and. &
                bc%west%bc_type /= OBC_CLAMPED) then
               allocate (bc%rx_west(bc%ny_total, bc%nz_ml), source=0.0_wp)
               allocate (bc%u_prev_west(bc%ny_total, bc%nz_ml), source=0.0_wp)
            end if
            ! East
            if (bc%east%bc_type /= OBC_WALL .and. bc%east%bc_type /= OBC_PERIODIC .and. &
                bc%east%bc_type /= OBC_CLAMPED) then
               allocate (bc%rx_east(bc%ny_total, bc%nz_ml), source=0.0_wp)
               allocate (bc%u_prev_east(bc%ny_total, bc%nz_ml), source=0.0_wp)
            end if
            ! South
            if (bc%south%bc_type /= OBC_WALL .and. bc%south%bc_type /= OBC_PERIODIC .and. &
                bc%south%bc_type /= OBC_CLAMPED) then
               allocate (bc%rx_south(bc%nx_total, bc%nz_ml), source=0.0_wp)
               allocate (bc%u_prev_south(bc%nx_total, bc%nz_ml), source=0.0_wp)
            end if
            ! North
            if (bc%north%bc_type /= OBC_WALL .and. bc%north%bc_type /= OBC_PERIODIC .and. &
                bc%north%bc_type /= OBC_CLAMPED) then
               allocate (bc%rx_north(bc%nx_total, bc%nz_ml), source=0.0_wp)
               allocate (bc%u_prev_north(bc%nx_total, bc%nz_ml), source=0.0_wp)
            end if

            if (compute_rank == 0) then
               call logger%info("OBC Orlanski: rx_max=" &
                                //to_string(bc%orlanski_rx_max) &
                                //" gamma="//to_string(bc%orlanski_gamma))
            end if
         end if

         ! ---- Full Flather with exterior velocity (§4, v2) ----
         select case (trim(bc_cfg%flather_form))
         case ("full")
            bc%use_full_flather = .true.
         case default   ! "legacy" or anything else
            bc%use_full_flather = .false.
         end select
         bc%ext_u_west = bc_cfg%west_ext_u
         bc%ext_u_east = bc_cfg%east_ext_u
         bc%ext_v_south = bc_cfg%south_ext_v
         bc%ext_v_north = bc_cfg%north_ext_v

         ! ---- Periodic validation (also refreshes periodic_x/y flags) ----
         ! Only need to call if any axis could be periodic.
         if (bc%west%bc_type == OBC_PERIODIC .or. bc%east%bc_type == OBC_PERIODIC .or. &
             bc%south%bc_type == OBC_PERIODIC .or. bc%north%bc_type == OBC_PERIODIC) then
            if (present(ierr)) then
               call ocean_bc_validate_periodic(bc, ierr=local_ierr)
               if (local_ierr /= 0) then
                  ierr = local_ierr
                  return
               end if
            else
               call ocean_bc_validate_periodic(bc)
            end if
         else
            ! Derive periodic flags from final tags even for non-periodic runs.
            bc%periodic_x = .false.
            bc%periodic_y = .false.
         end if

         ! ---- Tripolar north-fold validation (also refreshes north_fold) ----
         ! Always call: it enforces "fold only on north" for every tag set,
         ! and "fold requires periodic w/e + nghost>=3" when the tag is set.
         if (present(ierr)) then
            call ocean_bc_validate_fold(bc, ierr=local_ierr)
            if (local_ierr /= 0) then
               ierr = local_ierr
               return
            end if
         else
            call ocean_bc_validate_fold(bc)
         end if

         if (compute_rank == 0) then
            call logger%info("OBC:  west="//trim(bc_cfg%west)// &
                             " east="//trim(bc_cfg%east)// &
                             " south="//trim(bc_cfg%south)// &
                             " north="//trim(bc_cfg%north))
         end if
      end associate
      if (present(ierr)) ierr = OCEAN_STATUS_OK
   end subroutine configure_ocean_bc

   subroutine configure_ocean_sponge(cfg, ocean_state, grid, compute_rank)
      !! Populate `ocean_state%sponge`'s per-cell `idamp_h`/`idamp_u`/
      !! `idamp_v` maps from `cfg%ocean%sponge` (`&ocean_sponge_nml`) + the
      !! already-configured `ocean_state%bc` edge tags. Run AFTER
      !! `configure_ocean_bc` (reads the edge tags + `has_*` flags) and
      !! BEFORE `ocean_state_enter_data`.
      !!
      !! `damp_source = "band"` (the only value implemented in v1): cosine
      !! ramp from every `OBC_SPONGE`-tagged edge, at the exact cell/u-face/
      !! v-face offsets the legacy `rdb_ocean_sponge::ocean_sponge_apply{,
      !! _tracers}` kernels use (§3.2 of the plan — reproduces today's band
      !! so `enable=.true., damp_source="band"` is *physically* the same
      !! sponge as today). Overlapping edges (corners) SUM their rates — the
      !! exact continuation of the legacy kernel's SEQUENTIAL per-edge decay
      !! (`exp(-a*dt)*exp(-b*dt) = exp(-(a+b)*dt)`).
      !!
      !! Per-edge `west_width`/`west_strength` etc. (sentinel `< 0` ⇒
      !! inherit the single global `&ocean_bc_nml sponge_width`/
      !! `sponge_strength` already fanned out onto `bc%<edge>%sponge_*` by
      !! `configure_ocean_bc`) close the audit gap that the per-edge tag
      !! fields were previously unreachable (no config path could ever make
      !! them differ).
      !!
      !! Does NOT snapshot the reference state — `ocean_sponge_snapshot_reference`
      !! does that, called MUCH earlier from the driver (right after the IC
      !! seed, before any restart read — see that subroutine's docstring for
      !! why the ordering is load-bearing).
      !!
      !! No-op when `.not. ocean_state%sponge%enable` (default) — the maps
      !! stay at their `init`-time zero fill (or unallocated, if `enable`
      !! was false at `init` time too).
      type(config_t), intent(in) :: cfg
      type(ocean_state_t), intent(inout) :: ocean_state
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: compute_rank

      integer :: i0, i1, j0, j1, band, ramp
      real(wp) :: strength
      integer :: nonzero_h, nonzero_u, nonzero_v

      associate (sp => ocean_state%sponge, bc => ocean_state%bc, s_cfg => cfg%ocean%sponge)

         if (.not. sp%enable) return

         sp%relax_uv = s_cfg%relax_uv
         sp%relax_tracers = s_cfg%relax_tracers
         sp%relax_h = s_cfg%relax_h
         sp%damp_source = s_cfg%damp_source
         sp%target_source = s_cfg%target_source
         sp%lin_t_ref = s_cfg%lin_t_ref
         sp%lin_dt_dz = s_cfg%lin_dt_dz
         sp%lin_s_ref = s_cfg%lin_s_ref
         sp%lin_ds_dz = s_cfg%lin_ds_dz
         ramp = merge(SPONGE_RAMP_LINEAR, SPONGE_RAMP_COSINE, &
                      trim(s_cfg%ramp) == "linear")

         ! Latch the two registry indices the analytic `linear_z` refresh
         ! needs, so the per-step kernel never walks the tracer registry
         ! (the array-of-derived-types device-indirection rule).
         sp%idx_t = ocean_state%multilayer%idx_temperature
         sp%idx_s = ocean_state%multilayer%idx_salinity

         ! Latch the geopotential depth of the column TOP.  The draft is
         ! static by design, so this is a configure-time copy and the
         ! sponge slot stays self-contained at run time.  No cavity ⇒
         ! `z_draft` is the (1,1) placeholder and `z_top` keeps its
         ! init-time zero (the open-ocean free-surface datum).
         if (ocean_state%metrics%use_cavity .and. &
             size(ocean_state%metrics%z_draft, 1) == size(sp%z_top, 1) .and. &
             size(ocean_state%metrics%z_draft, 2) == size(sp%z_top, 2)) then
            sp%z_top = ocean_state%metrics%z_draft
         else
            sp%z_top = 0.0_wp
         end if

         i0 = grid%nghost + 1
         i1 = grid%nghost + grid%nx_phys
         j0 = grid%nghost + 1
         j1 = grid%nghost + grid%ny_phys

         if (trim(sp%damp_source) == "band") then
            ! ---- West edge ----
            if (bc%west%bc_type == OBC_SPONGE .and. bc%has_west) then
               band = merge(s_cfg%west_width, bc%west%sponge_width, s_cfg%west_width >= 0)
               strength = merge(s_cfg%west_strength, bc%west%sponge_strength, &
                                s_cfg%west_strength >= 0.0_wp)
               if (band > 0) then
                  call sponge_add_band_x(sp%idamp_h, sp%idamp_u, sp%idamp_v, &
                                         wall_face=grid%nghost + 1, band=band, &
                                         strength=strength, side=+1, j0=j0, j1=j1, &
                                         ramp=ramp)
               end if
            end if
            ! ---- East edge ----
            if (bc%east%bc_type == OBC_SPONGE .and. bc%has_east) then
               band = merge(s_cfg%east_width, bc%east%sponge_width, s_cfg%east_width >= 0)
               strength = merge(s_cfg%east_strength, bc%east%sponge_strength, &
                                s_cfg%east_strength >= 0.0_wp)
               if (band > 0) then
                  call sponge_add_band_x(sp%idamp_h, sp%idamp_u, sp%idamp_v, &
                                         wall_face=grid%nghost + grid%nx_phys + 1, band=band, &
                                         strength=strength, side=-1, j0=j0, j1=j1, &
                                         ramp=ramp)
               end if
            end if
            ! ---- South edge ----
            if (bc%south%bc_type == OBC_SPONGE .and. bc%has_south) then
               band = merge(s_cfg%south_width, bc%south%sponge_width, s_cfg%south_width >= 0)
               strength = merge(s_cfg%south_strength, bc%south%sponge_strength, &
                                s_cfg%south_strength >= 0.0_wp)
               if (band > 0) then
                  call sponge_add_band_y(sp%idamp_h, sp%idamp_u, sp%idamp_v, &
                                         wall_face=grid%nghost + 1, band=band, &
                                         strength=strength, side=+1, i0=i0, i1=i1, &
                                         ramp=ramp)
               end if
            end if
            ! ---- North edge ----
            if (bc%north%bc_type == OBC_SPONGE .and. bc%has_north) then
               band = merge(s_cfg%north_width, bc%north%sponge_width, s_cfg%north_width >= 0)
               strength = merge(s_cfg%north_strength, bc%north%sponge_strength, &
                                s_cfg%north_strength >= 0.0_wp)
               if (band > 0) then
                  call sponge_add_band_y(sp%idamp_h, sp%idamp_u, sp%idamp_v, &
                                         wall_face=grid%nghost + grid%ny_phys + 1, band=band, &
                                         strength=strength, side=-1, i0=i0, i1=i1, &
                                         ramp=ramp)
               end if
            end if
         end if

         if (compute_rank == 0) then
            nonzero_h = count(sp%idamp_h > 0.0_wp)
            nonzero_u = count(sp%idamp_u > 0.0_wp)
            nonzero_v = count(sp%idamp_v > 0.0_wp)
            call logger%info("Sponge (map-driven): damp_source="//trim(sp%damp_source)// &
                             " target_source="//trim(sp%target_source)// &
                             " idamp_h max="//to_string(maxval(sp%idamp_h))//" 1/s"// &
                             " nonzero(h/u/v)="//to_string(nonzero_h)//"/"// &
                             to_string(nonzero_u)//"/"//to_string(nonzero_v))
         end if
      end associate
   end subroutine configure_ocean_sponge

   pure subroutine sponge_add_band_x(idamp_h, idamp_u, idamp_v, wall_face, band, &
                                     strength, side, j0, j1, ramp)
      !! Add a cosine-ramp `Idamp` band for a west/east (x-normal) sponge
      !! edge into the three maps, SUMMING onto whatever is already there
      !! (§3.2 corner composition). Offsets reproduce
      !! `rdb_ocean_sponge::ocean_sponge_apply{,_tracers}` exactly:
      !! `idamp_h`/`idamp_v` share the cell-column offset; `idamp_u` (the
      !! x-normal face) sits one further column in for the west edge (`side
      !! = +1`) and shares the offset for the east edge (`side = -1`) — the
      !! legacy kernel's own asymmetry (see Risk 2 of the plan), reproduced
      !! verbatim so `damp_source="band"` matches today's band exactly.
      real(wp), intent(inout) :: idamp_h(:, :), idamp_u(:, :), idamp_v(:, :)
      integer, intent(in) :: wall_face, band, side, j0, j1
      integer, intent(in) :: ramp
         !! `SPONGE_RAMP_COSINE` (default, the legacy shape, bit-identical)
         !! or `SPONGE_RAMP_LINEAR` (ISOMIP+ Eq. 20 at cell centres).
      real(wp), intent(in) :: strength

      integer :: d, j, i_h, i_u
      real(wp) :: rate, alpha

      do d = 0, band - 1
         ! `sponge_band_alpha` holds both shapes; its cosine branch keeps
         ! the legacy kernel's own `acos(-1.0_wp)` spelling of pi
         ! bit-for-bit (rather than `rdb_constants::PI`, a `4*atan(1)`
         ! formulation) — test 9.4 pins the band to 1e-14.
         alpha = sponge_band_alpha(d, band, ramp)
         rate = strength*alpha
         if (side > 0) then
            i_h = wall_face + d
            i_u = wall_face + d + 1
         else
            i_h = wall_face - d - 1
            i_u = wall_face - d - 1
         end if
         if (i_h >= 1 .and. i_h <= size(idamp_h, 1)) then
            do j = j0, j1
               idamp_h(i_h, j) = idamp_h(i_h, j) + rate
            end do
            do j = j0, j1 + 1
               idamp_v(i_h, j) = idamp_v(i_h, j) + rate
            end do
         end if
         if (i_u >= 1 .and. i_u <= size(idamp_u, 1)) then
            do j = j0, j1
               idamp_u(i_u, j) = idamp_u(i_u, j) + rate
            end do
         end if
      end do
   end subroutine sponge_add_band_x

   pure subroutine sponge_add_band_y(idamp_h, idamp_u, idamp_v, wall_face, band, &
                                     strength, side, i0, i1, ramp)
      !! Mirror of `sponge_add_band_x` for a south/north (y-normal) sponge
      !! edge: `idamp_h`/`idamp_u` share the row offset; `idamp_v` (the
      !! y-normal face) sits one further row in for the south edge (`side =
      !! +1`) and shares the offset for the north edge (`side = -1`).
      real(wp), intent(inout) :: idamp_h(:, :), idamp_u(:, :), idamp_v(:, :)
      integer, intent(in) :: wall_face, band, side, i0, i1
      integer, intent(in) :: ramp
         !! See `sponge_add_band_x`.
      real(wp), intent(in) :: strength

      integer :: d, i, j_h, j_v
      real(wp) :: rate, alpha

      do d = 0, band - 1
         alpha = sponge_band_alpha(d, band, ramp)
         rate = strength*alpha
         if (side > 0) then
            j_h = wall_face + d
            j_v = wall_face + d + 1
         else
            j_h = wall_face - d - 1
            j_v = wall_face - d - 1
         end if
         if (j_h >= 1 .and. j_h <= size(idamp_h, 2)) then
            do i = i0, i1
               idamp_h(i, j_h) = idamp_h(i, j_h) + rate
            end do
            do i = i0, i1 + 1
               idamp_u(i, j_h) = idamp_u(i, j_h) + rate
            end do
         end if
         if (j_v >= 1 .and. j_v <= size(idamp_v, 2)) then
            do i = i0, i1
               idamp_v(i, j_v) = idamp_v(i, j_v) + rate
            end do
         end if
      end do
   end subroutine sponge_add_band_y

   subroutine configure_obc_edge_nodal(face, f_all, u_all, v_all, edge_name, compute_rank, ierr)
      !! Bake the C3 nodal/astronomical correction into one OBC edge, then
      !! fail loud (rank-0 error + `error stop`) if any of the edge's
      !! constituent frequencies matches no tide-catalog entry, else log the
      !! resolved constituent → (name, f_c, V_c+u_c) map on rank 0.  Host-side
      !! setup wrapper around the `pure` `obc_tide_nodal_fill` — the pure
      !! helper signals the failure via `fill_ierr`; the loud policy lives
      !! here.
      type(ocean_bc_face_tag_t), intent(inout) :: face
      real(wp), intent(in) :: f_all(TIDES_CATALOG_SIZE)
      real(wp), intent(in) :: u_all(TIDES_CATALOG_SIZE)
      real(wp), intent(in) :: v_all(TIDES_CATALOG_SIZE)
      character(len=*), intent(in) :: edge_name
      integer, intent(in) :: compute_rank
      integer, intent(out), optional :: ierr
         !! Non-zero when an OBC tidal constituent matches no tide-catalog
         !! entry, when present; absent behaves as today (`error stop`).

      integer :: fill_ierr, nc, ic

      call obc_tide_nodal_fill(face, f_all, u_all, v_all, fill_ierr)
      if (fill_ierr /= 0) then
         call fail("OBC tides: "//trim(edge_name)//" edge constituent " &
                   //to_string(fill_ierr)//" omega="//to_string(face%tidal_omega(fill_ierr)) &
                   //" rad/s matches no tide-catalog entry (rel tol 1e-4)", ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if
      if (compute_rank == 0) then
         do nc = 1, face%n_tidal_constituents
            ic = obc_match_constituent(face%tidal_omega(nc))
            call logger%info("OBC tide nodal ["//trim(edge_name)//"]: "//TIDE_NAME(ic) &
                             //" f="//to_string(face%tidal_fnodal(nc)) &
                             //" V+u="//to_string(face%tidal_arg(nc))//" rad")
         end do
      end if
      if (present(ierr)) ierr = OCEAN_STATUS_OK
   end subroutine configure_obc_edge_nodal

end module rdb_ocean_setup
