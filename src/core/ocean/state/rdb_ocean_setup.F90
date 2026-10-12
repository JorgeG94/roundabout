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
   use rdb_constants, only: H_VANISHED, VCOORD_ZSTAR_FULL, VCOORD_ZSTAR
   use rdb_constants, only: FRHAT_ARITHMETIC, FRHAT_HYBRID, parse_frhat_scheme
   use rdb_config, only: config_t, ocean_bt_forcing_visc_rem_on, &
                         ocean_bt_renorm_visc_rem_on, ocean_bt_rem_from_visc_rem_on, &
                         ocean_bt_visc_rem_producer_on
   use rdb_grid, only: hgrid_t
   use rdb_decomp, only: decomp_t
   use rdb_ocean_state, only: ocean_state_t, ocean_state_seed_land_cells, &
                              pgf_nonoverlap_gate_on
   use rdb_ocean_pressure_force, only: OPGF_VARIANT_GPRIME, OPGF_VARIANT_FV_MOM6, &
                                       OPGF_VARIANT_FV_WRIGHT, parse_opgf_variant, &
                                       nonoverlap_vanish_tol_for
   use rdb_eos, only: parse_eos_variant, eos_validate, &
                      EOS_VARIANT_ROQUET_SPV, EOS_VARIANT_ROQUET_RHO, &
                      parse_tfreeze_set, eos_apply_tfreeze_set, &
                      TFREEZE_SET_SEAICE
   use rdb_ocean_vcoord, only: parse_ocean_vcoord_type, VCOORD_EULERIAN_Z, &
                               VCOORD_RHO, VCOORD_HYCOM, VCOORD_LAGRANGIAN, &
                               VCOORD_Z_FIXED, ocean_vcoord_z_fixed_target, &
                               ocean_vcoord_closed_face_masks, &
                               ocean_vcoord_eta0_target, &
                               ocean_vcoord_k_top_from_target, &
                               ocean_vcoord_k_bot_from_target, &
                               ocean_vcoord_set_z_fixed_profile, &
                               ocean_vcoord_count_ledges, &
                               ocean_vcoord_count_bed_steps, ocean_vcoord_t
   use rdb_vcoord, only: parse_remap_method, parse_z_fixed_profile, z_fixed_nominal_dz, &
                         ZFIXED_PROFILE_UNIFORM, ZFIXED_PROFILE_INVALID, ZFIXED_DZ_OK
   use rdb_ocean_bottom_drag, only: parse_bdrag_variant
   use rdb_ocean_vdiff, only: vdiff_bbl_configure, BBL_FORM_LINEAR, BBL_FORM_QUADRATIC
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
                                       OBC_TRIPOLAR_FOLD, ocean_bc_outer_face_tag, &
                                       OBC_MAX_TIDAL_CONSTITUENTS, &
                                       ocean_bc_face_tag_t, &
                                       obc_tide_nodal_fill, obc_match_constituent
   use rdb_ocean_sponge, only: sponge_band_alpha, SPONGE_RAMP_COSINE, SPONGE_RAMP_LINEAR
   use rdb_ocean_halo, only: ocean_halo_centre, ocean_halo_is_decomposed, &
                             ocean_halo_is_decomposed_x, ocean_halo_is_decomposed_y, &
                             ocean_halo_face_x, ocean_halo_face_y, ocean_halo_is_init
   use rdb_ocean_halo_state, only: ocean_seam_refresh_surface_stress
   use rdb_ocean_fold_apply, only: ocean_fold_wrap_eta_2d
   use rdb_ocean_fold_exchange, only: ocean_fold_is_distributed
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
   use rdb_halo, only: halo_allreduce_min, halo_allreduce_sum
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
   public :: configure_ocean_k_bot
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
   ! R2 split (submodules): these stay resident in this parent module (real
   ! bodies, host-associated into the submodules below) but were file-private
   ! before the split.  gfortran -O3 gives a module-contained procedure
   ! INTERNAL linkage whenever it is private and has no caller in the SAME
   ! translation unit -- it has no way to know a submodule calls it via host
   ! association, so the body is silently dropped (link error) or demoted to
   ! a local symbol.  Making them public forces external linkage; no argument
   ! list, intent, or behaviour changes -- a pure widening of the API.
   public :: tags_periodic_x
   public :: tags_north_fold
   public :: metrics_dx_min
   public :: fill_coriolis_corner
   public :: fill_coriolis_centre
   public :: fill_f_corner_seam_ghosts
   public :: hycom_top_dz
   public :: sponge_add_band_x
   public :: sponge_add_band_y

   interface
      module subroutine configure_ocean_metrics(cfg, ocean_state, grid, compute_rank, ierr)
            !! See submodule [[rdb_ocean_setup_grid]].
         implicit none
         type(config_t), intent(in) :: cfg
         type(ocean_state_t), intent(inout) :: ocean_state
         type(hgrid_t), intent(in) :: grid
         integer, intent(in) :: compute_rank
         integer, intent(out), optional :: ierr
      end subroutine configure_ocean_metrics

      module subroutine configure_ocean_land_mask(cfg, ocean_state, grid, compute_rank, warm_restart)
            !! See submodule [[rdb_ocean_setup_grid]].
         implicit none
         type(config_t), intent(in) :: cfg
         type(ocean_state_t), intent(inout) :: ocean_state
         type(hgrid_t), intent(in) :: grid
         integer, intent(in) :: compute_rank
         logical, intent(in), optional :: warm_restart
      end subroutine configure_ocean_land_mask

      module subroutine configure_ocean_forcing(cfg, ocean_state, grid, compute_rank, decomp)
            !! See submodule [[rdb_ocean_setup_grid]].
         implicit none
         type(config_t), intent(in) :: cfg
         type(ocean_state_t), intent(inout) :: ocean_state
         type(hgrid_t), intent(in) :: grid
         integer, intent(in) :: compute_rank
         type(decomp_t), intent(in), optional :: decomp
      end subroutine configure_ocean_forcing

      module subroutine configure_ocean_z_fixed_profile(cfg, ocean_state, compute_rank, log_it)
            !! See submodule [[rdb_ocean_setup_grid]].
         implicit none
         type(config_t), intent(in) :: cfg
         type(ocean_state_t), intent(inout) :: ocean_state
         integer, intent(in) :: compute_rank
         logical, intent(in) :: log_it
      end subroutine configure_ocean_z_fixed_profile

      module subroutine configure_ocean_k_top(cfg, ocean_state, grid, compute_rank)
            !! See submodule [[rdb_ocean_setup_grid]].
         implicit none
         type(config_t), intent(in) :: cfg
         type(ocean_state_t), intent(inout) :: ocean_state
         type(hgrid_t), intent(in) :: grid
         integer, intent(in) :: compute_rank
      end subroutine configure_ocean_k_top

      module subroutine configure_ocean_k_bot(ocean_state, grid, compute_rank)
            !! See submodule [[rdb_ocean_setup_grid]].
         implicit none
         type(ocean_state_t), intent(inout) :: ocean_state
         type(hgrid_t), intent(in) :: grid
         integer, intent(in) :: compute_rank
      end subroutine configure_ocean_k_bot

      module subroutine configure_ocean_closed_faces(cfg, ocean_state, grid, compute_rank, ierr)
            !! See submodule [[rdb_ocean_setup_grid]].
         implicit none
         type(config_t), intent(in) :: cfg
         type(ocean_state_t), intent(inout) :: ocean_state
         type(hgrid_t), intent(in) :: grid
         integer, intent(in) :: compute_rank
         integer, intent(out), optional :: ierr
      end subroutine configure_ocean_closed_faces

      module subroutine refuse_open_zfixed_staircase(ocean_state, grid, ierr)
            !! See submodule [[rdb_ocean_setup_grid]].
         implicit none
         type(ocean_state_t), intent(in) :: ocean_state
         type(hgrid_t), intent(in) :: grid
         integer, intent(out), optional :: ierr
      end subroutine refuse_open_zfixed_staircase

      module subroutine configure_ocean_porous(cfg, ocean_state, grid, compute_rank, ierr)
            !! See submodule [[rdb_ocean_setup_grid]].
         implicit none
         type(config_t), intent(in) :: cfg
         type(ocean_state_t), intent(inout) :: ocean_state
         type(hgrid_t), intent(in) :: grid
         integer, intent(in) :: compute_rank
         integer, intent(out), optional :: ierr
      end subroutine configure_ocean_porous

      module subroutine configure_ocean_wetdry(cfg, ocean_state, grid, compute_rank)
            !! See submodule [[rdb_ocean_setup_grid]].
         implicit none
         type(config_t), intent(in) :: cfg
         type(ocean_state_t), intent(inout) :: ocean_state
         type(hgrid_t), intent(in) :: grid
         integer, intent(in) :: compute_rank
      end subroutine configure_ocean_wetdry

      module subroutine configure_ocean_drag(cfg, ocean_state, compute_rank, ierr)
            !! See submodule [[rdb_ocean_setup_vmix]].
         implicit none
         type(config_t), intent(in) :: cfg
         type(ocean_state_t), intent(inout) :: ocean_state
         integer, intent(in) :: compute_rank
         integer, intent(out), optional :: ierr
      end subroutine configure_ocean_drag

      module subroutine configure_ocean_hdiff(cfg, ocean_state, compute_rank)
            !! See submodule [[rdb_ocean_setup_vmix]].
         implicit none
         type(config_t), intent(in) :: cfg
         type(ocean_state_t), intent(inout) :: ocean_state
         integer, intent(in) :: compute_rank
      end subroutine configure_ocean_hdiff

      module subroutine configure_ocean_vmix(cfg, ocean_state, compute_rank, ierr)
            !! See submodule [[rdb_ocean_setup_vmix]].
         implicit none
         type(config_t), intent(in) :: cfg
         type(ocean_state_t), intent(inout) :: ocean_state
         integer, intent(in) :: compute_rank
         integer, intent(out), optional :: ierr
      end subroutine configure_ocean_vmix

      module subroutine configure_ocean_tracers(cfg, ocean_state, compute_rank)
            !! See submodule [[rdb_ocean_setup_vmix]].
         implicit none
         type(config_t), intent(in) :: cfg
         type(ocean_state_t), intent(inout) :: ocean_state
         integer, intent(in) :: compute_rank
      end subroutine configure_ocean_tracers

      module subroutine configure_ocean_epbl(cfg, ocean_state, grid, compute_rank, ierr)
            !! See submodule [[rdb_ocean_setup_vmix]].
         implicit none
         type(config_t), intent(in) :: cfg
         type(ocean_state_t), intent(inout) :: ocean_state
         type(hgrid_t), intent(in) :: grid
         integer, intent(in) :: compute_rank
         integer, intent(out), optional :: ierr
      end subroutine configure_ocean_epbl

      module subroutine configure_ocean_foxkemper(cfg, ocean_state, compute_rank, ierr)
            !! See submodule [[rdb_ocean_setup_vmix]].
         implicit none
         type(config_t), intent(in) :: cfg
         type(ocean_state_t), intent(inout) :: ocean_state
         integer, intent(in) :: compute_rank
         integer, intent(out), optional :: ierr
      end subroutine configure_ocean_foxkemper

      module subroutine configure_ocean_kappa_shear(cfg, ocean_state, grid, compute_rank, ierr)
            !! See submodule [[rdb_ocean_setup_vmix]].
         implicit none
         type(config_t), intent(in) :: cfg
         type(ocean_state_t), intent(inout) :: ocean_state
         type(hgrid_t), intent(in) :: grid
         integer, intent(in) :: compute_rank
         integer, intent(out), optional :: ierr
      end subroutine configure_ocean_kappa_shear

      module subroutine configure_ocean_tidal_mixing(cfg, ocean_state, compute_rank, ierr)
            !! See submodule [[rdb_ocean_setup_vmix]].
         implicit none
         type(config_t), intent(in) :: cfg
         type(ocean_state_t), intent(inout) :: ocean_state
         integer, intent(in) :: compute_rank
         integer, intent(out), optional :: ierr
      end subroutine configure_ocean_tidal_mixing

      module subroutine configure_ocean_conv(cfg, ocean_state, compute_rank, ierr)
            !! See submodule [[rdb_ocean_setup_vmix]].
         implicit none
         type(config_t), intent(in) :: cfg
         type(ocean_state_t), intent(inout) :: ocean_state
         integer, intent(in) :: compute_rank
         integer, intent(out), optional :: ierr
      end subroutine configure_ocean_conv

      module subroutine configure_ocean_ddiff(cfg, ocean_state, compute_rank, ierr)
            !! See submodule [[rdb_ocean_setup_vmix]].
         implicit none
         type(config_t), intent(in) :: cfg
         type(ocean_state_t), intent(inout) :: ocean_state
         integer, intent(in) :: compute_rank
         integer, intent(out), optional :: ierr
      end subroutine configure_ocean_ddiff

      module subroutine configure_ocean_lateral(cfg, ocean_state, grid, compute_rank, ierr)
            !! See submodule [[rdb_ocean_setup_lateral]].
         implicit none
         type(config_t), intent(in) :: cfg
         type(ocean_state_t), intent(inout) :: ocean_state
         type(hgrid_t), intent(in) :: grid
         integer, intent(in) :: compute_rank
         integer, intent(out), optional :: ierr
      end subroutine configure_ocean_lateral

      module subroutine configure_rho_target(cfg, rho_target, nz_ml)
            !! See submodule [[rdb_ocean_setup_lateral]].
         implicit none
         type(config_t), intent(in) :: cfg
         real(wp), intent(inout) :: rho_target(0:)
         integer, intent(in) :: nz_ml
      end subroutine configure_rho_target

      module subroutine configure_ocean_wavespeed(cfg, ocean_state, grid, ierr)
            !! See submodule [[rdb_ocean_setup_lateral]].
         implicit none
         type(config_t), intent(in) :: cfg
         type(ocean_state_t), intent(inout) :: ocean_state
         type(hgrid_t), intent(in) :: grid
         integer, intent(out), optional :: ierr
      end subroutine configure_ocean_wavespeed

      module subroutine configure_ocean_varmix(cfg, ocean_state, grid)
            !! See submodule [[rdb_ocean_setup_lateral]].
         implicit none
         type(config_t), intent(in) :: cfg
         type(ocean_state_t), intent(inout) :: ocean_state
         type(hgrid_t), intent(in) :: grid
      end subroutine configure_ocean_varmix

      module subroutine configure_ocean_meke(cfg, ocean_state, grid)
            !! See submodule [[rdb_ocean_setup_lateral]].
         implicit none
         type(config_t), intent(in) :: cfg
         type(ocean_state_t), intent(inout) :: ocean_state
         type(hgrid_t), intent(in) :: grid
      end subroutine configure_ocean_meke

      module subroutine configure_ocean_reference_density(ocean_state, geo)
            !! See submodule [[rdb_ocean_setup_lateral]].
         implicit none
         type(ocean_state_t), intent(inout) :: ocean_state
         type(ocean_geothermal_t), intent(inout), optional :: geo
      end subroutine configure_ocean_reference_density

      module subroutine configure_ocean_pgf(cfg, ocean_state, compute_rank, ierr)
            !! See submodule [[rdb_ocean_setup_lateral]].
         implicit none
         type(config_t), intent(in) :: cfg
         type(ocean_state_t), intent(inout) :: ocean_state
         integer, intent(in) :: compute_rank
         integer, intent(out), optional :: ierr
      end subroutine configure_ocean_pgf

      module subroutine configure_ocean_tides(cfg, ocean_state, grid, compute_rank, ierr)
            !! See submodule [[rdb_ocean_setup_bt]].
         implicit none
         type(config_t), intent(in) :: cfg
         type(ocean_state_t), intent(inout) :: ocean_state
         type(hgrid_t), intent(in) :: grid
         integer, intent(in) :: compute_rank
         integer, intent(out), optional :: ierr
      end subroutine configure_ocean_tides

      module subroutine configure_ocean_p_surf(cfg, ocean_state, grid, compute_rank)
            !! See submodule [[rdb_ocean_setup_bt]].
         implicit none
         type(config_t), intent(in) :: cfg
         type(ocean_state_t), intent(inout) :: ocean_state
         type(hgrid_t), intent(in) :: grid
         integer, intent(in) :: compute_rank
      end subroutine configure_ocean_p_surf

      module subroutine configure_ocean_wave_drag(cfg, ocean_state, grid, compute_rank, ierr)
            !! See submodule [[rdb_ocean_setup_bt]].
         implicit none
         type(config_t), intent(in) :: cfg
         type(ocean_state_t), intent(inout) :: ocean_state
         type(hgrid_t), intent(in) :: grid
         integer, intent(in) :: compute_rank
         integer, intent(out), optional :: ierr
      end subroutine configure_ocean_wave_drag

      module subroutine wave_drag_roughness_proxy(b, wet_T, nx, ny, nghost, kappa, n_bot, h2_max, r_h)
            !! See submodule [[rdb_ocean_setup_bt]].
         implicit none
         integer, intent(in) :: nx, ny, nghost
         real(wp), intent(in) :: b(nx, ny), wet_T(nx, ny)
         real(wp), intent(in) :: kappa, n_bot, h2_max
         real(wp), intent(inout) :: r_h(nx, ny)
      end subroutine wave_drag_roughness_proxy

      module subroutine parse_constituent_list(list, cat_idx, nconst, ierr)
            !! See submodule [[rdb_ocean_setup_bt]].
         implicit none
         character(len=*), intent(in) :: list
         integer, intent(out) :: cat_idx(:)
         integer, intent(out) :: nconst
         integer, intent(out), optional :: ierr
      end subroutine parse_constituent_list

      module subroutine configure_ocean_bt(cfg, ocean_state, grid, compute_rank)
            !! See submodule [[rdb_ocean_setup_bt]].
         implicit none
         type(config_t), intent(in) :: cfg
         type(ocean_state_t), intent(inout) :: ocean_state
         type(hgrid_t), intent(in) :: grid
         integer, intent(in) :: compute_rank
      end subroutine configure_ocean_bt

      module subroutine configure_ocean_bt_split(cfg, ocean_state, grid, compute_rank)
            !! See submodule [[rdb_ocean_setup_bt]].
         implicit none
         type(config_t), intent(inout) :: cfg
         type(ocean_state_t), intent(inout) :: ocean_state
         type(hgrid_t), intent(in) :: grid
         integer, intent(in) :: compute_rank
      end subroutine configure_ocean_bt_split

      module subroutine configure_ocean_bc(cfg, ocean_state, compute_rank, ierr)
            !! See submodule [[rdb_ocean_setup_bc]].
         implicit none
         type(config_t), intent(inout) :: cfg
         type(ocean_state_t), intent(inout) :: ocean_state
         integer, intent(in) :: compute_rank
         integer, intent(out), optional :: ierr
      end subroutine configure_ocean_bc

      module subroutine configure_ocean_sponge(cfg, ocean_state, grid, compute_rank)
            !! See submodule [[rdb_ocean_setup_bc]].
         implicit none
         type(config_t), intent(in) :: cfg
         type(ocean_state_t), intent(inout) :: ocean_state
         type(hgrid_t), intent(in) :: grid
         integer, intent(in) :: compute_rank
      end subroutine configure_ocean_sponge

      module subroutine configure_obc_edge_nodal(face, f_all, u_all, v_all, edge_name, compute_rank, ierr)
            !! See submodule [[rdb_ocean_setup_bc]].
         implicit none
         type(ocean_bc_face_tag_t), intent(inout) :: face
         real(wp), intent(in) :: f_all(TIDES_CATALOG_SIZE)
         real(wp), intent(in) :: u_all(TIDES_CATALOG_SIZE)
         real(wp), intent(in) :: v_all(TIDES_CATALOG_SIZE)
         character(len=*), intent(in) :: edge_name
         integer, intent(in) :: compute_rank
         integer, intent(out), optional :: ierr
      end subroutine configure_obc_edge_nodal

      module subroutine configure_ocean_cavity(cfg, ocean_state, grid, compute_rank, ierr)
            !! See submodule [[rdb_ocean_setup_cavity]].
         implicit none
         type(config_t), intent(in) :: cfg
         type(ocean_state_t), intent(inout) :: ocean_state
         type(hgrid_t), intent(in) :: grid
         integer, intent(in) :: compute_rank
         integer, intent(out), optional :: ierr
      end subroutine configure_ocean_cavity

      module subroutine configure_ocean_cavity_melt(cfg, ocean_state, grid, compute_rank, ierr)
            !! See submodule [[rdb_ocean_setup_cavity]].
         implicit none
         type(config_t), intent(in) :: cfg
         type(ocean_state_t), intent(inout) :: ocean_state
         type(hgrid_t), intent(in) :: grid
         integer, intent(in) :: compute_rank
         integer, intent(out), optional :: ierr
      end subroutine configure_ocean_cavity_melt

      module subroutine configure_ocean_top_drag(cfg, ocean_state, grid, compute_rank)
            !! See submodule [[rdb_ocean_setup_cavity]].
         implicit none
         type(config_t), intent(in) :: cfg
         type(ocean_state_t), intent(inout) :: ocean_state
         type(hgrid_t), intent(in) :: grid
         integer, intent(in) :: compute_rank
      end subroutine configure_ocean_top_drag
   end interface

contains

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

   subroutine fill_f_corner_seam_ghosts(grid, f_corner)
      !! Periodic-x wrap + north-fold (scalar copy) of the static corner
      !! Coriolis array for a tripolar grid.  f is reflection-invariant
      !! (same latitude at the conjugate corner), so negate=.false.
      !!
      !! Decomposed: only an axis this tile holds WHOLE is wrapped, and only
      !! the tile that owns the north edge folds (the `seed_wrap_static_2d`
      !! rule).  Every other seam ghost was already filled by the generator
      !! from the global-slice geography (`metrics_fill_tripolar`), which is
      !! exactly the neighbour's corner value.
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(inout) :: f_corner(:, :)
      integer :: ng, ni, i, j
      ng = grid%nghost
      ni = grid%nx_phys
      ! px > 1: x ghosts come from the generator; the fold runs later, in
      ! engine_setup, through the distributed exchange (needs its plan).
      if (grid%nx_phys /= grid%nx_global) return
      ! Periodic-x ghost columns (corner/Bu face-type: physical i=ng+1..ng+ni+1).
      do j = 1, size(f_corner, 2)
         do i = 1, ng
            f_corner(i, j) = f_corner(i + ni, j)
            f_corner(ng + ni + 1 + i, j) = f_corner(ng + 1 + i, j)
         end do
      end do
      ! North fold (copy) — the north-edge tile only.
      if (grid%j_offset_global + grid%ny_phys /= grid%ny_global) return
      call fold_north_corner(f_corner, grid%nx_total + 1, grid%ny_total + 1, &
                             grid%nx_phys, grid%ny_phys, grid%nghost, negate=.false.)
   end subroutine fill_f_corner_seam_ghosts

   pure function hycom_top_dz(vcoord) result(dz)
      !! Surface-layer thickness of the `hycom` z* nominal floor (m), for
      !! the configure banner: the stretched profile's top entry, else the
      !! uniform `z_fixed_h_ref/nz`, else 0 (the unconfigured sigma
      !! fallback of `ocean_vcoord_rho_target_column`).
      type(ocean_vcoord_t), intent(in) :: vcoord
      real(wp) :: dz
      dz = 0.0_wp
      if (vcoord%z_fixed_use_profile) then
         dz = vcoord%z_fixed_dz(vcoord%nz_ml)
      else if (vcoord%z_fixed_h_ref > 0.0_wp .and. vcoord%nz_ml > 0) then
         dz = vcoord%z_fixed_h_ref/real(vcoord%nz_ml, wp)
      end if
   end function hycom_top_dz

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

end module rdb_ocean_setup
