!! Bundled test driver — fpm gcc build entry point.
!!
!! Compiles every `tests/test_*.F90` module into one executable that
!! runs every testdrive suite in sequence.
!! CMake builds the same modules as N standalone executables (one per
!! suite) via a `configure_file` template in `tests/CMakeLists.txt`, so
!! per-test isolation under gdb / sanitizers is preserved on that path.
!!
!! Each `use` line renames `collect_<name>_tests` to a unique local
!! identifier (suffix `_collect`) so the bundled program doesn't trip
!! Fortran's name-collision rules when two test modules export the same
!! `collect_*` symbol (e.g. `collect_dam_break_tests` exists in both
!! `test_dam_break` and `test_dam_break_unstr`).
program rdb_tests
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use rdb_comm_env, only: comm_env_finalize
   use test_mem_report, only: mem_report_collect => collect_mem_report_tests
   use test_efp, only: efp_collect => collect_efp_tests
   use test_ocean_console_stats_efp, only: ocean_console_stats_efp_collect => &
                                           collect_ocean_console_stats_efp_tests
   use test_config, only: config_collect => collect_config_tests
   use test_state_roundtrip, only: state_roundtrip_collect => collect_state_roundtrip_tests
   use test_continuity_barotropic, only: continuity_barotropic_collect => collect_continuity_barotropic_tests
   use test_continuity_multilayer, only: continuity_multilayer_collect => collect_continuity_multilayer_tests
   use test_ocean_isopycnal_floor, only: ocean_isopycnal_floor_collect => collect_ocean_isopycnal_floor_tests
   use test_ocean_isopycnal_vanished_vel, only: ocean_isopycnal_vanished_vel_collect => collect_ocean_isopycnal_vanished_vel_tests
   use test_ocean_isopycnal_cfl, only: ocean_isopycnal_cfl_collect => collect_ocean_isopycnal_cfl_tests
   use test_ocean_eos, only: ocean_eos_collect => collect_ocean_eos_tests
   use test_ocean_eos_handle, only: ocean_eos_handle_collect => collect_ocean_eos_handle_tests
   use test_ocean_eos_roquet, only: ocean_eos_roquet_collect => collect_ocean_eos_roquet_tests
   use test_ocean_eos_p_top, only: ocean_eos_p_top_collect => collect_ocean_eos_p_top_tests
   use test_ocean_linear_eos_knobs, only: ocean_linear_eos_knobs_collect => &
                                          collect_ocean_linear_eos_knobs_tests
   use test_ocean_salinity_ic, only: ocean_salinity_ic_collect => &
                                     collect_ocean_salinity_ic_tests
   use test_ocean_freezing_point, only: ocean_freezing_point_collect => collect_ocean_freezing_point_tests
   use test_ocean_cavity_melt, only: ocean_cavity_melt_collect => &
                                     collect_ocean_cavity_melt_tests
   use test_ocean_wright_eos, only: ocean_wright_eos_collect => collect_ocean_wright_eos_tests
   use test_ocean_pgf, only: ocean_pgf_collect => collect_ocean_pgf_tests
   use test_ocean_pgf_fv, only: ocean_pgf_fv_collect => collect_ocean_pgf_fv_tests
   use test_ocean_pgf_grounded, only: ocean_pgf_grounded_collect => collect_ocean_pgf_grounded_tests
   use test_ocean_pgf_fv_mom6, only: ocean_pgf_fv_mom6_collect => collect_ocean_pgf_fv_mom6_tests
   use test_ocean_pgf_reconstruct, only: ocean_pgf_reconstruct_collect => &
                                         collect_ocean_pgf_reconstruct_tests
   use test_ocean_pgf_sigma_rest, only: ocean_pgf_sigma_rest_collect => &
                                        collect_ocean_pgf_sigma_rest_tests
   use test_ocean_pgf_p_top_bc, only: ocean_pgf_p_top_bc_collect => &
                                      collect_ocean_pgf_p_top_bc_tests
   use test_ocean_cavity_draft, only: ocean_cavity_draft_collect => &
                                      collect_ocean_cavity_draft_tests
   use test_ocean_cavity_equivalence, only: ocean_cavity_equivalence_collect => &
                                            collect_ocean_cavity_equivalence_tests
   use test_ocean_cavity_load, only: ocean_cavity_load_collect => &
                                     collect_ocean_cavity_load_tests
   use test_ocean_cavity_flux, only: ocean_cavity_flux_collect => &
                                     collect_ocean_cavity_flux_tests
   use test_ocean_pgf_wright, only: ocean_pgf_wright_collect => collect_ocean_pgf_wright_tests
   use test_ocean_pgf_rho_ref, only: ocean_pgf_rho_ref_collect => collect_ocean_pgf_rho_ref_tests
   use test_ocean_forcing_rho_ref, only: ocean_forcing_rho_ref_collect => &
                                         collect_ocean_forcing_rho_ref_tests
   use test_ocean_vertical_advection, only: ocean_vertical_advection_collect => &
                                            collect_ocean_vertical_advection_tests
   use test_ocean_hdiff_tracer, only: ocean_hdiff_tracer_collect => &
                                      collect_ocean_hdiff_tracer_tests
   use test_ocean_vdiff, only: ocean_vdiff_collect => collect_ocean_vdiff_tests
   use test_ocean_pp81, only: ocean_pp81_collect => collect_ocean_pp81_tests
   use test_ocean_foxkemper, only: ocean_foxkemper_collect => collect_ocean_foxkemper_tests
   use test_ocean_kpp, only: ocean_kpp_collect => collect_ocean_kpp_tests
   use test_ocean_dt_tracer_advect, only: ocean_dt_tracer_advect_collect => &
                                          collect_ocean_dt_tracer_advect_tests
   use test_ocean_tracer_advect_window, only: ocean_tracer_advect_window_collect => &
                                              collect_ocean_tracer_advect_window_tests
   use test_ocean_tracer_weno, only: ocean_tracer_weno_collect => &
                                     collect_ocean_tracer_weno_tests
   use test_ocean_wave_speed, only: ocean_wave_speed_collect => collect_ocean_wave_speed_tests
   use test_ocean_metrics, only: ocean_metrics_collect => collect_ocean_metrics_tests
   use test_ocean_metrics_conservation, only: ocean_metrics_cons_collect => &
                                              collect_ocean_metrics_conservation_tests
   use test_ocean_bipolar, only: ocean_bipolar_collect => collect_ocean_bipolar_tests
   use test_ocean_status_returns, only: ocean_status_returns_collect => &
                                        collect_ocean_status_returns_tests
#ifndef RDB_NO_NETCDF
   use test_ocean_supergrid, only: ocean_supergrid_collect => collect_ocean_supergrid_tests
   use test_ocean_zinit, only: ocean_zinit_collect => collect_ocean_zinit_tests
   use test_ocean_data_input, only: ocean_data_input_collect => collect_ocean_data_input_tests
   use test_ocean_data_forcing, only: ocean_data_forcing_collect => collect_ocean_data_forcing_tests
#endif
   use test_ocean_barotropic_substep, only: ocean_barotropic_substep_collect => collect_ocean_barotropic_substep_tests
   use test_ocean_cor_ref_seiche, only: ocean_cor_ref_seiche_collect => &
                                        collect_ocean_cor_ref_seiche_tests
   use test_ocean_bt_slow_forcing, only: ocean_bt_slow_forcing_collect => &
                                         collect_ocean_bt_slow_forcing_tests
   use test_ocean_tides_astronomy, only: ocean_tides_astronomy_collect => collect_ocean_tides_astronomy_tests
   use test_ocean_tides_disabled_bitident, only: ocean_tides_bitident_collect => &
                                                 collect_ocean_tides_disabled_bitident_tests
   use test_ocean_tidal_forcing, only: ocean_tidal_forcing_collect => collect_ocean_tidal_forcing_tests
   use test_ocean_tides_sal, only: ocean_tides_sal_collect => collect_ocean_tides_sal_tests
   use test_ocean_obc_tide_nodal, only: ocean_obc_tide_nodal_collect => collect_ocean_obc_tide_nodal_tests
   use test_ocean_dyn_split, only: ocean_dyn_split_collect => collect_ocean_dyn_split_tests
   use test_ocean_sponge, only: ocean_sponge_collect => collect_ocean_sponge_tests
   use test_ocean_surface_flux, only: ocean_surface_flux_collect => collect_ocean_surface_flux_tests
   use test_ocean_surface_forcing_type, only: ocean_surface_forcing_type_collect => &
                                              collect_ocean_surface_forcing_type_tests
   use test_ocean_frazil, only: ocean_frazil_collect => collect_ocean_frazil_tests
   use test_ocean_ice_enthalpy, only: ocean_ice_enthalpy_collect => collect_ocean_ice_enthalpy_tests
   use test_ocean_ice_column, only: ocean_ice_column_collect => collect_ocean_ice_column_tests
   use test_ocean_ice_coupling, only: ocean_ice_coupling_collect => collect_ocean_ice_coupling_tests
   use test_ocean_ice_driver_column, only: ocean_ice_driver_column_collect => &
                                           collect_ocean_ice_driver_column_tests
   use test_ocean_ice_snowfall, only: ocean_ice_snowfall_collect => &
                                      collect_ocean_ice_snowfall_tests
   use test_ocean_ice_itd, only: ocean_ice_itd_collect => collect_ocean_ice_itd_tests
   use test_ocean_ice_transport, only: ocean_ice_transport_collect => &
                                       collect_ocean_ice_transport_tests
   use test_ocean_ice_evp, only: ocean_ice_evp_collect => collect_ocean_ice_evp_tests
   use test_ocean_ice_diags, only: ocean_ice_diags_collect => collect_ocean_ice_diags_tests
   use test_ocean_ice_init, only: ocean_ice_init_collect => collect_ocean_ice_init_tests
   use test_ocean_sw_penetration, only: ocean_sw_penetration_collect => collect_ocean_sw_penetration_tests
   use test_ocean_restore, only: ocean_restore_collect => collect_ocean_restore_tests
   use test_ocean_geothermal, only: ocean_geothermal_collect => collect_ocean_geothermal_tests
   use test_ocean_bkgnd_mixing, only: ocean_bkgnd_mixing_collect => collect_ocean_bkgnd_mixing_tests
   use test_ocean_kpp_convective, only: ocean_kpp_convective_collect => collect_ocean_kpp_convective_tests
   use test_ocean_buoyancy_flux, only: ocean_buoyancy_flux_collect => collect_ocean_buoyancy_flux_tests
   use test_ocean_vcoord, only: ocean_vcoord_collect => collect_ocean_vcoord_tests
   use test_ocean_vcoord_rho, only: ocean_vcoord_rho_collect => collect_ocean_vcoord_rho_tests
   use test_ocean_vcoord_hycom, only: ocean_vcoord_hycom_collect => collect_ocean_vcoord_hycom_tests
   use test_ocean_remap, only: ocean_remap_collect => collect_ocean_remap_tests
   use test_ocean_regrid_refine, only: ocean_regrid_refine_collect => collect_ocean_regrid_refine_tests
   use test_ocean_remap_e2e, only: ocean_remap_e2e_collect => collect_ocean_remap_e2e_tests
   use test_ocean_ppm_h4_remap, only: ocean_ppm_h4_remap_collect => collect_ocean_ppm_h4_remap_tests
   use test_safe_math, only: safe_math_collect => collect_safe_math_tests
   use test_ocean_tracer_adv, only: ocean_tracer_adv_collect => collect_ocean_tracer_adv_tests
   use test_ocean_conservation_salt_heat, only: ocean_conservation_salt_heat_collect => &
                                                collect_ocean_conservation_salt_heat_tests
   use test_ocean_dyn_multilayer, only: ocean_dyn_multilayer_collect => collect_ocean_dyn_multilayer_tests
   use test_ocean_hvisc, only: ocean_hvisc_collect => collect_ocean_hvisc_tests
   use test_ocean_hvisc_resoln, only: ocean_hvisc_resoln_collect => collect_ocean_hvisc_resoln_tests
   use test_ocean_leith, only: ocean_leith_collect => collect_ocean_leith_tests
   use test_ocean_isopycnal_slopes, only: ocean_isopycnal_slopes_collect => &
                                          collect_ocean_isopycnal_slopes_tests
   use test_ocean_gm, only: ocean_gm_collect => collect_ocean_gm_tests
   use test_ocean_redi, only: ocean_redi_collect => collect_ocean_redi_tests
   use test_ocean_varmix, only: ocean_varmix_collect => collect_ocean_varmix_tests
   use test_ocean_meke, only: ocean_meke_collect => collect_ocean_meke_tests
   use test_ocean_meke_backscatter, only: ocean_meke_backscatter_collect => collect_ocean_meke_backscatter_tests
   use test_ocean_smag_ah, only: ocean_smag_ah_collect => collect_ocean_smag_ah_tests
   use test_ocean_hvisc_bih_cfl, only: ocean_hvisc_bih_cfl_collect => collect_ocean_hvisc_bih_cfl_tests
   use test_ocean_hvisc_leith_biharm, only: ocean_hvisc_leith_biharm_collect => &
                                            collect_ocean_hvisc_leith_biharm_tests
   use test_ocean_hvisc_aniso, only: ocean_hvisc_aniso_collect => collect_ocean_hvisc_aniso_tests
   use test_ocean_fail_loud_dispatch, only: ocean_fail_loud_dispatch_collect => &
                                            collect_ocean_fail_loud_dispatch_tests
   use test_ocean_ideal_age, only: ocean_ideal_age_collect => collect_ocean_ideal_age_tests
   use test_ocean_diag, only: ocean_diag_collect => collect_ocean_diag_tests
   use test_ocean_diag_reduce, only: ocean_diag_reduce_collect => collect_ocean_diag_reduce_tests
   use test_ocean_porous, only: ocean_porous_collect => collect_ocean_porous_tests
   use test_ocean_diag_remap, only: ocean_diag_remap_collect => collect_ocean_diag_remap_tests
#ifndef RDB_NO_NETCDF
   use test_ocean_diag_netcdf, only: ocean_diag_netcdf_collect => collect_ocean_diag_netcdf_tests
   use test_ocean_restart, only: ocean_restart_collect => collect_ocean_restart_tests
#endif
   use test_ocean_bottom_drag, only: ocean_bottom_drag_collect => collect_ocean_bottom_drag_tests
   use test_ocean_cfl_trunc, only: ocean_cfl_trunc_collect => collect_ocean_cfl_trunc_tests
   use test_ocean_periodic, only: ocean_periodic_collect => collect_ocean_periodic_tests
   use test_ocean_fold, only: ocean_fold_collect => collect_ocean_fold_tests
   use test_ocean_tripolar, only: ocean_tripolar_collect => collect_ocean_tripolar_tests
   use test_ocean_stability_audit, only: ocean_stability_audit_collect => &
                                         collect_ocean_stability_audit_tests
   use test_ocean_obc_baroclinic, only: ocean_obc_baroclinic_collect => collect_ocean_obc_baroclinic_tests
   use test_ocean_obc_eta_ghost, only: ocean_obc_eta_ghost_collect => collect_ocean_obc_eta_ghost_tests
   use test_ocean_surface_stress, only: ocean_surface_stress_collect => collect_ocean_surface_stress_tests
   use test_coriolis_inertial, only: coriolis_inertial_collect => collect_coriolis_inertial_tests
   use test_coriolis_multilayer, only: coriolis_multilayer_collect => collect_coriolis_multilayer_tests
   use test_ocean_dyn_step, only: ocean_dyn_step_collect => collect_ocean_dyn_step_tests
   use test_decomp, only: decomp_collect => collect_decomp_tests
   use test_driver_ocean, only: driver_ocean_collect => collect_driver_ocean_tests
   use test_remap, only: remap_collect => collect_remap_tests
   use test_remap_pqm, only: remap_pqm_collect => collect_remap_pqm_tests
   use test_vcoord, only: vcoord_collect => collect_vcoord_tests
   use test_profiler, only: profiler_collect => collect_profiler_tests
   use test_vcoord_target, only: vcoord_target_collect => collect_vcoord_target_tests
   use test_vcoord_zstar_full, only: vcoord_zstar_full_collect => collect_vcoord_zstar_full_tests
   use test_decomp_extras, only: decomp_extras_collect => collect_decomp_extras_tests
   use test_rdb_ocean_api, only: rdb_ocean_api_collect => collect_rdb_ocean_api_tests
   use test_ocean_api_p2, only: ocean_api_p2_collect => collect_ocean_api_p2_tests
   implicit none

   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mem_report", mem_report_collect), &
                new_testsuite("efp", efp_collect), &
                new_testsuite("ocean_console_stats_efp", ocean_console_stats_efp_collect), &
                new_testsuite("config", config_collect), &
                new_testsuite("state_roundtrip", state_roundtrip_collect), &
                new_testsuite("continuity_barotropic", continuity_barotropic_collect), &
                new_testsuite("continuity_multilayer", continuity_multilayer_collect), &
                new_testsuite("ocean_isopycnal_floor", ocean_isopycnal_floor_collect), &
                new_testsuite("ocean_isopycnal_vanished_vel", ocean_isopycnal_vanished_vel_collect), &
                new_testsuite("ocean_isopycnal_cfl", ocean_isopycnal_cfl_collect), &
                new_testsuite("ocean_eos", ocean_eos_collect), &
                new_testsuite("ocean_eos_handle", ocean_eos_handle_collect), &
                new_testsuite("ocean_eos_roquet", ocean_eos_roquet_collect), &
                new_testsuite("ocean_eos_p_top", ocean_eos_p_top_collect), &
                new_testsuite("ocean_linear_eos_knobs", ocean_linear_eos_knobs_collect), &
                new_testsuite("ocean_salinity_ic", ocean_salinity_ic_collect), &
                new_testsuite("ocean_freezing_point", ocean_freezing_point_collect), &
                new_testsuite("ocean_cavity_melt", ocean_cavity_melt_collect), &
                new_testsuite("ocean_wright_eos", ocean_wright_eos_collect), &
                new_testsuite("ocean_pgf", ocean_pgf_collect), &
                new_testsuite("ocean_dt_tracer_advect", ocean_dt_tracer_advect_collect), &
                new_testsuite("ocean_tracer_advect_window", ocean_tracer_advect_window_collect), &
                new_testsuite("ocean_tracer_weno", ocean_tracer_weno_collect), &
                new_testsuite("ocean_pgf_fv", ocean_pgf_fv_collect), &
                new_testsuite("ocean_pgf_grounded", ocean_pgf_grounded_collect), &
                new_testsuite("ocean_pgf_fv_mom6", ocean_pgf_fv_mom6_collect), &
                new_testsuite("ocean_pgf_reconstruct", ocean_pgf_reconstruct_collect), &
                new_testsuite("ocean_pgf_sigma_rest", ocean_pgf_sigma_rest_collect), &
                new_testsuite("ocean_pgf_p_top_bc", ocean_pgf_p_top_bc_collect), &
                new_testsuite("ocean_cavity_draft", ocean_cavity_draft_collect), &
                new_testsuite("ocean_cavity_equivalence", ocean_cavity_equivalence_collect), &
                new_testsuite("ocean_cavity_load", ocean_cavity_load_collect), &
                new_testsuite("ocean_cavity_flux", ocean_cavity_flux_collect), &
                new_testsuite("ocean_pgf_wright", ocean_pgf_wright_collect), &
                new_testsuite("ocean_pgf_rho_ref", ocean_pgf_rho_ref_collect), &
                new_testsuite("ocean_forcing_rho_ref", ocean_forcing_rho_ref_collect), &
                new_testsuite("ocean_vertical_advection", &
                              ocean_vertical_advection_collect), &
                new_testsuite("ocean_hdiff_tracer", ocean_hdiff_tracer_collect), &
                new_testsuite("ocean_vdiff", ocean_vdiff_collect), &
                new_testsuite("ocean_pp81", ocean_pp81_collect), &
                new_testsuite("ocean_foxkemper", ocean_foxkemper_collect), &
                new_testsuite("ocean_kpp", ocean_kpp_collect), &
                new_testsuite("ocean_wave_speed", ocean_wave_speed_collect), &
                new_testsuite("ocean_metrics", ocean_metrics_collect), &
                new_testsuite("ocean_metrics_conservation", ocean_metrics_cons_collect), &
                new_testsuite("ocean_bipolar", ocean_bipolar_collect), &
                new_testsuite("ocean_status_returns", ocean_status_returns_collect), &
#ifndef RDB_NO_NETCDF
                new_testsuite("ocean_supergrid", ocean_supergrid_collect), &
                new_testsuite("ocean_zinit", ocean_zinit_collect), &
                new_testsuite("ocean_data_input", ocean_data_input_collect), &
                new_testsuite("ocean_data_forcing", ocean_data_forcing_collect), &
#endif
                new_testsuite("ocean_barotropic_substep", ocean_barotropic_substep_collect), &
                new_testsuite("ocean_cor_ref_seiche", ocean_cor_ref_seiche_collect), &
                new_testsuite("ocean_bt_slow_forcing", ocean_bt_slow_forcing_collect), &
                new_testsuite("ocean_tides_astronomy", ocean_tides_astronomy_collect), &
                new_testsuite("ocean_tides_disabled_bitident", ocean_tides_bitident_collect), &
                new_testsuite("ocean_tidal_forcing", ocean_tidal_forcing_collect), &
                new_testsuite("ocean_tides_sal", ocean_tides_sal_collect), &
                new_testsuite("ocean_obc_tide_nodal", ocean_obc_tide_nodal_collect), &
                new_testsuite("ocean_dyn_split", ocean_dyn_split_collect), &
                new_testsuite("ocean_sponge", ocean_sponge_collect), &
                new_testsuite("ocean_surface_flux", ocean_surface_flux_collect), &
                new_testsuite("ocean_surface_forcing_type", ocean_surface_forcing_type_collect), &
                new_testsuite("ocean_frazil", ocean_frazil_collect), &
                new_testsuite("ocean_ice_enthalpy", ocean_ice_enthalpy_collect), &
                new_testsuite("ocean_ice_column", ocean_ice_column_collect), &
                new_testsuite("ocean_ice_coupling", ocean_ice_coupling_collect), &
                new_testsuite("ocean_ice_driver_column", ocean_ice_driver_column_collect), &
                new_testsuite("ocean_ice_snowfall", ocean_ice_snowfall_collect), &
                new_testsuite("ocean_ice_itd", ocean_ice_itd_collect), &
                new_testsuite("ocean_ice_transport", ocean_ice_transport_collect), &
                new_testsuite("ocean_ice_evp", ocean_ice_evp_collect), &
                new_testsuite("ocean_ice_diags", ocean_ice_diags_collect), &
                new_testsuite("ocean_ice_init", ocean_ice_init_collect), &
                new_testsuite("ocean_sw_penetration", ocean_sw_penetration_collect), &
                new_testsuite("ocean_restore", ocean_restore_collect), &
                new_testsuite("ocean_geothermal", ocean_geothermal_collect), &
                new_testsuite("ocean_bkgnd_mixing", ocean_bkgnd_mixing_collect), &
                new_testsuite("ocean_kpp_convective", ocean_kpp_convective_collect), &
                new_testsuite("ocean_buoyancy_flux", ocean_buoyancy_flux_collect), &
                new_testsuite("ocean_vcoord", ocean_vcoord_collect), &
                new_testsuite("ocean_vcoord_rho", ocean_vcoord_rho_collect), &
                new_testsuite("ocean_vcoord_hycom", ocean_vcoord_hycom_collect), &
                new_testsuite("ocean_remap", ocean_remap_collect), &
                new_testsuite("ocean_regrid_refine", ocean_regrid_refine_collect), &
                new_testsuite("ocean_remap_e2e", ocean_remap_e2e_collect), &
                new_testsuite("ocean_ppm_h4_remap", ocean_ppm_h4_remap_collect), &
                new_testsuite("safe_math", safe_math_collect), &
                new_testsuite("ocean_tracer_adv", ocean_tracer_adv_collect), &
                new_testsuite("ocean_conservation_salt_heat", ocean_conservation_salt_heat_collect), &
                new_testsuite("ocean_dyn_multilayer", ocean_dyn_multilayer_collect), &
                new_testsuite("ocean_hvisc", ocean_hvisc_collect), &
                new_testsuite("ocean_hvisc_resoln", ocean_hvisc_resoln_collect), &
                new_testsuite("ocean_leith", ocean_leith_collect), &
                new_testsuite("ocean_isopycnal_slopes", ocean_isopycnal_slopes_collect), &
                new_testsuite("ocean_gm", ocean_gm_collect), &
                new_testsuite("ocean_redi", ocean_redi_collect), &
                new_testsuite("ocean_varmix", ocean_varmix_collect), &
                new_testsuite("ocean_meke", ocean_meke_collect), &
                new_testsuite("ocean_meke_backscatter", ocean_meke_backscatter_collect), &
                new_testsuite("ocean_smag_ah", ocean_smag_ah_collect), &
                new_testsuite("ocean_hvisc_bih_cfl", ocean_hvisc_bih_cfl_collect), &
                new_testsuite("ocean_hvisc_leith_biharm", ocean_hvisc_leith_biharm_collect), &
                new_testsuite("ocean_hvisc_aniso", ocean_hvisc_aniso_collect), &
                new_testsuite("ocean_ideal_age", ocean_ideal_age_collect), &
                new_testsuite("ocean_fail_loud_dispatch", ocean_fail_loud_dispatch_collect), &
                new_testsuite("ocean_diag", ocean_diag_collect), &
                new_testsuite("ocean_diag_reduce", ocean_diag_reduce_collect), &
                new_testsuite("ocean_porous", ocean_porous_collect), &
                new_testsuite("ocean_diag_remap", ocean_diag_remap_collect), &
#ifndef RDB_NO_NETCDF
                new_testsuite("ocean_diag_netcdf", ocean_diag_netcdf_collect), &
                new_testsuite("ocean_restart", ocean_restart_collect), &
#endif
                new_testsuite("ocean_bottom_drag", ocean_bottom_drag_collect), &
                new_testsuite("ocean_cfl_trunc", ocean_cfl_trunc_collect), &
                new_testsuite("ocean_periodic", ocean_periodic_collect), &
                new_testsuite("ocean_fold", ocean_fold_collect), &
                new_testsuite("ocean_tripolar", ocean_tripolar_collect), &
                new_testsuite("ocean_stability_audit", ocean_stability_audit_collect), &
                new_testsuite("ocean_obc_baroclinic", ocean_obc_baroclinic_collect), &
                new_testsuite("ocean_obc_eta_ghost", ocean_obc_eta_ghost_collect), &
                new_testsuite("ocean_surface_stress", ocean_surface_stress_collect), &
                new_testsuite("coriolis_inertial", coriolis_inertial_collect), &
                new_testsuite("coriolis_multilayer", coriolis_multilayer_collect), &
                new_testsuite("ocean_dyn_step", ocean_dyn_step_collect), &
                new_testsuite("decomp", decomp_collect), &
                new_testsuite("driver_ocean", driver_ocean_collect), &
                new_testsuite("remap", remap_collect), &
                new_testsuite("remap_pqm", remap_pqm_collect), &
                new_testsuite("vcoord", vcoord_collect), &
                new_testsuite("profiler", profiler_collect), &
                new_testsuite("vcoord_target", vcoord_target_collect), &
                new_testsuite("vcoord_zstar_full", vcoord_zstar_full_collect), &
                new_testsuite("decomp_extras", decomp_extras_collect), &
                new_testsuite("rdb_ocean_api", rdb_ocean_api_collect), &
                new_testsuite("ocean_api_p2", ocean_api_p2_collect), &
                ]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      ! `parallel=.false.` is REQUIRED, not a preference.  test-drive's
      ! `run_testsuite` defaults to `parallel=.true.` and wraps the per-test-case
      ! loop in `!$omp parallel do`.  gfortran and nvfortran build this tree
      ! without any OpenMP flag, so that directive is inert there — but ifx gets
      ! `-qopenmp` unconditionally (it is how ifx maps `do concurrent` onto
      ! threads), so every test CASE in a suite ran concurrently in one process.
      ! The suites are not written for that: kernels keep module-level persistent
      ! workspaces, the profiler and the error ring are module state, and several
      ! cases share a `save`d config.  Running them serially is the contract the
      ! tests were written against.
      call run_testsuite(testsuites(is)%collect, error_unit, stat, parallel=.false.)
   end do

   ! Finalise the comm-env if a test initialised it (no-op otherwise / on
   ! serial builds) so the bundled runner exits cleanly with or without MPI.
   call comm_env_finalize()

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop 1
   end if

end program rdb_tests
