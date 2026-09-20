# --------------------------------------------------------------------------
# Test build
#
# Each test source file (`tests/test_*.F90`) is *module-only* — it exports a
# `collect_*_tests` subroutine for the testdrive framework. The matching program
# driver is auto-generated from `tester_template.F90.in` so we preserve the
# per-test executable model (one binary per suite, easy to run under gdb /
# sanitizers / `ctest -R`) without checking in 48 nearly identical driver files.
#
# The same modules are bundled into a single `fpm test` executable via
# `tests/tester.F90`. The flat `tests/` layout is required because fpm's
# `[[test]]` source-dir doesn't recurse subdirs (verified against fpm 0.13).
# Keep the two paths in sync when adding a new test: - drop the module file in
# `tests/` - append it to `RDB_TESTS` below (CMake) - append it to
# `tests/tester.F90` (fpm)
# --------------------------------------------------------------------------

# Helper modules used by some unit tests (Ritter dam-break + Thacker bowl
# analytical solutions). Built as a static lib and linked into every test.
add_library(rdb_analytical STATIC ritter_dambreak.F90 thacker_bowl.F90
                                  ocean_test_metrics.F90)
set_target_properties(rdb_analytical PROPERTIES Fortran_MODULE_DIRECTORY
                                                ${CMAKE_BINARY_DIR}/modules)
target_include_directories(rdb_analytical PRIVATE "${CMAKE_BINARY_DIR}/modules")
# `ocean_test_metrics` uses `rdb_ocean_metrics`, whose own `use` chain reaches
# `rdb_io_netcdf` — so the compiler must open `netcdf.mod` while walking the
# .mod tree, even though nothing here calls NetCDF. Every test target picks the
# dir up transitively from ${main_lib}; this target links no core (see below),
# so take the include dirs alone and keep the no-link property intact.
if(RDB_ENABLE_NETCDF)
  target_include_directories(
    rdb_analytical
    PRIVATE $<TARGET_PROPERTY:netCDF::netcdff,INTERFACE_INCLUDE_DIRECTORIES>)
endif()
# Same reasoning for pic. `rdb_constants` re-exports `pic_types`, so a compiler
# that opens the TRANSITIVE .mod files while walking the chain needs pic's
# module dir on the include path -- nvfortran does, and fails with "Unable to
# open MODULE file pic_types.mod"; gfortran does not, because its .mod files are
# self-contained, which is why this was invisible on gfortran.
#
# It was also invisible before pic was pinned: pic <= v0.8.1 emitted its .mod
# files into ${CMAKE_BINARY_DIR}/modules, the same dir this target already adds,
# so the include path was satisfied by coincidence. v0.9.0 emits them into its
# own _deps/pic-build/modules instead. Take the include dirs alone, keeping the
# no-link property the block above depends on.
target_include_directories(
  rdb_analytical
  PRIVATE $<TARGET_PROPERTY:pic::pic,INTERFACE_INCLUDE_DIRECTORIES>)
# Depend on the core for build ordering + the shared .mod dir, but do NOT LINK
# it: a static lib propagates its PRIVATE link deps to consumers, so linking the
# static core here would drag it into EVERY test and defeat the shared-core
# split below. Each test links its own core; analytical's core symbols resolve
# against whichever one the test picks. (Archiving a static lib does no symbol
# resolution, so the missing link is harmless at this stage.)
#
# THIS TARGET IS WHY THE BUILD REQUIRES THE MAKEFILES GENERATOR.  These sources
# `use rdb_constants`, so they must not compile until the core's .mod files
# exist -- but `add_dependencies` does NOT express that.  It only chains the
# object-ORDER phonies (cmake_object_order_depends_target_*), which never depend
# on the core's actual outputs.  Makefiles happens to serialise it anyway
# because make recurses directories in order; Ninja's single global graph does
# not, so the analytical objects race ahead of the core and the build dies on
# "Unable to open MODULE file rdb_constants.mod".
#
# It is NOT fixable by tightening the dependency: CMake feeds the Ninja Fortran
# dyndep scanner only the module dirs of targets this one LINKS, so
# `add_dependencies` and `$<COMPILE_ONLY:${core_objs}>` are both invisible to it
# (verified: neither produces a .tdi entry, both still race).  Making it
# Ninja-clean means giving this target a real link edge to a core, which is
# exactly what the paragraph above forbids -- a test that links the shared core
# would then drag the static one in too.  Switching generators therefore needs
# that no-link constraint revisited first; don't just flip CMAKE_GENERATOR.
add_dependencies(rdb_analytical ${main_lib})

# Per-test metadata.
#
# Each row encodes one test suite as a single string, fields separated by `|`:
# <test-name>|<module-name>|<collect-fn>|<suite-display>|<regime>
#
# Where: <test-name>     — file stem (also the ctest name and executable name)
# <module-name>   — Fortran module name to `use` (matches <test-name> in every
# current row; column kept so a future test with a different module name can
# plug in trivially) <collect-fn>    — testdrive collector subroutine in that
# module <suite-display> — string passed to testdrive's `new_testsuite()`
# <regime>        — ctest LABEL for the regime the test belongs to: `ocean`
# (C-grid dyn-core + sea-ice), `core` (shared infrastructure), or `ocean+core`
# (both labels).  Mandatory — configure fails loud on a missing/invalid regime,
# and the `test-regime-labels` pre-commit hook enforces it too.  Run a single
# regime suite with `ctest -L ocean` / `-L core`.
set(RDB_TESTS
    "test_mem_report|test_mem_report|collect_mem_report_tests|mem_report|core"
    "test_efp|test_efp|collect_efp_tests|efp|core"
    "test_ocean_console_stats_efp|test_ocean_console_stats_efp|collect_ocean_console_stats_efp_tests|ocean_console_stats_efp|ocean"
    "test_config|test_config|collect_config_tests|config|ocean+core"
    "test_nml_schema|test_nml_schema|collect_nml_schema_tests|nml_schema|core"
    "test_config_schema|test_config_schema|collect_config_schema_tests|config_schema|core"
    "test_ocean_knob_retirement|test_ocean_knob_retirement|collect_ocean_knob_retirement_tests|ocean_knob_retirement|ocean"
    "test_state_roundtrip|test_state_roundtrip|collect_state_roundtrip_tests|state_roundtrip|ocean"
    "test_continuity_barotropic|test_continuity_barotropic|collect_continuity_barotropic_tests|continuity_barotropic|ocean"
    "test_continuity_multilayer|test_continuity_multilayer|collect_continuity_multilayer_tests|continuity_multilayer|ocean"
    "test_ocean_isopycnal_floor|test_ocean_isopycnal_floor|collect_ocean_isopycnal_floor_tests|ocean_isopycnal_floor|ocean"
    "test_ocean_isopycnal_vanished_vel|test_ocean_isopycnal_vanished_vel|collect_ocean_isopycnal_vanished_vel_tests|ocean_isopycnal_vanished_vel|ocean"
    "test_ocean_isopycnal_cfl|test_ocean_isopycnal_cfl|collect_ocean_isopycnal_cfl_tests|ocean_isopycnal_cfl|ocean"
    "test_ocean_eos|test_ocean_eos|collect_ocean_eos_tests|ocean_eos|ocean"
    "test_ocean_eos_handle|test_ocean_eos_handle|collect_ocean_eos_handle_tests|ocean_eos_handle|ocean"
    "test_ocean_eos_roquet|test_ocean_eos_roquet|collect_ocean_eos_roquet_tests|ocean_eos_roquet|ocean"
    "test_ocean_eos_p_top|test_ocean_eos_p_top|collect_ocean_eos_p_top_tests|ocean_eos_p_top|ocean"
    "test_ocean_linear_eos_knobs|test_ocean_linear_eos_knobs|collect_ocean_linear_eos_knobs_tests|ocean_linear_eos_knobs|ocean"
    "test_ocean_salinity_ic|test_ocean_salinity_ic|collect_ocean_salinity_ic_tests|ocean_salinity_ic|ocean"
    "test_ocean_freezing_point|test_ocean_freezing_point|collect_ocean_freezing_point_tests|ocean_freezing_point|ocean"
    "test_ocean_cavity_melt|test_ocean_cavity_melt|collect_ocean_cavity_melt_tests|ocean_cavity_melt|ocean"
    "test_ocean_wright_eos|test_ocean_wright_eos|collect_ocean_wright_eos_tests|ocean_wright_eos|ocean"
    "test_ocean_pgf|test_ocean_pgf|collect_ocean_pgf_tests|ocean_pgf|ocean"
    "test_ocean_pgf_fv|test_ocean_pgf_fv|collect_ocean_pgf_fv_tests|ocean_pgf_fv|ocean"
    "test_ocean_pgf_grounded|test_ocean_pgf_grounded|collect_ocean_pgf_grounded_tests|ocean_pgf_grounded|ocean"
    "test_ocean_pgf_mont|test_ocean_pgf_mont|collect_ocean_pgf_mont_tests|ocean_pgf_mont|ocean"
    "test_ocean_pgf_fv_mom6|test_ocean_pgf_fv_mom6|collect_ocean_pgf_fv_mom6_tests|ocean_pgf_fv_mom6|ocean"
    "test_ocean_pgf_reconstruct|test_ocean_pgf_reconstruct|collect_ocean_pgf_reconstruct_tests|ocean_pgf_reconstruct|ocean"
    "test_ocean_pgf_p_top_bc|test_ocean_pgf_p_top_bc|collect_ocean_pgf_p_top_bc_tests|ocean_pgf_p_top_bc|ocean"
    "test_ocean_pgf_wright|test_ocean_pgf_wright|collect_ocean_pgf_wright_tests|ocean_pgf_wright|ocean"
    "test_ocean_pgf_gprime|test_ocean_pgf_gprime|collect_ocean_pgf_gprime_tests|ocean_pgf_gprime|ocean"
    "test_ocean_pgf_gfs_scale|test_ocean_pgf_gfs_scale|collect_ocean_pgf_gfs_scale_tests|ocean_pgf_gfs_scale|ocean"
    "test_ocean_pgf_rho_ref|test_ocean_pgf_rho_ref|collect_ocean_pgf_rho_ref_tests|ocean_pgf_rho_ref|ocean"
    "test_ocean_vertical_advection|test_ocean_vertical_advection|collect_ocean_vertical_advection_tests|ocean_vertical_advection|ocean"
    "test_ocean_hdiff_tracer|test_ocean_hdiff_tracer|collect_ocean_hdiff_tracer_tests|ocean_hdiff_tracer|ocean"
    "test_ocean_vdiff|test_ocean_vdiff|collect_ocean_vdiff_tests|ocean_vdiff|ocean"
    "test_ocean_vdiff_implicit_stress_drag|test_ocean_vdiff_implicit_stress_drag|collect_ocean_vdiff_implicit_stress_drag_tests|ocean_vdiff_implicit_stress_drag|ocean"
    "test_ocean_visc_rem|test_ocean_visc_rem|collect_ocean_visc_rem_tests|ocean_visc_rem|ocean"
    "test_ocean_vdiff_bbl|test_ocean_vdiff_bbl|collect_ocean_vdiff_bbl_tests|ocean_vdiff_bbl|ocean"
    "test_ocean_pp81|test_ocean_pp81|collect_ocean_pp81_tests|ocean_pp81|ocean"
    "test_ocean_kv_ml_invz2|test_ocean_kv_ml_invz2|collect_ocean_kv_ml_invz2_tests|ocean_kv_ml_invz2|ocean"
    "test_ocean_dt_therm|test_ocean_dt_therm|collect_ocean_dt_therm_tests|ocean_dt_therm|ocean"
    "test_ocean_dt_tracer_advect|test_ocean_dt_tracer_advect|collect_ocean_dt_tracer_advect_tests|ocean_dt_tracer_advect|ocean"
    "test_ocean_tracer_advect_window|test_ocean_tracer_advect_window|collect_ocean_tracer_advect_window_tests|ocean_tracer_advect_window|ocean"
    "test_ocean_tracer_weno|test_ocean_tracer_weno|collect_ocean_tracer_weno_tests|ocean_tracer_weno|ocean+core"
    "test_ocean_bathymetry_fill|test_ocean_bathymetry_fill|collect_ocean_bathymetry_fill_tests|ocean_bathymetry_fill|ocean"
    "test_ocean_neverworld2|test_ocean_neverworld2|collect_ocean_neverworld2_tests|ocean_neverworld2|ocean"
    "test_ocean_layer_rho_init|test_ocean_layer_rho_init|collect_ocean_layer_rho_init_tests|ocean_layer_rho_init|ocean"
    "test_ocean_thickness_config|test_ocean_thickness_config|collect_ocean_thickness_config_tests|ocean_thickness_config|ocean"
    "test_ocean_kpp|test_ocean_kpp|collect_ocean_kpp_tests|ocean_kpp|ocean"
    "test_ocean_metrics|test_ocean_metrics|collect_ocean_metrics_tests|ocean_metrics|ocean"
    "test_ocean_land_mask|test_ocean_land_mask|collect_ocean_land_mask_tests|ocean_land_mask|ocean"
    "test_ocean_metrics_conservation|test_ocean_metrics_conservation|collect_ocean_metrics_conservation_tests|ocean_metrics_conservation|ocean"
    "test_ocean_supergrid|test_ocean_supergrid|collect_ocean_supergrid_tests|ocean_supergrid|ocean"
    # P0.1 review F1/F3: the ierr-present branches of the solver-creation status
    # contract, including the F1 regression gate (missing supergrid file must
    # return ERR_IO, not error stop). NetCDF I/O.
    "test_ocean_status_returns|test_ocean_status_returns|collect_ocean_status_returns_tests|ocean_status_returns|ocean"
    "test_ocean_bipolar|test_ocean_bipolar|collect_ocean_bipolar_tests|ocean_bipolar|ocean"
    "test_ocean_zinit|test_ocean_zinit|collect_ocean_zinit_tests|ocean_zinit|ocean"
    "test_ocean_data_input|test_ocean_data_input|collect_ocean_data_input_tests|ocean_data_input|ocean"
    "test_ocean_data_forcing|test_ocean_data_forcing|collect_ocean_data_forcing_tests|ocean_data_forcing|ocean"
    "test_ocean_epbl|test_ocean_epbl|collect_ocean_epbl_tests|ocean_epbl|ocean"
    "test_ocean_wave_speed|test_ocean_wave_speed|collect_ocean_wave_speed_tests|ocean_wave_speed|ocean"
    "test_ocean_foxkemper|test_ocean_foxkemper|collect_ocean_foxkemper_tests|ocean_foxkemper|ocean"
    "test_ocean_kappa_shear|test_ocean_kappa_shear|collect_ocean_kappa_shear_tests|ocean_kappa_shear|ocean"
    "test_ocean_isopycnal_slopes|test_ocean_isopycnal_slopes|collect_ocean_isopycnal_slopes_tests|ocean_isopycnal_slopes|ocean"
    "test_ocean_gm|test_ocean_gm|collect_ocean_gm_tests|ocean_gm|ocean"
    "test_ocean_redi|test_ocean_redi|collect_ocean_redi_tests|ocean_redi|ocean"
    "test_ocean_varmix|test_ocean_varmix|collect_ocean_varmix_tests|ocean_varmix|ocean"
    "test_ocean_meke|test_ocean_meke|collect_ocean_meke_tests|ocean_meke|ocean"
    "test_ocean_meke_backscatter|test_ocean_meke_backscatter|collect_ocean_meke_backscatter_tests|ocean_meke_backscatter|ocean"
    "test_ocean_tidal_mixing|test_ocean_tidal_mixing|collect_ocean_tidal_mixing_tests|ocean_tidal_mixing|ocean"
    "test_ocean_convection|test_ocean_convection|collect_ocean_convection_tests|ocean_convection|ocean"
    "test_ocean_double_diffusion|test_ocean_double_diffusion|collect_ocean_double_diffusion_tests|ocean_double_diffusion|ocean"
    "test_ocean_tides_astronomy|test_ocean_tides_astronomy|collect_ocean_tides_astronomy_tests|ocean_tides_astronomy|ocean"
    "test_ocean_tides_disabled_bitident|test_ocean_tides_disabled_bitident|collect_ocean_tides_disabled_bitident_tests|ocean_tides_disabled_bitident|ocean"
    "test_ocean_tidal_forcing|test_ocean_tidal_forcing|collect_ocean_tidal_forcing_tests|ocean_tidal_forcing|ocean"
    "test_ocean_tides_sal|test_ocean_tides_sal|collect_ocean_tides_sal_tests|ocean_tides_sal|ocean"
    "test_ocean_p_surf|test_ocean_p_surf|collect_ocean_p_surf_tests|ocean_p_surf|ocean"
    "test_ocean_p_surf_bitident|test_ocean_p_surf_bitident|collect_ocean_p_surf_bitident_tests|ocean_p_surf_bitident|ocean"
    "test_ocean_obc_tide_nodal|test_ocean_obc_tide_nodal|collect_ocean_obc_tide_nodal_tests|ocean_obc_tide_nodal|ocean"
    "test_ocean_massless|test_ocean_massless|collect_ocean_massless_tests|ocean_massless|ocean"
    "test_ocean_vmix_assembly|test_ocean_vmix_assembly|collect_ocean_vmix_assembly_tests|ocean_vmix_assembly|ocean"
    "test_ocean_bkgnd_mixing|test_ocean_bkgnd_mixing|collect_ocean_bkgnd_mixing_tests|ocean_bkgnd_mixing|ocean"
    "test_ocean_barotropic_substep|test_ocean_barotropic_substep|collect_ocean_barotropic_substep_tests|ocean_barotropic_substep|ocean"
    "test_ocean_bt_cor_ref|test_ocean_bt_cor_ref|collect_ocean_bt_cor_ref_tests|ocean_bt_cor_ref|ocean"
    "test_ocean_bt_substep_zeta_ke|test_ocean_bt_substep_zeta_ke|collect_bt_substep_zeta_ke_tests|ocean_bt_substep_zeta_ke|ocean"
    "test_ocean_pred_corr|test_ocean_pred_corr|collect_ocean_pred_corr_tests|ocean_pred_corr|ocean"
    "test_ocean_cor_ref_seiche|test_ocean_cor_ref_seiche|collect_ocean_cor_ref_seiche_tests|ocean_cor_ref_seiche|ocean"
    "test_ocean_bt_cfl|test_ocean_bt_cfl|collect_ocean_bt_cfl_tests|ocean_bt_cfl|ocean"
    "test_ocean_bt_corrector_hweight|test_ocean_bt_corrector_hweight|collect_bt_corrector_hweight_tests|ocean_bt_corrector_hweight|ocean"
    "test_ocean_bt_cont_type|test_ocean_bt_cont_type|collect_ocean_bt_cont_type_tests|ocean_bt_cont_type|ocean"
    "test_ocean_bt_upstream_h_face|test_ocean_bt_upstream_h_face|collect_ocean_bt_upstream_h_face_tests|ocean_bt_upstream_h_face|ocean"
    "test_ocean_dyn_split|test_ocean_dyn_split|collect_ocean_dyn_split_tests|ocean_dyn_split|ocean"
    "test_ocean_sponge|test_ocean_sponge|collect_ocean_sponge_tests|ocean_sponge|ocean"
    "test_ocean_surface_flux|test_ocean_surface_flux|collect_ocean_surface_flux_tests|ocean_surface_flux|ocean"
    "test_ocean_surface_forcing_type|test_ocean_surface_forcing_type|collect_ocean_surface_forcing_type_tests|ocean_surface_forcing_type|ocean"
    "test_ocean_forcing_rho_ref|test_ocean_forcing_rho_ref|collect_ocean_forcing_rho_ref_tests|ocean_forcing_rho_ref|ocean"
    "test_ocean_frazil|test_ocean_frazil|collect_ocean_frazil_tests|ocean_frazil|ocean"
    "test_ocean_ice_enthalpy|test_ocean_ice_enthalpy|collect_ocean_ice_enthalpy_tests|ocean_ice_enthalpy|ocean"
    "test_ocean_ice_column|test_ocean_ice_column|collect_ocean_ice_column_tests|ocean_ice_column|ocean"
    "test_ocean_ice_coupling|test_ocean_ice_coupling|collect_ocean_ice_coupling_tests|ocean_ice_coupling|ocean"
    "test_ocean_ice_driver_column|test_ocean_ice_driver_column|collect_ocean_ice_driver_column_tests|ocean_ice_driver_column|ocean"
    "test_ocean_ice_snowfall|test_ocean_ice_snowfall|collect_ocean_ice_snowfall_tests|ocean_ice_snowfall|ocean"
    "test_ocean_ice_itd|test_ocean_ice_itd|collect_ocean_ice_itd_tests|ocean_ice_itd|ocean"
    "test_ocean_ice_transport|test_ocean_ice_transport|collect_ocean_ice_transport_tests|ocean_ice_transport|ocean"
    "test_ocean_ice_evp|test_ocean_ice_evp|collect_ocean_ice_evp_tests|ocean_ice_evp|ocean"
    "test_ocean_ice_restart|test_ocean_ice_restart|collect_ocean_ice_restart_tests|ocean_ice_restart|ocean"
    "test_ocean_ice_diags|test_ocean_ice_diags|collect_ocean_ice_diags_tests|ocean_ice_diags|ocean"
    "test_ocean_ice_init|test_ocean_ice_init|collect_ocean_ice_init_tests|ocean_ice_init|ocean"
    "test_ocean_ice_stress_mag|test_ocean_ice_stress_mag|collect_ocean_ice_stress_mag_tests|ocean_ice_stress_mag|ocean"
    "test_ocean_sw_penetration|test_ocean_sw_penetration|collect_ocean_sw_penetration_tests|ocean_sw_penetration|ocean"
    "test_ocean_sw_bl_coupling|test_ocean_sw_bl_coupling|collect_ocean_sw_bl_coupling_tests|ocean_sw_bl_coupling|ocean"
    "test_ocean_restore|test_ocean_restore|collect_ocean_restore_tests|ocean_restore|ocean"
    "test_ocean_geothermal|test_ocean_geothermal|collect_ocean_geothermal_tests|ocean_geothermal|ocean"
    "test_ocean_kpp_convective|test_ocean_kpp_convective|collect_ocean_kpp_convective_tests|ocean_kpp_convective|ocean"
    "test_ocean_kpp_nonlocal|test_ocean_kpp_nonlocal|collect_ocean_kpp_nonlocal_tests|ocean_kpp_nonlocal|ocean"
    "test_ocean_validation|test_ocean_validation|collect_ocean_validation_tests|ocean_validation|ocean"
    "test_ocean_vcoord|test_ocean_vcoord|collect_ocean_vcoord_tests|ocean_vcoord|ocean"
    "test_ocean_vcoord_rho|test_ocean_vcoord_rho|collect_ocean_vcoord_rho_tests|ocean_vcoord_rho|ocean"
    "test_ocean_vcoord_hycom|test_ocean_vcoord_hycom|collect_ocean_vcoord_hycom_tests|ocean_vcoord_hycom|ocean"
    "test_vcoord_parse|test_vcoord_parse|collect_vcoord_parse_tests|vcoord_parse|ocean"
    "test_ocean_remap|test_ocean_remap|collect_ocean_remap_tests|ocean_remap|ocean"
    "test_ocean_regrid_refine|test_ocean_regrid_refine|collect_ocean_regrid_refine_tests|ocean_regrid_refine|ocean"
    "test_ocean_remap_e2e|test_ocean_remap_e2e|collect_ocean_remap_e2e_tests|ocean_remap_e2e|ocean"
    "test_ocean_ppm_h4_remap|test_ocean_ppm_h4_remap|collect_ocean_ppm_h4_remap_tests|ocean_ppm_h4_remap|ocean+core"
    "test_safe_math|test_safe_math|collect_safe_math_tests|safe_math|core"
    "test_scratch_3d_device|test_scratch_3d_device|collect_scratch_3d_device_tests|scratch_3d_device|core"
    "test_remap|test_remap|collect_remap_tests|remap|core"
    "test_remap_pqm|test_remap_pqm|collect_remap_pqm_tests|remap_pqm|core"
    "test_ocean_tracer_adv|test_ocean_tracer_adv|collect_ocean_tracer_adv_tests|ocean_tracer_adv|ocean"
    "test_ocean_conservation_salt_heat|test_ocean_conservation_salt_heat|collect_ocean_conservation_salt_heat_tests|ocean_conservation_salt_heat|ocean"
    "test_ocean_min_thickness|test_ocean_min_thickness|collect_ocean_min_thickness_tests|ocean_min_thickness|ocean"
    "test_ocean_continuity_positive|test_ocean_continuity_positive|collect_ocean_continuity_positive_tests|ocean_continuity_positive|ocean"
    "test_ocean_dyn_multilayer|test_ocean_dyn_multilayer|collect_ocean_dyn_multilayer_tests|ocean_dyn_multilayer|ocean"
    "test_ocean_hvisc|test_ocean_hvisc|collect_ocean_hvisc_tests|ocean_hvisc|ocean"
    "test_ocean_hvisc_resoln|test_ocean_hvisc_resoln|collect_ocean_hvisc_resoln_tests|ocean_hvisc_resoln|ocean"
    "test_ocean_leith|test_ocean_leith|collect_ocean_leith_tests|ocean_leith|ocean"
    "test_ocean_smag_ah|test_ocean_smag_ah|collect_ocean_smag_ah_tests|ocean_smag_ah|ocean"
    "test_ocean_hvisc_bih_cfl|test_ocean_hvisc_bih_cfl|collect_ocean_hvisc_bih_cfl_tests|ocean_hvisc_bih_cfl|ocean"
    "test_ocean_hvisc_kh_bound|test_ocean_hvisc_kh_bound|collect_ocean_hvisc_kh_bound_tests|ocean_hvisc_kh_bound|ocean"
    "test_ocean_hvisc_leith_biharm|test_ocean_hvisc_leith_biharm|collect_ocean_hvisc_leith_biharm_tests|ocean_hvisc_leith_biharm|ocean"
    "test_ocean_hvisc_aniso|test_ocean_hvisc_aniso|collect_ocean_hvisc_aniso_tests|ocean_hvisc_aniso|ocean"
    "test_ocean_hvisc_stress_biharm|test_ocean_hvisc_stress_biharm|collect_ocean_hvisc_stress_biharm_tests|ocean_hvisc_stress_biharm|ocean"
    "test_ocean_ideal_age|test_ocean_ideal_age|collect_ocean_ideal_age_tests|ocean_ideal_age|ocean"
    "test_ocean_fail_loud_dispatch|test_ocean_fail_loud_dispatch|collect_ocean_fail_loud_dispatch_tests|ocean_fail_loud_dispatch|ocean"
    "test_ocean_pseudo_salt|test_ocean_pseudo_salt|collect_ocean_pseudo_salt_tests|ocean_pseudo_salt|ocean"
    "test_ocean_eady_ic|test_ocean_eady_ic|collect_ocean_eady_ic_tests|ocean_eady_ic|ocean"
    "test_ocean_geostrophic_adjustment_ic|test_ocean_geostrophic_adjustment_ic|collect_ocean_geostrophic_adjustment_ic_tests|ocean_geostrophic_adjustment_ic|ocean"
    "test_ocean_diag|test_ocean_diag|collect_ocean_diag_tests|ocean_diag|ocean"
    "test_ocean_diag_reduce|test_ocean_diag_reduce|collect_ocean_diag_reduce_tests|ocean_diag_reduce|ocean"
    "test_ocean_porous|test_ocean_porous|collect_ocean_porous_tests|ocean_porous|ocean"
    "test_ocean_diag_remap|test_ocean_diag_remap|collect_ocean_diag_remap_tests|ocean_diag_remap|ocean"
    "test_ocean_diag_netcdf|test_ocean_diag_netcdf|collect_ocean_diag_netcdf_tests|ocean_diag_netcdf|ocean"
    "test_ocean_diag_precision|test_ocean_diag_precision|collect_ocean_diag_precision_tests|ocean_diag_precision|ocean"
    "test_ocean_boundary|test_ocean_boundary|collect_ocean_boundary_tests|ocean_boundary|ocean"
    "test_ocean_periodic|test_ocean_periodic|collect_ocean_periodic_tests|ocean_periodic|ocean"
    "test_ocean_fold|test_ocean_fold|collect_ocean_fold_tests|ocean_fold|ocean"
    "test_ocean_tripolar|test_ocean_tripolar|collect_ocean_tripolar_tests|ocean_tripolar|ocean"
    "test_ocean_stability_audit|test_ocean_stability_audit|collect_ocean_stability_audit_tests|ocean_stability_audit|ocean"
    "test_ocean_obc_baroclinic|test_ocean_obc_baroclinic|collect_ocean_obc_baroclinic_tests|ocean_obc_baroclinic|ocean"
    "test_ocean_obc_eta_ghost|test_ocean_obc_eta_ghost|collect_ocean_obc_eta_ghost_tests|ocean_obc_eta_ghost|ocean"
    "test_ocean_driver_diag|test_ocean_driver_diag|collect_ocean_driver_diag_tests|ocean_driver_diag|ocean"
    "test_ocean_restart|test_ocean_restart|collect_ocean_restart_tests|ocean_restart|ocean"
    "test_ocean_double_gyre_helpers|test_ocean_double_gyre_helpers|collect_ocean_double_gyre_helpers_tests|ocean_double_gyre_helpers|ocean"
    "test_ocean_analytical|test_ocean_analytical|collect_ocean_analytical_tests|ocean_analytical|ocean"
    # Split per-subtest (same source) so the heavy 20000-step geostrophic run is
    # its own ctest binary, focusable/skippable in isolation.
    "test_ocean_baroclinic_geostrophic|test_ocean_baroclinic_longrun|collect_ocean_baroclinic_geostrophic_tests|ocean_baroclinic_geostrophic|ocean"
    "test_ocean_baroclinic_stommel|test_ocean_baroclinic_longrun|collect_ocean_baroclinic_stommel_tests|ocean_baroclinic_stommel|ocean"
    "test_ocean_baroclinic_budget|test_ocean_baroclinic_longrun|collect_ocean_baroclinic_budget_tests|ocean_baroclinic_budget|ocean"
    "test_ocean_bottom_drag|test_ocean_bottom_drag|collect_ocean_bottom_drag_tests|ocean_bottom_drag|ocean"
    "test_ppm_limit_pos|test_ppm_limit_pos|collect_ppm_limit_pos_tests|ppm_limit_pos|ocean"
    "test_ocean_bt_corrector_bc_pgf|test_ocean_bt_corrector_bc_pgf|collect_bt_corrector_bc_pgf_tests|ocean_bt_corrector_bc_pgf|ocean"
    "test_ocean_bt_substep_drag|test_ocean_bt_substep_drag|collect_bt_substep_drag_tests|ocean_bt_substep_drag|ocean"
    "test_ocean_wave_drag|test_ocean_wave_drag|collect_ocean_wave_drag_tests|ocean_wave_drag|ocean"
    "test_ocean_wetdry|test_ocean_wetdry|collect_ocean_wetdry_tests|ocean_wetdry|ocean"
    "test_ocean_wetdry_driver|test_ocean_wetdry_driver|collect_ocean_wetdry_driver_tests|ocean_wetdry_driver|ocean"
    "test_ocean_surface_stress|test_ocean_surface_stress|collect_ocean_surface_stress_tests|ocean_surface_stress|ocean"
    "test_coriolis_inertial|test_coriolis_inertial|collect_coriolis_inertial_tests|coriolis_inertial|ocean"
    "test_coriolis_multilayer|test_coriolis_multilayer|collect_coriolis_multilayer_tests|coriolis_multilayer|ocean"
    "test_ocean_dyn_step|test_ocean_dyn_step|collect_ocean_dyn_step_tests|ocean_dyn_step|ocean"
    "test_ocean_cfl_trunc|test_ocean_cfl_trunc|collect_ocean_cfl_trunc_tests|ocean_cfl_trunc|ocean"
    "test_ocean_chksum|test_ocean_chksum|collect_ocean_chksum_tests|ocean_chksum|ocean"
    "test_dc_reduce|test_dc_reduce|collect_dc_reduce_tests|dc_reduce|core"
    "test_dc_transfer_sentinel|test_dc_transfer_sentinel|collect_dc_transfer_sentinel_tests|dc_transfer_sentinel|core"
    "test_ocean_accel_visc_rem|test_ocean_accel_visc_rem|collect_ocean_accel_visc_rem_tests|ocean_accel_visc_rem|ocean"
    "test_decomp|test_decomp|collect_decomp_tests|decomp|core"
    "test_driver_ocean|test_driver_ocean|collect_driver_ocean_tests|driver_ocean|ocean"
    # Per-test isolation for the rdb_api tests — each create/destroy of a solver
    # handle leaves residual state in the OpenMP-target present table on Aurora
    # ifx that poisons the next test's enter_data. One ctest binary per scenario
    # keeps each in a fresh process. The ${MODULE_NAME}.F90 substitution below
    # lets four binaries share tests/test_rdb_api.F90 as their source.
    "test_vcoord|test_vcoord|collect_vcoord_tests|vcoord|core"
    "test_profiler|test_profiler|collect_profiler_tests|profiler|core"
    "test_vcoord_target|test_vcoord_target|collect_vcoord_target_tests|vcoord_target|core"
    "test_vcoord_zstar_full|test_vcoord_zstar_full|collect_vcoord_zstar_full_tests|vcoord_zstar_full|core"
    # Per-test isolation — module-scope kernel workspaces (PPM/hdiff scratch)
    # accumulate state across these four tests and the host-side allocate trips
    # "array already allocated" on the second test. One ctest binary per
    # scenario sidesteps it.
    "test_decomp_extras|test_decomp_extras|collect_decomp_extras_tests|decomp_extras|core"
    # Per-test isolation — same workspace-pollution reason as above.
    "test_ocean_halo_single_rank|test_ocean_halo_single_rank|collect_ocean_halo_single_rank_tests|ocean_halo_single_rank|ocean"
    # Ocean C ABI lifecycle skeleton (Python runtime API plan P1). Own binary:
    # the single-live-handle module state in rdb_ocean_api (g_handle_live) must
    # start from a fresh process, same reasoning as every other
    # handle/module-state-isolated suite in this list.
    "test_rdb_ocean_api|test_rdb_ocean_api|collect_rdb_ocean_api_tests|rdb_ocean_api|ocean"
    "test_ocean_api_p2|test_ocean_api_p2|collect_ocean_api_p2_tests|ocean_api_p2|ocean"
    # P2.5 pre-create geometry injection. Own binary: same single-live-handle
    # module-state reasoning as test_rdb_ocean_api/test_ocean_api_p2.
    "test_ocean_geometry_inject|test_ocean_geometry_inject|collect_ocean_geometry_inject_tests|ocean_geometry_inject|ocean"
)

# Tests that MUST link the STATIC core.
#
# The shared core (.so) is the default test link (RDB_TESTS_LINK_SHARED, on by
# default) because it skips the per-exe static-archive pull + RDC device link —
# the slow part of a GPU test build.
#
# A test needs the static core ONLY when its own `do concurrent` (a device
# kernel) CALLS a core `!$acc routine seq`: that call is a device→device
# reference that can only be resolved by device-linking WITHIN one binary, never
# across a .so boundary. A self-contained `do concurrent`, or a host-level call
# to an acc-seq routine, does NOT trigger it — those link fine against the .so.
#
# Empirically that is just test_remap_pqm (its on-device `do concurrent` run
# calls remap_column_pqm). This list is the authority: if a GPU build's device
# link fails with an unresolved device symbol for some test, add it here. On a
# non-GPU (gfortran/ifx) build the distinction is moot — nothing generates
# device code — so a stale entry only costs that test the shared-link speedup.
set(RDB_STATIC_CORE_TESTS test_remap_pqm)

# Tests that pull in NetCDF/HDF5 — must run single-threaded because HDF5 is not
# built thread-safe (OMP_NUM_THREADS=1 prevents crashes when the OpenMP runtime
# is active via -mp=multicore or -fopenmp).
set(NETCDF_TESTS
    test_ocean_diag_netcdf
    test_ocean_diag_precision
    test_ocean_supergrid
    test_ocean_status_returns
    test_ocean_zinit
    test_ocean_data_input
    test_ocean_data_forcing
    test_ocean_driver_diag
    test_ocean_restart
    test_ocean_ice_evp
    test_ocean_ice_restart
    test_driver_ocean)

# Tests that directly `use` modules from the NetCDF I/O subsystem (rdb_forcing,
# rdb_gauges, rdb_output, ...).  Skipped when -DRDB_ENABLE_NETCDF=OFF because
# those modules aren't built.
set(NETCDF_REQUIRED_TESTS
    test_ocean_diag_netcdf
    test_ocean_diag_precision
    test_ocean_supergrid
    # `use rdb_driver`, which src/CMakeLists.txt only compiles with NetCDF on
    test_ocean_diag_remap
    # one case writes a restart, and ocean_state_restart_write is fail-loud
    # without NetCDF
    test_ocean_ice_init
    test_ocean_zinit
    test_ocean_data_input
    test_ocean_data_forcing
    test_ocean_driver_diag
    test_ocean_restart
    test_ocean_ice_restart
    test_driver_ocean)

# Template-generated tests that drive the real comm-env lifecycle (driver_run ->
# comm_env_init).  In MPI builds, launch these under the MPI runner with 1 rank
# so MPI starts up cleanly under any implementation (a bare singleton run is not
# universally supported); in serial builds they run directly.  The comm-env is
# finalised by the shared per-test main.
set(MPI_LAUNCH_TESTS test_driver_ocean test_ocean_halo_single_rank
                     test_ocean_console_stats_efp)

# ---------------------------------------------------------------------------
# LFortran 0.64 (alpha) known failures.  Two categories, handled differently:
#
# (A) UNBUILDABLE — do not COMPILE under lfortran, so they must be excluded from
# the build (otherwise `cmake --build` itself fails; they can't be "tests that
# fail"): test_config   - LLVM-backend codegen ICE ("x is nullptr") when a unit
# consumes the .mod of a module that (transitively) calls rdb_config's
# read_config_from_string. test_profiler - VOLATILE attribute statement,
# rejected by lfortran 0.64 ("Attribute declaration not supported").
#
# (B) RUNTIME FAILURES — BUILD fine but crash when run.  Deliberately LEFT IN
# the suite as live regression markers: they will start passing once lfortran
# fixes the bug, so re-running the suite after an upgrade tells us for free.
# test_nml_schema, test_config_schema - SIGSEGV on allocate(dst, source=src) of
# a type with a character- pointer component (the nml_string/nml_enum key
# `tgt`).
#
# MREs for every bug live in tmp_local_artifacts/lf_nz_mre/; track upstream at
# lfortran/lfortran.  Extend either list as more surface.
# ---------------------------------------------------------------------------
if(CMAKE_Fortran_COMPILER_ID STREQUAL "LFortran")
  set(RDB_LFORTRAN_UNBUILDABLE test_config test_profiler)
  set(RDB_LFORTRAN_RUNTIME_FAILURES test_nml_schema test_config_schema)

  # (A) remove the unbuildable ones from the build
  foreach(_row ${RDB_TESTS})
    string(REPLACE "|" ";" _fields "${_row}")
    list(GET _fields 0 _tname)
    list(FIND RDB_LFORTRAN_UNBUILDABLE "${_tname}" _skip_idx)
    if(NOT _skip_idx EQUAL -1)
      list(REMOVE_ITEM RDB_TESTS "${_row}")
    endif()
  endforeach()
  # (B) runtime failures stay in RDB_TESTS on purpose (regression markers).

  message(STATUS "")
  message(
    STATUS "LFortran 0.64 known failures (re-check after any lfortran upgrade;")
  message(STATUS "  MREs in tmp_local_artifacts/lf_nz_mre/):")
  message(STATUS "    excluded from build (do not compile):")
  foreach(_t ${RDB_LFORTRAN_UNBUILDABLE})
    message(STATUS "      ${_t}")
  endforeach()
  message(STATUS "    built, expected to FAIL at runtime (kept as markers):")
  foreach(_t ${RDB_LFORTRAN_RUNTIME_FAILURES})
    message(STATUS "      ${_t}")
  endforeach()
  message(STATUS "")
endif()

# Valid values for the <regime> field of a RDB_TESTS row (see below).
# `ocean+core` expands to BOTH labels (cross-regime tests that stay closed-side
# but also exercise shared infrastructure).
set(RDB_TEST_REGIMES ocean core ocean+core)

# Build per-test executables from the shared template.
foreach(row ${RDB_TESTS})
  string(REPLACE "|" ";" fields "${row}")
  list(LENGTH fields n_fields)
  if(NOT n_fields EQUAL 5)
    message(
      FATAL_ERROR
        "RDB_TESTS row '${row}' has ${n_fields} fields; expected 5: "
        "<test-name>|<module-name>|<collect-fn>|<suite-display>|<regime> "
        "with <regime> one of: ${RDB_TEST_REGIMES}")
  endif()
  list(GET fields 0 test_name)
  list(GET fields 1 MODULE_NAME)
  list(GET fields 2 COLLECT_FN)
  list(GET fields 3 SUITE_NAME)
  list(GET fields 4 TEST_REGIME)
  if(NOT TEST_REGIME IN_LIST RDB_TEST_REGIMES)
    message(
      FATAL_ERROR "RDB_TESTS row '${row}': invalid regime '${TEST_REGIME}' "
                  "(must be one of: ${RDB_TEST_REGIMES})")
  endif()
  # `+`-separated regimes become multiple ctest LABELS.
  string(REPLACE "+" ";" TEST_REGIME_LABELS "${TEST_REGIME}")

  if(NOT RDB_ENABLE_NETCDF AND ${test_name} IN_LIST NETCDF_REQUIRED_TESTS)
    continue()
  endif()

  set(TESTER_NAME "tester_${test_name}")
  set(driver_path "${CMAKE_CURRENT_BINARY_DIR}/generated/${TESTER_NAME}.F90")
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/tester_template.F90.in
                 ${driver_path} @ONLY)

  # Source file is named after the MODULE, not the binary -- this lets multiple
  # per-test binaries share one source file when present-table isolation
  # requires running each test in its own process (e.g. OpenMP target offload
  # with module-scope allocatables). For existing rows where binary == module
  # name this is a no-op.
  add_executable(${test_name} ${CMAKE_CURRENT_SOURCE_DIR}/${MODULE_NAME}.F90
                              ${driver_path})
  # Per-target module dir so multiple binaries sharing one source (rdb_api
  # split, zstar_full_tracer_conservation split, etc.) don't race to write the
  # same .mod file under parallel make. Each target gets its own scratch dir;
  # the shared modules/ dir is still on the include path so the test source can
  # use core_rdb / analytical.
  set_target_properties(
    ${test_name} PROPERTIES Fortran_MODULE_DIRECTORY
                            "${CMAKE_CURRENT_BINARY_DIR}/modules/${test_name}")
  target_include_directories(
    ${test_name} PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/modules/${test_name}"
                         "${CMAKE_BINARY_DIR}/modules")
  # Core selection: the SHARED core is the default fast link; only tests in the
  # RDB_STATIC_CORE_TESTS authority list (device code in their own source) fall
  # back to the static core. When the shared core wasn't built
  # (RDB_TESTS_LINK_SHARED=OFF) every test links static.
  set(_core_for_test ${main_lib})
  if(TARGET ${main_lib}_shared AND NOT (${test_name} IN_LIST
                                        RDB_STATIC_CORE_TESTS))
    set(_core_for_test ${main_lib}_shared)
  endif()
  target_link_libraries(${test_name} PRIVATE ${_core_for_test} rdb_analytical
                                             test-drive::test-drive)
  if(RDB_ENABLE_MPI AND ${test_name} IN_LIST MPI_LAUNCH_TESTS)
    add_test(NAME rdb_${test_name}
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1
                     $<TARGET_FILE:${test_name}>)
  else()
    add_test(NAME rdb_${test_name} COMMAND ${test_name})
  endif()
  set_tests_properties(rdb_${test_name} PROPERTIES LABELS
                                                   "${TEST_REGIME_LABELS}")
  if(${test_name} IN_LIST NETCDF_TESTS)
    set_tests_properties(rdb_${test_name} PROPERTIES ENVIRONMENT
                                                     "OMP_NUM_THREADS=1")
  endif()
  # SENTINEL: nvfortran 26.5 miscompiles `transfer()` inside `do concurrent`
  # (~95% of elements wrong, a shifted-index bug; the sentinel source documents
  # it in full). The test asserts the CORRECT behaviour, so it fails there
  # today; invert the verdict for NVHPC only, so the suite is quiet while the
  # bug exists AND goes red the moment NVIDIA fixes it.  When that happens: drop
  # this block, drop the `dc-transfer` pre-commit lint, and rdb_ocean_chksum can
  # move to `do concurrent`.  Every other compiler runs it as an ordinary test.
  #
  # Scoped to builds that actually get `-stdpar`.  The miscompile needs `do
  # concurrent` to be LOWERED to parallel code, which on NVHPC means
  # `-stdpar=gpu` (RDB_ENABLE_GPU) or `-stdpar=multicore` (RDB_ENABLE_THREADS).
  # A plain serial NVHPC build gets neither -- just `-Mfree -Mbackslash` -- so
  # TRANSFER is correct there, the sentinel passes, and an unconditional
  # WILL_FAIL would INVERT that into a red. Verified on nvfortran 26.5: serial
  # build, sentinel exits 0.
  if(${test_name} STREQUAL "test_dc_transfer_sentinel"
     AND CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC"
     AND (RDB_ENABLE_GPU OR RDB_ENABLE_THREADS))
    set_tests_properties(rdb_${test_name} PROPERTIES WILL_FAIL TRUE)
  endif()
endforeach()

# MPI integration tests (only when MPI is enabled). Per-test rank counts because
# the structured halo test wants a 2x2 grid (4 ranks) while the unstructured
# halo and lake-at-rest tests run on 2 ranks.
if(RDB_ENABLE_MPI)
  set(MPI_TESTS test_halo_mpi test_halo_ocean_mpi test_ocean_dyn_mpi
                test_ocean_bt_cfl_mpi test_efp_mpi)
  set(MPI_TEST_NPROCS_test_halo_mpi 4)
  set(MPI_TEST_NPROCS_test_halo_ocean_mpi 4)
  set(MPI_TEST_NPROCS_test_ocean_dyn_mpi 2)
  set(MPI_TEST_NPROCS_test_ocean_bt_cfl_mpi 2)
  set(MPI_TEST_NPROCS_test_efp_mpi 4)
  # Per-test regime label (same values as the RDB_TESTS 5th field;
  # classification is by the source path of the modules each test imports).
  set(MPI_TEST_REGIME_test_halo_mpi core)
  set(MPI_TEST_REGIME_test_halo_ocean_mpi ocean)
  set(MPI_TEST_REGIME_test_ocean_dyn_mpi ocean)
  set(MPI_TEST_REGIME_test_ocean_bt_cfl_mpi ocean)
  set(MPI_TEST_REGIME_test_efp_mpi core)
  set(MPI_NETCDF_TESTS)
  foreach(test ${MPI_TESTS})
    add_executable(${test} mpi/${test}.F90)
    target_include_directories(${test} PRIVATE "${CMAKE_BINARY_DIR}/modules")
    target_link_libraries(${test} PRIVATE ${main_lib} test-drive::test-drive)
    add_test(NAME rdb_${test}
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG}
                     ${MPI_TEST_NPROCS_${test}} $<TARGET_FILE:${test}>)
    if(NOT DEFINED MPI_TEST_REGIME_${test} OR NOT MPI_TEST_REGIME_${test}
                                              IN_LIST RDB_TEST_REGIMES)
      message(
        FATAL_ERROR "MPI test '${test}' has no valid MPI_TEST_REGIME_${test} "
                    "(must be one of: ${RDB_TEST_REGIMES})")
    endif()
    set_tests_properties(rdb_${test} PROPERTIES LABELS
                                                "${MPI_TEST_REGIME_${test}}")
    if(${test} IN_LIST MPI_NETCDF_TESTS)
      set_tests_properties(rdb_${test} PROPERTIES ENVIRONMENT
                                                  "OMP_NUM_THREADS=1")
    endif()
  endforeach()
  # test_ocean_dyn_mpi uses ocean_test_metrics from rdb_analytical.
  target_link_libraries(test_ocean_dyn_mpi PRIVATE rdb_analytical)
  # test_ocean_bt_cfl_mpi uses make_spherical_metrics from rdb_analytical.
  target_link_libraries(test_ocean_bt_cfl_mpi PRIVATE rdb_analytical)
  # Extra rank-count legs for the ocean halo test (same binary, the test
  # branches on nprocs): 2-rank periodic-x wrap-vs-reference bit-exactness and
  # the 1-rank D0 gate (halo == wrap kernel / non-periodic no-op).
  add_test(NAME rdb_test_halo_ocean_mpi_2rank
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                   $<TARGET_FILE:test_halo_ocean_mpi>)
  add_test(NAME rdb_test_halo_ocean_mpi_1rank
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1
                   $<TARGET_FILE:test_halo_ocean_mpi>)
  set_tests_properties(rdb_test_halo_ocean_mpi_2rank
                       rdb_test_halo_ocean_mpi_1rank PROPERTIES LABELS "ocean")
  # O2 gate: 2-rank and 1-rank legs for the dyn-core MPI test. The 2-rank leg is
  # registered by the foreach above (MPI_TEST_NPROCS=2). Add the 1-rank leg here
  # so both are in the test suite.
  add_test(NAME rdb_test_ocean_dyn_mpi_1rank
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1
                   $<TARGET_FILE:test_ocean_dyn_mpi>)
  # O4: 4-rank leg — px=4 x-split (the test decomposes over the launched rank
  # count); covers multi-seam (two interior seams per row) + the 4-rank
  # periodic-x ring on the same three physics legs.
  add_test(NAME rdb_test_ocean_dyn_mpi_4rank
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
                   $<TARGET_FILE:test_ocean_dyn_mpi>)
  set_tests_properties(rdb_test_ocean_dyn_mpi_1rank rdb_test_ocean_dyn_mpi_4rank
                       PROPERTIES LABELS "ocean")
  # Spherical meridional-split n_inner consistency regression (the
  # MPI_ERR_TRUNCATE bug).  The 2-rank leg is registered by the foreach above;
  # add a 1-rank (identity-reduction) and 4-rank (stronger latitude spread) leg.
  add_test(NAME rdb_test_ocean_bt_cfl_mpi_1rank
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1
                   $<TARGET_FILE:test_ocean_bt_cfl_mpi>)
  add_test(NAME rdb_test_ocean_bt_cfl_mpi_4rank
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
                   $<TARGET_FILE:test_ocean_bt_cfl_mpi>)
  set_tests_properties(
    rdb_test_ocean_bt_cfl_mpi_1rank rdb_test_ocean_bt_cfl_mpi_4rank
    PROPERTIES LABELS "ocean")
endif()

# --------------------------------------------------------------------------
# Python package tests (P3, python/tests/) -- wired into ctest.
#
# Two independent prerequisites, NEITHER guaranteed on a given tree, so this
# degrades to a clear configure-time SKIP (no test registered, no build failure)
# rather than failing configure or making the Fortran build depend on Python: 1.
# A built shared core (RDB_BUILD_SHARED=ON -> librdb_core.so), which is what
# python/rdb/_ffi.py dlopens. 2. A Python 3 interpreter that ALREADY has pytest
# importable -- this repo never runs `pip install`, so "find one" means "find
# one that already has it", not "install it". RDB_PYTHON overrides the search
# with an explicit interpreter path.
# --------------------------------------------------------------------------
set(_rdb_python_candidates "")
if(DEFINED ENV{RDB_PYTHON})
  list(APPEND _rdb_python_candidates "$ENV{RDB_PYTHON}")
endif()
find_package(
  Python3
  COMPONENTS Interpreter
  QUIET)
if(Python3_EXECUTABLE)
  list(APPEND _rdb_python_candidates "${Python3_EXECUTABLE}")
endif()
list(APPEND _rdb_python_candidates "python3")

set(RDB_PYTHON_EXECUTABLE "")
foreach(_cand IN LISTS _rdb_python_candidates)
  if(NOT _cand)
    continue()
  endif()
  execute_process(
    COMMAND "${_cand}" -c "import pytest"
    RESULT_VARIABLE _rdb_pytest_rc
    OUTPUT_QUIET ERROR_QUIET)
  if(_rdb_pytest_rc EQUAL 0)
    set(RDB_PYTHON_EXECUTABLE "${_cand}")
    break()
  endif()
endforeach()

if(NOT RDB_ENABLE_NETCDF)
  # The suite exercises the diagnostics path, and engine_configure_diag is
  # fail-loud on a build without NetCDF. Skipping keeps the NetCDF-off
  # portability configuration green; the suite runs in every NetCDF build.
  message(
    STATUS "Python tests (python/tests): SKIPPED -- RDB_ENABLE_NETCDF=OFF")
elseif(NOT TARGET ${main_lib}_shared)
  message(STATUS "Python tests (python/tests): SKIPPED -- "
                 "RDB_BUILD_SHARED is OFF, no librdb_core.so to test "
                 "against. Configure with -DRDB_BUILD_SHARED=ON to " "enable.")
elseif(NOT RDB_PYTHON_EXECUTABLE)
  message(
    STATUS "Python tests (python/tests): SKIPPED -- no Python 3 "
           "interpreter with pytest importable was found (tried: "
           "${_rdb_python_candidates}). Set RDB_PYTHON to an "
           "interpreter that already has pytest (never `pip install` "
           "in this repo) to enable.")
else()
  message(STATUS "Python tests (python/tests): ENABLED, using "
                 "${RDB_PYTHON_EXECUTABLE}")
  add_test(
    NAME rdb_python_tests
    COMMAND "${RDB_PYTHON_EXECUTABLE}" -m pytest
            "${PROJECT_SOURCE_DIR}/python/tests" -q
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}")
  set_tests_properties(
    rdb_python_tests PROPERTIES LABELS "ocean;python" ENVIRONMENT
                                "RDB_LIB=$<TARGET_FILE:${main_lib}_shared>")
endif()
