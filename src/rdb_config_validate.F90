submodule(rdb_config) rdb_config_validate
   !! validate_config and its helper predicates.
   implicit none

contains

   module subroutine validate_config(cfg, ierr)
      !! Validate configuration parameters after reading
      !!
      !! Logs errors for invalid values that would crash the solver,
      !! and warnings for suspicious but non-fatal settings.
      !!
      !! `ierr`, when present, returns `OCEAN_STATUS_ERR_CONFIG_VALIDATE`
      !! instead of `error stop`-ing on the first accumulated failure
      !! (every individual check still logs via `global_logger%error`
      !! unchanged); absent behaves as today.
      use rdb_ocean_lateral_mix, only: parse_lateral_closure, &
                                       lateral_closure_is_implemented, &
                                       lateral_closure_conflicts_smag_ah, &
                                       has_biharmonic_backstop, &
                                       leith_biharm_is_inert
      use rdb_ocean_horizontal_viscosity, only: aniso_mode_is_implemented
      use rdb_vcoord, only: parse_vcoord_type, vcoord_h_min_is_coherent, &
                            parse_z_fixed_profile, z_fixed_nominal_dz, &
                            ZFIXED_PROFILE_INVALID, ZFIXED_PROFILE_UNIFORM, &
                            ZFIXED_PROFILE_LIST, ZFIXED_DZ_OK, ZFIXED_DZ_ERR_COUNT, &
                            ZFIXED_DZ_ERR_TOO_DEEP
      use rdb_constants, only: VCOORD_SIGMA, VCOORD_ZSTAR, VCOORD_EULERIAN_Z, &
                               VCOORD_ZSIGMA, VCOORD_LAGRANGIAN, VCOORD_ZSTAR_SIGMA, &
                               VCOORD_ZSTAR_FULL, VCOORD_Z_FIXED, VCOORD_RHO, VCOORD_HYCOM, &
                               H_VANISHED
      use rdb_ocean_boundary_types, only: ocean_bc_type_from_string, OBC_PERIODIC, OBC_WALL, &
                                          OBC_INVALID, OBC_TRIPOLAR_FOLD
      use rdb_coriolis_adv, only: parse_pv_variant, pv_variant_is_implemented, &
                                  parse_pv_adv_scheme, pv_adv_scheme_is_implemented, &
                                  pv_adv_required_nghost
      use rdb_ocean_bottom_drag, only: parse_bdrag_variant, bdrag_variant_is_implemented
      use rdb_ocean_pressure_force, only: parse_opgf_variant, gprime_nz_is_supported, &
                                          OPGF_VARIANT_FV_MOM6
      use rdb_ocean_tidal_mixing, only: tidal_mixing_is_inert
      use rdb_ocean_pseudo_salt, only: pseudo_salt_conflicts_restore, &
                                       pseudo_salt_conflicts_ice, &
                                       pseudo_salt_needs_thermo_warning
      use rdb_ocean_surface_flux, only: sw_source_is_implemented
      use rdb_ocean_vmix, only: kpp_sw_method_is_implemented, &
                                bkgnd_henyey_conflicts_profile
      use rdb_eos, only: parse_tfreeze_set, TFREEZE_SET_INVALID
      use rdb_ocean_top_drag, only: parse_tdrag_variant, tdrag_variant_is_implemented, &
                                    TDRAG_LINEAR, TDRAG_QUADRATIC
      use rdb_ocean_cavity_melt, only: parse_cavity_exchange_law, parse_cavity_ice_mode, &
                                       CAVITY_LAW_INVALID, CAVITY_LAW_CONST_GAMMA, &
                                       CAVITY_LAW_HJ99, CAVITY_LAW_YUNG25, &
                                       CAVITY_ICE_INVALID, CAVITY_ICE_INSULATING, &
                                       CAVITY_ICE_ADV_DIFF, &
                                       parse_cavity_freshwater, parse_cavity_volume_comp, &
                                       CAVITY_FW_INVALID, CAVITY_FW_VIRTUAL, CAVITY_FW_MASS, &
                                       CAVITY_VC_INVALID, CAVITY_VC_NONE, CAVITY_VC_UNIFORM_OPEN
      use rdb_ocean_cavity, only: parse_cavity_draft_sign, CAVITY_SIGN_INVALID
      use rdb_calendar, only: date_t, parse_date, calendar_name_is_implemented
      type(config_t), intent(in) :: cfg
      integer, intent(out), optional :: ierr
         !! Non-zero on any cross-knob semantic validation failure when
         !! present; absent behaves as today (`error stop`).

      logical :: has_error
      integer :: atm_code
      integer :: lateral_closure_code
      type(date_t) :: start_date_parsed
      logical :: start_date_ok

      has_error = .false.

      ! Simulation regime
      if (trim(cfg%sim_type) /= "ocean") then
         call logger%error("Invalid sim_type = '"//trim(cfg%sim_type)// &
                           "': must be 'ocean' (the coastal regime was split out into "// &
                           "its own repository)")
         has_error = .true.
      end if

      ! Grid parameters
      if (cfg%nx < 1) then
         call logger%error("Invalid nx = "//to_string(cfg%nx)//": must be >= 1")
         has_error = .true.
      end if
      if (cfg%ny < 1) then
         call logger%error("Invalid ny = "//to_string(cfg%ny)//": must be >= 1")
         has_error = .true.
      end if
      if (cfg%nghost < 1) then
         call logger%error("Invalid nghost = "//to_string(cfg%nghost)//": must be >= 1")
         has_error = .true.
      end if
      if (cfg%dx <= 0.0_wp) then
         call logger%error("Invalid dx = "//to_string(cfg%dx)//": must be > 0")
         has_error = .true.
      end if
      if (cfg%dy <= 0.0_wp) then
         call logger%error("Invalid dy = "//to_string(cfg%dy)//": must be > 0")
         has_error = .true.
      end if

      ! Time parameters
      if (cfg%cfl <= 0.0_wp .or. cfg%cfl > 1.0_wp) then
         call logger%error("Invalid cfl = "//to_string(cfg%cfl)//": must be in (0, 1]")
         has_error = .true.
      end if
      if (cfg%t_end <= 0.0_wp) then
         call logger%error("Invalid t_end = "//to_string(cfg%t_end)//": must be > 0")
         has_error = .true.
      end if
      if (.not. calendar_name_is_implemented(cfg%calendar)) then
         call logger%error("Invalid calendar = '"//trim(cfg%calendar)// &
                           "': only gregorian/standard/proleptic_gregorian are "// &
                           "implemented (noleap/360_day are a different day-count rule)")
         has_error = .true.
      end if
      if (len_trim(cfg%start_date) > 0) then
         call parse_date(cfg%start_date, start_date_parsed, start_date_ok)
         if (.not. start_date_ok) then
            call logger%error("Invalid start_date = '"//trim(cfg%start_date)// &
                              "': want 'YYYY-MM-DD' or 'YYYY-MM-DD hh:mm:ss'")
            has_error = .true.
         end if
      end if

      ! Multilayer parameters
      if (cfg%use_multilayer) then
         if (cfg%nz_layers < 1) then
            call logger%error("Invalid nz_layers = "//to_string(cfg%nz_layers)// &
                              ": must be >= 1 for multilayer")
            has_error = .true.
         end if
      end if

      ! Per-column stack-workspace bound.  FAIL LOUD, both regimes.
      !
      ! Every layered kernel (BPG, ALE remap, kappa-shear, Redi, the
      ! diag-remap fills) carries fixed-size `NZ_STACK_MAX` thread-local
      ! column arrays.  Overrunning them corrupts thread-local storage
      ! SILENTLY — wrong answers, no crash — so this must refuse the run,
      ! not warn.  It replaces a warning-only check that additionally
      ! never fired on the ocean path: no ocean namelist sets
      ! `use_multilayer` (it defaults .false.), so the ocean regime had no
      ! nz guard at all beyond the wavespeed- and Redi-specific ones.
      !
      ! The requirement is `nz + 1` — see `nz_stack_required`.
      if (cfg%use_multilayer .or. trim(cfg%sim_type) == "ocean") then
         if (.not. nz_stack_is_sufficient(cfg%nz_layers)) then
            call logger%error("nz_layers = "//to_string(cfg%nz_layers)// &
                              " needs NZ_STACK_MAX >= "// &
                              to_string(nz_stack_required(cfg%nz_layers))// &
                              " but this binary was compiled with NZ_STACK_MAX = "// &
                              to_string(NZ_STACK_MAX)// &
                              ". Per-column stack kernels (BPG, ALE remap, kappa-shear, "// &
                              "Redi, diag remap) would overrun thread-local storage and "// &
                              "silently produce wrong answers. Rebuild with "// &
                              "-DRDB_NZ_STACK_MAX="// &
                              to_string(nz_stack_required(cfg%nz_layers))//" (or larger).")
            has_error = .true.
         end if
      end if

      ! Ocean windowed-drain tracer_recon guard (Q6).  A SEPARATE knob from
      ! the coastal tracer_recon above — it lives in `&ocean_vmix_nml` and
      ! only takes effect on sim_type='ocean' (ignored on coastal, like all
      ! ocean knobs).  Reuses the coastal parse + per-rung nghost helpers
      ! but NOT the coastal support-status predicate (which rejects ocean by
      ! design).  Fail loud on pqm/unknown and on an insufficient nghost for
      ! the requested rung (weno5→3, weno7→4, weno9→5).  The default nghost
      ! is 3 (bumped from 2 — see the grid config default), so ppm/weno5 are
      ! satisfied by default and only weno7/9 need an explicit extra column.
      if (trim(cfg%sim_type) == "ocean") then
         block
            use rdb_recon_weno, only: parse_tracer_recon, &
                                      tracer_recon_required_nghost, &
                                      TRACER_RECON_PPM
            integer :: ocn_recon_code
            ocn_recon_code = parse_tracer_recon(trim(cfg%ocean%vmix%tracer_recon))
            if (ocn_recon_code == -2) then
               call logger%error("&ocean_vmix_nml tracer_recon='pqm' is not implemented; "// &
                                 "use ppm, weno5, weno7 or weno9")
               has_error = .true.
            else if (ocn_recon_code < 0) then
               call logger%error("&ocean_vmix_nml tracer_recon='"// &
                                 trim(cfg%ocean%vmix%tracer_recon)// &
                                 "' not recognised; valid: ppm, weno5, weno7, weno9")
               has_error = .true.
            else if (ocn_recon_code /= TRACER_RECON_PPM) then
               if (cfg%nghost < tracer_recon_required_nghost(ocn_recon_code)) then
                  call logger%error("&ocean_vmix_nml tracer_recon='"// &
                                    trim(cfg%ocean%vmix%tracer_recon)// &
                                    "' requires nghost >= "// &
                                    to_string(tracer_recon_required_nghost(ocn_recon_code))// &
                                    " but nghost = "//to_string(cfg%nghost))
                  has_error = .true.
               end if
            end if
         end block
      end if

      ! (PR-21) Shortwave source + boundary-layer coupling selectors.
      ! Fail-loud on unknown strings (PR-6 idiom); q_sw source requires
      ! the PR-12 component set (makes the host-side source branch total);
      ! the net-heat source with penetrating SW warns about the
      ! night-time negative-I0 hazard (§2.4 of the plan).
      if (trim(cfg%sim_type) == "ocean") then
         if (.not. sw_source_is_implemented(trim(cfg%ocean%thermo%sw_source))) then
            call logger%error("&ocean_thermo_nml sw_source='"// &
                              trim(cfg%ocean%thermo%sw_source)// &
                              "' not recognised; valid: net_heat, q_sw")
            has_error = .true.
         end if
         if (.not. kpp_sw_method_is_implemented(trim(cfg%ocean%thermo%kpp_sw_method))) then
            call logger%error("&ocean_thermo_nml kpp_sw_method='"// &
                              trim(cfg%ocean%thermo%kpp_sw_method)// &
                              "' not recognised; valid: all_sw, mxl_sw, lv1_sw")
            has_error = .true.
         end if
         if (trim(cfg%ocean%thermo%sw_source) == "q_sw" .and. &
             .not. cfg%ocean%forcing%enable_components) then
            call logger%error("&ocean_thermo_nml sw_source='q_sw' requires "// &
                              "&ocean_forcing_nml enable_components=.true. "// &
                              "(q_sw is allocated only with the component set)")
            has_error = .true.
         end if
         if (trim(cfg%ocean%thermo%sw_source) == "net_heat" .and. &
             cfg%ocean%thermo%sw_pen_frac > 0.0_wp) then
            call logger%warning("shortwave penetration is reading the NET heat flux "// &
                                "(&ocean_thermo_nml sw_source='net_heat', sw_pen_frac>0): "// &
                                "I0 < 0 wherever Q_heat < 0, which drives unphysical "// &
                                "negative irradiance down the two-band profile. Set "// &
                                "sw_source='q_sw' with &ocean_forcing_nml "// &
                                "enable_components=.true.")
         end if
      end if

      ! Sponge width
      if (cfg%sponge_width < 0) then
         call logger%error("Invalid sponge_width = "//to_string(cfg%sponge_width)// &
                           ": must be >= 0")
         has_error = .true.
      end if

      ! Warn about unknown BC strings (they silently default to WALL)
      call warn_unknown_bc(cfg%bc_west, "bc_west")
      call warn_unknown_bc(cfg%bc_east, "bc_east")
      call warn_unknown_bc(cfg%bc_south, "bc_south")
      call warn_unknown_bc(cfg%bc_north, "bc_north")

      ! `thickness_config` is consumed only by the ocean state builder's IC
      ! seed.  The coastal regimes lay their layers on a separate path and
      ! would silently ignore a non-default value, so fail loud instead.
      if (trim(cfg%thickness_config) /= "sigma" .and. &
          trim(cfg%sim_type) /= "ocean") then
         call logger%error("thickness_config = '"//trim(cfg%thickness_config)// &
                           "' is an ocean-path knob (sim_type = 'ocean'); the coastal "// &
                           "regimes seed their layers on a separate path and would "// &
                           "silently ignore it")
         has_error = .true.
      end if

      ! Retired coastal-legacy `&tracer_nml` linear-EOS quartet.  These
      ! four keys reach `tracer_t%eos_coeff`/`eos_ref` and nothing else;
      ! the C-grid ocean path's linear EOS reads `&ocean_ic_nml` alone.
      ! A knob that validates and silently does nothing is the bug, so
      ! moving any of them off its historical default is fatal and names
      ! the live replacement rather than being quietly ignored.
      if (cfg%alpha_T /= LEGACY_TRACER_ALPHA_T) then
         call logger%error("&tracer_nml alpha_T is RETIRED on the ocean path: it only ever "// &
                           "reached tracer_t%eos_coeff, which no ocean kernel reads, so "// &
                           "setting it changed nothing. Use &ocean_ic_nml alpha_T "// &
                           "(kg/m^3 per degC) instead, or delete the key.")
         has_error = .true.
      end if
      if (cfg%beta_S /= LEGACY_TRACER_BETA_S) then
         call logger%error("&tracer_nml beta_S is RETIRED on the ocean path: it only ever "// &
                           "reached tracer_t%eos_coeff, which no ocean kernel reads, so "// &
                           "setting it changed nothing. Use &ocean_ic_nml beta_S "// &
                           "(kg/m^3 per PSU) instead, or delete the key.")
         has_error = .true.
      end if
      if (cfg%T_ref /= LEGACY_TRACER_T_REF) then
         call logger%error("&tracer_nml T_ref is RETIRED on the ocean path: it only ever "// &
                           "reached tracer_t%eos_ref, which no ocean kernel reads, so "// &
                           "setting it changed nothing. Use &ocean_ic_nml T_ref (degC) "// &
                           "instead, or delete the key.")
         has_error = .true.
      end if
      if (cfg%S_ref /= LEGACY_TRACER_S_REF) then
         call logger%error("&tracer_nml S_ref is RETIRED on the ocean path: it only ever "// &
                           "reached tracer_t%eos_ref, which no ocean kernel reads, so "// &
                           "setting it changed nothing. Use &ocean_ic_nml S_ref (PSU) "// &
                           "instead, or delete the key.")
         has_error = .true.
      end if

      ! Linear-EOS stratified salinity IC: both ends must be set (the
      ! seed gates on `/= 0` for BOTH, mirroring T_init_surface/bottom),
      ! so exactly one non-zero is a silently-uniform column — the same
      ! class of bug the retired knobs above were.
      if ((cfg%S_init_surface /= 0.0_wp) .neqv. (cfg%S_init_bottom /= 0.0_wp)) then
         call logger%error("&tracer_nml S_init_surface / S_init_bottom must BOTH be "// &
                           "non-zero to build the linear S(z) IC (they gate together, "// &
                           "exactly as T_init_surface/T_init_bottom do); one alone is "// &
                           "ignored and the column stays uniform at initial_salinity.")
         has_error = .true.
      end if

      ! Freezing-point (liquidus) coefficient set.  Belt-and-braces on top
      ! of the `nml_enum allowed=` list, in the style of the
      ! `conc_config` check above — and NOT optional: `parse_tfreeze_set`
      ! deliberately has no default fallback, so a string that reaches
      ! `eos_apply_tfreeze_set` unrecognised would silently leave the
      ! sea-ice set in place.  At S = 34.5 the two shipped sets differ by
      ! ~0.03 degC, which is enough to flip the sign of an ice-shelf
      ! basal melt rate — a mistyped liquidus must stop the run.
      if (parse_tfreeze_set(cfg%ocean%eos%tfreeze_set) == TFREEZE_SET_INVALID) then
         call logger%error("&ocean_eos_nml tfreeze_set = '"// &
                           trim(cfg%ocean%eos%tfreeze_set)// &
                           "' is not recognised (allowed: seaice, isomip). "// &
                           "'seaice' is the SIS2/MOM6 sea-ice liquidus "// &
                           "(-0.054*S - 7.53e-8*p); 'isomip' is the ISOMIP+ "// &
                           "ice-shelf-cavity set (-0.0573*S + 0.0832 - 7.53e-8*p, "// &
                           "Asay-Davis et al. 2016 Table 4).")
         has_error = .true.
      end if

      ! Ocean diag manager
      if (trim(cfg%sim_type) == "ocean") then
         if (trim(cfg%ocean%diag%vgrid) /= "layer" .and. &
             trim(cfg%ocean%diag%vgrid) /= "z_fixed" .and. &
             trim(cfg%ocean%diag%vgrid) /= "sigma" .and. &
             trim(cfg%ocean%diag%vgrid) /= "zstar" .and. &
             trim(cfg%ocean%diag%vgrid) /= "density") then
            call logger%error("Invalid ocean_diag_vgrid = '"// &
                              trim(cfg%ocean%diag%vgrid)// &
                              "': must be 'layer', 'z_fixed', 'sigma', 'zstar', or 'density'")
            has_error = .true.
         end if
         if (trim(cfg%ocean%diag%vgrid) == "z_fixed" .and. &
             cfg%ocean%diag%n_z_levels <= 0) then
            call logger%error("ocean_diag_vgrid = 'z_fixed' requires "// &
                              "ocean_diag_n_z_levels > 0")
            has_error = .true.
         end if
         if (cfg%ocean%diag%n_z_levels < 0 .or. &
             cfg%ocean%diag%n_z_levels > MAX_OCEAN_DIAG_Z_LEVELS) then
            call logger%error("ocean_diag_n_z_levels = "// &
                              to_string(cfg%ocean%diag%n_z_levels)// &
                              " out of range [0, "// &
                              to_string(MAX_OCEAN_DIAG_Z_LEVELS)//"]")
            has_error = .true.
         end if
         ! Density bins have no auto-fill (unlike sigma/z*), so the global
         ! vgrid='density' selection OR any per-diagnostic ':density'/':rho'
         ! attribute in `diags` requires a valid rho_levels list.
         if (.not. diag_density_levels_ok(trim(cfg%ocean%diag%vgrid), &
                                          trim(cfg%ocean%diag%diags), &
                                          cfg%ocean%diag%n_rho_levels, &
                                          cfg%ocean%diag%rho_levels)) then
            call logger%error("&ocean_diag_nml vgrid='density' (or a ':density'/':rho' "// &
                              "entry in diags) requires n_rho_levels > 0, "// &
                              "rho_levels(1:n_rho_levels) strictly increasing, and "// &
                              "n_rho_levels <= "//to_string(MAX_OCEAN_DIAG_Z_LEVELS))
            has_error = .true.
         end if
         if (cfg%ocean%diag%enabled .and. cfg%ocean%diag%dt_out <= 0.0_wp) then
            call logger%error("ocean_diag_dt_out must be > 0 when "// &
                              "ocean_diag_enabled = .true.")
            has_error = .true.
         end if
         if (cfg%ocean%diag%enabled .and. cfg%ocean%diag%dt_out > cfg%t_end) then
            call logger%warning("ocean_diag_dt_out ("// &
                                to_string(cfg%ocean%diag%dt_out)//" s) > t_end ("// &
                                to_string(cfg%t_end)//" s): the diagnostic will "// &
                                "fire at most once — is dt_out in your time_unit "// &
                                "(e.g. days), not seconds?")
         end if
         if (trim(cfg%ocean%topo%topo_config) /= "flat" .and. &
             trim(cfg%ocean%topo%topo_config) /= "spoon" .and. &
             trim(cfg%ocean%topo%topo_config) /= "seamount" .and. &
             trim(cfg%ocean%topo%topo_config) /= "neverworld2" .and. &
             trim(cfg%ocean%topo%topo_config) /= "island" .and. &
             trim(cfg%ocean%topo%topo_config) /= "double_drake" .and. &
             trim(cfg%ocean%topo%topo_config) /= "isomip_plus" .and. &
             trim(cfg%ocean%topo%topo_config) /= "file") then
            call logger%error("Invalid topo_config = '"//trim(cfg%ocean%topo%topo_config)// &
                              "': must be 'flat', 'spoon', 'seamount', 'neverworld2', "// &
                              "'island', 'double_drake', 'isomip_plus', or 'file'")
            has_error = .true.
         end if
         if (trim(cfg%ocean%ic%ic_config) /= "" .and. &
             trim(cfg%ocean%ic%ic_config) /= "eady" .and. &
             trim(cfg%ocean%ic%ic_config) /= "geostrophic_adjustment" .and. &
             trim(cfg%ocean%ic%ic_config) /= "baroclinic_jet") then
            call logger%error("Invalid ic_config = '"//trim(cfg%ocean%ic%ic_config)// &
                              "': must be '', 'eady', 'geostrophic_adjustment', or "// &
                              "'baroclinic_jet'")
            has_error = .true.
         end if
         if (trim(cfg%ocean%topo%topo_config) == "file" .and. &
             len_trim(cfg%bathymetry_file) == 0) then
            call logger%error("topo_config = 'file' requires bathymetry_file "// &
                              "to be set in &output_nml")
            has_error = .true.
         end if
         if (trim(cfg%ocean%topo%wind_config) /= "constant" .and. &
             trim(cfg%ocean%topo%wind_config) /= "2gyre" .and. &
             trim(cfg%ocean%topo%wind_config) /= "neverworld2") then
            call logger%error("Invalid wind_config = '"//trim(cfg%ocean%topo%wind_config)// &
                              "': must be 'constant', '2gyre', or 'neverworld2'")
            has_error = .true.
         end if
         if (cfg%ocean%topo%max_depth <= 0.0_wp) then
            call logger%error("ocean_max_depth must be > 0")
            has_error = .true.
         end if
         if (trim(cfg%ocean%topo%topo_config) == "spoon") then
            if (cfg%ocean%topo%edge_depth <= 0.0_wp) then
               call logger%error("ocean_edge_depth must be > 0 for spoon bathymetry")
               has_error = .true.
            end if
            if (cfg%ocean%topo%edge_depth >= cfg%ocean%topo%max_depth) then
               call logger%error("ocean_edge_depth must be < ocean_max_depth")
               has_error = .true.
            end if
            if (cfg%ocean%topo%slope_scale <= 0.0_wp) then
               call logger%error("ocean_slope_scale must be > 0 for spoon bathymetry")
               has_error = .true.
            end if
         end if

         ! `uniform_z` seeds collapsed layers at the isopycnal `angstrom_h`
         ! floor, whereas the wet/dry path pins the emerged-column invariant
         ! `Sum h_layer = nz*2*H_VANISHED` and floors `bt_h` to exactly that
         ! same value.  Mixing the two would seed a bt_h that disagrees with
         ! Sum h_layer on every emerged column, so refuse the combination.
         if (trim(cfg%thickness_config) == "uniform_z" .and. &
             cfg%ocean%wetdry%enable) then
            call logger%error("thickness_config = 'uniform_z' is incompatible with "// &
                              "&ocean_wetdry_nml enable = .true.: the emerged-column "// &
                              "seed invariant (bt_h = nz*2*H_VANISHED) assumes the "// &
                              "'sigma' even split")
            has_error = .true.
         end if

         ! Stretched z* nominal profile.  Default "uniform" ⇒ no check
         ! fires and nothing downstream changes.  Anything else must be on
         ! a family that reads it — `z_fixed` (its levels), `zstar` (MOM6 z*
         ! levels) or `hycom` (its z* nominal floor); silently ignoring it
         ! would be the bug — and must build: list length = nz_layers, tanh
         ! parameters in range and leaving room to stretch.
         block
            integer :: zf_code, zf_ierr
            real(wp), allocatable :: zf_dz(:)
            zf_code = parse_z_fixed_profile(cfg%z_fixed_profile)
            if (zf_code == ZFIXED_PROFILE_INVALID) then
               call logger%error("&vcoord_nml z_fixed_profile = '"// &
                                 trim(cfg%z_fixed_profile)//"' is not one of "// &
                                 "'uniform', 'list', 'tanh'")
               has_error = .true.
            else if (zf_code /= ZFIXED_PROFILE_UNIFORM) then
               if (parse_vcoord_type(cfg%vcoord_type, VCOORD_EULERIAN_Z) /= VCOORD_Z_FIXED .and. &
                   parse_vcoord_type(cfg%vcoord_type, VCOORD_EULERIAN_Z) /= VCOORD_ZSTAR .and. &
                   parse_vcoord_type(cfg%vcoord_type, VCOORD_EULERIAN_Z) /= VCOORD_HYCOM) then
                  call logger%error("&vcoord_nml z_fixed_profile = '"// &
                                    trim(cfg%z_fixed_profile)//"' is only read by "// &
                                    "vcoord_type = 'z_fixed', 'zstar' or 'hycom' (got '"// &
                                    trim(cfg%vcoord_type)//"'); it would be silently ignored")
                  has_error = .true.
               else if (cfg%nz_layers >= 1) then
                  allocate (zf_dz(cfg%nz_layers))
                  call z_fixed_nominal_dz(zf_code, cfg%nz_layers, cfg%ocean%topo%max_depth, &
                                          cfg%z_fixed_dz, cfg%z_fixed_dz_top, &
                                          cfg%z_fixed_tanh_center, cfg%z_fixed_tanh_width, &
                                          zf_dz, zf_ierr)
                  if (zf_ierr == ZFIXED_DZ_ERR_COUNT) then
                     call logger%error("&vcoord_nml z_fixed_profile = 'list' needs exactly "// &
                                       "nz_layers = "//to_string(cfg%nz_layers)// &
                                       " leading positive z_fixed_dz entries (surface first, "// &
                                       "no gaps), got "//to_string(count(cfg%z_fixed_dz > 0.0_wp)))
                     has_error = .true.
                  else if (zf_ierr == ZFIXED_DZ_ERR_TOO_DEEP) then
                     call logger%error("&vcoord_nml z_fixed_profile = 'tanh': nz_layers * "// &
                                       "z_fixed_dz_top = "// &
                                       to_string(real(cfg%nz_layers, wp)*cfg%z_fixed_dz_top)// &
                                       " m is not below &ocean_topo_nml max_depth = "// &
                                       to_string(cfg%ocean%topo%max_depth)// &
                                       " m — there is no depth left to stretch into")
                     has_error = .true.
                  else if (zf_ierr /= ZFIXED_DZ_OK) then
                     call logger%error("&vcoord_nml z_fixed_profile = 'tanh' needs "// &
                                       "&ocean_topo_nml max_depth > 0, z_fixed_dz_top > 0, "// &
                                       "z_fixed_tanh_width > 0 and 0 <= z_fixed_tanh_center <= 1")
                     has_error = .true.
                  else if (zf_code == ZFIXED_PROFILE_LIST .and. &
                           sum(zf_dz) < cfg%ocean%topo%max_depth) then
                     call logger%warning("&vcoord_nml z_fixed_dz sums to "// &
                                         to_string(sum(zf_dz))//" m, shallower than "// &
                                         "&ocean_topo_nml max_depth = "// &
                                         to_string(cfg%ocean%topo%max_depth)// &
                                         " m: columns deeper than the profile carry the "// &
                                         "excess in their bed layer")
                  end if
               end if
            end if
         end block

         ! Target-density profile of the density families.  Default
         ! "uniform" ⇒ the light→dense linspace, nothing checked.  "list"
         ! must be on rho/hycom (nothing else reads it) and give exactly
         ! nz_layers+1 strictly increasing interface densities, light first.
         block
            integer :: n_rho, kr, vc_code
            logical :: rho_ok
            vc_code = parse_vcoord_type(cfg%vcoord_type, VCOORD_EULERIAN_Z)
            if (trim(cfg%rho_target_profile) /= "uniform" .and. &
                trim(cfg%rho_target_profile) /= "list") then
               call logger%error("&vcoord_nml rho_target_profile = '"// &
                                 trim(cfg%rho_target_profile)//"' is not one of "// &
                                 "'uniform', 'list'")
               has_error = .true.
            else if (trim(cfg%rho_target_profile) == "list") then
               if (vc_code /= VCOORD_RHO .and. vc_code /= VCOORD_HYCOM) then
                  call logger%error("&vcoord_nml rho_target_profile = 'list' is only "// &
                                    "read by vcoord_type = 'rho' or 'hycom' (got '"// &
                                    trim(cfg%vcoord_type)//"'); it would be silently ignored")
                  has_error = .true.
               else
                  ! Leading positive entries, no gaps.
                  n_rho = 0
                  do kr = 1, size(cfg%rho_target_list)
                     if (cfg%rho_target_list(kr) <= 0.0_wp) exit
                     n_rho = n_rho + 1
                  end do
                  rho_ok = n_rho == cfg%nz_layers + 1 .and. &
                           count(cfg%rho_target_list > 0.0_wp) == n_rho
                  if (rho_ok) then
                     do kr = 2, n_rho
                        if (cfg%rho_target_list(kr) <= cfg%rho_target_list(kr - 1)) then
                           rho_ok = .false.
                        end if
                     end do
                  end if
                  if (.not. rho_ok) then
                     call logger%error("&vcoord_nml rho_target_profile = 'list' needs "// &
                                       "exactly nz_layers+1 = "// &
                                       to_string(cfg%nz_layers + 1)// &
                                       " leading positive, strictly increasing "// &
                                       "rho_target_list entries (lightest first, no gaps), "// &
                                       "got "//to_string(n_rho))
                     has_error = .true.
                  end if
               end if
            else if (count(cfg%rho_target_list > 0.0_wp) > 0) then
               call logger%error("&vcoord_nml rho_target_list is set but "// &
                                 "rho_target_profile = 'uniform' ignores it; set "// &
                                 "rho_target_profile = 'list'")
               has_error = .true.
            end if
         end block

         ! `VCOORD_ZSIGMA` is NOT a working coordinate on the ocean path.
         ! Its deep branch reads `z_ref_global` as a table of absolute
         ! reference depths in METRES
         ! (`z_top_k = min(z_ref_global(nz-k), column_total)`,
         ! `rdb_ocean_vcoord :: ocean_vcoord_compute_target_h_impl`), but the
         ! ONLY writer of that array anywhere in `src/` is the DIMENSIONLESS
         ! `z_ref_global(k) = k/nz` init in `ocean_vcoord_init` — nothing on
         ! the namelist path or the Python path ever replaces it with metres.
         ! So every z-level interval is `1/nz` metres, every layer collapses,
         ! and the whole column is dumped into `target_h(:,:,1)` (the BED
         ! layer) by the deficit line.  `Sum target_h = H + eta` still holds,
         ! which is exactly why no conservation test ever caught it; the
         ! placement is measured interface-by-interface in
         ! `test_ocean_vcoord_interface_depths ::
         ! documents_zsigma_dimensionless_zref_collapse` (nine 0.1 m layers
         ! in the top 90 cm of a 1000 m column).
         !
         ! Refuse it rather than silently running a broken coordinate.  No
         ! shipped namelist selects it.  Follow-up: fill `z_ref_global` in
         ! metres (the family is the natural seat for a sigma-near-the-top /
         ! z-below HYBRID), then delete this refusal and the `documents_*`
         ! test with it.  `VCOORD_ZSTAR_SIGMA` consumes the same table
         ! FRACTIONALLY (rescaled by `z_ref_global(nz)`) so it is unaffected
         ! by the units — but note that with a uniform table its deep branch
         ! is numerically indistinguishable from SIGMA.
         if (parse_vcoord_type(cfg%vcoord_type, default_code=VCOORD_EULERIAN_Z) &
             == VCOORD_ZSIGMA) then
            call logger%error("&vcoord_nml vcoord_type = 'zsigma' is refused on the "// &
                              "ocean path: its deep branch reads `z_ref_global` as "// &
                              "absolute depths in METRES, but the only writer of that "// &
                              "table is the dimensionless `k/nz` init, so every "// &
                              "z-level interval is 1/nz metres and the whole column "// &
                              "collapses into the bed layer (the column sum is still "// &
                              "exact, which is why it looked healthy).  Use 'sigma', "// &
                              "'zstar', 'zstar_sigma' or 'zstar_full'; ZSIGMA returns "// &
                              "when `z_ref_global` is filled in metres.")
            has_error = .true.
         end if

         ! `&vcoord_nml zstar_h_min` carries TWO different contracts, picked
         ! by the coordinate family rather than by the value (see
         ! `rdb_vcoord :: vcoord_h_min_role`).  On the GEOMETRIC families
         ! (ZSTAR_FULL / Z_FIXED) it is the anti-zero thickness of below-bed
         ! FILLER layers that are meant to read as vanished downstream, so it
         ! belongs at or below the D4 skip/merge marker `H_VANISHED`; on the
         ! DENSITY families (RHO / HYCOM) the collapsed layers carry real
         ! tracer mass and that path deliberately floors at
         ! `max(zstar_h_min, 2*H_VANISHED)` instead.  Nothing else pins the
         ! knob, so an INERT-role run with `zstar_h_min > H_VANISHED` silently
         ! promotes its below-bed filler to LIVE layers (real EOS density off
         ! ghost T/S, a PGF column entry, a remap-drain concentration, a vdiff
         ! interface) while the coordinate still treats them as throwaway.
         !
         ! This was a WARNING until the rigid-top work made it blocking.  The
         ! reason it is now an ERROR: a coordinate that vanishes layers
         ! against the TOP of the column (an ice-shelf cavity) puts its
         ! fillers where the surface fluxes, the pressure gradient, the melt
         ! sampler and the tracer budgets all read — so a filler that is not
         ! skipped is not a cosmetic slip, it is a conservation hole.  The
         ! only configuration in the tree that sat in the band was the repo's
         ! own Python worked example (`python/tests/test_worked_example.py`,
         ! `ZStarFull(h_min=1.0e-3)`), corrected in the same change.
         !
         ! NOTE the boundary this must NOT move: five shipped namelists set
         ! `zstar_h_min = 1.5e-4`, H_VANISHED EXACTLY, which is legal —
         ! `vcoord_h_min_is_coherent` is a strict `>` and every downstream
         ! vanish test is a strict `> H_VANISHED`, so a layer on the marker
         ! reads as vanished.  Relaxing either to `>=` would refuse the
         ! canonical double-gyre reference; `test_ocean_vcoord_hygiene ::
         ! h_min_on_the_marker_is_accepted` guards that.
         !
         ! A non-positive floor was already refused and still is — it defeats
         ! the knob's single documented purpose.
         block
            integer :: hmin_vcoord_code
            hmin_vcoord_code = parse_vcoord_type(cfg%vcoord_type, &
                                                 default_code=VCOORD_EULERIAN_Z)
            if (.not. vcoord_h_min_is_coherent(hmin_vcoord_code, cfg%zstar_h_min)) then
               if (cfg%zstar_h_min <= 0.0_wp) then
                  call logger%error("&vcoord_nml zstar_h_min = "// &
                                    to_string(cfg%zstar_h_min)//" must be > 0: it exists "// &
                                    "so a vanishing layer's target thickness is never "// &
                                    "exactly zero (kernels that divide by h_layer)")
               else
                  call logger%error("&vcoord_nml zstar_h_min = "// &
                                    to_string(cfg%zstar_h_min)//" m exceeds H_VANISHED = "// &
                                    to_string(H_VANISHED)//" m under vcoord_type = '"// &
                                    trim(cfg%vcoord_type)//"': that family uses the knob "// &
                                    "as an anti-zero floor for filler layers that are "// &
                                    "MEANT to stay vanished — above H_VANISHED they "// &
                                    "become dynamically live (EOS/PGF/remap-drain/vdiff) "// &
                                    "while the coordinate still treats them as "// &
                                    "throwaway.  Use a value <= "// &
                                    to_string(H_VANISHED)//"; for a genuinely live "// &
                                    "minimum layer thickness use &ocean_isopycnal_nml "// &
                                    "angstrom_h (the D4 floor knob); the rho/hycom "// &
                                    "regrid has its own keep-alive floor")
               end if
               has_error = .true.
            end if
            ! NOTE the boundary, deliberately NOT warned about at runtime:
            ! the five shipped namelists set `zstar_h_min = 1.5e-4`, which is
            ! H_VANISHED EXACTLY.  That is legal — every downstream vanish
            ! test is a strict `> H_VANISHED`, so a layer sitting on the
            ! marker still reads as vanished — but it carries zero margin,
            ! and any gate relaxed to `>=` would change those runs' answers.
            ! A per-run warning here would fire on the canonical double-gyre
            ! reference forever while recommending nothing an operator can do
            ! without changing answers, so the fact lives in the docs and at
            ! the gate (`rdb_ocean_remap :: H_FLOOR`) instead.
         end block

         ! GM thickness diffusion consumes the stored isopycnal slope, so
         ! the slopes slot MUST be enabled.  Loud invariant (the kernel is
         ! pure/device and cannot fail loud); a silent no-op would hide a
         ! misconfigured GM run.
         if (cfg%ocean%gm%enable .and. .not. cfg%ocean%slopes%enable) then
            call logger%error("&ocean_gm_nml enable=.true. requires "// &
                              "&ocean_slopes_nml enable=.true. (GM reads the stored slope)")
            has_error = .true.
         end if
         ! khth_slope_max enters as 1/slope_max^2 in the safe-streamfunction
         ! limiter; a zero/negative value divides by zero.
         if (cfg%ocean%gm%enable .and. cfg%ocean%gm%khth_slope_max <= 0.0_wp) then
            call logger%error("&ocean_gm_nml khth_slope_max must be > 0 "// &
                              "(enters the safe-streamfunction limiter as 1/slope_max^2)")
            has_error = .true.
         end if

         ! Redi (capability [3]) ships the CONTINUOUS variant only; the
         ! discontinuous (regula-falsi) path is deferred (R4).  Loud
         ! invariant — the device kernel cannot fail loud.
         if (cfg%ocean%redi%enable .and. .not. cfg%ocean%redi%continuous) then
            call logger%error("&ocean_redi_nml continuous=.false. (discontinuous "// &
                              "variant) is not implemented yet (deferred R4)")
            has_error = .true.
         end if

         ! VarMix (capability [4]) needs the stored slope + N² (slopes slot)
         ! for the Eady term AND the first-mode cg1 (wavespeed slot) for the
         ! resolution function.  Loud invariant — the kernel is pure/device.
         if (cfg%ocean%varmix%enable .and. .not. cfg%ocean%slopes%enable) then
            call logger%error("&ocean_varmix_nml enable=.true. requires "// &
                              "&ocean_slopes_nml enable=.true. (VarMix reads the stored slope + N^2)")
            has_error = .true.
         end if
         if (cfg%ocean%varmix%enable .and. .not. cfg%ocean%wavespeed%enable) then
            call logger%error("&ocean_varmix_nml enable=.true. requires "// &
                              "&ocean_wavespeed_nml enable=.true. (the resolution function needs cg1)")
            has_error = .true.
         end if
         ! The resolution-function power is an integer >= 1 (Hallberg 2013
         ! form `s^p/(s^p + (alpha*cg1)^p)`); p < 1 has no meaning and would
         ! divide by zero or invert the limits.  Refused regardless of enable.
         if (cfg%ocean%varmix%kh_res_fn_power < 1) then
            call logger%error("&ocean_varmix_nml kh_res_fn_power must be >= 1 "// &
                              "(the resolution-function power p; 2 is the production form)")
            has_error = .true.
         end if

         ! Resolution-scaled momentum viscosity (Gap 1) reads the VarMix
         ! resolution function `Res_fn`, so VarMix MUST be enabled.  Loud
         ! invariant — the lateral-mix kernel is pure/device.
         if (cfg%ocean%hvisc%resoln_scaled_visc .and. .not. cfg%ocean%varmix%enable) then
            call logger%error("&ocean_hvisc_nml resoln_scaled_visc=.true. requires "// &
                              "&ocean_varmix_nml enable=.true. (the resolution function lives in VarMix)")
            has_error = .true.
         end if

         ! MEKE (capability [5]) sources its eddy energy from the GM PE
         ! release (`gm%gm_src`), so GM MUST be enabled.  Loud invariant —
         ! the kernel is pure/device and cannot fail loud.  (The VarMix
         ! feedback seam is OPTIONAL: with VarMix off, MEKE still evolves E
         ! but `meke%kh` has no face accumulator to feed — documented, not an
         ! error.)
         if (cfg%ocean%meke%enable .and. .not. cfg%ocean%gm%enable) then
            call logger%error("&ocean_meke_nml enable=.true. requires "// &
                              "&ocean_gm_nml enable=.true. (MEKE sources from gm_src)")
            has_error = .true.
         end if

         ! ---- Lateral-mixing closure: fail loud on unimplemented tags ----
         ! `parse_lateral_closure` maps the &ocean_hvisc_nml string to an
         ! LMIX_* code; a mistyped/garbage string yields LMIX_INVALID and a
         ! tag with no dispatcher path is not implemented.  Either case must
         ! ABORT here — never silently fall through to background-only
         ! viscosity (the original silent-wrong-answer foot-gun).
         lateral_closure_code = parse_lateral_closure(cfg%ocean%hvisc%lateral_closure)
         if (.not. lateral_closure_is_implemented(lateral_closure_code)) then
            call logger%error("Invalid/unimplemented &ocean_hvisc_nml lateral_closure = '"// &
                              trim(cfg%ocean%hvisc%lateral_closure)//"': must be one of "// &
                              "'none', 'leith', 'smagorinsky'/'smag', 'biharmonic', "// &
                              "'leith_biharm' — a closure with no kernel must not "// &
                              "silently disable lateral viscosity")
            has_error = .true.
         end if
         ! `leith_biharm` and `smag_ah` are BOTH flow-aware biharmonic
         ! closures that fill the same `nu4_face_*` arrays; the dispatcher
         ! runs the closure first, then `compute_smag_ah` unconditionally,
         ! so an enabled `smag_ah` would silently OVERWRITE the
         ! Leith-biharmonic fill.  Fail loud rather than let one closure
         ! silently win (MOM6 max-combines them; we do not).
         if (lateral_closure_conflicts_smag_ah(lateral_closure_code, &
                                               cfg%ocean%hvisc%smag_ah)) then
            call logger%error("&ocean_hvisc_nml lateral_closure='leith_biharm' and "// &
                              "smag_ah=.true. are both flow-aware biharmonic closures "// &
                              "filling nu4_face — smag_ah would overwrite the "// &
                              "Leith-biharmonic fill; enable only one")
            has_error = .true.
         end if
         ! Only anisotropy mode 0 (grid-relative `aniso_dir`) has an
         ! implemented direction tensor; a requested-but-unimplemented mode
         ! must abort, not silently fall back to the grid-i default.
         if (.not. aniso_mode_is_implemented(cfg%ocean%hvisc%aniso_mode)) then
            call logger%error("&ocean_hvisc_nml aniso_mode must be 0 "// &
                              "(grid-relative aniso_dir) — other modes are not "// &
                              "implemented and must not silently fall back to grid-i")
            has_error = .true.
         end if
         ! MEKE backscatter injects a NEGATIVE harmonic viscosity that
         ! amplifies grid-scale modes; only a POSITIVE biharmonic
         ! (nu_4 / smag_ah / leith_biharm) can dissipate them.  Abort if
         ! backscatter is on without a biharmonic backstop configured —
         ! the CFL floor bounds the growth rate, it does not make a
         ! negative harmonic operator stable on its own.
         if (cfg%ocean%meke%backscatter .and. &
             .not. has_biharmonic_backstop(cfg%ocean%hvisc%nu_4, &
                                           cfg%ocean%hvisc%smag_ah, &
                                           cfg%ocean%hvisc%smag_bi_const, &
                                           lateral_closure_code, &
                                           cfg%ocean%hvisc%c_leith_bi, &
                                           cfg%ocean%hvisc%nu_4_bg)) then
            call logger%error("&ocean_meke_nml backscatter=.true. requires a "// &
                              "biharmonic backstop with a NON-ZERO coefficient "// &
                              "(&ocean_hvisc_nml nu_4>0 with no flow-aware "// &
                              "closure selected, or smag_ah=.true. with "// &
                              "smag_bi_const>0, or lateral_closure='leith_biharm' "// &
                              "with c_leith_bi>0 — or nu_4_bg>0 as a floor in "// &
                              "either flow-aware case) — the negative harmonic "// &
                              "backscatter is unstable without one")
            has_error = .true.
         end if

         ! ---- PR-6 fail-loud dispatch pack ----
         ! Each of these string→enum / config→kernel seams previously
         ! accepted a value the schema (or a hand-written reader) advertised,
         ! then silently ran DIFFERENT physics than the name promised.  Each
         ! guard consumes a `pure` predicate next to its dispatcher and
         ! aborts once via `has_error` — same idiom as the lateral-closure
         ! guard above.  Every guard is inert on valid config (bit-identical).
         ! Coriolis form: a typo (→ PV_VARIANT_INVALID) or the
         ! reserved-but-unwired `al81` (→ PV_VARIANT_AL81) must NOT silently
         ! run enstrophy-only Sadourny — a different conservation law.
         if (.not. pv_variant_is_implemented( &
             parse_pv_variant(cfg%ocean%coriolis%form))) then
            call logger%error("&ocean_coriolis_nml form = '"// &
                              trim(cfg%ocean%coriolis%form)//"' is not implemented — "// &
                              "must be 'sadourny', 'sadourny_hk' or 'sadourny_energy' "// &
                              "('al81' is reserved but its Arakawa-Lamb kernel is not "// &
                              "yet wired)")
            has_error = .true.
         end if
         ! PV face-interpolation scheme: a typo (→ PV_ADV_INVALID) must not
         ! silently fall back to centred.
         if (.not. pv_adv_scheme_is_implemented( &
             parse_pv_adv_scheme(cfg%ocean%coriolis%pv_adv_scheme))) then
            call logger%error("&ocean_coriolis_nml pv_adv_scheme = '"// &
                              trim(cfg%ocean%coriolis%pv_adv_scheme)//"' is not recognised — "// &
                              "must be 'centered', 'weno3', 'weno5' or 'weno7'")
            has_error = .true.
         end if
         ! weno5/weno7 (stencil radius 3/4) need a halo of radius + 1 = 4/5
         ! (at radius alone a decomposed run is not bit-identical to one rank;
         ! see `pv_adv_required_nghost`) — fail-loud, mirroring the
         ! tracer-WENO ladder's per-rung nghost gate.
         if (cfg%nghost < pv_adv_required_nghost( &
             parse_pv_adv_scheme(cfg%ocean%coriolis%pv_adv_scheme))) then
            call logger%error("&ocean_coriolis_nml pv_adv_scheme='"// &
                              trim(cfg%ocean%coriolis%pv_adv_scheme)//"' requires nghost >= "// &
                              to_string(pv_adv_required_nghost( &
                                        parse_pv_adv_scheme(cfg%ocean%coriolis%pv_adv_scheme)))// &
                              " but nghost = "//to_string(cfg%nghost))
            has_error = .true.
         end if
         ! WENO PV reconstruction is wired only into the Sadourny enstrophy
         ! path; pairing it with the hk/energy forms would silently run the
         ! centred interpolation those kernels hard-code.
         if (trim(cfg%ocean%coriolis%pv_adv_scheme) /= "centered" .and. &
             trim(cfg%ocean%coriolis%form) /= "sadourny") then
            call logger%error("&ocean_coriolis_nml pv_adv_scheme='"// &
                              trim(cfg%ocean%coriolis%pv_adv_scheme)// &
                              "' is only wired into form='sadourny'; got form='"// &
                              trim(cfg%ocean%coriolis%form)//"'")
            has_error = .true.
         end if
         ! Mass-consistent CorAdCalc (use_state_fluxes): only the
         ! sadourny_energy transport form consumes uh/vh, and only the
         ! mom6 corrector has a predictor continuity solve to be
         ! consistent WITH — any other combination would silently run
         ! the recompute path while the namelist promised otherwise.
         if (cfg%ocean%coriolis%use_state_fluxes) then
            if (trim(cfg%ocean%coriolis%form) /= "sadourny_energy") then
               call logger%error("&ocean_coriolis_nml use_state_fluxes requires "// &
                                 "form='sadourny_energy' (the transport form is what "// &
                                 "consumes uh/vh); got form='"// &
                                 trim(cfg%ocean%coriolis%form)//"'")
               has_error = .true.
            end if
            if (trim(cfg%ocean%bt%split_scheme) /= "pred_corr") then
               call logger%error("&ocean_coriolis_nml use_state_fluxes requires "// &
                                 "&ocean_bt_nml split_scheme='pred_corr' (the corrector "// &
                                 "consumes the predictor chain's fluxes); got "// &
                                 "split_scheme='"//trim(cfg%ocean%bt%split_scheme)//"'")
               has_error = .true.
            end if
         end if
         ! BOUND_CORIOLIS clamps the ENERGY-scheme CAu (q·vh) to the
         ! (f+ζ)·v velocity-form range; it is wired into the energy impl
         ! only (MOM6 also applies it to the enstrophy scheme, but Roundabout's
         ! enstrophy path is a separate kernel — v1 covers the audit target).
         ! Silently ignoring it on another form would be a positivity-request
         ! foot-gun, so fail loud.
         if (cfg%ocean%coriolis%bound_coriolis .and. &
             trim(cfg%ocean%coriolis%form) /= "sadourny_energy") then
            call logger%error("&ocean_coriolis_nml bound_coriolis is implemented for "// &
                              "form='sadourny_energy' only (the energy-scheme q·vh blow-up "// &
                              "it cures); got form='"//trim(cfg%ocean%coriolis%form)//"'")
            has_error = .true.
         end if
         ! corner_h="mom6_area" is wired into the energy impl's Pass 2 only
         ! (MOM6 shares the q construction across schemes, but Roundabout's
         ! enstrophy/hk are separate kernels — v1 covers the audit target).
         if (trim(cfg%ocean%coriolis%corner_h) == "mom6_area" .and. &
             trim(cfg%ocean%coriolis%form) /= "sadourny_energy") then
            call logger%error("&ocean_coriolis_nml corner_h='mom6_area' is implemented for "// &
                              "form='sadourny_energy' only; got form='"// &
                              trim(cfg%ocean%coriolis%form)//"'")
            has_error = .true.
         end if
         ! Bottom-drag form: linear (τ=ρ·r·u) and quadratic (τ=ρ·C_d·|u|·u)
         ! obey different physics; a typo must abort, not silently pick the
         ! quadratic default.
         if (.not. bdrag_variant_is_implemented( &
             parse_bdrag_variant(cfg%ocean%bdrag%form))) then
            call logger%error("&ocean_bdrag_nml form = '"// &
                              trim(cfg%ocean%bdrag%form)//"' is not implemented — "// &
                              "must be 'linear'/'rayleigh' or 'quadratic'/'cd'")
            has_error = .true.
         end if
         ! gprime PGF is a 2-layer reduced-gravity form — it hard-writes
         ! ONLY the top+bottom layer, so nz/=2 silently zeros the PGF on the
         ! other layers.  Config-time check against nz_layers (not nz_ml).
         if (.not. gprime_nz_is_supported( &
             parse_opgf_variant(cfg%ocean%pgf%form), cfg%nz_layers)) then
            call logger%error("&ocean_pgf_nml form = 'gprime' requires "// &
                              "&nonhydrostatic_nml nz_layers = 2 (the reduced-gravity "// &
                              "form writes only the top + bottom layer; use "// &
                              "the default form='mont' for a general-nz PGF)")
            has_error = .true.
         end if
         ! Leith-biharmonic with c_leith_bi<=0 is a provable no-op (ν₄ is
         ! linear in c_leith_bi) — the user asked for biharmonic dissipation
         ! and got none.  Promoted from a configure-time warning to an abort.
         if (leith_biharm_is_inert(lateral_closure_code, cfg%ocean%hvisc%c_leith_bi)) then
            call logger%error("&ocean_hvisc_nml lateral_closure='leith_biharm' with "// &
                              "c_leith_bi <= 0 is inert (ν₄ is linear in c_leith_bi) — "// &
                              "set c_leith_bi > 0 or choose a different closure")
            has_error = .true.
         end if
         ! Tidal mixing enabled with no energy source (e_uniform<=0 and
         ! e_compute off) makes Kd ≡ 0 — the whole flux sweep runs for
         ! nothing.
         if (tidal_mixing_is_inert(cfg%ocean%tidal_mixing%enable, &
                                   cfg%ocean%tidal_mixing%e_uniform, &
                                   cfg%ocean%tidal_mixing%e_compute)) then
            call logger%error("&ocean_tidal_mixing_nml enable=.true. with e_uniform <= 0 "// &
                              "and e_compute=.false. supplies no energy — Kd is "// &
                              "identically zero; set e_uniform > 0 or e_compute=.true.")
            has_error = .true.
         end if
         ! OBC edge strings bypass the schema (external group), so a typo
         ! silently closed the boundary to a wall.  Reject INVALID per edge,
         ! naming which edge (the guard runs before configure_ocean_bc, so an
         ! invalid tag never reaches the setup path).
         if (ocean_bc_type_from_string(cfg%ocean%bc%west) == OBC_INVALID) then
            call logger%error("&ocean_bc_nml west = '"//trim(cfg%ocean%bc%west)// &
                              "' is not a recognised boundary type (wall/open/tidal/"// &
                              "nested/inflow/discharge/clamped/sponge/chapman/periodic/"// &
                              "tripolar_fold)")
            has_error = .true.
         end if
         if (ocean_bc_type_from_string(cfg%ocean%bc%east) == OBC_INVALID) then
            call logger%error("&ocean_bc_nml east = '"//trim(cfg%ocean%bc%east)// &
                              "' is not a recognised boundary type")
            has_error = .true.
         end if
         if (ocean_bc_type_from_string(cfg%ocean%bc%south) == OBC_INVALID) then
            call logger%error("&ocean_bc_nml south = '"//trim(cfg%ocean%bc%south)// &
                              "' is not a recognised boundary type")
            has_error = .true.
         end if
         if (ocean_bc_type_from_string(cfg%ocean%bc%north) == OBC_INVALID) then
            call logger%error("&ocean_bc_nml north = '"//trim(cfg%ocean%bc%north)// &
                              "' is not a recognised boundary type")
            has_error = .true.
         end if

         ! Along-coordinate tracer Laplacian (`&ocean_hdiff_nml kappa_h`):
         ! the explicit forward-Euler stability bound
         ! kappa_h*dt_therm*(1/dx^2+1/dy^2) <= 0.5 used to be checked HERE
         ! with the nominal `cfg%dx`/`cfg%dy` — broken on non-Cartesian
         ! grids (DEGREES there, not metres) and blind to the true minimum
         ! cell even on Cartesian.  MOVED to `ocean_stability_audit`
         ! (`rdb_ocean_stability_audit.F90`), which runs AFTER the real
         ! per-cell `ocean_metrics_t` arrays are built and works on every
         ! grid type.  See that module's docstring, or
         ! docs/CAPABILITIES_AND_LIMITATIONS.md, for the full check list.
      end if

      ! ---- Implicit vdiff stress/drag folding (mutual exclusions) ----
      ! `implicit_drag` (fold into the vdiff bed diagonal) is mutually
      ! exclusive with the legacy `&ocean_bdrag_nml implicit` split-apply
      ! path (both would damp the bed velocity ⇒ double drag), and the
      ! bed-only fold cannot represent HBBL-distributed drag (`hbbl > 0`)
      ! — that needs a per-layer rate (follow-up PR).  Fail loud.  (Not a
      ! `k = 1` problem: the fold and the HBBL band both start at the
      ! face's first live layer `k_bot_u/v`; what is missing is a 3-D
      ! `lambda_bot` the vdiff interior rows can add to their diagonal.)
      if (cfg%ocean%vdiff%implicit_drag .and. cfg%ocean%bdrag%implicit) then
         call logger%error("ocean_vdiff implicit_drag is mutually exclusive with "// &
                           "ocean_bdrag implicit (split-apply): set only one")
         has_error = .true.
      end if
      if (cfg%ocean%vdiff%implicit_drag .and. cfg%ocean%bdrag%hbbl > 0.0_wp &
          .and. .not. bbl_glue_is_effective(cfg)) then
         ! Exempt only while the glue is EFFECTIVE (it carries the band):
         ! a requested glue that setup turns off (no drag coefficient)
         ! would otherwise leave the unsupported fold + band running.
         call logger%error("ocean_vdiff implicit_drag does not yet support "// &
                           "HBBL-distributed drag (ocean_bdrag hbbl > 0); use the "// &
                           "split-apply path (ocean_bdrag implicit) for HBBL, or "// &
                           "bbl_glue with a nonzero bottom drag (quadratic cd > 0, "// &
                           "or linear r > 0 with bg_vel > 0)")
         has_error = .true.
      end if
      ! `implicit_top_drag` is the `k = nz` twin of the rule above, plus
      ! one of its own: the surface row is the row the WIND stress owns
      ! as a Neumann RHS, so the fold both adds a diagonal term and
      ! masks that RHS off on the faces the ice covers.  That is only
      ! meaningful with a top drag configured.
      if (cfg%ocean%vdiff%implicit_top_drag) then
         if (.not. cfg%ocean%tdrag%enable) then
            call logger%error("&ocean_vdiff_nml implicit_top_drag=.true. requires "// &
                              "&ocean_tdrag_nml enable=.true.  The fold consumes the "// &
                              "top-drag slot's lambda_top_u/v and its face cover "// &
                              "masks; with the slot disabled those are placeholder "// &
                              "arrays and the knob would silently do nothing except "// &
                              "look like a top drag was configured.")
            has_error = .true.
         end if
         if (cfg%ocean%tdrag%implicit) then
            call logger%error("&ocean_vdiff_nml implicit_top_drag is mutually "// &
                              "exclusive with &ocean_tdrag_nml implicit: both damp "// &
                              "the top layer, so running both is a DOUBLE COUNT, not "// &
                              "a stronger drag.  Pick one — the vdiff fold if the "// &
                              "column also has real vertical viscosity to couple "// &
                              "against, the in-kernel backward-Euler form otherwise.")
            has_error = .true.
         end if
         if (cfg%ocean%tdrag%htbl > 0.0_wp) then
            call logger%error("&ocean_vdiff_nml implicit_top_drag does not support "// &
                              "the HTBL-distributed top drag (&ocean_tdrag_nml htbl "// &
                              "> 0): the fold is a SINGLE k = nz Rayleigh rate on the "// &
                              "diagonal and cannot represent a band spread over "// &
                              "several layers (the mirror of the implicit_drag/HBBL "// &
                              "restriction).  Use &ocean_tdrag_nml implicit for a "// &
                              "distributed top drag.")
            has_error = .true.
         end if
      end if
      ! `bbl_glue` (MOM6 BOTTOMDRAGLAW) needs the height-above-bed stack
      ! that only the hvel_mom6 path accumulates.  It composes with every
      ! drag form and fold: its piston IS the bed sink (the explicit apply,
      ! the `&ocean_bdrag_nml implicit` split-apply and the `implicit_drag`
      ! fold are all skipped on the layers; the explicit tendency still
      ! feeds the barotropic F_slow, as under `implicit_drag`), which is
      ! also why `implicit_drag` + `hbbl > 0` is accepted under it.
      if (cfg%ocean%vdiff%bbl_glue) then
         if (.not. cfg%ocean%vdiff%hvel_mom6) then
            call logger%error("ocean_vdiff bbl_glue requires hvel_mom6=.true. — the "// &
                              "botfn glue reads the height-above-bed stack that only "// &
                              "the hvel_mom6 face-thickness build accumulates")
            has_error = .true.
         end if
      end if
      ! `implicit_stress` injects the wind stress at the surface row (k=nz)
      ! only — it cannot represent the DIRECT_STRESS distribution of stress
      ! over the top `hmix_stress` metres (`vmix%direct_stress`).  Enabling
      ! both would silently drop the distribution (the fold wins).  Symmetric
      ! to the implicit_drag/HBBL guard above.  Fail loud.
      if (cfg%ocean%vdiff%implicit_stress .and. cfg%ocean%vmix%direct_stress) then
         call logger%error("ocean_vdiff implicit_stress is incompatible with "// &
                           "distributed wind stress (ocean_vmix direct_stress): the "// &
                           "surface-row fold cannot spread stress over hmix_stress")
         has_error = .true.
      end if

      ! RETIRED `correction_h_weighted`: the h-weighted barotropic-
      ! correction fold is not energy-conserving on any column whose
      ! (open) layers differ in thickness — which is every column of a
      ! stretched z_fixed stack — and MOM6 has no such fold.  Refused,
      ! never silently ignored, so a namelist that relied on it learns
      ! that its answer changes.  See `apply_bt_correction`.
      if (cfg%ocean%bt%correction_h_weighted) then
         ! Pushed to the error ring as well as logged, so a C/Python
         ! caller (and `test_ocean_bt_correction_weight`) reads the
         ! specific reason, not only the generic validation rollup.
         block
            character(len=*), parameter :: msg = &
                                           "&ocean_bt_nml correction_h_weighted is RETIRED: the h-weighted "// &
                                           "barotropic-correction fold was energy-non-conserving (beyond the "// &
                                           "barotropic KE change it adds a positive source 0.5*D^2*H*(kappa-1), "// &
                                           "kappa = sum(h^3)sum(h)/sum(h^2)^2 >= 1, plus a shear feedback that grew "// &
                                           "stretched z_fixed runs non-finite). MOM6 has no h-weighted fold; the "// &
                                           "uniform fold is the default and MOM6's own one (accel_layer_u applies "// &
                                           "the BT acceleration uniformly). For drag-aware damping use "// &
                                           "visc_rem_chain=.true. (with &ocean_vdiff_nml implicit_drag or "// &
                                           "bbl_glue) -- it does NOT re-weight this fold (correction_visc_rem, "// &
                                           "the knob that used to, is itself retired), it feeds visc_rem into "// &
                                           "bt_rem/av_rem instead. Else delete the key."
            call error_ring_push(msg)
            call logger%error(msg)
         end block
         has_error = .true.
      end if
      ! The vdiff operator's row sums are exactly 1 (a no-flux-top/no-flux-
      ! bottom viscous operator cannot remove a uniform acceleration)
      ! UNLESS the implicit-drag fold breaks the k=1 row sum.  So without
      ! `implicit_drag` (or `bbl_glue`), the remnant producer still runs
      ! but returns gamma ≡ 1 identically, and every `*_visc_rem` consumer
      ! degenerates to its own no-op (the forcing weight to the plain
      ! h-mean, the continuity renormaliser to the uniform `du`, av_rem to
      ! 1).  Legal and mathematically correct — merely inert.  Warn (not
      ! error): same "enabled but inert" precedent as the tidal-mixing
      ! e_uniform=0 warning.  PR-3 (D1) follow-up: `forcing_visc_rem`/
      ! `renorm_visc_rem`/`bt_rem_from_visc_rem` are each now
      ! SELF-SUFFICIENT — the producer is decoupled from the retired
      ! `correction_visc_rem` weighted fold and instead runs whenever ANY
      ! of them (or `visc_rem_chain`) is on, so none of them "requires" a
      ! separate producer knob any more.
      if (ocean_bt_visc_rem_producer_on(cfg) .and. .not. &
          (cfg%ocean%vdiff%implicit_drag .or. cfg%ocean%vdiff%bbl_glue)) then
         call logger%warning("&ocean_bt_nml forcing_visc_rem/renorm_visc_rem/"// &
                             "bt_rem_from_visc_rem/visc_rem_chain is on but neither "// &
                             "&ocean_vdiff_nml implicit_drag nor bbl_glue is: the vdiff "// &
                             "operator then carries no drag, so visc_rem = 1 identically "// &
                             "and every consumer is a no-op")
      end if
      ! The bc-PGF retro-correction (MOM6 btstep_layer_accel) builds its
      ! per-layer `pbce` from the FV_MOM6 interface-height stack
      ! `pgf%e_face`, which no other form fills.  It used to be ACCEPTED here
      ! and `error stop` inside step 1 (`compute_pbce`); refuse it at
      ! configure, naming the same requirement.  Not generalised: the
      ! MONT / FV_LITE / FV_WRIGHT forms are surface-relative (they carry no
      ! free-surface term, `g_pf = 0`), so the response `pbce = dp_k/deta`
      ! the correction redistributes is not the one their PGF sees, and
      ! GPRIME runs the fast loop at the reduced `g_FS`.
      if (cfg%ocean%bt%correction_bc_pgf .and. &
          parse_opgf_variant(cfg%ocean%pgf%form) /= OPGF_VARIANT_FV_MOM6) then
         call logger%error("&ocean_bt_nml correction_bc_pgf=.true. requires "// &
                           "&ocean_pgf_nml form='fv_mom6' (got '"// &
                           trim(adjustl(cfg%ocean%pgf%form))//"'). compute_pbce "// &
                           "builds the per-layer pressure response from the FV_MOM6 "// &
                           "interface-height stack (pgf%e_face), which no other PGF "// &
                           "form fills.")
         has_error = .true.
      end if
      ! pred_corr (SPEC S4) envelope: any ALE / Lagrangian-within-step
      ! vcoord (lagrangian, sigma, zstar, zsigma, zstar_sigma,
      ! zstar_full) — MOM6's own model: the dynamics step is always
      ! Lagrangian, regridding is orthogonal, and the pc corrector flows
      ! into the same thermo-cadence ALE remap the ssp path uses.  The
      ! LEGACY pure Eulerian-z path stays excluded: its per-stage
      ! vertical advection + BT-fold h-rescale would compose differently
      ! inside the restructured loop (the fold's rescale lands BEFORE the
      ! deferred continuity — a double eta application) — untested, fail
      ! loud.  Also still excluded: dynamic wet/dry (per-stage masking
      ! composes with the ssp stages only) and windowed tracer advection
      ! (the predictor's TR_MODE_NONE + window accounting is unexercised).
      if (trim(cfg%ocean%bt%split_scheme) == "pred_corr") then
         if (parse_vcoord_type(cfg%vcoord_type, &
                               default_code=VCOORD_EULERIAN_Z) == VCOORD_EULERIAN_Z) then
            call logger%error("&ocean_bt_nml split_scheme='pred_corr' requires an ALE "// &
                              "vertical coordinate (lagrangian/sigma/zstar/zsigma/"// &
                              "zstar_sigma/zstar_full) — the legacy 'eulerian_z' "// &
                              "per-stage vertical-advection + h-rescale path is not "// &
                              "wired through the restructured pc loop; got '"// &
                              trim(cfg%vcoord_type)//"'. pred_corr is the DEFAULT, so "// &
                              "set split_scheme='ssp_rk2' to keep this configuration.")
            has_error = .true.
         end if
         if (cfg%ocean%wetdry%enable) then
            call logger%error("&ocean_bt_nml split_scheme='pred_corr' is incompatible "// &
                              "with &ocean_wetdry_nml enable (v1 envelope). pred_corr "// &
                              "is the DEFAULT, so set split_scheme='ssp_rk2' to keep "// &
                              "this configuration.")
            has_error = .true.
         end if
         if (cfg%ocean%vmix%dt_tracer_advect_ratio > 1) then
            call logger%error("&ocean_bt_nml split_scheme='pred_corr' requires "// &
                              "dt_tracer_advect_ratio = 1 (v1 envelope). pred_corr is "// &
                              "the DEFAULT, so set split_scheme='ssp_rk2' to keep this "// &
                              "configuration.")
            has_error = .true.
         end if
         ! z-level closed faces need their SOLID WALLS to be walls.  With
         ! `mask_wall_velocity = .false.` a wall face keeps a free per-layer
         ! velocity that carries no mass (continuity zeroes the wall flux)
         ! and whose depth mean the BT fold resets to zero every stage —
         ! but whose BAROCLINIC part nothing restores under pred_corr: the
         ! Coriolis reads `u_av`, which the transport renormaliser never
         ! writes at a skipped wall face, so the wall velocity never sees
         ! its own rotation and integrates the layer Coriolis of its
         ! interior neighbour without bound.  The closed-face mask is what
         ! makes that neighbour baroclinic (a closed bed layer beside open
         ! ones).  Measured on the rotating ledge basin of
         ! tests/test_ocean_zfixed_cor_ref.F90: KE+PE x1783 over 4000
         ! outer steps with the physical interior at x1.24 — the whole
         ! growth sits on the unmasked wall faces — against x0.74 with
         ! the walls masked.  The mask is the DEFAULT; refused rather than
         ! silently overridden.
         if (cfg%zfixed_closed_faces .and. .not. cfg%ocean%bc%mask_wall_velocity) then
            call logger%error("&ocean_bt_nml split_scheme='pred_corr' with "// &
                              "&vcoord_nml zfixed_closed_faces=.true. requires "// &
                              "&ocean_bc_nml mask_wall_velocity=.true. (the default): "// &
                              "an unmasked solid-wall face carries a baroclinic layer "// &
                              "velocity that pred_corr's u_av-evaluated Coriolis never "// &
                              "rotates, so it grows without bound. Set "// &
                              "mask_wall_velocity=.true., or split_scheme='ssp_rk2'.")
            has_error = .true.
         end if
      end if
      ! RETIRED `accel_visc_rem`: PR-3's audit found no MOM6 state-update
      ! equivalent — `btstep_layer_accel` and the split-explicit RK2
      ! corrector's `up`/`vp` update both apply the depth-mean
      ! barotropic acceleration `u_accel_bt` UNIFORMLY across every
      ! layer — no `visc_rem` weight anywhere in that path (the only
      ! `visc_rem x u_accel_bt` products in MOM6 are a diagnostic-only
      ! `id_u_BT_accel_visc_rem` post-product, never fed back into state,
      ! and `RESCALE_STRONG_DRAG`'s depth-MEAN, not per-layer, rescale —
      ! already its own knob).  The real MOM6 mechanisms that multiply a
      ! velocity correction by `visc_rem` are `renorm_visc_rem` (continuity
      ! `u_cor = u + du*visc_rem`) and `rescale_strong_drag`.  Refused,
      ! never silently ignored.
      if (cfg%ocean%vdiff%accel_visc_rem) then
         block
            character(len=*), parameter :: msg = &
                                           "&ocean_vdiff_nml accel_visc_rem is RETIRED: PR-3's audit found no "// &
                                           "MOM6 state-update equivalent -- MOM6's btstep_layer_accel applies "// &
                                           "the depth-mean barotropic acceleration uniformly across every "// &
                                           "layer, with no visc_rem weight. The real MOM6 mechanisms that "// &
                                           "multiply a velocity correction by visc_rem are "// &
                                           "&ocean_bt_nml renorm_visc_rem (continuity u_cor = u + du*visc_rem) "// &
                                           "and rescale_strong_drag. Use one of those, or delete the key."
            call error_ring_push(msg)
            call logger%error(msg)
         end block
         has_error = .true.
      end if
      ! RETIRED `correction_visc_rem`: MOM6's `accel_layer_u` gives every
      ! layer the SAME `u_accel_bt` plus only the depth-mean-zero `pbce`
      ! baroclinic-pressure term — NO `visc_rem` weight — and that
      ! unweighted acceleration is folded into `up` BEFORE `vertvisc`
      ! (in the split-explicit RK2 corrector), so the glue's implicit
      ! friction (which already includes the BBL
      ! piston) is what then distributes it across layers — ONCE, not
      ! twice.  roundabout's `correction_visc_rem` applies
      ! `apply_bt_correction` BEFORE that same implicit friction
      ! (`vmix_apply_in_stage`) and weights it by `visc_rem_k/⟨visc_rem⟩_h`
      ! on top — a SECOND, unbounded-ratio damping that concentrates the
      ! correction into whichever layers the glue left least damped.
      ! Measured on the 1-degree Southern Ocean z* OPEN-step probe
      ! (`probes/pr3_producer_only`): `hvel_mom6`+`bbl_glue`+
      ! `implicit_drag` alone runs clean; adding ONLY `correction_visc_rem`
      ! NaNs at step 38.  Under D1 ("exactly MOM6's set") this fold is not
      ! part of the chain and has no standalone namelist path of its own;
      ! refused, never silently ignored.  The kernel itself
      ! (`apply_bt_correction`'s `use_visc_rem` dispatch) and its direct
      ! unit tests (`tests/test_ocean_bt_correction_weight.F90`
      ! `visc_rem_unity_is_uniform`/`visc_rem_biases_against_bed`/
      ! `closed_faces_open_column`) are untouched — they call it with a
      ! raw logical, not through `cfg`.
      if (cfg%ocean%bt%correction_visc_rem) then
         block
            character(len=*), parameter :: msg = &
                                           "&ocean_bt_nml correction_visc_rem is RETIRED: MOM6's accel_layer_u "// &
                                           "applies the barotropic acceleration UNIFORMLY across every layer, "// &
                                           "before vertvisc distributes it via "// &
                                           "the SAME implicit friction the glue uses -- this fold re-weights it "// &
                                           "a second time by visc_rem/<visc_rem>_h, an unbounded ratio that "// &
                                           "NaNs the 1-degree Southern Ocean z* open-step case under bbl_glue "// &
                                           "at step ~40. Use visc_rem_chain (producer + bt_rem_from_av_rem + "// &
                                           "wt_u forcing + renorm_visc_rem, uniform BT-correction fold) instead, "// &
                                           "or delete the key."
            call error_ring_push(msg)
            call logger%error(msg)
         end block
         has_error = .true.
      end if
      ! D2: `substep_drag`'s linear piston and the visc_rem chain both put
      ! bed drag into bt_rem — composing them double-counts it (once via
      ! the glue/implicit_drag fold inside visc_rem, once via the piston).
      if (ocean_bt_rem_from_visc_rem_on(cfg) .and. cfg%ocean%bt%substep_drag) then
         call logger%error("&ocean_bt_nml bt_rem_from_visc_rem=.true. (or "// &
                           "visc_rem_chain=.true.) is mutually exclusive with "// &
                           "substep_drag=.true. (bed drag would be double-counted: once "// &
                           "inside the visc_rem producer via the bbl_glue/implicit_drag "// &
                           "fold, once again via the linear piston law)")
         has_error = .true.
      end if
      ! Decomposition: av_rem/bt_rem are built on the SAME normal-width
      ! face stencil visc_rem occupies (valid after PR-1's halo refresh);
      ! the wide-halo BT clone's metrics_w/halo-widened arrays carry no
      ! av_rem/visc_rem ghost width yet — same posture as porous.
      if (ocean_bt_rem_from_visc_rem_on(cfg) .and. cfg%ocean%bt%bt_halo > 0) then
         call logger%error("&ocean_bt_nml bt_rem_from_visc_rem=.true. (or "// &
                           "visc_rem_chain=.true.) is mutually exclusive with "// &
                           "bt_halo > 0 (the wide-halo BT clone carries no "// &
                           "av_rem/visc_rem ghost-width statistics, like porous)")
         has_error = .true.
      end if
      ! D3: BT_STRONG_DRAG / RESCALE_STRONG_DRAG are refinements OF the
      ! av_rem chain, not independent knobs.
      if (cfg%ocean%bt%strong_drag .and. .not. ocean_bt_rem_from_visc_rem_on(cfg)) then
         call logger%error("&ocean_bt_nml strong_drag=.true. requires "// &
                           "bt_rem_from_visc_rem=.true. (or visc_rem_chain=.true.) — "// &
                           "the rational-approximation bt_rem form is only defined in "// &
                           "terms of av_rem")
         has_error = .true.
      end if
      if (cfg%ocean%bt%rescale_strong_drag .and. .not. cfg%ocean%bt%strong_drag) then
         call logger%error("&ocean_bt_nml rescale_strong_drag=.true. requires "// &
                           "strong_drag=.true. — the rescale corrects for the rational "// &
                           "form's bt_rem**n_inner /= av_rem gap, which the plain power "// &
                           "form does not have")
         has_error = .true.
      end if
      if (substep_drag_ignores_bdrag_form(cfg)) then
         call logger%warning("&ocean_bt_nml substep_drag=.true. with "// &
                             "&ocean_bdrag_nml form='"// &
                             trim(adjustl(cfg%ocean%bdrag%form))//"': the "// &
                             "barotropic substep damping is built from the "// &
                             "LINEAR coefficient &ocean_bdrag_nml r (times hbbl) "// &
                             "only, so it does NOT follow this bottom drag — "// &
                             "with r = 0 (the default) substep_drag is a no-op, "// &
                             "otherwise it damps the barotropic mode with a "// &
                             "linear drag the slow step never applies.  Use "// &
                             "form='linear', or drop substep_drag")
      end if

      ! Equilibrium tide (C1) requires lat/lon — meaningless on a
      ! cartesian grid (geolatT/geolonT stay 0).  Fail loud rather than
      ! silently forcing a flat basin.
      if (cfg%ocean%tides%enable .and. &
          trim(cfg%ocean%grid%grid_config) == "cartesian") then
         call logger%error("ocean_tides_nml: enable=.true. requires a "// &
                           "non-cartesian grid_config (spherical/tripolar/"// &
                           "supergrid) — the equilibrium tide needs lat/lon")
         has_error = .true.
      end if
      if (cfg%ocean%tides%enable .and. len_trim(cfg%ocean%tides%constituents) == 0) then
         call logger%error("ocean_tides_nml: enable=.true. but constituents "// &
                           "list is empty")
         has_error = .true.
      end if

      ! C7 Henyey latitude factor: copy of the equilibrium-tide cartesian
      ! guard above.  `geolatT` is identically zero on a cartesian grid (only
      ! the spherical/supergrid/tripolar fills populate it), so EVERY column
      ! would take the equatorial factor L(0 deg) = 0 and the background
      ! would collapse to a uniform `bkgnd_kd_min` everywhere while the
      ! console cheerfully logged "Henyey IGW latitude factor ON".  The
      ! `kd_min` floor makes that less catastrophic than the unfloored
      ! zeroing it used to be, but it is still a latitude parameterisation
      ! on a grid with no meaningful latitude — strictly worse than leaving
      ! the knob off, so refuse at configure.
      if (cfg%ocean%vmix%bkgnd_henyey .and. &
          trim(cfg%ocean%grid%grid_config) == "cartesian") then
         call logger%error("&ocean_vmix_nml bkgnd_henyey=.true. requires a "// &
                           "non-cartesian &ocean_grid_nml grid_config (spherical/"// &
                           "tripolar/supergrid) — geolatT is meaningless on cartesian, "// &
                           "so every column would take the equatorial L(0 deg)=0 and "// &
                           "the background would flatten to a uniform bkgnd_kd_min")
         has_error = .true.
      end if
      ! Bryan-Lewis and Henyey are MUTUALLY EXCLUSIVE background schemes,
      ! matching the reference code, which FATALs when a second background
      ! scheme is selected.  Henyey scales the SCALAR background; letting it
      ! also multiply the Bryan-Lewis deep asymptote would suppress abyssal /
      ! internal-tide mixing at the equator, which is not what the Henyey
      ! scaling describes.
      if (bkgnd_henyey_conflicts_profile(cfg%ocean%vmix%bkgnd_henyey, &
                                         cfg%ocean%vmix%bkgnd_profile)) then
         call logger%error("&ocean_vmix_nml bkgnd_henyey and bkgnd_profile are "// &
                           "mutually exclusive background schemes — pick ONE "// &
                           "(Bryan-Lewis depth profile, or the Henyey latitude "// &
                           "factor on the scalar background)")
         has_error = .true.
      end if
      ! Surface-pressure loading / inverse barometer (PR-17).  The seam it
      ! folds into (`eta_forcing`) lives in the barotropic substep, so the
      ! feature is split-solver only; the `p_surf` field it reads is
      ! allocated only with the PR-12 component set.  Every violation fails
      ! loud rather than silently reading an unallocated array (a silent
      ! stale device read on the mem:separate GPU build).
      if (cfg%ocean%psurf%enable) then
         if (.not. cfg%ocean%forcing%enable_components) then
            call logger%error("&ocean_psurf_nml enable=.true. requires "// &
                              "&ocean_forcing_nml enable_components=.true. "// &
                              "(p_surf is allocated only with the component set)")
            has_error = .true.
         end if
         if (cfg%ocean%bt%n_inner < 1 .and. .not. cfg%ocean%bt%auto_n_inner) then
            call logger%error("&ocean_psurf_nml enable=.true. requires the "// &
                              "split-explicit solver (&ocean_bt_nml n_inner >= 1 "// &
                              "or auto_n_inner=.true.) — the eta_forcing seam "// &
                              "lives in the barotropic substep")
            has_error = .true.
         end if
         ! A UNIFORM surface pressure has no gradient => provably inert
         ! (gauge invariance).  A user enabling p_surf with only a uniform
         ! constant and no file/override path gets nothing — say so rather
         ! than repeat the &ocean_tidal_mixing_nml e_uniform silent-no-op.
         if (cfg%ocean%psurf%p_surf_const /= 0.0_wp) then
            if (cfg%ocean%psurf%in_eos) then
               ! With in_eos the gauge argument does NOT apply to the EOS:
               ! it is nonlinear in pressure, so a spatially UNIFORM load
               ! still changes the IN-SITU densities it reaches.  Say so
               ! instead of the (then wrong) "provably inert" warning.
               call logger%info("&ocean_psurf_nml in_eos=.true.: a uniform "// &
                                "p_surf_const is inert on the barotropic seam "// &
                                "(only grad(p_surf) is physical there) but NOT "// &
                                "in the IN-SITU EOS pressure — it compresses "// &
                                "the water. The potential density ms%rho_layer "// &
                                "is unaffected either way (uniform p_ref).")
            else
               call logger%warning("&ocean_psurf_nml enable=.true. with a uniform "// &
                                   "p_surf_const and no file/override path: only "// &
                                   "grad(p_surf) is physical, so a spatially "// &
                                   "constant load is inert (gauge invariance). "// &
                                   "File-driven p_surf lands in PR-14/PR-15.")
            end if
         end if
      end if
      ! ---- Top-of-column pressure in the IN-SITU EOS arguments (E3) ----
      ! `in_eos` offsets the EOS's IN-SITU pressure arguments by
      ! `ms%p_top`.  v1 ports exactly one: the FV_WRIGHT Picard column
      ! sweep.  It does NOT touch `ms%rho_layer`, which is a potential
      ! density at the horizontally uniform `&ocean_eos_nml p_ref` and must
      ! stay that way.  Every OTHER in-situ consumer builds its own
      ! surface-relative hydrostatic pressure starting at 0 Pa at the free
      ! surface and has NOT been ported; running them against a loaded
      ! column would mix two incompatible pressure conventions inside one
      ! time step with no symptom.  Refuse fail-loud rather than be
      ! silently inconsistent — each line below is a named follow-up, not a
      ! permanent limit.
      if (cfg%ocean%psurf%in_eos) then
         ! Inert-configuration warning, not a refusal: with no in-situ
         ! consumer selected the knob legitimately does nothing, and the
         ! house rule (cf. &ocean_tidal_mixing_nml e_uniform) is to SAY so
         ! rather than let a user believe a cavity load reached the EOS.
         ! E4 adds a THIRD ported consumer: `buoyancy_coeffs="eos"` seeds
         ! the KPP B_0 coefficients at `p_top` and the double-diffusion
         ! interface stack from it, so the knob is no longer inert when
         ! that is selected.
         if (trim(adjustl(cfg%ocean%pgf%form)) /= "fv_wright" .and. &
             .not. cfg%ocean%epbl%enable .and. &
             trim(adjustl(cfg%ocean%vmix%buoyancy_coeffs)) /= "eos") then
            call logger%warning("&ocean_psurf_nml in_eos=.true. is INERT for "// &
                                "&ocean_pgf_nml form='"// &
                                trim(adjustl(cfg%ocean%pgf%form))//"' with no "// &
                                "other ported in-situ consumer: the ported ones "// &
                                "are the FV_WRIGHT Picard column sweep, the "// &
                                "EPBL column stack (&ocean_epbl_nml) and the "// &
                                "EOS-derived buoyancy coefficients "// &
                                "(&ocean_vmix_nml buoyancy_coeffs='eos'). The "// &
                                "other PGF forms read the POTENTIAL density "// &
                                "ms%rho_layer, which is referenced to the "// &
                                "uniform &ocean_eos_nml p_ref BY DESIGN and is "// &
                                "not offset by the load.")
         end if
         if (.not. cfg%ocean%psurf%enable) then
            call logger%error("&ocean_psurf_nml in_eos=.true. requires "// &
                              "enable=.true. (p_top is filled from the "// &
                              "assembled sf%p_surf the seam owns)")
            has_error = .true.
         end if
         ! EPBL was refused here until Phase 4b.  It is PORTED now:
         ! `epbl_column_kernel` seeds its stack at `ms%p_top(i,j)` when
         ! `in_eos`, which moves BOTH consumers of that stack (the
         ! in-situ `eos_specvol_derivs` argument and the PE weight
         ! `dmass*p_mid*dsv`) together.  Gate:
         ! `test_ocean_bl_under_ice` (`epbl_p_top_*`).
         if (cfg%ocean%kshear%enable) then
            call logger%error("&ocean_psurf_nml in_eos=.true. is not supported "// &
                              "with &ocean_kappa_shear_nml enable=.true. — "// &
                              "ks_solve_column builds its own interface pressure "// &
                              "from 0 Pa at the surface; not yet ported to p_top")
            has_error = .true.
         end if
         if (cfg%ocean%tidal_mixing%enable) then
            call logger%error("&ocean_psurf_nml in_eos=.true. is not supported "// &
                              "with &ocean_tidal_mixing_nml enable=.true. — "// &
                              "tidal_mixing_column_kernel builds its own interface "// &
                              "pressure from 0 Pa at the surface; not yet ported")
            has_error = .true.
         end if
         if (cfg%ocean%redi%enable) then
            call logger%error("&ocean_psurf_nml in_eos=.true. is not supported "// &
                              "with &ocean_redi_nml enable=.true. — "// &
                              "redi_build_column seeds Pint(1)=0 at the surface; "// &
                              "not yet ported to p_top")
            has_error = .true.
         end if
         if (cfg%ocean%slopes%enable) then
            call logger%error("&ocean_psurf_nml in_eos=.true. is not supported "// &
                              "with &ocean_slopes_nml enable=.true. — "// &
                              "pressure_above_x sums g*rho0*h down from 0 Pa at "// &
                              "the surface; not yet ported to p_top")
            has_error = .true.
         end if
         if (cfg%ocean%pgf%reconstruct_for_pressure) then
            call logger%error("&ocean_psurf_nml in_eos=.true. is not supported "// &
                              "with &ocean_pgf_nml reconstruct_for_pressure=.true. "// &
                              "— boole_dpa_intz_layer builds the EOS pressure as "// &
                              "p = -g*rho0*z from the surface-relative interface "// &
                              "height; not yet ported to p_top")
            has_error = .true.
         end if
         if (trim(adjustl(cfg%ocean%pgf%form)) == "fv_mom6" .and. &
             .not. cfg%ocean%pgf%reconstruct_for_pressure .and. &
             cfg%ocean%pgf%insitu_density .and. &
             trim(adjustl(cfg%ocean%eos%eos)) /= "linear") then
            call logger%error("&ocean_psurf_nml in_eos=.true. is not supported "// &
                              "with the FV_MOM6 in-situ PCM density (&ocean_pgf_nml "// &
                              "insitu_density=.true., the default, under a "// &
                              "pressure-dependent EOS) — boole_dpa_intz_layer builds "// &
                              "the EOS pressure as p = -g*rho0*z; not yet ported to "// &
                              "p_top.  Set insitu_density=.false. for the legacy "// &
                              "potential-density integral")
            has_error = .true.
         end if
         if (cfg%ocean%ice%enable) then
            call logger%error("&ocean_psurf_nml in_eos=.true. is not supported "// &
                              "with &ocean_ice_nml enable=.true. — the freezing "// &
                              "point is evaluated at p=0 in the frazil / basal-flux "// &
                              "kernels, and the liquidus pressure depression is "// &
                              "exactly what an ice-shelf load changes; not yet ported")
            has_error = .true.
         end if
      end if
      ! ---- Z-level T/S initial-condition overlay (`&ocean_zinit_nml`) ----
      if (cfg%ocean%zinit%enable) then
         select case (trim(adjustl(cfg%ocean%zinit%source)))
         case ("file")
            if (len_trim(cfg%ocean%zinit%file) == 0) then
               call logger%error("&ocean_zinit_nml enable=.true. with "// &
                                 "source='file' requires a non-blank file= path.")
               has_error = .true.
            end if
         case ("linear")
            ! An analytic profile plus a file path is ambiguous: the file
            ! would be silently ignored.  Say so rather than pick one.
            if (len_trim(cfg%ocean%zinit%file) > 0) then
               call logger%error("&ocean_zinit_nml source='linear' takes the "// &
                                 "analytic lin_* profile and opens NOTHING, so the "// &
                                 "file='"//trim(adjustl(cfg%ocean%zinit%file))// &
                                 "' you also set would be silently ignored.  Pick "// &
                                 "one: drop file=, or set source='file'.")
               has_error = .true.
            end if
         case default
            call logger%error("&ocean_zinit_nml source='"// &
                              trim(adjustl(cfg%ocean%zinit%source))// &
                              "' is not a known profile source; expected 'file' "// &
                              "(pre-regridded NetCDF) or 'linear' (analytic "// &
                              "affine T(z)/S(z)).")
            has_error = .true.
         end select
      end if
      ! ---- Top-of-column load in the PGF surface boundary condition (P5.0) ----
      ! `pa(nz+1) = rho_ref*g*eta_geo + ms%p_top`.  Only the FV_MOM6 family
      ! builds a `pa` stack at all — `mont` hard-zeroes `M(nz)` and
      ! `fv_lite`/`fv_wright` seed `p_edge(nz+1) = 0` — so there is
      ! literally no boundary condition to inject anywhere else.  Refuse
      ! rather than accept a knob that would silently do nothing on a form
      ! the user believes is carrying an ice load.
      if (cfg%ocean%pgf%p_top_in_bc) then
         if (trim(adjustl(cfg%ocean%pgf%form)) /= "fv_mom6") then
            call logger%error("&ocean_pgf_nml p_top_in_bc=.true. requires "// &
                              "form='fv_mom6' (got '"// &
                              trim(adjustl(cfg%ocean%pgf%form))//"'). Only the "// &
                              "FV_MOM6 family builds the pa(nz+1) pressure-stack "// &
                              "boundary condition the load is injected into; mont "// &
                              "hard-zeroes M(nz) and fv_lite/fv_wright seed "// &
                              "p_edge(nz+1)=0.")
            has_error = .true.
         end if
         ! Inert-configuration warning, not a refusal (the house rule, cf.
         ! &ocean_psurf_nml in_eos and &ocean_tidal_mixing_nml e_uniform).
         ! `ms%p_top = metrics%p_ice_ref + sf%p_surf` has TWO producers,
         ! and the warning must name both or it is a lie: the
         ! &ocean_psurf_nml seam (the atmospheric half, `sf%p_surf`) and
         ! the &ocean_cavity_dyn_nml ice-shelf load (the static half,
         ! `p_ice_ref = rho_ref*g*z_draft`, assembled in
         ! `configure_ocean_cavity`).  Under a cavity `p_top` carries the
         ! ice load and `p_top_in_bc` is not merely live — it is REQUIRED
         ! for a varying draft, refused above and again at configure.  The
         ! warning fires only when NEITHER producer is on, which is the
         ! one case in which `p_top` really is the zero array it was
         ! allocated as.
         if (.not. p_top_has_producer(cfg)) then
            call logger%warning("&ocean_pgf_nml p_top_in_bc=.true. is INERT "// &
                                "without a producer for ms%p_top: enable "// &
                                "&ocean_psurf_nml (the atmospheric surface-pressure "// &
                                "seam) or &ocean_cavity_dyn_nml (the static "// &
                                "ice-shelf load), else p_top is the zero array and "// &
                                "pa(nz+1) is unchanged.")
         end if
      end if
      ! ---- Static ice-shelf cavity geometry (&ocean_cavity_dyn_nml, P5.1) ----
      ! The draft is absorbed into the barotropic DATUM (bt_H_ref =
      ! b - z_draft), so every consumer of the water-column thickness
      ! D = bt_H_ref + bt_eta is correct with no cavity branch of its own.
      ! What that buys is paid for by a narrow envelope, and EVERY
      ! restriction below fails loud naming the knob and the reason: a
      ! cavity that silently runs outside it looks plausible and is wrong
      ! (a coordinate anchored at z = 0 under 500 m of ice, a second
      ! un-reconciled surface load, a wide-halo BT clone with no draft).
      if (cfg%ocean%cavity_dyn%trim_ic_for_p_surf .and. &
          .not. cfg%ocean%cavity_dyn%enable) then
         call logger%error("&ocean_cavity_dyn_nml trim_ic_for_p_surf=.true. "// &
                           "requires enable=.true. (there is no ice load to "// &
                           "trim the initial column against)")
         has_error = .true.
      end if
      if (cfg%ocean%cavity_dyn%enable) then
         if (trim(cfg%sim_type) /= "ocean") then
            call logger%error("&ocean_cavity_dyn_nml enable=.true. requires "// &
                              "sim_type='ocean'")
            has_error = .true.
         end if
         ! --- geometry source envelope ---
         select case (trim(adjustl(cfg%ocean%cavity_dyn%draft_config)))
         case ("none", "flat", "linear")
            continue
         case ("file")
            ! Ships single-rank, through the PR-14 static-2-D reader
            ! (`ocean_data_input_load_static_2d`), which DOES apply the
            ! global offset — but the cavity as a whole is single-rank
            ! fenced below, and the loader re-asserts it.
            if (len_trim(cfg%ocean%cavity_dyn%draft_file) == 0) then
               call logger%error("&ocean_cavity_dyn_nml draft_config='file' requires "// &
                                 "draft_file")
               has_error = .true.
            end if
            if (len_trim(cfg%ocean%cavity_dyn%draft_var) == 0) then
               call logger%error("&ocean_cavity_dyn_nml draft_config='file' requires "// &
                                 "draft_var (the 2-D variable name; ISOMIP+ ships "// &
                                 "'iceDraft')")
               has_error = .true.
            end if
            if (parse_cavity_draft_sign(cfg%ocean%cavity_dyn%draft_sign) == &
                CAVITY_SIGN_INVALID) then
               call logger%error("&ocean_cavity_dyn_nml draft_sign='"// &
                                 trim(adjustl(cfg%ocean%cavity_dyn%draft_sign))// &
                                 "' is not recognised (depth|positive_down|"// &
                                 "elevation|positive_up).  There is no default that "// &
                                 "guesses from the data: the ISOMIP+ file carries an "// &
                                 "ELEVATION (z_d <= 0) and a depth file carries "// &
                                 "z_draft >= 0, and the two differ by the whole load.")
               has_error = .true.
            end if
         case default
            call logger%error("&ocean_cavity_dyn_nml draft_config='"// &
                              trim(adjustl(cfg%ocean%cavity_dyn%draft_config))// &
                              "' is not recognised (none|flat|linear|file)")
            has_error = .true.
         end select
         select case (trim(adjustl(cfg%ocean%cavity_dyn%draft_source)))
         case ("draft", "thickness")
            continue
         case ("in_situ")
            call logger%error("&ocean_cavity_dyn_nml draft_source='in_situ' (true "// &
                              "isostasy, p_ice = g*int(rho)) is not implemented: it "// &
                              "needs a per-column root find and does NOT admit exact "// &
                              "discrete rest in the split solver.  Use 'draft' "// &
                              "(the Boussinesq-isostatic flotation load ISOMIP+ "// &
                              "prescribes) or 'thickness'.")
            has_error = .true.
         case default
            call logger%error("&ocean_cavity_dyn_nml draft_source='"// &
                              trim(adjustl(cfg%ocean%cavity_dyn%draft_source))// &
                              "' is not recognised (draft|thickness|in_situ)")
            has_error = .true.
         end select
         ! `"linear"` measures its profile from `draft_x0`, so an
         ! unbounded (sentinel) anchor would make `draft_depth` meaningless
         ! and the whole shelf depth an artefact of 1e30*slope.
         if (trim(adjustl(cfg%ocean%cavity_dyn%draft_config)) == "linear" .and. &
             abs(cfg%ocean%cavity_dyn%draft_x0) >= 1.0e29_wp) then
            call logger%error("&ocean_cavity_dyn_nml draft_config='linear' requires "// &
                              "a finite draft_x0: it is the ANCHOR of the profile "// &
                              "(draft_depth is the draft AT draft_x0), not just the "// &
                              "western edge of the box.")
            has_error = .true.
         end if
         if (cfg%ocean%cavity_dyn%draft_depth < 0.0_wp) then
            call logger%error("&ocean_cavity_dyn_nml draft_depth must be >= 0 "// &
                              "(it is a DEPTH below z = 0, positive down)")
            has_error = .true.
         end if
         if (cfg%ocean%cavity_dyn%h_min_cavity <= 0.0_wp) then
            call logger%error("&ocean_cavity_dyn_nml h_min_cavity must be > 0 "// &
                              "(the grounding cutoff; 0 would admit a zero-thickness "// &
                              "water column under the ice)")
            has_error = .true.
         end if
         if (cfg%ocean%cavity_dyn%grounded_max_frac <= 0.0_wp .or. &
             cfg%ocean%cavity_dyn%grounded_max_frac > 1.0_wp) then
            call logger%error("&ocean_cavity_dyn_nml grounded_max_frac must be in "// &
                              "(0, 1] (the fraction of interior columns allowed to "// &
                              "ground)")
            has_error = .true.
         end if
         if (trim(adjustl(cfg%ocean%cavity_dyn%draft_source)) == "thickness" .and. &
             cfg%ocean%cavity_dyn%rho_ice <= 0.0_wp) then
            call logger%error("&ocean_cavity_dyn_nml draft_source='thickness' "// &
                              "requires rho_ice > 0")
            has_error = .true.
         end if
         ! MOM6 TRIM_IC_FOR_P_SURF.  The trim depth solves
         ! g*int_{-s}^{0} rho dz = p_ice_ref in CLOSED FORM, which needs a
         ! density that is affine in z above the ice base: the linear EOS
         ! over the analytic linear zinit profile.  A nonlinear EOS or a
         ! file profile would need a per-column root find against the
         ! column's own extrapolated T/S (MOM6 cut_off_column_top) and is
         ! not wired; a uniform_z seed lays interfaces from z = 0, not from
         ! the (trimmed) ice base.
         if (cfg%ocean%cavity_dyn%trim_ic_for_p_surf) then
            if (trim(adjustl(cfg%ocean%eos%eos)) /= "linear") then
               call logger%error("&ocean_cavity_dyn_nml trim_ic_for_p_surf=.true. "// &
                                 "requires &ocean_eos_nml eos='linear' (the trim "// &
                                 "depth is the closed-form root for a density "// &
                                 "affine in z; a nonlinear-EOS trim is not wired)")
               has_error = .true.
            end if
            if (.not. cfg%ocean%zinit%enable .or. &
                trim(adjustl(cfg%ocean%zinit%source)) /= "linear") then
               call logger%error("&ocean_cavity_dyn_nml trim_ic_for_p_surf=.true. "// &
                                 "requires &ocean_zinit_nml enable=.true., "// &
                                 "source='linear': the analytic T(z)/S(z) profile "// &
                                 "is what defines the density of the water the "// &
                                 "ice displaces")
               has_error = .true.
            end if
            if (trim(cfg%thickness_config) == "uniform_z") then
               call logger%error("&ocean_cavity_dyn_nml trim_ic_for_p_surf=.true. "// &
                                 "is incompatible with thickness_config='uniform_z'")
               has_error = .true.
            end if
         end if
         ! --- the one atmospheric-forcing path the cover mask does NOT
         !     reach (P2c) ---
         ! Every static forcing field is masked: the wind pair and the
         ! scalar q_heat/q_salt once at configure, the component bands
         ! every thermo step in the assembler.  The FILE-DRIVEN override
         ! is the exception: `ocean_data_forcing_apply` rewrites `tau_x`/
         ! `tau_y` (and, with `heat_to_component=.false.`, `Q_heat`
         ! itself) from the next time bracket with no access to
         ! `metrics%cover_frac` — it is handed `ss`, `sf`, `grid` and
         ! `bc`, and nothing else.  Re-masking per bracket means
         ! threading the metrics slot through the reader, which is a
         ! separate change.  Refused rather than half-wired: a cavity run
         ! whose wind is silently restored to its unmasked file value on
         ! the first bracket read looks entirely plausible and is wrong.
         if (cfg%ocean%dataovr%enable) then
            call logger%error("&ocean_cavity_dyn_nml enable=.true. is mutually "// &
                              "exclusive with &ocean_dataovr_nml enable=.true.  The "// &
                              "ice-cover mask on the atmospheric forcing is applied "// &
                              "to the wind pair at configure and to the surface-flux "// &
                              "components in the assembler; the data-override reader "// &
                              "rewrites tau_x/tau_y (and Q_heat, unless "// &
                              "heat_to_component=.true.) per time bracket without "// &
                              "the cover, which would restore the unmasked "// &
                              "atmosphere under the shelf.  Follow-up: thread "// &
                              "cover_frac through ocean_data_forcing_apply and the "// &
                              "ocean_seam_refresh_surface_stress seam.")
            has_error = .true.
         end if
         ! --- pressure-gradient envelope ---
         if (trim(adjustl(cfg%ocean%pgf%form)) /= "fv_mom6") then
            call logger%error("&ocean_cavity_dyn_nml enable=.true. requires "// &
                              "&ocean_pgf_nml form='fv_mom6' (got '"// &
                              trim(adjustl(cfg%ocean%pgf%form))//"'). Only the "// &
                              "FV_MOM6 family builds the pa(nz+1) pressure-stack "// &
                              "boundary condition the ice load is injected into; "// &
                              "mont hard-zeroes M(nz) and fv_lite/fv_wright seed "// &
                              "p_edge(nz+1)=0.")
            has_error = .true.
         end if
         ! P5.2 — the load must have a consumer once it has a GRADIENT.
         ! `p_top_in_bc` is the only route by which the isostatic load
         ! rho_ref*g*z_draft reaches the FV_MOM6 pa(nz+1) surface BC;
         ! without it a varying draft leaves the pressure stack ~5e6 Pa
         ! off its anomaly scale, the unsplit driver feels a raw
         ! g*grad(z_draft), and a non-uniform BT-correction weight turns the
         ! uncancelled depth-uniform force into a real per-layer shear.
         ! REFUSED rather than auto-enabled: an answer-changing knob that
         ! a second namelist group switches on behind the user's back is
         ! exactly the class of silent coupling this file exists to
         ! prevent.  A UNIFORM draft is exempt — a load with no gradient
         ! is bit-identically inert in the top BC (the theorem in
         ! `compute_fv_mom6_impl`'s docstring) — which is what keeps the
         ! flat-lid datum-equivalence gate expressible.  `draft_config`
         ! "none" is uniform (identically zero); "flat" is uniform only
         ! when no box bound clips it, else the calving front is a step.
         if (.not. cfg%ocean%pgf%p_top_in_bc) then
            block
               character(len=:), allocatable :: dcfg
               dcfg = trim(adjustl(cfg%ocean%cavity_dyn%draft_config))
               if (.not. cavity_draft_is_uniform(cfg)) then
                  call logger%error("&ocean_cavity_dyn_nml enable=.true. with "// &
                                    "draft_config='"//dcfg//"' requires "// &
                                    "&ocean_pgf_nml p_top_in_bc=.true.  That knob is "// &
                                    "the ONLY route by which the isostatic load "// &
                                    "rho_ref*g*z_draft reaches the FV_MOM6 pa(nz+1) "// &
                                    "surface boundary condition; without it a draft "// &
                                    "that VARIES leaves the pressure stack ~5e6 Pa "// &
                                    "off its anomaly scale and the column out of "// &
                                    "hydrostatic balance.  Only a draft that is "// &
                                    "uniform over the whole domain is exempt (a load "// &
                                    "with no gradient is provably inert there).")
                  has_error = .true.
               end if
            end block
         end if
         if (cfg%ocean%pgf%gfs_scale /= 1.0_wp) then
            call logger%error("&ocean_cavity_dyn_nml enable=.true. requires "// &
                              "&ocean_pgf_nml gfs_scale=1: the datum "// &
                              "(bt_H_ref, which the barotropic substep feels "// &
                              "through g_bt = gfs_scale*GRAVITY) and the load "// &
                              "(rho_ref*GRAVITY*z_draft, which the PGF feels "// &
                              "through GRAVITY) would then sit on two different "// &
                              "gravities and drift apart.")
            has_error = .true.
         end if
         ! --- vertical-coordinate envelope ---
         ! sigma (and zstar-lite, which shares its ocean branch) rescales
         ! the live column and so follows the draft for free; `z_fixed`
         ! is the FIRST family taught about the ice base explicitly
         ! (P6.2) — it reads `vcoord%z_top`, keeps its nominal interface
         ! depths GEOPOTENTIAL, vanishes the layers that outcrop into the
         ! ice to the inert filler and cuts the first live layer at the
         ! draft (Yung, Hallberg, Adcroft & Morrison 2026, JAMES 18,
         ! e2025MS005645, Fig. 1b).  Its own envelope is fenced below.
         !
         ! The refusal used to be a two-value whitelist whose message
         ! enumerated SIX families and named neither `lagrangian` nor
         ! `zstar_sigma` — both of which it refused.  A refusal that does
         ! not name what it refused, or gives a reason that is not the
         ! real one, sends the operator to fix the wrong thing.  So the
         ! accept test stays a whitelist (nothing else has been validated
         ! under a shelf) but the message now carries the offending
         ! family's OWN reason, and the reasons are not all "anchors at
         ! z = 0": three of the seven refused families are geometrically
         ! datum-safe and are refused for want of validation, which is a
         ! different, and recoverable, kind of no.  Measured per family in
         ! `test_ocean_vcoord_interface_depths`.
         block
            integer :: cav_vcoord_code
            real(wp) :: cav_h_nominal
            character(len=:), allocatable :: cav_reason
            cav_vcoord_code = parse_vcoord_type(cfg%vcoord_type, &
                                                default_code=VCOORD_EULERIAN_Z)
            if (.not. (cav_vcoord_code == VCOORD_SIGMA .or. &
                       cav_vcoord_code == VCOORD_Z_FIXED)) then
               select case (cav_vcoord_code)
               case (VCOORD_ZSTAR)
                  cav_reason = "'zstar' is MOM6 z*: the fixed z_fixed nominal "// &
                               "profile dilated by (H + eta)/H from the FREE "// &
                               "SURFACE.  It has no rigid-top branch (MOM6's "// &
                               "build_zstar_column z_rigid_top path is not "// &
                               "ported), so under a draft its fine near-surface "// &
                               "levels would hang from the ice base; it was "// &
                               "'sigma' under another name until the z* slice, "// &
                               "and is refused under a cavity until that branch "// &
                               "lands.  Use 'z_fixed' (the z-like family taught "// &
                               "the ice base) or 'sigma'"
               case (VCOORD_LAGRANGIAN)
                  cav_reason = "'lagrangian' is geometrically datum-FREE (the target "// &
                               "IS the live h_layer and the remap is a no-op), so the "// &
                               "placement objection does not apply to it.  It is "// &
                               "refused in v1 only because no cavity run has been "// &
                               "validated on an isopycnal coordinate; it is the "// &
                               "intended isopycnal control leg of the coordinate study"
               case (VCOORD_EULERIAN_Z)
                  cav_reason = "'eulerian_z' is a stretched SIGMA with the free "// &
                               "surface DROPPED (it is not a geopotential coordinate "// &
                               "despite the name), and the dropped eta — not a z "// &
                               "anchor — is why it cannot carry an ice-base "// &
                               "displacement"
               case (VCOORD_ZSIGMA)
                  cav_reason = "'zsigma' is refused on EVERY path, cavity or not: its "// &
                               "deep branch reads z_ref_global as metres while the "// &
                               "only writer fills it with the dimensionless k/nz"
               case (VCOORD_ZSTAR_SIGMA)
                  cav_reason = "'zstar_sigma' is a purely FRACTIONAL rescale of the "// &
                               "global reference table, hence datum-invariant and "// &
                               "geometrically safe under a draft.  It is refused in "// &
                               "v1 only because it follows BOTH boundaries (so it "// &
                               "solves nothing a cavity needs) and no cavity run has "// &
                               "been validated on it"
               case (VCOORD_ZSTAR_FULL)
                  cav_reason = "'zstar_full' builds its per-column reference table "// &
                               "from the TRUE bed while the target walk sees only the "// &
                               "live thickness, so under a draft it keeps the table's "// &
                               "SHALLOW entries: the fine near-surface band lands "// &
                               "against the ice base and the deep water is carried by "// &
                               "shallow-ocean spacing.  It does not merely anchor at "// &
                               "z = 0, it inverts which half of the column is resolved"
               case default
                  cav_reason = "'"//trim(cfg%vcoord_type)//"' is a DENSITY-space "// &
                               "coordinate with no geometric anchor, so the placement "// &
                               "objection does not apply.  It is refused in v1 because "// &
                               "cavity x isopycnal is unvalidated; note also that "// &
                               "hycom's z* nominal-floor band accumulates from the "// &
                               "column top, so under a shelf that band is "// &
                               "draft-FOLLOWING and must not be quoted as z-like"
               end select
               call logger%error("&ocean_cavity_dyn_nml enable=.true. accepts "// &
                                 "vcoord_type='sigma' or 'z_fixed' ONLY — the family "// &
                                 "that rescales the live column and so follows the "// &
                                 "ice base for free, plus the one that has been "// &
                                 "TAUGHT the ice base "// &
                                 "(z_fixed reads vcoord%z_top, keeps its nominal "// &
                                 "interface depths geopotential, and vanishes the "// &
                                 "layers that outcrop into the ice).  Got '"// &
                                 trim(cfg%vcoord_type)//"': "//cav_reason//".")
               has_error = .true.
            end if

            if (cav_vcoord_code == VCOORD_Z_FIXED) then
               ! ---- The staircase, and why this is a WARNING ----
               !
               ! A z-like coordinate removes the interior sigma tilt but
               ! replaces it with the ice-base STAIRCASE: where the draft
               ! crosses a nominal level the two columns' filler counts
               ! differ by one, the interface offset across that face
               ! jumps by up to `h_nominal`, and the FV_MOM6 acceleration
               ! in a vanished layer is h-INDEPENDENT.  The corrections
               ! that arrest it — top-side mass weighting (MWIPG) and the
               ! interior reference interface, Yung, Hallberg, Adcroft &
               ! Morrison (2026), JAMES 18, e2025MS005645, §3.3.2 and
               ! §3.2 — are NOT implemented.  A UNIFORM draft has no
               ! staircase (every column vanishes the same layers and
               ! cuts at the same depth, so every interface offset is
               ! identically zero) and is validated bit-zero to 30 days;
               ! a VARYING draft is not validated at all, and measured on
               ! `cavity_sloping_lid_rest_zfixed.nml` it is 47x the sigma
               ! leg's resting pressure-gradient residual at step 1 and
               ! ends in a non-finite state on day 18 (gfortran) or a
               ! saturated En = 3.8E-04 (nvfortran GPU).  That is the
               ! knob-OFF state.  `&vcoord_nml zfixed_closed_faces`
               ! (Adcroft, Hill & Marshall 1997 partial steps) closes the
               ! staircase faces, and with it the same file completes 30
               ! days and ISOMIP+ Ocean0 runs 30 days — but the residual on
               ! the faces left open is still ~2 decades above the sigma
               ! leg, so the warning stays, reworded for each state.
               !
               ! WARNING, not a refusal, deliberately: those corrections
               ! are the next slices and they need this configuration to
               ! be runnable to be developed and measured against.  What
               ! the user must not do is walk into it silently.
               if (.not. cavity_draft_is_uniform(cfg)) then
                  if (cfg%zfixed_closed_faces) then
                     call logger%warning("&vcoord_nml vcoord_type='z_fixed' under "// &
                                         "&ocean_cavity_dyn_nml with a draft that VARIES "// &
                                         "(draft_config='"// &
                                         trim(adjustl(cfg%ocean%cavity_dyn%draft_config))// &
                                         "') is EXPERIMENTAL.  &vcoord_nml "// &
                                         "zfixed_closed_faces=.true. closes the ice-base "// &
                                         "STAIRCASE faces, which is what makes this "// &
                                         "configuration runnable: the sloping-lid rest "// &
                                         "case completes 30 days and ISOMIP+ Ocean0 runs "// &
                                         "30 days.  It does not remove the staircase "// &
                                         "residual on the faces that stay open — the "// &
                                         "top-side mass weighting (MWIPG) and the interior "// &
                                         "reference interface, Yung, Hallberg, Adcroft & "// &
                                         "Morrison (2026), JAMES 18, e2025MS005645, "// &
                                         "sections 3.3.2 and 3.2, are not implemented — so "// &
                                         "the resting sloping-lid case still carries about "// &
                                         "two decades more spurious energy than its sigma "// &
                                         "leg.  Only a UNIFORM draft is bit-zero at rest.")
                  else
                     call logger%warning("&vcoord_nml vcoord_type='z_fixed' under "// &
                                         "&ocean_cavity_dyn_nml with a draft that VARIES "// &
                                         "(draft_config='"// &
                                         trim(adjustl(cfg%ocean%cavity_dyn%draft_config))// &
                                         "') and &vcoord_nml zfixed_closed_faces=.false. "// &
                                         "is NOT VALIDATED.  Only a UNIFORM draft is: "// &
                                         "there every column vanishes the same layers and "// &
                                         "cuts at the same depth, so the answer is "// &
                                         "bit-zero at rest.  Where the draft crosses a "// &
                                         "nominal level the filler count changes column "// &
                                         "to column and the FV pressure gradient across "// &
                                         "the open ice-base STAIRCASE drives a spurious "// &
                                         "flow: measured at rest, 47x the sigma leg's "// &
                                         "step-1 residual, and the run does not survive "// &
                                         "30 days.  Set zfixed_closed_faces=.true. (the "// &
                                         "partial-step face closure, with which a varying "// &
                                         "draft does survive 30 days), or use "// &
                                         "vcoord_type='sigma' for a sloping lid.")
                  end if
               end if
               ! ---- z_fixed × cavity, v1 envelope ----
               !
               ! Under a quasi-geopotential coordinate the ice base cuts
               ! the nominal stack, so on an ICE-COVERED column `k = nz`
               ! is an inert filler (`h <= H_VANISHED`), NOT the
               ! ice-adjacent live layer.  This is a state no consumer in
               ! the tree has ever seen: no family on this branch
               ! vanishes a layer against the TOP.  The shared
               ! `k_top(i,j)` (the first live layer, counting down) that
               ! fixes them is the NEXT slice — P6.3 for the tracer/flux
               ! consumers, P6.4 for momentum and the boundary-layer
               ! schemes — so until it lands the only thing this
               ! combination may run is ADIABATIC DYNAMICS.
               !
               ! What is NOT refused, and why: OPEN-OCEAN columns have
               ! `z_top = 0`, so `k = nz` is live there and behaves
               ! exactly as today; and on a COVERED column the binary
               ! cover mask (`cover_frac`) already zeroes every
               ! ATMOSPHERIC input at source — wind stress (and with it
               ! `stress_mag`, hence u*), the scalar and assembled
               ! surface heat/salt fluxes, the shortwave deposit and both
               ! restoring increments.  Those therefore compose with a
               ! vanished `k = nz` and are left alone.  What follows is
               ! everything that acts ON the ice-adjacent layer itself,
               ! and so is not masked by anything.
               ! NOTE `&ocean_cavity_melt_nml enable` and
               ! `&ocean_tdrag_nml enable` used to be refused HERE, and
               ! are not any more: both are routed through the shared
               ! first-live-layer index `ms%k_top` (and its two face
               ! twins) — P6.3/P6.4.  The melt heat/salt deposit and the
               ! `freshwater="mass"` column source land on
               ! `k_top(i,j)`; the ice-ocean drag's band walk starts at
               ! `k_top_u/v`, gates on `H_VANISHED` instead of on zero,
               ! and captures its implicit-fold rate on the same row the
               ! vdiff diagonal adds it to.  `k_top ≡ nz` off a rigid
               ! top, so nothing else moved.
               if (cfg%ocean%vmix%use_kpp) then
                  call logger%error("&ocean_vmix_nml use_kpp=.true. is refused with "// &
                                    "vcoord_type='z_fixed' under a cavity: KPP's "// &
                                    "surface reference column is k = nz "// &
                                    "(b_ref = -g*rho_layer(nz)/rho_0, "// &
                                    "d_centre_ref = h(nz)/2), which on a covered "// &
                                    "column is the filler — a spurious buoyancy jump "// &
                                    "at the very first interface and a reference "// &
                                    "depth of ~1e-4 m.  Set use_kpp=.false. (an "// &
                                    "adiabatic cavity run is the v1 envelope).  "// &
                                    "Needs the shared k_top: follow-up P6.4.")
                  has_error = .true.
               end if
               if (cfg%ocean%epbl%enable) then
                  call logger%error("&ocean_epbl_nml enable=.true. is refused with "// &
                                    "vcoord_type='z_fixed' under a cavity: the EPBL "// &
                                    "column captures the surface at k = nz and gates "// &
                                    "on `h > 0` rather than on H_VANISHED, so a "// &
                                    "filler top layer reaches the specific-volume "// &
                                    "derivatives as an absurd concentration.  Needs "// &
                                    "the shared k_top: follow-up P6.4.")
                  has_error = .true.
               end if
               if (cfg%ocean%tracers%enable_ideal_age) then
                  call logger%error("&ocean_tracers_nml enable_ideal_age=.true. is "// &
                                    "refused with vcoord_type='z_fixed' under a "// &
                                    "cavity: the young-band reset writes "// &
                                    "hTr_age(:,:,nz) = young*h(:,:,nz), which on a "// &
                                    "covered column ventilates a filler (and there "// &
                                    "is no ventilation under a shelf at all).  "// &
                                    "Follow-up P6.3.")
                  has_error = .true.
               end if
               if (cfg%ocean%gm%enable) then
                  call logger%error("&ocean_gm_nml enable=.true. is refused with "// &
                                    "vcoord_type='z_fixed' under a cavity: without "// &
                                    "zfixed_closed_faces the non-divergence closure "// &
                                    "dumps the residual streamfunction transport "// &
                                    "into k = nz (a filler under the ice); with it "// &
                                    "the closure is open-column, but the slopes "// &
                                    "slot GM needs is itself refused under a "// &
                                    "cavity (below).  GM x cavity is unvalidated "// &
                                    "on any coordinate; revisited with the "// &
                                    "coordinate study.")
                  has_error = .true.
               end if
               if (cfg%ocean%redi%enable .or. cfg%ocean%slopes%enable) then
                  call logger%error("&ocean_redi_nml / &ocean_slopes_nml enable="// &
                                    ".true. is refused with vcoord_type='z_fixed' "// &
                                    "under a cavity: the isopycnal-slope surface "// &
                                    "fill is built from h(:,:,nz), the filler on a "// &
                                    "covered column.  Unvalidated with a cavity on "// &
                                    "any coordinate; revisited with the coordinate "// &
                                    "study.")
                  has_error = .true.
               end if
               ! ---- kappa_h: refused ONLY with the face mask off ----
               !
               ! The objection is unchanged where it still applies: the
               ! along-coordinate tracer diffusion gates its T = hTr/h
               ! division on `h > 0` (1/0 armour) and not on H_VANISHED,
               ! weights the face flux by the ARITHMETIC mean thickness
               ! — a filler beside a 20 m cell is weighted by 10 m — and
               ! carries no mass-availability limiter, so it can drive a
               ! vanished cell's hTr strongly negative in one step.
               !
               ! But `&vcoord_nml zfixed_closed_faces` already answers
               ! it, from the other end: `ocean_hdiff_tracer_step`
               ! multiplies every face flux by `metrics%open_u/open_v`,
               ! which is exactly zero wherever the layer is a filler on
               ! EITHER side.  A filler cell then has all four of its
               ! own-layer faces closed, so its hTr divergence is
               ! identically zero and the garbage concentration the
               ! `h > 0` gate computes for it never leaves the cell.
               ! Flux-zero, not flux-limited, so the scheme stays
               ! conservative by construction.  That is the P6.5 gate,
               ! reached through the mask instead of through a new
               ! threshold — and it is a strictly stronger statement,
               ! because it also closes the partial⇄filler face the
               ! threshold alone would leave open on the thick side.
               if (cfg%ocean%hdiff%kappa_h /= 0.0_wp .and. &
                   .not. cfg%zfixed_closed_faces) then
                  call logger%error("&ocean_hdiff_nml kappa_h /= 0 is refused with "// &
                                    "vcoord_type='z_fixed' under a cavity UNLESS "// &
                                    "&vcoord_nml zfixed_closed_faces=.true.: the "// &
                                    "along-coordinate tracer diffusion gates its "// &
                                    "T = hTr/h division on `h > 0` (1/0 armour), "// &
                                    "not on H_VANISHED, weights the face flux by the "// &
                                    "ARITHMETIC mean thickness — so a filler beside "// &
                                    "a 20 m cell is weighted by 10 m — and has no "// &
                                    "mass-availability limiter, so it can drive a "// &
                                    "vanished cell's hTr strongly negative in one "// &
                                    "step.  With the partial-step face mask on, "// &
                                    "every face touching a filler is CLOSED and the "// &
                                    "flux is exactly zero, which answers all three.  "// &
                                    "Set zfixed_closed_faces=.true., or kappa_h=0.")
                  has_error = .true.
               end if
               if (cfg%ocean%kshear%enable) then
                  call logger%error("&ocean_kappa_shear_nml enable=.true. is refused "// &
                                    "with vcoord_type='z_fixed' under a cavity: the "// &
                                    "JHL08 column solve closes its SURFACE row on "// &
                                    "k = nz (u_c/v_c/t_c/s_c(nz)), which on a "// &
                                    "covered column is an inert filler inside the "// &
                                    "ice draft, and its own massless-merge helper is "// &
                                    "off by default.  Not covered by the shared "// &
                                    "k_top (P6.3/P6.4), which routes the FORCING "// &
                                    "sites; a column solver wants the compacted "// &
                                    "column rdb_massless already builds.")
                  has_error = .true.
               end if
               if (cfg%ocean%tidal_mixing%enable) then
                  call logger%error("&ocean_tidal_mixing_nml enable=.true. is "// &
                                    "refused with vcoord_type='z_fixed' under a "// &
                                    "cavity: the N^2 column sets its top boundary "// &
                                    "at k = nz and forms the k = nz-1 interface "// &
                                    "spacing as 0.5*(h(nz-1) + h(nz)), which on a "// &
                                    "covered column HALVES that spacing against a "// &
                                    "filler and inflates the buoyancy frequency at "// &
                                    "the first live interface.  Unvalidated under a "// &
                                    "rigid top.")
                  has_error = .true.
               end if
               if (cfg%regrid_time_scale > 0.0_wp) then
                  call logger%error("&vcoord_nml regrid_time_scale > 0 is refused "// &
                                    "with vcoord_type='z_fixed' under a cavity: the "// &
                                    "grid time-filter is a per-layer convex blend "// &
                                    "that does not respect the vanish marker, so a "// &
                                    "filler relaxing toward h_min from a live "// &
                                    "thickness passes THROUGH H_VANISHED and the "// &
                                    "layer oscillates live/dead on successive steps "// &
                                    "— the remap drain deleting its content on the "// &
                                    "step it reads dead.  Excluding the fillers from "// &
                                    "the blend is its own change.")
                  has_error = .true.
               end if
               ! ---- In-layer T/S reconstruction for the PGF: REFUSED ----
               !
               ! Measured on ISOMIP+ Ocean0 idealised (z_fixed + closed
               ! faces, melt off), V100, on the v0.1.0 defaults (bebt =
               ! 0.1, renorm_consistent_flux, I1'): from rest the knob
               ! spins up En = 3.4E-04 m2/s2 in the first 3 hours
               ! (MaxCFL 0.43) and holds ~1E-03 (MaxCFL up to 0.85) for
               ! 30 days — ~19 000x the layer-mean-density run's 5.6E-08
               ! and ~5x the melt-ON circulation.  (Before I1' and bebt =
               ! 0.1 it went non-finite at outer step 118, day 0.41.)  The
               ! PLM/PPM edge build reads the inert top-side FILLERS as
               ! neighbouring water at the partial top cell.  The
               ! filler-aware reconstruction that skips them is a held
               ! slice, so until it lands this is a refusal, not a
               ! warning: the spurious flow is there within 3 hours.
               ! (`cavity_rest_growth_diagnosis.md` §Q.0 item 9.)
               if (cfg%ocean%pgf%reconstruct_for_pressure) then
                  call logger%error("&ocean_pgf_nml reconstruct_for_pressure=.true. "// &
                                    "is refused with vcoord_type='z_fixed' under a "// &
                                    "cavity: the in-layer PLM/PPM T/S edge build "// &
                                    "reads the layer means of the inert top-side "// &
                                    "FILLERS as neighbouring water at the partial top "// &
                                    "cell — ISOMIP+ Ocean0 (melt off, from rest) spins "// &
                                    "up a spurious En ~ 1E-03 m2/s2 within hours, "// &
                                    "19 000x the layer-mean-density run.  Pending the "// &
                                    "filler-aware reconstruction that skips vanished "// &
                                    "layers; set reconstruct_for_pressure=.false. "// &
                                    "(the layer-mean PCM density) until it lands.")
                  has_error = .true.
               end if
               ! ---- Lateral viscosity lower envelope: a WARNING ----
               !
               ! Not a refusal: the vcoord stability matrix runs this
               ! combination inviscid on purpose, to measure the mode.
               ! But a user who sets nu_h below the envelope should be
               ! told.  See `zfixed_cavity_nu_h_below_envelope`.
               if (zfixed_cavity_nu_h_below_envelope(cfg)) then
                  call logger%warning("&ocean_hvisc_nml nu_h = "// &
                                      to_string(cfg%ocean%hvisc%nu_h)// &
                                      " m2/s is below the "// &
                                      to_string(ZFIXED_CAVITY_NU_H_MIN)// &
                                      " m2/s lower envelope of vcoord_type='z_fixed' "// &
                                      "under a cavity.  Measured on ISOMIP+ Ocean0 "// &
                                      "(2 km, melt off): nu_h = 0 carries an INVISCID "// &
                                      "mode that grows EXPONENTIALLY at 0.18 /day "// &
                                      "(5.6-day e-folding, accelerating; En 1.4E-06 at "// &
                                      "day 30, 26x the protocol leg) while nu_h = 2 "// &
                                      "already decelerates.  Flow-aware closures do "// &
                                      "not replace it at these speeds (Smagorinsky "// &
                                      "adds ~1 m2/s).  The ISOMIP+ Table-4 value is "// &
                                      "6; nu_h >= 30 also removes the calving-front "// &
                                      "partial-cell jet (cavity_rest_growth_diagnosis "// &
                                      "section Q).  Proceeding, as configured.")
               end if
               ! ISOMIP+ (Asay-Davis et al. 2016 §3.1.5): "the minimum
               ! thickness is likely to be approximately two grid cells
               ! (~40 m if z levels are equally spaced)".  Under a
               ! terrain-following coordinate a just-afloat column still
               ! carries all nz layers; under a z-like one a column
               ! thinner than 2*h_nominal carries ONE partial live layer
               ! and nz-1 fillers, and every column-walking consumer then
               ! operates on a single cell.  `h_nominal` is
               ! `max_depth/nz` — Z_FIXED's spacing is written from
               ! `&ocean_topo_nml max_depth`, deliberately with no second
               ! spelling.
               cav_h_nominal = 0.0_wp
               if (cfg%nz_layers > 0) then
                  cav_h_nominal = cfg%ocean%topo%max_depth/real(cfg%nz_layers, wp)
               end if
               ! A stretched profile has no single `h_nominal`: whether a
               ! thin cavity column spans two layers depends on the DEPTH
               ! of its draft, which configure does not see.  The check
               ! is then made against the thickest nominal layer and
               ! downgraded to a warning — conservative, and it never
               ! refuses a configuration the geometry might satisfy.
               if (parse_z_fixed_profile(cfg%z_fixed_profile) /= ZFIXED_PROFILE_UNIFORM .and. &
                   parse_z_fixed_profile(cfg%z_fixed_profile) /= ZFIXED_PROFILE_INVALID .and. &
                   cfg%nz_layers > 0) then
                  block
                     real(wp), allocatable :: cav_dz(:)
                     integer :: cav_ierr
                     allocate (cav_dz(cfg%nz_layers))
                     call z_fixed_nominal_dz(parse_z_fixed_profile(cfg%z_fixed_profile), &
                                             cfg%nz_layers, cfg%ocean%topo%max_depth, &
                                             cfg%z_fixed_dz, cfg%z_fixed_dz_top, &
                                             cfg%z_fixed_tanh_center, &
                                             cfg%z_fixed_tanh_width, cav_dz, cav_ierr)
                     if (cav_ierr == ZFIXED_DZ_OK .and. &
                         cfg%ocean%cavity_dyn%h_min_cavity < 2.0_wp*maxval(cav_dz)) then
                        call logger%warning("&ocean_cavity_dyn_nml h_min_cavity = "// &
                                            to_string(cfg%ocean%cavity_dyn%h_min_cavity)// &
                                            " m is below twice the thickest nominal "// &
                                            "z_fixed layer ("//to_string(maxval(cav_dz))// &
                                            " m) of the stretched profile: a cavity column "// &
                                            "thinner than two nominal layers AT ITS DEPTH "// &
                                            "carries a single partial live layer (ISOMIP+ "// &
                                            "Asay-Davis et al. 2016, §3.1.5).  Not refused: "// &
                                            "under a stretched profile it depends on the draft.")
                     end if
                  end block
                  cav_h_nominal = 0.0_wp
               end if
               if (cav_h_nominal > 0.0_wp .and. &
                   cfg%ocean%cavity_dyn%h_min_cavity < 2.0_wp*cav_h_nominal) then
                  call logger%error("&ocean_cavity_dyn_nml h_min_cavity = "// &
                                    to_string(cfg%ocean%cavity_dyn%h_min_cavity)// &
                                    " m is below 2*h_nominal = "// &
                                    to_string(2.0_wp*cav_h_nominal)//" m "// &
                                    "(h_nominal = &ocean_topo_nml max_depth / "// &
                                    "nz_layers = "//to_string(cav_h_nominal)// &
                                    " m), which vcoord_type='z_fixed' requires under "// &
                                    "a cavity: a thinner water column carries a "// &
                                    "single partial live layer and nz-1 inert "// &
                                    "fillers.  This is ISOMIP+'s own rule "// &
                                    "(Asay-Davis et al. 2016, GMD 9, 2471, "// &
                                    "§3.1.5).  NOTE it moves the grounding line "// &
                                    "relative to a sigma leg — a confound the "// &
                                    "coordinate study must control for.")
                  has_error = .true.
               end if
            end if
         end block
         if (trim(cfg%thickness_config) == "uniform_z") then
            call logger%error("&ocean_cavity_dyn_nml enable=.true. is mutually "// &
                              "exclusive with thickness_config='uniform_z': that "// &
                              "seed lays uniform z interfaces from z = 0 down, so "// &
                              "under a draft it would fill the ice with water.  Use "// &
                              "the default sigma-style split.")
            has_error = .true.
         end if
         ! NOTE: `&ocean_zinit_nml enable` used to be refused here — the
         ! overlay measured depth from the COLUMN TOP, which under a
         ! draft is `z_draft` metres below `z = 0`, so a geopotential
         ! profile landed systematically too shallow.  `build_z_ctr` now
         ! takes the column-top depth and the seed passes
         ! `metrics%z_draft`, so the two compose and the refusal is gone.
         ! --- solver envelope ---
         if (cfg%ocean%bt%n_inner < 1 .and. .not. cfg%ocean%bt%auto_n_inner) then
            call logger%error("&ocean_cavity_dyn_nml enable=.true. requires the "// &
                              "SPLIT solver (&ocean_bt_nml n_inner >= 1, or "// &
                              "auto_n_inner): the unsplit driver has neither the "// &
                              "barotropic correction nor the eta_forcing seam, so it "// &
                              "carries no barotropic response to a surface load and "// &
                              "the cavity is unvalidated there.")
            has_error = .true.
         end if
         if (cfg%ocean%bt%bt_halo > 0) then
            call logger%error("&ocean_cavity_dyn_nml enable=.true. is mutually "// &
                              "exclusive with &ocean_bt_nml bt_halo > 0: the "// &
                              "wide-halo BT clone rebuilds its own metrics from the "// &
                              "grid formula and carries no z_draft, so its reference "// &
                              "depth would be the bed and its solve would ignore the "// &
                              "ice.")
            has_error = .true.
         end if
         if (cfg%px*cfg%py > 1) then
            call logger%error("&ocean_cavity_dyn_nml enable=.true. is single-rank in "// &
                              "v1 (px*py = "//to_string(cfg%px*cfg%py)//").  The "// &
                              "draft halo itself is two lines, but the grounding "// &
                              "statistics are global reductions the v1 configure does "// &
                              "not take, and draft_config='file' inherits the "// &
                              "local-nx reader.  Lifting the fence is its own PR.")
            has_error = .true.
         end if
         ! --- mutually exclusive capabilities ---
         if (cfg%ocean%wetdry%enable) then
            call logger%error("&ocean_cavity_dyn_nml enable=.true. is mutually "// &
                              "exclusive with &ocean_wetdry_nml enable: both decide "// &
                              "whether a column can carry water, from different "// &
                              "thresholds, and wet/dry also forces split_scheme="// &
                              "'ssp_rk2'.")
            has_error = .true.
         end if
         if (cfg%ocean%porous%enable) then
            call logger%error("&ocean_cavity_dyn_nml enable=.true. is mutually "// &
                              "exclusive with &ocean_porous_nml enable: both narrow "// &
                              "the same faces from static geometry and the "// &
                              "combination is unvalidated.")
            has_error = .true.
         end if
         if (cfg%ocean%ice%enable) then
            call logger%error("&ocean_cavity_dyn_nml enable=.true. is mutually "// &
                              "exclusive with &ocean_ice_nml enable: sea ice is a "// &
                              "SECOND surface load on the same column, and the two "// &
                              "are not reconciled in v1.")
            has_error = .true.
         end if
         if (cfg%ocean%tides%enable .and. cfg%ocean%tides%use_sal) then
            call logger%error("&ocean_cavity_dyn_nml enable=.true. is mutually "// &
                              "exclusive with &ocean_tides_nml use_sal: the scalar "// &
                              "SAL elevation is beta_sal*bt_eta, and under the cavity "// &
                              "datum bt_eta is the departure from the LOADED "// &
                              "equilibrium, not the sea-surface elevation the SAL "// &
                              "response is defined on.  The body tide itself "// &
                              "(enable alone) is datum-independent and allowed.")
            has_error = .true.
         end if
         ! --- honest-inertness warnings (the house rule: warn, do not refuse) ---
         if (trim(adjustl(cfg%ocean%cavity_dyn%draft_config)) == "none") then
            call logger%warning("&ocean_cavity_dyn_nml enable=.true. with "// &
                                "draft_config='none': z_draft is identically zero, so "// &
                                "bt_H_ref = b and the run is bit-identical to a "// &
                                "cavity-free one.")
         end if
         ! Design Q4: WARN, do not refuse, on cavity x in_eos=.false.  The
         ! load is depth-uniform in `pa`, so the rest and equivalence
         ! gates pass either way; what is wrong without `in_eos` is the
         ! THERMOBARICITY, and only a pressure-dependent EOS has any.  A
         ! linear EOS is pressure-blind, so there the knob is honestly
         ! inert and there is nothing to warn about.
         if (.not. cfg%ocean%psurf%in_eos .and. &
             trim(adjustl(cfg%ocean%eos%eos)) /= "linear") then
            call logger%warning("&ocean_cavity_dyn_nml enable=.true. with a "// &
                                "NONLINEAR equation of state (&ocean_eos_nml eos='"// &
                                trim(adjustl(cfg%ocean%eos%eos))//"') but without "// &
                                "&ocean_psurf_nml in_eos=.true.: the in-situ EOS "// &
                                "pressure still starts at 0 Pa at the ice base, so "// &
                                "the EOS ignores up to ~5e6 Pa of ice load (wrong "// &
                                "thermobaricity, wrong freezing point).  Harmless "// &
                                "for the rest and equivalence gates; not for a "// &
                                "production cavity.")
         end if
      end if
      ! ---- Ice-shelf basal melt, v1 scope (Phase 2b) ----
      ! Every restriction fails loud.  One of them exists because a piece
      ! of the coupling is NOT in this slice: there is no per-cell cover
      ! mask on the atmospheric forcing yet, and leaving that silently
      ! half-wired would run an atmosphere through several hundred metres
      ! of solid ice.  The interface pressure is NOT such a gap — the
      ! cavity assembles `ms%p_top = p_ice_ref + sf%p_surf` (P5.2) and
      ! the melt liquidus is its third consumer.
      if (cfg%ocean%cavity_melt%enable) then
         if (.not. cfg%ocean%cavity_dyn%enable) then
            call logger%error("&ocean_cavity_melt_nml enable=.true. requires "// &
                              "&ocean_cavity_dyn_nml enable=.true.  The melt "// &
                              "interface stands on that group's geometry: without a "// &
                              "draft there is no ice base, cover_frac is identically "// &
                              "zero and the kernel would never be called.")
            has_error = .true.
         end if
         if (trim(adjustl(cfg%ocean%eos%tfreeze_set)) /= "isomip") then
            call logger%error("&ocean_cavity_melt_nml enable=.true. requires "// &
                              "&ocean_eos_nml tfreeze_set='isomip' (got '"// &
                              trim(adjustl(cfg%ocean%eos%tfreeze_set))//"').  The "// &
                              "two shipped liquidi differ by ~0.03 degC at S = 34.5, "// &
                              "which is enough to flip the SIGN of the basal melt "// &
                              "rate over a 0.03 degC band of far-field temperature — "// &
                              "the sea-ice set is not a defensible default for a "// &
                              "cavity.")
            has_error = .true.
         end if
         if (.not. cfg%ocean%forcing%enable_components) then
            call logger%error("&ocean_cavity_melt_nml enable=.true. requires "// &
                              "&ocean_forcing_nml enable_components=.true.  The melt "// &
                              "fluxes are delivered as the OWNED components "// &
                              "heat_cavity/salt_cavity, which are allocated only "// &
                              "with the component set, and only "// &
                              "ocean_surface_flux_assemble folds them into "// &
                              "Q_heat/Q_salt.  Writing Q_heat/Q_salt directly is "// &
                              "refused by the fill contract — the sea-ice coupler "// &
                              "full-overwrites those.")
            has_error = .true.
         end if
         block
            integer :: melt_law_code, melt_ice_code
            melt_law_code = parse_cavity_exchange_law(cfg%ocean%cavity_melt%exchange_law)
            select case (melt_law_code)
            case (CAVITY_LAW_CONST_GAMMA, CAVITY_LAW_HJ99, CAVITY_LAW_YUNG25)
               continue
            case (CAVITY_LAW_INVALID)
               call logger%error("&ocean_cavity_melt_nml exchange_law = '"// &
                                 trim(adjustl(cfg%ocean%cavity_melt%exchange_law))// &
                                 "' is not recognised (shipped: const_gamma, hj99, "// &
                                 "yung25).")
               has_error = .true.
            case default
               call logger%error("&ocean_cavity_melt_nml exchange_law = '"// &
                                 trim(adjustl(cfg%ocean%cavity_melt%exchange_law))// &
                                 "' is RESERVED, not implemented "// &
                                 "(CAVITY_MELT_NOT_IMPLEMENTED).  Its enum value is "// &
                                 "nailed down so adding it later is not a "// &
                                 "renumbering, and the Python prototype has it, but "// &
                                 "no Fortran physics ships.  Shipped: const_gamma, "// &
                                 "hj99, yung25.")
               has_error = .true.
            end select
            ! Holland & Jenkins (1999) eq. (15) p. 1792 takes
            ! `ln(u* xi_N eta*^2 / (|f| h_nu))` and eq. (18) divides by
            ! `f L_O`: the law does not exist on the equator.  The
            ! per-column check over the covered cells is taken at
            ! configure (`configure_ocean_cavity_melt`), where f_centre
            ! exists; here we can only catch the whole-domain f = 0 case,
            ! and catching it early is worth the duplication.
            if (melt_law_code == CAVITY_LAW_HJ99) then
               if (cfg%coriolis_f == 0.0_wp .and. &
                   cfg%ocean%topo%coriolis_beta == 0.0_wp) then
                  call logger%error("&ocean_cavity_melt_nml exchange_law='hj99' on an "// &
                                    "f = 0 grid (coriolis_f = 0, coriolis_beta = 0). "// &
                                    "Holland & Jenkins (1999) eq. (15) takes "// &
                                    "ln(.../|f| h_nu) and eq. (18) divides by f*L_O, "// &
                                    "so the law has no value there — the same "// &
                                    "fail-loud stance &ocean_vmix_nml bkgnd_henyey "// &
                                    "takes on a cartesian grid.  Use "// &
                                    "exchange_law='const_gamma' or set a Coriolis "// &
                                    "parameter.")
                  has_error = .true.
               end if
            end if
            melt_ice_code = parse_cavity_ice_mode(cfg%ocean%cavity_melt%ice_conduction)
            select case (melt_ice_code)
            case (CAVITY_ICE_INSULATING, CAVITY_ICE_ADV_DIFF)
               continue
            case (CAVITY_ICE_INVALID)
               call logger%error("&ocean_cavity_melt_nml ice_conduction = '"// &
                                 trim(adjustl(cfg%ocean%cavity_melt%ice_conduction))// &
                                 "' is not recognised (shipped: insulating, "// &
                                 "adv_diff).")
               has_error = .true.
            case default
               call logger%error("&ocean_cavity_melt_nml ice_conduction = '"// &
                                 trim(adjustl(cfg%ocean%cavity_melt%ice_conduction))// &
                                 "' is RESERVED, not implemented "// &
                                 "(CAVITY_MELT_NOT_IMPLEMENTED).  The steady "// &
                                 "diffusive form makes q_ice independent of m_mass, "// &
                                 "so sign(m) = sign(T*) no longer holds and the "// &
                                 "pre-solve melt/freeze branch has to be revisited — "// &
                                 "it is not a coefficient change.")
               has_error = .true.
            end select
         end block
         if (cfg%ocean%cavity_melt%gamma_t <= 0.0_wp) then
            call logger%error("&ocean_cavity_melt_nml gamma_t must be > 0 (it "// &
                              "multiplies u* to give the heat exchange velocity; "// &
                              "zero is identically zero melt, which is what "// &
                              "enable=.false. is for).")
            has_error = .true.
         end if
         if (cfg%ocean%cavity_melt%gamma_s >= 0.0_wp .and. &
             cfg%ocean%cavity_melt%gamma_s <= 0.0_wp) then
            call logger%error("&ocean_cavity_melt_nml gamma_s = 0 is refused: the "// &
                              "three-equation form divides by gamma_s.  Leave it "// &
                              "NEGATIVE (the default) to take the ISOMIP+ "// &
                              "gamma_t/35, or set a positive value.")
            has_error = .true.
         end if
         if (cfg%ocean%cavity_melt%far_field_depth <= 0.0_wp) then
            call logger%error("&ocean_cavity_melt_nml far_field_depth must be > 0 m "// &
                              "(it is the thickness the far-field T/S/u are averaged "// &
                              "over; zero would sample nothing).")
            has_error = .true.
         end if
         ! --- Phase 3: real freshwater MASS, and its sea-level partner ---
         ! `freshwater="mass"` moves the meltwater as a REAL Boussinesq
         ! volume on the top layer.  Two envelope holes are refused BY
         ! NAME rather than half-wired, because each would put the mass
         ! source and the machinery that owns `h_layer(:,:,nz)` out of
         ! step with one another:
         !
         !   * dynamic wet/dry re-decides every step which columns carry
         !     water, and its positive-definite outflow limiter is the
         !     other writer of a top-layer thickness source.  Composing
         !     the two needs the limiter to SEE the melt volume;
         !   * a windowed tracer-advection ratio > 1 freezes `hTr` for
         !     `ratio` steps while `h` keeps moving, so the dilution the
         !     mass form relies on would be applied to a tracer load that
         !     is deliberately stale.
         block
            integer :: fw_code, vc_code
            fw_code = parse_cavity_freshwater(cfg%ocean%cavity_melt%freshwater)
            vc_code = parse_cavity_volume_comp(cfg%ocean%cavity_melt%volume_compensation)
            if (fw_code == CAVITY_FW_INVALID) then
               call logger%error("&ocean_cavity_melt_nml freshwater = '"// &
                                 trim(adjustl(cfg%ocean%cavity_melt%freshwater))// &
                                 "' is not recognised (shipped: virtual, mass).")
               has_error = .true.
            end if
            if (vc_code == CAVITY_VC_INVALID) then
               call logger%error("&ocean_cavity_melt_nml volume_compensation = '"// &
                                 trim(adjustl(cfg%ocean%cavity_melt%volume_compensation))// &
                                 "' is not recognised (shipped: none, "// &
                                 "uniform_open_ocean).")
               has_error = .true.
            end if
            if (fw_code == CAVITY_FW_MASS) then
               if (cfg%ocean%wetdry%enable) then
                  call logger%error("&ocean_cavity_melt_nml freshwater='mass' is "// &
                                    "refused with &ocean_wetdry_nml enable=.true.  "// &
                                    "Wet/dry owns the other top-layer thickness "// &
                                    "source (its positive-definite outflow limiter) "// &
                                    "and re-decides per step which columns hold "// &
                                    "water; composing the two needs the limiter to "// &
                                    "see the melt volume.  Follow-up: "// &
                                    "'cavity real freshwater under wet/dry'.")
                  has_error = .true.
               end if
               if (cfg%ocean%vmix%dt_tracer_advect_ratio > 1) then
                  call logger%error("&ocean_cavity_melt_nml freshwater='mass' is "// &
                                    "refused with &ocean_vmix_nml "// &
                                    "dt_tracer_advect_ratio > 1.  The windowed drain "// &
                                    "holds hTr fixed for the window while h keeps "// &
                                    "moving, so the dilution the mass form relies on "// &
                                    "would act on a deliberately stale tracer load.  "// &
                                    "Follow-up: 'cavity real freshwater under the "// &
                                    "windowed tracer-advect drain'.")
                  has_error = .true.
               end if
            end if
            if (vc_code == CAVITY_VC_UNIFORM_OPEN .and. fw_code /= CAVITY_FW_MASS) then
               call logger%error("&ocean_cavity_melt_nml volume_compensation="// &
                                 "'uniform_open_ocean' requires freshwater='mass'.  "// &
                                 "The virtual form adds no volume, so there is "// &
                                 "nothing to compensate and the sink would be a "// &
                                 "pure, unexplained mass loss.")
               has_error = .true.
            end if
         end block
         ! --- the cover mask SHIPS (P2c) ---
         ! Wind stress, surface restoring, shortwave penetration and the
         ! uniform scalar q_heat/q_salt were each refused here while
         ! there was no per-cell ice-COVER MASK on the atmospheric
         ! forcing.  There is one now:
         !
         !   * the wind-stress PAIR is masked face-wise at configure
         !     (`ocean_surface_stress_apply_cover`, called from
         !     `configure_ocean_cavity`), which also silences the
         !     implicit vdiff stress fold and the MLE front sampler —
         !     both read `ss%tau_x` raw — and refreshes `stress_mag`, so
         !     KPP/EPBL `u*` is zero-wind under cover;
         !   * `q_heat`/`q_salt`, the radiative/turbulent bands, the
         !     mass-flux enthalpies and `salt_flux` are masked in
         !     `ocean_surface_flux_assemble` (cover-aware twin), which is
         !     the one place they are still separable from the cavity's
         !     own `heat_cavity`/`salt_cavity` — and which makes the
         !     `Q_heat`/`Q_salt` KPP/EPBL read for `B_0` the MASKED
         !     values;
         !   * shortwave penetration and surface restoring carry the
         !     factor themselves (they do not route through `Q_*`).
         !
         ! `&ocean_psurf_nml enable` was never refused here: the cavity
         ! is the `ms%p_top` producer (P5.2: `p_ice_ref + sf%p_surf`), so
         ! the seam composes with the ice load rather than clobbering it.
         if (cfg%ocean%ice%enable) then
            call logger%error("&ocean_cavity_melt_nml enable=.true. is mutually "// &
                              "exclusive with &ocean_ice_nml enable=.true.  Sea ice "// &
                              "full-overwrites heat_added/salt_flux and runs its own "// &
                              "instant-relaxation basal flux with a DIFFERENT "// &
                              "liquidus set; two interface thermodynamics in one "// &
                              "column is not a configuration.")
            has_error = .true.
         end if
         ! (E4) Constant α/β under a NONLINEAR EOS in a cavity — a WARNING,
         ! deliberately, not a refusal.  ISOMIP+ prescribes the LINEAR EOS
         ! (Asay-Davis et al. 2016 Table 4), and there the constants ARE
         ! that EOS's exact coefficients, so every shipped cavity namelist
         ! is unaffected and must keep running untouched.  Under Wright or
         ! Roquet the pair is a constant stand-in for a coefficient that
         ! collapses toward zero at the freezing point and roughly doubles
         ! by 1000 dbar — exactly the corner a cavity sits in.
         if (trim(adjustl(cfg%ocean%vmix%buoyancy_coeffs)) == "constant" .and. &
             trim(adjustl(cfg%ocean%eos%eos)) /= "linear") then
            call logger%warning("&ocean_cavity_melt_nml enable=.true. with a "// &
                                "NONLINEAR EOS (&ocean_eos_nml eos='"// &
                                trim(adjustl(cfg%ocean%eos%eos))//"') but "// &
                                "&ocean_vmix_nml buoyancy_coeffs='constant'.  The "// &
                                "KPP surface buoyancy flux B_0 and the "// &
                                "double-diffusion density ratio will size the "// &
                                "melt-driven buoyancy with the SCALAR "// &
                                "&ocean_ic_nml alpha_T/beta_S, not with the "// &
                                "derivatives of the density this run actually "// &
                                "integrates.  Near the freezing point thermal "// &
                                "expansion is several times smaller than at 10 degC "// &
                                "and grows strongly with pressure, so the boundary "// &
                                "layer under the shelf can be mis-sized (and in the "// &
                                "cold-fresh corner mis-signed).  Set "// &
                                "buoyancy_coeffs='eos' unless you are reproducing a "// &
                                "constant-coefficient reference.")
         end if
      end if
      ! ---- Ice-shelf TOP drag (&ocean_tdrag_nml, Phase 4a) ----
      ! Default off ⇒ this whole block is skipped and every path is
      ! bit-identical.  Every restriction below fails loud: a top drag
      ! that is silently inert (no cover, no coefficient) is worse than
      ! no top drag, because a cavity circulation would then be
      ! frictionless at the ice and LOOK like it was damped.
      if (cfg%ocean%tdrag%enable) then
         if (.not. cfg%ocean%cavity_dyn%enable) then
            call logger%error("&ocean_tdrag_nml enable=.true. requires "// &
                              "&ocean_cavity_dyn_nml enable=.true.  The top drag "// &
                              "stands on that group's geometry: without a draft "// &
                              "there is no ice base, cover_frac is identically zero, "// &
                              "and every face mask would be zero — an inert kernel "// &
                              "with a cost.  An OPEN surface's momentum boundary "// &
                              "condition is the wind stress (&ocean_topo_nml "// &
                              "wind_config), not a drag.")
            has_error = .true.
         end if
         block
            integer :: tdrag_code
            tdrag_code = parse_tdrag_variant(cfg%ocean%tdrag%form)
            if (.not. tdrag_variant_is_implemented(tdrag_code)) then
               call logger%error("&ocean_tdrag_nml form = '"// &
                                 trim(adjustl(cfg%ocean%tdrag%form))// &
                                 "' is not recognised (quadratic|linear).  The two "// &
                                 "forms take coefficients of DIFFERENT dimensions "// &
                                 "(cd dimensionless, r in 1/s), so a typo cannot be "// &
                                 "defaulted.")
               has_error = .true.
            else if (tdrag_code == TDRAG_QUADRATIC .and. cfg%ocean%tdrag%cd <= 0.0_wp) then
               call logger%error("&ocean_tdrag_nml form='quadratic' requires cd > 0 "// &
                                 "(ISOMIP+ prescribes 2.5e-3).  cd = 0 is an "// &
                                 "identically zero drag, which is what enable=.false. "// &
                                 "is for.")
               has_error = .true.
            else if (tdrag_code == TDRAG_LINEAR .and. cfg%ocean%tdrag%r <= 0.0_wp) then
               call logger%error("&ocean_tdrag_nml form='linear' requires r > 0 "// &
                                 "(1/s).  r = 0 is an identically zero drag, which "// &
                                 "is what enable=.false. is for.")
               has_error = .true.
            end if
            ! ---- ONE drag coefficient for momentum and melt ----
            ! MOM6 carries two independent top-drag coefficients (one in
            ! the momentum BC, one in the melt u*); we deliberately do
            ! not.  A cavity in which the ice base takes momentum out of
            ! the flow at one C_d and reports a friction velocity built
            ! on another is not a closure, it is two closures sharing a
            ! boundary.  So: when both groups are on, the melt slot TAKES
            ! its C_d from this group (`configure_ocean_cavity_melt`), and
            ! a user who set both to DIFFERENT values is told rather than
            ! silently overridden.
            if (cfg%ocean%cavity_melt%enable .and. tdrag_code == TDRAG_QUADRATIC) then
               if (cfg%ocean%cavity_melt%cdrag_top /= cfg%ocean%tdrag%cd) then
                  call logger%error("&ocean_tdrag_nml cd = "// &
                                    to_string(cfg%ocean%tdrag%cd)//" and "// &
                                    "&ocean_cavity_melt_nml cdrag_top = "// &
                                    to_string(cfg%ocean%cavity_melt%cdrag_top)// &
                                    " disagree.  There is ONE ice-base drag "// &
                                    "coefficient in this model: the same C_d sets "// &
                                    "the momentum sink and the melt friction "// &
                                    "velocity u* = sqrt(C_d*(U^2 + u_tide^2)).  Set "// &
                                    "them equal (or leave cdrag_top at its default "// &
                                    "and set &ocean_tdrag_nml cd alone).")
                  has_error = .true.
               end if
            end if
            if (cfg%ocean%cavity_melt%enable .and. tdrag_code == TDRAG_LINEAR) then
               call logger%warning("&ocean_tdrag_nml form='linear' with "// &
                                   "&ocean_cavity_melt_nml enable=.true.: the melt "// &
                                   "friction velocity is quadratic by construction "// &
                                   "(u* = sqrt(C_d*(U^2 + u_tide^2))), so it keeps "// &
                                   "its own cdrag_top and the momentum sink uses r. "// &
                                   "The two boundary conditions are then NOT the "// &
                                   "same closure — intended for analytic work only.")
            end if
         end block
      end if
      ! ---- Dynamic wetting/drying v1 scope (docs/ocean_wetdry_plan.md §6) ----
      ! Every restriction fails loud: silently running wet/dry outside its
      ! validated envelope is the coastal ZSTAR_FULL salt-leak foot-gun class.
      if (cfg%ocean%wetdry%enable) then
         if (.not. (cfg%ocean%wetdry%rewet_depth > cfg%ocean%wetdry%dry_depth &
                    .and. cfg%ocean%wetdry%dry_depth > 0.0_wp)) then
            call logger%error("&ocean_wetdry_nml requires rewet_depth > dry_depth > 0 "// &
                              "(hysteresis band; equal thresholds flip-flop the front)")
            has_error = .true.
         end if
         ! v1 vertical-coordinate restriction: sigma only.  On sigma the
         ! layer partition scales all layers to zero TOGETHER as D -> 0, so
         ! a dry column is just "all layers vanished" (a state the
         ! H_VANISHED gates already handle).  ZSTAR_FULL bed layers pinch
         ! independently of the surface — the argument fails and the
         ! documented 1-2%/cycle intertidal salt leak would compound.
         ! 'zstar' was accepted while it was sigma under another name; as
         ! MOM6 z* it carries zstar_h_min bed fillers (below H_VANISHED)
         ! where the wet/dry seed and floor assume an emerged column of
         ! nz*2*H_VANISHED, and its dilation is floored at a dry column —
         ! unvalidated, so refused.
         ! Parse (not string-compare): the vcoord string has aliases
         ! ("z-star", "SIGMA", ...); same ocean default as rdb_ocean_vcoord.
         block
            integer :: wd_vcoord_code
            wd_vcoord_code = parse_vcoord_type(cfg%vcoord_type, &
                                               default_code=VCOORD_EULERIAN_Z)
            if (.not. (wd_vcoord_code == VCOORD_SIGMA)) then
               call logger%error("&ocean_wetdry_nml enable=.true. supports "// &
                                 "vcoord_type='sigma' only — got '"// &
                                 trim(cfg%vcoord_type)//"' (ZSTAR_FULL bed layers "// &
                                 "pinch independently; 'zstar' is MOM6 z* with bed "// &
                                 "fillers below the wet/dry emerged-column floor; "// &
                                 "other coords unvalidated)")
               has_error = .true.
            end if
         end block
         ! Positive-definite continuity composes its own per-layer outflux
         ! limiter; wet/dry owns the barotropic drying limiter.  The two
         ! positivity machineries have not been reconciled (v1) — fail loud
         ! rather than run both over the same faces.
         if (cfg%ocean%continuity%positive_definite) then
            call logger%error("&ocean_wetdry_nml enable=.true. is mutually "// &
                              "exclusive with &ocean_continuity_nml "// &
                              "positive_definite=.true. (wet/dry owns its own "// &
                              "barotropic limiter; composition deferred)")
            has_error = .true.
         end if
         ! Layer-side positivity is the continuity-PPM limiter's job — the BT
         ! limiter bounds the column, ppm_limit_pos bounds each layer.
         if (.not. cfg%ocean%continuity%ppm_limit_pos) then
            call logger%error("&ocean_wetdry_nml enable=.true. requires "// &
                              "&ocean_continuity_nml ppm_limit_pos=.true. "// &
                              "(per-layer positivity under drying transports)")
            has_error = .true.
         end if
         ! The BT_cont / upstream-h Pass-1 flux branches are not composed
         ! with the wet/dry upwind+limiter branch in v1 (each replaces the
         ! same face-thickness logic).
         if (cfg%ocean%bt%use_cont_type .or. cfg%ocean%bt%upstream_h_face) then
            call logger%error("&ocean_wetdry_nml enable=.true. is mutually "// &
                              "exclusive with &ocean_bt_nml use_cont_type / "// &
                              "upstream_h_face (all three replace the BT Pass-1 "// &
                              "face thickness; composition deferred)")
            has_error = .true.
         end if
         ! v1 single-rank: the theta/mask ghost exchange lands with the
         ! C-grid MPI halo work.
         if (cfg%px*cfg%py > 1) then
            call logger%error("&ocean_wetdry_nml enable=.true. is single-rank in "// &
                              "v1 (px*py = 1); the wet-mask/limiter halo exchange "// &
                              "is deferred to the C-grid MPI work")
            has_error = .true.
         end if
         ! Wet/dry needs the split solver: the limiter lives in the BT substep.
         if (cfg%ocean%bt%n_inner < 1 .and. .not. cfg%ocean%bt%auto_n_inner) then
            call logger%error("&ocean_wetdry_nml enable=.true. requires the "// &
                              "split-explicit solver (&ocean_bt_nml n_inner >= 1 "// &
                              "or auto_n_inner=.true.)")
            has_error = .true.
         end if
         ! Surface-flux composition: only the MAIN heat/salt deposit is
         ! dynamic-mask gated in v1.  sw_pen removes the surface deposit
         ! and re-adds a distributed profile — on a masked (dry) column
         ! it would remove heat that was never added; the restoring
         ! piston rate lambda = piston/h_top and geothermal's bed deposit
         ! blow up on a ~0-thickness sliver.  All three mis-compose:
         ! fail loud, defer the gated variants.
         if (cfg%ocean%thermo%sw_pen_frac > 0.0_wp) then
            call logger%error("&ocean_wetdry_nml enable=.true. does not yet "// &
                              "compose with shortwave penetration "// &
                              "(&ocean_thermo_nml sw_pen_frac > 0): the additive "// &
                              "correction would withdraw a surface deposit the "// &
                              "dry-column mask never made")
            has_error = .true.
         end if
         if (cfg%ocean%restore%enable_restore_temp .or. &
             cfg%ocean%restore%enable_restore_salt) then
            call logger%error("&ocean_wetdry_nml enable=.true. does not yet "// &
                              "compose with surface restoring (&ocean_restore_nml): "// &
                              "the piston rate lambda = piston/h_top diverges on a "// &
                              "drying column")
            has_error = .true.
         end if
         if (cfg%ocean%geothermal%enable) then
            call logger%error("&ocean_wetdry_nml enable=.true. does not yet "// &
                              "compose with geothermal heating "// &
                              "(&ocean_geothermal_nml): the bed deposit blows up "// &
                              "a dried column's vanished bed layer")
            has_error = .true.
         end if
      end if
      ! ---- Pseudo-salt: fail-loud exclusions ----
      ! Pseudo-salt's contract is "receives EXACTLY the operators salinity
      ! receives" (PR-28).  SSS piston restoring and sea-ice frazil/basal
      ! salt exchange are salinity sources this PR does NOT mirror into
      ! pseudo-salt — an un-mirrored source would make the deviation
      ! diagnostic (pseudo_salt - S) measure that source instead of the
      ! passive-transport-path error it exists to isolate.  Abort rather
      ! than silently ship a lying diagnostic.
      if (pseudo_salt_conflicts_restore(cfg%ocean%tracers%enable_pseudo_salt, &
                                        cfg%ocean%restore%enable_restore_salt)) then
         call logger%error("&ocean_tracers_nml enable_pseudo_salt=.true. does not "// &
                           "compose with &ocean_restore_nml enable_restore_salt=.true. "// &
                           "(SSS piston restoring): the deviation diagnostic would "// &
                           "measure the un-mirrored restoring term, not the "// &
                           "passive-transport-path error")
         has_error = .true.
      end if
      if (pseudo_salt_conflicts_ice(cfg%ocean%tracers%enable_pseudo_salt, &
                                    cfg%ocean%ice%enable)) then
         call logger%error("&ocean_tracers_nml enable_pseudo_salt=.true. does not "// &
                           "compose with &ocean_ice_nml enable=.true.: sea-ice "// &
                           "frazil/basal salt exchange is an un-mirrored salinity "// &
                           "source that would corrupt the deviation diagnostic")
         has_error = .true.
      end if
      if (pseudo_salt_needs_thermo_warning(cfg%ocean%tracers%enable_pseudo_salt, &
                                           cfg%ocean%thermo%enable_thermodynamics)) then
         call logger%warning("&ocean_tracers_nml enable_pseudo_salt=.true. with "// &
                             "&ocean_thermo_nml enable_thermodynamics=.false.: "// &
                             "pseudo-salt will never receive the surface salt flux "// &
                             "or KPP non-local mirrors, so it measures pure "// &
                             "advection/diffusion transport only — a legitimate "// &
                             "but narrower use of the diagnostic")
      end if
      ! ---- Sea-ice: whole-slot envelope (applies to enable=.true.) ----
      ! The sea ice runs on the ocean's decomposition: the category state,
      ! the EVP ice velocity and the blended surface stress are
      ! halo-exchanged (`engine_step_ice`) and, on a tripolar north fold,
      ! folded across the seam through the same px-aware fold dispatcher
      ! as the ocean state (A7, #143), so a folded grid runs on any px.
      if (cfg%ocean%ice%enable) then
         ! Sea ice is tuned for global, ice-covered basins — not for
         ! regional/coastal runs.  It runs every outer step and is NOT
         ! performance-tuned; expect a large throughput hit.  Warn loudly
         ! so it is never mistaken for a free rider on a coastal run.
         call logger%warning("=====================================================")
         call logger%warning("&ocean_ice_nml enable=.true.: SEA ICE IS ON.")
         call logger%warning("Sea ice is tuned for GLOBAL, ice-covered simulations "// &
                             "and is NOT performance-tuned. Expect it to be VERY "// &
                             "SLOW on regional/coastal domains where most of the "// &
                             "grid is ice-free.")
         call logger%warning("=====================================================")
      end if
      ! ---- Sea-ice PR 4b: category ice transport envelope ----
      ! `transport=.true.` is a request to run `ice_transport_step` every
      ! thermo window; the v1 envelope excludes configurations the kernel
      ! does not (yet) handle — fail loud rather than silently produce a
      ! wrong/unstable answer.
      if (cfg%ocean%ice%transport) then
         if (.not. cfg%ocean%ice%enable) then
            call logger%error("&ocean_ice_nml transport=.true. requires "// &
                              "enable=.true. (the ice slot must be live)")
            has_error = .true.
         end if
         if (cfg%ocean%ice%ncat <= 1) then
            call logger%error("&ocean_ice_nml transport=.true. requires ncat > 1 "// &
                              "(ncat==1 is the frozen legacy lumped mode — no ITD, "// &
                              "nothing to transport by category)")
            has_error = .true.
         end if
         ! The driver only reaches `ice_transport_step` inside the
         ! `enable_thermodynamics .and. is_thermo_step()` ice block, so
         ! `transport=.true.` with thermo off is a WORDLESS no-op — fail
         ! loud (house rule) instead.
         if (.not. cfg%ocean%thermo%enable_thermodynamics) then
            call logger%error("&ocean_ice_nml transport=.true. requires "// &
                              "&ocean_thermo_nml enable_thermodynamics=.true. "// &
                              "(the transport call fires only on the thermo cadence)")
            has_error = .true.
         end if
         ! Periodic edges and more than one rank are both fine: the CAS
         ! masses and riding tracers are halo-exchanged (or, on one rank,
         ! wrapped) at the top of every advective substep, and a seam face
         ! is never zeroed as a wall (`rdb_ice_transport`).
         if (cfg%ocean%ice%adv_substeps < 1) then
            call logger%error("&ocean_ice_nml adv_substeps must be >= 1")
            has_error = .true.
         end if
      end if

      ! ---- Sea-ice PR 5: C-grid EVP dynamics envelope ----
      ! `dynamics=.true.` requests `ice_evp_step` every outer step; the v1
      ! envelope excludes configurations the kernel does not (yet) handle —
      ! fail loud rather than silently produce a wrong/unstable answer.
      ! EVP allows periodic edges (the ghost-wrap machinery in `rdb_ice_evp`
      ! mirrors the ocean's own periodic-wrap contract, as transport now
      ! does too); its edge envelope is its own, not a re-use of the
      ! transport block above.
      ! Multi-rank: the ice velocity is halo-exchanged every subcycle.
      if (cfg%ocean%ice%dynamics) then
         if (.not. cfg%ocean%ice%enable) then
            call logger%error("&ocean_ice_nml dynamics=.true. requires "// &
                              "enable=.true. (the ice slot must be live)")
            has_error = .true.
         end if
         ! v1 envelope: every edge must be WALL or PERIODIC — no OBC/tidal/
         ! sponge/clamped/Chapman edges for ice (the momentum solve has no
         ! notion of an open ice boundary yet).
         if (.not. any(ocean_bc_type_from_string(cfg%ocean%bc%west) == &
                       [OBC_WALL, OBC_PERIODIC]) .or. &
             .not. any(ocean_bc_type_from_string(cfg%ocean%bc%east) == &
                       [OBC_WALL, OBC_PERIODIC]) .or. &
             .not. any(ocean_bc_type_from_string(cfg%ocean%bc%south) == &
                       [OBC_WALL, OBC_PERIODIC]) .or. &
             .not. any(ocean_bc_type_from_string(cfg%ocean%bc%north) == &
                       [OBC_WALL, OBC_PERIODIC, OBC_TRIPOLAR_FOLD])) then
            call logger%error("&ocean_ice_nml dynamics=.true. requires every "// &
                              "&ocean_bc_nml edge to be 'wall' or 'periodic' (the north "// &
                              "edge may also be 'tripolar_fold'; v1: no "// &
                              "OBC/tidal/sponge/clamped/Chapman edges for ice dynamics)")
            has_error = .true.
         end if
         ! Tripolar north fold (A6, #164; multi-rank A7, #143): EVP folds the
         ! ice velocity every subcycle and recomputes its stencil one row past
         ! the seam, which needs two ghost rows.
         if (trim(cfg%ocean%bc%north) == "tripolar_fold" .and. cfg%nghost < 2) then
            call logger%error("&ocean_ice_nml dynamics=.true. with north='tripolar_fold' "// &
                              "needs nghost >= 2 (the EVP stencil reads one row past "// &
                              "the fold line); nghost = "//to_string(cfg%nghost))
            has_error = .true.
         end if
         if (cfg%ocean%ice%evp_sub_steps < 1) then
            call logger%error("&ocean_ice_nml evp_sub_steps must be >= 1")
            has_error = .true.
         end if
         ! Physical-parameter positivity guards — fail loud on a typo/garbage
         ! value that would silently produce a wrong or NaN rheology. `ec`
         ! allows 0 (the documented cavitating-fluid mode: `i_ec2` is forced
         ! to 0 in `rdb_ice_evp`), so it is guarded `>= 0`, not `> 0`.
         ! `tdamp` is legitimately sign-free (`< 0` selects the
         ! fraction-of-slow-step form) — no guard.
         if (cfg%ocean%ice%p0 <= 0.0_wp) then
            call logger%error("&ocean_ice_nml p0 (ice-strength P*) must be > 0")
            has_error = .true.
         end if
         if (cfg%ocean%ice%c0 <= 0.0_wp) then
            call logger%error("&ocean_ice_nml c0 (ice-strength C*) must be > 0")
            has_error = .true.
         end if
         if (cfg%ocean%ice%ec < 0.0_wp) then
            call logger%error("&ocean_ice_nml ec (yield ellipticity) must be >= 0 "// &
                              "(0 = cavitating-fluid rheology)")
            has_error = .true.
         end if
         if (cfg%ocean%ice%cdw <= 0.0_wp) then
            call logger%error("&ocean_ice_nml cdw (ice-ocean drag coefficient) must be > 0")
            has_error = .true.
         end if
         if (cfg%ocean%ice%rho_ocean <= 0.0_wp) then
            call logger%error("&ocean_ice_nml rho_ocean (ice-drag reference density) "// &
                              "must be > 0")
            has_error = .true.
         end if
         if (cfg%ocean%ice%del_sh_min_scale <= 0.0_wp) then
            call logger%error("&ocean_ice_nml del_sh_min_scale (viscosity-floor scale) "// &
                              "must be > 0")
            has_error = .true.
         end if
         ! Drift-only mode (dynamics without transport) has no CFL backstop
         ! on u_ice: transport's positivity abort is the only ice-velocity
         ! guard. Warn (NOT error — the envelope is unchanged); this
         ! combination is unvalidated.
         if (.not. cfg%ocean%ice%transport) then
            call logger%warning("&ocean_ice_nml dynamics=.true. without "// &
                                "transport=.true.: drift-only mode has no CFL backstop "// &
                                "on u_ice (transport's positivity abort is the only "// &
                                "ice-velocity guard) and is unvalidated")
         end if
         ! PR 36: the per-iteration clip needs a ceiling to clip to — SIS2's
         ! own triple gate (`do_trunc_its = CFL_check_its .and. CFL_trunc>0
         ! .and. dt_slow>0`) makes this combination a silent no-op; fail
         ! loud instead of silently reproducing SIS2's silence.
         if (cfg%ocean%ice%cfl_trunc_dyn_its .and. cfg%ocean%ice%cfl_trunc <= 0.0_wp) then
            call logger%error("&ocean_ice_nml cfl_trunc_dyn_its=.true. requires "// &
                              "cfl_trunc > 0 (the per-iteration check needs a ceiling)")
            has_error = .true.
         end if
         ! SIS2 permits cfl_trunc > 1 but calls it unwise ("instability can
         ! occur past 0.5") — warn, do not fail loud.
         if (cfg%ocean%ice%cfl_trunc > 1.0_wp) then
            call logger%warning("&ocean_ice_nml cfl_trunc > 1.0 is permitted but unwise "// &
                                "(SIS2: 'instability can occur past 0.5')")
         end if
      end if

      ! ---- Sea-ice PR 36: cfl_trunc / cfl_trunc_dyn_its / project_ci ----
      ! require dynamics; a knob that only has meaning inside a disabled
      ! feature is a wordless no-op — fail loud instead (same pattern as
      ! a_face_stress below). Deliberately OUTSIDE the `if (dynamics)`
      ! envelope above: it must fire precisely when `dynamics` is false.
      if (cfg%ocean%ice%cfl_trunc > 0.0_wp .and. .not. cfg%ocean%ice%dynamics) then
         call logger%error("&ocean_ice_nml cfl_trunc > 0 requires dynamics=.true. "// &
                           "(the knob only reaches the EVP velocity solve)")
         has_error = .true.
      end if
      if (cfg%ocean%ice%cfl_trunc_dyn_its .and. .not. cfg%ocean%ice%dynamics) then
         call logger%error("&ocean_ice_nml cfl_trunc_dyn_its=.true. requires "// &
                           "dynamics=.true. (the knob only reaches the EVP velocity solve)")
         has_error = .true.
      end if
      if (cfg%ocean%ice%project_ci .and. .not. cfg%ocean%ice%dynamics) then
         call logger%error("&ocean_ice_nml project_ci=.true. requires dynamics=.true. "// &
                           "(the knob only reaches the EVP subcycle loop)")
         has_error = .true.
      end if

      ! ---- Map-driven sponge (&ocean_sponge_nml, PR-23) ----
      if (cfg%ocean%sponge%enable) then
         if (.not. sponge_source_is_implemented(cfg%ocean%sponge%damp_source)) then
            call logger%error("&ocean_sponge_nml enable=.true. with unimplemented "// &
                              "damp_source = '"//trim(cfg%ocean%sponge%damp_source)// &
                              "': must be 'band' ('file' is PR-23b, needs the PR-14 reader)")
            has_error = .true.
         end if
         if (.not. sponge_target_is_implemented(cfg%ocean%sponge%target_source)) then
            call logger%error("&ocean_sponge_nml enable=.true. with unimplemented "// &
                              "target_source = '"//trim(cfg%ocean%sponge%target_source)// &
                              "': must be 'ic' or 'linear_z' "// &
                              "('file' is PR-23b, needs the PR-14 reader)")
            has_error = .true.
         end if
         if (.not. sponge_ramp_is_valid(cfg%ocean%sponge%ramp)) then
            call logger%error("&ocean_sponge_nml ramp = '"//trim(cfg%ocean%sponge%ramp)// &
                              "': must be 'cosine' (default, the legacy band shape) "// &
                              "or 'linear' (ISOMIP+ Eq. 20)")
            has_error = .true.
         end if
         ! A `linear_z` target with both gradients AND both references left
         ! at their defaults would relax toward T = 0 degC / S = 35 PSU
         ! everywhere — almost certainly not what was meant, and silent.
         if (trim(cfg%ocean%sponge%target_source) == "linear_z" .and. &
             cfg%ocean%sponge%lin_t_ref == 0.0_wp .and. &
             cfg%ocean%sponge%lin_dt_dz == 0.0_wp .and. &
             cfg%ocean%sponge%lin_ds_dz == 0.0_wp) then
            call logger%warning("&ocean_sponge_nml target_source='linear_z' with "// &
                                "lin_t_ref = lin_dt_dz = lin_ds_dz = 0: the sponge will "// &
                                "relax toward a uniform T = 0 degC column. Set the "// &
                                "lin_* profile knobs.")
         end if
         ! relax_h (interior-interface thickness damping) is NOT implemented
         ! in PR-23 v1 — deferred to PR-23b alongside the file targets (the
         ! plan's §14 Q1 recommendation). Abort unconditionally rather than
         ! silently ignore the knob.
         if (cfg%ocean%sponge%relax_h) then
            call logger%error("&ocean_sponge_nml relax_h=.true. is not implemented "// &
                              "in PR-23 v1 (interior-interface thickness damping is "// &
                              "deferred to PR-23b; see docs/plans/PLAN_PR23_real_sponge.md "// &
                              "§14 Q1) — leave relax_h=.false.")
            has_error = .true.
         end if
         ! The sponge is dispatched only from run_stage_split (the split
         ! barotropic solver); the unsplit run_stage never receives bc or sp.
         ! With auto_n_inner=.true. n_inner resolves later at setup, so only
         ! a HARD n_inner < 1 with auto off is an error here (mirrors the
         ! &ocean_wetdry_nml n_inner precedent above).
         if (cfg%ocean%bt%n_inner < 1 .and. .not. cfg%ocean%bt%auto_n_inner) then
            call logger%error("&ocean_sponge_nml enable=.true. requires the "// &
                              "split-explicit solver (&ocean_bt_nml n_inner >= 1 "// &
                              "or auto_n_inner=.true.) — the sponge is dispatched only "// &
                              "from the split barotropic path")
            has_error = .true.
         end if
         ! A "band" idamp source with no OBC_SPONGE-tagged edge would build
         ! an all-zero map — the sponge would silently do nothing.
         if (trim(cfg%ocean%sponge%damp_source) == "band") then
            if (trim(cfg%ocean%bc%west) /= "sponge" .and. &
                trim(cfg%ocean%bc%east) /= "sponge" .and. &
                trim(cfg%ocean%bc%south) /= "sponge" .and. &
                trim(cfg%ocean%bc%north) /= "sponge") then
               call logger%error("&ocean_sponge_nml enable=.true., damp_source='band' "// &
                                 "requires at least one &ocean_bc_nml edge = 'sponge' "// &
                                 "(otherwise idamp_h/u/v would build all-zero)")
               has_error = .true.
            end if
         end if
         ! Per-edge width/strength overrides: sentinel is < 0 (inherit); a
         ! non-sentinel width must be non-negative, and a positive width
         ! paired with a zero strength is a silent no-op sponge.
         if (cfg%ocean%sponge%west_width < -1 .or. cfg%ocean%sponge%east_width < -1 .or. &
             cfg%ocean%sponge%south_width < -1 .or. cfg%ocean%sponge%north_width < -1) then
            call logger%error("&ocean_sponge_nml *_width overrides must be >= -1 "// &
                              "(-1 = inherit &ocean_bc_nml sponge_width)")
            has_error = .true.
         end if
         if ((cfg%ocean%sponge%west_width > 0 .and. cfg%ocean%sponge%west_strength == 0.0_wp) .or. &
             (cfg%ocean%sponge%east_width > 0 .and. cfg%ocean%sponge%east_strength == 0.0_wp) .or. &
             (cfg%ocean%sponge%south_width > 0 .and. cfg%ocean%sponge%south_strength == 0.0_wp) .or. &
             (cfg%ocean%sponge%north_width > 0 .and. cfg%ocean%sponge%north_strength == 0.0_wp)) then
            call logger%error("&ocean_sponge_nml a positive *_width override with "// &
                              "*_strength == 0.0 would build an all-zero band on that "// &
                              "edge — set a nonzero *_strength or leave *_width at -1")
            has_error = .true.
         end if
      end if
      ! ---- Barotropic linear wave drag (&ocean_bt_nml wave_drag) ----
      ! `wave_drag_form` is an `nml_enum` (allowed = uniform/roughness_proxy/
      ! file), so a garbage tag is already rejected at nml-apply time; the
      ! one form-specific case left to catch here is "file", which is
      ! REGISTERED but has no reader yet (PR-14) — fail loud rather than
      ! silently no-op (the audit's trap #12: a named-but-dead option).
      if (trim(cfg%ocean%bt%wave_drag_form) == "file") then
         call logger%error("&ocean_bt_nml wave_drag_form='file' is not "// &
                           "implemented — the NetCDF map reader is PR-14; "// &
                           "use wave_drag_form='uniform' or 'roughness_proxy'")
         has_error = .true.
      end if
      if (cfg%ocean%bt%wave_drag_scale < 0.0_wp) then
         call logger%error("&ocean_bt_nml wave_drag_scale must be >= 0 "// &
                           "(negative r_H makes bt_rem > 1, exponentially "// &
                           "growing the barotropic mode)")
         has_error = .true.
      end if
      if (cfg%ocean%bt%wave_drag_r_uniform < 0.0_wp) then
         call logger%error("&ocean_bt_nml wave_drag_r_uniform must be >= 0")
         has_error = .true.
      end if
      if (cfg%ocean%bt%wave_drag_kappa <= 0.0_wp) then
         call logger%error("&ocean_bt_nml wave_drag_kappa must be > 0")
         has_error = .true.
      end if
      if (cfg%ocean%bt%wave_drag_n_bot < 0.0_wp) then
         call logger%error("&ocean_bt_nml wave_drag_n_bot must be >= 0")
         has_error = .true.
      end if
      if (cfg%ocean%bt%wave_drag_h2_max <= 0.0_wp) then
         call logger%error("&ocean_bt_nml wave_drag_h2_max must be > 0")
         has_error = .true.
      end if
      if (cfg%ocean%bt%wave_drag) then
         if (trim(cfg%ocean%bt%wave_drag_form) == "uniform" .and. &
             cfg%ocean%bt%wave_drag_r_uniform <= 0.0_wp) then
            call logger%warning("&ocean_bt_nml wave_drag=.true. with "// &
                                "wave_drag_form='uniform' and wave_drag_r_uniform <= 0 "// &
                                "=> inert (r_H == 0 everywhere)")
         end if
         if (trim(cfg%ocean%bt%wave_drag_form) == "roughness_proxy") then
            call logger%warning("&ocean_bt_nml wave_drag_form='roughness_proxy' derives "// &
                                "<h^2> from RESOLVED bathymetry variance — it is a "// &
                                "placeholder for a real subgrid roughness map "// &
                                "(PR-14/PR-30), not a substitute")
         end if
      end if
      ! ---- Sea-ice PR-58: ITD category-bound override envelope ----
      ! `hlim=...` is consumed entirely inside `ice%init` (called from
      ! `ocean_state_init_from_config`, BEFORE the thermo cadence), so —
      ! unlike transport/dynamics — there is no `enable_thermodynamics`
      ! clause to copy here: a wordless no-op would require the value to
      ! be read on a cadence it never reaches, and `hlim` has no such gate.
      if (ice_hlim_count(cfg%ocean%ice%hlim) > 0) then
         if (.not. cfg%ocean%ice%enable) then
            call logger%error("&ocean_ice_nml hlim=... requires enable=.true. "// &
                              "(the ice slot must be live for the override to apply)")
            has_error = .true.
         end if
         block
            logical :: hlim_ok
            character(len=:), allocatable :: reason
            call ice_hlim_spec_is_valid(cfg%ocean%ice%hlim, cfg%ocean%ice%ncat, hlim_ok, reason)
            if (.not. hlim_ok) then
               call logger%error("&ocean_ice_nml hlim: "//reason)
               has_error = .true.
            end if
         end block
      end if
      ! ---- Sea-ice PR 24: analytic initial-condition envelope ----
      ! `conc_config /= "zero"` is a request to seed a live ice pack
      ! before the first step; the envelope excludes configurations the
      ! seeder cannot handle consistently — fail loud rather than
      ! silently seed a meaningless or self-melting pack.
      if (trim(cfg%ocean%ice_ic%conc_config) /= "zero") then
         ! V1: mirrors "transport=.true. requires enable=.true." above.
         if (.not. cfg%ocean%ice%enable) then
            call logger%error("&ocean_ice_ic_nml conc_config /= 'zero' requires "// &
                              "&ocean_ice_nml enable=.true. (the ice slot must be live)")
            has_error = .true.
         end if
         ! V2: an IC with no thickness is a wordless no-op; h_ice > 0
         ! also keeps the pack above h_lim(1) = 1e-10 m so the ITD cat-1
         ! compress cannot fire on the very first thermo window.
         if (cfg%ocean%ice_ic%h_ice <= 0.0_wp) then
            call logger%error("&ocean_ice_ic_nml conc_config /= 'zero' requires "// &
                              "h_ice > 0 (an IC with no thickness is a wordless no-op, "// &
                              "and h_ice must clear the category-1 floor)")
            has_error = .true.
         end if
         ! V3: mirrors the "wordless no-op => fail loud" argument used
         ! for transport above.
         if (trim(cfg%ocean%ice_ic%conc_config) == "uniform" .and. &
             cfg%ocean%ice_ic%conc <= 0.0_wp) then
            call logger%error("&ocean_ice_ic_nml conc_config='uniform' requires "// &
                              "conc > 0 (conc <= 0 is a wordless open-water no-op)")
            has_error = .true.
         end if
         ! V4: copy of the equilibrium-tide cartesian-grid guard above —
         ! geolatT/geolonT stay 0 on a cartesian grid (only
         ! metrics_init_spherical fills them), so a latitudes IC would
         ! silently ice the whole basin or none of it.
         if (trim(cfg%ocean%ice_ic%conc_config) == "latitudes" .and. &
             trim(cfg%ocean%grid%grid_config) == "cartesian") then
            call logger%error("&ocean_ice_ic_nml conc_config='latitudes' requires a "// &
                              "non-cartesian &ocean_grid_nml grid_config (spherical/"// &
                              "tripolar/supergrid) — geolatT is meaningless on cartesian")
            has_error = .true.
         end if
         ! V5: the legacy lumped mode (ncat==1) has no partial-cover
         ! representation — ice_cell_concentration_impl returns a BINARY
         ! ci, so a fractional conc would be silently reinterpreted.
         if (cfg%ocean%ice%ncat == 1 .and. &
             trim(cfg%ocean%ice_ic%conc_config) == "uniform" .and. &
             cfg%ocean%ice_ic%conc /= 1.0_wp) then
            call logger%error("&ocean_ice_ic_nml conc_config='uniform' with "// &
                              "&ocean_ice_nml ncat==1 requires conc == 1.0 (the legacy "// &
                              "lumped mode has no partial-cover representation — "// &
                              "ice_cell_concentration_impl returns a BINARY ci)")
            has_error = .true.
         end if
         ! V6: at/above the liquidus, ice_enth_from_ts returns the
         ! LIQUID-WATER branch — the pack would carry m_ice > 0 whose
         ! enthalpy says "this is water" and melt out on the first
         ! thermo step, silently.
         if (cfg%ocean%ice_ic%t_ice >= ice_t_freeze(cfg%ocean%ice_ic%s_ice)) then
            call logger%error("&ocean_ice_ic_nml conc_config /= 'zero' requires "// &
                              "t_ice < the liquidus at s_ice (t_ice >= t_freeze(s_ice) "// &
                              "silently yields LIQUID-WATER enthalpy — the pack would "// &
                              "melt out on the first thermo step)")
            has_error = .true.
         end if
         ! V7: orphan snow. The column entry gate is on m_ice, and
         ! ice_adjust_categories' massless-cleanup step would silently
         ! absorb snow-with-no-ice into the "massless" taxonomy.
         if (cfg%ocean%ice_ic%h_snow > 0.0_wp .and. cfg%ocean%ice_ic%h_ice <= 0.0_wp) then
            call logger%error("&ocean_ice_ic_nml h_snow > 0 requires h_ice > 0 "// &
                              "(orphan snow with no ice is silently absorbed as "// &
                              "'massless' by the ITD restore)")
            has_error = .true.
         end if
         ! V8: belt-and-braces on top of nml_enum allowed= — catches an
         ! unrecognised conc_config string even if the schema layer is
         ! ever bypassed.
         if (ice_ic_parse_conc_config(cfg%ocean%ice_ic%conc_config) == ICE_IC_CONC_INVALID) then
            call logger%error("&ocean_ice_ic_nml conc_config='"// &
                              trim(cfg%ocean%ice_ic%conc_config)// &
                              "' is not recognised (allowed: zero, uniform, latitudes)")
            has_error = .true.
         end if
         ! V9: the IC is skipped on resume (a restart REPLACES it, never
         ! merges) — warn, don't abort, so the legitimate "same nml, cold
         ! then warm" workflow still works.
         if (len_trim(cfg%restart_file) > 0) then
            call logger%warning("&ocean_ice_ic_nml conc_config /= 'zero' with "// &
                                "restart_file set: the ice IC is SKIPPED on a warm "// &
                                "restart (the restart replaces it, it does not merge "// &
                                "with it) — this namelist's ice IC knobs have no effect "// &
                                "on this run")
         end if
      end if
      ! ---- Sea-ice PR 26: snowfall source-term envelope ----
      ! `snowfall > 0` is a request for `ice_atm_forcing_restoring` +
      ! `ice_thermo_driver_step`/`ice_snowfall_ocean_share` to fire every
      ! thermo window; both live inside the `enable_thermodynamics .and.
      ! is_thermo_step()` ice block in the driver, so a config with the
      ! ice slot or thermodynamics off would leave `snowfall` a wordless
      ! no-op — fail loud instead (house rule, mirrors the transport
      ! envelope above).
      if (cfg%ocean%ice%snowfall > 0.0_wp) then
         if (.not. cfg%ocean%ice%enable) then
            call logger%error("&ocean_ice_nml snowfall > 0 requires "// &
                              "enable=.true. (the ice slot must be live)")
            has_error = .true.
         end if
         if (.not. cfg%ocean%thermo%enable_thermodynamics) then
            call logger%error("&ocean_ice_nml snowfall > 0 requires "// &
                              "&ocean_thermo_nml enable_thermodynamics=.true. "// &
                              "(the atmospheric-forcing fill and the column "// &
                              "snow add both fire only on the thermo cadence)")
            has_error = .true.
         end if
      end if

      ! ---- Sea-ice PR 62: a_face_stress requires dynamics ----
      ! `a_face_stress` only reaches the EVP momentum solve — with
      ! `dynamics=.false.` it would be a wordless no-op.  Fail loud instead
      ! (house rule).  Deliberately OUTSIDE the `if (dynamics)` envelope
      ! above: it must fire precisely when `dynamics` is false.
      if (cfg%ocean%ice%a_face_stress .and. .not. cfg%ocean%ice%dynamics) then
         call logger%error("&ocean_ice_nml a_face_stress=.true. requires dynamics=.true. "// &
                           "(the knob only reaches the EVP momentum solve)")
         has_error = .true.
      end if

      ! ---- Sea-ice PR 27: Archimedes snow-ice flooding envelope ----
      ! `snow_ice=.true.` is a request for `ice_snow_ice_flood` to fire
      ! inside `ice_column_step`'s resize stage every thermo window; that
      ! call lives inside the `enable_thermodynamics .and. is_thermo_
      ! step()` ice block in the driver, so a config with the ice slot or
      ! thermodynamics off would leave `snow_ice` a wordless no-op — fail
      ! loud instead (house rule, mirrors the transport/snowfall envelopes
      ! above).
      if (cfg%ocean%ice%snow_ice) then
         if (.not. cfg%ocean%ice%enable) then
            call logger%error("&ocean_ice_nml snow_ice=.true. requires "// &
                              "enable=.true. (the ice slot must be live)")
            has_error = .true.
         end if
         if (.not. cfg%ocean%thermo%enable_thermodynamics) then
            call logger%error("&ocean_ice_nml snow_ice=.true. requires "// &
                              "&ocean_thermo_nml enable_thermodynamics=.true. "// &
                              "(the column runs only on the thermo cadence)")
            has_error = .true.
         end if
      end if

      ! ---- BT wide-halo march-in (bt_halo > 0) exclusions ----
      if (cfg%ocean%bt%bt_halo > 0) then
         if (cfg%ocean%wetdry%enable) then
            call logger%error("&ocean_bt_nml bt_halo > 0 is mutually exclusive "// &
                              "with &ocean_wetdry_nml enable (wd arrays not widened in v1)")
            has_error = .true.
         end if
         if (cfg%ocean%bt%use_cont_type) then
            call logger%error("&ocean_bt_nml bt_halo > 0 is mutually exclusive "// &
                              "with use_cont_type (BTCL arrays not widened in v1)")
            has_error = .true.
         end if
         if (cfg%ocean%bt%upstream_h_face) then
            call logger%error("&ocean_bt_nml bt_halo > 0 is mutually exclusive "// &
                              "with upstream_h_face (upstream h-face arrays not widened in v1)")
            has_error = .true.
         end if
         if (cfg%ocean%tides%enable) then
            call logger%error("&ocean_bt_nml bt_halo > 0 is mutually exclusive "// &
                              "with &ocean_tides_nml enable (wide eta_forcing copy deferred)")
            has_error = .true.
         end if
         if (cfg%ocean%psurf%enable) then
            call logger%error("&ocean_bt_nml bt_halo > 0 is mutually exclusive "// &
                              "with &ocean_psurf_nml enable (wide eta_forcing copy deferred)")
            has_error = .true.
         end if
         if (cfg%ocean%porous%enable) then
            ! `bt_wide` carries its OWN `metrics_w`, re-filled from the
            ! grid formula on the widened grid.  Nothing fills its porous
            ! statistics, so `bt_wide_substep` would hand the fast loop the
            ! UN-narrowed `metrics_w%dy_cu` / `dx_cv` — silently reverting
            ! the porous-aware barotropic solve and giving an answer that
            ! is neither the porous one nor the baseline.  Fail loud until
            ! the wide shadow carries `dy_cu_bt` / `dx_cv_bt` too.
            call logger%error("&ocean_bt_nml bt_halo > 0 is mutually exclusive "// &
                              "with &ocean_porous_nml enable (the wide-halo "// &
                              "metrics shadow carries no porous statistics, so the "// &
                              "BT solve would silently transport on un-narrowed widths)")
            has_error = .true.
         end if
         if (cfg%ocean%cavity_dyn%enable) then
            ! Same failure mode as porous, one level up: `metrics_w` is
            ! re-filled from the grid formula and carries no `z_draft`, so
            ! the wide fast loop would take the BED as its reference depth
            ! and solve a column that is twice as deep as the cavity's.
            ! (The cavity block above refuses this from its own side too —
            ! the knob that is "wrong" depends on which one the user meant.)
            call logger%error("&ocean_bt_nml bt_halo > 0 is mutually exclusive "// &
                              "with &ocean_cavity_dyn_nml enable (the wide-halo "// &
                              "metrics shadow carries no ice draft, so the BT solve "// &
                              "would reference the bed instead of the ice base)")
            has_error = .true.
         end if
         if (trim(cfg%ocean%grid%grid_config) == "supergrid" .or. &
             trim(cfg%ocean%grid%grid_config) == "tripolar") then
            call logger%error("&ocean_bt_nml bt_halo > 0 is mutually exclusive "// &
                              "with grid_config='"//trim(cfg%ocean%grid%grid_config)// &
                              "' (file-based metrics cannot be re-filled on a wide grid; "// &
                              "cartesian/spherical formula fills supported)")
            has_error = .true.
         end if
         if (ocean_bt_rem_from_visc_rem_on(cfg)) then
            ! Same check as above, from the bt_halo side — belt-and-
            ! braces so the error fires regardless of which knob a reader
            ! finds first in the namelist.  PR-3 (D1): reads through the
            ! helper, so visc_rem_chain = .true. is caught too.
            call logger%error("&ocean_bt_nml bt_halo > 0 is mutually exclusive "// &
                              "with bt_rem_from_visc_rem=.true. (or "// &
                              "visc_rem_chain=.true.; the wide-halo BT clone carries "// &
                              "no av_rem/visc_rem ghost-width statistics, like porous)")
            has_error = .true.
         end if
      end if
      ! ---- Porous barriers (&ocean_porous_nml) exclusions ----
      if (cfg%ocean%porous%enable) then
         if (cfg%ocean%wetdry%enable) then
            ! A drying column drives every layer's face thickness below
            ! H_VANISHED, where the open fractions all go to zero and the
            ! barotropic width goes with them.  That is self-consistent but
            ! wholly untested against the wet/dry outflow limiter, and the
            ! two schemes each own part of the same "can this face carry
            ! transport" decision.  Refuse the combination rather than ship
            ! an unvalidated interaction.
            call logger%error("&ocean_porous_nml enable is mutually exclusive with "// &
                              "&ocean_wetdry_nml enable (the vanishing-column "// &
                              "interaction between the open fractions and the "// &
                              "wet/dry outflow limiter is unvalidated)")
            has_error = .true.
         end if
      end if

      ! ---- Ideal-age tracer knob sanity (PR-7, fail-quiet class) ----
      ! Neither combination is unsafe, so both are warnings, not aborts —
      ! and the default (0/0) must not trip either guard.
      if (cfg%ocean%tracers%ideal_age_sfc_growth_rate /= 0.0_wp .and. &
          cfg%ocean%tracers%ideal_age_young_val == 0.0_wp) then
         call logger%warning("&ocean_tracers_nml ideal_age_sfc_growth_rate is nonzero but "// &
                             "ideal_age_young_val = 0: young_val*exp(...) is identically 0, "// &
                             "so the growth-rate knob is silently inert (MOM6 seeds its "// &
                             "vintage tracer with young_val = 1e-20 for exactly this reason)")
      end if
      if (.not. cfg%ocean%tracers%enable_ideal_age .and. &
          (cfg%ocean%tracers%ideal_age_young_val /= 0.0_wp .or. &
           cfg%ocean%tracers%ideal_age_sfc_growth_rate /= 0.0_wp)) then
         call logger%warning("&ocean_tracers_nml ideal_age_young_val / "// &
                             "ideal_age_sfc_growth_rate are set but enable_ideal_age = "// &
                             ".false.: no age tracer is registered, so both knobs are "// &
                             "silently inert")
      end if

      ! ---- Surface-flux component set (PR-12) ----
      ! The component set only feeds the thermo path (the assembler
      ! derives Q_heat/Q_salt, which apply_tracers stamps onto the
      ! tracer hTr); with thermo off the whole thing is inert. Fail loud
      ! rather than let a user carry the (real, gated-off) device memory
      ! for nothing.
      if (forcing_components_need_thermo(cfg%ocean%forcing%enable_components, &
                                         cfg%ocean%thermo%enable_thermodynamics)) then
         call logger%error("&ocean_forcing_nml enable_components=.true. requires "// &
                           "&ocean_thermo_nml enable_thermodynamics=.true. "// &
                           "(the component set feeds the tracer path; it is inert "// &
                           "under enable_thermodynamics=.false.)")
         has_error = .true.
      end if

      ! ---- &ocean_dataovr_nml (PR-15 file-backed surface forcing) ----
      !
      ! Validated HERE, not at registration, on the "a broken config must
      ! not run at all" principle: `ocean_data_forcing_configure` runs
      ! deep in the driver's setup phase, AFTER the grid, bathymetry,
      ! initial condition and a dozen slots are built.  Catching a
      ! namelist typo there means the user pays for all of that first.
      ! Everything below is decidable from `cfg` alone, so it costs
      ! nothing to decide it before any work happens.
      !
      ! The fail-loud guards inside `rdb_ocean_data_forcing` STAY: they
      ! are the gate for a programmatic caller (a future Python driver
      ! registers fields directly and never parses a namelist), and this
      ! block is unreachable for that path.
      if (cfg%ocean%dataovr%enable) then
#ifdef RDB_NO_NETCDF
         call logger%error("&ocean_dataovr_nml enable=.true. requires a NetCDF build "// &
                           "(RDB_ENABLE_NETCDF=ON); the reader is compiled out here")
         has_error = .true.
#endif
         ! D1: a named file with no variable name.  No guessed fallbacks
         ! for forcing — an unintended variable is worse than a stop.
         if (.not. dataovr_entry_is_valid(cfg%ocean%dataovr%tau_x)) then
            call logger%error("&ocean_dataovr_nml tau_x_file is set but tau_x_var is blank")
            has_error = .true.
         end if
         if (.not. dataovr_entry_is_valid(cfg%ocean%dataovr%tau_y)) then
            call logger%error("&ocean_dataovr_nml tau_y_file is set but tau_y_var is blank")
            has_error = .true.
         end if
         if (.not. dataovr_entry_is_valid(cfg%ocean%dataovr%heat)) then
            call logger%error("&ocean_dataovr_nml heat_file is set but heat_var is blank")
            has_error = .true.
         end if
         if (.not. dataovr_entry_is_valid(cfg%ocean%dataovr%evap)) then
            call logger%error("&ocean_dataovr_nml evap_file is set but evap_var is blank")
            has_error = .true.
         end if
         if (.not. dataovr_entry_is_valid(cfg%ocean%dataovr%lprec)) then
            call logger%error("&ocean_dataovr_nml lprec_file is set but lprec_var is blank")
            has_error = .true.
         end if
         if (.not. dataovr_entry_is_valid(cfg%ocean%dataovr%salt)) then
            call logger%error("&ocean_dataovr_nml salt_file is set but salt_var is blank")
            has_error = .true.
         end if

         ! D2: cyclic climatology with no period is undecidable, not
         ! defaultable — there is no sane fallback for "one year".
         if (.not. dataovr_time_is_valid(cfg%ocean%dataovr)) then
            call logger%error("&ocean_dataovr_nml time_mode='cyclic' requires cycle_period > 0")
            has_error = .true.
         end if

         ! D3: the freshwater tags exist ONLY in the surface-flux
         ! component set; without it there is no slot to write into.
         if (dataovr_freshwater_needs_components(cfg%ocean%dataovr, &
                                                 cfg%ocean%forcing%enable_components)) then
            call logger%error("&ocean_dataovr_nml evap/lprec file forcing requires "// &
                              "&ocean_forcing_nml enable_components=.true. "// &
                              "(there is no non-component freshwater slot to write into)")
            has_error = .true.
         end if

         ! D4: a named file that does not exist.  A typo'd path is the
         ! single most common way to get this group wrong, and without
         ! this check it surfaces only once registration opens the file
         ! — deep in setup, as a bare "NetCDF operation failed" that
         ! does not even name the path.  `inquire` needs no NetCDF, so
         ! this works in every build.  It is a genuine
         ! time-of-check/time-of-use race (the file could vanish before
         ! registration), which is fine: registration still fails loud,
         ! this only moves the COMMON case earlier and makes it legible.
         call check_dataovr_file(cfg%ocean%dataovr%tau_x, "tau_x", has_error)
         call check_dataovr_file(cfg%ocean%dataovr%tau_y, "tau_y", has_error)
         call check_dataovr_file(cfg%ocean%dataovr%heat, "heat", has_error)
         call check_dataovr_file(cfg%ocean%dataovr%evap, "evap", has_error)
         call check_dataovr_file(cfg%ocean%dataovr%lprec, "lprec", has_error)
         call check_dataovr_file(cfg%ocean%dataovr%salt, "salt", has_error)

         ! D5: enabled but nothing driven is a wordless no-op — the
         ! house "silent no-op => fail loud" rule (cf. V2/V3 above).
         if (.not. dataovr_any_tag_set(cfg%ocean%dataovr)) then
            call logger%error("&ocean_dataovr_nml enable=.true. but no <tag>_file is set — "// &
                              "either name a forcing file or set enable=.false.")
            has_error = .true.
         end if
      end if

      ! ---- &ocean_atm_forcing_nml (PR-C2c: JRA55-do prescribed atm state) ----
      if (cfg%ocean%atm_forcing%enable) then
#ifdef RDB_NO_NETCDF
         call logger%error("&ocean_atm_forcing_nml enable=.true. requires a NetCDF "// &
                           "build (RDB_ENABLE_NETCDF=ON); the reader is compiled out here")
         has_error = .true.
#endif
         atm_code = atm_forcing_prereq_code(cfg)
         if (atm_code /= ATM_PREREQ_OK) then
            call logger%error(atm_forcing_prereq_message(atm_code))
            has_error = .true.
         end if
      end if
      ! C2c-3a wires both the provider (`atm_state_update`) and the ocean
      ! bulk-flux consumer (`bulk_forcing_apply`) into the engine -- the
      ! generic "not yet wired" refusal that used to block this group is
      ! gone (`atm_forcing_prereq_code` above already requires
      ! `&ocean_bulk_flux_nml enable` whenever `&ocean_atm_forcing_nml
      ! enable` is set). `&ocean_bulk_flux_nml enable` on its OWN has no
      ! atmospheric state to act on -- the ocean consumer reads `atm%tas`
      ! etc, which stay unallocated unless the provider is configured --
      ! so that combination is still refused, honestly, rather than
      ! silently running as a no-op.
      if (cfg%ocean%bulk_flux%enable .and. .not. cfg%ocean%atm_forcing%enable) then
         call logger%error("&ocean_bulk_flux_nml enable=.true. requires "// &
                           "&ocean_atm_forcing_nml enable=.true. -- the bulk-flux "// &
                           "kernel has no atmospheric state to act on standalone")
         has_error = .true.
      end if

      if (has_error) then
         ! NOTE: this is a terminal ROLLUP over many independent checks above
         ! (most of which log their own specific reason via `logger%error`
         ! but do not individually push to the error ring — instrumenting
         ! all ~160 of them is out of scope for this pass). This message is
         ! therefore GENERIC by construction; a caller wanting the specific
         ! reason should also look at ring index 1 (the most recently
         ! instrumented check, if any fired) rather than only index 0.
         call fail("Configuration validation failed — see errors above", ierr, &
                   OCEAN_STATUS_ERR_CONFIG_VALIDATE)
         return
      end if
      if (present(ierr)) ierr = OCEAN_STATUS_OK

   end subroutine validate_config

   ! ======================================================================
   ! &ocean_dataovr_nml predicates (PR-15).
   !
   ! Pure and non-aborting by design: a `pure` function cannot `error
   ! stop`, which is exactly what makes these testable from both sides
   ! without a subprocess/death-test harness — the house idiom, see
   ! `ice_hlim_spec_is_valid` / `zinit_dims_ok` / `data_input_dims_ok`.
   ! `validate_config` turns a `.false.` into the abort.
   ! ======================================================================

   module subroutine check_dataovr_file(e, tag, has_error)
      !! Existence check for one tag's file.  NOT `pure` — `inquire` is
      !! an I/O statement, which is exactly why this is a subroutine
      !! setting `has_error` rather than another predicate: the pure
      !! ones above stay pure and independently testable, and this
      !! keeps the filesystem dependency in one visible place.
      !! Silent for a tag that names no file.
      type(dataovr_entry_config_t), intent(in) :: e
      character(len=*), intent(in) :: tag
      logical, intent(inout) :: has_error
      logical :: exists

      if (len_trim(e%file) == 0) return
      inquire (file=trim(e%file), exist=exists)
      if (.not. exists) then
         call logger%error("&ocean_dataovr_nml "//tag//"_file not found: '"// &
                           trim(e%file)//"'")
         has_error = .true.
      end if
   end subroutine check_dataovr_file

   pure module function dataovr_entry_is_valid(e) result(ok)
      !! `.true.` unless a tag names a file without naming the variable
      !! inside it.  A blank `file` (tag not file-driven) is valid, and
      !! a `var` with no `file` is harmless — the tag simply never
      !! registers — so only the one combination is rejected.
      type(dataovr_entry_config_t), intent(in) :: e
      logical :: ok
      ok = .not. (len_trim(e%file) > 0 .and. len_trim(e%var) == 0)
   end function dataovr_entry_is_valid

   pure module function dataovr_time_is_valid(cfg) result(ok)
      !! `.true.` unless `time_mode='cyclic'` was requested without a
      !! positive `cycle_period`.  Checked only for the cyclic mode;
      !! `cycle_period` is ignored by `linear`/`static`.
      type(ocean_dataovr_config_t), intent(in) :: cfg
      logical :: ok
      ok = .true.
      if (trim(cfg%time_mode) == "cyclic") ok = (cfg%cycle_period > 0.0_wp)
   end function dataovr_time_is_valid

   pure module function dataovr_freshwater_needs_components(cfg, use_components) result(bad)
      !! `.true.` when `evap`/`lprec` are file-driven but the
      !! surface-flux component set they write into is switched off.
      !! Note the polarity: this reports the PROBLEM, not validity —
      !! named for how it reads at the call site.
      type(ocean_dataovr_config_t), intent(in) :: cfg
      logical, intent(in) :: use_components
      logical :: bad
      bad = (len_trim(cfg%evap%file) > 0 .or. len_trim(cfg%lprec%file) > 0) &
            .and. (.not. use_components)
   end function dataovr_freshwater_needs_components

   pure module function atm_forcing_prereq_code(cfg) result(code)
      !! First failing prerequisite of `&ocean_atm_forcing_nml enable` as an
      !! `ATM_PREREQ_*` code (`ATM_PREREQ_OK` when all hold). Independent of
      !! the NetCDF build and of the engine-wiring refusal, so each reason
      !! is testable on every toolchain.
      type(config_t), intent(in) :: cfg
      integer :: code

      code = ATM_PREREQ_OK
      if (len_trim(cfg%start_date) == 0) then
         code = ATM_PREREQ_START_DATE
      else if (.not. cfg%ocean%forcing%enable_components) then
         code = ATM_PREREQ_COMPONENTS
      else if (trim(cfg%ocean%grid%grid_config) == "cartesian") then
         code = ATM_PREREQ_CARTESIAN
      else if (len_trim(cfg%ocean%atm_forcing%filelist_tas) == 0 .or. &
               len_trim(cfg%ocean%atm_forcing%filelist_huss) == 0 .or. &
               len_trim(cfg%ocean%atm_forcing%filelist_uas) == 0 .or. &
               len_trim(cfg%ocean%atm_forcing%filelist_vas) == 0 .or. &
               len_trim(cfg%ocean%atm_forcing%filelist_psl) == 0 .or. &
               len_trim(cfg%ocean%atm_forcing%filelist_rlds) == 0 .or. &
               len_trim(cfg%ocean%atm_forcing%filelist_rsds) == 0 .or. &
               len_trim(cfg%ocean%atm_forcing%filelist_prra) == 0 .or. &
               len_trim(cfg%ocean%atm_forcing%filelist_prsn) == 0) then
         code = ATM_PREREQ_FILELIST
      else if (cfg%ocean%thermo%sw_pen_frac > 0.0_wp .and. &
               trim(cfg%ocean%thermo%sw_source) /= "q_sw") then
         code = ATM_PREREQ_SW_SOURCE
      else if (len_trim(cfg%ocean%dataovr%tau_x%file) > 0 .or. &
               len_trim(cfg%ocean%dataovr%tau_y%file) > 0 .or. &
               len_trim(cfg%ocean%dataovr%heat%file) > 0 .or. &
               len_trim(cfg%ocean%dataovr%evap%file) > 0 .or. &
               len_trim(cfg%ocean%dataovr%lprec%file) > 0) then
         code = ATM_PREREQ_DATAOVR
      else if (.not. cfg%ocean%bulk_flux%enable) then
         code = ATM_PREREQ_BULK_FLUX
      else if (.not. bulk_flux_config_n_iter_ok(cfg%ocean%bulk_flux)) then
         code = ATM_PREREQ_N_ITER
      else if (cfg%ocean%atm_forcing%orb_year /= 1990) then
         code = ATM_PREREQ_ORB_YEAR
      end if
   end function atm_forcing_prereq_code

   pure module function atm_forcing_prereq_message(code) result(msg)
      !! The fail-loud message for an `ATM_PREREQ_*` code.
      integer, intent(in) :: code
      character(len=:), allocatable :: msg

      select case (code)
      case (ATM_PREREQ_START_DATE)
         msg = "&ocean_atm_forcing_nml enable=.true. requires &time_nml start_date "// &
               "(absolute file-time axes need a model-clock epoch)"
      case (ATM_PREREQ_COMPONENTS)
         msg = "&ocean_atm_forcing_nml enable=.true. requires &ocean_forcing_nml "// &
               "enable_components=.true. (the bulk-flux outputs land on the component set)"
      case (ATM_PREREQ_CARTESIAN)
         msg = "&ocean_atm_forcing_nml enable=.true. requires a non-cartesian "// &
               "&ocean_grid_nml grid_config (spherical/supergrid/tripolar) -- a "// &
               "cartesian grid has no geographic lon/lat to regrid the atmosphere onto"
      case (ATM_PREREQ_FILELIST)
         msg = "&ocean_atm_forcing_nml enable=.true. requires all nine filelist_<field> "// &
               "paths (tas/huss/uas/vas/psl/rlds/rsds/prra/prsn) -- a blank one is a "// &
               "fail-loud setup, never a silently-skipped field"
      case (ATM_PREREQ_SW_SOURCE)
         msg = "&ocean_atm_forcing_nml enable=.true. with &ocean_thermo_nml "// &
               "sw_pen_frac > 0 requires sw_source='q_sw' (the atmosphere writes "// &
               "q_sw, not net_heat)"
      case (ATM_PREREQ_DATAOVR)
         msg = "&ocean_atm_forcing_nml enable=.true. conflicts with &ocean_dataovr_nml "// &
               "tau_x/tau_y/heat/evap/lprec -- both would write the same surface-flux "// &
               "components (a salt_file MAY coexist: SSS restoring is independent)"
      case (ATM_PREREQ_BULK_FLUX)
         msg = "&ocean_atm_forcing_nml enable=.true. requires &ocean_bulk_flux_nml "// &
               "enable=.true. (set it explicitly)"
      case (ATM_PREREQ_N_ITER)
         msg = "&ocean_bulk_flux_nml n_iter must be >= 1"
      case (ATM_PREREQ_ORB_YEAR)
         msg = "&ocean_atm_forcing_nml orb_year is not wired yet (the coszen kernel "// &
               "uses fixed Spencer (1971) coefficients); leave it at the default 1990"
      case default
         msg = ""
      end select
   end function atm_forcing_prereq_message

   pure module function dataovr_any_tag_set(cfg) result(any_set)
      !! `.true.` when at least one tag names a file, i.e. enabling the
      !! group would actually do something.
      type(ocean_dataovr_config_t), intent(in) :: cfg
      logical :: any_set
      any_set = len_trim(cfg%tau_x%file) > 0 .or. len_trim(cfg%tau_y%file) > 0 .or. &
                len_trim(cfg%heat%file) > 0 .or. len_trim(cfg%evap%file) > 0 .or. &
                len_trim(cfg%lprec%file) > 0 .or. len_trim(cfg%salt%file) > 0
   end function dataovr_any_tag_set

   pure module function ice_hlim_count(hlim) result(n)
      !! Length of the LEADING run of non-sentinel (>= 0) entries in
      !! `&ocean_ice_nml hlim` (PR-58). Contiguity matters here (unlike
      !! `apply_layer_rho_init`'s bare `count(>= 0)`): a value AFTER the
      !! first sentinel is silently dropped by this count, so
      !! `ice_hlim_spec_is_valid` checks contiguity explicitly rather than
      !! trusting the count alone.
      real(wp), intent(in) :: hlim(:)
      integer :: n
      integer :: k

      n = 0
      do k = 1, size(hlim)
         if (hlim(k) < 0.0_wp) exit
         n = n + 1
      end do
   end function ice_hlim_count

   pure module subroutine ice_hlim_spec_is_valid(hlim, ncat, ok, reason)
      !! Shape/monotonicity predicate for a supplied `&ocean_ice_nml hlim`
      !! list (PR-58). A SUBROUTINE, not a function — a pure FUNCTION may
      !! not carry an `intent(out)` dummy (`reason`); pure SUBROUTINES
      !! can. `n = ice_hlim_count(hlim)`; `ok` iff ALL of:
      !!   (a) n >= 2                          -- one edge gives no width
      !!                                          to extrapolate (SIS2
      !!                                          divergence D1: Roundabout
      !!                                          fails loud here instead
      !!                                          of SIS2's silent
      !!                                          fallback to the default
      !!                                          table).
      !!   (b) n <= ncat + 1                   -- more edges than the ITD
      !!                                          has bins.
      !!   (c) all(hlim(n+1:) < 0)             -- contiguity: no value
      !!                                          after the first
      !!                                          sentinel.
      !!   (d) hlim(1) > 0                     -- mh_lim(1) > 0 is a live
      !!                                          gate (`rdb_ice_itd`
      !!                                          `ice_adjust_categories`).
      !!   (e) hlim(k+1) > hlim(k), k = 1..n-1  -- strictly increasing.
      !! Checking (d)+(e) on the SUPPLIED list is sufficient for the FULL
      !! `h_lim(1..ncat+1)`: the constant-width extrapolation preserves
      !! both properties by induction
      !! (`h_lim(k)-h_lim(k-1) = h_lim(k-1)-h_lim(k-2)`), so no redundant
      !! post-extrapolation check is needed.
      real(wp), intent(in) :: hlim(:)
      integer, intent(in) :: ncat
      logical, intent(out) :: ok
      character(len=:), allocatable, intent(out) :: reason
      integer :: n, k

      n = ice_hlim_count(hlim)
      ok = .false.
      reason = ""

      if (n < 2) then
         reason = "hlim must supply >= 2 entries (a single edge has no "// &
                  "width to extrapolate); got "//to_string(n)
         return
      end if
      if (n > ncat + 1) then
         reason = "hlim supplies "//to_string(n)//" entries, more than "// &
                  "ncat+1 ("//to_string(ncat + 1)//")"
         return
      end if
      if (.not. all(hlim(n + 1:) < 0.0_wp)) then
         reason = "hlim is not a contiguous leading run: a value follows "// &
                  "the first sentinel (entry "//to_string(n + 1)//")"
         return
      end if
      if (hlim(1) <= 0.0_wp) then
         reason = "hlim(1) must be > 0"
         return
      end if
      do k = 1, n - 1
         if (.not. (hlim(k + 1) > hlim(k))) then
            reason = "hlim must be strictly increasing (violated at entry "// &
                     to_string(k + 1)//")"
            return
         end if
      end do

      ok = .true.
      reason = "ok"
   end subroutine ice_hlim_spec_is_valid

   pure module function resolve_bt_halo(requested, compute_size, exclusion_active) result(width)
      !! Resolve the `&ocean_bt_nml bt_halo` sentinel to a concrete march-in
      !! width.
      !!   requested == BT_HALO_AUTO_SENTINEL (-1, the default) => AUTO: 0.
      !!     The march-in used to switch itself on (BT_HALO_AUTO_WIDTH) on
      !!     every compatible multi-rank run.  It is not bit-reproducible
      !!     against the serial run over variable bathymetry or with open
      !!     boundaries (`test_ocean_decomp_bitid_mpi`), so a default
      !!     multi-rank run must not pick it: AUTO is off, serial or not.
      !!   requested >= 0 (user set it explicitly) => returned UNCHANGED; the
      !!     fail-loud exclusion checks in validate_config police an explicit
      !!     bt_halo > 0 against an incompatible feature (user asked for the
      !!     impossible), so the auto-resolution never overrides an explicit 0
      !!     or an explicit width.
      !! `compute_size` / `exclusion_active` no longer change the answer; they
      !! stay so the call site (and its exclusion log) keeps its shape for the
      !! day the march-in is exact again.
      integer, intent(in) :: requested
         !! The namelist value: BT_HALO_AUTO_SENTINEL (-1) for auto, else >= 0.
      integer, intent(in) :: compute_size
         !! Number of compute ranks (1 = serial).  Unused by AUTO (see above).
      logical, intent(in) :: exclusion_active
         !! .true. iff any march-in exclusion feature is on (see
         !! `bt_halo_auto_exclusion`).  Unused by AUTO (see above).
      integer :: width
      width = 0
      if (requested >= 0) width = requested
      if (compute_size < 0 .and. exclusion_active) width = 0   ! never taken; keeps the dummies referenced
   end function resolve_bt_halo

   pure module subroutine bt_halo_auto_exclusion(cfg, excluded, reason)
      !! Single source of truth for the BT march-in exclusion set (mirrors the
      !! explicit-`bt_halo > 0` fail-loud checks in `validate_config`).  Reports
      !! whether ANY exclusion is active and names the first one (for the AUTO
      !! resolution log).  When none is active, `reason` is "".
      type(config_t), intent(in) :: cfg
      logical, intent(out) :: excluded
      character(len=:), allocatable, intent(out) :: reason
      reason = ""
      if (cfg%ocean%wetdry%enable) then
         reason = "wetdry enable"
      else if (cfg%ocean%bt%use_cont_type) then
         reason = "use_cont_type"
      else if (cfg%ocean%bt%upstream_h_face) then
         reason = "upstream_h_face"
      else if (cfg%ocean%tides%enable) then
         reason = "tides enable"
      else if (cfg%ocean%psurf%enable) then
         ! Surface-pressure loading rides the SAME `eta_forcing` seam as
         ! the body tide, and the wide BT clone carries no copy of it:
         ! `run_stage_split` `error stop`s on `bt_halo > 0` with an
         ! `eta_forcing` actual.  Must exclude alongside tides.
         reason = "psurf enable"
      else if (cfg%ocean%porous%enable) then
         reason = "porous enable"
      else if (cfg%ocean%cavity_dyn%enable) then
         ! The wide BT clone rebuilds `metrics_w` from the grid formula and
         ! carries no `z_draft`, so its reference depth would be the BED —
         ! it would solve a 1000 m ocean where the cavity has 500 m of
         ! water under 500 m of ice.  `validate_config` refuses an EXPLICIT
         ! `bt_halo > 0`; AUTO must resolve to 0 here or a multi-rank
         ! cavity run manufactures a width the user never asked for and
         ! then trips that abort.
         reason = "cavity_dyn enable"
      else if (trim(cfg%ocean%grid%grid_config) == "supergrid" .or. &
               trim(cfg%ocean%grid%grid_config) == "tripolar") then
         reason = "grid_config='"//trim(cfg%ocean%grid%grid_config)//"'"
      end if
      excluded = len_trim(reason) > 0
   end subroutine bt_halo_auto_exclusion

   pure module function p_top_has_producer(cfg) result(has)
      !! Is there anything in this configuration that WRITES
      !! `multilayer_state_t%p_top`?
      !!
      !! `p_top = metrics%p_ice_ref + sf%p_surf` has exactly two
      !! producers, and both halves count:
      !!
      !!   * `&ocean_psurf_nml enable` — the atmospheric surface-pressure
      !!     seam, which fills `sf%p_surf` and refreshes `p_top` once per
      !!     outer step in `ocean_dyn_step_split`;
      !!   * `&ocean_cavity_dyn_nml enable` — the STATIC ice-shelf load
      !!     `p_ice_ref = rho_ref*GRAVITY*z_draft`, assembled into
      !!     `p_top` by `configure_ocean_cavity`.  Static is not the same
      !!     as absent: the draft never changes, so the configure-time
      !!     assembly is the final value and there is nothing to refresh.
      !!
      !! Used by the `&ocean_pgf_nml p_top_in_bc` inert-knob warning,
      !! which must name both — a cavity run is precisely the case where
      !! `p_top_in_bc` is not merely live but REQUIRED (for a draft that
      !! varies), so warning that it is inert there was telling the user
      !! the opposite of the truth.
      type(config_t), intent(in) :: cfg
      logical :: has
      has = cfg%ocean%psurf%enable .or. cfg%ocean%cavity_dyn%enable
   end function p_top_has_producer

   pure module function bbl_glue_is_effective(cfg) result(on)
      !! Will `&ocean_vdiff_nml bbl_glue` actually be ON after setup?
      !!
      !! `vdiff_bbl_configure` (`rdb_ocean_vdiff.F90`) turns a requested
      !! glue OFF when its effective drag coefficient or BBL thickness is
      !! not positive: quadratic `cd <= 0`, linear `bg_vel <= 0` (no
      !! `cd = r·hbbl/bg_vel`) or `r <= 0`, or `hbbl` and its `hbbl_visc`
      !! fallback both `<= 0`.  This is the same predicate on `cfg`, so a
      !! configure check that only holds while the glue is on (the
      !! `implicit_drag` + `hbbl > 0` exception) cannot accept a namelist
      !! the setup then runs without it.  Keep the two in step.
      type(config_t), intent(in) :: cfg
      logical :: on
      real(wp) :: hbbl, cd_eff
      on = .false.
      if (.not. cfg%ocean%vdiff%bbl_glue) return
      hbbl = cfg%ocean%bdrag%hbbl
      if (hbbl <= 0.0_wp) hbbl = cfg%ocean%vdiff%hbbl_visc
      if (trim(adjustl(cfg%ocean%bdrag%form)) == "linear") then
         cd_eff = 0.0_wp
         if (cfg%ocean%bdrag%bg_vel > 0.0_wp) then
            cd_eff = cfg%ocean%bdrag%r*hbbl/cfg%ocean%bdrag%bg_vel
         end if
      else
         cd_eff = cfg%ocean%bdrag%cd
      end if
      on = cd_eff > 0.0_wp .and. hbbl > 0.0_wp
   end function bbl_glue_is_effective

   pure module function substep_drag_ignores_bdrag_form(cfg) result(ignores)
      !! Is `&ocean_bt_nml substep_drag` blind to the configured bottom
      !! drag?
      !!
      !! `compute_bt_rem` builds the barotropic substep damping
      !! `Htot/(Htot + r·hbbl·dt_inner)` from the LINEAR coefficient
      !! `&ocean_bdrag_nml r` alone.  Under any `form` other than
      !! `"linear"` the slow bottom drag is not that operator, so the
      !! knob either does nothing at all (`r = 0`, the default — the
      !! factor is identically 1) or damps the barotropic mode with a
      !! linear drag the slow step never applies.  Neither is what the
      !! user asked for; the configure WARNING that names both knobs
      !! reads this predicate.  A warning, not a refusal: the knob is
      !! harmless-if-useless in the default case, and a deliberate
      !! BT-only linear sponge is a legitimate (if unusual) request.
      type(config_t), intent(in) :: cfg
      logical :: ignores
      ignores = cfg%ocean%bt%substep_drag .and. &
                trim(adjustl(cfg%ocean%bdrag%form)) /= "linear"
   end function substep_drag_ignores_bdrag_form

   pure module function ocean_bt_forcing_visc_rem_on(cfg) result(on)
      !! PR-3 (D1): is the MOM6 `wt_u` BT-forcing weight on, either
      !! directly (`forcing_visc_rem`) or via `visc_rem_chain`?
      type(config_t), intent(in) :: cfg
      logical :: on
      on = cfg%ocean%bt%visc_rem_chain .or. cfg%ocean%bt%forcing_visc_rem
   end function ocean_bt_forcing_visc_rem_on

   pure module function ocean_bt_renorm_visc_rem_on(cfg) result(on)
      !! PR-3 (D1): is the MOM6 continuity `u_cor = u + du*visc_rem`
      !! transport-matching renormaliser on, either directly
      !! (`renorm_visc_rem`) or via `visc_rem_chain`?
      type(config_t), intent(in) :: cfg
      logical :: on
      on = cfg%ocean%bt%visc_rem_chain .or. cfg%ocean%bt%renorm_visc_rem
   end function ocean_bt_renorm_visc_rem_on

   pure module function ocean_bt_rem_from_visc_rem_on(cfg) result(on)
      !! PR-3 (D1): is `bt_rem` built from `av_rem` (MOM6's barotropic
      !! viscous-remnant depth mean) on, either directly
      !! (`bt_rem_from_visc_rem`) or via `visc_rem_chain`?
      type(config_t), intent(in) :: cfg
      logical :: on
      on = cfg%ocean%bt%visc_rem_chain .or. cfg%ocean%bt%bt_rem_from_visc_rem
   end function ocean_bt_rem_from_visc_rem_on

   pure module function ocean_bt_visc_rem_producer_on(cfg) result(on)
      !! D1 follow-up: is the visc_rem PRODUCER needed, independent of
      !! the (retired) weighted BT-correction fold?  `.true.` whenever
      !! ANY real consumer is on — `forcing_visc_rem`/`renorm_visc_rem`/
      !! `bt_rem_from_visc_rem` (each already `.or.`-ed with
      !! `visc_rem_chain` by their own helper) — or the legacy
      !! `correction_visc_rem` field itself, so a test that constructs
      !! `cfg` directly and sets that field alone (bypassing the nml
      !! retirement check) still gets a live producer.  This is what
      !! `configure_ocean_bt` wires into `bt_work%bt_visc_rem_producer`,
      !! which `vmix_apply_in_stage`'s `do_remnant` reads — NOT
      !! `bt_work%bt_correction_visc_rem`, which now drives ONLY the
      !! weighted-fold dispatch in `apply_bt_correction` (and is never
      !! set by `visc_rem_chain`).
      type(config_t), intent(in) :: cfg
      logical :: on
      on = cfg%ocean%bt%correction_visc_rem .or. ocean_bt_forcing_visc_rem_on(cfg) &
           .or. ocean_bt_renorm_visc_rem_on(cfg) .or. ocean_bt_rem_from_visc_rem_on(cfg)
   end function ocean_bt_visc_rem_producer_on

   pure module function cavity_draft_is_uniform(cfg) result(uniform)
      !! Is the configured ice-shelf draft UNIFORM over the whole array?
      !!
      !! The one predicate, so the two rules that turn on it cannot drift
      !! apart: the `&ocean_pgf_nml p_top_in_bc` REFUSAL (a load with a
      !! gradient must have a consumer — a uniform load is bit-identically
      !! inert in the FV_MOM6 top BC, which is the theorem in
      !! `compute_fv_mom6_impl`'s docstring), and the `z_fixed` x cavity
      !! staircase WARNING (only a uniform draft is validated on a
      !! quasi-geopotential coordinate).
      !!
      !! `draft_config = "none"` is uniform because it is identically
      !! zero.  `"flat"` is uniform only when NO box bound clips it —
      !! a clipped flat draft has a calving front, which is a step, and a
      !! step is the largest gradient in the domain.  `"linear"` and
      !! `"file"` are never assumed uniform: this is a NAMELIST-level
      !! predicate and cannot see the filled array (the configure-time
      !! twin in `configure_ocean_cavity` tests `maxval /= minval` on the
      !! field itself, which is the stricter check and runs later).
      type(config_t), intent(in) :: cfg
      logical :: uniform
      character(len=:), allocatable :: dcfg
      dcfg = trim(adjustl(cfg%ocean%cavity_dyn%draft_config))
      uniform = (dcfg == "none")
      if (dcfg == "flat") then
         uniform = abs(cfg%ocean%cavity_dyn%draft_x0) >= 1.0e29_wp .and. &
                   abs(cfg%ocean%cavity_dyn%draft_x1) >= 1.0e29_wp .and. &
                   abs(cfg%ocean%cavity_dyn%draft_y0) >= 1.0e29_wp .and. &
                   abs(cfg%ocean%cavity_dyn%draft_y1) >= 1.0e29_wp
      end if
   end function cavity_draft_is_uniform

   pure module function zfixed_cavity_nu_h_below_envelope(cfg) result(below)
      !! Is the constant harmonic viscosity below the `z_fixed` x cavity
      !! lower envelope `ZFIXED_CAVITY_NU_H_MIN`?  The predicate behind
      !! `validate_config`'s WARNING (never a refusal — the vcoord
      !! stability matrix runs this combination inviscid on purpose).
      !!
      !! The envelope is measured, not derived.  ISOMIP+ Ocean0 idealised
      !! at 2 km, `z_fixed` + `zfixed_closed_faces`, melt off, 30 days
      !! (`cavity_rest_growth_diagnosis.md` §Q.3/§Q.0 item 7):
      !! `nu_h = 0` grows EXPONENTIALLY at `0.180 /day` (d15-30) and
      !! accelerating (0.142 -> 0.220 /day, d15-20 -> d25-30;
      !! `En(30 d) = 1.44E-06`), `nu_h = 2` decelerates (0.092 -> 0.068
      !! /day, d5-10 -> d25-30; `1.12E-07`), the protocol's 6 gives
      !! `5.63E-08` (V100, v0.1.0 defaults: bebt = 0.1,
      !! renorm_consistent_flux, I1').  The
      !! inviscid mode is arrested somewhere in `0 < nu_h < 2`; the
      !! bracket below 2 was not refined, so 2 is the bound of record.
      !!
      !! `.false.` whenever the cavity is off or the coordinate is not
      !! `z_fixed`, so it cannot fire on any other configuration.
      use rdb_vcoord, only: parse_vcoord_type
      use rdb_constants, only: VCOORD_EULERIAN_Z, VCOORD_Z_FIXED
      type(config_t), intent(in) :: cfg
      logical :: below
      below = .false.
      if (.not. cfg%ocean%cavity_dyn%enable) return
      if (parse_vcoord_type(cfg%vcoord_type, default_code=VCOORD_EULERIAN_Z) &
          /= VCOORD_Z_FIXED) return
      below = cfg%ocean%hvisc%nu_h < ZFIXED_CAVITY_NU_H_MIN
   end function zfixed_cavity_nu_h_below_envelope

   module subroutine warn_unknown_bc(bc_str, param_name)
      !! Warn if a BC string does not match any known type
      character(len=*), intent(in) :: bc_str
      character(len=*), intent(in) :: param_name

      select case (trim(bc_str))
      case ("wall", "open", "tidal", "nested", "inflow", &
            "discharge", "clamped", "sponge", "chapman")
         continue
      case default
         call logger%warning("Unknown boundary type '"//trim(bc_str)// &
                             "' for "//param_name//", defaulting to wall")
      end select

   end subroutine warn_unknown_bc

   pure logical module function forcing_components_need_thermo(enable_components, &
                                                               enable_thermodynamics) result(bad)
      !! `&ocean_forcing_nml enable_components=.true.` with
      !! `&ocean_thermo_nml enable_thermodynamics=.false.` is a
      !! configuration that validates but does nothing (the component
      !! set only feeds the assembler->apply_tracers thermo path) —
      !! named predicate per the `lateral_closure_is_implemented` idiom
      !! so `validate_config` reads as one line.
      logical, intent(in) :: enable_components, enable_thermodynamics
      bad = enable_components .and. .not. enable_thermodynamics
   end function forcing_components_need_thermo

   pure logical module function diag_density_levels_ok(vgrid, diags, n_rho_levels, rho_levels) result(ok)
      !! .true. iff a density-coordinate diagnostic selection (global
      !! `vgrid='density'`, or a per-diagnostic `:density`/`:rho` attribute
      !! anywhere in `diags`) has a usable `rho_levels` axis: at least one
      !! entry, strictly increasing (`invert_density_targets` assumes a
      !! monotone light->dense target list), and within the declared
      !! `MAX_OCEAN_DIAG_Z_LEVELS` array bound.  Unlike sigma/z*, density
      !! bins have no auto-fill, so an unset axis must abort rather than
      !! silently size to zero.  When density is not requested at all this
      !! is unconditionally `.true.` (no constraint).
      character(len=*), intent(in) :: vgrid
      character(len=*), intent(in) :: diags
      integer, intent(in) :: n_rho_levels
      real(wp), intent(in) :: rho_levels(:)
      logical :: density_requested
      integer :: k

      density_requested = (trim(vgrid) == "density") .or. &
                          (index(diags, ":density") > 0) .or. &
                          (index(diags, ":rho") > 0)
      if (.not. density_requested) then
         ok = .true.
         return
      end if

      ok = .true.
      if (n_rho_levels <= 0) then
         ok = .false.
         return
      end if
      if (n_rho_levels > size(rho_levels)) then
         ok = .false.
         return
      end if
      do k = 2, n_rho_levels
         if (rho_levels(k) <= rho_levels(k - 1)) then
            ok = .false.
            return
         end if
      end do
   end function diag_density_levels_ok

   pure logical module function sponge_source_is_implemented(tag) result(ok)
      !! .true. iff `&ocean_sponge_nml damp_source` names a source that is
      !! actually filled. `"file"` is a recognised name (PR-23b, needs the
      !! PR-14 reader) but has no kernel yet — a source with no filler
      !! must abort, not silently build an all-zero map (the
      !! `lateral_closure_is_implemented` idiom).
      !!
      !! Kept under its original name because it is the DAMP axis; the
      !! TARGET axis has its own predicate below. They were one function
      !! until `"linear_z"` arrived, at which point a shared list would
      !! have accepted `damp_source="linear_z"` — a spelling with no
      !! meaning on that axis — as valid.
      character(len=*), intent(in) :: tag
      select case (trim(tag))
      case ("band")
         ok = .true.
      case default
         ok = .false.
      end select
   end function sponge_source_is_implemented

   pure logical module function sponge_target_is_implemented(tag) result(ok)
      !! .true. iff `&ocean_sponge_nml target_source` names a reference
      !! state that is actually filled: `"ic"` (snapshot of the seeded
      !! initial condition) or `"linear_z"` (analytic affine geopotential
      !! profile, re-evaluated on the live layer geometry). `"file"` is
      !! recognised but has no reader yet (PR-23b).
      character(len=*), intent(in) :: tag
      select case (trim(tag))
      case ("ic", "linear_z")
         ok = .true.
      case default
         ok = .false.
      end select
   end function sponge_target_is_implemented

   pure logical module function sponge_ramp_is_valid(tag) result(ok)
      !! .true. iff `&ocean_sponge_nml ramp` names a band shape the
      !! `damp_source="band"` filler implements.
      character(len=*), intent(in) :: tag
      select case (trim(tag))
      case ("cosine", "linear")
         ok = .true.
      case default
         ok = .false.
      end select
   end function sponge_ramp_is_valid

end submodule rdb_config_validate
