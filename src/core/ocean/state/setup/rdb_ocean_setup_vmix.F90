submodule(rdb_ocean_setup) rdb_ocean_setup_vmix
   !! R2 split: vmix-group configure_ocean_* bodies, moved verbatim
   !! from rdb_ocean_setup.F90 (host-only configure code; see
   !! src/core/ocean/state/rdb_ocean_setup.F90 for the public interface).
#ifndef RDB_NO_NETCDF
   use rdb_bathymetry, only: load_static_2d_field
#endif
   use rdb_ocean_tidal_mixing, only: tidal_mixing_fold_static_fields
   use rdb_ocean_status, only: OCEAN_STATUS_ERR_IO
   implicit none

contains

   module subroutine configure_ocean_drag(cfg, ocean_state, compute_rank, ierr)
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
      ! T/S convention + liquidus form: same site and reasons as the
      ! liquidus set above (flat POD, host assignment before any copy and
      ! before `enter_data`).  Invalid strings are fatal in `validate_config`.
      call eos_apply_ts_convention(ocean_state%eos, &
                                   parse_ts_convention(cfg%ocean%eos%ts_convention), &
                                   parse_tfreeze_form(cfg%ocean%eos%tfreeze_form), &
                                   cfg%ocean%eos%tfreeze_sat_frac)
      ocean_state%eos%c_p = cfg%ocean%eos%c_p
      if (compute_rank == 0) then
         call logger%info("EOS variant:      "//trim(cfg%ocean%eos%eos))
         if (ocean_state%eos%ts_convention == TS_CONS_ABS) then
            call logger%info("T/S convention:   cons_abs (T = Conservative Temperature, "// &
                             "S = Absolute Salinity)")
         end if
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
            ! Bed-only mode drags the face's first LIVE layer `k_bot_u/v`.
            ! Under z_fixed that index is filled at configure
            ! (`configure_ocean_k_bot`) and skips the static bed fillers, so
            ! the bed-only mode reaches water there.  `k_bot ≡ 1` under
            ! zstar_full (its bed-side layers vanish DYNAMICALLY, no static
            ! pattern) and, until `configure_ocean_k_bot` learns it, under
            ! zstar (static z_fixed fillers): k = 1 is an inert filler in
            ! every column shallower than the deepest nominal interface and
            ! almost the whole domain runs with NO bottom drag.  HBBL mode
            ! accumulates thickness from the bed up, skips those fillers and
            ! reaches the live bottom layer.
            if ((cfg%ocean%bdrag%cd > 0.0_wp .or. cfg%ocean%bdrag%r > 0.0_wp) .and. &
                (parse_ocean_vcoord_type(cfg%vcoord_type) == VCOORD_ZSTAR .or. &
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

   module subroutine configure_ocean_hdiff(cfg, ocean_state, compute_rank)
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

   module subroutine configure_ocean_vmix(cfg, ocean_state, compute_rank, ierr)
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
      ocean_state%vdiff%hvel_harmonic = cfg%ocean%vdiff%hvel_harmonic
      ! MOM6 bottomdraglaw coupling parity: kv_bbl botfn glue + piston bed
      ! drag (PGF_BUG.md §9), with MOM6's per-face `set_viscous_BBL`
      ! (`vdiff_set_viscous_bbl`) built from the bottom-drag configuration:
      ! the drag law, `cd`/`r`, `bg_vel`, `bbl_thick_min`, and `HBBL` =
      ! `&ocean_bdrag_nml hbbl`, falling back to `&ocean_vdiff_nml
      ! hbbl_visc` for a drag configured bed-only (`hbbl = 0`) — MOM6 has
      ! ONE HBBL and no bed-only form.  No drag configured ⇒ the latch
      ! leaves the glue OFF.  `hvel_mom6` is validated at configure.
      ocean_state%vdiff%bbl_glue = cfg%ocean%vdiff%bbl_glue
      ocean_state%vdiff%bbl_piston = cfg%ocean%vdiff%bbl_piston
      ocean_state%vdiff%hvel_upwind = cfg%ocean%vdiff%hvel_upwind
      block
         integer :: bbl_form
         real(wp) :: bbl_hbbl
         bbl_form = BBL_FORM_QUADRATIC
         if (trim(cfg%ocean%bdrag%form) == "linear") bbl_form = BBL_FORM_LINEAR
         bbl_hbbl = cfg%ocean%bdrag%hbbl
         if (bbl_hbbl <= 0.0_wp) bbl_hbbl = cfg%ocean%vdiff%hbbl_visc
         call vdiff_bbl_configure(ocean_state%vdiff, &
                                  size(ocean_state%multilayer%h_layer, 1), &
                                  size(ocean_state%multilayer%h_layer, 2), &
                                  size(ocean_state%multilayer%h_layer, 3), &
                                  bbl_form, cfg%ocean%bdrag%cd, cfg%ocean%bdrag%r, &
                                  bbl_hbbl, cfg%ocean%bdrag%bg_vel, &
                                  cfg%ocean%bdrag%bbl_thick_min, &
                                  cfg%ocean%kshear%enable, ocean_state%eos%rho0, &
                                  cfg%ocean%vmix%pp81_nu_bg)
      end block
      if (compute_rank == 0 .and. cfg%ocean%vdiff%hvel_mom6 &
          .and. .not. cfg%ocean%vdiff%hvel_upwind) then
         call logger%info("vdiff: hvel near-bed upwind blend OFF (pure harmonic hvel)")
      end if
      if (compute_rank == 0 .and. cfg%ocean%vdiff%hvel_mom6) then
         if (cfg%ocean%vdiff%hvel_harmonic) then
            call logger%info("vdiff: MOM6 hvel, HARMONIC_VISC=True branch (harmonic + "// &
                             "near-bed upwind blend, arithmetic h_shear) ON")
         else
            call logger%info("vdiff: MOM6 hvel, HARMONIC_VISC=False branch (arithmetic, "// &
                             "z_clear near-bed harmonic blend, arithmetic h_shear) ON")
         end if
      end if
      if (compute_rank == 0 .and. ocean_state%vdiff%bbl_glue) then
         call logger%info("vdiff: MOM6 bottomdraglaw BBL glue ON (per-face "// &
                          "set_viscous_BBL kv_bbl/bbl_thick, "//trim(cfg%ocean%bdrag%form)// &
                          " law; botfn interfaces + piston bed row replace the "// &
                          "explicit bed-drag apply)")
      else if (compute_rank == 0 .and. cfg%ocean%vdiff%bbl_glue) then
         call logger%info("vdiff: bbl_glue requested but no bottom drag is configured "// &
                          "(&ocean_bdrag_nml cd/r = 0, or linear with bg_vel = 0) — glue OFF")
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

   module subroutine configure_ocean_tracers(cfg, ocean_state, compute_rank)
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

   module subroutine configure_ocean_epbl(cfg, ocean_state, grid, compute_rank, ierr)
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
                      trim(ecfg%combine)//"' (add/max/max_shear)", ierr, OCEAN_STATUS_ERR_SETUP)
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
         if (epbl%combine_mode == EPBL_COMBINE_MAX_SHEAR .and. cfg%ocean%vmix%use_kpp) then
            call fail("ocean_epbl_nml: combine='max_shear' is refused with use_kpp=.true. "// &
                      "(the max-with-kappa-shear combine is defined against the "// &
                      "ePBL-only surface scheme); set use_kpp=.false.", ierr, OCEAN_STATUS_ERR_SETUP)
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

   module subroutine configure_ocean_foxkemper(cfg, ocean_state, compute_rank, ierr)
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

   module subroutine configure_ocean_kappa_shear(cfg, ocean_state, grid, compute_rank, ierr)
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

   module subroutine configure_ocean_tidal_mixing(cfg, ocean_state, grid, compute_rank, ierr)
      !! Copy the `&ocean_tidal_mixing_nml` knobs onto the tidal-mixing
      !! slot — every field, end to end (no dead-config gaps).  Seeds the
      !! prescribed bottom energy field `e_in` (v1 uniform path), copies
      !! the shared EOS hookup for the N^2 buoyancy derivatives, and
      !! validates.  No mutual exclusion: tidal mixing is an interior
      !! closure that coexists with KPP / EPBL / PP81 / kappa-shear.
      type(config_t), intent(in) :: cfg
      type(ocean_state_t), intent(inout) :: ocean_state
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: compute_rank
      integer, intent(out), optional :: ierr
         !! Non-zero on a tidal-mixing configuration conflict when present;
         !! absent behaves as today (`error stop`).
      integer :: local_ierr

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

         ! OM3 refinements (all default to the legacy behaviour).
         tm%n2_bot_rough_avg = tcfg%n2_bot_rough_avg
         tm%limiter_reinject = (trim(tcfg%energy_limiter) == "reinject")
         tm%bed_lump = (trim(tcfg%bed_share) == "lump")

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

         ! Static 2-D tide fields.  "uniform" spreads the scalar; "file"
         ! reads the time-less 2-D variable (all ranks read their band).
         tm%tideamp(:, :) = tm%utide
         tm%h2(:, :) = tm%h2_rough
         tm%use_fields = (trim(tcfg%tideamp_source) == "file") .or. &
                         (trim(tcfg%h2_source) == "file")
         if (tm%use_fields .and. .not. tm%e_compute) then
            call fail("ocean_tidal_mixing_nml: tideamp_source/h2_source = 'file' "// &
                      "requires e_compute=.true. (the fields feed the "// &
                      "state-dependent energy input only)", ierr, OCEAN_STATUS_ERR_SETUP)
            return
         end if
         if (trim(tcfg%tideamp_source) == "file") then
            call load_tide_field(tcfg%tideamp_file, tcfg%tideamp_var, tm%tideamp, &
                                 "tideamp_source", local_ierr)
            if (local_ierr /= OCEAN_STATUS_OK) then
               if (present(ierr)) ierr = local_ierr
               if (.not. present(ierr)) error stop "ocean_tidal_mixing_nml: tideamp load failed"
               return
            end if
         end if
         if (trim(tcfg%h2_source) == "file") then
            call load_tide_field(tcfg%h2_file, tcfg%h2_var, tm%h2, "h2_source", local_ierr)
            if (local_ierr /= OCEAN_STATUS_OK) then
               if (present(ierr)) ierr = local_ierr
               if (.not. present(ierr)) error stop "ocean_tidal_mixing_nml: h2 load failed"
               return
            end if
         end if
         if (tm%use_fields) then
            block
               real(wp), allocatable :: depth(:, :)
               allocate (depth(size(tm%h2, 1), size(tm%h2, 2)))
               depth(:, :) = sum(ocean_state%multilayer%h_layer, dim=3)
               call tidal_mixing_fold_static_fields(tm, size(tm%h2, 1), size(tm%h2, 2), &
                                                    depth, ocean_state%multilayer%wet_mask)
            end block
         end if

         if (compute_rank == 0) then
            call logger%info("Tidal mixing:     interior closure enabled "// &
                             "(zeta="//to_string(tcfg%zeta)//" m  "// &
                             "E="//to_string(tcfg%e_uniform)//" W/m^2)")
         end if
         tm%enable = .true.
      end associate
      if (present(ierr)) ierr = OCEAN_STATUS_OK
   contains

      subroutine load_tide_field(fname, vname, arr, what, lierr)
         !! Read one static 2-D tide field into `arr`; fail loud (via the
         !! error ring) on a blank path, a missing variable or a grid
         !! mismatch.
         character(len=*), intent(in) :: fname, vname, what
         real(wp), intent(inout) :: arr(:, :)
         integer, intent(out) :: lierr

         lierr = OCEAN_STATUS_OK
         if (len_trim(fname) == 0) then
            call fail("ocean_tidal_mixing_nml: "//what//" = 'file' needs the "// &
                      "matching *_file path", lierr, OCEAN_STATUS_ERR_SETUP)
            lierr = OCEAN_STATUS_ERR_SETUP
            return
         end if
#ifndef RDB_NO_NETCDF
         call load_static_2d_field(trim(fname), trim(vname), arr, grid, ierr=lierr)
         if (lierr /= OCEAN_STATUS_OK) then
            call fail("ocean_tidal_mixing_nml: could not read '"//trim(vname)// &
                      "' from '"//trim(fname)//"' (missing file/variable or "// &
                      "grid mismatch)", lierr, OCEAN_STATUS_ERR_IO)
            lierr = OCEAN_STATUS_ERR_IO
         end if
#else
         call fail("ocean_tidal_mixing_nml: "//what//" = 'file' needs "// &
                   "RDB_ENABLE_NETCDF=ON", lierr, OCEAN_STATUS_ERR_IO)
         lierr = OCEAN_STATUS_ERR_IO
#endif
      end subroutine load_tide_field
   end subroutine configure_ocean_tidal_mixing

   module subroutine configure_ocean_conv(cfg, ocean_state, compute_rank, ierr)
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

   module subroutine configure_ocean_ddiff(cfg, ocean_state, compute_rank, ierr)
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

end submodule rdb_ocean_setup_vmix
