submodule(rdb_ocean_setup) rdb_ocean_setup_lateral
   !! R2 split: lateral-group configure_ocean_* bodies, moved verbatim
   !! from rdb_ocean_setup.F90 (host-only configure code; see
   !! src/core/ocean/state/rdb_ocean_setup.F90 for the public interface).
   implicit none

contains

   module subroutine configure_ocean_lateral(cfg, ocean_state, grid, compute_rank, ierr)
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
                             to_string(cfg%rho_ref_pressure)//" Pa, targets "// &
                             to_string(ocean_state%vcoord%rho_target(0))//" .. "// &
                             to_string(ocean_state%vcoord%rho_target( &
                                       ocean_state%vcoord%nz_ml))//" kg/m^3 ("// &
                             trim(cfg%rho_target_profile)//"))")
         else if (ocean_state%vcoord%coord_type == VCOORD_HYCOM) then
            ! HYCOM = the RHO density-space inversion + a z* nominal-floor
            ! sweep (fixed-resolution surface band, isopycnal interior).
            ! Reuses the RHO rho_target + EOS; the floor is the z* nominal
            ! profile in metres (z_fixed_profile).  Needs S+T.
            call logger%info("Vertical coord:  hycom (hybrid z*/isopycnal, ALE "// &
                             "remap; needs S+T + EOS; rho_ref_p="// &
                             to_string(cfg%rho_ref_pressure)//" Pa, targets "// &
                             to_string(ocean_state%vcoord%rho_target(0))//" .. "// &
                             to_string(ocean_state%vcoord%rho_target( &
                                       ocean_state%vcoord%nz_ml))//" kg/m^3 ("// &
                             trim(cfg%rho_target_profile)//"); z* floor "// &
                             trim(cfg%z_fixed_profile)//", top dz "// &
                             to_string(hycom_top_dz(ocean_state%vcoord))//" m)")
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
      ! MANDATORY re-derive: kv_bg/kt_bg/ks_bg and the kt/ks/kd_bg arrays
      ! were set from the TYPE-DEFAULT pp81_* at `init` time (before this
      ! configure call runs); without this, `vmix_assemble`'s background
      ! floor would silently clamp against the old default even though
      ! the user just set a new pp81_nu_bg/pp81_kappa_bg, and a
      ! Bryan-Lewis/Henyey background config would lose `kd_bg` entirely
      ! on every warm restart (`kd_bg` is set ONLY here) — see
      ! `vmix_seed_backgrounds`'s docstring.  `kv` alone is skipped when
      ! `engine_setup`'s restart read actually found it in the checkpoint
      ! (PR-2, bt-rem-from-av-rem: `vmix%kv_from_restart`, set by
      ! `ocean_state_restart_read`) -- it is the one CARRIED field here
      ! (`visc_rem_precompute` reads the PREVIOUS stage's `kv`), so an
      ! unconditional reseed would stomp the just-restored value back to
      ! the cold background on every resume; a cold start or an older
      ! checkpoint without `vmix_kv` leaves `kv_from_restart = .false.`
      ! and `kv` still reseeds normally.  Runs before
      ! `ocean_state_enter_data` (driver), so no `!$acc update` is needed
      ! either way -- the host array this leaves in place (restored or
      ! freshly seeded) is what the later `copyin` maps.
      call ocean_state%vmix%seed_backgrounds(skip_kv=ocean_state%vmix%kv_from_restart)

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

   module subroutine configure_rho_target(cfg, rho_target, nz_ml)
      !! Populate the isopycnal `rho_target(0:nz_ml)` interface densities
      !! for `VCOORD_RHO` / `VCOORD_HYCOM`.  `rho_target(0)` is the surface
      !! (lightest) interface, `rho_target(nz_ml)` the bed (densest).
      !!
      !!   * `&vcoord_nml rho_target_profile = "uniform"` (default) — a
      !!     uniform light→dense linspace from `rho_target_light` /
      !!     `rho_target_dense` (MOM6 `ALE_COORDINATE_CONFIG=UNIFORM`);
      !!   * `"list"` — the `nz_ml+1` interface densities of
      !!     `rho_target_list`, light first (MOM6 target densities from a
      !!     list or a file, e.g. `HYBRID:file,sigma2,dz`).
      !!
      !! `validate_config` has already refused a list of the wrong length
      !! or one that is not strictly increasing; a short list reaching here
      !! is an `error stop`.  Harmless on other families.
      type(config_t), intent(in) :: cfg
      real(wp), intent(inout) :: rho_target(0:)
      integer, intent(in) :: nz_ml
      integer :: k
      real(wp) :: frac

      if (trim(cfg%rho_target_profile) == "list") then
         ! validate_config refuses a short list; reaching here with one is
         ! a bug, so stop rather than fall back to the linspace.
         if (count(cfg%rho_target_list > 0.0_wp) < nz_ml + 1) then
            error stop "configure_rho_target: rho_target_list shorter than nz_ml+1"
         end if
         do k = 0, nz_ml
            rho_target(k) = cfg%rho_target_list(k + 1)
         end do
         return
      end if
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

   module subroutine configure_ocean_wavespeed(cfg, ocean_state, grid, ierr)
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

   module subroutine configure_ocean_varmix(cfg, ocean_state, grid)
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

   module subroutine configure_ocean_meke(cfg, ocean_state, grid)
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

   module subroutine configure_ocean_reference_density(ocean_state, geo)
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

   module subroutine configure_ocean_pgf(cfg, ocean_state, compute_rank, ierr)
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
      ! Grounded = non-overlapping AND vanished on one side: the gate never
      ! zeroes the PGF of a layer that is massive in both columns.
      ocean_state%pressure_force%nonoverlap_vanish_tol = &
         nonoverlap_vanish_tol_for(cfg%ocean%isopycnal%angstrom_h)
      if (ocean_state%pressure_force%skip_nonoverlap) then
         if (compute_rank == 0) then
            call logger%info("Isopycnal PGF:    pgf_skip_nonoverlap = ON "// &
                             "(grounded-layer face PGF zeroed where layers do not overlap "// &
                             "and one side is <= "// &
                             to_string(ocean_state%pressure_force%nonoverlap_vanish_tol)//" m)")
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
      ocean_state%pressure_force%insitu_density = cfg%ocean%pgf%insitu_density
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
      ! The bc-PGF retro-correction reads `pgf%e_face`, which only FV_MOM6
      ! fills.  Mirrors the `validate_config` refusal so an in-memory /
      ! API caller fails at configure, not with an `error stop` in step 1.
      if (cfg%ocean%bt%correction_bc_pgf .and. &
          ocean_state%pressure_force%variant /= OPGF_VARIANT_FV_MOM6) then
         call fail("configure_ocean_pgf: &ocean_bt_nml correction_bc_pgf=.true. "// &
                   "requires form='fv_mom6'. compute_pbce builds the per-layer "// &
                   "pressure response from the FV_MOM6 interface-height stack "// &
                   "(pgf%e_face), which no other PGF form fills.", &
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
      ! FV_MOM6: BT gravity is ALWAYS gfs_scale*GRAVITY (the Pass-5
      ! Montgomery correction supplies the complementary share when
      ! gfs_scale < 1). Gate must NOT condition on nz_layers (NK=1 SSH
      ! bug, 2026-05-26).
      !
      ! UNCONDITIONAL on gfs_scale (fixed 2026-09-28, was gated on
      ! `gfs_scale < 1.0 - 1e-12`): `pgf_free_surface_gravity` builds the
      ! slow PGF's shed free-surface term as
      ! `gfs_scale*GRAVITY*rho_ref/rho0` unconditionally (GRAVITY =
      ! 9.80665, `rdb_constants`), but with the old gate `g_bt` stayed at
      ! `barotropic_workstate_t`'s struct-default literal `9.81_wp`
      ! whenever `gfs_scale >= 1` (the common default) -- a ~3.4e-4
      ! relative mismatch between the two gravities. `set_fast_forcing_eta_pf`
      ! subtracts `g_pf*grad(eta)` from the depth-mean slow PGF so the
      ! fast loop's own live `-g_bt*grad(eta)` can re-supply exactly that
      ! term; with `g_pf /= g_bt` the two do not cancel and the fast loop
      ! carries a real leftover forcing `(g_bt - g_pf)*grad(eta)`.
      ! Invisible in an ordinary run (eta gradients are generated by real
      ! dynamics, and the bias is a tiny wave-speed error), but exposed
      ! bit-for-bit whenever `eta` starts with a real static gradient at
      ! rest -- exactly what `&ocean_cavity_dyn_nml trim_ic_for_p_surf`
      ! seeds under a sloping ice draft. Root cause of the
      ! `vcm_lid_slope_sigma_unstrat` rest-state blip (peak En 9.997e-13);
      ! see `tests/test_ocean_cavity_load.F90::test_trim_uniform_rho_rest`
      ! (the end-to-end regression; the raw PGF face force was ALREADY
      ! round-off before this fix -- see the sibling
      ! `test_trim_balances_pfu_uniform_rho`).
      if (ocean_state%pressure_force%variant == OPGF_VARIANT_FV_MOM6) then
         ocean_state%dyn%bt_work%g_bt = ocean_state%pressure_force%gfs_scale*GRAVITY
      end if
      if (present(ierr)) ierr = OCEAN_STATUS_OK
   end subroutine configure_ocean_pgf

end submodule rdb_ocean_setup_lateral
