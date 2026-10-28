submodule(rdb_ocean_setup) rdb_ocean_setup_grid
   !! R2 split: grid-group configure_ocean_* bodies, moved verbatim
   !! from rdb_ocean_setup.F90 (host-only configure code; see
   !! src/core/ocean/state/rdb_ocean_setup.F90 for the public interface).
   implicit none

contains

   module subroutine configure_ocean_metrics(cfg, ocean_state, grid, compute_rank, ierr)
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

   module subroutine configure_ocean_land_mask(cfg, ocean_state, grid, compute_rank, warm_restart)
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
      !! literal no-op ⇒ bit-identical to a no-mask build — except that a
      !! face of ZERO width (`dy_cu`/`dx_cv = 0`, a node-aligned tripolar
      !! cap's pole columns) is closed even between wet cells.  This call
      !! is the only place that closure happens, so every geometry path
      !! must run it.
      type(config_t), intent(in) :: cfg
      type(ocean_state_t), intent(inout) :: ocean_state
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: compute_rank
      logical, intent(in), optional :: warm_restart
         !! `.true.` when the prognostic state was just read from a
         !! checkpoint: the land-state seed below is then SKIPPED (the
         !! masks and metrics are still derived).  Absent ⇒ cold start.

      logical :: seed_land

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
      ! Distributed tripolar fold (px > 1, plan site S2):
      ! `metrics_apply_land_mask` folds its working copy with the LOCAL
      ! kernel, which is exact only on a tile holding the whole fold row.
      ! Fold the STORED wet mask here instead, through the owner-routed
      ! exchange (collective over the north rank row, host mode), and tell
      ! the land mask not to fold again.  Folding the stored mask, not a
      ! copy, also closes the seed-time gap for it: the seed derives
      ! `wet_mask` elementwise from the water column, whose north ghosts
      ! the single-rank seed folded first (`seed_wrap_static_2d`) and the
      ! px > 1 seed could not, so the stored ghosts now equal the serial
      ! run's bit for bit (an elementwise map commutes with the fold copy).
      if (ocean_fold_is_distributed()) then
         call ocean_fold_wrap_eta_2d(grid, ocean_state%bc, ocean_state%multilayer%wet_mask, &
                                     device_resident=.false.)
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
      ! The periodic flags select a LOCAL wrap of the working wet mask's
      ! ghosts, which is right only on an axis this rank holds whole: on a
      ! decomposed periodic axis the wrap ghosts are an MPI seam the
      ! exchange above already filled, and a local wrap overwrote them with
      ! the tile's own opposite edge (every rank of a px >= 3 periodic run
      ! with land near a seam then masked the wrong ghost faces; the
      ! decomposed-vs-serial bit-identity test caught it on 4x1).  The wall
      ! flags below carry the edge tags, so nothing else reads these two.
      call metrics_apply_land_mask(ocean_state%metrics, &
                                   ocean_state%multilayer%wet_mask, grid, &
                                   ocean_state%bc%periodic_x .and. .not. ocean_halo_is_decomposed_x(), &
                                   ocean_state%bc%periodic_y .and. .not. ocean_halo_is_decomposed_y(), &
                                   ocean_state%bc%north_fold .and. .not. ocean_fold_is_distributed(), &
                                   mask_wall_velocity=cfg%ocean%bc%mask_wall_velocity, &
                                   wall_west=(ocean_bc_outer_face_tag(ocean_state%bc%west%bc_type) == OBC_WALL &
                                              .and. ocean_state%bc%has_west), &
                                   wall_east=(ocean_bc_outer_face_tag(ocean_state%bc%east%bc_type) == OBC_WALL &
                                              .and. ocean_state%bc%has_east), &
                                   wall_south=(ocean_bc_outer_face_tag(ocean_state%bc%south%bc_type) == OBC_WALL &
                                               .and. ocean_state%bc%has_south), &
                                   wall_north=(ocean_bc_outer_face_tag(ocean_state%bc%north%bc_type) == OBC_WALL &
                                               .and. ocean_state%bc%has_north))

      ! Hold land T-cells at finite reference values so masked arithmetic
      ! never multiplies 0 by a NaN (0*NaN = NaN).  Land h_layer is floored
      ! to H_VANISHED (never 0), tracers held at the IC reference, layer +
      ! face velocities zeroed on land.
      !
      ! Cold start only.  The seed is the land state's INITIAL value, not an
      ! invariant the step maintains: the ALE remap regrids a land column
      ! like any other (a `z_fixed` land column leaves step 1 as `nz-1`
      ! `zstar_h_min` fillers over a bed cell holding the rest of its
      ! `nz·H_VANISHED`), and the checkpoint saves that state verbatim.
      ! Re-seeding after a restart read rewrote every land column back to
      ! uniform `H_VANISHED` -- column total unchanged, layout not -- so the
      ! resumed run was not the run that wrote the file (1/4-degree Southern
      ! Ocean, 2026-10-01: 56 046 land columns x 50 layers per rank).
      seed_land = .true.
      if (present(warm_restart)) seed_land = .not. warm_restart
      if (seed_land) call ocean_state_seed_land_cells(ocean_state, grid)

      if (compute_rank == 0) then
         associate (unused => cfg%ocean%grid%grid_config)
         end associate
      end if
   end subroutine configure_ocean_land_mask

   module subroutine configure_ocean_forcing(cfg, ocean_state, grid, compute_rank, decomp)
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
      ! HK pair floor on the coordinates with STATIC bed fillers, closed
      ! faces or not (see `coriolis_adv_t%hk_pair_floor`).
      ocean_state%coriolis_adv%hk_pair_floor = &
         any(parse_ocean_vcoord_type(cfg%vcoord_type) == &
             [VCOORD_Z_FIXED, VCOORD_ZSTAR, VCOORD_ZSTAR_FULL])
      ! Coastal lateral BC (land mask, C1): free-slip default, shared knob.
      ocean_state%coriolis_adv%no_slip = cfg%ocean%hvisc%no_slip
      ! Mass-consistent CorAdCalc (MOM6 parity); config fail-loud
      ! guarantees form="sadourny_energy" + split_scheme="pred_corr".
      ocean_state%coriolis_adv%state_fluxes = cfg%ocean%coriolis%use_state_fluxes
      ! BOUND_CORIOLIS velocity-form clamp (energy form, or enstrophy form
      ! with centred PV interpolation; config fail-loud guarantees it).
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
            call logger%info("Coriolis bound:   BOUND_CORIOLIS ON — "//trim(cfg%ocean%coriolis%form)// &
                             " accel clamped to the (f+zeta)*v velocity-form range (MOM6)")
         end if
         if (trim(cfg%ocean%coriolis%corner_h) == "mom6_area") then
            call logger%info("Coriolis corner:  PV corner-h = mom6_area (MOM6 "// &
                             "Area_q/(hArea_q+vol_neglect); round-off-identical to "// &
                             "cell_mean above H_MIN_PV)")
         end if
      end if
   end subroutine configure_ocean_forcing

   module subroutine configure_ocean_z_fixed_profile(cfg, ocean_state, compute_rank, log_it)
      !! Resolve the `VCOORD_Z_FIXED` nominal layering onto the vcoord slot —
      !! and the `VCOORD_ZSTAR` one, which is the same nominal profile
      !! (MOM6 z* dilates it per column; `ocean_vcoord_zstar_target`):
      !! `z_fixed_h_ref = &ocean_topo_nml max_depth` (the uniform
      !! `max_depth/nz` spacing — the default, byte-identical), or, under
      !! `&vcoord_nml z_fixed_profile = "list" | "tanh"`, the stretched
      !! per-layer tables `z_fixed_zi` / `z_fixed_dz` built by
      !! `rdb_vcoord :: z_fixed_nominal_dz` (`z_fixed_h_ref` then becomes
      !! the profile's total depth).
      !!
      !! The same tables are the z* NOMINAL FLOOR of `VCOORD_HYCOM` (MOM6
      !! HYCOM1 floors its interfaces at the z* `coordinateResolution` its
      !! ALE_COORDINATE_CONFIG defines), so this also runs for `hycom`.
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
      ! Three readers of the z* nominal resolution: `z_fixed` (its levels),
      ! `zstar` (MOM6 z*, its levels) and `hycom` (its z* nominal floor,
      ! MOM6 HYCOM1 coordinateResolution).
      if (ocean_state%vcoord%coord_type /= VCOORD_Z_FIXED .and. &
          ocean_state%vcoord%coord_type /= VCOORD_ZSTAR .and. &
          ocean_state%vcoord%coord_type /= VCOORD_HYCOM) return
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
         call logger%info(trim(cfg%vcoord_type)//" z* profile:  "//trim(cfg%z_fixed_profile)// &
                          " — nominal dz "//to_string(dz(1))//" m (surface) … "// &
                          to_string(dz(nz))//" m (bed), total "// &
                          to_string(ocean_state%vcoord%z_fixed_h_ref)//" m over "// &
                          to_string(nz)//" layers")
      end if
   end subroutine configure_ocean_z_fixed_profile

   module subroutine configure_ocean_k_top(cfg, ocean_state, grid, compute_rank)
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

   module subroutine configure_ocean_k_bot(ocean_state, grid, compute_rank)
      !! Fill `ms%k_bot` / `k_bot_u` / `k_bot_v` — the shared index of the
      !! first LIVE layer counting UP from the bed, and the field every
      !! bed-side consumer reads instead of spelling `1`.  The bed-side
      !! mirror of `configure_ocean_k_top`.
      !!
      !! Under `vcoord_type = "z_fixed"` (with a resolved `z_fixed_h_ref`)
      !! every column shallower than the nominal stack carries inert
      !! FILLERS at `zstar_h_min` below its partial bed cell, so `k = 1`
      !! is not the bed-adjacent layer: a bottom drag, an implicit-drag
      !! diagonal, a geothermal deposit or a tidal-mixing bed anchor put
      !! on `k = 1` lands on a layer the vdiff solve has cut out of the
      !! column.  Unlike `k_top` this is NOT gated on a cavity — bed
      !! fillers exist on every `z_fixed` run with topography.
      !!
      !! **Static.** Read ONCE from `ocean_vcoord_z_fixed_target` at
      !! `eta = 0` — the same kernel, `bt_H_ref`, `z_top` and profile the
      !! closed-face mask and `k_top` are built from (one definition of
      !! "live").  The bed is static and `eta` is absorbed by the first
      !! live layer at the TOP, so the bed-side pattern never moves.  It
      !! is derived from bathymetry + the vcoord config, so it is rebuilt
      !! on every start (cold or warm) and is not restart state.
      !!
      !! **Every other coordinate is a literal no-op**: the arrays were
      !! allocated at `source = 1` in `multilayer_state_init`, which IS
      !! the answer wherever nothing vanishes against the bed.
      !!
      !! **Seams.** The centre index is built from the ghost-filled
      !! `bt_H_ref`, so it is seam-correct by construction; the OUTERMOST
      !! face of each face twin cannot be evaluated locally (it needs a
      !! cell beyond the array) and takes its one column.  A plain face
      !! halo exchange (integer → real → integer: the halo layer moves
      !! `real(wp)`, and an integer < 2**53 round-trips exactly) makes
      !! that ghost face carry its owner's value on a tile seam or a
      !! periodic wrap.  Single-rank non-periodic: a no-op.
      !!
      !! **Ordering.** Same as `configure_ocean_k_top`: after the cavity
      !! draft, after `configure_ocean_bt_split` (`bt_H_ref`), after the
      !! periodic-wrap / halo pass, BEFORE `ocean_state_enter_data` (the
      !! host fill is what the `copyin` captures).
      type(ocean_state_t), intent(inout) :: ocean_state
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: compute_rank

      integer :: nx, ny, nz, n_filler
      real(wp) :: h_nominal, h_min
      real(wp), allocatable :: tgt(:, :, :), eta0(:, :), fx(:, :), fy(:, :)

      if (.not. ocean_state%multilayer%is_init) return
      if (ocean_state%vcoord%coord_type /= VCOORD_Z_FIXED) return
      if (ocean_state%vcoord%z_fixed_h_ref <= 0.0_wp) return
      if (.not. ocean_state%dyn%bt_work%is_init) return
      if (maxval(ocean_state%dyn%bt_work%bt_H_ref) <= 0.0_wp) return

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
      call ocean_vcoord_k_bot_from_target(ocean_state%multilayer%k_bot, &
                                          ocean_state%multilayer%k_bot_u, &
                                          ocean_state%multilayer%k_bot_v, &
                                          tgt, nx, ny, nz, H_VANISHED)
      deallocate (tgt, eta0)

      if (ocean_halo_is_init()) then
         allocate (fx(nx + 1, ny), fy(nx, ny + 1))
         fx = real(ocean_state%multilayer%k_bot_u, wp)
         fy = real(ocean_state%multilayer%k_bot_v, wp)
         call ocean_halo_face_x(fx, device_resident=.false.)
         call ocean_halo_face_y(fy, device_resident=.false.)
         ocean_state%multilayer%k_bot_u = nint(fx)
         ocean_state%multilayer%k_bot_v = nint(fy)
         deallocate (fx, fy)
      end if

      n_filler = count(ocean_state%multilayer%k_bot > 1)
      if (compute_rank == 0) then
         call logger%info("z_fixed k_bot: "//to_string(n_filler)//"/"// &
                          to_string(nx*ny)//" columns carry bed-side "// &
                          "fillers (k_bot > 1); max k_bot = "// &
                          to_string(maxval(ocean_state%multilayer%k_bot))// &
                          ", nz = "//to_string(nz))
      end if
   end subroutine configure_ocean_k_bot

   module subroutine configure_ocean_closed_faces(cfg, ocean_state, grid, compute_rank, ierr)
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
      !! the coordinate's target at `η = 0` (`ocean_vcoord_eta0_target`)
      !! — the same kernel the ALE regrid and the IC seed use, so there
      !! is no second definition of "live".  The mask is STATIC: the bed
      !! and the draft are static, and under `z_fixed` `η` is absorbed by
      !! the first LIVE layer (the partial cell), so the live/filler
      !! pattern does not move to first order in `η/h_partial`.
      !!
      !! **`zstar_full` too.**  `build_zref_full` lays a z-level fine
      !! zone (`h_surf` layers from the surface) over a terrain-following
      !! coarse zone; a column shallower than the fine zone ends in a
      !! partial cell and every layer below it is a `zstar_h_min` filler —
      !! the same staircase as `z_fixed`'s bed, with the same
      !! thickness-independent FV-PGF defect across it.  Its target puts
      !! `η ≥ 0` in the surface layer (pattern EXACTLY static) and clips
      !! `η < 0` from the bed (a bed-most live layer thinner than `|η|`
      !! flips — fewer columns than `z_fixed` flips at the same `|η|` on
      !! the 1-degree Southern Ocean).  So the mask is built from the
      !! `ZSTAR_FULL` target exactly as it is from the `z_fixed` one, and
      !! every consumer is reused unchanged.  It closes FILLER faces
      !! only; the terrain-following coarse zone's open faces keep the
      !! sigma PGF error.
      !!
      !! **`zstar` (MOM6 z*).**  `ocean_vcoord_zstar_target` dilates the
      !! `z_fixed` nominal profile by the column's free-surface stretching
      !! and decides every layer's liveness at `η = 0`, so its pattern is
      !! EXACTLY static for `η` of either sign (MOM6 `build_zstar_column`:
      !! the dilation keeps the ratios) and its `η = 0` target is the
      !! `z_fixed` one bit for bit — the same mask, every consumer reused.
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
      real(wp) :: h_face, sum_all, sum_open
      real(wp), allocatable :: tgt(:, :, :)
      character(len=16) :: vcoord_label

      if (present(ierr)) ierr = OCEAN_STATUS_OK
      if (.not. cfg%zfixed_closed_faces) then
         call refuse_open_zfixed_staircase(ocean_state, grid, ierr)
         return
      end if

      if (.not. ocean_state%multilayer%is_init) then
         call fail("&vcoord_nml zfixed_closed_faces requires the ocean "// &
                   "multilayer path (the mask is per-layer)", ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if
      if (ocean_state%vcoord%coord_type /= VCOORD_Z_FIXED .and. &
          ocean_state%vcoord%coord_type /= VCOORD_ZSTAR .and. &
          ocean_state%vcoord%coord_type /= VCOORD_ZSTAR_FULL) then
         call fail("&vcoord_nml zfixed_closed_faces is only defined for "// &
                   "vcoord_type='z_fixed', 'zstar' and 'zstar_full': the "// &
                   "live/filler staircase it closes is made by a GEOMETRIC "// &
                   "coordinate that vanishes bed-side layers at fixed "// &
                   "reference depths. On sigma / zstar_sigma every layer is live on "// &
                   "every wet face, and rho / hycom vanish layers by DENSITY, "// &
                   "so their pattern is not static", &
                   ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if
      if ((ocean_state%vcoord%coord_type == VCOORD_Z_FIXED .or. &
           ocean_state%vcoord%coord_type == VCOORD_ZSTAR) .and. &
          ocean_state%vcoord%z_fixed_h_ref <= 0.0_wp) then
         call fail("&vcoord_nml zfixed_closed_faces needs a resolved "// &
                   "z_fixed_h_ref (set &ocean_topo_nml max_depth): without "// &
                   "it the z_fixed / zstar target degenerates to uniform sigma, "// &
                   "there are no fillers, and the mask would close nothing", &
                   ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if
      if (ocean_state%vcoord%coord_type == VCOORD_ZSTAR_FULL .and. &
          ocean_state%vcoord%zstar_h_surf_target <= 0.0_wp) then
         call fail("&vcoord_nml zfixed_closed_faces under vcoord_type="// &
                   "'zstar_full' needs zstar_h_surf_target > 0: without it "// &
                   "build_zref_full lays uniform sigma, there are no fillers, "// &
                   "and the mask would close nothing", &
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
      ! GM composes: `gm_column_x/y` build the streamfunction on each face's
      ! OPEN column (`metrics%open_u/open_v` x live on both sides) and the
      ! slopes slot masks slope / N^2 to it, so `uhD`/`vhD` are zero on
      ! every closed face-layer by construction (`test_ocean_gm_zfixed`).
      ! Redi composes: each face pairs only its OPEN WINDOW (the contiguous
      ! open + live-both-sides run), so the neutral-surface sweep, the PPM
      ! reconstruction and the flux scatter never touch a closed face-layer
      ! or a filler (`test_ocean_redi_zfixed`).
      ! Fox-Kemper MLE composes: the ML walk skips fillers and each face
      ! builds its overturning on its OPEN column (masked face thickness,
      ! `H_vel` clamped to it), so `uhml`/`vhml` — folded into
      ! mass_flux_*_layer after the per-layer mask — are zero on every
      ! closed face-layer and filler and still sum to zero per face
      ! (`test_ocean_mle_zfixed`).
      ! The velocity-form BIHARMONIC (scalar `nu_4` and the flow-aware
      ! `smag_ah` / `leith_biharm` `nu4_face_*`) composes: both chained
      ! Laplacians gate every difference by `open(a)*open(b)` (free-slip,
      ! `hvisc_biharm_lap_closed`).  The stress-tensor assembly (and
      ! `kh_aniso`, which only it reads) does not: its T-cell tension and
      ! corner shear are built from the face velocities with the 2-D
      ! `wet_*` masks only.
      if (cfg%ocean%hvisc%stress_tensor) then
         call fail("&vcoord_nml zfixed_closed_faces does not yet compose "// &
                   "with &ocean_hvisc_nml stress_tensor (nor kh_aniso, which "// &
                   "only the stress-tensor path reads): its T-cell tension and "// &
                   "corner shear are masked by the 2-D wet_* fields only, so a "// &
                   "closed face-layer's zeroed velocity enters the strain as a "// &
                   "Dirichlet value (no-slip at every staircase step).  The "// &
                   "velocity-Laplacian harmonic kernels and both biharmonic "// &
                   "paths (nu_4, smag_ah, leith_biharm) carry the free-slip "// &
                   "closure; use those", &
                   ierr, OCEAN_STATUS_ERR_SETUP)
         return
      end if
      ! ---- Barotropic paths still on FULL-column weights ----
      !
      ! `derive_bt_from_layers`, `face_depth_mean_*`, `set_cor_ref_velocity`,
      ! the `apply_bt_correction` fold, `compute_h_face_upstream`
      ! (`upstream_h_face`), `compute_bt_rem` (`substep_drag`) and
      ! `compute_bt_rem_wave_drag` (`wave_drag`) weight by `h_face·open`
      ! under this knob.  The bc-PGF correction does NOT: it sums
      ! `0.5·(h_L + h_R)` over EVERY layer.  Its answer is not wrong by a
      ! round-off — the closed layers' thickness is O(h_nominal) at a
      ! staircase face — so it is refused until it is ported, rather than
      ! left to run on the wrong column.
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

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ocean_state%multilayer%nz_ml
      vcoord_label = "z_fixed"
      if (ocean_state%vcoord%coord_type == VCOORD_ZSTAR_FULL) vcoord_label = "zstar_full"
      if (ocean_state%vcoord%coord_type == VCOORD_ZSTAR) vcoord_label = "zstar"

      call metrics_closed_faces_alloc(ocean_state%metrics, grid, nz)

      ! The eta = 0 target of the running coordinate — `z_fixed`, `zstar`
      ! or `zstar_full` — through the SAME kernel the ALE regrid dispatches
      ! to, so "live" has one definition.  Under `zstar_full` it walks the
      ! per-column `z_ref` table the engine rebuilt from the
      ! halo-exchanged bathymetry just before the static-geometry pass.
      allocate (tgt(nx, ny, nz), source=0.0_wp)
      call ocean_vcoord_eta0_target(ocean_state%vcoord, tgt, &
                                    ocean_state%dyn%bt_work%bt_H_ref, nx, ny, nz)
      call ocean_vcoord_closed_face_masks(ocean_state%metrics%open_u, &
                                          ocean_state%metrics%open_v, &
                                          tgt, nx, ny, nz, H_VANISHED)
      ! The builder cannot evaluate the OUTERMOST face of the array (it
      ! needs a cell beyond it) and leaves it open.  On a tile that face
      ! is a SEAM ghost, three cells from the owned region, and on a
      ! periodic edge it is the wrap of an interior face — so without
      ! this exchange the ghost band of the mask is decomposition-
      ! dependent.  The harmonic kernels never read that deep; the
      ! biharmonic's chained stencil does (its ghost-band tendency at
      ! depth 1 reads the mask at depth 3), which broke 1x4 bit-identity
      ! in `test_ocean_decomp_bitid_mpi` (`file_readers`).  A plain face
      ! exchange (no sign flip: a 0/1 mask, not a vector) makes every
      ! ghost face carry its owner's value.  Single-rank non-periodic:
      ! a no-op.
      if (ocean_halo_is_init()) then
         call ocean_halo_face_x(ocean_state%metrics%open_u, nz, device_resident=.false.)
         call ocean_halo_face_y(ocean_state%metrics%open_v, nz, device_resident=.false.)
      end if

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

      deallocate (tgt)

      ! Flip the switches LAST — every consumer branches on them, and the
      ! mask has to be in place before any of them can read a closed face.
      ocean_state%metrics%use_closed_faces = .true.
      ocean_state%vcoord%zfixed_closed_faces = .true.
      ocean_state%vdiff%zlevel_faces = .true.
      ! #178: `ocean_bottom_drag_t%zlevel_faces` mirrors `vdiff`'s so the two
      ! modules' `kb_live` searches (shared `rdb_blf_is_live` criterion) agree
      ! on which row is the live bed row of a given face.
      ocean_state%bdrag%zlevel_faces = .true.

      if (compute_rank == 0) then
         call logger%info(trim(vcoord_label)//" closed faces: ON (partial steps, "// &
                          "Adcroft/Hill/Marshall 1997) — u closed "// &
                          to_string(n_closed_u)//"/"// &
                          to_string(n_closed_u + n_open_u)//", v closed "// &
                          to_string(n_closed_v)//"/"// &
                          to_string(n_closed_v + n_open_v)// &
                          ", isolated ledge cells "//to_string(n_ledge))
         if (n_ledge > 0) then
            call logger%warning(trim(vcoord_label)//" closed faces: "//to_string(n_ledge)// &
                                " LIVE cells have all four own-layer faces "// &
                                "closed — the mask has isolated water.  They "// &
                                "are inert (no flux in or out, velocity zeroed "// &
                                "every stage) but a one-cell spike in the bed "// &
                                "or the draft is worth looking at.")
         end if
      end if
      if (present(ierr)) ierr = OCEAN_STATUS_OK
   end subroutine configure_ocean_closed_faces

   module subroutine refuse_open_zfixed_staircase(ocean_state, grid, ierr)
      !! The `zfixed_closed_faces = .false.` leg of
      !! `configure_ocean_closed_faces`: refuse `vcoord_type = "z_fixed"`
      !! with OPEN staircase faces over a STEPPED bed.
      !!
      !! A face where the two columns' bed falls in different nominal
      !! `z_fixed` layers pairs a live layer with a bed filler.  Left
      !! open, the FV pressure gradient across that step drives
      !! `|ρ′|·g·Δz_step/(ρ₀·dx)` out of rest whatever the filler
      !! thickness — measured on the 1° Southern Ocean: u ≈ 250 m/s by
      !! step 3, 4e9 m/s on the filler faces, `h → −2e9`, the first
      !! non-finite in `post_bt`, with or without any closure.  That is
      !! the case the closed-face mask exists for, so the combination
      !! fails loud here instead of blowing up three steps in.
      !!
      !! A flat or step-free bed is ACCEPTED untouched (no state is
      !! written, so it stays bit-identical), and so is every coordinate
      !! but `z_fixed`.  Top-side (ice-draft) steps are not counted.
      !!
      !! MPI: each rank counts its OWNED faces off the same GHOST-FILLED
      !! `bt_H_ref` / `z_top` the mask builder reads (so a tile-seam
      !! step is seen by exactly one rank) and the count is summed over
      !! the compute communicator, so every rank takes the same decision.
      !! Every early return is on rank-uniform state, so either all
      !! ranks reach the collective or none does.
      type(ocean_state_t), intent(in) :: ocean_state
      type(hgrid_t), intent(in) :: grid
      integer, intent(out), optional :: ierr

      integer :: nx, ny, nz, ng, n_local
      real(wp) :: h_nominal, steps_local, steps_global
      real(wp), allocatable :: tgt(:, :, :), eta0(:, :)

      if (present(ierr)) ierr = OCEAN_STATUS_OK
      if (.not. ocean_state%multilayer%is_init) return
      if (ocean_state%vcoord%coord_type /= VCOORD_Z_FIXED) return
      if (ocean_state%vcoord%z_fixed_h_ref <= 0.0_wp) return
      if (.not. ocean_state%dyn%bt_work%is_init) return

      nx = grid%nx_total
      ny = grid%ny_total
      ng = grid%nghost
      nz = ocean_state%multilayer%nz_ml
      h_nominal = ocean_state%vcoord%z_fixed_h_ref/real(nz, wp)

      allocate (tgt(nx, ny, nz), source=0.0_wp)
      allocate (eta0(nx, ny), source=0.0_wp)
      call ocean_vcoord_z_fixed_target(tgt, ocean_state%dyn%bt_work%bt_H_ref, &
                                       eta0, ocean_state%vcoord%z_top, &
                                       nx, ny, nz, h_nominal, &
                                       ocean_state%vcoord%z_fixed_use_profile, &
                                       ocean_state%vcoord%z_fixed_zi, &
                                       ocean_state%vcoord%z_fixed_dz, &
                                       ocean_state%vcoord%zstar_h_min)
      n_local = ocean_vcoord_count_bed_steps(tgt, ocean_state%dyn%bt_work%bt_H_ref, &
                                             ocean_state%metrics%dy_cu, &
                                             ocean_state%metrics%dx_cv, &
                                             nx, ny, nz, ng + 1, ng + grid%nx_phys, &
                                             ng + 1, ng + grid%ny_phys, H_VANISHED)
      deallocate (tgt, eta0)

      ! An integer count rides the real64 sum exactly (far below 2**53).
      steps_local = real(n_local, wp)
      call halo_allreduce_sum(steps_local, steps_global)
      if (steps_global < 0.5_wp) return

      call fail("&vcoord_nml zfixed_closed_faces = .false. is refused with "// &
                "vcoord_type='z_fixed' over a stepped bed: "// &
                to_string(nint(steps_global))//" wet faces join columns whose "// &
                "bed lies in different nominal z_fixed layers, so an open face "// &
                "pairs a live layer with a bed filler and the pressure gradient "// &
                "across the open staircase step drives flow from rest (blows "// &
                "up within steps).  Set &vcoord_nml zfixed_closed_faces = .true.", &
                ierr, OCEAN_STATUS_ERR_SETUP)
   end subroutine refuse_open_zfixed_staircase

   module subroutine configure_ocean_porous(cfg, ocean_state, grid, compute_rank, ierr)
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

   module subroutine configure_ocean_wetdry(cfg, ocean_state, grid, compute_rank)
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

end submodule rdb_ocean_setup_grid
