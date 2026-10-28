submodule(rdb_ocean_setup) rdb_ocean_setup_bt
   !! R2 split: bt-group configure_ocean_* bodies, moved verbatim
   !! from rdb_ocean_setup.F90 (host-only configure code; see
   !! src/core/ocean/state/rdb_ocean_setup.F90 for the public interface).
   implicit none

contains

   module subroutine configure_ocean_tides(cfg, ocean_state, grid, compute_rank, ierr)
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

   module subroutine configure_ocean_p_surf(cfg, ocean_state, grid, compute_rank)
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

   module subroutine configure_ocean_wave_drag(cfg, ocean_state, grid, compute_rank, ierr)
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

   module subroutine wave_drag_roughness_proxy(b, wet_T, nx, ny, nghost, kappa, n_bot, h2_max, r_h)
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

   module subroutine parse_constituent_list(list, cat_idx, nconst, ierr)
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

   module subroutine configure_ocean_bt(cfg, ocean_state, grid, compute_rank)
      !! Barotropic-substep correction knobs (MOM6 frhatu h-weighting, bc-PGF
      !! retro-correction, bt_rem_u drag damping, visc_rem joint weight), the
      !! BT_cont_type / upstream-PPM h_face workspace allocations, and the
      !! rank-0 PGF/BT configuration log lines.
      type(config_t), intent(in) :: cfg
      type(ocean_state_t), intent(inout) :: ocean_state
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: compute_rank

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

      ! PR-2 (bt-rem-from-av-rem): bt_rem_u/v from the SAME viscous
      ! remnant the layered momentum solve uses (MOM6's barotropic
      ! solver). Mutually exclusive with substep_drag and bt_halo > 0
      ! (validated at configure).
      ocean_state%dyn%bt_work%bt_rem_from_visc_rem = ocean_bt_rem_from_visc_rem_on(cfg)
      ocean_state%dyn%bt_work%bt_strong_drag = cfg%ocean%bt%strong_drag
      ocean_state%dyn%bt_work%bt_rescale_strong_drag = cfg%ocean%bt%rescale_strong_drag

      ! frhat port (&ocean_bt_nml frhat_scheme): the per-layer face-thickness
      ! closure every barotropic depth mean reads (MOM6 btcalc HYBRID vs the
      ! plain two-abutting-cell arithmetic mean). Default "arithmetic" ==
      ! bit-identical to the pre-port tree.
      ocean_state%dyn%bt_work%frhat_scheme = parse_frhat_scheme(cfg%ocean%bt%frhat_scheme, &
                                                                FRHAT_ARITHMETIC)
      if (compute_rank == 0 .and. ocean_state%dyn%bt_work%frhat_scheme == FRHAT_HYBRID) then
         call logger%info("BT frhat:         hybrid ON "// &
                          "(MOM6 btcalc HVEL_SCHEME=HYBRID face-thickness closure)")
      end if
      if (compute_rank == 0 .and. ocean_bt_rem_from_visc_rem_on(cfg)) then
         block
            character(len=:), allocatable :: bt_rem_msg
            bt_rem_msg = "BT substep:       bt_rem = av_rem**(1/n_inner) ON "// &
                         "(MOM6 av_rem/bt_rem from visc_rem"
            if (cfg%ocean%bt%strong_drag) bt_rem_msg = bt_rem_msg//", strong_drag"
            if (cfg%ocean%bt%rescale_strong_drag) bt_rem_msg = bt_rem_msg//", rescale_strong_drag"
            bt_rem_msg = bt_rem_msg//")"
            call logger%info(bt_rem_msg)
         end block
      end if

      ! MOM6 planetary-only fast loop — drop live ζ_bt/∇KE from the substeps
      ! (they stay frozen inside F_bt_*_fast); Cor_ref reduces to f·v̄.
      ocean_state%dyn%bt_work%substep_zeta_ke = cfg%ocean%bt%substep_zeta_ke
      if (compute_rank == 0 .and. .not. cfg%ocean%bt%substep_zeta_ke) then
         call logger%info("BT substep:       planetary Coriolis only "// &
                          "(substep_zeta_ke=.false.; MOM6 q=f/D parity — "// &
                          "live ζ_bt/∇KE off, frozen copies stay in F_bt)")
      end if

      ! RETIRED visc_rem-weighted BT-corrector fold (wt = visc_rem/<visc_rem>_h):
      ! `correction_visc_rem` is now fail-loud at configure (D1 follow-up —
      ! MOM6's accel_layer_u never weights this fold), so this field is
      ! wired straight from the RAW (deprecated) knob, never from
      ! `visc_rem_chain` — a test that constructs `cfg` directly and sets
      ! the raw field (bypassing the nml refusal) still gets the weighted
      ! dispatch the kernel itself supports; everyone else reads 1.
      ocean_state%dyn%bt_work%bt_correction_visc_rem = cfg%ocean%bt%correction_visc_rem
      if (compute_rank == 0 .and. cfg%ocean%bt%correction_visc_rem) then
         call logger%info("BT correction:    visc_rem/<visc_rem>_h weight ON "// &
                          "(visc_rem produced by vdiff)")
      end if
      ! D1 follow-up: the visc_rem PRODUCER is decoupled from the retired
      ! weighted fold above — it runs whenever ANY real consumer needs
      ! visc_rem_u/v fresh (forcing/renorm/bt_rem_from, or visc_rem_chain,
      ! which implies all three).
      ocean_state%dyn%bt_work%bt_visc_rem_producer = ocean_bt_visc_rem_producer_on(cfg)
      ! MOM6 wt_u parity for the BT FORCING assembly (PGF_BUG.md §9):
      ! friction-damped layers stop forcing the fast loop.  Self-sufficient
      ! (see `bt_visc_rem_producer` above).
      ocean_state%dyn%bt_work%bt_forcing_visc_rem = ocean_bt_forcing_visc_rem_on(cfg)
      if (compute_rank == 0 .and. ocean_bt_forcing_visc_rem_on(cfg)) then
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
      ocean_state%dyn%bt_work%bt_renorm_visc_rem = ocean_bt_renorm_visc_rem_on(cfg)
      if (compute_rank == 0 .and. ocean_bt_renorm_visc_rem_on(cfg)) then
         call logger%info("BT renormaliser:  gamma-weighted du + u_cor (MOM6 "// &
                          "continuity inversion parity)")
      end if
      ! PR-8 / D1 follow-up: every `*_visc_rem` consumer couples against
      ! `visc_rem_u/v`, which are only ever filled away from their
      ! `source=1.0` default by the implicit bottom-drag diagonal fold
      ! (`&ocean_vdiff_nml implicit_drag`) or `bbl_glue` breaking vdiff's
      ! per-row unit-sum normalisation — without either the chain is
      ! exactly gamma==1 (a x1.0). Warn rather than abort (tidal_mixing
      ! e_uniform=0 precedent): PR-19 (visc_rem) is the named owner that
      ! fills the arrays for other configurations.  (`validate_config` also
      ! carries this check at configure; repeated here for the rank-0 log.)
      if (compute_rank == 0 .and. ocean_bt_visc_rem_producer_on(cfg) .and. &
          .not. (cfg%ocean%vdiff%implicit_drag .or. ocean_state%vdiff%bbl_glue)) then
         call logger%warning("&ocean_bt_nml forcing_visc_rem/renorm_visc_rem/"// &
                             "bt_rem_from_visc_rem/visc_rem_chain has no effect without "// &
                             "&ocean_vdiff_nml implicit_drag=.true. or bbl_glue=.true. "// &
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

   module subroutine configure_ocean_bt_split(cfg, ocean_state, grid, compute_rank)
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

      ! `&ocean_bt_nml dtbt`: OM3-style step definition (see
      ! `rdb_barotropic_dtbt`).  Initial n_inner from the seeded bathymetry
      ! `b` (Sigma h_layer = b at t = 0); `dtbt_reset_period >= 0` re-derives
      ! it from the live thicknesses every step / period in `engine_step`.
      ! Host loop: the state is not device-mapped yet.
      if (bt_dtbt_is_set(cfg%ocean%bt)) then
         block
            integer :: nxb, nyb, ng_, n_inner_derived
            real(wp), allocatable :: b3(:, :, :)
            real(wp) :: dt_max_local, dt_max_glob
            nxb = size(ocean_state%barotropic%b, 1)
            nyb = size(ocean_state%barotropic%b, 2)
            ng_ = grid%nghost
            allocate (b3(nxb, nyb, 1))
            b3(:, :, 1) = ocean_state%barotropic%b
            call bt_dtbt_dt_max(nxb, nyb, 1, ng_ + 1, ng_ + grid%nx_phys, &
                                ng_ + 1, ng_ + grid%ny_phys, b3, &
                                ocean_state%multilayer%wet_mask, &
                                ocean_state%metrics%dy_cu_bt, ocean_state%metrics%dx_cv_bt, &
                                ocean_state%metrics%idxCu, ocean_state%metrics%idyCv, &
                                ocean_state%metrics%iareaT, &
                                ocean_state%coriolis_adv%f_corner, &
                                ocean_state%dyn%bt_work%g_bt, cfg%ocean%bt%bebt, &
                                .false., dt_max_local)
            call halo_allreduce_min(dt_max_local, dt_max_glob)
            n_inner_derived = bt_dtbt_n_inner(cfg%dt_fixed, cfg%ocean%bt%dtbt, dt_max_glob)
            if (compute_rank == 0) then
               call logger%info("dtbt: dt_max = "//to_string(dt_max_glob)// &
                                " s, dtbt = "//to_string(cfg%ocean%bt%dtbt)// &
                                " (reset period "//to_string(cfg%ocean%bt%dtbt_reset_period)// &
                                " s) → n_inner = "//to_string(n_inner_derived)// &
                                " (was "//to_string(cfg%ocean%bt%n_inner)//")")
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

end submodule rdb_ocean_setup_bt
