submodule(rdb_ocean_setup) rdb_ocean_setup_cavity
   !! R2 split: cavity-group configure_ocean_* bodies, moved verbatim
   !! from rdb_ocean_setup.F90 (host-only configure code; see
   !! src/core/ocean/state/rdb_ocean_setup.F90 for the public interface).
   implicit none

contains

   module subroutine configure_ocean_cavity(cfg, ocean_state, grid, compute_rank, ierr)
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
      ! and a non-uniform barotropic-correction weight turns the
      ! uncancelled depth-uniform force into a real per-layer shear.
      ! Refused, not auto-enabled: the namelist should say what the run does.
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

   module subroutine configure_ocean_cavity_melt(cfg, ocean_state, grid, compute_rank, ierr)
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

   module subroutine configure_ocean_top_drag(cfg, ocean_state, grid, compute_rank)
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

end submodule rdb_ocean_setup_cavity
