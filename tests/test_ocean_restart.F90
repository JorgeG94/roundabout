!! Bit-exact restart roundtrip gate for the ocean dynamical core
!! (ROADMAP A1).  The acceptance gate (contract): a write/kill/read
!! cycle must reproduce step (N+1) bit-for-bit.
!!
!! Strategy:
!!   1. Build a full `ocean_state_t` (state A), seed it, run N steps.
!!   2. Write a per-rank restart from A via the registry.
!!   3. Build a FRESH `ocean_state_t` (state B), seed it (cold), then
!!      read the restart into B (overwriting the seed).
!!   4. Advance BOTH A and B one more step.
!!   5. Assert every prognostic field (h_layer, u/v_face_x/y_layer,
!!      every tracers(:)%hTr, barotropic h/u_face_x/v_face_y) is
!!      byte-identical between A and B.
!!
!! Variants:
!!   * closed-wall (default) — full-local-array semantics: the read
!!     restores the WHOLE local array (interior + wall ghosts), so the
!!     post-restart step matches bit-for-bit with no wall-mirror op.
!!   * periodic-x — seam ghost-refresh proof: the driver's wrap
!!     re-establishes the redundant seam ghosts from the restored
!!     interior so the post-restart step still matches.
!!   * physics-rich — KPP overlay ON (lagged bl_depth), nonzero surface
!!     heat flux, dt_therm_ratio=2 (rho_layer carried across thermo-skip
!!     steps).  This variant exercises vmix_bl_depth + ml_rho_layer +
!!     the closure transient checkpoints — it FAILS if any of those is
!!     dropped from the registry.
!!
!! Negative test: a decomposition-metadata mismatch is detected both via
!! the status-returning check path AND through the production read
!! wrapper (`ocean_state_restart_read` with `ierr` present returns
!! non-zero rather than error-stopping).
module test_ocean_restart
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, H_VANISHED
   use rdb_config, only: config_t
   use rdb_grid, only: hgrid_t
   use rdb_decomp, only: decomp_t, decomp_init
   use rdb_ocean_metrics, only: metrics_fill_cartesian, metrics_finalize
   use rdb_ocean_state, only: ocean_state_t, ocean_state_enter_data, &
                              ocean_state_exit_data, ocean_state_seed_from_cfg, &
                              ocean_state_restart_write, ocean_state_restart_read
   use rdb_ocean_restart, only: restart_registry_t
   use rdb_ocean_restart_io, only: ocean_restart_check_decomp
   use rdb_state, only: register_default_tracers
   use rdb_ocean_dyn, only: ocean_dyn_step_split
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t
   use rdb_ocean_periodic, only: ocean_periodic_wrap_state, ocean_periodic_wrap_centre_2d
   implicit none
   private

   public :: collect_ocean_restart_tests

   integer, parameter :: NX = 8, NY = 8, NZ = 3, NGHOST = 2
   real(wp), parameter :: DX = 10.0e3_wp, H_TOTAL = 1000.0_wp
   real(wp), parameter :: DT = 300.0_wp
   integer, parameter :: N_INNER = 20
   integer, parameter :: N_STEPS = 5
   ! Wetdry variant: 2x2 interior patch overwritten to a near-dry depth
   ! (b = 0.01 m < dry_depth = 0.05 m) inside a 0.3 m basin.  Both the A
   ! (write) and B (read) builds must patch b identically — b is NOT in
   ! the restart registry, so the bit-exact step-(N+1) gate requires the
   ! static fields to agree by construction.
   real(wp), parameter :: WD_B_PATCH = 0.01_wp
   integer, parameter :: WD_PATCH_LO = NGHOST + 2, WD_PATCH_HI = NGHOST + 3
   ! Fix-6 adversarial probe cells for the wd_wet_dyn round-trip.  Both
   ! carry a HISTORY-dependent mask value that a naive "re-seed from D"
   ! would misclassify, so the round-trip must restore the SAVED value:
   !   * HELD-WET-IN-BAND (WD_WET_I/J): a column at 0.07 m depth
   !     (dry_depth=0.05 < D=0.07 < rewet_depth=0.10) that history holds WET
   !     (wd_wet_dyn=1); a re-seed from D<rewet would flip it dry.  B cold-
   !     seeds it 0 ⇒ the read must restore 1.
   !   * EMERGED (WD_EMRG_I/J): a genuinely emerged column (b<0, held dry
   !     wd_wet_dyn=0 in A); B cold-seeds it 1 ⇒ the read must restore 0.
   real(wp), parameter :: WD_HELD_DEPTH = 0.07_wp    ! in the hysteresis band
   real(wp), parameter :: WD_EMRG_B = -0.5_wp        ! emerged bed
   integer, parameter :: WD_WET_I = NGHOST + 6, WD_WET_J = NGHOST + 6
   integer, parameter :: WD_EMRG_I = NGHOST + 6, WD_EMRG_J = NGHOST + 2

contains

   subroutine collect_ocean_restart_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("restart_bit_exact_closed_wall", test_bit_exact_closed), &
                  new_unittest("restart_bit_exact_periodic_x", test_bit_exact_periodic), &
                  new_unittest("restart_bit_exact_physics_rich", test_bit_exact_physics), &
                  new_unittest("restart_decomp_mismatch_errors", test_decomp_mismatch), &
                  new_unittest("restart_decomp_mismatch_read_wrapper", &
                               test_decomp_mismatch_read), &
                  new_unittest("restart_bit_exact_wetdry", test_bit_exact_wetdry) &
                  ]
   end subroutine collect_ocean_restart_tests

   subroutine make_cfg(cfg, physics_rich, wetdry)
      type(config_t), intent(out) :: cfg
      logical, intent(in), optional :: physics_rich
      logical, intent(in), optional :: wetdry
      cfg%sim_type = "ocean"
      cfg%nx = NX
      cfg%ny = NY
      cfg%dx = DX
      cfg%dy = DX
      cfg%nz_layers = NZ
      cfg%nghost = NGHOST
      cfg%ocean%topo%max_depth = H_TOTAL
      cfg%initial_temperature = 15.0_wp
      cfg%initial_salinity = 35.0_wp
      cfg%coriolis_f = 1.0e-4_wp
      ! Small SSH bump to drive non-trivial (but stable) dynamics.
      cfg%ocean%ic%ic_config = "geostrophic_adjustment"
      cfg%ocean%ic%ga_length_scale = 3.0_wp*DX
      cfg%ocean%ic%ga_eta_amp = 0.1_wp
      if (present(physics_rich)) then
         if (physics_rich) then
            ! KPP overlay (lagged bl_depth) + surface heat flux +
            ! thermo subcycling (rho_layer carried across skipped steps).
            cfg%ocean%vmix%use_closure = .true.
            cfg%ocean%vmix%use_kpp = .true.
            cfg%ocean%vmix%dt_therm_ratio = 2
            cfg%ocean%thermo%q_heat = 200.0_wp
         end if
      end if
      if (present(wetdry)) then
         if (wetdry) then
            ! Use a shallow domain so the 2x2 dry patch (0.01 m) is only a
            ! 30x contrast to normal cells (~0.3 m), not 100 000x.  A tiny
            ! SSH amplitude avoids negative h_layer at 0.3 m depth.
            cfg%ocean%topo%max_depth = 0.3_wp  ! override H_TOTAL for this variant
            cfg%ocean%ic%ga_eta_amp = 0.001_wp
            ! Dynamic wet/dry: enable + hysteresis depths + land_margin.
            ! ppm_limit_pos required by wet/dry (§4.3 of the plan).
            cfg%ocean%wetdry%enable = .true.
            cfg%ocean%wetdry%dry_depth = 0.05_wp
            cfg%ocean%wetdry%rewet_depth = 0.10_wp
            cfg%ocean%wetdry%land_margin = 5.0_wp
            cfg%ocean%continuity%ppm_limit_pos = .true.
         end if
      end if
   end subroutine make_cfg

   subroutine setup_state(cfg, grid, state, sf, periodic_x, physics_rich, wetdry)
      !! Mirror the driver's ocean init path: init, seed, metrics, wrap
      !! (when periodic), enter_data.  Leaves the state device-resident.
      !!
      !! wetdry variant: after seed, overwrites the 2x2 interior patch
      !! (WD_PATCH_LO:WD_PATCH_HI square) to b = WD_B_PATCH = 0.01 m,
      !! rescaling h_layer to b/NZ and each tracer's hTr by the thickness
      !! ratio (concentrations unchanged — no EOS surprises).  Sets
      !! bt_H_ref = b (what configure_ocean_bt_split does) so the
      !! bed-blocking gate sees real bed elevations, and gives the PGF its
      !! bathymetry copy (what configure_ocean_pressure does) so the
      !! patched b is not misread as a 0.29 m SSH jump under the default
      !! flat-bed assumption.  Then mirrors configure_ocean_wetdry:
      !! allocates the wd_* workspaces, seeds wd_wet_dyn from depth
      !! (patch D = 0.01 < dry_depth = 0.05 ⇒ 0), and sets
      !! continuity%use_ppm_limit_pos.  With eta_wet ≈ ga_eta_amp =
      !! 0.001 m < zb_patch + dry_depth = 0.04 m the patch faces stay
      !! CLOSED, so the dry patch survives N_STEPS — the hysteresis mask
      !! is genuinely nontrivial at restart-write time.
      type(config_t), intent(in) :: cfg
      type(hgrid_t), intent(inout) :: grid
      type(ocean_state_t), intent(inout) :: state
      type(ocean_surface_flux_t), intent(inout) :: sf
      logical, intent(in) :: periodic_x
      logical, intent(in), optional :: physics_rich
      logical, intent(in), optional :: wetdry

      integer :: nx_w, ny_w, i, j
      logical :: do_wd

      do_wd = .false.
      if (present(wetdry)) do_wd = wetdry

      call grid%init(NX, NY, NGHOST, DX, DX)
      state%multilayer%nz_ml = NZ
      call state%init(grid)
      call metrics_fill_cartesian(state%metrics, grid, grid%dx, grid%dy)
      call metrics_finalize(state%metrics)
      call ocean_state_seed_from_cfg(state, grid, cfg)
      call register_default_tracers( &
         state%multilayer%tracers(state%multilayer%idx_salinity), &
         state%multilayer%tracers(state%multilayer%idx_temperature), cfg)
      ! Physics knobs applied directly (the test does not run the full
      ! configure_ocean_* path; these mirror what it would set).
      if (present(physics_rich)) then
         if (physics_rich) then
            state%vmix%use_closure = .true.
            state%vmix%use_kpp = .true.
            state%dyn%dt_therm_ratio = 2
         end if
      end if
      if (do_wd) then
         ! Step 1: overwrite the 2x2 interior patch to b = 0.01 m
         ! (< dry_depth = 0.05 m).  Rescale h_layer to b/NZ and each
         ! tracer's hTr by the per-layer thickness ratio so concentrations
         ! (T = 15°C, S = 35 PSU) are unchanged — the EOS sees normal
         ! values (the earlier unscaled-hTr attempt produced 1.5e6 °C and
         ! NaN'd the run).
         block
            real(wp), parameter :: H_PATCH_LAYER = WD_B_PATCH/real(NZ, wp)
            real(wp) :: h_ratio
            integer :: it_t, k
            do j = WD_PATCH_LO, WD_PATCH_HI
               do i = WD_PATCH_LO, WD_PATCH_HI
                  do k = 1, NZ
                     h_ratio = H_PATCH_LAYER/state%multilayer%h_layer(i, j, k)
                     do it_t = 1, size(state%multilayer%tracers)
                        state%multilayer%tracers(it_t)%hTr(i, j, k) = &
                           state%multilayer%tracers(it_t)%hTr(i, j, k)*h_ratio
                     end do
                     state%multilayer%h_layer(i, j, k) = H_PATCH_LAYER
                  end do
                  state%barotropic%b(i, j) = WD_B_PATCH
                  state%barotropic%h(i, j) = WD_B_PATCH
               end do
            end do
         end block
         ! Step 1b (Fix-6): two extra probe cells with history-dependent
         ! masks a naive re-seed-from-D would misclassify.  Patch b + rescale
         ! h_layer/hTr so their STATIC bathymetry agrees A↔B (b is NOT in the
         ! restart; the step-(N+1) gate needs it identical by construction).
         !   HELD-WET-IN-BAND: b = 0.07 m (dry<D<rewet) — the depth seed gives
         !     wd_wet_dyn=1 (D>dry_depth); the ADVERSARY (B) cold-seeds it 0.
         !   EMERGED: b = -0.5 m (< 0) — the wetdry seed floors h_layer to
         !     2·H_VANISHED, wd_wet_dyn seeds 0 (dry); B cold-seeds it 1.
         block
            real(wp) :: h_new, h_ratio
            integer :: it_t, k
            ! held-wet-in-band cell.
            do k = 1, NZ
               h_new = WD_HELD_DEPTH/real(NZ, wp)
               h_ratio = h_new/state%multilayer%h_layer(WD_WET_I, WD_WET_J, k)
               do it_t = 1, size(state%multilayer%tracers)
                  state%multilayer%tracers(it_t)%hTr(WD_WET_I, WD_WET_J, k) = &
                     state%multilayer%tracers(it_t)%hTr(WD_WET_I, WD_WET_J, k)*h_ratio
               end do
               state%multilayer%h_layer(WD_WET_I, WD_WET_J, k) = h_new
            end do
            state%barotropic%b(WD_WET_I, WD_WET_J) = WD_HELD_DEPTH
            state%barotropic%h(WD_WET_I, WD_WET_J) = WD_HELD_DEPTH
            ! emerged cell: floor h_layer to 2·H_VANISHED (the production seed
            ! floor), rescale tracers to preserve concentration.
            do k = 1, NZ
               h_new = 2.0_wp*H_VANISHED
               h_ratio = h_new/state%multilayer%h_layer(WD_EMRG_I, WD_EMRG_J, k)
               do it_t = 1, size(state%multilayer%tracers)
                  state%multilayer%tracers(it_t)%hTr(WD_EMRG_I, WD_EMRG_J, k) = &
                     state%multilayer%tracers(it_t)%hTr(WD_EMRG_I, WD_EMRG_J, k)*h_ratio
               end do
               state%multilayer%h_layer(WD_EMRG_I, WD_EMRG_J, k) = h_new
            end do
            state%barotropic%b(WD_EMRG_I, WD_EMRG_J) = WD_EMRG_B
            state%barotropic%h(WD_EMRG_I, WD_EMRG_J) = &
               real(NZ, wp)*2.0_wp*H_VANISHED
         end block
         ! Step 2: mode-split reference depth bt_H_ref = b (what
         ! configure_ocean_bt_split does).  derive_bt_from_layers then gives
         ! bt_eta = Σh − b ≈ 0 everywhere (exactly 0 at the patch, ±ga_eta_amp
         ! elsewhere), so the bed-blocking gate sees real bed elevations:
         ! flooding the patch needs eta_wet > zb_patch + dry_depth =
         ! −0.01 + 0.05 = 0.04 m, while eta_wet ≈ 0.001 m — patch faces stay
         ! CLOSED and the patch stays dry through N_STEPS.
         state%dyn%bt_work%bt_H_ref = state%barotropic%b
         ! Step 3: PGF bathymetry copy (what configure_ocean_pressure does).
         ! Without it the Montgomery PGF assumes a flat bed and the patched b
         ! masquerades as a 0.29 m surface jump at the patch faces.
         call state%pressure_force%set_bathymetry(state%barotropic%b)
         ! Step 4: mirror configure_ocean_wetdry — knobs + lazy wd_*
         ! allocations + depth-based hysteresis-mask seed (eta = 0 at seed
         ! so D0 = b: patch D = 0.01 < dry_depth = 0.05 ⇒ 0).
         nx_w = grid%nx_total
         ny_w = grid%ny_total
         state%dyn%bt_work%wetdry_enable = .true.
         state%dyn%bt_work%wd_dry_depth = cfg%ocean%wetdry%dry_depth
         state%dyn%bt_work%wd_rewet_depth = cfg%ocean%wetdry%rewet_depth
         allocate (state%dyn%bt_work%wd_wet_dyn(nx_w, ny_w), source=1.0_wp)
         allocate (state%dyn%bt_work%wd_theta(nx_w, ny_w), source=1.0_wp)
         allocate (state%dyn%bt_work%wd_flux_x(nx_w + 1, ny_w), source=0.0_wp)
         allocate (state%dyn%bt_work%wd_flux_y(nx_w, ny_w + 1), source=0.0_wp)
         allocate (state%dyn%bt_work%wd_open_u(nx_w + 1, ny_w), source=1.0_wp)
         allocate (state%dyn%bt_work%wd_open_v(nx_w, ny_w + 1), source=1.0_wp)
         do j = 1, ny_w
            do i = 1, nx_w
               if (state%barotropic%b(i, j) < cfg%ocean%wetdry%dry_depth) then
                  state%dyn%bt_work%wd_wet_dyn(i, j) = 0.0_wp
               end if
            end do
         end do
         ! ppm_limit_pos required by wet/dry (§4.3).
         state%continuity%use_ppm_limit_pos = .true.
      end if
      if (periodic_x) then
         state%bc%periodic_x = .true.
         call wrap_host(grid, state)
      end if
      call sf%init(grid)
      if (present(physics_rich)) then
         if (physics_rich) then
            call sf%set_surface_flux_const(cfg%ocean%thermo%q_heat, 0.0_wp)
         end if
      end if
      call ocean_state_enter_data(state)
      call sf%enter_data()
   end subroutine setup_state

   subroutine wrap_host(grid, state)
      !! Re-establish periodic-x ghosts on the host (driver §1.5 order).
      type(hgrid_t), intent(in) :: grid
      type(ocean_state_t), intent(inout) :: state
      call ocean_periodic_wrap_centre_2d(state%barotropic%b, &
                                         grid%nx_total, grid%ny_total, &
                                         grid%nx_phys, grid%ny_phys, grid%nghost, &
                                         .true., .false.)
      call ocean_periodic_wrap_state(grid, state%bc, state%multilayer)
   end subroutine wrap_host

   subroutine advance(grid, state, sf)
      type(hgrid_t), intent(in) :: grid
      type(ocean_state_t), intent(inout) :: state
      type(ocean_surface_flux_t), intent(in) :: sf
      call ocean_dyn_step_split(grid, state%metrics, state%dyn, state%eos, &
                                state%coriolis_adv, state%continuity, &
                                state%pressure_force, state%hvisc, &
                                state%bdrag, state%surface_stress, &
                                state%vert_advect, state%hdiff_tracer, &
                                state%vdiff, state%vmix, &
                                state%multilayer, DT, N_INNER, sf=sf, &
                                vcoord=state%vcoord, bc=state%bc, &
                                lateral_mix=state%lateral_mix, &
                                epbl=state%epbl, kshear=state%kshear)
   end subroutine advance

   subroutine teardown(state, sf)
      type(ocean_state_t), intent(inout) :: state
      type(ocean_surface_flux_t), intent(inout) :: sf
      call sf%exit_data()
      call ocean_state_exit_data(state)
      call sf%destroy()
      call state%destroy()
   end subroutine teardown

   subroutine roundtrip(error, periodic_x, physics_rich, wetdry)
      !! The bit-exact gate, parameterised on the periodic-x + physics
      !! + wetdry flags.
      type(error_type), allocatable, intent(out) :: error
      logical, intent(in) :: periodic_x
      logical, intent(in), optional :: physics_rich
      logical, intent(in), optional :: wetdry

      type(config_t) :: cfg
      type(hgrid_t) :: gA, gB
      type(ocean_state_t) :: A, B
      type(ocean_surface_flux_t) :: sfA, sfB
      type(decomp_t) :: decomp
      character(len=*), parameter :: FN = "test_ocean_restart_rt.nc"
      real(wp) :: t_read
      integer :: step_read, s, it
      logical :: ok, phys, do_wd
      real(wp), allocatable :: bld_A(:, :), rho_A(:, :, :)
      real(wp), allocatable :: wd_A(:, :)

      phys = .false.
      if (present(physics_rich)) phys = physics_rich
      do_wd = .false.
      if (present(wetdry)) do_wd = wetdry

      call make_cfg(cfg, physics_rich=phys, wetdry=do_wd)
      call decomp_init(decomp, NX, NY, 1, 1, 0)
      call cleanup_file(FN)

      ! State A: seed, run N steps, write restart.
      call setup_state(cfg, gA, A, sfA, periodic_x, physics_rich=phys, wetdry=do_wd)
      ! Wetdry precondition: the seed placed 0s at the shallow 2x2 patch.
      if (do_wd) then
         call check(error, &
                    minval(A%dyn%bt_work%wd_wet_dyn(NGHOST + 1:NGHOST + NX, NGHOST + 1:NGHOST + NY)) < 0.5_wp, &
                    "wetdry restart: seed did not place wd_wet_dyn==0 at shallow patch")
         if (allocated(error)) then
            call teardown(A, sfA)
            call cleanup_file(FN)
            return
         end if
      end if
      do s = 1, N_STEPS
         call advance(gA, A, sfA)
      end do
      ! Wetdry gate: the hysteresis mask must still be NONTRIVIAL at write
      ! time — the dry patch survived N_STEPS of dynamics (bed-blocking gate
      ! keeps its faces closed: eta_wet ≈ 0.001 m < zb_patch + dry_depth =
      ! 0.04 m).  Without this the round-trip compare below would pass even
      ! if the registry silently dropped wd_wet_dyn (all-1 == all-1).
      if (do_wd) then
         !$acc update self(A%dyn%bt_work%wd_wet_dyn)
         call check(error, &
                    minval(A%dyn%bt_work%wd_wet_dyn(NGHOST + 1:NGHOST + NX, NGHOST + 1:NGHOST + NY)) < 0.5_wp, &
                    "wetdry restart: hysteresis mask trivial at write time (patch flooded?)")
         if (allocated(error)) then
            call teardown(A, sfA)
            call cleanup_file(FN)
            return
         end if
         ! Fix-6 adversarial HISTORY overrides (host, post-run, pre-write):
         ! stamp the two probe cells with mask values a naive re-seed-from-D
         ! would get WRONG, so the round-trip must carry the SAVED value.
         ! HELD-WET-IN-BAND: force wd_wet_dyn=1 (history holds it wet) even
         ! though D=0.07 sits below rewet_depth; EMERGED: force wd_wet_dyn=0.
         A%dyn%bt_work%wd_wet_dyn(WD_WET_I, WD_WET_J) = 1.0_wp
         A%dyn%bt_work%wd_wet_dyn(WD_EMRG_I, WD_EMRG_J) = 0.0_wp
         !$acc update device(A%dyn%bt_work%wd_wet_dyn)
         wd_A = A%dyn%bt_work%wd_wet_dyn
      end if
      call ocean_state_restart_write(A, gA, decomp, FN, real(N_STEPS, wp)*DT, N_STEPS)

      ! Snapshot A's checkpointed transients at write time (the write
      ! already pulled them host-ward via update self).  These are the
      ! exact values the restart must restore into a fresh state, so a
      ! direct compare after the read proves the round-trip of bl_depth
      ! (review #1) + rho_layer (review #2) independent of whether they
      ! feed back into the next step's prognostics.
      if (phys) then
         bld_A = A%vmix%bl_depth
         rho_A = A%multilayer%rho_layer
      end if

      ! State B: cold seed, then read restart (overwrites seed), wrap,
      ! enter_data.  Read runs on the host BEFORE enter_data.
      call gB%init(NX, NY, NGHOST, DX, DX)
      B%multilayer%nz_ml = NZ
      call B%init(gB)
      call metrics_fill_cartesian(B%metrics, gB, gB%dx, gB%dy)
      call metrics_finalize(B%metrics)
      call ocean_state_seed_from_cfg(B, gB, cfg)
      call register_default_tracers( &
         B%multilayer%tracers(B%multilayer%idx_salinity), &
         B%multilayer%tracers(B%multilayer%idx_temperature), cfg)
      if (phys) then
         B%vmix%use_closure = .true.
         B%vmix%use_kpp = .true.
         B%dyn%dt_therm_ratio = 2
      end if
      ! Wetdry B: same STATIC fields as A (b patch, bt_H_ref = b, PGF bathy —
      ! none of these are in the restart, so the step-(N+1) gate needs them to
      ! agree by construction) + wd_* workspace allocations (needed by the
      ! read to accept the optional wd_wet_dyn field and by enter_data).
      ! B's h_layer / hTr / bt_h at the patch are overwritten by the read, so
      ! only the static fields are patched here.
      if (do_wd) then
         block
            integer :: nx_w, ny_w, i_b, j_b
            do j_b = WD_PATCH_LO, WD_PATCH_HI
               do i_b = WD_PATCH_LO, WD_PATCH_HI
                  B%barotropic%b(i_b, j_b) = WD_B_PATCH
               end do
            end do
            ! Fix-6 probe cells: SAME static b as A (bt_H_ref/PGF bathy are
            ! derived from b and are NOT in the restart, so must agree A↔B).
            B%barotropic%b(WD_WET_I, WD_WET_J) = WD_HELD_DEPTH
            B%barotropic%b(WD_EMRG_I, WD_EMRG_J) = WD_EMRG_B
            B%dyn%bt_work%bt_H_ref = B%barotropic%b
            call B%pressure_force%set_bathymetry(B%barotropic%b)
            nx_w = gB%nx_total
            ny_w = gB%ny_total
            B%dyn%bt_work%wetdry_enable = .true.
            B%dyn%bt_work%wd_dry_depth = cfg%ocean%wetdry%dry_depth
            B%dyn%bt_work%wd_rewet_depth = cfg%ocean%wetdry%rewet_depth
            ! ADVERSARIAL cold seed: wd_wet_dyn = ALL ONES, deliberately
            ! different from A's saved mask (0 at the patch).  The restart
            ! read MUST overwrite these 1s with A's 0s — the direct compare
            ! below genuinely proves restoration (like the physics-rich
            ! bl_depth check, where B's cold value differs from A's saved
            ! one).  If the read silently drops the field, the compare fails.
            allocate (B%dyn%bt_work%wd_wet_dyn(nx_w, ny_w), source=1.0_wp)
            ! Fix-6: cold-seed the HELD-WET-IN-BAND probe cell to 0 — the
            ! ADVERSARY value (A saved 1 from history).  A naive re-seed-from-D
            ! (D=0.07 < rewet_depth=0.10) would indeed classify it dry(0), so
            ! the read genuinely must restore A's saved 1.  (The all-ones
            ! source already makes the EMERGED cell adversarial: A saved 0,
            ! B cold 1 ⇒ read must restore 0.)
            B%dyn%bt_work%wd_wet_dyn(WD_WET_I, WD_WET_J) = 0.0_wp
            allocate (B%dyn%bt_work%wd_theta(nx_w, ny_w), source=1.0_wp)
            allocate (B%dyn%bt_work%wd_flux_x(nx_w + 1, ny_w), source=0.0_wp)
            allocate (B%dyn%bt_work%wd_flux_y(nx_w, ny_w + 1), source=0.0_wp)
            allocate (B%dyn%bt_work%wd_open_u(nx_w + 1, ny_w), source=1.0_wp)
            allocate (B%dyn%bt_work%wd_open_v(nx_w, ny_w + 1), source=1.0_wp)
            B%continuity%use_ppm_limit_pos = .true.
         end block
      end if
      if (periodic_x) B%bc%periodic_x = .true.
      call ocean_state_restart_read(B, gB, decomp, FN, t_read, step_read)
      call check(error, step_read == N_STEPS, "restart step count mismatch")
      if (allocated(error)) return
      call check(error, abs(t_read - real(N_STEPS, wp)*DT) < 1.0e-9_wp, &
                 "restart time mismatch")
      if (allocated(error)) return
      ! Direct round-trip proof: B's freshly-read host transients must
      ! equal A's snapshot.  Seeded B starts these at 0 / EOS-fallback,
      ! so this fails outright if bl_depth / rho_layer is dropped from
      ! the registry (review #1 / #2).
      if (phys) then
         ok = arrays_identical_2d(bld_A, B%vmix%bl_depth)
         call check(error, ok, "restored vmix bl_depth /= A (dropped from registry?)")
         if (allocated(error)) return
         ok = arrays_identical_3d(rho_A, B%multilayer%rho_layer)
         call check(error, ok, "restored rho_layer /= A (dropped from registry?)")
         if (allocated(error)) return
      end if
      ! Wetdry direct round-trip: wd_wet_dyn must be bit-exact after read.
      ! (wd_* flux/theta/open are per-substep transient — NOT round-tripped.)
      ! Verdict: BIT-EXACT because wd_wet_dyn values are 0.0/1.0 (exactly
      ! representable in IEEE 754 double), and the restart carries the full
      ! local double-precision array for every registered optional field.
      if (do_wd) then
         ok = arrays_identical_2d(wd_A, B%dyn%bt_work%wd_wet_dyn)
         call check(error, ok, &
                    "wd_wet_dyn not bit-identical after restart (dropped from registry?)")
         if (allocated(error)) return
      end if
      ! Ghost refresh: the read restored owned cells only; re-establish
      ! periodic ghosts on the host before mapping (the driver does this).
      if (periodic_x) call wrap_host(gB, B)
      call sfB%init(gB)
      if (phys) call sfB%set_surface_flux_const(cfg%ocean%thermo%q_heat, 0.0_wp)
      call ocean_state_enter_data(B)
      call sfB%enter_data()

      ! Advance both one more step; pull state down and compare.  This is
      ! the acceptance gate: step (N+1) must reproduce bit-for-bit.
      call advance(gA, A, sfA)
      call advance(gB, B, sfB)

      !$acc update self(A%multilayer%h_layer, A%multilayer%u_face_x_layer, &
      !$acc&            A%multilayer%v_face_y_layer, A%barotropic%h, &
      !$acc&            A%barotropic%u_face_x, A%barotropic%v_face_y)
      !$acc update self(B%multilayer%h_layer, B%multilayer%u_face_x_layer, &
      !$acc&            B%multilayer%v_face_y_layer, B%barotropic%h, &
      !$acc&            B%barotropic%u_face_x, B%barotropic%v_face_y)

      ok = arrays_identical_3d(A%multilayer%h_layer, B%multilayer%h_layer)
      call check(error, ok, "h_layer not bit-identical after restart")
      if (allocated(error)) return
      ok = arrays_identical_3d(A%multilayer%u_face_x_layer, B%multilayer%u_face_x_layer)
      call check(error, ok, "u_face_x_layer not bit-identical after restart")
      if (allocated(error)) return
      ok = arrays_identical_3d(A%multilayer%v_face_y_layer, B%multilayer%v_face_y_layer)
      call check(error, ok, "v_face_y_layer not bit-identical after restart")
      if (allocated(error)) return
      ok = arrays_identical_2d(A%barotropic%h, B%barotropic%h)
      call check(error, ok, "barotropic h not bit-identical after restart")
      if (allocated(error)) return
      ok = arrays_identical_2d(A%barotropic%u_face_x, B%barotropic%u_face_x)
      call check(error, ok, "barotropic u_face_x not bit-identical after restart")
      if (allocated(error)) return

      do it = 1, size(A%multilayer%tracers)
         !$acc update self(A%multilayer%tracers(it)%hTr, B%multilayer%tracers(it)%hTr)
         ok = arrays_identical_3d(A%multilayer%tracers(it)%hTr, &
                                  B%multilayer%tracers(it)%hTr)
         call check(error, ok, "tracer hTr not bit-identical after restart")
         if (allocated(error)) exit
      end do
      if (allocated(error)) then
         call teardown(A, sfA)
         call teardown(B, sfB)
         call cleanup_file(FN)
         return
      end if

      ! Physics-rich: assert the checkpointed transients themselves match.
      ! A dropped registration (bl_depth / rho_layer) shows up here AND in
      ! the prognostics above, but the direct check localizes the failure.
      if (phys) then
         !$acc update self(A%vmix%bl_depth, B%vmix%bl_depth)
         !$acc update self(A%multilayer%rho_layer, B%multilayer%rho_layer)
         ok = arrays_identical_2d(A%vmix%bl_depth, B%vmix%bl_depth)
         call check(error, ok, "vmix bl_depth not bit-identical after restart")
         if (.not. allocated(error)) then
            ok = arrays_identical_3d(A%multilayer%rho_layer, B%multilayer%rho_layer)
            call check(error, ok, "rho_layer not bit-identical after restart")
         end if
      end if

      call teardown(A, sfA)
      call teardown(B, sfB)
      call cleanup_file(FN)
   end subroutine roundtrip

   subroutine test_bit_exact_closed(error)
      type(error_type), allocatable, intent(out) :: error
      call roundtrip(error, periodic_x=.false.)
   end subroutine test_bit_exact_closed

   subroutine test_bit_exact_wetdry(error)
      !! Bit-exact restart round-trip with dynamic wet/dry enabled.
      !!
      !! Geometry: 8x8x3, max_depth=0.3 m (wetdry cfg override).  A 2x2
      !! interior patch is overwritten to b=0.01 m < dry_depth=0.05 m so
      !! wd_wet_dyn seeds 0 there.  With bt_H_ref=b the bed-blocking gate
      !! keeps the patch faces closed (eta_wet ~ 0.001 m < zb_patch +
      !! dry_depth = 0.04 m), so the patch stays dry through N_STEPS=5
      !! and the hysteresis mask is genuinely nontrivial at write time
      !! (asserted).
      !!
      !! Fix-6 adversarial probe cells (mask values a naive re-seed-from-D
      !! would MISCLASSIFY, so the round-trip must carry the saved value):
      !!   * HELD-WET-IN-BAND (b=0.07 m, dry<D<rewet): A saves wd_wet_dyn=1
      !!     (history holds it wet); B cold-seeds 0 ⇒ read must restore 1.
      !!   * EMERGED (b=-0.5 m<0, h_layer floored to 2·H_VANISHED): A saves
      !!     wd_wet_dyn=0; B cold-seeds 1 ⇒ read must restore 0.
      !! Both are covered by the interior-wide wd_wet_dyn bit-exact compare.
      !!
      !! Gates:
      !!   1. Nontrivial front: >=1 interior wd_wet_dyn==0 after N_STEPS.
      !!   2. wd_wet_dyn round-trip: B%wd_wet_dyn == A%wd_wet_dyn (interior)
      !!      bit-for-bit immediately after the restart read (covers the dry
      !!      patch + both adversarial probe cells).
      !!   3. Bit-exact continuation: step (N+1) from B equals step (N+1)
      !!      from A for all prognostics (h_layer, u/v faces, tracers,
      !!      barotropic).
      !!
      !! Precision verdict: BIT-EXACT.  wd_wet_dyn holds 0.0/1.0 values
      !! (exactly representable in IEEE 754 double precision), and the
      !! restart file carries the full local double-precision array for
      !! every registered optional field.  Any drift would indicate the
      !! field was dropped from the registry or the float representation
      !! was corrupted — both are hard bugs to fix, not tolerances to loosen.
      type(error_type), allocatable, intent(out) :: error
      call roundtrip(error, periodic_x=.false., wetdry=.true.)
   end subroutine test_bit_exact_wetdry

   subroutine test_bit_exact_physics(error)
      type(error_type), allocatable, intent(out) :: error
      call roundtrip(error, periodic_x=.false., physics_rich=.true.)
   end subroutine test_bit_exact_physics

   subroutine test_bit_exact_periodic(error)
      type(error_type), allocatable, intent(out) :: error
      call roundtrip(error, periodic_x=.true.)
   end subroutine test_bit_exact_periodic

   subroutine test_decomp_mismatch(error)
      !! Write a restart with a 1x1 decomp, then check it against a
      !! 2x1 decomp — must report a mismatch (ierr /= 0).
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(hgrid_t) :: g
      type(ocean_state_t) :: A
      type(ocean_surface_flux_t) :: sf
      type(decomp_t) :: d_write, d_other
      character(len=*), parameter :: FN = "test_ocean_restart_mismatch.nc"
      integer :: ierr

      call make_cfg(cfg)
      call decomp_init(d_write, NX, NY, 1, 1, 0)
      ! A 2x1 decomposition of the same global grid: different px / nx_local.
      call decomp_init(d_other, NX, NY, 2, 1, 0)
      call cleanup_file(FN)

      call setup_state(cfg, g, A, sf, periodic_x=.false.)
      call ocean_state_restart_write(A, g, d_write, FN, 0.0_wp, 0)
      call teardown(A, sf)

      call ocean_restart_check_decomp(FN, d_other, ierr)
      call check(error, ierr /= 0, &
                 "decomp mismatch should have been detected (ierr /= 0)")
      call cleanup_file(FN)
   end subroutine test_decomp_mismatch

   subroutine test_decomp_mismatch_read(error)
      !! Route a decomp mismatch through the PRODUCTION read wrapper
      !! (`ocean_state_restart_read` with `ierr` present): it must return
      !! ierr /= 0 rather than error-stopping (review #9b / #4).
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(hgrid_t) :: gA, gB
      type(ocean_state_t) :: A, B
      type(ocean_surface_flux_t) :: sfA
      type(decomp_t) :: d_write, d_other
      character(len=*), parameter :: FN = "test_ocean_restart_mismatch_rd.nc"
      real(wp) :: t_read
      integer :: step_read, ierr

      call make_cfg(cfg)
      call decomp_init(d_write, NX, NY, 1, 1, 0)
      call decomp_init(d_other, NX, NY, 2, 1, 0)
      call cleanup_file(FN)

      ! Write with a 1x1 decomp.
      call setup_state(cfg, gA, A, sfA, periodic_x=.false.)
      call ocean_state_restart_write(A, gA, d_write, FN, 0.0_wp, 0)
      call teardown(A, sfA)

      ! Fresh host state, read with the mismatching 2x1 decomp + ierr.
      call gB%init(NX, NY, NGHOST, DX, DX)
      B%multilayer%nz_ml = NZ
      call B%init(gB)
      call metrics_fill_cartesian(B%metrics, gB, gB%dx, gB%dy)
      call metrics_finalize(B%metrics)
      call ocean_state_seed_from_cfg(B, gB, cfg)
      call register_default_tracers( &
         B%multilayer%tracers(B%multilayer%idx_salinity), &
         B%multilayer%tracers(B%multilayer%idx_temperature), cfg)
      ierr = 0
      call ocean_state_restart_read(B, gB, d_other, FN, t_read, step_read, ierr=ierr)
      call check(error, ierr /= 0, &
                 "read wrapper should return ierr /= 0 on decomp mismatch")
      call B%destroy()
      call cleanup_file(FN)
   end subroutine test_decomp_mismatch_read

   logical function arrays_identical_2d(a, b) result(same)
      !! Compare OWNED (interior) cells only.  The restart restores owned
      !! cells; ghosts re-establish from BC/wrap each step, so an
      !! interior match is the bit-exact gate (ROADMAP A1).
      real(wp), intent(in) :: a(:, :), b(:, :)
      integer :: i0, i1, j0, j1
      i0 = NGHOST + 1
      i1 = NGHOST + NX
      j0 = NGHOST + 1
      j1 = NGHOST + NY
      same = all(a(i0:i1, j0:j1) == b(i0:i1, j0:j1))
   end function arrays_identical_2d

   logical function arrays_identical_3d(a, b) result(same)
      real(wp), intent(in) :: a(:, :, :), b(:, :, :)
      integer :: i0, i1, j0, j1
      i0 = NGHOST + 1
      i1 = NGHOST + NX
      j0 = NGHOST + 1
      j1 = NGHOST + NY
      same = all(a(i0:i1, j0:j1, :) == b(i0:i1, j0:j1, :))
   end function arrays_identical_3d

   subroutine cleanup_file(fn)
      character(len=*), intent(in) :: fn
      integer :: u, ios
      character(len=256) :: msg
      logical :: exists
      inquire (file=fn, exist=exists)
      if (exists) then
         open (newunit=u, file=fn, status="old", action="readwrite", &
               iostat=ios, iomsg=msg)
         if (ios == 0) close (u, status="delete")
      end if
   end subroutine cleanup_file

end module test_ocean_restart
