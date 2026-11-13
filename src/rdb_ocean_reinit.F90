!! Re-derive the setup-time snapshots of a LIVE model after its initial state
!! was overwritten through the C API (`rdb_ocean_set_h/u/v/tracer`).
module rdb_ocean_reinit
   !! `engine_setup` builds the built-in initial condition and THEN takes a
   !! series of snapshots and derivations of it (sponge reference, pseudo-salt
   !! seed, land-face masking, seam ghosts, `rho_layer`).  A caller that
   !! overwrites `h`/`u`/`v`/tracers on the finished handle leaves all of them
   !! describing the OLD initial state.  `ocean_rebaseline_apply` redoes the
   !! ones that matter physically, in the order `engine_setup` runs them, on
   !! the current prognostics.  (The preferred route is the staged path,
   !! `rdb_ocean_ic_inject`, which gets this right by construction.)
   !!
   !! ## Device discipline (`-gpu=...,mem:separate`)
   !! The caller has already drained the windowed tracer accumulator and
   !! refreshed the HOST copies of `h`, `u`, `v` and every `c` from the
   !! device, so the host is authoritative on entry.  Every host write below
   !! is followed by an `!$acc update device(<leaf array>)`; every value a
   !! device kernel produced (`rho_layer`, `bt_eta`, the seam-wrapped state)
   !! is pulled with `!$acc update self(<leaf array>)` before the host reads
   !! it.  Never an aggregate: only leaf component arrays are named.  All of
   !! it is inert on a host-only build.
   !!
   !! ## Known limit
   !! The tidal-mixing static fields are folded with `depth = sum(h_layer)` at
   !! setup and are NOT re-folded here (the error is O(eta/H) after a `set_h`).
   !! The wet/dry hysteresis mask is likewise left to the first step.
   use rdb_grid, only: hgrid_t
   use rdb_config, only: config_t
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   use rdb_ocean_balance, only: ocean_balance_velocity
   use rdb_ocean_state, only: ocean_state_t, ocean_state_seed_land_cells
   use rdb_ocean_ic_inject, only: ocean_ic_fill_state_ghosts
   use rdb_ocean_eos_compute, only: ocean_eos_compute
   use rdb_ocean_pseudo_salt, only: ocean_pseudo_salt_seed
   use rdb_ocean_sponge, only: ocean_sponge_snapshot_reference
   use rdb_ocean_halo, only: ocean_halo_is_decomposed_x, ocean_halo_is_decomposed_y
   use rdb_ocean_halo_state, only: ocean_halo_exchange_ml_state
   use rdb_ocean_periodic, only: ocean_periodic_wrap_state
   use rdb_ocean_fold_apply, only: ocean_fold_wrap_state
   use rdb_barotropic_coupling, only: derive_bt_from_layers
   implicit none
   private

   public :: ocean_rebaseline_apply
   public :: ocean_balance_apply

contains

   subroutine ocean_rebaseline_apply(state, grid, thermo_active, reseed_pseudo_salt)
      !! Steps 3-8 of `rdb_ocean_rebaseline` (the API wrapper does the
      !! precondition checks, the windowed-tracer flush and the host refresh
      !! first).  Decomposition-safe (exercised decomposed by the `balance_*`
      !! cases of `test_ocean_decomp_bitid_mpi`); the C API handle that reaches
      !! it is single-rank.
      type(ocean_state_t), intent(inout) :: state
      type(hgrid_t), intent(in) :: grid
      logical, intent(in) :: thermo_active
         !! `&ocean_thermo_nml enable_thermodynamics`: gates the EOS exactly
         !! as the step does, so a prescribed `layer_rho_init` survives.
      logical, intent(in) :: reseed_pseudo_salt
         !! `.false.` when the caller set pseudo-salt itself.

      ! ---- 3. seam ghosts ----
      call refresh_seam_ghosts(state, grid)

      ! ---- 4. land / closed faces (host) ----
      ! Zeroes u/v/hu/hv on closed faces, holds land tracer content at 0 and
      ! floors land h: a nonzero user u on a closed face would otherwise feed
      ! the enstrophy Coriolis, the viscosity and the fast-loop reference.
      call ocean_state_seed_land_cells(state, grid)

      ! ---- 6a. pseudo-salt follows S (host; whole array, ghosts included) ----
      if (reseed_pseudo_salt) call ocean_pseudo_salt_seed(state%multilayer)

      ! One push of everything the host steps above may have written.
      call push_state(state)
      associate (ms => state%multilayer)
         ! Host hu/hv are current: step 2 pulled them, step 4 only zeroed closed faces.
         !$acc update device(ms%hu_face_x_layer, ms%hv_face_y_layer)
      end associate

      ! ---- 5. rho_layer from the final T/S/h (device kernel), host mirror ----
      call ocean_eos_compute(state%eos, state%multilayer, active=thermo_active)
      associate (ms => state%multilayer)
         !$acc update self(ms%rho_layer)
      end associate

      ! ---- 7. sponge reference snapshot (host), then its device copy ----
      if (state%sponge%enable .and. state%sponge%is_init) then
         call ocean_sponge_snapshot_reference(state%sponge, grid, state%multilayer)
         associate (sp => state%sponge)
            !$acc update device(sp%ref_tracer, sp%u_ref, sp%v_ref)
         end associate
      end if

      ! ---- 8. barotropic state, so the eta getter is current ----
      if (state%dyn%bt_work%is_init) then
         call derive_bt_from_layers(grid, state%dyn%bt_work, state%multilayer, state%metrics)
         associate (bw => state%dyn%bt_work)
            !$acc update self(bw%bt_eta)
         end associate
      end if
   end subroutine ocean_rebaseline_apply

   subroutine ocean_balance_apply(state, grid, cfg, mode, ref_level, thermo_active, &
                                  reseed_pseudo_salt, ierr)
      !! `ocean_balance_velocity` (device-resident state, host authoritative on
      !! entry) followed by `ocean_rebaseline_apply`, so the sponge reference,
      !! land mask, pseudo-salt and eta describe the balanced state.  On a
      !! refusal `ierr /= OCEAN_STATUS_OK`, the error ring carries the reason
      !! and nothing was modified.
      type(ocean_state_t), intent(inout) :: state
      type(hgrid_t), intent(in) :: grid
      type(config_t), intent(in) :: cfg
      integer, intent(in) :: mode
         !! `BALANCE_*`.
      character(len=*), intent(in) :: ref_level
         !! `"keep_h"` or `"bottom_rest"`.
      logical, intent(in) :: thermo_active
      logical, intent(in) :: reseed_pseudo_salt
      integer, intent(out) :: ierr

      ! The balance stencil reads ghost thickness/density at the seams: a live
      ! caller has only written the interior.
      call refresh_seam_ghosts(state, grid)
      call ocean_balance_velocity(grid, cfg, state, mode, ierr, ref_level=ref_level)
      if (ierr /= OCEAN_STATUS_OK) return
      call ocean_rebaseline_apply(state, grid, thermo_active, reseed_pseudo_salt)
   end subroutine ocean_balance_apply

   subroutine refresh_seam_ghosts(state, grid)
      !! Make the host AND device ghost cells of `h`, `u`, `v`, `c` consistent
      !! with their (host-authoritative) interior.
      type(ocean_state_t), intent(inout) :: state
      type(hgrid_t), intent(in) :: grid

      integer :: it

      ! Host: periodic image / zero-gradient (the staged path's fill), so the
      ! ghosts of an interior-only overwrite are consistent with the new
      ! interior.  Push, then run the SAME exchange + periodic wrap + north
      ! fold the step runs at stage entry (device kernels), and pull the
      ! result back so the host copy carries the seam-correct ghosts for the
      ! host steps below.
      call ocean_ic_fill_state_ghosts(state, grid, state%bc%periodic_x, state%bc%periodic_y)
      call push_state(state)
      call ocean_halo_exchange_ml_state(state%multilayer)
      call ocean_periodic_wrap_state(grid, state%bc, state%multilayer, &
                                     skip_x=ocean_halo_is_decomposed_x(), &
                                     skip_y=ocean_halo_is_decomposed_y())
      call ocean_fold_wrap_state(grid, state%bc, state%multilayer)
      associate (ms => state%multilayer)
         !$acc update self(ms%h_layer, ms%u_face_x_layer, ms%v_face_y_layer)
         if (allocated(ms%tracers)) then
            do it = 1, size(ms%tracers)
               if (.not. allocated(ms%tracers(it)%c)) cycle
               !$acc update self(ms%tracers(it)%c)
            end do
         end if
      end associate
   end subroutine refresh_seam_ghosts

   subroutine push_state(state)
      !! Host -> device for the prognostic leaves this module writes: `h`,
      !! `u`, `v` and every registered tracer's `c`.  Leaf arrays only.
      type(ocean_state_t), intent(inout) :: state
      integer :: it

      associate (ms => state%multilayer)
         !$acc update device(ms%h_layer, ms%u_face_x_layer, ms%v_face_y_layer)
         if (allocated(ms%tracers)) then
            do it = 1, size(ms%tracers)
               if (.not. allocated(ms%tracers(it)%c)) cycle
               !$acc update device(ms%tracers(it)%c)
            end do
         end if
      end associate
   end subroutine push_state

end module rdb_ocean_reinit
