!! Stage hooks of the CICE-port ice step bound to the engine's state (roundabout-original).
module rdb_ice_cice_hooks
   !! The production implementations of the ten hook stages of `ice_step_cice`
   !! (`rdb_ice_step_cice`), each a thin adapter from the engine's state to an existing
   !! driver, plus the entry `ice_cice_step` the engine calls once per outer step and the
   !! cold-start seeder `ice_cice_ic_apply`.  Nothing here computes physics except the
   !! small copy / wrap / bit kernels named below; no test double is installed.
   !!
   !! **Binding.**  A hook is `subroutine(ctx, status)` and carries no engine handle, so the
   !! module holds pointers to the engine's state and grid (`ice_cice_bind`), re-bound at the
   !! top of EVERY step by the engine (the engine dummy has the `target` attribute only for
   !! the duration of a call).  The ice slots themselves arrive through `ctx`.
   !!
   !! **Stage by stage** (contract text of `rdb_ice_step_cice` honoured one to one):
   !!  1. `ocean_import`: the melt-potential provider (`ice_melt_potential`, depth `hfreeze`,
   !!     the boundary-layer depth of EPBL or KPP when one is active) is called EVERY ice step
   !!     with the ice `dt` on the ocean-side frazil bank `ice%frazil_heat` (the SIS2 thermo-
   !!     cadence provider call and `ice_frazil_uptake` never run under this scheme: the bank is
   !!     consumed once, by new ice).  `ice_import_ocean` follows with the OPEN averager window
   !!     (the melt cap).  `work%sst, Tf, frzmlt` are copies of the import; `work%bank_ocean`
   !!     <- the ocean bank, `work%bank` <- `bank_ocean` (spent by the new-ice stage).
   !!  2. `velocity_to_T` (`highfreq` only): `ice_atm_velocity_T` on the face velocities.
   !!  3. `atm_prepare`: `ice_atm_prepare_cell_inputs` from the prescribed atmospheric state
   !!     (`atm%enable`) or the uniform atmosphere of the run slot; the reference heights;
   !!     `work%sw_down` is the surface downward shortwave.
   !!  9. and 15. `halo_state`, `halo_state_post`: `ice_halo_category_state` (the ocean's exchange, the
   !!     periodic wrap and the tripolar fold of the category state and `aice0`; a wall edge keeps its
   !!     ghosts: the C-grid metrics pin them to land).  The scheme runs on several ranks on walls and
   !!     periodic edges and across the tripolar fold (bitwise: `cice_ice*` of `test_ocean_decomp_bitid_mpi`,
   !!     the fold as a symmetry in `test_cice_tripole`; the fold of a scalar category field is a copy).
   !! 11. `stress_to_faces`: the T wind stress (`strair_x_T/y_T`, written on owned cells by the
   !!     thermodynamics) gets its ghost ring (`ice_halo_centre_pair`: exchange on a split axis, the
   !!     local wrap otherwise, the fold on the fold row with the VECTOR sign: both components negate;
   !!     2 exchanges), then
   !!     `ice_atm_stress_to_faces` averages it onto the dynamics slot's `strairxE/yN` (the dynamics
   !!     keeps them: `wind_on_faces`).
   !!
   !! **Decomposition contract (what each exchange feeds).**  Everything the thermodynamics, ITD and
   !! ridging compute runs on OWNED cells from owned inputs, so a rank needs no halo for them.  What
   !! a later stage reads across a seam is refreshed before it: the category state ring (`halo_state`,
   !! `halo_state_post`: dynamics strength and transport read it), the T-cell stress (`stress_to_faces`),
   !! the face velocities (the EVP's own per-subcycle exchange).  The ghosts of the assembled `Q_heat`/`Q_salt`
   !! (and `sw_net`) that KPP and EPBL read are refreshed for every path by `ocean_seam_refresh_surface_flux`
   !! (`engine_step_finalize`), not here.  Every decision (a hook's
   !! status, the thermodynamics' verdict, the error word) is taken on a reduced value
   !! (`ice_transport_seam_bind_reductions`, `ice_err_bind_comm`), so every rank stops together.
   !! `ice_mr_test_skip` is the test seam of the negative controls: each exchange above (stress, state
   !! halo) is skipped in turn and the decomposition test must see the difference.
   !! 12. `dynamics`: the scheme-neutral aggregates are refreshed from the haloed state, then the
   !!     C-grid EVP of `ice_dyn_c_step` runs with the CATEGORY strength (`kstrength`) from the
   !!     post-thermodynamics state; `work%tau_x/y` (the persistent centre ice-ocean stress) and
   !!     `work%rdg_conv/rdg_shear` are copied from the dynamics slot.
   !! 13. `transport`: `ice_step_transport` with the post mode of the context (NONE); the transport's
   !!     decomposition seam (`rdb_ice_transport_seam`: reductions, fold-row flux pairing) is bound at
   !!     every call, a no-op on one rank without a fold.
   !! 14. `ridging`: `ice_step_ridge`; cells that needed more than one pass get the counted
   !!     `ICE_ERR_RIDGE_PASSES` bit.
   !! 17. `ocean_handoff`: accumulate the step's fluxes into the window mean; the face stress
   !!     blend (`ice_ocean_stress_flux_c`, its own restart mirror); at the ocean thermo step the
   !!     adapters apply the window mean (`ice_cice_apply_window`, virtual-salt form, with the seam
   !!     halo exchange) and reset the window; the spent bank is written back to the ocean-side
   !!     bank; the scheme-neutral aggregates are refreshed for every consumer (diagnostics,
   !!     console totals, the open-water weight of the bulk forcing).
   !!
   !! **Atmosphere and stress.**  With `&ocean_atm_forcing_nml enable` the `atm_prepare` hook builds
   !! the ice's inputs from the SAME `atm_state_t` arrays and time level the bulk forcing used at
   !! the top of this outer step (`tas, huss, uas, vas, psl, rlds, prra, prsn`, `rsds` for the ice
   !! radiation, the one wind rotation); otherwise from the uniform atmosphere of the run slot.  The
   !! ocean-side stress blend (`ice_ocean_stress_flux_c`) reads `ice%tau_a_x/y`: under the prescribed
   !! atmosphere the bulk forcing REFRESHES them every step with the unweighted wind stress of the
   !! step (and applies the blend itself at the top of the step, with the concentration and the ice
   !! stress of the previous ice step); with the constant wind of `&physics_nml` they are the
   !! configure-time snapshot.  Either way the hand-off blend of this stage uses the wind the ice
   !! just felt.  The `ice_cice_stress_blend` of the adapters (a different formulation: per-ice-area
   !! stress) is not used, the face blend of the C-grid dynamics conserves the momentum of the pair
   !! exactly.  Rain through the ice (`rain_through_impl`) is removed from the window mean of `fresh`
   !! before the hand-over (unless `&ocean_forcing_nml freshwater = "virtual"`, where the open
   !! fraction's rain dilutes too and the term stays): see the module docstring of `rdb_ice_cice_adapters`.
   use, intrinsic :: iso_fortran_env, only: int32
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_config, only: config_t
   use rdb_ocean_state, only: ocean_state_t
   use rdb_ocean_halo, only: ocean_halo_nranks
   use rdb_ice_state, only: ice_aggregates_fill
   use rdb_ice_cice_state, only: ice_cice_state_t, ice_update_state, ice_state_seed_uniform
   use rdb_ice_step_cice, only: ice_step_cice, ice_step_hooks_t, ice_step_opts_t, ice_step_result_t, &
                                ice_step_ctx_t, ice_stage_name, ICE_STEP_OK, ICE_STEP_FATAL, &
                                ICE_STAGE_OCEAN_IMPORT, ICE_STAGE_VELOCITY_TO_T, ICE_STAGE_ATM_PREPARE, &
                                ICE_STAGE_HALO_STATE, ICE_STAGE_STRESS_TO_FACES, ICE_STAGE_DYNAMICS, &
                                ICE_STAGE_TRANSPORT, ICE_STAGE_RIDGING, ICE_STAGE_HALO_STATE_POST, &
                                ICE_STAGE_OCEAN_HANDOFF
   use rdb_ice_errors, only: ice_err_enforce, ice_err_set_combine, ice_err_combine_iface, ICE_ERR_RIDGE_PASSES, &
                             ICE_ERR_FATAL_MASK
   use rdb_ice_melt_potential, only: ice_melt_potential
   use rdb_ice_import_ocean, only: ice_import_ocean, ICE_IMPORT_OK
   use rdb_ice_atm_inputs, only: ice_atm_prepare_cell_inputs, ice_atm_heights, ice_atm_velocity_T, &
                                 ice_atm_stress_to_faces
   use rdb_ice_dyn_c, only: ice_dyn_c_step, ice_dyn_c_cat_t
   use rdb_ice_transport_cice, only: ice_step_transport, transport_result_t, TRANSPORT_CFG_OK
   use rdb_ice_transport_seam, only: seam_max_one
   use rdb_ice_transport_seam_comm, only: ice_transport_seam_bind, ice_halo_category_state, ice_halo_centre_pair, &
                                          ice_transport_seam_bind_reductions, ice_err_bind_comm
   use rdb_ice_cice_run, only: ICE_LCUM_NCLIP, ICE_LCUM_AREA, ICE_LCUM_VICE, ICE_LCUM_VSNO, ICE_LCUM_EICE, ICE_LCUM_ESNO
   use rdb_ice_ridge, only: ice_step_ridge
   use rdb_ice_flux_avg, only: ice_flux_avg_accumulate, ice_flux_avg_stash_pending, ICE_AVG_OK
   use rdb_ice_cice_adapters, only: ice_cice_apply_window, ICE_ADAPT_OK
   use rdb_ice_ocean_coupler, only: ice_ocean_stress_flux_c
   use rdb_ice_init, only: ice_ic_params_t, ICE_IC_CONC_UNIFORM, ICE_IC_CONC_LATITUDES
   use pic_logger, only: logger => global_logger
   use pic_strings, only: to_string
   implicit none
   private

   public :: ice_cice_bind, ice_cice_step, ice_cice_step_result, ice_cice_ic_apply, ice_cice_post_enter, ice_cice_ledger_report
   public :: ice_cice_unbind, ice_cice_hooks_manifest

   integer, save, public :: ice_mr_test_skip = 0
      !! TEST SEAM of the decomposition tests' negative controls, `0` in every run: bit 1 skips the exchange of the
      !! T-cell wind stress, bit 4 the category-state halo, bit 8 folds that stress pair as a SCALAR (the wrong sign), bit 16
      !! skips its fold, bit 32 switches the transport's fold-row flux pairing off.  A
      !! decomposed run with one bit set must differ from the serial one (the test proves each exchange is needed).
   integer, parameter, public :: ICE_SKIP_STRESS_HALO = 1, ICE_SKIP_STATE_HALO = 4, &
                                 ICE_SKIP_FOLD_SIGN = 8, ICE_SKIP_FOLD_PAIR = 16, ICE_SKIP_FOLD_FLUX = 32

   type(ocean_state_t), pointer, save :: S => null()
   type(hgrid_t), pointer, save :: G => null()
   type(ice_step_hooks_t), save :: HOOKS
   logical, save :: hooks_ready = .false.

contains

   ! ------------------------------------------------------------------
   ! entry points
   ! ------------------------------------------------------------------

   subroutine ice_cice_bind(state, grid)
      !! Point the hooks at the engine's state and grid (every step).
      type(ocean_state_t), intent(inout), target :: state
      type(hgrid_t), intent(in), target :: grid

      S => state
      G => grid
      if (.not. hooks_ready) then
         call HOOKS%reset()
         call HOOKS%install(ICE_STAGE_OCEAN_IMPORT, hook_ocean_import)
         call HOOKS%install(ICE_STAGE_VELOCITY_TO_T, hook_velocity_to_T)
         call HOOKS%install(ICE_STAGE_ATM_PREPARE, hook_atm_prepare)
         call HOOKS%install(ICE_STAGE_HALO_STATE, hook_halo_state)
         call HOOKS%install(ICE_STAGE_STRESS_TO_FACES, hook_stress_to_faces)
         call HOOKS%install(ICE_STAGE_DYNAMICS, hook_dynamics)
         call HOOKS%install(ICE_STAGE_TRANSPORT, hook_transport)
         call HOOKS%install(ICE_STAGE_RIDGING, hook_ridging)
         call HOOKS%install(ICE_STAGE_HALO_STATE_POST, hook_halo_state)
         call HOOKS%install(ICE_STAGE_OCEAN_HANDOFF, hook_ocean_handoff)
         hooks_ready = .true.
      end if
   end subroutine ice_cice_bind

   subroutine ice_cice_unbind()
      !! Release what `ice_cice_bind` and the step bound: the state and grid pointers and the cross-rank error
      !! combination (a process-global procedure pointer: left bound it would make the NEXT engine's error
      !! reductions collective).  Called by `engine_exit_data`.
      procedure(ice_err_combine_iface), pointer :: p

      p => null()
      call ice_err_set_combine(p)
      S => null()
      G => null()
   end subroutine ice_cice_unbind

   subroutine ice_cice_hooks_manifest(installed)
      !! Which hook slots hold an implementation (for the tests): `installed(stage)`.
      logical, intent(out) :: installed(:)
      integer :: k

      do k = 1, size(installed)
         installed(k) = HOOKS%installed(k)
      end do
   end subroutine ice_cice_hooks_manifest

   subroutine ice_cice_step_result(dt, thermo_step, istep, res)
      !! One ice step, returning the sequencer's result instead of stopping on a failure (`ice_cice_step` is the
      !! production entry that stops).  On a decomposed run every decision inside is taken on a reduced value, so
      !! `res%status`, `res%stage` and `res%summary` are the same on every rank.
      real(wp), intent(in) :: dt
      logical, intent(in) :: thermo_step
      integer, intent(in) :: istep
      type(ice_step_result_t), intent(out) :: res
      type(ice_step_opts_t) :: opts
      integer :: nlen, nstat

      opts%ndtd = 1
      opts%nranks = ocean_halo_nranks()
      ! the decomposition seam of the step: every decision (a hook's status, the thermodynamics' failure, the
      ! error word) is taken on a reduced value, so every rank stops together; pointer assignments, no
      ! collective, and nothing on one rank
      call ice_transport_seam_bind_reductions(opts%seam, G, S%bc)
      call ice_err_bind_comm()
      opts%dt_outer = dt
      opts%fatal_mask = ICE_ERR_FATAL_MASK
      if (S%icec%skip_thermo) then
         ! a test-harness diagnostic (it bypasses the stage manifest): only with the harness's environment switch
         call get_environment_variable("RDB_ICE_ALLOW_SKIP_THERMO", length=nlen, status=nstat)
         if (nstat /= 0 .or. nlen < 1) then
            error stop "ice_cice_step: skip_thermo is a test-harness diagnostic (set RDB_ICE_ALLOW_SKIP_THERMO)"
         end if
      end if
      opts%skip_thermo = S%icec%skip_thermo
      if (S%icec%ledger%on) opts%ledger => S%icec%ledger
      call ice_step_cice(S%ice%cice, S%icec%flux, S%icec%rad, S%icec%work, HOOKS, opts, dt, thermo_step, istep, res)
      S%icec%n_steps = S%icec%n_steps + 1
   end subroutine ice_cice_step_result

   subroutine ice_cice_step(dt, thermo_step, istep)
      !! One ice step.  A refusal or a stage failure stops the run with the stage named; a fatal cell error bit
      !! aborts through `ice_err_enforce` (first offender dumped) on EVERY rank: the summary is the combined one.
      real(wp), intent(in) :: dt
      logical, intent(in) :: thermo_step
      integer, intent(in) :: istep
      type(ice_step_result_t) :: res

      call ice_cice_step_result(dt, thermo_step, istep, res)
      if (res%status == ICE_STEP_OK) return
      if (res%status == ICE_STEP_FATAL) call ice_err_enforce(res%summary, S%ice%cice%nx_global, ICE_ERR_FATAL_MASK)
      call logger%error("ice_step_cice: step "//to_string(istep)//" stopped, status "//to_string(res%status)// &
                        ", stage '"//trim(ice_stage_name(res%stage))//"' (raw status "// &
                        to_string(res%stage_status)//"): "//trim(res%message))
      error stop "ice_step_cice: the ice step failed (see the log)"
   end subroutine ice_cice_step

   subroutine ice_cice_ledger_report(state, compute_rank)
      !! The console line of the named-term ice ledger (`&ocean_ice_nml ledger`): the cumulative
      !! residual of each line (water, salt, heat) relative to its cumulative scale, with the stage
      !! that carries the largest one.  Called at the status cadence next to the ocean's console
      !! block; the sums are EFP reproducing sums, so the numbers do not depend on the backend.  Inert
      !! (prints nothing, touches nothing) under `scheme = "sis2"` and with the ledger off.  Single
      !! rank: the ledger refuses more, so no collective is needed.
      type(ocean_state_t), intent(in) :: state
      integer, intent(in) :: compute_rank

      if (compute_rank /= 0) return
      if (.not. state%icec%is_init) return
      if (.not. state%icec%ledger%on .or. state%icec%ledger%nsteps < 0.5_wp) return
      call logger%info(trim(state%icec%ledger%console_line()))
   end subroutine ice_cice_ledger_report

   subroutine ice_cice_post_enter(state)
      !! After `enter_data` (device): the aggregates are not restart state, so rebuild them over
      !! the ghosted extent and publish them to the scheme-neutral slot.
      type(ocean_state_t), intent(inout) :: state

      if (.not. state%ice%cice%is_init) return
      call ice_update_state(state%ice%cice, .true.)
      call ice_aggregates_fill(state%ice)
   end subroutine ice_cice_post_enter

   subroutine ice_cice_ic_apply(state, par, hm, geolatT, ierr)
      !! Cold-start seed from `&ocean_ice_ic_nml` (host, BEFORE `enter_data`): the analytic
      !! concentration of the SIS2 seeder (`uniform`, or the `latitudes` polar caps) as ONE
      !! slab per cell (thickness `h_ice`, snow `h_snow`, surface temperature `t_ice`; the
      !! linear profile to the freezing point `-depressT*S_ref` of the EOS reference salinity)
      !! in the category that holds it; land and open-water cells empty.  `zero`: empty ocean.
      type(ocean_state_t), intent(inout), target :: state
      type(ice_ic_params_t), intent(in) :: par
      real(wp), intent(in) :: hm(:, :), geolatT(:, :)
      integer, intent(out) :: ierr
      type(ice_cice_state_t), pointer :: c
      real(wp) :: tf, aice_in
      integer :: i, j, n

      c => state%ice%cice
      tf = -c%p%c%depressT*state%eos%S_ref
      select case (par%conc_config)
      case (ICE_IC_CONC_UNIFORM)
         aice_in = par%conc
      case (ICE_IC_CONC_LATITUDES)
         aice_in = 1.0_wp
      case default
         call ice_state_seed_uniform(c, 0.0_wp, 0.0_wp, 0.0_wp, min(tf, c%p%c%Timelt), tf, ierr)
         return
      end select
      call ice_state_seed_uniform(c, aice_in, par%h_ice, par%h_snow, min(par%t_ice, c%p%c%Timelt), tf, ierr)
      if (ierr /= 0) return
      do j = 1, c%ny
         do i = 1, c%nx
            if (hm(i, j) > 0.5_wp) then
               if (par%conc_config /= ICE_IC_CONC_LATITUDES) cycle
               if (geolatT(i, j) > par%arctic_edge .or. geolatT(i, j) < par%antarctic_edge) cycle
            end if
            do n = 1, c%ncat
               c%aicen(i, j, n) = 0.0_wp
               c%vicen(i, j, n) = 0.0_wp
               c%vsnon(i, j, n) = 0.0_wp
               c%Tsfcn(i, j, n) = tf
            end do
            c%aice(i, j) = 0.0_wp
            c%vice(i, j) = 0.0_wp
            c%vsno(i, j) = 0.0_wp
            c%aice0(i, j) = 1.0_wp
         end do
      end do
   end subroutine ice_cice_ic_apply

   ! ------------------------------------------------------------------
   ! small flat kernels
   ! ------------------------------------------------------------------

   subroutine copy2_impl(nx, ny, src, dst)
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: src(nx, ny)
      real(wp), intent(inout) :: dst(nx, ny)
      integer :: i, j

      do concurrent(j=1:ny, i=1:nx)
         dst(i, j) = src(i, j)
      end do
   end subroutine copy2_impl

   subroutine import_copy_impl(nx, ny, sst_i, tf_i, frz_i, bank_o, sst, tf, frz, bank_oc, bank)
      !! The ocean import into the work slot, and the frazil bank into the spendable bank.
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: sst_i(nx, ny), tf_i(nx, ny), frz_i(nx, ny), bank_o(nx, ny)
      real(wp), intent(inout) :: sst(nx, ny), tf(nx, ny), frz(nx, ny), bank_oc(nx, ny), bank(nx, ny)
      integer :: i, j

      do concurrent(j=1:ny, i=1:nx)
         sst(i, j) = sst_i(i, j)
         tf(i, j) = tf_i(i, j)
         frz(i, j) = frz_i(i, j)
         bank_oc(i, j) = bank_o(i, j)
         bank(i, j) = bank_o(i, j)
      end do
   end subroutine import_copy_impl

   subroutine rain_through_impl(nx, ny, is, ie, js, je, w, frain, aice_init, acc_fresh)
      !! Take this step's rain-through-the-ice term `frain*aice_init` out of the running window mean
      !! of `fresh` (the mean moved by `w = dt/t_acc` of this step's value).
      integer, intent(in) :: nx, ny, is, ie, js, je
      real(wp), intent(in) :: w
      real(wp), intent(in) :: frain(nx, ny), aice_init(nx, ny)
      real(wp), intent(inout) :: acc_fresh(nx, ny)
      integer :: i, j

      do concurrent(j=js:je, i=is:ie)
         acc_fresh(i, j) = acc_fresh(i, j) - w*(frain(i, j)*aice_init(i, j))
      end do
   end subroutine rain_through_impl

   subroutine ridge_bits_impl(nx, ny, is, ie, js, je, niter, ice_err)
      !! Counted bit for the cells that needed more than one ridging pass.
      integer, intent(in) :: nx, ny, is, ie, js, je
      integer(int32), intent(in) :: niter(nx, ny)
      integer(int32), intent(inout) :: ice_err(nx, ny)
      integer :: i, j

      do concurrent(j=js:je, i=is:ie)
         if (niter(i, j) > 1_int32) ice_err(i, j) = ior(ice_err(i, j), ICE_ERR_RIDGE_PASSES)
      end do
   end subroutine ridge_bits_impl

   subroutine wrap_category_state(c)
      !! The halo refresh of the category state and `aice0` (module docstring): the ocean's two-pass
      !! exchange (the periodic wrap on an axis that is not split; a wall edge keeps its ghosts) and the
      !! tripolar fold, through `ice_halo_category_state` (`rdb_ice_transport_seam_comm`).
      type(ice_cice_state_t), intent(inout) :: c

      call ice_halo_category_state(G, S%bc, c%nx, c%ny, c%ncat, c%nilyr, c%nslyr, c%aicen, c%vicen, c%vsnon, &
                                   c%Tsfcn, c%qicen, c%sicen, c%qsnon, c%aice0)
   end subroutine wrap_category_state

   ! ------------------------------------------------------------------
   ! hooks
   ! ------------------------------------------------------------------

   subroutine hook_ocean_import(ctx, status)
      type(ice_step_ctx_t), intent(inout) :: ctx
      integer, intent(out) :: status
      integer :: ierr

      status = 0
      if (S%epbl%enable .and. allocated(S%epbl%mld)) then
         call ice_melt_potential(G, S%eos, S%multilayer, S%ice%frazil_heat, S%icec%hfreeze, ctx%dt, &
                                 S%ice%melt_potential, S%ice%frzmlt, h_bl=S%epbl%mld)
      else if (S%vmix%use_kpp .and. allocated(S%vmix%bl_depth)) then
         call ice_melt_potential(G, S%eos, S%multilayer, S%ice%frazil_heat, S%icec%hfreeze, ctx%dt, &
                                 S%ice%melt_potential, S%ice%frzmlt, h_bl=S%vmix%bl_depth)
      else
         call ice_melt_potential(G, S%eos, S%multilayer, S%ice%frazil_heat, S%icec%hfreeze, ctx%dt, &
                                 S%ice%melt_potential, S%ice%frzmlt)
      end if
      call ice_import_ocean(G, S%eos, S%multilayer, S%ice%frzmlt, ctx%cice, S%icec%oin, ierr, &
                            S%icec%avg, ctx%dt, pend_applied=ctx%thermo_step)
      if (ierr /= ICE_IMPORT_OK) then
         status = ierr
         return
      end if
      call import_copy_impl(ctx%work%nx, ctx%work%ny, S%icec%oin%sst, S%icec%oin%Tf, S%icec%oin%frzmlt, &
                            S%ice%frazil_heat, ctx%work%sst, ctx%work%Tf, ctx%work%frzmlt, &
                            ctx%work%bank_ocean, ctx%work%bank)
   end subroutine hook_ocean_import

   subroutine hook_velocity_to_T(ctx, status)
      !! Runs before `ice_save_init`, so `aice == aice_init` here (the start-of-step concentration).
      type(ice_step_ctx_t), intent(inout) :: ctx
      integer, intent(out) :: status
      integer :: is, ie, js, je

      status = 0
      is = ctx%cice%nghost + 1
      ie = ctx%cice%nghost + ctx%cice%nx_phys
      js = ctx%cice%nghost + 1
      je = ctx%cice%nghost + ctx%cice%ny_phys
      associate (cm => S%ice_cgrid)
         call ice_atm_velocity_T(ctx%cice%nx, ctx%cice%ny, is, ie, js, je, ctx%cice%p%c%puny, &
                                 S%ice%u_ice, S%ice%v_ice, cm%earea, cm%narea, cm%uarea, cm%uvm, cm%hm, &
                                 ctx%cice%aice, ctx%work%uvel, ctx%work%vvel, ctx%cice%ice_err)
      end associate
   end subroutine hook_velocity_to_T

   subroutine hook_atm_prepare(ctx, status)
      type(ice_step_ctx_t), intent(inout) :: ctx
      integer, intent(out) :: status
      integer :: nx, ny, is, ie, js, je
      real(wp) :: zl, zs
      logical :: ok

      status = 0
      nx = ctx%cice%nx
      ny = ctx%cice%ny
      is = ctx%cice%nghost + 1
      ie = ctx%cice%nghost + ctx%cice%nx_phys
      js = ctx%cice%nghost + 1
      je = ctx%cice%nghost + ctx%cice%ny_phys
      call ice_atm_heights(S%icec%z_wind, zl, zs, ok)
      if (.not. ok) then
         status = 1
         return
      end if
      ctx%work%zlvl = zl
      ctx%work%zlvs = zs
      ! The atmosphere is the prescribed state the bulk forcing used THIS step (same arrays, same
      ! time level: `atm_state_update` ran at the top of the outer step) when
      ! `&ocean_atm_forcing_nml enable`, else the uniform one of the run slot (setting both is refused
      ! at configure).  The two paths are the same builder with different arrays; `tair_offset` is 0
      ! because the state is in kelvin on both.
      if (S%atm%enable) then
         call ice_atm_prepare_cell_inputs(nx, ny, is, ie, js, je, 0.0_wp, S%ice_cgrid%hm, S%icec%cosA, &
                                          S%icec%sinA, S%atm%tas, S%atm%huss, S%atm%uas, &
                                          S%atm%vas, S%atm%psl, S%atm%rlds, S%atm%prra, &
                                          S%atm%prsn, &
                                          ctx%work%potT, ctx%work%rhoa, ctx%work%wind, ctx%work%uatm, &
                                          ctx%work%vatm, ctx%work%Qa, ctx%work%flw, ctx%work%fsnow, &
                                          ctx%work%frain, ctx%cice%ice_err)
         call copy2_impl(nx, ny, S%atm%rsds, ctx%work%sw_down)
      else
         call ice_atm_prepare_cell_inputs(nx, ny, is, ie, js, je, 0.0_wp, S%ice_cgrid%hm, S%icec%cosA, &
                                          S%icec%sinA, S%icec%u_tas, S%icec%u_huss, S%icec%u_uas, &
                                          S%icec%u_vas, S%icec%u_psl, S%icec%u_rlds, S%icec%u_prra, &
                                          S%icec%u_prsn, &
                                          ctx%work%potT, ctx%work%rhoa, ctx%work%wind, ctx%work%uatm, &
                                          ctx%work%vatm, ctx%work%Qa, ctx%work%flw, ctx%work%fsnow, &
                                          ctx%work%frain, ctx%cice%ice_err)
         call copy2_impl(nx, ny, S%icec%u_rsds, ctx%work%sw_down)
      end if
   end subroutine hook_atm_prepare

   subroutine hook_halo_state(ctx, status)
      !! Used for both `halo_state` and `halo_state_post`.
      type(ice_step_ctx_t), intent(inout) :: ctx
      integer, intent(out) :: status

      status = 0
      if (iand(ice_mr_test_skip, ICE_SKIP_STATE_HALO) /= 0) return
      call wrap_category_state(ctx%cice)
   end subroutine hook_halo_state

   subroutine hook_stress_to_faces(ctx, status)
      type(ice_step_ctx_t), intent(inout) :: ctx
      integer, intent(out) :: status
      integer :: is, ie, js, je

      status = 0
      is = ctx%cice%nghost + 1
      ie = ctx%cice%nghost + ctx%cice%nx_phys
      js = ctx%cice%nghost + 1
      je = ctx%cice%nghost + ctx%cice%ny_phys
      ! the face hand-off averages the T-cell stress across the seam: refresh its ghost ring (exchange on a split
      ! axis, the local wrap otherwise, the vector fold on the fold row)
      if (iand(ice_mr_test_skip, ICE_SKIP_STRESS_HALO) == 0) then
         call ice_halo_centre_pair(G, S%bc, ctx%cice%nx, ctx%cice%ny, ctx%work%strair_x_T, ctx%work%strair_y_T, &
                                   test_fold=merge(2, merge(1, 0, iand(ice_mr_test_skip, ICE_SKIP_FOLD_SIGN) /= 0), &
                                                   iand(ice_mr_test_skip, ICE_SKIP_FOLD_PAIR) /= 0))
      end if
      call ice_atm_stress_to_faces(ctx%cice%nx, ctx%cice%ny, is, ie, js, je, ctx%work%strair_x_T, &
                                   ctx%work%strair_y_T, S%ice_cgrid%tarea, S%ice_cgrid%earea, &
                                   S%ice_cgrid%narea, S%ice%dync%strairxE, S%ice%dync%strairyN)
   end subroutine hook_stress_to_faces

   subroutine hook_dynamics(ctx, status)
      type(ice_step_ctx_t), intent(inout) :: ctx
      integer, intent(out) :: status
      type(ice_dyn_c_cat_t) :: cat
      integer :: nx, ny

      status = 0
      nx = ctx%cice%nx
      ny = ctx%cice%ny
      call ice_aggregates_fill(S%ice)
      cat%ncat = ctx%cice%ncat
      cat%wind_on_faces = .true.
      cat%aice_init => ctx%cice%aice_init
      cat%aicen => ctx%cice%aicen
      cat%vicen => ctx%cice%vicen
      call ice_dyn_c_step(G, S%metrics, S%coriolis_adv%f_corner, S%barotropic%b, S%ice_cgrid, S%ice, &
                          S%multilayer, S%bc, ctx%dt, cat=cat)
      call copy2_impl(nx, ny, S%ice%dync%strocnT_x, ctx%work%tau_x)
      call copy2_impl(nx, ny, S%ice%dync%strocnT_y, ctx%work%tau_y)
      call copy2_impl(nx, ny, S%ice%dync%w%rdg_conv, ctx%work%rdg_conv)
      call copy2_impl(nx, ny, S%ice%dync%w%rdg_shear, ctx%work%rdg_shear)
   end subroutine hook_dynamics

   subroutine hook_transport(ctx, status)
      type(ice_step_ctx_t), intent(inout) :: ctx
      integer, intent(out) :: status
      type(transport_result_t) :: tres
      integer :: is, ie, js, je

      status = 0
      is = ctx%cice%nghost + 1
      ie = ctx%cice%nghost + ctx%cice%nx_phys
      js = ctx%cice%nghost + 1
      je = ctx%cice%nghost + ctx%cice%ny_phys
      ! the decomposition seam of the transport (collectives, fold pairing), bound ONCE (the first step: it takes a
      ! min allreduce of the category batch size, which never changes); one rank and no fold is the default
      if (S%icec%tws%seam%nbatch == 0) then
         call ice_transport_seam_bind(S%icec%tws%seam, G, S%bc, S%icec%tws%fw%nbatch, S%icec%tws%ntrace)
      end if
      ! test harness only (`ice_mr_test_skip`): the fold-row flux pairing off, a negative control of the decomposition test
      S%icec%tws%seam%pair_fold = iand(ice_mr_test_skip, ICE_SKIP_FOLD_FLUX) == 0
      associate (cm => S%ice_cgrid, c => ctx%cice, w => ctx%work, f => ctx%flux)
         call ice_step_transport(S%icec%tws, c%p, ctx%dt, ctx%transport_post, is, ie, js, je, &
                                 S%ice%u_ice, S%ice%v_ice, cm%HTE, cm%HTN, cm%tarea, cm%dxT, cm%dyT, &
                                 cm%epm, cm%npm, cm%hm, w%Tf, c%hin_max, c%aice0, c%aicen, c%vicen, &
                                 c%vsnon, c%Tsfcn, c%qicen, c%sicen, c%qsnon, c%aice, c%vice, c%vsno, &
                                 w%Tsfc_c, w%qice_c, w%sice_c, w%qsno_c, f%fresh, f%fsalt, f%fhocn, &
                                 c%ice_err, tres, dxU=cm%dxU, dyU=cm%dyU, earea=cm%earea, narea=cm%narea, uvm=cm%uvm)
      end associate
      ! the step's clip/drop ledger of the remap update, kept for the conservation ledger (module slot)
      S%icec%tres = tres
      S%icec%ledger_cum(ICE_LCUM_NCLIP) = S%icec%ledger_cum(ICE_LCUM_NCLIP) + real(tres%ledger%n_clip, wp)
      S%icec%ledger_cum(ICE_LCUM_AREA) = S%icec%ledger_cum(ICE_LCUM_AREA) + tres%ledger%clip_area
      S%icec%ledger_cum(ICE_LCUM_VICE) = S%icec%ledger_cum(ICE_LCUM_VICE) + tres%ledger%drop_vice
      S%icec%ledger_cum(ICE_LCUM_VSNO) = S%icec%ledger_cum(ICE_LCUM_VSNO) + tres%ledger%drop_vsno
      S%icec%ledger_cum(ICE_LCUM_EICE) = S%icec%ledger_cum(ICE_LCUM_EICE) + tres%ledger%drop_eice
      S%icec%ledger_cum(ICE_LCUM_ESNO) = S%icec%ledger_cum(ICE_LCUM_ESNO) + tres%ledger%drop_esno
      ! A refused step (cell error bits) is reported by the sequencer's post-transport enforcement
      ! with the first offender; only a configuration refusal is a hook failure.
      if (tres%cfg /= TRANSPORT_CFG_OK) status = tres%cfg
   end subroutine hook_transport

   subroutine hook_ridging(ctx, status)
      type(ice_step_ctx_t), intent(inout) :: ctx
      integer, intent(out) :: status
      logical :: ok
      integer :: is, ie, js, je

      status = 0
      is = ctx%cice%nghost + 1
      ie = ctx%cice%nghost + ctx%cice%nx_phys
      js = ctx%cice%nghost + 1
      je = ctx%cice%nghost + ctx%cice%ny_phys
      associate (c => ctx%cice, w => ctx%work, f => ctx%flux)
         call ice_step_ridge(c%nx, c%ny, c%ncat, c%nilyr, c%nslyr, is, ie, js, je, c%p, c%hin_max, ctx%dt, 1, &
                             w%wet, w%Tf, w%rdg_conv, w%rdg_shear, c%aicen, c%vicen, c%vsnon, c%Tsfcn, &
                             c%qicen, c%sicen, c%qsnon, c%aice, c%vice, c%vsno, c%aice0, w%Tsfc_c, &
                             w%qice_c, w%sice_c, w%qsno_c, f%fresh, f%fsalt, f%fhocn, f%dardg1dt, &
                             f%dardg2dt, f%dvirdgdt, f%opening, w%niter, c%ice_err, ok)
         if (.not. ok) then
            status = 1
            return
         end if
         call ridge_bits_impl(c%nx, c%ny, is, ie, js, je, w%niter, c%ice_err)
      end associate
   end subroutine hook_ridging

   subroutine hook_ocean_handoff(ctx, status)
      type(ice_step_ctx_t), intent(inout) :: ctx
      integer, intent(out) :: status
      integer :: ierr

      status = 0
      call ice_flux_avg_accumulate(S%icec%avg, ctx%flux, ctx%dt, ierr)
      ! the exchanges below (stress seam, window apply) are collective: decide on the reduced status.  (Its
      ! failures are shape / not-initialised ones, rank-uniform by construction; the reduction makes it a fact.)
      call seam_max_one(ctx%seam, ierr)
      if (ierr /= ICE_AVG_OK) then
         status = ierr
         return
      end if
      ! Rain that falls through the ice cover is in `fresh` (the ice step books it so its own water
      ! ledger closes) but the ocean's own precipitation on the open fraction carries ENTHALPY ONLY
      ! (the bulk forcing writes no salt or mass term), so the virtual-salt dilution of the hand-over
      ! must not count it either.  Prescribed atmosphere only: with the uniform atmosphere the rain is
      ! zero and nothing runs.  Under `&ocean_forcing_nml freshwater = "virtual"` the open fraction's
      ! rain DOES dilute (`salt_fw`), so the ice-covered share must too: the term stays in `fresh`.
      if (S%atm%enable .and. .not. S%surface_flux%fw_virtual) then
         call rain_through_impl(ctx%work%nx, ctx%work%ny, ctx%cice%nghost + 1, &
                                ctx%cice%nghost + ctx%cice%nx_phys, ctx%cice%nghost + 1, &
                                ctx%cice%nghost + ctx%cice%ny_phys, ctx%dt/S%icec%avg%t_acc, &
                                ctx%work%frain, ctx%cice%aice_init, S%icec%avg%acc_fresh)
      end if
      call ice_ocean_stress_flux_c(S%metrics, S%surface_stress, S%ice, G, S%bc)
      if (ctx%thermo_step .and. S%dyn%enable_thermodynamics) then
         call ice_flux_avg_stash_pending(S%icec%avg)
         ! `apply_window` returns before its exchanges only on a shape error, an empty window or a window length
         ! different from the thermo window: all functions of the configuration and of the (rank-uniform) step
         ! schedule, never of a rank's data, so every rank takes the same branch.
         call ice_cice_apply_window(G, S%bc, S%icec%avg, S%surface_flux, S%icec%oin%sss_ocn, ierr, &
                                    window_expected=S%dyn%therm_dt(ctx%dt))
         if (ierr /= ICE_ADAPT_OK) then
            status = 100 + ierr
            return
         end if
      end if
      call copy2_impl(ctx%work%nx, ctx%work%ny, ctx%work%bank, ctx%work%bank_ocean)
      call copy2_impl(ctx%work%nx, ctx%work%ny, ctx%work%bank_ocean, S%ice%frazil_heat)
      call ice_aggregates_fill(S%ice)
   end subroutine hook_ocean_handoff

end module rdb_ice_cice_hooks
