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
   !!  9. and 15. `halo_state`, `halo_state_post`: single rank, so a LOCAL refresh: the periodic
   !!     wrap of the category state and `aice0` (a wall edge keeps its ghosts: the C-grid metrics pin
   !!     them to land).  No tripolar fold, no exchange: refused at configure.
   !! 11. `stress_to_faces`: periodic wrap of the T wind stress, then `ice_atm_stress_to_faces`
   !!     into the dynamics slot's `strairxE/yN` (the dynamics keeps them: `wind_on_faces`).
   !! 12. `dynamics`: the scheme-neutral aggregates are refreshed from the haloed state, then the
   !!     C-grid EVP of `ice_dyn_c_step` runs with the CATEGORY strength (`kstrength`) from the
   !!     post-thermodynamics state; `work%tau_x/y` (the persistent centre ice-ocean stress) and
   !!     `work%rdg_conv/rdg_shear` are copied from the dynamics slot.
   !! 13. `transport`: `ice_step_transport` with the post mode of the context (NONE).
   !! 14. `ridging`: `ice_step_ridge`; cells that needed more than one pass get the counted
   !!     `ICE_ERR_RIDGE_PASSES` bit.
   !! 17. `ocean_handoff`: accumulate the step's fluxes into the window mean; the face stress
   !!     blend (`ice_ocean_stress_flux_c`, its own restart mirror); at the ocean thermo step the
   !!     adapters apply the window mean (`ice_cice_apply_window`, virtual-salt form, with the seam
   !!     halo exchange) and reset the window; the spent bank is written back to the ocean-side
   !!     bank; the scheme-neutral aggregates are refreshed for every consumer (diagnostics,
   !!     console totals, the open-water weight of the bulk forcing).
   !!
   !! **Interim choices.**  The stress blend uses the configure-time wind-stress snapshot of the
   !! ocean (`ice%tau_a_x/y`, as the SIS2 coupler does): correct because the only wind this scheme
   !! accepts is the constant one of `&physics_nml` (time-varying forcing is refused at configure),
   !! while the ice itself feels the wind of the atmosphere boundary layer; the `ice_cice_stress_blend` of the adapters (a different
   !! formulation: per-ice-area stress) is not used, the face blend of the C-grid dynamics
   !! conserves the momentum of the pair exactly.
   use, intrinsic :: iso_fortran_env, only: int32
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_config, only: config_t
   use rdb_ocean_state, only: ocean_state_t
   use rdb_ocean_periodic, only: ocean_periodic_wrap_centre_2d, ocean_periodic_wrap_centre_3d
   use rdb_ice_state, only: ice_aggregates_fill
   use rdb_ice_cice_state, only: ice_cice_state_t, ice_update_state, ice_state_seed_uniform
   use rdb_ice_step_cice, only: ice_step_cice, ice_step_hooks_t, ice_step_opts_t, ice_step_result_t, &
                                ice_step_ctx_t, ice_stage_name, ICE_STEP_OK, ICE_STEP_FATAL, &
                                ICE_STAGE_OCEAN_IMPORT, ICE_STAGE_VELOCITY_TO_T, ICE_STAGE_ATM_PREPARE, &
                                ICE_STAGE_HALO_STATE, ICE_STAGE_STRESS_TO_FACES, ICE_STAGE_DYNAMICS, &
                                ICE_STAGE_TRANSPORT, ICE_STAGE_RIDGING, ICE_STAGE_HALO_STATE_POST, &
                                ICE_STAGE_OCEAN_HANDOFF
   use rdb_ice_errors, only: ice_err_enforce, ICE_ERR_RIDGE_PASSES, ICE_ERR_FATAL_MASK
   use rdb_ice_melt_potential, only: ice_melt_potential
   use rdb_ice_import_ocean, only: ice_import_ocean, ICE_IMPORT_OK
   use rdb_ice_atm_inputs, only: ice_atm_prepare_cell_inputs, ice_atm_heights, ice_atm_velocity_T, &
                                 ice_atm_stress_to_faces
   use rdb_ice_dyn_c, only: ice_dyn_c_step, ice_dyn_c_cat_t
   use rdb_ice_transport_cice, only: ice_step_transport, transport_result_t, TRANSPORT_CFG_OK
   use rdb_ice_ridge, only: ice_step_ridge
   use rdb_ice_flux_avg, only: ice_flux_avg_accumulate, ice_flux_avg_stash_pending, ICE_AVG_OK
   use rdb_ice_cice_adapters, only: ice_cice_apply_window, ICE_ADAPT_OK
   use rdb_ice_ocean_coupler, only: ice_ocean_stress_flux_c
   use rdb_ice_init, only: ice_ic_params_t, ICE_IC_CONC_UNIFORM, ICE_IC_CONC_LATITUDES
   use pic_logger, only: logger => global_logger
   use pic_strings, only: to_string
   implicit none
   private

   public :: ice_cice_bind, ice_cice_step, ice_cice_ic_apply, ice_cice_post_enter
   public :: ice_cice_hooks_manifest

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

   subroutine ice_cice_hooks_manifest(installed)
      !! Which hook slots hold an implementation (for the tests): `installed(stage)`.
      logical, intent(out) :: installed(:)
      integer :: k

      do k = 1, size(installed)
         installed(k) = HOOKS%installed(k)
      end do
   end subroutine ice_cice_hooks_manifest

   subroutine ice_cice_step(dt, thermo_step, istep)
      !! One ice step.  A refusal or a stage failure stops the run with the stage named; a fatal
      !! cell error bit aborts through `ice_err_enforce` (first offender dumped).
      real(wp), intent(in) :: dt
      logical, intent(in) :: thermo_step
      integer, intent(in) :: istep
      type(ice_step_opts_t) :: opts
      type(ice_step_result_t) :: res

      opts%ndtd = 1
      opts%nranks = 1
      opts%dt_outer = dt
      opts%fatal_mask = ICE_ERR_FATAL_MASK
      call ice_step_cice(S%ice%cice, S%icec%flux, S%icec%rad, S%icec%work, HOOKS, opts, dt, thermo_step, istep, res)
      S%icec%n_steps = S%icec%n_steps + 1
      if (res%status == ICE_STEP_OK) return
      if (res%status == ICE_STEP_FATAL) call ice_err_enforce(res%summary, G%nx_phys, opts%fatal_mask)
      call logger%error("ice_step_cice: step "//to_string(istep)//" stopped, status "//to_string(res%status)// &
                        ", stage '"//trim(ice_stage_name(res%stage))//"' (raw status "// &
                        to_string(res%stage_status)//"): "//trim(res%message))
      error stop "ice_step_cice: the ice step failed (see the log)"
   end subroutine ice_cice_step

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

   subroutine wrap4_impl(nx, ny, nc, nl, nxp, nyp, ng, wrap_x, wrap_y, q)
      !! Periodic ghost wrap of a `(nx, ny, ncat, nlayer)` field (x then y, as the ocean's).
      integer, intent(in) :: nx, ny, nc, nl, nxp, nyp, ng
      logical, intent(in) :: wrap_x, wrap_y
      real(wp), intent(inout) :: q(nx, ny, nc, nl)
      integer :: i, j, c, l, g

      if (wrap_x) then
         do concurrent(l=1:nl, c=1:nc, j=1:ny, g=1:ng)
            q(g, j, c, l) = q(g + nxp, j, c, l)
            q(nxp + ng + g, j, c, l) = q(ng + g, j, c, l)
         end do
      end if
      if (wrap_y) then
         do concurrent(l=1:nl, c=1:nc, g=1:ng, i=1:nx)
            q(i, g, c, l) = q(i, g + nyp, c, l)
            q(i, nyp + ng + g, c, l) = q(i, ng + g, c, l)
         end do
      end if
   end subroutine wrap4_impl

   subroutine wrap_category_state(c)
      !! The local halo refresh of the category state and `aice0` (module docstring).
      type(ice_cice_state_t), intent(inout) :: c
      logical :: wx, wy

      wx = S%bc%periodic_x
      wy = S%bc%periodic_y
      if (.not. (wx .or. wy)) return
      call ocean_periodic_wrap_centre_3d(c%aicen, c%nx, c%ny, c%ncat, c%nx_phys, c%ny_phys, c%nghost, wx, wy)
      call ocean_periodic_wrap_centre_3d(c%vicen, c%nx, c%ny, c%ncat, c%nx_phys, c%ny_phys, c%nghost, wx, wy)
      call ocean_periodic_wrap_centre_3d(c%vsnon, c%nx, c%ny, c%ncat, c%nx_phys, c%ny_phys, c%nghost, wx, wy)
      call ocean_periodic_wrap_centre_3d(c%Tsfcn, c%nx, c%ny, c%ncat, c%nx_phys, c%ny_phys, c%nghost, wx, wy)
      call wrap4_impl(c%nx, c%ny, c%ncat, c%nilyr, c%nx_phys, c%ny_phys, c%nghost, wx, wy, c%qicen)
      call wrap4_impl(c%nx, c%ny, c%ncat, c%nilyr, c%nx_phys, c%ny_phys, c%nghost, wx, wy, c%sicen)
      call wrap4_impl(c%nx, c%ny, c%ncat, c%nslyr, c%nx_phys, c%ny_phys, c%nghost, wx, wy, c%qsnon)
      call ocean_periodic_wrap_centre_2d(c%aice0, c%nx, c%ny, c%nx_phys, c%ny_phys, c%nghost, wx, wy)
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
      ! The atmosphere is the uniform one of the run slot: a prescribed atmospheric state / bulk
      ! forcing is refused at configure under this scheme (the ocean-side stress blend needs the
      ! unweighted wind stress of the step).
      call ice_atm_prepare_cell_inputs(nx, ny, is, ie, js, je, 0.0_wp, S%ice_cgrid%hm, S%icec%cosA, &
                                       S%icec%sinA, S%icec%u_tas, S%icec%u_huss, S%icec%u_uas, &
                                       S%icec%u_vas, S%icec%u_psl, S%icec%u_rlds, S%icec%u_prra, &
                                       S%icec%u_prsn, &
                                       ctx%work%potT, ctx%work%rhoa, ctx%work%wind, ctx%work%uatm, &
                                       ctx%work%vatm, ctx%work%Qa, ctx%work%flw, ctx%work%fsnow, &
                                       ctx%work%frain, ctx%cice%ice_err)
      call copy2_impl(nx, ny, S%icec%u_rsds, ctx%work%sw_down)
   end subroutine hook_atm_prepare

   subroutine hook_halo_state(ctx, status)
      !! Used for both `halo_state` and `halo_state_post`.
      type(ice_step_ctx_t), intent(inout) :: ctx
      integer, intent(out) :: status

      status = 0
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
      call ocean_periodic_wrap_centre_2d(ctx%work%strair_x_T, ctx%cice%nx, ctx%cice%ny, ctx%cice%nx_phys, &
                                         ctx%cice%ny_phys, ctx%cice%nghost, S%bc%periodic_x, S%bc%periodic_y)
      call ocean_periodic_wrap_centre_2d(ctx%work%strair_y_T, ctx%cice%nx, ctx%cice%ny, ctx%cice%nx_phys, &
                                         ctx%cice%ny_phys, ctx%cice%nghost, S%bc%periodic_x, S%bc%periodic_y)
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
      associate (cm => S%ice_cgrid, c => ctx%cice, w => ctx%work, f => ctx%flux)
         call ice_step_transport(S%icec%tws, c%p, ctx%dt, ctx%transport_post, is, ie, js, je, &
                                 S%ice%u_ice, S%ice%v_ice, cm%HTE, cm%HTN, cm%tarea, cm%dxT, cm%dyT, &
                                 cm%epm, cm%npm, cm%hm, w%Tf, c%hin_max, c%aice0, c%aicen, c%vicen, &
                                 c%vsnon, c%Tsfcn, c%qicen, c%sicen, c%qsnon, c%aice, c%vice, c%vsno, &
                                 w%Tsfc_c, w%qice_c, w%sice_c, w%qsno_c, f%fresh, f%fsalt, f%fhocn, &
                                 c%ice_err, tres)
      end associate
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
      if (ierr /= ICE_AVG_OK) then
         status = ierr
         return
      end if
      call ice_ocean_stress_flux_c(S%metrics, S%surface_stress, S%ice, G, S%bc)
      if (ctx%thermo_step .and. S%dyn%enable_thermodynamics) then
         call ice_flux_avg_stash_pending(S%icec%avg)
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
