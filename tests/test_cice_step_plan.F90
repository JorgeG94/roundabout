!! The ice step sequencer (`ice_step_cice`, D11.3): stage presence, ordering, contracts.
!!
!! The manifest is asserted against a literal list written from the design (not read back from
!! the module), with each ordering fact as a predicate; a missing required stage refuses with its
!! name and leaves a seeded state bitwise untouched; recording doubles prove the call order and
!! the cadence; a fatal bit stops the step before the next group while a counted one does not;
!! the frazil bank is frozen once per window; the ice-ocean stress survives a restart; and two
!! 50-step smoke runs (a dark one and a daylight one with shortwave, a bank, snow, rain and a melt
!! potential) with the REAL thermodynamic stages and no-motion doubles for the dynamics, transport
!! and ridging close the shipped conservation ledger.  The doubles act through `ctx` (the hook
!! contract), never through a module global.  Links the main library (the slots need the grid).
module test_cice_step_plan
   use, intrinsic :: iso_fortran_env, only: int32, int64, output_unit
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_decomp, only: decomp_t, decomp_init
   use rdb_ocean_restart, only: restart_registry_t
   use rdb_ocean_restart_io, only: ocean_restart_write_local, ocean_restart_read_local, ocean_restart_metadata_t
   use rdb_ice_params, only: ice_params_t, ice_params_default, ICE_CONSTANTS_ICEPACK, ice_q_from_T, ice_qsno_from_T, &
                             ice_salinity_profile, ICE_ADVECT_UPWIND
   use rdb_ice_cice_state, only: ice_cice_state_t, ice_update_state
   use rdb_ice_flux_state, only: ice_flux_t
   use rdb_ice_rad_state, only: ice_rad_state_t
   use rdb_ice_sw_wiring, only: ice_sw_wiring_cleanup
   use rdb_ice_ledger, only: ice_ledger_mirror_t, ICE_LEDGER_WATER, ICE_LEDGER_SALT, ICE_LEDGER_HEAT
   use rdb_ice_thickness, only: ICE_FERRMAX
   use rdb_ice_errors, only: ICE_ERR_DYN_VEL, ICE_ERR_RIDGE_PASSES, ICE_ERR_NONFINITE, ICE_ERR_TRANSPORT_CFL
   use rdb_ice_transport_cice, only: ICE_TRANSPORT_POST_NONE, ICE_TRANSPORT_POST_AGGREGATE
   use rdb_ice_step_cice, only: ice_step_cice, ice_step_plan_build, ice_step_plan_t, ice_step_hooks_t, ice_step_opts_t, &
                                ice_step_result_t, ice_step_ctx_t, ice_step_work_t, ice_stage_name, &
                                ICE_NSTAGE, ICE_STAGE_OCEAN_IMPORT, ICE_STAGE_FLUX_ZERO, ICE_STAGE_SAVE_INIT, &
                                ICE_STAGE_VELOCITY_TO_T, ICE_STAGE_ATM_PREPARE, ICE_STAGE_RADIATION, ICE_STAGE_THERM1, &
                                ICE_STAGE_THERM2, ICE_STAGE_HALO_STATE, ICE_STAGE_AGGREGATE_GHOST, &
                                ICE_STAGE_STRESS_TO_FACES, ICE_STAGE_DYNAMICS, ICE_STAGE_TRANSPORT, &
                                ICE_STAGE_RIDGING, ICE_STAGE_HALO_STATE_POST, ICE_STAGE_AGGREGATE, ICE_STAGE_OCEAN_HANDOFF, &
                                ICE_GROUP_PREP, ICE_GROUP_THERMO, ICE_GROUP_DYNAMICS, ICE_GROUP_TRANSPORT, ICE_GROUP_CLOSE, &
                                ICE_GROUP_OUT, ICE_STEP_OK, ICE_STEP_REFUSED_MISSING, ICE_STEP_REFUSED_DOUBLE, &
                                ICE_STEP_REFUSED_ENVELOPE, ICE_STEP_STAGE_FAILED, ICE_STEP_FATAL, ICE_HOOK_MISSING
   implicit none
   private

   public :: collect_cice_step_plan_tests

   integer, parameter :: NG = 3, NXP = 4, NYP = 3, NCAT = 5, NIL = 4, IC = 4, JC = 4, NLOGMAX = 512
   real(wp), parameter :: AREA = 1.0e6_wp, DT = 3600.0_wp, TFRZ = -1.8_wp

   type :: world_t
      type(hgrid_t) :: g
      type(ice_params_t) :: p
      type(ice_cice_state_t) :: st
      type(ice_flux_t) :: fx
      type(ice_rad_state_t) :: rad
      type(ice_step_work_t) :: wk
      real(wp), allocatable :: area(:, :)
      logical :: mapped = .false.
   end type world_t

   ! the recorder (an observer: it only logs what the hooks were handed)
   integer, save :: nlog = 0
   integer, save :: log_stage(NLOGMAX) = 0
   logical, save :: log_thermo(NLOGMAX) = .false.
   integer, save :: log_post(NLOGMAX) = -1
   integer, save :: log_istep(NLOGMAX) = -1
   real(wp), save :: log_dt(NLOGMAX) = 0.0_wp
   ! the scenario and the fault injection the doubles apply (inputs, not world state)
   integer, save :: inject_stage = 0
   integer(int32), save :: inject_bit = 0_int32
   integer, save :: fail_stage = 0
   logical, save :: nan_sw = .false.
   real(wp), save :: scn_sw = 0.0_wp
   logical, save :: scn_melt = .false.
   logical, save :: scn_precip = .false.
   logical, save :: scn_ghost = .false.

contains

   subroutine collect_cice_step_plan_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("plan_is_the_literal_sequence", test_plan_literal), &
                  new_unittest("plan_ordering_predicates", test_plan_predicates), &
                  new_unittest("missing_required_stage_refuses_untouched", test_missing_refuses), &
                  new_unittest("test_doubles_refused_outside_tests", test_doubles_refused), &
                  new_unittest("envelope_refuses_untouched", test_envelope), &
                  new_unittest("call_order_and_cadence", test_order_cadence), &
                  new_unittest("highfreq_adds_velocity_stage", test_highfreq), &
                  new_unittest("fatal_bit_stops_before_next_group", test_fatal), &
                  new_unittest("transport_fatal_stops_before_ridging", test_transport_fatal), &
                  new_unittest("counted_bit_does_not_stop", test_counted), &
                  new_unittest("stage_failure_stops_the_step", test_stage_failure), &
                  new_unittest("halo_then_ghost_aggregate_before_dynamics", test_halo_ghost), &
                  new_unittest("bank_is_spent_once_over_a_window", test_bank_once), &
                  new_unittest("stress_state_survives_a_restart", test_stress_restart), &
                  new_unittest("smoke_dark_50_steps_ledger_closes", test_smoke_dark), &
                  new_unittest("smoke_daylight_50_steps_ledger_closes", test_smoke_daylight) &
                  ]
   end subroutine collect_cice_step_plan_tests

   ! ------------------------------------------------------------------
   ! device helpers (one-level components, so every directive names a plain dummy's component)
   ! ------------------------------------------------------------------

   subroutine push_err(st)
      type(ice_cice_state_t), intent(inout) :: st
      !$acc update device(st%ice_err)
   end subroutine push_err

   subroutine push_ocn(wk)
      type(ice_step_work_t), intent(inout) :: wk
      !$acc update device(wk%Tf, wk%sst, wk%frzmlt, wk%bank)
   end subroutine push_ocn

   subroutine push_atm(wk)
      type(ice_step_work_t), intent(inout) :: wk
      !$acc update device(wk%uatm, wk%vatm, wk%wind, wk%potT, wk%rhoa, wk%Qa, wk%flw)
      !$acc update device(wk%fsnow, wk%frain, wk%sw_down)
   end subroutine push_atm

   subroutine push_tau(wk)
      type(ice_step_work_t), intent(inout) :: wk
      !$acc update device(wk%tau_x, wk%tau_y)
   end subroutine push_tau

   subroutine push_bank_ocean(wk)
      type(ice_step_work_t), intent(inout) :: wk
      !$acc update device(wk%bank_ocean)
   end subroutine push_bank_ocean

   subroutine pull_bank(wk)
      type(ice_step_work_t), intent(inout) :: wk
      !$acc update self(wk%bank, wk%bank_ocean, wk%frzmlt)
   end subroutine pull_bank

   subroutine pull_work(wk)
      type(ice_step_work_t), intent(inout) :: wk
      !$acc update self(wk%tau_x, wk%tau_y, wk%flw, wk%bank, wk%bank_ocean, wk%frzmlt)
   end subroutine pull_work

   subroutine push_state_ring(st)
      type(ice_cice_state_t), intent(inout) :: st
      !$acc update device(st%aicen, st%aice0)
   end subroutine push_state_ring

   subroutine pull_state(st, fx)
      type(ice_cice_state_t), intent(inout) :: st
      type(ice_flux_t), intent(inout) :: fx
      !$acc update self(st%aicen, st%vicen, st%vsnon, st%Tsfcn, st%qicen, st%qsnon, st%aice, st%vice)
      !$acc update self(st%vsno, st%aice0, st%aicen_init, st%vicen_init, st%aice_init, st%ice_err, st%sicen)
      !$acc update self(fx%fresh, fx%fsalt, fx%fhocn, fx%fswthru, fx%flwout, fx%freshn, fx%fhocnn)
      !$acc update self(fx%frazil_vol, fx%meltl)
   end subroutine pull_state

   subroutine push_restart(st, wk)
      !! After a restart read (host I/O): every registered array to the device.
      type(ice_cice_state_t), intent(inout) :: st
      type(ice_step_work_t), intent(inout) :: wk
      !$acc update device(st%aicen, st%vicen, st%vsnon, st%Tsfcn, st%qicen, st%sicen, st%qsnon)
      !$acc update device(wk%tau_x, wk%tau_y)
   end subroutine push_restart

   subroutine map_area(a)
      real(wp), intent(inout) :: a(:, :)
      !$acc enter data copyin(a)
   end subroutine map_area

   subroutine unmap_area(a)
      real(wp), intent(inout) :: a(:, :)
      !$acc exit data delete(a)
   end subroutine unmap_area

   subroutine pull_ledger_inputs(st, wk)
      type(ice_cice_state_t), intent(inout) :: st
      type(ice_step_work_t), intent(inout) :: wk
      !$acc update self(st%aice_init, st%aice)
      !$acc update self(wk%flw, wk%fsnow, wk%frain, wk%bank, wk%bank_ocean, wk%frzmlt)
   end subroutine pull_ledger_inputs

   ! ------------------------------------------------------------------
   ! doubles (they act through ctx)
   ! ------------------------------------------------------------------

   subroutine record(ctx, stage, status)
      type(ice_step_ctx_t), intent(inout) :: ctx
      integer, intent(in) :: stage
      integer, intent(out) :: status

      if (nlog < NLOGMAX) then
         nlog = nlog + 1
         log_stage(nlog) = stage
         log_thermo(nlog) = ctx%thermo_step
         log_post(nlog) = ctx%transport_post
         log_istep(nlog) = ctx%istep
         log_dt(nlog) = ctx%dt
      end if
      status = 0
      if (stage == fail_stage) status = 7
      if (stage == inject_stage) then
         ctx%cice%ice_err = 0_int32
         ctx%cice%ice_err(IC, JC) = inject_bit
         call push_err(ctx%cice)
      end if
   end subroutine record

   subroutine dbl_ocean_import(ctx, status)
      !! The bank contract: `frzmlt` from the CURRENT ocean-side bank every step.
      type(ice_step_ctx_t), intent(inout) :: ctx
      integer, intent(out) :: status
      integer :: i, j, ng_

      call pull_bank(ctx%work)
      ng_ = ctx%cice%nghost
      ctx%work%Tf = TFRZ
      ctx%work%sst = TFRZ
      ctx%work%bank = ctx%work%bank_ocean
      ctx%work%frzmlt = 0.0_wp
      do j = 1, ctx%cice%ny
         do i = 1, ctx%cice%nx
            if (ctx%work%bank(i, j) > 0.0_wp) then
               ctx%work%frzmlt(i, j) = ctx%work%bank(i, j)/ctx%dt
            else if (scn_melt .and. mod(i - ng_, 2) == 0) then
               ctx%work%frzmlt(i, j) = -60.0_wp
               ctx%work%sst(i, j) = TFRZ + 0.3_wp
            end if
         end do
      end do
      call push_ocn(ctx%work)
      call record(ctx, ICE_STAGE_OCEAN_IMPORT, status)
   end subroutine dbl_ocean_import

   subroutine dbl_velocity(ctx, status)
      type(ice_step_ctx_t), intent(inout) :: ctx
      integer, intent(out) :: status

      call record(ctx, ICE_STAGE_VELOCITY_TO_T, status)
   end subroutine dbl_velocity

   subroutine dbl_atm(ctx, status)
      !! A cold atmosphere over every cell; shortwave, snow and rain per the scenario.
      type(ice_step_ctx_t), intent(inout) :: ctx
      integer, intent(out) :: status
      real(wp), parameter :: T = 268.0_wp, SIGMA = 5.67e-8_wp

      ctx%work%uatm = 2.0_wp
      ctx%work%vatm = 0.0_wp
      ctx%work%wind = 2.0_wp
      ctx%work%potT = T
      ctx%work%rhoa = 1.3_wp
      ctx%work%Qa = 0.8_wp*11637800.0_wp*exp(-5897.8_wp/T)/1.3_wp
      ctx%work%flw = 0.8_wp*SIGMA*T**4
      ctx%work%fsnow = 0.0_wp
      ctx%work%frain = 0.0_wp
      if (scn_precip .and. ctx%istep <= 10) ctx%work%fsnow = 2.0e-4_wp
      if (scn_precip .and. ctx%istep > 10 .and. ctx%istep <= 20) ctx%work%frain = 1.0e-4_wp
      ctx%work%sw_down = scn_sw
      if (nan_sw) ctx%work%sw_down(IC, JC) = ieee_value(1.0_wp, ieee_quiet_nan)
      ctx%work%zlvl = 10.0_wp
      ctx%work%zlvs = 10.0_wp
      call push_atm(ctx%work)
      call record(ctx, ICE_STAGE_ATM_PREPARE, status)
   end subroutine dbl_atm

   subroutine dbl_halo(ctx, status)
      !! A ghost ring that carries ice (what an exchange would deliver); the sequencer's ghost
      !! aggregate must then make the ghost aggregates consistent.
      type(ice_step_ctx_t), intent(inout) :: ctx
      integer, intent(out) :: status

      if (scn_ghost) then
         ctx%cice%aicen(1, 1, :) = 0.1_wp
         call push_state_ring(ctx%cice)
      end if
      call record(ctx, ICE_STAGE_HALO_STATE, status)
   end subroutine dbl_halo

   subroutine dbl_halo_post(ctx, status)
      type(ice_step_ctx_t), intent(inout) :: ctx
      integer, intent(out) :: status

      call record(ctx, ICE_STAGE_HALO_STATE_POST, status)
   end subroutine dbl_halo_post

   subroutine dbl_stress(ctx, status)
      type(ice_step_ctx_t), intent(inout) :: ctx
      integer, intent(out) :: status

      call record(ctx, ICE_STAGE_STRESS_TO_FACES, status)
   end subroutine dbl_stress

   subroutine dbl_dynamics(ctx, status)
      !! No motion; the persistent stress of this step (a known function of the step counter).
      type(ice_step_ctx_t), intent(inout) :: ctx
      integer, intent(out) :: status

      ctx%work%tau_x = 0.005_wp*real(ctx%istep, wp)
      ctx%work%tau_y = 0.002_wp*real(ctx%istep, wp)
      call push_tau(ctx%work)
      call record(ctx, ICE_STAGE_DYNAMICS, status)
   end subroutine dbl_dynamics

   subroutine dbl_transport(ctx, status)
      type(ice_step_ctx_t), intent(inout) :: ctx
      integer, intent(out) :: status

      call record(ctx, ICE_STAGE_TRANSPORT, status)
   end subroutine dbl_transport

   subroutine dbl_ridging(ctx, status)
      type(ice_step_ctx_t), intent(inout) :: ctx
      integer, intent(out) :: status

      call record(ctx, ICE_STAGE_RIDGING, status)
   end subroutine dbl_ridging

   subroutine dbl_handoff(ctx, status)
      !! The single write-back owner of the bank: the spent bank goes back to the ocean side.
      type(ice_step_ctx_t), intent(inout) :: ctx
      integer, intent(out) :: status

      call pull_bank(ctx%work)
      ctx%work%bank_ocean = ctx%work%bank
      call push_bank_ocean(ctx%work)
      call record(ctx, ICE_STAGE_OCEAN_HANDOFF, status)
   end subroutine dbl_handoff

   subroutine all_doubles(h, with_velocity)
      type(ice_step_hooks_t), intent(inout) :: h
      logical, intent(in) :: with_velocity

      call h%install(ICE_STAGE_OCEAN_IMPORT, dbl_ocean_import, is_double=.true.)
      if (with_velocity) call h%install(ICE_STAGE_VELOCITY_TO_T, dbl_velocity, is_double=.true.)
      call h%install(ICE_STAGE_ATM_PREPARE, dbl_atm, is_double=.true.)
      call h%install(ICE_STAGE_HALO_STATE, dbl_halo, is_double=.true.)
      call h%install(ICE_STAGE_STRESS_TO_FACES, dbl_stress, is_double=.true.)
      call h%install(ICE_STAGE_DYNAMICS, dbl_dynamics, is_double=.true.)
      call h%install(ICE_STAGE_TRANSPORT, dbl_transport, is_double=.true.)
      call h%install(ICE_STAGE_RIDGING, dbl_ridging, is_double=.true.)
      call h%install(ICE_STAGE_HALO_STATE_POST, dbl_halo_post, is_double=.true.)
      call h%install(ICE_STAGE_OCEAN_HANDOFF, dbl_handoff, is_double=.true.)
   end subroutine all_doubles

   subroutine reset_log()
      nlog = 0
      log_stage = 0
      log_thermo = .false.
      log_post = -1
      log_istep = -1
      log_dt = 0.0_wp
      inject_stage = 0
      inject_bit = 0_int32
      fail_stage = 0
      nan_sw = .false.
      scn_sw = 0.0_wp
      scn_melt = .false.
      scn_precip = .false.
      scn_ghost = .false.
   end subroutine reset_log

   ! ------------------------------------------------------------------
   ! the world
   ! ------------------------------------------------------------------

   subroutine seed_cell(wd, i, j, scale)
      type(world_t), intent(inout) :: wd
      integer, intent(in) :: i, j
      real(wp), intent(in) :: scale
      real(wp), parameter :: A(NCAT) = [0.05_wp, 0.15_wp, 0.30_wp, 0.25_wp, 0.10_wp]
      real(wp), parameter :: HS(NCAT) = [0.0_wp, 0.03_wp, 0.08_wp, 0.15_wp, 0.2_wp]
      real(wp) :: Sp(NIL + 1), Tz, hlo, hhi
      integer :: n, k

      call ice_salinity_profile(NIL, wd%p%c%saltmax, wd%p%c%min_salin, Sp)
      do n = 1, NCAT
         hlo = wd%st%hin_max(n - 1)
         hhi = wd%st%hin_max(n)
         wd%st%aicen(i, j, n) = A(n)*scale
         wd%st%vicen(i, j, n) = A(n)*scale*(hlo + 0.5_wp*min(hhi - hlo, 1.0_wp))
         wd%st%vsnon(i, j, n) = A(n)*scale*HS(n)
         wd%st%Tsfcn(i, j, n) = -9.0_wp
         do k = 1, NIL
            wd%st%sicen(i, j, n, k) = Sp(k)
            Tz = -7.0_wp + 5.0_wp*real(k - 1, wp)/real(NIL - 1, wp)
            wd%st%qicen(i, j, n, k) = ice_q_from_T(Tz, Sp(k), wd%p%c, wd%p%thermo%l_brine)
         end do
         wd%st%qsnon(i, j, n, 1) = ice_qsno_from_T(-4.5_wp, wd%p%c)
      end do
   end subroutine seed_cell

   subroutine make_world(wd, ng, map, set_wet, error)
      !! A seeded world: cells differ (the scale varies), one owned cell is ice free.
      type(world_t), intent(inout) :: wd
      integer, intent(in) :: ng
      logical, intent(in) :: map, set_wet
      type(error_type), allocatable, intent(out) :: error
      integer :: ierr, i, j
      logical, allocatable :: wet(:, :)

      call wd%g%init(NXP, NYP, ng, 1000.0_wp, 1000.0_wp)
      wd%p = ice_params_default(ICE_CONSTANTS_ICEPACK, NCAT, NIL, 1)
      wd%p%transport%advection = ICE_ADVECT_UPWIND
      call wd%st%init(wd%g, wd%p, ierr)
      call check(error, ierr, 0, "state init")
      if (allocated(error)) return
      call wd%fx%init(wd%g, NCAT, ierr)
      call check(error, ierr, 0, "flux init")
      if (allocated(error)) return
      call wd%rad%init(wd%st%nx, wd%st%ny, NCAT, NIL, ierr)
      call check(error, ierr, 0, "radiation init")
      if (allocated(error)) return
      call wd%rad%configure(wd%p%sw%ahmax, wd%p%c%snowpatch)
      call wd%wk%init(wd%st, ierr)
      call check(error, ierr, 0, "work init")
      if (allocated(error)) return
      do j = ng + 1, ng + NYP
         do i = ng + 1, ng + NXP
            if (i == ng + 2 .and. j == ng + 2) cycle      ! open water
            call seed_cell(wd, i, j, 0.6_wp + 0.1_wp*real(i - ng, wp) + 0.05_wp*real(j - ng, wp))
         end do
      end do
      wd%st%aice = sum(wd%st%aicen, dim=3)
      wd%st%vice = sum(wd%st%vicen, dim=3)
      wd%st%vsno = sum(wd%st%vsnon, dim=3)
      wd%st%aice0 = max(1.0_wp - wd%st%aice, 0.0_wp)
      if (set_wet) then
         allocate (wet(wd%st%nx, wd%st%ny), source=.true.)
         call wd%wk%set_wet(wet)
      end if
      allocate (wd%area(wd%st%nx, wd%st%ny), source=AREA)
      wd%mapped = map
      if (map) then
         call wd%st%enter_data()
         call wd%fx%enter_data()
         call wd%rad%enter_data()
         call wd%wk%enter_data()
         call map_area(wd%area)
      end if
   end subroutine make_world

   subroutine free_world(wd)
      type(world_t), intent(inout) :: wd

      if (wd%mapped) then
         call unmap_area(wd%area)
         call wd%wk%exit_data()
         call wd%rad%exit_data()
         call wd%fx%exit_data()
         call wd%st%exit_data()
      end if
      wd%mapped = .false.
      if (allocated(wd%area)) deallocate (wd%area)
      call wd%wk%destroy()
      call wd%rad%destroy()
      call wd%fx%destroy()
      call wd%st%destroy()
      call ice_sw_wiring_cleanup()
   end subroutine free_world

   subroutine step(wd, h, opts, thermo, istep, res)
      type(world_t), intent(inout) :: wd
      type(ice_step_hooks_t), intent(in) :: h
      type(ice_step_opts_t), intent(in) :: opts
      logical, intent(in) :: thermo
      integer, intent(in) :: istep
      type(ice_step_result_t), intent(out) :: res

      call ice_step_cice(wd%st, wd%fx, wd%rad, wd%wk, h, opts, DT, thermo, istep, res)
   end subroutine step

   function pack_world(wd, with_work) result(v)
      !! Every state and flux value the sequencer could write (host copy after a device pull); with
      !! `with_work` also the persistent stress.
      type(world_t), intent(inout) :: wd
      logical, intent(in) :: with_work
      real(wp), allocatable :: v(:)

      call pull_state(wd%st, wd%fx)
      call pull_work(wd%wk)
      v = [reshape(wd%st%aicen, [size(wd%st%aicen)]), reshape(wd%st%vicen, [size(wd%st%vicen)]), &
           reshape(wd%st%vsnon, [size(wd%st%vsnon)]), reshape(wd%st%Tsfcn, [size(wd%st%Tsfcn)]), &
           reshape(wd%st%qicen, [size(wd%st%qicen)]), reshape(wd%st%qsnon, [size(wd%st%qsnon)]), &
           reshape(wd%st%sicen, [size(wd%st%sicen)]), &
           reshape(wd%st%aice, [size(wd%st%aice)]), reshape(wd%st%vice, [size(wd%st%vice)]), &
           reshape(wd%st%vsno, [size(wd%st%vsno)]), reshape(wd%st%aice0, [size(wd%st%aice0)]), &
           reshape(wd%st%aicen_init, [size(wd%st%aicen_init)]), reshape(wd%st%vicen_init, [size(wd%st%vicen_init)]), &
           reshape(wd%st%aice_init, [size(wd%st%aice_init)]), real(reshape(wd%st%ice_err, [size(wd%st%ice_err)]), wp), &
           reshape(wd%fx%fresh, [size(wd%fx%fresh)]), reshape(wd%fx%fsalt, [size(wd%fx%fsalt)]), &
           reshape(wd%fx%fhocn, [size(wd%fx%fhocn)]), reshape(wd%fx%fswthru, [size(wd%fx%fswthru)]), &
           reshape(wd%fx%flwout, [size(wd%fx%flwout)]), reshape(wd%fx%freshn, [size(wd%fx%freshn)]), &
           reshape(wd%fx%fhocnn, [size(wd%fx%fhocnn)])]
      if (with_work) v = [v, reshape(wd%wk%tau_x, [size(wd%wk%tau_x)]), reshape(wd%wk%tau_y, [size(wd%wk%tau_y)])]
   end function pack_world

   function has_text(res, needle) result(found)
      type(ice_step_result_t), intent(in) :: res
      character(len=*), intent(in) :: needle
      logical :: found

      found = index(res%message, needle) > 0
   end function has_text

   ! ------------------------------------------------------------------
   ! the manifest
   ! ------------------------------------------------------------------

   subroutine test_plan_literal(error)
      !! The stage list is the one of the design, written out here.
      type(error_type), allocatable, intent(out) :: error
      type(ice_step_hooks_t) :: h
      type(ice_step_plan_t) :: plan
      character(len=16), parameter :: WANT(17) = [character(len=16) :: "ocean_import", "velocity_to_T", "atm_prepare", &
                    "flux_zero", "save_init", "radiation", "therm1", "therm2", "halo_state", "aggregate_ghost", "stress_to_faces", &
                                                "dynamics", "transport", "ridging", "halo_state_post", "aggregate", "ocean_handoff"]
      integer, parameter :: GROUP(17) = [1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 5, 5, 5, 6]
      logical, parameter :: HOOK(17) = [.true., .true., .true., .false., .false., .false., .false., .false., .true., &
                                        .false., .true., .true., .true., .true., .true., .false., .true.]
      logical, parameter :: ENFORCE(17) = [.false., .false., .false., .false., .false., .false., .false., .true., &
                                           .false., .false., .false., .true., .true., .false., .false., .true., .false.]
      integer :: k

      call all_doubles(h, .true.)
      plan = ice_step_plan_build(h, .true.)
      call check(error, plan%n, 17, "seventeen stages")
      if (allocated(error)) return
      call check(error, ICE_NSTAGE, 17, "ICE_NSTAGE")
      if (allocated(error)) return
      do k = 1, 17
         call check(error, trim(ice_stage_name(plan%id(k))), trim(WANT(k)), "stage name at position "//trim(WANT(k)))
         if (allocated(error)) return
         call check(error, plan%group(k), GROUP(k), "group of "//trim(WANT(k)))
         if (allocated(error)) return
         call check(error, plan%is_hook(k) .eqv. HOOK(k), "hook-or-driver of "//trim(WANT(k)))
         if (allocated(error)) return
         call check(error, plan%enforce_after(k) .eqv. ENFORCE(k), "error enforcement point at "//trim(WANT(k)))
         if (allocated(error)) return
         call check(error, plan%present(k), "present: "//trim(WANT(k)))
         if (allocated(error)) return
         call check(error, plan%required(k), "required (highfreq on): "//trim(WANT(k)))
         if (allocated(error)) return
      end do
      call check(error, plan%group(1) == ICE_GROUP_PREP .and. plan%group(6) == ICE_GROUP_THERMO .and. &
                 plan%group(12) == ICE_GROUP_DYNAMICS .and. plan%group(13) == ICE_GROUP_TRANSPORT .and. &
                 plan%group(14) == ICE_GROUP_CLOSE .and. plan%group(17) == ICE_GROUP_OUT, "group constants")
      if (allocated(error)) return
      ! with no hook installed only the driver stages are present
      call h%reset()
      plan = ice_step_plan_build(h, .false.)
      do k = 1, 17
         call check(error, plan%present(k) .eqv. (.not. HOOK(k)), "default slot is 'missing': "//trim(WANT(k)))
         if (allocated(error)) return
      end do
      call check(error,.not. plan%required(ICE_STAGE_VELOCITY_TO_T), "velocity_to_T is not required without highfreq")
   end subroutine test_plan_literal

   subroutine test_plan_predicates(error)
      !! The binding ordering facts, each as a predicate on positions in the manifest.
      type(error_type), allocatable, intent(out) :: error
      type(ice_step_hooks_t) :: h
      type(ice_step_plan_t) :: plan
      integer :: pos(ICE_NSTAGE), k

      call all_doubles(h, .true.)
      plan = ice_step_plan_build(h, .true.)
      do k = 1, plan%n
         pos(plan%id(k)) = k
      end do
      ! (a) zero < save_init < radiation < therm1 < therm2
      call check(error, pos(ICE_STAGE_FLUX_ZERO) < pos(ICE_STAGE_SAVE_INIT), "zero < save_init")
      if (allocated(error)) return
      call check(error, pos(ICE_STAGE_SAVE_INIT) < pos(ICE_STAGE_RADIATION), "save_init < radiation")
      if (allocated(error)) return
      call check(error, pos(ICE_STAGE_RADIATION) < pos(ICE_STAGE_THERM1), "radiation < therm1")
      if (allocated(error)) return
      call check(error, pos(ICE_STAGE_THERM1) < pos(ICE_STAGE_THERM2), "therm1 < therm2")
      if (allocated(error)) return
      ! every input is in before the first write
      call check(error, pos(ICE_STAGE_OCEAN_IMPORT) < pos(ICE_STAGE_FLUX_ZERO) .and. &
                 pos(ICE_STAGE_VELOCITY_TO_T) < pos(ICE_STAGE_FLUX_ZERO) .and. &
                 pos(ICE_STAGE_ATM_PREPARE) < pos(ICE_STAGE_FLUX_ZERO), "inputs before the first write")
      if (allocated(error)) return
      ! thermodynamics before dynamics; halo then ghost aggregate before anything reads the ring
      call check(error, pos(ICE_STAGE_THERM2) < pos(ICE_STAGE_HALO_STATE), "thermodynamics before the halo")
      if (allocated(error)) return
      call check(error, pos(ICE_STAGE_HALO_STATE) < pos(ICE_STAGE_AGGREGATE_GHOST), "halo < ghost aggregate")
      if (allocated(error)) return
      call check(error, pos(ICE_STAGE_AGGREGATE_GHOST) < pos(ICE_STAGE_STRESS_TO_FACES), "ghost aggregate < stress_to_faces")
      if (allocated(error)) return
      call check(error, pos(ICE_STAGE_STRESS_TO_FACES) < pos(ICE_STAGE_DYNAMICS), "stress_to_faces < dynamics")
      if (allocated(error)) return
      ! the strength is the dynamics' own, from the haloed post-thermodynamics state
      call check(error, plan%dynamics_computes_strength, "the dynamics computes the strength")
      if (allocated(error)) return
      call check(error, trim(plan%strength_routine), "ice_compute_strength", "and the routine is ice_compute_strength")
      if (allocated(error)) return
      call check(error, pos(ICE_STAGE_HALO_STATE) < pos(ICE_STAGE_DYNAMICS) .and. &
                 pos(ICE_STAGE_AGGREGATE_GHOST) < pos(ICE_STAGE_DYNAMICS), "the strength sees the haloed, aggregated state")
      if (allocated(error)) return
      ! (c) dynamics < transport < ridging < halo < aggregate, nothing between transport and ridging
      call check(error, pos(ICE_STAGE_DYNAMICS) < pos(ICE_STAGE_TRANSPORT), "dynamics < transport")
      if (allocated(error)) return
      call check(error, pos(ICE_STAGE_RIDGING) - pos(ICE_STAGE_TRANSPORT) == 1, "ridging directly after transport")
      if (allocated(error)) return
      call check(error, pos(ICE_STAGE_RIDGING) < pos(ICE_STAGE_HALO_STATE_POST), "ridging < halo_post")
      if (allocated(error)) return
      call check(error, pos(ICE_STAGE_HALO_STATE_POST) < pos(ICE_STAGE_AGGREGATE), "halo_post < final aggregate")
      if (allocated(error)) return
      call check(error, plan%transport_post, ICE_TRANSPORT_POST_NONE, "transport post mode is NONE")
      if (allocated(error)) return
      call check(error, plan%transport_post /= ICE_TRANSPORT_POST_AGGREGATE, "transport does not aggregate")
      if (allocated(error)) return
      call check(error, plan%ridging_uses_advected_aice0, "ridging consumes the advected open water")
      if (allocated(error)) return
      ! (e) the ocean hand-off is last, the import first
      call check(error, pos(ICE_STAGE_OCEAN_HANDOFF) == ICE_NSTAGE .and. pos(ICE_STAGE_OCEAN_IMPORT) == 1, &
                 "import first, hand-off last")
      if (allocated(error)) return
      ! (f) the error word is enforced after therm2, the dynamics, transport (alone) and the final aggregate
      call check(error, plan%enforce_after(pos(ICE_STAGE_THERM2)) .and. plan%enforce_after(pos(ICE_STAGE_DYNAMICS)) &
                 .and. plan%enforce_after(pos(ICE_STAGE_TRANSPORT)) .and. plan%enforce_after(pos(ICE_STAGE_AGGREGATE)), &
                 "enforcement after thermo, dynamics, transport, close")
      if (allocated(error)) return
      call check(error, count(plan%enforce_after), 4, "and nowhere else")
   end subroutine test_plan_predicates

   ! ------------------------------------------------------------------
   ! refusals
   ! ------------------------------------------------------------------

   subroutine test_missing_refuses(error)
      !! A missing required stage refuses with its name; nothing is called, nothing is written.
      type(error_type), allocatable, intent(out) :: error
      type(world_t) :: wd
      type(ice_step_hooks_t) :: h
      type(ice_step_opts_t) :: opts
      type(ice_step_result_t) :: res
      real(wp), allocatable :: before(:), after(:)

      call make_world(wd, NG, .true., .true., error)
      if (allocated(error)) then
         call free_world(wd)
         return
      end if
      call reset_log()
      call h%reset()
      call h%install(ICE_STAGE_ATM_PREPARE, dbl_atm, is_double=.true.)
      call h%install(ICE_STAGE_HALO_STATE, dbl_halo, is_double=.true.)
      call h%install(ICE_STAGE_STRESS_TO_FACES, dbl_stress, is_double=.true.)
      call h%install(ICE_STAGE_DYNAMICS, dbl_dynamics, is_double=.true.)
      call h%install(ICE_STAGE_TRANSPORT, dbl_transport, is_double=.true.)
      call h%install(ICE_STAGE_HALO_STATE_POST, dbl_halo_post, is_double=.true.)
      call h%install(ICE_STAGE_OCEAN_HANDOFF, dbl_handoff, is_double=.true.)
      ! ocean_import and ridging are missing
      opts%allow_test_doubles = .true.
      before = pack_world(wd, .true.)
      call step(wd, h, opts, .true., 1, res)
      after = pack_world(wd, .true.)
      call check(error, res%status, ICE_STEP_REFUSED_MISSING, "refused: a required stage is missing")
      if (.not. allocated(error)) call check(error, res%stage, ICE_STAGE_OCEAN_IMPORT, "first missing stage in plan order")
      if (.not. allocated(error)) call check(error, has_text(res, "ocean_import") .and. has_text(res, "ridging"), &
                                             "both missing stages are named: "//trim(res%message))
      if (.not. allocated(error)) call check(error, nlog, 0, "no hook was called")
      if (.not. allocated(error)) call check(error, res%n_exec, 0, "no stage executed")
      if (.not. allocated(error)) call check(error, all(before == after), "state and fluxes bitwise untouched")
      if (.not. allocated(error)) call check(error, ICE_HOOK_MISSING, -1, "default hook status")
      call free_world(wd)
   end subroutine test_missing_refuses

   subroutine test_doubles_refused(error)
      type(error_type), allocatable, intent(out) :: error
      type(world_t) :: wd
      type(ice_step_hooks_t) :: h
      type(ice_step_opts_t) :: opts
      type(ice_step_result_t) :: res
      real(wp), allocatable :: before(:), after(:)

      call make_world(wd, NG, .true., .true., error)
      if (allocated(error)) then
         call free_world(wd)
         return
      end if
      call reset_log()
      call all_doubles(h, .false.)
      before = pack_world(wd, .true.)
      call step(wd, h, opts, .false., 1, res)
      after = pack_world(wd, .true.)
      call check(error, res%status, ICE_STEP_REFUSED_DOUBLE, "a test double is refused by default")
      if (.not. allocated(error)) call check(error, nlog, 0, "nothing called")
      if (.not. allocated(error)) call check(error, all(before == after), "state untouched")
      if (.not. allocated(error)) then
         opts%allow_test_doubles = .true.
         call step(wd, h, opts, .false., 1, res)
         call check(error, res%status, ICE_STEP_OK, "accepted when the test switch is on")
      end if
      call free_world(wd)
   end subroutine test_doubles_refused

   subroutine test_envelope(error)
      !! Every envelope violation is refused before a state or flux value is written.
      type(error_type), allocatable, intent(out) :: error
      type(world_t) :: wd
      type(ice_step_hooks_t) :: h
      type(ice_step_opts_t) :: opts
      type(ice_step_result_t) :: res
      real(wp), allocatable :: before(:), after(:)
      integer :: icase

      call make_world(wd, NG, .true., .true., error)
      if (allocated(error)) then
         call free_world(wd)
         return
      end if
      call reset_log()
      call all_doubles(h, .false.)
      before = pack_world(wd, .false.)
      do icase = 1, 4
         opts = ice_step_opts_t()
         opts%allow_test_doubles = .true.
         select case (icase)
         case (1)
            opts%ndtd = 2
         case (2)
            opts%nranks = 2
         case (3)
            opts%dt_outer = 0.5_wp*DT          ! the ice step is not the outer step
         case default
            opts%dt_outer = -1.0_wp
         end select
         if (icase == 4) then
            call ice_step_cice(wd%st, wd%fx, wd%rad, wd%wk, h, opts, ieee_value(1.0_wp, ieee_quiet_nan), .false., 1, res)
         else
            call step(wd, h, opts, .false., 1, res)
         end if
         call check(error, res%status, ICE_STEP_REFUSED_ENVELOPE, "envelope refusal, case "//achar(48 + icase))
         if (allocated(error)) exit
      end do
      if (.not. allocated(error)) then
         after = pack_world(wd, .false.)
         call check(error, all(before == after), "state untouched by every refusal")
      end if
      if (.not. allocated(error)) call check(error, nlog, 0, "no hook called")
      ! non-positive forcing height: refused after the input hooks, still before the first state or flux write
      if (.not. allocated(error)) then
         opts = ice_step_opts_t()
         opts%allow_test_doubles = .true.
         call h%install(ICE_STAGE_ATM_PREPARE, dbl_atm_bad_height, is_double=.true.)
         call step(wd, h, opts, .false., 1, res)
         call check(error, res%status, ICE_STEP_REFUSED_ENVELOPE, "zlvl <= 0 refused")
         if (.not. allocated(error)) call check(error, res%n_exec, 2, "after the input hooks, before any write")
         if (.not. allocated(error)) then
            after = pack_world(wd, .false.)
            call check(error, all(before == after), "state and fluxes untouched by that refusal")
         end if
      end if
      call free_world(wd)
      if (allocated(error)) return

      ! nghost = 2: refused (the dynamics needs 3)
      call make_world(wd, 2, .true., .true., error)
      if (.not. allocated(error)) then
         opts = ice_step_opts_t()
         opts%allow_test_doubles = .true.
         call all_doubles(h, .false.)
         call step(wd, h, opts, .false., 1, res)
         call check(error, res%status, ICE_STEP_REFUSED_ENVELOPE, "nghost = 2 refused")
         if (.not. allocated(error)) call check(error, index(res%message, "nghost") > 0, "the message names nghost")
      end if
      call free_world(wd)
      if (allocated(error)) return

      ! unmapped slots: refused
      call make_world(wd, NG, .false., .true., error)
      if (.not. allocated(error)) then
         call step(wd, h, opts, .false., 1, res)
         call check(error, res%status, ICE_STEP_REFUSED_ENVELOPE, "unmapped slots refused")
         if (.not. allocated(error)) call check(error, index(res%message, "mapped") > 0, "the message says unmapped")
      end if
      call free_world(wd)
      if (allocated(error)) return

      ! wet mask never set: refused (no default of all wet)
      call make_world(wd, NG, .true., .false., error)
      if (.not. allocated(error)) then
         opts = ice_step_opts_t()
         opts%allow_test_doubles = .true.
         call all_doubles(h, .false.)
         call reset_log()
         call step(wd, h, opts, .false., 1, res)
         call check(error, res%status, ICE_STEP_REFUSED_ENVELOPE, "a wet mask never set is refused")
         if (.not. allocated(error)) call check(error, index(res%message, "wet") > 0, "the message names the wet mask")
         if (.not. allocated(error)) call check(error, nlog, 0, "before any hook")
      end if
      call free_world(wd)
   end subroutine test_envelope

   subroutine dbl_atm_bad_height(ctx, status)
      type(ice_step_ctx_t), intent(inout) :: ctx
      integer, intent(out) :: status

      call dbl_atm(ctx, status)
      ctx%work%zlvl = 0.0_wp
   end subroutine dbl_atm_bad_height

   ! ------------------------------------------------------------------
   ! order, cadence, enforcement
   ! ------------------------------------------------------------------

   subroutine run_steps(wd, h, opts, nsteps, thermo_every, res)
      type(world_t), intent(inout) :: wd
      type(ice_step_hooks_t), intent(in) :: h
      type(ice_step_opts_t), intent(in) :: opts
      integer, intent(in) :: nsteps, thermo_every
      type(ice_step_result_t), intent(out) :: res
      integer :: s

      do s = 1, nsteps
         call step(wd, h, opts, mod(s, thermo_every) == 0, s, res)
         if (res%status /= ICE_STEP_OK) return
      end do
   end subroutine run_steps

   subroutine test_order_cadence(error)
      !! The recording doubles see the manifest's order, every step, with the thermo flag at the right cadence.
      type(error_type), allocatable, intent(out) :: error
      type(world_t) :: wd
      type(ice_step_hooks_t) :: h
      type(ice_step_opts_t) :: opts
      type(ice_step_result_t) :: res
      integer, parameter :: NH = 9
      integer, parameter :: HOOKS_WANT(NH) = [ICE_STAGE_OCEAN_IMPORT, ICE_STAGE_ATM_PREPARE, ICE_STAGE_HALO_STATE, &
                                              ICE_STAGE_STRESS_TO_FACES, ICE_STAGE_DYNAMICS, ICE_STAGE_TRANSPORT, &
                                              ICE_STAGE_RIDGING, ICE_STAGE_HALO_STATE_POST, ICE_STAGE_OCEAN_HANDOFF]
      integer, parameter :: NE = 16
      integer, parameter :: EXEC_WANT(NE) = [ICE_STAGE_OCEAN_IMPORT, ICE_STAGE_ATM_PREPARE, ICE_STAGE_FLUX_ZERO, &
                                             ICE_STAGE_SAVE_INIT, ICE_STAGE_RADIATION, ICE_STAGE_THERM1, ICE_STAGE_THERM2, &
                                             ICE_STAGE_HALO_STATE, ICE_STAGE_AGGREGATE_GHOST, ICE_STAGE_STRESS_TO_FACES, &
                                             ICE_STAGE_DYNAMICS, ICE_STAGE_TRANSPORT, ICE_STAGE_RIDGING, &
                                             ICE_STAGE_HALO_STATE_POST, ICE_STAGE_AGGREGATE, ICE_STAGE_OCEAN_HANDOFF]
      integer :: s, k, m

      call make_world(wd, NG, .true., .true., error)
      if (allocated(error)) then
         call free_world(wd)
         return
      end if
      call reset_log()
      call check(error, wd%st%p%atm%highfreq, .false., "highfreq is off by default")
      if (allocated(error)) then
         call free_world(wd)
         return
      end if
      call all_doubles(h, .true.)         ! velocity_to_T installed but not required: must NOT run
      opts%allow_test_doubles = .true.
      call run_steps(wd, h, opts, 4, 2, res)  ! thermo steps are 2 and 4
      call check(error, res%status, ICE_STEP_OK, "four steps run: "//trim(res%message))
      if (.not. allocated(error)) call check(error, nlog, 4*NH, "nine hooks per step, velocity_to_T not among them")
      if (.not. allocated(error)) then
         call check(error, res%n_exec, NE, "16 stages executed per step")
         if (.not. allocated(error)) call check(error, all(res%executed(1:NE) == EXEC_WANT), &
                                                "the executed order is the manifest order")
      end if
      do s = 1, 4
         if (allocated(error)) exit
         do k = 1, NH
            m = (s - 1)*NH + k
            call check(error, log_stage(m), HOOKS_WANT(k), "hook order in step")
            if (allocated(error)) exit
            call check(error, log_istep(m), s, "step counter reaches the hook")
            if (allocated(error)) exit
            call check(error, log_thermo(m) .eqv. (mod(s, 2) == 0), "thermo flag follows the cadence")
            if (allocated(error)) exit
            call check(error, log_dt(m) == DT, "dt reaches the hook")
            if (allocated(error)) exit
            if (HOOKS_WANT(k) == ICE_STAGE_TRANSPORT) then
               call check(error, log_post(m), ICE_TRANSPORT_POST_NONE, "transport is told to post nothing")
               if (allocated(error)) exit
            end if
         end do
      end do
      call free_world(wd)
   end subroutine test_order_cadence

   subroutine test_highfreq(error)
      !! `highfreq` makes velocity_to_T required, among the inputs; missing it refuses.
      type(error_type), allocatable, intent(out) :: error
      type(world_t) :: wd
      type(ice_step_hooks_t) :: h
      type(ice_step_opts_t) :: opts
      type(ice_step_result_t) :: res
      integer, parameter :: WANT(5) = [ICE_STAGE_OCEAN_IMPORT, ICE_STAGE_VELOCITY_TO_T, ICE_STAGE_ATM_PREPARE, &
                                       ICE_STAGE_FLUX_ZERO, ICE_STAGE_SAVE_INIT]

      call make_world(wd, NG, .true., .true., error)
      if (allocated(error)) then
         call free_world(wd)
         return
      end if
      call reset_log()
      wd%st%p%atm%highfreq = .true.
      opts%allow_test_doubles = .true.
      call all_doubles(h, .false.)
      call step(wd, h, opts, .false., 1, res)
      call check(error, res%status, ICE_STEP_REFUSED_MISSING, "highfreq without velocity_to_T refuses")
      if (.not. allocated(error)) call check(error, res%stage, ICE_STAGE_VELOCITY_TO_T, "the missing stage is velocity_to_T")
      if (.not. allocated(error)) then
         call all_doubles(h, .true.)
         call step(wd, h, opts, .false., 1, res)
         call check(error, res%status, ICE_STEP_OK, "runs with it: "//trim(res%message))
      end if
      if (.not. allocated(error)) call check(error, all(res%executed(1:5) == WANT), &
                                             "velocity_to_T sits among the inputs, before the first write")
      if (.not. allocated(error)) call check(error, count(log_stage(1:nlog) == ICE_STAGE_VELOCITY_TO_T), 1, "called once")
      call free_world(wd)
   end subroutine test_highfreq

   subroutine test_fatal(error)
      !! A fatal bit stops the step before the next stage group; the stages of the next group never run.
      type(error_type), allocatable, intent(out) :: error
      type(world_t) :: wd
      type(ice_step_hooks_t) :: h
      type(ice_step_opts_t) :: opts
      type(ice_step_result_t) :: res

      call make_world(wd, NG, .true., .true., error)
      if (allocated(error)) then
         call free_world(wd)
         return
      end if
      opts%allow_test_doubles = .true.
      call all_doubles(h, .false.)

      ! (1) the dynamics raises DYN_VEL: transport, ridging, the closing group and the hand-off never run
      call reset_log()
      inject_stage = ICE_STAGE_DYNAMICS
      inject_bit = ICE_ERR_DYN_VEL
      call step(wd, h, opts, .true., 1, res)
      call check(error, res%status, ICE_STEP_FATAL, "fatal bit stops the step")
      if (.not. allocated(error)) call check(error, res%stage, ICE_STAGE_DYNAMICS, "at the dynamics group")
      if (.not. allocated(error)) call check(error, iand(res%summary%mask, ICE_ERR_DYN_VEL) /= 0_int32, "the summary carries it")
      if (.not. allocated(error)) call check(error, res%summary%first_idx > 0_int64, "and names a cell")
      if (.not. allocated(error)) call check(error, count(log_stage(1:nlog) == ICE_STAGE_TRANSPORT) &
                                             + count(log_stage(1:nlog) == ICE_STAGE_RIDGING) &
                                             + count(log_stage(1:nlog) == ICE_STAGE_OCEAN_HANDOFF), 0, "no later hook ran")
      if (.not. allocated(error)) call check(error, res%executed(res%n_exec), ICE_STAGE_DYNAMICS, &
                                             "the last stage run is the dynamics")

      ! (1b) a fatal bit raised by the ocean import (stage 1) is not wiped by the step-start clear
      if (.not. allocated(error)) then
         call reset_log()
         inject_stage = ICE_STAGE_OCEAN_IMPORT
         inject_bit = ICE_ERR_DYN_VEL
         call step(wd, h, opts, .true., 3, res)
         call check(error, res%status, ICE_STEP_FATAL, "a bit raised by the import hook stops the step")
         if (.not. allocated(error)) call check(error, count(log_stage(1:nlog) == ICE_STAGE_STRESS_TO_FACES), 0, &
                                                "before the dynamics group")
      end if

      ! (2) the REAL thermodynamic stages raise a fatal bit (a NaN shortwave): stopped before the halo
      if (.not. allocated(error)) then
         call reset_log()
         nan_sw = .true.
         call step(wd, h, opts, .true., 2, res)
         call check(error, res%status, ICE_STEP_FATAL, "a NaN shortwave is fatal")
         if (.not. allocated(error)) call check(error, res%stage, ICE_STAGE_THERM2, "reported after the thermodynamics group")
         if (.not. allocated(error)) call check(error, iand(res%summary%mask, ICE_ERR_NONFINITE) /= 0_int32, "as non-finite")
         if (.not. allocated(error)) call check(error, count(log_stage(1:nlog) == ICE_STAGE_HALO_STATE) &
                                                + count(log_stage(1:nlog) == ICE_STAGE_DYNAMICS), 0, &
                                                "dynamics group never ran")
      end if
      call free_world(wd)
   end subroutine test_fatal

   subroutine test_transport_fatal(error)
      !! The error word is enforced after TRANSPORT alone: the ridging never runs on a state transport flagged,
      !! and the reported stage is the transport.
      type(error_type), allocatable, intent(out) :: error
      type(world_t) :: wd
      type(ice_step_hooks_t) :: h
      type(ice_step_opts_t) :: opts
      type(ice_step_result_t) :: res

      call make_world(wd, NG, .true., .true., error)
      if (allocated(error)) then
         call free_world(wd)
         return
      end if
      opts%allow_test_doubles = .true.
      call all_doubles(h, .false.)
      call reset_log()
      inject_stage = ICE_STAGE_TRANSPORT
      inject_bit = ICE_ERR_TRANSPORT_CFL
      call step(wd, h, opts, .false., 1, res)
      call check(error, res%status, ICE_STEP_FATAL, "a transport bit is fatal")
      if (.not. allocated(error)) call check(error, res%stage, ICE_STAGE_TRANSPORT, "and the stage is the transport")
      if (.not. allocated(error)) call check(error, count(log_stage(1:nlog) == ICE_STAGE_RIDGING), 0, "ridging never ran")
      call free_world(wd)
   end subroutine test_transport_fatal

   subroutine test_counted(error)
      !! A COUNTED bit is reported in the mask and does not stop the step.
      type(error_type), allocatable, intent(out) :: error
      type(world_t) :: wd
      type(ice_step_hooks_t) :: h
      type(ice_step_opts_t) :: opts
      type(ice_step_result_t) :: res

      call make_world(wd, NG, .true., .true., error)
      if (allocated(error)) then
         call free_world(wd)
         return
      end if
      opts%allow_test_doubles = .true.
      call all_doubles(h, .false.)
      call reset_log()
      inject_stage = ICE_STAGE_RIDGING
      inject_bit = ICE_ERR_RIDGE_PASSES
      call step(wd, h, opts, .true., 1, res)
      call check(error, res%status, ICE_STEP_OK, "a counted bit does not stop the step")
      if (.not. allocated(error)) call check(error, iand(res%err_mask, ICE_ERR_RIDGE_PASSES) /= 0_int32, "but it is reported")
      if (.not. allocated(error)) call check(error, res%executed(res%n_exec), ICE_STAGE_OCEAN_HANDOFF, &
                                             "the step ran to the hand-off")
      call free_world(wd)
   end subroutine test_counted

   subroutine test_stage_failure(error)
      !! A hook that returns a non-zero status stops the step at that stage.
      type(error_type), allocatable, intent(out) :: error
      type(world_t) :: wd
      type(ice_step_hooks_t) :: h
      type(ice_step_opts_t) :: opts
      type(ice_step_result_t) :: res

      call make_world(wd, NG, .true., .true., error)
      if (allocated(error)) then
         call free_world(wd)
         return
      end if
      opts%allow_test_doubles = .true.
      call all_doubles(h, .false.)
      call reset_log()
      fail_stage = ICE_STAGE_TRANSPORT
      call step(wd, h, opts, .false., 1, res)
      call check(error, res%status, ICE_STEP_STAGE_FAILED, "a failing hook stops the step")
      if (.not. allocated(error)) call check(error, res%stage, ICE_STAGE_TRANSPORT, "at that stage")
      if (.not. allocated(error)) call check(error, res%stage_status, 7, "with its status")
      if (.not. allocated(error)) call check(error, count(log_stage(1:nlog) == ICE_STAGE_RIDGING), 0, "ridging never ran")
      call free_world(wd)
   end subroutine test_stage_failure

   subroutine test_halo_ghost(error)
      !! The halo hook delivers a ghost ring that carries ice; the sequencer's ghosted aggregate makes the
      !! ghost aggregates consistent with it before anything reads them (and again at the end of the step).
      type(error_type), allocatable, intent(out) :: error
      type(world_t) :: wd
      type(ice_step_hooks_t) :: h
      type(ice_step_opts_t) :: opts
      type(ice_step_result_t) :: res
      real(wp), allocatable :: va(:, :)

      call make_world(wd, NG, .true., .true., error)
      if (allocated(error)) then
         call free_world(wd)
         return
      end if
      opts%allow_test_doubles = .true.
      call all_doubles(h, .false.)
      call reset_log()
      scn_ghost = .true.
      call step(wd, h, opts, .false., 1, res)
      call check(error, res%status, ICE_STEP_OK, "step runs: "//trim(res%message))
      if (.not. allocated(error)) then
         call pull_state(wd%st, wd%fx)
         va = sum(wd%st%aicen, dim=3)
         call check(error, wd%st%aice(1, 1) > 0.4_wp, "the ghost aggregate was rebuilt from the delivered ghost state")
         if (.not. allocated(error)) call check(error, all(abs(va - wd%st%aice) < 1.0e-14_wp), &
                                                "aggregate area equals the category sum on the whole ghosted array")
         if (.not. allocated(error)) call check(error, all(abs(wd%st%aice0 - max(1.0_wp - wd%st%aice, 0.0_wp)) < 1.0e-14_wp), &
                                                "and so does the open water")
      end if
      call free_world(wd)
   end subroutine test_halo_ghost

   ! ------------------------------------------------------------------
   ! contracts: the bank, the stress
   ! ------------------------------------------------------------------

   subroutine test_bank_once(error)
      !! The bank is formed into `frzmlt` every step from the CURRENT ocean-side bank, spent by the new-ice stage,
      !! and written back by the hand-off hook: over a 3-step thermo window it is frozen exactly once.
      type(error_type), allocatable, intent(out) :: error
      type(world_t) :: wd
      type(ice_step_hooks_t) :: h
      type(ice_step_opts_t) :: opts
      type(ice_step_result_t) :: res
      real(wp), parameter :: BANK0 = 5.0e5_wp
      real(wp) :: vol1, expect
      integer :: s

      call make_world(wd, NG, .true., .true., error)
      if (allocated(error)) then
         call free_world(wd)
         return
      end if
      opts%allow_test_doubles = .true.
      call all_doubles(h, .false.)
      call reset_log()
      wd%wk%bank_ocean(IC, JC) = BANK0
      call push_bank_ocean(wd%wk)
      do s = 1, 3
         call step(wd, h, opts, s == 3, s, res)
         call check(error, res%status, ICE_STEP_OK, "step runs: "//trim(res%message))
         if (allocated(error)) exit
         call pull_state(wd%st, wd%fx)
         call pull_work(wd%wk)
         if (s == 1) then
            vol1 = wd%fx%frazil_vol(IC, JC)
            expect = BANK0/(wd%p%c%rhoi*wd%p%c%Lfresh)
            call check(error, vol1 > 0.0_wp, "new ice formed in the first step")
            if (.not. allocated(error)) call check(error, abs(vol1 - expect) <= 1.0e-6_wp*expect, &
                                                   "its volume is the bank's (bank/(rhoi Lfresh))")
            if (.not. allocated(error)) call check(error, wd%wk%bank(IC, JC) == 0.0_wp, "the new-ice stage spent the bank")
            if (.not. allocated(error)) call check(error, wd%wk%bank_ocean(IC, JC) == 0.0_wp, &
                                                   "the hand-off wrote the spent bank back to the ocean side")
         else
            call check(error, wd%wk%frzmlt(IC, JC) == 0.0_wp, "no freezing potential from a spent bank")
            if (.not. allocated(error)) call check(error, sum(abs(wd%fx%frazil_vol)) == 0.0_wp, &
                                                   "no new ice from it in a later step of the window")
         end if
         if (allocated(error)) exit
      end do
      call free_world(wd)
   end subroutine test_bank_once

   subroutine test_stress_restart(error)
      !! The ice-ocean stress is persistent state: a run restarted from a file in the middle resumes with the same
      !! stress input the thermodynamics reads, and ends bitwise equal to the straight run.
      type(error_type), allocatable, intent(out) :: error
      type(world_t), target :: wa, wb
      type(ice_step_hooks_t) :: h
      type(ice_step_opts_t) :: opts
      type(ice_step_result_t) :: res
      type(restart_registry_t) :: ra, rb
      type(decomp_t) :: d
      type(ocean_restart_metadata_t) :: meta
      character(len=*), parameter :: FN = "test_cice_step_plan_restart.nc"
      real(wp), allocatable :: va(:), vb(:)
      real(wp) :: t
      integer :: s, ierr, step_r, osc, lun

      call make_world(wa, NG, .true., .true., error)
      if (allocated(error)) then
         call free_world(wa)
         return
      end if
      opts%allow_test_doubles = .true.
      call all_doubles(h, .false.)
      call reset_log()
      do s = 1, 3
         call step(wa, h, opts, .false., s, res)
         call check(error, res%status, ICE_STEP_OK, "first half runs")
         if (allocated(error)) then
            call free_world(wa)
            return
         end if
      end do
      call pull_state(wa%st, wa%fx)      ! the writer is host I/O: pull first
      call pull_work(wa%wk)
      call wa%st%register_restart(ra)
      call wa%wk%register_restart(ra)
      call decomp_init(d, NXP, NYP, 1, 1, 0)
      meta%nghost = NG
      open (newunit=lun, file=FN, status="replace")
      close (lun, status="delete")
      call ocean_restart_write_local(FN, ra, d, meta, 3.0_wp*DT, 3, 0, ierr)
      call check(error, ierr, 0, "restart write")
      if (allocated(error)) then
         call free_world(wa)
         return
      end if
      do s = 4, 6
         call step(wa, h, opts, .false., s, res)
         call check(error, res%status, ICE_STEP_OK, "second half runs (straight)")
         if (allocated(error)) then
            call free_world(wa)
            return
         end if
      end do
      va = pack_world(wa, .true.)

      ! the restarted world: a fresh seed overwritten by the file, then the same last three steps
      call make_world(wb, NG, .true., .true., error)
      if (allocated(error)) then
         call free_world(wa)
         call free_world(wb)
         return
      end if
      wb%st%aicen = 0.0_wp
      wb%wk%tau_x = 0.0_wp
      wb%wk%tau_y = 0.0_wp
      call wb%st%register_restart(rb)
      call wb%wk%register_restart(rb)
      call ocean_restart_read_local(FN, rb, d, meta, t, step_r, osc, ierr)
      call check(error, ierr, 0, "restart read")
      if (.not. allocated(error)) then
         call push_restart(wb%st, wb%wk)
         call ice_update_state(wb%st, .true.)       ! the aggregates are not restart state
         do s = 4, 6
            call step(wb, h, opts, .false., s, res)
            call check(error, res%status, ICE_STEP_OK, "second half runs (restarted)")
            if (allocated(error)) exit
         end do
      end if
      if (.not. allocated(error)) then
         vb = pack_world(wb, .true.)
         call check(error, size(va), size(vb), "same shape")
         if (.not. allocated(error)) call check(error, all(va == vb), &
                                                "the restarted run equals the straight run bitwise (state, fluxes, stress)")
      end if
      open (newunit=lun, file=FN, status="old")
      close (lun, status="delete")
      call free_world(wa)
      call free_world(wb)
   end subroutine test_stress_restart

   ! ------------------------------------------------------------------
   ! smoke runs
   ! ------------------------------------------------------------------

   subroutine test_smoke_dark(error)
      type(error_type), allocatable, intent(out) :: error

      call run_smoke(.false., error)
   end subroutine test_smoke_dark

   subroutine test_smoke_daylight(error)
      type(error_type), allocatable, intent(out) :: error

      call run_smoke(.true., error)
   end subroutine test_smoke_daylight

   subroutine run_smoke(daylight, error)
      !! 50 steps with the real radiation, therm1 and therm2 and no-motion doubles: clean, and the shipped ledger
      !! closes (water, salt, heat) on every step.  The daylight case adds shortwave, a frazil bank on the odd
      !! columns, a melt potential with a warm ocean on the even ones, snowfall (steps 1-10) and rain (11-20).
      logical, intent(in) :: daylight
      type(error_type), allocatable, intent(out) :: error
      type(world_t) :: wd
      type(ice_step_hooks_t) :: h
      type(ice_step_opts_t) :: opts
      type(ice_step_result_t) :: res
      type(ice_ledger_mirror_t) :: led
      integer, parameter :: NSTEP = 50
      integer :: s, i, j
      real(wp) :: worst_w, worst_s, worst_h, worst_ratio, w_ext, h_ext, bound, rw, rs, rh, vice0, vice1, e_frz, aref
      real(wp) :: fthr_max, meltl_sum, frazil_sum, bank_start, bank_end, a_sum
      real(wp), allocatable :: va(:, :)

      call make_world(wd, NG, .true., .true., error)
      if (allocated(error)) then
         call free_world(wd)
         return
      end if
      opts%allow_test_doubles = .true.
      opts%dt_outer = DT
      call all_doubles(h, .false.)
      call reset_log()
      bank_start = 0.0_wp
      if (daylight) then
         scn_sw = 250.0_wp
         scn_melt = .true.
         scn_precip = .true.
         wd%st%p%thermo%update_ocn_f = .true.     ! the new ice's mass leaves the ocean through `fresh`
         do j = NG + 1, NG + NYP
            do i = NG + 1, NG + NXP
               if (mod(i - NG, 2) == 1) then
                  wd%wk%bank_ocean(i, j) = 5.0e5_wp
                  bank_start = bank_start + wd%wk%bank_ocean(i, j)
               end if
            end do
         end do
         call push_bank_ocean(wd%wk)
      end if
      call pull_state(wd%st, wd%fx)
      vice0 = sum(wd%st%vice)
      worst_w = 0.0_wp
      worst_s = 0.0_wp
      worst_h = 0.0_wp
      worst_ratio = 0.0_wp
      fthr_max = 0.0_wp
      meltl_sum = 0.0_wp
      frazil_sum = 0.0_wp
      do s = 1, NSTEP
         call led%open(wd%st, wd%area, DT, 1)
         call step(wd, h, opts, mod(s, 8) == 0, s, res)
         call check(error, res%status, ICE_STEP_OK, "step runs: "//trim(res%message))
         if (allocated(error)) exit
         call check(error, res%err_mask == 0_int32, "no error bit, counted or fatal")
         if (allocated(error)) exit
         call pull_ledger_inputs(wd%st, wd%wk)
         ! the external sources of the ledger over this step (the forcing the doubles set, read back)
         a_sum = sum(wd%st%aice_init(NG + 1:NG + NXP, NG + 1:NG + NYP))
         aref = a_sum
         w_ext = DT*AREA*sum(wd%wk%fsnow(NG + 1:NG + NXP, NG + 1:NG + NYP) &
                             *wd%st%aice_init(NG + 1:NG + NXP, NG + 1:NG + NYP) &
                             + wd%wk%frain(NG + 1:NG + NXP, NG + 1:NG + NYP)*wd%st%aice_init(NG + 1:NG + NXP, NG + 1:NG + NYP))
         h_ext = DT*AREA*sum((wd%wk%flw(NG + 1:NG + NXP, NG + 1:NG + NYP) &
                              - wd%p%c%Lfresh*wd%wk%fsnow(NG + 1:NG + NXP, NG + 1:NG + NYP)) &
                             *wd%st%aice_init(NG + 1:NG + NXP, NG + 1:NG + NYP))
         ! the bank energy spent by new ice (the bank formed into a positive freezing potential)
         e_frz = 0.0_wp
         do j = NG + 1, NG + NYP
            do i = NG + 1, NG + NXP
               if (wd%wk%frzmlt(i, j) > 0.0_wp) e_frz = e_frz + AREA*wd%wk%frzmlt(i, j)*DT
            end do
         end do
         h_ext = h_ext - e_frz      ! the new ice's enthalpy is NEGATIVE: freezing the bank lowers the stored energy
         call led%close(wd%st, wd%fx, wd%area, w_ext, h_ext)
         if (.not. led%valid) then
            call test_failed(error, "the ledger did not close")
            exit
         end if
         rw = led%rel_residual(ICE_LEDGER_WATER)
         rs = led%rel_residual(ICE_LEDGER_SALT)
         rh = led%rel_residual(ICE_LEDGER_HEAT)
         if (.not. (ieee_is_finite(rw) .and. ieee_is_finite(rs) .and. ieee_is_finite(rh))) then
            call test_failed(error, "non-finite ledger residual")
            exit
         end if
         worst_w = max(worst_w, rw)
         worst_s = max(worst_s, rs)
         worst_h = max(worst_h, rh)
         bound = 1.1_wp*ICE_FERRMAX*DT*AREA*aref + 1.0e-15_wp*led%scale(ICE_LEDGER_HEAT)
         worst_ratio = max(worst_ratio, abs(led%res(ICE_LEDGER_HEAT))/bound)
         call pull_state(wd%st, wd%fx)
         fthr_max = max(fthr_max, maxval(wd%fx%fswthru(NG + 1:NG + NXP, NG + 1:NG + NYP)))
         meltl_sum = meltl_sum + sum(wd%fx%meltl(NG + 1:NG + NXP, NG + 1:NG + NYP))
         frazil_sum = frazil_sum + sum(wd%fx%frazil_vol(NG + 1:NG + NXP, NG + 1:NG + NYP))
      end do
      if (.not. allocated(error)) then
         call pull_state(wd%st, wd%fx)
         call pull_work(wd%wk)
         vice1 = sum(wd%st%vice)
         bank_end = sum(wd%wk%bank_ocean)
         write (output_unit, '(a,l1,a,3es11.3,a,es10.3,a,es10.3)') "step smoke (daylight=", daylight, &
            "), 50 steps, worst relative residual water/salt/heat: ", worst_w, worst_s, worst_h, &
            "; worst heat residual / BL99 bound ", worst_ratio, "; total ice volume change (m) ", vice1 - vice0
         call check(error, abs(vice1 - vice0) > 1.0e-3_wp*abs(vice0), "the ice actually evolved (not a vacuous run)")
         if (.not. allocated(error)) call check(error, worst_w < 1.0e-9_wp, "water closes")
         if (.not. allocated(error)) call check(error, worst_s < 1.0e-9_wp, "salt closes")
         if (.not. allocated(error)) call check(error, worst_ratio <= 1.0_wp, "heat closes within the BL99 solver bound")
         if (daylight .and. .not. allocated(error)) then
            write (output_unit, '(a,es10.3,a,es10.3,a,es10.3)') "  daylight: max fswthru ", fthr_max, &
               " W/m2, lateral melt (m, summed) ", meltl_sum, ", new ice (m, summed) ", frazil_sum
            call check(error, fthr_max > 0.0_wp, "shortwave went through the ice (fswthru > 0)")
            if (.not. allocated(error)) call check(error, bank_start > 0.0_wp .and. bank_end == 0.0_wp, &
                                                   "the whole bank was spent (and written back as zero)")
            if (.not. allocated(error)) call check(error, frazil_sum > 0.0_wp, "new ice formed from it")
            if (.not. allocated(error)) call check(error, meltl_sum > 0.0_wp, "lateral melt ran")
         end if
         ! every aggregate is consistent with the categories over the whole ghosted array
         if (.not. allocated(error)) then
            allocate (va(wd%st%nx, wd%st%ny))
            va = sum(wd%st%aicen, dim=3)
            call check(error, all(abs(va - wd%st%aice) < 1.0e-14_wp), "aggregate area equals the category sum")
         end if
      end if
      call free_world(wd)
   end subroutine run_smoke

end module test_cice_step_plan
