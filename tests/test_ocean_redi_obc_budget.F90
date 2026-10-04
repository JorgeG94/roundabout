!! Redi neutral diffusion with a Flather OPEN edge must keep the closed salt
!! and heat budgets closed (found by the pairwise compatibility matrix, row
!! `redi_obc_salt_budget`).
!!
!! ## The bug
!!
!! `redi_apply_flux` lets the along-isopycnal flux cross an OPEN physical
!! face, read against the OBC-filled ghost column (MOM6 does the same:
!! `neutral_diffusion` gates its faces on `G%mask2dCu`, which stays 1 on an
!! open segment's normal face).  That exchange changes the domain's salt
!! and heat content, but it was booked in no budget accumulator, so the
!! console's closed residual `((Q - Q0) + out - src)/Q0` would have printed
!! the Redi boundary flux as a leak.  The driver avoided that by demoting the
!! WHOLE Salt/Heat `Error` column to raw drift whenever Redi met an open edge
!! (`redi_with_open_edge`).  Raw drift cannot subtract the ADVECTIVE exchange
!! through the open face, so every such run printed a -4e-5 "leak" in 24
!! steps.  That happened even with `khtr = 0`.
!!
!! Fix: `redi_apply_flux` books its realised increment into
!! `ms%salt_budget_hdiff` / `heat_budget_hdiff` (the accumulator
!! `tracer_hdiff` already uses), and the fall-back gate is removed.
!!
!! ## The test
!!
!! A wind-driven stratified channel over a seamount (sigma, so isopycnals
!! cross the layers and Redi acts), open (Flather) east edge, walls
!! elsewhere.  After `N_STEPS` it forms the closed residual exactly as the
!! console does, from the model's own budget (`budget_out`) and the salt and
!! heat totals.  The residual must be at round-off.  Before the fix,
!! the Redi exchange through the east face was missing from `out`. The
!! residual was then 2.0e-6 for salt and 6e-5 to 8e-5 for heat on both
!! schemes. Both split schemes run, because the budget stage weight differs
!! between them.
module test_ocean_redi_obc_budget
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_constants, only: wp, RHO_WATER
   use rdb_config, only: config_t, read_config_from_string, validate_config
   use rdb_ocean_engine, only: ocean_engine_t, engine_setup, engine_enter_data, &
                               engine_step, engine_step_finalize, &
                               engine_exit_data, engine_teardown
   use rdb_console_stats, only: console_stats_t, conservation_budget_t
   use rdb_ocean_console_stats, only: ocean_console_stats_report, ocean_budget_stage_weight
   use rdb_ocean_dyn, only: SPLIT_SCHEME_PRED_CORR
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   use rdb_comm_env, only: comm_env_init, comm_env_setup_roles
   implicit none
   private

   public :: collect_ocean_redi_obc_budget_tests

   integer, parameter :: NX_G = 24
   integer, parameter :: NY_G = 12
   integer, parameter :: N_STEPS = 24
   real(wp), parameter :: DT = 900.0_wp
   real(wp), parameter :: RESID_TOL = 1.0e-12_wp
      !! Relative closed-budget residual bound.  The fixed code sits at
      !! round-off here, while the missing Redi open-face term was 2e-6 to 8e-5.
   real(wp), parameter :: REDI_SIGNAL_MIN = 1.0e-11_wp
      !! Positive control: Redi ON vs OFF must move the salt total by at
      !! least this relative amount, otherwise the test proves nothing.

   logical :: comm_inited = .false.

contains

   subroutine collect_ocean_redi_obc_budget_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("redi_open_edge_budget_closes_pred_corr", test_pred_corr), &
                  new_unittest("redi_open_edge_budget_closes_ssp_rk2", test_ssp_rk2) &
                  ]
   end subroutine collect_ocean_redi_obc_budget_tests

   subroutine ensure_comm()
      !! On an `RDB_ENABLE_MPI=ON` build the first collective in
      !! `engine_setup` hits `MPI_Comm_f2c` before `MPI_Init` without this.
      if (.not. comm_inited) then
         call comm_env_init()
         call comm_env_setup_roles(.false.)
         comm_inited = .true.
      end if
   end subroutine ensure_comm

   function channel_nml(scheme, khtr) result(nml)
      !! Stratified seamount channel, Flather open east edge, Redi on.
      character(len=*), intent(in) :: scheme
      character(len=*), intent(in) :: khtr
      character(len=:), allocatable :: nml
      character(len=*), parameter :: NL = new_line("a")
      character(len=16) :: snx, sny

      write (snx, '(i0)') NX_G
      write (sny, '(i0)') NY_G
      nml = "&sim_nml sim_type = 'ocean' /"//NL// &
            "&time_nml t_end = 21600.0, dt_fixed = 900.0 /"//NL// &
            "&nonhydrostatic_nml nz_layers = 6 /"//NL// &
            "&tracer_nml initial_temperature = 12.0, initial_salinity = 35.0, "// &
            "T_init_surface = 20.0, T_init_bottom = 4.0, "// &
            "S_init_surface = 34.5, S_init_bottom = 35.2 /"//NL// &
            "&ocean_thermo_nml enable_thermodynamics = .true. /"//NL// &
            "&ocean_eos_nml eos = 'wright' /"//NL// &
            "&ocean_bt_nml auto_n_inner = .true., split_scheme = '"//scheme//"' /"//NL// &
            "&ocean_hvisc_nml nu_h = 200.0, lateral_closure = 'smagorinsky', "// &
            "smag_ah = .true. /"//NL// &
            "&ocean_diag_nml enabled = .false. /"//NL// &
            "&grid_nml nx = "//trim(snx)//", ny = "//trim(sny)//", nghost = 3, "// &
            "dx = 20000.0, dy = 20000.0 /"//NL// &
            "&physics_nml coriolis_f = 1.0e-4, wind_stress_x = 0.1 /"//NL// &
            "&vcoord_nml vcoord_type = 'sigma' /"//NL// &
            "&ocean_topo_nml topo_config = 'seamount', max_depth = 2000.0, "// &
            "edge_depth = 800.0, slope_scale = 80000.0 /"//NL// &
            "&ocean_bc_nml east = 'open' /"//NL// &
            "&ocean_slopes_nml enable = .true. /"//NL// &
            "&ocean_redi_nml enable = .true., khtr = "//khtr//" /"//NL// &
            "&output_nml output_to_file = .false. /"//NL
   end function channel_nml

   subroutine run_channel(scheme, khtr, salt_res, heat_res, salt1, ok)
      !! Step the channel `N_STEPS` and form the closed salt/heat residuals
      !! `((Q1 - Q0)*rho + out - src) / (Q0*rho)` the console prints.
      character(len=*), intent(in) :: scheme
      character(len=*), intent(in) :: khtr
      real(wp), intent(out) :: salt_res, heat_res, salt1
      logical, intent(out) :: ok
      type(ocean_engine_t) :: engine
      type(config_t) :: cfg
      type(console_stats_t) :: cstats
      type(conservation_budget_t) :: bud
      integer :: ierr, n
      real(wp) :: t, salt0, heat0, heat1

      ok = .false.
      salt_res = huge(1.0_wp)
      heat_res = huge(1.0_wp)
      salt1 = 0.0_wp
      call ensure_comm()
      call read_config_from_string(channel_nml(scheme, khtr), cfg, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call validate_config(cfg, ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call engine_setup(engine, cfg, ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call engine_enter_data(engine, cfg)

      call physical_totals(engine, salt0, heat0)
      t = 0.0_wp
      ierr = OCEAN_STATUS_OK
      do n = 1, N_STEPS
         call engine_step(engine, DT, t, ierr=ierr)
         if (ierr /= OCEAN_STATUS_OK) exit
         call engine_step_finalize(engine, DT, t, ierr=ierr)
         if (ierr /= OCEAN_STATUS_OK) exit
         t = t + DT
      end do
      if (ierr == OCEAN_STATUS_OK) then
         ! The budget accumulators run from t = 0, so one report at the end
         ! carries the whole window.  Its latch of the reference totals is
         ! irrelevant here: the residual is formed below from our own t = 0
         ! snapshot.
         call ocean_console_stats_report(cstats, engine%grid, engine%state%metrics, &
                                         engine%state%multilayer, t, DT, N_STEPS, &
                                         compute_rank=0, reproducing_sums=.true., &
                                         budget_stage_weight=ocean_budget_stage_weight( &
                                         engine%state%dyn%split_scheme == SPLIT_SCHEME_PRED_CORR), &
                                         budget_out=bud)
         call physical_totals(engine, salt1, heat1)
         ok = bud%salt_active .and. bud%heat_active
         salt_res = ((salt1 - salt0)*RHO_WATER + bud%salt_out - bud%salt_src)/(salt0*RHO_WATER)
         heat_res = ((heat1 - heat0)*RHO_WATER + bud%heat_out - bud%heat_src)/(heat0*RHO_WATER)
         ok = ok .and. ieee_is_finite(salt_res) .and. ieee_is_finite(heat_res)
      end if
      call engine_exit_data(engine)
      call engine_teardown(engine)
   end subroutine run_channel

   subroutine physical_totals(engine, salt, heat)
      !! `Σ hTr·areaT` over the physical cells for salinity and temperature
      !! (host, after pulling the two tracer payloads off the device).
      type(ocean_engine_t), intent(inout) :: engine
      real(wp), intent(out) :: salt, heat
      integer :: is, it, ng, nx, ny, nz, i, j, k

      is = engine%state%multilayer%idx_salinity
      it = engine%state%multilayer%idx_temperature
      !$acc update self(engine%state%multilayer%tracers(is)%hTr)
      !$acc update self(engine%state%multilayer%tracers(it)%hTr)
      ng = engine%grid%nghost
      nx = engine%grid%nx_total
      ny = engine%grid%ny_total
      nz = engine%state%multilayer%nz_ml
      salt = 0.0_wp
      heat = 0.0_wp
      do k = 1, nz
         do j = ng + 1, ny - ng
            do i = ng + 1, nx - ng
               salt = salt + engine%state%multilayer%tracers(is)%hTr(i, j, k)* &
                      engine%state%metrics%areaT(i, j)
               heat = heat + engine%state%multilayer%tracers(it)%hTr(i, j, k)* &
                      engine%state%metrics%areaT(i, j)
            end do
         end do
      end do
   end subroutine physical_totals

   subroutine check_scheme(error, scheme)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: scheme
      real(wp) :: salt_res, heat_res, salt_on, salt_off, r0, r1
      logical :: ok
      character(len=96) :: vals

      call run_channel(scheme, "0.0", r0, r1, salt_off, ok)
      call check(error, ok, scheme//": Redi-off control run failed or budget inactive")
      if (allocated(error)) return
      call run_channel(scheme, "500.0", salt_res, heat_res, salt_on, ok)
      call check(error, ok, scheme//": Redi run failed, or the closed budget is inactive")
      if (allocated(error)) return

      ! Positive control: Redi must actually move salt through the open edge.
      call check(error, abs(salt_on - salt_off) > REDI_SIGNAL_MIN*abs(salt_off), &
                 scheme//": Redi left the salt total unchanged -- the test has no signal")
      if (allocated(error)) return

      write (vals, '(a,es10.3,a,es10.3,a,es10.3)') " salt_res =", salt_res, &
         "  heat_res =", heat_res, "  redi_dS =", (salt_on - salt_off)/salt_off
      call check(error, abs(salt_res) <= RESID_TOL, &
                 scheme//": closed SALT residual with Redi + open edge is not at round-off;"// &
                 trim(vals))
      if (allocated(error)) return
      call check(error, abs(heat_res) <= RESID_TOL, &
                 scheme//": closed HEAT residual with Redi + open edge is not at round-off;"// &
                 trim(vals))
   end subroutine check_scheme

   subroutine test_pred_corr(error)
      type(error_type), allocatable, intent(out) :: error
      call check_scheme(error, "pred_corr")
   end subroutine test_pred_corr

   subroutine test_ssp_rk2(error)
      type(error_type), allocatable, intent(out) :: error
      call check_scheme(error, "ssp_rk2")
   end subroutine test_ssp_rk2

end module test_ocean_redi_obc_budget
