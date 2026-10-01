!! Bit-exact warm restart through the PRODUCTION engine path.
!!
!! `test_ocean_restart` builds its two states by hand and never runs
!! `engine_setup`, so it cannot see what the configure chain does to a
!! state AFTER the restart read.  Three defects lived exactly there, and
!! all three made a 1/4-degree Southern Ocean resume diverge from the run
!! that wrote the checkpoint (2026-10-01):
!!   1. the `pred_corr` predictor's carried viscous tendency
!!      (`hvisc%du_visc`/`dv_visc`) was not checkpointed, and the scratch
!!      attach zero-filled it even once it was;
!!   2. `configure_ocean_land_mask` re-seeded every land column to uniform
!!      `H_VANISHED`, while the running model (the ALE remap regrids land
!!      columns like any other) had carried a different layout;
!!   3. the init-time periodic wrap, the host halo exchange and the
!!      device warm-up exchange re-derived the checkpointed ghosts, and at
!!      a step boundary the duplicated seam faces of `u` are not a pure
!!      function of the owned interior.
!!
!! Gate: engine A steps N_WRITE steps, writes a checkpoint, steps
!! N_AFTER more.  Engine B resumes from the checkpoint and steps N_AFTER.
!! Every prognostic array -- FULL local arrays, ghosts included -- must be
!! bitwise identical.  The configuration carries each ingredient: land (an
!! island) under `z_fixed` (so the remap reshapes the land columns), a
!! periodic axis (seam ghosts), and lateral viscosity under the default
!! `pred_corr` (the carried tendency).
module test_ocean_restart_engine
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_config, only: config_t, read_config_from_string, validate_config
   use rdb_ocean_engine, only: ocean_engine_t, engine_setup, engine_enter_data, &
                               engine_step, engine_step_finalize, &
                               engine_exit_data, engine_teardown
   use rdb_ocean_state, only: ocean_state_restart_write
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   use rdb_comm_env, only: comm_env_init, comm_env_setup_roles
   implicit none
   private

   public :: collect_ocean_restart_engine_tests

   integer, parameter :: N_WRITE = 6
   integer, parameter :: N_AFTER = 3
   real(wp), parameter :: DT = 600.0_wp
   character(len=*), parameter :: FN = "test_ocean_restart_engine_rt.nc"

   logical :: comm_inited = .false.

   type :: snapshot_t
      !! Host copies of the prognostic arrays, full local extent.
      real(wp), allocatable :: h(:, :, :), u(:, :, :), v(:, :, :)
      real(wp), allocatable :: s(:, :, :), t(:, :, :)
      real(wp), allocatable :: bt_h(:, :), bt_u(:, :), bt_v(:, :)
   end type snapshot_t

contains

   subroutine collect_ocean_restart_engine_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("restart_engine_bit_exact_island_periodic_zfixed", &
                               test_engine_bit_exact) &
                  ]
   end subroutine collect_ocean_restart_engine_tests

   subroutine ensure_comm()
      !! See `test_ocean_budget_periodic_sponge_serial::ensure_comm`.
      if (.not. comm_inited) then
         call comm_env_init()
         call comm_env_setup_roles(.false.)
         comm_inited = .true.
      end if
   end subroutine ensure_comm

   function case_nml() result(nml)
      !! Periodic channel with an island, `z_fixed` (tanh) layers, wind,
      !! Smagorinsky + a Laplacian floor, KPP (default), `pred_corr`
      !! (default).
      character(len=:), allocatable :: nml
      character(len=*), parameter :: NL = new_line("a")
      nml = "&sim_nml sim_type = 'ocean' /"//NL// &
            "&time_nml t_end = 86400.0, dt_fixed = 600.0 /"//NL// &
            "&nonhydrostatic_nml nz_layers = 6 /"//NL// &
            "&tracer_nml initial_temperature = 12.0, initial_salinity = 35.0, "// &
            "T_init_surface = 20.0, T_init_bottom = 4.0 /"//NL// &
            "&ocean_bt_nml auto_n_inner = .true. /"//NL// &
            "&ocean_hvisc_nml nu_h = 200.0, lateral_closure = 'smagorinsky', "// &
            "smag_ah = .true. /"//NL// &
            "&ocean_diag_nml enabled = .false. /"//NL// &
            "&grid_nml nx = 24, ny = 16, nghost = 3, dx = 10000.0, dy = 10000.0 /"//NL// &
            "&physics_nml coriolis_f = 1.0e-4, wind_stress_x = 0.08 /"//NL// &
            "&vcoord_nml vcoord_type = 'z_fixed', z_fixed_profile = 'tanh', "// &
            "z_fixed_dz_top = 20.0, z_fixed_tanh_center = 0.5, "// &
            "z_fixed_tanh_width = 0.25 /"//NL// &
            "&ocean_topo_nml topo_config = 'island', max_depth = 1000.0, "// &
            "slope_scale = 0.25 /"//NL// &
            "&ocean_bc_nml west = 'periodic', east = 'periodic', south = 'wall', "// &
            "north = 'wall' /"//NL// &
            "&output_nml output_to_file = .false. /"//NL
   end function case_nml

   subroutine make_engine(engine, cfg, ok, restart_file, t0, step0)
      type(ocean_engine_t), intent(inout) :: engine
      type(config_t), intent(inout) :: cfg
      logical, intent(out) :: ok
      character(len=*), intent(in), optional :: restart_file
      real(wp), intent(out), optional :: t0
      integer, intent(out), optional :: step0
      integer :: ierr

      ok = .false.
      call ensure_comm()
      call read_config_from_string(case_nml(), cfg, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call validate_config(cfg, ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      if (present(restart_file)) then
         call engine_setup(engine, cfg, ierr, restart_file=restart_file, &
                           t_restart=t0, step_restart=step0)
      else
         call engine_setup(engine, cfg, ierr)
      end if
      if (ierr /= OCEAN_STATUS_OK) return
      call engine_enter_data(engine, cfg)
      ok = .true.
   end subroutine make_engine

   subroutine advance(engine, t, nsteps, ok)
      type(ocean_engine_t), intent(inout) :: engine
      real(wp), intent(inout) :: t
      integer, intent(in) :: nsteps
      logical, intent(out) :: ok
      integer :: n, ierr
      ok = .false.
      do n = 1, nsteps
         call engine_step(engine, DT, t, ierr=ierr)
         if (ierr /= OCEAN_STATUS_OK) return
         call engine_step_finalize(engine, DT, t, ierr=ierr)
         if (ierr /= OCEAN_STATUS_OK) return
         t = t + DT
      end do
      ok = .true.
   end subroutine advance

   subroutine take_snapshot(engine, snap)
      !! Pull the device copies down COMPONENT-wise (never the aggregate
      !! derived type -- see CLAUDE.md), then copy the full arrays.
      type(ocean_engine_t), intent(inout) :: engine
      type(snapshot_t), intent(out) :: snap
      integer :: is, it
      is = engine%state%multilayer%idx_salinity
      it = engine%state%multilayer%idx_temperature
      !$acc update self(engine%state%multilayer%h_layer)
      !$acc update self(engine%state%multilayer%u_face_x_layer)
      !$acc update self(engine%state%multilayer%v_face_y_layer)
      !$acc update self(engine%state%multilayer%tracers(is)%hTr)
      !$acc update self(engine%state%multilayer%tracers(it)%hTr)
      !$acc update self(engine%state%barotropic%h)
      !$acc update self(engine%state%barotropic%u_face_x)
      !$acc update self(engine%state%barotropic%v_face_y)
      snap%h = engine%state%multilayer%h_layer
      snap%u = engine%state%multilayer%u_face_x_layer
      snap%v = engine%state%multilayer%v_face_y_layer
      snap%s = engine%state%multilayer%tracers(is)%hTr
      snap%t = engine%state%multilayer%tracers(it)%hTr
      snap%bt_h = engine%state%barotropic%h
      snap%bt_u = engine%state%barotropic%u_face_x
      snap%bt_v = engine%state%barotropic%v_face_y
   end subroutine take_snapshot

   subroutine test_engine_bit_exact(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t) :: ea, eb
      type(config_t) :: cfg_a, cfg_b
      type(snapshot_t) :: sa, sb
      real(wp) :: t_a, t_b
      integer :: step_b, ierr
      logical :: ok

      ! ---- A: N_WRITE steps, checkpoint, N_AFTER more ----
      call make_engine(ea, cfg_a, ok)
      call check(error, ok, "engine A setup failed")
      if (allocated(error)) return
      t_a = 0.0_wp
      call advance(ea, t_a, N_WRITE, ok)
      call check(error, ok, "engine A failed before the checkpoint")
      if (allocated(error)) return
      call ocean_state_restart_write(ea%state, ea%grid, ea%decomp, FN, t_a, N_WRITE, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "checkpoint write failed")
      if (allocated(error)) return
      call advance(ea, t_a, N_AFTER, ok)
      call check(error, ok, "engine A failed after the checkpoint")
      if (allocated(error)) return
      call take_snapshot(ea, sa)
      call engine_exit_data(ea)
      call engine_teardown(ea)

      ! ---- B: resume from the checkpoint, N_AFTER steps ----
      call make_engine(eb, cfg_b, ok, restart_file=FN, t0=t_b, step0=step_b)
      call check(error, ok, "engine B (warm restart) setup failed")
      if (allocated(error)) return
      call check(error, step_b == N_WRITE, "restart step count not restored")
      if (allocated(error)) return
      call advance(eb, t_b, N_AFTER, ok)
      call check(error, ok, "engine B failed after the restart")
      if (allocated(error)) return
      call take_snapshot(eb, sb)
      call engine_exit_data(eb)
      call engine_teardown(eb)
      call delete_file(FN)

      call check(error, all(sa%h == sb%h), "h_layer differs after a warm restart")
      if (allocated(error)) return
      call check(error, all(sa%u == sb%u), "u_face_x_layer differs after a warm restart")
      if (allocated(error)) return
      call check(error, all(sa%v == sb%v), "v_face_y_layer differs after a warm restart")
      if (allocated(error)) return
      call check(error, all(sa%s == sb%s), "salinity hTr differs after a warm restart")
      if (allocated(error)) return
      call check(error, all(sa%t == sb%t), "temperature hTr differs after a warm restart")
      if (allocated(error)) return
      call check(error, all(sa%bt_h == sb%bt_h) .and. all(sa%bt_u == sb%bt_u) .and. &
                 all(sa%bt_v == sb%bt_v), "barotropic state differs after a warm restart")
   end subroutine test_engine_bit_exact

   subroutine delete_file(fname)
      character(len=*), intent(in) :: fname
      integer :: u, ios
      logical :: exists
      inquire (file=fname, exist=exists)
      if (.not. exists) return
      open (newunit=u, file=fname, status="old", iostat=ios)
      if (ios == 0) close (u, status="delete")
   end subroutine delete_file

end module test_ocean_restart_engine
