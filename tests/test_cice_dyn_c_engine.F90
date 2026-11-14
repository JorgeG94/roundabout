!! The CICE6 C-grid EVP dynamics under the SIS2-derived ice state, through the
!! PRODUCTION engine path (D9.2): the switch `&ocean_ice_dyn_nml dyn_scheme`,
!! the slot, its coupling to the ocean (tilt from the surface height, ocean
!! surface current, stress hand-off), restart and the configure-time refusals.
!!
!! Every case builds an engine from a namelist (`engine_setup` ->
!! `engine_enter_data`) and steps it the driver's way (`engine_step` ->
!! `engine_step_ice` -> `engine_step_finalize`); the analytic cases call
!! `ice_dyn_c_step` on the engine's own state with the ocean held frozen, which
!! is how a closed form can be asked of a coupled step.
!!
!! Cases:
!!   * `default_slot_absent`    `dyn_scheme = "sis2_c"` (the default): the C-grid slot
!!                              is not initialised, counts 0 bytes, registers no
!!                              restart tag and the SIS2 workspace is there.
!!   * `cice_c_slot_present`    `"cice_c"`: slot, metrics and restart tags exist,
!!                              the SIS2 EVP workspace does not.
!!   * `refusals`               every refusal of `ice_config_check` names its cause.
!!   * `free_drift_tilt`        uniform ice, no wind, no Coriolis, the rheology off, a
!!                              sea-surface slope built into the layer thicknesses:
!!                              the coupled step converges to `rho_w C_w |u| u =
!!                              -m g grad(eta)` (closed form, both components).
!!   * `new_ice_ocean_start`    faces that become ice-covered start at the ocean
!!                              current and stay on it with no force (u and v).
!!   * `closed_basin_wind`      walls and wind: no ice velocity through a wall, the
!!                              ocean stress is `(1-ai) tau_a - strocn` face by face and
!!                              `strocnT` is the flux-conserving map of the face stress.
!!   * `velocity_guard`         a NaN velocity and a speed above `vel_cap` raise
!!                              `ICE_ERR_DYN_VEL` (nothing is clamped).
!!   * `restart_round_trip`     straight N+M steps against N + checkpoint + warm
!!                              restart + M steps, every registry field bitwise,
!!                              ghosts included (periodic channel, thermodynamics,
!!                              transport on the new velocities).
!!   * `vs_sis2_free_drift` / `vs_sis2_basin`  the two schemes on a uniform pack and in
!!                              a closed basin (numbers are printed: they are
!!                              different discretisations).
module test_cice_dyn_c_engine
   use, intrinsic :: iso_fortran_env, only: int32, int64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_config, only: config_t, read_config_from_string
   use rdb_ocean_engine, only: ocean_engine_t, engine_setup, engine_enter_data, &
                               engine_step, engine_step_ice, engine_step_finalize, &
                               engine_exit_data, engine_teardown
   use rdb_ocean_state, only: ocean_state_restart_write, ocean_state_build_restart_registry
   use rdb_ocean_restart, only: restart_registry_t
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   use rdb_comm_env, only: comm_env_init, comm_env_setup_roles
   use rdb_ice_config, only: ice_config_check, ICE_CONFIG_MSG_LEN, ICE_CONFIG_MAX_MSGS
   use rdb_ice_dyn_c, only: ice_dyn_c_step
   use rdb_ice_evp, only: ice_evp_step
   use rdb_ice_errors, only: ice_err_summary_t, ICE_ERR_DYN_VEL
   use rdb_ice_column, only: ICE_RHO_ICE
   implicit none
   private

   public :: collect_cice_dyn_c_engine_tests

   real(wp), parameter :: DT = 600.0_wp
   character(len=*), parameter :: FN = "test_cice_dyn_c_engine_rt.nc"
   integer, parameter :: MAXF = 200
   logical :: comm_inited = .false.
   integer, parameter :: MAXEX = 8
   integer :: n_ex = 0
   character(len=40) :: ex_group(MAXEX) = "", ex_key(MAXEX) = "", ex_val(MAXEX) = ""
      !! Experiment overrides applied to every `case_nml`: (group, key, value).

   type :: field_t
      character(len=64) :: tag = ""
      real(wp), allocatable :: a(:, :, :)
   end type field_t

   type :: snapshot_t
      integer :: n = 0
      type(field_t) :: f(MAXF)
   end type snapshot_t

contains

   subroutine collect_cice_dyn_c_engine_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("default_slot_absent", test_default_absent), &
                  new_unittest("cice_c_slot_present", test_slot_present), &
                  new_unittest("refusals", test_refusals), &
                  new_unittest("free_drift_tilt", test_free_drift_tilt), &
                  new_unittest("new_ice_ocean_start", test_ocean_start), &
                  new_unittest("closed_basin_wind", test_closed_basin), &
                  new_unittest("velocity_guard", test_guard), &
                  new_unittest("restart_round_trip", test_restart), &
                  new_unittest("vs_sis2_free_drift", test_vs_sis2_free_drift), &
                  new_unittest("vs_sis2_basin", test_vs_sis2_basin), &
                  new_unittest("exp1_basin_no_strength", test_exp1), &
                  new_unittest("exp2_uniform_steady", test_exp2), &
                  new_unittest("exp3_basin_viscosity_floor", test_exp3), &
                  new_unittest("restart_refused_across_schemes", test_restart_cross), &
                  new_unittest("step_throughput", test_throughput) &
                  ]
   end subroutine collect_cice_dyn_c_engine_tests

   ! ------------------------------------------------------------------
   ! namelists and engine plumbing
   ! ------------------------------------------------------------------

   function case_nml(label) result(nml)
      !! `label`: `<geometry>[_cice]`, geometry one of `basin` (closed, flat,
      !! 100 % ice, no thermodynamics, wind), `free` (the same, no wind, no
      !! Coriolis, the rheology dropped: pure free drift), `uniform` (doubly
      !! periodic uniform pack, wind), `channel` (periodic x, thermodynamics,
      !! transport: the restart case).
      character(len=*), intent(in) :: label
      character(len=:), allocatable :: nml
      character(len=*), parameter :: NL = new_line("a")
      character(len=:), allocatable :: geom, ice, bc, topo, phys, thermo, dynnml
      logical :: cice
      integer :: k

      cice = index(label, "_cice") > 0
      geom = label
      if (cice) geom = label(1:index(label, "_cice") - 1)
      dynnml = ""
      if (cice) dynnml = "&ocean_ice_dyn_nml dyn_scheme = 'cice_c', ndte = 20 /"//NL// &
                         "&ocean_ice_ridge_nml kstrength = 0 /"//NL
      ice = "&ocean_ice_nml enable = .true., ncat = 5, nk_ice = 2, dynamics = .true., "// &
            "evp_sub_steps = 20, a_face_stress = .true., cdw = 5.36e-3, rho_ocean = 1026.0 /"//NL// &
            "&ocean_ice_ic_nml conc_config = 'uniform', h_ice = 1.0, conc = 1.0 /"//NL
      bc = "&ocean_bc_nml west = 'wall', east = 'wall', south = 'wall', north = 'wall' /"//NL
      phys = "&physics_nml coriolis_f = 1.0e-4, wind_stress_x = 0.1, wind_stress_y = 0.03 /"//NL
      thermo = ""
      topo = "&ocean_topo_nml topo_config = 'flat', max_depth = 200.0 /"//NL
      select case (geom)
      case ("basin")
         continue
      case ("free")
         phys = "&physics_nml coriolis_f = 0.0 /"//NL
         dynnml = "&ocean_ice_dyn_nml rheo_area_min = 1.1 /"//NL
         if (cice) dynnml = "&ocean_ice_dyn_nml dyn_scheme = 'cice_c', ndte = 20, rheo_area_min = 1.1 /"//NL// &
                            "&ocean_ice_ridge_nml kstrength = 0 /"//NL
      case ("uniform")
         bc = "&ocean_bc_nml west = 'periodic', east = 'periodic', south = 'periodic', "// &
              "north = 'periodic' /"//NL
         phys = "&physics_nml coriolis_f = 1.4e-4, wind_stress_x = 0.1, wind_stress_y = 0.0 /"//NL
      case ("channel")
         bc = "&ocean_bc_nml west = 'periodic', east = 'periodic', south = 'wall', north = 'wall' /"//NL
         phys = "&physics_nml coriolis_f = 1.4e-4, wind_stress_x = 0.1, wind_stress_y = 0.04 /"//NL
         topo = "&ocean_topo_nml topo_config = 'seamount', max_depth = 1000.0, edge_depth = 800.0, "// &
                "slope_scale = 60000.0 /"//NL
         thermo = "&ocean_thermo_nml enable_thermodynamics = .true., q_heat = -100.0 /"//NL// &
                  "&ocean_vmix_nml dt_therm_ratio = 2 /"//NL
         ice = "&ocean_ice_nml enable = .true., ncat = 5, nk_ice = 2, air_temp = -20.0, "// &
               "restore_lambda = 20.0, sw_down = 40.0, snowfall = 2.0e-5, dynamics = .true., "// &
               "evp_sub_steps = 20, transport = .true., adv_substeps = 2 /"//NL// &
               "&ocean_ice_ic_nml conc_config = 'uniform', h_ice = 1.0, conc = 0.7 /"//NL
      case default
         error stop "test_cice_dyn_c_engine: unknown case"
      end select
      nml = "&sim_nml sim_type = 'ocean' /"//NL// &
            "&time_nml t_end = 86400.0, dt_fixed = 600.0 /"//NL// &
            "&nonhydrostatic_nml nz_layers = 3 /"//NL// &
            "&tracer_nml initial_temperature = -1.83, initial_salinity = 34.0, "// &
            "T_init_surface = -1.82, T_init_bottom = -1.835 /"//NL// &
            "&ocean_bt_nml auto_n_inner = .true. /"//NL// &
            "&ocean_hvisc_nml nu_h = 200.0 /"//NL// &
            "&ocean_diag_nml enabled = .false. /"//NL// &
            "&grid_nml nx = 24, ny = 16, nghost = 3, dx = 20000.0, dy = 20000.0 /"//NL// &
            phys// &
            "&vcoord_nml vcoord_type = 'sigma' /"//NL// &
            topo//thermo//ice//dynnml//bc// &
            "&output_nml output_to_file = .false. /"//NL
      do k = 1, n_ex
         nml = set_nml(nml, trim(ex_group(k)), trim(ex_key(k)), trim(ex_val(k)))
      end do
   end function case_nml

   function set_nml(nml, group, key, val) result(r)
      !! `key = val` in the namelist group `group` of `nml`: replaced when present, added
      !! otherwise (a group that is absent is appended).
      character(len=*), intent(in) :: nml, group, key, val
      character(len=:), allocatable :: r
      integer :: g, e, k, ve
      character(len=*), parameter :: NL = new_line("a")

      g = index(nml, "&"//group//" ")
      if (g == 0) then
         r = nml//"&"//group//" "//key//" = "//val//" /"//NL
         return
      end if
      e = g + index(nml(g:), " /") - 1          ! the closing " /" of this group
      k = index(nml(g:e), " "//key//" =")
      if (k == 0) then
         r = nml(1:e - 1)//", "//key//" = "//val//nml(e:)
      else
         k = g + k - 1 + len(" "//key//" = ")    ! first character of the old value
         ve = k + scan(nml(k:e), ", /") - 2
         r = nml(1:k - 1)//val//nml(ve + 1:)
      end if
   end function set_nml

   subroutine ex_clear()
      n_ex = 0
   end subroutine ex_clear

   subroutine ex_add(group, key, val)
      character(len=*), intent(in) :: group, key, val
      n_ex = n_ex + 1
      ex_group(n_ex) = group
      ex_key(n_ex) = key
      ex_val(n_ex) = val
   end subroutine ex_add

   subroutine ensure_comm()
      if (.not. comm_inited) then
         call comm_env_init()
         call comm_env_setup_roles(.false.)
         comm_inited = .true.
      end if
   end subroutine ensure_comm

   subroutine make_engine(engine, cfg, nml, ok, restart_file, t0, step0)
      type(ocean_engine_t), intent(inout) :: engine
      type(config_t), intent(inout) :: cfg
      character(len=*), intent(in) :: nml
      logical, intent(out) :: ok
      character(len=*), intent(in), optional :: restart_file
      real(wp), intent(out), optional :: t0
      integer, intent(out), optional :: step0
      integer :: ierr

      ok = .false.
      call ensure_comm()
      call read_config_from_string(nml, cfg, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      if (present(restart_file)) then
         call engine_setup(engine, cfg, ierr, restart_file=restart_file, t_restart=t0, step_restart=step0)
      else
         call engine_setup(engine, cfg, ierr)
      end if
      if (ierr /= OCEAN_STATUS_OK) return
      call engine_enter_data(engine, cfg)
      ok = .true.
   end subroutine make_engine

   subroutine advance(engine, cfg, t, nsteps, ok)
      type(ocean_engine_t), intent(inout) :: engine
      type(config_t), intent(inout) :: cfg
      real(wp), intent(inout) :: t
      integer, intent(in) :: nsteps
      logical, intent(out) :: ok
      integer :: n, ierr
      ok = .false.
      do n = 1, nsteps
         call engine_step(engine, DT, t, ierr=ierr)
         if (ierr /= OCEAN_STATUS_OK) return
         call engine_step_ice(engine, cfg, DT, t, ierr=ierr)
         if (ierr /= OCEAN_STATUS_OK) return
         call engine_step_finalize(engine, DT, t, ierr=ierr)
         if (ierr /= OCEAN_STATUS_OK) return
         t = t + DT
      end do
      ok = .true.
   end subroutine advance

   subroutine finish_engine(engine)
      type(ocean_engine_t), intent(inout) :: engine
      call engine_exit_data(engine)
      call engine_teardown(engine)
   end subroutine finish_engine

   subroutine pull_ice(engine)
      !! Device -> host of the component arrays the cases read.
      type(ocean_engine_t), intent(inout) :: engine
      !$acc update self(engine%state%ice%u_ice, engine%state%ice%v_ice)
      !$acc update self(engine%state%ice%tau_a_x, engine%state%ice%tau_a_y)
      !$acc update self(engine%state%surface_stress%tau_x, engine%state%surface_stress%tau_y)
   end subroutine pull_ice

   subroutine pull_dync(engine)
      type(ocean_engine_t), intent(inout) :: engine
      !$acc update self(engine%state%ice%dync%w%strocnxE, engine%state%ice%dync%w%strocnyN)
      !$acc update self(engine%state%ice%dync%w%aiE, engine%state%ice%dync%w%aiN)
      !$acc update self(engine%state%ice%dync%strocnT_x, engine%state%ice%dync%strocnT_y)
      !$acc update self(engine%state%ice%dync%wE, engine%state%ice%dync%wN)
   end subroutine pull_dync

   ! ------------------------------------------------------------------
   ! tests: the switch
   ! ------------------------------------------------------------------

   subroutine test_default_absent(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      type(restart_registry_t) :: reg
      logical :: ok
      integer :: e, n_dync

      call make_engine(engine, cfg, case_nml("basin"), ok)
      call check(error, ok, "engine setup failed")
      if (allocated(error)) return
      call check(error,.not. engine%state%ice%dync%is_init, "the C-grid slot is not initialised by default")
      if (allocated(error)) return
      call check(error, engine%state%ice%dync%bytes() == 0_int64, "the absent slot counts 0 bytes")
      if (allocated(error)) return
      call check(error,.not. engine%state%ice_cgrid%is_init, "no C-grid metric set by default")
      if (allocated(error)) return
      call check(error, engine%state%ice%evp_ws%bytes() > 0_int64, "the SIS2 EVP workspace is there")
      if (allocated(error)) return
      call ocean_state_build_restart_registry(engine%state, engine%grid, reg)
      n_dync = 0
      do e = 1, reg%n
         if (index(reg%entries(e)%tag, "ice_dync") == 1) n_dync = n_dync + 1
      end do
      call check(error, n_dync == 0, "no C-grid restart tag by default")
      call finish_engine(engine)
   end subroutine test_default_absent

   subroutine test_slot_present(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      type(restart_registry_t) :: reg
      logical :: ok
      integer :: e, n_dync

      call make_engine(engine, cfg, case_nml("basin_cice"), ok)
      call check(error, ok, "engine setup failed")
      if (allocated(error)) return
      call check(error, engine%state%ice%dync%is_init, "the C-grid slot is initialised")
      if (allocated(error)) return
      call check(error, engine%state%ice%dync%bytes() > 0_int64, "the slot counts its arrays")
      if (allocated(error)) return
      call check(error, engine%state%ice%bytes() > engine%state%ice%dync%bytes(), &
                 "the ice slot total includes the C-grid slot")
      if (allocated(error)) return
      call check(error, engine%state%ice_cgrid%is_init .and. engine%state%ice_cgrid%mapped, &
                 "the C-grid metric set is built and mapped")
      if (allocated(error)) return
      call check(error, engine%state%ice%evp_ws%bytes() == 0_int64, "the SIS2 EVP workspace is not allocated")
      if (allocated(error)) return
      call ocean_state_build_restart_registry(engine%state, engine%grid, reg)
      n_dync = 0
      do e = 1, reg%n
         if (index(reg%entries(e)%tag, "ice_dync") == 1) n_dync = n_dync + 1
      end do
      call check(error, n_dync == 7, "seven C-grid restart tags (3 stresses, 2 masks, strocnT x/y)")
      call finish_engine(engine)
   end subroutine test_slot_present

   subroutine refusal_text(nml, text)
      !! All refusal messages of `nml` joined into one string.
      character(len=*), intent(in) :: nml
      character(len=:), allocatable, intent(out) :: text
      type(config_t) :: cfg
      character(len=ICE_CONFIG_MSG_LEN) :: msgs(ICE_CONFIG_MAX_MSGS)
      integer :: n, i, ierr

      text = ""
      call read_config_from_string(nml, cfg, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) then
         text = "PARSE-REFUSED"
         return
      end if
      call ice_config_check(cfg, msgs, n)
      do i = 1, n
         text = text//" | "//trim(msgs(i))
      end do
   end subroutine refusal_text

   subroutine test_refusals(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: t, base
      character(len=*), parameter :: NL = new_line("a")

      call ensure_comm()
      base = case_nml("basin_cice")
      call refusal_text(base, t)
      call check(error, len_trim(t) == 0, "a valid cice_c configuration is accepted: "//t)
      if (allocated(error)) return
      call refusal_text(case_nml("basin"), t)
      call check(error, len_trim(t) == 0, "the default configuration is accepted: "//t)
      if (allocated(error)) return
      call refusal_text(replace(base, "dyn_scheme = 'cice_c'", "dyn_scheme = 'nope'"), t)
      call check(error, index(t, "PARSE-REFUSED") > 0 .or. index(t, "dyn_scheme") > 0, "an unknown scheme is refused")
      if (allocated(error)) return
      call refusal_text(replace(base, "dynamics = .true.", "dynamics = .false."), t)
      call check(error, index(t, "dynamics=.true.") > 0, "cice_c needs dynamics=.true.: "//t)
      if (allocated(error)) return
      call refusal_text(replace(base, "nghost = 3", "nghost = 2"), t)
      call check(error, index(t, "nghost >= 3") > 0, "cice_c needs nghost >= 3: "//t)
      if (allocated(error)) return
      call refusal_text(base//"&mpi_nml px = 2, py = 1 /"//NL, t)
      call check(error, index(t, "single-rank") > 0, "cice_c is single-rank: "//t)
      if (allocated(error)) return
      call refusal_text(replace(base, "north = 'wall'", "north = 'tripolar_fold'"), t)
      call check(error, index(t, "tripolar") > 0, "cice_c refuses the tripolar fold: "//t)
   end subroutine test_refusals

   ! ------------------------------------------------------------------
   ! tests: the coupled step against closed forms (ocean frozen)
   ! ------------------------------------------------------------------

   subroutine test_free_drift_tilt(error)
      !! No wind, no Coriolis, rheology dropped, a uniform sea-surface slope
      !! `grad(eta) = (sx, sy)`: the free-drift balance of the drag against the
      !! tilt force, `rho_w C_w |u| u = -m g grad(eta)`, is the fixed point of the
      !! coupled step.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      logical :: ok
      integer :: i, j, n, nx, ny, nz, ng, is, ie, js, je, nbad
      real(wp), parameter :: SX = 1.0e-5_wp, SY = -0.5e-5_wp
      real(wp) :: dx, dy, g, rhow, cw, m, gmag, vmag, ue, ve, worst, mass_ice
      real(wp), allocatable :: dh(:, :)
      type(ice_err_summary_t) :: summ

      call make_engine(engine, cfg, case_nml("free_cice"), ok)
      call check(error, ok, "engine setup failed")
      if (allocated(error)) return
      associate (st => engine%state, gr => engine%grid)
         nx = gr%nx_total
         ny = gr%ny_total
         nz = st%multilayer%nz_ml
         ng = gr%nghost
         is = ng + 1
         ie = ng + gr%nx_phys
         js = ng + 1
         je = ng + gr%ny_phys
         dx = 20000.0_wp
         dy = 20000.0_wp
         allocate (dh(nx, ny))
         do j = 1, ny
            do i = 1, nx
               dh(i, j) = SX*dx*real(i, wp) + SY*dy*real(j, wp)
            end do
         end do
         ! the slope goes into the top layer thickness (ssh = sum h - b)
         !$acc update self(st%multilayer%h_layer)
         do j = 1, ny
            do i = 1, nx
               st%multilayer%h_layer(i, j, nz) = st%multilayer%h_layer(i, j, nz) + dh(i, j)
            end do
         end do
         !$acc update device(st%multilayer%h_layer)
         g = st%ice%dync%p%c%gravit
         rhow = st%ice%dync%p%c%rhow
         cw = st%ice%dync%p%dyn%cdn_ocn
         mass_ice = ICE_RHO_ICE*1.0_wp
         m = mass_ice
         gmag = sqrt(SX*SX + SY*SY)
         vmag = sqrt(m*g*gmag/(rhow*cw))
         ue = -vmag*SX/gmag
         ve = -vmag*SY/gmag
         do n = 1, 150
            call ice_dyn_c_step(gr, st%metrics, st%coriolis_adv%f_corner, st%barotropic%b, st%ice_cgrid, &
                                st%ice, st%multilayer, st%bc, DT, err_summary=summ)
            call check(error, summ%mask == 0_int32, "free-drift step raised an error bit")
            if (allocated(error)) return
         end do
         call pull_ice(engine)
         worst = 0.0_wp
         nbad = 0
         ! Away from the walls: the no-slip coast perturbs the tangential velocity of the
         ! wall-adjacent faces, and that perturbation decays over the next few cells.
         do j = js + 4, je - 4
            do i = is + 4, ie - 4
               ! interior faces: x-face i+1 of E(i,j), y-face j+1 of N(i,j)
               worst = max(worst, abs(st%ice%u_ice(i + 1, j) - ue)/vmag)
               worst = max(worst, abs(st%ice%v_ice(i, j + 1) - ve)/vmag)
            end do
         end do
         if (.not. ieee_is_finite(worst)) nbad = 1
         write (*, '(a,es10.3,a,2es12.4)') "free_drift_tilt: worst |u-u_exact|/|u_exact| = ", worst, &
            "  (u_exact, v_exact) =", ue, ve
         call check(error, nbad == 0, "non-finite velocity")
         if (allocated(error)) return
         call check(error, worst < 1.5e-6_wp, "free drift on a tilted surface is the closed form")
      end associate
      call finish_engine(engine)
      if (allocated(dh)) deallocate (dh)
   end subroutine test_free_drift_tilt

   subroutine test_ocean_start(error)
      !! The ocean surface current is uniform (0.2, 0.1) m/s, no force at all: a
      !! face that becomes ice-covered starts at that current (along-face
      !! component), and with the tangential component interpolated from the
      !! four neighbours the relative velocity is zero, so it stays there.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      logical :: ok
      integer :: i, j, nx, ny, nz, ng, is, ie, js, je
      real(wp), parameter :: UO = 0.2_wp, VO = 0.1_wp
      real(wp) :: worst
      type(ice_err_summary_t) :: summ

      call make_engine(engine, cfg, case_nml("free_cice"), ok)
      call check(error, ok, "engine setup failed")
      if (allocated(error)) return
      associate (st => engine%state, gr => engine%grid)
         nx = gr%nx_total
         ny = gr%ny_total
         nz = st%multilayer%nz_ml
         ng = gr%nghost
         is = ng + 1
         ie = ng + gr%nx_phys
         js = ng + 1
         je = ng + gr%ny_phys
         st%multilayer%u_face_x_layer = UO
         st%multilayer%v_face_y_layer = VO
         !$acc update device(st%multilayer%u_face_x_layer, st%multilayer%v_face_y_layer)
         call ice_dyn_c_step(gr, st%metrics, st%coriolis_adv%f_corner, st%barotropic%b, st%ice_cgrid, &
                             st%ice, st%multilayer, st%bc, DT, err_summary=summ)
         call check(error, summ%mask == 0_int32, "step raised an error bit")
         if (allocated(error)) return
         call pull_ice(engine)
         worst = 0.0_wp
         do j = js + 1, je - 1
            do i = is + 1, ie - 1
               worst = max(worst, abs(st%ice%u_ice(i + 1, j) - UO), abs(st%ice%v_ice(i, j + 1) - VO))
            end do
         end do
         write (*, '(a,es10.3)') "new_ice_ocean_start: worst |u_ice - u_ocean| = ", worst
         call check(error, worst < 1.0e-12_wp, "newly covered faces start at, and stay on, the ocean current")
      end associate
      call finish_engine(engine)
   end subroutine test_ocean_start

   ! ------------------------------------------------------------------
   ! tests: closed basin, guard
   ! ------------------------------------------------------------------

   subroutine test_closed_basin(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      logical :: ok
      real(wp) :: t, worst, scale, umax, usum, tsum, fsum
      integer :: i, j, nx, ny, ng, is, ie, js, je, nwall
      real(wp) :: a_u, resid

      call make_engine(engine, cfg, case_nml("basin_cice"), ok)
      call check(error, ok, "engine setup failed")
      if (allocated(error)) return
      t = 0.0_wp
      call advance(engine, cfg, t, 12, ok)
      call check(error, ok, "12 coupled steps failed")
      if (allocated(error)) return
      call pull_ice(engine)
      call pull_dync(engine)
      associate (st => engine%state, gr => engine%grid)
         nx = gr%nx_total
         ny = gr%ny_total
         ng = gr%nghost
         is = ng + 1
         ie = ng + gr%nx_phys
         js = ng + 1
         je = ng + gr%ny_phys
         ! ---- no velocity through a wall ----
         nwall = 0
         do j = 1, ny
            if (st%ice%u_ice(is, j) /= 0.0_wp) nwall = nwall + 1
            if (st%ice%u_ice(ie + 1, j) /= 0.0_wp) nwall = nwall + 1
         end do
         do i = 1, nx
            if (st%ice%v_ice(i, js) /= 0.0_wp) nwall = nwall + 1
            if (st%ice%v_ice(i, je + 1) /= 0.0_wp) nwall = nwall + 1
         end do
         call check(error, nwall == 0, "no ice velocity through a wall face")
         if (allocated(error)) return
         umax = 0.0_wp
         usum = 0.0_wp
         do j = js, je
            do i = is + 1, ie
               if (.not. ieee_is_finite(st%ice%u_ice(i, j))) then
                  call check(error, .false., "non-finite ice velocity")
                  return
               end if
               umax = max(umax, abs(st%ice%u_ice(i, j)))
               usum = usum + st%ice%u_ice(i, j)
            end do
         end do
         call check(error, umax > 1.0e-6_wp, "the wind moves the ice")
         if (allocated(error)) return
         call check(error, usum > 0.0_wp, "the ice drifts with the (eastward) wind")
         if (allocated(error)) return
         ! ---- the ocean receives equal and opposite stress, face by face ----
         worst = 0.0_wp
         scale = 0.0_wp
         do j = js, je
            do i = is + 1, ie
               ! x-face i is E(i-1,j)
               a_u = st%ice%dync%w%aiE(i - 1, j)
               resid = st%surface_stress%tau_x(i, j) - (1.0_wp - a_u)*st%ice%tau_a_x(i, j) &
                       + st%ice%dync%w%strocnxE(i - 1, j)
               worst = max(worst, abs(resid))
               scale = max(scale, abs(st%ice%dync%w%strocnxE(i - 1, j)))
            end do
         end do
         write (*, '(a,es10.3,a,es10.3)') "closed_basin_wind: ocean stress residual ", worst, &
            "  max |strocnxE| ", scale
         call check(error, scale > 1.0e-6_wp, "the ice-ocean stress is not vacuous")
         if (allocated(error)) return
         call check(error, worst <= 1.0e-13_wp*max(scale, 1.0e-3_wp), "ocean stress = (1-ai) tau_a - strocn")
         if (allocated(error)) return
         ! ---- strocnT: flux-conserving map of the face stress per ice area ----
         tsum = 0.0_wp
         fsum = 0.0_wp
         do j = js, je
            do i = is, ie
               tsum = tsum + st%ice_cgrid%tarea(i, j)*st%ice%dync%strocnT_x(i, j)
               fsum = fsum + st%ice_cgrid%earea(i, j)*st%ice%dync%wE(i, j)
            end do
         end do
         call check(error, abs(tsum - fsum) <= 1.0e-10_wp*max(abs(fsum), 1.0e-30_wp) + 1.0e-12_wp*abs(tsum), &
                    "sum(A_T strocnT_x) = sum(A_E strocn/ai)")
      end associate
      call finish_engine(engine)
   end subroutine test_closed_basin

   subroutine test_guard(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      logical :: ok
      integer :: ng
      type(ice_err_summary_t) :: summ

      call make_engine(engine, cfg, case_nml("basin_cice"), ok)
      call check(error, ok, "engine setup failed")
      if (allocated(error)) return
      associate (st => engine%state, gr => engine%grid)
         ng = gr%nghost
         call ice_dyn_c_step(gr, st%metrics, st%coriolis_adv%f_corner, st%barotropic%b, st%ice_cgrid, &
                             st%ice, st%multilayer, st%bc, DT, err_summary=summ)
         call check(error, iand(summ%mask, ICE_ERR_DYN_VEL) == 0_int32, "a clean step raises no velocity bit")
         if (allocated(error)) return
         ! a NaN on an interior ice face
         st%ice%u_ice(ng + 8, ng + 5) = ieee_value(1.0_wp, ieee_quiet_nan)
         !$acc update device(st%ice%u_ice)
         call ice_dyn_c_step(gr, st%metrics, st%coriolis_adv%f_corner, st%barotropic%b, st%ice_cgrid, &
                             st%ice, st%multilayer, st%bc, DT, err_summary=summ)
         call check(error, iand(summ%mask, ICE_ERR_DYN_VEL) /= 0_int32, "a NaN velocity raises ICE_ERR_DYN_VEL")
         if (allocated(error)) return
         call check(error, summ%first_idx > 0_int64, "the first offender is reported")
      end associate
      call finish_engine(engine)
   end subroutine test_guard

   ! ------------------------------------------------------------------
   ! restart
   ! ------------------------------------------------------------------

   subroutine take_snapshot(engine, snap)
      type(ocean_engine_t), intent(inout), target :: engine
      type(snapshot_t), intent(inout) :: snap
      type(restart_registry_t) :: reg
      integer :: e

      call ocean_state_build_restart_registry(engine%state, engine%grid, reg)
      snap%n = 0
      do e = 1, reg%n
         associate (en => reg%entries(e))
            if (en%device_mapped .and. en%rank == 2) then
               !$acc update self(en%p2)
            else if (en%device_mapped .and. en%rank == 3) then
               !$acc update self(en%p3)
            end if
            select case (en%rank)
            case (0)
               call push(snap, en%tag, reshape([en%p0], [1, 1, 1]))
            case (2)
               call push(snap, en%tag, reshape(en%p2, [size(en%p2, 1), size(en%p2, 2), 1]))
            case default
               call push(snap, en%tag, en%p3)
            end select
         end associate
      end do
   end subroutine take_snapshot

   subroutine push(snap, tag, a)
      type(snapshot_t), intent(inout) :: snap
      character(len=*), intent(in) :: tag
      real(wp), intent(in) :: a(:, :, :)
      if (snap%n >= MAXF) error stop "test_cice_dyn_c_engine: raise MAXF"
      snap%n = snap%n + 1
      snap%f(snap%n)%tag = tag
      snap%f(snap%n)%a = a
   end subroutine push

   subroutine compare(error, sa, sb, n_dync, when)
      type(error_type), allocatable, intent(out) :: error
      type(snapshot_t), intent(in) :: sa, sb
      character(len=*), intent(in) :: when
      integer, intent(out) :: n_dync
      integer :: e, i, j, k, nbad_fields, nb
      character(len=:), allocatable :: msg

      n_dync = 0
      nbad_fields = 0
      msg = ""
      call check(error, sa%n == sb%n, "registry field count differs "//when)
      if (allocated(error)) return
      do e = 1, sa%n
         if (sa%f(e)%tag /= sb%f(e)%tag) then
            call check(error, .false., "registry field order differs: "//trim(sa%f(e)%tag))
            return
         end if
         if (index(sa%f(e)%tag, "ice_dync") == 1) n_dync = n_dync + 1
         if (any(shape(sa%f(e)%a) /= shape(sb%f(e)%a))) then
            nbad_fields = nbad_fields + 1
            msg = msg//" "//trim(sa%f(e)%tag)//"(shape)"
            cycle
         end if
         nb = 0
         do k = 1, size(sa%f(e)%a, 3)
            do j = 1, size(sa%f(e)%a, 2)
               do i = 1, size(sa%f(e)%a, 1)
                  if (transfer(sa%f(e)%a(i, j, k), 0_int64) /= transfer(sb%f(e)%a(i, j, k), 0_int64)) then
                     if (nb == 0) write (*, '(a,a,3(1x,i0),a,es24.16,a,es24.16)') when//": ", &
                        trim(sa%f(e)%tag)//" first at", i, j, k, ": ", sa%f(e)%a(i, j, k), " vs ", &
                        sb%f(e)%a(i, j, k)
                     nb = nb + 1
                  end if
               end do
            end do
         end do
         if (nb > 0) then
            nbad_fields = nbad_fields + 1
            msg = msg//" "//trim(sa%f(e)%tag)
         end if
      end do
      call check(error, nbad_fields == 0, "fields differ "//when//":"//msg)
   end subroutine compare

   subroutine test_restart(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: ea, eb
      type(config_t) :: cfg_a, cfg_b
      type(snapshot_t) :: sa, sb, sa0, sb0
      real(wp) :: t_a, t_b
      integer :: step_b, ierr, n_dync, nw, na
      logical :: ok
      real(wp) :: umax

      nw = 5          ! odd: the checkpoint falls inside a thermodynamics window (ratio 2)
      na = 4
      call make_engine(ea, cfg_a, case_nml("channel_cice"), ok)
      call check(error, ok, "engine A setup failed")
      if (allocated(error)) return
      t_a = 0.0_wp
      call advance(ea, cfg_a, t_a, nw, ok)
      call check(error, ok, "engine A failed before the checkpoint")
      if (allocated(error)) return
      call ocean_state_restart_write(ea%state, ea%grid, ea%decomp, FN, t_a, nw, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "checkpoint write failed")
      if (allocated(error)) return
      call take_snapshot(ea, sa0)
      call advance(ea, cfg_a, t_a, na, ok)
      call check(error, ok, "engine A failed after the checkpoint")
      if (allocated(error)) return
      call take_snapshot(ea, sa)
      !$acc update self(ea%state%ice%u_ice)
      umax = maxval(abs(ea%state%ice%u_ice))
      call finish_engine(ea)
      call check(error, umax > 1.0e-6_wp, "the ice moves in the restart case")
      if (allocated(error)) return

      call make_engine(eb, cfg_b, case_nml("channel_cice"), ok, restart_file=FN, t0=t_b, step0=step_b)
      call check(error, ok, "engine B (warm restart) setup failed")
      if (allocated(error)) return
      call check(error, step_b == nw .and. eb%warm_restart, "warm restart taken at the written step")
      if (allocated(error)) return
      call take_snapshot(eb, sb0)
      call compare(error, sa0, sb0, n_dync, "at the resume point")
      if (allocated(error)) then
         call finish_engine(eb)
         return
      end if
      call check(error, n_dync == 7, "the seven C-grid tags are part of the comparison")
      if (allocated(error)) return
      call advance(eb, cfg_b, t_b, na, ok)
      call check(error, ok, "engine B failed after the restart")
      if (allocated(error)) return
      call take_snapshot(eb, sb)
      call finish_engine(eb)
      call delete_file(FN)
      call compare(error, sa, sb, n_dync, "after the resumed steps")
   end subroutine test_restart

   function replace(s, a, b) result(r)
      character(len=*), intent(in) :: s, a, b
      character(len=:), allocatable :: r
      integer :: p
      p = index(s, a)
      if (p == 0) then
         r = s
      else
         r = s(1:p - 1)//b//s(p + len(a):)
      end if
   end function replace

   subroutine delete_file(fn_)
      character(len=*), intent(in) :: fn_
      integer :: u, stat
      open (newunit=u, file=fn_, status="old", iostat=stat)
      if (stat == 0) close (u, status="delete")
   end subroutine delete_file

   ! ------------------------------------------------------------------
   ! against the SIS2 core
   ! ------------------------------------------------------------------

   subroutine run_stats(label, nsteps, uavg, vavg, ke, ok)
      !! Mean ice velocity over the interior and total kinetic energy of the ice.
      character(len=*), intent(in) :: label
      integer, intent(in) :: nsteps
      real(wp), intent(out) :: uavg, vavg, ke
      logical, intent(out) :: ok
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      real(wp) :: t, nsum
      integer :: i, j, ng, is, ie, js, je

      uavg = 0.0_wp
      vavg = 0.0_wp
      ke = 0.0_wp
      call make_engine(engine, cfg, case_nml(label), ok)
      if (.not. ok) return
      t = 0.0_wp
      call advance(engine, cfg, t, nsteps, ok)
      if (.not. ok) return
      call pull_ice(engine)
      !$acc update self(engine%state%ice%vice, engine%state%ice%aice)
      ng = engine%grid%nghost
      is = ng + 1
      ie = ng + engine%grid%nx_phys
      js = ng + 1
      je = ng + engine%grid%ny_phys
      nsum = 0.0_wp
      associate (ice => engine%state%ice, m => engine%state%metrics)
         do j = js + 2, je - 2
            do i = is + 2, ie - 2
               uavg = uavg + 0.5_wp*(ice%u_ice(i, j) + ice%u_ice(i + 1, j))
               vavg = vavg + 0.5_wp*(ice%v_ice(i, j) + ice%v_ice(i, j + 1))
               nsum = nsum + 1.0_wp
            end do
         end do
         uavg = uavg/nsum
         vavg = vavg/nsum
         do j = js, je
            do i = is, ie
               ke = ke + 0.5_wp*ICE_RHO_ICE*ice%vice(i, j)*m%areaT(i, j)* &
                    (0.25_wp*(ice%u_ice(i, j) + ice%u_ice(i + 1, j))**2 + &
                     0.25_wp*(ice%v_ice(i, j) + ice%v_ice(i, j + 1))**2)
            end do
         end do
      end associate
      call finish_engine(engine)
   end subroutine run_stats

   ! ------------------------------------------------------------------
   ! experiments: why the two cores differ (frozen ocean, direct steps)
   ! ------------------------------------------------------------------

   subroutine frozen_stats(label, nsteps, margin, uavg, vavg, ke, ok)
      !! `nsteps` ice steps of the core selected by `label` on a FROZEN ocean
      !! (the ocean is never stepped, so both cores see the same forcing and the
      !! only difference is the ice dynamics).  Mean ice velocity of the
      !! interior (`margin` cells off the walls) and the total ice kinetic energy.
      character(len=*), intent(in) :: label
      integer, intent(in) :: nsteps, margin
      real(wp), intent(out) :: uavg, vavg, ke
      logical, intent(out) :: ok
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      type(ice_err_summary_t) :: summ
      integer :: n, i, j, ng, is, ie, js, je
      real(wp) :: nsum

      uavg = 0.0_wp
      vavg = 0.0_wp
      ke = 0.0_wp
      call make_engine(engine, cfg, case_nml(label), ok)
      if (.not. ok) return
      associate (st => engine%state, gr => engine%grid)
         do n = 1, nsteps
            if (st%ice%dync%is_init) then
               call ice_dyn_c_step(gr, st%metrics, st%coriolis_adv%f_corner, st%barotropic%b, &
                                   st%ice_cgrid, st%ice, st%multilayer, st%bc, DT, err_summary=summ)
               ok = summ%mask == 0_int32
               if (.not. ok) return
            else
               call ice_evp_step(gr, st%metrics, st%coriolis_adv%f_corner, st%ice, st%multilayer, DT, &
                                 engine%evp_params, st%bc)
            end if
         end do
         call pull_ice(engine)
         !$acc update self(st%ice%vice)
         ng = gr%nghost
         is = ng + 1
         ie = ng + gr%nx_phys
         js = ng + 1
         je = ng + gr%ny_phys
         nsum = 0.0_wp
         do j = js + margin, je - margin
            do i = is + margin, ie - margin
               uavg = uavg + 0.5_wp*(st%ice%u_ice(i, j) + st%ice%u_ice(i + 1, j))
               vavg = vavg + 0.5_wp*(st%ice%v_ice(i, j) + st%ice%v_ice(i, j + 1))
               nsum = nsum + 1.0_wp
            end do
         end do
         uavg = uavg/nsum
         vavg = vavg/nsum
         do j = js, je
            do i = is, ie
               ke = ke + 0.5_wp*ICE_RHO_ICE*st%ice%vice(i, j)*st%metrics%areaT(i, j)* &
                    (0.25_wp*(st%ice%u_ice(i, j) + st%ice%u_ice(i + 1, j))**2 + &
                     0.25_wp*(st%ice%v_ice(i, j) + st%ice%v_ice(i, j + 1))**2)
            end do
         end do
      end associate
      call finish_engine(engine)
   end subroutine frozen_stats

   function analytic_free_drift(tau, f, m, rc) result(vel)
      !! Steady free drift `tau = rc |u| u + m f k x u` (`tau` along x): the speed `s`
      !! solves `s^2 ((rc s)^2 + (m f)^2) = tau^2` (bisection), then
      !! `u = tau rc s / ((rc s)^2 + (m f)^2)`, `v = -m f u / (rc s)`.
      real(wp), intent(in) :: tau, f, m, rc
      real(wp) :: vel(2)
      real(wp) :: lo, hi, s, d
      integer :: it

      lo = 0.0_wp
      hi = 10.0_wp
      do it = 1, 200
         s = 0.5_wp*(lo + hi)
         if (s*s*((rc*s)**2 + (m*f)**2) < tau*tau) then
            lo = s
         else
            hi = s
         end if
      end do
      s = 0.5_wp*(lo + hi)
      d = (rc*s)**2 + (m*f)**2
      vel(1) = tau*rc*s/d
      vel(2) = -m*f*vel(1)/(rc*s)
   end function analytic_free_drift

   subroutine test_exp1(error)
      !! E1: closed basin with the ice strength switched off in BOTH cores:
      !! free drift against walls, which the two must agree on.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: ua, va, ka, ub, vb, kb
      logical :: ok

      call ex_clear()
      call ex_add("ocean_ice_nml", "p0", "1.0e-3")
      call ex_add("ocean_ice_ridge_nml", "pstar", "1.0e-3")
      call frozen_stats("basin", 72, 3, ua, va, ka, ok)
      call check(error, ok, "sis2_c run failed")
      if (allocated(error)) return
      call frozen_stats("basin_cice", 72, 3, ub, vb, kb, ok)
      call ex_clear()
      call check(error, ok, "cice_c run failed")
      if (allocated(error)) return
      write (*, '(a,2es13.5,a,2es13.5,a,es10.3,a,es10.3)') "E1 no strength: sis2_c (u,v) =", ua, va, &
         "  cice_c =", ub, vb, "  |du|/|u| =", sqrt((ua - ub)**2 + (va - vb)**2)/sqrt(ua**2 + va**2), &
         "  KE ratio -1 =", kb/ka - 1.0_wp
      ! Measured 1.8e-4 (velocity) and 4.5e-3 (kinetic energy): with no strength the pack is in
      ! free drift against the walls and the two cores agree, so the coupling (wind, drag,
      ! tilt, masks) is the same in both.  A factor 3 here would be a coupling bug.
      call check(error, ieee_is_finite(kb/ka) .and. ka > 0.0_wp, "finite, non-zero energy")
      if (allocated(error)) return
      call check(error, sqrt((ua - ub)**2 + (va - vb)**2)/sqrt(ua**2 + va**2) < 2.0e-3_wp, &
                 "E1: free drift against walls, velocity")
      if (allocated(error)) return
      call check(error, abs(kb/ka - 1.0_wp) < 2.0e-2_wp, "E1: free drift against walls, kinetic energy")
   end subroutine test_exp1

   subroutine test_exp2(error)
      !! E2: doubly periodic uniform pack, frozen ocean, 200 steps (steady state),
      !! without and with Coriolis, against the analytic free drift.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: ua, va, ka, ub, vb, kb, ref(2), f
      real(wp) :: rc, m, e_sis, e_cice
      integer :: icase
      logical :: ok
      character(len=*), parameter :: NL = new_line("a")

      rc = 1026.0_wp*5.36e-3_wp
      m = ICE_RHO_ICE
      do icase = 1, 2
         f = merge(0.0_wp, 1.4e-4_wp, icase == 1)
         call ex_clear()
         call ex_add("physics_nml", "coriolis_f", merge("0.0   ", "1.4e-4", icase == 1))
         call frozen_stats("uniform", 200, 0, ua, va, ka, ok)
         call check(error, ok, "sis2_c run failed")
         if (allocated(error)) return
         call frozen_stats("uniform_cice", 200, 0, ub, vb, kb, ok)
         call check(error, ok, "cice_c run failed")
         if (allocated(error)) return
         ref = analytic_free_drift(0.1_wp, f, m, rc)
         write (*, '(a,i2,a,2es13.5,a,2es13.5,a,2es13.5)') "E2 case", icase, " analytic", ref, "  sis2_c", ua, va, &
            "  cice_c", ub, vb
         e_sis = sqrt((ua - ref(1))**2 + (va - ref(2))**2)/sqrt(sum(ref**2))
         e_cice = sqrt((ub - ref(1))**2 + (vb - ref(2))**2)/sqrt(sum(ref**2))
         write (*, '(a,2es10.3)') "E2   relative error |u-ref|/|ref|: sis2_c, cice_c", e_sis, e_cice
         ! Measured: f = 0: 8e-16 both; f /= 0: cice_c 7e-15, sis2_c 3.7e-4 (its own Coriolis
         ! time discretisation; it is converged: 200 steps).
         call check(error, e_cice < 1.0e-12_wp, "E2: cice_c is on the analytic free drift")
         if (allocated(error)) return
         if (icase == 1) then
            call check(error, e_sis < 1.0e-12_wp, "E2: sis2_c is on the analytic free drift (f = 0)")
         else
            call check(error, e_sis < 2.0e-3_wp, "E2: sis2_c is near the analytic free drift (f /= 0)")
         end if
         if (allocated(error)) return
         ! the 6-step transient of `vs_sis2_free_drift`
         call frozen_stats("uniform", 6, 0, ua, va, ka, ok)
         call frozen_stats("uniform_cice", 6, 0, ub, vb, kb, ok)
         write (*, '(a,i2,a,es10.3)') "E2 case", icase, ": 6-step relative difference sis2_c/cice_c ", &
            sqrt((ua - ub)**2 + (va - vb)**2)/sqrt(ua**2 + va**2)
      end do
      call ex_clear()
      call check(error, ieee_is_finite(ub), "finite")
   end subroutine test_exp2

   subroutine test_exp3(error)
      !! E3: closed basin, rheology on, frozen ocean: the viscosity floor and the
      !! elastic damping of the new core moved toward the existing core's
      !! (floor `del_sh_min_pr p0/rho = 2 s dt_sub^2/(Tdamp dxharm^2) p0/rho`).
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: ua, va, ka, ub, vb, kb, dmin, k120, kv(6), r120
      logical :: ok
      character(len=*), parameter :: NL = new_line("a")
      integer :: ivar
      character(len=64) :: lab

      ! existing core: scale 2, dt_sub = 600/20, Tdamp = max(0.2*600, 3 dt_sub) = 120 s,
      ! dxharm = 20 km, p0/rho_ice = 27500/905
      dmin = 2.0_wp*2.0_wp*30.0_wp**2/(120.0_wp*(2.0e4_wp)**2)*(2.75e4_wp/ICE_RHO_ICE)
      write (*, '(a,es11.4)') "E3 existing core's viscosity floor [1/s] =", dmin
      call ex_clear()
      call frozen_stats("basin", 72, 3, ua, va, ka, ok)
      call check(error, ok, "sis2_c run failed")
      if (allocated(error)) return
      write (*, '(a,es12.4)') "E3 sis2_c (20 subcycles) KE =", ka
      do ivar = 1, 6
         select case (ivar)
         case (1)
            call ex_clear()
            lab = "default"
         case (2)
            write (lab, '(a,es9.2)') "deltaminEVP =", dmin
            call ex_clear()
            call ex_add("ocean_ice_dyn_nml", "deltaminEVP", trim(adjustl(num(dmin))))
         case (3)
            lab = "deltamin + elasticDamp 0.2 (Tdamp 120 s)"
            call ex_clear()
            call ex_add("ocean_ice_dyn_nml", "deltaminEVP", trim(adjustl(num(dmin))))
            call ex_add("ocean_ice_dyn_nml", "elasticDamp", "0.2")
         case (4)
            lab = "deltamin + elasticDamp 0.2 + ndte 120"
            call ex_clear()
            call ex_add("ocean_ice_dyn_nml", "deltaminEVP", trim(adjustl(num(dmin))))
            call ex_add("ocean_ice_dyn_nml", "elasticDamp", "0.2")
            call ex_add("ocean_ice_dyn_nml", "ndte", "120")
         case (5)
            lab = "deltamin + elasticDamp 0.2 + ndte 120, sis2_c also 120"
            call ex_clear()
            call ex_add("ocean_ice_dyn_nml", "deltaminEVP", trim(adjustl(num(dmin))))
            call ex_add("ocean_ice_dyn_nml", "elasticDamp", "0.2")
            call ex_add("ocean_ice_dyn_nml", "ndte", "120")
            call ex_add("ocean_ice_nml", "evp_sub_steps", "120")
         case (6)
            lab = "ndte 120, floor of the existing core at 120 subcycles, elasticDamp 0.2; sis2_c also 120"
            call ex_clear()
            call ex_add("ocean_ice_dyn_nml", "deltaminEVP", trim(adjustl(num(dmin/36.0_wp))))
            call ex_add("ocean_ice_dyn_nml", "elasticDamp", "0.2")
            call ex_add("ocean_ice_dyn_nml", "ndte", "120")
            call ex_add("ocean_ice_nml", "evp_sub_steps", "120")
         end select
         call frozen_stats("basin_cice", 72, 3, ub, vb, kb, ok)
         call check(error, ok, "cice_c run failed")
         if (allocated(error)) return
         kv(ivar) = kb
         write (*, '(a,a,a,es14.7,a,f8.4)') "E3 cice_c ", trim(lab), ": KE =", kb, "  KE/KE(sis2_c 20) =", kb/ka
         if (ivar == 5 .or. ivar == 6) then
            call frozen_stats("basin", 72, 3, ua, va, k120, ok)
            write (*, '(a,es14.7,a,f8.4)') "E3   sis2_c with 120 subcycles: KE =", k120, "  cice_c/that =", kb/k120
            if (ivar == 6) r120 = kb/k120
         end if
      end do
      call ex_clear()
      ! Measured: default 0.0899 (the new core is close to rigid at deltaminEVP = 1e-11); with the
      ! existing core's floor 1.137; ndte (20 -> 120) and elasticDamp move the kinetic energy by
      ! < 1e-6 relative (a steady state does not depend on them, so unconverged subcycling is NOT
      ! the cause); floors matched at 120 subcycles: 1.37 (a residual of 14-37 % remains: the two
      ! regularisations are not the same function of the strain rate).
      call check(error, kv(1)/ka > 0.05_wp .and. kv(1)/ka < 0.15_wp, "E3: default floor is far too stiff")
      if (allocated(error)) return
      call check(error, kv(2)/ka > 0.9_wp .and. kv(2)/ka < 1.4_wp, "E3: matched floor brings the energy to the existing core's")
      if (allocated(error)) return
      call check(error, abs(kv(3)/kv(2) - 1.0_wp) < 1.0e-6_wp .and. abs(kv(4)/kv(2) - 1.0_wp) < 1.0e-6_wp, &
                 "E3: elasticDamp and ndte do not change the steady state")
      if (allocated(error)) return
      call check(error, r120 > 1.0_wp .and. r120 < 1.8_wp, "E3: matched floors at 120 subcycles")
   contains
      function num(x) result(r)
         real(wp), intent(in) :: x
         character(len=32) :: r
         write (r, '(es14.6)') x
      end function num
   end subroutine test_exp3

   subroutine test_restart_cross(error)
      !! A checkpoint written by one core is refused by the other, both ways; a
      !! file without the scheme id is the existing core's.
      type(error_type), allocatable, intent(out) :: error

      call cross_one(error, "channel_cice", "channel")
      if (allocated(error)) return
      call cross_one(error, "channel", "channel_cice")
   end subroutine test_restart_cross

   subroutine cross_one(error, lab_w, lab_r)
      !! Write a checkpoint with the core of `lab_w`, try to resume it with that of `lab_r`
      !! (fresh engines and configurations each call).
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: lab_w, lab_r
      type(ocean_engine_t), target :: ea, eb
      type(config_t) :: cfg_a, cfg_b
      real(wp) :: t_a, t_b
      integer :: ierr, step_b
      logical :: ok

      call make_engine(ea, cfg_a, case_nml(lab_w), ok)
      call check(error, ok, "writer setup failed")
      if (allocated(error)) return
      t_a = 0.0_wp
      call advance(ea, cfg_a, t_a, 2, ok)
      call check(error, ok, "writer steps failed")
      if (allocated(error)) return
      call ocean_state_restart_write(ea%state, ea%grid, ea%decomp, FN, t_a, 2, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "checkpoint write failed")
      if (allocated(error)) return
      call finish_engine(ea)
      call ensure_comm()
      call read_config_from_string(case_nml(lab_r), cfg_b, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "reader config failed")
      if (allocated(error)) return
      call engine_setup(eb, cfg_b, ierr, restart_file=FN, t_restart=t_b, step_restart=step_b)
      call check(error, ierr /= OCEAN_STATUS_OK, "a checkpoint of the other dynamics core was NOT refused: "//lab_w//" -> "//lab_r)
      call delete_file(FN)
      ! No `engine_teardown`: the refused setup stopped at the restart read, before the
      ! halo, the device map and the diagnostics exist, and tearing those down crashes.
   end subroutine cross_one

   subroutine test_throughput(error)
      !! Opt-in (`RDB_CICE_EVP_BENCH=1`): wall time of `ice_dyn_c_step` on a 512 x 512
      !! closed basin with 100 % ice (ncat = 5) at the upstream default of 120
      !! subcycles; prints ms per step and the cell-subcycle rate.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      character(len=8) :: flag
      character(len=:), allocatable :: nml
      logical :: ok
      integer :: n, nstep
      integer(int64) :: c0, c1, crate
      real(wp) :: ms
      type(ice_err_summary_t) :: summ

      call get_environment_variable("RDB_CICE_EVP_BENCH", flag)
      if (trim(flag) /= "1") return
      nml = replace(replace(case_nml("basin_cice"), "nx = 24, ny = 16", "nx = 512, ny = 512"), "ndte = 20", "ndte = 120")
      call make_engine(engine, cfg, nml, ok)
      call check(error, ok, "engine setup failed")
      if (allocated(error)) return
      nstep = 5
      associate (st => engine%state, gr => engine%grid)
         call ice_dyn_c_step(gr, st%metrics, st%coriolis_adv%f_corner, st%barotropic%b, st%ice_cgrid, &
                             st%ice, st%multilayer, st%bc, DT, err_summary=summ)
         call system_clock(c0, crate)
         do n = 1, nstep
            call ice_dyn_c_step(gr, st%metrics, st%coriolis_adv%f_corner, st%barotropic%b, st%ice_cgrid, &
                                st%ice, st%multilayer, st%bc, DT, err_summary=summ)
         end do
         call system_clock(c1)
         ms = 1.0e3_wp*real(c1 - c0, wp)/real(crate, wp)/real(nstep, wp)
         write (*, '(a,f10.3,a,es11.3,a)') "step_throughput 512x512, ndte=120: ", ms, " ms/step, ", &
            512.0_wp*512.0_wp*120.0_wp/(ms*1.0e-3_wp), " cell-subcycles/s"
      end associate
      call finish_engine(engine)
   end subroutine test_throughput

   subroutine test_vs_sis2_free_drift(error)
      !! A doubly periodic uniform pack under a steady wind: no strain, so the
      !! two cores differ only in how they discretise the drag, the Coriolis
      !! term and the subcycling.  The numbers are printed.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: ua, va, ka, ub, vb, kb, rel
      logical :: ok

      call run_stats("uniform", 6, ua, va, ka, ok)
      call check(error, ok, "sis2_c run failed")
      if (allocated(error)) return
      call run_stats("uniform_cice", 6, ub, vb, kb, ok)
      call check(error, ok, "cice_c run failed")
      if (allocated(error)) return
      rel = sqrt((ua - ub)**2 + (va - vb)**2)/sqrt(ua**2 + va**2)
      write (*, '(a,2es13.5,a,2es13.5,a,es10.3)') "vs_sis2_free_drift: sis2_c (u,v) =", ua, va, &
         "  cice_c (u,v) =", ub, vb, "  relative difference ", rel
      call check(error, ieee_is_finite(rel), "finite comparison")
      if (allocated(error)) return
      ! Measured 3.7e-3 after 6 steps: a spin-up transient (the Coriolis and drag time
      ! discretisations differ).  At steady state (E2, 200 steps) the cice_c core is on the
      ! analytic free drift to 7e-15 and the SIS2 core to 3.7e-4.
      call check(error, rel < 1.0e-2_wp, "the two C-grid cores agree on a uniform pack under wind")
   end subroutine test_vs_sis2_free_drift

   subroutine test_vs_sis2_basin(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: ua, va, ka, ub, vb, kb, ratio
      logical :: ok

      call run_stats("basin", 72, ua, va, ka, ok)
      call check(error, ok, "sis2_c run failed")
      if (allocated(error)) return
      call run_stats("basin_cice", 72, ub, vb, kb, ok)
      call check(error, ok, "cice_c run failed")
      if (allocated(error)) return
      ratio = kb/ka
      write (*, '(a,2es13.5,a,2es13.5,a,2es12.4,a,f8.4)') "vs_sis2_basin: sis2_c mean (u,v) =", ua, va, &
         "  cice_c =", ub, vb, "  KE (sis2_c, cice_c) =", ka, kb, "  ratio ", ratio
      call check(error, ieee_is_finite(ratio) .and. ka > 0.0_wp, "finite, non-zero kinetic energy")
      if (allocated(error)) return
      ! Pinned to the measured 0.097, and the reason is E3 (`exp3_basin_viscosity_floor`): the
      ! CICE core regularises the viscosities at deltaminEVP = 1e-11 1/s, the SIS2 core at
      ! 2 dt_sub^2 s/(Tdamp dxharm^2) p0/rho = 2.3e-6 1/s in this case, so the compact pack
      ! creeps in one and is close to rigid in the other.  With the floors matched the ratio
      ! is 1.14.  Both cores have a no-slip coast.
      call check(error, ratio > 0.05_wp .and. ratio < 0.15_wp, "basin kinetic energy ratio (floor-dominated, E3)")
   end subroutine test_vs_sis2_basin

end module test_cice_dyn_c_engine
