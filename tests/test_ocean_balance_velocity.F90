!! Tests for `ocean_balance_velocity` (geostrophic mode, PR-B1).
!!
!! The reference is the balanced-vortex seed (`ic_config = "balanced_vortex"`).
!! The same vortex is built two ways through the PRODUCTION engine:
!!   * `vortex_balance = "geostrophic"`: the seed's analytic velocity;
!!   * `vortex_balance = "none"`: the same density and SSH at rest, then
!!     `ocean_balance_velocity(BALANCE_GEOSTROPHIC)` from the model's own PGF.
!! They are two discretisations of one balance, so they agree to the
!! difference between the seed's rest-thickness hydrostatic integral and the
!! model's Montgomery PGF on the stretched columns, which is O(eta/H) and is
!! pinned below by showing it scales LINEARLY with the vortex amplitude.
!!
!! GPU contract (`-gpu=mem:separate`): the routine runs after `enter_data`;
!! inputs are host-modified then pushed by the routine, and the results are
!! read back with `update self` on LEAF arrays, so on a GPU build the
!! comparison reads what the device holds.
module test_ocean_balance_velocity
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_config, only: config_t, read_config_from_string, validate_config
   use rdb_ocean_engine, only: ocean_engine_t, engine_setup, engine_enter_data, &
                               engine_exit_data, engine_teardown, engine_step, &
                               engine_step_finalize
   use rdb_ocean_balance, only: ocean_balance_velocity, BALANCE_NONE, &
                                BALANCE_GEOSTROPHIC, BALANCE_GRADIENT
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_IC_SEED
   use rdb_comm_env, only: comm_env_init, comm_env_setup_roles
   implicit none
   private

   public :: collect_ocean_balance_velocity_tests

   integer, parameter :: NXP = 32
   integer, parameter :: NG = 4
   real(wp), parameter :: DX = 10000.0_wp
   real(wp), parameter :: DTDZ = 0.03_wp
      !! Background dT/dz (K/m): statically stable for every |dT| <= 12 K here.
   logical, save :: comm_inited = .false.

contains

   subroutine collect_ocean_balance_velocity_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("balance_fplane_matches_vortex_seed", test_fplane), &
                  new_unittest("balance_betaplane_matches_vortex_seed", test_betaplane), &
                  new_unittest("balance_difference_is_interface_tilt", test_stratification_scaling), &
                  new_unittest("balance_is_steadier_than_seed", test_steadier), &
                  new_unittest("balance_bottom_rest_recovers_seed_ssh", test_bottom_rest), &
                  new_unittest("balance_refusals_leave_state_untouched", test_refusals) &
                  ]
   end subroutine collect_ocean_balance_velocity_tests

   subroutine ensure_comm()
      if (.not. comm_inited) then
         call comm_env_init()
         call comm_env_setup_roles(.false.)
         comm_inited = .true.
      end if
   end subroutine ensure_comm

   function vortex_nml(beta, balance, dT, dtdz) result(nml)
      real(wp), intent(in) :: beta, dT, dtdz
      character(len=*), intent(in) :: balance
      character(len=:), allocatable :: nml
      character(len=*), parameter :: NL = new_line("a")
      character(len=32) :: s_dt, s_beta, s_yref, s_dtdz
      write (s_dtdz, "(es14.6)") dtdz
      write (s_dt, "(es14.6)") dT
      write (s_beta, "(es14.6)") beta
      write (s_yref, "(es14.6)") 0.5_wp*real(NXP, wp)*DX
      nml = "&sim_nml sim_type = 'ocean' /"//NL// &
            "&grid_nml nx = 32, ny = 32, dx = 10000.0, dy = 10000.0, nghost = 4 /"//NL// &
            "&nonhydrostatic_nml nz_layers = 6 /"//NL// &
            "&time_nml t_end = 1.0e9, dt_fixed = 1200.0 /"//NL// &
            "&physics_nml coriolis_f = 1.0e-4, wind_stress_x = 0.0, wind_stress_y = 0.0 /"//NL// &
            "&tracer_nml initial_salinity = 35.0 /"//NL// &
            "&ocean_topo_nml topo_config = 'flat', max_depth = 2000.0, "// &
            "wind_config = 'constant', coriolis_beta = "//trim(adjustl(s_beta))// &
            ", coriolis_y_ref = "//trim(adjustl(s_yref))//" /"//NL// &
            "&ocean_thermo_nml enable_thermodynamics = .true. /"//NL// &
            "&ocean_eos_nml eos = 'linear' /"//NL// &
            "&ocean_ic_nml ic_config = 'balanced_vortex', alpha_T = 0.17, rho_0 = 1035.0, "// &
            "vortex_radius = 40000.0, vortex_hv = 500.0, vortex_dT_dz = "//trim(adjustl(s_dtdz))//", vortex_dT = "// &
            trim(adjustl(s_dt))//", vortex_balance = '"//balance//"' /"//NL// &
            "&ocean_hvisc_nml nu_h = 350.0 /"//NL// &
            "&ocean_vmix_nml use_closure = .false. /"//NL// &
            "&ocean_bt_nml split_scheme = 'pred_corr', auto_n_inner = .true. /"//NL// &
            "&vcoord_nml vcoord_type = 'sigma' /"//NL// &
            "&ocean_diag_nml enabled = .false. /"//NL// &
            "&output_nml output_to_file = .false. /"//NL
   end function vortex_nml

   subroutine setup_engine(engine, cfg, beta, balance, dT, dtdz, ok)
      type(ocean_engine_t), intent(inout) :: engine
      type(config_t), intent(out) :: cfg
      real(wp), intent(in) :: beta, dT, dtdz
      character(len=*), intent(in) :: balance
      logical, intent(out) :: ok
      integer :: ierr
      ok = .false.
      call ensure_comm()
      call read_config_from_string(vortex_nml(beta, balance, dT, dtdz), cfg, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call validate_config(cfg, ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call engine_setup(engine, cfg, ierr)
      ok = (ierr == OCEAN_STATUS_OK)
   end subroutine setup_engine

   subroutine pull_uvh(engine)
      !! Device -> host on LEAF arrays: what the device holds is what is read.
      type(ocean_engine_t), intent(inout) :: engine
      associate (uf => engine%state%multilayer%u_face_x_layer, &
                 vf => engine%state%multilayer%v_face_y_layer, &
                 hl => engine%state%multilayer%h_layer)
         !$omp target update from(uf, vf, hl)
      end associate
   end subroutine pull_uvh

   subroutine seed_fields(beta, dT, dtdz, u, v, h)
      !! The analytic seed's fields (host, straight after `engine_setup`).
      real(wp), intent(in) :: beta, dT, dtdz
      real(wp), allocatable, intent(out) :: u(:, :, :), v(:, :, :), h(:, :, :)
      type(ocean_engine_t) :: engine
      type(config_t) :: cfg
      logical :: ok
      call setup_engine(engine, cfg, beta, "geostrophic", dT, dtdz, ok)
      if (.not. ok) then
         allocate (u(1, 1, 1), v(1, 1, 1), h(1, 1, 1), source=-1.0_wp)
         return
      end if
      u = engine%state%multilayer%u_face_x_layer
      v = engine%state%multilayer%v_face_y_layer
      h = engine%state%multilayer%h_layer
      call engine_teardown(engine)
   end subroutine seed_fields

   subroutine balanced_fields(beta, dT, dtdz, ref_level, flatten, u, v, h, ierr)
      !! The `none` seed, mapped to the device, then `ocean_balance_velocity`.
      !! `flatten`: set every column to `sum h = b` (flat SSH) BEFORE the call,
      !! scaling `hTr` with `h` (concentration kept).
      real(wp), intent(in) :: beta, dT, dtdz
      character(len=*), intent(in) :: ref_level
      logical, intent(in) :: flatten
      real(wp), allocatable, intent(out) :: u(:, :, :), v(:, :, :), h(:, :, :)
      integer, intent(out) :: ierr
      type(ocean_engine_t) :: engine
      type(config_t) :: cfg
      logical :: ok
      integer :: i, j, it, nx, ny, nz
      real(wp) :: s

      ierr = -1
      call setup_engine(engine, cfg, beta, "none", dT, dtdz, ok)
      if (.not. ok) then
         allocate (u(1, 1, 1), v(1, 1, 1), h(1, 1, 1), source=-1.0_wp)
         return
      end if
      call engine_enter_data(engine, cfg)
      if (flatten) then
         nx = engine%grid%nx_total
         ny = engine%grid%ny_total
         nz = engine%state%multilayer%nz_ml
         do j = 1, ny
            do i = 1, nx
               s = engine%state%barotropic%b(i, j)/sum(engine%state%multilayer%h_layer(i, j, :))
               engine%state%multilayer%h_layer(i, j, :) = s*engine%state%multilayer%h_layer(i, j, :)
               do it = 1, size(engine%state%multilayer%tracers)
                  engine%state%multilayer%tracers(it)%hTr(i, j, :) = &
                     s*engine%state%multilayer%tracers(it)%hTr(i, j, :)
               end do
            end do
         end do
      end if
      call ocean_balance_velocity(engine%grid, cfg, engine%state, BALANCE_GEOSTROPHIC, ierr, &
                                  ref_level=ref_level)
      call pull_uvh(engine)
      u = engine%state%multilayer%u_face_x_layer
      v = engine%state%multilayer%v_face_y_layer
      h = engine%state%multilayer%h_layer
      call engine_exit_data(engine)
      call engine_teardown(engine)
   end subroutine balanced_fields

   pure function rel_diff(a, b) result(r)
      !! max |a - b| / max |b| over the interior, one cell in from the walls.
      real(wp), intent(in) :: a(:, :, :), b(:, :, :)
      real(wp) :: r
      integer :: lo, hi
      lo = NG + 2
      hi = NG + NXP - 1
      r = maxval(abs(a(lo:hi, lo:hi, :) - b(lo:hi, lo:hi, :)))/maxval(abs(b(lo:hi, lo:hi, :)))
   end function rel_diff

   subroutine compare(error, beta, dT, dtdz, tol, label)
      type(error_type), allocatable, intent(out) :: error
      real(wp), intent(in) :: beta, dT, dtdz, tol
      character(len=*), intent(in) :: label
      real(wp), allocatable :: us(:, :, :), vs(:, :, :), hs(:, :, :)
      real(wp), allocatable :: ub(:, :, :), vb(:, :, :), hb(:, :, :)
      integer :: ierr
      real(wp) :: eu, ev, umax
      character(len=200) :: msg

      call seed_fields(beta, dT, dtdz, us, vs, hs)
      call balanced_fields(beta, dT, dtdz, "keep_h", .false., ub, vb, hb, ierr)
      call check(error, ierr, OCEAN_STATUS_OK, "balance returned an error status")
      if (allocated(error)) return
      umax = maxval(abs(us))
      eu = rel_diff(ub, us)
      ev = rel_diff(vb, vs)
      write (*, "(3a,es10.3,a,es10.3,a,es10.3)") "  [balance] ", label, ": max|u_seed| = ", umax, &
         " m/s; rel diff u = ", eu, ", v = ", ev
      call check(error, umax > 0.01_wp, "seed velocity is a real mesoscale signal")
      if (allocated(error)) return
      write (msg, "(2a,es10.3,a,es10.3)") label, ": u/v differ from the seed by ", max(eu, ev), &
         " > tol ", tol
      call check(error, max(eu, ev) <= tol, trim(msg))
      if (allocated(error)) return
      ! keep_h must not touch the thickness
      call check(error, maxval(abs(hb - hs)) == 0.0_wp, "keep_h leaves h_layer untouched")
   end subroutine compare

   subroutine test_fplane(error)
      type(error_type), allocatable, intent(out) :: error
      call compare(error, 0.0_wp, 1.0_wp, DTDZ, 1.0e-2_wp, "f-plane")
   end subroutine test_fplane

   subroutine test_betaplane(error)
      type(error_type), allocatable, intent(out) :: error
      call compare(error, 2.0e-11_wp, 1.0_wp, DTDZ, 1.0e-2_wp, "beta-plane")
   end subroutine test_betaplane

   function diff_for(dT, dtdz) result(e)
      !! Seed-vs-routine relative velocity difference at one (dT, dT/dz).
      real(wp), intent(in) :: dT, dtdz
      real(wp) :: e
      real(wp), allocatable :: us(:, :, :), vs(:, :, :), ub(:, :, :), vb(:, :, :), hb(:, :, :)
      integer :: ierr
      call seed_fields(0.0_wp, dT, dtdz, us, vs, hb)
      call balanced_fields(0.0_wp, dT, dtdz, "keep_h", .false., ub, vb, hb, ierr)
      e = rel_diff(ub, us)
   end function diff_for

   subroutine test_stratification_scaling(error)
      !! WHY the routine and the seed differ by ~0.5 %, pinned by scaling.
      !! The seed's potential integrates the density at the REST layer
      !! thickness, as if the layers were flat.  In the model the layers
      !! stretch with the SSH, so a column with eta > 0 has its interfaces
      !! displaced by O(eta); across them the BACKGROUND stratification
      !! contributes a real interface-tilt term, g/rho0 (rho_k - rho_{k+1})
      !! grad(z_interface), that the model's PGF carries and the seed omits.
      !! That term is proportional to the SSH (so to dT, like the velocity
      !! itself: the RELATIVE difference is independent of dT) and to the
      !! background density jump (so to dT/dz: halving dT/dz halves it, up
      !! to a small remainder from the anomaly's own stretching).
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: e_ref, e_half_dT, e_half_N, r_dT, r_N
      character(len=160) :: msg

      e_ref = diff_for(1.0_wp, DTDZ)
      e_half_dT = diff_for(0.5_wp, DTDZ)
      e_half_N = diff_for(1.0_wp, 0.5_wp*DTDZ)
      r_dT = e_ref/e_half_dT
      r_N = e_ref/e_half_N
      write (*, "(a,es10.3,a,es10.3,a,es10.3,a,f5.2,a,f5.2)") "  [balance] rel diff (dT=1, dTdz): ", &
         e_ref, "  (dT/2): ", e_half_dT, "  (dTdz/2): ", e_half_N, "  ratios ", r_dT, ", ", r_N
      write (msg, "(a,f5.2)") "relative difference must not depend on dT (ratio ~1); ratio ", r_dT
      call check(error, r_dT > 0.85_wp .and. r_dT < 1.15_wp, trim(msg))
      if (allocated(error)) return
      write (msg, "(a,f5.2)") "halving dT/dz must cut the difference ~2x; ratio ", r_N
      call check(error, r_N > 1.6_wp .and. r_N < 2.6_wp, trim(msg))
   end subroutine test_stratification_scaling

   subroutine test_bottom_rest(error)
      !! From a FLAT sea surface, `bottom_rest` recovers the seed's SSH
      !! (zero bottom-pressure anomaly) and velocities.
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: us(:, :, :), vs(:, :, :), hs(:, :, :)
      real(wp), allocatable :: ub(:, :, :), vb(:, :, :), hb(:, :, :)
      integer :: ierr, i, j, lo, hi
      real(wp) :: eta_err, eta_peak, eu, ev
      character(len=200) :: msg

      call seed_fields(0.0_wp, 1.0_wp, DTDZ, us, vs, hs)
      call balanced_fields(0.0_wp, 1.0_wp, DTDZ, "bottom_rest", .true., ub, vb, hb, ierr)
      call check(error, ierr, OCEAN_STATUS_OK, "bottom_rest returned an error status")
      if (allocated(error)) return
      lo = NG + 2
      hi = NG + NXP - 1
      eta_err = 0.0_wp
      eta_peak = 0.0_wp
      do j = lo, hi
         do i = lo, hi
            eta_err = max(eta_err, abs(sum(hb(i, j, :)) - sum(hs(i, j, :))))
            eta_peak = max(eta_peak, abs(sum(hs(i, j, :)) - 2000.0_wp))
         end do
      end do
      eu = rel_diff(ub, us)
      ev = rel_diff(vb, vs)
      write (*, "(a,es10.3,a,es10.3,a,es10.3,a,es10.3)") "  [balance] bottom_rest: eta peak ", eta_peak, &
         " m, max|eta - eta_seed| ", eta_err, " m; rel diff u ", eu, " v ", ev
      write (msg, "(a,es10.3,a,es10.3)") "SSH differs from the seed by ", eta_err, &
         " m (peak eta ", eta_peak
      call check(error, eta_peak > 0.01_wp, "seed eta is non-trivial")
      if (allocated(error)) return
      call check(error, eta_err <= 1.0e-3_wp*eta_peak, trim(msg))
      if (allocated(error)) return
      call check(error, max(eu, ev) <= 1.0e-2_wp, "bottom_rest velocities match the seed")
   end subroutine test_bottom_rest

   function drift_after_steps(balance, use_routine, nsteps, umax_ref) result(d)
      !! max |u(n) - u(0)| / umax_ref over the interior and all layers after
      !! `nsteps` full model steps, from the `balance` seed; `use_routine`
      !! applies `ocean_balance_velocity` (keep_h) to it before stepping.
      character(len=*), intent(in) :: balance
      logical, intent(in) :: use_routine
      integer, intent(in) :: nsteps
      real(wp), intent(inout) :: umax_ref
         !! Normaliser; set from this run's initial |u| when <= 0.
      real(wp) :: d
      real(wp), parameter :: DT_STEP = 1200.0_wp
      type(ocean_engine_t) :: engine
      type(config_t) :: cfg
      logical :: ok
      integer :: ierr, n, lo, hi
      real(wp) :: t, umax
      real(wp), allocatable :: u0(:, :, :)

      d = huge(1.0_wp)
      call setup_engine(engine, cfg, 0.0_wp, balance, 1.0_wp, DTDZ, ok)
      if (.not. ok) return
      call engine_enter_data(engine, cfg)
      if (use_routine) then
         call ocean_balance_velocity(engine%grid, cfg, engine%state, BALANCE_GEOSTROPHIC, ierr, &
                                     ref_level="keep_h")
         if (ierr /= OCEAN_STATUS_OK) return
      end if
      call pull_uvh(engine)
      u0 = engine%state%multilayer%u_face_x_layer
      t = 0.0_wp
      do n = 1, nsteps
         call engine_step(engine, DT_STEP, t, ierr=ierr)
         if (ierr == OCEAN_STATUS_OK) call engine_step_finalize(engine, DT_STEP, t, ierr=ierr)
         if (ierr /= OCEAN_STATUS_OK) return
         t = t + DT_STEP
      end do
      call pull_uvh(engine)
      lo = NG + 2
      hi = NG + NXP - 1
      umax = umax_ref
      if (umax <= 0.0_wp) umax = maxval(abs(u0(lo:hi, lo:hi, :)))
      umax_ref = umax
      d = maxval(abs(engine%state%multilayer%u_face_x_layer(lo:hi, lo:hi, :) - u0(lo:hi, lo:hi, :)))/umax
      call engine_exit_data(engine)
      call engine_teardown(engine)
   end function drift_after_steps

   subroutine test_steadier(error)
      !! The discriminating control: stepping the model, the routine's state
      !! must drift no more than the analytic seed's (it balances the model's
      !! OWN discrete PGF, the seed balances an idealised one), and both must
      !! drift far less than the same density at rest (which radiates
      !! inertia-gravity waves).
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: d_seed, d_routine, d_rest, umax
      integer :: ns
      character(len=160) :: msg

      do ns = 1, 6, 5
         umax = 0.0_wp
         d_seed = drift_after_steps("geostrophic", .false., ns, umax)
         d_routine = drift_after_steps("none", .true., ns, umax)
         d_rest = drift_after_steps("none", .false., ns, umax)
         write (*, "(a,i0,a,es10.3,a,es10.3,a,es10.3)") "  [balance] velocity drift after ", ns, &
            " steps: seed ", d_seed, ", routine ", d_routine, ", unbalanced rest ", d_rest
      end do
      write (msg, "(a,es10.3,a,es10.3)") "routine drift ", d_routine, " exceeds the seed's ", d_seed
      call check(error, d_routine <= 1.02_wp*d_seed, trim(msg))
      if (allocated(error)) return
      call check(error, d_routine < 0.5_wp*d_rest, "balanced start must drift far less than rest")
   end subroutine test_steadier

   subroutine test_refusals(error)
      !! Every refusal returns its status and leaves u, v (and h) untouched.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t) :: engine
      type(config_t) :: cfg, cfg_bad
      logical :: ok
      integer :: ierr
      real(wp), allocatable :: u0(:, :, :), v0(:, :, :), h0(:, :, :)

      call setup_engine(engine, cfg, 0.0_wp, "geostrophic", 1.0_wp, DTDZ, ok)
      call check(error, ok, "engine_setup failed")
      if (allocated(error)) return
      call engine_enter_data(engine, cfg)
      u0 = engine%state%multilayer%u_face_x_layer
      v0 = engine%state%multilayer%v_face_y_layer
      h0 = engine%state%multilayer%h_layer

      checks: block
         ! gradient-wind: reserved
         call ocean_balance_velocity(engine%grid, cfg, engine%state, BALANCE_GRADIENT, ierr)
         call check(error, ierr, OCEAN_STATUS_ERR_IC_SEED, "gradient mode must fail loud")
         if (allocated(error)) exit checks
         call untouched(error, engine, u0, v0, h0, "gradient")
         if (allocated(error)) exit checks

         ! unknown mode
         call ocean_balance_velocity(engine%grid, cfg, engine%state, 7, ierr)
         call check(error, ierr, OCEAN_STATUS_ERR_IC_SEED, "unknown mode must fail loud")
         if (allocated(error)) exit checks

         ! bad ref_level
         call ocean_balance_velocity(engine%grid, cfg, engine%state, BALANCE_GEOSTROPHIC, ierr, &
                                     ref_level="nonsense")
         call check(error, ierr, OCEAN_STATUS_ERR_IC_SEED, "bad ref_level must fail loud")
         if (allocated(error)) exit checks
         call untouched(error, engine, u0, v0, h0, "ref_level")
         if (allocated(error)) exit checks

         ! equatorial: f_min above the domain's f
         call ocean_balance_velocity(engine%grid, cfg, engine%state, BALANCE_GEOSTROPHIC, ierr, &
                                     f_min=1.0e-3_wp)
         call check(error, ierr, OCEAN_STATUS_ERR_IC_SEED, "|f| < f_min must fail loud")
         if (allocated(error)) exit checks
         call untouched(error, engine, u0, v0, h0, "f_min")
         if (allocated(error)) exit checks

         ! unsupported grids
         cfg_bad = cfg
         cfg_bad%ocean%grid%grid_config = "tripolar"
         call ocean_balance_velocity(engine%grid, cfg_bad, engine%state, BALANCE_GEOSTROPHIC, ierr)
         call check(error, ierr, OCEAN_STATUS_ERR_IC_SEED, "tripolar must be refused")
         if (allocated(error)) exit checks
         call untouched(error, engine, u0, v0, h0, "tripolar")
         if (allocated(error)) exit checks
         cfg_bad%ocean%grid%grid_config = "supergrid"
         call ocean_balance_velocity(engine%grid, cfg_bad, engine%state, BALANCE_GEOSTROPHIC, ierr)
         call check(error, ierr, OCEAN_STATUS_ERR_IC_SEED, "supergrid must be refused")
         if (allocated(error)) exit checks

         ! wet/dry
         cfg_bad = cfg
         cfg_bad%ocean%wetdry%enable = .true.
         call ocean_balance_velocity(engine%grid, cfg_bad, engine%state, BALANCE_GEOSTROPHIC, ierr)
         call check(error, ierr, OCEAN_STATUS_ERR_IC_SEED, "wet/dry must be refused")
         if (allocated(error)) exit checks
         call untouched(error, engine, u0, v0, h0, "wetdry")
         if (allocated(error)) exit checks

         ! porous barriers / cavity: the latches the gate reads
         engine%state%metrics%use_porous = .true.
         call ocean_balance_velocity(engine%grid, cfg, engine%state, BALANCE_GEOSTROPHIC, ierr)
         engine%state%metrics%use_porous = .false.
         call check(error, ierr, OCEAN_STATUS_ERR_IC_SEED, "porous must be refused")
         if (allocated(error)) exit checks
         call untouched(error, engine, u0, v0, h0, "porous")
         if (allocated(error)) exit checks
         engine%state%metrics%use_cavity = .true.
         call ocean_balance_velocity(engine%grid, cfg, engine%state, BALANCE_GEOSTROPHIC, ierr)
         engine%state%metrics%use_cavity = .false.
         call check(error, ierr, OCEAN_STATUS_ERR_IC_SEED, "cavity must be refused")
         if (allocated(error)) exit checks
         call untouched(error, engine, u0, v0, h0, "cavity")
         if (allocated(error)) exit checks

         ! BALANCE_NONE is a successful no-op
         call ocean_balance_velocity(engine%grid, cfg, engine%state, BALANCE_NONE, ierr)
         call check(error, ierr, OCEAN_STATUS_OK, "BALANCE_NONE is a no-op")
         if (allocated(error)) exit checks
         call untouched(error, engine, u0, v0, h0, "none")
      end block checks

      call engine_exit_data(engine)
      call engine_teardown(engine)
   end subroutine test_refusals

   subroutine untouched(error, engine, u0, v0, h0, what)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), intent(inout) :: engine
      real(wp), intent(in) :: u0(:, :, :), v0(:, :, :), h0(:, :, :)
      character(len=*), intent(in) :: what
      call pull_uvh(engine)
      call check(error, all(engine%state%multilayer%u_face_x_layer == u0) .and. &
                 all(engine%state%multilayer%v_face_y_layer == v0) .and. &
                 all(engine%state%multilayer%h_layer == h0), &
                 what//": refusal must leave u, v, h untouched")
   end subroutine untouched

end module test_ocean_balance_velocity
