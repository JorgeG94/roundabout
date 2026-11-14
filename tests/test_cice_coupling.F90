!! Ocean-ice coupling library of the CICE-port ice (D11.2): the ocean import, the window-mean
!! flux averager (+ restart tags), the `Q_*` adapters (virtual-salt form), the open-fraction
!! hook and the face stress blend.
!!
!! "Closure B on a stub ice": a stub emits known, time-varying fluxes with a varying ice
!! step; over several thermo windows (ratio 1, 2, 3, both component modes) the heat, salt
!! and fresh water the OCEAN-side fields receive, times the window, equal the time integral
!! of what the stub emitted.  The residual of each line is printed (`closure_B ...`).
!! Device rules: every slot is mapped with `enter_data`, host-set inputs are pushed with
!! `!$acc update device`, outputs pulled with `update self` (component arrays only), and
!! every owned cell is compared (cells all differ).  Finiteness is tested explicitly.
module test_cice_coupling
   use, intrinsic :: iso_fortran_env, only: int32, int64, error_unit
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
#ifdef RDB_GPU_OFFLOAD
   use openacc, only: acc_present_dump
#endif
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, RHO_WATER
   use lcg_deterministic, only: lcg_next, lcg_unit
   use rdb_grid, only: hgrid_t
   use rdb_config, only: config_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_eos, only: eos_t, TS_CONS_ABS
   use rdb_teos10_convert, only: pt_from_ct, sp_from_sr
   use rdb_decomp, only: decomp_t, decomp_init
   use rdb_ocean_restart, only: restart_registry_t
   use rdb_ocean_restart_io, only: ocean_restart_write_local, ocean_restart_read_local, &
                                   ocean_restart_metadata_t
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t
   use rdb_ocean_metrics, only: metrics_fill_cartesian, metrics_finalize
   use rdb_ocean_state, only: ocean_state_t, ocean_state_enter_data, ocean_state_exit_data, &
                              ocean_state_seed_from_cfg, ocean_state_build_restart_registry
   use rdb_state, only: register_default_tracers
   use rdb_ice_state, only: ocean_sea_ice_t
   use rdb_ice_params, only: ice_params_t, ice_params_default, ICE_CONSTANTS_ICEPACK
   use rdb_ice_errors, only: ICE_ERR_NONFINITE
   use rdb_ice_cice_state, only: ice_cice_state_t
   use rdb_ice_flux_state, only: ice_flux_t
   use rdb_ice_flux_avg, only: ice_flux_avg_t, ice_flux_avg_accumulate, ice_flux_avg_reset, &
                               ICE_AVG_OK, ICE_AVG_ERR_DT, ICE_AVG_ERR_SHAPE
   use rdb_ice_import_ocean, only: ice_ocean_in_t, ice_import_ocean, ICE_IMPORT_OK
   use rdb_ice_melt_potential, only: ice_melt_potential
   use rdb_ocean_halo, only: ocean_halo_init, ocean_halo_destroy
   use rdb_comm_env, only: comm_env_init, comm_env_setup_roles
   use rdb_ocean_boundary_types, only: ocean_bc_state_t, ocean_bc_state_init, OBC_PERIODIC, &
                                       ocean_bc_validate_periodic
   use rdb_ocean_surface_flux, only: ocean_surface_flux_assemble, ocean_surface_flux_apply_tracers
   use rdb_ice_import_ocean, only: ice_melt_cap
   use rdb_ice_cice_adapters, only: ice_cice_stress_mirror_t, ice_cice_stress_resume_apply, &
                                    ICE_ADAPT_ERR_NO_COMPONENT, ICE_ADAPT_ERR_WINDOW, &
                                    ice_cice_to_surface_flux, ice_cice_apply_window, ice_cice_stress_blend, &
                                    ICE_ADAPT_OK, ICE_ADAPT_ERR_EMPTY, PSU_PER_MASS_FRACTION
   implicit none
   private

   public :: collect_cice_coupling_tests

   integer, parameter :: NXP = 7, NYP = 5, NG = 3, NCAT = 2, NILYR = 7
   integer, parameter :: NX = NXP + 2*NG, NY = NYP + 2*NG
   integer, parameter :: IS = NG + 1, IE = NG + NXP, JS = NG + 1, JE = NG + NYP
   real(wp), parameter :: SREF = 5.0_wp
   real(wp), parameter :: RHOI = 917.0_wp
   logical, save :: halo_up = .false.
   logical, save :: comm_up = .false.

contains

   subroutine collect_cice_coupling_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("avg_lifecycle_and_default_absent", test_lifecycle), &
                  new_unittest("avg_scheme_default_registers_nothing", test_sis2_untouched), &
                  new_unittest("avg_refuses_bad_dt_and_shape", test_avg_refuses), &
                  new_unittest("closure_B_stub_ice", test_closure_b), &
                  new_unittest("ratio1_bitwise_direct_handoff", test_ratio1_bitwise), &
                  new_unittest("avg_mid_window_restart_bitwise", test_restart), &
                  new_unittest("adapters_virtual_salt_closed_forms", test_adapters), &
                  new_unittest("adapters_update_ocn_f_agnostic", test_update_ocn_f), &
                  new_unittest("apply_window_refuses_empty", test_empty_window), &
                  new_unittest("import_closed_forms_vanished_top", test_import), &
                  new_unittest("import_cons_abs_conversion", test_import_cons), &
                  new_unittest("import_nan_in", test_import_nan), &
                  new_unittest("import_frzmlt_cadence", test_import_cadence), &
                  new_unittest("open_fraction_hook", test_open_fraction), &
                  new_unittest("stress_blend_hook", test_stress_blend), &
                  new_unittest("stress_blend_needs_unweighted_wind", test_stress_unweighted), &
                  new_unittest("stress_mirror_restart_roundtrip", test_stress_restart), &
                  new_unittest("ice_components_summed_not_overwritten", test_ice_components_end_to_end), &
                  new_unittest("seam_ghosts_equal_owners_periodic", test_ghosts), &
                  new_unittest("melt_cap_window_bound", test_melt_cap), &
                  new_unittest("melt_cap_ratio1_bitwise", test_melt_cap_ratio1), &
                  new_unittest("avg_required_tags_and_window_check", test_required_and_window), &
                  new_unittest("nothing_left_mapped", test_nothing_left) &
                  ]
   end subroutine collect_cice_coupling_tests

   ! ---------------------------------------------------------------- helpers

   function rnd(s) result(u)
      !! Next uniform number in [0, 1) from the LCG state `s`.
      integer(int64), intent(inout) :: s
      real(wp) :: u
      s = lcg_next(s)
      u = lcg_unit(s)
   end function rnd

   subroutine make_flux(g, fx, ierr)
      type(hgrid_t), intent(out) :: g
      type(ice_flux_t), intent(inout) :: fx
      integer, intent(out) :: ierr
      call g%init(NXP, NYP, NG, 1000.0_wp, 1000.0_wp)
      call fx%init(g, NCAT, ierr)
   end subroutine make_flux

   subroutine emit(fx, base, step)
      !! Stub ice: every cell and variable has its own base value, modulated by a
      !! step-dependent triangle; ghosts get a large sentinel the averager must ignore.
      type(ice_flux_t), intent(inout) :: fx
      real(wp), intent(in) :: base(:, :, :)
      integer, intent(in) :: step
      real(wp) :: m
      integer :: i, j
      m = 1.0_wp + 0.25_wp*real(modulo(step*7, 11) - 5, wp)/5.0_wp
      do j = 1, NY
         do i = 1, NX
            if (i >= IS .and. i <= IE .and. j >= JS .and. j <= JE) then
               fx%fresh(i, j) = base(i, j, 1)*m*1.0e-5_wp
               fx%fsalt(i, j) = base(i, j, 2)*m*1.0e-7_wp
               fx%fhocn(i, j) = base(i, j, 3)*m*30.0_wp
               fx%fswthru(i, j) = base(i, j, 4)*m*10.0_wp
               fx%fswthru_vis(i, j) = base(i, j, 4)*m*10.0_wp
            else
               fx%fresh(i, j) = 1.0e9_wp
               fx%fsalt(i, j) = 1.0e9_wp
               fx%fhocn(i, j) = 1.0e9_wp
               fx%fswthru(i, j) = 1.0e9_wp
               fx%fswthru_vis(i, j) = 1.0e9_wp
            end if
         end do
      end do
      !$acc update device(fx%fresh, fx%fsalt, fx%fhocn, fx%fswthru, fx%fswthru_vis)
   end subroutine emit

   subroutine make_base(base, seed)
      real(wp), intent(out) :: base(:, :, :)
      integer(int64), intent(in) :: seed
      integer(int64) :: s
      integer :: i, j, v
      s = seed
      do v = 1, 4
         do j = 1, NY
            do i = 1, NX
               base(i, j, v) = 2.0_wp*(rnd(s) - 0.5_wp)
            end do
         end do
      end do
   end subroutine make_base

   subroutine make_halo(grid, bc, periodic)
      !! Boundary state + single-rank halo: all four edges periodic when `periodic`.
      type(hgrid_t), intent(in) :: grid
      type(ocean_bc_state_t), intent(out) :: bc
      logical, intent(in) :: periodic
      type(decomp_t) :: d

      call ocean_bc_state_init(bc, grid, nz_ml=1, n_tracers=0)
      if (periodic) then
         bc%west%bc_type = OBC_PERIODIC
         bc%east%bc_type = OBC_PERIODIC
         bc%south%bc_type = OBC_PERIODIC
         bc%north%bc_type = OBC_PERIODIC
         call ocean_bc_validate_periodic(bc)
      end if
      if (.not. comm_up) then
         call comm_env_init()
         call comm_env_setup_roles(.false.)
         comm_up = .true.
      end if
      if (halo_up) call ocean_halo_destroy()
      call decomp_init(d, grid%nx_phys, grid%ny_phys, 1, 1, 0)
      call ocean_halo_init(d, grid%nghost, bc%periodic_x, bc%periodic_y)
      halo_up = .true.
   end subroutine make_halo

   function owned_max(a) result(m)
      real(wp), intent(in) :: a(:, :)
      real(wp) :: m
      m = maxval(abs(a(IS:IE, JS:JE)))
   end function owned_max

   ! ------------------------------------------------------------ lifecycle

   subroutine test_lifecycle(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: g
      type(ice_flux_avg_t) :: av, none
      type(ice_ocean_in_t) :: oin, none_in
      type(restart_registry_t) :: reg
      integer(int64) :: b0

      call check(error, none%bytes() == 0_int64 .and. none_in%bytes() == 0_int64 .and. &
                 .not. allocated(none%acc_fresh) .and. .not. allocated(none_in%sst), &
                 "a default slot holds nothing (no arrays, 0 bytes)")
      call none%register_restart(reg)
      call check(error, reg%n == 0, "an uninitialised averager registers no restart tag")
      if (allocated(error)) return
      call g%init(NXP, NYP, NG, 1000.0_wp, 1000.0_wp)
      call av%init(g)
      call oin%init(g)
      b0 = 5_int64*int(NX*NY, int64)*storage_size(1.0_wp)/8
      call check(error, av%bytes() == b0 .and. oin%bytes() == b0, "bytes() counts the five arrays of each slot")
      call av%enter_data()
      call oin%enter_data()
      call check(error, av%mapped .and. oin%mapped, "mapped")
      call av%register_restart(reg)
      call check(error, reg%n == 7, "five arrays + two window scalars registered")
      call check(error, all(reg%entries(1:7)%optional), "every averager tag is optional")
      call check(error, trim(reg%entries(1)%tag) == "ice_avg_fresh" .and. trim(reg%entries(7)%tag) == "ice_avg_t", &
                 "tags")
      call av%destroy()
      call oin%destroy()
      call check(error, av%bytes() == 0_int64 .and. .not. av%is_init .and. .not. allocated(av%acc_fresh), "destroy")
   end subroutine test_lifecycle

   subroutine test_sis2_untouched(error)
      !! The default scheme: the parent ocean state composes none of the new slots, so its
      !! restart registry holds no `ice_avg_*` tag.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      type(restart_registry_t) :: reg
      integer :: e

      cfg%sim_type = "ocean"
      cfg%nx = NXP
      cfg%ny = NYP
      cfg%dx = 10.0e3_wp
      cfg%dy = 10.0e3_wp
      cfg%nz_layers = 2
      cfg%nghost = NG
      cfg%ocean%topo%max_depth = 1000.0_wp
      cfg%initial_temperature = 10.0_wp
      cfg%initial_salinity = 34.0_wp
      call grid%init(NXP, NYP, NG, 10.0e3_wp, 10.0e3_wp)
      state%multilayer%nz_ml = 2
      state%ice%enable = .true.
      state%ice%ncat = 2
      state%ice%nk_ice = 2
      call state%init(grid)
      call metrics_fill_cartesian(state%metrics, grid, grid%dx, grid%dy)
      call metrics_finalize(state%metrics)
      call ocean_state_seed_from_cfg(state, grid, cfg)
      call register_default_tracers( &
         state%multilayer%tracers(state%multilayer%idx_salinity), &
         state%multilayer%tracers(state%multilayer%idx_temperature), cfg)
      call ocean_state_enter_data(state)
      call ocean_state_build_restart_registry(state, grid, reg)
      do e = 1, reg%n
         call check(error, index(reg%entries(e)%tag, "ice_avg_") == 0, &
                    "no averager tag registered: "//trim(reg%entries(e)%tag))
         if (allocated(error)) exit
      end do
      call ocean_state_exit_data(state)
      call state%destroy()
   end subroutine test_sis2_untouched

   subroutine test_avg_refuses(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: g, g2
      type(ice_flux_t) :: fx, fx2
      type(ice_flux_avg_t) :: av
      integer :: ierr

      call make_flux(g, fx, ierr)
      call av%init(g)
      call fx%enter_data()
      call av%enter_data()
      call ice_flux_avg_accumulate(av, fx, 0.0_wp, ierr)
      call check(error, ierr == ICE_AVG_ERR_DT, "dt = 0 refused")
      call ice_flux_avg_accumulate(av, fx, ieee_value(1.0_wp, ieee_quiet_nan), ierr)
      call check(error, ierr == ICE_AVG_ERR_DT, "dt = NaN refused")
      call ice_flux_avg_accumulate(av, fx, -1.0_wp, ierr)
      call check(error, ierr == ICE_AVG_ERR_DT, "dt < 0 refused")
      call check(error, av%n_acc == 0.0_wp .and. av%t_acc == 0.0_wp, "a refused call leaves the window untouched")
      call g2%init(NXP + 1, NYP, NG, 1000.0_wp, 1000.0_wp)
      call fx2%init(g2, NCAT, ierr)
      call ice_flux_avg_accumulate(av, fx2, 900.0_wp, ierr)
      call check(error, ierr == ICE_AVG_ERR_SHAPE, "shape mismatch refused")
      call av%exit_data()
      call fx%exit_data()
   end subroutine test_avg_refuses

   ! ---------------------------------------------------- closure B (stub ice)

   subroutine test_closure_b(error)
      !! Stub ice -> averager -> adapters -> ocean surface-flux fields, ratio 1, 2, 3, both
      !! component modes.  Per cell, over `NWIN` windows, with ice steps of 900, 1000, 1100 s:
      !!   heat : sum_win (Q_heat - Q_const)*window       = sum_steps (fhocn + fswthru)*dt
      !!   salt : sum_win (Q_salt - Q_const)*window       = sum_steps (1000 fsalt - S fresh)*dt
      !!   water: sum_win fresh_mass_out*window           = sum_steps fresh*dt
      !! (S constant in time, so the hand-off S is the S of every step; components on: heat
      !! = heat_added + q_sw).  The salt constant is 0 here: subtracting a constant ~1e3 times
      !! the flux would put the cancellation error of the TEST (3e-13) into the residual; the
      !! constants are pinned in `test_adapters`.  Tolerance 5e-15 of the sum of absolute
      !! contributions (measured worst 3.9e-16, printed per case).
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NWIN = 4
      integer :: ratio, mode
      logical :: comps

      do ratio = 1, 3
         do mode = 1, 2
            comps = mode == 2
            call closure_case(ratio, NWIN, comps, error)
            if (allocated(error)) return
         end do
      end do
   end subroutine test_closure_b

   subroutine closure_case(ratio, nwin, comps, error)
      integer, intent(in) :: ratio, nwin
      logical, intent(in) :: comps
      type(error_type), allocatable, intent(inout) :: error
      type(hgrid_t) :: g
      type(ice_flux_t) :: fx
      type(ice_flux_avg_t) :: av
      type(ocean_surface_flux_t) :: sf
      type(ocean_bc_state_t) :: bc
      real(wp), allocatable :: base(:, :, :), sss(:, :), fw(:, :)
      real(wp), allocatable :: i_h(:, :), i_s(:, :), i_w(:, :), o_h(:, :), o_s(:, :), o_w(:, :)
      real(wp), allocatable :: a_h(:, :), a_s(:, :), a_w(:, :)
      real(wp) :: dt, win, qhc, qsc, r_h, r_s, r_w, tolr
      integer :: ierr, step, w, k, i, j
      integer(int64) :: s

      call make_flux(g, fx, ierr)
      call av%init(g)
      call sf%init(g)
      if (comps) call sf%set_components(g, .true.)
      if (comps) call sf%set_ice_components(g, .true.)
      call make_halo(g, bc, .false.)
      qhc = 3.0_wp
      qsc = 0.0_wp
      sf%Q_heat_const = qhc
      sf%Q_salt_const = qsc
      allocate (base(NX, NY, 4), sss(NX, NY), fw(NX, NY))
      allocate (i_h(NX, NY), i_s(NX, NY), i_w(NX, NY), o_h(NX, NY), o_s(NX, NY), o_w(NX, NY))
      allocate (a_h(NX, NY), a_s(NX, NY), a_w(NX, NY))
      call make_base(base, 424242_int64 + 7_int64*ratio)
      s = 99_int64
      do j = 1, NY
         do i = 1, NX
            sss(i, j) = 30.0_wp + 6.0_wp*rnd(s)
         end do
      end do
      fw = 0.0_wp
      i_h = 0.0_wp
      i_s = 0.0_wp
      i_w = 0.0_wp
      o_h = 0.0_wp
      o_s = 0.0_wp
      o_w = 0.0_wp
      a_h = 0.0_wp
      a_s = 0.0_wp
      a_w = 0.0_wp
      !$acc enter data copyin(sss, fw)
      call fx%enter_data()
      call av%enter_data()
      call sf%enter_data()
      step = 0
      do w = 1, nwin
         win = 0.0_wp
         do k = 1, ratio
            step = step + 1
            dt = 900.0_wp + 100.0_wp*real(modulo(step, 3), wp)
            call emit(fx, base, step)
            call ice_flux_avg_accumulate(av, fx, dt, ierr)
            call check(error, ierr, ICE_AVG_OK, "accumulate")
            if (allocated(error)) goto 900
            win = win + dt
            i_h = i_h + (fx%fhocn + fx%fswthru)*dt
            i_s = i_s + (PSU_PER_MASS_FRACTION*fx%fsalt - sss*fx%fresh)*dt
            i_w = i_w + fx%fresh*dt
            a_h = a_h + (abs(fx%fhocn) + abs(fx%fswthru))*dt
            a_s = a_s + (abs(PSU_PER_MASS_FRACTION*fx%fsalt) + abs(sss*fx%fresh))*dt
            a_w = a_w + abs(fx%fresh)*dt
         end do
         call ice_cice_apply_window(g, bc, av, sf, sss, ierr, fresh_mass_out=fw, window_expected=win)
         call check(error, ierr, ICE_ADAPT_OK, "apply_window")
         if (allocated(error)) goto 900
         call check(error, av%n_acc == 0.0_wp .and. av%t_acc == 0.0_wp .and. owned_max(av%acc_fhocn) == 0.0_wp, &
                    "the window is reset after the hand-off")
         if (comps) then
            !$acc update self(sf%ice_heat, sf%ice_sw, sf%ice_salt, fw)
            o_h = o_h + (sf%ice_heat + sf%ice_sw)*win
            o_s = o_s + sf%ice_salt*win
         else
            !$acc update self(sf%Q_heat, sf%Q_salt, fw)
            o_h = o_h + (sf%Q_heat - qhc)*win
            o_s = o_s + (sf%Q_salt - qsc)*win
         end if
         o_w = o_w + fw*win
      end do
      r_h = maxval(abs(o_h(IS:IE, JS:JE) - i_h(IS:IE, JS:JE))/a_h(IS:IE, JS:JE))
      r_s = maxval(abs(o_s(IS:IE, JS:JE) - i_s(IS:IE, JS:JE))/a_s(IS:IE, JS:JE))
      r_w = maxval(abs(o_w(IS:IE, JS:JE) - i_w(IS:IE, JS:JE))/a_w(IS:IE, JS:JE))
      write (error_unit, "(a,i0,a,l1,a,3es10.2)") "closure_B ratio=", ratio, " components=", comps, &
         "  max rel residual heat/salt/water =", r_h, r_s, r_w
      tolr = 5.0e-15_wp
      call check(error, ieee_is_finite(r_h) .and. ieee_is_finite(r_s) .and. ieee_is_finite(r_w), &
                 "closure residuals are finite")
      if (allocated(error)) goto 900
      call check(error, r_h <= tolr, "closure B heat: ocean-applied == ice-emitted")
      call check(error, r_s <= tolr, "closure B salt (virtual form): ocean-applied == ice-emitted")
      call check(error, r_w <= tolr, "closure B water: mass flux handed on == ice-emitted")
      if (allocated(error)) goto 900
900   continue
      call sf%exit_data()
      call av%exit_data()
      call fx%exit_data()
      !$acc exit data delete(sss, fw)
   end subroutine closure_case

   subroutine test_ratio1_bitwise(error)
      !! Ratio 1: the averaged hand-off is bit for bit the direct hand-off of the same step's
      !! fluxes (adapters on `fx%*` versus accumulate + apply_window), over several steps
      !! with different dt; `==` on every owned cell of every output.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: g
      type(ice_flux_t) :: fx
      type(ice_flux_avg_t) :: av
      type(ocean_surface_flux_t) :: sa, sb
      type(ocean_bc_state_t) :: bc
      real(wp), allocatable :: base(:, :, :), sss(:, :), fw(:, :)
      integer :: ierr, step
      real(wp) :: dt

      call make_flux(g, fx, ierr)
      call av%init(g)
      call sa%init(g)
      call sb%init(g)
      call make_halo(g, bc, .false.)
      sa%Q_heat_const = 3.0_wp
      sb%Q_heat_const = 3.0_wp
      allocate (base(NX, NY, 4), sss(NX, NY), fw(NX, NY))
      call make_base(base, 5150_int64)
      sss = 33.3_wp
      fw = 0.0_wp
      !$acc enter data copyin(sss, fw)
      call fx%enter_data()
      call av%enter_data()
      call sa%enter_data()
      call sb%enter_data()
      do step = 1, 5
         dt = 600.0_wp + 137.0_wp*real(step, wp)
         call emit(fx, base, step)
         call ice_cice_to_surface_flux(g, sa, fx%fresh, fx%fsalt, fx%fhocn, fx%fswthru, sss, ierr)
         call check(error, ierr, ICE_ADAPT_OK, "direct")
         call ice_flux_avg_accumulate(av, fx, dt, ierr)
         call ice_cice_apply_window(g, bc, av, sb, sss, ierr, fresh_mass_out=fw)
         call check(error, ierr, ICE_ADAPT_OK, "averaged")
         if (allocated(error)) goto 900
         !$acc update self(sa%Q_heat, sa%Q_salt, sb%Q_heat, sb%Q_salt, fw)
         call check(error, all(sa%Q_heat(IS:IE, JS:JE) == sb%Q_heat(IS:IE, JS:JE)), &
                    "Q_heat bitwise identical to the direct hand-off")
         call check(error, all(sa%Q_salt(IS:IE, JS:JE) == sb%Q_salt(IS:IE, JS:JE)), &
                    "Q_salt bitwise identical to the direct hand-off")
         call check(error, all(fw(IS:IE, JS:JE) == fx%fresh(IS:IE, JS:JE)), "fresh bitwise identical")
         if (allocated(error)) goto 900
      end do
900   continue
      call sb%exit_data()
      call sa%exit_data()
      call av%exit_data()
      call fx%exit_data()
      !$acc exit data delete(sss, fw)
   end subroutine test_ratio1_bitwise

   subroutine test_restart(error)
      !! A checkpoint in the middle of a 3-step window: write the averager through the registry
      !! after step 2, read into a second averager, finish step 3 on both: the means and the
      !! window counters are bitwise equal, ghosts included.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: g
      type(ice_flux_t) :: fx
      type(ice_flux_avg_t) :: a, b
      type(restart_registry_t) :: ra, rb
      type(decomp_t) :: d
      type(ocean_restart_metadata_t) :: meta
      character(len=*), parameter :: FN = "test_cice_coupling_restart.nc"
      real(wp), allocatable :: base(:, :, :)
      real(wp) :: t
      integer :: ierr, step, osc, lun, rstep
      real(wp) :: dts(3)

      dts = [900.0_wp, 1000.0_wp, 1100.0_wp]
      call make_flux(g, fx, ierr)
      call a%init(g)
      call b%init(g)
      allocate (base(NX, NY, 4))
      call make_base(base, 777_int64)
      call fx%enter_data()
      call a%enter_data()
      do step = 1, 2
         call emit(fx, base, step)
         call ice_flux_avg_accumulate(a, fx, dts(step), ierr)
      end do
      ! The registry writer is host-only: like the state's own restart path, the caller pulls
      ! the device-mapped entries first.
      !$acc update self(a%acc_fresh, a%acc_fsalt, a%acc_fhocn, a%acc_fswthru, a%acc_fswthru_vis)
      call a%register_restart(ra)
      call b%register_restart(rb)
      call decomp_init(d, NXP, NYP, 1, 1, 0)
      meta%nghost = NG
      open (newunit=lun, file=FN, status="replace")
      close (lun, status="delete")
      call ocean_restart_write_local(FN, ra, d, meta, 3.5_wp, 7, 2, ierr)
      call check(error, ierr, 0, "write ok")
      call ocean_restart_read_local(FN, rb, d, meta, t, rstep, osc, ierr)
      call check(error, ierr, 0, "read ok")
      open (newunit=lun, file=FN, status="old")
      close (lun, status="delete")
      if (allocated(error)) goto 900
      call check(error, all(rb%entries(1:7)%found), "every averager entry found in the file")
      call check(error, b%n_acc == 2.0_wp .and. b%t_acc == dts(1) + dts(2), "window counters restored")
      if (allocated(error)) goto 900
      call b%enter_data()
      call emit(fx, base, 3)
      call ice_flux_avg_accumulate(a, fx, dts(3), ierr)
      call ice_flux_avg_accumulate(b, fx, dts(3), ierr)
      !$acc update self(a%acc_fresh, a%acc_fsalt, a%acc_fhocn, a%acc_fswthru, a%acc_fswthru_vis)
      !$acc update self(b%acc_fresh, b%acc_fsalt, b%acc_fhocn, b%acc_fswthru, b%acc_fswthru_vis)
      call check(error, all(a%acc_fresh == b%acc_fresh) .and. all(a%acc_fsalt == b%acc_fsalt) .and. &
                 all(a%acc_fhocn == b%acc_fhocn) .and. all(a%acc_fswthru == b%acc_fswthru) .and. &
                 all(a%acc_fswthru_vis == b%acc_fswthru_vis), "resumed window mean is bitwise the uninterrupted one")
      call check(error, a%n_acc == b%n_acc .and. a%t_acc == b%t_acc, "counters equal after the resume")
900   continue
      call a%exit_data()
      call b%exit_data()
      call fx%exit_data()
   end subroutine test_restart

   ! -------------------------------------------------------------- adapters

   subroutine adapt_one(comps, fresh, fsalt, fhocn, fswthru, s, qh, qs, qsw, hadd, saltf)
      !! One 7x5 grid with the same flux everywhere, run through the shim; returns an owned
      !! cell of each written field.
      logical, intent(in) :: comps
      real(wp), intent(in) :: fresh, fsalt, fhocn, fswthru, s
      real(wp), intent(out) :: qh, qs, qsw, hadd, saltf
      type(hgrid_t) :: g
      type(ocean_surface_flux_t) :: sf
      real(wp), allocatable :: a(:, :), b(:, :), c(:, :), d(:, :), e(:, :)
      integer :: ierr

      call g%init(NXP, NYP, NG, 1000.0_wp, 1000.0_wp)
      call sf%init(g)
      if (comps) call sf%set_components(g, .true.)
      if (comps) call sf%set_ice_components(g, .true.)
      sf%Q_heat_const = 3.0_wp
      sf%Q_salt_const = -0.5_wp
      allocate (a(NX, NY), b(NX, NY), c(NX, NY), d(NX, NY), e(NX, NY))
      a = fresh
      b = fsalt
      c = fhocn
      d = fswthru
      e = s
      !$acc enter data copyin(a, b, c, d, e)
      call sf%enter_data()
      call ice_cice_to_surface_flux(g, sf, a, b, c, d, e, ierr)
      qh = -1.0e30_wp
      qs = -1.0e30_wp
      qsw = -1.0e30_wp
      hadd = -1.0e30_wp
      saltf = -1.0e30_wp
      if (comps) then
         !$acc update self(sf%ice_heat, sf%ice_sw, sf%ice_salt)
         hadd = sf%ice_heat(IS + 1, JS + 1)
         qsw = sf%ice_sw(IS + 1, JS + 1)
         saltf = sf%ice_salt(IS + 1, JS + 1)
         qh = merge(1.0_wp, 0.0_wp, sf%has_heat .and. sf%has_salt)
         if (any(sf%ice_heat(1:NG, :) /= 0.0_wp)) qh = -1.0_wp
      else
         !$acc update self(sf%Q_heat, sf%Q_salt)
         qh = sf%Q_heat(IS + 1, JS + 1)
         qs = sf%Q_salt(IS + 1, JS + 1)
         if (any(sf%Q_heat(1:NG, :) /= 0.0_wp) .or. any(sf%Q_salt(:, NY - NG + 1:NY) /= 0.0_wp)) qs = -2.0e30_wp
      end if
      !$acc exit data delete(a, b, c, d, e)
      call sf%exit_data()
   end subroutine adapt_one

   subroutine test_adapters(error)
      !! Closed forms of the virtual-salt form: ice growth of `dhi` m in `dt` s (new ice with
      !! `update_ocn_f`): fresh = -rhoi dhi/dt, fsalt = -Sref 1e-3 rhoi dhi/dt  =>
      !! Q_salt = rhoi dhi (S - Sref)/dt > 0 (S > Sref: brine rejection salinifies); melt flips.
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: DT = 900.0_wp, S = 34.0_wp, DHI = 0.02_wp
      real(wp) :: qh, qs, qsw, hadd, saltf, fresh, fsalt, want

      fresh = -RHOI*DHI/DT
      fsalt = -SREF*1.0e-3_wp*RHOI*DHI/DT
      call adapt_one(.false., fresh, fsalt, 12.5_wp, 7.25_wp, S, qh, qs, qsw, hadd, saltf)
      want = RHOI*DHI*(S - SREF)/DT
      call check(error, qs - (-0.5_wp) > 0.0_wp, "growth salinifies (Q_salt - const > 0)")
      call check(error, abs((qs + 0.5_wp) - want) <= 1.0e-14_wp*want, "growth: Q_salt = rhoi dhi (S - Sref)/dt")
      call check(error, qh == (3.0_wp + 12.5_wp) + 7.25_wp, "Q_heat = (const + fhocn) + fswthru, same operand order")
      if (allocated(error)) return
      call adapt_one(.false., -fresh, -fsalt, -4.0_wp, 0.0_wp, S, qh, qs, qsw, hadd, saltf)
      call check(error, qs - (-0.5_wp) < 0.0_wp, "melt freshens (Q_salt - const < 0)")
      call check(error, abs((qs + 0.5_wp) + want) <= 1.0e-14_wp*want, "melt: Q_salt = -rhoi dhi (S - Sref)/dt")
      call check(error, qh == (3.0_wp - 4.0_wp) + 0.0_wp, "melt heat: fhocn < 0 cools")
      call check(error, qs > -1.0e30_wp, "ghost cells of the ocean fields are untouched")
      if (allocated(error)) return
      call adapt_one(.true., fresh, fsalt, 12.5_wp, 7.25_wp, S, qh, qs, qsw, hadd, saltf)
      call check(error, hadd == 12.5_wp .and. qsw == 7.25_wp, "components on: ice_heat = fhocn, ice_sw = fswthru")
      call check(error, abs(saltf - want) <= 1.0e-14_wp*want, "components on: ice_salt has no const term")
      call check(error, qh == 1.0_wp, "components on: has_heat, has_salt latched; ghosts untouched")
   end subroutine test_adapters

   subroutine test_update_ocn_f(error)
      !! `update_ocn_f` on versus off emits the frazil term `-rhoi V/dt` into fresh and
      !! `-Sref 1e-3 rhoi V/dt` into fsalt, or nothing.  The adapters are agnostic: the ocean
      !! Q_salt differs between the two by exactly `rhoi V (S - Sref)/dt` (the water and salt
      !! the ice gained and the ocean did NOT lose when the knob is off), and the heat is the
      !! same.
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: DT = 900.0_wp, S = 33.0_wp, V = 0.03_wp
      real(wp) :: qh1, qs1, qh2, qs2, x1, x2, x3, x4, f0, s0, dfr, dfs

      f0 = 2.0e-6_wp
      s0 = 1.0e-9_wp
      dfr = -RHOI*V/DT
      dfs = -SREF*1.0e-3_wp*RHOI*V/DT
      call adapt_one(.false., f0 + dfr, s0 + dfs, 5.0_wp, 1.0_wp, S, qh1, qs1, x1, x2, x3)
      call adapt_one(.false., f0, s0, 5.0_wp, 1.0_wp, S, qh2, qs2, x1, x2, x4)
      call check(error, qh1 == qh2, "heat does not depend on update_ocn_f")
      call check(error, abs((qs1 - qs2) - RHOI*V*(S - SREF)/DT) <= 1.0e-12_wp*RHOI*V*(S - SREF)/DT, &
                 "salt differs by the unpaid frazil water and salt")
   end subroutine test_update_ocn_f

   subroutine test_empty_window(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: g
      type(ice_flux_avg_t) :: av
      type(ocean_surface_flux_t) :: sf
      type(ocean_bc_state_t) :: bc
      real(wp), allocatable :: sss(:, :)
      integer :: ierr

      call g%init(NXP, NYP, NG, 1000.0_wp, 1000.0_wp)
      call make_halo(g, bc, .false.)
      call av%init(g)
      call sf%init(g)
      sf%Q_heat = 11.0_wp
      allocate (sss(NX, NY))
      sss = 35.0_wp
      !$acc enter data copyin(sss)
      call av%enter_data()
      call sf%enter_data()
      call ice_cice_apply_window(g, bc, av, sf, sss, ierr)
      call check(error, ierr, ICE_ADAPT_ERR_EMPTY, "an empty window is refused")
      !$acc update self(sf%Q_heat)
      call check(error, all(sf%Q_heat == 11.0_wp), "and the ocean field is untouched")
      !$acc exit data delete(sss)
      call av%exit_data()
      call sf%exit_data()
   end subroutine test_empty_window

   ! ---------------------------------------------------------------- import

   subroutine import_setup(grid, ms, eos, cice, oin, frz, nz)
      type(hgrid_t), intent(out) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(eos_t), intent(out) :: eos
      type(ice_cice_state_t), intent(inout) :: cice
      type(ice_ocean_in_t), intent(inout) :: oin
      real(wp), allocatable, intent(out) :: frz(:, :)
      integer, intent(in) :: nz
      type(ice_params_t) :: p
      integer :: ierr, k, i, j
      integer(int64) :: s

      call grid%init(NXP, NYP, NG, 1000.0_wp, 1000.0_wp)
      ms%nz_ml = nz
      call ms%init(grid)
      call eos%init(grid)
      p = ice_params_default(ICE_CONSTANTS_ICEPACK, NCAT, NILYR, 1)
      call cice%init(grid, p, ierr)
      call oin%init(grid)
      allocate (frz(NX, NY))
      s = 31337_int64
      do j = 1, NY
         do i = 1, NX
            frz(i, j) = 200.0_wp*(rnd(s) - 0.5_wp)
            do k = 1, nz
               ms%h_layer(i, j, k) = 2.0_wp + 8.0_wp*rnd(s)
               ms%tracers(ms%idx_temperature)%c(i, j, k) = -1.5_wp + 6.0_wp*rnd(s)
               ms%tracers(ms%idx_salinity)%c(i, j, k) = 30.0_wp + 6.0_wp*rnd(s)
            end do
         end do
      end do
   end subroutine import_setup

   subroutine import_map(ms, cice, oin, frz)
      type(multilayer_state_t), intent(inout) :: ms
      type(ice_cice_state_t), intent(inout) :: cice
      type(ice_ocean_in_t), intent(inout) :: oin
      real(wp), intent(inout) :: frz(:, :)
      !$acc enter data copyin(ms)
      call ms%enter_data()
      call cice%enter_data()
      call oin%enter_data()
      !$acc enter data copyin(frz)
      ! A `copyin` of an array that a previous test left present keeps the device copy, so
      ! the host-set inputs are pushed explicitly (CLAUDE.md, "Writing GPU tests").
      !$acc update device(frz)
      !$acc update device(ms%wet_mask, ms%h_layer, ms%k_top)
      !$acc update device(ms%tracers(ms%idx_temperature)%c, ms%tracers(ms%idx_salinity)%c)
   end subroutine import_map

   subroutine import_unmap(ms, cice, oin, frz)
      type(multilayer_state_t), intent(inout) :: ms
      type(ice_cice_state_t), intent(inout) :: cice
      type(ice_ocean_in_t), intent(inout) :: oin
      real(wp), intent(inout) :: frz(:, :)
      !$acc exit data delete(frz)
      call oin%exit_data()
      call cice%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms)
   end subroutine import_unmap

   subroutine test_import(error)
      !! Every owned cell against an independent host reference, on a random state where the
      !! top layer is vanished in ~25 % of the columns and some columns are land: sst, sss,
      !! sss_ocn, Tf = -depressT sss and frzmlt come from the FIRST LIVE layer; land gives 0.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 4
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ice_cice_state_t) :: cice
      type(ice_ocean_in_t) :: oin
      real(wp), allocatable :: frz(:, :)
      integer :: i, j, k, kk, ierr, nbad, nvan, nland
      integer(int64) :: s
      real(wp) :: r, t, sa, dep

      call import_setup(grid, ms, eos, cice, oin, frz, NZ)
      s = 4711_int64
      nvan = 0
      nland = 0
      do j = 1, NY
         do i = 1, NX
            r = rnd(s)
            if (r < 0.25_wp) ms%h_layer(i, j, NZ) = 0.0_wp
            r = rnd(s)
            if (r < 0.2_wp) then
               ms%wet_mask(i, j) = 0.0_wp
            end if
         end do
      end do
      dep = cice%p%c%depressT
      call check(error, abs(dep - 0.054_wp) < 1.0e-15_wp, "depressT is the CICE 0.054")
      call import_map(ms, cice, oin, frz)
      call ice_import_ocean(grid, eos, ms, frz, cice, oin, ierr)
      call check(error, ierr, ICE_IMPORT_OK, "import ran")
      !$acc update self(oin%sst, oin%sss, oin%sss_ocn, oin%Tf, oin%frzmlt, cice%ice_err)
      nbad = 0
      do j = JS, JE
         do i = IS, IE
            if (ms%wet_mask(i, j) <= 0.5_wp) then
               nland = nland + 1
               if (oin%sst(i, j) /= 0.0_wp .or. oin%sss(i, j) /= 0.0_wp .or. oin%Tf(i, j) /= 0.0_wp .or. &
                   oin%frzmlt(i, j) /= 0.0_wp .or. oin%sss_ocn(i, j) /= 0.0_wp) nbad = nbad + 1
               cycle
            end if
            kk = NZ
            do k = NZ, 1, -1
               if (ms%h_layer(i, j, k) > 1.5e-4_wp) then
                  kk = k
                  exit
               end if
            end do
            if (kk /= NZ) nvan = nvan + 1
            t = ms%tracers(ms%idx_temperature)%c(i, j, kk)
            sa = ms%tracers(ms%idx_salinity)%c(i, j, kk)
            if (.not. all(ieee_is_finite([oin%sst(i, j), oin%sss(i, j), oin%Tf(i, j), oin%frzmlt(i, j)]))) then
               nbad = nbad + 1
            else if (oin%sst(i, j) /= t .or. oin%sss(i, j) /= sa .or. oin%sss_ocn(i, j) /= sa .or. &
                     abs(oin%Tf(i, j) + dep*sa) > 1.0e-14_wp*abs(dep*sa) .or. oin%frzmlt(i, j) /= frz(i, j)) then
               nbad = nbad + 1
            end if
            if (cice%ice_err(i, j) /= 0_int32) nbad = nbad + 1
         end do
      end do
      call check(error, nbad == 0, "every owned cell matches the host reference")
      call check(error, nvan > 3 .and. nland > 3, "the random state exercises vanished tops and land")
      call check(error, all(oin%sst(1:NG, :) == 0.0_wp) .and. all(oin%Tf(:, NY - NG + 1:NY) == 0.0_wp), &
                 "ghosts are not written")
      call import_unmap(ms, cice, oin, frz)
   end subroutine test_import

   subroutine test_import_cons(error)
      !! `cons_abs`: sst/sss are converted to potential temperature / practical salinity (the
      !! ice-boundary convention), sss_ocn stays as the ocean carries it, Tf follows the
      !! converted salinity.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ice_cice_state_t) :: cice
      type(ice_ocean_in_t) :: oin
      real(wp), allocatable :: frz(:, :)
      integer :: i, j, ierr, nbad
      real(wp) :: t, sa

      call import_setup(grid, ms, eos, cice, oin, frz, 3)
      eos%ts_convention = TS_CONS_ABS
      call import_map(ms, cice, oin, frz)
      call ice_import_ocean(grid, eos, ms, frz, cice, oin, ierr)
      !$acc update self(oin%sst, oin%sss, oin%sss_ocn, oin%Tf)
      nbad = 0
      do j = JS, JE
         do i = IS, IE
            t = ms%tracers(ms%idx_temperature)%c(i, j, 3)
            sa = ms%tracers(ms%idx_salinity)%c(i, j, 3)
            if (abs(oin%sst(i, j) - pt_from_ct(sa, t)) > 1.0e-10_wp) nbad = nbad + 1
            if (abs(oin%sss(i, j) - sp_from_sr(sa)) > 1.0e-13_wp*sa) nbad = nbad + 1
            if (oin%sss_ocn(i, j) /= sa) nbad = nbad + 1
            if (abs(oin%Tf(i, j) + 0.054_wp*oin%sss(i, j)) > 1.0e-13_wp) nbad = nbad + 1
         end do
      end do
      call check(error, nbad == 0, "cons_abs conversion on every owned cell")
      call import_unmap(ms, cice, oin, frz)
   end subroutine test_import_cons

   subroutine test_import_nan(error)
      !! A NaN in the surface T, S, thickness or the frzmlt source of a wet cell sets the
      !! cell's ice error bit and NaNs all five outputs; the neighbours are unaffected.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ice_cice_state_t) :: cice
      type(ice_ocean_in_t) :: oin
      real(wp), allocatable :: frz(:, :)
      integer :: ierr, c, i0, j0, nbad
      real(wp) :: nan
      logical :: nonfin, bit

      nan = ieee_value(1.0_wp, ieee_quiet_nan)
      do c = 1, 4
         call import_setup(grid, ms, eos, cice, oin, frz, 3)
         i0 = IS + 2
         j0 = JS + 1
         select case (c)
         case (1)
            ms%tracers(ms%idx_temperature)%c(i0, j0, 3) = nan
         case (2)
            ms%tracers(ms%idx_salinity)%c(i0, j0, 3) = nan
         case (3)
            ms%h_layer(i0, j0, 3) = nan
         case default
            frz(i0, j0) = nan
         end select
         call import_map(ms, cice, oin, frz)
         call ice_import_ocean(grid, eos, ms, frz, cice, oin, ierr)
         !$acc update self(oin%sst, oin%sss, oin%sss_ocn, oin%Tf, oin%frzmlt, cice%ice_err)
         nonfin = .not. any([ieee_is_finite(oin%sst(i0, j0)), ieee_is_finite(oin%sss(i0, j0)), &
                             ieee_is_finite(oin%sss_ocn(i0, j0)), ieee_is_finite(oin%Tf(i0, j0)), &
                             ieee_is_finite(oin%frzmlt(i0, j0))])
         bit = iand(cice%ice_err(i0, j0), ICE_ERR_NONFINITE) /= 0_int32
         call check(error, nonfin, "case: all five outputs are NaN, never a plausible number")
         call check(error, bit, "case: ICE_ERR_NONFINITE is set")
         nbad = 0
         if (any(cice%ice_err(IS:IE, JS:JE) /= 0_int32)) nbad = count(cice%ice_err(IS:IE, JS:JE) /= 0_int32)
         call check(error, nbad == 1, "only the offending cell is flagged")
         call import_unmap(ms, cice, oin, frz)
         call oin%destroy()
         call cice%destroy()
         call eos%destroy()
         call ms%destroy()
         if (allocated(error)) return
      end do
   end subroutine test_import_nan

   subroutine test_import_cadence(error)
      !! The provider -> import chain at the ICE step: with bank B and the ice dt the imported
      !! frzmlt is B/dt_ice; once the bank is spent (zeroed, as `add_new_ice` does) and the
      !! provider re-run, the next import reads 0, so the same bank is never given to the
      !! ice twice.  (A provider run at thermo cadence with the window dt would give
      !! B/(ratio dt) and then keep returning it: the wiring must run it every ice step.)
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: DT = 900.0_wp, BANK = 2.7e6_wp
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ice_cice_state_t) :: cice
      type(ice_ocean_in_t) :: oin
      type(ocean_sea_ice_t) :: ice
      real(wp), allocatable :: frz(:, :)
      integer :: ierr
      real(wp) :: want, got, got2

      call import_setup(grid, ms, eos, cice, oin, frz, 3)
      ms%tracers(ms%idx_temperature)%c = -1.9_wp
      ms%tracers(ms%idx_salinity)%c = 34.0_wp
      ice%enable = .true.
      ice%melt_potential_enable = .true.
      call ice%init(grid)
      ice%frazil_heat = BANK
      call import_map(ms, cice, oin, frz)
      call ice%enter_data()
      call ice_melt_potential(grid, eos, ms, ice%frazil_heat, 10.0_wp, DT, ice%melt_potential, ice%frzmlt)
      call ice_import_ocean(grid, eos, ms, ice%frzmlt, cice, oin, ierr)
      !$acc update self(oin%frzmlt)
      want = BANK/DT
      got = oin%frzmlt(IS + 1, JS + 1)
      call check(error, abs(got - want) <= 1.0e-14_wp*want, "imported frzmlt = bank/dt_ice")
      ice%frazil_heat = 0.0_wp
      !$acc update device(ice%frazil_heat)
      call ice_melt_potential(grid, eos, ms, ice%frazil_heat, 10.0_wp, DT, ice%melt_potential, ice%frzmlt)
      call ice_import_ocean(grid, eos, ms, ice%frzmlt, cice, oin, ierr)
      !$acc update self(oin%frzmlt)
      got2 = oin%frzmlt(IS + 1, JS + 1)
      call check(error, got2 <= 0.0_wp .and. ieee_is_finite(got2), "spent bank: no freezing potential next step")
      call ice%exit_data()
      call ice%destroy()
      call import_unmap(ms, cice, oin, frz)
   end subroutine test_import_cadence

   ! ------------------------------------------------------------------ hooks

   subroutine test_open_fraction(error)
      !! `ice_aggregates_fill` copies the CICE aggregates into the scheme-neutral slot arrays
      !! when the CICE slot is initialised (every cell, ghosts included), and is the SIS2 fill
      !! (aice0 = part_size(0)) otherwise.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: g
      type(ocean_sea_ice_t) :: ice, ice2
      type(ice_params_t) :: p
      integer :: ierr, i, j
      real(wp) :: v

      call g%init(NXP, NYP, NG, 1000.0_wp, 1000.0_wp)
      ice%enable = .true.
      ice%ncat = 2
      ice%nk_ice = 2
      call ice%init(g)
      p = ice_params_default(ICE_CONSTANTS_ICEPACK, NCAT, NILYR, 1)
      call ice%cice%init(g, p, ierr)
      do j = 1, NY
         do i = 1, NX
            v = real(i + 17*j, wp)/101.0_wp
            ice%cice%aice(i, j) = v
            ice%cice%aice0(i, j) = 1.0_wp - v
            ice%cice%vice(i, j) = 2.0_wp*v
            ice%cice%vsno(i, j) = 3.0_wp*v
         end do
      end do
      call ice%enter_data()
      !$acc update self(ice%aice, ice%aice0, ice%vice, ice%vsno)
      call check(error, all(ice%aice == ice%cice%aice) .and. all(ice%aice0 == ice%cice%aice0) .and. &
                 all(ice%vice == ice%cice%vice) .and. all(ice%vsno == ice%cice%vsno), &
                 "aggregates are the CICE aggregates on every cell")
      call ice%exit_data()
      call ice%destroy()
      ice2%enable = .true.
      ice2%ncat = 2
      ice2%nk_ice = 2
      call ice2%init(g)
      call ice2%enter_data()
      !$acc update self(ice2%aice0)
      call check(error, all(ice2%aice0 == ice2%part_size(:, :, 0)), "SIS2 path unchanged: aice0 = part_size(0)")
      call ice2%exit_data()
      call ice2%destroy()
   end subroutine test_open_fraction

   subroutine test_stress_blend(error)
      !! Face blend with `a = 0.6` everywhere: interior faces tau = 0.4 tau_atm + 0.6 tau_ice; the
      !! magnitude is re-derived.  a = 0 gives the wind, a = 1 the ice stress.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: g
      type(ice_cice_state_t) :: cice
      type(ice_params_t) :: p
      type(ocean_surface_stress_t) :: stress
      type(ice_cice_stress_mirror_t) :: mirror
      real(wp), allocatable :: ta_x(:, :), ta_y(:, :), ti_x(:, :), ti_y(:, :)
      integer :: ierr, k
      real(wp) :: a, ex, ey

      call g%init(NXP, NYP, NG, 1000.0_wp, 1000.0_wp)
      p = ice_params_default(ICE_CONSTANTS_ICEPACK, NCAT, NILYR, 1)
      call cice%init(g, p, ierr)
      call stress%init(g, nz_ml=1)
      call mirror%init(g)
      call mirror%enter_data()
      call stress%set_wind_stress_const(0.0_wp, 0.0_wp)
      allocate (ta_x(NX + 1, NY), ta_y(NX, NY + 1), ti_x(NX + 1, NY), ti_y(NX, NY + 1))
      ta_x = 0.3_wp
      ta_y = 0.1_wp
      ti_x = 0.7_wp
      ti_y = -0.2_wp
      call cice%enter_data()
      call stress%enter_data()
      !$acc enter data copyin(ta_x, ta_y, ti_x, ti_y)
      do k = 1, 3
         select case (k)
         case (1)
            a = 0.6_wp
         case (2)
            a = 0.0_wp
         case default
            a = 1.0_wp
         end select
         cice%aice = a
         !$acc update device(cice%aice)
         call ice_cice_stress_blend(stress, ta_x, ta_y, ti_x, ti_y, cice, mirror)
         !$acc update self(stress%tau_x, stress%tau_y, stress%stress_mag)
         ex = (1.0_wp - a)*0.3_wp + a*0.7_wp
         ey = (1.0_wp - a)*0.1_wp + a*(-0.2_wp)
         call check(error, all(abs(stress%tau_x(IS:IE + 1, JS:JE) - ex) <= 1.0e-14_wp) .and. &
                    all(abs(stress%tau_y(IS:IE, JS:JE + 1) - ey) <= 1.0e-14_wp), "blend closed form")
         call check(error, all(abs(stress%stress_mag(IS:IE, JS:JE) - sqrt(ex*ex + ey*ey)) <= 1.0e-12_wp), &
                    "stress_mag re-derived from the blended pair")
         if (allocated(error)) exit
      end do
      !$acc exit data delete(ta_x, ta_y, ti_x, ti_y)
      call mirror%exit_data()
      call stress%exit_data()
      call cice%exit_data()
   end subroutine test_stress_blend

   ! ------------------------------------------------- fix-round tests

   subroutine stress_setup(g, cice, stress, mirror, a, ta, ti)
      !! Shared setup of the blend tests: uniform cover `a`, uniform unweighted wind `ta`
      !! (both components), uniform ice stress `ti`.
      type(hgrid_t), intent(out) :: g
      type(ice_cice_state_t), intent(inout) :: cice
      type(ocean_surface_stress_t), intent(inout) :: stress
      type(ice_cice_stress_mirror_t), intent(inout) :: mirror
      real(wp), intent(in) :: a, ta, ti
      type(ice_params_t) :: p
      integer :: ierr

      call g%init(NXP, NYP, NG, 1000.0_wp, 1000.0_wp)
      p = ice_params_default(ICE_CONSTANTS_ICEPACK, NCAT, NILYR, 1)
      call cice%init(g, p, ierr)
      cice%aice = a
      call stress%init(g, nz_ml=1)
      call stress%set_wind_stress_const(0.0_wp, 0.0_wp)
      call mirror%init(g)
      call cice%enter_data()
      call stress%enter_data()
      call mirror%enter_data()
      !$acc update device(cice%aice)
   end subroutine stress_setup

   subroutine test_stress_unweighted(error)
      !! The blend's `tau_a` is the UNWEIGHTED wind: with cover a = 0.6 and wind 0.5 the ocean
      !! gets 0.4*0.5 + 0.6*ti.  Feeding the bulk forcing's already open-water-weighted wind
      !! (0.4*0.5) instead weights it by 0.4 twice: 0.4*(0.4*0.5) + 0.6*ti, wrong by 0.12.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: g
      type(ice_cice_state_t) :: cice
      type(ocean_surface_stress_t) :: stress
      type(ice_cice_stress_mirror_t) :: mirror
      real(wp), allocatable :: ta(:, :), tw(:, :), ti(:, :), tay(:, :), tiy(:, :)
      real(wp) :: right, double

      call stress_setup(g, cice, stress, mirror, 0.6_wp, 0.5_wp, 0.2_wp)
      allocate (ta(NX + 1, NY), tw(NX + 1, NY), ti(NX + 1, NY), tay(NX, NY + 1), tiy(NX, NY + 1))
      ta = 0.5_wp
      tw = 0.4_wp*0.5_wp
      ti = 0.2_wp
      tay = 0.0_wp
      tiy = 0.0_wp
      !$acc enter data copyin(ta, tw, ti, tay, tiy)
      call ice_cice_stress_blend(stress, ta, tay, ti, tiy, cice, mirror)
      !$acc update self(stress%tau_x)
      right = stress%tau_x(IS + 1, JS + 1)
      call ice_cice_stress_blend(stress, tw, tay, ti, tiy, cice, mirror)
      !$acc update self(stress%tau_x)
      double = stress%tau_x(IS + 1, JS + 1)
      call check(error, abs(right - (0.4_wp*0.5_wp + 0.6_wp*0.2_wp)) <= 1.0e-14_wp, "unweighted wind: 0.4 tau_a + 0.6 tau_ice")
      call check(error, abs((right - double) - 0.4_wp*0.5_wp*0.6_wp) <= 1.0e-14_wp, &
                 "a pre-weighted wind is double-weighted (off by 0.12): the contract is the unweighted stress")
      !$acc exit data delete(ta, tw, ti, tay, tiy)
      call mirror%exit_data()
      call stress%exit_data()
      call cice%exit_data()
   end subroutine test_stress_unweighted

   subroutine test_stress_restart(error)
      !! The blended face stress survives a restart: write the mirror through the registry, read it
      !! into a fresh mirror, apply the resume to a re-seeded (wind-only) stress slot: the faces equal
      !! the blend bit for bit.  A fresh run (valid = 0) keeps the wind.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: g
      type(ice_cice_state_t) :: cice
      type(ocean_surface_stress_t) :: stress, s2, s3
      type(ice_cice_stress_mirror_t) :: mirror, m2
      type(restart_registry_t) :: ra, rb
      type(decomp_t) :: d
      type(ocean_restart_metadata_t) :: meta
      character(len=*), parameter :: FN = "test_cice_coupling_stress.nc"
      real(wp), allocatable :: ta(:, :), ti(:, :), tay(:, :), tiy(:, :)
      real(wp) :: t
      integer :: ierr, step, osc, lun, i, j
      integer(int64) :: sd

      call stress_setup(g, cice, stress, mirror, 0.35_wp, 0.5_wp, 0.2_wp)
      allocate (ta(NX + 1, NY), ti(NX + 1, NY), tay(NX, NY + 1), tiy(NX, NY + 1))
      sd = 321_int64
      do j = 1, NY
         do i = 1, NX + 1
            ta(i, j) = 0.5_wp*rnd(sd)
            ti(i, j) = 0.5_wp*rnd(sd)
         end do
      end do
      do j = 1, NY + 1
         do i = 1, NX
            tay(i, j) = 0.5_wp*rnd(sd)
            tiy(i, j) = 0.5_wp*rnd(sd)
         end do
      end do
      !$acc enter data copyin(ta, ti, tay, tiy)
      call ice_cice_stress_blend(stress, ta, tay, ti, tiy, cice, mirror)
      !$acc update self(stress%tau_x, stress%tau_y, mirror%tau_x, mirror%tau_y)
      call check(error, mirror%valid == 1.0_wp .and. all(mirror%tau_x == stress%tau_x) .and. &
                 all(mirror%tau_y == stress%tau_y), "the mirror is the blend, bit for bit")
      call mirror%register_restart(ra)
      call m2%init(g)
      call m2%register_restart(rb)
      call decomp_init(d, NXP, NYP, 1, 1, 0)
      meta%nghost = 0
      open (newunit=lun, file=FN, status="replace")
      close (lun, status="delete")
      call ocean_restart_write_local(FN, ra, d, meta, 1.0_wp, 1, 1, ierr)
      call check(error, ierr, 0, "write")
      call ocean_restart_read_local(FN, rb, d, meta, t, step, osc, ierr)
      call check(error, ierr, 0, "read")
      open (newunit=lun, file=FN, status="old")
      close (lun, status="delete")
      if (allocated(error)) goto 900
      call s2%init(g, nz_ml=1)
      call s2%set_wind_stress_const(0.7_wp, 0.0_wp)
      call ice_cice_stress_resume_apply(s2, m2)
      call check(error, m2%valid == 1.0_wp, "valid flag restored")
      call check(error, all(s2%tau_x == stress%tau_x) .and. all(s2%tau_y == stress%tau_y), &
                 "resumed faces are the blended stress, not the wind")
      call s3%init(g, nz_ml=1)
      call s3%set_wind_stress_const(0.7_wp, 0.0_wp)
      call mirror%init(g)
      call ice_cice_stress_resume_apply(s3, mirror)
      call check(error, all(s3%tau_x == 0.7_wp), "a fresh run keeps the wind")
900   continue
      !$acc exit data delete(ta, ti, tay, tiy)
      call stress%exit_data()
      call cice%exit_data()
   end subroutine test_stress_restart

   subroutine stub_growth(fx, ratio_scale, step)
      !! Growing-ice stub: dhi = 1 cm * cell factor per step; fresh = -rhoi dhi/dt, fsalt = -Sref 1e-3 rhoi dhi/dt
      !! (update_ocn_f form), a heat flux and an under-ice shortwave.
      type(ice_flux_t), intent(inout) :: fx
      real(wp), intent(in) :: ratio_scale
      integer, intent(in) :: step
      integer :: i, j
      real(wp) :: dhi, m

      m = 1.0_wp + 0.2_wp*real(modulo(step, 3), wp)
      do j = 1, NY
         do i = 1, NX
            dhi = 0.01_wp*(1.0_wp + 0.1_wp*real(i, wp) + 0.05_wp*real(j, wp))*m*ratio_scale
            fx%fresh(i, j) = -RHOI*dhi/900.0_wp
            fx%fsalt(i, j) = -SREF*1.0e-3_wp*RHOI*dhi/900.0_wp
            fx%fhocn(i, j) = 4.0_wp + 0.5_wp*real(i, wp)
            fx%fswthru(i, j) = 3.0_wp + 0.25_wp*real(j, wp)
            fx%fswthru_vis(i, j) = fx%fswthru(i, j)
         end do
      end do
      !$acc update device(fx%fresh, fx%fsalt, fx%fhocn, fx%fswthru, fx%fswthru_vis)
   end subroutine stub_growth

   subroutine test_ice_components_end_to_end(error)
      !! B1 + S5.  Open-water shortwave (`q_sw`, as the bulk forcing writes it) and a reader's
      !! `heat_added` are present when the ice hands off; the assembler sums the ice components
      !! beside them, and the REAL surface-flux apply changes the ocean's content by exactly what
      !! the ice emitted plus the other writers' share: sum(h*S) rises (growing ice), and
      !! `rho0 cp d(hT) = (q_sw + heat_added) W + sum((fhocn + fswthru) dt)`, per cell, ratio 1 and 3.
      type(error_type), allocatable, intent(out) :: error
      integer :: ratio, ierr, step, k, i, j, it, is_
      type(hgrid_t) :: g
      type(ice_flux_t) :: fx
      type(ice_flux_avg_t) :: av
      type(ocean_surface_flux_t) :: sf
      type(multilayer_state_t) :: ms
      type(ocean_bc_state_t) :: bc
      real(wp), allocatable :: sss(:, :), s0(:, :, :), t0(:, :, :), eh(:, :), es(:, :)
      real(wp) :: dt, w, dh, ds, dtt, rho_cp, worst_h, worst_s, tot, htop
      integer(int64) :: sd

      do ratio = 1, 3, 2
         call make_flux(g, fx, ierr)
         call av%init(g)
         call make_halo(g, bc, .false.)
         ms%nz_ml = 3
         call ms%init(g)
         call sf%init(g)
         call sf%set_components(g, .true.)
         call sf%set_ice_components(g, .true.)
         it = ms%idx_temperature
         is_ = ms%idx_salinity
         ms%h_layer = 10.0_wp
         ms%tracers(it)%c = 5.0_wp
         allocate (sss(NX, NY), s0(NX, NY, 3), t0(NX, NY, 3), eh(NX, NY), es(NX, NY))
         sd = 11_int64
         do j = 1, NY
            do i = 1, NX
               ms%tracers(is_)%c(i, j, :) = 30.0_wp + 6.0_wp*rnd(sd)
            end do
         end do
         s0 = ms%tracers(is_)%c
         t0 = ms%tracers(it)%c
         sss = s0(:, :, 3)
         sf%q_sw = 50.0_wp
         sf%heat_added = 2.0_wp
         sf%has_q_sw = .true.
         !$acc enter data copyin(ms)
         call ms%enter_data()
         !$acc update device(ms%wet_mask, ms%h_layer, ms%k_top)
         !$acc update device(ms%tracers(it)%c, ms%tracers(is_)%c)
         !$acc enter data copyin(sss)
         call fx%enter_data()
         call av%enter_data()
         call sf%enter_data()
         eh = 0.0_wp
         es = 0.0_wp
         w = 0.0_wp
         do step = 1, ratio
            dt = 900.0_wp
            call stub_growth(fx, 1.0_wp, step)
            call ice_flux_avg_accumulate(av, fx, dt, ierr)
            w = w + dt
            eh = eh + (fx%fhocn + fx%fswthru)*dt
            es = es + (PSU_PER_MASS_FRACTION*fx%fsalt - sss*fx%fresh)*dt
         end do
         call ice_cice_apply_window(g, bc, av, sf, sss, ierr, window_expected=w)
         call check(error, ierr, ICE_ADAPT_OK, "hand-off")
         call ocean_surface_flux_assemble(g, sf, ms)
         call ocean_surface_flux_apply_tracers(g, sf, ms, w)
         !$acc update self(ms%tracers(it)%c, ms%tracers(is_)%c)
         worst_h = 0.0_wp
         worst_s = 0.0_wp
         tot = 0.0_wp
         rho_cp = sf%rho0*sf%cp
         do j = JS, JE
            do i = IS, IE
               htop = ms%h_layer(i, j, 3)
               dtt = (ms%tracers(it)%c(i, j, 3) - t0(i, j, 3))*htop
               ds = (ms%tracers(is_)%c(i, j, 3) - s0(i, j, 3))*htop
               dh = (50.0_wp + 2.0_wp)*w + eh(i, j)
               worst_h = max(worst_h, abs(dtt*rho_cp - dh)/abs(dh))
               worst_s = max(worst_s, abs(ds*sf%rho0 - es(i, j))/abs(es(i, j)))
               tot = tot + ds
               if (.not. ieee_is_finite(dtt + ds)) worst_h = ieee_value(1.0_wp, ieee_quiet_nan)
            end do
         end do
         write (error_unit, "(a,i0,a,2es10.2)") "ocean content closure ratio=", ratio, " rel residual heat/salt =", &
            worst_h, worst_s
         call check(error, ieee_is_finite(worst_h) .and. ieee_is_finite(worst_s), "content changes are finite")
         if (allocated(error)) goto 800
         call check(error, tot > 0.0_wp .and. all(es(IS:IE, JS:JE) > 0.0_wp), "growing ice salinifies the ocean (sum(h*S) rises)")
         call check(error, worst_h <= 3.0e-12_wp, "ocean heat = open-water SW + reader heat + ice heat + under-ice SW")
         call check(error, worst_s <= 1.0e-12_wp, "ocean salt change = what the ice emitted (virtual-salt form)")
         if (allocated(error)) goto 800
800      continue
         call sf%exit_data()
         call av%exit_data()
         call fx%exit_data()
         !$acc exit data delete(sss)
         call ms%exit_data()
         !$acc exit data delete(ms)
         deallocate (sss, s0, t0, eh, es)
         call ms%destroy()
         call sf%destroy()
         call av%destroy()
         if (allocated(error)) return
      end do
   end subroutine test_ice_components_end_to_end

   subroutine test_ghosts(error)
      !! Periodic domain: after the hand-off every seam ghost of the ice components, and of the
      !! returned mass flux, equals its owner (two-pass, corners included).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: g
      type(ice_flux_t) :: fx
      type(ice_flux_avg_t) :: av
      type(ocean_surface_flux_t) :: sf
      type(ocean_bc_state_t) :: bc
      real(wp), allocatable :: sss(:, :), fw(:, :)
      integer :: ierr, i, j, io, jo, nbad

      call make_flux(g, fx, ierr)
      call av%init(g)
      call sf%init(g)
      call sf%set_components(g, .true.)
      call sf%set_ice_components(g, .true.)
      call make_halo(g, bc, .true.)
      allocate (sss(NX, NY), fw(NX, NY))
      sss = 34.0_wp
      fw = -1.0_wp
      !$acc enter data copyin(sss, fw)
      call fx%enter_data()
      call av%enter_data()
      call sf%enter_data()
      call stub_growth(fx, 1.0_wp, 1)
      call ice_flux_avg_accumulate(av, fx, 900.0_wp, ierr)
      call ice_cice_apply_window(g, bc, av, sf, sss, ierr, fresh_mass_out=fw)
      call check(error, ierr, ICE_ADAPT_OK, "hand-off")
      !$acc update self(sf%ice_heat, sf%ice_sw, sf%ice_salt, fw)
      nbad = 0
      do j = 1, NY
         do i = 1, NX
            io = IS + modulo(i - IS, NXP)
            jo = JS + modulo(j - JS, NYP)
            if (sf%ice_heat(i, j) /= sf%ice_heat(io, jo) .or. sf%ice_sw(i, j) /= sf%ice_sw(io, jo) .or. &
                sf%ice_salt(i, j) /= sf%ice_salt(io, jo) .or. fw(i, j) /= fw(io, jo)) nbad = nbad + 1
            if (.not. ieee_is_finite(sf%ice_salt(i, j))) nbad = nbad + 1
         end do
      end do
      call check(error, nbad == 0, "every seam ghost equals its owner")
      call check(error, sf%ice_heat(1, 1) /= 0.0_wp .and. sf%ice_salt(1, JS) /= 0.0_wp, "and the ghosts are not just zero")
      call sf%exit_data()
      call av%exit_data()
      call fx%exit_data()
      !$acc exit data delete(sss, fw)
   end subroutine test_ghosts

   subroutine test_melt_cap(error)
      !! Ratio 3, strong melt: provider -> capped import -> a stub ice that takes everything offered
      !! (fhocn = frzmlt), accumulated in the open window.  Total heat taken over the window never
      !! exceeds the layer's available heat MP (equals it: the first step takes all); the same chain
      !! without the cap draws 3 MP.  Plus closed forms of `ice_melt_cap`.
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: DT = 900.0_wp
      type(hgrid_t) :: grid, g
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ice_cice_state_t) :: cice
      type(ice_ocean_in_t) :: oin
      type(ice_flux_t) :: fx
      type(ice_flux_avg_t) :: av
      type(ocean_sea_ice_t) :: ice
      real(wp), allocatable :: frz(:, :), taken(:, :), taken_nocap(:, :)
      integer :: ierr, step, mode, i, j
      real(wp) :: worst, mp, nan

      nan = ieee_value(1.0_wp, ieee_quiet_nan)
      call check(error, abs(ice_melt_cap(-100.0_wp, -20.0_wp, 1800.0_wp, 900.0_wp, nan) + 60.0_wp) <= 1.0e-12_wp, &
                 "cap closed form: 90000 J offered, 36000 taken, 54000 left over 900 s")
      call check(error, ice_melt_cap(-100.0_wp, -60.0_wp, 1800.0_wp, 900.0_wp, nan) == 0.0_wp, "exhausted: nothing offered")
      call check(error, ice_melt_cap(40.0_wp, -60.0_wp, 1800.0_wp, 900.0_wp, nan) == 40.0_wp, "freezing side untouched")
      call check(error, ice_melt_cap(-100.0_wp, 5.0_wp, 1800.0_wp, 900.0_wp, nan) == -100.0_wp, &
                 "heat given back to the ocean does not raise the offer")
      call check(error,.not. ieee_is_finite(ice_melt_cap(-100.0_wp, nan, 1800.0_wp, 900.0_wp, nan)), "NaN window poisons")
      if (allocated(error)) return
      call import_setup(grid, ms, eos, cice, oin, frz, 3)
      ms%tracers(ms%idx_temperature)%c = 3.0_wp
      ms%tracers(ms%idx_salinity)%c = 34.0_wp
      ms%h_layer = 10.0_wp
      ice%enable = .true.
      ice%melt_potential_enable = .true.
      call ice%init(grid)
      call fx%init(grid, NCAT, ierr)
      call av%init(grid)
      call import_map(ms, cice, oin, frz)
      call ice%enter_data()
      call fx%enter_data()
      call av%enter_data()
      allocate (taken(NX, NY), taken_nocap(NX, NY))
      do mode = 1, 2
         call ice_flux_avg_reset(av)
         taken = 0.0_wp
         do step = 1, 3
            call ice_melt_potential(grid, eos, ms, ice%frazil_heat, 10.0_wp, DT, ice%melt_potential, ice%frzmlt)
            if (mode == 1) then
               call ice_import_ocean(grid, eos, ms, ice%frzmlt, cice, oin, ierr, avg=av, dt=DT)
            else
               call ice_import_ocean(grid, eos, ms, ice%frzmlt, cice, oin, ierr)
            end if
            !$acc update self(oin%frzmlt)
            fx%fresh = 0.0_wp
            fx%fsalt = 0.0_wp
            fx%fswthru = 0.0_wp
            fx%fswthru_vis = 0.0_wp
            fx%fhocn = oin%frzmlt
            !$acc update device(fx%fresh, fx%fsalt, fx%fhocn, fx%fswthru, fx%fswthru_vis)
            taken = taken - oin%frzmlt*DT
            call ice_flux_avg_accumulate(av, fx, DT, ierr)
         end do
         if (mode == 1) then
            !$acc update self(ice%melt_potential)
            worst = 0.0_wp
            do j = JS, JE
               do i = IS, IE
                  mp = ice%melt_potential(i, j)
                  if (.not. (mp > 0.0_wp)) worst = ieee_value(1.0_wp, ieee_quiet_nan)
                  if (.not. ieee_is_finite(taken(i, j))) then
                     worst = ieee_value(1.0_wp, ieee_quiet_nan)
                  else
                     worst = max(worst, taken(i, j)/mp - 1.0_wp)
                  end if
               end do
            end do
            write (error_unit, "(a,es10.2)") "melt cap ratio 3: max (taken/MP - 1) =", worst
            call check(error, ieee_is_finite(worst), "finite")
            call check(error, worst <= 1.0e-12_wp, "capped: total heat taken over the window never exceeds MP")
            call check(error, minval(taken(IS:IE, JS:JE)/ice%melt_potential(IS:IE, JS:JE)) >= 1.0_wp - 1.0e-12_wp, &
                       "and the ice does take all of it")
         else
            taken_nocap = taken
            call check(error, minval(taken_nocap(IS:IE, JS:JE)/ice%melt_potential(IS:IE, JS:JE)) > 2.99_wp, &
                       "uncapped: the same upper-ocean heat is drawn three times")
         end if
         if (allocated(error)) goto 900
      end do
900   continue
      call av%exit_data()
      call fx%exit_data()
      call ice%exit_data()
      call import_unmap(ms, cice, oin, frz)
   end subroutine test_melt_cap

   subroutine test_melt_cap_ratio1(error)
      !! At the first step of a window (t_acc = 0: every ratio-1 step) the capped import is bit for
      !! bit the uncapped one, melt and freezing cells alike; the cap also changes nothing when the
      !! window has only given heat back (fhocn > 0).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ice_cice_state_t) :: cice
      type(ice_ocean_in_t) :: oin
      type(ice_flux_avg_t) :: av
      real(wp), allocatable :: frz(:, :), base(:, :)
      integer :: ierr

      call import_setup(grid, ms, eos, cice, oin, frz, 3)
      call av%init(grid)
      call import_map(ms, cice, oin, frz)
      call av%enter_data()
      call ice_import_ocean(grid, eos, ms, frz, cice, oin, ierr)
      !$acc update self(oin%frzmlt)
      base = oin%frzmlt
      call ice_import_ocean(grid, eos, ms, frz, cice, oin, ierr, avg=av, dt=900.0_wp)
      !$acc update self(oin%frzmlt)
      call check(error, all(oin%frzmlt == base), "empty window: capped == uncapped, bitwise")
      call check(error, any(base(IS:IE, JS:JE) < 0.0_wp) .and. any(base(IS:IE, JS:JE) > 0.0_wp), &
                 "the state has both melt and freezing cells")
      call av%exit_data()
      call import_unmap(ms, cice, oin, frz)
   end subroutine test_melt_cap_ratio1

   subroutine test_required_and_window(error)
      !! Restart tags are REQUIRED when asked (file written by the new scheme); the hand-off checks the
      !! accumulated time against the thermo-window length and refuses without the ice components.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: g
      type(ice_flux_t) :: fx
      type(ice_flux_avg_t) :: av
      type(ocean_surface_flux_t) :: sf, sf2
      type(ocean_bc_state_t) :: bc
      type(restart_registry_t) :: r1, r2
      real(wp), allocatable :: sss(:, :)
      real(wp) :: tw
      integer :: ierr, e

      call make_flux(g, fx, ierr)
      call av%init(g)
      call av%register_restart(r1)
      call av%register_restart(r2, required=.true.)
      call check(error, all(r1%entries(1:r1%n)%optional) .and. .not. any(r2%entries(1:r2%n)%optional) .and. &
                 r1%n == 7 .and. r2%n == 7, "required=.true. makes all seven tags non-optional")
      call make_halo(g, bc, .false.)
      call sf%init(g)
      call sf%set_components(g, .true.)
      call sf%set_ice_components(g, .true.)
      call sf2%init(g)
      call sf2%set_components(g, .true.)
      allocate (sss(NX, NY))
      sss = 34.0_wp
      !$acc enter data copyin(sss)
      call fx%enter_data()
      call av%enter_data()
      call sf%enter_data()
      call sf2%enter_data()
      call stub_growth(fx, 1.0_wp, 1)
      call ice_flux_avg_accumulate(av, fx, 900.0_wp, ierr)
      call ice_cice_apply_window(g, bc, av, sf, sss, ierr, window_expected=1800.0_wp, t_window=tw)
      call check(error, ierr, ICE_ADAPT_ERR_WINDOW, "a window of 900 s against an expected 1800 s is refused")
      call check(error, tw == 900.0_wp .and. av%n_acc == 1.0_wp, "and left untouched; the accumulated time is reported")
      call ice_cice_apply_window(g, bc, av, sf2, sss, ierr)
      call check(error, ierr, ICE_ADAPT_ERR_NO_COMPONENT, "components on without the ice components is refused")
      call ice_cice_apply_window(g, bc, av, sf, sss, ierr, window_expected=900.0_wp, t_window=tw)
      call check(error, ierr, ICE_ADAPT_OK, "the right window passes")
      e = 0
      call sf2%exit_data()
      call sf%exit_data()
      call av%exit_data()
      call fx%exit_data()
      !$acc exit data delete(sss)
   end subroutine test_required_and_window

   subroutine test_nothing_left(error)
      !! Last test of the file: with RDB_ACC_PRESENT_DUMP set the OpenACC present table is dumped
      !! (offload build); every earlier test unmaps what it mapped, so nothing of this file's
      !! slots may appear in it.  Always passes (the dump is read by eye / by the repeat-run gate).
      type(error_type), allocatable, intent(out) :: error
      character(len=8) :: v
      integer :: st

      call get_environment_variable("RDB_ACC_PRESENT_DUMP", v, status=st)
#ifdef RDB_GPU_OFFLOAD
      if (st == 0) call acc_present_dump()
#endif
      call check(error, .true., "present table dumped on request")
   end subroutine test_nothing_left

end module test_cice_coupling
