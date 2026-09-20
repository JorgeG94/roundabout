!! Unit tests for namelist configuration reading
module test_config
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, DEG2RAD, NZ_STACK_MAX, &
                            nz_stack_required, nz_stack_is_sufficient
   use rdb_config, only: config_t, read_config, read_config_from_string, &
                         resolve_bt_halo, bt_halo_auto_exclusion, &
                         BT_HALO_AUTO_SENTINEL, BT_HALO_AUTO_WIDTH, &
                         diag_density_levels_ok, MAX_OCEAN_DIAG_Z_LEVELS, &
                         ice_hlim_count, ice_hlim_spec_is_valid, MAX_ICE_HLIM_VALS, &
                         validate_config
   implicit none
   private

   public :: collect_config_tests

contains

   subroutine portable_sleep(seconds)
      !! Portable sleep using system_clock (works on gfortran + nvhpc)
      integer, intent(in) :: seconds
      integer :: count_start, count_end, count_rate

      call system_clock(count_start, count_rate)
      do
         call system_clock(count_end)
         if ((count_end - count_start) >= seconds*count_rate) exit
      end do
   end subroutine portable_sleep

   subroutine wait_for_file(filename)
      !! Wait until file is visible and readable (parallel FS flush delay)
      character(len=*), intent(in) :: filename
      logical :: exists
      integer :: i, sz

      call portable_sleep(1)

      do i = 1, 100
         inquire (file=filename, exist=exists, size=sz)
         if (exists .and. sz > 0) return
         call portable_sleep(1)
      end do
   end subroutine wait_for_file

   subroutine collect_config_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("defaults", test_defaults), &
                  new_unittest("read_grid_nml", test_read_grid), &
                  new_unittest("read_all_namelists", test_read_all), &
                  new_unittest("missing_file_keeps_defaults", test_missing_file), &
                  new_unittest("read_ocean_subnmls", test_read_ocean_subnmls), &
                  new_unittest("ocean_diag_dt_out_time_unit", test_ocean_diag_time_unit), &
                  new_unittest("cartesian_degrees_extent", test_cartesian_degrees), &
                  new_unittest("cartesian_meters_default_bitident", test_cartesian_meters_default), &
                  new_unittest("bt_halo_default_is_auto_sentinel", test_bt_halo_default_sentinel), &
                  new_unittest("resolve_bt_halo", test_resolve_bt_halo), &
                  new_unittest("bt_halo_auto_exclusion", test_bt_halo_auto_exclusion), &
                  new_unittest("diag_density_requires_rho_levels", test_diag_density_requires_rho_levels), &
                  new_unittest("ice_hlim_count_leading_run", test_ice_hlim_count_leading_run), &
                  new_unittest("ice_hlim_spec_validity", test_ice_hlim_spec_validity), &
                  new_unittest("nz_stack_guard_refuses_oversized_nz", test_nz_stack_guard) &
                  ]
   end subroutine collect_config_tests

   subroutine test_nz_stack_guard(error)
      !! Per-column stack-workspace guard.  Every layered kernel (BPG, ALE
      !! remap, kappa-shear, Redi, diag remap) carries fixed-size
      !! `NZ_STACK_MAX` thread-local column arrays; overrunning them
      !! corrupts thread-local storage SILENTLY (wrong answers, no crash).
      !! `validate_config` therefore refuses such a run outright.  Tested
      !! at the predicate level (house pattern — see
      !! `test_diag_density_requires_rho_levels`) rather than by exercising
      !! the fatal `error stop`.
      !!
      !! The requirement is `nz + 1`, NOT `2*nz + 2`: Redi's neutral-surface
      !! locals are declared `2*NZ_STACK_MAX + 2` and scale with the constant.
      type(error_type), allocatable, intent(out) :: error

      checks: block
         ! The rule itself: nz + 1.
         call check(error, nz_stack_required(50) == 51, &
                    "nz_stack_required(50) must be 51 (nz + 1)")
         if (allocated(error)) exit checks

         ! A comfortably-sized run is accepted.
         call check(error, nz_stack_is_sufficient(1), &
                    "nz = 1 must always fit")
         if (allocated(error)) exit checks

         ! The exact boundary: nz = NZ_STACK_MAX - 1 needs NZ_STACK_MAX, fits.
         call check(error, nz_stack_is_sufficient(NZ_STACK_MAX - 1), &
                    "nz = NZ_STACK_MAX - 1 needs exactly NZ_STACK_MAX => must fit")
         if (allocated(error)) exit checks

         ! One past it: nz = NZ_STACK_MAX needs NZ_STACK_MAX + 1 => REFUSED.
         ! This is the off-by-one that the `+1`-indexed arrays (k-eps S2/N2,
         ! si_closures kapS/kapT, ml_settling_tracer) would overrun.
         call check(error,.not. nz_stack_is_sufficient(NZ_STACK_MAX), &
                    "nz = NZ_STACK_MAX needs nz+1 slots => must be REFUSED")
         if (allocated(error)) exit checks

         ! Grossly oversized is refused.
         call check(error,.not. nz_stack_is_sufficient(NZ_STACK_MAX + 500), &
                    "nz far above NZ_STACK_MAX must be REFUSED")
         if (allocated(error)) exit checks

         ! The shipped envelope must actually be buildable with the
         ! compiled default: nz = 90 is the largest shipped namelist
         ! (acc_channel_kitchensink_xl) and nz = 100 the largest test.
         ! If either of these ever fails, the default was lowered too far.
         call check(error, nz_stack_is_sufficient(90), &
                    "largest shipped namelist (nz=90) must fit the compiled NZ_STACK_MAX")
         if (allocated(error)) exit checks

         call check(error, nz_stack_is_sufficient(100), &
                    "largest test case (nz=100) must fit the compiled NZ_STACK_MAX")
         if (allocated(error)) exit checks
      end block checks
   end subroutine test_nz_stack_guard

   subroutine test_cartesian_degrees(error)
      !! MOM6-style Cartesian sizing: len_lon/len_lat in degrees over nx/ny
      !! cells derive uniform dx/dy in metres via the arc-length formula
      !! (dx = rad_earth·len·π/180/n, no cos(lat)).  Mirrors the MOM6
      !! double_gyre: 22°×20° over 300×300 ⇒ dx≈8163.3, dy≈7421.1 m.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      real(wp) :: dx_exp, dy_exp

      call read_config_from_string( &
         "&sim_nml"//new_line('a')// &
         '   sim_type = "ocean"'//new_line('a')// &
         "/"//new_line('a')// &
         "&grid_nml"//new_line('a')// &
         "   nx = 300"//new_line('a')// &
         "   ny = 300"//new_line('a')// &
         "/"//new_line('a')// &
         "&ocean_grid_nml"//new_line('a')// &
         '   grid_config = "cartesian"'//new_line('a')// &
         '   axis_units = "degrees"'//new_line('a')// &
         "   len_lon = 22.0"//new_line('a')// &
         "   len_lat = 20.0"//new_line('a')// &
         "/"//new_line('a'), cfg)

      ! rad_earth defaults to 6.378e6 (MOM6); dx = rad_earth·len·DEG2RAD/n.
      dx_exp = 6.378e6_wp*22.0_wp*DEG2RAD/300.0_wp
      dy_exp = 6.378e6_wp*20.0_wp*DEG2RAD/300.0_wp
      call check(error, abs(cfg%dx - dx_exp) < 1.0e-6_wp, &
                 "cartesian degrees: dx should be "//" derived from len_lon")
      if (allocated(error)) return
      call check(error, abs(cfg%dy - dy_exp) < 1.0e-6_wp, &
                 "cartesian degrees: dy should be derived from len_lat")
      if (allocated(error)) return
      ! Sanity: the concrete MOM6 double_gyre numbers.
      call check(error, abs(cfg%dx - 8163.26_wp) < 0.1_wp, "dx ~ 8163.3 m")
      if (allocated(error)) return
      call check(error, abs(cfg%dy - 7421.14_wp) < 0.1_wp, "dy ~ 7421.1 m")
   end subroutine test_cartesian_degrees

   subroutine test_cartesian_meters_default(error)
      !! Default axis_units="meters" with len_lon unset leaves &grid_nml
      !! dx/dy untouched — bit-identical to the pre-feature path.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg

      call read_config_from_string( &
         "&sim_nml"//new_line('a')// &
         '   sim_type = "ocean"'//new_line('a')// &
         "/"//new_line('a')// &
         "&grid_nml"//new_line('a')// &
         "   nx = 300"//new_line('a')// &
         "   ny = 300"//new_line('a')// &
         "   dx = 8163.3"//new_line('a')// &
         "   dy = 7421.1"//new_line('a')// &
         "/"//new_line('a')// &
         "&ocean_grid_nml"//new_line('a')// &
         '   grid_config = "cartesian"'//new_line('a')// &
         "/"//new_line('a'), cfg)

      call check(error, abs(cfg%dx - 8163.3_wp) < 1.0e-9_wp, &
                 "meters default: dx untouched")
      if (allocated(error)) return
      call check(error, abs(cfg%dy - 7421.1_wp) < 1.0e-9_wp, &
                 "meters default: dy untouched")
   end subroutine test_cartesian_meters_default

   subroutine test_bt_halo_default_sentinel(error)
      !! The shipping default must be the AUTO sentinel (-1), not a concrete
      !! width — resolution to 0/8 happens later in the driver.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      call check(error, cfg%ocean%bt%bt_halo == BT_HALO_AUTO_SENTINEL, &
                 "bt_halo default is not the AUTO sentinel")
      if (allocated(error)) return
      call check(error, BT_HALO_AUTO_SENTINEL == -1, "AUTO sentinel changed from -1")
      if (allocated(error)) return
      call check(error, BT_HALO_AUTO_WIDTH == 8, "AUTO width changed from 8")
   end subroutine test_bt_halo_default_sentinel

   subroutine test_resolve_bt_halo(error)
      !! The sentinel-resolution truth table (compute_size / exclusion /
      !! explicit-value combinations).
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: AUTO = BT_HALO_AUTO_SENTINEL

      ! AUTO, serial => 0 (bit-identical to the historical default).
      call check(error, resolve_bt_halo(AUTO, 1, .false.) == 0, "auto serial clean should be 0")
      if (allocated(error)) return
      ! AUTO, multi-rank, clean => the validated width (8).
      call check(error, resolve_bt_halo(AUTO, 2, .false.) == BT_HALO_AUTO_WIDTH, &
                 "auto multi-rank clean should be 8")
      if (allocated(error)) return
      ! AUTO, multi-rank, but an exclusion is active => 0.
      call check(error, resolve_bt_halo(AUTO, 2, .true.) == 0, &
                 "auto multi-rank with exclusion should be 0")
      if (allocated(error)) return
      ! AUTO, serial, exclusion active => 0 (serial dominates anyway).
      call check(error, resolve_bt_halo(AUTO, 1, .true.) == 0, "auto serial excluded should be 0")
      if (allocated(error)) return
      ! Explicit 0 is never overridden (stays the v1 per-substep path).
      call check(error, resolve_bt_halo(0, 4, .false.) == 0, "explicit 0 must stay 0")
      if (allocated(error)) return
      ! Explicit width is left untouched — even on a serial run.
      call check(error, resolve_bt_halo(8, 1, .false.) == 8, "explicit 8 must stay 8 (serial)")
      if (allocated(error)) return
      call check(error, resolve_bt_halo(8, 2, .false.) == 8, "explicit 8 must stay 8 (multi)")
      if (allocated(error)) return
      ! Explicit width + exclusion is left untouched here; validate_config is
      ! the fail-loud gate for that impossible combination, not resolution.
      call check(error, resolve_bt_halo(4, 2, .true.) == 4, "explicit width unchanged by exclusion")
   end subroutine test_resolve_bt_halo

   subroutine test_bt_halo_auto_exclusion(error)
      !! The exclusion detector mirrors validate_config's explicit-bt_halo>0 set.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      logical :: excluded
      character(len=:), allocatable :: reason

      ! Clean default config => no exclusion.
      call bt_halo_auto_exclusion(cfg, excluded, reason)
      call check(error,.not. excluded, "clean config wrongly flagged as excluded")
      if (allocated(error)) return
      call check(error, len_trim(reason) == 0, "clean config should have empty reason")
      if (allocated(error)) return

      ! Each excluded feature must flip the flag.
      cfg%ocean%wetdry%enable = .true.
      call bt_halo_auto_exclusion(cfg, excluded, reason)
      call check(error, excluded, "wetdry enable not detected as exclusion")
      if (allocated(error)) return
      cfg%ocean%wetdry%enable = .false.

      cfg%ocean%tides%enable = .true.
      call bt_halo_auto_exclusion(cfg, excluded, reason)
      call check(error, excluded, "tides enable not detected as exclusion")
      if (allocated(error)) return
      cfg%ocean%tides%enable = .false.

      ! Surface-pressure loading rides the SAME `eta_forcing` seam as the
      ! body tide, and the wide-halo BT clone carries no copy of it — so
      ! `run_stage_split` `error stop`s on `bt_halo > 0` with an
      ! `eta_forcing` actual.  validate_config only polices that under an
      ! EXPLICIT `bt_halo > 0`; AUTO must resolve to 0 here or a default
      ! multi-rank psurf run passes validation and dies at the runtime
      ! backstop.
      cfg%ocean%psurf%enable = .true.
      call bt_halo_auto_exclusion(cfg, excluded, reason)
      call check(error, excluded, "psurf enable not detected as exclusion")
      if (allocated(error)) return
      call check(error, resolve_bt_halo(BT_HALO_AUTO_SENTINEL, 4, excluded) == 0, &
                 "psurf + auto bt_halo on 4 ranks must resolve to 0, not the "// &
                 "wide-halo march-in")
      if (allocated(error)) return
      cfg%ocean%psurf%enable = .false.

      ! Porous barriers: `bt_wide` re-fills its own `metrics_w` from the
      ! grid formula and nothing gives it the porous statistics, so the
      ! wide fast loop would silently transport on the UN-narrowed
      ! `dy_cu`/`dx_cv` — an answer that is neither the porous one nor the
      ! baseline.  AUTO must therefore resolve to 0 with porous on.
      cfg%ocean%porous%enable = .true.
      call bt_halo_auto_exclusion(cfg, excluded, reason)
      call check(error, excluded, "porous enable not detected as exclusion")
      if (allocated(error)) return
      call check(error, resolve_bt_halo(BT_HALO_AUTO_SENTINEL, 4, excluded) == 0, &
                 "porous + auto bt_halo on 4 ranks must resolve to 0, not the "// &
                 "wide-halo march-in")
      if (allocated(error)) return
      cfg%ocean%porous%enable = .false.

      ! Ice-shelf cavity: `bt_wide`'s own `metrics_w` is re-filled from the
      ! grid formula and carries no `z_draft`, so the wide fast loop would
      ! take the BED as its reference depth — a 1000 m ocean where the
      ! cavity has 500 m of water under 500 m of ice.  AUTO must resolve
      ! to 0 rather than manufacture a width that then trips the explicit
      ! abort the user never asked for.
      cfg%ocean%cavity_dyn%enable = .true.
      call bt_halo_auto_exclusion(cfg, excluded, reason)
      call check(error, excluded, "cavity_dyn enable not detected as exclusion")
      if (allocated(error)) return
      call check(error, resolve_bt_halo(BT_HALO_AUTO_SENTINEL, 4, excluded) == 0, &
                 "cavity + auto bt_halo on 4 ranks must resolve to 0, not the "// &
                 "wide-halo march-in")
      if (allocated(error)) return
      cfg%ocean%cavity_dyn%enable = .false.

      cfg%ocean%grid%grid_config = "tripolar"
      call bt_halo_auto_exclusion(cfg, excluded, reason)
      call check(error, excluded, "tripolar grid_config not detected as exclusion")
   end subroutine test_bt_halo_auto_exclusion

   subroutine test_ice_hlim_count_leading_run(error)
      !! PR-58 case 7 (PLAN_PR58_ice_hlim_override.md §9): the count
      !! derivation, including the non-contiguous case that feeds case 8
      !! (c). The `0` case matters most — it is the default-off dispatch
      !! gate (`ocean_state_init_from_config`'s latch reads
      !! `ice_hlim_count(...) > 0`).
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: hlim(MAX_ICE_HLIM_VALS)

      ! Sentinel-filled => 0.
      hlim = -1.0_wp
      call check(error, ice_hlim_count(hlim) == 0, "all-sentinel hlim must count 0")
      if (allocated(error)) return

      ! [0.1, 0.3, -1, -1, ...] => 2.
      hlim = -1.0_wp
      hlim(1) = 0.1_wp
      hlim(2) = 0.3_wp
      call check(error, ice_hlim_count(hlim) == 2, "leading run of 2 must count 2")
      if (allocated(error)) return

      ! [0.1, -5.0, 0.3, -1, ...] => 1 (leading run only, non-contiguous
      ! trailing value dropped from the count by design).
      hlim = -1.0_wp
      hlim(1) = 0.1_wp
      hlim(2) = -5.0_wp
      hlim(3) = 0.3_wp
      call check(error, ice_hlim_count(hlim) == 1, "non-contiguous list must count only the leading run")
   end subroutine test_ice_hlim_count_leading_run

   subroutine test_ice_hlim_spec_validity(error)
      !! PR-58 case 8 (PLAN_PR58_ice_hlim_override.md §9): each rejected
      !! case is a real failure downstream, not hygiene — (d)/(e)/(f)
      !! are what keep every `mh_lim(c+1)` divisor in `rdb_ice_itd`/
      !! `rdb_ice_transport` positive and the ladder ordered. Testing the
      !! predicate directly is the only way to reach these:
      !! `validate_config` `error stop`s (test_ocean_ice_transport.F90's
      !! stated constraint for this exact envelope class).
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: hlim(MAX_ICE_HLIM_VALS)
      logical :: ok
      character(len=:), allocatable :: reason
      integer, parameter :: NCAT = 5

      ! Valid: [0.05, 0.5, 1.2, -1, ...] at ncat=5.
      hlim = -1.0_wp
      hlim(1) = 0.05_wp
      hlim(2) = 0.5_wp
      hlim(3) = 1.2_wp
      call ice_hlim_spec_is_valid(hlim, NCAT, ok, reason)
      call check(error, ok, "a valid 3-entry list at ncat=5 must pass: "//reason)
      if (allocated(error)) return

      ! (a) n=1: [0.5, -1, ...] -- D1, the deliberate divergence from
      ! SIS2's silent fallback.
      hlim = -1.0_wp
      hlim(1) = 0.5_wp
      call ice_hlim_spec_is_valid(hlim, NCAT, ok, reason)
      call check(error,.not. ok, "a 1-element hlim must be rejected (D1)")
      if (allocated(error)) return

      ! (b) n > ncat+1: 7 entries at ncat=5 (ncat+1=6).
      hlim = -1.0_wp
      hlim(1:7) = [0.05_wp, 0.2_wp, 0.4_wp, 0.6_wp, 0.8_wp, 1.0_wp, 1.2_wp]
      call ice_hlim_spec_is_valid(hlim, NCAT, ok, reason)
      call check(error,.not. ok, "more than ncat+1 entries must be rejected")
      if (allocated(error)) return

      ! (c) non-contiguous: [0.1, -5.0, 0.3, ...].
      hlim = -1.0_wp
      hlim(1) = 0.1_wp
      hlim(2) = -5.0_wp
      hlim(3) = 0.3_wp
      call ice_hlim_spec_is_valid(hlim, NCAT, ok, reason)
      call check(error,.not. ok, "a non-contiguous list must be rejected")
      if (allocated(error)) return

      ! (d) hlim(1) == 0.0.
      hlim = -1.0_wp
      hlim(1) = 0.0_wp
      hlim(2) = 0.5_wp
      call ice_hlim_spec_is_valid(hlim, NCAT, ok, reason)
      call check(error,.not. ok, "hlim(1) == 0 must be rejected")
      if (allocated(error)) return

      ! (e) non-increasing: [0.1, 0.5, 0.4, ...].
      hlim = -1.0_wp
      hlim(1) = 0.1_wp
      hlim(2) = 0.5_wp
      hlim(3) = 0.4_wp
      call ice_hlim_spec_is_valid(hlim, NCAT, ok, reason)
      call check(error,.not. ok, "a non-increasing list must be rejected")
      if (allocated(error)) return

      ! (f) repeated edge: [0.1, 0.5, 0.5, ...] -- STRICT, not >=.
      hlim = -1.0_wp
      hlim(1) = 0.1_wp
      hlim(2) = 0.5_wp
      hlim(3) = 0.5_wp
      call ice_hlim_spec_is_valid(hlim, NCAT, ok, reason)
      call check(error,.not. ok, "a repeated edge must be rejected (strict, not >=)")
   end subroutine test_ice_hlim_spec_validity

   subroutine test_read_ocean_subnmls(error)
      !! Round-trip the 11 ocean sub-namelists through read_config and
      !! assert each `cfg%ocean%group%field` parsed.  Locks in the
      !! PR-B2 split + prefix-drop.  In particular the three distinct
      !! `form` keys (coriolis / pgf / bdrag) verify the per-sub-nml
      !! local-var scoping doesn't leak across groups.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      integer :: io_unit

      open (newunit=io_unit, file="test_ocean.nml", status="replace", action="write")
      write (io_unit, '(a)') "&sim_nml"
      write (io_unit, '(a)') '   sim_type = "ocean"'
      write (io_unit, '(a)') "/"
      write (io_unit, '(a)') "&ocean_coriolis_nml"
      write (io_unit, '(a)') '   form = "sadourny_energy"'
      write (io_unit, '(a)') "/"
      write (io_unit, '(a)') "&ocean_thermo_nml"
      write (io_unit, '(a)') "   enable_thermodynamics = .false."
      write (io_unit, '(a)') "/"
      write (io_unit, '(a)') "&ocean_tracers_nml"
      write (io_unit, '(a)') "   enable_ideal_age = .true."
      write (io_unit, '(a)') "/"
      write (io_unit, '(a)') "&ocean_bt_nml"
      write (io_unit, '(a)') "   n_inner = 60"
      write (io_unit, '(a)') "   use_cont_type = .true."
      write (io_unit, '(a)') "   correction_h_weighted = .true."
      write (io_unit, '(a)') "/"
      write (io_unit, '(a)') "&ocean_pgf_nml"
      write (io_unit, '(a)') '   form = "fv_mom6"'
      write (io_unit, '(a)') "   gfs_scale = 0.98"
      write (io_unit, '(a)') "/"
      write (io_unit, '(a)') "&ocean_bdrag_nml"
      write (io_unit, '(a)') '   form = "linear"'
      write (io_unit, '(a)') "   r = 2.5e-5"
      write (io_unit, '(a)') "   hbbl = 10.0"
      write (io_unit, '(a)') "/"
      write (io_unit, '(a)') "&ocean_hvisc_nml"
      write (io_unit, '(a)') "   nu_h = 10000.0"
      write (io_unit, '(a)') "   smag_ah = .true."
      write (io_unit, '(a)') "/"
      write (io_unit, '(a)') "&ocean_vmix_nml"
      write (io_unit, '(a)') "   use_closure = .false."
      write (io_unit, '(a)') "   dt_therm_ratio = 2"
      write (io_unit, '(a)') "/"
      write (io_unit, '(a)') "&ocean_continuity_nml"
      write (io_unit, '(a)') "   ppm_limit_pos = .true."
      write (io_unit, '(a)') "   h_min = 2.0"
      write (io_unit, '(a)') "/"
      write (io_unit, '(a)') "&ocean_topo_nml"
      write (io_unit, '(a)') '   topo_config = "spoon"'
      write (io_unit, '(a)') "   max_depth = 4000.0"
      write (io_unit, '(a)') "   coriolis_beta = 1.76e-11"
      write (io_unit, '(a)') "/"
      write (io_unit, '(a)') "&ocean_ic_nml"
      write (io_unit, '(a)') '   ic_config = "eady"'
      write (io_unit, '(a)') "   alpha_T = 2.0e-4"
      write (io_unit, '(a)') "/"
      write (io_unit, '(a)') "&ocean_diag_nml"
      write (io_unit, '(a)') "   enabled = .true."
      write (io_unit, '(a)') '   vgrid = "z_fixed"'
      write (io_unit, '(a)') "/"
      close (io_unit)

      call wait_for_file("test_ocean.nml")
      call read_config("test_ocean.nml", cfg)

      ! coriolis / pgf / bdrag all have a `form` key — the three must
      ! land in their own sub-structs without cross-contamination.
      call check(error, trim(cfg%ocean%coriolis%form) == "sadourny_energy", &
                 "ocean%coriolis%form should be sadourny_energy")
      if (allocated(error)) return
      call check(error, trim(cfg%ocean%pgf%form) == "fv_mom6", &
                 "ocean%pgf%form should be fv_mom6")
      if (allocated(error)) return
      call check(error, trim(cfg%ocean%bdrag%form) == "linear", &
                 "ocean%bdrag%form should be linear")
      if (allocated(error)) return
      ! thermo
      call check(error,.not. cfg%ocean%thermo%enable_thermodynamics, &
                 "ocean%thermo%enable_thermodynamics should be .false.")
      if (allocated(error)) return
      ! tracers
      call check(error, cfg%ocean%tracers%enable_ideal_age, &
                 "ocean%tracers%enable_ideal_age should be .true.")
      if (allocated(error)) return
      ! bt
      call check(error, cfg%ocean%bt%n_inner == 60, "ocean%bt%n_inner should be 60")
      if (allocated(error)) return
      call check(error, cfg%ocean%bt%use_cont_type, "ocean%bt%use_cont_type should be .true.")
      if (allocated(error)) return
      call check(error, cfg%ocean%bt%correction_h_weighted, &
                 "ocean%bt%correction_h_weighted should be .true.")
      if (allocated(error)) return
      ! pgf / bdrag scalars
      call check(error, abs(cfg%ocean%pgf%gfs_scale - 0.98_wp) < 1.0e-12_wp, &
                 "ocean%pgf%gfs_scale should be 0.98")
      if (allocated(error)) return
      call check(error, abs(cfg%ocean%bdrag%r - 2.5e-5_wp) < 1.0e-12_wp, &
                 "ocean%bdrag%r should be 2.5e-5")
      if (allocated(error)) return
      call check(error, abs(cfg%ocean%bdrag%hbbl - 10.0_wp) < 1.0e-12_wp, &
                 "ocean%bdrag%hbbl should be 10.0")
      if (allocated(error)) return
      ! hvisc
      call check(error, abs(cfg%ocean%hvisc%nu_h - 10000.0_wp) < 1.0e-9_wp, &
                 "ocean%hvisc%nu_h should be 10000")
      if (allocated(error)) return
      call check(error, cfg%ocean%hvisc%smag_ah, "ocean%hvisc%smag_ah should be .true.")
      if (allocated(error)) return
      ! vmix
      call check(error,.not. cfg%ocean%vmix%use_closure, &
                 "ocean%vmix%use_closure should be .false.")
      if (allocated(error)) return
      call check(error, cfg%ocean%vmix%dt_therm_ratio == 2, &
                 "ocean%vmix%dt_therm_ratio should be 2")
      if (allocated(error)) return
      ! continuity
      call check(error, cfg%ocean%continuity%ppm_limit_pos, &
                 "ocean%continuity%ppm_limit_pos should be .true.")
      if (allocated(error)) return
      call check(error, abs(cfg%ocean%continuity%h_min - 2.0_wp) < 1.0e-12_wp, &
                 "ocean%continuity%h_min should be 2.0")
      if (allocated(error)) return
      ! topo
      call check(error, trim(cfg%ocean%topo%topo_config) == "spoon", &
                 "ocean%topo%topo_config should be spoon")
      if (allocated(error)) return
      call check(error, abs(cfg%ocean%topo%max_depth - 4000.0_wp) < 1.0e-9_wp, &
                 "ocean%topo%max_depth should be 4000")
      if (allocated(error)) return
      call check(error, abs(cfg%ocean%topo%coriolis_beta - 1.76e-11_wp) < 1.0e-20_wp, &
                 "ocean%topo%coriolis_beta should be 1.76e-11")
      if (allocated(error)) return
      ! ic
      call check(error, trim(cfg%ocean%ic%ic_config) == "eady", &
                 "ocean%ic%ic_config should be eady")
      if (allocated(error)) return
      call check(error, abs(cfg%ocean%ic%alpha_T - 2.0e-4_wp) < 1.0e-12_wp, &
                 "ocean%ic%alpha_T should be 2.0e-4")
      if (allocated(error)) return
      ! diag
      call check(error, cfg%ocean%diag%enabled, "ocean%diag%enabled should be .true.")
      if (allocated(error)) return
      call check(error, trim(cfg%ocean%diag%vgrid) == "z_fixed", &
                 "ocean%diag%vgrid should be z_fixed")
      if (allocated(error)) return

      ! A knob NOT set in the file keeps its default (bdrag%bed_factor = 1.0).
      call check(error, abs(cfg%ocean%bdrag%bed_factor - 1.0_wp) < 1.0e-12_wp, &
                 "unset ocean%bdrag%bed_factor should keep default 1.0")

      open (newunit=io_unit, file="test_ocean.nml", status="old")
      close (io_unit, status="delete")
   end subroutine test_read_ocean_subnmls

   subroutine test_ocean_diag_time_unit(error)
      !! `ocean%diag%dt_out` must pick up the `&time_nml time_unit`
      !! multiplier (regression: PR-B2 initially left the cascade
      !! operating on a dead local instead of cfg%ocean%diag%dt_out).
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      integer :: io_unit

      open (newunit=io_unit, file="test_diag_tu.nml", status="replace", action="write")
      write (io_unit, '(a)') "&sim_nml"
      write (io_unit, '(a)') '   sim_type = "ocean"'
      write (io_unit, '(a)') "/"
      write (io_unit, '(a)') "&time_nml"
      write (io_unit, '(a)') "   t_end     = 1.0"
      write (io_unit, '(a)') '   time_unit = "day"'
      write (io_unit, '(a)') "/"
      write (io_unit, '(a)') "&ocean_diag_nml"
      write (io_unit, '(a)') "   dt_out = 0.25"   ! quarter-day → 21600 s
      write (io_unit, '(a)') "/"
      close (io_unit)

      call wait_for_file("test_diag_tu.nml")
      call read_config("test_diag_tu.nml", cfg)

      ! 0.25 day × 86400 s/day = 21600 s.
      call check(error, abs(cfg%ocean%diag%dt_out - 21600.0_wp) < 1.0e-6_wp, &
                 "ocean%diag%dt_out (0.25 day) should convert to 21600 s")

      open (newunit=io_unit, file="test_diag_tu.nml", status="old")
      close (io_unit, status="delete")
   end subroutine test_ocean_diag_time_unit

   subroutine test_defaults(error)
      !! Default config_t values should be sensible
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg

      call check(error, cfg%nx == 200, "default nx should be 200")
      if (allocated(error)) return
      call check(error, cfg%ny == 1, "default ny should be 1")
      if (allocated(error)) return
      call check(error, abs(cfg%cfl - 0.45_wp) < 1.0e-15_wp, "default cfl should be 0.45")
      if (allocated(error)) return
      call check(error, cfg%nghost == 3, "default nghost should be 3")
      if (allocated(error)) return
      call check(error, trim(cfg%bc_west) == "wall", "default bc_west should be wall")
   end subroutine test_defaults

   subroutine test_read_grid(error)
      !! Reading a namelist file should override grid defaults
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      integer :: io_unit

      ! Write a temporary namelist file
      open (newunit=io_unit, file="test_grid.nml", status="replace", action="write")
      write (io_unit, '(a)') "&grid_nml"
      write (io_unit, '(a)') "  nx = 500"
      write (io_unit, '(a)') "  ny = 100"
      write (io_unit, '(a)') "  dx = 50.0"
      write (io_unit, '(a)') "  dy = 25.0"
      write (io_unit, '(a)') "/"
      close (io_unit)
      call wait_for_file("test_grid.nml")

      call read_config("test_grid.nml", cfg)

      call check(error, cfg%nx == 500, "nx should be 500")
      if (allocated(error)) return
      call check(error, cfg%ny == 100, "ny should be 100")
      if (allocated(error)) return
      call check(error, abs(cfg%dx - 50.0_wp) < 1.0e-12_wp, "dx should be 50.0")
      if (allocated(error)) return
      call check(error, abs(cfg%dy - 25.0_wp) < 1.0e-12_wp, "dy should be 25.0")
      if (allocated(error)) return
      ! Unset values should keep defaults
      call check(error, abs(cfg%cfl - 0.45_wp) < 1.0e-15_wp, &
                 "cfl should keep default 0.45")

      ! Clean up
      open (newunit=io_unit, file="test_grid.nml", status="old")
      close (io_unit, status="delete")
   end subroutine test_read_grid

   subroutine test_read_all(error)
      !! Reading all namelist groups should populate all fields
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      integer :: io_unit

      open (newunit=io_unit, file="test_all.nml", status="replace", action="write")
      write (io_unit, '(a)') "&grid_nml"
      write (io_unit, '(a)') "  nx = 300"
      write (io_unit, '(a)') "  ny = 200"
      write (io_unit, '(a)') "  dx = 10.0"
      write (io_unit, '(a)') "  dy = 10.0"
      write (io_unit, '(a)') "/"
      write (io_unit, '(a)') "&time_nml"
      write (io_unit, '(a)') "  t_end = 86400.0"
      write (io_unit, '(a)') "  cfl = 0.3"
      write (io_unit, '(a)') "  dt_max = 5.0"
      write (io_unit, '(a)') "/"
      write (io_unit, '(a)') "&physics_nml"
      write (io_unit, '(a)') "  manning_n = 0.025"
      write (io_unit, '(a)') "  wind_stress_x = 0.1"
      write (io_unit, '(a)') "/"
      write (io_unit, '(a)') "&boundary_nml"
      write (io_unit, '(a)') "  bc_west = 'tidal'"
      write (io_unit, '(a)') "  bc_east = 'open'"
      write (io_unit, '(a)') "/"
      close (io_unit)
      call wait_for_file("test_all.nml")

      call read_config("test_all.nml", cfg)

      call check(error, cfg%nx == 300, "nx should be 300")
      if (allocated(error)) return
      call check(error, abs(cfg%t_end - 86400.0_wp) < 1.0e-10_wp, &
                 "t_end should be 86400")
      if (allocated(error)) return
      call check(error, abs(cfg%cfl - 0.3_wp) < 1.0e-15_wp, "cfl should be 0.3")
      if (allocated(error)) return
      call check(error, abs(cfg%dt_max - 5.0_wp) < 1.0e-15_wp, "dt_max should be 5.0")
      if (allocated(error)) return
      call check(error, abs(cfg%manning_n - 0.025_wp) < 1.0e-15_wp, &
                 "manning_n should be 0.025")
      if (allocated(error)) return
      call check(error, abs(cfg%wind_stress_x - 0.1_wp) < 1.0e-15_wp, &
                 "wind_stress_x should be 0.1")
      if (allocated(error)) return
      call check(error, trim(cfg%bc_west) == "tidal", "bc_west should be tidal")
      if (allocated(error)) return
      call check(error, trim(cfg%bc_east) == "open", "bc_east should be open")

      open (newunit=io_unit, file="test_all.nml", status="old")
      close (io_unit, status="delete")
   end subroutine test_read_all

   subroutine test_missing_file(error)
      !! Missing file should not crash, defaults should be preserved
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg

      call read_config("nonexistent_file_12345.nml", cfg)

      call check(error, cfg%nx == 200, "nx should keep default on missing file")
      if (allocated(error)) return
      call check(error, abs(cfg%cfl - 0.45_wp) < 1.0e-15_wp, &
                 "cfl should keep default on missing file")
   end subroutine test_missing_file

   subroutine test_diag_density_requires_rho_levels(error)
      !! PR-9 §9.5 / §5B guard.  Density bins have no auto-fill (unlike
      !! sigma/z*), so the config-time predicate `diag_density_levels_ok`
      !! must reject vgrid='density' (or a ':density'/':rho' diags entry)
      !! with no/non-monotone/oversized rho_levels, and accept a valid
      !! increasing list.  `validate_config` calls this predicate and
      !! `error stop`s on failure — tested here at the predicate level
      !! (the house pattern: `si_bc_supported`/`lateral_closure_is_implemented`)
      !! rather than by exercising the fatal `error stop` path.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: levels(MAX_OCEAN_DIAG_Z_LEVELS)
      real(wp) :: oversized(MAX_OCEAN_DIAG_Z_LEVELS + 1)

      checks: block
         ! Density not requested at all -> always ok, regardless of levels.
         levels = -1.0_wp
         call check(error, diag_density_levels_ok("layer", "", 0, levels), &
                    "density not requested => no constraint")
         if (allocated(error)) exit checks

         ! vgrid='density' with n_rho_levels = 0 -> reject.
         call check(error,.not. diag_density_levels_ok("density", "", 0, levels), &
                    "vgrid='density' with n_rho_levels=0 must be rejected")
         if (allocated(error)) exit checks

         ! Per-diagnostic ':density' attribute (vgrid stays 'layer') with
         ! n_rho_levels = 0 -> also reject.
         call check(error,.not. diag_density_levels_ok("layer", "temperature:density", 0, levels), &
                    "diags ':density' with n_rho_levels=0 must be rejected")
         if (allocated(error)) exit checks

         ! ':rho' spelling triggers the same requirement.
         call check(error,.not. diag_density_levels_ok("layer", "salinity:rho", 0, levels), &
                    "diags ':rho' with n_rho_levels=0 must be rejected")
         if (allocated(error)) exit checks

         ! Non-monotone rho_levels -> reject.
         levels(1:3) = [1024.0_wp, 1023.0_wp, 1025.0_wp]
         call check(error,.not. diag_density_levels_ok("density", "", 3, levels), &
                    "non-monotone rho_levels must be rejected")
         if (allocated(error)) exit checks

         ! Oversized n_rho_levels (> MAX_OCEAN_DIAG_Z_LEVELS) -> reject.
         oversized = 0.0_wp
         oversized(1) = 1020.0_wp
         call check(error,.not. diag_density_levels_ok("density", "", &
                                                       MAX_OCEAN_DIAG_Z_LEVELS + 1, oversized), &
                    "n_rho_levels > MAX_OCEAN_DIAG_Z_LEVELS must be rejected")
         if (allocated(error)) exit checks

         ! Valid strictly-increasing list -> accept.
         levels(1:4) = [1021.0_wp, 1022.5_wp, 1024.0_wp, 1027.0_wp]
         call check(error, diag_density_levels_ok("density", "", 4, levels), &
                    "valid strictly-increasing rho_levels must be accepted")
      end block checks
   end subroutine test_diag_density_requires_rho_levels

end module test_config
