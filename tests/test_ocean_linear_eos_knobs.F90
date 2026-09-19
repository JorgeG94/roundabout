!! Analytical + contract tests for the `&ocean_ic_nml` linear-EOS
!! reference state (`alpha_T`, `beta_S`, `T_ref`, `S_ref`, `rho_0`).
module test_ocean_linear_eos_knobs
   !! Before this capability only `alpha_T` and `rho_0` reached the ocean
   !! path's linear EOS; `beta_S`, `T_ref` and `S_ref` were fixed
   !! `eos_t` component defaults that no namelist could move, and the
   !! coastal-legacy `&tracer_nml` spellings of all four validated and
   !! then did nothing.  This suite pins all three halves of the fix:
   !!
   !!   (a) **Configure path** — a namelist string carrying the full
   !!       quintet lands on `ocean_state%eos` through the PRODUCTION
   !!       route (`read_config_from_string` -> `init_from_config`), and
   !!       `eos_density_point` then reproduces a hand-computed density.
   !!       The case is the ISOMIP+ protocol reference state
   !!       (Asay-Davis et al. 2016) because it is the one that exposes
   !!       the units trap below.
   !!
   !!   (b) **Units convention** — Roundabout's linear EOS is the
   !!       DENSITY-ANOMALY form
   !!
   !!         rho = rho_0 + beta_S*(S - S_ref) - alpha_T*(T - T_ref)
   !!
   !!       so `alpha_T`/`beta_S` are DIMENSIONAL (kg/m^3 per degC /
   !!       per PSU).  Protocols quote the FRACTIONAL coefficients of
   !!       the equivalent
   !!
   !!         rho = rho_0*(1 - alpha*(T - T_ref) + beta*(S - S_ref))
   !!
   !!       so `alpha_T = rho_0*alpha`, `beta_S = rho_0*beta`.  The
   !!       density check is written against the FRACTIONAL form, i.e.
   !!       independently of the expression the kernel evaluates, so it
   !!       fails if the conversion documented on the knobs is wrong.
   !!
   !!   (c) **Retirement** — the coastal-legacy `&tracer_nml` quartet
   !!       (`alpha_T`, `beta_S`, `T_ref`, `S_ref`) reached
   !!       `tracer_t%eos_coeff`/`eos_ref`, which no ocean kernel reads.
   !!       Moving any of them off its historical default is now a
   !!       fail-loud `validate_config` error naming the live
   !!       `&ocean_ic_nml` replacement, instead of being accepted and
   !!       silently ignored.
   !!
   !!   (d) **Bit-identity** — with `&ocean_ic_nml` untouched, all five
   !!       EOS members equal the `eos_t` component defaults they had
   !!       before the knobs existed, so every shipped namelist is
   !!       unchanged.
   !!
   !! Host-only by construction: `eos_t` is a flat POD passed BY VALUE
   !! into the `!$acc routine seq` point routines, and the configure-time
   !! assignment under test runs before `ocean_state_enter_data`, so
   !! there is no device-resident array in the loop and no `mem:separate`
   !! mapping to arrange.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_config, only: config_t, read_config_from_string, validate_config
   use rdb_grid, only: hgrid_t
   use rdb_ocean_state, only: ocean_state_t
   use rdb_eos, only: eos_t, eos_density_point
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_CONFIG_VALIDATE
   implicit none
   private

   public :: collect_ocean_linear_eos_knobs_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NX_PHYS = 6
   integer, parameter :: NY_PHYS = 5

   ! ISOMIP+ (Asay-Davis et al. 2016) linear-EOS reference state, in the
   ! protocol's own FRACTIONAL units.  Test inputs only — the numbers are
   ! carried here to exercise the conversion, not to certify the protocol.
   real(wp), parameter :: ISOMIP_ALPHA_FRAC = 3.733e-5_wp   ! 1/degC
   real(wp), parameter :: ISOMIP_BETA_FRAC = 7.843e-4_wp    ! 1/PSU
   real(wp), parameter :: ISOMIP_T_REF = -1.0_wp            ! degC
   real(wp), parameter :: ISOMIP_S_REF = 34.2_wp            ! PSU
   real(wp), parameter :: ISOMIP_RHO_0 = 1027.51_wp         ! kg/m^3

   ! The DIMENSIONAL knob values `&ocean_ic_nml` actually takes.
   real(wp), parameter :: ISOMIP_ALPHA_T = ISOMIP_RHO_0*ISOMIP_ALPHA_FRAC
   real(wp), parameter :: ISOMIP_BETA_S = ISOMIP_RHO_0*ISOMIP_BETA_FRAC

   ! Sample water mass: a shelf-water column a little warmer and a little
   ! saltier than the reference point, so BOTH anomaly terms are live and
   ! carry opposite signs (a sign slip cannot cancel).
   real(wp), parameter :: SAMPLE_T = 0.5_wp
   real(wp), parameter :: SAMPLE_S = 34.5_wp

   ! Hand-computed density for (SAMPLE_T, SAMPLE_S):
   !   alpha_T = 1027.51 * 3.733e-5  = 0.0383569483 kg/m^3/degC
   !   beta_S  = 1027.51 * 7.843e-4  = 0.805876093  kg/m^3/PSU
   !   dT = 0.5 - (-1.0) = 1.5 ;  dS = 34.5 - 34.2 = 0.3
   !   rho = 1027.51 + 0.805876093*0.3 - 0.0383569483*1.5
   !       = 1027.51 + 0.2417628279 - 0.0575354225
   !       = 1027.6942274055
   real(wp), parameter :: SAMPLE_RHO_HAND = 1027.6942274055_wp

   ! `eos_t` component defaults as they stood before the knobs existed —
   ! the bit-identity contract for every namelist that stays silent.
   real(wp), parameter :: LEGACY_EOS_ALPHA_T = 1.7e-4_wp
   real(wp), parameter :: LEGACY_EOS_BETA_S = 7.6e-4_wp
   real(wp), parameter :: LEGACY_EOS_T_REF = 10.0_wp
   real(wp), parameter :: LEGACY_EOS_S_REF = 35.0_wp
   real(wp), parameter :: LEGACY_EOS_RHO0 = 1035.0_wp

contains

   subroutine collect_ocean_linear_eos_knobs_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("linear_eos_isomip_reference_state", test_isomip_reference_state), &
                  new_unittest("linear_eos_defaults_bit_identical", test_defaults_bit_identical), &
                  new_unittest("linear_eos_tracer_alpha_T_retired", test_tracer_alpha_t_retired), &
                  new_unittest("linear_eos_tracer_beta_S_retired", test_tracer_beta_s_retired), &
                  new_unittest("linear_eos_tracer_T_ref_retired", test_tracer_t_ref_retired), &
                  new_unittest("linear_eos_tracer_S_ref_retired", test_tracer_s_ref_retired) &
                  ]
   end subroutine collect_ocean_linear_eos_knobs_tests

   subroutine build_state(cfg, grid, state, nml, ierr)
      !! Drive the PRODUCTION configure path: parse the namelist text
      !! exactly as `rdb_ocean_create` / the driver does, then hand it to
      !! `ocean_state_t%init_from_config`, which is the single site that
      !! copies `&ocean_ic_nml` onto `state%eos`.
      type(config_t), intent(out) :: cfg
      type(hgrid_t), intent(out) :: grid
      type(ocean_state_t), intent(out) :: state
      character(len=*), intent(in) :: nml
      integer, intent(out) :: ierr

      call read_config_from_string(nml, cfg, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call grid%init(NX_PHYS, NY_PHYS, NGHOST, 1000.0_wp, 1000.0_wp)
      call state%init_from_config(cfg, grid)
   end subroutine build_state

   subroutine test_isomip_reference_state(error)
      !! (a) + (b): the ISOMIP+ quintet survives the configure path and
      !! the linear EOS reproduces the hand-computed density.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      character(len=:), allocatable :: nml
      real(wp) :: rho, rho_fractional_form
      integer :: ierr

      ! DIMENSIONAL values: alpha_T = rho_0*alpha, beta_S = rho_0*beta.
      nml = "&sim_nml sim_type = 'ocean' /"//new_line("a")// &
            "&grid_nml nx = 6, ny = 5, nghost = 2, dx = 1000.0, dy = 1000.0 /"//new_line("a")// &
            "&ocean_eos_nml eos = 'linear' /"//new_line("a")// &
            "&ocean_ic_nml"//new_line("a")// &
            "   alpha_T = 3.83569483e-2"//new_line("a")// &
            "   beta_S  = 8.05876093e-1"//new_line("a")// &
            "   T_ref   = -1.0"//new_line("a")// &
            "   S_ref   = 34.2"//new_line("a")// &
            "   rho_0   = 1027.51"//new_line("a")// &
            "/"//new_line("a")

      call build_state(cfg, grid, state, nml, ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "ISOMIP+ namelist must parse")
      if (allocated(error)) return

      checks: block
         ! Every member reached the EOS handle — not just alpha_T/rho_0.
         call check(error, abs(state%eos%rho0 - ISOMIP_RHO_0) < 1.0e-12_wp, &
                    "&ocean_ic_nml rho_0 must reach eos%rho0")
         if (allocated(error)) exit checks
         call check(error, abs(state%eos%T_ref - ISOMIP_T_REF) < 1.0e-12_wp, &
                    "&ocean_ic_nml T_ref must reach eos%T_ref")
         if (allocated(error)) exit checks
         call check(error, abs(state%eos%S_ref - ISOMIP_S_REF) < 1.0e-12_wp, &
                    "&ocean_ic_nml S_ref must reach eos%S_ref")
         if (allocated(error)) exit checks
         call check(error, abs(state%eos%alpha_T - ISOMIP_ALPHA_T) < 1.0e-9_wp, &
                    "&ocean_ic_nml alpha_T must reach eos%alpha_T")
         if (allocated(error)) exit checks
         call check(error, abs(state%eos%beta_S - ISOMIP_BETA_S) < 1.0e-9_wp, &
                    "&ocean_ic_nml beta_S must reach eos%beta_S")
         if (allocated(error)) exit checks

         rho = eos_density_point(state%eos, SAMPLE_T, SAMPLE_S, 0.0_wp)

         ! Hand-computed literal (see the parameter's derivation above).
         call check(error, abs(rho - SAMPLE_RHO_HAND) < 1.0e-8_wp, &
                    "linear EOS must reproduce the hand-computed ISOMIP+ density")
         if (allocated(error)) exit checks

         ! Same number via the protocol's OWN fractional form.  This is
         ! the units contract: if the documented `alpha_T = rho_0*alpha`
         ! conversion were wrong, these two would not agree.
         rho_fractional_form = ISOMIP_RHO_0*(1.0_wp &
                                             - ISOMIP_ALPHA_FRAC*(SAMPLE_T - ISOMIP_T_REF) &
                                             + ISOMIP_BETA_FRAC*(SAMPLE_S - ISOMIP_S_REF))
         call check(error, abs(rho - rho_fractional_form) < 1.0e-8_wp, &
                    "dimensional knobs must equal rho_0 * the fractional coefficients")
         if (allocated(error)) exit checks

         ! Sanity: the salinity term dominates here, so a run that fed
         ! the FRACTIONAL numbers straight in (the trap the knob doc
         ! warns about) would be ~1000x too weak.  Assert we are not.
         call check(error, rho - ISOMIP_RHO_0 > 0.1_wp, &
                    "density anomaly must be the dimensional-magnitude one")
      end block checks

      call state%destroy()
   end subroutine test_isomip_reference_state

   subroutine test_defaults_bit_identical(error)
      !! (d): a namelist that never mentions the EOS knobs leaves all
      !! five members at the pre-capability `eos_t` defaults.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      character(len=:), allocatable :: nml
      integer :: ierr

      nml = "&sim_nml sim_type = 'ocean' /"//new_line("a")// &
            "&grid_nml nx = 6, ny = 5, nghost = 2, dx = 1000.0, dy = 1000.0 /"//new_line("a")

      call build_state(cfg, grid, state, nml, ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "bare namelist must parse")
      if (allocated(error)) return

      checks: block
         call check(error, state%eos%alpha_T == LEGACY_EOS_ALPHA_T, &
                    "default alpha_T must stay 1.7e-4")
         if (allocated(error)) exit checks
         call check(error, state%eos%beta_S == LEGACY_EOS_BETA_S, &
                    "default beta_S must stay 7.6e-4")
         if (allocated(error)) exit checks
         call check(error, state%eos%T_ref == LEGACY_EOS_T_REF, &
                    "default T_ref must stay 10.0")
         if (allocated(error)) exit checks
         call check(error, state%eos%S_ref == LEGACY_EOS_S_REF, &
                    "default S_ref must stay 35.0")
         if (allocated(error)) exit checks
         call check(error, state%eos%rho0 == LEGACY_EOS_RHO0, &
                    "default rho_0 must stay 1035.0")
      end block checks

      call state%destroy()
   end subroutine test_defaults_bit_identical

   subroutine assert_default_config_valid(error)
      !! Baseline for the retirement checks below: an untouched
      !! `config_t` must validate cleanly, so a failure in those tests
      !! can only come from the one key they move.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      integer :: ierr

      ierr = -999
      call validate_config(cfg, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK, &
                 "default config must validate cleanly (retirement baseline)")
   end subroutine assert_default_config_valid

   subroutine test_tracer_alpha_t_retired(error)
      !! (c): `&tracer_nml alpha_T` only ever reached
      !! `tracer_t%eos_coeff`; the live spelling is `&ocean_ic_nml
      !! alpha_T`.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      integer :: ierr

      call assert_default_config_valid(error)
      if (allocated(error)) return

      cfg%alpha_T = 0.0_wp
      ierr = -999
      call validate_config(cfg, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_ERR_CONFIG_VALIDATE, &
                 "&tracer_nml alpha_T off its default must fail loud")
   end subroutine test_tracer_alpha_t_retired

   subroutine test_tracer_beta_s_retired(error)
      !! (c): live spelling is `&ocean_ic_nml beta_S`.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      integer :: ierr

      call assert_default_config_valid(error)
      if (allocated(error)) return

      cfg%beta_S = 0.0_wp
      ierr = -999
      call validate_config(cfg, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_ERR_CONFIG_VALIDATE, &
                 "&tracer_nml beta_S off its default must fail loud")
   end subroutine test_tracer_beta_s_retired

   subroutine test_tracer_t_ref_retired(error)
      !! (c): live spelling is `&ocean_ic_nml T_ref`.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      integer :: ierr

      call assert_default_config_valid(error)
      if (allocated(error)) return

      cfg%T_ref = ISOMIP_T_REF
      ierr = -999
      call validate_config(cfg, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_ERR_CONFIG_VALIDATE, &
                 "&tracer_nml T_ref off its default must fail loud")
   end subroutine test_tracer_t_ref_retired

   subroutine test_tracer_s_ref_retired(error)
      !! (c): live spelling is `&ocean_ic_nml S_ref`.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      integer :: ierr

      call assert_default_config_valid(error)
      if (allocated(error)) return

      cfg%S_ref = ISOMIP_S_REF
      ierr = -999
      call validate_config(cfg, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_ERR_CONFIG_VALIDATE, &
                 "&tracer_nml S_ref off its default must fail loud")
   end subroutine test_tracer_s_ref_retired

end module test_ocean_linear_eos_knobs
