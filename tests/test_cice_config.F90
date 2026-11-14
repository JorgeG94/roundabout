!! Sea-ice (CICE port) configuration: analytical test A20.
!!
!! Every refusal fires with its message, the defaults are the upstream code
!! defaults (transcribed from the upstream parameter tables, not computed by
!! the module under test), both constants sets select correctly and
!! `ice_params_from_config` round-trips each knob.  This test links the main
!! library because the config types live there; it is still a `cice` test.
module test_cice_config
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_config, only: config_t, read_config_from_string, validate_config
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_CONFIG_VALIDATE
   use rdb_ice_dims, only: ICE_NCAT_MAX, ICE_NILYR_MAX
   use rdb_ice_params, only: ICE_ADVECT_UPWIND, ICE_ADVECT_REMAP, ice_params_t, ICE_CONSTANTS_ICEPACK, ICE_CONSTANTS_CESM, &
                             ICE_TFRZ_LINEAR_SALT, ICE_TFRZ_CONSTANT, ICE_TFRZ_MINUS1P8, &
                             ICE_CONDUCT_BUBBLY, ICE_CONDUCT_MU71, ICE_VISC_AVG_ZETA, ICE_VISC_AVG_STRENGTH
   use rdb_ice_shortwave, only: ice_albedo_fhtan, ice_albedos_ccsm3, ice_absorbed_ccsm3, &
                                ICE_ALBICEV_DEFAULT, ICE_ALBICEI_DEFAULT, ICE_ALBSNOWV_DEFAULT, &
                                ICE_ALBSNOWI_DEFAULT, ICE_AHMAX_DEFAULT, ICE_SNOWPATCH_DEFAULT
   use rdb_ice_config, only: ice_params_from_config, ice_config_check, &
                             ICE_CONFIG_MSG_LEN, ICE_CONFIG_MAX_MSGS
   implicit none
   private

   public :: collect_cice_config_tests

   character(len=*), parameter :: BASE_NML = &
                                  "&sim_nml sim_type = 'ocean' /"//new_line('a')// &
                                  "&grid_nml nx = 24, ny = 16, nghost = 3, dx = 20000.0, dy = 20000.0 /"//new_line('a')

contains

   subroutine collect_cice_config_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("defaults_are_upstream", test_defaults), &
                  new_unittest("default_sis2_validates", test_default_validates), &
                  new_unittest("scheme_cice_refused_not_runnable", test_scheme_cice), &
                  new_unittest("structure_refusals", test_structure_refusals), &
                  new_unittest("structure_limits", test_structure_limits), &
                  new_unittest("unsupported_option_refusals", test_option_refusals), &
                  new_unittest("ahmax_nonpositive_refused", test_ahmax), &
                  new_unittest("all_reasons_reported", test_all_reasons), &
                  new_unittest("constants_sets_select", test_constants_sets), &
                  new_unittest("params_round_trip", test_round_trip), &
                  new_unittest("ridge_atm_floe_defaults_are_upstream", test_ridge_atm_defaults), &
                  new_unittest("ridge_refusals", test_ridge_refusals), &
                  new_unittest("atm_floe_refusals", test_atm_floe_refusals), &
                  new_unittest("atm_cesm_constants", test_atm_cesm_constants), &
                  new_unittest("ridge_atm_all_reasons_reported", test_ridge_atm_all_reasons), &
                  new_unittest("ridge_atm_round_trip_om3", test_ridge_atm_round_trip), &
                  new_unittest("sw_defaults_are_upstream", test_sw_defaults), &
                  new_unittest("sw_refusals", test_sw_refusals), &
                  new_unittest("sw_edge_values_accepted", test_sw_edges), &
                  new_unittest("sw_all_reasons_reported", test_sw_all_reasons), &
                  new_unittest("sw_round_trip_om3", test_sw_round_trip), &
                  new_unittest("sw_kernel_default_params_bitwise", test_sw_kernel_defaults), &
                  new_unittest("transport_defaults_are_upstream", test_transport_defaults), &
                  new_unittest("transport_refusals", test_transport_refusals), &
                  new_unittest("transport_round_trip", test_transport_round_trip), &
                  new_unittest("dyn_defaults_are_upstream", test_dyn_defaults), &
                  new_unittest("dyn_refusals", test_dyn_refusals), &
                  new_unittest("dyn_all_reasons_reported", test_dyn_all_reasons), &
                  new_unittest("dyn_round_trip", test_dyn_round_trip), &
                  new_unittest("step_envelope_defaults", test_step_defaults), &
                  new_unittest("step_envelope_refusals", test_step_refusals), &
                  new_unittest("step_envelope_all_reasons", test_step_all_reasons), &
                  new_unittest("unrepresentable_options_rejected", test_unrepresentable) &
                  ]
   end subroutine collect_cice_config_tests

   subroutine parse(extra, cfg, error)
      character(len=*), intent(in) :: extra
      type(config_t), intent(out) :: cfg
      type(error_type), allocatable, intent(out) :: error
      integer :: ierr

      ierr = -999
      call read_config_from_string(BASE_NML//extra, cfg, ierr=ierr)
      call check(error, ierr, OCEAN_STATUS_OK, "namelist text must parse: "//extra)
   end subroutine parse

   function parse_fails(extra) result(failed)
      !! The namelist text is rejected by the parser.
      character(len=*), intent(in) :: extra
      logical :: failed
      type(config_t) :: cfg
      integer :: ierr

      ierr = OCEAN_STATUS_OK
      call read_config_from_string(BASE_NML//extra, cfg, ierr=ierr)
      failed = ierr /= OCEAN_STATUS_OK
   end function parse_fails

   function count_msgs(cfg) result(n)
      type(config_t), intent(in) :: cfg
      integer :: n
      character(len=ICE_CONFIG_MSG_LEN) :: msgs(ICE_CONFIG_MAX_MSGS)
      call ice_config_check(cfg, msgs, n)
   end function count_msgs

   function has_msg(cfg, needle) result(found)
      type(config_t), intent(in) :: cfg
      character(len=*), intent(in) :: needle
      logical :: found
      character(len=ICE_CONFIG_MSG_LEN) :: msgs(ICE_CONFIG_MAX_MSGS)
      integer :: n, i

      call ice_config_check(cfg, msgs, n)
      found = .false.
      do i = 1, n
         if (index(msgs(i), needle) > 0) found = .true.
      end do
   end function has_msg

   subroutine expect_refusal(extra, needle, error)
      !! `extra` parses, and the check reports a message containing `needle`.
      character(len=*), intent(in) :: extra, needle
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      character(len=ICE_CONFIG_MSG_LEN) :: msgs(ICE_CONFIG_MAX_MSGS)
      integer :: n, i
      logical :: found

      call parse(extra, cfg, error)
      if (allocated(error)) return
      call ice_config_check(cfg, msgs, n)
      found = .false.
      do i = 1, n
         if (index(msgs(i), needle) > 0) found = .true.
      end do
      call check(error, found, "refusal containing '"//needle//"' for: "//extra)
   end subroutine expect_refusal

   subroutine test_defaults(error)
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(ice_params_t) :: p
      real(wp), parameter :: tol = 1.0e-15_wp

      cfg = config_t()
      call check(error, count_msgs(cfg), 0, "default config has no refusal")
      if (allocated(error)) return
      call check(error, trim(cfg%ocean%ice%scheme), "sis2")
      if (allocated(error)) return
      p = ice_params_from_config(cfg)
      ! Upstream CICE code defaults (not OM3 values).
      call check(error, p%constants_set, ICE_CONSTANTS_ICEPACK)
      if (allocated(error)) return
      call check(error, p%ncat, 5)
      if (allocated(error)) return
      call check(error, p%nilyr, 7)
      if (allocated(error)) return
      call check(error, p%nslyr, 1)
      if (allocated(error)) return
      call check(error, p%itd%kcatbound, 1)
      if (allocated(error)) return
      call check(error, p%itd%kitd, 1)
      if (allocated(error)) return
      call check(error, p%thermo%ktherm, 1)
      if (allocated(error)) return
      call check(error, p%thermo%conduct, ICE_CONDUCT_BUBBLY)
      if (allocated(error)) return
      call check(error, p%thermo%tfrz_option, ICE_TFRZ_LINEAR_SALT)
      if (allocated(error)) return
      call check(error, p%thermo%congel_one_step, .false.)
      if (allocated(error)) return
      call check(error, p%thermo%update_ocn_f, .false.)
      if (allocated(error)) return
      call check(error, p%thermo%calc_Tsfc, .true.)
      if (allocated(error)) return
      call check(error, abs(p%thermo%ice_ref_salinity - 4.0_wp) < tol, "ice_ref_salinity 4")
      if (allocated(error)) return
      call check(error, abs(p%thermo%ustar_min - 0.005_wp) < tol, "ustar_min 0.005")
      if (allocated(error)) return
      call check(error, abs(p%thermo%emissivity - 0.985_wp) < tol, "emissivity 0.985")
      if (allocated(error)) return
      call check(error, abs(p%thermo%hfrazilmin - 0.05_wp) < tol, "hfrazilmin 0.05")
      if (allocated(error)) return
      call check(error, abs(p%thermo%hi_min - 0.01_wp) < tol, "hi_min 0.01")
      if (allocated(error)) return
      call check(error, abs(p%thermo%floediam - 300.0_wp) < tol, "floediam 300")
      if (allocated(error)) return
      call check(error, abs(p%c%ksno - 0.3_wp) < tol, "ksno 0.3")
      if (allocated(error)) return
      call check(error, abs(p%c%snowpatch - 0.02_wp) < tol, "snowpatch follows icepack set")
      if (allocated(error)) return
      call check(error, abs(p%sw%albicev - 0.78_wp) < tol, "albicev 0.78")
      if (allocated(error)) return
      call check(error, abs(p%sw%albicei - 0.36_wp) < tol, "albicei 0.36")
      if (allocated(error)) return
      call check(error, abs(p%sw%albsnowv - 0.98_wp) < tol, "albsnowv 0.98")
      if (allocated(error)) return
      call check(error, abs(p%sw%albsnowi - 0.70_wp) < tol, "albsnowi 0.70")
      if (allocated(error)) return
      call check(error, abs(p%sw%ahmax - 0.3_wp) < tol, "ahmax 0.3")
   end subroutine test_defaults

   subroutine test_default_validates(error)
      !! The default scheme and the new knobs at their defaults leave
      !! `validate_config` accepting a minimal ocean namelist.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      integer :: ierr

      call parse("", cfg, error)
      if (allocated(error)) return
      ierr = -999
      call validate_config(cfg, ierr=ierr)
      call check(error, ierr, OCEAN_STATUS_OK, "default sis2 namelist validates")
   end subroutine test_default_validates

   subroutine test_scheme_cice(error)
      !! The parser accepts scheme='cice'; configure refuses it, and the
      !! full `validate_config` path (what --validate-only runs) reports it.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      integer :: ierr

      call expect_refusal("&ocean_ice_nml enable = .true., scheme = 'cice' /", &
                          "scheme='cice' needs enable_components", error)
      if (allocated(error)) return
      call parse("&ocean_ice_nml enable = .true., scheme = 'cice' /", cfg, error)
      if (allocated(error)) return
      ierr = -999
      call validate_config(cfg, ierr=ierr)
      call check(error, ierr, OCEAN_STATUS_ERR_CONFIG_VALIDATE, &
                 "validate_config refuses a scheme='cice' namelist outside the envelope")
      if (allocated(error)) return
      call check(error, parse_fails("&ocean_ice_nml scheme = 'bogus' /"), &
                 "an unknown scheme is a parse error")
   end subroutine test_scheme_cice

   subroutine test_structure_refusals(error)
      type(error_type), allocatable, intent(out) :: error

      call expect_refusal("&ocean_ice_nml scheme='cice', nslyr=2 /", "nslyr /= 1", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_nml scheme='cice', nilyr=1 /", "nilyr >= 2", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_nml scheme='cice', ncat=1 /", "ncat >= 2", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_nml scheme='cice', kcatbound=-1 /", "kcatbound must be 0..3", &
                          error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_nml scheme='cice', kcatbound=4 /", "kcatbound must be 0..3", &
                          error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_nml scheme='cice', kcatbound=2, ncat=4 /", &
                          "kcatbound=2", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_nml scheme='cice', kcatbound=2, ncat=8 /", &
                          "kcatbound=2", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_nml scheme='cice', hlim=0.0,0.6,1.4 /", "hlim", error)
   end subroutine test_structure_refusals

   subroutine test_structure_limits(error)
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      character(len=16) :: s1, s2

      write (s1, "(i0)") ICE_NILYR_MAX + 1
      write (s2, "(i0)") ICE_NCAT_MAX + 1
      call expect_refusal("&ocean_ice_nml scheme='cice', nilyr="//trim(s1)//" /", &
                          "ICE_NILYR_MAX", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_nml scheme='cice', ncat="//trim(s2)//" /", &
                          "ICE_NCAT_MAX", error)
      if (allocated(error)) return
      ! The maxima themselves are legal: only the 'not runnable' reason remains.
      write (s1, "(i0)") ICE_NILYR_MAX
      write (s2, "(i0)") ICE_NCAT_MAX
      call parse("&ocean_ice_nml scheme='cice', nilyr="//trim(s1)//", ncat="//trim(s2)//" /"// &
                 new_line('a')//"&ocean_ice_dyn_nml advection='upwind' /", cfg, error)
      if (allocated(error)) return
 call check(error, count_msgs(cfg), 5, "only the envelope reasons (enable, components, hfreeze, dynamics, transport) at the maxima")
      if (allocated(error)) return
      ! Under scheme='sis2' the structure keys are inert.
      call parse("&ocean_ice_nml nilyr=1, nslyr=3, kcatbound=9 /", cfg, error)
      if (allocated(error)) return
      call check(error, count_msgs(cfg), 0, "structure keys inert under sis2")
   end subroutine test_structure_limits

   subroutine test_option_refusals(error)
      type(error_type), allocatable, intent(out) :: error

      call expect_refusal("&ocean_ice_thermo_nml ktherm=2 /", "ktherm /= 1", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_thermo_nml ktherm=0 /", "ktherm /= 1", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_thermo_nml tfrz_option='mushy' /", "mushy", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_thermo_nml saltflux_option='prognostic' /", &
                          "saltflux_option", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_thermo_nml fbot_xfer_type='Cdn_ocn' /", &
                          "fbot_xfer_type", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_thermo_nml cpl_frazil='internal' /", "cpl_frazil", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_thermo_nml cpl_frazil='external' /", "cpl_frazil", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_thermo_nml floediam=0.0 /", "floediam must be > 0", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_thermo_nml kitd=0 /", "kitd /= 1", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_thermo_nml calc_Tsfc=.false. /", "calc_Tsfc", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_thermo_nml semi_implicit_Tsfc=.true. /", &
                          "semi_implicit_Tsfc", error)
   end subroutine test_option_refusals

   subroutine test_ahmax(error)
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg

      call expect_refusal("&ocean_ice_sw_nml ahmax=0.0 /", "ahmax must be > 0", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_sw_nml ahmax=-0.1 /", "ahmax must be > 0", error)
      if (allocated(error)) return
      call parse("&ocean_ice_sw_nml ahmax=1.0e-6 /", cfg, error)
      if (allocated(error)) return
      call check(error, count_msgs(cfg), 0, "a tiny positive ahmax is accepted")
   end subroutine test_ahmax

   subroutine test_all_reasons(error)
      !! Every reason is reported, not just the first.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg

      call parse("&ocean_ice_nml scheme='cice', nilyr=1, nslyr=2 /"//new_line('a')// &
                 "&ocean_ice_thermo_nml ktherm=2, kitd=0 /"//new_line('a')// &
                 "&ocean_ice_sw_nml ahmax=0.0 /"//new_line('a')// &
                 "&ocean_ice_dyn_nml advection='upwind' /", cfg, error)
      if (allocated(error)) return
      call check(error, count_msgs(cfg), 10, "envelope x5 + nilyr + nslyr + ktherm + kitd + ahmax")
      if (allocated(error)) return
      call check(error, has_msg(cfg, "needs enable=.true."), "envelope reason")
      if (allocated(error)) return
      call check(error, has_msg(cfg, "nilyr >= 2"), "nilyr reason")
      if (allocated(error)) return
      call check(error, has_msg(cfg, "nslyr /= 1"), "nslyr reason")
      if (allocated(error)) return
      call check(error, has_msg(cfg, "ktherm /= 1"), "ktherm reason")
      if (allocated(error)) return
      call check(error, has_msg(cfg, "kitd /= 1"), "kitd reason")
      if (allocated(error)) return
      call check(error, has_msg(cfg, "ahmax must be > 0"), "ahmax reason")
   end subroutine test_all_reasons

   subroutine test_constants_sets(error)
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(ice_params_t) :: p
      real(wp), parameter :: tol = 1.0e-12_wp

      call parse("&ocean_ice_nml constants_set='cesm' /", cfg, error)
      if (allocated(error)) return
      p = ice_params_from_config(cfg)
      call check(error, p%constants_set, ICE_CONSTANTS_CESM)
      if (allocated(error)) return
      call check(error, abs(p%c%cp_ocn - 3996.0_wp) < tol, "cesm cp_ocn")
      if (allocated(error)) return
      call check(error, abs(p%c%snowpatch - 0.005_wp) < tol, "cesm snowpatch follows the set")
      if (allocated(error)) return
      ! An explicit snowpatch wins over either set.
      call parse("&ocean_ice_nml constants_set='cesm' /"//new_line('a')// &
                 "&ocean_ice_sw_nml snowpatch=0.013 /", cfg, error)
      if (allocated(error)) return
      p = ice_params_from_config(cfg)
      call check(error, abs(p%c%snowpatch - 0.013_wp) < tol, "explicit snowpatch (cesm)")
      if (allocated(error)) return
      call parse("&ocean_ice_nml constants_set='icepack' /"//new_line('a')// &
                 "&ocean_ice_sw_nml snowpatch=0.013 /", cfg, error)
      if (allocated(error)) return
      p = ice_params_from_config(cfg)
      call check(error, p%constants_set, ICE_CONSTANTS_ICEPACK)
      if (allocated(error)) return
      call check(error, abs(p%c%cp_ocn - 4218.0_wp) < tol, "icepack cp_ocn")
      if (allocated(error)) return
      call check(error, abs(p%c%snowpatch - 0.013_wp) < tol, "explicit snowpatch (icepack)")
      if (allocated(error)) return
      call check(error, parse_fails("&ocean_ice_nml constants_set='bogus' /"), &
                 "an unknown constants_set is a parse error")
   end subroutine test_constants_sets

   subroutine test_round_trip(error)
      !! Every non-default knob reaches `ice_params_t` unchanged (OM3 values).
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(ice_params_t) :: p
      real(wp), parameter :: tol = 1.0e-15_wp

      call parse("&ocean_ice_nml ncat=5, nilyr=4, nslyr=1, kcatbound=0, "// &
                 "constants_set='cesm' /"//new_line('a')// &
                 "&ocean_ice_thermo_nml conduct='MU71', tfrz_option='constant', "// &
                 "ice_ref_salinity=5.0, congel_freeze='one-step', update_ocn_f=.true., "// &
                 "hfrazilmin=0.07, hi_min=0.02, ustar_min=5.0e-4, floediam=250.0, "// &
                 "ksno=0.31, emissivity=0.95 /"//new_line('a')// &
                 "&ocean_ice_sw_nml albicev=0.86, albicei=0.44, albsnowv=0.97, "// &
                 "albsnowi=0.69, ahmax=0.1, snowpatch=0.011 /", cfg, error)
      if (allocated(error)) return
      call check(error, count_msgs(cfg), 0, "OM3-style values need no refusal under sis2")
      if (allocated(error)) return
      p = ice_params_from_config(cfg)
      call check(error, p%ncat, 5)
      if (allocated(error)) return
      call check(error, p%nilyr, 4)
      if (allocated(error)) return
      call check(error, p%nslyr, 1)
      if (allocated(error)) return
      call check(error, p%itd%kcatbound, 0)
      if (allocated(error)) return
      call check(error, p%constants_set, ICE_CONSTANTS_CESM)
      if (allocated(error)) return
      call check(error, p%thermo%conduct, ICE_CONDUCT_MU71)
      if (allocated(error)) return
      call check(error, p%thermo%tfrz_option, ICE_TFRZ_CONSTANT)
      if (allocated(error)) return
      call check(error, p%thermo%congel_one_step, .true.)
      if (allocated(error)) return
      call check(error, p%thermo%update_ocn_f, .true.)
      if (allocated(error)) return
      call check(error, abs(p%thermo%ice_ref_salinity - 5.0_wp) < tol, "ice_ref_salinity")
      if (allocated(error)) return
      call check(error, abs(p%thermo%hfrazilmin - 0.07_wp) < tol, "hfrazilmin")
      if (allocated(error)) return
      call check(error, abs(p%thermo%hi_min - 0.02_wp) < tol, "hi_min")
      if (allocated(error)) return
      call check(error, abs(p%thermo%ustar_min - 5.0e-4_wp) < tol, "ustar_min")
      if (allocated(error)) return
      call check(error, abs(p%thermo%floediam - 250.0_wp) < tol, "floediam")
      if (allocated(error)) return
      call check(error, abs(p%c%ksno - 0.31_wp) < tol, "ksno")
      if (allocated(error)) return
      call check(error, abs(p%thermo%emissivity - 0.95_wp) < tol, "emissivity")
      if (allocated(error)) return
      call check(error, abs(p%sw%albicev - 0.86_wp) < tol, "albicev")
      if (allocated(error)) return
      call check(error, abs(p%sw%albicei - 0.44_wp) < tol, "albicei")
      if (allocated(error)) return
      call check(error, abs(p%sw%albsnowv - 0.97_wp) < tol, "albsnowv")
      if (allocated(error)) return
      call check(error, abs(p%sw%albsnowi - 0.69_wp) < tol, "albsnowi")
      if (allocated(error)) return
      call check(error, abs(p%sw%ahmax - 0.1_wp) < tol, "ahmax")
      if (allocated(error)) return
      call check(error, abs(p%c%snowpatch - 0.011_wp) < tol, "snowpatch")
      if (allocated(error)) return
      ! The other two tfrz spellings.
      call parse("&ocean_ice_thermo_nml tfrz_option='minus1p8' /", cfg, error)
      if (allocated(error)) return
      p = ice_params_from_config(cfg)
      call check(error, p%thermo%tfrz_option, ICE_TFRZ_MINUS1P8)
   end subroutine test_round_trip

   subroutine test_ridge_atm_defaults(error)
      !! Upstream CICE code defaults, transcribed from its parameter table.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(ice_params_t) :: p
      real(wp), parameter :: tol = 1.0e-15_wp

      cfg = config_t()
      p = ice_params_from_config(cfg)
      call check(error, p%ridge%kstrength, 1)
      if (allocated(error)) return
      call check(error, p%ridge%krdg_partic, 1)
      if (allocated(error)) return
      call check(error, p%ridge%krdg_redist, 1)
      if (allocated(error)) return
      call check(error, abs(p%ridge%mu_rdg - 3.0_wp) < tol, "mu_rdg 3")
      if (allocated(error)) return
      call check(error, abs(p%ridge%Cf - 17.0_wp) < tol, "Cf 17")
      if (allocated(error)) return
      call check(error, abs(p%ridge%Pstar - 2.75e4_wp) < tol*1.0e4_wp, "Pstar 2.75e4")
      if (allocated(error)) return
      call check(error, abs(p%ridge%Cstar - 20.0_wp) < tol, "Cstar 20")
      if (allocated(error)) return
      call check(error, p%atm%highfreq, .false.)
      if (allocated(error)) return
      call check(error, p%atm%natmiter, 5)
      if (allocated(error)) return
      call check(error, abs(p%atm%atmiter_conv) < tol, "atmiter_conv 0")
      if (allocated(error)) return
      call check(error, abs(p%atm%iceruf - 5.0e-4_wp) < tol, "iceruf 5e-4")
      if (allocated(error)) return
      ! rdn = kappa/ln(zref/iceruf) = 0.4/ln(20000) written out from the closed form.
      call check(error, abs(p%atm%rdn - 0.4_wp/log(10.0_wp/5.0e-4_wp)) < 1.0e-12_wp, &
                 "neutral drag finalised at the default roughness")
      if (allocated(error)) return
      call check(error, abs(p%thermo%floeshape - 0.66_wp) < tol, "floeshape 0.66")
      if (allocated(error)) return
      call check(error, abs(p%thermo%floediam - 300.0_wp) < tol, "floediam 300")
      if (allocated(error)) return
      call check(error, count_msgs(cfg), 0, "defaults need no refusal")
   end subroutine test_ridge_atm_defaults

   subroutine test_ridge_refusals(error)
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg

      call expect_refusal("&ocean_ice_ridge_nml kstrength=2 /", "kstrength must be 0", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_ridge_nml kstrength=-1 /", "kstrength must be 0", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_ridge_nml krdg_partic=2 /", "krdg_partic must be 0", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_ridge_nml krdg_redist=3 /", "krdg_redist must be 0", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_ridge_nml mu_rdg=0.0 /", "mu_rdg must be > 0", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_ridge_nml mu_rdg=-1.0 /", "mu_rdg must be > 0", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_ridge_nml Cf=0.0 /", "Cf must be > 0", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_ridge_nml Cf=-17.0 /", "Cf must be > 0", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_ridge_nml Pstar=-2.75e4 /", "Pstar must be > 0", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_ridge_nml Pstar=0.0 /", "Pstar must be > 0", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_ridge_nml Cstar=-20.0 /", "Cstar must be > 0", error)
      if (allocated(error)) return
      ! Every legal option combination is accepted (upstream documents no
      ! discouraged pairing, so there is no warning either).
      call parse("&ocean_ice_ridge_nml kstrength=1, krdg_partic=0, krdg_redist=0 /", cfg, error)
      if (allocated(error)) return
      call check(error, count_msgs(cfg), 0, "legal option combination accepted")
   end subroutine test_ridge_refusals

   subroutine test_atm_floe_refusals(error)
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg

      call expect_refusal("&ocean_ice_atm_nml natmiter=-1 /", "natmiter must be >= 0", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_atm_nml atmiter_conv=-1.0e-3 /", "atmiter_conv must be >= 0", &
                          error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_atm_nml iceruf=0.0 /", "iceruf must be > 0 and below", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_atm_nml iceruf=-1.0e-3 /", "iceruf must be > 0 and below", &
                          error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_atm_nml iceruf=10.0 /", "iceruf must be > 0 and below", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_atm_nml iceruf=12.5 /", "iceruf must be > 0 and below", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_atm_nml iceruf=12.5 /", "zref=10", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_atm_nml atmbndy='constant' /", "atmbndy /= 'similarity'", &
                          error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_atm_nml atmbndy='mixed' /", "atmbndy /= 'similarity'", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_atm_nml formdrag=.true. /", "formdrag=.true.", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_atm_nml calc_strair=.false. /", "calc_strair=.false.", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_thermo_nml floeshape=0.0 /", "floeshape must be > 0", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_thermo_nml floeshape=-0.66 /", "floeshape must be > 0", error)
      if (allocated(error)) return
      call check(error, parse_fails("&ocean_ice_thermo_nml floediam=-5.0 /"), "negative floediam is a parse error (schema minimum)")
      if (allocated(error)) return
      ! natmiter = 0 (neutral drag only) and a roughness just below zref are legal.
      call parse("&ocean_ice_atm_nml natmiter=0, iceruf=9.99 /", cfg, error)
      if (allocated(error)) return
      call check(error, count_msgs(cfg), 0, "natmiter=0 and iceruf<zref accepted")
   end subroutine test_atm_floe_refusals

   subroutine test_atm_cesm_constants(error)
      !! The iceruf bound and the neutral drag follow the selected constants set
      !! (zref is 10 m in both sets; rdn is the closed form for the set's kappa).
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(ice_params_t) :: p

      call parse("&ocean_ice_nml constants_set='cesm' /"//new_line('a')// &
                 "&ocean_ice_atm_nml iceruf=2.0e-3 /", cfg, error)
      if (allocated(error)) return
      call check(error, count_msgs(cfg), 0, "cesm set with iceruf=2e-3 accepted")
      if (allocated(error)) return
      p = ice_params_from_config(cfg)
      call check(error, p%constants_set, ICE_CONSTANTS_CESM)
      if (allocated(error)) return
      call check(error, abs(p%atm%rdn - p%c%vonkar/log(p%c%zref/2.0e-3_wp)) < 1.0e-12_wp, &
                 "rdn = kappa/ln(zref/iceruf) with the cesm set's constants")
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_nml constants_set='cesm' /"//new_line('a')// &
                          "&ocean_ice_atm_nml iceruf=11.0 /", "iceruf must be > 0 and below", error)
   end subroutine test_atm_cesm_constants

   subroutine test_ridge_atm_all_reasons(error)
      !! Every new refusal is reported together, not just the first.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg

      call parse("&ocean_ice_ridge_nml kstrength=2, krdg_partic=5, krdg_redist=-1, "// &
                 "mu_rdg=0.0, Cf=0.0, Pstar=0.0, Cstar=0.0 /"//new_line('a')// &
                 "&ocean_ice_atm_nml natmiter=-2, atmiter_conv=-1.0, iceruf=20.0, "// &
                 "atmbndy='mixed', formdrag=.true., calc_strair=.false. /"//new_line('a')// &
                 "&ocean_ice_thermo_nml floediam=0.0, floeshape=0.0 /", cfg, error)
      if (allocated(error)) return
      call check(error, count_msgs(cfg), 15, "7 ridge + 6 atm + 2 floe reasons")
      if (allocated(error)) return
      call check(error, has_msg(cfg, "kstrength must be 0"), "kstrength reason")
      if (allocated(error)) return
      call check(error, has_msg(cfg, "krdg_partic must be 0"), "krdg_partic reason")
      if (allocated(error)) return
      call check(error, has_msg(cfg, "krdg_redist must be 0"), "krdg_redist reason")
      if (allocated(error)) return
      call check(error, has_msg(cfg, "mu_rdg must be > 0"), "mu_rdg reason")
      if (allocated(error)) return
      call check(error, has_msg(cfg, "Cf must be > 0"), "Cf reason")
      if (allocated(error)) return
      call check(error, has_msg(cfg, "Pstar must be > 0"), "Pstar reason")
      if (allocated(error)) return
      call check(error, has_msg(cfg, "Cstar must be > 0"), "Cstar reason")
      if (allocated(error)) return
      call check(error, has_msg(cfg, "natmiter must be >= 0"), "natmiter reason")
      if (allocated(error)) return
      call check(error, has_msg(cfg, "atmiter_conv must be >= 0"), "atmiter_conv reason")
      if (allocated(error)) return
      call check(error, has_msg(cfg, "iceruf must be > 0 and below"), "iceruf reason")
      if (allocated(error)) return
      call check(error, has_msg(cfg, "atmbndy /= 'similarity'"), "atmbndy reason")
      if (allocated(error)) return
      call check(error, has_msg(cfg, "formdrag=.true."), "formdrag reason")
      if (allocated(error)) return
      call check(error, has_msg(cfg, "calc_strair=.false."), "calc_strair reason")
      if (allocated(error)) return
      call check(error, has_msg(cfg, "floediam must be > 0"), "floediam reason")
      if (allocated(error)) return
      call check(error, has_msg(cfg, "floeshape must be > 0"), "floeshape reason")
   end subroutine test_ridge_atm_all_reasons

   subroutine test_ridge_atm_round_trip(error)
      !! OM3 sets highfreq=.true. and no ridge key; the rest are non-default
      !! values chosen to be distinguishable.  `atm%rdn` follows `iceruf`.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(ice_params_t) :: p
      real(wp), parameter :: tol = 1.0e-15_wp

      call parse("&ocean_ice_atm_nml highfreq=.true. /", cfg, error)
      if (allocated(error)) return
      call check(error, count_msgs(cfg), 0, "OM3 atmosphere settings accepted")
      if (allocated(error)) return
      p = ice_params_from_config(cfg)
      call check(error, p%atm%highfreq, .true.)
      if (allocated(error)) return

      call parse("&ocean_ice_ridge_nml kstrength=0, krdg_partic=0, krdg_redist=0, "// &
                 "mu_rdg=4.5, Cf=12.0, Pstar=3.0e4, Cstar=15.0 /"//new_line('a')// &
                 "&ocean_ice_atm_nml highfreq=.true., natmiter=3, atmiter_conv=1.0e-3, "// &
                 "iceruf=1.0e-3, atmbndy='similarity', formdrag=.false., "// &
                 "calc_strair=.true. /"//new_line('a')// &
                 "&ocean_ice_thermo_nml floediam=250.0, floeshape=0.5 /", cfg, error)
      if (allocated(error)) return
      call check(error, count_msgs(cfg), 0, "non-default legal values accepted")
      if (allocated(error)) return
      p = ice_params_from_config(cfg)
      call check(error, p%ridge%kstrength, 0)
      if (allocated(error)) return
      call check(error, p%ridge%krdg_partic, 0)
      if (allocated(error)) return
      call check(error, p%ridge%krdg_redist, 0)
      if (allocated(error)) return
      call check(error, abs(p%ridge%mu_rdg - 4.5_wp) < tol, "mu_rdg")
      if (allocated(error)) return
      call check(error, abs(p%ridge%Cf - 12.0_wp) < tol, "Cf")
      if (allocated(error)) return
      call check(error, abs(p%ridge%Pstar - 3.0e4_wp) < tol*1.0e4_wp, "Pstar")
      if (allocated(error)) return
      call check(error, abs(p%ridge%Cstar - 15.0_wp) < tol, "Cstar")
      if (allocated(error)) return
      call check(error, p%atm%natmiter, 3)
      if (allocated(error)) return
      call check(error, abs(p%atm%atmiter_conv - 1.0e-3_wp) < tol, "atmiter_conv")
      if (allocated(error)) return
      call check(error, abs(p%atm%iceruf - 1.0e-3_wp) < tol, "iceruf")
      if (allocated(error)) return
      ! rdn finalised from the configured roughness: kappa/ln(zref/iceruf).
      call check(error, abs(p%atm%rdn - 0.4_wp/log(10.0_wp/1.0e-3_wp)) < 1.0e-12_wp, "rdn")
      if (allocated(error)) return
      call check(error, abs(p%thermo%floeshape - 0.5_wp) < tol, "floeshape")
      if (allocated(error)) return
      call check(error, abs(p%thermo%floediam - 250.0_wp) < tol, "floediam")
   end subroutine test_ridge_atm_round_trip

   subroutine test_sw_defaults(error)
      !! Upstream code defaults of the shortwave knobs (the values are written here from
      !! the upstream parameter tables, not read back from the module under test).
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(ice_params_t) :: p
      real(wp), parameter :: tol = 1.0e-15_wp

      cfg = config_t()
      call check(error, trim(cfg%ocean%ice_sw%shortwave), "ccsm3", "shortwave default ccsm3")
      if (allocated(error)) return
      call check(error, trim(cfg%ocean%ice_sw%albedo_type), "ccsm3", "albedo_type default ccsm3")
      if (allocated(error)) return
      p = ice_params_from_config(cfg)
      call check(error, p%sw%sw_redist, .false., "sw_redist upstream default .false.")
      if (allocated(error)) return
      call check(error, abs(p%sw%sw_frac - 0.9_wp) < tol, "sw_frac upstream default 0.9")
      if (allocated(error)) return
      call check(error, abs(p%sw%sw_dtemp - 0.02_wp) < tol, "sw_dtemp upstream default 0.02")
      if (allocated(error)) return
      ! Interim roundabout choice, not an upstream value: CICE standalone 0.28 + 0.24.
      call check(error, abs(p%sw%sw_vis_frac - 0.52_wp) < tol, "sw_vis_frac interim 0.52")
   end subroutine test_sw_defaults

   subroutine test_sw_refusals(error)
      type(error_type), allocatable, intent(out) :: error

      call expect_refusal("&ocean_ice_sw_nml shortwave='dEdd' /", "shortwave='dEdd' is not ported", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_sw_nml shortwave='dEdd_snicar_ad' /", "is not ported", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_sw_nml albedo_type='constant' /", &
                          "albedo_type='constant' is not ported", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_sw_nml albicev=1.01 /", "albicev must lie in [0, 1]", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_sw_nml albicei=-0.01 /", "albicei must lie in [0, 1]", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_sw_nml albsnowv=1.5 /", "albsnowv must lie in [0, 1]", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_sw_nml albsnowi=-1.0 /", "albsnowi must lie in [0, 1]", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_sw_nml sw_frac=0.0 /", "sw_frac must lie in (0, 1]", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_sw_nml sw_frac=1.2 /", "sw_frac must lie in (0, 1]", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_sw_nml sw_dtemp=-0.01 /", "sw_dtemp must be >= 0", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_sw_nml sw_vis_frac=-0.1 /", "sw_vis_frac must lie in [0, 1]", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_sw_nml sw_vis_frac=1.1 /", "sw_vis_frac must lie in [0, 1]", error)
   end subroutine test_sw_refusals

   subroutine test_sw_edges(error)
      !! The closed ends of every range are legal.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg

      call parse("&ocean_ice_sw_nml albicev=0.0, albicei=1.0, albsnowv=1.0, albsnowi=0.0, "// &
                 "sw_frac=1.0, sw_dtemp=0.0, sw_vis_frac=0.0 /", cfg, error)
      if (allocated(error)) return
      call check(error, count_msgs(cfg), 0, "closed range ends accepted (vis_frac 0)")
      if (allocated(error)) return
      call parse("&ocean_ice_sw_nml sw_vis_frac=1.0, sw_frac=1.0e-6, shortwave='ccsm3', "// &
                 "albedo_type='ccsm3' /", cfg, error)
      if (allocated(error)) return
      call check(error, count_msgs(cfg), 0, "sw_vis_frac 1 and tiny sw_frac accepted")
   end subroutine test_sw_edges

   subroutine test_sw_all_reasons(error)
      !! Every shortwave reason is reported together, none hidden by the schema.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg

      call parse("&ocean_ice_sw_nml shortwave='dEdd', albedo_type='constant', albicev=2.0, "// &
                 "albicei=-1.0, albsnowv=3.0, albsnowi=-3.0, ahmax=0.0, sw_frac=0.0, "// &
                 "sw_dtemp=-1.0, sw_vis_frac=2.0 /", cfg, error)
      if (allocated(error)) return
      call check(error, count_msgs(cfg), 10, "ten shortwave reasons at once")
      if (allocated(error)) return
      call check(error, has_msg(cfg, "shortwave='dEdd'"), "scheme reason")
      if (allocated(error)) return
      call check(error, has_msg(cfg, "albedo_type='constant'"), "albedo_type reason")
      if (allocated(error)) return
      call check(error, has_msg(cfg, "sw_vis_frac"), "band fraction reason")
   end subroutine test_sw_all_reasons

   subroutine test_sw_round_trip(error)
      !! OM3 values (`ice_in`: albicev .86, albicei .44, albsnowv .98, albsnowi .7, ahmax .1,
      !! sw_redist .true.; sw_frac and sw_dtemp unset so upstream) reach `ice_sw_params_t`.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(ice_params_t) :: p
      real(wp), parameter :: tol = 1.0e-15_wp

      call parse("&ocean_ice_sw_nml shortwave='ccsm3', albedo_type='ccsm3', albicev=0.86, "// &
                 "albicei=0.44, albsnowv=0.98, albsnowi=0.7, ahmax=0.1, sw_redist=.true., "// &
                 "sw_frac=0.8, sw_dtemp=0.03, sw_vis_frac=0.6 /", cfg, error)
      if (allocated(error)) return
      call check(error, count_msgs(cfg), 0, "OM3 shortwave values accepted")
      if (allocated(error)) return
      p = ice_params_from_config(cfg)
      call check(error, abs(p%sw%albicev - 0.86_wp) < tol, "albicev")
      if (allocated(error)) return
      call check(error, abs(p%sw%albicei - 0.44_wp) < tol, "albicei")
      if (allocated(error)) return
      call check(error, abs(p%sw%albsnowv - 0.98_wp) < tol, "albsnowv")
      if (allocated(error)) return
      call check(error, abs(p%sw%albsnowi - 0.7_wp) < tol, "albsnowi")
      if (allocated(error)) return
      call check(error, abs(p%sw%ahmax - 0.1_wp) < tol, "ahmax")
      if (allocated(error)) return
      call check(error, p%sw%sw_redist, .true., "sw_redist")
      if (allocated(error)) return
      call check(error, abs(p%sw%sw_frac - 0.8_wp) < tol, "sw_frac")
      if (allocated(error)) return
      call check(error, abs(p%sw%sw_dtemp - 0.03_wp) < tol, "sw_dtemp")
      if (allocated(error)) return
      call check(error, abs(p%sw%sw_vis_frac - 0.6_wp) < tol, "sw_vis_frac")
   end subroutine test_sw_round_trip

   subroutine test_sw_kernel_defaults(error)
      !! The D7.2 kernels fed with the parameters of the default config give bit-for-bit the
      !! result they give with their own module defaults.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(ice_params_t) :: p
      real(wp), parameter :: hi(4) = [0.05_wp, 0.4_wp, 1.3_wp, 3.0_wp]
      real(wp), parameter :: hs(4) = [0.0_wp, 0.03_wp, 0.2_wp, 0.0_wp]
      real(wp), parameter :: ts(4) = [-15.0_wp, -0.5_wp, 0.0_wp, -3.0_wp]
      real(wp) :: a(8), b(8), fa(3), fb(3), ia(4), ib(4), fhp, fhd
      integer :: m

      cfg = config_t()
      p = ice_params_from_config(cfg)
      fhp = ice_albedo_fhtan(p%sw%ahmax)
      fhd = ice_albedo_fhtan(ICE_AHMAX_DEFAULT)
      do m = 1, 4
         call ice_albedos_ccsm3(hi(m), hs(m), ts(m), p%sw%albicev, p%sw%albicei, p%sw%albsnowv, &
                                p%sw%albsnowi, fhp, p%c%snowpatch, a(1), a(2), a(3), a(4), a(5), a(6), a(7))
         call ice_albedos_ccsm3(hi(m), hs(m), ts(m), ICE_ALBICEV_DEFAULT, ICE_ALBICEI_DEFAULT, &
                                ICE_ALBSNOWV_DEFAULT, ICE_ALBSNOWI_DEFAULT, fhd, ICE_SNOWPATCH_DEFAULT, &
                                b(1), b(2), b(3), b(4), b(5), b(6), b(7))
         call check(error, all(a(1:7) == b(1:7)), "albedos bitwise equal")
         if (allocated(error)) return
         call ice_absorbed_ccsm3(4, hi(m), a(5), 100.0_wp*p%sw%sw_vis_frac, 100.0_wp*(1.0_wp - p%sw%sw_vis_frac), &
                                 a(1), a(2), a(3), a(4), fa(1), fa(2), fa(3), ia)
         call ice_absorbed_ccsm3(4, hi(m), b(5), 100.0_wp*0.52_wp, 100.0_wp*(1.0_wp - 0.52_wp), b(1), b(2), b(3), b(4), &
                                 fb(1), fb(2), fb(3), ib)
         call check(error, all(fa == fb) .and. all(ia == ib), "absorbed solar bitwise equal")
         if (allocated(error)) return
      end do
   end subroutine test_sw_kernel_defaults

   subroutine test_transport_defaults(error)
      !! `advection = 'remap'` (upstream default) and `transport_check = .false.`;
      !! the default configuration (scheme sis2) needs no refusal.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(ice_params_t) :: p

      cfg = config_t()
      p = ice_params_from_config(cfg)
      call check(error, p%transport%advection, ICE_ADVECT_REMAP, "advection default is remap")
      if (allocated(error)) return
      call check(error, p%transport%check, .false., "transport_check default off")
      if (allocated(error)) return
      call check(error, count_msgs(cfg), 0, "defaults need no refusal under sis2")
   end subroutine test_transport_defaults

   subroutine test_transport_refusals(error)
      !! Every refusal of the transport keys, with its message.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg

      ! remap is parsed and accepted under sis2 (inert), refused under cice
      call parse("&ocean_ice_dyn_nml advection='remap' /", cfg, error)
      if (allocated(error)) return
      call check(error, has_msg(cfg, "advection='remap' is not runnable"), .false., "remap inert under sis2")
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_nml scheme='cice' /"//new_line('a')// &
                          "&ocean_ice_dyn_nml advection='remap' /", "advection='remap' is not runnable", error)
      if (allocated(error)) return
      ! upwind is accepted under cice (the scheme itself is still refused: not runnable)
      call parse("&ocean_ice_nml scheme='cice' /"//new_line('a')// &
                 "&ocean_ice_dyn_nml advection='upwind' /", cfg, error)
      if (allocated(error)) return
      call check(error, has_msg(cfg, "advection="), .false., "upwind accepted under cice")
      if (allocated(error)) return
      ! an unknown spelling fails the strict schema
      call check(error, parse_fails("&ocean_ice_dyn_nml advection='ppm' /"), "unknown advection spelling rejected")
      if (allocated(error)) return
      call check(error, parse_fails("&ocean_ice_dyn_nml transport_check = 3 /"), "transport_check is a logical")
   end subroutine test_transport_refusals

   subroutine test_transport_round_trip(error)
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(ice_params_t) :: p

      call parse("&ocean_ice_dyn_nml advection='upwind', transport_check=.true. /", cfg, error)
      if (allocated(error)) return
      p = ice_params_from_config(cfg)
      call check(error, p%transport%advection, ICE_ADVECT_UPWIND, "advection round trip")
      if (allocated(error)) return
      call check(error, p%transport%check, .true., "transport_check round trip")
   end subroutine test_transport_round_trip

   subroutine test_dyn_defaults(error)
      !! `&ocean_ice_dyn_nml` defaults are the upstream CICE code defaults
      !! (written out from its parameter table).
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(ice_params_t) :: p
      real(wp), parameter :: tol = 1.0e-15_wp

      cfg = config_t()
      p = ice_params_from_config(cfg)
      call check(error, p%dyn%ndte, 120)
      if (allocated(error)) return
      call check(error, abs(p%dyn%elasticDamp - 0.36_wp) < tol, "elasticDamp 0.36")
      if (allocated(error)) return
      call check(error, p%dyn%revised_evp, .false.)
      if (allocated(error)) return
      call check(error, abs(p%dyn%arlx - 300.0_wp) < tol .and. abs(p%dyn%brlx - 300.0_wp) < tol, "arlx = brlx = 300")
      if (allocated(error)) return
      call check(error, abs(p%dyn%e_yieldcurve - 2.0_wp) < tol .and. abs(p%dyn%e_plasticpot - 2.0_wp) < tol, "e = 2")
      if (allocated(error)) return
      call check(error, abs(p%dyn%Ktens) < tol, "Ktens 0")
      if (allocated(error)) return
      call check(error, abs(p%dyn%deltaminEVP - 1.0e-11_wp) < 1.0e-26_wp, "deltaminEVP 1e-11")
      if (allocated(error)) return
      call check(error, p%dyn%capping_max, .true., "capping max")
      if (allocated(error)) return
      call check(error, abs(p%dyn%dyn_area_min - 1.0e-11_wp) < 1.0e-26_wp, "dyn_area_min 1e-11")
      if (allocated(error)) return
      call check(error, abs(p%dyn%dyn_mass_min - 1.0e-10_wp) < 1.0e-25_wp, "dyn_mass_min 1e-10")
      if (allocated(error)) return
      call check(error, abs(p%dyn%rheo_area_min - 1.0e-3_wp) < 1.0e-18_wp, "rheo_area_min 1e-3")
      if (allocated(error)) return
      call check(error, p%dyn%free_drift_approx, .true., "free_drift_approx")
      if (allocated(error)) return
      call check(error, p%dyn%visc_method, ICE_VISC_AVG_ZETA)
      if (allocated(error)) return
      call check(error, abs(p%dyn%cdn_ocn - 0.00536_wp) < tol, "cdn_ocn 0.00536")
      if (allocated(error)) return
      call check(error, abs(p%dyn%vel_cap - 10.0_wp) < tol, "vel_cap 10")
      if (allocated(error)) return
      call check(error, abs(p%dyn%cosw - 1.0_wp) < tol .and. abs(p%dyn%sinw) < tol, "no turning angle")
   end subroutine test_dyn_defaults

   subroutine test_dyn_refusals(error)
      !! Every illegal value is accepted by the parser and refused by the check
      !! with a message naming the knob.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg

      call expect_refusal("&ocean_ice_dyn_nml ndte=0 /", "ndte must be >= 1", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_dyn_nml ndte=-5 /", "ndte must be >= 1", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_dyn_nml elasticDamp=0.0 /", "elasticDamp must be > 0", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_dyn_nml revised_evp=.true., arlx=0.0 /", "arlx must be > 0", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_dyn_nml revised_evp=.true., brlx=-1.0 /", "brlx must be > 0", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_dyn_nml e_yieldcurve=0.0 /", "e_yieldcurve must be > 0", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_dyn_nml e_plasticpot=-2.0 /", "e_plasticpot must be > 0", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_dyn_nml Ktens=1.0 /", "Ktens must be in [0, 1)", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_dyn_nml Ktens=-0.1 /", "Ktens must be in [0, 1)", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_dyn_nml deltaminEVP=0.0 /", "deltaminEVP must be > 0", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_dyn_nml capping_method='min' /", "capping_method must be", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_dyn_nml dyn_area_min=-1.0e-3 /", "dyn_area_min must be >= 0", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_dyn_nml dyn_mass_min=-1.0 /", "dyn_mass_min must be >= 0", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_dyn_nml rheo_area_min=-1.0 /", "rheo_area_min must be >= 0", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_dyn_nml visc_method='avg_viscosity' /", "visc_method must be", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_dyn_nml cdn_ocn=-0.1 /", "cdn_ocn must be >= 0", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_dyn_nml vel_cap=0.0 /", "vel_cap must be > 0", error)
      if (allocated(error)) return
      ! legal alternatives are accepted
      call parse("&ocean_ice_dyn_nml revised_evp=.true., capping_method='sum', visc_method='avg_strength', "// &
                 "Ktens=0.2, ndte=240 /", cfg, error)
      if (allocated(error)) return
      call check(error, count_msgs(cfg), 0, "legal alternatives accepted")
      if (allocated(error)) return
      ! revised-only relaxation parameters are not checked under the classic EVP
      call parse("&ocean_ice_dyn_nml arlx=-3.0, brlx=0.0 /", cfg, error)
      if (allocated(error)) return
      call check(error, count_msgs(cfg), 0, "arlx/brlx unused (classic EVP) are not refused")
   end subroutine test_dyn_refusals

   subroutine test_dyn_all_reasons(error)
      !! Every dynamics refusal is reported together (and with the other groups').
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg

      call parse("&ocean_ice_dyn_nml ndte=0, elasticDamp=-1.0, e_yieldcurve=0.0, e_plasticpot=0.0, Ktens=2.0, "// &
                 "deltaminEVP=0.0, capping_method='x', dyn_area_min=-1.0, dyn_mass_min=-1.0, "// &
                 "rheo_area_min=-1.0, visc_method='x', cdn_ocn=-1.0, vel_cap=-1.0 /"//new_line('a')// &
                 "&ocean_ice_ridge_nml Cf=0.0 /", cfg, error)
      if (allocated(error)) return
      call check(error, count_msgs(cfg), 14, "13 dynamics reasons + 1 ridge reason")
      if (allocated(error)) return
      call check(error, has_msg(cfg, "ndte must be >= 1") .and. has_msg(cfg, "elasticDamp must be > 0") .and. &
                 has_msg(cfg, "Ktens must be in") .and. has_msg(cfg, "capping_method must be") .and. &
                 has_msg(cfg, "visc_method must be") .and. has_msg(cfg, "vel_cap must be > 0") .and. &
                 has_msg(cfg, "Cf must be > 0"), "reasons present")
      if (allocated(error)) return
      call parse("&ocean_ice_dyn_nml revised_evp=.true., arlx=0.0, brlx=0.0 /", cfg, error)
      if (allocated(error)) return
      call check(error, count_msgs(cfg), 2, "revised EVP: arlx and brlx both reported")
   end subroutine test_dyn_all_reasons

   subroutine test_step_defaults(error)
      !! The step-envelope knobs default to the supported values, and a plain `scheme='cice'` namelist trips none of
      !! them (only the three envelope reasons a bare namelist lacks).
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg

      cfg = config_t()
      call check(error, cfg%ocean%ice_dyn%ndtd, 1, "ndtd default")
      if (allocated(error)) return
      call check(error, cfg%ocean%ice%dt_ice <= 0.0_wp, "dt_ice default: the outer step")
      if (allocated(error)) return
      call parse("&ocean_ice_nml scheme='cice' /"//new_line('a')//"&ocean_ice_dyn_nml advection='upwind' /", cfg, error)
      if (allocated(error)) return
      call check(error, count_msgs(cfg), 5, "only the envelope refusals (enable, components, hfreeze, dynamics, transport)")
   end subroutine test_step_defaults

   subroutine test_step_refusals(error)
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      character(len=*), parameter :: CICE = "&ocean_ice_nml scheme='cice' /"//new_line('a')// &
                                     "&ocean_ice_dyn_nml advection='upwind' /"//new_line('a')

      call parse(CICE, cfg, error)
      if (allocated(error)) return
      cfg%nghost = 2
      call check(error, has_msg(cfg, "nghost >= 3"), "nghost = 2 refused")
      if (allocated(error)) return
      call check(error, has_msg(cfg, "superseded"), "and the message says the old minimum of 2 is superseded")
      if (allocated(error)) return
      cfg%nghost = 3
      cfg%px = 2
      call check(error, has_msg(cfg, "single rank"), "px*py > 1 refused")
      if (allocated(error)) return
      cfg%px = 1
      call expect_refusal("&ocean_ice_nml scheme='cice' /"//new_line('a')// &
                          "&ocean_ice_dyn_nml advection='upwind', ndtd = 2 /", "ndtd /= 1", error)
      if (allocated(error)) return
      call expect_refusal("&ocean_ice_nml scheme='cice', dt_ice = 450.0 /"//new_line('a')// &
                          "&ocean_ice_dyn_nml advection='upwind' /"//new_line('a')//"&time_nml dt_fixed = 900.0 /", &
                          "dt_ice must equal the outer step", error)
      if (allocated(error)) return
      ! the equal value is accepted
      call parse("&ocean_ice_nml scheme='cice', dt_ice = 900.0 /"//new_line('a')// &
                 "&ocean_ice_dyn_nml advection='upwind' /"//new_line('a')//"&time_nml dt_fixed = 900.0 /", cfg, error)
      if (allocated(error)) return
      call check(error,.not. has_msg(cfg, "dt_ice"), "dt_ice = dt_fixed accepted")
      if (allocated(error)) return
      ! a non-finite dt_ice is refused
      call parse(CICE, cfg, error)
      if (allocated(error)) return
      cfg%ocean%ice%dt_ice = ieee_value(1.0_wp, ieee_quiet_nan)
      call check(error, has_msg(cfg, "dt_ice must be finite"), "a NaN dt_ice is refused")
      if (allocated(error)) return
      ! none of the envelope refusals applies under the default scheme
      call parse("&ocean_ice_dyn_nml ndtd = 2 /", cfg, error)
      if (allocated(error)) return
      cfg%nghost = 2
      call check(error, count_msgs(cfg), 0, "scheme='sis2' untouched by the cice envelope")
   end subroutine test_step_refusals

   subroutine test_step_all_reasons(error)
      !! Every envelope reason is reported together with the scheme refusal.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg

      call parse("&ocean_ice_nml scheme='cice', dt_ice = 100.0 /"//new_line('a')// &
                 "&ocean_ice_dyn_nml advection='upwind', ndtd = 3 /"//new_line('a')// &
                 "&time_nml dt_fixed = 900.0 /", cfg, error)
      if (allocated(error)) return
      cfg%nghost = 2
      cfg%px = 2
      call check(error, count_msgs(cfg), 9, "envelope x5 + nghost + ndtd + dt_ice + ranks")
      if (allocated(error)) return
      call check(error, has_msg(cfg, "needs enable=.true."), "envelope reason")
      if (allocated(error)) return
      call check(error, has_msg(cfg, "nghost >= 3") .and. has_msg(cfg, "ndtd /= 1") .and. has_msg(cfg, "dt_ice must equal") &
                 .and. has_msg(cfg, "single rank"), "every envelope reason")
   end subroutine test_step_all_reasons

   subroutine test_unrepresentable(error)
      !! The options the port excludes entirely (melt ponds, floe-size distribution, wave fracture, biogeochemistry,
      !! aerosols, snow physics, EAP/VP rheology) have NO knob: the namelist parser rejects each one, so they cannot be
      !! requested; those that do have a knob (mushy, delta-Eddington, `calc_Tsfc`, ...) are refused by the earlier tests.
      type(error_type), allocatable, intent(out) :: error

      call check(error, parse_fails("&ocean_ice_nml tr_pond_lvl = .true. /"), "melt ponds have no knob")
      if (allocated(error)) return
      call check(error, parse_fails("&ocean_ice_nml tr_fsd = .true. /"), "the floe size distribution has no knob")
      if (allocated(error)) return
      call check(error, parse_fails("&ocean_ice_nml wave_spec = .true. /"), "wave fracture has no knob")
      if (allocated(error)) return
      call check(error, parse_fails("&ocean_ice_nml tr_aero = .true. /"), "aerosols have no knob")
      if (allocated(error)) return
      call check(error, parse_fails("&ocean_ice_dyn_nml kdyn = 2 /"), "VP/EAP rheology has no knob")
      if (allocated(error)) return
      call check(error, parse_fails("&ponds_nml tr_pond_topo = .true. /"), "no ponds namelist")
   end subroutine test_unrepresentable

   subroutine test_dyn_round_trip(error)
      !! Each knob reaches `ice_params_t%dyn`.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      type(ice_params_t) :: p
      real(wp), parameter :: tol = 1.0e-15_wp

      call parse("&ocean_ice_dyn_nml ndte=240, elasticDamp=0.5, revised_evp=.true., arlx=150.0, brlx=250.0, "// &
                 "e_yieldcurve=1.5, e_plasticpot=2.5, Ktens=0.2, deltaminEVP=2.0e-9, capping_method='sum', "// &
                 "dyn_area_min=1.0e-6, dyn_mass_min=1.0e-4, rheo_area_min=5.0e-3, free_drift_approx=.false., "// &
                 "visc_method='avg_strength', cdn_ocn=0.0045, vel_cap=3.0 /", cfg, error)
      if (allocated(error)) return
      p = ice_params_from_config(cfg)
      call check(error, p%dyn%ndte, 240)
      if (allocated(error)) return
      call check(error, abs(p%dyn%elasticDamp - 0.5_wp) < tol, "elasticDamp")
      if (allocated(error)) return
      call check(error, p%dyn%revised_evp, .true., "revised_evp")
      if (allocated(error)) return
      call check(error, abs(p%dyn%arlx - 150.0_wp) < tol .and. abs(p%dyn%brlx - 250.0_wp) < tol, "arlx brlx")
      if (allocated(error)) return
      call check(error, abs(p%dyn%e_yieldcurve - 1.5_wp) < tol .and. abs(p%dyn%e_plasticpot - 2.5_wp) < tol, "e")
      if (allocated(error)) return
      call check(error, abs(p%dyn%Ktens - 0.2_wp) < tol, "Ktens")
      if (allocated(error)) return
      call check(error, abs(p%dyn%deltaminEVP - 2.0e-9_wp) < 1.0e-24_wp, "deltaminEVP")
      if (allocated(error)) return
      call check(error, p%dyn%capping_max, .false., "capping sum")
      if (allocated(error)) return
      call check(error, abs(p%dyn%dyn_area_min - 1.0e-6_wp) < 1.0e-21_wp, "dyn_area_min")
      if (allocated(error)) return
      call check(error, abs(p%dyn%dyn_mass_min - 1.0e-4_wp) < 1.0e-19_wp, "dyn_mass_min")
      if (allocated(error)) return
      call check(error, abs(p%dyn%rheo_area_min - 5.0e-3_wp) < 1.0e-18_wp, "rheo_area_min")
      if (allocated(error)) return
      call check(error, p%dyn%free_drift_approx, .false., "free_drift_approx")
      if (allocated(error)) return
      call check(error, p%dyn%visc_method, ICE_VISC_AVG_STRENGTH)
      if (allocated(error)) return
      call check(error, abs(p%dyn%cdn_ocn - 0.0045_wp) < tol .and. abs(p%dyn%vel_cap - 3.0_wp) < tol, "cdn_ocn vel_cap")
   end subroutine test_dyn_round_trip

end module test_cice_config
