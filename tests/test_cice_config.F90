!! Sea-ice (CICE port) configuration: analytical test A20.
!!
!! Every refusal fires with its message, the defaults are the upstream code
!! defaults (transcribed from the upstream parameter tables, not computed by
!! the module under test), both constants sets select correctly and
!! `ice_params_from_config` round-trips each knob.  This test links the main
!! library because the config types live there; it is still a `cice` test.
module test_cice_config
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_config, only: config_t, read_config_from_string, validate_config
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_CONFIG_VALIDATE
   use rdb_ice_dims, only: ICE_NCAT_MAX, ICE_NILYR_MAX
   use rdb_ice_params, only: ice_params_t, ICE_CONSTANTS_ICEPACK, ICE_CONSTANTS_CESM, &
                             ICE_TFRZ_LINEAR_SALT, ICE_TFRZ_CONSTANT, ICE_TFRZ_MINUS1P8, &
                             ICE_CONDUCT_BUBBLY, ICE_CONDUCT_MU71
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
                  new_unittest("params_round_trip", test_round_trip) &
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
                          "not runnable yet", error)
      if (allocated(error)) return
      call parse("&ocean_ice_nml enable = .true., scheme = 'cice' /", cfg, error)
      if (allocated(error)) return
      ierr = -999
      call validate_config(cfg, ierr=ierr)
      call check(error, ierr, OCEAN_STATUS_ERR_CONFIG_VALIDATE, &
                 "validate_config refuses scheme='cice'")
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
      call parse("&ocean_ice_nml scheme='cice', nilyr="//trim(s1)//", ncat="//trim(s2)//" /", &
                 cfg, error)
      if (allocated(error)) return
      call check(error, count_msgs(cfg), 1, "only the not-runnable reason at the maxima")
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
                 "&ocean_ice_sw_nml ahmax=0.0 /", cfg, error)
      if (allocated(error)) return
      call check(error, count_msgs(cfg), 6, "scheme + nilyr + nslyr + ktherm + kitd + ahmax")
      if (allocated(error)) return
      call check(error, has_msg(cfg, "not runnable yet"), "scheme reason")
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

end module test_cice_config
