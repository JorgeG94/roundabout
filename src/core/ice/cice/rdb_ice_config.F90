!! Sea-ice (CICE port) configuration glue: the configure-time refusals and the
!! translation of the namelist into the device-friendly `ice_params_t`.
!!
!! This module needs the config types, so it lives in the MAIN library, not in
!! `rdb_ice_core_objs` (the ice core depends on the shared infrastructure only).
!! The namelist groups are `&ocean_ice_nml` (structure: `scheme`, `nilyr`,
!! `nslyr`, `kcatbound`, `constants_set`), `&ocean_ice_thermo_nml` (thermo,
!! ITD and floe-geometry knobs), `&ocean_ice_sw_nml` (albedo parameters),
!! `&ocean_ice_ridge_nml` (ridging and strength) and `&ocean_ice_atm_nml`
!! (atmosphere-ice boundary layer; the option keys the port supports one value
!! of are refused here with a message naming the value); every default is the
!! upstream CICE code default EXCEPT `tfrz_option = "linear_salt"` (upstream's
!! code default is `mushy`, which is out of scope and refused; linear_salt is
!! OM3's value).  Every new knob is registered DEAD until the PR that consumes
!! it wires it (only `scheme` is live: it drives a refusal).
!!
!! **Refusals.**  `ice_config_check` collects every reason a configuration is
!! not supported (all of them, not just the first); `ice_config_validate`
!! logs them for `validate_config`, so each one is reported by
!! `./rdb --validate-only`.  Value-level checks (an unsupported thermodynamics
!! option, a non-positive `ahmax`) apply under either scheme; structure checks
!! (layer and category counts, `hlim`) apply under `scheme="cice"`.
!! `scheme="cice"` itself is refused until the port is runnable.
module rdb_ice_config
   use rdb_constants, only: wp
   use rdb_config, only: config_t
   use rdb_ice_dims, only: ICE_NCAT_MAX, ICE_NILYR_MAX
   use rdb_ice_params, only: ice_params_t, ice_params_default, ice_constants, ice_constants_t, &
                             ice_atm_params_finalize, ice_ridge_params_t, ice_ridge_params_default, &
                             ice_ridge_params_valid, &
                             ICE_CONSTANTS_ICEPACK, ICE_CONSTANTS_CESM, &
                             ICE_TFRZ_LINEAR_SALT, ICE_TFRZ_CONSTANT, ICE_TFRZ_MINUS1P8, &
                             ICE_CONDUCT_BUBBLY, ICE_CONDUCT_MU71
   use pic_logger, only: logger => global_logger
   implicit none
   private

   public :: ice_params_from_config, ice_config_check, ice_config_validate
   public :: ice_scheme_is_cice, ice_constants_set_from_name
   public :: ICE_CONFIG_MSG_LEN, ICE_CONFIG_MAX_MSGS

   integer, parameter :: ICE_CONFIG_MSG_LEN = 240
      !! Length of one refusal message.
   integer, parameter :: ICE_CONFIG_MAX_MSGS = 48
      !! Capacity of the refusal list (more than every check can raise).

contains

   pure function ice_scheme_is_cice(cfg) result(is_cice)
      !! `&ocean_ice_nml scheme = "cice"`.
      type(config_t), intent(in) :: cfg
      logical :: is_cice
      is_cice = trim(adjustl(cfg%ocean%ice%scheme)) == "cice"
   end function ice_scheme_is_cice

   pure function ice_constants_set_from_name(name) result(set)
      !! `ICE_CONSTANTS_*` for a `constants_set` spelling; the upstream set
      !! for anything unrecognised (the check reports it).
      character(len=*), intent(in) :: name
      integer :: set
      select case (trim(adjustl(name)))
      case ("cesm")
         set = ICE_CONSTANTS_CESM
      case default
         set = ICE_CONSTANTS_ICEPACK
      end select
   end function ice_constants_set_from_name

   pure subroutine add_msg(msgs, n, text)
      character(len=ICE_CONFIG_MSG_LEN), intent(inout) :: msgs(:)
      integer, intent(inout) :: n
      character(len=*), intent(in) :: text
      if (n < size(msgs)) then
         n = n + 1
         msgs(n) = text
      end if
   end subroutine add_msg

   pure subroutine ice_config_check(cfg, msgs, n)
      !! Every refusal for `cfg`, in `msgs(1:n)` (empty list: accepted).
      type(config_t), intent(in) :: cfg
      character(len=ICE_CONFIG_MSG_LEN), intent(out) :: msgs(ICE_CONFIG_MAX_MSGS)
      integer, intent(out) :: n
      character(len=24) :: s
      integer :: nilyr, nslyr, ncat, kcb
      real(wp) :: zref
      type(ice_ridge_params_t) :: rp
      character(len=24) :: szref
      type(ice_constants_t) :: cst

      msgs = ""
      n = 0

      ! ---- structure: scheme and the not-yet-runnable refusal ----
      if (ice_scheme_is_cice(cfg)) then
         call add_msg(msgs, n, "&ocean_ice_nml scheme='cice' is not runnable yet "// &
                      "(CICE6/Icepack port in progress); use scheme='sis2'")
         nilyr = cfg%ocean%ice%nilyr
         nslyr = cfg%ocean%ice%nslyr
         ncat = cfg%ocean%ice%ncat
         kcb = cfg%ocean%ice%kcatbound
         if (nilyr < 2) then
            call add_msg(msgs, n, "&ocean_ice_nml scheme='cice' requires nilyr >= 2 "// &
                         "(a single ice layer is not supported)")
         end if
         if (nilyr > ICE_NILYR_MAX) then
            write (s, "(i0)") ICE_NILYR_MAX
            call add_msg(msgs, n, "&ocean_ice_nml nilyr exceeds the compile-time "// &
                         "maximum ICE_NILYR_MAX="//trim(s))
         end if
         if (nslyr /= 1) then
            call add_msg(msgs, n, "&ocean_ice_nml nslyr /= 1 is not supported "// &
                         "(one snow layer only)")
         end if
         if (ncat < 2) then
            call add_msg(msgs, n, "&ocean_ice_nml scheme='cice' requires ncat >= 2 "// &
                         "(a single category is not supported)")
         end if
         if (ncat > ICE_NCAT_MAX) then
            write (s, "(i0)") ICE_NCAT_MAX
            call add_msg(msgs, n, "&ocean_ice_nml ncat exceeds the compile-time "// &
                         "maximum ICE_NCAT_MAX="//trim(s))
         end if
         if (kcb < 0 .or. kcb > 3) then
            call add_msg(msgs, n, "&ocean_ice_nml kcatbound must be 0..3 "// &
                         "(-1 is the single-category formula and needs ncat=1, "// &
                         "which is not supported)")
         end if
         if (kcb == 2 .and. (ncat < 5 .or. ncat > 7)) then
            call add_msg(msgs, n, "&ocean_ice_nml kcatbound=2 (WMO classes) requires "// &
                         "ncat of 5, 6 or 7")
         end if
         if (any(cfg%ocean%ice%hlim >= 0.0_wp)) then
            call add_msg(msgs, n, "&ocean_ice_nml hlim is a scheme='sis2' knob; "// &
                         "scheme='cice' takes its category bounds from kcatbound")
         end if
      end if

      ! ---- thermodynamics / ITD options the port does not carry ----
      if (cfg%ocean%ice_thermo%ktherm /= 1) then
         call add_msg(msgs, n, "&ocean_ice_thermo_nml ktherm /= 1 is not supported "// &
                      "(only Bitz-Lipscomb 1999; mushy thermodynamics is excluded)")
      end if
      if (trim(adjustl(cfg%ocean%ice_thermo%tfrz_option)) == "mushy") then
         call add_msg(msgs, n, "&ocean_ice_thermo_nml tfrz_option='mushy' is not "// &
                      "supported (mushy thermodynamics is excluded)")
      end if
      if (trim(adjustl(cfg%ocean%ice_thermo%saltflux_option)) /= "constant") then
         call add_msg(msgs, n, "&ocean_ice_thermo_nml saltflux_option /= 'constant' "// &
                      "is not supported (prognostic ice salinity is excluded)")
      end if
      if (trim(adjustl(cfg%ocean%ice_thermo%fbot_xfer_type)) /= "constant") then
         call add_msg(msgs, n, "&ocean_ice_thermo_nml fbot_xfer_type /= 'constant' "// &
                      "is not supported (only the constant transfer coefficient)")
      end if
      if (trim(adjustl(cfg%ocean%ice_thermo%cpl_frazil)) /= "fresh_ice_correction") then
         call add_msg(msgs, n, "&ocean_ice_thermo_nml cpl_frazil /= "// &
                      "'fresh_ice_correction' is not supported")
      end if
      if (cfg%ocean%ice_thermo%kitd /= 1) then
         call add_msg(msgs, n, "&ocean_ice_thermo_nml kitd /= 1 is not supported "// &
                      "(only the linear ITD remap; the delta-function ITD is excluded)")
      end if
      if (.not. cfg%ocean%ice_thermo%calc_tsfc) then
         call add_msg(msgs, n, "&ocean_ice_thermo_nml calc_Tsfc=.false. is not "// &
                      "supported (the surface temperature is always solved)")
      end if
      if (cfg%ocean%ice_thermo%semi_implicit_tsfc) then
         call add_msg(msgs, n, "&ocean_ice_thermo_nml semi_implicit_Tsfc=.true. is "// &
                      "not supported")
      end if
      if (.not. (cfg%ocean%ice_thermo%floediam > 0.0_wp)) then
         call add_msg(msgs, n, "&ocean_ice_thermo_nml floediam must be > 0 "// &
                      "(it divides the lateral-melt rate)")
      end if
      if (.not. (cfg%ocean%ice_thermo%floeshape > 0.0_wp)) then
         call add_msg(msgs, n, "&ocean_ice_thermo_nml floeshape must be > 0 "// &
                      "(it divides the lateral-melt rate)")
      end if

      ! ---- ridging and ice strength ----
      rp = ice_ridge_params_default()
      rp%kstrength = cfg%ocean%ice_ridge%kstrength
      if (.not. ice_ridge_params_valid(rp)) then
         call add_msg(msgs, n, "&ocean_ice_ridge_nml kstrength must be 0 (Hibler 1979) "// &
                      "or 1 (Rothrock 1975)")
      end if
      rp = ice_ridge_params_default()
      rp%krdg_partic = cfg%ocean%ice_ridge%krdg_partic
      if (.not. ice_ridge_params_valid(rp)) then
         call add_msg(msgs, n, "&ocean_ice_ridge_nml krdg_partic must be 0 (Thorndike "// &
                      "et al. 1975) or 1 (Lipscomb et al. 2007)")
      end if
      rp = ice_ridge_params_default()
      rp%krdg_redist = cfg%ocean%ice_ridge%krdg_redist
      if (.not. ice_ridge_params_valid(rp)) then
         call add_msg(msgs, n, "&ocean_ice_ridge_nml krdg_redist must be 0 (Hibler 1980) "// &
                      "or 1 (Lipscomb et al. 2007)")
      end if
      if (.not. (cfg%ocean%ice_ridge%mu_rdg > 0.0_wp)) then
         call add_msg(msgs, n, "&ocean_ice_ridge_nml mu_rdg must be > 0 (the "// &
                      "e-folding scale of the ridged-ice thickness)")
      end if
      if (.not. (cfg%ocean%ice_ridge%cf > 0.0_wp)) then
         call add_msg(msgs, n, "&ocean_ice_ridge_nml Cf must be > 0 (it divides "// &
                      "the ridging work)")
      end if
      if (.not. (cfg%ocean%ice_ridge%pstar > 0.0_wp)) then
         call add_msg(msgs, n, "&ocean_ice_ridge_nml Pstar must be > 0 (Hibler "// &
                      "strength constant)")
      end if
      if (.not. (cfg%ocean%ice_ridge%cstar > 0.0_wp)) then
         call add_msg(msgs, n, "&ocean_ice_ridge_nml Cstar must be > 0 (Hibler "// &
                      "concentration constant)")
      end if

      ! ---- atmosphere-ice boundary layer ----
      if (cfg%ocean%ice_atm%natmiter < 0) then
         call add_msg(msgs, n, "&ocean_ice_atm_nml natmiter must be >= 0")
      end if
      if (.not. (cfg%ocean%ice_atm%atmiter_conv >= 0.0_wp)) then
         call add_msg(msgs, n, "&ocean_ice_atm_nml atmiter_conv must be >= 0")
      end if
      cst = ice_constants(ice_constants_set_from_name(cfg%ocean%ice%constants_set))
      zref = cst%zref
      if (.not. (cfg%ocean%ice_atm%iceruf > 0.0_wp .and. cfg%ocean%ice_atm%iceruf < zref)) then
         write (szref, "(g0)") zref
         call add_msg(msgs, n, "&ocean_ice_atm_nml iceruf must be > 0 and below the "// &
                      "reference height zref="//trim(szref)//" m (the neutral drag "// &
                      "takes ln(zref/iceruf))")
      end if
      if (trim(adjustl(cfg%ocean%ice_atm%atmbndy)) /= "similarity") then
         call add_msg(msgs, n, "&ocean_ice_atm_nml atmbndy /= 'similarity' is not "// &
                      "supported (only the similarity-theory boundary layer)")
      end if
      if (cfg%ocean%ice_atm%formdrag) then
         call add_msg(msgs, n, "&ocean_ice_atm_nml formdrag=.true. is not supported "// &
                      "(form drag needs melt ponds and level ice, which are excluded)")
      end if
      if (.not. cfg%ocean%ice_atm%calc_strair) then
         call add_msg(msgs, n, "&ocean_ice_atm_nml calc_strair=.false. is not "// &
                      "supported (the wind stress is always computed in the ice model)")
      end if

      ! ---- shortwave albedo parameters ----
      if (cfg%ocean%ice_sw%ahmax <= 0.0_wp) then
         call add_msg(msgs, n, "&ocean_ice_sw_nml ahmax must be > 0 (it divides "// &
                      "the thin-ice albedo ramp)")
      end if
   end subroutine ice_config_check

   subroutine ice_config_validate(cfg, has_error)
      !! Log every refusal from `ice_config_check` and set `has_error`.
      type(config_t), intent(in) :: cfg
      logical, intent(inout) :: has_error
      character(len=ICE_CONFIG_MSG_LEN) :: msgs(ICE_CONFIG_MAX_MSGS)
      integer :: n, i

      call ice_config_check(cfg, msgs, n)
      do i = 1, n
         call logger%error(trim(msgs(i)))
      end do
      if (n > 0) has_error = .true.
   end subroutine ice_config_validate

   pure function ice_params_from_config(cfg) result(p)
      !! The ice parameter bundle for `cfg`.  Does not validate: call
      !! `ice_config_check` first.  Unrecognised enum spellings fall back to the
      !! upstream default.  A negative `snowpatch` keeps the value of the
      !! selected `constants_set`.
      type(config_t), intent(in) :: cfg
      type(ice_params_t) :: p
      integer :: set

      set = ice_constants_set_from_name(cfg%ocean%ice%constants_set)
      p = ice_params_default(set, cfg%ocean%ice%ncat, cfg%ocean%ice%nilyr, &
                             cfg%ocean%ice%nslyr)

      p%thermo%ktherm = cfg%ocean%ice_thermo%ktherm
      select case (trim(adjustl(cfg%ocean%ice_thermo%conduct)))
      case ("MU71")
         p%thermo%conduct = ICE_CONDUCT_MU71
      case default
         p%thermo%conduct = ICE_CONDUCT_BUBBLY
      end select
      select case (trim(adjustl(cfg%ocean%ice_thermo%tfrz_option)))
      case ("constant")
         p%thermo%tfrz_option = ICE_TFRZ_CONSTANT
      case ("minus1p8")
         p%thermo%tfrz_option = ICE_TFRZ_MINUS1P8
      case default
         p%thermo%tfrz_option = ICE_TFRZ_LINEAR_SALT
      end select
      p%thermo%ice_ref_salinity = cfg%ocean%ice_thermo%ice_ref_salinity
      p%thermo%congel_one_step = trim(adjustl(cfg%ocean%ice_thermo%congel_freeze)) == "one-step"
      p%thermo%update_ocn_f = cfg%ocean%ice_thermo%update_ocn_f
      p%thermo%hfrazilmin = cfg%ocean%ice_thermo%hfrazilmin
      p%thermo%hi_min = cfg%ocean%ice_thermo%hi_min
      p%thermo%ustar_min = cfg%ocean%ice_thermo%ustar_min
      p%thermo%floediam = cfg%ocean%ice_thermo%floediam
      p%thermo%floeshape = cfg%ocean%ice_thermo%floeshape
      p%thermo%emissivity = cfg%ocean%ice_thermo%emissivity
      p%thermo%calc_Tsfc = cfg%ocean%ice_thermo%calc_tsfc
      p%c%ksno = cfg%ocean%ice_thermo%ksno

      p%itd%kitd = cfg%ocean%ice_thermo%kitd
      p%itd%kcatbound = cfg%ocean%ice%kcatbound

      p%ridge%kstrength = cfg%ocean%ice_ridge%kstrength
      p%ridge%krdg_partic = cfg%ocean%ice_ridge%krdg_partic
      p%ridge%krdg_redist = cfg%ocean%ice_ridge%krdg_redist
      p%ridge%mu_rdg = cfg%ocean%ice_ridge%mu_rdg
      p%ridge%Cf = cfg%ocean%ice_ridge%cf
      p%ridge%Pstar = cfg%ocean%ice_ridge%pstar
      p%ridge%Cstar = cfg%ocean%ice_ridge%cstar

      p%atm%highfreq = cfg%ocean%ice_atm%highfreq
      p%atm%natmiter = cfg%ocean%ice_atm%natmiter
      p%atm%atmiter_conv = cfg%ocean%ice_atm%atmiter_conv
      p%atm%iceruf = cfg%ocean%ice_atm%iceruf
      call ice_atm_params_finalize(p%c, p%atm)

      p%sw%albicev = cfg%ocean%ice_sw%albicev
      p%sw%albicei = cfg%ocean%ice_sw%albicei
      p%sw%albsnowv = cfg%ocean%ice_sw%albsnowv
      p%sw%albsnowi = cfg%ocean%ice_sw%albsnowi
      p%sw%ahmax = cfg%ocean%ice_sw%ahmax
      if (cfg%ocean%ice_sw%snowpatch >= 0.0_wp) p%c%snowpatch = cfg%ocean%ice_sw%snowpatch
   end function ice_params_from_config

end module rdb_ice_config
