!! Sea-ice constants, parameter bundle, tracer index model, salinity profile
!! and freezing point (analytical test A1; category bounds are D5.1's).
!!
!! Expected values are transcribed from the Bitz and Lipscomb (1999) profile,
!! from the documented upstream parameter defaults, and from the CESM shared
!! constants, never computed by the module under test.
module test_cice_params
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_ice_dims, only: ICE_NCAT_MAX, ICE_NILYR_MAX, ICE_NSLYR_MAX, ICE_NTRCR_MAX
   use rdb_ice_params, only: ice_constants_t, ice_params_t, ice_constants, ice_params_default, &
                             ICE_CONSTANTS_ICEPACK, ICE_CONSTANTS_CESM, ICE_TFRZ_LINEAR_SALT, &
                             ICE_TFRZ_CONSTANT, ICE_TFRZ_MINUS1P8, ICE_CONDUCT_BUBBLY, &
                             ice_nt_Tsfc, ice_nt_qice, ice_nt_qsno, ice_nt_sice, ice_ntrcr, &
                             ice_trcr_depend, ice_trcr_base, ice_intensive, ice_salinity_profile, &
                             ice_sea_freezing_temperature
   use test_cice_util, only: ice_golden_check
   implicit none
   private

   public :: collect_cice_params_tests

contains

   subroutine collect_cice_params_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("dims_bounds", test_dims), &
                  new_unittest("constants_upstream_defaults", test_constants_icepack), &
                  new_unittest("constants_cesm_overrides", test_constants_cesm), &
                  new_unittest("params_default_bundle", test_params_default), &
                  new_unittest("tracer_indices", test_tracer_indices), &
                  new_unittest("tracer_depend_and_base", test_trcr_depend), &
                  new_unittest("intensive_from_extensive", test_intensive), &
                  new_unittest("salinity_profile_closed_form", test_salinity_profile), &
                  new_unittest("salinity_profile_properties", test_salinity_properties), &
                  new_unittest("sea_freezing_temperature", test_tfrz), &
                  new_unittest("refusals_not_guesses", test_refusals) &
                  ]
   end subroutine collect_cice_params_tests

   subroutine test_dims(error)
      type(error_type), allocatable, intent(out) :: error

      call check(error, ICE_NCAT_MAX, 8, "ncat max")
      if (allocated(error)) return
      call check(error, ICE_NILYR_MAX, 8, "nilyr max")
      if (allocated(error)) return
      call check(error, ICE_NSLYR_MAX, 1, "nslyr max")
      if (allocated(error)) return
      call check(error, ICE_NTRCR_MAX, 18, "ntrcr max")
      if (allocated(error)) return
      call check(error, ice_ntrcr(ICE_NILYR_MAX, ICE_NSLYR_MAX), ICE_NTRCR_MAX, "ntrcr(max) == max")
   end subroutine test_dims

   subroutine test_constants_icepack(error)
      !! Upstream code default values (documented in the Icepack parameter table).
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      real(wp) :: got(14), ref(14)

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      got = [c%rhoi, c%rhos, c%rhosi, c%rhow, c%rhofresh, c%cp_ice, c%cp_ocn, c%Lsub, c%Lvap, &
             c%depressT, c%Tffresh, c%Tocnfrz, c%kice, c%snowpatch]
      ref = [917.0_wp, 330.0_wp, 940.0_wp, 1026.0_wp, 1000.0_wp, 2106.0_wp, 4218.0_wp, 2.835e6_wp, &
             2.501e6_wp, 0.054_wp, 273.15_wp, -1.8_wp, 2.03_wp, 0.02_wp]
      call ice_golden_check(error, "icepack constants", got, ref, 0.0_wp, 0.0_wp)
      if (allocated(error)) return
      ! Derived: Lfresh = Lsub - Lvap = 3.34e5 J/kg.
      call check(error, abs(c%Lfresh - 3.34e5_wp) < 1.0e-6_wp, "Lfresh 3.34e5")
      if (allocated(error)) return
      ! BL99 salinity constants.
      call check(error, c%saltmax == 3.2_wp .and. c%min_salin == 0.1_wp, "saltmax/min_salin")
      if (allocated(error)) return
      call check(error, c%puny == 1.0e-11_wp .and. c%hs_min == 1.0e-4_wp .and. c%ksno == 0.3_wp, &
                 "puny/hs_min/ksno")
   end subroutine test_constants_icepack

   subroutine test_constants_cesm(error)
      !! The CESM-coupled set differs where the CESM shared constants differ.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c, u

      c = ice_constants(ICE_CONSTANTS_CESM)
      u = ice_constants(ICE_CONSTANTS_ICEPACK)
      call check(error, c%cp_ocn == 3996.0_wp, "cesm cp_ocn")
      if (allocated(error)) return
      call check(error, c%cp_ice == 2117.27_wp, "cesm cp_ice")
      if (allocated(error)) return
      call check(error, abs(c%Lfresh - 3.337e5_wp) < 1.0e-6_wp, "cesm Lfresh = 3.337e5")
      if (allocated(error)) return
      call check(error, abs(c%Tocnfrz + 34.0_wp*0.054_wp) < 1.0e-15_wp, "cesm Tocnfrz")
      if (allocated(error)) return
      call check(error, c%snowpatch == 0.005_wp, "cesm snowpatch")
      if (allocated(error)) return
      call check(error, abs(c%zvir - 0.60779_wp) < 1.0e-4_wp, "cesm zvir ~ 0.6078")
      if (allocated(error)) return
      ! Unchanged: densities, depressT, saltmax.
      call check(error, c%rhoi == u%rhoi .and. c%rhow == u%rhow .and. c%rhos == u%rhos .and. &
                 c%depressT == u%depressT .and. c%saltmax == u%saltmax, "shared values untouched")
      if (allocated(error)) return
      call check(error, abs(c%Lfresh - (c%Lsub - c%Lvap)) == 0.0_wp, "Lfresh derived")
   end subroutine test_constants_cesm

   subroutine test_params_default(error)
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p

      p = ice_params_default(ICE_CONSTANTS_ICEPACK, 5, 4, 1)
      call check(error, p%ncat == 5 .and. p%nilyr == 4 .and. p%nslyr == 1, "structure")
      if (allocated(error)) return
      call check(error, p%thermo%ktherm == 1 .and. p%thermo%conduct == ICE_CONDUCT_BUBBLY .and. &
                 p%thermo%tfrz_option == ICE_TFRZ_LINEAR_SALT, "thermo enums")
      if (allocated(error)) return
      call check(error, p%thermo%ice_ref_salinity == 4.0_wp .and. p%thermo%ustar_min == 0.005_wp .and. &
                 p%thermo%hi_min == 0.01_wp .and. p%thermo%hfrazilmin == 0.05_wp .and. &
                 p%thermo%emissivity == 0.985_wp .and. p%thermo%floediam == 300.0_wp, "upstream defaults")
      if (allocated(error)) return
      call check(error, p%thermo%l_brine .and. p%thermo%calc_Tsfc .and. .not. p%thermo%update_ocn_f, &
                 "thermo flags")
      if (allocated(error)) return
      call check(error, p%itd%kcatbound == 1 .and. p%itd%kitd == 1, "itd defaults")
      if (allocated(error)) return
      p = ice_params_default(ICE_CONSTANTS_CESM, 5, 4, 1)
      call check(error, p%constants_set == ICE_CONSTANTS_CESM .and. p%c%cp_ocn == 3996.0_wp, &
                 "cesm selector reaches the bundle")
   end subroutine test_params_default

   subroutine test_tracer_indices(error)
      !! Packed column Tsfc | qice(4) | qsno(1) | sice(4): 10 tracers.
      type(error_type), allocatable, intent(out) :: error

      call check(error, ice_nt_Tsfc(), 1, "nt_Tsfc")
      if (allocated(error)) return
      call check(error, ice_nt_qice(), 2, "nt_qice")
      if (allocated(error)) return
      call check(error, ice_nt_qsno(4), 6, "nt_qsno")
      if (allocated(error)) return
      call check(error, ice_nt_sice(4, 1), 7, "nt_sice")
      if (allocated(error)) return
      call check(error, ice_ntrcr(4, 1), 10, "ntrcr")
      if (allocated(error)) return
      call check(error, ice_ntrcr(7, 1), 16, "ntrcr nilyr=7")
   end subroutine test_tracer_indices

   subroutine test_trcr_depend(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: it, expect(10)
      real(wp) :: b(3)

      expect = [0, 1, 1, 1, 1, 2, 1, 1, 1, 1]
      do it = 1, 10
         call check(error, ice_trcr_depend(it, 4, 1), expect(it), "depend")
         if (allocated(error)) return
      end do
      call check(error, ice_trcr_depend(0, 4, 1), -1, "below range")
      if (allocated(error)) return
      call check(error, ice_trcr_depend(11, 4, 1), -1, "above range")
      if (allocated(error)) return
      b = ice_trcr_base(0)
      call check(error, all(b == [1.0_wp, 0.0_wp, 0.0_wp]), "base area")
      if (allocated(error)) return
      b = ice_trcr_base(1)
      call check(error, all(b == [0.0_wp, 1.0_wp, 0.0_wp]), "base ice volume")
      if (allocated(error)) return
      b = ice_trcr_base(2)
      call check(error, all(b == [0.0_wp, 0.0_wp, 1.0_wp]), "base snow volume")
   end subroutine test_trcr_depend

   subroutine test_intensive(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: PUNY = 1.0e-11_wp

      ! Area tracer: extensive a*Tsfc / a.
      call check(error, ice_intensive(-0.5_wp*0.4_wp, 0, 0.4_wp, 1.0_wp, 0.2_wp, PUNY, -1.8_wp) == -0.5_wp, &
                 "area tracer")
      if (allocated(error)) return
      ! Empty category: area tracer falls back to the supplied value (Tf for Tsfc).
      call check(error, ice_intensive(0.0_wp, 0, 0.0_wp, 0.0_wp, 0.0_wp, PUNY, -1.8_wp) == -1.8_wp, &
                 "empty area tracer -> fallback")
      if (allocated(error)) return
      ! Ice-volume tracer: q*v / v.
      call check(error, ice_intensive(-3.0e8_wp*0.5_wp, 1, 0.4_wp, 0.5_wp, 0.2_wp, PUNY, 0.0_wp) == -3.0e8_wp, &
                 "ice-volume tracer")
      if (allocated(error)) return
      ! Snow-volume tracer divides by the SNOW volume only.
      call check(error, ice_intensive(-1.0e8_wp*0.2_wp, 2, 0.4_wp, 0.5_wp, 0.2_wp, PUNY, 0.0_wp) == -1.0e8_wp, &
                 "snow-volume tracer")
      if (allocated(error)) return
      ! Volume tracers with a zero base give 0 (upstream tests > 0, no puny).
      call check(error, ice_intensive(5.0_wp, 1, 0.4_wp, 0.0_wp, 0.2_wp, PUNY, 9.0_wp) == 0.0_wp, &
                 "empty volume tracer -> 0")
   end subroutine test_intensive

   subroutine test_salinity_profile(error)
      !! Bitz and Lipscomb (1999) profile, evaluated here with the intrinsic
      !! functions as an independent reference (tolerance covers the
      !! reproducible library's ulp-level differences).
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NL = 4
      real(wp) :: s(NL + 1), ref(NL + 1), zn
      real(wp), parameter :: PI = 3.14159265358979323846_wp
      integer :: k

      call ice_salinity_profile(NL, 3.2_wp, 0.1_wp, s)
      do k = 1, NL
         zn = (real(k, wp) - 0.5_wp)/real(NL, wp)
         ref(k) = max(1.6_wp*(1.0_wp - cos(PI*zn**(0.407_wp/(0.573_wp + zn)))), 0.1_wp)
      end do
      ref(NL + 1) = 3.2_wp
      call ice_golden_check(error, "BL99 salinity (4 layers)", s, ref, 1.0e-13_wp, 1.0e-14_wp)
   end subroutine test_salinity_profile

   subroutine test_salinity_properties(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: s(ICE_NILYR_MAX + 1), z(5)
      integer :: nl, k

      do nl = 1, ICE_NILYR_MAX
         call ice_salinity_profile(nl, 3.2_wp, 0.1_wp, s(1:nl + 1))
         call check(error, s(nl + 1) == 3.2_wp, "base value is saltmax")
         if (allocated(error)) return
         do k = 1, nl
            call check(error, s(k) >= 0.1_wp .and. s(k) <= 3.2_wp, "bounded by [min_salin, saltmax]")
            if (allocated(error)) return
            if (k > 1) then
               call check(error, s(k) >= s(k - 1), "monotone with depth")
               if (allocated(error)) return
            end if
         end do
      end do
      ! saltmax <= min_salin: brine pockets off, profile identically zero.
      call ice_salinity_profile(4, 0.1_wp, 0.1_wp, z(1:5))
      call check(error, all(z == 0.0_wp), "brine off -> zero profile")
   end subroutine test_salinity_properties

   subroutine test_tfrz(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: tf
      logical :: bad

      call ice_sea_freezing_temperature(35.0_wp, ICE_TFRZ_LINEAR_SALT, 0.054_wp, -1.8_wp, tf, bad)
      call check(error,.not. bad .and. abs(tf - (-1.89_wp)) < 1.0e-14_wp, "linear_salt at 35 ppt = -1.89")
      if (allocated(error)) return
      call ice_sea_freezing_temperature(0.0_wp, ICE_TFRZ_LINEAR_SALT, 0.054_wp, -1.8_wp, tf, bad)
      call check(error, tf == 0.0_wp, "fresh water freezes at 0")
      if (allocated(error)) return
      call ice_sea_freezing_temperature(35.0_wp, ICE_TFRZ_CONSTANT, 0.054_wp, -1.836_wp, tf, bad)
      call check(error, tf == -1.836_wp, "constant option returns Tocnfrz")
      if (allocated(error)) return
      call ice_sea_freezing_temperature(35.0_wp, ICE_TFRZ_MINUS1P8, 0.054_wp, -1.836_wp, tf, bad)
      call check(error, tf == -1.8_wp, "minus1p8 option")
      if (allocated(error)) return
      call ice_sea_freezing_temperature(35.0_wp, 99, 0.054_wp, -1.836_wp, tf, bad)
      call check(error, bad, "unknown tfrz_option is refused, not guessed")
   end subroutine test_tfrz

   subroutine test_refusals(error)
      !! Unknown selectors and out-of-range structure set the flag.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      type(ice_constants_t) :: c

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      call check(error,.not. c%bad, "upstream set accepted")
      if (allocated(error)) return
      c = ice_constants(ICE_CONSTANTS_CESM)
      call check(error,.not. c%bad, "cesm set accepted")
      if (allocated(error)) return
      c = ice_constants(7)
      call check(error, c%bad, "unknown constants_set flagged")
      if (allocated(error)) return
      p = ice_params_default(ICE_CONSTANTS_ICEPACK, 5, 4, 1)
      call check(error,.not. p%bad, "OM3 structure accepted")
      if (allocated(error)) return
      p = ice_params_default(ICE_CONSTANTS_ICEPACK, 9, 4, 1)
      call check(error, p%bad, "ncat > ICE_NCAT_MAX flagged")
      if (allocated(error)) return
      p = ice_params_default(ICE_CONSTANTS_ICEPACK, 5, 4, 2)
      call check(error, p%bad, "nslyr > ICE_NSLYR_MAX flagged")
      if (allocated(error)) return
      p = ice_params_default(5, 5, 4, 1)
      call check(error, p%bad, "unknown constants_set flagged in the bundle")
   end subroutine test_refusals

end module test_cice_params
