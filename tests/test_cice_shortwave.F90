!! Analytic tests of the CCSM3 shortwave scalar kernels (`rdb_ice_shortwave`; MASTER_PLAN A14-1..3).
!!
!!   * `albedo_limits`         -- known albedos: cold/melting bare ice, cold/melting snow,
!!                                thin-ice ramp, floor, linear melt ramp, `ahmax`, code defaults.
!!   * `partition_closure`     -- property test over (hi, hs, Tsfc, SW, nilyr): `SW_in = reflected +
!!                                fswsfc + fswint + fswthru`, `sum Iswabs = fswint`, signs, night.
!!   * `beer_attenuation`      -- closed-form Beer law: layer ratios, refinement invariance, limits.
!!   * `four_band_reference`   -- an in-test 4-band (vdr, vdf, idr, idf) transcription of the
!!                                formulas against the 2-band production kernels on pseudo-random
!!                                inputs (proves the direct == diffuse collapse).
!!   * `oracle_golden`         -- 52 rows of the Icepack oracle (40 upstream-constant rows, 12 CESM-coupled
!!                                rows with snowpatch = 0.005), to a few ulp.
!!   * `device_kernel_closure` -- the kernels driven from a `do concurrent` loop over mapped
!!                                arrays (a GPU build maps them; inert on host), closure + host match.
!!
!! Expected values are closed forms written here with the INTRINSIC `atan`/`exp`, independent of
!! the safe-math bodies the kernels use (<= 2.5 ulp), hence tolerances of a few 1e-15 relative.
!! `oracle_golden` pins a subset of the Icepack column-driver oracle tables (embedded literals).
module test_cice_shortwave
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_ice_shortwave, only: ice_albedo_fhtan, ice_albedos_ccsm3, ice_absorbed_ccsm3, &
                                ICE_ALBOCN, ICE_KAPPAV, ICE_I0VIS, ICE_SNOWPATCH_DEFAULT, &
                                ICE_ALBICEV_DEFAULT, ICE_ALBICEI_DEFAULT, ICE_ALBSNOWV_DEFAULT, &
                                ICE_ALBSNOWI_DEFAULT, ICE_AHMAX_DEFAULT
   implicit none
   private

   public :: collect_cice_shortwave_tests

   ! OM3 shipped values (the code defaults are tested separately)
   real(wp), parameter :: OM3_V = 0.86_wp, OM3_I = 0.44_wp, OM3_SV = 0.98_wp, OM3_SI = 0.70_wp
   real(wp), parameter :: OM3_AHMAX = 0.1_wp
   real(wp), parameter :: TOL = 1.0e-14_wp

contains

   subroutine collect_cice_shortwave_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("albedo_limits", test_albedo_limits), &
                  new_unittest("partition_closure", test_partition_closure), &
                  new_unittest("beer_attenuation", test_beer_attenuation), &
                  new_unittest("four_band_reference", test_four_band_reference), &
                  new_unittest("oracle_golden", test_oracle_golden), &
                  new_unittest("device_kernel_closure", test_device_kernel_closure) &
                  ]
   end subroutine collect_cice_shortwave_tests

   subroutine albedos(hi, hs, tsfc, av, ai, sv, si, ahmax, alv_i, ali_i, alv_s, ali_s, asnow, alv, ali)
      !! Thin wrapper: computes `fhtan` from `ahmax` like a host caller would.
      real(wp), intent(in) :: hi, hs, tsfc, av, ai, sv, si, ahmax
      real(wp), intent(out) :: alv_i, ali_i, alv_s, ali_s, asnow, alv, ali
      call ice_albedos_ccsm3(hi, hs, tsfc, av, ai, sv, si, ice_albedo_fhtan(ahmax), ICE_SNOWPATCH_DEFAULT, &
                             alv_i, ali_i, alv_s, ali_s, asnow, alv, ali)
   end subroutine albedos

   subroutine test_albedo_limits(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: vi, ii, vs, is_, as, v, i, fh, e_v, e_i, a
      integer :: n

      ! cold, thick bare ice: the baseline albedos, exactly
      call albedos(1.0_wp, 0.0_wp, -10.0_wp, OM3_V, OM3_I, OM3_SV, OM3_SI, OM3_AHMAX, vi, ii, vs, is_, as, v, i)
      call check(error, vi, OM3_V, thr=TOL, message="cold bare ice vis")
      if (allocated(error)) return
      call check(error, ii, OM3_I, thr=TOL, message="cold bare ice nir")
      if (allocated(error)) return
      call check(error, as, 0.0_wp, thr=TOL, message="no snow -> asnow 0")
      if (allocated(error)) return
      call check(error, v, OM3_V, thr=TOL, message="category vis = ice vis without snow")
      if (allocated(error)) return

      ! melting bare ice (Tsfc = 0): 0.075 less in both bands
      call albedos(1.0_wp, 0.0_wp, 0.0_wp, OM3_V, OM3_I, OM3_SV, OM3_SI, OM3_AHMAX, vi, ii, vs, is_, as, v, i)
      call check(error, vi, 0.785_wp, thr=TOL, message="melting bare ice vis")
      if (allocated(error)) return
      call check(error, ii, 0.365_wp, thr=TOL, message="melting bare ice nir")
      if (allocated(error)) return

      ! the temperature ramp is linear between -1 and 0 degC: midpoint = half the reduction
      call albedos(1.0_wp, 0.0_wp, -0.5_wp, OM3_V, OM3_I, OM3_SV, OM3_SI, OM3_AHMAX, vi, ii, vs, is_, as, v, i)
      call check(error, vi, OM3_V - 0.0375_wp, thr=TOL, message="melt ramp midpoint vis")
      if (allocated(error)) return
      call albedos(1.0_wp, 0.0_wp, -1.0_wp, OM3_V, OM3_I, OM3_SV, OM3_SI, OM3_AHMAX, vi, ii, vs, is_, as, v, i)
      call check(error, vi, OM3_V, thr=TOL, message="ramp starts at -1 degC")
      if (allocated(error)) return
      call albedos(1.0_wp, 0.0_wp, 0.5_wp, OM3_V, OM3_I, OM3_SV, OM3_SI, OM3_AHMAX, vi, ii, vs, is_, as, v, i)
      call check(error, vi, OM3_V - 0.1125_wp, thr=TOL, message="ramp continues linearly above 0 degC (upstream form, no clamp)")
      if (allocated(error)) return

      ! cold snow on thick ice, hs = 0.5 m
      call albedos(1.0_wp, 0.5_wp, -10.0_wp, OM3_V, OM3_I, OM3_SV, OM3_SI, OM3_AHMAX, vi, ii, vs, is_, as, v, i)
      a = 0.5_wp/(0.5_wp + ICE_SNOWPATCH_DEFAULT)
      call check(error, as, a, thr=TOL, message="asnow = hs/(hs+snowpatch)")
      if (allocated(error)) return
      call check(error, vs, OM3_SV, thr=TOL, message="cold snow vis")
      if (allocated(error)) return
      call check(error, v, OM3_V*(1.0_wp - a) + OM3_SV*a, thr=TOL, message="cat vis with snow")
      if (allocated(error)) return
      call check(error, v, 0.9754_wp, thr=1.0e-4_wp, message="cat vis with snow (hand number)")
      if (allocated(error)) return
      call check(error, i, OM3_I*(1.0_wp - a) + OM3_SI*a, thr=TOL, message="cat nir with snow")
      if (allocated(error)) return

      ! melting snow, deep: asnow -> 1, vis 0.98-0.1, nir 0.70-0.15
      call albedos(1.0_wp, 1.0e3_wp, 0.0_wp, OM3_V, OM3_I, OM3_SV, OM3_SI, OM3_AHMAX, vi, ii, vs, is_, as, v, i)
      call check(error, vs, 0.88_wp, thr=TOL, message="melting snow vis")
      if (allocated(error)) return
      call check(error, is_, 0.55_wp, thr=TOL, message="melting snow nir")
      if (allocated(error)) return
      call check(error, v, 0.88_wp, thr=1.0e-4_wp, message="deep melting snow cat vis")
      if (allocated(error)) return

      ! thin ice ramp: hi = 0.05
      call albedos(0.05_wp, 0.0_wp, -10.0_wp, OM3_V, OM3_I, OM3_SV, OM3_SI, OM3_AHMAX, vi, ii, vs, is_, as, v, i)
      fh = atan(4.0_wp*0.05_wp)/atan(4.0_wp*OM3_AHMAX)
      e_v = OM3_V*fh + ICE_ALBOCN*(1.0_wp - fh)
      e_i = OM3_I*fh + ICE_ALBOCN*(1.0_wp - fh)
      call check(error, vi, e_v, thr=TOL, message="thin ice vis (closed form)")
      if (allocated(error)) return
      call check(error, ii, e_i, thr=TOL, message="thin ice nir (closed form)")
      if (allocated(error)) return
      call check(error, vi, 0.4751_wp, thr=1.0e-4_wp, message="thin ice vis (hand number)")
      if (allocated(error)) return
      call check(error, ii, 0.2572_wp, thr=1.0e-4_wp, message="thin ice nir (hand number)")
      if (allocated(error)) return

      ! hi -> 0: open-ocean albedo; melting thin ice never below it
      call albedos(0.0_wp, 0.0_wp, -10.0_wp, OM3_V, OM3_I, OM3_SV, OM3_SI, OM3_AHMAX, vi, ii, vs, is_, as, v, i)
      call check(error, vi, ICE_ALBOCN, thr=TOL, message="hi=0 vis -> albocn")
      if (allocated(error)) return
      call check(error, ii, ICE_ALBOCN, thr=TOL, message="hi=0 nir -> albocn")
      if (allocated(error)) return
      do n = 0, 20
         call albedos(1.0e-3_wp*real(n, wp), 0.0_wp, 0.0_wp, OM3_V, OM3_I, OM3_SV, OM3_SI, OM3_AHMAX, &
                      vi, ii, vs, is_, as, v, i)
         if (vi < ICE_ALBOCN .or. ii < ICE_ALBOCN) then
            call check(error, .false., message="melting thin ice albedo fell below albocn")
            return
         end if
      end do
      call albedos(1.0e-3_wp, 0.0_wp, 0.0_wp, OM3_V, OM3_I, OM3_SV, OM3_SI, OM3_AHMAX, vi, ii, vs, is_, as, v, i)
      call check(error, ii, ICE_ALBOCN, thr=TOL, message="thin melting nir is held at the floor")
      if (allocated(error)) return

      ! ahmax moves only the thin-ice ramp: thick ice unchanged, hi == ahmax gives fh = 1
      call albedos(2.0_wp, 0.0_wp, -10.0_wp, OM3_V, OM3_I, OM3_SV, OM3_SI, 0.5_wp, vi, ii, vs, is_, as, v, i)
      call check(error, vi, OM3_V, thr=TOL, message="thick ice independent of ahmax")
      if (allocated(error)) return
      call albedos(0.5_wp, 0.0_wp, -10.0_wp, OM3_V, OM3_I, OM3_SV, OM3_SI, 0.5_wp, vi, ii, vs, is_, as, v, i)
      call check(error, vi, OM3_V, thr=TOL, message="hi == ahmax -> fh = 1 exactly")
      if (allocated(error)) return
      call albedos(0.2_wp, 0.0_wp, -10.0_wp, OM3_V, OM3_I, OM3_SV, OM3_SI, 0.5_wp, vi, ii, vs, is_, as, v, i)
      fh = atan(0.8_wp)/atan(2.0_wp)
      call check(error, vi, OM3_V*fh + ICE_ALBOCN*(1.0_wp - fh), thr=TOL, message="ahmax = 0.5 ramp")
      if (allocated(error)) return

      ! upstream code defaults (0.78 / 0.36 / 0.98 / 0.70 / 0.3)
      call albedos(0.1_wp, 0.0_wp, -10.0_wp, ICE_ALBICEV_DEFAULT, ICE_ALBICEI_DEFAULT, ICE_ALBSNOWV_DEFAULT, &
                   ICE_ALBSNOWI_DEFAULT, ICE_AHMAX_DEFAULT, vi, ii, vs, is_, as, v, i)
      fh = atan(0.4_wp)/atan(4.0_wp*0.3_wp)
      call check(error, vi, 0.78_wp*fh + ICE_ALBOCN*(1.0_wp - fh), thr=TOL, message="code defaults vis")
      if (allocated(error)) return
      call check(error, ii, 0.36_wp*fh + ICE_ALBOCN*(1.0_wp - fh), thr=TOL, message="code defaults nir")
      if (allocated(error)) return
   end subroutine test_albedo_limits

   subroutine test_partition_closure(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: HI(7) = [0.0_wp, 0.01_wp, 0.1_wp, 0.5_wp, 1.0_wp, 3.0_wp, 8.0_wp]
      real(wp), parameter :: HS(4) = [0.0_wp, 0.01_wp, 0.1_wp, 0.5_wp]
      real(wp), parameter :: TS(3) = [-10.0_wp, -0.5_wp, 0.0_wp]
      real(wp), parameter :: SWV(4) = [100.0_wp, 0.0_wp, 250.0_wp, 0.0_wp]
      real(wp), parameter :: SWN(4) = [80.0_wp, 0.0_wp, 0.0_wp, 200.0_wp]
      integer, parameter :: NLY(4) = [1, 2, 4, 8]
      real(wp) :: vi, ii, vs, is_, as, v, i, fsfc, fint, fthru, isw(8), refl, sw_in, fpen
      real(wp) :: eps, res
      integer :: a, b, c, d, e, k, nl, nbad
      eps = epsilon(1.0_wp)
      nbad = 0
      do a = 1, size(HI)
         do b = 1, size(HS)
            do c = 1, size(TS)
               call albedos(HI(a), HS(b), TS(c), OM3_V, OM3_I, OM3_SV, OM3_SI, OM3_AHMAX, &
                            vi, ii, vs, is_, as, v, i)
               do d = 1, size(SWV)
                  do e = 1, size(NLY)
                     nl = NLY(e)
                     call ice_absorbed_ccsm3(nl, HI(a), as, SWV(d), SWN(d), vi, ii, vs, is_, &
                                             fsfc, fint, fthru, isw(1:nl))
                     sw_in = SWV(d) + SWN(d)
                     refl = SWV(d)*(vi*(1.0_wp - as) + vs*as) + SWN(d)*(ii*(1.0_wp - as) + is_*as)
                     fpen = SWV(d)*(1.0_wp - vi)*(1.0_wp - as)*ICE_I0VIS
                     res = abs(sw_in - refl - fsfc - fint - fthru)
                     ! 8 eps (plan: 4): the residual sums ~8 terms of magnitude SW_in, each with up to
                     ! 0.5 ulp rounding, so 4 eps is not a safe bound; 8 eps still catches any real leak.
                     ! Every comparison is written so that NaN FAILS (a bare `res > tol` is false for NaN).
                     if (.not. (res <= 8.0_wp*eps*sw_in)) nbad = nbad + 1
                     if (.not. (abs(sum(isw(1:nl)) - fint) <= 8.0_wp*eps*fpen)) nbad = nbad + 1
                     if (.not. (fsfc >= 0.0_wp .and. fint >= 0.0_wp .and. fthru >= 0.0_wp)) nbad = nbad + 1
                     do k = 1, nl
                        if (.not. (isw(k) >= 0.0_wp)) nbad = nbad + 1
                     end do
                     if (sw_in == 0.0_wp) then
                        if (.not. (fsfc == 0.0_wp .and. fint == 0.0_wp .and. fthru == 0.0_wp &
                                   .and. all(isw(1:nl) == 0.0_wp))) nbad = nbad + 1
                     end if
                     ! snow-covered fraction transmits nothing; NIR-only input transmits nothing
                     if (as >= 1.0_wp .and. .not. (fthru == 0.0_wp)) nbad = nbad + 1
                     if (SWV(d) == 0.0_wp .and. .not. (fthru == 0.0_wp)) nbad = nbad + 1
                  end do
               end do
            end do
         end do
      end do
      call check(error, nbad, 0, "energy partition closure / sign / night violations")
   end subroutine test_partition_closure

   subroutine test_beer_attenuation(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: fsfc, fint, fthru, isw(8), fpen, thru_ref(4), fint_ref(4), a
      integer :: nl, k, e
      integer, parameter :: NLY(4) = [1, 2, 4, 8]
      real(wp), parameter :: SWV = 100.0_wp
      fpen = SWV*(1.0_wp - OM3_V)*ICE_I0VIS
      call check(error, fpen, 9.8_wp, thr=1.0e-13_wp, message="fswpen = 100*0.14*0.7")
      if (allocated(error)) return
      do e = 1, 4
         nl = NLY(e)
         call ice_absorbed_ccsm3(nl, 1.0_wp, 0.0_wp, SWV, 0.0_wp, OM3_V, OM3_I, OM3_SV, OM3_SI, &
                                 fsfc, fint, fthru, isw(1:nl))
         thru_ref(e) = fthru
         fint_ref(e) = fint
         call check(error, fthru, 9.8_wp*exp(-ICE_KAPPAV), thr=1.0e-13_wp, message="fswthru = fswpen exp(-kappav hi)")
         if (allocated(error)) return
         ! layer closed form and constant ratio exp(-kappav hi/nilyr)
         do k = 1, nl
            a = 9.8_wp*(exp(-ICE_KAPPAV*real(k - 1, wp)/real(nl, wp)) - exp(-ICE_KAPPAV*real(k, wp)/real(nl, wp)))
            call check(error, isw(k), a, thr=1.0e-13_wp, message="Iswabs(k) closed form")
            if (allocated(error)) return
            if (k > 1) then
               call check(error, isw(k)/isw(k - 1), exp(-ICE_KAPPAV/real(nl, wp)), thr=1.0e-12_wp, &
                          message="layer-to-layer ratio")
               if (allocated(error)) return
            end if
         end do
         call check(error, sum(isw(1:nl)), 9.8_wp*(1.0_wp - exp(-ICE_KAPPAV)), thr=1.0e-13_wp, &
                    message="layers sum to the continuous integral")
         if (allocated(error)) return
      end do
      ! refinement: nilyr 1,2,4,8 give the same fswthru and fswint
      do e = 2, 4
         call check(error, thru_ref(e), thru_ref(1), thr=1.0e-13_wp, message="fswthru independent of nilyr")
         if (allocated(error)) return
         call check(error, fint_ref(e), fint_ref(1), thr=1.0e-13_wp, message="fswint independent of nilyr")
         if (allocated(error)) return
      end do
      ! hi -> 0: everything penetrating is transmitted
      call ice_absorbed_ccsm3(4, 1.0e-12_wp, 0.0_wp, SWV, 0.0_wp, OM3_V, OM3_I, OM3_SV, OM3_SI, &
                              fsfc, fint, fthru, isw(1:4))
      call check(error, fthru, fpen, thr=1.0e-10_wp, message="hi -> 0: fswthru -> fswpen")
      if (allocated(error)) return
      ! thick ice: nearly nothing transmitted
      call ice_absorbed_ccsm3(4, 8.0_wp, 0.0_wp, SWV, 0.0_wp, OM3_V, OM3_I, OM3_SV, OM3_SI, &
                              fsfc, fint, fthru, isw(1:4))
      if (.not. (fthru < 1.0e-4_wp*fpen)) then
         call check(error, .false., message="thick ice transmits more than 1e-4 of the penetrating flux")
         return
      end if
      ! snow-covered (asnow = 1) and NIR-only: no transmission; surface takes all
      call ice_absorbed_ccsm3(4, 1.0_wp, 1.0_wp, SWV, 0.0_wp, OM3_V, OM3_I, 0.9_wp, OM3_SI, &
                              fsfc, fint, fthru, isw(1:4))
      call check(error, fthru, 0.0_wp, thr=0.0_wp, message="snow-covered ice transmits nothing")
      if (allocated(error)) return
      call check(error, fsfc, SWV*0.1_wp, thr=1.0e-13_wp, message="snow-covered: all absorbed at the surface")
      if (allocated(error)) return
      call ice_absorbed_ccsm3(4, 1.0_wp, 0.0_wp, 0.0_wp, 50.0_wp, OM3_V, OM3_I, OM3_SV, OM3_SI, &
                              fsfc, fint, fthru, isw(1:4))
      call check(error, fthru, 0.0_wp, thr=0.0_wp, message="NIR-only transmits nothing")
      if (allocated(error)) return
      call check(error, fsfc, 50.0_wp*(1.0_wp - OM3_I), thr=1.0e-13_wp, message="NIR all absorbed at the surface")
   end subroutine test_beer_attenuation

   subroutine test_four_band_reference(error)
      !! A 4-band transcription of the formulas (one albedo set per band, direct and diffuse
      !! computed separately) against the 2-band production kernels.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NSAMP = 400
      real(wp) :: hi, hs, ts, vdr, vdf, idr, idf, v_sw, n_sw
      real(wp) :: vi, ii, vs, is_, as, v, i, fsfc, fint, fthru, isw(4)
      real(wp) :: r_alvdri, r_alvdfi, r_alidri, r_alidfi, r_alvdrs, r_alvdfs, r_alidrs, r_alidfs
      real(wp) :: r_as, r_swabs, r_pen, r_sfc, r_thru, r_int, r_tran, r_isw(4), r_cat_v
      integer :: s, k, nbad, nlyr
      nbad = 0
      nlyr = 4
      do s = 1, NSAMP
         hi = 4.0_wp*frac(s, 1)**2
         hs = merge(0.0_wp, 0.6_wp*frac(s, 2), frac(s, 2) < 0.25_wp)
         ts = merge(0.0_wp, -8.0_wp*frac(s, 3), frac(s, 3) > 0.9_wp)
         vdr = 150.0_wp*frac(s, 4)
         vdf = 90.0_wp*frac(s, 5)
         idr = 120.0_wp*frac(s, 6)
         idf = 70.0_wp*frac(s, 7)
         v_sw = vdr + vdf
         n_sw = idr + idf
         ! --- production (2-band) ---
         call albedos(hi, hs, ts, OM3_V, OM3_I, OM3_SV, OM3_SI, OM3_AHMAX, vi, ii, vs, is_, as, v, i)
         call ice_absorbed_ccsm3(nlyr, hi, as, v_sw, n_sw, vi, ii, vs, is_, fsfc, fint, fthru, isw)
         ! --- reference (4-band, intrinsics) ---
         call ref_albedo(hi, ts, OM3_V, r_alvdri)
         call ref_albedo(hi, ts, OM3_V, r_alvdfi)
         call ref_albedo(hi, ts, OM3_I, r_alidri)
         call ref_albedo(hi, ts, OM3_I, r_alidfi)
         r_alvdrs = OM3_SV - (-0.1_wp)*ref_ft(ts)
         r_alvdfs = r_alvdrs
         r_alidrs = OM3_SI - (-0.15_wp)*ref_ft(ts)
         r_alidfs = r_alidrs
         if (hs > 1.0e-11_wp) then
            r_as = hs/(hs + 0.02_wp)
         else
            r_as = 0.0_wp
            r_alvdrs = 0.06_wp; r_alvdfs = 0.06_wp; r_alidrs = 0.06_wp; r_alidfs = 0.06_wp
         end if
         r_swabs = vdr*((1.0_wp - r_alvdri)*(1.0_wp - r_as) + (1.0_wp - r_alvdrs)*r_as) &
                   + vdf*((1.0_wp - r_alvdfi)*(1.0_wp - r_as) + (1.0_wp - r_alvdfs)*r_as) &
                   + idr*((1.0_wp - r_alidri)*(1.0_wp - r_as) + (1.0_wp - r_alidrs)*r_as) &
                   + idf*((1.0_wp - r_alidfi)*(1.0_wp - r_as) + (1.0_wp - r_alidfs)*r_as)
         r_pen = vdr*(1.0_wp - r_alvdri)*(1.0_wp - r_as)*0.7_wp + vdf*(1.0_wp - r_alvdfi)*(1.0_wp - r_as)*0.7_wp
         r_sfc = r_swabs - r_pen
         r_tran = 1.0_wp
         do k = 1, nlyr
            r_isw(k) = r_pen*(r_tran - exp(-1.4_wp*hi/real(nlyr, wp)*real(k, wp)))
            r_tran = exp(-1.4_wp*hi/real(nlyr, wp)*real(k, wp))
         end do
         r_thru = r_pen*r_tran
         r_int = r_pen - r_thru
         r_cat_v = r_alvdri*(1.0_wp - r_as) + r_alvdrs*r_as
         if (.not. (abs(v - r_cat_v) <= 1.0e-14_wp)) nbad = nbad + 1
         if (.not. (abs(fsfc - r_sfc) <= 1.0e-12_wp*max(1.0_wp, r_swabs))) nbad = nbad + 1
         if (.not. (abs(fint - r_int) <= 1.0e-12_wp*max(1.0_wp, r_pen))) nbad = nbad + 1
         if (.not. (abs(fthru - r_thru) <= 1.0e-12_wp*max(1.0_wp, r_pen))) nbad = nbad + 1
         do k = 1, nlyr
            if (.not. (abs(isw(k) - r_isw(k)) <= 1.0e-12_wp*max(1.0_wp, r_pen))) nbad = nbad + 1
         end do
      end do
      call check(error, nbad, 0, "2-band kernels disagree with the 4-band reference")
   contains
      function frac(s_, d_) result(f)
         !! Deterministic quasi-random number in [0,1) (additive recurrence).
         integer, intent(in) :: s_, d_
         real(wp) :: f
         real(wp), parameter :: PHI(7) = [0.6180339887498949_wp, 0.7548776662466927_wp, 0.5698402909980532_wp, &
                                          0.4301597090019468_wp, 0.8191725133961645_wp, 0.3247179572447460_wp, &
                                          0.2207440846057595_wp]
         f = real(s_, wp)*PHI(d_) + 0.5_wp*real(d_, wp)
         f = f - real(floor(f), wp)
      end function frac
      function ref_ft(t) result(f)
         real(wp), intent(in) :: t
         real(wp) :: f
         f = min((0.0_wp - t)/1.0_wp - 1.0_wp, 0.0_wp)
      end function ref_ft
      subroutine ref_albedo(h, t, cold, al)
         real(wp), intent(in) :: h, t, cold
         real(wp), intent(out) :: al
         real(wp) :: fh
         fh = min(atan(4.0_wp*h)/atan(4.0_wp*OM3_AHMAX), 1.0_wp)
         al = cold*fh + 0.06_wp*(1.0_wp - fh)
         al = al - (-0.075_wp)*ref_ft(t)
         al = max(al, 0.06_wp)
      end subroutine ref_albedo
   end subroutine test_four_band_reference

   subroutine oracle_rows(g)
      !! 52 rows of the Icepack column-driver oracle (compute_albedos + absorbed_solar, nilyr = 4).
      !! Row 1-40: upstream constants (snowpatch 0.02); rows 41-52: CESM-coupled set (snowpatch 0.005).
      !! Columns: 1 par (1 = albicev/i 0.86/0.44, ahmax 0.1; 2 = 0.78/0.36, ahmax 0.3; snow 0.98/0.70),
      !! 2 hi, 3 hs, 4 Tsfc, 5-8 swvdr swvdf swidr swidf, 9-10 alvdfni alidfni, 11-12 alvdfns alidfns
      !! (-999 where hs <= puny), 13-14 alvdfn alidfn, 15-17 fswsfc fswint fswthru, 18-21 Iswabs(1:4).
      real(wp), intent(out) :: g(21, 52)
      g(:, 1) = [ &
                2.0_wp, 0.10000000000000001_wp, 0.0_wp, -30.0_wp, &
                100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                0.372724243940037_wp, 0.19030176830834875_wp, -999.0_wp, -999.0_wp, &
                0.372724243940037_wp, 0.19030176830834875_wp, 109.19723219186346_wp, 8.604583224638283_wp, &
                57.259371161657832_wp, 2.2653632947472846_wp, 2.1874470671990589_wp, 2.1122107358641475_wp, &
                2.039562126827795_wp]
      g(:, 2) = [ &
                2.0_wp, 0.5_wp, 9.9999999999999998e-13_wp, -0.5_wp, &
                100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                0.74250000000000005_wp, 0.32250000000000001_wp, -999.0_wp, -999.0_wp, &
                0.74250000000000005_wp, 0.32250000000000001_wp, 79.337500000000006_wp, 13.611074848739761_wp, &
                13.426425151260233_wp, 4.3406808009525548_wp, 3.6438149732777285_wp, 3.0588260617019505_wp, &
                2.5677530128075277_wp]
      g(:, 3) = [ &
                2.0_wp, 0.0_wp, 9.9999999999999998e-13_wp, -30.0_wp, &
                100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                0.059999999999999998_wp, 0.059999999999999998_wp, -999.0_wp, -999.0_wp, &
                0.059999999999999998_wp, 0.059999999999999998_wp, 136.30000000000001_wp, 0.0_wp, &
                98.699999999999989_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                0.0_wp]
      g(:, 4) = [ &
                1.0_wp, 0.050000000000000003_wp, 9.9999999999999998e-13_wp, -1.0_wp, &
                100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                0.47501656050634677_wp, 0.25713286624051473_wp, -999.0_wp, -999.0_wp, &
                0.47501656050634677_wp, 0.25713286624051473_wp, 97.910968153162912_wp, 3.7266731204642767_wp, &
                51.396588026369308_wp, 0.95626534374789562_wp, 0.93967627792290509_wp, 0.92337499530237466_wp, &
                0.90735650349110031_wp]
      g(:, 5) = [ &
                1.0_wp, 0.050000000000000003_wp, 0.0_wp, 0.0_wp, &
                100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                0.40001656050634676_wp, 0.18213286624051472_wp, -999.0_wp, -999.0_wp, &
                0.40001656050634676_wp, 0.18213286624051472_wp, 108.78596815316291_wp, 4.2590717887049365_wp, &
                58.739189358128648_wp, 1.0928789878854441_wp, 1.0739199807570206_wp, 1.0552898700162756_wp, &
                1.0369829500461925_wp]
      g(:, 6) = [ &
                1.0_wp, 1.0_wp, 9.9999999999999998e-13_wp, -0.5_wp, &
                100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                0.82250000000000001_wp, 0.40250000000000002_wp, -999.0_wp, -999.0_wp, &
                0.82250000000000001_wp, 0.40250000000000002_wp, 67.737499999999997_wp, 14.041549084538309_wp, &
                4.5959509154616907_wp, 5.5038757278674781_wp, 3.8785156727201264_wp, 2.7331438003532353_wp, &
                1.9260138835974685_wp]
      g(:, 7) = [ &
                2.0_wp, 0.01_wp, 9.9999999999999998e-13_wp, -30.0_wp, &
                100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                0.092857017533386044_wp, 0.073690423972244179_wp, -999.0_wp, -999.0_wp, &
                0.092857017533386044_wp, 0.073690423972244179_wp, 133.4523918137732_wp, 1.3242090919047484_wp, &
                93.925804067089715_wp, 0.3327923197714589_wp, 0.33162958262921394_wp, 0.33047090795350514_wp, &
                0.32931628155056497_wp]
      g(:, 8) = [ &
                1.0_wp, 0.01_wp, 0.0_wp, -5.0_wp, &
                100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                0.14405364961646189_wp, 0.099925483567819398_wp, -999.0_wp, -999.0_wp, &
                0.14405364961646189_wp, 0.099925483567819398_wp, 128.5250374104773_wp, 1.2494744062050671_wp, &
                88.624892384066413_wp, 0.31401044493505204_wp, 0.31291332944984618_wp, 0.31182004715684086_wp, &
                0.31073058466332543_wp]
      g(:, 9) = [ &
                1.0_wp, 0.0_wp, 0.0_wp, -0.5_wp, &
                100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                0.059999999999999998_wp, 0.059999999999999998_wp, -999.0_wp, -999.0_wp, &
                0.059999999999999998_wp, 0.059999999999999998_wp, 136.30000000000001_wp, 0.0_wp, &
                98.699999999999989_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                0.0_wp]
      g(:, 10) = [ &
                 1.0_wp, 0.050000000000000003_wp, 9.9999999999999998e-13_wp, -0.5_wp, &
                 100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                 0.4375165605063468_wp, 0.21963286624051473_wp, -999.0_wp, -999.0_wp, &
                 0.4375165605063468_wp, 0.21963286624051473_wp, 103.34846815316291_wp, 3.9928724545846066_wp, &
                 55.067888692248978_wp, 1.02457216581667_wp, 1.0067981293399628_wp, 0.98933243265932513_wp, &
                 0.97216972676864644_wp]
      g(:, 11) = [ &
                 1.0_wp, 0.10000000000000001_wp, 9.9999999999999998e-13_wp, 0.0_wp, &
                 100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                 0.78500000000000003_wp, 0.36499999999999999_wp, -999.0_wp, -999.0_wp, &
                 0.78500000000000003_wp, 0.36499999999999999_wp, 73.175000000000011_wp, 2.9492378358719584_wp, &
                 19.625762164128037_wp, 0.776457727985437_wp, 0.7497517876377805_wp, 0.72396438699183607_wp, &
                 0.69906393325690375_wp]
      g(:, 12) = [ &
                 2.0_wp, 0.050000000000000003_wp, 0.0_wp, -0.10000000000000001_wp, &
                 100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                 0.15473217513367743_wp, 0.060096739639032226_wp, -999.0_wp, -999.0_wp, &
                 0.15473217513367743_wp, 0.060096739639032226_wp, 132.02737815508129_wp, 6.0002595235400946_wp, &
                 82.752862087423765_wp, 1.5396682376960995_wp, 1.5129584359545847_wp, 1.4867119895591008_wp, &
                 1.4609208603303119_wp]
      g(:, 13) = [ &
                 2.0_wp, 0.5_wp, 0.0_wp, -0.5_wp, &
                 100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                 0.74250000000000005_wp, 0.32250000000000001_wp, -999.0_wp, -999.0_wp, &
                 0.74250000000000005_wp, 0.32250000000000001_wp, 79.337500000000006_wp, 13.611074848739761_wp, &
                 13.426425151260233_wp, 4.3406808009525548_wp, 3.6438149732777285_wp, 3.0588260617019505_wp, &
                 2.5677530128075277_wp]
      g(:, 14) = [ &
                 2.0_wp, 0.10000000000000001_wp, 0.0_wp, -5.0_wp, &
                 100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                 0.372724243940037_wp, 0.19030176830834875_wp, -999.0_wp, -999.0_wp, &
                 0.372724243940037_wp, 0.19030176830834875_wp, 109.19723219186346_wp, 8.604583224638283_wp, &
                 57.259371161657832_wp, 2.2653632947472846_wp, 2.1874470671990589_wp, 2.1122107358641475_wp, &
                 2.039562126827795_wp]
      g(:, 15) = [ &
                 1.0_wp, 0.0_wp, 9.9999999999999998e-13_wp, -30.0_wp, &
                 100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                 0.059999999999999998_wp, 0.059999999999999998_wp, -999.0_wp, -999.0_wp, &
                 0.059999999999999998_wp, 0.059999999999999998_wp, 136.30000000000001_wp, 0.0_wp, &
                 98.699999999999989_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                 0.0_wp]
      g(:, 16) = [ &
                 2.0_wp, 1.0_wp, 0.0_wp, -0.10000000000000001_wp, &
                 100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                 0.71250000000000002_wp, 0.29249999999999998_wp, -999.0_wp, -999.0_wp, &
                 0.71250000000000002_wp, 0.29249999999999998_wp, 83.6875_wp, 22.743354151012749_wp, &
                 7.4441458489872439_wp, 8.9147282916163366_wp, 6.2821028501804852_wp, 4.4269230569101694_wp, &
                 3.1195999523057583_wp]
      g(:, 17) = [ &
                 2.0_wp, 1.0_wp, 9.9999999999999998e-13_wp, 0.0_wp, &
                 100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                 0.70500000000000007_wp, 0.28499999999999998_wp, -999.0_wp, -999.0_wp, &
                 0.70500000000000007_wp, 0.28499999999999998_wp, 84.774999999999991_wp, 23.336659041908735_wp, &
                 7.63834095809126_wp, 9.1472864209628497_wp, 6.4459837940982379_wp, 4.5424080062208692_wp, &
                 3.200980820626778_wp]
      g(:, 18) = [ &
                 2.0_wp, 0.5_wp, 0.0_wp, -1.0_wp, &
                 100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                 0.78000000000000003_wp, 0.35999999999999999_wp, -999.0_wp, -999.0_wp, &
                 0.78000000000000003_wp, 0.35999999999999999_wp, 73.900000000000006_wp, 11.628879482418437_wp, &
                 11.471120517581557_wp, 3.7085428202313091_wp, 3.1131623072664087_wp, 2.6133659556288507_wp, &
                 2.1938083992918682_wp]
      g(:, 19) = [ &
                 2.0_wp, 0.5_wp, 9.9999999999999998e-13_wp, -1.0_wp, &
                 100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                 0.78000000000000003_wp, 0.35999999999999999_wp, -999.0_wp, -999.0_wp, &
                 0.78000000000000003_wp, 0.35999999999999999_wp, 73.900000000000006_wp, 11.628879482418437_wp, &
                 11.471120517581557_wp, 3.7085428202313091_wp, 3.1131623072664087_wp, 2.6133659556288507_wp, &
                 2.1938083992918682_wp]
      g(:, 20) = [ &
                 2.0_wp, 0.0_wp, 0.0_wp, -0.5_wp, &
                 100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                 0.059999999999999998_wp, 0.059999999999999998_wp, -999.0_wp, -999.0_wp, &
                 0.059999999999999998_wp, 0.059999999999999998_wp, 136.30000000000001_wp, 0.0_wp, &
                 98.699999999999989_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                 0.0_wp]
      g(:, 21) = [ &
                 1.0_wp, 0.0_wp, 0.0_wp, -5.0_wp, &
                 100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                 0.059999999999999998_wp, 0.059999999999999998_wp, -999.0_wp, -999.0_wp, &
                 0.059999999999999998_wp, 0.059999999999999998_wp, 136.30000000000001_wp, 0.0_wp, &
                 98.699999999999989_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                 0.0_wp]
      g(:, 22) = [ &
                 1.0_wp, 0.01_wp, 9.9999999999999998e-13_wp, -1.0_wp, &
                 100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                 0.14405364961646189_wp, 0.099925483567819398_wp, -999.0_wp, -999.0_wp, &
                 0.14405364961646189_wp, 0.099925483567819398_wp, 128.5250374104773_wp, 1.2494744062050671_wp, &
                 88.624892384066413_wp, 0.31401044493505204_wp, 0.31291332944984618_wp, 0.31182004715684086_wp, &
                 0.31073058466332543_wp]
      g(:, 23) = [ &
                 2.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, &
                 100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                 0.70500000000000007_wp, 0.28499999999999998_wp, -999.0_wp, -999.0_wp, &
                 0.70500000000000007_wp, 0.28499999999999998_wp, 84.774999999999991_wp, 23.336659041908735_wp, &
                 7.63834095809126_wp, 9.1472864209628497_wp, 6.4459837940982379_wp, 4.5424080062208692_wp, &
                 3.200980820626778_wp]
      g(:, 24) = [ &
                 2.0_wp, 0.0_wp, 9.9999999999999998e-13_wp, -0.10000000000000001_wp, &
                 100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                 0.059999999999999998_wp, 0.059999999999999998_wp, -999.0_wp, -999.0_wp, &
                 0.059999999999999998_wp, 0.059999999999999998_wp, 136.30000000000001_wp, 0.0_wp, &
                 98.699999999999989_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                 0.0_wp]
      g(:, 25) = [ &
                 1.0_wp, 1.0_wp, 9.9999999999999998e-13_wp, -5.0_wp, &
                 100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                 0.85999999999999999_wp, 0.44_wp, -999.0_wp, -999.0_wp, &
                 0.85999999999999999_wp, 0.44_wp, 62.299999999999997_wp, 11.075024630058385_wp, &
                 3.6249753699416156_wp, 4.3410850811349126_wp, 3.0591109531313676_wp, 2.1557190537997353_wp, &
                 1.5191095419923697_wp]
      g(:, 26) = [ &
                 1.0_wp, 0.0_wp, 0.0_wp, -0.10000000000000001_wp, &
                 100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                 0.059999999999999998_wp, 0.059999999999999998_wp, -999.0_wp, -999.0_wp, &
                 0.059999999999999998_wp, 0.059999999999999998_wp, 136.30000000000001_wp, 0.0_wp, &
                 98.699999999999989_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                 0.0_wp]
      g(:, 27) = [ &
                 2.0_wp, 0.29999999999999999_wp, 9.9999999999999998e-13_wp, -1.0_wp, &
                 100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                 0.78000000000000003_wp, 0.35999999999999999_wp, -999.0_wp, -999.0_wp, &
                 0.78000000000000003_wp, 0.35999999999999999_wp, 73.900000000000006_wp, 7.9222184622721876_wp, &
                 15.177781537727807_wp, 2.3025035282572639_wp, 2.0730003898314133_wp, 1.8663730862961088_wp, &
                 1.680341457887401_wp]
      g(:, 28) = [ &
                 1.0_wp, 2.5_wp, 0.0_wp, 0.0_wp, &
                 100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                 0.78500000000000003_wp, 0.36499999999999999_wp, -999.0_wp, -999.0_wp, &
                 0.78500000000000003_wp, 0.36499999999999999_wp, 73.175000000000011_wp, 21.893294069241154_wp, &
                 0.68170593075884001_wp, 13.16433990575767_wp, 5.4877133208485267_wp, 2.2876192583455719_wp, &
                 0.95362158428938626_wp]
      g(:, 29) = [ &
                 2.0_wp, 1.0_wp, 0.0_wp, -1.0_wp, &
                 100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                 0.78000000000000003_wp, 0.35999999999999999_wp, -999.0_wp, -999.0_wp, &
                 0.78000000000000003_wp, 0.35999999999999999_wp, 73.900000000000006_wp, 17.403610132948884_wp, &
                 5.6963898670511082_wp, 6.8217051274977178_wp, 4.8071743549207193_wp, 3.3875585131138686_wp, &
                 2.38717213741658_wp]
      g(:, 30) = [ &
                 1.0_wp, 0.5_wp, 9.9999999999999998e-13_wp, -0.5_wp, &
                 100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                 0.82250000000000001_wp, 0.40250000000000002_wp, -999.0_wp, -999.0_wp, &
                 0.82250000000000001_wp, 0.40250000000000002_wp, 67.737499999999997_wp, 9.3823914005876041_wp, &
                 9.2551085994123952_wp, 2.9921197754138977_wp, 2.5117559524535804_wp, 2.1085111687460052_wp, &
                 1.7700045039741215_wp]
      g(:, 31) = [ &
                 1.0_wp, 0.29999999999999999_wp, 9.9999999999999998e-13_wp, -5.0_wp, &
                 100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                 0.85999999999999999_wp, 0.44_wp, -999.0_wp, -999.0_wp, &
                 0.85999999999999999_wp, 0.44_wp, 62.299999999999997_wp, 5.0414117487186658_wp, &
                 9.6585882512813352_wp, 1.4652295179818957_wp, 1.3191820662563545_wp, 1.187691964006615_wp, &
                 1.069308200473801_wp]
      g(:, 32) = [ &
                 2.0_wp, 0.01_wp, 0.0_wp, 0.0_wp, &
                 100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                 0.059999999999999998_wp, 0.059999999999999998_wp, -999.0_wp, -999.0_wp, &
                 0.059999999999999998_wp, 0.059999999999999998_wp, 136.30000000000001_wp, 1.3721723812555382_wp, &
                 97.32782761874445_wp, 0.3448461671770518_wp, 0.34364131531264303_wp, 0.34244067305864601_wp, &
                 0.34124422570719032_wp]
      g(:, 33) = [ &
                 1.0_wp, 0.10000000000000001_wp, 9.9999999999999998e-13_wp, -30.0_wp, &
                 100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                 0.85999999999999999_wp, 0.44_wp, -999.0_wp, -999.0_wp, &
                 0.85999999999999999_wp, 0.44_wp, 62.299999999999997_wp, 1.9204339396375545_wp, &
                 12.779566060362447_wp, 0.50560038101377303_wp, 0.4882104663687874_wp, 0.47141867059933523_wp, &
                 0.45520442165565839_wp]
      g(:, 34) = [ &
                 1.0_wp, 2.5_wp, 0.0_wp, -30.0_wp, &
                 100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                 0.85999999999999999_wp, 0.44_wp, -999.0_wp, -999.0_wp, &
                 0.85999999999999999_wp, 0.44_wp, 62.299999999999997_wp, 14.256098463691918_wp, &
                 0.44390153630808199_wp, 8.5721283107259261_wp, 3.57339472055253_wp, 1.4896125403180474_wp, &
                 0.62096289209541444_wp]
      g(:, 35) = [ &
                 2.0_wp, 2.5_wp, 9.9999999999999998e-13_wp, -30.0_wp, &
                 100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                 0.78000000000000003_wp, 0.35999999999999999_wp, -999.0_wp, -999.0_wp, &
                 0.78000000000000003_wp, 0.35999999999999999_wp, 73.900000000000006_wp, 22.402440442944439_wp, &
                 0.69755955705555717_wp, 13.470487345426452_wp, 5.6153345608682601_wp, 2.3408197062140736_wp, &
                 0.97579883043565108_wp]
      g(:, 36) = [ &
                 1.0_wp, 0.5_wp, 9.9999999999999998e-13_wp, -1.0_wp, &
                 100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                 0.85999999999999999_wp, 0.44_wp, -999.0_wp, -999.0_wp, &
                 0.85999999999999999_wp, 0.44_wp, 62.299999999999997_wp, 7.4001960342662807_wp, &
                 7.2998039657337204_wp, 2.359981794692652_wp, 1.9811032864422609_wp, 1.6630510626729056_wp, &
                 1.396059890458462_wp]
      g(:, 37) = [ &
                 2.0_wp, 0.29999999999999999_wp, 0.0_wp, -0.10000000000000001_wp, &
                 100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                 0.71250000000000002_wp, 0.29249999999999998_wp, -999.0_wp, -999.0_wp, &
                 0.71250000000000002_wp, 0.29249999999999998_wp, 83.6875_wp, 10.352899126832973_wp, &
                 19.83460087316702_wp, 3.0089534744271065_wp, 2.70903460034787_wp, 2.4390102832278693_wp, &
                 2.1959007688301262_wp]
      g(:, 38) = [ &
                 1.0_wp, 0.10000000000000001_wp, 0.0_wp, -5.0_wp, &
                 100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                 0.85999999999999999_wp, 0.44_wp, -999.0_wp, -999.0_wp, &
                 0.85999999999999999_wp, 0.44_wp, 62.299999999999997_wp, 1.9204339396375545_wp, &
                 12.779566060362447_wp, 0.50560038101377303_wp, 0.4882104663687874_wp, 0.47141867059933523_wp, &
                 0.45520442165565839_wp]
      g(:, 39) = [ &
                 2.0_wp, 0.050000000000000003_wp, 0.0_wp, -0.5_wp, &
                 100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                 0.18473217513367743_wp, 0.090096739639032225_wp, -999.0_wp, -999.0_wp, &
                 0.18473217513367743_wp, 0.090096739639032225_wp, 127.67737815508129_wp, 5.7873000562438364_wp, &
                 79.815821554720031_wp, 1.4850227800410802_wp, 1.4592609548209385_wp, 1.4339460396735406_wp, &
                 1.4090702817082752_wp]
      g(:, 40) = [ &
                 1.0_wp, 0.0_wp, 0.0_wp, -30.0_wp, &
                 100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                 0.059999999999999998_wp, 0.059999999999999998_wp, -999.0_wp, -999.0_wp, &
                 0.059999999999999998_wp, 0.059999999999999998_wp, 136.30000000000001_wp, 0.0_wp, &
                 98.699999999999989_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                 0.0_wp]
      g(:, 41) = [ &
                 1.0_wp, 0.0_wp, 9.9999999999999998e-13_wp, -0.10000000000000001_wp, &
                 100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                 0.059999999999999998_wp, 0.059999999999999998_wp, -999.0_wp, -999.0_wp, &
                 0.059999999999999998_wp, 0.059999999999999998_wp, 136.30000000000001_wp, 0.0_wp, &
                 98.699999999999989_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                 0.0_wp]
      g(:, 42) = [ &
                 1.0_wp, 0.5_wp, 9.9999999999999998e-13_wp, 0.0_wp, &
                 100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                 0.78500000000000003_wp, 0.36499999999999999_wp, -999.0_wp, -999.0_wp, &
                 0.78500000000000003_wp, 0.36499999999999999_wp, 73.175000000000011_wp, 11.364586766908928_wp, &
                 11.210413233091067_wp, 3.624257756135143_wp, 3.0424086184648997_wp, 2.5539712748191046_wp, &
                 2.1439491174897807_wp]
      g(:, 43) = [ &
                 1.0_wp, 0.050000000000000003_wp, 9.9999999999999998e-13_wp, 0.0_wp, &
                 100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                 0.40001656050634676_wp, 0.18213286624051472_wp, -999.0_wp, -999.0_wp, &
                 0.40001656050634676_wp, 0.18213286624051472_wp, 108.78596815316291_wp, 4.2590717887049365_wp, &
                 58.739189358128648_wp, 1.0928789878854441_wp, 1.0739199807570206_wp, 1.0552898700162756_wp, &
                 1.0369829500461925_wp]
      g(:, 44) = [ &
                 2.0_wp, 0.01_wp, 9.9999999999999998e-13_wp, -0.5_wp, &
                 100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                 0.059999999999999998_wp, 0.059999999999999998_wp, -999.0_wp, -999.0_wp, &
                 0.059999999999999998_wp, 0.059999999999999998_wp, 136.30000000000001_wp, 1.3721723812555382_wp, &
                 97.32782761874445_wp, 0.3448461671770518_wp, 0.34364131531264303_wp, 0.34244067305864601_wp, &
                 0.34124422570719032_wp]
      g(:, 45) = [ &
                 1.0_wp, 0.0_wp, 9.9999999999999998e-13_wp, -1.0_wp, &
                 100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                 0.059999999999999998_wp, 0.059999999999999998_wp, -999.0_wp, -999.0_wp, &
                 0.059999999999999998_wp, 0.059999999999999998_wp, 136.30000000000001_wp, 0.0_wp, &
                 98.699999999999989_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                 0.0_wp]
      g(:, 46) = [ &
                 2.0_wp, 0.10000000000000001_wp, 9.9999999999999998e-13_wp, -5.0_wp, &
                 100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                 0.372724243940037_wp, 0.19030176830834875_wp, -999.0_wp, -999.0_wp, &
                 0.372724243940037_wp, 0.19030176830834875_wp, 109.19723219186346_wp, 8.604583224638283_wp, &
                 57.259371161657832_wp, 2.2653632947472846_wp, 2.1874470671990589_wp, 2.1122107358641475_wp, &
                 2.039562126827795_wp]
      g(:, 47) = [ &
                 1.0_wp, 0.29999999999999999_wp, 9.9999999999999998e-13_wp, -0.10000000000000001_wp, &
                 100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                 0.79249999999999998_wp, 0.3725_wp, -999.0_wp, -999.0_wp, &
                 0.79249999999999998_wp, 0.3725_wp, 72.087500000000006_wp, 7.4720924132794497_wp, &
                 14.315407586720548_wp, 2.1716794641517381_wp, 1.9552162767728107_wp, 1.7603291609383755_wp, &
                 1.5848675114165263_wp]
      g(:, 48) = [ &
                 2.0_wp, 0.01_wp, 9.9999999999999998e-13_wp, -30.0_wp, &
                 100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                 0.092857017533386044_wp, 0.073690423972244179_wp, -999.0_wp, -999.0_wp, &
                 0.092857017533386044_wp, 0.073690423972244179_wp, 133.4523918137732_wp, 1.3242090919047484_wp, &
                 93.925804067089715_wp, 0.3327923197714589_wp, 0.33162958262921394_wp, 0.33047090795350514_wp, &
                 0.32931628155056497_wp]
      g(:, 49) = [ &
                 1.0_wp, 2.5_wp, 9.9999999999999998e-13_wp, 0.0_wp, &
                 100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                 0.78500000000000003_wp, 0.36499999999999999_wp, -999.0_wp, -999.0_wp, &
                 0.78500000000000003_wp, 0.36499999999999999_wp, 73.175000000000011_wp, 21.893294069241154_wp, &
                 0.68170593075884001_wp, 13.16433990575767_wp, 5.4877133208485267_wp, 2.2876192583455719_wp, &
                 0.95362158428938626_wp]
      g(:, 50) = [ &
                 2.0_wp, 1.0_wp, 9.9999999999999998e-13_wp, -0.5_wp, &
                 100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                 0.74250000000000005_wp, 0.32250000000000001_wp, -999.0_wp, -999.0_wp, &
                 0.74250000000000005_wp, 0.32250000000000001_wp, 79.337500000000006_wp, 20.370134587428812_wp, &
                 6.6673654125711836_wp, 7.9844957742302833_wp, 5.6265790745094781_wp, 3.9649832596673691_wp, &
                 2.7940764790216792_wp]
      g(:, 51) = [ &
                 2.0_wp, 0.0_wp, 9.9999999999999998e-13_wp, -0.5_wp, &
                 100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                 0.059999999999999998_wp, 0.059999999999999998_wp, -999.0_wp, -999.0_wp, &
                 0.059999999999999998_wp, 0.059999999999999998_wp, 136.30000000000001_wp, 0.0_wp, &
                 98.699999999999989_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                 0.0_wp]
      g(:, 52) = [ &
                 2.0_wp, 0.0_wp, 9.9999999999999998e-13_wp, -1.0_wp, &
                 100.0_wp, 50.0_wp, 60.0_wp, 40.0_wp, &
                 0.059999999999999998_wp, 0.059999999999999998_wp, -999.0_wp, -999.0_wp, &
                 0.059999999999999998_wp, 0.059999999999999998_wp, 136.30000000000001_wp, 0.0_wp, &
                 98.699999999999989_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                 0.0_wp]
   end subroutine oracle_rows

   subroutine test_oracle_golden(error)
      !! Compare the kernels with the oracle rows.  Tolerance 8e-15 relative (about 36 eps): the
      !! oracle uses libm `atan`/`exp`, the port the library bodies (<= 2.5 ulp each), and the
      !! absorbed fluxes accumulate a few such errors; the reviewer measured <= 1.7e-15.
      !! Rows with hs <= puny have undefined snow albedos upstream (-999): there the port is
      !! asserted to return `albocn` and `asnow = 0` instead.
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: RTOL = 8.0e-15_wp
      real(wp) :: g(21, 52), av, ai, ahmax, sp, vi, ii, vs, is_, as, v, i, fsfc, fint, fthru, isw(4)
      integer :: r, k, nbad
      call oracle_rows(g)
      nbad = 0
      do r = 1, 52
         if (nint(g(1, r)) == 1) then
            av = 0.86_wp; ai = 0.44_wp; ahmax = 0.1_wp
         else
            av = 0.78_wp; ai = 0.36_wp; ahmax = 0.3_wp
         end if
         sp = merge(0.02_wp, 0.005_wp, r <= 40)
         call ice_albedos_ccsm3(g(2, r), g(3, r), g(4, r), av, ai, 0.98_wp, 0.70_wp, ice_albedo_fhtan(ahmax), sp, &
                                vi, ii, vs, is_, as, v, i)
         call ice_absorbed_ccsm3(4, g(2, r), as, g(5, r) + g(6, r), g(7, r) + g(8, r), vi, ii, vs, is_, &
                                 fsfc, fint, fthru, isw)
         if (.not. near(vi, g(9, r)) .or. .not. near(ii, g(10, r))) nbad = nbad + 1
         if (g(3, r) > 1.0e-11_wp) then
            if (.not. near(vs, g(11, r)) .or. .not. near(is_, g(12, r))) nbad = nbad + 1
         else
            if (.not. (vs == ICE_ALBOCN .and. is_ == ICE_ALBOCN .and. as == 0.0_wp)) nbad = nbad + 1
         end if
         if (.not. near(v, g(13, r)) .or. .not. near(i, g(14, r))) nbad = nbad + 1
         if (.not. near(fsfc, g(15, r)) .or. .not. near(fint, g(16, r)) .or. .not. near(fthru, g(17, r))) &
            nbad = nbad + 1
         do k = 1, 4
            if (.not. near(isw(k), g(17 + k, r))) nbad = nbad + 1
         end do
      end do
      call check(error, nbad, 0, "kernels disagree with the Icepack oracle rows")
   contains
      logical function near(x, ref)
         real(wp), intent(in) :: x, ref
         ! relative with an absolute floor for exact zeros; NaN -> .false.
         near = abs(x - ref) <= RTOL*max(abs(ref), 1.0_wp) .or. x == ref
      end function near
   end subroutine test_oracle_golden

   subroutine test_device_kernel_closure(error)
      !! The kernels from a `do concurrent` loop over mapped arrays: every array the kernel touches
      !! is device-present (`mem:separate`), results come back by `update self` of each array.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: N = 2000, NL = 4
      real(wp) :: hi(N), hs(N), ts(N), swv(N), swn(N)
      real(wp) :: fsfc(N), fint(N), fthru(N), isw(NL, N), alv(N)
      real(wp) :: fhtan, reflected, sw_in, tolh
      integer :: nbad_close, nbad_layers, nbad_host
      real(wp) :: vi, ii, vs, is_, as, v, i, hf, hn, ht, h_isw(NL)
      integer :: m
      do m = 1, N
         hi(m) = 5.0_wp*mod(real(m, wp)*0.6180339887498949_wp, 1.0_wp)**2
         hs(m) = merge(0.0_wp, 0.5_wp*mod(real(m, wp)*0.7548776662466927_wp, 1.0_wp), mod(m, 3) == 0)
         ts(m) = -6.0_wp*mod(real(m, wp)*0.5698402909980532_wp, 1.0_wp)
         swv(m) = 200.0_wp*mod(real(m, wp)*0.4301597090019468_wp, 1.0_wp)
         swn(m) = 150.0_wp*mod(real(m, wp)*0.8191725133961645_wp, 1.0_wp)
      end do
      fsfc = 0.0_wp; fint = 0.0_wp; fthru = 0.0_wp; isw = 0.0_wp; alv = 0.0_wp
      fhtan = ice_albedo_fhtan(OM3_AHMAX)
      !$acc enter data copyin(hi, hs, ts, swv, swn) create(fsfc, fint, fthru, isw, alv)
      call run_kernel(N, NL, fhtan, hi, hs, ts, swv, swn, fsfc, fint, fthru, isw, alv)
      !$acc update self(fsfc, fint, fthru, isw, alv)
      !$acc exit data delete(hi, hs, ts, swv, swn, fsfc, fint, fthru, isw, alv)
      nbad_close = 0
      nbad_layers = 0
      nbad_host = 0
      do m = 1, N
         call ice_albedos_ccsm3(hi(m), hs(m), ts(m), OM3_V, OM3_I, OM3_SV, OM3_SI, fhtan, &
                                ICE_SNOWPATCH_DEFAULT, vi, ii, vs, is_, as, v, i)
         sw_in = swv(m) + swn(m)
         reflected = swv(m)*(vi*(1.0_wp - as) + vs*as) + swn(m)*(ii*(1.0_wp - as) + is_*as)
         ! NaN-proof: every test is written so a NaN counts as a failure
         if (.not. (abs(sw_in - reflected - fsfc(m) - fint(m) - fthru(m)) <= 1.0e-13_wp*max(sw_in, 1.0_wp))) &
            nbad_close = nbad_close + 1
         if (.not. (abs(sum(isw(:, m)) - fint(m)) <= 1.0e-13_wp*max(sw_in, 1.0_wp))) nbad_layers = nbad_layers + 1
         call ice_absorbed_ccsm3(NL, hi(m), as, swv(m), swn(m), vi, ii, vs, is_, hf, hn, ht, h_isw)
         tolh = 1.0e-12_wp*max(sw_in, 1.0_wp)
         if (.not. (abs(hf - fsfc(m)) <= tolh .and. abs(hn - fint(m)) <= tolh .and. &
                    abs(ht - fthru(m)) <= tolh .and. abs(v - alv(m)) <= 1.0e-12_wp)) nbad_host = nbad_host + 1
      end do
      call check(error, nbad_close, 0, "device-kernel energy partition closure")
      if (allocated(error)) return
      call check(error, nbad_layers, 0, "device-kernel sum(Iswabs) == fswint")
      if (allocated(error)) return
      call check(error, nbad_host, 0, "device-kernel results match host scalar calls")
   end subroutine test_device_kernel_closure

   subroutine run_kernel(n, nl, fhtan, hi, hs, ts, swv, swn, fsfc, fint, fthru, isw, alv)
      !! One `do concurrent` over the cells, explicit-shape dummies (project GPU rule).
      integer, intent(in) :: n, nl
      real(wp), intent(in) :: fhtan
      real(wp), intent(in) :: hi(n), hs(n), ts(n), swv(n), swn(n)
      real(wp), intent(inout) :: fsfc(n), fint(n), fthru(n), isw(nl, n), alv(n)
      real(wp) :: vi, ii, vs, is_, as, v, i
      integer :: m
      do concurrent(m=1:n) local(vi, ii, vs, is_, as, v, i)
         call ice_albedos_ccsm3(hi(m), hs(m), ts(m), OM3_V, OM3_I, OM3_SV, OM3_SI, fhtan, ICE_SNOWPATCH_DEFAULT, &
                                vi, ii, vs, is_, as, v, i)
         alv(m) = v
         call ice_absorbed_ccsm3(nl, hi(m), as, swv(m), swn(m), vi, ii, vs, is_, &
                                 fsfc(m), fint(m), fthru(m), isw(1, m))
      end do
   end subroutine run_kernel

end module test_cice_shortwave
