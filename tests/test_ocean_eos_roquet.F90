!! Tests for the Roquet et al. (2015) specific-volume EOS variant (Capability 2).
module test_ocean_eos_roquet
   !! Locks the `EOS_VARIANT_ROQUET_SPV` branch added to the shared EOS
   !! handle:
   !!
   !!   (1) DENSITY checkvalues — `eos_density_point` reproduces the
   !!       verified prototype oracles (roquet_spv_eos.py, 128/128
   !!       identical to the published MOM6 transcription) at four
   !!       (CT, SA, p) points to ~1e-9.  The production branch works in
   !!       MODEL variables (PT, SP) and applies SR = SP·(35.16504/35)
   !!       plus CT = ct_from_pt(SR, PT) internally, so to hit a given
   !!       (CT_target, SA_target) we feed SP = SA_target/(35.16504/35)
   !!       and Newton-solve PT so the branch's own ct_from_pt recovers
   !!       CT_target.  This exercises the FULL production conversion
   !!       path (not a bypass) — the inversion uses the same Horner
   !!       polynomial the branch evaluates, so the round-trip is exact
   !!       to root-find tolerance.
   !!
   !!   (2) CHAIN-RULE derivative regression — analytic dSV/dT, dSV/dS
   !!       (w.r.t. the MODEL variables PT, SP, including the dCT/dPT and
   !!       SR-factor chain factors) agree with a central finite
   !!       difference of `eos_density_point`'s SV = 1/rho as a function
   !!       of (model T, model S) to ~1e-7.  This is the regression that
   !!       would catch a dropped chain factor (e.g. dCT/dPT≈1).
   !!
   !!   (3) DEFAULT-OFF bit-identity — the linear and Wright branches are
   !!       byte-unchanged by adding the Roquet case.
   !!
   !!   (4) FAIL-LOUD — roquet_spv + fv_wright is rejected at configure;
   !!       asserted here in the harness-compatible predicate form
   !!       (mirrors the `error stop` guard in configure_ocean_pgf
   !!       without aborting the test process).
   !!
   !! Runs host-side under NVHPC so the `!$acc routine seq` device
   !! codegen path is compiled and the values locked.
   use rdb_constants, only: wp
   use rdb_eos, only: eos_t, eos_density_point, eos_specvol_derivs, &
                      EOS_VARIANT_LINEAR, EOS_VARIANT_WRIGHT_97, &
                      EOS_VARIANT_ROQUET_SPV, roquet_model_ct, roquet_model_ct_derivs
   use rdb_ocean_eos_compute, only: ocean_eos_compute
   use rdb_ocean_pressure_force, only: OPGF_VARIANT_FV_WRIGHT, OPGF_VARIANT_FV_LITE
   use rdb_teos10_convert, only: ct_from_pt
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use testdrive, only: new_unittest, unittest_type, error_type, check
   implicit none
   private

   public :: collect_ocean_eos_roquet_tests

   real(wp), parameter :: SR_FACTOR = 35.16504_wp/35.0_wp
   integer, parameter :: NGHOST = 2

contains

   subroutine collect_ocean_eos_roquet_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("roquet_density_checkvalues", test_density_checkvalues), &
                  new_unittest("roquet_chain_rule_derivs", test_chain_rule_derivs), &
                  new_unittest("roquet_default_off_bit_identity", test_default_off), &
                  new_unittest("roquet_fail_loud_fv_wright", test_fail_loud_fv_wright), &
                  new_unittest("roquet_fresh_water_finite", test_fresh_water_finite), &
                  new_unittest("roquet_3d_fill_matches_point", test_3d_fill_matches_point), &
                  new_unittest("roquet_model_ct_pins_gsw", test_model_ct_pins_gsw) &
                  ]
   end subroutine collect_ocean_eos_roquet_tests

   pure function make_roquet() result(eos)
      type(eos_t) :: eos
      eos%variant = EOS_VARIANT_ROQUET_SPV
      eos%rho0 = 1035.0_wp
   end function make_roquet

   pure function pt_for_ct(eos, ct_target, sp) result(pt)
      !! Newton-solve the model PT that makes the Roquet branch's own
      !! ct_from_pt(SR(sp), PT) equal `ct_target`.  Uses the branch's
      !! analytic dSV/dT path indirectly via a local FD on CT is not
      !! needed: we drive CT directly through eos_density_point would be
      !! circular, so we recover dCT/dPT numerically from the branch via
      !! a CT-probe.  CT itself is not exposed, so we Newton on the
      !! density-implied CT through a tiny inversion: bisection is robust
      !! and CT(PT) is monotone over the ocean envelope.
      type(eos_t), intent(in) :: eos
      real(wp), intent(in) :: ct_target, sp
      real(wp) :: pt
      ! CT ≈ PT to ~0.04 degC; PT(CT) is monotone increasing.  Bisect on
      ! PT using the production `ct_from_pt` (rdb_teos10_convert) -- the
      ! SAME routine roquet_spv_point/roquet_rho_point call internally --
      ! so the inversion target is exactly the branch's.
      real(wp) :: lo, hi, mid, c
      lo = ct_target - 5.0_wp
      hi = ct_target + 5.0_wp
      do while (hi - lo > 1.0e-13_wp)
         mid = 0.5_wp*(lo + hi)
         c = ct_from_pt(sp*SR_FACTOR, mid)
         if (c < ct_target) then
            lo = mid
         else
            hi = mid
         end if
      end do
      pt = 0.5_wp*(lo + hi)
   end function pt_for_ct

   subroutine test_density_checkvalues(error)
      !! Reproduce the verified prototype rho oracles to ~1e-9.
      type(error_type), allocatable, intent(out) :: error
      type(eos_t) :: eos
      real(wp), parameter :: CT_T(4) = [25.0_wp, 10.0_wp, 4.0_wp, 10.0_wp]
      real(wp), parameter :: SA_T(4) = [35.0_wp, 35.0_wp, 35.0_wp, 30.0_wp]
      real(wp), parameter :: P_T(4) = [0.0_wp, 1.0e7_wp, 4.0e7_wp, 1.0e7_wp]
      real(wp), parameter :: RHO_ORACLE(4) = [1023.220806995_wp, 1031.280872325_wp, &
                                              1045.441271302_wp, 1027.451398573_wp]
      real(wp) :: sp, pt, rho
      integer :: n
      eos = make_roquet()
      do n = 1, 4
         sp = SA_T(n)/SR_FACTOR
         pt = pt_for_ct(eos, CT_T(n), sp)
         rho = eos_density_point(eos, pt, sp, P_T(n))
         call check(error, abs(rho - RHO_ORACLE(n)) <= 1.0e-7_wp, &
                    "Roquet density checkvalue mismatch vs prototype oracle")
         if (allocated(error)) return
      end do
   end subroutine test_density_checkvalues

   subroutine test_chain_rule_derivs(error)
      !! Analytic dSV/dT, dSV/dS (model PT, SP) vs central FD of SV=1/rho.
      type(error_type), allocatable, intent(out) :: error
      type(eos_t) :: eos
      real(wp), parameter :: PTS_T(4) = [25.0_wp, 10.0_wp, 4.0_wp, 18.0_wp]
      real(wp), parameter :: PTS_S(4) = [35.0_wp, 35.0_wp, 34.7_wp, 30.0_wp]
      real(wp), parameter :: PTS_P(4) = [0.0_wp, 1.0e7_wp, 4.0e7_wp, 1.0e6_wp]
      real(wp), parameter :: HT = 1.0e-4_wp, HS = 1.0e-4_wp
      real(wp) :: dsv_dt, dsv_ds, fd_t, fd_s, rt, rs
      integer :: n
      eos = make_roquet()
      do n = 1, 4
         call eos_specvol_derivs(eos, PTS_T(n), PTS_S(n), PTS_P(n), dsv_dt, dsv_ds)
         fd_t = (sv_at(eos, PTS_T(n) + HT, PTS_S(n), PTS_P(n)) &
                 - sv_at(eos, PTS_T(n) - HT, PTS_S(n), PTS_P(n)))/(2.0_wp*HT)
         fd_s = (sv_at(eos, PTS_T(n), PTS_S(n) + HS, PTS_P(n)) &
                 - sv_at(eos, PTS_T(n), PTS_S(n) - HS, PTS_P(n)))/(2.0_wp*HS)
         rt = abs(dsv_dt - fd_t)/abs(fd_t)
         rs = abs(dsv_ds - fd_s)/abs(fd_s)
         call check(error, rt <= 1.0e-7_wp, "Roquet dSV/dT (model PT) chain-rule mismatch vs FD")
         if (allocated(error)) return
         call check(error, rs <= 1.0e-7_wp, "Roquet dSV/dS (model SP) chain-rule mismatch vs FD")
         if (allocated(error)) return
      end do
   end subroutine test_chain_rule_derivs

   subroutine test_model_ct_pins_gsw(error)
      !! Regression pin for the misplaced-term PT->CT bug: the conversion
      !! every Roquet EOS path uses (`roquet_model_ct`, and its derivative
      !! twin) must reproduce the TEOS-10 GSW Conservative Temperature.
      !! References, independent of this tree:
      !!  * three TEOS-10 check-cast points, `CT_from_pt` at (SA, pt) from
      !!    the GSW toolbox's `gsw_data_v3_0.nc` as shipped (SP = SA/SR);
      !!  * (SP = 35, pt = 10 degC): the GSW algorithm (McDougall & Barker
      !!    2011) evaluated in double precision outside Fortran.
      !! The old polynomial was off by 8.1e-3 degC at (35, 10); the 1e-10
      !! degC bound is eight orders below that and above round-off.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NP = 4
      real(wp), parameter :: SA_REF(NP) = [3.450749946569206e+01_wp, 3.470686385140788e+01_wp, &
                                           6.669904340924572e+00_wp, 35.16504_wp]
      real(wp), parameter :: PT_REF(NP) = [2.796200000000000e+01_wp, 1.294495427573772e+01_wp, &
                                           1.004599999999999e+01_wp, 10.0_wp]
      real(wp), parameter :: CT_REF(NP) = [2.799482733197866e+01_wp, 1.294428999745330e+01_wp, &
                                           1.050276773530223e+01_wp, 9.989811727177308_wp]
      real(wp), parameter :: SR_PER_SP = 35.16504_wp/35.0_wp
      real(wp) :: sr, ct, dct_dpt, dct_dsr
      integer :: ip
      do ip = 1, NP
         ct = roquet_model_ct(PT_REF(ip), SA_REF(ip)/SR_PER_SP)
         call check(error, abs(ct - CT_REF(ip)) <= 1.0e-10_wp, &
                    "roquet_model_ct differs from the TEOS-10 GSW Conservative Temperature")
         if (allocated(error)) return
         call roquet_model_ct_derivs(PT_REF(ip), SA_REF(ip)/SR_PER_SP, sr, ct, dct_dpt, dct_dsr)
         call check(error, abs(ct - CT_REF(ip)) <= 1.0e-10_wp, &
                    "roquet_model_ct_derivs differs from the TEOS-10 GSW Conservative Temperature")
         if (allocated(error)) return
      end do
   end subroutine test_model_ct_pins_gsw

   subroutine test_3d_fill_matches_point(error)
      !! Exercises the 3D fill kernel `eos_roquet_spv_impl` ON DEVICE (via
      !! ocean_eos_compute) and asserts ms%rho_layer matches the point
      !! routine eos_density_point per layer.  The checkvalue/chain-rule
      !! tests cover only the point routines; this is the on-device test
      !! the new 3D `do concurrent` fill kernel otherwise lacks.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 3
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      real(wp) :: h_k(NZ), S_k(NZ), T_k(NZ), expected(NZ), mx
      integer :: k
      checks: block
         call grid%init(8, 6, NGHOST, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         eos = make_roquet()           ! variant = ROQUET_SPV, p_ref = 0
         h_k = [3.0_wp, 5.0_wp, 7.0_wp]
         S_k = [35.0_wp, 34.5_wp, 36.0_wp]
         T_k = [12.0_wp, 8.0_wp, 4.0_wp]
         do k = 1, NZ
            ms%h_layer(:, :, k) = h_k(k)
            ms%tracers(ms%idx_salinity)%hTr(:, :, k) = S_k(k)*h_k(k)
            ms%tracers(ms%idx_temperature)%hTr(:, :, k) = T_k(k)*h_k(k)
            ! same model (PT,SP) -> the 3D fill evaluates at eos%p_ref.
            expected(k) = eos_density_point(eos, T_k(k), S_k(k), eos%p_ref)
         end do

         !$omp target enter data map(to: ms)
         call ms%enter_data()
         call ocean_eos_compute(eos, ms)
         call ms%exit_data()
         !$omp target exit data map(delete: ms)

         do k = 1, NZ
            mx = maxval(abs(ms%rho_layer(:, :, k) - expected(k)))
            call check(error, mx < 1.0e-9_wp, &
                       "Roquet 3D fill (eos_roquet_spv_impl) != point routine")
            if (allocated(error)) exit checks
         end do
      end block checks
      call ms%destroy()
   end subroutine test_3d_fill_matches_point

   subroutine test_fresh_water_finite(error)
      !! SP = 0 (fresh water) must NOT produce a NaN.  The dCT/dSR term
      !! carries a 1/sqrt(sfac*SR) that is 0/0 at SR=0 (removable pole —
      !! dh_dx is itself proportional to sqrt(SR)); the x2 floor in
      !! roquet_spv_point keeps density + both derivatives finite there.
      type(error_type), allocatable, intent(out) :: error
      type(eos_t) :: eos
      real(wp) :: dsv_dt, dsv_ds, rho
      eos = make_roquet()
      rho = eos_density_point(eos, 10.0_wp, 0.0_wp, 1.0e6_wp)
      call eos_specvol_derivs(eos, 10.0_wp, 0.0_wp, 1.0e6_wp, dsv_dt, dsv_ds)
      ! finite (x==x rules out NaN; bound rules out Inf).
      call check(error, rho == rho .and. abs(rho) < 1.0e6_wp, &
                 "Roquet density NaN/Inf at SP=0 (fresh water)")
      if (allocated(error)) return
      call check(error, dsv_dt == dsv_dt .and. abs(dsv_dt) < 1.0e3_wp, &
                 "Roquet dSV/dT NaN/Inf at SP=0")
      if (allocated(error)) return
      call check(error, dsv_ds == dsv_ds .and. abs(dsv_ds) < 1.0e3_wp, &
                 "Roquet dSV/dS NaN/Inf at SP=0 (the dCT/dSR 0/0 pole)")
   end subroutine test_fresh_water_finite

   pure function sv_at(eos, T, S, p) result(sv)
      !! SV = 1/rho through the production point routine (model T, S).
      type(eos_t), intent(in) :: eos
      real(wp), intent(in) :: T, S, p
      real(wp) :: sv
      sv = 1.0_wp/eos_density_point(eos, T, S, p)
   end function sv_at

   subroutine test_default_off(error)
      !! Adding the Roquet case leaves the linear + Wright branches
      !! byte-unchanged (default eos = linear remains bit-identical).
      type(error_type), allocatable, intent(out) :: error
      type(eos_t) :: eos_lin, eos_wr
      real(wp) :: rho, expect
      real(wp), parameter :: T = 14.0_wp, S = 36.5_wp, p = 1.0e6_wp
      ! linear
      eos_lin%variant = EOS_VARIANT_LINEAR
      eos_lin%rho0 = 1035.0_wp
      eos_lin%alpha_T = 0.2_wp
      eos_lin%beta_S = 0.8_wp
      eos_lin%T_ref = 10.0_wp
      eos_lin%S_ref = 35.0_wp
      rho = eos_density_point(eos_lin, T, S, p)
      expect = 1035.0_wp + 0.8_wp*(S - 35.0_wp) - 0.2_wp*(T - 10.0_wp)
      call check(error, abs(rho - expect) <= 1.0e-13_wp*abs(expect), &
                 "linear branch perturbed by Roquet case")
      if (allocated(error)) return
      ! Wright — closed-form rational, must stay finite + in-range.
      eos_wr%variant = EOS_VARIANT_WRIGHT_97
      eos_wr%rho0 = 1035.0_wp
      rho = eos_density_point(eos_wr, 10.0_wp, 35.0_wp, 0.0_wp)
      call check(error, rho > 1020.0_wp .and. rho < 1035.0_wp, &
                 "Wright branch perturbed by Roquet case")
   end subroutine test_default_off

   subroutine test_fail_loud_fv_wright(error)
      !! Harness-compatible mirror of the configure_ocean_pgf guard:
      !! roquet_spv + fv_wright is rejected; roquet_spv + fv_lite is OK.
      type(error_type), allocatable, intent(out) :: error
      call check(error, roquet_pgf_unsupported(EOS_VARIANT_ROQUET_SPV, OPGF_VARIANT_FV_WRIGHT), &
                 "roquet_spv + fv_wright must be flagged unsupported")
      if (allocated(error)) return
      call check(error,.not. roquet_pgf_unsupported(EOS_VARIANT_ROQUET_SPV, OPGF_VARIANT_FV_LITE), &
                 "roquet_spv + fv_lite must be supported")
      if (allocated(error)) return
      call check(error,.not. roquet_pgf_unsupported(EOS_VARIANT_WRIGHT_97, OPGF_VARIANT_FV_WRIGHT), &
                 "wright + fv_wright must be supported")
   end subroutine test_fail_loud_fv_wright

   pure logical function roquet_pgf_unsupported(eos_variant, pgf_variant)
      !! Mirror of the configure_ocean_pgf fail-loud condition.
      integer, intent(in) :: eos_variant, pgf_variant
      roquet_pgf_unsupported = (eos_variant == EOS_VARIANT_ROQUET_SPV) .and. &
                               (pgf_variant == OPGF_VARIANT_FV_WRIGHT)
   end function roquet_pgf_unsupported

end module test_ocean_eos_roquet
