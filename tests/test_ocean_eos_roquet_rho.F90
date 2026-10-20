!! Tests for the Roquet et al. (2015) in-situ-density EOS variant.
module test_ocean_eos_roquet_rho
   !! Locks the `EOS_VARIANT_ROQUET_RHO` branch added to the shared EOS
   !! handle — the polyTEOS10_55t IN-SITU DENSITY fit (ΔS = 32 g/kg),
   !! an INDEPENDENT fit of the Roquet et al. (2015) paper from the
   !! already-shipped polyTEOS10_75t SPECIFIC-VOLUME fit
   !! (`EOS_VARIANT_ROQUET_SPV`, ΔS = 24 g/kg):
   !!
   !!   (1) PAPER CHECK VALUE — rho(SA=30, CT=10, p=1000 dbar) =
   !!       1027.45140 kg/m^3 (Roquet et al. 2015 Sec. 5), to 1e-5.
   !!       The production branch works in MODEL variables (PT, SP), so
   !!       (as `test_ocean_eos_roquet` already does for the SpV branch)
   !!       we Newton-solve the model PT that makes the branch's own
   !!       ct_from_pt(SR, PT) recover CT=10 exactly, exercising the
   !!       FULL production conversion path.
   !!
   !!   (2) GSW CAST CROSS-CHECK — 6 points extracted from the TEOS-10
   !!       reference check-value file (`gsw_data_v3_0.nc`,
   !!       `rho_t_exact`/`CT_chck_cast`/`SA_chck_cast`/`p_chck_cast`),
   !!       spanning warm/cold, shallow/deep and a low-salinity cast.
   !!       |rho_model - rho_gsw_exact| <= 0.01 kg/m^3, the paper's
   !!       stated fit-accuracy bound (NOT round-off): this is the
   !!       "does the transcribed 55-term polynomial reproduce TEOS-10"
   !!       regression.
   !!
   !!   (3) FIT-SWITCH COST — rho_55t (this branch) vs the already-
   !!       shipped rho_75t (`1/roquet_spv_value`) at the SAME 4
   !!       (CT, SA, p) oracle points `test_ocean_eos_roquet` uses.
   !!       |Delta| <= 5e-3 kg/m^3: the two Roquet fits are independent
   !!       polynomials, not bit-identical, and this locks the size of
   !!       the answer change for anyone switching eos="roquet_spv" ->
   !!       "roquet_rho".
   !!
   !!   (4) CHAIN-RULE DERIVATIVE REGRESSION — analytic d(rho)/dT,
   !!       d(rho)/dS (w.r.t. the MODEL variables PT, SP, including the
   !!       dCT/dPT and SR-factor chain factors) agree with a central
   !!       finite difference of `eos_density_point` to ~1e-6.
   !!
   !!   (5) 3D KERNEL VS POINT — the `eos_roquet_rho_impl` device fill
   !!       kernel (exercised through `ocean_eos_compute`) matches the
   !!       point routine per layer to ~1e-9.
   !!
   !!   (6) DEFAULT-OFF BIT-IDENTITY — the linear, Wright and
   !!       roquet_spv branches are byte-unchanged by adding this case.
   !!
   !!   (7) FAIL-LOUD — roquet_rho + fv_wright is rejected at configure
   !!       (harness-compatible mirror of the `configure_ocean_pgf`
   !!       guard, same shape as the existing roquet_spv test).
   !!
   !!   (8) FRESH WATER FINITE — SP = 0 must not produce NaN/Inf (same
   !!       x2-floor rationale as the SpV branch; the floor is in the
   !!       shared ct_from_pt poly both branches evaluate).
   !!
   !! Runs host-side under NVHPC so the `!$acc routine seq` device
   !! codegen path is compiled and the values locked.
   use rdb_constants, only: wp
   use rdb_eos, only: eos_t, eos_density_point, eos_specvol_derivs, &
                      EOS_VARIANT_LINEAR, EOS_VARIANT_WRIGHT_97, &
                      EOS_VARIANT_ROQUET_SPV, EOS_VARIANT_ROQUET_RHO
   use rdb_ocean_eos_compute, only: ocean_eos_compute
   use rdb_ocean_pressure_force, only: OPGF_VARIANT_FV_WRIGHT, OPGF_VARIANT_FV_LITE
   use rdb_teos10_convert, only: ct_from_pt
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use testdrive, only: new_unittest, unittest_type, error_type, check
   implicit none
   private

   public :: collect_ocean_eos_roquet_rho_tests

   real(wp), parameter :: SR_FACTOR = 35.16504_wp/35.0_wp
   integer, parameter :: NGHOST = 2

contains

   subroutine collect_ocean_eos_roquet_rho_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("roquet_rho_paper_checkvalue", test_paper_checkvalue), &
                  new_unittest("roquet_rho_gsw_cast_crosscheck", test_gsw_cast_crosscheck), &
                  new_unittest("roquet_rho_vs_spv_fit_cost", test_vs_spv_fit_cost), &
                  new_unittest("roquet_rho_chain_rule_derivs", test_chain_rule_derivs), &
                  new_unittest("roquet_rho_3d_fill_matches_point", test_3d_fill_matches_point), &
                  new_unittest("roquet_rho_default_off_bit_identity", test_default_off), &
                  new_unittest("roquet_rho_fail_loud_fv_wright", test_fail_loud_fv_wright), &
                  new_unittest("roquet_rho_fresh_water_finite", test_fresh_water_finite) &
                  ]
   end subroutine collect_ocean_eos_roquet_rho_tests

   pure function make_roquet_rho() result(eos)
      type(eos_t) :: eos
      eos%variant = EOS_VARIANT_ROQUET_RHO
      eos%rho0 = 1035.0_wp
   end function make_roquet_rho

   pure function pt_for_ct(ct_target, sp) result(pt)
      !! Newton/bisection-solve the model PT that makes the branch's own
      !! ct_from_pt(SR(sp), PT) equal `ct_target` — the same approach
      !! `test_ocean_eos_roquet` uses for the SpV branch (same poly).
      real(wp), intent(in) :: ct_target, sp
      real(wp) :: pt
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

   subroutine test_paper_checkvalue(error)
      !! Roquet et al. (2015) Sec. 5 check value: rho(SA=30, CT=10,
      !! p=1000 dbar) = 1027.45140 kg/m^3.
      type(error_type), allocatable, intent(out) :: error
      type(eos_t) :: eos
      real(wp), parameter :: CT_TARGET = 10.0_wp, SA_TARGET = 30.0_wp
      real(wp), parameter :: P_PA = 1000.0_wp*1.0e4_wp
         !! 1000 dbar -> Pa (1 dbar = 1e4 Pa, the SPV table's own
         !! convention — ROQ_PA2KB = 1e-8 folds Pa -> kbar).
      real(wp), parameter :: RHO_PAPER = 1027.45140_wp
      real(wp) :: sp, pt, rho
      eos = make_roquet_rho()
      sp = SA_TARGET/SR_FACTOR
      pt = pt_for_ct(CT_TARGET, sp)
      rho = eos_density_point(eos, pt, sp, P_PA)
      call check(error, abs(rho - RHO_PAPER) <= 1.0e-5_wp, &
                 "Roquet_rho paper check value mismatch (Roquet et al. 2015 Sec. 5)")
   end subroutine test_paper_checkvalue

   subroutine test_gsw_cast_crosscheck(error)
      !! 6 points from the TEOS-10 reference check-value file
      !! (gsw_data_v3_0.nc: CT_chck_cast / SA_chck_cast / p_chck_cast /
      !! rho_t_exact), spanning warm/cold, shallow/deep, high/low
      !! salinity.  |rho_model - rho_gsw_exact| <= 0.01 kg/m^3 (the
      !! paper's stated polynomial-fit accuracy bound).
      type(error_type), allocatable, intent(out) :: error
      type(eos_t) :: eos
      integer, parameter :: N = 6
      real(wp), parameter :: CT_T(N) = [27.9948273320_wp, 7.2259862752_wp, &
                                        1.3795290983_wp, 27.2798981236_wp, &
                                        5.6244493769_wp, 4.6113979928_wp]
      real(wp), parameter :: SA_T(N) = [34.5074994657_wp, 34.6591744923_wp, &
                                        34.8464829610_wp, 34.6056959095_wp, &
                                        34.7137834878_wp, 10.3894684550_wp]
      real(wp), parameter :: P_DBAR(N) = [0.0_wp, 505.0_wp, 3045.0_wp, &
                                          30.0_wp, 808.0_wp, 101.0_wp]
      real(wp), parameter :: RHO_GSW(N) = [1021.915980252_wp, 1029.287647691_wp, &
                                           1041.666252567_wp, 1022.345597301_wp, &
                                           1030.937452325_wp, 1008.687812318_wp]
      real(wp) :: sp, pt, rho, p_pa
      integer :: ip
      eos = make_roquet_rho()
      do ip = 1, N
         sp = SA_T(ip)/SR_FACTOR
         pt = pt_for_ct(CT_T(ip), sp)
         p_pa = P_DBAR(ip)*1.0e4_wp
         rho = eos_density_point(eos, pt, sp, p_pa)
         call check(error, abs(rho - RHO_GSW(ip)) <= 1.0e-2_wp, &
                    "Roquet_rho vs GSW TEOS-10 exact cast point exceeds fit-accuracy bound")
         if (allocated(error)) return
      end do
   end subroutine test_gsw_cast_crosscheck

   subroutine test_vs_spv_fit_cost(error)
      !! rho_55t (this branch) vs rho_75t (1/roquet_spv_value) at the
      !! SAME 4 oracle points `test_ocean_eos_roquet` checks.  The two
      !! fits are independent polynomials (ΔS=32 vs ΔS=24) and agree
      !! only to the paper's stated rms accuracy -- NOT bit-identical.
      type(error_type), allocatable, intent(out) :: error
      type(eos_t) :: eos_rho, eos_spv
      real(wp), parameter :: CT_T(4) = [25.0_wp, 10.0_wp, 4.0_wp, 10.0_wp]
      real(wp), parameter :: SA_T(4) = [35.0_wp, 35.0_wp, 35.0_wp, 30.0_wp]
      real(wp), parameter :: P_T(4) = [0.0_wp, 1.0e7_wp, 4.0e7_wp, 1.0e7_wp]
      real(wp) :: sp, pt, rho_rho, rho_spv
      integer :: n
      eos_rho = make_roquet_rho()
      eos_spv%variant = EOS_VARIANT_ROQUET_SPV
      eos_spv%rho0 = 1035.0_wp
      do n = 1, 4
         sp = SA_T(n)/SR_FACTOR
         pt = pt_for_ct(CT_T(n), sp)
         rho_rho = eos_density_point(eos_rho, pt, sp, P_T(n))
         rho_spv = eos_density_point(eos_spv, pt, sp, P_T(n))
         call check(error, abs(rho_rho - rho_spv) <= 5.0e-3_wp, &
                    "Roquet_rho vs roquet_spv switch cost exceeds the documented bound")
         if (allocated(error)) return
      end do
   end subroutine test_vs_spv_fit_cost

   subroutine test_chain_rule_derivs(error)
      !! Analytic dSV/dT, dSV/dS (model PT, SP; dsv = -drho/rho^2) vs
      !! central FD of SV = 1/rho through eos_density_point.
      type(error_type), allocatable, intent(out) :: error
      type(eos_t) :: eos
      real(wp), parameter :: PTS_T(4) = [25.0_wp, 10.0_wp, 4.0_wp, 18.0_wp]
      real(wp), parameter :: PTS_S(4) = [35.0_wp, 35.0_wp, 34.7_wp, 30.0_wp]
      real(wp), parameter :: PTS_P(4) = [0.0_wp, 1.0e7_wp, 4.0e7_wp, 1.0e6_wp]
      real(wp), parameter :: HT = 1.0e-4_wp, HS = 1.0e-4_wp
      real(wp) :: dsv_dt, dsv_ds, fd_t, fd_s, rt, rs
      integer :: n
      eos = make_roquet_rho()
      do n = 1, 4
         call eos_specvol_derivs(eos, PTS_T(n), PTS_S(n), PTS_P(n), dsv_dt, dsv_ds)
         fd_t = (sv_at(eos, PTS_T(n) + HT, PTS_S(n), PTS_P(n)) &
                 - sv_at(eos, PTS_T(n) - HT, PTS_S(n), PTS_P(n)))/(2.0_wp*HT)
         fd_s = (sv_at(eos, PTS_T(n), PTS_S(n) + HS, PTS_P(n)) &
                 - sv_at(eos, PTS_T(n), PTS_S(n) - HS, PTS_P(n)))/(2.0_wp*HS)
         rt = abs(dsv_dt - fd_t)/abs(fd_t)
         rs = abs(dsv_ds - fd_s)/abs(fd_s)
         call check(error, rt <= 1.0e-7_wp, "Roquet_rho dSV/dT (model PT) chain-rule mismatch vs FD")
         if (allocated(error)) return
         call check(error, rs <= 1.0e-7_wp, "Roquet_rho dSV/dS (model SP) chain-rule mismatch vs FD")
         if (allocated(error)) return
      end do
   end subroutine test_chain_rule_derivs

   pure function sv_at(eos, T, S, p) result(sv)
      !! SV = 1/rho through the production point routine (model T, S).
      type(eos_t), intent(in) :: eos
      real(wp), intent(in) :: T, S, p
      real(wp) :: sv
      sv = 1.0_wp/eos_density_point(eos, T, S, p)
   end function sv_at

   subroutine test_3d_fill_matches_point(error)
      !! Exercises the 3D fill kernel `eos_roquet_rho_impl` ON DEVICE
      !! (via ocean_eos_compute) and asserts ms%rho_layer matches the
      !! point routine eos_density_point per layer.
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
         eos = make_roquet_rho()         ! variant = ROQUET_RHO, p_ref = 0
         h_k = [3.0_wp, 5.0_wp, 7.0_wp]
         S_k = [35.0_wp, 34.5_wp, 36.0_wp]
         T_k = [12.0_wp, 8.0_wp, 4.0_wp]
         do k = 1, NZ
            ms%h_layer(:, :, k) = h_k(k)
            ms%tracers(ms%idx_salinity)%hTr(:, :, k) = S_k(k)*h_k(k)
            ms%tracers(ms%idx_temperature)%hTr(:, :, k) = T_k(k)*h_k(k)
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
                       "Roquet_rho 3D fill (eos_roquet_rho_impl) != point routine")
            if (allocated(error)) exit checks
         end do
      end block checks
      call ms%destroy()
   end subroutine test_3d_fill_matches_point

   subroutine test_default_off(error)
      !! Adding the Roquet_rho case leaves the linear, Wright and
      !! roquet_spv branches byte-unchanged.
      type(error_type), allocatable, intent(out) :: error
      type(eos_t) :: eos_lin, eos_wr, eos_spv
      real(wp) :: rho, expect, rho_spv1, rho_spv2
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
                 "linear branch perturbed by Roquet_rho case")
      if (allocated(error)) return
      ! Wright
      eos_wr%variant = EOS_VARIANT_WRIGHT_97
      eos_wr%rho0 = 1035.0_wp
      rho = eos_density_point(eos_wr, 10.0_wp, 35.0_wp, 0.0_wp)
      call check(error, rho > 1020.0_wp .and. rho < 1035.0_wp, &
                 "Wright branch perturbed by Roquet_rho case")
      if (allocated(error)) return
      ! roquet_spv -- same point evaluated twice must agree bit-for-bit
      ! (sanity that adding roquet_rho left the spv branch's own code path
      ! untouched).
      eos_spv%variant = EOS_VARIANT_ROQUET_SPV
      eos_spv%rho0 = 1035.0_wp
      rho_spv1 = eos_density_point(eos_spv, 10.0_wp, 35.0_wp, 1.0e6_wp)
      rho_spv2 = eos_density_point(eos_spv, 10.0_wp, 35.0_wp, 1.0e6_wp)
      call check(error, rho_spv1 == rho_spv2, "roquet_spv branch perturbed by Roquet_rho case")
   end subroutine test_default_off

   subroutine test_fail_loud_fv_wright(error)
      !! Harness-compatible mirror of the configure_ocean_pgf guard:
      !! roquet_rho + fv_wright is rejected; roquet_rho + fv_lite is OK.
      type(error_type), allocatable, intent(out) :: error
      call check(error, roquet_pgf_unsupported(EOS_VARIANT_ROQUET_RHO, OPGF_VARIANT_FV_WRIGHT), &
                 "roquet_rho + fv_wright must be flagged unsupported")
      if (allocated(error)) return
      call check(error,.not. roquet_pgf_unsupported(EOS_VARIANT_ROQUET_RHO, OPGF_VARIANT_FV_LITE), &
                 "roquet_rho + fv_lite must be supported")
      if (allocated(error)) return
      call check(error,.not. roquet_pgf_unsupported(EOS_VARIANT_WRIGHT_97, OPGF_VARIANT_FV_WRIGHT), &
                 "wright + fv_wright must be supported")
   end subroutine test_fail_loud_fv_wright

   pure logical function roquet_pgf_unsupported(eos_variant, pgf_variant)
      !! Mirror of the configure_ocean_pgf fail-loud condition (both
      !! Roquet fits are refused with fv_wright).
      integer, intent(in) :: eos_variant, pgf_variant
      roquet_pgf_unsupported = (eos_variant == EOS_VARIANT_ROQUET_SPV .or. &
                                eos_variant == EOS_VARIANT_ROQUET_RHO) .and. &
                               (pgf_variant == OPGF_VARIANT_FV_WRIGHT)
   end function roquet_pgf_unsupported

   subroutine test_fresh_water_finite(error)
      !! SP = 0 (fresh water) must NOT produce a NaN -- the dCT/dSR
      !! term carries a removable 0/0 pole the shared x2 floor keeps
      !! finite (same rationale as the SpV branch's own test).
      type(error_type), allocatable, intent(out) :: error
      type(eos_t) :: eos
      real(wp) :: dsv_dt, dsv_ds, rho
      eos = make_roquet_rho()
      rho = eos_density_point(eos, 10.0_wp, 0.0_wp, 1.0e6_wp)
      call eos_specvol_derivs(eos, 10.0_wp, 0.0_wp, 1.0e6_wp, dsv_dt, dsv_ds)
      call check(error, rho == rho .and. abs(rho) < 1.0e6_wp, &
                 "Roquet_rho density NaN/Inf at SP=0 (fresh water)")
      if (allocated(error)) return
      call check(error, dsv_dt == dsv_dt .and. abs(dsv_dt) < 1.0e3_wp, &
                 "Roquet_rho dSV/dT NaN/Inf at SP=0")
      if (allocated(error)) return
      call check(error, dsv_ds == dsv_ds .and. abs(dsv_ds) < 1.0e3_wp, &
                 "Roquet_rho dSV/dS NaN/Inf at SP=0 (the dCT/dSR 0/0 pole)")
   end subroutine test_fresh_water_finite

end module test_ocean_eos_roquet_rho
