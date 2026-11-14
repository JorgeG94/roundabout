!! Tests of `rdb_ice_ridge_itd` (participation and redistribution preparation, ice
!! strength) and the `ice_compute_strength` driver: analytical tests A12 (the
!! participation fractions sum to one, the redistribution of each ridging category
!! conserves area-weighted volume) and A13 (Rothrock strength of the plan's reference
!! column, 5.887e4 N/m, and the Hibler-strength window), plus the Icepack oracle tables.
!!
!! Expected values come from closed forms written here, from the Icepack oracle (literals
!! in `test_cice_ridge_itd_golden.inc`, upstream compiled -O0 without FMA), or from brute
!! force; never from the port.  The oracle tables are pure arithmetic (and `sqrt`) except
!! where the participation function or the Hibler strength evaluate `exp`: those go
!! through the reproducible `rdb_sm_exp` (<= 1.2 ulp from libm), so they carry a
!! measured tolerance.  The arithmetic-only tables compare exactly on CPU builds (FMA off)
!! and to 1e-12 on the offload build, where FMA is on.
module test_cice_ridge_itd
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use, intrinsic :: iso_fortran_env, only: int32, int64
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use rdb_constants, only: wp
   use rdb_ice_dims, only: ICE_NCAT_MAX
   use rdb_ice_errors, only: ICE_ERR_NONFINITE, ICE_ERR_RIDGE_DEGENERATE
   use rdb_ice_params, only: ice_params_t, ice_params_default, ice_ridge_params_t, ice_constants_t, ice_constants, &
                             ice_ridge_params_valid, &
                             ICE_CONSTANTS_ICEPACK, ICE_CONSTANTS_CESM
   use rdb_ice_ridge_itd, only: ice_asum_cell, ridge_itd_cell, ice_strength_cell
   use rdb_ice_ridge, only: ice_compute_strength
   use lcg_deterministic, only: lcg_next, lcg_unit
   use test_cice_util, only: ice_golden_check
   implicit none
   private

   public :: collect_cice_ridge_itd_tests

   real(wp), parameter :: PUNY = 1.0e-11_wp   !! `puny` of both constant sets
   real(wp), parameter :: QNAN = transfer(int(z'7FF8000000000000', int64), 1.0_wp)
   integer, parameter :: NCAT = 5
   integer, parameter :: NCOL = 40
   real(wp), parameter :: HIN_MAX(0:NCAT) = [0.0_wp, 0.6445072168194258_wp, 1.3914334975763036_wp, &
                                             2.4701793819598885_wp, 4.567287918850491_wp, 999.9_wp]
      !! Category bounds (kcatbound = 0, ncat = 5; the D5.1 oracle table).
#ifdef RDB_GPU_OFFLOAD
   ! FMA is on for host code too on the offload build: scale-aware tolerances, relative plus an
   ! absolute floor of the O(1) participation scale.  Measured worst on the offload build (V100,
   ! nvfortran): 1.9e-15 relative on the ncat 1/7 oracle tables.
   real(wp), parameter :: TOL_ARITH = 1.0e-12_wp
   real(wp), parameter :: TOL_DEV = 1.0e-12_wp
   real(wp), parameter :: TOL_EXP = 1.0e-12_wp
   real(wp), parameter :: ATOL_EXP = 1.0e-12_wp
   real(wp), parameter :: ATOL_ARITH = 1.0e-12_wp
#else
   real(wp), parameter :: TOL_ARITH = 0.0_wp       ! pure arithmetic + sqrt: bitwise on CPU
   real(wp), parameter :: TOL_DEV = 1.0e-13_wp
   real(wp), parameter :: TOL_EXP = 1.0e-14_wp
      !! Tables through `exp`: `rdb_sm_exp` vs the oracle's libm.  Measured worst relative difference
      !! 8.7e-16 (apartic entries), 3.2e-16 (strength), 2.5e-16 (Hibler) on gfortran 15; tolerance ~10x.
   real(wp), parameter :: ATOL_EXP = 1.0e-15_wp
   real(wp), parameter :: ATOL_ARITH = 0.0_wp
#endif

   ! Oracle tables: GI inputs, GA aksum + apartic(0:5), GK krdg, GH hrmin/hrmax/hrexp, GS strength.
   ! Option-pair index k = 1 + 2*krdg_partic + krdg_redist; GS(5, :) is kstrength = 0.
   real(wp), save :: GI(13, NCOL), GA(7, 4, NCOL), GK(5, 4, NCOL), GH(15, 4, NCOL), GS(5, NCOL)
   ! ncat = 1 and ncat = 7 oracle tables: inputs aice0, aicen(n), vicen(n), aice, vice; aksum + apartic(0:n); strength
   real(wp), save :: GI1(5, NCOL), GA1(3, 4, NCOL), GS1(5, NCOL)
   real(wp), save :: GI7(17, NCOL), GA7(9, 4, NCOL), GS7(5, NCOL)

contains

   subroutine collect_cice_ridge_itd_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("partic_closed_forms", test_partic_closed_forms), &
                  new_unittest("partic_sums_to_one_random", test_partic_sums), &
                  new_unittest("redistribution_telescopes", test_redistribution_telescopes), &
                  new_unittest("rothrock_known_answer", test_rothrock_known), &
                  new_unittest("hibler_closed_form_and_window", test_hibler_window), &
                  new_unittest("oracle_ridge_tables", test_oracle_ridge), &
                  new_unittest("oracle_strength_tables", test_oracle_strength), &
                  new_unittest("oracle_ncat_1_and_7", test_oracle_ncat), &
                  new_unittest("option_validation", test_valid), &
                  new_unittest("empty_and_single_category_cells", test_empty_single), &
                  new_unittest("nonfinite_input_flagged", test_nonfinite), &
                  new_unittest("compute_strength_driver", test_driver), &
                  new_unittest("strength_driver_on_device", test_device) &
                  ]
   end subroutine collect_cice_ridge_itd_tests

   subroutine fill_golden()
      !! Fill the oracle tables (one statement per column).
#include "test_cice_ridge_itd_golden.inc"
   end subroutine fill_golden

   ! ------------------------------------------------------------------
   ! helpers
   ! ------------------------------------------------------------------

   subroutine random_cell(st, aice0, aicen, vicen)
      !! Pseudo-random consistent cell: ~30 % empty categories, mean thickness 0.05 .. 3.5 m,
      !! total ice area 0.1 .. 1.2 (so `aice0 = 0` and `asum > 1` occur), `aice0 = max(1-aice, 0)`.
      integer(int64), intent(inout) :: st
      real(wp), intent(out) :: aice0, aicen(NCAT), vicen(NCAT)
      integer :: n
      real(wp) :: u1, u2, u3, tot, tgt

      tot = 0.0_wp
      do n = 1, NCAT
         st = lcg_next(st); u1 = lcg_unit(st)
         st = lcg_next(st); u2 = lcg_unit(st)
         st = lcg_next(st); u3 = lcg_unit(st)
         aicen(n) = 0.01_wp + u2
         vicen(n) = aicen(n)*(0.05_wp + 3.45_wp*u3*u3)
         if (u1 < 0.3_wp) then
            aicen(n) = 0.0_wp
            vicen(n) = 0.0_wp
         end if
         tot = tot + aicen(n)
      end do
      st = lcg_next(st); tgt = 0.1_wp + 1.1_wp*lcg_unit(st)
      if (tot > 0.0_wp) then
         aicen = aicen*(tgt/tot)
         vicen = vicen*(tgt/tot)
      end if
      aice0 = max(1.0_wp - sum(aicen), 0.0_wp)
   end subroutine random_cell

   function cst() result(c)
      type(ice_constants_t) :: c

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
   end function cst

   function params_for(kstrength, kpartic, kredist) result(rp)
      integer, intent(in) :: kstrength, kpartic, kredist
      type(ice_ridge_params_t) :: rp
      type(ice_params_t) :: p

      p = ice_params_default(ICE_CONSTANTS_ICEPACK, NCAT, 4, 1)
      rp = p%ridge
      rp%kstrength = kstrength
      rp%krdg_partic = kpartic
      rp%krdg_redist = kredist
   end function params_for

   function strength_of(rp, aice0, aicen, vicen, ierr) result(s)
      type(ice_ridge_params_t), intent(in) :: rp
      real(wp), intent(in) :: aice0, aicen(NCAT), vicen(NCAT)
      integer(int32), intent(inout) :: ierr
      real(wp) :: s

      call ice_strength_cell(NCAT, rp, cst(), sum(aicen), sum(vicen), aice0, aicen, vicen, s, ierr)
   end function strength_of

   subroutine ridge_of(rp, aice0, aicen, vicen, aksum, apartic, hrmin, hrmax, hrexp, krdg)
      type(ice_ridge_params_t), intent(in) :: rp
      real(wp), intent(in) :: aice0, aicen(NCAT), vicen(NCAT)
      real(wp), intent(out) :: aksum, apartic(0:NCAT), hrmin(NCAT), hrmax(NCAT), hrexp(NCAT), krdg(NCAT)

      call ridge_itd_cell(NCAT, rp, PUNY, aice0, aicen, vicen, aksum, apartic, hrmin, hrmax, hrexp, krdg)
   end subroutine ridge_of

   function rel_diff(x, y) result(d)
      real(wp), intent(in) :: x, y
      real(wp) :: d

      d = abs(x - y)/max(abs(y), tiny(1.0_wp))
   end function rel_diff

   ! ------------------------------------------------------------------
   ! A12: participation
   ! ------------------------------------------------------------------

   subroutine test_partic_closed_forms(error)
      !! Hand-evaluated participation functions.  Thorndike: G* = 0.15 truncation branch and
      !! the open-water-saturated case; Lipscomb: the telescoped exponential differences.
      type(error_type), allocatable, intent(out) :: error
      type(ice_ridge_params_t) :: rp
      real(wp) :: ak, ap(0:NCAT), hrmin(NCAT), hrmax(NCAT), hrexp(NCAT), krdg(NCAT)
      real(wp) :: G(0:2), e(0:2), norm, aice0, aicen(NCAT), vicen(NCAT), ex(0:2), eps

      ! open water 0.05, two categories 0.45 and 0.5: G = (0.05, 0.5, 1)
      aice0 = 0.05_wp
      aicen = [0.45_wp, 0.5_wp, 0.0_wp, 0.0_wp, 0.0_wp]
      vicen = [0.45_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
      G = [0.05_wp, 0.5_wp, 1.0_wp]

      rp = params_for(1, 0, 1)
      call ridge_of(rp, aice0, aicen, vicen, ak, ap, hrmin, hrmax, hrexp, krdg)
      ! apartic(0) = (1/G*)(0.05)(2 - 0.05/G*); apartic(1) = (1/G*)(G* - 0.05)(2 - (G* + 0.05)/G*)
      ex = [(1.0_wp/0.15_wp)*0.05_wp*(2.0_wp - 0.05_wp/0.15_wp), &
            (1.0_wp/0.15_wp)*(0.15_wp - 0.05_wp)*(2.0_wp - (0.15_wp + 0.05_wp)/0.15_wp), 0.0_wp]
      eps = 1.0e-14_wp
      call ice_golden_check(error, "Thorndike apartic(0:2)", ap(0:2), ex, eps, eps)
      if (allocated(error)) return
      call ice_golden_check(error, "Thorndike apartic(3:5) zero", ap(3:5), [0.0_wp, 0.0_wp, 0.0_wp], 0.0_wp, 0.0_wp)
      if (allocated(error)) return
      call check(error, ap(0) + ap(1), 1.0_wp, "Thorndike sum", thr=1.0e-14_wp)
      if (allocated(error)) return

      ! Thorndike, open water above G*: all participation in open water
      aice0 = 0.3_wp
      aicen = [0.4_wp, 0.3_wp, 0.0_wp, 0.0_wp, 0.0_wp]
      vicen = [0.4_wp, 0.6_wp, 0.0_wp, 0.0_wp, 0.0_wp]
      call ridge_of(rp, aice0, aicen, vicen, ak, ap, hrmin, hrmax, hrexp, krdg)
      call ice_golden_check(error, "Thorndike open water saturates", ap, [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp], &
                            1.0e-14_wp, 1.0e-14_wp)
      if (allocated(error)) return

      ! Lipscomb: apartic(n) = (exp(-G(n-1)/a*) - exp(-G(n)/a*))/(1 - exp(-1/a*)), a* = 0.05
      aice0 = 0.05_wp
      aicen = [0.45_wp, 0.5_wp, 0.0_wp, 0.0_wp, 0.0_wp]
      vicen = [0.45_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
      rp = params_for(1, 1, 1)
      call ridge_of(rp, aice0, aicen, vicen, ak, ap, hrmin, hrmax, hrexp, krdg)
      norm = 1.0_wp - exp(-1.0_wp/0.05_wp)
      e = [exp(-G(0)/0.05_wp), exp(-G(1)/0.05_wp), exp(-G(2)/0.05_wp)]
      ex = [(1.0_wp - e(0))/norm, (e(0) - e(1))/norm, (e(1) - e(2))/norm]
      ! open water is category 0: apartic(0) = e(-1) - e(0) with G(-1) = 0
      call ice_golden_check(error, "Lipscomb apartic(0:2)", ap(0:2), ex, 1.0e-13_wp, 1.0e-15_wp)
      if (allocated(error)) return
      call check(error, ap(0) + ap(1) + ap(2), 1.0_wp, "Lipscomb sum", thr=1.0e-14_wp)
   end subroutine test_partic_closed_forms

   subroutine test_partic_sums(error)
      !! A12: 4000 random cells x both options: the participation fractions sum to 1
      !! (open water included), are zero where `aicen <= puny` and non-negative, and `aksum`
      !! lies in (0, 1].
      type(error_type), allocatable, intent(out) :: error
      type(ice_ridge_params_t) :: rp
      integer(int64) :: st
      integer :: c, n, kp, kr
      real(wp) :: aice0, aicen(NCAT), vicen(NCAT), ak, ap(0:NCAT), hrmin(NCAT), hrmax(NCAT), hrexp(NCAT), krdg(NCAT)
      real(wp) :: worst, s, tiny_area(NCAT)

      worst = 0.0_wp
      st = lcg_next(20261113_int64)
      do c = 1, 1000
         call random_cell(st, aice0, aicen, vicen)
         ! one category below the puny threshold: it must not participate
         tiny_area = aicen
         if (aicen(2) > 0.0_wp) then
            tiny_area(2) = 0.5_wp*PUNY
         end if
         do kp = 0, 1
            do kr = 0, 1
               rp = params_for(1, kp, kr)
               call ridge_of(rp, aice0, tiny_area, vicen, ak, ap, hrmin, hrmax, hrexp, krdg)
               s = sum(ap)
               if (.not. (ieee_is_finite(s) .and. ieee_is_finite(ak))) then
                  call test_failed(error, "non-finite participation sum or aksum")
                  return
               end if
               if (.not. (abs(s - 1.0_wp) <= worst)) worst = abs(s - 1.0_wp)
               if (.not. (abs(s - 1.0_wp) <= 1.0e-14_wp)) then
                  call test_failed(error, "participation fractions do not sum to 1")
                  return
               end if
               if (.not. (ak > 0.0_wp .and. ak <= 1.0_wp + 1.0e-14_wp)) then
                  call test_failed(error, "aksum outside (0, 1]")
                  return
               end if
               do n = 0, NCAT
                  if (.not. (ap(n) >= -1.0e-15_wp)) then
                     call test_failed(error, "negative participation")
                     return
                  end if
               end do
               do n = 1, NCAT
                  if (tiny_area(n) <= PUNY .and. ap(n) /= 0.0_wp) then
                     call test_failed(error, "empty category participates")
                     return
                  end if
               end do
            end do
         end do
      end do
      call check(error, worst < 1.0e-14_wp, .true., "sum(apartic) - 1 worst below 1e-14")
   end subroutine test_partic_sums

   ! ------------------------------------------------------------------
   ! A12: redistribution
   ! ------------------------------------------------------------------

   subroutine test_redistribution_telescopes(error)
      !! A12: for every ridging category the area and volume fractions each ridge puts into
      !! the `ncat` thickness categories sum to 1, and the mean ridge thickness they imply
      !! equals `krdg*hi` of the port (Hibler 1980: uniform on [hrmin, hrmax];
      !! Lipscomb et al. 2007: exponential beyond hrmin).  The fractions are written here
      !! from the plan's formulas, independent of the port.
      type(error_type), allocatable, intent(out) :: error
      type(ice_ridge_params_t) :: rp
      integer(int64) :: st
      integer :: c, n, m, kr
      real(wp) :: aice0, aicen(NCAT), vicen(NCAT), ak, ap(0:NCAT), hrmin(NCAT), hrmax(NCAT), hrexp(NCAT), krdg(NCAT)
      real(wp) :: hi, hL, hR, fa(NCAT), fv(NCAT), mean, lam, eL, eR, hm(NCAT), worst
      integer :: nchecked

      worst = 0.0_wp
      nchecked = 0
      st = lcg_next(777_int64)
      do c = 1, 600
         call random_cell(st, aice0, aicen, vicen)
         do kr = 0, 1
            rp = params_for(1, 1, kr)
            call ridge_of(rp, aice0, aicen, vicen, ak, ap, hrmin, hrmax, hrexp, krdg)
            hm = HIN_MAX(1:NCAT)
            hm(NCAT) = 1.0e8_wp
            do n = 1, NCAT
               if (.not. (aicen(n) > PUNY)) cycle
               hi = vicen(n)/aicen(n)
               if (kr == 1) hi = max(hi, PUNY)
               fa = 0.0_wp
               fv = 0.0_wp
               do m = 1, NCAT
                  if (kr == 0) then
                     hL = max(hrmin(n), HIN_MAX(m - 1))
                     hR = min(hrmax(n), hm(m))
                     if (hrmin(n) >= hm(m) .or. hrmax(n) <= HIN_MAX(m - 1)) cycle
                     fa(m) = (hR - hL)/(hrmax(n) - hrmin(n))
                     fv(m) = (hR*hR - hL*hL)/(hrmax(n)*hrmax(n) - hrmin(n)*hrmin(n))
                  else
                     lam = hrexp(n)
                     if (hrmin(n) >= hm(m)) cycle
                     hL = max(hrmin(n), HIN_MAX(m - 1))
                     eL = exp(-(hL - hrmin(n))/lam)
                     if (m < NCAT) then
                        eR = exp(-(HIN_MAX(m) - hrmin(n))/lam)
                        fa(m) = eL - eR
                        fv(m) = ((hL + lam)*eL - (HIN_MAX(m) + lam)*eR)/(hrmin(n) + lam)
                     else
                        fa(m) = eL
                        fv(m) = (hL + lam)*eL/(hrmin(n) + lam)
                     end if
                  end if
               end do
               nchecked = nchecked + 1
               if (.not. (abs(sum(fa) - 1.0_wp) <= 1.0e-13_wp .and. abs(sum(fv) - 1.0_wp) <= 1.0e-13_wp)) then
                  call test_failed(error, "ridge area / volume fractions do not sum to 1")
                  return
               end if
               ! mean ridge thickness implied by the distribution = krdg*hi
               if (kr == 0) then
                  mean = 0.5_wp*(hrmax(n)*hrmax(n) - hrmin(n)*hrmin(n))*sum(fv)/((hrmax(n) - hrmin(n))*sum(fa))
               else
                  mean = (hrmin(n) + lam)*sum(fv)/sum(fa)
               end if
               if (.not. (rel_diff(mean, krdg(n)*hi) <= worst)) worst = rel_diff(mean, krdg(n)*hi)
               if (.not. (rel_diff(mean, krdg(n)*hi) <= 1.0e-12_wp)) then
                  call test_failed(error, "mean ridge thickness differs from krdg*hi")
                  return
               end if
            end do
         end do
      end do
      call check(error, nchecked > 1000, .true., "enough ridging categories exercised")
   end subroutine test_redistribution_telescopes

   ! ------------------------------------------------------------------
   ! A13: strength
   ! ------------------------------------------------------------------

   subroutine test_rothrock_known(error)
      !! A13: the reference column `aice0 = 0`, one category `aicen = 1`, `hi = 1`:
      !! `krdg_redist = 1`: `hrmin = 2`, `lam = 3`, `krdg = 5`, `h2rdg = 4 + 12 + 18 = 34`,
      !! `dh2rdg = -1 + 34/5 = 5.8`, `aksum = 1 - 1/5 = 0.8`, `P = Cf Cp 5.8/0.8` = 5.887e4 N/m.
      !! `krdg_redist = 0`: `hrmax = 2 sqrt(25) = 10`, `krdg = 6`, `h2rdg = (10^3 - 2^3)/(3*8)`.
      !! Both constant sets (identical for these constants).
      type(error_type), allocatable, intent(out) :: error
      type(ice_ridge_params_t) :: rp
      type(ice_params_t) :: p
      integer(int32) :: ierr
      real(wp) :: s, ex, cp, aicen(NCAT), vicen(NCAT)
      integer :: kset, kp

      aicen = [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
      vicen = [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
      do kset = ICE_CONSTANTS_ICEPACK, ICE_CONSTANTS_CESM
         p = ice_params_default(kset, NCAT, 4, 1)
         cp = 0.5_wp*p%c%gravit*(p%c%rhow - p%c%rhoi)*p%c%rhoi/p%c%rhow    ! (g/2)(rhow-rhoi) rhoi/rhow
         do kp = 0, 1
            ! redist = 1 (OM3 default)
            rp = p%ridge
            rp%kstrength = 1
            rp%krdg_partic = kp
            rp%krdg_redist = 1
            ierr = 0_int32
            s = strength_of(rp, 0.0_wp, aicen, vicen, ierr)
            ex = 17.0_wp*cp*((-1.0_wp + 34.0_wp/5.0_wp)/0.8_wp)
            if (.not. ieee_is_finite(s)) then
               call test_failed(error, "Rothrock strength not finite")
               return
            end if
            call check(error, rel_diff(s, ex) < 1.0e-13_wp, .true., "Rothrock strength, redist=1, closed form")
            if (allocated(error)) return
            call check(error, rel_diff(s, 5.887e4_wp) < 1.0e-3_wp, .true., "Rothrock strength is 5.887e4 N/m")
            if (allocated(error)) return
            call check(error, ierr, 0_int32, "no error bits")
            if (allocated(error)) return
            ! redist = 0
            rp%krdg_redist = 0
            s = strength_of(rp, 0.0_wp, aicen, vicen, ierr)
            ex = 17.0_wp*cp*((-1.0_wp + ((1000.0_wp - 8.0_wp)/(3.0_wp*(10.0_wp - 2.0_wp)))/6.0_wp)/(1.0_wp - 1.0_wp/6.0_wp))
            call check(error, rel_diff(s, ex) < 1.0e-13_wp, .true., "Rothrock strength, redist=0, closed form")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_rothrock_known

   subroutine test_hibler_window(error)
      !! A13: `kstrength = 0` against `Pstar vice exp(-Cstar(1-aice))` (closed form), then
      !! the Rothrock / Hibler ratio for one-category columns of concentration {0.9, 1} and
      !! thickness {0.5, 1, 2, 3} lies in [0.5, 5] (measured 1.36 .. 4.00); Rothrock strength increases with
      !! concentration at fixed ITD shape and the thin-ice scaling is `h^1.5` (Rothrock) versus
      !! `h` (Hibler).
      type(error_type), allocatable, intent(out) :: error
      type(ice_ridge_params_t) :: rh, rr
      integer(int32) :: ierr
      integer :: ia, ih, i
      real(wp), parameter :: AS(2) = [0.9_wp, 1.0_wp], HS(4) = [0.5_wp, 1.0_wp, 2.0_wp, 3.0_wp]
      real(wp) :: a, h, sh, sr, aicen(NCAT), vicen(NCAT), worst, ratio, prev, slope_r, slope_h

      rh = params_for(0, 1, 1)
      rr = params_for(1, 1, 1)
      ierr = 0_int32
      worst = 0.0_wp
      do ia = 1, 20
         do ih = 1, 20
            a = 0.05_wp*real(ia, wp)
            h = 0.15_wp*real(ih, wp)
            aicen = [a, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
            vicen = [a*h, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
            sh = strength_of(rh, 1.0_wp - a, aicen, vicen, ierr)
            if (.not. (rel_diff(sh, 2.75e4_wp*(a*h)*exp(-20.0_wp*(1.0_wp - a))) <= worst)) &
               worst = rel_diff(sh, 2.75e4_wp*(a*h)*exp(-20.0_wp*(1.0_wp - a)))
         end do
      end do
      call check(error, worst < 1.0e-14_wp, .true., "Hibler strength matches Pstar vice exp(-Cstar(1-aice))")
      if (allocated(error)) return
      call check(error, ierr, 0_int32, "no error bits")
      if (allocated(error)) return
      ! the plan's value at a = 1, h = 1
      aicen = [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
      vicen = aicen
      sh = strength_of(rh, 0.0_wp, aicen, vicen, ierr)
      call check(error, rel_diff(sh, 2.75e4_wp) < 1.0e-14_wp, .true., "Hibler strength is 27 500 N/m at a=1, h=1")
      if (allocated(error)) return

      ! window
      do ia = 1, 2
         do ih = 1, 4
            a = AS(ia)
            h = HS(ih)
            aicen = [a, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
            vicen = [a*h, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
            sh = strength_of(rh, 1.0_wp - a, aicen, vicen, ierr)
            sr = strength_of(rr, 1.0_wp - a, aicen, vicen, ierr)
            ratio = sr/sh
            if (.not. (ratio >= 0.5_wp .and. ratio <= 5.0_wp)) then
               call test_failed(error, "Rothrock / Hibler ratio outside [0.5, 5]")
               return
            end if
         end do
      end do
      ! monotone in concentration at fixed shape (h = 1 m)
      prev = -1.0_wp
      do i = 1, 20
         a = 0.05_wp*real(i, wp)
         aicen = [a, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
         vicen = [a, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
         sr = strength_of(rr, 1.0_wp - a, aicen, vicen, ierr)
         if (.not. (sr > prev)) then
            call test_failed(error, "Rothrock strength not increasing with concentration")
            return
         end if
         prev = sr
      end do
      ! thin-ice scaling between h = 0.1 and 0.4 m at a = 1
      aicen = [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
      vicen = [0.1_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
      sr = strength_of(rr, 0.0_wp, aicen, vicen, ierr)
      sh = strength_of(rh, 0.0_wp, aicen, vicen, ierr)
      vicen = [0.4_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
      slope_r = log(strength_of(rr, 0.0_wp, aicen, vicen, ierr)/sr)/log(4.0_wp)
      slope_h = log(strength_of(rh, 0.0_wp, aicen, vicen, ierr)/sh)/log(4.0_wp)
      call check(error, abs(slope_h - 1.0_wp) < 1.0e-12_wp, .true., "Hibler thin-ice slope is 1")
      if (allocated(error)) return
      call check(error, slope_r > 1.4_wp .and. slope_r < 1.65_wp, .true., "Rothrock thin-ice slope near 1.54")
   end subroutine test_hibler_window

   ! ------------------------------------------------------------------
   ! Icepack oracle tables
   ! ------------------------------------------------------------------

   subroutine test_oracle_ridge(error)
      !! 40 oracle columns x 4 option pairs: `aksum`, `apartic(0:5)`, `krdg`, `hrmin`, `hrmax`,
      !! `hrexp`.  `krdg_partic = 0` rows are arithmetic only (exact on CPU builds);
      !! `krdg_partic = 1` rows carry the `exp` tolerance.
      type(error_type), allocatable, intent(out) :: error
      type(ice_ridge_params_t) :: rp
      integer :: c, kp, kr, k
      real(wp) :: ak, ap(0:NCAT), hrmin(NCAT), hrmax(NCAT), hrexp(NCAT), krdg(NCAT), tol
      character(len=48) :: lab

      call fill_golden()
      do c = 1, NCOL
         do kp = 0, 1
            do kr = 0, 1
               k = 1 + 2*kp + kr
               rp = params_for(1, kp, kr)
               call ridge_of(rp, GI(1, c), GI(2:6, c), GI(7:11, c), ak, ap, hrmin, hrmax, hrexp, krdg)
               tol = merge(TOL_EXP, TOL_ARITH, kp == 1)
               write (lab, '(a,i0,a,i0,a,i0)') "oracle ridge col ", c - 1, " partic ", kp, " redist ", kr
         call ice_golden_check(error, trim(lab)//" aksum+apartic", [ak, ap], GA(:, k, c), tol, merge(ATOL_EXP, ATOL_ARITH, kp == 1))
               if (allocated(error)) return
               ! krdg/hrmin/hrmax/hrexp depend on the participation only through aksum, not themselves
               call ice_golden_check(error, trim(lab)//" krdg", krdg, GK(:, k, c), TOL_ARITH, ATOL_ARITH)
               if (allocated(error)) return
            call ice_golden_check(error, trim(lab)//" hrmin/hrmax/hrexp", [hrmin, hrmax, hrexp], GH(:, k, c), TOL_ARITH, ATOL_ARITH)
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_oracle_ridge

   subroutine test_oracle_strength(error)
      !! 40 oracle columns: Rothrock strength for the four option pairs and Hibler strength.
      type(error_type), allocatable, intent(out) :: error
      type(ice_ridge_params_t) :: rp
      type(ice_params_t) :: p
      integer(int32) :: ierr
      integer :: c, kp, kr, k, kset
      real(wp) :: s
      character(len=48) :: lab

      call fill_golden()
      ierr = 0_int32
      do kset = ICE_CONSTANTS_ICEPACK, ICE_CONSTANTS_CESM
         ! both constant sets give the same tables (nothing ridging reads differs)
         p = ice_params_default(kset, NCAT, 4, 1)
         do c = 1, NCOL
            do kp = 0, 1
               do kr = 0, 1
                  k = 1 + 2*kp + kr
                  rp = p%ridge
                  rp%kstrength = 1
                  rp%krdg_partic = kp
                  rp%krdg_redist = kr
                  s = strength_of(rp, GI(1, c), GI(2:6, c), GI(7:11, c), ierr)
                  write (lab, '(a,i0,a,i0,a,i0)') "oracle strength col ", c - 1, " partic ", kp, " redist ", kr
                  call ice_golden_check(error, trim(lab), [s], [GS(k, c)], merge(TOL_EXP, TOL_ARITH, kp == 1), &
                                        merge(ATOL_EXP, ATOL_ARITH, kp == 1))
                  if (allocated(error)) return
               end do
            end do
            rp = p%ridge
            rp%kstrength = 0
            s = strength_of(rp, GI(1, c), GI(2:6, c), GI(7:11, c), ierr)
            write (lab, '(a,i0)') "oracle Hibler strength col ", c - 1
            call ice_golden_check(error, trim(lab), [s], [GS(5, c)], TOL_EXP, ATOL_EXP)
            if (allocated(error)) return
         end do
      end do
      call check(error, ierr, 0_int32, "no error bits on consistent cells")
   end subroutine test_oracle_strength

   ! ------------------------------------------------------------------
   ! edge cases
   ! ------------------------------------------------------------------

   subroutine oracle_n(error, n, gi, ga, gs)
      !! Oracle columns of one category count: `aksum`, `apartic(0:n)` and both strengths, every
      !! option pair.  `gi` rows: aice0, aicen(n), vicen(n), ...; `ga` rows: aksum, apartic(0:n).
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: n
      real(wp), intent(in) :: gi(:, :), ga(:, :, :), gs(:, :)
      type(ice_ridge_params_t) :: rp
      integer :: c, kp, kr, k
      real(wp) :: ak, ap(0:n), hrmin(n), hrmax(n), hrexp(n), krdg(n), s, aice, vice
      integer(int32) :: ierr
      character(len=48) :: lab

      ierr = 0_int32
      do c = 1, NCOL
         aice = sum(gi(2:n + 1, c))
         vice = sum(gi(n + 2:2*n + 1, c))
         do kp = 0, 1
            do kr = 0, 1
               k = 1 + 2*kp + kr
               rp = params_for(1, kp, kr)
               call ridge_itd_cell(n, rp, PUNY, gi(1, c), gi(2:n + 1, c), gi(n + 2:2*n + 1, c), ak, ap, hrmin, hrmax, hrexp, krdg)
               write (lab, '(a,i0,a,i0,a,i0,a,i0)') "ncat ", n, " col ", c - 1, " partic ", kp, " redist ", kr
               call ice_golden_check(error, trim(lab)//" aksum+apartic", [ak, ap], ga(:, k, c), &
                                     merge(TOL_EXP, TOL_ARITH, kp == 1), merge(ATOL_EXP, ATOL_ARITH, kp == 1))
               if (allocated(error)) return
               call ice_strength_cell(n, rp, cst(), aice, vice, gi(1, c), gi(2:n + 1, c), gi(n + 2:2*n + 1, c), s, ierr)
               call ice_golden_check(error, trim(lab)//" strength", [s], [gs(k, c)], merge(TOL_EXP, TOL_ARITH, kp == 1), &
                                     merge(ATOL_EXP, ATOL_ARITH, kp == 1))
               if (allocated(error)) return
            end do
         end do
         rp = params_for(0, 1, 1)
         call ice_strength_cell(n, rp, cst(), aice, vice, gi(1, c), gi(2:n + 1, c), gi(n + 2:2*n + 1, c), s, ierr)
         call ice_golden_check(error, "ncat Hibler strength", [s], [gs(5, c)], TOL_EXP, ATOL_EXP)
         if (allocated(error)) return
      end do
      call check(error, ierr, 0_int32, "no error bits on consistent cells")
   end subroutine oracle_n

   subroutine test_oracle_ncat(error)
      !! Oracle tables at ncat = 1 and ncat = 7 (the category count enters the cumulative sums, the
      !! last-category handling and the local array extents).
      type(error_type), allocatable, intent(out) :: error

      call fill_golden()
      call oracle_n(error, 1, GI1, GA1, GS1)
      if (allocated(error)) return
      call oracle_n(error, 7, GI7, GA7, GS7)
   end subroutine test_oracle_ncat

   subroutine test_valid(error)
      !! Option validation: 0/1 accepted, anything else refused (the drivers stop on it; the
      !! config PR must refuse it at configure).
      type(error_type), allocatable, intent(out) :: error
      type(ice_ridge_params_t) :: rp

      rp = params_for(1, 1, 1)
      call check(error, ice_ridge_params_valid(rp), .true., "defaults valid")
      if (allocated(error)) return
      rp%kstrength = 2
      call check(error, ice_ridge_params_valid(rp), .false., "kstrength = 2 refused")
      if (allocated(error)) return
      rp = params_for(1, -1, 1)
      call check(error, ice_ridge_params_valid(rp), .false., "krdg_partic = -1 refused")
      if (allocated(error)) return
      rp = params_for(0, 0, 3)
      call check(error, ice_ridge_params_valid(rp), .false., "krdg_redist = 3 refused")
   end subroutine test_valid

   subroutine test_empty_single(error)
      !! Empty open-ocean cell (`aice0 = 1`): strength 0 and no error for both strength options;
      !! inconsistent empty cell (`aice0 = 0`, no ice): Rothrock gives 0 with
      !! `ICE_ERR_RIDGE_DEGENERATE`, Hibler 0 and no error; single-category cell: only that
      !! category participates; categories at or below `puny` do not.
      type(error_type), allocatable, intent(out) :: error
      type(ice_ridge_params_t) :: rp
      integer(int32) :: ierr
      real(wp) :: s, z(NCAT), ak, ap(0:NCAT), hrmin(NCAT), hrmax(NCAT), hrexp(NCAT), krdg(NCAT)
      real(wp) :: aicen(NCAT), vicen(NCAT)
      integer :: ks

      z = 0.0_wp
      do ks = 0, 1
         rp = params_for(ks, 1, 1)
         ierr = 0_int32
         s = strength_of(rp, 1.0_wp, z, z, ierr)
         call check(error, s, 0.0_wp, "empty open-ocean cell has zero strength")
         if (allocated(error)) return
         call check(error, ierr, 0_int32, "empty open-ocean cell is clean")
         if (allocated(error)) return
      end do
      rp = params_for(1, 1, 1)
      ierr = 0_int32
      s = strength_of(rp, 0.0_wp, z, z, ierr)
      call check(error, s, 0.0_wp, "inconsistent empty cell: strength 0")
      if (allocated(error)) return
      call check(error, iand(ierr, ICE_ERR_RIDGE_DEGENERATE) /= 0_int32, .true., "inconsistent empty cell: degenerate bit")
      if (allocated(error)) return
      rp = params_for(0, 1, 1)
      ierr = 0_int32
      s = strength_of(rp, 0.0_wp, z, z, ierr)
      call check(error, s == 0.0_wp .and. ierr == 0_int32, .true., "Hibler on an empty full-cover cell is 0 and clean")
      if (allocated(error)) return

      ! single category, the third one: only it (and open water) participates
      aicen = [0.0_wp, 0.0_wp, 0.6_wp, 0.0_wp, 0.0_wp]
      vicen = [0.0_wp, 0.0_wp, 1.2_wp, 0.0_wp, 0.0_wp]
      rp = params_for(1, 1, 1)
      call ridge_of(rp, 0.4_wp, aicen, vicen, ak, ap, hrmin, hrmax, hrexp, krdg)
      call check(error, ap(1) == 0.0_wp .and. ap(2) == 0.0_wp .and. ap(4) == 0.0_wp .and. ap(5) == 0.0_wp, .true., &
                 "only the occupied category participates")
      if (allocated(error)) return
      call check(error, abs(sum(ap) - 1.0_wp) < 1.0e-14_wp, .true., "single category: participation sums to 1")
      if (allocated(error)) return
      call check(error, hrmin(3) == 3.0_wp .and. krdg(1) == 1.0_wp, .true., "hi = 2: hrmin = min(4, 3) = 3 (rafting cap)")
      if (allocated(error)) return
      ! area budget: ice plus open water may exceed 1 after transport
      call check(error, ice_asum_cell(NCAT, [0.6_wp, 0.5_wp, 0.0_wp, 0.0_wp, 0.0_wp], 0.1_wp), 1.2_wp, "asum = aice0 + sum aicen", &
                 thr=1.0e-15_wp)
   end subroutine test_empty_single

   subroutine test_nonfinite(error)
      !! NaN in any input: non-finite strength and the non-finite bit, for both strength options;
      !! the bit is ORed, never cleared.
      type(error_type), allocatable, intent(out) :: error
      type(ice_ridge_params_t) :: rp
      integer(int32) :: ierr
      real(wp) :: s, aicen(NCAT), vicen(NCAT), aice0
      integer :: ks, which

      do ks = 0, 1
         rp = params_for(ks, 1, 1)
         do which = 1, 3
            aicen = [0.3_wp, 0.2_wp, 0.1_wp, 0.0_wp, 0.0_wp]
            vicen = [0.3_wp, 0.4_wp, 0.3_wp, 0.0_wp, 0.0_wp]
            aice0 = 0.4_wp
            select case (which)
            case (1)
               aicen(2) = QNAN
            case (2)
               vicen(3) = QNAN
            case default
               aice0 = QNAN
            end select
            ierr = 4_int32
            s = strength_of(rp, aice0, aicen, vicen, ierr)
            if (ieee_is_finite(s)) then
               call test_failed(error, "NaN input must not give a finite strength")
               return
            end if
            call check(error, iand(ierr, ICE_ERR_NONFINITE) /= 0_int32 .and. iand(ierr, 4_int32) /= 0_int32, .true., &
                       "NaN input sets the non-finite bit and keeps the old bits")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_nonfinite

   ! ------------------------------------------------------------------
   ! driver
   ! ------------------------------------------------------------------

   subroutine test_driver(error)
      !! `ice_compute_strength` on a 6x5 grid: dry cells and empty open-ocean cells are exactly 0;
      !! wet ice cells match the cell kernel (`TOL_DEV`: FMA on the offload build); only the index window is written; error bits are
      !! ORed into `ice_err`; a NaN-area cell is not skipped and is flagged.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nx = 6, ny = 5
      type(ice_ridge_params_t) :: rp
      type(ice_constants_t) :: cc
      integer(int64) :: st
      real(wp) :: aice(nx, ny), vice(nx, ny), aice0(nx, ny), aicen(nx, ny, NCAT), vicen(nx, ny, NCAT), strength(nx, ny)
      real(wp) :: a0, an(NCAT), vn(NCAT), s
      logical :: wet(nx, ny)
      integer(int32) :: ice_err(nx, ny), ierr
      integer :: i, j

      st = lcg_next(99_int64)
      do j = 1, ny
         do i = 1, nx
            call random_cell(st, a0, an, vn)
            aicen(i, j, :) = an
            vicen(i, j, :) = vn
            aice(i, j) = sum(an)
            vice(i, j) = sum(vn)
            aice0(i, j) = a0
         end do
      end do
      wet = .true.
      wet(2, 2) = .false.
      aicen(3, 3, :) = 0.0_wp   ! empty open-ocean cell
      vicen(3, 3, :) = 0.0_wp
      aice(3, 3) = 0.0_wp
      vice(3, 3) = 0.0_wp
      aice0(3, 3) = 1.0_wp
      aice(4, 4) = QNAN         ! corrupted cell
      strength = -7.0_wp
      ice_err = 0_int32
      ice_err(5, 2) = 8_int32
      rp = params_for(1, 1, 1)
      cc = cst()

      !$acc enter data copyin(rp, cc, wet, aice, vice, aice0, aicen, vicen, strength, ice_err)
      call ice_compute_strength(nx, ny, NCAT, 2, nx - 1, 2, ny - 1, rp, cc, wet, aice, vice, aice0, aicen, vicen, &
                                strength, ice_err)
      !$acc update self(strength, ice_err)
      !$acc exit data delete(rp, cc, wet, aice, vice, aice0, aicen, vicen, strength, ice_err)

      call check(error, strength(1, 1) == -7.0_wp .and. strength(nx, ny) == -7.0_wp .and. strength(1, ny) == -7.0_wp, .true., &
                 "cells outside the window are untouched")
      if (allocated(error)) return
      call check(error, strength(2, 2) == 0.0_wp, .true., "dry cell is exactly 0")
      if (allocated(error)) return
      call check(error, strength(3, 3) == 0.0_wp .and. ice_err(3, 3) == 0_int32, .true., "empty open-ocean cell is exactly 0")
      if (allocated(error)) return
      call check(error, iand(ice_err(4, 4), ICE_ERR_NONFINITE) /= 0_int32 .and. .not. ieee_is_finite(strength(4, 4)), .true., &
                 "NaN-area cell is flagged, not skipped")
      if (allocated(error)) return
      call check(error, ice_err(5, 2), 8_int32, "pre-existing error bits survive")
      if (allocated(error)) return
      do j = 2, ny - 1
         do i = 2, nx - 1
            if (.not. wet(i, j) .or. (i == 3 .and. j == 3) .or. (i == 4 .and. j == 4)) cycle
            ierr = 0_int32
            s = strength_of(rp, aice0(i, j), aicen(i, j, :), vicen(i, j, :), ierr)
            if (.not. (rel_diff(strength(i, j), s) <= TOL_DEV)) then
               call test_failed(error, "driver differs from the cell kernel")
               return
            end if
         end do
      end do
   end subroutine test_driver

   subroutine test_device(error)
      !! The driver on the offload build (a host `do concurrent` otherwise): every array and the
      !! params bundle are mapped (`mem:separate`), results compared with the host evaluation of
      !! the cell kernel (`TOL_DEV`; FMA is on there).  Columns: random, a dry cell, an empty cell
      !! and a NaN cell; both strength options.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nx = 40, ny = 30
      type(ice_ridge_params_t) :: rp
      type(ice_constants_t) :: cc
      integer(int64) :: st
      real(wp), allocatable :: aice(:, :), vice(:, :), aice0(:, :), aicen(:, :, :), vicen(:, :, :), strength(:, :)
      logical, allocatable :: wet(:, :)
      integer(int32), allocatable :: ice_err(:, :)
      real(wp) :: a0, an(NCAT), vn(NCAT), s, worst
      integer(int32) :: ierr
      integer :: i, j, ks

      allocate (aice(nx, ny), vice(nx, ny), aice0(nx, ny), aicen(nx, ny, NCAT), vicen(nx, ny, NCAT), strength(nx, ny))
      allocate (wet(nx, ny), ice_err(nx, ny))
      st = lcg_next(31337_int64)
      do j = 1, ny
         do i = 1, nx
            call random_cell(st, a0, an, vn)
            aicen(i, j, :) = an
            vicen(i, j, :) = vn
            aice(i, j) = sum(an)
            vice(i, j) = sum(vn)
            aice0(i, j) = a0
         end do
      end do
      wet = .true.
      wet(7, 7) = .false.
      aicen(8, 8, :) = 0.0_wp
      vicen(8, 8, :) = 0.0_wp
      aice(8, 8) = 0.0_wp
      vice(8, 8) = 0.0_wp
      aice0(8, 8) = 1.0_wp
      aice(9, 9) = QNAN

      do ks = 0, 1
         rp = params_for(ks, 1, 1)
         cc = cst()
         strength = -1.0_wp
         ice_err = 0_int32
         !$acc enter data copyin(rp, cc, wet, aice, vice, aice0, aicen, vicen, strength, ice_err)
         call ice_compute_strength(nx, ny, NCAT, 1, nx, 1, ny, rp, cc, wet, aice, vice, aice0, aicen, vicen, strength, ice_err)
         !$acc update self(strength, ice_err)
         !$acc exit data delete(rp, cc, wet, aice, vice, aice0, aicen, vicen, strength, ice_err)

         worst = 0.0_wp
         do j = 1, ny
            do i = 1, nx
               if (i == 9 .and. j == 9) then
                  if (ieee_is_finite(strength(i, j)) .or. iand(ice_err(i, j), ICE_ERR_NONFINITE) == 0_int32) then
                     call test_failed(error, "device: NaN cell flagged")
                     return
                  end if
                  cycle
               end if
               if (.not. wet(i, j) .or. (i == 8 .and. j == 8)) then
                  if (strength(i, j) /= 0.0_wp) then
                     call test_failed(error, "device: dry / empty cell is exactly 0")
                     return
                  end if
                  cycle
               end if
               ierr = 0_int32
               s = strength_of(rp, aice0(i, j), aicen(i, j, :), vicen(i, j, :), ierr)
               if (.not. (rel_diff(strength(i, j), s) <= worst)) worst = rel_diff(strength(i, j), s)
               if (ice_err(i, j) /= 0_int32) then
                  call test_failed(error, "device: clean cell has error bits")
                  return
               end if
            end do
         end do
         if (.not. (worst <= TOL_DEV)) then
            call test_failed(error, "device strength differs from the host evaluation")
            return
         end if
      end do
   end subroutine test_device

end module test_cice_ridge_itd
