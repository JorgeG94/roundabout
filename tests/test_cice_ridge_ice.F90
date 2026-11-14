!! Tests of the ridging step (`rdb_ice_ridge_shift`, `ice_step_ridge`): the plan's A12 cases 3-9 and 12-14.
!!
!! * `pure_convergence_known_answer` (3), `single_category_full_cover` (4, both redistribution forms and the
!!   shear-only fraction), `pure_divergence_opens_water` (5), `limiter_forces_second_pass` (6, with the
!!   non-convergence bit), `snow_flux_signs` (8), `no_ridging_identity` (9), `empty_cell_skip` (13): closed forms
!!   written HERE from the formulas of the papers (libm `exp`), never from the port;
!! * `exact_conservation_random_columns` (7): 100 000 random cells for nilyr = 1, 4, 7 over both participation
!!   and both redistribution forms: ice volume, snow volume plus the snow thrown into the ocean, ice and snow
!!   energy plus `fhocn dt`, salt, and total area = 1 after the step, worst residuals printed;
!! * `err_flags` (12): `ICE_ERR_NONFINITE`, `ICE_ERR_RIDGE_AICE0`, `ICE_ERR_RIDGE_DEGENERATE`,
!!   `ICE_ERR_RIDGE_NOCONV` (`ICE_ERR_RIDGE_OVERSHOOT` is a defensive guard the limiter makes unreachable, see
!!   the module note), a NaN is refused and never comes out finite;
!! * the Icepack ORACLE: the real upstream `ridge_ice`, 12 cells per constant set per participation x
!!   redistribution combination (literals in `test_cice_ridge_ice_golden.inc`; both constant sets, nilyr 7, 4, 1,
!!   multi-pass cells, `ndtd = 2`).  The arithmetic-only combination (Thorndike participation, Hibler
!!   redistribution) agrees bitwise on CPU builds; the others evaluate `exp` through the reproducible
!!   `rdb_sm_exp` and carry a measured tolerance; the offload build is scale-aware;
!! * `device_sweep_and_throughput` (14): `ice_step_ridge` over an INTERIOR range of a random state against a host
!!   reference (ridge cell + clean-up cell), decisions exact, the ring untouched, then the throughput on a
!!   random state and on a spatially coherent state with ~85 % ice-free cells.
module test_cice_ridge_ice
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan, ieee_value, ieee_quiet_nan
   use, intrinsic :: iso_fortran_env, only: int32, int64, output_unit
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use rdb_constants, only: wp
   use rdb_ice_dims, only: ICE_NTRCR_MAX, ICE_NCAT_MAX
use rdb_ice_errors, only: ICE_ERR_FATAL_MASK, ICE_ERR_NONFINITE, ICE_ERR_RIDGE_NOCONV, ICE_ERR_RIDGE_AICE0, ICE_ERR_RIDGE_DEGENERATE
   use rdb_ice_params, only: ice_params_t, ice_params_default, ice_ntrcr, ice_trcr_get, ice_trcr_set, ice_qsno_init_from_T, &
                             ICE_CONSTANTS_ICEPACK, ICE_CONSTANTS_CESM
   use rdb_ice_ridge_itd, only: RIDGE_NITERMAX
   use rdb_ice_ridge_shift, only: ice_ridge_ice_cell
   use rdb_ice_ridge, only: ice_step_ridge
   use rdb_ice_itd_cleanup, only: ice_cleanup_itd_cell, ice_cleanup_itd_step
   use lcg_deterministic, only: lcg_next, lcg_unit
   implicit none
   private

   public :: collect_cice_ridge_ice_tests

   real(wp), parameter :: PUNY = 1.0e-11_wp
   integer, parameter :: NCAT = 5
   integer, parameter :: NR = 60
   real(wp), parameter :: HIN(0:NCAT) = [0.0_wp, 0.6445072168194258_wp, 1.3914334975763036_wp, &
                                         2.4701793819598885_wp, 4.567287918850491_wp, 999.9_wp]
#ifdef RDB_GPU_OFFLOAD
   logical, parameter :: GPU_BUILD = .true.
   real(wp), parameter :: TOL_GOLD = 1.0e-12_wp
   real(wp), parameter :: TOL_EXT = 2.0e-14_wp
   real(wp), parameter :: TOL_INT = 5.0e-14_wp
#else
   logical, parameter :: GPU_BUILD = .false.
   real(wp), parameter :: TOL_GOLD = 0.0_wp
   real(wp), parameter :: TOL_EXT = 2.0e-14_wp
   real(wp), parameter :: TOL_INT = 5.0e-14_wp
#endif
   real(wp), parameter :: TOL_EXP = 4.0e-15_wp
      !! Oracle rows that evaluate `exp` (`rdb_sm_exp` vs libm, a few ulp per exponential, carried through the
      !! cancellations of the redistribution and the next pass).
   real(wp), parameter :: CLOSE_TOL = 1.0e-13_wp
   real(wp), parameter :: SNOW0_RHOS = 330.0_wp

   type :: cell_t
      !! One packed cell (host test container).
      real(wp) :: a(NCAT), v(NCAT), s(NCAT)
      real(wp) :: tr(ICE_NTRCR_MAX, NCAT)
      real(wp) :: a0, fresh, fhocn, d1, d2, dv, op
      integer(int32) :: ierr
      integer :: niter
   end type cell_t

   integer, save :: GMETA(5, NR)
   real(wp), save :: GHD(5, NR), GII(100, NR), GIC(100, NR), GOI(110, NR), GOC(110, NR)

contains

   subroutine collect_cice_ridge_ice_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("pure_convergence_known_answer", test_pure_convergence), &
                  new_unittest("single_category_full_cover", test_single_category), &
                  new_unittest("shear_only_fraction", test_shear_only), &
                  new_unittest("pure_divergence_opens_water", test_divergence), &
                  new_unittest("limiter_forces_second_pass", test_limiter), &
                  new_unittest("snow_flux_signs", test_snow_signs), &
                  new_unittest("no_ridging_identity", test_no_ridging), &
                  new_unittest("empty_cell_skip", test_empty_skip), &
                  new_unittest("err_flags", test_err_flags), &
                  new_unittest("exact_conservation_random_columns_nilyr7", test_conserve_7), &
                  new_unittest("exact_conservation_random_columns_nilyr4", test_conserve_4), &
                  new_unittest("exact_conservation_random_columns_nilyr1", test_conserve_1), &
                  new_unittest("oracle_ridge_icepack", test_oracle_icepack), &
                  new_unittest("oracle_ridge_cesm", test_oracle_cesm), &
                  new_unittest("device_sweep_and_throughput", test_device) &
                  ]
   end subroutine collect_cice_ridge_ice_tests

   subroutine fill_golden()
      GMETA = 0; GHD = 0.0_wp; GII = 0.0_wp; GIC = 0.0_wp; GOI = 0.0_wp; GOC = 0.0_wp
#include "test_cice_ridge_ice_golden.inc"
   end subroutine fill_golden

   ! ------------------------------------------------------------------
   ! helpers
   ! ------------------------------------------------------------------

   subroutine blank(c)
      !! Zero every field (no default initialisers in the type: nvfortran does not apply them reliably to a host
      !! structure constructor).
      type(cell_t), intent(out) :: c

      c%a = 0.0_wp; c%v = 0.0_wp; c%s = 0.0_wp; c%tr = 0.0_wp
      c%a0 = 0.0_wp; c%fresh = 0.0_wp; c%fhocn = 0.0_wp; c%d1 = 0.0_wp; c%d2 = 0.0_wp; c%dv = 0.0_wp; c%op = 0.0_wp
      c%ierr = 0_int32; c%niter = 0
   end subroutine blank

   pure logical function same(a, b, tol, floor)
      !! Bitwise (`tol = 0`) or relative agreement against `max(|b|, floor)`; NaN matches NaN only.
      real(wp), intent(in) :: a, b, tol, floor

      if (ieee_is_nan(a) .or. ieee_is_nan(b)) then
         same = ieee_is_nan(a) .and. ieee_is_nan(b)
      else if (tol == 0.0_wp) then
         same = (a == b)
      else
         same = abs(a - b) <= tol*max(abs(b), floor)
      end if
   end function same

   subroutine params(kset, nilyr, kpartic, kredist, p)
      integer, intent(in) :: kset, nilyr, kpartic, kredist
      type(ice_params_t), intent(out) :: p

      p = ice_params_default(kset, NCAT, nilyr, 1)
      p%ridge%krdg_partic = kpartic
      p%ridge%krdg_redist = kredist
      p%ridge%mu_rdg = 3.0_wp
   end subroutine params

   subroutine base_tracers(nilyr, p, Tsn, tr)
      !! Packed tracers: Tsfc -5, qice(k) -(2e8 + 1e7 k), sice(k) 3 + 0.2 k, snow enthalpy of `Tsn` degC.
      integer, intent(in) :: nilyr
      type(ice_params_t), intent(in) :: p
      real(wp), intent(in) :: Tsn
      real(wp), intent(out) :: tr(ICE_NTRCR_MAX, NCAT)
      integer :: n, k

      tr = 0.0_wp
      do n = 1, NCAT
         tr(1, n) = -5.0_wp
         do k = 1, nilyr
            tr(1 + k, n) = -(2.0e8_wp + 1.0e7_wp*real(k, wp))
            tr(2 + nilyr + k, n) = 3.0_wp + 0.2_wp*real(k, wp)
         end do
         tr(2 + nilyr, n) = ice_qsno_init_from_T(Tsn, p%c)
      end do
   end subroutine base_tracers

   subroutine run_cell(p, nilyr, dt, ndtd, itermax, tfx, conv, shear, c)
      !! The cell kernel on a 1 x 1 grid, the state in and out through the packed layout; the fluxes and
      !! diagnostics of `c` are added to, the error word ORed.
      type(ice_params_t), intent(in) :: p
      integer, intent(in) :: nilyr, ndtd, itermax
      real(wp), intent(in) :: dt, tfx, conv, shear
      type(cell_t), intent(inout) :: c
      real(wp) :: a3(1, 1, NCAT), v3(1, 1, NCAT), s3(1, 1, NCAT), T3(1, 1, NCAT), q4(1, 1, NCAT, nilyr), si4(1, 1, NCAT, nilyr)
      real(wp) :: qs4(1, 1, NCAT, 1), a0(1, 1), fr(1, 1), fh(1, 1), d1(1, 1), d2(1, 1), dv(1, 1), op(1, 1)
      integer(int32) :: ie(1, 1)
      integer :: n, it, nit

      do n = 1, NCAT
         a3(1, 1, n) = c%a(n); v3(1, 1, n) = c%v(n); s3(1, 1, n) = c%s(n)
         do it = 1, ice_ntrcr(nilyr, 1)
            call ice_trcr_set(it, 1, 1, n, 1, 1, NCAT, nilyr, 1, c%tr(it, n), T3, q4, si4, qs4)
         end do
      end do
      a0 = c%a0; fr = c%fresh; fh = c%fhocn; d1 = c%d1; d2 = c%d2; dv = c%dv; op = c%op; ie = c%ierr
      call ice_ridge_ice_cell(1, 1, 1, 1, NCAT, nilyr, 1, p%ridge, p%c%puny, p%c%rhos, dt, ndtd, itermax, HIN, tfx, &
                              conv, shear, a3, v3, s3, T3, q4, si4, qs4, a0, fr, fh, d1, d2, dv, op, ie, nit)
      do n = 1, NCAT
         c%a(n) = a3(1, 1, n); c%v(n) = v3(1, 1, n); c%s(n) = s3(1, 1, n)
         do it = 1, ice_ntrcr(nilyr, 1)
            c%tr(it, n) = ice_trcr_get(it, 1, 1, n, 1, 1, NCAT, nilyr, 1, T3, q4, si4, qs4)
         end do
      end do
      c%a0 = a0(1, 1); c%fresh = fr(1, 1); c%fhocn = fh(1, 1); c%d1 = d1(1, 1); c%d2 = d2(1, 1); c%dv = dv(1, 1)
      c%op = op(1, 1); c%ierr = ie(1, 1); c%niter = nit
   end subroutine run_cell

   pure function area_sum(c) result(asum)
      type(cell_t), intent(in) :: c
      real(wp) :: asum
      asum = c%a0 + sum(c%a)
   end function area_sum

   pure subroutine lipscomb_fracs(hrmin, lam, m, farea, fvol)
      !! Closed form (Lipscomb et al. 2007) of the fractions of a ridge `exp(-(h - hrmin)/lam)` beyond `hrmin` that
      !! land in category `m` of `HIN`, libm exponentials.
      real(wp), intent(in) :: hrmin, lam
      integer, intent(in) :: m
      real(wp), intent(out) :: farea, fvol
      real(wp) :: hl, hr, el, er

      if (hrmin >= HIN(m) .and. m < NCAT) then
         farea = 0.0_wp; fvol = 0.0_wp
         return
      end if
      hl = max(hrmin, HIN(m - 1))
      el = exp(-(hl - hrmin)/lam)
      if (m < NCAT) then
         hr = HIN(m)
         er = exp(-(hr - hrmin)/lam)
         farea = el - er
         fvol = ((hl + lam)*el - (hr + lam)*er)/(hrmin + lam)
      else
         farea = el
         fvol = (hl + lam)*el/(hrmin + lam)
      end if
   end subroutine lipscomb_fracs

   pure subroutine hibler_fracs(hrmin, hrmax, m, farea, fvol)
      !! Closed form (Hibler 1980, modified): ridges uniform on `[hrmin, hrmax]`.
      real(wp), intent(in) :: hrmin, hrmax
      integer, intent(in) :: m
      real(wp), intent(out) :: farea, fvol
      real(wp) :: hl, hr, hb

      hb = HIN(m)
      if (m == NCAT) hb = 1.0e8_wp
      if (hrmin >= hb .or. hrmax <= HIN(m - 1)) then
         farea = 0.0_wp; fvol = 0.0_wp
      else
         hl = max(hrmin, HIN(m - 1)); hr = min(hrmax, hb)
         farea = (hr - hl)/(hrmax - hrmin)
         fvol = (hr*hr - hl*hl)/(hrmax*hrmax - hrmin*hrmin)
      end if
   end subroutine hibler_fracs

   subroutine totals(nilyr, c, vice, vsno, eice, esno, sice, mag_e, mag_es)
      !! Ice volume, snow volume, ice energy `sum q v/nilyr`, snow energy, salt `sum S v/nilyr`, and the sums of
      !! magnitudes of the two energies.
      integer, intent(in) :: nilyr
      type(cell_t), intent(in) :: c
      real(wp), intent(out) :: vice, vsno, eice, esno, sice, mag_e, mag_es
      integer :: n, k

      vice = 0.0_wp; vsno = 0.0_wp; eice = 0.0_wp; esno = 0.0_wp; sice = 0.0_wp; mag_e = 0.0_wp; mag_es = 0.0_wp
      do n = 1, NCAT
         vice = vice + c%v(n)
         vsno = vsno + c%s(n)
         do k = 1, nilyr
            eice = eice + c%tr(1 + k, n)*c%v(n)/real(nilyr, wp)
            mag_e = mag_e + abs(c%tr(1 + k, n)*c%v(n)/real(nilyr, wp))
            sice = sice + c%tr(2 + nilyr + k, n)*c%v(n)/real(nilyr, wp)
         end do
         esno = esno + c%tr(2 + nilyr, n)*c%s(n)
         mag_es = mag_es + abs(c%tr(2 + nilyr, n)*c%s(n))
      end do
   end subroutine totals

   ! ------------------------------------------------------------------
   ! closed forms
   ! ------------------------------------------------------------------

   subroutine test_pure_convergence(error)
      !! A12 case 3: `aice0 = 0.02`, `aicen = (0.6, 0.4)`, `hi = (0.5, 1.5)`, no deformation, `dt = 900`:
      !! `asum = 1.02`, `closing_net = 0.02/dt`, `opning = 0`; participation and ridge geometry from the formulas
      !! of the papers, one pass, area balance exact, ice volume conserved.  Both participation functions.
      type(error_type), allocatable, intent(out) :: error
      integer :: kp
      type(ice_params_t) :: p
      type(cell_t) :: c, c0
      real(wp) :: dt, G(0:2), e(0:2), ap(0:2), norm, hi(2), hrmin(2), lam(2), krdg(2), aksum, cn, cg, a1n(2), a2n(2), v1n(2)
      real(wp) :: fa(2, NCAT), fv(2, NCAT), ea(NCAT), ev(NCAT), a0_exp, d1_exp, d2_exp, dv_exp, worst
      integer :: n, m
      real(wp), parameter :: GSTAR = 0.15_wp, ASTAR = 0.05_wp

      dt = 900.0_wp
      worst = 0.0_wp
      do kp = 0, 1
         call params(ICE_CONSTANTS_ICEPACK, 7, kp, 1, p)
         call blank(c)
         c%a = [0.6_wp, 0.4_wp, 0.0_wp, 0.0_wp, 0.0_wp]
         hi = [0.5_wp, 1.5_wp]
         c%v(1:2) = c%a(1:2)*hi
         c%a0 = 0.02_wp
         call base_tracers(7, p, -10.0_wp, c%tr)
         c0 = c
         call run_cell(p, 7, dt, 1, RIDGE_NITERMAX, -1.8_wp, 0.0_wp, 0.0_wp, c)
         call check(error, c%ierr, 0_int32, "no error bit")
         if (allocated(error)) return
         call check(error, c%niter, 1, "one pass")
         if (allocated(error)) return
         ! G(n) = cumulative area over open water and ice, normalised
         G(0) = 0.02_wp/1.02_wp; G(1) = 0.62_wp/1.02_wp; G(2) = 1.0_wp
         if (kp == 1) then
            norm = 1.0_wp - exp(-1.0_wp/ASTAR)
            e = exp(-G/ASTAR)
            ap(0) = (1.0_wp - e(0))/norm
            ap(1) = (e(0) - e(1))/norm
            ap(2) = (e(1) - e(2))/norm
         else
            ap = 0.0_wp
            ap(0) = (1.0_wp/GSTAR)*G(0)*(2.0_wp - G(0)/GSTAR)
            ap(1) = (1.0_wp/GSTAR)*(GSTAR - G(0))*(2.0_wp - (GSTAR + G(0))/GSTAR)
         end if
         lam = 3.0_wp*sqrt(hi)
         hrmin = [min(2.0_wp*hi(1), hi(1) + 1.0_wp), min(2.0_wp*hi(2), hi(2) + 1.0_wp)]
         krdg = (hrmin + lam)/hi
         aksum = ap(0) + ap(1)*(1.0_wp - 1.0_wp/krdg(1)) + ap(2)*(1.0_wp - 1.0_wp/krdg(2))
         cn = 0.02_wp/dt
         cg = cn/aksum
         call check(error, ap(0)*cg*dt < 0.02_wp, "open-water limiter not active in the hand case")
         if (allocated(error)) return
         do n = 1, 2
            a1n(n) = ap(n)*cg*dt
            a2n(n) = a1n(n)/krdg(n)
            v1n(n) = c0%v(n)*(a1n(n)/c0%a(n))
         end do
         do n = 1, 2
            do m = 1, NCAT
               call lipscomb_fracs(hrmin(n), lam(n), m, fa(n, m), fv(n, m))
            end do
         end do
         ea = c0%a; ev = c0%v
         ea(1) = ea(1) - a1n(1); ea(2) = ea(2) - a1n(2); ev(1) = ev(1) - v1n(1); ev(2) = ev(2) - v1n(2)
         do m = 1, NCAT
            ea(m) = ea(m) + fa(1, m)*a2n(1) + fa(2, m)*a2n(2)
            ev(m) = ev(m) + fv(1, m)*v1n(1) + fv(2, m)*v1n(2)
         end do
         a0_exp = 0.02_wp - ap(0)*cg*dt
         d1_exp = (a1n(1) + a1n(2))/dt; d2_exp = (a2n(1) + a2n(2))/dt; dv_exp = (v1n(1) + v1n(2))/dt
         if (kp == 0) then
            ! Thorndike: only the thinnest category and the open water participate; the hand case above
            ! contains the second category with zero participation
            call check(error, ap(2) == 0.0_wp, "Thorndike participation zero above G*")
            if (allocated(error)) return
         end if
         call check(error, same(c%a0, a0_exp, 1.0e-13_wp, 1.0e-300_wp), "open water after the pass (closed form)")
         if (allocated(error)) return
         do m = 1, NCAT
            call check(error, same(c%a(m), ea(m), 1.0e-13_wp, 1.0e-300_wp), "category area after the pass (closed form)")
            if (allocated(error)) return
            call check(error, same(c%v(m), ev(m), 1.0e-13_wp, 1.0e-300_wp), "category volume after the pass (closed form)")
            if (allocated(error)) return
         end do
         call check(error, same(c%d1, d1_exp, 1.0e-13_wp, 1.0e-300_wp), "dardg1dt")
         if (allocated(error)) return
         call check(error, same(c%d2, d2_exp, 1.0e-13_wp, 1.0e-300_wp), "dardg2dt")
         if (allocated(error)) return
         call check(error, same(c%dv, dv_exp, 1.0e-13_wp, 1.0e-300_wp), "dvirdgdt")
         if (allocated(error)) return
         call check(error, same(c%op, 0.0_wp, 0.0_wp, 0.0_wp), "opening zero (the opening of the pass is 0)")
         if (allocated(error)) return
         ! the asserts of the plan
         call check(error, abs(area_sum(c) - 1.0_wp) < 1.0e-14_wp, "asum = 1 after the pass")
         if (allocated(error)) return
         call check(error, same(sum(c%v), sum(c0%v), 1.0e-14_wp, 1.0e-300_wp), "ice volume conserved")
         if (allocated(error)) return
         call check(error, c%a0 >= 0.0_wp, "open water not negative")
         if (allocated(error)) return
         worst = max(worst, abs(area_sum(c) - 1.0_wp))
      end do
      write (output_unit, '(a,es10.2)') "pure convergence: worst |asum - 1| ", worst
   end subroutine test_pure_convergence

   subroutine single_cover_case(error, kredist, conv, shear, c_out, a1n_out)
      !! The hand case of the plan: `aice0 = 0`, `aicen(1) = 1`, `hi = 1`, closing only through `rdg_conv`/`rdg_shear`.
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: kredist
      real(wp), intent(in) :: conv, shear
      type(cell_t), intent(out) :: c_out
      real(wp), intent(out) :: a1n_out
      type(ice_params_t) :: p
      type(cell_t) :: c, c0
      real(wp) :: dt, cn, cdt, cg, a1, a2, v1, krdg, hrmin, hrmax, lam, fa(NCAT), fv(NCAT), ea(NCAT), ev(NCAT), s1, vs
      integer :: m, k
      real(wp), parameter :: TOLX = 2.0e-14_wp

      dt = 3600.0_wp
      call params(ICE_CONSTANTS_ICEPACK, 7, 1, kredist, p)
      call blank(c)
      c%a(1) = 1.0_wp; c%v(1) = 1.0_wp; c%s(1) = 0.2_wp
      c%a0 = 0.0_wp
      call base_tracers(7, p, -10.0_wp, c%tr)
      do k = 1, 7
         c%tr(1 + k, 1) = -(3.0e8_wp + 2.0e7_wp*real(k, wp))     ! the donor's own enthalpy profile
      end do
      c%tr(1, 1) = -12.0_wp
      c0 = c
      call run_cell(p, 7, dt, 1, RIDGE_NITERMAX, -1.8_wp, conv, shear, c)
      call check(error, c%ierr, 0_int32, "no error bit")
      if (allocated(error)) return
      call check(error, c%niter, 1, "converges in one pass")
      if (allocated(error)) return
      ! closed form: apartic(1) = 1, hrmin = 2, krdg = 5 (redist 1: lam = 3) or 6 (redist 0: hrmax = 10)
      cn = conv + 0.25_wp*shear
      cdt = cn*dt
      hrmin = 2.0_wp
      if (kredist == 1) then
         lam = 3.0_wp
         krdg = (hrmin + lam)/1.0_wp
      else
         hrmax = 2.0_wp*sqrt(25.0_wp*1.0_wp)
         krdg = 0.5_wp*(hrmin + hrmax)/1.0_wp
      end if
      cg = cn/(1.0_wp - 1.0_wp/krdg)
      a1 = cg*dt; a2 = a1/krdg; v1 = a1; vs = 0.2_wp*a1
      a1n_out = a1
      do m = 1, NCAT
         if (kredist == 1) then
            call lipscomb_fracs(hrmin, lam, m, fa(m), fv(m))
         else
            call hibler_fracs(hrmin, hrmax, m, fa(m), fv(m))
         end if
      end do
      ea = fa*a2; ev = fv*v1
      ea(1) = ea(1) + (1.0_wp - a1); ev(1) = ev(1) + (1.0_wp - v1)
      call check(error, same(c%a0, cdt, 1.0e-15_wp, 1.0e-300_wp), "new open water = C dt")
      if (allocated(error)) return
      call check(error, same(c%d1, a1/dt, 1.0e-14_wp, 1.0e-300_wp), "ardg1 = cg dt")
      if (allocated(error)) return
      call check(error, same(c%d2, a2/dt, 1.0e-14_wp, 1.0e-300_wp), "ardg2 = ardg1/krdg")
      if (allocated(error)) return
      call check(error, same(c%dv, v1/dt, 1.0e-14_wp, 1.0e-300_wp), "virdg = ardg1 hi")
      if (allocated(error)) return
      call check(error, same(c%op, cdt/dt, 1.0e-15_wp, 1.0e-300_wp), "opening = opning")
      if (allocated(error)) return
      do m = 1, NCAT
         call check(error, same(c%a(m), ea(m), TOLX, 1.0e-300_wp), "category area (closed form)")
         if (allocated(error)) return
         call check(error, same(c%v(m), ev(m), TOLX, 1.0e-300_wp), "category volume (closed form)")
         if (allocated(error)) return
      end do
      call check(error, abs(area_sum(c) - 1.0_wp) < 1.0e-15_wp, "area balance exact")
      if (allocated(error)) return
      ! categories 1 and 2 receive nothing in the exponential form (hrmin = 2 >= 1.3914)
      if (kredist == 1) then
         call check(error, fa(2) == 0.0_wp .and. fv(2) == 0.0_wp .and. c%a(2) == 0.0_wp, "category 2 receives nothing")
         if (allocated(error)) return
         call check(error, same(sum(fa), 1.0_wp, 1.0e-15_wp, 1.0_wp), "closed-form fractions sum to 1")
         if (allocated(error)) return
      end if
      ! every receiving category carries the donor's own tracers; an empty one the freezing point / zero
      do m = 1, NCAT
         if (ev(m) > 0.0_wp .and. m /= 1 .or. m == 1) then
            do k = 1, 7
               call check(error, same(c%tr(1 + k, m), c0%tr(1 + k, 1), 4.0e-15_wp, 1.0_wp), "volume tracer follows the donor")
               if (allocated(error)) return
            end do
            call check(error, same(c%tr(1, m), c0%tr(1, 1), 4.0e-15_wp, 1.0_wp), "area tracer follows the donor")
            if (allocated(error)) return
         else
            call check(error, c%tr(1, m) == -1.8_wp .and. c%tr(2, m) == 0.0_wp, "empty category: Tf and zero")
            if (allocated(error)) return
         end if
      end do
      s1 = vs
      call check(error, same(c%fresh, SNOW0_RHOS*s1*0.5_wp/dt, 1.0e-14_wp, 1.0e-300_wp), "snow thrown into the ocean")
      if (allocated(error)) return
      c_out = c
   end subroutine single_cover_case

   subroutine test_single_category(error)
      !! A12 case 4: apartic(1) = 1, hrmin = 2, krdg = 5 (Lipscomb: lam = 3) and the Hibler form (hrmax = 10,
      !! krdg = 6); `cg = C/aksum`, `ardg1 = cg dt`, `ardg2 = ardg1/krdg`, the exponential receipts.
      type(error_type), allocatable, intent(out) :: error
      type(cell_t) :: c
      real(wp) :: a1

      call single_cover_case(error, 1, 2.0e-6_wp, 0.0_wp, c, a1)
      if (allocated(error)) return
      call check(error, same(a1, 1.25_wp*2.0e-6_wp*3600.0_wp, 1.0e-14_wp, 1.0e-300_wp), "ardg1 = 1.25 C dt (aksum = 0.8)")
      if (allocated(error)) return
      call single_cover_case(error, 0, 2.0e-6_wp, 0.0_wp, c, a1)
      if (allocated(error)) return
      call check(error, same(a1, 1.2_wp*2.0e-6_wp*3600.0_wp, 1.0e-14_wp, 1.0e-300_wp), "ardg1 = 1.2 C dt (aksum = 5/6)")
   end subroutine test_single_category

   subroutine test_shear_only(error)
      !! Shear only: `closing_net = Cs rdg_shear`, `Cs = 0.25` (Flato and Hibler 1995): the cell given
      !! `rdg_shear = 8e-6` is, bit for bit, the cell given `rdg_conv = 2e-6`; and a mix adds linearly.
      type(error_type), allocatable, intent(out) :: error
      type(cell_t) :: c1, c2, c3
      real(wp) :: a1

      call single_cover_case(error, 1, 2.0e-6_wp, 0.0_wp, c1, a1)
      if (allocated(error)) return
      call single_cover_case(error, 1, 0.0_wp, 8.0e-6_wp, c2, a1)
      if (allocated(error)) return
      call check(error, c1%a0 == c2%a0 .and. all(c1%a == c2%a) .and. all(c1%v == c2%v) .and. c1%d1 == c2%d1, &
                 "shear 8e-6 == convergence 2e-6, bitwise")
      if (allocated(error)) return
      call single_cover_case(error, 1, 1.0e-6_wp, 4.0e-6_wp, c3, a1)
      if (allocated(error)) return
      call check(error, c1%a0 == c3%a0 .and. all(c1%a == c3%a), "convergence 1e-6 + shear 4e-6 == 2e-6, bitwise")
   end subroutine test_shear_only

   subroutine test_divergence(error)
      !! A12 case 5: `asum = 0.97`, no rates: `closing_net = 0`, no ridging, `aice0 += 0.03`; area and volume of
      !! every category bit for bit, tracers re-derived (round-off), one pass, no flux.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      type(cell_t) :: c, c0
      integer :: n, it
      real(wp) :: dt

      dt = 1800.0_wp
      call params(ICE_CONSTANTS_ICEPACK, 7, 1, 1, p)
      call blank(c)
      c%a = [0.4_wp, 0.3_wp, 0.2_wp, 0.0_wp, 0.0_wp]
      c%v = [0.2_wp, 0.45_wp, 0.4_wp, 0.0_wp, 0.0_wp]
      c%s = [0.05_wp, 0.01_wp, 0.0_wp, 0.0_wp, 0.0_wp]
      c%a0 = 0.07_wp
      call base_tracers(7, p, -8.0_wp, c%tr)
      c0 = c
      call run_cell(p, 7, dt, 1, RIDGE_NITERMAX, -1.8_wp, 0.0_wp, 0.0_wp, c)
      call check(error, c%ierr, 0_int32, "no error bit")
      if (allocated(error)) return
      call check(error, c%niter, 1, "one pass")
      if (allocated(error)) return
      call check(error, same(c%a0, 0.07_wp + 0.03_wp, 1.0e-14_wp, 1.0_wp), "aice0 += 1 - asum")
      if (allocated(error)) return
      call check(error, all(c%a == c0%a) .and. all(c%v == c0%v) .and. all(c%s == c0%s), "areas and volumes untouched bit for bit")
      if (allocated(error)) return
      do n = 1, NCAT
         do it = 1, ice_ntrcr(7, 1)
            ! a tracer over an empty carrier (no snow) is normalised to 0, as upstream does
            if (c0%a(n) > 0.0_wp .and. .not. (it == 9 .and. c0%s(n) == 0.0_wp)) then
               call check(error, same(c%tr(it, n), c0%tr(it, n), 4.0e-16_wp, 1.0_wp), "tracers re-derived to round-off")
               if (allocated(error)) return
            end if
         end do
      end do
      call check(error, c%fresh == 0.0_wp .and. c%fhocn == 0.0_wp .and. c%d1 == 0.0_wp .and. c%d2 == 0.0_wp .and. c%dv == 0.0_wp, &
                 "no flux and no ridging diagnostics")
      if (allocated(error)) return
      call check(error, same(c%op, 0.03_wp/dt, 1.0e-13_wp, 1.0e-300_wp), "opening = (1 - asum)/dt")
   end subroutine test_divergence

   subroutine limiter_case(p, itermax, c0, c)
      type(ice_params_t), intent(in) :: p
      integer, intent(in) :: itermax
      type(cell_t), intent(out) :: c0, c

      call blank(c)
      c%a = [0.5_wp, 0.52_wp, 0.0_wp, 0.0_wp, 0.0_wp]
      c%v = [0.15_wp, 0.8_wp, 0.0_wp, 0.0_wp, 0.0_wp]
      c%s = [0.02_wp, 0.05_wp, 0.0_wp, 0.0_wp, 0.0_wp]
      c%a0 = 0.004_wp
      call base_tracers(7, p, -6.0_wp, c%tr)
      c0 = c
      call run_cell(p, 7, 3600.0_wp, 1, itermax, -1.8_wp, 5.0e-5_wp, 0.0_wp, c)
   end subroutine limiter_case

   subroutine test_limiter(error)
      !! A12 case 6: a small `aice0` and a large `rdg_conv` so that `apartic(0) cg dt > aice0`: the first pass is
      !! cut by `tmpfac` (closed form from the pass-1 state with `itermax = 1`), the misfit `(1 - asum)(1 - tmpfac)`
      !! needs a second pass, area and volume conserved; `itermax = 1` returns `ICE_ERR_RIDGE_NOCONV`.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      type(cell_t) :: c0, c1, c2
      real(wp) :: asum0, cn, dvv, G(0:2), e(0:2), ap0, norm, hi(2), hrmin(2), lam(2), krdg(2), aksum, cg, tmpfac, opn
      real(wp), parameter :: ASTAR = 0.05_wp

      call params(ICE_CONSTANTS_ICEPACK, 7, 1, 1, p)
      call limiter_case(p, 1, c0, c1)
      call check(error, iand(c1%ierr, ICE_ERR_RIDGE_NOCONV) /= 0_int32, "itermax = 1: ICE_ERR_RIDGE_NOCONV set")
      if (allocated(error)) return
      call check(error, c1%niter, 1, "one pass run")
      if (allocated(error)) return
      ! pass-1 closed form
      asum0 = c0%a0 + sum(c0%a)
      cn = max(5.0e-5_wp, (asum0 - 1.0_wp)/3600.0_wp)
      opn = cn + (1.0_wp - asum0)/3600.0_wp
      G(0) = c0%a0/asum0; G(1) = (c0%a0 + c0%a(1))/asum0; G(2) = 1.0_wp
      norm = 1.0_wp - exp(-1.0_wp/ASTAR)
      e = exp(-G/ASTAR)
      ap0 = (1.0_wp - e(0))/norm
      hi = [c0%v(1)/c0%a(1), c0%v(2)/c0%a(2)]
      hrmin = [min(2.0_wp*hi(1), hi(1) + 1.0_wp), min(2.0_wp*hi(2), hi(2) + 1.0_wp)]
      lam = 3.0_wp*sqrt(hi)
      krdg = (hrmin + lam)/hi
      aksum = ap0 + ((e(0) - e(1))/norm)*(1.0_wp - 1.0_wp/krdg(1)) + ((e(1) - e(2))/norm)*(1.0_wp - 1.0_wp/krdg(2))
      cg = cn/aksum
      call check(error, ap0*cg*3600.0_wp > c0%a0, "the open-water limiter fires in this case")
      if (allocated(error)) return
      tmpfac = c0%a0/(ap0*cg*3600.0_wp)
      call check(error, same(c1%a0, opn*tmpfac*3600.0_wp, 1.0e-12_wp, 1.0e-300_wp), "pass 1: aice0 = tmpfac opning dt")
      if (allocated(error)) return
      dvv = abs(sum(c1%v) - sum(c0%v))
      call check(error, dvv <= 1.0e-14_wp*sum(c0%v), "pass 1 conserves ice volume")
      if (allocated(error)) return
      call check(error, same(area_sum(c1) - 1.0_wp, (asum0 - 1.0_wp)*(1.0_wp - tmpfac), 1.0e-9_wp, 1.0e-300_wp), &
                 "pass 1 misfit = (asum - 1)(1 - tmpfac)")
      if (allocated(error)) return
      ! full run: converges in two passes
      call limiter_case(p, RIDGE_NITERMAX, c0, c2)
      call check(error, c2%ierr, 0_int32, "full run: no error bit")
      if (allocated(error)) return
      call check(error, c2%niter, 2, "the cut limiter forces a second pass")
      if (allocated(error)) return
      call check(error, abs(area_sum(c2) - 1.0_wp) < PUNY, "asum = 1 within puny after the passes")
      if (allocated(error)) return
      call check(error, c2%a0 >= 0.0_wp, "aice0 not negative")
      if (allocated(error)) return
      call check(error, same(sum(c2%v), sum(c0%v), 1.0e-14_wp, 1.0e-300_wp), "ice volume conserved over both passes")
      if (allocated(error)) return
      write (output_unit, '(a,es10.2,a,es10.2)') "limiter: tmpfac ", tmpfac, ", pass-1 misfit ", area_sum(c1) - 1.0_wp
   end subroutine test_limiter

   subroutine test_snow_signs(error)
      !! A12 case 8: snow at 0 degC (`qsno = -rhos Lfresh`): `fresh > 0`, `fhocn = -Lfresh fresh` to round-off;
      !! half of the ridged snow volume lands on the ridges (`fsnowrdg = 0.5`), the other half is the flux.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      type(cell_t) :: c, c0
      real(wp) :: dt, lf, vsr, dsno

      dt = 3600.0_wp
      call params(ICE_CONSTANTS_ICEPACK, 4, 1, 1, p)
      lf = p%c%Lfresh
      call blank(c)
      c%a = [0.3_wp, 0.4_wp, 0.35_wp, 0.0_wp, 0.0_wp]
      c%v = [0.12_wp, 0.5_wp, 0.7_wp, 0.0_wp, 0.0_wp]
      c%s = [0.06_wp, 0.12_wp, 0.07_wp, 0.0_wp, 0.0_wp]
      c%a0 = 0.0_wp
      call base_tracers(4, p, 0.0_wp, c%tr)
      c0 = c
      call run_cell(p, 4, dt, 1, RIDGE_NITERMAX, -1.8_wp, 0.0_wp, 0.0_wp, c)
      call check(error, c%ierr, 0_int32, "no error bit")
      if (allocated(error)) return
      call check(error, c%fresh > 0.0_wp, "fresh > 0: snow into the ocean")
      if (allocated(error)) return
      call check(error, c%fhocn < 0.0_wp, "fhocn < 0: the ocean cools to melt it")
      if (allocated(error)) return
      call check(error, same(c%fhocn, -lf*c%fresh, 1.0e-14_wp, 1.0e-300_wp), "fhocn = -Lfresh fresh")
      if (allocated(error)) return
      ! half of the ridged snow stays: sum(snow before) - sum(snow after) = (1 - fsnowrdg) ridged volume = rhos^-1 dt fresh
      dsno = sum(c0%s) - sum(c%s)
      vsr = c%fresh*dt/p%c%rhos
      call check(error, same(dsno, vsr, 1.0e-12_wp, 1.0e-300_wp), "snow volume lost = fresh dt/rhos (the thrown half)")
      if (allocated(error)) return
      call check(error, all(c%s >= 0.0_wp), "snow volumes stay non-negative")
   end subroutine test_snow_signs

   subroutine test_no_ridging(error)
      !! A12 case 9: `asum = 1` exactly and no rates: state, fluxes and diagnostics bit for bit, `niter = 0`.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      type(cell_t) :: c, c0

      call params(ICE_CONSTANTS_ICEPACK, 7, 1, 1, p)
      call blank(c)
      c%a = [0.25_wp, 0.25_wp, 0.5_wp, 0.0_wp, 0.0_wp]
      c%v = [0.1_wp, 0.4_wp, 1.5_wp, 0.0_wp, 0.0_wp]
      c%s = [0.01_wp, 0.0_wp, 0.2_wp, 0.0_wp, 0.0_wp]
      c%a0 = 0.0_wp
      c%fresh = 0.3_wp; c%fhocn = -1.25_wp
      call base_tracers(7, p, -8.0_wp, c%tr)
      c%tr(3, 2) = -2.7000000000000001e8_wp   ! a non-trivial value: a re-derived tracer could change by an ulp
      c0 = c
      call run_cell(p, 7, 3600.0_wp, 1, RIDGE_NITERMAX, -1.8_wp, 0.0_wp, 0.0_wp, c)
      call check(error, c%niter, 0, "early-out: no pass")
      if (allocated(error)) return
      call check(error, c%ierr, 0_int32, "no error bit")
      if (allocated(error)) return
      call check(error, all(c%a == c0%a) .and. all(c%v == c0%v) .and. all(c%s == c0%s) .and. all(c%tr == c0%tr) .and. &
                 c%a0 == c0%a0, "state bit for bit")
      if (allocated(error)) return
      call check(error, c%fresh == c0%fresh .and. c%fhocn == c0%fhocn .and. c%d1 == 0.0_wp .and. c%op == 0.0_wp, &
                 "fluxes and diagnostics untouched")
   end subroutine test_no_ridging

   subroutine test_empty_skip(error)
      !! A12 case 13: an ice-free cell (`aice0 = 1`) with deformation rates is left bit for bit (upstream would move
      !! `aice0` by an ulp), also with junk tracers; an ice-free cell with `asum /= 1` is NOT skipped and opens to 1.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      type(cell_t) :: c, c0

      call params(ICE_CONSTANTS_ICEPACK, 7, 1, 1, p)
      call blank(c)
      c%a0 = 1.0_wp
      c%tr = 123.456_wp
      c%fresh = 0.1_wp
      c0 = c
      call run_cell(p, 7, 3600.0_wp, 1, RIDGE_NITERMAX, -1.8_wp, 3.0e-6_wp, 1.0e-6_wp, c)
      call check(error, c%niter, 0, "ice-free cell: no pass")
      if (allocated(error)) return
      call check(error, c%a0 == 1.0_wp .and. all(c%a == 0.0_wp) .and. all(c%tr == c0%tr) .and. c%fresh == 0.1_wp .and. &
                 c%ierr == 0_int32, "ice-free cell untouched")
      if (allocated(error)) return
      call check(error, same(c%op, 3.0e-6_wp + 0.25_wp*1.0e-6_wp, 1.0e-13_wp, 1.0e-300_wp), "the opening diagnostic is booked: min(c, aice0/dt)")
      if (allocated(error)) return
      call blank(c)
      c%a0 = 0.99_wp
      c0 = c
      call run_cell(p, 7, 3600.0_wp, 1, RIDGE_NITERMAX, -1.8_wp, 0.0_wp, 0.0_wp, c)
      call check(error, c%niter, 1, "ice-free cell with asum /= 1 runs its pass")
      if (allocated(error)) return
      call check(error, abs(c%a0 - 1.0_wp) < 1.0e-15_wp, "open water goes to 1")
   end subroutine test_empty_skip

   subroutine test_err_flags(error)
      !! A12 case 12: the error bits and the NaN rule.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      type(cell_t) :: c, c0
      integer :: n
      real(wp) :: nan
      real(wp) :: a3(1, 1, 9), v3(1, 1, 9), s3(1, 1, 9), T3(1, 1, 9), q4(1, 1, 9, 1), si4(1, 1, 9, 1), qs4(1, 1, 9, 1)
      real(wp) :: h9(0:9), a0(1, 1), f1(1, 1), f2(1, 1), f3(1, 1), f4(1, 1), f5(1, 1), f6(1, 1)
      integer(int32) :: ie(1, 1)
      integer :: nit

      nan = ieee_value(1.0_wp, ieee_quiet_nan)
      call params(ICE_CONSTANTS_ICEPACK, 7, 1, 1, p)
      ! a ridging cell
      call blank(c0)
      c0%a = [0.5_wp, 0.5_wp, 0.0_wp, 0.0_wp, 0.0_wp]
      c0%v = [0.2_wp, 0.9_wp, 0.0_wp, 0.0_wp, 0.0_wp]
      c0%s = [0.01_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
      call base_tracers(7, p, -8.0_wp, c0%tr)
      ! NaN in an area, a volume, a snow volume, the open water, a rate, a tracer of the ridging cell: refused, untouched
      do n = 1, 6
         c = c0
         select case (n)
         case (1); c%a(2) = nan
         case (2); c%v(1) = nan
         case (3); c%s(1) = nan
         case (4); c%a0 = nan
         case (5); c%tr(4, 2) = nan
         case (6); c%tr(1, 1) = nan
         end select
         call run_cell(p, 7, 3600.0_wp, 1, RIDGE_NITERMAX, -1.8_wp, 3.0e-6_wp, 0.0_wp, c)
         call check(error, iand(c%ierr, ICE_ERR_NONFINITE) /= 0_int32, "NaN in the state raises ICE_ERR_NONFINITE")
         if (allocated(error)) return
         call check(error, c%niter <= 1 .and. c%fresh == 0.0_wp, "no flux booked")
         if (allocated(error)) return
         select case (n)
         case (1); call check(error, ieee_is_nan(c%a(2)) .and. c%a(1) == c0%a(1), "the NaN area is not laundered")
         case (2); call check(error, ieee_is_nan(c%v(1)) .and. c%v(2) == c0%v(2), "the NaN volume is not laundered")
         case (3); call check(error, ieee_is_nan(c%s(1)), "the NaN snow is not laundered")
         case (4); call check(error, ieee_is_nan(c%a0), "the NaN open water is not laundered")
         case (5); call check(error, ieee_is_nan(c%tr(4, 2)) .and. all(c%a == c0%a) .and. c%a0 == c0%a0, &
 "a NaN tracer: refused before anything is written")
         case (6); call check(error, ieee_is_nan(c%tr(1, 1)) .and. all(c%v == c0%v), "a NaN Tsfc: refused")
         end select
         if (allocated(error)) return
      end do
      c = c0
      call run_cell(p, 7, 3600.0_wp, 1, RIDGE_NITERMAX, -1.8_wp, nan, 0.0_wp, c)
      call check(error, iand(c%ierr, ICE_ERR_NONFINITE) /= 0_int32 .and. all(c%a == c0%a), "NaN rate refused")
      if (allocated(error)) return
      c = c0
      call run_cell(p, 7, 3600.0_wp, 1, RIDGE_NITERMAX, nan, 3.0e-6_wp, 0.0_wp, c)
      call check(error, iand(c%ierr, ICE_ERR_NONFINITE) /= 0_int32 .and. all(c%a == c0%a), "NaN freezing point refused")
      if (allocated(error)) return
      c = c0
      call run_cell(p, 7, 0.0_wp, 1, RIDGE_NITERMAX, -1.8_wp, 3.0e-6_wp, 0.0_wp, c)
      call check(error, iand(c%ierr, ICE_ERR_NONFINITE) /= 0_int32 .and. all(c%a == c0%a), "dt = 0 refused")
      if (allocated(error)) return
      ! ncat above the fixed local size
      a3 = 0.0_wp; v3 = 0.0_wp; s3 = 0.0_wp; T3 = 0.0_wp; q4 = 0.0_wp; si4 = 0.0_wp; qs4 = 0.0_wp; h9 = 0.0_wp
      a0 = 1.0_wp; f1 = 0.0_wp; f2 = 0.0_wp; f3 = 0.0_wp; f4 = 0.0_wp; f5 = 0.0_wp; f6 = 0.0_wp; ie = 0_int32
      call ice_ridge_ice_cell(1, 1, 1, 1, 9, 1, 1, p%ridge, PUNY, 330.0_wp, 3600.0_wp, 1, 20, h9, -1.8_wp, 1.0e-6_wp, 0.0_wp, &
                              a3, v3, s3, T3, q4, si4, qs4, a0, f1, f2, f3, f4, f5, f6, ie, nit)
      call check(error, iand(ie(1, 1), ICE_ERR_NONFINITE) /= 0_int32, "ncat > ICE_NCAT_MAX refused")
      if (allocated(error)) return
      ! open water driven negative: a negative net closing (not produced by the dynamics) with asum = 1
      c = c0
      call run_cell(p, 7, 3600.0_wp, 1, RIDGE_NITERMAX, -1.8_wp, -1.0e-3_wp, 0.0_wp, c)
      call check(error, iand(c%ierr, ICE_ERR_RIDGE_AICE0) /= 0_int32, "ICE_ERR_RIDGE_AICE0")
      if (allocated(error)) return
      call check(error, all(c%a == c0%a) .and. c%fresh == 0.0_wp, "AICE0: nothing written, nothing booked")
      if (allocated(error)) return
      ! inconsistent empty cell (nothing at all) with a closing rate: degenerate, untouched
      call blank(c)
      c%a0 = 0.0_wp
      call base_tracers(7, p, -8.0_wp, c%tr)
      c0 = c
      call run_cell(p, 7, 3600.0_wp, 1, RIDGE_NITERMAX, -1.8_wp, 1.0e-6_wp, 0.0_wp, c)
      call check(error, iand(c%ierr, ICE_ERR_RIDGE_DEGENERATE) /= 0_int32, "ICE_ERR_RIDGE_DEGENERATE")
      if (allocated(error)) return
      call check(error, c%a0 == 0.0_wp .and. all(c%a == 0.0_wp) .and. c%niter == 1, "degenerate cell untouched")
      if (allocated(error)) return
      ! the same empty cell without a closing rate simply opens to 1
      c = c0
      call run_cell(p, 7, 3600.0_wp, 1, RIDGE_NITERMAX, -1.8_wp, 0.0_wp, 0.0_wp, c)
      call check(error, c%ierr == 0_int32 .and. abs(c%a0 - 1.0_wp) < 1.0e-15_wp, "empty cell, no closing: opens to open water")
      if (allocated(error)) return
      ! non-convergence: tested in limiter_forces_second_pass (itermax = 1)
      call limiter_case(p, 1, c0, c)
      call check(error, iand(c%ierr, ICE_ERR_RIDGE_NOCONV) /= 0_int32 .and. c%fresh == 0.0_wp, "NOCONV: no flux booked")
   end subroutine test_err_flags

   ! ------------------------------------------------------------------
   ! random states
   ! ------------------------------------------------------------------

   subroutine fill_random(nx, ny, nilyr, seed, p, coherent, a, v, s, Tsfcn, qicen, sicen, qsnon, a0, conv, shear)
      !! A valid post-transport state on every cell: ~25 % of the categories empty, mean thickness 0.05-4 m,
      !! total ice area 0.85-1.08 (open water `max(1 - tot, 0)`, perturbed in 20 % of the cells so that `asum /= 1`),
      !! snow up to 0.4 m on 70 % of the categories, deformation rates in 60 % of the cells.  `coherent`: the global
      !! case, ~85 % of the cells ice free in 32 x 32 blocks, ice cells with `asum` 1 +- 1 %.
      integer, intent(in) :: nx, ny, nilyr
      integer(int64), intent(in) :: seed
      type(ice_params_t), intent(in) :: p
      logical, intent(in) :: coherent
      real(wp), intent(out) :: a(nx, ny, NCAT), v(nx, ny, NCAT), s(nx, ny, NCAT), Tsfcn(nx, ny, NCAT)
      real(wp), intent(out) :: qicen(nx, ny, NCAT, nilyr), sicen(nx, ny, NCAT, nilyr), qsnon(nx, ny, NCAT, 1)
      real(wp), intent(out) :: a0(nx, ny), conv(nx, ny), shear(nx, ny)
      integer(int64) :: st
      integer :: i, j, n, k
      real(wp) :: u1, u2, u3, tot, tgt, hmax, T
      logical :: ice

      do j = 1, ny
         do i = 1, nx
            ice = .true.
            if (coherent) then
               st = lcg_next(seed + int((i - 1)/32, int64)*7919_int64 + int((j - 1)/32, int64)*104729_int64*65537_int64)
               st = lcg_next(st)
               ice = lcg_unit(st) < 0.15_wp
            end if
            st = lcg_next(seed + int(i, int64)*31_int64 + int(j, int64)*104729_int64*131_int64)
            st = lcg_next(st)
            tot = 0.0_wp
            hmax = 4.0_wp
            do n = 1, NCAT
               a(i, j, n) = 0.0_wp; v(i, j, n) = 0.0_wp; s(i, j, n) = 0.0_wp
               Tsfcn(i, j, n) = -1.8_wp; qicen(i, j, n, :) = 0.0_wp; sicen(i, j, n, :) = 0.0_wp; qsnon(i, j, n, 1) = 0.0_wp
               if (ice) then
                  st = lcg_next(st); u1 = lcg_unit(st)
                  st = lcg_next(st); u2 = lcg_unit(st)
                  st = lcg_next(st); u3 = lcg_unit(st)
                  if (u1 < 0.75_wp) then
                     a(i, j, n) = 0.01_wp + u2
                     v(i, j, n) = a(i, j, n)*(0.05_wp + hmax*u3*u3)
                     tot = tot + a(i, j, n)
                     st = lcg_next(st); u1 = lcg_unit(st)
                     st = lcg_next(st); u2 = lcg_unit(st)
                     if (u1 < 0.7_wp) s(i, j, n) = a(i, j, n)*0.4_wp*u2
                     st = lcg_next(st); Tsfcn(i, j, n) = -30.0_wp*lcg_unit(st)
                     do k = 1, nilyr
                        st = lcg_next(st); qicen(i, j, n, k) = -(2.0e8_wp + 2.0e8_wp*lcg_unit(st))
                        st = lcg_next(st); sicen(i, j, n, k) = 1.0_wp + 5.0_wp*lcg_unit(st)
                     end do
                     st = lcg_next(st); T = -1.0_wp - 24.0_wp*lcg_unit(st)
                     qsnon(i, j, n, 1) = ice_qsno_init_from_T(T, p%c)
                  end if
               end if
            end do
            st = lcg_next(st); u1 = lcg_unit(st)
            st = lcg_next(st); u2 = lcg_unit(st)
            st = lcg_next(st); u3 = lcg_unit(st)
            if (coherent) then
               tgt = 0.99_wp + 0.02_wp*u1
            else
               tgt = 0.85_wp + 0.23_wp*u1
            end if
            if (tot > 0.0_wp) then
               a(i, j, :) = a(i, j, :)*(tgt/tot); v(i, j, :) = v(i, j, :)*(tgt/tot); s(i, j, :) = s(i, j, :)*(tgt/tot)
               tot = sum(a(i, j, :))
            end if
            a0(i, j) = max(1.0_wp - tot, 0.0_wp)
            if (.not. coherent .and. u2 < 0.2_wp) a0(i, j) = a0(i, j)*(0.6_wp + 0.6_wp*u3)
            conv(i, j) = 0.0_wp; shear(i, j) = 0.0_wp
            st = lcg_next(st); u1 = lcg_unit(st)
            st = lcg_next(st); u2 = lcg_unit(st)
            st = lcg_next(st); u3 = lcg_unit(st)
            if (tot > 0.0_wp .and. u1 < 0.6_wp) then
               conv(i, j) = 3.0e-5_wp*u2*u2
               shear(i, j) = 2.0e-5_wp*u3*u3
            end if
         end do
      end do
   end subroutine fill_random

   subroutine cell_to_packed(i, j, nx, ny, nilyr, a, v, s, Tsfcn, qicen, sicen, qsnon, a0, c)
      integer, intent(in) :: i, j, nx, ny, nilyr
      real(wp), intent(in) :: a(nx, ny, NCAT), v(nx, ny, NCAT), s(nx, ny, NCAT), Tsfcn(nx, ny, NCAT)
      real(wp), intent(in) :: qicen(nx, ny, NCAT, nilyr), sicen(nx, ny, NCAT, nilyr), qsnon(nx, ny, NCAT, 1), a0(nx, ny)
      type(cell_t), intent(out) :: c
      integer :: n, it

      call blank(c)
      do n = 1, NCAT
         c%a(n) = a(i, j, n); c%v(n) = v(i, j, n); c%s(n) = s(i, j, n)
         do it = 1, ice_ntrcr(nilyr, 1)
            c%tr(it, n) = ice_trcr_get(it, i, j, n, nx, ny, NCAT, nilyr, 1, Tsfcn, qicen, sicen, qsnon)
         end do
      end do
      c%a0 = a0(i, j)
   end subroutine cell_to_packed

   ! ------------------------------------------------------------------
   ! conservation on random columns
   ! ------------------------------------------------------------------

   subroutine conserve_run(error, nilyr)
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: nilyr
      integer, parameter :: nx = 500, ny = 200
      real(wp), allocatable :: a(:, :, :), v(:, :, :), s(:, :, :), Tsfcn(:, :, :), qicen(:, :, :, :), sicen(:, :, :, :)
      real(wp), allocatable :: qsnon(:, :, :, :), a0(:, :), conv(:, :), shear(:, :), fr(:, :), fh(:, :), d1(:, :), d2(:, :)
      real(wp), allocatable :: dv(:, :), op(:, :)
      real(wp), allocatable :: a_0(:, :, :), v_0(:, :, :), s_0(:, :, :), T_0(:, :, :), q_0(:, :, :, :), si_0(:, :, :, :)
      real(wp), allocatable :: qs_0(:, :, :, :), a0_0(:, :)
      integer(int32), allocatable :: ierr(:, :)
      type(ice_params_t) :: p4(0:3)
      type(cell_t) :: c0, c1
      real(wp) :: dt, vi0, vs0, ei0, es0, si0, me0, mes0, vi1, vs1, ei1, es1, si1, me1, mes1, dfr, dfh
      real(wp) :: w_v, w_s, w_e, w_sn, w_es, w_area, r, amin, a0min
      integer :: i, j, k, nit, nerr, nmulti, nridge, nquiet, nmax, ncnt(0:RIDGE_NITERMAX)
      character(len=160) :: msg

      dt = 3600.0_wp
      do k = 0, 3
         call params(ICE_CONSTANTS_ICEPACK, nilyr, k/2, mod(k, 2), p4(k))
      end do
      allocate (a(nx, ny, NCAT), v(nx, ny, NCAT), s(nx, ny, NCAT), Tsfcn(nx, ny, NCAT), qicen(nx, ny, NCAT, nilyr), &
                sicen(nx, ny, NCAT, nilyr), qsnon(nx, ny, NCAT, 1), a0(nx, ny), conv(nx, ny), shear(nx, ny), fr(nx, ny), &
                fh(nx, ny), d1(nx, ny), d2(nx, ny), dv(nx, ny), op(nx, ny), ierr(nx, ny))
      allocate (a_0(nx, ny, NCAT), v_0(nx, ny, NCAT), s_0(nx, ny, NCAT), T_0(nx, ny, NCAT), q_0(nx, ny, NCAT, nilyr), &
                si_0(nx, ny, NCAT, nilyr), qs_0(nx, ny, NCAT, 1), a0_0(nx, ny))
      call fill_random(nx, ny, nilyr, 90210_int64 + int(nilyr, int64), p4(0), .false., a, v, s, Tsfcn, qicen, sicen, qsnon, a0, &
                       conv, shear)
      a_0 = a; v_0 = v; s_0 = s; T_0 = Tsfcn; q_0 = qicen; si_0 = sicen; qs_0 = qsnon; a0_0 = a0
      fr = 0.0_wp; fh = 0.0_wp; d1 = 0.0_wp; d2 = 0.0_wp; dv = 0.0_wp; op = 0.0_wp; ierr = 0_int32
      w_v = 0.0_wp; w_s = 0.0_wp; w_e = 0.0_wp; w_sn = 0.0_wp; w_es = 0.0_wp; w_area = 0.0_wp
      nerr = 0; nmulti = 0; nridge = 0; nquiet = 0; nmax = 0; ncnt = 0; amin = 0.0_wp; a0min = 1.0_wp
      do j = 1, ny
         do i = 1, nx
            k = mod(i + j, 4)
            call ice_ridge_ice_cell(i, j, nx, ny, NCAT, nilyr, 1, p4(k)%ridge, p4(k)%c%puny, p4(k)%c%rhos, dt, 1, &
                                    RIDGE_NITERMAX, HIN, -1.8_wp, conv(i, j), shear(i, j), a, v, s, Tsfcn, qicen, sicen, qsnon, &
                                    a0, fr, fh, d1, d2, dv, op, ierr, nit)
            if (ierr(i, j) /= 0_int32) then
               nerr = nerr + 1
               if (nerr <= 3) write (output_unit, '(a,2i5,a,i0)') "   cell with error word: ", i, j, " ierr ", ierr(i, j)
               cycle
            end if
            ncnt(nit) = ncnt(nit) + 1
            if (nit >= 2) nmulti = nmulti + 1
            if (nit >= 1) nridge = nridge + 1
            if (nit == 0) nquiet = nquiet + 1
            nmax = max(nmax, nit)
            call cell_to_packed(i, j, nx, ny, nilyr, a_0, v_0, s_0, T_0, q_0, si_0, qs_0, a0_0, c0)
            call cell_to_packed(i, j, nx, ny, nilyr, a, v, s, Tsfcn, qicen, sicen, qsnon, a0, c1)
            call totals(nilyr, c0, vi0, vs0, ei0, es0, si0, me0, mes0)
            call totals(nilyr, c1, vi1, vs1, ei1, es1, si1, me1, mes1)
            if (.not. (ieee_is_finite(vi1 + vs1 + ei1 + es1 + si1 + a0(i, j)))) then
               call test_failed(error, "non-finite state after the ridging step")
               return
            end if
            dfr = fr(i, j); dfh = fh(i, j)
            r = abs(vi1 - vi0)/max(sum(abs(c0%v)), tiny(1.0_wp)); w_v = max(w_v, r)
            r = abs(si1 - si0)/max(sum(abs(c0%v))*7.0_wp, tiny(1.0_wp)); w_s = max(w_s, r)
            r = abs(ei1 - ei0)/max(me0, tiny(1.0_wp)); w_e = max(w_e, r)
            r = abs((vs0 - vs1) - dfr*dt/p4(k)%c%rhos)/max(vs0, tiny(1.0_wp)); w_sn = max(w_sn, r)
            r = abs((es0 - es1) - dfh*dt)/max(mes0, tiny(1.0_wp)); w_es = max(w_es, r)
            if (nit >= 1) w_area = max(w_area, abs(c1%a0 + sum(c1%a) - 1.0_wp))
            amin = min(amin, minval(c1%a))
            a0min = min(a0min, c1%a0)
            if (.not. (a0(i, j) >= 0.0_wp)) then
               call test_failed(error, "open water negative")
               return
            end if
         end do
      end do
      write (output_unit, '(a,i0,a,i0,a,i0,a,i0,a,i0,a)') "ridging conservation, nilyr ", nilyr, ": cells ", nx*ny, ", ridged ", &
         nridge, ", multi-pass ", nmulti, ", quiescent/ice-free ", nquiet, " (max passes "
      write (output_unit, '(i0,a,*(i0,1x))') nmax, ")  pass histogram 0..: ", ncnt(0:min(nmax, 6))
      write (output_unit, '(a,5es10.2)') "  worst rel: ice vol, salt, ice energy, snow vol vs fresh, snow energy vs fhocn: ", &
         w_v, w_s, w_e, w_sn, w_es
      write (output_unit, '(a,es10.2,a,es10.2,a,es10.2)') "  worst |asum - 1| ", w_area, "; min category area ", amin, "; min aice0 ", a0min
      write (msg, '(a,i0)') "no error bit on any random cell (fatal-bit cells: ", nerr
      call check(error, nerr == 0, trim(msg)//")")
      if (allocated(error)) return
      call check(error, nmulti > 100, "the random set contains multi-pass cells")
      if (allocated(error)) return
      call check(error, w_v <= CLOSE_TOL .and. w_s <= CLOSE_TOL, "ice volume and salt conserved")
      if (allocated(error)) return
      call check(error, w_e <= CLOSE_TOL, "ice energy conserved")
      if (allocated(error)) return
      call check(error, w_sn <= CLOSE_TOL, "snow volume: before = after + dt fresh/rhos")
      if (allocated(error)) return
      call check(error, w_es <= CLOSE_TOL, "snow energy: before = after + dt fhocn")
      if (allocated(error)) return
      call check(error, w_area < PUNY, "total area = 1 within puny after the step")
      if (allocated(error)) return
      call check(error, amin >= -PUNY, "no category area below -puny")
   end subroutine conserve_run

   subroutine test_conserve_7(error)
      type(error_type), allocatable, intent(out) :: error
      call conserve_run(error, 7)
   end subroutine test_conserve_7

   subroutine test_conserve_4(error)
      type(error_type), allocatable, intent(out) :: error
      call conserve_run(error, 4)
   end subroutine test_conserve_4

   subroutine test_conserve_1(error)
      type(error_type), allocatable, intent(out) :: error
      call conserve_run(error, 1)
   end subroutine test_conserve_1

   ! ------------------------------------------------------------------
   ! oracle
   ! ------------------------------------------------------------------

   subroutine oracle_run(error, kset, gin, gold)
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: kset
      real(wp), intent(in) :: gin(100, NR), gold(110, NR)
      type(ice_params_t) :: p
      type(cell_t) :: c
      integer :: r, nilyr, ntr, kp, kr, ndtd, it, n, nbad, nmulti(0:3), ncomb, nsnow, nthrown, nquiet
      real(wp) :: dt, tfx, tol, got(110), worst(0:3), d, fl
      character(len=240) :: msg

      call fill_golden()
      worst = 0.0_wp
      nmulti = 0
      nsnow = 0
      nthrown = 0
      nquiet = 0
      do r = 1, NR
         nilyr = GMETA(4, r); ndtd = GMETA(5, r); kp = GMETA(2, r); kr = GMETA(3, r)
         ntr = ice_ntrcr(nilyr, 1)
         call params(kset, nilyr, kp, kr, p)
         dt = GHD(1, r); tfx = GHD(2, r)
         call blank(c)
         c%a0 = GHD(3, r)
         c%a = gin(1:5, r); c%v = gin(6:10, r); c%s = gin(11:15, r)
         c%tr(1:ntr, :) = reshape(gin(16:15 + 5*ntr, r), [ntr, NCAT])
         c%fresh = 0.25_wp; c%fhocn = -3.5_wp
         call run_cell(p, nilyr, dt, ndtd, RIDGE_NITERMAX, tfx, GHD(4, r), GHD(5, r), c)
         call check(error, c%ierr, 0_int32, "oracle case: no error bit")
         if (allocated(error)) return
         ncomb = 2*kp + kr
         if (c%niter >= 2) nmulti(ncomb) = nmulti(ncomb) + 1
         if (c%fresh /= 0.25_wp) nthrown = nthrown + 1
         got = 0.0_wp
         got(1) = c%a0; got(2:6) = c%a; got(7:11) = c%v; got(12:16) = c%s
         got(17) = c%fresh; got(18) = c%fhocn; got(19) = c%d1; got(20) = c%d2; got(21) = c%dv; got(22) = c%op
         got(23:22 + 5*ntr) = reshape(c%tr(1:ntr, :), [ntr*NCAT])
         tol = TOL_GOLD
         if (.not. GPU_BUILD .and. .not. (kp == 0 .and. kr == 0)) tol = TOL_EXP
         nbad = 0
         do it = 1, 22 + 5*ntr
            ! micro categories keep garbage intensive tracers (a tiny remainder over a tiny volume): judged only when
            ! healthy.  A QUIESCENT cell (no pass: nothing opens or closes) is left alone by the port while upstream
            ! re-derives every tracer (a round-off `x a/a`) and normalises the empty categories: tracers of healthy
            ! categories to 1e-15, the rest not compared (areas, volumes and fluxes are bitwise as everywhere).
            if (it > 22) then
               n = (it - 23)/ntr + 1
              if ((tol /= 0.0_wp .or. c%niter == 0) .and. .not. (gold(1 + n, r) > 1.0e-6_wp .and. gold(6 + n, r) > 1.0e-6_wp)) cycle
               if (c%niter == 0) then
                  nquiet = nquiet + 1
                  if (.not. same(got(it), gold(it, r), max(tol, 1.0e-15_wp), 1.0_wp)) then
                     call test_failed(error, "oracle quiescent cell: tracer of a healthy category differs by more than 1e-15")
                     return
                  end if
                  cycle
               end if
            end if
            if (.not. (ieee_is_finite(got(it)) .and. ieee_is_finite(gold(it, r)))) then
               write (msg, '(a,i0,a,i0)') "oracle row ", r, ": non-finite value in output ", it
               call test_failed(error, trim(msg))
               return
            end if
            d = 0.0_wp
            if (gold(it, r) /= 0.0_wp) d = abs(got(it) - gold(it, r))/abs(gold(it, r))
            worst(ncomb) = max(worst(ncomb), d)
            fl = 1.0e-300_wp
            if (it <= 16 .or. it == 17 .or. it == 18) fl = 1.0e-3_wp
            if (.not. same(got(it), gold(it, r), tol, fl)) then
               nbad = nbad + 1
               if (nbad <= 6) write (output_unit, '(a,i0,a,i0,a,i0,a,i0,a,es24.16,a,es24.16)') "   oracle row ", r, " (kp kr ", &
                  ncomb, ") passes ", c%niter, " output ", it, " got ", got(it), " expected ", gold(it, r)
               if (nbad == 1) then
                  write (msg, '(a,i0,a,i0,a,i0,a,es24.16,a,es24.16)') "oracle row ", r, " (kp kr ", ncomb, ") output ", it, &
                     " got ", got(it), " expected ", gold(it, r)
                  if (.not. allocated(error)) call test_failed(error, trim(msg))
               end if
            end if
         end do
      end do
      write (output_unit, '(a,4es10.2)') "oracle worst relative difference (partic x redist 00 01 10 11), every output: ", worst
      write (output_unit, '(a,4i4,a,i0,a,i0)') "oracle multi-pass cells per combination: ", nmulti, "; cells that throw snow: ", nthrown, &
         "; tracers of quiescent cells compared: ", nquiet
      call check(error, all(nmulti >= 1), "every combination has a multi-pass cell in the oracle set")
      if (allocated(error)) return
      call check(error, nthrown >= 24, "the oracle set exercises the snow flux")
   end subroutine oracle_run

   subroutine test_oracle_icepack(error)
      type(error_type), allocatable, intent(out) :: error
      call fill_golden()
      call oracle_run(error, ICE_CONSTANTS_ICEPACK, GII, GOI)
   end subroutine test_oracle_icepack

   subroutine test_oracle_cesm(error)
      type(error_type), allocatable, intent(out) :: error
      call fill_golden()
      call oracle_run(error, ICE_CONSTANTS_CESM, GIC, GOC)
   end subroutine test_oracle_cesm

   ! ------------------------------------------------------------------
   ! device sweep and throughput
   ! ------------------------------------------------------------------

   subroutine test_device(error)
      !! `ice_step_ridge` over an INTERIOR range of a random 64 x 48 state against a host reference (the ridging cell
      !! kernel then the clean-up cell kernel), EVERY output of EVERY interior cell, the ring untouched.  Decisions
      !! (error word, pass count, the zero pattern of every category) are compared EXACTLY, extensive values and
      !! flux increments relatively, healthy-category tracers at their own tolerance.  Production parameters
      !! (Lipscomb x Lipscomb) and the arithmetic-only pair.  Then the throughput.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p

      call params(ICE_CONSTANTS_ICEPACK, 7, 1, 1, p)
      call device_sweep(error, p)
      if (allocated(error)) return
      call params(ICE_CONSTANTS_ICEPACK, 7, 0, 0, p)
      call device_sweep(error, p)
      if (allocated(error)) return
      call params(ICE_CONSTANTS_ICEPACK, 7, 1, 1, p)
      call throughput(error, p, 1)
      if (allocated(error)) return
      call throughput(error, p, 2)
   end subroutine test_device

   subroutine device_sweep(error, p)
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t), intent(in) :: p
      integer, parameter :: nx = 64, ny = 48, nilyr = 7
      real(wp), allocatable :: a(:, :, :), v(:, :, :), s(:, :, :), Tsfcn(:, :, :), qicen(:, :, :, :), sicen(:, :, :, :)
      real(wp), allocatable :: qsnon(:, :, :, :), aice(:, :), vice(:, :), vsno(:, :), aice0(:, :), Tc(:, :), qc(:, :, :)
      real(wp), allocatable :: sc(:, :, :), qsc(:, :, :), fr(:, :), fs(:, :), fh(:, :), d1(:, :), d2(:, :), dv(:, :), op(:, :)
      real(wp), allocatable :: Tf2(:, :), conv(:, :), shear(:, :)
      real(wp), allocatable :: a_0(:, :, :), v_0(:, :, :), s_0(:, :, :), T_0(:, :, :), q_0(:, :, :, :), si_0(:, :, :, :)
      real(wp), allocatable :: qs_0(:, :, :, :), a0_0(:, :)
      real(wp), allocatable :: a_r(:, :, :), v_r(:, :, :), s_r(:, :, :), T_r(:, :, :), q_r(:, :, :, :), si_r(:, :, :, :)
      real(wp), allocatable :: qs_r(:, :, :, :), a0_r(:, :), ai_r(:, :), vi_r(:, :), vs_r(:, :), Tc_r(:, :), qc_r(:, :, :)
      real(wp), allocatable :: sc_r(:, :, :), qsc_r(:, :, :), fr_r(:, :), fs_r(:, :), fh_r(:, :), d1_r(:, :), d2_r(:, :)
      real(wp), allocatable :: dv_r(:, :), op_r(:, :)
      integer(int32), allocatable :: ierr(:, :), ierr_r(:, :), niter(:, :), niter_r(:, :)
      logical, allocatable :: wet(:, :)
      real(wp) :: hin_loc(0:NCAT), dt
      integer :: i, j, n, it, nflip, nring, nridge, nmulti, nit, wi(4), wj(4), nambig, nskip
      logical :: inside, ambig, ok
      real(wp) :: wv(4), asum_h, mscale
      character(len=16) :: wn(4)
      type(cell_t) :: c

      allocate (a(nx, ny, NCAT), v(nx, ny, NCAT), s(nx, ny, NCAT), Tsfcn(nx, ny, NCAT), qicen(nx, ny, NCAT, nilyr), &
                sicen(nx, ny, NCAT, nilyr), qsnon(nx, ny, NCAT, 1), aice(nx, ny), vice(nx, ny), vsno(nx, ny), aice0(nx, ny), &
                Tc(nx, ny), qc(nx, ny, nilyr), sc(nx, ny, nilyr), qsc(nx, ny, 1), fr(nx, ny), fs(nx, ny), fh(nx, ny), &
                d1(nx, ny), d2(nx, ny), dv(nx, ny), op(nx, ny), Tf2(nx, ny), conv(nx, ny), shear(nx, ny), ierr(nx, ny), &
                niter(nx, ny), wet(nx, ny))
      allocate (a_0(nx, ny, NCAT), v_0(nx, ny, NCAT), s_0(nx, ny, NCAT), T_0(nx, ny, NCAT), q_0(nx, ny, NCAT, nilyr), &
                si_0(nx, ny, NCAT, nilyr), qs_0(nx, ny, NCAT, 1), a0_0(nx, ny))
      allocate (a_r(nx, ny, NCAT), v_r(nx, ny, NCAT), s_r(nx, ny, NCAT), T_r(nx, ny, NCAT), q_r(nx, ny, NCAT, nilyr), &
                si_r(nx, ny, NCAT, nilyr), qs_r(nx, ny, NCAT, 1), a0_r(nx, ny), ai_r(nx, ny), vi_r(nx, ny), vs_r(nx, ny), &
                Tc_r(nx, ny), qc_r(nx, ny, nilyr), sc_r(nx, ny, nilyr), qsc_r(nx, ny, 1), fr_r(nx, ny), fs_r(nx, ny), &
                fh_r(nx, ny), d1_r(nx, ny), d2_r(nx, ny), dv_r(nx, ny), op_r(nx, ny), ierr_r(nx, ny), niter_r(nx, ny))
      call fill_random(nx, ny, nilyr, 424242_int64, p, .false., a, v, s, Tsfcn, qicen, sicen, qsnon, aice0, conv, shear)
      a_0 = a; v_0 = v; s_0 = s; T_0 = Tsfcn; q_0 = qicen; si_0 = sicen; qs_0 = qsnon; a0_0 = aice0
      hin_loc = HIN
      Tf2 = -1.8_wp
      dt = 2400.0_wp
      do j = 1, ny
         do i = 1, nx
            wet(i, j) = mod(i + 2*j, 17) /= 0     ! a scattering of land cells
            if (i >= 2 .and. i <= nx - 1 .and. j >= 2 .and. j <= ny - 1) then
               fr(i, j) = 0.0_wp; fs(i, j) = 0.0_wp; fh(i, j) = 0.0_wp
               d1(i, j) = 0.0_wp; d2(i, j) = 0.0_wp; dv(i, j) = 0.0_wp; op(i, j) = 0.0_wp
            else
               fr(i, j) = 0.01_wp*real(i, wp) - 0.02_wp*real(j, wp) + 0.5_wp
               fs(i, j) = 0.03_wp*real(j, wp) - 0.5_wp
               fh(i, j) = 0.1_wp*real(i + j, wp)
               d1(i, j) = 1.0e-7_wp*real(i, wp); d2(i, j) = 2.0e-7_wp*real(j, wp); dv(i, j) = 3.0e-7_wp; op(i, j) = 4.0e-7_wp
            end if
         end do
      end do
      aice = -777.0_wp; vice = -777.0_wp; vsno = -777.0_wp
      Tc = -777.0_wp; qc = -777.0_wp; sc = -777.0_wp; qsc = -777.0_wp
      ierr = 0_int32; niter = -777_int32
      do j = 2, ny - 1
         do i = 2, nx - 1
            if (.not. wet(i, j)) Tf2(i, j) = ieee_value(1.0_wp, ieee_quiet_nan)   ! land: its freezing point is never read
            if (wet(i, j) .and. mod(7*i + 3*j, 29) == 0) ierr(i, j) = ICE_ERR_RIDGE_NOCONV   ! a cell that arrives fatal
         end do
      end do
      ! host reference on copies
      a_r = a; v_r = v; s_r = s; T_r = Tsfcn; q_r = qicen; si_r = sicen; qs_r = qsnon; a0_r = aice0
      ai_r = aice; vi_r = vice; vs_r = vsno; Tc_r = Tc; qc_r = qc; sc_r = sc; qsc_r = qsc
      fr_r = fr; fs_r = fs; fh_r = fh; d1_r = d1; d2_r = d2; dv_r = dv; op_r = op; ierr_r = ierr; niter_r = niter
      nridge = 0; nmulti = 0
      do j = 2, ny - 1
         do i = 2, nx - 1
            niter_r(i, j) = 0
            if (wet(i, j) .and. iand(ierr_r(i, j), ICE_ERR_FATAL_MASK) == 0_int32) then
               call ice_ridge_ice_cell(i, j, nx, ny, NCAT, nilyr, 1, p%ridge, p%c%puny, p%c%rhos, dt, 1, RIDGE_NITERMAX, hin_loc, &
                                       Tf2(i, j), conv(i, j), shear(i, j), a_r, v_r, s_r, T_r, q_r, si_r, qs_r, a0_r, fr_r, fh_r, &
                                       d1_r, d2_r, dv_r, op_r, ierr_r, nit)
               niter_r(i, j) = nit
               if (nit >= 1) nridge = nridge + 1
               if (nit >= 2) nmulti = nmulti + 1
               if (iand(ierr_r(i, j), ICE_ERR_FATAL_MASK) == 0_int32) &
                  call ice_cleanup_itd_cell(i, j, nx, ny, NCAT, nilyr, 1, hin_loc, dt, Tf2(i, j), p%c%rhoi, p%c%rhos, p%c%Lfresh, &
                                            p%c%cp_ice, p%c%hs_min, p%c%puny, p%c%itd_area_min, p%c%itd_mass_min, &
                                            p%thermo%ice_ref_salinity, .true., .true., .true., a_r, v_r, s_r, T_r, q_r, si_r, &
                                            qs_r, ai_r, vi_r, vs_r, a0_r, Tc_r, qc_r, sc_r, qsc_r, fr_r, fs_r, fh_r, ierr_r)
            end if
         end do
      end do
      !$acc enter data copyin(hin_loc, wet, Tf2, conv, shear, a, v, s, Tsfcn, qicen, sicen, qsnon, aice, vice, vsno, aice0, Tc, qc, sc, qsc, fr, fs, fh, d1, d2, dv, op, niter, ierr)
      call ice_step_ridge(nx, ny, NCAT, nilyr, 1, 2, nx - 1, 2, ny - 1, p, hin_loc, dt, 1, wet, Tf2, conv, shear, a, v, s, &
                          Tsfcn, qicen, sicen, qsnon, aice, vice, vsno, aice0, Tc, qc, sc, qsc, fr, fs, fh, d1, d2, dv, op, &
                          niter, ierr, ok)
      !$acc update self(a, v, s, Tsfcn, qicen, sicen, qsnon, aice, vice, vsno, aice0, Tc, qc, sc, qsc, fr, fs, fh, d1, d2, dv, op, niter, ierr)
      !$acc exit data delete(hin_loc, wet, Tf2, conv, shear, a, v, s, Tsfcn, qicen, sicen, qsnon, aice, vice, vsno, aice0, Tc, qc, sc, qsc, fr, fs, fh, d1, d2, dv, op, niter, ierr)
      call check(error, ok, "the driver accepts the configuration")
      if (allocated(error)) return
      nring = 0; nflip = 0; nambig = 0; nskip = 0
      wv = 0.0_wp; wi = 0; wj = 0; wn = ''
      do j = 1, ny
         do i = 1, nx
            inside = i >= 2 .and. i <= nx - 1 .and. j >= 2 .and. j <= ny - 1
            if (.not. inside) then
               do n = 1, NCAT
                  if (a(i, j, n) /= a_0(i, j, n) .or. v(i, j, n) /= v_0(i, j, n) .or. s(i, j, n) /= s_0(i, j, n)) nring = nring + 1
                  if (Tsfcn(i, j, n) /= T_0(i, j, n) .or. any(qicen(i, j, n, :) /= q_0(i, j, n, :)) .or. &
                      any(sicen(i, j, n, :) /= si_0(i, j, n, :)) .or. qsnon(i, j, n, 1) /= qs_0(i, j, n, 1)) nring = nring + 1
               end do
               if (aice0(i, j) /= a0_0(i, j) .or. aice(i, j) /= -777.0_wp .or. vice(i, j) /= -777.0_wp .or. &
                   vsno(i, j) /= -777.0_wp .or. Tc(i, j) /= -777.0_wp .or. any(qc(i, j, :) /= -777.0_wp) .or. &
                   any(sc(i, j, :) /= -777.0_wp) .or. qsc(i, j, 1) /= -777.0_wp) nring = nring + 1
               if (fr(i, j) /= fr_r(i, j) .or. fs(i, j) /= fs_r(i, j) .or. fh(i, j) /= fh_r(i, j) .or. d1(i, j) /= d1_r(i, j) .or. &
                   d2(i, j) /= d2_r(i, j) .or. dv(i, j) /= dv_r(i, j) .or. op(i, j) /= op_r(i, j) .or. ierr(i, j) /= 0_int32 .or. &
                   niter(i, j) /= -777_int32) nring = nring + 1
               cycle
            end if
            if (.not. wet(i, j) .or. iand(ierr_r(i, j), ICE_ERR_FATAL_MASK) /= 0_int32) then
               ! land (its Tf is NaN) and cells that arrive with a fatal bit: nothing touched, no new bit, no pass
               nskip = nskip + 1
               if (any(a(i, j, :) /= a_0(i, j, :)) .or. any(v(i, j, :) /= v_0(i, j, :)) .or. any(s(i, j, :) /= s_0(i, j, :)) .or. &
                 any(Tsfcn(i, j, :) /= T_0(i, j, :)) .or. any(qicen(i, j, :, :) /= q_0(i, j, :, :)) .or. aice0(i, j) /= a0_0(i, j) &
                   .or. aice(i, j) /= -777.0_wp .or. fr(i, j) /= 0.0_wp .or. fs(i, j) /= 0.0_wp .or. fh(i, j) /= 0.0_wp .or. &
                   niter(i, j) /= 0_int32 .or. ierr(i, j) /= ierr_r(i, j)) nflip = nflip + 1
               cycle
            end if
            if (ierr(i, j) /= ierr_r(i, j)) then
               nflip = nflip + 1
               write (output_unit, '(a,2i4,2i9)') "   FLIP error word, cell ", i, j, ierr_r(i, j), ierr(i, j)
            end if
            if (ierr_r(i, j) /= 0_int32) cycle
            ! a cell whose open water plus ice area is within a few ulp of 1 and which has no deformation sits ON the
            ! decision edge `asum == 1` of the quiescent early-out: the device (relaxed floating point, reassociated
            ! sums) and the host may land on either side, and the two outcomes differ by ~1e-16 in area but can
            ! book a round-off-excess flux of the clean-up that one side does not.  Such cells are not compared on
            ! pass count or emptiness; everything else about them (areas, volumes, aggregates, fluxes, diagnostics,
            ! healthy tracers) is compared against the cell scales like any other cell.
            asum_h = a0_0(i, j)
            do n = 1, NCAT
               asum_h = asum_h + a_0(i, j, n)
            end do
            ambig = abs(asum_h - 1.0_wp) <= 1.0e-15_wp .and. conv(i, j) == 0.0_wp .and. shear(i, j) == 0.0_wp
            if (ambig) nambig = nambig + 1
            if (.not. ambig) then
               if (niter(i, j) /= niter_r(i, j)) then
                  nflip = nflip + 1
                  if (nflip <= 5) write (output_unit, '(a,2i4,2i5)') "   FLIP pass count, cell ", i, j, niter_r(i, j), niter(i, j)
               end if
               do n = 1, NCAT
                  if ((a_r(i, j, n) == 0.0_wp) .neqv. (a(i, j, n) == 0.0_wp) .or. &
                      (v_r(i, j, n) == 0.0_wp) .neqv. (v(i, j, n) == 0.0_wp) .or. &
                      (s_r(i, j, n) == 0.0_wp) .neqv. (s(i, j, n) == 0.0_wp)) then
                     nflip = nflip + 1
                     if (nflip <= 5) write (output_unit, '(a,3i4,a,3es12.4,a,3es12.4)') "   FLIP zero pattern, cell/cat ", i, j, &
                        n, " host a,v,s ", a_r(i, j, n), v_r(i, j, n), s_r(i, j, n), " device ", a(i, j, n), v(i, j, n), s(i, j, n)
                  end if
               end do
            end if
            do n = 1, NCAT
               call scl("area", a_r(i, j, n), a(i, j, n), 1.0_wp, 3, i, j)
               call scl("ice volume", v_r(i, j, n), v(i, j, n), 1.0_wp, 3, i, j)
               call scl("snow volume", s_r(i, j, n), s(i, j, n), 1.0_wp, 3, i, j)
            end do
            call scl("aice", ai_r(i, j), aice(i, j), 1.0_wp, 3, i, j)
            call scl("vice", vi_r(i, j), vice(i, j), 1.0_wp, 3, i, j)
            call scl("vsno", vs_r(i, j), vsno(i, j), 1.0_wp, 3, i, j)
            call scl("aice0", a0_r(i, j), aice0(i, j), 1.0_wp, 3, i, j)
            ! booked fluxes against the cell's flux scale: the clean-up books a round-off-excess share (aice - 1)/aice
            ! of ~1e-16 on one side and nothing on the other when FMA / reassociation moves aice across 1
            mscale = p%c%rhoi*sum(v_0(i, j, :)) + p%c%rhos*sum(s_0(i, j, :))
            if (mscale > 0.0_wp) then
               call scl("fresh", fr_r(i, j), fr(i, j), mscale/dt, 4, i, j)
               call scl("fsalt", fs_r(i, j), fs(i, j), 0.01_wp*mscale/dt, 4, i, j)
               call scl("fhocn", fh_r(i, j), fh(i, j), 4.0e5_wp*mscale/dt, 4, i, j)
            else
               call scl("fresh", fr_r(i, j), fr(i, j), 1.0_wp, 4, i, j)
               call scl("fsalt", fs_r(i, j), fs(i, j), 1.0_wp, 4, i, j)
               call scl("fhocn", fh_r(i, j), fh(i, j), 1.0_wp, 4, i, j)
            end if
            ! the ridging diagnostics: sums over categories of terms through differences of exponentials, judged against
            ! the cell's total area (rates) and ice volume (volume rate) per step
            call scl("dardg1dt", d1_r(i, j), d1(i, j), max(sum(a_0(i, j, :)), tiny(1.0_wp))/dt, 4, i, j)
            call scl("dardg2dt", d2_r(i, j), d2(i, j), max(sum(a_0(i, j, :)), tiny(1.0_wp))/dt, 4, i, j)
            call scl("opening", op_r(i, j), op(i, j), max(sum(a_0(i, j, :)), 1.0_wp)/dt, 4, i, j)
            call scl("dvirdgdt", dv_r(i, j), dv(i, j), max(sum(v_0(i, j, :)), tiny(1.0_wp))/dt, 4, i, j)
            do n = 1, NCAT
               ! healthy: not drained by the ridging to below 0.1 % of its own area, ice or snow volume (a tracer is
               ! a ratio of a near-cancelling difference, conditioned by the inverse of that fraction) and not carried
               ! by a round-off residue: a category that had no snow can receive ~1e-19 m of it from a ridge through
               ! the redistribution, whose tracer (a ratio of two round-off numbers, or the stale value of an empty
               ! carrier) a decision-edge cell, which the device may leave quiescent, legitimately does not reproduce
               if (a_r(i, j, n) > 1.0e3_wp*PUNY .and. a_r(i, j, n) > 1.0e-3_wp*a_0(i, j, n) .and. &
           v_r(i, j, n) > max(1.0e-3_wp*v_0(i, j, n), 1.0e-12_wp) .and. s_r(i, j, n) > max(1.0e-3_wp*s_0(i, j, n), 1.0e-12_wp)) then
                  do it = 1, ice_ntrcr(nilyr, 1)
                     call rel(2, "tracer", ice_trcr_get(it, i, j, n, nx, ny, NCAT, nilyr, 1, T_r, q_r, si_r, qs_r), &
                              ice_trcr_get(it, i, j, n, nx, ny, NCAT, nilyr, 1, Tsfcn, qicen, sicen, qsnon), i, j)
                  end do
               end if
            end do
         end do
      end do
      write (output_unit, '(a,2i3,a,i0,a,i0,a,i0,a,i0)') "ridge device sweep (kpartic, kredist ", p%ridge%krdg_partic, &
   p%ridge%krdg_redist, "): ridged cells ", nridge, ", multi-pass ", nmulti, ", decision flips ", nflip, ", ring violations ", nring
   write (output_unit, '(a,es10.2,a,a,2i4)') "  worst healthy-tracer relative difference ", wv(2), " at ", trim(wn(2)), wi(2), wj(2)
      write (output_unit, '(a,es10.2,a,a,2i4)') "  worst absolute difference of areas, volumes, aggregates ", wv(3), " at ", trim(wn(3)), wi(3), wj(3)
    write (output_unit, '(a,es10.2,a,a,2i4)') "  worst scaled flux/diagnostic difference ", wv(4), " at ", trim(wn(4)), wi(4), wj(4)
      write (output_unit, '(a,i0,a,i0)') "  decision-edge cells (asum == 1) compared on values only: ", nambig, "; land / pre-fatal cells skipped: ", nskip
      call check(error, nring, 0, "the ring outside the owned range is untouched (all outputs)")
      if (allocated(error)) return
 call check(error, nflip, 0, "pass counts, error words, emptiness; land and fatal-bit cells untouched, no error from a NaN land Tf")
      if (allocated(error)) return
call check(error, nmulti >= 3 .and. nridge >= 500 .and. nskip >= 20, "the sweep contains ridging, multi-pass, land and fatal cells")
      if (allocated(error)) return
      call check(error, wv(2) <= TOL_INT, "intensive tracers of healthy categories agree relatively")
      if (allocated(error)) return
      call check(error, wv(3) <= 1.0e-14_wp, "areas, volumes and aggregates agree to round-off of O(1) values")
      if (allocated(error)) return
      call check(error, wv(4) <= TOL_EXT, "fluxes and ridging diagnostics agree within the cell scales")
      if (allocated(error)) return
      ! ndtd /= 1 is refused with nothing touched
      call ice_step_ridge(nx, ny, NCAT, nilyr, 1, 2, nx - 1, 2, ny - 1, p, hin_loc, dt, 2, wet, Tf2, conv, shear, a, v, s, &
                          Tsfcn, qicen, sicen, qsnon, aice, vice, vsno, aice0, Tc, qc, sc, qsc, fr, fs, fh, d1, d2, dv, op, &
                          niter, ierr, ok)
      call check(error,.not. ok, "ndtd = 2 is refused (ok = .false.)")
   contains
      subroutine scl(name, ref, got, scale_f, cls, ci, cj)
         !! Absolute difference over `scale_f`, with an explicit finiteness test before anything is accumulated.
         character(len=*), intent(in) :: name
         real(wp), intent(in) :: ref, got, scale_f
         integer, intent(in) :: cls, ci, cj
         real(wp) :: dd

         if (.not. (ieee_is_finite(ref) .and. ieee_is_finite(got))) then
            nflip = nflip + 1
            return
         end if
         dd = abs(ref - got)/scale_f
         if (dd > wv(cls)) then
            wv(cls) = dd; wn(cls) = name; wi(cls) = ci; wj(cls) = cj
         end if
      end subroutine scl

      subroutine rel(cls, name, ref, got, ci, cj)
         integer, intent(in) :: cls, ci, cj
         character(len=*), intent(in) :: name
         real(wp), intent(in) :: ref, got
         real(wp) :: dd

         if (.not. (ieee_is_finite(ref) .and. ieee_is_finite(got))) then
            nflip = nflip + 1
            return
         end if
         if (ref == got) return
         dd = abs(ref - got)/max(abs(ref), abs(got))
         if (dd > wv(cls)) then
            wv(cls) = dd; wn(cls) = name; wi(cls) = ci; wj(cls) = cj
         end if
      end subroutine rel
   end subroutine device_sweep

   subroutine throughput(error, p, kind)
      !! 512 x 512 cells x 5 categories x 7 layers, repeated sweeps of `ice_step_ridge` on re-seeded data.
      !! `kind = 1`: the random state; `kind = 2`: ~85 % ice-free cells in coherent blocks.  The clean-up alone
      !! (same state) is timed too: the difference is the ridging launch.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t), intent(in) :: p
      integer, intent(in) :: kind
      integer, parameter :: nx = 512, ny = 512, nilyr = 7, nrep = 5
      real(wp), allocatable :: a(:, :, :), v(:, :, :), s(:, :, :), Tsfcn(:, :, :), qicen(:, :, :, :), sicen(:, :, :, :)
      real(wp), allocatable :: qsnon(:, :, :, :), aice(:, :), vice(:, :), vsno(:, :), aice0(:, :), Tc(:, :), qc(:, :, :)
      real(wp), allocatable :: sc(:, :, :), qsc(:, :, :), fr(:, :), fs(:, :), fh(:, :), d1(:, :), d2(:, :), dv(:, :), op(:, :)
      real(wp), allocatable :: Tf2(:, :), conv(:, :), shear(:, :)
      real(wp), allocatable :: a_s(:, :, :), v_s(:, :, :), s_s(:, :, :), T_s(:, :, :), q_s(:, :, :, :), si_s(:, :, :, :)
      real(wp), allocatable :: qs_s(:, :, :, :), a0_s(:, :)
      integer(int32), allocatable :: ierr(:, :), niter(:, :), ierr_c(:, :)
      logical, allocatable :: wet(:, :)
      real(wp) :: hin_loc(0:NCAT), secs, secs_c
      integer(int64) :: t0, t1, rate
      integer :: k, i, j, nice, nmulti
      logical :: coherent, ok

      coherent = kind == 2
      allocate (a(nx, ny, NCAT), v(nx, ny, NCAT), s(nx, ny, NCAT), Tsfcn(nx, ny, NCAT), qicen(nx, ny, NCAT, nilyr), &
                sicen(nx, ny, NCAT, nilyr), qsnon(nx, ny, NCAT, 1), aice(nx, ny), vice(nx, ny), vsno(nx, ny), aice0(nx, ny), &
                Tc(nx, ny), qc(nx, ny, nilyr), sc(nx, ny, nilyr), qsc(nx, ny, 1), fr(nx, ny), fs(nx, ny), fh(nx, ny), &
                d1(nx, ny), d2(nx, ny), dv(nx, ny), op(nx, ny), Tf2(nx, ny), conv(nx, ny), shear(nx, ny), ierr(nx, ny), &
                niter(nx, ny), wet(nx, ny), ierr_c(nx, ny))
      allocate (a_s(nx, ny, NCAT), v_s(nx, ny, NCAT), s_s(nx, ny, NCAT), T_s(nx, ny, NCAT), q_s(nx, ny, NCAT, nilyr), &
                si_s(nx, ny, NCAT, nilyr), qs_s(nx, ny, NCAT, 1), a0_s(nx, ny))
      call fill_random(nx, ny, nilyr, merge(818181_int64, 616161_int64, coherent), p, coherent, a, v, s, Tsfcn, qicen, sicen, &
                       qsnon, aice0, conv, shear)
      a_s = a; v_s = v; s_s = s; T_s = Tsfcn; q_s = qicen; si_s = sicen; qs_s = qsnon; a0_s = aice0
      hin_loc = HIN
      Tf2 = -1.8_wp
      wet = .true.
      fr = 0.0_wp; fs = 0.0_wp; fh = 0.0_wp; d1 = 0.0_wp; d2 = 0.0_wp; dv = 0.0_wp; op = 0.0_wp
      aice = 0.0_wp; vice = 0.0_wp; vsno = 0.0_wp; Tc = 0.0_wp; qc = 0.0_wp; sc = 0.0_wp; qsc = 0.0_wp
      ierr = 0_int32; niter = 0_int32; ierr_c = 0_int32
      !$acc enter data copyin(hin_loc, wet, Tf2, conv, shear, a, v, s, Tsfcn, qicen, sicen, qsnon, aice, vice, vsno, aice0, Tc, qc, sc, qsc, fr, fs, fh, d1, d2, dv, op, niter, ierr, ierr_c, a_s, v_s, s_s, T_s, q_s, si_s, qs_s, a0_s)
      secs = 0.0_wp
      secs_c = 0.0_wp
      call system_clock(count_rate=rate)
      do k = 0, nrep
         !$acc kernels present(a, v, s, Tsfcn, qicen, sicen, qsnon, aice0, a_s, v_s, s_s, T_s, q_s, si_s, qs_s, a0_s)
         a = a_s; v = v_s; s = s_s; Tsfcn = T_s; qicen = q_s; sicen = si_s; qsnon = qs_s; aice0 = a0_s
         !$acc end kernels
         call system_clock(t0)
         call ice_step_ridge(nx, ny, NCAT, nilyr, 1, 1, nx, 1, ny, p, hin_loc, 3600.0_wp, 1, wet, Tf2, conv, shear, a, v, s, &
                             Tsfcn, qicen, sicen, qsnon, aice, vice, vsno, aice0, Tc, qc, sc, qsc, fr, fs, fh, d1, d2, dv, op, &
                             niter, ierr, ok)
         !$acc wait
         call system_clock(t1)
         if (k > 0) secs = secs + real(t1 - t0, wp)/real(rate, wp)
         !$acc kernels present(a, v, s, Tsfcn, qicen, sicen, qsnon, aice0, a_s, v_s, s_s, T_s, q_s, si_s, qs_s, a0_s)
         a = a_s; v = v_s; s = s_s; Tsfcn = T_s; qicen = q_s; sicen = si_s; qsnon = qs_s; aice0 = a0_s
         !$acc end kernels
         call system_clock(t0)
         call ice_cleanup_itd_step(nx, ny, NCAT, nilyr, 1, 1, nx, 1, ny, p, hin_loc, 3600.0_wp, Tf2, .true., .true., .true., a, v, &
                                   s, Tsfcn, qicen, sicen, qsnon, aice, vice, vsno, aice0, Tc, qc, sc, qsc, fr, fs, fh, ierr_c)
         !$acc wait
         call system_clock(t1)
         if (k > 0) secs_c = secs_c + real(t1 - t0, wp)/real(rate, wp)
      end do
      !$acc update self(ierr, niter)
      !$acc exit data delete(hin_loc, wet, Tf2, conv, shear, a, v, s, Tsfcn, qicen, sicen, qsnon, aice, vice, vsno, aice0, Tc, qc, sc, qsc, fr, fs, fh, d1, d2, dv, op, niter, ierr, ierr_c, a_s, v_s, s_s, T_s, q_s, si_s, qs_s, a0_s)
      nice = 0; nmulti = 0
      do j = 1, ny
         do i = 1, nx
            if (niter(i, j) >= 1) nice = nice + 1
            if (niter(i, j) >= 2) nmulti = nmulti + 1
         end do
      end do
      write (output_unit, '(a,i0,a,f10.3,a,es11.3,a,f10.3,a,i0,a,i0)') "ridge step 512x512x5x7, state ", kind, ": ", &
         secs/real(nrep, wp)*1.0e3_wp, " ms/sweep, ", real(nx*ny, wp)*real(nrep, wp)/secs, " cells/s (clean-up alone ", &
         secs_c/real(nrep, wp)*1.0e3_wp, " ms); ridged cells ", nice, ", multi-pass ", nmulti
      if (any(ierr /= 0_int32)) then
         do j = 1, ny
            do i = 1, nx
               if (ierr(i, j) /= 0_int32) then
                  write (output_unit, '(a,2i5,a,i0)') "   first cell with an error word: ", i, j, " ierr ", ierr(i, j)
                  exit
               end if
            end do
            if (any(ierr(:, j) /= 0_int32)) exit
         end do
      end if
      call check(error, all(ierr == 0_int32), "no error bit")
   end subroutine throughput

end module test_cice_ridge_ice
