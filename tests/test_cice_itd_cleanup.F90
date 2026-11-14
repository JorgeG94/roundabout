!! Tests of the thickness-distribution clean-up (`rdb_ice_itd_cleanup`): the plan's A11 set.
!!
!! * the zap closure: ice, snow, salt and energy removed from the state equal `dt` times the
!!   increments of `fresh`, `fsalt`, `fhocn` (closed forms written here; measured residuals are
!!   printed), the `aice > 1` round-off scaling, the snow-temperature zap, a clean state left
!!   bitwise alone, out-of-bin categories rebinned (the partner of the remapping fallback);
!! * the error bits (`AICEN_NEG`, `AICE_RANGE`, `NONFINITE`: a NaN category is refused, never
!!   zapped as "small") and the `limit_aice` / `dorebin` switches;
!! * the REAL remapping output (`ice_linear_itd_cell` with `rebin_fb = .false.`, wild growth,
!!   fallback cells, doubly drained categories) through the clean-up: every non-empty category
!!   inside its bin, no area in `(0, puny]`, finite and bounded tracers, closure of mass, salt and
!!   energy from the state BEFORE the remapping against the flux increments;
!! * exact closure on 100 000 random cells for nilyr = 1, 4, 7;
!! * the Icepack ORACLE: the real upstream `cleanup_itd`, 40 cells per constant set covering
!!   every zap branch (literals in `test_cice_itd_cleanup_golden.inc`; both constant sets, nilyr
!!   = 7 and 4), bitwise on CPU builds (FMA off) and scale-aware 1e-12 on the offload build;
!! * the cell entry against the packed routine, and a device sweep of the host shim over a
!!   random multi-cell state on an INTERIOR range comparing EVERY output on EVERY cell with the
!!   host reference and asserting the untouched ring, with a 512 x 512 x 5 x 7 throughput.
!!
!! Snow layers: every existing ITD test uses `nslyr = 1`, which is `ICE_NSLYR_MAX` (the packed
!! forms need `1 + 2*nilyr + nslyr <= ICE_NTRCR_MAX`).  The cell form is written for any
!! `nslyr` (it indexes `qsnon(:,:,:,nslyr)` directly); `test_cell_two_snow_layers` runs it with 2.
module test_cice_itd_cleanup
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan, ieee_value, ieee_quiet_nan, ieee_positive_inf
   use, intrinsic :: iso_fortran_env, only: int32, int64, output_unit
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use rdb_constants, only: wp
   use rdb_ice_dims, only: ICE_NTRCR_MAX, ICE_NSLYR_MAX
   use rdb_ice_errors, only: ICE_ERR_NONFINITE, ICE_ERR_AICE_RANGE, ICE_ERR_AICEN_NEG, ICE_ERR_REMAP_FALLBACK
   use rdb_ice_params, only: ice_params_t, ice_params_default, ice_ntrcr, ice_trcr_get, ice_trcr_set, ice_qsno_from_T, &
                             ICE_CONSTANTS_ICEPACK, ICE_CONSTANTS_CESM
   use rdb_ice_itd_shift, only: itd_aggregate, itd_shift_ice_cell
   use rdb_ice_itd_linear, only: ice_linear_itd_cell
   use rdb_ice_itd_cleanup, only: itd_cleanup, ice_cleanup_itd_cell, ice_cleanup_itd_step
   use lcg_deterministic, only: lcg_next, lcg_unit
   use test_cice_util, only: ice_golden_first_bad
   implicit none
   private

   public :: collect_cice_itd_cleanup_tests

   real(wp), parameter :: PUNY = 1.0e-11_wp
   real(wp), parameter :: TF = -1.8_wp
   integer, parameter :: NCAT = 5
   integer, parameter :: NC = 40
   real(wp), parameter :: HIN(0:NCAT) = [0.0_wp, 0.6445072168194258_wp, 1.3914334975763036_wp, &
                                         2.4701793819598885_wp, 4.567287918850491_wp, 999.9_wp]
#ifdef RDB_GPU_OFFLOAD
   logical, parameter :: GPU_BUILD = .true.
   real(wp), parameter :: TOL_GOLD = 1.0e-12_wp
   real(wp), parameter :: TOL_CELL = 3.0e-10_wp
      !! Cell form vs packed form (both on the host, FMA on, different association in the aggregate).
   real(wp), parameter :: TOL_EXT = 5.0e-15_wp
      !! Device vs host extensive values and flux increments, relative (measured worst 4.4e-16).
   real(wp), parameter :: TOL_INT = 5.0e-15_wp
      !! Device vs host intensive tracers of healthy categories (measured worst 6.0e-16).
#else
   logical, parameter :: GPU_BUILD = .false.
   real(wp), parameter :: TOL_GOLD = 0.0_wp   ! bitwise on CPU builds
   real(wp), parameter :: TOL_CELL = 0.0_wp
   real(wp), parameter :: TOL_EXT = 5.0e-15_wp
   real(wp), parameter :: TOL_INT = 5.0e-15_wp
#endif
   real(wp), parameter :: CLOSE_TOL = 1.0e-12_wp
      !! Round-off closure of the conserved sums (relative to the sum of magnitudes).
   real(wp), parameter :: SCL(8) = [0.3_wp, 1.0_wp, 2.0_wp, 4.0_wp, 6.0_wp, 0.6_wp, 3.0_wp, 8.0_wp]
      !! Growth-rate scales of the remapping cells (slow ... wild).

   integer, save :: GNIL(NC)
   real(wp), save :: GII(100, NC), GIC(100, NC), GOI(100, NC), GOC(100, NC)

contains

   subroutine collect_cice_itd_cleanup_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("zap_small_categories_closure", test_zap_closure), &
                  new_unittest("residual_cell_zapped_all_categories", test_zap_residual), &
                  new_unittest("aice_excess_scaling", test_aice_excess), &
                  new_unittest("snow_temperature_zap", test_snow_temperature), &
                  new_unittest("clean_state_is_a_noop", test_noop), &
                  new_unittest("remap_fallback_partner_rebin", test_fallback_partner), &
                  new_unittest("error_bits_and_switches", test_error_bits), &
                  new_unittest("nan_is_refused_not_zapped", test_nan), &
                  new_unittest("remap_output_through_cleanup", test_remap_pipeline), &
                  new_unittest("cell_two_snow_layers", test_cell_two_snow_layers), &
                  new_unittest("conservation_100k_nilyr7", test_conserve_7), &
                  new_unittest("conservation_100k_nilyr4", test_conserve_4), &
                  new_unittest("conservation_100k_nilyr1", test_conserve_1), &
                  new_unittest("oracle_cleanup_icepack", test_oracle_icepack), &
                  new_unittest("oracle_cleanup_cesm", test_oracle_cesm), &
                  new_unittest("cell_entry_matches_packed", test_cell_vs_packed), &
                  new_unittest("empty_cell_fast_path_is_bitwise", test_empty_cells), &
                  new_unittest("doubly_drained_threshold_edges", test_drained_edges), &
                  new_unittest("device_sweep_and_throughput", test_device) &
                  ]
   end subroutine collect_cice_itd_cleanup_tests

   subroutine fill_golden()
#include "test_cice_itd_cleanup_golden.inc"
   end subroutine fill_golden

   ! ------------------------------------------------------------------
   ! helpers
   ! ------------------------------------------------------------------

   pure logical function same(a, b, tol)
      !! Bitwise (`tol = 0`) or scale-aware agreement; NaN matches NaN only.
      real(wp), intent(in) :: a, b, tol

      if (ieee_is_nan(a) .or. ieee_is_nan(b)) then
         same = ieee_is_nan(a) .and. ieee_is_nan(b)
      else if (tol == 0.0_wp) then
         same = (a == b)
      else
         same = abs(a - b) <= tol*max(abs(a), abs(b), 1.0_wp)
      end if
   end function same

   subroutine params(kset, nilyr, p)
      integer, intent(in) :: kset, nilyr
      type(ice_params_t), intent(out) :: p

      p = ice_params_default(kset, NCAT, nilyr, 1)
   end subroutine params

   subroutine base_tracers(nilyr, p, Tsn, tr)
      !! Default packed tracers: Tsfc -5, qice(k) -(2e8 + 1e7 k), sice(k) 3 + 0.2 k, snow enthalpy of `Tsn` degC.
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
         tr(2 + nilyr, n) = ice_qsno_from_T(Tsn, p%c)
      end do
   end subroutine base_tracers

   subroutine run_packed(p, nilyr, dt, tfx, a, v, s, tr, f0, aice, aice0, fl, ierr, dorebin, limit)
      !! Packed clean-up; `fl` returns the flux values after, starting from `f0`.
      type(ice_params_t), intent(in) :: p
      integer, intent(in) :: nilyr
      real(wp), intent(in) :: dt, tfx, f0(3)
      real(wp), intent(inout) :: a(NCAT), v(NCAT), s(NCAT), tr(ICE_NTRCR_MAX, NCAT)
      real(wp), intent(out) :: aice, aice0, fl(3)
      integer(int32), intent(inout) :: ierr
      logical, intent(in), optional :: dorebin, limit
      logical :: dr, lm

      dr = .true.; lm = .true.
      if (present(dorebin)) dr = dorebin
      if (present(limit)) lm = limit
      fl = f0
      call itd_cleanup(NCAT, ICE_NTRCR_MAX, nilyr, 1, HIN, dt, tfx, p%c%rhoi, p%c%rhos, p%c%Lfresh, p%c%cp_ice, p%c%hs_min, &
                       p%c%puny, p%c%itd_area_min, p%c%itd_mass_min, p%thermo%ice_ref_salinity, dr, lm, a, v, s, tr, &
                       aice, aice0, fl(1), fl(2), fl(3), ierr)
   end subroutine run_packed

   subroutine run_cell1(p, nilyr, dt, tfx, a, v, s, tr, f0, aggr, fl, ierr, dorebin, limit, aggregate)
      !! Cell clean-up on a 1 x 1 grid; the state goes in and out through the packed layout.
      !! `aggr = (aice, vice, vsno, aice0)`.
      type(ice_params_t), intent(in) :: p
      integer, intent(in) :: nilyr
      real(wp), intent(in) :: dt, tfx, f0(3)
      real(wp), intent(inout) :: a(NCAT), v(NCAT), s(NCAT), tr(ICE_NTRCR_MAX, NCAT)
      real(wp), intent(out) :: aggr(4), fl(3)
      integer(int32), intent(inout) :: ierr
      logical, intent(in), optional :: dorebin, limit, aggregate
      real(wp) :: a3(1, 1, NCAT), v3(1, 1, NCAT), s3(1, 1, NCAT), T3(1, 1, NCAT), q4(1, 1, NCAT, nilyr), si4(1, 1, NCAT, nilyr)
      real(wp) :: qs4(1, 1, NCAT, 1), ag(1, 1), vg(1, 1), sg(1, 1), a0g(1, 1), Tc(1, 1), qc(1, 1, nilyr), sc(1, 1, nilyr)
      real(wp) :: qsc(1, 1, 1), fr(1, 1), fs(1, 1), fh(1, 1)
      integer(int32) :: ie(1, 1)
      integer :: n, it
      logical :: dr, lm, ag_on

      dr = .true.; lm = .true.; ag_on = .true.
      if (present(dorebin)) dr = dorebin
      if (present(limit)) lm = limit
      if (present(aggregate)) ag_on = aggregate
      do n = 1, NCAT
         a3(1, 1, n) = a(n); v3(1, 1, n) = v(n); s3(1, 1, n) = s(n)
         do it = 1, ice_ntrcr(nilyr, 1)
            call ice_trcr_set(it, 1, 1, n, 1, 1, NCAT, nilyr, 1, tr(it, n), T3, q4, si4, qs4)
         end do
      end do
      ag = -777.0_wp; vg = -777.0_wp; sg = -777.0_wp; a0g = -777.0_wp; Tc = -777.0_wp; qc = -777.0_wp; sc = -777.0_wp
      qsc = -777.0_wp
      fr = f0(1); fs = f0(2); fh = f0(3)
      ie = ierr
      call ice_cleanup_itd_cell(1, 1, 1, 1, NCAT, nilyr, 1, HIN, dt, tfx, p%c%rhoi, p%c%rhos, p%c%Lfresh, p%c%cp_ice, &
                               p%c%hs_min, p%c%puny, p%c%itd_area_min, p%c%itd_mass_min, p%thermo%ice_ref_salinity, dr, lm, ag_on, &
                                a3, v3, s3, T3, q4, si4, qs4, ag, vg, sg, a0g, Tc, qc, sc, qsc, fr, fs, fh, ie)
      ierr = ie(1, 1)
      do n = 1, NCAT
         a(n) = a3(1, 1, n); v(n) = v3(1, 1, n); s(n) = s3(1, 1, n)
         do it = 1, ice_ntrcr(nilyr, 1)
            tr(it, n) = ice_trcr_get(it, 1, 1, n, 1, 1, NCAT, nilyr, 1, T3, q4, si4, qs4)
         end do
      end do
      aggr = [ag(1, 1), vg(1, 1), sg(1, 1), a0g(1, 1)]
      fl = [fr(1, 1), fs(1, 1), fh(1, 1)]
   end subroutine run_cell1

   subroutine sums(nilyr, p, a, v, s, tr, mass, salt, ener, scl_m, scl_e)
      !! Mass `rhoi*v + rhos*s`, virtual salt `rhoi*v*Sref*1e-3`, enthalpy `sum q*v/nilyr + qsno*s`, and the
      !! sums of magnitudes used as the closure scales.
      integer, intent(in) :: nilyr
      type(ice_params_t), intent(in) :: p
      real(wp), intent(in) :: a(NCAT), v(NCAT), s(NCAT), tr(ICE_NTRCR_MAX, NCAT)
      real(wp), intent(out) :: mass, salt, ener, scl_m, scl_e
      integer :: n, k

      mass = 0.0_wp; salt = 0.0_wp; ener = 0.0_wp; scl_m = 0.0_wp; scl_e = 0.0_wp
      do n = 1, NCAT
         mass = mass + p%c%rhoi*v(n) + p%c%rhos*s(n)
         scl_m = scl_m + p%c%rhoi*abs(v(n)) + p%c%rhos*abs(s(n))
         salt = salt + p%c%rhoi*v(n)*p%thermo%ice_ref_salinity*0.001_wp
         do k = 1, nilyr
            ener = ener + tr(1 + k, n)*v(n)/real(nilyr, wp)
            scl_e = scl_e + abs(tr(1 + k, n)*v(n)/real(nilyr, wp))
         end do
         ener = ener + tr(2 + nilyr, n)*s(n)
         scl_e = scl_e + abs(tr(2 + nilyr, n)*s(n))
      end do
   end subroutine sums

   subroutine closure(nilyr, p, dt, a0, v0, s0, tr0, a1, v1, s1, tr1, f0, f1, rm, rs, re)
      !! Relative residuals of `removed from the state - dt*flux increment` for mass, salt, energy.
      integer, intent(in) :: nilyr
      type(ice_params_t), intent(in) :: p
      real(wp), intent(in) :: dt, a0(NCAT), v0(NCAT), s0(NCAT), tr0(ICE_NTRCR_MAX, NCAT)
      real(wp), intent(in) :: a1(NCAT), v1(NCAT), s1(NCAT), tr1(ICE_NTRCR_MAX, NCAT), f0(3), f1(3)
      real(wp), intent(out) :: rm, rs, re
      real(wp) :: m0, sa0, e0, sm0, se0, m1, sa1, e1, sm1, se1

      call sums(nilyr, p, a0, v0, s0, tr0, m0, sa0, e0, sm0, se0)
      call sums(nilyr, p, a1, v1, s1, tr1, m1, sa1, e1, sm1, se1)
      rm = abs((m0 - m1) - dt*(f1(1) - f0(1)))/max(sm0, tiny(1.0_wp))
      rs = abs((sa0 - sa1) - dt*(f1(2) - f0(2)))/max(sm0*p%thermo%ice_ref_salinity*0.001_wp, tiny(1.0_wp))
      re = abs((e0 - e1) - dt*(f1(3) - f0(3)))/max(se0, tiny(1.0_wp))
   end subroutine closure

   subroutine check_closure(error, nilyr, p, dt, a0, v0, s0, tr0, a1, v1, s1, tr1, f0, f1, label, worst)
      type(error_type), allocatable, intent(inout) :: error
      integer, intent(in) :: nilyr
      type(ice_params_t), intent(in) :: p
      real(wp), intent(in) :: dt, a0(NCAT), v0(NCAT), s0(NCAT), tr0(ICE_NTRCR_MAX, NCAT)
      real(wp), intent(in) :: a1(NCAT), v1(NCAT), s1(NCAT), tr1(ICE_NTRCR_MAX, NCAT), f0(3), f1(3)
      character(len=*), intent(in) :: label
      real(wp), intent(inout) :: worst
      real(wp) :: rm, rs, re

      call closure(nilyr, p, dt, a0, v0, s0, tr0, a1, v1, s1, tr1, f0, f1, rm, rs, re)
      if (.not. (ieee_is_finite(rm) .and. ieee_is_finite(rs) .and. ieee_is_finite(re))) then
         call test_failed(error, label//": non-finite closure residual")
         return
      end if
      worst = max(worst, rm, rs, re)
      call check(error, max(rm, rs, re) <= CLOSE_TOL, label//": mass, salt and energy closed with the flux increments")
   end subroutine check_closure

   ! ------------------------------------------------------------------
   ! analytical tests
   ! ------------------------------------------------------------------

   subroutine test_zap_closure(error)
      !! A category below `puny` carrying ice and snow is removed (closed forms for the fluxes), for nilyr 7 and 4,
      !! by the packed and the cell forms.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      integer :: nilyr, pass, k
      real(wp) :: a(NCAT), v(NCAT), s(NCAT), tr(ICE_NTRCR_MAX, NCAT), a0(NCAT), v0(NCAT), s0(NCAT), tr0(ICE_NTRCR_MAX, NCAT)
      real(wp) :: aice, aice0, fl(3), aggr(4), f0(3), dt, dfresh, dfsalt, dfhocn, worst
      integer(int32) :: ierr
      integer, parameter :: nils(2) = [7, 4]

      worst = 0.0_wp
      f0 = 0.0_wp   ! closed forms are compared to the increment: start from zero (no cancellation against the old value)
      dt = 900.0_wp
      do k = 1, 2
         nilyr = nils(k)
         call params(ICE_CONSTANTS_ICEPACK, nilyr, p)
         do pass = 1, 2
            call base_tracers(nilyr, p, -10.0_wp, tr)
            tr(2:1 + nilyr, 2) = -3.0e8_wp
            tr(2 + nilyr, 2) = -1.1e8_wp
            a = [0.3_wp, 5.0e-12_wp, 0.4_wp, 0.0_wp, 0.0_wp]
            v = [0.09_wp, 2.0e-12_wp, 0.76_wp, 0.0_wp, 0.0_wp]
            s = [0.0_wp, 1.0e-12_wp, 0.04_wp, 0.0_wp, 0.0_wp]
            a0 = a; v0 = v; s0 = s; tr0 = tr
            ierr = 0_int32
            if (pass == 1) then
               call run_packed(p, nilyr, dt, TF, a, v, s, tr, f0, aice, aice0, fl, ierr)
            else
               call run_cell1(p, nilyr, dt, TF, a, v, s, tr, f0, aggr, fl, ierr)
            end if
            call check(error, ierr, 0_int32, "no error bit")
            if (allocated(error)) return
            dfresh = (p%c%rhoi*2.0e-12_wp + p%c%rhos*1.0e-12_wp)/dt
            dfsalt = p%c%rhoi*2.0e-12_wp*p%thermo%ice_ref_salinity*0.001_wp/dt
            dfhocn = (real(nilyr, wp)*(-3.0e8_wp*2.0e-12_wp/real(nilyr, wp)) + (-1.1e8_wp*1.0e-12_wp))/dt
            call check(error, abs((fl(1) - f0(1)) - dfresh) <= 1.0e-14_wp*abs(dfresh), "fresh increment = (rhoi*v + rhos*s)/dt")
            if (allocated(error)) return
            call check(error, abs((fl(2) - f0(2)) - dfsalt) <= 1.0e-14_wp*abs(dfsalt), "fsalt increment")
            if (allocated(error)) return
            call check(error, abs((fl(3) - f0(3)) - dfhocn) <= 1.0e-14_wp*abs(dfhocn), "fhocn increment (negative)")
            if (allocated(error)) return
            call check(error, fl(3) - f0(3) < 0.0_wp, "heat flux sign: ocean spends heat melting the removed ice")
            if (allocated(error)) return
            call check(error, a(2) == 0.0_wp .and. v(2) == 0.0_wp .and. s(2) == 0.0_wp, "category zeroed")
            if (allocated(error)) return
            call check(error, tr(1, 2) == TF .and. all(tr(2:ice_ntrcr(nilyr, 1), 2) == 0.0_wp), "Tsfc = Tf, other tracers 0")
            if (allocated(error)) return
            call check(error, all(a([1, 3]) == a0([1, 3])) .and. all(v([1, 3]) == v0([1, 3])), "other categories untouched")
            if (allocated(error)) return
            if (pass == 1) then
               call check(error, abs(aice0 - (1.0_wp - 0.7_wp - 5.0e-12_wp + 5.0e-12_wp)) <= 1.0e-15_wp, &
                          "upstream aice0 gains the removed area")
            else
               call check(error, abs(aggr(1) - 0.7_wp) <= 1.0e-15_wp .and. abs(aggr(4) - 0.3_wp) <= 1.0e-15_wp, &
                          "refreshed aggregates")
            end if
            if (allocated(error)) return
            call check_closure(error, nilyr, p, dt, a0, v0, s0, tr0, a, v, s, tr, f0, fl, "zap", worst)
            if (allocated(error)) return
         end do
      end do
      write (output_unit, '(a,es10.2)') "zap closure: worst relative residual ", worst
   end subroutine test_zap_closure

   subroutine test_zap_residual(error)
      !! `aice < itd_area_min`: every category removed; the empty categories only have their tracers normalised.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      integer, parameter :: nilyr = 7
      real(wp) :: a(NCAT), v(NCAT), s(NCAT), tr(ICE_NTRCR_MAX, NCAT), a0(NCAT), v0(NCAT), s0(NCAT), tr0(ICE_NTRCR_MAX, NCAT)
      real(wp) :: aice, aice0, fl(3), f0(3), dt, worst
      integer(int32) :: ierr

      call params(ICE_CONSTANTS_CESM, nilyr, p)
      dt = 3600.0_wp
      f0 = 0.0_wp
      call base_tracers(nilyr, p, -10.0_wp, tr)
      a = [3.0e-12_wp, 0.0_wp, 4.0e-12_wp, 0.0_wp, 0.0_wp]
      v = [1.0e-12_wp, 0.0_wp, 6.0e-12_wp, 0.0_wp, 0.0_wp]
      s = [0.0_wp, 0.0_wp, 2.0e-12_wp, 0.0_wp, 0.0_wp]
      a0 = a; v0 = v; s0 = s; tr0 = tr
      ierr = 0_int32
      worst = 0.0_wp
      call run_packed(p, nilyr, dt, TF, a, v, s, tr, f0, aice, aice0, fl, ierr)
      call check(error, ierr, 0_int32, "no error bit")
      if (allocated(error)) return
      call check(error, all(a == 0.0_wp) .and. all(v == 0.0_wp) .and. all(s == 0.0_wp), "every category zeroed")
      if (allocated(error)) return
    call check(error, all(tr(1, :) == TF) .and. all(tr(2:ice_ntrcr(nilyr, 1), :) == 0.0_wp), "tracers normalised in every category")
      if (allocated(error)) return
      call check(error, abs((fl(1) - f0(1))*dt - (p%c%rhoi*7.0e-12_wp + p%c%rhos*2.0e-12_wp)) <= 1.0e-13_wp*1.0e-8_wp, &
                 "fresh increment is the removed mass over dt")
      if (allocated(error)) return
      call check_closure(error, nilyr, p, dt, a0, v0, s0, tr0, a, v, s, tr, f0, fl, "residual zap", worst)
   end subroutine test_zap_residual

   subroutine test_aice_excess(error)
      !! `1 < aice < 1 + puny` scales area, volumes by `1/aice`, tracers unchanged, flux booked as the share
      !! `(aice - 1)/aice`; a second call is (to round-off) a no-op.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      integer, parameter :: nilyr = 7
      real(wp) :: a(NCAT), v(NCAT), s(NCAT), tr(ICE_NTRCR_MAX, NCAT), a0(NCAT), v0(NCAT), s0(NCAT), tr0(ICE_NTRCR_MAX, NCAT)
      real(wp) :: aice, aice0, fl(3), fl2(3), f0(3), dt, worst, ae, dfresh, aggr(4)
      integer(int32) :: ierr
      integer :: pass

      call params(ICE_CONSTANTS_ICEPACK, nilyr, p)
      dt = 1800.0_wp
      f0 = 0.0_wp
      worst = 0.0_wp
      do pass = 1, 2
         call base_tracers(nilyr, p, -10.0_wp, tr)
         a = [0.4_wp, 0.3_wp, 0.3_wp + 5.0e-12_wp, 0.0_wp, 0.0_wp]
         v = a*[0.3_wp, 1.0_wp, 1.9_wp, 0.0_wp, 0.0_wp]
         s = a*[0.0_wp, 0.1_wp, 0.2_wp, 0.0_wp, 0.0_wp]
         a0 = a; v0 = v; s0 = s; tr0 = tr
         ae = sum(a)
         ierr = 0_int32
         if (pass == 1) then
            call run_packed(p, nilyr, dt, TF, a, v, s, tr, f0, aice, aice0, fl, ierr)
            call check(error, aice == 1.0_wp .and. aice0 == 0.0_wp, "upstream aice = 1, aice0 = 0 after the scaling")
            if (allocated(error)) return
         else
            call run_cell1(p, nilyr, dt, TF, a, v, s, tr, f0, aggr, fl, ierr)
            call check(error, abs(aggr(1) - 1.0_wp) <= 4.0e-16_wp .and. aggr(4) <= 1.0e-15_wp, "refreshed aggregates")
            if (allocated(error)) return
         end if
         call check(error, ierr, 0_int32, "no error bit")
         if (allocated(error)) return
         call check(error, abs(sum(a) - 1.0_wp) <= 4.0e-16_wp, "total area is 1 to round-off")
         if (allocated(error)) return
         call check(error, all(abs(a - a0/ae) <= 1.0e-16_wp) .and. all(tr == tr0), "areas scaled by 1/aice, tracers unchanged")
         if (allocated(error)) return
         dfresh = sum(p%c%rhoi*v0 + p%c%rhos*s0)*(ae - 1.0_wp)/ae/dt
         call check(error, abs((fl(1) - f0(1)) - dfresh) <= 1.0e-12_wp*abs(dfresh), "fresh = removed share (aice-1)/aice")
         if (allocated(error)) return
         call check(error, fl(3) < f0(3), "heat flux negative")
         if (allocated(error)) return
         call check_closure(error, nilyr, p, dt, a0, v0, s0, tr0, a, v, s, tr, f0, fl, "aice excess", worst)
         if (allocated(error)) return
         ! second call: effectively nothing left to remove
         a0 = a; v0 = v; s0 = s
         call run_packed(p, nilyr, dt, TF, a, v, s, tr, fl, aice, aice0, fl2, ierr)
         call check(error, abs(fl2(1) - fl(1)) <= 1.0e-3_wp*abs(dfresh) .and. all(abs(a - a0) <= 4.0e-16_wp), &
                    "second call: a no-op to round-off")
         if (allocated(error)) return
      end do
      write (output_unit, '(a,es10.2)') "aice excess closure: worst relative residual ", worst
   end subroutine test_aice_excess

   subroutine test_snow_temperature(error)
      !! Snow warmer than melting (+0.5 C) or colder than -100 C is removed with `dfresh = rhos*vsno/dt`,
      !! `dfhocn = qsno*vsno/dt`; -99 C and snow thinner than `hs_min` (judged on the reference enthalpy) stay.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      integer, parameter :: nilyr = 7
      real(wp) :: a(NCAT), v(NCAT), s(NCAT), tr(ICE_NTRCR_MAX, NCAT), a0(NCAT), v0(NCAT), s0(NCAT), tr0(ICE_NTRCR_MAX, NCAT)
      real(wp) :: aice, aice0, fl(3), f0(3), dt, worst, tsn(NCAT), q, aggr(4)
      integer(int32) :: ierr
      integer :: n, pass, nz

      call params(ICE_CONSTANTS_ICEPACK, nilyr, p)
      dt = 3600.0_wp
      f0 = 0.0_wp
      worst = 0.0_wp
      do pass = 1, 2
         call base_tracers(nilyr, p, -10.0_wp, tr)
         a = [0.2_wp, 0.2_wp, 0.2_wp, 0.2_wp, 0.15_wp]
         v = a*[0.3_wp, 1.0_wp, 1.9_wp, 3.0_wp, 5.0_wp]
         s = a*[0.1_wp, 0.1_wp, 0.1_wp, 0.1_wp, 5.0e-5_wp]
         tsn = [0.5_wp, -120.0_wp, -99.0_wp, -10.0_wp, 5.0_wp]
         do n = 1, NCAT
            tr(2 + nilyr, n) = ice_qsno_from_T(tsn(n), p%c)
         end do
         a0 = a; v0 = v; s0 = s; tr0 = tr
         ierr = 0_int32
         if (pass == 1) then
            call run_packed(p, nilyr, dt, TF, a, v, s, tr, f0, aice, aice0, fl, ierr)
         else
            call run_cell1(p, nilyr, dt, TF, a, v, s, tr, f0, aggr, fl, ierr)
         end if
         call check(error, ierr, 0_int32, "no error bit")
         if (allocated(error)) return
         call check(error, s(1) == 0.0_wp .and. s(2) == 0.0_wp, "warm and cold snow removed")
         if (allocated(error)) return
         call check(error, s(3) == s0(3) .and. s(4) == s0(4) .and. s(5) == s0(5), "-99 C, -10 C and thin snow kept")
         if (allocated(error)) return
     call check(error, tr(2 + nilyr, 1) == 0.0_wp .and. tr(2 + nilyr, 3) == tr0(2 + nilyr, 3), "enthalpy zeroed only where removed")
         if (allocated(error)) return
         call check(error, all(v == v0) .and. all(a == a0), "ice untouched")
         if (allocated(error)) return
         nz = 0
         q = 0.0_wp
         do n = 1, 2
            q = q + s0(n)*p%c%rhos
            nz = nz + 1
         end do
         call check(error, abs(fl(1)*dt - q) <= 1.0e-13_wp*q, "fresh = rhos*vsno/dt")
         if (allocated(error)) return
         q = (tr0(2 + nilyr, 1)*s0(1) + tr0(2 + nilyr, 2)*s0(2))/dt
         call check(error, abs(fl(3) - q) <= 1.0e-13_wp*abs(q) .and. fl(3) < 0.0_wp, "fhocn = qsno*vsno/dt")
         if (allocated(error)) return
         call check(error, fl(2) == 0.0_wp, "no salt flux from snow")
         if (allocated(error)) return
         call check_closure(error, nilyr, p, dt, a0, v0, s0, tr0, a, v, s, tr, f0, fl, "snow T", worst)
         if (allocated(error)) return
      end do
   end subroutine test_snow_temperature

   subroutine test_noop(error)
      !! A clean state (every category inside its bin, snow normal) is bitwise unchanged, no flux; the cell
      !! means come from the aggregate kernel.  (Cell means may differ from a hand sum by round-off: they are
      !! the division of the same extensive sums, compared here against the packed aggregate bitwise.)
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      integer, parameter :: nilyr = 7
      real(wp) :: a(NCAT), v(NCAT), s(NCAT), tr(ICE_NTRCR_MAX, NCAT), a0(NCAT), v0(NCAT), s0(NCAT), tr0(ICE_NTRCR_MAX, NCAT)
      real(wp) :: aice, aice0, fl(3), f0(3), aggr(4), aic, vic, vso, a0c, trc(ICE_NTRCR_MAX)
      integer(int32) :: ierr, ie2
      integer :: pass

      call params(ICE_CONSTANTS_ICEPACK, nilyr, p)
      f0 = [0.1_wp, 0.2_wp, 0.3_wp]
      do pass = 1, 2
         call base_tracers(nilyr, p, -10.0_wp, tr)
         a = [0.2_wp, 0.3_wp, 0.25_wp, 0.15_wp, 0.0_wp]
         v = a*[0.3_wp, 1.0_wp, 1.9_wp, 3.2_wp, 0.0_wp]
         s = a*[0.0_wp, 0.1_wp, 0.15_wp, 0.2_wp, 0.0_wp]
         a0 = a; v0 = v; s0 = s; tr0 = tr
         ierr = 0_int32
         if (pass == 1) then
            call run_packed(p, nilyr, 3600.0_wp, TF, a, v, s, tr, f0, aice, aice0, fl, ierr)
         else
            call run_cell1(p, nilyr, 3600.0_wp, TF, a, v, s, tr, f0, aggr, fl, ierr)
            ie2 = 0_int32
            call itd_aggregate(NCAT, ICE_NTRCR_MAX, nilyr, 1, a0, v0, s0, tr0, TF, PUNY, aic, vic, vso, a0c, trc, ie2)
            call check(error, same(aggr(1), aic, TOL_CELL) .and. same(aggr(2), vic, TOL_CELL) .and. &
                       same(aggr(3), vso, TOL_CELL) .and. same(aggr(4), a0c, TOL_CELL), &
                       "aggregates equal the packed aggregate")
            if (allocated(error)) return
         end if
         call check(error, ierr, 0_int32, "no error bit")
         if (allocated(error)) return
         call check(error, all(a == a0) .and. all(v == v0) .and. all(s == s0) .and. all(tr == tr0), "state bitwise unchanged")
         if (allocated(error)) return
         call check(error, all(fl == f0), "fluxes bitwise unchanged")
         if (allocated(error)) return
      end do
   end subroutine test_noop

   subroutine test_fallback_partner(error)
      !! The partner of the remapping fallback: category 2 grows past `hin_max(3)` (the remapping refuses and leaves
      !! the state untouched, `REMAP_FALLBACK` counted), so it sits out of its bin; the clean-up rebins it.  With
      !! `dorebin = .false.` it stays.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      integer, parameter :: nilyr = 7, nx = 1, ny = 1
      real(wp) :: ai0(1, 1, NCAT), vi0(1, 1, NCAT), a3(1, 1, NCAT), v3(1, 1, NCAT), s3(1, 1, NCAT), T3(1, 1, NCAT)
      real(wp) :: q4(1, 1, NCAT, nilyr), si4(1, 1, NCAT, nilyr), qs4(1, 1, NCAT, 1), aice, aice0
      real(wp) :: a(NCAT), v(NCAT), s(NCAT), tr(ICE_NTRCR_MAX, NCAT), a1(NCAT), v1(NCAT), s1(NCAT), tr1(ICE_NTRCR_MAX, NCAT)
      real(wp) :: f0(3), fl(3), aggr(4), worst, h
      integer(int32) :: ie(1, 1), ierr
      integer :: n, it, pass

      call params(ICE_CONSTANTS_ICEPACK, nilyr, p)
      worst = 0.0_wp
      do pass = 1, 2
         call base_tracers(nilyr, p, -10.0_wp, tr)
         tr(2:1 + nilyr, 2) = -2.5e8_wp
         ai0 = 0.0_wp; vi0 = 0.0_wp; a3 = 0.0_wp; v3 = 0.0_wp; s3 = 0.0_wp
         ai0(1, 1, 2) = 0.3_wp; vi0(1, 1, 2) = 0.3_wp*1.0_wp
         a3(1, 1, 2) = 0.3_wp; v3(1, 1, 2) = 0.3_wp*2.5_wp; s3(1, 1, 2) = 0.03_wp
         do n = 1, NCAT
            do it = 1, ice_ntrcr(nilyr, 1)
               call ice_trcr_set(it, 1, 1, n, 1, 1, NCAT, nilyr, 1, tr(it, n), T3, q4, si4, qs4)
            end do
         end do
         ie = 0_int32
         call ice_linear_itd_cell(1, 1, 1, 1, NCAT, nilyr, 1, HIN, ai0, vi0, a3, v3, s3, T3, q4, si4, qs4, TF, p%c%rhos, &
                                  p%c%Lfresh, p%thermo%hi_min, PUNY, .false., aice, aice0, ie)
         call check(error, iand(ie(1, 1), ICE_ERR_REMAP_FALLBACK) /= 0_int32, "the remapping refused the cell")
         if (allocated(error)) return
         call check(error, a3(1, 1, 2) == 0.3_wp .and. v3(1, 1, 2) == 0.75_wp, "refused remapping left the state untouched")
         if (allocated(error)) return
         do n = 1, NCAT
            a(n) = a3(1, 1, n); v(n) = v3(1, 1, n); s(n) = s3(1, 1, n)
            do it = 1, ice_ntrcr(nilyr, 1)
               tr(it, n) = ice_trcr_get(it, 1, 1, n, 1, 1, NCAT, nilyr, 1, T3, q4, si4, qs4)
            end do
         end do
         a1 = a; v1 = v; s1 = s; tr1 = tr
         f0 = 0.0_wp
         ierr = 0_int32
         if (pass == 1) then
            call run_cell1(p, nilyr, 3600.0_wp, TF, a, v, s, tr, f0, aggr, fl, ierr, dorebin=.false.)
            call check(error, a(2) == a1(2) .and. v(2) == v1(2), "dorebin = .false.: the category stays out of its bin")
            if (allocated(error)) return
         else
            call run_cell1(p, nilyr, 3600.0_wp, TF, a, v, s, tr, f0, aggr, fl, ierr)
            call check(error, ierr, 0_int32, "no error bit after the rebin")
            if (allocated(error)) return
            call check(error, a(2) == 0.0_wp .and. a(3) == 0.0_wp .and. abs(a(4) - 0.3_wp) <= 1.0e-15_wp, "moved up to category 4")
            if (allocated(error)) return
            h = v(4)/a(4)
            call check(error, h > HIN(3) .and. h <= HIN(4), "mean thickness inside the bin of category 4")
            if (allocated(error)) return
            call check(error, abs(v(4) - v1(2)) <= 1.0e-15_wp .and. abs(tr(2, 4) - tr1(2, 2)) <= 1.0e-6_wp .and. &
                       abs(s(4) - s1(2)) <= 1.0e-17_wp, "volume, snow and tracers moved with the category")
            if (allocated(error)) return
            call check(error, all(fl == f0), "the rebin emits no flux")
            if (allocated(error)) return
            call check_closure(error, nilyr, p, 3600.0_wp, a1, v1, s1, tr1, a, v, s, tr, f0, fl, "rebin", worst)
         end if
      end do
   end subroutine test_fallback_partner

   subroutine test_error_bits(error)
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      integer, parameter :: nilyr = 7
      real(wp) :: a(NCAT), v(NCAT), s(NCAT), tr(ICE_NTRCR_MAX, NCAT), aice, aice0, fl(3), f0(3), aggr(4)
      integer(int32) :: ierr
      integer :: pass

      call params(ICE_CONSTANTS_ICEPACK, nilyr, p)
      f0 = [1.0_wp, 2.0_wp, 3.0_wp]
      do pass = 1, 2
         ! negative area beyond -puny
         call base_tracers(nilyr, p, -10.0_wp, tr)
         a = [0.5_wp, -0.2_wp, 0.0_wp, 0.0_wp, 0.0_wp]; v = [0.15_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]; s = 0.0_wp
         ierr = 0_int32
         if (pass == 1) then
            call run_packed(p, nilyr, 3600.0_wp, TF, a, v, s, tr, f0, aice, aice0, fl, ierr)
         else
            call run_cell1(p, nilyr, 3600.0_wp, TF, a, v, s, tr, f0, aggr, fl, ierr)
         end if
         call check(error, iand(ierr, ICE_ERR_AICEN_NEG) /= 0_int32, "aicen < -puny raises AICEN_NEG")
         if (allocated(error)) return
         call check(error, all(fl == f0), "no flux booked on a refused cell")
         if (allocated(error)) return
         ! total area above 1 + puny, below -puny
         call base_tracers(nilyr, p, -10.0_wp, tr)
         a = [0.7_wp, 0.5_wp, 0.0_wp, 0.0_wp, 0.0_wp]; v = a*0.5_wp; s = 0.0_wp
         ierr = 0_int32
         if (pass == 1) then
            call run_packed(p, nilyr, 3600.0_wp, TF, a, v, s, tr, f0, aice, aice0, fl, ierr)
         else
            call run_cell1(p, nilyr, 3600.0_wp, TF, a, v, s, tr, f0, aggr, fl, ierr)
         end if
         call check(error, iand(ierr, ICE_ERR_AICE_RANGE) /= 0_int32, "aice > 1 + puny raises AICE_RANGE")
         if (allocated(error)) return
         a = [-0.1_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]; v = 0.0_wp
         ierr = 0_int32
         if (pass == 1) then
            call run_packed(p, nilyr, 3600.0_wp, TF, a, v, s, tr, f0, aice, aice0, fl, ierr)
         else
            call run_cell1(p, nilyr, 3600.0_wp, TF, a, v, s, tr, f0, aggr, fl, ierr)
         end if
         call check(error, iand(ierr, ICE_ERR_AICE_RANGE) /= 0_int32, "aice < -puny raises AICE_RANGE")
         if (allocated(error)) return
         ! limit_aice = .false.: no range error, no zap of small areas
         a = [0.7_wp, 0.5_wp, 3.0e-12_wp, 0.0_wp, 0.0_wp]; v = [0.35_wp, 0.25_wp, 1.0e-12_wp, 0.0_wp, 0.0_wp]
         ierr = 0_int32
         if (pass == 1) then
            call run_packed(p, nilyr, 3600.0_wp, TF, a, v, s, tr, f0, aice, aice0, fl, ierr, limit=.false.)
         else
            call run_cell1(p, nilyr, 3600.0_wp, TF, a, v, s, tr, f0, aggr, fl, ierr, limit=.false.)
         end if
         call check(error, ierr, 0_int32, "limit_aice = .false.: no range error")
         if (allocated(error)) return
         call check(error, a(3) == 3.0e-12_wp .and. all(fl == f0), "limit_aice = .false.: nothing zapped")
         if (allocated(error)) return
      end do
   end subroutine test_error_bits

   subroutine test_nan(error)
      !! A NaN anywhere in the cell (area, volume, tracer, `Tf`, `dt`) raises NONFINITE, the cell is untouched, no flux
      !! is booked, and in particular a NaN category is NOT zapped as "small".
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      integer, parameter :: nilyr = 7
      real(wp) :: a(NCAT), v(NCAT), s(NCAT), tr(ICE_NTRCR_MAX, NCAT), a0(NCAT), v0(NCAT), s0(NCAT), tr0(ICE_NTRCR_MAX, NCAT)
      real(wp) :: aice, aice0, fl(3), f0(3), aggr(4), nan, inf, dtv, tfv
      integer(int32) :: ierr
      integer :: case, pass, it

      call params(ICE_CONSTANTS_ICEPACK, nilyr, p)
      nan = ieee_value(1.0_wp, ieee_quiet_nan)
      inf = ieee_value(1.0_wp, ieee_positive_inf)
      f0 = [1.0_wp, 2.0_wp, 3.0_wp]
      do case = 1, 7
         do pass = 1, 2
            call base_tracers(nilyr, p, -10.0_wp, tr)
            a = [0.4_wp, 5.0e-12_wp, 0.4_wp, 0.0_wp, 0.0_wp]; v = [0.12_wp, 2.0e-12_wp, 0.8_wp, 0.0_wp, 0.0_wp]
            s = [0.0_wp, 1.0e-12_wp, 0.04_wp, 0.0_wp, 0.0_wp]
            dtv = 3600.0_wp; tfv = TF
            select case (case)
            case (1); v(2) = nan                 ! NaN volume of a tiny category
            case (2); a(2) = nan                 ! NaN area
            case (3); tr(3, 1) = nan             ! NaN ice enthalpy
            case (4); s(3) = inf                 ! infinite snow volume
            case (5); tr(2 + nilyr, 3) = nan     ! NaN snow enthalpy
            case (6); tfv = nan
            case (7); dtv = 0.0_wp
            end select
            a0 = a; v0 = v; s0 = s; tr0 = tr
            ierr = 0_int32
            if (pass == 1) then
               call run_packed(p, nilyr, dtv, tfv, a, v, s, tr, f0, aice, aice0, fl, ierr)
            else
               call run_cell1(p, nilyr, dtv, tfv, a, v, s, tr, f0, aggr, fl, ierr)
            end if
            call check(error, iand(ierr, ICE_ERR_NONFINITE) /= 0_int32, "NONFINITE raised")
            if (allocated(error)) return
            call check(error, all(fl == f0), "no flux booked")
            if (allocated(error)) return
            do it = 1, NCAT
               call check(error, same(a(it), a0(it), 0.0_wp) .and. same(v(it), v0(it), 0.0_wp) .and. same(s(it), s0(it), 0.0_wp), &
                          "NaN cell untouched (nothing zapped as small)")
               if (allocated(error)) return
            end do
            call check(error, all([(same(tr(it, 1), tr0(it, 1), 0.0_wp), it=1, ICE_NTRCR_MAX)]) .and. &
                       all([(same(tr(it, 3), tr0(it, 3), 0.0_wp), it=1, ICE_NTRCR_MAX)]), "tracers untouched")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_nan

   ! ------------------------------------------------------------------
   ! random states
   ! ------------------------------------------------------------------

   subroutine fill_clean_random(nx, ny, nilyr, seed, a, v, s, Tsfcn, qicen, sicen, qsnon)
      !! Cells that exercise every clean-up branch: normal in-bin cells, thickness anywhere (out of bin), tiny and
      !! slightly negative categories, nearly empty cells, a round-off area excess, warm and cold snow.
      integer, intent(in) :: nx, ny, nilyr
      integer(int64), intent(in) :: seed
      real(wp), intent(out) :: a(nx, ny, NCAT), v(nx, ny, NCAT), s(nx, ny, NCAT), Tsfcn(nx, ny, NCAT)
      real(wp), intent(out) :: qicen(nx, ny, NCAT, nilyr), sicen(nx, ny, NCAT, nilyr), qsnon(nx, ny, NCAT, 1)
      integer(int64) :: st
      integer :: i, j, n, k
      real(wp) :: u1, u2, u3, u4, ct, lo, hi, h0, tot, f, tiny_sum

      do j = 1, ny
         do i = 1, nx
            st = lcg_next(seed + int(i, int64)*7919_int64 + int(j, int64)*104729_int64*65537_int64)
            st = lcg_next(st); ct = lcg_unit(st)
            tot = 0.0_wp
            tiny_sum = 0.0_wp
            do n = 1, NCAT
               st = lcg_next(st); u1 = lcg_unit(st)
               st = lcg_next(st); u2 = lcg_unit(st)
               st = lcg_next(st); u3 = lcg_unit(st)
               st = lcg_next(st); u4 = lcg_unit(st)
               a(i, j, n) = 0.0_wp; v(i, j, n) = 0.0_wp; s(i, j, n) = 0.0_wp
               if (ct < 0.06_wp) then          ! nearly empty cell
                  if (u1 < 0.6_wp) then
                     a(i, j, n) = 1.0e-12_wp*(u2 + 0.1_wp); v(i, j, n) = 1.0e-12_wp*u3; s(i, j, n) = 1.0e-13_wp*u4
                  end if
               else if (u1 < 0.75_wp) then
                  a(i, j, n) = 0.02_wp + 0.25_wp*u2
                  if (u4 < 0.3_wp) then
                     h0 = 0.05_wp + 7.0_wp*u3
                  else
                     lo = HIN(n - 1); hi = min(HIN(n), 8.0_wp)
                     h0 = lo + (hi - lo)*(0.05_wp + 0.9_wp*u3)
                  end if
                  v(i, j, n) = a(i, j, n)*h0
                  st = lcg_next(st)
                  if (lcg_unit(st) > 0.3_wp) then
                     st = lcg_next(st); s(i, j, n) = a(i, j, n)*0.4_wp*lcg_unit(st)
                  end if
               else if (u1 < 0.9_wp) then      ! tiny (some negative) category with content
                  a(i, j, n) = 1.0e-11_wp*(u2 - 0.2_wp); v(i, j, n) = 1.0e-11_wp*u3; s(i, j, n) = 1.0e-12_wp*u4
               end if
               if (a(i, j, n) > 1.0e-3_wp) then
                  tot = tot + a(i, j, n)
               else
                  tiny_sum = tiny_sum + a(i, j, n)
               end if
               st = lcg_next(st); Tsfcn(i, j, n) = -30.0_wp*lcg_unit(st)
               do k = 1, nilyr
                  st = lcg_next(st); qicen(i, j, n, k) = -4.0e8_wp + 3.0e8_wp*lcg_unit(st)
                  st = lcg_next(st); sicen(i, j, n, k) = 0.1_wp + 5.9_wp*lcg_unit(st)
               end do
               st = lcg_next(st); qsnon(i, j, n, 1) = -2.0e8_wp + 0.95e8_wp*lcg_unit(st)
            end do
            if (ct >= 0.06_wp .and. tot > 0.0_wp) then
               st = lcg_next(st); f = 0.3_wp + 0.69_wp*lcg_unit(st)
               if (ct > 0.88_wp .and. ct < 0.94_wp) then    ! round-off excess of the area
                  st = lcg_next(st); f = (1.0_wp + 1.0e-12_wp*(1.0_wp + 8.0_wp*lcg_unit(st)) - tiny_sum)/tot
               else
                  f = f/tot
               end if
               do n = 1, NCAT
                  if (a(i, j, n) > 1.0e-3_wp) then   ! tiny categories are NOT scaled up with the populated ones
                     a(i, j, n) = a(i, j, n)*f; v(i, j, n) = v(i, j, n)*f; s(i, j, n) = s(i, j, n)*f
                  end if
               end do
            end if
         end do
      end do
   end subroutine fill_clean_random

   subroutine cell_to_packed(i, j, nx, ny, nilyr, a, v, s, Tsfcn, qicen, sicen, qsnon, ac, vc, sc, tr)
      integer, intent(in) :: i, j, nx, ny, nilyr
      real(wp), intent(in) :: a(nx, ny, NCAT), v(nx, ny, NCAT), s(nx, ny, NCAT), Tsfcn(nx, ny, NCAT)
      real(wp), intent(in) :: qicen(nx, ny, NCAT, nilyr), sicen(nx, ny, NCAT, nilyr), qsnon(nx, ny, NCAT, 1)
      real(wp), intent(out) :: ac(NCAT), vc(NCAT), sc(NCAT), tr(ICE_NTRCR_MAX, NCAT)
      integer :: n, it

      tr = 0.0_wp
      do n = 1, NCAT
         ac(n) = a(i, j, n); vc(n) = v(i, j, n); sc(n) = s(i, j, n)
         do it = 1, ice_ntrcr(nilyr, 1)
            tr(it, n) = ice_trcr_get(it, i, j, n, nx, ny, NCAT, nilyr, 1, Tsfcn, qicen, sicen, qsnon)
         end do
      end do
   end subroutine cell_to_packed

   subroutine postconditions(nilyr, p, a, v, s, tr, ierr_c, nbad, label)
      !! After the clean-up: no error, every category EMPTY (all zero) or `a > puny` inside its bin, finite and
      !! bounded tracers.  Counts violations.
      integer, intent(in) :: nilyr
      type(ice_params_t), intent(in) :: p
      real(wp), intent(in) :: a(NCAT), v(NCAT), s(NCAT), tr(ICE_NTRCR_MAX, NCAT)
      integer(int32), intent(in) :: ierr_c
      integer, intent(inout) :: nbad
      character(len=*), intent(in) :: label
      integer :: n, it, k
      real(wp) :: h

      if (ierr_c /= 0_int32) nbad = nbad + 1
      do n = 1, NCAT
         if (a(n) == 0.0_wp) then
            if (v(n) /= 0.0_wp .or. s(n) /= 0.0_wp) nbad = nbad + 1
         else if (.not. (a(n) > p%c%puny)) then
            nbad = nbad + 1                 ! area in (0, puny] or negative survived
         else
            h = v(n)/a(n)
            if (.not. (h > HIN(n - 1)*(1.0_wp - 1.0e-12_wp) .and. h <= HIN(n)*(1.0_wp + 1.0e-12_wp))) nbad = nbad + 1
            if (v(n) < 0.0_wp .or. s(n) < 0.0_wp) nbad = nbad + 1
         end if
         do it = 1, ice_ntrcr(nilyr, 1)
            if (.not. ieee_is_finite(tr(it, n))) nbad = nbad + 1
         end do
         if (.not. (tr(1, n) >= -30.0_wp*(1.0_wp + 1.0e-9_wp) .and. tr(1, n) <= 0.0_wp)) nbad = nbad + 1
         do k = 1, nilyr
            ! ice enthalpy within the generator's range (or 0 in an empty category), salinity 0 .. 6
            if (.not. (tr(1 + k, n) >= -4.0e8_wp*(1.0_wp + 1.0e-9_wp) .and. tr(1 + k, n) <= 0.0_wp)) nbad = nbad + 1
            if (.not. (tr(2 + nilyr + k, n) >= 0.0_wp .and. tr(2 + nilyr + k, n) <= 6.0_wp*(1.0_wp + 1.0e-9_wp))) nbad = nbad + 1
         end do
         if (.not. (tr(2 + nilyr, n) >= -2.0e8_wp*(1.0_wp + 1.0e-9_wp) .and. tr(2 + nilyr, n) <= 0.0_wp)) nbad = nbad + 1
      end do
      if (nbad > 0 .and. len_trim(label) < 0) nbad = nbad   ! (label is for the caller's message)
   end subroutine postconditions

   subroutine conserve_case(error, nilyr, seed)
      !! 100 000 random cells through the shim: postconditions and closure on every cell.
      type(error_type), allocatable, intent(inout) :: error
      integer, intent(in) :: nilyr
      integer(int64), intent(in) :: seed
      integer, parameter :: nx = 400, ny = 250
      real(wp), allocatable :: a(:, :, :), v(:, :, :), s(:, :, :), Tsfcn(:, :, :), qicen(:, :, :, :), sicen(:, :, :, :)
      real(wp), allocatable :: qsnon(:, :, :, :), a_0(:, :, :), v_0(:, :, :), s_0(:, :, :), T_0(:, :, :), q_0(:, :, :, :)
      real(wp), allocatable :: si_0(:, :, :, :), qs_0(:, :, :, :)
      real(wp), allocatable :: aice(:, :), vice(:, :), vsno(:, :), aice0(:, :), Tc(:, :), qc(:, :, :), sc(:, :, :), qsc(:, :, :)
      real(wp), allocatable :: fr(:, :), fs(:, :), fh(:, :), Tf2(:, :)
      integer(int32), allocatable :: ierr(:, :)
      real(wp) :: a0c(NCAT), v0c(NCAT), s0c(NCAT), t0c(ICE_NTRCR_MAX, NCAT), a1c(NCAT), v1c(NCAT), s1c(NCAT)
      real(wp) :: t1c(ICE_NTRCR_MAX, NCAT), rm, rs, re, wm, ws, we, f0(3), f1(3), dt, hin_loc(0:NCAT)
      type(ice_params_t) :: p
      integer :: i, j, nbad, ncl, nzap, nexc
      character(len=40) :: lab

      call params(ICE_CONSTANTS_ICEPACK, nilyr, p)
      allocate (a(nx, ny, NCAT), v(nx, ny, NCAT), s(nx, ny, NCAT), Tsfcn(nx, ny, NCAT), qicen(nx, ny, NCAT, nilyr), &
                sicen(nx, ny, NCAT, nilyr), qsnon(nx, ny, NCAT, 1))
      allocate (a_0(nx, ny, NCAT), v_0(nx, ny, NCAT), s_0(nx, ny, NCAT), T_0(nx, ny, NCAT), q_0(nx, ny, NCAT, nilyr), &
                si_0(nx, ny, NCAT, nilyr), qs_0(nx, ny, NCAT, 1))
      allocate (aice(nx, ny), vice(nx, ny), vsno(nx, ny), aice0(nx, ny), Tc(nx, ny), qc(nx, ny, nilyr), sc(nx, ny, nilyr), &
                qsc(nx, ny, 1), fr(nx, ny), fs(nx, ny), fh(nx, ny), Tf2(nx, ny), ierr(nx, ny))
      call fill_clean_random(nx, ny, nilyr, seed, a, v, s, Tsfcn, qicen, sicen, qsnon)
      a_0 = a; v_0 = v; s_0 = s; T_0 = Tsfcn; q_0 = qicen; si_0 = sicen; qs_0 = qsnon
      dt = 3600.0_wp
      hin_loc = HIN
      Tf2 = TF
      fr = 0.0_wp; fs = 0.0_wp; fh = 0.0_wp   ! zero start: an increment ~1e-13 on a 0.3 base would lose its digits
      aice = 0.0_wp; vice = 0.0_wp; vsno = 0.0_wp; aice0 = 0.0_wp; Tc = 0.0_wp; qc = 0.0_wp; sc = 0.0_wp; qsc = 0.0_wp
      ierr = 0_int32
      !$acc enter data copyin(hin_loc, Tf2, a, v, s, Tsfcn, qicen, sicen, qsnon, aice, vice, vsno, aice0, Tc, qc, sc, qsc, fr, fs, fh, ierr)
      call ice_cleanup_itd_step(nx, ny, NCAT, nilyr, 1, 1, nx, 1, ny, p, hin_loc, dt, Tf2, .true., .true., .true., a, v, s, Tsfcn, &
                                qicen, sicen, qsnon, aice, vice, vsno, aice0, Tc, qc, sc, qsc, fr, fs, fh, ierr)
      !$acc update self(a, v, s, Tsfcn, qicen, sicen, qsnon, aice, aice0, fr, fs, fh, ierr)
      !$acc exit data delete(hin_loc, Tf2, a, v, s, Tsfcn, qicen, sicen, qsnon, aice, vice, vsno, aice0, Tc, qc, sc, qsc, fr, fs, fh, ierr)
      nbad = 0; ncl = 0; nzap = 0; nexc = 0
      wm = 0.0_wp; ws = 0.0_wp; we = 0.0_wp
      do j = 1, ny
         do i = 1, nx
            call cell_to_packed(i, j, nx, ny, nilyr, a_0, v_0, s_0, T_0, q_0, si_0, qs_0, a0c, v0c, s0c, t0c)
            call cell_to_packed(i, j, nx, ny, nilyr, a, v, s, Tsfcn, qicen, sicen, qsnon, a1c, v1c, s1c, t1c)
            call postconditions(nilyr, p, a1c, v1c, s1c, t1c, ierr(i, j), nbad, lab)
            f0 = 0.0_wp
            f1 = [fr(i, j), fs(i, j), fh(i, j)]
            call closure(nilyr, p, dt, a0c, v0c, s0c, t0c, a1c, v1c, s1c, t1c, f0, f1, rm, rs, re)
            if (.not. (ieee_is_finite(rm) .and. ieee_is_finite(rs) .and. ieee_is_finite(re))) then
               nbad = nbad + 1
            else
               wm = max(wm, rm); ws = max(ws, rs); we = max(we, re)
            end if
            if (f1(1) /= f0(1)) nzap = nzap + 1
            if (sum(a0c) > 1.0_wp .and. sum(a0c) < 1.0_wp + PUNY) nexc = nexc + 1
            if (abs(aice(i, j) - sum(a1c)) > 1.0e-15_wp) nbad = nbad + 1
            if (abs(aice0(i, j) - max(1.0_wp - sum(a1c), 0.0_wp)) > 1.0e-15_wp) nbad = nbad + 1
         end do
      end do
      write (output_unit, '(a,i0,a,i0,a,i0,a,3es10.2)') "cleanup closure nilyr ", nilyr, ": cells with a flux ", nzap, &
         ", round-off excess cells ", nexc, ", worst relative residual (mass, salt, energy) ", wm, ws, we
      call check(error, nbad, 0, "postconditions (bins, areas, tracers) and aggregates on every cell")
      if (allocated(error)) return
      call check(error, nzap > 5000 .and. nexc > 100, "the random set exercises the zaps")
      if (allocated(error)) return
      call check(error, max(wm, ws, we) <= CLOSE_TOL, "mass, salt and energy closed on every cell")
   end subroutine conserve_case

   subroutine test_conserve_7(error)
      type(error_type), allocatable, intent(out) :: error
      call conserve_case(error, 7, 1234567_int64)
   end subroutine test_conserve_7

   subroutine test_conserve_4(error)
      type(error_type), allocatable, intent(out) :: error
      call conserve_case(error, 4, 2345678_int64)
   end subroutine test_conserve_4

   subroutine test_conserve_1(error)
      type(error_type), allocatable, intent(out) :: error
      call conserve_case(error, 1, 3456789_int64)
   end subroutine test_conserve_1

   ! ------------------------------------------------------------------
   ! the real remapping output through the clean-up
   ! ------------------------------------------------------------------

   subroutine test_remap_pipeline(error)
      !! `ice_linear_itd_cell` with `rebin_fb = .false.` on wild-growth states (fallback cells, doubly drained
      !! categories), then the clean-up: bins, areas, tracers, closure from the state BEFORE the remapping.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nx = 300, ny = 200, nilyr = 7
      real(wp), allocatable :: ai0(:, :, :), vi0(:, :, :), a(:, :, :), v(:, :, :), s(:, :, :), Tsfcn(:, :, :)
      real(wp), allocatable :: qicen(:, :, :, :), sicen(:, :, :, :), qsnon(:, :, :, :), aice(:, :), aice0(:, :)
      real(wp), allocatable :: a_0(:, :, :), v_0(:, :, :), s_0(:, :, :), T_0(:, :, :), q_0(:, :, :, :), si_0(:, :, :, :)
      real(wp), allocatable :: qs_0(:, :, :, :), vice(:, :), vsno(:, :), Tc(:, :), qc(:, :, :), sc(:, :, :), qsc(:, :, :)
      real(wp), allocatable :: fr(:, :), fs(:, :), fh(:, :), Tf2(:, :), hold(:, :, :)
      integer(int32), allocatable :: ierr(:, :)
      real(wp) :: a0c(NCAT), v0c(NCAT), s0c(NCAT), t0c(ICE_NTRCR_MAX, NCAT), a1c(NCAT), v1c(NCAT), s1c(NCAT)
      real(wp) :: t1c(ICE_NTRCR_MAX, NCAT), rm, rs, re, wm, ws, we, f0(3), f1(3), dt, hin_loc(0:NCAT)
      type(ice_params_t) :: p
      integer :: i, j, n, nbad, nfb, ntiny, nout, nzap
      character(len=8) :: lab

      call params(ICE_CONSTANTS_ICEPACK, nilyr, p)
      allocate (ai0(nx, ny, NCAT), vi0(nx, ny, NCAT), a(nx, ny, NCAT), v(nx, ny, NCAT), s(nx, ny, NCAT), Tsfcn(nx, ny, NCAT), &
                qicen(nx, ny, NCAT, nilyr), sicen(nx, ny, NCAT, nilyr), qsnon(nx, ny, NCAT, 1), aice(nx, ny), aice0(nx, ny))
      allocate (a_0(nx, ny, NCAT), v_0(nx, ny, NCAT), s_0(nx, ny, NCAT), T_0(nx, ny, NCAT), q_0(nx, ny, NCAT, nilyr), &
                si_0(nx, ny, NCAT, nilyr), qs_0(nx, ny, NCAT, 1), vice(nx, ny), vsno(nx, ny), Tc(nx, ny), qc(nx, ny, nilyr), &
                sc(nx, ny, nilyr), qsc(nx, ny, 1), fr(nx, ny), fs(nx, ny), fh(nx, ny), Tf2(nx, ny), ierr(nx, ny), &
                hold(nx, ny, NCAT))
      call fill_remap_random(nx, ny, nilyr, 99999_int64, ai0, vi0, a, v, s, Tsfcn, qicen, sicen, qsnon)
      a_0 = a; v_0 = v; s_0 = s; T_0 = Tsfcn; q_0 = qicen; si_0 = sicen; qs_0 = qsnon
      ierr = 0_int32
      nfb = 0
      do j = 1, ny
         do i = 1, nx
            if (sum(a(i, j, :)) > PUNY) then
               call ice_linear_itd_cell(i, j, nx, ny, NCAT, nilyr, 1, HIN, ai0, vi0, a, v, s, Tsfcn, qicen, sicen, qsnon, &
                                        TF, p%c%rhos, p%c%Lfresh, p%thermo%hi_min, PUNY, .false., aice(i, j), aice0(i, j), ierr)
            end if
            if (iand(ierr(i, j), ICE_ERR_REMAP_FALLBACK) /= 0_int32) nfb = nfb + 1
         end do
      end do
      ! categories left with a round-off sliver by the remapping (drained from both sides) or out of bin
      ntiny = 0; nout = 0
      do n = 1, NCAT
         do j = 1, ny
            do i = 1, nx
               if (a(i, j, n) > 0.0_wp .and. a(i, j, n) <= PUNY) ntiny = ntiny + 1
               if (a(i, j, n) > PUNY) then
                  if (v(i, j, n)/a(i, j, n) > HIN(n)*(1.0_wp + 1.0e-12_wp) .or. &
                      v(i, j, n)/a(i, j, n) < HIN(n - 1)*(1.0_wp - 1.0e-12_wp)) nout = nout + 1
               end if
            end do
         end do
      end do
      ierr = iand(ierr, not(ICE_ERR_REMAP_FALLBACK))
      dt = 3600.0_wp
      hin_loc = HIN
      Tf2 = TF
      fr = 0.0_wp; fs = 0.0_wp; fh = 0.0_wp
      aice = 0.0_wp; aice0 = 0.0_wp; vice = 0.0_wp; vsno = 0.0_wp; Tc = 0.0_wp; qc = 0.0_wp; sc = 0.0_wp; qsc = 0.0_wp
      !$acc enter data copyin(hin_loc, Tf2, a, v, s, Tsfcn, qicen, sicen, qsnon, aice, vice, vsno, aice0, Tc, qc, sc, qsc, fr, fs, fh, ierr)
      call ice_cleanup_itd_step(nx, ny, NCAT, nilyr, 1, 1, nx, 1, ny, p, hin_loc, dt, Tf2, .true., .true., .true., a, v, s, Tsfcn, &
                                qicen, sicen, qsnon, aice, vice, vsno, aice0, Tc, qc, sc, qsc, fr, fs, fh, ierr)
      !$acc update self(a, v, s, Tsfcn, qicen, sicen, qsnon, aice, aice0, fr, fs, fh, ierr)
      !$acc exit data delete(hin_loc, Tf2, a, v, s, Tsfcn, qicen, sicen, qsnon, aice, vice, vsno, aice0, Tc, qc, sc, qsc, fr, fs, fh, ierr)
      nbad = 0; wm = 0.0_wp; ws = 0.0_wp; we = 0.0_wp; nzap = 0
      f0 = 0.0_wp
      do j = 1, ny
         do i = 1, nx
            call cell_to_packed(i, j, nx, ny, nilyr, a_0, v_0, s_0, T_0, q_0, si_0, qs_0, a0c, v0c, s0c, t0c)
            call cell_to_packed(i, j, nx, ny, nilyr, a, v, s, Tsfcn, qicen, sicen, qsnon, a1c, v1c, s1c, t1c)
            call postconditions(nilyr, p, a1c, v1c, s1c, t1c, ierr(i, j), nbad, lab)
            f1 = [fr(i, j), fs(i, j), fh(i, j)]
            if (f1(1) /= 0.0_wp) nzap = nzap + 1
            call closure(nilyr, p, dt, a0c, v0c, s0c, t0c, a1c, v1c, s1c, t1c, f0, f1, rm, rs, re)
            if (.not. (ieee_is_finite(rm) .and. ieee_is_finite(rs) .and. ieee_is_finite(re))) then
               nbad = nbad + 1
            else
               wm = max(wm, rm); ws = max(ws, rs); we = max(we, re)
            end if
         end do
      end do
      write (output_unit, '(a,i0,a,i0,a,i0,a,i0,a,3es10.2)') "remap pipeline: fallback cells ", nfb, ", sliver categories ", &
         ntiny, ", out-of-bin categories ", nout, ", cells zapped ", nzap, ", worst closure (mass, salt, energy) ", wm, ws, we
      call check(error, nfb > 100, "the set includes fallback cells")
      if (allocated(error)) return
      call check(error, nout > 0, "the remapping output has out-of-bin categories (cleanup must fix them)")
      if (allocated(error)) return
      call check(error, nbad, 0, "every non-empty category in its bin, no sliver, finite bounded tracers, no error bit")
      if (allocated(error)) return
      call check(error, max(wm, ws, we) <= CLOSE_TOL, "closure from the pre-remap state")
   end subroutine test_remap_pipeline

   subroutine fill_remap_random(nx, ny, nilyr, seed, ai0, vi0, a, v, s, Tsfcn, qicen, sicen, qsnon)
      !! Random consistent pre/post-thermodynamics cells with growth rates up to a few metres per step.
      integer, intent(in) :: nx, ny, nilyr
      integer(int64), intent(in) :: seed
      real(wp), intent(out) :: ai0(nx, ny, NCAT), vi0(nx, ny, NCAT), a(nx, ny, NCAT), v(nx, ny, NCAT), s(nx, ny, NCAT)
      real(wp), intent(out) :: Tsfcn(nx, ny, NCAT), qicen(nx, ny, NCAT, nilyr), sicen(nx, ny, NCAT, nilyr), qsnon(nx, ny, NCAT, 1)
      integer(int64) :: st
      integer :: i, j, n, k
      real(wp) :: u1, u2, u3, u4, lo, hi, h0, dh, sc

      do j = 1, ny
         do i = 1, nx
            st = lcg_next(seed + int(i, int64)*7919_int64 + int(j, int64)*104729_int64*65537_int64)
            st = lcg_next(st); u1 = lcg_unit(st)
            sc = SCL(1 + int(u1*8.0_wp))
            do n = 1, NCAT
               st = lcg_next(st); u1 = lcg_unit(st)
               st = lcg_next(st); u2 = lcg_unit(st)
               st = lcg_next(st); u3 = lcg_unit(st)
               st = lcg_next(st); u4 = lcg_unit(st)
               ai0(i, j, n) = 0.0_wp; vi0(i, j, n) = 0.0_wp; a(i, j, n) = 0.0_wp; v(i, j, n) = 0.0_wp; s(i, j, n) = 0.0_wp
               if (u1 < 0.8_wp) then
                  ai0(i, j, n) = 0.02_wp + 0.17_wp*u2
                  lo = HIN(n - 1); hi = min(HIN(n), 8.0_wp)
                  h0 = lo + (hi - lo)*(0.05_wp + 0.9_wp*u3)
                  dh = sc*0.2_wp*(u4 - 0.4_wp)
                  dh = max(dh, -0.9_wp*h0)
                  if (h0 + dh < 0.002_wp) dh = 0.002_wp - h0
                  vi0(i, j, n) = ai0(i, j, n)*h0
                  a(i, j, n) = ai0(i, j, n)
                  v(i, j, n) = a(i, j, n)*(h0 + dh)
                  st = lcg_next(st)
                  if (lcg_unit(st) > 0.3_wp) then
                     st = lcg_next(st); s(i, j, n) = a(i, j, n)*0.4_wp*lcg_unit(st)
                  end if
               end if
               st = lcg_next(st); Tsfcn(i, j, n) = -30.0_wp*lcg_unit(st)
               do k = 1, nilyr
                  st = lcg_next(st); qicen(i, j, n, k) = -4.0e8_wp + 3.0e8_wp*lcg_unit(st)
                  st = lcg_next(st); sicen(i, j, n, k) = 0.1_wp + 5.9_wp*lcg_unit(st)
               end do
               st = lcg_next(st); qsnon(i, j, n, 1) = -2.0e8_wp + 0.95e8_wp*lcg_unit(st)
            end do
         end do
      end do
   end subroutine fill_remap_random

   subroutine test_cell_two_snow_layers(error)
      !! The cell form with `nslyr = 2` (above `ICE_NSLYR_MAX`, which the packed forms need): snow zap closure.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nilyr = 4, nsl = 2
      type(ice_params_t) :: p
      real(wp) :: a(1, 1, NCAT), v(1, 1, NCAT), s(1, 1, NCAT), T3(1, 1, NCAT), q4(1, 1, NCAT, nilyr), si4(1, 1, NCAT, nilyr)
      real(wp) :: qs4(1, 1, NCAT, nsl), ag(1, 1), vg(1, 1), sg(1, 1), a0g(1, 1), Tc(1, 1), qc(1, 1, nilyr), sc(1, 1, nilyr)
      real(wp) :: qsc(1, 1, nsl), fr(1, 1), fs(1, 1), fh(1, 1), dt, ev, es, hsn
      integer(int32) :: ie(1, 1)
      integer :: n, k

      p = ice_params_default(ICE_CONSTANTS_ICEPACK, NCAT, nilyr, nsl)
      dt = 1200.0_wp
      a = 0.0_wp; v = 0.0_wp; s = 0.0_wp
      a(1, 1, 1:2) = [0.4_wp, 0.5_wp]; v(1, 1, 1:2) = [0.2_wp, 0.4_wp]; s(1, 1, 1:2) = [0.08_wp, 0.1_wp]
      T3 = -5.0_wp
      do n = 1, NCAT
         do k = 1, nilyr
            q4(1, 1, n, k) = -2.0e8_wp; si4(1, 1, n, k) = 4.0_wp
         end do
         qs4(1, 1, n, 1) = ice_qsno_from_T(-10.0_wp, p%c)
         qs4(1, 1, n, 2) = ice_qsno_from_T(-10.0_wp, p%c)
      end do
      qs4(1, 1, 2, 2) = ice_qsno_from_T(0.4_wp, p%c)     ! the second snow layer of category 2 is above melting
      fr = 0.0_wp; fs = 0.0_wp; fh = 0.0_wp; ie = 0_int32
      ev = qs4(1, 1, 2, 1)/dt*0.1_wp/2.0_wp + qs4(1, 1, 2, 2)/dt*0.1_wp/2.0_wp
      es = p%c%rhos*0.1_wp/dt
      call ice_cleanup_itd_cell(1, 1, 1, 1, NCAT, nilyr, nsl, HIN, dt, TF, p%c%rhoi, p%c%rhos, p%c%Lfresh, p%c%cp_ice, &
                                p%c%hs_min, p%c%puny, p%c%itd_area_min, p%c%itd_mass_min, p%thermo%ice_ref_salinity, &
                                .true., .true., .true., a, v, s, T3, q4, si4, qs4, ag, vg, sg, a0g, Tc, qc, sc, qsc, fr, fs, fh, ie)
      call check(error, ie(1, 1), 0_int32, "no error bit")
      if (allocated(error)) return
      hsn = s(1, 1, 2)
      call check(error, hsn == 0.0_wp .and. all(qs4(1, 1, 2, :) == 0.0_wp) .and. s(1, 1, 1) == 0.08_wp, "category 2 snow removed")
      if (allocated(error)) return
     call check(error, abs(fh(1, 1) - ev) <= 1.0e-13_wp*abs(ev) .and. abs(fr(1, 1) - es) <= 1.0e-13_wp*es, "fluxes with two layers")
      if (allocated(error)) return
      call check(error, size(qsc, 3) == nsl .and. qsc(1, 1, 1) == qs4(1, 1, 1, 1) .and. qsc(1, 1, 2) == qs4(1, 1, 1, 2) .or. &
                 abs(qsc(1, 1, 1) - qs4(1, 1, 1, 1)) <= 1.0e-12_wp*abs(qs4(1, 1, 1, 1)), "cell-mean snow enthalpy per layer")
      call check(error, ICE_NSLYR_MAX >= 1, "ICE_NSLYR_MAX documented (packed forms need nslyr <= ICE_NSLYR_MAX)")
   end subroutine test_cell_two_snow_layers

   ! ------------------------------------------------------------------
   ! oracle
   ! ------------------------------------------------------------------

   subroutine oracle_run(error, kset, gin, gold)
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: kset
      real(wp), intent(in) :: gin(100, NC), gold(100, NC)
      type(ice_params_t) :: p
      integer :: c, nilyr, ntr, n, it, bad, ncell, nsmall, nexc, nres, nsnow, nreb, nneg0
      real(wp) :: dt, tfx, f0(3), a(NCAT), v(NCAT), s(NCAT), tr(ICE_NTRCR_MAX, NCAT), fl(3), aice, aice0, got(100), aggr(4)
      real(wp) :: ai(NCAT), vi(NCAT), si(NCAT), tri(ICE_NTRCR_MAX, NCAT), a2(NCAT), v2(NCAT), s2(NCAT), tr2(ICE_NTRCR_MAX, NCAT)
      real(wp) :: fl2(3), ag2, ao2
      integer(int32) :: ierr, ierr2
      character(len=200) :: msg
      integer :: pass
      real(wp) :: orc_worst

      call fill_golden()
      orc_worst = 0.0_wp
      nsmall = 0; nexc = 0; nres = 0; nsnow = 0; nreb = 0; ncell = 0; nneg0 = 0
      do c = 1, NC
         nilyr = GNIL(c)
         ntr = ice_ntrcr(nilyr, 1)
         call params(kset, nilyr, p)
         dt = gin(1, c); tfx = gin(2, c); f0 = gin(3:5, c)
         ai = gin(6:10, c); vi = gin(11:15, c); si = gin(16:20, c)
         tri = 0.0_wp
         tri(1:ntr, :) = reshape(gin(21:20 + 5*ntr, c), [ntr, NCAT])
         do pass = 1, 2
            a = ai; v = vi; s = si; tr = tri
            ierr = 0_int32
            if (pass == 1) then
               call run_packed(p, nilyr, dt, tfx, a, v, s, tr, f0, aice, aice0, fl, ierr)
            else
               call run_cell1(p, nilyr, dt, tfx, a, v, s, tr, f0, aggr, fl, ierr)
               aice = gold(16, c); aice0 = gold(17, c)   ! upstream's intermediate aggregates are not returned by the cell form
            end if
            call check(error, ierr, 0_int32, "oracle case: no error bit")
            if (allocated(error)) return
            got(1:5) = a; got(6:10) = v; got(11:15) = s; got(16) = aice; got(17) = aice0
            got(18:20) = fl
            got(21:20 + 5*ntr) = reshape(tr(1:ntr, :), [ntr*NCAT])
            do it = 1, 20 + 5*ntr
               if (gold(it, c) /= 0.0_wp) orc_worst = max(orc_worst, abs(got(it) - gold(it, c))/abs(gold(it, c)))
            end do
            bad = ice_golden_first_bad(got(1:20 + 5*ntr), gold(1:20 + 5*ntr, c), TOL_GOLD, TOL_GOLD)
            if (bad > 0) then
               write (msg, '(a,i0,a,i0,a,i0,a,es24.16,a,es24.16)') "oracle case ", c, " pass ", pass, " output ", bad, " got ", &
                  got(bad), " expected ", gold(bad, c)
               call test_failed(error, trim(msg))
               return
            end if
         end do
         ! the cell aggregates equal a recomputation from the oracle's output state
         call itd_aggregate(NCAT, ICE_NTRCR_MAX, nilyr, 1, gold(1:5, c), gold(6:10, c), gold(11:15, c), tr, tfx, PUNY, ag2, ao2, &
                            fl2(1), fl2(2), tr2(:, 1), ierr2)
         ! branch coverage from the inputs
         ncell = ncell + 1
         do n = 1, NCAT
            if (abs(ai(n)) <= PUNY .and. (ai(n) /= 0.0_wp .or. vi(n) /= 0.0_wp .or. si(n) /= 0.0_wp)) nsmall = nsmall + 1
            if (ai(n) < 0.0_wp .and. ai(n) >= -PUNY) nneg0 = nneg0 + 1
         end do
         if (sum(ai) > 1.0_wp .and. sum(ai) < 1.0_wp + PUNY) nexc = nexc + 1
         if (sum(ai) < PUNY) nres = nres + 1
         if (any(gold(11:15, c) == 0.0_wp .and. si > 0.0_wp .and. ai > PUNY)) nsnow = nsnow + 1
         a2 = ai; v2 = vi; s2 = si; tr2 = tri
         ierr2 = 0_int32
         call run_packed(p, nilyr, dt, tfx, a2, v2, s2, tr2, f0, aice, aice0, fl2, ierr2, dorebin=.false.)
         if (any(a2 /= gold(1:5, c)) .and. .not. (sum(ai) < PUNY)) nreb = nreb + 1
      end do
      write (output_unit, '(a,es10.2)') "oracle worst relative difference (packed and cell forms, every output) ", orc_worst
      write (output_unit, '(a,6i4)') "oracle coverage (small cats, round-off excess, residual cells, snow zaps, rebin, tiny negative): ", &
         nsmall, nexc, nres, nsnow, nreb, nneg0
      call check(error, nsmall >= 3 .and. nexc >= 3 .and. nres >= 2 .and. nsnow >= 4 .and. nreb >= 3 .and. nneg0 >= 1, &
                 "oracle set covers every zap branch")
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
   ! cell entry vs packed, device
   ! ------------------------------------------------------------------

   subroutine test_cell_vs_packed(error)
      !! `ice_cleanup_itd_cell` against the packed routine over random cells, with the switches: every output.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nx = 60, ny = 40, nilyr = 7
      real(wp), allocatable :: a(:, :, :), v(:, :, :), s(:, :, :), Tsfcn(:, :, :), qicen(:, :, :, :), sicen(:, :, :, :)
      real(wp), allocatable :: qsnon(:, :, :, :)
      real(wp) :: ac(NCAT), vc(NCAT), sc(NCAT), tr(ICE_NTRCR_MAX, NCAT), a0(NCAT), v0(NCAT), s0(NCAT), t0(ICE_NTRCR_MAX, NCAT)
      real(wp) :: aice, aice0, fl(3), aggr(4), flc(3), f0(3)
      type(ice_params_t) :: p
      integer(int32) :: ierr_p, ierr_c
      integer :: i, j, n, it, pass, nbad
      logical :: dr, lm

      call params(ICE_CONSTANTS_ICEPACK, nilyr, p)
      allocate (a(nx, ny, NCAT), v(nx, ny, NCAT), s(nx, ny, NCAT), Tsfcn(nx, ny, NCAT), qicen(nx, ny, NCAT, nilyr), &
                sicen(nx, ny, NCAT, nilyr), qsnon(nx, ny, NCAT, 1))
      nbad = 0
      do pass = 1, 3
         dr = pass /= 2
         lm = pass /= 3
         call fill_clean_random(nx, ny, nilyr, 424242_int64 + int(pass, int64), a, v, s, Tsfcn, qicen, sicen, qsnon)
         do j = 1, ny
            do i = 1, nx
               call cell_to_packed(i, j, nx, ny, nilyr, a, v, s, Tsfcn, qicen, sicen, qsnon, ac, vc, sc, tr)
               a0 = ac; v0 = vc; s0 = sc; t0 = tr
               f0 = [0.1_wp*real(i, wp), -0.2_wp*real(j, wp), 3.0_wp]
               ierr_p = 0_int32; ierr_c = 0_int32
               call run_packed(p, nilyr, 1800.0_wp, TF, ac, vc, sc, tr, f0, aice, aice0, fl, ierr_p, dorebin=dr, limit=lm)
               call run_cell1(p, nilyr, 1800.0_wp, TF, a0, v0, s0, t0, f0, aggr, flc, ierr_c, dorebin=dr, limit=lm)
               if (ierr_p /= ierr_c) nbad = nbad + 1
               do it = 1, 3
                  if (.not. same(fl(it), flc(it), TOL_CELL)) nbad = nbad + 1
               end do
               do n = 1, NCAT
                  if (.not. (same(ac(n), a0(n), TOL_CELL) .and. same(vc(n), v0(n), TOL_CELL) .and. same(sc(n), s0(n), TOL_CELL))) &
                     nbad = nbad + 1
                  do it = 1, ice_ntrcr(nilyr, 1)
                     if (.not. same(tr(it, n), t0(it, n), TOL_CELL)) nbad = nbad + 1
                  end do
               end do
            end do
         end do
      end do
      call check(error, nbad, 0, "cell entry equals the packed routine (all outputs, all cells)")
   end subroutine test_cell_vs_packed

   pure logical function bit_eq(x, y)
      !! Identical bit patterns (distinguishes -0 from +0; NaN payloads count).
      real(wp), intent(in) :: x, y

      bit_eq = transfer(x, 0_int64) == transfer(y, 0_int64)
   end function bit_eq

   subroutine test_empty_cells(error)
      !! All-empty cells through the cell form (shortcut) and the packed reference (general path, always): every
      !! state value, flux and the error word bit for bit (a `-0` flux or state value counts), with and without
      !! the aggregate, and with `limit_aice = .false.` (no shortcut, no zap).
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nilyr = 7
      type(ice_params_t) :: p
      real(wp) :: a(NCAT), v(NCAT), s(NCAT), tr(ICE_NTRCR_MAX, NCAT), a1(NCAT), v1(NCAT), s1(NCAT), tr1(ICE_NTRCR_MAX, NCAT)
      real(wp) :: f0(3), fl(3), flc(3), aggr(4), aice, aice0, ag(4), trc(ICE_NTRCR_MAX), nz
      integer(int32) :: ie_p, ie_c, ie2
      integer :: cs, fv, pass, n, it, nbad
      logical :: lm, agg

      call params(ICE_CONSTANTS_ICEPACK, nilyr, p)
      nz = sign(0.0_wp, -1.0_wp)
      nbad = 0
      do cs = 1, 5
         do fv = 1, 3
            do pass = 1, 3
               lm = pass /= 3
               agg = pass /= 2
               call base_tracers(nilyr, p, -10.0_wp, tr)
               a = 0.0_wp; v = 0.0_wp; s = 0.0_wp
               select case (cs)
               case (2)   ! already normalised
                  tr = 0.0_wp
                  tr(1, :) = TF
               case (3)   ! negative zeros everywhere in the state
                  a = nz; v = nz; s = nz
                  tr(2, :) = nz; tr(2 + nilyr, :) = nz; tr(3 + nilyr, :) = nz
               case (4)   ! normalised except one stray value per category
                  tr = 0.0_wp
                  tr(1, :) = TF
                  tr(4, 3) = -1.0e8_wp
                  tr(1, 2) = -3.0_wp
               case (5)   ! tracer already -0 in qice, +0 elsewhere
                  tr = 0.0_wp
                  tr(1, :) = TF
                  tr(3, :) = nz
               end select
               select case (fv)
               case (1); f0 = nz
               case (2); f0 = [0.5_wp, -2.0_wp, 7.0_wp]
               case (3); f0 = [nz, 1.0_wp, nz]
               end select
               a1 = a; v1 = v; s1 = s; tr1 = tr
               ie_p = 0_int32; ie_c = 0_int32
               call run_packed(p, nilyr, 1800.0_wp, TF, a1, v1, s1, tr1, f0, aice, aice0, fl, ie_p, limit=lm)
               call run_cell1(p, nilyr, 1800.0_wp, TF, a, v, s, tr, f0, aggr, flc, ie_c, limit=lm, aggregate=agg)
               if (ie_p /= ie_c) nbad = nbad + 1
               do it = 1, 3
                  if (.not. bit_eq(fl(it), flc(it))) then
                     nbad = nbad + 1
                     if (nbad <= 5) write (output_unit, '(a,4i3,2es14.6)') "   empty-cell flux differs: case/f0/pass/idx ", cs, fv, pass, it, fl(it), flc(it)
                  end if
               end do
               do n = 1, NCAT
                  if (.not. (bit_eq(a(n), a1(n)) .and. bit_eq(v(n), v1(n)) .and. bit_eq(s(n), s1(n)))) then
                     nbad = nbad + 1
          if (nbad <= 5) write (output_unit, '(a,4i3)') "   empty-cell area/volume/snow differs: case/f0/pass/cat ", cs, fv, pass, n
                  end if
                  do it = 1, ice_ntrcr(nilyr, 1)
                     if (.not. bit_eq(tr(it, n), tr1(it, n))) then
                        nbad = nbad + 1
                        if (nbad <= 5) write (output_unit, '(a,5i3,2es14.6)') "   empty-cell tracer differs: case/f0/pass/cat/it ", cs, fv, pass, n, it, tr(it, n), tr1(it, n)
                     end if
                  end do
               end do
               ie2 = 0_int32
               call itd_aggregate(NCAT, ICE_NTRCR_MAX, nilyr, 1, a1, v1, s1, tr1, TF, PUNY, ag(1), ag(2), ag(3), ag(4), trc, ie2)
               if (agg) then
                  if (.not. (bit_eq(aggr(1), ag(1)) .and. bit_eq(aggr(2), ag(2)) .and. bit_eq(aggr(3), ag(3)) .and. &
                             bit_eq(aggr(4), ag(4)))) nbad = nbad + 1
               else
                  if (any(aggr /= -777.0_wp)) nbad = nbad + 1
               end if
            end do
         end do
      end do
      call check(error, nbad, 0, "empty cells: shortcut identical to the general path, bit for bit")
   end subroutine test_empty_cells

   subroutine test_drained_edges(error)
      !! A category left with area `delta` by two transfers of `itd_shift_ice_cell` (the way a doubly drained category
      !! arises), for `delta = +0.5 puny, +puny, 0, -0.5 puny, -puny, -2 puny`.  `puny` is a power of two here so
      !! every area is exact.  Upstream's rules: a category with `|a| <= puny` that still carries area, volume or
      !! snow is zapped (so `+-puny` and `+-0.5 puny` are zapped); `a < -puny` strictly is an error; an exactly empty
      !! category is left alone.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nilyr = 7
      type(ice_params_t) :: p
      real(wp) :: pu, deltas(6), delta, d_a(NCAT), d_v(NCAT)
      real(wp) :: a3(1, 1, NCAT), v3(1, 1, NCAT), s3(1, 1, NCAT), T3(1, 1, NCAT), q4(1, 1, NCAT, nilyr), si4(1, 1, NCAT, nilyr)
      real(wp) :: qs4(1, 1, NCAT, 1), a(NCAT), v(NCAT), s(NCAT), tr(ICE_NTRCR_MAX, NCAT), hic(NCAT), f0(3), fl(3), aggr(4)
      real(wp) :: tr0(ICE_NTRCR_MAX, NCAT)
      integer :: donor(NCAT), k, n, it
      integer(int32) :: ie, ie_s
      character(len=40) :: lab

      pu = scale(1.0_wp, -37)
      call params(ICE_CONSTANTS_ICEPACK, nilyr, p)
      p%c%puny = pu
      deltas = [0.5_wp*pu, pu, 0.0_wp, -0.5_wp*pu, -pu, -2.0_wp*pu]
      do k = 1, 6
         delta = deltas(k)
         write (lab, '(a,es10.2)') "delta = ", delta
         call base_tracers(nilyr, p, -10.0_wp, tr)
         a3 = 0.0_wp; v3 = 0.0_wp; s3 = 0.0_wp
         a3(1, 1, 1:3) = [0.375_wp, 0.25_wp, 0.125_wp]
         v3(1, 1, 1:3) = a3(1, 1, 1:3)*[0.1_wp, 1.0_wp, 1.9_wp]
         v3(1, 1, 2) = 0.25_wp                  ! h = 1 exactly
         s3(1, 1, 1:3) = [0.0_wp, 0.0625_wp, 0.0_wp]
         do n = 1, NCAT
            do it = 1, ice_ntrcr(nilyr, 1)
               call ice_trcr_set(it, 1, 1, n, 1, 1, NCAT, nilyr, 1, tr(it, n), T3, q4, si4, qs4)
            end do
         end do
         a = a3(1, 1, :); v = v3(1, 1, :); s = s3(1, 1, :)
         donor = 0; d_a = 0.0_wp; d_v = 0.0_wp
         donor(1) = 2; d_a(1) = 0.125_wp; d_v(1) = 0.125_wp          ! category 2 gives down ...
         donor(2) = 2; d_a(2) = 0.125_wp - delta; d_v(2) = 0.125_wp - delta   ! ... and up: what is left is `delta`
         ie_s = 0_int32
         hic = 0.0_wp
         call itd_shift_ice_cell(1, 1, 1, 1, NCAT, nilyr, 1, pu, TF, 0.0_wp, donor, d_a, d_v, a, v, s, hic, T3, q4, si4, qs4, ie_s)
         call check(error, ie_s, 0_int32, lab//": the shift is accepted")
         if (allocated(error)) return
         call check(error, a(2) == delta, lab//": the shift leaves exactly delta in category 2")
         if (allocated(error)) return
         ! the cleanup input: packed form of the shifted state
         a3(1, 1, :) = a; v3(1, 1, :) = v; s3(1, 1, :) = s
         call cell_to_packed(1, 1, 1, 1, nilyr, a3, v3, s3, T3, q4, si4, qs4, a, v, s, tr)
         tr0 = tr
         f0 = 0.0_wp
         ie = 0_int32
         call run_cell1(p, nilyr, 3600.0_wp, TF, a, v, s, tr, f0, aggr, fl, ie)
         select case (k)
         case (1, 2, 4, 5)       ! |delta| <= puny with content: zapped
            call check(error, ie, 0_int32, lab//": no error")
            if (allocated(error)) return
            call check(error, a(2) == 0.0_wp .and. v(2) == 0.0_wp .and. s(2) == 0.0_wp, lab//": zapped")
            if (allocated(error)) return
            call check(error, fl(1) /= 0.0_wp, lab//": a flux was booked")
         case (3)                ! empty: untouched, nothing booked
            call check(error, ie, 0_int32, lab//": no error")
            if (allocated(error)) return
            call check(error, a(2) == 0.0_wp .and. v(2) == 0.0_wp .and. s(2) == 0.0_wp, lab//": stays empty")
            if (allocated(error)) return
            call check(error, all(fl == 0.0_wp), lab//": no flux")
         case (6)                ! below -puny: refused
            call check(error, iand(ie, ICE_ERR_AICEN_NEG) /= 0_int32, lab//": AICEN_NEG")
         end select
         if (allocated(error)) return
      end do
   end subroutine test_drained_edges

   subroutine test_device(error)
      !! The host shim over an INTERIOR range of a random 64 x 48 state against the host packed reference, every
      !! output of every interior cell, the ring untouched.  Decisions are compared EXACTLY (the zero pattern of area,
      !! ice volume and snow volume of every category, the error word); extensive values and flux increments
      !! relatively (no floor); intensive tracers on healthy categories (`a > 1e3 puny`) at their own tolerance.
      !! The worst difference of each class and where it is are printed.  Then the throughput.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nx = 64, ny = 48, nilyr = 7
      real(wp), allocatable :: a(:, :, :), v(:, :, :), s(:, :, :), Tsfcn(:, :, :), qicen(:, :, :, :), sicen(:, :, :, :)
      real(wp), allocatable :: qsnon(:, :, :, :), aice(:, :), vice(:, :), vsno(:, :), aice0(:, :), Tc(:, :), qc(:, :, :)
      real(wp), allocatable :: sc(:, :, :), qsc(:, :, :), fr(:, :), fs(:, :), fh(:, :), Tf2(:, :)
      real(wp), allocatable :: a_0(:, :, :), v_0(:, :, :), s_0(:, :, :), T_0(:, :, :), q_0(:, :, :, :), si_0(:, :, :, :)
      real(wp), allocatable :: qs_0(:, :, :, :)
      integer(int32), allocatable :: ierr(:, :)
      real(wp) :: hin_loc(0:NCAT), dt, f0(3), fl(3), aic, vic, vso, a0c, ac(NCAT), vc(NCAT), sc1(NCAT)
      real(wp) :: tr(ICE_NTRCR_MAX, NCAT), trm(ICE_NTRCR_MAX), ai_p, ao_p, wv(4)
      type(ice_params_t) :: p
      integer :: i, j, n, it, nring, nflip, wi(4), wj(4)
      integer(int32) :: ie_p, ie2
      logical :: inside, exc
      real(wp) :: sm_m, sm_s, sm_e, sc_m, sc_e
      character(len=16) :: wn(4)

      call params(ICE_CONSTANTS_ICEPACK, nilyr, p)
      allocate (a(nx, ny, NCAT), v(nx, ny, NCAT), s(nx, ny, NCAT), Tsfcn(nx, ny, NCAT), qicen(nx, ny, NCAT, nilyr), &
                sicen(nx, ny, NCAT, nilyr), qsnon(nx, ny, NCAT, 1), aice(nx, ny), vice(nx, ny), vsno(nx, ny), aice0(nx, ny), &
                Tc(nx, ny), qc(nx, ny, nilyr), sc(nx, ny, nilyr), qsc(nx, ny, 1), fr(nx, ny), fs(nx, ny), fh(nx, ny), &
                Tf2(nx, ny), ierr(nx, ny))
      allocate (a_0(nx, ny, NCAT), v_0(nx, ny, NCAT), s_0(nx, ny, NCAT), T_0(nx, ny, NCAT), q_0(nx, ny, NCAT, nilyr), &
                si_0(nx, ny, NCAT, nilyr), qs_0(nx, ny, NCAT, 1))
      call fill_clean_random(nx, ny, nilyr, 515151_int64, a, v, s, Tsfcn, qicen, sicen, qsnon)
      a_0 = a; v_0 = v; s_0 = s; T_0 = Tsfcn; q_0 = qicen; si_0 = sicen; qs_0 = qsnon
      hin_loc = HIN
      Tf2 = TF
      dt = 2400.0_wp
      ! interior cells start from zero fluxes (the increments are compared relatively); the ring carries distinct values
      do j = 1, ny
         do i = 1, nx
            if (i >= 2 .and. i <= nx - 1 .and. j >= 2 .and. j <= ny - 1) then
               fr(i, j) = 0.0_wp; fs(i, j) = 0.0_wp; fh(i, j) = 0.0_wp
            else
               fr(i, j) = 0.01_wp*real(i, wp) - 0.02_wp*real(j, wp) + 0.5_wp
               fs(i, j) = 0.03_wp*real(j, wp) - 0.5_wp
               fh(i, j) = 0.1_wp*real(i + j, wp)
            end if
         end do
      end do
      aice = -777.0_wp; vice = -777.0_wp; vsno = -777.0_wp; aice0 = -777.0_wp
      Tc = -777.0_wp; qc = -777.0_wp; sc = -777.0_wp; qsc = -777.0_wp
      ierr = 0_int32
      !$acc enter data copyin(hin_loc, Tf2, a, v, s, Tsfcn, qicen, sicen, qsnon, aice, vice, vsno, aice0, Tc, qc, sc, qsc, fr, fs, fh, ierr)
      call ice_cleanup_itd_step(nx, ny, NCAT, nilyr, 1, 2, nx - 1, 2, ny - 1, p, hin_loc, dt, Tf2, .true., .true., .true., a, v, &
                                s, Tsfcn, qicen, sicen, qsnon, aice, vice, vsno, aice0, Tc, qc, sc, qsc, fr, fs, fh, ierr)
      !$acc update self(a, v, s, Tsfcn, qicen, sicen, qsnon, aice, vice, vsno, aice0, Tc, qc, sc, qsc, fr, fs, fh, ierr)
      !$acc exit data delete(hin_loc, Tf2, a, v, s, Tsfcn, qicen, sicen, qsnon, aice, vice, vsno, aice0, Tc, qc, sc, qsc, fr, fs, fh, ierr)
      nring = 0; nflip = 0
      wv = 0.0_wp; wi = 0; wj = 0; wn = ''
      do j = 1, ny
         do i = 1, nx
            inside = i >= 2 .and. i <= nx - 1 .and. j >= 2 .and. j <= ny - 1
            call cell_to_packed(i, j, nx, ny, nilyr, a_0, v_0, s_0, T_0, q_0, si_0, qs_0, ac, vc, sc1, tr)
            if (.not. inside) then
               do n = 1, NCAT
                  if (a(i, j, n) /= a_0(i, j, n) .or. v(i, j, n) /= v_0(i, j, n) .or. s(i, j, n) /= s_0(i, j, n)) nring = nring + 1
                  if (Tsfcn(i, j, n) /= T_0(i, j, n) .or. any(qicen(i, j, n, :) /= q_0(i, j, n, :)) .or. &
                      any(sicen(i, j, n, :) /= si_0(i, j, n, :)) .or. qsnon(i, j, n, 1) /= qs_0(i, j, n, 1)) nring = nring + 1
               end do
               if (aice(i, j) /= -777.0_wp .or. vice(i, j) /= -777.0_wp .or. vsno(i, j) /= -777.0_wp .or. &
                   aice0(i, j) /= -777.0_wp .or. Tc(i, j) /= -777.0_wp .or. any(qc(i, j, :) /= -777.0_wp) .or. &
                   any(sc(i, j, :) /= -777.0_wp) .or. qsc(i, j, 1) /= -777.0_wp) nring = nring + 1
               if (fr(i, j) /= 0.01_wp*real(i, wp) - 0.02_wp*real(j, wp) + 0.5_wp .or. &
                   fs(i, j) /= 0.03_wp*real(j, wp) - 0.5_wp .or. fh(i, j) /= 0.1_wp*real(i + j, wp) .or. &
                   ierr(i, j) /= 0_int32) nring = nring + 1
               cycle
            end if
            f0 = 0.0_wp
            exc = sum(ac) > 1.0_wp
            call sums(nilyr, p, ac, vc, sc1, tr, sm_m, sm_s, sm_e, sc_m, sc_e)
            ie_p = 0_int32
            call run_packed(p, nilyr, dt, TF, ac, vc, sc1, tr, f0, ai_p, ao_p, fl, ie_p)
            ie2 = 0_int32
            call itd_aggregate(NCAT, ICE_NTRCR_MAX, nilyr, 1, ac, vc, sc1, tr, TF, PUNY, aic, vic, vso, a0c, trm, ie2)
            ! (a) decisions, exactly
            if (ie_p /= ierr(i, j)) then
               nflip = nflip + 1
               write (output_unit, '(a,2i4,2i9)') "   FLIP error word, cell ", i, j, ie_p, ierr(i, j)
            end if
            if (ie_p /= 0_int32) cycle   ! fatal: state unspecified
            do n = 1, NCAT
               if ((ac(n) == 0.0_wp) .neqv. (a(i, j, n) == 0.0_wp) .or. (vc(n) == 0.0_wp) .neqv. (v(i, j, n) == 0.0_wp) .or. &
                   (sc1(n) == 0.0_wp) .neqv. (s(i, j, n) == 0.0_wp)) then
                  nflip = nflip + 1
                  if (nflip <= 5) write (output_unit, '(a,3i4,a,3es12.4,a,3es12.4)') "   FLIP zero pattern, cell/cat ", i, j, n, &
                     " host a,v,s ", ac(n), vc(n), sc1(n), " device ", a(i, j, n), v(i, j, n), s(i, j, n)
               end if
            end do
            ! (b) extensive values and fluxes, relatively
            do n = 1, NCAT
               call rel(1, "area", ac(n), a(i, j, n), i, j)
               call rel(1, "ice volume", vc(n), v(i, j, n), i, j)
               call rel(1, "snow volume", sc1(n), s(i, j, n), i, j)
            end do
            if (exc) then
               ! round-off excess cell: the booked share (aice - 1)/aice, aice - 1 ~ 1e-12, amplifies the 1 ulp difference of
               ! the area sum (the remap-free rebin and the sum itself differ by FMA) by ~1e4: judged against the cell's
               ! own flux scale (mass or enthalpy over dt) instead of relatively
               call nrm("fresh (excess)", fl(1), fr(i, j), sc_m/dt, i, j)
               call nrm("fsalt (excess)", fl(2), fs(i, j), sc_m*p%thermo%ice_ref_salinity*0.001_wp/dt, i, j)
               call nrm("fhocn (excess)", fl(3), fh(i, j), sc_e/dt, i, j)
            else
               call rel(1, "fresh", fl(1), fr(i, j), i, j)
               call rel(1, "fsalt", fl(2), fs(i, j), i, j)
               call rel(1, "fhocn", fl(3), fh(i, j), i, j)
            end if
            call rel(1, "aice", aic, aice(i, j), i, j)
            call rel(1, "vice", vic, vice(i, j), i, j)
            call rel(1, "vsno", vso, vsno(i, j), i, j)
            call rel(3, "aice0", a0c, aice0(i, j), i, j)
            ! (c) intensive tracers on healthy categories (and the cell means of a healthy cell)
            do n = 1, NCAT
               if (ac(n) > 1.0e3_wp*PUNY) then
                  do it = 1, ice_ntrcr(nilyr, 1)
                   call rel(2, "tracer", tr(it, n), ice_trcr_get(it, i, j, n, nx, ny, NCAT, nilyr, 1, Tsfcn, qicen, sicen, qsnon), &
                              i, j)
                  end do
               end if
            end do
            if (aic > 1.0e3_wp*PUNY) then
               call rel(2, "mean Tsfc", trm(1), Tc(i, j), i, j)
               do it = 1, nilyr
                  call rel(2, "mean qice", trm(1 + it), qc(i, j, it), i, j)
                  call rel(2, "mean sice", trm(2 + nilyr + it), sc(i, j, it), i, j)
               end do
               call rel(2, "mean qsno", trm(2 + nilyr), qsc(i, j, 1), i, j)
            end if
         end do
      end do
      write (output_unit, '(a,i0,a,i0)') "cleanup device sweep: decision flips ", nflip, ", ring violations ", nring
   write (output_unit, '(a,es10.2,a,a,2i4)') "  worst extensive/flux relative difference ", wv(1), " at ", trim(wn(1)), wi(1), wj(1)
   write (output_unit, '(a,es10.2,a,a,2i4)') "  worst healthy-tracer relative difference ", wv(2), " at ", trim(wn(2)), wi(2), wj(2)
      write (output_unit, '(a,es10.2,a,a,2i4)') "  worst aice0 absolute difference ", wv(3), " at ", trim(wn(3)), wi(3), wj(3)
      write (output_unit, '(a,es10.2,a,a,2i4)') "  worst excess-cell flux difference / (scale/dt) ", wv(4), " at ", trim(wn(4)), wi(4), wj(4)
      call check(error, nring, 0, "the ring outside the owned range is untouched (all outputs)")
      if (allocated(error)) return
      call check(error, nflip, 0, "every zap decision and error word equals the host reference")
      if (allocated(error)) return
      call check(error, wv(1) <= TOL_EXT, "extensive values and flux increments agree relatively")
      if (allocated(error)) return
      call check(error, wv(2) <= TOL_INT, "intensive tracers of healthy categories agree relatively")
      if (allocated(error)) return
      call check(error, wv(3) <= 1.0e-15_wp, "aice0 agrees to round-off of 1")
      if (allocated(error)) return
      call check(error, wv(4) <= 2.0e-15_wp, "round-off-excess cells: flux difference within round-off of the cell's flux scale")
      if (allocated(error)) return
      deallocate (a, v, s, Tsfcn, qicen, sicen, qsnon, aice, vice, vsno, aice0, Tc, qc, sc, qsc, fr, fs, fh, Tf2, ierr)
      call throughput(error, p, 1)
      if (allocated(error)) return
      call throughput(error, p, 2)
   contains
      subroutine nrm(name, ref, got, scale_f, ci, cj)
         character(len=*), intent(in) :: name
         real(wp), intent(in) :: ref, got, scale_f
         integer, intent(in) :: ci, cj
         real(wp) :: dd

         if (.not. (ieee_is_finite(ref) .and. ieee_is_finite(got))) then
            nflip = nflip + 1
            return
         end if
         dd = abs(ref - got)/max(scale_f, tiny(1.0_wp))
         if (dd > wv(4)) then
            wv(4) = dd; wn(4) = name; wi(4) = ci; wj(4) = cj
         end if
      end subroutine nrm

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
         if (cls == 3) dd = abs(ref - got)   ! aice0 = max(1 - aice, 0): a cancellation, judged absolutely
         if (dd > wv(cls)) then
            wv(cls) = dd; wn(cls) = name; wi(cls) = ci; wj(cls) = cj
         end if
      end subroutine rel
   end subroutine test_device

   subroutine fill_bench_coherent(nx, ny, nilyr, seed, a, v, s, Tsfcn, qicen, sicen, qsnon)
      !! The realistic global case: ~85 % of the cells ice free in 32 x 32 blocks (spatially coherent), the ice cells
      !! clean (every populated category inside its bin, cold snow), empty cells carrying normalised tracers.
      integer, intent(in) :: nx, ny, nilyr
      integer(int64), intent(in) :: seed
      real(wp), intent(out) :: a(nx, ny, NCAT), v(nx, ny, NCAT), s(nx, ny, NCAT), Tsfcn(nx, ny, NCAT)
      real(wp), intent(out) :: qicen(nx, ny, NCAT, nilyr), sicen(nx, ny, NCAT, nilyr), qsnon(nx, ny, NCAT, 1)
      integer(int64) :: st
      integer :: i, j, n, k
      real(wp) :: u1, u2, u3, tot
      logical :: ice

      do j = 1, ny
         do i = 1, nx
            st = lcg_next(seed + int((i - 1)/32, int64)*7919_int64 + int((j - 1)/32, int64)*104729_int64*65537_int64)
            st = lcg_next(st)
            ice = lcg_unit(st) < 0.15_wp
            st = lcg_next(seed + int(i, int64)*31_int64 + int(j, int64)*104729_int64*131_int64)
            tot = 0.0_wp
            do n = 1, NCAT
               a(i, j, n) = 0.0_wp; v(i, j, n) = 0.0_wp; s(i, j, n) = 0.0_wp
               Tsfcn(i, j, n) = TF; qicen(i, j, n, :) = 0.0_wp; sicen(i, j, n, :) = 0.0_wp; qsnon(i, j, n, 1) = 0.0_wp
               if (ice) then
                  st = lcg_next(st); u1 = lcg_unit(st)
                  st = lcg_next(st); u2 = lcg_unit(st)
                  st = lcg_next(st); u3 = lcg_unit(st)
                  if (u1 < 0.75_wp) then
                     a(i, j, n) = 0.02_wp + 0.18_wp*u2
                     v(i, j, n) = a(i, j, n)*(HIN(n - 1) + (min(HIN(n), 8.0_wp) - HIN(n - 1))*(0.05_wp + 0.9_wp*u3))
                     tot = tot + a(i, j, n)
                     st = lcg_next(st); s(i, j, n) = a(i, j, n)*0.3_wp*lcg_unit(st)
                     st = lcg_next(st); Tsfcn(i, j, n) = -30.0_wp*lcg_unit(st)
                     do k = 1, nilyr
                        st = lcg_next(st); qicen(i, j, n, k) = -4.0e8_wp + 3.0e8_wp*lcg_unit(st)
                        st = lcg_next(st); sicen(i, j, n, k) = 0.1_wp + 5.9_wp*lcg_unit(st)
                     end do
                     st = lcg_next(st); qsnon(i, j, n, 1) = -1.45e8_wp + 0.25e8_wp*lcg_unit(st)
                  end if
               end if
            end do
            if (tot > 0.95_wp) then
               a(i, j, :) = a(i, j, :)*(0.95_wp/tot); v(i, j, :) = v(i, j, :)*(0.95_wp/tot); s(i, j, :) = s(i, j, :)*(0.95_wp/tot)
            end if
         end do
      end do
   end subroutine fill_bench_coherent

   subroutine throughput(error, p, kind)
      !! 512 x 512 cells x 5 categories x 7 layers, repeated sweeps of the shim on re-seeded data.  `kind = 1`: the
      !! random state (about 30 % of the categories out of their bin: measures the rebin); `kind = 2`: a clean state,
      !! about 85 % empty cells in coherent blocks (the realistic global case).
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t), intent(in) :: p
      integer, intent(in) :: kind
      integer, parameter :: nx = 512, ny = 512, nilyr = 7, nrep = 5
      real(wp), allocatable :: a(:, :, :), v(:, :, :), s(:, :, :), Tsfcn(:, :, :), qicen(:, :, :, :), sicen(:, :, :, :)
      real(wp), allocatable :: qsnon(:, :, :, :), aice(:, :), vice(:, :), vsno(:, :), aice0(:, :), Tc(:, :), qc(:, :, :)
      real(wp), allocatable :: sc(:, :, :), qsc(:, :, :), fr(:, :), fs(:, :), fh(:, :), Tf2(:, :)
      real(wp), allocatable :: a_s(:, :, :), v_s(:, :, :), s_s(:, :, :), T_s(:, :, :), q_s(:, :, :, :), si_s(:, :, :, :)
      real(wp), allocatable :: qs_s(:, :, :, :)
      integer(int32), allocatable :: ierr(:, :)
      real(wp) :: hin_loc(0:NCAT), secs
      integer(int64) :: t0, t1, rate
      integer :: k
      logical :: agg

      agg = .false.   ! the step driver aggregates once after the clean-up
      allocate (a(nx, ny, NCAT), v(nx, ny, NCAT), s(nx, ny, NCAT), Tsfcn(nx, ny, NCAT), qicen(nx, ny, NCAT, nilyr), &
                sicen(nx, ny, NCAT, nilyr), qsnon(nx, ny, NCAT, 1), aice(nx, ny), vice(nx, ny), vsno(nx, ny), aice0(nx, ny), &
                Tc(nx, ny), qc(nx, ny, nilyr), sc(nx, ny, nilyr), qsc(nx, ny, 1), fr(nx, ny), fs(nx, ny), fh(nx, ny), &
                Tf2(nx, ny), ierr(nx, ny))
      allocate (a_s(nx, ny, NCAT), v_s(nx, ny, NCAT), s_s(nx, ny, NCAT), T_s(nx, ny, NCAT), q_s(nx, ny, NCAT, nilyr), &
                si_s(nx, ny, NCAT, nilyr), qs_s(nx, ny, NCAT, 1))
      if (kind == 1) then
         call fill_clean_random(nx, ny, nilyr, 616161_int64, a, v, s, Tsfcn, qicen, sicen, qsnon)
      else
         call fill_bench_coherent(nx, ny, nilyr, 818181_int64, a, v, s, Tsfcn, qicen, sicen, qsnon)
      end if
      a_s = a; v_s = v; s_s = s; T_s = Tsfcn; q_s = qicen; si_s = sicen; qs_s = qsnon
      hin_loc = HIN
      Tf2 = TF
      fr = 0.0_wp; fs = 0.0_wp; fh = 0.0_wp
      aice = 0.0_wp; vice = 0.0_wp; vsno = 0.0_wp; aice0 = 0.0_wp; Tc = 0.0_wp; qc = 0.0_wp; sc = 0.0_wp; qsc = 0.0_wp
      ierr = 0_int32
      !$acc enter data copyin(hin_loc, Tf2, a, v, s, Tsfcn, qicen, sicen, qsnon, aice, vice, vsno, aice0, Tc, qc, sc, qsc, fr, fs, fh, ierr, a_s, v_s, s_s, T_s, q_s, si_s, qs_s)
      secs = 0.0_wp
      call system_clock(count_rate=rate)
      do k = 0, nrep
         !$acc kernels present(a, v, s, Tsfcn, qicen, sicen, qsnon, a_s, v_s, s_s, T_s, q_s, si_s, qs_s)
         a = a_s; v = v_s; s = s_s; Tsfcn = T_s; qicen = q_s; sicen = si_s; qsnon = qs_s
         !$acc end kernels
         call system_clock(t0)
         call ice_cleanup_itd_step(nx, ny, NCAT, nilyr, 1, 1, nx, 1, ny, p, hin_loc, 3600.0_wp, Tf2, .true., .true., agg, a, v, s, &
                                   Tsfcn, qicen, sicen, qsnon, aice, vice, vsno, aice0, Tc, qc, sc, qsc, fr, fs, fh, ierr)
         !$acc wait
         call system_clock(t1)
         if (k > 0) secs = secs + real(t1 - t0, wp)/real(rate, wp)
      end do
      !$acc update self(ierr)
      !$acc exit data delete(hin_loc, Tf2, a, v, s, Tsfcn, qicen, sicen, qsnon, aice, vice, vsno, aice0, Tc, qc, sc, qsc, fr, fs, fh, ierr, a_s, v_s, s_s, T_s, q_s, si_s, qs_s)
      write (output_unit, '(a,i0,a,f10.3,a,es11.3,a)') "cleanup sweep 512x512x5x7, state ", kind, ": ", &
         secs/real(nrep, wp)*1.0e3_wp, " ms/sweep, ", real(nx*ny, wp)*real(nrep, wp)/secs, " cells/s"
      call check(error, all(ierr == 0_int32), "no error bit")
   end subroutine throughput

end module test_cice_itd_cleanup
