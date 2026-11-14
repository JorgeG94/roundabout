!! Tests of the linear remapping of the thickness distribution (`rdb_ice_itd_linear`), the
!! category shift and the rebin fallback (`rdb_ice_itd_shift`): the plan's A10 set.
!!
!! * a two-category known answer worked by hand (area, volume, snow and tracers moved across
!!   the upper bound of a category that grew uniformly), identity when nothing changed, the
!!   fallback to the rebin, the `hi_min` floor, the melt-area loss of category 1;
!! * each error bit (`ITD_SHIFT_RANGE`, `ITD_HBNEW_ORDER`, `REMAP_FALLBACK`, `NONFINITE`) and the
!!   tolerance trims of the shift;
!! * exact closure on 100 000 random cells for ncat = 5 and nilyr = 1, 4, 7: area (only the two
!!   documented open-water losses remove any), ice volume, snow volume, ice and snow energy, salt;
!! * the Icepack ORACLE: the real upstream `linear_itd`, 27 cells per constant set (literals in
!!   `test_cice_itd_linear_golden.inc`; both constant sets, nilyr = 7 and 4), bitwise on CPU
!!   builds (FMA off) and scale-aware 1e-12 on the offload build;
!! * the cell entry against the packed routine, and a device sweep of the host shim over a
!!   random multi-cell state comparing EVERY output on EVERY cell with the host reference, with
!!   a 512 x 512 x 5 x 7 throughput figure.
!!
!! `ICE_NSLYR_MAX` is 1, so every case here has `nslyr = 1`; a multi-snow-layer case needs that
!! limit lifted first (the kernels loop `1..nslyr` through the accessors and are written for it).
!!
!! Expected values come from closed forms written here or from the oracle literals, never from
!! the port.  The oracle regenerates its tracer inputs with the LCG of `gen_tracers` (same draw
!! order as the oracle driver), checked through a stored checksum.
module test_cice_itd_linear
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
   use, intrinsic :: iso_fortran_env, only: int32, int64, output_unit
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use rdb_constants, only: wp
   use rdb_ice_dims, only: ICE_NTRCR_MAX
   use rdb_ice_errors, only: ICE_ERR_NONFINITE, ICE_ERR_ITD_HBNEW_ORDER, ICE_ERR_ITD_SHIFT_RANGE, &
                             ICE_ERR_REMAP_FALLBACK
   use rdb_ice_params, only: ice_params_t, ice_params_default, ice_ntrcr, ice_trcr_get, &
                             ICE_CONSTANTS_ICEPACK, ICE_CONSTANTS_CESM
   use rdb_ice_itd_shift, only: itd_shift_ice, itd_shift_ice_cell, itd_rebin, itd_rebin_cell
   use rdb_ice_itd_linear, only: itd_linear_plan, itd_linear_itd, ice_linear_itd_cell, ice_linear_itd_step
   use lcg_deterministic, only: lcg_next, lcg_unit
   use test_cice_util, only: ice_golden_first_bad
   implicit none
   private

   public :: collect_cice_itd_linear_tests

   real(wp), parameter :: PUNY = 1.0e-11_wp
   real(wp), parameter :: TF = -1.8_wp
   integer, parameter :: NCAT = 5
   integer, parameter :: NC = 27
   integer, parameter :: NOUT = 97
   real(wp), parameter :: HIN(0:NCAT) = [0.0_wp, 0.6445072168194258_wp, 1.3914334975763036_wp, &
                                         2.4701793819598885_wp, 4.567287918850491_wp, 999.9_wp]
      !! Category bounds (kcatbound = 0, ncat = 5, kitd = 1; the D5.1 oracle table).
   real(wp), parameter :: SCL(16) = [0.1_wp, 0.3_wp, 0.7_wp, 1.5_wp, 3.0_wp, 0.5_wp, 0.2_wp, 1.0_wp, 2.0_wp, &
                                     0.4_wp, 0.15_wp, 1.2_wp, 2.5_wp, 0.6_wp, 0.25_wp, 1.8_wp]
      !! Growth-rate scales of the random cells (slow ... wild).
#ifdef RDB_GPU_OFFLOAD
   logical, parameter :: GPU_BUILD = .true.
   real(wp), parameter :: TOL_GOLD = 1.0e-12_wp
      !! Oracle and cell-vs-packed tolerance: FMA is on for host code on the offload build.
   real(wp), parameter :: TOL_DEV = 1.0e-12_wp
   real(wp), parameter :: TOL_TRC = 1.0e-10_wp
      !! Device tracers of a nearly drained category (area <= puny or volume < 1e-6 of the cell) divide
      !! a small difference of two large extensive products (conditioning ~1e4).  MEASURED on the
      !! V100 (see the GPU-build line printed by the test): area/volume/aggregates 3.3e-16, healthy
      !! tracers 5.5e-14, drained tracers 8.4e-12 (plain remap), 1.3e-12 (with fallback), rebin food 0; the GPU-build oracle
      !! worst (linear_itd and rebin, both constant sets) is 6.6e-16 scaled.
#else
   logical, parameter :: GPU_BUILD = .false.
   real(wp), parameter :: TOL_GOLD = 0.0_wp   ! bitwise on CPU builds
   real(wp), parameter :: TOL_DEV = 1.0e-13_wp
   real(wp), parameter :: TOL_TRC = 1.0e-13_wp
#endif
   real(wp), parameter :: CLOSE_TOL = 1.0e-12_wp
      !! Round-off closure of the conserved sums (relative to the sum of magnitudes).

   ! oracle literals
   integer, save :: GNIL(NC)
   real(wp), save :: GCHK(NC), GIN(25, NC), GOI(NOUT, NC), GOC(NOUT, NC), GRI(NOUT, NC), GRC(NOUT, NC)

contains

   subroutine collect_cice_itd_linear_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("known_answer_two_categories", test_known_answer), &
                  new_unittest("identity_when_nothing_changed", test_identity), &
                  new_unittest("cat1_melt_area_loss", test_cat1_melt), &
                  new_unittest("hi_min_floor", test_hi_min), &
                  new_unittest("cat1_zero_volume_melted_away", test_cat1_zero_volume), &
                  new_unittest("fallback_then_rebin", test_fallback), &
                  new_unittest("shift_range_and_trims", test_shift_range), &
                  new_unittest("error_bits_order_and_nonfinite", test_error_bits), &
                  new_unittest("rebin_conserves_and_bins", test_rebin_random), &
                  new_unittest("conservation_100k_nilyr7", test_conserve_7), &
                  new_unittest("conservation_100k_nilyr4", test_conserve_4), &
                  new_unittest("conservation_100k_nilyr1", test_conserve_1), &
                  new_unittest("oracle_linear_itd_icepack", test_oracle_icepack), &
                  new_unittest("oracle_linear_itd_cesm", test_oracle_cesm), &
                  new_unittest("cell_entry_matches_packed", test_cell_vs_packed), &
                  new_unittest("device_sweep_and_throughput", test_device) &
                  ]
   end subroutine collect_cice_itd_linear_tests

   subroutine fill_golden()
#include "test_cice_itd_linear_golden.inc"
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

   subroutine lcg_tracers(c, nilyr, ntr, trcrn)
      !! Tracer inputs of oracle case `c` (the draw order of the oracle driver's `gen_tracers`).
      integer, intent(in) :: c, nilyr, ntr
      real(wp), intent(out) :: trcrn(ntr, NCAT)
      integer(int64) :: st
      integer :: n, k

      st = lcg_next(int(c, int64)*104729_int64 + 17_int64)
      do n = 1, NCAT
         st = lcg_next(st); trcrn(1, n) = -30.0_wp*lcg_unit(st)
         do k = 1, nilyr
            st = lcg_next(st); trcrn(1 + k, n) = -4.0e8_wp + 3.0e8_wp*lcg_unit(st)
         end do
         st = lcg_next(st); trcrn(2 + nilyr, n) = -1.2e8_wp + 0.3e8_wp*lcg_unit(st)
         do k = 1, nilyr
            st = lcg_next(st); trcrn(2 + nilyr + k, n) = 6.0_wp*lcg_unit(st)
         end do
      end do
   end subroutine lcg_tracers

   subroutine call_packed(kset, nilyr, ai0, vi0, a, v, s, trcrn, rebin_fb, aice, aice0, ierr, hi_min_in)
      !! Packed `linear_itd` of one column with the constants of `kset` (`trcrn(ntr, ncat)`
      !! is embedded in the `ICE_NTRCR_MAX`-leading local the routine expects).
      integer, intent(in) :: kset, nilyr
      real(wp), intent(in) :: ai0(NCAT), vi0(NCAT)
      real(wp), intent(inout) :: a(NCAT), v(NCAT), s(NCAT), trcrn(:, :)
      logical, intent(in) :: rebin_fb
      real(wp), intent(out) :: aice, aice0
      integer(int32), intent(inout) :: ierr
      real(wp), intent(in), optional :: hi_min_in
      type(ice_params_t) :: p
      real(wp) :: tr(ICE_NTRCR_MAX, NCAT), hm
      integer :: ntr

      call params(kset, nilyr, p)
      hm = p%thermo%hi_min
      if (present(hi_min_in)) hm = hi_min_in
      ntr = ice_ntrcr(nilyr, 1)
      tr = 0.0_wp
      tr(1:ntr, :) = trcrn(1:ntr, :)
      call itd_linear_itd(NCAT, ICE_NTRCR_MAX, nilyr, 1, HIN, ai0, vi0, a, v, s, tr, TF, p%c%rhos, p%c%Lfresh, &
                          hm, PUNY, rebin_fb, aice, aice0, ierr)
      trcrn(1:ntr, :) = tr(1:ntr, :)
   end subroutine call_packed

   subroutine totals(nilyr, a, v, s, trcrn, va, vi, vs, ei, es, si)
      !! Column sums: area, ice and snow volume, ice and snow energy, salt (`ntr` packed tracers).
      integer, intent(in) :: nilyr
      real(wp), intent(in) :: a(NCAT), v(NCAT), s(NCAT), trcrn(:, :)
      real(wp), intent(out) :: va, vi, vs, ei, es, si
      integer :: n, k

      va = 0.0_wp; vi = 0.0_wp; vs = 0.0_wp; ei = 0.0_wp; es = 0.0_wp; si = 0.0_wp
      do n = 1, NCAT
         va = va + a(n)
         vi = vi + v(n)
         vs = vs + s(n)
         do k = 1, nilyr
            ei = ei + trcrn(1 + k, n)*v(n)/real(nilyr, wp)
            si = si + trcrn(2 + nilyr + k, n)*v(n)/real(nilyr, wp)
         end do
         es = es + trcrn(2 + nilyr, n)*s(n)
      end do
   end subroutine totals

   ! ------------------------------------------------------------------
   ! analytical tests
   ! ------------------------------------------------------------------

   subroutine test_known_answer(error)
      !! Category 2 alone, area 0.4, mean thickness 1.30 m, grown uniformly by 0.1 m: the new
      !! boundaries are `hin_max(n) + 0.1`, the mean (1.4) lies in the upper third of
      !! `[hbnew(1), hbnew(2)]`, so `g` lives on `[3*1.4 - 2*hbnew(2), hbnew(2)]`, and the ice above
      !! `hin_max(2)` moves up.  `daice`, `dvice` are the integrals of `g` and `h*g`, written out
      !! here; snow moves with the area, every tracer with the donor's value.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nilyr = 4, ntr = 10
      real(wp) :: ai0(NCAT), vi0(NCAT), a(NCAT), v(NCAT), s(NCAT), trcrn(ntr, NCAT), aice, aice0
      real(wp) :: hb1, hb2, hL, hR, dhr, wk1, wk2, g0, g1, e1, e2, x1, x2, x3, da, dv, ds
      real(wp) :: q0(ntr), tr0(ntr, NCAT), va0, vi0t, vs0, ei0, es0, si0, va1, vi1, vs1, ei1, es1, si1
      type(ice_params_t) :: pk
      real(wp) :: rhosL
      integer(int32) :: ierr
      integer :: it, n

      call params(ICE_CONSTANTS_ICEPACK, nilyr, pk)
      rhosL = pk%c%rhos*pk%c%Lfresh
      ai0 = [0.0_wp, 0.4_wp, 0.0_wp, 0.0_wp, 0.0_wp]
      vi0 = [0.0_wp, 0.52_wp, 0.0_wp, 0.0_wp, 0.0_wp]
      a = ai0
      v = [0.0_wp, 0.56_wp, 0.0_wp, 0.0_wp, 0.0_wp]
      s = [0.0_wp, 0.06_wp, 0.0_wp, 0.0_wp, 0.0_wp]
      q0 = [-7.0_wp, -3.0e8_wp, -2.8e8_wp, -2.6e8_wp, -2.4e8_wp, -1.1e8_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp]
      trcrn = 123.0_wp
      trcrn(:, 2) = q0
      tr0 = trcrn
      ierr = 0_int32
      call totals(nilyr, a, v, s, trcrn, va0, vi0t, vs0, ei0, es0, si0)
      call call_packed(ICE_CONSTANTS_ICEPACK, nilyr, ai0, vi0, a, v, s, trcrn, .false., aice, aice0, ierr)
      call check(error, ierr, 0_int32, "clean")
      if (allocated(error)) return

      ! closed form
      hb1 = HIN(1) + 0.1_wp
      hb2 = HIN(2) + 0.1_wp
      call check(error, 1.4_wp > (hb1 + 2.0_wp*hb2)/3.0_wp, "mean in the upper third of the new category")
      if (allocated(error)) return
      hL = 3.0_wp*1.4_wp - 2.0_wp*hb2
      hR = hb2
      dhr = 1.0_wp/(hR - hL)
      wk1 = 6.0_wp*0.4_wp*dhr
      wk2 = (1.4_wp - hL)*dhr
      g0 = wk1*(2.0_wp/3.0_wp - wk2)
      g1 = 2.0_wp*dhr*wk1*(wk2 - 0.5_wp)
      e1 = max(HIN(2), hL) - hL
      e2 = min(hb2, hR) - hL
      x1 = e2 - e1
      x2 = 0.5_wp*(e2*e2 - e1*e1)
      x3 = (e2*e2*e2 - e1*e1*e1)/3.0_wp
      da = g1*x2 + g0*x1
      dv = g1*x3 + g0*x2 + da*hL
      ds = 0.06_wp*da/0.4_wp
      call check(error, da > 0.0_wp .and. da < 0.4_wp, "transferred area inside (0, aicen)")
      if (allocated(error)) return
      call check(error, dv/da > HIN(2) .and. dv/da < hb2, "transferred mean thickness between the two bounds")
      if (allocated(error)) return

      call check(error, abs(a(3) - da) <= 1.0e-14_wp, "area received")
      if (allocated(error)) return
      call check(error, abs(a(2) - (0.4_wp - da)) <= 1.0e-14_wp, "area left")
      if (allocated(error)) return
      call check(error, abs(v(3) - dv) <= 1.0e-14_wp .and. abs(v(2) - (0.56_wp - dv)) <= 1.0e-14_wp, "ice volume moved")
      if (allocated(error)) return
      call check(error, abs(s(3) - ds) <= 1.0e-14_wp .and. abs(s(2) - (0.06_wp - ds)) <= 1.0e-14_wp, &
                 "snow moves in the area fraction da/aicen_init")
      if (allocated(error)) return
      call check(error, all(a([1, 4, 5]) == 0.0_wp) .and. all(v([1, 4, 5]) == 0.0_wp), "other categories stay empty")
      if (allocated(error)) return
      do it = 2, ntr   ! the receiver carries the DONOR's intensive tracers (it was empty)
         call check(error, abs(trcrn(it, 3) - q0(it)) <= 1.0e-13_wp*abs(q0(it)), "receiver tracer = donor tracer")
         if (allocated(error)) return
         call check(error, abs(trcrn(it, 2) - q0(it)) <= 1.0e-13_wp*abs(q0(it)), "donor tracer unchanged")
         if (allocated(error)) return
      end do
      call check(error, abs(trcrn(1, 3) - q0(1)) <= 1.0e-13_wp*abs(q0(1)), "Tsfc moves with the area")
      if (allocated(error)) return
      do n = 1, NCAT   ! empty categories: Tsfc = Tf, ice tracers 0, snow enthalpy = -rhos*Lfresh (upstream's offset)
         if (a(n) <= PUNY) then
            call check(error, trcrn(1, n) == TF .and. all(trcrn(2:5, n) == 0.0_wp) .and. all(trcrn(7:10, n) == 0.0_wp) &
                       .and. abs(trcrn(6, n) + rhosL) <= 1.0e-6_wp, "empty category tracers")
            if (allocated(error)) return
         end if
      end do
      call totals(nilyr, a, v, s, trcrn, va1, vi1, vs1, ei1, es1, si1)
      call check(error, abs(va1 - va0) <= 1.0e-15_wp .and. abs(vi1 - vi0t) <= 1.0e-15_wp .and. &
                 abs(vs1 - vs0) <= 1.0e-15_wp, "area and volumes conserved")
      if (allocated(error)) return
      call check(error, abs(ei1 - ei0) <= 1.0e-13_wp*abs(ei0) .and. abs(es1 - es0) <= 1.0e-13_wp*abs(es0) &
                 .and. abs(si1 - si0) <= 1.0e-13_wp*abs(si0), "energy and salt conserved")
      if (allocated(error)) return
      call check(error, abs(aice - 0.4_wp) <= 1.0e-15_wp .and. abs(aice0 - 0.6_wp) <= 1.0e-15_wp, "aggregates")
   end subroutine test_known_answer

   subroutine test_identity(error)
      !! Nothing grew or melted: area, volume and snow are bitwise unchanged; the tracers
      !! (recomputed by the shift's division, as upstream) agree to round-off.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nilyr = 7, ntr = 16
      real(wp) :: ai0(NCAT), vi0(NCAT), a(NCAT), v(NCAT), s(NCAT), a1(NCAT), v1(NCAT), s1(NCAT)
      real(wp) :: trcrn(ntr, NCAT), tr1(ntr, NCAT), aice, aice0
      type(ice_params_t) :: pk
      integer(int32) :: ierr
      integer :: n, it

      call params(ICE_CONSTANTS_ICEPACK, nilyr, pk)
      ai0 = [0.1_wp, 0.15_wp, 0.2_wp, 0.12_wp, 0.08_wp]
      vi0 = ai0*[0.3_wp, 1.0_wp, 1.9_wp, 3.5_wp, 6.0_wp]
      a = ai0; v = vi0
      s = ai0*[0.0_wp, 0.1_wp, 0.2_wp, 0.1_wp, 0.3_wp]
      call lcg_tracers(1, nilyr, ntr, trcrn)
      a1 = a; v1 = v; s1 = s; tr1 = trcrn
      ierr = 0_int32
      call call_packed(ICE_CONSTANTS_ICEPACK, nilyr, ai0, vi0, a, v, s, trcrn, .false., aice, aice0, ierr)
      call check(error, ierr, 0_int32, "clean")
      if (allocated(error)) return
      call check(error, all(a == a1) .and. all(v == v1) .and. all(s == s1), "area and volumes bitwise unchanged")
      if (allocated(error)) return
      do n = 1, NCAT
         if (s(n) > 0.0_wp) then
            call check(error, all(abs(trcrn(:, n) - tr1(:, n)) <= 1.0e-14_wp*abs(tr1(:, n))), "tracers unchanged to round-off")
         else   ! snow-free category: the snow enthalpy is the melting-point value -rhos*Lfresh (upstream)
            call check(error, all(abs(trcrn([(it, it=1, 8), (it, it=10, 16)], n) - tr1([(it, it=1, 8), (it, it=10, 16)], n)) &
                                  <= 1.0e-14_wp*abs(tr1([(it, it=1, 8), (it, it=10, 16)], n))) .and. &
                       abs(trcrn(9, n) + pk%c%rhos*pk%c%Lfresh) <= 1.0e-6_wp, "snow-free category: qsno at the melting point")
         end if
         if (allocated(error)) return
      end do
      call check(error, abs(aice - sum(a1)) <= 1.0e-15_wp, "aice")
   end subroutine test_identity

   subroutine test_cat1_melt(error)
      !! Category 1 alone (area 0.3, mean 0.3 m) melts by 0.1 m: the area under `g(h)` below
      !! `|dh0|` (a sliver of the uniform-middle-third fit on `[0, hin_max(1)]`) leaves for open
      !! water, the volume stays, so the thickness rises by `aicen/(aicen - da0)`; no transfer.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nilyr = 4, ntr = 10
      real(wp) :: ai0(NCAT), vi0(NCAT), a(NCAT), v(NCAT), s(NCAT), trcrn(ntr, NCAT), aice, aice0
      real(wp) :: dhr, wk1, wk2, g0, g1, da0, damax, hnew
      integer(int32) :: ierr

      ai0 = [0.3_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
      vi0 = [0.09_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
      a = ai0
      v = [0.06_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
      s = 0.0_wp
      call lcg_tracers(2, nilyr, ntr, trcrn)
      ierr = 0_int32
      call call_packed(ICE_CONSTANTS_ICEPACK, nilyr, ai0, vi0, a, v, s, trcrn, .false., aice, aice0, ierr)
      call check(error, ierr, 0_int32, "clean")
      if (allocated(error)) return

      call check(error, 0.3_wp > (2.0_wp*0.0_wp + HIN(1))/3.0_wp .and. 0.3_wp < (0.0_wp + 2.0_wp*HIN(1))/3.0_wp, &
                 "mean in the central third: g spans the whole category")
      if (allocated(error)) return
      dhr = 1.0_wp/HIN(1)
      wk1 = 6.0_wp*0.3_wp*dhr
      wk2 = 0.3_wp*dhr
      g0 = wk1*(2.0_wp/3.0_wp - wk2)
      g1 = 2.0_wp*dhr*wk1*(wk2 - 0.5_wp)
      da0 = g1*0.5_wp*0.1_wp*0.1_wp + g0*0.1_wp      ! area under g between 0 and |dh0| = 0.1
      damax = 0.3_wp*(1.0_wp - 0.2_wp/0.3_wp)
      call check(error, da0 > 0.0_wp .and. da0 < damax, "sliver smaller than the cap")
      if (allocated(error)) return
      hnew = 0.2_wp*0.3_wp/(0.3_wp - da0)
      call check(error, abs(a(1) - (0.3_wp - da0)) <= 1.0e-15_wp, "area removed by the sliver")
      if (allocated(error)) return
      call check(error, v(1) == 0.06_wp, "volume conserved")
      if (allocated(error)) return
      call check(error, abs(v(1)/a(1) - hnew) <= 1.0e-14_wp, "thickness rises by aicen/(aicen - da0)")
      if (allocated(error)) return
      call check(error, abs(aice0 - (1.0_wp - (0.3_wp - da0))) <= 1.0e-15_wp, "the area went to open water")
   end subroutine test_cat1_melt

   subroutine test_hi_min(error)
      !! Category 1 thinner than `hi_min` gives up area until its thickness is `hi_min`, volume kept;
      !! `hi_min = 0` switches the floor off.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nilyr = 4, ntr = 10
      real(wp) :: ai0(NCAT), vi0(NCAT), a(NCAT), v(NCAT), s(NCAT), trcrn(ntr, NCAT), aice, aice0
      integer(int32) :: ierr

      ai0 = [0.5_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
      vi0 = [0.002_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
      s = 0.0_wp
      a = ai0; v = vi0
      call lcg_tracers(3, nilyr, ntr, trcrn)
      ierr = 0_int32
      call call_packed(ICE_CONSTANTS_ICEPACK, nilyr, ai0, vi0, a, v, s, trcrn, .false., aice, aice0, ierr)
      call check(error, ierr, 0_int32, "clean")
      if (allocated(error)) return
      call check(error, abs(a(1) - 0.5_wp*(0.004_wp/0.01_wp)) <= 1.0e-16_wp, "area aicen*h/hi_min")
      if (allocated(error)) return
      call check(error, v(1) == 0.002_wp, "volume unchanged")
      if (allocated(error)) return
      call check(error, abs(v(1)/a(1) - 0.01_wp) <= 1.0e-15_wp, "thickness = hi_min")
      if (allocated(error)) return
      call check(error, abs(aice - a(1)) <= 1.0e-16_wp .and. abs(aice0 - (1.0_wp - a(1))) <= 1.0e-16_wp, "aggregates")
      if (allocated(error)) return

      a = ai0; v = vi0
      call lcg_tracers(3, nilyr, ntr, trcrn)
      call call_packed(ICE_CONSTANTS_ICEPACK, nilyr, ai0, vi0, a, v, s, trcrn, .false., aice, aice0, ierr, 0.0_wp)
      call check(error, a(1) == 0.5_wp, "hi_min = 0: no floor")
   end subroutine test_hi_min

   subroutine test_cat1_zero_volume(error)
      !! Regression: category 1 keeps area above `puny` but has no volume left (everything melted:
      !! `hicen = 0`, `damax = aicen`).  The melt-area loss then divides `0*aicen/(aicen - da0)`,
      !! which is 0/0 if `da0` reaches `aicen`.  Scan initial mean thicknesses over the whole bin
      !! (including one just below its top, where the fitted support ends at the bound) and two
      !! areas: no output may be non-finite, nothing may leave the cell but the area lost to open water.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nilyr = 4, ntr = 10
      real(wp) :: ai0(NCAT), vi0(NCAT), a(NCAT), v(NCAT), s(NCAT), trcrn(ntr, NCAT), aice, aice0, h0, a0
      integer(int32) :: ierr
      integer :: k, m

      do m = 1, 2
         a0 = merge(0.3_wp, 1.0e-9_wp, m == 1)
         do k = 1, 40
            h0 = HIN(1)*real(k, wp)/40.0_wp - 1.0e-7_wp
            ai0 = [a0, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
            vi0 = [a0*h0, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
            a = ai0
            v = 0.0_wp
            s = 0.0_wp
            call lcg_tracers(10 + k, nilyr, ntr, trcrn)
            ierr = 0_int32
            call call_packed(ICE_CONSTANTS_ICEPACK, nilyr, ai0, vi0, a, v, s, trcrn, .false., aice, aice0, ierr)
            call check(error, iand(ierr, not(ICE_ERR_REMAP_FALLBACK)) == 0_int32, "no error bit")
            if (allocated(error)) return
            call check(error, all(ieee_is_finite(a)) .and. all(ieee_is_finite(v)) .and. all(ieee_is_finite(s)) .and. &
                       all(ieee_is_finite(trcrn)) .and. ieee_is_finite(aice) .and. ieee_is_finite(aice0), &
                       "outputs finite when category 1 has no volume")
            if (allocated(error)) return
            call check(error, all(v == 0.0_wp) .and. all(a >= 0.0_wp) .and. a(1) <= a0, "no volume appears, area only shrinks")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_cat1_zero_volume

   subroutine test_fallback(error)
      !! A boundary that jumps past the mean thickness of its category abandons the remap: the
      !! state is untouched, `REMAP_FALLBACK` set.  The rebin then moves the whole category and
      !! merges it with the receiver by closed forms (area-, volume- and snow-volume-weighted).
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nilyr = 4, ntr = 10
      real(wp) :: ai0(NCAT), vi0(NCAT), a(NCAT), v(NCAT), s(NCAT), trcrn(ntr, NCAT), tr0(ntr, NCAT), aice, aice0
      real(wp) :: a1(NCAT), v1(NCAT), s1(NCAT), tr1(ntr, NCAT), w2, w3
      integer(int32) :: ierr
      integer :: it

      ai0 = [0.0_wp, 0.3_wp, 0.3_wp, 0.0_wp, 0.0_wp]
      vi0 = ai0*[0.0_wp, 1.0_wp, 2.0_wp, 0.0_wp, 0.0_wp]
      a = ai0
      v = ai0*[0.0_wp, 2.5_wp, 2.1_wp, 0.0_wp, 0.0_wp]
      s = [0.0_wp, 0.06_wp, 0.03_wp, 0.0_wp, 0.0_wp]
      call lcg_tracers(4, nilyr, ntr, trcrn)
      tr0 = trcrn
      a1 = a; v1 = v; s1 = s
      ierr = 0_int32
      call call_packed(ICE_CONSTANTS_ICEPACK, nilyr, ai0, vi0, a, v, s, trcrn, .false., aice, aice0, ierr)
      call check(error, iand(ierr, ICE_ERR_REMAP_FALLBACK) /= 0_int32, "fallback flagged")
      if (allocated(error)) return
      call check(error, iand(ierr, not(ICE_ERR_REMAP_FALLBACK)) == 0_int32, "nothing else flagged")
      if (allocated(error)) return
      call check(error, all(a == a1) .and. all(v == v1) .and. all(s == s1) .and. all(trcrn == tr0), "state untouched")
      if (allocated(error)) return
      call check(error, abs(aice - 0.6_wp) <= 1.0e-15_wp, "aggregate of the untouched state")
      if (allocated(error)) return

      ! with the rebin: category 2 (2.5 m > hin_max(2)) moves wholly into category 3
      ierr = 0_int32
      call call_packed(ICE_CONSTANTS_ICEPACK, nilyr, ai0, vi0, a, v, s, trcrn, .true., aice, aice0, ierr)
      call check(error, a(2) == 0.0_wp .and. v(2) == 0.0_wp .and. s(2) == 0.0_wp, "category 2 emptied")
      if (allocated(error)) return
      call check(error, abs(a(3) - 0.6_wp) <= 1.0e-15_wp .and. abs(v(3) - (0.75_wp + 0.63_wp)) <= 1.0e-15_wp &
                 .and. abs(s(3) - 0.09_wp) <= 1.0e-15_wp, "area, volume, snow merged in category 3")
      if (allocated(error)) return
      call check(error, v(3)/a(3) > HIN(2) .and. v(3)/a(3) <= HIN(3), "category 3 inside its bin")
      if (allocated(error)) return
      w2 = 0.3_wp/0.6_wp; w3 = 0.3_wp/0.6_wp
      tr1 = tr0
      call check(error, abs(trcrn(1, 3) - (tr1(1, 2)*w2 + tr1(1, 3)*w3)) <= 1.0e-13_wp*abs(trcrn(1, 3)), "Tsfc area-weighted")
      if (allocated(error)) return
      do it = 2, 5   ! qice layers: ice-volume-weighted
         call check(error, abs(trcrn(it, 3) - (tr1(it, 2)*0.75_wp + tr1(it, 3)*0.63_wp)/1.38_wp) <= &
                    1.0e-13_wp*abs(trcrn(it, 3)), "qice volume-weighted")
         if (allocated(error)) return
      end do
      call check(error, abs(trcrn(6, 3) - (tr1(6, 2)*0.06_wp + tr1(6, 3)*0.03_wp)/0.09_wp) <= 1.0e-13_wp*abs(trcrn(6, 3)), &
                 "qsno snow-volume-weighted")
      if (allocated(error)) return
      do it = 7, 10
         call check(error, abs(trcrn(it, 3) - (tr1(it, 2)*0.75_wp + tr1(it, 3)*0.63_wp)/1.38_wp) <= &
                    1.0e-13_wp*abs(trcrn(it, 3)), "sice volume-weighted")
         if (allocated(error)) return
      end do
   end subroutine test_fallback

   subroutine test_shift_range(error)
      !! Increments outside the tolerance band refuse the shift (`ITD_SHIFT_RANGE`, nothing
      !! modified); inside it they are trimmed (negative round-off to zero, `1 -+ puny` of the
      !! donor to the whole category).
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nilyr = 4, ntr = 10
      real(wp) :: a(NCAT), v(NCAT), s(NCAT), trcrn(ICE_NTRCR_MAX, NCAT), a0(NCAT), v0(NCAT), tr0(ICE_NTRCR_MAX, NCAT)
      real(wp) :: daice(NCAT), dvice(NCAT), hicen(NCAT)
      integer :: donor(NCAT), k
      integer(int32) :: ierr

      donor = [1, 0, 0, 0, 0]
      do k = 1, 4
         a = [0.3_wp, 0.2_wp, 0.0_wp, 0.0_wp, 0.0_wp]
         v = [0.15_wp, 0.3_wp, 0.0_wp, 0.0_wp, 0.0_wp]
         s = [0.03_wp, 0.01_wp, 0.0_wp, 0.0_wp, 0.0_wp]
         call lcg_tracers(5, nilyr, ntr, trcrn(1:ntr, :))
         trcrn(ntr + 1:, :) = 0.0_wp
         a0 = a; v0 = v; tr0 = trcrn
         select case (k)
         case (1)
            daice = [0.3_wp*(1.0_wp + 1.0e-9_wp), 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]   ! daice > aicen(1 + puny)
            dvice = [0.1_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
         case (2)
            daice = [0.1_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
            dvice = [0.15_wp*(1.0_wp + 1.0e-9_wp), 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]   ! dvice > vicen(1 + puny)
         case (3)
            daice = [-0.01_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]                        ! daice < -puny*aicen
            dvice = [0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
         case default
            daice = [0.1_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
            dvice = [-0.01_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]                        ! dvice < -puny*vicen
         end select
         ierr = 0_int32
         call itd_shift_ice(NCAT, ICE_NTRCR_MAX, nilyr, 1, PUNY, TF, 0.0_wp, donor, daice, dvice, a, v, s, hicen, trcrn, ierr)
         call check(error, iand(ierr, ICE_ERR_ITD_SHIFT_RANGE) /= 0_int32, "range violation flagged")
         if (allocated(error)) return
         call check(error, all(a == a0) .and. all(v == v0) .and. all(trcrn == tr0), "refused shift leaves the state")
         if (allocated(error)) return
      end do

      ! trims: round-off negative -> no shift; 1 - 0.5 puny -> whole category
      a = [0.3_wp, 0.2_wp, 0.0_wp, 0.0_wp, 0.0_wp]
      v = [0.15_wp, 0.3_wp, 0.0_wp, 0.0_wp, 0.0_wp]
      s = [0.03_wp, 0.01_wp, 0.0_wp, 0.0_wp, 0.0_wp]
      call lcg_tracers(5, nilyr, ntr, trcrn(1:ntr, :))
      a0 = a; v0 = v
      daice = [-0.5e-11_wp*0.3_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
      dvice = 0.0_wp
      ierr = 0_int32
      call itd_shift_ice(NCAT, ICE_NTRCR_MAX, nilyr, 1, PUNY, TF, 0.0_wp, donor, daice, dvice, a, v, s, hicen, trcrn, ierr)
      call check(error, ierr == 0_int32 .and. all(a == a0) .and. all(v == v0), "negative round-off: shift nothing")
      if (allocated(error)) return
      daice = [0.3_wp*(1.0_wp - 0.5e-11_wp), 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
      dvice = [0.15_wp*(1.0_wp - 0.5e-11_wp), 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
      call itd_shift_ice(NCAT, ICE_NTRCR_MAX, nilyr, 1, PUNY, TF, 0.0_wp, donor, daice, dvice, a, v, s, hicen, trcrn, ierr)
      call check(error, ierr == 0_int32 .and. a(1) == 0.0_wp .and. v(1) == 0.0_wp .and. abs(a(2) - 0.5_wp) <= 1.0e-16_wp, &
                 "within puny of the donor: the whole category moves")
   end subroutine test_shift_range

   subroutine test_error_bits(error)
      !! `HBNEW_ORDER` (initial thicknesses of two populated neighbours not increasing; upstream
      !! aborts), non-finite inputs (area, volume, tracer, `Tf`): the cell is untouched and flagged.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nilyr = 4, ntr = 10
      real(wp) :: ai0(NCAT), vi0(NCAT), a(NCAT), v(NCAT), s(NCAT), trcrn(ntr, NCAT), aice, aice0
      real(wp) :: a1(NCAT), v1(NCAT), s1(NCAT), tr1(ntr, NCAT), tfx
      type(ice_params_t) :: p
      integer(int32) :: ierr
      integer :: k
      real(wp), parameter :: QNAN = transfer(int(z'7FF8000000000000', int64), 1.0_wp)

      ! thicker initial mean in category 2 than in category 3
      ai0 = [0.0_wp, 0.3_wp, 0.3_wp, 0.0_wp, 0.0_wp]
      vi0 = ai0*[0.0_wp, 1.3_wp, 1.2_wp, 0.0_wp, 0.0_wp]
      a = ai0; v = vi0*1.01_wp; s = 0.0_wp
      call lcg_tracers(6, nilyr, ntr, trcrn)
      a1 = a; v1 = v; s1 = s; tr1 = trcrn
      ierr = 0_int32
      call call_packed(ICE_CONSTANTS_ICEPACK, nilyr, ai0, vi0, a, v, s, trcrn, .false., aice, aice0, ierr)
      call check(error, iand(ierr, ICE_ERR_ITD_HBNEW_ORDER) /= 0_int32 .and. iand(ierr, ICE_ERR_REMAP_FALLBACK) /= 0_int32, &
                 "order violation: HBNEW_ORDER and the fallback")
      if (allocated(error)) return
      call check(error, all(a == a1) .and. all(v == v1) .and. all(trcrn == tr1), "state untouched")
      if (allocated(error)) return

      ! non-finite inputs
      call params(ICE_CONSTANTS_ICEPACK, nilyr, p)
      do k = 1, 5
         ai0 = [0.1_wp, 0.15_wp, 0.2_wp, 0.12_wp, 0.08_wp]
         vi0 = ai0*[0.3_wp, 1.0_wp, 1.9_wp, 3.5_wp, 6.0_wp]
         a = ai0; v = vi0*1.02_wp; s = 0.0_wp
         call lcg_tracers(7, nilyr, ntr, trcrn)
         tfx = TF
         select case (k)
         case (1); a(3) = QNAN
         case (2); v(2) = QNAN
         case (3); trcrn(4, 5) = QNAN
         case (4); tfx = QNAN
         case default; ai0(1) = QNAN
         end select
         a1 = a; v1 = v; s1 = s; tr1 = trcrn
         ierr = 0_int32
         call run_one_packed(nilyr, p, tfx, ai0, vi0, a, v, s, trcrn, aice, aice0, ierr)
         call check(error, iand(ierr, ICE_ERR_NONFINITE) /= 0_int32, "non-finite input flagged")
         if (allocated(error)) return
         call check(error,.not. ieee_is_finite(aice) .and. .not. ieee_is_finite(aice0), "aggregates are NaN, not laundered")
         if (allocated(error)) return
         call check(error, all(same_arr(a, a1)) .and. all(same_arr(v, v1)) .and. all(same_arr(reshape(trcrn, [ntr*NCAT]), &
                                                                                              reshape(tr1, [ntr*NCAT]))), &
                    "input untouched")
         if (allocated(error)) return
      end do
   end subroutine test_error_bits

   elemental logical function same_arr(a, b)
      real(wp), intent(in) :: a, b

      same_arr = same(a, b, 0.0_wp)
   end function same_arr

   subroutine run_one_packed(nilyr, p, tfx, ai0, vi0, a, v, s, trcrn, aice, aice0, ierr)
      integer, intent(in) :: nilyr
      type(ice_params_t), intent(in) :: p
      real(wp), intent(in) :: tfx, ai0(NCAT), vi0(NCAT)
      real(wp), intent(inout) :: a(NCAT), v(NCAT), s(NCAT), trcrn(:, :)
      real(wp), intent(out) :: aice, aice0
      integer(int32), intent(inout) :: ierr
      real(wp) :: tr(ICE_NTRCR_MAX, NCAT)
      integer :: ntr

      ntr = ice_ntrcr(nilyr, 1)
      tr = 0.0_wp
      tr(1:ntr, :) = trcrn(1:ntr, :)
      call itd_linear_itd(NCAT, ICE_NTRCR_MAX, nilyr, 1, HIN, ai0, vi0, a, v, s, tr, tfx, p%c%rhos, p%c%Lfresh, &
                          p%thermo%hi_min, PUNY, .false., aice, aice0, ierr)
      trcrn(1:ntr, :) = tr(1:ntr, :)
   end subroutine run_one_packed

   ! ------------------------------------------------------------------
   ! random state
   ! ------------------------------------------------------------------

   subroutine fill_random(nx, ny, nilyr, seed, ai0, vi0, a, v, s, Tsfcn, qicen, sicen, qsnon, mode)
      !! Random consistent cells: each category populated with probability 0.8 (area 0.02 - 0.27),
      !! mean thickness inside its bin (up to 8 m), growth `dh = scale*0.2*(u - 0.4)` with the
      !! scale cycling slow ... wild per cell; a quarter of the cells also perturb the current
      !! area against the initial one.  `mode = 1` (default): as above; `mode = 2`: mean
      !! thicknesses drawn anywhere 0.05 - 8 m (outside their bins: rebin food), no growth.
      integer, intent(in) :: nx, ny, nilyr
      integer(int64), intent(in) :: seed
      real(wp), intent(out) :: ai0(nx, ny, NCAT), vi0(nx, ny, NCAT), a(nx, ny, NCAT), v(nx, ny, NCAT), s(nx, ny, NCAT)
      real(wp), intent(out) :: Tsfcn(nx, ny, NCAT), qicen(nx, ny, NCAT, nilyr), sicen(nx, ny, NCAT, nilyr), qsnon(nx, ny, NCAT, 1)
      integer, intent(in), optional :: mode
      integer(int64) :: st
      integer :: i, j, n, k, md
      real(wp) :: u1, u2, u3, u4, lo, hi, h0, dh, sc, pert

      md = 1
      if (present(mode)) md = mode
      do j = 1, ny
         do i = 1, nx
            st = lcg_next(seed + int(i, int64)*7919_int64 + int(j, int64)*104729_int64*65537_int64)
            st = lcg_next(st); u1 = lcg_unit(st)
            sc = SCL(1 + int(u1*16.0_wp))
            st = lcg_next(st); pert = 0.0_wp
            if (lcg_unit(st) < 0.25_wp) pert = 0.1_wp
            do n = 1, NCAT
               st = lcg_next(st); u1 = lcg_unit(st)
               st = lcg_next(st); u2 = lcg_unit(st)
               st = lcg_next(st); u3 = lcg_unit(st)
               st = lcg_next(st); u4 = lcg_unit(st)
               ai0(i, j, n) = 0.0_wp; vi0(i, j, n) = 0.0_wp; a(i, j, n) = 0.0_wp; v(i, j, n) = 0.0_wp; s(i, j, n) = 0.0_wp
               if (u1 < 0.8_wp) then
                  ai0(i, j, n) = 0.02_wp + 0.25_wp*u2
                  if (md == 1) then
                     lo = HIN(n - 1); hi = min(HIN(n), 8.0_wp)
                     h0 = lo + (hi - lo)*(0.05_wp + 0.9_wp*u3)
                     dh = sc*0.2_wp*(u4 - 0.4_wp)
                     dh = max(dh, -0.9_wp*h0)
                     if (h0 + dh < 0.002_wp) dh = 0.002_wp - h0
                  else
                     h0 = 0.05_wp + 7.95_wp*u3
                     dh = 0.0_wp
                  end if
                  vi0(i, j, n) = ai0(i, j, n)*h0
                  st = lcg_next(st)
                  a(i, j, n) = ai0(i, j, n)*(1.0_wp + pert*(lcg_unit(st) - 0.5_wp))
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
               st = lcg_next(st); qsnon(i, j, n, 1) = -1.2e8_wp + 0.3e8_wp*lcg_unit(st)
            end do
         end do
      end do
   end subroutine fill_random

   subroutine gather(i, j, nx, ny, nilyr, ai0, vi0, a, v, s, Tsfcn, qicen, sicen, qsnon, ai0c, vi0c, ac, vc, sc, tr)
      !! Packed column of cell `(i, j)` (`ICE_NTRCR_MAX`-leading tracer block).
      integer, intent(in) :: i, j, nx, ny, nilyr
      real(wp), intent(in) :: ai0(nx, ny, NCAT), vi0(nx, ny, NCAT), a(nx, ny, NCAT), v(nx, ny, NCAT), s(nx, ny, NCAT)
      real(wp), intent(in) :: Tsfcn(nx, ny, NCAT), qicen(nx, ny, NCAT, nilyr), sicen(nx, ny, NCAT, nilyr), qsnon(nx, ny, NCAT, 1)
      real(wp), intent(out) :: ai0c(NCAT), vi0c(NCAT), ac(NCAT), vc(NCAT), sc(NCAT), tr(ICE_NTRCR_MAX, NCAT)
      integer :: n, it

      tr = 0.0_wp
      do n = 1, NCAT
         ai0c(n) = ai0(i, j, n); vi0c(n) = vi0(i, j, n)
         ac(n) = a(i, j, n); vc(n) = v(i, j, n); sc(n) = s(i, j, n)
         do it = 1, ice_ntrcr(nilyr, 1)
            tr(it, n) = ice_trcr_get(it, i, j, n, nx, ny, NCAT, nilyr, 1, Tsfcn, qicen, sicen, qsnon)
         end do
      end do
   end subroutine gather

   subroutine conserve_case(error, nilyr, seed)
      !! 100 000 random cells (500 x 200) through the host shim: every closure per cell.
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: nilyr
      integer(int64), intent(in) :: seed
      integer, parameter :: nx = 500, ny = 200
      real(wp), allocatable :: ai0(:, :, :), vi0(:, :, :), a(:, :, :), v(:, :, :), s(:, :, :), Tsfcn(:, :, :)
      real(wp), allocatable :: qicen(:, :, :, :), sicen(:, :, :, :), qsnon(:, :, :, :)
      real(wp), allocatable :: Tf2(:, :), aice(:, :), aice0(:, :)
real(wp), allocatable :: a_b(:, :, :), v_b(:, :, :), s_b(:, :, :), T_b(:, :, :), q_b(:, :, :, :), si_b(:, :, :, :), qs_b(:, :, :, :)
      integer(int32), allocatable :: ierr(:, :)
      type(ice_params_t) :: p
      real(wp) :: hin_loc(0:NCAT), wa, wv, ws, we, wq, wsl, worst_a, da, cell_a0, h
      real(wp) :: va0, vi0t, vs0, ei0, es0, si0, va1, vi1, vs1, ei1, es1, si1, sc_e, sc_s
      real(wp) :: tr0(ICE_NTRCR_MAX, NCAT), tr1(ICE_NTRCR_MAX, NCAT), Ta0, Ta1
      integer :: i, j, n, k, nfb, nlossy, nbin, nact, nrb
      logical :: noloss_cell

      call params(ICE_CONSTANTS_ICEPACK, nilyr, p)
      allocate (ai0(nx, ny, NCAT), vi0(nx, ny, NCAT), a(nx, ny, NCAT), v(nx, ny, NCAT), s(nx, ny, NCAT), Tsfcn(nx, ny, NCAT), &
                qicen(nx, ny, NCAT, nilyr), sicen(nx, ny, NCAT, nilyr), qsnon(nx, ny, NCAT, 1), Tf2(nx, ny), aice(nx, ny), &
                aice0(nx, ny), ierr(nx, ny))
      allocate (a_b(nx, ny, NCAT), v_b(nx, ny, NCAT), s_b(nx, ny, NCAT), T_b(nx, ny, NCAT), q_b(nx, ny, NCAT, nilyr), &
                si_b(nx, ny, NCAT, nilyr), qs_b(nx, ny, NCAT, 1))
      call fill_random(nx, ny, nilyr, seed, ai0, vi0, a, v, s, Tsfcn, qicen, sicen, qsnon)
      Tf2 = TF
      ierr = 0_int32
      aice = 0.0_wp
      aice0 = 0.0_wp
      a_b = a; v_b = v; s_b = s; T_b = Tsfcn; q_b = qicen; si_b = sicen; qs_b = qsnon
      hin_loc = HIN

      !$acc enter data copyin(hin_loc, Tf2, ai0, vi0, a, v, s, Tsfcn, qicen, sicen, qsnon, ierr, aice, aice0)
      call ice_linear_itd_step(nx, ny, NCAT, nilyr, 1, 1, nx, 1, ny, p, hin_loc, Tf2, ai0, vi0, a, v, s, Tsfcn, qicen, sicen, &
                               qsnon, .false., aice, aice0, ierr)
      !$acc update self(a, v, s, Tsfcn, qicen, sicen, qsnon, ierr, aice, aice0)
      !$acc exit data delete(hin_loc, Tf2, ai0, vi0, a, v, s, Tsfcn, qicen, sicen, qsnon, ierr, aice, aice0)

      worst_a = 0.0_wp
      wv = 0.0_wp; ws = 0.0_wp; we = 0.0_wp; wq = 0.0_wp; wsl = 0.0_wp
      nfb = 0; nlossy = 0; nbin = 0; nact = 0; nrb = 0
      do j = 1, ny
         do i = 1, nx
            if (.not. all(ieee_is_finite([aice(i, j), aice0(i, j)]))) then
               call test_failed(error, "non-finite aggregate")
               return
            end if
            if (iand(ierr(i, j), not(ICE_ERR_REMAP_FALLBACK)) /= 0_int32) then
               call test_failed(error, "unexpected error bit")
               return
            end if
            if (iand(ierr(i, j), ICE_ERR_REMAP_FALLBACK) /= 0_int32) then
               nfb = nfb + 1
               cycle
            end if
            nact = nact + 1
            do n = 1, NCAT
               do k = 0, 0
                  if (.not. (ieee_is_finite(a(i, j, n)) .and. ieee_is_finite(v(i, j, n)) .and. ieee_is_finite(s(i, j, n)) &
                             .and. ieee_is_finite(Tsfcn(i, j, n)) .and. all(ieee_is_finite(qicen(i, j, n, :))) &
                             .and. all(ieee_is_finite(sicen(i, j, n, :))) .and. ieee_is_finite(qsnon(i, j, n, 1)))) then
                     call test_failed(error, "non-finite state")
                     return
                  end if
               end do
            end do
            ! before/after sums (before: the saved copies; the packed tracer layout of `totals`)
            call pack_tr(i, j, nx, ny, nilyr, T_b, q_b, si_b, qs_b, tr0)
            call pack_tr(i, j, nx, ny, nilyr, Tsfcn, qicen, sicen, qsnon, tr1)
            call totals(nilyr, a_b(i, j, :), v_b(i, j, :), s_b(i, j, :), tr0(1:ntr_of(nilyr), :), va0, vi0t, vs0, ei0, es0, si0)
            call totals(nilyr, a(i, j, :), v(i, j, :), s(i, j, :), tr1(1:ntr_of(nilyr), :), va1, vi1, vs1, ei1, es1, si1)
            if (.not. ieee_is_finite(va0 + vi0t + vs0 + ei0 + es0 + si0 + va1 + vi1 + vs1 + ei1 + es1 + si1)) then
               call test_failed(error, "non-finite column sum")
               return
            end if
            sc_e = 0.0_wp; sc_s = 0.0_wp
            do n = 1, NCAT
               sc_e = sc_e + sum(abs(q_b(i, j, n, :)))*v_b(i, j, n)/real(nilyr, wp)
               sc_s = sc_s + sum(abs(si_b(i, j, n, :)))*v_b(i, j, n)/real(nilyr, wp)
            end do
            wv = max(wv, abs(vi1 - vi0t)/max(vi0t, 1.0e-30_wp))
            ws = max(ws, abs(vs1 - vs0)/max(vs0, 1.0e-30_wp))
            we = max(we, abs(ei1 - ei0)/max(sc_e, 1.0e-30_wp))
            wq = max(wq, abs(es1 - es0)/max(sum(abs(qs_b(i, j, :, 1))*s_b(i, j, :)), 1.0e-30_wp))
            wsl = max(wsl, abs(si1 - si0)/max(sc_s, 1.0e-30_wp))
            ! area: only the category-1 melt loss and the hi_min floor remove any
            da = va0 - va1
            cell_a0 = a_b(i, j, 1)
            h = 0.0_wp
            if (ai0(i, j, 1) > PUNY) h = vi0(i, j, 1)/ai0(i, j, 1)
            noloss_cell = cell_a0 <= PUNY .or. (v_b(i, j, 1)/max(a_b(i, j, 1), 1.0e-30_wp) >= vi0(i, j, 1)/max(ai0(i, j, 1), &
                                                                                                               1.0e-30_wp) &
                                                .and. v_b(i, j, 1)/max(a_b(i, j, 1), 1.0e-30_wp) >= 0.0101_wp*1.0_wp)
            if (noloss_cell) then
               worst_a = max(worst_a, abs(da))
               ! area-weighted surface temperature is conserved when no area leaves
               Ta0 = sum(a_b(i, j, :)*tr0(1, :))
               Ta1 = sum(a(i, j, :)*tr1(1, :))
               if (abs(Ta1 - Ta0) > 1.0e-12_wp*sum(a_b(i, j, :)*abs(tr0(1, :))) + 1.0e-15_wp) then
                  call test_failed(error, "area-weighted Tsfc not conserved in a cell without area loss")
                  return
               end if
            else
               nlossy = nlossy + 1
               if (da < -1.0e-13_wp .or. da > cell_a0 + 1.0e-13_wp) then
                  call test_failed(error, "area loss outside [0, category-1 area]")
                  return
               end if
            end if
            if (abs(aice(i, j) - va1) > 1.0e-15_wp) then
               call test_failed(error, "aice does not match the state")
               return
            end if
            if (abs(aice0(i, j) - max(1.0_wp - va1, 0.0_wp)) > 1.0e-15_wp) then
               call test_failed(error, "aice0 does not match the state")
               return
            end if
            ! every populated category inside its bin once the remap ran
            do n = 1, NCAT
               if (a(i, j, n) > PUNY) then
                  h = v(i, j, n)/a(i, j, n)
                  if (h < HIN(n - 1) - 1.0e-10_wp .or. h > HIN(n) + 1.0e-10_wp) nbin = nbin + 1
               end if
            end do
         end do
      end do
      write (output_unit, '(a,i0,a,i0,a,i0,a,i0,a,i0)') "itd_linear conservation nilyr=", nilyr, ": remapped ", nact, &
         ", fallback ", nfb, ", with area loss ", nlossy, ", categories out of bin ", nbin
      write (output_unit, '(a,6es10.2)') "   worst rel. closure: area(no-loss cells), ice vol, snow vol, ice energy, snow energy, salt: ", &
         worst_a, wv, ws, we, wq, wsl
      call check(error, nact > 90000, "most cells remapped")
      if (allocated(error)) return
      call check(error, ieee_is_finite(worst_a + wv + ws + we + wq + wsl), "closure ledgers finite")
      if (allocated(error)) return
      call check(error, worst_a <= CLOSE_TOL, "area conserved where no mechanism removes any")
      if (allocated(error)) return
      call check(error, wv <= CLOSE_TOL .and. ws <= CLOSE_TOL, "ice and snow volume conserved")
      if (allocated(error)) return
      call check(error, we <= CLOSE_TOL .and. wq <= CLOSE_TOL, "ice and snow energy conserved")
      if (allocated(error)) return
      call check(error, wsl <= CLOSE_TOL, "salt conserved")
      if (allocated(error)) return
      call check(error, nbin, 0, "every populated category inside its bin after a remap")
      nrb = 0
   end subroutine conserve_case

   pure integer function ntr_of(nilyr)
      integer, intent(in) :: nilyr

      ntr_of = ice_ntrcr(nilyr, 1)
   end function ntr_of

   subroutine pack_tr(i, j, nx, ny, nilyr, Tsfcn, qicen, sicen, qsnon, tr)
      integer, intent(in) :: i, j, nx, ny, nilyr
      real(wp), intent(in) :: Tsfcn(nx, ny, NCAT), qicen(nx, ny, NCAT, nilyr), sicen(nx, ny, NCAT, nilyr), qsnon(nx, ny, NCAT, 1)
      real(wp), intent(out) :: tr(ICE_NTRCR_MAX, NCAT)
      integer :: n, it

      tr = 0.0_wp
      do n = 1, NCAT
         do it = 1, ice_ntrcr(nilyr, 1)
            tr(it, n) = ice_trcr_get(it, i, j, n, nx, ny, NCAT, nilyr, 1, Tsfcn, qicen, sicen, qsnon)
         end do
      end do
   end subroutine pack_tr

   subroutine test_conserve_7(error)
      type(error_type), allocatable, intent(out) :: error

      call conserve_case(error, 7, 111111_int64)
   end subroutine test_conserve_7

   subroutine test_conserve_4(error)
      type(error_type), allocatable, intent(out) :: error

      call conserve_case(error, 4, 222222_int64)
   end subroutine test_conserve_4

   subroutine test_conserve_1(error)
      type(error_type), allocatable, intent(out) :: error

      call conserve_case(error, 1, 333333_int64)
   end subroutine test_conserve_1

   subroutine test_rebin_random(error)
      !! Rebin of states whose category means are anywhere: conserves everything, leaves every
      !! populated category inside its bin (up to the category-1 hin_max(0) floor, inert here),
      !! and the cell form is bitwise the packed form (CPU builds).
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nx = 40, ny = 25, nilyr = 7
      real(wp), allocatable :: ai0(:, :, :), vi0(:, :, :), a(:, :, :), v(:, :, :), s(:, :, :), Tsfcn(:, :, :)
      real(wp), allocatable :: qicen(:, :, :, :), sicen(:, :, :, :), qsnon(:, :, :, :)
      real(wp) :: ai0c(NCAT), vi0c(NCAT), ac(NCAT), vc(NCAT), sc(NCAT), tr(ICE_NTRCR_MAX, NCAT), trc(ICE_NTRCR_MAX, NCAT)
      real(wp) :: va0, vi0t, vs0, ei0, es0, si0, va1, vi1, vs1, ei1, es1, si1, h, w, wcell
      integer(int32) :: ierr_p, ierr_c(nx, ny)
      integer :: i, j, n, it, nbad

      allocate (ai0(nx, ny, NCAT), vi0(nx, ny, NCAT), a(nx, ny, NCAT), v(nx, ny, NCAT), s(nx, ny, NCAT), Tsfcn(nx, ny, NCAT), &
                qicen(nx, ny, NCAT, nilyr), sicen(nx, ny, NCAT, nilyr), qsnon(nx, ny, NCAT, 1))
      call fill_random(nx, ny, nilyr, 444444_int64, ai0, vi0, a, v, s, Tsfcn, qicen, sicen, qsnon, mode=2)
      ierr_c = 0_int32
      w = 0.0_wp
      nbad = 0
      do j = 1, ny
         do i = 1, nx
            call gather(i, j, nx, ny, nilyr, ai0, vi0, a, v, s, Tsfcn, qicen, sicen, qsnon, ai0c, vi0c, ac, vc, sc, tr)
            trc = tr
            call totals(nilyr, ac, vc, sc, tr(1:ice_ntrcr(nilyr, 1), :), va0, vi0t, vs0, ei0, es0, si0)
            ierr_p = 0_int32
            call itd_rebin(NCAT, ICE_NTRCR_MAX, nilyr, 1, HIN, PUNY, TF, ac, vc, sc, tr, ierr_p)
            call itd_rebin_cell(i, j, nx, ny, NCAT, nilyr, 1, HIN, PUNY, TF, a, v, s, Tsfcn, qicen, sicen, qsnon, ierr_c)
            call check(error, ierr_p == 0_int32 .and. ierr_c(i, j) == 0_int32, "no error")
            if (allocated(error)) return
            ! explicit finiteness of every output of both forms (same() and max() hide NaN)
            if (.not. (all(ieee_is_finite(ac)) .and. all(ieee_is_finite(vc)) .and. all(ieee_is_finite(sc)) .and. &
                       all(ieee_is_finite(tr(1:ice_ntrcr(nilyr, 1), :))) .and. all(ieee_is_finite(a(i, j, :))) .and. &
                       all(ieee_is_finite(v(i, j, :))) .and. all(ieee_is_finite(s(i, j, :))) .and. &
                       all(ieee_is_finite(Tsfcn(i, j, :))) .and. all(ieee_is_finite(qicen(i, j, :, :))) .and. &
                       all(ieee_is_finite(sicen(i, j, :, :))) .and. all(ieee_is_finite(qsnon(i, j, :, :))))) then
               call test_failed(error, "non-finite state after the rebin")
               return
            end if
            ! cell form == packed form
            do n = 1, NCAT
               if (.not. (same(ac(n), a(i, j, n), TOL_GOLD) .and. same(vc(n), v(i, j, n), TOL_GOLD) .and. &
                          same(sc(n), s(i, j, n), TOL_GOLD))) nbad = nbad + 1
               do it = 1, ice_ntrcr(nilyr, 1)
                  if (.not. same(tr(it, n), ice_trcr_get(it, i, j, n, nx, ny, NCAT, nilyr, 1, Tsfcn, qicen, sicen, qsnon), &
                                 TOL_GOLD)) nbad = nbad + 1
               end do
            end do
            call totals(nilyr, ac, vc, sc, tr(1:ice_ntrcr(nilyr, 1), :), va1, vi1, vs1, ei1, es1, si1)
            wcell = max(abs(va1 - va0)/max(va0, 1.0e-30_wp), abs(vi1 - vi0t)/max(vi0t, 1.0e-30_wp), &
                        abs(ei1 - ei0)/max(abs(ei0), 1.0e-30_wp), abs(si1 - si0)/max(abs(si0), 1.0e-30_wp))
            if (vs0 > 0.0_wp) wcell = max(wcell, abs(vs1 - vs0)/vs0, abs(es1 - es0)/max(abs(es0), 1.0e-30_wp))
            if (.not. ieee_is_finite(wcell)) then
               call test_failed(error, "non-finite rebin closure")
               return
            end if
            w = max(w, wcell)
            do n = 1, NCAT
               if (ac(n) > PUNY) then
                  h = vc(n)/ac(n)
                  if (h < HIN(n - 1) - 1.0e-12_wp .or. h > HIN(n) + 1.0e-12_wp) nbad = nbad + 1000
               end if
            end do
         end do
      end do
      call check(error, w <= CLOSE_TOL, "rebin conserves area, volumes, energy and salt")
      if (allocated(error)) return
      call check(error, nbad, 0, "cell form equals packed form; every populated category inside its bin")
   end subroutine test_rebin_random

   ! ------------------------------------------------------------------
   ! oracle
   ! ------------------------------------------------------------------

   subroutine oracle_run(error, kset, gold, goldr)
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: kset
      real(wp), intent(in) :: gold(NOUT, NC), goldr(NOUT, NC)
      integer :: c, nilyr, ntr, n, it, nfb, nup, ndown, nwhole, ncat1, nhim, ncheck, bad, nrbchg
      real(wp) :: worst_o, worst_r, trb(ICE_NTRCR_MAX, NCAT), ar(NCAT), vr(NCAT), sr(NCAT), gotr(15 + ICE_NTRCR_MAX*NCAT)
      integer(int32) :: ierr_r
      real(wp) :: ai0(NCAT), vi0(NCAT), a(NCAT), v(NCAT), s(NCAT), a_in(NCAT), v_in(NCAT)
      real(wp), allocatable :: trcrn(:, :), tr_in(:, :), got(:), tr_nohm(:, :)
      real(wp) :: aice, aice0, chk, a_nohm(NCAT), v_nohm(NCAT), s_nohm(NCAT), ai_nohm, a0_nohm, hicen(NCAT), dai(NCAT), dvi(NCAT)
      real(wp) :: a_pl(NCAT)
      integer :: donor(NCAT)
      logical :: remap
      integer(int32) :: ierr, ierr_p
      character(len=200) :: msg

      call fill_golden()
      nfb = 0; nup = 0; ndown = 0; nwhole = 0; ncat1 = 0; nhim = 0; nrbchg = 0
      worst_o = 0.0_wp; worst_r = 0.0_wp
      do c = 1, NC
         nilyr = GNIL(c)
         ntr = ice_ntrcr(nilyr, 1)
         allocate (trcrn(ntr, NCAT), tr_in(ntr, NCAT), tr_nohm(ntr, NCAT), got(17 + ntr*NCAT))
         ai0 = GIN(1:5, c); vi0 = GIN(6:10, c); a = GIN(11:15, c); v = GIN(16:20, c); s = GIN(21:25, c)
         call lcg_tracers(c, nilyr, ntr, trcrn)
         ! input regeneration check (the LCG chain is exact; the formulas differ only under FMA)
         chk = 0.0_wp
         do n = 1, NCAT
            do it = 1, ntr
               chk = chk + trcrn(it, n)*real(it + 10*n, wp)*1.0e-8_wp
            end do
         end do
         call check(error, abs(chk - GCHK(c)) <= 1.0e-12_wp*abs(GCHK(c)), "tracer input checksum")
         if (allocated(error)) return
         tr_in = trcrn; a_in = a; v_in = v
         ierr = 0_int32
         call call_packed(kset, nilyr, ai0, vi0, a, v, s, trcrn, .false., aice, aice0, ierr)
         call check(error, iand(ierr, not(ICE_ERR_REMAP_FALLBACK)) == 0_int32, "oracle case: no error bit but the fallback")
         if (allocated(error)) return
         got(1:5) = a; got(6:10) = v; got(11:15) = s; got(16) = aice; got(17) = aice0
         got(18:) = reshape(trcrn, [ntr*NCAT])
         bad = ice_golden_first_bad(got, gold(1:17 + ntr*NCAT, c), TOL_GOLD, TOL_GOLD)
         if (bad > 0) then
            write (msg, '(a,i0,a,i0,a,es24.16,a,es24.16)') "oracle case ", c, " output ", bad, " got ", got(bad), &
               " expected ", gold(bad, c)
            call test_failed(error, trim(msg))
            return
         end if
         do it = 1, 17 + ntr*NCAT
            worst_o = max(worst_o, abs(got(it) - gold(it, c))/max(abs(gold(it, c)), 1.0_wp))
         end do
         ! upstream rebin applied to the linear_itd output (the fallback-ON path)
         ar = a; vr = v; sr = s
         trb = 0.0_wp
         trb(1:ntr, :) = trcrn
         ierr_r = 0_int32
         call itd_rebin(NCAT, ICE_NTRCR_MAX, nilyr, 1, HIN, PUNY, TF, ar, vr, sr, trb, ierr_r)
         call check(error, ierr_r == 0_int32, "rebin clean")
         if (allocated(error)) return
         gotr(1:5) = ar; gotr(6:10) = vr; gotr(11:15) = sr
         gotr(16:15 + ntr*NCAT) = reshape(trb(1:ntr, :), [ntr*NCAT])
         bad = ice_golden_first_bad(gotr(1:15 + ntr*NCAT), goldr(1:15 + ntr*NCAT, c), TOL_GOLD, TOL_GOLD)
         if (bad > 0) then
            write (msg, '(a,i0,a,i0,a,es24.16,a,es24.16)') "oracle rebin case ", c, " output ", bad, " got ", gotr(bad), &
               " expected ", goldr(bad, c)
            call test_failed(error, trim(msg))
            return
         end if
         do it = 1, 15 + ntr*NCAT
            worst_r = max(worst_r, abs(gotr(it) - goldr(it, c))/max(abs(goldr(it, c)), 1.0_wp))
         end do
         if (any(ar /= a) .or. any(vr /= v)) nrbchg = nrbchg + 1
         ! branch coverage, from the port's plan on the same inputs
         a_pl = a_in
         ierr_p = 0_int32
         call itd_linear_plan(NCAT, HIN, PUNY, ai0, vi0, a_pl, v_in, hicen, remap, donor, dai, dvi, ierr_p)
         if (.not. remap) then
            nfb = nfb + 1
         else
            if (a_pl(1) < a_in(1)) ncat1 = ncat1 + 1
            do n = 1, NCAT - 1
               if (donor(n) == n .and. dai(n) > 0.0_wp) nup = nup + 1
               if (donor(n) == n + 1 .and. dai(n) > 0.0_wp) ndown = ndown + 1
               if (donor(n) > 0 .and. dai(n) > 0.0_wp) then
                  if (dai(n) == a_pl(donor(n))) nwhole = nwhole + 1
               end if
            end do
            ! hi_min active iff switching it off changes category 1
            a_nohm = a_in; v_nohm = v_in; s_nohm = GIN(21:25, c); tr_nohm = tr_in
            ierr_p = 0_int32
            call call_packed(kset, nilyr, ai0, vi0, a_nohm, v_nohm, s_nohm, tr_nohm, .false., ai_nohm, a0_nohm, ierr_p, 0.0_wp)
            if (a_nohm(1) /= a(1)) nhim = nhim + 1
         end if
         deallocate (trcrn, tr_in, tr_nohm, got)
      end do
      ncheck = nfb
      write (output_unit, '(a,es10.2,a,es10.2)') "oracle worst scaled difference: linear_itd ", worst_o, ", rebin ", worst_r
      call check(error, ieee_is_finite(worst_o + worst_r), "oracle comparison finite")
      if (allocated(error)) return
      call check(error, nrbchg >= 3, "the oracle rebin rows include states the rebin changes (fallback ON)")
      if (allocated(error)) return
      write (output_unit, '(a,6i4)') "oracle coverage (fallback, up, down, whole, cat1 loss, hi_min): ", nfb, nup, ndown, &
         nwhole, ncat1, nhim
      call check(error, nfb >= 3, "oracle set covers the fallback (veto types)")
      if (allocated(error)) return
      call check(error, nup >= 4 .and. ndown >= 4, "oracle set covers both transfer directions")
      if (allocated(error)) return
      call check(error, nwhole >= 1, "oracle set covers a whole-category shift")
      if (allocated(error)) return
      call check(error, ncat1 >= 2, "oracle set covers the category-1 melt-area loss")
      if (allocated(error)) return
      call check(error, nhim >= 1, "oracle set covers the hi_min floor")
      ncheck = ncheck + 0
   end subroutine oracle_run

   subroutine test_oracle_icepack(error)
      type(error_type), allocatable, intent(out) :: error

      call fill_golden()
      call oracle_run(error, ICE_CONSTANTS_ICEPACK, GOI, GRI)
   end subroutine test_oracle_icepack

   subroutine test_oracle_cesm(error)
      type(error_type), allocatable, intent(out) :: error

      call fill_golden()
      call oracle_run(error, ICE_CONSTANTS_CESM, GOC, GRC)
   end subroutine test_oracle_cesm

   ! ------------------------------------------------------------------
   ! cell entry vs packed, device
   ! ------------------------------------------------------------------

   subroutine test_cell_vs_packed(error)
      !! `ice_linear_itd_cell` against the packed routine over random cells, with and without the
      !! rebin fallback: every output, bitwise on CPU builds.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nx = 60, ny = 40, nilyr = 7
      real(wp), allocatable :: ai0(:, :, :), vi0(:, :, :), a(:, :, :), v(:, :, :), s(:, :, :), Tsfcn(:, :, :)
      real(wp), allocatable :: qicen(:, :, :, :), sicen(:, :, :, :), qsnon(:, :, :, :)
      real(wp) :: ai0c(NCAT), vi0c(NCAT), ac(NCAT), vc(NCAT), sc(NCAT), tr(ICE_NTRCR_MAX, NCAT)
      real(wp) :: aice_p, aice0_p, aice_c, aice0_c
      type(ice_params_t) :: p
      integer(int32) :: ierr_p, ierr_c(nx, ny)
      integer :: i, j, n, it, pass, nbad, nfb, nup, ndown, nboth, nwhole, ncat1
      logical :: fb, remap
      real(wp) :: a_pl(NCAT), hicen(NCAT), dai(NCAT), dvi(NCAT)
      integer :: donor(NCAT)
      integer(int32) :: ierr_t

      call params(ICE_CONSTANTS_ICEPACK, nilyr, p)
      allocate (ai0(nx, ny, NCAT), vi0(nx, ny, NCAT), a(nx, ny, NCAT), v(nx, ny, NCAT), s(nx, ny, NCAT), Tsfcn(nx, ny, NCAT), &
                qicen(nx, ny, NCAT, nilyr), sicen(nx, ny, NCAT, nilyr), qsnon(nx, ny, NCAT, 1))
      nbad = 0
      nfb = 0; nup = 0; ndown = 0; nboth = 0; nwhole = 0; ncat1 = 0
      do pass = 1, 3
         fb = pass >= 2
         call fill_random(nx, ny, nilyr, 555555_int64 + int(pass, int64), ai0, vi0, a, v, s, Tsfcn, qicen, sicen, qsnon, &
                          mode=merge(2, 1, pass == 3))
         ierr_c = 0_int32
         do j = 1, ny
            do i = 1, nx
               call gather(i, j, nx, ny, nilyr, ai0, vi0, a, v, s, Tsfcn, qicen, sicen, qsnon, ai0c, vi0c, ac, vc, sc, tr)
               ! branch census from the port's plan on the same (pre-remap) inputs
               a_pl = ac
               ierr_t = 0_int32
               call itd_linear_plan(NCAT, HIN, PUNY, ai0c, vi0c, a_pl, vc, hicen, remap, donor, dai, dvi, ierr_t)
               if (remap .and. pass == 1) then
                  if (a_pl(1) < ac(1)) ncat1 = ncat1 + 1
                  do n = 1, NCAT - 1
                     if (donor(n) == n .and. dai(n) > 0.0_wp) nup = nup + 1
                     if (donor(n) == n + 1 .and. dai(n) > 0.0_wp) ndown = ndown + 1
                     if (donor(n) > 0 .and. dai(n) > 0.0_wp) then
                        if (dai(n) == a_pl(donor(n))) nwhole = nwhole + 1
                     end if
                  end do
                  do n = 2, NCAT - 1   ! category n donating across BOTH of its boundaries
                     if (donor(n) == n .and. dai(n) > 0.0_wp .and. donor(n - 1) == n .and. dai(n - 1) > 0.0_wp) nboth = nboth + 1
                  end do
               end if
               ierr_p = 0_int32
               call itd_linear_itd(NCAT, ICE_NTRCR_MAX, nilyr, 1, HIN, ai0c, vi0c, ac, vc, sc, tr, TF, p%c%rhos, p%c%Lfresh, &
                                   p%thermo%hi_min, PUNY, fb, aice_p, aice0_p, ierr_p)
               call ice_linear_itd_cell(i, j, nx, ny, NCAT, nilyr, 1, HIN, ai0, vi0, a, v, s, Tsfcn, qicen, sicen, qsnon, &
                                        TF, p%c%rhos, p%c%Lfresh, p%thermo%hi_min, PUNY, fb, aice_c, aice0_c, ierr_c)
               if (ierr_p /= ierr_c(i, j)) nbad = nbad + 1
               if (iand(ierr_p, ICE_ERR_REMAP_FALLBACK) /= 0_int32) nfb = nfb + 1
               if (.not. (same(aice_p, aice_c, TOL_GOLD) .and. same(aice0_p, aice0_c, TOL_GOLD))) nbad = nbad + 1
               do n = 1, NCAT
                  if (.not. (same(ac(n), a(i, j, n), TOL_GOLD) .and. same(vc(n), v(i, j, n), TOL_GOLD) .and. &
                             same(sc(n), s(i, j, n), TOL_GOLD))) nbad = nbad + 1
                  do it = 1, ice_ntrcr(nilyr, 1)
                     if (.not. same(tr(it, n), ice_trcr_get(it, i, j, n, nx, ny, NCAT, nilyr, 1, Tsfcn, qicen, sicen, qsnon), &
                                    TOL_GOLD)) nbad = nbad + 1
                  end do
               end do
               if (.not. (all(ieee_is_finite(ac)) .and. all(ieee_is_finite(vc)) .and. all(ieee_is_finite(sc)) .and. &
                          all(ieee_is_finite(tr(1:ice_ntrcr(nilyr, 1), :))) .and. ieee_is_finite(aice_p) .and. &
                          all(ieee_is_finite(a(i, j, :))) .and. all(ieee_is_finite(qicen(i, j, :, :))) .and. &
                          all(ieee_is_finite(sicen(i, j, :, :))) .and. all(ieee_is_finite(qsnon(i, j, :, :))) .and. &
                          ieee_is_finite(aice_c))) nbad = nbad + 1000
            end do
         end do
      end do
      write (output_unit, '(a,6i6)') "cell-vs-packed census (fallback, up, down, both-bounds, whole, cat1 loss): ", nfb, nup, &
         ndown, nboth, nwhole, ncat1
      call check(error, nfb > 0, "random cells include fallbacks")
      if (allocated(error)) return
      call check(error, nup > 0 .and. ndown > 0, "transfers in both directions")
      if (allocated(error)) return
      call check(error, nboth > 0, "a category donating across both of its boundaries")
      if (allocated(error)) return
      call check(error, nwhole > 0, "a whole-category shift")
      if (allocated(error)) return
      call check(error, ncat1 > 0, "the category-1 melt loss")
      if (allocated(error)) return
      call check(error, nbad, 0, "cell entry equals the packed routine (all outputs, all cells)")
   end subroutine test_cell_vs_packed

   subroutine test_device(error)
      !! The host shim `ice_linear_itd_step` (one `do concurrent` over the cells, the cell kernel
      !! inside) over random 64 x 48 states, EVERY output of EVERY cell against the host packed
      !! reference: the plain remap, the remap with the rebin fallback ON (the device rebin path),
      !! and rebin food (`mode = 2`, category means anywhere); then 512 x 512 x 5 x 7 throughput.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p

      call params(ICE_CONSTANTS_ICEPACK, 7, p)
      call device_case(error, "remap", 666666_int64, 1, .false.)
      if (allocated(error)) return
      call device_case(error, "remap + rebin fallback", 888888_int64, 1, .true.)
      if (allocated(error)) return
      call device_case(error, "rebin food + fallback", 999999_int64, 2, .true.)
      if (allocated(error)) return
      call throughput(error, p)
   end subroutine test_device

   subroutine device_case(error, label, seed, mode, rebin_fb)
      !! One device sweep against the host reference.  Measured worst scaled differences on the
      !! V100 (nvfortran 26.5, FMA on) are recorded in the comment of `TOL_DEV` / `TOL_TRC`.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: label
      integer(int64), intent(in) :: seed
      integer, intent(in) :: mode
      logical, intent(in) :: rebin_fb
      integer, parameter :: nx = 64, ny = 48, nilyr = 7
      real(wp), allocatable :: ai0(:, :, :), vi0(:, :, :), a(:, :, :), v(:, :, :), s(:, :, :), Tsfcn(:, :, :)
      real(wp), allocatable :: qicen(:, :, :, :), sicen(:, :, :, :), qsnon(:, :, :, :), Tf2(:, :), aice(:, :), aice0(:, :)
      integer(int32), allocatable :: ierr(:, :)
      real(wp) :: hin_loc(0:NCAT), w_main, w_hl, w_dr
      type(ice_params_t) :: p
      integer :: nbad, ndr

      call params(ICE_CONSTANTS_ICEPACK, nilyr, p)
      allocate (ai0(nx, ny, NCAT), vi0(nx, ny, NCAT), a(nx, ny, NCAT), v(nx, ny, NCAT), s(nx, ny, NCAT), Tsfcn(nx, ny, NCAT), &
                qicen(nx, ny, NCAT, nilyr), sicen(nx, ny, NCAT, nilyr), qsnon(nx, ny, NCAT, 1), Tf2(nx, ny), aice(nx, ny), &
                aice0(nx, ny), ierr(nx, ny))
      call fill_random(nx, ny, nilyr, seed, ai0, vi0, a, v, s, Tsfcn, qicen, sicen, qsnon, mode=mode)
      hin_loc = HIN
      Tf2 = TF
      aice = 0.0_wp; aice0 = 0.0_wp; ierr = 0_int32
      !$acc enter data copyin(hin_loc, Tf2, ai0, vi0, a, v, s, Tsfcn, qicen, sicen, qsnon, ierr, aice, aice0)
      call ice_linear_itd_step(nx, ny, NCAT, nilyr, 1, 1, nx, 1, ny, p, hin_loc, Tf2, ai0, vi0, a, v, s, Tsfcn, qicen, sicen, &
                               qsnon, rebin_fb, aice, aice0, ierr)
      !$acc update self(a, v, s, Tsfcn, qicen, sicen, qsnon, ierr, aice, aice0)
      !$acc exit data delete(hin_loc, Tf2, ai0, vi0, a, v, s, Tsfcn, qicen, sicen, qsnon, ierr, aice, aice0)
      ! the reference regenerates the INPUT state
      call compare_with_host(nx, ny, nilyr, seed, mode, rebin_fb, p, a, v, s, Tsfcn, qicen, sicen, qsnon, aice, aice0, ierr, &
                             nbad, w_main, w_hl, w_dr, ndr)
      write (output_unit, '(a,a,a,i0,a,3es10.2,a,i0)') "linear_itd device sweep [", label, "]: mismatches ", nbad, &
         "; worst scaled diff area/volume/aggregates, healthy tracers, drained tracers: ", w_main, w_hl, w_dr, &
         "; drained categories ", ndr
      call check(error, ieee_is_finite(w_main + w_hl + w_dr), "device comparison finite")
      if (allocated(error)) return
      call check(error, nbad, 0, "every output of every cell equals the host reference")
   end subroutine device_case

   subroutine compare_with_host(nx, ny, nilyr, seed, mode, rebin_fb, p, a, v, s, Tsfcn, qicen, sicen, qsnon, aice, aice0, ierr, &
                                nbad, w_main, w_hl, w_dr, ndr)
      !! Regenerate the input state and run the packed host reference cell by cell against the
      !! shim's outputs.  Tracers of a category that is nearly drained on the host (area <= puny
      !! or volume < 1e-6 of the cell volume) divide a small difference of two large products and
      !! take the loose `TOL_TRC` (also a category that lost more than 99.9 % of its entry area or volume); every other output takes `TOL_DEV`.
      integer, intent(in) :: nx, ny, nilyr, mode
      integer(int64), intent(in) :: seed
      logical, intent(in) :: rebin_fb
      type(ice_params_t), intent(in) :: p
      real(wp), intent(in) :: a(nx, ny, NCAT), v(nx, ny, NCAT), s(nx, ny, NCAT), Tsfcn(nx, ny, NCAT)
      real(wp), intent(in) :: qicen(nx, ny, NCAT, nilyr), sicen(nx, ny, NCAT, nilyr), qsnon(nx, ny, NCAT, 1)
      real(wp), intent(in) :: aice(nx, ny), aice0(nx, ny)
      integer(int32), intent(in) :: ierr(nx, ny)
      integer, intent(out) :: nbad, ndr
      real(wp), intent(out) :: w_main, w_hl, w_dr
      real(wp), allocatable :: ai0(:, :, :), vi0(:, :, :), a0(:, :, :), v0(:, :, :), s0(:, :, :), T0(:, :, :)
      real(wp), allocatable :: q0(:, :, :, :), si0(:, :, :, :), qs0(:, :, :, :)
      real(wp) :: ai0c(NCAT), vi0c(NCAT), ac(NCAT), vc(NCAT), sc(NCAT), tr(ICE_NTRCR_MAX, NCAT), aice_p, aice0_p, tot
      real(wp) :: ain(NCAT), vin(NCAT)
      integer(int32) :: ierr_p
      integer :: i, j, n, it
      logical :: drained

      allocate (ai0(nx, ny, NCAT), vi0(nx, ny, NCAT), a0(nx, ny, NCAT), v0(nx, ny, NCAT), s0(nx, ny, NCAT), T0(nx, ny, NCAT), &
                q0(nx, ny, NCAT, nilyr), si0(nx, ny, NCAT, nilyr), qs0(nx, ny, NCAT, 1))
      call fill_random(nx, ny, nilyr, seed, ai0, vi0, a0, v0, s0, T0, q0, si0, qs0, mode=mode)
      nbad = 0
      ndr = 0
      w_main = 0.0_wp
      w_hl = 0.0_wp
      w_dr = 0.0_wp
      do j = 1, ny
         do i = 1, nx
            tot = sum(a0(i, j, :))
            call gather(i, j, nx, ny, nilyr, ai0, vi0, a0, v0, s0, T0, q0, si0, qs0, ai0c, vi0c, ac, vc, sc, tr)
            ain = ac
            vin = vc
            ierr_p = 0_int32
            if (tot > PUNY) then
               call itd_linear_itd(NCAT, ICE_NTRCR_MAX, nilyr, 1, HIN, ai0c, vi0c, ac, vc, sc, tr, TF, p%c%rhos, p%c%Lfresh, &
                                   p%thermo%hi_min, PUNY, rebin_fb, aice_p, aice0_p, ierr_p)
            else
               aice_p = tot
               aice0_p = max(1.0_wp - tot, 0.0_wp)
            end if
            if (ierr_p /= ierr(i, j)) nbad = nbad + 1
            call cmp(aice_p, aice(i, j), TOL_DEV, w_main)
            call cmp(aice0_p, aice0(i, j), TOL_DEV, w_main)
            do n = 1, NCAT
               call cmp(ac(n), a(i, j, n), TOL_DEV, w_main)
               call cmp(vc(n), v(i, j, n), TOL_DEV, w_main)
               call cmp(sc(n), s(i, j, n), TOL_DEV, w_main)
               drained = ac(n) <= PUNY .or. vc(n) < 1.0e-6_wp*sum(vc) .or. ac(n) < 1.0e-3_wp*ain(n) .or. vc(n) < 1.0e-3_wp*vin(n)
               if (drained) ndr = ndr + 1
               do it = 1, ice_ntrcr(nilyr, 1)
                  if (drained) then
                     call cmp(tr(it, n), ice_trcr_get(it, i, j, n, nx, ny, NCAT, nilyr, 1, Tsfcn, qicen, sicen, qsnon), &
                              TOL_TRC, w_dr)
                  else
                     call cmp(tr(it, n), ice_trcr_get(it, i, j, n, nx, ny, NCAT, nilyr, 1, Tsfcn, qicen, sicen, qsnon), &
                              TOL_DEV, w_hl)
                  end if
               end do
            end do
         end do
      end do
   contains
      subroutine cmp(ref, got, tl, worst_l)
         real(wp), intent(in) :: ref, got, tl
         real(wp), intent(inout) :: worst_l
         real(wp) :: dd

         if (.not. (ieee_is_finite(ref) .and. ieee_is_finite(got))) then
            nbad = nbad + 1
            return
         end if
         dd = abs(ref - got)/max(abs(ref), abs(got), 1.0_wp)
         worst_l = max(worst_l, dd)
         if (dd > tl) then
            nbad = nbad + 1
            if (nbad <= 5) write (output_unit, '(a,2es24.16)') "   mismatch ref/got ", ref, got
         end if
      end subroutine cmp
   end subroutine compare_with_host

   subroutine throughput(error, p)
      !! 512 x 512 cells x 5 categories x 7 layers, repeated sweeps of the shim on re-seeded data.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t), intent(in) :: p
      integer, parameter :: nx = 512, ny = 512, nilyr = 7, nrep = 5
      real(wp), allocatable :: ai0(:, :, :), vi0(:, :, :), a(:, :, :), v(:, :, :), s(:, :, :), Tsfcn(:, :, :)
      real(wp), allocatable :: qicen(:, :, :, :), sicen(:, :, :, :), qsnon(:, :, :, :), Tf2(:, :), aice(:, :), aice0(:, :)
      real(wp), allocatable :: a_s(:, :, :), v_s(:, :, :), s_s(:, :, :), T_s(:, :, :), q_s(:, :, :, :), si_s(:, :, :, :)
      real(wp), allocatable :: qs_s(:, :, :, :)
      integer(int32), allocatable :: ierr(:, :)
      real(wp) :: hin_loc(0:NCAT), secs
      integer(int64) :: t0, t1, rate
      integer :: k

      allocate (ai0(nx, ny, NCAT), vi0(nx, ny, NCAT), a(nx, ny, NCAT), v(nx, ny, NCAT), s(nx, ny, NCAT), Tsfcn(nx, ny, NCAT), &
                qicen(nx, ny, NCAT, nilyr), sicen(nx, ny, NCAT, nilyr), qsnon(nx, ny, NCAT, 1), Tf2(nx, ny), aice(nx, ny), &
                aice0(nx, ny), ierr(nx, ny))
      allocate (a_s(nx, ny, NCAT), v_s(nx, ny, NCAT), s_s(nx, ny, NCAT), T_s(nx, ny, NCAT), q_s(nx, ny, NCAT, nilyr), &
                si_s(nx, ny, NCAT, nilyr), qs_s(nx, ny, NCAT, 1))
      call fill_random(nx, ny, nilyr, 777777_int64, ai0, vi0, a, v, s, Tsfcn, qicen, sicen, qsnon)
      a_s = a; v_s = v; s_s = s; T_s = Tsfcn; q_s = qicen; si_s = sicen; qs_s = qsnon
      hin_loc = HIN
      Tf2 = TF
      aice = 0.0_wp; aice0 = 0.0_wp; ierr = 0_int32
      !$acc enter data copyin(hin_loc, Tf2, ai0, vi0, a, v, s, Tsfcn, qicen, sicen, qsnon, ierr, aice, aice0, a_s, v_s, s_s, T_s, q_s, si_s, qs_s)
      secs = 0.0_wp
      call system_clock(count_rate=rate)
      do k = 0, nrep
         ! restore the input (device-side copy), time the sweep only
         !$acc kernels present(a, v, s, Tsfcn, qicen, sicen, qsnon, a_s, v_s, s_s, T_s, q_s, si_s, qs_s)
         a = a_s; v = v_s; s = s_s; Tsfcn = T_s; qicen = q_s; sicen = si_s; qsnon = qs_s
         !$acc end kernels
         call system_clock(t0)
         call ice_linear_itd_step(nx, ny, NCAT, nilyr, 1, 1, nx, 1, ny, p, hin_loc, Tf2, ai0, vi0, a, v, s, Tsfcn, qicen, sicen, &
                                  qsnon, .false., aice, aice0, ierr)
         !$acc wait
         call system_clock(t1)
         if (k > 0) secs = secs + real(t1 - t0, wp)/real(rate, wp)
      end do
      !$acc update self(ierr, aice, aice0)
      !$acc exit data delete(hin_loc, Tf2, ai0, vi0, a, v, s, Tsfcn, qicen, sicen, qsnon, ierr, aice, aice0, a_s, v_s, s_s, T_s, q_s, si_s, qs_s)
      write (output_unit, '(a,f10.3,a,es11.3,a)') "linear_itd sweep 512x512x5x7: ", secs/real(nrep, wp)*1.0e3_wp, " ms/sweep, ", &
         real(nx*ny, wp)*real(nrep, wp)/secs, " cells/s"
      call check(error, all(iand(ierr, not(ICE_ERR_REMAP_FALLBACK)) == 0_int32), "no error bit but the counted fallback")
      if (allocated(error)) return
      call check(error, all(ieee_is_finite(aice)) .and. all(ieee_is_finite(aice0)), "finite aggregates")
   end subroutine throughput

end module test_cice_itd_linear
