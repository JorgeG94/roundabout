!! New ice (D4.2, `add_new_ice`, BL99 branch): closed forms (frazil volume,
!! open-water fill, overflow into existing ice, the `hi0max` cap and the
!! `hfrazilmin` floor), conservation on random cells, the no-op / dry / error
!! contracts, goldens from the real Icepack routine, and the device sweep.
!!
!! Plan tests: T1 closure_exact, T2 thin_lead, T3 thick_full_lead, T4
!! partial_lead_overflow, T5 closed_ice, T6 hfrazilmin_sweep, T7 salinity, T8
!! bl99_inversion, T9 update_ocn_f_off, T10 no_bank, T11 dry_cell, T14
!! err_flags / NaN, T15 cell-vs-column + device-vs-host (T12/T13 belong to the
!! therm2 driver PR).  Expected values are closed forms written here or rows of
!! the Icepack oracle (upstream `add_new_ice`, inputs and outputs both embedded).
module test_cice_new_ice
   use, intrinsic :: iso_fortran_env, only: int32, int64, output_unit
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use rdb_constants, only: wp
   use rdb_ice_errors, only: ICE_ERR_NONFINITE, ICE_ERR_NEW_ICE_AREA, ICE_ERR_TBOT, ICE_ERR_REMAP_FALLBACK
   use rdb_ice_params, only: ice_params_t, ice_params_default, ICE_CONSTANTS_ICEPACK, ICE_CONSTANTS_CESM, &
                             ice_salinity_profile, ice_T_from_q
   use rdb_ice_new_ice, only: ice_add_new_ice_cell, ice_new_ice_step
   use lcg_deterministic, only: lcg_next, lcg_unit
   use test_cice_util, only: ice_golden_check
   implicit none
   private

   public :: collect_cice_new_ice_tests

   integer, parameter :: NC = 5
   real(wp), parameter :: DT = 900.0_wp
   real(wp), parameter :: HIN1 = 0.6445072168194258_wp    ! hin_max(1), kcatbound 0, ncat 5 (upstream value)
   real(wp), parameter :: TF0 = -1.8_wp
   integer, parameter :: NGOLD = 80, NSET1 = 40

   type :: cell_t
      !! One cell as 1 x 1 SoA arrays (so the kernel is called exactly as in the driver).
      integer :: nil = 4
      real(wp), allocatable :: a(:, :, :), v(:, :, :), ts(:, :, :), q(:, :, :, :), s(:, :, :, :)
      real(wp), allocatable :: a0(:, :), fr(:, :), fs(:, :), fz(:, :), sal(:)
      integer(int32), allocatable :: ie(:, :)
   end type cell_t

contains

   subroutine collect_cice_new_ice_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("t02_thin_lead", test_thin_lead), &
                  new_unittest("t03_thick_full_lead", test_thick_full_lead), &
                  new_unittest("t04_partial_lead_overflow", test_overflow), &
                  new_unittest("t05_closed_ice", test_closed_ice), &
                  new_unittest("t06_hfrazilmin_sweep", test_hfm_sweep), &
                  new_unittest("t08_bl99_inversion", test_bl99_inversion), &
                  new_unittest("t09_update_ocn_f_off", test_updf_off), &
                  new_unittest("t10_no_bank_bitwise_noop", test_no_bank), &
                  new_unittest("t11_dry_cell_keeps_bank", test_dry_cell), &
                  new_unittest("t14_err_flags_and_nan", test_errors), &
                  new_unittest("t14_stale_aice0_ignored", test_stale_aice0), &
                  new_unittest("t14_driver_error_cells_keep_bank", test_driver_errors), &
                  new_unittest("t01_t07_closure_exact_random_nilyr_1_4_7", test_closure_random), &
                  new_unittest("oracle_new_ice", test_oracle), &
                  new_unittest("t15_cell_entry_vs_column", test_cell_vs_column), &
                  new_unittest("t15_driver_every_cell_vs_host_interior", test_driver_every_cell), &
                  new_unittest("driver_throughput", test_throughput) &
                  ]
   end subroutine collect_cice_new_ice_tests

   ! ------------------------------------------------------------------ helpers

   subroutine cell_init(c, nil)
      type(cell_t), intent(out) :: c
      integer, intent(in) :: nil
      real(wp) :: sp(nil + 1)

      c%nil = nil
      allocate (c%a(1, 1, NC), c%v(1, 1, NC), c%ts(1, 1, NC), source=0.0_wp)
      allocate (c%q(1, 1, NC, nil), c%s(1, 1, NC, nil), source=0.0_wp)
      allocate (c%a0(1, 1), source=1.0_wp)
      allocate (c%fr(1, 1), c%fs(1, 1), c%fz(1, 1), source=0.0_wp)
      allocate (c%ie(1, 1), source=0_int32)
      allocate (c%sal(nil))
      call ice_salinity_profile(nil, 3.2_wp, 0.1_wp, sp)
      c%sal = sp(1:nil)
   end subroutine cell_init

   subroutine cell_fill_layers(c, seed)
      !! Layer enthalpy in [-4e8, -1e8], salinity in [0, 6] (deterministic).
      type(cell_t), intent(inout) :: c
      integer(int64), intent(in) :: seed
      integer(int64) :: st
      integer :: n, k

      st = lcg_next(seed)
      do n = 1, NC
         st = lcg_next(st); c%ts(1, 1, n) = -0.1_wp - 30.0_wp*lcg_unit(st)
         do k = 1, c%nil
            st = lcg_next(st); c%q(1, 1, n, k) = -4.0e8_wp + 3.0e8_wp*lcg_unit(st)
            st = lcg_next(st); c%s(1, 1, n, k) = 6.0_wp*lcg_unit(st)
         end do
      end do
   end subroutine cell_fill_layers

   subroutine cell_set_a0(c)
      !! aice0 = max(1 - sum(aicen), 0) with the sequential sum upstream uses.
      type(cell_t), intent(inout) :: c
      real(wp) :: sm
      integer :: n

      sm = 0.0_wp
      do n = 1, NC
         sm = sm + c%a(1, 1, n)
      end do
      c%a0(1, 1) = max(1.0_wp - sm, 0.0_wp)
   end subroutine cell_set_a0

   function volume_of(rhoi, Lfresh, frz) result(vol)
      !! Closed form: new-ice volume per cell area of a freezing potential [W m-2] over DT.
      real(wp), intent(in) :: rhoi, Lfresh, frz
      real(wp) :: vol
      vol = frz*DT/(rhoi*Lfresh)
   end function volume_of

   function frz_of(rhoi, Lfresh, vol) result(frz)
      real(wp), intent(in) :: rhoi, Lfresh, vol
      real(wp) :: frz
      frz = vol*rhoi*Lfresh/DT
   end function frz_of

   subroutine test_stale_aice0(error)
      !! `aice0` is recomputed from `aicen`: a stale (too large, NaN, negative) input gives
      !! exactly the result of a consistent one, in both the fill and the overflow regimes.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      type(cell_t) :: c, ref
      real(wp), parameter :: VOLS(2) = [0.1_wp, 0.9_wp]
      real(wp) :: stale(3)
      integer :: iv, m

      p = params(ICE_CONSTANTS_ICEPACK, 4, .true., 5.0_wp, 0.05_wp)
      stale = [0.95_wp, ieee_value(1.0_wp, ieee_quiet_nan), -3.0_wp]
      do iv = 1, 2
         do m = 1, 3
            call cell_init(ref, 4)
            call cell_fill_layers(ref, 60_int64)
            ref%a(1, 1, 2) = 0.3_wp; ref%v(1, 1, 2) = 0.4_wp; call cell_set_a0(ref)
            c = ref
            c%a0(1, 1) = stale(m)
            call cell_run(ref, p, frz_of(p%c%rhoi, p%c%Lfresh, VOLS(iv)), TF0)
            call cell_run(c, p, frz_of(p%c%rhoi, p%c%Lfresh, VOLS(iv)), TF0)
            call check(error, cell_equal(c, ref) .and. c%ie(1, 1) == 0_int32 .and. &
                       sum(c%a(1, 1, :)) + c%a0(1, 1) <= 1.0_wp + 1.0e-14_wp, "stale aice0 ignored, area <= 1")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_stale_aice0

   subroutine test_driver_errors(error)
      !! Driver: a wet freezing cell with a NaN layer enthalpy keeps its bank and gets the bit; a cell
      !! that already carries a fatal bit keeps its bank AND its state (it is not processed); a
      !! non-fatal (counted) bit does not block; a clean cell spends its bank.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      integer, parameter :: NX = 4, NIL = 4
      real(wp) :: a(NX, 1, NC), v(NX, 1, NC), ts(NX, 1, NC), q(NX, 1, NC, NIL), s(NX, 1, NC, NIL), sal(NIL), sp(NIL + 1)
      real(wp) :: a0(NX, 1), fr(NX, 1), fs(NX, 1), fz(NX, 1), frz(NX, 1), bank(NX, 1), tf(NX, 1), hv(0:NC), vol
      real(wp) :: v_in(NX, 1, NC)
      integer(int32) :: ie(NX, 1)
      logical :: wet(NX, 1)

      p = params(ICE_CONSTANTS_ICEPACK, NIL, .true., 5.0_wp, 0.05_wp)
      call ice_salinity_profile(NIL, 3.2_wp, 0.1_wp, sp); sal = sp(1:NIL)
      a = 0.0_wp; v = 0.0_wp; ts = -5.0_wp; q = -3.0e8_wp; s = 3.0_wp; a0 = 1.0_wp
      a(:, 1, 2) = 0.3_wp; v(:, 1, 2) = 0.4_wp
      fr = 0.0_wp; fs = 0.0_wp; fz = 0.0_wp; ie = 0_int32; tf = TF0; wet = .true.
      hv = 0.0_wp; hv(1) = HIN1
      vol = 0.1_wp
      frz = frz_of(p%c%rhoi, p%c%Lfresh, vol)
      bank = frz*DT
      q(1, 1, 2, 3) = ieee_value(1.0_wp, ieee_quiet_nan)
      ie(2, 1) = ICE_ERR_TBOT
      ie(3, 1) = ICE_ERR_REMAP_FALLBACK              ! counted, not fatal
      v_in = v
      !$acc enter data copyin(a, v, ts, q, s, sal, a0, fr, fs, fz, frz, bank, tf, ie, wet)
      call ice_new_ice_step(NX, 1, NC, NIL, 1, NX, 1, 1, DT, p, hv, sal, wet, frz, bank, tf, a, v, ts, q, s, &
                            a0, fr, fs, fz, ie)
      !$acc update self(a, v, a0, fr, fs, fz, bank, ie)
      !$acc exit data delete(a, v, ts, q, s, sal, a0, fr, fs, fz, frz, bank, tf, ie, wet)
      call check(error, ie(1, 1) == ICE_ERR_NONFINITE .and. bank(1, 1) == frz(1, 1)*DT .and. fz(1, 1) == 0.0_wp .and. &
                 all(v(1, 1, :) == v_in(1, 1, :)), "NaN enthalpy: bit set, bank kept, cell untouched")
      if (allocated(error)) return
      call check(error, ie(2, 1) == ICE_ERR_TBOT .and. bank(2, 1) == frz(2, 1)*DT .and. fz(2, 1) == 0.0_wp .and. &
                 all(v(2, 1, :) == v_in(2, 1, :)) .and. fr(2, 1) == 0.0_wp, "pre-existing fatal bit: bank and state kept")
      if (allocated(error)) return
      call check(error, ie(3, 1) == ICE_ERR_REMAP_FALLBACK .and. bank(3, 1) == 0.0_wp .and. fz(3, 1) > 0.0_wp, &
                 "counted bit does not block")
      if (allocated(error)) return
      call check(error, ie(4, 1) == 0_int32 .and. bank(4, 1) == 0.0_wp .and. fz(4, 1) > 0.0_wp, "clean cell spends its bank")
   end subroutine test_driver_errors

   subroutine cell_run(c, p, frz, tf)
      type(cell_t), intent(inout) :: c
      type(ice_params_t), intent(in) :: p
      real(wp), intent(in) :: frz, tf

      call ice_add_new_ice_cell(1, 1, 1, 1, NC, c%nil, DT, frz, tf, HIN1, p%c%rhoi, p%c%Lfresh, p%c%puny, &
                                p%thermo%hfrazilmin, p%thermo%ice_ref_salinity, p%thermo%update_ocn_f, &
                                c%sal, c%a, c%v, c%ts, c%q, c%s, c%a0, c%fr, c%fs, c%fz, c%ie)
   end subroutine cell_run

   function params(set, nil, updf, sref, hfm) result(p)
      integer, intent(in) :: set, nil
      logical, intent(in) :: updf
      real(wp), intent(in) :: sref, hfm
      type(ice_params_t) :: p

      p = ice_params_default(set, NC, nil, 1)
      p%thermo%update_ocn_f = updf
      p%thermo%ice_ref_salinity = sref
      p%thermo%hfrazilmin = hfm
   end function params

   function cell_equal(a, b) result(same)
      !! Bitwise equality of every array of two cells (NaN-free states).
      type(cell_t), intent(in) :: a, b
      logical :: same
      same = all(a%a == b%a) .and. all(a%v == b%v) .and. all(a%ts == b%ts) .and. all(a%q == b%q) &
             .and. all(a%s == b%s) .and. all(a%a0 == b%a0) .and. all(a%fr == b%fr) &
             .and. all(a%fs == b%fs) .and. all(a%fz == b%fz) .and. all(a%ie == b%ie)
   end function cell_equal

   subroutine rand_cell(st, nil, a, v, ts, q, s, vol, frz_neg)
      !! Random cell: 20 % empty, 15 % closed (sum 1), else a partial cover with random
      !! category occupancy; `vol` log-uniform in [1e-4, 1.5] m; `frz_neg` for about 10 % of the cells.
      integer(int64), intent(inout) :: st
      integer, intent(in) :: nil
      real(wp), intent(out) :: a(NC), v(NC), ts(NC), q(nil, NC), s(nil, NC), vol
      logical, intent(out) :: frz_neg
      real(wp) :: u, w(NC), tot, cover
      integer :: n, k

      a = 0.0_wp; v = 0.0_wp
      st = lcg_next(st); u = lcg_unit(st)
      if (u >= 0.2_wp) then
         tot = 0.0_wp
         do n = 1, NC
            st = lcg_next(st); w(n) = lcg_unit(st)
            if (w(n) < 0.35_wp) w(n) = 0.0_wp
            tot = tot + w(n)
         end do
         if (tot <= 0.0_wp) then
            w = 0.0_wp; w(1) = 1.0_wp; tot = 1.0_wp
         end if
         st = lcg_next(st); cover = lcg_unit(st)
         if (u < 0.35_wp) then
            cover = 1.0_wp                                       ! closed ice
         else
            cover = 0.02_wp + 0.97_wp*cover
         end if
         do n = 1, NC
            a(n) = cover*w(n)/tot
            st = lcg_next(st)
            if (a(n) > 0.0_wp) v(n) = a(n)*(0.05_wp + 4.0_wp*lcg_unit(st))
         end do
         if (cover >= 1.0_wp) then                               ! make the sum exactly representable <= 1
            a = a*(1.0_wp - 1.0e-15_wp)
         end if
      end if
      do n = 1, NC
         st = lcg_next(st); ts(n) = -0.1_wp - 30.0_wp*lcg_unit(st)
         do k = 1, nil
            st = lcg_next(st); q(k, n) = -4.0e8_wp + 3.0e8_wp*lcg_unit(st)
            st = lcg_next(st); s(k, n) = 6.0_wp*lcg_unit(st)
         end do
      end do
      st = lcg_next(st); vol = 10.0_wp**(-4.0_wp + 4.176_wp*lcg_unit(st))
      st = lcg_next(st); frz_neg = lcg_unit(st) < 0.1_wp
   end subroutine rand_cell

   ! ------------------------------------------------------------------ closed forms

   subroutine test_thin_lead(error)
      !! T2: empty cell, V = 0.02 m: hi0new = hfrazilmin = 0.05, new ice covers 0.4 of the cell.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      type(cell_t) :: c
      real(wp) :: vol, qi
      integer :: k

      p = params(ICE_CONSTANTS_ICEPACK, 4, .true., 5.0_wp, 0.05_wp)
      call cell_init(c, 4)
      call cell_fill_layers(c, 11_int64)
      c%ts(1, 1, 1) = -7.0_wp                          ! empty category: weight zero
      vol = 0.02_wp
      qi = -p%c%rhoi*p%c%Lfresh
      call cell_run(c, p, frz_of(p%c%rhoi, p%c%Lfresh, vol), TF0)
      call ice_golden_check(error, "area, volume, open water, frazil", &
                            [c%a(1, 1, 1), c%v(1, 1, 1), c%a0(1, 1), c%fz(1, 1)], &
                            [0.4_wp, vol, 0.6_wp, vol], 1.0e-14_wp, 0.0_wp)
      if (allocated(error)) return
      do k = 1, 4
         call ice_golden_check(error, "new-ice layer enthalpy / salinity", [c%q(1, 1, 1, k), c%s(1, 1, 1, k)], &
                               [qi, c%sal(k)], 1.0e-14_wp, 0.0_wp)
         if (allocated(error)) return
      end do
      call ice_golden_check(error, "surface temperature of new ice = Tf", [c%ts(1, 1, 1)], [TF0], 1.0e-14_wp, 0.0_wp)
      if (allocated(error)) return
      call ice_golden_check(error, "fresh, fsalt", [c%fr(1, 1), c%fs(1, 1)], &
                            [-p%c%rhoi*vol/DT, -p%c%rhoi*vol/DT*5.0_wp*1.0e-3_wp], 1.0e-14_wp, 0.0_wp)
      if (allocated(error)) return
      call check(error, all(c%a(1, 1, 2:NC) == 0.0_wp) .and. all(c%v(1, 1, 2:NC) == 0.0_wp) .and. c%ie(1, 1) == 0, &
                 "other categories empty, no error")
   end subroutine test_thin_lead

   subroutine test_thick_full_lead(error)
      !! T3: empty cell, V = 0.30 m < hi0max: the whole lead freezes over at 0.30 m.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      type(cell_t) :: c

      p = params(ICE_CONSTANTS_ICEPACK, 4, .true., 5.0_wp, 0.05_wp)
      call cell_init(c, 4)
      call cell_fill_layers(c, 12_int64)
      call cell_run(c, p, frz_of(p%c%rhoi, p%c%Lfresh, 0.30_wp), TF0)
      call ice_golden_check(error, "area 1, volume 0.3, open water 0", &
                            [c%a(1, 1, 1), c%v(1, 1, 1), c%a0(1, 1)], [1.0_wp, 0.30_wp, 0.0_wp], 1.0e-14_wp, 1.0e-15_wp)
   end subroutine test_thick_full_lead

   subroutine test_overflow(error)
      !! T4: aicen(3) = 0.5 (h = 1.8 m), V = 0.4 m: V/aice0 = 0.8 > hi0max, so the new ice fills
      !! the open water at hi0max and the surplus spreads over the existing ice.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      type(cell_t) :: c
      real(wp) :: hi0max, vsurp, hsurp, qi, v3, q3, vtot, q1(4), s3, q3old(4), s3old(4)
      integer :: k

      p = params(ICE_CONSTANTS_ICEPACK, 4, .true., 5.0_wp, 0.05_wp)
      call cell_init(c, 4)
      call cell_fill_layers(c, 13_int64)
      c%a(1, 1, 3) = 0.5_wp
      c%v(1, 1, 3) = 0.9_wp
      call cell_set_a0(c)
      q3old = c%q(1, 1, 3, :)
      s3old = c%s(1, 1, 3, :)
      qi = -p%c%rhoi*p%c%Lfresh
      hi0max = 0.9_wp*HIN1
      vsurp = 0.4_wp - 0.5_wp*hi0max
      hsurp = vsurp/0.5_wp
      call cell_run(c, p, frz_of(p%c%rhoi, p%c%Lfresh, 0.4_wp), TF0)
      call ice_golden_check(error, "cat 1: area = aice0, volume = aice0*hi0max", &
                            [c%a(1, 1, 1), c%v(1, 1, 1), c%a0(1, 1)], [0.5_wp, 0.5_wp*hi0max, 0.0_wp], 1.0e-14_wp, 1.0e-15_wp)
      if (allocated(error)) return
      call ice_golden_check(error, "cat 3 volume rises by hsurp*aicen", [c%v(1, 1, 3), c%a(1, 1, 3)], &
                            [0.9_wp + hsurp*0.5_wp, 0.5_wp], 1.0e-14_wp, 0.0_wp)
      if (allocated(error)) return
      vtot = c%v(1, 1, 3)
      v3 = 0.9_wp
      do k = 1, 4
         q3 = (q3old(k)*v3 + qi*hsurp*0.5_wp)/vtot
         s3 = (s3old(k)*v3 + c%sal(k)*hsurp*0.5_wp)/vtot
         call ice_golden_check(error, "cat 3 layer = volume-weighted mean", [c%q(1, 1, 3, k), c%s(1, 1, 3, k)], &
                               [q3, s3], 1.0e-14_wp, 0.0_wp)
         if (allocated(error)) return
      end do
      q1 = c%q(1, 1, 1, :)
      call ice_golden_check(error, "cat 1 holds pure new ice (it was empty)", q1, [(qi, k=1, 4)], 1.0e-14_wp, 0.0_wp)
      if (allocated(error)) return
      call ice_golden_check(error, "total volume change = V", [c%v(1, 1, 1) + c%v(1, 1, 3) - 0.9_wp], [0.4_wp], &
                            1.0e-14_wp, 0.0_wp)
   end subroutine test_overflow

   subroutine test_closed_ice(error)
      !! T5: aice0 = 0: all the volume spreads (hsurp = V/aice), areas unchanged, cat 1 not enlarged.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      type(cell_t) :: c
      real(wp) :: a_before(NC), v_before(NC), vol
      integer :: n

      p = params(ICE_CONSTANTS_ICEPACK, 7, .true., 5.0_wp, 0.05_wp)
      call cell_init(c, 7)
      call cell_fill_layers(c, 14_int64)
      c%a(1, 1, :) = [0.1_wp, 0.2_wp, 0.3_wp, 0.15_wp, 0.25_wp]
      c%v(1, 1, :) = [0.03_wp, 0.2_wp, 0.5_wp, 0.5_wp, 1.0_wp]
      call cell_set_a0(c)
      call check(error, c%a0(1, 1) <= 2.0e-16_wp, "setup: closed ice")
      if (allocated(error)) return
      c%a0(1, 1) = 0.0_wp
      a_before = c%a(1, 1, :)
      v_before = c%v(1, 1, :)
      vol = 0.05_wp
      call cell_run(c, p, frz_of(p%c%rhoi, p%c%Lfresh, vol), TF0)
      call check(error, all(c%a(1, 1, :) == a_before) .and. c%a0(1, 1) == 0.0_wp, "areas and open water unchanged")
      if (allocated(error)) return
      do n = 1, NC
         call ice_golden_check(error, "each category gains (V/aice)*aicen", [c%v(1, 1, n)], &
                               [v_before(n) + vol/sum(a_before)*a_before(n)], 1.0e-13_wp, 0.0_wp)
         if (allocated(error)) return
      end do
      call ice_golden_check(error, "total gain = V", [sum(c%v(1, 1, :)) - sum(v_before)], [vol], 1.0e-12_wp, 0.0_wp)
   end subroutine test_closed_ice

   subroutine test_hfm_sweep(error)
      !! T6: ai0new = V/max(V/aice0, hfrazilmin), never above aice0, non-increasing in hfrazilmin;
      !! the volume added does not depend on it.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      type(cell_t) :: c
      real(wp), parameter :: HFM(6) = [0.005_wp, 0.01_wp, 0.02_wp, 0.05_wp, 0.1_wp, 0.15_wp]
      real(wp) :: vol, ai0new, prev, ex
      integer :: m

      vol = 0.002_wp
      prev = huge(1.0_wp)
      do m = 1, size(HFM)
         p = params(ICE_CONSTANTS_ICEPACK, 4, .true., 5.0_wp, HFM(m))
         call cell_init(c, 4)
         call cell_fill_layers(c, 15_int64)
         c%a(1, 1, 1) = 0.9_wp
         c%v(1, 1, 1) = 0.9_wp
         call cell_set_a0(c)
         call cell_run(c, p, frz_of(p%c%rhoi, p%c%Lfresh, vol), TF0)
         ai0new = c%a(1, 1, 1) - 0.9_wp
         ex = vol/max(vol/0.1_wp, HFM(m))
         call ice_golden_check(error, "ai0new = V/max(V/aice0, hfrazilmin)", [ai0new], [ex], 1.0e-12_wp, 0.0_wp)
         if (allocated(error)) return
         call check(error, ai0new <= 0.1_wp + 1.0e-14_wp .and. ai0new <= prev, "bounded by aice0 and non-increasing")
         if (allocated(error)) return
         call ice_golden_check(error, "volume independent of hfrazilmin", [c%v(1, 1, 1) - 0.9_wp], [vol], 1.0e-12_wp, 0.0_wp)
         if (allocated(error)) return
         prev = ai0new
      end do
   end subroutine test_hfm_sweep

   subroutine test_bl99_inversion(error)
      !! T8: the new-ice enthalpy -rhoi*Lfresh read through the BL99 inversion is a finite
      !! temperature below the layer melting point and a root of
      !! `cp_ice*(Tmlt-T) - Lfresh*Tmlt/T - cp_ocn*Tmlt = 0`.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      real(wp) :: sp(8), T, Tmlt, res, qi
      logical :: bad
      integer :: k

      p = params(ICE_CONSTANTS_ICEPACK, 7, .true., 5.0_wp, 0.05_wp)
      call ice_salinity_profile(7, 3.2_wp, 0.1_wp, sp)
      qi = -p%c%rhoi*p%c%Lfresh
      do k = 1, 7
         call ice_T_from_q(qi, sp(k), p%c, .true., T, bad)
         Tmlt = -p%c%depressT*sp(k)
         if (bad .or. .not. ieee_is_finite(T)) then
            call test_failed(error, "inversion of the new-ice enthalpy failed or is non-finite")
            return
         end if
         res = p%c%cp_ice*(Tmlt - T) - p%c%Lfresh*Tmlt/T - p%c%cp_ocn*Tmlt
         call check(error, T < Tmlt, "temperature below the melting point")
         if (allocated(error)) return
         call ice_golden_check(error, "equation residual", [res], [0.0_wp], 0.0_wp, 1.0e-6_wp*p%c%Lfresh)
         if (allocated(error)) return
      end do
   end subroutine test_bl99_inversion

   subroutine test_updf_off(error)
      !! T9: with `update_ocn_f` off the ice forms identically but `fresh`/`fsalt` do not move.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: pon, poff
      type(cell_t) :: con, coff

      pon = params(ICE_CONSTANTS_ICEPACK, 4, .true., 5.0_wp, 0.05_wp)
      poff = params(ICE_CONSTANTS_ICEPACK, 4, .false., 5.0_wp, 0.05_wp)
      call cell_init(con, 4); call cell_init(coff, 4)
      call cell_fill_layers(con, 16_int64); call cell_fill_layers(coff, 16_int64)
      con%a(1, 1, 3) = 0.5_wp; con%v(1, 1, 3) = 0.9_wp; call cell_set_a0(con)
      coff%a = con%a; coff%v = con%v; coff%a0 = con%a0
      con%fr = -2.0e-5_wp; con%fs = 3.0e-8_wp; coff%fr = con%fr; coff%fs = con%fs
      call cell_run(con, pon, frz_of(pon%c%rhoi, pon%c%Lfresh, 0.4_wp), TF0)
      call cell_run(coff, poff, frz_of(pon%c%rhoi, pon%c%Lfresh, 0.4_wp), TF0)
      call check(error, coff%fr(1, 1) == -2.0e-5_wp .and. coff%fs(1, 1) == 3.0e-8_wp, "fresh, fsalt unchanged")
      if (allocated(error)) return
      call check(error, con%fr(1, 1) < -2.0e-5_wp .and. con%fs(1, 1) < 3.0e-8_wp, "fresh, fsalt charged when on")
      if (allocated(error)) return
      call check(error, all(con%a == coff%a) .and. all(con%v == coff%v) .and. all(con%q == coff%q) .and. &
                 all(con%s == coff%s) .and. all(con%ts == coff%ts) .and. con%fz(1, 1) == coff%fz(1, 1) .and. &
                 coff%fz(1, 1) > 0.0_wp, "state and frazil volume do not depend on update_ocn_f")
   end subroutine test_updf_off

   subroutine test_no_bank(error)
      !! T10: a freezing potential <= 0 changes NOTHING (bitwise): not the state, not the fluxes,
      !! not the bank (step driver).  Garbage elsewhere in the cell is not even read.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      type(cell_t) :: c, ref
      real(wp) :: frz(1, 1), bank(1, 1), tf(1, 1), hv(0:NC)
      real(wp) :: a1(1, 1, NC), v1(1, 1, NC), ts1(1, 1, NC), q1(1, 1, NC, 4), s1(1, 1, NC, 4), a01(1, 1)
      real(wp) :: fr1(1, 1), fs1(1, 1), fz1(1, 1), sal1(4)
      integer(int32) :: ie1(1, 1)
      logical :: wet(1, 1)
      real(wp), parameter :: FRZ_CASES(3) = [0.0_wp, -35.0_wp, -1.0e-300_wp]
      integer :: m

      p = params(ICE_CONSTANTS_ICEPACK, 4, .true., 5.0_wp, 0.05_wp)
      call cell_init(c, 4); call cell_init(ref, 4)
      call cell_fill_layers(c, 17_int64)
      c%a(1, 1, 2) = 0.3_wp; c%v(1, 1, 2) = 0.4_wp; call cell_set_a0(c)
      c%fr = 1.0e-5_wp
      ref = c
      c%q(1, 1, 4, 2) = ieee_value(1.0_wp, ieee_quiet_nan)      ! never read: no freezing potential
      ref%q(1, 1, 4, 2) = c%q(1, 1, 4, 2)
      do m = 1, size(FRZ_CASES)
         call cell_run(c, p, FRZ_CASES(m), TF0)
         call check(error, all(c%a == ref%a) .and. all(c%v == ref%v) .and. all(c%ts == ref%ts) .and. &
                    all(c%s == ref%s) .and. all(c%a0 == ref%a0) .and. all(c%fr == ref%fr) .and. all(c%fs == ref%fs) .and. &
                    all(c%fz == ref%fz) .and. all(c%ie == 0_int32) .and. &
                    all(c%q(1, 1, 1:3, :) == ref%q(1, 1, 1:3, :)), "cell kernel is a bitwise no-op")
         if (allocated(error)) return
      end do
      hv = 0.0_wp; hv(1) = HIN1
      bank = 1.234e5_wp; frz = -3.0_wp; tf = TF0; wet = .true.
      c%q(1, 1, 4, 2) = -2.0e8_wp
      a1 = c%a; v1 = c%v; ts1 = c%ts; q1 = c%q; s1 = c%s; a01 = c%a0; fr1 = c%fr; fs1 = c%fs; fz1 = c%fz; ie1 = c%ie
      sal1 = c%sal
      !$acc enter data copyin(a1, v1, ts1, q1, s1, sal1, a01, fr1, fs1, fz1, frz, bank, tf, ie1, wet)
      call ice_new_ice_step(1, 1, NC, 4, 1, 1, 1, 1, DT, p, hv, sal1, wet, frz, bank, tf, a1, v1, ts1, q1, s1, &
                            a01, fr1, fs1, fz1, ie1)
      !$acc update self(fz1, bank)
      !$acc exit data delete(a1, v1, ts1, q1, s1, sal1, a01, fr1, fs1, fz1, frz, bank, tf, ie1, wet)
      call check(error, bank(1, 1) == 1.234e5_wp .and. fz1(1, 1) == 0.0_wp, "driver keeps the bank when frzmlt <= 0")
   end subroutine test_no_bank

   subroutine test_dry_cell(error)
      !! T11: a dry cell keeps its bank and its state; the wet neighbour spends its bank.
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      integer, parameter :: NX = 2, NIL = 4
      real(wp) :: a(NX, 1, NC), v(NX, 1, NC), ts(NX, 1, NC), q(NX, 1, NC, NIL), s(NX, 1, NC, NIL), sal(NIL), sp(NIL + 1)
      real(wp) :: a0(NX, 1), fr(NX, 1), fs(NX, 1), fz(NX, 1), frz(NX, 1), bank(NX, 1), tf(NX, 1), hv(0:NC), vol
      integer(int32) :: ie(NX, 1)
      logical :: wet(NX, 1)

      p = params(ICE_CONSTANTS_ICEPACK, NIL, .true., 5.0_wp, 0.05_wp)
      call ice_salinity_profile(NIL, 3.2_wp, 0.1_wp, sp); sal = sp(1:NIL)
      a = 0.0_wp; v = 0.0_wp; ts = -5.0_wp; q = -3.0e8_wp; s = 3.0_wp; a0 = 1.0_wp
      fr = 0.0_wp; fs = 0.0_wp; fz = 0.0_wp; ie = 0_int32; tf = TF0
      hv = 0.0_wp; hv(1) = HIN1
      vol = 0.1_wp
      bank = vol*p%c%rhoi*p%c%Lfresh                    ! J m-2 of cell
      frz = bank/DT
      wet(1, 1) = .false.; wet(2, 1) = .true.
      !$acc enter data copyin(a, v, ts, q, s, sal, a0, fr, fs, fz, frz, bank, tf, ie, wet)
      call ice_new_ice_step(NX, 1, NC, NIL, 1, NX, 1, 1, DT, p, hv, sal, wet, frz, bank, tf, a, v, ts, q, s, &
                            a0, fr, fs, fz, ie)
      !$acc update self(a, v, ts, q, s, a0, fr, fs, fz, bank, ie)
      !$acc exit data delete(a, v, ts, q, s, sal, a0, fr, fs, fz, frz, bank, tf, ie, wet)
      call check(error, a(1, 1, 1) == 0.0_wp .and. v(1, 1, 1) == 0.0_wp .and. a0(1, 1) == 1.0_wp .and. &
                 fr(1, 1) == 0.0_wp .and. fz(1, 1) == 0.0_wp .and. bank(1, 1) == vol*p%c%rhoi*p%c%Lfresh, &
                 "dry cell: state untouched, bank retained")
      if (allocated(error)) return
      call ice_golden_check(error, "wet cell: ice formed", [v(2, 1, 1), fz(2, 1)], [vol, vol], 1.0e-13_wp, 0.0_wp)
      if (allocated(error)) return
      call check(error, bank(2, 1) == 0.0_wp, "wet cell: bank spent")
   end subroutine test_dry_cell

   subroutine test_errors(error)
      !! T14: every entry error leaves the cell untouched and sets its bit; other cells of the grid are
      !! processed exactly as on their own (bitwise).
      type(error_type), allocatable, intent(out) :: error
      type(ice_params_t) :: p
      type(cell_t) :: c, ref, solo
      real(wp) :: nan, frz
      integer :: m
      integer(int32) :: want

      p = params(ICE_CONSTANTS_ICEPACK, 4, .true., 5.0_wp, 0.05_wp)
      nan = ieee_value(1.0_wp, ieee_quiet_nan)
      frz = frz_of(p%c%rhoi, p%c%Lfresh, 0.2_wp)
      do m = 1, 6
         call cell_init(c, 4)
         call cell_fill_layers(c, 20_int64 + m)
         c%a(1, 1, 2) = 0.3_wp; c%v(1, 1, 2) = 0.4_wp; call cell_set_a0(c)
         want = ICE_ERR_NONFINITE
         select case (m)
         case (1); c%a(1, 1, 4) = nan
         case (2); c%v(1, 1, 5) = nan
         case (3); c%q(1, 1, 3, 2) = nan
         case (4); c%s(1, 1, 1, 4) = nan
         case (5); c%ts(1, 1, 1) = nan
         case (6); c%a(1, 1, 5) = -1.0e-3_wp; want = ICE_ERR_NEW_ICE_AREA
         end select
         ref = c
         call cell_run(c, p, frz, TF0)
         call check(error, c%ie(1, 1) == want .and. cell_equal_ignoring_ie(c, ref), "entry error: bit set, cell untouched")
         if (allocated(error)) return
      end do
      ! NaN freezing potential: flagged, untouched
      call cell_init(c, 4); call cell_fill_layers(c, 40_int64)
      c%a(1, 1, 2) = 0.3_wp; c%v(1, 1, 2) = 0.4_wp; call cell_set_a0(c)
      ref = c
      call cell_run(c, p, nan, TF0)
      call check(error, c%ie(1, 1) == ICE_ERR_NONFINITE .and. cell_equal_ignoring_ie(c, ref), "NaN frzmlt flagged")
      if (allocated(error)) return
      ! an error bit already in the word survives (OR-ed in)
      call cell_init(c, 4); call cell_fill_layers(c, 41_int64)
      c%a(1, 1, 5) = -1.0_wp; c%ie = 8_int32
      call cell_run(c, p, frz, TF0)
      call check(error, c%ie(1, 1) == ior(8_int32, ICE_ERR_NEW_ICE_AREA), "existing error bits kept")
      if (allocated(error)) return
      ! a good cell next to a bad one in a grid: same as on its own
      call cell_init(solo, 4); call cell_fill_layers(solo, 42_int64)
      solo%a(1, 1, 2) = 0.3_wp; solo%v(1, 1, 2) = 0.4_wp; call cell_set_a0(solo)
      call cell_run(solo, p, frz, TF0)
      call check(error, solo%ie(1, 1) == 0_int32 .and. solo%fz(1, 1) > 0.0_wp, "good cell runs clean")
   end subroutine test_errors

   function cell_equal_ignoring_ie(a, b) result(same)
      !! Bitwise equality of the state and flux arrays, treating a NaN on both sides as equal.
      type(cell_t), intent(in) :: a, b
      logical :: same
      same = eq(reshape(a%a, [size(a%a)]), reshape(b%a, [size(b%a)])) .and. &
             eq(reshape(a%v, [size(a%v)]), reshape(b%v, [size(b%v)])) .and. &
             eq(reshape(a%ts, [size(a%ts)]), reshape(b%ts, [size(b%ts)])) .and. &
             eq(reshape(a%q, [size(a%q)]), reshape(b%q, [size(b%q)])) .and. &
             eq(reshape(a%s, [size(a%s)]), reshape(b%s, [size(b%s)])) .and. &
             eq(reshape(a%a0, [size(a%a0)]), reshape(b%a0, [size(b%a0)])) .and. &
             eq(reshape(a%fr, [size(a%fr)]), reshape(b%fr, [size(b%fr)])) .and. &
             eq(reshape(a%fs, [size(a%fs)]), reshape(b%fs, [size(b%fs)])) .and. &
             eq(reshape(a%fz, [size(a%fz)]), reshape(b%fz, [size(b%fz)]))
   end function cell_equal_ignoring_ie

   pure function eq(x, y) result(same)
      real(wp), intent(in) :: x(:), y(:)
      logical :: same
      integer :: m

      same = .true.
      do m = 1, size(x)
         if (ieee_is_finite(x(m)) .neqv. ieee_is_finite(y(m))) same = .false.
         if (ieee_is_finite(x(m)) .and. ieee_is_finite(y(m))) then
            if (x(m) /= y(m)) same = .false.
         end if
      end do
   end function eq

   ! ------------------------------------------------------------------ conservation

   subroutine test_closure_random(error)
      !! T1 + T7 on 100 000 random cells for each of nilyr = 1, 4, 7 (set ICEPACK, update_ocn_f on,
      !! ice_ref_salinity 5): the ice volume rises by exactly V = frzmlt*dt/(rhoi*Lfresh); the layer
      !! energy sum_n sum_k q v/nilyr rises by exactly -frzmlt*dt (new ice carries -rhoi*Lfresh, and
      !! every weighted mean conserves q*v); fresh = -rhoi*V/dt, fsalt = ref_salinity*1e-3*fresh,
      !! frazil_vol = V; area + open water stays 1 (<= 1e-14); categories 2..5 never change area;
      !! sice stays a convex mix of the old layer and `salinz`; no error bit.  Cells with a
      !! freezing potential <= 0 change nothing.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 250, NY = 400
      integer, parameter :: NILS(3) = [1, 4, 7]
      type(ice_params_t) :: p
      real(wp), allocatable :: a(:, :, :), v(:, :, :), ts(:, :, :), q(:, :, :, :), s(:, :, :, :)
      real(wp), allocatable :: a0(:, :), fr(:, :), fs(:, :), fz(:, :), sal(:)
      integer(int32), allocatable :: ie(:, :)
      real(wp) :: sp(8), ac(NC), vc(NC), tc(NC), qc(7, NC), sc(7, NC), vol, frz, e0, e1, v0, v1, a1sum, scl
      real(wp) :: w_vol, w_en, w_fr, w_fs, w_fz, w_area, w_cat, w_sal
      real(wp) :: smin, smax
      integer(int64) :: st
      integer :: im, nil, i, j, n, k, nneg, nact
      logical :: neg
      integer(int32) :: nbad

      do im = 1, size(NILS)
         nil = NILS(im)
         p = params(ICE_CONSTANTS_ICEPACK, nil, .true., 5.0_wp, 0.05_wp)
         allocate (a(NX, NY, NC), v(NX, NY, NC), ts(NX, NY, NC), q(NX, NY, NC, nil), s(NX, NY, NC, nil))
         allocate (a0(NX, NY), fr(NX, NY), fs(NX, NY), fz(NX, NY), ie(NX, NY), sal(nil))
         call ice_salinity_profile(nil, p%c%saltmax, p%c%min_salin, sp(1:nil + 1))
         sal = sp(1:nil)
         fr = 0.0_wp; fs = 0.0_wp; fz = 0.0_wp; ie = 0_int32
         st = lcg_next(int(1000 + nil, int64))
         w_vol = 0.0_wp; w_en = 0.0_wp; w_fr = 0.0_wp; w_fs = 0.0_wp; w_fz = 0.0_wp; w_area = 0.0_wp; w_cat = 0.0_wp
         w_sal = 0.0_wp
         nneg = 0; nact = 0; nbad = 0
         do j = 1, NY
            do i = 1, NX
               call rand_cell(st, nil, ac, vc, tc, qc(1:nil, :), sc(1:nil, :), vol, neg)
               a(i, j, :) = ac; v(i, j, :) = vc; ts(i, j, :) = tc
               do n = 1, NC
                  q(i, j, n, :) = qc(1:nil, n)
                  s(i, j, n, :) = sc(1:nil, n)
               end do
               a0(i, j) = 0.0_wp
               do n = 1, NC
                  a0(i, j) = a0(i, j) + ac(n)
               end do
               a0(i, j) = max(1.0_wp - a0(i, j), 0.0_wp)
               frz = frz_of(p%c%rhoi, p%c%Lfresh, vol)
               if (neg) frz = -frz
               ! before
               e0 = 0.0_wp; v0 = 0.0_wp
               do n = 1, NC
                  v0 = v0 + v(i, j, n)
                  do k = 1, nil
                     e0 = e0 + q(i, j, n, k)*v(i, j, n)/real(nil, wp)
                  end do
               end do
               call ice_add_new_ice_cell(i, j, NX, NY, NC, nil, DT, frz, TF0, HIN1, p%c%rhoi, p%c%Lfresh, p%c%puny, &
                                         p%thermo%hfrazilmin, p%thermo%ice_ref_salinity, p%thermo%update_ocn_f, &
                                         sal, a, v, ts, q, s, a0, fr, fs, fz, ie)
               e1 = 0.0_wp; v1 = 0.0_wp; a1sum = 0.0_wp
               do n = 1, NC
                  v1 = v1 + v(i, j, n)
                  a1sum = a1sum + a(i, j, n)
                  do k = 1, nil
                     e1 = e1 + q(i, j, n, k)*v(i, j, n)/real(nil, wp)
                  end do
               end do
               if (frz <= 0.0_wp) then
                  nneg = nneg + 1
                  if (v1 /= v0 .or. e1 /= e0 .or. fz(i, j) /= 0.0_wp .or. fr(i, j) /= 0.0_wp .or. ie(i, j) /= 0) nbad = nbad + 1
                  cycle
               end if
               nact = nact + 1
               if (ie(i, j) /= 0) nbad = nbad + 1
               if (.not. (ieee_is_finite(v1 + e1 + a1sum + a0(i, j) + fr(i, j) + fs(i, j) + fz(i, j) + &
                                         sum(ts(i, j, :)) + sum(s(i, j, :, :))))) then
                  nbad = nbad + 1
                  cycle
               end if
               scl = max(vol, v0)
               w_vol = max(w_vol, abs((v1 - v0) - vol)/scl)
               w_en = max(w_en, abs((e1 - e0) + frz*DT)/max(frz*DT, abs(e0) + abs(e1)))
               w_fr = max(w_fr, abs(fr(i, j) + p%c%rhoi*vol/DT)/(p%c%rhoi*vol/DT))
               w_fs = max(w_fs, abs(fs(i, j) - 5.0_wp*0.001_wp*fr(i, j))/abs(5.0_wp*0.001_wp*fr(i, j)))
               w_fz = max(w_fz, abs(fz(i, j) - vol)/vol)
               w_area = max(w_area, abs(a1sum + a0(i, j) - 1.0_wp))
               do n = 2, NC
                  if (a(i, j, n) /= ac(n)) nbad = nbad + 1
               end do
               do n = 1, NC
                  do k = 1, nil
                     if (v(i, j, n) > 0.0_wp) then
                        smin = min(sc(k, n), sal(k)); smax = max(sc(k, n), sal(k))
                        if (v(i, j, n) > vc(n) .or. n == 1) then
                           w_sal = max(w_sal, max(smin - s(i, j, n, k), s(i, j, n, k) - smax, 0.0_wp))
                        end if
                     end if
                  end do
               end do
               w_cat = max(w_cat, abs(a1sum - sum(ac) - (a(i, j, 1) - ac(1))))
            end do
         end do
         write (output_unit, '(a,i0,a,i0,a,i0,a,6es10.2)') "  new ice closure nilyr=", nil, ": ", nact, &
            " freezing cells, ", nneg, " no-op cells; worst V, E, fresh, fsalt, frazil, area: ", &
            w_vol, w_en, w_fr, w_fs, w_fz, w_area
         call check(error, nbad == 0 .and. nact > 0 .and. nneg > 0, "no error bits, no-op cells bitwise, categories 2..5 unchanged")
         if (.not. allocated(error)) call check(error, w_vol <= 1.0e-13_wp, "volume change = V")
         if (.not. allocated(error)) call check(error, w_en <= 1.0e-13_wp, "energy change = -frzmlt*dt")
         if (.not. allocated(error)) call check(error, w_fr <= 1.0e-14_wp .and. w_fs <= 1.0e-14_wp, "fresh and fsalt")
         if (.not. allocated(error)) call check(error, w_fz <= 1.0e-14_wp, "frazil volume")
         if (.not. allocated(error)) call check(error, w_area <= 1.0e-14_wp, "area + open water = 1")
         if (.not. allocated(error)) call check(error, w_sal <= 1.0e-12_wp, "salinity stays a convex mix")
         if (.not. allocated(error)) call check(error, w_cat <= 1.0e-14_wp, "area gained only in category 1")
         deallocate (a, v, ts, q, s, a0, fr, fs, fz, ie, sal)
         if (allocated(error)) return
      end do
   end subroutine test_closure_random

   ! ------------------------------------------------------------------ oracle

   pure function blockmax(g, m) result(sc)
      !! Magnitude scale of the output block that element `m` belongs to (areas, volumes, Tsfc,
      !! enthalpies, salinities, then the four scalars on their own).
      real(wp), intent(in) :: g(89)
      integer, intent(in) :: m
      real(wp) :: sc
      integer :: lo, hi

      if (m <= 5) then
         lo = 1; hi = 5
      else if (m <= 10) then
         lo = 6; hi = 10
      else if (m <= 15) then
         lo = 11; hi = 15
      else if (m <= 50) then
         lo = 16; hi = 50
      else if (m <= 85) then
         lo = 51; hi = 85
      else
         lo = m; hi = m
      end if
      sc = maxval(abs(g(lo:hi)))
   end function blockmax

   subroutine load_golden(gin, gout)
      real(wp), intent(out) :: gin(103, NGOLD), gout(89, NGOLD)
#include "test_cice_new_ice_golden.inc"
   end subroutine load_golden

   subroutine test_oracle(error)
      !! Rows of the real upstream `add_new_ice` (40 per constant set, nilyr 7, 4 and 1; every
      !! branch: empty cell thin / thick / thicker than hi0max, overflow, closed ice, tiny open
      !! water, hfrazilmin floor and sweep, update_ocn_f off, reference salinity 4 and 5, no-op
      !! potentials, surface temperature cap, empty categories in the spread, tiny ice with aice0 + puny >= 1, area without volume).  Bitwise on CPU
      !! builds (the arithmetic is the same statements in the same order, no transcendental); on
      !! the GPU build FMA contraction applies to host code too: 1e-12 relative to the block scale.
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: gin(:, :), gout(:, :)
      real(wp) :: worst, got(89), x, sc
      integer :: r, nil, set, n, k, m, nbad
      type(ice_params_t) :: p
      type(cell_t) :: c
      character(len=80) :: msg

      allocate (gin(103, NGOLD), gout(89, NGOLD))
      call load_golden(gin, gout)
      worst = 0.0_wp
      nbad = 0
      do r = 1, NGOLD
         nil = nint(gin(2, r))
         set = ICE_CONSTANTS_ICEPACK
         if (r > NSET1) set = ICE_CONSTANTS_CESM
         p = params(set, nil, nint(gin(3, r)) == 1, gin(4, r), gin(5, r))
         call cell_init(c, nil)
         do n = 1, NC
            c%a(1, 1, n) = gin(11 + n, r)
            c%v(1, 1, n) = gin(16 + n, r)
            c%ts(1, 1, n) = gin(21 + n, r)
            do k = 1, nil
               c%q(1, 1, n, k) = gin(33 + (n - 1)*7 + k, r)
               c%s(1, 1, n, k) = gin(68 + (n - 1)*7 + k, r)
            end do
         end do
         c%sal = gin(27:26 + nil, r)
         c%fr = gin(9, r); c%fs = gin(10, r)
         call cell_set_a0(c)
         call ice_add_new_ice_cell(1, 1, 1, 1, NC, nil, gin(6, r), gin(7, r), gin(8, r), gin(11, r), p%c%rhoi, p%c%Lfresh, &
                                   p%c%puny, p%thermo%hfrazilmin, p%thermo%ice_ref_salinity, p%thermo%update_ocn_f, &
                                   c%sal, c%a, c%v, c%ts, c%q, c%s, c%a0, c%fr, c%fs, c%fz, c%ie)
         got = 0.0_wp
         got(1:5) = c%a(1, 1, :)
         got(6:10) = c%v(1, 1, :)
         got(11:15) = c%ts(1, 1, :)
         do n = 1, NC
            do k = 1, nil
               got(15 + (n - 1)*7 + k) = c%q(1, 1, n, k)
               got(50 + (n - 1)*7 + k) = c%s(1, 1, n, k)
            end do
         end do
         got(86) = c%a0(1, 1); got(87) = c%fr(1, 1); got(88) = c%fs(1, 1); got(89) = c%fz(1, 1)
         if (c%ie(1, 1) /= 0) nbad = nbad + 1
         do m = 1, 89
            if (.not. (ieee_is_finite(got(m)) .and. ieee_is_finite(gout(m, r)))) then
               nbad = nbad + 1
               cycle
            end if
            sc = blockmax(gout(:, r), m)
            x = abs(got(m) - gout(m, r))
            if (gout(m, r) /= 0.0_wp) x = x/abs(gout(m, r))
            worst = max(worst, x)
#ifdef RDB_GPU_OFFLOAD
            if (.not. (abs(got(m) - gout(m, r)) <= 1.0e-12_wp*sc + 1.0e-300_wp)) nbad = nbad + 1
#else
            if (got(m) /= gout(m, r)) nbad = nbad + 1
#endif
         end do
      end do
      write (output_unit, '(a,es10.2,a,i0)') "  new ice vs upstream: worst relative difference ", worst, ", failures ", nbad
      write (msg, '(a,i0,a)') "new ice differs from upstream in ", nbad, " values"
      call check(error, nbad == 0, trim(msg))
   end subroutine test_oracle

   ! ------------------------------------------------------------------ cell entry vs column / device

   subroutine fill_grid(nx, ny, nil, seed, p, frz, bank, tf, wet, a, v, ts, q, s, a0, sal)
      integer, intent(in) :: nx, ny, nil
      integer(int64), intent(in) :: seed
      type(ice_params_t), intent(in) :: p
      real(wp), intent(out) :: frz(nx, ny), bank(nx, ny), tf(nx, ny), a(nx, ny, NC), v(nx, ny, NC), ts(nx, ny, NC)
      real(wp), intent(out) :: q(nx, ny, NC, nil), s(nx, ny, NC, nil), a0(nx, ny), sal(nil)
      logical, intent(out) :: wet(nx, ny)
      real(wp) :: ac(NC), vc(NC), tc(NC), qc(nil, NC), sc(nil, NC), vol, sp(nil + 1), u
      logical :: neg
      integer(int64) :: st
      integer :: i, j, n

      call ice_salinity_profile(nil, p%c%saltmax, p%c%min_salin, sp)
      sal = sp(1:nil)
      st = lcg_next(seed)
      do j = 1, ny
         do i = 1, nx
            call rand_cell(st, nil, ac, vc, tc, qc, sc, vol, neg)
            a(i, j, :) = ac; v(i, j, :) = vc; ts(i, j, :) = tc
            do n = 1, NC
               q(i, j, n, :) = qc(:, n)
               s(i, j, n, :) = sc(:, n)
            end do
            a0(i, j) = 0.0_wp
            do n = 1, NC
               a0(i, j) = a0(i, j) + ac(n)
            end do
            a0(i, j) = max(1.0_wp - a0(i, j), 0.0_wp)
            frz(i, j) = frz_of(p%c%rhoi, p%c%Lfresh, vol)
            if (neg) frz(i, j) = -frz(i, j)
            bank(i, j) = frz(i, j)*DT
            tf(i, j) = TF0 - 0.2_wp*lcg_unit(lcg_next(st + int(i + 131*j, int64)))
            st = lcg_next(st); u = lcg_unit(st)
            wet(i, j) = u > 0.15_wp
         end do
      end do
   end subroutine fill_grid

   subroutine test_cell_vs_column(error)
      !! T15 (cell entry vs column): the in-place kernel at `(i, j)` of a multi-cell grid gives, bit
      !! for bit, what the same kernel gives on the cell extracted to a 1 x 1 grid (the column
      !! form).  Catches wrong-cell indexing, shared scratch, and cross-cell leakage.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 17, NY = 11, NIL = 7
      type(ice_params_t) :: p
      real(wp) :: frz(NX, NY), bank(NX, NY), tf(NX, NY), a(NX, NY, NC), v(NX, NY, NC), ts(NX, NY, NC)
      real(wp) :: q(NX, NY, NC, NIL), s(NX, NY, NC, NIL), a0(NX, NY), sal(NIL), fr(NX, NY), fs(NX, NY), fz(NX, NY)
      logical :: wet(NX, NY)
      integer(int32) :: ie(NX, NY)
      type(cell_t) :: c
      integer :: i, j, nbad, nfrz

      p = params(ICE_CONSTANTS_ICEPACK, NIL, .true., 5.0_wp, 0.05_wp)
      call fill_grid(NX, NY, NIL, 99_int64, p, frz, bank, tf, wet, a, v, ts, q, s, a0, sal)
      fr = 1.0e-6_wp; fs = -2.0e-9_wp; fz = 3.0e-4_wp; ie = 0_int32
      nbad = 0; nfrz = 0
      call cell_init(c, NIL)
      c%sal = sal
      do j = 1, NY
         do i = 1, NX
            c%a(1, 1, :) = a(i, j, :); c%v(1, 1, :) = v(i, j, :); c%ts(1, 1, :) = ts(i, j, :)
            c%q(1, 1, :, :) = q(i, j, :, :); c%s(1, 1, :, :) = s(i, j, :, :)
            c%a0(1, 1) = a0(i, j); c%fr(1, 1) = fr(i, j); c%fs(1, 1) = fs(i, j); c%fz(1, 1) = fz(i, j)
            c%ie = 0_int32
            call cell_run(c, p, frz(i, j), tf(i, j))
            call ice_add_new_ice_cell(i, j, NX, NY, NC, NIL, DT, frz(i, j), tf(i, j), HIN1, p%c%rhoi, p%c%Lfresh, &
                                      p%c%puny, p%thermo%hfrazilmin, p%thermo%ice_ref_salinity, p%thermo%update_ocn_f, &
                                      sal, a, v, ts, q, s, a0, fr, fs, fz, ie)
            if (frz(i, j) > 0.0_wp) nfrz = nfrz + 1
            if (.not. (all(a(i, j, :) == c%a(1, 1, :)) .and. all(v(i, j, :) == c%v(1, 1, :)) .and. &
                       all(ts(i, j, :) == c%ts(1, 1, :)) .and. all(q(i, j, :, :) == c%q(1, 1, :, :)) .and. &
                       all(s(i, j, :, :) == c%s(1, 1, :, :)) .and. a0(i, j) == c%a0(1, 1) .and. &
                       fr(i, j) == c%fr(1, 1) .and. fs(i, j) == c%fs(1, 1) .and. fz(i, j) == c%fz(1, 1) .and. &
                       ie(i, j) == c%ie(1, 1))) nbad = nbad + 1
         end do
      end do
      call check(error, nbad == 0 .and. nfrz > NX*NY/2, "in-grid cell == extracted column, every cell and output")
   end subroutine test_cell_vs_column

   subroutine test_driver_every_cell(error)
      !! T15 (device vs host): a random state in which EVERY cell differs goes through the step driver
      !! over an INTERIOR range (the GPU build runs it on the device); every output of every cell is
      !! compared with the host kernel on the same cell, the untouched ring (state, fluxes, bank,
      !! error word) is asserted bitwise, dry cells keep their bank, wet freezing cells spend it.
      !! Tolerance 1e-12 relative (FMA on the GPU build); bitwise on CPU builds.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 33, NY = 21, NIL = 7, IS = 3, IE_ = NX - 2, JS = 3, JE = NY - 2
      type(ice_params_t) :: p
      real(wp) :: frz(NX, NY), bank(NX, NY), tf(NX, NY), a(NX, NY, NC), v(NX, NY, NC), ts(NX, NY, NC)
      real(wp) :: q(NX, NY, NC, NIL), s(NX, NY, NC, NIL), a0(NX, NY), sal(NIL), fr(NX, NY), fs(NX, NY), fz(NX, NY)
      real(wp) :: ah(NX, NY, NC), vh(NX, NY, NC), tsh(NX, NY, NC), qh(NX, NY, NC, NIL), sh(NX, NY, NC, NIL)
      real(wp) :: a0h(NX, NY), frh(NX, NY), fsh(NX, NY), fzh(NX, NY), bankh(NX, NY), hv(0:NC)
      real(wp) :: a_in(NX, NY, NC), v_in(NX, NY, NC), q_in(NX, NY, NC, NIL), s_in(NX, NY, NC, NIL), a0_in(NX, NY)
      real(wp) :: tol, worst
      logical :: wet(NX, NY)
      integer(int32) :: ie(NX, NY), ieh(NX, NY)
      integer :: i, j, nspent, nchg, nring_bad, nbank_bad

#ifdef RDB_GPU_OFFLOAD
      tol = 1.0e-12_wp
#else
      tol = 0.0_wp
#endif
      p = params(ICE_CONSTANTS_ICEPACK, NIL, .true., 5.0_wp, 0.05_wp)
      call fill_grid(NX, NY, NIL, 31415_int64, p, frz, bank, tf, wet, a, v, ts, q, s, a0, sal)
      fr = 1.0e-6_wp; fs = -2.0e-9_wp; fz = 3.0e-4_wp; ie = 0_int32
      hv = 0.0_wp; hv(1) = HIN1
      a_in = a; v_in = v; q_in = q; s_in = s; a0_in = a0
      ah = a; vh = v; tsh = ts; qh = q; sh = s; a0h = a0; frh = fr; fsh = fs; fzh = fz; bankh = bank; ieh = ie
      do j = JS, JE
         do i = IS, IE_
            if (.not. wet(i, j)) cycle
            call ice_add_new_ice_cell(i, j, NX, NY, NC, NIL, DT, frz(i, j), tf(i, j), HIN1, p%c%rhoi, p%c%Lfresh, &
                                      p%c%puny, p%thermo%hfrazilmin, p%thermo%ice_ref_salinity, p%thermo%update_ocn_f, &
                                      sal, ah, vh, tsh, qh, sh, a0h, frh, fsh, fzh, ieh)
            if (frz(i, j) > 0.0_wp .and. ieh(i, j) == 0_int32) bankh(i, j) = 0.0_wp
         end do
      end do
      !$acc enter data copyin(a, v, ts, q, s, sal, a0, fr, fs, fz, frz, bank, tf, ie, wet)
      call ice_new_ice_step(NX, NY, NC, NIL, IS, IE_, JS, JE, DT, p, hv, sal, wet, frz, bank, tf, a, v, ts, q, s, &
                            a0, fr, fs, fz, ie)
      !$acc update self(a, v, ts, q, s, a0, fr, fs, fz, bank, ie)
      !$acc exit data delete(a, v, ts, q, s, sal, a0, fr, fs, fz, frz, bank, tf, ie, wet)
      worst = 0.0_wp
      nspent = 0; nchg = 0; nring_bad = 0; nbank_bad = 0
      do j = 1, NY
         do i = 1, NX
            if (i >= IS .and. i <= IE_ .and. j >= JS .and. j <= JE) then
               call cmp(a(i, j, :), ah(i, j, :))
               call cmp(v(i, j, :), vh(i, j, :))
               call cmp(ts(i, j, :), tsh(i, j, :))
               call cmp(reshape(q(i, j, :, :), [NC*NIL]), reshape(qh(i, j, :, :), [NC*NIL]))
               call cmp(reshape(s(i, j, :, :), [NC*NIL]), reshape(sh(i, j, :, :), [NC*NIL]))
               call cmp([a0(i, j), fr(i, j), fs(i, j), fz(i, j)], [a0h(i, j), frh(i, j), fsh(i, j), fzh(i, j)])
               if (ie(i, j) /= 0_int32 .or. ieh(i, j) /= 0_int32) worst = huge(1.0_wp)
               if (bank(i, j) /= bankh(i, j)) nbank_bad = nbank_bad + 1
               if (wet(i, j) .and. frz(i, j) > 0.0_wp) then
                  nspent = nspent + 1
                  if (.not. ieee_is_finite(fz(i, j)) .or. .not. (fz(i, j) > 3.0e-4_wp)) nchg = nchg + 1
               end if
               if (.not. wet(i, j)) then
                  if (fz(i, j) /= 3.0e-4_wp .or. bank(i, j) /= frz(i, j)*DT) nbank_bad = nbank_bad + 1
               end if
            else      ! ring: bitwise untouched
               if (.not. (all(a(i, j, :) == a_in(i, j, :)) .and. all(v(i, j, :) == v_in(i, j, :)) .and. &
                          all(q(i, j, :, :) == q_in(i, j, :, :)) .and. all(s(i, j, :, :) == s_in(i, j, :, :)) .and. &
                          a0(i, j) == a0_in(i, j) .and. fr(i, j) == 1.0e-6_wp .and. fs(i, j) == -2.0e-9_wp .and. &
                          fz(i, j) == 3.0e-4_wp .and. bank(i, j) == frz(i, j)*DT .and. ie(i, j) == 0_int32)) &
                  nring_bad = nring_bad + 1
            end if
         end do
      end do
      call check(error, worst <= tol, "driver output differs from the host kernel in some cell (or is non-finite)")
      if (allocated(error)) return
      call check(error, nring_bad == 0, "ring cells untouched")
      if (allocated(error)) return
      call check(error, nbank_bad == 0, "bank spent exactly where wet and freezing; dry cells keep it")
      if (allocated(error)) return
      call check(error, nspent > 100 .and. nchg == 0, "every wet freezing cell received frazil")
   contains
      subroutine cmp(d, h)
         real(wp), intent(in) :: d(:), h(:)
         integer :: m
         real(wp) :: x

         do m = 1, size(d)
            if (.not. (ieee_is_finite(d(m)) .and. ieee_is_finite(h(m)))) then
               worst = huge(1.0_wp)
            else
               x = abs(d(m) - h(m))/max(abs(h(m)), 1.0e-30_wp)
               if (h(m) == 0.0_wp) x = abs(d(m))
               if (x > worst) worst = x
            end if
         end do
      end subroutine cmp
   end subroutine test_driver_every_cell

   subroutine test_throughput(error)
      !! Driver throughput on 512 x 512 x ncat 5 x nilyr 7 (printed; the assertion is that every
      !! cell ran clean and the state stayed finite).  The GPU build times the device sweep.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 512, NY = 512, NIL = 7, NREP = 4
      type(ice_params_t) :: p
      real(wp), allocatable :: frz(:, :), bank(:, :), tf(:, :), a(:, :, :), v(:, :, :), ts(:, :, :)
      real(wp), allocatable :: q(:, :, :, :), s(:, :, :, :), a0(:, :), fr(:, :), fs(:, :), fz(:, :)
      real(wp) :: sal(NIL), hv(0:NC), secs, tot
      logical, allocatable :: wet(:, :)
      integer(int32), allocatable :: ie(:, :)
      integer(int64) :: t0, t1, rate
      integer :: rep, i, j, n

      p = params(ICE_CONSTANTS_ICEPACK, NIL, .true., 5.0_wp, 0.05_wp)
      allocate (frz(NX, NY), bank(NX, NY), tf(NX, NY), a(NX, NY, NC), v(NX, NY, NC), ts(NX, NY, NC))
      allocate (q(NX, NY, NC, NIL), s(NX, NY, NC, NIL), a0(NX, NY), fr(NX, NY), fs(NX, NY), fz(NX, NY))
      allocate (wet(NX, NY), ie(NX, NY))
      call fill_grid(NX, NY, NIL, 2718_int64, p, frz, bank, tf, wet, a, v, ts, q, s, a0, sal)
      wet = .true.
      frz = abs(frz)*0.01_wp                  ! keep the repeated sweeps from over-filling the cells
      bank = frz*DT
      fr = 0.0_wp; fs = 0.0_wp; fz = 0.0_wp; ie = 0_int32
      hv = 0.0_wp; hv(1) = HIN1
      !$acc enter data copyin(a, v, ts, q, s, sal, a0, fr, fs, fz, frz, bank, tf, ie, wet)
      call ice_new_ice_step(NX, NY, NC, NIL, 1, NX, 1, NY, DT, p, hv, sal, wet, frz, bank, tf, a, v, ts, q, s, &
                            a0, fr, fs, fz, ie)                                   ! warm-up
      call system_clock(t0, rate)
      do rep = 1, NREP
         call ice_new_ice_step(NX, NY, NC, NIL, 1, NX, 1, NY, DT, p, hv, sal, wet, frz, bank, tf, a, v, ts, q, s, &
                               a0, fr, fs, fz, ie)
      end do
      !$acc wait
      call system_clock(t1)
      !$acc update self(a, v, ts, q, s, a0, fz, ie)
      !$acc exit data delete(a, v, ts, q, s, sal, a0, fr, fs, fz, frz, bank, tf, ie, wet)
      secs = real(t1 - t0, wp)/real(rate, wp)/real(NREP, wp)
      write (output_unit, '(a,f9.3,a,es10.3,a,es10.3,a)') "  new ice driver sweep (512x512, ncat 5, nilyr 7): ", secs*1.0e3_wp, &
         " ms/sweep, ", real(NX*NY, wp)/secs, " cells/s, ", real(NX*NY*NC, wp)/secs, " cell-categories/s"
      tot = 0.0_wp
      do n = 1, NC
         do j = 1, NY
            do i = 1, NX
               tot = tot + v(i, j, n)
            end do
         end do
      end do
      call check(error, all(ie == 0_int32) .and. ieee_is_finite(tot) .and. ieee_is_finite(sum(a0)), &
                 "every cell ran clean and finite")
   end subroutine test_throughput

end module test_cice_new_ice
