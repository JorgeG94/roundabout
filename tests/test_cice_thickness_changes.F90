module test_cice_thickness_changes
   !! Tests of `ice_thickness_changes`, `ice_conservation_check_vthermo` and their
   !! cell entries (A8, A9): closed forms written here (Stefan top melt of snow and
   !! the melt cascade, congelation rate from `fcondbot - fbot`, melt-out returning
   !! the residual energy, snow melted from below, sublimation of snow then ice,
   !! condensation, snowfall, snow-ice, the fresh-ice option), the no-flux
   !! no-op, every error bit and NaN-in before any write, the exact per-column
   !! closure of energy, thickness, water and salt over 1e5 random columns per
   !! layer count and constant set (the ocean fluxes composed with the hand-off
   !! formulae of the thermodynamic driver), goldens from the real upstream routine
   !! on 66 columns per constant set, and a category-parallel `do concurrent`
   !! device test over mapped `(nx, ny, ncat, nlyr)` arrays comparing EVERY
   !! output on EVERY cell with a host sweep (with the throughput).
   use, intrinsic :: iso_fortran_env, only: int32, int64, output_unit
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan, ieee_positive_inf
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use rdb_constants, only: wp
   use lcg_deterministic, only: lcg_next, lcg_unit
   use rdb_ice_errors, only: ICE_ERR_NONFINITE, ICE_ERR_THERMO_CONSERV, ICE_ERR_TBOT, ICE_ERR_NEG_THICK
   use rdb_ice_params, only: ice_constants_t, ice_constants, ICE_CONSTANTS_ICEPACK, ICE_CONSTANTS_CESM, &
                             ice_q_from_T, ice_qsno_from_T
   use rdb_ice_thickness, only: ice_thickness_changes, ice_thickness_changes_cell, ice_conservation_check_vthermo, &
                                ice_conservation_check_vthermo_cell, ICE_FERRMAX
   use test_cice_util, only: ice_ledger_t
   implicit none
   private

   public :: collect_cice_thickness_changes_tests

   real(wp), parameter :: DT = 900.0_wp
   real(wp), parameter :: SREF = 4.0_wp     !! reference salinity of the salt-flux hand-off [ppt]
   real(wp), parameter :: HPOW2 = 0.5_wp !! 2**-1: layer sums and ratios are exact in binary
   integer, parameter :: NIL_SET(5) = [1, 2, 4, 7, 8]

#ifdef RDB_GPU_OFFLOAD
   real(wp), parameter :: CRT = 1.0e-11_wp
      !! Closed-form relative tolerance of the leftover-energy bookkeeping, GPU build (FMA on).
   real(wp), parameter :: GRTOL = 1.0e-12_wp
      !! Golden and device-vs-host relative tolerance, GPU build (FMA on, host code too).
   real(wp), parameter :: GATOL = 1.0e-300_wp
   real(wp), parameter :: DRTOL = 1.0e-12_wp
#else
   real(wp), parameter :: CRT = 1.0e-13_wp
      !! Closed-form relative tolerance of the leftover-energy bookkeeping, CPU builds.
   real(wp), parameter :: GRTOL = 0.0_wp
      !! CPU builds: exact against the oracle (measured worst: stated in the report).
   real(wp), parameter :: GATOL = 0.0_wp
   real(wp), parameter :: DRTOL = 0.0_wp
#endif

   type :: column_t
      !! One column with its inputs, the entry bookkeeping and the outputs.
      integer :: nilyr = 4
      logical :: l_brine = .true.
      real(wp) :: dt = DT, Tbot = -1.8_wp, fbot = 0.0_wp, flatn = 0.0_wp, fsurfn = 0.0_wp
      real(wp) :: fcondtopn = 0.0_wp, fcondbotn = 0.0_wp, fsnow = 0.0_wp, zSin_bot = 0.0_wp
      real(wp) :: hin = 0.0_wp, hsn = 0.0_wp, hilyr = 0.0_wp, hslyr = 0.0_wp
      real(wp) :: zqin(8) = 0.0_wp, zqsn(1) = 0.0_wp
      real(wp) :: meltt = 0.0_wp, melts = 0.0_wp, meltb = 0.0_wp, congel = 0.0_wp, snoice = 0.0_wp, dsnow = 0.0_wp
      real(wp) :: hsn_new = 0.0_wp, fhocnn = 0.0_wp, evapn = 0.0_wp, evapsn = 0.0_wp, evapin = 0.0_wp, efinal = 0.0_wp
      real(wp) :: hin0 = 0.0_wp, hsn0 = 0.0_wp, einter = 0.0_wp   !! entry thicknesses and column energy
      integer(int32) :: ierr = 0_int32
   end type column_t

contains

   subroutine collect_cice_thickness_changes_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("snow_top_melt_stefan", test_snow_top_melt), &
                  new_unittest("top_melt_cascade_through_snow_into_ice", test_top_cascade), &
                  new_unittest("congelation_rate_and_enthalpy", test_growth), &
                  new_unittest("bottom_melt_rate", test_bottom_melt), &
                  new_unittest("melt_out_returns_residual_energy_to_ocean", test_melt_out), &
                  new_unittest("snow_melted_from_below", test_snow_from_below), &
                  new_unittest("sublimation_snow_then_ice", test_sublimation), &
                  new_unittest("condensation_on_snow_and_on_ice", test_condensation), &
                  new_unittest("snowfall_and_snow_ice_formation", test_snowfall), &
                  new_unittest("no_flux_column_is_a_noop_but_relayering", test_noop), &
                  new_unittest("fresh_ice_option_l_brine_false", test_fresh_option), &
                  new_unittest("error_bits", test_error_bits), &
                  new_unittest("raw_snow_free_category_is_valid", test_raw_snow_free), &
                  new_unittest("nan_in_flagged_before_any_write", test_nan_in), &
                  new_unittest("conservation_check_threshold", test_check_threshold), &
                  new_unittest("random_columns_closure_1e5", test_closure_random), &
                  new_unittest("goldens_vs_upstream_both_sets", test_golden), &
                  new_unittest("cell_entry_matches_column", test_cell_vs_column), &
                  new_unittest("device_cell_sweep_matches_host", test_device) &
                  ]
   end subroutine collect_cice_thickness_changes_tests

   ! ------------------------------------------------------------------
   ! helpers
   ! ------------------------------------------------------------------

   function energy(col) result(e)
      !! Column energy [J/m2]: snow then ice, the order of the upstream `efinal`.
      type(column_t), intent(in) :: col
      real(wp) :: e
      integer :: k
      e = 0.0_wp
      e = e + col%hslyr*col%zqsn(1)
      do k = 1, col%nilyr
         e = e + col%hilyr*col%zqin(k)
      end do
   end function energy

   subroutine base_col(col, c, nilyr, hin, hsn, Tice, Tsnow)
      !! A column of uniform layers at `Tice`, zero salinity, snow at `Tsnow`, no fluxes.
      type(column_t), intent(out) :: col
      type(ice_constants_t), intent(in) :: c
      integer, intent(in) :: nilyr
      real(wp), intent(in) :: hin, hsn, Tice, Tsnow
      integer :: k
      col%nilyr = nilyr
      col%hin = hin
      col%hsn = hsn
      col%hilyr = hin/real(nilyr, wp)
      col%hslyr = hsn
      do k = 1, nilyr
         col%zqin(k) = ice_q_from_T(Tice, 0.0_wp, c, .true.)
      end do
      col%zqsn(1) = ice_qsno_from_T(Tsnow, c)
   end subroutine base_col

   subroutine run(col, c)
      type(column_t), intent(inout) :: col
      type(ice_constants_t), intent(in) :: c
      col%hin0 = col%hin
      col%hsn0 = col%hsn
      col%einter = energy(col)
      call ice_thickness_changes(col%nilyr, 1, c, col%l_brine, col%dt, col%Tbot, col%fbot, col%flatn, col%fsurfn, &
                                 col%fcondtopn, col%fcondbotn, col%fsnow, col%zSin_bot, col%hin, col%hilyr, col%hsn, &
                                 col%hslyr, col%zqin(1:col%nilyr), col%zqsn, col%meltt, col%melts, col%meltb, &
                                 col%congel, col%snoice, col%dsnow, col%hsn_new, col%fhocnn, col%evapn, col%evapsn, &
                                 col%evapin, col%efinal, col%ierr)
   end subroutine run

   function energy_residual(col, c) result(r)
      !! `(efinal - einter) - (fsurfn - fcondtopn - flatn + fcondbotn - fhocnn - fsnow Lfresh) dt` [J/m2]
      type(column_t), intent(in) :: col
      type(ice_constants_t), intent(in) :: c
      real(wp) :: r
      r = (col%efinal - col%einter) - (col%fsurfn - col%fcondtopn - col%flatn + col%fcondbotn - col%fhocnn &
                                       - col%fsnow*c%Lfresh)*col%dt
   end function energy_residual

   subroutine hand_off(col, c, fresh, fsalt)
      !! The ocean fluxes the thermodynamic driver books after the thickness change.
      type(column_t), intent(in) :: col
      type(ice_constants_t), intent(in) :: c
      real(wp), intent(out) :: fresh, fsalt
      real(wp) :: dhi, dhs
      dhi = col%hin - col%hin0
      dhs = col%hsn - col%hsn0 - col%hsn_new
      fresh = col%evapn - (c%rhoi*dhi + c%rhos*dhs)/col%dt
      fsalt = -c%rhoi*dhi*SREF*0.001_wp/col%dt
   end subroutine hand_off

   function near(a, b, rtol, atol) result(ok)
      real(wp), intent(in) :: a, b, rtol, atol
      logical :: ok
      ok = .false.
      if (.not. (ieee_is_finite(a) .and. ieee_is_finite(b))) return
      ok = abs(a - b) <= atol + rtol*max(abs(a), abs(b))
   end function near

   function uniform(st) result(u)
      integer(int64), intent(inout) :: st
      real(wp) :: u
      st = lcg_next(st)
      u = lcg_unit(st)
   end function uniform

   ! ------------------------------------------------------------------
   ! closed forms
   ! ------------------------------------------------------------------

   subroutine test_snow_top_melt(error)
      !! Stefan: `etop_mlt = (fsurfn - fcondtopn) dt` melts `etop/|zqsn|` of snow, nothing else.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      type(column_t) :: col
      integer :: iset, ni, s
      real(wp) :: etop, expect

      do iset = ICE_CONSTANTS_ICEPACK, ICE_CONSTANTS_CESM
         c = ice_constants(iset)
         do s = 1, 2
            ni = merge(1, 4, s == 1)
            call base_col(col, c, ni, 1.0_wp, 0.2_wp, -10.0_wp, -5.0_wp)
            col%fsurfn = 100.0_wp
            col%fcondtopn = 0.0_wp
            etop = 100.0_wp*DT
            expect = etop/(-col%zqsn(1))
            call run(col, c)
            call check(error, col%ierr == 0_int32, "no error flag")
            if (allocated(error)) return
            call check(error, near(col%melts, expect, 1.0e-14_wp, 0.0_wp), "melts = etop/|q_snow|")
            if (allocated(error)) return
            call check(error, near(col%hsn, 0.2_wp - expect, 1.0e-14_wp, 0.0_wp), "snow thickness")
            if (allocated(error)) return
            call check(error, abs(col%meltt) < 1.0e-12_wp .and. col%meltb == 0.0_wp .and. col%congel == 0.0_wp, "ice untouched")
            if (allocated(error)) return
            call check(error, near(col%hin, 1.0_wp, 1.0e-14_wp, 0.0_wp), "ice thickness unchanged")
            if (allocated(error)) return
            call check(error, abs(col%fhocnn) < 1.0e-9_wp, "no energy left over for the ocean")
            if (allocated(error)) return
            call check(error, abs(energy_residual(col, c)) < 1.0e-5_wp, "energy closes")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_snow_top_melt

   subroutine test_top_cascade(error)
      !! Melt energy consumed layer by layer: all the snow, then 1.5 layers of ice, which keeps its enthalpy.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      type(column_t) :: col
      integer :: iset, k
      real(wp) :: etop, qi

      do iset = ICE_CONSTANTS_ICEPACK, ICE_CONSTANTS_CESM
         c = ice_constants(iset)
         call base_col(col, c, 4, 0.8_wp, 0.1_wp, -8.0_wp, -5.0_wp)
         qi = col%zqin(1)
         etop = 0.1_wp*(-col%zqsn(1)) + 1.5_wp*0.2_wp*(-qi)
         col%fsurfn = etop/DT
         call run(col, c)
         call check(error, col%ierr == 0_int32, "no error flag")
         if (allocated(error)) return
         call check(error, near(col%melts, 0.1_wp, 1.0e-13_wp, 0.0_wp), "all snow melted")
         if (allocated(error)) return
         call check(error, col%hsn == 0.0_wp, "snow thickness exactly zero")
         if (allocated(error)) return
         call check(error, near(col%meltt, 0.3_wp, 1.0e-12_wp, 0.0_wp), "1.5 layers of ice melted from the top")
         if (allocated(error)) return
         call check(error, near(col%hin, 0.5_wp, 1.0e-12_wp, 0.0_wp), "ice thickness")
         if (allocated(error)) return
         do k = 1, 4
            call check(error, near(col%zqin(k), qi, 1.0e-13_wp, 0.0_wp), "remaining layers keep their enthalpy")
            if (allocated(error)) return
         end do
         call check(error, abs(col%fhocnn) < 1.0e-6_wp, "no residual energy")
         if (allocated(error)) return
         call check(error, abs(energy_residual(col, c)) < 1.0e-4_wp, "energy closes")
         if (allocated(error)) return
      end do
   end subroutine test_top_cascade

   subroutine test_growth(error)
      !! Congelation: `dhi = (fcondbot - fbot) dt/qbot` with the BL99 growth enthalpy `qbot`.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      type(column_t) :: col
      integer :: iset, s, ni
      real(wp) :: qbot, dhi, Tmlts, fresh, fsalt, qmax

      do iset = ICE_CONSTANTS_ICEPACK, ICE_CONSTANTS_CESM
         c = ice_constants(iset)
         do s = 1, 6
            ni = NIL_SET(1 + mod(s - 1, 5))
            call base_col(col, c, ni, 0.6_wp, 0.0_wp, -8.0_wp, -5.0_wp)
            col%hslyr = 0.0_wp
            col%fcondbotn = -15.0_wp
            col%fbot = -5.0_wp
            col%Tbot = -1.8_wp
            col%zSin_bot = 0.0_wp
            if (s == 6) then
               col%zSin_bot = 3.2_wp
               col%Tbot = -1.7_wp
            end if
            Tmlts = -col%zSin_bot*c%depressT
            qbot = -c%rhoi*(c%cp_ice*(Tmlts - col%Tbot) + c%Lfresh*(1.0_wp - Tmlts/col%Tbot) - c%cp_ocn*Tmlts)
            dhi = (col%fcondbotn - col%fbot)*DT/qbot
            call run(col, c)
            call check(error, col%ierr == 0_int32, "no error flag")
            if (allocated(error)) return
            call check(error, near(col%congel, dhi, 1.0e-14_wp, 0.0_wp) .and. dhi > 0.0_wp, "congel = growth rate x dt")
            if (allocated(error)) return
            call check(error, near(col%hin, 0.6_wp + dhi, 1.0e-14_wp, 0.0_wp), "ice thickness")
            if (allocated(error)) return
            call check(error, near(col%fhocnn, col%fbot, 0.0_wp, 1.0e-12_wp), "fhocn = fbot (nothing left over)")
            if (allocated(error)) return
            call check(error, near(col%efinal - col%einter, dhi*qbot, 1.0e-12_wp, 1.0e-3_wp), "new ice enters at qbot")
            if (allocated(error)) return
            call check(error, abs(energy_residual(col, c)) < 1.0e-4_wp, "energy closes")
            if (allocated(error)) return
            call hand_off(col, c, fresh, fsalt)
            call check(error, near(fresh, -c%rhoi*dhi/DT, 1.0e-9_wp, 0.0_wp), "fresh = -rhoi dhi/dt")
            if (allocated(error)) return
            call check(error, near(fsalt, -c%rhoi*dhi*SREF*0.001_wp/DT, 1.0e-9_wp, 0.0_wp), "fsalt")
            if (allocated(error)) return
         end do
         ! The growth enthalpy is capped at -rhoi Lfresh/2 when the base is barely below its melting point.
         call base_col(col, c, 4, 0.6_wp, 0.0_wp, -8.0_wp, -5.0_wp)
         col%hslyr = 0.0_wp
         col%zSin_bot = 3.2_wp
         col%Tbot = -3.2_wp*c%depressT - 1.0e-5_wp
         col%fcondbotn = -10.0_wp
         qmax = -0.5_wp*c%rhoi*c%Lfresh
         dhi = col%fcondbotn*DT/qmax
         call run(col, c)
         call check(error, col%ierr == 0_int32, "capped: no error flag")
         if (allocated(error)) return
         call check(error, near(col%congel, dhi, 1.0e-14_wp, 0.0_wp), "capped growth enthalpy")
         if (allocated(error)) return
      end do
   end subroutine test_growth

   subroutine test_bottom_melt(error)
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      type(column_t) :: col
      integer :: iset, ni
      real(wp) :: dhi

      do iset = ICE_CONSTANTS_ICEPACK, ICE_CONSTANTS_CESM
         c = ice_constants(iset)
         do ni = 1, 4, 3
            call base_col(col, c, ni, 0.8_wp, 0.0_wp, -8.0_wp, -5.0_wp)
            col%hslyr = 0.0_wp
            col%fbot = -200.0_wp
            dhi = 200.0_wp*DT/(-col%zqin(1))
            call run(col, c)
            call check(error, col%ierr == 0_int32, "no error flag")
            if (allocated(error)) return
            call check(error, near(col%meltb, dhi, 1.0e-13_wp, 0.0_wp), "meltb = ocean heat/|q|")
            if (allocated(error)) return
            call check(error, near(col%hin, 0.8_wp - dhi, 1.0e-13_wp, 0.0_wp), "ice thickness")
            if (allocated(error)) return
            call check(error, abs(col%fhocnn - col%fbot) < 1.0e-8_wp, "fhocn = fbot")
            if (allocated(error)) return
            call check(error, abs(energy_residual(col, c)) < 1.0e-4_wp, "energy closes")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_bottom_melt

   subroutine test_melt_out(error)
      !! Ice thinner than the melt energy: nothing is left, and the rest goes to the ocean.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      type(column_t) :: col
      integer :: iset, ni, k
      real(wp) :: eice, expect

      do iset = ICE_CONSTANTS_ICEPACK, ICE_CONSTANTS_CESM
         c = ice_constants(iset)
         do ni = 1, 4, 3
            call base_col(col, c, ni, 0.001_wp, 0.0_wp, -8.0_wp, -5.0_wp)
            col%hslyr = 0.0_wp
            col%fbot = -1000.0_wp
            eice = 0.001_wp*(-col%zqin(1))
            expect = col%fbot + (1000.0_wp*DT - eice)/DT
            call run(col, c)
            call check(error, col%ierr == 0_int32, "no error flag")
            if (allocated(error)) return
            call check(error, col%hin == 0.0_wp .and. col%hilyr == 0.0_wp, "ice thickness exactly zero")
            if (allocated(error)) return
            call check(error, near(col%meltb, 0.001_wp, 1.0e-13_wp, 0.0_wp), "all the ice melted from below")
            if (allocated(error)) return
            call check(error, near(col%fhocnn, expect, 1.0e-12_wp, 1.0e-9_wp), "residual energy returns to the ocean")
            if (allocated(error)) return
            call check(error, col%fhocnn > col%fbot, "fhocn above fbot")
            if (allocated(error)) return
            do k = 1, ni
               call check(error, col%zqin(k) == 0.0_wp, "no enthalpy left in an empty column")
               if (allocated(error)) return
            end do
            call check(error, abs(energy_residual(col, c)) < 1.0e-5_wp, "energy closes")
            if (allocated(error)) return
         end do
         ! Top melt through the whole column: the same bookkeeping from above.
         call base_col(col, c, 4, 0.001_wp, 0.0_wp, -8.0_wp, -5.0_wp)
         col%hslyr = 0.0_wp
         eice = 0.001_wp*(-col%zqin(1))
         col%fsurfn = (eice + 5000.0_wp)/DT
         call run(col, c)
         call check(error, col%hin == 0.0_wp .and. near(col%meltt, 0.001_wp, 1.0e-13_wp, 0.0_wp), "top melt-out")
         if (allocated(error)) return
         call check(error, near(col%fhocnn, 5000.0_wp/DT, 1.0e-9_wp, 1.0e-9_wp), "residual top-melt energy to the ocean")
         if (allocated(error)) return
         call check(error, abs(energy_residual(col, c)) < 1.0e-5_wp, "top melt-out energy closes")
         if (allocated(error)) return
      end do
   end subroutine test_melt_out

   subroutine test_snow_from_below(error)
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      type(column_t) :: col
      integer :: iset
      real(wp) :: eice, esno, dhsn

      do iset = ICE_CONSTANTS_ICEPACK, ICE_CONSTANTS_CESM
         c = ice_constants(iset)
         call base_col(col, c, 4, 0.001_wp, 0.01_wp, -8.0_wp, -5.0_wp)
         eice = 0.001_wp*(-col%zqin(1))
         esno = 0.01_wp*(-col%zqsn(1))
         col%fbot = -(eice + 0.6_wp*esno)/DT
         call run(col, c)
         call check(error, col%ierr == 0_int32, "no error flag")
         if (allocated(error)) return
         call check(error, near(col%melts, 0.006_wp, 1.0e-10_wp, 0.0_wp), "60 % of the snow melted from below")
         if (allocated(error)) return
         ! upstream quirk kept: the 4 mm of snow left on a column with no ice is below the waterline and becomes ice
         dhsn = 0.004_wp*c%rhoi/c%rhow
call check(error, near(col%snoice, dhsn*c%rhos/c%rhoi, 1.0e-8_wp, 0.0_wp) .and. near(col%hsn, 0.004_wp - dhsn, 1.0e-8_wp, 0.0_wp), &
                    "the snow left floats on no ice and is converted (upstream behaviour)")
         if (allocated(error)) return
         call check(error, abs(col%fhocnn - col%fbot) < 1.0e-6_wp, "no leftover")
         if (allocated(error)) return
         call check(error, abs(energy_residual(col, c)) < 1.0e-4_wp, "energy closes")
         if (allocated(error)) return
      end do
   end subroutine test_snow_from_below

   subroutine test_sublimation(error)
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      type(column_t) :: col
      integer :: iset
      real(wp) :: qsub_s, qsub_i, dhs, dhi, esub, fresh, fsalt, esub_snow

      do iset = ICE_CONSTANTS_ICEPACK, ICE_CONSTANTS_CESM
         c = ice_constants(iset)
         ! snow only
         call base_col(col, c, 4, 1.0_wp, 0.2_wp, -8.0_wp, -5.0_wp)
         qsub_s = col%zqsn(1) - c%rhos*c%Lvap
         col%flatn = -50.0_wp
         dhs = 50.0_wp*DT/qsub_s
         call run(col, c)
         call check(error, col%ierr == 0_int32, "no error flag")
         if (allocated(error)) return
         call check(error, near(col%hsn, 0.2_wp + dhs, 1.0e-13_wp, 0.0_wp), "snow thinned by esub/(q - rhos Lvap)")
         if (allocated(error)) return
         call check(error, near(col%evapsn, c%rhos*dhs/DT, CRT, 0.0_wp) .and. col%evapn == col%evapsn .and. &
                    abs(col%evapin) < 1.0e-12_wp*abs(col%evapn) .and. col%evapn < 0.0_wp, "evapn over snow")
         if (allocated(error)) return
         call check(error, abs(energy_residual(col, c)) < 1.0e-3_wp, "energy closes")
         if (allocated(error)) return
         call hand_off(col, c, fresh, fsalt)
         call check(error, abs(fresh) < 1.0e-10_wp*abs(col%evapn), "pure sublimation books no fresh water")
         if (allocated(error)) return
         ! snow, then ice
         call base_col(col, c, 4, 1.0_wp, 0.01_wp, -8.0_wp, -5.0_wp)
         qsub_s = col%zqsn(1) - c%rhos*c%Lvap
         qsub_i = col%zqin(1) - c%rhoi*c%Lvap
         esub_snow = 0.01_wp*(-qsub_s)
         esub = esub_snow + 0.05_wp*(-qsub_i)
         col%flatn = -esub/DT
         dhi = -0.05_wp
         call run(col, c)
         call check(error, col%ierr == 0_int32, "snow+ice: no error flag")
         if (allocated(error)) return
         call check(error, col%hsn == 0.0_wp, "snow sublimated first")
         if (allocated(error)) return
         call check(error, near(col%hin, 1.0_wp + dhi, 1.0e-12_wp, 0.0_wp), "the remainder sublimates ice")
         if (allocated(error)) return
         call check(error, near(col%evapin, c%rhoi*dhi/DT, 1.0e-10_wp, 0.0_wp) .and. &
                    near(col%evapsn, -c%rhos*0.01_wp/DT, 1.0e-12_wp, 0.0_wp), "evapin, evapsn")
         if (allocated(error)) return
         call check(error, abs(energy_residual(col, c)) < 1.0e-3_wp, "snow+ice: energy closes")
         if (allocated(error)) return
         call hand_off(col, c, fresh, fsalt)
         call check(error, abs(fresh) < 1.0e-8_wp*abs(col%evapn), "snow+ice: no fresh water")
         if (allocated(error)) return
         ! bare ice
         call base_col(col, c, 7, 1.0_wp, 0.0_wp, -8.0_wp, -5.0_wp)
         col%hslyr = 0.0_wp
         qsub_i = col%zqin(1) - c%rhoi*c%Lvap
         col%flatn = -30.0_wp
         dhi = 30.0_wp*DT/qsub_i
         call run(col, c)
         call check(error, col%ierr == 0_int32 .and. near(col%hin, 1.0_wp + dhi, 1.0e-13_wp, 0.0_wp), "bare-ice sublimation")
         if (allocated(error)) return
         call check(error, abs(energy_residual(col, c)) < 1.0e-3_wp, "bare ice: energy closes")
         if (allocated(error)) return
      end do
   end subroutine test_sublimation

   subroutine test_condensation(error)
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      type(column_t) :: col
      integer :: iset
      real(wp) :: dhs, dhi, fresh, fsalt

      do iset = ICE_CONSTANTS_ICEPACK, ICE_CONSTANTS_CESM
         c = ice_constants(iset)
         call base_col(col, c, 4, 1.0_wp, 0.2_wp, -8.0_wp, -5.0_wp)
         col%flatn = 20.0_wp
         dhs = -20.0_wp*DT/(col%zqsn(1) - c%rhos*c%Lvap)
         call run(col, c)
         call check(error, col%ierr == 0_int32, "no error flag")
         if (allocated(error)) return
         call check(error, dhs > 0.0_wp .and. near(col%hsn, 0.2_wp + dhs, 1.0e-13_wp, 0.0_wp), "snow grows")
         if (allocated(error)) return
         call check(error, col%evapn > 0.0_wp .and. near(col%evapsn, c%rhos*dhs/DT, 1.0e-13_wp, 0.0_wp), "evapn > 0")
         if (allocated(error)) return
         call check(error, abs(energy_residual(col, c)) < 1.0e-3_wp, "energy closes")
         if (allocated(error)) return
         call hand_off(col, c, fresh, fsalt)
         call check(error, abs(fresh) < 1.0e-10_wp*abs(col%evapn), "condensation books no fresh water")
         if (allocated(error)) return
         call base_col(col, c, 4, 1.0_wp, 0.0_wp, -8.0_wp, -5.0_wp)
         col%hslyr = 0.0_wp
         col%flatn = 20.0_wp
         dhi = -20.0_wp*DT/(col%zqin(1) - c%rhoi*c%Lvap)
         call run(col, c)
         call check(error, col%ierr == 0_int32 .and. near(col%hin, 1.0_wp + dhi, 1.0e-13_wp, 0.0_wp), "ice grows")
         if (allocated(error)) return
         call check(error, near(col%evapin, c%rhoi*dhi/DT, 1.0e-13_wp, 0.0_wp) .and. col%evapsn == 0.0_wp, "evapin")
         if (allocated(error)) return
         call check(error, abs(energy_residual(col, c)) < 1.0e-3_wp, "bare ice: energy closes")
         if (allocated(error)) return
      end do
   end subroutine test_condensation

   subroutine test_snowfall(error)
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      type(column_t) :: col
      integer :: iset
      real(wp) :: hnew, zq, wk1, dhsn, dhin, fresh, fsalt

      do iset = ICE_CONSTANTS_ICEPACK, ICE_CONSTANTS_CESM
         c = ice_constants(iset)
         ! snow on thick ice: new snow at 0 degC, mixed by thickness
         call base_col(col, c, 4, 2.0_wp, 0.1_wp, -8.0_wp, -5.0_wp)
         col%fsnow = 1.0e-4_wp
         hnew = col%fsnow*DT/c%rhos
         zq = (0.1_wp*col%zqsn(1) + hnew*(-c%rhos*c%Lfresh))/(0.1_wp + hnew)
         call run(col, c)
         call check(error, col%ierr == 0_int32, "no error flag")
         if (allocated(error)) return
         call check(error, near(col%hsn_new, hnew, 1.0e-14_wp, 0.0_wp) .and. near(col%hsn, 0.1_wp + hnew, 1.0e-14_wp, 0.0_wp), &
                    "hsn_new = fsnow dt/rhos")
         if (allocated(error)) return
         call check(error, near(col%zqsn(1), zq, 1.0e-13_wp, 0.0_wp), "snow enthalpy mixed with the new snow at 0 degC")
         if (allocated(error)) return
         call check(error, near(col%dsnow, hnew, 1.0e-12_wp, 1.0e-15_wp), "dsnow")
         if (allocated(error)) return
         call check(error, abs(energy_residual(col, c)) < 1.0e-3_wp, "energy closes with the snowfall latent term")
         if (allocated(error)) return
         call hand_off(col, c, fresh, fsalt)
         call check(error, abs(fresh) < 1.0e-10_wp*c%rhos*hnew/DT, "new snow is excluded from the ocean fresh water")
         if (allocated(error)) return
         ! snow-ice: 0.136 m of new snow on 0.3 m of ice sinks the surface below the waterline
         call base_col(col, c, 4, 0.3_wp, 0.0_wp, -8.0_wp, -5.0_wp)
         col%hslyr = 0.0_wp
         col%fsnow = 0.05_wp
         hnew = col%fsnow*DT/c%rhos
         wk1 = hnew - 0.3_wp*(c%rhow - c%rhoi)/c%rhos
         dhsn = wk1*c%rhoi/c%rhow
         dhin = dhsn*c%rhos/c%rhoi
         call run(col, c)
         call check(error, col%ierr == 0_int32, "snow-ice: no error flag")
         if (allocated(error)) return
         call check(error, wk1 > 0.0_wp .and. near(col%snoice, dhin, 1.0e-13_wp, 0.0_wp), "snoice = dhin")
         if (allocated(error)) return
         call check(error, near(col%hin, 0.3_wp + dhin, 1.0e-13_wp, 0.0_wp) .and. near(col%hsn, hnew - dhsn, 1.0e-12_wp, 0.0_wp), &
                    "thicknesses after freeboard")
         if (allocated(error)) return
         call check(error, abs(col%hsn - col%hin*(c%rhow - c%rhoi)/c%rhos) < 1.0e-13_wp, "surface on the waterline")
         if (allocated(error)) return
         call check(error, abs(energy_residual(col, c)) < 1.0e-3_wp, "snow-ice: energy closes")
         if (allocated(error)) return
         call hand_off(col, c, fresh, fsalt)
         call check(error, abs(fresh) < 1.0e-10_wp*c%rhos*hnew/DT, "snow-ice books no fresh water")
         if (allocated(error)) return
         call check(error, near(fsalt, -c%rhoi*dhin*SREF*0.001_wp/DT, 1.0e-10_wp, 0.0_wp), "snow-ice books salt leaving (upstream)")
         if (allocated(error)) return
      end do
   end subroutine test_snowfall

   subroutine test_noop(error)
      !! Zero fluxes: the thicknesses and the accumulators are bitwise unchanged, and with layers of a power of
      !! two the relayering is exact too; for general layers the enthalpies agree to a few ulp.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      type(column_t) :: col, ref
      integer :: iset, ni, k, s

      do iset = ICE_CONSTANTS_ICEPACK, ICE_CONSTANTS_CESM
         c = ice_constants(iset)
         do s = 1, 6
            ni = NIL_SET(1 + mod(s - 1, 5))
            if (s == 6) ni = 7
            call base_col(col, c, ni, 0.9_wp, 0.05_wp, -8.0_wp, -5.0_wp)
            if (s <= 5) then
               col%hilyr = HPOW2
               col%hin = real(ni, wp)*HPOW2
            end if
            do k = 1, ni
               col%zqin(k) = col%zqin(k) + real(k, wp)*1.0e6_wp
            end do
            col%zSin_bot = 3.2_wp
            ref = col
            call run(col, c)
            call check(error, col%ierr == 0_int32, "no error flag")
            if (allocated(error)) return
            if (s <= 5) then
               call check(error, col%hin == ref%hin .and. col%hsn == ref%hsn .and. col%hilyr == ref%hilyr .and. &
                          col%hslyr == ref%hslyr, "thicknesses bitwise unchanged")
            else
               call check(error, near(col%hin, ref%hin, 1.0e-15_wp, 0.0_wp) .and. col%hsn == ref%hsn .and. &
                          near(col%hilyr, ref%hilyr, 1.0e-15_wp, 0.0_wp), "thicknesses unchanged to an ulp")
            end if
            if (allocated(error)) return
            call check(error, col%meltt == 0.0_wp .and. col%melts == 0.0_wp .and. col%meltb == 0.0_wp .and. &
                       col%congel == 0.0_wp .and. col%snoice == 0.0_wp .and. col%dsnow == 0.0_wp, "no history")
            if (allocated(error)) return
            call check(error, col%fhocnn == 0.0_wp .and. col%evapn == 0.0_wp .and. col%hsn_new == 0.0_wp, "no fluxes")
            if (allocated(error)) return
            do k = 1, ni
               if (s <= 5) then
                  call check(error, col%zqin(k) == ref%zqin(k), "enthalpy bitwise unchanged (power-of-two layers)")
               else
                  call check(error, near(col%zqin(k), ref%zqin(k), 1.0e-15_wp, 0.0_wp), "enthalpy unchanged to a few ulp")
               end if
               if (allocated(error)) return
            end do
            call check(error, col%zqsn(1) == ref%zqsn(1), "snow enthalpy bitwise unchanged")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_noop

   subroutine test_fresh_option(error)
      !! `l_brine = .false.`: layers above 0 degC melt back to 0 degC (thickness), energy conserved.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      type(column_t) :: col
      integer :: iset
      real(wp) :: Ti, dhi, Ts, dhs, qbot, dzi

      do iset = ICE_CONSTANTS_ICEPACK, ICE_CONSTANTS_CESM
         c = ice_constants(iset)
         call base_col(col, c, 4, 1.0_wp, 0.2_wp, -8.0_wp, -5.0_wp)
         col%l_brine = .false.
         col%zqin(1) = -c%rhoi*c%Lfresh + c%rhoi*c%cp_ice*0.5_wp
         col%zqsn(1) = -c%rhos*c%Lfresh + c%rhos*c%cp_ice*0.3_wp
         Ti = (c%Lfresh + col%zqin(1)/c%rhoi)/c%cp_ice
         Ts = (c%Lfresh + col%zqsn(1)/c%rhos)/c%cp_ice
         dhi = c%cp_ice*Ti*0.25_wp/c%Lfresh
         dhs = c%cp_ice*Ts*0.2_wp/c%Lfresh
         call run(col, c)
         call check(error, col%ierr == 0_int32, "no error flag")
         if (allocated(error)) return
         call check(error, near(col%hin, 1.0_wp - dhi, 1.0e-13_wp, 0.0_wp), "warm ice layer thinned to 0 degC")
         if (allocated(error)) return
         call check(error, near(col%melts, dhs, 1.0e-13_wp, 0.0_wp) .and. near(col%hsn, 0.2_wp - dhs, 1.0e-13_wp, 0.0_wp), &
                    "warm snow melted back to 0 degC")
         if (allocated(error)) return
         call check(error, abs(energy_residual(col, c)) < 1.0e-3_wp, "energy closes")
         if (allocated(error)) return
         ! growth enthalpy without brine: -rhoi (Lfresh - cp_ice Tbot)
         call base_col(col, c, 4, 0.6_wp, 0.0_wp, -8.0_wp, -5.0_wp)
         col%l_brine = .false.
         col%hslyr = 0.0_wp
         col%fcondbotn = -10.0_wp
         col%Tbot = -1.8_wp
         qbot = -c%rhoi*(-c%cp_ice*col%Tbot + c%Lfresh)
         dzi = col%fcondbotn*DT/qbot
         call run(col, c)
         call check(error, col%ierr == 0_int32 .and. near(col%congel, dzi, 1.0e-14_wp, 0.0_wp), "fresh-ice growth rate")
         if (allocated(error)) return
         call check(error, abs(energy_residual(col, c)) < 1.0e-3_wp, "growth: energy closes")
         if (allocated(error)) return
      end do
   end subroutine test_fresh_option

   ! ------------------------------------------------------------------
   ! error bits
   ! ------------------------------------------------------------------

   subroutine test_error_bits(error)
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      type(column_t) :: col

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      ! Tbot at the freezing point of fresh water: the growth enthalpy divides by it
      call base_col(col, c, 4, 1.0_wp, 0.1_wp, -8.0_wp, -5.0_wp)
      col%Tbot = 0.0_wp
      call run(col, c)
      call check(error, col%ierr == ICE_ERR_TBOT, "Tbot = 0 sets ICE_ERR_TBOT only")
      if (allocated(error)) return
      call base_col(col, c, 4, 1.0_wp, 0.1_wp, -8.0_wp, -5.0_wp)
      col%Tbot = -0.5_wp*c%puny
      call run(col, c)
      call check(error, col%ierr == ICE_ERR_TBOT, "Tbot above -puny sets ICE_ERR_TBOT")
      if (allocated(error)) return
      call base_col(col, c, 4, 1.0_wp, 0.1_wp, -8.0_wp, -5.0_wp)
      col%l_brine = .false.
      col%Tbot = 0.0_wp
      call run(col, c)
      call check(error, col%ierr == 0_int32, "fresh-ice growth does not divide by Tbot")
      if (allocated(error)) return
      ! negative layer thickness
      call base_col(col, c, 4, 1.0_wp, 0.1_wp, -8.0_wp, -5.0_wp)
      col%hilyr = -1.0e-3_wp
      call run(col, c)
      call check(error, col%ierr == ICE_ERR_NEG_THICK, "negative layer thickness sets ICE_ERR_NEG_THICK")
      if (allocated(error)) return
      call base_col(col, c, 4, 1.0_wp, 0.1_wp, -8.0_wp, -5.0_wp)
      col%hslyr = -1.0e-3_wp
      call run(col, c)
      call check(error, col%ierr == ICE_ERR_NEG_THICK, "negative snow layer sets ICE_ERR_NEG_THICK")
      if (allocated(error)) return
      ! snow enthalpy that is not a snow state (the melt divides by it), dt <= 0, layer count
      call base_col(col, c, 4, 1.0_wp, 0.1_wp, -8.0_wp, -5.0_wp)
      col%zqsn(1) = 0.0_wp
      call run(col, c)
      call check(error, col%ierr == ICE_ERR_NONFINITE, "zero snow enthalpy is refused")
      if (allocated(error)) return
      call base_col(col, c, 4, 1.0_wp, 0.1_wp, -8.0_wp, -5.0_wp)
      col%dt = 0.0_wp
      call run(col, c)
      call check(error, col%ierr == ICE_ERR_NONFINITE, "dt = 0 is refused")
      if (allocated(error)) return
      call base_col(col, c, 4, 1.0_wp, 0.1_wp, -8.0_wp, -5.0_wp)
      col%nilyr = 9
      call run_unchecked_nilyr(col, c)
      call check(error, col%ierr == ICE_ERR_NONFINITE, "layer count above the bound is refused")
      if (allocated(error)) return
      ! the error word accumulates: a prior bit survives
      call base_col(col, c, 4, 1.0_wp, 0.1_wp, -8.0_wp, -5.0_wp)
      col%ierr = ICE_ERR_THERMO_CONSERV
      call run(col, c)
      call check(error, col%ierr == ICE_ERR_THERMO_CONSERV, "a clean column keeps the prior bits")
   end subroutine test_error_bits

   subroutine test_raw_snow_free(error)
      !! A bare-ice category with `qsno = 0` straight from the state runs clean and equals the run with
      !! `qsno = -rhos Lfresh`; real snow with `qsno = 0` is still refused.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      type(column_t) :: a, b
      integer :: iset, k

      do iset = ICE_CONSTANTS_ICEPACK, ICE_CONSTANTS_CESM
         c = ice_constants(iset)
         call base_col(a, c, 4, 1.0_wp, 0.0_wp, -8.0_wp, -5.0_wp)
         a%hslyr = 0.0_wp
         a%fsurfn = 20.0_wp
         a%fbot = -30.0_wp
         a%fsnow = 1.0e-5_wp
         b = a
         a%zqsn(1) = 0.0_wp
         b%zqsn(1) = -c%rhos*c%Lfresh
         call run(a, c)
         call run(b, c)
         call check(error, a%ierr == 0_int32 .and. b%ierr == 0_int32, "raw qsno = 0 runs clean")
         if (allocated(error)) return
         call check(error, same_bits(a%hin, b%hin) .and. same_bits(a%hsn, b%hsn) .and. same_bits(a%zqsn(1), b%zqsn(1)) .and. &
                    same_bits(a%fhocnn, b%fhocnn) .and. same_bits(a%efinal, b%efinal), "same result as qsno = -rhos Lfresh")
         if (allocated(error)) return
         do k = 1, 4
            call check(error, same_bits(a%zqin(k), b%zqin(k)), "same ice enthalpies")
            if (allocated(error)) return
         end do
         call base_col(a, c, 4, 1.0_wp, 0.1_wp, -8.0_wp, -5.0_wp)
         a%zqsn(1) = 0.0_wp
         call run(a, c)
         call check(error, a%ierr == ICE_ERR_NONFINITE, "real snow with qsno = 0 is refused")
         if (allocated(error)) return
      end do
   end subroutine test_raw_snow_free

   subroutine run_unchecked_nilyr(col, c)
      type(column_t), intent(inout) :: col
      type(ice_constants_t), intent(in) :: c
      real(wp) :: zq(9)
      zq = col%zqin(1)
      call ice_thickness_changes(9, 1, c, col%l_brine, col%dt, col%Tbot, col%fbot, col%flatn, col%fsurfn, col%fcondtopn, &
                                 col%fcondbotn, col%fsnow, col%zSin_bot, col%hin, col%hilyr, col%hsn, col%hslyr, zq, &
                                 col%zqsn, col%meltt, col%melts, col%meltb, col%congel, col%snoice, col%dsnow, &
                                 col%hsn_new, col%fhocnn, col%evapn, col%evapsn, col%evapin, col%efinal, col%ierr)
   end subroutine run_unchecked_nilyr

   subroutine test_nan_in(error)
      !! A non-finite input is flagged and NOTHING is written to the state.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      type(column_t) :: col, ref
      real(wp) :: nan
      integer :: which, k
      character(len=24) :: label

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      nan = ieee_value(1.0_wp, ieee_quiet_nan)
      do which = 1, 23
         call base_col(col, c, 4, 0.5_wp, 0.1_wp, -8.0_wp, -5.0_wp)
         col%fsnow = 1.0e-5_wp
         col%fbot = -10.0_wp
         write (label, '(a,i0)') "input ", which
         select case (which)
         case (1); col%dt = nan
         case (2); col%Tbot = nan
         case (3); col%fbot = nan
         case (4); col%flatn = nan
         case (5); col%fsurfn = nan
         case (6); col%fcondtopn = nan
         case (7); col%fcondbotn = nan
         case (8); col%fsnow = nan
         case (9); col%zSin_bot = nan
         case (10); col%hin = nan
         case (11); col%hsn = nan
         case (12); col%hilyr = nan
         case (13); col%hslyr = nan
         case (14); col%zqin(1) = nan
         case (15); col%zqin(4) = nan
         case (16); col%zqsn(1) = nan
         case (17); col%meltt = nan
         case (18); col%melts = nan
         case (19); col%meltb = nan
         case (20); col%congel = nan
         case (21); col%snoice = nan
         case (22); col%dsnow = nan
         case (23); col%dt = ieee_value(1.0_wp, ieee_positive_inf)
         end select
         ref = col
         call run(col, c)
         call check(error, iand(col%ierr, ICE_ERR_NONFINITE) /= 0_int32, trim(label)//": ICE_ERR_NONFINITE flagged")
         if (allocated(error)) return
         call check(error, same_bits(col%hin, ref%hin) .and. same_bits(col%hsn, ref%hsn) .and. same_bits(col%hilyr, ref%hilyr) &
                    .and. same_bits(col%hslyr, ref%hslyr) .and. same_bits(col%zqsn(1), ref%zqsn(1)) .and. &
                    same_bits(col%meltt, ref%meltt) .and. same_bits(col%melts, ref%melts) .and. same_bits(col%meltb, ref%meltb) &
                    .and. same_bits(col%congel, ref%congel) .and. same_bits(col%snoice, ref%snoice) .and. &
                    same_bits(col%dsnow, ref%dsnow), trim(label)//": state untouched")
         if (allocated(error)) return
         do k = 1, 4
            call check(error, same_bits(col%zqin(k), ref%zqin(k)), trim(label)//": layer enthalpies untouched")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_nan_in

   function same_bits(a, b) result(ok)
      real(wp), intent(in) :: a, b
      logical :: ok
      ok = transfer(a, 0_int64) == transfer(b, 0_int64)
   end function same_bits

   subroutine test_check_threshold(error)
      !! `|efinal - einit - einp|/dt > 1.1 ferrmax` sets the bit, exactly at the threshold it does not.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      real(wp) :: dt, ferr, einit, einp, efinal, nan
      integer(int32) :: ierr
      integer :: s
      real(wp), parameter :: FS(7) = [0.0_wp, 0.5_wp, 1.0_wp, 1.0999_wp, 1.1001_wp, 5.0_wp, -2.0_wp]
      logical, parameter :: EXPECT(7) = [.false., .false., .false., .false., .true., .true., .true.]

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      dt = 900.0_wp
      einit = -2.0e8_wp
      einp = (30.0_wp - (-5.0_wp) + 7.0_wp - (-20.0_wp) - 1.0e-4_wp*c%Lfresh)*dt
      do s = 1, 7
         efinal = einit + einp + FS(s)*ICE_FERRMAX*dt
         ierr = 0_int32
         call ice_conservation_check_vthermo(c, dt, 30.0_wp, -5.0_wp, -20.0_wp, 7.0_wp, 1.0e-4_wp, einit, efinal, &
                                             ICE_FERRMAX, ferr, ierr)
         call check(error, (iand(ierr, ICE_ERR_THERMO_CONSERV) /= 0_int32) .eqv. EXPECT(s), "threshold at 1.1 ferrmax")
         if (allocated(error)) return
         call check(error, near(ferr, abs(FS(s))*ICE_FERRMAX, 1.0e-4_wp, 1.0e-9_wp), "ferr value [W/m2]")
         if (allocated(error)) return
      end do
      nan = ieee_value(1.0_wp, ieee_quiet_nan)
      ierr = 0_int32
      call ice_conservation_check_vthermo(c, dt, nan, -5.0_wp, -20.0_wp, 7.0_wp, 1.0e-4_wp, einit, einit, ICE_FERRMAX, ferr, ierr)
      call check(error, ierr == ICE_ERR_NONFINITE, "NaN input is a non-finite error, not a conservation error")
   end subroutine test_check_threshold

   ! ------------------------------------------------------------------
   ! random columns: exact closure (A9)
   ! ------------------------------------------------------------------

   subroutine draw_column(st, c, nilyr, col, fswint, einit, dt_fixed)
      !! A random but physical column: ice 0.02-3 m (thin-biased), snow 0-0.5 m (30 % bare), layer
      !! temperatures below their melting points, fluxes whose energies span the column's own.
      integer(int64), intent(inout) :: st
      type(ice_constants_t), intent(in) :: c
      integer, intent(in) :: nilyr
      type(column_t), intent(out) :: col
      real(wp), intent(out) :: fswint, einit
      real(wp), intent(in), optional :: dt_fixed
      integer :: k
      real(wp) :: u, S, tm, T, etot, e_i, e_s, mag, ebot, etop, esub

      col%nilyr = nilyr
      col%dt = merge(900.0_wp, 3600.0_wp, uniform(st) < 0.7_wp)
      if (present(dt_fixed)) col%dt = dt_fixed
      u = uniform(st)
      col%hin = 0.02_wp + 3.0_wp*u*u
      col%hilyr = col%hin/real(nilyr, wp)
      col%hsn = 0.0_wp
      if (uniform(st) > 0.3_wp) col%hsn = 0.5_wp*uniform(st)
      col%hslyr = col%hsn
      do k = 1, nilyr
         S = 0.1_wp + 5.9_wp*uniform(st)
         tm = -c%depressT*S
         T = tm - 0.01_wp - 30.0_wp*uniform(st)**2
         col%zqin(k) = ice_q_from_T(T, S, c, .true.)
         if (k == nilyr) col%zSin_bot = S
      end do
      col%zqsn(1) = ice_qsno_from_T(-0.05_wp - 30.0_wp*uniform(st), c)
      col%Tbot = -1.9_wp + 0.8_wp*uniform(st)
      e_i = 0.0_wp
      do k = 1, nilyr
         e_i = e_i - col%hilyr*col%zqin(k)
      end do
      e_s = -col%hsn*col%zqsn(1)
      etot = e_i + e_s
      col%fcondtopn = -15.0_wp + 30.0_wp*uniform(st)
      col%fcondbotn = -25.0_wp + 30.0_wp*uniform(st)
      fswint = 30.0_wp*uniform(st)
      ! top melt energy: none, partial, or beyond the whole column
      etop = 0.0_wp
      u = uniform(st)
      if (u > 0.45_wp) then
         mag = 10.0_wp**(3.0_wp + log10(max(3.0_wp*etot, 2.0e3_wp)/1.0e3_wp)*uniform(st))
         etop = mag
      end if
      ! base: melt (including through the column) or growth
      ebot = 0.0_wp
      u = uniform(st)
      if (u < 0.4_wp) then
         ebot = 10.0_wp**(3.0_wp + log10(max(3.0_wp*etot, 2.0e3_wp)/1.0e3_wp)*uniform(st))
      else if (u < 0.85_wp) then
         ebot = -10.0_wp**(3.0_wp + 3.0_wp*uniform(st))
      end if
      ! latent: sublimation, condensation or none, up to the column's own latent energy
      esub = 0.0_wp
      u = uniform(st)
      if (u < 0.3_wp) then
         esub = 10.0_wp**(3.0_wp + 4.0_wp*uniform(st))
      else if (u < 0.5_wp) then
         esub = -10.0_wp**(3.0_wp + 3.0_wp*uniform(st))
      end if
      col%fsurfn = col%fcondtopn + etop/col%dt
      col%fbot = col%fcondbotn - ebot/col%dt
      col%flatn = -esub/col%dt
      col%fsnow = 0.0_wp
      if (uniform(st) < 0.4_wp) col%fsnow = 2.0e-3_wp*uniform(st)
      einit = energy(col) - (col%fcondtopn - col%fcondbotn + fswint)*col%dt   ! energy before the temperature solve
   end subroutine draw_column

   subroutine test_closure_random(error)
      !! Energy, thickness, water and salt close per column to round-off, and the check never trips.
      type(error_type), allocatable, intent(out) :: error
      type(ice_constants_t) :: c
      type(column_t) :: col
      type(ice_ledger_t) :: led
      integer(int64) :: st
      integer :: iset, is, ni, i, nfail, nmelt_out, ncong, nsnoice, nsub
      real(wp) :: fswint, einit, fresh, fsalt, ferr, scale, dhi_exp, dhs_exp, dm, ds
      integer(int32) :: ierr
      integer, parameter :: NCOL = 100000

      do iset = ICE_CONSTANTS_ICEPACK, ICE_CONSTANTS_CESM
         c = ice_constants(iset)
         led = ice_ledger_t()
         nfail = 0
         nmelt_out = 0
         ncong = 0
         nsnoice = 0
         nsub = 0
         do is = 1, size(NIL_SET)
            ni = NIL_SET(is)
            st = int(7919*is + 104729*iset + 17, int64)
            do i = 1, NCOL
               call draw_column(st, c, ni, col, fswint, einit)
               call run(col, c)
               if (col%ierr /= 0_int32) then
                  nfail = nfail + 1
                  cycle
               end if
               if (col%hin == 0.0_wp) nmelt_out = nmelt_out + 1
               if (col%congel > 0.0_wp) ncong = ncong + 1
               if (col%snoice > 0.0_wp) nsnoice = nsnoice + 1
               if (col%evapn /= 0.0_wp) nsub = nsub + 1
               scale = abs(col%efinal) + abs(col%einter) + abs(col%fhocnn)*col%dt + 1.0_wp
               call led%add("energy", energy_residual(col, c), scale)
               ierr = 0_int32
               call ice_conservation_check_vthermo(c, col%dt, col%fsurfn, col%flatn, col%fhocnn, fswint, col%fsnow, einit, &
                                                   col%efinal, ICE_FERRMAX, ferr, ierr)
               if (ierr /= 0_int32) nfail = nfail + 1
               call led%add("check_ferr_over_ferrmax", ferr, ICE_FERRMAX)
               ! thickness ledger: ice and snow changes against the accumulators
               dhi_exp = col%congel - col%meltt - col%meltb + col%snoice + col%evapin*col%dt/c%rhoi
               dhs_exp = col%hsn_new - col%melts + col%evapsn*col%dt/c%rhos - col%snoice*c%rhoi/c%rhos
               call led%add("ice_thickness", (col%hin - col%hin0) - dhi_exp, col%hin0 + col%hin + 1.0e-3_wp)
               call led%add("snow_thickness", (col%hsn - col%hsn0) - dhs_exp, col%hsn0 + col%hsn + col%hsn_new + 1.0e-3_wp)
               ! water: dM = (fsnow + evap - fresh) dt with the booked fresh water; salt: dS_ice = -fsalt dt
               call hand_off(col, c, fresh, fsalt)
               dm = c%rhoi*(col%hin - col%hin0) + c%rhos*(col%hsn - col%hsn0)
               call led%add("water", dm - (col%fsnow + col%evapn - fresh)*col%dt, c%rhoi*abs(col%hin - col%hin0) &
                            + c%rhos*abs(col%hsn - col%hsn0) + col%fsnow*col%dt + 1.0e-3_wp)
               ds = c%rhoi*(col%hin - col%hin0)*SREF*0.001_wp
               call led%add("salt", ds + fsalt*col%dt, abs(ds) + 1.0e-6_wp)
            end do
         end do
         write (output_unit, '(a,i0,a,6(1x,es10.3))') "closure set ", iset, " worst:", led%worst_of("energy"), &
            led%worst_of("ice_thickness"), led%worst_of("snow_thickness"), led%worst_of("water"), led%worst_of("salt"), &
            led%worst_of("check_ferr_over_ferrmax")
         call check(error, nfail == 0, "no error flags and no conservation trips over 5e5 columns")
         if (allocated(error)) return
         call check(error, nmelt_out > 100 .and. ncong > 1000 .and. nsnoice > 100 .and. nsub > 1000, &
                    "the random set exercises melt-out, growth, snow-ice and sublimation")
         if (allocated(error)) return
         call check(error, led%worst_of("energy") <= 1.0e-14_wp, "energy closes to round-off (measured worst 8.3e-16)")
         if (allocated(error)) return
         call check(error, led%worst_of("ice_thickness") <= 1.0e-14_wp, "ice thickness closes")
         if (allocated(error)) return
         call check(error, led%worst_of("snow_thickness") <= 1.0e-14_wp, "snow thickness closes")
         if (allocated(error)) return
         call check(error, led%worst_of("water") <= 4.0e-13_wp, "water closes (measured worst 3.5e-14)")
         if (allocated(error)) return
         call check(error, led%worst_of("salt") <= 2.0e-15_wp, "salt closes")
         if (allocated(error)) return
         call check(error, led%worst_of("check_ferr_over_ferrmax") <= 2.0e-5_wp, "the upstream energy check passes with 5e4 to spare (measured worst ferr 1.3e-9 W/m2)")
         if (allocated(error)) return
         call check(error,.not. led%overflow, "ledger did not overflow")
         if (allocated(error)) return
      end do
   end subroutine test_closure_random

   ! ------------------------------------------------------------------
   ! goldens from the upstream routine
   ! ------------------------------------------------------------------

   subroutine load_goldens(set, GT, GC, NT, NC)
      integer, intent(in) :: set
      integer, intent(in) :: NT, NC
      real(wp), intent(out) :: GT(55, NT), GC(12, NC)
#include "test_cice_thickness_changes_golden.inc"
   end subroutine load_goldens

   subroutine test_golden(error)
      !! 66 columns per constant set (top melt through several layers, bottom melt, growth with and without the
      !! enthalpy cap, melt-out, snow melted from below, snow-ice, sublimation, condensation, snowfall, the
      !! fresh-ice option; nilyr 1, 2, 4, 7, 8) against the real `thickness_changes`, every output; and the
      !! conservation-check trip decision against the real `conservation_check_vthermo`.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NT = 66, NC = 20
      real(wp) :: GT(55, NT), GC(12, NC)
      type(ice_constants_t) :: c
      type(column_t) :: col
      integer :: iset, r, k, ni, io
      real(wp) :: got(25), want(25), ferr, e0, worst
      integer(int32) :: ierr
      character(len=64) :: names(25)

      names = [character(len=64) :: "hin", "hsn", "hilyr", "hslyr", "zqin1", "zqin2", "zqin3", "zqin4", "zqin5", "zqin6", &
               "zqin7", "zqin8", "zqsn", "hsn_new", "fhocnn", "evapn", "evapsn", "evapin", "efinal", &
               "meltt", "melts", "meltb", "congel", "snoice", "dsnow"]
      worst = 0.0_wp
      do iset = ICE_CONSTANTS_ICEPACK, ICE_CONSTANTS_CESM
         c = ice_constants(iset)
         call load_goldens(iset, GT, GC, NT, NC)
         do r = 1, NT
            ni = nint(GT(1, r))
            col = column_t()
            col%nilyr = ni
            col%l_brine = GT(2, r) > 0.5_wp
            col%dt = GT(3, r); col%Tbot = GT(4, r); col%fbot = GT(5, r); col%flatn = GT(6, r); col%fsurfn = GT(7, r)
            col%fcondtopn = GT(8, r); col%fcondbotn = GT(9, r); col%fsnow = GT(10, r); col%zSin_bot = GT(11, r)
            col%hin = GT(12, r); col%hsn = GT(13, r); col%hilyr = GT(14, r); col%hslyr = GT(15, r)
            col%zqin = GT(16:23, r); col%zqsn(1) = GT(24, r)
            col%meltt = GT(25, r); col%melts = GT(26, r); col%meltb = GT(27, r); col%congel = GT(28, r)
            col%snoice = GT(29, r); col%dsnow = GT(30, r)
            call run(col, c)
            call check(error, col%ierr == 0_int32, "golden column raises no flag")
            if (allocated(error)) return
            got = [col%hin, col%hsn, col%hilyr, col%hslyr, col%zqin(1:8), col%zqsn(1), col%hsn_new, col%fhocnn, col%evapn, &
                   col%evapsn, col%evapin, col%efinal, col%meltt, col%melts, col%meltb, col%congel, col%snoice, col%dsnow]
            want = GT(31:55, r)
            do k = 1, 25
               if (k >= 5 .and. k <= 12 .and. k - 4 > ni) cycle   ! unused layers
               if (.not. ieee_is_finite(got(k))) then
                  call test_failed(error, "non-finite "//trim(names(k))//" in golden row")
                  return
               end if
               worst = max(worst, abs(got(k) - want(k))/max(abs(want(k)), golden_atol(k, want)*1.0e12_wp, 1.0e-300_wp))
               if (.not. near(got(k), want(k), GRTOL, golden_atol(k, want))) then
                  call golden_fail(error, iset, r, trim(names(k)), got(k), want(k))
                  return
               end if
            end do
         end do
         do r = 1, NC
            ierr = 0_int32
            e0 = GC(7, r)
            call ice_conservation_check_vthermo(c, GC(1, r), GC(2, r), GC(3, r), GC(4, r), GC(5, r), GC(6, r), e0, GC(8, r), &
                                                ICE_FERRMAX, ferr, ierr)
            io = merge(1, 0, iand(ierr, ICE_ERR_THERMO_CONSERV) /= 0_int32)
            if (io /= nint(GC(12, r))) then
               call test_failed(error, "conservation-check decision differs from upstream")
               return
            end if
         end do
      end do
      write (output_unit, '(a,es10.3)') "golden worst difference (relative to the output, or to its group scale where the output is ~0), 25 outputs x 66 columns x 2 sets: ", worst
   end subroutine test_golden

   function golden_atol(k, want) result(a)
      !! Absolute floor of the GPU-build golden comparison: 1e-12 of the scale of the output's group (a
      !! residual of the leftover-energy bookkeeping, not a relative error); 0 on CPU builds.
      integer, intent(in) :: k
      real(wp), intent(in) :: want(25)
      real(wp) :: a
      a = GATOL
      if (GRTOL == 0.0_wp) return
      select case (k)
      case (1:4, 14, 20:25)   ! thicknesses and their accumulators [m]
         a = 1.0e-12_wp*max(1.0_wp, want(1))
      case (16:18)            ! evaporation fluxes [kg/m2/s]
         a = 1.0e-12_wp*max(abs(want(16)), abs(want(17)), abs(want(18)))
      case (15)               ! fhocnn [W/m2]
         a = 1.0e-12_wp*max(abs(want(15)), 1.0e2_wp)
      case default
         a = GATOL
      end select
   end function golden_atol

   subroutine golden_fail(error, iset, row, label, got, want)
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: iset, row
      character(len=*), intent(in) :: label
      real(wp), intent(in) :: got, want
      character(len=200) :: msg
      write (msg, '(a,i0,a,i0,a,a,a,es24.16,a,es24.16)') "golden set ", iset, " row ", row, " ", label, ": got ", got, &
         " want ", want
      call test_failed(error, trim(msg))
   end subroutine golden_fail

   ! ------------------------------------------------------------------
   ! cell entry, device sweep
   ! ------------------------------------------------------------------

   subroutine fill_state(nx, ny, ncat, nilyr, c, seed, Tbot, fbot, fsnow, flatn, fsurfn, fcondtopn, fcondbotn, sicen, &
                         hin, hsn, hilyr, hslyr, qicen, qsnon, fswint, einit)
      !! Random state where every cell differs: a draw per (i, j, n); a few cells are given non-finite fluxes,
      !! `Tbot = 0` or a negative layer so that the error words differ between cells too.
      integer, intent(in) :: nx, ny, ncat, nilyr
      type(ice_constants_t), intent(in) :: c
      integer(int64), intent(in) :: seed
      real(wp), intent(out) :: Tbot(nx, ny), fbot(nx, ny), fsnow(nx, ny)
      real(wp), intent(out) :: flatn(nx, ny, ncat), fsurfn(nx, ny, ncat), fcondtopn(nx, ny, ncat), fcondbotn(nx, ny, ncat)
      real(wp), intent(out) :: sicen(nx, ny, ncat, nilyr), hin(nx, ny, ncat), hsn(nx, ny, ncat), hilyr(nx, ny, ncat)
      real(wp), intent(out) :: hslyr(nx, ny, ncat), qicen(nx, ny, ncat, nilyr), qsnon(nx, ny, ncat, 1)
      real(wp), intent(out) :: fswint(nx, ny, ncat), einit(nx, ny, ncat)
      integer(int64) :: st
      integer :: i, j, n, k
      type(column_t) :: col
      real(wp) :: fsw, ein, u

      st = seed
      do n = 1, ncat
         do j = 1, ny
            do i = 1, nx
               call draw_column(st, c, nilyr, col, fsw, ein, DT)
               Tbot(i, j) = col%Tbot
               fbot(i, j) = col%fbot
               fsnow(i, j) = col%fsnow
               flatn(i, j, n) = col%flatn
               fsurfn(i, j, n) = col%fsurfn
               fcondtopn(i, j, n) = col%fcondtopn
               fcondbotn(i, j, n) = col%fcondbotn
               hin(i, j, n) = col%hin
               hsn(i, j, n) = col%hsn
               hilyr(i, j, n) = col%hilyr
               hslyr(i, j, n) = col%hslyr
               do k = 1, nilyr
                  qicen(i, j, n, k) = col%zqin(k)
                  sicen(i, j, n, k) = col%zSin_bot*real(k, wp)/real(nilyr, wp)
               end do
               sicen(i, j, n, nilyr) = col%zSin_bot
               qsnon(i, j, n, 1) = col%zqsn(1)
               fswint(i, j, n) = fsw
               einit(i, j, n) = ein
               u = uniform(st)
               if (u < 0.002_wp) then
                  flatn(i, j, n) = ieee_value(1.0_wp, ieee_quiet_nan)
               else if (u < 0.004_wp) then
                  Tbot(i, j) = 0.0_wp
               else if (u < 0.006_wp) then
                  hilyr(i, j, n) = -1.0e-3_wp
               end if
            end do
         end do
      end do
      ! the fluxes of a cell are shared by its categories: Tbot, fbot, fsnow come from the last draw, as in the model
   end subroutine fill_state

   subroutine test_cell_vs_column(error)
      !! The cell entry (gather, column routine, scatter) equals the column routine on the same data, on a
      !! small grid where cells and categories differ, for all outputs.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 7, NY = 5, NCAT = 3, NIL = 4
      type(ice_constants_t) :: c
      real(wp) :: Tbot(NX, NY), fbot(NX, NY), fsnow(NX, NY)
      real(wp) :: flatn(NX, NY, NCAT), fsurfn(NX, NY, NCAT), fcondtopn(NX, NY, NCAT), fcondbotn(NX, NY, NCAT)
      real(wp) :: sicen(NX, NY, NCAT, NIL), hin(NX, NY, NCAT), hsn(NX, NY, NCAT), hilyr(NX, NY, NCAT)
      real(wp) :: hslyr(NX, NY, NCAT), qicen(NX, NY, NCAT, NIL), qsnon(NX, NY, NCAT, 1), fswint(NX, NY, NCAT)
      real(wp) :: einit(NX, NY, NCAT)
      real(wp) :: acc(NX, NY, NCAT, 6), outv(NX, NY, NCAT, 6)
      integer(int32) :: ierr(NX, NY, NCAT)
      real(wp) :: hin0(NX, NY, NCAT), hsn0(NX, NY, NCAT), hil0(NX, NY, NCAT), hsl0(NX, NY, NCAT)
      real(wp) :: qi0(NX, NY, NCAT, NIL), qs0(NX, NY, NCAT, 1)
      type(column_t) :: col
      integer :: i, j, n, k

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      call fill_state(NX, NY, NCAT, NIL, c, 99_int64, Tbot, fbot, fsnow, flatn, fsurfn, fcondtopn, fcondbotn, sicen, hin, hsn, &
                      hilyr, hslyr, qicen, qsnon, fswint, einit)
      hin0 = hin; hsn0 = hsn; hil0 = hilyr; hsl0 = hslyr; qi0 = qicen; qs0 = qsnon
      acc = 0.0_wp
      ierr = 0_int32
      do n = 1, NCAT
         do j = 1, NY
            do i = 1, NX
               call ice_thickness_changes_cell(i, j, n, NX, NY, NCAT, NIL, 1, c, .true., DT, Tbot, fbot, fsnow, flatn, fsurfn, &
                                               fcondtopn, fcondbotn, sicen, hin, hsn, hilyr, hslyr, qicen, qsnon, &
                                               acc(:, :, :, 1), acc(:, :, :, 2), acc(:, :, :, 3), acc(:, :, :, 4), &
                                               acc(:, :, :, 5), acc(:, :, :, 6), outv(:, :, :, 1), outv(:, :, :, 2), &
                                               outv(:, :, :, 3), outv(:, :, :, 4), outv(:, :, :, 5), outv(:, :, :, 6), ierr)
            end do
         end do
      end do
      do n = 1, NCAT
         do j = 1, NY
            do i = 1, NX
               col = column_t()
               col%nilyr = NIL
               col%Tbot = Tbot(i, j); col%fbot = fbot(i, j); col%fsnow = fsnow(i, j)
               col%flatn = flatn(i, j, n); col%fsurfn = fsurfn(i, j, n); col%fcondtopn = fcondtopn(i, j, n)
               col%fcondbotn = fcondbotn(i, j, n); col%zSin_bot = sicen(i, j, n, NIL)
               col%hin = hin0(i, j, n); col%hsn = hsn0(i, j, n); col%hilyr = hil0(i, j, n); col%hslyr = hsl0(i, j, n)
               col%zqin(1:NIL) = qi0(i, j, n, :); col%zqsn(1) = qs0(i, j, n, 1)
               call run(col, c)
               call check(error, col%ierr == ierr(i, j, n), "error word equal")
               if (allocated(error)) return
               call check(error, same_bits(col%hin, hin(i, j, n)) .and. same_bits(col%hsn, hsn(i, j, n)) .and. &
                          same_bits(col%hilyr, hilyr(i, j, n)) .and. same_bits(col%hslyr, hslyr(i, j, n)) .and. &
                          same_bits(col%zqsn(1), qsnon(i, j, n, 1)), "state equal")
               if (allocated(error)) return
               do k = 1, NIL
                  call check(error, same_bits(col%zqin(k), qicen(i, j, n, k)), "enthalpies equal")
                  if (allocated(error)) return
               end do
               call check(error, same_bits(col%meltt, acc(i, j, n, 1)) .and. same_bits(col%melts, acc(i, j, n, 2)) .and. &
                          same_bits(col%meltb, acc(i, j, n, 3)) .and. same_bits(col%congel, acc(i, j, n, 4)) .and. &
                          same_bits(col%snoice, acc(i, j, n, 5)) .and. same_bits(col%dsnow, acc(i, j, n, 6)), "accumulators equal")
               if (allocated(error)) return
               call check(error, same_bits(col%hsn_new, outv(i, j, n, 1)) .and. same_bits(col%fhocnn, outv(i, j, n, 2)) .and. &
                          same_bits(col%evapn, outv(i, j, n, 3)) .and. same_bits(col%evapsn, outv(i, j, n, 4)) .and. &
                          same_bits(col%evapin, outv(i, j, n, 5)) .and. same_bits(col%efinal, outv(i, j, n, 6)), "outputs equal")
               if (allocated(error)) return
            end do
         end do
      end do
      call check(error, einit(1, 1, 1) == einit(1, 1, 1) .and. fswint(1, 1, 1) >= 0.0_wp, "check inputs generated")
   end subroutine test_cell_vs_column

   subroutine sweep(nx, ny, ncat, nilyr, c, lb, Tbot, fbot, fsnow, flatn, fsurfn, fcondtopn, fcondbotn, sicen, hin, hsn, hilyr, &
                    hslyr, qicen, qsnon, meltt, melts, meltb, congel, snoice, dsnow, hsn_new, fhocnn, evapn, evapsn, evapin, &
                    efinal, fswint, einit, ferr, ierr)
      !! The category-parallel driver shape: ONE `do concurrent (n, j, i)` calling the two cell kernels (bare
      !! index arguments, the form the by-value dummies exist for).  Every array is device-resident on the GPU
      !! build (`mem:separate`), mapped by the caller.
      integer, intent(in) :: nx, ny, ncat, nilyr
      type(ice_constants_t), intent(in) :: c
      logical, intent(in) :: lb
      real(wp), intent(in) :: Tbot(nx, ny), fbot(nx, ny), fsnow(nx, ny)
      real(wp), intent(in) :: flatn(nx, ny, ncat), fsurfn(nx, ny, ncat), fcondtopn(nx, ny, ncat), fcondbotn(nx, ny, ncat)
      real(wp), intent(in) :: sicen(nx, ny, ncat, nilyr), fswint(nx, ny, ncat), einit(nx, ny, ncat)
      real(wp), intent(inout) :: hin(nx, ny, ncat), hsn(nx, ny, ncat), hilyr(nx, ny, ncat), hslyr(nx, ny, ncat)
      real(wp), intent(inout) :: qicen(nx, ny, ncat, nilyr), qsnon(nx, ny, ncat, 1)
      real(wp), intent(inout) :: meltt(nx, ny, ncat), melts(nx, ny, ncat), meltb(nx, ny, ncat)
      real(wp), intent(inout) :: congel(nx, ny, ncat), snoice(nx, ny, ncat), dsnow(nx, ny, ncat)
      real(wp), intent(inout) :: hsn_new(nx, ny, ncat), fhocnn(nx, ny, ncat), evapn(nx, ny, ncat)
      real(wp), intent(inout) :: evapsn(nx, ny, ncat), evapin(nx, ny, ncat), efinal(nx, ny, ncat), ferr(nx, ny, ncat)
      integer(int32), intent(inout) :: ierr(nx, ny, ncat)
      integer :: i, j, n

      do concurrent(n=1:ncat, j=1:ny, i=1:nx)
         call ice_thickness_changes_cell(i, j, n, nx, ny, ncat, nilyr, 1, c, lb, DT, Tbot, fbot, fsnow, flatn, fsurfn, &
                                         fcondtopn, fcondbotn, sicen, hin, hsn, hilyr, hslyr, qicen, qsnon, meltt, melts, &
                                         meltb, congel, snoice, dsnow, hsn_new, fhocnn, evapn, evapsn, evapin, efinal, ierr)
         call ice_conservation_check_vthermo_cell(i, j, n, nx, ny, ncat, c, DT, fsurfn, flatn, fhocnn, fswint, fsnow, einit, &
                                                  efinal, ICE_FERRMAX, ferr, ierr)
      end do
   end subroutine sweep

   subroutine test_device(error)
      type(error_type), allocatable, intent(out) :: error
      call run_device_case(error, .true.)
      if (allocated(error)) return
      call run_device_case(error, .false.)
   end subroutine test_device

   subroutine run_device_case(error, lb)
      type(error_type), allocatable, intent(out) :: error
      logical, intent(in) :: lb
      integer, parameter :: NX = 512, NY = 512, NCAT = 5, NIL = 7
      type(ice_constants_t) :: c
      real(wp), allocatable :: Tbot(:, :), fbot(:, :), fsnow(:, :)
      real(wp), allocatable :: flatn(:, :, :), fsurfn(:, :, :), fcondtopn(:, :, :), fcondbotn(:, :, :), fswint(:, :, :)
      real(wp), allocatable :: einit(:, :, :), sicen(:, :, :, :)
      real(wp), allocatable :: hin(:, :, :), hsn(:, :, :), hilyr(:, :, :), hslyr(:, :, :), qicen(:, :, :, :), qsnon(:, :, :, :)
      real(wp), allocatable :: hin_h(:, :, :), hsn_h(:, :, :), hilyr_h(:, :, :), hslyr_h(:, :, :)
      real(wp), allocatable :: qicen_h(:, :, :, :), qsnon_h(:, :, :, :)
      real(wp), allocatable :: o_d(:, :, :, :), o_h(:, :, :, :)   ! 13 per-category outputs, last index
      integer(int32), allocatable :: ierr_d(:, :, :), ierr_h(:, :, :)
      integer :: i, j, n, k, nbad_err
      integer(int64) :: t0, t1, rate
      real(wp) :: tdev, thost
      character(len=40) :: names(13)

      c = ice_constants(ICE_CONSTANTS_ICEPACK)
      allocate (Tbot(NX, NY), fbot(NX, NY), fsnow(NX, NY), flatn(NX, NY, NCAT), fsurfn(NX, NY, NCAT), fcondtopn(NX, NY, NCAT))
      allocate (fcondbotn(NX, NY, NCAT), fswint(NX, NY, NCAT), einit(NX, NY, NCAT), sicen(NX, NY, NCAT, NIL))
      allocate (hin(NX, NY, NCAT), hsn(NX, NY, NCAT), hilyr(NX, NY, NCAT), hslyr(NX, NY, NCAT))
      allocate (qicen(NX, NY, NCAT, NIL), qsnon(NX, NY, NCAT, 1))
      allocate (hin_h(NX, NY, NCAT), hsn_h(NX, NY, NCAT), hilyr_h(NX, NY, NCAT), hslyr_h(NX, NY, NCAT))
      allocate (qicen_h(NX, NY, NCAT, NIL), qsnon_h(NX, NY, NCAT, 1))
      allocate (o_d(NX, NY, NCAT, 13), o_h(NX, NY, NCAT, 13), ierr_d(NX, NY, NCAT), ierr_h(NX, NY, NCAT))
      names = [character(len=40) :: "meltt", "melts", "meltb", "congel", "snoice", "dsnow", "hsn_new", "fhocnn", "evapn", &
               "evapsn", "evapin", "efinal", "ferr"]
      call fill_state(NX, NY, NCAT, NIL, c, 424242_int64, Tbot, fbot, fsnow, flatn, fsurfn, fcondtopn, fcondbotn, sicen, &
                      hin, hsn, hilyr, hslyr, qicen, qsnon, fswint, einit)
      ! accumulators start non-zero and different per cell
      do n = 1, NCAT
         do j = 1, NY
            do i = 1, NX
               o_d(i, j, n, 1:6) = 1.0e-4_wp*real(1 + mod(i + 3*j + 5*n, 17), wp)
            end do
         end do
      end do
      o_d(:, :, :, 7:13) = -1.0_wp   ! sentinel: every one of these must be overwritten
      o_h = o_d
      ierr_d = 0_int32
      ierr_h = 0_int32
      hin_h = hin; hsn_h = hsn; hilyr_h = hilyr; hslyr_h = hslyr; qicen_h = qicen; qsnon_h = qsnon

      ! host sweep (plain loops over the same cell kernels)
      call system_clock(t0, rate)
      call host_sweep()
      call system_clock(t1)
      thost = real(t1 - t0, wp)/real(rate, wp)

      ! device sweep
      !$acc enter data copyin(c, Tbot, fbot, fsnow, flatn, fsurfn, fcondtopn, fcondbotn, sicen, fswint, einit, hin, hsn, hilyr, &
      !$acc&  hslyr, qicen, qsnon, o_d, ierr_d)
      call system_clock(t0)
      call sweep(NX, NY, NCAT, NIL, c, lb, Tbot, fbot, fsnow, flatn, fsurfn, fcondtopn, fcondbotn, sicen, hin, hsn, hilyr, hslyr, &
                 qicen, qsnon, o_d(:, :, :, 1), o_d(:, :, :, 2), o_d(:, :, :, 3), o_d(:, :, :, 4), o_d(:, :, :, 5), &
                 o_d(:, :, :, 6), o_d(:, :, :, 7), o_d(:, :, :, 8), o_d(:, :, :, 9), o_d(:, :, :, 10), o_d(:, :, :, 11), &
                 o_d(:, :, :, 12), fswint, einit, o_d(:, :, :, 13), ierr_d)
      !$acc wait
      call system_clock(t1)
      tdev = real(t1 - t0, wp)/real(rate, wp)
      !$acc update self(hin, hsn, hilyr, hslyr, qicen, qsnon, o_d, ierr_d)
      !$acc exit data delete(c, Tbot, fbot, fsnow, flatn, fsurfn, fcondtopn, fcondbotn, sicen, fswint, einit, hin, hsn, hilyr, &
      !$acc&  hslyr, qicen, qsnon, o_d, ierr_d)
      write (output_unit, '(a,i0,a,i0,a,i0,a,i0,a,f8.3,a,es10.3,a,f8.3,a,es10.3,a)') "thickness_changes sweep ", NX, "x", NY, "x", &
         NCAT, "x", NIL, ": sweep ", tdev, " s (", real(NX*NY*NCAT, wp)/max(tdev, 1.0e-9_wp), " cell-cat/s; host loop ", thost, &
         " s, ", real(NX*NY*NCAT, wp)/max(thost, 1.0e-9_wp), " cell-cat/s)"

      nbad_err = count(ierr_h /= 0_int32)
      write (output_unit, '(a,i0,a,i0,a,i0,a,i0)') "error words: nonzero ", nbad_err, ", nonfinite ", &
         count(iand(ierr_h, ICE_ERR_NONFINITE) /= 0_int32), ", Tbot ", count(iand(ierr_h, ICE_ERR_TBOT) /= 0_int32), &
         ", conservation ", count(iand(ierr_h, ICE_ERR_THERMO_CONSERV) /= 0_int32)
      ! (with l_brine off the brine-form enthalpies of the random state are not fresh-ice states, so the energy
      ! check trips on many columns; the sweep compares device and host, not physics)
      call check(error, nbad_err > 100 .and. (nbad_err < NX*NY*NCAT/50 .or. .not. lb), "the state exercises the error words")
      if (allocated(error)) return
      call check(error, all(ierr_d == ierr_h), "every error word equals the host's")
      if (allocated(error)) return
      call compare3(error, "hin", hin, hin_h); if (allocated(error)) return
      call compare3(error, "hsn", hsn, hsn_h); if (allocated(error)) return
      call compare3(error, "hilyr", hilyr, hilyr_h); if (allocated(error)) return
      call compare3(error, "hslyr", hslyr, hslyr_h); if (allocated(error)) return
      do k = 1, NIL
         call compare3(error, "qicen", qicen(:, :, :, k), qicen_h(:, :, :, k)); if (allocated(error)) return
      end do
      call compare3(error, "qsnon", qsnon(:, :, :, 1), qsnon_h(:, :, :, 1)); if (allocated(error)) return
      do k = 1, 13
         call compare3(error, trim(names(k)), o_d(:, :, :, k), o_h(:, :, :, k)); if (allocated(error)) return
      end do

   contains

      subroutine host_sweep()
         integer :: ii, jj, nn
         do nn = 1, NCAT
            do jj = 1, NY
               do ii = 1, NX
                  call ice_thickness_changes_cell(ii, jj, nn, NX, NY, NCAT, NIL, 1, c, lb, DT, Tbot, fbot, fsnow, flatn, &
                                                  fsurfn, fcondtopn, fcondbotn, sicen, hin_h, hsn_h, hilyr_h, hslyr_h, qicen_h, &
                                                  qsnon_h, o_h(:, :, :, 1), o_h(:, :, :, 2), o_h(:, :, :, 3), o_h(:, :, :, 4), &
                                                  o_h(:, :, :, 5), o_h(:, :, :, 6), o_h(:, :, :, 7), o_h(:, :, :, 8), &
                                                  o_h(:, :, :, 9), o_h(:, :, :, 10), o_h(:, :, :, 11), o_h(:, :, :, 12), ierr_h)
                  call ice_conservation_check_vthermo_cell(ii, jj, nn, NX, NY, NCAT, c, DT, fsurfn, flatn, o_h(:, :, :, 8), &
                                                           fswint, fsnow, einit, o_h(:, :, :, 12), ICE_FERRMAX, o_h(:, :, :, 13), &
                                                           ierr_h)
               end do
            end do
         end do
      end subroutine host_sweep

   end subroutine run_device_case

   subroutine compare3(error, label, a, b)
      !! Cell-by-cell comparison, NaN-proof: a non-finite value on either side is itself a failure.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: label
      real(wp), intent(in) :: a(:, :, :), b(:, :, :)
      integer :: i, j, n
      character(len=200) :: msg
      real(wp) :: atol
      atol = DRTOL*maxval(abs(b))   ! scale-aware floor: 1e-12 of the field's largest value on the GPU build, 0 on CPU
      do n = 1, size(a, 3)
         do j = 1, size(a, 2)
            do i = 1, size(a, 1)
               if (.not. (ieee_is_finite(a(i, j, n)) .and. ieee_is_finite(b(i, j, n)))) then
                  write (msg, '(a,a,3(1x,i0))') "non-finite ", label, i, j, n
                  call test_failed(error, trim(msg))
                  return
               end if
               if (.not. near(a(i, j, n), b(i, j, n), DRTOL, atol)) then
                  write (msg, '(a,a,3(1x,i0),2es24.16)') "device != host ", label, i, j, n, a(i, j, n), b(i, j, n)
                  call test_failed(error, trim(msg))
                  return
               end if
            end do
         end do
      end do
   end subroutine compare3

end module test_cice_thickness_changes
