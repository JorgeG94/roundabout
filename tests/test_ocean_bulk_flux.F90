!! Unit tests for the Large & Yeager bulk air-sea flux kernel
!! (`rdb_ocean_bulk_flux`) — C3 of the OM3 forcing plan.
!!
!! Cases:
!!   * T1 neutral stability — a hand-constructed case with
!!     `SST == T_air` and `q_air == q_sat(SST)` (so `dtemp = dqr = 0`
!!     identically) and `z_wind = 10 m` (so the height-adjustment log
!!     terms vanish): the whole Monin-Obukhov iteration collapses to a
!!     fixed point on its very first pass, and the converged `cd/ch/ce`
!!     must equal the NEUTRAL 10 m coefficients computed directly from
!!     `bulk_flux_cdn10_neutral` to round-off. Also exercises "SST ==
!!     air T and saturation => sensible ~= 0" (here: exactly 0).
!!   * T2 zero relative wind — wind equals the surface current exactly:
!!     `taux`/`tauy` are EXACTLY zero (the stress is proportional to the
!!     un-floored relative-wind component), while `sensible` stays
!!     non-zero whenever `SST /= T_air` — the free-convection floor.
!!   * T3 minimum-wind floor — relative wind of exactly 0 and of exactly
!!     `BULK_WIND_MIN` give IDENTICAL transfer coefficients/sensible/
!!     latent (both evaluate the bulk formulae at the SAME floored wind
!!     speed) but DIFFERENT stress (only the second case has a non-zero
!!     relative-wind vector to project the stress onto).
!!   * T4 latent-heat sign — a sub-saturated air column over a warm sea
!!     evaporates: `latent < 0` (cooling the ocean) and the column
!!     function's raw `evap > 0` (meteorological convention, mass
!!     leaving upward), while the driver's `evap_massflux < 0` (MOM6
!!     convention, after the sign flip).
!!   * T5 stability ordering — an unstable (warm sea, `SST > T_air`)
!!     column has a LARGER final `ch` than an otherwise-identical stable
!!     (`SST < T_air`) column at the same `|SST - T_air|`.
!!   * T6 driver == column — `bulk_flux_driver_2d` on a small 2D grid
!!     reproduces `bulk_flux_column` cell-by-cell, including through a
!!     `!$acc` device data region (inert on a non-OpenACC build) —
!!     the `mem:separate` GPU parity gate.
!!   * T7 land cells — a land cell whose raw inputs are a mix of 0 and
!!     NaN (`ieee_value(..., ieee_quiet_nan)`) must produce ALL EIGHT
!!     driver outputs exactly 0.0 and finite (`bulk_flux_column` is
!!     never even called there), while an adjacent wet cell's outputs
!!     are unaffected.
!!   * T8 longwave — `bulk_flux_longwave_up` at `SST = 300 K`,
!!     `emissivity = 1` matches `sigma*300**4` to 1e-12 relative, and
!!     the result scales linearly with `emissivity`.
module test_ocean_bulk_flux
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_ocean_bulk_flux, only: bulk_flux_column, bulk_flux_driver_2d, &
                                  bulk_flux_cdn10_neutral, bulk_flux_longwave_up, &
                                  BULK_WIND_MIN, BULK_T0_KELVIN, BULK_STEFAN_BOLTZMANN
   implicit none
   private

   public :: collect_ocean_bulk_flux_tests

   real(wp), parameter :: SLP_STD = 101325.0_wp
      !! Standard sea-level pressure (Pa).
   real(wp), parameter :: Z_WIND = 10.0_wp
   real(wp), parameter :: Z_TA = 2.0_wp
   integer, parameter :: N_ITER = 5

contains

   subroutine collect_ocean_bulk_flux_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("bulk_flux_neutral_matches_hand_computed", &
                               test_neutral_matches_hand_computed), &
                  new_unittest("bulk_flux_zero_relative_wind_zero_stress", &
                               test_zero_wind_zero_stress), &
                  new_unittest("bulk_flux_min_wind_floor", test_min_wind_floor), &
                  new_unittest("bulk_flux_latent_sign_evaporation", &
                               test_latent_sign_evaporation), &
                  new_unittest("bulk_flux_unstable_ch_exceeds_stable", &
                               test_stability_ordering), &
                  new_unittest("bulk_flux_driver_matches_column_gpu", &
                               test_driver_matches_column), &
                  new_unittest("bulk_flux_land_cells_zero_with_nan_inputs", &
                               test_land_cells_zero_with_nan_inputs), &
                  new_unittest("bulk_flux_longwave_up_matches_stefan_boltzmann", &
                               test_longwave_up_matches_stefan_boltzmann) &
                  ]
   end subroutine collect_ocean_bulk_flux_tests

   ! -----------------------------------------------------------------
   ! T1 — neutral stability reproduces the hand-computed neutral coeffs
   ! -----------------------------------------------------------------
   subroutine test_neutral_matches_hand_computed(error)
      !! Two sub-cases, both built on `dtemp = SST_K - T_air = 0`
      !! exactly (bit-identical subtraction of equal values):
      !!
      !!   (1) `n_iter = 0`: the Monin-Obukhov loop never runs, so
      !!       `cd/ch/ce` are exactly the PRE-LOOP neutral-at-10m values
      !!       (`bulk_flux_cdn10_neutral` + the fixed 18.0/32.7 Stanton
      !!       and 34.6 Dalton multipliers).  `stab` at this point is
      !!       `merge(0.0_wp, 1.0_wp, dtemp > 0.0_wp)` -- `dtemp` is
      !!       exactly `+0.0` here (`x - x`), which is NOT `> 0`, so
      !!       `stab = 1.0` deterministically on every compiler, landing
      !!       on the "18.0" branch (see the comment at the call site) —
      !!       the "hand-computed case with neutral stability" the work
      !!       package asks for.
      !!       NOT testing this at the knife-edge of the FULL iteration:
      !!       `dqr` is not pinned to exact zero there (no closed-form,
      !!       bit-exact saturation humidity is available), and at true
      !!       neutrality the buoyancy term is itself ~0, so ITS sign
      !!       (and hence which of the 18.0/32.7 branches a later pass
      !!       takes) is a coin-flip on round-off — a genuine property
      !!       of this empirical stability split, not a test bug.
      !!   (2) `n_iter = N_ITER` (5, the OM3 value): `tstar =
      !!       (ctn/sqrt(cdn))*(-dtemp)` is `(...)*0 = 0` on every pass
      !!       (multiplying by an EXACT zero is exact regardless of the
      !!       other factor), so `ta_use` never moves off `T_air` and
      !!       `dtemp` STAYS exactly 0 through all 5 iterations —
      !!       independent of humidity. Hence `sensible = -rho_air*cp*
      !!       ch*dtemp*wv` is exactly 0 too, for ANY `q_air`: the
      !!       robust form of "SST == air T => sensible ~= 0".
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: U10 = 8.0_wp
      real(wp), parameter :: SST = 20.0_wp
      real(wp), parameter :: QA = 0.010_wp
      real(wp) :: ta_k, taux, tauy, sensible, latent, evap
      real(wp) :: cd, ch, ce, ustar
      real(wp) :: cdn10_expected, ch_expected, ce_expected

      ta_k = SST + BULK_T0_KELVIN

      ! (1) n_iter = 0.
      call bulk_flux_column(U10, 0.0_wp, 0.0_wp, 0.0_wp, ta_k, QA, SLP_STD, SST, &
                            Z_WIND, Z_TA, 0, &
                            taux, tauy, sensible, latent, evap, cd, ch, ce, ustar)

      cdn10_expected = bulk_flux_cdn10_neutral(U10)
      ! `dtemp = sst_k - t_air` is `x - x`, i.e. EXACTLY `+0.0`: exact
      ! neutrality.  The kernel's stability split lands this on the
      ! STABLE (18.0e-3) branch deterministically, on every compiler
      ! (`merge` on `dtemp > 0.0_wp`, never `sign` of a value that can
      ! carry a signed zero).
      ch_expected = 18.0e-3_wp*sqrt(cdn10_expected)
      ce_expected = 34.6e-3_wp*sqrt(cdn10_expected)

      call check(error, abs(cd - cdn10_expected) < 1.0e-12_wp*cdn10_expected, &
                 "neutral (n_iter=0) cd must match the hand-computed cdn10")
      if (allocated(error)) return
      call check(error, abs(ch - ch_expected) < 1.0e-12_wp*ch_expected, &
                 "neutral (n_iter=0) ch must match the hand-computed 18e-3*sqrt(cdn10)")
      if (allocated(error)) return
      call check(error, abs(ce - ce_expected) < 1.0e-12_wp*ce_expected, &
                 "neutral (n_iter=0) ce must match the hand-computed 34.6e-3*sqrt(cdn10)")
      if (allocated(error)) return
      call check(error, sensible == 0.0_wp, &
                 "sensible heat must be exactly 0 at SST==airT with n_iter=0")
      if (allocated(error)) return

      ! (2) n_iter = N_ITER (5) -- the dtemp == 0 invariant survives
      ! the full iteration regardless of humidity.
      call bulk_flux_column(U10, 0.0_wp, 0.0_wp, 0.0_wp, ta_k, QA, SLP_STD, SST, &
                            Z_WIND, Z_TA, N_ITER, &
                            taux, tauy, sensible, latent, evap, cd, ch, ce, ustar)
      call check(error, sensible == 0.0_wp, &
                 "sensible heat must stay exactly 0 at SST==airT through the full iteration")
   end subroutine test_neutral_matches_hand_computed

   ! -----------------------------------------------------------------
   ! T2 — zero relative wind => zero stress, non-zero free convection
   ! -----------------------------------------------------------------
   subroutine test_zero_wind_zero_stress(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: SST = 25.0_wp
      real(wp), parameter :: TA_K = 15.0_wp + BULK_T0_KELVIN
         !! Air much colder than SST -- a strongly unstable column, so
         !! the free-convection sensible flux is unambiguously non-zero.
      real(wp), parameter :: QA = 0.006_wp
      real(wp) :: taux, tauy, sensible, latent, evap, cd, ch, ce, ustar

      ! Wind == current (both zero) => relative wind vector is exactly
      ! zero; the wind-speed MAGNITUDE used by the transfer coefficients
      ! is still floored at BULK_WIND_MIN, but the stress is the floored
      ! magnitude TIMES the (zero) relative-wind component, hence 0.
      call bulk_flux_column(0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, TA_K, QA, SLP_STD, SST, &
                            Z_WIND, Z_TA, N_ITER, &
                            taux, tauy, sensible, latent, evap, cd, ch, ce, ustar)

      call check(error, taux == 0.0_wp, "taux must be exactly zero at zero relative wind")
      if (allocated(error)) return
      call check(error, tauy == 0.0_wp, "tauy must be exactly zero at zero relative wind")
      if (allocated(error)) return
      ! SST (25 degC) is much warmer than the air (15 degC): the ocean
      ! LOSES sensible heat even at dead calm (the free-convection
      ! floor keeps the exchange coefficients finite) -- "positive
      ! down" means this is a sizeable NEGATIVE number, not zero.
      call check(error, sensible < -1.0_wp, &
                 "sensible heat must stay non-zero (free-convection floor) at zero wind")
   end subroutine test_zero_wind_zero_stress

   ! -----------------------------------------------------------------
   ! T3 — the BULK_WIND_MIN floor: coefficients agree, stress differs
   ! -----------------------------------------------------------------
   subroutine test_min_wind_floor(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: SST = 18.0_wp
      real(wp), parameter :: TA_K = 16.0_wp + BULK_T0_KELVIN
      real(wp), parameter :: QA = 0.009_wp
      real(wp) :: taux0, tauy0, sens0, lat0, evap0, cd0, ch0, ce0, ustar0
      real(wp) :: taux1, tauy1, sens1, lat1, evap1, cd1, ch1, ce1, ustar1

      call bulk_flux_column(0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, TA_K, QA, SLP_STD, SST, &
                            Z_WIND, Z_TA, N_ITER, &
                            taux0, tauy0, sens0, lat0, evap0, cd0, ch0, ce0, ustar0)
      call bulk_flux_column(BULK_WIND_MIN, 0.0_wp, 0.0_wp, 0.0_wp, TA_K, QA, SLP_STD, SST, &
                            Z_WIND, Z_TA, N_ITER, &
                            taux1, tauy1, sens1, lat1, evap1, cd1, ch1, ce1, ustar1)

      call check(error, abs(cd0 - cd1) < 1.0e-14_wp, &
                 "cd at zero wind must equal cd at the floored wind speed")
      if (allocated(error)) return
      call check(error, abs(sens0 - sens1) < 1.0e-12_wp, &
                 "sensible at zero wind must equal sensible at the floored wind speed")
      if (allocated(error)) return
      call check(error, taux0 == 0.0_wp, "taux must be zero with no relative-wind component")
      if (allocated(error)) return
      call check(error, taux1 > 0.0_wp, &
                 "taux must be non-zero once the relative wind has a component")
   end subroutine test_min_wind_floor

   ! -----------------------------------------------------------------
   ! T4 — latent heat sign (sub-saturated air over a warm sea)
   ! -----------------------------------------------------------------
   subroutine test_latent_sign_evaporation(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: SST = 28.0_wp
      real(wp), parameter :: TA_K = 26.0_wp + BULK_T0_KELVIN
      real(wp), parameter :: QA = 0.012_wp
         !! Well below q_sat(28 degC) (~0.024 kg/kg) -- the air is dry
         !! relative to the warm sea, so the ocean evaporates.
      real(wp) :: taux, tauy, sensible, latent, evap, cd, ch, ce, ustar
      real(wp) :: wet(1, 1), taux_c(1, 1), tauy_c(1, 1), qsens(1, 1), qlat(1, 1)
      real(wp) :: evapm(1, 1), cd2(1, 1), ch2(1, 1), ce2(1, 1)
      real(wp) :: u_w(1, 1), v_w(1, 1), u_c(1, 1), v_c(1, 1)
      real(wp) :: ta2(1, 1), qa2(1, 1), slp2(1, 1), sst2(1, 1)

      call bulk_flux_column(6.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, TA_K, QA, SLP_STD, SST, &
                            Z_WIND, Z_TA, N_ITER, &
                            taux, tauy, sensible, latent, evap, cd, ch, ce, ustar)

      call check(error, latent < 0.0_wp, &
                 "latent heat must be negative (cooling) when the sea evaporates")
      if (allocated(error)) return
      call check(error, evap > 0.0_wp, &
                 "column evap must be positive (meteorological convention: mass leaving up)")
      if (allocated(error)) return

      ! Driver boundary: the same case through the 2D kernel must flip
      ! evap to the MOM6 "<= 0, mass leaving" convention.
      u_w = 6.0_wp; v_w = 0.0_wp; u_c = 0.0_wp; v_c = 0.0_wp
      ta2 = TA_K; qa2 = QA; slp2 = SLP_STD; sst2 = SST; wet = 1.0_wp
      call bulk_flux_driver_2d(1, 1, u_w, v_w, u_c, v_c, ta2, qa2, slp2, sst2, wet, &
                               Z_WIND, Z_TA, N_ITER, &
                               taux_c, tauy_c, qsens, qlat, evapm, cd2, ch2, ce2)
      call check(error, evapm(1, 1) < 0.0_wp, &
                 "driver evap_massflux must be negative (MOM6 convention) when the sea evaporates")
   end subroutine test_latent_sign_evaporation

   ! -----------------------------------------------------------------
   ! T5 — unstable (warm sea) gives a larger ch than stable (cold sea)
   ! -----------------------------------------------------------------
   subroutine test_stability_ordering(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: U10 = 7.0_wp
      real(wp), parameter :: TA_K = 20.0_wp + BULK_T0_KELVIN
      real(wp), parameter :: QA = 0.010_wp
      real(wp) :: taux, tauy, sensible, latent, evap, ustar
      real(wp) :: cd_u, ch_u, ce_u, cd_s, ch_s, ce_s

      ! Unstable: SST 5 degC warmer than the air.
      call bulk_flux_column(U10, 0.0_wp, 0.0_wp, 0.0_wp, TA_K, QA, SLP_STD, 25.0_wp, &
                            Z_WIND, Z_TA, N_ITER, &
                            taux, tauy, sensible, latent, evap, cd_u, ch_u, ce_u, ustar)
      ! Stable: SST 5 degC colder than the air (same |difference|).
      call bulk_flux_column(U10, 0.0_wp, 0.0_wp, 0.0_wp, TA_K, QA, SLP_STD, 15.0_wp, &
                            Z_WIND, Z_TA, N_ITER, &
                            taux, tauy, sensible, latent, evap, cd_s, ch_s, ce_s, ustar)

      call check(error, ch_u > ch_s, &
                 "unstable (warm sea) ch must exceed stable (cold sea) ch at the same |dT|")
   end subroutine test_stability_ordering

   ! -----------------------------------------------------------------
   ! T6 — driver matches the column function, including a device
   !      data region (inert without OpenACC, a real device round
   !      trip with it).
   ! -----------------------------------------------------------------
   subroutine test_driver_matches_column(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 4, NY = 3
      real(wp) :: u_w(NX, NY), v_w(NX, NY), u_c(NX, NY), v_c(NX, NY)
      real(wp) :: ta(NX, NY), qa(NX, NY), slp(NX, NY), sst(NX, NY), wet(NX, NY)
      real(wp) :: taux_c(NX, NY), tauy_c(NX, NY), qsens(NX, NY), qlat(NX, NY)
      real(wp) :: evapm(NX, NY), cd(NX, NY), ch(NX, NY), ce(NX, NY)
      real(wp) :: taux_h, tauy_h, sens_h, lat_h, evap_h, cd_h, ch_h, ce_h, ustar_h
      integer :: i, j
      real(wp) :: tol

      do j = 1, NY
         do i = 1, NX
            u_w(i, j) = 2.0_wp + 0.7_wp*real(i, wp) - 0.3_wp*real(j, wp)
            v_w(i, j) = -1.0_wp + 0.2_wp*real(i, wp)
            u_c(i, j) = 0.05_wp*real(j, wp)
            v_c(i, j) = -0.03_wp*real(i, wp)
            ta(i, j) = 290.0_wp + 0.5_wp*real(i, wp) - 0.2_wp*real(j, wp)
            qa(i, j) = 0.008_wp + 0.0005_wp*real(i, wp)
            slp(i, j) = SLP_STD - 50.0_wp*real(j, wp)
            sst(i, j) = 18.0_wp + 0.3_wp*real(i, wp)
            wet(i, j) = merge(0.0_wp, 1.0_wp, i == 1 .and. j == 1)
         end do
      end do

      !$omp target enter data map(to: u_w, v_w, u_c, v_c, ta, qa, slp, sst, wet) &
      !$omp&   map(alloc: taux_c, tauy_c, qsens, qlat, evapm, cd, ch, ce)
      call bulk_flux_driver_2d(NX, NY, u_w, v_w, u_c, v_c, ta, qa, slp, sst, wet, &
                               Z_WIND, Z_TA, N_ITER, &
                               taux_c, tauy_c, qsens, qlat, evapm, cd, ch, ce)
      !$omp target update from(taux_c, tauy_c, qsens, qlat, evapm, cd, ch, ce)
      !$omp target exit data map(delete: u_w, v_w, u_c, v_c, ta, qa, slp, sst, wet, &
      !$omp&   taux_c, tauy_c, qsens, qlat, evapm, cd, ch, ce)

      do j = 1, NY
         do i = 1, NX
            call bulk_flux_column(u_w(i, j), v_w(i, j), u_c(i, j), v_c(i, j), &
                                  ta(i, j), qa(i, j), slp(i, j), sst(i, j), &
                                  Z_WIND, Z_TA, N_ITER, &
                                  taux_h, tauy_h, sens_h, lat_h, evap_h, &
                                  cd_h, ch_h, ce_h, ustar_h)
            if (i == 1 .and. j == 1) then
               call check(error, taux_c(i, j) == 0.0_wp, "masked cell must be exactly zero")
               if (allocated(error)) return
               cycle
            end if
            tol = 1.0e-12_wp*max(abs(taux_h), 1.0e-30_wp)
            call check(error, abs(taux_c(i, j) - taux_h) < tol, &
                       "driver taux must match the column function")
            if (allocated(error)) return
            tol = 1.0e-12_wp*max(abs(sens_h), 1.0e-30_wp)
            call check(error, abs(qsens(i, j) - sens_h) < tol, &
                       "driver q_sens must match the column function")
            if (allocated(error)) return
            tol = 1.0e-12_wp*max(abs(evap_h), 1.0e-30_wp)
            call check(error, abs(evapm(i, j) - (-evap_h)) < tol, &
                       "driver evap_massflux must match -1 * the column function's evap")
            if (allocated(error)) return
            tol = 1.0e-12_wp*max(abs(ch_h), 1.0e-30_wp)
            call check(error, abs(ch(i, j) - ch_h) < tol, &
                       "driver ch must match the column function")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_driver_matches_column

   ! -----------------------------------------------------------------
   ! T7 — land cells: explicit zero outputs, never a masked NaN
   ! -----------------------------------------------------------------
   subroutine test_land_cells_zero_with_nan_inputs(error)
      !! A land cell (`wet_mask == 0`) with raw inputs that are a mix
      !! of 0 and NaN must come back with ALL EIGHT driver outputs
      !! exactly 0.0 and finite -- `bulk_flux_driver_2d` must branch on
      !! `wet_mask` and skip the `bulk_flux_column` call outright
      !! rather than multiply a (possibly NaN) column result by a zero
      !! mask (`NaN*0.0 == NaN`).  The adjacent wet cell's outputs must
      !! be unaffected and must match a direct `bulk_flux_column` call.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 2, NY = 1
      real(wp) :: u_w(NX, NY), v_w(NX, NY), u_c(NX, NY), v_c(NX, NY)
      real(wp) :: ta(NX, NY), qa(NX, NY), slp(NX, NY), sst(NX, NY), wet(NX, NY)
      real(wp) :: taux_c(NX, NY), tauy_c(NX, NY), qsens(NX, NY), qlat(NX, NY)
      real(wp) :: evapm(NX, NY), cd(NX, NY), ch(NX, NY), ce(NX, NY)
      real(wp) :: taux_h, tauy_h, sens_h, lat_h, evap_h, cd_h, ch_h, ce_h, ustar_h
      real(wp) :: nan_val
      real(wp), parameter :: U10 = 7.0_wp
      real(wp), parameter :: SST_WET = 19.0_wp
      real(wp), parameter :: TA_WET = 17.0_wp + BULK_T0_KELVIN
      real(wp), parameter :: QA_WET = 0.009_wp

      nan_val = ieee_value(1.0_wp, ieee_quiet_nan)

      ! Cell 1: land, inputs a mix of 0 and NaN.
      u_w(1, 1) = nan_val; v_w(1, 1) = 0.0_wp
      u_c(1, 1) = nan_val; v_c(1, 1) = 0.0_wp
      ta(1, 1) = nan_val; qa(1, 1) = 0.0_wp
      slp(1, 1) = nan_val; sst(1, 1) = nan_val
      wet(1, 1) = 0.0_wp

      ! Cell 2: wet, ordinary inputs.
      u_w(2, 1) = U10; v_w(2, 1) = 0.0_wp
      u_c(2, 1) = 0.0_wp; v_c(2, 1) = 0.0_wp
      ta(2, 1) = TA_WET; qa(2, 1) = QA_WET
      slp(2, 1) = SLP_STD; sst(2, 1) = SST_WET
      wet(2, 1) = 1.0_wp

      call bulk_flux_driver_2d(NX, NY, u_w, v_w, u_c, v_c, ta, qa, slp, sst, wet, &
                               Z_WIND, Z_TA, N_ITER, &
                               taux_c, tauy_c, qsens, qlat, evapm, cd, ch, ce)

      call check(error, taux_c(1, 1) == 0.0_wp .and. ieee_is_finite(taux_c(1, 1)), &
                 "land taux_cell must be exactly 0 and finite")
      if (allocated(error)) return
      call check(error, tauy_c(1, 1) == 0.0_wp .and. ieee_is_finite(tauy_c(1, 1)), &
                 "land tauy_cell must be exactly 0 and finite")
      if (allocated(error)) return
      call check(error, qsens(1, 1) == 0.0_wp .and. ieee_is_finite(qsens(1, 1)), &
                 "land q_sens must be exactly 0 and finite")
      if (allocated(error)) return
      call check(error, qlat(1, 1) == 0.0_wp .and. ieee_is_finite(qlat(1, 1)), &
                 "land q_lat must be exactly 0 and finite")
      if (allocated(error)) return
      call check(error, evapm(1, 1) == 0.0_wp .and. ieee_is_finite(evapm(1, 1)), &
                 "land evap_massflux must be exactly 0 and finite")
      if (allocated(error)) return
      call check(error, cd(1, 1) == 0.0_wp .and. ieee_is_finite(cd(1, 1)), &
                 "land cd must be exactly 0 and finite")
      if (allocated(error)) return
      call check(error, ch(1, 1) == 0.0_wp .and. ieee_is_finite(ch(1, 1)), &
                 "land ch must be exactly 0 and finite")
      if (allocated(error)) return
      call check(error, ce(1, 1) == 0.0_wp .and. ieee_is_finite(ce(1, 1)), &
                 "land ce must be exactly 0 and finite")
      if (allocated(error)) return

      ! The wet cell must be unaffected by the land neighbour.
      call bulk_flux_column(U10, 0.0_wp, 0.0_wp, 0.0_wp, TA_WET, QA_WET, SLP_STD, SST_WET, &
                            Z_WIND, Z_TA, N_ITER, &
                            taux_h, tauy_h, sens_h, lat_h, evap_h, cd_h, ch_h, ce_h, ustar_h)
      call check(error, abs(taux_c(2, 1) - taux_h) < 1.0e-12_wp*max(abs(taux_h), 1.0e-30_wp), &
                 "wet cell taux must be unaffected by a land neighbour")
      if (allocated(error)) return
      call check(error, abs(ch(2, 1) - ch_h) < 1.0e-12_wp*max(abs(ch_h), 1.0e-30_wp), &
                 "wet cell ch must be unaffected by a land neighbour")
   end subroutine test_land_cells_zero_with_nan_inputs

   ! -----------------------------------------------------------------
   ! T8 — longwave: Stefan-Boltzmann law + linear emissivity scaling
   ! -----------------------------------------------------------------
   subroutine test_longwave_up_matches_stefan_boltzmann(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: SST_K = 300.0_wp
      real(wp) :: lw_full, lw_half, expected

      lw_full = bulk_flux_longwave_up(SST_K, 1.0_wp)
      expected = BULK_STEFAN_BOLTZMANN*SST_K**4
      call check(error, abs(lw_full - expected) < 1.0e-12_wp*expected, &
                 "longwave_up at emissivity=1 must match sigma*SST**4 to 1e-12 relative")
      if (allocated(error)) return

      lw_half = bulk_flux_longwave_up(SST_K, 0.5_wp)
      call check(error, abs(lw_half - 0.5_wp*lw_full) < 1.0e-12_wp*lw_full, &
                 "longwave_up must scale linearly with emissivity")
   end subroutine test_longwave_up_matches_stefan_boltzmann

end module test_ocean_bulk_flux
