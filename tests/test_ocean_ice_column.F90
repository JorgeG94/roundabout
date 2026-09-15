!! Analytic tests for the single-column Winton NkIce-layer sea-ice
!! thermodynamics (`rdb_ice_column`, sea-ice PR 3a) — the SIS2 port
!! (`ice_temp_SIS2`/`laytemp_SIS2`/`update_lay_enth`) plus the mass
!! bookkeeping (`rdb_ice_mass`) and shortwave optics (`rdb_ice_optics`)
!! it composes.
!!
!! Cases (SPEC_ice-pr3a-column.md §7):
!!   * `stefan_growth` — the TRAP-#1/#2/#3/#5 physics gate: growth of a
!!     conduction-limited ice cap under a cold stiff-pinned surface
!!     matches the classical Stefan analytic solution to < 1% at days
!!     1/5/10/30, plus a mass-bookkeeping identity (net ocean<->ice
!!     water flux equals the ice mass change).
!!   * `conduction_solve_energy_closure` (formerly `column_energy_closure`
!!     — renamed, PR-60: the old name read as "the column conserves
!!     energy", which it does not test) — `ice_temp_sis2`'s own internal
!!     energy-closure identity holds to the 1e-9 SIS2 budget every step
!!     over a 30-day run. Covers the CONDUCTION SOLVE ONLY — not the
!!     resize/peel/rebalance chain `ice_column_step` runs after it
!!     (`column_step_energy_closes_cold` below is the end-to-end case);
!!     deliberately excludes the liq-lim clamp (SIS2 checks it, this
!!     port does not — `rdb_ice_column.F90:310-316`); the identity omits
!!     the `-e_extra_sum` term and is exact only on a snow-free,
!!     non-pinning column (every case in this file), which is why it
!!     never surfaces.
!!   * `column_step_energy_closes_cold` / `_with_sw` / `_thick` (PR-60,
!!     RESUME item b) — the END-TO-END cold-column total-energy budget:
!!     drives the full `ice_column_step` (optics + conduction + resize +
!!     rebalance) and pins `Delta(E_col) ==
!!     dtt*((1-albedo)*sw_dn - sw_thru + fb) + h2o_ocn_to_ice*enth_ocean`
!!     to SIS2's 1e-9 sum-of-magnitudes budget, at sw_dn in {0,200},
!!     fb in {0,5}, h0 in {1,3} m — the property `conduction_solve_
!!     energy_closure`'s name used to (falsely) claim.
!!   * `equilibrium_fixed_point` — a column seeded exactly at the
!!     analytic equilibrium thickness (conductive supply = constant
!!     ocean heat flux) drifts by < 1e-6 m over 10 days.
!!   * `rebalance_conserves` — `ice_rebalance_layers` on an unequal
!!     2-layer column conserves total ice mass, enthalpy, and salt.
!!   * `optics_csim4` — `ice_optics_csim4` closed forms (bare cold ice,
!!     snow-covered, melting bare ice) plus the shortwave-partition
!!     identity.
!!   * `laytemp_quadratic` — `laytemp_sis2` and `update_lay_enth` satisfy
!!     their own implicit energy-balance identities to round-off.
!!   * `snow_ice_archimedes` — PR 27: `ice_snow_ice_flood` on a flooded
!!     column hits the Archimedes equilibrium EXACTLY in one call, and
!!     conserves mass/enthalpy/salt to round-off (the physics gate).
!!   * `snow_ice_threshold_and_clamp` — PR 27: below-threshold columns
!!     are bit-unchanged; the SIS2 `min(...)` clamp is provably dead
!!     code in this pond-free model (the branch gate).
!!   * `snow_ice_column_device` — PR 27: full `ice_thermo_columns` with
!!     snow present, `do_snow_ice=.true.` — the GPU `mem:separate`
!!     canary for `snow_to_ice` (the wiring gate).
!!   * `snow_ice_disabled_bitident` — PR 27: `do_snow_ice=.false.` is a
!!     bit-for-bit no-op against the same trajectory (the default-off
!!     gate).
module test_ocean_ice_column
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ice_state, only: ocean_sea_ice_t
   use rdb_ice_enthalpy, only: ICE_LAT_FUS, ICE_CP_ICE, ICE_CP_WATER, ICE_DTF_DS, &
                               ice_enth_from_ts, ice_t_freeze, ice_enthalpy_liquid, &
                               ice_temp_from_en_s
   use rdb_ice_column, only: ICE_K_ICE, ICE_RHO_ICE, ICE_RHO_SNOW, ICE_RHO_OCEAN, &
                             ICE_BULK_SALINITY, ICE_NK_MAX, &
                             ICE_TEMP_RANGE_EST, &
                             laytemp_sis2, update_lay_enth, ice_temp_sis2, &
                             ice_thermo_columns, ice_column_step
   use rdb_ice_mass, only: ice_bottom_freeze, ice_top_melt_peel, ice_bottom_melt_peel, &
                           ice_rebalance_layers, ice_snow_ice_flood
   use rdb_ice_optics, only: ice_optics_csim4, ICE_ALB_SNOW, ICE_ALB_ICE, &
                             ICE_PEN_ICE, ICE_OPT_DEP_ICE, ICE_T_RANGE_MELT, &
                             ICE_SNOW_PATCH
   implicit none
   private

   public :: collect_ocean_ice_column_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NX_PHYS = 6
   integer, parameter :: NY_PHYS = 4
   integer, parameter :: NK = 2
   real(wp), parameter :: DAY_S = 86400.0_wp
   real(wp), parameter :: S_INIT = 35.0_wp
      !! Sea-surface salinity (PSU) fed to `ice_column_step`'s TRAP-#1
      !! liquid-ocean enthalpy in the cold-energy-closure cases (PR-60).

contains

   subroutine collect_ocean_ice_column_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("stefan_growth", test_stefan_growth), &
                  new_unittest("melt_ablation", test_melt_ablation), &
                  new_unittest("conduction_solve_energy_closure", &
                               test_conduction_solve_energy_closure), &
                  new_unittest("equilibrium_fixed_point", test_equilibrium_fixed_point), &
                  new_unittest("rebalance_conserves", test_rebalance_conserves), &
                  new_unittest("optics_csim4", test_optics_csim4), &
                  new_unittest("laytemp_quadratic", test_laytemp_quadratic), &
                  new_unittest("snow_ice_archimedes", test_snow_ice_archimedes), &
                  new_unittest("snow_ice_threshold_and_clamp", test_snow_ice_threshold_and_clamp), &
                  new_unittest("snow_ice_column_device", test_snow_ice_column_device), &
                  new_unittest("snow_ice_disabled_bitident", test_snow_ice_disabled_bitident), &
                  new_unittest("column_step_energy_closes_cold", &
                               test_column_step_energy_closes_cold), &
                  new_unittest("column_step_energy_closes_cold_with_sw", &
                               test_column_step_energy_closes_cold_with_sw), &
                  new_unittest("column_step_energy_closes_cold_thick", &
                               test_column_step_energy_closes_cold_thick) &
                  ]
   end subroutine collect_ocean_ice_column_tests

   ! -----------------------------------------------------------------
   ! Shared helpers
   ! -----------------------------------------------------------------

   pure function stefan_thickness_analytic(t_seconds, delta_t_mag) result(h)
      !! Classic Stefan solution: h(t) = sqrt(2*k_ice*dT*t/(rho_ice*Lat_fus)).
      real(wp), intent(in) :: t_seconds, delta_t_mag
      real(wp) :: h
      h = sqrt(2.0_wp*ICE_K_ICE*delta_t_mag*t_seconds/(ICE_RHO_ICE*ICE_LAT_FUS))
   end function stefan_thickness_analytic

   subroutine make_initial_column(h_ice0, tfw, tsurf_target, has_target, &
                                  m_ice, m_snow, enth_ice_bu, enth_snow_pt, sal_ice_bu)
      !! Seed a snow-free ice column of thickness `h_ice0` (m) at bulk
      !! salinity `ICE_BULK_SALINITY`, BOTTOM-UP state order. When
      !! `has_target`, seeds a LINEAR temperature profile from
      !! `tsurf_target` (top) to `tfw` (base) — the quasi-steady profile
      !! the Stefan solution itself assumes (an isothermal IC creates a
      !! spurious one-step transient the growth-rate diagnostic would
      !! misread as extra melt/freeze; prototype `make_initial_column`
      !! docstring). Without a target, seeds isothermal at `tfw`.
      real(wp), intent(in) :: h_ice0, tfw, tsurf_target
      logical, intent(in) :: has_target
      real(wp), intent(out) :: m_ice, m_snow, enth_snow_pt
      real(wp), intent(out) :: enth_ice_bu(NK), sal_ice_bu(NK)

      integer :: k, k_td
      real(wp) :: frac, t_k

      m_ice = h_ice0*ICE_RHO_ICE
      m_snow = 0.0_wp
      sal_ice_bu(:) = ICE_BULK_SALINITY
      enth_snow_pt = ice_enth_from_ts(0.0_wp, 0.0_wp)

      do k = 1, NK
         ! k is BOTTOM-UP (1 = bed); the seeding profile is defined
         ! top-down, so convert via k_td = NK+1-k (k_td=1 is the top
         ! layer, matching the make_initial_column prototype).
         k_td = NK + 1 - k
         if (has_target) then
            frac = (real(k_td, wp) - 0.5_wp)/real(NK, wp)
            t_k = tsurf_target + frac*(tfw - tsurf_target)
         else
            t_k = tfw
         end if
         enth_ice_bu(k) = ice_enth_from_ts(t_k, ICE_BULK_SALINITY)
      end do
   end subroutine make_initial_column

   ! -----------------------------------------------------------------
   ! Cases
   ! -----------------------------------------------------------------

   subroutine test_stefan_growth(error)
      !! TRAP-#1/#2/#3/#5 physics gate. Seed a thin ice cap at the
      !! Stefan-analytic thickness for t0=1 day, then advance hourly to
      !! day 32 under a cold stiff-pinned surface (tsurf=-20 degC,
      !! tfw=-1.9 degC, fb=0) via the FULL `ice_thermo_columns` device
      !! driver. At the first step crossing each of days {1,5,10,30},
      !! assert the modeled thickness matches the analytic Stefan
      !! solution (evaluated at the actual elapsed time) to < 1%
      !! relative. Also assert the accumulated ocn<->ice water-mass
      !! identity: Σh2o_ocn_to_ice - Σh2o_ice_to_ocn == Δm_ice, to 1e-9
      !! relative.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_sea_ice_t) :: ice
      real(wp), parameter :: TFW = -1.9_wp
      real(wp), parameter :: TSURF_TARGET = -20.0_wp
      real(wp), parameter :: DELTA_T_MAG = abs(TSURF_TARGET - TFW)
      real(wp), parameter :: DTT = 3600.0_wp
      integer, parameter :: N_STEPS = int(32.0_wp*DAY_S/DTT)
      integer, parameter :: N_DAYS_CHECK = 4
      integer, parameter :: DAYS_CHECK(N_DAYS_CHECK) = [1, 5, 10, 30]

      real(wp), allocatable :: wet_mask(:, :), sf_0(:, :), dsf_dt(:, :), sw_dn(:, :), fprec(:, :)
      real(wp), allocatable :: tfw_arr(:, :), fb_arr(:, :), sst_arr(:, :), s_surf_arr(:, :)
      real(wp), allocatable :: tsurf_out(:, :, :), h2o_ocn_to_ice(:, :, :), &
                               h2o_ice_to_ocn(:, :, :), heat_to_ocn(:, :, :), sw_thru(:, :, :), &
                               snow_to_ice(:, :, :)
      real(wp) :: m_ice0, m_ice_seed, m_snow_seed, enth_snow_seed
      real(wp) :: enth_ice_seed(NK), sal_ice_seed(NK)
      real(wp) :: t_elapsed, h_model, h_analytic, rel_err
      real(wp) :: sum_ocn_to_ice, sum_ice_to_ocn, mass_identity, tol
      integer :: ip, jp, nx, ny, step, d, idx
      logical :: recorded(N_DAYS_CHECK)
      character(len=64) :: msg

      checks: block

         call grid%init(NX_PHYS, NY_PHYS, NGHOST, 1.0_wp, 1.0_wp)
         nx = grid%nx_total
         ny = grid%ny_total
         ice%enable = .true.
         ice%ncat = 1
         ice%nk_ice = NK
         call ice%init(grid)

         allocate (wet_mask(nx, ny), source=1.0_wp)
         allocate (sf_0(nx, ny), dsf_dt(nx, ny), sw_dn(nx, ny), fprec(nx, ny))
         allocate (tfw_arr(nx, ny), fb_arr(nx, ny), sst_arr(nx, ny), s_surf_arr(nx, ny))
         allocate (tsurf_out(nx, ny, 1), h2o_ocn_to_ice(nx, ny, 1), &
                   h2o_ice_to_ocn(nx, ny, 1), heat_to_ocn(nx, ny, 1), sw_thru(nx, ny, 1), &
                   snow_to_ice(nx, ny, 1))

         ! Stiff-SEB pin: SF(T) = sf_0 + dsf_dt*T, upward-positive,
         ! increasing in T, so tsurf tracks tsurf_target to O(1/stiff).
         dsf_dt = 1.0e6_wp
         sf_0 = -1.0e6_wp*TSURF_TARGET
         sw_dn = 0.0_wp
         fprec = 0.0_wp
         tfw_arr = TFW
         fb_arr = 0.0_wp
         sst_arr = TFW
         s_surf_arr = ICE_BULK_SALINITY

         call make_initial_column(stefan_thickness_analytic(DAY_S, DELTA_T_MAG), TFW, &
                                  TSURF_TARGET, .true., &
                                  m_ice_seed, m_snow_seed, enth_ice_seed, &
                                  enth_snow_seed, sal_ice_seed)
         m_ice0 = m_ice_seed

         ip = grid%nghost + 2
         jp = grid%nghost + 2

         !$acc enter data copyin(ice)
         call ice%enter_data()
         !$acc enter data copyin(wet_mask, sf_0, dsf_dt, sw_dn, fprec, tfw_arr, fb_arr, &
         !$acc                   sst_arr, s_surf_arr, tsurf_out, h2o_ocn_to_ice, &
         !$acc                   h2o_ice_to_ocn, heat_to_ocn, sw_thru, snow_to_ice)

         ice%m_ice(ip, jp, 1) = m_ice_seed
         ice%m_snow(ip, jp, 1) = m_snow_seed
         ice%enth_ice(ip, jp, 1, :) = enth_ice_seed
         ice%enth_snow(ip, jp, 1, 1) = enth_snow_seed
         ice%sal_ice(ip, jp, 1, :) = sal_ice_seed
         !$acc update device(ice%m_ice, ice%m_snow, ice%enth_ice, ice%enth_snow, ice%sal_ice)

         t_elapsed = DAY_S
         recorded = .false.
         sum_ocn_to_ice = 0.0_wp
         sum_ice_to_ocn = 0.0_wp

         do step = 1, N_STEPS
            call ice_thermo_columns(grid%nghost, nx, ny, 1, NK, DTT, .false., wet_mask, &
                                    ice%m_ice, ice%m_snow, ice%enth_ice, ice%enth_snow, &
                                    ice%sal_ice, sf_0, dsf_dt, sw_dn, fprec, tfw_arr, fb_arr, &
                                    sst_arr, s_surf_arr, &
                                    tsurf_out, h2o_ocn_to_ice, h2o_ice_to_ocn, &
                                    heat_to_ocn, sw_thru, snow_to_ice)
            t_elapsed = t_elapsed + DTT

            !$acc update self(ice%m_ice, h2o_ocn_to_ice, h2o_ice_to_ocn)
            sum_ocn_to_ice = sum_ocn_to_ice + h2o_ocn_to_ice(ip, jp, 1)
            sum_ice_to_ocn = sum_ice_to_ocn + h2o_ice_to_ocn(ip, jp, 1)

            do idx = 1, N_DAYS_CHECK
               d = DAYS_CHECK(idx)
               if (.not. recorded(idx) .and. t_elapsed >= real(d, wp)*DAY_S - 1.0e-6_wp) then
                  recorded(idx) = .true.
                  h_model = ice%m_ice(ip, jp, 1)/ICE_RHO_ICE
                  h_analytic = stefan_thickness_analytic(t_elapsed, DELTA_T_MAG)
                  rel_err = abs(h_model - h_analytic)/h_analytic
                  write (msg, "(A,I0,A,ES12.5,A,ES12.5)") &
                     "stefan day ", d, ": rel_err=", rel_err, " h_model=", h_model
                  call check(error, rel_err < 0.01_wp, trim(msg))
                  if (allocated(error)) exit checks
               end if
            end do
         end do

         call check(error, all(recorded), "all four Stefan checkpoints must be recorded")
         if (allocated(error)) exit checks

         ! Mass-bookkeeping identity: net ocn<->ice water flux == Δm_ice.
         mass_identity = sum_ocn_to_ice - sum_ice_to_ocn
         tol = 1.0e-9_wp*max(abs(ice%m_ice(ip, jp, 1) - m_ice0), 1.0_wp)
         call check(error, abs(mass_identity - (ice%m_ice(ip, jp, 1) - m_ice0)) <= tol, &
                    "accumulated ocn<->ice water flux must equal the ice mass change")

      end block checks

      !$acc exit data delete(wet_mask, sf_0, dsf_dt, sw_dn, fprec, tfw_arr, fb_arr, &
      !$acc                  sst_arr, s_surf_arr, tsurf_out, h2o_ocn_to_ice, &
      !$acc                  h2o_ice_to_ocn, heat_to_ocn, sw_thru, snow_to_ice)
      call ice%exit_data()
      !$acc exit data delete(ice)
      call ice%destroy()
   end subroutine test_stefan_growth

   subroutine test_melt_ablation(error)
      !! Melt-side coverage gate (complements the cold Stefan growth
      !! test, which never exercises the warm branches). Seed a 1.5 m cap
      !! ISOTHERMAL at tfw (so the fixed-tfw base BC drives no spurious
      !! basal freezing), then force a WARM surface (tsurf_target=+2 degC,
      !! above the melt point, so the tsurf>tsf melting-pin branch banks
      !! the excess into surface melt) plus a positive ocean base flux
      !! (fb=5 W/m^2, basal melt). Over a 12-day hourly run assert: the
      !! ice thins (net m_ice decrease), ice->ocean meltwater is produced
      !! (`ice_top_melt_peel`/`ice_bottom_melt_peel` fire), and the
      !! ocn<->ice water-mass identity holds. Without this case the entire
      !! warm/melt half of the port runs in no test.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_sea_ice_t) :: ice
      real(wp), parameter :: TFW = -1.9_wp
      real(wp), parameter :: TSURF_TARGET = 2.0_wp
      real(wp), parameter :: H_INIT = 1.5_wp
      real(wp), parameter :: FB = 5.0_wp
      real(wp), parameter :: DTT = 3600.0_wp
      integer, parameter :: N_STEPS = int(12.0_wp*DAY_S/DTT)

      real(wp), allocatable :: wet_mask(:, :), sf_0(:, :), dsf_dt(:, :), sw_dn(:, :), fprec(:, :)
      real(wp), allocatable :: tfw_arr(:, :), fb_arr(:, :), sst_arr(:, :), s_surf_arr(:, :)
      real(wp), allocatable :: tsurf_out(:, :, :), h2o_ocn_to_ice(:, :, :), &
                               h2o_ice_to_ocn(:, :, :), heat_to_ocn(:, :, :), sw_thru(:, :, :), &
                               snow_to_ice(:, :, :)
      real(wp) :: m_ice0, m_ice_seed, m_snow_seed, enth_snow_seed
      real(wp) :: enth_ice_seed(NK), sal_ice_seed(NK)
      real(wp) :: m_ice_end, sum_ocn_to_ice, sum_ice_to_ocn, mass_identity, tol
      integer :: ip, jp, nx, ny, step

      checks: block

         call grid%init(NX_PHYS, NY_PHYS, NGHOST, 1.0_wp, 1.0_wp)
         nx = grid%nx_total
         ny = grid%ny_total
         ice%enable = .true.
         ice%ncat = 1
         ice%nk_ice = NK
         call ice%init(grid)

         allocate (wet_mask(nx, ny), source=1.0_wp)
         allocate (sf_0(nx, ny), dsf_dt(nx, ny), sw_dn(nx, ny), fprec(nx, ny))
         allocate (tfw_arr(nx, ny), fb_arr(nx, ny), sst_arr(nx, ny), s_surf_arr(nx, ny))
         allocate (tsurf_out(nx, ny, 1), h2o_ocn_to_ice(nx, ny, 1), &
                   h2o_ice_to_ocn(nx, ny, 1), heat_to_ocn(nx, ny, 1), sw_thru(nx, ny, 1), &
                   snow_to_ice(nx, ny, 1))

         ! Warm stiff-SEB pin: tsurf tracks +2 degC; the tsurf>tsf melting
         ! branch pins to tsf and banks the excess into surface melt.
         dsf_dt = 1.0e6_wp
         sf_0 = -1.0e6_wp*TSURF_TARGET
         sw_dn = 0.0_wp
         fprec = 0.0_wp
         tfw_arr = TFW
         fb_arr = FB
         sst_arr = TFW
         s_surf_arr = ICE_BULK_SALINITY

         ! Isothermal-at-tfw cap (no cold interior => no basal freezing).
         call make_initial_column(H_INIT, TFW, TFW, .true., &
                                  m_ice_seed, m_snow_seed, enth_ice_seed, &
                                  enth_snow_seed, sal_ice_seed)
         m_ice0 = m_ice_seed

         ip = grid%nghost + 2
         jp = grid%nghost + 2

         !$acc enter data copyin(ice)
         call ice%enter_data()
         !$acc enter data copyin(wet_mask, sf_0, dsf_dt, sw_dn, fprec, tfw_arr, fb_arr, &
         !$acc                   sst_arr, s_surf_arr, tsurf_out, h2o_ocn_to_ice, &
         !$acc                   h2o_ice_to_ocn, heat_to_ocn, sw_thru, snow_to_ice)

         ice%m_ice(ip, jp, 1) = m_ice_seed
         ice%m_snow(ip, jp, 1) = m_snow_seed
         ice%enth_ice(ip, jp, 1, :) = enth_ice_seed
         ice%enth_snow(ip, jp, 1, 1) = enth_snow_seed
         ice%sal_ice(ip, jp, 1, :) = sal_ice_seed
         !$acc update device(ice%m_ice, ice%m_snow, ice%enth_ice, ice%enth_snow, ice%sal_ice)

         sum_ocn_to_ice = 0.0_wp
         sum_ice_to_ocn = 0.0_wp

         do step = 1, N_STEPS
            call ice_thermo_columns(grid%nghost, nx, ny, 1, NK, DTT, .false., wet_mask, &
                                    ice%m_ice, ice%m_snow, ice%enth_ice, ice%enth_snow, &
                                    ice%sal_ice, sf_0, dsf_dt, sw_dn, fprec, tfw_arr, fb_arr, &
                                    sst_arr, s_surf_arr, &
                                    tsurf_out, h2o_ocn_to_ice, h2o_ice_to_ocn, &
                                    heat_to_ocn, sw_thru, snow_to_ice)

            !$acc update self(h2o_ocn_to_ice, h2o_ice_to_ocn)
            sum_ocn_to_ice = sum_ocn_to_ice + h2o_ocn_to_ice(ip, jp, 1)
            sum_ice_to_ocn = sum_ice_to_ocn + h2o_ice_to_ocn(ip, jp, 1)
         end do

         !$acc update self(ice%m_ice)
         m_ice_end = ice%m_ice(ip, jp, 1)

         call check(error, m_ice_end < m_ice0 - 1.0_wp, &
                    "warm forcing must thin the ice cap (net m_ice decrease)")
         if (allocated(error)) exit checks
         call check(error, sum_ice_to_ocn > 0.0_wp, &
                    "surface/basal melt must produce ice->ocean meltwater")
         if (allocated(error)) exit checks

         ! Water-mass identity (same bookkeeping as the Stefan gate).
         mass_identity = sum_ocn_to_ice - sum_ice_to_ocn
         tol = 1.0e-9_wp*max(abs(m_ice_end - m_ice0), 1.0_wp)
         call check(error, abs(mass_identity - (m_ice_end - m_ice0)) <= tol, &
                    "accumulated ocn<->ice water flux must equal the ice mass change")

      end block checks

      !$acc exit data delete(wet_mask, sf_0, dsf_dt, sw_dn, fprec, tfw_arr, fb_arr, &
      !$acc                  sst_arr, s_surf_arr, tsurf_out, h2o_ocn_to_ice, &
      !$acc                  h2o_ice_to_ocn, heat_to_ocn, sw_thru, snow_to_ice)
      call ice%exit_data()
      !$acc exit data delete(ice)
      call ice%destroy()
   end subroutine test_melt_ablation

   subroutine test_conduction_solve_energy_closure(error)
      !! Formerly `column_energy_closure` — renamed (PR-60): the old name
      !! read as "the column conserves energy", a claim this test does
      !! NOT establish. Three verified exclusions:
      !! (a) covers `ice_temp_sis2` (the conduction solve) only, not the
      !!     resize/peel/rebalance chain `ice_column_step` runs after it
      !!     (bottom-freeze, top/bottom melt-peel, rebalance) — this test
      !!     *runs* those branches below only to ADVANCE the trajectory,
      !!     asserting nothing about them;
      !! (b) deliberately excludes the liq-lim clamp
      !!     (`rdb_ice_column.F90:310-316`: `col_enth_out` is measured
      !!     BEFORE it) — SIS2's own `col_check` DOES check that clamp
      !!     (`SIS2_ice_thm.F90:526-537`), this port does not;
      !! (c) the identity omits a `-e_extra_sum` term
      !!     (`rdb_ice_column.F90:494,511,529,544,561` — accumulated,
      !!     never read) and is exact only on a snow-free, non-pinning
      !!     column: EVERY column in this file is seeded `m_snow=0` and
      !!     cold, so the omission never surfaces here.
      !! `column_step_energy_closes_cold` (below) is the end-to-end
      !! answer to "does the column conserve energy" this test's old name
      !! implied.
      !!
      !! Host-side drive of `ice_temp_sis2` directly (the same routine
      !! the device calls). Every hourly step over a 30-day run, form
      !! `residual = (col_enth_out - col_enth_in) - (sum_sol +
      !! tflux_sfc + tflux_bot)` and assert `|residual| <=
      !! 1e-9*(sum of magnitudes)` — SIS2's own tfb_resid_err
      !! normalization convention (sum-of-magnitudes, not the
      !! near-cancelling net heat, which would manufacture fake error
      !! out of O(1e-15) roundoff). Advances resize/rebalance host-side
      !! after each conduction step so the run stays on the production
      !! trajectory.
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: TFW = -1.9_wp
      real(wp), parameter :: TSURF_TARGET = -20.0_wp
      real(wp), parameter :: DTT = 3600.0_wp
      integer, parameter :: N_STEPS = 30*24
      real(wp), parameter :: STIFF = 1.0e6_wp

      real(wp) :: m_ice, m_snow, enth_snow_pt
      real(wp) :: enth_ice_bu(NK), sal_ice_bu(NK)
      real(wp) :: enth_loc(0:NK), sal_loc(NK), sol(0:NK)
      real(wp) :: sf_0, dsf_dt, tsurf, tmelt, bmelt
      real(wp) :: col_enth_in, col_enth_out, sum_sol, tflux_sfc, tflux_bot
      real(wp) :: residual, norm, max_rel_imbalance, rel
      real(wp) :: m_lay(0:NK), enthalpy(0:NK), salin(0:NK)
      real(wp) :: enth_ocean, salin_freeze, mtot_ice
      real(wp) :: h2o_dummy, heat_dummy
      integer :: step, k

      checks: block

         call make_initial_column(1.0_wp, TFW, TSURF_TARGET, .true., &
                                  m_ice, m_snow, enth_ice_bu, enth_snow_pt, sal_ice_bu)

         dsf_dt = STIFF
         sf_0 = -STIFF*TSURF_TARGET
         sol = 0.0_wp

         max_rel_imbalance = 0.0_wp

         do step = 1, N_STEPS
            ! Gather (bottom-up -> top-down, TRAP #2 — mirrors the
            ! ice_column_step gather).
            do k = 1, NK
               enth_loc(k) = enth_ice_bu(NK + 1 - k)
               sal_loc(k) = sal_ice_bu(NK + 1 - k)
            end do
            enth_loc(0) = enth_snow_pt

            tmelt = 0.0_wp
            bmelt = 0.0_wp
            call ice_temp_sis2(NK, m_snow, m_ice, sal_loc, enth_loc, &
                               sf_0, dsf_dt, sol, TFW, 0.0_wp, DTT, &
                               tsurf, tmelt, bmelt, &
                               col_enth_in, col_enth_out, sum_sol, tflux_sfc, tflux_bot)

            residual = (col_enth_out - col_enth_in) - (sum_sol + tflux_sfc + tflux_bot)
            norm = 1.0e-9_wp*(abs(col_enth_out) + abs(col_enth_in) + abs(sum_sol) &
                              + abs(tflux_sfc) + abs(tflux_bot))
            rel = 0.0_wp
            if (norm > 0.0_wp) rel = abs(residual)/norm
            max_rel_imbalance = max(max_rel_imbalance, rel)

            ! Advance resize + rebalance host-side (production trajectory).
            m_lay(0) = m_snow
            m_lay(1:NK) = m_ice/real(NK, wp)
            salin(0) = 0.0_wp
            salin(1:NK) = sal_loc(1:NK)
            enthalpy(0:NK) = enth_loc(0:NK)
            salin_freeze = ICE_BULK_SALINITY
            enth_ocean = ice_enthalpy_liquid(TFW, salin(NK))  ! TRAP #1

            if (tmelt < 0.0_wp) then
               bmelt = bmelt + tmelt
               tmelt = 0.0_wp
            end if

            h2o_dummy = 0.0_wp
            heat_dummy = 0.0_wp
            call ice_bottom_freeze(NK, m_lay, enthalpy, salin, bmelt, &
                                   enth_ocean, salin_freeze, h2o_dummy)
            call ice_top_melt_peel(NK, m_lay, enthalpy, salin, tmelt, &
                                   heat_dummy, h2o_dummy)
            call ice_bottom_melt_peel(NK, m_lay, enthalpy, salin, bmelt, &
                                      heat_dummy, h2o_dummy)
            call ice_rebalance_layers(NK, m_lay, enthalpy, salin, mtot_ice)

            m_ice = mtot_ice
            m_snow = m_lay(0)
            enth_snow_pt = enthalpy(0)
            do k = 1, NK
               enth_ice_bu(NK + 1 - k) = enthalpy(k)
               sal_ice_bu(NK + 1 - k) = salin(k)
            end do
         end do

         call check(error, max_rel_imbalance <= 1.0_wp, &
                    "energy-closure residual must stay within the 1e-9 budget every step")

      end block checks
   end subroutine test_conduction_solve_energy_closure

   subroutine test_equilibrium_fixed_point(error)
      !! `h_eq = ICE_K_ICE*|tsurf-tfw|/fb` with tsurf=-15, tfw=-1.9,
      !! fb=8 => h_eq = 3.324125 m. Seed AT h_eq with the linear
      !! profile; run 10 days at dtt=21600s (40 steps) via
      !! `ice_thermo_columns`; assert |h(end) - h_eq| < 1e-6 m. A
      !! mis-signed/mis-wired `fb` (the failure this test guards) moves
      !! h by ~5.7e-4 m in a single step — the 1e-6 gate sits 3 decades
      !! above the prototype's measured 5.7e-10 m round-off drift and 3
      !! below the failure signal.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_sea_ice_t) :: ice
      real(wp), parameter :: TFW = -1.9_wp
      real(wp), parameter :: TSURF_TARGET = -15.0_wp
      real(wp), parameter :: FB = 8.0_wp
      real(wp), parameter :: H_EQ = ICE_K_ICE*15.0_wp*13.1_wp/(15.0_wp*8.0_wp)
      ! Written to avoid an intrinsic-name collision; equals
      ! ICE_K_ICE*abs(TSURF_TARGET-TFW)/FB = 2.03*13.1/8 = 3.324125.
      real(wp), parameter :: DTT = 21600.0_wp
      integer, parameter :: N_STEPS = int(10.0_wp*DAY_S/DTT)

      real(wp), allocatable :: wet_mask(:, :), sf_0(:, :), dsf_dt(:, :), sw_dn(:, :), fprec(:, :)
      real(wp), allocatable :: tfw_arr(:, :), fb_arr(:, :), sst_arr(:, :), s_surf_arr(:, :)
      real(wp), allocatable :: tsurf_out(:, :, :), h2o_ocn_to_ice(:, :, :), &
                               h2o_ice_to_ocn(:, :, :), heat_to_ocn(:, :, :), sw_thru(:, :, :), &
                               snow_to_ice(:, :, :)
      real(wp) :: m_ice_seed, m_snow_seed, enth_snow_seed
      real(wp) :: enth_ice_seed(NK), sal_ice_seed(NK)
      integer :: ip, jp, nx, ny, step
      real(wp) :: h_end

      checks: block

         call check(error, abs(H_EQ - 3.324125_wp) < 1.0e-9_wp, &
                    "H_EQ literal must equal the analytic 3.324125 m")
         if (allocated(error)) exit checks

         call grid%init(NX_PHYS, NY_PHYS, NGHOST, 1.0_wp, 1.0_wp)
         nx = grid%nx_total
         ny = grid%ny_total
         ice%enable = .true.
         ice%ncat = 1
         ice%nk_ice = NK
         call ice%init(grid)

         allocate (wet_mask(nx, ny), source=1.0_wp)
         allocate (sf_0(nx, ny), dsf_dt(nx, ny), sw_dn(nx, ny), fprec(nx, ny))
         allocate (tfw_arr(nx, ny), fb_arr(nx, ny), sst_arr(nx, ny), s_surf_arr(nx, ny))
         allocate (tsurf_out(nx, ny, 1), h2o_ocn_to_ice(nx, ny, 1), &
                   h2o_ice_to_ocn(nx, ny, 1), heat_to_ocn(nx, ny, 1), sw_thru(nx, ny, 1), &
                   snow_to_ice(nx, ny, 1))

         dsf_dt = 1.0e6_wp
         sf_0 = -1.0e6_wp*TSURF_TARGET
         sw_dn = 0.0_wp
         fprec = 0.0_wp
         tfw_arr = TFW
         fb_arr = FB
         sst_arr = TFW
         s_surf_arr = ICE_BULK_SALINITY

         call make_initial_column(H_EQ, TFW, TSURF_TARGET, .true., &
                                  m_ice_seed, m_snow_seed, enth_ice_seed, &
                                  enth_snow_seed, sal_ice_seed)

         ip = grid%nghost + 2
         jp = grid%nghost + 2

         !$acc enter data copyin(ice)
         call ice%enter_data()
         !$acc enter data copyin(wet_mask, sf_0, dsf_dt, sw_dn, fprec, tfw_arr, fb_arr, &
         !$acc                   sst_arr, s_surf_arr, tsurf_out, h2o_ocn_to_ice, &
         !$acc                   h2o_ice_to_ocn, heat_to_ocn, sw_thru, snow_to_ice)

         ice%m_ice(ip, jp, 1) = m_ice_seed
         ice%m_snow(ip, jp, 1) = m_snow_seed
         ice%enth_ice(ip, jp, 1, :) = enth_ice_seed
         ice%enth_snow(ip, jp, 1, 1) = enth_snow_seed
         ice%sal_ice(ip, jp, 1, :) = sal_ice_seed
         !$acc update device(ice%m_ice, ice%m_snow, ice%enth_ice, ice%enth_snow, ice%sal_ice)

         do step = 1, N_STEPS
            call ice_thermo_columns(grid%nghost, nx, ny, 1, NK, DTT, .false., wet_mask, &
                                    ice%m_ice, ice%m_snow, ice%enth_ice, ice%enth_snow, &
                                    ice%sal_ice, sf_0, dsf_dt, sw_dn, fprec, tfw_arr, fb_arr, &
                                    sst_arr, s_surf_arr, &
                                    tsurf_out, h2o_ocn_to_ice, h2o_ice_to_ocn, &
                                    heat_to_ocn, sw_thru, snow_to_ice)
         end do

         !$acc update self(ice%m_ice)
         h_end = ice%m_ice(ip, jp, 1)/ICE_RHO_ICE
         call check(error, abs(h_end - H_EQ) < 1.0e-6_wp, &
                    "equilibrium thickness must drift by < 1e-6 m over 10 days")

      end block checks

      !$acc exit data delete(wet_mask, sf_0, dsf_dt, sw_dn, fprec, tfw_arr, fb_arr, &
      !$acc                  sst_arr, s_surf_arr, tsurf_out, h2o_ocn_to_ice, &
      !$acc                  h2o_ice_to_ocn, heat_to_ocn, sw_thru, snow_to_ice)
      call ice%exit_data()
      !$acc exit data delete(ice)
      call ice%destroy()
   end subroutine test_equilibrium_fixed_point

   subroutine test_rebalance_conserves(error)
      !! Host-direct `ice_rebalance_layers` on an unequal 2-layer
      !! top-down column: m_lay = [0, 300, 100], distinct enthalpies
      !! and salinities. Assert mtot_ice = 400, m_lay(1) == m_lay(2) ==
      !! 200 to 1e-12 relative, total enthalpy and total salt conserved
      !! to 1e-12 relative, snow slot untouched.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: m_lay(0:NK), enthalpy(0:NK), salin(0:NK)
      real(wp) :: enth_before, salt_before, enth_after, salt_after, mtot_ice
      real(wp) :: snow_enth_before

      checks: block

         m_lay = [0.0_wp, 300.0_wp, 100.0_wp]
         enthalpy(0) = ice_enth_from_ts(0.0_wp, 0.0_wp)
         enthalpy(1) = ice_enth_from_ts(-5.0_wp, ICE_BULK_SALINITY)
         enthalpy(2) = ice_enth_from_ts(-1.0_wp, ICE_BULK_SALINITY)
         salin = [0.0_wp, 4.0_wp, 6.0_wp]
         snow_enth_before = enthalpy(0)

         enth_before = m_lay(1)*enthalpy(1) + m_lay(2)*enthalpy(2)
         salt_before = m_lay(1)*salin(1) + m_lay(2)*salin(2)

         call ice_rebalance_layers(NK, m_lay, enthalpy, salin, mtot_ice)

         call check(error, abs(mtot_ice - 400.0_wp) < 1.0e-9_wp, "mtot_ice must equal 400")
         if (allocated(error)) exit checks
         call check(error, abs(m_lay(1) - 200.0_wp) < 1.0e-12_wp*200.0_wp, &
                    "m_lay(1) must equal mtot/2")
         if (allocated(error)) exit checks
         call check(error, abs(m_lay(2) - 200.0_wp) < 1.0e-12_wp*200.0_wp, &
                    "m_lay(2) must equal mtot/2")
         if (allocated(error)) exit checks

         enth_after = m_lay(1)*enthalpy(1) + m_lay(2)*enthalpy(2)
         salt_after = m_lay(1)*salin(1) + m_lay(2)*salin(2)
         call check(error, abs(enth_after - enth_before) < 1.0e-12_wp*abs(enth_before), &
                    "total enthalpy must be conserved")
         if (allocated(error)) exit checks
         call check(error, abs(salt_after - salt_before) < 1.0e-12_wp*abs(salt_before), &
                    "total salt must be conserved")
         if (allocated(error)) exit checks
         call check(error, abs(enthalpy(0) - snow_enth_before) < 1.0e-14_wp, &
                    "snow slot must be untouched")

      end block checks
   end subroutine test_rebalance_conserves

   subroutine test_optics_csim4(error)
      !! `ice_optics_csim4` closed-form checks (SIS_optics.F90:371-409):
      !! bare thick cold ice, snow-covered, melting bare ice; each case
      !! asserts the shortwave-partition identity `abs_sfc + abs_snow +
      !! sum(abs_ice_lay) + abs_ocn == 1` to 1e-12.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: albedo, abs_sfc, abs_snow, abs_ocn, abs_int, pen
      real(wp) :: abs_ice_lay(NK)
      real(wp) :: fh, expect_albedo, snow_cover, t_fr, ai

      checks: block

         ! ---- Case 1: bare thick cold ice ----
         call ice_optics_csim4(NK, 0.0_wp, 2.0_wp, -10.0_wp, ICE_BULK_SALINITY, &
                               albedo, abs_sfc, abs_snow, abs_ice_lay, abs_ocn, abs_int, pen)
         fh = min(atan(5.0_wp*2.0_wp)/atan(5.0_wp*0.5_wp), 1.0_wp)
         t_fr = ice_t_freeze(ICE_BULK_SALINITY)
         call check(error, -10.0_wp + ICE_T_RANGE_MELT <= t_fr, &
                    "case-1 setup must be below the melt-reduction threshold")
         if (allocated(error)) exit checks
         expect_albedo = fh*ICE_ALB_ICE + (1.0_wp - fh)*0.06_wp
         call check(error, abs(albedo - expect_albedo) < 1.0e-12_wp, &
                    "bare cold ice albedo must match the closed form")
         if (allocated(error)) exit checks
         call check(error, abs(pen - ICE_PEN_ICE) < 1.0e-12_wp, "pen must equal ICE_PEN_ICE")
         if (allocated(error)) exit checks
         call check(error, abs(abs_sfc - (1.0_wp - ICE_PEN_ICE)) < 1.0e-12_wp, &
                    "abs_sfc must equal 1-pen")
         if (allocated(error)) exit checks
         call check(error, abs(abs_snow) < 1.0e-14_wp, "abs_snow must be 0 (CSIM4 branch)")
         if (allocated(error)) exit checks
         call check(error, abs(abs_sfc + abs_snow + sum(abs_ice_lay) + abs_ocn - 1.0_wp) &
                    < 1.0e-12_wp, "case-1 SW partition must sum to 1")
         if (allocated(error)) exit checks

         ! ---- Case 2: snow-covered ----
         call ice_optics_csim4(NK, 0.5_wp, 2.0_wp, -10.0_wp, ICE_BULK_SALINITY, &
                               albedo, abs_sfc, abs_snow, abs_ice_lay, abs_ocn, abs_int, pen)
         snow_cover = 0.5_wp/(0.5_wp + ICE_SNOW_PATCH)
         call check(error, abs(pen - (1.0_wp - snow_cover)*ICE_PEN_ICE) < 1.0e-12_wp, &
                    "snow-covered pen must equal (1-snow_cover)*ICE_PEN_ICE")
         if (allocated(error)) exit checks
         ! Closed form: snow_cover*alb_snow + (1-snow_cover)*bare-cold-ice
         ! albedo (ts=-10 is below the melt-reduction threshold, so the
         ! ice component equals case-1's expect_albedo). Pinning the exact
         ! blend rejects a swapped snow_cover*ai + (1-snow_cover)*as fill.
         call check(error, abs(albedo - (snow_cover*ICE_ALB_SNOW &
                                         + (1.0_wp - snow_cover)*expect_albedo)) < 1.0e-12_wp, &
                    "snow-covered albedo must match the exact snow_cover blend")
         if (allocated(error)) exit checks
         call check(error, abs(abs_sfc + abs_snow + sum(abs_ice_lay) + abs_ocn - 1.0_wp) &
                    < 1.0e-12_wp, "case-2 SW partition must sum to 1")
         if (allocated(error)) exit checks

         ! ---- Case 3: melting bare ice (ts == t_fr, ramp saturated) ----
         call ice_optics_csim4(NK, 0.0_wp, 2.0_wp, t_fr, ICE_BULK_SALINITY, &
                               albedo, abs_sfc, abs_snow, abs_ice_lay, abs_ocn, abs_int, pen)
         ai = ICE_ALB_ICE - 0.075_wp
         ai = fh*ai + (1.0_wp - fh)*0.06_wp
         call check(error, abs(albedo - ai) < 1.0e-12_wp, &
                    "melting bare ice albedo must match the closed form")
         if (allocated(error)) exit checks
         call check(error, abs(abs_sfc + abs_snow + sum(abs_ice_lay) + abs_ocn - 1.0_wp) &
                    < 1.0e-12_wp, "case-3 SW partition must sum to 1")

      end block checks
   end subroutine test_optics_csim4

   subroutine test_laytemp_quadratic(error)
      !! `laytemp_sis2`: for a grid of salty-quadratic inputs, the
      !! result T satisfies the per-layer implicit balance
      !! `m*(E(T)-E(tp)) - dtt*(qf-bf*T) == 0` (mushy energy E(x) =
      !! Cp_ice*(x-t_fr) - Lat_fus*(1-t_fr/x)) to <= 1e-9 relative to
      !! m*Lat_fus, skipping the pinned-to-freezing branch (identity
      !! does not apply there).
      !!
      !! `update_lay_enth`: for a small case matrix including an
      !! m_lay=0 diagnostic branch and a pin-to-max branch, the
      !! conservation identity `m_lay*(enth_new-enth_in) + extra_heat
      !! == dtt*((ht_body+ftop_new)-fbot_new)` holds to <= 1e-9
      !! relative, across all three flux-bookkeeping branches (driven
      !! via small vs. large dftop_dt/dfbot_dt magnitudes).
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: M = 452.5_wp
      real(wp), parameter :: DTT = 3600.0_wp
      real(wp) :: t_fr, tp, qf, bf, t_new, e_t, e_tp, res, tol
      integer :: i_tp, i_qf, i_bf
      real(wp), parameter :: TP_VALS(2) = [-8.0_wp, -3.0_wp]
      real(wp), parameter :: QF_VALS(2) = [-50.0_wp, 20.0_wp]
      real(wp), parameter :: BF_VALS(2) = [5.0_wp, 25.0_wp]
      real(wp) :: hf_err_rat

      checks: block

         t_fr = ice_t_freeze(ICE_BULK_SALINITY)

         do i_tp = 1, 2
            do i_qf = 1, 2
               do i_bf = 1, 2
                  tp = TP_VALS(i_tp)
                  qf = QF_VALS(i_qf)
                  bf = BF_VALS(i_bf)
                  t_new = laytemp_sis2(M, t_fr, qf, bf, tp, DTT)
                  if (t_new >= t_fr - 1.0e-12_wp) cycle  ! pinned branch, skip
                  e_t = ICE_CP_ICE*(t_new - t_fr) - ICE_LAT_FUS*(1.0_wp - t_fr/t_new)
                  e_tp = ICE_CP_ICE*(tp - t_fr) - ICE_LAT_FUS*(1.0_wp - t_fr/tp)
                  res = M*(e_t - e_tp) - DTT*(qf - bf*t_new)
                  tol = 1.0e-9_wp*M*ICE_LAT_FUS
                  call check(error, abs(res) <= tol, &
                             "laytemp_sis2 must satisfy its own implicit energy balance")
                  if (allocated(error)) exit checks
               end do
            end do
         end do

         ! ---- update_lay_enth conservation identity, 4 branches ----
         hf_err_rat = 0.7071_wp*DTT*ICE_TEMP_RANGE_EST/ &
                      (ICE_TEMP_RANGE_EST*ICE_CP_ICE + ICE_LAT_FUS)

         ! (a) both-explicit
         call run_update_lay_enth_case(M, ICE_BULK_SALINITY, &
                                       ice_enth_from_ts(-5.0_wp, ICE_BULK_SALINITY), &
                                       30.0_wp, 5.0_wp, -20.0_wp, -0.5_wp, 0.5_wp, DTT, &
                                       hf_err_rat, .false., 0.0_wp, error)
         if (allocated(error)) exit checks

         ! (b) top-inverted (large |dftop_dt|)
         call run_update_lay_enth_case(M, ICE_BULK_SALINITY, &
                                       ice_enth_from_ts(-5.0_wp, ICE_BULK_SALINITY), &
                                       30.0_wp, 5.0_wp, -20.0_wp, -1.0e6_wp, 0.5_wp, DTT, &
                                       hf_err_rat, .false., 0.0_wp, error)
         if (allocated(error)) exit checks

         ! (c) massless diagnostic branch
         call run_update_lay_enth_case(0.0_wp, ICE_BULK_SALINITY, &
                                       ice_enth_from_ts(-5.0_wp, ICE_BULK_SALINITY), &
                                       30.0_wp, 5.0_wp, -20.0_wp, -0.5_wp, 0.5_wp, DTT, &
                                       hf_err_rat, .false., 0.0_wp, error)
         if (allocated(error)) exit checks

         ! (d) pin-to-max branch (huge body heating pins to max_temp = t_fr)
         call run_update_lay_enth_case(M, ICE_BULK_SALINITY, &
                                       ice_enth_from_ts(-0.5_wp, ICE_BULK_SALINITY), &
                                       3000.0_wp, 500.0_wp, -2000.0_wp, -0.5_wp, 0.5_wp, DTT, &
                                       hf_err_rat, .false., 0.0_wp, error)

      end block checks
   end subroutine test_laytemp_quadratic

   subroutine run_update_lay_enth_case(m_lay, sice, enth_in, ftop, ht_body, fbot, &
                                       dftop_dt, dfbot_dt, dtt, hf_err_rat, &
                                       has_temp_max, temp_max, error)
      !! Drive `update_lay_enth` once and assert its own conservation
      !! identity `m_lay*(enth_new-enth_in) + extra_heat ==
      !! dtt*((ht_body+ftop_new)-fbot_new)` to <= 1e-9 relative.
      real(wp), intent(in) :: m_lay, sice, enth_in, ftop, ht_body, fbot
      real(wp), intent(in) :: dftop_dt, dfbot_dt, dtt, hf_err_rat, temp_max
      logical, intent(in) :: has_temp_max
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: enth, ftop_io, fbot_io, extra_heat, new_temp, lhs, rhs, mag_sum

      enth = enth_in
      ftop_io = ftop
      fbot_io = fbot
      call update_lay_enth(m_lay, sice, enth, ftop_io, ht_body, fbot_io, &
                           dftop_dt, dfbot_dt, dtt, hf_err_rat, &
                           extra_heat, new_temp, has_temp_max, temp_max)

      lhs = m_lay*(enth - enth_in) + extra_heat
      rhs = dtt*((ht_body + ftop_io) - fbot_io)
      mag_sum = abs(lhs) + abs(rhs) + 1.0_wp
      call check(error, abs(lhs - rhs) <= 1.0e-9_wp*mag_sum, &
                 "update_lay_enth must satisfy its own conservation identity")
   end subroutine run_update_lay_enth_case

   ! -----------------------------------------------------------------
   ! PR 27: Archimedes freeboard snow-ice flooding
   ! -----------------------------------------------------------------

   subroutine test_snow_ice_archimedes(error)
      !! PR 27 physics gate (PLAN_PR27_snow_ice_flooding.md §9.1). The
      !! plan's worked golden: h_i=1.0 m, h_s=0.5 m, nk=2 => m_lay =
      !! [165.0, 452.5, 452.5]. Distinct enthalpies (snow -10 degC, ice
      !! layers -5/-1 degC), salin = [0, 4, 4]. Asserts snow_to_ice
      !! equals the analytic 35.145631... kg/m^2 EXACTLY (to 1e-9
      !! relative), the post-flood column sits EXACTLY at the
      !! Archimedes waterline (to 1e-12 relative -- the property that
      !! proves the physics, not just "some snow converts"), mass/
      !! enthalpy/salt conservation to round-off, and that the
      !! snow's cold/fresh mass landed in the TOP ice layer (local
      !! index 1 -- TRAP #2 in rdb_ice_column's module docstring) while
      !! layer 2 is untouched.
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: RHO_RATIO = ICE_RHO_ICE/ICE_RHO_OCEAN
      real(wp), parameter :: SNOW_TO_ICE_GOLDEN = 35.145631067961_wp
      real(wp) :: m_lay(0:NK), enthalpy(0:NK), salin(0:NK)
      real(wp) :: m_lay0(0:NK), enthalpy0(0:NK), salin0(0:NK)
      real(wp) :: snow_to_ice, m_i_before, m_s_before, m_i_after, m_s_after
      real(wp) :: enth_before, enth_after, salt_before, salt_after, equilibrium_resid

      checks: block

         m_lay = [165.0_wp, 452.5_wp, 452.5_wp]
         enthalpy(0) = ice_enth_from_ts(-10.0_wp, 0.0_wp)
         enthalpy(1) = ice_enth_from_ts(-5.0_wp, 4.0_wp)
         enthalpy(2) = ice_enth_from_ts(-1.0_wp, 4.0_wp)
         salin = [0.0_wp, 4.0_wp, 4.0_wp]
         m_lay0 = m_lay
         enthalpy0 = enthalpy
         salin0 = salin

         m_i_before = m_lay(1) + m_lay(2)
         m_s_before = m_lay(0)
         enth_before = sum(m_lay(0:NK)*enthalpy(0:NK))
         salt_before = m_lay(1)*salin(1) + m_lay(2)*salin(2)

         call ice_snow_ice_flood(NK, m_lay, enthalpy, salin, RHO_RATIO, snow_to_ice)

         call check(error, abs(snow_to_ice - SNOW_TO_ICE_GOLDEN) < 1.0e-9_wp*SNOW_TO_ICE_GOLDEN, &
                    "snow_to_ice must equal the analytic Archimedes worked golden")
         if (allocated(error)) exit checks

         m_i_after = m_lay(1) + m_lay(2)
         m_s_after = m_lay(0)

         equilibrium_resid = (m_i_after + m_s_after)*RHO_RATIO - m_i_after
         call check(error, abs(equilibrium_resid) < 1.0e-12_wp*m_i_after, &
                    "post-flood snow-ice interface must sit exactly on the waterline")
         if (allocated(error)) exit checks

         call check(error, abs((m_s_after + m_i_after) - (m_s_before + m_i_before)) < &
                    1.0e-12_wp*(m_s_before + m_i_before), "total mass must be conserved")
         if (allocated(error)) exit checks

         enth_after = sum(m_lay(0:NK)*enthalpy(0:NK))
         call check(error, abs(enth_after - enth_before) < 1.0e-12_wp*abs(enth_before), &
                    "total column enthalpy must be conserved")
         if (allocated(error)) exit checks

         salt_after = m_lay(1)*salin(1) + m_lay(2)*salin(2)
         call check(error, abs(salt_after - salt_before) < 1.0e-12_wp*abs(salt_before), &
                    "total salt must be conserved")
         if (allocated(error)) exit checks

         call check(error, enthalpy(1) < enthalpy0(1), &
                    "top ice layer enthalpy must strictly decrease (snow's cold mass mixed in)")
         if (allocated(error)) exit checks
         call check(error, salin(1) < salin0(1), &
                    "top ice layer salinity must strictly decrease (diluted by zero-salinity snow)")
         if (allocated(error)) exit checks

         call check(error, enthalpy(2) == enthalpy0(2), "layer 2 enthalpy must be untouched")
         if (allocated(error)) exit checks
         call check(error, salin(2) == salin0(2), "layer 2 salinity must be untouched")
         if (allocated(error)) exit checks
         call check(error, m_lay(2) == m_lay0(2), "layer 2 mass must be untouched")

      end block checks
   end subroutine test_snow_ice_archimedes

   subroutine test_snow_ice_threshold_and_clamp(error)
      !! PR 27 branch gate (plan §9.2). Three sub-cases on
      !! `ice_snow_ice_flood` directly.
      !! (a) Below threshold (h_i=1.0, h_s=0.30 < 0.37879 m): asserts
      !!     `snow_to_ice == 0.0` EXACTLY and the column bit-unchanged.
      !! (b) At threshold (h_s = h_i*(rho_w-rho_i)/rho_s exactly, so
      !!     m_s = 125.0 kg/m^2 when h_i=1.0): asserts
      !!     `|snow_to_ice| < 1e-9` -- already exactly at the waterline.
      !! (c) Deep snow, a ~20-pair sweep over h_i in [0.01,5], h_s in
      !!     [0,3]: asserts the SIS2 `min(...)` clamp is PROVABLY DEAD
      !!     CODE in this pond-free model (`snow_to_ice < m_lay(0)`
      !!     strictly whenever snow is present, so a genuine flood never
      !!     consumes all the snow), and that every flooded column hits
      !!     the exact Archimedes equilibrium.
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: RHO_RATIO = ICE_RHO_ICE/ICE_RHO_OCEAN
      integer, parameter :: N_HI = 5, N_HS = 4
      real(wp) :: m_lay(0:NK), enthalpy(0:NK), salin(0:NK)
      real(wp) :: m_lay0(0:NK), enthalpy0(0:NK), salin0(0:NK)
      real(wp) :: snow_to_ice, h_i, h_s, m_i, m_s, resid
      integer :: ii, jj

      checks: block

         ! ---- (a) below threshold ----
         m_lay = [0.30_wp*ICE_RHO_SNOW, 452.5_wp, 452.5_wp]
         enthalpy(0) = ice_enth_from_ts(-10.0_wp, 0.0_wp)
         enthalpy(1) = ice_enth_from_ts(-5.0_wp, 4.0_wp)
         enthalpy(2) = ice_enth_from_ts(-1.0_wp, 4.0_wp)
         salin = [0.0_wp, 4.0_wp, 4.0_wp]
         m_lay0 = m_lay
         enthalpy0 = enthalpy
         salin0 = salin

         call ice_snow_ice_flood(NK, m_lay, enthalpy, salin, RHO_RATIO, snow_to_ice)

         call check(error, snow_to_ice == 0.0_wp, "below threshold: snow_to_ice must be exactly 0")
         if (allocated(error)) exit checks
         call check(error, all(m_lay == m_lay0) .and. all(enthalpy == enthalpy0) .and. &
                    all(salin == salin0), "below threshold: column must be bit-unchanged")
         if (allocated(error)) exit checks

         ! ---- (b) at threshold: m_s = h_i*(rho_w-rho_i) exactly (h_i=1.0) ----
         m_lay = [ICE_RHO_OCEAN - ICE_RHO_ICE, 452.5_wp, 452.5_wp]
         call ice_snow_ice_flood(NK, m_lay, enthalpy, salin, RHO_RATIO, snow_to_ice)
         call check(error, abs(snow_to_ice) < 1.0e-9_wp, &
                    "at threshold: column already at the waterline, snow_to_ice ~ 0")
         if (allocated(error)) exit checks

         ! ---- (c) deep snow: the min(...) clamp must never bind ----
         do ii = 1, N_HI
            h_i = 0.01_wp + (real(ii, wp) - 1.0_wp)*(5.0_wp - 0.01_wp)/real(N_HI - 1, wp)
            do jj = 1, N_HS
               h_s = (real(jj, wp) - 1.0_wp)*3.0_wp/real(N_HS - 1, wp)
               m_i = h_i*ICE_RHO_ICE
               m_s = h_s*ICE_RHO_SNOW
               m_lay(0) = m_s
               m_lay(1) = 0.5_wp*m_i
               m_lay(2) = 0.5_wp*m_i
               enthalpy(0) = ice_enth_from_ts(-10.0_wp, 0.0_wp)
               enthalpy(1) = ice_enth_from_ts(-5.0_wp, 4.0_wp)
               enthalpy(2) = ice_enth_from_ts(-1.0_wp, 4.0_wp)
               salin = [0.0_wp, 4.0_wp, 4.0_wp]

               call ice_snow_ice_flood(NK, m_lay, enthalpy, salin, RHO_RATIO, snow_to_ice)

               if (m_s > 0.0_wp) then
                  call check(error, snow_to_ice < m_s, &
                             "the SIS2 min(...) clamp must never bind in a pond-free model")
                  if (allocated(error)) exit checks
               end if
               if (snow_to_ice > 0.0_wp) then
                  call check(error, m_lay(0) > 0.0_wp, &
                             "a genuine flood must never consume all the snow")
                  if (allocated(error)) exit checks
                  resid = (m_lay(1) + m_lay(2) + m_lay(0))*RHO_RATIO - (m_lay(1) + m_lay(2))
                  call check(error, abs(resid) < 1.0e-10_wp*max(m_lay(1) + m_lay(2), 1.0_wp), &
                             "flooded columns must sit exactly at Archimedes equilibrium")
                  if (allocated(error)) exit checks
               end if
            end do
         end do

      end block checks
   end subroutine test_snow_ice_threshold_and_clamp

   subroutine test_snow_ice_column_device(error)
      !! PR 27 wiring gate (plan §9.3). Full `ice_thermo_columns` device
      !! driver on ONE physical cell/category, snow-carrying (m_ice=905,
      !! m_snow=165 kg/m^2 -- the §9.1 golden), forcing held at a
      !! GENUINE numerical fixed point: isothermal AT tfw across snow
      !! AND ice, stiff SEB pin AT tfw, sw_dn=0, fb=0 -- the same
      !! established idiom `test_snow_arrives_at_layer_enthalpy` (part
      !! b, `test_ocean_ice_snowfall.F90`) uses to guarantee conduction
      !! moves no mass this step, so only the Archimedes flood acts.
      !! `do_snow_ice=.true.`.
      !!
      !! Asserts `snow_to_ice(ip,jp,1) > 0` -- the `mem:separate`
      !! canary (CLAUDE.md): if `snow_to_ice` were not mapped in
      !! `enter_data`/`exit_data`, the device write would land in
      !! unmapped memory and this would read a stale host 0.0 with no
      !! crash -- meaningful ONLY on the GPU build, per CLAUDE.md's own
      !! warning that a green multicore run proves nothing here. Also
      !! asserts m_ice+m_snow conservation, the mass moving the right
      !! direction, column enthalpy conservation (looser tolerance --
      !! the conduction solve legitimately runs), and that every OTHER
      !! cell (including the ghost band) is exactly untouched.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_sea_ice_t) :: ice
      real(wp), parameter :: TFW = -1.9_wp
      real(wp), parameter :: DTT = 3600.0_wp
      real(wp), parameter :: M_ICE0 = ICE_RHO_ICE*1.0_wp
      real(wp), parameter :: M_SNOW0 = ICE_RHO_SNOW*0.5_wp

      real(wp), allocatable :: wet_mask(:, :), sf_0(:, :), dsf_dt(:, :), sw_dn(:, :), fprec(:, :)
      real(wp), allocatable :: tfw_arr(:, :), fb_arr(:, :), sst_arr(:, :), s_surf_arr(:, :)
      real(wp), allocatable :: tsurf_out(:, :, :), h2o_ocn_to_ice(:, :, :), &
                               h2o_ice_to_ocn(:, :, :), heat_to_ocn(:, :, :), sw_thru(:, :, :), &
                               snow_to_ice(:, :, :)
      real(wp) :: enth_ice_seed(NK), sal_ice_seed(NK), enth_snow_seed
      real(wp) :: m_ice_before, m_snow_before, col_enth_before, col_enth_after
      logical :: only_seeded_nonzero
      integer :: ip, jp, nx, ny, i, j

      checks: block

         call grid%init(NX_PHYS, NY_PHYS, NGHOST, 1.0_wp, 1.0_wp)
         nx = grid%nx_total
         ny = grid%ny_total
         ice%enable = .true.
         ice%ncat = 1
         ice%nk_ice = NK
         call ice%init(grid)

         allocate (wet_mask(nx, ny), source=1.0_wp)
         allocate (sf_0(nx, ny), dsf_dt(nx, ny), sw_dn(nx, ny), fprec(nx, ny))
         allocate (tfw_arr(nx, ny), fb_arr(nx, ny), sst_arr(nx, ny), s_surf_arr(nx, ny))
         allocate (tsurf_out(nx, ny, 1), h2o_ocn_to_ice(nx, ny, 1), &
                   h2o_ice_to_ocn(nx, ny, 1), heat_to_ocn(nx, ny, 1), sw_thru(nx, ny, 1), &
                   snow_to_ice(nx, ny, 1))

         ! Stiff SEB pin AT tfw, combined with an isothermal-at-tfw seed
         ! (snow AND ice): a true fixed point, conduction moves nothing.
         dsf_dt = 1.0e6_wp
         sf_0 = -1.0e6_wp*TFW
         sw_dn = 0.0_wp
         fprec = 0.0_wp
         tfw_arr = TFW
         fb_arr = 0.0_wp
         sst_arr = TFW
         s_surf_arr = ICE_BULK_SALINITY

         enth_ice_seed(:) = ice_enth_from_ts(TFW, ICE_BULK_SALINITY)
         sal_ice_seed(:) = ICE_BULK_SALINITY
         enth_snow_seed = ice_enth_from_ts(TFW, 0.0_wp)

         ip = grid%nghost + 2
         jp = grid%nghost + 2

         !$acc enter data copyin(ice)
         call ice%enter_data()
         !$acc enter data copyin(wet_mask, sf_0, dsf_dt, sw_dn, fprec, tfw_arr, fb_arr, &
         !$acc                   sst_arr, s_surf_arr, tsurf_out, h2o_ocn_to_ice, &
         !$acc                   h2o_ice_to_ocn, heat_to_ocn, sw_thru, snow_to_ice)

         ice%m_ice(ip, jp, 1) = M_ICE0
         ice%m_snow(ip, jp, 1) = M_SNOW0
         ice%enth_ice(ip, jp, 1, :) = enth_ice_seed
         ice%enth_snow(ip, jp, 1, 1) = enth_snow_seed
         ice%sal_ice(ip, jp, 1, :) = sal_ice_seed
         !$acc update device(ice%m_ice, ice%m_snow, ice%enth_ice, ice%enth_snow, ice%sal_ice)

         m_ice_before = M_ICE0
         m_snow_before = M_SNOW0
         col_enth_before = M_SNOW0*enth_snow_seed + sum(enth_ice_seed(:))*M_ICE0/real(NK, wp)

         call ice_thermo_columns(grid%nghost, nx, ny, 1, NK, DTT, .true., wet_mask, &
                                 ice%m_ice, ice%m_snow, ice%enth_ice, ice%enth_snow, &
                                 ice%sal_ice, sf_0, dsf_dt, sw_dn, fprec, tfw_arr, fb_arr, &
                                 sst_arr, s_surf_arr, &
                                 tsurf_out, h2o_ocn_to_ice, h2o_ice_to_ocn, &
                                 heat_to_ocn, sw_thru, snow_to_ice)

         !$acc update self(ice%m_ice, ice%m_snow, ice%enth_ice, ice%enth_snow, snow_to_ice)

         call check(error, snow_to_ice(ip, jp, 1) > 0.0_wp, &
                    "snow_to_ice must be nonzero on the flooded cell -- the mem:separate canary")
         if (allocated(error)) exit checks

         call check(error, abs((ice%m_ice(ip, jp, 1) + ice%m_snow(ip, jp, 1)) - &
                               (m_ice_before + m_snow_before)) < &
                    1.0e-12_wp*(m_ice_before + m_snow_before), &
                    "m_ice + m_snow must be conserved across the step")
         if (allocated(error)) exit checks

         call check(error, ice%m_snow(ip, jp, 1) < m_snow_before, "m_snow must decrease")
         if (allocated(error)) exit checks
         call check(error, ice%m_ice(ip, jp, 1) > m_ice_before, "m_ice must increase")
         if (allocated(error)) exit checks

         col_enth_after = ice%m_snow(ip, jp, 1)*ice%enth_snow(ip, jp, 1, 1) + &
                          sum(ice%enth_ice(ip, jp, 1, :))*ice%m_ice(ip, jp, 1)/real(NK, wp)
         call check(error, abs(col_enth_after - col_enth_before) < 1.0e-9_wp*abs(col_enth_before), &
                    "column enthalpy (snow + nk equal-mass ice layers) must be conserved")
         if (allocated(error)) exit checks

         only_seeded_nonzero = .true.
         do j = 1, ny
            do i = 1, nx
               if (i == ip .and. j == jp) cycle
               if (snow_to_ice(i, j, 1) /= 0.0_wp) only_seeded_nonzero = .false.
            end do
         end do
         call check(error, only_seeded_nonzero, &
                    "snow_to_ice must be exactly 0 everywhere except the seeded cell "// &
                    "(ghost band + every other physical cell)")

      end block checks

      !$acc exit data delete(wet_mask, sf_0, dsf_dt, sw_dn, fprec, tfw_arr, fb_arr, &
      !$acc                  sst_arr, s_surf_arr, tsurf_out, h2o_ocn_to_ice, &
      !$acc                  h2o_ice_to_ocn, heat_to_ocn, sw_thru, snow_to_ice)
      call ice%exit_data()
      !$acc exit data delete(ice)
      call ice%destroy()
   end subroutine test_snow_ice_column_device

   subroutine test_snow_ice_disabled_bitident(error)
      !! PR 27 default-off gate (plan §9.4) -- MANDATORY. A snow-carrying
      !! column (m_ice=905, m_snow=165 kg/m^2, well past the flood
      !! threshold, so a wrongly-wired gate would show up immediately)
      !! run for 12 hourly steps through `ice_thermo_columns` TWICE from
      !! the identical seed and forcing, both times with
      !! `do_snow_ice=.false.`. Asserts (a) `snow_to_ice` is EXACTLY
      !! 0.0 on the seeded cell every step in BOTH runs (the gate
      !! actually gates -- with snow present and well past threshold, a
      !! broken gate would flood and this would catch it), and (b) the
      !! two independent runs are BIT-FOR-BIT identical in every
      !! prognostic (m_ice, m_snow, enth_ice, enth_snow, sal_ice) --
      !! `do_snow_ice=.false.` is a reproducible, deterministic no-op.
      !! The companion half of the acceptance criterion -- the EXISTING
      !! `test_stefan_growth`/`test_melt_ablation` passing at their
      !! ORIGINAL tolerances with the new `do_snow_ice`/`snow_to_ice`
      !! arguments threaded through, no tolerance edited -- is exercised
      !! by those two tests directly (see this file's other cases).
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: m_ice_a, m_snow_a, enth_snow_a
      real(wp) :: enth_ice_a(NK), sal_ice_a(NK)
      real(wp) :: m_ice_b, m_snow_b, enth_snow_b
      real(wp) :: enth_ice_b(NK), sal_ice_b(NK)
      logical :: zero_a, zero_b

      checks: block

         call run_snowy_trajectory_no_flood(m_ice_a, m_snow_a, enth_snow_a, enth_ice_a, &
                                            sal_ice_a, zero_a)
         call run_snowy_trajectory_no_flood(m_ice_b, m_snow_b, enth_snow_b, enth_ice_b, &
                                            sal_ice_b, zero_b)

         call check(error, zero_a, "snow_to_ice must be exactly 0 every step (run A)")
         if (allocated(error)) exit checks
         call check(error, zero_b, "snow_to_ice must be exactly 0 every step (run B)")
         if (allocated(error)) exit checks

         call check(error, m_ice_a == m_ice_b, &
                    "m_ice must be bit-identical across two do_snow_ice=.false. runs")
         if (allocated(error)) exit checks
         call check(error, m_snow_a == m_snow_b, "m_snow must be bit-identical")
         if (allocated(error)) exit checks
         call check(error, enth_snow_a == enth_snow_b, "enth_snow must be bit-identical")
         if (allocated(error)) exit checks
         call check(error, all(enth_ice_a == enth_ice_b), "enth_ice must be bit-identical")
         if (allocated(error)) exit checks
         call check(error, all(sal_ice_a == sal_ice_b), "sal_ice must be bit-identical")

      end block checks
   end subroutine test_snow_ice_disabled_bitident

   subroutine run_snowy_trajectory_no_flood(m_ice_out, m_snow_out, enth_snow_out, &
                                            enth_ice_out, sal_ice_out, snow_to_ice_all_zero)
      !! Host helper for `test_snow_ice_disabled_bitident` (not a
      !! test-drive case itself): seeds a snow-carrying column (same
      !! §9.1 golden masses) under `test_melt_ablation`'s warm-forcing
      !! recipe (TSURF_TARGET=+2 degC, FB=5 W/m^2 -- exercises real
      !! melt/freeze, not a fixed point) and steps it 12 hours through
      !! `ice_thermo_columns` with `do_snow_ice=.false.`, checking
      !! `snow_to_ice` stays exactly 0 every step.
      real(wp), intent(out) :: m_ice_out, m_snow_out, enth_snow_out
      real(wp), intent(out) :: enth_ice_out(NK), sal_ice_out(NK)
      logical, intent(out) :: snow_to_ice_all_zero

      type(hgrid_t) :: grid
      type(ocean_sea_ice_t) :: ice
      real(wp), parameter :: TFW = -1.9_wp
      real(wp), parameter :: TSURF_TARGET = 2.0_wp
      real(wp), parameter :: FB = 5.0_wp
      real(wp), parameter :: DTT = 3600.0_wp
      real(wp), parameter :: M_ICE0 = ICE_RHO_ICE*1.0_wp
      real(wp), parameter :: M_SNOW0 = ICE_RHO_SNOW*0.5_wp
      integer, parameter :: N_STEPS = 12

      real(wp), allocatable :: wet_mask(:, :), sf_0(:, :), dsf_dt(:, :), sw_dn(:, :), fprec(:, :)
      real(wp), allocatable :: tfw_arr(:, :), fb_arr(:, :), sst_arr(:, :), s_surf_arr(:, :)
      real(wp), allocatable :: tsurf_out(:, :, :), h2o_ocn_to_ice(:, :, :), &
                               h2o_ice_to_ocn(:, :, :), heat_to_ocn(:, :, :), sw_thru(:, :, :), &
                               snow_to_ice(:, :, :)
      integer :: ip, jp, nx, ny, step

      call grid%init(NX_PHYS, NY_PHYS, NGHOST, 1.0_wp, 1.0_wp)
      nx = grid%nx_total
      ny = grid%ny_total
      ice%enable = .true.
      ice%ncat = 1
      ice%nk_ice = NK
      call ice%init(grid)

      allocate (wet_mask(nx, ny), source=1.0_wp)
      allocate (sf_0(nx, ny), dsf_dt(nx, ny), sw_dn(nx, ny), fprec(nx, ny))
      allocate (tfw_arr(nx, ny), fb_arr(nx, ny), sst_arr(nx, ny), s_surf_arr(nx, ny))
      allocate (tsurf_out(nx, ny, 1), h2o_ocn_to_ice(nx, ny, 1), &
                h2o_ice_to_ocn(nx, ny, 1), heat_to_ocn(nx, ny, 1), sw_thru(nx, ny, 1), &
                snow_to_ice(nx, ny, 1))

      dsf_dt = 1.0e6_wp
      sf_0 = -1.0e6_wp*TSURF_TARGET
      sw_dn = 0.0_wp
      fprec = 0.0_wp
      tfw_arr = TFW
      fb_arr = FB
      sst_arr = TFW
      s_surf_arr = ICE_BULK_SALINITY

      ip = grid%nghost + 2
      jp = grid%nghost + 2

      !$acc enter data copyin(ice)
      call ice%enter_data()
      !$acc enter data copyin(wet_mask, sf_0, dsf_dt, sw_dn, fprec, tfw_arr, fb_arr, &
      !$acc                   sst_arr, s_surf_arr, tsurf_out, h2o_ocn_to_ice, &
      !$acc                   h2o_ice_to_ocn, heat_to_ocn, sw_thru, snow_to_ice)

      ice%m_ice(ip, jp, 1) = M_ICE0
      ice%m_snow(ip, jp, 1) = M_SNOW0
      ice%enth_ice(ip, jp, 1, :) = ice_enth_from_ts(TFW, ICE_BULK_SALINITY)
      ice%enth_snow(ip, jp, 1, 1) = ice_enth_from_ts(-10.0_wp, 0.0_wp)
      ice%sal_ice(ip, jp, 1, :) = ICE_BULK_SALINITY
      !$acc update device(ice%m_ice, ice%m_snow, ice%enth_ice, ice%enth_snow, ice%sal_ice)

      snow_to_ice_all_zero = .true.

      do step = 1, N_STEPS
         call ice_thermo_columns(grid%nghost, nx, ny, 1, NK, DTT, .false., wet_mask, &
                                 ice%m_ice, ice%m_snow, ice%enth_ice, ice%enth_snow, &
                                 ice%sal_ice, sf_0, dsf_dt, sw_dn, fprec, tfw_arr, fb_arr, &
                                 sst_arr, s_surf_arr, &
                                 tsurf_out, h2o_ocn_to_ice, h2o_ice_to_ocn, &
                                 heat_to_ocn, sw_thru, snow_to_ice)
         !$acc update self(snow_to_ice)
         if (snow_to_ice(ip, jp, 1) /= 0.0_wp) snow_to_ice_all_zero = .false.
      end do

      !$acc update self(ice%m_ice, ice%m_snow, ice%enth_ice, ice%enth_snow, ice%sal_ice)
      m_ice_out = ice%m_ice(ip, jp, 1)
      m_snow_out = ice%m_snow(ip, jp, 1)
      enth_snow_out = ice%enth_snow(ip, jp, 1, 1)
      enth_ice_out(:) = ice%enth_ice(ip, jp, 1, :)
      sal_ice_out(:) = ice%sal_ice(ip, jp, 1, :)

      !$acc exit data delete(wet_mask, sf_0, dsf_dt, sw_dn, fprec, tfw_arr, fb_arr, &
      !$acc                  sst_arr, s_surf_arr, tsurf_out, h2o_ocn_to_ice, &
      !$acc                  h2o_ice_to_ocn, heat_to_ocn, sw_thru, snow_to_ice)
      call ice%exit_data()
      !$acc exit data delete(ice)
      call ice%destroy()
   end subroutine run_snowy_trajectory_no_flood

   ! Cold-ice total-energy closure (PR-60, RESUME item b)
   ! -----------------------------------------------------------------

   pure function col_enth(nk, m_snow, m_ice, enth_snow, enth_ice_bu) result(e)
      !! `E_col = m_snow*enth_snow + Sum_k (m_ice/nk)*enth_ice_bu(k)`
      !! [J/m²], from the BOTTOM-UP state. Flip-invariant (PR-60 plan
      !! §7/§11.5): `ice_column_step` gathers top-down internally, but
      !! every layer carries EQUAL mass (`m_ice/nk`,
      !! `rdb_ice_column.F90:702-704`), so this sum over ALL layers
      !! cannot see the flip — it is the same number computed from the
      !! top-down or bottom-up ordering.
      integer, intent(in) :: nk
      real(wp), intent(in) :: m_snow, m_ice, enth_snow
      real(wp), intent(in) :: enth_ice_bu(nk)
      real(wp) :: e
      integer :: k

      e = m_snow*enth_snow
      do k = 1, nk
         e = e + (m_ice/real(nk, wp))*enth_ice_bu(k)
      end do
   end function col_enth

   subroutine test_column_step_energy_closes_cold(error)
      !! The RESUME's item b, in its sharpest form (PR-60 plan §3.2).
      !! Seed via `make_initial_column(1.0, TFW=-1.9, TSURF_TARGET=-20,
      !! has_target=.true.)`. Force a column CLOSED TO HEAT: `sf_0=0`,
      !! `dsf_dt=0` (adiabatic top — zero SEB flux regardless of surface
      !! temperature), `sw_dn=0`, `fb=0`; `sst=tfw`, `s_surf=S_INIT`.
      !! Drive `ice_column_step` (the FULL end-to-end orchestrator —
      !! optics + conduction + resize + rebalance, not just the
      !! conduction solve `conduction_solve_energy_closure` covers) for
      !! 24 hourly steps. Every step: (a) `h2o_ice_to_ocn == 0` and
      !! `heat_to_ocn == 0` EXACTLY — the cold-branch precondition and
      !! scope fence (no melt anywhere, base only freezes: Stefan
      !! growth); (b) `sw_thru == 0` exactly (no shortwave in); (c) the
      !! closure: with a column exchanging energy with the outside ONLY
      !! through mass, `Delta(E_col) == h2o_ocn_to_ice*enth_ocean`
      !! (`enth_ocean = ice_enthalpy_liquid(sst,s_surf)`, the SAME public
      !! liquid-ocean formula the kernel itself uses — TRAP #1's exact
      !! gate: substituting the mushy `ice_enth_from_ts(tfw,sice)` here
      !! (~8-9x more negative) would break this by ~8x, not by the ~3%
      !! rate-error margin `stefan_growth` catches it at), bounded by
      !! SIS2's own `IMBALANCE_TOLERANCE` (1e-9, sum-of-magnitudes
      !! normalisation — same convention as
      !! `conduction_solve_energy_closure`, reused verbatim).
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: TFW = -1.9_wp
      real(wp), parameter :: TSURF_TARGET = -20.0_wp
      real(wp), parameter :: DTT = 3600.0_wp
      real(wp), parameter :: SW_DN = 0.0_wp
      real(wp), parameter :: FB = 0.0_wp
      integer, parameter :: N_STEPS = 24

      real(wp) :: m_ice, m_snow, enth_snow_pt
      real(wp) :: enth_ice_bu(NK), sal_ice_bu(NK)
      real(wp) :: tsurf, h2o_ocn_to_ice, h2o_ice_to_ocn, heat_to_ocn, sw_thru
      real(wp) :: snow_to_ice
      real(wp) :: e_before, e_after, enth_ocean, residual, norm, rel, max_rel_imbalance
      integer :: step

      checks: block

         call make_initial_column(1.0_wp, TFW, TSURF_TARGET, .true., &
                                  m_ice, m_snow, enth_ice_bu, enth_snow_pt, sal_ice_bu)
         max_rel_imbalance = 0.0_wp

         do step = 1, N_STEPS
            e_before = col_enth(NK, m_snow, m_ice, enth_snow_pt, enth_ice_bu)
            enth_ocean = ice_enthalpy_liquid(TFW, S_INIT)

            call ice_column_step(NK, m_snow, m_ice, enth_snow_pt, enth_ice_bu, sal_ice_bu, &
                                 0.0_wp, 0.0_wp, SW_DN, TFW, FB, TFW, S_INIT, DTT, &
                                 .false., 0.0_wp, &
                                 tsurf, h2o_ocn_to_ice, h2o_ice_to_ocn, heat_to_ocn, sw_thru, &
                                 snow_to_ice)

            call check(error, h2o_ice_to_ocn == 0.0_wp, &
                       "h2o_ice_to_ocn must be exactly zero on a cold step (scope fence)")
            if (allocated(error)) exit checks
            call check(error, heat_to_ocn == 0.0_wp, &
                       "heat_to_ocn must be exactly zero on a cold step (scope fence)")
            if (allocated(error)) exit checks
            call check(error, sw_thru == 0.0_wp, &
                       "sw_thru must be exactly zero at sw_dn=0")
            if (allocated(error)) exit checks

            e_after = col_enth(NK, m_snow, m_ice, enth_snow_pt, enth_ice_bu)
            residual = (e_after - e_before) - h2o_ocn_to_ice*enth_ocean
            norm = 1.0e-9_wp*(abs(e_before) + abs(e_after) + abs(h2o_ocn_to_ice*enth_ocean))
            rel = 0.0_wp
            if (norm > 0.0_wp) rel = abs(residual)/norm
            max_rel_imbalance = max(max_rel_imbalance, rel)
         end do

         call check(error, max_rel_imbalance <= 1.0_wp, &
                    "cold-column energy closure must stay within the 1e-9 budget every step")

      end block checks
   end subroutine test_column_step_energy_closes_cold

   subroutine test_column_step_energy_closes_cold_with_sw(error)
      !! As `column_step_energy_closes_cold` but `sw_dn=200 W/m^2`,
      !! `fb=5 W/m^2` (still cold overall — the base still nets freeze,
      !! not melt: `h2o_ice_to_ocn==0`/`heat_to_ocn==0` hold every step,
      !! asserted not assumed). Closes the loop through the OPTICS
      !! partition and the `fb` post-hoc fold (TRAP #3), neither
      !! exercised at `sw_dn=fb=0`. `albedo`/`sw_thru` are obtained by
      !! calling the PUBLIC `ice_optics_csim4` on the PRE-step state,
      !! exactly as `ice_column_step` does internally (`rdb_ice_column.
      !! F90:678-688`) — snow-free here, so the top-down layer-1 temp
      !! (== the bottom-up top layer, `enth_ice_bu(NK)`) feeds
      !! `ice_temp_from_en_s`. Full budget:
      !! `Delta(E_col) == DTT*((1-albedo)*sw_dn - sw_thru + fb) +
      !! h2o_ocn_to_ice*enth_ocean`. Also `sw_thru > 0` (the term is
      !! live, not incidentally zero) — this is also the analytical
      !! statement of the `sw_thru`-to-ocean leak PR-31 fixes: the
      !! magnitude subtracted from the ice's own budget here and given
      !! to nobody downstream is now an asserted number.
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: TFW = -1.9_wp
      real(wp), parameter :: TSURF_TARGET = -20.0_wp
      real(wp), parameter :: DTT = 3600.0_wp
      real(wp), parameter :: SW_DN = 200.0_wp
      real(wp), parameter :: FB = 5.0_wp
      integer, parameter :: N_STEPS = 24

      real(wp) :: m_ice, m_snow, enth_snow_pt
      real(wp) :: enth_ice_bu(NK), sal_ice_bu(NK)
      real(wp) :: tsurf, h2o_ocn_to_ice, h2o_ice_to_ocn, heat_to_ocn, sw_thru
      real(wp) :: snow_to_ice
      real(wp) :: e_before, e_after, enth_ocean, residual, norm, rel, max_rel_imbalance
      real(wp) :: ts_opt, albedo, abs_sfc, abs_snow, abs_ocn, abs_int, pen
      real(wp) :: abs_ice_lay(NK), expected
      integer :: step

      checks: block

         call make_initial_column(1.0_wp, TFW, TSURF_TARGET, .true., &
                                  m_ice, m_snow, enth_ice_bu, enth_snow_pt, sal_ice_bu)
         max_rel_imbalance = 0.0_wp

         do step = 1, N_STEPS
            e_before = col_enth(NK, m_snow, m_ice, enth_snow_pt, enth_ice_bu)
            enth_ocean = ice_enthalpy_liquid(TFW, S_INIT)

            ! Pre-step optics, snow-free: ts_opt == the top (bottom-up
            ! index NK) ice layer's temperature (rdb_ice_column.F90:
            ! 671-680 with m_snow==0).
            ts_opt = ice_temp_from_en_s(enth_ice_bu(NK), sal_ice_bu(NK))
            call ice_optics_csim4(NK, 0.0_wp, m_ice/ICE_RHO_ICE, ts_opt, sal_ice_bu(NK), &
                                  albedo, abs_sfc, abs_snow, abs_ice_lay, abs_ocn, &
                                  abs_int, pen)

            call ice_column_step(NK, m_snow, m_ice, enth_snow_pt, enth_ice_bu, sal_ice_bu, &
                                 0.0_wp, 0.0_wp, SW_DN, TFW, FB, TFW, S_INIT, DTT, &
                                 .false., 0.0_wp, &
                                 tsurf, h2o_ocn_to_ice, h2o_ice_to_ocn, heat_to_ocn, sw_thru, &
                                 snow_to_ice)

            call check(error, h2o_ice_to_ocn == 0.0_wp, &
                       "h2o_ice_to_ocn must be exactly zero on a cold step (scope fence)")
            if (allocated(error)) exit checks
            call check(error, heat_to_ocn == 0.0_wp, &
                       "heat_to_ocn must be exactly zero on a cold step (scope fence)")
            if (allocated(error)) exit checks
            call check(error, sw_thru > 0.0_wp, &
                       "sw_thru must be strictly positive (the optics term is live)")
            if (allocated(error)) exit checks

            e_after = col_enth(NK, m_snow, m_ice, enth_snow_pt, enth_ice_bu)
            expected = DTT*((1.0_wp - albedo)*SW_DN - sw_thru + FB) + h2o_ocn_to_ice*enth_ocean
            residual = (e_after - e_before) - expected
            norm = 1.0e-9_wp*(abs(e_before) + abs(e_after) + abs(h2o_ocn_to_ice*enth_ocean) &
                              + abs(expected))
            rel = 0.0_wp
            if (norm > 0.0_wp) rel = abs(residual)/norm
            max_rel_imbalance = max(max_rel_imbalance, rel)
         end do

         call check(error, max_rel_imbalance <= 1.0_wp, &
                    "cold-column energy closure (with SW+fb) must stay within the 1e-9 &
                    &budget every step")

      end block checks
   end subroutine test_column_step_energy_closes_cold_with_sw

   subroutine test_column_step_energy_closes_cold_thick(error)
      !! Sensitivity guard: `column_step_energy_closes_cold` at
      !! `h0=3.0 m` instead of `1.0 m`, same 24-step trajectory. A
      !! conduction-limited regime (`kk = K_ICE/hl_ice_eff` ~3x smaller,
      !! `m_freeze` per step ~3x smaller) — the closure residual must
      !! stay inside the same RELATIVE budget while the absolute terms
      !! move a decade. A test that only ever runs at one thickness
      !! cannot distinguish "conserves" from "the residual happens to be
      !! below tolerance at this one scale".
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: TFW = -1.9_wp
      real(wp), parameter :: TSURF_TARGET = -20.0_wp
      real(wp), parameter :: DTT = 3600.0_wp
      real(wp), parameter :: SW_DN = 0.0_wp
      real(wp), parameter :: FB = 0.0_wp
      real(wp), parameter :: H0 = 3.0_wp
      integer, parameter :: N_STEPS = 24

      real(wp) :: m_ice, m_snow, enth_snow_pt
      real(wp) :: enth_ice_bu(NK), sal_ice_bu(NK)
      real(wp) :: tsurf, h2o_ocn_to_ice, h2o_ice_to_ocn, heat_to_ocn, sw_thru
      real(wp) :: snow_to_ice
      real(wp) :: e_before, e_after, enth_ocean, residual, norm, rel, max_rel_imbalance
      integer :: step

      checks: block

         call make_initial_column(H0, TFW, TSURF_TARGET, .true., &
                                  m_ice, m_snow, enth_ice_bu, enth_snow_pt, sal_ice_bu)
         max_rel_imbalance = 0.0_wp

         do step = 1, N_STEPS
            e_before = col_enth(NK, m_snow, m_ice, enth_snow_pt, enth_ice_bu)
            enth_ocean = ice_enthalpy_liquid(TFW, S_INIT)

            call ice_column_step(NK, m_snow, m_ice, enth_snow_pt, enth_ice_bu, sal_ice_bu, &
                                 0.0_wp, 0.0_wp, SW_DN, TFW, FB, TFW, S_INIT, DTT, &
                                 .false., 0.0_wp, &
                                 tsurf, h2o_ocn_to_ice, h2o_ice_to_ocn, heat_to_ocn, sw_thru, &
                                 snow_to_ice)

            call check(error, h2o_ice_to_ocn == 0.0_wp, &
                       "h2o_ice_to_ocn must be exactly zero on a cold step (scope fence)")
            if (allocated(error)) exit checks
            call check(error, heat_to_ocn == 0.0_wp, &
                       "heat_to_ocn must be exactly zero on a cold step (scope fence)")
            if (allocated(error)) exit checks
            call check(error, sw_thru == 0.0_wp, &
                       "sw_thru must be exactly zero at sw_dn=0")
            if (allocated(error)) exit checks

            e_after = col_enth(NK, m_snow, m_ice, enth_snow_pt, enth_ice_bu)
            residual = (e_after - e_before) - h2o_ocn_to_ice*enth_ocean
            norm = 1.0e-9_wp*(abs(e_before) + abs(e_after) + abs(h2o_ocn_to_ice*enth_ocean))
            rel = 0.0_wp
            if (norm > 0.0_wp) rel = abs(residual)/norm
            max_rel_imbalance = max(max_rel_imbalance, rel)
         end do

         call check(error, max_rel_imbalance <= 1.0_wp, &
                    "cold-column energy closure (h0=3m) must stay within the 1e-9 budget &
                    &every step")

      end block checks
   end subroutine test_column_step_energy_closes_cold_thick

end module test_ocean_ice_column
