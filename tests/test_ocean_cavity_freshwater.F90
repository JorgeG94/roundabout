!! REAL ICE-SHELF MELTWATER MASS: the column grows, and every gram of it
!! is accounted for.
module test_ocean_cavity_freshwater
   !! ### What this suite gates
   !!
   !! `&ocean_cavity_melt_nml freshwater = "mass"` stops emulating the
   !! meltwater and adds it: `dh = m*dt/rho_0` on the top layer, `d(hS) =
   !! dh*s_ice`, `d(hT) = dh*T_b`, with the virtual salt flux removed
   !! again from the tracer so the dilution is not counted twice.  The
   !! derivation is in `rdb_ocean_cavity_flux`'s module docstring; this
   !! suite asserts the CONSEQUENCES, analytically where it can and
   !! through the production split solver where it must.
   !!
   !! Three things are easy to get wrong and each has its own gate here.
   !!
   !! 1. **Double counting.**  Applying both the volume and the virtual
   !!    salt flux gives a dilution twice as strong, and it still looks
   !!    plausible.  `top_salinity_follows_the_exact_dilution_law` pins
   !!    `S(t) = S0*h0/(h0 + m t/rho_0)` — the closed form, not a
   !!    tendency — so a double count is off by a factor the tolerance
   !!    cannot absorb.
   !!
   !! 2. **"Agrees with the old form".**  A double count ALSO agrees with
   !!    the virtual form to leading order, so "the two are close" proves
   !!    nothing.  The two equivalence gates therefore assert the EXACT
   !!    difference — `(S0 - s_ice)*eps^2/(1 + eps)` for salinity and
   !!    `-eps*(A + T0 - T_b)/(1 + eps)` for temperature, `eps = m dt/
   !!    (rho_0 h0)` — and separately that it is NOT zero.
   !!
   !! 3. **An untracked source.**  The maintainer's budget principle is
   !!    that the residual sits at round-off at EVERY step from step 0,
   !!    and that a change in the domain total is a TRACKED SOURCE.  Mass
   !!    now changes, so `mass_src` has to name it.  The 3-D gates assert
   !!    all three relative residuals at `1e-12` at every step, and that
   !!    the mass total actually GREW (a closed budget with a zero source
   !!    is the same statement as "nothing happened").
   !!
   !! ### `mem:separate` discipline
   !!
   !! The direct-kernel cases map every local array they hand a device
   !! kernel with `!$acc enter data copyin(...)`, push host-set values,
   !! read results back with `!$acc update self(...)` on the COMPONENT
   !! arrays (never an aggregate), and release with `!$acc exit data
   !! delete(...)`.  The engine cases use `engine_enter_data` and read
   !! through the PRODUCTION reducers, which carry their own `present(...)`
   !! clauses.  All of it is inert on the host build.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_constants, only: wp, RHO_WATER, H_VANISHED
   use rdb_config, only: config_t, read_config_from_string, validate_config
   use rdb_ocean_engine, only: ocean_engine_t, engine_setup, engine_enter_data, &
                               engine_step, engine_step_ice, engine_step_finalize, &
                               engine_exit_data, engine_teardown
   use rdb_ocean_console_stats, only: compute_total_h, compute_total_tracer, &
                                      ocean_budget_src, ocean_budget_out, &
                                      ocean_budget_stage_weight, &
                                      ocean_salt_src_sum, ocean_heat_src_sum
   use rdb_ocean_dyn, only: SPLIT_SCHEME_PRED_CORR
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   use rdb_comm_env, only: comm_env_init, comm_env_setup_roles
   use rdb_ocean_cavity_flux, only: cavity_mass_apply_impl, cavity_mass_salt_mirror_impl, &
                                    cavity_mass_totals_impl, cavity_comp_apply_impl, &
                                    cavity_comp_scale_tracer_impl, &
                                    cavity_comp_withdrawal, cavity_mass_thin_is_fatal
   use rdb_ocean_cavity_melt, only: parse_cavity_freshwater, parse_cavity_volume_comp, &
                                    CAVITY_FW_VIRTUAL, CAVITY_FW_MASS, CAVITY_FW_INVALID, &
                                    CAVITY_VC_NONE, CAVITY_VC_UNIFORM_OPEN, CAVITY_VC_INVALID
   implicit none
   private

   public :: collect_ocean_cavity_freshwater_tests

   logical, save :: comm_inited = .false.
      !! One-shot guard for `ensure_comm` (test-drive runs every case of a
      !! suite in ONE process, so the comm env must come up exactly once).

   ! ---- The single-column analytic setup --------------------------------
   real(wp), parameter :: RHO0 = 1027.51_wp
      !! ISOMIP+ `rho_ref` (Asay-Davis et al. 2016 Table 4) — deliberately
      !! NOT equal to `RHO_WATER`, so any place that confuses the
      !! Boussinesq reference with the console's mass scale shows up.
   real(wp), parameter :: CP = 3992.0_wp
      !! `SEAWATER_CP`; only the virtual heat stamp needs it.
   real(wp), parameter :: H0 = 20.0_wp
      !! Top-layer thickness (m).
   real(wp), parameter :: S0 = 34.5_wp
      !! Far-field / top-layer salinity (g/kg).
   real(wp), parameter :: T0 = -1.0_wp
      !! Top-layer temperature (degC).
   real(wp), parameter :: TB = -2.05_wp
      !! Interface temperature (degC), on the ISOMIP+ liquidus at depth.
   real(wp), parameter :: MELT = 1.0e-3_wp
      !! Melt mass flux (kg/m^2/s) — ~30 m/yr of ice, the ISOMIP+ Ocean0
      !! band.
   real(wp), parameter :: QOC = 340.0_wp
      !! Turbulent ocean -> interface heat flux (W/m^2), roughly `m*L`.
   real(wp), parameter :: DTC = 300.0_wp
      !! Column-test timestep (s).
   integer, parameter :: NSTEP_COL = 40
      !! Steps the column cases integrate for.

   ! ---- The 3-D engine setup --------------------------------------------
   integer, parameter :: N_STEPS = 20
      !! Steps of the FULL split solver every budget gate runs.
   real(wp), parameter :: DT = 300.0_wp
   real(wp), parameter :: BUDGET_TOL = 1.0e-12_wp
      !! The gate: a RELATIVE residual, so one number serves mass, salt
      !! and heat.  Not tuned — a closed budget in double precision
      !! drifts at `n_steps * eps * cancellation`, which for this domain
      !! is 1e-15..1e-13.  Raising it would mean something stopped being
      !! tracked.
   real(wp), parameter :: MELT_BED = 800.0_wp
      !! Flat bed (m) under the test cavity.
   real(wp), parameter :: DXC = 2000.0_wp
      !! Cell size (m).
   integer, parameter :: NXC = 40
      !! Cells in x.  Wide enough that the linear draft can taper to ZERO
      !! inside the domain and leave real OPEN OCEAN east of it — which
      !! the compensation gate needs (a fully ice-covered domain has
      !! nowhere to put the volume, and `cavity_comp_withdrawal` then
      !! correctly does nothing).
   integer, parameter :: NYC = 12
      !! Cells in y.

contains

   subroutine collect_ocean_cavity_freshwater_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("freshwater_and_compensation_parse", test_parse), &
                  new_unittest("mass_source_grows_the_column_by_m_dt_over_rho0", &
                               test_column_grows), &
                  new_unittest("top_salinity_follows_the_exact_dilution_law", &
                               test_dilution_law), &
                  new_unittest("top_temperature_follows_the_enthalpy_mixing_law", &
                               test_enthalpy_law), &
                  new_unittest("virtual_and_mass_differ_at_second_order_in_salinity", &
                               test_salt_second_order), &
                  new_unittest("virtual_and_mass_differ_by_the_derived_amount_in_heat", &
                               test_heat_first_order), &
                  new_unittest("freezing_removes_mass_and_concentrates_salt", &
                               test_freezing), &
                  new_unittest("freezing_through_h_vanished_is_clamped_counted_and_fatal", &
                               test_freezing_clamp), &
                  new_unittest("passive_tracer_is_diluted_by_the_added_volume", &
                               test_passive_dilution), &
                  new_unittest("pseudo_salt_mirror_matches_salinity_exactly", &
                               test_pseudo_salt_mirror), &
                  new_unittest("compensation_withdrawal_is_volume_over_open_area", &
                               test_comp_withdrawal), &
                  new_unittest("compensation_leaves_open_ocean_concentrations_alone", &
                               test_comp_concentrations), &
                  new_unittest("mass_budget_closes_from_step_zero_with_a_live_source", &
                               test_budget_mass), &
                  new_unittest("compensation_holds_the_domain_mass_constant", &
                               test_budget_compensated), &
                  new_unittest("virtual_default_is_bit_identical_to_the_explicit_knob", &
                               test_default_bit_identical), &
                  new_unittest("b0_salt_forcing_is_the_same_under_virtual_and_mass", &
                               test_b0_equivalence), &
                  new_unittest("mass_form_refusals_are_by_name", test_refusals) &
                  ]
   end subroutine collect_ocean_cavity_freshwater_tests

   ! ======================================================================
   ! Parsing
   ! ======================================================================

   subroutine test_parse(error)
      !! The two new enums, including the "reserved vs typo" distinction
      !! every other melt selector keeps.
      type(error_type), allocatable, intent(out) :: error
      call check(error, parse_cavity_freshwater("virtual") == CAVITY_FW_VIRTUAL, &
                 "virtual parses")
      if (allocated(error)) return
      call check(error, parse_cavity_freshwater("mass") == CAVITY_FW_MASS, "mass parses")
      if (allocated(error)) return
      call check(error, parse_cavity_freshwater("nope") == CAVITY_FW_INVALID, &
                 "a typo is INVALID, not a silent default")
      if (allocated(error)) return
      call check(error, parse_cavity_volume_comp("none") == CAVITY_VC_NONE, "none parses")
      if (allocated(error)) return
      call check(error, parse_cavity_volume_comp("uniform_open_ocean") == &
                 CAVITY_VC_UNIFORM_OPEN, "uniform_open_ocean parses")
      if (allocated(error)) return
      call check(error, parse_cavity_volume_comp("uniform") == CAVITY_VC_INVALID, &
                 "a near-miss is INVALID")
   end subroutine test_parse

   ! ======================================================================
   ! Single-column driving
   ! ======================================================================

   subroutine column_step(h, hs, ht, hps, sb, hb, melt_val, do_mass, n_thin)
      !! ONE outer step of the top layer of one column, exactly as the
      !! production stage runs it: the surface-flux stamp
      !! (`apply_surface_src_2d_impl`'s arithmetic, spelt out because the
      !! whole point is that the undo is its exact negation), then — when
      !! `do_mass` — the real-mass kernel.
      !!
      !! The `(1,1,1)` arrays are mapped, pushed, called and pulled every
      !! step: the cost is irrelevant at this size and it keeps the
      !! `mem:separate` contract honest for a kernel the GPU build will
      !! run over a whole domain.
      real(wp), intent(inout) :: h, hs, ht, hps, sb, hb
         !! Top-layer thickness, `h*S`, `h*T`, `h*pseudo-salt`, and the
         !! salt / heat surface budget contributors.
      real(wp), intent(in) :: melt_val
         !! Melt mass flux (kg/m^2/s).
      logical, intent(in) :: do_mass
         !! `.true.` ⇒ the real-mass form; `.false.` ⇒ virtual only.
      integer, intent(out) :: n_thin
         !! Clamped-withdrawal count this step.

      real(wp) :: a_h(1, 1), a_melt(1, 1), a_sfar(1, 1), a_tb(1, 1)
      real(wp) :: a_hl(1, 1, 1), a_hs(1, 1, 1), a_ht(1, 1, 1), a_hps(1, 1, 1)
      real(wp) :: a_sb(1, 1, 1), a_hb(1, 1, 1)
      real(wp) :: salt_cav, stamp_s, stamp_t

      ! --- the virtual stamp, on the host, exactly as the assembler +
      ! apply pair would leave it on a fully covered column (open-water
      ! factor 0 annihilates every atmospheric band, so Q_salt IS
      ! salt_cavity and Q_heat IS -q_ocean).
      salt_cav = -melt_val*(hs/h - 0.0_wp)
      stamp_s = (DTC/RHO0)*salt_cav
      stamp_t = (DTC/(RHO0*CP))*(-QOC)
      hs = hs + stamp_s
      sb = sb + stamp_s
      ht = ht + stamp_t
      hb = hb + stamp_t
      hps = hps + stamp_s
      n_thin = 0
      if (.not. do_mass) return

      a_h(1, 1) = 1.0_wp
      a_melt(1, 1) = melt_val
      a_sfar(1, 1) = (hs - stamp_s)/h
      a_tb(1, 1) = TB
      a_hl(1, 1, 1) = h
      a_hs(1, 1, 1) = hs
      a_ht(1, 1, 1) = ht
      a_hps(1, 1, 1) = hps
      a_sb(1, 1, 1) = sb
      a_hb(1, 1, 1) = hb

      !$acc enter data copyin(a_h, a_melt, a_sfar, a_tb, a_hl, a_hs, a_ht, &
      !$acc&                  a_hps, a_sb, a_hb)
      !$acc update device(a_h, a_melt, a_sfar, a_tb, a_hl, a_hs, a_ht, &
      !$acc&              a_hps, a_sb, a_hb)
      call cavity_mass_apply_impl(1, 1, 1, DTC/RHO0, DTC/RHO0, 0.0_wp, &
                                  a_h, a_melt, a_sfar, a_tb, &
                                  a_hl, a_hs, a_ht, a_sb, a_hb, n_thin)
      call cavity_mass_salt_mirror_impl(1, 1, 1, DTC/RHO0, DTC/RHO0, 0.0_wp, &
                                        a_h, a_melt, a_sfar, a_hps)
      !$acc update self(a_hl, a_hs, a_ht, a_hps, a_sb, a_hb)
      !$acc exit data delete(a_h, a_melt, a_sfar, a_tb, a_hl, a_hs, a_ht, &
      !$acc&                 a_hps, a_sb, a_hb)

      h = a_hl(1, 1, 1)
      hs = a_hs(1, 1, 1)
      ht = a_ht(1, 1, 1)
      hps = a_hps(1, 1, 1)
      sb = a_sb(1, 1, 1)
      hb = a_hb(1, 1, 1)
   end subroutine column_step

   ! ======================================================================
   ! The analytic column gates
   ! ======================================================================

   subroutine test_column_grows(error)
      !! Constant `m` for `N` steps ⇒ the column is thicker by exactly
      !! `N*m*dt/rho_0`.  "Exactly" to the accumulation bound: `N` adds,
      !! each rounding at `eps*h`.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: h, hs, ht, hps, sb, hb, dh, tol
      integer :: n, n_thin

      call seed(h, hs, ht, hps, sb, hb)
      dh = MELT*(DTC/RHO0)
      do n = 1, NSTEP_COL
         call column_step(h, hs, ht, hps, sb, hb, MELT, .true., n_thin)
      end do
      tol = 4.0_wp*real(NSTEP_COL, wp)*epsilon(1.0_wp)*(H0 + NSTEP_COL*dh)
      call check(error, abs(h - (H0 + real(NSTEP_COL, wp)*dh)) <= tol, &
                 "thickness grew by exactly N*m*dt/rho_0")
      if (allocated(error)) return
      call check(error, dh > 0.0_wp .and. h > H0, "and it actually grew")
   end subroutine test_column_grows

   subroutine test_dilution_law(error)
      !! `s_ice = 0`, no mixing ⇒ the salinity follows the CLOSED FORM
      !! `S(t) = S0*h0/(h0 + m t/rho_0)`.  A double count (volume AND
      !! virtual flux) halves the numerator's survival and misses this by
      !! ~`eps` of the signal, which the bound below cannot absorb.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: h, hs, ht, hps, sb, hb, dh, s_exact, tol
      integer :: n, n_thin

      call seed(h, hs, ht, hps, sb, hb)
      dh = MELT*(DTC/RHO0)
      do n = 1, NSTEP_COL
         call column_step(h, hs, ht, hps, sb, hb, MELT, .true., n_thin)
      end do
      s_exact = S0*H0/(H0 + real(NSTEP_COL, wp)*dh)
      ! The undo is the exact negation of the stamp, so `hs` only carries
      ! `N` roundings of `(hs + X) - X`, each at `eps*|hs|`.
      tol = 8.0_wp*real(NSTEP_COL, wp)*epsilon(1.0_wp)*S0
      call check(error, abs(hs/h - s_exact) <= tol, "S follows the exact dilution law")
      if (allocated(error)) return
      call check(error, s_exact < S0 - 1.0e-4_wp, &
                 "and the dilution is a real, resolvable freshening")
   end subroutine test_dilution_law

   subroutine test_enthalpy_law(error)
      !! With the turbulent flux ALSO on, the closed form is the mixing
      !! law `T(t) = (h0 T0 + A h0 + (m t/rho_0) T_b)/(h0 + m t/rho_0)`
      !! per step, accumulated — asserted here as the equivalent
      !! statement that `h*T` grew by exactly `N*(A*h0_stamp + dh*T_b)`,
      !! where the stamp is the same every step because `Q_heat` is.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: h, hs, ht, hps, sb, hb, dh, stamp_t, ht_exact, tol
      integer :: n, n_thin

      call seed(h, hs, ht, hps, sb, hb)
      dh = MELT*(DTC/RHO0)
      stamp_t = (DTC/(RHO0*CP))*(-QOC)
      do n = 1, NSTEP_COL
         call column_step(h, hs, ht, hps, sb, hb, MELT, .true., n_thin)
      end do
      ht_exact = H0*T0 + real(NSTEP_COL, wp)*(stamp_t + dh*TB)
      tol = 8.0_wp*real(NSTEP_COL, wp)*epsilon(1.0_wp)*abs(H0*T0)
      call check(error, abs(ht - ht_exact) <= tol, &
                 "h*T grew by the turbulent stamp plus the meltwater enthalpy dh*T_b")
      if (allocated(error)) return
      ! The meltwater arrives COLDER than the column, so the enthalpy
      ! term cools it on top of the turbulent loss — both signs negative.
      call check(error, dh*TB < 0.0_wp .and. ht < H0*T0, &
                 "and the added water cools the column")
   end subroutine test_enthalpy_law

   subroutine test_salt_second_order(error)
      !! ONE step from the same state under both forms.  The difference
      !! is EXACTLY `(S0 - s_ice)*eps^2/(1 + eps)` — second order — and
      !! it is NOT zero.  Asserting only "they agree" would also pass on
      !! a double count, which agrees to leading order too.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: hm, hsm, htm, hpsm, sbm, hbm
      real(wp) :: hv, hsv, htv, hpsv, sbv, hbv
      real(wp) :: eps, ds_m, ds_v, diff, predicted
      integer :: n_thin

      call seed(hm, hsm, htm, hpsm, sbm, hbm)
      call seed(hv, hsv, htv, hpsv, sbv, hbv)
      call column_step(hm, hsm, htm, hpsm, sbm, hbm, MELT, .true., n_thin)
      call column_step(hv, hsv, htv, hpsv, sbv, hbv, MELT, .false., n_thin)

      eps = MELT*(DTC/RHO0)/H0
      ds_m = hsm/hm - S0
      ds_v = hsv/hv - S0
      diff = ds_m - ds_v
      predicted = (S0 - 0.0_wp)*eps*eps/(1.0_wp + eps)
      call check(error, abs(diff - predicted) <= 1.0e-9_wp*abs(predicted) + &
                 8.0_wp*epsilon(1.0_wp)*S0, &
                 "the salinity difference is exactly (S0-s_ice)*eps^2/(1+eps)")
      if (allocated(error)) return
      call check(error, abs(diff) > 0.5_wp*predicted .and. predicted > 0.0_wp, &
                 "and it is NOT zero — the two forms are not the same answer")
      if (allocated(error)) return
      call check(error, abs(diff) < 1.0e-3_wp*abs(ds_v), &
                 "but it is a small correction to the first-order term")
   end subroutine test_salt_second_order

   subroutine test_heat_first_order(error)
      !! The temperature twin: `dT_mass - dT_virtual = -eps*(A + T0 -
      !! T_b)/(1 + eps)` with `A = -q_ocean*dt/(rho_0 cp h0)`.  Unlike
      !! the salt case this is FIRST order in `eps` — the virtual form
      !! deliberately drops the `-m*c_w*(T_w - T_b)` dilution term — and
      !! the gate asserts the size, not that it vanishes.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: hm, hsm, htm, hpsm, sbm, hbm
      real(wp) :: hv, hsv, htv, hpsv, sbv, hbv
      real(wp) :: eps, a_term, dt_m, dt_v, diff, predicted
      integer :: n_thin

      call seed(hm, hsm, htm, hpsm, sbm, hbm)
      call seed(hv, hsv, htv, hpsv, sbv, hbv)
      call column_step(hm, hsm, htm, hpsm, sbm, hbm, MELT, .true., n_thin)
      call column_step(hv, hsv, htv, hpsv, sbv, hbv, MELT, .false., n_thin)

      eps = MELT*(DTC/RHO0)/H0
      a_term = -QOC*DTC/(RHO0*CP*H0)
      dt_m = htm/hm - T0
      dt_v = htv/hv - T0
      diff = dt_m - dt_v
      predicted = -eps*(a_term + T0 - TB)/(1.0_wp + eps)
      call check(error, abs(diff - predicted) <= 1.0e-9_wp*abs(predicted) + &
                 8.0_wp*epsilon(1.0_wp)*abs(T0), &
                 "the temperature difference is exactly -eps*(A + T0 - T_b)/(1+eps)")
      if (allocated(error)) return
      call check(error, abs(diff) > 0.0_wp, "and it is not zero")
   end subroutine test_heat_first_order

   subroutine test_freezing(error)
      !! `m < 0` withdraws mass and CONCENTRATES the salt — the same
      !! closed form with a negative `dh`, so `S > S0`.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: h, hs, ht, hps, sb, hb, dh, s_exact
      integer :: n, n_thin

      call seed(h, hs, ht, hps, sb, hb)
      dh = (-MELT)*(DTC/RHO0)
      do n = 1, NSTEP_COL
         call column_step(h, hs, ht, hps, sb, hb, -MELT, .true., n_thin)
         call check(error, n_thin == 0, "no clamping at this thickness")
         if (allocated(error)) return
      end do
      s_exact = S0*H0/(H0 + real(NSTEP_COL, wp)*dh)
      call check(error, h < H0 .and. s_exact > S0, "freezing thins and concentrates")
      if (allocated(error)) return
      call check(error, abs(hs/h - s_exact) <= &
                 8.0_wp*real(NSTEP_COL, wp)*epsilon(1.0_wp)*S0, &
                 "and it follows the same closed form with dh < 0")
   end subroutine test_freezing

   subroutine test_freezing_clamp(error)
      !! A freezing column whose top layer is thinner than `|m|dt/rho_0`
      !! must NOT be driven through `H_VANISHED`.  The withdrawal is
      !! clamped to what is there, the column is COUNTED, and the count
      !! is FATAL (`cavity_mass_thin_is_fatal`) because a clamped
      !! withdrawal no longer matches the tracked mass source.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: a_h(1, 1), a_melt(1, 1), a_sfar(1, 1), a_tb(1, 1)
      real(wp) :: a_hl(1, 1, 1), a_hs(1, 1, 1), a_ht(1, 1, 1)
      real(wp) :: a_sb(1, 1, 1), a_hb(1, 1, 1)
      real(wp) :: h_small, big_melt
      integer :: n_thin

      ! A 1 mm top layer against a freezing flux that would take a metre.
      h_small = 1.0e-3_wp
      big_melt = -1.0_wp*RHO0/DTC
      a_h(1, 1) = 1.0_wp
      a_melt(1, 1) = big_melt
      a_sfar(1, 1) = S0
      a_tb(1, 1) = TB
      a_hl(1, 1, 1) = h_small
      a_hs(1, 1, 1) = h_small*S0
      a_ht(1, 1, 1) = h_small*T0
      a_sb(1, 1, 1) = 0.0_wp
      a_hb(1, 1, 1) = 0.0_wp

      !$acc enter data copyin(a_h, a_melt, a_sfar, a_tb, a_hl, a_hs, a_ht, a_sb, a_hb)
      !$acc update device(a_h, a_melt, a_sfar, a_tb, a_hl, a_hs, a_ht, a_sb, a_hb)
      call cavity_mass_apply_impl(1, 1, 1, DTC/RHO0, DTC/RHO0, 0.0_wp, &
                                  a_h, a_melt, a_sfar, a_tb, &
                                  a_hl, a_hs, a_ht, a_sb, a_hb, n_thin)
      !$acc update self(a_hl)
      !$acc exit data delete(a_h, a_melt, a_sfar, a_tb, a_hl, a_hs, a_ht, a_sb, a_hb)

      call check(error, n_thin == 1, "the starved column is counted")
      if (allocated(error)) return
      call check(error, a_hl(1, 1, 1) >= H_VANISHED, &
                 "and clamped AT the vanish marker, never through it")
      if (allocated(error)) return
      call check(error, a_hl(1, 1, 1) == H_VANISHED, &
                 "pinned EXACTLY at it, never a ulp below")
      if (allocated(error)) return
      call check(error, cavity_mass_thin_is_fatal(1), "any clamped column is fatal")
      if (allocated(error)) return
      call check(error,.not. cavity_mass_thin_is_fatal(0), "zero is not")
   end subroutine test_freezing_clamp

   subroutine test_passive_dilution(error)
      !! A passive tracer with ZERO concentration in the meltwater needs
      !! no kernel of its own: `h` grows, `h*C` does not, so `C` falls by
      !! exactly the dilution factor.  The gate is that the mass kernel
      !! leaves an untouched `h*C` alone — i.e. that nothing secretly
      !! adds tracer with the volume.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: h, hs, ht, hps, sb, hb, hc, dh, c_exact
      integer :: n, n_thin

      call seed(h, hs, ht, hps, sb, hb)
      hc = H0*1.0_wp
      dh = MELT*(DTC/RHO0)
      do n = 1, NSTEP_COL
         call column_step(h, hs, ht, hps, sb, hb, MELT, .true., n_thin)
      end do
      c_exact = 1.0_wp*H0/(H0 + real(NSTEP_COL, wp)*dh)
      call check(error, abs(hc/h - c_exact) <= &
                 8.0_wp*real(NSTEP_COL, wp)*epsilon(1.0_wp), &
                 "an untouched h*C is diluted by exactly the added volume")
      if (allocated(error)) return
      call check(error, c_exact < 1.0_wp, "and the dilution is real")
   end subroutine test_passive_dilution

   subroutine test_pseudo_salt_mirror(error)
      !! Pseudo-salt is "given S's surface salt flux".  Under the mass
      !! form that flux is `dh*s_ice` and the virtual increment must be
      !! taken back from it too — otherwise its deviation from `S` would
      !! measure the BOOKKEEPING rather than the transport path, which is
      !! the one thing it exists to isolate.  Seeded to `S`, it must stay
      !! bit-for-bit equal to it.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: h, hs, ht, hps, sb, hb
      integer :: n, n_thin

      call seed(h, hs, ht, hps, sb, hb)
      do n = 1, NSTEP_COL
         call column_step(h, hs, ht, hps, sb, hb, MELT, .true., n_thin)
      end do
      call check(error, hps == hs, &
                 "pseudo-salt tracks salinity bit-for-bit under the mass form")
   end subroutine test_pseudo_salt_mirror

   subroutine test_comp_withdrawal(error)
      !! The uniform sink's arithmetic, and its no-open-ocean fallback —
      !! a fully ice-covered domain has nowhere to put the volume, and
      !! the honest answer is to leave it in (and let the mass budget
      !! report the growth) rather than divide by zero.
      type(error_type), allocatable, intent(out) :: error
      call check(error, abs(cavity_comp_withdrawal(1000.0_wp, 250.0_wp) - 4.0_wp) <= &
                 1.0e-14_wp, "dw = volume/open area")
      if (allocated(error)) return
      call check(error, cavity_comp_withdrawal(1000.0_wp, 0.0_wp) == 0.0_wp, &
                 "no open ocean ⇒ no withdrawal, no division by zero")
   end subroutine test_comp_withdrawal

   subroutine test_comp_concentrations(error)
      !! The sink removes VOLUME carrying each cell's own `T` and `S`, so
      !! no concentration anywhere changes — including every passive
      !! tracer, which is why `comp_scale` is published rather than the
      !! thickness change.  A covered cell is left completely alone.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 4, NY = 3
      real(wp) :: wet(NX, NY), cover(NX, NY), scale(NX, NY)
      real(wp) :: hl(NX, NY, 1), hs(NX, NY, 1), ht(NX, NY, 1), hc(NX, NY, 1)
      real(wp) :: sb(NX, NY, 1), hb(NX, NY, 1)
      real(wp) :: dw, s_before, t_before, c_before
      integer :: i, j, n_thin

      wet = 1.0_wp
      cover = 0.0_wp
      cover(1, :) = 1.0_wp          ! one covered column: untouched
      wet(NX, :) = 0.0_wp           ! one land column: untouched
      do j = 1, NY
         do i = 1, NX
            hl(i, j, 1) = 50.0_wp + real(i, wp)
            hs(i, j, 1) = hl(i, j, 1)*S0
            ht(i, j, 1) = hl(i, j, 1)*T0
            hc(i, j, 1) = hl(i, j, 1)*7.0_wp
         end do
      end do
      sb = 0.0_wp
      hb = 0.0_wp
      dw = 0.25_wp
      s_before = hs(2, 2, 1)/hl(2, 2, 1)
      t_before = ht(2, 2, 1)/hl(2, 2, 1)
      c_before = hc(2, 2, 1)/hl(2, 2, 1)

      !$acc enter data copyin(wet, cover, scale, hl, hs, ht, hc, sb, hb)
      !$acc update device(wet, cover, hl, hs, ht, hc, sb, hb)
      call cavity_comp_apply_impl(NX, NY, 1, dw, wet, cover, hl, hs, ht, sb, hb, &
                                  scale, n_thin)
      call cavity_comp_scale_tracer_impl(NX, NY, 1, scale, hc)
      !$acc update self(hl, hs, ht, hc, sb, hb, scale)
      !$acc exit data delete(wet, cover, scale, hl, hs, ht, hc, sb, hb)

      call check(error, n_thin == 0, "no clamping on a 50 m top layer")
      if (allocated(error)) return
      call check(error, abs(hl(2, 2, 1) - (52.0_wp - dw)) <= 1.0e-12_wp, &
                 "an open-ocean column loses exactly dw")
      if (allocated(error)) return
      call check(error, abs(hs(2, 2, 1)/hl(2, 2, 1) - s_before) <= 1.0e-12_wp*S0, &
                 "its salinity is unchanged")
      if (allocated(error)) return
      call check(error, abs(ht(2, 2, 1)/hl(2, 2, 1) - t_before) <= 1.0e-12_wp*abs(T0), &
                 "its temperature is unchanged")
      if (allocated(error)) return
      call check(error, abs(hc(2, 2, 1)/hl(2, 2, 1) - c_before) <= 1.0e-12_wp*c_before, &
                 "and so is every passive tracer's concentration")
      if (allocated(error)) return
      call check(error, abs(sb(2, 2, 1) + dw*S0) <= 1.0e-9_wp, &
                 "the salt the parcel carried out is TRACKED, not lost")
      if (allocated(error)) return
      call check(error, hl(1, 2, 1) == 51.0_wp .and. scale(1, 2) == 1.0_wp, &
                 "an ICE-COVERED column is left alone, scale exactly 1")
      if (allocated(error)) return
      call check(error, hl(NX, 2, 1) == 50.0_wp + real(NX, wp) .and. &
                 scale(NX, 2) == 1.0_wp, "and so is a LAND column")
   end subroutine test_comp_concentrations

   pure subroutine seed(h, hs, ht, hps, sb, hb)
      !! The single-column initial state, in one place so the paired
      !! virtual/mass gates cannot start from different numbers.
      real(wp), intent(out) :: h, hs, ht, hps, sb, hb
      h = H0
      hs = H0*S0
      ht = H0*T0
      hps = H0*S0
      sb = 0.0_wp
      hb = 0.0_wp
   end subroutine seed

   ! ======================================================================
   ! Namelists for the 3-D gates
   ! ======================================================================

   function melt_nml(freshwater, compensation, extra) result(nml)
      !! A partly grounded sloping cavity over a `MELT_BED` basin — the
      !! same geometry `test_ocean_cavity_grounded_budget` uses for its
      !! melt gate, because it is the one that carries BOTH grounded
      !! columns and covered wet ones, and the melt source only exists on
      !! the latter.  No wind, no restoring, walls all round: the only
      !! things that may move a total are the ones the budget names.
      character(len=*), intent(in) :: freshwater
         !! `&ocean_cavity_melt_nml freshwater`; empty ⇒ the key is
         !! omitted entirely (the default path).
      character(len=*), intent(in) :: compensation
         !! `volume_compensation`; empty ⇒ omitted.
      character(len=*), intent(in) :: extra
         !! Extra namelist groups appended verbatim (the refusal cases).
      character(len=:), allocatable :: nml
      character(len=32) :: bed_s, dx_s, nx_s, ny_s
      character(len=:), allocatable :: melt_keys

      write (bed_s, '(F12.2)') MELT_BED
      write (dx_s, '(F12.2)') DXC
      write (nx_s, '(I0)') NXC
      write (ny_s, '(I0)') NYC

      melt_keys = ""
      if (len_trim(freshwater) > 0) then
         melt_keys = melt_keys//", freshwater = '"//trim(freshwater)//"'"
      end if
      if (len_trim(compensation) > 0) then
         melt_keys = melt_keys//", volume_compensation = '"//trim(compensation)//"'"
      end if

      nml = "&sim_nml sim_type = 'ocean' /"//new_line("a")// &
            "&grid_nml nx = "//trim(adjustl(nx_s))//", ny = "// &
            trim(adjustl(ny_s))//", nghost = 2, dx = "//trim(adjustl(dx_s))// &
            ", dy = "//trim(adjustl(dx_s))//" /"//new_line("a")// &
            "&nonhydrostatic_nml nz_layers = 4 /"//new_line("a")// &
            "&time_nml t_end = 1.0e9, dt_fixed = 300.0 /"//new_line("a")// &
            "&physics_nml coriolis_f = -1.409e-4, wind_stress_x = 0.0, "// &
            "wind_stress_y = 0.0 /"//new_line("a")// &
            "&ocean_topo_nml topo_config = 'flat', max_depth = "// &
            trim(adjustl(bed_s))//", wind_config = 'constant', "// &
            "taux_magnitude = 0.0 /"//new_line("a")// &
            "&tracer_nml initial_salinity = 34.5, T_init_bottom = 1.0, "// &
            "T_init_surface = 3.0 /"//new_line("a")// &
            "&ocean_pgf_nml form = 'fv_mom6', p_top_in_bc = .true. /"//new_line("a")// &
            "&ocean_eos_nml eos = 'linear', tfreeze_set = 'isomip' /"//new_line("a")// &
            "&vcoord_nml vcoord_type = 'sigma' /"//new_line("a")// &
            "&ocean_bc_nml west = 'wall', east = 'wall', south = 'wall', "// &
            "north = 'wall' /"//new_line("a")// &
            ! `auto_n_inner`, not a pinned count: the open ocean east of
            ! the calving taper carries the FULL 800 m bed, so the
            ! gravity-wave CFL there is ~13 substeps and a hand-pinned 12
            ! silently runs the barotropic loop unstable.
            "&ocean_bt_nml split_scheme = 'pred_corr', auto_n_inner = .true. /"// &
            new_line("a")// &
            ! A draft of 750 m at the west wall thinning at 1e-2 m/m:
            ! the western five columns GROUND (their water column falls
            ! below `h_min_cavity`) and become land, the middle of the
            ! domain is covered ocean where the melt kernel delivers, and
            ! the formula clips at zero near x = 75 km so the last
            ! columns are genuinely OPEN ocean.  20 m of ice base per
            ! 2 km cell, the same order as the neighbouring
            ! grounded-budget gate's slope -- a steeper lid is a
            ! barotropic shock, not a cavity.
            "&ocean_cavity_dyn_nml enable = .true., draft_config = 'linear', "// &
            "draft_depth = 750.0, draft_slope = -1.0e-2, draft_x0 = 0.0, "// &
            "h_min_cavity = 150.0 /"//new_line("a")// &
            "&ocean_forcing_nml enable_components = .true. /"//new_line("a")// &
            "&ocean_cavity_melt_nml enable = .true., exchange_law = 'const_gamma', "// &
            "gamma_t = 2.2e-2, gamma_s = -1.0, cdrag_top = 2.5e-3, "// &
            "ice_conduction = 'insulating', s_ice = 0.0, far_field_depth = 10.0"// &
            melt_keys//" /"//new_line("a")// &
            "&ocean_tdrag_nml enable = .true., form = 'quadratic', cd = 2.5e-3 /"// &
            new_line("a")// &
            "&ocean_hvisc_nml nu_h = 100.0 /"//new_line("a")// &
            "&ocean_diag_nml enabled = .false. /"//new_line("a")// &
            "&output_nml output_to_file = .false. /"//new_line("a")//extra
   end function melt_nml

   ! ======================================================================
   ! Engine driving + budget bookkeeping (the console's own terms)
   ! ======================================================================

   subroutine ensure_comm()
      !! Bring the MPI comm env up once per process.  Every engine case
      !! below reaches the communicator — `engine_setup` builds the decomp
      !! and the budget reducers in `totals` are collectives — so on an
      !! `RDB_ENABLE_MPI=ON` build the first of them hits `MPI_Comm_f2c`
      !! before `MPI_INIT` and OpenMPI aborts the process (CI: SEGFAULT).
      !! No-op-equivalent on the single-rank backend.  Finalised by the
      !! shared per-test main, so there is no teardown here.  Same pattern
      !! as `test_ocean_console_stats_efp` / `test_driver_ocean`.
      if (.not. comm_inited) then
         call comm_env_init()
         call comm_env_setup_roles(.false.)
         comm_inited = .true.
      end if
   end subroutine ensure_comm

   subroutine start_engine(nml, engine, cfg, ok)
      character(len=*), intent(in) :: nml
      type(ocean_engine_t), intent(inout) :: engine
      type(config_t), intent(inout) :: cfg
      logical, intent(out) :: ok
      integer :: ierr
      ok = .false.
      call ensure_comm()
      call read_config_from_string(nml, cfg, ierr=ierr)
      if (ierr /= 0) return
      call validate_config(cfg, ierr)
      if (ierr /= 0) return
      call engine_setup(engine, cfg, ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call engine_enter_data(engine, cfg)
      ok = .true.
   end subroutine start_engine

   subroutine advance(engine, cfg, t, ierr)
      !! One outer step in `driver_run_ocean`'s MANDATED ORDER.  The melt
      !! package writes its owned components inside `engine_step_finalize`,
      !! so a gate that skipped it would never see a source at all.
      type(ocean_engine_t), intent(inout) :: engine
      type(config_t), intent(in) :: cfg
      real(wp), intent(in) :: t
      integer, intent(out) :: ierr
      call engine_step(engine, DT, t, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call engine_step_ice(engine, cfg, DT, t, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call engine_step_finalize(engine, DT, t, ierr=ierr)
   end subroutine advance

   subroutine totals(engine, total_mass, total_salt, total_heat, &
                     src_salt, src_heat, out_salt, out_heat, mass_out, mass_src)
      !! Exactly what `ocean_console_stats_report` forms, through the
      !! PRODUCTION reducers and the production `ocean_budget_*` helpers.
      type(ocean_engine_t), intent(inout) :: engine
      real(wp), intent(out) :: total_mass, total_salt, total_heat
      real(wp), intent(out) :: src_salt, src_heat, out_salt, out_heat
      real(wp), intent(out) :: mass_out, mass_src
      real(wp) :: total_h, bud_w

      associate (ms => engine%state%multilayer, mt => engine%state%metrics, &
                 ng => engine%grid%nghost)
         bud_w = ocean_budget_stage_weight(engine%state%dyn%split_scheme == &
                                           SPLIT_SCHEME_PRED_CORR)
         total_h = compute_total_h(ms%h_layer, mt%areaT, ng)
         total_mass = total_h*RHO_WATER
         total_salt = compute_total_tracer(ms%tracers(ms%idx_salinity)%hTr, &
                                           mt%areaT, ng)*RHO_WATER
         total_heat = compute_total_tracer(ms%tracers(ms%idx_temperature)%hTr, &
                                           mt%areaT, ng)*RHO_WATER
         src_salt = ocean_budget_src( &
                    ocean_salt_src_sum( &
                    compute_total_tracer(ms%salt_budget_surface, mt%areaT, ng), &
                    compute_total_tracer(ms%salt_budget_sponge, mt%areaT, ng)), &
                    stage_weight=bud_w)
         out_salt = ocean_budget_out( &
                    compute_total_tracer(ms%salt_budget_horiz_adv, mt%areaT, ng), &
                    compute_total_tracer(ms%salt_budget_hdiff, mt%areaT, ng), &
                    stage_weight=bud_w)
         src_heat = ocean_budget_src( &
                    ocean_heat_src_sum( &
                    compute_total_tracer(ms%heat_budget_surface, mt%areaT, ng), &
                    compute_total_tracer(ms%heat_budget_geothermal, mt%areaT, ng), &
                    compute_total_tracer(ms%heat_budget_sponge, mt%areaT, ng)), &
                    stage_weight=bud_w)
         out_heat = ocean_budget_out( &
                    compute_total_tracer(ms%heat_budget_horiz_adv, mt%areaT, ng), &
                    compute_total_tracer(ms%heat_budget_hdiff, mt%areaT, ng), &
                    stage_weight=bud_w)
         mass_out = ms%mass_out
         mass_src = ms%mass_src
      end associate
   end subroutine totals

   pure function residual(total, ref, out_term, src_term) result(rel)
      !! `rdb_console_stats::emit_drift_line`'s residual, relative — what
      !! the `Error` column prints.
      real(wp), intent(in) :: total, ref, out_term, src_term
      real(wp) :: rel
      if (abs(ref) > 1.0e-30_wp) then
         rel = ((total - ref) + out_term - src_term)/ref
      else
         rel = 0.0_wp
      end if
   end function residual

   ! ======================================================================
   ! The 3-D gates
   ! ======================================================================

   subroutine budget_run(error, nml, label, expect_growth)
      !! Step the cavity for `N_STEPS` and assert the three RELATIVE
      !! residuals at EVERY step, step 0 included — and that the mass
      !! total moved the way the configuration says it should.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: nml
      character(len=*), intent(in) :: label
      logical, intent(in) :: expect_growth
         !! `.true.` ⇒ uncompensated: the domain must GAIN mass, and by
         !! exactly the integrated melt volume.  `.false.` ⇒ compensated:
         !! the domain mass must stay put to round-off.

      type(ocean_engine_t) :: engine
      type(config_t) :: cfg
      logical :: ok
      integer :: n, ierr
      real(wp) :: t, m0, s0_, h0_
      real(wp) :: mt_, st_, ht_, src_s, src_h, out_s, out_h, m_out, m_src
      real(wp) :: r_m, r_s, r_h, worst_m, worst_s, worst_h, vol_acc

      call start_engine(nml, engine, cfg, ok)
      call check(error, ok, label//": engine starts")
      if (allocated(error)) return

      gate: block
         call totals(engine, m0, s0_, h0_, src_s, src_h, out_s, out_h, m_out, m_src)
         call check(error, ieee_is_finite(m0) .and. abs(m0) > 0.0_wp, &
                    label//": step-0 mass total is finite and non-zero")
         if (allocated(error)) exit gate
         call check(error, m_src == 0.0_wp, label//": no source before any step")
         if (allocated(error)) exit gate

         worst_m = 0.0_wp
         worst_s = 0.0_wp
         worst_h = 0.0_wp
         vol_acc = 0.0_wp
         t = 0.0_wp
         do n = 1, N_STEPS
            call advance(engine, cfg, t, ierr)
            call check(error, ierr == OCEAN_STATUS_OK, label//": step succeeds")
            if (allocated(error)) exit gate
            t = t + DT
            ! The melt volume the slot reports for THIS thermo step.
            ! Accumulated here so the gate can compare the tracked source
            ! against `sum(m*area*dt)/rho_0` independently of the
            ! accumulator that produced it.
            vol_acc = vol_acc + engine%state%cavity_flux%melt_volume_step

            call totals(engine, mt_, st_, ht_, src_s, src_h, out_s, out_h, m_out, m_src)
            r_m = residual(mt_, m0, m_out, m_src)
            r_s = residual(st_, s0_, out_s, src_s)
            r_h = residual(ht_, h0_, out_h, src_h)
            worst_m = max(worst_m, abs(r_m))
            worst_s = max(worst_s, abs(r_s))
            worst_h = max(worst_h, abs(r_h))
            call check(error, abs(r_m) <= BUDGET_TOL, label//": mass residual")
            if (allocated(error)) exit gate
            call check(error, abs(r_s) <= BUDGET_TOL, label//": salt residual")
            if (allocated(error)) exit gate
            call check(error, abs(r_h) <= BUDGET_TOL, label//": heat residual")
            if (allocated(error)) exit gate
         end do

         call check(error, vol_acc > 0.0_wp, &
                    label//": the cavity actually melted (a closed budget with a "// &
                    "zero source proves nothing)")
         if (allocated(error)) exit gate

         if (expect_growth) then
            ! The tracked source IS the integrated melt volume scaled to
            ! the console's own mass measure.
            call check(error, abs(m_src - RHO_WATER*vol_acc) <= &
                       1.0e-12_wp*abs(m_src), &
                       label//": mass_src = RHO_WATER * sum(m*area*dt)/rho_0")
            if (allocated(error)) exit gate
            call check(error, mt_ - m0 > 0.0_wp, label//": the domain GAINED mass")
            if (allocated(error)) exit gate
            ! `mt_ - m0` is a difference of two ~1e13 kg totals, so its
            ! own cancellation floor is `eps*m0` — the bound is therefore
            ! relative to `m0`, the same denominator the console's Error
            ! column uses.  The source itself is many decades above that
            ! floor (asserted next), so this is not a vacuous bound.
            call check(error, abs((mt_ - m0) - m_src) <= BUDGET_TOL*abs(m0), &
                       label//": and gained exactly the tracked source")
            if (allocated(error)) exit gate
            call check(error, m_src > 1.0e-9_wp*abs(m0), &
                       label//": the growth is resolvable, not a round-off artefact")
         else
            call check(error, abs(m_src) <= 1.0e-9_wp*abs(m0), &
                       label//": source and sink cancel in the accumulator")
            if (allocated(error)) exit gate
            call check(error, abs(mt_ - m0) <= 1.0e-12_wp*abs(m0), &
                       label//": the domain mass is constant to round-off")
         end if
      end block gate

      call engine_exit_data(engine)
      call engine_teardown(engine)
   end subroutine budget_run

   subroutine test_budget_mass(error)
      type(error_type), allocatable, intent(out) :: error
      call budget_run(error, melt_nml("mass", "", ""), "mass", .true.)
   end subroutine test_budget_mass

   subroutine test_budget_compensated(error)
      type(error_type), allocatable, intent(out) :: error
      call budget_run(error, melt_nml("mass", "uniform_open_ocean", ""), &
                      "mass+compensation", .false.)
   end subroutine test_budget_compensated

   subroutine test_default_bit_identical(error)
      !! The knob default must change NOTHING: an omitted `freshwater`
      !! and an explicit `freshwater='virtual'` must produce the same
      !! state, bit for bit, after five steps of the full solver.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t) :: e_a, e_b
      type(config_t) :: c_a, c_b
      logical :: ok_a, ok_b
      integer :: n, ierr
      real(wp) :: t
      real(wp), allocatable :: h_a(:, :, :), h_b(:, :, :)
      real(wp), allocatable :: s_a(:, :, :), s_b(:, :, :)

      call start_engine(melt_nml("", "", ""), e_a, c_a, ok_a)
      call start_engine(melt_nml("virtual", "", ""), e_b, c_b, ok_b)
      call check(error, ok_a .and. ok_b, "both engines start")
      if (allocated(error)) return

      gate: block
         t = 0.0_wp
         do n = 1, 5
            call advance(e_a, c_a, t, ierr)
            call check(error, ierr == OCEAN_STATUS_OK, "default run steps")
            if (allocated(error)) exit gate
            call advance(e_b, c_b, t, ierr)
            call check(error, ierr == OCEAN_STATUS_OK, "explicit-virtual run steps")
            if (allocated(error)) exit gate
            t = t + DT
         end do
         call pull_h_s(e_a, h_a, s_a)
         call pull_h_s(e_b, h_b, s_b)
         call check(error, all(h_a == h_b), "h_layer is bit-identical")
         if (allocated(error)) exit gate
         call check(error, all(s_a == s_b), "h*S is bit-identical")
         if (allocated(error)) exit gate
         call check(error, e_a%state%multilayer%mass_src == 0.0_wp, &
                    "and the virtual form tracks no mass source")
      end block gate

      call engine_exit_data(e_a)
      call engine_teardown(e_a)
      call engine_exit_data(e_b)
      call engine_teardown(e_b)
   end subroutine test_default_bit_identical

   subroutine test_b0_equivalence(error)
      !! The surface buoyancy forcing KPP and EPBL read (`Q_salt`, which
      !! is `beta*Q_salt/rho_0` in `B_0`) must be the SAME under both
      !! forms — that is the whole reason the virtual flux is still
      !! assembled under `"mass"` and removed from the TRACER instead.
      !!
      !! Decisive version: after ONE outer step the two runs have had
      !! identical state everywhere the mass kernel could act (melt is
      !! still zero during step 1), so the melt solve at the end of step
      !! 1 sees identical inputs and `Q_salt` must come out BIT-IDENTICAL.
      !! It is also asserted NON-ZERO, so the gate cannot pass on two
      !! runs that both forgot to assemble anything.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t) :: e_v, e_m
      type(config_t) :: c_v, c_m
      logical :: ok_v, ok_m
      integer :: ierr
      real(wp), allocatable :: qs_v(:, :), qs_m(:, :), sc_v(:, :), sc_m(:, :)

      call start_engine(melt_nml("virtual", "", ""), e_v, c_v, ok_v)
      call start_engine(melt_nml("mass", "", ""), e_m, c_m, ok_m)
      call check(error, ok_v .and. ok_m, "both engines start")
      if (allocated(error)) return

      gate: block
         call advance(e_v, c_v, 0.0_wp, ierr)
         call check(error, ierr == OCEAN_STATUS_OK, "virtual run steps")
         if (allocated(error)) exit gate
         call advance(e_m, c_m, 0.0_wp, ierr)
         call check(error, ierr == OCEAN_STATUS_OK, "mass run steps")
         if (allocated(error)) exit gate

         call pull_q_salt(e_v, qs_v, sc_v)
         call pull_q_salt(e_m, qs_m, sc_m)
         call check(error, maxval(abs(sc_v)) > 0.0_wp, &
                    "the virtual run assembled a non-zero salt_cavity")
         if (allocated(error)) exit gate
         call check(error, all(sc_v == sc_m), &
                    "salt_cavity — the meltwater buoyancy signal — is identical")
         if (allocated(error)) exit gate
         call check(error, all(qs_v == qs_m), &
                    "so the Q_salt KPP/EPBL build B_0 from is identical")
      end block gate

      call engine_exit_data(e_v)
      call engine_teardown(e_v)
      call engine_exit_data(e_m)
      call engine_teardown(e_m)
   end subroutine test_b0_equivalence

   subroutine pull_h_s(engine, h, hs)
      !! Host copies of `h_layer` and `h*S`.  COMPONENT arrays only, and
      !! the registry indirection is resolved through `associate` before
      !! the `update self` names it.
      type(ocean_engine_t), intent(inout) :: engine
      real(wp), allocatable, intent(out) :: h(:, :, :), hs(:, :, :)
      associate (ms => engine%state%multilayer)
         associate (hl => ms%h_layer, hsal => ms%tracers(ms%idx_salinity)%hTr)
            !$acc update self(hl, hsal) if_present
            allocate (h, source=hl)
            allocate (hs, source=hsal)
         end associate
      end associate
   end subroutine pull_h_s

   subroutine pull_q_salt(engine, q_salt, salt_cavity)
      !! Host copies of the assembled `Q_salt` and the cavity's own
      !! component.  COMPONENT arrays only.
      type(ocean_engine_t), intent(inout) :: engine
      real(wp), allocatable, intent(out) :: q_salt(:, :), salt_cavity(:, :)
      associate (sf => engine%state%surface_flux)
         associate (qs => sf%Q_salt, sc => sf%salt_cavity)
            !$acc update self(qs, sc) if_present
            allocate (q_salt, source=qs)
            allocate (salt_cavity, source=sc)
         end associate
      end associate
   end subroutine pull_q_salt

   ! ======================================================================
   ! Refusals
   ! ======================================================================

   subroutine test_refusals(error)
      !! Every envelope hole is refused BY NAME at configure, and each
      !! case is asserted to PARSE first — otherwise the test would pass
      !! on a malformed namelist rather than on the rule it names.
      type(error_type), allocatable, intent(out) :: error
      call expect_refused(error, melt_nml("mass", "", &
                                          "&ocean_wetdry_nml enable = .true. /"//new_line("a")), &
                          "mass x wet/dry")
      if (allocated(error)) return
      call expect_refused(error, melt_nml("mass", "", &
                                          "&ocean_vmix_nml dt_tracer_advect_ratio = 2 /"// &
                                          new_line("a")), "mass x windowed tracer advect")
      if (allocated(error)) return
      call expect_refused(error, melt_nml("virtual", "uniform_open_ocean", ""), &
                          "compensation without the mass form")
      if (allocated(error)) return
      call expect_accepted(error, melt_nml("mass", "uniform_open_ocean", ""), &
                           "mass + compensation")
   end subroutine test_refusals

   subroutine expect_refused(error, nml, label)
      type(error_type), allocatable, intent(inout) :: error
      character(len=*), intent(in) :: nml, label
      type(config_t) :: cfg
      integer :: ierr
      call read_config_from_string(nml, cfg, ierr=ierr)
      call check(error, ierr == 0, label//": the namelist PARSES")
      if (allocated(error)) return
      call validate_config(cfg, ierr)
      call check(error, ierr /= 0, label//": and validate_config refuses it")
   end subroutine expect_refused

   subroutine expect_accepted(error, nml, label)
      type(error_type), allocatable, intent(inout) :: error
      character(len=*), intent(in) :: nml, label
      type(config_t) :: cfg
      integer :: ierr
      call read_config_from_string(nml, cfg, ierr=ierr)
      call check(error, ierr == 0, label//": parses")
      if (allocated(error)) return
      call validate_config(cfg, ierr)
      call check(error, ierr == 0, label//": and is accepted")
   end subroutine expect_accepted

end module test_ocean_cavity_freshwater
