!! Tests for the ice-shelf basal-melt kernel (`rdb_ocean_cavity_melt`).
module test_ocean_cavity_melt
   !! Five gates on the cavity three-equation kernel — which is KERNEL
   !! ONLY at this point: nothing here runs inside the ocean step.
   !!
   !!   (a) ORACLE — 48 golden cases replayed against
   !!       `cavity_solve_melt`, covering `const_gamma` / `hj99` /
   !!       `yung25`, both liquidus SETS (applied through the real
   !!       `eos_apply_tfreeze_set` handle path, never a local copy of
   !!       the coefficients), both shipped ice-conduction modes, both
   !!       melting and freezing, three ice salinities, a pressure
   !!       (ice-pump) sweep and a friction-velocity sweep.
   !!   (b) RESIDUALS — the three equations are satisfied to round-off by
   !!       the returned state, independently of any golden number.
   !!   (c) LIMITS — the four analytic limits the prototype established,
   !!       the sign conventions and the ice pump.
   !!   (d) GUARDS — NaN / Inf / negative salinity / zero `u*` / `f = 0`
   !!       return the documented status AND the safe state, never a
   !!       finite-but-wrong melt rate.
   !!   (e) DEVICE — the whole solver runs inside a `do concurrent` loop
   !!       under the `mem:separate` data directives, so the GPU build
   !!       actually exercises the `!$acc routine seq` path.
   !!
   !! ============ PROVENANCE OF THE GOLDEN NUMBERS ============
   !! Source: `python_prototypes/ice_shelf_melt/oracle.json`
   !!   schema          `ice_shelf_melt/oracle/v1`
   !!   generator       `ice_shelf_melt/make_oracle.py`
   !!   generating commit  38cbe556772edb96a7b55fb4756329d9aed7eed7
   !!                      (the generator's parent commit, as the oracle
   !!                       records it; `git_dirty = true`)
   !!   86 cases at 17 significant digits; the 48 embedded here are
   !!   every case whose exchange law, liquidus set, ice mode and
   !!   constants bundle this Fortran kernel implements:
   !!
   !!     law_const_isomip__{warm_deep_fast,mild_deep_mid,
   !!                        quiescent_warm,freeze_mild}
   !!     law_const_jenkins10__{same four}
   !!     law_hj99_south__{same four}
   !!     law_yung25__{same four}
   !!     ice_insulating__{mild_deep_mid,warm_deep_fast,freeze_mild}
   !!     ice_adv_diff__{same three}
   !!     ice_adv_diff_warm__{same three}
   !!     liq_roundabout__{mild_deep_mid,cold_deep_slow}
   !!     liq_isomip__{mild_deep_mid,cold_deep_slow}
   !!     sice_{0,1,4}__mild_deep_mid
   !!     pump_{0,1,2,3,4}
   !!     ustar_yung25_{0,1,2,3,4}
   !!     ustar_hj99_south_{0,1,2,3,4}
   !!     jenkins10_native__mild_deep_mid
   !!
   !!   DELIBERATELY NOT EMBEDDED, and why:
   !!     * `law_{jenkins91,rosevear22,vt19,mk18,burchard22}_*`,
   !!       `ustar_rosevear22_*`, `mk18_45deg_*`,
   !!       `rosevear_below_floor_*`, `burchard_{hard_switch,native}_*`
   !!       — RESERVED laws; the dispatcher returns
   !!       `CAVITY_MELT_NOT_IMPLEMENTED` and `test_reserved_laws_refuse`
   !!       pins that instead.
   !!     * `ice_diffusive__*` — RESERVED ice mode (same treatment).
   !!     * `liq_hj99__*`, `liq_yung25__*`, `hj99_native__*` — those
   !!       liquidus coefficient sets are not among the named sets
   !!       `eos_apply_tfreeze_set` ships, and this suite refuses to
   !!       hand-inject coefficients the EOS handle would not produce.
   !!
   !! TOLERANCE: relative 1e-12 with a per-field absolute floor, never
   !! bit-equality — GPU/FMA builds differ from gfortran in the last ulp,
   !! the liquidus is evaluated in a different association order here
   !! than in the Python (`(l1*S + l3*p) + l2` vs `l1*S + l2 + l3*p`),
   !! and the `hj99`/`yung25` answers come through a bracketed bisection
   !! whose endpoints inherit those ulps.
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
   use rdb_constants, only: wp, RHO_WATER
   use rdb_eos, only: eos_t, eos_freezing_point, eos_apply_tfreeze_set, &
                      TFREEZE_SET_SEAICE, TFREEZE_SET_ISOMIP
   use rdb_ocean_surface_flux, only: SEAWATER_CP
   use rdb_ocean_cavity_melt, only: ocean_cavity_const_t, ocean_cavity_exchange_t, &
                                    ocean_cavity_ice_t, ocean_cavity_solution_t, &
                                    parse_cavity_exchange_law, parse_cavity_ice_mode, &
                                    cavity_ustar, cavity_exchange_velocities, &
                                    cavity_three_equation, cavity_two_equation, &
                                    cavity_solve_melt, cavity_heat_fluxes, &
                                    cavity_salt_fluxes, cavity_m_weq_from_mass, &
                                    cavity_melt_columns, &
                                    CAVITY_LAW_INVALID, CAVITY_LAW_CONST_GAMMA, &
                                    CAVITY_LAW_HJ99, CAVITY_LAW_YUNG25, &
                                    CAVITY_LAW_JENKINS91, CAVITY_LAW_ROSEVEAR22, &
                                    CAVITY_LAW_VT19, CAVITY_LAW_MK18, &
                                    CAVITY_LAW_BURCHARD22, CAVITY_LAW_JENKINS21, &
                                    CAVITY_ICE_INVALID, CAVITY_ICE_INSULATING, &
                                    CAVITY_ICE_ADV_DIFF, CAVITY_ICE_DIFFUSIVE, &
                                    CAVITY_MELT_OK, CAVITY_MELT_NONFINITE_INPUT, &
                                    CAVITY_MELT_BAD_INPUT, CAVITY_MELT_NOT_IMPLEMENTED, &
                                    CAVITY_MELT_LAW_INVALID, CAVITY_MELT_NO_CORIOLIS, &
                                    CAVITY_GAMMA_T_ISOMIP, CAVITY_GAMMA_S_ISOMIP, &
                                    CAVITY_CD_ISOMIP, CAVITY_U_TIDE_ISOMIP, &
                                    CAVITY_USTAR_MIN_YUNG25, CAVITY_L_PLUS_NEUTRAL
   use testdrive, only: new_unittest, unittest_type, error_type, check
   implicit none
   private

   public :: collect_ocean_cavity_melt_tests

   real(wp), parameter :: RTOL = 1.0e-12_wp
      !! Relative tolerance against the golden numbers.  Explicitly NOT
      !! bit-equality: see the module header.

   ! Per-field absolute floors, each chosen well below the smallest
   ! non-zero magnitude the embedded cases carry (min |value| in the
   ! oracle subset is given in brackets).
   real(wp), parameter :: ATOL_T = 1.0e-12_wp      !! degC          [1.29]
   real(wp), parameter :: ATOL_S = 1.0e-12_wp      !! g/kg          [17.6]
   real(wp), parameter :: ATOL_M = 1.0e-18_wp      !! kg/m^2/s      [1.4e-6]
   real(wp), parameter :: ATOL_GAMMA = 1.0e-20_wp  !! m/s           [4.4e-9]
   real(wp), parameter :: ATOL_Q = 1.0e-12_wp      !! W/m^2         [0.47]
   real(wp), parameter :: ATOL_B = 1.0e-22_wp      !! m^2/s^3       [2.4e-10]
   real(wp), parameter :: ATOL_LP = 1.0e-14_wp     !! dimensionless [0.022]

   ! ---- The golden cases, generated from oracle.json (see the header) ----
   integer, parameter :: NCASE = 48
   character(len=35), parameter :: CASE_ID(NCASE) = [character(len=35) :: &
                                                     "law_const_isomip__warm_deep_fast", "law_const_isomip__mild_deep_mid", &
                                                     "law_const_isomip__quiescent_warm", "law_const_isomip__freeze_mild", &
                                                     "law_const_jenkins10__warm_deep_fast", "law_const_jenkins10__mild_deep_mid", &
                                                     "law_const_jenkins10__quiescent_warm", "law_const_jenkins10__freeze_mild", &
                                                     "law_hj99_south__warm_deep_fast", "law_hj99_south__mild_deep_mid", &
                                                     "law_hj99_south__quiescent_warm", "law_hj99_south__freeze_mild", &
                                                     "law_yung25__warm_deep_fast", "law_yung25__mild_deep_mid", &
                                                     "law_yung25__quiescent_warm", "law_yung25__freeze_mild", &
                                                     "ice_insulating__mild_deep_mid", "ice_insulating__warm_deep_fast", &
                                                     "ice_insulating__freeze_mild", "ice_adv_diff__mild_deep_mid", &
                                                     "ice_adv_diff__warm_deep_fast", "ice_adv_diff__freeze_mild", &
                                                     "ice_adv_diff_warm__mild_deep_mid", "ice_adv_diff_warm__warm_deep_fast", &
                                                     "ice_adv_diff_warm__freeze_mild", "liq_roundabout__mild_deep_mid", &
                                                     "liq_roundabout__cold_deep_slow", "liq_isomip__mild_deep_mid", &
                                                     "liq_isomip__cold_deep_slow", "sice_0__mild_deep_mid", &
                                                     "sice_1__mild_deep_mid", "sice_4__mild_deep_mid", &
                                                     "pump_0", "pump_1", &
                                                     "pump_2", "pump_3", &
                                                     "pump_4", "ustar_yung25_0", &
                                                     "ustar_yung25_1", "ustar_yung25_2", &
                                                     "ustar_yung25_3", "ustar_yung25_4", &
                                                     "ustar_hj99_south_0", "ustar_hj99_south_1", &
                                                     "ustar_hj99_south_2", "ustar_hj99_south_3", &
                                                     "ustar_hj99_south_4", "jenkins10_native__mild_deep_mid" &
                                                     ]
   integer, parameter :: CASE_LAW(NCASE) = [ &
                         1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, &
                         4, 4, 4, 4, 1, 1, 1, 1, 1, 1, 1, 1, &
                         1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
                         1, 4, 4, 4, 4, 4, 2, 2, 2, 2, 2, 1 &
                         ]
   integer, parameter :: CASE_ICE_MODE(NCASE) = [ &
                         1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
                         1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, &
                         2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
                         1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 &
                         ]
   integer, parameter :: CASE_TFREEZE_SET(NCASE) = [ &
                         1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
                         1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
                         1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, &
                         1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2 &
                         ]
   integer, parameter :: CASE_CONST_SET(NCASE) = [ &
                         1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
                         1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
                         1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
                         1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2 &
                         ]
   real(wp), parameter :: CASE_T_W(NCASE) = [ &
                          1.0_wp, -1.0_wp, 0.5_wp, -2.35_wp, &
                          1.0_wp, -1.0_wp, 0.5_wp, -2.35_wp, &
                          1.0_wp, -1.0_wp, 0.5_wp, -2.35_wp, &
                          1.0_wp, -1.0_wp, 0.5_wp, -2.35_wp, &
                          -1.0_wp, 1.0_wp, -2.35_wp, -1.0_wp, &
                          1.0_wp, -2.35_wp, -1.0_wp, 1.0_wp, &
                          -2.35_wp, -1.0_wp, -1.9_wp, -1.0_wp, &
                          -1.9_wp, -1.0_wp, -1.0_wp, -1.0_wp, &
                          -2.05_wp, -2.05_wp, -2.05_wp, -2.05_wp, &
                          -2.05_wp, 0.5_wp, 0.5_wp, 0.5_wp, &
                          0.5_wp, 0.5_wp, 0.5_wp, 0.5_wp, &
                          0.5_wp, 0.5_wp, 0.5_wp, -1.0_wp &
                          ]
   real(wp), parameter :: CASE_S_W(NCASE) = [ &
                          34.7_wp, 34.5_wp, 34.6_wp, 34.5_wp, &
                          34.7_wp, 34.5_wp, 34.6_wp, 34.5_wp, &
                          34.7_wp, 34.5_wp, 34.6_wp, 34.5_wp, &
                          34.7_wp, 34.5_wp, 34.6_wp, 34.5_wp, &
                          34.5_wp, 34.7_wp, 34.5_wp, 34.5_wp, &
                          34.7_wp, 34.5_wp, 34.5_wp, 34.7_wp, &
                          34.5_wp, 34.5_wp, 34.3_wp, 34.5_wp, &
                          34.3_wp, 34.5_wp, 34.5_wp, 34.5_wp, &
                          34.5_wp, 34.5_wp, 34.5_wp, 34.5_wp, &
                          34.5_wp, 34.6_wp, 34.6_wp, 34.6_wp, &
                          34.6_wp, 34.6_wp, 34.6_wp, 34.6_wp, &
                          34.6_wp, 34.6_wp, 34.6_wp, 34.5_wp &
                          ]
   real(wp), parameter :: CASE_P_B(NCASE) = [ &
                          4500000.0_wp, 4500000.0_wp, 4500000.0_wp, 4500000.0_wp, &
                          4500000.0_wp, 4500000.0_wp, 4500000.0_wp, 4500000.0_wp, &
                          4500000.0_wp, 4500000.0_wp, 4500000.0_wp, 4500000.0_wp, &
                          4500000.0_wp, 4500000.0_wp, 4500000.0_wp, 4500000.0_wp, &
                          4500000.0_wp, 4500000.0_wp, 4500000.0_wp, 4500000.0_wp, &
                          4500000.0_wp, 4500000.0_wp, 4500000.0_wp, 4500000.0_wp, &
                          4500000.0_wp, 4500000.0_wp, 8100000.0_wp, 4500000.0_wp, &
                          8100000.0_wp, 4500000.0_wp, 4500000.0_wp, 4500000.0_wp, &
                          450000.0_wp, 1800000.0_wp, 4500000.0_wp, 7200000.0_wp, &
                          9000000.0_wp, 4500000.0_wp, 4500000.0_wp, 4500000.0_wp, &
                          4500000.0_wp, 4500000.0_wp, 4500000.0_wp, 4500000.0_wp, &
                          4500000.0_wp, 4500000.0_wp, 4500000.0_wp, 4500000.0_wp &
                          ]
   real(wp), parameter :: CASE_U_STAR(NCASE) = [ &
                          0.01_wp, 0.002_wp, 0.0001_wp, 0.002_wp, &
                          0.01_wp, 0.002_wp, 0.0001_wp, 0.002_wp, &
                          0.01_wp, 0.002_wp, 0.0001_wp, 0.002_wp, &
                          0.01_wp, 0.002_wp, 0.0001_wp, 0.002_wp, &
                          0.002_wp, 0.01_wp, 0.002_wp, 0.002_wp, &
                          0.01_wp, 0.002_wp, 0.002_wp, 0.01_wp, &
                          0.002_wp, 0.002_wp, 0.0005_wp, 0.002_wp, &
                          0.0005_wp, 0.002_wp, 0.002_wp, 0.002_wp, &
                          0.002_wp, 0.002_wp, 0.002_wp, 0.002_wp, &
                          0.002_wp, 0.0001_wp, 0.0005_wp, 0.002_wp, &
                          0.01_wp, 0.05_wp, 0.0001_wp, 0.0005_wp, &
                          0.002_wp, 0.01_wp, 0.05_wp, 0.002_wp &
                          ]
   real(wp), parameter :: CASE_S_I(NCASE) = [ &
                          0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                          0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                          0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                          0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                          0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                          0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                          0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                          0.0_wp, 0.0_wp, 1.0_wp, 4.0_wp, &
                          0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                          0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                          0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                          0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp &
                          ]
   real(wp), parameter :: CASE_T_ICE(NCASE) = [ &
                          0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                          0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                          0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                          0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                          0.0_wp, 0.0_wp, 0.0_wp, -25.0_wp, &
                          -25.0_wp, -25.0_wp, -10.0_wp, -10.0_wp, &
                          -10.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                          0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                          0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                          0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                          0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                          0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp &
                          ]
   real(wp), parameter :: CASE_GAMMA_T_COEFF(NCASE) = [ &
                          0.022_wp, 0.022_wp, 0.022_wp, 0.022_wp, &
                          0.011_wp, 0.011_wp, 0.011_wp, 0.011_wp, &
                          0.022_wp, 0.022_wp, 0.022_wp, 0.022_wp, &
                          0.022_wp, 0.022_wp, 0.022_wp, 0.022_wp, &
                          0.022_wp, 0.022_wp, 0.022_wp, 0.022_wp, &
                          0.022_wp, 0.022_wp, 0.022_wp, 0.022_wp, &
                          0.022_wp, 0.022_wp, 0.022_wp, 0.022_wp, &
                          0.022_wp, 0.022_wp, 0.022_wp, 0.022_wp, &
                          0.022_wp, 0.022_wp, 0.022_wp, 0.022_wp, &
                          0.022_wp, 0.022_wp, 0.022_wp, 0.022_wp, &
                          0.022_wp, 0.022_wp, 0.022_wp, 0.022_wp, &
                          0.022_wp, 0.022_wp, 0.022_wp, 0.011_wp &
                          ]
   real(wp), parameter :: CASE_GAMMA_S_COEFF(NCASE) = [ &
                          0.0006285714285714285_wp, 0.0006285714285714285_wp, 0.0006285714285714285_wp, 0.0006285714285714285_wp, &
                          0.00031_wp, 0.00031_wp, 0.00031_wp, 0.00031_wp, &
                          0.0006285714285714285_wp, 0.0006285714285714285_wp, 0.0006285714285714285_wp, 0.0006285714285714285_wp, &
                          0.0006285714285714285_wp, 0.0006285714285714285_wp, 0.0006285714285714285_wp, 0.0006285714285714285_wp, &
                          0.0006285714285714285_wp, 0.0006285714285714285_wp, 0.0006285714285714285_wp, 0.0006285714285714285_wp, &
                          0.0006285714285714285_wp, 0.0006285714285714285_wp, 0.0006285714285714285_wp, 0.0006285714285714285_wp, &
                          0.0006285714285714285_wp, 0.0006285714285714285_wp, 0.0006285714285714285_wp, 0.0006285714285714285_wp, &
                          0.0006285714285714285_wp, 0.0006285714285714285_wp, 0.0006285714285714285_wp, 0.0006285714285714285_wp, &
                          0.0006285714285714285_wp, 0.0006285714285714285_wp, 0.0006285714285714285_wp, 0.0006285714285714285_wp, &
                          0.0006285714285714285_wp, 0.0006285714285714285_wp, 0.0006285714285714285_wp, 0.0006285714285714285_wp, &
                          0.0006285714285714285_wp, 0.0006285714285714285_wp, 0.0006285714285714285_wp, 0.0006285714285714285_wp, &
                          0.0006285714285714285_wp, 0.0006285714285714285_wp, 0.0006285714285714285_wp, 0.00031_wp &
                          ]
   real(wp), parameter :: CASE_F_COR(NCASE) = [ &
                          -0.00014_wp, -0.00014_wp, -0.00014_wp, -0.00014_wp, &
                          -0.00014_wp, -0.00014_wp, -0.00014_wp, -0.00014_wp, &
                          -0.00014_wp, -0.00014_wp, -0.00014_wp, -0.00014_wp, &
                          -0.00014_wp, -0.00014_wp, -0.00014_wp, -0.00014_wp, &
                          -0.00014_wp, -0.00014_wp, -0.00014_wp, -0.00014_wp, &
                          -0.00014_wp, -0.00014_wp, -0.00014_wp, -0.00014_wp, &
                          -0.00014_wp, -0.00014_wp, -0.00014_wp, -0.00014_wp, &
                          -0.00014_wp, -0.00014_wp, -0.00014_wp, -0.00014_wp, &
                          -0.00014_wp, -0.00014_wp, -0.00014_wp, -0.00014_wp, &
                          -0.00014_wp, -0.00014_wp, -0.00014_wp, -0.00014_wp, &
                          -0.00014_wp, -0.00014_wp, -0.00014_wp, -0.00014_wp, &
                          -0.00014_wp, -0.00014_wp, -0.00014_wp, -0.00014_wp &
                          ]
   real(wp), parameter :: WANT_T_B(NCASE) = [ &
                          -1.2966288520010791_wp, -1.755839190691419_wp, -1.385467845392689_wp, -2.2678453014161155_wp, &
                          -1.291289355792308_wp, -1.7526209979011487_wp, -1.3803863426766836_wp, -2.268357889705781_wp, &
                          -1.4588528370810727_wp, -1.8471822608048039_wp, -1.7220700635717148_wp, -2.257787906252856_wp, &
                          -1.3467369646372465_wp, -1.8378880669105255_wp, -1.7629833973749456_wp, -2.2630808829103537_wp, &
                          -1.755839190691419_wp, -1.2966288520010791_wp, -2.2678453014161155_wp, -1.786113521967163_wp, &
                          -1.3483610734966445_wp, -2.2678453014161155_wp, -1.7671133151653051_wp, -1.3164344382340127_wp, &
                          -2.2678453014161155_wp, -1.755839190691419_wp, -2.2352324859798873_wp, -1.758138660405212_wp, &
                          -2.245137020877307_wp, -1.755839190691419_wp, -1.7648508594890104_wp, -1.7929837361776209_wp, &
                          -1.9651362346385384_wp, -2.0211748791454243_wp, -2.1368335954402933_wp, -2.2571869195412626_wp, &
                          -2.3399720926983667_wp, -1.7629833973749456_wp, -1.6457789061756818_wp, -1.5343395799256767_wp, &
                          -1.432978402296661_wp, -1.432978402296661_wp, -1.7220700635717148_wp, -1.6132173003436983_wp, &
                          -1.55536733524022_wp, -1.536477790573384_wp, -1.5494570719083318_wp, -1.75478076506155_wp &
                          ]
   real(wp), parameter :: WANT_S_B(NCASE) = [ &
                          17.736645407427392_wp, 26.24054056835961_wp, 19.38181195171646_wp, 35.722135211409544_wp, &
                          17.637765848005703_wp, 26.180944405576827_wp, 19.287710049568215_wp, 35.73162758714409_wp, &
                          20.740793279279124_wp, 27.93207890379266_wp, 25.61518636243916_wp, 35.53588715283067_wp, &
                          18.664573419208267_wp, 27.759964202046767_wp, 26.37284069212862_wp, 35.633905239080626_wp, &
                          26.24054056835961_wp, 17.736645407427392_wp, 35.722135211409544_wp, 26.80117633272524_wp, &
                          18.694649509197117_wp, 35.722135211409544_wp, 26.44932065120935_wp, 18.103415522852085_wp, &
                          35.722135211409544_wp, 26.24054056835961_wp, 30.098194184812726_wp, 26.221442589968795_wp, &
                          29.989651324211295_wp, 26.24054056835961_wp, 26.40742332387056_wp, 26.928402521807794_wp, &
                          35.76391175256553_wp, 34.91916442861897_wp, 33.29599250815358_wp, 31.75975776928264_wp, &
                          30.782816531451232_wp, 26.37284069212862_wp, 24.202387151401513_wp, 22.138695924549566_wp, &
                          20.261637079567794_wp, 20.261637079567794_wp, 25.61518636243916_wp, 23.599394450809225_wp, &
                          22.528098800744814_wp, 22.178292418025627_wp, 22.418649479783923_wp, 26.162840576990398_wp &
                          ]
   real(wp), parameter :: WANT_M_MASS(NCASE) = [ &
                          0.006179993351987645_wp, 0.0004067771916715866_wp, 5.0736011348810716e-5_wp, -4.4213978296108075e-5_wp, &
                          0.0030828126568946076_wp, 0.0002025226130701427_wp, 2.529963665392638e-5_wp, -2.1969056881825025e-5_wp, &
                          0.0030237749887498775_wp, 0.0002114177652402276_wp, 1.5233235600714158e-5_wp, -2.6352925940988435e-5_wp, &
                          0.0034444522681428267_wp, 0.00011649602749917285_wp, 1.403987687260532e-6_wp, -2.5515322297800084e-5_wp, &
                          0.0004067771916715866_wp, 0.006179993351987645_wp, -4.4213978296108075e-5_wp, 0.0003712344432679312_wp, &
                          0.005532171216336799_wp, -4.4213978296108075e-5_wp, 0.0003933650352974713_wp, 0.005923875904780385_wp, &
                          -4.4213978296108075e-5_wp, 0.0004067771916715866_wp, 4.510381668594096e-5_wp, 0.0004080147192356906_wp, &
                          4.64364212068511e-5_wp, 0.0004067771916715866_wp, 0.00041162708748395576_wp, 0.00042676762625724614_wp, &
                        -4.5671942621595244e-5_wp, -1.5513090423496048e-5_wp, 4.6732064876971965e-5_wp, 0.00011150376207005869_wp, &
                          0.00015605704888505976_wp, 1.403987687260532e-6_wp, 2.2407965391907115e-5_wp, 0.0002433731696603681_wp, &
                          0.002837153009643337_wp, 0.014185765048216687_wp, 1.5233235600714158e-5_wp, 0.00010356654240507306_wp, &
                          0.0004804072903788798_wp, 0.0025168721615574066_wp, 0.012187925964215956_wp, 0.0002034989278731592_wp &
                          ]
   real(wp), parameter :: WANT_GAMMA_T(NCASE) = [ &
                          0.00021999999999999998_wp, 4.4e-5_wp, 2.2e-6_wp, 4.4e-5_wp, &
                          0.00010999999999999999_wp, 2.2e-5_wp, 1.1e-6_wp, 2.2e-5_wp, &
                          0.00010054081999007261_wp, 2.0402816038130124e-5_wp, 5.604797160539437e-7_wp, 2.336503563874024e-5_wp, &
                          0.00012_wp, 1.136712407102279e-5_wp, 5.072329068968142e-8_wp, 2.4e-5_wp, &
                          4.4e-5_wp, 0.00021999999999999998_wp, 4.4e-5_wp, 4.4e-5_wp, &
                          0.00021999999999999998_wp, 4.4e-5_wp, 4.4e-5_wp, 0.00021999999999999998_wp, &
                          4.4e-5_wp, 4.4e-5_wp, 1.1e-5_wp, 4.4e-5_wp, &
                          1.1e-5_wp, 4.4e-5_wp, 4.4e-5_wp, 4.4e-5_wp, &
                          4.4e-5_wp, 4.4e-5_wp, 4.4e-5_wp, 4.4e-5_wp, &
                          4.4e-5_wp, 5.072329068968142e-8_wp, 8.53774034720185e-7_wp, 9.78080811878929e-6_wp, &
                          0.00012_wp, 0.0006000000000000001_wp, 5.604797160539437e-7_wp, 4.006829008979823e-6_wp, &
                          1.9109337398948267e-5_wp, 0.00010104317131827508_wp, 0.0004862016975477157_wp, 2.2e-5_wp &
                          ]
   real(wp), parameter :: WANT_GAMMA_S(NCASE) = [ &
                          6.285714285714286e-6_wp, 1.2571428571428571e-6_wp, 6.285714285714286e-8_wp, 1.2571428571428571e-6_wp, &
                          3.1e-6_wp, 6.2e-7_wp, 3.1e-8_wp, 6.2e-7_wp, &
                          4.370397911857091e-6_wp, 8.746286446302056e-7_wp, 4.224618444966326e-8_wp, 8.794080605780259e-7_wp, &
                          3.9e-6_wp, 4.6673965420745825e-7_wp, 4.378014432526135e-9_wp, 7.8e-7_wp, &
                          1.2571428571428571e-6_wp, 6.285714285714286e-6_wp, 1.2571428571428571e-6_wp, 1.2571428571428571e-6_wp, &
                          6.285714285714286e-6_wp, 1.2571428571428571e-6_wp, 1.2571428571428571e-6_wp, 6.285714285714286e-6_wp, &
                          1.2571428571428571e-6_wp, 1.2571428571428571e-6_wp, 3.142857142857143e-7_wp, 1.2571428571428571e-6_wp, &
                          3.142857142857143e-7_wp, 1.2571428571428571e-6_wp, 1.2571428571428571e-6_wp, 1.2571428571428571e-6_wp, &
                          1.2571428571428571e-6_wp, 1.2571428571428571e-6_wp, 1.2571428571428571e-6_wp, 1.2571428571428571e-6_wp, &
                          1.2571428571428571e-6_wp, 4.378014432526135e-9_wp, 5.0738061627248945e-8_wp, 4.205988930194123e-7_wp, &
                          3.9e-6_wp, 1.95e-5_wp, 4.224618444966326e-8_wp, 2.1612779851272622e-7_wp, &
                          8.720981083978105e-7_wp, 4.37134261201938e-6_wp, 2.181979666981811e-5_wp, 6.2e-7_wp &
                          ]
   real(wp), parameter :: WANT_Q_OCEAN(NCASE) = [ &
                          2064.1177795638732_wp, 135.86358201830984_wp, 16.945827790502783_wp, -14.767468750900067_wp, &
                          1029.6594274027989_wp, 67.64255276542767_wp, 8.450078642411407_wp, -7.33766499852947_wp, &
                          1009.9408462424582_wp, 70.61353359023609_wp, 5.08790069063853_wp, -8.801877264290168_wp, &
                          1150.447057559704_wp, 38.90967318472374_wp, 0.4689318875450174_wp, -8.522117647465294_wp, &
                          135.86358201830984_wp, 2064.1177795638732_wp, -14.767468750900067_wp, 141.3054526450088_wp, &
                          2110.612622678072_wp, -14.767468750900067_wp, 137.89012807996784_wp, 2081.918245077691_wp, &
                          -14.767468750900067_wp, 135.86358201830984_wp, 15.064674773104294_wp, 136.27691622472054_wp, &
                          15.509764683088266_wp, 135.86358201830984_wp, 137.4834472196411_wp, 142.54038716992014_wp, &
                          -15.254428835612927_wp, -5.1813722014477985_wp, 15.608509668908574_wp, 37.24225653139964_wp, &
                          52.12305432760988_wp, 0.4689318875450174_wp, 7.484260440896973_wp, 81.28663866656295_wp, &
                          947.6091052208742_wp, 4738.045526104372_wp, 5.08790069063853_wp, 34.591225163294425_wp, &
                          160.45603498654583_wp, 840.6353019601743_wp, 4070.767272048129_wp, 67.96864190963521_wp &
                          ]
   real(wp), parameter :: WANT_Q_ICE(NCASE) = [ &
                          0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                          0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                          0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                          0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                          0.0_wp, 0.0_wp, 0.0_wp, 17.31314859351977_wp, &
                          262.86743642158126_wp, 0.0_wp, 6.506206290612583_wp, 103.3436928810429_wp, &
                          0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                          0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                          0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                          0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                          0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                          0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp &
                          ]
   real(wp), parameter :: WANT_Q_LATENT(NCASE) = [ &
                          2064.1177795638732_wp, 135.86358201830993_wp, 16.94582779050278_wp, -14.767468750900097_wp, &
                          1029.6594274027989_wp, 67.64255276542767_wp, 8.45007864241141_wp, -7.337664998529558_wp, &
                          1009.9408462424591_wp, 70.61353359023602_wp, 5.087900690638529_wp, -8.801877264290138_wp, &
                          1150.447057559704_wp, 38.90967318472373_wp, 0.4689318875450177_wp, -8.522117647465228_wp, &
                          135.86358201830993_wp, 2064.1177795638732_wp, -14.767468750900097_wp, 123.99230405148903_wp, &
                          1847.7451862564908_wp, -14.767468750900097_wp, 131.3839217893554_wp, 1978.5745521966485_wp, &
                          -14.767468750900097_wp, 135.86358201830993_wp, 15.064674773104281_wp, 136.27691622472065_wp, &
                          15.50976468308827_wp, 135.86358201830993_wp, 137.48344721964122_wp, 142.54038716992022_wp, &
                          -15.254428835612812_wp, -5.18137220144768_wp, 15.608509668908637_wp, 37.242256531399605_wp, &
                          52.123054327609964_wp, 0.4689318875450177_wp, 7.484260440896977_wp, 81.28663866656295_wp, &
                          947.6091052208747_wp, 4738.0455261043735_wp, 5.087900690638529_wp, 34.5912251632944_wp, &
                          160.45603498654586_wp, 840.6353019601738_wp, 4070.767272048129_wp, 67.96864190963517_wp &
                          ]
   real(wp), parameter :: WANT_T_STAR(NCASE) = [ &
                          3.21265_wp, 1.2018499999999999_wp, 2.70725_wp, -0.14815000000000023_wp, &
                          3.21265_wp, 1.2018499999999999_wp, 2.70725_wp, -0.14815000000000023_wp, &
                          3.21265_wp, 1.2018499999999999_wp, 2.70725_wp, -0.14815000000000023_wp, &
                          3.21265_wp, 1.2018499999999999_wp, 2.70725_wp, -0.14815000000000023_wp, &
                          1.2018499999999999_wp, 3.21265_wp, -0.14815000000000023_wp, 1.2018499999999999_wp, &
                          3.21265_wp, -0.14815000000000023_wp, 1.2018499999999999_wp, 3.21265_wp, &
                          -0.14815000000000023_wp, 1.2018499999999999_wp, 0.5621300000000002_wp, 1.2325_wp, &
                          0.59212_wp, 1.2018499999999999_wp, 1.2018499999999999_wp, 1.2018499999999999_wp, &
                          -0.1531149999999999_wp, -0.05145999999999984_wp, 0.15185000000000004_wp, 0.35516000000000014_wp, &
                          0.49070000000000036_wp, 2.70725_wp, 2.70725_wp, 2.70725_wp, &
                          2.70725_wp, 2.70725_wp, 2.70725_wp, 2.70725_wp, &
                          2.70725_wp, 2.70725_wp, 2.70725_wp, 1.2325_wp &
                          ]
   real(wp), parameter :: WANT_S_STAR(NCASE) = [ &
                          16.96335459257261_wp, 8.25945943164039_wp, 15.21818804828354_wp, -1.2221352114095438_wp, &
                          17.0622341519943_wp, 8.319055594423173_wp, 15.312289950431786_wp, -1.2316275871440894_wp, &
                          13.959206720720879_wp, 6.567921096207339_wp, 8.98481363756084_wp, -1.035887152830668_wp, &
                          16.035426580791736_wp, 6.740035797953233_wp, 8.227159307871382_wp, -1.133905239080626_wp, &
                          8.25945943164039_wp, 16.96335459257261_wp, -1.2221352114095438_wp, 7.698823667274759_wp, &
                          16.005350490802886_wp, -1.2221352114095438_wp, 8.050679348790649_wp, 16.596584477147918_wp, &
                          -1.2221352114095438_wp, 8.25945943164039_wp, 4.201805815187271_wp, 8.278557410031205_wp, &
                          4.310348675788703_wp, 8.25945943164039_wp, 8.09257667612944_wp, 7.571597478192206_wp, &
                          -1.2639117525655266_wp, -0.41916442861896996_wp, 1.2040074918464185_wp, 2.74024223071736_wp, &
                          3.717183468548768_wp, 8.227159307871382_wp, 10.397612848598488_wp, 12.461304075450435_wp, &
                          14.338362920432207_wp, 14.338362920432207_wp, 8.98481363756084_wp, 11.000605549190777_wp, &
                          12.071901199255187_wp, 12.421707581974374_wp, 12.181350520216078_wp, 8.337159423009602_wp &
                          ]
   real(wp), parameter :: WANT_B_FLUX(NCASE) = [ &
                          -6.353554936045811e-7_wp, -6.771016237194087e-8_wp, -5.840808727151461e-9_wp, 1.0497255720928493e-8_wp, &
                          -3.146577147377587e-7_wp, -3.362060015969593e-8_wp, -2.8947151296011303e-9_wp, 5.217440855949678e-9_wp, &
                          -3.788570688649693e-7_wp, -3.786816526752866e-8_wp, -2.464352828632766e-9_wp, 6.219960693787233e-9_wp, &
                        -3.7804050317480286e-7_wp, -2.0716158086699036e-8_wp, -2.3509120859477056e-10_wp, 6.040983189309075e-9_wp, &
                          -6.771016237194087e-8_wp, -6.353554936045811e-7_wp, 1.0497255720928493e-8_wp, -6.179963498085352e-8_wp, &
                          -5.848564901144008e-7_wp, 1.0497255720928493e-8_wp, -6.550908909615159e-8_wp, -6.160220426259157e-7_wp, &
                          1.0497255720928493e-8_wp, -6.771016237194087e-8_wp, -8.810013798660682e-9_wp, -6.785783470709607e-8_wp, &
                          -9.03258401775636e-9_wp, -6.771016237194087e-8_wp, -6.595079352635232e-8_wp, -6.045834724426584e-8_wp, &
                          1.0857685110686188e-8_wp, 3.58987834755012e-9_wp, -1.0246524885642078e-8_wp, -2.316638902914681e-8_wp, &
                        -3.128186287588762e-8_wp, -2.3509120859477056e-10_wp, -3.3881015079142593e-9_wp, -3.303918173517655e-8_wp, &
                          -3.453000739751746e-7_wp, -1.726500369875873e-6_wp, -2.464352828632766e-9_wp, -1.519193943032761e-8_wp, &
                          -6.661792952918356e-8_wp, -3.424244627489573e-7_wp, -1.6801119849115234e-6_wp, -3.368956000530119e-8_wp &
                          ]
   real(wp), parameter :: WANT_L_PLUS(NCASE) = [ &
                          20178.48739731174_wp, 302.9503961331666_wp, 0.021949893275764458_wp, -1954.112680319284_wp, &
                          40744.31428194176_wp, 610.1265419232789_wp, 0.044289376489628496_wp, -3931.586591811355_wp, &
                          33839.97257573282_wp, 541.6903715272026_wp, 0.05202385255696524_wp, -3297.901951906161_wp, &
                          33913.06675566587_wp, 990.184590548689_wp, 0.5453420779597007_wp, -3395.6096002920062_wp, &
                          302.9503961331666_wp, 20178.48739731174_wp, -1954.112680319284_wp, 331.9246225188175_wp, &
                          21920.78405081059_wp, -1954.112680319284_wp, 313.12938091251164_wp, 20811.776094671633_wp, &
                          -1954.112680319284_wp, 302.9503961331666_wp, 9.095128221068894_wp, 302.2911149664997_wp, &
                          8.871016861917715_wp, 302.9503961331666_wp, 311.0322016766045_wp, 339.28847624537156_wp, &
                          -1889.2443742572434_wp, -5714.071209911417_wp, 2001.9295070042783_wp, 885.4561013808536_wp, &
                          655.7416543319742_wp, 0.5453420779597007_wp, 23.64988325793401_wp, 620.8634547078025_wp, &
                          37128.61301452982_wp, 4641076.6268163_wp, 0.05202385255696524_wp, 5.274389454729194_wp, &
                          307.91741289159205_wp, 37440.41158038357_wp, 4769218.114495227_wp, 608.8776614949184_wp &
                          ]

contains

   subroutine collect_ocean_cavity_melt_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("cavity_melt_oracle_const_gamma", test_oracle_const_gamma), &
                  new_unittest("cavity_melt_oracle_hj99", test_oracle_hj99), &
                  new_unittest("cavity_melt_oracle_yung25", test_oracle_yung25), &
                  new_unittest("cavity_melt_residuals", test_residuals), &
                  new_unittest("cavity_melt_limit_two_equation", test_limit_two_equation), &
                  new_unittest("cavity_melt_limit_gamma_ratio_zero", test_limit_ratio_zero), &
                  new_unittest("cavity_melt_relaxation_identity", test_relaxation_identity), &
                  new_unittest("cavity_melt_sign_conventions", test_sign_conventions), &
                  new_unittest("cavity_melt_ice_pump", test_ice_pump), &
                  new_unittest("cavity_melt_ustar_law", test_ustar_law), &
                  new_unittest("cavity_melt_parse_strings", test_parse_strings), &
                  new_unittest("cavity_melt_reserved_laws_refuse", test_reserved_laws_refuse), &
                  new_unittest("cavity_melt_guards", test_guards), &
                  new_unittest("cavity_melt_do_concurrent", test_do_concurrent) &
                  ]
   end subroutine collect_ocean_cavity_melt_tests

   ! ======================================================================
   ! Helpers
   ! ======================================================================

   pure function make_eos(tfreeze_set) result(eos)
      !! EOS handle carrying the requested liquidus SET, written by the
      !! production `eos_apply_tfreeze_set` — the same call the ocean
      !! configure path makes.  No liquidus coefficient is spelled out in
      !! this test file.
      integer, intent(in) :: tfreeze_set
      type(eos_t) :: eos
      call eos_apply_tfreeze_set(eos, tfreeze_set)
   end function make_eos

   pure function make_const(const_set) result(con)
      !! The two constants bundles the embedded oracle cases use.
      !! `const_set = 1` is ISOMIP+ (the type defaults); `2` is the
      !! Jenkins, Nicholls & Corr (2010) Table 1 p. 2300 bundle, which
      !! differs from ISOMIP+ only in the two densities.
      integer, intent(in) :: const_set
      type(ocean_cavity_const_t) :: con
      if (const_set == 2) then
         con%rho_w = 1030.0_wp
         con%rho_i = 916.0_wp
      end if
   end function make_const

   subroutine check_close(error, got, want, atol, what, case_id)
      !! Relative comparison with an absolute floor — never bit-equality.
      type(error_type), allocatable, intent(out) :: error
      real(wp), intent(in) :: got, want, atol
      character(len=*), intent(in) :: what, case_id
      character(len=64) :: got_str, want_str
      write (got_str, '(es24.17)') got
      write (want_str, '(es24.17)') want
      call check(error, abs(got - want) <= RTOL*abs(want) + atol, &
                 trim(case_id)//": "//what//" = "//trim(adjustl(got_str))// &
                 ", oracle "//trim(adjustl(want_str)))
   end subroutine check_close

   subroutine run_oracle_law(error, law_code)
      !! Replay every embedded oracle case for one exchange law.
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: law_code
      type(ocean_cavity_exchange_t) :: par
      type(ocean_cavity_ice_t) :: ice
      type(ocean_cavity_const_t) :: con
      type(ocean_cavity_solution_t) :: sol
      type(eos_t) :: eos
      integer :: c, ierr, n_run

      n_run = 0
      do c = 1, NCASE
         if (CASE_LAW(c) /= law_code) cycle
         n_run = n_run + 1
         par%law = CASE_LAW(c)
         par%gamma_t_coeff = CASE_GAMMA_T_COEFF(c)
         par%gamma_s_coeff = CASE_GAMMA_S_COEFF(c)
         par%f_cor = CASE_F_COR(c)
         ice%mode = CASE_ICE_MODE(c)
         ice%T_ice = CASE_T_ICE(c)
         con = make_const(CASE_CONST_SET(c))
         if (CASE_TFREEZE_SET(c) == 2) then
            eos = make_eos(TFREEZE_SET_ISOMIP)
         else
            eos = make_eos(TFREEZE_SET_SEAICE)
         end if

         call cavity_solve_melt(CASE_T_W(c), CASE_S_W(c), CASE_P_B(c), CASE_U_STAR(c), &
                                CASE_S_I(c), par, ice, eos, con, sol, ierr)
         call check(error, ierr == CAVITY_MELT_OK, &
                    trim(CASE_ID(c))//": solver did not return OK")
         if (allocated(error)) return

         call check_close(error, sol%T_b, WANT_T_B(c), ATOL_T, "T_b", CASE_ID(c))
         if (allocated(error)) return
         call check_close(error, sol%S_b, WANT_S_B(c), ATOL_S, "S_b", CASE_ID(c))
         if (allocated(error)) return
         call check_close(error, sol%m_mass, WANT_M_MASS(c), ATOL_M, "m_mass", CASE_ID(c))
         if (allocated(error)) return
         call check_close(error, sol%gamma_t, WANT_GAMMA_T(c), ATOL_GAMMA, "gamma_t", &
                          CASE_ID(c))
         if (allocated(error)) return
         call check_close(error, sol%gamma_s, WANT_GAMMA_S(c), ATOL_GAMMA, "gamma_s", &
                          CASE_ID(c))
         if (allocated(error)) return
         call check_close(error, sol%q_ocean, WANT_Q_OCEAN(c), ATOL_Q, "q_ocean", &
                          CASE_ID(c))
         if (allocated(error)) return
         call check_close(error, sol%q_ice, WANT_Q_ICE(c), ATOL_Q, "q_ice", CASE_ID(c))
         if (allocated(error)) return
         call check_close(error, sol%q_latent, WANT_Q_LATENT(c), ATOL_Q, "q_latent", &
                          CASE_ID(c))
         if (allocated(error)) return
         call check_close(error, sol%T_star, WANT_T_STAR(c), ATOL_T, "T_star", CASE_ID(c))
         if (allocated(error)) return
         call check_close(error, sol%S_star, WANT_S_STAR(c), ATOL_S, "S_star", CASE_ID(c))
         if (allocated(error)) return
         call check_close(error, sol%b_flux, WANT_B_FLUX(c), ATOL_B, "b_flux", CASE_ID(c))
         if (allocated(error)) return
         call check_close(error, sol%l_plus, WANT_L_PLUS(c), ATOL_LP, "l_plus", CASE_ID(c))
         if (allocated(error)) return
      end do

      call check(error, n_run > 0, "no oracle cases ran for this law")
   end subroutine run_oracle_law

   ! ======================================================================
   ! (a) Oracle
   ! ======================================================================

   subroutine test_oracle_const_gamma(error)
      !! `gamma = Gamma*u*` (Jenkins et al. 2010 / ISOMIP+) over both
      !! liquidus sets, both ice modes, three ice salinities and the
      !! pressure sweep — 36 of the 48 embedded cases.
      type(error_type), allocatable, intent(out) :: error
      call run_oracle_law(error, CAVITY_LAW_CONST_GAMMA)
   end subroutine test_oracle_const_gamma

   subroutine test_oracle_hj99(error)
      !! Holland & Jenkins (1999) eqs. (14)-(18), including the outer
      !! stratification bisection and the `eta* := 1` destabilising
      !! branch on the freezing case.
      type(error_type), allocatable, intent(out) :: error
      call run_oracle_law(error, CAVITY_LAW_HJ99)
   end subroutine test_oracle_hj99

   subroutine test_oracle_yung25(error)
      !! Yung et al. (2025) StratFeedback, eqs. (7)-(8), including the
      !! `min()` caps and the freezing-branch fallback to ConstCoeff.
      type(error_type), allocatable, intent(out) :: error
      call run_oracle_law(error, CAVITY_LAW_YUNG25)
   end subroutine test_oracle_yung25

   ! ======================================================================
   ! (b) Residuals — the three equations, independent of the oracle
   ! ======================================================================

   subroutine test_residuals(error)
      !! For a sweep over far-field state, pressure, friction velocity,
      !! law and ice mode, assert that the RETURNED state satisfies
      !!   (E1) `T_b = T_f(S_b, p_b)`
      !!   (E2) `q_ocean - q_ice - q_latent = 0`
      !!   (E3) `rho_w*gamma_s*(S_w - S_b) = m_mass*(S_b - S_i)`
      !! to round-off.  This gate owes nothing to the golden numbers: a
      !! kernel whose transcendentals differ in the last bits still has
      !! the three equations right, or it does not.
      !!
      !! NORMALISATION, and why it is not just `max|flux|`.  `q_ocean` is
      !! `rho_w*c_w*gamma_t*(T_w - T_b)` and `f_turb` is
      !! `rho_w*gamma_s*(S_w - S_b)` — both DIFFERENCES.  Near neutrality
      !! (`T* -> 0`) those differences cancel to a small fraction of
      !! their operands, so a residual divided by the RESULTING flux
      !! reports the cancellation, not the solver: measured 2.8e-13 at
      !! `T* = -6.4e-3` even though every step is correctly rounded.
      !! Each residual is therefore normalised by the largest INTERMEDIATE
      !! that entered it, which is the honest conditioning scale — and on
      !! that scale the whole sweep sits at 1e-15, i.e. round-off.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_cavity_exchange_t) :: par
      type(ocean_cavity_ice_t) :: ice
      type(ocean_cavity_const_t) :: con
      type(ocean_cavity_solution_t) :: sol
      type(eos_t) :: eos
      real(wp), parameter :: T_SWEEP(4) = [-2.4_wp, -1.9_wp, -0.5_wp, 1.2_wp]
      real(wp), parameter :: S_SWEEP(3) = [33.8_wp, 34.5_wp, 34.9_wp]
      real(wp), parameter :: P_SWEEP(3) = [0.0_wp, 2.0e6_wp, 8.0e6_wp]
      real(wp), parameter :: U_SWEEP(3) = [5.0e-4_wp, 3.0e-3_wp, 2.0e-2_wp]
      integer, parameter :: LAWS(3) = [CAVITY_LAW_CONST_GAMMA, CAVITY_LAW_HJ99, &
                                       CAVITY_LAW_YUNG25]
      integer, parameter :: MODES(2) = [CAVITY_ICE_INSULATING, CAVITY_ICE_ADV_DIFF]
      real(wp) :: r_liq, r_heat, r_salt, tf, scale, f_turb, f_phase
      integer :: it, is, ip, iu, il, im, ierr

      eos = make_eos(TFREEZE_SET_ISOMIP)
      do il = 1, size(LAWS)
         do im = 1, size(MODES)
            do it = 1, size(T_SWEEP)
               do is = 1, size(S_SWEEP)
                  do ip = 1, size(P_SWEEP)
                     do iu = 1, size(U_SWEEP)
                        par%law = LAWS(il)
                        ice%mode = MODES(im)
                        ice%T_ice = -20.0_wp
                        call cavity_solve_melt(T_SWEEP(it), S_SWEEP(is), P_SWEEP(ip), &
                                               U_SWEEP(iu), 0.0_wp, par, ice, eos, con, &
                                               sol, ierr)
                        call check(error, ierr == CAVITY_MELT_OK, &
                                   "residual sweep: solver did not return OK")
                        if (allocated(error)) return

                        tf = eos_freezing_point(eos, sol%S_b, P_SWEEP(ip))
                        r_liq = (sol%T_b - tf)/max(1.0_wp, abs(tf))
                        call check(error, abs(r_liq) <= 1.0e-13_wp, &
                                   "(E1) liquidus residual too large")
                        if (allocated(error)) return

                        scale = max(abs(sol%q_ocean), abs(sol%q_ice), &
                                    abs(sol%q_latent), &
                                    con%rho_w*con%c_w*sol%gamma_t* &
                                    max(1.0_wp, abs(T_SWEEP(it)), abs(sol%T_b)), &
                                    1.0e-30_wp)
                        r_heat = (sol%q_ocean - sol%q_ice - sol%q_latent)/scale
                        call check(error, abs(r_heat) <= 1.0e-13_wp, &
                                   "(E2) heat residual too large")
                        if (allocated(error)) return

                        call cavity_salt_fluxes(con, S_SWEEP(is), sol%S_b, sol%m_mass, &
                                                sol%gamma_s, 0.0_wp, f_turb, f_phase)
                        scale = max(abs(f_turb), abs(f_phase), &
                                    con%rho_w*sol%gamma_s* &
                                    max(1.0_wp, S_SWEEP(is), abs(sol%S_b)), &
                                    1.0e-30_wp)
                        r_salt = (f_turb - f_phase)/scale
                        call check(error, abs(r_salt) <= 1.0e-13_wp, &
                                   "(E3) salt residual too large")
                        if (allocated(error)) return
                     end do
                  end do
               end do
            end do
         end do
      end do
   end subroutine test_residuals

   ! ======================================================================
   ! (c) Limits, signs, ice pump
   ! ======================================================================

   subroutine test_limit_two_equation(error)
      !! `gamma_S/gamma_T -> infinity` is the TWO-equation form: the
      !! interface salinity goes to the far-field salinity and the melt
      !! rate to `q_ocean/L_f` evaluated at `T_f(S_w, p)`.  (Note this is
      !! the limit of LARGE salt exchange, not small — the prototype
      !! found the opposite claim in the planning document.)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_cavity_ice_t) :: ice
      type(ocean_cavity_const_t) :: con
      type(eos_t) :: eos
      real(wp), parameter :: T_W = -0.5_wp, S_W = 34.6_wp, P_B = 4.0e6_wp
      real(wp), parameter :: G_T = 2.2e-4_wp
      real(wp) :: T3, S3, M3, T2, S2, M2
      integer :: ierr

      eos = make_eos(TFREEZE_SET_ISOMIP)
      call cavity_three_equation(T_W, S_W, P_B, G_T, 1.0e8_wp*G_T, 0.0_wp, ice, eos, &
                                 con, T3, S3, M3, ierr)
      call check(error, ierr == CAVITY_MELT_OK, "three-equation solve failed")
      if (allocated(error)) return
      call cavity_two_equation(T_W, S_W, P_B, G_T, ice, eos, con, T2, S2, M2, ierr)
      call check(error, ierr == CAVITY_MELT_OK, "two-equation solve failed")
      if (allocated(error)) return

      call check(error, S2 == S_W, "two-equation S_b is not exactly S_w")
      if (allocated(error)) return
      call check(error, abs(S3 - S_W) <= 1.0e-6_wp*S_W, &
                 "gamma_s/gamma_t -> infinity did not drive S_b to S_w")
      if (allocated(error)) return
      call check(error, abs(M3 - M2) <= 1.0e-6_wp*abs(M2), &
                 "gamma_s/gamma_t -> infinity did not reproduce the two-equation melt")
      if (allocated(error)) return
      call check(error, abs(T3 - T2) <= 1.0e-6_wp*abs(T2), &
                 "gamma_s/gamma_t -> infinity did not reproduce the two-equation T_b")
   end subroutine test_limit_two_equation

   subroutine test_limit_ratio_zero(error)
      !! `gamma_S/gamma_T -> 0` sends `T_b -> T_w` and the melt rate to
      !! zero: the interface becomes salty enough that its OWN freezing
      !! point is the far-field temperature, so the thermal driving
      !! across the sublayer vanishes.  The interface state at fixed
      !! ratio is independent of the MAGNITUDE of gamma — only the fluxes
      !! scale — which is why the melt rate is exactly proportional to
      !! `gamma_T` at fixed ratio (asserted here too).
      type(error_type), allocatable, intent(out) :: error
      type(ocean_cavity_ice_t) :: ice
      type(ocean_cavity_const_t) :: con
      type(eos_t) :: eos
      real(wp), parameter :: T_W = -0.5_wp, S_W = 34.6_wp, P_B = 4.0e6_wp
      real(wp), parameter :: G_T = 2.2e-4_wp
      real(wp) :: T_b, S_b, m_mass, T_b2, S_b2, m2
      integer :: ierr

      eos = make_eos(TFREEZE_SET_ISOMIP)
      call cavity_three_equation(T_W, S_W, P_B, G_T, 1.0e-10_wp*G_T, 0.0_wp, ice, eos, &
                                 con, T_b, S_b, m_mass, ierr)
      call check(error, ierr == CAVITY_MELT_OK, "three-equation solve failed")
      if (allocated(error)) return
      call check(error, abs(T_b - T_W) <= 1.0e-6_wp*max(1.0_wp, abs(T_W)), &
                 "gamma_s/gamma_t -> 0 did not drive T_b to T_w")
      if (allocated(error)) return
      call check(error, abs(m_mass) <= 1.0e-6_wp, &
                 "gamma_s/gamma_t -> 0 did not drive the melt rate to zero")
      if (allocated(error)) return

      ! m proportional to gamma_T at fixed ratio.
      call cavity_three_equation(T_W, S_W, P_B, G_T, G_T/35.0_wp, 0.0_wp, ice, eos, &
                                 con, T_b, S_b, m_mass, ierr)
      if (allocated(error)) return
      call cavity_three_equation(T_W, S_W, P_B, 4.0_wp*G_T, 4.0_wp*G_T/35.0_wp, 0.0_wp, &
                                 ice, eos, con, T_b2, S_b2, m2, ierr)
      call check(error, abs(S_b2 - S_b) <= 1.0e-12_wp*S_b, &
                 "interface state moved when only the magnitude of gamma changed")
      if (allocated(error)) return
      call check(error, abs(m2 - 4.0_wp*m_mass) <= 1.0e-12_wp*abs(m2), &
                 "melt rate is not proportional to gamma_T at fixed ratio")
   end subroutine test_limit_ratio_zero

   subroutine test_relaxation_identity(error)
      !! The shipped sea-ice basal-flux law
      !! (`rdb_ice_basal_flux`: `fb = rho*c_p*max(0, SST - T_f)*h/dt`)
      !! IS the two-equation form at `gamma_T*dt = h`.  Asserted here as
      !! an exact identity on the ocean-heat-flux side, with the sea-ice
      !! kernel's own constants (`RHO_WATER`, `SEAWATER_CP`) and its
      !! surface-pressure convention (`p = 0`).
      type(error_type), allocatable, intent(out) :: error
      type(ocean_cavity_ice_t) :: ice
      type(ocean_cavity_const_t) :: con
      type(eos_t) :: eos
      real(wp), parameter :: SST = -1.2_wp, S_SURF = 34.0_wp
      real(wp), parameter :: H_TOP = 12.5_wp, DT_THERM = 900.0_wp
      real(wp) :: gamma_t, T_b, S_b, m_mass, q_oc, q_ic, q_lat, fb_seaice, tfw
      integer :: ierr

      eos = make_eos(TFREEZE_SET_SEAICE)
      con%rho_w = RHO_WATER
      con%c_w = SEAWATER_CP
      gamma_t = H_TOP/DT_THERM

      call cavity_two_equation(SST, S_SURF, 0.0_wp, gamma_t, ice, eos, con, &
                               T_b, S_b, m_mass, ierr)
      call check(error, ierr == CAVITY_MELT_OK, "two-equation solve failed")
      if (allocated(error)) return
      call cavity_heat_fluxes(SST, T_b, m_mass, gamma_t, ice, con, q_oc, q_ic, q_lat)

      tfw = eos_freezing_point(eos, S_SURF, 0.0_wp)
      fb_seaice = RHO_WATER*SEAWATER_CP*max(0.0_wp, SST - tfw)*H_TOP/DT_THERM
      call check(error, abs(q_oc - fb_seaice) <= 1.0e-12_wp*abs(fb_seaice), &
                 "two-equation q_ocean at gamma_t*dt = h is not the sea-ice basal flux")
      if (allocated(error)) return
      call check(error, q_oc > 0.0_wp, "above-freezing surface gave a non-positive flux")
   end subroutine test_relaxation_identity

   subroutine test_sign_conventions(error)
      !! Warm water melts, cools the ocean and freshens it; supercooled
      !! water freezes.  And for the two shipped (insulating /
      !! advective-diffusive) ice modes `sign(m) = sign(T*)` EXACTLY,
      !! which is what lets the conduction branch be decided before the
      !! quadratic is solved.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_cavity_exchange_t) :: par
      type(ocean_cavity_ice_t) :: ice
      type(ocean_cavity_const_t) :: con
      type(ocean_cavity_solution_t) :: sol
      type(eos_t) :: eos
      real(wp), parameter :: S_W = 34.5_wp, P_B = 4.5e6_wp, U_STAR = 5.0e-3_wp
      real(wp) :: tf, T_w
      integer :: k, ierr

      eos = make_eos(TFREEZE_SET_ISOMIP)
      tf = eos_freezing_point(eos, S_W, P_B)

      ! Warm: melt, ocean cools (q_ocean > 0), interface fresher.
      call cavity_solve_melt(tf + 0.75_wp, S_W, P_B, U_STAR, 0.0_wp, par, ice, eos, &
                             con, sol, ierr)
      call check(error, ierr == CAVITY_MELT_OK, "warm solve failed")
      if (allocated(error)) return
      call check(error, sol%m_mass > 0.0_wp, "warm water did not melt")
      if (allocated(error)) return
      call check(error, sol%q_ocean > 0.0_wp, "melting did not draw heat from the ocean")
      if (allocated(error)) return
      call check(error, sol%S_b < S_W, "melting did not freshen the interface")
      if (allocated(error)) return
      call check(error, cavity_m_weq_from_mass(con, sol%m_mass) > 0.0_wp, &
                 "freshwater-equivalent conversion flipped the sign")
      if (allocated(error)) return

      ! Supercooled: freeze, interface saltier.
      call cavity_solve_melt(tf - 0.05_wp, S_W, P_B, U_STAR, 0.0_wp, par, ice, eos, &
                             con, sol, ierr)
      call check(error, ierr == CAVITY_MELT_OK, "supercooled solve failed")
      if (allocated(error)) return
      call check(error, sol%m_mass < 0.0_wp, "supercooled water did not freeze")
      if (allocated(error)) return
      call check(error, sol%S_b > S_W, "freezing did not salinify the interface")
      if (allocated(error)) return

      ! sign(m) = sign(T*) across the transition, insulating AND adv-diff.
      ice%mode = CAVITY_ICE_ADV_DIFF
      ice%T_ice = -25.0_wp
      do k = -3, 3
         T_w = tf + 0.2_wp*real(k, wp)
         call cavity_solve_melt(T_w, S_W, P_B, U_STAR, 0.0_wp, par, ice, eos, con, &
                                sol, ierr)
         call check(error, ierr == CAVITY_MELT_OK, "sign sweep solve failed")
         if (allocated(error)) return
         if (k == 0) then
            call check(error, abs(sol%m_mass) <= 1.0e-18_wp, &
                       "zero thermal driving did not give zero melt")
         else
            call check(error, sol%m_mass*sol%T_star > 0.0_wp, &
                       "sign(m) /= sign(T*) for an insulating/adv-diff interface")
         end if
         if (allocated(error)) return
      end do
   end subroutine test_sign_conventions

   subroutine test_ice_pump(error)
      !! THE ICE PUMP: identical water melts MORE at higher pressure,
      !! because the liquidus pressure coefficient depresses the freezing
      !! point with depth and so raises the thermal driving.  Monotone in
      !! `p_b` over the full ice-shelf range.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_cavity_exchange_t) :: par
      type(ocean_cavity_ice_t) :: ice
      type(ocean_cavity_const_t) :: con
      type(ocean_cavity_solution_t) :: sol
      type(eos_t) :: eos
      real(wp), parameter :: T_W = -1.6_wp, S_W = 34.5_wp, U_STAR = 5.0e-3_wp
      real(wp) :: prev, p_b
      integer :: k, ierr

      eos = make_eos(TFREEZE_SET_ISOMIP)
      prev = -huge(1.0_wp)
      do k = 0, 8
         p_b = 1.0e6_wp*real(k, wp)
         call cavity_solve_melt(T_W, S_W, p_b, U_STAR, 0.0_wp, par, ice, eos, con, &
                                sol, ierr)
         call check(error, ierr == CAVITY_MELT_OK, "ice-pump solve failed")
         if (allocated(error)) return
         call check(error, sol%m_mass > prev, &
                    "melt rate did not increase with interface pressure (no ice pump)")
         if (allocated(error)) return
         prev = sol%m_mass
      end do
   end subroutine test_ice_pump

   subroutine test_ustar_law(error)
      !! `u* = max(sqrt(cd*(u^2+v^2+u_tide^2)), u*_min)`, and the
      !! ISOMIP+ at-rest value: `cd = 2.5e-3`, `u_tide = 0.01 m/s` gives
      !! exactly 5e-4 m/s with no flow.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: u_star
      integer :: ierr

      call cavity_ustar(0.0_wp, 0.0_wp, CAVITY_CD_ISOMIP, CAVITY_U_TIDE_ISOMIP, &
                        0.0_wp, u_star, ierr)
      call check(error, ierr == CAVITY_MELT_OK, "ustar failed on the at-rest case")
      if (allocated(error)) return
      call check(error, abs(u_star - 5.0e-4_wp) <= 1.0e-16_wp, &
                 "ISOMIP+ at-rest friction velocity is not 5e-4 m/s")
      if (allocated(error)) return

      ! The floor engages below itself and not above it.
      call cavity_ustar(0.0_wp, 0.0_wp, CAVITY_CD_ISOMIP, 0.0_wp, &
                        CAVITY_USTAR_MIN_YUNG25, u_star, ierr)
      call check(error, u_star == CAVITY_USTAR_MIN_YUNG25, &
                 "the friction-velocity floor did not engage at rest")
      if (allocated(error)) return
      call cavity_ustar(0.4_wp, 0.3_wp, CAVITY_CD_ISOMIP, 0.0_wp, &
                        CAVITY_USTAR_MIN_YUNG25, u_star, ierr)
      call check(error, abs(u_star - sqrt(CAVITY_CD_ISOMIP*0.25_wp)) <= 1.0e-16_wp, &
                 "the floor overrode a resolved flow")
   end subroutine test_ustar_law

   subroutine test_parse_strings(error)
      !! Namelist-string parsing, including the fail-loud INVALID codes.
      type(error_type), allocatable, intent(out) :: error
      call check(error, parse_cavity_exchange_law("const_gamma") == CAVITY_LAW_CONST_GAMMA, &
                 "const_gamma did not parse")
      if (allocated(error)) return
      call check(error, parse_cavity_exchange_law("  hj99 ") == CAVITY_LAW_HJ99, &
                 "hj99 did not parse (leading/trailing blanks)")
      if (allocated(error)) return
      call check(error, parse_cavity_exchange_law("yung25") == CAVITY_LAW_YUNG25, &
                 "yung25 did not parse")
      if (allocated(error)) return
      call check(error, parse_cavity_exchange_law("mk18") == CAVITY_LAW_MK18, &
                 "a reserved law must still parse (so it can be refused as such)")
      if (allocated(error)) return
      call check(error, parse_cavity_exchange_law("hj_99") == CAVITY_LAW_INVALID, &
                 "a mistyped law must parse to INVALID, never to a default")
      if (allocated(error)) return
      call check(error, parse_cavity_ice_mode("insulating") == CAVITY_ICE_INSULATING, &
                 "insulating did not parse")
      if (allocated(error)) return
      call check(error, parse_cavity_ice_mode("adv_diff") == CAVITY_ICE_ADV_DIFF, &
                 "adv_diff did not parse")
      if (allocated(error)) return
      call check(error, parse_cavity_ice_mode("conductive") == CAVITY_ICE_INVALID, &
                 "a mistyped ice mode must parse to INVALID")
   end subroutine test_parse_strings

   subroutine test_reserved_laws_refuse(error)
      !! Every reserved law and the reserved ice mode return
      !! NOT_IMPLEMENTED — not a wrong answer, not a silent fallback to
      !! the default law — and an INVALID code returns LAW_INVALID.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_cavity_exchange_t) :: par
      type(ocean_cavity_ice_t) :: ice
      type(ocean_cavity_const_t) :: con
      type(ocean_cavity_solution_t) :: sol
      type(eos_t) :: eos
      integer, parameter :: RESERVED(6) = [CAVITY_LAW_JENKINS91, CAVITY_LAW_ROSEVEAR22, &
                                           CAVITY_LAW_VT19, CAVITY_LAW_MK18, &
                                           CAVITY_LAW_BURCHARD22, CAVITY_LAW_JENKINS21]
      integer :: k, ierr

      eos = make_eos(TFREEZE_SET_ISOMIP)
      do k = 1, size(RESERVED)
         par%law = RESERVED(k)
         call cavity_solve_melt(-0.5_wp, 34.5_wp, 4.0e6_wp, 5.0e-3_wp, 0.0_wp, par, &
                                ice, eos, con, sol, ierr)
         call check(error, ierr == CAVITY_MELT_NOT_IMPLEMENTED, &
                    "a reserved exchange law did not report NOT_IMPLEMENTED")
         if (allocated(error)) return
         call check(error, sol%m_mass == 0.0_wp, &
                    "a refused law returned a non-zero melt rate")
         if (allocated(error)) return
      end do

      par%law = CAVITY_LAW_INVALID
      call cavity_solve_melt(-0.5_wp, 34.5_wp, 4.0e6_wp, 5.0e-3_wp, 0.0_wp, par, ice, &
                             eos, con, sol, ierr)
      call check(error, ierr == CAVITY_MELT_LAW_INVALID, &
                 "an invalid law code did not report LAW_INVALID")
      if (allocated(error)) return

      par%law = CAVITY_LAW_CONST_GAMMA
      ice%mode = CAVITY_ICE_DIFFUSIVE
      call cavity_solve_melt(-0.5_wp, 34.5_wp, 4.0e6_wp, 5.0e-3_wp, 0.0_wp, par, ice, &
                             eos, con, sol, ierr)
      call check(error, ierr == CAVITY_MELT_NOT_IMPLEMENTED, &
                 "the reserved DIFFUSIVE ice mode did not report NOT_IMPLEMENTED")
      if (allocated(error)) return
      call check(error, sol%m_mass == 0.0_wp, &
                 "a refused ice mode returned a non-zero melt rate")
   end subroutine test_reserved_laws_refuse

   ! ======================================================================
   ! (d) Guards
   ! ======================================================================

   subroutine test_guards(error)
      !! Non-finite and out-of-domain inputs return the documented status
      !! AND the safe state — `m_mass` EXACTLY zero, never a
      !! finite-but-plausible melt rate.  The friction-velocity guard is
      !! the pointed one: under nvfortran's relaxed FP a NaN reaching the
      !! `max(., u*_min)` clamp would come back out as `u*_min`, i.e. as
      !! a small, plausible, permanently wrong melt rate.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_cavity_exchange_t) :: par
      type(ocean_cavity_ice_t) :: ice
      type(ocean_cavity_const_t) :: con
      type(ocean_cavity_solution_t) :: sol
      type(eos_t) :: eos
      real(wp) :: nan, inf, u_star, T_b, S_b, m_mass, gamma_t, gamma_s
      integer :: ierr

      nan = ieee_value(1.0_wp, ieee_quiet_nan)
      inf = ieee_value(1.0_wp, ieee_positive_inf)
      eos = make_eos(TFREEZE_SET_ISOMIP)

      ! --- u* : NaN must NOT be laundered into the floor ---
      call cavity_ustar(nan, 0.0_wp, CAVITY_CD_ISOMIP, CAVITY_U_TIDE_ISOMIP, &
                        CAVITY_USTAR_MIN_YUNG25, u_star, ierr)
      call check(error, ierr == CAVITY_MELT_NONFINITE_INPUT, &
                 "a NaN velocity did not report NONFINITE_INPUT")
      if (allocated(error)) return
      call check(error, u_star == 0.0_wp, &
                 "a NaN velocity was laundered into a plausible friction velocity")
      if (allocated(error)) return
      call cavity_ustar(0.0_wp, 0.0_wp, -1.0e-3_wp, 0.0_wp, 0.0_wp, u_star, ierr)
      call check(error, ierr == CAVITY_MELT_BAD_INPUT .and. u_star == 0.0_wp, &
                 "a negative drag coefficient was accepted")
      if (allocated(error)) return

      ! --- solve_melt: non-finite far field, non-finite pressure ---
      call cavity_solve_melt(nan, 34.5_wp, 4.0e6_wp, 5.0e-3_wp, 0.0_wp, par, ice, eos, &
                             con, sol, ierr)
      call check(error, ierr == CAVITY_MELT_NONFINITE_INPUT .and. sol%m_mass == 0.0_wp, &
                 "a NaN far-field temperature did not return the safe state")
      if (allocated(error)) return
      call cavity_solve_melt(-0.5_wp, 34.5_wp, inf, 5.0e-3_wp, 0.0_wp, par, ice, eos, &
                             con, sol, ierr)
      call check(error, ierr == CAVITY_MELT_NONFINITE_INPUT .and. sol%m_mass == 0.0_wp, &
                 "an infinite interface pressure did not return the safe state")
      if (allocated(error)) return

      ! --- zero and negative friction velocity ---
      call cavity_solve_melt(-0.5_wp, 34.5_wp, 4.0e6_wp, 0.0_wp, 0.0_wp, par, ice, eos, &
                             con, sol, ierr)
      call check(error, ierr == CAVITY_MELT_BAD_INPUT .and. sol%m_mass == 0.0_wp, &
                 "u* = 0 did not return BAD_INPUT and zero melt")
      if (allocated(error)) return

      ! --- salinity domain ---
      call cavity_solve_melt(-0.5_wp, -1.0_wp, 4.0e6_wp, 5.0e-3_wp, 0.0_wp, par, ice, &
                             eos, con, sol, ierr)
      call check(error, ierr == CAVITY_MELT_BAD_INPUT .and. sol%m_mass == 0.0_wp, &
                 "a negative far-field salinity was accepted")
      if (allocated(error)) return
      call cavity_solve_melt(-0.5_wp, 34.5_wp, 4.0e6_wp, 5.0e-3_wp, -1.0_wp, par, ice, &
                             eos, con, sol, ierr)
      call check(error, ierr == CAVITY_MELT_BAD_INPUT .and. sol%m_mass == 0.0_wp, &
                 "a negative ice salinity was accepted")
      if (allocated(error)) return
      call cavity_solve_melt(-0.5_wp, 34.5_wp, 4.0e6_wp, 5.0e-3_wp, 40.0_wp, par, ice, &
                             eos, con, sol, ierr)
      call check(error, ierr == CAVITY_MELT_BAD_INPUT .and. sol%m_mass == 0.0_wp, &
                 "S_i >= S_w was accepted (the root bracketing needs S_w > S_i)")
      if (allocated(error)) return

      ! --- gamma_s = 0 is not a limit of the three-equation form ---
      call cavity_three_equation(-0.5_wp, 34.5_wp, 4.0e6_wp, 1.0e-4_wp, 0.0_wp, 0.0_wp, &
                                 ice, eos, con, T_b, S_b, m_mass, ierr)
      call check(error, ierr == CAVITY_MELT_BAD_INPUT .and. m_mass == 0.0_wp, &
                 "gamma_s = 0 was accepted by the three-equation form")
      if (allocated(error)) return
      call check(error, S_b == 34.5_wp, "the safe state did not set S_b = S_w")
      if (allocated(error)) return
      call check(error, T_b == eos_freezing_point(eos, 34.5_wp, 4.0e6_wp), &
                 "the safe state did not set T_b = T_f(S_w, p_b)")
      if (allocated(error)) return

      ! --- hj99 on the equator ---
      par%law = CAVITY_LAW_HJ99
      par%f_cor = 0.0_wp
      call cavity_exchange_velocities(par, con, 5.0e-3_wp, CAVITY_L_PLUS_NEUTRAL, &
                                      -0.5_wp, 34.5_wp, -2.0_wp, 34.0_wp, &
                                      gamma_t, gamma_s, ierr)
      call check(error, ierr == CAVITY_MELT_NO_CORIOLIS, &
                 "hj99 at f = 0 did not report NO_CORIOLIS")
      if (allocated(error)) return
      call check(error, gamma_t == 0.0_wp .and. gamma_s == 0.0_wp, &
                 "hj99 at f = 0 returned non-zero exchange velocities")
      if (allocated(error)) return
      call cavity_solve_melt(-0.5_wp, 34.5_wp, 4.0e6_wp, 5.0e-3_wp, 0.0_wp, par, ice, &
                             eos, con, sol, ierr)
      call check(error, ierr == CAVITY_MELT_NO_CORIOLIS .and. sol%m_mass == 0.0_wp, &
                 "hj99 at f = 0 did not propagate NO_CORIOLIS with zero melt")
      if (allocated(error)) return

      ! --- a NaN exchange coefficient must not reach the multiply ---
      par%law = CAVITY_LAW_CONST_GAMMA
      par%gamma_t_coeff = nan
      call cavity_exchange_velocities(par, con, 5.0e-3_wp, CAVITY_L_PLUS_NEUTRAL, &
                                      -0.5_wp, 34.5_wp, -2.0_wp, 34.0_wp, &
                                      gamma_t, gamma_s, ierr)
      call check(error, ierr == CAVITY_MELT_NONFINITE_INPUT .and. gamma_t == 0.0_wp, &
                 "a NaN Gamma_T was not caught")
   end subroutine test_guards

   ! ======================================================================
   ! (e) Device path
   ! ======================================================================

   subroutine test_do_concurrent(error)
      !! Run the whole solver over an array of columns through
      !! `cavity_melt_columns` — the module's own `do concurrent` driver
      !! — with the `mem:separate` data directives the GPU build needs
      !! (inert comments on gfortran/ifx).  This is what makes the
      !! `!$acc routine seq` chain — `cavity_melt_point` ->
      !! `cavity_solve_melt` -> `cavity_state_at_x` ->
      !! `cavity_exchange_velocities` -> `cavity_three_equation` ->
      !! `eos_freezing_point` — actually compile and run on device,
      !! including the outer bisection loop.  The answers must match the
      !! scalar host path.
      !!
      !! The loop deliberately lives in the MODULE, not here: on the GPU
      !! build nvlink cannot resolve a `routine seq` device symbol out of
      !! `librdb_core.so` into a `do concurrent` compiled in this file
      !! (measured — `nvlink error: Undefined reference to
      !! 'rdb_ocean_cavity_melt_cavity_solve_melt_'`, GPU toolchain only;
      !! gfortran linked the same source without complaint).  Writing the
      !! loop here would therefore be a test that cannot exist on the
      !! toolchain it is meant to gate.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NCOL = 6
      type(ocean_cavity_exchange_t) :: par
      type(ocean_cavity_ice_t) :: ice
      type(ocean_cavity_const_t) :: con
      type(ocean_cavity_solution_t) :: sol
      type(eos_t) :: eos
      real(wp) :: t_w(NCOL), s_w(NCOL), p_b(NCOL), u_s(NCOL), s_i(NCOL)
      real(wp) :: t_b(NCOL), s_b(NCOL), m_mass(NCOL), q_oc(NCOL)
      integer :: status(NCOL)
      integer :: i, ierr

      eos = make_eos(TFREEZE_SET_ISOMIP)
      par%law = CAVITY_LAW_YUNG25
      ice%mode = CAVITY_ICE_ADV_DIFF
      ice%T_ice = -25.0_wp
      do i = 1, NCOL
         t_w(i) = -2.0_wp + 0.5_wp*real(i - 1, wp)
         s_w(i) = 34.0_wp + 0.1_wp*real(i - 1, wp)
         p_b(i) = 1.0e6_wp*real(i, wp)
         u_s(i) = 1.0e-3_wp*real(i, wp)
         s_i(i) = 0.0_wp
         t_b(i) = 0.0_wp
         s_b(i) = 0.0_wp
         m_mass(i) = 0.0_wp
         q_oc(i) = 0.0_wp
         status(i) = -1
      end do

      !$acc enter data copyin(t_w, s_w, p_b, u_s, s_i) &
      !$acc            create(t_b, s_b, m_mass, q_oc, status)
      call cavity_melt_columns(NCOL, t_w, s_w, p_b, u_s, s_i, par, ice, eos, con, &
                               t_b, s_b, m_mass, q_oc, status)
      !$acc update self(t_b, s_b, m_mass, q_oc, status)
      !$acc exit data delete(t_w, s_w, p_b, u_s, s_i, t_b, s_b, m_mass, q_oc, status)

      do i = 1, NCOL
         call check(error, status(i) == CAVITY_MELT_OK, &
                    "a do-concurrent column did not solve")
         if (allocated(error)) return
         call cavity_solve_melt(t_w(i), s_w(i), p_b(i), u_s(i), 0.0_wp, par, ice, eos, &
                                con, sol, ierr)
         call check_close(error, t_b(i), sol%T_b, ATOL_T, "T_b", "do_concurrent")
         if (allocated(error)) return
         call check_close(error, s_b(i), sol%S_b, ATOL_S, "S_b", "do_concurrent")
         if (allocated(error)) return
         call check_close(error, m_mass(i), sol%m_mass, ATOL_M, "m_mass", "do_concurrent")
         if (allocated(error)) return
         call check_close(error, q_oc(i), sol%q_ocean, ATOL_Q, "q_ocean", "do_concurrent")
         if (allocated(error)) return
      end do
   end subroutine test_do_concurrent

end module test_ocean_cavity_melt
