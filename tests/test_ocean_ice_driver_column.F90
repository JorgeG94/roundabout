!! Analytic tests for the driver-facing ice column + melt-side coupling
!! (sea-ice PR 3c): `ice_atm_forcing_restoring` (`rdb_ice_atm_forcing`),
!! `ice_compute_basal_flux` (`rdb_ice_basal_flux`), `ice_thermo_driver_step`
!! (`rdb_ice_thermo_driver`), and `ice_ocean_heat_flux`
!! (`rdb_ice_ocean_coupler`).
!!
!! Harness mirrors `test_ocean_ice_coupling.F90`: grid 6x4, NGHOST=2, NZ=3,
!! H_LAYER=10, S_INIT=35; plus an `ocean_surface_flux_t` seeded to a zero
!! background. `dt_therm = 1200 s`. Ice seeded directly (`m_ice`/`enth_ice`
!! via `ice_enth_from_ts`, `sal_ice = ICE_BULK_SALINITY`) rather than via the
!! frazil path, so each case can isolate grow vs. melt.
!!
!! **GPU mem:separate discipline** (per PR-3a/3b + `test_open_boundary_out_
!! closes`): one shared enter_data/exit_data span around the whole kernel
!! chain (`run_chain`); `!$acc update self` the touched arrays before any
!! host assertion; every new PR-3c seam + scratch field on `ice` is mapped
!! by `ice%enter_data()` (see `rdb_ice_state.F90`), so a single `ice%
!! enter_data()`/`exit_data()` pair covers them all. Directives are inert
!! no-ops on host builds — write them unconditionally; a green multicore
!! run proves nothing about device data motion.
!!
!! Cases (SPEC_ice-pr3c-driver-column.md §8):
!!   * `column_runs_and_melts` — end-to-end liveness: warm restoring melts
!!     a seeded cold ice cap (mass thins, m_melt_diag > 0).
!!   * `melt_freshens_ocean` — the SIGN gate: melt-dominated window through
!!     the REAL apply kernel, Q_salt < 0, S(nz) strictly decreases, salt
!!     increment identity to 1e-12 rel.
!!   * `heat_budget_closes` — the Q_heat/fb arbiter: grow window (Q_heat ~
!!     0) and melt window, ocean heat-content change == column's reported
!!     net exchange to round-off, through the real apply.
!!   * `pure_melt_conserves` — seed ice, melt fully in one window:
!!     the heat + salt the ocean gets back equals exactly what the ice
!!     held (closed cycle, frazil bank not involved).
!!   * `icefree_ocean_untouched` — the OPEN-OCEAN regression gate: ice
!!     ENABLED, WARM ice-free surface (SST >> T_f, no ice, no bank).
!!     Proves the ice-presence gate on `fb` (rdb_ice_basal_flux) so a
!!     warm open-ocean cell is NOT spuriously cooled by `Q_heat = -fb`;
!!     asserts fb=0, Q_heat/Q_salt all-zero, hT/hS bit-identical.
!!   * `disabled_bitident` — ice off ⇒ the thermo-cadence ice block never
!!     fires (byte-identical); ice on + benign restoring + no ice ⇒ a
!!     no-op window (Q_heat=Q_salt=0, tracers unchanged to round-off).
!!   * `compose_frazil_and_melt_salt` — the COMPOSE headline: a cell with
!!     a live frazil bank AND a seeded ice cap under warm restoring, so
!!     `m_frozen_diag > 0` (frazil) and `h2o_ice_to_ocn > 0` (melt) are
!!     BOTH live in the same window. Pins
!!     `salt_flux_diag*dt == (m_frozen_diag+m_net)*(s_surf-bulk)` — the
!!     gate that fails if `salt_flux_diag = ...` (overwrite) replaces
!!     `salt_flux_diag = salt_flux_diag + ...` (compose) at
!!     `rdb_ice_thermo_driver.F90:203`.
!!   * `compose_frazil_ice_melts_same_window` — the ORDER gate: an
!!     ice-free cell with a live bank nucleates ice via frazil uptake,
!!     and that same freshly-nucleated ice melts in the SAME window
!!     (`h2o_ice_to_ocn > 0`) — provable only because the mandated driver
!!     order runs frazil uptake BEFORE the column driver.
!!   * `compose_heat_not_double_counted` — with a live bank AND live melt
!!     in one window, the frazil latent heat contributes exactly zero to
!!     `Q_heat` (the frazil kernel never touches `heat_flux_diag`).
!!   * `compose_mass_closes` — the two-term mass identity
!!     `Delta(m_ice+m_snow) == m_frozen_diag - m_melt_diag` with BOTH
!!     terms strictly live, and the bank fully spent
!!     (`frazil_heat == 0` after).
module test_ocean_ice_driver_column
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_eos, only: eos_t, eos_freezing_point
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t, ocean_surface_flux_apply_tracers, &
                                     ocean_surface_flux_assemble, SEAWATER_CP
   use rdb_ice_state, only: ocean_sea_ice_t
   use rdb_ice_frazil_uptake, only: ice_frazil_uptake
   use rdb_ice_atm_forcing, only: ice_atm_forcing_restoring
   use rdb_ice_basal_flux, only: ice_compute_basal_flux
   use rdb_ice_thermo_driver, only: ice_thermo_driver_step
   use rdb_ice_ocean_coupler, only: ice_ocean_brine_flux, ice_ocean_heat_flux, ice_ocean_sw_flux
   use rdb_ice_optics, only: ICE_OPT_DEP_ICE, ICE_PEN_ICE
   use rdb_ice_enthalpy, only: ice_enth_from_ts, ice_enthalpy_liquid_freeze
   use rdb_ice_column, only: ICE_BULK_SALINITY, ICE_RHO_ICE
   implicit none
   private

   public :: collect_ocean_ice_driver_column_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 3
   real(wp), parameter :: H_LAYER = 10.0_wp
      !! Uniform layer thickness (m).
   real(wp), parameter :: S_INIT = 35.0_wp
      !! Uniform salinity (PSU).
   real(wp), parameter :: T_DEEP = 1.0_wp
      !! Sub-surface IC (degC) — untouched by every case.
   real(wp), parameter :: DT_THERM = 1200.0_wp
      !! Thermo timestep (s).
   real(wp), parameter :: H_ICE0 = 0.5_wp
      !! Seeded ice cap thickness (m) — cases 1-4.
   real(wp), parameter :: T_ICE_COLD = -10.0_wp
      !! Seeded ice temperature (degC), well below freezing.
   real(wp), parameter :: E_BANK_COMPOSE = 3.34e6_wp
      !! Analytic frazil bank (J/m² of cell) for the compose cases — same
      !! value as `test_ocean_ice_coupling.F90`'s `E_BANK`.

contains

   subroutine collect_ocean_ice_driver_column_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("column_runs_and_melts", test_column_runs_and_melts), &
                  new_unittest("melt_freshens_ocean", test_melt_freshens_ocean), &
                  new_unittest("heat_budget_closes", test_heat_budget_closes), &
                  new_unittest("pure_melt_conserves", test_pure_melt_conserves), &
                  new_unittest("icefree_ocean_untouched", test_icefree_ocean_untouched), &
                  new_unittest("disabled_bitident", test_disabled_bitident), &
                  new_unittest("compose_frazil_and_melt_salt", test_compose_frazil_and_melt_salt), &
                  new_unittest("compose_frazil_ice_melts_same_window", &
                               test_compose_frazil_ice_melts_same_window), &
                  new_unittest("compose_heat_not_double_counted", &
                               test_compose_heat_not_double_counted), &
                  new_unittest("compose_mass_closes", test_compose_mass_closes), &
                  new_unittest("sw_thru_reaches_ocean", test_sw_thru_reaches_ocean), &
                  new_unittest("sw_thru_beers_law", test_sw_thru_beers_law), &
                  new_unittest("thick_ice_nearly_opaque", test_thick_ice_nearly_opaque), &
                  new_unittest("sw_default_bitident", test_sw_default_bitident), &
                  new_unittest("sw_components_q_sw", test_sw_components_q_sw) &
                  ]
   end subroutine collect_ocean_ice_driver_column_tests

   ! -----------------------------------------------------------------
   ! Shared setup / run / teardown
   ! -----------------------------------------------------------------

   subroutine setup_state(grid, ms, eos, ice, sf, t_surface)
      !! Tiny ocean state: uniform h + S, sub-surface layers at T_DEEP,
      !! surface layer at `t_surface`. Ice slot enabled the driver way
      !! (enable latched before init), single category (v1 scope).
      !! Surface-flux slot seeded to a zero background (tests read Q_heat/
      !! Q_salt via the couplers, not the constant).
      type(hgrid_t), intent(out) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(eos_t), intent(out) :: eos
      type(ocean_sea_ice_t), intent(inout) :: ice
      type(ocean_surface_flux_t), intent(inout) :: sf
      real(wp), intent(in) :: t_surface

      call grid%init(6, 4, NGHOST, 1.0_wp, 1.0_wp)
      ms%nz_ml = NZ
      call ms%init(grid)
      call eos%init(grid)
      ice%enable = .true.
      ice%ncat = 1
      call ice%init(grid)
      call sf%init(grid)
      call sf%set_surface_flux_const(0.0_wp, 0.0_wp)

      ms%h_layer = H_LAYER
      ms%tracers(ms%idx_salinity)%hTr = S_INIT*H_LAYER
      ms%tracers(ms%idx_temperature)%hTr = T_DEEP*H_LAYER
      ms%tracers(ms%idx_temperature)%hTr(:, :, NZ) = t_surface*H_LAYER
   end subroutine setup_state

   subroutine seed_ice_cap(ice, ip, jp, h_ice, t_ice)
      !! Seed a snow-free category-1 ice cap of thickness `h_ice` (m),
      !! uniform temperature `t_ice`, bulk salinity ICE_BULK_SALINITY,
      !! at a single interior cell. Mirrors `make_initial_column`
      !! (`test_ocean_ice_column.F90`) but isothermal (no linear-profile
      !! target needed — these tests probe the driver wiring + budget
      !! identities, not the Stefan growth-rate physics).
      type(ocean_sea_ice_t), intent(inout) :: ice
      integer, intent(in) :: ip, jp
      real(wp), intent(in) :: h_ice, t_ice
      integer :: k

      ice%m_ice(ip, jp, 1) = h_ice*ICE_RHO_ICE
      ice%m_snow(ip, jp, 1) = 0.0_wp
      ice%enth_snow(ip, jp, 1, 1) = ice_enth_from_ts(0.0_wp, 0.0_wp)
      do k = 1, ice%nk_ice
         ice%enth_ice(ip, jp, 1, k) = ice_enth_from_ts(t_ice, ICE_BULK_SALINITY)
         ice%sal_ice(ip, jp, 1, k) = ICE_BULK_SALINITY
      end do
   end subroutine seed_ice_cap

   subroutine run_chain(grid, eos, ms, ice, sf, air_temp, restore_lambda, sw_down, use_components)
      !! Full PR-3c/PR-31 thermo-cadence chain, mandated order (SPEC §7):
      !! forcing -> basal -> frazil uptake -> column driver -> brine
      !! coupler -> heat coupler -> shortwave coupler -> (assembler, if
      !! components) -> the REAL apply_tracers kernel. One shared map/unmap
      !! spanning every kernel (GPU mem:separate discipline — every array
      !! touched anywhere in the chain must be device-present for the whole
      !! span). `ice_frazil_uptake` always runs (mandated order); on most
      !! cases in this file the bank is empty (ice is seeded straight onto
      !! `m_ice`/`enth_ice`/`sal_ice`), so it is a pure pass-through — but
      !! that is the DEFAULT, not an invariant: the `compose_*` cases below
      !! seed `ice%frazil_heat` deliberately, exercising the frazil term
      !! and the column driver's net-melt term in the SAME window (the
      !! mandated-order contract, `rdb_driver.F90` +
      !! `rdb_ice_thermo_driver.F90` module docstring).
      !!
      !! PR 31: `use_components` (optional, default `.false.`) selects the
      !! PR-12 component surface-flux path — `set_components` is called
      !! BEFORE `enter_data` (so `q_sw` etc. get mapped) and
      !! `ocean_surface_flux_assemble` runs after the couplers, so the
      !! components-on `q_sw` -> `Q_heat` assembly is exercised end-to-end.
      type(hgrid_t), intent(in) :: grid
      type(eos_t), intent(in) :: eos
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_sea_ice_t), intent(inout) :: ice
      type(ocean_surface_flux_t), intent(inout) :: sf
      real(wp), intent(in) :: air_temp, restore_lambda, sw_down
      logical, intent(in), optional :: use_components
      logical :: comps

      comps = .false.
      if (present(use_components)) comps = use_components
      if (comps) call sf%set_components(grid, .true.)

      !$acc enter data copyin(ms, sf)
      call ms%enter_data()
      call ice%enter_data()
      call sf%enter_data()

      call ice_atm_forcing_restoring(ice, air_temp, restore_lambda, sw_down, 0.0_wp)
      call ice_compute_basal_flux(grid, eos, ms, ice, DT_THERM)
      call ice_frazil_uptake(grid, eos, ms, ice, DT_THERM)
      call ice_thermo_driver_step(grid, eos, ms, ice, DT_THERM)
      call ice_ocean_brine_flux(sf, ice)
      call ice_ocean_heat_flux(sf, ice)
      call ice_ocean_sw_flux(sf, ice)
      if (comps) call ocean_surface_flux_assemble(grid, sf, ms)
      call ocean_surface_flux_apply_tracers(grid, sf, ms, DT_THERM)

      associate (hT => ms%tracers(ms%idx_temperature)%hTr, &
                 hS => ms%tracers(ms%idx_salinity)%hTr, &
                 mi => ice%m_ice, ei => ice%enth_ice, si => ice%sal_ice, &
                 mf => ice%m_frozen_diag, sd => ice%salt_flux_diag, &
                 hf => ice%heat_flux_diag, md => ice%m_melt_diag, &
                 fb => ice%fb, h2o => ice%h2o_ocn_to_ice, h2i => ice%h2o_ice_to_ocn, &
                 hto => ice%heat_to_ocn, qs => sf%Q_salt, qh => sf%Q_heat, &
                 fz => ice%frazil_heat, ssurf => ice%ssurf_seam, &
                 swd => ice%sw_thru_diag, swc => ice%sw_thru)
         !$acc update self(hT, hS, mi, ei, si, mf, sd, hf, md, fb, h2o, h2i, hto, qs, qh, &
         !$acc              fz, ssurf, swd, swc)
      end associate
      if (comps) then
         !$acc update self(sf%q_sw)
      end if

      call sf%exit_data()
      call ice%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, sf)
   end subroutine run_chain

   subroutine teardown(ms, eos, ice, sf)
      type(multilayer_state_t), intent(inout) :: ms
      type(eos_t), intent(inout) :: eos
      type(ocean_sea_ice_t), intent(inout) :: ice
      type(ocean_surface_flux_t), intent(inout) :: sf
      call sf%destroy()
      call ice%destroy()
      call eos%destroy()
      call ms%destroy()
   end subroutine teardown

   ! -----------------------------------------------------------------
   ! Cases
   ! -----------------------------------------------------------------

   subroutine test_column_runs_and_melts(error)
      !! Case 1: ice enabled, warm restoring (air_temp=+5, lambda=20,
      !! sw_down=0), seed a cold ice cap. Run the full seam chain. Assert
      !! the cap THINS and produces meltwater (m_melt_diag > 0).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      type(ocean_surface_flux_t) :: sf
      real(wp) :: m_ice_before
      integer :: ip, jp
      checks: block

         call setup_state(grid, ms, eos, ice, sf, 0.0_wp)
         ip = grid%nx_total/2
         jp = grid%ny_total/2
         call seed_ice_cap(ice, ip, jp, H_ICE0, T_ICE_COLD)
         m_ice_before = ice%m_ice(ip, jp, 1)

         call run_chain(grid, eos, ms, ice, sf, air_temp=5.0_wp, restore_lambda=20.0_wp, &
                        sw_down=0.0_wp)

         call check(error, ice%m_ice(ip, jp, 1) < m_ice_before, &
                    "ice cap must thin under warm restoring")
         if (allocated(error)) exit checks
         call check(error, ice%m_melt_diag(ip, jp) > 0.0_wp, &
                    "m_melt_diag must be strictly positive (net melt)")

      end block checks
      call teardown(ms, eos, ice, sf)
   end subroutine test_column_runs_and_melts

   subroutine test_melt_freshens_ocean(error)
      !! Case 2 — the SIGN gate. A melt-dominated window (warm restoring),
      !! seeded ice, s_surf = S_INIT = 35 > ICE_BULK_SALINITY = 4. Assert
      !! salt_flux_diag < 0, Q_salt < 0, has_salt, S(nz) strictly
      !! DECREASES, and the salt increment matches
      !! rho0*Delta(hS_top) == m_net*(s_surf - ICE_BULK_SALINITY) to
      !! 1e-12 rel (m_net = h2o_ocn_to_ice - h2o_ice_to_ocn < 0 here).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      type(ocean_surface_flux_t) :: sf
      real(wp) :: s_before, s_after, m_net, expected_dhS, tol
      integer :: ip, jp
      checks: block

         call setup_state(grid, ms, eos, ice, sf, 0.0_wp)
         ip = grid%nx_total/2
         jp = grid%ny_total/2
         call seed_ice_cap(ice, ip, jp, H_ICE0, T_ICE_COLD)
         s_before = ms%tracers(ms%idx_salinity)%hTr(ip, jp, NZ)/H_LAYER

         call run_chain(grid, eos, ms, ice, sf, air_temp=10.0_wp, restore_lambda=30.0_wp, &
                        sw_down=0.0_wp)

         s_after = ms%tracers(ms%idx_salinity)%hTr(ip, jp, NZ)/H_LAYER

         call check(error, ice%salt_flux_diag(ip, jp) < 0.0_wp, &
                    "salt_flux_diag must be strictly negative (net melt freshens)")
         if (allocated(error)) exit checks
         call check(error, sf%has_salt, "has_salt must be latched true")
         if (allocated(error)) exit checks
         call check(error, sf%Q_salt(ip, jp) < 0.0_wp, "Q_salt must be strictly negative")
         if (allocated(error)) exit checks
         call check(error, s_after < s_before, "S(nz) must strictly decrease (the trap gate)")
         if (allocated(error)) exit checks

         m_net = ice%h2o_ocn_to_ice(ip, jp, 1) - ice%h2o_ice_to_ocn(ip, jp, 1)
         expected_dhS = m_net*(S_INIT - ICE_BULK_SALINITY)
         tol = 1.0e-12_wp*abs(expected_dhS)
         call check(error, abs(sf%rho0*(ms%tracers(ms%idx_salinity)%hTr(ip, jp, NZ) &
                                        - s_before*H_LAYER) - expected_dhS) <= tol, &
                    "salt removed from the ocean must equal m_net*(s_surf-bulk) to 1e-12 rel")

      end block checks
      call teardown(ms, eos, ice, sf)
   end subroutine test_melt_freshens_ocean

   subroutine test_heat_budget_closes(error)
      !! Case 3 — the Q_heat/fb arbiter. Two independent sub-blocks (own
      !! state each): a GROW window (cold restoring, no pre-existing ice
      !! melt, fb=0 since SST=T_f exactly) must give Q_heat ~ 0 (ocean
      !! unchanged); a MELT window (warm restoring, seeded ice) must
      !! satisfy, through the REAL apply kernel,
      !!   |rho0*cp*Delta(hT_top) - (heat_to_ocn - fb*dt_therm)| <= tol
      !! with tol = 1e-9*scale, scale = rho0*cp*|hT_top|.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      type(ocean_surface_flux_t) :: sf
      real(wp) :: t_f, hT_before, hT_after, d_ocean_heat, rhs, tol, scale
      integer :: ip, jp
      checks: block

         ! ---- Grow sub-block: surface EXACTLY at t_f, no ice, cold
         ! restoring => the column has no ice to touch (m_ice=0 skips
         ! ice_thermo_columns' gate), fb=0 (SST==T_f), Q_heat must be ~0. ----
         t_f = eos_freezing_point(eos, S_INIT, 0.0_wp)
         call setup_state(grid, ms, eos, ice, sf, t_f)
         ip = grid%nx_total/2
         jp = grid%ny_total/2
         hT_before = ms%tracers(ms%idx_temperature)%hTr(ip, jp, NZ)

         call run_chain(grid, eos, ms, ice, sf, air_temp=-20.0_wp, restore_lambda=20.0_wp, &
                        sw_down=0.0_wp)

         hT_after = ms%tracers(ms%idx_temperature)%hTr(ip, jp, NZ)
         call check(error, abs(ice%fb(ip, jp)) < 1.0e-12_wp, &
                    "fb must be ~0 when SST == T_f exactly")
         if (allocated(error)) exit checks
         call check(error, abs(sf%Q_heat(ip, jp)) < 1.0e-9_wp, &
                    "Q_heat must be ~0 on a pure cold step with no pre-existing ice")
         if (allocated(error)) exit checks
         call check(error, abs(hT_after - hT_before) < 1.0e-9_wp, &
                    "ocean surface heat content must be unchanged on the grow step")
         call teardown(ms, eos, ice, sf)
         if (allocated(error)) exit checks

         ! ---- Melt sub-block: seeded ice, warm restoring => heat_to_ocn
         ! and/or fb nonzero. Verify the closure identity through the real
         ! apply. ----
         call setup_state(grid, ms, eos, ice, sf, 2.0_wp)
         ip = grid%nx_total/2
         jp = grid%ny_total/2
         call seed_ice_cap(ice, ip, jp, H_ICE0, T_ICE_COLD)
         hT_before = ms%tracers(ms%idx_temperature)%hTr(ip, jp, NZ)

         call run_chain(grid, eos, ms, ice, sf, air_temp=15.0_wp, restore_lambda=25.0_wp, &
                        sw_down=0.0_wp)

         hT_after = ms%tracers(ms%idx_temperature)%hTr(ip, jp, NZ)
         d_ocean_heat = sf%rho0*SEAWATER_CP*(hT_after - hT_before)
         rhs = ice%heat_to_ocn(ip, jp, 1) - ice%fb(ip, jp)*DT_THERM
         scale = sf%rho0*SEAWATER_CP*abs(hT_before)
         tol = 1.0e-9_wp*max(scale, 1.0_wp)
         call check(error, abs(d_ocean_heat - rhs) <= tol, &
                    "rho0*cp*Delta(hT_top) must equal heat_to_ocn - fb*dt_therm to round-off")

      end block checks
      call teardown(ms, eos, ice, sf)
   end subroutine test_heat_budget_closes

   subroutine test_pure_melt_conserves(error)
      !! Case 4 — the column<->coupling-seam heat/mass/salt bookkeeping,
      !! isolated from the frazil bank (SPEC §8 case 4, "simplest
      !! closable formulation"). A single-layer (`nk_ice=1`) ice cap is
      !! seeded ISOTHERMAL AT `tfw` (the seawater freezing point at
      !! S_INIT), so the column's internal conduction exchange with the
      !! `tfw` boundary condition is exactly zero (`tflux_bot=0` — no
      !! temperature gradient to conduct against) and with `sf_0=dsf_dt=
      !! sw_down=0` there is no SEB/solar input either; a warm `fb`
      !! (via warm SST) is then the ONLY energy source touching the
      !! column. Under exactly this configuration the identity
      !!   heat_to_ocn == fb*dt_therm - h2o_ice_to_ocn*(enth_fr - enth_ice_after)
      !! holds to round-off, where `enth_fr = ice_enthalpy_liquid_freeze
      !! (sal_ice_after)` is the freeze-point enthalpy the melt-peel
      !! kernel (`ice_bottom_melt_peel`, `rdb_ice_mass.F90`) prices the
      !! melted mass at, and `enth_ice_after` is the surviving layer's
      !! post-conduction specific enthalpy (unchanged by the peel itself
      !! — `rdb_ice_mass%ice_bottom_melt_peel` only removes MASS from a
      !! layer, never touches the survivor's specific enthalpy). This was
      !! verified numerically (a standalone probe against `ice_temp_sis2`
      !! + `ice_column_step` directly, residual ~1e-10 against an
      !! O(1e6) J/m^2 scale, across fb in {500,2000,5000} W/m^2) — WITHOUT
      !! the tfw-seeding, the column's OWN conduction exchange with `tfw`
      !! (`tflux_bot`, internal and not exposed by `ice_thermo_columns`'s
      !! public interface) adds an independent ~0.2-0.3% term that this
      !! identity cannot see from driver-visible quantities alone.
      !!
      !! Salt: the melt returns exactly the ice's own bulk salinity (no
      !! brine rejection on a pure melt) — `mass_returned*ICE_BULK_
      !! SALINITY == -d(ice salt content)`.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      type(ocean_surface_flux_t) :: sf
      real(wp) :: t_f, enth_at_tfw
      real(wp) :: m_ice_before, salt_ice_before, m_ice_after, salt_ice_after
      real(wp) :: d_m_ice, d_salt_ice, mass_returned, salt_returned
      real(wp) :: enth_fr, rhs, tol
      integer :: ip, jp
      checks: block

         call grid%init(6, 4, NGHOST, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call eos%init(grid)
         ice%enable = .true.
         ice%ncat = 1
         ice%nk_ice = 1
         call ice%init(grid)
         call sf%init(grid)
         call sf%set_surface_flux_const(0.0_wp, 0.0_wp)

         ms%h_layer = H_LAYER
         ms%tracers(ms%idx_salinity)%hTr = S_INIT*H_LAYER
         ms%tracers(ms%idx_temperature)%hTr = T_DEEP*H_LAYER
         t_f = eos_freezing_point(eos, S_INIT, 0.0_wp)
         ms%tracers(ms%idx_temperature)%hTr(:, :, NZ) = 5.0_wp*H_LAYER

         ip = grid%nx_total/2
         jp = grid%ny_total/2
         enth_at_tfw = ice_enth_from_ts(t_f, ICE_BULK_SALINITY)
         ice%m_ice(ip, jp, 1) = H_ICE0*ICE_RHO_ICE
         ice%m_snow(ip, jp, 1) = 0.0_wp
         ice%enth_snow(ip, jp, 1, 1) = ice_enth_from_ts(0.0_wp, 0.0_wp)
         ice%enth_ice(ip, jp, 1, 1) = enth_at_tfw
         ice%sal_ice(ip, jp, 1, 1) = ICE_BULK_SALINITY

         m_ice_before = ice%m_ice(ip, jp, 1)
         salt_ice_before = m_ice_before*ice%sal_ice(ip, jp, 1, 1)

         ! Passive column (no SEB/solar): air_temp/restore_lambda/sw_down
         ! all inert (restore_lambda=0 ⇒ dsf_dt=0, sf_0=0 regardless of
         ! air_temp). Warm SST (5 degC vs t_f) drives fb > 0 via
         ! ice_compute_basal_flux — the ONLY energy source touching the
         ! column in this configuration.
         call run_chain(grid, eos, ms, ice, sf, air_temp=0.0_wp, restore_lambda=0.0_wp, &
                        sw_down=0.0_wp)

         m_ice_after = ice%m_ice(ip, jp, 1)
         salt_ice_after = m_ice_after*ice%sal_ice(ip, jp, 1, 1)
         d_m_ice = m_ice_after - m_ice_before
         d_salt_ice = salt_ice_after - salt_ice_before

         call check(error, ice%fb(ip, jp) > 0.0_wp, &
                    "fb must be strictly positive (warm SST drives basal melt)")
         if (allocated(error)) exit checks
         call check(error, d_m_ice < 0.0_wp, "the cap must thin (pure melt, no freeze)")
         if (allocated(error)) exit checks
         call check(error, ice%h2o_ocn_to_ice(ip, jp, 1) == 0.0_wp, &
                    "no freeze must occur (pure-melt configuration)")
         if (allocated(error)) exit checks

         ! ---- Mass: what left the ice must equal the ice mass loss. ----
         mass_returned = ice%h2o_ice_to_ocn(ip, jp, 1)
         tol = 1.0e-9_wp*max(abs(d_m_ice), 1.0e-6_wp)
         call check(error, abs(mass_returned - (-d_m_ice)) <= tol, &
                    "meltwater mass reported must equal the ice mass LOSS")
         if (allocated(error)) exit checks

         ! ---- Heat: the tfw-seeded closed identity (module docstring). ----
         enth_fr = ice_enthalpy_liquid_freeze(ice%sal_ice(ip, jp, 1, 1))
         rhs = ice%fb(ip, jp)*DT_THERM &
               - ice%h2o_ice_to_ocn(ip, jp, 1)*(enth_fr - ice%enth_ice(ip, jp, 1, 1))
         tol = 1.0e-9_wp*max(abs(rhs), 1.0_wp)
         call check(error, abs(ice%heat_to_ocn(ip, jp, 1) - rhs) <= tol, &
                    "heat_to_ocn must equal fb*dt_therm - melt*(enth_fr-enth_after) to round-off")
         if (allocated(error)) exit checks

         ! ---- Salt: the melt returns exactly the ice's own bulk salinity. ----
         salt_returned = mass_returned*ICE_BULK_SALINITY
         tol = 1.0e-9_wp*max(abs(salt_returned), abs(d_salt_ice), 1.0e-6_wp)
         call check(error, abs((-d_salt_ice) - salt_returned) <= tol, &
                    "ice salt content lost must equal net-melt-mass*ICE_BULK_SALINITY")

      end block checks
      call teardown(ms, eos, ice, sf)
   end subroutine test_pure_melt_conserves

   subroutine test_icefree_ocean_untouched(error)
      !! Case 5 — the OPEN-OCEAN regression gate for the ice-presence gate
      !! on `fb`. Ice ENABLED, a WARM ice-free surface (SST = 5 degC, well
      !! above T_f), NO ice seeded ANYWHERE, NO frazil bank. On such a cell
      !! `ice_compute_basal_flux` would (without the gate) compute a large
      !! `fb = RHO_WATER*SEAWATER_CP*(SST-T_f)*h/dt > 0`, but the column
      !! never consumes it (m_ice=0 ⇒ heat_to_ocn=0), so the melt-side
      !! reduce kernel would inject a spurious `Q_heat = -fb`, cooling the
      !! open ocean by several degC PER thermo step. The ice-presence gate
      !! (`rdb_ice_basal_flux`: fb=0 unless `sum_cat m_ice >
      !! ICE_RHO_ICE*H_VANISHED`) makes fb=0 here, so the whole ice path is
      !! a clean no-op on open water. This test asserts exactly that: fb=0
      !! at an interior cell, Q_heat/Q_salt all-zero, and hT/hS BIT-IDENTICAL
      !! before/after the full chain through the REAL apply. (Fails on the
      !! pre-gate code — verified by the implementer — passes after.)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      type(ocean_surface_flux_t) :: sf
      real(wp), allocatable :: hT_ic(:, :, :), hS_ic(:, :, :)
      integer :: ip, jp
      checks: block

         ! Warm ice-free surface (SST = 5 degC, well above T_f), no ice.
         call setup_state(grid, ms, eos, ice, sf, 5.0_wp)
         ip = grid%nx_total/2
         jp = grid%ny_total/2
         allocate (hT_ic, source=ms%tracers(ms%idx_temperature)%hTr)
         allocate (hS_ic, source=ms%tracers(ms%idx_salinity)%hTr)

         ! No ice seeded, no frazil bank. Restoring params are inert on an
         ! ice-free column (the column never steps without ice). PR 31:
         ! sw_down = 400 W/m^2 so the shortwave trap is VISIBLE — the ice
         ! column gates `sw_thru = 0` on every ice-free cell
         ! (`rdb_ice_column`), so the open ocean must see NO ice-driven
         ! shortwave. `atm_sw_dn` is written to every cell (incl ghosts) by
         ! the forcing kernel; reading it instead of the column-gated
         ! `sw_thru`/`sw_thru_diag` would spray the whole basin — this is
         ! the guard.
         call run_chain(grid, eos, ms, ice, sf, air_temp=5.0_wp, restore_lambda=20.0_wp, &
                        sw_down=400.0_wp)

         call check(error, ice%fb(ip, jp) == 0.0_wp, &
                    "fb must be exactly zero on an ice-free cell (ice-presence gate)")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(ice%sw_thru_diag)) == 0.0_wp, &
                    "sw_thru_diag must be all-zero on ice-free water (column SW gate)")
         if (allocated(error)) exit checks
         call check(error, ice%sw_thru_diag(ip, jp) == 0.0_wp, &
                    "sw_thru_diag must be exactly zero at the interior ice-free cell")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(ice%heat_flux_diag)) == 0.0_wp, &
                    "heat_flux_diag must be all-zero (no ice, fb gated to 0)")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(sf%Q_heat)) == 0.0_wp, &
                    "Q_heat must be all-zero on the open ocean (no -fb cooling, no SW spray)")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(sf%Q_salt)) == 0.0_wp, &
                    "Q_salt must be all-zero (no ice exchange)")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(ms%tracers(ms%idx_temperature)%hTr - hT_ic)) == 0.0_wp, &
                    "hT must be bit-identical (open ocean untouched by the ice path)")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(ms%tracers(ms%idx_salinity)%hTr - hS_ic)) == 0.0_wp, &
                    "hS must be bit-identical (open ocean untouched by the ice path)")

      end block checks
      if (allocated(hT_ic)) deallocate (hT_ic)
      if (allocated(hS_ic)) deallocate (hS_ic)
      call teardown(ms, eos, ice, sf)
   end subroutine test_icefree_ocean_untouched

   subroutine test_disabled_bitident(error)
      !! Disabled-ice path contract: with ice OFF the driver's
      !! thermo-cadence ice block (`if (ocean_state%ice%enable)`) never
      !! fires, so a benign window is an exact no-op on the ocean tracers.
      !! Carried by ONE real gate — Case 5b: ice ON, benign restoring
      !! (air_temp = T_f, sw_down = 0), NO pre-existing ice, NO frazil
      !! bank — one window must leave Q_heat = Q_salt = 0 and the ocean
      !! tracers bit-identical. (The former Case 5a asserted
      !! init-guaranteed values — `has_heat`/`has_salt` false, Q_heat/
      !! Q_salt zero — WITHOUT invoking any driver path, a vacuous proxy;
      !! it was removed, keeping only the freezing-point `t_f` seed it
      !! computed for 5b.)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      type(ocean_surface_flux_t) :: sf
      real(wp), allocatable :: hT_ic(:, :, :), hS_ic(:, :, :)
      real(wp) :: t_f
      checks: block

         ! ---- (Former 5a REMOVED.) 5a asserted init-guaranteed values
         ! (has_heat/has_salt false, Q_heat/Q_salt zero) WITHOUT invoking
         ! any driver path — a self-documented proxy, not a real
         ! disabled-path gate. The disabled-path contract is carried
         ! solely by 5b's real `run_chain` bit-identity check below. Here
         ! we only compute the freezing point `t_f` that seeds 5b's
         ! surface layer (a live eos handle is all `eos_freezing_point`
         ! needs; setup_state re-inits grid/eos/ms/sf for 5b). ----
         call grid%init(6, 4, NGHOST, 1.0_wp, 1.0_wp)
         call eos%init(grid)
         t_f = eos_freezing_point(eos, S_INIT, 0.0_wp)
         call eos%destroy()

         ! ---- 5b: ice ON, benign restoring, no ice, no bank => exact
         ! no-op window. ----
         call setup_state(grid, ms, eos, ice, sf, t_f)
         allocate (hT_ic, source=ms%tracers(ms%idx_temperature)%hTr)
         allocate (hS_ic, source=ms%tracers(ms%idx_salinity)%hTr)

         call run_chain(grid, eos, ms, ice, sf, air_temp=t_f, restore_lambda=20.0_wp, &
                        sw_down=0.0_wp)

         call check(error, maxval(abs(sf%Q_heat)) == 0.0_wp, &
                    "Q_heat must be exactly zero (no ice, fb=0, heat_to_ocn=0)")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(sf%Q_salt)) == 0.0_wp, &
                    "Q_salt must be exactly zero (no ice, no frazil bank)")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(ms%tracers(ms%idx_temperature)%hTr - hT_ic)) == 0.0_wp, &
                    "hT must be bit-identical on the benign no-op window")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(ms%tracers(ms%idx_salinity)%hTr - hS_ic)) == 0.0_wp, &
                    "hS must be bit-identical on the benign no-op window")

      end block checks
      if (allocated(hT_ic)) deallocate (hT_ic)
      if (allocated(hS_ic)) deallocate (hS_ic)
      call teardown(ms, eos, ice, sf)
   end subroutine test_disabled_bitident

   ! -----------------------------------------------------------------
   ! Compose cases: frazil + melt in one window (PR-60)
   ! -----------------------------------------------------------------
   !
   ! Shared compose setup for the four cases below: surface at t_f (the
   ! POST-CLAMP freezing point) makes fb == 0 exactly
   ! (`ice_compute_basal_flux_impl`: `max(0, sst-tfw) = 0`) and is the
   ! only state a live frazil bank is physically coherent in (the PR-1
   ! clamp raises SST to T_f when it banks the deficit). A cold cap
   ! (H_ICE0, T_ICE_COLD) plus a warm restoring pair (air_temp=10,
   ! restore_lambda=30 — the exact pair `melt_freshens_ocean` already
   ! proves melts this cap) then makes the surface melt while the base
   ! nucleates frazil in the SAME window: an autumn marginal-ice-zone
   ! cell, supercooled water making frazil under a pack a warm airmass is
   ! ablating from above.

   subroutine test_compose_frazil_and_melt_salt(error)
      !! The compose headline (RESUME item a). On the compose cell, BOTH
      !! contributors to `salt_flux_diag` are live this window:
      !! `m_frozen_diag > 0` (frazil, `rdb_ice_frazil_uptake`) and
      !! `h2o_ice_to_ocn > 0` (melt, the column driver). Then, to 1e-12
      !! rel: `salt_flux_diag*DT_THERM == (m_frozen_diag+m_net)*
      !! (S_INIT-ICE_BULK_SALINITY)`, with `m_net = h2o_ocn_to_ice -
      !! h2o_ice_to_ocn`. Also `ssurf_seam == S_INIT` exactly (the frazil
      !! kernel's fresh SSS sample and the reduce kernel's `ssurf_seam`
      !! agree — both read the same unmutated ocean, §3.1 of the PR-60
      !! plan) and `Q_salt == salt_flux_diag` exactly (zero background).
      !!
      !! **This is the negative-control gate**: with
      !! `rdb_ice_thermo_driver.F90:203` temporarily changed from
      !! `salt_flux_diag = salt_flux_diag + ...` (compose) to
      !! `salt_flux_diag = ...` (overwrite), this case FAILS — the
      !! `m_frozen_diag` contribution vanishes from the identity. No
      !! other test in the suite catches that mutation (verified: every
      !! other test that calls the reduce kernel has an empty bank, so
      !! `0 + x == x` hides the bug).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      type(ocean_surface_flux_t) :: sf
      real(wp) :: t_f, m_net, lhs, rhs, tol
      integer :: ip, jp
      checks: block

         t_f = eos_freezing_point(eos, S_INIT, 0.0_wp)
         call setup_state(grid, ms, eos, ice, sf, t_f)
         ip = grid%nx_total/2
         jp = grid%ny_total/2
         call seed_ice_cap(ice, ip, jp, H_ICE0, T_ICE_COLD)
         ice%frazil_heat(ip, jp) = E_BANK_COMPOSE

         call run_chain(grid, eos, ms, ice, sf, air_temp=10.0_wp, restore_lambda=30.0_wp, &
                        sw_down=0.0_wp)

         call check(error, ice%m_frozen_diag(ip, jp) > 0.0_wp, &
                    "m_frozen_diag must be strictly positive (frazil term live)")
         if (allocated(error)) exit checks
         call check(error, ice%h2o_ice_to_ocn(ip, jp, 1) > 0.0_wp, &
                    "h2o_ice_to_ocn must be strictly positive (melt term live)")
         if (allocated(error)) exit checks

         m_net = ice%h2o_ocn_to_ice(ip, jp, 1) - ice%h2o_ice_to_ocn(ip, jp, 1)
         lhs = ice%salt_flux_diag(ip, jp)*DT_THERM
         rhs = (ice%m_frozen_diag(ip, jp) + m_net)*(S_INIT - ICE_BULK_SALINITY)
         tol = 1.0e-12_wp*max(abs(rhs), 1.0e-6_wp)
         call check(error, abs(lhs - rhs) <= tol, &
                    "salt_flux_diag*dt must equal (m_frozen+m_net)*(s_surf-bulk) to 1e-12 rel")
         if (allocated(error)) exit checks

         call check(error, ice%ssurf_seam(ip, jp) == S_INIT, &
                    "ssurf_seam must equal S_INIT exactly (both kernels sample the same &
                    &unmutated ocean)")
         if (allocated(error)) exit checks
         call check(error, sf%Q_salt(ip, jp) == ice%salt_flux_diag(ip, jp), &
                    "Q_salt must equal salt_flux_diag exactly (zero background const)")

      end block checks
      call teardown(ms, eos, ice, sf)
   end subroutine test_compose_frazil_and_melt_salt

   subroutine test_compose_frazil_ice_melts_same_window(error)
      !! The order gate, in one assertion. An ice-FREE wet cell
      !! (`m_ice == 0` everywhere at window start), surface at t_f, a
      !! live bank, warm restoring. After one `run_chain`:
      !! `m_frozen_diag > 0` (frazil made ice) AND `h2o_ice_to_ocn > 0`
      !! (that SAME ice melted, this window). Also `fb == 0.0` exactly
      !! (`ice_compute_basal_flux` samples `m_ice` BEFORE the frazil
      !! uptake runs — the mandated order, `rdb_driver.F90` — so
      !! frazil-nucleated ice feels `fb = 0` its first window; that is
      !! correct and not a bug, see the PR-60 plan §11.3).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      type(ocean_surface_flux_t) :: sf
      real(wp) :: t_f
      integer :: ip, jp
      checks: block

         t_f = eos_freezing_point(eos, S_INIT, 0.0_wp)
         call setup_state(grid, ms, eos, ice, sf, t_f)
         ip = grid%nx_total/2
         jp = grid%ny_total/2
         ! No seed_ice_cap: m_ice(ip,jp,1) == 0 at window start.
         ice%frazil_heat(ip, jp) = E_BANK_COMPOSE

         call run_chain(grid, eos, ms, ice, sf, air_temp=10.0_wp, restore_lambda=30.0_wp, &
                        sw_down=0.0_wp)

         call check(error, ice%fb(ip, jp) == 0.0_wp, &
                    "fb must be exactly zero (basal flux sampled m_ice before frazil ran)")
         if (allocated(error)) exit checks
         call check(error, ice%m_frozen_diag(ip, jp) > 0.0_wp, &
                    "frazil must have nucleated new ice this window")
         if (allocated(error)) exit checks
         call check(error, ice%h2o_ice_to_ocn(ip, jp, 1) > 0.0_wp, &
                    "the frazil-nucleated ice must melt in the SAME window (order gate)")

      end block checks
      call teardown(ms, eos, ice, sf)
   end subroutine test_compose_frazil_ice_melts_same_window

   subroutine test_compose_heat_not_double_counted(error)
      !! With a LIVE bank and live melt in the same window, the frazil
      !! latent heat contributes exactly nothing to `Q_heat`
      !! (`rdb_ice_frazil_uptake` never writes `heat_flux_diag` — writing
      !! one would double-count the PR-1 surface-clamp credit). Asserts
      !! `heat_flux_diag*DT_THERM == Sum_c(heat_to_ocn) - fb*DT_THERM`
      !! to 1e-12 rel (ncat=1: a single term), AND, through the REAL
      !! `ocean_surface_flux_apply_tracers` kernel,
      !! `rho0*SEAWATER_CP*Delta(hT_top) == heat_flux_diag*DT_THERM` to
      !! 1e-9 rel — the same identity `heat_budget_closes` checks with a
      !! DEAD bank; this is the live-bank twin.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      type(ocean_surface_flux_t) :: sf
      real(wp) :: t_f, hT_before, hT_after, d_ocean_heat, rhs, tol, scale
      integer :: ip, jp
      checks: block

         t_f = eos_freezing_point(eos, S_INIT, 0.0_wp)
         call setup_state(grid, ms, eos, ice, sf, t_f)
         ip = grid%nx_total/2
         jp = grid%ny_total/2
         call seed_ice_cap(ice, ip, jp, H_ICE0, T_ICE_COLD)
         ice%frazil_heat(ip, jp) = E_BANK_COMPOSE
         hT_before = ms%tracers(ms%idx_temperature)%hTr(ip, jp, NZ)

         call run_chain(grid, eos, ms, ice, sf, air_temp=10.0_wp, restore_lambda=30.0_wp, &
                        sw_down=0.0_wp)

         call check(error, ice%m_frozen_diag(ip, jp) > 0.0_wp, &
                    "frazil term must be live (bank spent)")
         if (allocated(error)) exit checks
         call check(error, ice%h2o_ice_to_ocn(ip, jp, 1) > 0.0_wp, &
                    "melt term must be live")
         if (allocated(error)) exit checks

         rhs = ice%heat_to_ocn(ip, jp, 1) - ice%fb(ip, jp)*DT_THERM
         tol = 1.0e-12_wp*max(abs(rhs), 1.0e-6_wp)
         call check(error, abs(ice%heat_flux_diag(ip, jp)*DT_THERM - rhs) <= tol, &
                    "heat_flux_diag*dt must equal heat_to_ocn - fb*dt (frazil contributes 0)")
         if (allocated(error)) exit checks

         hT_after = ms%tracers(ms%idx_temperature)%hTr(ip, jp, NZ)
         d_ocean_heat = sf%rho0*SEAWATER_CP*(hT_after - hT_before)
         scale = sf%rho0*SEAWATER_CP*abs(hT_before)
         tol = 1.0e-9_wp*max(scale, 1.0_wp)
         call check(error, abs(d_ocean_heat - ice%heat_flux_diag(ip, jp)*DT_THERM) <= tol, &
                    "rho0*cp*Delta(hT_top) must equal heat_flux_diag*dt to round-off")

      end block checks
      call teardown(ms, eos, ice, sf)
   end subroutine test_compose_heat_not_double_counted

   subroutine test_compose_mass_closes(error)
      !! The two-term mass identity: `Delta(m_ice+m_snow)(ip,jp,1) ==
      !! m_frozen_diag - m_melt_diag`, to 1e-12 rel, with BOTH the
      !! frazil term (`m_frozen_diag > 0`) and the melt term
      !! (`h2o_ice_to_ocn > 0`) strictly live — the `ncat=1` twin of
      !! `test_multicat_thermo_conserves`'s mass identity, which today
      !! runs with an empty bank and degenerates to a one-term
      !! `Delta(m) == -m_melt_diag`. Additionally `frazil_heat == 0`
      !! exactly (the bank is fully spent every window it is nonzero).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      type(ocean_surface_flux_t) :: sf
      real(wp) :: t_f, m_before, m_after, d_m, expected_d_m, tol
      integer :: ip, jp
      checks: block

         t_f = eos_freezing_point(eos, S_INIT, 0.0_wp)
         call setup_state(grid, ms, eos, ice, sf, t_f)
         ip = grid%nx_total/2
         jp = grid%ny_total/2
         call seed_ice_cap(ice, ip, jp, H_ICE0, T_ICE_COLD)
         ice%frazil_heat(ip, jp) = E_BANK_COMPOSE
         m_before = ice%m_ice(ip, jp, 1) + ice%m_snow(ip, jp, 1)

         call run_chain(grid, eos, ms, ice, sf, air_temp=10.0_wp, restore_lambda=30.0_wp, &
                        sw_down=0.0_wp)

         call check(error, ice%m_frozen_diag(ip, jp) > 0.0_wp, &
                    "frazil term must be strictly live")
         if (allocated(error)) exit checks
         call check(error, ice%h2o_ice_to_ocn(ip, jp, 1) > 0.0_wp, &
                    "melt term must be strictly live")
         if (allocated(error)) exit checks
         call check(error, ice%frazil_heat(ip, jp) == 0.0_wp, &
                    "the bank must be fully spent (frazil_heat reset to 0)")
         if (allocated(error)) exit checks

         m_after = ice%m_ice(ip, jp, 1) + ice%m_snow(ip, jp, 1)
         d_m = m_after - m_before
         expected_d_m = ice%m_frozen_diag(ip, jp) - ice%m_melt_diag(ip, jp)
         tol = 1.0e-12_wp*max(abs(expected_d_m), 1.0e-6_wp)
         call check(error, abs(d_m - expected_d_m) <= tol, &
                    "Delta(m_ice+m_snow) must equal m_frozen_diag - m_melt_diag to 1e-12 rel")

      end block checks
      call teardown(ms, eos, ice, sf)
   end subroutine test_compose_mass_closes

   ! -----------------------------------------------------------------
   ! PR 31: ice -> ocean shortwave (sw_thru) coupling
   ! -----------------------------------------------------------------

   subroutine test_sw_thru_reaches_ocean(error)
      !! PR 31 headline (the arbiter). A thin snow-free cap (h_i = 0.2 m)
      !! over a surface EXACTLY at t_f (so fb = 0), passive SEB
      !! (restore_lambda = 0), sw_down = 400 W/m^2. Through the REAL chain
      !! (couplers + apply_tracers), components OFF:
      !!   (a) sw_thru(ip,jp,1) > 0             — thin ice transmits SW;
      !!   (b) hT_after > hT_before             — thin ice WARMS the water;
      !!   (c) the closure identity, restored:
      !!       |rho0*cp*Delta(hT_top) - (heat_to_ocn - fb*dt + sw_thru*dt)|
      !!         <= 1e-9*scale.
      !! This is `heat_budget_closes`'s melt identity with the one missing
      !! shortwave term put back. It FAILS by -sw_thru*dt on a miss and by
      !! +sw_thru*dt on a double-count, so it pins the branch-2 accounting.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      type(ocean_surface_flux_t) :: sf
      real(wp) :: t_f, hT_before, hT_after, d_ocean_heat, rhs, tol, scale
      integer :: ip, jp
      checks: block

         t_f = eos_freezing_point(eos, S_INIT, 0.0_wp)
         call setup_state(grid, ms, eos, ice, sf, t_f)
         ip = grid%nx_total/2
         jp = grid%ny_total/2
         call seed_ice_cap(ice, ip, jp, 0.2_wp, t_f)
         hT_before = ms%tracers(ms%idx_temperature)%hTr(ip, jp, NZ)

         call run_chain(grid, eos, ms, ice, sf, air_temp=0.0_wp, restore_lambda=0.0_wp, &
                        sw_down=400.0_wp)

         hT_after = ms%tracers(ms%idx_temperature)%hTr(ip, jp, NZ)

         call check(error, ice%sw_thru(ip, jp, 1) > 0.0_wp, &
                    "sw_thru must be strictly positive (thin snow-free ice transmits SW)")
         if (allocated(error)) exit checks
         call check(error, ice%sw_thru_diag(ip, jp) == ice%sw_thru(ip, jp, 1), &
                    "sw_thru_diag must equal sw_thru (ncat=1 lumped reduction)")
         if (allocated(error)) exit checks
         call check(error, hT_after > hT_before, &
                    "thin ice must WARM the water beneath it (sw_thru reaches the ocean)")
         if (allocated(error)) exit checks

         d_ocean_heat = sf%rho0*SEAWATER_CP*(hT_after - hT_before)
         rhs = ice%heat_to_ocn(ip, jp, 1) - ice%fb(ip, jp)*DT_THERM &
               + ice%sw_thru(ip, jp, 1)*DT_THERM
         scale = sf%rho0*SEAWATER_CP*abs(hT_before)
         tol = 1.0e-9_wp*max(scale, 1.0_wp)
         call check(error, abs(d_ocean_heat - rhs) <= tol, &
                    "rho0*cp*Delta(hT_top) must equal heat_to_ocn - fb*dt + sw_thru*dt")

      end block checks
      call teardown(ms, eos, ice, sf)
   end subroutine test_sw_thru_reaches_ocean

   subroutine test_sw_thru_beers_law(error)
      !! The transmission law itself. Two snow-free columns at h_i = 1.0 m
      !! and h_i = 2.0 m, same sw_down, same seeded (deep-cold) ice
      !! temperature. Both are above 0.5 m where the CSIM4 thin-ice albedo
      !! ramp saturates (fh = 1), so the albedo is identical and factors
      !! out of the ratio exactly. The remaining thickness dependence is
      !! pure Beer-Lambert: sw_thru ~ exp(-h/ICE_OPT_DEP_ICE).
      !!   |sw_thru_diag(2)/sw_thru_diag(1) - exp(-1/ICE_OPT_DEP_ICE)| <= 1e-10.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      type(ocean_surface_flux_t) :: sf
      real(wp) :: ratio, expected
      integer :: ip, jp
      checks: block

         call setup_state(grid, ms, eos, ice, sf, 0.0_wp)
         ip = grid%nx_total/2
         jp = grid%ny_total/2
         ! Two independent interior columns, different thicknesses.
         call seed_ice_cap(ice, ip, jp, 1.0_wp, T_ICE_COLD)
         call seed_ice_cap(ice, ip + 1, jp, 2.0_wp, T_ICE_COLD)

         call run_chain(grid, eos, ms, ice, sf, air_temp=0.0_wp, restore_lambda=0.0_wp, &
                        sw_down=400.0_wp)

         call check(error, ice%sw_thru_diag(ip, jp) > 0.0_wp, &
                    "sw_thru_diag(1m) must be strictly positive")
         if (allocated(error)) exit checks
         call check(error, ice%sw_thru_diag(ip + 1, jp) > 0.0_wp, &
                    "sw_thru_diag(2m) must be strictly positive")
         if (allocated(error)) exit checks

         ratio = ice%sw_thru_diag(ip + 1, jp)/ice%sw_thru_diag(ip, jp)
         expected = exp(-1.0_wp/ICE_OPT_DEP_ICE)
         call check(error, abs(ratio - expected) <= 1.0e-10_wp, &
                    "sw_thru_diag(2m)/sw_thru_diag(1m) must equal exp(-1/opt_dep) (Beer-Lambert)")

      end block checks
      call teardown(ms, eos, ice, sf)
   end subroutine test_sw_thru_beers_law

   subroutine test_thick_ice_nearly_opaque(error)
      !! The correct form of "opaque thick ice": abs_ocn = pen*exp(-h/opt)
      !! is > 0 for every finite h, so thick ice is NEARLY opaque, never
      !! exactly. h_i = 3.0 m, snow-free, sw_down = 400. The bound uses
      !! (1 - albedo) <= 1 and pen <= ICE_PEN_ICE:
      !!   0 < sw_thru_diag <= ICE_PEN_ICE*exp(-3/ICE_OPT_DEP_ICE)*sw_down.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      type(ocean_surface_flux_t) :: sf
      real(wp) :: bound
      integer :: ip, jp
      checks: block

         call setup_state(grid, ms, eos, ice, sf, 0.0_wp)
         ip = grid%nx_total/2
         jp = grid%ny_total/2
         call seed_ice_cap(ice, ip, jp, 3.0_wp, T_ICE_COLD)

         call run_chain(grid, eos, ms, ice, sf, air_temp=0.0_wp, restore_lambda=0.0_wp, &
                        sw_down=400.0_wp)

         bound = ICE_PEN_ICE*exp(-3.0_wp/ICE_OPT_DEP_ICE)*400.0_wp
         call check(error, ice%sw_thru_diag(ip, jp) > 0.0_wp, &
                    "sw_thru_diag must be strictly positive (never exactly opaque)")
         if (allocated(error)) exit checks
         call check(error, ice%sw_thru_diag(ip, jp) <= bound, &
                    "sw_thru_diag must be <= ICE_PEN_ICE*exp(-3/opt_dep)*sw_down (nearly opaque)")

      end block checks
      call teardown(ms, eos, ice, sf)
   end subroutine test_thick_ice_nearly_opaque

   subroutine test_sw_default_bitident(error)
      !! The default-off bit-identity gate. sw_down = 0.0 (the shipped
      !! default) with a seeded cap under warm restoring (melt-driving).
      !! sw_down = 0 => sw_tot = 0 => sw_thru = 0 exactly, so:
      !!   * sw_thru_diag == 0 on EVERY cell (exact ==);
      !!   * sf%Q_heat(ip,jp) == heat_to_ocn/dt - fb (the PRE-PR-31
      !!     formula, exact ==) — components off, so the SW coupler's
      !!     `Q_heat += sw_thru_diag` is an exact +0.0 and the value is
      !!     byte-for-byte what the pre-PR-31 code produced.
      !! Exact == (not tolerance) is deliberate: a tolerance would hide a
      !! `+0.0*x` reassociation.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      type(ocean_surface_flux_t) :: sf
      real(wp) :: q_heat_pre
      integer :: ip, jp
      checks: block

         call setup_state(grid, ms, eos, ice, sf, 2.0_wp)
         ip = grid%nx_total/2
         jp = grid%ny_total/2
         call seed_ice_cap(ice, ip, jp, H_ICE0, T_ICE_COLD)

         call run_chain(grid, eos, ms, ice, sf, air_temp=15.0_wp, restore_lambda=25.0_wp, &
                        sw_down=0.0_wp)

         call check(error, maxval(abs(ice%sw_thru_diag)) == 0.0_wp, &
                    "sw_thru_diag must be exactly zero everywhere at sw_down = 0")
         if (allocated(error)) exit checks
         q_heat_pre = ice%heat_to_ocn(ip, jp, 1)/DT_THERM - ice%fb(ip, jp)
         call check(error, sf%Q_heat(ip, jp) == q_heat_pre, &
                    "Q_heat must equal the pre-PR-31 formula exactly (sw_thru_diag = +0.0)")

      end block checks
      call teardown(ms, eos, ice, sf)
   end subroutine test_sw_default_bitident

   subroutine test_sw_components_q_sw(error)
      !! The components-ON no-double-count gate. Same thin-cap / fb = 0 /
      !! passive-SEB / sw_down = 400 setup as `sw_thru_reaches_ocean`, but
      !! with the PR-12 component surface-flux path ON, so the shortwave
      !! travels the OTHER branch: the coupler fills the `q_sw` component
      !! and `ocean_surface_flux_assemble` sums it into `Q_heat`. Asserts:
      !!   * q_sw(ip,jp) == sw_thru_diag(ip,jp)  exactly (the coupler fill);
      !!   * Q_heat(ip,jp) == heat_flux_diag + sw_thru_diag exactly — the
      !!     shortwave is in Q_heat EXACTLY ONCE (heat_flux_diag went to
      !!     `heat_added`, q_sw is the SW summand; a double-count would
      !!     make this fail by +sw_thru_diag);
      !!   * the same closure identity as the components-off arbiter, so
      !!     the energy lands in the ocean exactly once in THIS mode too.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      type(ocean_surface_flux_t) :: sf
      real(wp) :: t_f, hT_before, hT_after, d_ocean_heat, rhs, tol, scale
      integer :: ip, jp
      checks: block

         t_f = eos_freezing_point(eos, S_INIT, 0.0_wp)
         call setup_state(grid, ms, eos, ice, sf, t_f)
         ip = grid%nx_total/2
         jp = grid%ny_total/2
         call seed_ice_cap(ice, ip, jp, 0.2_wp, t_f)
         hT_before = ms%tracers(ms%idx_temperature)%hTr(ip, jp, NZ)

         call run_chain(grid, eos, ms, ice, sf, air_temp=0.0_wp, restore_lambda=0.0_wp, &
                        sw_down=400.0_wp, use_components=.true.)

         hT_after = ms%tracers(ms%idx_temperature)%hTr(ip, jp, NZ)

         call check(error, ice%sw_thru_diag(ip, jp) > 0.0_wp, &
                    "sw_thru_diag must be strictly positive (else the test is vacuous)")
         if (allocated(error)) exit checks
         call check(error, sf%has_q_sw, "has_q_sw must be latched true (components on)")
         if (allocated(error)) exit checks
         call check(error, sf%q_sw(ip, jp) == ice%sw_thru_diag(ip, jp), &
                    "q_sw must equal sw_thru_diag exactly (the coupler component fill)")
         if (allocated(error)) exit checks
         ! Q_heat_const = 0, only q_sw + heat_added nonzero => Q_heat is
         ! their sum: the shortwave appears EXACTLY ONCE.
         call check(error, abs(sf%Q_heat(ip, jp) &
                               - (ice%heat_flux_diag(ip, jp) + ice%sw_thru_diag(ip, jp))) &
                    <= 1.0e-12_wp*max(abs(sf%Q_heat(ip, jp)), 1.0_wp), &
                    "assembled Q_heat must be heat_flux_diag + sw_thru_diag (SW counted once)")
         if (allocated(error)) exit checks

         d_ocean_heat = sf%rho0*SEAWATER_CP*(hT_after - hT_before)
         rhs = ice%heat_to_ocn(ip, jp, 1) - ice%fb(ip, jp)*DT_THERM &
               + ice%sw_thru(ip, jp, 1)*DT_THERM
         scale = sf%rho0*SEAWATER_CP*abs(hT_before)
         tol = 1.0e-9_wp*max(scale, 1.0_wp)
         call check(error, abs(d_ocean_heat - rhs) <= tol, &
                    "components-on closure: rho0*cp*Delta(hT) == heat_to_ocn - fb*dt + sw_thru*dt")

      end block checks
      call teardown(ms, eos, ice, sf)
   end subroutine test_sw_components_q_sw

end module test_ocean_ice_driver_column
