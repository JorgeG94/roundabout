!! Analytic tests for the sea-ice snowfall source term (`m_snow`, PR 26):
!! `ice_snow_accumulate` (`rdb_ice_mass`), the `fprec`/`snow` threading
!! through `rdb_ice_column`, the atmospheric seam fill (`rdb_ice_atm_
!! forcing`), the pre-column ocean-share snapshot
!! (`ice_snow_part_ocn_fill_impl`, `rdb_ice_thermo_driver`), and the
!! ice-free-share ocean contributor (`ice_snowfall_ocean_share`,
!! `rdb_ice_snow`).
!!
!! Harness mirrors `test_ocean_ice_driver_column.F90`: grid 6x4, NGHOST=2,
!! NZ=3, H_LAYER=10, S_INIT=35, DT_THERM=1200. Cases 1-4 call
!! `ice_column_step`/`ice_snow_accumulate` directly (no GPU mapping
!! needed -- pure scalar/local-array arguments). Cases 5-9 drive the full
!! thermo-cadence chain through `run_chain_snow`, the PR-26 twin of that
!! file's `run_chain` (adds the has_snowfall-gated `ice_snow_part_ocn_
!! fill_impl` call inside `ice_thermo_driver_step` + the `ice_snowfall_
!! ocean_share` contributor, mandated order per `rdb_driver.F90`).
!!
!! **GPU mem:separate discipline** (per PR-3a/3b/3c + `test_open_
!! boundary_out_closes`): one shared enter_data/exit_data span around the
!! whole kernel chain; `!$acc update self` the touched arrays before any
!! host assertion; `atm_fprec`/`snow_part_ocn`/`fprec_ocn_diag` are mapped
!! by `ice%enter_data()` (`rdb_ice_state.F90`), so the existing `ice%
!! enter_data()`/`exit_data()` pair covers them too. Directives are inert
!! no-ops on host builds -- write them unconditionally; a green multicore
!! run proves nothing about device data motion.
!!
!! Cases (PLAN_PR26_snowfall.md S9):
!!   * `column_snow_accumulates_linearly` -- the roadmap's acceptance
!!     test: m_snow after N cold, non-melting windows == N*fprec*dt_therm
!!     to 1e-12 rel.
!!   * `snow_arrives_at_layer_enthalpy` -- SIS2's "the snow enthalpy
!!     should not have changed" convention, both as an exact algebraic
!!     identity on `ice_snow_accumulate` directly and as a loose
!!     (numerical-noise-tolerant) sanity check through the full column.
!!   * `new_snow_melts_before_ice` -- new snow is a real thermodynamic
!!     buffer IN THE SAME WINDOW it falls (the add is before the melt
!!     peels, not after).
!!   * `snow_reaches_snow_albedo` -- the roadmap's literal "done when":
!!     snow albedo (0.85 vs bare 0.5826) is reachable, and the melt point
!!     switches from T_f(S_ice) to 0 once snow is present.
!!   * `mass_budget_closes_ncat1` / `mass_budget_closes_ncat5` -- the
!!     roadmap's second acceptance test (snow mass + ocean-received mass
!!     = precipitated mass to round-off), at both area conventions
!!     (ncat==1 per-CELL, ncat>1 per-ICE-AREA).
!!   * `no_orphan_snow_created` -- guards the `rdb_ice_transport`
!!     fail-loud `mca_snow>0` where `mca_ice<=0` reduction.
!!   * `ocean_energy_and_salt_close` -- no silent energy leak (the
!!     `sw_thru` failure mode) through the REAL apply-tracers kernel.
!!   * `snowfall_disabled_bitident` -- the house default-off gate,
!!     against an independent reference chain that never references the
!!     new PR-26 procedures at all.
module test_ocean_ice_snowfall
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_eos, only: eos_t, eos_freezing_point
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t, ocean_surface_flux_apply_tracers, &
                                     SEAWATER_CP
   use rdb_ice_state, only: ocean_sea_ice_t
   use rdb_ice_frazil_uptake, only: ice_frazil_uptake
   use rdb_ice_atm_forcing, only: ice_atm_forcing_restoring
   use rdb_ice_basal_flux, only: ice_compute_basal_flux
   use rdb_ice_thermo_driver, only: ice_thermo_driver_step
   use rdb_ice_snow, only: ice_snowfall_ocean_share
   use rdb_ice_ocean_coupler, only: ice_ocean_brine_flux, ice_ocean_heat_flux
   use rdb_ice_mass, only: ice_snow_accumulate
   use rdb_ice_column, only: ice_column_step, ICE_BULK_SALINITY, ICE_RHO_ICE, ICE_RHO_SNOW
   use rdb_ice_enthalpy, only: ice_enth_from_ts, ice_t_freeze, ICE_LAT_FUS, ICE_CP_WATER
   use rdb_ice_optics, only: ice_optics_csim4
   implicit none
   private

   public :: collect_ocean_ice_snowfall_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NX_PHYS = 6
   integer, parameter :: NY_PHYS = 4
   integer, parameter :: NZ = 3
   integer, parameter :: NK = 2
   real(wp), parameter :: H_LAYER = 10.0_wp
      !! Uniform layer thickness (m).
   real(wp), parameter :: S_INIT = 35.0_wp
      !! Uniform salinity (PSU).
   real(wp), parameter :: T_DEEP = 1.0_wp
      !! Sub-surface IC (degC), untouched by every driver-chain case.
   real(wp), parameter :: DT_THERM = 1200.0_wp
      !! Thermo timestep (s).
   real(wp), parameter :: FPREC = 2.0e-5_wp
      !! Uniform test snowfall rate (kg/m^2/s) for the driver-chain cases
      !! (~0.024 kg/m^2 per DT_THERM window -- comfortably resolvable in
      !! double precision, well inside the &ocean_ice_nml snowfall min=0
      !! envelope).

contains

   subroutine collect_ocean_ice_snowfall_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("column_snow_accumulates_linearly", &
                               test_column_snow_accumulates_linearly), &
                  new_unittest("snow_arrives_at_layer_enthalpy", &
                               test_snow_arrives_at_layer_enthalpy), &
                  new_unittest("new_snow_melts_before_ice", test_new_snow_melts_before_ice), &
                  new_unittest("snow_reaches_snow_albedo", test_snow_reaches_snow_albedo), &
                  new_unittest("mass_budget_closes_ncat1", test_mass_budget_closes_ncat1), &
                  new_unittest("mass_budget_closes_ncat5", test_mass_budget_closes_ncat5), &
                  new_unittest("no_orphan_snow_created", test_no_orphan_snow_created), &
                  new_unittest("ocean_energy_and_salt_close", test_ocean_energy_and_salt_close), &
                  new_unittest("snowfall_disabled_bitident", test_snowfall_disabled_bitident) &
                  ]
   end subroutine collect_ocean_ice_snowfall_tests

   ! -----------------------------------------------------------------
   ! Shared helpers -- driver-chain cases (5-9)
   ! -----------------------------------------------------------------

   subroutine setup_state_snow(grid, ms, eos, ice, sf, ncat, has_snowfall, t_surface)
      !! Tiny ocean state, PR-26 twin of `test_ocean_ice_driver_column`'s
      !! `setup_state`: adds `ncat`/`has_snowfall` so the mass-budget
      !! cases can drive the SIS2 ITD mode (ncat>1) and the disabled-
      !! bitident case can latch the config gate directly.
      type(hgrid_t), intent(out) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(eos_t), intent(out) :: eos
      type(ocean_sea_ice_t), intent(inout) :: ice
      type(ocean_surface_flux_t), intent(inout) :: sf
      integer, intent(in) :: ncat
      logical, intent(in) :: has_snowfall
      real(wp), intent(in) :: t_surface

      call grid%init(NX_PHYS, NY_PHYS, NGHOST, 1.0_wp, 1.0_wp)
      ms%nz_ml = NZ
      call ms%init(grid)
      call eos%init(grid)
      ice%enable = .true.
      ice%ncat = ncat
      ice%has_snowfall = has_snowfall
      call ice%init(grid)
      call sf%init(grid)
      call sf%set_surface_flux_const(0.0_wp, 0.0_wp)

      ms%h_layer = H_LAYER
      ms%tracers(ms%idx_salinity)%hTr = S_INIT*H_LAYER
      ms%tracers(ms%idx_temperature)%hTr = T_DEEP*H_LAYER
      ms%tracers(ms%idx_temperature)%hTr(:, :, NZ) = t_surface*H_LAYER
   end subroutine setup_state_snow

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

   subroutine run_chain_snow(grid, eos, ms, ice, sf, air_temp, restore_lambda, sw_down, snowfall)
      !! Full thermo-cadence chain INCLUDING PR 26, mandated order
      !! (rdb_driver.F90): forcing -> basal -> frazil uptake -> column
      !! driver (which internally fills `snow_part_ocn` when
      !! `has_snowfall`) -> snowfall ocean share (gated on
      !! `has_snowfall`) -> brine coupler -> heat coupler -> the REAL
      !! apply_tracers kernel. One shared map/unmap spanning every kernel
      !! (GPU mem:separate discipline).
      type(hgrid_t), intent(in) :: grid
      type(eos_t), intent(in) :: eos
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_sea_ice_t), intent(inout) :: ice
      type(ocean_surface_flux_t), intent(inout) :: sf
      real(wp), intent(in) :: air_temp, restore_lambda, sw_down, snowfall

      !$acc enter data copyin(ms, sf)
      call ms%enter_data()
      call ice%enter_data()
      call sf%enter_data()

      call ice_atm_forcing_restoring(ice, air_temp, restore_lambda, sw_down, snowfall)
      call ice_compute_basal_flux(grid, eos, ms, ice, DT_THERM)
      call ice_frazil_uptake(grid, eos, ms, ice, DT_THERM)
      call ice_thermo_driver_step(grid, eos, ms, ice, DT_THERM)
      if (ice%has_snowfall) then
         call ice_snowfall_ocean_share(grid, ice, DT_THERM)
      end if
      call ice_ocean_brine_flux(sf, ice)
      call ice_ocean_heat_flux(sf, ice)
      call ocean_surface_flux_apply_tracers(grid, sf, ms, DT_THERM)

      associate (hT => ms%tracers(ms%idx_temperature)%hTr, &
                 hS => ms%tracers(ms%idx_salinity)%hTr, &
                 mi => ice%m_ice, msn => ice%m_snow, ei => ice%enth_ice, si => ice%sal_ice, &
                 mf => ice%m_frozen_diag, sd => ice%salt_flux_diag, &
                 hf => ice%heat_flux_diag, md => ice%m_melt_diag, &
                 fb => ice%fb, h2o => ice%h2o_ocn_to_ice, h2i => ice%h2o_ice_to_ocn, &
                 hto => ice%heat_to_ocn, qs => sf%Q_salt, qh => sf%Q_heat, &
                 af => ice%atm_fprec, spo => ice%snow_part_ocn, fod => ice%fprec_ocn_diag, &
                 psz => ice%part_size)
         !$acc update self(hT, hS, mi, msn, ei, si, mf, sd, hf, md, fb, h2o, h2i, hto, qs, qh, &
         !$acc              af, spo, fod, psz)
      end associate

      call sf%exit_data()
      call ice%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, sf)
   end subroutine run_chain_snow

   subroutine run_chain_pre26(grid, eos, ms, ice, sf, air_temp, restore_lambda, sw_down)
      !! Independent reference chain that NEVER references
      !! `ice_snowfall_ocean_share` or the has_snowfall branch at all --
      !! the literal pre-PR-26 driver sequence (with the now-mandatory
      !! `snowfall=0.0` argument to `ice_atm_forcing_restoring`, the only
      !! textual difference). Used by `test_snowfall_disabled_bitident`
      !! as a true independent baseline: if a future change accidentally
      !! removes the `has_snowfall` gate in `run_chain_snow`, THIS
      !! subroutine still gives a different answer to diff against.
      type(hgrid_t), intent(in) :: grid
      type(eos_t), intent(in) :: eos
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_sea_ice_t), intent(inout) :: ice
      type(ocean_surface_flux_t), intent(inout) :: sf
      real(wp), intent(in) :: air_temp, restore_lambda, sw_down

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
      call ocean_surface_flux_apply_tracers(grid, sf, ms, DT_THERM)

      associate (hT => ms%tracers(ms%idx_temperature)%hTr, &
                 hS => ms%tracers(ms%idx_salinity)%hTr, &
                 mi => ice%m_ice, msn => ice%m_snow, &
                 qs => sf%Q_salt, qh => sf%Q_heat)
         !$acc update self(hT, hS, mi, msn, qs, qh)
      end associate

      call sf%exit_data()
      call ice%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, sf)
   end subroutine run_chain_pre26

   subroutine seed_cold_pack(ice, ip, jp)
      !! Seed an ice cap WITH pre-existing snow, so the disabled-bitident
      !! comparison is non-trivial (proves existing snow rides through
      !! the chain unperturbed too, not just that no NEW snow appears).
      type(ocean_sea_ice_t), intent(inout) :: ice
      integer, intent(in) :: ip, jp
      ice%m_ice(ip, jp, 1) = 0.5_wp*ICE_RHO_ICE
      ice%m_snow(ip, jp, 1) = 1.0_wp
      ice%enth_snow(ip, jp, 1, 1) = ice_enth_from_ts(-5.0_wp, 0.0_wp)
      ice%enth_ice(ip, jp, 1, :) = ice_enth_from_ts(-10.0_wp, ICE_BULK_SALINITY)
      ice%sal_ice(ip, jp, 1, :) = ICE_BULK_SALINITY
   end subroutine seed_cold_pack

   ! -----------------------------------------------------------------
   ! Cases 1-4: direct ice_column_step / ice_snow_accumulate calls
   ! -----------------------------------------------------------------

   subroutine test_column_snow_accumulates_linearly(error)
      !! The roadmap's acceptance test: m_snow after N calls on a cold,
      !! non-melting column == N*snow_step to 1e-12 rel. Cold stiff-SEB
      !! pin (tsurf_target=-10, dsf_dt=1e6) + fb=0 + sst=tfw guarantees
      !! tmelt=bmelt=0 identically (mirrors the proven `test_stefan_
      !! growth`/`equilibrium_fixed_point` cold-pin pattern in
      !! `test_ocean_ice_column.F90`), so no resize branch can destroy
      !! the accumulated snow.
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: TSURF_TARGET = -10.0_wp
      integer, parameter :: N_WINDOWS = 10
      real(wp) :: m_snow, m_ice_tot, enth_snow_pt
      real(wp) :: enth_ice_bu(NK), sal_ice_bu(NK)
      real(wp) :: sf_0, dsf_dt, sw_dn, tfw, fb, sst, s_surf, snow_step
      real(wp) :: tsurf, h2o_ocn_to_ice, h2o_ice_to_ocn, heat_to_ocn, sw_thru, snow_to_ice
      real(wp) :: expected
      integer :: step

      sf_0 = -1.0e6_wp*TSURF_TARGET
      dsf_dt = 1.0e6_wp
      sw_dn = 0.0_wp
      fb = 0.0_wp
      tfw = ice_t_freeze(S_INIT)
      sst = tfw
      s_surf = S_INIT

      m_snow = 0.0_wp
      m_ice_tot = 0.5_wp*ICE_RHO_ICE
      enth_snow_pt = ice_enth_from_ts(0.0_wp, 0.0_wp)
      enth_ice_bu(:) = ice_enth_from_ts(TSURF_TARGET, ICE_BULK_SALINITY)
      sal_ice_bu(:) = ICE_BULK_SALINITY

      snow_step = FPREC*DT_THERM

      do step = 1, N_WINDOWS
         call ice_column_step(NK, m_snow, m_ice_tot, enth_snow_pt, enth_ice_bu, sal_ice_bu, &
                              sf_0, dsf_dt, sw_dn, tfw, fb, sst, s_surf, DT_THERM, .false., &
                              snow_step, tsurf, h2o_ocn_to_ice, h2o_ice_to_ocn, &
                              heat_to_ocn, sw_thru, snow_to_ice)
      end do

      expected = real(N_WINDOWS, wp)*snow_step
      call check(error, abs(m_snow - expected) < 1.0e-12_wp*expected, &
                 "m_snow must accumulate linearly to 1e-12 rel over N non-melting windows")
   end subroutine test_column_snow_accumulates_linearly

   subroutine test_snow_arrives_at_layer_enthalpy(error)
      !! SIS2's "the snow enthalpy should not have changed" convention.
      !! Part (a): exact algebraic identity directly on `ice_snow_
      !! accumulate` -- the routine does not take enthalpy as an
      !! argument at all, so this pins the CONTRACT (a future "helpful"
      !! rewrite that mixes in an atmospheric temperature would have to
      !! change the signature, breaking every call site). Part (b): one
      !! full `ice_column_step` call on a genuinely isothermal, zero-flux
      !! column (sf_0=dsf_dt=sw_dn=fb=0, uniform T everywhere including
      !! the pre-existing snow) -- a true numerical fixed point, so
      !! `enth_snow_pt` returning unchanged is a meaningful integration
      !! check, not just unit-level algebra.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: m_lay(0:NK), enthalpy(0:NK)
      real(wp) :: enth_snow0, snow_amt, col_before, col_after
      real(wp) :: t_iso
      real(wp) :: m_snow, m_ice_tot, enth_snow_pt
      real(wp) :: enth_ice_bu(NK), sal_ice_bu(NK)
      real(wp) :: sf_0, dsf_dt, sw_dn, tfw, fb, sst, s_surf
      real(wp) :: tsurf, h2o_ocn_to_ice, h2o_ice_to_ocn, heat_to_ocn, sw_thru, snow_to_ice
      integer :: k

      ! ---- (a) direct, exact ----
      enth_snow0 = ice_enth_from_ts(-5.0_wp, 0.0_wp)
      m_lay(0) = 3.0_wp
      enthalpy(0) = enth_snow0
      do k = 1, NK
         m_lay(k) = 100.0_wp
         enthalpy(k) = ice_enth_from_ts(-8.0_wp, ICE_BULK_SALINITY)
      end do
      snow_amt = FPREC*DT_THERM
      col_before = sum(m_lay(0:NK)*enthalpy(0:NK))

      call ice_snow_accumulate(NK, m_lay, snow_amt)

      call check(error, enthalpy(0) == enth_snow0, &
                 "ice_snow_accumulate must not touch enthalpy(0) at all")
      if (allocated(error)) return

      col_after = sum(m_lay(0:NK)*enthalpy(0:NK))
      call check(error, abs((col_after - col_before) - snow_amt*enth_snow0) < &
                 1.0e-13_wp*abs(snow_amt*enth_snow0), &
                 "column enthalpy must grow by exactly snow*enth_snow0")
      if (allocated(error)) return

      ! ---- (b) integration, loose tolerance (numerical fixed point) ----
      ! The column's bottom BC is ALWAYS Dirichlet-pinned at tfw (not the
      ! seeded sst) -- so a genuine zero-flux fixed point needs the WHOLE
      ! column (snow + ice) isothermal AT tfw, not an arbitrary cold
      ! temperature (an isothermal seed away from tfw creates a real,
      ! if small, one-step basal conduction transient -- verified
      ! numerically while developing this test).
      sf_0 = 0.0_wp
      dsf_dt = 0.0_wp
      sw_dn = 0.0_wp
      fb = 0.0_wp
      tfw = ice_t_freeze(S_INIT)
      t_iso = tfw
      sst = t_iso
      s_surf = S_INIT

      m_snow = 4.0_wp
      enth_snow_pt = ice_enth_from_ts(t_iso, 0.0_wp)
      m_ice_tot = 0.5_wp*ICE_RHO_ICE
      enth_ice_bu(:) = ice_enth_from_ts(t_iso, ICE_BULK_SALINITY)
      sal_ice_bu(:) = ICE_BULK_SALINITY

      call check(error, enth_snow_pt /= 0.0_wp, "sanity: seeded enth_snow_pt must be nonzero")
      if (allocated(error)) return

      block
         real(wp) :: enth_snow_before
         enth_snow_before = enth_snow_pt
         call ice_column_step(NK, m_snow, m_ice_tot, enth_snow_pt, enth_ice_bu, sal_ice_bu, &
                              sf_0, dsf_dt, sw_dn, tfw, fb, sst, s_surf, DT_THERM, .false., &
                              snow_amt, tsurf, h2o_ocn_to_ice, h2o_ice_to_ocn, &
                              heat_to_ocn, sw_thru, snow_to_ice)
         call check(error, abs(enth_snow_pt - enth_snow_before) < 1.0e-7_wp*abs(enth_snow_before), &
                    "on an isothermal (at tfw) zero-flux column, snow enthalpy must stay "// &
                    "(near-)unchanged through the full column step")
      end block
   end subroutine test_snow_arrives_at_layer_enthalpy

   subroutine test_new_snow_melts_before_ice(error)
      !! New snow is meltable in the SAME window it falls: a warm-restoring,
      !! snow-free-entering column with a generous new snowfall (50 kg/m2)
      !! must lose SOME of that new snow (proving the add precedes the
      !! melt peels) while leaving m_ice untouched (the peel loop must
      !! stop inside the snow layer). Uses a REALISTIC restoring
      !! coefficient (lambda=20 W/m2/K, the driver's own convention,
      !! SF(T)=lambda*(T-air_temp)) -- a huge "stiff pin" dsf_dt (as used
      !! elsewhere to hold tsurf tightly for a NON-melting experiment)
      !! turns out to inject an effectively unbounded tmelt once the
      !! clamp fires, melting through everything in one window; that is
      !! a modelling artifact of the pin technique, not something this
      !! case should exercise. The ice column is seeded ISOTHERMAL AT
      !! tfw (the column's bottom BC) so the base carries no spurious
      !! one-step conduction transient -- verified numerically while
      !! developing this test (an isothermal seed away from tfw leaves a
      !! small but real bottom-melt residual unrelated to the top-melt
      !! question this case is about).
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: AIR_TEMP = 20.0_wp
      real(wp), parameter :: LAMBDA = 20.0_wp
      real(wp), parameter :: SNOW_AMT = 50.0_wp
      real(wp) :: m_snow, m_ice_tot, m_ice_before, enth_snow_pt
      real(wp) :: enth_ice_bu(NK), sal_ice_bu(NK)
      real(wp) :: sf_0, dsf_dt, sw_dn, tfw, fb, sst, s_surf
      real(wp) :: tsurf, h2o_ocn_to_ice, h2o_ice_to_ocn, heat_to_ocn, sw_thru, snow_to_ice

      sf_0 = -LAMBDA*AIR_TEMP
      dsf_dt = LAMBDA
      sw_dn = 0.0_wp
      fb = 0.0_wp
      tfw = ice_t_freeze(S_INIT)
      sst = tfw
      s_surf = S_INIT

      m_snow = 0.0_wp
      m_ice_tot = 0.5_wp*ICE_RHO_ICE
      m_ice_before = m_ice_tot
      enth_snow_pt = ice_enth_from_ts(0.0_wp, 0.0_wp)
      enth_ice_bu(:) = ice_enth_from_ts(tfw, ICE_BULK_SALINITY)
      sal_ice_bu(:) = ICE_BULK_SALINITY

      call ice_column_step(NK, m_snow, m_ice_tot, enth_snow_pt, enth_ice_bu, sal_ice_bu, &
                           sf_0, dsf_dt, sw_dn, tfw, fb, sst, s_surf, DT_THERM, .false., &
                           SNOW_AMT, tsurf, h2o_ocn_to_ice, h2o_ice_to_ocn, &
                           heat_to_ocn, sw_thru, snow_to_ice)

      call check(error, m_snow < SNOW_AMT, &
                 "some of the new snow must melt this window (proves add-before-peel)")
      if (allocated(error)) return
      call check(error, m_snow > 0.0_wp, &
                 "snow must not fully melt through into ice (test margin too tight)")
      if (allocated(error)) return
      call check(error, abs(m_ice_tot - m_ice_before) < 1.0e-6_wp*m_ice_before, &
                 "ice mass must stay (numerically) unchanged -- the peel must stop "// &
                 "inside the snow layer")
   end subroutine test_new_snow_melts_before_ice

   subroutine test_snow_reaches_snow_albedo(error)
      !! The roadmap's literal "done when": snow albedo (0.85 vs bare
      !! 0.5826) is reachable, and the melt point switches to 0 once snow
      !! is present. Part A accumulates a thick snowpack (5 kg/m2 x 10
      !! cold windows = 50 kg/m2, well past the ICE_SNOW_PATCH=0.02 m
      !! masking depth) then checks the CSIM4 albedo directly. Part B
      !! restores BOTH a bare and a snow-covered column, isothermal AT
      !! tfw (the column's bottom BC -- eliminates a spurious one-step
      !! basal transient, same rationale as `test_new_snow_melts_
      !! before_ice`), toward a target (-0.1 degC) BETWEEN bare ice's
      !! T_f(S_ice)=-0.216 and snow's 0: the bare column must melt, the
      !! snow-covered column must not -- the observable signature of tsf
      !! switching from T_f(S_ice) to 0. Realistic restoring (lambda=300
      !! W/m2/K; see `test_new_snow_melts_before_ice`'s docstring for why
      !! a "stiff pin" dsf_dt is the wrong tool once a clamp can fire).
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: TSURF_COLD = -10.0_wp
      real(wp), parameter :: AIR_TEMP_NEAR0 = -0.1_wp
      real(wp), parameter :: LAMBDA_NEAR0 = 300.0_wp
      real(wp), parameter :: SNOW_STEP = 5.0_wp
      integer, parameter :: N_WINDOWS = 10
      real(wp) :: m_snow, m_ice_tot, enth_snow_pt
      real(wp) :: enth_ice_bu(NK), sal_ice_bu(NK)
      real(wp) :: sf_0, dsf_dt, sw_dn, tfw, fb, sst, s_surf
      real(wp) :: tsurf, h2o_ocn_to_ice, h2o_ice_to_ocn, heat_to_ocn, sw_thru, snow_to_ice
      real(wp) :: albedo_snow, albedo_bare, abs_sfc, abs_snow, abs_ocn, abs_int, pen
      real(wp) :: abs_ice_lay(NK)
      real(wp) :: m_snow_accum
      real(wp) :: h2o_ice_to_ocn_bare, h2o_ice_to_ocn_snow
      integer :: step

      ! ---- Part A: accumulate a thick snowpack under cold forcing ----
      sf_0 = -1.0e6_wp*TSURF_COLD
      dsf_dt = 1.0e6_wp
      sw_dn = 0.0_wp
      fb = 0.0_wp
      tfw = ice_t_freeze(S_INIT)
      sst = tfw
      s_surf = S_INIT

      m_snow = 0.0_wp
      m_ice_tot = 0.5_wp*ICE_RHO_ICE
      enth_snow_pt = ice_enth_from_ts(0.0_wp, 0.0_wp)
      enth_ice_bu(:) = ice_enth_from_ts(TSURF_COLD, ICE_BULK_SALINITY)
      sal_ice_bu(:) = ICE_BULK_SALINITY

      do step = 1, N_WINDOWS
         call ice_column_step(NK, m_snow, m_ice_tot, enth_snow_pt, enth_ice_bu, sal_ice_bu, &
                              sf_0, dsf_dt, sw_dn, tfw, fb, sst, s_surf, DT_THERM, .false., &
                              SNOW_STEP, tsurf, h2o_ocn_to_ice, h2o_ice_to_ocn, &
                              heat_to_ocn, sw_thru, snow_to_ice)
      end do
      m_snow_accum = m_snow

      call ice_optics_csim4(NK, m_snow_accum/ICE_RHO_SNOW, m_ice_tot/ICE_RHO_ICE, TSURF_COLD, &
                            ICE_BULK_SALINITY, albedo_snow, abs_sfc, abs_snow, abs_ice_lay, &
                            abs_ocn, abs_int, pen)
      call ice_optics_csim4(NK, 0.0_wp, m_ice_tot/ICE_RHO_ICE, TSURF_COLD, &
                            ICE_BULK_SALINITY, albedo_bare, abs_sfc, abs_snow, abs_ice_lay, &
                            abs_ocn, abs_int, pen)

      call check(error, albedo_snow > 0.80_wp, "snow-covered albedo must exceed 0.80")
      if (allocated(error)) return
      call check(error, albedo_bare < 0.60_wp, "bare-ice control albedo must be below 0.60")
      if (allocated(error)) return

      ! ---- Part B: the snow-covered melt point (tsf=0) vs bare ----
      sf_0 = -LAMBDA_NEAR0*AIR_TEMP_NEAR0
      dsf_dt = LAMBDA_NEAR0

      ! BARE: seeded isothermal AT tfw (the bottom BC) -- no snow.
      m_snow = 0.0_wp
      m_ice_tot = 0.5_wp*ICE_RHO_ICE
      enth_snow_pt = ice_enth_from_ts(0.0_wp, 0.0_wp)
      enth_ice_bu(:) = ice_enth_from_ts(tfw, ICE_BULK_SALINITY)
      sal_ice_bu(:) = ICE_BULK_SALINITY
      call ice_column_step(NK, m_snow, m_ice_tot, enth_snow_pt, enth_ice_bu, sal_ice_bu, &
                           sf_0, dsf_dt, sw_dn, tfw, fb, sst, s_surf, DT_THERM, .false., &
                           0.0_wp, tsurf, h2o_ocn_to_ice, h2o_ice_to_ocn, &
                           heat_to_ocn, sw_thru, snow_to_ice)
      h2o_ice_to_ocn_bare = h2o_ice_to_ocn

      ! SNOW-COVERED: same isothermal-at-tfw ice, pre-existing snow cover
      ! (also at tfw) -- only difference from BARE is m_snow > 0 on entry.
      m_snow = m_snow_accum
      m_ice_tot = 0.5_wp*ICE_RHO_ICE
      enth_snow_pt = ice_enth_from_ts(tfw, 0.0_wp)
      enth_ice_bu(:) = ice_enth_from_ts(tfw, ICE_BULK_SALINITY)
      sal_ice_bu(:) = ICE_BULK_SALINITY
      call ice_column_step(NK, m_snow, m_ice_tot, enth_snow_pt, enth_ice_bu, sal_ice_bu, &
                           sf_0, dsf_dt, sw_dn, tfw, fb, sst, s_surf, DT_THERM, .false., &
                           0.0_wp, tsurf, h2o_ocn_to_ice, h2o_ice_to_ocn, &
                           heat_to_ocn, sw_thru, snow_to_ice)
      h2o_ice_to_ocn_snow = h2o_ice_to_ocn

      call check(error, h2o_ice_to_ocn_bare > 1.0e-4_wp, &
                 "bare ice must melt: restoring target (-0.1) exceeds T_f(S_ice) (-0.216)")
      if (allocated(error)) return
      call check(error, h2o_ice_to_ocn_snow < 1.0e-6_wp, &
                 "snow-covered ice must NOT (meaningfully) melt: the melt point has "// &
                 "switched to 0 degC")
   end subroutine test_snow_reaches_snow_albedo

   ! -----------------------------------------------------------------
   ! Cases 5-9: full thermo-cadence driver chain
   ! -----------------------------------------------------------------

   subroutine test_mass_budget_closes_ncat1(error)
      !! ncat==1 (legacy lumped, per-CELL area). One iced cell, one
      !! ice-free wet cell, one land cell, under a fully passive thermo
      !! window (restore_lambda=0 => SF(T)=0 identically; sst=T_f =>
      !! fb=0) so the ONLY thing that happens is the snow source term.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      type(ocean_surface_flux_t) :: sf
      real(wp) :: t_f, snow_sum, ocn_sum, expected, tol
      integer :: ip_ice, jp_ice, ip_free, jp_free, ip_land, jp_land
      integer :: i, j, i_lo, i_hi, j_lo, j_hi, n_wet

      checks: block
         call setup_state_snow(grid, ms, eos, ice, sf, 1, .true., 0.0_wp)
         t_f = eos_freezing_point(eos, S_INIT, 0.0_wp)
         ms%tracers(ms%idx_temperature)%hTr(:, :, NZ) = t_f*H_LAYER

         ip_ice = NGHOST + 2; jp_ice = NGHOST + 2
         ip_free = NGHOST + 3; jp_free = NGHOST + 2
         ip_land = NGHOST + 4; jp_land = NGHOST + 2
         ms%wet_mask(ip_land, jp_land) = 0.0_wp

         ice%m_ice(ip_ice, jp_ice, 1) = 0.5_wp*ICE_RHO_ICE
         ice%m_snow(ip_ice, jp_ice, 1) = 0.0_wp
         ice%enth_snow(ip_ice, jp_ice, 1, 1) = ice_enth_from_ts(0.0_wp, 0.0_wp)
         ice%enth_ice(ip_ice, jp_ice, 1, :) = ice_enth_from_ts(-5.0_wp, ICE_BULK_SALINITY)
         ice%sal_ice(ip_ice, jp_ice, 1, :) = ICE_BULK_SALINITY

         call run_chain_snow(grid, eos, ms, ice, sf, air_temp=0.0_wp, restore_lambda=0.0_wp, &
                             sw_down=0.0_wp, snowfall=FPREC)

         expected = FPREC*DT_THERM
         tol = 1.0e-12_wp*expected
         call check(error, abs(ice%m_snow(ip_ice, jp_ice, 1) - expected) < tol, &
                    "iced cell m_snow must grow by exactly FPREC*DT_THERM")
         if (allocated(error)) exit checks
         call check(error, ice%fprec_ocn_diag(ip_ice, jp_ice) == 0.0_wp, &
                    "fully iced (ncat=1) cell must receive zero ocean-bound share")
         if (allocated(error)) exit checks
         call check(error, abs(ice%fprec_ocn_diag(ip_free, jp_free) - FPREC) < 1.0e-14_wp*FPREC, &
                    "ice-free cell must receive the full FPREC as ocean-bound share")
         if (allocated(error)) exit checks
         call check(error, ice%fprec_ocn_diag(ip_land, jp_land) == 0.0_wp, &
                    "land cell must receive zero ocean-bound share")
         if (allocated(error)) exit checks

         i_lo = NGHOST + 1; i_hi = grid%nx_total - NGHOST
         j_lo = NGHOST + 1; j_hi = grid%ny_total - NGHOST
         snow_sum = 0.0_wp
         ocn_sum = 0.0_wp
         n_wet = 0
         do j = j_lo, j_hi
            do i = i_lo, i_hi
               snow_sum = snow_sum + ice%m_snow(i, j, 1)
               ocn_sum = ocn_sum + ice%fprec_ocn_diag(i, j)*DT_THERM
               if (ms%wet_mask(i, j) > 0.5_wp) n_wet = n_wet + 1
            end do
         end do
         expected = FPREC*DT_THERM*real(n_wet, wp)
         tol = 1.0e-10_wp*expected
         call check(error, abs((snow_sum + ocn_sum) - expected) < tol, &
                    "domain identity: snow accumulated + ocean share == FPREC*DT_THERM*n_wet")
      end block checks
      call teardown(ms, eos, ice, sf)
   end subroutine test_mass_budget_closes_ncat1

   subroutine test_mass_budget_closes_ncat5(error)
      !! ncat==5 (SIS2 ITD, per-ICE-AREA). part_size(0)=0.4 (open water),
      !! part_size(2)=0.6 (the only occupied category); cats 1/3/4/5
      !! empty. The column must add fprec*dt UNWEIGHTED to cat 2's
      !! m_snow (per-ICE-AREA, not part-weighted); only the CELL-TOTAL
      !! bookkeeping takes the part_size weight.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      type(ocean_surface_flux_t) :: sf
      real(wp) :: t_f, expected, tol, cell_total
      integer :: ip, jp, cat
      integer, parameter :: NCAT5 = 5

      checks: block
         call setup_state_snow(grid, ms, eos, ice, sf, NCAT5, .true., 0.0_wp)
         t_f = eos_freezing_point(eos, S_INIT, 0.0_wp)
         ms%tracers(ms%idx_temperature)%hTr(:, :, NZ) = t_f*H_LAYER

         ip = NGHOST + 2; jp = NGHOST + 2

         ice%part_size(ip, jp, 0) = 0.4_wp
         ice%part_size(ip, jp, 1) = 0.0_wp
         ice%part_size(ip, jp, 2) = 0.6_wp
         ice%part_size(ip, jp, 3) = 0.0_wp
         ice%part_size(ip, jp, 4) = 0.0_wp
         ice%part_size(ip, jp, 5) = 0.0_wp
         ice%m_ice(ip, jp, 2) = 0.5_wp*ICE_RHO_ICE
         ice%m_snow(ip, jp, :) = 0.0_wp
         ice%enth_snow(ip, jp, 2, 1) = ice_enth_from_ts(0.0_wp, 0.0_wp)
         ice%enth_ice(ip, jp, 2, :) = ice_enth_from_ts(-5.0_wp, ICE_BULK_SALINITY)
         ice%sal_ice(ip, jp, 2, :) = ICE_BULK_SALINITY

         call run_chain_snow(grid, eos, ms, ice, sf, air_temp=0.0_wp, restore_lambda=0.0_wp, &
                             sw_down=0.0_wp, snowfall=FPREC)

         expected = FPREC*DT_THERM
         tol = 1.0e-12_wp*expected
         call check(error, abs(ice%m_snow(ip, jp, 2) - expected) < tol, &
                    "occupied category's m_snow must grow by FPREC*DT_THERM UNWEIGHTED")
         if (allocated(error)) exit checks

         do cat = 1, NCAT5
            if (cat == 2) cycle
            call check(error, ice%m_snow(ip, jp, cat) == 0.0_wp, &
                       "empty categories must receive exactly zero snow")
            if (allocated(error)) exit checks
         end do

         tol = 1.0e-13_wp*(0.4_wp*FPREC)
         call check(error, abs(ice%fprec_ocn_diag(ip, jp) - 0.4_wp*FPREC) < tol, &
                    "ocean-bound share must be the 0.4 open-water fraction of FPREC")
         if (allocated(error)) exit checks

         cell_total = 0.0_wp
         do cat = 1, NCAT5
            cell_total = cell_total + ice%part_size(ip, jp, cat)*ice%m_snow(ip, jp, cat)
         end do
         cell_total = cell_total + ice%fprec_ocn_diag(ip, jp)*DT_THERM
         tol = 1.0e-10_wp*expected
         call check(error, abs(cell_total - expected) < tol, &
                    "cell-total identity: part-weighted snow + ocean share == FPREC*DT_THERM")
      end block checks
      call teardown(ms, eos, ice, sf)
   end subroutine test_mass_budget_closes_ncat5

   subroutine test_no_orphan_snow_created(error)
      !! ncat=5, category 3 has part_size>0 (SIS2 would gate on this and
      !! add snow -- the orphan-snow trap) but ZERO ice mass -- fails the
      !! column's own entry gate. Guards a hard `error stop`:
      !! `rdb_ice_transport.F90` aborts the run on `mca_snow>0` where
      !! `mca_ice<=0`.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      type(ocean_surface_flux_t) :: sf
      real(wp) :: t_f
      integer :: ip, jp
      integer, parameter :: NCAT5 = 5

      checks: block
         call setup_state_snow(grid, ms, eos, ice, sf, NCAT5, .true., 0.0_wp)
         t_f = eos_freezing_point(eos, S_INIT, 0.0_wp)
         ms%tracers(ms%idx_temperature)%hTr(:, :, NZ) = t_f*H_LAYER

         ip = NGHOST + 2; jp = NGHOST + 2

         ! Category 3: part_size > 0 (SIS2's own gate) but m_ice == 0
         ! (already the init default) -- fails the column's entry gate.
         ice%part_size(ip, jp, 0) = 0.8_wp
         ice%part_size(ip, jp, 1) = 0.0_wp
         ice%part_size(ip, jp, 2) = 0.0_wp
         ice%part_size(ip, jp, 3) = 0.2_wp
         ice%part_size(ip, jp, 4) = 0.0_wp
         ice%part_size(ip, jp, 5) = 0.0_wp

         call run_chain_snow(grid, eos, ms, ice, sf, air_temp=0.0_wp, restore_lambda=0.0_wp, &
                             sw_down=0.0_wp, snowfall=FPREC)

         call check(error, ice%m_snow(ip, jp, 3) == 0.0_wp, &
                    "category with part_size>0 but m_ice==0 must receive ZERO snow "// &
                    "(the SIS2 orphan-snow trap)")
         if (allocated(error)) exit checks
         call check(error, abs(ice%fprec_ocn_diag(ip, jp) - FPREC) < 1.0e-13_wp*FPREC, &
                    "the whole cell's snowfall (incl. the would-be-orphaned 0.2 share) "// &
                    "must reach the ocean, not vanish")
      end block checks
      call teardown(ms, eos, ice, sf)
   end subroutine test_no_orphan_snow_created

   subroutine test_ocean_energy_and_salt_close(error)
      !! No silent energy leak (the sw_thru failure mode) through the
      !! REAL apply-tracers kernel, on an ice-free wet cell. (a) and (b)
      !! are asserted SEPARATELY on purpose: ICE_CP_WATER (4200) and
      !! SEAWATER_CP (3992) are deliberately distinct constants
      !! (rdb_ice_enthalpy docstring, "do not unify them").
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      type(ocean_surface_flux_t) :: sf
      real(wp), parameter :: SST = 5.0_wp
      real(wp) :: m_ocn, rhs, tol
      real(wp) :: hT_before, hT_after, hS_before, hS_after
      integer :: ip, jp

      checks: block
         call setup_state_snow(grid, ms, eos, ice, sf, 1, .true., SST)
         ip = NGHOST + 2; jp = NGHOST + 2

         hT_before = sum(ms%tracers(ms%idx_temperature)%hTr(ip, jp, :))
         hS_before = sum(ms%tracers(ms%idx_salinity)%hTr(ip, jp, :))

         call run_chain_snow(grid, eos, ms, ice, sf, air_temp=0.0_wp, restore_lambda=0.0_wp, &
                             sw_down=0.0_wp, snowfall=FPREC)

         hT_after = sum(ms%tracers(ms%idx_temperature)%hTr(ip, jp, :))
         hS_after = sum(ms%tracers(ms%idx_salinity)%hTr(ip, jp, :))

         m_ocn = ice%fprec_ocn_diag(ip, jp)*DT_THERM
         call check(error, m_ocn > 0.0_wp, "sanity: the ice-free cell must receive ocean-bound snow")
         if (allocated(error)) exit checks

         ! (a) heat_flux_diag*dt == -m_ocn*(Cp_water*sst + L_f)
         rhs = -m_ocn*(ICE_CP_WATER*SST + ICE_LAT_FUS)
         tol = 1.0e-12_wp*abs(rhs)
         call check(error, abs(ice%heat_flux_diag(ip, jp)*DT_THERM - rhs) < tol, &
                    "heat_flux_diag*dt must equal -m_ocn*(Cp_water*sst + L_f) to 1e-12 rel")
         if (allocated(error)) exit checks

         ! (b) rho0*SEAWATER_CP*Delta(hT) == heat_flux_diag*dt. Tolerance
         ! 1e-9 rel (not 1e-12): this identity round-trips through
         ! Q_heat = heat_flux_diag then hT += Q_heat*dt/(rho0*cp) then
         ! back through rho0*cp -- a few extra float ops beyond (a)'s
         ! direct identity, so it closes to ~1e-12 REL not exact bitwise;
         ! 1e-9 keeps 3 orders of margin while still proving no leak.
         rhs = ice%heat_flux_diag(ip, jp)*DT_THERM
         tol = 1.0e-9_wp*max(abs(rhs), 1.0_wp)
         call check(error, abs(sf%rho0*SEAWATER_CP*(hT_after - hT_before) - rhs) < tol, &
                    "ocean heat-content change must equal heat_flux_diag*dt to 1e-9 rel")
         if (allocated(error)) exit checks

         ! (c) rho0*Delta(hS) == -m_ocn*s_surf (same round-trip margin as (b))
         rhs = -m_ocn*S_INIT
         tol = 1.0e-9_wp*abs(rhs)
         call check(error, abs(sf%rho0*(hS_after - hS_before) - rhs) < tol, &
                    "ocean salt-content change must equal -m_ocn*s_surf to 1e-9 rel")
         if (allocated(error)) exit checks

         call check(error, sf%Q_heat(ip, jp) < 0.0_wp, "Q_heat must be strictly negative (cools)")
         if (allocated(error)) exit checks
         call check(error, sf%Q_salt(ip, jp) < 0.0_wp, "Q_salt must be strictly negative (freshens)")
      end block checks
      call teardown(ms, eos, ice, sf)
   end subroutine test_ocean_energy_and_salt_close

   subroutine test_snowfall_disabled_bitident(error)
      !! The house default-off gate: `snowfall=0` (has_snowfall=.false.)
      !! must be BIT-IDENTICAL to an independent reference chain that
      !! never references the new PR-26 procedures. Cold pack (some
      !! pre-existing snow + ice) under nontrivial forcing (not a
      !! passive no-op window), so real thermodynamics are exercised on
      !! both sides.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid_a, grid_b
      type(multilayer_state_t) :: ms_a, ms_b
      type(eos_t) :: eos_a, eos_b
      type(ocean_sea_ice_t) :: ice_a, ice_b
      type(ocean_surface_flux_t) :: sf_a, sf_b
      integer :: ip, jp

      checks: block
         call setup_state_snow(grid_a, ms_a, eos_a, ice_a, sf_a, 1, .false., -5.0_wp)
         call setup_state_snow(grid_b, ms_b, eos_b, ice_b, sf_b, 1, .false., -5.0_wp)

         ip = NGHOST + 2; jp = NGHOST + 2
         call seed_cold_pack(ice_a, ip, jp)
         call seed_cold_pack(ice_b, ip, jp)

         call run_chain_snow(grid_a, eos_a, ms_a, ice_a, sf_a, air_temp=-10.0_wp, &
                             restore_lambda=20.0_wp, sw_down=0.0_wp, snowfall=0.0_wp)
         call run_chain_pre26(grid_b, eos_b, ms_b, ice_b, sf_b, air_temp=-10.0_wp, &
                              restore_lambda=20.0_wp, sw_down=0.0_wp)

         call check(error,.not. ice_a%has_snowfall, "has_snowfall must be false by default")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(ice_a%atm_fprec)) == 0.0_wp, &
                    "atm_fprec must be exactly 0 (snowfall=0 filler)")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(ice_a%snow_part_ocn)) == 0.0_wp, &
                    "snow_part_ocn must be exactly 0 (never filled -- has_snowfall false)")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(ice_a%fprec_ocn_diag)) == 0.0_wp, &
                    "fprec_ocn_diag must be exactly 0 (never filled -- has_snowfall false)")
         if (allocated(error)) exit checks

         call check(error, all(ice_a%m_snow == ice_b%m_snow), &
                    "m_snow must be bit-identical to the PR-26-free reference chain")
         if (allocated(error)) exit checks
         call check(error, all(ice_a%m_ice == ice_b%m_ice), &
                    "m_ice must be bit-identical to the PR-26-free reference chain")
         if (allocated(error)) exit checks
         call check(error, all(sf_a%Q_heat == sf_b%Q_heat), &
                    "Q_heat must be bit-identical to the PR-26-free reference chain")
         if (allocated(error)) exit checks
         call check(error, all(sf_a%Q_salt == sf_b%Q_salt), &
                    "Q_salt must be bit-identical to the PR-26-free reference chain")
         if (allocated(error)) exit checks
         call check(error, all(ms_a%tracers(ms_a%idx_temperature)%hTr == &
                               ms_b%tracers(ms_b%idx_temperature)%hTr), &
                    "hT must be bit-identical to the PR-26-free reference chain")
         if (allocated(error)) exit checks
         call check(error, all(ms_a%tracers(ms_a%idx_salinity)%hTr == &
                               ms_b%tracers(ms_b%idx_salinity)%hTr), &
                    "hS must be bit-identical to the PR-26-free reference chain")
      end block checks
      call teardown(ms_a, eos_a, ice_a, sf_a)
      call teardown(ms_b, eos_b, ice_b, sf_b)
   end subroutine test_snowfall_disabled_bitident

end module test_ocean_ice_snowfall
