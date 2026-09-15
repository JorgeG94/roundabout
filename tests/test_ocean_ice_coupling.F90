!! Analytic tests for the frazil->ice uptake + brine-rejection coupling
!! (sea-ice PR 3b): `ice_frazil_uptake` (`rdb_ice_frazil_uptake`) and
!! `ice_ocean_brine_flux` (`rdb_ice_ocean_coupler`).
!!
!! Harness mirrors `test_ocean_frazil.F90`: grid 6x4, NGHOST=2, NZ=3,
!! H_LAYER=10, S_INIT=35; plus an `ocean_surface_flux_t` initialised from
!! the grid.  `dt_therm = 1200 s`.
!!
!! Cases:
!!   * frazil_uptake_conserves — bank spend at t_f: mass identity, the
!!     analytic per-layer freeze mass, ENERGY CLOSURE (banked energy ==
!!     the enthalpy drop from seawater to new ice), salinity invariant,
!!     ocean tracers untouched by the uptake itself.
!!   * brine_rejection_salinifies — the SIGN gate: Q_salt > 0, run the
!!     REAL production apply kernel, assert S strictly increases and the
!!     salt increment matches the rejected salt exactly.
!!   * salt_budget_closes — full chain from genuine supercooling: the
!!     Boussinesq virtual-flux identity ocean-gain + ice-salt ==
!!     m_frozen*s_surf.
!!   * disabled_and_zero_bank_noop — an all-zero bank is an exact no-op
!!     through the whole chain.
!!   * uptake_skips_ghosts_dry_and_vanished — ghosts / land / vanished
!!     columns keep their banks and stay at zero ice mass.
module test_ocean_ice_coupling
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, H_VANISHED
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_eos, only: eos_t, eos_freezing_point
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t, ocean_surface_flux_apply_tracers
   use rdb_ice_state, only: ocean_sea_ice_t
   use rdb_ice_frazil, only: ice_frazil_accumulate
   use rdb_ice_frazil_uptake, only: ice_frazil_uptake, ICE_FRAZIL_T_OFFSET
   use rdb_ice_ocean_coupler, only: ice_ocean_brine_flux
   use rdb_ice_enthalpy, only: ICE_LAT_FUS, ICE_LIQ_LIM, ice_enthalpy_liquid, &
                               ice_enth_from_ts, ice_t_freeze
   use rdb_ice_column, only: ICE_BULK_SALINITY
   implicit none
   private

   public :: collect_ocean_ice_coupling_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 3
   real(wp), parameter :: H_LAYER = 10.0_wp
      !! Uniform layer thickness (m).
   real(wp), parameter :: S_INIT = 35.0_wp
      !! Uniform salinity (PSU).
   real(wp), parameter :: T_DEEP = 1.0_wp
      !! Sub-surface IC (°C) — untouched by every case.
   real(wp), parameter :: T_WARM = 2.0_wp
      !! Above-freezing surface IC (°C) — case 4 (zero bank, no-op).
   real(wp), parameter :: T_COLD = -3.0_wp
      !! Genuinely supercooled surface IC (°C) — case 3.
   real(wp), parameter :: DT_THERM = 1200.0_wp
      !! Thermo timestep (s).
   real(wp), parameter :: E_BANK = 3.34e6_wp
      !! Analytic frazil bank (J/m² of cell) seeded directly — case 1/2/5.

contains

   subroutine collect_ocean_ice_coupling_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("frazil_uptake_conserves", test_uptake_conserves), &
                  new_unittest("brine_rejection_salinifies", test_brine_sign), &
                  new_unittest("salt_budget_closes", test_salt_budget), &
                  new_unittest("disabled_and_zero_bank_noop", test_zero_bank_noop), &
                  new_unittest("uptake_skips_ghosts_dry_and_vanished", test_skip_gates) &
                  ]
   end subroutine collect_ocean_ice_coupling_tests

   ! -----------------------------------------------------------------
   ! Shared setup / run / teardown
   ! -----------------------------------------------------------------

   subroutine setup_state(grid, ms, eos, ice, sf, t_surface)
      !! Tiny ocean state: uniform h + S, sub-surface layers at T_DEEP,
      !! surface layer at `t_surface`.  Ice slot enabled the way the
      !! driver path does it (enable latched before init).  Surface-flux
      !! slot seeded to a zero background (tests override Q_salt via the
      !! coupler, not the constant).
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
      call ice%init(grid)
      call sf%init(grid)
      call sf%set_surface_flux_const(0.0_wp, 0.0_wp)

      ms%h_layer = H_LAYER
      ms%tracers(ms%idx_salinity)%hTr = S_INIT*H_LAYER
      ms%tracers(ms%idx_temperature)%hTr = T_DEEP*H_LAYER
      ms%tracers(ms%idx_temperature)%hTr(:, :, NZ) = t_surface*H_LAYER
   end subroutine setup_state

   subroutine run_uptake(grid, eos, ms, ice)
      !! GPU mem:separate discipline: map, run the uptake kernel only,
      !! pull the touched arrays host-ward, unmap.  Host-seeded values
      !! (bank, tracers) must be set BEFORE this call so the copyin
      !! carries them.
      type(hgrid_t), intent(in) :: grid
      type(eos_t), intent(in) :: eos
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_sea_ice_t), intent(inout) :: ice

      !$acc enter data copyin(ms)
      call ms%enter_data()
      call ice%enter_data()
      call ice_frazil_uptake(grid, eos, ms, ice, DT_THERM)
      associate (hT => ms%tracers(ms%idx_temperature)%hTr, &
                 hS => ms%tracers(ms%idx_salinity)%hTr, &
                 fz => ice%frazil_heat, mi => ice%m_ice, ei => ice%enth_ice, &
                 si => ice%sal_ice, mf => ice%m_frozen_diag, sd => ice%salt_flux_diag)
         !$acc update self(hT, hS, fz, mi, ei, si, mf, sd)
      end associate
      call ice%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms)
   end subroutine run_uptake

   subroutine run_chain(grid, eos, ms, ice, sf, do_accumulate)
      !! Full chain: (optional) frazil accumulate -> uptake -> brine
      !! coupler -> the REAL apply_tracers kernel.  One shared map/unmap
      !! spanning all four kernels (GPU mem:separate discipline — every
      !! array touched anywhere in the chain must be device-present for
      !! the whole span).
      type(hgrid_t), intent(in) :: grid
      type(eos_t), intent(in) :: eos
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_sea_ice_t), intent(inout) :: ice
      type(ocean_surface_flux_t), intent(inout) :: sf
      logical, intent(in) :: do_accumulate

      !$acc enter data copyin(ms, sf)
      call ms%enter_data()
      call ice%enter_data()
      call sf%enter_data()

      if (do_accumulate) then
         call ice_frazil_accumulate(grid, eos, ms, ice%frazil_heat, ice%heat_budget_frazil)
      end if
      call ice_frazil_uptake(grid, eos, ms, ice, DT_THERM)
      call ice_ocean_brine_flux(sf, ice)
      call ocean_surface_flux_apply_tracers(grid, sf, ms, DT_THERM)

      associate (hT => ms%tracers(ms%idx_temperature)%hTr, &
                 hS => ms%tracers(ms%idx_salinity)%hTr, &
                 fz => ice%frazil_heat, mi => ice%m_ice, ei => ice%enth_ice, &
                 si => ice%sal_ice, mf => ice%m_frozen_diag, sd => ice%salt_flux_diag, &
                 qs => sf%Q_salt)
         !$acc update self(hT, hS, fz, mi, ei, si, mf, sd, qs)
      end associate

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

   subroutine test_uptake_conserves(error)
      !! Case 1: seed surface T = t_f (post-clamp state) and a uniform
      !! bank E_BANK at all physical cells with empty ice.  After the
      !! uptake: bank zeroed, mass identity, analytic mass, ENERGY
      !! CLOSURE, salinity invariant, ocean tracers untouched.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      type(ocean_surface_flux_t) :: sf
      real(wp) :: t_f, enth_ocean, min_denth, t_frazil, enth_frazil, m_frozen_analytic
      real(wp) :: e_ice, tol
      real(wp), allocatable :: hT_ic(:, :, :), hS_ic(:, :, :)
      integer :: ip, jp, k
      checks: block

         call setup_state(grid, ms, eos, ice, sf, 0.0_wp)
         t_f = eos_freezing_point(eos, S_INIT, 0.0_wp)
         ! Re-seed the surface at t_f exactly (post-clamp state).
         ms%tracers(ms%idx_temperature)%hTr(:, :, NZ) = t_f*H_LAYER
         ice%frazil_heat = E_BANK

         allocate (hT_ic, source=ms%tracers(ms%idx_temperature)%hTr)
         allocate (hS_ic, source=ms%tracers(ms%idx_salinity)%hTr)

         call run_uptake(grid, eos, ms, ice)

         ip = grid%nx_total/2
         jp = grid%ny_total/2

         call check(error, ice%frazil_heat(ip, jp) == 0.0_wp, &
                    "bank must be zeroed exactly after the uptake")
         if (allocated(error)) exit checks

         call check(error, ice%m_ice(ip, jp, 1) > 0.0_wp, &
                    "category-1 ice mass must be strictly positive")
         if (allocated(error)) exit checks
         call check(error, abs(ice%m_ice(ip, jp, 1) - ice%m_frozen_diag(ip, jp)) &
                    <= 1.0e-14_wp*ice%m_ice(ip, jp, 1), &
                    "m_ice must equal m_frozen_diag to 1e-14 rel (empty IC)")
         if (allocated(error)) exit checks

         ! ---- Analytic mass (uniform bulk salinity: per-layer m_frazil identical) ----
         enth_ocean = ice_enthalpy_liquid(t_f, S_INIT)
         min_denth = ICE_LAT_FUS*(1.0_wp - ICE_LIQ_LIM)
         t_frazil = min(t_f, ice_t_freeze(ICE_BULK_SALINITY) - ICE_FRAZIL_T_OFFSET)
         enth_frazil = min(ice_enth_from_ts(t_frazil, ICE_BULK_SALINITY), enth_ocean - min_denth)
         m_frozen_analytic = E_BANK/(enth_ocean - enth_frazil)
         call check(error, abs(ice%m_frozen_diag(ip, jp) - m_frozen_analytic) &
                    <= 1.0e-12_wp*m_frozen_analytic, &
                    "m_frozen must match the analytic per-layer freeze mass to 1e-12 rel")
         if (allocated(error)) exit checks

         ! ---- ENERGY CLOSURE: banked energy == enthalpy drop to new ice ----
         e_ice = 0.0_wp
         do k = 1, ice%nk_ice
            e_ice = e_ice + (ice%m_ice(ip, jp, 1)/real(ice%nk_ice, wp))*ice%enth_ice(ip, jp, 1, k)
         end do
         tol = 1.0e-10_wp*E_BANK
         call check(error, abs(ice%m_frozen_diag(ip, jp)*enth_ocean - e_ice - E_BANK) <= tol, &
                    "banked energy must equal the seawater->ice enthalpy drop")
         if (allocated(error)) exit checks

         ! ---- Salinity invariant (started bulk, mixed with bulk) ----
         do k = 1, ice%nk_ice
            call check(error, abs(ice%sal_ice(ip, jp, 1, k) - ICE_BULK_SALINITY) < 1.0e-12_wp, &
                       "sal_ice must stay at ICE_BULK_SALINITY")
            if (allocated(error)) exit checks
         end do

         ! ---- Sanity: enthalpies finite/negative ----
         do k = 1, ice%nk_ice
            call check(error, ice%enth_ice(ip, jp, 1, k) < 0.0_wp, &
                       "new-ice enthalpy must be negative")
            if (allocated(error)) exit checks
         end do

         ! ---- Ocean tracers untouched by the uptake itself ----
         call check(error, maxval(abs(ms%tracers(ms%idx_temperature)%hTr - hT_ic)) < 1.0e-14_wp, &
                    "uptake must not touch hTr_T")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(ms%tracers(ms%idx_salinity)%hTr - hS_ic)) < 1.0e-14_wp, &
                    "uptake must not touch hTr_S")

      end block checks
      if (allocated(hT_ic)) deallocate (hT_ic)
      if (allocated(hS_ic)) deallocate (hS_ic)
      call teardown(ms, eos, ice, sf)
   end subroutine test_uptake_conserves

   subroutine test_brine_sign(error)
      !! Case 2 — the SIGN gate.  Continue from a case-1-style uptake
      !! (m_frozen > 0), then run the REAL `ice_ocean_brine_flux` +
      !! `ocean_surface_flux_apply_tracers` (no reimplementation).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      type(ocean_surface_flux_t) :: sf
      real(wp) :: t_f, s_before, s_after, q_salt_analytic, expected_dhS, tol
      integer :: ip, jp, k
      checks: block

         call setup_state(grid, ms, eos, ice, sf, 0.0_wp)
         t_f = eos_freezing_point(eos, S_INIT, 0.0_wp)
         ms%tracers(ms%idx_temperature)%hTr(:, :, NZ) = t_f*H_LAYER
         ice%frazil_heat = E_BANK

         ip = grid%nx_total/2
         jp = grid%ny_total/2
         s_before = ms%tracers(ms%idx_salinity)%hTr(ip, jp, NZ)/H_LAYER

         call run_chain(grid, eos, ms, ice, sf, do_accumulate=.false.)

         s_after = ms%tracers(ms%idx_salinity)%hTr(ip, jp, NZ)/H_LAYER

         call check(error, sf%has_salt, "has_salt must be latched true")
         if (allocated(error)) exit checks

         q_salt_analytic = ice%m_frozen_diag(ip, jp)*(S_INIT - ICE_BULK_SALINITY)/DT_THERM
         tol = 1.0e-12_wp*abs(q_salt_analytic)
         call check(error, abs(sf%Q_salt(ip, jp) - q_salt_analytic) <= tol, &
                    "Q_salt must equal m_frozen*(35-ICE_BULK_SALINITY)/dt_therm to 1e-12 rel")
         if (allocated(error)) exit checks
         call check(error, sf%Q_salt(ip, jp) > 0.0_wp, "Q_salt must be strictly positive")
         if (allocated(error)) exit checks

         call check(error, s_after > s_before, "S(nz) must strictly increase (the trap gate)")
         if (allocated(error)) exit checks

         expected_dhS = ice%m_frozen_diag(ip, jp)*(S_INIT - ICE_BULK_SALINITY)
         tol = 1.0e-12_wp*abs(expected_dhS)
         call check(error, abs(sf%rho0*(ms%tracers(ms%idx_salinity)%hTr(ip, jp, NZ) &
                                        - s_before*H_LAYER) - expected_dhS) <= tol, &
                    "salt added to the ocean must equal the salt rejected")
         if (allocated(error)) exit checks

         do k = 1, NZ - 1
            call check(error, abs(ms%tracers(ms%idx_salinity)%hTr(ip, jp, k) &
                                  - S_INIT*H_LAYER) < 1.0e-12_wp, &
                       "sub-surface salinity must be untouched")
            if (allocated(error)) exit checks
         end do

      end block checks
      call teardown(ms, eos, ice, sf)
   end subroutine test_brine_sign

   subroutine test_salt_budget(error)
      !! Case 3: full chain from genuine supercooling (T = -3.0).  The
      !! Boussinesq virtual-flux identity: total(ocean+ice) salt grows by
      !! exactly the salt content of the virtually-frozen seawater.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      type(ocean_surface_flux_t) :: sf
      real(wp), allocatable :: hS_ic(:, :, :)
      real(wp) :: s_surf_pre, t_f, d_ocean, ice_salt, resid, tol
      integer :: ip, jp, k
      checks: block

         call setup_state(grid, ms, eos, ice, sf, T_COLD)
         allocate (hS_ic, source=ms%tracers(ms%idx_salinity)%hTr)

         call run_chain(grid, eos, ms, ice, sf, do_accumulate=.true.)

         ip = grid%nx_total/2
         jp = grid%ny_total/2

         ! s_surf_pre: what the uptake saw.  The accumulate clamp raises T
         ! to t_f but never touches S, so s_surf_pre == the IC value.
         t_f = eos_freezing_point(eos, S_INIT, 0.0_wp)
         s_surf_pre = S_INIT
         call check(error, abs(hS_ic(ip, jp, NZ)/H_LAYER - s_surf_pre) < 1.0e-12_wp, &
                    "the clamp must not move S (s_surf_pre == IC value)")
         if (allocated(error)) exit checks

         call check(error, ice%m_frozen_diag(ip, jp) > 0.0_wp, &
                    "genuine supercooling must freeze a strictly positive mass")
         if (allocated(error)) exit checks

         d_ocean = 0.0_wp
         do k = 1, NZ
            d_ocean = d_ocean + sf%rho0*(ms%tracers(ms%idx_salinity)%hTr(ip, jp, k) &
                                         - hS_ic(ip, jp, k))
         end do

         ice_salt = 0.0_wp
         do k = 1, ice%nk_ice
            ice_salt = ice_salt + (ice%m_ice(ip, jp, 1)/real(ice%nk_ice, wp)) &
                       *ice%sal_ice(ip, jp, 1, k)
         end do

         resid = d_ocean + ice_salt - ice%m_frozen_diag(ip, jp)*s_surf_pre
         tol = 1.0e-10_wp*abs(ice%m_frozen_diag(ip, jp)*s_surf_pre)
         call check(error, abs(resid) <= tol, &
                    "d_ocean + ice_salt must equal m_frozen*s_surf_pre (virtual-flux closure)")

      end block checks
      if (allocated(hS_ic)) deallocate (hS_ic)
      call teardown(ms, eos, ice, sf)
   end subroutine test_salt_budget

   subroutine test_zero_bank_noop(error)
      !! Case 4: ice enabled, bank all-zero, warm surface (T=2).  The
      !! whole chain (with Q_salt_const=0) must be an exact no-op.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      type(ocean_surface_flux_t) :: sf
      real(wp), allocatable :: hT_ic(:, :, :), hS_ic(:, :, :)
      real(wp), allocatable :: m_ice_ic(:, :, :), enth_ice_ic(:, :, :, :)
      checks: block

         call setup_state(grid, ms, eos, ice, sf, T_WARM)
         ! Bank already all-zero from ice%init.
         allocate (hT_ic, source=ms%tracers(ms%idx_temperature)%hTr)
         allocate (hS_ic, source=ms%tracers(ms%idx_salinity)%hTr)
         allocate (m_ice_ic, source=ice%m_ice)
         allocate (enth_ice_ic, source=ice%enth_ice)

         call run_chain(grid, eos, ms, ice, sf, do_accumulate=.true.)

         call check(error, maxval(abs(ms%tracers(ms%idx_temperature)%hTr - hT_ic)) == 0.0_wp, &
                    "hT must be bit-identical")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(ms%tracers(ms%idx_salinity)%hTr - hS_ic)) == 0.0_wp, &
                    "hS must be bit-identical")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(ice%m_ice - m_ice_ic)) == 0.0_wp, &
                    "m_ice must be bit-identical")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(ice%enth_ice - enth_ice_ic)) == 0.0_wp, &
                    "enth_ice must be bit-identical")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(ice%m_frozen_diag)) == 0.0_wp, &
                    "m_frozen_diag must be all zero")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(ice%salt_flux_diag)) == 0.0_wp, &
                    "salt_flux_diag must be all zero")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(sf%Q_salt)) == 0.0_wp, &
                    "Q_salt must be all zero (Q_salt_const=0, no freezing)")

      end block checks
      if (allocated(hT_ic)) deallocate (hT_ic)
      if (allocated(hS_ic)) deallocate (hS_ic)
      if (allocated(m_ice_ic)) deallocate (m_ice_ic)
      if (allocated(enth_ice_ic)) deallocate (enth_ice_ic)
      call teardown(ms, eos, ice, sf)
   end subroutine test_zero_bank_noop

   subroutine test_skip_gates(error)
      !! Case 5: bank seeded everywhere INCLUDING ghosts; one land column
      !! and one vanished column.  Ghost/land/vanished banks must be
      !! untouched; wet cells must spend normally.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      type(ocean_surface_flux_t) :: sf
      real(wp) :: t_f
      integer :: i_land, j_land, i_van, j_van, ip, jp
      checks: block

         call setup_state(grid, ms, eos, ice, sf, 0.0_wp)
         t_f = eos_freezing_point(eos, S_INIT, 0.0_wp)
         ms%tracers(ms%idx_temperature)%hTr(:, :, NZ) = t_f*H_LAYER
         ice%frazil_heat = E_BANK   ! seeded everywhere, INCLUDING ghosts

         i_land = NGHOST + 2
         j_land = NGHOST + 2
         ms%wet_mask(i_land, j_land) = 0.0_wp

         i_van = NGHOST + 3
         j_van = NGHOST + 3
         ms%h_layer(i_van, j_van, NZ) = 0.5_wp*H_VANISHED

         call run_uptake(grid, eos, ms, ice)

         ! Ghost corner (1,1): untouched.
         call check(error, ice%frazil_heat(1, 1) == E_BANK, "ghost bank must be untouched")
         if (allocated(error)) exit checks
         call check(error, ice%m_ice(1, 1, 1) == 0.0_wp, "ghost m_ice must stay 0")
         if (allocated(error)) exit checks
         call check(error, ice%m_frozen_diag(1, 1) == 0.0_wp, "ghost diag must stay 0")
         if (allocated(error)) exit checks

         ! Land column: untouched.
         call check(error, ice%frazil_heat(i_land, j_land) == E_BANK, &
                    "land bank must be untouched")
         if (allocated(error)) exit checks
         call check(error, ice%m_ice(i_land, j_land, 1) == 0.0_wp, "land m_ice must stay 0")
         if (allocated(error)) exit checks
         call check(error, ice%m_frozen_diag(i_land, j_land) == 0.0_wp, &
                    "land diag must stay 0")
         if (allocated(error)) exit checks

         ! Vanished column: untouched.
         call check(error, ice%frazil_heat(i_van, j_van) == E_BANK, &
                    "vanished-column bank must be untouched")
         if (allocated(error)) exit checks
         call check(error, ice%m_ice(i_van, j_van, 1) == 0.0_wp, &
                    "vanished-column m_ice must stay 0")
         if (allocated(error)) exit checks
         call check(error, ice%m_frozen_diag(i_van, j_van) == 0.0_wp, &
                    "vanished-column diag must stay 0")
         if (allocated(error)) exit checks

         ! A normal wet interior cell must spend normally.
         ip = grid%nx_total/2
         jp = grid%ny_total/2
         call check(error, ice%frazil_heat(ip, jp) == 0.0_wp, &
                    "wet cell bank must be spent")
         if (allocated(error)) exit checks
         call check(error, ice%m_ice(ip, jp, 1) > 0.0_wp, &
                    "wet cell must freeze a strictly positive mass")

      end block checks
      call teardown(ms, eos, ice, sf)
   end subroutine test_skip_gates

end module test_ocean_ice_coupling
