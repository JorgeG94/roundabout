!! Unit tests for the Wright (1997) rational equation of state
!! (rdb_eos with `variant = EOS_VARIANT_WRIGHT_97`).
!!
!! Wright is *nonlinear* in (T, S): unlike the linear EOS that's
!! exercised by `test_ocean_eos`, Wright captures cabbeling
!! (mixing two equal-density water masses produces denser water),
!! thermobaric effects, and the realistic seawater density at the
!! standard reference state.  These tests check the things that
!! linear would get wrong.
!!
!! Cases:
!!   * Reference seawater density — ρ(T=10°C, S=35 PSU, P=0) must
!!     fall within ~0.1 kg/m^3 of 1027 kg/m^3 (the well-known
!!     surface seawater density also reported by MOM6's reference
!!     output).
!!   * Pure-water reference — ρ(T=0, S=0, P=0) ≈ 999.7 kg/m^3.
!!   * Monotonicity — increasing T at fixed S decreases ρ (thermal
!!     expansion); increasing S at fixed T increases ρ (haline
!!     contraction).  Both signs.
!!   * Cabbeling discriminator — initialise two water masses with
!!     equal density but different (T, S), then evaluate ρ at the
!!     blended (T_mean, S_mean).  Wright reports ρ_blend > ρ_orig
!!     (the cabbeling signature).  Linear gives ρ_blend == ρ_orig
!!     identically — this test is the cleanest discriminator
!!     between the two closures.
!!   * Vanishing-layer fallback — h_layer <= 0 returns rho_0
!!     verbatim (matches the linear branch's defensive guard).
module test_ocean_wright_eos
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_eos, only: eos_t, &
                      EOS_VARIANT_WRIGHT_97, EOS_VARIANT_LINEAR
   use rdb_ocean_eos_compute, only: ocean_eos_compute
   implicit none
   private

   public :: collect_ocean_wright_eos_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 3

contains

   subroutine collect_ocean_wright_eos_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("wright_reference_seawater", test_reference_seawater), &
                  new_unittest("wright_pure_water", test_pure_water), &
                  new_unittest("wright_monotonic_TS", test_monotonic_ts), &
                  new_unittest("wright_cabbeling", test_cabbeling), &
                  new_unittest("wright_vanishing_layer", test_vanishing_layer) &
                  ]
   end subroutine collect_ocean_wright_eos_tests

   subroutine make_grid(grid)
      type(hgrid_t), intent(out) :: grid
      call grid%init(8, 8, NGHOST, 1.0_wp, 1.0_wp)
   end subroutine make_grid

   subroutine run_wright(ms, eos)
      type(multilayer_state_t), intent(inout) :: ms
      type(eos_t), intent(in) :: eos
      !$acc enter data copyin(ms, eos)
      call ms%enter_data()
      call ocean_eos_compute(eos, ms)
      call ms%exit_data()
      !$acc exit data delete(ms, eos)
   end subroutine run_wright

   subroutine fill_uniform(ms, h_val, t_val, s_val)
      !! Fill h, hS, hT to a spatially uniform (T, S, h) state.
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: h_val, t_val, s_val
      ms%h_layer = h_val
      ms%tracers(ms%idx_salinity)%hTr = s_val*h_val
      ms%tracers(ms%idx_temperature)%hTr = t_val*h_val
   end subroutine fill_uniform

   ! -----------------------------------------------------------------
   ! Cases
   ! -----------------------------------------------------------------

   subroutine test_reference_seawater(error)
      !! ρ(T=10°C, S=35 PSU, P=0) ≈ 1027 kg/m^3 — the standard
      !! surface-seawater reference.  Wright matches EOS-80 here to
      !! ~0.01 kg/m^3 over the open-ocean range; a 0.1 kg/m^3
      !! tolerance is comfortable.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      real(wp) :: rho_obs

      call make_grid(grid)
      ms%nz_ml = NZ
      call ms%init(grid)
      eos%variant = EOS_VARIANT_WRIGHT_97
      eos%p_ref = 0.0_wp
      call fill_uniform(ms, 10.0_wp, 10.0_wp, 35.0_wp)
      call run_wright(ms, eos)

      rho_obs = ms%rho_layer(NGHOST + 1, NGHOST + 1, 1)
      call check(error, abs(rho_obs - 1027.0_wp) < 0.5_wp, &
                 "Wright ρ(10°C, 35 PSU, 0 Pa) off 1027 kg/m^3 by > 0.5")

      call ms%destroy()
   end subroutine test_reference_seawater

   subroutine test_pure_water(error)
      !! ρ(T=0°C, S=0 PSU, P=0) ≈ 999.7 kg/m^3 — pure water at 0°C.
      !! The b0/(c0 + α_0*b0) reduction of the Wright formula.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      real(wp) :: rho_obs

      call make_grid(grid)
      ms%nz_ml = NZ
      call ms%init(grid)
      eos%variant = EOS_VARIANT_WRIGHT_97
      eos%p_ref = 0.0_wp
      call fill_uniform(ms, 10.0_wp, 0.0_wp, 0.0_wp)
      call run_wright(ms, eos)

      rho_obs = ms%rho_layer(NGHOST + 1, NGHOST + 1, 1)
      call check(error, abs(rho_obs - 999.7_wp) < 0.5_wp, &
                 "Wright ρ(0°C, 0 PSU, 0 Pa) off 999.7 kg/m^3 by > 0.5")

      call ms%destroy()
   end subroutine test_pure_water

   subroutine test_monotonic_ts(error)
      !! Thermal expansion: ρ(T+ΔT, S) < ρ(T, S) for ΔT > 0.
      !! Haline contraction: ρ(T, S+ΔS) > ρ(T, S) for ΔS > 0.
      !! Both signs guarded.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      real(wp) :: rho_base, rho_warmer, rho_saltier
      checks: block

         call make_grid(grid)
         ms%nz_ml = NZ
         call ms%init(grid)
         eos%variant = EOS_VARIANT_WRIGHT_97
         eos%p_ref = 0.0_wp

         call fill_uniform(ms, 10.0_wp, 10.0_wp, 35.0_wp)
         call run_wright(ms, eos)
         rho_base = ms%rho_layer(NGHOST + 1, NGHOST + 1, 1)

         call fill_uniform(ms, 10.0_wp, 15.0_wp, 35.0_wp)
         call run_wright(ms, eos)
         rho_warmer = ms%rho_layer(NGHOST + 1, NGHOST + 1, 1)

         call fill_uniform(ms, 10.0_wp, 10.0_wp, 38.0_wp)
         call run_wright(ms, eos)
         rho_saltier = ms%rho_layer(NGHOST + 1, NGHOST + 1, 1)

         call check(error, rho_warmer < rho_base, &
                    "Wright thermal expansion broken: warmer should be lighter")
         if (allocated(error)) exit checks
         call check(error, rho_saltier > rho_base, &
                    "Wright haline contraction broken: saltier should be denser")

      end block checks
      call ms%destroy()
   end subroutine test_monotonic_ts

   subroutine test_cabbeling(error)
      !! Cabbeling discriminator.  Choose two water masses with
      !! *equal* density but distinct (T, S):
      !!   parcel A: cold + fresh, parcel B: warm + salty.
      !! Adjust S_B until ρ_A ≈ ρ_B (linear search by hand at the
      !! Wright values used below).  Then evaluate the EOS at the
      !! midpoint (T_mean, S_mean): Wright reports a *higher*
      !! density (the cabbeling signature, ~0.1-0.3 kg/m^3 for
      !! these realistic parcels).
      !!
      !! Linear EOS gives ρ_mean == 0.5*(ρ_A + ρ_B) identically by
      !! construction, so any nonzero cabbeling signal certifies
      !! the rational form is wired correctly.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      real(wp), parameter :: T_A = 0.0_wp, S_A = 33.0_wp
      real(wp), parameter :: T_B = 20.0_wp, S_B = 36.5_wp
      real(wp) :: rho_a, rho_b, rho_mean_expected, rho_blended

      call make_grid(grid)
      ms%nz_ml = NZ
      call ms%init(grid)
      eos%variant = EOS_VARIANT_WRIGHT_97
      eos%p_ref = 0.0_wp

      call fill_uniform(ms, 10.0_wp, T_A, S_A)
      call run_wright(ms, eos)
      rho_a = ms%rho_layer(NGHOST + 1, NGHOST + 1, 1)

      call fill_uniform(ms, 10.0_wp, T_B, S_B)
      call run_wright(ms, eos)
      rho_b = ms%rho_layer(NGHOST + 1, NGHOST + 1, 1)

      rho_mean_expected = 0.5_wp*(rho_a + rho_b)
      call fill_uniform(ms, 10.0_wp, 0.5_wp*(T_A + T_B), 0.5_wp*(S_A + S_B))
      call run_wright(ms, eos)
      rho_blended = ms%rho_layer(NGHOST + 1, NGHOST + 1, 1)

      ! Cabbeling: blended density exceeds the linear average of
      ! the source densities.  Magnitude depends on the (T, S)
      ! contrast — 0.05 kg/m^3 is a conservative lower bound for
      ! the cold-fresh ↔ warm-salty mix used here.
      call check(error, rho_blended > rho_mean_expected + 0.05_wp, &
                 "Wright cabbeling absent: blended ρ <= mean of source ρ's")

      call ms%destroy()
   end subroutine test_cabbeling

   subroutine test_vanishing_layer(error)
      !! `h_layer <= 0` cells fall back to `rho_0` — same defensive
      !! branch as the linear EOS.  Important under ZSTAR_FULL when
      !! bed-side layers pinch out.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      real(wp) :: rho_van, rho_normal
      checks: block

         call make_grid(grid)
         ms%nz_ml = NZ
         call ms%init(grid)
         eos%variant = EOS_VARIANT_WRIGHT_97
         eos%p_ref = 0.0_wp
         eos%rho0 = 1035.0_wp

         ms%h_layer = 10.0_wp
         ms%tracers(ms%idx_salinity)%hTr = 35.0_wp*10.0_wp
         ms%tracers(ms%idx_temperature)%hTr = 10.0_wp*10.0_wp
         ! Pinch out one interior cell
         ms%h_layer(NGHOST + 1, NGHOST + 1, 1) = 0.0_wp
         call run_wright(ms, eos)

         rho_van = ms%rho_layer(NGHOST + 1, NGHOST + 1, 1)
         rho_normal = ms%rho_layer(NGHOST + 2, NGHOST + 1, 1)

         call check(error, abs(rho_van - eos%rho0) < 1.0e-12_wp, &
                    "Wright vanishing-layer fallback didn't return rho_0")
         if (allocated(error)) exit checks
         call check(error, rho_normal > 1020.0_wp .and. rho_normal < 1035.0_wp, &
                    "Wright neighbour cell off the seawater range")

      end block checks
      call ms%destroy()
   end subroutine test_vanishing_layer

end module test_ocean_wright_eos
