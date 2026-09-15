!! Analytic tests for the ocean-side frazil accumulator
!! (`ice_frazil_accumulate` in `rdb_ice_frazil` + the `frazil_heat` /
!! `heat_budget_frazil` arrays on `ocean_sea_ice_t`) — sea-ice PR 1.
!!
!! The kernel clamps the PHYSICAL surface layer (k = nz, bottom-up
!! convention) at the freezing point `T_f = eos_freezing_point(S, 0)`
!! and banks the removed supercooling deficit
!!   deficit = ρ·Cp·h·(T_f − T)⁺   [J/m²]
!! into `frazil_heat`, mirroring the tracer increment `h·(T_f − T)`
!! (K·m) into `heat_budget_frazil`.
!!
!! Cases:
!!   * Supercooled cell — surface T clamped EXACTLY to T_f; deficit
!!     matches ρ·Cp·h·(T_f − T₀) to round-off; sub-surface layers and
!!     salinity untouched.
!!   * Warm cell — above-freezing water is untouched and banks zero.
!!   * Energy closure — Σ ρ·Cp·ΔhTr_T == Σ frazil_heat over the domain
!!     (the sensible heat ADDED to the ocean equals the latent bank).
!!   * Ghost / dry cells — ghosts (physical-cells-only contract) and
!!     `wet_mask = 0` land columns are untouched.
!!   * Idempotence — a second call on the clamped state is a no-op.
module test_ocean_frazil
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, RHO_WATER
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_eos, only: eos_t, eos_freezing_point
   use rdb_ocean_surface_flux, only: SEAWATER_CP
   use rdb_ice_state, only: ocean_sea_ice_t
   use rdb_ice_frazil, only: ice_frazil_accumulate
   implicit none
   private

   public :: collect_ocean_frazil_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 3
   real(wp), parameter :: H_LAYER = 10.0_wp
      !! Uniform layer thickness (m).
   real(wp), parameter :: S_INIT = 35.0_wp
      !! Uniform salinity (PSU) — T_f = −1.89 °C.
   real(wp), parameter :: T_COLD = -3.0_wp
      !! Supercooled surface IC (°C), well below T_f.
   real(wp), parameter :: T_WARM = 2.0_wp
      !! Above-freezing surface IC (°C).
   real(wp), parameter :: T_DEEP = 1.0_wp
      !! Sub-surface IC (°C) — below-freezing checks must not touch it.

contains

   subroutine collect_ocean_frazil_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("frazil_clamps_to_freezing_point", test_clamp), &
                  new_unittest("frazil_banks_analytic_deficit", test_bank), &
                  new_unittest("frazil_warm_cell_untouched", test_warm_noop), &
                  new_unittest("frazil_energy_conservation", test_energy), &
                  new_unittest("frazil_skips_ghosts_and_dry", test_ghost_dry), &
                  new_unittest("frazil_second_call_idempotent", test_idempotent) &
                  ]
   end subroutine collect_ocean_frazil_tests

   subroutine setup_state(grid, ms, eos, ice, t_surface)
      !! Tiny ocean state: uniform h + S, sub-surface layers at T_DEEP,
      !! the surface layer at `t_surface`.  The ice slot is enabled the
      !! way the driver path does it (enable latched before init).
      type(hgrid_t), intent(out) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(eos_t), intent(out) :: eos
      type(ocean_sea_ice_t), intent(inout) :: ice
      real(wp), intent(in) :: t_surface

      call grid%init(6, 4, NGHOST, 1.0_wp, 1.0_wp)
      ms%nz_ml = NZ
      call ms%init(grid)
      call eos%init(grid)
      ice%enable = .true.
      call ice%init(grid)

      ms%h_layer = H_LAYER
      ms%tracers(ms%idx_salinity)%hTr = S_INIT*H_LAYER
      ms%tracers(ms%idx_temperature)%hTr = T_DEEP*H_LAYER
      ms%tracers(ms%idx_temperature)%hTr(:, :, NZ) = t_surface*H_LAYER
   end subroutine setup_state

   subroutine run_frazil(grid, eos, ms, ice)
      !! Device-wrapped kernel invocation (mirrors `run_apply` in
      !! test_ocean_surface_flux): map state, run, pull hTr + the ice
      !! arrays host-ward, unmap.
      type(hgrid_t), intent(in) :: grid
      type(eos_t), intent(in) :: eos
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_sea_ice_t), intent(inout) :: ice

      !$acc enter data copyin(ms)
      call ms%enter_data()
      call ice%enter_data()
      call ice_frazil_accumulate(grid, eos, ms, ice%frazil_heat, &
                                 ice%heat_budget_frazil)
      associate (hT => ms%tracers(ms%idx_temperature)%hTr, &
                 fz => ice%frazil_heat, fb => ice%heat_budget_frazil)
         !$acc update self(hT, fz, fb)
      end associate
      call ice%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms)
   end subroutine run_frazil

   subroutine teardown(ms, eos, ice)
      type(multilayer_state_t), intent(inout) :: ms
      type(eos_t), intent(inout) :: eos
      type(ocean_sea_ice_t), intent(inout) :: ice
      call ice%destroy()
      call eos%destroy()
      call ms%destroy()
   end subroutine teardown

   ! -----------------------------------------------------------------
   ! Cases
   ! -----------------------------------------------------------------

   subroutine test_clamp(error)
      !! (a) A supercooled surface cell lands EXACTLY at T_f; salinity
      !! and the sub-surface layers are untouched.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      real(wp) :: t_f, t_after, drift
      integer :: ip, jp, k
      checks: block

         call setup_state(grid, ms, eos, ice, T_COLD)
         t_f = eos_freezing_point(eos, S_INIT, 0.0_wp)
         call run_frazil(grid, eos, ms, ice)

         ip = grid%nx_total/2
         jp = grid%ny_total/2
         t_after = ms%tracers(ms%idx_temperature)%hTr(ip, jp, NZ)/H_LAYER
         call check(error, abs(t_after - t_f) < 1.0e-12_wp, &
                    "surface T must be clamped exactly to T_f")
         if (allocated(error)) exit checks

         drift = 0.0_wp
         do k = 1, NZ - 1
            drift = max(drift, abs(ms%tracers(ms%idx_temperature)%hTr(ip, jp, k) &
                                   - T_DEEP*H_LAYER))
         end do
         call check(error, drift < 1.0e-14_wp, &
                    "sub-surface layers must be untouched")
         if (allocated(error)) exit checks
         call check(error, abs(ms%tracers(ms%idx_salinity)%hTr(ip, jp, NZ) &
                               - S_INIT*H_LAYER) < 1.0e-14_wp, &
                    "salinity must be untouched")

      end block checks
      call teardown(ms, eos, ice)
   end subroutine test_clamp

   subroutine test_bank(error)
      !! (b) The banked deficit matches ρ·Cp·h·(T_f − T₀) to round-off,
      !! and the budget accumulator mirrors h·(T_f − T₀).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      real(wp) :: t_f, expected, tol
      integer :: ip, jp
      checks: block

         call setup_state(grid, ms, eos, ice, T_COLD)
         t_f = eos_freezing_point(eos, S_INIT, 0.0_wp)
         call run_frazil(grid, eos, ms, ice)

         ip = grid%nx_total/2
         jp = grid%ny_total/2
         expected = RHO_WATER*SEAWATER_CP*H_LAYER*(t_f - T_COLD)
         tol = 1.0e-14_wp*abs(expected)
         call check(error, abs(ice%frazil_heat(ip, jp) - expected) <= tol, &
                    "frazil_heat must equal rho*Cp*h*(T_f - T0) to round-off")
         if (allocated(error)) exit checks
         call check(error, abs(ice%heat_budget_frazil(ip, jp, 1) &
                               - H_LAYER*(t_f - T_COLD)) < 1.0e-12_wp, &
                    "budget accumulator must mirror h*(T_f - T0)")

      end block checks
      call teardown(ms, eos, ice)
   end subroutine test_bank

   subroutine test_warm_noop(error)
      !! (c) A cell already above freezing is untouched and banks zero.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      integer :: ip, jp
      checks: block

         call setup_state(grid, ms, eos, ice, T_WARM)
         call run_frazil(grid, eos, ms, ice)

         ip = grid%nx_total/2
         jp = grid%ny_total/2
         call check(error, abs(ms%tracers(ms%idx_temperature)%hTr(ip, jp, NZ) &
                               - T_WARM*H_LAYER) < 1.0e-14_wp, &
                    "above-freezing surface T must be untouched")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(ice%frazil_heat)) < 1.0e-14_wp, &
                    "warm ocean must bank zero frazil heat")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(ice%heat_budget_frazil)) < 1.0e-14_wp, &
                    "warm ocean must leave the budget accumulator at zero")

      end block checks
      call teardown(ms, eos, ice)
   end subroutine test_warm_noop

   subroutine test_energy(error)
      !! (d) Energy closure: the sensible heat GAINED by the ocean
      !! (ρ·Cp·Σ ΔhTr_T) equals the frazil heat banked (Σ frazil_heat),
      !! cell-by-cell summed over the whole domain — banked, not
      !! discarded.  Uses a mixed IC (west supercooled, east warm) so
      !! the closure is not trivially uniform.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      real(wp), allocatable :: hT_ic(:, :, :)
      real(wp) :: ocean_gain, banked, tol
      integer :: i, j, k
      checks: block

         call setup_state(grid, ms, eos, ice, T_COLD)
         ! East half warm — those columns must not participate.
         do j = 1, grid%ny_total
            do i = grid%nx_total/2 + 1, grid%nx_total
               ms%tracers(ms%idx_temperature)%hTr(i, j, NZ) = T_WARM*H_LAYER
            end do
         end do
         allocate (hT_ic, source=ms%tracers(ms%idx_temperature)%hTr)

         call run_frazil(grid, eos, ms, ice)

         ! Ocean sensible-heat gain over ALL layers (uniform cell area
         ! ⇒ the area factor cancels between the two sums).
         ocean_gain = 0.0_wp
         do k = 1, NZ
            do j = 1, grid%ny_total
               do i = 1, grid%nx_total
                  ocean_gain = ocean_gain + RHO_WATER*SEAWATER_CP* &
                               (ms%tracers(ms%idx_temperature)%hTr(i, j, k) - hT_ic(i, j, k))
               end do
            end do
         end do
         banked = sum(ice%frazil_heat)

         call check(error, banked > 0.0_wp, &
                    "mixed IC must bank a strictly positive deficit")
         if (allocated(error)) exit checks
         tol = 1.0e-13_wp*banked
         call check(error, abs(ocean_gain - banked) <= tol, &
                    "ocean sensible-heat gain must equal frazil heat banked")

      end block checks
      if (allocated(hT_ic)) deallocate (hT_ic)
      call teardown(ms, eos, ice)
   end subroutine test_energy

   subroutine test_ghost_dry(error)
      !! Ghost rows (physical-cells-only contract) and land columns
      !! (`wet_mask = 0`) stay supercooled and bank nothing.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      integer :: i_dry, j_dry
      checks: block

         call setup_state(grid, ms, eos, ice, T_COLD)
         ! One interior land column.
         i_dry = NGHOST + 2
         j_dry = NGHOST + 2
         ms%wet_mask(i_dry, j_dry) = 0.0_wp

         call run_frazil(grid, eos, ms, ice)

         call check(error, abs(ms%tracers(ms%idx_temperature)%hTr(1, 1, NZ) &
                               - T_COLD*H_LAYER) < 1.0e-14_wp, &
                    "ghost-corner surface T must be untouched")
         if (allocated(error)) exit checks
         call check(error, abs(ice%frazil_heat(1, 1)) < 1.0e-14_wp, &
                    "ghost cells must bank nothing")
         if (allocated(error)) exit checks
         call check(error, abs(ms%tracers(ms%idx_temperature)%hTr(i_dry, j_dry, NZ) &
                               - T_COLD*H_LAYER) < 1.0e-14_wp, &
                    "land column surface T must be untouched")
         if (allocated(error)) exit checks
         call check(error, abs(ice%frazil_heat(i_dry, j_dry)) < 1.0e-14_wp, &
                    "land column must bank nothing")

      end block checks
      call teardown(ms, eos, ice)
   end subroutine test_ghost_dry

   subroutine test_idempotent(error)
      !! A second pass over the already-clamped state must change
      !! nothing (T == T_f is NOT < T_f) — the bank only grows while
      !! genuine supercooling keeps being produced.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      real(wp), allocatable :: bank_1(:, :), hT_1(:, :, :)
      checks: block

         call setup_state(grid, ms, eos, ice, T_COLD)
         call run_frazil(grid, eos, ms, ice)
         allocate (bank_1, source=ice%frazil_heat)
         allocate (hT_1, source=ms%tracers(ms%idx_temperature)%hTr)

         call run_frazil(grid, eos, ms, ice)

         call check(error, maxval(abs(ice%frazil_heat - bank_1)) < 1.0e-14_wp, &
                    "second call must not grow the bank")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(ms%tracers(ms%idx_temperature)%hTr - hT_1)) &
                    < 1.0e-14_wp, "second call must not move the tracer")

      end block checks
      if (allocated(bank_1)) deallocate (bank_1)
      if (allocated(hT_1)) deallocate (hT_1)
      call teardown(ms, eos, ice)
   end subroutine test_idempotent

end module test_ocean_frazil
