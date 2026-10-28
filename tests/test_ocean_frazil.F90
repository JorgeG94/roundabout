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
   use rdb_constants, only: wp, RHO_WATER, GRAVITY
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_eos, only: eos_t, eos_freezing_point, eos_apply_ts_convention, &
                      TS_CONS_ABS, TFREEZE_FORM_TEOS_POLY
   use rdb_teos10_convert, only: ct_freezing_poly
   use rdb_ocean_surface_flux, only: SEAWATER_CP
   use rdb_ice_state, only: ocean_sea_ice_t
   use rdb_ice_frazil, only: ice_frazil_accumulate, ice_frazil_accumulate_column
   use rdb_ocean_surface_flux, only: bound_salinity_top_impl
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
                  new_unittest("frazil_second_call_idempotent", test_idempotent), &
                  new_unittest("frazil_pot_prac_bit_identical", test_pot_prac_bits), &
                  new_unittest("frazil_cons_abs_onset_and_cp", test_cons_abs_onset), &
                  new_unittest("frazil_column_pressure_dependence", test_col_pressure), &
                  new_unittest("frazil_column_carry_up", test_col_carry_up), &
                  new_unittest("frazil_column_no_absorber_hands_off", test_col_no_absorber), &
                  new_unittest("frazil_column_reclaim_branches", test_col_reclaim), &
                  new_unittest("frazil_bound_salinity", test_bound_salinity) &
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

   subroutine test_pot_prac_bits(error)
      !! (d) pot_prac default: the bank is EXACTLY rho*3992*(h*(T_f - T)) --
      !! the `eos%c_p` read is bit-identical to the old SEAWATER_CP constant.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      real(wp) :: t_f, t_sfc, expected
      integer :: ip, jp
      checks: block
         call setup_state(grid, ms, eos, ice, T_COLD)
         t_f = eos_freezing_point(eos, S_INIT, 0.0_wp)
         t_sfc = (T_COLD*H_LAYER)/H_LAYER
         expected = 0.0_wp + RHO_WATER*SEAWATER_CP*(H_LAYER*(t_f - t_sfc))
         call run_frazil(grid, eos, ms, ice)
         ip = grid%nx_total/2
         jp = grid%ny_total/2
         call check(error, ice%frazil_heat(ip, jp) == expected, &
                    "pot_prac bank must be bit-identical to rho*SEAWATER_CP*h*(Tf-T)")
      end block checks
      call teardown(ms, eos, ice)
   end subroutine test_pot_prac_bits

   subroutine test_cons_abs_onset(error)
      !! (a)+(c) cons_abs + teos_poly at SA = 35.16504, p = 0: the onset is
      !! the CT liquidus ct_freezing_poly = -1.9165 degC (below the linear
      !! -1.899, so a cell at -1.905 is NOT supercooled), a cell at -3 banks
      !! rho*c_p*h*(CT_f - CT) with c_p taken from the knob (4100, not 3992),
      !! and the cell is clamped to CT_f.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      real(wp), parameter :: SA = 35.16504_wp
      real(wp), parameter :: CP_TEST = 4100.0_wp
      real(wp) :: ct_f, t_lin, expected
      integer :: ip, jp
      checks: block
         ! --- no frazil between the two liquidi (-1.905 < -1.899 linear) ---
         call setup_state(grid, ms, eos, ice, -1.905_wp)
         t_lin = eos_freezing_point(eos, SA, 0.0_wp)
         call eos_apply_ts_convention(eos, TS_CONS_ABS, TFREEZE_FORM_TEOS_POLY, 0.0_wp)
         eos%c_p = CP_TEST
         ms%tracers(ms%idx_salinity)%hTr = SA*H_LAYER
         ct_f = eos_freezing_point(eos, SA, 0.0_wp)
         call check(error, abs(ct_f - ct_freezing_poly(SA, 0.0_wp, 0.0_wp)) < 1.0e-14_wp, &
                    "cons_abs onset must be ct_freezing_poly(SA, 0, 0)")
         if (allocated(error)) exit checks
         call check(error, abs(ct_f + 1.92_wp) < 0.01_wp .and. ct_f < t_lin - 0.01_wp, &
                    "CT_f ~ -1.9165 degC, shifted DOWN against the linear -1.899")
         if (allocated(error)) exit checks
         call run_frazil(grid, eos, ms, ice)
         ip = grid%nx_total/2
         jp = grid%ny_total/2
         call check(error, ice%frazil_heat(ip, jp) == 0.0_wp, &
                    "-1.905 is above CT_f: no frazil under cons_abs")
         if (allocated(error)) exit checks
         call teardown(ms, eos, ice)

         ! --- supercooled: bank uses c_p from the knob, clamp to CT_f ---
         call setup_state(grid, ms, eos, ice, T_COLD)
         call eos_apply_ts_convention(eos, TS_CONS_ABS, TFREEZE_FORM_TEOS_POLY, 0.0_wp)
         eos%c_p = CP_TEST
         ms%tracers(ms%idx_salinity)%hTr = SA*H_LAYER
         call run_frazil(grid, eos, ms, ice)
         expected = RHO_WATER*CP_TEST*H_LAYER*(ct_f - T_COLD)
         call check(error, abs(ice%frazil_heat(ip, jp) - expected) < 1.0e-9_wp*expected, &
                    "cons_abs bank must be rho*c_p(knob)*h*(CT_f - CT)")
         if (allocated(error)) exit checks
         call check(error, abs(ms%tracers(ms%idx_temperature)%hTr(ip, jp, NZ)/H_LAYER - ct_f) &
                    < 1.0e-12_wp, "surface CT clamped to CT_f")
      end block checks
      call teardown(ms, eos, ice)
   end subroutine test_cons_abs_onset

   subroutine run_column(grid, eos, ms, ice, pdep, reclaim)
      !! Device-wrapped `ice_frazil_accumulate_column` (same mapping as
      !! `run_frazil`); `ice%frazil_heat` set on the host beforehand rides
      !! the `copyin` of `ice%enter_data`.
      type(hgrid_t), intent(in) :: grid
      type(eos_t), intent(in) :: eos
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_sea_ice_t), intent(inout) :: ice
      logical, intent(in) :: pdep, reclaim

      !$acc enter data copyin(ms)
      call ms%enter_data()
      call ice%enter_data()
      call ice_frazil_accumulate_column(grid, eos, ms, ice%frazil_heat, &
                                        ice%heat_budget_frazil, pdep, reclaim)
      associate (hT => ms%tracers(ms%idx_temperature)%hTr, &
                 fz => ice%frazil_heat, fb => ice%heat_budget_frazil)
         !$acc update self(hT, fz, fb)
      end associate
      call ice%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms)
   end subroutine run_column

   subroutine test_col_pressure(error)
      !! Linear liquidus T_f = -0.054 S - 7.53e-8 p (the `seaice` set) falls
      !! with depth, so a uniform T between T_f(layer 3 centre) and
      !! T_f(layer 2 centre) is supercooled ONLY in the surface layer when
      !! pressure dependence is on (3 x 300 m layers, centres at 150, 450,
      !! 750 m).  With it off (p = 0, T_f = -1.89) all three freeze.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      real(wp), parameter :: HL = 300.0_wp
      real(wp) :: gr, tf3, tf2, t0, expected
      integer :: ip, jp
      checks: block
         call setup_state(grid, ms, eos, ice, T_DEEP)
         gr = GRAVITY*eos%rho0
         tf3 = eos_freezing_point(eos, S_INIT, 0.5_wp*gr*HL)
         tf2 = eos_freezing_point(eos, S_INIT, 0.5_wp*gr*HL + 0.5_wp*gr*(HL + HL))
         t0 = 0.5_wp*(tf3 + tf2)
         ms%h_layer = HL
         ms%tracers(ms%idx_salinity)%hTr = S_INIT*HL
         ms%tracers(ms%idx_temperature)%hTr = t0*HL
         ip = grid%nx_total/2
         jp = grid%ny_total/2
         call run_column(grid, eos, ms, ice, .true., .true.)
         call check(error, abs(ms%tracers(ms%idx_temperature)%hTr(ip, jp, NZ)/HL - tf3) &
                    < 1.0e-12_wp, "pressure-dependent: surface layer must reach T_f(p_3)")
         if (allocated(error)) exit checks
         call check(error, ms%tracers(ms%idx_temperature)%hTr(ip, jp, 2) == t0*HL .and. &
                    ms%tracers(ms%idx_temperature)%hTr(ip, jp, 1) == t0*HL, &
                    "pressure-dependent: deeper layers (T_f lower than T) must NOT freeze")
         if (allocated(error)) exit checks
         expected = RHO_WATER*eos%c_p*HL*(tf3 - t0)
         call check(error, abs(ice%frazil_heat(ip, jp) - expected) <= 1.0e-12_wp*expected, &
                    "pressure-dependent: F = rho c_p h (T_f(p_3) - T0), surface layer only")
         if (allocated(error)) exit checks
         call teardown(ms, eos, ice)

         ! Control: pressure dependence OFF => p = 0, all three layers freeze.
         call setup_state(grid, ms, eos, ice, T_DEEP)
         ms%h_layer = HL
         ms%tracers(ms%idx_salinity)%hTr = S_INIT*HL
         ms%tracers(ms%idx_temperature)%hTr = t0*HL
         call run_column(grid, eos, ms, ice, .false., .true.)
         expected = 3.0_wp*RHO_WATER*eos%c_p*HL*(eos_freezing_point(eos, S_INIT, 0.0_wp) - t0)
         call check(error, abs(ice%frazil_heat(ip, jp) - expected) <= 1.0e-12_wp*expected, &
                    "p = 0 control: all three layers freeze at T_f(S, 0)")
      end block checks
      call teardown(ms, eos, ice)
   end subroutine test_col_pressure

   subroutine test_col_carry_up(error)
      !! Bed layer at -4 degC (supercooled), a warm layer (+1) above it,
      !! a warmer surface (+2).  The bed freezes to T_f and carries
      !! c = hc (T_f + 4) upward; the warm layer absorbs it and cools by
      !! exactly c/hc.  F stays 0, the surface is untouched, and the column
      !! heat sum(hTr) is unchanged.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      real(wp) :: tf, c_over_hc, tot0, tot1
      integer :: ip, jp, k
      checks: block
         call setup_state(grid, ms, eos, ice, 2.0_wp)
         ms%tracers(ms%idx_temperature)%hTr(:, :, 1) = -4.0_wp*H_LAYER
         ms%tracers(ms%idx_temperature)%hTr(:, :, 2) = 1.0_wp*H_LAYER
         tf = eos_freezing_point(eos, S_INIT, 0.0_wp)
         c_over_hc = tf + 4.0_wp
         ip = grid%nx_total/2
         jp = grid%ny_total/2
         tot0 = 0.0_wp
         do k = 1, NZ
            tot0 = tot0 + ms%tracers(ms%idx_temperature)%hTr(ip, jp, k)
         end do
         call run_column(grid, eos, ms, ice, .false., .true.)
         tot1 = 0.0_wp
         do k = 1, NZ
            tot1 = tot1 + ms%tracers(ms%idx_temperature)%hTr(ip, jp, k)
         end do
         call check(error, abs(ms%tracers(ms%idx_temperature)%hTr(ip, jp, 1)/H_LAYER - tf) &
                    < 1.0e-14_wp, "bed layer clamped to T_f")
         if (allocated(error)) exit checks
         call check(error, abs(ms%tracers(ms%idx_temperature)%hTr(ip, jp, 2)/H_LAYER &
                               - (1.0_wp - c_over_hc)) < 1.0e-14_wp, &
                    "warm layer cooled by exactly c/hc")
         if (allocated(error)) exit checks
         call check(error, ms%tracers(ms%idx_temperature)%hTr(ip, jp, NZ) == 2.0_wp*H_LAYER, &
                    "surface layer above the absorber is untouched")
         if (allocated(error)) exit checks
         call check(error, ice%frazil_heat(ip, jp) == 0.0_wp, "nothing reaches the ice: F = 0")
         if (allocated(error)) exit checks
         call check(error, abs(tot1 - tot0) <= 1.0e-14_wp*abs(tot0) + 1.0e-14_wp, &
                    "column heat sum(hTr) unchanged by the carry")
         if (allocated(error)) exit checks
         call check(error, abs(ice%heat_budget_frazil(ip, jp, 1)) < 1.0e-12_wp, &
                    "budget accumulator sees no net heat")
      end block checks
      call teardown(ms, eos, ice)
   end subroutine test_col_carry_up

   subroutine test_col_no_absorber(error)
      !! Same bed deficit but every layer above is already at T_f: the
      !! carry cannot be absorbed (c + hc*0 > 0) and reaches the top.
      !! F = hc (T_f + 4) and the budget accumulator equals h (T_f + 4).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      real(wp) :: tf, expected
      integer :: ip, jp
      checks: block
         call setup_state(grid, ms, eos, ice, T_DEEP)
         tf = eos_freezing_point(eos, S_INIT, 0.0_wp)
         ms%tracers(ms%idx_temperature)%hTr(:, :, 1) = -4.0_wp*H_LAYER
         ms%tracers(ms%idx_temperature)%hTr(:, :, 2) = tf*H_LAYER
         ms%tracers(ms%idx_temperature)%hTr(:, :, NZ) = tf*H_LAYER
         ip = grid%nx_total/2
         jp = grid%ny_total/2
         call run_column(grid, eos, ms, ice, .false., .true.)
         expected = RHO_WATER*eos%c_p*H_LAYER*(tf + 4.0_wp)
         call check(error, abs(ice%frazil_heat(ip, jp) - expected) <= 1.0e-12_wp*expected, &
                    "unabsorbed carry is handed to the ice: F = hc (T_f + 4)")
         if (allocated(error)) exit checks
         call check(error, abs(ice%heat_budget_frazil(ip, jp, 1) - H_LAYER*(tf + 4.0_wp)) &
                    < 1.0e-12_wp, "budget accumulator = h (T_f + 4)")
      end block checks
      call teardown(ms, eos, ice)
   end subroutine test_col_no_absorber

   subroutine test_col_reclaim(error)
      !! Standing bank F on a warm surface (T = 2): (a) F < hc (T - T_f)
      !! cools the surface by F/hc and zeroes F; (b) F > hc (T - T_f)
      !! clamps the surface to T_f and leaves F - hc (T - T_f); (c) with
      !! reclaim off the bank is untouched.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      real(wp) :: tf, hc, room
      integer :: ip, jp
      checks: block
         call setup_state(grid, ms, eos, ice, 2.0_wp)
         tf = eos_freezing_point(eos, S_INIT, 0.0_wp)
         hc = RHO_WATER*eos%c_p*H_LAYER
         room = hc*(2.0_wp - tf)
         ip = grid%nx_total/2
         jp = grid%ny_total/2
         ice%frazil_heat = 0.25_wp*room
         call run_column(grid, eos, ms, ice, .false., .true.)
         call check(error, ice%frazil_heat(ip, jp) == 0.0_wp, "(a) bank fully reclaimed")
         if (allocated(error)) exit checks
         call check(error, abs(ms%tracers(ms%idx_temperature)%hTr(ip, jp, NZ)/H_LAYER &
                               - (2.0_wp - 0.25_wp*room/hc)) < 1.0e-14_wp, &
                    "(a) surface cooled by exactly F/hc")
         if (allocated(error)) exit checks
         call check(error, abs(ice%heat_budget_frazil(ip, jp, 1) &
                               + 0.25_wp*room/(RHO_WATER*eos%c_p)) < 1.0e-12_wp, &
                    "(a) budget accumulator = -F/(rho c_p)")
         if (allocated(error)) exit checks
         call teardown(ms, eos, ice)

         call setup_state(grid, ms, eos, ice, 2.0_wp)
         ice%frazil_heat = 3.0_wp*room
         call run_column(grid, eos, ms, ice, .false., .true.)
         call check(error, abs(ice%frazil_heat(ip, jp) - 2.0_wp*room) <= 1.0e-12_wp*room, &
                    "(b) bank left with F - hc (T - T_f)")
         if (allocated(error)) exit checks
         call check(error, abs(ms%tracers(ms%idx_temperature)%hTr(ip, jp, NZ)/H_LAYER - tf) &
                    < 1.0e-14_wp, "(b) surface clamped to T_f")
         if (allocated(error)) exit checks
         call teardown(ms, eos, ice)

         call setup_state(grid, ms, eos, ice, 2.0_wp)
         ice%frazil_heat = 3.0_wp*room
         call run_column(grid, eos, ms, ice, .false., .false.)
         call check(error, ice%frazil_heat(ip, jp) == 3.0_wp*room .and. &
                    ms%tracers(ms%idx_temperature)%hTr(ip, jp, NZ) == 2.0_wp*H_LAYER, &
                    "(c) reclaim off: bank and surface untouched")
      end block checks
      call teardown(ms, eos, ice)
   end subroutine test_col_reclaim

   subroutine test_bound_salinity(error)
      !! A negative top-layer S is raised to `min_salinity`; positive S and
      !! every sub-surface layer are untouched.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 3, NY = 2, NK = 2
      real(wp) :: hTr(NX, NY, NK), h(NX, NY, NK), wet(NX, NY)
      integer :: kt(NX, NY)
      real(wp), parameter :: SMIN = 0.5_wp
      h = 4.0_wp
      wet = 1.0_wp
      kt = NK
      hTr = 30.0_wp*4.0_wp
      hTr(1, 1, NK) = -2.0_wp*4.0_wp
      hTr(2, 1, NK) = 35.0_wp*4.0_wp
      hTr(1, 2, 1) = -1.0_wp
      wet(3, 2) = 0.0_wp
      hTr(3, 2, NK) = -8.0_wp
      !$acc enter data copyin(hTr, h, wet, kt)
      call bound_salinity_top_impl(hTr, h, wet, SMIN, kt, NK, NX, NY)
      !$acc update self(hTr)
      !$acc exit data delete(hTr, h, wet, kt)
      call check(error, abs(hTr(1, 1, NK)/4.0_wp - SMIN) < 1.0e-14_wp, &
                 "negative S raised to min_salinity")
      if (allocated(error)) return
      call check(error, hTr(2, 1, NK) == 140.0_wp .and. hTr(1, 2, NK) == 120.0_wp, &
                 "positive S untouched")
      if (allocated(error)) return
      call check(error, hTr(1, 2, 1) == -1.0_wp, "sub-surface layer untouched")
      if (allocated(error)) return
      call check(error, hTr(3, 2, NK) == -8.0_wp, "dry column untouched")
   end subroutine test_bound_salinity

end module test_ocean_frazil
