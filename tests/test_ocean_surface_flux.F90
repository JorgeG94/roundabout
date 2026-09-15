!! Unit tests for the surface heat + salt flux apply kernel
!! (`ocean_surface_flux_apply_tracers` in
!! `rdb_ocean_surface_flux`).
!!
!! The kernel adds `dt · Q_heat(i,j) / (ρ₀·cp)` to the top layer's
!! `hT_temperature` and `dt · Q_salt(i,j) / ρ₀` to the top layer's
!! `hT_salinity`.  Bottom-up convention: top layer = `k = nz`.
!!
!! Cases:
!!   * Zero-flux no-op — Q_heat = Q_salt = 0, run N steps, every
!!     tracer field stays at its IC to round-off.
!!   * Constant heat flux warms the top layer — Q_heat > 0, after
!!     N steps the top-layer T has increased by N·dt·Q_heat /
!!     (ρ₀·cp·h_top) within tolerance.  Other layers untouched.
!!   * Constant salt flux salinifies the top layer — symmetric.
!!   * Non-uniform 2D Q_heat field: half the domain heated, half
!!     at zero; only the heated half warms.  Validates that the
!!     per-column field read is not broadcast from a scalar.
module test_ocean_surface_flux
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_eos, only: eos_t
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t, &
                                     ocean_surface_flux_apply_tracers, &
                                     SEAWATER_CP
   use rdb_ocean_budgets, only: ocean_budgets_t, BUDGET_HEAT_TOTAL, BUDGET_SALT_TOTAL
   implicit none
   private

   public :: collect_ocean_surface_flux_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 4

contains

   subroutine collect_ocean_surface_flux_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("surface_flux_zero_noop", test_zero_noop), &
                  new_unittest("surface_heat_flux_warms_top", test_heat_warms_top), &
                  new_unittest("surface_salt_flux_salinifies_top", test_salt_top), &
                  new_unittest("surface_flux_contributors_match_analytic", &
                               test_surface_contributors_analytic), &
                  new_unittest("nonuniform_q_heat_field_per_column", &
                               test_nonuniform_qheat) &
                  ]
   end subroutine collect_ocean_surface_flux_tests

   subroutine make_grid(grid, nx_phys, ny_phys)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      call grid%init(nx_phys, ny_phys, NGHOST, 1.0_wp, 1.0_wp)
   end subroutine make_grid

   subroutine run_apply(grid, sf, ms, dt)
      type(hgrid_t), intent(in) :: grid
      type(ocean_surface_flux_t), intent(inout) :: sf
         !! inout: enter_data / exit_data update device state
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: dt
      !$acc enter data copyin(ms, sf)
      call ms%enter_data()
      call sf%enter_data()
      call ocean_surface_flux_apply_tracers(grid, sf, ms, dt)
      associate (hT => ms%tracers(ms%idx_temperature)%hTr, &
                 hS => ms%tracers(ms%idx_salinity)%hTr)
         !$acc update self(hT, hS)
      end associate
      call sf%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, sf)
   end subroutine run_apply

   ! -----------------------------------------------------------------
   ! Cases
   ! -----------------------------------------------------------------

   subroutine test_zero_noop(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_surface_flux_t) :: sf
      type(eos_t) :: eos
      real(wp), parameter :: H_LAYER = 10.0_wp
      real(wp), parameter :: DT = 1.0_wp
      real(wp), allocatable :: hT_ic(:, :, :), hS_ic(:, :, :)
      checks: block

         call make_grid(grid, 6, 4)
         ms%nz_ml = NZ
         call ms%init(grid)
         call eos%init(grid)
         call sf%init(grid)

         ms%h_layer = H_LAYER
         ms%tracers(ms%idx_temperature)%hTr = eos%T_ref*H_LAYER
         ms%tracers(ms%idx_salinity)%hTr = eos%S_ref*H_LAYER
         allocate (hT_ic, source=ms%tracers(ms%idx_temperature)%hTr)
         allocate (hS_ic, source=ms%tracers(ms%idx_salinity)%hTr)

         call sf%set_surface_flux_const(0.0_wp, 0.0_wp)
         call run_apply(grid, sf, ms, DT)

         call check(error, maxval(abs(ms%tracers(ms%idx_temperature)%hTr - hT_ic)) < 1.0e-14_wp, &
                    "zero-flux apply modified hT")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(ms%tracers(ms%idx_salinity)%hTr - hS_ic)) < 1.0e-14_wp, &
                    "zero-flux apply modified hS")

      end block checks
      deallocate (hT_ic, hS_ic)
      call sf%destroy(); call eos%destroy(); call ms%destroy()
   end subroutine test_zero_noop

   subroutine test_heat_warms_top(error)
      !! Q_heat = 100 W/m² downward, no salt flux, h_top = 10 m,
      !! one apply step at dt = 3600 s.  Expected ΔT at the top:
      !!   ΔT = dt · Q / (ρ₀·cp·h_top) ≈ 9.06e-3 °C
      !! Convert to ΔhT = ΔT · h_top ≈ 9.06e-2 °C·m.  Other layers
      !! must stay at their IC.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_surface_flux_t) :: sf
      type(eos_t) :: eos
      real(wp), parameter :: H_LAYER = 10.0_wp
      real(wp), parameter :: DT = 3600.0_wp
      real(wp), parameter :: Q_HEAT = 100.0_wp
      real(wp) :: expected_dhT, max_top_err, max_other_drift
      integer :: i_probe, j_probe, k
      checks: block

         call make_grid(grid, 6, 4)
         ms%nz_ml = NZ
         call ms%init(grid)
         call eos%init(grid)
         call sf%init(grid)

         ms%h_layer = H_LAYER
         ms%tracers(ms%idx_temperature)%hTr = eos%T_ref*H_LAYER
         ms%tracers(ms%idx_salinity)%hTr = eos%S_ref*H_LAYER

         call sf%set_surface_flux_const(Q_HEAT, 0.0_wp)
         call run_apply(grid, sf, ms, DT)

         i_probe = grid%nx_total/2
         j_probe = grid%ny_total/2
         expected_dhT = DT*Q_HEAT/(sf%rho0*sf%cp)
         max_top_err = abs(ms%tracers(ms%idx_temperature)%hTr(i_probe, j_probe, NZ) &
                           - (eos%T_ref*H_LAYER + expected_dhT))
         max_other_drift = 0.0_wp
         do k = 1, NZ - 1
            max_other_drift = max(max_other_drift, &
                                  abs(ms%tracers(ms%idx_temperature)%hTr(i_probe, j_probe, k) - eos%T_ref*H_LAYER))
         end do

         call check(error, max_top_err < 1.0e-10_wp, &
                    "surface heat: top hT didn't match analytic")
         if (allocated(error)) exit checks
         call check(error, max_other_drift < 1.0e-14_wp, &
                    "surface heat: layers below the surface drifted")

      end block checks
      call sf%destroy(); call eos%destroy(); call ms%destroy()
   end subroutine test_heat_warms_top

   subroutine test_salt_top(error)
      !! Symmetric salt-flux test.  Q_salt > 0 → top hS rises by
      !! dt · Q_salt / ρ₀.  Other layers untouched.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_surface_flux_t) :: sf
      type(eos_t) :: eos
      real(wp), parameter :: H_LAYER = 10.0_wp
      real(wp), parameter :: DT = 3600.0_wp
      real(wp), parameter :: Q_SALT = 1.0e-5_wp   ! kg salt / m² / s
      real(wp) :: expected_dhS, max_top_err, max_other_drift
      integer :: i_probe, j_probe, k
      checks: block

         call make_grid(grid, 6, 4)
         ms%nz_ml = NZ
         call ms%init(grid)
         call eos%init(grid)
         call sf%init(grid)

         ms%h_layer = H_LAYER
         ms%tracers(ms%idx_temperature)%hTr = eos%T_ref*H_LAYER
         ms%tracers(ms%idx_salinity)%hTr = eos%S_ref*H_LAYER

         call sf%set_surface_flux_const(0.0_wp, Q_SALT)
         call run_apply(grid, sf, ms, DT)

         i_probe = grid%nx_total/2
         j_probe = grid%ny_total/2
         expected_dhS = DT*Q_SALT/sf%rho0
         max_top_err = abs(ms%tracers(ms%idx_salinity)%hTr(i_probe, j_probe, NZ) &
                           - (eos%S_ref*H_LAYER + expected_dhS))
         max_other_drift = 0.0_wp
         do k = 1, NZ - 1
            max_other_drift = max(max_other_drift, &
                                  abs(ms%tracers(ms%idx_salinity)%hTr(i_probe, j_probe, k) - eos%S_ref*H_LAYER))
         end do

         call check(error, max_top_err < 1.0e-10_wp, &
                    "surface salt: top hS didn't match analytic")
         if (allocated(error)) exit checks
         call check(error, max_other_drift < 1.0e-14_wp, &
                    "surface salt: layers below the surface drifted")

      end block checks
      call sf%destroy(); call eos%destroy(); call ms%destroy()
   end subroutine test_salt_top

   subroutine test_surface_contributors_analytic(error)
      !! Phase D v2 first sign-defined source: register the surface
      !! heat + salt budget contributors against a local
      !! `ocean_budgets_t`, apply uniform Q over N steps, drain, and
      !! verify the contributor totals match the closed-form
      !!   ∫ Q · dt / (ρ₀·cp) · dA · N_steps
      !! for heat and the analogous expression for salt to FP.  This
      !! is the validation gate that proves the contributor RHS
      !! actually closes against an independently-computed LHS for a
      !! sign-defined source — closed-basin internal redistribution
      !! would always pass trivially with everything = 0.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_surface_flux_t) :: sf
      type(eos_t) :: eos
      type(ocean_budgets_t) :: budgets
      real(wp), parameter :: H_LAYER = 10.0_wp
      real(wp), parameter :: DT = 3600.0_wp
      real(wp), parameter :: Q_HEAT = 100.0_wp
      real(wp), parameter :: Q_SALT = 1.0e-5_wp
      integer, parameter :: N_STEPS = 3
      real(wp) :: area, expected_heat, expected_salt
      integer :: step, idx_heat, idx_salt
      checks: block

         call make_grid(grid, 6, 4)
         ms%nz_ml = NZ
         call ms%init(grid)
         call eos%init(grid)
         call sf%init(grid)
         call budgets%init(grid)
         call budgets%register_contributor("surface_heat", BUDGET_HEAT_TOTAL, &
                                           ms%heat_budget_surface, &
                                           device_resident=.true.)
         idx_heat = budgets%n_contributors
         call budgets%register_contributor("surface_salt", BUDGET_SALT_TOTAL, &
                                           ms%salt_budget_surface, &
                                           device_resident=.true.)
         idx_salt = budgets%n_contributors

         ms%h_layer = H_LAYER
         ms%tracers(ms%idx_temperature)%hTr = eos%T_ref*H_LAYER
         ms%tracers(ms%idx_salinity)%hTr = eos%S_ref*H_LAYER

         call sf%set_surface_flux_const(Q_HEAT, Q_SALT)

         !$acc enter data copyin(ms, sf)
         call ms%enter_data()
         call sf%enter_data()
         do step = 1, N_STEPS
            call ocean_surface_flux_apply_tracers(grid, sf, ms, DT)
         end do
         ! Drain before exit_data — `acc update self if_present` needs
         ! the device buffer to still be mapped to pull the kernel's
         ! writes back.
         call budgets%drain_contributors()
         call sf%exit_data()
         call ms%exit_data()
         !$acc exit data delete(ms, sf)

         ! Surface flux fires at every interior cell on the surface
         ! layer.  `apply_tracers` writes at i=1..nx_total, j=1..ny_total
         ! (full grid including ghost), so the integral covers the full
         ! grid area.  Expected hTr·area integral after N_STEPS:
         area = real(grid%nx_total*grid%ny_total, wp)*grid%dx*grid%dy
         expected_heat = real(N_STEPS, wp)*DT*Q_HEAT/(sf%rho0*sf%cp)*area
         expected_salt = real(N_STEPS, wp)*DT*Q_SALT/sf%rho0*area

         call check(error, abs(budgets%contributors(idx_heat)%total_integrated &
                               - expected_heat) < 1.0e-10_wp*abs(expected_heat), &
                    "surface_heat contributor should match N · dt · Q/(ρ₀·cp) · area to FP")
         if (allocated(error)) exit checks
         call check(error, abs(budgets%contributors(idx_salt)%total_integrated &
                               - expected_salt) < 1.0e-10_wp*abs(expected_salt), &
                    "surface_salt contributor should match N · dt · Q/ρ₀ · area to FP")

      end block checks
      call budgets%destroy()
      call sf%destroy(); call eos%destroy(); call ms%destroy()
   end subroutine test_surface_contributors_analytic

   subroutine test_nonuniform_qheat(error)
      !! Non-uniform 2D Q_heat gate: fill only the west half of the domain
      !! with Q_heat = 100 W/m² (i <= nx_half) and the east half with zero.
      !! After one apply step only the west-half top layer must warm; the
      !! east half must stay at its IC.  This catches any path that ignores
      !! Q_heat(i,j) and broadcasts Q_heat_const uniformly.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_surface_flux_t) :: sf
      type(eos_t) :: eos
      real(wp), parameter :: H_LAYER = 10.0_wp
      real(wp), parameter :: DT = 3600.0_wp
      real(wp), parameter :: Q_HEATED = 100.0_wp
      integer :: i, j, nx, nx_half, i_east
      real(wp) :: hT_ic_east, hT_final_east, hT_ic_west, hT_final_west
      real(wp) :: expected_dhT
      checks: block

         call make_grid(grid, 6, 4)
         ms%nz_ml = NZ
         call ms%init(grid)
         call eos%init(grid)
         call sf%init(grid)

         ms%h_layer = H_LAYER
         ms%tracers(ms%idx_temperature)%hTr = eos%T_ref*H_LAYER
         ms%tracers(ms%idx_salinity)%hTr = eos%S_ref*H_LAYER

         ! Seed scalar fill (non-zero so apply_tracers doesn't early-exit)
         ! then overwrite the east half with zero to create the split.
         call sf%set_surface_flux_const(Q_HEATED, 0.0_wp)
         nx = grid%nx_total
         nx_half = nx/2
         do j = 1, grid%ny_total
            do i = nx_half + 1, nx
               sf%Q_heat(i, j) = 0.0_wp
            end do
         end do
         ! Q_salt stays zero; Q_heat_const != 0 keeps the guard alive.

         call run_apply(grid, sf, ms, DT)

         ! Probe representative west (heated) and east (unheated) cells.
         i_east = nx_half + 2
         expected_dhT = DT*Q_HEATED/(sf%rho0*sf%cp)

         ! West half (i = 2): hT should have risen by expected_dhT
         hT_ic_west = eos%T_ref*H_LAYER
         hT_final_west = ms%tracers(ms%idx_temperature)%hTr(2, grid%ny_total/2, NZ)
         call check(error, abs(hT_final_west - (hT_ic_west + expected_dhT)) < 1.0e-10_wp, &
                    "nonuniform: heated west cell top layer must warm by expected amount")
         if (allocated(error)) exit checks

         ! East half (i_east): hT must be unchanged (zero Q_heat there)
         hT_ic_east = eos%T_ref*H_LAYER
         hT_final_east = ms%tracers(ms%idx_temperature)%hTr(i_east, grid%ny_total/2, NZ)
         call check(error, abs(hT_final_east - hT_ic_east) < 1.0e-14_wp, &
                    "nonuniform: unheated east cell must not warm")

      end block checks
      call sf%destroy(); call eos%destroy(); call ms%destroy()
   end subroutine test_nonuniform_qheat

end module test_ocean_surface_flux
