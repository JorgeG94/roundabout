!! Unit tests for the geothermal bottom-heat-flux apply kernel
!! (`ocean_geothermal_apply_tracers` in `rdb_ocean_geothermal`).
!!
!! The kernel adds `dt · Q_geo / (ρ₀·cp)` to the lowest *massive*
!! layer's `hT_temperature` (bottom-up convention: bed layer = `k = 1`).
!! No salt analogue.
!!
!! Cases:
!!   * T1 energy/conservation — constant Q_geo, one column, one step;
!!     the column heat integral `Σ_k(ρ₀·cp·hTr_T)·area` gains exactly
!!     `Q_geo·dt·area` to round-off.
!!   * T2 bed-layer-only — multi-layer column; only k=1 warms,
!!     k=2..nz hTr_T bit-unchanged.
!!   * T3 bit-identity (default off) — enable=.false. + q_geo=0 ⇒
!!     hTr_T and the budget array bit-identical to a no-geothermal run.
!!   * T4 vanishing bed layer — h_layer(:,:,1) < h_min, h_layer(:,:,2)
!!     massive ⇒ the increment lands in k=2, column integral closes.
module test_ocean_geothermal
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_eos, only: eos_t
   use rdb_ocean_geothermal, only: ocean_geothermal_t, &
                                   ocean_geothermal_apply_tracers
   implicit none
   private

   public :: collect_ocean_geothermal_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 4

contains

   subroutine collect_ocean_geothermal_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("geothermal_energy_conservation", test_energy_conservation), &
                  new_unittest("geothermal_bed_layer_only", test_bed_layer_only), &
                  new_unittest("geothermal_default_off_bit_identity", test_default_off), &
                  new_unittest("geothermal_vanishing_bed_layer", test_vanishing_bed) &
                  ]
   end subroutine collect_ocean_geothermal_tests

   subroutine make_grid(grid, nx_phys, ny_phys)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      call grid%init(nx_phys, ny_phys, NGHOST, 1.0_wp, 1.0_wp)
   end subroutine make_grid

   subroutine run_apply(grid, geo, ms, dt)
      !! Map the state + slot to the device, apply one geothermal step,
      !! pull hTr + budget back, tear down.  Mirrors the surface-flux
      !! test harness so the kernel is exercised on the GPU build.
      type(hgrid_t), intent(in) :: grid
      type(ocean_geothermal_t), intent(in) :: geo
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: dt
      !$acc enter data copyin(ms, geo)
      call ms%enter_data()
      call ocean_geothermal_apply_tracers(grid, geo, ms, dt)
      associate (hT => ms%tracers(ms%idx_temperature)%hTr, &
                 bgeo => ms%heat_budget_geothermal)
         !$acc update self(hT, bgeo)
      end associate
      call ms%exit_data()
      !$acc exit data delete(ms, geo)
   end subroutine run_apply

   ! -----------------------------------------------------------------
   ! T1 — energy / conservation
   ! -----------------------------------------------------------------
   subroutine test_energy_conservation(error)
      !! Q_geo = 80 W/m², no other forcing, all layers massive, one
      !! apply step at dt = 3600 s.  The column heat integral
      !!   E = Σ_k ρ₀·cp·hTr_T(i,j,k)
      !! must gain exactly Q_geo·dt per unit area (the kernel deposits
      !! ΔhT = Q_geo·dt/(ρ₀·cp) into k=1, so ρ₀·cp·ΔhT = Q_geo·dt).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_geothermal_t) :: geo
      type(eos_t) :: eos
      real(wp), parameter :: H_LAYER = 25.0_wp
      real(wp), parameter :: DT = 3600.0_wp
      real(wp), parameter :: Q_GEO = 80.0_wp
      real(wp) :: e_before, e_after, expected_gain
      integer :: ip, jp, k
      checks: block

         call make_grid(grid, 6, 4)
         ms%nz_ml = NZ
         call ms%init(grid)
         call eos%init(grid)
         call geo%init(grid)

         ms%h_layer = H_LAYER
         ms%tracers(ms%idx_temperature)%hTr = eos%T_ref*H_LAYER
         ms%tracers(ms%idx_salinity)%hTr = eos%S_ref*H_LAYER

         ip = grid%nx_total/2
         jp = grid%ny_total/2
         e_before = 0.0_wp
         do k = 1, NZ
            e_before = e_before + geo%rho0*geo%cp*ms%tracers(ms%idx_temperature)%hTr(ip, jp, k)
         end do

         geo%enable = .true.
         geo%q_geo_const = Q_GEO
         call run_apply(grid, geo, ms, DT)

         e_after = 0.0_wp
         do k = 1, NZ
            e_after = e_after + geo%rho0*geo%cp*ms%tracers(ms%idx_temperature)%hTr(ip, jp, k)
         end do
         ! Heat (J/m²) gained = Q_geo · dt.  Probe cell is interior (wet).
         expected_gain = Q_GEO*DT

         call check(error, abs((e_after - e_before) - expected_gain) < 1.0e-8_wp*abs(expected_gain), &
                    "geothermal: column heat integral didn't gain Q_geo·dt")

      end block checks
      call geo%destroy(); call eos%destroy(); call ms%destroy()
   end subroutine test_energy_conservation

   ! -----------------------------------------------------------------
   ! T2 — bed-layer-only
   ! -----------------------------------------------------------------
   subroutine test_bed_layer_only(error)
      !! Multi-layer column, all layers massive.  Only k=1 warms by the
      !! analytic ΔhT = Q_geo·dt/(ρ₀·cp); k=2..nz bit-unchanged.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_geothermal_t) :: geo
      type(eos_t) :: eos
      real(wp), parameter :: H_LAYER = 25.0_wp
      real(wp), parameter :: DT = 3600.0_wp
      real(wp), parameter :: Q_GEO = 80.0_wp
      real(wp) :: expected_dhT, max_bed_err, max_other_drift
      integer :: ip, jp, k
      checks: block

         call make_grid(grid, 6, 4)
         ms%nz_ml = NZ
         call ms%init(grid)
         call eos%init(grid)
         call geo%init(grid)

         ms%h_layer = H_LAYER
         ms%tracers(ms%idx_temperature)%hTr = eos%T_ref*H_LAYER
         ms%tracers(ms%idx_salinity)%hTr = eos%S_ref*H_LAYER

         geo%enable = .true.
         geo%q_geo_const = Q_GEO
         call run_apply(grid, geo, ms, DT)

         ip = grid%nx_total/2
         jp = grid%ny_total/2
         expected_dhT = DT*Q_GEO/(geo%rho0*geo%cp)
         ! Bed = k=1 (bottom-up convention).
         max_bed_err = abs(ms%tracers(ms%idx_temperature)%hTr(ip, jp, 1) &
                           - (eos%T_ref*H_LAYER + expected_dhT))
         max_other_drift = 0.0_wp
         do k = 2, NZ
            max_other_drift = max(max_other_drift, &
                                  abs(ms%tracers(ms%idx_temperature)%hTr(ip, jp, k) - eos%T_ref*H_LAYER))
         end do

         call check(error, max_bed_err < 1.0e-10_wp, &
                    "geothermal: bed hT didn't match analytic")
         if (allocated(error)) exit checks
         call check(error, max_other_drift < 1.0e-14_wp, &
                    "geothermal: layers above the bed drifted")

      end block checks
      call geo%destroy(); call eos%destroy(); call ms%destroy()
   end subroutine test_bed_layer_only

   ! -----------------------------------------------------------------
   ! T3 — bit-identity (default off)
   ! -----------------------------------------------------------------
   subroutine test_default_off(error)
      !! enable=.false. and q_geo=0 ⇒ no-op: every tracer field and the
      !! geothermal budget array stay bit-identical to their IC.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_geothermal_t) :: geo
      type(eos_t) :: eos
      real(wp), parameter :: H_LAYER = 25.0_wp
      real(wp), parameter :: DT = 3600.0_wp
      real(wp), allocatable :: hT_ic(:, :, :), budget_ic(:, :, :)
      checks: block

         call make_grid(grid, 6, 4)
         ms%nz_ml = NZ
         call ms%init(grid)
         call eos%init(grid)
         call geo%init(grid)

         ms%h_layer = H_LAYER
         ms%tracers(ms%idx_temperature)%hTr = eos%T_ref*H_LAYER
         ms%tracers(ms%idx_salinity)%hTr = eos%S_ref*H_LAYER
         allocate (hT_ic, source=ms%tracers(ms%idx_temperature)%hTr)
         allocate (budget_ic, source=ms%heat_budget_geothermal)

         geo%enable = .false.
         geo%q_geo_const = 0.0_wp
         call run_apply(grid, geo, ms, DT)

         call check(error, maxval(abs(ms%tracers(ms%idx_temperature)%hTr - hT_ic)) < 1.0e-14_wp, &
                    "default-off geothermal modified hT")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(ms%heat_budget_geothermal - budget_ic)) < 1.0e-14_wp, &
                    "default-off geothermal modified the budget array")

      end block checks
      if (allocated(hT_ic)) deallocate (hT_ic)
      if (allocated(budget_ic)) deallocate (budget_ic)
      call geo%destroy(); call eos%destroy(); call ms%destroy()
   end subroutine test_default_off

   ! -----------------------------------------------------------------
   ! T4 — vanishing bed layer
   ! -----------------------------------------------------------------
   subroutine test_vanishing_bed(error)
      !! h_layer(:,:,1) < h_min (pinched bed under ZSTAR_FULL),
      !! h_layer(:,:,2) massive.  The increment must land in k=2 (the
      !! lowest massive layer); k=1 unchanged; the column heat integral
      !! still gains Q_geo·dt per unit area.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_geothermal_t) :: geo
      type(eos_t) :: eos
      real(wp), parameter :: H_THICK = 25.0_wp
      real(wp), parameter :: H_THIN = 1.0e-6_wp   ! below h_min = 1e-3
      real(wp), parameter :: DT = 3600.0_wp
      real(wp), parameter :: Q_GEO = 80.0_wp
      real(wp) :: expected_dhT, k1_drift, k2_err, e_before, e_after
      integer :: ip, jp, k
      checks: block

         call make_grid(grid, 6, 4)
         ms%nz_ml = NZ
         call ms%init(grid)
         call eos%init(grid)
         call geo%init(grid)

         ! Bed layer (k=1) pinched; k=2..nz massive.
         ms%h_layer = H_THICK
         ms%h_layer(:, :, 1) = H_THIN
         ! hTr = T·h so the thin layer carries a tiny hTr (consistent IC).
         ms%tracers(ms%idx_temperature)%hTr = eos%T_ref*ms%h_layer
         ms%tracers(ms%idx_salinity)%hTr = eos%S_ref*ms%h_layer

         ip = grid%nx_total/2
         jp = grid%ny_total/2
         e_before = 0.0_wp
         do k = 1, NZ
            e_before = e_before + geo%rho0*geo%cp*ms%tracers(ms%idx_temperature)%hTr(ip, jp, k)
         end do

         geo%enable = .true.
         geo%q_geo_const = Q_GEO
         call run_apply(grid, geo, ms, DT)

         expected_dhT = DT*Q_GEO/(geo%rho0*geo%cp)
         ! k=1 (thin) must be untouched.
         k1_drift = abs(ms%tracers(ms%idx_temperature)%hTr(ip, jp, 1) - eos%T_ref*H_THIN)
         ! k=2 (lowest massive) absorbs the increment.
         k2_err = abs(ms%tracers(ms%idx_temperature)%hTr(ip, jp, 2) &
                      - (eos%T_ref*H_THICK + expected_dhT))
         e_after = 0.0_wp
         do k = 1, NZ
            e_after = e_after + geo%rho0*geo%cp*ms%tracers(ms%idx_temperature)%hTr(ip, jp, k)
         end do

         call check(error, k1_drift < 1.0e-14_wp, &
                    "geothermal vanishing-bed: pinched k=1 was modified")
         if (allocated(error)) exit checks
         call check(error, k2_err < 1.0e-10_wp, &
                    "geothermal vanishing-bed: increment didn't land in k=2")
         if (allocated(error)) exit checks
         call check(error, abs((e_after - e_before) - Q_GEO*DT) < 1.0e-8_wp*abs(Q_GEO*DT), &
                    "geothermal vanishing-bed: column heat integral didn't close")

      end block checks
      call geo%destroy(); call eos%destroy(); call ms%destroy()
   end subroutine test_vanishing_bed

end module test_ocean_geothermal
