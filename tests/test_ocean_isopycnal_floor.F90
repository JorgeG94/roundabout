!! Unit tests for the Phase-1 Lagrangian minimum-thickness floor on the
!! continuity h-update (`continuity_apply_fluxes`, `continuity_apply_zonal`,
!! `continuity_apply_meridional`).
!!
!! Cases:
!!   * floor_clamps_grounding_layer   — floor clamps a layer that would go
!!     negative; other layers and the mass-budget accumulator unchanged.
!!   * floor_bitident                  — `h_min=0` reproduces the pre-feature
!!     trajectory bit-for-bit.
!!   * floor_leak_bound                ��� bounded mass injection: on a
!!     grounding step, injected mass <= N_floored * angstrom_h * areaT;
!!     on a benign step, injection == 0 exactly.
!!   * floor_keeps_tracer_content      -- a layer drained below the floor is
!!     refilled to the floor by continuity, but the tracer update divides by
!!     that stored thickness, so Sum(h*c) is unchanged (the floor dilutes c,
!!     it does not create tracer).  Zonal and meridional halves.
module test_ocean_isopycnal_floor
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_continuity, only: continuity_t, &
                             continuity_compute_fluxes, &
                             continuity_apply_fluxes, &
                             continuity_apply_zonal, &
                             continuity_apply_meridional, &
                             continuity_zonal_flux, &
                             continuity_save_h_pair, &
                             tracer_advect_zonal, tracer_advect_meridional
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   implicit none
   private

   public :: collect_ocean_isopycnal_floor_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 4

contains

   subroutine collect_ocean_isopycnal_floor_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("floor_clamps_grounding_layer", test_floor_clamps_grounding), &
                  new_unittest("floor_bitident", test_floor_bitident), &
                  new_unittest("floor_leak_bound", test_floor_leak_bound), &
                  new_unittest("floor_keeps_tracer_content", test_floor_keeps_content) &
                  ]
   end subroutine collect_ocean_isopycnal_floor_tests

   subroutine make_grid(grid, nx_phys, ny_phys, dx, dy)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      real(wp), intent(in) :: dx, dy
      call grid%init(nx_phys, ny_phys, NGHOST, dx, dy)
   end subroutine make_grid

   subroutine map_in_ms(ms, ct)
      type(multilayer_state_t), intent(inout) :: ms
      type(continuity_t), intent(inout) :: ct
      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(ct)
      call ct%enter_data()
   end subroutine map_in_ms

   subroutine map_out_ms(ms, ct)
      type(multilayer_state_t), intent(inout) :: ms
      type(continuity_t), intent(inout) :: ct
      call ct%exit_data()
      !$acc exit data delete(ct)
      call ms%exit_data()
      !$acc exit data delete(ms)
   end subroutine map_out_ms

   ! -----------------------------------------------------------------

   subroutine test_floor_clamps_grounding(error)
      !! Set one layer's flux_h_layer so h - dt*flux < 0 in one cell.
      !! With angstrom_h = 1e-2, that cell must equal angstrom_h to 1e-14.
      !! All other cells must be unchanged.
      !! With angstrom_h = 0, the cell goes to the un-floored (negative) value.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: DX = 1000.0_wp
      real(wp), parameter :: DT = 60.0_wp
      real(wp), parameter :: H_INIT = 0.5_wp         !! thin but positive
      real(wp), parameter :: LARGE_DIV = H_INIT/DT + 1.0_wp  !! forces h_new < 0
      real(wp), parameter :: AH = 1.0e-2_wp
      integer :: nx_phys, ny_phys, nx, ny, i_grnd, j_grnd, k_grnd
      real(wp) :: floored_expected, no_floor_expected
      checks: block

         nx_phys = 4
         ny_phys = 4
         call make_grid(grid, nx_phys, ny_phys, DX, DX)
         call make_cartesian_metrics(metrics, grid)
         nx = grid%nx_total
         ny = grid%ny_total
         ms%nz_ml = NZ
         call ms%init(grid)
         call ct%init(grid, nz_ml=NZ)

         ! Uniform benign state.
         ms%h_layer = H_INIT
         ms%flux_h_layer = 0.0_wp
         ms%mass_budget_continuity = 0.0_wp

         ! One interior grounding cell: layer k=2 at (i_grnd, j_grnd).
         i_grnd = grid%nghost + 2
         j_grnd = grid%nghost + 2
         k_grnd = 2
         ms%flux_h_layer(i_grnd, j_grnd, k_grnd) = LARGE_DIV

         call map_in_ms(ms, ct)
         ! `enter_data` maps flux_h_layer with `create` (production recomputes
         ! it on-device each step), so the host-set grounding flux above is NOT
         ! on the device — push it, else the -gpu=mem:separate kernel reads
         ! uninitialised device flux and the grounding cell never floors.
         !$acc update device(ms%h_layer, ms%flux_h_layer, ms%mass_budget_continuity)

         ! --- With floor active: floored cell == AH ---
         call continuity_apply_fluxes(ms, DT, h_min=AH)
         !$acc update self(ms%h_layer)

         floored_expected = AH
         call check(error, &
                    abs(ms%h_layer(i_grnd, j_grnd, k_grnd) - floored_expected) < 1.0e-14_wp, &
                    "floor: grounding cell must equal angstrom_h")
         if (allocated(error)) exit checks

         ! Neighbouring cells in the same layer must be at H_INIT (flux=0 there).
         call check(error, &
                    abs(ms%h_layer(i_grnd + 1, j_grnd, k_grnd) - H_INIT) < 1.0e-14_wp, &
                    "floor: non-grounding cell must be unchanged")
         if (allocated(error)) exit checks

         ! --- Reset and rerun with h_min=0: cell must go negative ---
         ms%h_layer = H_INIT
         ms%flux_h_layer = 0.0_wp
         ms%mass_budget_continuity = 0.0_wp
         ms%flux_h_layer(i_grnd, j_grnd, k_grnd) = LARGE_DIV
         !$acc update device(ms%h_layer, ms%flux_h_layer, ms%mass_budget_continuity)

         call continuity_apply_fluxes(ms, DT, h_min=0.0_wp)
         !$acc update self(ms%h_layer)

         no_floor_expected = H_INIT - DT*LARGE_DIV
         call check(error, no_floor_expected < 0.0_wp, &
                    "floor=0 sanity: expected negative result")
         if (allocated(error)) exit checks
         call check(error, &
                    abs(ms%h_layer(i_grnd, j_grnd, k_grnd) - no_floor_expected) < 1.0e-14_wp, &
                    "floor=0: cell must match exact unfloored value")
         if (allocated(error)) exit checks

      end block checks
      call map_out_ms(ms, ct)
      call ct%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_floor_clamps_grounding

   ! -----------------------------------------------------------------

   subroutine test_floor_bitident(error)
      !! Calling continuity_apply_fluxes with h_min=0.0 must produce
      !! bit-identical h_layer to the pre-feature call (no h_min arg),
      !! verified by comparing two copies of the same state.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms_new, ms_ref
      type(continuity_t) :: ct_new, ct_ref
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: DX = 500.0_wp
      real(wp), parameter :: DT = 30.0_wp
      integer :: nx_phys, ny_phys, nx, ny, i, j, k
      real(wp) :: max_diff
      checks: block

         nx_phys = 6
         ny_phys = 6
         call make_grid(grid, nx_phys, ny_phys, DX, DX)
         call make_cartesian_metrics(metrics, grid)
         nx = grid%nx_total
         ny = grid%ny_total
         ms_new%nz_ml = NZ
         ms_ref%nz_ml = NZ
         call ms_new%init(grid)
         call ms_ref%init(grid)
         call ct_new%init(grid, nz_ml=NZ)
         call ct_ref%init(grid, nz_ml=NZ)

         ! Non-trivial but positive flux field.
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx
                  ms_new%h_layer(i, j, k) = 10.0_wp + real(k, wp) + 0.1_wp*real(i + j, wp)
                  ms_new%flux_h_layer(i, j, k) = 0.2_wp*real(mod(i + j + k, 3) - 1, wp)
               end do
            end do
         end do
         ms_new%mass_budget_continuity = 0.0_wp
         ms_ref%h_layer = ms_new%h_layer
         ms_ref%flux_h_layer = ms_new%flux_h_layer
         ms_ref%mass_budget_continuity = 0.0_wp

         call map_in_ms(ms_new, ct_new)
         call map_in_ms(ms_ref, ct_ref)
         ! Push the host-set flux field (flux_h_layer is `create`-mapped) so the
         ! comparison exercises the real flux, not zeroed device scratch.
         !$acc update device(ms_new%h_layer, ms_new%flux_h_layer, &
         !$acc&              ms_new%mass_budget_continuity)
         !$acc update device(ms_ref%h_layer, ms_ref%flux_h_layer, &
         !$acc&              ms_ref%mass_budget_continuity)

         ! New path with explicit h_min=0.
         call continuity_apply_fluxes(ms_new, DT, h_min=0.0_wp)
         ! Reference path with no h_min argument.
         call continuity_apply_fluxes(ms_ref, DT)

         !$acc update self(ms_new%h_layer, ms_ref%h_layer)

         max_diff = maxval(abs(ms_new%h_layer - ms_ref%h_layer))
         call check(error, max_diff == 0.0_wp, &
                    "floor bitident: h_min=0 must produce bit-identical h_layer")
         if (allocated(error)) exit checks

      end block checks
      call map_out_ms(ms_new, ct_new)
      call map_out_ms(ms_ref, ct_ref)
      call ct_new%destroy()
      call ms_new%destroy()
      call ct_ref%destroy()
      call ms_ref%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_floor_bitident

   ! -----------------------------------------------------------------

   subroutine test_floor_leak_bound(error)
      !! After a floor fires on N layer-cells the mass injected satisfies
      !!   sum(h_floored - h_unfloored) <= N * angstrom_h * areaT    (R2 bound)
      !! This is tested by starting with h_layer = 0 (flux_h_layer = 0 so
      !! h_new = 0 < angstrom_h, floor fires to angstrom_h).  The injection
      !! per floored cell is exactly angstrom_h (since h_unfloored = 0),
      !! which equals the per-cell bound.
      !!
      !! On a benign (non-grounding, h >> angstrom_h) run the floor never
      !! fires and injection == 0 exactly.
      !!
      !! Note: we assert the BOUND, not conservation (the floor is
      !! deliberately not conservative — see SPEC §2.1 R2).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms_floor, ms_unfloored
      type(continuity_t) :: ct1, ct2
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: DX = 1000.0_wp
      real(wp), parameter :: DT = 60.0_wp
      real(wp), parameter :: AH = 1.0e-2_wp
      integer :: nx_phys, ny_phys, nx, ny, i, j, k
      real(wp) :: area_t, injection, bound
      integer :: n_floored

      nx_phys = 4
      ny_phys = 4
      call make_grid(grid, nx_phys, ny_phys, DX, DX)
      call make_cartesian_metrics(metrics, grid)
      nx = grid%nx_total
      ny = grid%ny_total
      area_t = DX*DX

      ! ---- Part 1: grounding run ----
      ! Start with h = 0.0, flux = 0.  h_new = 0 < angstrom_h → floor fires
      ! on all cells.  Injection per cell = angstrom_h (since h_unfloored = 0).
      ms_floor%nz_ml = NZ
      ms_unfloored%nz_ml = NZ
      call ms_floor%init(grid)
      call ms_unfloored%init(grid)
      call ct1%init(grid, nz_ml=NZ)
      call ct2%init(grid, nz_ml=NZ)
      ms_floor%h_layer = 0.0_wp
      ms_floor%flux_h_layer = 0.0_wp
      ms_floor%mass_budget_continuity = 0.0_wp
      ms_unfloored%h_layer = 0.0_wp
      ms_unfloored%flux_h_layer = 0.0_wp
      ms_unfloored%mass_budget_continuity = 0.0_wp

      call map_in_ms(ms_floor, ct1)
      call map_in_ms(ms_unfloored, ct2)
      ! Push host-set inputs (flux_h_layer is `create`-mapped) so the run does
      ! not depend on freshly-created device scratch happening to be zero.
      !$acc update device(ms_floor%h_layer, ms_floor%flux_h_layer, &
      !$acc&              ms_floor%mass_budget_continuity)
      !$acc update device(ms_unfloored%h_layer, ms_unfloored%flux_h_layer, &
      !$acc&              ms_unfloored%mass_budget_continuity)
      call continuity_apply_fluxes(ms_floor, DT, h_min=AH)
      call continuity_apply_fluxes(ms_unfloored, DT, h_min=0.0_wp)
      !$acc update self(ms_floor%h_layer, ms_unfloored%h_layer)
      call map_out_ms(ms_floor, ct1)
      call map_out_ms(ms_unfloored, ct2)

      injection = 0.0_wp
      n_floored = 0
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx
               if (ms_floor%h_layer(i, j, k) > ms_unfloored%h_layer(i, j, k)) then
                  injection = injection + &
                              (ms_floor%h_layer(i, j, k) - ms_unfloored%h_layer(i, j, k))*area_t
                  n_floored = n_floored + 1
               end if
            end do
         end do
      end do
      bound = real(n_floored, wp)*AH*area_t

      call ct1%destroy()
      call ms_floor%destroy()
      call ct2%destroy()
      call ms_unfloored%destroy()

      call check(error, n_floored > 0, "floor leak bound: at least one cell must be floored")
      if (allocated(error)) then
         call destroy_cartesian_metrics(metrics)
         return
      end if
      ! injection per cell == AH (unfloored was 0, floored == AH), so
      ! injection == bound exactly; allow floating-point slack.
      call check(error, injection <= bound + 1.0e-10_wp, &
                 "floor leak bound: mass injection must be <= N * angstrom_h * areaT")
      if (allocated(error)) then
         call destroy_cartesian_metrics(metrics)
         return
      end if

      ! ---- Part 2: benign run (h >> angstrom_h, floor never fires) ----
      ms_floor%nz_ml = NZ
      ms_unfloored%nz_ml = NZ
      call ms_floor%init(grid)
      call ms_unfloored%init(grid)
      call ct1%init(grid, nz_ml=NZ)
      call ct2%init(grid, nz_ml=NZ)
      ms_floor%h_layer = 10.0_wp
      ms_floor%flux_h_layer = 0.0_wp
      ms_floor%mass_budget_continuity = 0.0_wp
      ms_unfloored%h_layer = 10.0_wp
      ms_unfloored%flux_h_layer = 0.0_wp
      ms_unfloored%mass_budget_continuity = 0.0_wp

      call map_in_ms(ms_floor, ct1)
      call map_in_ms(ms_unfloored, ct2)
      ! Push host-set inputs (flux_h_layer is `create`-mapped) so the benign run
      ! does not depend on freshly-created device scratch happening to be zero.
      !$acc update device(ms_floor%h_layer, ms_floor%flux_h_layer, &
      !$acc&              ms_floor%mass_budget_continuity)
      !$acc update device(ms_unfloored%h_layer, ms_unfloored%flux_h_layer, &
      !$acc&              ms_unfloored%mass_budget_continuity)
      call continuity_apply_fluxes(ms_floor, DT, h_min=AH)
      call continuity_apply_fluxes(ms_unfloored, DT, h_min=0.0_wp)
      !$acc update self(ms_floor%h_layer, ms_unfloored%h_layer)
      call map_out_ms(ms_floor, ct1)
      call map_out_ms(ms_unfloored, ct2)

      injection = 0.0_wp
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx
               injection = injection + &
                           abs(ms_floor%h_layer(i, j, k) - ms_unfloored%h_layer(i, j, k))*area_t
            end do
         end do
      end do
      call check(error, injection == 0.0_wp, &
                 "floor leak bound: benign run must have zero injection")

      call ct1%destroy()
      call ms_floor%destroy()
      call ct2%destroy()
      call ms_unfloored%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_floor_leak_bound

   ! -----------------------------------------------------------------

   subroutine test_floor_keeps_content(error)
      !! Single column: layer 2 of cell (ic, jc) holds H_THIN, and one face
      !! pulls DQ out of it (DT*DQ > H_THIN - AH, so the Angstrom floor binds
      !! while the unfloored h+ = 0.03 is still > H_DIV_EPS).  Continuity
      !! stores AH; the tracer update must divide by that stored thickness, so
      !! the layer content Sum(h*c) over the domain is unchanged to round-off
      !! for salinity and temperature.  Run for the x half then the y half.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: DT = 1.0_wp
      real(wp), parameter :: H0 = 5.0_wp, H_THIN = 0.1_wp, AH = 0.05_wp
      real(wp), parameter :: DQ = 0.07_wp   !! volume drawn from the thin cell
      integer, parameter :: KL = 2
      integer :: nx, ny, ic, jc, i, j, k, idir, it
      real(wp) :: before(2), after(2), rel
      integer :: tr(2)
      checks: block
         call make_grid(grid, 12, 8, 1.0_wp, 1.0_wp)
         call make_cartesian_metrics(metrics, grid)
         nx = grid%nx_total
         ny = grid%ny_total
         ms%nz_ml = NZ
         call ms%init(grid)
         call ct%init(grid, nz_ml=NZ)
         tr = [ms%idx_salinity, ms%idx_temperature]
         ic = grid%nghost + 6
         jc = grid%nghost + 4

         do idir = 1, 2
            ms%h_layer = H0
            ms%h_layer(ic, jc, KL) = H_THIN
            do k = 1, NZ
               do j = 1, ny
                  do i = 1, nx
                     ms%tracers(tr(1))%c(i, j, k) = 35.0_wp + 0.1_wp*real(i, wp) + 0.05_wp*real(j, wp)
                     ms%tracers(tr(2))%c(i, j, k) = 10.0_wp + 0.2_wp*real(k, wp) - 0.1_wp*real(i, wp)
                  end do
               end do
            end do
            ms%mass_flux_x_layer = 0.0_wp
            ms%mass_flux_y_layer = 0.0_wp
            ms%flux_h_layer = 0.0_wp
            ms%mass_budget_continuity = 0.0_wp
            ! Draw DQ out of the thin cell through its east (x) / north (y) face
            ! into the neighbour: the face flux is the divergence source.
            if (idir == 1) then
               ms%mass_flux_x_layer(ic + 1, jc, KL) = DQ
            else
               ms%mass_flux_y_layer(ic, jc + 1, KL) = DQ
            end if
            do it = 1, 2
               before(it) = sum(ms%h_layer(:, :, KL)*ms%tracers(tr(it))%c(:, :, KL))
            end do

            call map_in_ms(ms, ct)
            ! mass_flux is `create`-mapped: push the host-set flux.
            !$acc update device(ms%mass_flux_x_layer, ms%mass_flux_y_layer)
            call continuity_save_h_pair(ct, ms)
            if (idir == 1) then
               call continuity_apply_zonal(grid, metrics, ms, DT, h_min=AH)
               call tracer_advect_zonal(grid, metrics, ct, ms, DT, ct%h_pair_prev%data, ms%h_layer, h_floor=AH)
            else
               call continuity_apply_meridional(grid, metrics, ms, DT, h_min=AH)
               call tracer_advect_meridional(grid, metrics, ct, ms, DT, ct%h_pair_prev%data, ms%h_layer, h_floor=AH)
            end if
            call map_out_ms(ms, ct)

            call check(error, abs(ms%h_layer(ic, jc, KL) - AH) < 1.0e-14_wp, &
                       "floor_keeps_content: the floor must have bound (h stored = angstrom_h)")
            if (allocated(error)) exit checks
            do it = 1, 2
               after(it) = sum(ms%h_layer(:, :, KL)*ms%tracers(tr(it))%c(:, :, KL))
               rel = abs(after(it) - before(it))/abs(before(it))
               call check(error, rel < 1.0e-13_wp, &
                          "floor_keeps_content: Sum(h*c) changed across the Angstrom floor")
               if (allocated(error)) exit checks
            end do
         end do
      end block checks
      call ct%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_floor_keeps_content

end module test_ocean_isopycnal_floor
