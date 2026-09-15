!! Default-off guard for the surface-pressure loading seam (PR-17).
!!
!! `&ocean_psurf_nml enable=.false.` must be byte-for-byte identical to the
!! path that never carried the knob at all.  This exercises the 3-way stage
!! dispatch in `ocean_dyn_step_split` at the level where PR-17 changed it —
!! ABOVE the barotropic substep kernel that `test_ocean_tides_disabled_bitident`
!! already covers.  The bug PR-17 can introduce (a mis-ordered 3-way branch,
!! an unmapped eta_seam) lives here, where the kernel-level test cannot see it.
!!
!! Invariant: `ocean_dyn_step_split` with `psurf` PRESENT-but-disabled
!! produces byte-for-byte identical bt_eta / u_face_x_layer / v_face_y_layer
!! to the same call with `psurf` ABSENT, over N steps from a non-trivial IC.
!!
!! GPU/mem:separate: state + scratch companions are enter_data'd; fields read
!! back are update self'd; directives are inert on multicore.  Verify on GPU.
module test_ocean_p_surf_bitident
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_continuity, only: continuity_t
   use rdb_coriolis_adv, only: coriolis_adv_t
   use rdb_eos, only: eos_t
   use rdb_ocean_pressure_force, only: ocean_pressure_force_t
   use rdb_ocean_horizontal_viscosity, only: ocean_horizontal_viscosity_t
   use rdb_ocean_bottom_drag, only: ocean_bottom_drag_t
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t
   use rdb_ocean_vertical_advection, only: ocean_vertical_advection_t
   use rdb_ocean_hdiff_tracer, only: ocean_hdiff_tracer_t
   use rdb_ocean_vdiff, only: ocean_vdiff_t
   use rdb_ocean_vmix, only: ocean_vmix_t
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t
   use rdb_ocean_dyn, only: ocean_dyn_t, ocean_dyn_step_split
   use rdb_ocean_p_surf, only: ocean_p_surf_t, p_surf_configure
   implicit none
   private

   public :: collect_ocean_p_surf_bitident_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NX_PHYS = 16, NY_PHYS = 16, NZ = 2
   real(wp), parameter :: DX = 50000.0_wp, DY = 50000.0_wp
   real(wp), parameter :: H_TOTAL = 4000.0_wp
   real(wp), parameter :: DT = 300.0_wp
   integer, parameter :: N_INNER = 16
   integer, parameter :: N_STEPS = 20

contains

   subroutine collect_ocean_p_surf_bitident_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("psurf_off_equals_absent", test_off_equals_absent)]
   end subroutine collect_ocean_p_surf_bitident_tests

   subroutine test_off_equals_absent(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: eta_off(:, :), u_off(:, :, :), v_off(:, :, :)
      real(wp), allocatable :: eta_abs(:, :), u_abs(:, :, :), v_abs(:, :, :)
      real(wp) :: dmax

      ! psurf present but enable=.false.
      call run(error, with_psurf=.true., eta_out=eta_off, u_out=u_off, v_out=v_off)
      if (allocated(error)) return
      ! psurf absent entirely.
      call run(error, with_psurf=.false., eta_out=eta_abs, u_out=u_abs, v_out=v_abs)
      if (allocated(error)) return

      dmax = maxval(abs(eta_off - eta_abs))
      dmax = max(dmax, maxval(abs(u_off - u_abs)))
      dmax = max(dmax, maxval(abs(v_off - v_abs)))
      call check(error, dmax == 0.0_wp, &
                 "psurf present-but-disabled must be byte-for-byte identical "// &
                 "to psurf absent (default-off bit-identity)")

      deallocate (eta_off, u_off, v_off, eta_abs, u_abs, v_abs)
   end subroutine test_off_equals_absent

   subroutine run(error, with_psurf, eta_out, u_out, v_out)
      type(error_type), allocatable, intent(out) :: error
      logical, intent(in) :: with_psurf
      real(wp), allocatable, intent(out) :: eta_out(:, :), u_out(:, :, :), v_out(:, :, :)
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(coriolis_adv_t) :: cor
      type(ocean_pressure_force_t) :: pgf
      type(ocean_horizontal_viscosity_t) :: hv
      type(ocean_bottom_drag_t) :: bd
      type(ocean_surface_stress_t) :: ss
      type(ocean_vertical_advection_t) :: va
      type(ocean_hdiff_tracer_t) :: hd
      type(ocean_vdiff_t) :: vd
      type(ocean_vmix_t) :: vmix
      type(eos_t) :: eos
      type(ocean_dyn_t) :: dyn
      type(ocean_p_surf_t) :: psurf
      integer :: nx, ny, i, j, k, step
      real(wp) :: tilt

      call grid%init(NX_PHYS, NY_PHYS, NGHOST, DX, DY)
      nx = grid%nx_total
      ny = grid%ny_total
      call make_cartesian_metrics(metrics, grid)

      ms%nz_ml = NZ
      call ms%init(grid)
      call ct%init(grid, nz_ml=NZ)
      cor%f_0 = 1.0e-4_wp
      call cor%init(grid, nz_ml=NZ)
      call pgf%init(grid, nz_ml=NZ)
      call hv%init(grid, nz_ml=NZ)
      call bd%init(grid, nz_ml=NZ)
      call ss%init(grid, nz_ml=NZ)
      call va%init(grid, nz_ml=NZ)
      call hd%init(grid, nz_ml=NZ)
      call vd%init(grid, nz_ml=NZ)
      call vmix%init(grid, nz_ml=NZ)
      call eos%init(grid)
      call dyn%init(grid, nz_ml=NZ)
      dyn%bt_work%bt_H_ref = H_TOTAL

      ! Non-trivial, dynamically active IC.
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx
               tilt = 0.05_wp*real(i - nx/2, wp)/real(nx, wp)
               ms%h_layer(i, j, k) = (H_TOTAL + tilt)/real(NZ, wp)
               ms%tracers(ms%idx_salinity)%hTr(i, j, k) = &
                  eos%S_ref*ms%h_layer(i, j, k)
               ms%tracers(ms%idx_temperature)%hTr(i, j, k) = &
                  eos%T_ref*ms%h_layer(i, j, k)
               ms%u_face_x_layer(i, j, k) = &
                  0.02_wp*sin(2.0_wp*acos(-1.0_wp)*real(j, wp)/real(ny, wp))
               ms%v_face_y_layer(i, j, k) = 0.0_wp
            end do
         end do
      end do

      ! Disabled: enable=.false. => the slot allocates nothing and the
      ! 3-way dispatch must fall through to the pre-PR-17 branch.
      psurf%enable = .false.
      call p_surf_configure(psurf, nx, ny)

      !$acc enter data copyin(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn, psurf)
      call ms%enter_data()
      call ct%enter_data()
      call cor%enter_data()
      call pgf%enter_data()
      call hv%enter_data()
      call bd%enter_data()
      call ss%enter_data()
      call va%enter_data()
      call hd%enter_data()
      call vd%enter_data()
      call vmix%enter_data()
      call dyn%enter_data()
      call psurf%enter_data()

      do step = 1, N_STEPS
         if (with_psurf) then
            call ocean_dyn_step_split(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
                                      va, hd, vd, vmix, ms, DT, N_INNER, psurf=psurf)
         else
            call ocean_dyn_step_split(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
                                      va, hd, vd, vmix, ms, DT, N_INNER)
         end if
      end do

      !$acc update self(ms%u_face_x_layer, ms%v_face_y_layer, dyn%bt_work%bt_eta)
      ! C-grid faces are staggered (u has nx+1 in dim 1, v has ny+1 in
      ! dim 2), so take the arrays' own shapes via source-only allocation.
      allocate (eta_out, source=dyn%bt_work%bt_eta)
      allocate (u_out, source=ms%u_face_x_layer)
      allocate (v_out, source=ms%v_face_y_layer)

      call psurf%exit_data()
      call dyn%exit_data()
      call vmix%exit_data()
      call vd%exit_data()
      call hd%exit_data()
      call va%exit_data()
      call ss%exit_data()
      call bd%exit_data()
      call hv%exit_data()
      call pgf%exit_data()
      call cor%exit_data()
      call ct%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn, psurf)
      call metrics%exit_data()

      call check(error, all(eta_out == eta_out), "eta_out is NaN")

      call psurf%destroy()
      call dyn%destroy()
      call eos%destroy()
      call vmix%destroy()
      call vd%destroy()
      call hd%destroy()
      call va%destroy()
      call ss%destroy()
      call bd%destroy()
      call hv%destroy()
      call pgf%destroy()
      call cor%destroy()
      call ct%destroy()
      call ms%destroy()
   end subroutine run

end module test_ocean_p_surf_bitident
