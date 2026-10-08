!! [C1] ACCESS-OM3 vertical-grid parity: run on OM3's exact 75-level z*
!! stack, and exercise the diag output-level bound raised for it.
!!
!! `tools/vgrid_to_dz.py` converts OM3's `ocean_vgrid.nc` (75 levels, only
!! available on Gadi -- see its own self-test for synthetic coverage) into
!! a `&vcoord_nml z_fixed_profile = "list", z_fixed_dz = ...` block. This
!! suite checks the TWO things that block needs from the engine once
!! pasted in:
!!
!!   * `om3_75_level_zstar_steps_cleanly` -- a 75-entry `z_fixed_dz` list
!!     under `vcoord_type = "zstar"` (OM3's coordinate) builds and steps.
!!     The profile is a synthetic surface-first stretch (10 m at the
!!     surface, +2 m per layer, summing to 6300 m over 75 layers -- the
!!     same profile `tools/vgrid_to_dz.py --self-test` checks its own
!!     conversion against), so a regression here and a regression in the
!!     tool's self-test would show up at the same numbers.
!!   * `om3_diag_z_levels_over_64_conserves` -- `MAX_OCEAN_DIAG_Z_LEVELS`
!!     was raised 64 -> 128 for this plan item; this registers a
!!     `DIAG_VGRID_Z_FIXED` diagnostic on 100 output levels (> the old
!!     bound) through the SAME device pipeline (`enter_data` / `step` /
!!     `exit_data`) `test_ocean_diag_remap.F90::test_on_device` uses, and
!!     checks the remap is finite and conserves the column integral --
!!     the `mem:separate` gotcha this project's CLAUDE.md warns about
!!     (every array the kernel touches must be device-present).
module test_ocean_vgrid_om3
   use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_null_ptr, c_f_pointer
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_state, only: ocean_state_t, ocean_state_enter_data, ocean_state_exit_data
   use rdb_ocean_diag, only: DIAG_OP_INSTANT, DIAG_VGRID_Z_FIXED
   use rdb_ocean_diag_fills, only: fill_temperature, remap_layer_to_z
   use rdb_ocean_metrics, only: metrics_fill_cartesian, metrics_finalize
   use rdb_ocean_api, only: rdb_ocean_create_from_string, rdb_ocean_step, &
                            rdb_ocean_destroy, rdb_ocean_refresh_host, &
                            rdb_ocean_get_grid_info, rdb_ocean_get_h_layer_ptr
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   use testdrive, only: error_type, check, new_unittest, unittest_type
   implicit none
   private

   public :: collect_ocean_vgrid_om3_tests

   integer, parameter :: NZ75 = 75
      !! ACCESS-OM3's layer count.
   real(wp), parameter :: MAX_DEPTH_75 = 6300.0_wp
      !! Sum of the synthetic 75-entry profile below (10 + 2*(k-1), k=1..75).

   integer, parameter :: NX_D = 5, NY_D = 4, NZ_D = 4
      !! Small grid for the diag-level-bound test (unrelated to NZ75).
   real(wp), parameter :: DX_D = 1.0_wp

contains

   subroutine collect_ocean_vgrid_om3_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("om3_75_level_zstar_steps_cleanly", test_om3_75_zstar), &
                  new_unittest("om3_diag_z_levels_over_64_conserves", test_diag_over64) &
                  ]
   end subroutine collect_ocean_vgrid_om3_tests

   ! ------------------------------------------------------------------
   ! (a) 75-entry z_fixed_dz list, vcoord_type="zstar", a few steps.
   ! ------------------------------------------------------------------

   pure function dz75_surface_first() result(dz)
      !! Surface-first synthetic OM3-like profile (fine near the surface,
      !! coarsening with depth): 10, 12, 14, ..., 158 m, sum = 6300 m.
      !! Matches `vgrid_to_dz._expected_dz(75)` in
      !! `tools/vgrid_to_dz.py` so both regressions are pinned to the same
      !! numbers.
      real(wp) :: dz(NZ75)
      integer :: k
      do k = 1, NZ75
         dz(k) = 10.0_wp + 2.0_wp*real(k - 1, wp)
      end do
   end function dz75_surface_first

   pure function dz_list_nml(dz) result(txt)
      !! Comma-separated namelist literal for a real array, e.g.
      !! "10.000, 12.000, ...".
      real(wp), intent(in) :: dz(:)
      character(len=:), allocatable :: txt
      character(len=24) :: tmp
      integer :: k
      txt = ""
      do k = 1, size(dz)
         write (tmp, "(f12.4)") dz(k)
         if (k < size(dz)) then
            txt = txt//trim(adjustl(tmp))//", "
         else
            txt = txt//trim(adjustl(tmp))
         end if
      end do
   end function dz_list_nml

   subroutine test_om3_75_zstar(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), parameter :: nl = new_line("a")
      character(len=:), allocatable :: nml, dzlist
      type(c_ptr) :: handle, ptr
      integer(c_int) :: status, nx_p, ny_p, nz_p, ng, nx, ny, nz, gen
      real(wp), pointer :: h(:, :, :)
      real(wp) :: col_total, err_total
      character(len=160) :: msg

      dzlist = dz_list_nml(dz75_surface_first())
      nml = '&sim_nml sim_type = "ocean" /'//nl// &
            "&grid_nml nx = 4, ny = 4, dx = 2000.0, dy = 2000.0, nghost = 2 /"//nl// &
            "&time_nml t_end = 900.0, dt_fixed = 300.0 /"//nl// &
            "&physics_nml coriolis_f = -1.0e-4 /"//nl// &
            "&nonhydrostatic_nml nz_layers = 75 /"//nl// &
            '&vcoord_nml vcoord_type = "zstar", z_fixed_profile = "list", '// &
            "z_fixed_dz = "//dzlist//" /"//nl// &
            '&ocean_topo_nml topo_config = "flat", max_depth = 6300.0 /'//nl// &
            "&ocean_vmix_nml use_closure = .false., use_kpp = .false. /"//nl// &
            "&ocean_bt_nml auto_n_inner = .true. /"//nl// &
            "&ocean_diag_nml enabled = .false. /"//nl// &
            "&output_nml output_to_file = .false. /"//nl

      handle = c_null_ptr
      status = rdb_ocean_create_from_string(nml, len(nml, kind=c_int), handle)
      call check(error, status == OCEAN_STATUS_OK, "75-level zstar case must build")
      if (allocated(error)) return

      body: block
         status = rdb_ocean_get_grid_info(handle, nx_p, ny_p, nz_p, ng)
         call check(error, nz_p == NZ75, "engine must resolve nz = 75")
         if (allocated(error)) exit body

         status = rdb_ocean_step(handle, 3_c_int)
         call check(error, status == OCEAN_STATUS_OK, "75-level zstar case must step")
         if (allocated(error)) exit body

         status = rdb_ocean_refresh_host(handle)
         status = rdb_ocean_get_h_layer_ptr(handle, ptr, nx, ny, nz, gen)
         call c_f_pointer(ptr, h, [nx, ny, nz])
         col_total = sum(h(ng + 2, ng + 2, :))
         err_total = abs(col_total - MAX_DEPTH_75)
         write (msg, "(a,es11.4,a,es11.4,a)") &
            "zstar column total should track max_depth=6300 m (got ", &
            col_total, ", |err| ", err_total, " m)"
         ! zstar tracks eta, so a small resting-state deviation is fine; a
         ! 1% band catches a scrambled profile/ordering, not round-off.
         call check(error, err_total < 0.01_wp*MAX_DEPTH_75, trim(msg))
      end block body
      status = rdb_ocean_destroy(handle)
   end subroutine test_om3_75_zstar

   ! ------------------------------------------------------------------
   ! (b) diag z-level output above the old 64-entry bound.
   ! ------------------------------------------------------------------

   subroutine setup_diag_state(grid, state)
      type(hgrid_t), intent(inout) :: grid
      type(ocean_state_t), intent(inout) :: state
      call grid%init(NX_D, NY_D, 1, DX_D, DX_D)
      state%multilayer%nz_ml = NZ_D
      call state%init(grid)
      call metrics_fill_cartesian(state%metrics, grid, grid%dx, grid%dy)
      call metrics_finalize(state%metrics)
      state%eos%rho0 = 1025.0_wp
      state%eos%alpha_T = 1.0_wp
      state%eos%beta_S = 0.0_wp
      state%eos%T_ref = 0.0_wp
      state%eos%S_ref = 0.0_wp
   end subroutine setup_diag_state

   subroutine set_uniform_TS_d(state, temp_by_layer)
      type(ocean_state_t), intent(inout) :: state
      real(wp), intent(in) :: temp_by_layer(:)
      integer :: it, is_, k
      it = state%multilayer%idx_temperature
      is_ = state%multilayer%idx_salinity
      state%multilayer%h_layer = 10.0_wp
      do k = 1, NZ_D
         state%multilayer%tracers(it)%hTr(:, :, k) = 10.0_wp*temp_by_layer(k)
         state%multilayer%tracers(is_)%hTr(:, :, k) = 0.0_wp
      end do
   end subroutine set_uniform_TS_d

   pure elemental function finite(x) result(ok)
      real(wp), intent(in) :: x
      logical :: ok
      ok = (x == x) .and. (abs(x) < huge(1.0_wp))
   end function finite

   subroutine test_diag_over64(error)
      !! 100 output z-levels (> the pre-C1 MAX_OCEAN_DIAG_Z_LEVELS=64
      !! bound, <= the post-C1 bound of 128) spanning a 40 m column
      !! (4 layers x 10 m) at 0.4 m spacing. Runs the full device
      !! pipeline (enter_data/step/exit_data): finite output + the
      !! thickness-weighted column integral of T is preserved.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NLEV = 100
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      real(wp) :: z_levels(NLEV)
      real(wp), parameter :: TEMP_BY_LAYER(NZ_D) = [0.0_wp, 1.0_wp, 2.0_wp, 3.0_wp]
      integer :: k, iv, it
      real(wp) :: out_int, src_int, dz_lev

      do k = 1, NLEV
         z_levels(k) = real(k, wp)*0.4_wp   ! 0.4, 0.8, ..., 40.0 m
      end do

      checks: block
         call setup_diag_state(grid, state)
         call set_uniform_TS_d(state, TEMP_BY_LAYER)

         call state%diag%set_output_z_levels(z_levels)
         call state%diag%register("T_z100", units="degC", fill=fill_temperature, &
                                  n1=NX_D, n2=NY_D, n3=NZ_D, &
                                  output_vgrid=DIAG_VGRID_Z_FIXED, &
                                  remap=remap_layer_to_z, &
                                  time_op=DIAG_OP_INSTANT, dt_out=1.0_wp)

         call ocean_state_enter_data(state)
         call state%diag%step(state, dt=2.0_wp, t=2.0_wp)
         call ocean_state_exit_data(state)

         it = 0
         do iv = 1, state%diag%nvars
            if (state%diag%vars(iv)%name == "T_z100") it = iv
         end do
         call check(error, it > 0, "T_z100 must register")
         if (allocated(error)) exit checks

         call check(error, size(state%diag%vars(it)%output_buffer, 3) == NLEV, &
                    "output_buffer must carry all 100 levels")
         if (allocated(error)) exit checks

         call check(error, all(finite(state%diag%vars(it)%output_buffer)), &
                    "100-level z remap output must be finite")
         if (allocated(error)) exit checks

         ! Column integral: source is piecewise-constant T=0,1,2,3 over
         ! 10 m layers -> integral = (0+1+2+3)*10 = 60. The 100-level
         ! target spans the whole 40 m column at uniform 0.4 m spacing,
         ! so the output integral must match.
         dz_lev = 0.4_wp
         out_int = sum(state%diag%vars(it)%output_buffer(2, 2, :))*dz_lev
         src_int = sum(TEMP_BY_LAYER)*10.0_wp
         call check(error, abs(out_int - src_int) < 1.0e-6_wp, &
                    "100-level z remap must preserve the column T integral (60)")
      end block checks
      call state%destroy()
   end subroutine test_diag_over64

end module test_ocean_vgrid_om3
