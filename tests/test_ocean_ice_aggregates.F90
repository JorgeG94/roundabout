!! D1.6 scheme-neutral ice aggregates: `aice`, `vice`, `vsno`, `aice0` of
!! `ocean_sea_ice_t` (filled by `ice_aggregates_fill`) must reproduce, BITWISE
!! on the CPU build, the inline category gathers the consumers used before the
!! retarget -- the old expressions are copied below as the reference:
!!   - diag `aice`/`ice_conc`, `hi`/`ice_thick` (`fill_ice_conc_thick_impl`)
!!   - diag `hs` (`fill_ice_hs_impl`)
!!   - console ice totals (`compute_ice_totals[_efp]`)
!!   - bulk-forcing open-water weight (`part_size(:,:,0)`)
!! Cells: empty, partially covered, full, over-saturated (sum > 1), land,
!! every cell of the array including the ghost ring.  `aice` (sum of
!! fractions + `min`) and `aice0` (a copy) are compared BITWISE on every
!! build.  `vice`/`vsno` use a reference computed by a `do concurrent` of the
!! test on the same device (same arithmetic form), so device-vs-device is
!! bitwise too; the console totals keep a 1e-12 relative tolerance on the GPU
!! build, where the reduction order differs.  An engine-level test pins the
!! refresh point: after EVERY step of a `dt_therm_ratio > 1` run the
!! aggregates equal a fresh recompute from the category state.
module test_ocean_ice_aggregates
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_state, only: ocean_state_t, ocean_state_enter_data, ocean_state_exit_data
   use rdb_ocean_metrics, only: metrics_fill_cartesian, metrics_finalize
   use rdb_ocean_console_stats, only: compute_ice_totals, compute_ice_totals_efp
   use rdb_efp, only: efp_t, efp_to_real
   use rdb_ice_state, only: ice_aggregates_fill
   use rdb_config, only: config_t, read_config_from_string, validate_config
   use rdb_ocean_engine, only: ocean_engine_t, engine_setup, engine_enter_data, &
                               engine_step, engine_step_ice, engine_step_finalize, &
                               engine_exit_data, engine_teardown
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   use rdb_comm_env, only: comm_env_init, comm_env_setup_roles
   use rdb_ice_column, only: ICE_RHO_ICE, ICE_RHO_SNOW
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use testdrive, only: error_type, check, new_unittest, unittest_type
   implicit none
   private

   public :: collect_ocean_ice_aggregates_tests

   integer, parameter :: NX = 9, NY = 7, NZ = 2, NG = 2

contains

   subroutine collect_ocean_ice_aggregates_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("aggregates_match_inline_ncat1", test_ncat1), &
                  new_unittest("aggregates_match_inline_ncat5", test_ncat5), &
                  new_unittest("aggregates_refresh_after_change", test_refresh), &
                  new_unittest("aice0_lumped_stays_one", test_aice0_lumped), &
                  new_unittest("aggregates_track_state_every_step_engine", test_engine_every_step) &
                  ]
   end subroutine collect_ocean_ice_aggregates_tests

   subroutine setup(grid, state, ncat)
      type(hgrid_t), intent(inout) :: grid
      type(ocean_state_t), intent(inout) :: state
      integer, intent(in) :: ncat
      integer :: i, j, c
      real(wp) :: r
      call grid%init(NX, NY, NG, 1.0_wp, 1.0_wp)
      state%multilayer%nz_ml = NZ
      state%ice%enable = .true.
      state%ice%ncat = ncat
      call state%init(grid)
      call metrics_fill_cartesian(state%metrics, grid, grid%dx, grid%dy)
      call metrics_finalize(state%metrics)
      ! Deterministic, cell-varying state over the WHOLE array (ghosts too).
      do j = 1, size(state%ice%m_ice, 2)
         do i = 1, size(state%ice%m_ice, 1)
            if (ncat == 1) then
               select case (mod(i + 2*j, 4))
               case (0)
                  state%ice%m_ice(i, j, 1) = 0.0_wp               ! empty
               case (1)
                  state%ice%m_ice(i, j, 1) = 300.0_wp + 7.0_wp*real(i, wp)
               case (2)
                  state%ice%m_ice(i, j, 1) = 1.0e-9_wp            ! tiny but > 0
               case default
                  state%ice%m_ice(i, j, 1) = 905.0_wp + real(j, wp)
               end select
               state%ice%m_snow(i, j, 1) = 12.5_wp*real(mod(i + j, 3), wp)
            else
               r = 0.0_wp
               do c = 1, ncat
                  select case (mod(i + 2*j, 5))
                  case (0)
                     state%ice%part_size(i, j, c) = 0.0_wp        ! empty cell
                  case (1)
                     state%ice%part_size(i, j, c) = 0.05_wp*real(c, wp) + 0.003_wp*real(i, wp)
                  case (2)
                     state%ice%part_size(i, j, c) = 1.0_wp/real(ncat, wp)   ! full
                  case (3)
                     state%ice%part_size(i, j, c) = 0.35_wp           ! sum > 1
                  case default
                     state%ice%part_size(i, j, c) = merge(0.1_wp, 0.0_wp, mod(c, 2) == 1)
                  end select
                  r = r + state%ice%part_size(i, j, c)
                  state%ice%m_ice(i, j, c) = 40.0_wp*real(c, wp) + 0.37_wp*real(i*j, wp)
                  state%ice%m_snow(i, j, c) = 3.3_wp*real(c, wp) + 0.011_wp*real(i + j, wp)
               end do
               state%ice%part_size(i, j, 0) = 1.0_wp - r
            end if
         end do
      end do
      ! A few land cells (aggregates are not wet-masked; consumers are).
      state%metrics%wet_T(4, 3) = 0.0_wp
      state%metrics%wet_T(6, 5) = 0.0_wp
   end subroutine setup

   subroutine ref_gather(ncat, nx, ny, part_size, m_ice, m_snow, ref_a, ref_v, ref_s)
      !! The OLD inline expressions (diag fills / console statistics before
      !! D1.6), verbatim operand order, as a `do concurrent` so the reference
      !! compiles through the same device path as the kernel under test.
      integer, intent(in) :: ncat, nx, ny
      real(wp), intent(in) :: part_size(nx, ny, 0:ncat), m_ice(nx, ny, ncat), m_snow(nx, ny, ncat)
      real(wp), intent(inout) :: ref_a(nx, ny), ref_v(nx, ny), ref_s(nx, ny)
      integer :: i, j, c
      real(wp) :: ci_sum, mice_val, hs_val, ci_val, w
      if (ncat == 1) then
         do concurrent(j=1:ny, i=1:nx) local(mice_val, ci_val, hs_val)
            if (m_ice(i, j, 1) > 0.0_wp) then
               mice_val = m_ice(i, j, 1)
               ci_val = 1.0_wp
            else
               mice_val = 0.0_wp
               ci_val = 0.0_wp
            end if
            ref_a(i, j) = ci_val
            ref_v(i, j) = mice_val/ICE_RHO_ICE
            hs_val = 0.0_wp
            hs_val = hs_val + 1.0_wp*m_snow(i, j, 1)
            ref_s(i, j) = hs_val/ICE_RHO_SNOW
         end do
      else
         do concurrent(j=1:ny, i=1:nx) local(c, mice_val, ci_sum, hs_val, w)
            mice_val = 0.0_wp
            ci_sum = 0.0_wp
            do c = 1, ncat
               mice_val = mice_val + part_size(i, j, c)*m_ice(i, j, c)
               ci_sum = ci_sum + part_size(i, j, c)
            end do
            ref_a(i, j) = min(1.0_wp, ci_sum)
            ref_v(i, j) = mice_val/ICE_RHO_ICE
            hs_val = 0.0_wp
            do c = 1, ncat
               w = part_size(i, j, c)
               hs_val = hs_val + w*m_snow(i, j, c)
            end do
            ref_s(i, j) = hs_val/ICE_RHO_SNOW
         end do
      end if
   end subroutine ref_gather

   subroutine agg_mismatches(state, ncat, nbad, ref_a, ref_v, ref_s)
      !! Count the cells (ghosts included) where an aggregate differs, BITWISE,
      !! from the device-computed reference of the current category state.
      !! The caller has `update self`-ed the aggregates already.
      type(ocean_state_t), intent(inout) :: state
      integer, intent(in) :: ncat
      integer, intent(out) :: nbad
      real(wp), allocatable, intent(out) :: ref_a(:, :), ref_v(:, :), ref_s(:, :)
      integer :: i, j, nx, ny
      nx = size(state%ice%m_ice, 1)
      ny = size(state%ice%m_ice, 2)
      allocate (ref_a(nx, ny), ref_v(nx, ny), ref_s(nx, ny), source=0.0_wp)
      !$acc enter data create(ref_a, ref_v, ref_s)
      call ref_gather(ncat, nx, ny, state%ice%part_size, state%ice%m_ice, state%ice%m_snow, &
                      ref_a, ref_v, ref_s)
      !$acc update self(ref_a, ref_v, ref_s, state%ice%part_size)
      !$acc exit data delete(ref_a, ref_v, ref_s)
      nbad = 0
      do j = 1, ny
         do i = 1, nx
            if (.not. (ieee_is_finite(state%ice%aice(i, j)) .and. ieee_is_finite(state%ice%vice(i, j)) &
                       .and. ieee_is_finite(state%ice%vsno(i, j)) &
                       .and. ieee_is_finite(state%ice%aice0(i, j)))) then
               nbad = nbad + 1
               cycle
            end if
            if (state%ice%aice(i, j) /= ref_a(i, j)) nbad = nbad + 1
            if (state%ice%vice(i, j) /= ref_v(i, j)) nbad = nbad + 1
            if (state%ice%vsno(i, j) /= ref_s(i, j)) nbad = nbad + 1
            if (state%ice%aice0(i, j) /= state%ice%part_size(i, j, 0)) nbad = nbad + 1
         end do
      end do
   end subroutine agg_mismatches

   subroutine compare_all(error, state, ncat)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_state_t), intent(inout) :: state
      integer, intent(in) :: ncat
      real(wp), allocatable :: ref_a(:, :), ref_v(:, :), ref_s(:, :)
      real(wp) :: sw, sc, sh, ew, ec, eh
      integer :: i, j, nx, ny, nbad
      type(efp_t) :: wefp, cefp, hefp

      nx = size(state%ice%m_ice, 1)
      ny = size(state%ice%m_ice, 2)
      call agg_mismatches(state, ncat, nbad, ref_a, ref_v, ref_s)
      call check(error, nbad, 0, "every aggregate in every cell (ghosts included) is BITWISE the old inline expression")
      if (allocated(error)) return

      ! ---- Console totals: old inline accumulation vs the retargeted routines.
      sw = 0.0_wp
      sc = 0.0_wp
      sh = 0.0_wp
      do j = NG + 1, ny - NG
         do i = NG + 1, nx - NG
            sw = sw + state%metrics%wet_T(i, j)*state%metrics%areaT(i, j)
            if (state%metrics%wet_T(i, j) > 0.5_wp) then
               sc = sc + ref_a(i, j)*state%metrics%areaT(i, j)
               sh = sh + ref_v(i, j)*state%metrics%areaT(i, j)
            end if
         end do
      end do
      call compute_ice_totals(state%metrics%wet_T, state%metrics%areaT, state%ice%aice, &
                              state%ice%vice, NG, ew, ec, eh)
      call check(error, abs(ew - sw) <= 1.0e-12_wp*abs(sw), "console wet area")
      if (allocated(error)) return
      call check(error, abs(ec - sc) <= 1.0e-12_wp*max(abs(sc), 1.0e-30_wp), "console ci area")
      if (allocated(error)) return
      call check(error, abs(eh - sh) <= 1.0e-12_wp*max(abs(sh), 1.0e-30_wp), "console hi area")
      if (allocated(error)) return
#ifndef RDB_GPU_OFFLOAD
      call check(error, ec == sc .and. eh == sh .and. ew == sw, &
                 "FP console totals bitwise equal on the CPU build")
      if (allocated(error)) return
#endif
      call compute_ice_totals_efp(state%metrics%wet_T, state%metrics%areaT, state%ice%aice, &
                                  state%ice%vice, NG, wefp, cefp, hefp)
      call check(error, abs(real(efp_to_real(cefp), wp) - sc) <= 1.0e-12_wp*max(abs(sc), 1.0e-30_wp), &
                 "EFP ci area")
      if (allocated(error)) return
      call check(error, abs(real(efp_to_real(hefp), wp) - sh) <= 1.0e-12_wp*max(abs(sh), 1.0e-30_wp), &
                 "EFP hi area")
   end subroutine compare_all

   subroutine run_case(error, ncat)
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: ncat
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      checks: block
         call setup(grid, state, ncat)
         call ocean_state_enter_data(state)       ! fills the aggregates on the device
         !$acc update self(state%ice%aice, state%ice%vice, state%ice%vsno, state%ice%aice0)
         call compare_all(error, state, ncat)
         call ocean_state_exit_data(state)
      end block checks
      call state%destroy()
   end subroutine run_case

   subroutine test_ncat1(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 1)
   end subroutine test_ncat1

   subroutine test_ncat5(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, 5)
   end subroutine test_ncat5

   subroutine test_refresh(error)
      !! The aggregates track the state through `ice_aggregates_fill` (the call
      !! the engine makes after every thermo window): change the category
      !! state on host and device, refill, and they follow.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      checks: block
         call setup(grid, state, 3)
         call ocean_state_enter_data(state)
         state%ice%part_size(:, :, 1) = 0.25_wp
         state%ice%part_size(:, :, 2) = 0.25_wp
         state%ice%part_size(:, :, 3) = 0.0_wp
         state%ice%part_size(:, :, 0) = 0.5_wp
         state%ice%m_ice(:, :, 1) = 450.0_wp
         state%ice%m_ice(:, :, 2) = 900.0_wp
         state%ice%m_snow(:, :, 1) = 33.0_wp
         state%ice%m_snow(:, :, 2) = 66.0_wp
         !$acc update device(state%ice%part_size, state%ice%m_ice, state%ice%m_snow)
         call ice_aggregates_fill(state%ice)
         !$acc update self(state%ice%aice, state%ice%vice, state%ice%vsno, state%ice%aice0)
         call compare_all(error, state, 3)
         if (allocated(error)) exit checks
         call check(error, abs(state%ice%aice(3, 3) - 0.5_wp) < 1.0e-14_wp, "aice follows the state")
         if (allocated(error)) exit checks
         call check(error, abs(state%ice%aice0(3, 3) - 0.5_wp) < 1.0e-14_wp, "aice0 follows the state")
         call ocean_state_exit_data(state)
      end block checks
      call state%destroy()
   end subroutine test_refresh

   subroutine test_aice0_lumped(error)
      !! Lumped mode (`ncat == 1`): `part_size` is not maintained, so `aice0`
      !! (a verbatim copy of the stored category 0) stays 1 even where
      !! `aice = 1`.  Inherited behaviour: lumped ice never reduces the
      !! bulk-forcing open-water weight.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      integer :: i, j, ncovered, nbad
      checks: block
         call setup(grid, state, 1)
         call ocean_state_enter_data(state)
         !$acc update self(state%ice%aice, state%ice%aice0)
         ncovered = 0
         nbad = 0
         do j = 1, size(state%ice%aice, 2)
            do i = 1, size(state%ice%aice, 1)
               if (state%ice%aice(i, j) == 1.0_wp) then
                  ncovered = ncovered + 1
                  if (state%ice%aice0(i, j) /= 1.0_wp) nbad = nbad + 1
               end if
            end do
         end do
         call check(error, ncovered > 0, "the lumped case has fully covered cells")
         if (allocated(error)) exit checks
         call check(error, nbad, 0, "aice0 stays 1 under full lumped cover")
         call ocean_state_exit_data(state)
      end block checks
      call state%destroy()
   end subroutine test_aice0_lumped

   subroutine test_engine_every_step(error)
      !! Pins the refresh point on the PRODUCTION path: ncat 5, category
      !! transport, `dt_therm_ratio = 3` (the category state changes only on
      !! every third step); after EVERY step the aggregates equal a fresh
      !! device recompute from the category state, bitwise.  A refresh that
      !! lagged the state by a step, or missed the thermo window, fails here.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t), target :: e
      type(config_t) :: cfg
      character(len=*), parameter :: NL = new_line("a")
      character(len=:), allocatable :: nml
      real(wp), allocatable :: ref_a(:, :), ref_v(:, :), ref_s(:, :)
      real(wp) :: t, dtl, amax
      integer :: n, ierr, nbad
      nml = "&sim_nml sim_type = 'ocean' /"//NL// &
            "&grid_nml nx = 20, ny = 16, dx = 10000.0, dy = 10000.0, nghost = 3 /"//NL// &
            "&time_nml t_end = 86400.0, dt_fixed = 600.0 /"//NL// &
            "&physics_nml coriolis_f = 1.435e-4 /"//NL// &
            "&nonhydrostatic_nml nz_layers = 4 /"//NL// &
            "&tracer_nml initial_salinity = 34.0, initial_temperature = -1.80 /"//NL// &
            "&ocean_ic_nml alpha_T = 1.7e-4, rho_0 = 1027.0 /"//NL// &
            "&ocean_thermo_nml enable_thermodynamics = .true., q_heat = -60.0 /"//NL// &
            "&ocean_topo_nml topo_config = 'flat', max_depth = 50.0 /"//NL// &
            "&ocean_hvisc_nml nu_h = 500.0 /"//NL// &
            "&ocean_vmix_nml kv_ml_invz2 = 0.01, hmix_fixed = 20.0, dt_therm_ratio = 3 /"//NL// &
            "&vcoord_nml vcoord_type = 'sigma' /"//NL// &
            "&ocean_bt_nml auto_n_inner = .true. /"//NL// &
            "&ocean_ice_nml enable = .true., ncat = 5, nk_ice = 2, air_temp = -20.0, "// &
            "restore_lambda = 20.0, sw_down = 0.0, transport = .true., dynamics = .false. /"//NL// &
            "&ocean_ice_ic_nml conc_config = 'uniform', h_ice = 1.0, conc = 0.7 /"//NL// &
            "&ocean_diag_nml enabled = .false. /"//NL// &
            "&output_nml output_to_file = .false. /"//NL
      checks: block
         call comm_env_init()
         call comm_env_setup_roles(.false.)
         call read_config_from_string(nml, cfg, ierr=ierr)
         call check(error, ierr == OCEAN_STATUS_OK, "config parse")
         if (allocated(error)) exit checks
         call validate_config(cfg, ierr)
         call check(error, ierr == OCEAN_STATUS_OK, "config validate")
         if (allocated(error)) exit checks
         call engine_setup(e, cfg, ierr)
         call check(error, ierr == OCEAN_STATUS_OK, "engine_setup")
         if (allocated(error)) exit checks
         call engine_enter_data(e, cfg)
         t = 0.0_wp
         dtl = 600.0_wp
         amax = 0.0_wp
         do n = 1, 9
            call engine_step(e, dtl, t, ierr=ierr)
            if (ierr == OCEAN_STATUS_OK) call engine_step_ice(e, cfg, dtl, t, ierr=ierr)
            if (ierr == OCEAN_STATUS_OK) call engine_step_finalize(e, dtl, t, ierr=ierr)
            call check(error, ierr == OCEAN_STATUS_OK, "engine step")
            if (allocated(error)) exit checks
            t = t + dtl
            !$acc update self(e%state%ice%aice, e%state%ice%vice, e%state%ice%vsno, e%state%ice%aice0)
            call agg_mismatches(e%state, 5, nbad, ref_a, ref_v, ref_s)
            call check(error, nbad, 0, "aggregates equal a fresh recompute after every step")
            if (allocated(error)) exit checks
            amax = max(amax, maxval(e%state%ice%vice))
         end do
         call check(error, amax > 0.0_wp, "the run carried live ice")
      end block checks
      call engine_exit_data(e)
      call engine_teardown(e)
   end subroutine test_engine_every_step

end module test_ocean_ice_aggregates
