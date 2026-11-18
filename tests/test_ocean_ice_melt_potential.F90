!! Tests of the melt-potential provider `ice_melt_potential`
!! (`rdb_ice_melt_potential`): the signed `frzmlt` [W m-2] the sea ice reads.
!!
!! Closed forms use the SAME ocean freezing point the frazil kernel clamps
!! to (`eos_freezing_point(eos, S, 0)`), so the tests pin the integral, the
!! depth cap, the sign rule and the guards, not the liquidus.
!!
!! Cases:
!!   * uniform warm column, `hfreeze` cutting a layer partially; `hfreeze`
!!     deeper than the column (stops at the bed);
!!   * a vanished top layer is skipped (integration starts on the first live
!!     layer, the partial last layer is taken from the live ones);
!!   * `h_bl` < `hfreeze` wins, `h_bl` > `hfreeze` and `h_bl <= 0` fall back;
!!   * sign exclusivity: bank > 0 wins over warm water, no bank + cold water
!!     is exactly +0, warm is negative;
!!   * land and ghosts; non-finite input gives NaN in both outputs, never 0,
!!     and a NaN below the integrated depth is not read;
!!   * `hfreeze` validation in `validate_config`;
!!   * every cell of a random state compared with an independent host
!!     reference, with every array device-mapped (`mem:separate`).
module test_ocean_ice_melt_potential
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, RHO_WATER, H_VANISHED
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_eos, only: eos_t, eos_freezing_point
   use rdb_ice_state, only: ocean_sea_ice_t
   use rdb_ice_melt_potential, only: ice_melt_potential
   use rdb_config, only: config_t, read_config_from_string, validate_config
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_CONFIG_VALIDATE
   use rdb_ocean_engine, only: ocean_engine_t, engine_setup, engine_enter_data, &
                               engine_step, engine_step_ice, engine_step_finalize, &
                               engine_exit_data, engine_teardown
   use rdb_comm_env, only: comm_env_init, comm_env_setup_roles
   implicit none
   private

   public :: collect_ocean_ice_melt_potential_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 4
   logical, save :: comm_inited = .false.
   real(wp), parameter :: H0 = 10.0_wp
   real(wp), parameter :: S0 = 35.0_wp
   real(wp), parameter :: DT = 900.0_wp

contains

   subroutine collect_ocean_ice_melt_potential_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("melt_potential_partial_layer", test_partial), &
                  new_unittest("melt_potential_deeper_than_column", test_deep), &
                  new_unittest("melt_potential_vanished_top_skipped", test_vanished), &
                  new_unittest("melt_potential_h_bl_cap", test_hbl), &
                  new_unittest("melt_potential_bed_fillers_and_strict_live", test_fillers), &
                  new_unittest("melt_potential_slot_off_by_default", test_slot_off), &
                  new_unittest("melt_potential_sign_exclusivity", test_sign), &
                  new_unittest("melt_potential_land_and_ghosts", test_land_ghost), &
                  new_unittest("melt_potential_nan_in", test_nan), &
                  new_unittest("melt_potential_validate_hfreeze", test_validate), &
                  new_unittest("melt_potential_engine_cadence_and_h_bl", test_engine), &
                  new_unittest("melt_potential_device_vs_host", test_device_vs_host) &
                  ]
   end subroutine collect_ocean_ice_melt_potential_tests

   subroutine setup(grid, ms, eos, ice, nx, ny, nz_in)
      type(hgrid_t), intent(out) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(eos_t), intent(out) :: eos
      type(ocean_sea_ice_t), intent(inout) :: ice
      integer, intent(in) :: nx, ny, nz_in
      call grid%init(nx, ny, NGHOST, 1.0_wp, 1.0_wp)
      ms%nz_ml = nz_in
      call ms%init(grid)
      call eos%init(grid)
      ice%enable = .true.
      ice%melt_potential_enable = .true.
      call ice%init(grid)
      ms%h_layer = H0
      ms%tracers(ms%idx_salinity)%c = S0
      ms%tracers(ms%idx_temperature)%c = 2.0_wp
   end subroutine setup

   subroutine teardown(ms, eos, ice)
      type(multilayer_state_t), intent(inout) :: ms
      type(eos_t), intent(inout) :: eos
      type(ocean_sea_ice_t), intent(inout) :: ice
      call ice%destroy()
      call eos%destroy()
      call ms%destroy()
   end subroutine teardown

   subroutine run(grid, eos, ms, ice, hfreeze, hbl)
      !! Map everything, run the provider, pull the outputs back.
      type(hgrid_t), intent(in) :: grid
      type(eos_t), intent(in) :: eos
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_sea_ice_t), intent(inout) :: ice
      real(wp), intent(in) :: hfreeze
      real(wp), intent(inout), optional :: hbl(:, :)
      !$acc enter data copyin(ms)
      call ms%enter_data()
      call ice%enter_data()
      if (present(hbl)) then
         !$acc enter data copyin(hbl)
         call ice_melt_potential(grid, eos, ms, ice%frazil_heat, hfreeze, DT, &
                                 ice%melt_potential, ice%frzmlt, h_bl=hbl)
         !$acc exit data delete(hbl)
      else
         call ice_melt_potential(grid, eos, ms, ice%frazil_heat, hfreeze, DT, &
                                 ice%melt_potential, ice%frzmlt)
      end if
      associate (mp => ice%melt_potential, fm => ice%frzmlt)
         !$acc update self(mp, fm)
      end associate
      call ice%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms)
   end subroutine run

   subroutine test_partial(error)
      !! 4 x 10 m at T = 2: hfreeze = 25 takes two layers and half of the
      !! third, MP = rho Cp 25 (T - Tf), frzmlt = -MP/dt.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      real(wp) :: tf, mp
      integer :: ip, jp
      call setup(grid, ms, eos, ice, 6, 4, NZ)
      tf = eos_freezing_point(eos, S0, 0.0_wp)
      ! Distinct T per layer so a wrong cut is visible: surface-first 3, 2, 1, 0.
      ms%tracers(ms%idx_temperature)%c(:, :, 4) = 3.0_wp
      ms%tracers(ms%idx_temperature)%c(:, :, 3) = 2.0_wp
      ms%tracers(ms%idx_temperature)%c(:, :, 2) = 1.0_wp
      ms%tracers(ms%idx_temperature)%c(:, :, 1) = 0.0_wp
      call run(grid, eos, ms, ice, 25.0_wp)
      ip = grid%nx_total/2
      jp = grid%ny_total/2
      mp = RHO_WATER*eos%c_p*(10.0_wp*(3.0_wp - tf) + 10.0_wp*(2.0_wp - tf) + 5.0_wp*(1.0_wp - tf))
      call check(error, abs(ice%melt_potential(ip, jp) - mp) <= 1.0e-12_wp*abs(mp), &
                 "partial last layer: melt potential closed form")
      if (allocated(error)) goto 100
      call check(error, abs(ice%frzmlt(ip, jp) + mp/DT) <= 1.0e-12_wp*abs(mp/DT), &
                 "frzmlt = -MP/dt")
      if (allocated(error)) goto 100
      call check(error, ice%frzmlt(ip, jp) < 0.0_wp, "melt potential is negative frzmlt")
100   call teardown(ms, eos, ice)
   end subroutine test_partial

   subroutine test_deep(error)
      !! hfreeze deeper than the 40 m column: the sum stops at the bed.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      real(wp) :: tf, mp
      integer :: ip, jp
      call setup(grid, ms, eos, ice, 6, 4, NZ)
      tf = eos_freezing_point(eos, S0, 0.0_wp)
      call run(grid, eos, ms, ice, 500.0_wp)
      ip = grid%nx_total/2
      jp = grid%ny_total/2
      mp = RHO_WATER*eos%c_p*40.0_wp*(2.0_wp - tf)
      call check(error, abs(ice%melt_potential(ip, jp) - mp) <= 1.0e-12_wp*abs(mp), &
                 "deeper than the column: whole column, not 500 m")
      call teardown(ms, eos, ice)
   end subroutine test_deep

   subroutine test_vanished(error)
      !! Top layer collapsed, hfreeze = 15: first live layer (k = 3, T = 3)
      !! in full, then 5 m of k = 2 (T = 1); the empty layer's T (99) is
      !! never read.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      real(wp) :: tf, mp
      integer :: ip, jp
      call setup(grid, ms, eos, ice, 6, 4, NZ)
      tf = eos_freezing_point(eos, S0, 0.0_wp)
      ms%h_layer(:, :, 4) = 0.0_wp
      ms%tracers(ms%idx_temperature)%c(:, :, 4) = 99.0_wp
      ms%tracers(ms%idx_temperature)%c(:, :, 3) = 3.0_wp
      ms%tracers(ms%idx_temperature)%c(:, :, 2) = 1.0_wp
      call run(grid, eos, ms, ice, 15.0_wp)
      ip = grid%nx_total/2
      jp = grid%ny_total/2
      mp = RHO_WATER*eos%c_p*(10.0_wp*(3.0_wp - tf) + 5.0_wp*(1.0_wp - tf))
      call check(error, abs(ice%melt_potential(ip, jp) - mp) <= 1.0e-12_wp*abs(mp), &
                 "vanished top layer skipped; exactly hfreeze metres integrated")
      call teardown(ms, eos, ice)
   end subroutine test_vanished

   subroutine test_fillers(error)
      !! z_fixed-style column: live layers on top, bed fillers at exactly
      !! H_VANISHED below (strict `>`: not live) and a layer a hair above it
      !! (live).  hfreeze reaches past the live water: the fillers add nothing
      !! (their T is a poison value) and the sum stops at the last live layer.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      real(wp) :: tf, mp
      integer :: ip, jp
      call setup(grid, ms, eos, ice, 6, 4, NZ)
      tf = eos_freezing_point(eos, S0, 0.0_wp)
      ms%h_layer(:, :, 4) = 10.0_wp
      ms%h_layer(:, :, 3) = 10.0_wp
      ms%h_layer(:, :, 2) = H_VANISHED
      ms%h_layer(:, :, 1) = H_VANISHED
      ms%tracers(ms%idx_temperature)%c(:, :, 2) = 500.0_wp
      ms%tracers(ms%idx_temperature)%c(:, :, 1) = 500.0_wp
      call run(grid, eos, ms, ice, 100.0_wp)
      ip = grid%nx_total/2
      jp = grid%ny_total/2
      mp = RHO_WATER*eos%c_p*20.0_wp*(2.0_wp - tf)
      call check(error, abs(ice%melt_potential(ip, jp) - mp) <= 1.0e-12_wp*abs(mp), &
                 "layers exactly at H_VANISHED are not live: only the 20 m of live water counts")
      if (allocated(error)) then
         call teardown(ms, eos, ice)
         return
      end if
      ! One ulp above H_VANISHED is live again (strict `>`).
      ms%h_layer(:, :, 2) = nearest(H_VANISHED, 1.0_wp)
      call run(grid, eos, ms, ice, 100.0_wp)
      call check(error, ice%melt_potential(ip, jp) > mp*(1.0_wp + 1.0e-4_wp), &
                 "a layer just above H_VANISHED is live (its T = 500 is read)")
      call teardown(ms, eos, ice)
   end subroutine test_fillers

   subroutine test_slot_off(error)
      !! Knob off: the slot does not allocate the two arrays.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_sea_ice_t) :: ice
      call grid%init(6, 4, NGHOST, 1.0_wp, 1.0_wp)
      ice%enable = .true.
      call ice%init(grid)
      call check(error,.not. allocated(ice%melt_potential) .and. .not. allocated(ice%frzmlt), &
                 "melt_potential / frzmlt are not allocated unless hfreeze > 0")
      call ice%destroy()
   end subroutine test_slot_off

   subroutine test_hbl(error)
      !! h_bl < hfreeze wins; h_bl > hfreeze and h_bl <= 0 give hfreeze.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      real(wp), allocatable :: hbl(:, :)
      real(wp) :: tf, rc
      integer :: ip, jp
      call setup(grid, ms, eos, ice, 6, 4, NZ)
      tf = eos_freezing_point(eos, S0, 0.0_wp)
      rc = RHO_WATER*eos%c_p
      ms%tracers(ms%idx_temperature)%c(:, :, 4) = 3.0_wp
      ms%tracers(ms%idx_temperature)%c(:, :, 3) = 2.0_wp
      allocate (hbl(grid%nx_total, grid%ny_total))
      ip = grid%nx_total/2
      jp = grid%ny_total/2
      checks: block
         hbl = 6.0_wp
         call run(grid, eos, ms, ice, 25.0_wp, hbl)
         call check(error, abs(ice%melt_potential(ip, jp) - rc*6.0_wp*(3.0_wp - tf)) &
                    <= 1.0e-12_wp*rc*6.0_wp*(3.0_wp - tf), "h_bl = 6 < hfreeze = 25")
         if (allocated(error)) exit checks
         hbl = 60.0_wp
         call run(grid, eos, ms, ice, 12.0_wp, hbl)
         call check(error, abs(ice%melt_potential(ip, jp) &
                               - rc*(10.0_wp*(3.0_wp - tf) + 2.0_wp*(2.0_wp - tf))) &
                    <= 1.0e-12_wp*rc*30.0_wp, "h_bl = 60 > hfreeze = 12 -> 12 m")
         if (allocated(error)) exit checks
         hbl = 0.0_wp
         call run(grid, eos, ms, ice, 12.0_wp, hbl)
         call check(error, ice%melt_potential(ip, jp) == 0.0_wp .and. ice%frzmlt(ip, jp) == 0.0_wp, &
                    "h_bl = 0: zero depth gives zero melt potential, NOT hfreeze")
         if (allocated(error)) exit checks
         hbl = -3.0_wp
         call run(grid, eos, ms, ice, 12.0_wp, hbl)
         call check(error, ice%melt_potential(ip, jp) == 0.0_wp, "h_bl < 0 clamps to zero depth")
         if (allocated(error)) exit checks
         ! A zero depth still passes a positive bank through.
         ice%frazil_heat(ip, jp) = 9.0e5_wp
         call run(grid, eos, ms, ice, 12.0_wp, hbl)
         call check(error, ice%frzmlt(ip, jp) == 9.0e5_wp/DT, "zero depth + bank: bank/dt")
      end block checks
      call teardown(ms, eos, ice)
   end subroutine test_hbl

   subroutine test_sign(error)
      !! Bank > 0 wins over warm water below; no bank + cold water is +0
      !! (not negative, not positive); warm water gives a negative frzmlt.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      real(wp) :: tf
      integer :: i0, j0
      call setup(grid, ms, eos, ice, 6, 4, NZ)
      tf = eos_freezing_point(eos, S0, 0.0_wp)
      i0 = NGHOST + 1
      j0 = NGHOST + 1
      ! cell (i0, j0): warm + bank; (i0+1, j0): cold, no bank; (i0+2, j0): warm
      ms%tracers(ms%idx_temperature)%c(i0 + 1, j0, :) = tf - 0.5_wp
      ice%frazil_heat(i0, j0) = 4.5e5_wp
      call run(grid, eos, ms, ice, 20.0_wp)
      checks: block
         call check(error, ice%frzmlt(i0, j0) == 4.5e5_wp/DT, "bank > 0: frzmlt = bank/dt")
         if (allocated(error)) exit checks
         call check(error, ice%melt_potential(i0, j0) > 0.0_wp, "warm water still diagnosed")
         if (allocated(error)) exit checks
         call check(error, ice%frzmlt(i0 + 1, j0) == 0.0_wp .and. &
                    .not. (ice%frzmlt(i0 + 1, j0) < 0.0_wp) .and. &
                    ice%melt_potential(i0 + 1, j0) < 0.0_wp, &
                    "cold water, no bank: frzmlt is 0, melt potential negative")
         if (allocated(error)) exit checks
         call check(error, sign(1.0_wp, ice%frzmlt(i0 + 1, j0)) > 0.0_wp, "cold-water zero is +0")
         if (allocated(error)) exit checks
         call check(error, ice%frzmlt(i0 + 2, j0) < 0.0_wp, "warm water, no bank: negative")
      end block checks
      call teardown(ms, eos, ice)
   end subroutine test_sign

   subroutine test_land_ghost(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      integer :: i0, j0
      call setup(grid, ms, eos, ice, 6, 4, NZ)
      i0 = NGHOST + 2
      j0 = NGHOST + 1
      ms%wet_mask(i0, j0) = 0.0_wp
      ice%frazil_heat(i0, j0) = 1.0e6_wp
      ice%melt_potential = -77.0_wp
      ice%frzmlt = -77.0_wp
      call run(grid, eos, ms, ice, 20.0_wp)
      checks: block
         call check(error, ice%melt_potential(i0, j0) == 0.0_wp .and. ice%frzmlt(i0, j0) == 0.0_wp, &
                    "land cell: both outputs zero even with a bank and warm water")
         if (allocated(error)) exit checks
         call check(error, ice%frzmlt(1, 1) == -77.0_wp .and. ice%melt_potential(1, 1) == -77.0_wp, &
                    "ghost cells untouched")
      end block checks
      call teardown(ms, eos, ice)
   end subroutine test_land_ghost

   subroutine test_nan(error)
      !! Non-finite T, S, h (inside the integrated depth), bank or h_bl gives
      !! NaN in BOTH outputs; a NaN below the integrated depth is not read.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      real(wp), allocatable :: hbl(:, :)
      real(wp) :: nan
      integer :: i0, j0
      nan = ieee_value(1.0_wp, ieee_quiet_nan)
      call setup(grid, ms, eos, ice, 8, 4, NZ)
      allocate (hbl(grid%nx_total, grid%ny_total))
      hbl = 50.0_wp
      i0 = NGHOST + 1
      j0 = NGHOST + 1
      ms%tracers(ms%idx_temperature)%c(i0, j0, 4) = nan
      ms%tracers(ms%idx_salinity)%c(i0 + 1, j0, 3) = nan
      ms%h_layer(i0 + 2, j0, 4) = nan
      ice%frazil_heat(i0 + 3, j0) = nan
      hbl(i0 + 4, j0) = nan
      ms%tracers(ms%idx_temperature)%c(i0, j0 + 1, 1) = nan   ! below the 20 m taken
      call run(grid, eos, ms, ice, 20.0_wp, hbl)
      checks: block
         call check(error,.not. ieee_is_finite(ice%frzmlt(i0, j0)) .and. &
                    .not. ieee_is_finite(ice%melt_potential(i0, j0)), "NaN T -> NaN")
         if (allocated(error)) exit checks
         call check(error,.not. ieee_is_finite(ice%frzmlt(i0 + 1, j0)) .and. &
                    .not. ieee_is_finite(ice%melt_potential(i0 + 1, j0)), "NaN S -> NaN")
         if (allocated(error)) exit checks
         call check(error,.not. ieee_is_finite(ice%frzmlt(i0 + 2, j0)) .and. &
                    .not. ieee_is_finite(ice%melt_potential(i0 + 2, j0)), "NaN h -> NaN")
         if (allocated(error)) exit checks
         call check(error,.not. ieee_is_finite(ice%frzmlt(i0 + 3, j0)) .and. &
                    .not. ieee_is_finite(ice%melt_potential(i0 + 3, j0)), "NaN bank -> NaN")
         if (allocated(error)) exit checks
         call check(error,.not. ieee_is_finite(ice%frzmlt(i0 + 4, j0)) .and. &
                    .not. ieee_is_finite(ice%melt_potential(i0 + 4, j0)), "NaN h_bl -> NaN")
         if (allocated(error)) exit checks
         call check(error, ieee_is_finite(ice%frzmlt(i0, j0 + 1)) .and. &
                    ieee_is_finite(ice%melt_potential(i0, j0 + 1)), &
                    "NaN below the integrated depth is not read")
      end block checks
      call teardown(ms, eos, ice)
   end subroutine test_nan

   subroutine test_validate(error)
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      integer :: ierr
      checks: block
         call parse_case(cfg, "&ocean_ice_nml enable = .true., hfreeze = 10.0 /")
         ierr = -999
         call validate_config(cfg, ierr=ierr)
         call check(error, ierr == OCEAN_STATUS_OK, "hfreeze = 10 with ice on accepted")
         if (allocated(error)) exit checks
         call parse_case(cfg, "&ocean_ice_nml enable = .true., hfreeze = 0.0 /")
         ierr = -999
         call validate_config(cfg, ierr=ierr)
         call check(error, ierr == OCEAN_STATUS_ERR_CONFIG_VALIDATE, "hfreeze = 0 refused")
         if (allocated(error)) exit checks
         call parse_case(cfg, "&ocean_ice_nml enable = .true., hfreeze = -5.0 /")
         ierr = -999
         call validate_config(cfg, ierr=ierr)
         call check(error, ierr == OCEAN_STATUS_ERR_CONFIG_VALIDATE, "hfreeze = -5 refused")
         if (allocated(error)) exit checks
         call parse_case(cfg, "&ocean_ice_nml enable = .false., hfreeze = 10.0 /")
         ierr = -999
         call validate_config(cfg, ierr=ierr)
         call check(error, ierr == OCEAN_STATUS_ERR_CONFIG_VALIDATE, "hfreeze > 0 without ice refused")
         if (allocated(error)) exit checks
         call parse_case(cfg, "&ocean_ice_nml enable = .true., hfreeze = 10.0 /"//new_line("a")// &
                         "&ocean_cavity_dyn_nml enable = .true. /")
         ierr = -999
         call validate_config(cfg, ierr=ierr)
         call check(error, ierr == OCEAN_STATUS_ERR_CONFIG_VALIDATE, &
                    "hfreeze > 0 with an ice-shelf cavity refused")
         if (allocated(error)) exit checks
         call parse_case(cfg, "&ocean_ice_nml enable = .true. /")
         call check(error, cfg%ocean%ice%hfreeze == -1.0_wp, "default hfreeze is -1 (off)")
      end block checks
   end subroutine test_validate

   subroutine parse_case(cfg, extra)
      type(config_t), intent(out) :: cfg
      character(len=*), intent(in) :: extra
      character(len=:), allocatable :: nml
      nml = '&sim_nml sim_type = "ocean" /'//new_line("a")// &
            "&grid_nml nx = 8, ny = 8, dx = 2000.0, dy = 2000.0 /"//new_line("a")// &
            "&nonhydrostatic_nml nz_layers = 3 /"//new_line("a")// &
            "&time_nml t_end = 3600.0, dt_fixed = 300.0 /"//new_line("a")// &
            "&ocean_forcing_nml enable_components = .true. /"//new_line("a")// &
            "&ocean_bt_nml n_inner = 8 /"//new_line("a")// &
            extra//new_line("a")
      call read_config_from_string(nml, cfg)
   end subroutine parse_case

   subroutine test_device_vs_host(error)
      !! Random state where every cell differs; every output of every owned
      !! cell compared with an independent host reference (interface-depth
      !! clipping formulation, not the kernel's remaining-depth loop).
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 40, NY = 30, NZR = 12
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      type(ocean_sea_ice_t) :: ice
      real(wp), allocatable :: hbl(:, :)
      real(wp) :: r, mp, fm, ztop, zbot, dsel, tol, wf, rc, tf, hfr
      integer :: i, j, k, n, nbad, nneg, npos, nzero
      integer(8) :: seed
      call setup(grid, ms, eos, ice, NX, NY, NZR)
      allocate (hbl(grid%nx_total, grid%ny_total))
      rc = RHO_WATER*eos%c_p
      hfr = 30.0_wp
      seed = 123456789_8
      do j = 1, grid%ny_total
         do i = 1, grid%nx_total
            do k = 1, NZR
               call rnd(seed, r)
               ms%h_layer(i, j, k) = merge(0.0_wp, 0.5_wp + 9.0_wp*r, r < 0.15_wp)
               call rnd(seed, r)
               ms%tracers(ms%idx_temperature)%c(i, j, k) = -2.5_wp + 8.0_wp*r
               call rnd(seed, r)
               ms%tracers(ms%idx_salinity)%c(i, j, k) = 28.0_wp + 9.0_wp*r
            end do
            call rnd(seed, r)
            ice%frazil_heat(i, j) = merge(1.0e5_wp*r, 0.0_wp, r > 0.6_wp)
            call rnd(seed, r)
            hbl(i, j) = merge(0.0_wp, 80.0_wp*r, r < 0.1_wp)
            call rnd(seed, r)
            ms%wet_mask(i, j) = merge(0.0_wp, 1.0_wp, r < 0.1_wp)
         end do
      end do
      call run(grid, eos, ms, ice, hfr, hbl)

      tol = 1.0e-12_wp
      nbad = 0
      nneg = 0
      npos = 0
      nzero = 0
      do j = NGHOST + 1, grid%ny_total - NGHOST
         do i = NGHOST + 1, grid%nx_total - NGHOST
            mp = 0.0_wp
            fm = 0.0_wp
            if (ms%wet_mask(i, j) > 0.5_wp) then
               dsel = min(hfr, max(hbl(i, j), 0.0_wp))
               ztop = 0.0_wp
               do k = NZR, 1, -1
                  if (ms%h_layer(i, j, k) <= 1.5e-4_wp) cycle
                  zbot = ztop + ms%h_layer(i, j, k)
                  tf = eos_freezing_point(eos, ms%tracers(ms%idx_salinity)%c(i, j, k), 0.0_wp)
                  wf = max(0.0_wp, min(zbot, dsel) - min(ztop, dsel))
                  mp = mp + wf*(ms%tracers(ms%idx_temperature)%c(i, j, k) - tf)
                  ztop = zbot
               end do
               mp = rc*mp
               if (ice%frazil_heat(i, j) > 0.0_wp) then
                  fm = ice%frazil_heat(i, j)/DT
               else if (mp > 0.0_wp) then
                  fm = -mp/DT
               end if
            end if
            if (.not. ieee_is_finite(ice%melt_potential(i, j)) .or. &
                .not. ieee_is_finite(ice%frzmlt(i, j))) then
               nbad = nbad + 1
               cycle
            end if
            if (abs(ice%melt_potential(i, j) - mp) > tol*(1.0e4_wp + abs(mp)) .or. &
                abs(ice%frzmlt(i, j) - fm) > tol*(1.0_wp + abs(fm))) nbad = nbad + 1
            if (fm < 0.0_wp) nneg = nneg + 1
            if (fm > 0.0_wp) npos = npos + 1
            if (fm == 0.0_wp) nzero = nzero + 1
         end do
      end do
      n = nbad
      call check(error, n == 0, "every owned cell matches the host reference")
      if (allocated(error)) goto 200
      call check(error, nneg > 50 .and. npos > 50 .and. nzero > 20, &
                 "random state exercises all three branches (counts "// &
                 itoa(nneg)//" "//itoa(npos)//" "//itoa(nzero)//")")
200   call teardown(ms, eos, ice)
   end subroutine test_device_vs_host

   subroutine test_engine(error)
      !! Engine-level: dt_therm_ratio = 3.  The provider runs only on the
      !! thermo step and divides by the thermo window (3 dt): a seeded bank B
      !! gives frzmlt = B/(3 dt) there and nothing is written on the steps in
      !! between.  Then, with no bank and warm water, the depth is the fixed
      !! hfreeze under PP81 only, the KPP bl_depth under KPP, the EPBL mld
      !! under EPBL.
      type(error_type), allocatable, intent(out) :: error
      integer :: kind
      character(len=5), parameter :: NAMES(3) = ["pp81 ", "kpp  ", "epbl "]
      do kind = 1, 3
         call engine_case(kind, NAMES(kind), error)
         if (allocated(error)) return
      end do
   end subroutine test_engine

   subroutine engine_case(kind, name, error)
      integer, intent(in) :: kind
      character(len=*), intent(in) :: name
      type(error_type), allocatable, intent(inout) :: error
      real(wp), parameter :: DTE = 600.0_wp
      real(wp), parameter :: BANK = 2.7e6_wp
      real(wp), parameter :: HF = 1000.0_wp
      character(len=1), parameter :: NL = new_line("a")
      character(len=:), allocatable :: nml, vm
      type(config_t) :: cfg
      type(ocean_engine_t) :: engine
      integer :: ierr, istep, ip, jp, ng, k
      logical :: changed(6)
      real(wp) :: tf, t_top, s_top, depth, expect, got, colh, rem, dh

      select case (kind)
      case (1)
         vm = "&ocean_vmix_nml use_kpp = .false., dt_therm_ratio = 3 /"//NL
      case (2)
         vm = "&ocean_vmix_nml use_kpp = .true., dt_therm_ratio = 3 /"//NL
      case default
         vm = "&ocean_vmix_nml use_kpp = .false., dt_therm_ratio = 3 /"//NL// &
              "&ocean_epbl_nml enable = .true. /"//NL
      end select
      nml = "&sim_nml sim_type = 'ocean' /"//NL// &
            "&time_nml t_end = 86400.0, dt_fixed = 600.0 /"//NL// &
            "&nonhydrostatic_nml nz_layers = 4 /"//NL// &
            "&tracer_nml initial_salinity = 34.0, T_init_surface = 3.0, T_init_bottom = 1.0 /"//NL// &
            "&grid_nml nx = 10, ny = 8, nghost = 3, dx = 10000.0, dy = 10000.0 /"//NL// &
            "&physics_nml coriolis_f = 1.435e-4 /"//NL// &
            "&ocean_topo_nml topo_config = 'flat', max_depth = 100.0 /"//NL// &
            "&ocean_thermo_nml enable_thermodynamics = .true., q_heat = 0.0 /"//NL// &
            "&ocean_bt_nml auto_n_inner = .true. /"//NL// &
            "&ocean_diag_nml enabled = .false. /"//NL// &
            "&ocean_ice_nml enable = .true., hfreeze = 1000.0, air_temp = 2.0, "// &
            "restore_lambda = 0.0, sw_down = 0.0 /"//NL// &
            vm
      if (.not. comm_inited) then
         call comm_env_init()
         call comm_env_setup_roles(.false.)
         comm_inited = .true.
      end if
      call read_config_from_string(nml, cfg, ierr=ierr)
      call check(error, ierr == 0, name//": namelist parses")
      if (allocated(error)) return
      call validate_config(cfg, ierr)
      call check(error, ierr == 0, name//": namelist validates")
      if (allocated(error)) return
      call engine_setup(engine, cfg, ierr)
      call check(error, ierr == OCEAN_STATUS_OK, name//": engine_setup")
      if (allocated(error)) return
      call engine_enter_data(engine, cfg)
      ng = engine%grid%nghost
      ip = ng + engine%grid%nx_total/2 - ng/2
      jp = ng + 2
      associate (fz => engine%state%ice%frazil_heat, fm => engine%state%ice%frzmlt)
         fz = BANK
         fm = -77.0_wp
         !$acc update device(fz, fm)
      end associate

      changed = .false.
      do istep = 1, 3
         call engine_step(engine, DTE, real(istep - 1, wp)*DTE, ierr=ierr)
         if (ierr == OCEAN_STATUS_OK) call engine_step_ice(engine, cfg, DTE, real(istep - 1, wp)*DTE, ierr=ierr)
         if (ierr == OCEAN_STATUS_OK) call engine_step_finalize(engine, DTE, real(istep - 1, wp)*DTE, ierr=ierr)
         call check(error, ierr == OCEAN_STATUS_OK, name//": outer step")
         if (allocated(error)) return
         associate (fm => engine%state%ice%frzmlt)
            !$acc update self(fm)
         end associate
         changed(istep) = engine%state%ice%frzmlt(ip, jp) /= -77.0_wp
      end do
      call check(error,.not. changed(1) .and. .not. changed(2) .and. changed(3), &
                 name//": provider runs on the thermo step only (step 3 of 3)")
      if (allocated(error)) return
      got = engine%state%ice%frzmlt(ip, jp)
      expect = BANK/(3.0_wp*DTE)
      call check(error, abs(got - expect) <= 1.0e-12_wp*expect, &
                 name//": frzmlt = bank / thermo window")
      if (allocated(error)) return

      ! Bank spent by the uptake; three more steps, warm water, no bank.
      do istep = 4, 6
         call engine_step(engine, DTE, real(istep - 1, wp)*DTE, ierr=ierr)
         if (ierr == OCEAN_STATUS_OK) call engine_step_ice(engine, cfg, DTE, real(istep - 1, wp)*DTE, ierr=ierr)
         if (ierr == OCEAN_STATUS_OK) call engine_step_finalize(engine, DTE, real(istep - 1, wp)*DTE, ierr=ierr)
         call check(error, ierr == OCEAN_STATUS_OK, name//": outer step (window 2)")
         if (allocated(error)) return
      end do
      associate (fm => engine%state%ice%frzmlt, mp => engine%state%ice%melt_potential, &
                 hl => engine%state%multilayer%h_layer, bl => engine%state%vmix%bl_depth, &
                 ml => engine%state%epbl%mld)
         !$acc update self(fm, mp, hl)
         if (kind == 2) then
            !$acc update self(bl)
         end if
         if (kind == 3) then
            !$acc update self(ml)
         end if
      end associate
      do k = 1, size(engine%state%multilayer%tracers)
         !$acc update self(engine%state%multilayer%tracers(k)%c)
      end do
      t_top = engine%state%multilayer%tracers(engine%state%multilayer%idx_temperature)%c(ip, jp, 4)
      s_top = engine%state%multilayer%tracers(engine%state%multilayer%idx_salinity)%c(ip, jp, 4)
      tf = eos_freezing_point(engine%state%eos, s_top, 0.0_wp)
      select case (kind)
      case (1)
         depth = HF
      case (2)
         depth = min(HF, max(engine%state%vmix%bl_depth(ip, jp), 0.0_wp))
      case default
         depth = min(HF, max(engine%state%epbl%mld(ip, jp), 0.0_wp))
      end select
      colh = sum(engine%state%multilayer%h_layer(ip, jp, :))
      depth = min(depth, colh)
      expect = 0.0_wp
      rem = depth
      do k = 4, 1, -1
         dh = min(engine%state%multilayer%h_layer(ip, jp, k), rem)
         expect = expect + dh*(engine%state%multilayer%tracers( &
                               engine%state%multilayer%idx_temperature)%c(ip, jp, k) &
                               - eos_freezing_point(engine%state%eos, &
                                                    engine%state%multilayer%tracers( &
                                                    engine%state%multilayer%idx_salinity)%c(ip, jp, k), 0.0_wp))
         rem = rem - dh
      end do
      expect = RHO_WATER*engine%state%eos%c_p*expect
      got = engine%state%ice%melt_potential(ip, jp)
      call check(error, abs(got - expect) <= 1.0e-9_wp*expect, &
                 name//": melt potential uses the active boundary-layer depth "// &
                 itoa_r(got)//itoa_r(expect)//itoa_r(depth)//itoa_r(t_top))
      if (allocated(error)) return
      if (kind /= 1) then
         call check(error, depth < colh - 1.0_wp, name//": the BL depth is the binding one (non-vacuous) "// &
                    itoa_r(depth)//itoa_r(colh))
         if (allocated(error)) return
      end if
      call check(error, abs(engine%state%ice%frzmlt(ip, jp) + got/(3.0_wp*DTE)) &
                 <= 1.0e-12_wp*abs(got/(3.0_wp*DTE)), name//": frzmlt = -MP/window")
      call engine_exit_data(engine)
      call engine_teardown(engine)
   end subroutine engine_case

   function itoa_r(x) result(s)
      real(wp), intent(in) :: x
      character(len=24) :: s
      write (s, "(es24.10)") x
   end function itoa_r

   function itoa(n) result(s)
      integer, intent(in) :: n
      character(len=12) :: s
      write (s, "(i12)") n
   end function itoa

   subroutine rnd(seed, r)
      !! Park-Miller minimal-standard generator, uniform in (0,1).
      integer(8), intent(inout) :: seed
      real(wp), intent(out) :: r
      seed = mod(seed*16807_8, 2147483647_8)
      r = real(seed, wp)/2147483647.0_wp
   end subroutine rnd

end module test_ocean_ice_melt_potential
