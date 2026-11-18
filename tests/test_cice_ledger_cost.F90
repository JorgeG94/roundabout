!! What the named-term conservation ledger costs per ice step (`&ocean_ice_nml ledger`), measured on the production
!! engine path: the same closed basin (remap transport, EVP, ridging, thermodynamics every step) run with the ledger off
!! and on, `ice` time per step reported (`engine_step_ice` only, device-synchronised).  The test asserts only what is
!! machine independent: both runs finish, the ledger-on state is bitwise the ledger-off one (observer) in the category
!! volumes, and the ledger report is finite; the timings are printed for the PR record.
module test_cice_ledger_cost
   use, intrinsic :: iso_fortran_env, only: int64, output_unit
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_config, only: config_t, read_config_from_string
   use rdb_ocean_engine, only: ocean_engine_t, engine_setup, engine_enter_data, engine_step, engine_step_ice, &
                               engine_step_finalize, engine_exit_data, engine_teardown
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   use rdb_comm_env, only: comm_env_init, comm_env_setup_roles
   implicit none
   private

   public :: collect_cice_ledger_cost_tests

   real(wp), parameter :: DT = 600.0_wp
   integer, parameter :: NX = 192, NY = 128, NWARM = 4, NMEAS = 20
   logical :: comm_inited = .false.

contains

   subroutine collect_cice_ledger_cost_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [new_unittest("ledger_step_cost", test_ledger_step_cost)]
   end subroutine collect_cice_ledger_cost_tests

   function basin_nml(ledger) result(nml)
      logical, intent(in) :: ledger
      character(len=:), allocatable :: nml
      character(len=*), parameter :: NL = new_line("a")
      character(len=16) :: snx, sny

      write (snx, '(i0)') NX
      write (sny, '(i0)') NY
      nml = "&sim_nml sim_type = 'ocean' /"//NL// &
            "&time_nml t_end = 86400.0, dt_fixed = 600.0 /"//NL// &
            "&nonhydrostatic_nml nz_layers = 3 /"//NL// &
            "&tracer_nml initial_temperature = -1.0, initial_salinity = 34.0, "// &
            "T_init_surface = -1.0, T_init_bottom = -1.0 /"//NL// &
            "&ocean_bt_nml auto_n_inner = .true. /"//NL// &
            "&ocean_hvisc_nml nu_h = 200.0 /"//NL// &
            "&ocean_diag_nml enabled = .false. /"//NL// &
            "&grid_nml nx = "//trim(snx)//", ny = "//trim(sny)//", nghost = 3, dx = 20000.0, dy = 20000.0 /"//NL// &
            "&physics_nml coriolis_f = 1.0e-4, wind_stress_x = 0.1, wind_stress_y = 0.03 /"//NL// &
            "&vcoord_nml vcoord_type = 'sigma' /"//NL// &
            "&ocean_topo_nml topo_config = 'flat', max_depth = 100.0 /"//NL// &
            "&ocean_thermo_nml enable_thermodynamics = .true. /"//NL// &
            "&ocean_forcing_nml enable_components = .true. /"//NL// &
            "&ocean_vmix_nml dt_therm_ratio = 1 /"//NL// &
            "&ocean_ice_nml enable = .true., scheme = 'cice', ncat = 5, nilyr = 4, hfreeze = 10.0, "// &
            "dynamics = .true., transport = .true., ledger = "//merge(".true. ", ".false.", ledger)//" /"//NL// &
            "&ocean_ice_thermo_nml update_ocn_f = .true. /"//NL// &
            "&ocean_ice_dyn_nml advection = 'remap', ndte = 20 /"//NL// &
            "&ocean_ice_ic_nml conc_config = 'uniform', h_ice = 1.0, conc = 0.9, h_snow = 0.1, t_ice = -8.0 /"//NL// &
            "&ocean_ice_atm_nml uni_tas = 255.0, uni_rlds = 180.0, uni_uas = 12.0, uni_vas = 3.0 /"//NL// &
            "&ocean_bc_nml west = 'wall', east = 'wall', south = 'wall', north = 'wall' /"//NL// &
            "&output_nml output_to_file = .false. /"//NL
   end function basin_nml

   subroutine run_basin(ledger, ms_per_step, vsum, finite, ok)
      !! `NWARM` steps, then `NMEAS` timed ice steps.  `vsum`: the ice volume summed over the categories (a state checksum).
      logical, intent(in) :: ledger
      real(wp), intent(out) :: ms_per_step, vsum
      logical, intent(out) :: finite, ok
      type(ocean_engine_t), target :: engine
      type(config_t) :: cfg
      integer :: ierr, n
      integer(int64) :: c0, c1, rate
      real(wp) :: t, tice

      ok = .false.
      finite = .true.
      ms_per_step = 0.0_wp
      vsum = 0.0_wp
      if (.not. comm_inited) then
         call comm_env_init()
         call comm_env_setup_roles(.false.)
         comm_inited = .true.
      end if
      call read_config_from_string(basin_nml(ledger), cfg, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call engine_setup(engine, cfg, ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call engine_enter_data(engine, cfg)
      t = 0.0_wp
      tice = 0.0_wp
      do n = 1, NWARM + NMEAS
         call engine_step(engine, DT, t, ierr=ierr)
         if (ierr /= OCEAN_STATUS_OK) exit
         !$acc wait
         call system_clock(c0, rate)
         call engine_step_ice(engine, cfg, DT, t, ierr=ierr)
         !$acc wait
         call system_clock(c1)
         if (ierr /= OCEAN_STATUS_OK) exit
         if (n > NWARM) tice = tice + real(c1 - c0, wp)/real(rate, wp)
         call engine_step_finalize(engine, DT, t, ierr=ierr)
         if (ierr /= OCEAN_STATUS_OK) exit
         t = t + DT
      end do
      if (ierr == OCEAN_STATUS_OK) then
         ms_per_step = 1.0e3_wp*tice/real(NMEAS, wp)
         !$acc update self(engine%state%ice%cice%vicen)
         vsum = sum(engine%state%ice%cice%vicen)
         if (ledger) finite = ieee_is_finite(engine%state%icec%ledger%cum_scale(3)) .and. engine%state%icec%ledger%finite
         ok = .true.
      end if
      call engine_exit_data(engine)
      call engine_teardown(engine)
   end subroutine run_basin

   subroutine test_ledger_step_cost(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: ms_off, ms_on, v_off, v_on
      logical :: fin, ok

      call run_basin(.false., ms_off, v_off, fin, ok)
      call check(error, ok, "the ledger-off run finished")
      if (allocated(error)) return
      call run_basin(.true., ms_on, v_on, fin, ok)
      call check(error, ok, "the ledger-on run finished")
      if (allocated(error)) return
      write (output_unit, "(a,i0,a,i0,a,f9.3,a,f9.3,a,f6.1,a)") "  ice step, ", NX, " x ", NY, " cells: ledger off ", ms_off, &
         " ms, ledger on ", ms_on, " ms (", 100.0_wp*(ms_on/max(ms_off, 1.0e-12_wp) - 1.0_wp), " % more)"
      call check(error, fin, "the ledger report is finite")
      if (.not. allocated(error)) call check(error, transfer(v_on, 0_int64) == transfer(v_off, 0_int64), &
                                             "ledger on == ledger off (category volume sum, bitwise)")
   end subroutine test_ledger_step_cost

end module test_cice_ledger_cost
