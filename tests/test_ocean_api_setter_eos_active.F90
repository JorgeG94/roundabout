module test_ocean_api_setter_eos_active
   !! `rdb_ocean_set_h` / `rdb_ocean_set_tracer` recompute `rho_layer` with
   !! `ocean_eos_compute`.  The time step passes
   !! `active = enable_thermodynamics`; the setters must too, or with
   !! `enable_thermodynamics = .false.` and prescribed layer densities
   !! (`&ocean_ic_nml layer_rho_init`) an API setter overwrites those
   !! densities with EOS values the step itself would never produce.
   !!
   !! Gates: thermodynamics OFF leaves `rho_layer` bit-identical across both
   !! setters; thermodynamics ON still recomputes it (unchanged behaviour).
   use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_null_ptr, c_f_pointer
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_ocean_api, only: rdb_ocean_create_from_string, rdb_ocean_destroy, &
                            rdb_ocean_get_rho_layer_ptr, rdb_ocean_get_grid_info, &
                            rdb_ocean_set_tracer, rdb_ocean_set_h
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   implicit none
   private

   public :: collect_ocean_api_setter_eos_active_tests

contains

   subroutine collect_ocean_api_setter_eos_active_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("setters_keep_prescribed_rho_thermo_off", test_thermo_off), &
                  new_unittest("setters_recompute_rho_thermo_on", test_thermo_on) &
                  ]
   end subroutine collect_ocean_api_setter_eos_active_tests

   function case_nml(thermo) result(nml)
      logical, intent(in) :: thermo
      character(len=:), allocatable :: nml
      character(len=8) :: flag

      flag = ".true."
      if (.not. thermo) flag = ".false."
      nml = "&sim_nml sim_type = 'ocean' /"//new_line("a")// &
            "&grid_nml nx = 8, ny = 6, nghost = 3, dx = 10000.0, dy = 10000.0 /"// &
            new_line("a")// &
            "&nonhydrostatic_nml nz_layers = 3 /"//new_line("a")// &
            "&time_nml t_end = 1.0e7, dt_fixed = 60.0 /"//new_line("a")// &
            "&ocean_topo_nml max_depth = 300.0, taux_magnitude = 0.0 /"//new_line("a")// &
            "&physics_nml coriolis_f = 0.0 /"//new_line("a")// &
            "&tracer_nml initial_salinity = 35.0, T_init_bottom = 10.0, "// &
            "T_init_surface = 10.0 /"//new_line("a")// &
            "&ocean_ic_nml alpha_T = 0.2, T_ref = 10.0, S_ref = 35.0, rho_0 = 1035.0, "// &
            "layer_rho_init = 1030.0, 1028.0, 1026.0 /"//new_line("a")// &
            "&ocean_thermo_nml enable_thermodynamics = "//trim(flag)//" /"//new_line("a")// &
            "&vcoord_nml vcoord_type = 'sigma' /"//new_line("a")// &
            "&ocean_bc_nml west = 'periodic', east = 'periodic' /"//new_line("a")// &
            "&ocean_bdrag_nml form = 'linear', r = 0.0 /"//new_line("a")// &
            "&ocean_bt_nml auto_n_inner = .false., n_inner = 10 /"//new_line("a")// &
            "&ocean_diag_nml enabled = .false. /"//new_line("a")// &
            "&output_nml output_to_file = .false. /"//new_line("a")
   end function case_nml

   subroutine run_case(thermo, changed_t, changed_h, ok)
      !! Snapshot `rho_layer`, call each setter, and report whether the
      !! density moved after the temperature setter and after the h setter.
      logical, intent(in) :: thermo
      logical, intent(out) :: changed_t, changed_h, ok

      type(c_ptr) :: handle, ptr
      integer(c_int) :: status, nx, ny, nz, gen, nxp, nyp, nzp, ngc
      real(wp), pointer :: r3(:, :, :)
      real(wp), allocatable :: tval(:, :, :), hval(:, :, :), r0(:, :, :)
      character(len=:), allocatable :: nml

      ok = .false.
      changed_t = .false.
      changed_h = .false.
      handle = c_null_ptr
      nml = case_nml(thermo)
      status = rdb_ocean_create_from_string(nml, len(nml, kind=c_int), handle)
      if (status /= OCEAN_STATUS_OK) return
      status = rdb_ocean_get_grid_info(handle, nxp, nyp, nzp, ngc)
      if (status /= OCEAN_STATUS_OK) then
         status = rdb_ocean_destroy(handle); return
      end if

      ! Getters hand back a host copy refreshed from the device.
      status = rdb_ocean_get_rho_layer_ptr(handle, ptr, nx, ny, nz, gen)
      if (status /= OCEAN_STATUS_OK) then
         status = rdb_ocean_destroy(handle); return
      end if
      call c_f_pointer(ptr, r3, [int(nx), int(ny), int(nz)])
      r0 = r3

      allocate (tval(int(nxp), int(nyp), int(nzp)))
      tval = 25.0_wp
      status = rdb_ocean_set_tracer(handle, "temperature", 11_c_int, tval, nxp, nyp, nzp)
      if (status /= OCEAN_STATUS_OK) then
         status = rdb_ocean_destroy(handle); return
      end if
      status = rdb_ocean_get_rho_layer_ptr(handle, ptr, nx, ny, nz, gen)
      call c_f_pointer(ptr, r3, [int(nx), int(ny), int(nz)])
      changed_t = any(r3 /= r0)

      r0 = r3
      allocate (hval(int(nxp), int(nyp), int(nzp)))
      hval = 90.0_wp
      status = rdb_ocean_set_h(handle, hval, nxp, nyp, nzp)
      if (status /= OCEAN_STATUS_OK) then
         status = rdb_ocean_destroy(handle); return
      end if
      status = rdb_ocean_get_rho_layer_ptr(handle, ptr, nx, ny, nz, gen)
      call c_f_pointer(ptr, r3, [int(nx), int(ny), int(nz)])
      changed_h = any(r3 /= r0)

      status = rdb_ocean_destroy(handle)
      ok = (status == OCEAN_STATUS_OK)
   end subroutine run_case

   subroutine test_thermo_off(error)
      type(error_type), allocatable, intent(out) :: error
      logical :: ct, ch, ok

      call run_case(.false., ct, ch, ok)
      call check(error, ok, "model must build and the setters must succeed")
      if (allocated(error)) return
      call check(error,.not. ct, &
                 "set_tracer must leave prescribed rho_layer bit-identical when "// &
                 "enable_thermodynamics = .false.")
      if (allocated(error)) return
      call check(error,.not. ch, &
                 "set_h must leave prescribed rho_layer bit-identical when "// &
                 "enable_thermodynamics = .false.")
   end subroutine test_thermo_off

   subroutine test_thermo_on(error)
      type(error_type), allocatable, intent(out) :: error
      logical :: ct, ch, ok

      call run_case(.true., ct, ch, ok)
      call check(error, ok, "model must build and the setters must succeed")
      if (allocated(error)) return
      call check(error, ct, "set_tracer must still recompute rho_layer with thermodynamics on")
   end subroutine test_thermo_on

end module test_ocean_api_setter_eos_active
