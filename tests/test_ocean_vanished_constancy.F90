!! Tracer CONSTANCY through vanished layers — the gate for invariant I1′.
!!
!! A uniform tracer must stay uniform: whatever the dynamics do to the
!! thicknesses, a column that starts at one concentration everywhere has no
!! way to make another.  That is a basic property of a conservative,
!! consistent transport scheme, and it is the one the previous vanished-layer
!! rule broke.
!!
!! The previous rule (I1, `h <= H_VANISHED ⇒ hTr = 0`) kept a filler's MASS but
!! zeroed its CONTENT.  Between two remaps the continuity step moves thickness
!! out of a filler (and nudges it a hair above the marker, so it reads live);
!! that thickness carries the filler's concentration, `hTr/h = 0`, into the
!! live layers — fresh, 0 °C water appearing in a uniform-S = 35, T = 15 run
!! (`double_gyre_mom6`, 10 days: 258 thin live cells at S = 0–34.43 and 1273
!! thick cells off 35 by up to 2e-5).  Salt was conserved; concentration was
!! not consistent.
!!
!! I1′ (`h <= H_VANISHED ⇒ hTr = h·c_live`, the donor live layer's
!! concentration) keeps it.  Both tests here run the FULL split solver through
!! the C API — continuity, the ALE remap and the enforcement point, ≥ 50 outer
!! steps with the flow forced — and assert that every wet cell, filler or
!! live, still holds the initial concentration to `1e-12` relative, for T, S
!! and (where the envelope allows one) a passive tracer:
!!
!!   * `ZSTAR_FULL`, the reduced `double_gyre_mom6` (spoon bed, two-gyre wind,
!!     `dt_therm_ratio = 2`): the bed layer vanishes over the shelf;
!!   * `Z_FIXED` under an ice-shelf cavity: a sloping draft vanishes the top
!!     layers into the ice (fillers ABOVE the live column, donor `k_top`) and a
!!     seamount vanishes the bed layers, so one column carries fillers at both
!!     ends.
!!
!! Each test also asserts it is not vacuous: fillers exist in wet columns,
!! and the flow is not at rest.
module test_ocean_vanished_constancy
   use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_null_ptr, c_f_pointer
   use rdb_constants, only: wp, H_VANISHED
   use rdb_ocean_api, only: rdb_ocean_create_from_string, rdb_ocean_step, &
                            rdb_ocean_destroy, rdb_ocean_refresh_host, &
                            rdb_ocean_get_grid_info, rdb_ocean_get_h_layer_ptr, &
                            rdb_ocean_get_wet_t_ptr, rdb_ocean_get_tracer_ptr, &
                            rdb_ocean_get_u_face_x_layer_ptr, rdb_ocean_set_u
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   use testdrive, only: error_type, check, new_unittest, unittest_type
   implicit none
   private

   public :: collect_ocean_vanished_constancy_tests

   real(wp), parameter :: S0 = 35.0_wp
      !! Uniform initial salinity (psu).
   real(wp), parameter :: T0 = 15.0_wp
      !! Uniform initial temperature (degC).
   real(wp), parameter :: REL_TOL = 1.0e-12_wp
      !! Constancy tolerance, relative.
   integer, parameter :: N_STEPS = 60
      !! Outer steps per run (the brief asks for >= 50).

contains

   subroutine collect_ocean_vanished_constancy_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("zstar_full_uniform_tracer_stays_uniform", test_zstar_full), &
                  new_unittest("z_fixed_cavity_uniform_tracer_stays_uniform", test_z_fixed_cavity) &
                  ]
   end subroutine collect_ocean_vanished_constancy_tests

   function nml_zstar_full() result(nml)
      !! `double_gyre_mom6.nml`, reduced to 22 x 20 cells: spoon bed
      !! (100–2000 m), two-gyre wind, 2 layers with a 1000 m surface target,
      !! so layer 1 vanishes to the filler over the shelf.  `dt_therm_ratio =
      !! 2` is the configuration that exposed the defect: the coordinate is
      !! only restored on thermo steps, so between two remaps continuity moves
      !! thickness through the fillers.  A passive pseudo-salt rides along.
      character(len=:), allocatable :: nml
      character(len=1), parameter :: nl = new_line("a")
      nml = "&sim_nml sim_type = 'ocean' /"//nl// &
            "&grid_nml nx = 22, ny = 20, nghost = 3 /"//nl// &
            "&ocean_grid_nml axis_units = 'degrees', len_lon = 22.0, len_lat = 20.0, "// &
            "rad_earth = 6.378e6 /"//nl// &
            "&time_nml t_end = 1.0e9, dt_fixed = 1200.0, cfl_interval = 1 /"//nl// &
            "&physics_nml coriolis_f = 9.4e-5 /"//nl// &
            "&nonhydrostatic_nml nz_layers = 2 /"//nl// &
            "&tracer_nml initial_temperature = 15.0, initial_salinity = 35.0 /"//nl// &
            "&ocean_ic_nml layer_rho_init = 1036.0, 1035.0 /"//nl// &
            "&ocean_coriolis_nml form = 'sadourny_energy' /"//nl// &
            "&ocean_thermo_nml enable_thermodynamics = .false. /"//nl// &
            "&ocean_topo_nml topo_config = 'spoon', max_depth = 2000.0, edge_depth = 100.0, "// &
            "slope_scale = 400000.0, wind_config = '2gyre', taux_magnitude = 0.1, "// &
            "coriolis_beta = 1.76e-11, coriolis_y_ref = 1113200.0 /"//nl// &
            "&ocean_pgf_nml form = 'gprime', gprime_gfs = 0.98, gprime_gint = 0.0098, "// &
            "maxvel = 6.0 /"//nl// &
            "&ocean_bdrag_nml form = 'linear', r = 2.5e-5, hbbl = 10.0, bg_vel = 0.1, "// &
            "bbl_thick_min = 0.1 /"//nl// &
            "&ocean_hvisc_nml nu_h = 10000.0, smag_ah = .true., lateral_closure = 'smagorinsky' /"//nl// &
            "&ocean_vmix_nml hmix_fixed = 20.0, hmix_stress = 20.0, harmonic_visc = .true., "// &
            "dt_therm_ratio = 2 /"//nl// &
            "&ocean_bt_nml auto_n_inner = .true., cfl_bt_safety = 0.65, bebt = 0.2 /"//nl// &
            "&vcoord_nml vcoord_type = 'zstar_full', zstar_h_surf_target = 1000.0, "// &
            "zstar_h_min = 1.5e-4, check_vanished_content = .true. /"//nl// &
            "&ocean_tracers_nml enable_pseudo_salt = .true. /"//nl// &
            "&ocean_diag_nml enabled = .false. /"//nl// &
            "&output_nml output_to_file = .false. /"//nl
   end function nml_zstar_full

   function nml_z_fixed_cavity() result(nml)
      !! `Z_FIXED` under an ice-shelf cavity with fillers at BOTH ends of the
      !! column: a linear draft (250 m at the west wall, shallowing to the
      !! calving front at 36 km) vanishes the top layers into the ice, and a
      !! Gaussian seamount (720 m → 450 m) vanishes the bed layers of the
      !! 15-layer, 48 m stack.  Uniform T/S under a linear EOS is a rest
      !! state, so the flow is stirred by an initial zonal jet written through
      !! the API.  Adiabatic (the `z_fixed` x cavity v1 envelope: no KPP, no
      !! melt), so nothing but transport acts on the tracers.
      character(len=:), allocatable :: nml
      character(len=1), parameter :: nl = new_line("a")
      nml = "&sim_nml sim_type = 'ocean' /"//nl// &
            "&grid_nml nx = 24, ny = 12, nghost = 2, dx = 2000.0, dy = 2000.0 /"//nl// &
            "&time_nml t_end = 1.0e9, dt_fixed = 300.0, cfl_interval = 1 /"//nl// &
            "&physics_nml coriolis_f = -1.409e-4 /"//nl// &
            "&nonhydrostatic_nml nz_layers = 15 /"//nl// &
            "&vcoord_nml vcoord_type = 'z_fixed', zfixed_closed_faces = .true., "// &
            "check_vanished_content = .true. /"//nl// &
            "&ocean_topo_nml topo_config = 'seamount', max_depth = 720.0, edge_depth = 450.0, "// &
            "slope_scale = 12000.0 /"//nl// &
            "&ocean_cavity_dyn_nml enable = .true., draft_config = 'linear', "// &
            "draft_depth = 250.0, draft_slope = -4.4e-3, draft_x0 = -2000.0, "// &
            "draft_x1 = 36000.0, h_min_cavity = 96.0 /"//nl// &
            "&ocean_pgf_nml form = 'fv_mom6', p_top_in_bc = .true., maxvel = 2.0 /"//nl// &
            "&ocean_eos_nml eos = 'linear' /"//nl// &
            "&ocean_ic_nml alpha_T = 3.8356948e-2, beta_S = 8.0587609e-1, T_ref = -1.0, "// &
            "S_ref = 34.2, rho_0 = 1027.51 /"//nl// &
            "&tracer_nml initial_temperature = 15.0, initial_salinity = 35.0 /"//nl// &
            "&ocean_coriolis_nml form = 'sadourny' /"//nl// &
            "&ocean_hvisc_nml nu_h = 20.0 /"//nl// &
            "&ocean_bdrag_nml form = 'quadratic', cd = 2.5e-3 /"//nl// &
            "&ocean_vmix_nml use_closure = .false., use_kpp = .false. /"//nl// &
            "&ocean_bt_nml auto_n_inner = .true. /"//nl// &
            "&ocean_diag_nml enabled = .false. /"//nl// &
            "&output_nml output_to_file = .false. /"//nl
   end function nml_z_fixed_cavity

   subroutine run_constancy(error, nml, what, seed_u, tracer_names, refs)
      !! Create, (optionally) stir, step `N_STEPS`, then check every WET
      !! physical cell of every named tracer against its uniform initial
      !! value.  Fillers are checked too — under I1′ their `hTr/h` IS the
      !! donor's concentration, i.e. the same uniform value.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: nml, what
      real(wp), intent(in) :: seed_u
         !! Initial depth-uniform zonal velocity (m/s); 0 ⇒ no stirring.
      character(len=*), intent(in) :: tracer_names(:)
      real(wp), intent(in) :: refs(:)

      type(c_ptr) :: handle, ptr
      integer(c_int) :: status, nx_p, ny_p, nz_p, ng, nx, ny, nz, gen
      real(wp), pointer :: h(:, :, :), q(:, :, :), wet(:, :), u(:, :, :)
      real(wp), allocatable :: ubuf(:, :, :)
      real(wp) :: worst_live, worst_fill, umax, c
      integer :: i, j, k, t, n_fill, n_live
      character(len=256) :: msg

      handle = c_null_ptr
      status = rdb_ocean_create_from_string(nml, len(nml, kind=c_int), handle)
      call check(error, status == OCEAN_STATUS_OK, what//": the namelist must build an ocean")
      if (allocated(error)) return

      body: block
         status = rdb_ocean_get_grid_info(handle, nx_p, ny_p, nz_p, ng)
         if (seed_u /= 0.0_wp) then
            allocate (ubuf(nx_p + 1, ny_p, nz_p), source=seed_u)
            status = rdb_ocean_set_u(handle, ubuf, nx_p, ny_p, nz_p)
            call check(error, status == OCEAN_STATUS_OK, what//": seeding the jet")
            if (allocated(error)) exit body
         end if

         status = rdb_ocean_step(handle, int(N_STEPS, c_int))
         call check(error, status == OCEAN_STATUS_OK, what//": the run must step cleanly")
         if (allocated(error)) exit body
         status = rdb_ocean_refresh_host(handle)

         status = rdb_ocean_get_h_layer_ptr(handle, ptr, nx, ny, nz, gen)
         call c_f_pointer(ptr, h, [nx, ny, nz])
         status = rdb_ocean_get_wet_t_ptr(handle, ptr, nx, ny, gen)
         call c_f_pointer(ptr, wet, [nx, ny])
         status = rdb_ocean_get_u_face_x_layer_ptr(handle, ptr, nx, ny, nz, gen)
         call c_f_pointer(ptr, u, [nx, ny, nz])

         ! Non-vacuity: fillers in wet columns, and a moving ocean.
         n_fill = 0
         n_live = 0
         umax = 0.0_wp
         do k = 1, nz_p
            do j = ng + 1, ng + ny_p
               do i = ng + 1, ng + nx_p
                  if (wet(i, j) <= 0.5_wp) cycle
                  if (h(i, j, k) <= H_VANISHED) then
                     n_fill = n_fill + 1
                  else
                     n_live = n_live + 1
                  end if
                  umax = max(umax, abs(u(i, j, k)))
               end do
            end do
         end do
         call check(error, n_fill > 0, what//": no filler in a wet column — the test is vacuous")
         if (allocated(error)) exit body
         call check(error, umax > 1.0e-4_wp, what//": the ocean is at rest — the test is vacuous")
         if (allocated(error)) exit body

         do t = 1, size(tracer_names)
            status = rdb_ocean_get_tracer_ptr(handle, trim(tracer_names(t)), &
                                              int(len_trim(tracer_names(t)), c_int), &
                                              ptr, nx, ny, nz, gen)
            call check(error, status == OCEAN_STATUS_OK, &
                       what//": tracer "//trim(tracer_names(t))//" must be registered")
            if (allocated(error)) exit body
            call c_f_pointer(ptr, q, [nx, ny, nz])
            worst_live = 0.0_wp
            worst_fill = 0.0_wp
            do k = 1, nz_p
               do j = ng + 1, ng + ny_p
                  do i = ng + 1, ng + nx_p
                     if (wet(i, j) <= 0.5_wp) cycle
                     if (h(i, j, k) <= 0.0_wp) cycle
                     c = q(i, j, k)/h(i, j, k)
                     if (h(i, j, k) > H_VANISHED) then
                        worst_live = max(worst_live, abs(c/refs(t) - 1.0_wp))
                     else
                        worst_fill = max(worst_fill, abs(c/refs(t) - 1.0_wp))
                     end if
                  end do
               end do
            end do
            write (msg, "(a,a,a,i0,a,i0,a,es10.3,a,es10.3)") what//": ", trim(tracer_names(t)), &
               " live=", n_live, " fillers=", n_fill, " worst rel (live)=", worst_live, &
               " (fillers)=", worst_fill
            call check(error, worst_live <= REL_TOL, trim(msg)//" — a LIVE cell left the uniform value")
            if (allocated(error)) exit body
            call check(error, worst_fill <= REL_TOL, trim(msg)//" — a FILLER does not carry c_live")
            if (allocated(error)) exit body
         end do
      end block body

      status = rdb_ocean_destroy(handle)
   end subroutine run_constancy

   subroutine test_zstar_full(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=16) :: names(3)
      names = [character(len=16) :: "salinity", "temperature", "pseudo_salt"]
      call run_constancy(error, nml_zstar_full(), "zstar_full double-gyre", 0.0_wp, &
                         names, [S0, T0, S0])
   end subroutine test_zstar_full

   subroutine test_z_fixed_cavity(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=16) :: names(2)
      names = [character(len=16) :: "salinity", "temperature"]
      call run_constancy(error, nml_z_fixed_cavity(), "z_fixed cavity", 0.05_wp, &
                         names, [S0, T0])
   end subroutine test_z_fixed_cavity

end module test_ocean_vanished_constancy
