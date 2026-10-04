!! A sigma-seeded LAGRANGIAN stack over a staircase shelf must stay bounded.
!!
!! The defect this pins
!! --------------------
!! `&ocean_isopycnal_nml pgf_skip_nonoverlap` (default ON under
!! `vcoord_type="lagrangian"`) zeroes the face PGF of a layer whose z-extents in
!! the two abutting columns do not overlap.  That is right for what it was
!! written for — an isopycnal layer GROUNDED against the bed, squeezed onto the
!! floor on one side and massive on the other — and wrong for a layer that is
!! MASSIVE on both sides and merely sits at different depths.  A sigma-seeded
!! stack (the default `thickness_config`) over a staircase is exactly that: at a
!! 100 m -> 250 m step a 10 m layer at 90-100 m faces a 25 m layer at
!! 225-250 m, and nine of the ten layers of every such face were gated (17 % of
!! all wet layer-faces on the compat-matrix domain).  Continuity kept moving
!! their mass across the face while the momentum equation no longer felt the
!! pressure gradient that mass flux works against, so the PGF-work / PE
!! exchange was broken: `En` grew ~200x between steps 50 and 150 and the run
!! hit the MaxCFL panic at step 178 (`pred_corr`) / 154 (`ssp_rk2`).  The same
!! geometry under `vcoord_type="sigma"` is bounded, and so is the Lagrangian
!! run once the gate requires the layer to be vanished on one side
!! (`ocean_pressure_force_t%nonoverlap_vanish_tol`).
!!
!! The case
!! --------
!! The compatibility matrix's base domain (24 x 16 x 10, dx = 20 km, f + beta,
!! staircase shelf with an along-shelf offset so steps face both ways, a 3 x 2
!! island), staged in memory through `rdb_ocean_stage_bathymetry`; analytic
!! linear T/S, wind 0.1 Pa, 60 W/m^2 cooling, KPP, Smagorinsky, `fv_mom6`.
!! The analytic IC tilts every layer's density with the bed, so the run spins up
!! a real shelf-break current in the first ~25 steps (En ~0.1 m^2/s^2, the same
!! as the sigma twin) and must then stay at that level.
!!
!! What is asserted, for BOTH split schemes, over `N_STEPS` outer steps:
!!   * every chunk steps cleanly and the state stays finite;
!!   * the layer KE at the end is at most `GROWTH_MAX` x its value at
!!     `N_SPINUP` (measured: 1.3x after the fix; on the defect it is > 100x and
!!     the CFL wall is reached).
module test_ocean_lagrangian_staircase
   use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double, c_null_ptr, c_f_pointer
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_constants, only: wp
   use rdb_ocean_api, only: rdb_ocean_create_pending, rdb_ocean_stage_bathymetry, &
                            rdb_ocean_create_finalize, rdb_ocean_step, &
                            rdb_ocean_destroy, rdb_ocean_refresh_host, &
                            rdb_ocean_get_grid_info, rdb_ocean_get_h_layer_ptr, &
                            rdb_ocean_get_wet_t_ptr, &
                            rdb_ocean_get_u_face_x_layer_ptr, &
                            rdb_ocean_get_v_face_y_layer_ptr
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   use testdrive, only: error_type, check, new_unittest, unittest_type
   implicit none
   private

   public :: collect_ocean_lagrangian_staircase_tests

   integer, parameter :: NXP = 24, NYP = 16, NZ = 10
   integer, parameter :: N_SPINUP = 40
      !! Outer steps before the reference KE sample (the IC spin-up is over).
   integer, parameter :: N_STEPS = 160
      !! Total outer steps (the defect is > 100x by step 150).
   integer, parameter :: CHUNK = 20
   real(wp), parameter :: GROWTH_MAX = 3.0_wp
   integer(c_int), parameter :: BATHY_DEPTH_POSITIVE_DOWN = 1_c_int
   real(wp), parameter :: STAIRCASE(8) = [100.0_wp, 100.0_wp, 250.0_wp, 250.0_wp, &
                                          500.0_wp, 800.0_wp, 1200.0_wp, 1600.0_wp]
   real(wp), parameter :: MAX_DEPTH = 2000.0_wp

contains

   subroutine collect_ocean_lagrangian_staircase_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("pred_corr_staircase_stays_bounded", test_pred_corr), &
                  new_unittest("ssp_rk2_staircase_stays_bounded", test_ssp_rk2) &
                  ]
   end subroutine collect_ocean_lagrangian_staircase_tests

   pure function staircase_depth() result(b)
      !! `compat_matrix.bathymetry()`: rows south (j=1) to north; columns
      !! i = 9..16 sit one row further south (an along-shelf offset, so the
      !! steps face both x and y); a 3 x 2 island at i = 11..13, j = 10..11.
      real(c_double) :: b(NXP, NYP)
      integer :: i, j, js
      do j = 1, NYP
         do i = 1, NXP
            js = j - 1
            if (i >= 9 .and. i <= 16) js = js - 1
            if (js < 0) then
               b(i, j) = STAIRCASE(1)
            else if (js < size(STAIRCASE)) then
               b(i, j) = STAIRCASE(js + 1)
            else
               b(i, j) = MAX_DEPTH
            end if
            if (i >= 11 .and. i <= 13 .and. j >= 10 .and. j <= 11) b(i, j) = 0.0_wp
         end do
      end do
   end function staircase_depth

   function staircase_nml(scheme) result(nml)
      character(len=*), intent(in) :: scheme
      character(len=:), allocatable :: nml
      character(len=1), parameter :: nl = new_line("a")
      nml = "&sim_nml sim_type = 'ocean' /"//nl// &
            "&grid_nml nx = 24, ny = 16, dx = 20000.0, dy = 20000.0, nghost = 4 /"//nl// &
            "&nonhydrostatic_nml nz_layers = 10 /"//nl// &
            "&time_nml t_end = 1.0e9, dt_fixed = 900.0 /"//nl// &
            "&physics_nml coriolis_f = 1.0e-4, wind_stress_x = 0.1, wind_stress_y = 0.0 /"//nl// &
            "&tracer_nml initial_temperature = 9.0, initial_salinity = 34.8, "// &
            "T_init_surface = 16.0, T_init_bottom = 3.0, "// &
            "S_init_surface = 34.6, S_init_bottom = 35.0 /"//nl// &
            "&ocean_topo_nml topo_config = 'flat', max_depth = 2000.0, "// &
            "wind_config = 'constant', coriolis_beta = 2.0e-11, coriolis_y_ref = 160000.0 /"//nl// &
            "&ocean_thermo_nml enable_thermodynamics = .true., q_heat = -60.0 /"//nl// &
            "&ocean_eos_nml eos = 'wright' /"//nl// &
            "&ocean_pgf_nml form = 'fv_mom6' /"//nl// &
            "&ocean_coriolis_nml form = 'sadourny_energy' /"//nl// &
            "&ocean_hvisc_nml nu_h = 200.0, c_smag = 0.15, lateral_closure = 'smagorinsky', "// &
            "smag_ah = .true. /"//nl// &
            "&ocean_bdrag_nml form = 'quadratic', cd = 0.003, bg_vel = 0.05, hbbl = 10.0 /"//nl// &
            "&ocean_vmix_nml use_closure = .true., use_kpp = .true. /"//nl// &
            "&ocean_bt_nml auto_n_inner = .true., split_scheme = '"//scheme//"' /"//nl// &
            "&vcoord_nml vcoord_type = 'lagrangian', check_vanished_content = .true. /"//nl// &
            "&ocean_diag_nml enabled = .false. /"//nl// &
            "&output_nml output_to_file = .false. /"//nl
   end function staircase_nml

   subroutine layer_ke(handle, ke, finite)
      !! Thickness-weighted mean kinetic energy (m^2/s^2) over the wet
      !! physical cells, from face velocities averaged to cell centres.
      type(c_ptr), intent(in) :: handle
      real(wp), intent(out) :: ke
      logical, intent(out) :: finite
      type(c_ptr) :: ptr
      integer(c_int) :: status, nx_p, ny_p, nz_p, ng, nx, ny, nz, gen
      real(wp), pointer :: h(:, :, :), u(:, :, :), v(:, :, :), wet(:, :)
      real(wp) :: uc, vc, num, den
      integer :: i, j, k

      status = rdb_ocean_refresh_host(handle)
      status = rdb_ocean_get_grid_info(handle, nx_p, ny_p, nz_p, ng)
      status = rdb_ocean_get_h_layer_ptr(handle, ptr, nx, ny, nz, gen)
      call c_f_pointer(ptr, h, [nx, ny, nz])
      status = rdb_ocean_get_u_face_x_layer_ptr(handle, ptr, nx, ny, nz, gen)
      call c_f_pointer(ptr, u, [nx, ny, nz])
      status = rdb_ocean_get_v_face_y_layer_ptr(handle, ptr, nx, ny, nz, gen)
      call c_f_pointer(ptr, v, [nx, ny, nz])
      status = rdb_ocean_get_wet_t_ptr(handle, ptr, nx, ny, gen)
      call c_f_pointer(ptr, wet, [nx, ny])

      num = 0.0_wp
      den = 0.0_wp
      do k = 1, nz_p
         do j = ng + 1, ng + ny_p
            do i = ng + 1, ng + nx_p
               if (wet(i, j) <= 0.5_wp) cycle
               uc = 0.5_wp*(u(i, j, k) + u(i + 1, j, k))
               vc = 0.5_wp*(v(i, j, k) + v(i, j + 1, k))
               num = num + 0.5_wp*h(i, j, k)*(uc*uc + vc*vc)
               den = den + h(i, j, k)
            end do
         end do
      end do
      ke = num/max(den, tiny(1.0_wp))
      finite = ieee_is_finite(ke) .and. ieee_is_finite(num)
   end subroutine layer_ke

   subroutine run_staircase(error, scheme)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: scheme
      type(c_ptr) :: handle
      integer(c_int) :: status
      character(len=:), allocatable :: nml
      real(c_double) :: bathy(NXP, NYP)
      real(wp) :: ke, ke_ref
      logical :: finite
      integer :: n_done
      character(len=200) :: msg

      nml = staircase_nml(scheme)
      bathy = staircase_depth()
      handle = c_null_ptr
      status = rdb_ocean_create_pending(nml, len(nml, kind=c_int), handle)
      call check(error, status == OCEAN_STATUS_OK, scheme//": create_pending")
      if (allocated(error)) return

      body: block
         status = rdb_ocean_stage_bathymetry(handle, bathy, int(NXP, c_int), int(NYP, c_int), &
                                             BATHY_DEPTH_POSITIVE_DOWN)
         call check(error, status == OCEAN_STATUS_OK, scheme//": stage_bathymetry")
         if (allocated(error)) exit body
         status = rdb_ocean_create_finalize(handle)
         call check(error, status == OCEAN_STATUS_OK, scheme//": create_finalize")
         if (allocated(error)) exit body

         status = rdb_ocean_step(handle, int(N_SPINUP, c_int))
         call check(error, status == OCEAN_STATUS_OK, scheme//": spin-up steps")
         if (allocated(error)) exit body
         call layer_ke(handle, ke_ref, finite)
         call check(error, finite .and. ke_ref > 1.0e-6_wp, &
                    scheme//": the spin-up must leave a finite, moving ocean (non-vacuity)")
         if (allocated(error)) exit body

         n_done = N_SPINUP
         do while (n_done < N_STEPS)
            status = rdb_ocean_step(handle, int(CHUNK, c_int))
            n_done = n_done + CHUNK
            call layer_ke(handle, ke, finite)
            write (msg, "(a,a,i0,a,es10.3,a,es10.3)") scheme, ": step ", n_done, &
               "  KE ", ke, "  KE(spin-up) ", ke_ref
            call check(error, status == OCEAN_STATUS_OK .and. finite, &
                       trim(msg)//" -- the run went non-finite / stopped")
            if (allocated(error)) exit body
            call check(error, ke <= GROWTH_MAX*ke_ref, &
                       trim(msg)//" -- layer KE is growing: the grounded-layer PGF "// &
                       "gate is zeroing massive non-overlapping layers again")
            if (allocated(error)) exit body
         end do
      end block body
      status = rdb_ocean_destroy(handle)
   end subroutine run_staircase

   subroutine test_pred_corr(error)
      type(error_type), allocatable, intent(out) :: error
      call run_staircase(error, "pred_corr")
   end subroutine test_pred_corr

   subroutine test_ssp_rk2(error)
      type(error_type), allocatable, intent(out) :: error
      call run_staircase(error, "ssp_rk2")
   end subroutine test_ssp_rk2

end module test_ocean_lagrangian_staircase
