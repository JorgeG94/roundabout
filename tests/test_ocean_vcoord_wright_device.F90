!! RHO / HYCOM density coordinates under the WRIGHT equation of state,
!! driven through the device data path.
!!
!! The isopycnal target builder (`ocean_vcoord_compute_target_h_rho`)
!! evaluates `eos_density_point(eos, T, S, p)` per layer inside its
!! column `do concurrent`.  The linear EOS never reads `p`; Wright does.
!! On the nvfortran GPU build (`-stdpar=gpu -gpu=mem:separate`) the
!! reference pressure used to reach the device callee as the HOST address
!! of `vcoord%rho_ref_pressure` (an `associate`-name over the scalar
!! component, passed by reference to the out-of-module `!$acc routine
!! seq` function), so every Wright read of `p` was an illegal address —
!! `CUDA_ERROR_ILLEGAL_ADDRESS` at the first regrid, even on a flat bed at
!! rest.  The sibling suites (`test_ocean_vcoord_rho`, `_hycom`) only use
!! the linear EOS, which is why they never saw it.
!!
!! Cases:
!!   1. rho_wright_target_on_device / 2. hycom_wright_target_on_device —
!!      the target builder alone on a two-layer warm-over-cold column
!!      whose interior target density lies INSIDE the density jump at the
!!      coordinate's reference pressure (2000 dbar) but DENSER than the
!!      whole column at the surface.  The analytic answer is "the new
!!      interface lands on the jump" (target_h = h_old); a reference
!!      pressure read as ~0 instead sends it to the bed, so the case
!!      catches a silently wrong `p` as well as a faulting one.
!!   3. rho_wright_engine_steps / 4. hycom_wright_engine_steps — a tiny
!!      flat-bed, salinity-stratified rest-state engine (the vcoord
!!      stability-matrix cell, shrunk) created from an in-memory namelist
!!      and stepped through several ALE regrids: layers stay finite and
!!      positive and mass is conserved.
!!
!! Mapping mirrors production: parent-first `copyin` of the vcoord
!! aggregate (production gets it from the root `copyin(state)`), then the
!! slot's own `enter_data` for the component arrays.  All directives are
!! inert on host builds, so the suite runs everywhere.
module test_ocean_vcoord_wright_device
   use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double, c_null_ptr, c_f_pointer
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_constants, only: wp, VCOORD_RHO, VCOORD_HYCOM
   use rdb_grid, only: hgrid_t
   use rdb_eos, only: eos_t, eos_density_point, EOS_VARIANT_WRIGHT_97
   use rdb_ocean_vcoord, only: ocean_vcoord_t
   use rdb_ocean_api, only: rdb_ocean_create_from_string, rdb_ocean_step, &
                            rdb_ocean_destroy, rdb_ocean_refresh_host, &
                            rdb_ocean_get_h_layer_ptr, rdb_ocean_get_total_mass
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   use testdrive, only: error_type, check, new_unittest, unittest_type
   implicit none
   private

   public :: collect_ocean_vcoord_wright_device_tests

   integer, parameter :: NZ2 = 2
      !! Layers in the target-builder column.
   real(wp), parameter :: H_TOP = 70.0_wp
      !! Warm upper layer thickness (m), state index k = 2.
   real(wp), parameter :: H_BED = 30.0_wp
      !! Cold bed layer thickness (m), state index k = 1.  Unequal to
      !! H_TOP so the HYCOM z* floor (at dsig·H = 50 m) is not binding and
      !! cannot mask a mis-placed interface.
   real(wp), parameter :: T_WARM = 14.0_wp, T_COLD = 2.0_wp, S_UNI = 35.0_wp
   real(wp), parameter :: P_REF = 2.0e7_wp
      !! Coordinate reference pressure (Pa) — sigma-2, the shipped default.

contains

   subroutine collect_ocean_vcoord_wright_device_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("rho_wright_target_on_device", test_rho_target), &
                  new_unittest("hycom_wright_target_on_device", test_hycom_target), &
                  new_unittest("rho_wright_engine_steps", test_rho_engine), &
                  new_unittest("hycom_wright_engine_steps", test_hycom_engine) &
                  ]
   end subroutine collect_ocean_vcoord_wright_device_tests

   subroutine test_rho_target(error)
      type(error_type), allocatable, intent(out) :: error
      call run_target_case(error, VCOORD_RHO, "rho")
   end subroutine test_rho_target

   subroutine test_hycom_target(error)
      type(error_type), allocatable, intent(out) :: error
      call run_target_case(error, VCOORD_HYCOM, "hycom")
   end subroutine test_hycom_target

   subroutine run_target_case(error, coord_type, label)
      !! Two-layer jump column, Wright EOS, target builder on the device.
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: coord_type
      character(len=*), intent(in) :: label
      type(hgrid_t) :: grid
      type(ocean_vcoord_t) :: vc
      type(eos_t) :: eos
      real(wp), allocatable :: eta(:, :)
      real(wp) :: rho_warm, rho_cold, rho_cold_surf, tol, err_max
      integer :: nx, ny

      call grid%init(3, 3, 1, 1.0_wp, 1.0_wp)
      nx = grid%nx_total
      ny = grid%ny_total

      eos%variant = EOS_VARIANT_WRIGHT_97
      eos%is_init = .true.

      call vc%init(grid, nz_ml=NZ2)
      vc%coord_type = coord_type
      vc%rho_ref_pressure = P_REF

      ! Bottom-up state: k = 1 cold bed layer, k = 2 warm upper layer.
      vc%remap_h_old(:, :, 1) = H_BED
      vc%remap_h_old(:, :, 2) = H_TOP
      vc%remap_conc_t(:, :, 1) = T_COLD
      vc%remap_conc_t(:, :, 2) = T_WARM
      vc%remap_conc_s(:, :, :) = S_UNI
      vc%remap_h_ref(:, :) = H_BED + H_TOP
      allocate (eta(nx, ny), source=0.0_wp)

      ! Interior target = the midpoint of the jump AT the reference
      ! pressure; the outer two bracket it (only 1..nz-1 are inverted).
      rho_warm = eos_density_point(eos, T_WARM, S_UNI, P_REF)
      rho_cold = eos_density_point(eos, T_COLD, S_UNI, P_REF)
      rho_cold_surf = eos_density_point(eos, T_COLD, S_UNI, 0.0_wp)
      vc%rho_target(0) = rho_warm - 1.0_wp
      vc%rho_target(1) = 0.5_wp*(rho_warm + rho_cold)
      vc%rho_target(2) = rho_cold + 1.0_wp

      checks: block
         ! Premise: the target is denser than the whole column at p = 0,
         ! so a reference pressure lost on the way to the device moves
         ! the interface to the bed instead of leaving it on the jump.
         call check(error, vc%rho_target(1) > rho_cold_surf, &
                    label//": premise — p matters for this target")
         if (allocated(error)) exit checks

         !$acc enter data copyin(vc)
         call vc%enter_data()
         !$acc enter data copyin(eta)
         call vc%compute_target_h_rho(vc%remap_h_ref, eta, vc%remap_conc_t, &
                                      vc%remap_conc_s, eos, &
                                      hybrid=(coord_type == VCOORD_HYCOM))
         !$acc update self(vc%target_h)
         !$acc exit data delete(eta)
         call vc%exit_data()
         !$acc exit data delete(vc)

         tol = 1.0e-9_wp*(H_BED + H_TOP)
         err_max = max(maxval(abs(vc%target_h(:, :, 1) - H_BED)), &
                       maxval(abs(vc%target_h(:, :, 2) - H_TOP)))
         call check(error, err_max <= tol, &
                    label//": interface lands on the density jump at p_ref")
      end block checks
      deallocate (eta)
      call vc%destroy()
   end subroutine run_target_case

   subroutine test_rho_engine(error)
      type(error_type), allocatable, intent(out) :: error
      call run_engine_case(error, "rho")
   end subroutine test_rho_engine

   subroutine test_hycom_engine(error)
      type(error_type), allocatable, intent(out) :: error
      call run_engine_case(error, "hycom")
   end subroutine test_hycom_engine

   function engine_nml(vcoord) result(txt)
      !! The vcoord stability-matrix rest-state cell (flat bed, salinity
      !! stratified 33.8 -> 34.55 PSU surface -> bed, Wright EOS, fv_mom6
      !! PGF) shrunk to 6 x 4 x 4 so it runs in a second on any toolchain.
      !! The stratification comes from the `&tracer_nml` linear-in-layer IC,
      !! not `&ocean_zinit_nml` (which needs a NetCDF build), so the suite
      !! also runs on the NetCDF-off legs.
      character(len=*), intent(in) :: vcoord
      character(len=:), allocatable :: txt
      character(len=*), parameter :: nl = new_line("a")
      txt = '&sim_nml sim_type = "ocean" /'//nl// &
            "&grid_nml nx = 6, ny = 4, dx = 2000.0, dy = 2000.0, nghost = 2 /"//nl// &
            "&time_nml t_end = 86400.0, dt_fixed = 600.0 /"//nl// &
            "&physics_nml coriolis_f = -1.409e-4 /"//nl// &
            "&nonhydrostatic_nml nz_layers = 4 /"//nl// &
            '&vcoord_nml vcoord_type = "'//vcoord//'", zstar_h_min = 1.0e-4, '// &
            "rho_target_light = 1027.20, rho_target_dense = 1027.90, "// &
            "rho_ref_pressure = 0.0 /"//nl// &
            '&ocean_topo_nml topo_config = "flat", max_depth = 1000.0 /'//nl// &
            '&ocean_pgf_nml form = "fv_mom6", reconstruct_for_pressure = .true. /'//nl// &
            '&ocean_eos_nml eos = "wright" /'//nl// &
            "&tracer_nml initial_temperature = -1.9, initial_salinity = 33.8, "// &
            "S_init_surface = 33.8, S_init_bottom = 34.55 /"//nl// &
            "&ocean_vmix_nml use_closure = .false., use_kpp = .false. /"//nl// &
            "&ocean_bt_nml auto_n_inner = .true. /"//nl// &
            "&ocean_diag_nml enabled = .false. /"//nl// &
            "&output_nml output_to_file = .false. /"//nl
   end function engine_nml

   subroutine run_engine_case(error, vcoord)
      !! Create -> step through several ALE regrids -> check -> destroy.
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: vcoord
      integer(c_int), parameter :: N_STEPS = 4_c_int
      type(c_ptr) :: handle, ptr
      character(len=:), allocatable :: nml
      integer(c_int) :: status, nx, ny, nz, gen
      real(c_double) :: mass0, mass1
      real(wp), pointer :: h(:, :, :)
      logical :: ok

      handle = c_null_ptr
      nml = engine_nml(vcoord)
      status = rdb_ocean_create_from_string(nml, len(nml, kind=c_int), handle)
      checks: block
         call check(error, status == OCEAN_STATUS_OK, vcoord//": create succeeds")
         if (allocated(error)) exit checks
         status = rdb_ocean_get_total_mass(handle, mass0)
         call check(error, status == OCEAN_STATUS_OK, vcoord//": initial mass")
         if (allocated(error)) exit checks

         status = rdb_ocean_step(handle, N_STEPS)
         call check(error, status == OCEAN_STATUS_OK, vcoord//": steps succeed")
         if (allocated(error)) exit checks

         status = rdb_ocean_get_total_mass(handle, mass1)
         call check(error, status == OCEAN_STATUS_OK, vcoord//": final mass")
         if (allocated(error)) exit checks
         call check(error, ieee_is_finite(real(mass1, wp)) .and. &
                    abs(mass1 - mass0) <= 1.0e-12_wp*abs(mass0), &
                    vcoord//": mass conserved through the regrids")
         if (allocated(error)) exit checks

         status = rdb_ocean_refresh_host(handle)
         if (status == OCEAN_STATUS_OK) &
            status = rdb_ocean_get_h_layer_ptr(handle, ptr, nx, ny, nz, gen)
         call check(error, status == OCEAN_STATUS_OK, vcoord//": h_layer readable")
         if (allocated(error)) exit checks
         call c_f_pointer(ptr, h, [int(nx), int(ny), int(nz)])
         ok = all(ieee_is_finite(h)) .and. all(h > 0.0_wp)
         call check(error, ok, vcoord//": layer thicknesses finite and positive")
      end block checks
      status = rdb_ocean_destroy(handle)
   end subroutine run_engine_case

end module test_ocean_vcoord_wright_device
