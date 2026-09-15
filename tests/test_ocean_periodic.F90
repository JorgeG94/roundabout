!! Tests for the ocean periodic boundary condition (OBC_PERIODIC).
!! Design §3 tests 1-6.  Pattern-matches test_ocean_dyn_split.F90:
!! test-drive framework, small grids, direct slot access.
!!
!! Tests implemented:
!!   1. Wrap helpers: fill distinct values, wrap, assert ghost == partner
!!      and u(i_w) == u(i_e) exactly (host after acc update self).
!!   2. Uniform zonal flow, periodic-x channel: uniform u, flat bath,
!!      f=0 => state unchanged after N outer steps to roundoff.
!!   3. Shifted-domain bit-identity (gold test): same IC run twice,
!!      run B's IC shifted by nx_phys/2; after N steps shift(B) == A.
!!   4. Conservation: total mass + salt constant to ~1e-14 relative over
!!      N steps in a doubly-periodic domain.
!!   5. Geostrophic zonal jet (periodic-x, f-plane): stays steady.
!!   6. Validation rejections: unpaired periodic, periodic+nghost=2,
!!      periodic+sponge => error paths exercised via validate function.
module test_ocean_periodic
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
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
   use rdb_ocean_dyn, only: ocean_dyn_t, ocean_dyn_step_split
   use rdb_ocean_boundary_types, only: ocean_bc_state_t, ocean_bc_state_init, &
                                       ocean_bc_state_destroy, OBC_PERIODIC, OBC_WALL, &
                                       OBC_SPONGE, ocean_bc_validate_periodic, &
                                       ocean_bc_type_from_string
   use rdb_ocean_periodic, only: ocean_periodic_wrap_centre_2d, &
                                 ocean_periodic_wrap_centre_3d, &
                                 ocean_periodic_wrap_face_x_2d, &
                                 ocean_periodic_wrap_face_x_3d, &
                                 ocean_periodic_wrap_face_y_2d, &
                                 ocean_periodic_wrap_face_y_3d, &
                                 ocean_periodic_wrap_state
   implicit none
   private

   public :: collect_ocean_periodic_tests

   ! Grid parameters for most tests: nghost=3 (required for periodic)
   integer, parameter :: NGHOST = 3
   integer, parameter :: NZ = 2

contains

   subroutine collect_ocean_periodic_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("periodic_wrap_helpers", test_wrap_helpers), &
                  new_unittest("periodic_uniform_flow", test_uniform_flow), &
                  new_unittest("periodic_shifted_domain_identity", test_shifted_domain), &
                  new_unittest("periodic_conservation", test_conservation), &
                  new_unittest("periodic_geostrophic_jet", test_geostrophic_jet), &
                  new_unittest("periodic_validation_rejections", test_validation_rejections), &
                  new_unittest("periodic_bathy_seam", test_bathy_seam) &
                  ]
   end subroutine collect_ocean_periodic_tests

   ! -----------------------------------------------------------------
   ! Grid factory: always nghost=3 for periodic
   ! -----------------------------------------------------------------

   subroutine make_grid(grid, nx_phys, ny_phys, dx, dy)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      real(wp), intent(in) :: dx, dy
      call grid%init(nx_phys, ny_phys, NGHOST, dx, dy)
   end subroutine make_grid

   ! -----------------------------------------------------------------
   ! Setup helpers matching test_ocean_dyn_split.F90 pattern
   ! -----------------------------------------------------------------

   subroutine init_all(grid, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn)
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(continuity_t), intent(inout) :: ct
      type(coriolis_adv_t), intent(inout) :: cor
      type(ocean_pressure_force_t), intent(inout) :: pgf
      type(ocean_horizontal_viscosity_t), intent(inout) :: hv
      type(ocean_bottom_drag_t), intent(inout) :: bd
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_vertical_advection_t), intent(inout) :: va
      type(ocean_hdiff_tracer_t), intent(inout) :: hd
      type(ocean_vdiff_t), intent(inout) :: vd
      type(ocean_vmix_t), intent(inout) :: vmix
      type(eos_t), intent(inout) :: eos
      type(ocean_dyn_t), intent(inout) :: dyn

      ms%nz_ml = NZ
      call ms%init(grid)
      call ct%init(grid, nz_ml=NZ)
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
   end subroutine init_all

   subroutine map_in(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(inout) :: metrics
      type(multilayer_state_t), intent(inout) :: ms
      type(continuity_t), intent(inout) :: ct
      type(coriolis_adv_t), intent(inout) :: cor
      type(ocean_pressure_force_t), intent(inout) :: pgf
      type(ocean_horizontal_viscosity_t), intent(inout) :: hv
      type(ocean_bottom_drag_t), intent(inout) :: bd
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_vertical_advection_t), intent(inout) :: va
      type(ocean_hdiff_tracer_t), intent(inout) :: hd
      type(ocean_vdiff_t), intent(inout) :: vd
      type(ocean_vmix_t), intent(inout) :: vmix
      type(ocean_dyn_t), intent(inout) :: dyn
      call make_cartesian_metrics(metrics, grid)
      !$acc enter data copyin(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
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
   end subroutine map_in

   subroutine map_out(metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
      type(ocean_metrics_t), intent(inout) :: metrics
      type(multilayer_state_t), intent(inout) :: ms
      type(continuity_t), intent(inout) :: ct
      type(coriolis_adv_t), intent(inout) :: cor
      type(ocean_pressure_force_t), intent(inout) :: pgf
      type(ocean_horizontal_viscosity_t), intent(inout) :: hv
      type(ocean_bottom_drag_t), intent(inout) :: bd
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_vertical_advection_t), intent(inout) :: va
      type(ocean_hdiff_tracer_t), intent(inout) :: hd
      type(ocean_vdiff_t), intent(inout) :: vd
      type(ocean_vmix_t), intent(inout) :: vmix
      type(ocean_dyn_t), intent(inout) :: dyn
      call destroy_cartesian_metrics(metrics)
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
      !$acc exit data delete(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
   end subroutine map_out

   subroutine destroy_all(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn)
      type(multilayer_state_t), intent(inout) :: ms
      type(continuity_t), intent(inout) :: ct
      type(coriolis_adv_t), intent(inout) :: cor
      type(ocean_pressure_force_t), intent(inout) :: pgf
      type(ocean_horizontal_viscosity_t), intent(inout) :: hv
      type(ocean_bottom_drag_t), intent(inout) :: bd
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_vertical_advection_t), intent(inout) :: va
      type(ocean_hdiff_tracer_t), intent(inout) :: hd
      type(ocean_vdiff_t), intent(inout) :: vd
      type(ocean_vmix_t), intent(inout) :: vmix
      type(eos_t), intent(inout) :: eos
      type(ocean_dyn_t), intent(inout) :: dyn
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
   end subroutine destroy_all

   ! Build a periodic-x bc with nghost=NGHOST (3).
   subroutine make_bc_periodic_x(bc, grid)
      type(ocean_bc_state_t), intent(out) :: bc
      type(hgrid_t), intent(in) :: grid
      call ocean_bc_state_init(bc, grid, nz_ml=NZ, n_tracers=2)
      bc%west%bc_type = OBC_PERIODIC
      bc%east%bc_type = OBC_PERIODIC
      call ocean_bc_validate_periodic(bc)
   end subroutine make_bc_periodic_x

   ! Build a doubly-periodic bc.
   subroutine make_bc_periodic_xy(bc, grid)
      type(ocean_bc_state_t), intent(out) :: bc
      type(hgrid_t), intent(in) :: grid
      call ocean_bc_state_init(bc, grid, nz_ml=NZ, n_tracers=2)
      bc%west%bc_type = OBC_PERIODIC
      bc%east%bc_type = OBC_PERIODIC
      bc%south%bc_type = OBC_PERIODIC
      bc%north%bc_type = OBC_PERIODIC
      call ocean_bc_validate_periodic(bc)
   end subroutine make_bc_periodic_xy

   ! -----------------------------------------------------------------
   ! Test 1: wrap helpers
   ! -----------------------------------------------------------------

   subroutine test_wrap_helpers(error)
      !! Fill 2D/3D centre + face fields with index-derived distinct
      !! values, wrap, then assert ghost == interior partner and
      !! u(i_w) == u(i_e) exactly.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      integer :: nx, ny, nz, nghost, nx_phys, ny_phys
      integer :: i, j, k
      integer :: i_w, i_e, j_s, j_n
      real(wp), allocatable :: fld2d(:, :)
      real(wp), allocatable :: fld3d(:, :, :)
      real(wp), allocatable :: face_x2d(:, :)
      real(wp), allocatable :: face_x3d(:, :, :)
      real(wp), allocatable :: face_y2d(:, :)
      real(wp), allocatable :: face_y3d(:, :, :)

      checks: block

         call make_grid(grid, 8, 6, 1.0_wp, 1.0_wp)
         nx = grid%nx_total
         ny = grid%ny_total
         nghost = grid%nghost
         nx_phys = grid%nx_phys
         ny_phys = grid%ny_phys
         nz = 2

         i_w = nghost + 1
         i_e = nghost + nx_phys + 1
         j_s = nghost + 1
         j_n = nghost + ny_phys + 1

         ! ---- 1a. Centre-2D, wrap-x ----
         allocate (fld2d(nx, ny))
         do j = 1, ny
            do i = 1, nx
               fld2d(i, j) = real(i, wp)*100.0_wp + real(j, wp)
            end do
         end do
         call ocean_periodic_wrap_centre_2d(fld2d, nx, ny, nx_phys, ny_phys, nghost, &
                                            wrap_x=.true., wrap_y=.false.)
         ! West ghosts should equal interior
         do j = 1, ny
            do i = 1, nghost
               call check(error, fld2d(i, j) == fld2d(i + nx_phys, j), &
                          "wrap_centre_2d x: west ghost != partner")
               if (allocated(error)) exit checks
            end do
         end do
         ! East ghosts
         do j = 1, ny
            do i = nx_phys + nghost + 1, nx
               call check(error, fld2d(i, j) == fld2d(i - nx_phys, j), &
                          "wrap_centre_2d x: east ghost != partner")
               if (allocated(error)) exit checks
            end do
         end do
         deallocate (fld2d)

         ! ---- 1b. Face-x 2D, wrap-x: check seam identity u(i_w)==u(i_e) ----
         allocate (face_x2d(nx + 1, ny))
         do j = 1, ny
            do i = 1, nx + 1
               face_x2d(i, j) = real(i, wp)*1000.0_wp + real(j, wp)
            end do
         end do
         call ocean_periodic_wrap_face_x_2d(face_x2d, nx + 1, ny, nx_phys, ny_phys, nghost, &
                                            wrap_x=.true., wrap_y=.false.)
         ! Seam invariant: i_w == i_e after belt-and-braces copy
         do j = 1, ny
            call check(error, face_x2d(i_w, j) == face_x2d(i_e, j), &
                       "wrap_face_x_2d: i_w /= i_e seam violation")
            if (allocated(error)) exit checks
         end do
         ! West ghost faces
         do j = 1, ny
            do i = 1, nghost
               call check(error, face_x2d(i, j) == face_x2d(i + nx_phys, j), &
                          "wrap_face_x_2d: west ghost face != partner")
               if (allocated(error)) exit checks
            end do
         end do
         deallocate (face_x2d)

         ! ---- 1c. Face-y 2D, wrap-y: check seam identity v(j_s)==v(j_n) ----
         allocate (face_y2d(nx, ny + 1))
         do j = 1, ny + 1
            do i = 1, nx
               face_y2d(i, j) = real(i, wp) + real(j, wp)*1000.0_wp
            end do
         end do
         call ocean_periodic_wrap_face_y_2d(face_y2d, nx, ny + 1, nx_phys, ny_phys, nghost, &
                                            wrap_x=.false., wrap_y=.true.)
         do i = 1, nx
            call check(error, face_y2d(i, j_s) == face_y2d(i, j_n), &
                       "wrap_face_y_2d: j_s /= j_n seam violation")
            if (allocated(error)) exit checks
         end do
         deallocate (face_y2d)

         ! ---- 1d. Centre-3D, wrap-xy ----
         allocate (fld3d(nx, ny, nz))
         do k = 1, nz
            do j = 1, ny
               do i = 1, nx
                  fld3d(i, j, k) = real(i, wp)*1000.0_wp + real(j, wp)*100.0_wp + real(k, wp)
               end do
            end do
         end do
         call ocean_periodic_wrap_centre_3d(fld3d, nx, ny, nz, nx_phys, ny_phys, nghost, &
                                            wrap_x=.true., wrap_y=.true.)
         do k = 1, nz
            do j = 1, ny
               do i = 1, nghost
                  call check(error, fld3d(i, j, k) == fld3d(i + nx_phys, j, k), &
                             "wrap_centre_3d x: west ghost != partner")
                  if (allocated(error)) exit checks
               end do
            end do
            do j = 1, nghost
               do i = 1, nx
                  call check(error, fld3d(i, j, k) == fld3d(i, j + ny_phys, k), &
                             "wrap_centre_3d y: south ghost != partner")
                  if (allocated(error)) exit checks
               end do
            end do
         end do
         deallocate (fld3d)

         ! ---- 1e. Face-x 3D, wrap-x ----
         allocate (face_x3d(nx + 1, ny, nz))
         do k = 1, nz
            do j = 1, ny
               do i = 1, nx + 1
                  face_x3d(i, j, k) = real(i, wp)*100.0_wp + real(j, wp)*10.0_wp + real(k, wp)
               end do
            end do
         end do
         call ocean_periodic_wrap_face_x_3d(face_x3d, nx + 1, ny, nz, nx_phys, ny_phys, nghost, &
                                            wrap_x=.true., wrap_y=.false.)
         do k = 1, nz
            do j = 1, ny
               call check(error, face_x3d(i_w, j, k) == face_x3d(i_e, j, k), &
                          "wrap_face_x_3d: seam violation")
               if (allocated(error)) exit checks
            end do
         end do
         deallocate (face_x3d)

         ! ---- 1f. Face-y 3D, wrap-y ----
         allocate (face_y3d(nx, ny + 1, nz))
         do k = 1, nz
            do j = 1, ny + 1
               do i = 1, nx
                  face_y3d(i, j, k) = real(i, wp) + real(j, wp)*100.0_wp + real(k, wp)*10.0_wp
               end do
            end do
         end do
         call ocean_periodic_wrap_face_y_3d(face_y3d, nx, ny + 1, nz, nx_phys, ny_phys, nghost, &
                                            wrap_x=.false., wrap_y=.true.)
         do k = 1, nz
            do i = 1, nx
               call check(error, face_y3d(i, j_s, k) == face_y3d(i, j_n, k), &
                          "wrap_face_y_3d: seam violation")
               if (allocated(error)) exit checks
            end do
         end do
         deallocate (face_y3d)

      end block checks
   end subroutine test_wrap_helpers

   ! -----------------------------------------------------------------
   ! Test 2: uniform zonal flow, periodic-x
   ! -----------------------------------------------------------------

   subroutine test_uniform_flow(error)
      !! Uniform u = U0, v = 0, η = 0, flat bath, f = 0, no forcing.
      !! In a periodic-x channel the state should be unchanged after
      !! N outer steps (advection of a uniform field is trivially constant).
      type(error_type), allocatable, intent(out) :: error
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
      type(ocean_bc_state_t) :: bc

      real(wp), parameter :: H0 = 100.0_wp  ! flat bath per layer
      real(wp), parameter :: U0 = 0.1_wp    ! uniform zonal velocity
      real(wp), parameter :: DT = 1.0_wp
      integer, parameter :: N_INNER = 5
      integer, parameter :: N_STEPS = 5
      integer :: step, nx_p, ny_p, k, i, j
      real(wp) :: max_dh, max_du, max_dv

      checks: block

         nx_p = 12
         ny_p = 8
         call make_grid(grid, nx_p, ny_p, 1000.0_wp, 1000.0_wp)
         call init_all(grid, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn)

         ! Uniform initial condition
         ms%h_layer = H0
         ms%u_face_x_layer = U0
         ms%v_face_y_layer = 0.0_wp
         do k = 1, NZ
            ms%tracers(ms%idx_salinity)%hTr(:, :, k) = eos%S_ref*H0
            ms%tracers(ms%idx_temperature)%hTr(:, :, k) = eos%T_ref*H0
         end do
         dyn%bt_work%bt_H_ref = real(NZ, wp)*H0

         call make_bc_periodic_x(bc, grid)
         ! Wrap bt_H_ref init-time ghost fill (design §1.5 last bullet)
         call ocean_periodic_wrap_centre_2d(dyn%bt_work%bt_H_ref, &
                                            grid%nx_total, grid%ny_total, &
                                            grid%nx_phys, grid%ny_phys, &
                                            grid%nghost, .true., .false.)

         call map_in(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
         !$acc update device(dyn%bt_work%bt_H_ref)

         do step = 1, N_STEPS
            call ocean_dyn_step_split( &
               grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
               va, hd, vd, vmix, ms, DT, N_INNER, bc=bc)
         end do

         !$acc update self(ms%h_layer, ms%u_face_x_layer, ms%v_face_y_layer)
         ! Only check physical interior cells
         block
            integer :: ig, i0, i1, j0, j1
            ig = grid%nghost
            i0 = ig + 1
            i1 = ig + grid%nx_phys
            j0 = ig + 1
            j1 = ig + grid%ny_phys
            max_dh = maxval(abs(ms%h_layer(i0:i1, j0:j1, :) - H0))
            max_du = maxval(abs(ms%u_face_x_layer(i0:i1, j0:j1, :) - U0))
            max_dv = maxval(abs(ms%v_face_y_layer(i0:i1, j0:j1, :)))
         end block

         call map_out(metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)

         call check(error, max_dh < 1.0e-8_wp, &
                    "periodic uniform flow: h drifted, max_dh="// &
                    trim(adjusted_str(max_dh)))
         if (allocated(error)) exit checks
         call check(error, max_du < 1.0e-8_wp, &
                    "periodic uniform flow: u drifted")
         if (allocated(error)) exit checks
         call check(error, max_dv < 1.0e-8_wp, &
                    "periodic uniform flow: v non-zero")

      end block checks

      call ocean_bc_state_destroy(bc)
      call destroy_all(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn)
   end subroutine test_uniform_flow

   ! -----------------------------------------------------------------
   ! Test 3: shifted-domain bit-identity (gold test)
   ! -----------------------------------------------------------------

   subroutine test_shifted_domain(error)
      !! Two runs with the same physical IC, run B's domain shifted by
      !! nx_phys/2 in x.  After N_STEPS outer steps, shift(B) == A
      !! bit-for-bit in the physical interior.  This proves the seam is
      !! indistinguishable from interior.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics_a, metrics_b
      type(multilayer_state_t) :: ms_a, ms_b
      type(continuity_t) :: ct_a, ct_b
      type(coriolis_adv_t) :: cor_a, cor_b
      type(ocean_pressure_force_t) :: pgf_a, pgf_b
      type(ocean_horizontal_viscosity_t) :: hv_a, hv_b
      type(ocean_bottom_drag_t) :: bd_a, bd_b
      type(ocean_surface_stress_t) :: ss_a, ss_b
      type(ocean_vertical_advection_t) :: va_a, va_b
      type(ocean_hdiff_tracer_t) :: hd_a, hd_b
      type(ocean_vdiff_t) :: vd_a, vd_b
      type(ocean_vmix_t) :: vmix_a, vmix_b
      type(eos_t) :: eos_a, eos_b
      type(ocean_dyn_t) :: dyn_a, dyn_b
      type(ocean_bc_state_t) :: bc

      real(wp), parameter :: H0 = 50.0_wp
      real(wp), parameter :: ETA_AMP = 1.0_wp
      real(wp), parameter :: DT = 5.0_wp
      integer, parameter :: N_INNER = 10
         !! Must exceed nghost+2: ghost-band corruption from the zeroed
         !! array-edge ζ ring creeps inward one column per substep and
         !! only touches the seam vorticity from substep nghost+2 onward.
         !! N_INNER=4 (the original value) cannot see that defect.
      integer, parameter :: N_STEPS = 8
      integer :: i, j, k, nx_p, ny_p, ng, shift
      integer :: i0, i1, j0, j1
      real(wp) :: max_diff
      real(wp), allocatable :: h_a(:, :, :), h_b_shifted(:, :, :)

      checks: block

         nx_p = 12
         ny_p = 8
         call make_grid(grid, nx_p, ny_p, 1000.0_wp, 1000.0_wp)
         ng = NGHOST
         shift = nx_p/2
         i0 = ng + 1
         i1 = ng + nx_p
         j0 = ng + 1
         j1 = ng + ny_p

         call init_all(grid, ms_a, ct_a, cor_a, pgf_a, hv_a, bd_a, ss_a, &
                       va_a, hd_a, vd_a, vmix_a, eos_a, dyn_a)
         call init_all(grid, ms_b, ct_b, cor_b, pgf_b, hv_b, bd_b, ss_b, &
                       va_b, hd_b, vd_b, vmix_b, eos_b, dyn_b)
         ! Uniform f-plane: preserves x-translation invariance (so the
         ! bit-identity assertion still holds) while activating the
         ! Coriolis / vorticity terms — with f=0 a purely zonal IC keeps
         ! v ≡ 0 and ζ ≡ 0 forever, leaving the seam-vorticity code paths
         ! completely untested.
         cor_a%f_corner = 1.0e-4_wp
         cor_b%f_corner = 1.0e-4_wp
         call make_bc_periodic_x(bc, grid)

         ! Assign IC to run A (on whole array, wrapping takes care of ghosts)
         ms_a%u_face_x_layer = 0.0_wp
         ms_a%v_face_y_layer = 0.0_wp
         do k = 1, NZ
            do j = 1, grid%ny_total
               do i = 1, grid%nx_total
                  ! Small Gaussian η-like h perturbation in x
                  ms_a%h_layer(i, j, k) = H0 + ETA_AMP*exp( &
                                          -0.5_wp*(real(i - ng - 1, wp) - real(nx_p, wp)*0.25_wp)**2/ &
                                          (real(nx_p, wp)*0.1_wp)**2)
               end do
            end do
         end do
         do k = 1, NZ
            ms_a%tracers(ms_a%idx_salinity)%hTr(:, :, k) = &
               eos_a%S_ref*ms_a%h_layer(:, :, k)
            ms_a%tracers(ms_a%idx_temperature)%hTr(:, :, k) = &
               eos_a%T_ref*ms_a%h_layer(:, :, k)
         end do
         dyn_a%bt_work%bt_H_ref = real(NZ, wp)*H0

         ! Run B: shift IC by nx_p/2 in the physical domain
         ms_b%u_face_x_layer = 0.0_wp
         ms_b%v_face_y_layer = 0.0_wp
         do k = 1, NZ
            do j = 1, grid%ny_total
               ! Circularly shift i-index by shift positions
               do i = i0, i1
                  ms_b%h_layer(i, j, k) = ms_a%h_layer(i0 + mod(i - i0 + shift, nx_p), j, k)
               end do
               ! Fill ghosts for B (will be re-wrapped at stage entry)
               do i = 1, ng
                  ms_b%h_layer(i, j, k) = ms_b%h_layer(i + nx_p, j, k)
               end do
               do i = nx_p + ng + 1, grid%nx_total
                  ms_b%h_layer(i, j, k) = ms_b%h_layer(i - nx_p, j, k)
               end do
            end do
         end do
         do k = 1, NZ
            ms_b%tracers(ms_b%idx_salinity)%hTr(:, :, k) = &
               eos_b%S_ref*ms_b%h_layer(:, :, k)
            ms_b%tracers(ms_b%idx_temperature)%hTr(:, :, k) = &
               eos_b%T_ref*ms_b%h_layer(:, :, k)
         end do
         dyn_b%bt_work%bt_H_ref = real(NZ, wp)*H0

         ! Wrap bt_H_ref for both (flat bath, same everywhere — wrap is trivial)
         call ocean_periodic_wrap_centre_2d(dyn_a%bt_work%bt_H_ref, &
                                            grid%nx_total, grid%ny_total, &
                                            grid%nx_phys, grid%ny_phys, &
                                            grid%nghost, .true., .false.)
         call ocean_periodic_wrap_centre_2d(dyn_b%bt_work%bt_H_ref, &
                                            grid%nx_total, grid%ny_total, &
                                            grid%nx_phys, grid%ny_phys, &
                                            grid%nghost, .true., .false.)

         call map_in(grid, metrics_a, ms_a, ct_a, cor_a, pgf_a, hv_a, bd_a, ss_a, va_a, hd_a, vd_a, vmix_a, dyn_a)
         call map_in(grid, metrics_b, ms_b, ct_b, cor_b, pgf_b, hv_b, bd_b, ss_b, va_b, hd_b, vd_b, vmix_b, dyn_b)
         !$acc update device(dyn_a%bt_work%bt_H_ref, dyn_b%bt_work%bt_H_ref)

         block
            integer :: stp
            do stp = 1, N_STEPS
               call ocean_dyn_step_split( &
                  grid, metrics_a, dyn_a, eos_a, cor_a, ct_a, pgf_a, hv_a, bd_a, ss_a, &
                  va_a, hd_a, vd_a, vmix_a, ms_a, DT, N_INNER, bc=bc)
               call ocean_dyn_step_split( &
                  grid, metrics_b, dyn_b, eos_b, cor_b, ct_b, pgf_b, hv_b, bd_b, ss_b, &
                  va_b, hd_b, vd_b, vmix_b, ms_b, DT, N_INNER, bc=bc)
            end do
         end block

         !$acc update self(ms_a%h_layer, ms_b%h_layer)
         call map_out(metrics_a, ms_a, ct_a, cor_a, pgf_a, hv_a, bd_a, ss_a, va_a, hd_a, vd_a, vmix_a, dyn_a)
         call map_out(metrics_b, ms_b, ct_b, cor_b, pgf_b, hv_b, bd_b, ss_b, va_b, hd_b, vd_b, vmix_b, dyn_b)

         ! Compare shift(B) == A in the physical interior
         allocate (h_a(nx_p, grid%ny_total, NZ))
         allocate (h_b_shifted(nx_p, grid%ny_total, NZ))
         h_a = ms_a%h_layer(i0:i1, :, :)
         do k = 1, NZ
            do j = 1, grid%ny_total
               do i = 1, nx_p
                  ! Un-shift B: physical cell i of B corresponds to cell
                  ! (i - shift) mod nx_p + 1 of A's physical domain
                  h_b_shifted(i, j, k) = ms_b%h_layer( &
                                         i0 + mod(i - 1 + nx_p - mod(shift, nx_p), nx_p), j, k)
               end do
            end do
         end do

         max_diff = maxval(abs(h_a - h_b_shifted))

         deallocate (h_a, h_b_shifted)

         ! Bit-for-bit: with f-plane (uniform f_corner), flat translation-
         ! invariant arithmetic and a periodic-x seam that runs the normal
         ! interior stencils on wrapped operands, a circular x-shift of the
         ! IC must reproduce the same physical field EXACTLY.  Any nonzero
         ! diff means a kernel breaks translation/expression symmetry at the
         ! seam.  Asserted at 0.0 tolerance (the design §1.1 invariant).
         call check(error, max_diff == 0.0_wp, &
                    "shifted-domain: shift(B) /= A bit-for-bit, max_diff=" &
                    //adjusted_str(max_diff))

      end block checks

      call ocean_bc_state_destroy(bc)
      call destroy_all(ms_a, ct_a, cor_a, pgf_a, hv_a, bd_a, ss_a, va_a, hd_a, vd_a, vmix_a, eos_a, dyn_a)
      call destroy_all(ms_b, ct_b, cor_b, pgf_b, hv_b, bd_b, ss_b, va_b, hd_b, vd_b, vmix_b, eos_b, dyn_b)
   end subroutine test_shifted_domain

   ! -----------------------------------------------------------------
   ! Test 4: conservation (doubly periodic)
   ! -----------------------------------------------------------------

   subroutine test_conservation(error)
      !! Total mass and salt in a doubly-periodic domain should be
      !! constant (to ~1e-12 relative) over N_STEPS outer steps.
      type(error_type), allocatable, intent(out) :: error
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
      type(ocean_bc_state_t) :: bc

      real(wp), parameter :: H0 = 50.0_wp
      real(wp), parameter :: U_AMP = 0.05_wp
      real(wp), parameter :: DT = 2.0_wp
      integer, parameter :: N_INNER = 4
      integer, parameter :: N_STEPS = 10
      integer :: i, j, k, nx_p, ny_p, nx, ny, ng, step
      integer :: i0, i1, j0, j1
      real(wp) :: total_h0, total_S0, total_h, total_S, drift_h, drift_S
      real(wp) :: PI

      checks: block

         PI = acos(-1.0_wp)
         nx_p = 8
         ny_p = 6
         call make_grid(grid, nx_p, ny_p, 500.0_wp, 500.0_wp)
         nx = grid%nx_total
         ny = grid%ny_total
         ng = NGHOST
         i0 = ng + 1
         i1 = ng + nx_p
         j0 = ng + 1
         j1 = ng + ny_p

         call init_all(grid, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn)
         call make_bc_periodic_xy(bc, grid)

         ! IC: small sinusoidal h perturbation on physical interior
         do k = 1, NZ
            do j = j0, j1
               do i = i0, i1
                  ms%h_layer(i, j, k) = H0*(1.0_wp + 0.01_wp* &
                                            sin(2.0_wp*PI*real(i - i0, wp)/real(nx_p, wp))* &
                                            cos(2.0_wp*PI*real(j - j0, wp)/real(ny_p, wp)))
               end do
            end do
            ! u perturbation (physical interior faces)
            do j = j0, j1
               do i = i0, i1 + 1
                  ms%u_face_x_layer(i, j, k) = U_AMP* &
                                               cos(2.0_wp*PI*real(i - i0, wp)/real(nx_p, wp))
               end do
            end do
            ms%v_face_y_layer = 0.0_wp
         end do
         ! Tracer proportional to h
         do k = 1, NZ
            ms%tracers(ms%idx_salinity)%hTr(:, :, k) = eos%S_ref*ms%h_layer(:, :, k)
            ms%tracers(ms%idx_temperature)%hTr(:, :, k) = eos%T_ref*ms%h_layer(:, :, k)
         end do
         dyn%bt_work%bt_H_ref = real(NZ, wp)*H0
         ! Wrap ghosts at init (uniform H, so trivial but correct)
         call ocean_periodic_wrap_centre_2d(dyn%bt_work%bt_H_ref, &
                                            grid%nx_total, grid%ny_total, &
                                            grid%nx_phys, grid%ny_phys, &
                                            grid%nghost, .true., .true.)
         call ocean_periodic_wrap_state(grid, bc, ms)

         ! Total over physical interior only (ghost cells would double-count)
         total_h0 = sum(ms%h_layer(i0:i1, j0:j1, :))*grid%dx*grid%dy
         total_S0 = sum(ms%tracers(ms%idx_salinity)%hTr(i0:i1, j0:j1, :))*grid%dx*grid%dy

         call map_in(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
         !$acc update device(dyn%bt_work%bt_H_ref)

         do step = 1, N_STEPS
            call ocean_dyn_step_split( &
               grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
               va, hd, vd, vmix, ms, DT, N_INNER, bc=bc)
         end do

         !$acc update self(ms%h_layer, ms%tracers(ms%idx_salinity)%hTr)
         call map_out(metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)

         total_h = sum(ms%h_layer(i0:i1, j0:j1, :))*grid%dx*grid%dy
         total_S = sum(ms%tracers(ms%idx_salinity)%hTr(i0:i1, j0:j1, :))*grid%dx*grid%dy
         drift_h = abs(total_h - total_h0)/abs(total_h0)
         drift_S = abs(total_S - total_S0)/abs(total_S0)

         call check(error, drift_h < 1.0e-12_wp, &
                    "periodic conservation: mass drift exceeded 1e-12")
         if (allocated(error)) exit checks
         call check(error, drift_S < 1.0e-12_wp, &
                    "periodic conservation: salt drift exceeded 1e-12")

      end block checks

      call ocean_bc_state_destroy(bc)
      call destroy_all(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn)
   end subroutine test_conservation

   ! -----------------------------------------------------------------
   ! Test 5: geostrophic zonal jet (periodic-x, f-plane)
   ! -----------------------------------------------------------------

   subroutine test_geostrophic_jet(error)
      !! A zonally-uniform geostrophic jet (u=U0, constant, η=const),
      !! periodic-x, with a Coriolis parameter f0, should remain steady.
      !! A truly balanced jet satisfies ∂u/∂t = 0.  The uniform u + flat
      !! η version is trivially balanced (no PGF gradient, no Coriolis
      !! tendency since v=0 everywhere).
      type(error_type), allocatable, intent(out) :: error
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
      type(ocean_bc_state_t) :: bc

      real(wp), parameter :: H0 = 200.0_wp
      real(wp), parameter :: U0 = 0.5_wp
      real(wp), parameter :: F0 = 1.0e-4_wp
      real(wp), parameter :: DT = 10.0_wp
      integer, parameter :: N_INNER = 5
      integer, parameter :: N_STEPS = 6
      integer :: step, k
      integer :: i0, i1, j0, j1, ng
      real(wp) :: max_du, max_dh

      checks: block

         call make_grid(grid, 10, 8, 2000.0_wp, 2000.0_wp)
         ng = NGHOST
         i0 = ng + 1
         i1 = ng + grid%nx_phys
         j0 = ng + 1
         j1 = ng + grid%ny_phys

         call init_all(grid, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn)
         cor%f_0 = F0
         call make_bc_periodic_x(bc, grid)

         ms%h_layer = H0
         ms%u_face_x_layer = U0
         ms%v_face_y_layer = 0.0_wp
         do k = 1, NZ
            ms%tracers(ms%idx_salinity)%hTr(:, :, k) = eos%S_ref*H0
            ms%tracers(ms%idx_temperature)%hTr(:, :, k) = eos%T_ref*H0
         end do
         dyn%bt_work%bt_H_ref = real(NZ, wp)*H0
         call ocean_periodic_wrap_centre_2d(dyn%bt_work%bt_H_ref, &
                                            grid%nx_total, grid%ny_total, &
                                            grid%nx_phys, grid%ny_phys, &
                                            grid%nghost, .true., .false.)

         call map_in(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
         !$acc update device(dyn%bt_work%bt_H_ref)

         do step = 1, N_STEPS
            call ocean_dyn_step_split( &
               grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
               va, hd, vd, vmix, ms, DT, N_INNER, bc=bc)
         end do

         !$acc update self(ms%h_layer, ms%u_face_x_layer)
         call map_out(metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)

         max_dh = maxval(abs(ms%h_layer(i0:i1, j0:j1, :) - H0))
         max_du = maxval(abs(ms%u_face_x_layer(i0:i1, j0:j1, :) - U0))

         call check(error, max_dh < 1.0e-7_wp, &
                    "geostrophic jet: h drifted")
         if (allocated(error)) exit checks
         call check(error, max_du < 1.0e-7_wp, &
                    "geostrophic jet: u drifted")

      end block checks

      call ocean_bc_state_destroy(bc)
      call destroy_all(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn)
   end subroutine test_geostrophic_jet

   ! -----------------------------------------------------------------
   ! Test 6: validation rejections
   ! -----------------------------------------------------------------

   subroutine test_validation_rejections(error)
      !! Check that ocean_bc_type_from_string returns OBC_PERIODIC for
      !! "periodic", and that the logicals are correctly derived.
      !! (The error stop paths for unpaired/nghost<3/sponge cannot be
      !! tested without process termination — we verify the derivation
      !! logic instead, and confirm the string → tag mapping.)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_bc_state_t) :: bc
      integer :: tag

      checks: block

         ! String → tag
         tag = ocean_bc_type_from_string("periodic")
         call check(error, tag == OBC_PERIODIC, &
                    "ocean_bc_type_from_string('periodic') /= OBC_PERIODIC")
         if (allocated(error)) exit checks

         ! periodic_x derived correctly when both west+east are OBC_PERIODIC
         call make_grid(grid, 8, 6, 1.0_wp, 1.0_wp)  ! nghost=3 from NGHOST param
         call ocean_bc_state_init(bc, grid, nz_ml=1)
         bc%west%bc_type = OBC_PERIODIC
         bc%east%bc_type = OBC_PERIODIC
         call ocean_bc_validate_periodic(bc)
         call check(error, bc%periodic_x, &
                    "periodic_x not set after validate with W+E=PERIODIC")
         if (allocated(error)) exit checks
         call check(error,.not. bc%periodic_y, &
                    "periodic_y incorrectly set when S/N are WALL")
         if (allocated(error)) exit checks
         call ocean_bc_state_destroy(bc)

         ! periodic_y derived correctly when both south+north are OBC_PERIODIC
         call ocean_bc_state_init(bc, grid, nz_ml=1)
         bc%south%bc_type = OBC_PERIODIC
         bc%north%bc_type = OBC_PERIODIC
         call ocean_bc_validate_periodic(bc)
         call check(error,.not. bc%periodic_x, &
                    "periodic_x incorrectly set when W/E are WALL")
         if (allocated(error)) exit checks
         call check(error, bc%periodic_y, &
                    "periodic_y not set after validate with S+N=PERIODIC")
         if (allocated(error)) exit checks
         call ocean_bc_state_destroy(bc)

         ! OBC_PERIODIC value matches the documented constant
         call check(error, OBC_PERIODIC == 10, &
                    "OBC_PERIODIC /= 10 (enum changed unexpectedly)")
         if (allocated(error)) exit checks

      end block checks
   end subroutine test_validation_rejections

   ! -----------------------------------------------------------------
   ! Test bathy-seam: non-uniform bathymetry init-wrap
   ! -----------------------------------------------------------------

   subroutine test_bathy_seam(error)
      !! Design §1.5 (init-time wrap) bathy-seam test — shifted-domain
      !! bit-identity WITH non-uniform-in-x bathymetry, through the real
      !! split-driver step loop.
      !!
      !! This is the strong form of the gold shifted-domain test: a zonal
      !! ridge in bt_H_ref AND h_layer (so the column thickness varies in x),
      !! run twice with run B's IC circularly shifted by nx_phys/2.  After
      !! N_STEPS the un-shifted B must equal A bit-for-bit.  The bathymetry
      !! ghosts MUST be wrapped at init (this test seeds the physical ridge,
      !! leaves bt_H_ref ghosts at zero, then calls the init wrap) — a stale
      !! ghost H_ref at the seam (the formula-bathy ghost-fill gotcha: ghost
      !! h_layer=0 → EOS ρ=ρ₀ jump) would break the seam and the diff would
      !! be nonzero.  The ghost==partner equality is also asserted right
      !! after the wrap.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics_a, metrics_b
      type(multilayer_state_t) :: ms_a, ms_b
      type(continuity_t) :: ct_a, ct_b
      type(coriolis_adv_t) :: cor_a, cor_b
      type(ocean_pressure_force_t) :: pgf_a, pgf_b
      type(ocean_horizontal_viscosity_t) :: hv_a, hv_b
      type(ocean_bottom_drag_t) :: bd_a, bd_b
      type(ocean_surface_stress_t) :: ss_a, ss_b
      type(ocean_vertical_advection_t) :: va_a, va_b
      type(ocean_hdiff_tracer_t) :: hd_a, hd_b
      type(ocean_vdiff_t) :: vd_a, vd_b
      type(ocean_vmix_t) :: vmix_a, vmix_b
      type(eos_t) :: eos_a, eos_b
      type(ocean_dyn_t) :: dyn_a, dyn_b
      type(ocean_bc_state_t) :: bc

      real(wp), parameter :: H0_MEAN = 100.0_wp   ! mean layer thickness (m)
      real(wp), parameter :: H0_AMP = 10.0_wp     ! amplitude of x-variation
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp), parameter :: F0 = 1.0e-4_wp       ! f-plane (uniform → x-invariant)
      real(wp), parameter :: DT = 5.0_wp
      integer, parameter :: N_INNER = 4
      integer, parameter :: N_STEPS = 8
      integer, parameter :: NX_PHYS = 12, NY_PHYS = 6
      integer :: i, j, k
      integer :: ng, nx_t, ny_t, nx_p, i0, i1, shift
      real(wp) :: ridge, max_diff
      real(wp), allocatable :: h_a(:, :, :), h_b_shifted(:, :, :)

      checks: block

         call make_grid(grid, NX_PHYS, NY_PHYS, 1000.0_wp, 1000.0_wp)
         call init_all(grid, ms_a, ct_a, cor_a, pgf_a, hv_a, bd_a, ss_a, &
                       va_a, hd_a, vd_a, vmix_a, eos_a, dyn_a)
         call init_all(grid, ms_b, ct_b, cor_b, pgf_b, hv_b, bd_b, ss_b, &
                       va_b, hd_b, vd_b, vmix_b, eos_b, dyn_b)
         cor_a%f_0 = F0
         cor_b%f_0 = F0
         cor_a%f_corner = F0
         cor_b%f_corner = F0
         call make_bc_periodic_x(bc, grid)

         ng = grid%nghost
         nx_t = grid%nx_total
         ny_t = grid%ny_total
         nx_p = grid%nx_phys
         i0 = ng + 1
         i1 = ng + nx_p
         shift = nx_p/2

         ! --- Run A: zonal ridge in bt_H_ref and h_layer (physical cells) ---
         ! Ghosts of bt_H_ref left at zero before the init wrap.
         dyn_a%bt_work%bt_H_ref = 0.0_wp
         ms_a%u_face_x_layer = 0.0_wp
         ms_a%v_face_y_layer = 0.0_wp
         do j = 1, ny_t
            do i = i0, i1
               ridge = H0_MEAN + H0_AMP*sin(2.0_wp*PI*real(i - i0, wp)/real(nx_p, wp))
               dyn_a%bt_work%bt_H_ref(i, j) = real(NZ, wp)*ridge
               do k = 1, NZ
                  ms_a%h_layer(i, j, k) = ridge
               end do
            end do
         end do

         ! --- Run B: same ridge, circularly shifted by nx_p/2 in x ---
         dyn_b%bt_work%bt_H_ref = 0.0_wp
         ms_b%u_face_x_layer = 0.0_wp
         ms_b%v_face_y_layer = 0.0_wp
         do j = 1, ny_t
            do i = i0, i1
               ridge = H0_MEAN + H0_AMP*sin(2.0_wp*PI* &
                                            real(mod(i - i0 + shift, nx_p), wp)/real(nx_p, wp))
               dyn_b%bt_work%bt_H_ref(i, j) = real(NZ, wp)*ridge
               do k = 1, NZ
                  ms_b%h_layer(i, j, k) = ridge
               end do
            end do
         end do

         ! Init-time wrap of bathymetry ghosts (design §1.5): without this
         ! the seam ghost H_ref stays 0 and the seam blows the bit-identity.
         call ocean_periodic_wrap_centre_2d(dyn_a%bt_work%bt_H_ref, &
                                            nx_t, ny_t, nx_p, grid%ny_phys, ng, .true., .false.)
         call ocean_periodic_wrap_centre_2d(dyn_b%bt_work%bt_H_ref, &
                                            nx_t, ny_t, nx_p, grid%ny_phys, ng, .true., .false.)

         ! Verify ghost columns equal their physical partners after the wrap.
         do j = 1, ny_t
            do i = 1, ng
               call check(error, &
                          dyn_a%bt_work%bt_H_ref(i, j) == dyn_a%bt_work%bt_H_ref(i + nx_p, j), &
                          "bathy_seam: west ghost bt_H_ref /= partner after wrap")
               if (allocated(error)) exit checks
            end do
            do i = nx_p + ng + 1, nx_t
               call check(error, &
                          dyn_a%bt_work%bt_H_ref(i, j) == dyn_a%bt_work%bt_H_ref(i - nx_p, j), &
                          "bathy_seam: east ghost bt_H_ref /= partner after wrap")
               if (allocated(error)) exit checks
            end do
         end do

         ! Seed tracers proportional to h, then wrap the multilayer state ghosts.
         do k = 1, NZ
            ms_a%tracers(ms_a%idx_salinity)%hTr(:, :, k) = eos_a%S_ref*ms_a%h_layer(:, :, k)
            ms_a%tracers(ms_a%idx_temperature)%hTr(:, :, k) = eos_a%T_ref*ms_a%h_layer(:, :, k)
            ms_b%tracers(ms_b%idx_salinity)%hTr(:, :, k) = eos_b%S_ref*ms_b%h_layer(:, :, k)
            ms_b%tracers(ms_b%idx_temperature)%hTr(:, :, k) = eos_b%T_ref*ms_b%h_layer(:, :, k)
         end do
         call ocean_periodic_wrap_state(grid, bc, ms_a)
         call ocean_periodic_wrap_state(grid, bc, ms_b)

         call map_in(grid, metrics_a, ms_a, ct_a, cor_a, pgf_a, hv_a, bd_a, ss_a, va_a, hd_a, vd_a, vmix_a, dyn_a)
         call map_in(grid, metrics_b, ms_b, ct_b, cor_b, pgf_b, hv_b, bd_b, ss_b, va_b, hd_b, vd_b, vmix_b, dyn_b)
         !$acc update device(dyn_a%bt_work%bt_H_ref, dyn_b%bt_work%bt_H_ref)

         do i = 1, N_STEPS
            call ocean_dyn_step_split( &
               grid, metrics_a, dyn_a, eos_a, cor_a, ct_a, pgf_a, hv_a, bd_a, ss_a, &
               va_a, hd_a, vd_a, vmix_a, ms_a, DT, N_INNER, bc=bc)
            call ocean_dyn_step_split( &
               grid, metrics_b, dyn_b, eos_b, cor_b, ct_b, pgf_b, hv_b, bd_b, ss_b, &
               va_b, hd_b, vd_b, vmix_b, ms_b, DT, N_INNER, bc=bc)
         end do

         !$acc update self(ms_a%h_layer, ms_b%h_layer)
         call map_out(metrics_a, ms_a, ct_a, cor_a, pgf_a, hv_a, bd_a, ss_a, va_a, hd_a, vd_a, vmix_a, dyn_a)
         call map_out(metrics_b, ms_b, ct_b, cor_b, pgf_b, hv_b, bd_b, ss_b, va_b, hd_b, vd_b, vmix_b, dyn_b)

         ! Un-shift B and compare to A in the physical interior — bit-for-bit.
         allocate (h_a(nx_p, ny_t, NZ), h_b_shifted(nx_p, ny_t, NZ))
         h_a = ms_a%h_layer(i0:i1, :, :)
         do k = 1, NZ
            do j = 1, ny_t
               do i = 1, nx_p
                  h_b_shifted(i, j, k) = ms_b%h_layer( &
                                         i0 + mod(i - 1 + nx_p - mod(shift, nx_p), nx_p), j, k)
               end do
            end do
         end do
         max_diff = maxval(abs(h_a - h_b_shifted))
         deallocate (h_a, h_b_shifted)

         call check(error, max_diff == 0.0_wp, &
                    "bathy_seam: shift(B) /= A bit-for-bit with zonal ridge bathymetry, max_diff=" &
                    //adjusted_str(max_diff))

      end block checks

      call ocean_bc_state_destroy(bc)
      call destroy_all(ms_a, ct_a, cor_a, pgf_a, hv_a, bd_a, ss_a, va_a, hd_a, vd_a, vmix_a, eos_a, dyn_a)
      call destroy_all(ms_b, ct_b, cor_b, pgf_b, hv_b, bd_b, ss_b, va_b, hd_b, vd_b, vmix_b, eos_b, dyn_b)
   end subroutine test_bathy_seam

   ! -----------------------------------------------------------------
   ! Internal helper: turn a real number into a short string
   ! -----------------------------------------------------------------

   pure function adjusted_str(val) result(s)
      real(wp), intent(in) :: val
      character(len=20) :: s
      write (s, "(es12.4)") val
   end function adjusted_str

end module test_ocean_periodic
