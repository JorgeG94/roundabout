!! Unit tests for the convective KPP boundary-layer extension
!! (`vmix_apply_kpp_overlay` with surface_flux Phase 2
!! pathway active).
!!
!! Phase 2 adds the convective velocity scale `w_*` on top of the
!! Phase 1 shear-driven `u_*`:
!!   w_*³ = max(0, -B_0) · h_b
!!   w_s = √(u_*² + w_*²)
!!   kv_kpp = h_b · w_s · G(σ)
!! where B_0 = (g/ρ₀)·(α·F_T - β·F_S) is the surface buoyancy flux.
!! B_0 < 0 (cooling / salting) → w_* > 0 (destabilizing → convective
!! mixing).
!!
!! Cases:
!!   * Stabilizing flux (Q_heat > 0) → B_0 > 0 → w_* = 0.  kv
!!     matches the Phase 1 shear-only result bit-for-bit (passing
!!     the same overlay with `sf` absent vs `sf` with stabilizing
!!     fluxes should produce the same kv).
!!   * Destabilizing cooling (Q_heat < 0) → B_0 < 0 → w_* > 0.
!!     With zero wind stress, the Phase 1 result is `kv ≈ 0` (just
!!     PP81 interior).  Phase 2 with cooling lifts kv well above
!!     the interior baseline.
!!   * Combined shear + convective vs shear-only — both active,
!!     verify kv > shear-only at the same h_b (w_s = √(u_*² + w_*²) > u_*).
module test_ocean_kpp_convective
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t
   use rdb_ocean_vmix, only: ocean_vmix_t, &
                             vmix_apply_kpp_overlay
   implicit none
   private

   public :: collect_ocean_kpp_convective_tests

   integer, parameter :: NGHOST = 2

contains

   subroutine collect_ocean_kpp_convective_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("convective_kpp_stabilizing_no_op", test_stabilizing_no_op), &
                  new_unittest("convective_kpp_cooling_deepens_mixing", test_cooling_deepens), &
                  new_unittest("convective_kpp_combined_exceeds_shear", test_combined_exceeds_shear) &
                  ]
   end subroutine collect_ocean_kpp_convective_tests

   subroutine make_grid(grid, nx_phys, ny_phys, dx, dy)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      real(wp), intent(in) :: dx, dy
      call grid%init(nx_phys, ny_phys, NGHOST, dx, dy)
   end subroutine make_grid

   subroutine setup_state(grid, ms, vmix, ss, sf, nz)
      !! Mixed-layer-over-thermocline IC + scalar wind stress.
      !! Heat / salt fluxes set per test by the caller AFTER this
      !! returns.
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_vmix_t), intent(inout) :: vmix
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_surface_flux_t), intent(inout) :: sf
      integer, intent(in) :: nz
      integer, parameter :: NZ_ML = 3
      real(wp), parameter :: H_LAYER = 25.0_wp
      real(wp), parameter :: U_ML = 0.5_wp
      integer :: k

      ms%nz_ml = nz
      call ms%init(grid)
      call vmix%init(grid, nz_ml=nz)
      call ss%init(grid)
      call sf%init(grid)

      ms%h_layer = H_LAYER
      ms%v_face_y_layer = 0.0_wp
      do k = 1, nz
         if (k > nz - NZ_ML) then
            ms%u_face_x_layer(:, :, k) = U_ML
            ms%rho_layer(:, :, k) = 1027.0_wp
         else
            ms%u_face_x_layer(:, :, k) = 0.0_wp
            ms%rho_layer(:, :, k) = 1030.0_wp - 0.5_wp*real(k - 1, wp)
         end if
      end do
   end subroutine setup_state

   subroutine run_kpp(grid, ms, vmix, ss, sf, use_sf)
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_vmix_t), intent(inout) :: vmix
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_surface_flux_t), intent(inout) :: sf
      logical, intent(in) :: use_sf
      !$acc enter data copyin(ms, vmix, ss, sf)
      call ms%enter_data()
      call vmix%enter_data()
      if (use_sf) then
         call vmix_apply_kpp_overlay(grid, vmix, ms, ss, sf)
      else
         block
            ! Zero-flux slot: sf is REQUIRED post-A7; zero fields = the
            ! old no-flux behaviour.
            type(ocean_surface_flux_t) :: sf_zero
            call sf_zero%init(grid)
            !$acc enter data copyin(sf_zero)
            call sf_zero%enter_data()
            call vmix_apply_kpp_overlay(grid, vmix, ms, ss, sf_zero)
            call sf_zero%exit_data()
            !$acc exit data delete(sf_zero)
            call sf_zero%destroy()
         end block
      end if
      !$acc update self(vmix%kv, vmix%kt, vmix%bl_depth)
      call vmix%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, vmix, ss, sf)
   end subroutine run_kpp

   ! -----------------------------------------------------------------
   ! Cases
   ! -----------------------------------------------------------------

   subroutine test_stabilizing_no_op(error)
      !! Heating with no salt flux → B_0 > 0 → max(0, -B_0) = 0 →
      !! w_* = 0.  The Phase 2 path with stabilizing fluxes must
      !! produce identical kv to the Phase 1 path.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix_p1, vmix_p2
      type(ocean_surface_stress_t) :: ss
      type(ocean_surface_flux_t) :: sf
      real(wp), allocatable :: kv_p1(:, :, :)
      integer, parameter :: NZ = 8
      real(wp) :: max_diff

      call make_grid(grid, 14, 12, 1.0_wp, 1.0_wp)
      call setup_state(grid, ms, vmix_p1, ss, sf, NZ)
      call ss%set_wind_stress_const(0.1_wp, 0.0_wp)

      ! Phase 1 reference (no surface flux argument)
      call run_kpp(grid, ms, vmix_p1, ss, sf, .false.)
      allocate (kv_p1, source=vmix_p1%kv)

      ! Phase 2 with strongly stabilizing flux
      call vmix_p2%init(grid, nz_ml=NZ)
      call sf%set_surface_flux_const(500.0_wp, 0.0_wp)     ! 500 W/m² downward (heating)
      call run_kpp(grid, ms, vmix_p2, ss, sf, .true.)

      max_diff = maxval(abs(vmix_p2%kv - kv_p1))
      call check(error, max_diff < 1.0e-12_wp, &
                 "stabilizing flux changed kv vs shear-only")

      deallocate (kv_p1)
      call sf%destroy(); call ss%destroy(); call vmix_p2%destroy()
      call vmix_p1%destroy(); call ms%destroy()
   end subroutine test_stabilizing_no_op

   subroutine test_cooling_deepens(error)
      !! Strong surface cooling (Q_heat < 0) over the same mixed-
      !! layer-over-thermocline IC.  With zero wind stress, the
      !! Phase 1 KPP overlay produces no overlay (u_* = 0 → w_s = 0
      !! → kv_kpp = 0).  Phase 2 with cooling drives w_* > 0 so
      !! kv inside the BL rises well above the interior PP81
      !! baseline.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      type(ocean_surface_stress_t) :: ss
      type(ocean_surface_flux_t) :: sf
      integer, parameter :: NZ = 8
      real(wp) :: max_kv_interior

      call make_grid(grid, 14, 12, 1.0_wp, 1.0_wp)
      call setup_state(grid, ms, vmix, ss, sf, NZ)
      call ss%set_wind_stress_const(0.0_wp, 0.0_wp)     ! no wind
      call sf%set_surface_flux_const(-500.0_wp, 0.0_wp)  ! 500 W/m² upward (cooling)

      call run_kpp(grid, ms, vmix, ss, sf, .true.)

      max_kv_interior = maxval(vmix%kv(grid%nghost + 1:grid%nx_total - grid%nghost, &
                                       grid%nghost + 1:grid%ny_total - grid%nghost, &
                                       2:NZ))

      ! Sanity: kv must exceed the bare PP81 background by at least
      ! 10x — convective mixing is the dominant signal.
      call check(error, max_kv_interior > 10.0_wp*vmix%pp81_nu_bg, &
                 "convective KPP didn't lift kv above 10x background")

      call sf%destroy(); call ss%destroy(); call vmix%destroy(); call ms%destroy()
   end subroutine test_cooling_deepens

   subroutine test_combined_exceeds_shear(error)
      !! Both shear (wind) and convection (cooling) active.  The
      !! combined w_s = √(u_*² + w_*²) is strictly larger than the
      !! shear-only u_*, so the Phase 2 kv must exceed the Phase 1
      !! kv at the same h_b.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix_p1, vmix_p2
      type(ocean_surface_stress_t) :: ss
      type(ocean_surface_flux_t) :: sf
      real(wp), allocatable :: kv_p1(:, :, :)
      integer, parameter :: NZ = 8
      real(wp) :: max_p1, max_p2

      call make_grid(grid, 14, 12, 1.0_wp, 1.0_wp)
      call setup_state(grid, ms, vmix_p1, ss, sf, NZ)
      call ss%set_wind_stress_const(0.1_wp, 0.0_wp)

      ! Phase 1 baseline (shear only)
      call run_kpp(grid, ms, vmix_p1, ss, sf, .false.)
      allocate (kv_p1, source=vmix_p1%kv)
      max_p1 = maxval(kv_p1(grid%nghost + 1:grid%nx_total - grid%nghost, &
                            grid%nghost + 1:grid%ny_total - grid%nghost, 2:NZ))

      ! Phase 2 with shear + cooling
      call vmix_p2%init(grid, nz_ml=NZ)
      call sf%set_surface_flux_const(-500.0_wp, 0.0_wp)
      call run_kpp(grid, ms, vmix_p2, ss, sf, .true.)
      max_p2 = maxval(vmix_p2%kv(grid%nghost + 1:grid%nx_total - grid%nghost, &
                                 grid%nghost + 1:grid%ny_total - grid%nghost, 2:NZ))

      call check(error, max_p2 > max_p1, &
                 "combined shear+convective kv didn't exceed shear-only")

      deallocate (kv_p1)
      call sf%destroy(); call ss%destroy(); call vmix_p2%destroy()
      call vmix_p1%destroy(); call ms%destroy()
   end subroutine test_combined_exceeds_shear

end module test_ocean_kpp_convective
