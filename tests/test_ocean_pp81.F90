!! Unit tests for the Pacanowski-Philander 1981 closure
!! (`vmix_compute_pp81` in `rdb_ocean_vmix`).
!!
!! PP81 produces 3D `kv` (momentum) and `kt` (tracer) eddy diffusivity
!! fields at layer interfaces from the local Richardson number
!! `Ri = N² / |∂u/∂z|²`:
!!
!!   kv = ν_bg + ν_0 / (1 + α·Ri)²
!!   kt = κ_bg + ν_0 / (1 + α·Ri)³
!!
!! Convective limit (Ri < 0): clip Ri = 0 → kv = ν_bg + ν_0.
!!
!! Cases:
!!   * Quiescent stratified column → kv = ν_bg + ν_0, kt = κ_bg + ν_0
!!     (Ri = 0 because shear² is floored, no stable factor applied).
!!   * Strong shear, no stratification → Ri = 0, full ν_0 / κ_0
!!     mixing.  Same kv as quiescent (both Ri = 0).
!!   * Strong stratification + weak shear → Ri large → kv → ν_bg,
!!     kt → κ_bg.  Catches the stability function clamping mixing.
!!   * Unstable column → Ri < 0 → clipped to Ri = 0 → saturates at
!!     ν_bg + ν_0 ≈ 1.01e-2 m²/s.  NOT convective mixing (that is
!!     `vmix_apply_convection`, `test_ocean_convection`) — this only
!!     demonstrates PP81's own Ri<0 clip.
!!   * Boundary interfaces (bed at k=1, surface at k=nz+1) → kv = 0,
!!     kt = 0 by construction (closed BCs).
module test_ocean_pp81
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_vmix, only: ocean_vmix_t, vmix_compute_pp81
   implicit none
   private

   public :: collect_ocean_pp81_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 4

contains

   subroutine collect_ocean_pp81_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("pp81_quiescent_stratified", test_quiescent), &
                  new_unittest("pp81_shear_no_strat_full_mixing", test_shear_no_strat), &
                  new_unittest("pp81_strong_strat_clamps_mixing", test_strong_strat), &
                  new_unittest("pp81_unstable_ri_clip", test_convective), &
                  new_unittest("pp81_boundary_interfaces_zero", test_boundary_zero), &
                  new_unittest("pp81_wall_kv_matches_interior", test_wall_symmetry) &
                  ]
   end subroutine collect_ocean_pp81_tests

   subroutine make_grid(grid, nx_phys, ny_phys, dx, dy)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      real(wp), intent(in) :: dx, dy
      call grid%init(nx_phys, ny_phys, NGHOST, dx, dy)
   end subroutine make_grid

   subroutine run_pp81(grid, ms, vmix)
      !! Map ms + vmix to the device, run compute_pp81, pull kv/kt
      !! back to the host for inspection.
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_vmix_t), intent(inout) :: vmix
      !$acc enter data copyin(ms, vmix)
      call ms%enter_data()
      call vmix%enter_data()
      call vmix_compute_pp81(grid, vmix, ms)
      !$acc update self(vmix%kv, vmix%kt)
      call vmix%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, vmix)
   end subroutine run_pp81

   ! -----------------------------------------------------------------
   ! Cases
   ! -----------------------------------------------------------------

   subroutine test_quiescent(error)
      !! Stratified column at rest.  Shear is zero so shear² is
      !! floored at `shear2_floor` and Ri ~ N²/eps → huge.  With α·Ri
      !! >> 1, kv → ν_bg and kt → κ_bg.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      real(wp), parameter :: H_LAYER = 10.0_wp
      integer :: i_probe, j_probe, k
      real(wp) :: max_kv_excess, max_kt_excess
      checks: block

         call make_grid(grid, 12, 10, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call vmix%init(grid, nz_ml=NZ)
         ms%h_layer = H_LAYER
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         ! Stable stratification: heavier at bed.
         do k = 1, NZ
            ms%rho_layer(:, :, k) = 1030.0_wp - 1.0_wp*real(k - 1, wp)
         end do

         call run_pp81(grid, ms, vmix)

         i_probe = grid%nx_total/2
         j_probe = grid%ny_total/2

         ! Sample interior interfaces only.  With shear floored, Ri is
         ! huge → factor ~ 0 → kv collapses to ν_bg, kt to κ_bg.
         max_kv_excess = 0.0_wp
         max_kt_excess = 0.0_wp
         do k = 2, NZ
            max_kv_excess = max(max_kv_excess, &
                                abs(vmix%kv(i_probe, j_probe, k) - vmix%pp81_nu_bg))
            max_kt_excess = max(max_kt_excess, &
                                abs(vmix%kt(i_probe, j_probe, k) - vmix%pp81_kappa_bg))
         end do

         call check(error, max_kv_excess < 1.0e-6_wp, &
                    "quiescent stratified: kv did not collapse to nu_bg")
         if (allocated(error)) exit checks
         call check(error, max_kt_excess < 1.0e-6_wp, &
                    "quiescent stratified: kt did not collapse to kappa_bg")

      end block checks
      call vmix%destroy(); call ms%destroy()
   end subroutine test_quiescent

   subroutine test_shear_no_strat(error)
      !! Strong shear in x, no stratification.  N² = 0 → Ri = 0 →
      !! factor = 1 → kv = ν_bg + ν_0.  Same kt = κ_bg + ν_0.  This
      !! is the upper limit of PP81 interior mixing.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      real(wp), parameter :: H_LAYER = 10.0_wp
      real(wp), parameter :: U_TOP = 1.0_wp
      real(wp), parameter :: TOL = 1.0e-8_wp
      real(wp) :: expected_kv, expected_kt
      integer :: i_probe, j_probe, k
      checks: block

         call make_grid(grid, 12, 10, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call vmix%init(grid, nz_ml=NZ)
         ms%h_layer = H_LAYER
         ms%rho_layer = 1030.0_wp   ! uniform → N² = 0
         ms%v_face_y_layer = 0.0_wp
         ! Linear u profile in z: u(k) = U_TOP * (k-1)/(NZ-1)
         do k = 1, NZ
            ms%u_face_x_layer(:, :, k) = U_TOP*real(k - 1, wp)/real(NZ - 1, wp)
         end do

         call run_pp81(grid, ms, vmix)

         i_probe = grid%nx_total/2
         j_probe = grid%ny_total/2
         expected_kv = vmix%pp81_nu_bg + vmix%pp81_nu0
         expected_kt = vmix%pp81_kappa_bg + vmix%pp81_nu0

         ! Interior interfaces (k = 2..NZ) should all hit the upper limit.
         do k = 2, NZ
            call check(error, &
                       abs(vmix%kv(i_probe, j_probe, k) - expected_kv) < TOL, &
                       "shear_no_strat: kv didn't hit nu_bg + nu_0")
            if (allocated(error)) exit checks
            call check(error, &
                       abs(vmix%kt(i_probe, j_probe, k) - expected_kt) < TOL, &
                       "shear_no_strat: kt didn't hit kappa_bg + nu_0")
            if (allocated(error)) exit checks
         end do

      end block checks
      call vmix%destroy(); call ms%destroy()
   end subroutine test_shear_no_strat

   subroutine test_strong_strat(error)
      !! Strong stratification + weak shear: Ri = N²/shear² >> 1/α
      !! → factor << 1 → kv well below ν_bg + ν_0.  Demonstrates
      !! the stability function clamping the closure away from full
      !! mixing.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      real(wp), parameter :: H_LAYER = 10.0_wp
      real(wp), parameter :: U_TOP = 1.0e-3_wp   ! weak shear
      real(wp), parameter :: DRHO = 5.0_wp       ! strong stratification
      real(wp) :: kv_max_interior
      integer :: k

      call make_grid(grid, 12, 10, 1.0_wp, 1.0_wp)
      ms%nz_ml = NZ
      call ms%init(grid)
      call vmix%init(grid, nz_ml=NZ)
      ms%h_layer = H_LAYER
      ms%v_face_y_layer = 0.0_wp
      do k = 1, NZ
         ms%u_face_x_layer(:, :, k) = U_TOP*real(k - 1, wp)/real(NZ - 1, wp)
         ms%rho_layer(:, :, k) = 1030.0_wp - DRHO*real(k - 1, wp)
      end do

      call run_pp81(grid, ms, vmix)

      ! At a central probe, kv must be much less than the no-strat
      ! upper limit (ν_bg + ν_0) — orders of magnitude smaller.
      kv_max_interior = maxval(vmix%kv(grid%nghost + 2:grid%nx_total - grid%nghost - 1, &
                                       grid%nghost + 2:grid%ny_total - grid%nghost - 1, &
                                       2:NZ))

      call check(error, kv_max_interior < 0.1_wp*(vmix%pp81_nu_bg + vmix%pp81_nu0), &
                 "strong_strat: kv did not drop below 10% of full mixing")

      call vmix%destroy(); call ms%destroy()
   end subroutine test_strong_strat

   subroutine test_convective(error)
      !! Heavier water sitting on top of lighter water → N² < 0,
      !! Ri < 0.  PP81 clips Ri to zero in the unstable regime,
      !! saturating at kv = ν_bg + ν_0 ≈ 1.01e-2 m²/s.  This is
      !! faithful PP81 and is NOT convective mixing (1.01e-2 m²/s is
      !! ~1% of a real convective diffusivity) — the convective
      !! response lives in `vmix_apply_convection`
      !! (`test_ocean_convection`, `kd_conv` default 1.0 m²/s).  This
      !! case only catches a sign mistake in the PP81 stability factor.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      real(wp), parameter :: H_LAYER = 10.0_wp
      real(wp), parameter :: TOL = 1.0e-8_wp
      integer :: i_probe, j_probe, k
      real(wp) :: expected_kv
      checks: block

         call make_grid(grid, 12, 10, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call vmix%init(grid, nz_ml=NZ)
         ms%h_layer = H_LAYER
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         ! Unstable: density INCREASES with height (heavier at top).
         do k = 1, NZ
            ms%rho_layer(:, :, k) = 1030.0_wp + 1.0_wp*real(k - 1, wp)
         end do

         call run_pp81(grid, ms, vmix)

         i_probe = grid%nx_total/2
         j_probe = grid%ny_total/2
         expected_kv = vmix%pp81_nu_bg + vmix%pp81_nu0

         do k = 2, NZ
            call check(error, &
                       abs(vmix%kv(i_probe, j_probe, k) - expected_kv) < TOL, &
                       "convective: kv didn't hit full-mixing upper limit")
            if (allocated(error)) exit checks
         end do

      end block checks
      call vmix%destroy(); call ms%destroy()
   end subroutine test_convective

   subroutine test_boundary_zero(error)
      !! Bed (k=1) and surface (k=nz+1) interfaces must be exactly
      !! zero — closed BCs feed into the Thomas solver via the
      !! `kv_centre(:, :, 1) = 0` and `kv_centre(:, :, nz+1) = 0`
      !! values that the impl assumes.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      real(wp), parameter :: H_LAYER = 10.0_wp
      real(wp) :: max_bed, max_surf
      integer :: k
      checks: block

         call make_grid(grid, 8, 6, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call vmix%init(grid, nz_ml=NZ)
         ms%h_layer = H_LAYER
         ms%u_face_x_layer = 0.5_wp
         ms%v_face_y_layer = 0.3_wp
         do k = 1, NZ
            ms%rho_layer(:, :, k) = 1030.0_wp - 0.5_wp*real(k - 1, wp)
         end do

         call run_pp81(grid, ms, vmix)

         max_bed = max(maxval(abs(vmix%kv(:, :, 1))), maxval(abs(vmix%kt(:, :, 1))))
         max_surf = max(maxval(abs(vmix%kv(:, :, NZ + 1))), &
                        maxval(abs(vmix%kt(:, :, NZ + 1))))

         call check(error, max_bed < 1.0e-15_wp, "PP81 bed interface kv/kt not zero")
         if (allocated(error)) exit checks
         call check(error, max_surf < 1.0e-15_wp, "PP81 surface interface kv/kt not zero")

      end block checks
      call vmix%destroy(); call ms%destroy()
   end subroutine test_boundary_zero

   subroutine test_wall_symmetry(error)
      !! Regression for the PP81 wall-fallback bug (commit 45421d4).
      !! Under quiescent stratified conditions, kv at the wall faces
      !! (i = 1, i = nx_total, j = 1, j = ny_total) must equal kv at
      !! the immediately adjacent interior column to within FP
      !! roundoff.  Pre-fix the wall block re-set kv to `pp81_nu_bg`
      !! exactly, while the interior formula gave
      !! `pp81_nu_bg + nu0/denom²` — a 1.6e-13 m²/s discontinuity at
      !! every wall column.  That seeded asymmetric vertical mixing
      !! and a 12-hr e-fold baroclinic instability in stratified
      !! closed basins.  Post-fix the interior `do concurrent` covers
      !! the full (i, j) extent, so walls and interior compute via
      !! the same formula.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vmix_t) :: vmix
      real(wp), parameter :: H_LAYER = 10.0_wp
      real(wp), parameter :: TOL = 1.0e-18_wp
      real(wp) :: max_diff
      integer :: k, nx, ny
      checks: block

         call make_grid(grid, 12, 10, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call vmix%init(grid, nz_ml=NZ)
         nx = grid%nx_total
         ny = grid%ny_total
         ms%h_layer = H_LAYER
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         ! Stable stratification, horizontally uniform.
         do k = 1, NZ
            ms%rho_layer(:, :, k) = 1030.0_wp - 1.0_wp*real(k - 1, wp)
         end do

         call run_pp81(grid, ms, vmix)

         ! At every interior interface k, wall-column kv must match
         ! the adjacent-interior-column kv to FP precision.  Pre-fix
         ! the diff was 1.6e-13 m²/s; post-fix it is identical bit-by-
         ! bit because both columns evaluate the same formula.
         max_diff = 0.0_wp
         do k = 2, NZ
            ! West wall vs first interior column
            max_diff = max(max_diff, &
                           abs(vmix%kv(1, 5, k) - vmix%kv(2, 5, k)))
            ! East wall vs last interior column
            max_diff = max(max_diff, &
                           abs(vmix%kv(nx, 5, k) - vmix%kv(nx - 1, 5, k)))
            ! South wall vs first interior row
            max_diff = max(max_diff, &
                           abs(vmix%kv(5, 1, k) - vmix%kv(5, 2, k)))
            ! North wall vs last interior row
            max_diff = max(max_diff, &
                           abs(vmix%kv(5, ny, k) - vmix%kv(5, ny - 1, k)))
         end do

         call check(error, max_diff < TOL, &
                    "PP81 wall kv differs from adjacent interior")

      end block checks
      call vmix%destroy(); call ms%destroy()
   end subroutine test_wall_symmetry

end module test_ocean_pp81
