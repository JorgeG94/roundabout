!! Phase 5a Leith lateral-closure tests.
!!
!! Covers:
!!   1. Zero-flow gives zero ζ → A_h = ah_bg everywhere.
!!   2. Uniform flow gives zero ζ → A_h = ah_bg everywhere.
!!   3. A spatially varying solenoidal flow with non-zero |∇ζ| gives
!!      A_h > ah_bg at some interior faces.
!!   4. ah_max clip is honoured even under a noisy flow that would
!!      otherwise produce huge A_h.
!!   5. End-to-end: horizontal-viscosity kernel with `lateral_mix`
!!      argument produces a Laplacian tendency using the per-face A_h.
!!   6. End-to-end: kernel called WITHOUT `lateral_mix` falls back
!!      bit-identically to the scalar-`nu_h` path.
!!
!! GPU/host pattern (mirrors test_ocean_hvisc):
!!   - All state types (ms, lmix, hv) are attached to the device
!!     via `enter_data` BEFORE the kernel is invoked.
!!   - Kernel runs on device.
!!   - For arrays whose `exit_data` does NOT copyout (lmix's
!!     `ah_face_*`, hv's `du_visc%data`), an explicit
!!     `!$acc update self(...)` syncs device→host before assertion.
!!   - exit_data tears the device attachment down.
module test_ocean_leith
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_lateral_mix, only: ocean_lateral_mix_t, &
                                    ocean_lateral_mix_compute_leith, &
                                    LMIX_NONE, LMIX_LEITH
   use rdb_ocean_horizontal_viscosity, only: ocean_horizontal_viscosity_t, &
                                             ocean_horizontal_viscosity_compute_tendencies
   use testdrive, only: error_type, check, new_unittest, unittest_type
   implicit none
   private

   public :: collect_ocean_leith_tests

contains

   subroutine collect_ocean_leith_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("leith_zero_flow_is_background", test_zero_flow), &
                  new_unittest("leith_uniform_flow_is_background", test_uniform_flow), &
                  new_unittest("leith_solenoidal_flow_excites_A_h", test_solenoidal), &
                  new_unittest("leith_ah_max_clip_honoured", test_ah_max_clip), &
                  new_unittest("hvisc_with_leith_yields_face_coefficient", test_hvisc_with_leith), &
                  new_unittest("hvisc_no_lateral_mix_is_scalar_nuh", test_hvisc_no_lateral_mix) &
                  ]
   end subroutine collect_ocean_leith_tests

   subroutine setup_state(grid, ms, nx, ny, nz, dx)
      type(hgrid_t), intent(inout) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: dx
      call grid%init(nx, ny, 1, dx, dx)
      ms%nz_ml = nz
      call ms%init(grid)
      ms%h_layer = 10.0_wp
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp
   end subroutine setup_state

   subroutine map_in(ms, metrics, grid, lmix, hv)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_metrics_t), intent(inout) :: metrics
      type(hgrid_t), intent(in) :: grid
      type(ocean_lateral_mix_t), intent(inout), optional :: lmix
      type(ocean_horizontal_viscosity_t), intent(inout), optional :: hv
      call make_cartesian_metrics(metrics, grid)
      !$acc enter data copyin(ms)
      call ms%enter_data()
      if (present(lmix)) then
         !$acc enter data copyin(lmix)
         call lmix%enter_data()
      end if
      if (present(hv)) then
         !$acc enter data copyin(hv)
         call hv%enter_data()
      end if
   end subroutine map_in

   subroutine map_out(ms, metrics, lmix, hv)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_metrics_t), intent(inout) :: metrics
      type(ocean_lateral_mix_t), intent(inout), optional :: lmix
      type(ocean_horizontal_viscosity_t), intent(inout), optional :: hv
      if (present(hv)) then
         call hv%exit_data()
         !$acc exit data delete(hv)
      end if
      if (present(lmix)) then
         call lmix%exit_data()
         !$acc exit data delete(lmix)
      end if
      call ms%exit_data()
      !$acc exit data delete(ms)
      call destroy_cartesian_metrics(metrics)
   end subroutine map_out

   subroutine test_zero_flow(error)
      !! Zero velocity → ζ = 0 → ∇ζ = 0 → A_h = ah_bg everywhere.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_lateral_mix_t) :: lmix
      integer, parameter :: NX = 6, NY = 6, NZ = 3
      real(wp), parameter :: DX = 1.0_wp, AH_BG = 1.5_wp
      real(wp) :: max_x, max_y

      call setup_state(grid, ms, NX, NY, NZ, DX)
      call lmix%init(grid, nz_ml=NZ)
      lmix%closure = LMIX_LEITH
      lmix%c_leith = 1.0_wp
      lmix%ah_bg = AH_BG
      lmix%ah_max = 1.0e6_wp
      ! Pre-fill ah_face_* on the host so the copyin captures the
      ! background floor — kernel will overwrite interior faces where
      ! ∇ζ ≠ 0, leave wall faces at ah_bg via its own explicit pass.
      lmix%ah_face_x = AH_BG
      lmix%ah_face_y = AH_BG

      call map_in(ms, metrics, grid, lmix=lmix)
      call ocean_lateral_mix_compute_leith(grid, metrics, lmix, ms)
      !$acc update self(lmix%ah_face_x, lmix%ah_face_y)
      call map_out(ms, metrics, lmix=lmix)

      max_x = maxval(abs(lmix%ah_face_x - AH_BG))
      max_y = maxval(abs(lmix%ah_face_y - AH_BG))
      call check(error, max_x < 1.0e-12_wp .and. max_y < 1.0e-12_wp, &
                 "zero flow should give ah_face = ah_bg everywhere")
      call lmix%destroy()
      call ms%destroy()
   end subroutine test_zero_flow

   subroutine test_uniform_flow(error)
      !! Constant u-velocity → ζ = 0 → A_h = ah_bg.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_lateral_mix_t) :: lmix
      integer, parameter :: NX = 5, NY = 5, NZ = 2
      real(wp), parameter :: DX = 1.0_wp, AH_BG = 0.5_wp
      real(wp) :: max_x, max_y

      call setup_state(grid, ms, NX, NY, NZ, DX)
      ms%u_face_x_layer = 0.25_wp        ! uniform eastward
      ms%v_face_y_layer = -0.1_wp        ! uniform southward
      call lmix%init(grid, nz_ml=NZ)
      lmix%closure = LMIX_LEITH
      lmix%ah_bg = AH_BG
      lmix%ah_max = 1.0e6_wp
      lmix%ah_face_x = AH_BG
      lmix%ah_face_y = AH_BG

      call map_in(ms, metrics, grid, lmix=lmix)
      call ocean_lateral_mix_compute_leith(grid, metrics, lmix, ms)
      !$acc update self(lmix%ah_face_x, lmix%ah_face_y)
      call map_out(ms, metrics, lmix=lmix)

      max_x = maxval(abs(lmix%ah_face_x - AH_BG))
      max_y = maxval(abs(lmix%ah_face_y - AH_BG))
      call check(error, max_x < 1.0e-12_wp .and. max_y < 1.0e-12_wp, &
                 "uniform flow has zero curl, ah_face should equal ah_bg")
      call lmix%destroy()
      call ms%destroy()
   end subroutine test_uniform_flow

   subroutine test_solenoidal(error)
      !! Stamp a spatially varying flow with curl that varies in i
      !! and j.  Some faces should see A_h > ah_bg.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_lateral_mix_t) :: lmix
      integer, parameter :: NX = 8, NY = 8, NZ = 2
      real(wp), parameter :: DX = 1.0_wp, AH_BG = 0.0_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp) :: max_x, max_y
      integer :: i, j, k, nx_tot, ny_tot
      checks: block

         call setup_state(grid, ms, NX, NY, NZ, DX)
         nx_tot = grid%nx_total
         ny_tot = grid%ny_total
         ! Sinusoidal velocity gives non-trivial ζ across the domain.
         do k = 1, NZ
            do j = 1, ny_tot
               do i = 1, nx_tot + 1
                  ms%u_face_x_layer(i, j, k) = 0.5_wp*sin(PI*real(i, wp)/real(nx_tot, wp))* &
                                               cos(PI*real(j, wp)/real(ny_tot, wp))
               end do
            end do
            do j = 1, ny_tot + 1
               do i = 1, nx_tot
                  ms%v_face_y_layer(i, j, k) = 0.3_wp*cos(PI*real(i, wp)/real(nx_tot, wp))* &
                                               sin(PI*real(j, wp)/real(ny_tot, wp))
               end do
            end do
         end do

         call lmix%init(grid, nz_ml=NZ)
         lmix%closure = LMIX_LEITH
         lmix%c_leith = 1.0_wp
         lmix%ah_bg = AH_BG
         lmix%ah_max = 1.0e6_wp

         call map_in(ms, metrics, grid, lmix=lmix)
         call ocean_lateral_mix_compute_leith(grid, metrics, lmix, ms)
         !$acc update self(lmix%ah_face_x, lmix%ah_face_y)
         call map_out(ms, metrics, lmix=lmix)

         ! Some interior face should be > 0 (the floor at ah_bg=0).
         max_x = maxval(lmix%ah_face_x)
         max_y = maxval(lmix%ah_face_y)
         call check(error, max_x > 0.0_wp .or. max_y > 0.0_wp, &
                    "solenoidal flow should produce non-zero A_h somewhere")
         if (allocated(error)) exit checks
         ! And no negative values.
         call check(error, minval(lmix%ah_face_x) >= 0.0_wp, &
                    "Leith A_h should be non-negative on x-faces")
         if (allocated(error)) exit checks
         call check(error, minval(lmix%ah_face_y) >= 0.0_wp, &
                    "Leith A_h should be non-negative on y-faces")

      end block checks
      call lmix%destroy()
      call ms%destroy()
   end subroutine test_solenoidal

   subroutine test_ah_max_clip(error)
      !! Stamp a wild high-vorticity flow and verify A_h never exceeds
      !! ah_max.  Uses a tiny ah_max value so the clip fires.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_lateral_mix_t) :: lmix
      integer, parameter :: NX = 6, NY = 6, NZ = 1
      real(wp), parameter :: DX = 1.0_wp
      real(wp), parameter :: AH_MAX = 0.1_wp
      integer :: i, j
      checks: block

         call setup_state(grid, ms, NX, NY, NZ, DX)
         ! Checkerboard u, v — grid-scale noise.  Produces large |∇ζ|.
         do j = 1, grid%ny_total
            do i = 1, grid%nx_total + 1
               ms%u_face_x_layer(i, j, 1) = real(mod(i + j, 2), wp)*5.0_wp
            end do
         end do
         do j = 1, grid%ny_total + 1
            do i = 1, grid%nx_total
               ms%v_face_y_layer(i, j, 1) = real(mod(i + j, 2), wp)*(-5.0_wp)
            end do
         end do

         call lmix%init(grid, nz_ml=NZ)
         lmix%closure = LMIX_LEITH
         lmix%c_leith = 10.0_wp        ! force big raw A
         lmix%ah_bg = 0.0_wp
         lmix%ah_max = AH_MAX

         call map_in(ms, metrics, grid, lmix=lmix)
         call ocean_lateral_mix_compute_leith(grid, metrics, lmix, ms)
         !$acc update self(lmix%ah_face_x, lmix%ah_face_y)
         call map_out(ms, metrics, lmix=lmix)

         call check(error, maxval(lmix%ah_face_x) <= AH_MAX + 1.0e-12_wp, &
                    "ah_face_x must respect ah_max clip")
         if (allocated(error)) exit checks
         call check(error, maxval(lmix%ah_face_y) <= AH_MAX + 1.0e-12_wp, &
                    "ah_face_y must respect ah_max clip")

      end block checks
      call lmix%destroy()
      call ms%destroy()
   end subroutine test_ah_max_clip

   subroutine test_hvisc_with_leith(error)
      !! The horizontal-viscosity kernel called with `lateral_mix=` and
      !! `closure=LMIX_LEITH` must read per-face A_h, not scalar nu_h.
      !! Test by setting nu_h to a "wrong" value and lateral_mix's
      !! face viscosity to known constants — the Laplacian tendency
      !! must reflect the face values.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_lateral_mix_t) :: lmix
      type(ocean_horizontal_viscosity_t) :: hv
      integer, parameter :: NX = 6, NY = 6, NZ = 2
      real(wp), parameter :: DX = 1.0_wp
      real(wp), parameter :: AH_FORCED = 25.0_wp
      real(wp), parameter :: NU_WRONG = -999.0_wp
      real(wp) :: expected_du_visc
      integer :: i, j, k

      call setup_state(grid, ms, NX, NY, NZ, DX)
      ! u(i, j, k) = i² → d²u/dx² = 2/dx² + d²u/dy² = 0.
      do k = 1, NZ
         do j = 1, grid%ny_total
            do i = 1, grid%nx_total + 1
               ms%u_face_x_layer(i, j, k) = real(i*i, wp)
            end do
         end do
      end do

      call lmix%init(grid, nz_ml=NZ)
      lmix%closure = LMIX_LEITH
      ! Stuff ah_face_* with a known constant value to bypass the
      ! compute_leith step — this test isolates the dispatch.
      lmix%ah_face_x = AH_FORCED
      lmix%ah_face_y = AH_FORCED

      call hv%init(grid, nz_ml=NZ)
      hv%nu_h = NU_WRONG    ! deliberately wrong to verify it's ignored

      call map_in(ms, metrics, grid, lmix=lmix, hv=hv)
      call ocean_horizontal_viscosity_compute_tendencies(grid, metrics, hv, ms, lmix)
      !$acc update self(hv%du_visc%data)
      call map_out(ms, metrics, lmix=lmix, hv=hv)

      ! Interior u-face Laplacian: u(i, j) = i² → d²u/dx² = 2/dx²,
      ! d²u/dy² = 0.  du_visc = AH_FORCED * 2 / (DX*DX) = 50 at DX=1.
      expected_du_visc = AH_FORCED*2.0_wp/(DX*DX)
      call check(error, abs(hv%du_visc%data(3, 3, 1) - expected_du_visc) < 1.0e-10_wp, &
                 "hvisc with Leith should use ah_face_x for the coefficient")

      call hv%destroy()
      call lmix%destroy()
      call ms%destroy()
   end subroutine test_hvisc_with_leith

   subroutine test_hvisc_no_lateral_mix(error)
      !! The horizontal-viscosity kernel called WITHOUT `lateral_mix`
      !! preserves the scalar-nu_h behaviour bit-identically.  Regression
      !! check that the closure dispatch doesn't change pre-existing
      !! call sites.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_horizontal_viscosity_t) :: hv
      integer, parameter :: NX = 6, NY = 6, NZ = 2
      real(wp), parameter :: DX = 1.0_wp
      real(wp), parameter :: NU_H = 7.0_wp
      real(wp) :: expected_du_visc
      integer :: i, j, k

      call setup_state(grid, ms, NX, NY, NZ, DX)
      do k = 1, NZ
         do j = 1, grid%ny_total
            do i = 1, grid%nx_total + 1
               ms%u_face_x_layer(i, j, k) = real(i*i, wp)
            end do
         end do
      end do

      call hv%init(grid, nz_ml=NZ)
      hv%nu_h = NU_H

      call map_in(ms, metrics, grid, hv=hv)
      call ocean_horizontal_viscosity_compute_tendencies(grid, metrics, hv, ms)
      !$acc update self(hv%du_visc%data)
      call map_out(ms, metrics, hv=hv)

      expected_du_visc = NU_H*2.0_wp/(DX*DX)
      call check(error, abs(hv%du_visc%data(3, 3, 1) - expected_du_visc) < 1.0e-10_wp, &
                 "hvisc without lateral_mix must equal nu_h * Laplacian")

      call hv%destroy()
      call ms%destroy()
   end subroutine test_hvisc_no_lateral_mix

end module test_ocean_leith
