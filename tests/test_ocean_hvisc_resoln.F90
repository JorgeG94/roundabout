!! Gap 1 — resolution-function scaling of the momentum horizontal
!! viscosity (MOM6 RESOLN_SCALED_KH analogue, Hallberg 2013).
!!
!! The dynamic LAPLACIAN lateral-viscosity coefficient `A_h` (Leith /
!! Smagorinsky_KH) is multiplied, BEFORE the `ah_max` clamp, by the
!! VarMix resolution function `Res_fn ∈ [0,1]` interpolated onto each
!! u/v face.  This turns the closure OFF where the deformation radius is
!! well resolved (`Res_fn → 0`) and leaves it ON where it is not
!! (`Res_fn → 1`).  The BIHARMONIC `nu_4` (Smagorinsky_AH / Leith_AH) is
!! deliberately NOT scaled — strict MOM6 parity: the resolution function
!! suppresses only the scale-non-selective Laplacian, while the ∝k⁴
!! biharmonic already spares the resolved scales.
!!
!! Covers:
!!   1. Knob OFF + res_fn passed ⇒ coefficient UNCHANGED from the
!!      no-scaling baseline (bit-identity; default off is byte-identical).
!!   2. Knob ON + Res_fn ≡ 1 (coarse / unresolved) ⇒ coefficient
!!      bit-identical to baseline (the identity scaling).
!!   3. Knob ON + Res_fn ≡ 0.3 (partially resolved) ⇒ the DYNAMIC part
!!      is suppressed: at faces where the baseline exceeds `ah_bg`, the
!!      scaled coefficient is strictly smaller and floors at `ah_bg`.
!!   4. Knob ON + Res_fn ≡ 0 (fully resolved) ⇒ coefficient suppressed
!!      all the way to the background floor `ah_bg` everywhere.
!!   5. Smagorinsky_AH biharmonic `nu_4` is NOT resolution-scaled —
!!      toggling the knob leaves it bit-identical (strict MOM6 parity).
!!
!! GPU/host pattern (mirrors test_ocean_leith): all state mapped to the
!! device via `enter_data` before the kernel; the per-face arrays are
!! `update self`'d back before assertion.  The `res_fn_*` arrays are
!! plain local copyins (they are dummy `intent(in)` to the kernel).
module test_ocean_hvisc_resoln
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_lateral_mix, only: ocean_lateral_mix_t, &
                                    ocean_lateral_mix_compute_leith, &
                                    ocean_lateral_mix_compute_smag_ah, &
                                    LMIX_LEITH
   use testdrive, only: error_type, check, new_unittest, unittest_type
   implicit none
   private

   public :: collect_ocean_hvisc_resoln_tests

contains

   subroutine collect_ocean_hvisc_resoln_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("resoln_knob_off_is_baseline", test_knob_off), &
                  new_unittest("resoln_unity_is_baseline", test_res_unity), &
                  new_unittest("resoln_partial_suppresses", test_res_partial), &
                  new_unittest("resoln_zero_collapses_to_bg", test_res_zero), &
                  new_unittest("resoln_does_not_scale_smag_ah_nu4", test_smag_ah) &
                  ]
   end subroutine collect_ocean_hvisc_resoln_tests

   subroutine setup_solenoidal(grid, ms, nx, ny, nz, dx)
      !! Build a sinusoidal (non-zero |∇ζ|, non-zero strain) flow so the
      !! dynamic closure produces A_h > ah_bg at interior faces.
      type(hgrid_t), intent(inout) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: dx
      real(wp), parameter :: PI = acos(-1.0_wp)
      integer :: i, j, k
      call grid%init(nx, ny, 1, dx, dx)
      ms%nz_ml = nz
      call ms%init(grid)
      ms%h_layer = 10.0_wp
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp
      do k = 1, nz
         do j = 1, grid%ny_total
            do i = 1, grid%nx_total + 1
               ms%u_face_x_layer(i, j, k) = 0.5_wp* &
                                            sin(PI*real(i, wp)/real(grid%nx_total, wp))* &
                                            cos(PI*real(j, wp)/real(grid%ny_total, wp))
            end do
         end do
         do j = 1, grid%ny_total + 1
            do i = 1, grid%nx_total
               ms%v_face_y_layer(i, j, k) = 0.3_wp* &
                                            cos(PI*real(i, wp)/real(grid%nx_total, wp))* &
                                            sin(PI*real(j, wp)/real(grid%ny_total, wp))
            end do
         end do
      end do
   end subroutine setup_solenoidal

   subroutine init_lmix(lmix, grid, nz, resoln_on)
      type(ocean_lateral_mix_t), intent(inout) :: lmix
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      logical, intent(in) :: resoln_on
      call lmix%init(grid, nz_ml=nz)
      lmix%closure = LMIX_LEITH
      lmix%c_leith = 1.0_wp
      lmix%ah_bg = 0.0_wp
      lmix%ah_max = 1.0e9_wp
      lmix%resoln_scaled_visc = resoln_on
   end subroutine init_lmix

   subroutine run_leith(grid, ms, metrics, lmix, res_u, res_v, use_res)
      !! Map state + (optionally) the res_fn fields to the device, run the
      !! Leith kernel (with or without the resolution scaling), pull the
      !! per-face coefficient back to the host.
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_metrics_t), intent(inout) :: metrics
      type(ocean_lateral_mix_t), intent(inout) :: lmix
      real(wp), intent(in) :: res_u(grid%nx_total + 1, grid%ny_total)
      real(wp), intent(in) :: res_v(grid%nx_total, grid%ny_total + 1)
      logical, intent(in) :: use_res
      call make_cartesian_metrics(metrics, grid)
      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(lmix)
      call lmix%enter_data()
      if (use_res) then
         !$acc enter data copyin(res_u, res_v)
         call ocean_lateral_mix_compute_leith(grid, metrics, lmix, ms, &
                                              res_fn_u=res_u, res_fn_v=res_v)
         !$acc exit data delete(res_u, res_v)
      else
         call ocean_lateral_mix_compute_leith(grid, metrics, lmix, ms)
      end if
      !$acc update self(lmix%ah_face_x, lmix%ah_face_y)
      call lmix%exit_data()
      !$acc exit data delete(lmix)
      call ms%exit_data()
      !$acc exit data delete(ms)
      call destroy_cartesian_metrics(metrics)
   end subroutine run_leith

   subroutine test_knob_off(error)
      !! Knob OFF but res_fn (≡ 0.3) supplied ⇒ the scaling must NOT fire ⇒
      !! coefficient bit-identical to the no-res baseline.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms_base, ms_off
      type(ocean_metrics_t) :: metrics
      type(ocean_lateral_mix_t) :: lmix_base, lmix_off
      integer, parameter :: NX = 8, NY = 8, NZ = 2
      real(wp), parameter :: DX = 1.0_wp
      real(wp), allocatable :: res_u(:, :), res_v(:, :)
      real(wp), allocatable :: base_x(:, :, :), base_y(:, :, :)

      ! Baseline (no res scaling).
      call setup_solenoidal(grid, ms_base, NX, NY, NZ, DX)
      call init_lmix(lmix_base, grid, NZ, resoln_on=.false.)
      allocate (res_u(grid%nx_total + 1, grid%ny_total), source=0.3_wp)
      allocate (res_v(grid%nx_total, grid%ny_total + 1), source=0.3_wp)
      call run_leith(grid, ms_base, metrics, lmix_base, res_u, res_v, use_res=.false.)
      base_x = lmix_base%ah_face_x
      base_y = lmix_base%ah_face_y

      ! Knob OFF, res_fn supplied — must equal baseline.
      call setup_solenoidal(grid, ms_off, NX, NY, NZ, DX)
      call init_lmix(lmix_off, grid, NZ, resoln_on=.false.)
      call run_leith(grid, ms_off, metrics, lmix_off, res_u, res_v, use_res=.true.)

      call check(error, maxval(abs(lmix_off%ah_face_x - base_x)) < 1.0e-13_wp .and. &
                 maxval(abs(lmix_off%ah_face_y - base_y)) < 1.0e-13_wp, &
                 "knob off must be bit-identical to the no-scaling baseline")
      call lmix_base%destroy()
      call lmix_off%destroy()
      call ms_base%destroy()
      call ms_off%destroy()
   end subroutine test_knob_off

   subroutine test_res_unity(error)
      !! Knob ON, Res_fn ≡ 1 (coarse / fully unresolved) ⇒ identity scaling
      !! ⇒ bit-identical to baseline.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms_base, ms_res
      type(ocean_metrics_t) :: metrics
      type(ocean_lateral_mix_t) :: lmix_base, lmix_res
      integer, parameter :: NX = 8, NY = 8, NZ = 2
      real(wp), parameter :: DX = 1.0_wp
      real(wp), allocatable :: res_u(:, :), res_v(:, :)
      real(wp), allocatable :: base_x(:, :, :), base_y(:, :, :)

      call setup_solenoidal(grid, ms_base, NX, NY, NZ, DX)
      call init_lmix(lmix_base, grid, NZ, resoln_on=.false.)
      allocate (res_u(grid%nx_total + 1, grid%ny_total), source=0.0_wp)
      allocate (res_v(grid%nx_total, grid%ny_total + 1), source=0.0_wp)
      call run_leith(grid, ms_base, metrics, lmix_base, res_u, res_v, use_res=.false.)
      base_x = lmix_base%ah_face_x
      base_y = lmix_base%ah_face_y

      res_u = 1.0_wp
      res_v = 1.0_wp
      call setup_solenoidal(grid, ms_res, NX, NY, NZ, DX)
      call init_lmix(lmix_res, grid, NZ, resoln_on=.true.)
      call run_leith(grid, ms_res, metrics, lmix_res, res_u, res_v, use_res=.true.)

      call check(error, maxval(abs(lmix_res%ah_face_x - base_x)) < 1.0e-13_wp .and. &
                 maxval(abs(lmix_res%ah_face_y - base_y)) < 1.0e-13_wp, &
                 "Res_fn=1 (unresolved) must leave the coefficient unchanged")
      call lmix_base%destroy()
      call lmix_res%destroy()
      call ms_base%destroy()
      call ms_res%destroy()
   end subroutine test_res_unity

   subroutine test_res_partial(error)
      !! Knob ON, Res_fn ≡ 0.3 ⇒ where the baseline A_h exceeds ah_bg the
      !! scaled value is 0.3× the baseline (ah_bg=0 ⇒ no floor interference)
      !! ⇒ strictly smaller, and ≈ 0.3× to round-off.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms_base, ms_res
      type(ocean_metrics_t) :: metrics
      type(ocean_lateral_mix_t) :: lmix_base, lmix_res
      integer, parameter :: NX = 8, NY = 8, NZ = 2
      real(wp), parameter :: DX = 1.0_wp, RES = 0.3_wp
      real(wp), allocatable :: res_u(:, :), res_v(:, :)
      real(wp), allocatable :: base_x(:, :, :)
      real(wp) :: max_base, max_diff
      checks: block

         call setup_solenoidal(grid, ms_base, NX, NY, NZ, DX)
         call init_lmix(lmix_base, grid, NZ, resoln_on=.false.)
         allocate (res_u(grid%nx_total + 1, grid%ny_total), source=0.0_wp)
         allocate (res_v(grid%nx_total, grid%ny_total + 1), source=0.0_wp)
         call run_leith(grid, ms_base, metrics, lmix_base, res_u, res_v, use_res=.false.)
         base_x = lmix_base%ah_face_x
         max_base = maxval(base_x)

         res_u = RES
         res_v = RES
         call setup_solenoidal(grid, ms_res, NX, NY, NZ, DX)
         call init_lmix(lmix_res, grid, NZ, resoln_on=.true.)
         call run_leith(grid, ms_res, metrics, lmix_res, res_u, res_v, use_res=.true.)

         ! There IS dynamic dissipation to suppress.
         call check(error, max_base > 1.0e-6_wp, &
                    "baseline should excite a non-trivial A_h to suppress")
         if (allocated(error)) exit checks
         ! Scaled <= baseline everywhere (Res_fn < 1 cannot increase A_h).
         call check(error, maxval(lmix_res%ah_face_x - base_x) <= 1.0e-12_wp, &
                    "Res_fn=0.3 must not increase the coefficient")
         if (allocated(error)) exit checks
         ! And exactly 0.3× the baseline (ah_bg=0 ⇒ the floor is inert).
         max_diff = maxval(abs(lmix_res%ah_face_x - RES*base_x))
         call check(error, max_diff < 1.0e-10_wp*max_base + 1.0e-12_wp, &
                    "Res_fn=0.3 must scale the coefficient to 0.3x the baseline")
         if (allocated(error)) exit checks
         ! Some face is strictly suppressed.
         call check(error, maxval(base_x - lmix_res%ah_face_x) > 1.0e-7_wp, &
                    "Res_fn=0.3 must strictly suppress the coefficient somewhere")

      end block checks
      call lmix_base%destroy()
      call lmix_res%destroy()
      call ms_base%destroy()
      call ms_res%destroy()
   end subroutine test_res_partial

   subroutine test_res_zero(error)
      !! Knob ON, Res_fn ≡ 0 (fully resolved) ⇒ the dynamic coefficient is
      !! suppressed all the way to the background floor `ah_bg` everywhere.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_lateral_mix_t) :: lmix
      integer, parameter :: NX = 8, NY = 8, NZ = 2
      real(wp), parameter :: DX = 1.0_wp, AH_BG = 2.0_wp
      real(wp), allocatable :: res_u(:, :), res_v(:, :)

      call setup_solenoidal(grid, ms, NX, NY, NZ, DX)
      call init_lmix(lmix, grid, NZ, resoln_on=.true.)
      lmix%ah_bg = AH_BG       ! non-zero floor: collapsed value must hit it
      allocate (res_u(grid%nx_total + 1, grid%ny_total), source=0.0_wp)
      allocate (res_v(grid%nx_total, grid%ny_total + 1), source=0.0_wp)
      call run_leith(grid, ms, metrics, lmix, res_u, res_v, use_res=.true.)

      call check(error, maxval(abs(lmix%ah_face_x - AH_BG)) < 1.0e-12_wp .and. &
                 maxval(abs(lmix%ah_face_y - AH_BG)) < 1.0e-12_wp, &
                 "Res_fn=0 (fully resolved) must collapse A_h to ah_bg")
      call lmix%destroy()
      call ms%destroy()
   end subroutine test_res_zero

   subroutine test_smag_ah(error)
      !! Strict MOM6 parity: the Smagorinsky_AH biharmonic `nu_4` is NOT
      !! resolution-scaled.  Toggling `resoln_scaled_visc` must leave the
      !! per-face nu4 bit-identical — and the closure must actually excite a
      !! non-trivial nu4 > nu4_bg, so the guard is not vacuous.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms_off, ms_on
      type(ocean_metrics_t) :: metrics
      type(ocean_lateral_mix_t) :: lmix_off, lmix_on
      integer, parameter :: NX = 8, NY = 8, NZ = 2
      ! nu4_bg = 0 so the dynamic Smag_AH coefficient shows through rather
      ! than being masked by the floor (with a non-zero floor the strain on
      ! this flow is far below it, pinning nu4 and making the test vacuous).
      real(wp), parameter :: DX = 1.0_wp, NU4_BG = 0.0_wp
      real(wp), allocatable :: ref_x(:, :, :)
      checks: block

         ! Knob OFF — baseline biharmonic nu4.
         call setup_solenoidal(grid, ms_off, NX, NY, NZ, DX)
         call init_smag_ah(lmix_off, grid, NZ, NU4_BG, resoln_on=.false.)
         call run_smag_ah(grid, ms_off, metrics, lmix_off)
         ref_x = lmix_off%nu4_face_x

         ! Closure must excite a non-trivial nu4 somewhere (non-vacuous guard).
         call check(error, maxval(ref_x) > 1.0e-6_wp, &
                    "Smag_AH should excite a non-trivial nu4 to test against")
         if (allocated(error)) exit checks

         ! Knob ON — must be bit-identical (biharmonic is never res-scaled).
         call setup_solenoidal(grid, ms_on, NX, NY, NZ, DX)
         call init_smag_ah(lmix_on, grid, NZ, NU4_BG, resoln_on=.true.)
         call run_smag_ah(grid, ms_on, metrics, lmix_on)
         call check(error, maxval(abs(lmix_on%nu4_face_x - ref_x)) < 1.0e-13_wp .and. &
                    maxval(abs(lmix_on%nu4_face_y - lmix_off%nu4_face_y)) < 1.0e-13_wp, &
                    "resoln_scaled_visc must NOT change the Smag_AH biharmonic nu4")

      end block checks
      call lmix_off%destroy()
      call lmix_on%destroy()
      call ms_off%destroy()
      call ms_on%destroy()
   end subroutine test_smag_ah

   subroutine init_smag_ah(lmix, grid, nz, nu4_bg, resoln_on)
      type(ocean_lateral_mix_t), intent(inout) :: lmix
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: nu4_bg
      logical, intent(in) :: resoln_on
      call lmix%init(grid, nz_ml=nz)
      lmix%smag_ah_active = .true.
      lmix%smag_bi_const = 0.06_wp
      lmix%nu4_bg = nu4_bg
      lmix%nu4_max = 1.0e15_wp
      lmix%resoln_scaled_visc = resoln_on
   end subroutine init_smag_ah

   subroutine run_smag_ah(grid, ms, metrics, lmix)
      !! Map state to the device, run Smag_AH, pull nu4 back, tear the
      !! device attachment down (clean per-call lifecycle).  The biharmonic
      !! kernel takes no res_fn — it is not resolution-scaled (MOM6 parity).
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_metrics_t), intent(inout) :: metrics
      type(ocean_lateral_mix_t), intent(inout) :: lmix
      call make_cartesian_metrics(metrics, grid)
      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(lmix)
      call lmix%enter_data()
      call ocean_lateral_mix_compute_smag_ah(grid, metrics, lmix, ms)
      !$acc update self(lmix%nu4_face_x, lmix%nu4_face_y)
      call lmix%exit_data()
      !$acc exit data delete(lmix)
      call ms%exit_data()
      !$acc exit data delete(ms)
      call destroy_cartesian_metrics(metrics)
   end subroutine run_smag_ah

end module test_ocean_hvisc_resoln
