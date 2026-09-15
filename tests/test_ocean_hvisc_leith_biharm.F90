!! Gap 3 — 2-D Leith biharmonic lateral closure (MOM6 `LEITH_AH`
!! analogue) + the fail-loud closure dispatcher audit.
!!
!! Covers:
!!   1. Zero flow → ζ = 0 → ∇²ζ = 0 → nu4_face = nu4_bg everywhere.
!!   2. Analytic ∇²ζ: a flow giving ζ quadratic in i (∇²ζ constant)
!!      reproduces nu4_face = C_lb · dx⁶ · inv_PI6 · |∇²ζ| to ~1e-8
!!      on both u- and v-faces in the interior.
!!   3. nu4_max clip honoured under a noisy flow.
!!   4. End-to-end: the biharmonic apply reads the Leith-biharm
!!      per-face nu4 (closure=LMIX_LEITH_BIHARM) and produces a
!!      tendency that differs from the scalar nu_4 path.
!!   5. Bit-identity: closure unset (LMIX_NONE) leaves the
!!      horizontal-viscosity tendency on the scalar nu_h path,
!!      bit-for-bit identical to omitting the lateral_mix argument.
!!   6. FAIL-LOUD audit: every namelist tag that `parse_lateral_closure`
!!      accepts maps to an IMPLEMENTED dispatcher code, and a garbage
!!      tag maps to LMIX_INVALID which `lateral_closure_is_implemented`
!!      reports false for — this is exactly what `validate_config`
!!      aborts on, so no tag can silently produce background-only
!!      viscosity.
!!
!! GPU/host pattern mirrors test_ocean_leith / test_ocean_smag_ah.
module test_ocean_hvisc_leith_biharm
   use rdb_constants, only: wp, PI
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_lateral_mix, only: ocean_lateral_mix_t, &
                                    ocean_lateral_mix_compute_leith_biharm, &
                                    parse_lateral_closure, &
                                    lateral_closure_is_implemented, &
                                    lateral_closure_conflicts_smag_ah, &
                                    LMIX_NONE, LMIX_LEITH, LMIX_SMAGORINSKY, &
                                    LMIX_BIHARMONIC, LMIX_LEITH_BIHARM, &
                                    LMIX_INVALID
   use rdb_ocean_horizontal_viscosity, only: ocean_horizontal_viscosity_t, &
                                             ocean_horizontal_viscosity_compute_tendencies
   use testdrive, only: error_type, check, new_unittest, unittest_type
   implicit none
   private

   public :: collect_ocean_hvisc_leith_biharm_tests

contains

   subroutine collect_ocean_hvisc_leith_biharm_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("leith_biharm_zero_flow_is_background", test_zero_flow), &
                  new_unittest("leith_biharm_matches_analytic_del2vort", test_analytic_del2vort), &
                  new_unittest("leith_biharm_nu4_max_clip_honoured", test_nu4_max_clip), &
                  new_unittest("hvisc_with_leith_biharm_uses_face_coefficient", &
                               test_hvisc_with_leith_biharm), &
                  new_unittest("leith_biharm_unset_is_bit_identical", test_bit_identity), &
                  new_unittest("lateral_closure_dispatch_fail_loud", test_fail_loud), &
                  new_unittest("leith_biharm_smag_ah_conflict_fail_loud", test_smag_ah_conflict) &
                  ]
   end subroutine collect_ocean_hvisc_leith_biharm_tests

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
      !! Zero velocity → ζ = 0 → ∇²ζ = 0 → nu4_face = nu4_bg.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_lateral_mix_t) :: lmix
      real(wp), parameter :: NU4_BG = 1.0e8_wp
      real(wp) :: max_diff_x, max_diff_y

      call setup_state(grid, ms, 8, 8, 2, 50000.0_wp)
      call lmix%init(grid, nz_ml=2)
      lmix%closure = LMIX_LEITH_BIHARM
      lmix%c_leith_bi = 8.0_wp
      lmix%nu4_bg = NU4_BG
      lmix%nu4_max = 1.0e12_wp

      call map_in(ms, metrics, grid, lmix)
      call ocean_lateral_mix_compute_leith_biharm(grid, metrics, lmix, ms)
      !$acc update self(lmix%nu4_face_x, lmix%nu4_face_y)
      call map_out(ms, metrics, lmix)

      max_diff_x = maxval(abs(lmix%nu4_face_x - NU4_BG))
      max_diff_y = maxval(abs(lmix%nu4_face_y - NU4_BG))

      call check(error, max_diff_x < 1.0e-9_wp .and. max_diff_y < 1.0e-9_wp, &
                 "Zero flow must leave nu4_face at nu4_bg")

      call lmix%destroy()
      call ms%destroy()
   end subroutine test_zero_flow

   subroutine test_analytic_del2vort(error)
      !! Construct a flow whose corner vorticity ζ is QUADRATIC in i
      !! (uniform in j), so ∇²ζ is a known constant, and check the
      !! per-face nu4 against the closed-form Leith-biharm coefficient.
      !!
      !! With u=0 and v(i,j) = C·i³, on uniform square metrics
      !! (dyCv=dx, iareaBu=1/dx²):
      !!     ζ(i,j) = (v(i)-v(i-1))/dx = C·(3i²-3i+1)/dx
      !! whose discrete x-Laplacian is the second i-difference of a
      !! quadratic, = 6C/dx, divided by dx2q=dx²:
      !!     ∇²ζ = 6C/dx³   (constant, j-independent).
      !! Hence
      !!     nu4_face = C_lb · dx⁶ · inv_PI6 · |∇²ζ|
      !!              = 6·C_lb·inv_PI6·|C|·dx³.
      !! Holds at interior faces away from the zeroed wall corners.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_lateral_mix_t) :: lmix
      integer, parameter :: NX = 14, NY = 14, NZ = 1
      real(wp), parameter :: DX = 40000.0_wp
      real(wp), parameter :: C_LB = 8.0_wp
      real(wp), parameter :: CV = 1.0e-15_wp   ! tiny so v stays physical
      integer, parameter :: IP = 7, JP = 7
      real(wp) :: inv_pi6, expected, obs_x, obs_y
      integer :: i, j

      call setup_state(grid, ms, NX, NY, NZ, DX)
      ! v(i,j) = CV · i³  (constant across j); u stays 0.
      do j = 1, grid%ny_total + 1
         do i = 1, grid%nx_total
            ms%v_face_y_layer(i, j, 1) = CV*real(i, wp)**3
         end do
      end do
      call lmix%init(grid, nz_ml=NZ)
      lmix%closure = LMIX_LEITH_BIHARM
      lmix%c_leith_bi = C_LB
      lmix%nu4_bg = 0.0_wp
      lmix%nu4_max = 1.0e30_wp   ! no clip

      call map_in(ms, metrics, grid, lmix)
      call ocean_lateral_mix_compute_leith_biharm(grid, metrics, lmix, ms)
      !$acc update self(lmix%nu4_face_x, lmix%nu4_face_y)
      call map_out(ms, metrics, lmix)

      inv_pi6 = (1.0_wp/PI)**6
      ! |∇²ζ| = 6·|CV|/dx³ ; nu4 = C_lb·dx⁶·inv_PI6·|∇²ζ|.
      expected = C_LB*(DX**6)*inv_pi6*(6.0_wp*abs(CV)/(DX**3))

      obs_x = lmix%nu4_face_x(IP, JP, 1)
      obs_y = lmix%nu4_face_y(IP, JP, 1)

      call check(error, abs(obs_x - expected) <= 1.0e-8_wp*expected + 1.0e-30_wp, &
                 "nu4_face_x must match analytic Leith-biharm coefficient")
      if (.not. allocated(error)) then
         call check(error, abs(obs_y - expected) <= 1.0e-8_wp*expected + 1.0e-30_wp, &
                    "nu4_face_y must match analytic Leith-biharm coefficient")
      end if
      ! Sanity: the closure actually produced a non-zero coefficient.
      if (.not. allocated(error)) then
         call check(error, expected > 0.0_wp .and. obs_x > 0.0_wp, &
                    "Leith-biharm must produce a non-zero nu4 for a sheared flow")
      end if

      call lmix%destroy()
      call ms%destroy()
   end subroutine test_analytic_del2vort

   subroutine test_nu4_max_clip(error)
      !! Noisy alternating flow → large |∇²ζ| → assert nu4_max cap held.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      type(ocean_lateral_mix_t) :: lmix
      real(wp), parameter :: NU4_MAX = 3.0e8_wp
      integer, parameter :: NX = 12, NY = 12, NZ = 1
      real(wp) :: peak_x, peak_y
      integer :: i, j

      call setup_state(grid, ms, NX, NY, NZ, 40000.0_wp)
      do j = 1, grid%ny_total
         do i = 1, grid%nx_total + 1
            ms%u_face_x_layer(i, j, 1) = merge(1.0_wp, -1.0_wp, mod(i + j, 2) == 0)
         end do
      end do
      do j = 1, grid%ny_total + 1
         do i = 1, grid%nx_total
            ms%v_face_y_layer(i, j, 1) = merge(-1.0_wp, 1.0_wp, mod(i + j, 2) == 0)
         end do
      end do
      call lmix%init(grid, nz_ml=NZ)
      lmix%closure = LMIX_LEITH_BIHARM
      lmix%c_leith_bi = 8.0_wp
      lmix%nu4_bg = 0.0_wp
      lmix%nu4_max = NU4_MAX

      call map_in(ms, metrics, grid, lmix)
      call ocean_lateral_mix_compute_leith_biharm(grid, metrics, lmix, ms)
      !$acc update self(lmix%nu4_face_x, lmix%nu4_face_y)
      call map_out(ms, metrics, lmix)

      peak_x = maxval(lmix%nu4_face_x)
      peak_y = maxval(lmix%nu4_face_y)

      call check(error, peak_x <= NU4_MAX + 1.0e-6_wp .and. peak_y <= NU4_MAX + 1.0e-6_wp, &
                 "nu4_max clip not honoured under noisy flow")

      call lmix%destroy()
      call ms%destroy()
   end subroutine test_nu4_max_clip

   subroutine test_hvisc_with_leith_biharm(error)
      !! End-to-end: hvisc compute_tendencies with closure
      !! LMIX_LEITH_BIHARM engages the per-face biharmonic apply and
      !! produces a tendency that differs from the scalar nu_4 path.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid_face, grid_scalar
      type(multilayer_state_t) :: ms_face, ms_scalar
      type(ocean_metrics_t) :: metrics_face, metrics_scalar
      type(ocean_lateral_mix_t) :: lmix_face
      type(ocean_horizontal_viscosity_t) :: hv_face, hv_scalar
      real(wp), parameter :: DX = 40000.0_wp, NU4_CONST = 1.0e9_wp
      integer, parameter :: NX = 14, NY = 14, NZ = 1
      real(wp) :: max_diff
      integer :: i, j

      call setup_state(grid_face, ms_face, NX, NY, NZ, DX)
      call setup_state(grid_scalar, ms_scalar, NX, NY, NZ, DX)
      do j = 1, grid_face%ny_total + 1
         do i = 1, grid_face%nx_total
            ms_face%v_face_y_layer(i, j, 1) = 0.01_wp*sin(0.7_wp*real(i, wp))
            ms_scalar%v_face_y_layer(i, j, 1) = ms_face%v_face_y_layer(i, j, 1)
         end do
      end do
      call hv_face%init(grid_face, nz_ml=NZ)
      call hv_scalar%init(grid_scalar, nz_ml=NZ)
      hv_face%nu_h = 0.0_wp
      hv_face%nu_4 = 0.0_wp     ! bypassed by the Leith-biharm face path
      hv_scalar%nu_h = 0.0_wp
      hv_scalar%nu_4 = NU4_CONST

      call lmix_face%init(grid_face, nz_ml=NZ)
      lmix_face%closure = LMIX_LEITH_BIHARM
      lmix_face%smag_ah_active = .false.
      lmix_face%c_leith_bi = 8.0_wp
      lmix_face%nu4_bg = 0.0_wp
      lmix_face%nu4_max = 1.0e30_wp

      call map_in(ms_face, metrics_face, grid_face, lmix_face, hv_face)
      call ocean_lateral_mix_compute_leith_biharm(grid_face, metrics_face, lmix_face, ms_face)
      call ocean_horizontal_viscosity_compute_tendencies(grid_face, metrics_face, hv_face, ms_face, &
                                                         lateral_mix=lmix_face)
      !$acc update self(hv_face%du_visc%data, hv_face%dv_visc%data)
      call map_out(ms_face, metrics_face, lmix_face, hv_face)

      call map_in(ms_scalar, metrics_scalar, grid_scalar, hv=hv_scalar)
      call ocean_horizontal_viscosity_compute_tendencies(grid_scalar, metrics_scalar, hv_scalar, ms_scalar)
      !$acc update self(hv_scalar%du_visc%data, hv_scalar%dv_visc%data)
      call map_out(ms_scalar, metrics_scalar, hv=hv_scalar)

      max_diff = maxval(abs(hv_face%dv_visc%data - hv_scalar%dv_visc%data))

      call check(error, max_diff > 1.0e-12_wp, &
                 "Leith-biharm face path must differ from scalar nu_4 tendency")

      call lmix_face%destroy()
      call hv_face%destroy()
      call hv_scalar%destroy()
      call ms_face%destroy()
      call ms_scalar%destroy()
   end subroutine test_hvisc_with_leith_biharm

   subroutine test_bit_identity(error)
      !! Closure unset (LMIX_NONE) with the lateral_mix argument present
      !! must give a tendency bit-for-bit identical to omitting the
      !! argument entirely (scalar nu_h Laplacian only).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid_a, grid_b
      type(multilayer_state_t) :: ms_a, ms_b
      type(ocean_metrics_t) :: metrics_a, metrics_b
      type(ocean_lateral_mix_t) :: lmix_off
      type(ocean_horizontal_viscosity_t) :: hv_a, hv_b
      real(wp), parameter :: DX = 40000.0_wp, NU_H = 500.0_wp
      integer, parameter :: NX = 12, NY = 12, NZ = 1
      real(wp) :: max_diff_u, max_diff_v
      integer :: i, j

      call setup_state(grid_a, ms_a, NX, NY, NZ, DX)
      call setup_state(grid_b, ms_b, NX, NY, NZ, DX)
      do j = 1, grid_a%ny_total
         do i = 1, grid_a%nx_total + 1
            ms_a%u_face_x_layer(i, j, 1) = 0.02_wp*sin(0.5_wp*real(j, wp))
            ms_b%u_face_x_layer(i, j, 1) = ms_a%u_face_x_layer(i, j, 1)
         end do
      end do
      call hv_a%init(grid_a, nz_ml=NZ)
      call hv_b%init(grid_b, nz_ml=NZ)
      hv_a%nu_h = NU_H
      hv_a%nu_4 = 0.0_wp
      hv_b%nu_h = NU_H
      hv_b%nu_4 = 0.0_wp

      call lmix_off%init(grid_a, nz_ml=NZ)
      lmix_off%closure = LMIX_NONE
      lmix_off%smag_ah_active = .false.

      ! Path A: lateral_mix present but closure=NONE.
      call map_in(ms_a, metrics_a, grid_a, lmix_off, hv_a)
      call ocean_horizontal_viscosity_compute_tendencies(grid_a, metrics_a, hv_a, ms_a, &
                                                         lateral_mix=lmix_off)
      !$acc update self(hv_a%du_visc%data, hv_a%dv_visc%data)
      call map_out(ms_a, metrics_a, lmix_off, hv_a)

      ! Path B: no lateral_mix argument.
      call map_in(ms_b, metrics_b, grid_b, hv=hv_b)
      call ocean_horizontal_viscosity_compute_tendencies(grid_b, metrics_b, hv_b, ms_b)
      !$acc update self(hv_b%du_visc%data, hv_b%dv_visc%data)
      call map_out(ms_b, metrics_b, hv=hv_b)

      max_diff_u = maxval(abs(hv_a%du_visc%data - hv_b%du_visc%data))
      max_diff_v = maxval(abs(hv_a%dv_visc%data - hv_b%dv_visc%data))

      call check(error, max_diff_u == 0.0_wp .and. max_diff_v == 0.0_wp, &
                 "closure=NONE must be bit-identical to omitting lateral_mix")

      call lmix_off%destroy()
      call hv_a%destroy()
      call hv_b%destroy()
      call ms_a%destroy()
      call ms_b%destroy()
   end subroutine test_bit_identity

   subroutine test_fail_loud(error)
      !! Audit the closure dispatcher: every accepted namelist string
      !! maps to an IMPLEMENTED code (no silent no-op), and a garbage
      !! string maps to LMIX_INVALID which is NOT implemented — this is
      !! exactly the predicate `validate_config` aborts on.
      type(error_type), allocatable, intent(out) :: error

      ! 1. All documented tags parse to implemented codes.
      call check(error, parse_lateral_closure("none") == LMIX_NONE .and. &
                 lateral_closure_is_implemented(LMIX_NONE), &
                 "'none' must parse + be implemented")
      if (allocated(error)) return
      call check(error, parse_lateral_closure("leith") == LMIX_LEITH .and. &
                 lateral_closure_is_implemented(LMIX_LEITH), &
                 "'leith' must parse + be implemented")
      if (allocated(error)) return
      call check(error, parse_lateral_closure("smagorinsky") == LMIX_SMAGORINSKY .and. &
                 lateral_closure_is_implemented(LMIX_SMAGORINSKY), &
                 "'smagorinsky' must parse + be implemented")
      if (allocated(error)) return
      call check(error, parse_lateral_closure("biharmonic") == LMIX_BIHARMONIC .and. &
                 lateral_closure_is_implemented(LMIX_BIHARMONIC), &
                 "'biharmonic' must parse + be implemented")
      if (allocated(error)) return
      call check(error, parse_lateral_closure("leith_biharm") == LMIX_LEITH_BIHARM .and. &
                 lateral_closure_is_implemented(LMIX_LEITH_BIHARM), &
                 "'leith_biharm' must parse + be implemented")
      if (allocated(error)) return

      ! 2. A garbage tag must NOT silently become a working closure: it
      !    parses to LMIX_INVALID and reports unimplemented, which makes
      !    validate_config fail loud rather than run background-only.
      call check(error, parse_lateral_closure("not_a_closure") == LMIX_INVALID, &
                 "garbage tag must parse to LMIX_INVALID")
      if (allocated(error)) return
      call check(error,.not. lateral_closure_is_implemented(LMIX_INVALID), &
                 "LMIX_INVALID must report unimplemented (drives fail-loud abort)")
      if (allocated(error)) return

      ! 3. Empty string is the explicit 'off' spelling, not garbage.
      call check(error, parse_lateral_closure("") == LMIX_NONE, &
                 "empty string must map to LMIX_NONE (off)")
   end subroutine test_fail_loud

   subroutine test_smag_ah_conflict(error)
      !! `leith_biharm` and `smag_ah` are both flow-aware biharmonic
      !! closures that fill `nu4_face_*`; with both on, `compute_smag_ah`
      !! silently overwrites the Leith-biharmonic fill.  This asserts the
      !! `lateral_closure_conflicts_smag_ah` predicate that `validate_config`
      !! aborts on — exercising the exact code path of the fail-loud guard
      !! (the `error stop` itself can't run in-process).
      type(error_type), allocatable, intent(out) :: error

      ! Only LMIX_LEITH_BIHARM + smag_ah is a conflict.
      call check(error, lateral_closure_conflicts_smag_ah(LMIX_LEITH_BIHARM, .true.), &
                 "leith_biharm + smag_ah must be flagged as a conflict")
      if (allocated(error)) return
      ! Leith-biharm alone is fine.
      call check(error,.not. lateral_closure_conflicts_smag_ah(LMIX_LEITH_BIHARM, .false.), &
                 "leith_biharm without smag_ah must not conflict")
      if (allocated(error)) return
      ! Smag_KH (Laplacian) + smag_ah (biharmonic) fill DIFFERENT arrays
      ! (ah_face vs nu4_face) — the intended combo, not a conflict.
      call check(error,.not. lateral_closure_conflicts_smag_ah(LMIX_SMAGORINSKY, .true.), &
                 "smagorinsky (Laplacian) + smag_ah must not conflict")
      if (allocated(error)) return
      ! No closure + smag_ah is the plain Smag_AH path — not a conflict.
      call check(error,.not. lateral_closure_conflicts_smag_ah(LMIX_NONE, .true.), &
                 "none + smag_ah must not conflict")
      if (allocated(error)) return
      call check(error,.not. lateral_closure_conflicts_smag_ah(LMIX_LEITH, .true.), &
                 "leith (Laplacian) + smag_ah must not conflict")
   end subroutine test_smag_ah_conflict

end module test_ocean_hvisc_leith_biharm
