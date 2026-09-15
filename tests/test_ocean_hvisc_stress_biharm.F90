!! PR-4 — regression + physics gate for the `stress_tensor` early-return
!! fix in `ocean_horizontal_viscosity_compute_tendencies`.
!!
!! Before this PR, `if (this%stress_tensor) then ... return` short-
!! circuited the procedure, so the biharmonic add-on below it (constant
!! `nu_4`, flow-aware `smag_ah`, `leith_biharm`) was dead code whenever
!! `stress_tensor = .true.` — exactly zero on `main`.  This PR converts
!! the `return` into an `else if` so the biharmonic composes with the
!! stress-divergence harmonic operator, matching MOM6 (`BIHARMONIC`
!! "may be used with `LAPLACIAN`").
!!
!! Covers:
!!   1. `stress_plus_nu4_dissipates` — the roadmap's named regression
!!      gate: a grid-scale checkerboard u-mode is measurably damped
!!      under `stress_tensor + nu_4`.  Exactly zero on `main`.
!!   2. `stress_biharm_matches_laplacian_biharm` — the combined
!!      `stress_tensor + nu_h + nu_4` tendency agrees with the
!!      `Laplacian + nu_h + nu_4` tendency to round-off (uniform grid,
!!      uniform h, all-wet, h >> H_VANISHED).
!!   3. `stress_biharm_analytic_k4_decay` — magnitude AND sign: a single
!!      Fourier mode decays at the discrete `-nu_4 * beta_sq**2`
!;      eigenvalue (beta_sq the discrete Laplacian eigenvalue), not just
!!      "non-zero".
!!   4. `stress_plus_smag_ah_engages` — the flow-aware face-biharmonic
!!      arm (Smag_AH) is also reachable under `stress_tensor`, not just
!!      the scalar `nu_4` arm.
!!   5. `stress_biharm_off_is_bit_identical` — an inert biharmonic
!!      config (`nu_4=0`, no flow-aware closure) is bit-for-bit (`==`)
!!      identical whether or not a `lateral_mix` argument is passed —
!!      pins the new dispatch reached via the deleted `return`.
!!   6. `aniso_plus_biharm_compose` — `kh_aniso` (off-axis direction) and
!!      `nu_4` are independent linear operators that superpose exactly
!!      under `stress_tensor` — the exclusion the roadmap asked to be
!!      guarded dissolves; this proves no guard is needed.
module test_ocean_hvisc_stress_biharm
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, PI
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_lateral_mix, only: ocean_lateral_mix_t, &
                                    ocean_lateral_mix_compute_smag_ah, &
                                    LMIX_NONE
   use rdb_ocean_horizontal_viscosity, only: &
      ocean_horizontal_viscosity_t, &
      ocean_horizontal_viscosity_compute_tendencies, &
      ocean_hvisc_set_aniso_direction
   implicit none
   private

   public :: collect_ocean_hvisc_stress_biharm_tests

contains

   subroutine collect_ocean_hvisc_stress_biharm_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("stress_plus_nu4_dissipates", &
                               test_stress_plus_nu4_dissipates), &
                  new_unittest("stress_biharm_matches_laplacian_biharm", &
                               test_stress_biharm_matches_laplacian_biharm), &
                  new_unittest("stress_biharm_analytic_k4_decay", &
                               test_stress_biharm_analytic_k4_decay), &
                  new_unittest("stress_plus_smag_ah_engages", &
                               test_stress_plus_smag_ah_engages), &
                  new_unittest("stress_biharm_off_is_bit_identical", &
                               test_stress_biharm_off_is_bit_identical), &
                  new_unittest("aniso_plus_biharm_compose", &
                               test_aniso_plus_biharm_compose) &
                  ]
   end subroutine collect_ocean_hvisc_stress_biharm_tests

   ! ------------------------------------------------------------------
   ! Shared helpers
   ! ------------------------------------------------------------------

   subroutine setup_state(grid, ms, nx, ny, nz, dx)
      type(hgrid_t), intent(inout) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: dx
      call grid%init(nx, ny, 2, dx, dx)
      ms%nz_ml = nz
      call ms%init(grid)
      ms%h_layer = 1.0e6_wp    ! h >> H_VANISHED: h_neglect floor negligible
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

   ! ------------------------------------------------------------------
   ! Test 1 — the roadmap's named regression gate.
   ! ------------------------------------------------------------------

   subroutine test_stress_plus_nu4_dissipates(error)
      !! Grid-scale checkerboard u-mode `u(i,j) = (-1)^(i+j)` under
      !! `stress_tensor=.true., nu_h=0, nu_4>0`.  Exactly zero on `main`
      !! (the `return` at the top of the stress branch skipped the
      !! biharmonic add-on entirely).  After the fix, the biharmonic
      !! must measurably damp this mode: `max|du_visc| > 0` and the
      !! tendency opposes the mode (`Σ u·du_visc < 0`, KE removed).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_horizontal_viscosity_t) :: hv
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: DX = 1000.0_wp, NU4 = 1.0e10_wp, DT = 1.0e-3_wp
      integer, parameter :: NX = 20, NY = 18, NZ = 1
      integer :: i, j, k, nx_t, ny_t
      real(wp) :: max_abs_du, sum_u_du
      checks: block

         call setup_state(grid, ms, NX, NY, NZ, DX)
         nx_t = grid%nx_total
         ny_t = grid%ny_total
         do k = 1, NZ
            do j = 1, ny_t
               do i = 1, nx_t + 1
                  ms%u_face_x_layer(i, j, k) = merge(1.0_wp, -1.0_wp, mod(i + j, 2) == 0)
               end do
            end do
         end do

         call hv%init(grid, nz_ml=NZ)
         hv%nu_h = 0.0_wp
         hv%nu_4 = NU4
         hv%stress_tensor = .true.

         call map_in(ms, metrics, grid, hv=hv)
         call ocean_horizontal_viscosity_compute_tendencies(grid, metrics, hv, ms, dt=DT)
         !$acc update self(hv%du_visc%data)
         call map_out(ms, metrics, hv=hv)

         ! Interior window well away from every wall (Pass 1/2 mirror BC
         ! at the immediate wall-adjacent rows/columns).
         max_abs_du = 0.0_wp
         sum_u_du = 0.0_wp
         do k = 1, NZ
            do j = 4, ny_t - 3
               do i = 4, nx_t - 2
                  max_abs_du = max(max_abs_du, abs(hv%du_visc%data(i, j, k)))
                  sum_u_du = sum_u_du + ms%u_face_x_layer(i, j, k)*hv%du_visc%data(i, j, k)
               end do
            end do
         end do

         call check(error, max_abs_du > 0.0_wp, &
                    "stress_tensor + nu_4: grid-scale checkerboard mode must be damped "// &
                    "(max|du_visc| == 0 means the biharmonic add-on is still dead code)")
         if (allocated(error)) exit checks
         call check(error, sum_u_du < 0.0_wp, &
                    "stress_tensor + nu_4: dissipation sign wrong — Sum(u*du_visc) must be < 0")

      end block checks
      call hv%destroy(); call ms%destroy()
   end subroutine test_stress_plus_nu4_dissipates

   ! ------------------------------------------------------------------
   ! Test 2 — additive composition matches the Laplacian + nu_4 path.
   ! ------------------------------------------------------------------

   subroutine test_stress_biharm_matches_laplacian_biharm(error)
      !! Two `ocean_horizontal_viscosity_t` slots, same field, uniform
      !! square grid + uniform h + all-wet, both with `nu_h=NU` and
      !! `nu_4=NU4`; one with `stress_tensor=.true.`, one without.
      !! Interior tendencies must agree to round-off — proves the
      !! composition is additive and correctly ordered (the harmonic
      !! reduction is already proven by `test_stress_reduces_to_laplacian`
      !! in `test_ocean_hvisc.F90`; any discrepancy here is a
      !! biharmonic-composition error).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_horizontal_viscosity_t) :: hv_lap, hv_m6
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: L = 1234.5_wp, NU = 7.5e3_wp, NU4 = 1.0e10_wp, DT = 1.0e-3_wp
      integer, parameter :: NX = 14, NY = 11, NZ = 3
      real(wp), allocatable :: du_lap(:, :, :), dv_lap(:, :, :)
      real(wp) :: max_du, max_dv, scale
      integer :: i, j, k, nx_t, ny_t
      checks: block

         call setup_state(grid, ms, NX, NY, NZ, L)
         nx_t = grid%nx_total
         ny_t = grid%ny_total
         do k = 1, NZ
            do j = 1, ny_t
               do i = 1, nx_t + 1
                  ms%u_face_x_layer(i, j, k) = &
                     sin(0.21_wp*real(i, wp))*cos(0.17_wp*real(j, wp)) + 0.3_wp*real(k, wp)
               end do
            end do
            do j = 1, ny_t + 1
               do i = 1, nx_t
                  ms%v_face_y_layer(i, j, k) = &
                     cos(0.13_wp*real(i, wp))*sin(0.29_wp*real(j, wp)) - 0.2_wp*real(k, wp)
               end do
            end do
         end do

         call hv_lap%init(grid, nz_ml=NZ)
         call hv_m6%init(grid, nz_ml=NZ)
         hv_lap%nu_h = NU
         hv_lap%nu_4 = NU4
         hv_m6%nu_h = NU
         hv_m6%nu_4 = NU4
         hv_m6%stress_tensor = .true.
         hv_m6%bound_coef = 0.8_wp

         !$acc enter data copyin(ms)
         call ms%enter_data()
         call make_cartesian_metrics(metrics, grid)
         !$acc enter data copyin(hv_lap)
         call hv_lap%enter_data()
         !$acc enter data copyin(hv_m6)
         call hv_m6%enter_data()

         call ocean_horizontal_viscosity_compute_tendencies(grid, metrics, hv_lap, ms, dt=DT)
         call ocean_horizontal_viscosity_compute_tendencies(grid, metrics, hv_m6, ms, dt=DT)
         !$acc update self(hv_lap%du_visc%data, hv_lap%dv_visc%data)
         !$acc update self(hv_m6%du_visc%data, hv_m6%dv_visc%data)

         allocate (du_lap, source=hv_lap%du_visc%data)
         allocate (dv_lap, source=hv_lap%dv_visc%data)

         call hv_m6%exit_data()
         !$acc exit data delete(hv_m6)
         call hv_lap%exit_data()
         !$acc exit data delete(hv_lap)
         call ms%exit_data()
         !$acc exit data delete(ms)
         call destroy_cartesian_metrics(metrics)

         max_du = 0.0_wp
         max_dv = 0.0_wp
         do k = 1, NZ
            do j = 2, ny_t - 1
               do i = 2, nx_t
                  max_du = max(max_du, abs(hv_m6%du_visc%data(i, j, k) - du_lap(i, j, k)))
               end do
            end do
            do j = 2, ny_t
               do i = 2, nx_t - 1
                  max_dv = max(max_dv, abs(hv_m6%dv_visc%data(i, j, k) - dv_lap(i, j, k)))
               end do
            end do
         end do

         scale = max(maxval(abs(du_lap)), maxval(abs(dv_lap)))
         call check(error, max_du < 1.0e-9_wp*scale, &
                    "stress_tensor + nu_4: u tendency differs from Laplacian + nu_4 "// &
                    "beyond round-off")
         if (allocated(error)) exit checks
         call check(error, max_dv < 1.0e-9_wp*scale, &
                    "stress_tensor + nu_4: v tendency differs from Laplacian + nu_4 "// &
                    "beyond round-off")

      end block checks
      if (allocated(du_lap)) deallocate (du_lap)
      if (allocated(dv_lap)) deallocate (dv_lap)
      call hv_m6%destroy(); call hv_lap%destroy(); call ms%destroy()
   end subroutine test_stress_biharm_matches_laplacian_biharm

   ! ------------------------------------------------------------------
   ! Test 3 — analytic k^4 decay rate, magnitude AND sign.
   ! ------------------------------------------------------------------

   subroutine test_stress_biharm_analytic_k4_decay(error)
      !! `stress_tensor=.true., nu_h=0, nu_4=NU4`, `u = cos(k*x)`, `v=0`.
      !! The harmonic contribution is exactly zero (A ≡ 0 under
      !! `nu_h=0`), so `du_visc` at an interior probe must match the
      !! discrete two-pass metric-Laplacian eigenvalue
      !!   `du_visc ≈ -NU4 * beta_sq**2 * u_ic`,
      !! `beta_sq = 2*(1-cos(k*dx))/dx**2` — the SAME discrete eigenvalue
      !! `test_sinusoid_decay` (`test_ocean_hvisc.F90`) uses for the
      !! single-pass Laplacian, applied twice here.  Fixes the
      !! *magnitude and sign*, not just non-zero-ness (closes audit G13
      !! for the stress path).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_horizontal_viscosity_t) :: hv
      type(ocean_metrics_t) :: metrics
      ! DX=1000 (not 1) so the per-face biharmonic CFL clamp
      ! (hvisc_nu4_cfl_bound, which scales as bound_coef/(dt*(1/dx^2 +
      ! 1/dy^2)^2)) has enormous headroom over NU4 — on a dx=1 grid the
      ! default bound_coef=0.8 clamp caps nu_4 at O(4), far below any
      ! interesting NU4, which silently invalidates the analytic
      ! comparison (confirmed empirically: dx=1 clamped NU4=5e3 down to
      ! ~4, an 1150x reduction).
      real(wp), parameter :: DX = 1000.0_wp, NU4 = 5.0e3_wp, DT = 1.0e-3_wp
      integer, parameter :: NX = 40, NY = 12, NZ = 1
      integer :: i, j, k, nx_t, ny_t, i_probe, j_probe
      real(wp) :: k_wave, beta_sq, ic_at_probe, expect, obs, rel_err
      checks: block

         call setup_state(grid, ms, NX, NY, NZ, DX)
         nx_t = grid%nx_total
         ny_t = grid%ny_total
         k_wave = 2.0_wp*PI/real(nx_t, wp)
         beta_sq = 2.0_wp*(1.0_wp - cos(k_wave*grid%dx))/(grid%dx*grid%dx)

         do k = 1, NZ
            do j = 1, ny_t
               do i = 1, nx_t + 1
                  ms%u_face_x_layer(i, j, k) = cos(k_wave*real(i - 1, wp)*grid%dx)
               end do
            end do
         end do

         call hv%init(grid, nz_ml=NZ)
         hv%nu_h = 0.0_wp
         hv%nu_4 = NU4
         hv%stress_tensor = .true.

         call map_in(ms, metrics, grid, hv=hv)
         call ocean_horizontal_viscosity_compute_tendencies(grid, metrics, hv, ms, dt=DT)
         !$acc update self(hv%du_visc%data)
         call map_out(ms, metrics, hv=hv)

         ! Probe well interior in i (mirror BC contaminates the 2 columns
         ! nearest each x-wall) and interior in j (no j-dependence in the
         ! field, so the mirror BC there is analytically exact anyway).
         i_probe = nx_t/2 + 1
         j_probe = ny_t/2
         ic_at_probe = cos(k_wave*real(i_probe - 1, wp)*grid%dx)
         expect = -NU4*beta_sq*beta_sq*ic_at_probe
         obs = hv%du_visc%data(i_probe, j_probe, NZ)
         rel_err = abs(obs - expect)/abs(expect)

         call check(error, rel_err < 1.0e-6_wp, &
                    "stress_tensor + nu_4: measured decay rate does not match the "// &
                    "discrete -nu_4*beta_sq**2 eigenvalue (magnitude/sign error in "// &
                    "the composition)")

      end block checks
      call hv%destroy(); call ms%destroy()
   end subroutine test_stress_biharm_analytic_k4_decay

   ! ------------------------------------------------------------------
   ! Test 4 — flow-aware face-biharmonic arm (Smag_AH) under stress_tensor.
   ! ------------------------------------------------------------------

   subroutine test_stress_plus_smag_ah_engages(error)
      !! `stress_tensor=.true.` + `lateral_mix%smag_ah_active=.true.`
      !! with `smag_bi_const > 0` on a sheared field must produce a
      !! `du_visc` different from the `smag_ah_active=.false.` run.
      !! Covers the `hvisc_compute_biharmonic_face_impl` arm (`:474` in
      !! the module), a DIFFERENT branch of the dispatch than the
      !! scalar-`nu_4` arm test 1 exercises.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid_on, grid_off
      type(multilayer_state_t) :: ms_on, ms_off
      type(ocean_metrics_t) :: metrics_on, metrics_off
      type(ocean_lateral_mix_t) :: lmix_on
      type(ocean_horizontal_viscosity_t) :: hv_on, hv_off
      real(wp), parameter :: DX = 50000.0_wp, NU_H = 500.0_wp
      integer, parameter :: NX = 12, NY = 12, NZ = 1
      real(wp) :: max_diff
      integer :: i, j

      call setup_state(grid_on, ms_on, NX, NY, NZ, DX)
      call setup_state(grid_off, ms_off, NX, NY, NZ, DX)
      do j = 1, grid_on%ny_total
         do i = 1, grid_on%nx_total + 1
            ms_on%u_face_x_layer(i, j, 1) = 0.01_wp*sin(2.0_wp*real(j, wp))
            ms_off%u_face_x_layer(i, j, 1) = ms_on%u_face_x_layer(i, j, 1)
         end do
      end do

      call hv_on%init(grid_on, nz_ml=NZ)
      call hv_off%init(grid_off, nz_ml=NZ)
      hv_on%nu_h = NU_H
      hv_on%nu_4 = 0.0_wp   ! bypassed by the Smag_AH face path
      hv_on%stress_tensor = .true.
      hv_off%nu_h = NU_H
      hv_off%nu_4 = 0.0_wp
      hv_off%stress_tensor = .true.

      call lmix_on%init(grid_on, nz_ml=NZ)
      lmix_on%closure = LMIX_NONE
      lmix_on%smag_ah_active = .true.
      lmix_on%smag_bi_const = 0.06_wp
      lmix_on%nu4_bg = 0.0_wp
      lmix_on%nu4_max = 1.0e12_wp

      call map_in(ms_on, metrics_on, grid_on, lmix_on, hv_on)
      call ocean_lateral_mix_compute_smag_ah(grid_on, metrics_on, lmix_on, ms_on)
      call ocean_horizontal_viscosity_compute_tendencies(grid_on, metrics_on, hv_on, ms_on, &
                                                         lateral_mix=lmix_on)
      !$acc update self(hv_on%du_visc%data, hv_on%dv_visc%data)
      call map_out(ms_on, metrics_on, lmix_on, hv_on)

      call map_in(ms_off, metrics_off, grid_off, hv=hv_off)
      call ocean_horizontal_viscosity_compute_tendencies(grid_off, metrics_off, hv_off, ms_off)
      !$acc update self(hv_off%du_visc%data, hv_off%dv_visc%data)
      call map_out(ms_off, metrics_off, hv=hv_off)

      max_diff = maxval(abs(hv_on%du_visc%data - hv_off%du_visc%data))

      call check(error, max_diff > 1.0e-12_wp, &
                 "stress_tensor + smag_ah_active: face-biharmonic arm must produce a "// &
                 "tendency different from smag_ah_active=.false.")

      call lmix_on%destroy()
      call hv_on%destroy()
      call hv_off%destroy()
      call ms_on%destroy()
      call ms_off%destroy()
   end subroutine test_stress_plus_smag_ah_engages

   ! ------------------------------------------------------------------
   ! Test 5 — inert-biharmonic bit-identity gate.
   ! ------------------------------------------------------------------

   subroutine test_stress_biharm_off_is_bit_identical(error)
      !! `stress_tensor=.true., nu_4=0`, no flow-aware closure.  Result
      !! must be bit-for-bit (`==`) identical whether the `lateral_mix`
      !! argument is OMITTED entirely, or PRESENT-but-inert
      !! (`is_init=.true.`, `closure=LMIX_NONE`, `smag_ah_active=.false.`)
      !! — both fall through to the (never-entered) scalar `nu_4 > 0`
      !! check.  This exercises the exact dispatch reached via the
      !! deleted `return`: the `if (present(lateral_mix))` / `if
      !! (lateral_mix%is_init .and. ...)` guards on the biharmonic
      !! add-on, now reached from the stress-tensor branch too.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid_a, grid_b
      type(multilayer_state_t) :: ms_a, ms_b
      type(ocean_metrics_t) :: metrics_a, metrics_b
      type(ocean_lateral_mix_t) :: lmix_inert
      type(ocean_horizontal_viscosity_t) :: hv_a, hv_b
      real(wp), parameter :: DX = 40000.0_wp, NU_H = 1000.0_wp
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
      hv_a%stress_tensor = .true.
      hv_b%nu_h = NU_H
      hv_b%nu_4 = 0.0_wp
      hv_b%stress_tensor = .true.

      call lmix_inert%init(grid_a, nz_ml=NZ)
      lmix_inert%closure = LMIX_NONE
      lmix_inert%smag_ah_active = .false.

      ! Path A: stress_tensor, lateral_mix present but inert.
      call map_in(ms_a, metrics_a, grid_a, lmix_inert, hv_a)
      call ocean_horizontal_viscosity_compute_tendencies(grid_a, metrics_a, hv_a, ms_a, &
                                                         lateral_mix=lmix_inert)
      !$acc update self(hv_a%du_visc%data, hv_a%dv_visc%data)
      call map_out(ms_a, metrics_a, lmix_inert, hv_a)

      ! Path B: stress_tensor, no lateral_mix argument at all.
      call map_in(ms_b, metrics_b, grid_b, hv=hv_b)
      call ocean_horizontal_viscosity_compute_tendencies(grid_b, metrics_b, hv_b, ms_b)
      !$acc update self(hv_b%du_visc%data, hv_b%dv_visc%data)
      call map_out(ms_b, metrics_b, hv=hv_b)

      max_diff_u = maxval(abs(hv_a%du_visc%data - hv_b%du_visc%data))
      max_diff_v = maxval(abs(hv_a%dv_visc%data - hv_b%dv_visc%data))

      call check(error, max_diff_u == 0.0_wp .and. max_diff_v == 0.0_wp, &
                 "stress_tensor with an inert biharmonic config must be bit-identical "// &
                 "regardless of whether lateral_mix is present-but-inert or omitted")

      call lmix_inert%destroy()
      call hv_a%destroy()
      call hv_b%destroy()
      call ms_a%destroy()
      call ms_b%destroy()
   end subroutine test_stress_biharm_off_is_bit_identical

   ! ------------------------------------------------------------------
   ! Test 6 — kh_aniso + nu_4 superpose exactly (the dissolved exclusion).
   ! ------------------------------------------------------------------

   subroutine test_aniso_plus_biharm_compose(error)
      !! `stress_tensor=.true.`, `kh_aniso > 0` with an off-axis
      !! `aniso_dir` (cross terms live) PLUS `nu_4 > 0`.  Both are
      !! independent LINEAR operators on `u`/`v`, so the combined
      !! tendency must equal the sum of the `kh_aniso`-only and the
      !! `nu_4`-only tendencies to round-off.  This is the assertion
      !! that would catch a future re-introduction of an early return
      !! (or any other coupling bug between the two operators) — it
      !! replaces the mutual-exclusion GUARD the roadmap asked for with
      !! a proof that no guard is needed.
      !!
      !! `nu_h = 0` throughout (in `run_case` below): isolates the two
      !! add-on operators under test.  With a shared nonzero `nu_h`
      !! baseline present in BOTH the "aniso-only" and "nu4-only" runs,
      !! naive addition would double-count the isotropic Laplacian term
      !! — the superposition principle only holds relative to a common
      !! baseline, and zeroing it out is the simplest way to keep this
      !! test an honest, direct additivity check.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: L = 1234.5_wp, KH_A = 5.0e3_wp
      real(wp), parameter :: NU4 = 1.0e10_wp, DT = 1.0e-3_wp
      real(wp), parameter :: N1 = 0.6_wp, N2 = 0.8_wp   ! off-axis direction
      integer, parameter :: NX = 14, NY = 11, NZ = 1
      real(wp), allocatable :: du_aniso(:, :, :), dv_aniso(:, :, :)
      real(wp), allocatable :: du_nu4(:, :, :), dv_nu4(:, :, :)
      real(wp), allocatable :: du_combo(:, :, :), dv_combo(:, :, :)
      real(wp) :: max_du, max_dv, scale
      integer :: i, j, k, nx_t, ny_t
      checks: block

         call run_case(KH_A, N1, N2, 0.0_wp, du_aniso, dv_aniso, nx_t, ny_t)
         call run_case(0.0_wp, N1, N2, NU4, du_nu4, dv_nu4, nx_t, ny_t)
         call run_case(KH_A, N1, N2, NU4, du_combo, dv_combo, nx_t, ny_t)

         max_du = 0.0_wp
         max_dv = 0.0_wp
         do k = 1, NZ
            do j = 2, ny_t - 1
               do i = 2, nx_t
                  max_du = max(max_du, abs(du_combo(i, j, k) - &
                                           (du_aniso(i, j, k) + du_nu4(i, j, k))))
               end do
            end do
            do j = 2, ny_t
               do i = 2, nx_t - 1
                  max_dv = max(max_dv, abs(dv_combo(i, j, k) - &
                                           (dv_aniso(i, j, k) + dv_nu4(i, j, k))))
               end do
            end do
         end do

         scale = max(maxval(abs(du_combo)), maxval(abs(dv_combo)))
         call check(error, max_du < 1.0e-9_wp*scale, &
                    "kh_aniso + nu_4 under stress_tensor: u tendency does not superpose "// &
                    "to round-off")
         if (allocated(error)) exit checks
         call check(error, max_dv < 1.0e-9_wp*scale, &
                    "kh_aniso + nu_4 under stress_tensor: v tendency does not superpose "// &
                    "to round-off")

      end block checks
      if (allocated(du_aniso)) deallocate (du_aniso, dv_aniso)
      if (allocated(du_nu4)) deallocate (du_nu4, dv_nu4)
      if (allocated(du_combo)) deallocate (du_combo, dv_combo)
   contains
      subroutine run_case(kh_aniso, n1, n2, nu4, du_out, dv_out, nx_t, ny_t)
         real(wp), intent(in) :: kh_aniso, n1, n2, nu4
         real(wp), allocatable, intent(out) :: du_out(:, :, :), dv_out(:, :, :)
         integer, intent(out) :: nx_t, ny_t
         type(hgrid_t) :: g
         type(multilayer_state_t) :: m
         type(ocean_horizontal_viscosity_t) :: hv
         type(ocean_metrics_t) :: met
         integer :: ii, jj, kk

         call setup_state(g, m, NX, NY, NZ, L)
         nx_t = g%nx_total
         ny_t = g%ny_total
         do kk = 1, NZ
            do jj = 1, ny_t
               do ii = 1, nx_t + 1
                  m%u_face_x_layer(ii, jj, kk) = &
                     sin(0.21_wp*real(ii, wp))*cos(0.17_wp*real(jj, wp))
               end do
            end do
            do jj = 1, ny_t + 1
               do ii = 1, nx_t
                  m%v_face_y_layer(ii, jj, kk) = &
                     cos(0.13_wp*real(ii, wp))*sin(0.29_wp*real(jj, wp))
               end do
            end do
         end do

         call hv%init(g, nz_ml=NZ)
         hv%nu_h = 0.0_wp
         hv%nu_4 = nu4
         hv%stress_tensor = .true.
         hv%bound_coef = 0.8_wp
         hv%kh_aniso = kh_aniso
         call ocean_hvisc_set_aniso_direction(hv, n1, n2)

         call map_in(m, met, g, hv=hv)
         call ocean_horizontal_viscosity_compute_tendencies(g, met, hv, m, dt=DT)
         !$acc update self(hv%du_visc%data, hv%dv_visc%data)
         call map_out(m, met, hv=hv)

         allocate (du_out, source=hv%du_visc%data)
         allocate (dv_out, source=hv%dv_visc%data)

         call hv%destroy(); call m%destroy()
      end subroutine run_case
   end subroutine test_aniso_plus_biharm_compose

end module test_ocean_hvisc_stress_biharm
