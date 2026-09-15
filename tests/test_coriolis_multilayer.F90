!! Unit tests for the multilayer Sadourny Coriolis-advection kernel
!! (rdb_coriolis_adv Phase 5b — per-layer lift of the Phase 3b
!! barotropic kernel).  Each k-slice runs the same ζ-at-corners,
!! KE-at-centres, (ζ+f)*v - grad(KE) pipeline; layers are
!! independent (no vertical coupling in this Phase).
!!
!! Cases:
!!   * Per-layer zero-velocity null check — every layer must stay at
!!     zero u, v after a Coriolis step.  Trivial guard against
!!     uninitialized branches injecting spurious tendency.
!!   * Per-layer inertial oscillation — each layer initialised with
!!     a *different* uniform U_k = U0 * k velocity in x, zero in y.
!!     After a quarter inertial period each layer's velocity vector
!!     must have rotated through f*t = π/2 with magnitude preserved
!!     within 2% (the FE tolerance — same as the Phase 3a barotropic
!!     test, but per-layer).  Confirms (a) the per-layer kernel
!!     produces the same dynamics as the barotropic, and (b) the
!!     layers don't bleed into each other.
!!   * Per-layer Sadourny solid-body discriminator — a different
!!     omega per layer; each layer's v-face tendency at an interior
!!     probe must match the discrete Sadourny formula `omega^2 +
!!     f*omega` times the v-face displacement, and must differ from
!!     plain Coriolis by the expected factor.
module test_coriolis_multilayer
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, H_DIV_EPS
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_coriolis_adv, only: coriolis_adv_t, &
                               coriolis_adv_compute_tendencies, &
                               coriolis_adv_compute_tendencies_hk, &
                               coriolis_adv_compute_tendencies_sadourny_energy, &
                               coriolis_adv_apply_tendencies, &
                               PV_VARIANT_SADOURNY_HK, &
                               PV_VARIANT_SADOURNY_ENERGY, &
                               CORNER_H_CELL_MEAN, CORNER_H_MOM6_AREA, &
                               PV_ADV_CENTERED, PV_ADV_WENO3, PV_ADV_WENO5, PV_ADV_WENO7, &
                               PV_ADV_INVALID, weno3_recon, weno5_recon, weno7_recon, &
                               parse_pv_adv_scheme, pv_adv_scheme_is_implemented, &
                               pv_adv_required_nghost
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics, &
                                 make_anisotropic_metrics
   implicit none
   private

   public :: collect_coriolis_multilayer_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 3
   real(wp), parameter :: BC_H_MIN_PV = 1.0e-12_wp
      !! Mirror of the kernel's CORIOLIS_H_MIN_PV (BOUND_CORIOLIS h_corner oracle).

contains

   subroutine collect_coriolis_multilayer_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("per_layer_zero_velocity", test_zero_velocity), &
                  new_unittest("per_layer_inertial_oscillation", &
                               test_inertial_oscillation), &
                  new_unittest("per_layer_sadourny_solid_body", &
                               test_sadourny_solid_body), &
                  new_unittest("beta_plane_f_profile", &
                               test_beta_plane_profile), &
                  new_unittest("beta_plane_rossby_drift", &
                               test_beta_plane_rossby_drift), &
                  new_unittest("variable_h_weights_coriolis", &
                               test_variable_h_weights_coriolis), &
                  new_unittest("hk_uniform_reduction", &
                               test_hk_uniform_reduction), &
                  new_unittest("hk_nonuniform_h_discriminator", &
                               test_hk_nonuniform_h_discriminator), &
                  new_unittest("hk_shear_flow_discriminator", &
                               test_hk_shear_flow_discriminator), &
                  new_unittest("hk_corner_h_area_weighted", &
                               test_hk_corner_h_area_weighted), &
                  new_unittest("energy_uniform_reduction", &
                               test_energy_uniform_reduction), &
                  new_unittest("energy_nonuniform_h_exact", &
                               test_energy_nonuniform_h_exact), &
                  new_unittest("energy_zero_velocity_null", &
                               test_energy_zero_velocity_null), &
                  new_unittest("energy_anisotropic_metrics_exact", &
                               test_energy_anisotropic_metrics_exact), &
                  new_unittest("energy_beta_plane_discriminator", &
                               test_energy_beta_plane_discriminator), &
                  new_unittest("energy_state_fluxes_gate", &
                               test_energy_state_fluxes_gate), &
                  new_unittest("bound_coriolis_allwet_inert", &
                               test_bound_coriolis_allwet_inert), &
                  new_unittest("bound_coriolis_masked_clamps", &
                               test_bound_coriolis_masked_clamps), &
                  new_unittest("corner_h_mom6_equiv_above_floor", &
                               test_corner_h_equiv), &
                  new_unittest("corner_h_mom6_floor_differs", &
                               test_corner_h_floor), &
                  new_unittest("weno3_recon_constant_exact", test_weno3_constant), &
                  new_unittest("weno3_recon_linear_exact", test_weno3_linear), &
                  new_unittest("weno3_recon_step_no_overshoot", test_weno3_step), &
                  new_unittest("weno3_recon_upwind_asymmetric", test_weno3_upwind), &
                  new_unittest("weno3_scheme_parse_and_gate", test_weno3_parse), &
                  new_unittest("weno_uniform_vort_matches_centered", test_weno_uniform), &
                  new_unittest("weno_sharp_vort_differs_from_centered", test_weno_sharp), &
                  new_unittest("weno5_recon_exact_and_eno", test_weno5_recon), &
                  new_unittest("weno7_recon_exact_and_eno", test_weno7_recon), &
                  new_unittest("weno5_kernel_uniform_matches_centered", test_weno5_kernel) &
                  ]
   end subroutine collect_coriolis_multilayer_tests

   subroutine make_grid(grid, nx_phys, ny_phys, dx, dy)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      real(wp), intent(in) :: dx, dy
      call grid%init(nx_phys, ny_phys, NGHOST, dx, dy)
   end subroutine make_grid

   subroutine map_in(ms, cor)
      type(multilayer_state_t), intent(inout) :: ms
      type(coriolis_adv_t), intent(inout) :: cor
      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(cor)
      call cor%enter_data()
   end subroutine map_in

   subroutine map_out(ms, cor)
      type(multilayer_state_t), intent(inout) :: ms
      type(coriolis_adv_t), intent(inout) :: cor
      call cor%exit_data()
      !$acc exit data delete(cor)
      call ms%exit_data()
      !$acc exit data delete(ms)
   end subroutine map_out

   ! -----------------------------------------------------------------
   ! Cases
   ! -----------------------------------------------------------------

   subroutine test_zero_velocity(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(coriolis_adv_t) :: cor
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: F_C = 1.0_wp
      real(wp), parameter :: DT = 0.01_wp
      real(wp) :: max_u, max_v
      checks: block

         call make_grid(grid, 8, 8, 1.0_wp, 1.0_wp)
         call make_cartesian_metrics(metrics, grid)
         ms%nz_ml = NZ
         call ms%init(grid)
         cor%f_0 = F_C
         call cor%init(grid, nz_ml=NZ)
         ms%h_layer = 1.0_wp
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp

         call map_in(ms, cor)
         call coriolis_adv_compute_tendencies(grid, metrics, cor, ms)
         call coriolis_adv_apply_tendencies(cor, ms, DT)
         call map_out(ms, cor)

         max_u = maxval(abs(ms%u_face_x_layer))
         max_v = maxval(abs(ms%v_face_y_layer))
         call check(error, max_u < 1.0e-14_wp, &
                    "multilayer Coriolis injected u from zero state")
         if (allocated(error)) exit checks
         call check(error, max_v < 1.0e-14_wp, &
                    "multilayer Coriolis injected v from zero state")

      end block checks
      call cor%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_zero_velocity

   subroutine test_inertial_oscillation(error)
      !! Per-layer inertial oscillation with a layer-dependent
      !! initial u amplitude.  Each layer must rotate at angular
      !! rate f (independent of u amplitude) and preserve its
      !! magnitude within the FE tolerance.  Tests independence of
      !! layers as well as correctness per layer.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(coriolis_adv_t) :: cor
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: F_C = 1.0_wp
      real(wp), parameter :: U0 = 1.0_wp
      real(wp), parameter :: DT = 0.01_wp
      integer, parameter :: N_STEPS = 157   ! ~ pi/(2*F_C*DT)
      real(wp) :: u_k, t_end, u_expected, v_expected, mag_expected
      real(wp) :: u_obs, v_obs, mag_obs, mag_err
      integer :: step, k, i_probe, j_probe, nx, ny
      checks: block

         call make_grid(grid, 16, 16, 1.0_wp, 1.0_wp)
         call make_cartesian_metrics(metrics, grid)
         ms%nz_ml = NZ
         call ms%init(grid)
         cor%f_0 = F_C
         call cor%init(grid, nz_ml=NZ)
         ms%h_layer = 1.0_wp
         nx = grid%nx_total
         ny = grid%ny_total

         ! Layer-dependent IC: U_k = U0 * k
         ms%v_face_y_layer = 0.0_wp
         do k = 1, NZ
            u_k = U0*real(k, wp)
            ms%u_face_x_layer(:, :, k) = u_k
         end do

         call map_in(ms, cor)
         do step = 1, N_STEPS
            call coriolis_adv_compute_tendencies(grid, metrics, cor, ms)
            call coriolis_adv_apply_tendencies(cor, ms, DT)
         end do
         call map_out(ms, cor)

         i_probe = nx/2
         j_probe = ny/2
         t_end = real(N_STEPS, wp)*DT

         do k = 1, NZ
            u_k = U0*real(k, wp)
            u_obs = ms%u_face_x_layer(i_probe, j_probe, k)
            v_obs = ms%v_face_y_layer(i_probe, j_probe, k)
            u_expected = u_k*cos(F_C*t_end)
            v_expected = -u_k*sin(F_C*t_end)
            mag_obs = sqrt(u_obs*u_obs + v_obs*v_obs)
            mag_expected = sqrt(u_expected*u_expected + v_expected*v_expected)
            mag_err = abs(mag_obs - mag_expected)/mag_expected

            call check(error, mag_err < 0.02_wp, &
                       "layer inertial-oscillation magnitude error > 2%")
            if (allocated(error)) exit checks
         end do

      end block checks
      call cor%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_inertial_oscillation

   subroutine test_sadourny_solid_body(error)
      !! Layer-dependent solid-body rotation: omega(k) = OMEGA0 * k.
      !! Each layer's expected discrete Sadourny dv/dt at the probe
      !! is (omega(k)^2 + f*omega(k)) * r_y; plain Coriolis would
      !! give omega(k) * r_y.  Catches cross-layer scratch reuse
      !! errors: a kernel that mistakenly read q_corner / ke_centre
      !! from the wrong layer would still satisfy layer 1 (where
      !! omega = OMEGA0) but fail layers 2 and 3.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(coriolis_adv_t) :: cor
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: OMEGA0 = 1.0_wp
      real(wp), parameter :: F_C = 1.0_wp
      real(wp), parameter :: DT = 0.01_wp
      integer, parameter :: I_PROBE_OFFSET = 3
      integer, parameter :: J_PROBE_OFFSET = 5
      integer :: i, j, k, nx, ny, i_c, j_c, i_probe, j_probe
      real(wp) :: omega_k, r_y
      real(wp) :: v_initial, v_final, dv_dt_obs
      real(wp) :: dv_dt_expected_sadourny, dv_dt_expected_plain
      checks: block

         call make_grid(grid, 16, 16, 1.0_wp, 1.0_wp)
         call make_cartesian_metrics(metrics, grid)
         ms%nz_ml = NZ
         call ms%init(grid)
         cor%f_0 = F_C
         call cor%init(grid, nz_ml=NZ)
         ms%h_layer = 1.0_wp
         nx = grid%nx_total
         ny = grid%ny_total
         i_c = nx/2
         j_c = ny/2

         ! Solid-body rotation per layer with omega(k) = OMEGA0 * k
         do k = 1, NZ
            omega_k = OMEGA0*real(k, wp)
            do j = 1, ny
               do i = 1, nx + 1
                  ms%u_face_x_layer(i, j, k) = -omega_k*(real(j, wp) - real(j_c, wp))
               end do
            end do
            do j = 1, ny + 1
               do i = 1, nx
                  ms%v_face_y_layer(i, j, k) = omega_k*(real(i, wp) - real(i_c, wp))
               end do
            end do
         end do

         i_probe = i_c + I_PROBE_OFFSET
         j_probe = j_c + J_PROBE_OFFSET
         r_y = real(j_probe, wp) - real(j_c, wp) - 0.5_wp

         call map_in(ms, cor)
         call coriolis_adv_compute_tendencies(grid, metrics, cor, ms)
         call coriolis_adv_apply_tendencies(cor, ms, DT)
         call map_out(ms, cor)

         do k = 1, NZ
            omega_k = OMEGA0*real(k, wp)
            v_initial = omega_k*(real(i_probe, wp) - real(i_c, wp))
            v_final = ms%v_face_y_layer(i_probe, j_probe, k)
            dv_dt_obs = (v_final - v_initial)/DT

            dv_dt_expected_sadourny = r_y*(omega_k**2 + F_C*omega_k)
            dv_dt_expected_plain = omega_k*r_y

            call check(error, abs(dv_dt_obs - dv_dt_expected_sadourny) < 1.0e-9_wp, &
                       "layer Sadourny dv/dt deviates from analytic")
            if (allocated(error)) exit checks
            call check(error, &
                       abs(dv_dt_obs - dv_dt_expected_plain) > 0.5_wp* &
                       abs(dv_dt_expected_sadourny - dv_dt_expected_plain), &
                       "layer dv/dt matched plain Coriolis (Sadourny missing)")
            if (allocated(error)) exit checks
         end do

      end block checks
      call cor%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_sadourny_solid_body

   subroutine test_beta_plane_profile(error)
      !! Verify `set_beta_plane` populates `f_corner` with the analytic
      !! linear-in-y profile `f(y) = f_0 + beta*(y - y_ref)`.  This is
      !! the algebra-level check: every interior + boundary corner must
      !! match the closed-form expression at machine precision.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(coriolis_adv_t) :: cor
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: F0 = 1.0e-4_wp     ! 30°N-ish
      real(wp), parameter :: BETA = 2.0e-11_wp  ! 1/(s·m)
      real(wp), parameter :: DX = 5.0e3_wp, DY = 5.0e3_wp
      real(wp), parameter :: TOL = 1.0e-14_wp
      integer, parameter :: NX_PHYS = 8, NY_PHYS = 6
      real(wp) :: y_ref, y, f_expected, max_err
      integer :: i, j
      checks: block

         call make_grid(grid, NX_PHYS, NY_PHYS, DX, DY)
         call make_cartesian_metrics(metrics, grid)
         call cor%init(grid, nz_ml=NZ)

         y_ref = 0.5_wp*real(NY_PHYS, wp)*DY
         call cor%set_beta_plane(grid, F0, BETA, y_ref)

         ! f_corner sits at C-grid corners (i, j) = (1..nx+1, 1..ny+1).
         ! Corner j has y = (j - 1 - nghost)*dy by `set_beta_plane`'s convention.
         max_err = 0.0_wp
         do j = 1, grid%ny_total + 1
            y = real(j - 1 - grid%nghost, wp)*grid%dy
            f_expected = F0 + BETA*(y - y_ref)
            do i = 1, grid%nx_total + 1
               max_err = max(max_err, abs(cor%f_corner(i, j) - f_expected))
            end do
         end do

         call check(error, max_err < TOL, &
                    "f_corner did not match analytic beta-plane profile")
         if (allocated(error)) exit checks
         call check(error, abs(cor%f_0 - F0) < TOL, "f_0 not stored")
         if (allocated(error)) exit checks
         call check(error, abs(cor%beta - BETA) < TOL, "beta not stored")

      end block checks
      call cor%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_beta_plane_profile

   subroutine test_beta_plane_rossby_drift(error)
      !! Sanity check that a non-zero beta produces a different per-face
      !! tendency than an f-plane with the same f_0.  Initial state is
      !! uniform u = U0 over a non-trivial grid; with beta > 0 the
      !! Coriolis parameter varies linearly across y, so the v-tendency
      !! (proportional to -f*u) must vary monotonically in y.  An
      !! f-plane control run produces uniform v-tendency.  This is a
      !! discriminator, not a quantitative comparison — confirms the
      !! refactor wired the corner field through the kernel rather than
      !! silently falling back to a constant.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(coriolis_adv_t) :: cor_beta, cor_fplane
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: F0 = 1.0e-4_wp
      real(wp), parameter :: BETA = 1.0e-8_wp   ! exaggerated to make drift visible
      real(wp), parameter :: U0 = 1.0_wp
      real(wp), parameter :: DX = 1.0e3_wp, DY = 1.0e3_wp
      real(wp) :: y_ref, dv_dt_south, dv_dt_north, dv_dt_fplane
      integer :: i_probe, j_south, j_north
      checks: block

         call make_grid(grid, 16, 16, DX, DY)
         call make_cartesian_metrics(metrics, grid)
         ms%nz_ml = NZ
         call ms%init(grid)
         y_ref = 0.5_wp*real(grid%ny_total, wp)*DY

         ! ---- Beta-plane run ----
         call cor_beta%init(grid, nz_ml=NZ)
         call cor_beta%set_beta_plane(grid, F0, BETA, y_ref)
         ms%h_layer = 1.0_wp
         ms%u_face_x_layer = U0
         ms%v_face_y_layer = 0.0_wp
         !$acc enter data copyin(ms)
         call ms%enter_data()
         !$acc enter data copyin(cor_beta)
         call cor_beta%enter_data()
         call coriolis_adv_compute_tendencies(grid, metrics, cor_beta, ms)
         !$acc update self(cor_beta%pv_flux_y%data)
         call cor_beta%exit_data()
         !$acc exit data delete(cor_beta)
         call ms%exit_data()
         !$acc exit data delete(ms)

         i_probe = grid%nx_total/2
         j_south = grid%nghost + 2
         j_north = grid%ny_total - grid%nghost - 1
         dv_dt_south = cor_beta%pv_flux_y%data(i_probe, j_south, 1)
         dv_dt_north = cor_beta%pv_flux_y%data(i_probe, j_north, 1)

         ! With u = U0 > 0 and f increasing northward, -f*u becomes more
         ! negative northward; the v-tendency at the north probe must be
         ! strictly less than at the south probe.
         call check(error, dv_dt_north < dv_dt_south, &
                    "beta-plane dv/dt failed to vary monotonically with y")
         if (allocated(error)) exit checks

         ! ---- f-plane control ----
         cor_fplane%f_0 = F0
         call cor_fplane%init(grid, nz_ml=NZ)
         !$acc enter data copyin(ms)
         call ms%enter_data()
         !$acc enter data copyin(cor_fplane)
         call cor_fplane%enter_data()
         call coriolis_adv_compute_tendencies(grid, metrics, cor_fplane, ms)
         !$acc update self(cor_fplane%pv_flux_y%data)
         call cor_fplane%exit_data()
         !$acc exit data delete(cor_fplane)
         call ms%exit_data()
         !$acc exit data delete(ms)

         dv_dt_fplane = cor_fplane%pv_flux_y%data(i_probe, j_south, 1)
         ! f-plane uniform: north and south probes must agree.
         call check(error, &
                    abs(cor_fplane%pv_flux_y%data(i_probe, j_north, 1) - dv_dt_fplane) < 1.0e-15_wp, &
                    "f-plane control showed y-variation (cor%init not zeroing field properly)")

      end block checks
      call cor_beta%destroy()
      call cor_fplane%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_beta_plane_rossby_drift

   subroutine test_variable_h_weights_coriolis(error)
      !! Thickness-weighted Coriolis discriminator.  At a v-face the
      !! kernel computes `u_at_v = Σ(u·h_face) / Σ(h_face)` over the
      !! four neighbouring u-faces, instead of the plain 4-point
      !! average.  Two runs with identical u, f, on the same grid but
      !! different `h_layer` patterns must therefore produce
      !! different dv/dt where u is non-uniform.  Setup: u varies
      !! across i (left half = 0, right half = U0), v = 0.  Run A
      !! uses uniform h; run B uses h heavily weighted toward the
      !! left (low-u) columns.  At an interior v-face straddling
      !! the u-jump, run B must shift `u_at_v` toward the lower-u
      !! side, so |dv/dt| in run B differs from run A.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(coriolis_adv_t) :: cor_A, cor_B
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: F0 = 1.0e-4_wp
      real(wp), parameter :: U0 = 1.0_wp
      real(wp), parameter :: DX = 1.0e3_wp, DY = 1.0e3_wp
      integer, parameter :: NX_PHYS = 8, NY_PHYS = 8
      real(wp) :: dv_dt_uniform, dv_dt_varh, abs_diff
      integer :: i, j, k, nx, ny, i_jump, j_probe

      call make_grid(grid, NX_PHYS, NY_PHYS, DX, DY)
      call make_cartesian_metrics(metrics, grid)
      ms%nz_ml = NZ
      call ms%init(grid)
      nx = grid%nx_total
      ny = grid%ny_total

      ! u: step jump across i = i_jump.  Left half = 0, right = U0.
      i_jump = nx/2
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx + 1
               if (i > i_jump) then
                  ms%u_face_x_layer(i, j, k) = U0
               else
                  ms%u_face_x_layer(i, j, k) = 0.0_wp
               end if
            end do
         end do
      end do
      ms%v_face_y_layer = 0.0_wp

      ! ---- Run A: uniform h ----
      cor_A%f_0 = F0
      call cor_A%init(grid, nz_ml=NZ)
      ms%h_layer = 1.0_wp
      !$acc enter data copyin(ms, cor_A)
      call ms%enter_data(); call cor_A%enter_data()
      call coriolis_adv_compute_tendencies(grid, metrics, cor_A, ms)
      !$acc update self(cor_A%pv_flux_y%data)
      call cor_A%exit_data(); call ms%exit_data()
      !$acc exit data delete(ms, cor_A)

      ! ---- Run B: h heavily weighted toward i <= i_jump ----
      cor_B%f_0 = F0
      call cor_B%init(grid, nz_ml=NZ)
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx
               if (i <= i_jump) then
                  ms%h_layer(i, j, k) = 10.0_wp
               else
                  ms%h_layer(i, j, k) = 1.0_wp
               end if
            end do
         end do
      end do
      !$acc enter data copyin(ms, cor_B)
      call ms%enter_data(); call cor_B%enter_data()
      call coriolis_adv_compute_tendencies(grid, metrics, cor_B, ms)
      !$acc update self(cor_B%pv_flux_y%data)
      call cor_B%exit_data(); call ms%exit_data()
      !$acc exit data delete(ms, cor_B)

      ! Probe a v-face straddling the u-jump.  i = i_jump puts the
      ! east-west neighbours one on each side of the discontinuity.
      j_probe = ny/2
      dv_dt_uniform = cor_A%pv_flux_y%data(i_jump, j_probe, 1)
      dv_dt_varh = cor_B%pv_flux_y%data(i_jump, j_probe, 1)
      abs_diff = abs(dv_dt_varh - dv_dt_uniform)

      ! Discriminator: the two must differ by at least 10% of the
      ! uniform-h magnitude.  (Thickness-weighting toward the u=0
      ! side reduces |u_at_v|, hence |dv/dt|.)
      call check(error, abs_diff > 0.10_wp*abs(dv_dt_uniform), &
                 "thickness-weighted Coriolis did not differ from uniform-h")

      call cor_A%destroy(); call cor_B%destroy(); call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_variable_h_weights_coriolis

   subroutine test_hk_uniform_reduction(error)
      !! Algebraic reduction property: under uniform `h_layer` and
      !! uniform face velocities the Arakawa-Hsu (1990) kernel must
      !! produce the same per-face tendencies as the Sadourny kernel
      !! to round-off.  Catches sign errors, weight typos, indexing
      !! flips in the HK stencil — the only way the two can agree at
      !! 1e-13 is if the algebra matches.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(coriolis_adv_t) :: cor_sad, cor_hk
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: F0 = 1.0e-4_wp
      real(wp), parameter :: U0 = 0.5_wp
      real(wp), parameter :: V0 = 0.3_wp
      real(wp), parameter :: DX = 1.0e3_wp, DY = 1.0e3_wp
      real(wp), parameter :: TOL = 1.0e-13_wp
      integer, parameter :: NX_PHYS = 8, NY_PHYS = 8
      real(wp) :: max_diff_u, max_diff_v

      call make_grid(grid, NX_PHYS, NY_PHYS, DX, DY)
      call make_cartesian_metrics(metrics, grid)
      ms%nz_ml = NZ
      call ms%init(grid)
      ms%h_layer = 1.0_wp
      ms%u_face_x_layer = U0
      ms%v_face_y_layer = V0

      ! ---- Sadourny run ----
      cor_sad%f_0 = F0
      call cor_sad%init(grid, nz_ml=NZ)
      !$acc enter data copyin(ms, cor_sad)
      call ms%enter_data(); call cor_sad%enter_data()
      call coriolis_adv_compute_tendencies(grid, metrics, cor_sad, ms)
      !$acc update self(cor_sad%pv_flux_x%data, cor_sad%pv_flux_y%data)
      call cor_sad%exit_data(); call ms%exit_data()
      !$acc exit data delete(ms, cor_sad)

      ! ---- HK run on the same state ----
      cor_hk%f_0 = F0
      cor_hk%pv_variant = PV_VARIANT_SADOURNY_HK
      call cor_hk%init(grid, nz_ml=NZ)
      !$acc enter data copyin(ms, cor_hk)
      call ms%enter_data(); call cor_hk%enter_data()
      call coriolis_adv_compute_tendencies_hk(grid, metrics, cor_hk, ms, &
                                              ms%u_face_x_layer, ms%v_face_y_layer, ms%h_layer)
      !$acc update self(cor_hk%pv_flux_x%data, cor_hk%pv_flux_y%data)
      call cor_hk%exit_data(); call ms%exit_data()
      !$acc exit data delete(ms, cor_hk)

      max_diff_u = maxval(abs(cor_hk%pv_flux_x%data - cor_sad%pv_flux_x%data))
      max_diff_v = maxval(abs(cor_hk%pv_flux_y%data - cor_sad%pv_flux_y%data))

      call check(error, max_diff_u < TOL, &
                 "HK and Sadourny u-tendencies differ under uniform conditions")
      if (allocated(error)) goto 99
      call check(error, max_diff_v < TOL, &
                 "HK and Sadourny v-tendencies differ under uniform conditions")

99    call cor_sad%destroy(); call cor_hk%destroy(); call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_hk_uniform_reduction

   subroutine test_hk_nonuniform_h_discriminator(error)
      !! Discriminator: under non-uniform `h_layer` the HK kernel reads
      !! per-mass PV q = (f + ζ) / h_at_corner — which varies across
      !! corners when h does — while Sadourny averages (ζ+f) at the
      !! face without normalising by h.  With ζ ≡ 0 and h(i, j)
      !! varying, both kernels see the same Coriolis force per unit
      !! mass but the HK weighting against mass-flux is 4-corner
      !! sensitive, so the tendencies must differ at faces straddling
      !! the h-jump.  Catches a regression where the dispatcher
      !! silently falls back to Sadourny.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(coriolis_adv_t) :: cor_sad, cor_hk
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: F0 = 1.0e-4_wp
      real(wp), parameter :: V0 = 0.5_wp
      real(wp), parameter :: DX = 1.0e3_wp, DY = 1.0e3_wp
      real(wp), parameter :: DISCRIM_TOL = 1.0e-8_wp
      integer, parameter :: NX_PHYS = 8, NY_PHYS = 8
      real(wp) :: max_abs_diff
      integer :: i, j, k, nx, ny

      call make_grid(grid, NX_PHYS, NY_PHYS, DX, DY)
      call make_cartesian_metrics(metrics, grid)
      ms%nz_ml = NZ
      call ms%init(grid)
      nx = grid%nx_total
      ny = grid%ny_total

      ! Non-uniform h: checkerboard pattern (high curvature in space).
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx
               if (mod(i + j, 2) == 0) then
                  ms%h_layer(i, j, k) = 2.0_wp
               else
                  ms%h_layer(i, j, k) = 1.0_wp
               end if
            end do
         end do
      end do
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = V0  ! uniform v, so ζ stays zero

      cor_sad%f_0 = F0
      call cor_sad%init(grid, nz_ml=NZ)
      !$acc enter data copyin(ms, cor_sad)
      call ms%enter_data(); call cor_sad%enter_data()
      call coriolis_adv_compute_tendencies(grid, metrics, cor_sad, ms)
      !$acc update self(cor_sad%pv_flux_x%data)
      call cor_sad%exit_data(); call ms%exit_data()
      !$acc exit data delete(ms, cor_sad)

      cor_hk%f_0 = F0
      cor_hk%pv_variant = PV_VARIANT_SADOURNY_HK
      call cor_hk%init(grid, nz_ml=NZ)
      !$acc enter data copyin(ms, cor_hk)
      call ms%enter_data(); call cor_hk%enter_data()
      call coriolis_adv_compute_tendencies_hk(grid, metrics, cor_hk, ms, &
                                              ms%u_face_x_layer, ms%v_face_y_layer, ms%h_layer)
      !$acc update self(cor_hk%pv_flux_x%data)
      call cor_hk%exit_data(); call ms%exit_data()
      !$acc exit data delete(ms, cor_hk)

      max_abs_diff = maxval(abs(cor_hk%pv_flux_x%data - cor_sad%pv_flux_x%data))
      call check(error, max_abs_diff > DISCRIM_TOL, &
                 "HK and Sadourny u-tendencies agreed under non-uniform h "// &
                 "(dispatcher fell back to Sadourny, or HK kernel ignores h)")

      call cor_sad%destroy(); call cor_hk%destroy(); call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_hk_nonuniform_h_discriminator

   subroutine test_hk_shear_flow_discriminator(error)
      !! Discriminator on a non-trivial vorticity field.  Under
      !! sinusoidal u(j) the relative vorticity ζ = -∂u/∂y is
      !! sinusoidal too, so ζ varies non-linearly across the 6 corners
      !! the HK stencil samples vs. the 2 corners Sadourny samples.
      !! The tendencies must differ wherever the ζ field has curvature.
      !! This exercises the part of the HK stencil that's responsible
      !! for suppressing the Hollingsworth instability at eddy-
      !! resolving resolutions — the bit-identical test would silently
      !! pass even if that branch were wrong, this one catches it.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(coriolis_adv_t) :: cor_sad, cor_hk
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: F0 = 1.0e-4_wp
      real(wp), parameter :: U_AMP = 1.0_wp
      real(wp), parameter :: DX = 1.0e3_wp, DY = 1.0e3_wp
      real(wp), parameter :: DISCRIM_TOL = 1.0e-12_wp
      real(wp), parameter :: PI = 3.141592653589793_wp
      integer, parameter :: NX_PHYS = 16, NY_PHYS = 16
      real(wp) :: max_abs_diff
      integer :: i, j, k, nx, ny
      real(wp) :: u_pattern

      call make_grid(grid, NX_PHYS, NY_PHYS, DX, DY)
      call make_cartesian_metrics(metrics, grid)
      ms%nz_ml = NZ
      call ms%init(grid)
      nx = grid%nx_total
      ny = grid%ny_total

      ! Sinusoidal u(i, j) → non-uniform ζ at corners.  Uniform h
      ! isolates the vorticity-stencil difference from any h
      ! weighting.
      do k = 1, NZ
         do j = 1, ny
            u_pattern = U_AMP*sin(2.0_wp*PI*real(j, wp)/real(ny, wp))
            do i = 1, nx + 1
               ms%u_face_x_layer(i, j, k) = u_pattern
            end do
         end do
      end do
      ms%v_face_y_layer = 0.0_wp
      ms%h_layer = 1.0_wp

      cor_sad%f_0 = F0
      call cor_sad%init(grid, nz_ml=NZ)
      !$acc enter data copyin(ms, cor_sad)
      call ms%enter_data(); call cor_sad%enter_data()
      call coriolis_adv_compute_tendencies(grid, metrics, cor_sad, ms)
      !$acc update self(cor_sad%pv_flux_y%data)
      call cor_sad%exit_data(); call ms%exit_data()
      !$acc exit data delete(ms, cor_sad)

      cor_hk%f_0 = F0
      cor_hk%pv_variant = PV_VARIANT_SADOURNY_HK
      call cor_hk%init(grid, nz_ml=NZ)
      !$acc enter data copyin(ms, cor_hk)
      call ms%enter_data(); call cor_hk%enter_data()
      call coriolis_adv_compute_tendencies_hk(grid, metrics, cor_hk, ms, &
                                              ms%u_face_x_layer, ms%v_face_y_layer, ms%h_layer)
      !$acc update self(cor_hk%pv_flux_y%data)
      call cor_hk%exit_data(); call ms%exit_data()
      !$acc exit data delete(ms, cor_hk)

      max_abs_diff = maxval(abs(cor_hk%pv_flux_y%data - cor_sad%pv_flux_y%data))
      call check(error, max_abs_diff > DISCRIM_TOL, &
                 "HK and Sadourny v-tendencies agreed under curved ζ field "// &
                 "(HK 6-corner stencil collapsed to Sadourny 2-corner)")

      call cor_sad%destroy(); call cor_hk%destroy(); call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_hk_shear_flow_discriminator

   subroutine test_hk_corner_h_area_weighted(error)
      !! Deliverable-1 gate: the HK PV corner thickness h_q is the
      !! AREA-WEIGHTED 4-cell mean
      !!   h_q = Σ areaT·h / Σ areaT
      !! NOT the plain arithmetic mean.  On anisotropic metrics (areaT
      !! varies with j) the two differ; we hand-compute the expected
      !! area-weighted h_q at a few interior corners and assert the
      !! kernel matches it AND differs from the plain mean.
      !!
      !! Probe trick: with ζ ≡ 0 (uniform v) and uniform f_corner = F0,
      !! the kernel writes q_corner = F0 / h_q.  So h_q = F0 / q_corner —
      !! recovered exactly from the kernel output, no internal access
      !! needed.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(coriolis_adv_t) :: cor_hk
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: F0 = 1.0e-4_wp
      real(wp), parameter :: V0 = 0.5_wp
      real(wp), parameter :: DX0 = 1.0e3_wp, DY = 1.0e3_wp, AMP = 0.5_wp
      integer, parameter :: NX_PHYS = 8, NY_PHYS = 8
      integer, parameter :: KPROBE = 2
      integer :: i, j, k, nx, ny, ng, iw, ie, js, jn
      real(wp) :: aSW, aSE, aNW, aNE
      real(wp) :: hSW, hSE, hNW, hNE
      real(wp) :: h_q_aw, h_q_plain, h_q_kernel
      logical :: differs

      call make_grid(grid, NX_PHYS, NY_PHYS, DX0, DY)
      ! Large amplitude (50%) so the area weighting is clearly visible.
      call make_anisotropic_metrics(metrics, grid, DX0, DY, AMP)
      ms%nz_ml = NZ
      call ms%init(grid)
      nx = grid%nx_total
      ny = grid%ny_total
      ng = grid%nghost

      ! Spatially-varying h (strictly positive), uniform v ⇒ ζ ≡ 0.
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx
               ms%h_layer(i, j, k) = 100.0_wp + 10.0_wp*real(i, wp) + &
                                     7.0_wp*real(j, wp) + 3.0_wp*real(k, wp)
            end do
         end do
      end do
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = V0   ! uniform ⇒ relative vorticity stays zero

      cor_hk%f_0 = F0
      cor_hk%pv_variant = PV_VARIANT_SADOURNY_HK
      call cor_hk%init(grid, nz_ml=NZ)
      ! f_corner stays uniform = F0 (init seeds it; no beta-plane call).
      !$acc enter data copyin(ms, cor_hk)
      call ms%enter_data(); call cor_hk%enter_data()
      call coriolis_adv_compute_tendencies_hk(grid, metrics, cor_hk, ms, &
                                              ms%u_face_x_layer, ms%v_face_y_layer, ms%h_layer)
      !$acc update self(cor_hk%q_corner%data)
      call cor_hk%exit_data(); call ms%exit_data()
      !$acc exit data delete(ms, cor_hk)

      ! Hand-check 3 interior corners against the area-weighted formula,
      ! and verify the kernel value differs from the plain 4-cell mean.
      differs = .false.
      k = KPROBE
      do j = ng + 2, ng + 4
         do i = ng + 2, ng + 4
            iw = max(1, i - 1); ie = min(nx, i)
            js = max(1, j - 1); jn = min(ny, j)
            aSW = metrics%areaT(iw, js); aSE = metrics%areaT(ie, js)
            aNW = metrics%areaT(iw, jn); aNE = metrics%areaT(ie, jn)
            hSW = ms%h_layer(iw, js, k); hSE = ms%h_layer(ie, js, k)
            hNW = ms%h_layer(iw, jn, k); hNE = ms%h_layer(ie, jn, k)
            h_q_aw = (aSW*hSW + aSE*hSE + aNW*hNW + aNE*hNE)/ &
                     (aSW + aSE + aNW + aNE)
            h_q_plain = 0.25_wp*(hSW + hSE + hNW + hNE)
            ! Recover the kernel's h_q from q_corner = F0 / h_q.
            h_q_kernel = F0/cor_hk%q_corner%data(i, j, k)

            call check(error, abs(h_q_kernel - h_q_aw) < 1.0e-9_wp*h_q_aw, &
                       "HK corner h_q /= area-weighted mean at interior corner")
            if (allocated(error)) exit
            if (abs(h_q_aw - h_q_plain) > 1.0e-6_wp*h_q_aw) differs = .true.
         end do
         if (allocated(error)) exit
      end do
      if (.not. allocated(error)) then
         call check(error, differs, &
                    "area-weighted h_q identical to plain mean — "// &
                    "test setup does not exercise the area weighting")
      end if

      call cor_hk%destroy(); call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_hk_corner_h_area_weighted

   subroutine test_energy_uniform_reduction(error)
      !! Reduction property (spec §5.1): under uniform `h_layer` and
      !! uniform face velocities the faithful SADOURNY75_ENERGY transport
      !! form (q·vh) must reproduce the enstrophy/velocity form `(f+ζ)·v`
      !! to round-off — the h_v/h_q = 1 collapse the Python prototype shows
      !! at diff = 0.0.  Catches sign/weight/index errors in the q·vh
      !! stencil (the only way to agree at 1e-13 is matching algebra).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(coriolis_adv_t) :: cor_sad, cor_en
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: F0 = 1.0e-4_wp
      real(wp), parameter :: U0 = 0.5_wp, V0 = 0.3_wp
      real(wp), parameter :: DX = 1.0e3_wp, DY = 1.0e3_wp
      real(wp), parameter :: TOL = 1.0e-13_wp
      integer, parameter :: NX_PHYS = 8, NY_PHYS = 8
      real(wp) :: max_diff_u, max_diff_v

      call make_grid(grid, NX_PHYS, NY_PHYS, DX, DY)
      call make_cartesian_metrics(metrics, grid)
      ms%nz_ml = NZ
      call ms%init(grid)
      ms%h_layer = 1.0_wp
      ms%u_face_x_layer = U0
      ms%v_face_y_layer = V0

      ! ---- Sadourny enstrophy run (default) ----
      cor_sad%f_0 = F0
      call cor_sad%init(grid, nz_ml=NZ)
      !$acc enter data copyin(ms, cor_sad)
      call ms%enter_data(); call cor_sad%enter_data()
      call coriolis_adv_compute_tendencies(grid, metrics, cor_sad, ms)
      !$acc update self(cor_sad%pv_flux_x%data, cor_sad%pv_flux_y%data)
      call cor_sad%exit_data(); call ms%exit_data()
      !$acc exit data delete(ms, cor_sad)

      ! ---- Energy run on the same state ----
      cor_en%f_0 = F0
      cor_en%pv_variant = PV_VARIANT_SADOURNY_ENERGY
      call cor_en%init(grid, nz_ml=NZ)
      !$acc enter data copyin(ms, cor_en)
      call ms%enter_data(); call cor_en%enter_data()
      call coriolis_adv_compute_tendencies_sadourny_energy(grid, metrics, cor_en, ms, &
                                                           ms%u_face_x_layer, ms%v_face_y_layer, ms%h_layer)
      !$acc update self(cor_en%pv_flux_x%data, cor_en%pv_flux_y%data)
      call cor_en%exit_data(); call ms%exit_data()
      !$acc exit data delete(ms, cor_en)

      max_diff_u = maxval(abs(cor_en%pv_flux_x%data - cor_sad%pv_flux_x%data))
      max_diff_v = maxval(abs(cor_en%pv_flux_y%data - cor_sad%pv_flux_y%data))

      call check(error, max_diff_u < TOL, &
                 "energy and enstrophy u-tendencies differ under uniform conditions")
      if (allocated(error)) goto 99
      call check(error, max_diff_v < TOL, &
                 "energy and enstrophy v-tendencies differ under uniform conditions")

99    call cor_sad%destroy(); call cor_en%destroy(); call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_energy_uniform_reduction

   subroutine test_energy_state_fluxes_gate(error)
      !! Gate test for the mass-consistent CorAdCalc path
      !! (`use_state_fluxes`): (a) EQUIVALENCE — with
      !! `ms%mass_flux_*_layer` filled to exactly the kernel's own
      !! `u·h_face·dy_cu` recompute, the state-fluxes path must
      !! reproduce the default path to round-off (proves the
      !! convention/shape/wall assumptions match); (b) ENGAGEMENT —
      !! doubling the state fluxes must change the PV part of the
      !! tendency (proves the flag actually reroutes the transport
      !! source and doesn't silently fall back to the recompute).
      !! Sheared v over a thickness step (the nonuniform-h regime where
      !! the transport form genuinely consumes vh).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(coriolis_adv_t) :: cor_en
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: F0 = 1.0e-4_wp, V0 = 0.5_wp
      real(wp), parameter :: DX = 1.0e3_wp, DY = 1.0e3_wp
      real(wp), parameter :: TOL = 1.0e-13_wp
      integer, parameter :: NX_PHYS = 8, NY_PHYS = 8
      integer :: i, j, k, nx, ny, nu, nv
      real(wp) :: h_face, max_diff_u, max_diff_v, engage_diff
      real(wp), allocatable :: ref_x(:, :, :), ref_y(:, :, :)

      call make_grid(grid, NX_PHYS, NY_PHYS, DX, DY)
      call make_cartesian_metrics(metrics, grid)
      ms%nz_ml = NZ
      call ms%init(grid)
      nx = grid%nx_total
      ny = grid%ny_total
      nu = size(ms%u_face_x_layer, 1)
      nv = size(ms%v_face_y_layer, 2)

      ! Sheared v (varying in j) over a thickness step in j — the regime
      ! where q·vh differs from the velocity form.
      ms%u_face_x_layer = 0.0_wp
      do k = 1, NZ
         do j = 1, nv
            ms%v_face_y_layer(:, j, k) = V0*(1.0_wp + 0.1_wp*real(j, wp))
         end do
         do j = 1, ny
            if (j <= ny/2) then
               ms%h_layer(:, j, k) = 10.0_wp
            else
               ms%h_layer(:, j, k) = 40.0_wp
            end if
         end do
      end do

      ! Fill the state fluxes with EXACTLY the recompute's convention.
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nu
               if (i == 1) then
                  h_face = ms%h_layer(1, j, k)
               else if (i == nu) then
                  h_face = ms%h_layer(nx, j, k)
               else
                  h_face = 0.5_wp*(ms%h_layer(i - 1, j, k) + ms%h_layer(i, j, k))
               end if
               ms%mass_flux_x_layer(i, j, k) = &
                  ms%u_face_x_layer(i, j, k)*h_face*metrics%dy_cu(i, j)
            end do
         end do
         do j = 1, nv
            do i = 1, nx
               if (j == 1) then
                  h_face = ms%h_layer(i, 1, k)
               else if (j == nv) then
                  h_face = ms%h_layer(i, ny, k)
               else
                  h_face = 0.5_wp*(ms%h_layer(i, j - 1, k) + ms%h_layer(i, j, k))
               end if
               ms%mass_flux_y_layer(i, j, k) = &
                  ms%v_face_y_layer(i, j, k)*h_face*metrics%dx_cv(i, j)
            end do
         end do
      end do

      cor_en%f_0 = F0
      cor_en%pv_variant = PV_VARIANT_SADOURNY_ENERGY
      call cor_en%init(grid, nz_ml=NZ)
      !$acc enter data copyin(ms, cor_en)
      call ms%enter_data(); call cor_en%enter_data()
      ! mem:separate: mass_flux_*_layer may be mapped `create` (the
      ! production step recomputes them on-device) — push the host
      ! values we just set.
      !$acc update device(ms%mass_flux_x_layer, ms%mass_flux_y_layer)

      ! (a) Reference: default recompute path.
      call coriolis_adv_compute_tendencies_sadourny_energy(grid, metrics, cor_en, ms, &
                                                           ms%u_face_x_layer, ms%v_face_y_layer, ms%h_layer)
      !$acc update self(cor_en%pv_flux_x%data, cor_en%pv_flux_y%data)
      ref_x = cor_en%pv_flux_x%data
      ref_y = cor_en%pv_flux_y%data

      ! (a) Equivalence: state fluxes == recompute values ⇒ identical.
      call coriolis_adv_compute_tendencies_sadourny_energy(grid, metrics, cor_en, ms, &
                                                           ms%u_face_x_layer, ms%v_face_y_layer, ms%h_layer, &
                                                           use_state_fluxes=.true.)
      !$acc update self(cor_en%pv_flux_x%data, cor_en%pv_flux_y%data)
      max_diff_u = maxval(abs(cor_en%pv_flux_x%data - ref_x))
      max_diff_v = maxval(abs(cor_en%pv_flux_y%data - ref_y))
      call check(error, max_diff_u < TOL .and. max_diff_v < TOL, &
                 "state-fluxes path with matching fluxes must reproduce the recompute")
      if (allocated(error)) goto 99

      ! (b) Engagement: doubled vh must change the u-tendency.
      ms%mass_flux_y_layer = 2.0_wp*ms%mass_flux_y_layer
      !$acc update device(ms%mass_flux_y_layer)
      call coriolis_adv_compute_tendencies_sadourny_energy(grid, metrics, cor_en, ms, &
                                                           ms%u_face_x_layer, ms%v_face_y_layer, ms%h_layer, &
                                                           use_state_fluxes=.true.)
      !$acc update self(cor_en%pv_flux_x%data)
      engage_diff = maxval(abs(cor_en%pv_flux_x%data - ref_x))
      call check(error, engage_diff > 1.0e-8_wp, &
                 "doubled state fluxes did not change the tendency — "// &
                 "flag is not rerouting the transport source")

99    call cor_en%exit_data(); call ms%exit_data()
      !$acc exit data delete(ms, cor_en)
      call cor_en%destroy(); call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_energy_state_fluxes_gate

   subroutine test_energy_nonuniform_h_exact(error)
      !! Exact analytic lock (spec §5.2) of the faithful SADOURNY75_ENERGY
      !! u-tendency, on uniform Cartesian metrics, in the regime where the
      !! transport form genuinely differs from the velocity form: SHEARED
      !! v (varying in j) over a strong thickness STEP in j.  With u ≡ 0 and
      !! v uniform in i, relative vorticity ζ ≡ 0 and the centre KE varies
      !! only in j, so ∂x KE = 0; the kernel reduces to the pure q·vh
      !! transport flux, hand-computed from `h_layer` and `v_face_y_layer`:
      !!   q(i,j)  = F0 / h_corner(i,j),   h_corner = ¼Σ(4 cell h)
      !!   vh(i,j) = v(i,j)·½(h(i,j-1)+h(i,j))·dx_cv
      !!   CAu(i,j)= ¼·( q_N·(vh_NW+vh_NE) + q_S·(vh_SW+vh_SE) )·idxCu
      !! Asserts (a) the kernel matches this oracle to ~1e-9 relative and
      !! (b) it DIFFERS from the velocity form `F0·v_at_u` (thickness-
      !! weighted), proving the transport weighting changed the answer.
      !! NOTE: under UNIFORM velocity the two forms are identically equal
      !! for any h (Σ face-h = 2·h_corner telescopes) — the shear is what
      !! makes this a real discriminator.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(coriolis_adv_t) :: cor_en
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: F0 = 1.0e-4_wp, V0 = 0.5_wp
      real(wp), parameter :: DX = 1.0e3_wp, DY = 1.0e3_wp
      integer, parameter :: NX_PHYS = 8, NY_PHYS = 8, KP = 2
      integer :: i, j, k, nx, ny, ip, jp
      real(wp) :: qN, qS, hcN, hcS, vS, vN
      real(wp) :: vh_SW, vh_SE, vh_NW, vh_NE
      real(wp) :: ca_oracle, ca_kernel, ca_velform

      call make_grid(grid, NX_PHYS, NY_PHYS, DX, DY)
      call make_cartesian_metrics(metrics, grid)
      ms%nz_ml = NZ
      call ms%init(grid)
      nx = grid%nx_total
      ny = grid%ny_total
      ip = nx/2
      jp = ny/2

      ! Thickness STEP in j (thin south of jp, thick north) ⇒ h_corner at
      ! the south corner ≪ the north corner, so transport weighting bites.
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx
               if (j > jp) then
                  ms%h_layer(i, j, k) = 100.0_wp + 3.0_wp*real(k, wp)
               else
                  ms%h_layer(i, j, k) = 10.0_wp + 3.0_wp*real(k, wp)
               end if
            end do
         end do
      end do
      ms%u_face_x_layer = 0.0_wp
      ! Sheared v: varies in j (⇒ different face transports), uniform in i
      ! (⇒ ζ ≡ 0).  KE varies only in j ⇒ ∂x KE = 0.
      do k = 1, NZ
         do j = 1, ny + 1
            do i = 1, nx
               ms%v_face_y_layer(i, j, k) = V0*real(j, wp)
            end do
         end do
      end do

      cor_en%f_0 = F0
      cor_en%pv_variant = PV_VARIANT_SADOURNY_ENERGY
      call cor_en%init(grid, nz_ml=NZ)
      !$acc enter data copyin(ms, cor_en)
      call ms%enter_data(); call cor_en%enter_data()
      call coriolis_adv_compute_tendencies_sadourny_energy(grid, metrics, cor_en, ms, &
                                                           ms%u_face_x_layer, ms%v_face_y_layer, ms%h_layer)
      !$acc update self(cor_en%pv_flux_x%data)
      call cor_en%exit_data(); call ms%exit_data()
      !$acc exit data delete(ms, cor_en)

      k = KP
      ! Corner PVs: q = F0 / h_corner (uniform area ⇒ plain 4-cell mean).
      hcS = 0.25_wp*(ms%h_layer(ip - 1, jp - 1, k) + ms%h_layer(ip, jp - 1, k) + &
                     ms%h_layer(ip - 1, jp, k) + ms%h_layer(ip, jp, k))
      hcN = 0.25_wp*(ms%h_layer(ip - 1, jp, k) + ms%h_layer(ip, jp, k) + &
                     ms%h_layer(ip - 1, jp + 1, k) + ms%h_layer(ip, jp + 1, k))
      qS = F0/hcS
      qN = F0/hcN
      ! v-transports vh(i,j) = v(i,j)·½(h(i,j-1)+h(i,j))·dx_cv ; dx_cv = DX.
      vh_SW = ms%v_face_y_layer(ip - 1, jp, k)*0.5_wp* &
              (ms%h_layer(ip - 1, jp - 1, k) + ms%h_layer(ip - 1, jp, k))*DX
      vh_SE = ms%v_face_y_layer(ip, jp, k)*0.5_wp* &
              (ms%h_layer(ip, jp - 1, k) + ms%h_layer(ip, jp, k))*DX
      vh_NW = ms%v_face_y_layer(ip - 1, jp + 1, k)*0.5_wp* &
              (ms%h_layer(ip - 1, jp, k) + ms%h_layer(ip - 1, jp + 1, k))*DX
      vh_NE = ms%v_face_y_layer(ip, jp + 1, k)*0.5_wp* &
              (ms%h_layer(ip, jp, k) + ms%h_layer(ip, jp + 1, k))*DX
      ! CAu = ¼( q_N·(vh_NW+vh_NE) + q_S·(vh_SW+vh_SE) )·idxCu ; idxCu = 1/DX.
      ca_oracle = 0.25_wp*(qN*(vh_NW + vh_NE) + qS*(vh_SW + vh_SE))/DX
      ca_kernel = cor_en%pv_flux_x%data(ip, jp, k)
      ! Velocity form: F0·v_at_u with v_at_u = thickness-weighted mean.
      ! Σ(face-h) = 2·h_corner per row ⇒ v_at_u = (vS·hcS + vN·hcN)/(hcS+hcN).
      vS = ms%v_face_y_layer(ip, jp, k)
      vN = ms%v_face_y_layer(ip, jp + 1, k)
      ca_velform = F0*(vS*hcS + vN*hcN)/(hcS + hcN)

      call check(error, abs(ca_kernel - ca_oracle) < 1.0e-9_wp*abs(ca_oracle), &
                 "energy-form CAu /= analytic q·vh oracle")
      if (allocated(error)) goto 99
      call check(error, abs(ca_oracle - ca_velform) > 1.0e-3_wp*abs(ca_velform), &
                 "energy-form CAu indistinguishable from velocity form "// &
                 "(variable-h transport weighting not exercised)")

99    call cor_en%destroy(); call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_energy_nonuniform_h_exact

   subroutine test_energy_zero_velocity_null(error)
      !! Null check for the energy variant via the public dispatcher
      !! (on-device): zero velocity ⇒ q = f/h_corner finite but every mass
      !! transport vanishes, so CAu/CAv ≡ 0 — no spurious tendency.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(coriolis_adv_t) :: cor
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: F_C = 1.0e-4_wp
      real(wp), parameter :: DT = 0.01_wp
      real(wp) :: max_u, max_v
      checks: block

         call make_grid(grid, 8, 8, 1.0e3_wp, 1.0e3_wp)
         call make_cartesian_metrics(metrics, grid)
         ms%nz_ml = NZ
         call ms%init(grid)
         cor%f_0 = F_C
         cor%pv_variant = PV_VARIANT_SADOURNY_ENERGY
         call cor%init(grid, nz_ml=NZ)
         ms%h_layer = 1.0_wp
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp

         call map_in(ms, cor)
         call coriolis_adv_compute_tendencies(grid, metrics, cor, ms)
         call coriolis_adv_apply_tendencies(cor, ms, DT)
         call map_out(ms, cor)

         max_u = maxval(abs(ms%u_face_x_layer))
         max_v = maxval(abs(ms%v_face_y_layer))
         call check(error, max_u < 1.0e-14_wp, &
                    "energy Coriolis injected u from zero state")
         if (allocated(error)) exit checks
         call check(error, max_v < 1.0e-14_wp, &
                    "energy Coriolis injected v from zero state")

      end block checks
      call cor%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_energy_zero_velocity_null

   subroutine test_energy_anisotropic_metrics_exact(error)
      !! Closes the "uniform-Cartesian-metrics only" coverage gap for the
      !! energy variant's NOVEL stencil (Pass 5/6): on ANISOTROPIC metrics
      !! (`dxT/areaT/dx_cv/idxCu` vary in j, uniform in i) the q·vh form
      !! exercises the area-weighted `h_corner`, the `dx_cv` factor in the
      !! mass transport, AND the `idxCu` factor that closes the flux — none
      !! of which vary on uniform Cartesian.  With u ≡ 0 and v uniform in i
      !! (= V0·j), ζ ≡ 0 (dyCv uniform in i) and KE is uniform in i so
      !! ∂x KE = 0; the kernel reduces to the pure q·vh transport flux.
      !! Asserts (a) the kernel matches a hand-computed oracle built from the
      !! ACTUAL metric arrays to ~1e-9 relative, and (b) that oracle DIFFERS
      !! from one built with a plain (un-area-weighted) corner-h mean — so
      !! the area weighting is genuinely exercised, not incidentally equal.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(coriolis_adv_t) :: cor_en
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: F0 = 1.0e-4_wp, V0 = 0.5_wp
      real(wp), parameter :: DX0 = 1.0e3_wp, DY = 1.0e3_wp, AMP = 0.5_wp
      integer, parameter :: NX_PHYS = 8, NY_PHYS = 8, KP = 2
      integer :: i, j, k, nx, ny, ip, jp
      real(wp) :: qN, qS, hcN, hcS, hcN_pl, hcS_pl
      real(wp) :: vh_SW, vh_SE, vh_NW, vh_NE, vhN, vhS
      real(wp) :: ca_oracle, ca_kernel, ca_plain

      call make_grid(grid, NX_PHYS, NY_PHYS, DX0, DY)
      call make_anisotropic_metrics(metrics, grid, DX0, DY, AMP)
      ms%nz_ml = NZ
      call ms%init(grid)
      nx = grid%nx_total
      ny = grid%ny_total
      ip = nx/2
      jp = ny/2

      ! Thickness step in j so the corner-h area weighting (areaT varies in
      ! j) genuinely differs from a plain 4-cell mean.
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx
               if (j > jp) then
                  ms%h_layer(i, j, k) = 100.0_wp + 3.0_wp*real(k, wp)
               else
                  ms%h_layer(i, j, k) = 10.0_wp + 3.0_wp*real(k, wp)
               end if
            end do
         end do
      end do
      ms%u_face_x_layer = 0.0_wp
      do k = 1, NZ
         do j = 1, ny + 1
            do i = 1, nx
               ms%v_face_y_layer(i, j, k) = V0*real(j, wp)   ! sheared in j, uniform in i
            end do
         end do
      end do

      cor_en%f_0 = F0
      cor_en%pv_variant = PV_VARIANT_SADOURNY_ENERGY
      call cor_en%init(grid, nz_ml=NZ)
      !$acc enter data copyin(ms, cor_en)
      call ms%enter_data(); call cor_en%enter_data()
      call coriolis_adv_compute_tendencies_sadourny_energy(grid, metrics, cor_en, ms, &
                                                           ms%u_face_x_layer, ms%v_face_y_layer, ms%h_layer)
      !$acc update self(cor_en%pv_flux_x%data)
      call cor_en%exit_data(); call ms%exit_data()
      !$acc exit data delete(ms, cor_en)

      k = KP
      ! Area-weighted corner h: q = F0 / [Σ(areaT·h)/Σ(areaT)] over 4 cells.
      hcS = (metrics%areaT(ip - 1, jp - 1)*ms%h_layer(ip - 1, jp - 1, k) + &
             metrics%areaT(ip, jp - 1)*ms%h_layer(ip, jp - 1, k) + &
             metrics%areaT(ip - 1, jp)*ms%h_layer(ip - 1, jp, k) + &
             metrics%areaT(ip, jp)*ms%h_layer(ip, jp, k))/ &
            (metrics%areaT(ip - 1, jp - 1) + metrics%areaT(ip, jp - 1) + &
             metrics%areaT(ip - 1, jp) + metrics%areaT(ip, jp))
      hcN = (metrics%areaT(ip - 1, jp)*ms%h_layer(ip - 1, jp, k) + &
             metrics%areaT(ip, jp)*ms%h_layer(ip, jp, k) + &
             metrics%areaT(ip - 1, jp + 1)*ms%h_layer(ip - 1, jp + 1, k) + &
             metrics%areaT(ip, jp + 1)*ms%h_layer(ip, jp + 1, k))/ &
            (metrics%areaT(ip - 1, jp) + metrics%areaT(ip, jp) + &
             metrics%areaT(ip - 1, jp + 1) + metrics%areaT(ip, jp + 1))
      qS = F0/hcS
      qN = F0/hcN
      ! v-transports vh(i,j) = v(i,j)·½(h(i,j-1)+h(i,j))·dx_cv(i,j).
      vh_SW = ms%v_face_y_layer(ip - 1, jp, k)*0.5_wp* &
              (ms%h_layer(ip - 1, jp - 1, k) + ms%h_layer(ip - 1, jp, k))*metrics%dx_cv(ip - 1, jp)
      vh_SE = ms%v_face_y_layer(ip, jp, k)*0.5_wp* &
              (ms%h_layer(ip, jp - 1, k) + ms%h_layer(ip, jp, k))*metrics%dx_cv(ip, jp)
      vh_NW = ms%v_face_y_layer(ip - 1, jp + 1, k)*0.5_wp* &
              (ms%h_layer(ip - 1, jp, k) + ms%h_layer(ip - 1, jp + 1, k))*metrics%dx_cv(ip - 1, jp + 1)
      vh_NE = ms%v_face_y_layer(ip, jp + 1, k)*0.5_wp* &
              (ms%h_layer(ip, jp, k) + ms%h_layer(ip, jp + 1, k))*metrics%dx_cv(ip, jp + 1)
      ! CAu = ¼( q_N·(vh_NW+vh_NE) + q_S·(vh_SW+vh_SE) )·idxCu(ip,jp).
      ca_oracle = 0.25_wp*(qN*(vh_NW + vh_NE) + qS*(vh_SW + vh_SE))*metrics%idxCu(ip, jp)
      ca_kernel = cor_en%pv_flux_x%data(ip, jp, k)
      ! Plain-mean corner-h control (NOT area-weighted) ⇒ different q ⇒
      ! different CAu, proving the area weighting is exercised.
      hcS_pl = 0.25_wp*(ms%h_layer(ip - 1, jp - 1, k) + ms%h_layer(ip, jp - 1, k) + &
                        ms%h_layer(ip - 1, jp, k) + ms%h_layer(ip, jp, k))
      hcN_pl = 0.25_wp*(ms%h_layer(ip - 1, jp, k) + ms%h_layer(ip, jp, k) + &
                        ms%h_layer(ip - 1, jp + 1, k) + ms%h_layer(ip, jp + 1, k))
      vhN = vh_NW + vh_NE
      vhS = vh_SW + vh_SE
      ca_plain = 0.25_wp*((F0/hcN_pl)*vhN + (F0/hcS_pl)*vhS)*metrics%idxCu(ip, jp)

      call check(error, abs(ca_kernel - ca_oracle) < 1.0e-9_wp*abs(ca_oracle), &
                 "energy-form CAu /= anisotropic-metric q·vh oracle "// &
                 "(area-weighted h_corner / dx_cv / idxCu mishandled)")
      if (allocated(error)) goto 99
      call check(error, abs(ca_oracle - ca_plain) > 1.0e-3_wp*abs(ca_oracle), &
                 "area-weighted h_corner identical to plain mean — "// &
                 "anisotropic weighting not exercised by this setup")

99    call cor_en%destroy(); call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_energy_anisotropic_metrics_exact

   subroutine test_energy_beta_plane_discriminator(error)
      !! Closes the "f-plane only" coverage gap for the energy variant: with
      !! a β-plane f(y) and uniform u = U0 (v = 0, uniform h), the energy
      !! v-tendency `CAv ≈ −q·uh ∝ −f(y)·U0` must vary monotonically in y
      !! (f increases northward).  An f-plane control with the same f_0 gives
      !! a y-uniform tendency.  Confirms the varying planetary f flows through
      !! the q·vh stencil (not silently dropped / constant-folded).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(coriolis_adv_t) :: cor_beta, cor_fplane
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: F0 = 1.0e-4_wp
      real(wp), parameter :: BETA = 1.0e-8_wp   ! exaggerated so drift is visible
      real(wp), parameter :: U0 = 1.0_wp
      real(wp), parameter :: DX = 1.0e3_wp, DY = 1.0e3_wp
      real(wp) :: y_ref, dv_dt_south, dv_dt_north, dv_dt_fp_s, dv_dt_fp_n
      integer :: i_probe, j_south, j_north
      checks: block

         call make_grid(grid, 16, 16, DX, DY)
         call make_cartesian_metrics(metrics, grid)
         ms%nz_ml = NZ
         call ms%init(grid)
         y_ref = 0.5_wp*real(grid%ny_total, wp)*DY
         ms%h_layer = 1.0_wp
         ms%u_face_x_layer = U0
         ms%v_face_y_layer = 0.0_wp

         ! ---- β-plane energy run ----
         cor_beta%pv_variant = PV_VARIANT_SADOURNY_ENERGY
         call cor_beta%init(grid, nz_ml=NZ)
         call cor_beta%set_beta_plane(grid, F0, BETA, y_ref)
         !$acc enter data copyin(ms, cor_beta)
         call ms%enter_data(); call cor_beta%enter_data()
         call coriolis_adv_compute_tendencies(grid, metrics, cor_beta, ms)
         !$acc update self(cor_beta%pv_flux_y%data)
         call cor_beta%exit_data(); call ms%exit_data()
         !$acc exit data delete(ms, cor_beta)

         i_probe = grid%nx_total/2
         j_south = grid%nghost + 2
         j_north = grid%ny_total - grid%nghost - 1
         dv_dt_south = cor_beta%pv_flux_y%data(i_probe, j_south, 1)
         dv_dt_north = cor_beta%pv_flux_y%data(i_probe, j_north, 1)
         ! u>0, f increasing northward ⇒ −f·u more negative northward.
         call check(error, dv_dt_north < dv_dt_south, &
                    "energy β-plane dv/dt failed to vary monotonically with y "// &
                    "(planetary f not flowing through q·vh)")
         if (allocated(error)) exit checks

         ! ---- f-plane control (same f_0) ----
         cor_fplane%f_0 = F0
         cor_fplane%pv_variant = PV_VARIANT_SADOURNY_ENERGY
         call cor_fplane%init(grid, nz_ml=NZ)
         !$acc enter data copyin(ms, cor_fplane)
         call ms%enter_data(); call cor_fplane%enter_data()
         call coriolis_adv_compute_tendencies(grid, metrics, cor_fplane, ms)
         !$acc update self(cor_fplane%pv_flux_y%data)
         call cor_fplane%exit_data(); call ms%exit_data()
         !$acc exit data delete(ms, cor_fplane)

         dv_dt_fp_s = cor_fplane%pv_flux_y%data(i_probe, j_south, 1)
         dv_dt_fp_n = cor_fplane%pv_flux_y%data(i_probe, j_north, 1)
         call check(error, abs(dv_dt_fp_n - dv_dt_fp_s) < 1.0e-15_wp, &
                    "energy f-plane control showed y-variation")

      end block checks
      call cor_beta%destroy(); call cor_fplane%destroy(); call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_energy_beta_plane_discriminator

   ! -----------------------------------------------------------------
   ! BOUND_CORIOLIS (MOM6 :896-909) — velocity-form clamp on the energy PV flux
   ! -----------------------------------------------------------------

   pure function recon_h_corner(ic, jc, k, nx, ny, h, wet_T, areaT) result(hc)
      !! Host oracle: the SAME wet-area-weighted 4-cell corner thickness the
      !! energy kernel's Pass 2 builds (so abs_vort = q_corner·hc matches the
      !! kernel's corner_abs_vort bit-for-bit).
      integer, intent(in) :: ic, jc, k, nx, ny
      real(wp), intent(in) :: h(:, :, :), wet_T(:, :), areaT(:, :)
      real(wp) :: hc
      integer :: iw, ie, js, jn
      real(wp) :: aSW, aSE, aNW, aNE, num, den
      iw = max(1, ic - 1); ie = min(nx, ic); js = max(1, jc - 1); jn = min(ny, jc)
      aSW = wet_T(iw, js)*areaT(iw, js); aSE = wet_T(ie, js)*areaT(ie, js)
      aNW = wet_T(iw, jn)*areaT(iw, jn); aNE = wet_T(ie, jn)*areaT(ie, jn)
      num = aSW*h(iw, js, k) + aSE*h(ie, js, k) + aNW*h(iw, jn, k) + aNE*h(ie, jn, k)
      den = aSW + aSE + aNW + aNE
      hc = max(num/max(den, H_DIV_EPS), BC_H_MIN_PV)
   end function recon_h_corner

   subroutine test_bound_coriolis_allwet_inert(error)
      !! For an ALL-WET column the energy-scheme PV flux is already a CONVEX
      !! combination of the four neighbouring `(f+ζ)·v` velocity-form estimates
      !! (Roundabout's cell-mean corner thickness is exactly the mean the two
      !! v-face thicknesses telescope into, so the per-corner coefficients are
      !! `0.25·hf/h_corner` and sum to 1).  BOUND_CORIOLIS therefore CANNOT
      !! change the result there — this pins that: knob ON vs OFF must agree to
      !! round-off on a non-uniform-h + sheared-velocity all-wet state.  (It is
      !! the healthy-flow no-op guard AND documents why the clamp only bites at
      !! land-masked / floored corners.)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(coriolis_adv_t) :: cor
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: F0 = 1.0e-4_wp, DX = 1.0e3_wp
      real(wp), parameter :: PI = 3.141592653589793_wp
      integer, parameter :: NX_PHYS = 8, NY_PHYS = 8
      integer :: i, j, k, nx, ny
      real(wp), allocatable :: fx_off(:, :, :), fy_off(:, :, :)
      real(wp) :: max_dx, max_dy, scal

      call make_grid(grid, NX_PHYS, NY_PHYS, DX, DX)
      call make_cartesian_metrics(metrics, grid)
      ms%nz_ml = NZ
      call ms%init(grid)
      nx = grid%nx_total
      ny = grid%ny_total
      ! Non-uniform h + sheared u,v (varies in BOTH i and j ⇒ nonzero ζ, wide
      ! fv range); everything wet.
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx
               ms%h_layer(i, j, k) = 100.0_wp + 30.0_wp*sin(0.7_wp*real(i, wp)) &
                                     + 20.0_wp*cos(0.5_wp*real(j, wp)) + 2.0_wp*real(k, wp)
            end do
         end do
      end do
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx + 1
               ms%u_face_x_layer(i, j, k) = 0.3_wp*sin(0.4_wp*real(i, wp) + 0.2_wp*real(j, wp))
            end do
         end do
      end do
      do k = 1, NZ
         do j = 1, ny + 1
            do i = 1, nx
               ms%v_face_y_layer(i, j, k) = 0.4_wp*cos(0.3_wp*real(i, wp)) + 0.2_wp*real(j, wp)*0.05_wp
            end do
         end do
      end do

      allocate (fx_off(nx + 1, ny, NZ), fy_off(nx, ny + 1, NZ))

      ! Knob OFF
      cor%f_0 = F0
      cor%pv_variant = PV_VARIANT_SADOURNY_ENERGY
      cor%bound_coriolis = .false.
      call cor%init(grid, nz_ml=NZ)
      !$acc enter data copyin(ms, cor)
      call ms%enter_data(); call cor%enter_data()
      call coriolis_adv_compute_tendencies_sadourny_energy(grid, metrics, cor, ms, &
                                                           ms%u_face_x_layer, ms%v_face_y_layer, ms%h_layer)
      !$acc update self(cor%pv_flux_x%data, cor%pv_flux_y%data)
      fx_off = cor%pv_flux_x%data
      fy_off = cor%pv_flux_y%data
      call cor%exit_data(); call ms%exit_data()
      !$acc exit data delete(ms, cor)
      call cor%destroy()

      ! Knob ON (same state)
      cor%f_0 = F0
      cor%pv_variant = PV_VARIANT_SADOURNY_ENERGY
      cor%bound_coriolis = .true.
      call cor%init(grid, nz_ml=NZ)
      !$acc enter data copyin(ms, cor)
      call ms%enter_data(); call cor%enter_data()
      call coriolis_adv_compute_tendencies_sadourny_energy(grid, metrics, cor, ms, &
                                                           ms%u_face_x_layer, ms%v_face_y_layer, ms%h_layer)
      !$acc update self(cor%pv_flux_x%data, cor%pv_flux_y%data)
      max_dx = maxval(abs(cor%pv_flux_x%data - fx_off))
      max_dy = maxval(abs(cor%pv_flux_y%data - fy_off))
      call cor%exit_data(); call ms%exit_data()
      !$acc exit data delete(ms, cor)

      ! Inert to round-off (the convex-combination property; the ~1-ulp abs_vort
      ! recovery in the clamp is the only source of any difference).
      scal = max(maxval(abs(fx_off)), maxval(abs(fy_off)))
      call check(error, max_dx <= 1.0e-12_wp*scal, &
                 "BOUND_CORIOLIS perturbed the u-tendency on an all-wet column "// &
                 "(should be convex-bounded ⇒ clamp inert)")
      if (.not. allocated(error)) then
         call check(error, max_dy <= 1.0e-12_wp*scal, &
                    "BOUND_CORIOLIS perturbed the v-tendency on an all-wet column")
      end if

      deallocate (fx_off, fy_off)
      call cor%destroy(); call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_bound_coriolis_allwet_inert

   subroutine test_bound_coriolis_masked_clamps(error)
      !! A LAND cell (wet_T=0) breaks the cell↔face thickness telescoping, so
      !! `h_corner` (which excludes the land cell) no longer equals the mean of
      !! the v-face thicknesses (which include it) — the energy PV flux can now
      !! leave the `(f+ζ)·v` velocity-form range, and BOUND_CORIOLIS clamps it
      !! back.  With the land cell made anomalously THICK relative to its wet
      !! neighbours the face thicknesses over-amplify `q·vh` (the thin-edge PV
      !! blow-up mechanism).  Asserts (i) the clamp ENGAGES (bound-off exceeds
      !! the velocity-form range on ≥1 face) and (ii) the clamp GUARANTEE holds
      !! (bound-on lands inside `[min fv, max fv]` on every face).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(coriolis_adv_t) :: cor
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: F0 = 1.0e-4_wp, DX = 1.0e3_wp
      real(wp), parameter :: RTOL = 1.0e-9_wp
      integer, parameter :: NX_PHYS = 8, NY_PHYS = 8, KPROBE = 2
      integer :: i, j, k, nx, ny, iL, jL, n_exceed, n_violate
      real(wp), allocatable :: fx_off(:, :, :), ke_off(:, :, :), q_off(:, :, :)
      real(wp), allocatable :: fx_on(:, :, :)
      real(wp) :: avN, avS, fv1, fv2, fv3, fv4, lo, hi, keg, ppart_off, ppart_on, sc

      call make_grid(grid, NX_PHYS, NY_PHYS, DX, DX)
      call make_cartesian_metrics(metrics, grid)
      ms%nz_ml = NZ
      nx = grid%nx_total
      ny = grid%ny_total
      iL = NGHOST + 4
      jL = NGHOST + 4
      ! Land cell: exclude it from the wet-area corner-h averages (Pass 2), but
      ! it still enters the v-face thickness (Pass 3b) ⇒ telescoping broken.
      metrics%wet_T(iL, jL) = 0.0_wp
      !$acc update device(metrics%wet_T)

      call ms%init(grid)
      ! Anomalously THICK land cell vs thin wet neighbours ⇒ over-amplified vh.
      ms%h_layer = 1.0_wp
      do k = 1, NZ
         ms%h_layer(iL, jL, k) = 200.0_wp
      end do
      ms%u_face_x_layer = 0.0_wp
      ! Sheared v (varies in i) ⇒ the four v-points around a u-face differ.
      do k = 1, NZ
         do j = 1, ny + 1
            do i = 1, nx
               ms%v_face_y_layer(i, j, k) = 0.5_wp + 0.3_wp*real(i, wp)
            end do
         end do
      end do

      allocate (fx_off(nx + 1, ny, NZ), ke_off(nx, ny, NZ), q_off(nx + 1, ny + 1, NZ))
      allocate (fx_on(nx + 1, ny, NZ))

      cor%f_0 = F0
      cor%pv_variant = PV_VARIANT_SADOURNY_ENERGY
      cor%bound_coriolis = .false.
      call cor%init(grid, nz_ml=NZ)
      !$acc enter data copyin(ms, cor)
      call ms%enter_data(); call cor%enter_data()
      call coriolis_adv_compute_tendencies_sadourny_energy(grid, metrics, cor, ms, &
                                                           ms%u_face_x_layer, ms%v_face_y_layer, ms%h_layer)
      !$acc update self(cor%pv_flux_x%data, cor%ke_centre%data, cor%q_corner%data)
      fx_off = cor%pv_flux_x%data
      ke_off = cor%ke_centre%data
      q_off = cor%q_corner%data
      call cor%exit_data(); call ms%exit_data()
      !$acc exit data delete(ms, cor)
      call cor%destroy()

      cor%f_0 = F0
      cor%pv_variant = PV_VARIANT_SADOURNY_ENERGY
      cor%bound_coriolis = .true.
      call cor%init(grid, nz_ml=NZ)
      !$acc enter data copyin(ms, cor)
      call ms%enter_data(); call cor%enter_data()
      call coriolis_adv_compute_tendencies_sadourny_energy(grid, metrics, cor, ms, &
                                                           ms%u_face_x_layer, ms%v_face_y_layer, ms%h_layer)
      !$acc update self(cor%pv_flux_x%data)
      fx_on = cor%pv_flux_x%data
      call cor%exit_data(); call ms%exit_data()
      !$acc exit data delete(ms, cor)

      ! Reconstruct pv_part = pv_flux + ke_grad (the pre-KE PV flux the clamp
      ! acts on) and the velocity-form range at every interior u-face.
      n_exceed = 0
      n_violate = 0
      k = KPROBE
      do j = 2, ny - 1
         do i = 3, nx - 1
            keg = (ke_off(i, j, k) - ke_off(i - 1, j, k))*metrics%idxCu(i, j)
            ppart_off = fx_off(i, j, k) + keg
            ppart_on = fx_on(i, j, k) + keg
            avN = q_off(i, j + 1, k)*recon_h_corner(i, j + 1, k, nx, ny, &
                                                    ms%h_layer, metrics%wet_T, metrics%areaT)
            avS = q_off(i, j, k)*recon_h_corner(i, j, k, nx, ny, &
                                                ms%h_layer, metrics%wet_T, metrics%areaT)
            fv1 = avN*ms%v_face_y_layer(i - 1, j + 1, k)
            fv2 = avN*ms%v_face_y_layer(i, j + 1, k)
            fv3 = avS*ms%v_face_y_layer(i - 1, j, k)
            fv4 = avS*ms%v_face_y_layer(i, j, k)
            hi = max(max(fv1, fv2), max(fv3, fv4))
            lo = min(min(fv1, fv2), min(fv3, fv4))
            sc = max(abs(hi), abs(lo)) + 1.0e-30_wp
            if (ppart_off > hi*(1.0_wp + RTOL) + RTOL*sc .or. &
                ppart_off < lo*(1.0_wp) - RTOL*sc) then
               ! bound-off left the velocity-form range on this face.
               if (ppart_off > hi + RTOL*sc .or. ppart_off < lo - RTOL*sc) n_exceed = n_exceed + 1
            end if
            ! bound-ON must be inside [lo, hi] (the clamp guarantee).
            if (ppart_on > hi + RTOL*sc .or. ppart_on < lo - RTOL*sc) n_violate = n_violate + 1
         end do
      end do

      call check(error, n_exceed > 0, &
                 "BOUND_CORIOLIS clamp never engaged (masked telescoping-break "// &
                 "did not push the PV flux outside the velocity-form range)")
      if (.not. allocated(error)) then
         call check(error, n_violate == 0, &
                    "BOUND_CORIOLIS bound-on PV flux left [min fv, max fv] "// &
                    "(clamp guarantee violated)")
      end if

      deallocate (fx_off, ke_off, q_off, fx_on)
      call cor%destroy(); call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_bound_coriolis_masked_clamps

   ! -----------------------------------------------------------------
   ! corner_h (MOM6 area-weighted PV corner thickness)
   ! -----------------------------------------------------------------

   subroutine run_energy_qcorner(corner_variant, grid, ms, metrics, q_out)
      !! Run the energy CorAdv with the given corner_h variant and return the
      !! q_corner field (Pass 2 writes q into q_corner; Pass 5/6 only read it).
      integer, intent(in) :: corner_variant
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_metrics_t), intent(in) :: metrics
      real(wp), allocatable, intent(out) :: q_out(:, :, :)
      type(coriolis_adv_t) :: cor
      cor%f_0 = 1.0e-4_wp
      cor%pv_variant = PV_VARIANT_SADOURNY_ENERGY
      cor%corner_h_variant = corner_variant
      call cor%init(grid, nz_ml=NZ)
      !$acc enter data copyin(ms, cor)
      call ms%enter_data(); call cor%enter_data()
      call coriolis_adv_compute_tendencies_sadourny_energy(grid, metrics, cor, ms, &
                                                           ms%u_face_x_layer, ms%v_face_y_layer, ms%h_layer)
      !$acc update self(cor%q_corner%data)
      allocate (q_out, source=cor%q_corner%data)
      call cor%exit_data(); call ms%exit_data()
      !$acc exit data delete(ms, cor)
      call cor%destroy()
   end subroutine run_energy_qcorner

   subroutine test_corner_h_equiv(error)
      !! ABOVE the floor, MOM6's `q = abs_vort·Area_q/(hArea_q + vol_neglect)`
      !! and Roundabout's `abs_vort/(hArea_q/Area_q)` are the SAME construction
      !! (identical numerator + denominator) — they agree to ROUND-OFF (only a
      !! floating-point operation-order difference + a negligible vol_neglect).
      !! This pins the algebraic equivalence that says corner_h is a no-op on
      !! any realistic (non-vanishing) column.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      real(wp), allocatable :: q_cell(:, :, :), q_mom6(:, :, :)
      real(wp), parameter :: DX = 1.0e3_wp
      integer, parameter :: NX_PHYS = 8, NY_PHYS = 8
      integer :: i, j, k, nx, ny
      real(wp) :: rel

      call make_grid(grid, NX_PHYS, NY_PHYS, DX, DX)
      call make_cartesian_metrics(metrics, grid)
      ms%nz_ml = NZ
      call ms%init(grid)
      nx = grid%nx_total
      ny = grid%ny_total
      ! Non-uniform, strictly-positive h + sheared velocity (nonzero ζ), all wet.
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx
               ms%h_layer(i, j, k) = 100.0_wp + 25.0_wp*sin(0.6_wp*real(i, wp)) &
                                     + 15.0_wp*cos(0.4_wp*real(j, wp))
            end do
         end do
      end do
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx + 1
               ms%u_face_x_layer(i, j, k) = 0.2_wp*sin(0.3_wp*real(j, wp))
            end do
         end do
      end do
      ms%v_face_y_layer = 0.0_wp

      call run_energy_qcorner(CORNER_H_CELL_MEAN, grid, ms, metrics, q_cell)
      call run_energy_qcorner(CORNER_H_MOM6_AREA, grid, ms, metrics, q_mom6)

      ! Max relative q difference over the interior corners.
      rel = 0.0_wp
      do k = 1, NZ
         do j = 2, ny
            do i = 2, nx
               if (abs(q_cell(i, j, k)) > 0.0_wp) then
                  rel = max(rel, abs(q_mom6(i, j, k) - q_cell(i, j, k))/abs(q_cell(i, j, k)))
               end if
            end do
         end do
      end do
      call check(error, rel <= 1.0e-12_wp, &
                 "corner_h: mom6_area and cell_mean must agree to round-off above the "// &
                 "floor (the constructions are algebraically identical)")

      deallocate (q_cell, q_mom6)
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_corner_h_equiv

   subroutine test_corner_h_floor(error)
      !! At a corner whose four cells are all far BELOW the H_MIN_PV thickness
      !! floor, the two constructions diverge exactly per the guard: cell_mean
      !! CAPS `q` at `abs_vort/H_MIN_PV`, while mom6_area's vol_neglect is pure
      !! 1/0 armor so `q ≈ abs_vort/h_thin` — much LARGER.  Hand oracle: with a
      !! block at `h_thin = 1e-15` (uniform area, ζ≈0 ⇒ abs_vort = f), the ratio
      !! `q_mom6/q_cell ≈ H_MIN_PV/h_thin = 1e3`.  (This is the ONLY regime the
      !! knob acts in — realistic layers never reach it, so corner_h cannot move
      !! a thin-LAYER (h≳1e-12) instability; if anything MOM6's floor is weaker.)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_metrics_t) :: metrics
      real(wp), allocatable :: q_cell(:, :, :), q_mom6(:, :, :)
      real(wp), parameter :: DX = 1.0e3_wp, H_THIN = 1.0e-15_wp
      real(wp), parameter :: H_MIN_PV = 1.0e-12_wp
      integer, parameter :: NX_PHYS = 8, NY_PHYS = 8, KPROBE = 1
      integer :: i, j, nx, ny, ic, jc
      real(wp) :: ratio, ratio_expect

      call make_grid(grid, NX_PHYS, NY_PHYS, DX, DX)
      call make_cartesian_metrics(metrics, grid)
      ms%nz_ml = NZ
      call ms%init(grid)
      nx = grid%nx_total
      ny = grid%ny_total
      ! Background 1 m, a 3×3 ultra-thin block centred on cell (ic0,jc0) so the
      ! corner (ic,jc) at the block interior has all four cells = H_THIN.
      ms%h_layer = 1.0_wp
      do j = NGHOST + 3, NGHOST + 5
         do i = NGHOST + 3, NGHOST + 5
            ms%h_layer(i, j, :) = H_THIN
         end do
      end do
      ms%u_face_x_layer = 0.0_wp   ! ζ = 0 ⇒ abs_vort = f_corner = F0
      ms%v_face_y_layer = 0.0_wp
      ! Corner (ic,jc) whose 4 cells (ic-1,jc-1),(ic,jc-1),(ic-1,jc),(ic,jc)
      ! are all inside the thin block.
      ic = NGHOST + 5
      jc = NGHOST + 5

      call run_energy_qcorner(CORNER_H_CELL_MEAN, grid, ms, metrics, q_cell)
      call run_energy_qcorner(CORNER_H_MOM6_AREA, grid, ms, metrics, q_mom6)

      ! cell_mean caps at abs_vort/H_MIN_PV; mom6_area ≈ abs_vort/H_THIN.
      ratio = q_mom6(ic, jc, KPROBE)/q_cell(ic, jc, KPROBE)
      ratio_expect = H_MIN_PV/H_THIN   ! = 1000
      call check(error, abs(ratio - ratio_expect) <= 1.0e-3_wp*ratio_expect, &
                 "corner_h floor: mom6_area q must exceed cell_mean q by ~H_MIN_PV/h_thin "// &
                 "(cell_mean caps, mom6_area does not)")

      deallocate (q_cell, q_mom6)
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_corner_h_floor

   ! -----------------------------------------------------------------
   ! WENO3 PV reconstruction (F1)
   ! -----------------------------------------------------------------

   subroutine test_weno3_constant(error)
      !! A locally constant vorticity field must reconstruct to the
      !! constant exactly (both upwind signs) -- the WENO-Z degenerate
      !! guard must not perturb it.
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: C = 3.14_wp
      call check(error, abs(weno3_recon(C, C, C, C, 1.0_wp) - C) < 1.0e-14_wp, &
                 "weno3 constant (u>0) not exact")
      if (allocated(error)) return
      call check(error, abs(weno3_recon(C, C, C, C, -1.0_wp) - C) < 1.0e-14_wp, &
                 "weno3 constant (u<0) not exact")
   end subroutine test_weno3_constant

   subroutine test_weno3_linear(error)
      !! On a linear ramp both candidate stencils are exact, so the
      !! reconstruction equals the linear face value for ANY weights --
      !! q = 2x+1 at x=-1,0,1,2 -> face at x=0.5 -> value 2.0.  Both signs.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: qf
      qf = weno3_recon(-1.0_wp, 1.0_wp, 3.0_wp, 5.0_wp, 1.0_wp)
      call check(error, abs(qf - 2.0_wp) < 1.0e-13_wp, "weno3 linear (u>0) not exact")
      if (allocated(error)) return
      qf = weno3_recon(-1.0_wp, 1.0_wp, 3.0_wp, 5.0_wp, -1.0_wp)
      call check(error, abs(qf - 2.0_wp) < 1.0e-13_wp, "weno3 linear (u<0) not exact")
   end subroutine test_weno3_linear

   subroutine test_weno3_step(error)
      !! A step [0,0,1,1] advected from the smooth (left) side must give
      !! an essentially-non-oscillatory value inside [0,1] with NO Gibbs
      !! overshoot, collapsing onto the smooth (left) stencil (~0).
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: qf
      qf = weno3_recon(0.0_wp, 0.0_wp, 1.0_wp, 1.0_wp, 1.0_wp)
      call check(error, qf >= -1.0e-12_wp .and. qf <= 1.0_wp + 1.0e-12_wp, &
                 "weno3 step overshoots [0,1] (not ENO)")
      if (allocated(error)) return
      call check(error, qf < 1.0e-6_wp, "weno3 step did not pick the smooth upwind stencil")
   end subroutine test_weno3_step

   subroutine test_weno3_upwind(error)
      !! Upwind bias: on a curved (quadratic) field q=j^2 sampled at
      !! j=-1,0,1,2 -> (1,0,1,4), the left stencil {1,0,1} is smooth
      !! (symmetric) while the right stencil carries the qp2=4 curvature, so
      !! the u>0 and u<0 reconstructions must differ materially.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: qp, qm
      qp = weno3_recon(1.0_wp, 0.0_wp, 1.0_wp, 4.0_wp, 1.0_wp)
      qm = weno3_recon(1.0_wp, 0.0_wp, 1.0_wp, 4.0_wp, -1.0_wp)
      call check(error, abs(qp - qm) > 1.0e-3_wp, &
                 "weno3 not upwind-biased (result independent of advecting sign)")
   end subroutine test_weno3_upwind

   subroutine test_weno3_parse(error)
      !! parse_pv_adv_scheme + the implemented gate + the per-rung nghost
      !! requirement: centered/weno3/weno5/weno7 all ship; a typo is INVALID;
      !! weno5 needs nghost>=3, weno7 needs nghost>=4.
      type(error_type), allocatable, intent(out) :: error
      call check(error, parse_pv_adv_scheme("centered") == PV_ADV_CENTERED .and. &
                 parse_pv_adv_scheme("weno3") == PV_ADV_WENO3 .and. &
                 parse_pv_adv_scheme("weno5") == PV_ADV_WENO5 .and. &
                 parse_pv_adv_scheme("weno7") == PV_ADV_WENO7 .and. &
                 parse_pv_adv_scheme("bogus") == PV_ADV_INVALID, "parse_pv_adv_scheme wrong")
      if (allocated(error)) return
      call check(error, pv_adv_scheme_is_implemented(PV_ADV_CENTERED) .and. &
                 pv_adv_scheme_is_implemented(PV_ADV_WENO3) .and. &
                 pv_adv_scheme_is_implemented(PV_ADV_WENO5) .and. &
                 pv_adv_scheme_is_implemented(PV_ADV_WENO7) .and. &
                 (.not. pv_adv_scheme_is_implemented(PV_ADV_INVALID)), &
                 "pv_adv_scheme_is_implemented gate wrong")
      if (allocated(error)) return
      call check(error, pv_adv_required_nghost(PV_ADV_CENTERED) == 2 .and. &
                 pv_adv_required_nghost(PV_ADV_WENO3) == 2 .and. &
                 pv_adv_required_nghost(PV_ADV_WENO5) == 3 .and. &
                 pv_adv_required_nghost(PV_ADV_WENO7) == 4, &
                 "pv_adv_required_nghost wrong")
   end subroutine test_weno3_parse

   subroutine run_coriolis_u(grid, metrics, scheme, uprof, uout)
      !! Set u = uprof(j) with a UNIFORM v (f-plane), run the Sadourny
      !! Coriolis-adv with the given pv_adv_scheme, forward-Euler apply,
      !! return the resulting u field.  The uniform v gives a non-zero
      !! v_at_u (so pv_flux_x = absvort*v_at_u is live) with no d/dx
      !! contribution to the vorticity, so the absolute vorticity varies
      !! only in j -- exactly the direction the u-face WENO stencil spans.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      integer, intent(in) :: scheme
      real(wp), intent(in) :: uprof(:)
      real(wp), allocatable, intent(out) :: uout(:, :, :)
      type(multilayer_state_t) :: ms
      type(coriolis_adv_t) :: cor
      real(wp), parameter :: DT = 0.01_wp
      integer :: j, ny
      ny = grid%ny_total
      ms%nz_ml = NZ; call ms%init(grid)
      cor%f_0 = 1.0e-4_wp
      call cor%init(grid, nz_ml=NZ)
      cor%pv_adv_scheme = scheme
      ms%h_layer = 1.0_wp
      ms%v_face_y_layer = 1.0e-3_wp
      do j = 1, ny
         ms%u_face_x_layer(:, j, :) = uprof(j)
      end do
      call map_in(ms, cor)
      call coriolis_adv_compute_tendencies(grid, metrics, cor, ms)
      call coriolis_adv_apply_tendencies(cor, ms, DT)
      call map_out(ms, cor)
      uout = ms%u_face_x_layer
      call cor%destroy(); call ms%destroy()
   end subroutine run_coriolis_u

   subroutine test_weno_uniform(error)
      !! A uniform-shear u = U0*j gives spatially CONSTANT absolute
      !! vorticity, which weno3 must reconstruct exactly -> the u-tendency
      !! is bit-identical to the centred scheme over the interior (away
      !! from the wall-zeroed corner band the WENO stencil never reaches).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      real(wp), allocatable :: uprof(:), u_c(:, :, :), u_w(:, :, :)
      real(wp) :: dmax
      integer :: j, ny, lo, hi
      call make_grid(grid, 12, 12, 1.0_wp, 1.0_wp)
      call make_cartesian_metrics(metrics, grid)
      ny = grid%ny_total
      allocate (uprof(ny))
      do j = 1, ny
         uprof(j) = 1.0e-3_wp*real(j, wp)
      end do
      call run_coriolis_u(grid, metrics, PV_ADV_CENTERED, uprof, u_c)
      call run_coriolis_u(grid, metrics, PV_ADV_WENO3, uprof, u_w)
      ! Interior band: >=3 cells from the j-array edges (WENO stencil radius 2
      ! + the centred fallback band) -> uniform vorticity, must match exactly.
      lo = 4; hi = ny - 3
      dmax = maxval(abs(u_c(:, lo:hi, :) - u_w(:, lo:hi, :)))
      call check(error, dmax < 1.0e-14_wp, &
                 "weno3 differs from centred on a uniform-vorticity interior")
      call destroy_cartesian_metrics(metrics)
   end subroutine test_weno_uniform

   subroutine test_weno_sharp(error)
      !! A sharp shear line (u steps 0 -> U0 at j=jc) makes the absolute
      !! vorticity strongly non-uniform in j, so the u-face weno3
      !! reconstruction MUST diverge from the centred 2-point average near
      !! the front -- proving the reconstruction is live (not silently the
      !! centred path).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      real(wp), allocatable :: uprof(:), u_c(:, :, :), u_w(:, :, :)
      real(wp) :: dmax
      integer :: j, ny, jc
      call make_grid(grid, 12, 12, 1.0_wp, 1.0_wp)
      call make_cartesian_metrics(metrics, grid)
      ny = grid%ny_total
      jc = ny/2
      allocate (uprof(ny))
      do j = 1, ny
         uprof(j) = merge(1.0e-2_wp, 0.0_wp, j >= jc)
      end do
      call run_coriolis_u(grid, metrics, PV_ADV_CENTERED, uprof, u_c)
      call run_coriolis_u(grid, metrics, PV_ADV_WENO3, uprof, u_w)
      dmax = maxval(abs(u_c - u_w))
      call check(error, dmax > 1.0e-12_wp, &
                 "weno3 identical to centred on a sharp front (reconstruction inert)")
      call destroy_cartesian_metrics(metrics)
   end subroutine test_weno_sharp

   subroutine test_weno5_recon(error)
      !! weno5_recon: exact on constant + linear (q=2x+1, x=-2..3 -> face@0.5
      !! = 2.0, both signs) and ENO (no overshoot) on a step.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: qf
      qf = weno5_recon(-3.0_wp, -1.0_wp, 1.0_wp, 3.0_wp, 5.0_wp, 7.0_wp, 1.0_wp)
      call check(error, abs(qf - 2.0_wp) < 1.0e-12_wp, "weno5 linear (u>0) not exact")
      if (allocated(error)) return
      qf = weno5_recon(-3.0_wp, -1.0_wp, 1.0_wp, 3.0_wp, 5.0_wp, 7.0_wp, -1.0_wp)
      call check(error, abs(qf - 2.0_wp) < 1.0e-12_wp, "weno5 linear (u<0) not exact")
      if (allocated(error)) return
      call check(error, abs(weno5_recon(5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 1.0_wp) &
                            - 5.0_wp) < 1.0e-14_wp, "weno5 constant not exact")
      if (allocated(error)) return
      qf = weno5_recon(0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp)
      call check(error, qf >= -1.0e-12_wp .and. qf <= 1.0_wp + 1.0e-12_wp, &
                 "weno5 step overshoots [0,1] (not ENO)")
   end subroutine test_weno5_recon

   subroutine test_weno7_recon(error)
      !! weno7_recon: exact on constant + linear (q=2x+1, x=-3..4 -> face@0.5
      !! = 2.0, both signs) and ENO (no overshoot) on a step.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: qf
      qf = weno7_recon(-5.0_wp, -3.0_wp, -1.0_wp, 1.0_wp, 3.0_wp, 5.0_wp, 7.0_wp, 9.0_wp, 1.0_wp)
      call check(error, abs(qf - 2.0_wp) < 1.0e-11_wp, "weno7 linear (u>0) not exact")
      if (allocated(error)) return
      qf = weno7_recon(-5.0_wp, -3.0_wp, -1.0_wp, 1.0_wp, 3.0_wp, 5.0_wp, 7.0_wp, 9.0_wp, -1.0_wp)
      call check(error, abs(qf - 2.0_wp) < 1.0e-11_wp, "weno7 linear (u<0) not exact")
      if (allocated(error)) return
      call check(error, abs(weno7_recon(5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, &
                                        5.0_wp, 1.0_wp) - 5.0_wp) < 1.0e-14_wp, "weno7 constant not exact")
      if (allocated(error)) return
      qf = weno7_recon(0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp)
      call check(error, qf >= -1.0e-12_wp .and. qf <= 1.0_wp + 1.0e-12_wp, &
                 "weno7 step overshoots [0,1] (not ENO)")
   end subroutine test_weno7_recon

   subroutine test_weno5_kernel(error)
      !! Kernel-level: weno5 (radius-3 stencil) reduces to centred on a
      !! uniform-vorticity column over the interior where the wider stencil
      !! never reaches the wall-zeroed corners -- proves the weno5 branch of
      !! the Coriolis kernel is wired + degrades to centred correctly.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      real(wp), allocatable :: uprof(:), u_c(:, :, :), u_w(:, :, :)
      integer :: j, ny, lo, hi
      call make_grid(grid, 16, 16, 1.0_wp, 1.0_wp)
      call make_cartesian_metrics(metrics, grid)
      ny = grid%ny_total
      allocate (uprof(ny))
      do j = 1, ny
         uprof(j) = 1.0e-3_wp*real(j, wp)
      end do
      call run_coriolis_u(grid, metrics, PV_ADV_CENTERED, uprof, u_c)
      call run_coriolis_u(grid, metrics, PV_ADV_WENO5, uprof, u_w)
      lo = 6; hi = ny - 5   ! inside the radius-3 stencil's reach of the edges
      call check(error, maxval(abs(u_c(:, lo:hi, :) - u_w(:, lo:hi, :))) < 1.0e-14_wp, &
                 "weno5 differs from centred on a uniform-vorticity interior")
      call destroy_cartesian_metrics(metrics)
   end subroutine test_weno5_kernel

end module test_coriolis_multilayer
