!! Unit tests for the ocean horizontal-viscosity kernel
!! (rdb_ocean_horizontal_viscosity).  Laplacian momentum diffusion at
!! the C-grid face velocities with closed-wall + free-slip boundary
!! conditions, applied per layer.
!!
!! Cases:
!!   * Uniform velocity — Laplacian vanishes exactly; the kernel must
!!     leave u, v bit-for-bit at the IC.  Constancy preservation.
!!   * Zero-coefficient short-circuit — `nu_h = 0` is a no-op even
!!     for non-uniform velocity.  Lets the driver wire the kernel
!!     unconditionally and gate via the namelist.
!!   * Sinusoid decay — initialise a single Fourier mode in u
!!     (cos(k*x)) and run several steps.  Analytic amplitude after
!!     N steps is (1 - dt * nu * k^2)^N.  Check the realised decay
!!     matches to better than 1% (centred 5-point Laplacian has a
!!     known wavenumber-dependent error, exact at low k).
!!   * KE monotone decay — non-trivial u, v field.  Total KE
!!     (sum 0.5 * (u^2 + v^2)) must decrease step-over-step.  No
!!     spurious energy injection regardless of state.
module test_ocean_hvisc
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t, metrics_fill_cartesian, &
                                metrics_finalize, metrics_apply_land_mask
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_horizontal_viscosity, only: &
      ocean_horizontal_viscosity_t, &
      ocean_horizontal_viscosity_compute_tendencies, &
      ocean_horizontal_viscosity_apply_tendencies, &
      ocean_horizontal_viscosity_compute_ke_diss
   implicit none
   private

   public :: collect_ocean_hvisc_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 3

contains

   subroutine collect_ocean_hvisc_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("uniform_velocity_no_op", test_uniform_no_op), &
                  new_unittest("zero_coefficient_no_op", test_zero_coeff_no_op), &
                  new_unittest("sinusoid_decay_rate", test_sinusoid_decay), &
                  new_unittest("kinetic_energy_monotone", test_ke_monotone), &
                  new_unittest("T5_metric_reduces_to_scalar_laplacian", &
                               test_metric_reduction), &
                  new_unittest("stress_tensor_reduces_to_laplacian", &
                               test_stress_reduces_to_laplacian), &
                  new_unittest("stress_tensor_momentum_conserving", &
                               test_stress_momentum_conservation), &
                  new_unittest("stress_tensor_cfl_clamp_finite", &
                               test_stress_cfl_clamp), &
                  new_unittest("stress_tensor_coast_mask", &
                               test_stress_coast_mask), &
                  new_unittest("ke_diss_rate_frictional_source", test_ke_diss) &
                  ]
   end subroutine collect_ocean_hvisc_tests

   ! ------------------------------------------------------------------
   ! `compute_ke_diss` fills the lateral-viscosity kinetic-energy
   ! dissipation rate that sources MEKE's frictional term (MOM6-inspired).
   ! Drives the REAL path (compute_tendencies fills du_visc), then checks
   ! the kernel against an INDEPENDENT host re-derivation:
   !   ke_diss(i,j) = Σ_k ρ_k h_k · [½(u_i·du_i + u_{i+1}·du_{i+1})
   !                                  + ½(v_j·dv_j + v_{j+1}·dv_{j+1})].
   ! For a sinusoidal field the viscous tendency opposes the velocity
   ! (du_visc = ν·∇²u = −ν·k²·u), so ke_diss ≤ 0 pointwise (dissipation).
   ! ------------------------------------------------------------------
   subroutine test_ke_diss(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_horizontal_viscosity_t) :: hv
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp), parameter :: NU = 0.5_wp, RHO = 1025.0_wp, HZ = 10.0_wp
      integer :: i, j, k, nx, ny
      real(wp) :: ref, ke, max_rel, max_pos
      checks: block
         call make_grid(grid, 16, 14, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call hv%init(grid, nz_ml=NZ)
         hv%nu_h = NU
         hv%compute_ke_diss = .true.
         nx = grid%nx_total
         ny = grid%ny_total
         ms%h_layer = HZ
         ms%rho_layer = RHO
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx + 1
                  ms%u_face_x_layer(i, j, k) = sin(2.0_wp*PI*real(i, wp)/real(nx, wp))
               end do
            end do
            do j = 1, ny + 1
               do i = 1, nx
                  ms%v_face_y_layer(i, j, k) = sin(2.0_wp*PI*real(j, wp)/real(ny, wp))
               end do
            end do
         end do
         call map_in(ms, hv, metrics, grid)
         ! Real path: fill du_visc/dv_visc, then the KE-dissipation rate.
         call ocean_horizontal_viscosity_compute_tendencies(grid, metrics, hv, ms)
         call ocean_horizontal_viscosity_compute_ke_diss(hv, ms)
         !$acc update self(hv%ke_diss, hv%du_visc%data, hv%dv_visc%data)
         !$acc update self(ms%u_face_x_layer, ms%v_face_y_layer)
         call map_out(ms, hv, metrics)
         ! Independent host re-derivation at every interior cell.
         max_rel = 0.0_wp
         max_pos = 0.0_wp
         do j = NGHOST + 1, ny - NGHOST
            do i = NGHOST + 1, nx - NGHOST
               ref = 0.0_wp
               do k = 1, NZ
                  ke = 0.5_wp*(ms%u_face_x_layer(i, j, k)*hv%du_visc%data(i, j, k) &
                               + ms%u_face_x_layer(i + 1, j, k)*hv%du_visc%data(i + 1, j, k)) &
                       + 0.5_wp*(ms%v_face_y_layer(i, j, k)*hv%dv_visc%data(i, j, k) &
                                 + ms%v_face_y_layer(i, j + 1, k)*hv%dv_visc%data(i, j + 1, k))
                  ref = ref + RHO*HZ*ke
               end do
               max_rel = max(max_rel, abs(hv%ke_diss(i, j) - ref))
               max_pos = max(max_pos, hv%ke_diss(i, j))
            end do
         end do
         call check(error, max_rel < 1.0e-8_wp, &
                    "ke_diss must equal the host Σ ρ·h·(u·du_visc) re-derivation")
         if (allocated(error)) exit checks
         call check(error, max_pos <= 1.0e-12_wp, &
                    "viscous KE dissipation ⇒ ke_diss ≤ 0 everywhere")
      end block checks
      call hv%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_ke_diss

   subroutine make_grid(grid, nx_phys, ny_phys, dx, dy)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      real(wp), intent(in) :: dx, dy
      call grid%init(nx_phys, ny_phys, NGHOST, dx, dy)
   end subroutine make_grid

   subroutine map_in(ms, hv, metrics, grid)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_horizontal_viscosity_t), intent(inout) :: hv
      type(ocean_metrics_t), intent(inout) :: metrics
      type(hgrid_t), intent(in) :: grid
      call make_cartesian_metrics(metrics, grid)
      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(hv)
      call hv%enter_data()
   end subroutine map_in

   subroutine map_out(ms, hv, metrics)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_horizontal_viscosity_t), intent(inout) :: hv
      type(ocean_metrics_t), intent(inout) :: metrics
      call hv%exit_data()
      !$acc exit data delete(hv)
      call ms%exit_data()
      !$acc exit data delete(ms)
      call destroy_cartesian_metrics(metrics)
   end subroutine map_out

   ! -----------------------------------------------------------------
   ! Cases
   ! -----------------------------------------------------------------

   subroutine test_uniform_no_op(error)
      !! Uniform u = U0, v = V0 across the domain.  Laplacian is
      !! exactly zero at every interior face, and wall faces are
      !! force-zeroed by the kernel.  After one apply step u and v
      !! must match IC bit-for-bit.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_horizontal_viscosity_t) :: hv
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: U0 = 0.5_wp, V0 = -0.3_wp
      real(wp), parameter :: DT = 0.1_wp
      real(wp) :: max_du, max_dv
      checks: block

         call make_grid(grid, 12, 10, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call hv%init(grid, nz_ml=NZ)
         hv%nu_h = 1.0_wp

         ms%h_layer = 10.0_wp
         ms%u_face_x_layer = U0
         ms%v_face_y_layer = V0

         call map_in(ms, hv, metrics, grid)
         call ocean_horizontal_viscosity_compute_tendencies(grid, metrics, hv, ms)
         call ocean_horizontal_viscosity_apply_tendencies(hv, ms, DT)
         call map_out(ms, hv, metrics)

         max_du = maxval(abs(ms%u_face_x_layer - U0))
         max_dv = maxval(abs(ms%v_face_y_layer - V0))

         call check(error, max_du < 1.0e-12_wp, &
                    "uniform u: hvisc tendency leaked")
         if (allocated(error)) exit checks
         call check(error, max_dv < 1.0e-12_wp, &
                    "uniform v: hvisc tendency leaked")

      end block checks
      call hv%destroy(); call ms%destroy()
   end subroutine test_uniform_no_op

   subroutine test_zero_coeff_no_op(error)
      !! `nu_h = 0` short-circuits the compute kernel.  Apply
      !! reads the zeroed tendency and leaves a non-trivial u, v
      !! pattern untouched.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_horizontal_viscosity_t) :: hv
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp), parameter :: DT = 0.1_wp
      real(wp), allocatable :: u_ic(:, :, :), v_ic(:, :, :)
      real(wp) :: max_du, max_dv
      integer :: i, j, k, nx, ny
      checks: block

         call make_grid(grid, 12, 10, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call hv%init(grid, nz_ml=NZ)
         hv%nu_h = 0.0_wp
         nx = grid%nx_total
         ny = grid%ny_total

         ms%h_layer = 10.0_wp
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx + 1
                  ms%u_face_x_layer(i, j, k) = sin(2.0_wp*PI*real(i, wp)/real(nx, wp))
               end do
            end do
            do j = 1, ny + 1
               do i = 1, nx
                  ms%v_face_y_layer(i, j, k) = cos(2.0_wp*PI*real(j, wp)/real(ny, wp))
               end do
            end do
         end do
         allocate (u_ic, source=ms%u_face_x_layer)
         allocate (v_ic, source=ms%v_face_y_layer)

         call map_in(ms, hv, metrics, grid)
         call ocean_horizontal_viscosity_compute_tendencies(grid, metrics, hv, ms)
         call ocean_horizontal_viscosity_apply_tendencies(hv, ms, DT)
         call map_out(ms, hv, metrics)

         max_du = maxval(abs(ms%u_face_x_layer - u_ic))
         max_dv = maxval(abs(ms%v_face_y_layer - v_ic))

         call check(error, max_du < 1.0e-12_wp, &
                    "nu_h=0: u must be unchanged")
         if (allocated(error)) exit checks
         call check(error, max_dv < 1.0e-12_wp, &
                    "nu_h=0: v must be unchanged")

      end block checks
      deallocate (u_ic, v_ic)
      call hv%destroy(); call ms%destroy()
   end subroutine test_zero_coeff_no_op

   subroutine test_sinusoid_decay(error)
      !! u(i, j) = cos(k * x_i) with k = 2*pi / Lx, v = 0.
      !! For the centred 5-point Laplacian
      !!   nabla^2 (cos(k*x)) = -beta^2 cos(k*x)   where
      !!   beta^2 = 2*(1 - cos(k*dx)) / dx^2
      !! so each forward-Euler step multiplies u by (1 - dt * nu * beta^2).
      !! After N steps amplitude is (1 - dt * nu * beta^2)^N.  We compare
      !! the realised peak (away from y-boundary rows where the kernel
      !! zeros the tendency) to the analytic prediction.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_horizontal_viscosity_t) :: hv
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp), parameter :: DT = 0.05_wp
      real(wp), parameter :: NU = 0.05_wp
      integer, parameter :: N_STEPS = 10
      integer :: i, j, k, step, nx, ny, j_probe, i_probe
      real(wp) :: k_wave, beta_sq, alpha
      real(wp) :: amp_expected, amp_obs, rel_err, ic_at_probe

      call make_grid(grid, 40, 12, 1.0_wp, 1.0_wp)
      ms%nz_ml = NZ
      call ms%init(grid)
      call hv%init(grid, nz_ml=NZ)
      hv%nu_h = NU
      nx = grid%nx_total
      ny = grid%ny_total

      k_wave = 2.0_wp*PI/real(nx, wp)
      beta_sq = 2.0_wp*(1.0_wp - cos(k_wave*grid%dx))/(grid%dx*grid%dx)
      alpha = 1.0_wp - DT*NU*beta_sq

      ms%h_layer = 10.0_wp
      ms%v_face_y_layer = 0.0_wp
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx + 1
               ms%u_face_x_layer(i, j, k) = cos(k_wave*real(i - 1, wp)*grid%dx)
            end do
         end do
      end do

      call map_in(ms, hv, metrics, grid)
      do step = 1, N_STEPS
         call ocean_horizontal_viscosity_compute_tendencies(grid, metrics, hv, ms)
         call ocean_horizontal_viscosity_apply_tendencies(hv, ms, DT)
      end do
      call map_out(ms, hv, metrics)

      ! Probe far from the x-walls so the wall-frozen boundary value
      ! (i=1 and i=nx+1 are forced to zero tendency, holding cos(0)=1
      ! and cos(k*(nx)*dx)=cos(2pi)=1) doesn't contaminate the local
      ! decay.  Diffusive length over N_STEPS is sqrt(N*dt*nu) ~ 0.5
      ! cells, so any probe more than a few cells from a wall is
      ! effectively pure periodic decay.  Probe at the cos-trough
      ! (k_wave*x = pi) to maximise signal.
      i_probe = nx/2 + 1
      j_probe = ny/2
      ic_at_probe = cos(k_wave*real(i_probe - 1, wp)*grid%dx)
      amp_expected = ic_at_probe*alpha**N_STEPS
      amp_obs = ms%u_face_x_layer(i_probe, j_probe, NZ/2)
      rel_err = abs(amp_obs - amp_expected)/abs(amp_expected)

      call check(error, rel_err < 5.0e-3_wp, &
                 "sinusoid decay: realised amplitude off analytic by > 5e-3")

      call hv%destroy(); call ms%destroy()
   end subroutine test_sinusoid_decay

   subroutine test_ke_monotone(error)
      !! Non-trivial u, v pattern; kinetic energy must decrease
      !! step-over-step.  Viscosity is a dissipative operator —
      !! KE budget under pure Laplacian friction satisfies
      !!   dKE/dt = -integral(nu * |grad u|^2 + nu * |grad v|^2) <= 0.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_horizontal_viscosity_t) :: hv
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp), parameter :: DT = 0.05_wp
      real(wp), parameter :: NU = 0.5_wp
      integer, parameter :: N_STEPS = 8
      integer :: i, j, k, step, nx, ny
      real(wp) :: ke_prev, ke_now, ke_initial
      checks: block

         call make_grid(grid, 24, 20, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call hv%init(grid, nz_ml=NZ)
         hv%nu_h = NU
         nx = grid%nx_total
         ny = grid%ny_total

         ms%h_layer = 10.0_wp
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx + 1
                  ms%u_face_x_layer(i, j, k) = sin(2.0_wp*PI*real(i, wp)/real(nx, wp))* &
                                               cos(PI*real(j, wp)/real(ny, wp))
               end do
            end do
            do j = 1, ny + 1
               do i = 1, nx
                  ms%v_face_y_layer(i, j, k) = cos(2.0_wp*PI*real(i, wp)/real(nx, wp))* &
                                               sin(PI*real(j, wp)/real(ny, wp))
               end do
            end do
         end do

         ke_initial = sum(ms%u_face_x_layer**2) + sum(ms%v_face_y_layer**2)
         ke_prev = ke_initial
         call map_in(ms, hv, metrics, grid)
         do step = 1, N_STEPS
            call ocean_horizontal_viscosity_compute_tendencies(grid, metrics, hv, ms)
            call ocean_horizontal_viscosity_apply_tendencies(hv, ms, DT)
            !$acc update self(ms%u_face_x_layer, ms%v_face_y_layer)
            ke_now = sum(ms%u_face_x_layer**2) + sum(ms%v_face_y_layer**2)
            if (ke_now >= ke_prev) exit
            ke_prev = ke_now
         end do
         call map_out(ms, hv, metrics)

         call check(error, ke_now < ke_prev + 1.0e-12_wp, &
                    "KE not monotone: viscosity injected energy")
         if (allocated(error)) exit checks
         call check(error, ke_now < ke_initial, &
                    "KE final >= initial: viscosity did not dissipate")

      end block checks
      call hv%destroy(); call ms%destroy()
   end subroutine test_ke_monotone

   subroutine test_metric_reduction(error)
      !! T5 reduction gate (design §5).  On a SQUARE uniform grid the
      !! curvilinear FV Laplacian must reproduce the OLD decoupled
      !! 5-point form `nu_h·(Δ²ₓu/dx² + Δ²ᵧu/dy²)` to round-off.  The old
      !! kernel has been replaced, so the expected value is computed
      !! analytically here from the same velocity field with the old
      !! grouping.  A square grid keeps every ratio in the bundle equal
      !! to 1, so the only departure is FP reassociation of the 3-point
      !! stencil — bounded well below the apply tolerance.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_horizontal_viscosity_t) :: hv
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: L = 1234.5_wp   ! square uniform spacing
      real(wp), parameter :: NU = 7.5e3_wp
      real(wp) :: inv_dx2, inv_dy2, lap, expect, max_du, max_dv
      real(wp), allocatable :: du_ref(:, :, :), dv_ref(:, :, :)
      integer :: i, j, k, nx, ny
      checks: block

         call make_grid(grid, 14, 11, L, L)
         ms%nz_ml = NZ
         call ms%init(grid)
         call hv%init(grid, nz_ml=NZ)
         hv%nu_h = NU
         nx = grid%nx_total
         ny = grid%ny_total
         inv_dx2 = 1.0_wp/(L*L)
         inv_dy2 = 1.0_wp/(L*L)

         ms%h_layer = 10.0_wp
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx + 1
                  ms%u_face_x_layer(i, j, k) = &
                     sin(0.21_wp*real(i, wp))*cos(0.17_wp*real(j, wp)) + 0.3_wp*real(k, wp)
               end do
            end do
            do j = 1, ny + 1
               do i = 1, nx
                  ms%v_face_y_layer(i, j, k) = &
                     cos(0.13_wp*real(i, wp))*sin(0.29_wp*real(j, wp)) - 0.2_wp*real(k, wp)
               end do
            end do
         end do

         ! Reference: OLD decoupled 5-point Laplacian × nu_h, host-side.
         allocate (du_ref(nx + 1, ny, NZ), source=0.0_wp)
         allocate (dv_ref(nx, ny + 1, NZ), source=0.0_wp)
         do k = 1, NZ
            do j = 2, ny - 1
               do i = 2, nx
                  lap = (ms%u_face_x_layer(i + 1, j, k) - 2.0_wp*ms%u_face_x_layer(i, j, k) + &
                         ms%u_face_x_layer(i - 1, j, k))*inv_dx2 + &
                        (ms%u_face_x_layer(i, j + 1, k) - 2.0_wp*ms%u_face_x_layer(i, j, k) + &
                         ms%u_face_x_layer(i, j - 1, k))*inv_dy2
                  du_ref(i, j, k) = NU*lap
               end do
            end do
            do j = 2, ny
               do i = 2, nx - 1
                  lap = (ms%v_face_y_layer(i + 1, j, k) - 2.0_wp*ms%v_face_y_layer(i, j, k) + &
                         ms%v_face_y_layer(i - 1, j, k))*inv_dx2 + &
                        (ms%v_face_y_layer(i, j + 1, k) - 2.0_wp*ms%v_face_y_layer(i, j, k) + &
                         ms%v_face_y_layer(i, j - 1, k))*inv_dy2
                  dv_ref(i, j, k) = NU*lap
               end do
            end do
         end do

         call map_in(ms, hv, metrics, grid)
         call ocean_horizontal_viscosity_compute_tendencies(grid, metrics, hv, ms)
         !$acc update self(hv%du_visc%data, hv%dv_visc%data)
         call map_out(ms, hv, metrics)

         ! Compare interior tendencies (the boundary rows are zeroed by
         ! both forms identically).
         max_du = 0.0_wp
         max_dv = 0.0_wp
         do k = 1, NZ
            do j = 2, ny - 1
               do i = 2, nx
                  max_du = max(max_du, abs(hv%du_visc%data(i, j, k) - du_ref(i, j, k)))
               end do
            end do
            do j = 2, ny
               do i = 2, nx - 1
                  max_dv = max(max_dv, abs(hv%dv_visc%data(i, j, k) - dv_ref(i, j, k)))
               end do
            end do
         end do

         ! Departure is FP reassociation only — scale by nu_h·|u|/dx².
         expect = NU/(L*L)
         call check(error, max_du < 1.0e-12_wp*expect, &
                    "T5: u metric Laplacian differs from scalar form beyond round-off")
         if (allocated(error)) exit checks
         call check(error, max_dv < 1.0e-12_wp*expect, &
                    "T5: v metric Laplacian differs from scalar form beyond round-off")

      end block checks
      if (allocated(du_ref)) deallocate (du_ref)
      if (allocated(dv_ref)) deallocate (dv_ref)
      call hv%destroy(); call ms%destroy()
   end subroutine test_metric_reduction

   ! -----------------------------------------------------------------
   ! MOM6 thickness-weighted stress-divergence path (spec PR2)
   ! -----------------------------------------------------------------

   subroutine test_stress_reduces_to_laplacian(error)
      !! Bit-identity hinge (spec test a): on a uniform square grid +
      !! uniform layer thickness + all-wet domain the MOM6 stress
      !! operator must reproduce the velocity-Laplacian tendency to
      !! round-off.  Run the SAME field through both paths and compare
      !! the interior `du_visc` / `dv_visc`.  The cross terms (`v_xy`)
      !! cancel discretely and the layer thickness cancels in the
      !! divide, so the only departure is FP reassociation.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_horizontal_viscosity_t) :: hv_lap, hv_m6
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: L = 1234.5_wp, NU = 7.5e3_wp, DT = 1.0e-3_wp
      real(wp), allocatable :: du_lap(:, :, :), dv_lap(:, :, :)
      real(wp) :: max_du, max_dv, expect
      integer :: i, j, k, nx, ny
      checks: block

         call make_grid(grid, 14, 11, L, L)
         ms%nz_ml = NZ
         call ms%init(grid)
         call hv_lap%init(grid, nz_ml=NZ)
         call hv_m6%init(grid, nz_ml=NZ)
         hv_lap%nu_h = NU
         hv_m6%nu_h = NU
         hv_m6%stress_tensor = .true.
         hv_m6%bound_coef = 0.8_wp     ! huge dt ⇒ clamp inactive
         nx = grid%nx_total
         ny = grid%ny_total

         ! Large uniform h so the `h_neglect` floor (H_VANISHED ~ 1.5e-4)
         ! is negligible vs the layer thickness — the reduction to the
         ! velocity Laplacian is exact only in the h ≫ h_neglect limit
         ! (MOM6 carries the same floor).  h = 1e6 m ⇒ h_neglect/h ~ 1e-10.
         ms%h_layer = 1.0e6_wp
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx + 1
                  ms%u_face_x_layer(i, j, k) = &
                     sin(0.21_wp*real(i, wp))*cos(0.17_wp*real(j, wp)) + 0.3_wp*real(k, wp)
               end do
            end do
            do j = 1, ny + 1
               do i = 1, nx
                  ms%v_face_y_layer(i, j, k) = &
                     cos(0.13_wp*real(i, wp))*sin(0.29_wp*real(j, wp)) - 0.2_wp*real(k, wp)
               end do
            end do
         end do

         call make_cartesian_metrics(metrics, grid)
         !$acc enter data copyin(ms)
         call ms%enter_data()
         !$acc enter data copyin(hv_lap)
         call hv_lap%enter_data()
         !$acc enter data copyin(hv_m6)
         call hv_m6%enter_data()

         call ocean_horizontal_viscosity_compute_tendencies(grid, metrics, hv_lap, ms)
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
            do j = 2, ny - 1
               do i = 2, nx
                  max_du = max(max_du, abs(hv_m6%du_visc%data(i, j, k) - du_lap(i, j, k)))
               end do
            end do
            do j = 2, ny
               do i = 2, nx - 1
                  max_dv = max(max_dv, abs(hv_m6%dv_visc%data(i, j, k) - dv_lap(i, j, k)))
               end do
            end do
         end do

         ! Departure is FP reassociation only — scale by nu_h·|u|/dx².
         expect = NU/(L*L)
         call check(error, max_du < 1.0e-10_wp*expect, &
                    "stress_tensor: u tendency differs from Laplacian beyond round-off")
         if (allocated(error)) exit checks
         call check(error, max_dv < 1.0e-10_wp*expect, &
                    "stress_tensor: v tendency differs from Laplacian beyond round-off")

      end block checks
      if (allocated(du_lap)) deallocate (du_lap)
      if (allocated(dv_lap)) deallocate (dv_lap)
      call hv_m6%destroy(); call hv_lap%destroy(); call ms%destroy()
   end subroutine test_stress_reduces_to_laplacian

   subroutine test_stress_momentum_conservation(error)
      !! Spec test b: the stress-divergence operator is
      !! momentum-conserving on a closed box.  The thickness-weighted
      !! momentum tendency `areaCu·h_u·diffu` is the discrete
      !! divergence of an interior stress flux; with a velocity field
      !! of COMPACT SUPPORT (zero in a band next to every wall) all the
      !! boundary stresses vanish, so the interior fluxes telescope and
      !! the total layer x- and y-momentum change must each vanish to
      !! round-off.  (A closed wall otherwise exerts a tangential
      !! free-slip stress, so compact support is what isolates the
      !! conservation property — equivalent to a periodic box with no
      !! flux through the seam.)  Large uniform h keeps the `h_neglect`
      !! floor negligible.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_horizontal_viscosity_t) :: hv
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp), parameter :: L = 5.0_wp, NU = 100.0_wp, DT = 1.0e-4_wp, H0 = 1.0e6_wp
      real(wp) :: sum_u, sum_v, scale, area, bump
      integer :: i, j, k, nx, ny
      checks: block

         call make_grid(grid, 24, 22, L, L)
         ms%nz_ml = NZ
         call ms%init(grid)
         call hv%init(grid, nz_ml=NZ)
         hv%nu_h = NU
         hv%stress_tensor = .true.
         nx = grid%nx_total
         ny = grid%ny_total
         area = L*L

         ! Compactly-supported velocity bump: smooth sin² hump centred
         ! in the interior, identically zero within 5 cells of every
         ! wall.  No stress reaches a boundary ⇒ pure interior flux
         ! divergence.
         ms%h_layer = H0
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         do k = 1, NZ
            do j = 8, ny - 7
               do i = 8, nx - 6
                  bump = sin(PI*real(i - 7, wp)/8.0_wp)**2*sin(PI*real(j - 7, wp)/8.0_wp)**2
                  ms%u_face_x_layer(i, j, k) = bump
               end do
            end do
            do j = 8, ny - 6
               do i = 8, nx - 7
                  bump = sin(PI*real(i - 7, wp)/8.0_wp)**2*sin(PI*real(j - 7, wp)/8.0_wp)**2
                  ms%v_face_y_layer(i, j, k) = 0.7_wp*bump
               end do
            end do
         end do

         call make_cartesian_metrics(metrics, grid)
         !$acc enter data copyin(ms)
         call ms%enter_data()
         !$acc enter data copyin(hv)
         call hv%enter_data()
         call ocean_horizontal_viscosity_compute_tendencies(grid, metrics, hv, ms, dt=DT)
         !$acc update self(hv%du_visc%data, hv%dv_visc%data)
         call hv%exit_data()
         !$acc exit data delete(hv)
         call ms%exit_data()
         !$acc exit data delete(ms)
         call destroy_cartesian_metrics(metrics)

         ! Thickness-weighted momentum change per face = h_u·area·diffu.
         ! Uniform h ⇒ h_u = H0; sum over all faces (wall faces are 0).
         ! `scale` = total gross momentum tendency (sum of |·|) — the
         ! residual is compared relative to it.
         sum_u = 0.0_wp
         sum_v = 0.0_wp
         scale = 0.0_wp
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx + 1
                  sum_u = sum_u + H0*area*hv%du_visc%data(i, j, k)
                  scale = scale + H0*area*abs(hv%du_visc%data(i, j, k))
               end do
            end do
            do j = 1, ny + 1
               do i = 1, nx
                  sum_v = sum_v + H0*area*hv%dv_visc%data(i, j, k)
                  scale = scale + H0*area*abs(hv%dv_visc%data(i, j, k))
               end do
            end do
         end do

         call check(error, abs(sum_u) < 1.0e-12_wp*scale, &
                    "stress_tensor: x-momentum not conserved over closed box")
         if (allocated(error)) exit checks
         call check(error, abs(sum_v) < 1.0e-12_wp*scale, &
                    "stress_tensor: y-momentum not conserved over closed box")

      end block checks
      call hv%destroy(); call ms%destroy()
   end subroutine test_stress_momentum_conservation

   subroutine test_stress_cfl_clamp(error)
      !! Spec test c: a viscosity coefficient that exceeds the explicit
      !! per-cell CFL bound is clamped, so the field stays finite even
      !! at a near-CFL time step.  Pick `nu_h` and `dt` so the
      !! unclamped viscous update would blow up
      !! (`nu·dt·(1/dx²+1/dy²) ≫ 0.5`); the clamp must hold the apply
      !! step bounded.  Compare against the same config WITHOUT the
      !! clamp budget (tiny bound_coef vs the value) to confirm the
      !! clamp is what tames it.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_horizontal_viscosity_t) :: hv
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp), parameter :: L = 1.0_wp, DT = 0.4_wp
      real(wp), parameter :: NU = 1.0e6_wp   ! wildly over the CFL bound
      integer, parameter :: N_STEPS = 40
      real(wp) :: max_abs
      integer :: i, j, k, step, nx, ny
      checks: block

         call make_grid(grid, 20, 18, L, L)
         ms%nz_ml = NZ
         call ms%init(grid)
         call hv%init(grid, nz_ml=NZ)
         hv%nu_h = NU
         hv%stress_tensor = .true.
         hv%bound_coef = 0.8_wp
         nx = grid%nx_total
         ny = grid%ny_total

         ms%h_layer = 10.0_wp
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx + 1
                  ms%u_face_x_layer(i, j, k) = sin(2.0_wp*PI*real(i, wp)/real(nx, wp))
               end do
            end do
            do j = 1, ny + 1
               do i = 1, nx
                  ms%v_face_y_layer(i, j, k) = cos(2.0_wp*PI*real(j, wp)/real(ny, wp))
               end do
            end do
         end do

         call make_cartesian_metrics(metrics, grid)
         !$acc enter data copyin(ms)
         call ms%enter_data()
         !$acc enter data copyin(hv)
         call hv%enter_data()
         do step = 1, N_STEPS
            call ocean_horizontal_viscosity_compute_tendencies(grid, metrics, hv, ms, dt=DT)
            call ocean_horizontal_viscosity_apply_tendencies(hv, ms, DT)
         end do
         !$acc update self(ms%u_face_x_layer, ms%v_face_y_layer)
         call hv%exit_data()
         !$acc exit data delete(hv)
         call ms%exit_data()
         !$acc exit data delete(ms)
         call destroy_cartesian_metrics(metrics)

         max_abs = max(maxval(abs(ms%u_face_x_layer)), maxval(abs(ms%v_face_y_layer)))
         ! With the CFL clamp the field can only DECAY: |u| <= IC peak (1).
         ! Without it, NU·dt/dx² ~ 4e5 ⇒ exponential blow-up in 1 step.
         call check(error, max_abs < 1.5_wp, &
                    "stress_tensor: CFL clamp failed — field grew past the IC amplitude")
         if (allocated(error)) exit checks
         call check(error, max_abs == max_abs, &   ! NaN guard
                    "stress_tensor: CFL clamp failed — field went NaN")

      end block checks
      call hv%destroy(); call ms%destroy()
   end subroutine test_stress_cfl_clamp

   subroutine test_stress_coast_mask(error)
      !! Spec test d: a block of land must not exchange viscous
      !! momentum with the surrounding ocean.  Seed a land block,
      !! derive the wet masks, set a strong velocity in the wet domain
      !! and zero velocity inside the land block, and confirm the
      !! stress operator produces ZERO tendency on the land-interior
      !! faces (no leakage across the coast).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_horizontal_viscosity_t) :: hv
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: L = 1.0_wp, NU = 50.0_wp, DT = 1.0e-6_wp
      real(wp), allocatable :: wet_mask(:, :)
      real(wp) :: max_land_du
      integer :: i, j, k, nx, ny
      integer :: il0, il1, jl0, jl1
      checks: block

         call make_grid(grid, 16, 14, L, L)
         ms%nz_ml = NZ
         call ms%init(grid)
         call hv%init(grid, nz_ml=NZ)
         hv%nu_h = NU
         hv%stress_tensor = .true.
         nx = grid%nx_total
         ny = grid%ny_total

         ! Land block (T-cell index range, interior, away from walls).
         il0 = 6; il1 = 9
         jl0 = 5; jl1 = 8

         ms%h_layer = 10.0_wp
         ! Strong, non-uniform wet-domain flow; zero inside land.
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx + 1
                  ms%u_face_x_layer(i, j, k) = 0.7_wp*real(i, wp) - 0.3_wp*real(j, wp)
               end do
            end do
            do j = 1, ny + 1
               do i = 1, nx
                  ms%v_face_y_layer(i, j, k) = -0.4_wp*real(i, wp) + 0.6_wp*real(j, wp)
               end do
            end do
         end do
         ! Zero velocity on faces interior to the land block.
         do k = 1, NZ
            do j = jl0, jl1
               do i = il0 + 1, il1
                  ms%u_face_x_layer(i, j, k) = 0.0_wp
               end do
            end do
            do j = jl0 + 1, jl1
               do i = il0, il1
                  ms%v_face_y_layer(i, j, k) = 0.0_wp
               end do
            end do
         end do

         ! Build masked metrics by hand (apply_land_mask must run before
         ! enter_data).
         call metrics%init(grid)
         call metrics_fill_cartesian(metrics, grid, grid%dx, grid%dy)
         call metrics_finalize(metrics)
         allocate (wet_mask(nx, ny), source=1.0_wp)
         do j = jl0, jl1
            do i = il0, il1
               wet_mask(i, j) = 0.0_wp
            end do
         end do
         call metrics_apply_land_mask(metrics, wet_mask, grid, .false., .false., .false.)
         !$acc enter data copyin(metrics)
         call metrics%enter_data()

         !$acc enter data copyin(ms)
         call ms%enter_data()
         !$acc enter data copyin(hv)
         call hv%enter_data()
         call ocean_horizontal_viscosity_compute_tendencies(grid, metrics, hv, ms, dt=DT)
         !$acc update self(hv%du_visc%data)
         call hv%exit_data()
         !$acc exit data delete(hv)
         call ms%exit_data()
         !$acc exit data delete(ms)
         call destroy_cartesian_metrics(metrics)

         ! u-faces strictly INSIDE the land block (both adjacent T-cells
         ! are land ⇒ wet_u==0) must carry zero viscous tendency — no
         ! momentum leaks across the coast into the dead interior.
         max_land_du = 0.0_wp
         do k = 1, NZ
            do j = jl0, jl1
               do i = il0 + 1, il1
                  max_land_du = max(max_land_du, abs(hv%du_visc%data(i, j, k)))
               end do
            end do
         end do

         call check(error, max_land_du < 1.0e-12_wp, &
                    "stress_tensor: viscous momentum leaked into the land interior")

      end block checks
      if (allocated(wet_mask)) deallocate (wet_mask)
      call hv%destroy(); call ms%destroy()
   end subroutine test_stress_coast_mask

end module test_ocean_hvisc
