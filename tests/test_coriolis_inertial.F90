!! Unit tests for the C-grid Coriolis force kernel
!! (rdb_coriolis_adv Phase 3a — plain Coriolis on the barotropic
!! momentum equation; PV-conserving Sadourny + Hollingsworth-Källén
!! land in Phase 3b/c).  See `docs/ROADMAP_OCEAN.md` Phase 3.
!!
!! Cases:
!!   * Zero-velocity null check — uniform f, zero u and v, the
!!     Coriolis step must leave the velocity field at zero
!!     bit-for-bit.  Catches sign errors in unused branches.
!!   * Inertial oscillation — uniform f, uniform u, zero v.  After
!!     N forward-Euler Coriolis steps the velocity vector must
!!     rotate by f*N*dt within a 5% magnitude tolerance and ±5%
!!     angle tolerance (FE has a small per-step amplitude growth
!!     and negligible phase error for f*dt << 1).
module test_coriolis_inertial
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_barotropic_state, only: barotropic_state_t
   use rdb_coriolis_adv, only: coriolis_adv_t, &
                               coriolis_adv_compute_tendencies_barotropic, &
                               coriolis_adv_apply_tendencies_barotropic
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   implicit none
   private

   public :: collect_coriolis_inertial_tests

   integer, parameter :: NX_PHYS = 16
   integer, parameter :: NY_PHYS = 16
   integer, parameter :: NGHOST = 2

contains

   subroutine collect_coriolis_inertial_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("zero_velocity_null_step", test_zero_velocity), &
                  new_unittest("inertial_oscillation_quarter_period", &
                               test_inertial_oscillation), &
                  new_unittest("sadourny_solid_body_rotation", &
                               test_sadourny_solid_body) &
                  ]
   end subroutine collect_coriolis_inertial_tests

   subroutine make_grid(grid)
      type(hgrid_t), intent(out) :: grid
      call grid%init(NX_PHYS, NY_PHYS, NGHOST, 1.0_wp, 1.0_wp)
   end subroutine make_grid

   subroutine map_in(grid, bs, cor, metrics)
      type(hgrid_t), intent(in) :: grid
      type(barotropic_state_t), intent(inout) :: bs
      type(coriolis_adv_t), intent(inout) :: cor
      type(ocean_metrics_t), intent(inout) :: metrics
      call make_cartesian_metrics(metrics, grid)
      !$acc enter data copyin(bs)
      call bs%enter_data()
      !$acc enter data copyin(cor)
      call cor%enter_data()
   end subroutine map_in

   subroutine map_out(bs, cor, metrics)
      type(barotropic_state_t), intent(inout) :: bs
      type(coriolis_adv_t), intent(inout) :: cor
      type(ocean_metrics_t), intent(inout) :: metrics
      call cor%exit_data()
      !$acc exit data delete(cor)
      call bs%exit_data()
      !$acc exit data delete(bs)
      call destroy_cartesian_metrics(metrics)
   end subroutine map_out

   ! -----------------------------------------------------------------
   ! Cases
   ! -----------------------------------------------------------------

   subroutine test_zero_velocity(error)
      !! Zero velocity in -> zero tendency out -> zero velocity
      !! after the step.  Confirms the kernel doesn't inject
      !! spurious velocity from any uninitialized branch.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(barotropic_state_t) :: bs
      type(coriolis_adv_t) :: cor
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: F_C = 1.0_wp
      real(wp), parameter :: DT = 0.01_wp
      real(wp) :: max_u, max_v

      call make_grid(grid)
      call bs%init(grid)
      cor%f_0 = F_C
      call cor%init(grid)
      bs%h = 1.0_wp
      bs%u_face_x = 0.0_wp
      bs%v_face_y = 0.0_wp

      call map_in(grid, bs, cor, metrics)
      call coriolis_adv_compute_tendencies_barotropic(grid, metrics, cor, bs)
      call coriolis_adv_apply_tendencies_barotropic(cor, bs, DT)
      call map_out(bs, cor, metrics)

      max_u = maxval(abs(bs%u_face_x))
      max_v = maxval(abs(bs%v_face_y))
      call check(error, max_u < 1.0e-14_wp, &
                 "Coriolis step injected u from zero initial state")
      if (allocated(error)) then
         call cor%destroy(); call bs%destroy(); return
      end if
      call check(error, max_v < 1.0e-14_wp, &
                 "Coriolis step injected v from zero initial state")

      call cor%destroy()
      call bs%destroy()
   end subroutine test_zero_velocity

   subroutine test_inertial_oscillation(error)
      !! Inertial oscillation over a quarter period.  Initial state:
      !! uniform u = U0, v = 0.  Continuous solution:
      !!   u(t) =  U0 * cos(f*t)
      !!   v(t) = -U0 * sin(f*t)
      !! At t = pi/(2f), u -> 0 and v -> -U0.
      !!
      !! Forward Euler is unstable for the oscillator (eigenvalues
      !! 1 ± i*f*dt have magnitude > 1), so we set f*dt = 0.01 to
      !! keep per-step amplitude growth at sqrt(1 + (f*dt)^2) ≈ 1 +
      !! 5e-5.  Over N=157 steps that's a ~0.8% amplitude rise.
      !! Tolerances allow up to 2% magnitude error and ±3% angle
      !! error (the latter is essentially zero for FE at f*dt << 1).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(barotropic_state_t) :: bs
      type(coriolis_adv_t) :: cor
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: F_C = 1.0_wp
      real(wp), parameter :: U0 = 1.0_wp
      real(wp), parameter :: DT = 0.01_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      integer, parameter :: N_STEPS = 157   ! ~ pi/(2*F_C*DT)
      real(wp) :: t_end, u_expected, v_expected
      real(wp) :: u_obs, v_obs, mag_obs, mag_expected
      real(wp) :: mag_err, angle_obs, angle_expected, angle_err
      integer :: step, i_probe, j_probe, nx, ny

      call make_grid(grid)
      call bs%init(grid)
      cor%f_0 = F_C
      call cor%init(grid)
      bs%h = 1.0_wp
      nx = grid%nx_total
      ny = grid%ny_total

      bs%u_face_x = U0
      bs%v_face_y = 0.0_wp

      call map_in(grid, bs, cor, metrics)
      do step = 1, N_STEPS
         call coriolis_adv_compute_tendencies_barotropic(grid, metrics, cor, bs)
         call coriolis_adv_apply_tendencies_barotropic(cor, bs, DT)
      end do
      call map_out(bs, cor, metrics)

      ! Probe an interior face well away from any wall stagger
      ! averaging asymmetry — middle of the domain.
      i_probe = nx/2
      j_probe = ny/2
      u_obs = bs%u_face_x(i_probe, j_probe)
      v_obs = bs%v_face_y(i_probe, j_probe)

      t_end = real(N_STEPS, wp)*DT
      u_expected = U0*cos(F_C*t_end)
      v_expected = -U0*sin(F_C*t_end)
      mag_obs = sqrt(u_obs*u_obs + v_obs*v_obs)
      mag_expected = sqrt(u_expected*u_expected + v_expected*v_expected)
      mag_err = abs(mag_obs - mag_expected)/mag_expected

      angle_obs = atan2(-v_obs, u_obs)   ! Angle rotated CCW from +u in (u, -v)
      angle_expected = F_C*t_end
      angle_err = abs(angle_obs - angle_expected)/angle_expected

      call check(error, mag_err < 0.02_wp, &
                 "inertial oscillation magnitude off by > 2%")
      if (allocated(error)) then
         call cor%destroy(); call bs%destroy(); return
      end if

      call check(error, angle_err < 0.03_wp, &
                 "inertial oscillation phase off by > 3%")

      call cor%destroy()
      call bs%destroy()
   end subroutine test_inertial_oscillation

   subroutine test_sadourny_solid_body(error)
      !! Solid-body rotation about cell-centre (i_c, j_c) provides a
      !! clean signal that distinguishes Sadourny from plain Coriolis:
      !!
      !!   u_face_x(i, j) = -omega * (j - j_c)
      !!   v_face_y(i, j) = +omega * (i - i_c)
      !!
      !! With this IC the discrete corner-vorticity is uniform
      !! zeta = 2*omega and the discrete KE at cell centres reduces
      !! to (omega^2/2) * ((i-i_c)^2 + (j-j_c)^2).  The Sadourny
      !! tendency at a v-face displaced (j - j_c - 1/2) cells from
      !! the rotation centre is
      !!
      !!   dv/dt|_Sadourny = (omega^2 + f*omega) * (j - j_c - 1/2)
      !!
      !! whereas plain Coriolis (Phase 3a) gives only
      !!
      !!   dv/dt|_plain = omega * (j - j_c - 1/2)
      !!
      !! For omega = f = 1, Sadourny is exactly 2x plain.  The test
      !! probes a v-face well inside the domain (so the corner
      !! stencil sees the full 4 cells and the zeta-fallback at
      !! outer corners doesn't pollute the result) and asserts:
      !!   (a) observed dv/dt matches Sadourny to ~1e-12
      !!   (b) observed dv/dt is far from plain Coriolis (factor 2),
      !!       i.e. the kernel actually does Sadourny.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(barotropic_state_t) :: bs
      type(coriolis_adv_t) :: cor
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: OMEGA = 1.0_wp
      real(wp), parameter :: F_C = 1.0_wp
      real(wp), parameter :: DT = 0.01_wp
      integer, parameter :: I_PROBE_OFFSET = 3
      integer, parameter :: J_PROBE_OFFSET = 5
      integer :: i, j, nx, ny, i_c, j_c, i_probe, j_probe
      real(wp) :: r_y, v_initial, v_final, dv_dt_obs
      real(wp) :: dv_dt_expected_sadourny, dv_dt_expected_plain

      call make_grid(grid)
      call bs%init(grid)
      cor%f_0 = F_C
      call cor%init(grid)
      bs%h = 1.0_wp
      nx = grid%nx_total
      ny = grid%ny_total
      i_c = nx/2
      j_c = ny/2

      ! Solid-body rotation IC
      do j = 1, ny
         do i = 1, nx + 1
            bs%u_face_x(i, j) = -OMEGA*(real(j, wp) - real(j_c, wp))
         end do
      end do
      do j = 1, ny + 1
         do i = 1, nx
            bs%v_face_y(i, j) = OMEGA*(real(i, wp) - real(i_c, wp))
         end do
      end do

      ! Probe a v-face inside the domain, far from the boundary
      ! zeta-fallback ring (which is the 1-cell-wide rim).
      i_probe = i_c + I_PROBE_OFFSET
      j_probe = j_c + J_PROBE_OFFSET
      v_initial = bs%v_face_y(i_probe, j_probe)

      ! Discrete v-face position offset from rotation centre.
      r_y = real(j_probe, wp) - real(j_c, wp) - 0.5_wp

      ! Sadourny (Phase 3b) and plain-Coriolis (Phase 3a) references
      dv_dt_expected_sadourny = r_y*(OMEGA**2 + F_C*OMEGA)
      dv_dt_expected_plain = OMEGA*r_y

      call map_in(grid, bs, cor, metrics)
      call coriolis_adv_compute_tendencies_barotropic(grid, metrics, cor, bs)
      call coriolis_adv_apply_tendencies_barotropic(cor, bs, DT)
      call map_out(bs, cor, metrics)

      v_final = bs%v_face_y(i_probe, j_probe)
      dv_dt_obs = (v_final - v_initial)/DT

      ! Sadourny match: tight, the discrete stencil reproduces the
      ! continuous formula exactly for the linear u/v + quadratic KE
      ! of a solid-body rotation.
      call check(error, abs(dv_dt_obs - dv_dt_expected_sadourny) < 1.0e-10_wp, &
                 "kernel dv/dt deviates from Sadourny analytic")
      if (allocated(error)) then
         call cor%destroy(); call bs%destroy(); return
      end if
      ! Not plain Coriolis: must differ by the full vorticity-flux
      ! + KE-gradient contribution.
      call check(error, &
                 abs(dv_dt_obs - dv_dt_expected_plain) > 0.5_wp* &
                 abs(dv_dt_expected_sadourny - dv_dt_expected_plain), &
                 "kernel dv/dt matched plain Coriolis (Sadourny terms missing)")

      call cor%destroy()
      call bs%destroy()
   end subroutine test_sadourny_solid_body

end module test_coriolis_inertial
