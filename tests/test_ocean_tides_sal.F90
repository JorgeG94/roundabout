!! Scalar self-attraction & loading (C2) test on the C1 tide seam.
!!
!! Two independent checks:
!!  (a) DIRECT relationship — with `use_sal=.true., beta_sal=0.1` and a
!!      known lagged SSH, the combined seam field satisfies
!!      `eta_forcing = eta_eq + 0.1*eta` exactly (the scalar-SAL fold);
!!  (b) EFFECTIVE-GRAVITY shift — the closed-basin M2 barotropic response
!!      slows with SAL on: `eta_sal = beta*eta` makes the surface term
!!      `-g(1-beta) grad(eta)`, so the free-surface gravity-wave speed
!!      drops by sqrt(1-beta) and the near-equilibrium SSH probe's RMS
!!      amplitude changes measurably versus `beta=0` (bit-identical to C1);
!!  (c) OFF ⇒ `eta_forcing` bit-identical to `eta_eq` (absent-SAL path).
module test_ocean_tides_sal
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_spherical_metrics
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
   use rdb_ocean_tides, only: ocean_tides_t, tides_configure_astronomy, &
                              tides_build_struct, tides_update_eta_eq, &
                              tides_update_eta_sal
   use rdb_ocean_tide_astro, only: TIDE_OMEGA, TIDE_AMP, TIDE_LOVE, TIDE_DEG2RAD
   implicit none
   private

   public :: collect_ocean_tides_sal_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NX_PHYS = 24, NY_PHYS = 24, NZ = 2
   real(wp), parameter :: LAT0 = 29.0_wp
   real(wp), parameter :: LON0 = 175.0_wp
   real(wp), parameter :: DLON = 0.5_wp
   real(wp), parameter :: DLAT = 0.1_wp
   real(wp), parameter :: PHI_MID = 30.2_wp
   real(wp), parameter :: RAD_EARTH = 6.371e6_wp
   real(wp), parameter :: H_TOTAL = 4000.0_wp
   real(wp), parameter :: DT = 300.0_wp
   integer, parameter :: N_INNER = 16
   real(wp), parameter :: BETA = 0.1_wp

contains

   subroutine collect_ocean_tides_sal_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("sal_direct_and_offident", test_sal_direct), &
                  new_unittest("sal_shifts_barotropic_response", test_sal_shift)]
   end subroutine collect_ocean_tides_sal_tests

   ! (a)+(c): unit-level check of the eta_forcing fold, no time-stepping.
   subroutine test_sal_direct(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ocean_tides_t) :: tides
      integer :: cat_idx(1), nx, ny, i, j
      real(wp), allocatable :: eta_ssh(:, :)
      real(wp) :: maxdev, expct, sref

      call grid%init(NX_PHYS, NY_PHYS, NGHOST, DLON, DLAT)
      nx = grid%nx_total
      ny = grid%ny_total
      call make_spherical_metrics(metrics, grid, LON0, LAT0, DLON, DLAT, RAD_EARTH)

      cat_idx(1) = 1     ! M2
      tides%enable = .true.
      call tides_configure_astronomy(tides, cat_idx, 1, 0.0_wp, 0.0_wp, .false., nx, ny)
      call tides_build_struct(tides, metrics%geolatT, metrics%geolonT, nx, ny)

      ! A known, spatially varying lagged SSH.
      allocate (eta_ssh(nx, ny))
      do j = 1, ny
         do i = 1, nx
            eta_ssh(i, j) = 0.01_wp*real(i, wp) - 0.02_wp*real(j, wp)
         end do
      end do

      ! Device-resident: the tide arrays + eta_ssh must be mapped for the
      ! stdpar=gpu kernels (the `update device` in tides_update_eta_eq and
      ! the fold in tides_update_eta_sal both run on device).
      call tides%enter_data()
      !$acc enter data copyin(eta_ssh)

      call tides_update_eta_eq(tides, 0.0_wp)

      ! --- (c) OFF: eta_forcing must be bit-identical to eta_eq. ---
      tides%use_sal = .false.
      tides%beta_sal = BETA           ! non-zero beta must be IGNORED when off
      call tides_update_eta_sal(tides, eta_ssh)
      !$acc update self(tides%eta_forcing, tides%eta_eq)
      maxdev = maxval(abs(tides%eta_forcing - tides%eta_eq))
      call check(error, maxdev == 0.0_wp, &
                 "use_sal=.false. must give eta_forcing bit-identical to eta_eq")
      if (allocated(error)) then
         !$acc exit data delete(eta_ssh)
         call tides%exit_data()
         call metrics%exit_data()
         deallocate (eta_ssh)
         call tides%destroy()
         return
      end if

      ! --- (a) ON physics relationship: eta_sal = beta*eta and
      ! eta_forcing = eta_eq + eta_sal.  A tight tolerance (not bit-exact):
      ! the reference `BETA*eta_ssh` re-contracts to an FMA against the
      ! kernel's stored product, a sub-ULP (~3e-18) comparison artefact on
      ! a floating multiply — physically irrelevant.  (The load-bearing
      ! bit-identity guarantee is the use_sal=.false. check above.)
      tides%use_sal = .true.
      tides%beta_sal = BETA
      call tides_update_eta_sal(tides, eta_ssh)
      !$acc update self(tides%eta_forcing, tides%eta_sal, tides%eta_eq)
      maxdev = 0.0_wp
      do j = 1, ny
         do i = 1, nx
            sref = BETA*eta_ssh(i, j)
            maxdev = max(maxdev, abs(tides%eta_sal(i, j) - sref))
            expct = tides%eta_eq(i, j) + tides%eta_sal(i, j)
            maxdev = max(maxdev, abs(tides%eta_forcing(i, j) - expct))
         end do
      end do
      call check(error, maxdev < 1.0e-14_wp, &
                 "eta_forcing must equal eta_eq + beta_sal*eta")

      !$acc exit data delete(eta_ssh)
      call tides%exit_data()
      call metrics%exit_data()
      deallocate (eta_ssh)
      call tides%destroy()
   end subroutine test_sal_direct

   ! (b): end-to-end — SAL on changes the barotropic response magnitude.
   subroutine test_sal_shift(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: rms_off, rms_on
      call run_response(error, .false., 0.0_wp, rms_off)
      if (allocated(error)) return
      call run_response(error, .true., BETA, rms_on)
      if (allocated(error)) return

      ! Both live + finite.
      call check(error, rms_off > 1.0e-4_wp .and. rms_on > 1.0e-4_wp, &
                 "SSH response inert in one of the runs")
      if (allocated(error)) return
      ! SAL folds beta*eta into the same PGF seam => the effective gravity
      ! drops (1-beta), a physically distinct barotropic response.  Assert a
      ! measurable, non-trivial change (>0.5% RMS) between the two runs.
      call check(error, abs(rms_on - rms_off)/rms_off > 5.0e-3_wp, &
                 "scalar SAL did not shift the barotropic response")
   end subroutine test_sal_shift

   subroutine run_response(error, use_sal, beta, rms)
      type(error_type), allocatable, intent(out) :: error
      logical, intent(in) :: use_sal
      real(wp), intent(in) :: beta
      real(wp), intent(out) :: rms
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
      type(ocean_tides_t) :: tides
      integer :: cat_idx(1)
      integer :: i, j, k, ip, jp, step, n_steps
      integer :: nx, ny
      real(wp) :: t, t_end, period_m2, ssh_probe, sumsq
      integer :: nsamp

      period_m2 = 2.0_wp*acos(-1.0_wp)/TIDE_OMEGA(1)
      t_end = 2.2_wp*period_m2
      n_steps = int(t_end/DT) + 1

      call grid%init(NX_PHYS, NY_PHYS, NGHOST, DLON, DLAT)
      nx = grid%nx_total
      ny = grid%ny_total
      ip = grid%nghost + NX_PHYS/4
      jp = grid%nghost + NY_PHYS/2

      call make_spherical_metrics(metrics, grid, LON0, LAT0, DLON, DLAT, RAD_EARTH)

      ms%nz_ml = NZ
      call ms%init(grid)
      call ct%init(grid, nz_ml=NZ)
      cor%f_0 = 2.0_wp*7.292115e-5_wp*sin(PHI_MID*TIDE_DEG2RAD)
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

      ms%h_layer = H_TOTAL/real(NZ, wp)
      ms%u_face_x_layer = 0.0_wp
      ms%v_face_y_layer = 0.0_wp
      do k = 1, NZ
         ms%tracers(ms%idx_salinity)%hTr(:, :, k) = eos%S_ref*(H_TOTAL/real(NZ, wp))
         ms%tracers(ms%idx_temperature)%hTr(:, :, k) = eos%T_ref*(H_TOTAL/real(NZ, wp))
      end do
      dyn%bt_work%bt_H_ref = H_TOTAL

      cat_idx(1) = 1     ! M2
      tides%enable = .true.
      tides%use_sal = use_sal
      tides%beta_sal = beta
      call tides_configure_astronomy(tides, cat_idx, 1, 0.0_wp, 0.0_wp, .false., nx, ny)
      call tides_build_struct(tides, metrics%geolatT, metrics%geolonT, nx, ny)

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
      call tides%enter_data()

      t = 0.0_wp
      sumsq = 0.0_wp
      nsamp = 0
      do step = 1, n_steps
         call ocean_dyn_step_split(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
                                   va, hd, vd, vmix, ms, DT, N_INNER, t=t, tides=tides)
         t = t + DT
         !$acc update self(ms%h_layer)
         ssh_probe = 0.0_wp
         do k = 1, NZ
            ssh_probe = ssh_probe + ms%h_layer(ip, jp, k)
         end do
         ssh_probe = ssh_probe - H_TOTAL
         if (step >= n_steps/4) then
            sumsq = sumsq + ssh_probe*ssh_probe
            nsamp = nsamp + 1
         end if
      end do
      rms = sqrt(sumsq/real(max(nsamp, 1), wp))

      call tides%exit_data()
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
      call metrics%exit_data()

      call check(error, rms == rms, "RMS is NaN")

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
      call tides%destroy()
   end subroutine run_response

end module test_ocean_tides_sal
