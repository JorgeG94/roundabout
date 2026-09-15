!! Unit tests for the ocean horizontal tracer-diffusion kernel
!! (rdb_ocean_hdiff_tracer).  Flux-form Laplacian on `T = hTr/h`,
!! constant scalar `kappa_h`, closed-wall BC.
!!
!! Cases:
!!   * Constancy — uniform T (any h_layer) ⇒ no change.  Strongest
!!     discriminator for flux-form bugs.
!!   * Zero-coefficient short-circuit — `kappa_h = 0` is a no-op.
!!   * Conservation — total tracer mass `sum(hTr)` preserved
!!     exactly under closed walls.
!!   * Sinusoid decay rate — initialise a single Fourier mode in
!!     T (with uniform h), run forward-Euler.  Realised amplitude
!!     after N steps must match `(1 - dt*kappa*beta^2)^N` where
!!     `beta^2 = 2*(1-cos(k*dx))/dx^2`.
!!   * Per-tracer flag — set `do_horizontal_diffusion = .false.`
!!     on one tracer and verify it stays bit-for-bit at IC while
!!     the others diffuse.
module test_ocean_hdiff_tracer
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_hdiff_tracer, only: ocean_hdiff_tracer_t, &
                                     tracer_hdiff
   use rdb_ocean_budgets, only: ocean_budgets_t, BUDGET_HEAT_TOTAL, BUDGET_SALT_TOTAL
   use rdb_ocean_boundary_types, only: ocean_bc_state_t
   use rdb_config, only: config_t
   use rdb_ocean_state, only: ocean_state_t
   use rdb_ocean_setup, only: configure_ocean_hdiff
   implicit none
   private

   public :: collect_ocean_hdiff_tracer_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 3

contains

   subroutine collect_ocean_hdiff_tracer_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("hdiff_constancy", test_constancy), &
                  new_unittest("hdiff_zero_coeff_noop", test_zero_coeff), &
                  new_unittest("hdiff_conservation", test_conservation), &
                  new_unittest("hdiff_sinusoid_decay_rate", test_sinusoid_decay), &
                  new_unittest("hdiff_per_tracer_flag", test_per_tracer_flag), &
                  new_unittest("hdiff_budget_contributor_telescopes", &
                               test_hdiff_budget_contributor), &
                  new_unittest("hdiff_wall_is_physical_edge", test_wall_physical_edge), &
                  new_unittest("hdiff_kappa_from_config", test_kappa_from_config) &
                  ]
   end subroutine collect_ocean_hdiff_tracer_tests

   subroutine make_grid(grid, nx_phys, ny_phys, dx, dy)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      real(wp), intent(in) :: dx, dy
      call grid%init(nx_phys, ny_phys, NGHOST, dx, dy)
   end subroutine make_grid

   subroutine map_in(ms, hd, metrics, grid)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_hdiff_tracer_t), intent(inout) :: hd
      type(ocean_metrics_t), intent(inout) :: metrics
      type(hgrid_t), intent(in) :: grid
      call make_cartesian_metrics(metrics, grid)
      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(hd)
      call hd%enter_data()
   end subroutine map_in

   subroutine map_out(ms, hd, metrics)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_hdiff_tracer_t), intent(inout) :: hd
      type(ocean_metrics_t), intent(inout) :: metrics
      call hd%exit_data()
      !$acc exit data delete(hd)
      call ms%exit_data()
      !$acc exit data delete(ms)
      call destroy_cartesian_metrics(metrics)
   end subroutine map_out

   ! -----------------------------------------------------------------
   ! Cases
   ! -----------------------------------------------------------------

   subroutine test_constancy(error)
      !! Uniform T everywhere (any `h_layer`).  Face fluxes are
      !! all zero by construction → hTr unchanged.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_hdiff_tracer_t) :: hd
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: T0 = 20.0_wp
      real(wp), parameter :: DT = 1.0_wp
      real(wp) :: max_T_dev, T_obs
      integer :: i, j, k, nx, ny

      call make_grid(grid, 12, 10, 1.0_wp, 1.0_wp)
      ms%nz_ml = NZ
      call ms%init(grid)
      call hd%init(grid, nz_ml=NZ)
      hd%kappa_h = 1.0_wp
      nx = grid%nx_total
      ny = grid%ny_total

      ! Non-uniform h_layer, but uniform T.  Verifies h-weighted
      ! flux still gives zero net flux at constancy.
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx
               ms%h_layer(i, j, k) = 10.0_wp + 0.5_wp*real(i + j, wp)
            end do
         end do
      end do
      ms%tracers(ms%idx_temperature)%hTr = T0*ms%h_layer
      ms%tracers(ms%idx_salinity)%hTr = 35.0_wp*ms%h_layer

      call map_in(ms, hd, metrics, grid)
      call tracer_hdiff(grid, metrics, hd, ms, DT)
      call map_out(ms, hd, metrics)

      max_T_dev = 0.0_wp
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx
               T_obs = ms%tracers(ms%idx_temperature)%hTr(i, j, k)/ms%h_layer(i, j, k)
               max_T_dev = max(max_T_dev, abs(T_obs - T0))
            end do
         end do
      end do
      call check(error, max_T_dev < 1.0e-12_wp, &
                 "hdiff: uniform T not preserved")

      call hd%destroy(); call ms%destroy()
   end subroutine test_constancy

   subroutine test_zero_coeff(error)
      !! `kappa_h = 0` short-circuits the registry walk.  Tracer
      !! pattern stays bit-for-bit at IC.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_hdiff_tracer_t) :: hd
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp), parameter :: DT = 0.5_wp
      real(wp), allocatable :: hT_ic(:, :, :)
      integer :: i, j, k, nx, ny

      call make_grid(grid, 8, 6, 1.0_wp, 1.0_wp)
      ms%nz_ml = NZ
      call ms%init(grid)
      call hd%init(grid, nz_ml=NZ)
      hd%kappa_h = 0.0_wp
      nx = grid%nx_total
      ny = grid%ny_total

      ms%h_layer = 10.0_wp
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx
               ms%tracers(ms%idx_temperature)%hTr(i, j, k) = &
                  10.0_wp*sin(2.0_wp*PI*real(i, wp)/real(nx, wp))* &
                  ms%h_layer(i, j, k)
               ms%tracers(ms%idx_salinity)%hTr(i, j, k) = 35.0_wp*ms%h_layer(i, j, k)
            end do
         end do
      end do
      allocate (hT_ic, source=ms%tracers(ms%idx_temperature)%hTr)

      call map_in(ms, hd, metrics, grid)
      call tracer_hdiff(grid, metrics, hd, ms, DT)
      call map_out(ms, hd, metrics)

      call check(error, &
                 maxval(abs(ms%tracers(ms%idx_temperature)%hTr - hT_ic)) < 1.0e-12_wp, &
                 "hdiff: kappa_h=0 not a no-op")

      deallocate (hT_ic)
      call hd%destroy(); call ms%destroy()
   end subroutine test_zero_coeff

   subroutine test_conservation(error)
      !! Non-uniform T, several diffusion steps, closed-wall (all default,
      !! single-rank) BC.  STRENGTHENED per PR-9 §9.6: the ghost band is
      !! POISONED to a value far from the interior pattern, and the sum is
      !! taken over the PHYSICAL domain ONLY
      !! (nghost+1:nghost+nx_phys, nghost+1:nghost+ny_phys) — summing the
      !! FULL array (including ghosts, the pre-PR-9 form of this test) is
      !! conserved by construction under array-edge-only wall closure and
      !! is therefore structurally incapable of detecting the physical-
      !! edge bug this PR fixes.  Against the pre-fix kernel (zeroing only
      !! the array-bound faces i=1/nx+1, not the physical wall at
      !! i=nghost+1/nghost+nx_phys+1) this assertion FAILS — the physical
      !! wall face carried a real flux from the poisoned ghost into the
      !! interior.  Against the fixed kernel it holds to < 1e-12.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_hdiff_tracer_t) :: hd
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp), parameter :: DT = 0.05_wp
      real(wp), parameter :: POISON = 1.0e4_wp
      integer, parameter :: N_STEPS = 8
      real(wp) :: total_before, total_after, drift
      integer :: i, j, k, step, nx, ny, i0, i1, j0, j1

      call make_grid(grid, 10, 8, 1.0_wp, 1.0_wp)
      ms%nz_ml = NZ
      call ms%init(grid)
      call hd%init(grid, nz_ml=NZ)
      hd%kappa_h = 0.1_wp
      nx = grid%nx_total
      ny = grid%ny_total
      i0 = grid%nghost + 1
      i1 = grid%nghost + grid%nx_phys
      j0 = grid%nghost + 1
      j1 = grid%nghost + grid%ny_phys

      ms%h_layer = 10.0_wp
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx
               ms%tracers(ms%idx_temperature)%hTr(i, j, k) = &
                  (5.0_wp + 3.0_wp*sin(2.0_wp*PI*real(i, wp)/real(nx, wp))* &
                   cos(PI*real(j, wp)/real(ny, wp)))*ms%h_layer(i, j, k)
               ms%tracers(ms%idx_salinity)%hTr(i, j, k) = 35.0_wp*ms%h_layer(i, j, k)
            end do
         end do
      end do
      ! Poison the ghost band with a value far outside the interior
      ! pattern's range — any physical-wall leak becomes visible.
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx
               if (i < i0 .or. i > i1 .or. j < j0 .or. j > j1) then
                  ms%tracers(ms%idx_temperature)%hTr(i, j, k) = POISON*ms%h_layer(i, j, k)
               end if
            end do
         end do
      end do
      total_before = sum(ms%tracers(ms%idx_temperature)%hTr(i0:i1, j0:j1, :))

      call map_in(ms, hd, metrics, grid)
      do step = 1, N_STEPS
         call tracer_hdiff(grid, metrics, hd, ms, DT)
      end do
      call map_out(ms, hd, metrics)
      total_after = sum(ms%tracers(ms%idx_temperature)%hTr(i0:i1, j0:j1, :))

      drift = abs(total_after - total_before)/abs(total_before)
      call check(error, drift < 1.0e-12_wp, &
                 "hdiff: PHYSICAL-domain tracer mass drift > 1e-12 (poisoned ghost band)")

      call hd%destroy(); call ms%destroy()
   end subroutine test_conservation

   subroutine test_sinusoid_decay(error)
      !! T(i, j) = cos(k_x * x).  For the 5-point Laplacian,
      !! analytic decay per step is (1 - dt*kappa*beta^2) where
      !! beta^2 = 2*(1 - cos(k*dx))/dx^2.  Compare interior point
      !! to analytic.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_hdiff_tracer_t) :: hd
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp), parameter :: DT = 0.05_wp
      real(wp), parameter :: KAPPA = 0.05_wp
      integer, parameter :: N_STEPS = 10
      real(wp), parameter :: H0 = 10.0_wp
      integer :: i, j, k, step, nx, ny, i_probe, j_probe
      real(wp) :: k_wave, beta_sq, alpha
      real(wp) :: T_ic_at_probe, T_expected, T_obs, rel_err

      call make_grid(grid, 40, 12, 1.0_wp, 1.0_wp)
      ms%nz_ml = NZ
      call ms%init(grid)
      call hd%init(grid, nz_ml=NZ)
      hd%kappa_h = KAPPA
      nx = grid%nx_total
      ny = grid%ny_total

      k_wave = 2.0_wp*PI/real(nx, wp)
      beta_sq = 2.0_wp*(1.0_wp - cos(k_wave*grid%dx))/(grid%dx*grid%dx)
      alpha = 1.0_wp - DT*KAPPA*beta_sq

      ms%h_layer = H0
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx
               ms%tracers(ms%idx_temperature)%hTr(i, j, k) = &
                  cos(k_wave*real(i, wp)*grid%dx)*H0
               ms%tracers(ms%idx_salinity)%hTr(i, j, k) = 35.0_wp*H0
            end do
         end do
      end do

      call map_in(ms, hd, metrics, grid)
      do step = 1, N_STEPS
         call tracer_hdiff(grid, metrics, hd, ms, DT)
      end do
      call map_out(ms, hd, metrics)

      ! Probe near a trough (k_wave*x ≈ pi) far from x-walls so
      ! the wall effect — frozen boundary T — hasn't propagated
      ! to the probe over the diffusion time.
      i_probe = nx/2 + 1
      j_probe = ny/2
      T_ic_at_probe = cos(k_wave*real(i_probe, wp)*grid%dx)
      T_expected = T_ic_at_probe*alpha**N_STEPS
      T_obs = ms%tracers(ms%idx_temperature)%hTr(i_probe, j_probe, NZ/2)/H0
      rel_err = abs(T_obs - T_expected)/abs(T_expected)

      call check(error, rel_err < 5.0e-3_wp, &
                 "hdiff: sinusoid decay rate off analytic by > 5e-3")

      call hd%destroy(); call ms%destroy()
   end subroutine test_sinusoid_decay

   subroutine test_per_tracer_flag(error)
      !! Disable `do_horizontal_diffusion` on salinity, leave it
      !! enabled on temperature.  Set up a non-trivial pattern in
      !! BOTH.  After diffusion, salinity must be bit-for-bit at
      !! IC; temperature must have evolved.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_hdiff_tracer_t) :: hd
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp), parameter :: DT = 0.1_wp
      real(wp), allocatable :: hS_ic(:, :, :), hT_ic(:, :, :)
      real(wp) :: max_S_dev, max_T_dev
      integer :: i, j, k, nx, ny
      checks: block

         call make_grid(grid, 12, 10, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call hd%init(grid, nz_ml=NZ)
         hd%kappa_h = 1.0_wp
         ms%tracers(ms%idx_salinity)%do_horizontal_diffusion = .false.
         nx = grid%nx_total
         ny = grid%ny_total

         ms%h_layer = 10.0_wp
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx
                  ms%tracers(ms%idx_temperature)%hTr(i, j, k) = &
                     cos(2.0_wp*PI*real(i, wp)/real(nx, wp))*ms%h_layer(i, j, k)
                  ms%tracers(ms%idx_salinity)%hTr(i, j, k) = &
                     (35.0_wp + sin(2.0_wp*PI*real(j, wp)/real(ny, wp)))* &
                     ms%h_layer(i, j, k)
               end do
            end do
         end do
         allocate (hS_ic, source=ms%tracers(ms%idx_salinity)%hTr)
         allocate (hT_ic, source=ms%tracers(ms%idx_temperature)%hTr)

         call map_in(ms, hd, metrics, grid)
         call tracer_hdiff(grid, metrics, hd, ms, DT)
         call map_out(ms, hd, metrics)

         max_S_dev = maxval(abs(ms%tracers(ms%idx_salinity)%hTr - hS_ic))
         max_T_dev = maxval(abs(ms%tracers(ms%idx_temperature)%hTr - hT_ic))

         call check(error, max_S_dev < 1.0e-12_wp, &
                    "per-tracer flag: salinity touched despite do_horizontal_diffusion=.false.")
         if (allocated(error)) exit checks
         call check(error, max_T_dev > 1.0e-3_wp, &
                    "per-tracer flag: temperature should have diffused")

      end block checks
      deallocate (hS_ic, hT_ic)
      call hd%destroy(); call ms%destroy()
   end subroutine test_per_tracer_flag

   subroutine test_hdiff_budget_contributor(error)
      !! Phase D v2 hdiff kernel patch: heat + salt contributor slots
      !! telescope to ~0 over the spatial integral for closed-wall
      !! Laplacian diffusion (wall fluxes forced to zero).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_hdiff_tracer_t) :: hd
      type(ocean_metrics_t) :: metrics
      type(ocean_budgets_t) :: budgets
      real(wp), parameter :: H0 = 10.0_wp, DT = 60.0_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp) :: total_T_before, total_S_before, total_T_after, total_S_after
      real(wp) :: rel_T, rel_S
      integer :: i, j, k, nx, ny, idx_heat, idx_salt
      checks: block

         call make_grid(grid, 8, 6, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call hd%init(grid, nz_ml=NZ)
         hd%kappa_h = 1.0e-2_wp
         call budgets%init(grid)
         call budgets%register_contributor("hdiff_heat", BUDGET_HEAT_TOTAL, &
                                           ms%heat_budget_hdiff, &
                                           device_resident=.true.)
         idx_heat = budgets%n_contributors
         call budgets%register_contributor("hdiff_salt", BUDGET_SALT_TOTAL, &
                                           ms%salt_budget_hdiff, &
                                           device_resident=.true.)
         idx_salt = budgets%n_contributors

         nx = grid%nx_total
         ny = grid%ny_total
         ms%h_layer = H0
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx
                  ms%tracers(ms%idx_temperature)%hTr(i, j, k) = &
                     (15.0_wp + 0.5_wp*sin(2.0_wp*PI*real(i, wp)/real(nx, wp))* &
                      sin(PI*real(j, wp)/real(ny, wp)))*H0
                  ms%tracers(ms%idx_salinity)%hTr(i, j, k) = &
                     (35.0_wp + 0.3_wp*cos(2.0_wp*PI*real(i, wp)/real(nx, wp)))*H0
               end do
            end do
         end do

         total_T_before = sum(ms%tracers(ms%idx_temperature)%hTr)*grid%dx*grid%dy
         total_S_before = sum(ms%tracers(ms%idx_salinity)%hTr)*grid%dx*grid%dy

         call map_in(ms, hd, metrics, grid)
         call tracer_hdiff(grid, metrics, hd, ms, DT)
         call budgets%drain_contributors()
         call map_out(ms, hd, metrics)

         total_T_after = sum(ms%tracers(ms%idx_temperature)%hTr)*grid%dx*grid%dy
         total_S_after = sum(ms%tracers(ms%idx_salinity)%hTr)*grid%dx*grid%dy

         rel_T = abs(budgets%contributors(idx_heat)%total_integrated)/abs(total_T_before)
         rel_S = abs(budgets%contributors(idx_salt)%total_integrated)/abs(total_S_before)
         call check(error, rel_T < 1.0e-12_wp, &
                    "hdiff heat contributor should telescope to FP")
         if (allocated(error)) exit checks
         call check(error, rel_S < 1.0e-12_wp, &
                    "hdiff salt contributor should telescope to FP")
         if (allocated(error)) exit checks

         call check(error, abs((total_T_after - total_T_before) - &
                               budgets%contributors(idx_heat)%total_integrated) &
                    < 1.0e-12_wp*abs(total_T_before), &
                    "heat LHS = RHS for hdiff")
         if (allocated(error)) exit checks
         call check(error, abs((total_S_after - total_S_before) - &
                               budgets%contributors(idx_salt)%total_integrated) &
                    < 1.0e-12_wp*abs(total_S_before), &
                    "salt LHS = RHS for hdiff")
      end block checks
      call budgets%destroy()
      call hd%destroy(); call ms%destroy()
   end subroutine test_hdiff_budget_contributor

   subroutine test_wall_physical_edge(error)
      !! PR-9 §9.7.  All-wall `bc` (default-constructed: every edge
      !! OBC_WALL, every has_* = .true.), ghost band poisoned far from
      !! the interior value.  A closed PHYSICAL wall must admit NO flux
      !! from the ghosts: the interior stays unchanged to round-off —
      !! this is the assertion the pre-fix array-edge-only zeroing
      !! (i=1/nx+1 only, not i=nghost+1/nghost+nx_phys+1) would FAIL,
      !! since with NGHOST=2 the physical wall sits strictly inside the
      !! array bound.  Then `bc%has_west = .false.` (simulated MPI seam):
      !! the west column DOES respond — an MPI seam is a real interior
      !! face, mirroring `continuity_zonal_flux`'s has_*/OBC_WALL gate.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms_wall, ms_seam
      type(ocean_hdiff_tracer_t) :: hd_wall, hd_seam
      type(ocean_metrics_t) :: metrics_wall, metrics_seam
      type(ocean_bc_state_t) :: bc_wall, bc_seam
      real(wp), parameter :: T0 = 10.0_wp, POISON = 1.0e4_wp
      real(wp), parameter :: DT = 0.1_wp
      integer, parameter :: N_STEPS = 5
      integer :: i, j, k, step, nx, ny, i0, i1, j0, j1
      real(wp) :: max_interior_dev, west_dev

      checks: block
         call make_grid(grid, 10, 8, 1.0_wp, 1.0_wp)
         nx = grid%nx_total
         ny = grid%ny_total
         i0 = grid%nghost + 1
         i1 = grid%nghost + grid%nx_phys
         j0 = grid%nghost + 1
         j1 = grid%nghost + grid%ny_phys

         ! ---- Phase 1: all-wall, default bc -> interior untouched ----
         ms_wall%nz_ml = NZ
         call ms_wall%init(grid)
         call hd_wall%init(grid, nz_ml=NZ)
         hd_wall%kappa_h = 1.0_wp
         ms_wall%h_layer = 10.0_wp
         ms_wall%tracers(ms_wall%idx_temperature)%hTr = T0*ms_wall%h_layer
         ms_wall%tracers(ms_wall%idx_salinity)%hTr = 35.0_wp*ms_wall%h_layer
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx
                  if (i < i0 .or. i > i1 .or. j < j0 .or. j > j1) then
                     ms_wall%tracers(ms_wall%idx_temperature)%hTr(i, j, k) = &
                        POISON*ms_wall%h_layer(i, j, k)
                  end if
               end do
            end do
         end do
         ! bc_wall default-constructed: every edge OBC_WALL, every
         ! has_* = .true. (see rdb_ocean_boundary_types type defaults).

         call map_in(ms_wall, hd_wall, metrics_wall, grid)
         do step = 1, N_STEPS
            call tracer_hdiff(grid, metrics_wall, hd_wall, ms_wall, DT, bc=bc_wall)
         end do
         call map_out(ms_wall, hd_wall, metrics_wall)

         max_interior_dev = 0.0_wp
         do k = 1, NZ
            do j = j0, j1
               do i = i0, i1
                  max_interior_dev = max(max_interior_dev, &
                                         abs(ms_wall%tracers(ms_wall%idx_temperature)%hTr(i, j, k)/ &
                                             ms_wall%h_layer(i, j, k) - T0))
               end do
            end do
         end do
         call check(error, max_interior_dev < 1.0e-10_wp, &
                    "hdiff: closed physical wall must admit no flux from a poisoned ghost band")
         call hd_wall%destroy(); call ms_wall%destroy()
         if (allocated(error)) exit checks

         ! ---- Phase 2: bc%has_west = .false. (simulated MPI seam) ----
         ! -> the local wall-position face is a real interior face and
         ! MUST respond to the (poisoned) neighbour value there.
         ms_seam%nz_ml = NZ
         call ms_seam%init(grid)
         call hd_seam%init(grid, nz_ml=NZ)
         hd_seam%kappa_h = 1.0_wp
         ms_seam%h_layer = 10.0_wp
         ms_seam%tracers(ms_seam%idx_temperature)%hTr = T0*ms_seam%h_layer
         ms_seam%tracers(ms_seam%idx_salinity)%hTr = 35.0_wp*ms_seam%h_layer
         do k = 1, NZ
            do j = 1, ny
               do i = 1, nx
                  if (i < i0 .or. i > i1 .or. j < j0 .or. j > j1) then
                     ms_seam%tracers(ms_seam%idx_temperature)%hTr(i, j, k) = &
                        POISON*ms_seam%h_layer(i, j, k)
                  end if
               end do
            end do
         end do
         bc_seam%has_west = .false.

         call map_in(ms_seam, hd_seam, metrics_seam, grid)
         do step = 1, N_STEPS
            call tracer_hdiff(grid, metrics_seam, hd_seam, ms_seam, DT, bc=bc_seam)
         end do
         call map_out(ms_seam, hd_seam, metrics_seam)

         ! West-most physical column (i=i0) must have moved off T0 — the
         ! seam face carried the flux from the poisoned west ghost.
         west_dev = 0.0_wp
         do k = 1, NZ
            do j = j0, j1
               west_dev = max(west_dev, &
                              abs(ms_seam%tracers(ms_seam%idx_temperature)%hTr(i0, j, k)/ &
                                  ms_seam%h_layer(i0, j, k) - T0))
            end do
         end do
         call check(error, west_dev > 1.0e-6_wp, &
                    "hdiff: an MPI-seam (has_west=.false.) edge must NOT be hard-zeroed "// &
                    "— the west column must respond to the neighbour ghost value")
         call hd_seam%destroy(); call ms_seam%destroy()
      end block checks
   end subroutine test_wall_physical_edge

   subroutine test_kappa_from_config(error)
      !! PR-9 §9.8.  `&ocean_hdiff_nml kappa_h` reaches
      !! `ocean_state%hdiff_tracer%kappa_h` via `configure_ocean_hdiff` —
      !! driven from the CONFIG, not by setting the field directly (the
      !! gap the audit flagged: the existing suite only ever pokes
      !! `hd%kappa_h` directly, which is exactly why nobody noticed the
      !! namelist key did not exist).  Bit-identity half: the DEFAULT
      !! config leaves kappa_h = 0.0, so the kernel's short-circuit stays
      !! intact.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg, cfg_default
      type(ocean_state_t) :: ocean_state, ocean_state_default

      checks: block
         cfg%ocean%hdiff%kappa_h = 1.0e2_wp
         call configure_ocean_hdiff(cfg, ocean_state, compute_rank=0)
         call check(error, ocean_state%hdiff_tracer%kappa_h == 1.0e2_wp, &
                    "configure_ocean_hdiff must copy &ocean_hdiff_nml kappa_h onto "// &
                    "ocean_state%hdiff_tracer%kappa_h")
         if (allocated(error)) exit checks

         call configure_ocean_hdiff(cfg_default, ocean_state_default, compute_rank=0)
         call check(error, ocean_state_default%hdiff_tracer%kappa_h == 0.0_wp, &
                    "default config must leave kappa_h == 0.0 (bit-identical short-circuit)")
      end block checks
   end subroutine test_kappa_from_config

end module test_ocean_hdiff_tracer
