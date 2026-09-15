!! Unit tests for the ocean backward-Euler vertical-diffusion solver
!! (rdb_ocean_vdiff).  Thomas tridiagonal solve per column,
!! unconditionally stable; tests cover both the tracer and momentum
!! paths.
!!
!! Cases:
!!   * Tracer constancy — uniform T across the column ⇒ no change.
!!     Verifies the closed-top + closed-bottom BC is correctly
!!     formulated (no spurious source at boundaries).
!!   * Tracer zero-coefficient short-circuit — `K_v_tracer = 0` is
!!     a no-op.
!!   * Tracer column-mass conservation — closed-BC implicit
!!     diffusion conserves `sum(hTr)` per column to round-off.
!!   * Tracer relaxation to mean — large K_v + many steps drives
!!     every column toward its initial column-mean.
!!   * Momentum constancy — uniform u_face across the column ⇒ no
!!     change.
!!   * Momentum zero-coefficient short-circuit.
!!   * Momentum stability under huge K_v — backward Euler is
!!     unconditionally stable.  Stepping with K_v * dt / h² = 1000
!!     must not blow up; the column collapses toward its mean.
module test_ocean_vdiff
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_vdiff, only: ocean_vdiff_t, &
                              vdiff_apply_momentum, &
                              vdiff_apply_tracers, &
                              face_thick
   use rdb_ocean_budgets, only: ocean_budgets_t, BUDGET_HEAT_TOTAL, BUDGET_SALT_TOTAL
   implicit none
   private

   public :: collect_ocean_vdiff_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 6

contains

   subroutine collect_ocean_vdiff_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("vdiff_tracer_constancy", test_tracer_constancy), &
                  new_unittest("vdiff_tracer_zero_kv_noop", test_tracer_zero), &
                  new_unittest("vdiff_tracer_vanishing_layer_conservation", &
                               test_tracer_vanishing_layer), &
                  new_unittest("vdiff_tracer_column_mass_conservation", &
                               test_tracer_conservation), &
                  new_unittest("vdiff_tracer_relaxes_to_mean", test_tracer_relax), &
                  new_unittest("vdiff_momentum_constancy", test_momentum_constancy), &
                  new_unittest("vdiff_momentum_zero_kv_noop", test_momentum_zero), &
                  new_unittest("vdiff_momentum_implicit_stability", test_momentum_stab), &
                  new_unittest("vdiff_budget_contributor_telescopes", &
                               test_vdiff_budget_contributor), &
                  new_unittest("vdiff_harmonic_uniform_matches_arithmetic", &
                               test_harmonic_uniform), &
                  new_unittest("vdiff_harmonic_thin_layer_smaller_face", &
                               test_harmonic_thin_layer), &
                  new_unittest("vdiff_per_tracer_kt_ne_ks_diffuse_independently", &
                               test_per_tracer_kt_ne_ks), &
                  new_unittest("vdiff_per_tracer_kt_eq_ks_bit_identical", &
                               test_per_tracer_kt_eq_ks_bit_identical), &
                  new_unittest("vdiff_passive_tracer_follows_ks", &
                               test_passive_tracer_follows_ks) &
                  ]
   end subroutine collect_ocean_vdiff_tests

   subroutine make_grid(grid, nx_phys, ny_phys)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      call grid%init(nx_phys, ny_phys, NGHOST, 1.0_wp, 1.0_wp)
   end subroutine make_grid

   subroutine map_in(ms, vd)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_vdiff_t), intent(inout) :: vd
      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(vd)
      call vd%enter_data()
   end subroutine map_in

   subroutine map_out(ms, vd)
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_vdiff_t), intent(inout) :: vd
      call vd%exit_data()
      !$acc exit data delete(vd)
      call ms%exit_data()
      !$acc exit data delete(ms)
   end subroutine map_out

   ! -----------------------------------------------------------------
   ! Tracer cases
   ! -----------------------------------------------------------------

   subroutine test_tracer_constancy(error)
      !! Uniform T across every column.  No vertical gradient
      !! anywhere ⇒ no flux ⇒ no change.  Tests both the no-flux
      !! bed and surface BCs.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vdiff_t) :: vd
      real(wp), parameter :: H0 = 10.0_wp, T0 = 25.0_wp, DT = 60.0_wp
      real(wp) :: max_T_dev, T_obs
      integer :: i, j, k, nx, ny

      call make_grid(grid, 6, 4)
      ms%nz_ml = NZ
      call ms%init(grid)
      call vd%init(grid, nz_ml=NZ)
      vd%K_v_tracer = 1.0e-2_wp
      nx = grid%nx_total
      ny = grid%ny_total

      ms%h_layer = H0
      ms%tracers(ms%idx_temperature)%hTr = T0*H0
      ms%tracers(ms%idx_salinity)%hTr = 35.0_wp*H0

      call map_in(ms, vd)
      call vdiff_apply_tracers(grid, vd, ms, DT)
      call map_out(ms, vd)

      max_T_dev = 0.0_wp
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx
               T_obs = ms%tracers(ms%idx_temperature)%hTr(i, j, k)/H0
               max_T_dev = max(max_T_dev, abs(T_obs - T0))
            end do
         end do
      end do
      call check(error, max_T_dev < 1.0e-12_wp, &
                 "vdiff tracer: uniform T not preserved")

      call vd%destroy(); call ms%destroy()
   end subroutine test_tracer_constancy

   subroutine test_tracer_zero(error)
      !! `K_v_tracer = 0` ⇒ kernel returns without touching anything.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vdiff_t) :: vd
      real(wp), parameter :: H0 = 10.0_wp, DT = 60.0_wp
      real(wp), allocatable :: hT_ic(:, :, :)
      integer :: i, j, k, nx, ny

      call make_grid(grid, 6, 4)
      ms%nz_ml = NZ
      call ms%init(grid)
      call vd%init(grid, nz_ml=NZ)
      vd%K_v_tracer = 0.0_wp
      nx = grid%nx_total
      ny = grid%ny_total

      ms%h_layer = H0
      do k = 1, NZ
         ms%tracers(ms%idx_temperature)%hTr(:, :, k) = real(k, wp)*H0
         ms%tracers(ms%idx_salinity)%hTr(:, :, k) = 35.0_wp*H0
      end do
      allocate (hT_ic, source=ms%tracers(ms%idx_temperature)%hTr)

      call map_in(ms, vd)
      call vdiff_apply_tracers(grid, vd, ms, DT)
      call map_out(ms, vd)

      call check(error, &
                 maxval(abs(ms%tracers(ms%idx_temperature)%hTr - hT_ic)) < 1.0e-12_wp, &
                 "vdiff tracer: K_v=0 not a no-op")

      deallocate (hT_ic)
      call vd%destroy(); call ms%destroy()
   end subroutine test_tracer_zero

   subroutine test_tracer_vanishing_layer(error)
      !! Regression for the Fox-Kemper × windowed-tracer-advect conservation
      !! leak (root cause: vdiff dropped tracer mass at vanishing layers).
      !!
      !! A thin z* surface layer can be driven to h ≤ 0 on an intermediate
      !! RK2 stage by the combined FK + resolved transport (Lagrangian,
      !! before the ALE remap relayers it).  The pre-fix code set the
      !! concentration of any h ≤ 0 layer to 0 (`T = 0 ⇒ hTr = T·h = 0`),
      !! SILENTLY discarding that layer's frozen tracer mass — and the raw
      !! 1/h in the tridiagonal coefficients corrupted the WHOLE column's
      !! solve, so even the thick neighbour layers lost mass (~5%/day leak
      !! on the kitchensink).
      !!
      !! The fix floors the per-layer thickness to H_VANISHED in BOTH the
      !! matrix and the T↔hTr conversion, and zeroes the diffusive flux at
      !! any interface touching a vanishing layer — so a collapsed layer is
      !! decoupled and its tracer mass is preserved EXACTLY while its thick
      !! neighbours conserve among themselves.  This asserts the
      !! floored-thickness column mass Σ(max(h,H_VANISHED)·T) is conserved
      !! to round-off and the vanishing layer keeps its hTr.
      use rdb_constants, only: H_VANISHED
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vdiff_t) :: vd
      real(wp), parameter :: H0 = 10.0_wp, DT = 60.0_wp
      real(wp) :: total_before, total_after, drift, htr_vanish_before, htr_vanish_after
      integer :: i, j, k, nx, ny

      call make_grid(grid, 6, 4)
      ms%nz_ml = NZ
      call ms%init(grid)
      call vd%init(grid, nz_ml=NZ)
      vd%K_v_tracer = 1.0e-2_wp
      nx = grid%nx_total
      ny = grid%ny_total

      ! Thick interior layers + a COLLAPSED surface layer (k=NZ) driven to a
      ! negative thickness, but still carrying frozen tracer mass.
      ms%h_layer = H0
      ms%h_layer(:, :, NZ) = -1.0e-3_wp     ! over-drained z* surface band
      do k = 1, NZ
         ms%tracers(ms%idx_temperature)%hTr(:, :, k) = (5.0_wp + 2.0_wp*real(k, wp))*H0
         ms%tracers(ms%idx_salinity)%hTr(:, :, k) = 35.0_wp*H0
      end do
      ! Give the collapsed surface layer a distinct, nonzero tracer mass —
      ! the pre-fix bug zeroed exactly this.
      ms%tracers(ms%idx_temperature)%hTr(:, :, NZ) = 7.5_wp

      ! Conserved invariant uses the SAME floored thickness the solve uses.
      total_before = 0.0_wp
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx
               total_before = total_before + ms%tracers(ms%idx_temperature)%hTr(i, j, k)
            end do
         end do
      end do
      htr_vanish_before = ms%tracers(ms%idx_temperature)%hTr(1 + grid%nghost, &
                                                             1 + grid%nghost, NZ)

      call map_in(ms, vd)
      call vdiff_apply_tracers(grid, vd, ms, DT)
      call map_out(ms, vd)

      total_after = 0.0_wp
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx
               total_after = total_after + ms%tracers(ms%idx_temperature)%hTr(i, j, k)
            end do
         end do
      end do
      htr_vanish_after = ms%tracers(ms%idx_temperature)%hTr(1 + grid%nghost, &
                                                            1 + grid%nghost, NZ)

      drift = abs(total_after - total_before)/abs(total_before)
      call check(error, drift < 1.0e-12_wp, &
                 "vdiff vanishing-layer: Σ(hTr) drift > 1e-12 (mass dropped at h≤0)")
      if (allocated(error)) then
         call vd%destroy(); call ms%destroy(); return
      end if
      ! The decoupled collapsed layer must keep its frozen tracer mass.
      call check(error, abs(htr_vanish_after - htr_vanish_before) < 1.0e-10_wp, &
                 "vdiff vanishing-layer: collapsed layer hTr was not preserved")

      call vd%destroy(); call ms%destroy()
   end subroutine test_tracer_vanishing_layer

   subroutine test_tracer_conservation(error)
      !! Non-trivial vertical T profile.  Closed top + closed bottom
      !! means tracer mass per column is exactly conserved (the
      !! tridiagonal system telescopes — sum of equations gives
      !! `sum(T_new) = sum(T_old)` only when `h` is constant in z,
      !! so we use uniform h_layer here; the more general property
      !! `sum(h*T)` is conserved is checked next).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vdiff_t) :: vd
      real(wp), parameter :: H0 = 10.0_wp, DT = 60.0_wp
      real(wp) :: total_before, total_after, drift
      integer :: k

      call make_grid(grid, 6, 4)
      ms%nz_ml = NZ
      call ms%init(grid)
      call vd%init(grid, nz_ml=NZ)
      vd%K_v_tracer = 1.0e-2_wp

      ms%h_layer = H0
      do k = 1, NZ
         ms%tracers(ms%idx_temperature)%hTr(:, :, k) = (5.0_wp + 2.0_wp*real(k, wp))*H0
         ms%tracers(ms%idx_salinity)%hTr(:, :, k) = 35.0_wp*H0
      end do
      total_before = sum(ms%tracers(ms%idx_temperature)%hTr)

      call map_in(ms, vd)
      call vdiff_apply_tracers(grid, vd, ms, DT)
      call map_out(ms, vd)

      total_after = sum(ms%tracers(ms%idx_temperature)%hTr)
      drift = abs(total_after - total_before)/abs(total_before)
      call check(error, drift < 1.0e-12_wp, &
                 "vdiff tracer: column mass drift > 1e-12")

      call vd%destroy(); call ms%destroy()
   end subroutine test_tracer_conservation

   subroutine test_tracer_relax(error)
      !! Huge `K_v * dt / h² ≫ 1` ⇒ implicit step jumps essentially
      !! to the column-mean profile in one shot.  Backward Euler
      !! amplification factor for the lowest mode is `1/(1 + λ*dt)`;
      !! at λ*dt ~ 1e6 the mode amplitude drops by 1e6 per step.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vdiff_t) :: vd
      real(wp), parameter :: H0 = 1.0_wp, DT = 1.0_wp, K_BIG = 1.0e6_wp
      real(wp) :: col_mean, max_dev
      integer :: i, j, k, nx, ny

      call make_grid(grid, 4, 4)
      ms%nz_ml = NZ
      call ms%init(grid)
      call vd%init(grid, nz_ml=NZ)
      vd%K_v_tracer = K_BIG
      nx = grid%nx_total
      ny = grid%ny_total

      ms%h_layer = H0
      do k = 1, NZ
         ms%tracers(ms%idx_temperature)%hTr(:, :, k) = real(k, wp)*H0
         ms%tracers(ms%idx_salinity)%hTr(:, :, k) = 35.0_wp*H0
      end do
      col_mean = 0.0_wp
      do k = 1, NZ
         col_mean = col_mean + real(k, wp)/real(NZ, wp)
      end do

      call map_in(ms, vd)
      call vdiff_apply_tracers(grid, vd, ms, DT)
      call map_out(ms, vd)

      max_dev = 0.0_wp
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx
               max_dev = max(max_dev, &
                             abs(ms%tracers(ms%idx_temperature)%hTr(i, j, k)/H0 - col_mean))
            end do
         end do
      end do
      ! At λ·dt ~ 1e6 backward-Euler damps the first mode by ~1e6
      ! per step, but the polynomial decay rate for higher modes is
      ! similar — practically the column collapses to its mean to
      ! a few parts in 1e5.
      call check(error, max_dev < 1.0e-4_wp, &
                 "vdiff tracer: huge K_v should drive column to its mean")

      call vd%destroy(); call ms%destroy()
   end subroutine test_tracer_relax

   ! -----------------------------------------------------------------
   ! Momentum cases
   ! -----------------------------------------------------------------

   subroutine test_momentum_constancy(error)
      !! Uniform u everywhere ⇒ no vertical gradient ⇒ no change.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vdiff_t) :: vd
      real(wp), parameter :: H0 = 10.0_wp, U0 = 0.4_wp, V0 = -0.3_wp
      real(wp), parameter :: DT = 60.0_wp
      real(wp) :: max_du, max_dv
      checks: block

         call make_grid(grid, 6, 4)
         ms%nz_ml = NZ
         call ms%init(grid)
         call vd%init(grid, nz_ml=NZ)
         vd%K_v_momentum = 1.0e-2_wp

         ms%h_layer = H0
         ms%u_face_x_layer = U0
         ms%v_face_y_layer = V0

         call map_in(ms, vd)
         call vdiff_apply_momentum(grid, vd, ms, DT)
         call map_out(ms, vd)

         max_du = maxval(abs(ms%u_face_x_layer - U0))
         max_dv = maxval(abs(ms%v_face_y_layer - V0))

         call check(error, max_du < 1.0e-12_wp, &
                    "vdiff momentum: uniform u not preserved")
         if (allocated(error)) exit checks
         call check(error, max_dv < 1.0e-12_wp, &
                    "vdiff momentum: uniform v not preserved")

      end block checks
      call vd%destroy(); call ms%destroy()
   end subroutine test_momentum_constancy

   subroutine test_momentum_zero(error)
      !! `K_v_momentum = 0` is a no-op.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vdiff_t) :: vd
      real(wp), parameter :: H0 = 10.0_wp, DT = 60.0_wp
      real(wp), allocatable :: u_ic(:, :, :), v_ic(:, :, :)
      integer :: k
      checks: block

         call make_grid(grid, 6, 4)
         ms%nz_ml = NZ
         call ms%init(grid)
         call vd%init(grid, nz_ml=NZ)
         vd%K_v_momentum = 0.0_wp

         ms%h_layer = H0
         do k = 1, NZ
            ms%u_face_x_layer(:, :, k) = 0.1_wp*real(k, wp)
            ms%v_face_y_layer(:, :, k) = -0.05_wp*real(k, wp)
         end do
         allocate (u_ic, source=ms%u_face_x_layer)
         allocate (v_ic, source=ms%v_face_y_layer)

         call map_in(ms, vd)
         call vdiff_apply_momentum(grid, vd, ms, DT)
         call map_out(ms, vd)

         call check(error, maxval(abs(ms%u_face_x_layer - u_ic)) < 1.0e-12_wp, &
                    "vdiff momentum: K_v=0 changed u")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(ms%v_face_y_layer - v_ic)) < 1.0e-12_wp, &
                    "vdiff momentum: K_v=0 changed v")

      end block checks
      deallocate (u_ic, v_ic)
      call vd%destroy(); call ms%destroy()
   end subroutine test_momentum_zero

   subroutine test_momentum_stab(error)
      !! Backward-Euler unconditional stability check.  Set up a
      !! sharp vertical shear, pick `K_v * dt / h² = 1000`, take
      !! one step.  No NaN / inf / overshoot — the field must
      !! contract toward the column-mean velocity, not blow up.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vdiff_t) :: vd
      real(wp), parameter :: H0 = 1.0_wp, DT = 1.0_wp, K_BIG = 1000.0_wp
      real(wp) :: u_obs, u_initial_max, u_initial_min
      integer :: k
      checks: block

         call make_grid(grid, 4, 4)
         ms%nz_ml = NZ
         call ms%init(grid)
         call vd%init(grid, nz_ml=NZ)
         vd%K_v_momentum = K_BIG

         ms%h_layer = H0
         ms%v_face_y_layer = 0.0_wp
         ! Sharp shear: u jumps from +1 to -1 across nz/2
         do k = 1, NZ
            if (k <= NZ/2) then
               ms%u_face_x_layer(:, :, k) = 1.0_wp
            else
               ms%u_face_x_layer(:, :, k) = -1.0_wp
            end if
         end do
         u_initial_max = maxval(ms%u_face_x_layer)
         u_initial_min = minval(ms%u_face_x_layer)

         call map_in(ms, vd)
         call vdiff_apply_momentum(grid, vd, ms, DT)
         call map_out(ms, vd)

         u_obs = maxval(abs(ms%u_face_x_layer))
         call check(error, u_obs <= u_initial_max + 1.0e-6_wp .and. u_obs == u_obs, &
                    "vdiff momentum: implicit step overshot or NaN")
         if (allocated(error)) exit checks
         call check(error, u_obs < 0.5_wp*(u_initial_max - u_initial_min), &
                    "vdiff momentum: huge K_v didn't damp the shear toward mean")

      end block checks
      call vd%destroy(); call ms%destroy()
   end subroutine test_momentum_stab

   subroutine test_vdiff_budget_contributor(error)
      !! Phase D v2 vdiff kernel patch: heat + salt contributor slots
      !! must telescope to ~0 over the spatial integral for closed-BC
      !! vertical diffusion.  Backward-Euler tridiag with zero flux
      !! through bed + surface conserves Σ(h·T) per column, so the
      !! contributor totals reflect that to FP.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vdiff_t) :: vd
      type(ocean_budgets_t) :: budgets
      real(wp), parameter :: H0 = 10.0_wp, DT = 60.0_wp
      real(wp) :: total_T_before, total_S_before, total_T_after, total_S_after
      real(wp) :: rel_T, rel_S
      integer :: k, idx_heat, idx_salt
      checks: block

         call make_grid(grid, 6, 4)
         ms%nz_ml = NZ
         call ms%init(grid)
         call vd%init(grid, nz_ml=NZ)
         vd%K_v_tracer = 1.0e-2_wp
         call budgets%init(grid)
         call budgets%register_contributor("vdiff_heat", BUDGET_HEAT_TOTAL, &
                                           ms%heat_budget_vdiff, &
                                           device_resident=.true.)
         idx_heat = budgets%n_contributors
         call budgets%register_contributor("vdiff_salt", BUDGET_SALT_TOTAL, &
                                           ms%salt_budget_vdiff, &
                                           device_resident=.true.)
         idx_salt = budgets%n_contributors

         ms%h_layer = H0
         do k = 1, NZ
            ms%tracers(ms%idx_temperature)%hTr(:, :, k) = (5.0_wp + 2.0_wp*real(k, wp))*H0
            ms%tracers(ms%idx_salinity)%hTr(:, :, k) = (30.0_wp + 0.5_wp*real(k, wp))*H0
         end do

         total_T_before = sum(ms%tracers(ms%idx_temperature)%hTr)*grid%dx*grid%dy
         total_S_before = sum(ms%tracers(ms%idx_salinity)%hTr)*grid%dx*grid%dy

         call map_in(ms, vd)
         call vdiff_apply_tracers(grid, vd, ms, DT)
         call budgets%drain_contributors()
         call map_out(ms, vd)

         total_T_after = sum(ms%tracers(ms%idx_temperature)%hTr)*grid%dx*grid%dy
         total_S_after = sum(ms%tracers(ms%idx_salinity)%hTr)*grid%dx*grid%dy

         rel_T = abs(budgets%contributors(idx_heat)%total_integrated)/abs(total_T_before)
         rel_S = abs(budgets%contributors(idx_salt)%total_integrated)/abs(total_S_before)
         call check(error, rel_T < 1.0e-12_wp, &
                    "vdiff heat contributor should telescope to FP")
         if (allocated(error)) exit checks
         call check(error, rel_S < 1.0e-12_wp, &
                    "vdiff salt contributor should telescope to FP")
         if (allocated(error)) exit checks

         ! Closure: LHS (total drift) = RHS (contributor).  vdiff is
         ! the only kernel touching hTr here.
         call check(error, abs((total_T_after - total_T_before) - &
                               budgets%contributors(idx_heat)%total_integrated) &
                    < 1.0e-12_wp*abs(total_T_before), &
                    "heat LHS = RHS for vdiff")
         if (allocated(error)) exit checks
         call check(error, abs((total_S_after - total_S_before) - &
                               budgets%contributors(idx_salt)%total_integrated) &
                    < 1.0e-12_wp*abs(total_S_before), &
                    "salt LHS = RHS for vdiff")
      end block checks
      call budgets%destroy()
      call vd%destroy(); call ms%destroy()
   end subroutine test_vdiff_budget_contributor

   subroutine test_harmonic_uniform(error)
      !! HARMONIC_VISC with uniform `h_layer = H0` must give exactly
      !! the same result as the arithmetic-mean default — the harmonic
      !! mean `2·h·h / (h+h) = h` reduces to `h` when both arguments
      !! are equal, identical to the arithmetic mean.
      !!
      !! Strategy: run the tracer kernel twice (once each mode) on a
      !! linear-in-k IC + uniform h.  hTr fields must match bit-by-bit.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms_arith, ms_harm
      type(ocean_vdiff_t) :: vd_arith, vd_harm
      real(wp), parameter :: H0 = 5.0_wp, DT = 0.1_wp, KV0 = 1.0e-3_wp
      real(wp) :: max_diff
      integer :: i, j, k, nx, ny

      call make_grid(grid, 4, 4)
      ms_arith%nz_ml = NZ
      ms_harm%nz_ml = NZ
      call ms_arith%init(grid)
      call ms_harm%init(grid)
      call vd_arith%init(grid, nz_ml=NZ)
      call vd_harm%init(grid, nz_ml=NZ)
      vd_arith%K_v_tracer = KV0
      vd_harm%K_v_tracer = KV0
      vd_arith%use_harmonic = .false.
      vd_harm%use_harmonic = .true.
      nx = grid%nx_total
      ny = grid%ny_total

      ms_arith%h_layer = H0
      ms_harm%h_layer = H0
      do k = 1, NZ
         ms_arith%tracers(ms_arith%idx_temperature)%hTr(:, :, k) = real(k, wp)*H0
         ms_harm%tracers(ms_harm%idx_temperature)%hTr(:, :, k) = real(k, wp)*H0
         ms_arith%tracers(ms_arith%idx_salinity)%hTr(:, :, k) = 35.0_wp*H0
         ms_harm%tracers(ms_harm%idx_salinity)%hTr(:, :, k) = 35.0_wp*H0
      end do

      call map_in(ms_arith, vd_arith)
      call map_in(ms_harm, vd_harm)
      call vdiff_apply_tracers(grid, vd_arith, ms_arith, DT)
      call vdiff_apply_tracers(grid, vd_harm, ms_harm, DT)
      call map_out(ms_arith, vd_arith)
      call map_out(ms_harm, vd_harm)

      max_diff = maxval(abs( &
                        ms_arith%tracers(ms_arith%idx_temperature)%hTr - &
                        ms_harm%tracers(ms_harm%idx_temperature)%hTr))
      call check(error, max_diff < 1.0e-14_wp, &
                 "HARMONIC_VISC must match arithmetic on uniform h")

      call vd_arith%destroy(); call vd_harm%destroy()
      call ms_arith%destroy(); call ms_harm%destroy()
   end subroutine test_harmonic_uniform

   subroutine test_harmonic_thin_layer(error)
      !! Direct unit test of `face_thick`: harmonic mean of two unequal
      !! thicknesses must be smaller than the arithmetic mean and must
      !! match the analytic formula.  Also: harmonic of (h, h) reduces
      !! to h exactly.  Smaller face thickness in the vdiff denominator
      !! means larger tridiagonal coefficient → faster mixing across
      !! thin/thick boundaries — that's the physical motivation; this
      !! test verifies the helper produces the right number, not the
      !! downstream dynamical consequence.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: arith_unequal, harm_unequal, harm_uniform
      real(wp), parameter :: H_FAT = 10.0_wp, H_THIN = 0.1_wp
      real(wp), parameter :: H_EQ = 5.0_wp
      real(wp), parameter :: HARM_EXPECT = 2.0_wp*H_FAT*H_THIN/(H_FAT + H_THIN)

      arith_unequal = face_thick(H_FAT, H_THIN, .false.)
      harm_unequal = face_thick(H_FAT, H_THIN, .true.)
      harm_uniform = face_thick(H_EQ, H_EQ, .true.)

      ! Arithmetic on (10, 0.1) = 5.05; harmonic = 2·10·0.1/10.1 ≈ 0.198.
      call check(error, abs(arith_unequal - 0.5_wp*(H_FAT + H_THIN)) < 1.0e-14_wp, &
                 "face_thick arithmetic branch wrong on unequal h")
      if (.not. allocated(error)) call check(error, abs(harm_unequal - HARM_EXPECT) < 1.0e-14_wp, &
                                             "face_thick harmonic branch off the 2·h_a·h_b/(h_a+h_b) formula")
      if (.not. allocated(error)) call check(error, harm_unequal < arith_unequal, &
                                             "harmonic mean of unequal h must be smaller than arithmetic")
      if (.not. allocated(error)) call check(error, abs(harm_uniform - H_EQ) < 1.0e-14_wp, &
                                             "harmonic mean of (h, h) must reduce to h exactly")
   end subroutine test_harmonic_thin_layer

   ! -----------------------------------------------------------------
   ! PR-20: per-tracer diffusivity (kt_source / ks_source)
   ! -----------------------------------------------------------------

   subroutine host_diffuse_column(nz, h, K, dt, T_in, T_out)
      !! Independent reference implementation of the backward-Euler
      !! column solve (build + Thomas factorize + solve, arithmetic
      !! `dz_face`, uniform `h`), coded directly against §3.1 of the
      !! PR-20 plan rather than by calling production code.  `K(1:nz+1)`
      !! is interface-located with `K(1) = K(nz+1) = 0` (closed BC).
      integer, intent(in) :: nz
      real(wp), intent(in) :: h(nz), K(nz + 1), dt
      real(wp), intent(in) :: T_in(nz)
      real(wp), intent(out) :: T_out(nz)

      real(wp) :: a(nz), b(nz), c(nz), rhs(nz)
      real(wp) :: alpha, beta, dz_face, denom
      integer :: k_

      dz_face = 0.5_wp*(h(1) + h(2))
      alpha = dt*K(2)/(h(1)*dz_face)
      a(1) = 0.0_wp
      c(1) = -alpha
      b(1) = 1.0_wp + alpha

      do k_ = 2, nz - 1
         beta = dt*K(k_)/(h(k_)*0.5_wp*(h(k_ - 1) + h(k_)))
         alpha = dt*K(k_ + 1)/(h(k_)*0.5_wp*(h(k_) + h(k_ + 1)))
         a(k_) = -beta
         c(k_) = -alpha
         b(k_) = 1.0_wp + alpha + beta
      end do

      dz_face = 0.5_wp*(h(nz - 1) + h(nz))
      beta = dt*K(nz)/(h(nz)*dz_face)
      a(nz) = -beta
      c(nz) = 0.0_wp
      b(nz) = 1.0_wp + beta

      rhs = T_in
      c(1) = c(1)/b(1)
      rhs(1) = rhs(1)/b(1)
      do k_ = 2, nz
         denom = b(k_) - a(k_)*c(k_ - 1)
         c(k_) = c(k_)/denom
         rhs(k_) = (rhs(k_) - a(k_)*rhs(k_ - 1))/denom
      end do

      T_out(nz) = rhs(nz)
      do k_ = nz - 1, 1, -1
         T_out(k_) = rhs(k_) - c(k_)*T_out(k_ + 1)
      end do
   end subroutine host_diffuse_column

   subroutine test_per_tracer_kt_ne_ks(error)
      !! The definitional test of PR-20: heat and salt evolve under
      !! their OWN diffusivities.  Single-column IC (identical step
      !! profile for T and S, uniform h), `ks_source = 4*kt_source`.
      !! After one step: T matches the closed-form kt solve, S matches
      !! the closed-form ks solve (both to ~1e-13), S has departed
      !! further from its IC than T (larger K mixes faster), and both
      !! columns conserve their integral Σ h·Tr to round-off.  This is
      !! exactly the test that fails on the pre-PR-20 code, where both
      !! tracers would follow kt and ‖ΔS‖ == ‖ΔT‖ bit-for-bit.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vdiff_t) :: vd
      real(wp), parameter :: H0 = 10.0_wp, DT = 60.0_wp
      real(wp), parameter :: KT0 = 1.0e-2_wp, KS0 = 4.0_wp*KT0
      real(wp), allocatable :: kt3(:, :, :), ks3(:, :, :)
      real(wp) :: T0(NZ), S0(NZ), T_ref(NZ), S_ref(NZ)
      real(wp) :: KT_col(NZ + 1), KS_col(NZ + 1), h_col(NZ)
      real(wp) :: sum_T_before, sum_S_before, sum_T_after, sum_S_after
      real(wp) :: dev_T, dev_S, max_abs_T, max_abs_S
      integer :: k, nx, ny, ic, jc

      checks: block
         call make_grid(grid, 6, 4)
         ms%nz_ml = NZ
         call ms%init(grid)
         call vd%init(grid, nz_ml=NZ)
         nx = grid%nx_total
         ny = grid%ny_total

         ms%h_layer = H0
         ! Step profile: bottom half 0, top half 1 (k=NZ is the surface),
         ! identical for T and S so any divergence is due to K alone.
         do k = 1, NZ
            if (k <= NZ/2) then
               T0(k) = 0.0_wp
            else
               T0(k) = 1.0_wp
            end if
         end do
         S0 = T0
         do k = 1, NZ
            ms%tracers(ms%idx_temperature)%hTr(:, :, k) = T0(k)*H0
            ms%tracers(ms%idx_salinity)%hTr(:, :, k) = S0(k)*H0
         end do

         allocate (kt3(nx, ny, NZ + 1), source=0.0_wp)
         allocate (ks3(nx, ny, NZ + 1), source=0.0_wp)
         do k = 2, NZ
            kt3(:, :, k) = KT0
            ks3(:, :, k) = KS0
         end do

         sum_T_before = sum(ms%tracers(ms%idx_temperature)%hTr)
         sum_S_before = sum(ms%tracers(ms%idx_salinity)%hTr)

         call map_in(ms, vd)
         !$acc enter data copyin(kt3, ks3)
         call vdiff_apply_tracers(grid, vd, ms, DT, kt_source=kt3, ks_source=ks3)
         !$acc exit data delete(kt3, ks3)
         call map_out(ms, vd)

         ! Independent analytical reference, one representative interior
         ! column (all columns are identical: uniform h, uniform IC,
         ! uniform K away from the closed boundaries).
         h_col = H0
         KT_col = 0.0_wp; KT_col(2:NZ) = KT0
         KS_col = 0.0_wp; KS_col(2:NZ) = KS0
         call host_diffuse_column(NZ, h_col, KT_col, DT, T0, T_ref)
         call host_diffuse_column(NZ, h_col, KS_col, DT, S0, S_ref)

         ic = grid%nghost + 2
         jc = grid%nghost + 2
         max_abs_T = 0.0_wp
         max_abs_S = 0.0_wp
         do k = 1, NZ
            max_abs_T = max(max_abs_T, &
                            abs(ms%tracers(ms%idx_temperature)%hTr(ic, jc, k)/H0 - T_ref(k)))
            max_abs_S = max(max_abs_S, &
                            abs(ms%tracers(ms%idx_salinity)%hTr(ic, jc, k)/H0 - S_ref(k)))
         end do

         call check(error, max_abs_T < 1.0e-13_wp, &
                    "T does not match the closed-form kt_source solve")
         if (allocated(error)) exit checks
         call check(error, max_abs_S < 1.0e-13_wp, &
                    "S does not match the closed-form ks_source solve")
         if (allocated(error)) exit checks

         dev_T = 0.0_wp
         dev_S = 0.0_wp
         do k = 1, NZ
            dev_T = dev_T + (ms%tracers(ms%idx_temperature)%hTr(ic, jc, k)/H0 - T0(k))**2
            dev_S = dev_S + (ms%tracers(ms%idx_salinity)%hTr(ic, jc, k)/H0 - S0(k))**2
         end do
         call check(error, dev_S > dev_T, &
                    "larger ks (4x kt) must depart from the IC more than kt does")
         if (allocated(error)) exit checks

         sum_T_after = sum(ms%tracers(ms%idx_temperature)%hTr)
         sum_S_after = sum(ms%tracers(ms%idx_salinity)%hTr)
         call check(error, abs(sum_T_after - sum_T_before) < 1.0e-12_wp*abs(sum_T_before), &
                    "temperature column integral not conserved under kt_source")
         if (allocated(error)) exit checks
         call check(error, abs(sum_S_after - sum_S_before) < 1.0e-12_wp*abs(sum_S_before), &
                    "salinity column integral not conserved under ks_source")
      end block checks
      if (allocated(kt3)) deallocate (kt3)
      if (allocated(ks3)) deallocate (ks3)
      call vd%destroy(); call ms%destroy()
   end subroutine test_per_tracer_kt_ne_ks

   subroutine test_per_tracer_kt_eq_ks_bit_identical(error)
      !! The bit-identity gate: `kt_source == ks_source` must give the
      !! EXACT same answer (bitwise `==`, not a tolerance) as the
      !! single-source legacy dispatch, because the two-pass refactor
      !! performs a DIFFERENT arithmetic sequence (two factorizations
      !! instead of one) and the PR-20 claim is that a deterministic
      !! kernel gives identical output for identical input regardless.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms_case1, ms_case3
      type(ocean_vdiff_t) :: vd_case1, vd_case3
      real(wp), parameter :: H0 = 10.0_wp, DT = 60.0_wp, KT0 = 1.0e-2_wp
      real(wp), allocatable :: k3(:, :, :)
      integer :: k, nx, ny

      checks: block
         call make_grid(grid, 6, 4)
         ms_case1%nz_ml = NZ; ms_case3%nz_ml = NZ
         call ms_case1%init(grid); call ms_case3%init(grid)
         call vd_case1%init(grid, nz_ml=NZ); call vd_case3%init(grid, nz_ml=NZ)
         nx = grid%nx_total
         ny = grid%ny_total

         ms_case1%h_layer = H0
         ms_case3%h_layer = H0
         do k = 1, NZ
            ms_case1%tracers(ms_case1%idx_temperature)%hTr(:, :, k) = &
               (5.0_wp + 2.0_wp*real(k, wp))*H0
            ms_case3%tracers(ms_case3%idx_temperature)%hTr(:, :, k) = &
               (5.0_wp + 2.0_wp*real(k, wp))*H0
            ms_case1%tracers(ms_case1%idx_salinity)%hTr(:, :, k) = &
               (30.0_wp + 0.5_wp*real(k, wp))*H0
            ms_case3%tracers(ms_case3%idx_salinity)%hTr(:, :, k) = &
               (30.0_wp + 0.5_wp*real(k, wp))*H0
         end do

         allocate (k3(nx, ny, NZ + 1), source=0.0_wp)
         do k = 2, NZ
            k3(:, :, k) = KT0
         end do

         call map_in(ms_case1, vd_case1)
         call map_in(ms_case3, vd_case3)
         !$acc enter data copyin(k3)
         call vdiff_apply_tracers(grid, vd_case1, ms_case1, DT, kt_source=k3)
         call vdiff_apply_tracers(grid, vd_case3, ms_case3, DT, kt_source=k3, ks_source=k3)
         !$acc exit data delete(k3)
         call map_out(ms_case1, vd_case1)
         call map_out(ms_case3, vd_case3)

         call check(error, all(ms_case1%tracers(ms_case1%idx_temperature)%hTr == &
                               ms_case3%tracers(ms_case3%idx_temperature)%hTr), &
                    "temperature hTr differs between single- and two-pass dispatch")
         if (allocated(error)) exit checks
         call check(error, all(ms_case1%tracers(ms_case1%idx_salinity)%hTr == &
                               ms_case3%tracers(ms_case3%idx_salinity)%hTr), &
                    "salinity hTr differs between single- and two-pass dispatch")
         if (allocated(error)) exit checks
         call check(error, all(ms_case1%heat_budget_vdiff == ms_case3%heat_budget_vdiff), &
                    "heat_budget_vdiff differs between single- and two-pass dispatch")
         if (allocated(error)) exit checks
         call check(error, all(ms_case1%salt_budget_vdiff == ms_case3%salt_budget_vdiff), &
                    "salt_budget_vdiff differs between single- and two-pass dispatch")
      end block checks
      if (allocated(k3)) deallocate (k3)
      call vd_case1%destroy(); call vd_case3%destroy()
      call ms_case1%destroy(); call ms_case3%destroy()
   end subroutine test_per_tracer_kt_eq_ks_bit_identical

   subroutine test_passive_tracer_follows_ks(error)
      !! MOM6's `Kd_salt` is "the diapycnal diffusivity of salt AND
      !! PASSIVE TRACERS" (`MOM_diabatic_driver.F90:579`).  Register the
      !! ideal-age tracer (the multilayer-C-grid registry's passive
      !! tracer) alongside T and S, seed all three with the same
      !! profile, run with `ks_source = 4*kt_source`: the passive must
      !! match salinity exactly and diverge from temperature.  Pins the
      !! convention against a future refactor flipping it.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vdiff_t) :: vd
      real(wp), parameter :: H0 = 10.0_wp, DT = 60.0_wp
      real(wp), parameter :: KT0 = 1.0e-2_wp, KS0 = 4.0_wp*KT0
      real(wp), allocatable :: kt3(:, :, :), ks3(:, :, :)
      real(wp) :: profile(NZ)
      integer :: k, nx, ny, ic, jc

      checks: block
         call make_grid(grid, 6, 4)
         ms%nz_ml = NZ
         call ms%init(grid, with_ideal_age=.true.)
         call vd%init(grid, nz_ml=NZ)
         nx = grid%nx_total
         ny = grid%ny_total

         ms%h_layer = H0
         do k = 1, NZ
            profile(k) = real(k, wp)
         end do
         do k = 1, NZ
            ms%tracers(ms%idx_temperature)%hTr(:, :, k) = profile(k)*H0
            ms%tracers(ms%idx_salinity)%hTr(:, :, k) = profile(k)*H0
            ms%tracers(ms%idx_age)%hTr(:, :, k) = profile(k)*H0
         end do

         allocate (kt3(nx, ny, NZ + 1), source=0.0_wp)
         allocate (ks3(nx, ny, NZ + 1), source=0.0_wp)
         do k = 2, NZ
            kt3(:, :, k) = KT0
            ks3(:, :, k) = KS0
         end do

         call map_in(ms, vd)
         !$acc enter data copyin(kt3, ks3)
         call vdiff_apply_tracers(grid, vd, ms, DT, kt_source=kt3, ks_source=ks3)
         !$acc exit data delete(kt3, ks3)
         call map_out(ms, vd)

         ic = grid%nghost + 2
         jc = grid%nghost + 2
         call check(error, all(ms%tracers(ms%idx_age)%hTr(ic, jc, :) == &
                               ms%tracers(ms%idx_salinity)%hTr(ic, jc, :)), &
                    "passive tracer (age) did not follow ks_source like salinity")
         if (allocated(error)) exit checks
         call check(error, any(ms%tracers(ms%idx_age)%hTr(ic, jc, :) /= &
                               ms%tracers(ms%idx_temperature)%hTr(ic, jc, :)), &
                    "passive tracer (age) unexpectedly matched kt_source (temperature)")
      end block checks
      if (allocated(kt3)) deallocate (kt3)
      if (allocated(ks3)) deallocate (ks3)
      call vd%destroy(); call ms%destroy()
   end subroutine test_passive_tracer_follows_ks

end module test_ocean_vdiff
