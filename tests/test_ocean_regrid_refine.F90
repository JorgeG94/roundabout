!! ALE-regrid refinement tests (P5): grid time-filter + KE-conserving
!! velocity rescale.  Both refinements default-off and bit-identical
!! when off; these tests drive them through `ocean_apply_ale_remap_step`
!! / the face-velocity remap.
!!
!! Cases (analytical — the high-leverage kind):
!!   1. regrid_tfilter_off_bit_identical — `regrid_time_scale = 0`
!!      ⇒ post-remap `h_layer` is byte-identical to the no-filter result.
!!   2. regrid_tfilter_relaxes — `τ > 0` ⇒ each layer sits a fraction
!!      `dt/(τ+dt)` of the way from old toward target (between the two),
!!      and the column total is conserved.
!!   3. regrid_tfilter_conserves — T·h / S·h conserved with the filter on.
!!   4. remap_vel_ke_off_bit_identical — `remap_vel_conserve_ke = .false.`
!!      ⇒ remapped face velocities unchanged vs the baseline remap.
!!   5. remap_vel_ke_conserves — flag on ⇒ per-face column KE preserved
!!      (within the 1.25× cap) and the barotropic mean `Σh·u/Σh`
!!      unchanged.
!!   6. regrid_refine_on_device — both knobs exercised through a device
!!      enter_data round-trip (GPU exercise): finite + conservative.
!!
!! All cases use VCOORD_SIGMA so the target grid is the simple uniform
!! `(H+η)·dsig` stencil — making the old/target spread an exact closed
!! form for the time-filter prediction.
module test_ocean_regrid_refine
   use rdb_constants, only: wp, REMAP_PPM, VCOORD_SIGMA
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_vcoord, only: ocean_vcoord_t
   use rdb_ocean_remap, only: ocean_apply_ale_remap_step, &
                              ocean_apply_ale_remap_faces
   use testdrive, only: error_type, check, new_unittest, unittest_type
   implicit none
   private

   public :: collect_ocean_regrid_refine_tests

contains

   subroutine collect_ocean_regrid_refine_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("regrid_tfilter_off_bit_identical", test_tfilter_off), &
                  new_unittest("regrid_tfilter_relaxes", test_tfilter_relaxes), &
                  new_unittest("regrid_tfilter_conserves", test_tfilter_conserves), &
                  new_unittest("remap_vel_ke_off_bit_identical", test_vel_ke_off), &
                  new_unittest("remap_vel_ke_conserves", test_vel_ke_conserves), &
                  new_unittest("regrid_refine_on_device", test_on_device) &
                  ]
   end subroutine collect_ocean_regrid_refine_tests

   ! -----------------------------------------------------------------
   ! Helpers
   ! -----------------------------------------------------------------

   subroutine setup_state(grid, ms, nx, ny, nz, layer_pattern, s_ref, t_ref)
      !! Build a uniform-over-(i,j) multilayer state from a bottom-up
      !! per-layer thickness pattern (k=1 bed .. k=nz surface) and
      !! uniform S/T concentrations.
      type(hgrid_t), intent(inout) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: layer_pattern(nz)
      real(wp), intent(in) :: s_ref, t_ref
      integer :: i, j, k
      call grid%init(nx, ny, 1, 1.0_wp, 1.0_wp)
      ms%nz_ml = nz
      call ms%init(grid)
      do k = 1, nz
         do j = 1, grid%ny_total
            do i = 1, grid%nx_total
               ms%h_layer(i, j, k) = layer_pattern(k)
               ms%tracers(ms%idx_salinity)%hTr(i, j, k) = s_ref*layer_pattern(k)
               ms%tracers(ms%idx_temperature)%hTr(i, j, k) = t_ref*layer_pattern(k)
            end do
         end do
      end do
   end subroutine setup_state

   ! -----------------------------------------------------------------
   ! 1. time-filter off (τ=0) ⇒ byte-identical to no-filter remap
   ! -----------------------------------------------------------------
   subroutine test_tfilter_off(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid_a, grid_b
      type(multilayer_state_t) :: ms_a, ms_b
      type(ocean_vcoord_t) :: vc_a, vc_b
      integer, parameter :: NX = 3, NY = 3, NZ = 4
      real(wp), parameter :: H_TOTAL = 400.0_wp, S0 = 35.0_wp, T0 = 12.0_wp
      real(wp), parameter :: DT = 1200.0_wp
      real(wp) :: pattern(NZ)
      real(wp), allocatable :: bt_eta_a(:, :), bt_H_ref_a(:, :)
      real(wp), allocatable :: bt_eta_b(:, :), bt_H_ref_b(:, :)
      real(wp) :: max_diff
      integer :: nx_t, ny_t
      checks: block
         pattern = [90.0_wp, 95.0_wp, 105.0_wp, 110.0_wp]

         ! Baseline: no dt argument (filter cannot fire).
         call setup_state(grid_a, ms_a, NX, NY, NZ, pattern, S0, T0)
         call vc_a%init(grid_a, nz_ml=NZ)
         vc_a%coord_type = VCOORD_SIGMA
         nx_t = grid_a%nx_total
         ny_t = grid_a%ny_total
         allocate (bt_eta_a(nx_t, ny_t), source=0.0_wp)
         allocate (bt_H_ref_a(nx_t, ny_t), source=H_TOTAL)
         call ocean_apply_ale_remap_step(grid_a, vc_a, ms_a, bt_eta_a, bt_H_ref_a, &
                                         method=REMAP_PPM)

         ! Filter present but τ=0 + dt passed: wtd=1 ⇒ jump to target.
         call setup_state(grid_b, ms_b, NX, NY, NZ, pattern, S0, T0)
         call vc_b%init(grid_b, nz_ml=NZ)
         vc_b%coord_type = VCOORD_SIGMA
         vc_b%regrid_time_scale = 0.0_wp
         allocate (bt_eta_b(nx_t, ny_t), source=0.0_wp)
         allocate (bt_H_ref_b(nx_t, ny_t), source=H_TOTAL)
         call ocean_apply_ale_remap_step(grid_b, vc_b, ms_b, bt_eta_b, bt_H_ref_b, &
                                         method=REMAP_PPM, dt=DT)

         max_diff = maxval(abs(ms_a%h_layer - ms_b%h_layer))
         call check(error, max_diff == 0.0_wp, &
                    "tfilter off: h_layer must be byte-identical to no-filter remap")
      end block checks
      if (allocated(bt_eta_a)) deallocate (bt_eta_a, bt_H_ref_a)
      if (allocated(bt_eta_b)) deallocate (bt_eta_b, bt_H_ref_b)
      call vc_a%destroy()
      call ms_a%destroy()
      call vc_b%destroy()
      call ms_b%destroy()
   end subroutine test_tfilter_off

   ! -----------------------------------------------------------------
   ! 2. time-filter relaxes: layer sits dt/(τ+dt) toward target
   ! -----------------------------------------------------------------
   subroutine test_tfilter_relaxes(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vcoord_t) :: vc
      integer, parameter :: NX = 3, NY = 3, NZ = 4
      real(wp), parameter :: H_TOTAL = 400.0_wp, S0 = 35.0_wp, T0 = 12.0_wp
      real(wp), parameter :: DT = 1200.0_wp, TAU = 3600.0_wp
      real(wp) :: pattern(NZ)
      real(wp) :: wtd, target_uniform, expected, col_sum
      real(wp), allocatable :: bt_eta(:, :), bt_H_ref(:, :)
      real(wp) :: max_err, lo, hi
      integer :: nx_t, ny_t, k
      logical :: between
      checks: block
         pattern = [90.0_wp, 95.0_wp, 105.0_wp, 110.0_wp]
         call setup_state(grid, ms, NX, NY, NZ, pattern, S0, T0)
         call vc%init(grid, nz_ml=NZ)
         vc%coord_type = VCOORD_SIGMA
         vc%regrid_time_scale = TAU
         nx_t = grid%nx_total
         ny_t = grid%ny_total
         allocate (bt_eta(nx_t, ny_t), source=0.0_wp)
         allocate (bt_H_ref(nx_t, ny_t), source=H_TOTAL)

         call ocean_apply_ale_remap_step(grid, vc, ms, bt_eta, bt_H_ref, &
                                         method=REMAP_PPM, dt=DT)

         ! SIGMA target is uniform H_TOTAL/NZ; blended h = old + wtd*(tgt-old).
         wtd = DT/(TAU + DT)
         target_uniform = H_TOTAL/real(NZ, wp)
         max_err = 0.0_wp
         between = .true.
         do k = 1, NZ
            expected = pattern(k) + wtd*(target_uniform - pattern(k))
            max_err = max(max_err, abs(ms%h_layer(1, 1, k) - expected))
            ! The relaxed interface must sit strictly between old and target.
            lo = min(pattern(k), target_uniform)
            hi = max(pattern(k), target_uniform)
            if (ms%h_layer(1, 1, k) < lo - 1.0e-9_wp .or. &
                ms%h_layer(1, 1, k) > hi + 1.0e-9_wp) between = .false.
         end do
         call check(error, max_err < 1.0e-9_wp, &
                    "tfilter relaxes: h_layer = old + dt/(tau+dt)*(target-old)")
         if (allocated(error)) exit checks
         call check(error, between, "tfilter relaxes: h_layer between old and target")
         if (allocated(error)) exit checks

         ! Column total conserved (convex blend of two grids summing to H).
         col_sum = sum(ms%h_layer(1, 1, :))
         call check(error, abs(col_sum - H_TOTAL) < 1.0e-9_wp, &
                    "tfilter relaxes: column total conserved")
      end block checks
      if (allocated(bt_eta)) deallocate (bt_eta, bt_H_ref)
      call vc%destroy()
      call ms%destroy()
   end subroutine test_tfilter_relaxes

   ! -----------------------------------------------------------------
   ! 3. time-filter conserves tracer mass
   ! -----------------------------------------------------------------
   subroutine test_tfilter_conserves(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vcoord_t) :: vc
      integer, parameter :: NX = 4, NY = 3, NZ = 4
      real(wp), parameter :: H_TOTAL = 400.0_wp, S0 = 35.0_wp, T0 = 12.0_wp
      real(wp), parameter :: DT = 1200.0_wp, TAU = 7200.0_wp
      real(wp) :: pattern(NZ)
      real(wp), allocatable :: bt_eta(:, :), bt_H_ref(:, :)
      real(wp) :: Th0, Sh0, Th1, Sh1
      integer :: nx_t, ny_t, k
      checks: block
         ! Non-uniform T/S per layer so the remap actually transports mass.
         pattern = [80.0_wp, 100.0_wp, 110.0_wp, 110.0_wp]
         call setup_state(grid, ms, NX, NY, NZ, pattern, S0, T0)
         do k = 1, NZ
            ms%tracers(ms%idx_temperature)%hTr(:, :, k) = &
               (T0 + 2.0_wp*real(k, wp))*pattern(k)
            ms%tracers(ms%idx_salinity)%hTr(:, :, k) = &
               (S0 - 0.5_wp*real(k, wp))*pattern(k)
         end do
         Th0 = sum(ms%tracers(ms%idx_temperature)%hTr(1, 1, :))
         Sh0 = sum(ms%tracers(ms%idx_salinity)%hTr(1, 1, :))

         call vc%init(grid, nz_ml=NZ)
         vc%coord_type = VCOORD_SIGMA
         vc%regrid_time_scale = TAU
         nx_t = grid%nx_total
         ny_t = grid%ny_total
         allocate (bt_eta(nx_t, ny_t), source=0.0_wp)
         allocate (bt_H_ref(nx_t, ny_t), source=H_TOTAL)
         call ocean_apply_ale_remap_step(grid, vc, ms, bt_eta, bt_H_ref, &
                                         method=REMAP_PPM, dt=DT)

         Th1 = sum(ms%tracers(ms%idx_temperature)%hTr(1, 1, :))
         Sh1 = sum(ms%tracers(ms%idx_salinity)%hTr(1, 1, :))
         call check(error, abs(Th1 - Th0) < 1.0e-7_wp, "tfilter conserves: T*h")
         if (allocated(error)) exit checks
         call check(error, abs(Sh1 - Sh0) < 1.0e-7_wp, "tfilter conserves: S*h")
      end block checks
      if (allocated(bt_eta)) deallocate (bt_eta, bt_H_ref)
      call vc%destroy()
      call ms%destroy()
   end subroutine test_tfilter_conserves

   ! -----------------------------------------------------------------
   ! Face-velocity test setup: sheared u/v on a variable-thickness column
   ! -----------------------------------------------------------------
   subroutine setup_faces(grid, h_old, h_new, u_face_x, v_face_y, nx, ny, nz)
      !! Variable old/new thicknesses (so the remap actually moves
      !! velocity) and a sheared face-velocity profile.
      type(hgrid_t), intent(inout) :: grid
      real(wp), allocatable, intent(out) :: h_old(:, :, :), h_new(:, :, :)
      real(wp), allocatable, intent(out) :: u_face_x(:, :, :), v_face_y(:, :, :)
      integer, intent(in) :: nx, ny, nz
      integer :: i, j, k, nx_t, ny_t
      call grid%init(nx, ny, 1, 1.0_wp, 1.0_wp)
      nx_t = grid%nx_total
      ny_t = grid%ny_total
      allocate (h_old(nx_t, ny_t, nz), h_new(nx_t, ny_t, nz))
      allocate (u_face_x(nx_t + 1, ny_t, nz))
      allocate (v_face_y(nx_t, ny_t + 1, nz))
      do k = 1, nz
         do j = 1, ny_t
            do i = 1, nx_t
               h_old(i, j, k) = 80.0_wp + 20.0_wp*real(k, wp)
               h_new(i, j, k) = 100.0_wp        ! uniform target
            end do
         end do
         ! Sheared velocities (k=1 bed .. k=nz surface).
         u_face_x(:, :, k) = 0.1_wp + 0.3_wp*real(k - 1, wp)
         v_face_y(:, :, k) = -0.2_wp + 0.15_wp*real(k - 1, wp)
      end do
   end subroutine setup_faces

   pure function col_anom_ke(nz, h, u) result(ke)
      !! Σ ½ h·(u - u_bar)² over a column, where u_bar = Σh·u/Σh.  This is
      !! the BAROCLINIC anomaly KE — the quantity the rescale conserves
      !! (the barotropic KE ½·Σh·u_bar² changes with the grid and is left
      !! to momentum conservation).
      integer, intent(in) :: nz
      real(wp), intent(in) :: h(nz), u(nz)
      real(wp) :: ke, ub, anom
      integer :: k
      ub = col_mean(nz, h, u)
      ke = 0.0_wp
      do k = 1, nz
         anom = u(k) - ub
         ke = ke + 0.5_wp*h(k)*anom*anom
      end do
   end function col_anom_ke

   pure function col_mean(nz, h, u) result(ub)
      !! Σ(h·u)/Σh — the barotropic / depth-mean velocity.
      integer, intent(in) :: nz
      real(wp), intent(in) :: h(nz), u(nz)
      real(wp) :: ub, hsum, msum
      integer :: k
      hsum = 0.0_wp
      msum = 0.0_wp
      do k = 1, nz
         hsum = hsum + h(k)
         msum = msum + h(k)*u(k)
      end do
      ub = msum/hsum
   end function col_mean

   ! -----------------------------------------------------------------
   ! 4. KE-conserve off ⇒ face velocities unchanged vs baseline
   ! -----------------------------------------------------------------
   subroutine test_vel_ke_off(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      integer, parameter :: NX = 3, NY = 3, NZ = 4
      real(wp), allocatable :: h_old(:, :, :), h_new(:, :, :)
      real(wp), allocatable :: ux_a(:, :, :), vy_a(:, :, :)
      real(wp), allocatable :: ux_b(:, :, :), vy_b(:, :, :)
      checks: block
         call setup_faces(grid, h_old, h_new, ux_a, vy_a, NX, NY, NZ)
         allocate (ux_b, source=ux_a)
         allocate (vy_b, source=vy_a)

         ! Baseline: no conserve_ke argument.
         call ocean_apply_ale_remap_faces(grid, h_old, h_new, ux_a, vy_a, method=REMAP_PPM)
         ! Explicit conserve_ke = .false. must match byte-for-byte.
         call ocean_apply_ale_remap_faces(grid, h_old, h_new, ux_b, vy_b, &
                                          method=REMAP_PPM, conserve_ke=.false.)

         call check(error, maxval(abs(ux_a - ux_b)) == 0.0_wp, &
                    "ke off: u_face_x byte-identical")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(vy_a - vy_b)) == 0.0_wp, &
                    "ke off: v_face_y byte-identical")
      end block checks
      if (allocated(h_old)) deallocate (h_old, h_new, ux_a, vy_a, ux_b, vy_b)
   end subroutine test_vel_ke_off

   ! -----------------------------------------------------------------
   ! 5. KE-conserve on ⇒ column KE preserved + barotropic mean unchanged
   ! -----------------------------------------------------------------
   subroutine test_vel_ke_conserves(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      integer, parameter :: NX = 3, NY = 3, NZ = 4
      real(wp), allocatable :: h_old(:, :, :), h_new(:, :, :)
      real(wp), allocatable :: ux(:, :, :), vy(:, :, :)         ! KE-conserve run
      real(wp), allocatable :: ux0(:, :, :), vy0(:, :, :)       ! pristine input
      real(wp), allocatable :: ux_b(:, :, :), vy_b(:, :, :)     ! momentum-only baseline
      real(wp) :: hold_face(NZ), hnew_face(NZ)
      real(wp) :: ke_old, ke_new, ke_base, ub_base, ub_ke
      integer :: I, J, k, nx_t, ny_t
      logical :: ke_ok, bar_ok
      checks: block
         call setup_faces(grid, h_old, h_new, ux, vy, NX, NY, NZ)
         allocate (ux0, source=ux)
         allocate (vy0, source=vy)
         allocate (ux_b, source=ux)
         allocate (vy_b, source=vy)
         nx_t = grid%nx_total
         ny_t = grid%ny_total

         ! Momentum-only baseline remap (conserve_ke off).
         call ocean_apply_ale_remap_faces(grid, h_old, h_new, ux_b, vy_b, method=REMAP_PPM)
         ! KE-conserving remap on the same input.
         call ocean_apply_ale_remap_faces(grid, h_old, h_new, ux, vy, &
                                          method=REMAP_PPM, conserve_ke=.true.)

         ! Per interior x-face:
         !   - the BAROCLINIC anomaly KE on the NEW grid matches the
         !     PRE-remap anomaly KE on the OLD grid (the rescale restores
         !     what the limiter drained), uncapped here because the
         !     profile is mild.
         !   - the post-remap barotropic mean (Σ h_new·u/Σ h_new) is
         !     IDENTICAL between the baseline and the KE run — the rescale
         !     never touches the mean (mode-split consistency).
         ke_ok = .true.
         bar_ok = .true.
         do J = 1, ny_t
            do I = 2, nx_t            ! interior faces (avg of I-1, I)
               do k = 1, NZ
                  hold_face(k) = 0.5_wp*(h_old(I - 1, J, k) + h_old(I, J, k))
                  hnew_face(k) = 0.5_wp*(h_new(I - 1, J, k) + h_new(I, J, k))
               end do
               ke_old = col_anom_ke(NZ, hold_face, ux0(I, J, :))      ! pre-remap target
               ke_base = col_anom_ke(NZ, hnew_face, ux_b(I, J, :))    ! drained (momentum-only)
               ke_new = col_anom_ke(NZ, hnew_face, ux(I, J, :))       ! after KE-rescale
               ub_base = col_mean(NZ, hnew_face, ux_b(I, J, :))
               ub_ke = col_mean(NZ, hnew_face, ux(I, J, :))
               ! The rescale restores anomaly KE TOWARD the pre-remap value:
               !   ke_base (drained) <= ke_new <= ke_old (target).
               ! Exact (ke_new==ke_old) when the needed factor <= the 1.25x cap;
               ! cap-limited (ke_base < ke_new < ke_old) on a steep column.
               ! Either way it must not overshoot the target nor fall below the
               ! drained baseline, and must strictly increase it (restored some).
               if (ke_new > ke_old*(1.0_wp + 1.0e-9_wp)) ke_ok = .false.
               if (ke_new < ke_base*(1.0_wp - 1.0e-9_wp)) ke_ok = .false.
               if (ke_old > ke_base + 1.0e-10_wp .and. &
                   ke_new <= ke_base + 1.0e-12_wp) ke_ok = .false.
               if (abs(ub_ke - ub_base) > 1.0e-12_wp*max(abs(ub_base), 1.0_wp)) bar_ok = .false.
            end do
         end do
         call check(error, ke_ok, "ke on: per-face baroclinic anomaly KE restored to pre-remap value")
         if (allocated(error)) exit checks
         call check(error, bar_ok, "ke on: barotropic mean unchanged vs momentum-only baseline")
         if (allocated(error)) exit checks

         ! Sanity: the rescale actually changed the velocities (i.e. the
         ! baseline PPM remap did drain KE), so the test is meaningful.
         call check(error, maxval(abs(ux - ux_b)) > 0.0_wp, &
                    "ke on: rescale must alter the remapped velocities")
      end block checks
      if (allocated(h_old)) deallocate (h_old, h_new, ux, vy, ux0, vy0, ux_b, vy_b)
   end subroutine test_vel_ke_conserves

   ! -----------------------------------------------------------------
   ! 6. Both knobs through a device enter_data round-trip
   ! -----------------------------------------------------------------
   subroutine test_on_device(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vcoord_t) :: vc
      integer, parameter :: NX = 3, NY = 3, NZ = 4
      real(wp), parameter :: H_TOTAL = 400.0_wp, S0 = 35.0_wp, T0 = 12.0_wp
      real(wp), parameter :: DT = 1200.0_wp, TAU = 3600.0_wp
      real(wp) :: pattern(NZ)
      real(wp) :: col_sum, Th0, Th1
      real(wp), allocatable :: bt_eta(:, :), bt_H_ref(:, :)
      integer :: nx_t, ny_t, i, j, k
      logical :: finite
      checks: block
         pattern = [80.0_wp, 100.0_wp, 110.0_wp, 110.0_wp]
         call setup_state(grid, ms, NX, NY, NZ, pattern, S0, T0)
         ! Sheared face velocities so the KE rescale has something to do.
         do k = 1, NZ
            ms%u_face_x_layer(:, :, k) = 0.1_wp + 0.3_wp*real(k - 1, wp)
            ms%v_face_y_layer(:, :, k) = -0.2_wp + 0.15_wp*real(k - 1, wp)
            ms%tracers(ms%idx_temperature)%hTr(:, :, k) = &
               (T0 + 2.0_wp*real(k, wp))*pattern(k)
         end do
         Th0 = sum(ms%tracers(ms%idx_temperature)%hTr(1, 1, :))

         call vc%init(grid, nz_ml=NZ)
         vc%coord_type = VCOORD_SIGMA
         vc%regrid_time_scale = TAU
         vc%remap_vel_conserve_ke = .true.
         nx_t = grid%nx_total
         ny_t = grid%ny_total
         allocate (bt_eta(nx_t, ny_t), source=0.0_wp)
         allocate (bt_H_ref(nx_t, ny_t), source=H_TOTAL)

         !$acc enter data copyin(ms, vc, bt_eta, bt_H_ref)
         call ms%enter_data()
         call vc%enter_data()
         call ocean_apply_ale_remap_step(grid, vc, ms, bt_eta, bt_H_ref, &
                                         method=REMAP_PPM, dt=DT)
         !$acc update self(ms%h_layer, ms%u_face_x_layer, ms%v_face_y_layer)
         !$acc update self(ms%tracers(ms%idx_temperature)%hTr)
         call vc%exit_data()
         call ms%exit_data()
         !$acc exit data delete(ms, vc, bt_eta, bt_H_ref)

         finite = .true.
         do k = 1, NZ
            if (ms%h_layer(1, 1, k) /= ms%h_layer(1, 1, k)) finite = .false.
            if (ms%h_layer(1, 1, k) < 0.0_wp) finite = .false.
         end do
         call check(error, finite, "device: thicknesses finite + non-negative")
         if (allocated(error)) exit checks
         col_sum = sum(ms%h_layer(1, 1, :))
         call check(error, abs(col_sum - H_TOTAL) < 1.0e-7_wp, &
                    "device: column total conserved")
         if (allocated(error)) exit checks
         Th1 = sum(ms%tracers(ms%idx_temperature)%hTr(1, 1, :))
         call check(error, abs(Th1 - Th0) < 1.0e-6_wp, "device: T*h conserved")
      end block checks
      if (allocated(bt_eta)) deallocate (bt_eta, bt_H_ref)
      call vc%destroy()
      call ms%destroy()
   end subroutine test_on_device

end module test_ocean_regrid_refine
