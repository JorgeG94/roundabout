!! Unit tests for the upstream-PPM h_face slot + builder + derive_bt
!! switch (steps 4-6).
!!
!! Covers:
!!   1. `builder_flat_bath`         — uniform h, sign-independent, sum
!!                                    equals nz·h_uniform.
!!   2. `builder_spoon_margin`      — west bed vanished, east bed=30;
!!                                    u sign chooses the upstream
!!                                    column, sums differ.
!!   3. `derive_bt_from_layers_upstream`
!!                                  — sheared flow with vanished west
!!                                    bed gives bt_ubt biased toward
!!                                    the surface (only the surface
!!                                    upstream-h carries weight).
!!   4. `builder_off_is_noop`       — knob off ⇒ slots stay
!!                                    unallocated, routine doesn't
!!                                    touch anything.
!!   5. `derive_bt_off_matches_centred`
!!                                  — knob off ⇒ derive_bt_from_layers
!!                                    produces the same bt_ubt as the
!!                                    pre-knob centred path.
!!
!! All tests use a 2×1 ocean (nx_phys=2, ny_phys=1, nghost=0) with the
!! single interior u-face at `i=2`.  Boundary u-faces `i=1` and `i=3`
!! fall back to single-cell h (matches `derive_bt_from_layers`).
module test_ocean_bt_upstream_h_face
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_barotropic_workstate, only: barotropic_workstate_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_barotropic_coupling, only: compute_h_face_upstream, &
                                      derive_bt_from_layers, &
                                      apply_bt_correction
   use rdb_barotropic_substep, only: barotropic_substep_linear
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   implicit none
   private

   public :: collect_ocean_bt_upstream_h_face_tests

contains

   subroutine collect_ocean_bt_upstream_h_face_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("builder_flat_bath", test_builder_flat_bath), &
                  new_unittest("builder_spoon_margin", test_builder_spoon_margin), &
                  new_unittest("derive_bt_from_layers_upstream", &
                               test_derive_bt_upstream), &
                  new_unittest("builder_off_is_noop", test_builder_off_is_noop), &
                  new_unittest("derive_bt_off_matches_centred", &
                               test_derive_bt_off_matches_centred), &
                  new_unittest("bt_substep_linear_uniform_matches_centred", &
                               test_substep_linear_uniform_match), &
                  new_unittest("apply_bt_correction_uniform_depth_mean", &
                               test_apply_bt_corr_depth_mean) &
                  ]
   end subroutine collect_ocean_bt_upstream_h_face_tests

   subroutine make_setup(grid, ms, bt_work, nz, h_west_per_k, h_east_per_k, &
                         u_face_per_k, enable_upstream)
      !! Spin up a 2×1×nz multilayer state with two cells (west at
      !! `i=1`, east at `i=2`) and a single interior u-face at `i=2`.
      !! Caller supplies per-layer h on each side and the per-layer
      !! u_face_x at i=2; the helper builds the grid, allocates the
      !! multilayer state, fills h_layer + u_face_x_layer, allocates
      !! bt_work, sets `use_upstream_h_face`, and lazy-allocates the
      !! `h_face_up_x/y` slots when the knob is on.
      type(hgrid_t), intent(out) :: grid
      type(multilayer_state_t), intent(out) :: ms
      type(barotropic_workstate_t), intent(out) :: bt_work
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_west_per_k(nz)
      real(wp), intent(in) :: h_east_per_k(nz)
      real(wp), intent(in) :: u_face_per_k(nz)
      logical, intent(in) :: enable_upstream

      integer :: nx_phys, ny_phys, nghost, k

      nx_phys = 2
      ny_phys = 1
      nghost = 0
      call grid%init(nx_phys, ny_phys, nghost, 1000.0_wp, 1000.0_wp)

      ms%nz_ml = nz
      call ms%init(grid)

      do k = 1, nz
         ms%h_layer(1, 1, k) = h_west_per_k(k)
         ms%h_layer(2, 1, k) = h_east_per_k(k)
         ! Interior u-face is at i=2 (between west cell 1 and east 2).
         ms%u_face_x_layer(2, 1, k) = u_face_per_k(k)
      end do

      call bt_work%init(grid, nz_ml=nz)
      bt_work%use_upstream_h_face = enable_upstream
      if (enable_upstream) then
         allocate (bt_work%h_face_up_x(grid%nx_total + 1, grid%ny_total), &
                   source=0.0_wp)
         allocate (bt_work%h_face_up_y(grid%nx_total, grid%ny_total + 1), &
                   source=0.0_wp)
      end if
   end subroutine make_setup

   subroutine test_builder_flat_bath(error)
      !! Uniform `h_layer = 1000` over NZ=3 layers in both columns.
      !! Expected at the interior u-face (i=2):
      !!   h_face_up_x(2, 1) = 3·1000 = 3000 independent of u sign.
      !! Boundary faces fall back to single-cell h: same value.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(barotropic_workstate_t) :: bt_work
      real(wp), parameter :: h_uniform = 1000.0_wp
      real(wp), parameter :: H_total = 3.0_wp*h_uniform

      ! Try u > 0 (upstream = west).
      call make_setup(grid, ms, bt_work, 3, &
                      [h_uniform, h_uniform, h_uniform], &
                      [h_uniform, h_uniform, h_uniform], &
                      [0.1_wp, 0.1_wp, 0.1_wp], &
                      enable_upstream=.true.)
      call compute_h_face_upstream(grid, bt_work, ms)
      call check(error, abs(bt_work%h_face_up_x(2, 1) - H_total) < 1.0e-10_wp, &
                 "h_face_up_x ≠ nz·h_uniform (u>0)")
      if (allocated(error)) return
      ! Boundary faces — fall back to single-cell column.
      call check(error, abs(bt_work%h_face_up_x(1, 1) - H_total) < 1.0e-10_wp, &
                 "h_face_up_x boundary i=1 ≠ nz·h_uniform")
      if (allocated(error)) return
      call check(error, abs(bt_work%h_face_up_x(3, 1) - H_total) < 1.0e-10_wp, &
                 "h_face_up_x boundary i=3 ≠ nz·h_uniform")
      if (allocated(error)) return
      call ms%destroy()
      call bt_work%destroy()

      ! Repeat with u < 0 (upstream = east) — must give the same answer.
      call make_setup(grid, ms, bt_work, 3, &
                      [h_uniform, h_uniform, h_uniform], &
                      [h_uniform, h_uniform, h_uniform], &
                      [-0.1_wp, -0.1_wp, -0.1_wp], &
                      enable_upstream=.true.)
      call compute_h_face_upstream(grid, bt_work, ms)
      call check(error, abs(bt_work%h_face_up_x(2, 1) - H_total) < 1.0e-10_wp, &
                 "h_face_up_x ≠ nz·h_uniform (u<0)")
      call ms%destroy()
      call bt_work%destroy()
   end subroutine test_builder_flat_bath

   subroutine test_builder_spoon_margin(error)
      !! West column has the bed layer vanished (h_bed=0, h_surf=1000);
      !! east column has both layers present (h_bed=30, h_surf=1000).
      !! Per-layer u_face signs chosen to exercise both branches in
      !! the same setup:
      !!   k=1 (bed):    u >  0  ⇒ upstream = west ⇒ h_face_k =    0
      !!   k=2 (surf):   u <  0  ⇒ upstream = east ⇒ h_face_k = 1000
      !! Expected h_face_up_x(2, 1) = 0 + 1000 = 1000.
      !!
      !! Then re-run with both u flipped to fully eastward (u>0 in
      !! both layers) and check the all-west result: h_face_up_x =
      !! 0 (bed) + 1000 (surf) = 1000.  And westward (u<0 in both):
      !! h_face_up_x = 30 + 1000 = 1030.  Sign matters.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(barotropic_workstate_t) :: bt_work

      ! Case A: mixed signs (bed eastward, surf westward).
      call make_setup(grid, ms, bt_work, 2, &
                      [0.0_wp, 1000.0_wp], &       ! west: bed vanished
                      [30.0_wp, 1000.0_wp], &      ! east: full
                      [0.1_wp, -0.1_wp], &         ! u_face per k
                      enable_upstream=.true.)
      call compute_h_face_upstream(grid, bt_work, ms)
      ! k=1 west (h=0) + k=2 east (h=1000) = 1000.
      call check(error, abs(bt_work%h_face_up_x(2, 1) - 1000.0_wp) < 1.0e-10_wp, &
                 "mixed-sign upstream pick ≠ 0 + 1000")
      if (allocated(error)) return
      call ms%destroy()
      call bt_work%destroy()

      ! Case B: all eastward (u > 0 both layers) ⇒ all west upstream.
      call make_setup(grid, ms, bt_work, 2, &
                      [0.0_wp, 1000.0_wp], &
                      [30.0_wp, 1000.0_wp], &
                      [0.1_wp, 0.1_wp], &
                      enable_upstream=.true.)
      call compute_h_face_upstream(grid, bt_work, ms)
      ! West column sum: 0 + 1000 = 1000.
      call check(error, abs(bt_work%h_face_up_x(2, 1) - 1000.0_wp) < 1.0e-10_wp, &
                 "all-eastward upstream pick ≠ west-column sum")
      if (allocated(error)) return
      call ms%destroy()
      call bt_work%destroy()

      ! Case C: all westward (u < 0 both layers) ⇒ all east upstream.
      call make_setup(grid, ms, bt_work, 2, &
                      [0.0_wp, 1000.0_wp], &
                      [30.0_wp, 1000.0_wp], &
                      [-0.1_wp, -0.1_wp], &
                      enable_upstream=.true.)
      call compute_h_face_upstream(grid, bt_work, ms)
      ! East column sum: 30 + 1000 = 1030.
      call check(error, abs(bt_work%h_face_up_x(2, 1) - 1030.0_wp) < 1.0e-10_wp, &
                 "all-westward upstream pick ≠ east-column sum")
      if (allocated(error)) return
      call ms%destroy()
      call bt_work%destroy()
   end subroutine test_builder_spoon_margin

   subroutine test_derive_bt_upstream(error)
      !! Sheared flow on a vanished-bed margin.  West: h_bed=0,
      !! h_surf=1000; east: h_bed=30, h_surf=1000.  Velocity:
      !! `u_bed = +0.30, u_surf = +0.05` (both eastward, so upstream
      !! is west in both layers).
      !!
      !! Centred-h depth-mean (NOT what we're testing — for reference):
      !!   h_face_bed  = 0.5·(0 + 30)    =   15
      !!   h_face_surf = 0.5·(1000+1000) = 1000
      !!   bt_ubt = (0.30·15 + 0.05·1000)/(15+1000) ≈ 0.0537
      !!
      !! Upstream-h depth-mean (the test):
      !!   h_face_bed  = h_west_bed  =    0  (vanished)
      !!   h_face_surf = h_west_surf = 1000
      !!   bt_ubt = (0.30·0 + 0.05·1000)/(0+1000) = 0.05 exactly.
      !!
      !! The bed-layer velocity is suppressed because the upstream
      !! west column has no bed-layer thickness to carry it.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(barotropic_workstate_t) :: bt_work
      type(ocean_metrics_t) :: metrics

      call make_setup(grid, ms, bt_work, 2, &
                      [0.0_wp, 1000.0_wp], &        ! west h: bed=0, surf=1000
                      [30.0_wp, 1000.0_wp], &       ! east h: bed=30, surf=1000
                      [0.30_wp, 0.05_wp], &         ! u_face: bed, surf
                      enable_upstream=.true.)
      ! `metrics` is a REQUIRED argument of `derive_bt_from_layers` — an
      ! optional one is how the closed-face weighting got omitted at a
      ! call site.  All-open metrics (`use_closed_faces = .false.`) ⇒ the
      ! original full-column branch, which is what this test asserts.
      call make_cartesian_metrics(metrics, grid)
      ! bt_H_ref doesn't affect bt_ubt; leave the default zeros.
      call derive_bt_from_layers(grid, bt_work, ms, metrics)
      call check(error, abs(bt_work%bt_ubt(2, 1) - 0.05_wp) < 1.0e-10_wp, &
                 "upstream-h bt_ubt should suppress vanished-bed contribution")
      call destroy_cartesian_metrics(metrics)
      call ms%destroy()
      call bt_work%destroy()
   end subroutine test_derive_bt_upstream

   subroutine test_builder_off_is_noop(error)
      !! Knob off ⇒ `h_face_up_x/y` stay unallocated and the builder
      !! returns without touching anything.  Calling on a workstate
      !! whose upstream slots aren't allocated must not crash.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(barotropic_workstate_t) :: bt_work

      call make_setup(grid, ms, bt_work, 2, &
                      [0.0_wp, 1000.0_wp], &
                      [30.0_wp, 1000.0_wp], &
                      [0.1_wp, -0.1_wp], &
                      enable_upstream=.false.)
      call check(error,.not. allocated(bt_work%h_face_up_x), &
                 "h_face_up_x must not be allocated when knob is off")
      if (allocated(error)) return
      ! Must not crash even with unallocated slots.
      call compute_h_face_upstream(grid, bt_work, ms)
      call check(error,.not. allocated(bt_work%h_face_up_x), &
                 "h_face_up_x must remain unallocated after no-op call")
      call ms%destroy()
      call bt_work%destroy()
   end subroutine test_builder_off_is_noop

   subroutine test_derive_bt_off_matches_centred(error)
      !! Knob off ⇒ `derive_bt_from_layers` produces the same bt_ubt
      !! as the pre-knob centred-h formula.  Bit-identity check for
      !! the existing path.  Same vanished-bed setup; centred-h
      !! reference value computed by hand:
      !!   bt_ubt = (0.30·15 + 0.05·1000)/(15 + 1000)
      !!         = (4.5 + 50)/1015 = 54.5/1015 ≈ 0.05369458...
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(barotropic_workstate_t) :: bt_work
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: expected = 54.5_wp/1015.0_wp

      call make_setup(grid, ms, bt_work, 2, &
                      [0.0_wp, 1000.0_wp], &
                      [30.0_wp, 1000.0_wp], &
                      [0.30_wp, 0.05_wp], &
                      enable_upstream=.false.)
      ! All-open metrics ⇒ the original full-column branch (see
      ! `test_derive_bt_upstream` for why `metrics` is not optional).
      call make_cartesian_metrics(metrics, grid)
      call derive_bt_from_layers(grid, bt_work, ms, metrics)
      call check(error, abs(bt_work%bt_ubt(2, 1) - expected) < 1.0e-12_wp, &
                 "centred-h bt_ubt drifted from analytic reference")
      call destroy_cartesian_metrics(metrics)
      call ms%destroy()
      call bt_work%destroy()
   end subroutine test_derive_bt_off_matches_centred

   subroutine test_substep_linear_uniform_match(error)
      !! Linear BT-substep bit-identity check: with uniform h_layer
      !! everywhere and `bt_eta = 0` initially, `h_face_up_x` =
      !! `0.5·(H_ref_W + H_ref_E)` = `nz·h_uniform` on every face.
      !! Running the substep with the knob ON should give the same
      !! `bt_eta`, `bt_ubt`, `bt_vbt` as running it with the knob OFF.
      !!
      !! Setup: 4×4 grid, NZ=2, h_uniform=500m everywhere, η seed
      !! at (2,2), 10 inner substeps at dt_inner=10s.  Compare
      !! off-path and on-path final time-mean states.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(barotropic_workstate_t) :: bt_off, bt_on
      type(ocean_metrics_t) :: metrics
      integer, parameter :: nx = 4, ny = 4, nz = 2, n_steps = 10
      real(wp), parameter :: h_uniform = 500.0_wp
      real(wp), parameter :: dt_inner = 10.0_wp
      real(wp), allocatable :: force_u(:, :), force_v(:, :)
      integer :: i, j
      real(wp) :: max_diff_eta, max_diff_ubt

      call grid%init(nx, ny, 0, 1000.0_wp, 1000.0_wp)
      call bt_off%init(grid, nz_ml=nz)
      call bt_on%init(grid, nz_ml=nz)
      bt_on%use_upstream_h_face = .true.
      allocate (bt_on%h_face_up_x(nx + 1, ny), source=real(nz, wp)*h_uniform)
      allocate (bt_on%h_face_up_y(nx, ny + 1), source=real(nz, wp)*h_uniform)
      allocate (force_u(nx + 1, ny), source=0.0_wp)
      allocate (force_v(nx, ny + 1), source=0.0_wp)

      ! Static reference column thickness — both paths share this.
      do j = 1, ny
         do i = 1, nx
            bt_off%bt_H_ref(i, j) = real(nz, wp)*h_uniform
            bt_on%bt_H_ref(i, j) = real(nz, wp)*h_uniform
         end do
      end do
      ! Seed a 1 mm gravity-wave perturbation in η.
      bt_off%bt_eta(2, 2) = 1.0e-3_wp
      bt_on%bt_eta(2, 2) = 1.0e-3_wp

      ! Map to device — `-gpu=mem:separate` requires explicit
      ! enter_data; the substep loop touches many bt_work members.
      call make_cartesian_metrics(metrics, grid)
      !$acc enter data copyin(force_u, force_v)
      call bt_off%enter_data()
      call bt_on%enter_data()
      call barotropic_substep_linear(grid, metrics, bt_off, force_u, force_v, n_steps, dt_inner)
      call barotropic_substep_linear(grid, metrics, bt_on, force_u, force_v, n_steps, dt_inner)
      !$acc update self(bt_off%bt_eta, bt_off%bt_ubt, bt_off%bt_vbt)
      !$acc update self(bt_on%bt_eta, bt_on%bt_ubt, bt_on%bt_vbt)
      call bt_off%exit_data()
      call bt_on%exit_data()
      call destroy_cartesian_metrics(metrics)
      !$acc exit data delete(force_u, force_v)

      max_diff_eta = 0.0_wp
      do j = 1, ny
         do i = 1, nx
            max_diff_eta = max(max_diff_eta, &
                               abs(bt_off%bt_eta(i, j) - bt_on%bt_eta(i, j)))
         end do
      end do
      call check(error, max_diff_eta < 1.0e-14_wp, &
                 "bt_eta: upstream-on must match off on uniform-h")
      if (allocated(error)) return
      max_diff_ubt = 0.0_wp
      do j = 1, ny
         do i = 1, nx + 1
            max_diff_ubt = max(max_diff_ubt, &
                               abs(bt_off%bt_ubt(i, j) - bt_on%bt_ubt(i, j)))
         end do
      end do
      call check(error, max_diff_ubt < 1.0e-14_wp, &
                 "bt_ubt: upstream-on must match off on uniform-h")

      deallocate (force_u, force_v)
      call bt_off%destroy()
      call bt_on%destroy()
   end subroutine test_substep_linear_uniform_match

   subroutine test_apply_bt_corr_depth_mean(error)
      !! Uniform-Δu corrector + upstream-h depth-mean projection.
      !! Default (non-h-weighted, non-visc_rem) corrector adds
      !! `delta_bar = bt_ubt_end - bt_ubt_at_n - dt·F_bt_u` uniformly
      !! to every layer.  With `bt_ubt_at_n = 0` and initial u-face
      !! seeded to ZERO, the post-correction upstream-weighted depth
      !! mean equals exactly `delta_bar = bt_ubt_end`.
      !!
      !! Setup: spoon margin (west bed=0/surf=1000, east bed=30/
      !! surf=1000).  Initial u_face = 0 on both layers (only the
      !! sign matters for the upstream pick, so we seed tiny
      !! distinct signs in each layer).  Target bt_ubt_end = 0.10.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(barotropic_workstate_t) :: bt_work
      type(ocean_metrics_t) :: metrics
      real(wp), parameter :: bt_ubt_target = 0.10_wp
      real(wp) :: u_bed, u_surf, h_up_bed, h_up_surf, depth_mean

      ! Seed u_face signs only (magnitudes 1e-12 so initial depth
      ! mean is ~0): bed u>0 (upstream=west, h=0), surf u<0
      ! (upstream=east, h=1000).
      call make_setup(grid, ms, bt_work, 2, &
                      [0.0_wp, 1000.0_wp], &
                      [30.0_wp, 1000.0_wp], &
                      [1.0e-12_wp, -1.0e-12_wp], &
                      enable_upstream=.true.)
      call compute_h_face_upstream(grid, bt_work, ms)
      ! Set the corrector target.  bt_ubt_at_n = 0 + dt·F_bt_u = 0
      ! ⇒ delta_bar = bt_ubt_end.
      bt_work%bt_ubt_end(2, 1) = bt_ubt_target
      ! `metrics` is REQUIRED; all-open (`use_closed_faces = .false.`)
      ! selects the full-column fold this test asserts.
      call make_cartesian_metrics(metrics, grid)
      call apply_bt_correction(bt_work, ms, 1.0_wp, metrics, &
                               skip_h_rescale=.true.)
      call destroy_cartesian_metrics(metrics)

      u_bed = ms%u_face_x_layer(2, 1, 1)
      u_surf = ms%u_face_x_layer(2, 1, 2)
      ! Both u signs still match the original setup (added a uniform
      ! +0.10 to ~0); upstream picks unchanged.
      h_up_bed = 0.0_wp     ! west bed = 0
      h_up_surf = 1000.0_wp ! east surf = 1000
      depth_mean = (u_bed*h_up_bed + u_surf*h_up_surf)/(h_up_bed + h_up_surf)
      call check(error, abs(depth_mean - bt_ubt_target) < 1.0e-10_wp, &
                 "post-correction upstream-weighted depth-mean ≠ bt_ubt_end")

      call ms%destroy()
      call bt_work%destroy()
   end subroutine test_apply_bt_corr_depth_mean

end module test_ocean_bt_upstream_h_face
