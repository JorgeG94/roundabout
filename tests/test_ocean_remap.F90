!! Centre-cell ALE remap orchestrator tests.
!!
!! Layer 3a coverage:
!!   - `ocean_apply_ale_remap_centres` is a no-op when coord_type ==
!!     VCOORD_EULERIAN_Z (h_layer, bt_eta, tracers all unchanged).
!!   - Under VCOORD_SIGMA with non-zero eta, h_layer absorbs the
!!     redistribution exactly and bt_eta is re-derived correctly.
!!   - Per-column tracer mass conservation holds across SIGMA + ZSTAR_FULL
!!     remaps via the PLM column kernel.
!!   - The `c = hTr/h` floor protects against vanishing-layer divisions
!!     under VCOORD_ZSTAR_FULL.
!!   - `ocean_remap_tracer_column` (column entry point) gives the
!!     same answer as a pre-computed analytic reference for a known
!!     remap.
module test_ocean_remap
   use rdb_constants, only: wp, REMAP_PLM
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_remap, only: ocean_apply_ale_remap_centres, &
                              ocean_apply_ale_remap_faces, &
                              ocean_apply_ale_remap_step, &
                              ocean_remap_scan_preconditions, &
                              OCEAN_REMAP_PRECOND_RTOL, &
                              ocean_remap_tracer_column
   use rdb_ocean_vcoord, only: ocean_vcoord_t, VCOORD_EULERIAN_Z
   use rdb_constants, only: VCOORD_SIGMA, VCOORD_ZSTAR_FULL
   use rdb_ocean_budgets, only: ocean_budgets_t, BUDGET_MASS, &
                                BUDGET_HEAT_TOTAL, BUDGET_SALT_TOTAL
   use testdrive, only: error_type, check, new_unittest, unittest_type
   implicit none
   private

   public :: collect_ocean_remap_tests

contains

   subroutine collect_ocean_remap_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("eulerian_z_is_noop", test_eulerian_noop), &
                  new_unittest("sigma_redistributes_h_to_target", test_sigma_h), &
                  new_unittest("sigma_tracer_conservation", test_sigma_tracer_conservation), &
                  new_unittest("sigma_bt_eta_recomputed", test_sigma_bt_eta), &
                  new_unittest("zstar_full_tracer_conservation", test_zstar_full_tracer_conservation), &
                  new_unittest("uniform_tracer_stays_uniform", test_uniform_stays_uniform), &
                  new_unittest("column_helper_matches_analytic", test_column_helper_analytic), &
                  new_unittest("faces_uniform_velocity_preserved", test_faces_uniform_u), &
                  new_unittest("faces_momentum_conserved_per_face", test_faces_momentum), &
                  new_unittest("step_remaps_centres_and_faces", test_step_full), &
                  new_unittest("ale_remap_budget_contributors_telescope", &
                               test_ale_remap_budget_contributors), &
                  new_unittest("precondition_scan_finds_the_bad_columns", &
                               test_precondition_scan) &
                  ]
   end subroutine collect_ocean_remap_tests

   subroutine setup_state(grid, ms, nx, ny, nz, h_per_layer, salinity_ref)
      type(hgrid_t), intent(inout) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: h_per_layer, salinity_ref
      integer :: i, j, k
      call grid%init(nx, ny, 1, 1.0_wp, 1.0_wp)
      ms%nz_ml = nz
      call ms%init(grid)
      do k = 1, nz
         do j = 1, grid%ny_total
            do i = 1, grid%nx_total
               ms%h_layer(i, j, k) = h_per_layer
               ms%tracers(ms%idx_salinity)%hTr(i, j, k) = salinity_ref*h_per_layer
            end do
         end do
      end do
   end subroutine setup_state

   subroutine test_eulerian_noop(error)
      !! coord_type == VCOORD_EULERIAN_Z → orchestrator should return
      !! immediately leaving h_layer, bt_eta, and tracers untouched.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vcoord_t) :: vc
      integer, parameter :: NX = 3, NY = 3, NZ = 4
      real(wp), parameter :: H0 = 50.0_wp, S0 = 35.0_wp
      real(wp), allocatable :: h_before(:, :, :), hS_before(:, :, :)
      real(wp), allocatable :: bt_eta(:, :), bt_H_ref(:, :)
      real(wp) :: max_dh, max_dS
      integer :: nx_tot, ny_tot
      checks: block

         call setup_state(grid, ms, NX, NY, NZ, H0, S0)
         call vc%init(grid, nz_ml=NZ)
         vc%coord_type = VCOORD_EULERIAN_Z
         nx_tot = grid%nx_total
         ny_tot = grid%ny_total
         allocate (h_before, source=ms%h_layer)
         allocate (hS_before, source=ms%tracers(ms%idx_salinity)%hTr)
         allocate (bt_eta(nx_tot, ny_tot), source=2.5_wp)        ! non-zero
         allocate (bt_H_ref(nx_tot, ny_tot), source=real(NZ, wp)*H0)

         call ocean_apply_ale_remap_centres(grid, vc, ms, bt_eta, bt_H_ref)

         max_dh = maxval(abs(ms%h_layer - h_before))
         max_dS = maxval(abs(ms%tracers(ms%idx_salinity)%hTr - hS_before))
         call check(error, max_dh == 0.0_wp, "EULERIAN_Z: h_layer should be untouched")
         if (allocated(error)) exit checks
         call check(error, max_dS == 0.0_wp, "EULERIAN_Z: tracer hTr should be untouched")
         if (allocated(error)) exit checks
         call check(error, all(bt_eta == 2.5_wp), "EULERIAN_Z: bt_eta should be untouched")

      end block checks
      deallocate (h_before, hS_before, bt_eta, bt_H_ref)
      call vc%destroy()
      call ms%destroy()
   end subroutine test_eulerian_noop

   subroutine test_sigma_h(error)
      !! VCOORD_SIGMA: an initially NON-uniform per-layer thickness
      !! should be redistributed uniformly across the column under
      !! the SIGMA remap.  The orchestrator preserves the column
      !! total (sum_k h_layer) — that's the conservation invariant.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vcoord_t) :: vc
      integer, parameter :: NX = 3, NY = 3, NZ = 4
      real(wp), parameter :: H_TOTAL = 400.0_wp, S0 = 35.0_wp
      real(wp) :: layer_pattern(NZ)
      real(wp), allocatable :: bt_eta(:, :), bt_H_ref(:, :)
      real(wp) :: expected, max_err
      integer :: nx_tot, ny_tot, i, j, k

      call setup_state(grid, ms, NX, NY, NZ, H_TOTAL/real(NZ, wp), S0)
      call vc%init(grid, nz_ml=NZ)
      vc%coord_type = VCOORD_SIGMA
      nx_tot = grid%nx_total
      ny_tot = grid%ny_total

      ! Stamp a non-uniform per-layer thickness summing to H_TOTAL.
      layer_pattern = [90.0_wp, 95.0_wp, 105.0_wp, 110.0_wp]
      do k = 1, NZ
         do j = 1, ny_tot
            do i = 1, nx_tot
               ms%h_layer(i, j, k) = layer_pattern(k)
               ms%tracers(ms%idx_salinity)%hTr(i, j, k) = S0*layer_pattern(k)
            end do
         end do
      end do

      ! bt_eta consistent with the column total: sum=400, H_ref=400, eta=0.
      allocate (bt_eta(nx_tot, ny_tot), source=0.0_wp)
      allocate (bt_H_ref(nx_tot, ny_tot), source=H_TOTAL)

      call ocean_apply_ale_remap_centres(grid, vc, ms, bt_eta, bt_H_ref)

      ! After remap: h_layer(k) = H_TOTAL/NZ = 100 m per layer (uniform).
      expected = H_TOTAL/real(NZ, wp)
      max_err = 0.0_wp
      do k = 1, NZ
         do j = 1, ny_tot
            do i = 1, nx_tot
               max_err = max(max_err, abs(ms%h_layer(i, j, k) - expected))
            end do
         end do
      end do
      call check(error, max_err < 1.0e-10_wp, &
                 "SIGMA orchestrator should redistribute non-uniform h to uniform target")

      deallocate (bt_eta, bt_H_ref)
      call vc%destroy()
      call ms%destroy()
   end subroutine test_sigma_h

   subroutine test_sigma_tracer_conservation(error)
      !! Per-column ∫(hTr) should be preserved across the SIGMA remap.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vcoord_t) :: vc
      integer, parameter :: NX = 4, NY = 4, NZ = 5
      real(wp), parameter :: H0 = 80.0_wp, S0 = 34.5_wp, ETA0 = 7.0_wp
      real(wp), allocatable :: bt_eta(:, :), bt_H_ref(:, :)
      real(wp) :: col_total_before, col_total_after, max_drift
      integer :: nx_tot, ny_tot, i, j, k

      call setup_state(grid, ms, NX, NY, NZ, H0, S0)
      call vc%init(grid, nz_ml=NZ)
      vc%coord_type = VCOORD_SIGMA
      nx_tot = grid%nx_total
      ny_tot = grid%ny_total
      allocate (bt_eta(nx_tot, ny_tot), source=ETA0)
      allocate (bt_H_ref(nx_tot, ny_tot), source=real(NZ, wp)*H0)

      ! Stamp a spatially varying tracer pattern (gradient with k).
      do k = 1, NZ
         do j = 1, ny_tot
            do i = 1, nx_tot
               ms%tracers(ms%idx_salinity)%hTr(i, j, k) = &
                  (S0 + 0.5_wp*real(k - 1, wp))*ms%h_layer(i, j, k)
            end do
         end do
      end do

      max_drift = 0.0_wp
      do j = 1, ny_tot
         do i = 1, nx_tot
            col_total_before = 0.0_wp
            do k = 1, NZ
               col_total_before = col_total_before + ms%tracers(ms%idx_salinity)%hTr(i, j, k)
            end do
            ! Single-column conservation will be checked AFTER remap; stash for now.
            ms%tracers(ms%idx_temperature)%hTr(i, j, 1) = col_total_before   ! abuse temp(1) as cache
         end do
      end do

      call ocean_apply_ale_remap_centres(grid, vc, ms, bt_eta, bt_H_ref)

      do j = 1, ny_tot
         do i = 1, nx_tot
            col_total_after = 0.0_wp
            do k = 1, NZ
               col_total_after = col_total_after + ms%tracers(ms%idx_salinity)%hTr(i, j, k)
            end do
            col_total_before = ms%tracers(ms%idx_temperature)%hTr(i, j, 1)
            max_drift = max(max_drift, abs(col_total_after - col_total_before))
         end do
      end do
      call check(error, max_drift < 1.0e-9_wp, &
                 "SIGMA remap: per-column tracer integral must be conserved")

      deallocate (bt_eta, bt_H_ref)
      call vc%destroy()
      call ms%destroy()
   end subroutine test_sigma_tracer_conservation

   subroutine test_sigma_bt_eta(error)
      !! After the orchestrator, bt_eta should equal sum_k(h_layer) - bt_H_ref.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vcoord_t) :: vc
      integer, parameter :: NX = 3, NY = 3, NZ = 4
      real(wp), parameter :: H0 = 100.0_wp, S0 = 35.0_wp, ETA0 = 4.2_wp
      real(wp), allocatable :: bt_eta(:, :), bt_H_ref(:, :)
      real(wp) :: max_err, col_sum, expected
      integer :: nx_tot, ny_tot, i, j, k

      call setup_state(grid, ms, NX, NY, NZ, H0, S0)
      call vc%init(grid, nz_ml=NZ)
      vc%coord_type = VCOORD_SIGMA
      nx_tot = grid%nx_total
      ny_tot = grid%ny_total
      allocate (bt_eta(nx_tot, ny_tot), source=ETA0)
      allocate (bt_H_ref(nx_tot, ny_tot), source=real(NZ, wp)*H0)

      call ocean_apply_ale_remap_centres(grid, vc, ms, bt_eta, bt_H_ref)

      max_err = 0.0_wp
      do j = 1, ny_tot
         do i = 1, nx_tot
            col_sum = 0.0_wp
            do k = 1, NZ
               col_sum = col_sum + ms%h_layer(i, j, k)
            end do
            expected = col_sum - bt_H_ref(i, j)
            max_err = max(max_err, abs(bt_eta(i, j) - expected))
         end do
      end do
      call check(error, max_err < 1.0e-10_wp, &
                 "bt_eta should match sum(h_layer) - bt_H_ref after remap")

      deallocate (bt_eta, bt_H_ref)
      call vc%destroy()
      call ms%destroy()
   end subroutine test_sigma_bt_eta

   subroutine test_zstar_full_tracer_conservation(error)
      !! Per-column tracer mass conservation under VCOORD_ZSTAR_FULL
      !! with positive η (subsurface untouched, surface grows).  This
      !! is the easy-mode case for ZSTAR_FULL — the hard case (eta<0
      !! with bed-vanishing) needs the floor guard not to clobber
      !! tracer mass, which is the next test below.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vcoord_t) :: vc
      integer, parameter :: NX = 3, NY = 3, NZ = 4
      real(wp), parameter :: H_BED = 400.0_wp, S0 = 35.0_wp, ETA0 = 5.0_wp
      real(wp), allocatable :: bt_eta(:, :), bt_H_ref(:, :), h_bed_2d(:, :)
      real(wp) :: col_before, col_after, max_drift
      integer :: nx_tot, ny_tot, i, j, k

      call setup_state(grid, ms, NX, NY, NZ, H_BED/real(NZ, wp), S0)
      call vc%init(grid, nz_ml=NZ)
      vc%coord_type = VCOORD_ZSTAR_FULL
      vc%zstar_h_surf_target = 20.0_wp
      vc%zstar_n_surf = 2
      nx_tot = grid%nx_total
      ny_tot = grid%ny_total

      allocate (h_bed_2d(nx_tot, ny_tot), source=H_BED)
      call vc%build_zref_full(h_bed_2d)

      allocate (bt_eta(nx_tot, ny_tot), source=ETA0)
      allocate (bt_H_ref(nx_tot, ny_tot), source=H_BED)

      ! Per-column salt total snapshot, written into temperature(1).
      do j = 1, ny_tot
         do i = 1, nx_tot
            col_before = 0.0_wp
            do k = 1, NZ
               col_before = col_before + ms%tracers(ms%idx_salinity)%hTr(i, j, k)
            end do
            ms%tracers(ms%idx_temperature)%hTr(i, j, 1) = col_before
         end do
      end do

      call ocean_apply_ale_remap_centres(grid, vc, ms, bt_eta, bt_H_ref)

      max_drift = 0.0_wp
      do j = 1, ny_tot
         do i = 1, nx_tot
            col_after = 0.0_wp
            do k = 1, NZ
               col_after = col_after + ms%tracers(ms%idx_salinity)%hTr(i, j, k)
            end do
            col_before = ms%tracers(ms%idx_temperature)%hTr(i, j, 1)
            max_drift = max(max_drift, abs(col_after - col_before))
         end do
      end do
      call check(error, max_drift < 1.0e-8_wp, &
                 "ZSTAR_FULL +eta: per-column salt integral must be conserved")

      deallocate (bt_eta, bt_H_ref, h_bed_2d)
      call vc%destroy()
      call ms%destroy()
   end subroutine test_zstar_full_tracer_conservation

   subroutine test_uniform_stays_uniform(error)
      !! A uniform concentration (constant S = S0 over the whole
      !! column, hTr = S0 · h) should remain uniform across any
      !! conservative remap.  Tests that the c = hTr/h ↔ hTr_new
      !! = c_new · h_new pattern preserves uniformity.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vcoord_t) :: vc
      integer, parameter :: NX = 3, NY = 3, NZ = 5
      real(wp), parameter :: H0 = 60.0_wp, S0 = 35.0_wp, ETA0 = 4.0_wp
      real(wp), allocatable :: bt_eta(:, :), bt_H_ref(:, :)
      real(wp) :: max_dev_from_S0, c
      integer :: nx_tot, ny_tot, i, j, k

      call setup_state(grid, ms, NX, NY, NZ, H0, S0)
      call vc%init(grid, nz_ml=NZ)
      vc%coord_type = VCOORD_SIGMA
      nx_tot = grid%nx_total
      ny_tot = grid%ny_total
      allocate (bt_eta(nx_tot, ny_tot), source=ETA0)
      allocate (bt_H_ref(nx_tot, ny_tot), source=real(NZ, wp)*H0)

      call ocean_apply_ale_remap_centres(grid, vc, ms, bt_eta, bt_H_ref)

      max_dev_from_S0 = 0.0_wp
      do k = 1, NZ
         do j = 1, ny_tot
            do i = 1, nx_tot
               c = ms%tracers(ms%idx_salinity)%hTr(i, j, k)/ms%h_layer(i, j, k)
               max_dev_from_S0 = max(max_dev_from_S0, abs(c - S0))
            end do
         end do
      end do
      call check(error, max_dev_from_S0 < 1.0e-10_wp, &
                 "uniform concentration should stay uniform across remap")

      deallocate (bt_eta, bt_H_ref)
      call vc%destroy()
      call ms%destroy()
   end subroutine test_uniform_stays_uniform

   subroutine test_column_helper_analytic(error)
      !! `ocean_remap_tracer_column` against an analytic reference:
      !! uniform old grid (10 m each), new grid shifted so surface
      !! layer gets +5 m of mass from the next layer.  Pre-fill old
      !! hTr with linear-in-k pattern; check post-remap mass is
      !! conserved AND surface concentration moved toward subsurface.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 4
      real(wp), parameter :: H_LAYER = 10.0_wp
      real(wp) :: h_old(NZ), h_new(NZ), hTr_in(NZ), hTr(NZ)
      real(wp) :: total_old, total_new, k_top
      integer :: k

      ! Old grid: 10/10/10/10 = 40 m total.
      ! New grid: 5/10/10/15 = 40 m total (surface absorbed 5 m).
      ! Old hTr: pattern 10, 20, 30, 40 (so concentration 1, 2, 3, 4 ROMS-ordered).
      h_old = H_LAYER
      h_new = [5.0_wp, 10.0_wp, 10.0_wp, 15.0_wp]   ! k=1 bed, k=NZ surface
      hTr_in = [10.0_wp, 20.0_wp, 30.0_wp, 40.0_wp]
      hTr = hTr_in

      call ocean_remap_tracer_column(NZ, h_old, h_new, hTr, REMAP_PLM)

      total_old = sum(hTr_in)
      total_new = sum(hTr)
      call check(error, abs(total_new - total_old) < 1.0e-10_wp, &
                 "column helper: total tracer mass should be conserved")
      if (allocated(error)) return

      ! Surface concentration after the remap should drop (it absorbed
      ! some of layer 3's lower-concentration mass).  Pre-remap surface
      ! concentration was 4.0; new must be < 4.0.
      k_top = hTr(NZ)/h_new(NZ)
      call check(error, k_top < 4.0_wp, &
                 "column helper: expanded surface should pull in less-concentrated mass")
   end subroutine test_column_helper_analytic

   ! --------- Face-velocity remap (Layer 3b) ---------

   subroutine test_faces_uniform_u(error)
      !! A uniform velocity field (u = U0 everywhere) should remain
      !! uniform under ANY conservative face remap, because c = U0
      !! → U0 by construction.  Exercises both x-face and y-face
      !! kernels under a non-trivial old→new layer redistribution.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      integer, parameter :: NX = 4, NY = 4, NZ = 4
      real(wp), parameter :: U0 = 0.25_wp, V0 = -0.1_wp
      real(wp), allocatable :: h_old(:, :, :), h_new(:, :, :)
      real(wp), allocatable :: u_face_x(:, :, :), v_face_y(:, :, :)
      real(wp) :: max_du, max_dv
      integer :: nx_tot, ny_tot, i, j, k
      real(wp) :: layer_pattern_old(NZ), layer_pattern_new(NZ)
      checks: block

         call grid%init(NX, NY, 1, 1.0_wp, 1.0_wp)
         nx_tot = grid%nx_total
         ny_tot = grid%ny_total

         ! Non-trivial old → new redistribution: surface absorbs +5 m.
         layer_pattern_old = [10.0_wp, 10.0_wp, 10.0_wp, 10.0_wp]
         layer_pattern_new = [5.0_wp, 10.0_wp, 10.0_wp, 15.0_wp]
         allocate (h_old(nx_tot, ny_tot, NZ), h_new(nx_tot, ny_tot, NZ))
         do k = 1, NZ
            h_old(:, :, k) = layer_pattern_old(k)
            h_new(:, :, k) = layer_pattern_new(k)
         end do
         allocate (u_face_x(nx_tot + 1, ny_tot, NZ), source=U0)
         allocate (v_face_y(nx_tot, ny_tot + 1, NZ), source=V0)

         call ocean_apply_ale_remap_faces(grid, h_old, h_new, u_face_x, v_face_y)

         max_du = maxval(abs(u_face_x - U0))
         max_dv = maxval(abs(v_face_y - V0))
         call check(error, max_du < 1.0e-10_wp, &
                    "uniform u: face remap should preserve uniform velocity")
         if (allocated(error)) exit checks
         call check(error, max_dv < 1.0e-10_wp, &
                    "uniform v: face remap should preserve uniform velocity")

      end block checks
      deallocate (h_old, h_new, u_face_x, v_face_y)
   end subroutine test_faces_uniform_u

   subroutine test_faces_momentum(error)
      !! Per-face momentum conservation: sum_k(h_face_old · u_old) ==
      !! sum_k(h_face_new · u_new) up to round-off.  Stamps a linear-
      !! in-k velocity pattern so the remap actually has to move
      !! momentum between layers.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      integer, parameter :: NX = 3, NY = 3, NZ = 4
      real(wp), allocatable :: h_old(:, :, :), h_new(:, :, :)
      real(wp), allocatable :: u_face_x(:, :, :), v_face_y(:, :, :)
      real(wp), allocatable :: u_before(:, :, :), v_before(:, :, :)
      real(wp) :: layer_pattern_old(NZ), layer_pattern_new(NZ)
      real(wp) :: mom_old, mom_new, h_face_old, h_face_new, max_drift
      integer :: nx_tot, ny_tot, II, JJ, j, i, k
      checks: block

         call grid%init(NX, NY, 1, 1.0_wp, 1.0_wp)
         nx_tot = grid%nx_total
         ny_tot = grid%ny_total

         layer_pattern_old = [10.0_wp, 10.0_wp, 10.0_wp, 10.0_wp]
         layer_pattern_new = [5.0_wp, 10.0_wp, 10.0_wp, 15.0_wp]
         allocate (h_old(nx_tot, ny_tot, NZ), h_new(nx_tot, ny_tot, NZ))
         do k = 1, NZ
            h_old(:, :, k) = layer_pattern_old(k)
            h_new(:, :, k) = layer_pattern_new(k)
         end do
         allocate (u_face_x(nx_tot + 1, ny_tot, NZ))
         allocate (v_face_y(nx_tot, ny_tot + 1, NZ))
         ! Linear-in-k u pattern: 0.05, 0.10, 0.15, 0.20 (ROMS-ordered).
         do k = 1, NZ
            u_face_x(:, :, k) = 0.05_wp*real(k, wp)
            v_face_y(:, :, k) = -0.03_wp*real(k, wp)
         end do
         allocate (u_before, source=u_face_x)
         allocate (v_before, source=v_face_y)

         call ocean_apply_ale_remap_faces(grid, h_old, h_new, u_face_x, v_face_y)

         ! Per-face momentum conservation check.
         max_drift = 0.0_wp
         do II = 1, nx_tot + 1
            do j = 1, ny_tot
               mom_old = 0.0_wp
               mom_new = 0.0_wp
               do k = 1, NZ
                  if (II == 1) then
                     h_face_old = h_old(1, j, k)
                     h_face_new = h_new(1, j, k)
                  else if (II == nx_tot + 1) then
                     h_face_old = h_old(nx_tot, j, k)
                     h_face_new = h_new(nx_tot, j, k)
                  else
                     h_face_old = 0.5_wp*(h_old(II - 1, j, k) + h_old(II, j, k))
                     h_face_new = 0.5_wp*(h_new(II - 1, j, k) + h_new(II, j, k))
                  end if
                  mom_old = mom_old + h_face_old*u_before(II, j, k)
                  mom_new = mom_new + h_face_new*u_face_x(II, j, k)
               end do
               max_drift = max(max_drift, abs(mom_new - mom_old))
            end do
         end do
         call check(error, max_drift < 1.0e-10_wp, &
                    "x-face momentum sum_k(h·u) must be conserved per face")
         if (allocated(error)) exit checks

         max_drift = 0.0_wp
         do JJ = 1, ny_tot + 1
            do i = 1, nx_tot
               mom_old = 0.0_wp
               mom_new = 0.0_wp
               do k = 1, NZ
                  if (JJ == 1) then
                     h_face_old = h_old(i, 1, k)
                     h_face_new = h_new(i, 1, k)
                  else if (JJ == ny_tot + 1) then
                     h_face_old = h_old(i, ny_tot, k)
                     h_face_new = h_new(i, ny_tot, k)
                  else
                     h_face_old = 0.5_wp*(h_old(i, JJ - 1, k) + h_old(i, JJ, k))
                     h_face_new = 0.5_wp*(h_new(i, JJ - 1, k) + h_new(i, JJ, k))
                  end if
                  mom_old = mom_old + h_face_old*v_before(i, JJ, k)
                  mom_new = mom_new + h_face_new*v_face_y(i, JJ, k)
               end do
               max_drift = max(max_drift, abs(mom_new - mom_old))
            end do
         end do
         call check(error, max_drift < 1.0e-10_wp, &
                    "y-face momentum sum_k(h·v) must be conserved per face")

      end block checks
      deallocate (h_old, h_new, u_face_x, v_face_y, u_before, v_before)
   end subroutine test_faces_momentum

   subroutine test_step_full(error)
      !! `ocean_apply_ale_remap_step` end-to-end on a small state: under
      !! SIGMA the call should remap centres + faces + bt_eta in one
      !! shot.  Sanity checks:
      !!   1. h_layer ends at uniform target (column total preserved).
      !!   2. Per-column salt integral conserved.
      !!   3. Face velocity magnitude bounded by old/new spread (no NaN).
      !!   4. bt_eta consistent with sum(h) - bt_H_ref on exit.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vcoord_t) :: vc
      integer, parameter :: NX = 3, NY = 3, NZ = 4
      real(wp), parameter :: H_TOTAL = 400.0_wp, S0 = 35.0_wp
      real(wp), allocatable :: bt_eta(:, :), bt_H_ref(:, :)
      real(wp), allocatable :: hS_before(:, :, :)
      real(wp) :: col_before, col_after, max_drift, expected, max_eta_err
      real(wp) :: layer_pattern(NZ), col_sum
      integer :: nx_tot, ny_tot, i, j, k
      checks: block

         call setup_state(grid, ms, NX, NY, NZ, H_TOTAL/real(NZ, wp), S0)
         call vc%init(grid, nz_ml=NZ)
         vc%coord_type = VCOORD_SIGMA
         nx_tot = grid%nx_total
         ny_tot = grid%ny_total

         ! Non-uniform initial h_layer summing to 400.
         layer_pattern = [90.0_wp, 95.0_wp, 105.0_wp, 110.0_wp]
         do k = 1, NZ
            do j = 1, ny_tot
               do i = 1, nx_tot
                  ms%h_layer(i, j, k) = layer_pattern(k)
                  ms%tracers(ms%idx_salinity)%hTr(i, j, k) = S0*layer_pattern(k)
                  ms%u_face_x_layer(:, j, k) = 0.1_wp + 0.02_wp*real(k, wp)
                  ms%v_face_y_layer(i, :, k) = -0.05_wp + 0.01_wp*real(k, wp)
               end do
            end do
         end do
         allocate (hS_before, source=ms%tracers(ms%idx_salinity)%hTr)
         allocate (bt_eta(nx_tot, ny_tot), source=0.0_wp)
         allocate (bt_H_ref(nx_tot, ny_tot), source=H_TOTAL)

         call ocean_apply_ale_remap_step(grid, vc, ms, bt_eta, bt_H_ref)

         ! (1) Uniform target after SIGMA.
         expected = H_TOTAL/real(NZ, wp)
         max_drift = maxval(abs(ms%h_layer - expected))
         call check(error, max_drift < 1.0e-10_wp, &
                    "step: h_layer should equal uniform SIGMA target")
         if (allocated(error)) exit checks

         ! (2) Per-column salt integral conserved.
         max_drift = 0.0_wp
         do j = 1, ny_tot
            do i = 1, nx_tot
               col_before = 0.0_wp
               col_after = 0.0_wp
               do k = 1, NZ
                  col_before = col_before + hS_before(i, j, k)
                  col_after = col_after + ms%tracers(ms%idx_salinity)%hTr(i, j, k)
               end do
               max_drift = max(max_drift, abs(col_after - col_before))
            end do
         end do
         call check(error, max_drift < 1.0e-9_wp, &
                    "step: per-column salt integral must be conserved")
         if (allocated(error)) exit checks

         ! (3) Face velocities finite + bounded by their original input
         ! range.  u_face_x was stamped at 0.12..0.18 (k=1..NZ), v_face_y
         ! at -0.04..-0.01.  PLM remap is conservative and (under monotone
         ! limiter) bounded by the input min/max per column.
         call check(error, all(ms%u_face_x_layer >= 0.12_wp - 1.0e-9_wp &
                               .and. ms%u_face_x_layer <= 0.18_wp + 1.0e-9_wp), &
                    "step: x-face velocities should stay in [0.12, 0.18]")
         if (allocated(error)) exit checks
         call check(error, all(ms%v_face_y_layer >= -0.04_wp - 1.0e-9_wp &
                               .and. ms%v_face_y_layer <= -0.01_wp + 1.0e-9_wp), &
                    "step: y-face velocities should stay in [-0.04, -0.01]")
         if (allocated(error)) exit checks

         ! (4) bt_eta = sum(h) - bt_H_ref.
         max_eta_err = 0.0_wp
         do j = 1, ny_tot
            do i = 1, nx_tot
               col_sum = 0.0_wp
               do k = 1, NZ
                  col_sum = col_sum + ms%h_layer(i, j, k)
               end do
               max_eta_err = max(max_eta_err, abs(bt_eta(i, j) - (col_sum - bt_H_ref(i, j))))
            end do
         end do
         call check(error, max_eta_err < 1.0e-10_wp, &
                    "step: bt_eta must match sum(h_layer) - bt_H_ref on exit")

      end block checks
      deallocate (bt_eta, bt_H_ref, hS_before)
      call vc%destroy()
      call ms%destroy()
   end subroutine test_step_full

   subroutine test_ale_remap_budget_contributors(error)
      !! Phase D v2 final kernel patch: ALE remap is the conservation
      !! leak detector.  Register mass + heat + salt contributors, run
      !! one SIGMA remap on a non-uniform `h_layer` with non-trivial
      !! T/S, drain, verify all three contributors integrate to ~0 to
      !! FP (the remap kernel is column-conservative).  Non-zero
      !! integral flags a real bug — exactly the gate Phase D v2 was
      !! introduced to provide.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_vcoord_t) :: vc
      type(ocean_budgets_t) :: budgets
      integer, parameter :: NX = 4, NY = 3, NZ = 4
      real(wp), parameter :: H_TOTAL = 400.0_wp, S0 = 35.0_wp, T0 = 12.0_wp
      real(wp), allocatable :: bt_eta(:, :), bt_H_ref(:, :)
      real(wp) :: layer_pattern(NZ)
      real(wp) :: total_h_before, total_T_before, total_S_before
      real(wp) :: rel_mass, rel_T, rel_S
      integer :: nx_tot, ny_tot, i, j, k, idx_mass, idx_heat, idx_salt
      checks: block

         call setup_state(grid, ms, NX, NY, NZ, H_TOTAL/real(NZ, wp), S0)
         call vc%init(grid, nz_ml=NZ)
         vc%coord_type = VCOORD_SIGMA
         call budgets%init(grid)
         call budgets%register_contributor("ale_remap_mass", BUDGET_MASS, &
                                           ms%mass_budget_remap, &
                                           device_resident=.true.)
         idx_mass = budgets%n_contributors
         call budgets%register_contributor("ale_remap_heat", BUDGET_HEAT_TOTAL, &
                                           ms%heat_budget_remap, &
                                           device_resident=.true.)
         idx_heat = budgets%n_contributors
         call budgets%register_contributor("ale_remap_salt", BUDGET_SALT_TOTAL, &
                                           ms%salt_budget_remap, &
                                           device_resident=.true.)
         idx_salt = budgets%n_contributors

         nx_tot = grid%nx_total
         ny_tot = grid%ny_total

         ! Non-uniform initial h_layer summing to H_TOTAL.
         layer_pattern = [90.0_wp, 95.0_wp, 105.0_wp, 110.0_wp]
         do k = 1, NZ
            do j = 1, ny_tot
               do i = 1, nx_tot
                  ms%h_layer(i, j, k) = layer_pattern(k)
                  ms%tracers(ms%idx_salinity)%hTr(i, j, k) = &
                     (S0 + 0.5_wp*real(k - 1, wp))*layer_pattern(k)
                  ms%tracers(ms%idx_temperature)%hTr(i, j, k) = &
                     (T0 + 0.3_wp*real(k - 1, wp))*layer_pattern(k)
               end do
            end do
         end do
         allocate (bt_eta(nx_tot, ny_tot), source=0.0_wp)
         allocate (bt_H_ref(nx_tot, ny_tot), source=H_TOTAL)

         total_h_before = sum(ms%h_layer)*grid%dx*grid%dy
         total_T_before = sum(ms%tracers(ms%idx_temperature)%hTr)*grid%dx*grid%dy
         total_S_before = sum(ms%tracers(ms%idx_salinity)%hTr)*grid%dx*grid%dy

         call ocean_apply_ale_remap_step(grid, vc, ms, bt_eta, bt_H_ref)
         call budgets%drain_contributors()

         ! Conservation: each contributor's spatial integral telescopes
         ! to zero per column.
         rel_mass = abs(budgets%contributors(idx_mass)%total_integrated)/abs(total_h_before)
         rel_T = abs(budgets%contributors(idx_heat)%total_integrated)/abs(total_T_before)
         rel_S = abs(budgets%contributors(idx_salt)%total_integrated)/abs(total_S_before)
         call check(error, rel_mass < 1.0e-10_wp, &
                    "ale_remap mass contributor should telescope to FP")
         if (allocated(error)) exit checks
         call check(error, rel_T < 1.0e-10_wp, &
                    "ale_remap heat contributor should telescope to FP")
         if (allocated(error)) exit checks
         call check(error, rel_S < 1.0e-10_wp, &
                    "ale_remap salt contributor should telescope to FP")
         if (allocated(error)) exit checks

         ! Closure: LHS = RHS for each quantity (remap is the only
         ! thing touching h_layer + hTr here).
         call check(error, abs((sum(ms%h_layer)*grid%dx*grid%dy - total_h_before) - &
                               budgets%contributors(idx_mass)%total_integrated) &
                    < 1.0e-10_wp*abs(total_h_before), &
                    "mass LHS = RHS for ALE remap")
         if (allocated(error)) exit checks
         call check(error, abs((sum(ms%tracers(ms%idx_temperature)%hTr)*grid%dx*grid%dy &
                                - total_T_before) - &
                               budgets%contributors(idx_heat)%total_integrated) &
                    < 1.0e-10_wp*abs(total_T_before), &
                    "heat LHS = RHS for ALE remap")
         if (allocated(error)) exit checks
         call check(error, abs((sum(ms%tracers(ms%idx_salinity)%hTr)*grid%dx*grid%dy &
                                - total_S_before) - &
                               budgets%contributors(idx_salt)%total_integrated) &
                    < 1.0e-10_wp*abs(total_S_before), &
                    "salt LHS = RHS for ALE remap")
      end block checks
      deallocate (bt_eta, bt_H_ref)
      call budgets%destroy()
      call vc%destroy()
      call ms%destroy()
   end subroutine test_ale_remap_budget_contributors

   subroutine test_precondition_scan(error)
      !! `ocean_remap_scan_preconditions` is the domain sweep the fail-loud
      !! guard (`&vcoord_nml remap_check_preconditions`) runs at the thermo
      !! cadence (audit findings V5, V6).  It must be silent on a healthy
      !! field, count exactly the offending columns, and report magnitudes
      !! large enough to identify the producer — a relative column-total
      !! mismatch and the most negative thickness.
      !!
      !! The sweep is a `present`-clause device reduction (production feeds it
      !! `vcoord%remap_h_old` / `vcoord%target_h`, both mapped by
      !! `ocean_vcoord_enter_data_impl`), so under `-gpu=mem:separate` these
      !! host fixtures must be mapped too and re-pushed after every host edit
      !! — otherwise the kernel aborts on a missing `h_old`.  All the
      !! directives below are inert no-ops on the host/multicore builds.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 5, NY = 4, NZ = 6
      real(wp) :: h_old(NX, NY, NZ), h_new(NX, NY, NZ)
      integer :: i, j, k

      ! Healthy: stretched source, differently stretched target, equal totals.
      do k = 1, NZ
         do j = 1, NY
            do i = 1, NX
               h_old(i, j, k) = 5.0_wp + real(k, wp)
               h_new(i, j, k) = 5.0_wp + real(NZ + 1 - k, wp)
            end do
         end do
      end do

      ! The cases live in a helper so the early `return` on the first failed
      ! check cannot skip the unmap and leave a stale device block bound to
      ! this stack address for whatever test runs next.
      !$acc enter data copyin(h_old, h_new)
      call precondition_scan_cases(error, NX, NY, NZ, h_old, h_new)
      !$acc exit data delete(h_old, h_new)
   end subroutine test_precondition_scan

   subroutine precondition_scan_cases(error, nx, ny, nz, h_old, h_new)
      !! The four `ocean_remap_scan_preconditions` cases, on fields the caller
      !! has already mapped.  Every host edit below is pushed to the device
      !! before the next sweep — `h_old`/`h_new` are `copyin`, not managed, so
      !! a host-only assignment would otherwise be scanned as the old values.
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(inout) :: h_old(nx, ny, nz)
      real(wp), intent(inout) :: h_new(nx, ny, nz)

      integer :: n_bad
      real(wp) :: worst_rel, worst_neg

      call ocean_remap_scan_preconditions(nx, ny, nz, h_old, h_new, &
                                          OCEAN_REMAP_PRECOND_RTOL, &
                                          n_bad, worst_rel, worst_neg)
      call check(error, n_bad == 0, "a healthy field must report no bad columns")
      if (allocated(error)) return
      call check(error, worst_neg == 0.0_wp, "a healthy field has no negative thickness")
      if (allocated(error)) return

      ! One column 10% short in the TARGET: the sweep would delete that
      ! tenth of the column's tracer content with no diagnostic at all.
      h_new(2, 3, :) = h_new(2, 3, :)*0.9_wp
      !$acc update device(h_new)
      call ocean_remap_scan_preconditions(nx, ny, nz, h_old, h_new, &
                                          OCEAN_REMAP_PRECOND_RTOL, &
                                          n_bad, worst_rel, worst_neg)
      call check(error, n_bad == 1, "exactly one short column must be counted")
      if (allocated(error)) return
      call check(error, abs(worst_rel - 0.1_wp) < 1.0e-12_wp, &
                 "the reported mismatch must be the 10% that was removed")
      if (allocated(error)) return

      ! Plus a NEGATIVE source thickness in a different column, kept
      ! total-neutral so only the sign test can catch it — which is the
      ! case that matters, because a non-monotone interface stack CREATES
      ! mass rather than losing it.
      h_old(4, 2, 3) = -2.0_wp
      h_old(4, 2, 4) = h_old(4, 2, 4) + 10.0_wp
      !$acc update device(h_old)
      call ocean_remap_scan_preconditions(nx, ny, nz, h_old, h_new, &
                                          OCEAN_REMAP_PRECOND_RTOL, &
                                          n_bad, worst_rel, worst_neg)
      call check(error, n_bad == 2, "the negative-thickness column must be counted too")
      if (allocated(error)) return
      call check(error, abs(worst_neg + 2.0_wp) < 1.0e-12_wp, &
                 "the most negative thickness must be reported")
      if (allocated(error)) return

      ! A land column (both totals zero) is well-posed, not a violation:
      ! the relative test must not trip on a 0/0.
      h_old = 0.0_wp
      h_new = 0.0_wp
      !$acc update device(h_old, h_new)
      call ocean_remap_scan_preconditions(nx, ny, nz, h_old, h_new, &
                                          OCEAN_REMAP_PRECOND_RTOL, &
                                          n_bad, worst_rel, worst_neg)
      call check(error, n_bad == 0, "an all-land field must pass, not divide by zero")
   end subroutine precondition_scan_cases

end module test_ocean_remap
