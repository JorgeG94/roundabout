!! Phase 0e scaffold test for `rdb_ocean_vcoord`.
!!
!! Covers what the shell delivers today:
!!   - `ocean_vcoord_t` init allocates `dsig(nz_ml)` to uniform fractions
!!     that sum to 1.0; default `coord_type == VCOORD_EULERIAN_Z`.
!!   - `destroy` clears `is_init` and releases the allocations.
!!   - `parse_ocean_vcoord_type` round-trips known namelist strings and
!!     falls back to `VCOORD_EULERIAN_Z` on unrecognised input.
!!
!! The per-step `target_h` / `z_ref` paths land in Phase 5g and will get
!! their own tests there.
module test_ocean_vcoord
   use, intrinsic :: iso_fortran_env, only: real64
   use rdb_constants, only: wp, &
                            VCOORD_SIGMA, VCOORD_ZSIGMA, VCOORD_ZSTAR, &
                            VCOORD_ZSTAR_FULL, VCOORD_ZSTAR_SIGMA, &
                            VCOORD_Z_FIXED
   use rdb_grid, only: hgrid_t
   use rdb_ocean_vcoord, only: ocean_vcoord_t, parse_ocean_vcoord_type, &
                               VCOORD_EULERIAN_Z, VCOORD_LAGRANGIAN, &
                               STRETCH_UNIFORM, STRETCH_LOG
   use testdrive, only: error_type, check, new_unittest, unittest_type
   implicit none
   private

   public :: collect_ocean_vcoord_tests

contains

   subroutine collect_ocean_vcoord_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("dsig_uniform_after_init", test_dsig_uniform), &
                  new_unittest("default_coord_is_eulerian_z", test_default_coord), &
                  new_unittest("destroy_releases_allocations", test_destroy), &
                  new_unittest("parse_string_to_code", test_parse), &
                  new_unittest("target_h_eulerian_z", test_target_h_eulerian), &
                  new_unittest("target_h_sigma_conservation", test_target_h_sigma), &
                  new_unittest("target_h_zstar_uses_eta", test_target_h_zstar_eta), &
                  new_unittest("target_h_layer_ratios_match_dsig", test_target_h_ratios), &
                  new_unittest("zsigma_shallow_is_pure_sigma", test_zsigma_shallow), &
                  new_unittest("zsigma_deep_branch_uses_zref", test_zsigma_deep), &
                  new_unittest("zsigma_blend_smoothstep_monotone", test_zsigma_blend_monotone), &
                  new_unittest("zsigma_conservation_per_column", test_zsigma_conservation), &
                  new_unittest("zstar_sigma_shallow_is_pure_sigma", test_zstar_sigma_shallow), &
                  new_unittest("zstar_sigma_deep_is_pure_zstar", test_zstar_sigma_deep), &
                  new_unittest("zstar_sigma_conservation_per_column", test_zstar_sigma_conservation), &
                  new_unittest("zstar_full_build_zref_uniform", test_zstar_full_build_uniform), &
                  new_unittest("zstar_full_build_zref_surface_anchor", test_zstar_full_build_surface), &
                  new_unittest("zstar_full_build_zref_log_monotone", test_zstar_full_build_log), &
                  new_unittest("zstar_full_target_at_reference", test_zstar_full_at_ref), &
                  new_unittest("zstar_full_positive_eta_grows_surface", test_zstar_full_pos_eta), &
                  new_unittest("zstar_full_negative_eta_vanishes_bed", test_zstar_full_neg_eta), &
                  new_unittest("zstar_full_conservation_per_column", test_zstar_full_conservation), &
                  new_unittest("parse_lagrangian_and_isopycnal", test_parse_lagrangian), &
                  new_unittest("target_h_lagrangian_is_no_op", test_target_h_lagrangian_noop), &
                  new_unittest("parse_z_fixed_aliases", test_parse_z_fixed), &
                  new_unittest("z_fixed_deep_locks_interfaces", test_z_fixed_deep), &
                  new_unittest("z_fixed_shallow_vanishes_bed", test_z_fixed_shallow), &
                  new_unittest("z_fixed_unset_falls_back_to_sigma", test_z_fixed_fallback) &
                  ]
   end subroutine collect_ocean_vcoord_tests

   subroutine test_dsig_uniform(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      integer, parameter :: NZ = 8
      real(wp) :: sum_dsig
      integer :: k
      checks: block
         call grid%init(4, 4, 1, 1.0_wp, 1.0_wp)
         call vc%init(grid, nz_ml=NZ)
         call check(error, vc%is_init, "vcoord%is_init should be true after init")
         if (allocated(error)) exit checks
         call check(error, allocated(vc%dsig), "dsig should be allocated")
         if (allocated(error)) exit checks
         call check(error, size(vc%dsig) == NZ, "dsig size should match nz_ml")
         if (allocated(error)) exit checks
         sum_dsig = 0.0_wp
         do k = 1, NZ
            sum_dsig = sum_dsig + vc%dsig(k)
         end do
         call check(error, abs(sum_dsig - 1.0_wp) < 1.0e-12_wp, &
                    "sum(dsig) should equal 1.0")
         if (allocated(error)) exit checks
         do k = 1, NZ
            call check(error, abs(vc%dsig(k) - 1.0_wp/real(NZ, wp)) < 1.0e-12_wp, &
                       "dsig should be uniform 1/nz_ml")
            if (allocated(error)) exit checks
         end do
      end block checks
      call vc%destroy()
   end subroutine test_dsig_uniform

   subroutine test_default_coord(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      call grid%init(2, 2, 1, 1.0_wp, 1.0_wp)
      call vc%init(grid, nz_ml=3)
      call check(error, vc%coord_type == VCOORD_EULERIAN_Z, &
                 "default coord_type should be VCOORD_EULERIAN_Z")
      call vc%destroy()
   end subroutine test_default_coord

   subroutine test_destroy(error)
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      call grid%init(2, 2, 1, 1.0_wp, 1.0_wp)
      call vc%init(grid, nz_ml=4)
      call vc%destroy()
      call check(error,.not. vc%is_init, "is_init should be false after destroy")
      if (allocated(error)) return
      call check(error,.not. allocated(vc%dsig), "dsig should be deallocated")
      if (allocated(error)) return
      call check(error,.not. allocated(vc%target_h), "target_h should be deallocated")
      if (allocated(error)) return
      call check(error,.not. allocated(vc%z_ref), "z_ref should be deallocated")
   end subroutine test_destroy

   subroutine test_parse(error)
      type(error_type), allocatable, intent(out) :: error
      call check(error, parse_ocean_vcoord_type("eulerian_z") == VCOORD_EULERIAN_Z, &
                 "parse 'eulerian_z' -> VCOORD_EULERIAN_Z")
      if (allocated(error)) return
      call check(error, parse_ocean_vcoord_type("sigma") == VCOORD_SIGMA, &
                 "parse 'sigma' -> VCOORD_SIGMA")
      if (allocated(error)) return
      call check(error, parse_ocean_vcoord_type("ZSIGMA") == VCOORD_ZSIGMA, &
                 "parse 'ZSIGMA' (upper) -> VCOORD_ZSIGMA")
      if (allocated(error)) return
      call check(error, parse_ocean_vcoord_type("zstar") == VCOORD_ZSTAR, &
                 "parse 'zstar' -> VCOORD_ZSTAR")
      if (allocated(error)) return
      call check(error, parse_ocean_vcoord_type("zstar_full") == VCOORD_ZSTAR_FULL, &
                 "parse 'zstar_full' -> VCOORD_ZSTAR_FULL")
      if (allocated(error)) return
      call check(error, parse_ocean_vcoord_type("zstar_sigma") == VCOORD_ZSTAR_SIGMA, &
                 "parse 'zstar_sigma' -> VCOORD_ZSTAR_SIGMA")
      if (allocated(error)) return
      call check(error, parse_ocean_vcoord_type("bogus") == VCOORD_EULERIAN_Z, &
                 "unrecognised name should fall back to VCOORD_EULERIAN_Z")
   end subroutine test_parse

   subroutine test_target_h_eulerian(error)
      !! VCOORD_EULERIAN_Z: target_h(:,:,k) = H · dsig(k), independent of η.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      integer, parameter :: NZ = 4
      integer, parameter :: NX = 4, NY = 4
      real(wp) :: total_h(NX + 2, NY + 2), eta(NX + 2, NY + 2)
      real(wp) :: max_err
      integer :: i, j, k, nx_tot, ny_tot

      call grid%init(NX, NY, 1, 1.0_wp, 1.0_wp)
      call vc%init(grid, nz_ml=NZ)
      vc%coord_type = VCOORD_EULERIAN_Z
      nx_tot = grid%nx_total
      ny_tot = grid%ny_total

      total_h = 1000.0_wp
      eta = 5.0_wp        ! deliberately non-zero — should NOT appear
      call vc%compute_target_h(total_h(1:nx_tot, 1:ny_tot), eta(1:nx_tot, 1:ny_tot))

      max_err = 0.0_wp
      do k = 1, NZ
         do j = 1, ny_tot
            do i = 1, nx_tot
               max_err = max(max_err, abs(vc%target_h(i, j, k) - 1000.0_wp/real(NZ, wp)))
            end do
         end do
      end do

      call check(error, max_err < 1.0e-10_wp, &
                 "VCOORD_EULERIAN_Z target_h should equal H/nz_ml regardless of eta")
      call vc%destroy()
   end subroutine test_target_h_eulerian

   subroutine test_target_h_sigma(error)
      !! VCOORD_SIGMA: sum_k target_h(i,j,k) = H(i,j) + η(i,j) per column.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      integer, parameter :: NZ = 5
      integer, parameter :: NX = 4, NY = 3
      real(wp) :: total_h(NX + 2, NY + 2), eta(NX + 2, NY + 2)
      real(wp) :: col_sum, expected, max_drift
      integer :: i, j, k, nx_tot, ny_tot

      call grid%init(NX, NY, 1, 1.0_wp, 1.0_wp)
      call vc%init(grid, nz_ml=NZ)
      vc%coord_type = VCOORD_SIGMA
      nx_tot = grid%nx_total
      ny_tot = grid%ny_total

      ! Non-trivial H and η — spatially varying both, so conservation
      ! has to hold per-column, not just in aggregate.
      do j = 1, ny_tot
         do i = 1, nx_tot
            total_h(i, j) = 100.0_wp + 5.0_wp*real(i, wp) - 2.0_wp*real(j, wp)
            eta(i, j) = 0.1_wp*real(i - j, wp)
         end do
      end do
      call vc%compute_target_h(total_h(1:nx_tot, 1:ny_tot), eta(1:nx_tot, 1:ny_tot))

      max_drift = 0.0_wp
      do j = 1, ny_tot
         do i = 1, nx_tot
            col_sum = 0.0_wp
            do k = 1, NZ
               col_sum = col_sum + vc%target_h(i, j, k)
            end do
            expected = total_h(i, j) + eta(i, j)
            max_drift = max(max_drift, abs(col_sum - expected))
         end do
      end do

      call check(error, max_drift < 1.0e-10_wp, &
                 "VCOORD_SIGMA target_h must sum to H + eta per column")
      call vc%destroy()
   end subroutine test_target_h_sigma

   subroutine test_target_h_zstar_eta(error)
      !! VCOORD_ZSTAR: the η-perturbation must show up — a pure-η
      !! difference at fixed H gives a target_h difference of η·dsig(k).
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      integer, parameter :: NZ = 4
      integer, parameter :: NX = 2, NY = 2
      real(wp) :: total_h(NX + 2, NY + 2), eta_a(NX + 2, NY + 2), eta_b(NX + 2, NY + 2)
      real(wp), allocatable :: target_a(:, :, :)
      real(wp) :: max_err, expected
      integer :: i, j, k, nx_tot, ny_tot

      call grid%init(NX, NY, 1, 1.0_wp, 1.0_wp)
      call vc%init(grid, nz_ml=NZ)
      vc%coord_type = VCOORD_ZSTAR
      nx_tot = grid%nx_total
      ny_tot = grid%ny_total

      total_h = 800.0_wp
      eta_a = 0.0_wp
      eta_b = 3.0_wp

      call vc%compute_target_h(total_h(1:nx_tot, 1:ny_tot), eta_a(1:nx_tot, 1:ny_tot))
      allocate (target_a(nx_tot, ny_tot, NZ), source=vc%target_h)
      call vc%compute_target_h(total_h(1:nx_tot, 1:ny_tot), eta_b(1:nx_tot, 1:ny_tot))

      max_err = 0.0_wp
      do k = 1, NZ
         expected = 3.0_wp*vc%dsig(k)
         do j = 1, ny_tot
            do i = 1, nx_tot
               max_err = max(max_err, abs((vc%target_h(i, j, k) - target_a(i, j, k)) - expected))
            end do
         end do
      end do

      call check(error, max_err < 1.0e-10_wp, &
                 "VCOORD_ZSTAR target_h shift per layer should equal Δη · dsig(k)")
      deallocate (target_a)
      call vc%destroy()
   end subroutine test_target_h_zstar_eta

   subroutine test_target_h_ratios(error)
      !! Non-uniform `dsig` (e.g., stretched z) must propagate through
      !! `compute_target_h`: target_h(k) / target_h(j) == dsig(k) / dsig(j).
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      integer, parameter :: NZ = 3
      real(wp) :: total_h(3, 3), eta(3, 3)
      real(wp) :: ratio_expected, ratio_actual, err_max
      integer :: nx_tot, ny_tot, i, j

      call grid%init(1, 1, 1, 1.0_wp, 1.0_wp)
      call vc%init(grid, nz_ml=NZ)
      vc%coord_type = VCOORD_SIGMA
      nx_tot = grid%nx_total
      ny_tot = grid%ny_total

      ! Stretched profile: surface (k=NZ) thinner than bed (k=1).
      vc%dsig(1) = 0.6_wp
      vc%dsig(2) = 0.3_wp
      vc%dsig(3) = 0.1_wp

      total_h = 500.0_wp
      eta = 2.0_wp
      call vc%compute_target_h(total_h(1:nx_tot, 1:ny_tot), eta(1:nx_tot, 1:ny_tot))

      ! Check k=3 / k=1 ratio at every interior cell.
      ratio_expected = vc%dsig(3)/vc%dsig(1)
      err_max = 0.0_wp
      do j = 1, ny_tot
         do i = 1, nx_tot
            ratio_actual = vc%target_h(i, j, 3)/vc%target_h(i, j, 1)
            err_max = max(err_max, abs(ratio_actual - ratio_expected))
         end do
      end do

      call check(error, err_max < 1.0e-10_wp, &
                 "VCOORD_SIGMA layer ratios should equal dsig ratios")
      call vc%destroy()
   end subroutine test_target_h_ratios

   ! --------- Layer-1 hybrid coords (ZSIGMA, ZSTAR_SIGMA) ---------

   subroutine setup_zref_for_hybrid_tests(vc, nz)
      !! Stuff a realistic z-level reference profile into `z_ref_global`
      !! for the hybrid-coord tests.  4 reference layers spaced
      !! 0/100/250/500/1000 m, surface-thin / bed-thick.
      type(ocean_vcoord_t), intent(inout) :: vc
      integer, intent(in) :: nz
      if (nz /= 4) return
      vc%z_ref_global(0) = 0.0_wp
      vc%z_ref_global(1) = 100.0_wp
      vc%z_ref_global(2) = 250.0_wp
      vc%z_ref_global(3) = 500.0_wp
      vc%z_ref_global(4) = 1000.0_wp
   end subroutine setup_zref_for_hybrid_tests

   subroutine test_zsigma_shallow(error)
      !! H + η ≤ depth_transition → ZSIGMA collapses to pure SIGMA
      !! identically (no z-level contribution).
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      integer, parameter :: NZ = 4
      integer, parameter :: NX = 3, NY = 3
      real(wp) :: total_h(NX + 2, NY + 2), eta(NX + 2, NY + 2)
      real(wp) :: max_err, expected
      integer :: i, j, k, nx_tot, ny_tot

      call grid%init(NX, NY, 1, 1.0_wp, 1.0_wp)
      call vc%init(grid, nz_ml=NZ)
      vc%coord_type = VCOORD_ZSIGMA
      vc%zsigma_depth_transition = 200.0_wp
      vc%zsigma_blend_width = 100.0_wp
      call setup_zref_for_hybrid_tests(vc, NZ)
      nx_tot = grid%nx_total
      ny_tot = grid%ny_total

      total_h = 150.0_wp        ! H + η = 150 < 200 = transition
      eta = 0.0_wp
      call vc%compute_target_h(total_h(1:nx_tot, 1:ny_tot), eta(1:nx_tot, 1:ny_tot))

      max_err = 0.0_wp
      do k = 1, NZ
         expected = vc%dsig(k)*150.0_wp
         do j = 1, ny_tot
            do i = 1, nx_tot
               max_err = max(max_err, abs(vc%target_h(i, j, k) - expected))
            end do
         end do
      end do
      call check(error, max_err < 1.0e-10_wp, &
                 "ZSIGMA below depth_transition should equal pure SIGMA")
      call vc%destroy()
   end subroutine test_zsigma_shallow

   subroutine test_zsigma_deep(error)
      !! H + η ≥ depth_transition + blend_width → ZSIGMA uses the
      !! z-level branch (α = 1).  The z-level intervals come from
      !! `z_ref_global`; clip + bed-side deficit keep sum = H + η.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      integer, parameter :: NZ = 4
      integer, parameter :: NX = 2, NY = 2
      real(wp) :: total_h(NX + 2, NY + 2), eta(NX + 2, NY + 2)
      real(wp) :: col_sum, surface_h, expected
      integer :: i, j, k, nx_tot, ny_tot
      checks: block

         call grid%init(NX, NY, 1, 1.0_wp, 1.0_wp)
         call vc%init(grid, nz_ml=NZ)
         vc%coord_type = VCOORD_ZSIGMA
         vc%zsigma_depth_transition = 200.0_wp
         vc%zsigma_blend_width = 100.0_wp
         call setup_zref_for_hybrid_tests(vc, NZ)
         nx_tot = grid%nx_total
         ny_tot = grid%ny_total

         ! Column total H+η = 1200 m > 200 + 100 = 300 (well past blend zone).
         ! z_ref_global is 0/100/250/500/1000, so the deepest interface is
         ! shallower than H+η — no clipping needed, layers match z_ref intervals.
         total_h = 1200.0_wp
         eta = 0.0_wp
         call vc%compute_target_h(total_h(1:nx_tot, 1:ny_tot), eta(1:nx_tot, 1:ny_tot))

         ! Surface layer (k=NZ) should match z_ref interval z_ref(0)..z_ref(1) = 100 m.
         ! Bed-side layer (k=1) carries the deficit (1200 - 1000 = 200 from clipping
         ! against H, plus the z_ref(3)..z_ref(4) = 500 m baseline = 700 m).
         expected = 100.0_wp          ! z_ref(0)..z_ref(1) = surface layer
         surface_h = 0.0_wp
         do j = 1, ny_tot
            do i = 1, nx_tot
               surface_h = max(surface_h, abs(vc%target_h(i, j, NZ) - expected))
            end do
         end do
         call check(error, surface_h < 1.0e-10_wp, &
                    "ZSIGMA deep: surface layer should match z_ref(0..1) interval")
         if (allocated(error)) exit checks

         ! Per-column conservation.
         do j = 1, ny_tot
            do i = 1, nx_tot
               col_sum = 0.0_wp
               do k = 1, NZ
                  col_sum = col_sum + vc%target_h(i, j, k)
               end do
               call check(error, abs(col_sum - 1200.0_wp) < 1.0e-9_wp, &
                          "ZSIGMA deep: sum(target_h) should equal column total")
               if (allocated(error)) exit checks
            end do
         end do
      end block checks
      call vc%destroy()
   end subroutine test_zsigma_deep

   subroutine test_zsigma_blend_monotone(error)
      !! In the blend zone (transition < H ≤ transition+blend_width)
      !! the smoothstep α grows monotonically with H, so the deep-z
      !! contribution must grow correspondingly.  Probe three depths
      !! and verify monotonic α via target_h dependence.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      integer, parameter :: NZ = 4
      real(wp) :: total_h(3, 3), eta(3, 3)
      real(wp), allocatable :: h_low(:, :, :), h_mid(:, :, :), h_high(:, :, :)
      real(wp) :: diff_lo, diff_hi
      integer :: nx_tot, ny_tot

      call grid%init(1, 1, 1, 1.0_wp, 1.0_wp)
      call vc%init(grid, nz_ml=NZ)
      vc%coord_type = VCOORD_ZSIGMA
      vc%zsigma_depth_transition = 200.0_wp
      vc%zsigma_blend_width = 200.0_wp        ! transition zone 200..400
      call setup_zref_for_hybrid_tests(vc, NZ)
      nx_tot = grid%nx_total
      ny_tot = grid%ny_total

      eta = 0.0_wp
      total_h = 250.0_wp          ! α ≈ 0.16
      call vc%compute_target_h(total_h(1:nx_tot, 1:ny_tot), eta(1:nx_tot, 1:ny_tot))
      allocate (h_low(nx_tot, ny_tot, NZ), source=vc%target_h)

      total_h = 300.0_wp          ! α = 0.5
      call vc%compute_target_h(total_h(1:nx_tot, 1:ny_tot), eta(1:nx_tot, 1:ny_tot))
      allocate (h_mid(nx_tot, ny_tot, NZ), source=vc%target_h)

      total_h = 380.0_wp          ! α ≈ 0.95
      call vc%compute_target_h(total_h(1:nx_tot, 1:ny_tot), eta(1:nx_tot, 1:ny_tot))
      allocate (h_high(nx_tot, ny_tot, NZ), source=vc%target_h)

      ! Surface layer thickness: ZSIGMA pulls toward 100 m as α grows.
      ! At α=0, sigma gives 100/4 ratio of H = 25 m at H=100. At α=1,
      ! it's exactly 100 m.  So surface-h grows monotonically with α
      ! when H is held variable but z_ref dominates.  Compare deltas
      ! after subtracting the sigma component.
      diff_lo = h_mid(1, 1, NZ) - h_low(1, 1, NZ)
      diff_hi = h_high(1, 1, NZ) - h_mid(1, 1, NZ)
      call check(error, diff_lo > 0.0_wp .and. diff_hi > 0.0_wp, &
                 "ZSIGMA surface-h should be monotonically increasing across the blend zone")
      deallocate (h_low, h_mid, h_high)
      call vc%destroy()
   end subroutine test_zsigma_blend_monotone

   subroutine test_zsigma_conservation(error)
      !! ZSIGMA must conserve sum_k target_h = H + η across every regime
      !! (pure sigma, blend zone, pure z-level + clipping).  Probe three
      !! H values spanning all three.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      integer, parameter :: NZ = 4
      real(wp) :: total_h(3, 3), eta(3, 3)
      real(wp), parameter :: HVALS(3) = [150.0_wp, 300.0_wp, 1200.0_wp]
      real(wp) :: col_sum, expected, max_drift
      integer :: ihval, nx_tot, ny_tot, i, j, k

      call grid%init(1, 1, 1, 1.0_wp, 1.0_wp)
      call vc%init(grid, nz_ml=NZ)
      vc%coord_type = VCOORD_ZSIGMA
      vc%zsigma_depth_transition = 200.0_wp
      vc%zsigma_blend_width = 200.0_wp
      call setup_zref_for_hybrid_tests(vc, NZ)
      nx_tot = grid%nx_total
      ny_tot = grid%ny_total

      eta = 0.5_wp
      max_drift = 0.0_wp
      do ihval = 1, 3
         total_h = HVALS(ihval)
         call vc%compute_target_h(total_h(1:nx_tot, 1:ny_tot), eta(1:nx_tot, 1:ny_tot))
         expected = HVALS(ihval) + 0.5_wp
         do j = 1, ny_tot
            do i = 1, nx_tot
               col_sum = 0.0_wp
               do k = 1, NZ
                  col_sum = col_sum + vc%target_h(i, j, k)
               end do
               max_drift = max(max_drift, abs(col_sum - expected))
            end do
         end do
      end do
      call check(error, max_drift < 1.0e-9_wp, &
                 "ZSIGMA: sum(target_h) should equal H + eta in every regime")
      call vc%destroy()
   end subroutine test_zsigma_conservation

   subroutine test_zstar_sigma_shallow(error)
      !! Shallow ZSTAR_SIGMA → pure SIGMA.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      integer, parameter :: NZ = 4
      integer, parameter :: NX = 3, NY = 3
      real(wp) :: total_h(NX + 2, NY + 2), eta(NX + 2, NY + 2)
      real(wp) :: max_err, expected
      integer :: i, j, k, nx_tot, ny_tot

      call grid%init(NX, NY, 1, 1.0_wp, 1.0_wp)
      call vc%init(grid, nz_ml=NZ)
      vc%coord_type = VCOORD_ZSTAR_SIGMA
      vc%zsigma_depth_transition = 200.0_wp
      vc%zsigma_blend_width = 100.0_wp
      call setup_zref_for_hybrid_tests(vc, NZ)
      nx_tot = grid%nx_total
      ny_tot = grid%ny_total

      total_h = 100.0_wp
      eta = 0.0_wp
      call vc%compute_target_h(total_h(1:nx_tot, 1:ny_tot), eta(1:nx_tot, 1:ny_tot))

      max_err = 0.0_wp
      do k = 1, NZ
         expected = vc%dsig(k)*100.0_wp
         do j = 1, ny_tot
            do i = 1, nx_tot
               max_err = max(max_err, abs(vc%target_h(i, j, k) - expected))
            end do
         end do
      end do
      call check(error, max_err < 1.0e-10_wp, &
                 "ZSTAR_SIGMA below depth_transition should equal pure SIGMA")
      call vc%destroy()
   end subroutine test_zstar_sigma_shallow

   subroutine test_zstar_sigma_deep(error)
      !! Deep ZSTAR_SIGMA → pure z*-lite, target_h(k) = (z_ref(nz-k+1) -
      !! z_ref(nz-k)) · column_total / z_ref(nz).  No clipping needed
      !! because the formula stretches per column.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      integer, parameter :: NZ = 4
      integer, parameter :: NX = 2, NY = 2
      real(wp) :: total_h(NX + 2, NY + 2), eta(NX + 2, NY + 2)
      real(wp) :: col_total, expected_k, max_err
      integer :: i, j, k, nx_tot, ny_tot

      call grid%init(NX, NY, 1, 1.0_wp, 1.0_wp)
      call vc%init(grid, nz_ml=NZ)
      vc%coord_type = VCOORD_ZSTAR_SIGMA
      vc%zsigma_depth_transition = 200.0_wp
      vc%zsigma_blend_width = 100.0_wp
      call setup_zref_for_hybrid_tests(vc, NZ)
      nx_tot = grid%nx_total
      ny_tot = grid%ny_total

      total_h = 2000.0_wp        ! H + η = 2000 >> 300 = transition + width
      eta = 0.0_wp
      call vc%compute_target_h(total_h(1:nx_tot, 1:ny_tot), eta(1:nx_tot, 1:ny_tot))

      col_total = 2000.0_wp
      max_err = 0.0_wp
      do k = 1, NZ
         expected_k = (vc%z_ref_global(NZ - k + 1) - vc%z_ref_global(NZ - k)) &
                      *col_total/vc%z_ref_global(NZ)
         do j = 1, ny_tot
            do i = 1, nx_tot
               max_err = max(max_err, abs(vc%target_h(i, j, k) - expected_k))
            end do
         end do
      end do
      call check(error, max_err < 1.0e-9_wp, &
                 "ZSTAR_SIGMA deep: target_h should match z*-lite stretched formula")
      call vc%destroy()
   end subroutine test_zstar_sigma_deep

   subroutine test_zstar_sigma_conservation(error)
      !! Both branches sum to H + η by construction so the blend must too.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      integer, parameter :: NZ = 4
      real(wp) :: total_h(3, 3), eta(3, 3)
      real(wp), parameter :: HVALS(3) = [120.0_wp, 280.0_wp, 1500.0_wp]
      real(wp) :: col_sum, expected, max_drift
      integer :: ihval, nx_tot, ny_tot, i, j, k

      call grid%init(1, 1, 1, 1.0_wp, 1.0_wp)
      call vc%init(grid, nz_ml=NZ)
      vc%coord_type = VCOORD_ZSTAR_SIGMA
      vc%zsigma_depth_transition = 200.0_wp
      vc%zsigma_blend_width = 200.0_wp
      call setup_zref_for_hybrid_tests(vc, NZ)
      nx_tot = grid%nx_total
      ny_tot = grid%ny_total

      eta = -0.3_wp
      max_drift = 0.0_wp
      do ihval = 1, 3
         total_h = HVALS(ihval)
         call vc%compute_target_h(total_h(1:nx_tot, 1:ny_tot), eta(1:nx_tot, 1:ny_tot))
         expected = HVALS(ihval) - 0.3_wp
         do j = 1, ny_tot
            do i = 1, nx_tot
               col_sum = 0.0_wp
               do k = 1, NZ
                  col_sum = col_sum + vc%target_h(i, j, k)
               end do
               max_drift = max(max_drift, abs(col_sum - expected))
            end do
         end do
      end do
      call check(error, max_drift < 1.0e-9_wp, &
                 "ZSTAR_SIGMA: sum(target_h) should equal H + eta in every regime")
      call vc%destroy()
   end subroutine test_zstar_sigma_conservation

   ! --------- Layer-2 VCOORD_ZSTAR_FULL ---------

   subroutine test_zstar_full_build_uniform(error)
      !! `build_zref_full` with default (auto) parameters and uniform
      !! stretching on a flat bath: z_ref is monotonically increasing
      !! per column, z_ref(:, :, nz) == h_bed, z_ref(:, :, 0) == 0.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      integer, parameter :: NZ = 6
      integer, parameter :: NX = 3, NY = 3
      real(wp) :: h_bed(NX + 2, NY + 2)
      logical :: monotone_ok, bed_ok, surface_ok
      integer :: i, j, k, nx_tot, ny_tot
      checks: block

         call grid%init(NX, NY, 1, 1.0_wp, 1.0_wp)
         call vc%init(grid, nz_ml=NZ)
         vc%coord_type = VCOORD_ZSTAR_FULL
         vc%zstar_h_surf_target = 10.0_wp
         vc%zstar_n_surf = 2
         vc%zstar_stretching = STRETCH_UNIFORM
         nx_tot = grid%nx_total
         ny_tot = grid%ny_total

         h_bed = 1000.0_wp
         call vc%build_zref_full(h_bed(1:nx_tot, 1:ny_tot))

         monotone_ok = .true.
         bed_ok = .true.
         surface_ok = .true.
         do j = 1, ny_tot
            do i = 1, nx_tot
               if (vc%z_ref(i, j, 0) /= 0.0_wp) surface_ok = .false.
               if (abs(vc%z_ref(i, j, NZ) - 1000.0_wp) > 1.0e-10_wp) bed_ok = .false.
               do k = 1, NZ
                  if (vc%z_ref(i, j, k) <= vc%z_ref(i, j, k - 1)) monotone_ok = .false.
               end do
            end do
         end do
         call check(error, surface_ok, "ZSTAR_FULL: z_ref(:,:,0) should be 0")
         if (allocated(error)) exit checks
         call check(error, bed_ok, "ZSTAR_FULL: z_ref(:,:,nz) should equal h_bed")
         if (allocated(error)) exit checks
         call check(error, monotone_ok, "ZSTAR_FULL: z_ref should be strictly monotonic per column")
      end block checks
      call vc%destroy()
   end subroutine test_zstar_full_build_uniform

   subroutine test_zstar_full_build_surface(error)
      !! Uniform stretching with `n_surf=2, h_surf_target=10`: the top
      !! two layers should be exactly 10 m thick each.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      integer, parameter :: NZ = 6
      real(wp) :: h_bed(3, 3)
      real(wp) :: surf1, surf2, max_err1, max_err2
      integer :: i, j, nx_tot, ny_tot

      call grid%init(1, 1, 1, 1.0_wp, 1.0_wp)
      call vc%init(grid, nz_ml=NZ)
      vc%coord_type = VCOORD_ZSTAR_FULL
      vc%zstar_h_surf_target = 10.0_wp
      vc%zstar_n_surf = 2
      vc%zstar_stretching = STRETCH_UNIFORM
      nx_tot = grid%nx_total
      ny_tot = grid%ny_total

      h_bed = 500.0_wp
      call vc%build_zref_full(h_bed(1:nx_tot, 1:ny_tot))

      ! z_ref top-down: z_ref(0)=0, z_ref(1)=10, z_ref(2)=20, then uniform.
      max_err1 = 0.0_wp
      max_err2 = 0.0_wp
      do j = 1, ny_tot
         do i = 1, nx_tot
            surf1 = vc%z_ref(i, j, 1) - vc%z_ref(i, j, 0)
            surf2 = vc%z_ref(i, j, 2) - vc%z_ref(i, j, 1)
            max_err1 = max(max_err1, abs(surf1 - 10.0_wp))
            max_err2 = max(max_err2, abs(surf2 - 10.0_wp))
         end do
      end do
      call check(error, max_err1 < 1.0e-10_wp .and. max_err2 < 1.0e-10_wp, &
                 "ZSTAR_FULL UNIFORM: top n_surf layers should equal h_surf_target")
      call vc%destroy()
   end subroutine test_zstar_full_build_surface

   subroutine test_zstar_full_build_log(error)
      !! LOG stretching gives layer thicknesses that grow with depth.
      !! Verify: top layer ~ h_surf_target; layer k thickness < layer
      !! k+1 thickness across the fine zone.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      integer, parameter :: NZ = 6
      real(wp) :: h_bed(3, 3)
      real(wp) :: dz_top, dz_next, total
      logical :: top_ok, growth_ok, total_ok
      integer :: i, j, k, nx_tot, ny_tot
      checks: block

         call grid%init(1, 1, 1, 1.0_wp, 1.0_wp)
         call vc%init(grid, nz_ml=NZ)
         vc%coord_type = VCOORD_ZSTAR_FULL
         vc%zstar_h_surf_target = 5.0_wp
         vc%zstar_n_surf = 3
         vc%zstar_stretching = STRETCH_LOG
         nx_tot = grid%nx_total
         ny_tot = grid%ny_total

         h_bed = 1000.0_wp
         call vc%build_zref_full(h_bed(1:nx_tot, 1:ny_tot))

         top_ok = .true.
         growth_ok = .true.
         total_ok = .true.
         do j = 1, ny_tot
            do i = 1, nx_tot
               dz_top = vc%z_ref(i, j, 1) - vc%z_ref(i, j, 0)
               dz_next = vc%z_ref(i, j, 2) - vc%z_ref(i, j, 1)
               if (abs(dz_top - 5.0_wp) > 1.0e-9_wp) top_ok = .false.
               if (dz_next <= dz_top) growth_ok = .false.
               total = 0.0_wp
               do k = 1, NZ
                  total = total + (vc%z_ref(i, j, k) - vc%z_ref(i, j, k - 1))
               end do
               if (abs(total - 1000.0_wp) > 1.0e-7_wp) total_ok = .false.
            end do
         end do
         call check(error, top_ok, "ZSTAR_FULL LOG: top layer should match h_surf_target")
         if (allocated(error)) exit checks
         call check(error, growth_ok, "ZSTAR_FULL LOG: layer thickness should grow with depth")
         if (allocated(error)) exit checks
         call check(error, total_ok, "ZSTAR_FULL LOG: layer thicknesses should sum to h_bed")
      end block checks
      call vc%destroy()
   end subroutine test_zstar_full_build_log

   subroutine test_zstar_full_at_ref(error)
      !! When H = h_bed (η=0), target_h(k) equals the reference layer
      !! thickness z_ref(nz-k+1) - z_ref(nz-k) for every k.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      integer, parameter :: NZ = 4
      integer, parameter :: NX = 2, NY = 2
      real(wp) :: h_bed(NX + 2, NY + 2)
      real(wp) :: total_h(NX + 2, NY + 2), eta(NX + 2, NY + 2)
      real(wp) :: max_err, expected_k
      integer :: i, j, k, nx_tot, ny_tot

      call grid%init(NX, NY, 1, 1.0_wp, 1.0_wp)
      call vc%init(grid, nz_ml=NZ)
      vc%coord_type = VCOORD_ZSTAR_FULL
      vc%zstar_h_surf_target = 20.0_wp
      vc%zstar_n_surf = 2
      vc%zstar_stretching = STRETCH_UNIFORM
      nx_tot = grid%nx_total
      ny_tot = grid%ny_total

      h_bed = 400.0_wp
      call vc%build_zref_full(h_bed(1:nx_tot, 1:ny_tot))

      ! H = h_bed exactly (η = 0).  Pass total_h = h_bed, eta = 0.
      total_h = h_bed
      eta = 0.0_wp
      call vc%compute_target_h(total_h(1:nx_tot, 1:ny_tot), eta(1:nx_tot, 1:ny_tot))

      max_err = 0.0_wp
      do k = 1, NZ
         do j = 1, ny_tot
            do i = 1, nx_tot
               expected_k = vc%z_ref(i, j, NZ - k + 1) - vc%z_ref(i, j, NZ - k)
               max_err = max(max_err, abs(vc%target_h(i, j, k) - expected_k))
            end do
         end do
      end do
      call check(error, max_err < 1.0e-10_wp, &
                 "ZSTAR_FULL at reference: target_h should match z_ref intervals ROMS-ordered")
      call vc%destroy()
   end subroutine test_zstar_full_at_ref

   subroutine test_zstar_full_pos_eta(error)
      !! η > 0: only the surface layer (k=NZ) grows.  Subsurface layers
      !! stay at reference thickness.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      integer, parameter :: NZ = 4
      real(wp) :: h_bed(3, 3), total_h(3, 3), eta(3, 3)
      real(wp), allocatable :: target_at_ref(:, :, :)
      real(wp) :: max_subsurface_drift, surface_growth
      integer :: i, j, k, nx_tot, ny_tot
      checks: block

         call grid%init(1, 1, 1, 1.0_wp, 1.0_wp)
         call vc%init(grid, nz_ml=NZ)
         vc%coord_type = VCOORD_ZSTAR_FULL
         vc%zstar_h_surf_target = 25.0_wp
         vc%zstar_n_surf = 2
         vc%zstar_stretching = STRETCH_UNIFORM
         nx_tot = grid%nx_total
         ny_tot = grid%ny_total

         h_bed = 500.0_wp
         call vc%build_zref_full(h_bed(1:nx_tot, 1:ny_tot))

         ! Snapshot at η=0.
         total_h = h_bed
         eta = 0.0_wp
         call vc%compute_target_h(total_h(1:nx_tot, 1:ny_tot), eta(1:nx_tot, 1:ny_tot))
         allocate (target_at_ref(nx_tot, ny_tot, NZ), source=vc%target_h)

         ! Now η = +7.5 m.
         eta = 7.5_wp
         call vc%compute_target_h(total_h(1:nx_tot, 1:ny_tot), eta(1:nx_tot, 1:ny_tot))

         max_subsurface_drift = 0.0_wp
         surface_growth = 0.0_wp
         do j = 1, ny_tot
            do i = 1, nx_tot
               do k = 1, NZ - 1
                  max_subsurface_drift = max(max_subsurface_drift, &
                                             abs(vc%target_h(i, j, k) - target_at_ref(i, j, k)))
               end do
               surface_growth = max(surface_growth, &
                                    abs((vc%target_h(i, j, NZ) - target_at_ref(i, j, NZ)) - 7.5_wp))
            end do
         end do
         call check(error, max_subsurface_drift < 1.0e-10_wp, &
                    "ZSTAR_FULL η>0: subsurface layers should keep reference thickness")
         if (allocated(error)) exit checks
         call check(error, surface_growth < 1.0e-10_wp, &
                    "ZSTAR_FULL η>0: surface layer should absorb the full η offset")
      end block checks
      deallocate (target_at_ref)
      call vc%destroy()
   end subroutine test_zstar_full_pos_eta

   subroutine test_zstar_full_neg_eta(error)
      !! η < 0 (column shallower than reference): bed-side layers vanish
      !! to h_min when their interface drops below H, surface trims so
      !! sum = H exactly.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      integer, parameter :: NZ = 4
      real(wp) :: h_bed(3, 3), total_h(3, 3), eta(3, 3)
      logical :: bed_vanishing
      real(wp) :: col_sum, max_drift
      integer :: i, j, k, nx_tot, ny_tot
      real(wp) :: H_eff
      checks: block

         call grid%init(1, 1, 1, 1.0_wp, 1.0_wp)
         call vc%init(grid, nz_ml=NZ)
         vc%coord_type = VCOORD_ZSTAR_FULL
         vc%zstar_h_surf_target = 20.0_wp
         vc%zstar_n_surf = 2
         vc%zstar_stretching = STRETCH_UNIFORM
         vc%zstar_h_min = 1.0e-4_wp
         nx_tot = grid%nx_total
         ny_tot = grid%ny_total

         h_bed = 400.0_wp
         call vc%build_zref_full(h_bed(1:nx_tot, 1:ny_tot))

         ! Drive the column far below reference: H_eff = 60 m << z_ref(NZ).
         ! z_ref pattern: 0/20/40/220/400 — at H_eff=60 the bottom layer
         ! (k=1, interval 220..400) should vanish to h_min.
         total_h = h_bed
         eta = -340.0_wp        ! H_eff = 400 - 340 = 60 m
         call vc%compute_target_h(total_h(1:nx_tot, 1:ny_tot), eta(1:nx_tot, 1:ny_tot))

         bed_vanishing = .true.
         max_drift = 0.0_wp
         H_eff = 60.0_wp
         do j = 1, ny_tot
            do i = 1, nx_tot
               ! k=1 corresponds to z_upper=z_ref(NZ-1)=220, z_lower=z_ref(NZ)=400 — fully below H_eff.
               if (abs(vc%target_h(i, j, 1) - vc%zstar_h_min) > 1.0e-10_wp) bed_vanishing = .false.
               col_sum = 0.0_wp
               do k = 1, NZ
                  col_sum = col_sum + vc%target_h(i, j, k)
               end do
               max_drift = max(max_drift, abs(col_sum - H_eff))
            end do
         end do
         call check(error, bed_vanishing, &
                    "ZSTAR_FULL η<0: bed-side layer below H should vanish to h_min")
         if (allocated(error)) exit checks
         call check(error, max_drift < 1.0e-9_wp, &
                    "ZSTAR_FULL η<0: surface trim should restore sum(target_h) = H")
      end block checks
      call vc%destroy()
   end subroutine test_zstar_full_neg_eta

   subroutine test_zstar_full_conservation(error)
      !! Sum conservation across a sweep of η values spanning the
      !! three regimes (η>0, η=0, η<0 with partial clipping, η<<0 with
      !! full vanishing).
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      integer, parameter :: NZ = 4
      real(wp) :: h_bed(3, 3), total_h(3, 3), eta(3, 3)
      real(wp), parameter :: ETAS(4) = [5.0_wp, 0.0_wp, -50.0_wp, -300.0_wp]
      real(wp) :: col_sum, expected, max_drift
      integer :: ie, nx_tot, ny_tot, i, j, k

      call grid%init(1, 1, 1, 1.0_wp, 1.0_wp)
      call vc%init(grid, nz_ml=NZ)
      vc%coord_type = VCOORD_ZSTAR_FULL
      vc%zstar_h_surf_target = 15.0_wp
      vc%zstar_n_surf = 2
      vc%zstar_stretching = STRETCH_UNIFORM
      vc%zstar_h_min = 1.0e-4_wp
      nx_tot = grid%nx_total
      ny_tot = grid%ny_total

      h_bed = 400.0_wp
      call vc%build_zref_full(h_bed(1:nx_tot, 1:ny_tot))
      total_h = h_bed

      max_drift = 0.0_wp
      do ie = 1, size(ETAS)
         eta = ETAS(ie)
         call vc%compute_target_h(total_h(1:nx_tot, 1:ny_tot), eta(1:nx_tot, 1:ny_tot))
         expected = max(400.0_wp + ETAS(ie), 0.0_wp)
         do j = 1, ny_tot
            do i = 1, nx_tot
               col_sum = 0.0_wp
               do k = 1, NZ
                  col_sum = col_sum + vc%target_h(i, j, k)
               end do
               max_drift = max(max_drift, abs(col_sum - expected))
            end do
         end do
      end do
      ! Tolerance allows for the "column too thin to honour the floor"
      ! corner case where the surface is pinned at h_min — that's a
      ! deliberate physics-guard, not a bug.  At ETA = -300 (H_eff =
      ! 100 m), the column still has plenty of headroom above h_min so
      ! the trim is exact.
      call check(error, max_drift < 1.0e-9_wp, &
                 "ZSTAR_FULL: sum(target_h) should equal max(H, 0) across all η regimes")
      call vc%destroy()
   end subroutine test_zstar_full_conservation

   subroutine test_parse_lagrangian(error)
      !! `VCOORD_LAGRANGIAN` is reached via two synonymous namelist
      !! strings — `"lagrangian"` (the mechanism: layers track flow)
      !! and `"isopycnal"` (the physics: interfaces follow density
      !! surfaces).  Both casings work via the existing case-folding
      !! pattern in `parse_ocean_vcoord_type`.  This test pins the
      !! mapping so future namespace tidying can't silently break the
      !! user-facing tokens.
      type(error_type), allocatable, intent(out) :: error
      call check(error, parse_ocean_vcoord_type("lagrangian") == VCOORD_LAGRANGIAN, &
                 "parse 'lagrangian' -> VCOORD_LAGRANGIAN")
      if (allocated(error)) return
      call check(error, parse_ocean_vcoord_type("LAGRANGIAN") == VCOORD_LAGRANGIAN, &
                 "parse 'LAGRANGIAN' (upper) -> VCOORD_LAGRANGIAN")
      if (allocated(error)) return
      call check(error, parse_ocean_vcoord_type("isopycnal") == VCOORD_LAGRANGIAN, &
                 "parse 'isopycnal' (alias) -> VCOORD_LAGRANGIAN")
      if (allocated(error)) return
      call check(error, parse_ocean_vcoord_type("ISOPYCNAL") == VCOORD_LAGRANGIAN, &
                 "parse 'ISOPYCNAL' (upper alias) -> VCOORD_LAGRANGIAN")
      if (allocated(error)) return
      ! VCOORD_LAGRANGIAN must be distinguishable from the other codes;
      ! its sentinel value (-1) is intentionally outside the positive
      ! VCOORD_* enum range used by the coastal path.
      call check(error, VCOORD_LAGRANGIAN /= VCOORD_EULERIAN_Z, &
                 "VCOORD_LAGRANGIAN must differ from VCOORD_EULERIAN_Z")
      if (allocated(error)) return
      call check(error, VCOORD_LAGRANGIAN /= VCOORD_SIGMA, &
                 "VCOORD_LAGRANGIAN must differ from VCOORD_SIGMA")
   end subroutine test_parse_lagrangian

   subroutine test_target_h_lagrangian_noop(error)
      !! `VCOORD_LAGRANGIAN` means layers float freely under the
      !! continuity step — there's no geometric `target_h` to
      !! restore them to.  `compute_target_h` must therefore
      !! early-return without writing to `vc%target_h`, so any
      !! caller-side accidental read of `target_h` will see its
      !! alloc-time value rather than a stale geometric layout
      !! that would mislead the ALE remap.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      integer, parameter :: NZ = 4
      integer, parameter :: NX = 4, NY = 4
      real(wp) :: total_h(NX + 2, NY + 2), eta(NX + 2, NY + 2)
      real(wp), allocatable :: target_h_before(:, :, :)
      real(wp) :: max_drift
      integer :: nx_tot, ny_tot

      call grid%init(NX, NY, 1, 1.0_wp, 1.0_wp)
      call vc%init(grid, nz_ml=NZ)
      vc%coord_type = VCOORD_LAGRANGIAN
      nx_tot = grid%nx_total
      ny_tot = grid%ny_total

      ! Snapshot target_h at its post-init state, then call compute_target_h
      ! with non-trivial inputs.  The Lagrangian early-return guarantees
      ! the snapshot equals the post-call value.
      allocate (target_h_before, source=vc%target_h)
      total_h = 1000.0_wp
      eta = 5.0_wp     ! deliberately non-zero — should NOT appear
      call vc%compute_target_h(total_h(1:nx_tot, 1:ny_tot), eta(1:nx_tot, 1:ny_tot))

      max_drift = maxval(abs(vc%target_h - target_h_before))
      call check(error, max_drift < 1.0e-15_wp, &
                 "VCOORD_LAGRANGIAN: compute_target_h should leave target_h untouched")

      deallocate (target_h_before)
      call vc%destroy()
   end subroutine test_target_h_lagrangian_noop

   ! -----------------------------------------------------------------
   ! VCOORD_Z_FIXED cases (MOM6 gprime-style interface anchoring)
   ! -----------------------------------------------------------------
   !
   ! `vcoord_type = "z_fixed"` (with "z_levels", "gprime" as aliases)
   ! sets layer interfaces at fixed depths `z = k · h_ref / nz_ml` from
   ! the surface.  In cells where the local water column is deep enough
   ! every layer takes its nominal thickness; in shallow cells the
   ! bed-side layers vanish to `zstar_h_min` and the surface layer
   ! absorbs the residual.  When `z_fixed_h_ref` is left at its default
   ! 0 the path falls back to uniform `H · dsig(k)` so untouched tests
   ! still get a sensible target.

   subroutine test_parse_z_fixed(error)
      !! All four casing/spelling variants — `"z_fixed"`, `"Z_FIXED"`,
      !! `"z_levels"`, `"Z_LEVELS"`, `"gprime"`, `"GPRIME"` — must map
      !! to `VCOORD_Z_FIXED`.  Pins the user-facing tokens so any
      !! future namespace cleanup can't silently break the namelist
      !! contract.
      type(error_type), allocatable, intent(out) :: error
      call check(error, parse_ocean_vcoord_type("z_fixed") == VCOORD_Z_FIXED, &
                 "parse 'z_fixed' -> VCOORD_Z_FIXED")
      if (allocated(error)) return
      call check(error, parse_ocean_vcoord_type("Z_FIXED") == VCOORD_Z_FIXED, &
                 "parse 'Z_FIXED' (upper) -> VCOORD_Z_FIXED")
      if (allocated(error)) return
      call check(error, parse_ocean_vcoord_type("z_levels") == VCOORD_Z_FIXED, &
                 "parse 'z_levels' (alias) -> VCOORD_Z_FIXED")
      if (allocated(error)) return
      call check(error, parse_ocean_vcoord_type("gprime") == VCOORD_Z_FIXED, &
                 "parse 'gprime' (alias) -> VCOORD_Z_FIXED")
      if (allocated(error)) return
      call check(error, parse_ocean_vcoord_type("GPRIME") == VCOORD_Z_FIXED, &
                 "parse 'GPRIME' (upper alias) -> VCOORD_Z_FIXED")
      if (allocated(error)) return
      call check(error, VCOORD_Z_FIXED /= VCOORD_EULERIAN_Z, &
                 "VCOORD_Z_FIXED must differ from VCOORD_EULERIAN_Z")
   end subroutine test_parse_z_fixed

   subroutine test_z_fixed_deep(error)
      !! Deep column: H = h_ref = 1000 m, nz=2 ⇒ h_nominal = 500 m.
      !! Both layers should land at exactly 500 m everywhere — the
      !! interface locks at z = -500 m.  This is the "gprime canonical"
      !! configuration MOM6 uses for the double_gyre reference.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      integer, parameter :: NZ = 2
      integer, parameter :: NX = 4, NY = 4
      real(wp), parameter :: H_REF = 1000.0_wp
      real(wp) :: total_h(NX + 2, NY + 2), eta(NX + 2, NY + 2)
      real(wp) :: max_err
      integer :: nx_tot, ny_tot

      call grid%init(NX, NY, 1, 1.0_wp, 1.0_wp)
      call vc%init(grid, nz_ml=NZ)
      vc%coord_type = VCOORD_Z_FIXED
      vc%z_fixed_h_ref = H_REF
      nx_tot = grid%nx_total
      ny_tot = grid%ny_total

      total_h = H_REF
      eta = 0.0_wp
      call vc%compute_target_h(total_h(1:nx_tot, 1:ny_tot), eta(1:nx_tot, 1:ny_tot))

      ! Both layers exactly 500 m
      max_err = max(maxval(abs(vc%target_h(:, :, 1) - 500.0_wp)), &
                    maxval(abs(vc%target_h(:, :, 2) - 500.0_wp)))
      call check(error, max_err < 1.0e-10_wp, &
                 "Z_FIXED deep: both layers should equal h_ref/nz = 500 m")
      call vc%destroy()
   end subroutine test_z_fixed_deep

   subroutine test_z_fixed_shallow(error)
      !! Truly shallow column: H = 400 m < h_nominal = 500 m for the
      !! nz=2 / h_ref=1000 setup.  The bed layer's nominal slot
      !! [-1000, -500] is entirely below the bed (z=-400), so the
      !! algorithm collapses the bed to `h_min` and lets the surface
      !! layer absorb the residual.  Demonstrates the "vanishing bed
      !! in shallow cells" MOM6 gprime semantics.  Sum of targets ≈ H
      !! (column conservation to within h_min round-off).
      !!
      !! At H = 600 m (between h_nominal and 2·h_nominal) the bed
      !! layer fits and takes the residual (100 m); only H < h_nominal
      !! triggers vanishing.  That fitting case is exercised by the
      !! deep_locks_interfaces test below at the boundary.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      integer, parameter :: NZ = 2
      integer, parameter :: NX = 4, NY = 4
      real(wp), parameter :: H_REF = 1000.0_wp
      real(wp), parameter :: H_SHALLOW = 400.0_wp
      real(wp) :: total_h(NX + 2, NY + 2), eta(NX + 2, NY + 2)
      real(wp) :: h_bed_centre, h_surf_centre, h_min_expected, sum_h

      call grid%init(NX, NY, 1, 1.0_wp, 1.0_wp)
      call vc%init(grid, nz_ml=NZ)
      vc%coord_type = VCOORD_Z_FIXED
      vc%z_fixed_h_ref = H_REF
      h_min_expected = vc%zstar_h_min   ! default 1.0e-4_wp

      total_h = H_SHALLOW
      eta = 0.0_wp
      call vc%compute_target_h(total_h(1:grid%nx_total, 1:grid%ny_total), &
                               eta(1:grid%nx_total, 1:grid%ny_total))

      ! Pick a deeply-interior cell to avoid any edge effects
      h_bed_centre = vc%target_h(grid%nghost + 1, grid%nghost + 1, 1)
      h_surf_centre = vc%target_h(grid%nghost + 1, grid%nghost + 1, 2)
      sum_h = h_bed_centre + h_surf_centre

      call check(error, abs(h_bed_centre - h_min_expected) < 1.0e-10_wp, &
                 "Z_FIXED shallow: bed layer should collapse to h_min")
      if (allocated(error)) goto 200
      ! Surface takes (column - h_min) = 400 - 1e-4 ≈ 400
      call check(error, abs(h_surf_centre - (H_SHALLOW - h_min_expected)) < 1.0e-9_wp, &
                 "Z_FIXED shallow: surface absorbs residual (≈ H_SHALLOW)")
      if (allocated(error)) goto 200
      ! Column conservation: sum of targets = original column total
      call check(error, abs(sum_h - H_SHALLOW) < 1.0e-9_wp, &
                 "Z_FIXED shallow: sum(target_h) should equal H (column conservation)")
200   call vc%destroy()
   end subroutine test_z_fixed_shallow

   subroutine test_z_fixed_fallback(error)
      !! With `z_fixed_h_ref = 0` (default, knob unset) the kernel
      !! must fall back to uniform `H · dsig(k)` so tests / callers
      !! that haven't wired the knob through don't get a meaningless
      !! all-vanished column.  Mirrors the sigma path bit-identically
      !! when `dsig` is uniform (which it is post-`init`).
      type(error_type), allocatable, intent(out) :: error
      type(ocean_vcoord_t) :: vc
      type(hgrid_t) :: grid
      integer, parameter :: NZ = 4
      integer, parameter :: NX = 4, NY = 4
      real(wp), parameter :: H = 800.0_wp
      real(wp) :: total_h(NX + 2, NY + 2), eta(NX + 2, NY + 2)
      real(wp) :: max_err
      integer :: i, j, k, nx_tot, ny_tot

      call grid%init(NX, NY, 1, 1.0_wp, 1.0_wp)
      call vc%init(grid, nz_ml=NZ)
      vc%coord_type = VCOORD_Z_FIXED
      ! Deliberately do NOT set z_fixed_h_ref — it should stay at the
      ! default 0.0, triggering the fallback path.
      nx_tot = grid%nx_total
      ny_tot = grid%ny_total

      total_h = H
      eta = 0.0_wp
      call vc%compute_target_h(total_h(1:nx_tot, 1:ny_tot), eta(1:nx_tot, 1:ny_tot))

      ! Expect uniform H/nz everywhere
      max_err = 0.0_wp
      do k = 1, NZ
         do j = 1, ny_tot
            do i = 1, nx_tot
               max_err = max(max_err, abs(vc%target_h(i, j, k) - H/real(NZ, wp)))
            end do
         end do
      end do
      call check(error, max_err < 1.0e-10_wp, &
                 "Z_FIXED with h_ref=0 should fall back to uniform H·dsig")
      call vc%destroy()
   end subroutine test_z_fixed_fallback

end module test_ocean_vcoord
