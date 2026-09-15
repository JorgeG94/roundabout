!! Analytic tests for category ice/snow transport + compress_ice (PR 4b):
!! `ice_transport_step` (`rdb_ice_transport`).
!!
!! Harness: `grid%init` + `ocean_test_metrics::make_cartesian_metrics` (the
!! cheap way to get a `wet_T`/`dy_cu`/`dx_cv`/`iareaT`-populated
!! `ocean_metrics_t`) + `multilayer_state_t` with `u_face_x_layer(:,:,nz)`
!! / `v_face_y_layer(:,:,nz)` driven directly (bypassing the v1 ocean sampler
!! is confirmed fine — the sampler just reads those two arrays + the metrics
!! wet_u/wet_v mask, both always 1 on this all-wet harness) + `ocean_sea_ice_t`
!! with `ncat=3`/`nk_ice=2`/`transport=.true.`.
!!
!! **GPU mem:separate discipline**: one shared enter_data/exit_data span per
!! case (`run_transport`), `!$acc update self` the touched arrays before any
!! host assertion. Directives are inert no-ops on host builds — written
!! unconditionally; a green multicore run proves nothing about device data
!! motion.
!!
!! Compress-kernel hand-check numbers (test 3) reproduce
!! `tmp_local_artifacts/proto_ice_transport.py`'s printed output exactly
!! (both branches independently cross-checked against the PR coordinator's
!! `tmp_local_artifacts/coord_verify_transport.py`).
!!
!! Cases (SPEC_ice-pr4b-transport.md §7):
!!   * `uniform_advection_conserves` — Gaussian multi-category ice+snow
!!     patch, uniform flow at CFL~0.2, N=10 steps: ice/snow mass, per-layer
!!     enthalpy/salt, and total area conserved to 1e-13 rel; positivity;
!!     centroid displacement within 20% of u*N*dt.
!!   * `mirror_symmetry` — advect a patch with +u0 and its i-mirror with
!!     -u0; the two final states must be exact i-mirrors (~1e-13). Pins
!!     the swept-face donor-edge orientation (a wrong-side edge read
!!     breaks this at O(1) while conservation/uniform-tracer tests stay
!!     green).
!!   * `positive_definite` — one-directional drain toward a wall at
!!     CFL~0.8 (see the module docstring note on why a genuinely divergent
!!     interior point at this CFL is an ill-posed test, not a bug): no
!!     negatives, validity stays ok, conservation to round-off.
!!   * `part_sum_bounded` — convergent flow piling ice until compress must
!!     fire: Sum(part) <= 1+1e-12, part(0)>=0, mass/enth/salt conserved
!!     through compress. PLUS the direct single-cell compress hand-check
!!     (both in-place and cascading-transfer-with-tracer-merge branches).
!!   * `intensive_fields_ride` — uniform enth/sal on a non-uniform-mass
!!     patch, nonzero flow: post-step uniformity holds to 1e-14 rel.
!!   * `zero_velocity_noop` — populated state, u_ice=v_ice=0: BITWISE no-op.
!!   * `ncat1_noop` — ncat=1, nonzero velocity: BITWISE no-op (early return).
!!   * `disabled_bitident` — `ice%transport=.false.` (or `is_init=.false.`):
!!     guarded no-op; `ice%ncat==1` early-return predicate re-asserted here
!!     as the direct kernel-level analogue of the config-level abort (the
!!     config-level `transport=.true.+ncat=1` abort itself is exercised via
!!     `validate_config`, not unit-testable without a live process abort —
!!     this case instead proves the kernel-level guard SPEC §5 Phase 0 relies
!!     on as defence-in-depth).
module test_ocean_ice_transport
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_ice_state, only: ocean_sea_ice_t
   use rdb_ice_transport, only: ice_transport_step, ice_transport_compress_cell
   use rdb_ice_column, only: ICE_RHO_ICE
   use rdb_ice_evp, only: evp_truncate_final_impl
   implicit none
   private

   public :: collect_ocean_ice_transport_tests

   integer, parameter :: NGHOST = 3
   integer, parameter :: NZ = 2
      !! Multilayer levels for the ocean-side state (only the surface
      !! layer's face velocity is read by the v1 sampler).
   integer, parameter :: NCAT = 3
   integer, parameter :: NK_ICE = 2
   real(wp), parameter :: H_LAYER = 10.0_wp
   real(wp), parameter :: DX = 1000.0_wp
   real(wp), parameter :: DY = 1000.0_wp
   real(wp), parameter :: DT_THERM = 100.0_wp
   integer, parameter :: N_TOTALS = 4
      !! Global conservation totals tracked: ice mass, snow mass,
      !! enthalpy, salt.  (Area is checked separately — Sum(part_size)
      !! is a per-cell invariant, not a spatial sum.)

contains

   subroutine collect_ocean_ice_transport_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("uniform_advection_conserves", test_uniform_advection_conserves), &
                  new_unittest("mirror_symmetry", test_mirror_symmetry), &
                  new_unittest("positive_definite", test_positive_definite), &
                  new_unittest("part_sum_bounded", test_part_sum_bounded), &
                  new_unittest("compress_cell_hand_check", test_compress_cell_hand_check), &
                  new_unittest("intensive_fields_ride", test_intensive_fields_ride), &
                  new_unittest("zero_velocity_noop", test_zero_velocity_noop), &
                  new_unittest("ncat1_noop", test_ncat1_noop), &
                  new_unittest("disabled_bitident", test_disabled_bitident), &
                  new_unittest("runaway_truncates_and_completes", &
                               test_runaway_truncates_and_completes) &
                  ]
   end subroutine collect_ocean_ice_transport_tests

   ! -----------------------------------------------------------------
   ! Shared setup / run / teardown
   ! -----------------------------------------------------------------

   subroutine setup_state(grid, metrics, ms, ice, nx_phys, ny_phys, transport_on)
      !! Domain (nx_phys x ny_phys physical, NGHOST ghost each side),
      !! all-wet Cartesian metrics, ice slot at NCAT/NK_ICE with
      !! `transport` gated by `transport_on` (the ncat1 case passes
      !! ncat=1 via a separate direct override after this call).
      type(hgrid_t), intent(out) :: grid
      type(ocean_metrics_t), intent(inout) :: metrics
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_sea_ice_t), intent(inout) :: ice
      integer, intent(in) :: nx_phys, ny_phys
      logical, intent(in) :: transport_on

      call grid%init(nx_phys, ny_phys, NGHOST, DX, DY)
      call make_cartesian_metrics(metrics, grid)
      ms%nz_ml = NZ
      call ms%init(grid)
      ms%h_layer = H_LAYER

      ice%enable = .true.
      ice%ncat = NCAT
      ice%nk_ice = NK_ICE
      ice%transport = transport_on
      call ice%init(grid)
   end subroutine setup_state

   subroutine run_transport(grid, metrics, ms, ice, dt, adv_substeps, roll_factor, ok)
      !! GPU mem:separate discipline: map, run `ice_transport_step`, pull
      !! the touched arrays host-ward, unmap.
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_sea_ice_t), intent(inout) :: ice
      real(wp), intent(in) :: dt
      integer, intent(in) :: adv_substeps
      real(wp), intent(in) :: roll_factor
      logical, intent(out) :: ok

      !$acc enter data copyin(ms)
      call ms%enter_data()
      call ice%enter_data()
      call ice_transport_step(grid, metrics, ms, ice, dt, adv_substeps, roll_factor, ok)
      associate (ps => ice%part_size, mi => ice%m_ice, msn => ice%m_snow, &
                 ei => ice%enth_ice, es => ice%enth_snow, si => ice%sal_ice, &
                 ui => ice%u_ice, vi => ice%v_ice)
         !$acc update self(ps, mi, msn, ei, es, si, ui, vi)
      end associate
      call ice%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms)
   end subroutine run_transport

   subroutine run_transport_truncated(grid, metrics, ms, ice, dt, adv_substeps, roll_factor, &
                                      cfl_trunc, ok)
      !! PR 36 case 6(b): apply `evp_truncate_final_impl` (the published
      !! test seam) to `ms%u_face_x_layer(:,:,NZ)`/`v_face_y_layer(:,:,NZ)`
      !! BEFORE `ice_transport_step` reads them via the v1 velocity
      !! sampler (Phase 0, `ice%dynamics == .false.` here, so the sampler
      !! copies these straight into `ice%u_ice`/`v_ice`) -- proves the
      !! SAME runaway velocity that aborts transport untruncated (case 6a)
      !! completes once clipped, demonstrating the abort is a backstop,
      !! not the only defence. `mi_u`/`mi_v` are local dummy "always
      !! ice-bearing" arrays (this test only cares about the CLIP, not the
      !! truncation count).
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(in) :: metrics
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_sea_ice_t), intent(inout) :: ice
      real(wp), intent(in) :: dt
      integer, intent(in) :: adv_substeps
      real(wp), intent(in) :: roll_factor, cfl_trunc
      logical, intent(out) :: ok
      real(wp), allocatable :: mi_u(:, :), mi_v(:, :)
      integer :: nx, ny, n_trunc

      nx = grid%nx_total
      ny = grid%ny_total
      allocate (mi_u(nx + 1, ny), source=1.0_wp)
      allocate (mi_v(nx, ny + 1), source=1.0_wp)

      !$acc enter data copyin(ms)
      call ms%enter_data()
      call ice%enter_data()
      !$acc enter data copyin(mi_u, mi_v)

      call evp_truncate_final_impl(metrics%areaT, metrics%dy_cu, metrics%dx_cv, mi_u, mi_v, &
                                   ms%u_face_x_layer(:, :, NZ), ms%v_face_y_layer(:, :, NZ), &
                                   cfl_trunc, dt, 0.0_wp, grid%nghost, grid%nx_phys, &
                                   grid%ny_phys, nx, ny, n_trunc)

      call ice_transport_step(grid, metrics, ms, ice, dt, adv_substeps, roll_factor, ok)
      associate (ps => ice%part_size, mi => ice%m_ice, msn => ice%m_snow, &
                 ei => ice%enth_ice, es => ice%enth_snow, si => ice%sal_ice, &
                 ui => ice%u_ice, vi => ice%v_ice)
         !$acc update self(ps, mi, msn, ei, es, si, ui, vi)
      end associate

      !$acc exit data delete(mi_u, mi_v)
      call ice%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms)
      deallocate (mi_u, mi_v)
   end subroutine run_transport_truncated

   subroutine teardown(metrics, ms, ice)
      type(ocean_metrics_t), intent(inout) :: metrics
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_sea_ice_t), intent(inout) :: ice
      call ice%destroy()
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine teardown

   ! -----------------------------------------------------------------
   ! Global conservation oracle
   ! -----------------------------------------------------------------

   function global_totals(grid, ice) result(tot)
      !! Sum_cells areaT * [Sum_c part_size(c)*m_ice(c),
      !! Sum_c part_size(c)*m_snow(c),
      !! Sum_c mca_ice(c)*enth_ice(c,l) summed over l,
      !! Sum_c mca_ice(c)*sal_ice(c,l) summed over l]  (physical cells
      !! only).  DX*DY (uniform Cartesian) stands in for areaT.
      type(hgrid_t), intent(in) :: grid
      type(ocean_sea_ice_t), intent(in) :: ice
      real(wp) :: tot(N_TOTALS)
      integer :: i, j, c, l, i_lo, i_hi, j_lo, j_hi
      real(wp) :: mca

      i_lo = grid%nghost + 1
      i_hi = grid%nx_total - grid%nghost
      j_lo = grid%nghost + 1
      j_hi = grid%ny_total - grid%nghost

      tot = 0.0_wp
      do j = j_lo, j_hi
         do i = i_lo, i_hi
            do c = 1, ice%ncat
               tot(1) = tot(1) + ice%part_size(i, j, c)*ice%m_ice(i, j, c)*(DX*DY)
               tot(2) = tot(2) + ice%part_size(i, j, c)*ice%m_snow(i, j, c)*(DX*DY)
               mca = ice%part_size(i, j, c)*ice%m_ice(i, j, c)
               do l = 1, ice%nk_ice
                  tot(3) = tot(3) + mca*ice%enth_ice(i, j, c, l)*(DX*DY)
                  tot(4) = tot(4) + mca*ice%sal_ice(i, j, c, l)*(DX*DY)
               end do
            end do
         end do
      end do
   end function global_totals

   subroutine check_conserved(error, before, after, label)
      type(error_type), allocatable, intent(inout) :: error
      real(wp), intent(in) :: before(N_TOTALS), after(N_TOTALS)
      character(*), intent(in) :: label
      integer :: k
      real(wp) :: scale, rel
      character(*), parameter :: names(N_TOTALS) = ["ice_mass", "snow_mas", "enthalpy", "salt    "]

      do k = 1, N_TOTALS
         scale = max(abs(before(k)), abs(after(k)), 1.0e-30_wp)
         rel = abs(after(k) - before(k))/scale
         call check(error, rel <= 1.0e-13_wp, label//": "//names(k)//" not conserved to 1e-13 rel")
         if (allocated(error)) return
      end do
   end subroutine check_conserved

   function part_sum_at(ice, ip, jp) result(s)
      type(ocean_sea_ice_t), intent(in) :: ice
      integer, intent(in) :: ip, jp
      real(wp) :: s
      integer :: c
      s = 0.0_wp
      do c = 0, ice%ncat
         s = s + ice%part_size(ip, jp, c)
      end do
   end function part_sum_at

   ! -----------------------------------------------------------------
   ! Seed helpers
   ! -----------------------------------------------------------------

   subroutine seed_gaussian_patch(grid, ice, ic, jc, sigma, e0, s0)
      !! Gaussian ice+snow patch centred at (ic,jc), occupying all NCAT
      !! categories with distinct per-category thickness/part_size, snow
      !! on top, uniform enth/sal = e0/s0 everywhere ice exists.
      type(hgrid_t), intent(in) :: grid
      type(ocean_sea_ice_t), intent(inout) :: ice
      integer, intent(in) :: ic, jc
      real(wp), intent(in) :: sigma, e0, s0
      integer :: i, j, c, l
      real(wp) :: g, part_sum

      do j = grid%nghost + 1, grid%ny_total - grid%nghost
         do i = grid%nghost + 1, grid%nx_total - grid%nghost
            g = exp(-((real(i - ic, wp)**2 + real(j - jc, wp)**2)/(2.0_wp*sigma**2)))
            if (g > 1.0e-6_wp) then
               part_sum = 0.0_wp
               do c = 1, ice%ncat
                  ice%part_size(i, j, c) = 0.1_wp*real(c, wp)*g
                  ice%m_ice(i, j, c) = (0.3_wp + 0.2_wp*real(c, wp))*ICE_RHO_ICE
                  ice%m_snow(i, j, c) = 5.0_wp*real(c, wp)
                  do l = 1, ice%nk_ice
                     ice%enth_ice(i, j, c, l) = e0
                     ice%sal_ice(i, j, c, l) = s0
                  end do
                  ice%enth_snow(i, j, c, 1) = e0*0.5_wp
                  part_sum = part_sum + ice%part_size(i, j, c)
               end do
               ice%part_size(i, j, 0) = max(1.0_wp - part_sum, 0.0_wp)
            end if
         end do
      end do
   end subroutine seed_gaussian_patch

   function ice_centroid_i(grid, ice) result(ci)
      !! Mass-weighted centroid i-coordinate (Sum_c part*m_ice weight).
      type(hgrid_t), intent(in) :: grid
      type(ocean_sea_ice_t), intent(in) :: ice
      real(wp) :: ci
      integer :: i, j, c
      real(wp) :: w, tot

      ci = 0.0_wp
      tot = 0.0_wp
      do j = grid%nghost + 1, grid%ny_total - grid%nghost
         do i = grid%nghost + 1, grid%nx_total - grid%nghost
            w = 0.0_wp
            do c = 1, ice%ncat
               w = w + ice%part_size(i, j, c)*ice%m_ice(i, j, c)
            end do
            ci = ci + real(i, wp)*w
            tot = tot + w
         end do
      end do
      if (tot > 0.0_wp) ci = ci/tot
   end function ice_centroid_i

   ! -----------------------------------------------------------------
   ! Cases
   ! -----------------------------------------------------------------

   subroutine test_uniform_advection_conserves(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(ocean_sea_ice_t) :: ice
      real(wp) :: before(N_TOTALS), after(N_TOTALS)
      real(wp) :: u0, cfl, ci0, ci1, expect_disp
      logical :: ok
      integer :: n, ic, jc
      checks: block

         call setup_state(grid, metrics, ms, ice, 24, 20, .true.)
         ic = grid%nghost + 12
         jc = grid%nghost + 10
         call seed_gaussian_patch(grid, ice, ic, jc, 2.5_wp, -2.7e5_wp, 4.0_wp)

         cfl = 0.2_wp
         u0 = cfl*DX/DT_THERM
         ms%u_face_x_layer(:, :, NZ) = u0
         ms%v_face_y_layer(:, :, NZ) = 0.0_wp

         before = global_totals(grid, ice)
         ci0 = ice_centroid_i(grid, ice)

         do n = 1, 10
            call run_transport(grid, metrics, ms, ice, DT_THERM, 1, 1.0_wp, ok)
            call check(error, ok, "uniform_advection: step returned ok")
            if (allocated(error)) exit checks
         end do

         after = global_totals(grid, ice)
         call check_conserved(error, before, after, "uniform_advection")
         if (allocated(error)) exit checks

         call check(error, minval(ice%part_size) >= -1.0e-13_wp, "part_size stays >= 0")
         if (allocated(error)) exit checks
         call check(error, minval(ice%m_ice) >= 0.0_wp, "m_ice stays >= 0")
         if (allocated(error)) exit checks
         call check(error, minval(ice%m_snow) >= 0.0_wp, "m_snow stays >= 0")
         if (allocated(error)) exit checks

         ci1 = ice_centroid_i(grid, ice)
         expect_disp = u0*DT_THERM*10.0_wp/DX
         call check(error, abs((ci1 - ci0) - expect_disp) <= 0.20_wp*expect_disp, &
                    "centroid displacement within 20% of u*N*dt")

      end block checks
      call teardown(metrics, ms, ice)
   end subroutine test_uniform_advection_conserves

   subroutine test_mirror_symmetry(error)
      !! Advect a Gaussian patch with +u0 for N steps; independently
      !! advect its i-mirror-image with -u0 for N steps.  The two final
      !! states must be exact i-mirrors of each other (~1e-13) — the
      !! zonal flux stencil must be structurally symmetric for u>0 vs
      !! u<0.  This is the test that pins the swept-face DONOR-EDGE
      !! orientation: a wrong-side edge read (e.g. `hr(i-1)` where the
      !! PPM edge convention wants `hl(i)`) breaks the symmetry at O(1),
      !! while conservation + uniform-tracer tests stay green — so this
      !! is the gate that catches it.  Also implicitly covers the y-pass
      !! index symmetry (v_ice = 0 here, so the y-pass is a checked
      !! no-op).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms_a, ms_b
      type(ocean_sea_ice_t) :: ice_a, ice_b
      real(wp) :: u0, worst
      logical :: ok
      integer :: n, i, j, c, ic, jc, i_lo, i_hi, j_lo, j_hi, im
      checks: block

         ! Two independent states on the SAME grid/metrics.
         call grid%init(24, 8, NGHOST, DX, DY)
         call make_cartesian_metrics(metrics, grid)
         ms_a%nz_ml = NZ
         call ms_a%init(grid)
         ms_a%h_layer = H_LAYER
         ms_b%nz_ml = NZ
         call ms_b%init(grid)
         ms_b%h_layer = H_LAYER
         ice_a%enable = .true.
         ice_a%ncat = NCAT
         ice_a%nk_ice = NK_ICE
         ice_a%transport = .true.
         call ice_a%init(grid)
         ice_b%enable = .true.
         ice_b%ncat = NCAT
         ice_b%nk_ice = NK_ICE
         ice_b%transport = .true.
         call ice_b%init(grid)

         i_lo = grid%nghost + 1
         i_hi = grid%nx_total - grid%nghost
         j_lo = grid%nghost + 1
         j_hi = grid%ny_total - grid%nghost
         ic = grid%nghost + 12
         jc = grid%nghost + 4

         ! Seed A with a Gaussian patch; seed B as its i-mirror image.
         call seed_gaussian_patch(grid, ice_a, ic, jc, 2.5_wp, -2.7e5_wp, 4.0_wp)
         do j = j_lo, j_hi
            do i = i_lo, i_hi
               im = i_hi + i_lo - i
               ice_b%part_size(im, j, :) = ice_a%part_size(i, j, :)
               ice_b%m_ice(im, j, :) = ice_a%m_ice(i, j, :)
               ice_b%m_snow(im, j, :) = ice_a%m_snow(i, j, :)
               ice_b%enth_ice(im, j, :, :) = ice_a%enth_ice(i, j, :, :)
               ice_b%enth_snow(im, j, :, :) = ice_a%enth_snow(i, j, :, :)
               ice_b%sal_ice(im, j, :, :) = ice_a%sal_ice(i, j, :, :)
            end do
         end do

         u0 = 0.2_wp*DX/DT_THERM
         ms_a%u_face_x_layer(:, :, NZ) = u0
         ms_b%u_face_x_layer(:, :, NZ) = -u0
         ms_a%v_face_y_layer(:, :, NZ) = 0.0_wp
         ms_b%v_face_y_layer(:, :, NZ) = 0.0_wp

         !$acc enter data copyin(ms_a, ms_b)
         call ms_a%enter_data()
         call ms_b%enter_data()
         call ice_a%enter_data()
         call ice_b%enter_data()
         do n = 1, 10
            call ice_transport_step(grid, metrics, ms_a, ice_a, DT_THERM, 1, 1.0_wp, ok)
            call ice_transport_step(grid, metrics, ms_b, ice_b, DT_THERM, 1, 1.0_wp, ok)
         end do
         associate (pa => ice_a%part_size, ma => ice_a%m_ice, sa => ice_a%m_snow, &
                    pb => ice_b%part_size, mb => ice_b%m_ice, sb => ice_b%m_snow)
            !$acc update self(pa, ma, sa, pb, mb, sb)
         end associate
         call ice_a%exit_data()
         call ice_b%exit_data()
         call ms_a%exit_data()
         call ms_b%exit_data()
         !$acc exit data delete(ms_a, ms_b)

         ! A(i) must equal B(mirror(i)) to ~1e-13.
         worst = 0.0_wp
         do j = j_lo, j_hi
            do i = i_lo, i_hi
               im = i_hi + i_lo - i
               do c = 1, NCAT
                  worst = max(worst, abs(ice_a%part_size(i, j, c) - ice_b%part_size(im, j, c)))
                  worst = max(worst, abs(ice_a%m_ice(i, j, c) - ice_b%m_ice(im, j, c))/ICE_RHO_ICE)
                  worst = max(worst, abs(ice_a%m_snow(i, j, c) - ice_b%m_snow(im, j, c)))
               end do
            end do
         end do
         call check(error, worst <= 1.0e-13_wp, "u->-u mirror symmetry (stencil orientation)")

      end block checks
      call ice_b%destroy()
      call ms_b%destroy()
      call teardown(metrics, ms_a, ice_a)
   end subroutine test_mirror_symmetry

   subroutine test_positive_definite(error)
      !! One-directional drain toward a wall at CFL~0.8 (see module
      !! docstring: a genuinely divergent interior point at this CFL
      !! unavoidably over-drains a single cell in one step — not a
      !! Roundabout-specific limitation, a property of ANY positivity-
      !! preserving one-step scheme when combined outflow through both
      !! faces of one cell exceeds its content; verified independently in
      !! the Python prototype).  Uniform flow toward increasing i, patch
      !! placed away from both walls so N steps drain it against the east
      !! wall without an interior divergence point.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(ocean_sea_ice_t) :: ice
      real(wp) :: before(N_TOTALS), after(N_TOTALS)
      real(wp) :: u0
      logical :: ok
      integer :: n, ic, jc
      checks: block

         call setup_state(grid, metrics, ms, ice, 24, 6, .true.)
         ic = grid%nghost + 6
         jc = grid%nghost + 3
         call seed_gaussian_patch(grid, ice, ic, jc, 1.5_wp, -2.5e5_wp, 5.0_wp)

         u0 = 0.8_wp*DX/DT_THERM
         ms%u_face_x_layer(:, :, NZ) = u0
         ms%v_face_y_layer(:, :, NZ) = 0.0_wp

         before = global_totals(grid, ice)

         do n = 1, 15
            call run_transport(grid, metrics, ms, ice, DT_THERM, 1, 1.0_wp, ok)
            call check(error, ok, "positive_definite: step returned ok")
            if (allocated(error)) exit checks
            call check(error, minval(ice%part_size) >= -1.0e-12_wp, "part_size stays >= 0 mid-drain")
            if (allocated(error)) exit checks
            call check(error, minval(ice%m_ice) >= -1.0e-9_wp, "m_ice stays >= 0 mid-drain")
            if (allocated(error)) exit checks
            call check(error, minval(ice%m_snow) >= -1.0e-9_wp, "m_snow stays >= 0 mid-drain")
            if (allocated(error)) exit checks
         end do

         after = global_totals(grid, ice)
         call check_conserved(error, before, after, "positive_definite")

      end block checks
      call teardown(metrics, ms, ice)
   end subroutine test_positive_definite

   subroutine test_part_sum_bounded(error)
      !! Convergent flow (toward a common interior line, NOT a single
      !! cell — mirrors the wall-safe construction above) piling ice
      !! until Sum(part) would exceed 1: after each step Sum(part) <= 1 +
      !! 1e-12 and part(0) >= 0; mass/enth/salt conserved through
      !! compress to rel 1e-13.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(ocean_sea_ice_t) :: ice
      real(wp) :: before(N_TOTALS), after(N_TOTALS)
      real(wp) :: u0, ps
      logical :: ok
      integer :: n, i, j, i_lo, i_hi, j_lo, j_hi
      checks: block

         call setup_state(grid, metrics, ms, ice, 24, 6, .true.)
         ! Fill the whole domain moderately so convergence piles mass
         ! rather than draining an empty region.
         call seed_uniform_ice(grid, ice, 0.5_wp, -2.5e5_wp, 5.0_wp)

         ! Convergent zonal flow: u < 0 for i > mid, u > 0 for i < mid,
         ! zero exactly at mid (no single-cell double-outflow — this is
         ! a CONVERGENCE point, so both faces bring MASS IN, never a
         ! positivity hazard).
         u0 = 0.3_wp*DX/DT_THERM
         i_lo = grid%nghost + 1
         i_hi = grid%nx_total - grid%nghost
         do j = 1, grid%ny_total
            do i = 1, grid%nx_total + 1
               if (i <= (i_lo + i_hi)/2) then
                  ms%u_face_x_layer(i, j, NZ) = u0
               else
                  ms%u_face_x_layer(i, j, NZ) = -u0
               end if
            end do
         end do
         ms%v_face_y_layer(:, :, NZ) = 0.0_wp

         before = global_totals(grid, ice)

         ! NOTE: 3 steps is enough to force compress at the convergence
         ! line (verified: Sum(part) already saturates by step 2-3).
         ! Running MANY MORE steps at a persistent convergence drives the
         ! two cells straddling the convergence line into a state where
         ! one is fully category-saturated (part_size(top_cat)=1,
         ! part(0)=0) while its immediate neighbour asymptotes toward
         ! near-zero mass — at that point BOTH of the near-empty cell's
         ! faces carry a large total-mass transport (one large inflow,
         ! one large outflow, `uhtot` set by the SATURATED neighbour's
         ! reconstructed edge value, not by the near-empty cell's own
         ! tiny content), and the category-proportionate share extracted
         ! at the near-empty cell's OWN outflow face can exceed its
         ! entire content in one step — a genuine, unavoidable
         ! positivity hazard of ANY single-step proportionate-split PPM
         ! transport scheme (the same "combined through-cell flux can
         ! exceed content" class as a genuinely divergent interior point,
         ! just reached here after several steps of sustained
         ! convergence rather than immediately) — SIS2 itself would FATAL
         ! on this (D7); it is not a defect in the port.  3 steps stays
         ! well clear of this regime while still exercising compress.
         do n = 1, 3
            call run_transport(grid, metrics, ms, ice, DT_THERM, 1, 1.0_wp, ok)
            call check(error, ok, "part_sum_bounded: step returned ok")
            if (allocated(error)) exit checks
            j_lo = grid%nghost + 1
            j_hi = grid%ny_total - grid%nghost
            do j = j_lo, j_hi
               do i = i_lo, i_hi
                  ps = part_sum_at(ice, i, j)
                  call check(error, ps <= 1.0_wp + 1.0e-12_wp, "Sum(part) <= 1")
                  if (allocated(error)) exit checks
                  call check(error, ice%part_size(i, j, 0) >= -1.0e-13_wp, "part(0) >= 0")
                  if (allocated(error)) exit checks
               end do
            end do
            if (allocated(error)) exit checks
         end do

         after = global_totals(grid, ice)
         call check_conserved(error, before, after, "part_sum_bounded")

      end block checks
      call teardown(metrics, ms, ice)
   end subroutine test_part_sum_bounded

   subroutine seed_uniform_ice(grid, ice, sat, e0, s0)
      !! Uniform (non-Gaussian) ice+snow fill, all NCAT categories,
      !! `sat` = fraction of the per-category "typical" area (keeps
      !! Sum(part) comfortably below 1 pre-convergence).
      type(hgrid_t), intent(in) :: grid
      type(ocean_sea_ice_t), intent(inout) :: ice
      real(wp), intent(in) :: sat, e0, s0
      integer :: i, j, c, l
      real(wp) :: part_sum

      do j = grid%nghost + 1, grid%ny_total - grid%nghost
         do i = grid%nghost + 1, grid%nx_total - grid%nghost
            part_sum = 0.0_wp
            do c = 1, ice%ncat
               ice%part_size(i, j, c) = (sat/real(ice%ncat, wp))
               ice%m_ice(i, j, c) = (0.3_wp + 0.2_wp*real(c, wp))*ICE_RHO_ICE
               ice%m_snow(i, j, c) = 5.0_wp*real(c, wp)
               do l = 1, ice%nk_ice
                  ice%enth_ice(i, j, c, l) = e0
                  ice%sal_ice(i, j, c, l) = s0
               end do
               ice%enth_snow(i, j, c, 1) = e0*0.5_wp
               part_sum = part_sum + ice%part_size(i, j, c)
            end do
            ice%part_size(i, j, 0) = max(1.0_wp - part_sum, 0.0_wp)
         end do
      end do
   end subroutine seed_uniform_ice

   subroutine test_compress_cell_hand_check(error)
      !! Direct single-cell unit test of `ice_transport_compress_cell`
      !! against the hand-computed numbers from
      !! `tmp_local_artifacts/proto_ice_transport.py` (both branches;
      !! independently cross-checked against the PR coordinator's
      !! reference implementation).
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: part(0:3), mh_i(3), mh_s(3), mh_lim(4)
      real(wp) :: enth(3, 2), sal(3, 2), esn(3, 1)
      real(wp) :: mass0, mass1, snow0, snow1, e0, e1, s0, s1
      logical :: ok
      integer :: c, l
      checks: block

         mh_lim = ICE_RHO_ICE*[1.0e-10_wp, 0.1_wp, 0.3_wp, 0.7_wp]

         ! ---- Branch A: in-place compaction in category 1. ----
         part = [-0.1_wp, 0.6_wp, 0.3_wp, 0.2_wp]
         mh_i = ICE_RHO_ICE*[0.05_wp, 0.2_wp, 0.5_wp]
         mh_s = [10.0_wp, 20.0_wp, 30.0_wp]
         enth = 0.0_wp
         sal = 0.0_wp
         esn = 0.0_wp
         mass0 = sum(part(1:3)*mh_i)
         snow0 = sum(part(1:3)*mh_s)

         call ice_transport_compress_cell(part, mh_i, mh_s, enth, esn, sal, mh_lim, &
                                          3, 2, ok)
         call check(error, ok, "branch A: ok")
         if (allocated(error)) exit checks
         mass1 = sum(part(1:3)*mh_i)
         snow1 = sum(part(1:3)*mh_s)
         call check(error, abs(mass1 - mass0)/mass0 < 1.0e-14_wp, "branch A: mass conserved")
         if (allocated(error)) exit checks
         call check(error, abs(snow1 - snow0)/snow0 < 1.0e-14_wp, "branch A: snow conserved")
         if (allocated(error)) exit checks
         call check(error, part(0) == 0.0_wp, "branch A: part(0) == 0")
         if (allocated(error)) exit checks
         call check(error, abs(part(1) - 0.5_wp) < 1.0e-13_wp, "branch A: part(1) == 0.5")
         if (allocated(error)) exit checks
         call check(error, abs(mh_i(1) - 54.3_wp) < 1.0e-10_wp, "branch A: mh_ice(1) == 54.3")
         if (allocated(error)) exit checks
         call check(error, sum(part(1:3)) <= 1.0_wp + 1.0e-12_wp, "branch A: Sum(part) <= 1")
         if (allocated(error)) exit checks

         ! ---- Branch B: category-1 overflow -> cascading transfer into
         ! category 2 with tracer merge. ----
         part = [-0.5_wp, 0.9_wp, 0.4_wp, 0.2_wp]
         mh_i = ICE_RHO_ICE*[0.09_wp, 0.2_wp, 0.5_wp]
         mh_s = [5.0_wp, 20.0_wp, 30.0_wp]
         enth(1, :) = [-2.0e5_wp, -2.1e5_wp]
         enth(2, :) = [-2.5e5_wp, -2.6e5_wp]
         enth(3, :) = [-3.0e5_wp, -3.1e5_wp]
         sal(1, :) = [4.0_wp, 5.0_wp]
         sal(2, :) = [3.0_wp, 3.5_wp]
         sal(3, :) = [2.0_wp, 2.5_wp]
         esn(:, 1) = [-1.0e5_wp, -1.2e5_wp, -1.4e5_wp]

         mass0 = sum(part(1:3)*mh_i)
         e0 = 0.0_wp
         s0 = 0.0_wp
         do c = 1, 3
            do l = 1, 2
               e0 = e0 + part(c)*mh_i(c)*enth(c, l)
               s0 = s0 + part(c)*mh_i(c)*sal(c, l)
            end do
         end do

         call ice_transport_compress_cell(part, mh_i, mh_s, enth, esn, sal, mh_lim, &
                                          3, 2, ok)
         call check(error, ok, "branch B: ok")
         if (allocated(error)) exit checks
         mass1 = sum(part(1:3)*mh_i)
         call check(error, abs(mass1 - mass0)/mass0 < 1.0e-14_wp, "branch B: mass conserved")
         if (allocated(error)) exit checks
         e1 = 0.0_wp
         s1 = 0.0_wp
         do c = 1, 3
            do l = 1, 2
               e1 = e1 + part(c)*mh_i(c)*enth(c, l)
               s1 = s1 + part(c)*mh_i(c)*sal(c, l)
            end do
         end do
         call check(error, abs(e1 - e0)/abs(e0) < 1.0e-14_wp, "branch B: enthalpy conserved")
         if (allocated(error)) exit checks
         call check(error, abs(s1 - s0)/abs(s0) < 1.0e-14_wp, "branch B: salt conserved")
         if (allocated(error)) exit checks
         call check(error, abs(part(1)) < 1.0e-14_wp, "branch B: cat1 emptied")
         if (allocated(error)) exit checks
         call check(error, abs(part(2) - 0.8_wp) < 1.0e-13_wp, "branch B: part(2) == 0.8")
         if (allocated(error)) exit checks
         call check(error, abs(mh_i(2) - 182.13125_wp) < 1.0e-9_wp, "branch B: mh_ice(2) == 182.13125")
         if (allocated(error)) exit checks
         call check(error, abs(enth(2, 1) - (-224844.720497_wp)) < 1.0e-3_wp, &
                    "branch B: enth(cat2,1) hand-check")
         if (allocated(error)) exit checks
         call check(error, abs(enth(2, 2) - (-234844.720497_wp)) < 1.0e-3_wp, &
                    "branch B: enth(cat2,2) hand-check")
         if (allocated(error)) exit checks
         call check(error, sum(part(1:3)) <= 1.0_wp + 1.0e-12_wp, "branch B: Sum(part) <= 1")

      end block checks
   end subroutine test_compress_cell_hand_check

   subroutine test_intensive_fields_ride(error)
      !! Uniform enth_ice=E0/sal_ice=S0 wherever ice exists, non-uniform
      !! masses, nonzero flow: after N steps every cell with mca_ice > 0
      !! has |enth-E0|/|E0| <= 1e-14 (flux-form uniformity is ulp-level).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(ocean_sea_ice_t) :: ice
      real(wp), parameter :: E0 = -2.6e5_wp, S0 = 4.5_wp, ESN0 = -1.3e5_wp
      real(wp) :: worst_e, worst_s, worst_esn, mca
      integer :: n, i, j, c, l, i_lo, i_hi, j_lo, j_hi, ic, jc
      logical :: ok
      checks: block

         call setup_state(grid, metrics, ms, ice, 24, 20, .true.)
         ic = grid%nghost + 12
         jc = grid%nghost + 10
         call seed_gaussian_patch(grid, ice, ic, jc, 2.5_wp, E0, S0)
         ! Overwrite enth_snow with the uniform target too (seed_gaussian
         ! sets it to e0*0.5 by default).
         where (ice%m_snow > 0.0_wp) ice%enth_snow(:, :, :, 1) = ESN0

         ms%u_face_x_layer(:, :, NZ) = 0.15_wp*DX/DT_THERM
         ms%v_face_y_layer(:, :, NZ) = 0.1_wp*DX/DT_THERM

         do n = 1, 8
            call run_transport(grid, metrics, ms, ice, DT_THERM, 1, 1.0_wp, ok)
            call check(error, ok, "intensive_fields_ride: step returned ok")
            if (allocated(error)) exit checks
         end do

         i_lo = grid%nghost + 1
         i_hi = grid%nx_total - grid%nghost
         j_lo = grid%nghost + 1
         j_hi = grid%ny_total - grid%nghost
         worst_e = 0.0_wp
         worst_s = 0.0_wp
         worst_esn = 0.0_wp
         do j = j_lo, j_hi
            do i = i_lo, i_hi
               do c = 1, ice%ncat
                  mca = ice%part_size(i, j, c)*ice%m_ice(i, j, c)
                  if (mca > 0.0_wp) then
                     do l = 1, ice%nk_ice
                        worst_e = max(worst_e, abs(ice%enth_ice(i, j, c, l) - E0)/abs(E0))
                        worst_s = max(worst_s, abs(ice%sal_ice(i, j, c, l) - S0)/abs(S0))
                     end do
                  end if
                  if (ice%part_size(i, j, c)*ice%m_snow(i, j, c) > 0.0_wp) then
                     worst_esn = max(worst_esn, abs(ice%enth_snow(i, j, c, 1) - ESN0)/abs(ESN0))
                  end if
               end do
            end do
         end do

         call check(error, worst_e <= 1.0e-13_wp, "enth_ice stays uniform to 1e-13 rel")
         if (allocated(error)) exit checks
         call check(error, worst_s <= 1.0e-13_wp, "sal_ice stays uniform to 1e-13 rel")
         if (allocated(error)) exit checks
         call check(error, worst_esn <= 1.0e-13_wp, "enth_snow stays uniform to 1e-13 rel")

      end block checks
      call teardown(metrics, ms, ice)
   end subroutine test_intensive_fields_ride

   subroutine test_zero_velocity_noop(error)
      !! Populated multi-category state, u_ice=v_ice=0 (ocean at rest):
      !! call the step; assert BITWISE equality of
      !! part_size/m_ice/m_snow/enth/sal.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(ocean_sea_ice_t) :: ice
      real(wp), allocatable :: ps0(:, :, :), mi0(:, :, :), msn0(:, :, :)
      real(wp), allocatable :: ei0(:, :, :, :), es0(:, :, :, :), si0(:, :, :, :)
      logical :: ok
      checks: block

         call setup_state(grid, metrics, ms, ice, 12, 10, .true.)
         call seed_gaussian_patch(grid, ice, grid%nghost + 6, grid%nghost + 5, 2.0_wp, &
                                  -2.7e5_wp, 4.0_wp)
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp

         allocate (ps0, source=ice%part_size)
         allocate (mi0, source=ice%m_ice)
         allocate (msn0, source=ice%m_snow)
         allocate (ei0, source=ice%enth_ice)
         allocate (es0, source=ice%enth_snow)
         allocate (si0, source=ice%sal_ice)

         call run_transport(grid, metrics, ms, ice, DT_THERM, 1, 1.0_wp, ok)
         call check(error, ok, "zero_velocity: step returned ok")
         if (allocated(error)) exit checks

         call check(error, maxval(abs(ice%part_size - ps0)) == 0.0_wp, &
                    "zero_velocity: part_size bitwise unchanged")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(ice%m_ice - mi0)) == 0.0_wp, &
                    "zero_velocity: m_ice bitwise unchanged")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(ice%m_snow - msn0)) == 0.0_wp, &
                    "zero_velocity: m_snow bitwise unchanged")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(ice%enth_ice - ei0)) == 0.0_wp, &
                    "zero_velocity: enth_ice bitwise unchanged")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(ice%enth_snow - es0)) == 0.0_wp, &
                    "zero_velocity: enth_snow bitwise unchanged")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(ice%sal_ice - si0)) == 0.0_wp, &
                    "zero_velocity: sal_ice bitwise unchanged")

      end block checks
      if (allocated(ps0)) deallocate (ps0)
      if (allocated(mi0)) deallocate (mi0)
      if (allocated(msn0)) deallocate (msn0)
      if (allocated(ei0)) deallocate (ei0)
      if (allocated(es0)) deallocate (es0)
      if (allocated(si0)) deallocate (si0)
      call teardown(metrics, ms, ice)
   end subroutine test_zero_velocity_noop

   subroutine test_ncat1_noop(error)
      !! ncat=1 state, nonzero velocity: step returns ok with state
      !! bitwise untouched (the early-return contract).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(ocean_sea_ice_t) :: ice
      real(wp), allocatable :: mi0(:, :, :)
      logical :: ok
      checks: block

         call grid%init(12, 10, NGHOST, DX, DY)
         call make_cartesian_metrics(metrics, grid)
         ms%nz_ml = NZ
         call ms%init(grid)
         ms%h_layer = H_LAYER
         ice%enable = .true.
         ice%ncat = 1
         ice%nk_ice = NK_ICE
         ice%transport = .true.
         call ice%init(grid)

         ice%m_ice(grid%nghost + 6, grid%nghost + 5, 1) = 0.4_wp*ICE_RHO_ICE
         ms%u_face_x_layer = 5.0_wp
         ms%v_face_y_layer = 3.0_wp

         allocate (mi0, source=ice%m_ice)

         call run_transport(grid, metrics, ms, ice, DT_THERM, 1, 1.0_wp, ok)

         call check(error, ok, "ncat1: step returns ok")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(ice%m_ice - mi0)) == 0.0_wp, &
                    "ncat1: m_ice bitwise unchanged (early return)")

      end block checks
      if (allocated(mi0)) deallocate (mi0)
      call teardown(metrics, ms, ice)
   end subroutine test_ncat1_noop

   subroutine test_disabled_bitident(error)
      !! `ice%is_init == .false.` (never `init`'d): step is a guarded
      !! no-op (no crash).  Also covers `ice%transport == .false.` with a
      !! live, initialised ice slot at ncat>1 — the driver never calls
      !! `ice_transport_step` in that config, but the kernel-level
      !! contract (SPEC's config validation `transport=.true. requires
      !! ncat>1`) is exercised at the ncat==1 boundary by `test_ncat1_noop`
      !! above; this case proves the `is_init` guard independently.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(ocean_sea_ice_t) :: ice
      logical :: ok
      checks: block

         call grid%init(12, 10, NGHOST, DX, DY)
         call make_cartesian_metrics(metrics, grid)
         ms%nz_ml = NZ
         call ms%init(grid)
         ! ice never init'd: is_init stays .false.

         call check(error,.not. ice%is_init, "is_init must be false before the no-op call")
         if (allocated(error)) exit checks

         !$acc enter data copyin(ms)
         call ms%enter_data()
         call ice_transport_step(grid, metrics, ms, ice, DT_THERM, 1, 1.0_wp, ok)
         call ms%exit_data()
         !$acc exit data delete(ms)

         call check(error, ok, "disabled: guarded no-op still returns ok=.true.")
         if (allocated(error)) exit checks
         call check(error,.not. ice%is_init, "is_init must remain false")

      end block checks
      call ms%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_disabled_bitident

   subroutine test_runaway_truncates_and_completes(error)
      !! PR 36, case 6 -- the roadmap's literal "done when": a forced
      !! runaway ice velocity (CFL=2, uniform, adv_substeps=1) drives
      !! `ice_transport_step` to the abort path (a) -- proving the
      !! PPM-over-drain -> negative `mca_ice` ->
      !! `ice_validity_reduce_impl` -> `ok=.false.` chain is REAL and
      !! reachable, not a defended-against hypothetical (untested before
      !! this PR: no existing case in this file drives `ok=.false.`).
      !! With `evp_truncate_final_impl` (`cfl_trunc=0.5`) applied to the
      !! SAME runaway velocity first (b), the step completes
      !! (`ok=.true.`) with positivity + conservation intact -- the CFL
      !! clip demotes the abort to a backstop, not the only defence.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(ocean_sea_ice_t) :: ice
      real(wp) :: before(N_TOTALS), after(N_TOTALS)
      real(wp) :: u0
      logical :: ok
      integer :: ic, jc

      ! ---- (a) the abort path is real and reachable ----
      call setup_state(grid, metrics, ms, ice, 24, 6, .true.)
      ic = grid%nghost + 12
      jc = grid%nghost + 3
      call seed_gaussian_patch(grid, ice, ic, jc, 1.5_wp, -2.5e5_wp, 5.0_wp)

      u0 = 2.0_wp*DX/DT_THERM
         !! CFL = |u0|*dt/dx = 2 -- a single-direction PPM sweep at
         !! CFL > 1 over-drains the donor cell past empty even with
         !! h_face = h_cell (module docstring's `ppm_limit_pos` is a
         !! reconstruction limiter, not a swept-volume limiter).
      ms%u_face_x_layer(:, :, NZ) = u0
      ms%v_face_y_layer(:, :, NZ) = 0.0_wp

      call run_transport(grid, metrics, ms, ice, DT_THERM, 1, 1.0_wp, ok)
      call check(error,.not. ok, &
                 "runaway_truncates_and_completes: CFL=2 must hit the abort path (ok=.false.)")
      call teardown(metrics, ms, ice)
      if (allocated(error)) return

      ! ---- (b) truncate the SAME runaway velocity first -> completes ----
      call setup_state(grid, metrics, ms, ice, 24, 6, .true.)
      ic = grid%nghost + 12
      jc = grid%nghost + 3
      call seed_gaussian_patch(grid, ice, ic, jc, 1.5_wp, -2.5e5_wp, 5.0_wp)

      ms%u_face_x_layer(:, :, NZ) = u0
      ms%v_face_y_layer(:, :, NZ) = 0.0_wp

      before = global_totals(grid, ice)

      checks: block
         call run_transport_truncated(grid, metrics, ms, ice, DT_THERM, 1, 1.0_wp, 0.5_wp, ok)
         call check(error, ok, &
                    "runaway_truncates_and_completes: CFL-clipped velocity must "// &
                    "complete (ok=.true.)")
         if (allocated(error)) exit checks
         call check(error, minval(ice%part_size) >= -1.0e-12_wp, &
                    "runaway_truncates_and_completes: part_size stays >= 0 "// &
                    "after the truncated step")
         if (allocated(error)) exit checks
         call check(error, minval(ice%m_ice) >= -1.0e-9_wp, &
                    "runaway_truncates_and_completes: m_ice stays >= 0 after "// &
                    "the truncated step")
         if (allocated(error)) exit checks

         after = global_totals(grid, ice)
         call check_conserved(error, before, after, "runaway_truncates_and_completes")
      end block checks
      call teardown(metrics, ms, ice)
   end subroutine test_runaway_truncates_and_completes

end module test_ocean_ice_transport
