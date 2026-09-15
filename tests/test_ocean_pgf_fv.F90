!! Unit tests for the FV_LITE variant of the ocean pressure-force
!! kernel (rdb_ocean_pressure_force with
!! `variant = OPGF_VARIANT_FV_LITE`).
!!
!! FV_LITE adds the z-position correction
!!   PGF_x = -(1/rho0)*[(p_R - p_L)/dx + g*rho_face*(z_R - z_L)/dx]
!! that differences pressure at constant *z* rather than at constant
!! layer index.  The correction cancels the spurious bottom-current
!! PGF that a bare layer-centre pressure difference incurs when layer
!! centres sit at different z-values in adjacent columns (sigma-style
!! coordinates, variable bathymetry, or — as exercised here —
!! artificially redistributed layer thicknesses with uniform total H).
!!
!! `mont` is used as the CROSS-CHECK here, not as the broken foil: the
!! Montgomery potential carries the geopotential in `M` and is a second,
!! independent route to the same acceleration.  The broken arithmetic
!! these cases measure against is written out explicitly below
!! (`host_no_geopotential_dpdx`) so it cannot quietly get fixed.
!!
!! Cases:
!!   * Reduces to Mont on aligned columns — uniform h_layer across
!!     cells, layered density.  The z-correction term must vanish
!!     and FV_LITE's PGF must equal Mont's to round-off.
!!   * Cancels the spurious uniform-density PGF — uniform rho_layer
!!     across the entire 3D field, but h_layer varies per column
!!     (totals match per column so total H is uniform).  Differencing
!!     layer-centre pressure with NO geopotential term reports a large
!!     spurious PGF in the interior layers; FV_LITE and `mont` must
!!     both report zero (to round-off).
!!   * Same field, density-gradient response unchanged — set a true
!!     horizontal density gradient in one layer.  Both Mont and
!!     FV_LITE must capture it; FV is allowed to differ from Mont
!!     by a small amount but not flip sign or vanish.
module test_ocean_pgf_fv
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_ocean_pressure_force, only: ocean_pressure_force_t, &
                                       ocean_pressure_force_compute, &
                                       OPGF_VARIANT_MONT, OPGF_VARIANT_FV_LITE
   implicit none
   private

   public :: collect_ocean_pgf_fv_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 3

contains

   subroutine collect_ocean_pgf_fv_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("fv_matches_mont_on_aligned_columns", test_aligned_match), &
                  new_unittest("fv_cancels_uniform_density_misalignment", test_cancel_spurious), &
                  new_unittest("fv_keeps_real_density_gradient_response", test_real_gradient) &
                  ]
   end subroutine collect_ocean_pgf_fv_tests

   subroutine make_grid(grid, nx_phys, ny_phys)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      call grid%init(nx_phys, ny_phys, NGHOST, 1.0_wp, 1.0_wp)
   end subroutine make_grid

   pure subroutine host_no_geopotential_dpdx(h_layer, rho_layer, dx, dpdx)
      !! The arithmetic FV_LITE's z-correction exists to repair: layer-centre
      !! pressure differenced horizontally along the coordinate surface, with
      !! no geopotential term at all.  Kept here as an explicit reference so
      !! the cancellation cases below have something that demonstrably FAILS
      !! to compare against — previously that role was played by `mont`, which
      !! now carries a real Montgomery potential and cancels correctly too.
      real(wp), intent(in) :: h_layer(:, :, :), rho_layer(:, :, :)
      real(wp), intent(in) :: dx
      real(wp), intent(out) :: dpdx(:, :, :)
      real(wp) :: p_edge(size(h_layer, 1), size(h_layer, 2), NZ + 1)
      real(wp) :: p_c_l, p_c_r
      integer :: i, j, k, nx, ny
      nx = size(h_layer, 1)
      ny = size(h_layer, 2)
      do j = 1, ny
         do i = 1, nx
            p_edge(i, j, NZ + 1) = 0.0_wp
            do k = NZ, 1, -1
               p_edge(i, j, k) = p_edge(i, j, k + 1) + &
                                 GRAVITY*rho_layer(i, j, k)*h_layer(i, j, k)
            end do
         end do
      end do
      dpdx = 0.0_wp
      do k = 1, NZ
         do j = 1, ny
            do i = 2, nx
               p_c_r = 0.5_wp*(p_edge(i, j, k) + p_edge(i, j, k + 1))
               p_c_l = 0.5_wp*(p_edge(i - 1, j, k) + p_edge(i - 1, j, k + 1))
               dpdx(i, j, k) = -(1.0_wp/1035.0_wp)*(p_c_r - p_c_l)/dx
            end do
         end do
      end do
   end subroutine host_no_geopotential_dpdx

   subroutine run_pgf(ms, pgf)
      !! Runs one compute pass and pulls the dpdx/dpdy buffers back
      !! to the host explicitly (they're scratch — create/delete by
      !! default — so a `!$acc update self` before exit_data is
      !! required for the host-side comparison).
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_pressure_force_t), intent(inout) :: pgf
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      call grid%init(size(ms%h_layer, 1) - 2*NGHOST, &
                     size(ms%h_layer, 2) - 2*NGHOST, &
                     NGHOST, 1.0_wp, 1.0_wp)
      call make_cartesian_metrics(metrics, grid)
      !$acc enter data copyin(ms, pgf)
      call ms%enter_data()
      call pgf%enter_data()
      call ocean_pressure_force_compute(grid, metrics, pgf, ms)
      !$acc update self(pgf%dpdx_face%data, pgf%dpdy_face%data)
      call pgf%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, pgf)
      call destroy_cartesian_metrics(metrics)
   end subroutine run_pgf

   ! -----------------------------------------------------------------
   ! Cases
   ! -----------------------------------------------------------------

   subroutine test_aligned_match(error)
      !! Uniform h_layer (aligned columns) + layered rho.  The
      !! z-correction term has z_R == z_L by construction, so it
      !! must vanish identically.  FV_LITE and Mont must produce
      !! the same dpdx/dpdy face values to round-off.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_pressure_force_t) :: pgf_mont, pgf_fv
      real(wp), parameter :: H0 = 10.0_wp
      real(wp), allocatable :: dpdx_mont(:, :, :), dpdy_mont(:, :, :)
      real(wp) :: max_dpdx_diff, max_dpdy_diff
      integer :: k
      checks: block

         call make_grid(grid, 12, 10)
         ms%nz_ml = NZ
         call ms%init(grid)
         call pgf_mont%init(grid, nz_ml=NZ)
         call pgf_fv%init(grid, nz_ml=NZ)
         pgf_mont%variant = OPGF_VARIANT_MONT
         pgf_fv%variant = OPGF_VARIANT_FV_LITE

         ms%h_layer = H0
         ! Stratified rho: heavier at the bed (k=1), lighter at the
         ! surface (k=NZ).  Spatially uniform within each layer.
         do k = 1, NZ
            ms%rho_layer(:, :, k) = 1035.0_wp - 0.5_wp*real(k - 1, wp)
         end do

         call run_pgf(ms, pgf_mont)
         allocate (dpdx_mont, source=pgf_mont%dpdx_face%data)
         allocate (dpdy_mont, source=pgf_mont%dpdy_face%data)

         call run_pgf(ms, pgf_fv)
         max_dpdx_diff = maxval(abs(pgf_fv%dpdx_face%data - dpdx_mont))
         max_dpdy_diff = maxval(abs(pgf_fv%dpdy_face%data - dpdy_mont))

         call check(error, max_dpdx_diff < 1.0e-12_wp, &
                    "aligned columns: FV dpdx differs from Mont dpdx")
         if (allocated(error)) exit checks
         call check(error, max_dpdy_diff < 1.0e-12_wp, &
                    "aligned columns: FV dpdy differs from Mont dpdy")

      end block checks
      deallocate (dpdx_mont, dpdy_mont)
      call pgf_fv%destroy(); call pgf_mont%destroy(); call ms%destroy()
   end subroutine test_aligned_match

   subroutine test_cancel_spurious(error)
      !! Uniform rho throughout *and* per-column total H uniform,
      !! but h_layer redistributed across layers so layer centres
      !! sit at different z in adjacent columns.  Algebra (see
      !! design notes): the no-geopotential layer-centre pressure
      !! difference reports a non-zero PGF in the redistributed layers
      !! (g/2 * dδ/dx in the simple two-layer split, similar in three);
      !! FV_LITE and `mont` both report zero.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_pressure_force_t) :: pgf_mont, pgf_fv
      real(wp), allocatable :: dpdx_broken(:, :, :)
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp), parameter :: H_THIRD = 10.0_wp
      real(wp), parameter :: DELTA_AMP = 1.0_wp
      real(wp) :: max_broken_pgf, max_mont_pgf, max_fv_pgf, delta_x
      integer :: i, j, nx, ny
      checks: block

         call make_grid(grid, 24, 8)
         ms%nz_ml = NZ
         call ms%init(grid)
         call pgf_mont%init(grid, nz_ml=NZ)
         call pgf_fv%init(grid, nz_ml=NZ)
         pgf_mont%variant = OPGF_VARIANT_MONT
         pgf_fv%variant = OPGF_VARIANT_FV_LITE
         nx = grid%nx_total
         ny = grid%ny_total

         ! Uniform density everywhere → trivial rest state physically.
         ms%rho_layer = 1035.0_wp

         ! Layer thicknesses redistribute by ±delta(x) between k=1 and
         ! k=2, leaving total per-column H = NZ * H_THIRD.
         do j = 1, ny
            do i = 1, nx
               delta_x = DELTA_AMP*sin(2.0_wp*PI*real(i, wp)/real(nx, wp))
               ms%h_layer(i, j, 1) = H_THIRD + delta_x
               ms%h_layer(i, j, 2) = H_THIRD - delta_x
               ms%h_layer(i, j, 3) = H_THIRD
            end do
         end do

         allocate (dpdx_broken(nx, ny, NZ))
         call host_no_geopotential_dpdx(ms%h_layer, ms%rho_layer, 1.0_wp, dpdx_broken)
         max_broken_pgf = maxval(abs(dpdx_broken))

         call run_pgf(ms, pgf_mont)
         max_mont_pgf = maxval(abs(pgf_mont%dpdx_face%data))

         call run_pgf(ms, pgf_fv)
         max_fv_pgf = maxval(abs(pgf_fv%dpdx_face%data))

         ! NON-VACUITY: the misalignment must actually be visible to a scheme
         ! that omits the geopotential term, or the two zeros below are free.
         call check(error, max_broken_pgf > 0.01_wp, &
                    "no-geopotential PGF too small — test setup didn't "// &
                    "generate the spurious signal")
         if (allocated(error)) exit checks
         call check(error, max_fv_pgf < 1.0e-10_wp, &
                    "FV_LITE failed to cancel the spurious uniform-density PGF")
         if (allocated(error)) exit checks
         call check(error, max_mont_pgf < 1.0e-10_wp, &
                    "the Montgomery form failed to cancel the spurious "// &
                    "uniform-density PGF — with uniform rho every rho_star "// &
                    "difference is zero, so M must be identically zero")

      end block checks
      if (allocated(dpdx_broken)) deallocate (dpdx_broken)
      call pgf_fv%destroy(); call pgf_mont%destroy(); call ms%destroy()
   end subroutine test_cancel_spurious

   subroutine test_real_gradient(error)
      !! Add a true horizontal density gradient in one layer.  Both
      !! Mont and FV must respond with a non-zero PGF in that layer.
      !! Sanity-check that FV doesn't accidentally cancel the *real*
      !! signal — it should agree with Mont in sign and order of
      !! magnitude (here aligned columns → agreement to round-off).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_pressure_force_t) :: pgf_mont, pgf_fv
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp), parameter :: H0 = 10.0_wp
      real(wp) :: max_diff, max_mont
      integer :: i, j, nx, ny
      checks: block

         call make_grid(grid, 24, 8)
         ms%nz_ml = NZ
         call ms%init(grid)
         call pgf_mont%init(grid, nz_ml=NZ)
         call pgf_fv%init(grid, nz_ml=NZ)
         pgf_mont%variant = OPGF_VARIANT_MONT
         pgf_fv%variant = OPGF_VARIANT_FV_LITE
         nx = grid%nx_total
         ny = grid%ny_total

         ms%h_layer = H0
         ms%rho_layer = 1035.0_wp
         do j = 1, ny
            do i = 1, nx
               ms%rho_layer(i, j, 1) = 1035.0_wp + &
                                       0.5_wp*sin(2.0_wp*PI*real(i, wp)/real(nx, wp))
            end do
         end do

         call run_pgf(ms, pgf_mont)
         max_mont = maxval(abs(pgf_mont%dpdx_face%data))

         call run_pgf(ms, pgf_fv)
         max_diff = maxval(abs(pgf_fv%dpdx_face%data - pgf_mont%dpdx_face%data))

         call check(error, max_mont > 1.0e-4_wp, &
                    "real density gradient: Mont didn't respond as expected")
         if (allocated(error)) exit checks
         call check(error, max_diff < 1.0e-12_wp, &
                    "aligned columns: FV should match Mont to round-off")

      end block checks
      call pgf_fv%destroy(); call pgf_mont%destroy(); call ms%destroy()
   end subroutine test_real_gradient

end module test_ocean_pgf_fv
