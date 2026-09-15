!! Analytical tests for the Montgomery-potential PGF
!! (`&ocean_pgf_nml form = "mont"` → `OPGF_VARIANT_MONT`).
!!
!! The defect these pin
!! --------------------
!! `form="mont"` used to compute its faces as
!!
!!     PGF_x = -(1/rho0) * (p_centre^R - p_centre^L) / dx
!!
!! and nothing else.  That is a horizontal difference of layer-centre
!! pressure ALONG the coordinate surface with NO geopotential term — the one
!! term that makes differencing a layer quantity legitimate in the first
!! place.  It is exact only where the two layer centres either side of a face
!! sit at the same z; the moment the bed slopes, or the layer thickness
!! varies across the face, what is left is a pressure gradient on a
!! MOTIONLESS ocean.  Every sloping shipped case NaN'd within a day.
!!
!! A genuine Boussinesq Montgomery potential
!!
!!     M = p/rho0 + rho_star*z,        rho_star = g*rho_layer/rho0
!!
!! is CONSTANT through a layer of uniform density (rising by `dz` costs
!! `-g*rho*dz/rho0` of `p/rho0` and gains exactly `rho_star*dz`), so ONE
!! horizontal difference of `M` IS the acceleration — `M` already carries
!! `g*z`.  Built bottom-up (k=1 bed, k=nz surface) from the free surface,
!! where the surface-relative interface height and the pressure both vanish:
!!
!!     e_edge(k) = z_centre(k) + 0.5*h_layer(k)          ! TOP of layer k
!!     M(nz)     = 0
!!     M(k)      = M(k+1) + (rho_star(k) - rho_star(k+1)) * e_edge(k+1)
!!
!!     PGF_x = -(M_R - M_L)*idxCu + (rho_star_R - rho_star_L)*z_eff*idxCu
!!     z_eff = (e_L*h_R + e_R*h_L - h_L*h_R) / (h_L + h_R)
!!
!! The second term is not a refinement.  `-dM/dx` alone is the PGF only where
!! density is horizontally uniform WITHIN the layer; where it is not, the
!! exact relation picks up `+ z * d(rho_star)/dx`, and `z_eff` is the
!! thickness-weighted height at which to evaluate it.  Without it the scheme
!! gets the sign of a horizontal density contrast wrong.
!!
!! Cases
!! -----
!!   * `mont_is_at_rest_over_a_sloping_bed` — the regression proper.  Flat
!!     interior isopycnals with the bathymetry absorbed entirely by the bed
!!     layer: an EXACT rest state (pressure horizontally uniform at every
!!     depth).  Montgomery must report zero.  NON-VACUITY: the OLD
!!     no-geopotential arithmetic is evaluated on the SAME state and asserted
!!     to produce a large spurious gradient, so the case cannot pass by being
!!     trivially quiet.
!!   * `mont_matches_fv_lite_on_aligned_columns` — on columns of equal layer
!!     thickness `z_eff` collapses to the mean layer centre and the whole
!!     expression reduces ALGEBRAICALLY to `fv_lite`.  Asserted against a live
!!     `fv_lite` compute over a horizontal density wave, to round-off.
!!   * `mont_density_term_is_load_bearing` — the same density wave, with
!!     `-dM/dx` evaluated ALONE.  Asserted to be badly wrong (and, in the
!!     deep layers, of the WRONG SIGN), which is what makes the previous
!!     case a real test of the `z_eff` term rather than of the recursion.
module test_ocean_pgf_mont
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

   public :: collect_ocean_pgf_mont_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 8
   integer, parameter :: NX = 32
   integer, parameter :: NY = 6
   real(wp), parameter :: DX = 8000.0_wp
   real(wp), parameter :: RHO0 = 1035.0_wp
      !! Same default the slot carries; pinned here so the host reference
      !! arithmetic below uses the identical divisor.
   real(wp), parameter :: RHO_LIGHTEST = 1035.0_wp
   real(wp), parameter :: RHO_RANGE = 2.0_wp
   real(wp), parameter :: H_UPPER = 100.0_wp
      !! Thickness of every layer above the bed layer in the sloping-bed
      !! case.  Uniform across columns, which is what keeps the interior
      !! interfaces FLAT while the bed slopes.

contains

   subroutine collect_ocean_pgf_mont_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("mont_is_at_rest_over_a_sloping_bed", test_mont_rest_over_slope), &
                  new_unittest("mont_matches_fv_lite_on_aligned_columns", test_mont_matches_fv_lite), &
                  new_unittest("mont_density_term_is_load_bearing", test_mont_density_term) &
                  ]
   end subroutine collect_ocean_pgf_mont_tests

   ! ---------------------------------------------------------------------
   ! Scaffolding
   ! ---------------------------------------------------------------------

   subroutine run_pgf(ms, pgf)
      !! One compute pass, with the face buffers pulled back to the host.
      !! `dpdx/dpdy` are scratch (create/delete), so the `!$acc update self`
      !! before `exit_data` is required for the host-side comparison — the
      !! GPU build maps `mem:separate` and would otherwise compare stale host
      !! memory.  Every directive is an inert no-op on host builds.
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_pressure_force_t), intent(inout) :: pgf
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      call grid%init(NX, NY, NGHOST, DX, DX)
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

   pure subroutine host_interface_heights(h_layer, e_edge)
      !! Surface-relative interface heights, `e_edge(k)` = TOP of layer k
      !! (so `e_edge(nz+1)` would be the bed and is never needed).  Built by
      !! summing DOWN from the free surface, the same reference the kernel
      !! uses — bottom-up indexing, k=1 bed, k=NZ surface.
      real(wp), intent(in) :: h_layer(:, :, :)
      real(wp), intent(out) :: e_edge(:, :, :)
      integer :: i, j, k
      do j = 1, size(h_layer, 2)
         do i = 1, size(h_layer, 1)
            e_edge(i, j, NZ) = 0.0_wp
            do k = NZ - 1, 1, -1
               e_edge(i, j, k) = e_edge(i, j, k + 1) - h_layer(i, j, k + 1)
            end do
         end do
      end do
   end subroutine host_interface_heights

   pure subroutine host_old_mont_dpdx(h_layer, rho_layer, dpdx)
      !! The ARITHMETIC THAT WAS THERE: layer-centre pressure differenced
      !! horizontally, no geopotential term.  Reproduced here (not imported)
      !! precisely so the regression can assert it is wrong.
      real(wp), intent(in) :: h_layer(:, :, :), rho_layer(:, :, :)
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
               dpdx(i, j, k) = -(1.0_wp/RHO0)*(p_c_r - p_c_l)/DX
            end do
         end do
      end do
   end subroutine host_old_mont_dpdx

   pure subroutine host_mont_dM_only_dpdx(h_layer, rho_layer, dpdx)
      !! The Montgomery recursion with the horizontal-density term DELETED —
      !! i.e. `-dM/dx` on its own, which is what a naive reading of the
      !! recursion gives.  Used to show the `z_eff` term is load-bearing.
      real(wp), intent(in) :: h_layer(:, :, :), rho_layer(:, :, :)
      real(wp), intent(out) :: dpdx(:, :, :)
      real(wp) :: e_edge(size(h_layer, 1), size(h_layer, 2), NZ)
      real(wp) :: mm(size(h_layer, 1), size(h_layer, 2), NZ)
      real(wp) :: g_over_rho0
      integer :: i, j, k, nx, ny
      nx = size(h_layer, 1)
      ny = size(h_layer, 2)
      g_over_rho0 = GRAVITY/RHO0
      call host_interface_heights(h_layer, e_edge)
      do j = 1, ny
         do i = 1, nx
            mm(i, j, NZ) = 0.0_wp
            do k = NZ - 1, 1, -1
               mm(i, j, k) = mm(i, j, k + 1) + &
                             g_over_rho0*(rho_layer(i, j, k) - rho_layer(i, j, k + 1))* &
                             e_edge(i, j, k)
            end do
         end do
      end do
      dpdx = 0.0_wp
      do k = 1, NZ
         do j = 1, ny
            do i = 2, nx
               dpdx(i, j, k) = -(mm(i, j, k) - mm(i - 1, j, k))/DX
            end do
         end do
      end do
   end subroutine host_mont_dM_only_dpdx

   subroutine build_sloping_bed_rest_state(ms, grid)
      !! An EXACT rest state over a sloping bed: every layer above the bed is
      !! `H_UPPER` thick in every column, so all interior interfaces are FLAT,
      !! and the bed layer alone absorbs the bathymetry (1300 m in the deep
      !! column, 100 m in the shallow one — a 13x contrast).  Per-layer
      !! density is horizontally uniform, so pressure is horizontally uniform
      !! at every depth and nothing may move.
      type(multilayer_state_t), intent(inout) :: ms
      type(hgrid_t), intent(in) :: grid
      integer :: i, j, k, nx, ny
      real(wp) :: frac, h_bed
      nx = grid%nx_total
      ny = grid%ny_total
      do j = 1, ny
         do i = 1, nx
            frac = real(i - 1, wp)/real(nx - 1, wp)
            h_bed = 1300.0_wp + frac*(100.0_wp - 1300.0_wp)
            ms%h_layer(i, j, 1) = h_bed
            do k = 2, NZ
               ms%h_layer(i, j, k) = H_UPPER
            end do
         end do
      end do
      do k = 1, NZ
         ms%rho_layer(:, :, k) = RHO_LIGHTEST + RHO_RANGE* &
                                 (real(NZ - k, wp) + 0.5_wp)/real(NZ, wp)
      end do
   end subroutine build_sloping_bed_rest_state

   subroutine build_aligned_density_wave(ms, grid)
      !! Columns of IDENTICAL layer thickness (250 m x 8 = 2000 m) carrying a
      !! horizontal density wave.  Aligned columns are exactly the regime
      !! where the Montgomery face expression must reduce to `fv_lite`, and
      !! the density wave is what switches the `z_eff` term on.
      !!
      !! The wave amplitude is DEPTH-WEIGHTED — full strength in the bed
      !! layer, zero at the surface, like a dense bottom-water tongue.  A
      !! depth-UNIFORM wave would be degenerate for the third case: it leaves
      !! every `rho_star` difference between adjacent layers untouched, so
      !! `dM/dx` vanishes identically and there is nothing left to compare the
      !! `z_eff` term against.
      type(multilayer_state_t), intent(inout) :: ms
      type(hgrid_t), intent(in) :: grid
      integer :: i, j, k, nx, ny
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp) :: weight
      nx = grid%nx_total
      ny = grid%ny_total
      ms%h_layer = 250.0_wp
      do k = 1, NZ
         weight = real(NZ - k, wp)/real(NZ - 1, wp)
         do j = 1, ny
            do i = 1, nx
               ms%rho_layer(i, j, k) = RHO_LIGHTEST + RHO_RANGE* &
                                       (real(NZ - k, wp) + 0.5_wp)/real(NZ, wp) &
                                       + 0.35_wp*weight* &
                                       sin(2.0_wp*PI*real(i, wp)/real(nx, wp))
            end do
         end do
      end do
   end subroutine build_aligned_density_wave

   ! ---------------------------------------------------------------------
   ! Cases
   ! ---------------------------------------------------------------------

   subroutine test_mont_rest_over_slope(error)
      !! THE regression.  Flat isopycnals over a bed that shoals by 1200 m:
      !! an exact rest state, so the Montgomery PGF must be zero — and the
      !! no-geopotential arithmetic it replaced must NOT be, or the assertion
      !! above proves nothing.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_pressure_force_t) :: pgf
      real(wp), allocatable :: dpdx_old(:, :, :)
      real(wp) :: max_mont, max_old
      checks: block

         call grid%init(NX, NY, NGHOST, DX, DX)
         ms%nz_ml = NZ
         call ms%init(grid)
         call pgf%init(grid, nz_ml=NZ)
         pgf%variant = OPGF_VARIANT_MONT
         pgf%rho0 = RHO0

         call build_sloping_bed_rest_state(ms, grid)

         allocate (dpdx_old(grid%nx_total, grid%ny_total, NZ))
         call host_old_mont_dpdx(ms%h_layer, ms%rho_layer, dpdx_old)
         max_old = maxval(abs(dpdx_old))

         call run_pgf(ms, pgf)
         max_mont = maxval(abs(pgf%dpdx_face%data))

         ! NON-VACUITY.  The arithmetic that was in the tree must visibly
         ! break on this state, otherwise the zero below is free.
         call check(error, max_old > 1.0e-3_wp, &
                    "the no-geopotential form produced no spurious gradient "// &
                    "over the sloping bed — this regression would pass vacuously")
         if (allocated(error)) exit checks

         ! A resting ocean.  Not "small": zero.
         call check(error, max_mont < 1.0e-14_wp, &
                    "flat isopycnals over a sloping bed are a rest state, but "// &
                    "the Montgomery PGF is non-zero — the geopotential term "// &
                    "is missing or mis-signed")
         if (allocated(error)) exit checks

         ! And the v-faces, which are uniform in j and must be zero too.
         call check(error, maxval(abs(pgf%dpdy_face%data)) < 1.0e-14_wp, &
                    "the j-uniform state produced a non-zero north-face "// &
                    "Montgomery PGF")

      end block checks
      if (allocated(dpdx_old)) deallocate (dpdx_old)
      call pgf%destroy()
      call ms%destroy()
   end subroutine test_mont_rest_over_slope

   subroutine test_mont_matches_fv_lite(error)
      !! On aligned columns `z_eff` is exactly the mean layer centre and the
      !! Montgomery face expression is ALGEBRAICALLY `fv_lite`.  Held to
      !! round-off against a live `fv_lite` compute, not to a loose tolerance
      !! — the two differ only in floating-point association.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_pressure_force_t) :: pgf_mont, pgf_fv
      real(wp), allocatable :: dpdx_mont(:, :, :)
      real(wp) :: scale_x, max_diff
      checks: block

         call grid%init(NX, NY, NGHOST, DX, DX)
         ms%nz_ml = NZ
         call ms%init(grid)
         call pgf_mont%init(grid, nz_ml=NZ)
         call pgf_fv%init(grid, nz_ml=NZ)
         pgf_mont%variant = OPGF_VARIANT_MONT
         pgf_fv%variant = OPGF_VARIANT_FV_LITE
         pgf_mont%rho0 = RHO0
         pgf_fv%rho0 = RHO0

         call build_aligned_density_wave(ms, grid)

         call run_pgf(ms, pgf_mont)
         dpdx_mont = pgf_mont%dpdx_face%data

         call run_pgf(ms, pgf_fv)

         scale_x = maxval(abs(pgf_fv%dpdx_face%data))
         ! Guard the setup: a zero reference would make the comparison free.
         call check(error, scale_x > 1.0e-6_wp, &
                    "the density wave produced no fv_lite PGF to compare "// &
                    "against — the comparison would pass vacuously")
         if (allocated(error)) exit checks

         max_diff = maxval(abs(dpdx_mont - pgf_fv%dpdx_face%data))
         ! 1e-9 relative, not bitwise: the two forms are the same expression
         ! rearranged, but `fv_lite` reaches it by differencing a ~2e7 Pa
         ! pressure stack to recover ~1e3 Pa, a condition number of ~1e4 that
         ! costs it four digits.  Montgomery differences DENSITIES and never
         ! forms the large number.  Measured agreement is ~8e-12 relative;
         ! the bound is two orders looser than that and eight orders tighter
         ! than any discretisation difference would be.
         call check(error, max_diff < 1.0e-9_wp*scale_x, &
                    "Montgomery and fv_lite disagree on ALIGNED columns, "// &
                    "where they are the same expression rearranged")

      end block checks
      if (allocated(dpdx_mont)) deallocate (dpdx_mont)
      call pgf_fv%destroy()
      call pgf_mont%destroy()
      call ms%destroy()
   end subroutine test_mont_matches_fv_lite

   subroutine test_mont_density_term(error)
      !! NON-VACUITY for the `z_eff` term.  `-dM/dx` alone — the recursion
      !! without the horizontal-density correction — is not a small error on
      !! the same density wave: in the deep layers it carries the OPPOSITE
      !! SIGN, because the `z * d(rho_star)/dx` it drops grows with depth
      !! while the term it keeps does not.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_pressure_force_t) :: pgf
      real(wp), allocatable :: dpdx_nobc(:, :, :)
      real(wp) :: scale_x, max_diff, a_full, a_nobc
      integer :: i, j
      logical :: sign_flipped
      checks: block

         call grid%init(NX, NY, NGHOST, DX, DX)
         ms%nz_ml = NZ
         call ms%init(grid)
         call pgf%init(grid, nz_ml=NZ)
         pgf%variant = OPGF_VARIANT_MONT
         pgf%rho0 = RHO0

         call build_aligned_density_wave(ms, grid)

         allocate (dpdx_nobc(grid%nx_total, grid%ny_total, NZ))
         call host_mont_dM_only_dpdx(ms%h_layer, ms%rho_layer, dpdx_nobc)

         call run_pgf(ms, pgf)
         scale_x = maxval(abs(pgf%dpdx_face%data))
         max_diff = 0.0_wp
         sign_flipped = .false.
         j = NGHOST + 1
         do i = 2, grid%nx_total
            a_full = pgf%dpdx_face%data(i, j, 1)
            a_nobc = dpdx_nobc(i, j, 1)
            max_diff = max(max_diff, abs(a_full - a_nobc))
            if (abs(a_full) > 0.1_wp*scale_x .and. a_full*a_nobc < 0.0_wp) then
               sign_flipped = .true.
            end if
         end do

         call check(error, max_diff > scale_x, &
                    "dropping the horizontal-density term changed the bed-layer "// &
                    "PGF by less than its own magnitude — the term under test "// &
                    "is not being exercised")
         if (allocated(error)) exit checks

         call check(error, sign_flipped, &
                    "dropping the horizontal-density term did not flip the sign "// &
                    "of the bed-layer PGF anywhere — the case no longer "// &
                    "demonstrates why the term is load-bearing")

      end block checks
      if (allocated(dpdx_nobc)) deallocate (dpdx_nobc)
      call pgf%destroy()
      call ms%destroy()
   end subroutine test_mont_density_term

end module test_ocean_pgf_mont
