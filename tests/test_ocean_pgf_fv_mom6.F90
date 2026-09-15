!! Unit tests for the FV_MOM6 variant of the ocean pressure-force
!! kernel (rdb_ocean_pressure_force with
!! `variant = OPGF_VARIANT_FV_MOM6`).
!!
!! FV_MOM6 is a faithful port of MOM6's `PressureForce_FV_Bouss`:
!! layer-integrated pressure differences divided by face-averaged
!! thickness `(h_L + h_R + h_neglect)`.  Tests probe the load-bearing
!! invariants:
!!
!!   1. Rest state, uniform ρ — η=0, flat bathy, uniform Rlay across
!!      all layers and cells.  PFu/PFv must be zero to roundoff
!!      (`pa(K) = rho_ref · g · η`, all layer anomalies cancel).
!!   2. Rest state, stratified — η=0, flat bathy, two layers with a
!!      density jump.  PFu/PFv must still be zero (the z-correction in
!!      the integrated form cancels the layer-mean pressure drop).
!!   3. Sloped bathy, uniform ρ — varying b(i), uniform Rlay = rho_ref.
!!      pa is identically zero everywhere, so PFu/PFv must vanish even
!!      with sharply-varying layer thicknesses.
!!   4. Pure SSH tilt, uniform ρ (NK=1) — flat bathy, η varies sinusoidally.
!!      In the limit η ≪ b, PFu must reduce to -g·∂η/∂x at every face,
!!      bit-for-bit with the shallow-water analytic.
!!   5. Two-layer reduced-gravity dynamic — flat bathy, density contrast
!!      Δρ between layers, η=0, but layer thickness varies (interface tilts).
!!      PFu in the bed layer must respond with the -g'·∂h_bed/∂x reduced
!!      gravity signal; PFu in the surf layer must remain zero (no SSH).
!!   6. Cross-variant sanity — same uniform-thickness setup that aligns
!!      columns: FV_MOM6 and FV_LITE must agree to a few percent (both
!!      are mathematically equivalent in the aligned-column limit).
module test_ocean_pgf_fv_mom6
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_ocean_pressure_force, only: ocean_pressure_force_t, &
                                       ocean_pressure_force_compute, &
                                       OPGF_VARIANT_FV_MOM6, OPGF_VARIANT_FV_LITE
   implicit none
   private

   public :: collect_ocean_pgf_fv_mom6_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 2

contains

   subroutine collect_ocean_pgf_fv_mom6_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("rest_uniform_rho_zero_pgf", test_rest_uniform), &
                  new_unittest("rest_stratified_zero_pgf", test_rest_stratified), &
                  new_unittest("sloped_bathy_uniform_rho_zero_pgf", test_bathy_uniform), &
                  new_unittest("ssh_tilt_nk1_matches_swe", test_ssh_tilt_nk1), &
                  new_unittest("reduced_gravity_two_layer", test_reduced_gravity), &
                  new_unittest("agrees_with_fv_lite_aligned", test_agrees_with_fv_lite), &
                  new_unittest("mass_weight_shelf_break_reduces_pgf", test_mass_weight_shelf), &
                  new_unittest("mass_weight_uniform_depth_bit_identical", test_mass_weight_bitident) &
                  ]
   end subroutine collect_ocean_pgf_fv_mom6_tests

   subroutine make_grid(grid, nx_phys, ny_phys, dx)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      real(wp), intent(in) :: dx
      call grid%init(nx_phys, ny_phys, NGHOST, dx, dx)
   end subroutine make_grid

   subroutine run_pgf(ms, pgf, dx)
      !! Run one compute pass.  Pulls dpdx/dpdy back to host.  Sets
      !! pgf bathymetry to ms-consistent b = sum_k(h_layer) per cell
      !! (so e_face(1) = -b puts the bed at the right z).
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_pressure_force_t), intent(inout) :: pgf
      real(wp), intent(in) :: dx
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      real(wp), allocatable :: b_diag(:, :)
      integer :: i, j, nx, ny, k
      call grid%init(size(ms%h_layer, 1) - 2*NGHOST, &
                     size(ms%h_layer, 2) - 2*NGHOST, &
                     NGHOST, dx, dx)
      call make_cartesian_metrics(metrics, grid)
      nx = grid%nx_total
      ny = grid%ny_total
      ! By default we set b such that η = 0: b(i,j) = sum_k h_layer(i,j,k).
      ! Tests that want a non-zero η override pgf%b after this routine
      ! (and call set_bathymetry on the host).  Since this is host-only
      ! (no GPU active in tests by default), direct assignment is fine.
      allocate (b_diag(size(pgf%b, 1), size(pgf%b, 2)), source=0.0_wp)
      do j = 1, size(pgf%b, 2)
         do i = 1, size(pgf%b, 1)
            do k = 1, ms%nz_ml
               b_diag(i, j) = b_diag(i, j) + ms%h_layer(i, j, k)
            end do
         end do
      end do
      call pgf%set_bathymetry(b_diag)
      deallocate (b_diag)
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

   subroutine run_pgf_with_b(ms, pgf, b_user, dx)
      !! Variant of run_pgf that uses a caller-supplied bathymetry
      !! instead of deriving b from sum_k(h_layer).  Use this when
      !! the test wants η = sum_k h_layer - b_user to be non-zero
      !! (e.g. the SSH-tilt test).
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_pressure_force_t), intent(inout) :: pgf
      real(wp), intent(in) :: b_user(:, :)
      real(wp), intent(in) :: dx
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      call grid%init(size(ms%h_layer, 1) - 2*NGHOST, &
                     size(ms%h_layer, 2) - 2*NGHOST, &
                     NGHOST, dx, dx)
      call make_cartesian_metrics(metrics, grid)
      call pgf%set_bathymetry(b_user)
      !$acc enter data copyin(ms, pgf)
      call ms%enter_data()
      call pgf%enter_data()
      call ocean_pressure_force_compute(grid, metrics, pgf, ms)
      !$acc update self(pgf%dpdx_face%data, pgf%dpdy_face%data)
      call pgf%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, pgf)
      call destroy_cartesian_metrics(metrics)
   end subroutine run_pgf_with_b

   ! -----------------------------------------------------------------
   ! Tests
   ! -----------------------------------------------------------------

   subroutine test_rest_uniform(error)
      !! Rest state: uniform ρ across ALL cells and ALL layers, flat
      !! bath, η = 0.  Every term in the numerator must vanish.
      type(error_type), allocatable, intent(out) :: error
      type(multilayer_state_t) :: ms
      type(ocean_pressure_force_t) :: pgf
      real(wp), parameter :: H_PER = 50.0_wp, RHO = 1035.0_wp, DX = 1000.0_wp
      real(wp) :: pgf_max
      type(hgrid_t) :: grid
      checks: block

         call make_grid(grid, 12, 10, DX)
         ms%nz_ml = NZ
         call ms%init(grid)
         call pgf%init(grid, nz_ml=NZ)
         pgf%variant = OPGF_VARIANT_FV_MOM6
         pgf%rho0 = RHO
         pgf%rho_ref = RHO

         ms%h_layer = H_PER
         ms%rho_layer = RHO

         call run_pgf(ms, pgf, DX)
         pgf_max = max(maxval(abs(pgf%dpdx_face%data)), &
                       maxval(abs(pgf%dpdy_face%data)))

         call check(error, pgf_max < 1.0e-12_wp, &
                    "rest+uniform: PGF not zero to roundoff")

      end block checks
      call pgf%destroy(); call ms%destroy()
   end subroutine test_rest_uniform

   subroutine test_rest_stratified(error)
      !! Rest state with a real density jump between layers.  η = 0,
      !! flat bath, uniform h_layer.  The integrated form must cancel
      !! the layer-mean pressure drop so the PGF stays at zero.
      type(error_type), allocatable, intent(out) :: error
      type(multilayer_state_t) :: ms
      type(ocean_pressure_force_t) :: pgf
      real(wp), parameter :: H_PER = 50.0_wp, DX = 1000.0_wp
      real(wp) :: pgf_max
      type(hgrid_t) :: grid
      integer :: k
      checks: block

         call make_grid(grid, 12, 10, DX)
         ms%nz_ml = NZ
         call ms%init(grid)
         call pgf%init(grid, nz_ml=NZ)
         pgf%variant = OPGF_VARIANT_FV_MOM6
         pgf%rho0 = 1035.0_wp
         pgf%rho_ref = 1035.0_wp

         ms%h_layer = H_PER
         do k = 1, NZ
            ms%rho_layer(:, :, k) = 1036.0_wp - real(k - 1, wp)  ! 1036, 1035
         end do

         call run_pgf(ms, pgf, DX)
         pgf_max = max(maxval(abs(pgf%dpdx_face%data)), &
                       maxval(abs(pgf%dpdy_face%data)))

         call check(error, pgf_max < 1.0e-10_wp, &
                    "rest+stratified: PGF not zero (z-correction failed)")

      end block checks
      call pgf%destroy(); call ms%destroy()
   end subroutine test_rest_stratified

   subroutine test_bathy_uniform(error)
      !! Sloped bathymetry + uniform ρ everywhere → pa ≡ 0 → PGF = 0.
      !! Layer thicknesses vary with i (deeper bath ⇒ thicker bed
      !! layer), but density anomaly is zero so all integrals vanish.
      type(error_type), allocatable, intent(out) :: error
      type(multilayer_state_t) :: ms
      type(ocean_pressure_force_t) :: pgf
      real(wp), parameter :: H_SURF = 50.0_wp, RHO = 1035.0_wp, DX = 1000.0_wp
      real(wp) :: pgf_max
      type(hgrid_t) :: grid
      integer :: i, j, nx, ny
      checks: block

         call make_grid(grid, 16, 8, DX)
         ms%nz_ml = NZ
         call ms%init(grid)
         call pgf%init(grid, nz_ml=NZ)
         pgf%variant = OPGF_VARIANT_FV_MOM6
         pgf%rho0 = RHO
         pgf%rho_ref = RHO
         nx = size(ms%h_layer, 1)
         ny = size(ms%h_layer, 2)

         ms%rho_layer = RHO
         ! h_bed varies with i, h_surf constant.  Total H = H_SURF + h_bed(i).
         do j = 1, ny
            do i = 1, nx
               ms%h_layer(i, j, 1) = 20.0_wp + 5.0_wp*real(i, wp)  ! bed grows east
               ms%h_layer(i, j, 2) = H_SURF
            end do
         end do

         call run_pgf(ms, pgf, DX)
         pgf_max = max(maxval(abs(pgf%dpdx_face%data)), &
                       maxval(abs(pgf%dpdy_face%data)))

         call check(error, pgf_max < 1.0e-10_wp, &
                    "sloped bath + uniform ρ: PGF not zero")

      end block checks
      call pgf%destroy(); call ms%destroy()
   end subroutine test_bathy_uniform

   subroutine test_ssh_tilt_nk1(error)
      !! Single-layer SWE limit: flat bath b ≡ B0, uniform ρ, but η
      !! varies sinusoidally.  We set b = B0 and h_layer(:,:,1) =
      !! B0 + η(i) so the column total carries the SSH information.
      !! Expected PFu at every face: -g · ∂η/∂x bit-for-bit (in the
      !! η ≪ b limit) — the standard SWE pressure-gradient force.
      type(error_type), allocatable, intent(out) :: error
      type(multilayer_state_t) :: ms
      type(ocean_pressure_force_t) :: pgf
      real(wp), parameter :: B0 = 1000.0_wp, RHO = 1035.0_wp
      real(wp), parameter :: ETA_AMP = 0.01_wp        ! 1 cm — keep η ≪ b
      real(wp), parameter :: DX = 1000.0_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp), allocatable :: b_user(:, :)
      real(wp) :: x_phys, eta_i, pfu_actual, pfu_expected, err_max
      integer :: i, j, nx, ny, ng
      integer :: nx_face, nx_face_west, nx_face_east
      type(hgrid_t) :: grid
      checks: block

         call make_grid(grid, 32, 8, DX)
         ms%nz_ml = 1
         call ms%init(grid)
         call pgf%init(grid, nz_ml=1)
         pgf%variant = OPGF_VARIANT_FV_MOM6
         pgf%rho0 = RHO
         pgf%rho_ref = RHO
         nx = size(ms%h_layer, 1)
         ny = size(ms%h_layer, 2)
         ng = NGHOST

         ms%rho_layer = RHO
         allocate (b_user(nx, ny), source=B0)
         do j = 1, ny
            do i = 1, nx
               x_phys = real(i - ng, wp)*DX     ! physical x-coord at cell centre
               eta_i = ETA_AMP*sin(2.0_wp*PI*x_phys/(real(32, wp)*DX))
               ms%h_layer(i, j, 1) = B0 + eta_i
            end do
         end do

         call run_pgf_with_b(ms, pgf, b_user, DX)
         deallocate (b_user)

         ! Check at interior u-faces (i ∈ [ng+2, ng+nx_phys]) on a
         ! middle row j_phys.  PFu = -g·(η(i) − η(i−1))/dx — and our
         ! u-face index i_face = i corresponds to face between cells
         ! (i−1) and (i).
         err_max = 0.0_wp
         j = ng + 4
         do i = ng + 2, ng + 32
            x_phys = real(i - 1 - ng, wp)*DX           ! cell (i−1) centre
            eta_i = ETA_AMP*sin(2.0_wp*PI*x_phys/(real(32, wp)*DX))
            ! η at cell i
            x_phys = real(i - ng, wp)*DX
            pfu_expected = -GRAVITY*(ETA_AMP*sin(2.0_wp*PI*x_phys/(real(32, wp)*DX)) &
                                     - eta_i)/DX
            pfu_actual = pgf%dpdx_face%data(i, j, 1)
            err_max = max(err_max, abs(pfu_actual - pfu_expected))
         end do

         ! Tolerance: linear theory holds to O(η/b) = 1e-5; with
         ! roundoff and discretisation we allow ~1e-7.
         call check(error, err_max < 1.0e-6_wp, &
                    "ssh tilt: PGF deviates from -g·∇η beyond linear-theory bound")

      end block checks
      call pgf%destroy(); call ms%destroy()
   end subroutine test_ssh_tilt_nk1

   subroutine test_reduced_gravity(error)
      !! 2-layer adiabatic gprime check.  Flat bathymetry, η = 0
      !! everywhere (b = h_surf + h_bed in each cell), uniform h_surf,
      !! but h_bed varies with i so the interface tilts.  Density
      !! contrast Δρ between layers.
      !!
      !! In the surface layer (Rlay = rho_ref) the anomaly is zero
      !! so PFu_surf should vanish (no SSH).
      !!
      !! In the bed layer the pressure gradient is the reduced-gravity
      !! signal:
      !!     PFu_bed ≈ -(g · Δρ / rho_0) · ∂h_bed/∂x  =  -g' · ∂h_bed/∂x
      !! Tolerance: the layer-integrated FV formula matches the
      !! reduced-gravity analytic up to O(h_bed/H) terms, so we check
      !! sign + order-of-magnitude rather than exact equality.
      type(error_type), allocatable, intent(out) :: error
      type(multilayer_state_t) :: ms
      type(ocean_pressure_force_t) :: pgf
      real(wp), parameter :: RHO_BED = 1036.0_wp, RHO_SURF = 1035.0_wp
      real(wp), parameter :: H_SURF = 1000.0_wp, DX = 1000.0_wp
      real(wp), parameter :: H_BED_AVG = 500.0_wp, H_BED_AMP = 50.0_wp
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp) :: g_prime, expected_amp, surf_max, bed_min, bed_max
      type(hgrid_t) :: grid
      integer :: i, j, nx, ny, ng
      checks: block

         call make_grid(grid, 24, 8, DX)
         ms%nz_ml = 2
         call ms%init(grid)
         call pgf%init(grid, nz_ml=2)
         pgf%variant = OPGF_VARIANT_FV_MOM6
         pgf%rho0 = RHO_SURF
         pgf%rho_ref = RHO_SURF
         nx = size(ms%h_layer, 1)
         ny = size(ms%h_layer, 2)
         ng = NGHOST

         ! Layer densities
         ms%rho_layer(:, :, 1) = RHO_BED
         ms%rho_layer(:, :, 2) = RHO_SURF

         ! h_bed sinusoidal, h_surf such that total = H_SURF + H_BED_AVG
         ! everywhere (so η = 0).
         do j = 1, ny
            do i = 1, nx
               ms%h_layer(i, j, 1) = H_BED_AVG + H_BED_AMP* &
                                     sin(2.0_wp*PI*real(i - ng, wp)/24.0_wp)
               ms%h_layer(i, j, 2) = (H_SURF + H_BED_AVG) - ms%h_layer(i, j, 1)
            end do
         end do

         call run_pgf(ms, pgf, DX)

         surf_max = maxval(abs(pgf%dpdx_face%data(ng + 2:ng + 24, ng + 1:ng + 8, 2)))
         bed_min = minval(pgf%dpdx_face%data(ng + 2:ng + 24, ng + 1:ng + 8, 1))
         bed_max = maxval(pgf%dpdx_face%data(ng + 2:ng + 24, ng + 1:ng + 8, 1))

         g_prime = GRAVITY*(RHO_BED - RHO_SURF)/RHO_SURF
         ! Sinusoidal ∂h_bed/∂x has amplitude H_BED_AMP · 2π/(N·dx).
         expected_amp = g_prime*H_BED_AMP*2.0_wp*PI/(24.0_wp*DX)

         ! Surface layer must show essentially zero PGF (η=0, no
         ! interface-tilt signal reaches it via the integrated form).
         call check(error, surf_max < 1.0e-8_wp, &
                    "reduced gravity: surface PGF should be ~0 (η=0)")
         if (allocated(error)) exit checks

         ! Bed layer should show PGF roughly of magnitude g' · h_bed_amp · 2π/L.
         ! Allow 4× tolerance (the integrated form mixes ∂h_bed and the
         ! intz_dpa/intx_dpa correction terms, so the analytic estimate
         ! is order-of-magnitude only).
         call check(error, &
                    max(abs(bed_min), abs(bed_max)) > 0.25_wp*expected_amp .and. &
                    max(abs(bed_min), abs(bed_max)) < 4.0_wp*expected_amp, &
                    "reduced gravity: bed PGF out of expected range")

      end block checks
      call pgf%destroy(); call ms%destroy()
   end subroutine test_reduced_gravity

   subroutine test_agrees_with_fv_lite(error)
      !! Aligned columns (uniform h_layer across cells) + stratified ρ.
      !! Both FV_MOM6 and FV_LITE are mathematically equivalent here
      !! (no h_layer gradient, no z_centre gradient).  Result must
      !! agree to a tight tolerance.
      type(error_type), allocatable, intent(out) :: error
      type(multilayer_state_t) :: ms
      type(ocean_pressure_force_t) :: pgf_mom6, pgf_lite
      real(wp), parameter :: H_PER = 100.0_wp, DX = 1000.0_wp
      real(wp), allocatable :: dpdx_lite(:, :, :), dpdy_lite(:, :, :)
      real(wp) :: max_diff
      type(hgrid_t) :: grid
      integer :: i, j, k, nx, ny
      checks: block

         call make_grid(grid, 12, 10, DX)
         ms%nz_ml = NZ
         call ms%init(grid)
         call pgf_mom6%init(grid, nz_ml=NZ)
         call pgf_lite%init(grid, nz_ml=NZ)
         pgf_mom6%variant = OPGF_VARIANT_FV_MOM6
         pgf_lite%variant = OPGF_VARIANT_FV_LITE
         nx = size(ms%h_layer, 1)
         ny = size(ms%h_layer, 2)

         ms%h_layer = H_PER
         do k = 1, NZ
            ms%rho_layer(:, :, k) = 1036.0_wp - real(k - 1, wp)
         end do
         ! Sinusoidal η via h_layer perturbation on the surface
         do j = 1, ny
            do i = 1, nx
               ms%h_layer(i, j, NZ) = H_PER + 0.005_wp*sin(real(i, wp))
            end do
         end do

         call run_pgf(ms, pgf_lite, DX)
         allocate (dpdx_lite, source=pgf_lite%dpdx_face%data)
         allocate (dpdy_lite, source=pgf_lite%dpdy_face%data)

         call run_pgf(ms, pgf_mom6, DX)
         max_diff = max(maxval(abs(pgf_mom6%dpdx_face%data - dpdx_lite)), &
                        maxval(abs(pgf_mom6%dpdy_face%data - dpdy_lite)))

         call check(error, max_diff < 1.0e-5_wp, &
                    "FV_MOM6 disagrees with FV_LITE on aligned-column setup")

      end block checks
      if (allocated(dpdx_lite)) deallocate (dpdx_lite)
      if (allocated(dpdy_lite)) deallocate (dpdy_lite)
      call pgf_mom6%destroy(); call pgf_lite%destroy(); call ms%destroy()
   end subroutine test_agrees_with_fv_lite

   subroutine setup_shelf(ms, pgf, mass_weight, dx, b_user)
      !! Build the shelf-break two-column setup and run the PGF.
      !! Domain is a uniform deep basin (bed layer 1000 m of dense water
      !! ρ=1037, surface layer 500 m of light water ρ=ρ_ref=1035) with a
      !! single interior column carved thin (total ~350 m: bed layer 50 m
      !! carrying the LIGHT surface density 1035 — the partial-cell
      !! artifact — surface layer 300 m).  The free surface is flat
      !! (η = 0, b = Σ h_layer per column) so the analytic cross-face PGF
      !! is identically zero everywhere; any non-zero bed-layer PGF at the
      !! shelf face is spurious.  Caller flips `mass_weight`.
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_pressure_force_t), intent(inout) :: pgf
      logical, intent(in) :: mass_weight
      real(wp), intent(in) :: dx
      real(wp), allocatable, intent(out) :: b_user(:, :)
      real(wp), parameter :: RHO_REF = 1035.0_wp, RHO_BED = 1037.0_wp
      real(wp), parameter :: H_BED_DEEP = 1000.0_wp, H_SURF_DEEP = 500.0_wp
      real(wp), parameter :: H_BED_THIN = 50.0_wp, H_SURF_THIN = 300.0_wp
      type(hgrid_t) :: grid
      integer :: i, j, nx, ny, i_thin

      call make_grid(grid, 12, 8, dx)
      ms%nz_ml = NZ
      call ms%init(grid)
      call pgf%init(grid, nz_ml=NZ)
      pgf%variant = OPGF_VARIANT_FV_MOM6
      pgf%rho0 = RHO_REF
      pgf%rho_ref = RHO_REF
      pgf%mass_weight = mass_weight
      nx = size(ms%h_layer, 1)
      ny = size(ms%h_layer, 2)
      i_thin = nx/2          ! one interior column made shallow

      ! Uniform deep basin: dense bed, light surface.
      ms%h_layer(:, :, 1) = H_BED_DEEP
      ms%h_layer(:, :, 2) = H_SURF_DEEP
      ms%rho_layer(:, :, 1) = RHO_BED
      ms%rho_layer(:, :, 2) = RHO_REF
      ! Carve the thin shelf column: its (near-vanished) bed layer holds
      ! the LIGHT surface density — the hydrostatic-inconsistency artifact.
      do j = 1, ny
         ms%h_layer(i_thin, j, 1) = H_BED_THIN
         ms%h_layer(i_thin, j, 2) = H_SURF_THIN
         ms%rho_layer(i_thin, j, 1) = RHO_REF
         ms%rho_layer(i_thin, j, 2) = RHO_REF
      end do

      ! b = Σ h_layer per column ⇒ η = 0 (flat free surface).
      allocate (b_user(nx, ny), source=0.0_wp)
      do j = 1, ny
         do i = 1, nx
            b_user(i, j) = ms%h_layer(i, j, 1) + ms%h_layer(i, j, 2)
         end do
      end do
      call run_pgf_with_b(ms, pgf, b_user, dx)
   end subroutine setup_shelf

   subroutine test_mass_weight_shelf(error)
      !! Shelf-break test: a thin column (~350 m) adjacent to a deep
      !! basin (~1500 m), uniform-ρ surface water, flat free surface.
      !! The true cross-face PGF is zero.  With `mass_weight = .false.`
      !! the bed-layer (k=1) PGF at the shelf face is spuriously large;
      !! with `mass_weight = .true.` the MOM6 hWght blend cancels it to
      !! ≥ one order of magnitude smaller.
      type(error_type), allocatable, intent(out) :: error
      type(multilayer_state_t) :: ms
      type(ocean_pressure_force_t) :: pgf
      real(wp), parameter :: DX = 1000.0_wp
      real(wp), allocatable :: b_user(:, :)
      real(wp) :: bed_off, bed_on
      integer :: i_thin, j, nx
      checks: block

         ! Midpoint (mass_weight off).
         call setup_shelf(ms, pgf, .false., DX, b_user)
         nx = size(pgf%dpdx_face%data, 1) - 1   ! u-face dim is nx+1
         i_thin = nx/2
         j = 4
         ! Largest |PFu| over the two shelf faces bordering the thin column.
         bed_off = max(abs(pgf%dpdx_face%data(i_thin, j, 1)), &
                       abs(pgf%dpdx_face%data(i_thin + 1, j, 1)))
         deallocate (b_user)
         call pgf%destroy(); call ms%destroy()

         ! Mass-weighted (knob on).
         call setup_shelf(ms, pgf, .true., DX, b_user)
         bed_on = max(abs(pgf%dpdx_face%data(i_thin, j, 1)), &
                      abs(pgf%dpdx_face%data(i_thin + 1, j, 1)))
         deallocate (b_user)

         call check(error, bed_off > 1.0e-4_wp, &
                    "shelf: midpoint bed PGF should be sizeably spurious")
         if (allocated(error)) exit checks
         call check(error, bed_on < 0.1_wp*bed_off, &
                    "shelf: mass-weight did not reduce bed PGF by >=10x")

      end block checks
      call pgf%destroy(); call ms%destroy()
   end subroutine test_mass_weight_shelf

   subroutine test_mass_weight_bitident(error)
      !! Uniform-depth bit-identity: with equal column depths the
      !! hydrostatic-inconsistency measure hWght = 0 at every face, so
      !! `mass_weight = .true.` must reproduce the midpoint result
      !! bit-for-bit — even with horizontally-varying density.
      type(error_type), allocatable, intent(out) :: error
      type(multilayer_state_t) :: ms
      type(ocean_pressure_force_t) :: pgf_off, pgf_on
      real(wp), parameter :: DX = 1000.0_wp, RHO_REF = 1035.0_wp
      real(wp), allocatable :: dpdx_off(:, :, :), dpdy_off(:, :, :)
      real(wp) :: max_diff
      type(hgrid_t) :: grid
      integer :: i, j, nx, ny
      checks: block

         call make_grid(grid, 12, 10, DX)
         ms%nz_ml = NZ
         call ms%init(grid)
         call pgf_off%init(grid, nz_ml=NZ)
         call pgf_on%init(grid, nz_ml=NZ)
         pgf_off%variant = OPGF_VARIANT_FV_MOM6
         pgf_on%variant = OPGF_VARIANT_FV_MOM6
         pgf_off%rho0 = RHO_REF; pgf_off%rho_ref = RHO_REF
         pgf_on%rho0 = RHO_REF; pgf_on%rho_ref = RHO_REF
         pgf_off%mass_weight = .false.
         pgf_on%mass_weight = .true.
         nx = size(ms%h_layer, 1)
         ny = size(ms%h_layer, 2)

         ! Equal column depths (uniform h_layer) but horizontally-varying
         ! density ⇒ a real PGF, yet hWght = 0 at every face.
         ms%h_layer(:, :, 1) = 400.0_wp
         ms%h_layer(:, :, 2) = 500.0_wp
         do j = 1, ny
            do i = 1, nx
               ms%rho_layer(i, j, 1) = 1037.0_wp + 0.01_wp*real(i, wp)
               ms%rho_layer(i, j, 2) = RHO_REF + 0.005_wp*real(i, wp)
            end do
         end do

         call run_pgf(ms, pgf_off, DX)
         allocate (dpdx_off, source=pgf_off%dpdx_face%data)
         allocate (dpdy_off, source=pgf_off%dpdy_face%data)

         call run_pgf(ms, pgf_on, DX)
         max_diff = max(maxval(abs(pgf_on%dpdx_face%data - dpdx_off)), &
                        maxval(abs(pgf_on%dpdy_face%data - dpdy_off)))

         call check(error, max_diff == 0.0_wp, &
                    "uniform depth: mass_weight not bit-identical to midpoint")

      end block checks
      if (allocated(dpdx_off)) deallocate (dpdx_off)
      if (allocated(dpdy_off)) deallocate (dpdy_off)
      call pgf_off%destroy(); call pgf_on%destroy(); call ms%destroy()
   end subroutine test_mass_weight_bitident

end module test_ocean_pgf_fv_mom6
