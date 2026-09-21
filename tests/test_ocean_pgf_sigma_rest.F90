!! The defining property of the FV (analytic finite-volume) horizontal
!! pressure-gradient force: with a LINEAR equation of state and T, S
!! LINEAR in z, a motionless column under a flat free surface must produce
!! ZERO acceleration TO ROUND-OFF, for ANY layer geometry — tilted
!! terrain-following layers included.
!!
!! References (cite the paper, never another model's source):
!!   * Adcroft, A., Hallberg, R., & Harrison, M. (2008). A finite volume
!!     discretization of the pressure gradient force using analytic
!!     integration. Ocean Modelling 22(3-4), 106-113.
!!   * White, L., Adcroft, A., & Hallberg, R. (2009). High-order
!!     regridding-remapping schemes for continuous isopycnal and generalized
!!     coordinates in ocean models. J. Comput. Phys. 228(23), 8665-8692.
!!   * Yung, C. K., Hallberg, R. W., Adcroft, A., & Morrison, A. K. (2026).
!!     Assessment of a finite volume discretization of the horizontal
!!     pressure gradient force beneath sloping ice shelves. JAMES 18,
!!     e2025MS005645.  §2.4 / §3.1: the along-face pressure reconstruction,
!!     and Figure 5b->5c — "the pressure gradient acceleration reduces to
!!     numerical precision with the linear stratification for all grid-cell
!!     geometries".
!!
!! WHY THIS TEST EXISTS.  `tests/test_ocean_pgf_reconstruct.F90` asserts the
!! reconstruction improves the per-column moment `intz_dpa`.  It does, and it
!! did so while the ASSEMBLED face acceleration over a sloping bed stayed at
!! 1e-8 m/s^2 — five decades above round-off — because two other pieces were
!! missing:
!!
!!   1. the BOUNDARY layers (k=1 bed, k=nz surface) flattened their T/S
!!      reconstruction to PCM, leaving the full terrain-following truncation
!!      error in exactly the layers adjacent to the tilted boundary; and
!!   2. the along-face integral `intx_dpa` was the two-column TRAPEZOID
!!      `0.5*(dpa_L + dpa_R)`, exact only for a pressure linear in x along the
!!      edge.  Under a tilted interface z is linear in x but p is quadratic in
!!      z, so the trapezoid left a curvature residual
!!      `g*(-drho/dz)*Delta_e^2/12` at EVERY interface — the sigma
!!      "pressure-gradient error of the second kind" (Haney 1991; Mellor,
!!      Ezer & Oey 1994).
!!
!! Those two residuals are invisible to a per-column assertion and to any
!! 2-D x-z check that only looks at the column integral.  This suite asserts
!! the assembled `PFu` / `PFv`, per layer, against a derived round-off bound.
!!
!! THE BOUND.  `pa` is built by subtracting `rho_ref` from an ABSOLUTE EOS
!! density, so the anomaly carries an absolute error of order
!! `eps*rho0`, not `eps*|rho - rho_ref|`.  Propagating that through
!! `pa ~ g*rho0*H`, the face numerator `~ pa*h`, and the divisor
!! `rho0*dx*h` gives
!!
!!      |PF| <~ C * eps * GRAVITY * H / dx
!!
!! with `C` an O(10) constant for the `nz`-deep recurrence.  `C = 100` here
!! leaves ~50x headroom over the measured gfortran value and still sits FIVE
!! decades below the pre-fix acceleration, so the test cannot pass by
!! accident and does not need a per-toolchain tolerance.
module test_ocean_pgf_sigma_rest
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_eos, only: eos_t, EOS_VARIANT_LINEAR
   use rdb_ocean_pressure_force, only: ocean_pressure_force_t, &
                                       ocean_pressure_force_compute, &
                                       OPGF_VARIANT_FV_MOM6
   implicit none
   private

   public :: collect_ocean_pgf_sigma_rest_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NXP = 12, NYP = 8
   integer, parameter :: NZTEST = 15
   real(wp), parameter :: DX = 2000.0_wp

   ! ISOMIP+-flavoured linear EOS: T uniform, all the density signal in S.
   real(wp), parameter :: RHO0 = 1027.0_wp
   real(wp), parameter :: RHO_REF = RHO0
   real(wp), parameter :: T_REF = -1.9_wp, S_REF = 34.2_wp
   real(wp), parameter :: ALPHA_T = 0.0383_wp     ! kg/m^3 per degC
   real(wp), parameter :: BETA_S = 0.80588_wp     ! kg/m^3 per PSU

   ! Linear stratification, z surface-relative and <= 0 (so S increases with
   ! depth).  A weak T gradient too, so the test is not blind to the T path.
   real(wp), parameter :: DSDZ = -1.25e-3_wp      ! PSU per m
   real(wp), parameter :: DTDZ = 4.0e-4_wp        ! degC per m

   ! Bed tilted in BOTH x and y, so PFu and PFv are each exercised.
   real(wp), parameter :: B0 = 226.0_wp
   real(wp), parameter :: BSX = 21.0_wp           ! m of depth per cell in x
   real(wp), parameter :: BSY = 13.0_wp           ! m of depth per cell in y

   ! Round-off bound constant (see the module header).
   real(wp), parameter :: BOUND_C = 100.0_wp
   real(wp), parameter :: EPS_WP = epsilon(1.0_wp)

contains

   subroutine collect_ocean_pgf_sigma_rest_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("sigma_rest_exact_plm", test_sigma_rest_plm), &
                  new_unittest("sigma_rest_exact_ppm", test_sigma_rest_ppm), &
                  new_unittest("stretched_rest_exact_plm", test_stretched_rest_plm), &
                  new_unittest("stretched_rest_exact_ppm", test_stretched_rest_ppm), &
                  new_unittest("pcm_baseline_is_not_exact", test_pcm_not_exact), &
                  new_unittest("boundary_layers_not_worse", test_boundary_layers) &
                  ]
   end subroutine collect_ocean_pgf_sigma_rest_tests

   pure function rho_lin(T, S) result(r)
      !! Linear EOS absolute density (mirror of `eos_density_point` LINEAR).
      real(wp), intent(in) :: T, S
      real(wp) :: r
      r = RHO0 + BETA_S*(S - S_REF) - ALPHA_T*(T - T_REF)
   end function rho_lin

   subroutine make_eos(eos)
      type(eos_t), intent(out) :: eos
      eos%variant = EOS_VARIANT_LINEAR
      eos%rho0 = RHO0
      eos%T_ref = T_REF
      eos%S_ref = S_REF
      eos%alpha_T = ALPHA_T
      eos%beta_S = BETA_S
      eos%is_init = .true.
   end subroutine make_eos

   pure function layer_fracs(stretched) result(frac)
      !! Thickness fractions of the local depth, k=1 bed .. k=nz surface.
      !! They sum to 1, so `sum_k h = b` in every column and the free
      !! surface is flat at z = 0 — there is no barotropic signal to hide
      !! behind.  `stretched = .false.` is uniform sigma; `.true.` is a
      !! surface-refined stack, so the assertion covers a geometry where
      !! `h` varies strongly WITHIN a column as well as across the face.
      logical, intent(in) :: stretched
      real(wp) :: frac(NZTEST)
      real(wp) :: w(NZTEST)
      real(wp) :: tot
      integer :: k
      if (.not. stretched) then
         frac = 1.0_wp/real(NZTEST, wp)
         return
      end if
      do k = 1, NZTEST
         ! Geometric-ish refinement toward the surface (k = nz).
         w(k) = 1.0_wp/(1.0_wp + 0.25_wp*real(NZTEST - k, wp))
      end do
      tot = sum(w)
      do k = 1, NZTEST
         frac(k) = w(k)/tot
      end do
   end function layer_fracs

   subroutine build_rest_state(grid, ms, b, stretched)
      !! Motionless, flat-surfaced, linearly stratified column stack over a
      !! bed tilted in x and y.  The layer mean of a linear profile is its
      !! mid-depth value, so the seeded `hTr` is EXACT — no discretisation
      !! of the initial condition is folded into the measurement.
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(out), allocatable :: b(:, :)
      logical, intent(in) :: stretched
      integer :: i, j, k, nx, ny
      real(wp) :: e_lo, z_mid, Tk, Sk, dz
      real(wp) :: frac(NZTEST)

      nx = size(ms%h_layer, 1)
      ny = size(ms%h_layer, 2)
      frac = layer_fracs(stretched)
      allocate (b(nx, ny))
      do j = 1, ny
         do i = 1, nx
            b(i, j) = B0 + BSX*real(i - NGHOST - 1, wp) + BSY*real(j - NGHOST - 1, wp)
         end do
      end do

      do j = 1, ny
         do i = 1, nx
            e_lo = -b(i, j)
            do k = 1, NZTEST
               dz = frac(k)*b(i, j)
               ms%h_layer(i, j, k) = dz
               z_mid = e_lo + 0.5_wp*dz
               Tk = T_REF + DTDZ*z_mid
               Sk = S_REF + DSDZ*z_mid
               ms%tracers(ms%idx_temperature)%hTr(i, j, k) = Tk*dz
               ms%tracers(ms%idx_salinity)%hTr(i, j, k) = Sk*dz
               ms%rho_layer(i, j, k) = rho_lin(Tk, Sk)
               e_lo = e_lo + dz
            end do
         end do
      end do
      if (grid%nx_total < 1) return   ! keep `grid` referenced; never taken
   end subroutine build_rest_state

   subroutine run_pgf(grid, ms, pgf, eos)
      !! One compute pass; pulls the face accelerations back to the host.
      !! Device rules (`mem:separate`): map BOTH `ms` and `pgf` before the
      !! kernel and pull the COMPONENT arrays, never the aggregate.
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_pressure_force_t), intent(inout) :: pgf
      type(eos_t), intent(in) :: eos
      type(ocean_metrics_t) :: metrics
      call make_cartesian_metrics(metrics, grid)
      !$acc enter data copyin(ms, pgf)
      call ms%enter_data()
      call pgf%enter_data()
      call ocean_pressure_force_compute(grid, metrics, pgf, ms, eos=eos)
      !$acc update self(pgf%dpdx_face%data, pgf%dpdy_face%data)
      call pgf%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, pgf)
      call destroy_cartesian_metrics(metrics)
   end subroutine run_pgf

   subroutine make_pgf(grid, pgf, b, recon, scheme)
      type(hgrid_t), intent(in) :: grid
      type(ocean_pressure_force_t), intent(out) :: pgf
      real(wp), intent(in) :: b(:, :)
      logical, intent(in) :: recon
      integer, intent(in) :: scheme
      call pgf%init(grid, nz_ml=NZTEST)
      pgf%variant = OPGF_VARIANT_FV_MOM6
      pgf%rho0 = RHO0
      pgf%rho_ref = RHO_REF
      pgf%reconstruct_for_pressure = recon
      pgf%recon_scheme = scheme
      call pgf%set_bathymetry(b)
   end subroutine make_pgf

   pure function interior_max(a, nface_x, nface_y) result(m)
      !! Max |a| over INTERIOR faces only (the wall faces are zeroed by
      !! convention and the ghost columns carry an extrapolated bed).
      real(wp), intent(in) :: a(:, :, :)
      integer, intent(in) :: nface_x, nface_y
      real(wp) :: m
      integer :: i, j, k
      m = 0.0_wp
      do k = 1, NZTEST
         do j = NGHOST + 1, NGHOST + nface_y
            do i = NGHOST + 2, NGHOST + nface_x
               m = max(m, abs(a(i, j, k)))
            end do
         end do
      end do
   end function interior_max

   pure function rest_bound(b) result(tol)
      !! `C * eps * g * H / dx` — see the module header.
      real(wp), intent(in) :: b(:, :)
      real(wp) :: tol
      tol = BOUND_C*EPS_WP*GRAVITY*maxval(b)/DX
   end function rest_bound

   subroutine rest_core(error, scheme, stretched)
      !! Headline assertion: the assembled PFu and PFv are at round-off in
      !! EVERY layer, for PLM and PPM, on uniform-sigma and on stretched
      !! layers.  FAILS on the pre-fix kernel at ~1e-8 m/s^2.
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: scheme
      logical, intent(in) :: stretched
      type(multilayer_state_t) :: ms
      type(ocean_pressure_force_t) :: pgf
      type(eos_t) :: eos
      type(hgrid_t) :: grid
      real(wp), allocatable :: b(:, :)
      real(wp) :: tol, pfu_max, pfv_max
      checks: block
         call grid%init(NXP, NYP, NGHOST, DX, DX)
         ms%nz_ml = NZTEST
         call ms%init(grid)
         call make_eos(eos)
         call build_rest_state(grid, ms, b, stretched)
         call make_pgf(grid, pgf, b, .true., scheme)
         call run_pgf(grid, ms, pgf, eos)

         tol = rest_bound(b)
         pfu_max = interior_max(pgf%dpdx_face%data, NXP, NYP)
         pfv_max = interior_max(pgf%dpdy_face%data, NXP, NYP)

         call check(error, pfu_max < tol, &
                    "sigma rest: PFu is not at round-off on tilted layers")
         if (allocated(error)) exit checks
         call check(error, pfv_max < tol, &
                    "sigma rest: PFv is not at round-off on tilted layers")
      end block checks
      if (allocated(b)) deallocate (b)
      call pgf%destroy()
      call ms%destroy()
   end subroutine rest_core

   subroutine test_sigma_rest_plm(error)
      type(error_type), allocatable, intent(out) :: error
      call rest_core(error, 1, .false.)
   end subroutine test_sigma_rest_plm

   subroutine test_sigma_rest_ppm(error)
      type(error_type), allocatable, intent(out) :: error
      call rest_core(error, 2, .false.)
   end subroutine test_sigma_rest_ppm

   subroutine test_stretched_rest_plm(error)
      type(error_type), allocatable, intent(out) :: error
      call rest_core(error, 1, .true.)
   end subroutine test_stretched_rest_plm

   subroutine test_stretched_rest_ppm(error)
      type(error_type), allocatable, intent(out) :: error
      call rest_core(error, 2, .true.)
   end subroutine test_stretched_rest_ppm

   subroutine test_pcm_not_exact(error)
      !! The control that keeps the headline honest: the DEFAULT
      !! (`reconstruct_for_pressure = .false.`) PCM density integral is NOT
      !! exact on this state — it carries the documented sigma truncation
      !! error, decades above the bound.  If this ever drops to round-off
      !! the default path changed and the bit-identity promise is broken.
      type(error_type), allocatable, intent(out) :: error
      type(multilayer_state_t) :: ms
      type(ocean_pressure_force_t) :: pgf
      type(eos_t) :: eos
      type(hgrid_t) :: grid
      real(wp), allocatable :: b(:, :)
      real(wp) :: tol, pfu_max
      checks: block
         call grid%init(NXP, NYP, NGHOST, DX, DX)
         ms%nz_ml = NZTEST
         call ms%init(grid)
         call make_eos(eos)
         call build_rest_state(grid, ms, b, .false.)
         call make_pgf(grid, pgf, b, .false., 1)
         call run_pgf(grid, ms, pgf, eos)
         tol = rest_bound(b)
         pfu_max = interior_max(pgf%dpdx_face%data, NXP, NYP)
         call check(error, pfu_max > 1.0e4_wp*tol, &
                    "pcm baseline: PCM PGF unexpectedly exact — default path changed?")
      end block checks
      if (allocated(b)) deallocate (b)
      call pgf%destroy()
      call ms%destroy()
   end subroutine test_pcm_not_exact

   subroutine test_boundary_layers(error)
      !! Regression guard on the BOUNDARY-layer edges specifically.  Before
      !! the fix, layers k=1 and k=nz flattened to PCM while the interior was
      !! reconstructed: the per-layer error profile went from smooth to a
      !! step, and the layer-to-layer JUMP in the spurious acceleration —
      !! the part that drives a spurious vertical shear, and so the part a
      !! resting sigma column actually feels — grew ~70x over plain PCM even
      !! though the column rms fell.  Assert the two boundary layers are no
      !! worse than the interior, by a wide factor.
      type(error_type), allocatable, intent(out) :: error
      type(multilayer_state_t) :: ms
      type(ocean_pressure_force_t) :: pgf
      type(eos_t) :: eos
      type(hgrid_t) :: grid
      real(wp), allocatable :: b(:, :)
      real(wp) :: bnd_max, int_max, aval
      integer :: i, j, k
      checks: block
         call grid%init(NXP, NYP, NGHOST, DX, DX)
         ms%nz_ml = NZTEST
         call ms%init(grid)
         call make_eos(eos)
         call build_rest_state(grid, ms, b, .false.)
         call make_pgf(grid, pgf, b, .true., 1)
         call run_pgf(grid, ms, pgf, eos)

         bnd_max = 0.0_wp
         int_max = 0.0_wp
         do k = 1, NZTEST
            do j = NGHOST + 1, NGHOST + NYP
               do i = NGHOST + 2, NGHOST + NXP
                  aval = abs(pgf%dpdx_face%data(i, j, k))
                  if (k == 1 .or. k == NZTEST) then
                     bnd_max = max(bnd_max, aval)
                  else
                     int_max = max(int_max, aval)
                  end if
               end do
            end do
         end do

         ! Both are at round-off, so the ratio is ~1; 3x is a loose guard
         ! that still catches the PCM flatten, which put the boundary
         ! layers 8x over the interior on this geometry (and 4 decades over
         ! it at the surface layer of a deeper, more strongly tilted stack).
         call check(error, bnd_max < 3.0_wp*int_max, &
                    "boundary layers: k=1/k=nz PGF error far above the interior")
      end block checks
      if (allocated(b)) deallocate (b)
      call pgf%destroy()
      call ms%destroy()
   end subroutine test_boundary_layers

end module test_ocean_pgf_sigma_rest
