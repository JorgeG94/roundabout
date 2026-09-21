!! The FV in-layer T/S reconstruction must be FILLER-AWARE: a vanished
!! layer (`h <= H_VANISHED`) carries no water and must contribute nothing
!! to the pressure-gradient force.
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
!!     e2025MS005645.  §2.4 steps 3-4 (the vertical + along-face Boole
!!     quadratures), §3.1 (the surface-pressure reconstruction beneath a
!!     sloping ice base), §3.3 (vanishing layers).
!!
!! WHY THIS SUITE EXISTS.  Under a z-like coordinate (`VCOORD_Z_FIXED`)
!! every ice-covered or shallow column carries a run of inert FILLER layers
!! at `zstar_h_min`, and the first LIVE layer beside that run is a PARTIAL
!! cell cut by the ice base or the bed.  `test_ocean_pgf_sigma_rest`
!! asserts the reconstruction is exact on tilted sigma layers — where no
!! layer is ever vanished — and says nothing about this geometry.  Before
!! the fix the filler's layer mean `hTr/H_VANISHED` (a cell the ALE drain
!! emptied, so `T = S = 0`) entered the centred PLM/PPM slope of the
!! adjacent partial cell, the pa stack every live layer below inherits, and
!! the along-face quadrature at the wedge layers.  Measured on the two
!! geometries below, `reconstruct_for_pressure = .true.` was then ~1000x
!! WORSE than the PCM default it is supposed to beat.
!!
!! THE THREE PROPERTIES ASSERTED
!!
!!   1. On a bed-side staircase under a FLAT free surface the reconstructed
!!      PGF is at ROUND-OFF at the partial-cell face, while the PCM default
!!      sits at the closed-form quadrature error
!!      `numer_err = (g*gamma/12)*[h_L^3 - h_R^3 + de_top^3 - de_bot^3]`
!!      (the `intz_dpa` mid-point first moment plus the `intx_pa*de` edge
!!      trapezoid).  The PCM value is asserted AGAINST that formula, so the
!!      test states the mechanism, not just a number.
!!   2. On an ice-base staircase (`p_top_in_bc`) every ALIGNED face is at
!!      round-off, and the PARTIAL-cell face is left with exactly the
!!      Yung et al. (2026) §3.1 residual — the `intx_pa(nz+1)` seed
!!      `0.5*(pa_L + pa_R)` is a TRAPEZOID over an ice base along which the
!!      pressure is quadratic, leaving `g*gamma*de_surf^2/12` that the
!!      `(h_R - h_L)` weight then picks up.  That term is NOT this
!!      reconstruction's to remove (it is the surface-pressure
!!      reconstruction, deliberately untouched here) and it is asserted in
!!      closed form so it cannot silently grow.
!!   3. A column with NO vanished layer is BIT-IDENTICAL to the
!!      unsegmented builder.  The reference is a verbatim copy of the
!!      pre-filler-awareness algorithm kept in this file, so the assertion
!!      is `==`, not a tolerance.
module test_ocean_pgf_filler_recon
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_eos, only: eos_t, EOS_VARIANT_LINEAR
   use rdb_ocean_pgf_reconstruct, only: plm_edges_column, ppm_edges_column
   use rdb_ocean_pressure_force, only: ocean_pressure_force_t, &
                                       ocean_pressure_force_compute, &
                                       OPGF_VARIANT_FV_MOM6
   implicit none
   private

   public :: collect_ocean_pgf_filler_recon_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NXP = 6, NYP = 3
   integer, parameter :: NZTEST = 8
   real(wp), parameter :: DX = 2000.0_wp

   real(wp), parameter :: HNOM = 20.0_wp
      !! Nominal z-level thickness of the synthetic `z_fixed` stack.
   real(wp), parameter :: HFILL = 1.0e-4_wp
      !! `zstar_h_min` — the inert filler thickness (< `H_VANISHED`).

   ! ISOMIP+-flavoured linear EOS: T uniform, all the density signal in S.
   real(wp), parameter :: RHO0 = 1027.0_wp
   real(wp), parameter :: RHO_REF = RHO0
   real(wp), parameter :: T_REF = -1.9_wp, S_REF = 34.2_wp
   real(wp), parameter :: ALPHA_T = 0.0383_wp     ! kg/m^3 per degC
   real(wp), parameter :: BETA_S = 0.80588_wp     ! kg/m^3 per PSU
   real(wp), parameter :: DSDZ = -1.0583e-3_wp    ! PSU per m (z up, <= 0)
   real(wp), parameter :: GAMMA = -BETA_S*DSDZ    ! -drho/dz, kg/m^3 per m

   ! Round-off bound `C*eps*g*H/dx` — the same derivation as
   ! `test_ocean_pgf_sigma_rest` (the anomaly carries an absolute error of
   ! order `eps*rho0`, not `eps*|rho - rho_ref|`).
   real(wp), parameter :: BOUND_C = 100.0_wp
   real(wp), parameter :: EPS_WP = epsilon(1.0_wp)
   real(wp), parameter :: REST_TOL = BOUND_C*EPS_WP*GRAVITY*160.0_wp/DX

contains

   subroutine collect_ocean_pgf_filler_recon_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("bed_staircase_pcm_matches_quadrature_error", &
                               test_bed_pcm_analytic), &
                  new_unittest("bed_staircase_reconstruct_exact_plm", &
                               test_bed_recon_plm), &
                  new_unittest("bed_staircase_reconstruct_exact_ppm", &
                               test_bed_recon_ppm), &
                  new_unittest("ice_staircase_aligned_faces_exact_plm", &
                               test_ice_recon_plm), &
                  new_unittest("ice_staircase_aligned_faces_exact_ppm", &
                               test_ice_recon_ppm), &
                  new_unittest("ice_staircase_partial_face_is_the_surface_seed", &
                               test_ice_partial_face_seed), &
                  new_unittest("flat_lid_over_fillers_is_bit_zero", &
                               test_flat_lid_bit_zero), &
                  new_unittest("no_filler_column_bit_identical_plm", &
                               test_bit_identical_plm), &
                  new_unittest("no_filler_column_bit_identical_ppm", &
                               test_bit_identical_ppm) &
                  ]
   end subroutine collect_ocean_pgf_filler_recon_tests

   ! ---------------------------------------------------------------------
   ! State builders
   ! ---------------------------------------------------------------------

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

   subroutine seed_column(ms, i, j, e_bed, h)
      !! One column from a thickness stack: T uniform, S linear in z, the
      !! layer mean of a linear profile being its mid-depth value (so the
      !! seeded `hTr` is EXACT and no IC discretisation is folded in).
      !! A FILLER is seeded the way the ALE drain leaves one — `hTr = 0`,
      !! and `rho_layer = rho_0` the way `eos_linear_impl`'s vanished-layer
      !! branch does — which is exactly the state that poisoned the
      !! reconstruction.
      type(multilayer_state_t), intent(inout) :: ms
      integer, intent(in) :: i, j
      real(wp), intent(in) :: e_bed, h(NZTEST)
      integer :: k
      real(wp) :: e_lo, zmid, Sk
      e_lo = e_bed
      do k = 1, NZTEST
         ms%h_layer(i, j, k) = h(k)
         zmid = e_lo + 0.5_wp*h(k)
         Sk = S_REF + DSDZ*zmid
         if (h(k) > 1.5e-4_wp) then
            ms%tracers(ms%idx_temperature)%hTr(i, j, k) = T_REF*h(k)
            ms%tracers(ms%idx_salinity)%hTr(i, j, k) = Sk*h(k)
            ms%rho_layer(i, j, k) = rho_lin(T_REF, Sk)
         else
            ms%tracers(ms%idx_temperature)%hTr(i, j, k) = 0.0_wp
            ms%tracers(ms%idx_salinity)%hTr(i, j, k) = 0.0_wp
            ms%rho_layer(i, j, k) = RHO0
         end if
         e_lo = e_lo + h(k)
      end do
   end subroutine seed_column

   subroutine run_pgf(grid, ms, b, p_top_in, recon, scheme, pf)
      !! One compute pass, returning `dpdx_face` on the host.  Device rules
      !! (`mem:separate`): map BOTH `ms` and `pgf` before the kernel and
      !! pull the COMPONENT array, never the aggregate.
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: b(:, :)
      logical, intent(in) :: p_top_in, recon
      integer, intent(in) :: scheme
      real(wp), allocatable, intent(out) :: pf(:, :, :)
      type(ocean_pressure_force_t) :: pgf
      type(ocean_metrics_t) :: metrics
      type(eos_t) :: eos
      call make_eos(eos)
      call pgf%init(grid, nz_ml=NZTEST)
      pgf%variant = OPGF_VARIANT_FV_MOM6
      pgf%rho0 = RHO0
      pgf%rho_ref = RHO_REF
      pgf%reconstruct_for_pressure = recon
      pgf%recon_scheme = scheme
      pgf%p_top_in_bc = p_top_in
      call pgf%set_bathymetry(b)
      call make_cartesian_metrics(metrics, grid)
      !$acc enter data copyin(ms, pgf)
      call ms%enter_data()
      call pgf%enter_data()
      call ocean_pressure_force_compute(grid, metrics, pgf, ms, eos=eos)
      !$acc update self(pgf%dpdx_face%data)
      allocate (pf, source=pgf%dpdx_face%data)
      call pgf%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, pgf)
      call destroy_cartesian_metrics(metrics)
      call pgf%destroy()
   end subroutine run_pgf

   subroutine build_bed_staircase(grid, ms, b, iface)
      !! FLAT free surface at z = 0, no ice.  West of `iface` the bed is
      !! -160 m (8 x 20 m, nothing vanished); east of it the bed is -130 m,
      !! which under `z_fixed` vanishes the k=1 level to a filler and cuts
      !! k=2 into a 10 m partial cell.  The staircase face is `iface`; the
      !! partial-cell face is layer k=2 and the faces above it are aligned
      !! 20 m <-> 20 m.
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), allocatable, intent(out) :: b(:, :)
      integer, intent(out) :: iface
      real(wp) :: hdeep(NZTEST), hshal(NZTEST)
      integer :: i, j, k, nxt, nyt
      nxt = size(ms%h_layer, 1)
      nyt = size(ms%h_layer, 2)
      allocate (b(nxt, nyt))
      hdeep = HNOM
      hshal(1) = HFILL
      hshal(2) = 10.0_wp - HFILL
      do k = 3, NZTEST
         hshal(k) = HNOM
      end do
      iface = NGHOST + 3
      do j = 1, nyt
         do i = 1, nxt
            if (i < iface) then
               b(i, j) = 160.0_wp
               call seed_column(ms, i, j, -160.0_wp, hdeep)
            else
               b(i, j) = 130.0_wp
               call seed_column(ms, i, j, -130.0_wp, hshal)
            end if
         end do
      end do
      ms%p_top = 0.0_wp
      if (grid%nx_total < 1) return   ! keep `grid` referenced; never taken
   end subroutine build_bed_staircase

   subroutine build_ice_staircase(grid, ms, b, iface, de_surf)
      !! ICE-BASE staircase over a flat bed at -160 m.  West of `iface` the
      !! column is open ocean (eta = 0, 8 x 20 m); east of it an ice base at
      !! -30 m vanishes k=8 to a filler and cuts k=7 into a 10 m partial
      !! cell.  `p_top` is the ISOSTATIC load evaluated at the MODEL column
      !! top, so the two columns are in exact hydrostatic rest and any
      !! acceleration is discretisation error.
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), allocatable, intent(out) :: b(:, :)
      integer, intent(out) :: iface
      real(wp), intent(out) :: de_surf
         !! Column-top offset across the staircase face (m) — the length
         !! the §3.1 surface trapezoid is taken over.
      real(wp) :: hopen(NZTEST), hice(NZTEST), eta_ice
      integer :: i, j, k, nxt, nyt
      nxt = size(ms%h_layer, 1)
      nyt = size(ms%h_layer, 2)
      allocate (b(nxt, nyt))
      hopen = HNOM
      do k = 1, 6
         hice(k) = HNOM
      end do
      hice(7) = 10.0_wp
      hice(8) = HFILL
      eta_ice = -160.0_wp + 6.0_wp*HNOM + 10.0_wp + HFILL
      de_surf = eta_ice
      iface = NGHOST + 3
      do j = 1, nyt
         do i = 1, nxt
            b(i, j) = 160.0_wp
            if (i < iface) then
               call seed_column(ms, i, j, -160.0_wp, hopen)
               ms%p_top(i, j) = 0.0_wp
            else
               call seed_column(ms, i, j, -160.0_wp, hice)
               ! pa_true(z) = 0.5*g*gamma*z^2 with the open-ocean datum;
               ! the kernel seeds pa(nz+1) = rho_ref*g*eta + p_top.
               ms%p_top(i, j) = 0.5_wp*GRAVITY*GAMMA*eta_ice**2 &
                                - RHO_REF*GRAVITY*eta_ice
            end if
         end do
      end do
      if (grid%nx_total < 1) return   ! keep `grid` referenced; never taken
   end subroutine build_ice_staircase

   ! ---------------------------------------------------------------------
   ! (1) Bed-side staircase
   ! ---------------------------------------------------------------------

   pure function bed_quadrature_error() result(pfu)
      !! `(g*gamma/12)*[h_L^3 - h_R^3 + de_top^3 - de_bot^3]`, divided by
      !! the face assembly's `rho0*dx*(h_L + h_R)/2`.  Here `h_L = 20`,
      !! `h_R = 10`, `de_top = 0` (the levels above the cut are aligned)
      !! and `de_bot = 10` (the bed step).
      real(wp) :: pfu
      pfu = (GRAVITY*GAMMA/12.0_wp)*(20.0_wp**3 - 10.0_wp**3 - 10.0_wp**3)
      pfu = pfu*2.0_wp/(RHO0*DX*30.0_wp)
   end function bed_quadrature_error

   subroutine test_bed_pcm_analytic(error)
      !! The PCM default carries the closed-form quadrature error, and this
      !! asserts the FORMULA — if the two ever disagree, either the face
      !! assembly changed or the derivation is wrong.
      type(error_type), allocatable, intent(out) :: error
      type(multilayer_state_t) :: ms
      type(hgrid_t) :: grid
      real(wp), allocatable :: b(:, :), pf(:, :, :)
      real(wp) :: got, want
      integer :: iface
      checks: block
         call grid%init(NXP, NYP, NGHOST, DX, DX)
         ms%nz_ml = NZTEST
         call ms%init(grid)
         call build_bed_staircase(grid, ms, b, iface)
         call run_pgf(grid, ms, b, .false., .false., 1, pf)
         got = pf(iface, NGHOST + 2, 2)
         want = bed_quadrature_error()
         call check(error, abs(got - want) < 1.0e-4_wp*abs(want), &
                    "bed staircase PCM: face PGF does not match "// &
                    "(g*gamma/12)[hL^3-hR^3-de_bot^3]")
      end block checks
      if (allocated(pf)) deallocate (pf)
      if (allocated(b)) deallocate (b)
      call ms%destroy()
   end subroutine test_bed_pcm_analytic

   subroutine bed_recon_core(error, scheme)
      !! With the reconstruction filler-aware, the partial-cell face and
      !! every aligned face above it are at round-off.  Before the fix the
      !! partial-cell face sat at -2.26e-08 (PLM and PPM) — 7 decades over
      !! the bound — because the bed filler's `hTr/H_VANISHED` mean entered
      !! the partial cell's edge build and its own `dpa`.
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: scheme
      type(multilayer_state_t) :: ms
      type(hgrid_t) :: grid
      real(wp), allocatable :: b(:, :), pf(:, :, :)
      real(wp) :: worst
      integer :: iface, k
      checks: block
         call grid%init(NXP, NYP, NGHOST, DX, DX)
         ms%nz_ml = NZTEST
         call ms%init(grid)
         call build_bed_staircase(grid, ms, b, iface)
         call run_pgf(grid, ms, b, .false., .true., scheme, pf)
         ! k = 1 is the one-side-vanished face the closed-face mask walls;
         ! it is excluded here and covered by its own bound below.
         worst = 0.0_wp
         do k = 2, NZTEST
            worst = max(worst, abs(pf(iface, NGHOST + 2, k)))
         end do
         call check(error, worst < REST_TOL, &
                    "bed staircase: reconstructed PGF is not at round-off "// &
                    "on the partial-cell column")
         if (allocated(error)) exit checks
         ! The filler-vs-live face is not a rest state (the two cells do
         ! not overlap in z at all), but with the filler flattened onto the
         ! live edge it must still be tiny -- it was 1.32e-03 before.
         call check(error, abs(pf(iface, NGHOST + 2, 1)) < 1.0e-9_wp, &
                    "bed staircase: the filler-vs-live face still carries "// &
                    "a large acceleration")
      end block checks
      if (allocated(pf)) deallocate (pf)
      if (allocated(b)) deallocate (b)
      call ms%destroy()
   end subroutine bed_recon_core

   subroutine test_bed_recon_plm(error)
      type(error_type), allocatable, intent(out) :: error
      call bed_recon_core(error, 1)
   end subroutine test_bed_recon_plm

   subroutine test_bed_recon_ppm(error)
      type(error_type), allocatable, intent(out) :: error
      call bed_recon_core(error, 2)
   end subroutine test_bed_recon_ppm

   ! ---------------------------------------------------------------------
   ! (2) Ice-base staircase
   ! ---------------------------------------------------------------------

   subroutine ice_recon_core(error, scheme)
      !! Every ALIGNED face (k = 1..6, 20 m <-> 20 m) must be at round-off.
      !! Before the fix these sat at 1.32e-08 — three decades ABOVE the PCM
      !! default's 1.22e-11 on the same state, because the top filler's
      !! drained mean rode the pa stack down into every layer below it.
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: scheme
      type(multilayer_state_t) :: ms
      type(hgrid_t) :: grid
      real(wp), allocatable :: b(:, :), pf(:, :, :)
      real(wp) :: worst, de_surf
      integer :: iface, k
      checks: block
         call grid%init(NXP, NYP, NGHOST, DX, DX)
         ms%nz_ml = NZTEST
         call ms%init(grid)
         call build_ice_staircase(grid, ms, b, iface, de_surf)
         call run_pgf(grid, ms, b, .true., .true., scheme, pf)
         worst = 0.0_wp
         do k = 1, 6
            worst = max(worst, abs(pf(iface, NGHOST + 2, k)))
         end do
         call check(error, worst < REST_TOL, &
                    "ice staircase: reconstructed PGF is not at round-off "// &
                    "on the aligned faces beneath the ice base")
      end block checks
      if (allocated(pf)) deallocate (pf)
      if (allocated(b)) deallocate (b)
      call ms%destroy()
   end subroutine ice_recon_core

   subroutine test_ice_recon_plm(error)
      type(error_type), allocatable, intent(out) :: error
      call ice_recon_core(error, 1)
   end subroutine test_ice_recon_plm

   subroutine test_ice_recon_ppm(error)
      type(error_type), allocatable, intent(out) :: error
      call ice_recon_core(error, 2)
   end subroutine test_ice_recon_ppm

   subroutine test_ice_partial_face_seed(error)
      !! What is LEFT at the partial-cell face under a sloping ice base is
      !! the Yung et al. (2026) §3.1 surface-pressure reconstruction, and
      !! nothing else.  `intx_pa(nz+1) = 0.5*(pa_L + pa_R)` is a trapezoid
      !! along an edge over which `pa` is QUADRATIC in z (linear
      !! stratification), so it carries `g*gamma*de_surf^2/12`; the Pass-3
      !! numerator picks it up with the weight `(h_R - h_L)`:
      !!
      !!     PFu = (h_R - h_L)*g*gamma*de_surf^2/12 * 2/(rho0*dx*(h_L+h_R))
      !!
      !! Asserted to 1 % — it is a closed form, not a fitted bound.  It is
      !! deliberately OUTSIDE this commit: the along-face `boole_dpa_face`
      !! quadrature re-anchors every INTERIOR interface but the surface
      !! seed is still the arithmetic mean, which is exact only where the
      !! column top is flat (it is, in the open ocean and on the bed-side
      !! staircase above — hence the round-off there).
      type(error_type), allocatable, intent(out) :: error
      type(multilayer_state_t) :: ms
      type(hgrid_t) :: grid
      real(wp), allocatable :: b(:, :), pf(:, :, :)
      real(wp) :: de_surf, want, got
      integer :: iface
      checks: block
         call grid%init(NXP, NYP, NGHOST, DX, DX)
         ms%nz_ml = NZTEST
         call ms%init(grid)
         call build_ice_staircase(grid, ms, b, iface, de_surf)
         call run_pgf(grid, ms, b, .true., .true., 1, pf)
         want = (10.0_wp - HNOM)*GRAVITY*GAMMA*de_surf*de_surf/12.0_wp
         want = want*2.0_wp/(RHO0*DX*(HNOM + 10.0_wp))
         got = pf(iface, NGHOST + 2, 7)
         call check(error, abs(got - want) < 1.0e-2_wp*abs(want), &
                    "ice staircase: the partial-cell residual is not the "// &
                    "surface-seed trapezoid g*gamma*de^2/12")
      end block checks
      if (allocated(pf)) deallocate (pf)
      if (allocated(b)) deallocate (b)
      call ms%destroy()
   end subroutine test_ice_partial_face_seed

   ! ---------------------------------------------------------------------
   ! (3) Flat lid over a filler stack
   ! ---------------------------------------------------------------------

   subroutine test_flat_lid_bit_zero(error)
      !! A UNIFORM draft gives every column the same stack — same filler
      !! run, same partial cell, same everything — so the Pass-3 numerator
      !! is a difference of IDENTICAL operands and the PGF is EXACTLY zero.
      !! `==` is legitimate here: it is `x - x`, not a cancellation between
      !! separately-computed quantities.  This is the gate that catches a
      !! filler rule which depends on anything other than the column.
      type(error_type), allocatable, intent(out) :: error
      type(multilayer_state_t) :: ms
      type(hgrid_t) :: grid
      real(wp), allocatable :: b(:, :), pf(:, :, :)
      real(wp) :: h(NZTEST), worst, eta_ice
      integer :: i, j, k, nxt, nyt
      checks: block
         call grid%init(NXP, NYP, NGHOST, DX, DX)
         ms%nz_ml = NZTEST
         call ms%init(grid)
         nxt = size(ms%h_layer, 1)
         nyt = size(ms%h_layer, 2)
         allocate (b(nxt, nyt))
         do k = 1, 6
            h(k) = HNOM
         end do
         h(7) = 10.0_wp
         h(8) = HFILL
         eta_ice = -160.0_wp + 6.0_wp*HNOM + 10.0_wp + HFILL
         do j = 1, nyt
            do i = 1, nxt
               b(i, j) = 160.0_wp
               call seed_column(ms, i, j, -160.0_wp, h)
               ms%p_top(i, j) = 0.5_wp*GRAVITY*GAMMA*eta_ice**2 &
                                - RHO_REF*GRAVITY*eta_ice
            end do
         end do
         call run_pgf(grid, ms, b, .true., .true., 2, pf)
         worst = 0.0_wp
         do k = 1, NZTEST
            do j = NGHOST + 1, NGHOST + NYP
               do i = NGHOST + 2, NGHOST + NXP
                  worst = max(worst, abs(pf(i, j, k)))
               end do
            end do
         end do
         call check(error, worst == 0.0_wp, &
                    "flat lid over fillers: the PGF is not exactly zero")
      end block checks
      if (allocated(pf)) deallocate (pf)
      if (allocated(b)) deallocate (b)
      call ms%destroy()
   end subroutine test_flat_lid_bit_zero

   ! ---------------------------------------------------------------------
   ! (4) Bit-identity on a column with no vanished layer
   ! ---------------------------------------------------------------------

   pure subroutine ref_boundary_edges(h_self, h_nbr, q_self, dq_up, q_t, q_b)
      !! Verbatim copy of `boundary_edges_linear` for the reference builders.
      real(wp), intent(in)  :: h_self, h_nbr, q_self, dq_up
      real(wp), intent(out) :: q_t, q_b
      real(wp), parameter :: H_TINY = 1.0e-30_wp
      real(wp) :: d
      d = dq_up*h_self/max(h_self + h_nbr, H_TINY)
      d = sign(min(abs(d), abs(dq_up)), d)
      q_t = q_self + d
      q_b = q_self - d
   end subroutine ref_boundary_edges

   pure subroutine ref_plm_edges(nz, h, q, q_t, q_b)
      !! The PLM edge build AS IT STOOD before the filler-awareness change:
      !! one unsegmented pass over the whole column.  Kept verbatim so the
      !! bit-identity assertion has something to be identical TO.
      integer, intent(in) :: nz
      real(wp), intent(in)  :: h(nz), q(nz)
      real(wp), intent(out) :: q_t(nz), q_b(nz)
      real(wp) :: slp(nz)
      real(wp) :: h_l, h_c, h_r, sig_c, sig_l, sig_r, slp_max
      real(wp) :: e_t, e_b, q_lo, q_hi
      integer  :: k
      if (nz <= 1) then
         q_t(1) = q(1)
         q_b(1) = q(1)
         return
      end if
      slp(1) = 0.0_wp
      slp(nz) = 0.0_wp
      do k = 2, nz - 1
         h_l = h(k - 1)
         h_c = h(k)
         h_r = h(k + 1)
         sig_l = q(k) - q(k - 1)
         sig_r = q(k + 1) - q(k)
         if (sig_l*sig_r <= 0.0_wp) then
            slp(k) = 0.0_wp
         else
            sig_c = 2.0_wp*(q(k + 1) - q(k - 1))*h_c/(h_l + 2.0_wp*h_c + h_r)
            slp_max = 2.0_wp*min(abs(sig_l), abs(sig_r))
            slp(k) = sign(min(abs(sig_c), slp_max), sig_c)
         end if
      end do
      call ref_boundary_edges(h(1), h(2), q(1), q(2) - q(1), q_t(1), q_b(1))
      call ref_boundary_edges(h(nz), h(nz - 1), q(nz), q(nz) - q(nz - 1), &
                              q_t(nz), q_b(nz))
      do k = 2, nz - 1
         e_t = q(k) + 0.5_wp*slp(k)
         e_b = q(k) - 0.5_wp*slp(k)
         q_lo = min(q(k), q(k + 1))
         q_hi = max(q(k), q(k + 1))
         q_t(k) = max(q_lo, min(q_hi, e_t))
         q_lo = min(q(k), q(k - 1))
         q_hi = max(q(k), q(k - 1))
         q_b(k) = max(q_lo, min(q_hi, e_b))
      end do
   end subroutine ref_plm_edges

   pure subroutine ref_ppm_edges(nz, h, q, q_t, q_b)
      !! The PPM edge build as it stood before the change — see
      !! `ref_plm_edges`.
      integer, intent(in) :: nz
      real(wp), intent(in)  :: h(nz), q(nz)
      real(wp), intent(out) :: q_t(nz), q_b(nz)
      real(wp) :: edge(nz)
      real(wp) :: q_lo, q_hi, ql, qr, qm, dq, dq_l, dq_r, q6
      real(wp) :: h0, h1, h2, h3, hf, h_sum
      real(wp) :: h01, h12, h23, h012, h123, h0123
      real(wp) :: f1, f2, f3, et1, et2, et3
      real(wp), parameter :: H_NEGLECT = 1.0e-30_wp
      real(wp), parameter :: H_MIN_FRAC = 1.0e-5_wp
      integer  :: k
      if (nz <= 1) then
         q_t(1) = q(1)
         q_b(1) = q(1)
         return
      end if
      if (nz == 2) then
         call ref_boundary_edges(h(1), h(2), q(1), q(2) - q(1), q_t(1), q_b(1))
         call ref_boundary_edges(h(2), h(1), q(2), q(2) - q(1), q_t(2), q_b(2))
         return
      end if
      do k = 2, nz - 2
         h0 = h(k - 1)
         h1 = h(k)
         h2 = h(k + 1)
         h3 = h(k + 2)
         h_sum = h0 + h1 + h2 + h3
         if (h0 + h1 <= 0.0_wp .or. h1 + h2 <= 0.0_wp .or. h2 + h3 <= 0.0_wp) then
            hf = H_MIN_FRAC*max(H_NEGLECT, h_sum)
            h0 = max(h0, hf)
            h1 = max(h1, hf)
            h2 = max(h2, hf)
            h3 = max(h3, hf)
         end if
         h01 = h0 + h1
         h12 = h1 + h2
         h23 = h2 + h3
         h012 = h0 + h1 + h2
         h123 = h1 + h2 + h3
         h0123 = h0 + h1 + h2 + h3
         f1 = h01*h23/h12
         f2 = h2*q(k) + h1*q(k + 1)
         f3 = 1.0_wp/h012 + 1.0_wp/h123
         et1 = f1*f2*f3
         et2 = (h2*h23/(h012*h01))*((h0 + 2.0_wp*h1)*q(k) - h1*q(k - 1))
         et3 = (h1*h01/(h123*h23))*((2.0_wp*h2 + h3)*q(k + 1) - h2*q(k + 2))
         edge(k) = (et1 + et2 + et3)/h0123
      end do
      edge(1) = (q(1)*h(2) + q(2)*h(1))/(h(1) + h(2))
      edge(nz - 1) = (q(nz - 1)*h(nz) + q(nz)*h(nz - 1))/(h(nz - 1) + h(nz))
      call ref_boundary_edges(h(1), h(2), q(1), q(2) - q(1), q_t(1), q_b(1))
      call ref_boundary_edges(h(nz), h(nz - 1), q(nz), q(nz) - q(nz - 1), &
                              q_t(nz), q_b(nz))
      do k = 2, nz - 1
         qm = q(k)
         ql = edge(k - 1)
         qr = edge(k)
         q_lo = min(q(k - 1), q(k), q(k + 1))
         q_hi = max(q(k - 1), q(k), q(k + 1))
         ql = max(q_lo, min(q_hi, ql))
         qr = max(q_lo, min(q_hi, qr))
         dq = qr - ql
         dq_l = qm - ql
         dq_r = qr - qm
         if (dq_l*dq_r <= 0.0_wp) then
            ql = qm
            qr = qm
         else
            q6 = 6.0_wp*qm - 3.0_wp*(ql + qr)
            if (abs(q6) > abs(dq)) then
               if (q6*dq > 0.0_wp) then
                  ql = 3.0_wp*qm - 2.0_wp*qr
               else
                  qr = 3.0_wp*qm - 2.0_wp*ql
               end if
            end if
         end if
         q_b(k) = ql
         q_t(k) = qr
      end do
   end subroutine ref_ppm_edges

   pure subroutine live_profile(icase, nz, h, q)
      !! Four all-live (`h` well above `H_VANISHED`) column shapes, with a
      !! non-monotone `q` in one of them so the limiter branches are
      !! exercised, not just the smooth path.
      integer, intent(in) :: icase, nz
      real(wp), intent(out) :: h(nz), q(nz)
      integer :: k
      real(wp) :: z
      select case (icase)
      case (1)                       ! uniform
         do k = 1, nz
            h(k) = 25.0_wp
         end do
      case (2)                       ! surface-refined
         do k = 1, nz
            h(k) = 5.0_wp + 8.0_wp*real(nz - k, wp)
         end do
      case (3)                       ! bed-refined, wide ratio
         do k = 1, nz
            h(k) = 0.5_wp + 30.0_wp*real(k - 1, wp)
         end do
      case default                   ! sawtooth thicknesses
         do k = 1, nz
            h(k) = 12.0_wp + 9.0_wp*real(modulo(k, 3), wp)
         end do
      end select
      z = 0.0_wp
      do k = 1, nz
         z = z - 0.5_wp*h(k)
         if (icase == 4) then
            ! Non-monotone: a mid-column inversion the limiters must see.
            q(k) = 34.0_wp - 1.2e-3_wp*z + 0.05_wp*sin(0.9_wp*real(k, wp))
         else
            q(k) = 34.0_wp - 1.2e-3_wp*z
         end if
         z = z - 0.5_wp*h(k)
      end do
   end subroutine live_profile

   subroutine bit_identical_core(error, parabolic)
      !! An all-live column must reconstruct BIT-IDENTICALLY to the
      !! unsegmented pre-change builder: the segmentation finds exactly one
      !! run, `1..nz`, and the filler pass is empty.  This is what keeps
      !! every sigma / z* answer — and the exactness numbers of the
      !! `fix/sigma-pgf-rest-state` work — from moving.
      type(error_type), allocatable, intent(out) :: error
      logical, intent(in) :: parabolic
      integer, parameter :: NZ_MAX = 24
      real(wp) :: h(NZ_MAX), q(NZ_MAX)
      real(wp) :: got_t(NZ_MAX), got_b(NZ_MAX)
      real(wp) :: ref_t(NZ_MAX), ref_b(NZ_MAX)
      integer :: icase, nz, k
      logical :: same
      same = .true.
      do icase = 1, 4
         do nz = 1, NZ_MAX
            call live_profile(icase, nz, h(1:nz), q(1:nz))
            if (parabolic) then
               call ppm_edges_column(nz, h(1:nz), q(1:nz), got_t(1:nz), got_b(1:nz))
               call ref_ppm_edges(nz, h(1:nz), q(1:nz), ref_t(1:nz), ref_b(1:nz))
            else
               call plm_edges_column(nz, h(1:nz), q(1:nz), got_t(1:nz), got_b(1:nz))
               call ref_plm_edges(nz, h(1:nz), q(1:nz), ref_t(1:nz), ref_b(1:nz))
            end if
            do k = 1, nz
               if (got_t(k) /= ref_t(k) .or. got_b(k) /= ref_b(k)) same = .false.
            end do
         end do
      end do
      call check(error, same, &
                 "all-live column: the filler-aware edge build is not "// &
                 "bit-identical to the unsegmented one")
   end subroutine bit_identical_core

   subroutine test_bit_identical_plm(error)
      type(error_type), allocatable, intent(out) :: error
      call bit_identical_core(error, .false.)
   end subroutine test_bit_identical_plm

   subroutine test_bit_identical_ppm(error)
      type(error_type), allocatable, intent(out) :: error
      call bit_identical_core(error, .true.)
   end subroutine test_bit_identical_ppm

end module test_ocean_pgf_filler_recon
