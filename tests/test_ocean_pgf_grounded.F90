!! Analytical tests for the grounded-layer PGF gate
!! (`&ocean_isopycnal_nml pgf_skip_nonoverlap` →
!! `ocean_pressure_force_t%skip_nonoverlap`).
!!
!! The defect these pin
!! --------------------
!! The face PGF is a two-point Jacobian,
!!
!!     PGF_x = -(1/rho0)*[ (p_c^R - p_c^L)/dx + g*rho_face*(z_c^R - z_c^L)/dx ]
!!
!! whose two terms cancel AT REST only while the two abutting layer centres
!! lie in a common z-interval whose AMBIENT density is the layer's own
!! `rho_layer`.  In an isopycnal (VCOORD_LAGRANGIAN) column a layer that has
!! grounded against sloping topography is squeezed onto the `angstrom_h` floor
!! on the shallow side while remaining massive one cell away: the two centres
!! are then hundreds of metres apart in z, the interval between them is filled
!! with OTHER density classes, and the cancellation leaves
!!
!!     g*(rho_layer - rho_ambient)*dz/dx
!!
!! of pressure-gradient acceleration on a MOTIONLESS ocean.  Nothing bounds it,
!! so it spins a basin-scale current up out of nothing — the shipped quiescent
!! case `validation_examples/ocean/seamount/seamount_conservative_floor.nml`
!! reached 6.7 cm/s (En 2.24e-3) and was still growing when the run stopped.
!!
!! The fix asserted here: where the layer's z-extents in the two abutting
!! columns do NOT overlap there is no common depth to difference the pressure
!! across, so the honest face value is zero.
!!
!! Variants
!! --------
!! The gate is variant-agnostic — it is a statement about GEOMETRY, not about
!! the quadrature — so every case below runs twice, once for `fv_lite` (the
!! layer-centre two-point Jacobian) and once for `fv_mom6` (the faithful
!! layer-integrated FV-Bouss form).  `fv_mom6` reads the bathymetry, so every
!! case sets `pgf%b` to the column-integrated thickness, which puts the free
!! surface at z = 0 and makes the seeded state a genuine rest state for it
!! too; `fv_lite` measures z from the surface and ignores `b`, so the same
!! call is inert there.
!!
!! Cases
!! -----
!!   * `flat_isopycnals_over_slope_are_at_rest` — the regression proper.  Seed
!!     the exact `thickness_config="uniform_z"` resting state (flat isopycnals,
!!     sub-bathymetry layers collapsed onto the floor) over a sloping bed with
!!     a layered density stack.  That state is an EXACT rest state: pressure is
!!     horizontally uniform at every depth.  Ungated, the PGF is large and
!!     O(g*drho*db/dx); gated, it must be zero to round-off in EVERY layer.
!!   * `gate_is_inert_on_aligned_columns` — no false positives.  Columns whose
!!     layers align in z must be bit-identical with the gate on and off.
!!   * `gate_preserves_a_real_density_gradient` — the gate must not eat physics:
!!     a genuine horizontal density contrast in fully overlapping layers still
!!     produces the same PGF with the gate on.
module test_ocean_pgf_grounded
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_ocean_pressure_force, only: ocean_pressure_force_t, &
                                       ocean_pressure_force_compute, &
                                       OPGF_VARIANT_FV_LITE, OPGF_VARIANT_FV_MOM6
   implicit none
   private

   public :: collect_ocean_pgf_grounded_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 8
   real(wp), parameter :: MAX_DEPTH = 2000.0_wp
      !! Basin depth the uniform-z interfaces are laid over.
   real(wp), parameter :: RHO_LIGHTEST = 1035.0_wp
   real(wp), parameter :: RHO_RANGE = 2.0_wp
   real(wp), parameter :: ANGSTROM_H = 1.0e-2_wp
      !! Same floor the shipped `seamount_conservative_floor.nml` uses.

contains

   subroutine collect_ocean_pgf_grounded_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("flat_isopycnals_over_slope_are_at_rest", test_grounded_rest_fv_lite), &
                  new_unittest("gate_is_inert_on_aligned_columns", test_gate_inert_fv_lite), &
                  new_unittest("gate_preserves_a_real_density_gradient", test_keeps_gradient_fv_lite), &
                  new_unittest("fv_mom6_flat_isopycnals_over_slope_are_at_rest", &
                               test_grounded_rest_fv_mom6), &
                  new_unittest("fv_mom6_gate_is_inert_on_aligned_columns", test_gate_inert_fv_mom6), &
                  new_unittest("fv_mom6_gate_preserves_a_real_density_gradient", &
                               test_keeps_gradient_fv_mom6) &
                  ]
   end subroutine collect_ocean_pgf_grounded_tests

   subroutine set_b_from_column(ms, pgf)
      !! Bathymetry consistent with the seeded column: `b = Σ_k h_layer`, so
      !! the FV_MOM6 interface stack `e(1) = -b`, `e(k+1) = e(k) + h` puts the
      !! free surface at exactly z = 0 in every column.  Without it `pgf%b`
      !! stays at its zero default, the whole water column sits ABOVE z = 0,
      !! and `pa(nz+1) = rho_ref·g·η` injects a free-surface slope that has
      !! nothing to do with the defect under test.  `fv_lite` never reads `b`
      !! (it measures z down from the surface), so the call is inert there and
      !! the two variants stay directly comparable.
      !! Must run BEFORE `enter_data` — `set_bathymetry` issues no
      !! `update device` of its own.
      type(multilayer_state_t), intent(in) :: ms
      type(ocean_pressure_force_t), intent(inout) :: pgf
      real(wp), allocatable :: b(:, :)
      integer :: i, j
      allocate (b(size(ms%h_layer, 1), size(ms%h_layer, 2)))
      do j = 1, size(b, 2)
         do i = 1, size(b, 1)
            b(i, j) = sum(ms%h_layer(i, j, :))
         end do
      end do
      call pgf%set_bathymetry(b)
      deallocate (b)
   end subroutine set_b_from_column

   subroutine run_pgf(ms, pgf, dx)
      !! One compute pass, with the face buffers pulled back to the host.
      !! `dpdx/dpdy` are scratch (create/delete), so the `!$acc update self`
      !! before `exit_data` is required for the host-side comparison — the
      !! GPU build maps `mem:separate` and would otherwise compare stale host
      !! memory (see `test_ocean_pgf_fv`).
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_pressure_force_t), intent(inout) :: pgf
      real(wp), intent(in) :: dx
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      call grid%init(size(ms%h_layer, 1) - 2*NGHOST, &
                     size(ms%h_layer, 2) - 2*NGHOST, &
                     NGHOST, dx, dx)
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

   pure subroutine seed_uniform_z_column(depth, h_col)
      !! `thickness_config = "uniform_z"` for ONE column, mirroring
      !! `seed_h_layer_uniform_z_impl`: uniform z interfaces over the GLOBAL
      !! `MAX_DEPTH`, clipped bottom-up against the local bed, with whatever
      !! will not fit collapsed onto `ANGSTROM_H`.  The telescoping sum is
      !! exactly `depth`, and the surviving interfaces are FLAT across columns
      !! — which is what makes the seeded state a genuine rest state.
      real(wp), intent(in) :: depth
      real(wp), intent(out) :: h_col(NZ)
      integer :: k
      real(wp) :: z_bot, z_top, z_top_target
      z_bot = -depth
      do k = 1, NZ
         z_top_target = -MAX_DEPTH*real(NZ - k, wp)/real(NZ, wp)
         z_top = max(z_top_target, z_bot + ANGSTROM_H)
         z_top = min(z_top, -real(NZ - k, wp)*ANGSTROM_H)
         h_col(k) = z_top - z_bot
         z_bot = z_top
      end do
   end subroutine seed_uniform_z_column

   subroutine build_grounded_state(ms, grid, bed_shallow)
      !! Stratified isopycnal stack at rest over a bed that shoals linearly in
      !! i from `MAX_DEPTH` to `bed_shallow`.  `k=1` is the bed (heaviest),
      !! `k=NZ` the surface (lightest) — Roundabout's bottom-up convention.
      type(multilayer_state_t), intent(inout) :: ms
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: bed_shallow
      integer :: i, j, k, nx, ny
      real(wp) :: h_col(NZ)
      real(wp) :: depth, frac
      nx = grid%nx_total
      ny = grid%ny_total
      do j = 1, ny
         do i = 1, nx
            frac = real(i - 1, wp)/real(nx - 1, wp)
            depth = MAX_DEPTH + frac*(bed_shallow - MAX_DEPTH)
            call seed_uniform_z_column(depth, h_col)
            do k = 1, NZ
               ms%h_layer(i, j, k) = h_col(k)
            end do
         end do
      end do
      ! Layer density is the coordinate: fixed per layer, horizontally uniform.
      do k = 1, NZ
         ms%rho_layer(:, :, k) = RHO_LIGHTEST + RHO_RANGE* &
                                 (real(NZ - k, wp) + 0.5_wp)/real(NZ, wp)
      end do
   end subroutine build_grounded_state

   ! -----------------------------------------------------------------

   subroutine test_grounded_rest_fv_lite(error)
      type(error_type), allocatable, intent(out) :: error
      call test_grounded_rest(error, OPGF_VARIANT_FV_LITE)
   end subroutine test_grounded_rest_fv_lite

   subroutine test_grounded_rest_fv_mom6(error)
      type(error_type), allocatable, intent(out) :: error
      call test_grounded_rest(error, OPGF_VARIANT_FV_MOM6)
   end subroutine test_grounded_rest_fv_mom6

   subroutine test_grounded_rest(error, variant)
      !! THE regression. Flat isopycnals over a seamount-scale slope: pressure
      !! is horizontally uniform at every depth, so every massive layer must
      !! feel exactly zero force. Ungated the kernel reports a large spurious
      !! gradient on the grounded faces; gated it must report zero everywhere.
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: variant
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_pressure_force_t) :: pgf_open, pgf_gated
      real(wp) :: max_open, max_gated
      checks: block

         call grid%init(32, 6, NGHOST, 8000.0_wp, 8000.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call pgf_open%init(grid, nz_ml=NZ)
         call pgf_gated%init(grid, nz_ml=NZ)
         pgf_open%variant = variant
         pgf_gated%variant = variant
         pgf_gated%skip_nonoverlap = .true.

         call build_grounded_state(ms, grid, bed_shallow=100.0_wp)
         call set_b_from_column(ms, pgf_open)
         call set_b_from_column(ms, pgf_gated)

         call run_pgf(ms, pgf_open, 8000.0_wp)
         max_open = maxval(abs(pgf_open%dpdx_face%data))

         call run_pgf(ms, pgf_gated, 8000.0_wp)
         max_gated = maxval(abs(pgf_gated%dpdx_face%data))

         ! Guard the setup itself: if the ungated kernel is already quiet the
         ! test would pass vacuously and stop protecting anything.
         call check(error, max_open > 1.0e-5_wp, &
                    "setup produced no spurious grounded-layer PGF to gate — "// &
                    "the regression would pass vacuously")
         if (allocated(error)) exit checks

         ! A resting ocean. Not "small": zero.
         call check(error, max_gated < 1.0e-14_wp, &
                    "flat isopycnals over a slope are a rest state, but the "// &
                    "gated PGF is non-zero — a quiescent isopycnal seamount "// &
                    "will spin up a spurious current")

      end block checks
      call pgf_gated%destroy()
      call pgf_open%destroy()
      call ms%destroy()
   end subroutine test_grounded_rest

   subroutine test_gate_inert_fv_lite(error)
      type(error_type), allocatable, intent(out) :: error
      call test_gate_inert(error, OPGF_VARIANT_FV_LITE)
   end subroutine test_gate_inert_fv_lite

   subroutine test_gate_inert_fv_mom6(error)
      type(error_type), allocatable, intent(out) :: error
      call test_gate_inert(error, OPGF_VARIANT_FV_MOM6)
   end subroutine test_gate_inert_fv_mom6

   subroutine test_gate_inert(error, variant)
      !! No false positives: where every layer is present in both columns the
      !! z-extents overlap, the gate must never fire, and the two results must
      !! agree BITWISE (not merely closely) — that is the property that keeps
      !! every non-Lagrangian vertical coordinate bit-identical.
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: variant
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_pressure_force_t) :: pgf_open, pgf_gated
      real(wp), allocatable :: dpdx_open(:, :, :), dpdy_open(:, :, :)
      integer :: i, j, k, nx, ny
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp) :: wiggle
      checks: block

         call grid%init(24, 6, NGHOST, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call pgf_open%init(grid, nz_ml=NZ)
         call pgf_gated%init(grid, nz_ml=NZ)
         pgf_open%variant = variant
         pgf_gated%variant = variant
         pgf_gated%skip_nonoverlap = .true.
         nx = grid%nx_total
         ny = grid%ny_total

         ! Thick layers that undulate mildly: centres shift between columns but
         ! the layers still overlap generously in z, as under sigma / z*.
         do j = 1, ny
            do i = 1, nx
               wiggle = 2.0_wp*sin(2.0_wp*PI*real(i, wp)/real(nx, wp))
               do k = 1, NZ
                  ms%h_layer(i, j, k) = 50.0_wp + wiggle*real(mod(k, 2), wp)
               end do
            end do
         end do
         do k = 1, NZ
            ms%rho_layer(:, :, k) = RHO_LIGHTEST + RHO_RANGE* &
                                    (real(NZ - k, wp) + 0.5_wp)/real(NZ, wp)
         end do

         call set_b_from_column(ms, pgf_open)
         call set_b_from_column(ms, pgf_gated)

         call run_pgf(ms, pgf_open, 1.0_wp)
         allocate (dpdx_open, source=pgf_open%dpdx_face%data)
         allocate (dpdy_open, source=pgf_open%dpdy_face%data)

         call run_pgf(ms, pgf_gated, 1.0_wp)

         call check(error, all(pgf_gated%dpdx_face%data == dpdx_open), &
                    "gate fired on overlapping layers: dpdx is not bit-identical")
         if (allocated(error)) exit checks
         call check(error, all(pgf_gated%dpdy_face%data == dpdy_open), &
                    "gate fired on overlapping layers: dpdy is not bit-identical")

      end block checks
      if (allocated(dpdx_open)) deallocate (dpdx_open)
      if (allocated(dpdy_open)) deallocate (dpdy_open)
      call pgf_gated%destroy()
      call pgf_open%destroy()
      call ms%destroy()
   end subroutine test_gate_inert

   subroutine test_keeps_gradient_fv_lite(error)
      type(error_type), allocatable, intent(out) :: error
      call test_keeps_gradient(error, OPGF_VARIANT_FV_LITE)
   end subroutine test_keeps_gradient_fv_lite

   subroutine test_keeps_gradient_fv_mom6(error)
      type(error_type), allocatable, intent(out) :: error
      call test_keeps_gradient(error, OPGF_VARIANT_FV_MOM6)
   end subroutine test_keeps_gradient_fv_mom6

   subroutine test_keeps_gradient(error, variant)
      !! The gate must not eat real physics. Give one fully-overlapping layer a
      !! genuine horizontal density contrast: the response must survive the
      !! gate unchanged (and be large enough that "unchanged" means something).
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: variant
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_pressure_force_t) :: pgf_open, pgf_gated
      real(wp), allocatable :: dpdx_open(:, :, :)
      integer :: i, j, k, nx, ny
      real(wp) :: max_signal
      checks: block

         call grid%init(24, 6, NGHOST, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         call pgf_open%init(grid, nz_ml=NZ)
         call pgf_gated%init(grid, nz_ml=NZ)
         pgf_open%variant = variant
         pgf_gated%variant = variant
         pgf_gated%skip_nonoverlap = .true.
         nx = grid%nx_total
         ny = grid%ny_total

         ms%h_layer = 50.0_wp
         do k = 1, NZ
            ms%rho_layer(:, :, k) = RHO_LIGHTEST
         end do
         ! A real front in the bed-most layer: denser to the west.
         do j = 1, ny
            do i = 1, nx
               ms%rho_layer(i, j, 1) = RHO_LIGHTEST + &
                                       1.0_wp*real(nx - i, wp)/real(nx - 1, wp)
            end do
         end do

         call set_b_from_column(ms, pgf_open)
         call set_b_from_column(ms, pgf_gated)

         call run_pgf(ms, pgf_open, 1.0_wp)
         allocate (dpdx_open, source=pgf_open%dpdx_face%data)
         max_signal = maxval(abs(dpdx_open))

         call run_pgf(ms, pgf_gated, 1.0_wp)

         call check(error, max_signal > 1.0e-4_wp, &
                    "density front produced no PGF — setup is not exercising "// &
                    "the response the gate must preserve")
         if (allocated(error)) exit checks
         call check(error, all(pgf_gated%dpdx_face%data == dpdx_open), &
                    "the gate altered the response to a real density gradient")

      end block checks
      if (allocated(dpdx_open)) deallocate (dpdx_open)
      call pgf_gated%destroy()
      call pgf_open%destroy()
      call ms%destroy()
   end subroutine test_keeps_gradient

end module test_ocean_pgf_grounded
