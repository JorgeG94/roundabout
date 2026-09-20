!! Top-of-column pressure in the EOS's IN-SITU pressure arguments (E3),
!! and the potential-density reference pressure that must NOT carry it.
!!
!! An in-situ EOS pressure on the ocean path used to be measured down from
!! **0 Pa at the free surface**.  Under an Antarctic ice shelf that is wrong
!! by the whole ice load (1e6-2e7 Pa), which with a nonlinear EOS is a
!! systematic several-kg/m^3 density error.  `multilayer_state_t%p_top` (Pa,
!! filled from the assembled `sf%p_surf` when `&ocean_psurf_nml in_eos`) is
!! the offset the ported IN-SITU builder now starts from.
!!
!! The load must reach IN-SITU pressures ONLY.  `ms%rho_layer` is a
!! POTENTIAL density at the horizontally uniform `&ocean_eos_nml p_ref`, and
!! its consumers difference it ALONG a layer (the Montgomery PGF, the
!! FV-lite / FV-MOM6-PCM integrands) and VERTICALLY (the vmix N^2 builders).
!! A per-column reference pressure would give two columns of IDENTICAL water
!! at the same geopotential depth densities differing by `dRho/dp * dp_top`
!! (~2 kg/m^3 across a calving front) — a large, entirely spurious
!! along-layer density gradient, hence a spurious pressure gradient force.
!! `rho_layer_independent_of_p_top` is the standing guard against
!! reintroducing exactly that.
!!
!! Cases:
!!   * `p_top_sign_and_monotone` — the MANDATORY sign test.  MOM6 handed its
!!     EOS a NEGATIVE pressure for years in one code path, so the pressure
!!     that actually reached the EOS is RECOVERED (bisection on the monotone
!!     `eos_density_point`, no coefficient duplication) and asserted
!!     `>= p_top > 0` and strictly increasing toward the bed (`k = 1`).
!!   * `p_top_column_sweep_analytic` — the FV-Wright Picard sweep's in-situ
!!     density equals `eos_density_point` at `p_top + p_above +
!!     0.5*g*rho_seed*h`, layer by layer, while its `p_edge` anomaly stack
!!     stays seeded at 0 (the PGF top boundary condition is deliberately NOT
!!     changed here).
!!   * `rho_layer_independent_of_p_top` — THE REGRESSION GUARD.  A
!!     horizontally VARYING `p_top` over a column of uniform (T, S, h) must
!!     leave `ms%rho_layer` exactly uniform.  Any future edit that routes
!!     the load into the potential-density reference fails here.
!!   * `p_ref_uniform_shifts_wright` — the now-live `&ocean_eos_nml p_ref`:
!!     `p_ref = 1e7` moves the Wright `rho_layer` to the scalar EOS value at
!!     1e7 (several kg/m^3 away from the surface answer), and `p_ref = 0` is
!!     BYTE-identical to the scalar EOS at 0.
!!   * `p_ref_linear_unchanged` — the linear EOS is pressure-independent, so
!!     `p_ref = 1e7` must be BYTE-identical to `p_ref = 0`.
!!   * `p_top_gate_off_leaves_p_top_zero` — knob-OFF bit-identity through the
!!     REAL gate: `ocean_dyn_step_split` with the psurf seam enabled and
!!     `sf%p_surf = 1e7 Pa` fills `ms%p_top` only when `in_eos = .true.`,
!!     and `ms%rho_layer` is BYTE-identical either way (the end-to-end
!!     statement that the load never reaches the potential density).
!!   * `p_top_validate_config` — the envelope: `in_eos` without `enable` and
!!     `in_eos` with a STILL-unported in-situ builder (kappa-shear) are
!!     REFUSED; `in_eos` with EPBL is ACCEPTED (ported in Phase 4b — see
!!     `test_ocean_bl_under_ice`); `in_eos` with a PGF form that has no
!!     in-situ pressure is ACCEPTED (documented inert, warning only), which
!!     is the deliberate no-op decision.
!!
!! GPU/mem:separate: every object and its scratch companion is
!! `enter_data`'d, host-set inputs (`p_top` included) are pushed with
!! `!$acc update device`, and read-back is `!$acc update self` on the
!! COMPONENT arrays, never the aggregate.  All directives are inert
!! no-ops on the host build — verify on the actual GPU build.
module test_ocean_eos_p_top
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, GRAVITY
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_continuity, only: continuity_t
   use rdb_coriolis_adv, only: coriolis_adv_t
   use rdb_ocean_pressure_force, only: ocean_pressure_force_t
   use rdb_ocean_horizontal_viscosity, only: ocean_horizontal_viscosity_t
   use rdb_ocean_bottom_drag, only: ocean_bottom_drag_t
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t
   use rdb_ocean_vertical_advection, only: ocean_vertical_advection_t
   use rdb_ocean_hdiff_tracer, only: ocean_hdiff_tracer_t
   use rdb_ocean_vdiff, only: ocean_vdiff_t
   use rdb_ocean_vmix, only: ocean_vmix_t
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t
   use rdb_ocean_dyn, only: ocean_dyn_t, ocean_dyn_step_split
   use rdb_ocean_p_surf, only: ocean_p_surf_t, p_surf_configure
   use rdb_eos, only: eos_t, eos_density_point, &
                      eos_wright_pgf_column_sweep_impl, &
                      EOS_VARIANT_LINEAR, EOS_VARIANT_WRIGHT_97
   use rdb_ocean_eos_compute, only: ocean_eos_compute
   use rdb_config, only: config_t, read_config_from_string, validate_config
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_CONFIG_VALIDATE
   implicit none
   private

   public :: collect_ocean_eos_p_top_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 6

   real(wp), parameter :: P_TOP = 1.0e7_wp
      !! Representative Antarctic ice-shelf load (~1000 m of ice).
   real(wp), parameter :: T_UNI = 1.5_wp
      !! degC — cold cavity water, well inside the Wright fit range.
   real(wp), parameter :: S_UNI = 34.6_wp
      !! PSU.
   real(wp), parameter :: DZ_LAYER = 200.0_wp
      !! Metres per layer in the analytical columns.

contains

   subroutine collect_ocean_eos_p_top_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("p_top_sign_and_monotone", test_sign_and_monotone), &
                  new_unittest("p_top_column_sweep_analytic", test_column_sweep_analytic), &
                  new_unittest("rho_layer_independent_of_p_top", test_rho_layer_independent), &
                  new_unittest("p_ref_uniform_shifts_wright", test_p_ref_shift), &
                  new_unittest("p_ref_linear_unchanged", test_p_ref_linear_unchanged), &
                  new_unittest("p_top_gate_off_leaves_p_top_zero", test_gate_off), &
                  new_unittest("p_top_validate_config", test_validate_config) &
                  ]
   end subroutine collect_ocean_eos_p_top_tests

   ! ------------------------------------------------------------------
   ! Helpers
   ! ------------------------------------------------------------------

   subroutine make_grid(grid)
      type(hgrid_t), intent(out) :: grid
      call grid%init(6, 5, NGHOST, 1000.0_wp, 1000.0_wp)
   end subroutine make_grid

   pure function recover_eos_pressure(eos, T, S, rho) result(p)
      !! Invert `rho = rho(T, S, p)` for `p` by bisection.  This is how the
      !! sign test sees the pressure a kernel ACTUALLY handed the EOS
      !! without exposing an internal or duplicating the Wright
      !! coefficients (a duplicated constant would pass the test while the
      !! kernel used a different one).  `rho` is strictly increasing in `p`
      !! for seawater, so bisection converges unconditionally; 200 halvings
      !! of [-1e9, 1e9] Pa reach far below the 1 Pa the assertions need.  A
      !! NEGATIVE pressure is representable in the bracket ON PURPOSE — the
      !! test must be able to SEE one.
      type(eos_t), intent(in) :: eos
      real(wp), intent(in) :: T, S, rho
      real(wp) :: p
      real(wp) :: lo, hi, mid
      integer :: it
      lo = -1.0e9_wp
      hi = 1.0e9_wp
      do it = 1, 200
         mid = 0.5_wp*(lo + hi)
         if (eos_density_point(eos, T, S, mid) < rho) then
            lo = mid
         else
            hi = mid
         end if
      end do
      p = 0.5_wp*(lo + hi)
   end function recover_eos_pressure

   subroutine build_uniform_column(grid, ms, eos)
      !! Uniform (T, S), uniform-thickness Wright column.  `ms%p_top` is
      !! left at the zero its `init` allocates it with; callers that want a
      !! load set it themselves BEFORE `compute_density` (which `copyin`s).
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(eos_t), intent(inout) :: eos
      ms%nz_ml = NZ
      call ms%init(grid)
      call eos%init(grid)
      eos%variant = EOS_VARIANT_WRIGHT_97
      ms%h_layer = DZ_LAYER
      ms%tracers(ms%idx_salinity)%hTr = S_UNI*DZ_LAYER
      ms%tracers(ms%idx_temperature)%hTr = T_UNI*DZ_LAYER
   end subroutine build_uniform_column

   subroutine compute_density(ms, eos)
      !! `mem:separate` contract: map the state, push the host-set `p_top`
      !! (mapped `copyin`, but pushed explicitly so the intent survives a
      !! future switch to `create`), run, pull `rho_layer` back as a
      !! COMPONENT array — never the aggregate.
      type(multilayer_state_t), intent(inout) :: ms
      type(eos_t), intent(in) :: eos
      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc update device(ms%p_top)
      call ocean_eos_compute(eos, ms)
      !$acc update self(ms%rho_layer)
      call ms%exit_data()
      !$acc exit data delete(ms)
   end subroutine compute_density

   subroutine run_column_sweep(grid, eos, rho_seed, p_top_val, p_edge, rho_insitu)
      !! Drive `eos_wright_pgf_column_sweep_impl` on a uniform column with a
      !! uniform `p_top`, under the `mem:separate` contract.
      type(hgrid_t), intent(in) :: grid
      type(eos_t), intent(in) :: eos
      real(wp), intent(in) :: rho_seed(:, :, :)
      real(wp), intent(in) :: p_top_val
      real(wp), allocatable, intent(out) :: p_edge(:, :, :), rho_insitu(:, :, :)
      real(wp), allocatable :: h_layer(:, :, :), hS(:, :, :), hT(:, :, :)
      real(wp), allocatable :: p_top_fld(:, :), seed(:, :, :)
      integer :: nx, ny

      nx = grid%nx_total
      ny = grid%ny_total
      allocate (h_layer(nx, ny, NZ), source=DZ_LAYER)
      allocate (hS(nx, ny, NZ), source=S_UNI*DZ_LAYER)
      allocate (hT(nx, ny, NZ), source=T_UNI*DZ_LAYER)
      allocate (seed, source=rho_seed)
      allocate (p_top_fld(nx, ny), source=p_top_val)
      allocate (p_edge(nx, ny, NZ + 1), source=-1.0_wp)
      allocate (rho_insitu(nx, ny, NZ), source=0.0_wp)

      !$acc enter data copyin(h_layer, hS, hT, seed, p_top_fld)
      !$acc enter data create(p_edge, rho_insitu)
      call eos_wright_pgf_column_sweep_impl(h_layer, hS, hT, seed, p_top_fld, &
                                            p_edge, rho_insitu, &
                                            GRAVITY, eos%rho0, nx, ny, NZ)
      !$acc update self(p_edge, rho_insitu)
      !$acc exit data delete(p_edge, rho_insitu)
      !$acc exit data delete(h_layer, hS, hT, seed, p_top_fld)
   end subroutine run_column_sweep

   ! ------------------------------------------------------------------
   ! Case 1 — sign + monotonicity of the in-situ pressure (mandatory)
   ! ------------------------------------------------------------------

   subroutine test_sign_and_monotone(error)
      !! The one ported in-situ builder must hand the EOS a pressure that is
      !! `>= p_top > 0` and STRICTLY increases toward the bed.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      real(wp), allocatable :: p_edge(:, :, :), rho_insitu(:, :, :)
      real(wp) :: p_rec, p_prev, min_p, worst_drop
      integer :: nx, ny, i, j, k, n_pairs

      checks: block
         call make_grid(grid)
         call build_uniform_column(grid, ms, eos)
         call compute_density(ms, eos)
         call run_column_sweep(grid, eos, ms%rho_layer, P_TOP, p_edge, rho_insitu)

         nx = grid%nx_total
         ny = grid%ny_total
         call check(error, P_TOP > 0.0_wp, &
                    "p_top must be positive for this test to mean anything")
         if (allocated(error)) exit checks

         ! `worst_drop` is the LARGEST (shallower - deeper) pressure
         ! difference over every adjacent pair.  Seeded at -huge, NOT at 0,
         ! so "every pair strictly increases downward" is `worst_drop < 0`
         ! rather than a check the seed would mask.
         min_p = huge(1.0_wp)
         worst_drop = -huge(1.0_wp)
         n_pairs = 0
         do j = 1, ny
            do i = 1, nx
               p_prev = -huge(1.0_wp)
               do k = NZ, 1, -1
                  p_rec = recover_eos_pressure(eos, T_UNI, S_UNI, rho_insitu(i, j, k))
                  min_p = min(min_p, p_rec)
                  if (p_prev > -huge(1.0_wp)) then
                     worst_drop = max(worst_drop, p_prev - p_rec)
                     n_pairs = n_pairs + 1
                  end if
                  p_prev = p_rec
               end do
            end do
         end do
         call check(error, n_pairs == nx*ny*(NZ - 1), &
                    "monotonicity check saw the wrong number of adjacent pairs")
         if (allocated(error)) exit checks
         ! The half-layer centre of the SURFACE layer already sits
         ! 0.5*g*rho*h below p_top, so the shallowest sample is strictly
         ! deeper than the load, not merely at it.
         call check(error, min_p > P_TOP, &
                    "FV-Wright column sweep handed the EOS a pressure at or "// &
                    "below p_top (it must be p_top + hydrostatic > 0)")
         if (allocated(error)) exit checks
         call check(error, worst_drop < 0.0_wp, &
                    "FV-Wright EOS pressure must STRICTLY increase toward the bed")
      end block checks

      call eos%destroy()
      call ms%destroy()
   end subroutine test_sign_and_monotone

   ! ------------------------------------------------------------------
   ! Case 2 — FV-Wright Picard sweep against the closed form
   ! ------------------------------------------------------------------

   subroutine test_column_sweep_analytic(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      real(wp), allocatable :: p_edge(:, :, :), rho_insitu(:, :, :)
      real(wp) :: p_above, p_seed, rho_expect, max_dev, max_edge
      integer :: nx, ny, i, j, k

      checks: block
         call make_grid(grid)
         call build_uniform_column(grid, ms, eos)
         call compute_density(ms, eos)
         call run_column_sweep(grid, eos, ms%rho_layer, P_TOP, p_edge, rho_insitu)

         nx = grid%nx_total
         ny = grid%ny_total
         max_dev = 0.0_wp
         do j = 1, ny
            do i = 1, nx
               p_above = 0.0_wp
               do k = NZ, 1, -1
                  p_seed = P_TOP + p_above + 0.5_wp*GRAVITY*ms%rho_layer(i, j, k)*DZ_LAYER
                  rho_expect = eos_density_point(eos, T_UNI, S_UNI, p_seed)
                  max_dev = max(max_dev, abs(rho_insitu(i, j, k) - rho_expect))
                  p_above = p_above + GRAVITY*rho_insitu(i, j, k)*DZ_LAYER
               end do
            end do
         end do
         call check(error, max_dev < 1.0e-9_wp, &
                    "FV-Wright in-situ density != eos_density_point at "// &
                    "p_top + hydrostatic")
         if (allocated(error)) exit checks

         ! SCOPE GUARD: p_top reaches the EOS argument ONLY.  The PGF
         ! anomaly stack must still be seeded at 0 at the top of the column
         ! -- the load's depth-uniform gradient is carried by the barotropic
         ! eta_forcing seam, and adding it here too would double-count.
         max_edge = maxval(abs(p_edge(:, :, NZ + 1)))
         call check(error, max_edge == 0.0_wp, &
                    "p_top must NOT enter the PGF p_edge top boundary condition")
      end block checks

      call eos%destroy()
      call ms%destroy()
   end subroutine test_column_sweep_analytic

   ! ------------------------------------------------------------------
   ! Case 3 — THE REGRESSION GUARD
   ! ------------------------------------------------------------------

   subroutine test_rho_layer_independent(error)
      !! `ms%rho_layer` is a POTENTIAL density at the horizontally uniform
      !! `eos%p_ref`.  Give a column of IDENTICAL water a `p_top` that ramps
      !! across the domain — the calving-front geometry — and the density
      !! must come out exactly uniform.
      !!
      !! If the load were ever routed into the potential-density reference,
      !! the two ends of this ramp would differ by `dRho/dp * 1e7` ~ 4.5
      !! kg/m^3 and every along-layer consumer (the Montgomery PGF, the
      !! FV-lite integrand) would see a density gradient in water that has
      !! none.  That is the bug this test exists to prevent.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(eos_t) :: eos
      real(wp) :: rho_uniform, ramp
      integer :: nx, ny, i, j

      checks: block
         call make_grid(grid)
         call build_uniform_column(grid, ms, eos)
         nx = grid%nx_total
         ny = grid%ny_total

         ! Sloping ice draft: 0 -> P_TOP across the domain in i.
         do j = 1, ny
            do i = 1, nx
               ramp = real(i - 1, wp)/real(nx - 1, wp)
               ms%p_top(i, j) = P_TOP*ramp
            end do
         end do
         call check(error, maxval(ms%p_top) - minval(ms%p_top) == P_TOP, &
                    "the p_top ramp is not spanning the intended range")
         if (allocated(error)) exit checks

         call compute_density(ms, eos)

         rho_uniform = eos_density_point(eos, T_UNI, S_UNI, eos%p_ref)
         call check(error, all(ms%rho_layer == rho_uniform), &
                    "ms%rho_layer moved with p_top: the surface load must NOT "// &
                    "reach the POTENTIAL-density reference pressure (it would "// &
                    "fabricate an along-layer density gradient)")
      end block checks

      call eos%destroy()
      call ms%destroy()
   end subroutine test_rho_layer_independent

   ! ------------------------------------------------------------------
   ! Case 4 — the now-live &ocean_eos_nml p_ref
   ! ------------------------------------------------------------------

   subroutine test_p_ref_shift(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms0, ms1
      type(eos_t) :: eos0, eos1
      real(wp) :: rho_ref0, rho_ref1, shift

      checks: block
         call make_grid(grid)
         call build_uniform_column(grid, ms0, eos0)
         call build_uniform_column(grid, ms1, eos1)
         eos1%p_ref = P_TOP
         call compute_density(ms0, eos0)
         call compute_density(ms1, eos1)

         rho_ref0 = eos_density_point(eos0, T_UNI, S_UNI, 0.0_wp)
         rho_ref1 = eos_density_point(eos1, T_UNI, S_UNI, P_TOP)

         call check(error, all(ms0%rho_layer == rho_ref0), &
                    "p_ref = 0 must be byte-identical to the scalar EOS at p = 0")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(ms1%rho_layer - rho_ref1)) < 1.0e-9_wp, &
                    "p_ref = 1e7 density != scalar eos_density_point(T, S, 1e7)")
         if (allocated(error)) exit checks

         ! A 1e7 Pa reference compresses seawater by several kg/m^3 — the
         ! thermobaric state the knob exists to select.  Bracket it rather
         ! than pin a digit.
         shift = rho_ref1 - rho_ref0
         call check(error, shift > 3.0_wp .and. shift < 10.0_wp, &
                    "a 1e7 Pa reference pressure must raise the Wright density "// &
                    "by several kg/m^3")
      end block checks

      call eos1%destroy()
      call eos0%destroy()
      call ms1%destroy()
      call ms0%destroy()
   end subroutine test_p_ref_shift

   subroutine test_p_ref_linear_unchanged(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms0, ms1
      type(eos_t) :: eos0, eos1

      checks: block
         call make_grid(grid)
         call build_uniform_column(grid, ms0, eos0)
         call build_uniform_column(grid, ms1, eos1)
         eos0%variant = EOS_VARIANT_LINEAR
         eos1%variant = EOS_VARIANT_LINEAR
         eos1%p_ref = P_TOP
         call compute_density(ms0, eos0)
         call compute_density(ms1, eos1)

         call check(error, all(ms1%rho_layer == ms0%rho_layer), &
                    "the linear EOS has no pressure dependence: p_ref must be "// &
                    "byte-identical to p_ref = 0")
      end block checks

      call eos1%destroy()
      call eos0%destroy()
      call ms1%destroy()
      call ms0%destroy()
   end subroutine test_p_ref_linear_unchanged

   ! ------------------------------------------------------------------
   ! Case 5 — the gate, end to end through the outer step
   ! ------------------------------------------------------------------

   subroutine test_gate_off(error)
      !! `sf%p_surf` non-zero, psurf seam ENABLED, over real outer steps on
      !! the DEFAULT (Montgomery) PGF — a form with no in-situ EOS pressure,
      !! so `in_eos` is the documented no-op there.  `ms%p_top` must be
      !! filled only with the knob on, and `ms%rho_layer` must be
      !! BYTE-identical either way: the end-to-end statement that the load
      !! never reaches the potential density.
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable :: rho_off(:, :, :), rho_on(:, :, :)
      real(wp), allocatable :: ptop_off(:, :), ptop_on(:, :)

      checks: block
         call run_step(error, in_eos=.false., rho_out=rho_off, ptop_out=ptop_off)
         if (allocated(error)) exit checks
         call run_step(error, in_eos=.true., rho_out=rho_on, ptop_out=ptop_on)
         if (allocated(error)) exit checks

         call check(error, maxval(abs(ptop_off)) == 0.0_wp, &
                    "in_eos=.false. must leave ms%p_top at exactly zero")
         if (allocated(error)) exit checks
         call check(error, minval(ptop_on) == P_TOP .and. maxval(ptop_on) == P_TOP, &
                    "in_eos=.true. must fill ms%p_top from sf%p_surf")
         if (allocated(error)) exit checks

         call check(error, all(rho_on == rho_off), &
                    "ms%rho_layer must be byte-identical with and without "// &
                    "in_eos: the potential density is referenced to the "// &
                    "uniform eos%p_ref, never to the surface load")
      end block checks
   end subroutine test_gate_off

   subroutine run_step(error, in_eos, rho_out, ptop_out)
      !! Minimal but REAL outer-step harness (shaped after
      !! `test_ocean_p_surf_bitident`): every slot the split driver touches
      !! is init'd, `enter_data`'d and `exit_data`'d, so the `mem:separate`
      !! contract holds for the p_top refresh too.
      type(error_type), allocatable, intent(out) :: error
      logical, intent(in) :: in_eos
      real(wp), allocatable, intent(out) :: rho_out(:, :, :)
      real(wp), allocatable, intent(out) :: ptop_out(:, :)
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(multilayer_state_t) :: ms
      type(continuity_t) :: ct
      type(coriolis_adv_t) :: cor
      type(ocean_pressure_force_t) :: pgf
      type(ocean_horizontal_viscosity_t) :: hv
      type(ocean_bottom_drag_t) :: bd
      type(ocean_surface_stress_t) :: ss
      type(ocean_vertical_advection_t) :: va
      type(ocean_hdiff_tracer_t) :: hd
      type(ocean_vdiff_t) :: vd
      type(ocean_vmix_t) :: vmix
      type(ocean_surface_flux_t) :: sf
      type(eos_t) :: eos
      type(ocean_dyn_t) :: dyn
      type(ocean_p_surf_t) :: psurf
      integer :: nx, ny, step
      integer, parameter :: N_STEPS = 3
      real(wp), parameter :: DT = 300.0_wp
      integer, parameter :: N_INNER = 8
      real(wp), parameter :: H_TOTAL = 1200.0_wp

      call grid%init(8, 8, NGHOST, 20000.0_wp, 20000.0_wp)
      nx = grid%nx_total
      ny = grid%ny_total
      call make_cartesian_metrics(metrics, grid)

      ms%nz_ml = NZ
      call ms%init(grid)
      call ct%init(grid, nz_ml=NZ)
      cor%f_0 = 1.0e-4_wp
      call cor%init(grid, nz_ml=NZ)
      call pgf%init(grid, nz_ml=NZ)
      call hv%init(grid, nz_ml=NZ)
      call bd%init(grid, nz_ml=NZ)
      call ss%init(grid, nz_ml=NZ)
      call va%init(grid, nz_ml=NZ)
      call hd%init(grid, nz_ml=NZ)
      call vd%init(grid, nz_ml=NZ)
      call vmix%init(grid, nz_ml=NZ)
      call sf%init(grid)
      call sf%set_components(grid, .true.)
      call eos%init(grid)
      eos%variant = EOS_VARIANT_WRIGHT_97
      call dyn%init(grid, nz_ml=NZ)
      dyn%bt_work%bt_H_ref = H_TOTAL

      ! Uniform (T, S) so the potential density has an exact closed form.
      ms%h_layer = H_TOTAL/real(NZ, wp)
      ms%tracers(ms%idx_salinity)%hTr = S_UNI*ms%h_layer
      ms%tracers(ms%idx_temperature)%hTr = T_UNI*ms%h_layer

      ! A real load on the assembled total, and the seam switched on so the
      ! OFF case is "the seam runs, the EOS does not see it" rather than
      ! "nothing happens".
      call sf%set_p_surf_const(P_TOP)
      psurf%enable = .true.
      psurf%in_eos = in_eos
      psurf%rho0 = eos%rho0
      call p_surf_configure(psurf, nx, ny)

      !$acc enter data copyin(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn, psurf, sf)
      call ms%enter_data()
      call ct%enter_data()
      call cor%enter_data()
      call pgf%enter_data()
      call hv%enter_data()
      call bd%enter_data()
      call ss%enter_data()
      call va%enter_data()
      call hd%enter_data()
      call vd%enter_data()
      call vmix%enter_data()
      call dyn%enter_data()
      call psurf%enter_data()
      call sf%enter_data()
      !$acc update device(ms%p_top, sf%p_surf)

      do step = 1, N_STEPS
         call ocean_dyn_step_split(grid, metrics, dyn, eos, cor, ct, pgf, hv, bd, ss, &
                                   va, hd, vd, vmix, ms, DT, N_INNER, sf=sf, psurf=psurf)
      end do

      !$acc update self(ms%rho_layer, ms%p_top)
      allocate (rho_out, source=ms%rho_layer)
      allocate (ptop_out, source=ms%p_top)

      call sf%exit_data()
      call psurf%exit_data()
      call dyn%exit_data()
      call vmix%exit_data()
      call vd%exit_data()
      call hd%exit_data()
      call va%exit_data()
      call ss%exit_data()
      call bd%exit_data()
      call hv%exit_data()
      call pgf%exit_data()
      call cor%exit_data()
      call ct%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn, psurf, sf)
      call metrics%exit_data()

      call check(error, all(rho_out == rho_out), "rho_out is NaN")

      call sf%destroy()
      call psurf%destroy()
      call dyn%destroy()
      call eos%destroy()
      call vmix%destroy()
      call vd%destroy()
      call hd%destroy()
      call va%destroy()
      call ss%destroy()
      call bd%destroy()
      call hv%destroy()
      call pgf%destroy()
      call cor%destroy()
      call ct%destroy()
      call ms%destroy()
   end subroutine run_step

   ! ------------------------------------------------------------------
   ! Case 6 — the envelope
   ! ------------------------------------------------------------------

   subroutine test_validate_config(error)
      !! The unported IN-SITU pressure builders must be REFUSED, not
      !! silently run against a loaded column on the old 0-Pa-at-the-surface
      !! convention.  A PGF form with NO in-situ pressure is ACCEPTED — the
      !! deliberate no-op-plus-warning decision, since the load's
      !! depth-uniform gradient already reaches the momentum through the
      !! barotropic `eta_forcing` seam.  `validate_config` with `ierr`
      !! present returns instead of aborting, which is what makes this
      !! testable.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      integer :: ierr

      checks: block
         ! Legal and ACTIVE: the ported in-situ consumer.
         call parse_case(cfg, '&ocean_psurf_nml enable = .true., in_eos = .true. /'// &
                         new_line("a")//'&ocean_pgf_nml form = "fv_wright" /')
         ierr = -999
         call validate_config(cfg, ierr=ierr)
         call check(error, ierr == OCEAN_STATUS_OK, &
                    "in_eos with the FV_WRIGHT PGF must be accepted")
         if (allocated(error)) exit checks

         ! Legal but INERT: no in-situ EOS pressure anywhere in this config.
         ! Accepted with a warning, NOT refused — see the docstring.
         call parse_case(cfg, "&ocean_psurf_nml enable = .true., in_eos = .true. /")
         ierr = -999
         call validate_config(cfg, ierr=ierr)
         call check(error, ierr == OCEAN_STATUS_OK, &
                    "in_eos with a PGF form that has no in-situ pressure must be "// &
                    "accepted as a documented no-op, not refused")
         if (allocated(error)) exit checks

         ! in_eos without the master switch: there is no p_surf to read.
         call parse_case(cfg, "&ocean_psurf_nml enable = .false., in_eos = .true. /")
         ierr = -999
         call validate_config(cfg, ierr=ierr)
         call check(error, ierr == OCEAN_STATUS_ERR_CONFIG_VALIDATE, &
                    "in_eos without enable must be refused")
         if (allocated(error)) exit checks

         ! in_eos with EPBL: PORTED in Phase 4b (`epbl_column_kernel`
         ! seeds its stack at `ms%p_top`), so it must now be ACCEPTED.
         ! The refusal that stood here is gone, and this is the assertion
         ! that keeps it gone.
         call parse_case(cfg, "&ocean_psurf_nml enable = .true., in_eos = .true. /"// &
                         new_line("a")//"&ocean_epbl_nml enable = .true. /")
         ierr = -999
         call validate_config(cfg, ierr=ierr)
         call check(error, ierr == OCEAN_STATUS_OK, &
                    "in_eos with EPBL must be ACCEPTED — its column pressure "// &
                    "stack is ported to ms%p_top")
         if (allocated(error)) exit checks

         ! ... and a builder that is STILL unported is still refused, so
         ! the lift above is scoped to EPBL and did not empty the list.
         call parse_case(cfg, "&ocean_psurf_nml enable = .true., in_eos = .true. /"// &
                         new_line("a")//"&ocean_kappa_shear_nml enable = .true. /")
         ierr = -999
         call validate_config(cfg, ierr=ierr)
         call check(error, ierr == OCEAN_STATUS_ERR_CONFIG_VALIDATE, &
                    "in_eos with kappa-shear (still an unported in-situ "// &
                    "pressure builder) must be refused")
      end block checks
   end subroutine test_validate_config

   subroutine parse_case(cfg, extra)
      type(config_t), intent(out) :: cfg
      character(len=*), intent(in) :: extra
      character(len=:), allocatable :: nml
      nml = '&sim_nml sim_type = "ocean" /'//new_line("a")// &
            "&grid_nml nx = 8, ny = 8, dx = 2000.0, dy = 2000.0 /"//new_line("a")// &
            "&nonhydrostatic_nml nz_layers = 3 /"//new_line("a")// &
            "&time_nml t_end = 3600.0, dt_fixed = 300.0 /"//new_line("a")// &
            "&ocean_forcing_nml enable_components = .true. /"//new_line("a")// &
            "&ocean_bt_nml n_inner = 8 /"//new_line("a")// &
            extra//new_line("a")
      call read_config_from_string(nml, cfg)
   end subroutine parse_case

end module test_ocean_eos_p_top
