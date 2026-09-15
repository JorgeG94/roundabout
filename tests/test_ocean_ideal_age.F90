!! Unit tests for the ideal-age kernels — the MOM6 USE_IDEAL_AGE_TRACER
!! analogue, split (PR-7) into an interior-aging entry point
!! (`ocean_ideal_age_age_step`) and a surface-reset entry point
!! (`ocean_ideal_age_reset_step`), plus the host-scalar
!! `ocean_ideal_age_young_val` evaluator.
!!
!! Cases:
!!   * Aging from zero: after one step, interior age = dt; surface
!!     (untouched by the age-only kernel) stays 0.
!!   * Surface reset: pre-existing top-layer age gets overwritten to
!!     `young_eff * h_layer` regardless of its initial value.
!!   * Linear-in-time: after N steps with no flow, interior age = N·dt.
!!   * Interior conservation: surface reset doesn't touch subsurface
!!     layers (bottom-up-convention guard — the single highest-
!!     probability transcription bug against the top-down MOM6
!!     reference).
!!   * young_val sets a concentration*thickness, not a raw
!!     concentration — the single most likely units bug in the split
!!     kernel.
!!   * The growth-rate evaluator is exponential in `t` with the right
!!     sign, and short-circuits to bitwise-zero `exp`-free at the
!!     shipped defaults.
!!   * Tracer registry: `with_ideal_age=.true.` registers idx_age=3
!!     with the expected metadata; default init leaves idx_age=0.
module test_ocean_ideal_age
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_ideal_age, only: ocean_ideal_age_age_step, ocean_ideal_age_reset_step, &
                                  ocean_ideal_age_young_val
   implicit none
   private

   public :: collect_ocean_ideal_age_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NX_PHYS = 4
   integer, parameter :: NY_PHYS = 4
   integer, parameter :: NZ = 5
   real(wp), parameter :: H_EACH = 100.0_wp
   real(wp), parameter :: DT = 1800.0_wp

contains

   subroutine collect_ocean_ideal_age_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("ideal_age_init_registers_at_index_3", &
                               test_init_registers), &
                  new_unittest("ideal_age_init_default_off_leaves_idx_zero", &
                               test_init_default_off), &
                  new_unittest("ideal_age_one_step_from_zero_gives_dt", &
                               test_one_step_from_zero), &
                  new_unittest("ideal_age_surface_layer_reset", &
                               test_surface_reset), &
                  new_unittest("ideal_age_n_steps_linear", &
                               test_n_steps_linear), &
                  new_unittest("ideal_age_surface_reset_isolates_top", &
                               test_surface_reset_isolates_top), &
                  new_unittest("ideal_age_young_val_sets_surface", &
                               test_young_val_sets_surface), &
                  new_unittest("ideal_age_growth_rate_is_exponential", &
                               test_growth_rate_is_exponential), &
                  new_unittest("ideal_age_default_knobs_short_circuit", &
                               test_default_knobs_short_circuit), &
                  new_unittest("age_index_survives_passive_registration", &
                               test_age_survives_passive_registration) &
                  ]
   end subroutine collect_ocean_ideal_age_tests

   subroutine make_grid(grid)
      type(hgrid_t), intent(out) :: grid
      call grid%init(NX_PHYS, NY_PHYS, NGHOST, 1000.0_wp, 1000.0_wp)
   end subroutine make_grid

   subroutine setup(grid, ms, with_age)
      type(hgrid_t), intent(inout) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      logical, intent(in) :: with_age
      call make_grid(grid)
      ms%nz_ml = NZ
      call ms%init(grid, with_ideal_age=with_age)
      ms%h_layer = H_EACH
   end subroutine setup

   subroutine test_init_registers(error)
      !! When with_ideal_age=.true., the registry holds 3 tracers and
      !! idx_age points at the LAST slot (not the literal 3 — PR-28
      !! generalises this to "whichever slot registration produced",
      !! since a passive tracer registered after age would otherwise
      !! make a literal `== 3` assertion pass for the wrong reason).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      call setup(grid, ms, with_age=.true.)
      call check(error, ms%idx_age == size(ms%tracers), &
                 "idx_age should be the last registry slot when ideal age registered")
      if (.not. allocated(error)) call check(error, size(ms%tracers) == 3, &
                                             "tracer registry should have 3 entries")
      if (.not. allocated(error)) call check(error, &
                                             trim(ms%tracers(ms%idx_age)%name) == "age", &
                                             "age tracer name should be 'age'")
      if (.not. allocated(error)) call check(error, &
                                             trim(ms%tracers(ms%idx_age)%units) == "s", &
                                             "age tracer units should be 's'")
      call ms%destroy()
   end subroutine test_init_registers

   subroutine test_init_default_off(error)
      !! Without with_ideal_age, idx_age stays 0 and registry has 2.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      call setup(grid, ms, with_age=.false.)
      call check(error, ms%idx_age == 0, "idx_age should be 0 by default")
      if (.not. allocated(error)) call check(error, size(ms%tracers) == 2, &
                                             "tracer registry should be S + T only")
      call ms%destroy()
   end subroutine test_init_default_off

   subroutine test_one_step_from_zero(error)
      !! All hTr_age starts at 0.  After one interior-aging-only step:
      !!   * Subsurface (k < nz): hTr = dt · h_layer
      !!   * Surface (k = nz): untouched (still 0) — the age-step
      !!     kernel now owns k=1..nz-1 exclusively; the reset is a
      !!     separate kernel (PR-7).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      real(wp) :: expected_sub, actual_sub, surface_val
      integer :: nx, ny

      call setup(grid, ms, with_age=.true.)
      ms%tracers(ms%idx_age)%hTr = 0.0_wp
      nx = grid%nx_total
      ny = grid%ny_total

      call ocean_ideal_age_age_step(ms%tracers(ms%idx_age)%hTr, ms%h_layer, &
                                    DT, nx, ny, NZ)

      expected_sub = DT*H_EACH
      actual_sub = ms%tracers(ms%idx_age)%hTr(nx/2, ny/2, 1)
      call check(error, abs(actual_sub - expected_sub) < 1.0e-12_wp, &
                 "subsurface age after one step should be dt · h_layer")

      if (.not. allocated(error)) then
         surface_val = maxval(abs(ms%tracers(ms%idx_age)%hTr(:, :, NZ)))
         call check(error, surface_val < 1.0e-12_wp, &
                    "surface layer must be untouched by the age-only kernel")
      end if

      call ms%destroy()
   end subroutine test_one_step_from_zero

   subroutine test_surface_reset(error)
      !! Pre-load the surface layer with a large positive value; the
      !! reset kernel (young_eff=0) must zero it.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      real(wp) :: surface_val
      integer :: nx, ny

      call setup(grid, ms, with_age=.true.)
      nx = grid%nx_total
      ny = grid%ny_total
      ms%tracers(ms%idx_age)%hTr = 0.0_wp
      ms%tracers(ms%idx_age)%hTr(:, :, NZ) = 5.0e8_wp   ! garbage to clear

      call ocean_ideal_age_reset_step(ms%tracers(ms%idx_age)%hTr, ms%h_layer, &
                                      0.0_wp, nx, ny, NZ)

      surface_val = maxval(abs(ms%tracers(ms%idx_age)%hTr(:, :, NZ)))
      call check(error, surface_val < 1.0e-12_wp, &
                 "surface layer must be reset to zero even from large IC")
      call ms%destroy()
   end subroutine test_surface_reset

   subroutine test_n_steps_linear(error)
      !! No advection / diffusion / vertical exchange in this unit
      !! test, so subsurface concentration accumulates linearly:
      !! age(N·dt) = N · dt.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      integer, parameter :: N_STEPS = 7
      real(wp) :: conc, expected
      integer :: n_iter, nx, ny

      call setup(grid, ms, with_age=.true.)
      ms%tracers(ms%idx_age)%hTr = 0.0_wp
      nx = grid%nx_total
      ny = grid%ny_total
      do n_iter = 1, N_STEPS
         call ocean_ideal_age_age_step(ms%tracers(ms%idx_age)%hTr, ms%h_layer, &
                                       DT, nx, ny, NZ)
      end do

      ! Concentration in seconds = hTr / h_layer at any subsurface cell.
      conc = ms%tracers(ms%idx_age)%hTr(nx/2, ny/2, 2)/ms%h_layer(nx/2, ny/2, 2)
      expected = real(N_STEPS, wp)*DT
      call check(error, abs(conc - expected) < 1.0e-9_wp, &
                 "subsurface age after N steps must equal N·dt")
      call ms%destroy()
   end subroutine test_n_steps_linear

   subroutine test_surface_reset_isolates_top(error)
      !! After stepping (age kernel) then resetting (reset kernel,
      !! young_eff=0), the only layer touched by the reset is the
      !! surface; sub-surface layers keep their aged values.  This is
      !! the bottom-up-convention guard: MOM6's reference resets
      !! `k=1..nkbl` (its surface, top-down); Roundabout's `k=1` is the
      !! BED.  Resetting k=1 here would zero the wrong layer.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      real(wp) :: sub_age, surface_age
      integer :: nx, ny

      call setup(grid, ms, with_age=.true.)
      nx = grid%nx_total
      ny = grid%ny_total
      ms%tracers(ms%idx_age)%hTr = 0.0_wp

      ! Two age-steps so subsurface ages by 2·DT, then one reset.
      call ocean_ideal_age_age_step(ms%tracers(ms%idx_age)%hTr, ms%h_layer, &
                                    DT, nx, ny, NZ)
      call ocean_ideal_age_age_step(ms%tracers(ms%idx_age)%hTr, ms%h_layer, &
                                    DT, nx, ny, NZ)
      call ocean_ideal_age_reset_step(ms%tracers(ms%idx_age)%hTr, ms%h_layer, &
                                      0.0_wp, nx, ny, NZ)

      sub_age = ms%tracers(ms%idx_age)%hTr(nx/2, ny/2, 1)/H_EACH
      surface_age = ms%tracers(ms%idx_age)%hTr(nx/2, ny/2, NZ)/H_EACH
      call check(error, abs(sub_age - 2.0_wp*DT) < 1.0e-9_wp, &
                 "subsurface (bed, k=1) ages 2·dt over 2 steps and is untouched by the reset")
      if (.not. allocated(error)) call check(error, abs(surface_age) < 1.0e-12_wp, &
                                             "surface (k=nz) age is reset to zero")
      call ms%destroy()
   end subroutine test_surface_reset_isolates_top

   subroutine test_young_val_sets_surface(error)
      !! `ocean_ideal_age_reset_step` writes a CONCENTRATION times
      !! thickness (`young_eff * h_layer`), not a raw concentration —
      !! the single most likely units bug in the new kernel (the old
      !! hard-coded `young_eff=0` reset never exercised the h_layer
      !! factor since 0*h == 0 either way).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      real(wp), parameter :: YOUNG_EFF = 42.0_wp
      real(wp) :: expected_hTr, actual_hTr, sub_before
      integer :: nx, ny

      call setup(grid, ms, with_age=.true.)
      nx = grid%nx_total
      ny = grid%ny_total
      ms%tracers(ms%idx_age)%hTr = 0.0_wp
      call ocean_ideal_age_age_step(ms%tracers(ms%idx_age)%hTr, ms%h_layer, &
                                    DT, nx, ny, NZ)
      sub_before = ms%tracers(ms%idx_age)%hTr(nx/2, ny/2, 1)

      call ocean_ideal_age_reset_step(ms%tracers(ms%idx_age)%hTr, ms%h_layer, &
                                      YOUNG_EFF, nx, ny, NZ)

      expected_hTr = YOUNG_EFF*H_EACH
      actual_hTr = ms%tracers(ms%idx_age)%hTr(nx/2, ny/2, NZ)
      call check(error, abs(actual_hTr - expected_hTr) < 1.0e-12_wp, &
                 "surface hTr must equal young_eff * h_layer exactly")
      if (.not. allocated(error)) then
         call check(error, &
                    abs(ms%tracers(ms%idx_age)%hTr(nx/2, ny/2, 1) - sub_before) < 1.0e-12_wp, &
                    "reset must not touch subsurface layers")
      end if
      call ms%destroy()
   end subroutine test_young_val_sets_surface

   subroutine test_growth_rate_is_exponential(error)
      !! Analytical: with young_val = 1e-20 (MOM6's vintage-tracer
      !! seed) and rate = 1.057e-9 s^-1 (MOM6's 1/30 yr^-1 in
      !! seconds), the ratio of two evaluations must equal
      !! exp(rate*(t2-t1)) to 1e-12 relative.  A sign slip in the
      !! exponent would silently decay the vintage tracer instead of
      !! growing it, and no dynamical test would notice.
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: YOUNG_VAL = 1.0e-20_wp
      real(wp), parameter :: RATE = 1.057e-9_wp
      real(wp), parameter :: YEAR_S = 30.0_wp*365.0_wp*86400.0_wp
      real(wp) :: y1, y2, ratio, expected

      y1 = ocean_ideal_age_young_val(YOUNG_VAL, RATE, 0.0_wp)
      y2 = ocean_ideal_age_young_val(YOUNG_VAL, RATE, YEAR_S)
      ratio = y2/y1
      expected = exp(RATE*YEAR_S)
      call check(error, abs(ratio - expected)/expected < 1.0e-12_wp, &
                 "young_val ratio over 30 yr must equal exp(rate*dt)")
   end subroutine test_growth_rate_is_exponential

   subroutine test_default_knobs_short_circuit(error)
      !! At the shipped defaults (young_val=0, sfc_growth_rate=0) the
      !! evaluator must return bitwise 0, not merely close-to-zero —
      !! pinning that no exp() contamination reaches the default path
      !! (the bit-identity contract).
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: young

      young = ocean_ideal_age_young_val(0.0_wp, 0.0_wp, 1.0e9_wp)
      call check(error, young == 0.0_wp, &
                 "default knobs (0,0) must give bitwise-zero young_eff")
   end subroutine test_default_knobs_short_circuit

   subroutine test_age_survives_passive_registration(error)
      !! Registering a passive tracer via the PR-28 API AFTER ideal-age
      !! must leave idx_salinity / idx_temperature / idx_age unchanged
      !! and land the new tracer at the (new) last slot — the
      !! `register_passive_tracer` contract ("S/T/age keep their
      !! indices").
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      integer :: idx_s0, idx_t0, idx_age0, idx_new

      call setup(grid, ms, with_age=.true.)
      idx_s0 = ms%idx_salinity
      idx_t0 = ms%idx_temperature
      idx_age0 = ms%idx_age

      call ms%register_passive_tracer(grid, "dye", "1", "Dye tracer", idx_new)

      call check(error, idx_new == 4, "the new passive tracer should land at slot 4")
      if (.not. allocated(error)) call check(error, ms%idx_salinity == idx_s0, &
                                             "idx_salinity must be unchanged")
      if (.not. allocated(error)) call check(error, ms%idx_temperature == idx_t0, &
                                             "idx_temperature must be unchanged")
      if (.not. allocated(error)) call check(error, ms%idx_age == idx_age0, &
                                             "idx_age must be unchanged")
      if (.not. allocated(error)) call check(error, idx_new == size(ms%tracers), &
                                             "the new tracer must be the last registry slot")
      call ms%destroy()
   end subroutine test_age_survives_passive_registration

end module test_ocean_ideal_age
