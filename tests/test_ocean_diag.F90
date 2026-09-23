!! Phase 6a diag-manager tests.
!!
!! Exercises the diagnostic registry end-to-end (PR-1 scope):
!!   * `register` grows the registry, populates the var record.
!!   * `register` grows capacity via doubling past the initial cap.
!!   * `step` does NOT fire below cadence (dt_accum accumulates,
!!     output_buffer stays at the initial fill).
!!   * `step` fires at cadence, calls the bound `fill`, resets accum.
!!   * `register_default_diags` registers the six canonical vars
!!     (SSH, T, S, u, v, KE) on a fully-initialised ocean state.
!!   * SSH fill matches `h - b` on a stamped barotropic state.
!!   * KE fill matches `0.5 * (u² + v²)` on uniform flow.
module test_ocean_diag
   use rdb_constants, only: wp, REMAP_PCM, REMAP_PPM, H_VANISHED
   use rdb_grid, only: hgrid_t
   use rdb_ocean_state, only: ocean_state_t, ocean_state_enter_data, ocean_state_exit_data
   use rdb_ocean_diag, only: ocean_diag_t, DIAG_OP_INSTANT, DIAG_OP_MEAN, &
                             DIAG_OP_MAX, DIAG_OP_MIN, DIAG_OP_INTEGRAL, &
                             DIAG_VGRID_Z_FIXED, DIAG_OP_UNSET, &
                             DIAG_VGRID_SIGMA, DIAG_VGRID_LAYER, DIAG_COORD_UNSET, &
                             DIAG_VGRID_DENSITY, &
                             diag_spec_t, parse_diag_spec
   use rdb_ocean_diag_fills, only: register_default_diags, fill_ssh, fill_ke, &
                                   fill_temperature, remap_layer_to_z, &
                                   set_diag_remap_method, is_canonical_diag_name, &
                                   set_diag_mask_vanished, canonical_diag_gate_hint
   use rdb_ocean_diag_derived, only: register_derived, apply_diag_selection, &
                                     derived_catalog_size, derived_catalog_name
   use rdb_ocean_diag_mask, only: diag_mask_t, diag_mask_global, &
                                  diag_mask_bbox, diag_mask_h_section, &
                                  diag_mask_v_section
   use rdb_ocean_metrics, only: metrics_fill_cartesian, metrics_finalize
   use rdb_ocean_budgets, only: ocean_budgets_t, BUDGET_MASS, BUDGET_SALT_TOTAL, &
                                BUDGET_HEAT_TOTAL, BUDGET_KE, &
                                budget_total_mass, budget_total_tracer, &
                                budget_total_ke, budget_contributor_t
   use testdrive, only: error_type, check, new_unittest, unittest_type
   implicit none
   private

   public :: collect_ocean_diag_tests

   integer, parameter :: NX = 6, NY = 4, NZ = 3
   real(wp), parameter :: DX = 1.0_wp

   ! Module-level counter for the stateful `counting_fill` test helper.
   ! Reset to zero at the top of each test that uses it; the fill
   ! increments and stamps the buffer with the new call count.
   integer :: counting_fill_count = 0

contains

   subroutine collect_ocean_diag_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("register_grows_registry", test_register_grows), &
                  new_unittest("register_doubles_capacity", test_register_grows_capacity), &
                  new_unittest("step_does_not_fire_below_cadence", test_step_below_cadence), &
                  new_unittest("step_fires_at_cadence", test_step_fires), &
                  new_unittest("register_default_diags_six_vars", test_default_diag_set), &
                  new_unittest("diags_spec_adds_derived", test_extra_diags), &
                  new_unittest("diags_spec_empty_noop", test_extra_diags_empty), &
                  new_unittest("diags_spec_off_skips_canonical", test_disable_diags), &
                  new_unittest("diags_spec_blank_keeps_canonical", test_disable_diags_empty), &
                  new_unittest("disabled_diag_skipped_in_step", test_disabled_skipped), &
                  new_unittest("parse_diag_spec_attributes", test_parse_spec), &
                  new_unittest("parse_diag_spec_coord", test_parse_spec_coord), &
                  new_unittest("spec_parses_integral_op", test_spec_parses_integral_op), &
                  new_unittest("spec_parses_density_coord", test_spec_parses_density_coord), &
                  new_unittest("diags_spec_cadence_op_override", test_spec_override), &
                  new_unittest("diags_spec_coord_routes_to_vgrid", test_spec_coord_routing), &
                  new_unittest("mask_vanished_sets_has_missing", test_mask_has_missing), &
                  new_unittest("is_canonical_diag_name_recognises", test_is_canonical), &
                  new_unittest("register_default_diags_adds_age_when_enabled", test_default_diag_age), &
                  new_unittest("fill_age_matches_htr_over_h", test_fill_age), &
                  new_unittest("fill_temperature_vanished_layer_is_nan", test_fill_vanished_nan), &
                  new_unittest("fill_ssh_matches_h_minus_b", test_fill_ssh), &
                  new_unittest("fill_ke_matches_half_u2_plus_v2", test_fill_ke), &
                  new_unittest("set_output_z_levels_stores_grid", test_set_z_levels), &
                  new_unittest("zfixed_register_allocates_layer_and_output", test_zfixed_alloc), &
                  new_unittest("remap_identity_at_layer_centres", test_remap_identity), &
                  new_unittest("remap_clamps_above_surface", test_remap_clamp_top), &
                  new_unittest("remap_clamps_below_bed", test_remap_clamp_bot), &
                  new_unittest("remap_step_end_to_end_temperature", test_remap_e2e_temp), &
                  new_unittest("mean_fills_every_step_and_averages", test_mean_basic), &
                  new_unittest("mean_resets_after_fire", test_mean_reset), &
                  new_unittest("max_tracks_max_sample", test_max), &
                  new_unittest("min_tracks_min_sample", test_min), &
                  new_unittest("instant_does_not_accumulate", test_instant_no_accum), &
                  new_unittest("mean_dt_weighted_varies_with_dt", test_mean_dt_weighted), &
                  new_unittest("integral_constant_rate", test_integral_constant_rate), &
                  new_unittest("derived_catalog_lists_every_entry", test_derived_catalog_size), &
                  new_unittest("derived_h_layer_direct_copy", test_derived_h_layer), &
                  new_unittest("derived_rho_layer_direct_copy", test_derived_rho_layer), &
                  new_unittest("derived_vorticity_z_rigid_rotation", test_derived_vorticity), &
                  new_unittest("derived_ke_total_uniform_flow", test_derived_ke_total), &
                  new_unittest("derived_transport_x_uniform_flow", test_derived_transport_x), &
                  new_unittest("derived_transport_y_uniform_flow", test_derived_transport_y), &
                  new_unittest("derived_mld_density_step_profile", test_derived_mld), &
                  new_unittest("derived_mld_density_ignores_vanished_filler", &
                               test_derived_mld_vanished), &
                  new_unittest("mask_global_covers_full_grid", test_mask_global), &
                  new_unittest("mask_bbox_covers_only_box", test_mask_bbox), &
                  new_unittest("mask_h_section_strip", test_mask_h_section), &
                  new_unittest("mask_v_section_strip", test_mask_v_section), &
                  new_unittest("mask_zeros_outside_in_fold", test_mask_zeros_outside), &
                  new_unittest("budgets_evaluate_total_mass", test_budget_total_mass), &
                  new_unittest("budgets_evaluate_total_tracer", test_budget_total_tracer), &
                  new_unittest("budgets_evaluate_total_ke", test_budget_total_ke), &
                  new_unittest("budgets_init_snapshot_records_values", test_budget_snapshot), &
                  new_unittest("budgets_step_emits_no_drift_on_rest", test_budget_no_drift), &
                  new_unittest("budgets_register_contributor_appends", test_budget_register_contributor), &
                  new_unittest("budgets_drain_integrates_then_zeroes", test_budget_drain_lifecycle), &
                  new_unittest("budgets_drain_accumulates_across_calls", test_budget_drain_accumulates), &
                  new_unittest("budgets_register_same_name_is_idempotent", test_budget_register_idempotent), &
                  new_unittest("diag_gated_canonical_is_not_registered", test_gated_canonical_not_registered), &
                  new_unittest("diag_gate_open_registers_and_does_not_warn", test_gate_open_registers), &
                  new_unittest("diag_off_on_gated_canonical_stays_silent", test_off_on_gated_stays_silent), &
                  new_unittest("diag_typo_still_routes_to_derived_not_warn", test_typo_routes_to_derived), &
                  new_unittest("diag_gate_hint_covers_every_gate", test_gate_hint_coverage) &
                  ]
   end subroutine collect_ocean_diag_tests

   subroutine dummy_fill(state_handle, buf)
      !! Trivial fill — writes a known constant into the buffer.  Used
      !! by the cadence / registry tests; doesn't touch the state.
      class(*), intent(in) :: state_handle
      real(wp), intent(inout) :: buf(:, :, :)
      buf = 42.0_wp
      if (.false.) then
         select type (state_handle)
         class default
         end select
      end if
   end subroutine dummy_fill

   subroutine counting_fill(state_handle, buf)
      !! Stateful fill — increments `counting_fill_count` on every
      !! call and stamps the buffer with the new count.  Used by the
      !! time-op tests to drive a known sequence of samples.
      class(*), intent(in) :: state_handle
      real(wp), intent(inout) :: buf(:, :, :)
      counting_fill_count = counting_fill_count + 1
      buf = real(counting_fill_count, wp)
      if (.false.) then
         select type (state_handle)
         class default
         end select
      end if
   end subroutine counting_fill

   subroutine setup_state(grid, state)
      type(hgrid_t), intent(inout) :: grid
      type(ocean_state_t), intent(inout) :: state
      call grid%init(NX, NY, 1, DX, DX)
      state%multilayer%nz_ml = NZ
      call state%init(grid)
      ! Fill the metrics slot (cartesian) so derived fills + console/
      ! budget integrals read real per-cell metrics rather than the
      ! zero-initialised arrays.  Mirrors configure_ocean_metrics.
      call metrics_fill_cartesian(state%metrics, grid, grid%dx, grid%dy)
      call metrics_finalize(state%metrics)
   end subroutine setup_state

   subroutine setup_state_age(grid, state)
      !! Same as setup_state but with the ideal-age tracer enabled, so the
      !! multilayer registry carries a third tracer (idx_age) and
      !! register_default_diags exposes an "age" field.
      type(hgrid_t), intent(inout) :: grid
      type(ocean_state_t), intent(inout) :: state
      call grid%init(NX, NY, 1, DX, DX)
      state%multilayer%nz_ml = NZ
      state%enable_ideal_age = .true.
      call state%init(grid)
      call metrics_fill_cartesian(state%metrics, grid, grid%dx, grid%dy)
      call metrics_finalize(state%metrics)
   end subroutine setup_state_age

   subroutine setup_state_ice(grid, state, ncat)
      !! setup_state + the sea-ice slot enabled (ncat categories), so
      !! `state%init` allocates the ice arrays and `register_default_diags`
      !! exposes ice_conc/ice_thick — used as the "gate open" control for
      !! the PR-64 diag-selection-warn tests.
      type(hgrid_t), intent(inout) :: grid
      type(ocean_state_t), intent(inout) :: state
      integer, intent(in) :: ncat
      call grid%init(NX, NY, 1, DX, DX)
      state%multilayer%nz_ml = NZ
      state%ice%enable = .true.
      state%ice%ncat = ncat        ! before state%init — sizes the arrays
      call state%init(grid)
      call metrics_fill_cartesian(state%metrics, grid, grid%dx, grid%dy)
      call metrics_finalize(state%metrics)
   end subroutine setup_state_ice

   subroutine test_register_grows(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      checks: block
         call setup_state(grid, state)

         call state%diag%register("foo", units="m", fill=dummy_fill, n1=NX, n2=NY, n3=1)
         call check(error, state%diag%nvars == 1, "nvars should be 1 after one register")
         if (allocated(error)) exit checks
         call check(error, state%diag%vars(1)%name == "foo", "name not stored")
         if (allocated(error)) exit checks
         call check(error, associated(state%diag%vars(1)%fill), "fill not bound")
         if (allocated(error)) exit checks
         call check(error, allocated(state%diag%vars(1)%output_buffer), &
                    "output_buffer not allocated")

      end block checks
      call state%destroy()
   end subroutine test_register_grows

   subroutine test_register_grows_capacity(error)
      !! Initial capacity is 16; register 20 to force a doubling.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      integer :: i
      character(len=8) :: name
      checks: block

         call setup_state(grid, state)
         do i = 1, 20
            write (name, "(A,I0)") "v", i
            call state%diag%register(trim(name), units="-", fill=dummy_fill, &
                                     n1=NX, n2=NY, n3=1)
         end do

         call check(error, state%diag%nvars == 20, "nvars should be 20")
         if (allocated(error)) exit checks
         call check(error, state%diag%nvars_max >= 20, "capacity should have grown")
         if (allocated(error)) exit checks
         call check(error, state%diag%vars(20)%name == "v20", "post-grow var name preserved")

      end block checks
      call state%destroy()
   end subroutine test_register_grows_capacity

   subroutine test_step_below_cadence(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      checks: block

         call setup_state(grid, state)
         call state%diag%register("foo", units="m", fill=dummy_fill, &
                                  n1=NX, n2=NY, n3=1, &
                                  time_op=DIAG_OP_INSTANT, dt_out=10.0_wp)

         call state%diag%step(state, dt=4.0_wp, t=4.0_wp)
         call check(error, abs(state%diag%vars(1)%dt_accum - 4.0_wp) < 1.0e-12_wp, &
                    "dt_accum should accumulate but not fire")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(state%diag%vars(1)%output_buffer)) < 1.0e-12_wp, &
                    "fill should NOT have run yet")
         if (allocated(error)) exit checks
         call check(error, state%diag%vars(1)%n_accum == 0, "n_accum should be 0")

      end block checks
      call state%destroy()
   end subroutine test_step_below_cadence

   subroutine test_step_fires(error)
      !! dt_out=10, step with dt=12 → fire, buffer = 42, dt_accum reset.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      checks: block

         call setup_state(grid, state)
         call state%diag%register("foo", units="m", fill=dummy_fill, &
                                  n1=NX, n2=NY, n3=1, &
                                  time_op=DIAG_OP_INSTANT, dt_out=10.0_wp)

         call state%diag%step(state, dt=12.0_wp, t=12.0_wp)
         call check(error, abs(state%diag%vars(1)%dt_accum) < 1.0e-12_wp, &
                    "dt_accum should reset on fire")
         if (allocated(error)) exit checks
         call check(error, abs(state%diag%vars(1)%output_buffer(2, 2, 1) - 42.0_wp) < 1.0e-12_wp, &
                    "dummy_fill should have stamped buffer with 42")

      end block checks
      call state%destroy()
   end subroutine test_step_fires

   subroutine test_default_diag_set(error)
      !! `register_default_diags` registers SSH + 2 tracers + u + v + KE.
      !! Default multilayer init carries idx_salinity + idx_temperature.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      checks: block

         call setup_state(grid, state)
         call register_default_diags(state, dt_out=3600.0_wp)

         call check(error, state%diag%nvars == 6, "should register 6 default vars")
         if (allocated(error)) exit checks
         call check(error, state%diag%vars(1)%name == "SSH", "first var should be SSH")
         if (allocated(error)) exit checks
         call check(error, state%diag%vars(6)%name == "KE", "last var should be KE")

      end block checks
      call state%destroy()
   end subroutine test_default_diag_set

   subroutine test_extra_diags(error)
      !! `apply_diag_selection` adds named derived-catalog diags on top of the
      !! canonical set; comma + space separators both parse; INSTANT time-op.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      integer :: n0
      checks: block
         call setup_state(grid, state)
         call register_default_diags(state, dt_out=3600.0_wp)
         n0 = state%diag%nvars
         call state%diag%destroy()
         call state%diag%init(grid)
         call apply_diag_selection(state, "h_layer, rho_layer", dt_out=3600.0_wp)
         call check(error, state%diag%nvars == n0 + 2, "should add 2 derived vars")
         if (allocated(error)) exit checks
         call check(error, state%diag%vars(n0 + 1)%name == "h_layer", "first derived = h_layer")
         if (allocated(error)) exit checks
         call check(error, state%diag%vars(n0 + 2)%name == "rho_layer", "second derived = rho_layer")
         if (allocated(error)) exit checks
         call check(error, state%diag%vars(n0 + 1)%time_op == DIAG_OP_INSTANT, &
                    "derived diags register as INSTANT by default")
      end block checks
      call state%destroy()
   end subroutine test_extra_diags

   subroutine test_extra_diags_empty(error)
      !! Empty / blank `diags` reproduces the canonical set (bit-identical).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      integer :: n0
      checks: block
         call setup_state(grid, state)
         call register_default_diags(state, dt_out=3600.0_wp)
         n0 = state%diag%nvars
         call state%diag%destroy()
         call state%diag%init(grid)
         call apply_diag_selection(state, "   ", dt_out=3600.0_wp)
         call check(error, state%diag%nvars == n0, "blank diags == canonical set")
      end block checks
      call state%destroy()
   end subroutine test_extra_diags_empty

   subroutine test_disable_diags(error)
      !! A `:off` entry skips registration of the named canonical diag
      !! entirely (it is absent from the registry), leaving the rest.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      integer :: i, n_full
      logical :: ke_present
      checks: block
         call setup_state(grid, state)
         call register_default_diags(state, dt_out=3600.0_wp)
         n_full = state%diag%nvars
         call state%diag%destroy()
         call state%diag%init(grid)
         call apply_diag_selection(state, "KE:off", dt_out=3600.0_wp)
         call check(error, state%diag%nvars == n_full - 1, "KE:off drops one var")
         if (allocated(error)) exit checks
         ke_present = .false.
         do i = 1, state%diag%nvars
            if (state%diag%vars(i)%name == "KE") ke_present = .true.
         end do
         call check(error,.not. ke_present, "KE must be absent after :off")
      end block checks
      call state%destroy()
   end subroutine test_disable_diags

   subroutine test_disable_diags_empty(error)
      !! Blank `diags` keeps the full canonical set (the bit-identical default).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      integer :: i
      logical :: all_on
      checks: block
         call setup_state(grid, state)
         call apply_diag_selection(state, "   ", dt_out=3600.0_wp)
         all_on = .true.
         do i = 1, state%diag%nvars
            all_on = all_on .and. state%diag%vars(i)%enabled
         end do
         call check(error, all_on, "blank diags must leave all enabled")
      end block checks
      call state%destroy()
   end subroutine test_disable_diags_empty

   subroutine test_parse_spec(error)
      !! `parse_diag_spec` classifies self-identifying colon attributes
      !! (name / off / op / cadence) order-free and case-insensitively.
      type(error_type), allocatable, intent(out) :: error
      type(diag_spec_t), allocatable :: s(:)
      checks: block
         s = parse_diag_spec("KE:off  temperature:6h:MEAN  vorticity_z")
         call check(error, size(s) == 3, "three entries parsed")
         if (allocated(error)) exit checks
         call check(error, trim(s(1)%name) == "KE" .and. s(1)%off, "KE:off")
         if (allocated(error)) exit checks
         call check(error, trim(s(2)%name) == "temperature", "name = temperature")
         if (allocated(error)) exit checks
         call check(error, abs(s(2)%dt_out - 21600.0_wp) < 1.0e-9_wp, "6h => 21600 s")
         if (allocated(error)) exit checks
         call check(error, s(2)%time_op == DIAG_OP_MEAN, "MEAN op (case-insensitive)")
         if (allocated(error)) exit checks
         call check(error, trim(s(3)%name) == "vorticity_z" .and. .not. s(3)%off &
                    .and. s(3)%time_op == DIAG_OP_UNSET .and. s(3)%dt_out < 0.0_wp, &
                    "bare name => no attributes, sentinels intact")
         if (allocated(error)) exit checks
         call check(error, size(parse_diag_spec("   ")) == 0, "blank => empty array")
      end block checks
   end subroutine test_parse_spec

   subroutine test_spec_override(error)
      !! A canonical-name spec entry overrides its cadence + time-op at
      !! registration (no buffer realloc), without dropping the var.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      integer :: i, ke_idx
      checks: block
         call setup_state(grid, state)
         call apply_diag_selection(state, "KE:1h:max", dt_out=3600.0_wp)
         ke_idx = 0
         do i = 1, state%diag%nvars
            if (state%diag%vars(i)%name == "KE") ke_idx = i
         end do
         call check(error, ke_idx > 0, "KE still registered")
         if (allocated(error)) exit checks
         call check(error, state%diag%vars(ke_idx)%time_op == DIAG_OP_MAX, &
                    "KE op overridden to MAX")
         if (allocated(error)) exit checks
         call check(error, abs(state%diag%vars(ke_idx)%dt_out - 3600.0_wp) < 1.0e-9_wp, &
                    "KE cadence overridden to 1h")
      end block checks
      call state%destroy()
   end subroutine test_spec_override

   subroutine test_parse_spec_coord(error)
      !! `parse_diag_spec` classifies the `:coord` attribute (layer/z/zstar/
      !! sigma) and leaves it UNSET when absent.
      type(error_type), allocatable, intent(out) :: error
      type(diag_spec_t), allocatable :: s(:)
      checks: block
         s = parse_diag_spec("temperature:sigma  KE  vorticity_z:z:1d")
         call check(error, s(1)%coord == DIAG_VGRID_SIGMA, "temperature:sigma -> SIGMA")
         if (allocated(error)) exit checks
         call check(error, s(2)%coord == DIAG_COORD_UNSET, "bare KE -> coord unset")
         if (allocated(error)) exit checks
         call check(error, s(3)%coord == DIAG_VGRID_Z_FIXED, "vorticity_z:z -> Z_FIXED")
         if (allocated(error)) exit checks
         call check(error, abs(s(3)%dt_out - 86400.0_wp) < 1.0e-9_wp, &
                    "coord + cadence parse together (1d)")
      end block checks
   end subroutine test_parse_spec_coord

   subroutine test_spec_parses_integral_op(error)
      !! PR-9 §9.3 wire test: `parse_diag_spec` recognises the `:integral`
      !! attribute (`rdb_ocean_diag.F90:parse_one_spec_token`) and maps it
      !! to `DIAG_OP_INTEGRAL`.  `test_integral_constant_rate` already
      !! proves the OP's numerics end-to-end via a direct `register(...
      !! time_op=DIAG_OP_INTEGRAL)` call; this proves the STRING reaches
      !! that same enum — before this PR `:integral` in a `diags` string
      !! was a hard `error stop`, not a selection.
      type(error_type), allocatable, intent(out) :: error
      type(diag_spec_t), allocatable :: s(:)
      checks: block
         s = parse_diag_spec("temperature:integral")
         call check(error, size(s) == 1, "one entry parsed")
         if (allocated(error)) exit checks
         call check(error, trim(s(1)%name) == "temperature", "name = temperature")
         if (allocated(error)) exit checks
         call check(error, s(1)%time_op == DIAG_OP_INTEGRAL, &
                    "':integral' must parse to DIAG_OP_INTEGRAL")
      end block checks
   end subroutine test_spec_parses_integral_op

   subroutine test_spec_parses_density_coord(error)
      !! PR-9 §9.4 (parser half): `parse_diag_spec` recognises the
      !! `:density`/`:rho` attribute and maps it to `DIAG_VGRID_DENSITY`.
      !! Before this PR neither string was recognised by
      !! `parse_one_spec_token` at all.
      type(error_type), allocatable, intent(out) :: error
      type(diag_spec_t), allocatable :: s(:)
      checks: block
         s = parse_diag_spec("temperature:density  salinity:rho")
         call check(error, size(s) == 2, "two entries parsed")
         if (allocated(error)) exit checks
         call check(error, s(1)%coord == DIAG_VGRID_DENSITY, &
                    "':density' must parse to DIAG_VGRID_DENSITY")
         if (allocated(error)) exit checks
         call check(error, s(2)%coord == DIAG_VGRID_DENSITY, &
                    "':rho' must parse to DIAG_VGRID_DENSITY")
      end block checks
   end subroutine test_spec_parses_density_coord

   subroutine test_spec_coord_routing(error)
      !! A `:sigma` coord on a layered canonical diag routes it to the SIGMA
      !! output vgrid with a remap attached; a 2D diag (SSH) ignores coord.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      integer :: i, t_idx, ssh_idx
      checks: block
         call setup_state(grid, state)
         call state%diag%set_output_sigma_levels([0.5_wp, 1.0_wp])
         call apply_diag_selection(state, "temperature:sigma  SSH:sigma", dt_out=3600.0_wp)
         t_idx = 0
         ssh_idx = 0
         do i = 1, state%diag%nvars
            if (state%diag%vars(i)%name == "temperature") t_idx = i
            if (state%diag%vars(i)%name == "SSH") ssh_idx = i
         end do
         call check(error, t_idx > 0 .and. ssh_idx > 0, "both diags registered")
         if (allocated(error)) exit checks
         call check(error, state%diag%vars(t_idx)%output_vgrid == DIAG_VGRID_SIGMA, &
                    "temperature routed to SIGMA vgrid")
         if (allocated(error)) exit checks
         call check(error, associated(state%diag%vars(t_idx)%remap), &
                    "temperature got a remap proc")
         if (allocated(error)) exit checks
         call check(error, state%diag%vars(ssh_idx)%output_vgrid == DIAG_VGRID_LAYER, &
                    "2D SSH ignores coord (stays LAYER)")
      end block checks
      call state%destroy()
   end subroutine test_spec_coord_routing

   subroutine test_mask_has_missing(error)
      !! With vanished-masking on, a remapped (non-layer) diag is tagged
      !! `has_missing` (so the writer emits _FillValue); a LAYER diag is not.
      !! Resets the module flag so it does not leak into other tests.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      integer :: i, t_idx, ssh_idx
      checks: block
         call setup_state(grid, state)
         call state%diag%set_output_sigma_levels([0.5_wp, 1.0_wp])
         call set_diag_mask_vanished(.true.)
         call apply_diag_selection(state, "temperature:sigma", dt_out=3600.0_wp)
         t_idx = 0
         ssh_idx = 0
         do i = 1, state%diag%nvars
            if (state%diag%vars(i)%name == "temperature") t_idx = i
            if (state%diag%vars(i)%name == "SSH") ssh_idx = i
         end do
         call check(error, t_idx > 0 .and. state%diag%vars(t_idx)%has_missing, &
                    "remapped temperature tagged has_missing under masking")
         if (allocated(error)) exit checks
         call check(error, ssh_idx > 0 .and. .not. state%diag%vars(ssh_idx)%has_missing, &
                    "LAYER SSH not tagged has_missing")
      end block checks
      call set_diag_mask_vanished(.false.)   ! reset module state
      call state%destroy()
   end subroutine test_mask_has_missing

   subroutine test_is_canonical(error)
      !! `is_canonical_diag_name` recognises canonical defaults and rejects
      !! derived-catalog / unknown names (routes spec entries correctly).
      type(error_type), allocatable, intent(out) :: error
      checks: block
         call check(error, is_canonical_diag_name("SSH"), "SSH is canonical")
         if (allocated(error)) exit checks
         call check(error, is_canonical_diag_name("KE"), "KE is canonical")
         if (allocated(error)) exit checks
         call check(error,.not. is_canonical_diag_name("vorticity_z"), &
                    "vorticity_z is derived, not canonical")
         if (allocated(error)) exit checks
         call check(error,.not. is_canonical_diag_name("nonsense"), &
                    "unknown name is not canonical")
      end block checks
   end subroutine test_is_canonical

   subroutine test_gated_canonical_not_registered(error)
      !! PR-64 §9 case 1 — the gap, stated as a fact. On a state with
      !! ice/EPBL off and no age tracer (setup_state's defaults),
      !! `apply_diag_selection` requesting "age ice_conc Kd_EPBL" leaves all
      !! three UNREGISTERED — the silent drop this PR adds a warning for.
      !! nvars stays at the unconditional-four-plus-two-tracers set (SSH,
      !! temperature, salinity, u, v, KE). This case passes before AND
      !! after the PR — it pins the behaviour the warn describes, so a
      !! later change cannot "fix" the silence by registering a
      !! zero-filled placeholder and calling it done.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      checks: block
         call setup_state(grid, state)
         call apply_diag_selection(state, "age ice_conc Kd_EPBL", dt_out=3600.0_wp)

         call check(error,.not. state%diag%is_registered("age"), &
                    "age should NOT be registered: idx_age == 0")
         if (allocated(error)) exit checks
         call check(error,.not. state%diag%is_registered("ice_conc"), &
                    "ice_conc should NOT be registered: ice%enable == .false.")
         if (allocated(error)) exit checks
         call check(error,.not. state%diag%is_registered("Kd_EPBL"), &
                    "Kd_EPBL should NOT be registered: epbl%enable == .false.")
         if (allocated(error)) exit checks
         call check(error, state%diag%nvars == 6, &
                    "only the unconditional set (SSH,T,S,u,v,KE) should register")
      end block checks
      call state%destroy()
   end subroutine test_gated_canonical_not_registered

   subroutine test_gate_open_registers(error)
      !! PR-64 §9 case 2 — the control. Same spec as case 1, but with the
      !! ice gate OPEN: ice_conc must register. Without this control, case
      !! 1 would pass equally on a build where `apply_diag_selection`
      !! registers nothing at all.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      checks: block
         call setup_state_ice(grid, state, ncat=1)
         call apply_diag_selection(state, "age ice_conc Kd_EPBL", dt_out=3600.0_wp)

         call check(error, state%diag%is_registered("ice_conc"), &
                    "ice_conc SHOULD register when ice%enable == .true.")
      end block checks
      call state%destroy()
   end subroutine test_gate_open_registers

   subroutine test_off_on_gated_stays_silent(error)
      !! PR-64 §9 case 3 — the discriminator. `"ice_conc:off"` with the ice
      !! gate CLOSED must stay silent: ice_conc is unregistered (same
      !! observable as case 1), reached via the pre-existing `:off` cycle
      !! at apply_diag_selection's :234, not the new warn branch. Confirms
      !! the `:off` flag round-trips through `parse_diag_spec` so a
      !! reviewer can see the guard is hit BEFORE the new check. This is
      !! the case that fails if the warn is "simplified" to fire before
      !! the `off` guard, which would spray a warning on every reusable
      !! namelist that defensively disables ice diags.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      type(diag_spec_t), allocatable :: specs(:)
      checks: block
         specs = parse_diag_spec("ice_conc:off")
         call check(error, size(specs) == 1 .and. trim(specs(1)%name) == "ice_conc", &
                    "parse_diag_spec should produce one 'ice_conc' entry")
         if (allocated(error)) exit checks
         call check(error, specs(1)%off, &
                    "parse_diag_spec should set off=.true. for ':off'")
         if (allocated(error)) exit checks

         call setup_state(grid, state)
         call apply_diag_selection(state, "ice_conc:off", dt_out=3600.0_wp)
         call check(error,.not. state%diag%is_registered("ice_conc"), &
                    "ice_conc:off on a gated-off ice run must stay silently unregistered")
      end block checks
      call state%destroy()
   end subroutine test_off_on_gated_stays_silent

   subroutine test_typo_routes_to_derived(error)
      !! PR-64 §9 case 4 — the regression fence. A plausible typo of `age`
      !! ("aeg") is NOT a canonical name, so `apply_diag_selection` must
      !! route it past the new "canonical but ungated" warn branch to the
      !! derived-catalog path (`register_derived`, unchanged by this PR),
      !! which still fails loud on an unknown name. `register_derived`'s
      !! `error stop` cannot be exercised in-process (test-drive would
      !! terminate with the whole ctest binary — see
      !! test_ocean_periodic::test_validation_rejections for the same
      !! convention); we instead pin the ROUTING decision that a
      !! well-meaning implementer could accidentally soften: a typo must
      !! never be classified as canonical.
      type(error_type), allocatable, intent(out) :: error
      checks: block
         call check(error,.not. is_canonical_diag_name("aeg"), &
                    "'aeg' (a plausible 'age' typo) must not be canonical, so it " &
                    //"still routes to register_derived's fail-loud unknown-name path")
      end block checks
   end subroutine test_typo_routes_to_derived

   subroutine test_gate_hint_coverage(error)
      !! PR-64 §9 case 5 — the anti-drift lock. Each of the eight
      !! conditionally-registered canonical names has a non-empty
      !! `canonical_diag_gate_hint`; each of the four unconditional names
      !! (SSH/u/v/KE) has an empty one. The count of non-empty hints across
      !! the full canonical set must equal 8 — so a future PR that adds a
      !! ninth gated canonical diagnostic to `register_default_diags`
      !! without adding its hint here fails THIS test (§13's lock-step
      !! contract).
      type(error_type), allocatable, intent(out) :: error
      character(len=16), parameter :: gated(*) = [character(len=16) :: &
                                                  "temperature", "salinity", "age", &
                                                  "MLD_EPBL", "Kd_EPBL", "Kd_KSHEAR", &
                                                  "ice_conc", "ice_thick"]
      character(len=16), parameter :: ungated(*) = [character(len=16) :: &
                                                    "SSH", "u", "v", "KE"]
      character(len=16), parameter :: all_canonical(*) = [character(len=16) :: &
                                                          "SSH", "temperature", "salinity", &
                                                          "age", "u", "v", "KE", &
                                                          "MLD_EPBL", "Kd_EPBL", "Kd_KSHEAR", &
                                                          "ice_conc", "ice_thick"]
      integer :: i, n_hinted
      checks: block
         do i = 1, size(gated)
            call check(error, len_trim(canonical_diag_gate_hint(trim(gated(i)))) > 0, &
                       "gated canonical '"//trim(gated(i))//"' must have a non-empty hint")
            if (allocated(error)) exit checks
         end do
         do i = 1, size(ungated)
            call check(error, len_trim(canonical_diag_gate_hint(trim(ungated(i)))) == 0, &
                       "unconditional canonical '"//trim(ungated(i))//"' must have an empty hint")
            if (allocated(error)) exit checks
         end do

         n_hinted = 0
         do i = 1, size(all_canonical)
            if (len_trim(canonical_diag_gate_hint(trim(all_canonical(i)))) > 0) then
               n_hinted = n_hinted + 1
            end if
         end do
         call check(error, n_hinted == 8, &
                    "exactly 8 of the 12 canonical names should carry a gate hint")
      end block checks
   end subroutine test_gate_hint_coverage

   subroutine test_disabled_skipped(error)
      !! A disabled var's `fill` is skipped by `step`: its buffer stays at the
      !! init value (0), unlike an enabled var which `dummy_fill` stamps 42.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      checks: block
         call setup_state(grid, state)
         call state%diag%register("foo", units="m", fill=dummy_fill, &
                                  n1=NX, n2=NY, n3=1, &
                                  time_op=DIAG_OP_INSTANT, dt_out=10.0_wp)
         call state%diag%disable("foo")
         call state%diag%step(state, dt=12.0_wp, t=12.0_wp)
         call check(error, abs(state%diag%vars(1)%output_buffer(2, 2, 1)) < 1.0e-12_wp, &
                    "disabled var must NOT be filled (buffer stays 0, not 42)")
      end block checks
      call state%destroy()
   end subroutine test_disabled_skipped

   subroutine test_default_diag_age(error)
      !! With the ideal-age tracer enabled, register_default_diags exposes
      !! a 7th variable "age" (units s).  Guards the wiring that makes the
      !! tracer observable in output; default-off keeps the 6-var set
      !! (test_default_diag_set) bit-identical.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      integer :: iv, iage
      checks: block

         call setup_state_age(grid, state)
         call register_default_diags(state, dt_out=3600.0_wp)

         call check(error, state%diag%nvars == 7, "age-enabled set should register 7 vars")
         if (allocated(error)) exit checks
         iage = 0
         do iv = 1, state%diag%nvars
            if (state%diag%vars(iv)%name == "age") iage = iv
         end do
         call check(error, iage > 0, "an 'age' diagnostic should be registered")
         if (allocated(error)) exit checks
         call check(error, trim(state%diag%vars(iage)%units) == "s", "age units should be s")

      end block checks
      call state%destroy()
   end subroutine test_default_diag_age

   subroutine test_fill_age(error)
      !! End-to-end: seed hTr_age + h_layer, fire the registered age diag
      !! via step, and confirm the buffer reads back the concentration
      !! age = hTr_age / h_layer (the fill runs on-device under enter_data).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      integer :: iv, iage, i, j, k
      real(wp), parameter :: H_COL = 50.0_wp, AGE_S = 1234.0_wp
      checks: block

         call setup_state_age(grid, state)
         do k = 1, NZ
            do j = 1, NY
               do i = 1, NX
                  state%multilayer%h_layer(i, j, k) = H_COL
                  state%multilayer%tracers(state%multilayer%idx_age)%hTr(i, j, k) = H_COL*AGE_S
               end do
            end do
         end do

         call register_default_diags(state, dt_out=1.0_wp)
         iage = 0
         do iv = 1, state%diag%nvars
            if (state%diag%vars(iv)%name == "age") iage = iv
         end do
         call check(error, iage > 0, "age diag must be registered")
         if (allocated(error)) exit checks

         call ocean_state_enter_data(state)
         call state%diag%step(state, dt=2.0_wp, t=2.0_wp)
         call ocean_state_exit_data(state)

         call check(error, abs(state%diag%vars(iage)%output_buffer(3, 2, 2) - AGE_S) < 1.0e-9_wp, &
                    "age readout should equal hTr_age / h_layer")

      end block checks
      call state%destroy()
   end subroutine test_fill_age

   subroutine test_fill_vanished_nan(error)
      !! P7 F5 regression: fill_tracer_impl's vanished-layer sentinel
      !! must be NaN, not 0 -- 0 degC is a LEGAL ocean value, so a
      !! vanished bed layer used to read back as plausible ice-point
      !! freshwater rather than as missing data
      !! (rdb_ocean_diag_fills.F90's fill_tracer_impl docstring).
      !! Vanishes ONE bed-layer (k=1) cell's h_layer/hTr below
      !! H_VANISHED and confirms fill_temperature reads NaN there while
      !! an untouched neighbour cell still reads the ordinary hTr/h.
      use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      real(wp) :: buf(NX, NY, NZ)
      integer :: it
      real(wp), parameter :: H_COL = 10.0_wp, T_REF = 12.0_wp
      checks: block
         call setup_state(grid, state)
         it = state%multilayer%idx_temperature
         call check(error, it > 0, "temperature tracer must be registered")
         if (allocated(error)) exit checks

         state%multilayer%h_layer = H_COL
         state%multilayer%tracers(it)%hTr = H_COL*T_REF
         ! Vanish one bed-layer cell (k=1): both h and hTr go to zero,
         ! matching what a real dynamic-vanish (e.g. ZSTAR_FULL) leaves
         ! behind -- an exact-zero thickness, not merely "small".
         state%multilayer%h_layer(3, 2, 1) = 0.0_wp
         state%multilayer%tracers(it)%hTr(3, 2, 1) = 0.0_wp

         call fill_temperature(state, buf)

         call check(error, ieee_is_nan(buf(3, 2, 1)), &
                    "vanished layer must read NaN, not 0")
         if (allocated(error)) exit checks
         call check(error, abs(buf(4, 2, 1) - T_REF) < 1.0e-9_wp, &
                    "an untouched neighbour cell must be unaffected")
      end block checks
      call state%destroy()
   end subroutine test_fill_vanished_nan

   subroutine test_fill_ssh(error)
      !! Stamp h_layer (column-sum equivalent to barotropic h) and b
      !! with known values; call SSH fill via step.  Tests the
      !! multilayer-aware fill path: SSH = sum_k(h_layer) - b.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      integer :: i, j, k
      real(wp) :: h_total
      checks: block

         call setup_state(grid, state)

         do j = 1, NY
            do i = 1, NX
               h_total = 10.0_wp + real(i, wp)
               state%barotropic%h(i, j) = h_total          ! barotropic mirror
               state%barotropic%b(i, j) = real(i, wp)
               do k = 1, NZ
                  state%multilayer%h_layer(i, j, k) = h_total/real(NZ, wp)
               end do
            end do
         end do

         call state%diag%register("SSH", units="m", fill=fill_ssh, &
                                  n1=NX, n2=NY, n3=1, &
                                  time_op=DIAG_OP_INSTANT, dt_out=1.0_wp)

         call ocean_state_enter_data(state)
         call state%diag%step(state, dt=2.0_wp, t=2.0_wp)
         call ocean_state_exit_data(state)

         ! SSH = h - b = 10 everywhere.
         call check(error, abs(state%diag%vars(1)%output_buffer(3, 2, 1) - 10.0_wp) < 1.0e-10_wp, &
                    "SSH at (3,2) should be h(3,2) - b(3,2) = 10")
         if (allocated(error)) exit checks
         call check(error, abs(minval(state%diag%vars(1)%output_buffer) - 10.0_wp) < 1.0e-10_wp, &
                    "SSH min should be 10")
         if (allocated(error)) exit checks
         call check(error, abs(maxval(state%diag%vars(1)%output_buffer) - 10.0_wp) < 1.0e-10_wp, &
                    "SSH max should be 10")

      end block checks
      call state%destroy()
   end subroutine test_fill_ssh

   subroutine test_fill_ke(error)
      !! Uniform u=2, v=3 → cell-centred KE = 0.5*(4 + 9) = 6.5 every cell.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      real(wp), parameter :: U0 = 2.0_wp, V0 = 3.0_wp
      real(wp) :: expected_ke

      call setup_state(grid, state)
      state%multilayer%u_face_x_layer = U0
      state%multilayer%v_face_y_layer = V0
      state%multilayer%h_layer = 10.0_wp

      call state%diag%register("KE", units="m2 s-2", fill=fill_ke, &
                               n1=NX, n2=NY, n3=NZ, &
                               time_op=DIAG_OP_INSTANT, dt_out=1.0_wp)

      call ocean_state_enter_data(state)
      call state%diag%step(state, dt=2.0_wp, t=2.0_wp)
      call ocean_state_exit_data(state)

      expected_ke = 0.5_wp*(U0*U0 + V0*V0)
      call check(error, &
                 abs(state%diag%vars(1)%output_buffer(3, 2, 2) - expected_ke) < 1.0e-10_wp, &
                 "KE should be 0.5*(u^2 + v^2) at every interior cell")

      call state%destroy()
   end subroutine test_fill_ke

   subroutine test_set_z_levels(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      real(wp), parameter :: Z(4) = [1.0_wp, 5.0_wp, 20.0_wp, 100.0_wp]
      checks: block

         call setup_state(grid, state)
         call state%diag%set_output_z_levels(Z)
         call check(error, state%diag%nz_out == 4, "nz_out should be 4")
         if (allocated(error)) exit checks
         call check(error, abs(state%diag%z_out(3) - 20.0_wp) < 1.0e-12_wp, &
                    "z_out(3) should be 20.0")

      end block checks
      call state%destroy()
   end subroutine test_set_z_levels

   subroutine test_zfixed_alloc(error)
      !! Z_FIXED register: both layer_buffer and output_buffer allocated,
      !! output_buffer has nz_out third dim.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      real(wp), parameter :: Z(2) = [5.0_wp, 50.0_wp]
      checks: block

         call setup_state(grid, state)
         call state%diag%set_output_z_levels(Z)
         call state%diag%register("foo", units="m", fill=dummy_fill, &
                                  n1=NX, n2=NY, n3=NZ, &
                                  output_vgrid=DIAG_VGRID_Z_FIXED, &
                                  remap=remap_layer_to_z)

         call check(error, allocated(state%diag%vars(1)%layer_buffer), &
                    "layer_buffer should be allocated for Z_FIXED")
         if (allocated(error)) exit checks
         call check(error, size(state%diag%vars(1)%layer_buffer, 3) == NZ, &
                    "layer_buffer should have nz_native depth")
         if (allocated(error)) exit checks
         call check(error, size(state%diag%vars(1)%output_buffer, 3) == 2, &
                    "output_buffer should have nz_out depth")
         if (allocated(error)) exit checks
         call check(error, associated(state%diag%vars(1)%remap), &
                    "remap should be bound")

      end block checks
      call state%destroy()
   end subroutine test_zfixed_alloc

   subroutine test_remap_identity(error)
      !! CONSERVATIVE z-remap.  Uniform h_layer=10 (NZ=3 ⇒ H=30); z_out
      !! are target INTERFACE depths [5,15,25] with an implicit surface
      !! at 0, so output cells span [0,5], [5,15], [15,25].  Source
      !! TOP-DOWN values are 30 (k=NZ), 20, 10 (k=1).  Cell averages:
      !!   [0,5]   ⊂ top layer            → 30
      !!   [5,15]  half top, half mid     → (5·30 + 5·20)/10 = 25
      !!   [15,25] half mid, half bottom  → (5·20 + 5·10)/10 = 15
      !! These partial-cell averages are the PCM (piecewise-constant)
      !! overlap semantics, so this test pins `diag_remap_method = REMAP_PCM`
      !! (and doubles as a "PCM is still selectable + correct" check); the
      !! production default is PPM, which reconstructs the in-layer profile
      !! and would return the higher-order linear averages instead.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      real(wp) :: layer_buf(NX, NY, NZ), output_buf(NX, NY, 3)
      real(wp), parameter :: H_LAYER = 10.0_wp
      real(wp), parameter :: Z_IFACE(3) = [5.0_wp, 15.0_wp, 25.0_wp]
      integer :: i, j, k
      checks: block

         call set_diag_remap_method(REMAP_PCM)
         call setup_state(grid, state)
         state%multilayer%h_layer = H_LAYER

         ! Layer values: layer(:,:,k) = real(k) so we can distinguish them.
         do k = 1, NZ
            layer_buf(:, :, k) = real(k, wp)*10.0_wp
         end do

         call remap_layer_to_z(state, Z_IFACE, layer_buf, output_buf, .false.)

         call check(error, abs(output_buf(2, 2, 1) - 30.0_wp) < 1.0e-10_wp, &
                    "cell [0,5] cell-average should be 30")
         if (allocated(error)) exit checks
         call check(error, abs(output_buf(2, 2, 2) - 25.0_wp) < 1.0e-10_wp, &
                    "cell [5,15] cell-average should be 25")
         if (allocated(error)) exit checks
         call check(error, abs(output_buf(2, 2, 3) - 15.0_wp) < 1.0e-10_wp, &
                    "cell [15,25] cell-average should be 15")

      end block checks
      call set_diag_remap_method(REMAP_PPM)   ! restore the production default
      call state%destroy()
   end subroutine test_remap_identity

   subroutine test_remap_clamp_top(error)
      !! Single zero-thickness target cell [0,0] (z_out=0) ⇒ no overlap
      !! ⇒ the conservative remap reads 0 (above-column = missing).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      real(wp) :: layer_buf(NX, NY, NZ), output_buf(NX, NY, 1)
      real(wp), parameter :: Z_OUT(1) = [0.0_wp]
      integer :: k

      call setup_state(grid, state)
      state%multilayer%h_layer = 10.0_wp
      do k = 1, NZ
         layer_buf(:, :, k) = real(k, wp)*10.0_wp
      end do

      call remap_layer_to_z(state, Z_OUT, layer_buf, output_buf, .false.)
      call check(error, abs(output_buf(2, 2, 1)) < 1.0e-10_wp, &
                 "zero-thickness target cell [0,0] should read 0")

      call state%destroy()
   end subroutine test_remap_clamp_top

   subroutine test_remap_clamp_bot(error)
      !! Single target cell [0,100] clipped to the column total H=30 ⇒
      !! spans the whole column, so the cell-average is the column mean
      !! (30 + 20 + 10)/3 = 20.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      real(wp) :: layer_buf(NX, NY, NZ), output_buf(NX, NY, 1)
      real(wp), parameter :: Z_OUT(1) = [100.0_wp]
      integer :: k

      call setup_state(grid, state)
      state%multilayer%h_layer = 10.0_wp
      do k = 1, NZ
         layer_buf(:, :, k) = real(k, wp)*10.0_wp
      end do

      call remap_layer_to_z(state, Z_OUT, layer_buf, output_buf, .false.)
      call check(error, abs(output_buf(2, 2, 1) - 20.0_wp) < 1.0e-10_wp, &
                 "whole-column cell-average should be 20")

      call state%destroy()
   end subroutine test_remap_clamp_bot

   subroutine test_remap_e2e_temp(error)
      !! End-to-end: configure z_out, register T with Z_FIXED + remap,
      !! step.  Verify output_buffer has the right shape and contains
      !! the remapped temperature.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      real(wp), parameter :: H_LAYER = 5.0_wp
      real(wp), parameter :: Z_OUT(2) = [2.5_wp, 12.5_wp]
      integer :: it, k
      checks: block

         call setup_state(grid, state)
         state%multilayer%h_layer = H_LAYER

         ! T values: layer k=3 = 20°C, k=2 = 10°C, k=1 = 5°C.  T is stored
         ! as hTr = h * T, so hTr = H_LAYER * T per layer.
         it = state%multilayer%idx_temperature
         do k = 1, NZ
            state%multilayer%tracers(it)%hTr(:, :, k) = H_LAYER*real(k*5, wp)
         end do

         call state%diag%set_output_z_levels(Z_OUT)
         call state%diag%register("T", units="degC", fill=fill_temperature, &
                                  n1=NX, n2=NY, n3=NZ, &
                                  output_vgrid=DIAG_VGRID_Z_FIXED, &
                                  remap=remap_layer_to_z, &
                                  time_op=DIAG_OP_INSTANT, dt_out=1.0_wp)

         call ocean_state_enter_data(state)
         call state%diag%step(state, dt=2.0_wp, t=2.0_wp)
         call ocean_state_exit_data(state)

         call check(error, size(state%diag%vars(1)%output_buffer, 3) == 2, &
                    "output_buffer should have nz_out=2 depth")
         if (allocated(error)) exit checks
         ! Conservative cells (implicit 0 surface): [0,2.5] ⊂ top layer
         ! (T=15) → 15.
         call check(error, abs(state%diag%vars(1)%output_buffer(2, 2, 1) - 15.0_wp) < 1.0e-10_wp, &
                    "cell [0,2.5] cell-average should be top-layer T=15")
         if (allocated(error)) exit checks
         ! Cell [2.5,12.5] spans 2.5 top + 5 mid + 2.5 bottom (clipped to
         ! H=15): (2.5·15 + 5·10 + 2.5·5)/10 = 10.
         call check(error, abs(state%diag%vars(1)%output_buffer(2, 2, 2) - 10.0_wp) < 1.0e-10_wp, &
                    "cell [2.5,12.5] cell-average should be 10")

      end block checks
      call state%destroy()
   end subroutine test_remap_e2e_temp

   subroutine test_mean_basic(error)
      !! MEAN with dt=1, dt_out=4 → fill called every step; on fire
      !! samples are [1, 2, 3, 4] → mean = 2.5.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      integer :: step
      checks: block

         counting_fill_count = 0
         call setup_state(grid, state)
         call state%diag%register("foo", units="-", fill=counting_fill, &
                                  n1=NX, n2=NY, n3=1, &
                                  time_op=DIAG_OP_MEAN, dt_out=4.0_wp)

         do step = 1, 4
            call state%diag%step(state, dt=1.0_wp, t=real(step, wp))
         end do

         ! After 4 steps: dt_accum reaches 4.0 (== dt_out) on the 4th call.
         ! Samples are 1, 2, 3, 4 → mean = 2.5.
         call check(error, abs(state%diag%vars(1)%output_buffer(2, 2, 1) - 2.5_wp) < 1.0e-10_wp, &
                    "mean of [1,2,3,4] should be 2.5")
         if (allocated(error)) exit checks
         ! Post-fire: n_accum and dt_accum should be reset.
         call check(error, state%diag%vars(1)%n_accum == 0, &
                    "n_accum should reset after fire")
         if (allocated(error)) exit checks
         call check(error, abs(state%diag%vars(1)%dt_accum) < 1.0e-12_wp, &
                    "dt_accum should reset after fire")

      end block checks
      call state%destroy()
   end subroutine test_mean_basic

   subroutine test_mean_reset(error)
      !! Run two consecutive windows; the second window's mean is
      !! over [5, 6, 7, 8] = 6.5, NOT the cumulative mean since start.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      integer :: step

      counting_fill_count = 0
      call setup_state(grid, state)
      call state%diag%register("foo", units="-", fill=counting_fill, &
                               n1=NX, n2=NY, n3=1, &
                               time_op=DIAG_OP_MEAN, dt_out=4.0_wp)

      do step = 1, 8
         call state%diag%step(state, dt=1.0_wp, t=real(step, wp))
      end do

      call check(error, abs(state%diag%vars(1)%output_buffer(2, 2, 1) - 6.5_wp) < 1.0e-10_wp, &
                 "second-window mean of [5,6,7,8] should be 6.5")

      call state%destroy()
   end subroutine test_mean_reset

   subroutine test_max(error)
      !! Samples [1, 2, 3, 4] → max = 4.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      integer :: step

      counting_fill_count = 0
      call setup_state(grid, state)
      call state%diag%register("foo", units="-", fill=counting_fill, &
                               n1=NX, n2=NY, n3=1, &
                               time_op=DIAG_OP_MAX, dt_out=4.0_wp)

      do step = 1, 4
         call state%diag%step(state, dt=1.0_wp, t=real(step, wp))
      end do

      call check(error, abs(state%diag%vars(1)%output_buffer(2, 2, 1) - 4.0_wp) < 1.0e-10_wp, &
                 "max of [1,2,3,4] should be 4")

      call state%destroy()
   end subroutine test_max

   subroutine test_min(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      integer :: step

      counting_fill_count = 0
      call setup_state(grid, state)
      call state%diag%register("foo", units="-", fill=counting_fill, &
                               n1=NX, n2=NY, n3=1, &
                               time_op=DIAG_OP_MIN, dt_out=4.0_wp)

      do step = 1, 4
         call state%diag%step(state, dt=1.0_wp, t=real(step, wp))
      end do

      call check(error, abs(state%diag%vars(1)%output_buffer(2, 2, 1) - 1.0_wp) < 1.0e-10_wp, &
                 "min of [1,2,3,4] should be 1")

      call state%destroy()
   end subroutine test_min

   subroutine test_instant_no_accum(error)
      !! INSTANT mode: fill must NOT run except on cadence-fire steps.
      !! After 3 sub-cadence steps, counter stays at 0.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      integer :: step
      checks: block

         counting_fill_count = 0
         call setup_state(grid, state)
         call state%diag%register("foo", units="-", fill=counting_fill, &
                                  n1=NX, n2=NY, n3=1, &
                                  time_op=DIAG_OP_INSTANT, dt_out=10.0_wp)

         do step = 1, 3
            call state%diag%step(state, dt=1.0_wp, t=real(step, wp))
         end do
         call check(error, counting_fill_count == 0, &
                    "INSTANT fill should NOT run below cadence")
         if (allocated(error)) exit checks

         ! One more step crosses dt=10 → fire, counter = 1.
         do step = 4, 10
            call state%diag%step(state, dt=1.0_wp, t=real(step, wp))
         end do
         call check(error, counting_fill_count == 1, &
                    "INSTANT fill should run exactly once at cadence")

      end block checks
      call state%destroy()
   end subroutine test_instant_no_accum

   subroutine test_mean_dt_weighted(error)
      !! Phase A regression: with variable dt, MEAN must be dt-weighted
      !! (the true time integral of the field divided by total time),
      !! not count-weighted (Σ samples / n samples).  Catches a return
      !! to the pre-Phase-A behaviour.
      !!
      !! Schedule: dt sequence {1, 2, 4, 8, 16}, sample at every step
      !! via `counting_fill` so samples are {1, 2, 3, 4, 5}.
      !!   * Count-weighted mean  = (1+2+3+4+5) / 5         = 3
      !!   * dt-weighted mean     = Σ(s·dt) / Σ(dt)
      !!                          = (1·1+2·2+3·4+4·8+5·16) / 31
      !!                          = 129 / 31  ≈ 4.161290…
      !! Test passes for dt-weighted, fails for count-weighted.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      real(wp), parameter :: DT_OUT = 31.0_wp
      real(wp), parameter :: EXPECTED = 129.0_wp/31.0_wp
      real(wp) :: dt_schedule(5), t_running
      integer :: step
      checks: block

         counting_fill_count = 0
         call setup_state(grid, state)
         call state%diag%register("foo", units="-", fill=counting_fill, &
                                  n1=NX, n2=NY, n3=1, &
                                  time_op=DIAG_OP_MEAN, dt_out=DT_OUT)

         dt_schedule = [1.0_wp, 2.0_wp, 4.0_wp, 8.0_wp, 16.0_wp]
         t_running = 0.0_wp
         do step = 1, 5
            t_running = t_running + dt_schedule(step)
            call state%diag%step(state, dt=dt_schedule(step), t=t_running)
         end do

         call check(error, abs(state%diag%vars(1)%output_buffer(2, 2, 1) - EXPECTED) < 1.0e-10_wp, &
                    "dt-weighted MEAN of {1,2,3,4,5} with dt={1,2,4,8,16} must be 129/31")

      end block checks
      call state%destroy()
   end subroutine test_mean_dt_weighted

   subroutine test_integral_constant_rate(error)
      !! Phase A regression: DIAG_OP_INTEGRAL emits the cumulative
      !! time integral of the sample over the window, without
      !! dividing by `dt_accum`.  For a constant-rate field R over a
      !! variable-dt window of length T, INTEGRAL must equal `R · T`,
      !! NOT `R · n_samples` (which a count-weighted bug produces).
      !!
      !! `dummy_fill` emits the constant 42.0 every step.  dt
      !! schedule {1, 2, 4, 8, 16} ⇒ total elapsed time 31.
      !!   * dt-weighted INTEGRAL  = 42 · 31              = 1302
      !!   * count-weighted (bug)  = 42 · 5               = 210
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      real(wp), parameter :: DT_OUT = 31.0_wp
      real(wp), parameter :: EXPECTED = 42.0_wp*31.0_wp
      real(wp) :: dt_schedule(5), t_running
      integer :: step
      checks: block

         call setup_state(grid, state)
         call state%diag%register("foo", units="-", fill=dummy_fill, &
                                  n1=NX, n2=NY, n3=1, &
                                  time_op=DIAG_OP_INTEGRAL, dt_out=DT_OUT)

         dt_schedule = [1.0_wp, 2.0_wp, 4.0_wp, 8.0_wp, 16.0_wp]
         t_running = 0.0_wp
         do step = 1, 5
            t_running = t_running + dt_schedule(step)
            call state%diag%step(state, dt=dt_schedule(step), t=t_running)
         end do

         call check(error, abs(state%diag%vars(1)%output_buffer(2, 2, 1) - EXPECTED) < 1.0e-10_wp, &
                    "INTEGRAL of constant 42 over dt-sequence summing to 31 must be 1302")

      end block checks
      call state%destroy()
   end subroutine test_integral_constant_rate

   ! ---------------------------------------------------------------------
   ! Phase B-partial — derived diagnostic library
   ! ---------------------------------------------------------------------

   subroutine test_derived_catalog_size(error)
      !! Catalog ships 23 entries: 7 base + 3 sea-ice velocity diags +
      !! 13 ice-shelf-cavity diags (11 melt-interface, gated on
      !! `&ocean_cavity_melt_nml`, and 2 geometry, gated on
      !! `&ocean_cavity_dyn_nml`).  Locks the count so we notice if an
      !! entry is accidentally dropped.
      type(error_type), allocatable, intent(out) :: error
      integer :: n
      n = derived_catalog_size()
      call check(error, n == 23, "derived catalog should hold 23 entries")
   end subroutine test_derived_catalog_size

   subroutine test_derived_h_layer(error)
      !! `h_layer` is a direct copy.  Stamp the state's h_layer with a
      !! known pattern, fire the diagnostic, verify the buffer matches.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      integer :: i, j, k
      real(wp) :: diff
      checks: block
         call setup_state(grid, state)
         do k = 1, NZ
            do j = 1, grid%ny_total
               do i = 1, grid%nx_total
                  state%multilayer%h_layer(i, j, k) = real(10*k + i, wp)
               end do
            end do
         end do
         call register_derived(state, "h_layer", time_op=DIAG_OP_INSTANT, dt_out=1.0_wp)
         call ocean_state_enter_data(state)
         call state%diag%step(state, dt=1.0_wp, t=1.0_wp)
         call ocean_state_exit_data(state)

         diff = maxval(abs(state%diag%vars(1)%output_buffer &
                           - state%multilayer%h_layer(:grid%nx_total, :grid%ny_total, :NZ)))
         call check(error, diff < 1.0e-12_wp, &
                    "h_layer derived diag should match h_layer state to FP")
      end block checks
      call state%destroy()
   end subroutine test_derived_h_layer

   subroutine test_derived_rho_layer(error)
      !! `rho_layer` direct copy from the EOS slot on every LIVE layer, and
      !! the NaN missing-data sentinel on a vanished one: the EOS writes
      !! `rho_0` into a filler on purpose (so the PGF column is not
      !! perturbed), which is a plausible density, not a measurement.
      use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      integer :: i, j, k
      real(wp) :: diff
      checks: block
         call setup_state(grid, state)
         state%multilayer%h_layer = 10.0_wp
         do k = 1, NZ
            do j = 1, grid%ny_total
               do i = 1, grid%nx_total
                  state%multilayer%rho_layer(i, j, k) = 1025.0_wp + 0.1_wp*real(k, wp)
               end do
            end do
         end do
         ! One vanished bed filler, exactly ON the marker.
         state%multilayer%h_layer(3, 2, 1) = H_VANISHED
         call register_derived(state, "rho_layer", time_op=DIAG_OP_INSTANT, dt_out=1.0_wp)
         call ocean_state_enter_data(state)
         call state%diag%step(state, dt=1.0_wp, t=1.0_wp)
         call ocean_state_exit_data(state)

         call check(error, ieee_is_nan(state%diag%vars(1)%output_buffer(3, 2, 1)), &
                    "rho_layer on a vanished layer must be missing (NaN), not rho_0")
         if (allocated(error)) exit checks
         diff = 0.0_wp
         do k = 1, NZ
            do j = 1, grid%ny_total
               do i = 1, grid%nx_total
                  if (i == 3 .and. j == 2 .and. k == 1) cycle
                  diff = max(diff, abs(state%diag%vars(1)%output_buffer(i, j, k) &
                                       - state%multilayer%rho_layer(i, j, k)))
               end do
            end do
         end do
         call check(error, diff < 1.0e-12_wp, &
                    "rho_layer derived diag should match rho_layer state to FP on live layers")
      end block checks
      call state%destroy()
   end subroutine test_derived_rho_layer

   subroutine test_derived_vorticity(error)
      !! Rigid-rotation analytical:
      !!   u(face_x) = -Ω · y_centre_at_face
      !!   v(face_y) = +Ω · x_centre_at_face
      !! ⇒ ζ = ∂v/∂x − ∂u/∂y = 2Ω at every cell centre, to FP.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      real(wp), parameter :: OMEGA = 0.5_wp
      integer :: i, j, k, i0, i1, j0, j1
      real(wp) :: y_cell, x_cell, max_err, expected
      checks: block
         call setup_state(grid, state)
         do k = 1, NZ
            do j = 1, grid%ny_total
               y_cell = (real(j, wp) - 0.5_wp)*grid%dx
               do i = 1, grid%nx_total + 1
                  state%multilayer%u_face_x_layer(i, j, k) = -OMEGA*y_cell
               end do
            end do
         end do
         do k = 1, NZ
            do j = 1, grid%ny_total + 1
               do i = 1, grid%nx_total
                  x_cell = (real(i, wp) - 0.5_wp)*grid%dx
                  state%multilayer%v_face_y_layer(i, j, k) = OMEGA*x_cell
               end do
            end do
         end do

         call register_derived(state, "vorticity_z", time_op=DIAG_OP_INSTANT, dt_out=1.0_wp)
         call ocean_state_enter_data(state)
         call state%diag%step(state, dt=1.0_wp, t=1.0_wp)
         call ocean_state_exit_data(state)

         expected = 2.0_wp*OMEGA
         ! The fill leaves a 1-cell ring of zeros at the edges of the
         ! interior buffer (insufficient stencil).  Check the inner
         ! portion only.
         i0 = 3
         i1 = grid%nx_total - 2
         j0 = 3
         j1 = grid%ny_total - 2
         max_err = maxval(abs(state%diag%vars(1)%output_buffer(i0:i1, j0:j1, :) - expected))
         call check(error, max_err < 1.0e-12_wp, &
                    "vorticity_z on rigid rotation should equal 2Ω to FP")
      end block checks
      call state%destroy()
   end subroutine test_derived_vorticity

   subroutine test_derived_ke_total(error)
      !! Uniform flow + uniform thickness:
      !!   KE_total = Σ_k 0.5 · h · (u² + v²) = 0.5 · H · (u² + v²)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      real(wp), parameter :: U0 = 0.2_wp, V0 = 0.3_wp, H_LAYER = 50.0_wp
      real(wp) :: expected, err
      checks: block
         call setup_state(grid, state)
         state%multilayer%h_layer = H_LAYER
         state%multilayer%u_face_x_layer = U0
         state%multilayer%v_face_y_layer = V0

         call register_derived(state, "ke_total", time_op=DIAG_OP_INSTANT, dt_out=1.0_wp)
         call ocean_state_enter_data(state)
         call state%diag%step(state, dt=1.0_wp, t=1.0_wp)
         call ocean_state_exit_data(state)

         expected = 0.5_wp*H_LAYER*real(NZ, wp)*(U0*U0 + V0*V0)
         err = abs(state%diag%vars(1)%output_buffer(3, 2, 1) - expected)
         call check(error, err < 1.0e-12_wp, &
                    "ke_total on uniform flow should equal 0.5 H (u² + v²) to FP")
      end block checks
      call state%destroy()
   end subroutine test_derived_ke_total

   subroutine test_derived_transport_x(error)
      !! Uniform u + uniform h ⇒ transport_x = u · H per column.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      real(wp), parameter :: U0 = 0.1_wp, H_LAYER = 100.0_wp
      real(wp) :: expected, err
      checks: block
         call setup_state(grid, state)
         state%multilayer%h_layer = H_LAYER
         state%multilayer%u_face_x_layer = U0
         state%multilayer%v_face_y_layer = 0.0_wp

         call register_derived(state, "transport_x", time_op=DIAG_OP_INSTANT, dt_out=1.0_wp)
         call ocean_state_enter_data(state)
         call state%diag%step(state, dt=1.0_wp, t=1.0_wp)
         call ocean_state_exit_data(state)

         expected = U0*H_LAYER*real(NZ, wp)
         err = abs(state%diag%vars(1)%output_buffer(3, 2, 1) - expected)
         call check(error, err < 1.0e-12_wp, &
                    "transport_x on uniform flow should equal u·H to FP")
      end block checks
      call state%destroy()
   end subroutine test_derived_transport_x

   subroutine test_derived_transport_y(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      real(wp), parameter :: V0 = 0.05_wp, H_LAYER = 200.0_wp
      real(wp) :: expected, err
      checks: block
         call setup_state(grid, state)
         state%multilayer%h_layer = H_LAYER
         state%multilayer%u_face_x_layer = 0.0_wp
         state%multilayer%v_face_y_layer = V0

         call register_derived(state, "transport_y", time_op=DIAG_OP_INSTANT, dt_out=1.0_wp)
         call ocean_state_enter_data(state)
         call state%diag%step(state, dt=1.0_wp, t=1.0_wp)
         call ocean_state_exit_data(state)

         expected = V0*H_LAYER*real(NZ, wp)
         err = abs(state%diag%vars(1)%output_buffer(3, 2, 1) - expected)
         call check(error, err < 1.0e-12_wp, &
                    "transport_y on uniform flow should equal v·H to FP")
      end block checks
      call state%destroy()
   end subroutine test_derived_transport_y

   subroutine test_derived_mld(error)
      !! 3-layer column, surface (k=NZ=3) and middle (k=2) at the
      !! same density, bottom (k=1) denser by Δρ > threshold.  MLD
      !! is the cumulative thickness from surface down to the first
      !! layer that crosses the threshold, i.e. h(k=3) + h(k=2).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      real(wp), parameter :: RHO_LIGHT = 1024.0_wp, RHO_DENSE = 1025.0_wp
      real(wp), parameter :: H_TOP = 20.0_wp, H_MID = 30.0_wp, H_BOT = 50.0_wp
      real(wp) :: expected, err
      checks: block
         call setup_state(grid, state)
         state%multilayer%h_layer(:, :, 1) = H_BOT
         state%multilayer%h_layer(:, :, 2) = H_MID
         state%multilayer%h_layer(:, :, 3) = H_TOP
         state%multilayer%rho_layer(:, :, 1) = RHO_DENSE
         state%multilayer%rho_layer(:, :, 2) = RHO_LIGHT
         state%multilayer%rho_layer(:, :, 3) = RHO_LIGHT

         call register_derived(state, "mld_density", time_op=DIAG_OP_INSTANT, dt_out=1.0_wp)
         call ocean_state_enter_data(state)
         call state%diag%step(state, dt=1.0_wp, t=1.0_wp)
         call ocean_state_exit_data(state)

         expected = H_TOP + H_MID
         err = abs(state%diag%vars(1)%output_buffer(3, 2, 1) - expected)
         call check(error, err < 1.0e-12_wp, &
                    "mld_density should be h(surface) + h(mid) for the configured step profile")
      end block checks
      call state%destroy()
   end subroutine test_derived_mld

   subroutine test_derived_mld_vanished(error)
      !! A vanished bed filler cannot mark the mixed-layer base.  The EOS
      !! puts the reference density `rho_0` into a filler, which here is
      !! denser than the (uniform, light) water above it by more than the
      !! threshold; a scan that trusted it reported MLD = h(top) + h(mid)
      !! — a pycnocline made of nothing.  With the filler ignored the
      !! column has no crossing, so MLD is the full column depth,
      !! filler thickness included.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      real(wp), parameter :: RHO_LIGHT = 1024.0_wp, RHO_FILLER = 1035.0_wp
      real(wp), parameter :: H_TOP = 20.0_wp, H_MID = 30.0_wp
      real(wp) :: expected, err
      checks: block
         call setup_state(grid, state)
         state%multilayer%h_layer(:, :, 1) = H_VANISHED
         state%multilayer%h_layer(:, :, 2) = H_MID
         state%multilayer%h_layer(:, :, 3) = H_TOP
         state%multilayer%rho_layer(:, :, 1) = RHO_FILLER
         state%multilayer%rho_layer(:, :, 2) = RHO_LIGHT
         state%multilayer%rho_layer(:, :, 3) = RHO_LIGHT

         call register_derived(state, "mld_density", time_op=DIAG_OP_INSTANT, dt_out=1.0_wp)
         call ocean_state_enter_data(state)
         call state%diag%step(state, dt=1.0_wp, t=1.0_wp)
         call ocean_state_exit_data(state)

         expected = H_TOP + H_MID + H_VANISHED
         err = abs(state%diag%vars(1)%output_buffer(3, 2, 1) - expected)
         call check(error, err < 1.0e-9_wp, &
                    "mld_density must not cross at a vanished filler (MLD = full column)")
      end block checks
      call state%destroy()
   end subroutine test_derived_mld_vanished

   ! ---------------------------------------------------------------------
   ! Phase C — region masks
   ! ---------------------------------------------------------------------

   subroutine test_mask_global(error)
      !! Global mask covers the full grid with weight=1 everywhere
      !! and total_area = nx · ny · dx · dy.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(diag_mask_t) :: m
      real(wp) :: expected_area
      checks: block
         call grid%init(NX, NY, 1, DX, DX)
         m = diag_mask_global(grid)
         call check(error, m%nx == grid%nx_total, "mask nx mismatch")
         if (allocated(error)) exit checks
         call check(error, m%ny == grid%ny_total, "mask ny mismatch")
         if (allocated(error)) exit checks
         call check(error, all(m%weight == 1.0_wp), "global mask weight must be 1 everywhere")
         if (allocated(error)) exit checks
         expected_area = real(grid%nx_total*grid%ny_total, wp)*grid%dx*grid%dx
         call check(error, abs(m%total_area - expected_area) < 1.0e-12_wp, &
                    "global mask total_area must equal nx · ny · dx · dy")
      end block checks
   end subroutine test_mask_global

   subroutine test_mask_bbox(error)
      !! Bbox mask: weight=1 inside the configured rectangle, 0
      !! outside.  Indices are 1-based and inclusive.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(diag_mask_t) :: m
      integer :: i0, i1, j0, j1, count_inside
      real(wp) :: expected_area
      checks: block
         call grid%init(NX, NY, 1, DX, DX)
         i0 = 2
         i1 = 4
         j0 = 2
         j1 = 3
         m = diag_mask_bbox(grid, i0, i1, j0, j1)
         count_inside = (i1 - i0 + 1)*(j1 - j0 + 1)
         call check(error, count(m%weight > 0.5_wp) == count_inside, &
                    "bbox mask should cover exactly (i1-i0+1)·(j1-j0+1) cells")
         if (allocated(error)) exit checks
         call check(error, m%weight(i0, j0) == 1.0_wp .and. m%weight(i1, j1) == 1.0_wp, &
                    "bbox mask should have weight=1 at the corners")
         if (allocated(error)) exit checks
         call check(error, m%weight(i0 - 1, j0) == 0.0_wp .and. &
                    m%weight(i1 + 1, j1) == 0.0_wp, &
                    "bbox mask should have weight=0 just outside")
         if (allocated(error)) exit checks
         expected_area = real(count_inside, wp)*grid%dx*grid%dx
         call check(error, abs(m%total_area - expected_area) < 1.0e-12_wp, &
                    "bbox mask total_area should equal n_cells · dx · dy")
      end block checks
   end subroutine test_mask_bbox

   subroutine test_mask_h_section(error)
      !! Horizontal section: single cell-row strip.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(diag_mask_t) :: m
      integer :: j_row, i0, i1
      checks: block
         call grid%init(NX, NY, 1, DX, DX)
         j_row = 3
         i0 = 2
         i1 = 5
         m = diag_mask_h_section(grid, j_row, i0, i1)
         call check(error, count(m%weight > 0.5_wp) == (i1 - i0 + 1), &
                    "h_section should cover exactly i1-i0+1 cells")
         if (allocated(error)) exit checks
         call check(error, all(m%weight(i0:i1, j_row) == 1.0_wp), &
                    "h_section should have weight=1 along the row")
         if (allocated(error)) exit checks
         call check(error, all(m%weight(:, j_row - 1) == 0.0_wp) .and. &
                    all(m%weight(:, j_row + 1) == 0.0_wp), &
                    "h_section should not bleed into adjacent rows")
      end block checks
   end subroutine test_mask_h_section

   subroutine test_mask_v_section(error)
      !! Vertical section: single cell-column strip.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(diag_mask_t) :: m
      integer :: i_col, j0, j1
      checks: block
         call grid%init(NX, NY, 1, DX, DX)
         i_col = 4
         j0 = 1
         j1 = 3
         m = diag_mask_v_section(grid, i_col, j0, j1)
         call check(error, count(m%weight > 0.5_wp) == (j1 - j0 + 1), &
                    "v_section should cover exactly j1-j0+1 cells")
         if (allocated(error)) exit checks
         call check(error, all(m%weight(i_col, j0:j1) == 1.0_wp), &
                    "v_section should have weight=1 along the column")
         if (allocated(error)) exit checks
         call check(error, all(m%weight(i_col - 1, :) == 0.0_wp) .and. &
                    all(m%weight(i_col + 1, :) == 0.0_wp), &
                    "v_section should not bleed into adjacent columns")
      end block checks
   end subroutine test_mask_v_section

   subroutine test_mask_zeros_outside(error)
      !! End-to-end mask test: register a derived diagnostic with a
      !! bbox mask, fire MEAN.  Cells inside the bbox match the
      !! state value; cells outside the bbox are zero in the
      !! accumulator (since they were multiplied by mask=0 at fold).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      type(diag_mask_t) :: m
      real(wp), parameter :: H_LAYER = 50.0_wp
      integer :: i0, i1, j0, j1
      real(wp) :: err_inside, err_outside
      checks: block
         call setup_state(grid, state)
         state%multilayer%h_layer = H_LAYER

         i0 = 3
         i1 = 5
         j0 = 2
         j1 = 3
         m = diag_mask_bbox(grid, i0, i1, j0, j1)
         ! INTEGRAL with one substep of dt=1 ⇒ accumulator value
         ! equals the masked sample directly (no averaging).
         call state%diag%register("h_masked", units="m", fill=fill_h_layer_local, &
                                  n1=grid%nx_total, n2=grid%ny_total, n3=NZ, &
                                  time_op=DIAG_OP_INTEGRAL, dt_out=1.0_wp, mask=m)
         call ocean_state_enter_data(state)
         call state%diag%step(state, dt=1.0_wp, t=1.0_wp)
         call ocean_state_exit_data(state)

         err_inside = maxval(abs(state%diag%vars(1)%output_buffer(i0:i1, j0:j1, :) - H_LAYER))
         err_outside = maxval(abs(state%diag%vars(1)%output_buffer(1, 1, :)))
         call check(error, err_inside < 1.0e-12_wp, &
                    "masked INTEGRAL inside bbox should equal H_LAYER · 1 = 50")
         if (allocated(error)) exit checks
         call check(error, err_outside < 1.0e-12_wp, &
                    "masked INTEGRAL outside bbox should be zero")
      end block checks
      call state%destroy()
   end subroutine test_mask_zeros_outside

   ! ---------------------------------------------------------------------
   ! Phase D v1 — global integral budgets
   ! ---------------------------------------------------------------------

   subroutine test_budget_total_mass(error)
      !! `budget_total_mass` returns Σ h_layer · dA over the masked
      !! interior.  Verify against a hand-summed reference.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      type(diag_mask_t) :: m
      real(wp), parameter :: H_LAYER = 50.0_wp
      real(wp) :: total, expected
      real(wp), allocatable :: areaT(:, :)
      checks: block
         call setup_state(grid, state)
         state%multilayer%h_layer = H_LAYER
         m = diag_mask_global(grid)
         allocate (areaT(grid%nx_total, grid%ny_total), source=grid%dx*grid%dy)
         total = budget_total_mass(state%multilayer, m, areaT)
         expected = real(grid%nx_total*grid%ny_total, wp)*grid%dx*grid%dx* &
                    H_LAYER*real(NZ, wp)
         call check(error, abs(total - expected) < 1.0e-9_wp*expected, &
                    "total mass should equal nx · ny · dx · dy · H · NZ")
      end block checks
      call state%destroy()
   end subroutine test_budget_total_mass

   subroutine test_budget_total_tracer(error)
      !! `budget_total_tracer` returns Σ hTr · dA over the masked
      !! interior.  Set salinity to a uniform value and verify the
      !! integral matches.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      type(diag_mask_t) :: m
      real(wp), parameter :: S_REF = 35.0_wp, H_LAYER = 10.0_wp
      real(wp) :: total, expected
      real(wp), allocatable :: areaT(:, :)
      integer :: it_S
      checks: block
         call setup_state(grid, state)
         it_S = state%multilayer%idx_salinity
         call check(error, it_S > 0, "salinity tracer must be registered")
         if (allocated(error)) exit checks
         state%multilayer%h_layer = H_LAYER
         state%multilayer%tracers(it_S)%hTr = S_REF*H_LAYER
         m = diag_mask_global(grid)
         allocate (areaT(grid%nx_total, grid%ny_total), source=grid%dx*grid%dy)
         total = budget_total_tracer(state%multilayer, it_S, m, areaT)
         expected = real(grid%nx_total*grid%ny_total, wp)*grid%dx*grid%dx* &
                    S_REF*H_LAYER*real(NZ, wp)
         call check(error, abs(total - expected) < 1.0e-9_wp*expected, &
                    "total salt should equal nx · ny · dx · dy · S · H · NZ")
      end block checks
      call state%destroy()
   end subroutine test_budget_total_tracer

   subroutine test_budget_total_ke(error)
      !! `budget_total_ke` = Σ 0.5 · h · (u² + v²) · dA at cell centres.
      !! With u, v uniform on all faces the cell-centre value equals
      !! the face value, so KE = 0.5 · H_total · (u² + v²) · area.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      type(diag_mask_t) :: m
      real(wp), parameter :: U0 = 0.3_wp, V0 = 0.4_wp, H_LAYER = 25.0_wp
      real(wp) :: total, expected, nx_phys_area
      real(wp), allocatable :: areaT(:, :)
      checks: block
         call setup_state(grid, state)
         state%multilayer%h_layer = H_LAYER
         state%multilayer%u_face_x_layer = U0
         state%multilayer%v_face_y_layer = V0
         m = diag_mask_global(grid)
         allocate (areaT(grid%nx_total, grid%ny_total), source=grid%dx*grid%dy)
         total = budget_total_ke(state%multilayer, m, areaT)
         ! The kernel iterates 1..nx_total, 1..ny_total — the size
         ! mins evaluate to the full grid given u-face is (nx+1, ny)
         ! and v-face is (nx, ny+1).
         nx_phys_area = real(grid%nx_total*grid%ny_total, wp)*grid%dx*grid%dx
         expected = nx_phys_area*H_LAYER*real(NZ, wp)*0.5_wp*(U0*U0 + V0*V0)
         call check(error, abs(total - expected) < 1.0e-9_wp*expected, &
                    "total KE should equal 0.5 · area · H · (u² + v²)")
      end block checks
      call state%destroy()
   end subroutine test_budget_total_ke

   subroutine test_budget_snapshot(error)
      !! `init_snapshot(state%multilayer)` records `values_init`.
      !! `ocean_budgets_t` carries no `ocean_state_t` slot (PR-8 — it
      !! is a verification instrument, not a production reporter), so
      !! the test constructs its own local instance, per the idiom in
      !! `tests/test_ocean_surface_flux.F90`.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      type(ocean_budgets_t) :: budgets
      real(wp), parameter :: H_LAYER = 10.0_wp
      checks: block
         call setup_state(grid, state)
         call budgets%init(grid)
         state%multilayer%h_layer = H_LAYER
         call budgets%init_snapshot(state%multilayer)
         call check(error, budgets%has_snapshot, &
                    "has_snapshot should be true after init_snapshot")
         if (allocated(error)) exit checks
         call check(error, budgets%values_init(BUDGET_MASS) > 0.0_wp, &
                    "mass snapshot should be positive")
      end block checks
      call budgets%destroy()
      call state%destroy()
   end subroutine test_budget_snapshot

   subroutine test_budget_no_drift(error)
      !! At-rest snapshot followed by no-op evaluations should leave
      !! `values - values_init` at zero to FP (mass / salt / heat).
      !! Confirms the snapshot path is consistent with the evaluation
      !! path (would catch a sign / units bug in one or the other).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      type(ocean_budgets_t) :: budgets
      real(wp), parameter :: H_LAYER = 30.0_wp, S_REF = 35.0_wp, T_REF = 10.0_wp
      integer :: it_S, it_T
      real(wp) :: drift_mass, drift_S, drift_T
      checks: block
         call setup_state(grid, state)
         call budgets%init(grid)
         it_S = state%multilayer%idx_salinity
         it_T = state%multilayer%idx_temperature
         state%multilayer%h_layer = H_LAYER
         if (it_S > 0) state%multilayer%tracers(it_S)%hTr = S_REF*H_LAYER
         if (it_T > 0) state%multilayer%tracers(it_T)%hTr = T_REF*H_LAYER

         call budgets%init_snapshot(state%multilayer)
         call budgets%evaluate(state%multilayer)

         drift_mass = budgets%values(BUDGET_MASS) - budgets%values_init(BUDGET_MASS)
         drift_S = budgets%values(BUDGET_SALT_TOTAL) - budgets%values_init(BUDGET_SALT_TOTAL)
         drift_T = budgets%values(BUDGET_HEAT_TOTAL) - budgets%values_init(BUDGET_HEAT_TOTAL)

         call check(error, abs(drift_mass) < 1.0e-9_wp*budgets%values_init(BUDGET_MASS), &
                    "mass drift should be FP after no-op evaluation")
         if (allocated(error)) exit checks
         call check(error, abs(drift_S) < 1.0e-9_wp*max(budgets%values_init(BUDGET_SALT_TOTAL), 1.0_wp), &
                    "salt drift should be FP after no-op evaluation")
         if (allocated(error)) exit checks
         call check(error, abs(drift_T) < 1.0e-9_wp*max(budgets%values_init(BUDGET_HEAT_TOTAL), 1.0_wp), &
                    "heat drift should be FP after no-op evaluation")
      end block checks
      call budgets%destroy()
      call state%destroy()
   end subroutine test_budget_no_drift

   subroutine test_budget_register_contributor(error)
      !! `register_contributor` appends a new slot, stores name +
      !! quantity, and binds the pointer to the caller's array.
      !! `ocean_state_init` auto-registers the continuity mass
      !! contributor, so the synthetic registration here lands at
      !! `n_before + 1` rather than slot 1.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      type(ocean_budgets_t) :: budgets
      real(wp), allocatable, target :: synth(:, :, :)
      integer :: n_before, idx
      checks: block
         call setup_state(grid, state)
         call budgets%init(grid)
         allocate (synth(grid%nx_total, grid%ny_total, NZ), source=0.0_wp)
         n_before = budgets%n_contributors

         call budgets%register_contributor("synthetic_src", BUDGET_MASS, synth)
         call check(error, budgets%n_contributors == n_before + 1, &
                    "register_contributor should append one slot")
         if (allocated(error)) exit checks
         idx = budgets%n_contributors
         call check(error, trim(budgets%contributors(idx)%name) == "synthetic_src", &
                    "contributor name should be stored")
         if (allocated(error)) exit checks
         call check(error, budgets%contributors(idx)%quantity == BUDGET_MASS, &
                    "contributor quantity should be BUDGET_MASS")
         if (allocated(error)) exit checks
         call check(error, associated(budgets%contributors(idx)%per_cell, synth), &
                    "per_cell pointer should be bound to the caller's array")
         if (allocated(error)) exit checks
         call check(error, budgets%contributors(idx)%is_active, &
                    "contributor should be marked active after register")
      end block checks
      call budgets%destroy()
      call state%destroy()
   end subroutine test_budget_register_contributor

   subroutine test_budget_drain_lifecycle(error)
      !! `drain_contributors` integrates `per_cell · weight · dA` summed
      !! over k into `total_integrated`, then zeroes `per_cell` so the
      !! kernel's next-step writes start from zero.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      type(ocean_budgets_t) :: budgets
      real(wp), allocatable, target :: synth(:, :, :)
      real(wp), parameter :: C0 = 0.25_wp
      real(wp) :: dA, expected, sum_weight
      integer :: i, j, idx
      checks: block
         call setup_state(grid, state)
         call budgets%init(grid)
         dA = grid%dx*grid%dy
         allocate (synth(grid%nx_total, grid%ny_total, NZ), source=C0)
         call budgets%register_contributor("synthetic_src", BUDGET_MASS, synth)
         idx = budgets%n_contributors
         call budgets%drain_contributors()

         sum_weight = 0.0_wp
         do j = 1, budgets%mask%ny
            do i = 1, budgets%mask%nx
               sum_weight = sum_weight + budgets%mask%weight(i, j)
            end do
         end do
         expected = C0*sum_weight*dA*real(NZ, wp)

         call check(error, abs(budgets%contributors(idx)%total_integrated - expected) &
                    < 1.0e-12_wp*max(abs(expected), 1.0_wp), &
                    "total_integrated should equal C0 · Σ(weight · dA) · NZ")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(synth)) < 1.0e-15_wp, &
                    "per_cell should be zeroed after drain")
      end block checks
      call budgets%destroy()
      call state%destroy()
   end subroutine test_budget_drain_lifecycle

   subroutine test_budget_drain_accumulates(error)
      !! Successive drains accumulate into `total_integrated` rather
      !! than overwriting it — the run-cumulative RHS invariant.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      type(ocean_budgets_t) :: budgets
      real(wp), allocatable, target :: synth(:, :, :)
      real(wp), parameter :: C1 = 0.10_wp, C2 = 0.40_wp
      real(wp) :: dA, expected, sum_weight
      integer :: i, j, idx
      checks: block
         call setup_state(grid, state)
         call budgets%init(grid)
         dA = grid%dx*grid%dy
         allocate (synth(grid%nx_total, grid%ny_total, NZ), source=0.0_wp)
         call budgets%register_contributor("synthetic_src", BUDGET_MASS, synth)
         idx = budgets%n_contributors

         synth = C1
         call budgets%drain_contributors()
         synth = C2
         call budgets%drain_contributors()

         sum_weight = 0.0_wp
         do j = 1, budgets%mask%ny
            do i = 1, budgets%mask%nx
               sum_weight = sum_weight + budgets%mask%weight(i, j)
            end do
         end do
         expected = (C1 + C2)*sum_weight*dA*real(NZ, wp)
         call check(error, abs(budgets%contributors(idx)%total_integrated - expected) &
                    < 1.0e-12_wp*max(abs(expected), 1.0_wp), &
                    "total_integrated should accumulate across drains")
      end block checks
      call budgets%destroy()
      call state%destroy()
   end subroutine test_budget_drain_accumulates

   subroutine test_budget_register_idempotent(error)
      !! Re-registering the same `name` overwrites the existing slot
      !! (and resets `total_integrated`) rather than allocating a new
      !! one — important when a kernel re-inits between cases.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      type(ocean_budgets_t) :: budgets
      real(wp), allocatable, target :: synth(:, :, :), synth2(:, :, :)
      integer :: n_before, idx
      checks: block
         call setup_state(grid, state)
         call budgets%init(grid)
         allocate (synth(grid%nx_total, grid%ny_total, NZ), source=0.1_wp)
         allocate (synth2(grid%nx_total, grid%ny_total, NZ), source=0.0_wp)
         n_before = budgets%n_contributors
         call budgets%register_contributor("synthetic_src", BUDGET_MASS, synth)
         idx = budgets%n_contributors
         call budgets%drain_contributors()
         call check(error, budgets%contributors(idx)%total_integrated > 0.0_wp, &
                    "first drain should populate total_integrated")
         if (allocated(error)) exit checks

         call budgets%register_contributor("synthetic_src", BUDGET_MASS, synth2)
         call check(error, budgets%n_contributors == n_before + 1, &
                    "re-registering same name should not append a new slot")
         if (allocated(error)) exit checks
         call check(error, abs(budgets%contributors(idx)%total_integrated) < 1.0e-15_wp, &
                    "re-register should reset total_integrated")
         if (allocated(error)) exit checks
         call check(error, associated(budgets%contributors(idx)%per_cell, synth2), &
                    "per_cell pointer should rebind to new array on re-register")
      end block checks
      call budgets%destroy()
      call state%destroy()
   end subroutine test_budget_register_idempotent

   subroutine fill_h_layer_local(state_handle, buf)
      !! Local copy of fill_h_layer.  Shim writes device-side via
      !! `fill_h_layer_local_impl` because the manager now folds
      !! `output_buffer` on device — a host fill would land in the host
      !! shadow only and the device fold would read zeros.
      class(*), intent(in) :: state_handle
      real(wp), intent(inout) :: buf(:, :, :)
      select type (state => state_handle)
      class is (ocean_state_t)
         call fill_h_layer_local_impl(state%multilayer%h_layer, buf)
      end select
   end subroutine fill_h_layer_local

   subroutine fill_h_layer_local_impl(src, buf)
      real(wp), intent(in)    :: src(:, :, :)
      real(wp), intent(inout) :: buf(:, :, :)
      integer :: i, j, k, nx, ny, nz
      nx = min(size(buf, 1), size(src, 1))
      ny = min(size(buf, 2), size(src, 2))
      nz = min(size(buf, 3), size(src, 3))
      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         buf(i, j, k) = src(i, j, k)
      end do
   end subroutine fill_h_layer_local_impl

end module test_ocean_diag
