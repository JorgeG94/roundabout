!! Vertical-coordinate hygiene: knob plumbing + the fail-loud envelope.
!!
!! Three things this suite pins, each of which was found by reading the
!! `&vcoord_nml` group end-to-end rather than by a failing run:
!!
!!   1. **`zstar_h_surf_target` / `zstar_h_min` actually reach the slot on
!!      the PRODUCTION namelist path.**  Both are copied in `engine_setup`
!!      (`rdb_ocean_engine`), BEFORE the IC seed — which is load-bearing,
!!      because the seed's tail calls `vcoord%build_zref_full(b)` and that
!!      reads `zstar_h_surf_target`.  `configure_ocean_lateral` runs later
!!      and rewrites the OTHER vcoord members (`coord_type`,
!!      `remap_method`, `z_fixed_h_ref`, the RHO targets,
!!      `regrid_time_scale`, `remap_vel_conserve_ke`) while deliberately
!!      leaving this pair alone.  A future tidy-up that "completes" that
!!      block by copying the pair there too would move the assignment to
!!      AFTER the seed and silently un-anchor every ZSTAR_FULL reference
!!      table; this test is the tripwire for that.
!!
!!   2. **`vcoord_type = "zsigma"` is refused** — its deep branch reads
!!      `z_ref_global` as metres while the only writer fills it with the
!!      dimensionless `k/nz`.  See
!!      `test_ocean_vcoord_interface_depths :: documents_zsigma_*` for the
!!      measured collapse.
!!
!!   3. **`zstar_h_min > H_VANISHED` under an INERT-role family is an
!!      error**, not a warning: a filler above the D4 skip/merge marker is
!!      a filler every h-gating kernel treats as LIVE while the coordinate
!!      still treats it as throwaway.
!!
!! Each refusal is exercised through `validate_config`'s `ierr` form (the
!! `error stop` form is untestable), the same harness
!! `test_ocean_cavity_draft` uses.
module test_ocean_vcoord_hygiene
   use rdb_constants, only: wp, VCOORD_ZSTAR_FULL
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   use rdb_config, only: config_t, read_config_from_string, validate_config
   use rdb_ocean_engine, only: ocean_engine_t, engine_setup, engine_teardown
   use testdrive, only: error_type, check, new_unittest, unittest_type
   implicit none
   private

   public :: collect_ocean_vcoord_hygiene_tests

   real(wp), parameter :: H_SURF_PROBE = 137.0_wp
      !! Deliberately unlike BOTH the `&vcoord_nml` default (0.0) and the
      !! slot default (5.0), so the assertion cannot pass by accident.
   real(wp), parameter :: H_MIN_PROBE = 1.25e-4_wp
      !! Likewise unlike the shared 1.0e-4 default — and still at or below
      !! `H_VANISHED = 1.5e-4`, so it stays inside the inert-filler
      !! contract and does not trip the refusal tested below.

contains

   subroutine collect_ocean_vcoord_hygiene_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("zstar_knobs_reach_the_slot_before_the_seed", test_knob_plumbing), &
                  new_unittest("zsigma_is_refused", test_zsigma_refused), &
                  new_unittest("zstar_full_is_accepted", test_zstar_full_accepted), &
                  new_unittest("h_min_above_h_vanished_is_refused", test_h_min_refused), &
                  new_unittest("h_min_on_the_marker_is_accepted", test_h_min_on_marker) &
                  ]
   end subroutine collect_ocean_vcoord_hygiene_tests

   function nml_with_vcoord(vcoord_body) result(nml)
      !! A minimal, in-envelope single-rank ocean namelist whose
      !! `&vcoord_nml` body the caller supplies.  Flat 1000 m bed, 4
      !! layers, split solver.  Not `pure`: builds a deferred-length
      !! result from concatenation, which is fine, but it is a test
      !! fixture and kept plain for readability.
      character(len=*), intent(in) :: vcoord_body
      character(len=:), allocatable :: nml
      nml = "&sim_nml sim_type = 'ocean' /"//new_line("a")// &
            "&grid_nml nx = 8, ny = 6, nghost = 2, dx = 1000.0, dy = 1000.0 /"// &
            new_line("a")// &
            "&nonhydrostatic_nml nz_layers = 4 /"//new_line("a")// &
            "&time_nml t_end = 3600.0, dt_fixed = 60.0 /"//new_line("a")// &
            "&ocean_topo_nml max_depth = 1000.0 /"//new_line("a")// &
            "&ocean_bt_nml auto_n_inner = .false., n_inner = 8 /"//new_line("a")// &
            "&vcoord_nml "//vcoord_body//" /"//new_line("a")// &
            "&ocean_diag_nml enabled = .false. /"//new_line("a")// &
            "&output_nml output_to_file = .false. /"//new_line("a")
   end function nml_with_vcoord

   subroutine expect_config(error, nml, want_ok, what)
      !! Run `read_config_from_string` + `validate_config` in their
      !! `ierr` forms and assert the verdict.  Not `pure`: `check`
      !! allocates `error`.
      type(error_type), allocatable, intent(inout) :: error
      character(len=*), intent(in) :: nml, what
      logical, intent(in) :: want_ok
      type(config_t) :: cfg
      integer :: ierr
      call read_config_from_string(nml, cfg, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) then
         call check(error, .false., what//": the namelist itself must parse")
         return
      end if
      call validate_config(cfg, ierr)
      if (want_ok) then
         call check(error, ierr == OCEAN_STATUS_OK, what//" must be ACCEPTED")
      else
         call check(error, ierr /= OCEAN_STATUS_OK, what//" must be REFUSED")
      end if
   end subroutine expect_config

   subroutine test_knob_plumbing(error)
      !! Drive the PRODUCTION configure path — `engine_setup`, the same
      !! entry `driver_run_ocean` and `rdb_ocean_create_from_string` both
      !! call — with non-default `zstar_h_surf_target` / `zstar_h_min`,
      !! then read the values back off `state%vcoord`.
      !!
      !! `z_ref(:, :, 1)` is checked too, because the copy being present
      !! is necessary but not sufficient: it must happen BEFORE
      !! `ocean_state_seed_from_cfg` runs `build_zref_full`, or the
      !! reference table would be built from the stale default.  With
      !! `h_surf_target = 137 m` on a 1000 m flat bed and the auto
      !! `n_surf = max(1, nz/3) = 1`, the first table interface must sit
      !! at exactly 137 m.
      type(error_type), allocatable, intent(out) :: error
      type(ocean_engine_t) :: engine
      type(config_t) :: cfg
      character(len=:), allocatable :: nml
      integer :: ierr
      checks: block
         nml = nml_with_vcoord("vcoord_type = 'zstar_full', zstar_h_surf_target = "// &
                               "137.0, zstar_h_min = 1.25e-4")
         call read_config_from_string(nml, cfg, ierr=ierr)
         call check(error, ierr == OCEAN_STATUS_OK, "the probe namelist must parse")
         if (allocated(error)) exit checks
         call validate_config(cfg, ierr)
         call check(error, ierr == OCEAN_STATUS_OK, "the probe namelist must validate")
         if (allocated(error)) exit checks
         call engine_setup(engine, cfg, ierr)
         call check(error, ierr == OCEAN_STATUS_OK, "engine_setup must succeed")
         if (allocated(error)) exit checks

         call check(error, engine%state%vcoord%coord_type == VCOORD_ZSTAR_FULL, &
                    "vcoord_type must reach the slot")
         if (allocated(error)) exit checks
         call check(error, abs(engine%state%vcoord%zstar_h_surf_target - H_SURF_PROBE) &
                    < 1.0e-12_wp, &
                    "&vcoord_nml zstar_h_surf_target must reach vcoord%zstar_h_surf_target")
         if (allocated(error)) exit checks
         call check(error, abs(engine%state%vcoord%zstar_h_min - H_MIN_PROBE) < 1.0e-12_wp, &
                    "&vcoord_nml zstar_h_min must reach vcoord%zstar_h_min")
         if (allocated(error)) exit checks
         ! The ordering half: the reference table must already carry the
         ! namelist value, i.e. the copy preceded the seed.
         call check(error, abs(engine%state%vcoord%z_ref(3, 3, 1) - H_SURF_PROBE) < 1.0e-9_wp, &
                    "build_zref_full must see the namelist zstar_h_surf_target "// &
                    "(the copy has to precede the IC seed)")
      end block checks
      call engine_teardown(engine)
   end subroutine test_knob_plumbing

   subroutine test_zsigma_refused(error)
      !! `zsigma` runs against a DIMENSIONLESS `z_ref_global` and collapses
      !! the whole column into the bed layer while keeping `Sum = H + eta`
      !! exact.  Refused until the table is filled in metres.
      type(error_type), allocatable, intent(out) :: error
      call expect_config(error, nml_with_vcoord("vcoord_type = 'zsigma'"), .false., &
                         "vcoord_type = 'zsigma'")
      if (allocated(error)) return
      ! `parse_vcoord_type` also accepts the aliases `z-sigma` / `z_sigma`,
      ! but `&vcoord_nml`'s schema enum lists only the `zsigma` spelling,
      ! so the alias never reaches `validate_config` on the namelist path
      ! — the reader rejects it first.  Nothing to assert here beyond the
      ! canonical spelling; if the enum ever grows the aliases, add them.
      ! ZSTAR_SIGMA consumes the same table FRACTIONALLY and is unaffected
      ! by its units, so it must stay accepted — 21 shipped namelists use it.
      call expect_config(error, nml_with_vcoord("vcoord_type = 'zstar_sigma'"), .true., &
                         "vcoord_type = 'zstar_sigma' (fractional, unaffected)")
   end subroutine test_zsigma_refused

   subroutine test_zstar_full_accepted(error)
      !! The regression net for the refusal above: every other geometric
      !! family still configures.
      type(error_type), allocatable, intent(out) :: error
      call expect_config(error, nml_with_vcoord("vcoord_type = 'sigma'"), .true., "sigma")
      if (allocated(error)) return
      call expect_config(error, nml_with_vcoord("vcoord_type = 'zstar'"), .true., "zstar")
      if (allocated(error)) return
      call expect_config(error, nml_with_vcoord("vcoord_type = 'zstar_full'"), .true., &
                         "zstar_full")
      if (allocated(error)) return
      call expect_config(error, nml_with_vcoord("vcoord_type = 'z_fixed'"), .true., "z_fixed")
      if (allocated(error)) return
      call expect_config(error, nml_with_vcoord("vcoord_type = 'lagrangian'"), .true., &
                         "lagrangian")
   end subroutine test_zstar_full_accepted

   subroutine test_h_min_refused(error)
      !! On the INERT-role families the filler must stay at or below
      !! `H_VANISHED = 1.5e-4`; above it the filler is dynamically live
      !! (EOS / PGF / remap-drain / vdiff) while the coordinate still
      !! treats it as throwaway.  Previously a warning — promoted to an
      !! error because a top-side filler that is not skipped participates
      !! in the pressure gradient, the melt sampler and the budgets.
      type(error_type), allocatable, intent(out) :: error
      call expect_config(error, nml_with_vcoord("vcoord_type = 'zstar_full', "// &
                                                "zstar_h_min = 1.0e-3"), .false., &
                         "zstar_full with zstar_h_min = 1e-3 (> H_VANISHED)")
      if (allocated(error)) return
      call expect_config(error, nml_with_vcoord("vcoord_type = 'z_fixed', "// &
                                                "zstar_h_min = 1.0e-3"), .false., &
                         "z_fixed with zstar_h_min = 1e-3 (> H_VANISHED)")
      if (allocated(error)) return
      ! A non-positive floor was already an error and stays one.
      call expect_config(error, nml_with_vcoord("vcoord_type = 'zstar_full', "// &
                                                "zstar_h_min = 0.0"), .false., &
                         "zstar_full with zstar_h_min = 0")
      if (allocated(error)) return
      ! The DENSITY families use the knob under the opposite contract
      ! (`max(zstar_h_min, 2*H_VANISHED)` keep-alive), so a large value
      ! there is meaningful rather than wrong and must stay accepted.
      call expect_config(error, nml_with_vcoord("vcoord_type = 'rho', "// &
                                                "zstar_h_min = 1.0e-3"), .true., &
                         "rho with zstar_h_min = 1e-3 (keep-alive role)")
   end subroutine test_h_min_refused

   subroutine test_h_min_on_marker(error)
      !! The boundary that five shipped namelists sit on: `zstar_h_min`
      !! EXACTLY `H_VANISHED`.  Legal — every downstream vanish test is a
      !! strict `> H_VANISHED`, so a layer on the marker reads as
      !! vanished — and the promotion above must not have moved the
      !! comparison to `>=`, which would refuse the canonical
      !! double-gyre reference.
      type(error_type), allocatable, intent(out) :: error
      call expect_config(error, nml_with_vcoord("vcoord_type = 'zstar_full', "// &
                                                "zstar_h_surf_target = 1000.0, "// &
                                                "zstar_h_min = 1.5e-4"), .true., &
                         "zstar_full with zstar_h_min exactly on H_VANISHED")
   end subroutine test_h_min_on_marker

end module test_ocean_vcoord_hygiene
