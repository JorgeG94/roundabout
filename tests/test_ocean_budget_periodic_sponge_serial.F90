!! SERIAL (real single-rank, no MPI) regression for the periodic-seam /
!! sponge-seam-ghost bug: a periodic axis relaxed by a sponge on the
!! orthogonal edge must complete on ONE rank with a finite, bounded
!! closed-budget `mass_out`.  It is a completes-and-finite smoke test
!! only: it does NOT fail with the fix reverted.  The regression gate is
!! the MPI-enabled twin (`tests/mpi/test_ocean_decomp_bitid_mpi.F90 ::
!! periodic_sponge`, `check_periodic_sponge_serial_out`) -- see "Why this
!! file is a smoke test, not a tight regression" below.
!!
!! ## The bug
!!
!! `rdb_ocean_dyn.F90`'s post-sponge seam-ghost refresh used to be gated on
!! `ocean_halo_is_decomposed_x() .or. ocean_halo_is_decomposed_y()`.  Both
!! sponges (legacy band + map-driven) relax a tile's PHYSICAL cells only;
!! the periodic-seam ghost copies of those same physical cells keep the
!! UN-relaxed value until something refreshes them, and the corrector's
!! advection reads them before the stage-end exchange.  On a real MPI
!! decomposition this refresh is unconditionally needed and the old gate
!! provided it (any valid `px x py` factorization of `nprocs >= 2` has
!! `px > 1 .or. py > 1`, so the gate was always true there).  On a
!! genuinely single-rank run with a periodic axis, the gate was FALSE -- no
!! MPI seam exists, so `ocean_halo_is_decomposed_x/y()` are both false --
!! and the refresh never ran, even though the periodic axis still needs its
!! LOCAL wrap redone after the sponge touched only the physical columns.
!! Measured on the Southern Ocean 1-degree cut: console `Salt out` 1.6e8 on
!! 1 rank vs -81 on 2 ranks, while `Error` (the closed-budget residual)
!! stayed at round-off on BOTH -- the stale-ghost flux is self-consistently
!! integrated into both the real state and the `out` accumulator, so the
!! closed-budget identity still balances; only the raw `out` term itself is
!! orders of magnitude off.
!!
!! Fix: the refresh in `rdb_ocean_dyn.F90` (~line 4948, guarded by
!! `sponge_seam`) now runs D0-unconditionally -- `ocean_halo_face_x`/`_y`
!! and `refresh_tracer_ghosts` already no-op on a single-rank non-periodic
!! axis, locally re-wrap a single-rank periodic one, and exchange messages
!! when decomposed, so there is nothing left for a decomposed-only gate to
!! usefully skip.
!!
!! ## Why this file is a smoke test, not a tight regression
!!
!! At the Southern Ocean's scale the bug is unmistakable (orders of
!! magnitude).  On THIS file's toy grid (26x18, `N_STEPS = 48`) it is not:
!! `mass_out` here is dominated by ordinary discretization noise that is
!! itself large enough, and platform/compiler-sensitive enough, to swamp
!! the bug's own signature.  Measured directly while calibrating this file
!! (gfortran default flags vs nvfortran `-gpu=cc70`, otherwise IDENTICAL
!! namelist and step count, EFP-reduced `mass_out` via
!! `ocean_console_stats_report`'s `budget_out`):
!!
!! | build              | fixed        | broken       |
!! |--------------------|-------------:|-------------:|
!! | gfortran (host)    | -1.445E-04   | -1.038E-04   |
!! | nvfortran (GPU)    | -7.900E-06   |  1.103E-04   |
!!
!! Neither the MAGNITUDE nor even the SIGN of "fixed vs broken" is stable
!! across toolchains at this scale, so no fixed threshold here can both (a)
!! pass on every toolchain this suite runs on and (b) reliably fail without
!! the fix -- a tight bound would be flaky CI, not a regression gate.  The
!! MPI-enabled companion case has no threshold to calibrate: it compares
!! the single-rank reference's closed-budget totals BIT FOR BIT against
!! every decomposed `px x py` run of the same problem (EFP reproducing
!! sums make that an exact identity), so the stale-ghost reference
!! differs from the decomposed runs whenever the fix is reverted.  It can
!! only do that with something to compare against, i.e. on its 2- and
!! 4-rank ctest legs (`rdb_test_ocean_decomp_bitid_mpi`, `..._4rank`;
!! `RDB_ENABLE_MPI=ON` only).  This file exists so the
!! `RDB_ENABLE_MPI=OFF` build (the default local/CI build per CLAUDE.md)
!! still runs the configuration.
!!
!! `periodic_sponge_seam_out_is_finite_and_bounded` therefore only asserts
!! the run completes and `mass_out` stays finite and under a generous,
!! toolchain-portable ceiling (`MASS_OUT_CEILING`, comfortably above every
!! value in the table above).  It catches a NaN or a crash, not the seam
!! bug itself.
module test_ocean_budget_periodic_sponge_serial
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rdb_constants, only: wp
   use rdb_config, only: config_t, read_config_from_string, validate_config
   use rdb_ocean_engine, only: ocean_engine_t, engine_setup, engine_enter_data, &
                               engine_step, engine_step_finalize, &
                               engine_exit_data, engine_teardown
   use rdb_console_stats, only: console_stats_t, conservation_budget_t
   use rdb_ocean_console_stats, only: ocean_console_stats_report, ocean_budget_stage_weight
   use rdb_ocean_dyn, only: SPLIT_SCHEME_PRED_CORR
   use rdb_ocean_status, only: OCEAN_STATUS_OK
   use rdb_comm_env, only: comm_env_init, comm_env_setup_roles
   implicit none
   private

   public :: collect_ocean_budget_periodic_sponge_serial_tests

   integer, parameter :: NX_G = 26
   integer, parameter :: NY_G = 18
   integer, parameter :: N_STEPS = 48
   real(wp), parameter :: DT = 900.0_wp
   real(wp), parameter :: MASS_OUT_CEILING = 1.0e-2_wp
      !! Generous, toolchain-portable ceiling -- see the module docstring's
      !! measured noise-floor table (both builds' fixed AND broken values
      !! sit two-plus orders of magnitude below this).

   logical :: comm_inited = .false.

contains

   subroutine collect_ocean_budget_periodic_sponge_serial_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("periodic_sponge_seam_out_is_finite_and_bounded", &
                               test_periodic_sponge_seam_out) &
                  ]
   end subroutine collect_ocean_budget_periodic_sponge_serial_tests

   subroutine ensure_comm()
      !! Same reasoning as `test_ocean_cavity_freshwater::ensure_comm`: on
      !! an `RDB_ENABLE_MPI=ON` build the first collective in `engine_setup`
      !! / `ocean_console_stats_report` hits `MPI_Comm_f2c` before
      !! `MPI_Init` without this. No-op-equivalent on the single-rank
      !! backend. Finalised by the shared per-test main.
      if (.not. comm_inited) then
         call comm_env_init()
         call comm_env_setup_roles(.false.)
         comm_inited = .true.
      end if
   end subroutine ensure_comm

   function channel_nml() result(nml)
      !! Periodic west/east, sponge-relaxed north edge, closed south --
      !! same shape as the MPI companion case `periodic_sponge`
      !! (`tests/mpi/test_ocean_decomp_bitid_mpi.F90`).
      character(len=:), allocatable :: nml
      character(len=*), parameter :: NL = new_line("a")
      character(len=16) :: snx, sny

      write (snx, '(i0)') NX_G
      write (sny, '(i0)') NY_G
      nml = "&sim_nml sim_type = 'ocean' /"//NL// &
            "&time_nml t_end = 86400.0, dt_fixed = 900.0 /"//NL// &
            "&nonhydrostatic_nml nz_layers = 4 /"//NL// &
            "&tracer_nml initial_temperature = 12.0, initial_salinity = 35.0, "// &
            "T_init_surface = 20.0, T_init_bottom = 4.0 /"//NL// &
            "&ocean_bt_nml auto_n_inner = .true., split_scheme = 'pred_corr' /"//NL// &
            "&ocean_hvisc_nml nu_h = 200.0, lateral_closure = 'smagorinsky', "// &
            "smag_ah = .true. /"//NL// &
            "&ocean_diag_nml enabled = .false. /"//NL// &
            "&grid_nml nx = "//trim(snx)//", ny = "//trim(sny)//", nghost = 3, "// &
            "dx = 20000.0, dy = 20000.0 /"//NL// &
            "&physics_nml coriolis_f = 1.0e-4, wind_stress_x = 0.05 /"//NL// &
            "&vcoord_nml vcoord_type = 'sigma' /"//NL// &
            "&ocean_topo_nml topo_config = 'seamount', max_depth = 2000.0, "// &
            "edge_depth = 1500.0, slope_scale = 60000.0 /"//NL// &
            "&ocean_bc_nml west = 'periodic', east = 'periodic', south = 'wall', "// &
            "north = 'sponge', sponge_width = 3, sponge_strength = 1.0e-4, "// &
            "sponge_relax_tracers = .true. /"//NL// &
            "&output_nml output_to_file = .false. /"//NL
   end function channel_nml

   subroutine run_and_get_mass_out(mass_out, ok)
      !! Configure a real single-rank engine (no compute_rank/compute_size
      !! override -- production default), step N_STEPS, and read the
      !! EFP-reduced `mass_out` exactly as the console does.
      real(wp), intent(out) :: mass_out
      logical, intent(out) :: ok
      type(ocean_engine_t) :: engine
      type(config_t) :: cfg
      type(console_stats_t) :: cstats
      type(conservation_budget_t) :: bud
      integer :: ierr, n
      real(wp) :: t

      ok = .false.
      mass_out = 0.0_wp
      call ensure_comm()
      call read_config_from_string(channel_nml(), cfg, ierr=ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call validate_config(cfg, ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call engine_setup(engine, cfg, ierr)
      if (ierr /= OCEAN_STATUS_OK) return
      call engine_enter_data(engine, cfg)

      t = 0.0_wp
      ierr = OCEAN_STATUS_OK
      do n = 1, N_STEPS
         call engine_step(engine, DT, t, ierr=ierr)
         if (ierr /= OCEAN_STATUS_OK) exit
         call engine_step_finalize(engine, DT, t, ierr=ierr)
         if (ierr /= OCEAN_STATUS_OK) exit
         t = t + DT
      end do
      if (ierr == OCEAN_STATUS_OK) then
         call ocean_console_stats_report(cstats, engine%grid, engine%state%metrics, &
                                         engine%state%multilayer, t, DT, N_STEPS, &
                                         compute_rank=0, reproducing_sums=.true., &
                                         budget_stage_weight=ocean_budget_stage_weight( &
                                         engine%state%dyn%split_scheme == SPLIT_SCHEME_PRED_CORR), &
                                         budget_out=bud)
         mass_out = bud%mass_out
         ok = ieee_is_finite(mass_out)
      end if
      call engine_exit_data(engine)
      call engine_teardown(engine)
   end subroutine run_and_get_mass_out

   subroutine test_periodic_sponge_seam_out(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: mass_out
      logical :: ok

      call run_and_get_mass_out(mass_out, ok)
      call check(error, ok, "periodic-sponge run failed to complete, or mass_out is non-finite")
      if (allocated(error)) return

      call check(error, abs(mass_out) <= MASS_OUT_CEILING, &
                 "periodic-channel mass_out exceeded the toolchain-portable smoke ceiling -- "// &
                 "see tests/mpi/test_ocean_decomp_bitid_mpi.F90 :: periodic_sponge for the "// &
                 "tight regression")
   end subroutine test_periodic_sponge_seam_out

end module test_ocean_budget_periodic_sponge_serial
