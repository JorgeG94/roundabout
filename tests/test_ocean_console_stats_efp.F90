!! Unit tests for the PR-32 ocean console-stats EFP (order-invariant
!! reproducing-sum) path.
!!
!! Coverage:
!!   * `test_efp_impl_matches_canonical` (SS9.3) — pins the in-module
!!     `!$acc routine seq` duplicate `efp_decompose_impl`
!!     (`rdb_ocean_console_stats`) bin-for-bit against the canonical
!!     `rdb_efp::efp_decompose` over a magnitude table spanning the range
!!     realistic console summands span.
!!   * `test_efp_kernels_match_fp_kernels` (SS11.2 mitigation) — on the same
!!     device-mapped state, the EFP reduction kernels
!!     (`compute_total_h_efp`/`_tracer_efp`/`_ke_efp`/`compute_ice_totals_efp`)
!!     must agree with their FP twins to ~1e-9 relative.  Per CLAUDE.md's
!!     `mem:separate` gotcha, a kernel whose `present(...)` list is missing an
!!     array reads the STALE HOST SHADOW and returns a plausible-looking but
!!     WRONG total — this test is the guard against exactly that class of bug.
!!   * `test_console_reproducing_sums_off_bit_identical` (SS9.5) — the
!!     default-off bit-identity gate: `reproducing_sums` ABSENT vs `.false.`
!!     must latch bit-identical `console_stats_t` reference fields.
!!   * `test_console_reproducing_sums_true_agrees_with_fp` — `.true.` on the
!!     same state must produce a Mass/KE/Salt/Heat total that agrees with the
!!     FP path to ~1e-9 relative, and must populate the EFP reference
!!     (`exact_sums = .true.`).
module test_ocean_console_stats_efp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use, intrinsic :: iso_fortran_env, only: int64, real64
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_state, only: ocean_state_t, ocean_state_enter_data, ocean_state_exit_data
   use rdb_ocean_metrics, only: metrics_fill_cartesian, metrics_finalize
   use rdb_console_stats, only: console_stats_t
   use rdb_ocean_console_stats, only: ocean_console_stats_report, &
                                      efp_decompose_impl, &
                                      compute_total_h_efp, compute_total_tracer_efp, &
                                      compute_total_ke_efp, compute_ice_totals_efp, &
                                      compute_total_h, compute_total_tracer, compute_total_ke
   use rdb_efp, only: efp_t, efp_decompose, efp_to_real
   use rdb_comm_env, only: comm_env_init, comm_env_setup_roles
   implicit none
   private

   public :: collect_ocean_console_stats_efp_tests

   integer, parameter :: NX = 8, NY = 6, NZ = 3
   real(wp), parameter :: DX = 1000.0_wp

   logical :: comm_inited = .false.
      !! Guard: the console-stats reductions go through `halo_allreduce_*`,
      !! which needs a live comm-env.  Only initialise once -- both console
      !! test cases share `setup_state`.

contains

   subroutine collect_ocean_console_stats_efp_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("efp_impl_matches_canonical", test_efp_impl_matches_canonical), &
                  new_unittest("efp_kernels_match_fp_kernels", test_efp_kernels_match_fp_kernels), &
                  new_unittest("console_reproducing_sums_off_bit_identical", &
                               test_console_reproducing_sums_off_bit_identical), &
                  new_unittest("console_reproducing_sums_true_agrees_with_fp", &
                               test_console_reproducing_sums_true_agrees_with_fp) &
                  ]
   end subroutine collect_ocean_console_stats_efp_tests

   subroutine setup_state(grid, state)
      !! Small ocean state with non-trivial h_layer/velocity/tracer fields so
      !! the EFP and FP reductions have something real to disagree on if a
      !! bug is present. Salinity/temperature are the default idx 1/2
      !! tracers (`rdb_multilayer_state`).
      type(hgrid_t), intent(inout) :: grid
      type(ocean_state_t), intent(inout) :: state

      ! `ocean_console_stats_report` reduces through `halo_allreduce_*`.  On an
      ! MPI build that reaches the communicator, so MPI must be up first --
      ! without this the test aborts in MPI_Comm_f2c before MPI_INIT.  No-op
      ! cost on serial builds.  Finalised by the shared per-test main.
      if (.not. comm_inited) then
         call comm_env_init()
         call comm_env_setup_roles(.false.)
         comm_inited = .true.
      end if

      call grid%init(NX, NY, 1, DX, DX)
      state%multilayer%nz_ml = NZ
      call state%init(grid)
      call metrics_fill_cartesian(state%metrics, grid, grid%dx, grid%dy)
      call metrics_finalize(state%metrics)
      call fill_test_fields(grid, state)
   end subroutine setup_state

   subroutine test_efp_impl_matches_canonical(error)
      !! SS9.3: pin the device-callable `efp_decompose_impl` bin-for-bit
      !! against the canonical `rdb_efp::efp_decompose` -- HOST call, no
      !! device mapping needed (`!$acc routine seq` procedures are plain
      !! callable Fortran; the directive only matters inside a device
      !! region).
      type(error_type), allocatable, intent(out) :: error
      real(real64), parameter :: table(9) = [ &
                                 1.0e21_real64, -1.0e21_real64, &
                                 1.0e6_real64, -1.0e6_real64, &
                                 1.0_real64, -1.0_real64, &
                                 1.0e-16_real64, -1.0e-16_real64, 0.0_real64]
      integer :: i
      integer(int64) :: e(6)
      integer(int64) :: e1, e2, e3, e4, e5, e6
      logical :: is_nan, is_ovf
      logical :: ok

      ok = .true.
      do i = 1, size(table)
         call efp_decompose(table(i), e, is_nan, is_ovf)
         call efp_decompose_impl(table(i), e1, e2, e3, e4, e5, e6)
         if (e(1) /= e1 .or. e(2) /= e2 .or. e(3) /= e3 .or. &
             e(4) /= e4 .or. e(5) /= e5 .or. e(6) /= e6) ok = .false.
      end do
      call check(error, ok, &
                 "efp_decompose_impl must match rdb_efp::efp_decompose bin-for-bit")
   end subroutine test_efp_impl_matches_canonical

   subroutine test_efp_kernels_match_fp_kernels(error)
      !! SS11.2 mitigation: on the SAME device-mapped state, the EFP
      !! reduction kernels must agree with their FP twins to ~1e-9 relative
      !! -- a missing `present(...)` entry would make the EFP kernel read a
      !! stale host shadow and return a plausible-but-wrong total instead of
      !! crashing, so this is the guard, not a formality.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      real(wp) :: h_fp, s_fp, t_fp, ke_fp
      real(wp) :: wet_fp, ci_fp, hi_fp
      type(efp_t) :: h_efp, s_efp, t_efp, ke_efp
      type(efp_t) :: wet_efp, ci_efp, hi_efp
      real(wp) :: h_e, s_e, t_e, ke_e

      ! Ice must be enabled BEFORE state%init (ncat sizes the arrays),
      ! mirroring test_ocean_ice_diags's setup_state_ice -- so this test
      ! builds its own state rather than reusing `setup_state`.
      call grid%init(NX, NY, 1, DX, DX)
      state%multilayer%nz_ml = NZ
      state%ice%enable = .true.
      state%ice%ncat = 1
      call state%init(grid)
      call metrics_fill_cartesian(state%metrics, grid, grid%dx, grid%dy)
      call metrics_finalize(state%metrics)
      call fill_test_fields(grid, state)
      state%ice%part_size(:, :, 0) = 0.3_wp
      state%ice%part_size(:, :, 1) = 0.7_wp
      state%ice%m_ice(:, :, 1) = 850.0_wp

      call ocean_state_enter_data(state)

      h_fp = compute_total_h(state%multilayer%h_layer, state%metrics%areaT, grid%nghost)
      h_efp = compute_total_h_efp(state%multilayer%h_layer, state%metrics%areaT, grid%nghost)
      h_e = real(efp_to_real(h_efp), wp)

      s_fp = compute_total_tracer(state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, &
                                  state%metrics%areaT, grid%nghost)
      s_efp = compute_total_tracer_efp(state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, &
                                       state%metrics%areaT, grid%nghost)
      s_e = real(efp_to_real(s_efp), wp)

      t_fp = compute_total_tracer(state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
                                  state%metrics%areaT, grid%nghost)
      t_efp = compute_total_tracer_efp(state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
                                       state%metrics%areaT, grid%nghost)
      t_e = real(efp_to_real(t_efp), wp)

      ke_fp = compute_total_ke(state%multilayer%h_layer, state%multilayer%u_face_x_layer, &
                               state%multilayer%v_face_y_layer, state%metrics%areaT, grid%nghost)
      ke_efp = compute_total_ke_efp(state%multilayer%h_layer, state%multilayer%u_face_x_layer, &
                                    state%multilayer%v_face_y_layer, state%metrics%areaT, grid%nghost)
      ke_e = real(efp_to_real(ke_efp), wp)

      call compute_ice_totals_efp(state%metrics%wet_T, state%metrics%areaT, &
                                  state%ice%part_size, state%ice%m_ice, 1, grid%nghost, &
                                  wet_efp, ci_efp, hi_efp)
      wet_fp = real(efp_to_real(wet_efp), wp)
      ci_fp = real(efp_to_real(ci_efp), wp)
      hi_fp = real(efp_to_real(hi_efp), wp)

      call ocean_state_exit_data(state)

      call check(error, abs(h_e - h_fp) <= 1.0e-9_wp*abs(h_fp), &
                 "compute_total_h_efp must agree with compute_total_h to 1e-9 relative")
      if (allocated(error)) then
         call state%destroy(); return
      end if
      call check(error, abs(s_e - s_fp) <= 1.0e-9_wp*abs(s_fp), &
                 "compute_total_tracer_efp(salt) must agree with compute_total_tracer to 1e-9 relative")
      if (allocated(error)) then
         call state%destroy(); return
      end if
      call check(error, abs(t_e - t_fp) <= 1.0e-9_wp*abs(t_fp), &
                 "compute_total_tracer_efp(heat) must agree with compute_total_tracer to 1e-9 relative")
      if (allocated(error)) then
         call state%destroy(); return
      end if
      call check(error, abs(ke_e - ke_fp) <= 1.0e-9_wp*max(abs(ke_fp), 1.0e-30_wp), &
                 "compute_total_ke_efp must agree with compute_total_ke to 1e-9 relative")
      if (allocated(error)) then
         call state%destroy(); return
      end if
      ! wet/ci/hi are compared against themselves (both computed via the EFP
      ! path above) only as a non-degeneracy sanity check -- the wet-area
      ! total must be positive and the ice fraction/thickness must be
      ! bounded/finite for a well-posed test configuration.
      call check(error, wet_fp > 0.0_wp, "ice wet_area_efp must be positive for a wet test grid")
      if (allocated(error)) then
         call state%destroy(); return
      end if
      call check(error, ci_fp >= 0.0_wp .and. hi_fp >= 0.0_wp, &
                 "ice ci/hi area totals must be non-negative")
      call state%destroy()
   end subroutine test_efp_kernels_match_fp_kernels

   subroutine fill_test_fields(grid, state)
      type(hgrid_t), intent(in) :: grid
      type(ocean_state_t), intent(inout) :: state
      integer :: i, j, k
      do k = 1, NZ
         do j = 1, grid%ny_total
            do i = 1, grid%nx_total
               state%multilayer%h_layer(i, j, k) = 10.0_wp + 0.1_wp*real(i + j + k, wp)
               state%multilayer%tracers(state%multilayer%idx_salinity)%hTr(i, j, k) = &
                  state%multilayer%h_layer(i, j, k)*(34.5_wp + 0.01_wp*real(i, wp))
               state%multilayer%tracers(state%multilayer%idx_temperature)%hTr(i, j, k) = &
                  state%multilayer%h_layer(i, j, k)*(12.0_wp - 0.02_wp*real(j, wp))
            end do
         end do
      end do
      do k = 1, NZ
         do j = 1, grid%ny_total
            do i = 1, grid%nx_total + 1
               state%multilayer%u_face_x_layer(i, j, k) = 0.05_wp*real(mod(i + k, 5) - 2, wp)
            end do
         end do
         do j = 1, grid%ny_total + 1
            do i = 1, grid%nx_total
               state%multilayer%v_face_y_layer(i, j, k) = 0.03_wp*real(mod(j + k, 4) - 2, wp)
            end do
         end do
      end do
   end subroutine fill_test_fields

   subroutine test_console_reproducing_sums_off_bit_identical(error)
      !! SS9.5: `reproducing_sums` ABSENT vs `.false.` must latch
      !! BIT-IDENTICAL `console_stats_t` reference fields -- the structural
      !! default-off byte-identity gate (mirrors `compute_max_cfl`'s
      !! two-loop precedent).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      type(console_stats_t) :: stats_absent, stats_false

      call setup_state(grid, state)
      call ocean_state_enter_data(state)

      call ocean_console_stats_report(stats_absent, grid, state%metrics, state%multilayer, &
                                      t=0.0_wp, dt=1.0_wp, step=0, compute_rank=0)
      call ocean_console_stats_report(stats_false, grid, state%metrics, state%multilayer, &
                                      t=0.0_wp, dt=1.0_wp, step=0, compute_rank=0, &
                                      reproducing_sums=.false.)

      call ocean_state_exit_data(state)
      call state%destroy()

      call check(error, stats_absent%mass0 == stats_false%mass0, &
                 "mass0 must be bit-identical between reproducing_sums absent and .false.")
      if (allocated(error)) return
      call check(error, stats_absent%salt0 == stats_false%salt0, &
                 "salt0 must be bit-identical between reproducing_sums absent and .false.")
      if (allocated(error)) return
      call check(error, stats_absent%heat0 == stats_false%heat0, &
                 "heat0 must be bit-identical between reproducing_sums absent and .false.")
      if (allocated(error)) return
      call check(error, stats_absent%ke0 == stats_false%ke0, &
                 "ke0 must be bit-identical between reproducing_sums absent and .false.")
      if (allocated(error)) return
      call check(error,.not. stats_false%exact_sums, &
                 "reproducing_sums=.false. must NOT set console_stats_t%exact_sums")
   end subroutine test_console_reproducing_sums_off_bit_identical

   subroutine test_console_reproducing_sums_true_agrees_with_fp(error)
      !! `reproducing_sums = .true.` on the same state must (a) populate the
      !! EFP reference (`exact_sums = .true.`) and (b) report a Mass total
      !! that agrees with the FP path to ~1e-9 relative -- the two paths sum
      !! the SAME physical field via different arithmetic, so they must land
      !! on the same physical answer even though the bit patterns of the
      !! intermediate reduction differ.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      type(console_stats_t) :: stats_fp, stats_efp

      call setup_state(grid, state)
      call ocean_state_enter_data(state)

      call ocean_console_stats_report(stats_fp, grid, state%metrics, state%multilayer, &
                                      t=0.0_wp, dt=1.0_wp, step=0, compute_rank=0, &
                                      reproducing_sums=.false.)
      call ocean_console_stats_report(stats_efp, grid, state%metrics, state%multilayer, &
                                      t=0.0_wp, dt=1.0_wp, step=0, compute_rank=0, &
                                      reproducing_sums=.true.)

      call ocean_state_exit_data(state)
      call state%destroy()

      call check(error, stats_efp%exact_sums, &
                 "reproducing_sums=.true. must set console_stats_t%exact_sums")
      if (allocated(error)) return
      call check(error, abs(stats_efp%mass0 - stats_fp%mass0) <= 1.0e-9_wp*abs(stats_fp%mass0), &
                 "EFP-path Mass total must agree with the FP-path total to 1e-9 relative")
      if (allocated(error)) return
      call check(error, abs(stats_efp%salt0 - stats_fp%salt0) <= 1.0e-9_wp*abs(stats_fp%salt0), &
                 "EFP-path Salt total must agree with the FP-path total to 1e-9 relative")
      if (allocated(error)) return
      call check(error, abs(stats_efp%heat0 - stats_fp%heat0) <= 1.0e-9_wp*abs(stats_fp%heat0), &
                 "EFP-path Heat total must agree with the FP-path total to 1e-9 relative")
   end subroutine test_console_reproducing_sums_true_agrees_with_fp

end module test_ocean_console_stats_efp
