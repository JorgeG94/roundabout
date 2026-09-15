!! Unit tests for the map-driven ocean sponge (PR-23,
!! `rdb_ocean_sponge::ocean_sponge_t` + `ocean_sponge_apply_maps`, plus the
!! `&ocean_sponge_nml damp_source="band"` map builder in
!! `rdb_ocean_setup::configure_ocean_sponge`).
!!
!! See `docs/plans/PLAN_PR23_real_sponge.md` §9 for the full test list and
!! the reasoning behind each. Cases:
!!   (1) `sponge_relaxes_tracer_to_3d_target_analytically` — the headline
!!       property missing from the legacy band (§2.2): a VERTICALLY VARYING
!!       reference is only expressible with a 3-D field.
!!   (2) `sponge_relaxes_momentum_to_u_ref_not_zero` — the legacy band always
!!       relaxes toward zero (§2.1); the map-driven path relaxes toward the
!!       configured `u_ref`.
!!   (3) `sponge_idamp_zero_is_exact_identity` — `Idamp = 0` IS the sponge
!!       mask; a cell outside it must be bit-for-bit untouched.
!!   (4) `sponge_band_builder_reproduces_cosine_ramp` — the analytic
!!       `damp_source="band"` filler reproduces the legacy kernel's cosine
!!       ramp at the exact cell/u-face/v-face offsets.
!!   (5) `sponge_band_corners_sum` — overlapping edge bands SUM their rates.
!!   (6) `sponge_source_closes_the_heat_budget` — the tracer relaxation's S/T
!!       source, mirrored into `ms%heat_budget_sponge`, closes the console
!!       heat budget via `ocean_heat_src_sum`.
!!   (7) `sponge_reference_is_the_ic_not_the_restart` — the reference
!!       snapshot freezes whatever state existed AT CALL TIME and is inert
!!       to any later mutation of `ms` — the property the driver's call-site
!!       ordering (BEFORE the restart read) relies on to make
!!       `target_source="ic"` mean the IC and not a resumed run's mid-run
!!       state (docs/plans/PLAN_PR23_real_sponge.md §13.1 item 4).
!!
!! Every kernel-level test maps `ms` and `sp` to the device (`!$acc enter
!! data` + bound `enter_data`, `!$acc update self` before host reads) per
!! CLAUDE.md's `mem:separate` contract (§6.6 of the plan) — a green
!! multicore run proves nothing about device data motion.
module test_ocean_sponge
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, RHO_WATER
   use rdb_grid, only: hgrid_t
   use rdb_config, only: config_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_state, only: ocean_state_t
   use rdb_ocean_boundary_types, only: OBC_SPONGE
   use rdb_ocean_setup, only: configure_ocean_sponge
   use rdb_ocean_sponge, only: ocean_sponge_t, ocean_sponge_apply_maps, &
                               ocean_sponge_snapshot_reference
   use rdb_ocean_console_stats, only: ocean_budget_src, ocean_heat_src_sum
   implicit none
   private

   public :: collect_ocean_sponge_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 4

contains

   subroutine collect_ocean_sponge_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("sponge_relaxes_tracer_to_3d_target_analytically", &
                               test_tracer_3d_target), &
                  new_unittest("sponge_relaxes_momentum_to_u_ref_not_zero", &
                               test_momentum_to_u_ref), &
                  new_unittest("sponge_idamp_zero_is_exact_identity", &
                               test_idamp_zero_identity), &
                  new_unittest("sponge_band_builder_reproduces_cosine_ramp", &
                               test_band_builder), &
                  new_unittest("sponge_band_corners_sum", &
                               test_band_corners_sum), &
                  new_unittest("sponge_source_closes_the_heat_budget", &
                               test_budget_closes), &
                  new_unittest("sponge_reference_is_the_ic_not_the_restart", &
                               test_reference_is_ic_not_restart) &
                  ]
   end subroutine collect_ocean_sponge_tests

   subroutine make_grid(grid, nx_phys, ny_phys, dx, dy)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      real(wp), intent(in) :: dx, dy
      call grid%init(nx_phys, ny_phys, NGHOST, dx, dy)
   end subroutine make_grid

   pure function area_sum(field, nx, ny, area_w) result(s)
      !! Host-side interior (ghost-excluded) area-weighted reduction —
      !! mirrors `compute_total_tracer`'s stencil on a uniform Cartesian
      !! grid, matching what `ocean_console_stats_report` feeds
      !! `ocean_heat_src_sum` / `ocean_budget_src`.
      real(wp), intent(in) :: field(:, :, :)
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: area_w
      real(wp) :: s
      s = sum(field(NGHOST + 1:nx - NGHOST, NGHOST + 1:ny - NGHOST, :))*area_w
   end function area_sum

   ! ------------------------------------------------------------------
   ! (1) 3-D tracer target
   ! ------------------------------------------------------------------
   subroutine test_tracer_3d_target(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_sponge_t) :: sp
      real(wp), parameter :: LAMBDA = 0.05_wp
      real(wp), parameter :: DT = 1.0_wp
      real(wp), parameter :: S_IC = 35.0_wp
      integer, parameter :: NSTEPS = 20
      integer, parameter :: NX_PHYS = 4, NY_PHYS = 4
      real(wp) :: ref_k, expected, max_rel_err, err
      integer :: i, j, k, step, i0, i1, j0, j1

      checks: block
         call make_grid(grid, NX_PHYS, NY_PHYS, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         ms%h_layer = 1.0_wp
         ms%tracers(ms%idx_salinity)%hTr = S_IC*ms%h_layer

         sp%enable = .true.
         call sp%init(grid, nz_ml=NZ, n_tracers=size(ms%tracers))
         sp%relax_tracers = .true.
         sp%relax_uv = .false.
         sp%idamp_h = LAMBDA
         do k = 1, NZ
            sp%ref_tracer(:, :, k, ms%idx_salinity) = 30.0_wp + real(k, wp)
         end do

         !$acc enter data copyin(ms, sp)
         call ms%enter_data()
         call sp%enter_data()
         do step = 1, NSTEPS
            call ocean_sponge_apply_maps(grid, sp, ms, DT)
         end do
         !$acc update self(ms%tracers(ms%idx_salinity)%hTr)
         call sp%exit_data()
         call ms%exit_data()
         !$acc exit data delete(ms, sp)

         i0 = grid%nghost + 1; i1 = grid%nghost + NX_PHYS
         j0 = grid%nghost + 1; j1 = grid%nghost + NY_PHYS
         max_rel_err = 0.0_wp
         do k = 1, NZ
            ref_k = 30.0_wp + real(k, wp)
            expected = ref_k + (S_IC - ref_k)*exp(-LAMBDA*real(NSTEPS, wp)*DT)
            do j = j0, j1
               do i = i0, i1
                  ! h_layer == 1 everywhere, so hTr IS the concentration.
                  err = abs(ms%tracers(ms%idx_salinity)%hTr(i, j, k) - expected)/abs(expected)
                  max_rel_err = max(max_rel_err, err)
               end do
            end do
         end do

         call check(error, max_rel_err < 1.0e-12_wp, &
                    "sponge tracer relaxation must match phi_ref + (phi0-phi_ref)*exp(-Idamp*t) "// &
                    "at every k (proves the reference is 3-D, not a scalar)")

      end block checks
      call sp%destroy()
      call ms%destroy()
   end subroutine test_tracer_3d_target

   ! ------------------------------------------------------------------
   ! (2) momentum relaxes toward u_ref, not zero
   ! ------------------------------------------------------------------
   subroutine test_momentum_to_u_ref(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_sponge_t) :: sp
      real(wp), parameter :: LAMBDA = 0.05_wp
      real(wp), parameter :: DT = 1.0_wp
      real(wp), parameter :: U0 = 0.5_wp
      real(wp), parameter :: U_REF = 0.3_wp
      integer, parameter :: NSTEPS = 20
      integer, parameter :: NX_PHYS = 4, NY_PHYS = 4
      real(wp) :: expected, computed
      integer :: step, i_probe, j_probe

      checks: block
         call make_grid(grid, NX_PHYS, NY_PHYS, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         ms%h_layer = 1.0_wp
         ms%u_face_x_layer = U0
         ms%v_face_y_layer = 0.0_wp

         sp%enable = .true.
         call sp%init(grid, nz_ml=NZ, n_tracers=size(ms%tracers))
         sp%relax_uv = .true.
         sp%relax_tracers = .false.
         sp%idamp_u = LAMBDA
         sp%u_ref = U_REF

         !$acc enter data copyin(ms, sp)
         call ms%enter_data()
         call sp%enter_data()
         do step = 1, NSTEPS
            call ocean_sponge_apply_maps(grid, sp, ms, DT)
         end do
         !$acc update self(ms%u_face_x_layer)
         call sp%exit_data()
         call ms%exit_data()
         !$acc exit data delete(ms, sp)

         i_probe = grid%nghost + 2
         j_probe = grid%nghost + 2
         computed = ms%u_face_x_layer(i_probe, j_probe, 1)
         expected = U_REF + (U0 - U_REF)*exp(-LAMBDA*real(NSTEPS, wp)*DT)

         call check(error, abs(computed - expected) < 1.0e-12_wp*U0, &
                    "sponge momentum relaxation must match U_REF + (U0-U_REF)*exp(-Idamp*t)")
         if (allocated(error)) exit checks

         ! Steady state -> u_ref, NOT zero: today's legacy kernel gives
         ! u -> 0 unconditionally (§2.1); the whole point of the map-driven
         ! path is that a non-zero target is reachable.
         call check(error, abs(U_REF) > 1.0e-6_wp, &
                    "sanity: U_REF must be nonzero for this test to be meaningful")

      end block checks
      call sp%destroy()
      call ms%destroy()
   end subroutine test_momentum_to_u_ref

   ! ------------------------------------------------------------------
   ! (3) Idamp = 0 is an exact (bit-for-bit) identity
   ! ------------------------------------------------------------------
   subroutine test_idamp_zero_identity(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_sponge_t) :: sp
      real(wp), parameter :: DT = 1.0_wp
      integer, parameter :: NSTEPS = 20
      integer, parameter :: NX_PHYS = 4, NY_PHYS = 4
      real(wp), allocatable :: u0(:, :, :), v0(:, :, :), s0(:, :, :), t0(:, :, :), h0(:, :, :)
      integer :: step

      checks: block
         call make_grid(grid, NX_PHYS, NY_PHYS, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         ms%h_layer = 3.0_wp
         ms%u_face_x_layer = 0.7_wp
         ms%v_face_y_layer = -0.4_wp
         ms%tracers(ms%idx_salinity)%hTr = 35.0_wp*ms%h_layer
         ms%tracers(ms%idx_temperature)%hTr = 12.0_wp*ms%h_layer

         sp%enable = .true.
         call sp%init(grid, nz_ml=NZ, n_tracers=size(ms%tracers))
         sp%relax_uv = .true.
         sp%relax_tracers = .true.
         ! idamp_h/u/v stay at their init-time zero fill — Idamp=0 IS the
         ! sponge mask (no explicit assignment needed).
         ! Arbitrary, deliberately "wrong" references — if the zero-Idamp
         ! guard is broken these would leak into the interior instantly.
         sp%ref_tracer = 0.0_wp
         sp%u_ref = 999.0_wp
         sp%v_ref = -999.0_wp

         allocate (u0, source=ms%u_face_x_layer)
         allocate (v0, source=ms%v_face_y_layer)
         allocate (s0, source=ms%tracers(ms%idx_salinity)%hTr)
         allocate (t0, source=ms%tracers(ms%idx_temperature)%hTr)
         allocate (h0, source=ms%h_layer)

         !$acc enter data copyin(ms, sp)
         call ms%enter_data()
         call sp%enter_data()
         do step = 1, NSTEPS
            call ocean_sponge_apply_maps(grid, sp, ms, DT)
         end do
         !$acc update self(ms%u_face_x_layer, ms%v_face_y_layer, ms%h_layer)
         !$acc update self(ms%tracers(ms%idx_salinity)%hTr)
         !$acc update self(ms%tracers(ms%idx_temperature)%hTr)
         !$acc update self(ms%salt_budget_sponge, ms%heat_budget_sponge)
         call sp%exit_data()
         call ms%exit_data()
         !$acc exit data delete(ms, sp)

         call check(error, all(ms%u_face_x_layer == u0), &
                    "u_face_x_layer must be bit-for-bit unchanged when idamp_u==0 everywhere")
         if (allocated(error)) exit checks
         call check(error, all(ms%v_face_y_layer == v0), &
                    "v_face_y_layer must be bit-for-bit unchanged when idamp_v==0 everywhere")
         if (allocated(error)) exit checks
         call check(error, all(ms%tracers(ms%idx_salinity)%hTr == s0), &
                    "salinity hTr must be bit-for-bit unchanged when idamp_h==0 everywhere")
         if (allocated(error)) exit checks
         call check(error, all(ms%tracers(ms%idx_temperature)%hTr == t0), &
                    "temperature hTr must be bit-for-bit unchanged when idamp_h==0 everywhere")
         if (allocated(error)) exit checks
         call check(error, all(ms%h_layer == h0), &
                    "h_layer must be untouched (v1 has no relax_h capability)")
         if (allocated(error)) exit checks
         call check(error, all(ms%salt_budget_sponge == 0.0_wp), &
                    "salt_budget_sponge must stay exactly 0 when idamp_h==0 everywhere")
         if (allocated(error)) exit checks
         call check(error, all(ms%heat_budget_sponge == 0.0_wp), &
                    "heat_budget_sponge must stay exactly 0 when idamp_h==0 everywhere")

      end block checks
      if (allocated(u0)) deallocate (u0)
      if (allocated(v0)) deallocate (v0)
      if (allocated(s0)) deallocate (s0)
      if (allocated(t0)) deallocate (t0)
      if (allocated(h0)) deallocate (h0)
      call sp%destroy()
      call ms%destroy()
   end subroutine test_idamp_zero_identity

   ! ------------------------------------------------------------------
   ! (4) analytic band builder reproduces the legacy cosine ramp
   ! ------------------------------------------------------------------
   subroutine test_band_builder(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(config_t) :: cfg
      type(ocean_state_t) :: state
      real(wp), parameter :: PI_T = acos(-1.0_wp)
      real(wp), parameter :: STRENGTH = 1.0e-3_wp
      integer, parameter :: BAND = 6
      integer, parameter :: NX_PHYS = 16, NY_PHYS = 8
      real(wp) :: expected_h, expected_u, expected_v
      integer :: d, j_probe

      checks: block
         call make_grid(grid, NX_PHYS, NY_PHYS, 1.0_wp, 1.0_wp)
         cfg%sim_type = "ocean"
         state%multilayer%nz_ml = NZ
         state%sponge%enable = .true.
         call state%init(grid)

         state%bc%west%bc_type = OBC_SPONGE
         state%bc%west%sponge_width = BAND
         state%bc%west%sponge_strength = STRENGTH

         cfg%ocean%sponge%enable = .true.
         cfg%ocean%sponge%damp_source = "band"
         cfg%ocean%sponge%target_source = "ic"

         call configure_ocean_sponge(cfg, state, grid, compute_rank=0)

         j_probe = grid%nghost + 2
         do d = 0, BAND - 1
            expected_h = STRENGTH*0.5_wp*(1.0_wp + cos(PI_T*real(d, wp)/real(BAND, wp)))
            if (abs(state%sponge%idamp_h(grid%nghost + 1 + d, j_probe) - expected_h) > 1.0e-14_wp) then
               call check(error, .false., "idamp_h band offset/rate mismatch (west edge)")
               exit checks
            end if
            expected_u = expected_h
            if (abs(state%sponge%idamp_u(grid%nghost + 2 + d, j_probe) - expected_u) > 1.0e-14_wp) then
               call check(error, .false., "idamp_u band offset/rate mismatch (west edge, +1 face offset)")
               exit checks
            end if
            expected_v = expected_h
            if (abs(state%sponge%idamp_v(grid%nghost + 1 + d, j_probe) - expected_v) > 1.0e-14_wp) then
               call check(error, .false., "idamp_v band offset/rate mismatch (west edge)")
               exit checks
            end if
         end do

         ! Zero outside the band (east interior) and in the ghosts.
         call check(error, state%sponge%idamp_h(grid%nghost + NX_PHYS, j_probe) == 0.0_wp, &
                    "idamp_h must be zero far from the sponge-tagged edge")
         if (allocated(error)) exit checks
         call check(error, state%sponge%idamp_h(1, 1) == 0.0_wp, &
                    "idamp_h must be zero in every ghost cell")
         if (allocated(error)) exit checks
         call check(error, .true.)

      end block checks
      call state%destroy()
   end subroutine test_band_builder

   ! ------------------------------------------------------------------
   ! (5) overlapping edge bands SUM their rates at corners
   ! ------------------------------------------------------------------
   subroutine test_band_corners_sum(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(config_t) :: cfg
      type(ocean_state_t) :: state
      real(wp), parameter :: PI_T = acos(-1.0_wp)
      real(wp), parameter :: STRENGTH = 2.0e-3_wp
      integer, parameter :: BAND = 4
      integer, parameter :: NX_PHYS = 16, NY_PHYS = 16
      real(wp) :: rate_w, rate_s, corner_val, single_edge_val
      integer :: i_corner, j_corner, i_single, j_single

      checks: block
         call make_grid(grid, NX_PHYS, NY_PHYS, 1.0_wp, 1.0_wp)
         cfg%sim_type = "ocean"
         state%multilayer%nz_ml = NZ
         state%sponge%enable = .true.
         call state%init(grid)

         state%bc%west%bc_type = OBC_SPONGE
         state%bc%west%sponge_width = BAND
         state%bc%west%sponge_strength = STRENGTH
         state%bc%south%bc_type = OBC_SPONGE
         state%bc%south%sponge_width = BAND
         state%bc%south%sponge_strength = STRENGTH

         cfg%ocean%sponge%enable = .true.
         cfg%ocean%sponge%damp_source = "band"

         call configure_ocean_sponge(cfg, state, grid, compute_rank=0)

         ! South-west corner cell: d_x = 1, d_y = 1 from each edge.
         i_corner = grid%nghost + 2
         j_corner = grid%nghost + 2
         rate_w = STRENGTH*0.5_wp*(1.0_wp + cos(PI_T*1.0_wp/real(BAND, wp)))
         rate_s = STRENGTH*0.5_wp*(1.0_wp + cos(PI_T*1.0_wp/real(BAND, wp)))
         corner_val = state%sponge%idamp_h(i_corner, j_corner)

         call check(error, abs(corner_val - (rate_w + rate_s)) < 1.0e-13_wp, &
                    "overlapping west+south bands must SUM their rates at a corner cell")
         if (allocated(error)) exit checks

         ! Single-edge cell (west band, far from the south band in y):
         ! equals the single west rate exactly.
         i_single = grid%nghost + 2
         j_single = grid%nghost + NY_PHYS - 1
         single_edge_val = state%sponge%idamp_h(i_single, j_single)
         call check(error, abs(single_edge_val - rate_w) < 1.0e-13_wp, &
                    "a single-edge cell must equal its single rate (not summed)")

      end block checks
      call state%destroy()
   end subroutine test_band_corners_sum

   ! ------------------------------------------------------------------
   ! (6) the tracer source closes the console heat budget
   ! ------------------------------------------------------------------
   subroutine test_budget_closes(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_sponge_t) :: sp
      real(wp), parameter :: LAMBDA = 0.01_wp
      real(wp), parameter :: DT = 10.0_wp
      real(wp), parameter :: T0 = 5.0_wp
      real(wp), parameter :: H0 = 5.0_wp
      real(wp), parameter :: T_REF = 15.0_wp
      real(wp), parameter :: TOL = 1.0e-12_wp
      integer, parameter :: NX_PHYS = 10, NY_PHYS = 6
      integer :: nx, ny
      real(wp) :: area_w, scale_v
      real(wp) :: ref_heat, total_heat, change_heat, src_heat, res_heat

      checks: block
         call make_grid(grid, NX_PHYS, NY_PHYS, 3000.0_wp, 3000.0_wp)
         nx = grid%nx_total
         ny = grid%ny_total
         area_w = grid%dx*grid%dy
         scale_v = area_w*RHO_WATER

         ms%nz_ml = NZ
         call ms%init(grid)
         ms%h_layer = H0
         ms%tracers(ms%idx_temperature)%hTr = T0*H0

         sp%enable = .true.
         call sp%init(grid, nz_ml=NZ, n_tracers=size(ms%tracers))
         sp%relax_tracers = .true.
         sp%relax_uv = .false.
         sp%idamp_h = LAMBDA
         sp%ref_tracer(:, :, :, ms%idx_temperature) = T_REF

         ref_heat = area_sum(ms%tracers(ms%idx_temperature)%hTr, nx, ny, area_w)*RHO_WATER

         ! Faithful two-stage RK2 (the sponge fires once per stage, like
         ! the surface flux): budget accumulates the raw sum over both
         ! stages via `+=`; RK2_STAGE_WEIGHT=0.5 (inside `ocean_budget_src`)
         ! is the correct weight (§3.5 of the plan).
         !$acc enter data copyin(ms, sp)
         call ms%enter_data()
         call sp%enter_data()
         call ocean_sponge_apply_maps(grid, sp, ms, DT)
         call ocean_sponge_apply_maps(grid, sp, ms, DT)
         !$acc update self(ms%tracers(ms%idx_temperature)%hTr, ms%heat_budget_sponge)
         call sp%exit_data()
         call ms%exit_data()
         !$acc exit data delete(ms, sp)

         ! Manual RK2 average (physics 0.5) — matches the production
         ! rk2_average step (§3.5: the accumulator itself is NOT averaged).
         ms%tracers(ms%idx_temperature)%hTr = &
            0.5_wp*(T0*H0 + ms%tracers(ms%idx_temperature)%hTr)

         total_heat = area_sum(ms%tracers(ms%idx_temperature)%hTr, nx, ny, area_w)*RHO_WATER
         change_heat = total_heat - ref_heat

         src_heat = ocean_budget_src(ocean_heat_src_sum(0.0_wp, 0.0_wp, &
                                                        area_sum(ms%heat_budget_sponge, nx, ny, area_w)))
         res_heat = change_heat - src_heat

         call check(error, abs(src_heat) > 0.0_wp, &
                    "sponge heat source must be strictly nonzero (T_ref != T)")
         if (allocated(error)) exit checks
         call check(error, abs(res_heat) <= TOL*abs(ref_heat), &
                    "change_heat - src_heat must close to round-off with the sponge on")
         if (allocated(error)) exit checks
         call check(error, abs(src_heat) > 1.0e3_wp*max(abs(res_heat), tiny(1.0_wp)), &
                    "the sponge source must dominate the residual (load-bearing, not noise)")

      end block checks
      call sp%destroy()
      call ms%destroy()
   end subroutine test_budget_closes

   ! ------------------------------------------------------------------
   ! (7) the reference snapshot is the IC, not whatever `ms` becomes later
   ! ------------------------------------------------------------------
   subroutine test_reference_is_ic_not_restart(error)
      !! Directly tests the property the driver's call-site ordering relies
      !! on (docs/plans/PLAN_PR23_real_sponge.md §13.1 item 4): call
      !! `ocean_sponge_snapshot_reference` once, on an `ms` holding the
      !! seeded IC (`T_IC`); THEN mutate `ms` to a different value (`T_
      !! RESTART`, standing in for what a warm-restart read would overwrite
      !! it with AFTER the snapshot ran, exactly as the driver orders it —
      !! see `rdb_driver.F90` around `ocean_state_seed_from_cfg` /
      !! `ocean_sponge_snapshot_reference` / `ocean_state_restart_read`).
      !! The snapshot must be inert to the later mutation: `sp%ref_tracer`
      !! must equal `T_IC`, bit-for-bit, never `T_RESTART`.  Under the
      !! ORIGINAL (rejected) design — snapshotting at `configure_ocean_sponge`
      !! time, which the driver runs AFTER the restart read — this assertion
      !! would fail against `T_RESTART`.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_sponge_t) :: sp
      real(wp), parameter :: T_IC = 8.0_wp
      real(wp), parameter :: T_RESTART = 25.0_wp
      real(wp), parameter :: H0 = 4.0_wp
      integer, parameter :: NX_PHYS = 4, NY_PHYS = 4
      integer :: idx_T

      checks: block
         call make_grid(grid, NX_PHYS, NY_PHYS, 1.0_wp, 1.0_wp)
         ms%nz_ml = NZ
         call ms%init(grid)
         ms%h_layer = H0
         idx_T = ms%idx_temperature
         ms%tracers(idx_T)%hTr = T_IC*H0

         sp%enable = .true.
         call sp%init(grid, nz_ml=NZ, n_tracers=size(ms%tracers))

         ! Snapshot while ms holds the IC — mirrors the driver calling this
         ! right after ocean_state_seed_from_cfg, BEFORE any restart read.
         call ocean_sponge_snapshot_reference(sp, grid, ms)

         ! Mutate ms AFTER the snapshot — stands in for a warm-restart read
         ! overwriting the seeded IC with a different mid-run state.
         ms%tracers(idx_T)%hTr = T_RESTART*H0

         call check(error, all(sp%ref_tracer(:, :, :, idx_T) == T_IC), &
                    "sponge reference must equal the IC (T_IC), bit-for-bit, "// &
                    "and must be UNAFFECTED by a later mutation of ms "// &
                    "(a restart read arriving after the snapshot)")
         if (allocated(error)) exit checks
         call check(error,.not. any(sp%ref_tracer(:, :, :, idx_T) == T_RESTART), &
                    "sponge reference must NEVER pick up the post-snapshot "// &
                    "(restart-like) mutation")

      end block checks
      call sp%destroy()
      call ms%destroy()
   end subroutine test_reference_is_ic_not_restart

end module test_ocean_sponge
