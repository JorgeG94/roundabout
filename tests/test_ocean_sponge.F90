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
                               ocean_sponge_snapshot_reference, &
                               ocean_sponge_refresh_target, sponge_band_alpha, &
                               SPONGE_RAMP_COSINE, SPONGE_RAMP_LINEAR
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
                               test_reference_is_ic_not_restart), &
                  new_unittest("linear_z_target_is_the_analytic_geopotential_profile", &
                               test_linear_z_profile), &
                  new_unittest("linear_z_target_tracks_the_live_layer_geometry", &
                               test_linear_z_tracks_layers), &
                  new_unittest("linear_z_target_only_touches_sponge_cells", &
                               test_linear_z_masked), &
                  new_unittest("linear_z_relaxes_at_the_analytic_rate", &
                               test_linear_z_relax_rate), &
                  new_unittest("linear_z_target_may_differ_from_the_ic", &
                               test_linear_z_target_not_ic), &
                  new_unittest("target_source_ic_is_bit_identical_under_refresh", &
                               test_refresh_noop_for_ic), &
                  new_unittest("band_ramp_linear_is_isomip_eq20_at_cell_centres", &
                               test_band_ramp_linear) &
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

   ! ==================================================================
   ! PR-23b: the ANALYTIC `target_source = "linear_z"` reference
   ! ==================================================================

   subroutine seed_linear_z_case(grid, ms, sp, nx_phys, ny_phys, lambda, &
                                 t_ref, dt_dz, s_ref, ds_dz, draft_slope)
      !! Common fixture: a `nx_phys x ny_phys x NZ` column stack with a
      !! non-uniform layer thickness (so a target built on layer INDEX
      !! could never pass) and a draft that slopes in x (so a target that
      !! forgets `z_top` could never pass either).  The sponge damps the
      !! whole interior at `lambda` and leaves the ghosts alone.
      type(hgrid_t), intent(inout) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_sponge_t), intent(inout) :: sp
      integer, intent(in) :: nx_phys, ny_phys
      real(wp), intent(in) :: lambda, t_ref, dt_dz, s_ref, ds_dz, draft_slope
      integer :: i, j, k

      call make_grid(grid, nx_phys, ny_phys, 1000.0_wp, 1000.0_wp)
      ms%nz_ml = NZ
      call ms%init(grid)
      ! Deliberately uneven layers: 10, 20, 30, 40 m bottom-up.
      do k = 1, NZ
         ms%h_layer(:, :, k) = 10.0_wp*real(k, wp)
      end do

      sp%enable = .true.
      call sp%init(grid, nz_ml=NZ, n_tracers=size(ms%tracers))
      sp%relax_tracers = .true.
      sp%relax_uv = .false.
      sp%target_source = "linear_z"
      sp%lin_t_ref = t_ref
      sp%lin_dt_dz = dt_dz
      sp%lin_s_ref = s_ref
      sp%lin_ds_dz = ds_dz
      sp%idx_t = ms%idx_temperature
      sp%idx_s = ms%idx_salinity
      sp%idamp_h = 0.0_wp
      do j = grid%nghost + 1, grid%nghost + ny_phys
         do i = grid%nghost + 1, grid%nghost + nx_phys
            sp%idamp_h(i, j) = lambda
         end do
      end do
      ! Ice base deepening toward +x: z_top = slope * x_centre.
      do j = 1, grid%ny_total
         do i = 1, grid%nx_total
            sp%z_top(i, j) = draft_slope*(real(i - grid%nghost, wp) - 0.5_wp)*grid%dx
         end do
      end do
   end subroutine seed_linear_z_case

   pure function expected_z_ctr(h_col, z_top, k, nz) result(z)
      !! Layer-centre geopotential DEPTH, rebuilt independently of the
      !! production recurrence: `z_top` plus every layer ABOVE k plus half
      !! of k.  Bottom-up storage, so "above k" is k+1 .. nz.
      integer, intent(in) :: k, nz
      real(wp), intent(in) :: h_col(nz)
      real(wp), intent(in) :: z_top
      real(wp) :: z
      integer :: m
      z = z_top + 0.5_wp*h_col(k)
      do m = k + 1, nz
         z = z + h_col(m)
      end do
   end function expected_z_ctr

   subroutine test_linear_z_profile(error)
      !! The refreshed target must equal `v_ref - dv_dz*z_ctr` at every
      !! layer centre, with `z_ctr` measured from the `z = 0` datum
      !! THROUGH the ice draft.  Uneven layers + a sloping lid mean a
      !! layer-index profile or a draft-blind profile both fail.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_sponge_t) :: sp
      ! ISOMIP+ WARM (Asay-Davis et al. 2016 Tables 4 + 6), in the
      ! z-positive-UP convention: T0 = -1.9 degC at z = 0, Tbot = 1.0 degC
      ! at z = z_b,deep = -720 m, so dT/dz = -(1.0 - (-1.9))/720 < 0 --
      ! WARM is thermally UNSTABLE (CDW at depth), stabilised by the
      ! salinity gradient.  Sbot = 34.7 gives dS/dz = -0.9/720.
      real(wp), parameter :: T_REF = -1.9_wp, DT_DZ = -4.027777778e-3_wp
      real(wp), parameter :: S_REF = 33.8_wp, DS_DZ = -1.25e-3_wp
      integer, parameter :: NX_PHYS = 5, NY_PHYS = 3
      real(wp) :: z, want_t, want_s, worst
      integer :: i, j, k

      checks: block
         call seed_linear_z_case(grid, ms, sp, NX_PHYS, NY_PHYS, 1.0e-4_wp, &
                                 T_REF, DT_DZ, S_REF, DS_DZ, 0.05_wp)

         !$acc enter data copyin(ms, sp)
         call ms%enter_data()
         call sp%enter_data()
         call ocean_sponge_refresh_target(grid, sp, ms)
         !$acc update self(sp%ref_tracer)
         call sp%exit_data()
         call ms%exit_data()
         !$acc exit data delete(ms, sp)

         worst = 0.0_wp
         do k = 1, NZ
            do j = grid%nghost + 1, grid%nghost + NY_PHYS
               do i = grid%nghost + 1, grid%nghost + NX_PHYS
                  z = expected_z_ctr(ms%h_layer(i, j, :), sp%z_top(i, j), k, NZ)
                  want_t = T_REF - DT_DZ*z
                  want_s = S_REF - DS_DZ*z
                  worst = max(worst, abs(sp%ref_tracer(i, j, k, ms%idx_temperature) - want_t))
                  worst = max(worst, abs(sp%ref_tracer(i, j, k, ms%idx_salinity) - want_s))
               end do
            end do
         end do
         call check(error, worst < 1.0e-12_wp, &
                    "linear_z target must be v_ref - dv_dz*z_ctr at every layer centre, "// &
                    "with z_ctr measured through the ice draft")
         if (allocated(error)) exit checks

         ! A draft-blind target would agree at i = 1 and disagree at
         ! i = NX_PHYS; assert the x-variation is real so the test cannot
         ! pass with z_top dropped.
         call check(error, abs(sp%ref_tracer(grid%nghost + 1, grid%nghost + 1, NZ, &
                                             ms%idx_temperature) &
                               - sp%ref_tracer(grid%nghost + NX_PHYS, grid%nghost + 1, NZ, &
                                               ms%idx_temperature)) > 1.0e-6_wp, &
                    "the sloping draft must make the target vary in x")
      end block checks
      call sp%destroy()
      call ms%destroy()
   end subroutine test_linear_z_profile

   subroutine test_linear_z_tracks_layers(error)
      !! The design decision, gated: the target is RE-EVALUATED on the
      !! live layer geometry, not frozen at t = 0.  Thicken every layer
      !! by 50 % (what an ALE regrid or a rising free surface does) and
      !! refresh again; the target must move with the new layer centres.
      !! A frozen target would be bit-identical and fail this.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_sponge_t) :: sp
      real(wp), parameter :: T_REF = 0.0_wp, DT_DZ = 1.0e-2_wp
      integer, parameter :: NX_PHYS = 3, NY_PHYS = 3
      real(wp) :: before, z, want
      integer :: ip, jp, k

      checks: block
         call seed_linear_z_case(grid, ms, sp, NX_PHYS, NY_PHYS, 1.0e-4_wp, &
                                 T_REF, DT_DZ, 35.0_wp, 0.0_wp, 0.0_wp)
         ip = grid%nghost + 2
         jp = grid%nghost + 2

         !$acc enter data copyin(ms, sp)
         call ms%enter_data()
         call sp%enter_data()
         call ocean_sponge_refresh_target(grid, sp, ms)
         !$acc update self(sp%ref_tracer)
         before = sp%ref_tracer(ip, jp, 1, ms%idx_temperature)

         ms%h_layer = 1.5_wp*ms%h_layer
         !$acc update device(ms%h_layer)
         call ocean_sponge_refresh_target(grid, sp, ms)
         !$acc update self(sp%ref_tracer)
         call sp%exit_data()
         call ms%exit_data()
         !$acc exit data delete(ms, sp)

         k = 1
         z = expected_z_ctr(ms%h_layer(ip, jp, :), sp%z_top(ip, jp), k, NZ)
         want = T_REF - DT_DZ*z
         call check(error, abs(sp%ref_tracer(ip, jp, k, ms%idx_temperature) - want) < 1.0e-12_wp, &
                    "the refreshed target must follow the LIVE layer centres")
         if (allocated(error)) exit checks
         call check(error, abs(sp%ref_tracer(ip, jp, k, ms%idx_temperature) - before) > 1.0e-6_wp, &
                    "a target frozen at t = 0 would not have moved — the whole point")
      end block checks
      call sp%destroy()
      call ms%destroy()
   end subroutine test_linear_z_tracks_layers

   subroutine test_linear_z_masked(error)
      !! `Idamp = 0` IS the sponge mask for the refresh too: a cell
      !! outside the band keeps whatever `ref_tracer` it already had
      !! (in production, the IC snapshot), bit-for-bit.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_sponge_t) :: sp
      real(wp), parameter :: SENTINEL = -12345.0_wp
      integer, parameter :: NX_PHYS = 4, NY_PHYS = 4

      checks: block
         call seed_linear_z_case(grid, ms, sp, NX_PHYS, NY_PHYS, 1.0e-4_wp, &
                                 0.0_wp, 1.0e-2_wp, 35.0_wp, 0.0_wp, 0.0_wp)
         ! Carve one interior cell OUT of the band and mark its target.
         sp%idamp_h(grid%nghost + 2, grid%nghost + 2) = 0.0_wp
         sp%ref_tracer(grid%nghost + 2, grid%nghost + 2, :, ms%idx_temperature) = SENTINEL
         sp%ref_tracer(1, 1, :, ms%idx_temperature) = SENTINEL

         !$acc enter data copyin(ms, sp)
         call ms%enter_data()
         call sp%enter_data()
         call ocean_sponge_refresh_target(grid, sp, ms)
         !$acc update self(sp%ref_tracer)
         call sp%exit_data()
         call ms%exit_data()
         !$acc exit data delete(ms, sp)

         call check(error, all(sp%ref_tracer(grid%nghost + 2, grid%nghost + 2, :, &
                                             ms%idx_temperature) == SENTINEL), &
                    "a cell with Idamp = 0 must be left bit-for-bit alone")
         if (allocated(error)) exit checks
         call check(error, all(sp%ref_tracer(1, 1, :, ms%idx_temperature) == SENTINEL), &
                    "ghost cells (Idamp = 0) must be left alone")
         if (allocated(error)) exit checks
         call check(error, sp%ref_tracer(grid%nghost + 1, grid%nghost + 1, NZ, &
                                         ms%idx_temperature) /= SENTINEL, &
                    "...while a banded cell IS refreshed")
      end block checks
      call sp%destroy()
      call ms%destroy()
   end subroutine test_linear_z_masked

   subroutine test_linear_z_relax_rate(error)
      !! A column displaced from the analytic target must relax toward it
      !! at exactly the discrete rate the kernel advertises.  With
      !! `phi <- phi*decay + tgt*(1-decay)` and `decay = exp(-Idamp*dt)`,
      !! the DEVIATION obeys `d_n = d_0 * decay**n = d_0 * exp(-n*dt/tau)`
      !! EXACTLY (not to O(dt)): the per-step map is affine with the same
      !! fixed point every step, so the errors compose as a pure power.
      !! Here `tau = 1/Idamp = 8640 s` (0.1 days, the ISOMIP+ value) and
      !! 20 steps of dt = 864 s is exactly 2 e-foldings.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_sponge_t) :: sp
      real(wp), parameter :: TAU = 8640.0_wp          ! 0.1 days
      real(wp), parameter :: LAMBDA = 1.0_wp/TAU
      real(wp), parameter :: DT = 864.0_wp
      integer, parameter :: NSTEPS = 20               ! n*dt/tau = 2
      real(wp), parameter :: T_OFFSET = 3.0_wp        ! IC minus target
      integer, parameter :: NX_PHYS = 3, NY_PHYS = 3
      real(wp) :: tgt, got, dev, want_dev, worst
      integer :: i, j, k, step

      checks: block
         call seed_linear_z_case(grid, ms, sp, NX_PHYS, NY_PHYS, LAMBDA, &
                                 -1.9_wp, -4.0e-3_wp, 34.2_wp, 0.0_wp, 0.02_wp)

         !$acc enter data copyin(ms, sp)
         call ms%enter_data()
         call sp%enter_data()
         call ocean_sponge_refresh_target(grid, sp, ms)
         !$acc update self(sp%ref_tracer)
         ! Seed the state exactly T_OFFSET above the analytic target.
         do k = 1, NZ
            ms%tracers(ms%idx_temperature)%hTr(:, :, k) = &
               (sp%ref_tracer(:, :, k, ms%idx_temperature) + T_OFFSET)*ms%h_layer(:, :, k)
         end do
         !$acc update device(ms%tracers(ms%idx_temperature)%hTr)

         do step = 1, NSTEPS
            call ocean_sponge_apply_maps(grid, sp, ms, DT)
         end do
         !$acc update self(ms%tracers(ms%idx_temperature)%hTr)
         call sp%exit_data()
         call ms%exit_data()
         !$acc exit data delete(ms, sp)

         want_dev = T_OFFSET*exp(-real(NSTEPS, wp)*DT/TAU)
         worst = 0.0_wp
         do k = 1, NZ
            do j = grid%nghost + 1, grid%nghost + NY_PHYS
               do i = grid%nghost + 1, grid%nghost + NX_PHYS
                  tgt = sp%ref_tracer(i, j, k, ms%idx_temperature)
                  got = ms%tracers(ms%idx_temperature)%hTr(i, j, k)/ms%h_layer(i, j, k)
                  dev = got - tgt
                  worst = max(worst, abs(dev - want_dev))
               end do
            end do
         end do
         call check(error, worst < 1.0e-12_wp*T_OFFSET, &
                    "deviation from the analytic target must decay as exp(-t/tau) exactly")
         if (allocated(error)) exit checks
         ! Sanity: two e-foldings really did most of the work.
         call check(error, want_dev < 0.14_wp*T_OFFSET .and. want_dev > 0.13_wp*T_OFFSET, &
                    "20 steps of dt = tau/10 is exactly 2 e-foldings (exp(-2) = 0.1353)")
      end block checks
      call sp%destroy()
      call ms%destroy()
   end subroutine test_linear_z_relax_rate

   subroutine test_linear_z_target_not_ic(error)
      !! The Ocean1 / Ocean2 property: the restoring profile is NOT the
      !! initial condition.  Seed the column uniformly COLD, target a
      !! WARM linear profile, and check the column moves toward the WARM
      !! target and AWAY from the IC — which `target_source="ic"` could
      !! never do, because there the target IS the IC and the state would
      !! be a fixed point.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_sponge_t) :: sp
      real(wp), parameter :: T_COLD = -1.9_wp        ! ISOMIP+ COLD: uniform
      real(wp), parameter :: T_WARM_SFC = -1.9_wp    ! ISOMIP+ WARM: T0
      real(wp), parameter :: T_WARM_DTDZ = -4.027777778e-3_wp  ! -(1.0-(-1.9))/720
      real(wp), parameter :: LAMBDA = 1.0e-3_wp, DT = 100.0_wp
      integer, parameter :: NSTEPS = 5
      integer, parameter :: NX_PHYS = 3, NY_PHYS = 3
      real(wp) :: t_bed_before, t_bed_after, tgt_bed
      integer :: ip, jp, step

      checks: block
         call seed_linear_z_case(grid, ms, sp, NX_PHYS, NY_PHYS, LAMBDA, &
                                 T_WARM_SFC, T_WARM_DTDZ, 34.7_wp, -1.25e-3_wp, 0.0_wp)
         ip = grid%nghost + 2
         jp = grid%nghost + 2
         ! Uniform COLD initial condition.
         ms%tracers(ms%idx_temperature)%hTr = T_COLD*ms%h_layer
         t_bed_before = T_COLD

         !$acc enter data copyin(ms, sp)
         call ms%enter_data()
         call sp%enter_data()
         call ocean_sponge_refresh_target(grid, sp, ms)
         do step = 1, NSTEPS
            call ocean_sponge_apply_maps(grid, sp, ms, DT)
         end do
         !$acc update self(sp%ref_tracer)
         !$acc update self(ms%tracers(ms%idx_temperature)%hTr)
         call sp%exit_data()
         call ms%exit_data()
         !$acc exit data delete(ms, sp)

         tgt_bed = sp%ref_tracer(ip, jp, 1, ms%idx_temperature)
         t_bed_after = ms%tracers(ms%idx_temperature)%hTr(ip, jp, 1)/ms%h_layer(ip, jp, 1)

         call check(error, tgt_bed > T_COLD + 0.1_wp, &
                    "the WARM linear target must be warmer than the COLD IC at the bed")
         if (allocated(error)) exit checks
         call check(error, t_bed_after > t_bed_before, &
                    "the bed layer must warm toward the target, away from the IC")
         if (allocated(error)) exit checks
         call check(error, t_bed_after < tgt_bed, &
                    "...without overshooting the target")
      end block checks
      call sp%destroy()
      call ms%destroy()
   end subroutine test_linear_z_target_not_ic

   subroutine test_refresh_noop_for_ic(error)
      !! Default-off contract: with `target_source = "ic"` (the default)
      !! the refresh is a bit-for-bit no-op, so every shipped namelist
      !! and every existing test is unaffected by this feature existing.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_sponge_t) :: sp
      real(wp), allocatable :: before(:, :, :, :)
      integer, parameter :: NX_PHYS = 4, NY_PHYS = 4

      checks: block
         call seed_linear_z_case(grid, ms, sp, NX_PHYS, NY_PHYS, 1.0e-3_wp, &
                                 0.0_wp, 1.0e-2_wp, 35.0_wp, -1.0e-3_wp, 0.01_wp)
         sp%target_source = "ic"
         sp%ref_tracer = 7.25_wp
         before = sp%ref_tracer

         !$acc enter data copyin(ms, sp)
         call ms%enter_data()
         call sp%enter_data()
         call ocean_sponge_refresh_target(grid, sp, ms)
         !$acc update self(sp%ref_tracer)
         call sp%exit_data()
         call ms%exit_data()
         !$acc exit data delete(ms, sp)

         call check(error, all(sp%ref_tracer == before), &
                    "target_source='ic' must make the refresh a bit-for-bit no-op")
      end block checks
      call sp%destroy()
      call ms%destroy()
   end subroutine test_refresh_noop_for_ic

   subroutine test_band_ramp_linear(error)
      !! `ramp = "linear"` is the CELL-CENTRE evaluation of ISOMIP+
      !! Eq. (20), `gamma(x) = gamma0*(x - x_r0)/(x_r1 - x_r0)`
      !! (Asay-Davis et al. 2016, their gamma0 = 10/day over
      !! 790 <= x <= 800 km).  On a 2 km grid that band is 5 cells, and
      !! the cell centres sit at x = 799, 797, 795, 793, 791 km, i.e.
      !! gamma/gamma0 = 0.9, 0.7, 0.5, 0.3, 0.1 running inward from the
      !! wall (d = 0 .. 4) — exactly (band - d - 0.5)/band with band = 5.
      !! The cosine branch must be untouched.
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: PI_T = acos(-1.0_wp)
      real(wp), parameter :: WANT(5) = [0.9_wp, 0.7_wp, 0.5_wp, 0.3_wp, 0.1_wp]
      integer, parameter :: BAND = 5
      real(wp) :: band_mean
      integer :: d

      checks: block
         do d = 0, BAND - 1
            if (abs(sponge_band_alpha(d, BAND, SPONGE_RAMP_LINEAR) - WANT(d + 1)) &
                > 1.0e-14_wp) then
               call check(error, .false., "linear ramp must be ISOMIP+ Eq. (20) at cell centres")
               exit checks
            end if
            if (abs(sponge_band_alpha(d, BAND, SPONGE_RAMP_COSINE) &
                    - 0.5_wp*(1.0_wp + cos(PI_T*real(d, wp)/real(BAND, wp)))) > 1.0e-14_wp) then
               call check(error, .false., "cosine ramp must stay the legacy shape bit-for-bit")
               exit checks
            end if
         end do
         ! The linear ramp integrates to exactly half the wall strength,
         ! which is what "ramps linearly to zero" means for a band mean.
         band_mean = 0.0_wp
         do d = 0, BAND - 1
            band_mean = band_mean + sponge_band_alpha(d, BAND, SPONGE_RAMP_LINEAR)
         end do
         band_mean = band_mean/real(BAND, wp)
         call check(error, abs(band_mean - 0.5_wp) < 1.0e-14_wp, &
                    "the linear band mean must be 0.5")
      end block checks
   end subroutine test_band_ramp_linear

end module test_ocean_sponge
