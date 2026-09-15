!! Sea-ice diagnostics tests (ice_conc / ice_thick / ice_speed / ice_u /
!! ice_v) — the diag-manager wiring for the ice concentration/thickness/
!! velocity read-outs (SPEC_ice-diags.md).  Pure read-out: no ice physics
!! is exercised beyond stamping the state's ice arrays directly, so these
!! tests pin the FILL math (the two-mode gather inlined into
!! `fill_ice_conc_thick_impl`) against the canonical
!! `ice_cell_concentration_impl` (`rdb_ice_state`), plus the registry-level
!! bit-identity guarantee (ice off ⇒ no ice vars registered).
module test_ocean_ice_diags
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_state, only: ocean_state_t, ocean_state_enter_data, ocean_state_exit_data
   use rdb_ocean_diag, only: DIAG_OP_INSTANT
   use rdb_ocean_diag_fills, only: register_default_diags, is_canonical_diag_name
   use rdb_ocean_diag_derived, only: register_derived, derived_catalog_size, &
                                     derived_catalog_name
   use rdb_ocean_metrics, only: metrics_fill_cartesian, metrics_finalize
   use rdb_ocean_console_stats, only: compute_ice_totals
   use rdb_ice_state, only: ice_cell_concentration_impl
   use rdb_ice_column, only: ICE_RHO_ICE
   use testdrive, only: error_type, check, new_unittest, unittest_type
   implicit none
   private

   public :: collect_ocean_ice_diags_tests

   integer, parameter :: NX = 6, NY = 4, NZ = 3
   real(wp), parameter :: DX = 1.0_wp
   real(wp), parameter :: TOL = 1.0e-12_wp

contains

   subroutine collect_ocean_ice_diags_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("ice_conc_sum_part_size", test_ice_conc_sum_part_size), &
                  new_unittest("ice_conc_lumped_ncat1", test_ice_conc_lumped_ncat1), &
                  new_unittest("ice_conc_supersaturated_clamp", test_ice_conc_supersaturated_clamp), &
                  new_unittest("ice_thick_mass_over_rho", test_ice_thick_mass_over_rho), &
                  new_unittest("ice_speed_face_average", test_ice_speed_face_average), &
                  new_unittest("ice_diags_land_mask", test_ice_diags_land_mask), &
                  new_unittest("ice_fill_matches_concentration_impl", &
                               test_ice_fill_matches_concentration_impl), &
                  new_unittest("console_ice_totals_reference", test_console_ice_totals_reference), &
                  new_unittest("ice_off_no_ice_diags", test_ice_off_no_ice_diags), &
                  new_unittest("ice_canonical_names", test_ice_canonical_names), &
                  new_unittest("ice_derived_catalog_has_velocity", &
                               test_ice_derived_catalog_has_velocity), &
                  new_unittest("ice_speed_zero_when_ice_off", test_ice_speed_zero_when_ice_off) &
                  ]
   end subroutine collect_ocean_ice_diags_tests

   subroutine setup_state(grid, state)
      !! Plain ice-OFF setup (mirrors test_ocean_diag's setup_state) — used
      !! by the registry-level bit-identity + ice-off derived-fill tests.
      type(hgrid_t), intent(inout) :: grid
      type(ocean_state_t), intent(inout) :: state
      call grid%init(NX, NY, 1, DX, DX)
      state%multilayer%nz_ml = NZ
      call state%init(grid)
      call metrics_fill_cartesian(state%metrics, grid, grid%dx, grid%dy)
      call metrics_finalize(state%metrics)
   end subroutine setup_state

   subroutine setup_state_ice(grid, state, ncat)
      !! setup_state + the sea-ice slot enabled (ncat categories), so
      !! `state%init` allocates + zero-inits the ice arrays and
      !! `register_default_diags` exposes ice_conc/ice_thick.
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

   function find_var(state, name) result(idx)
      !! Scan `state%diag%vars(1:nvars)%name` for `name`; 0 if absent.
      type(ocean_state_t), intent(in) :: state
      character(len=*), intent(in) :: name
      integer :: idx, i
      idx = 0
      do i = 1, state%diag%nvars
         if (trim(state%diag%vars(i)%name) == trim(name)) then
            idx = i
            return
         end if
      end do
   end function find_var

   subroutine test_ice_conc_sum_part_size(error)
      !! ncat=3 ITD mode: part_size(1)=0.3, part_size(2)=0.4 (rest 0) on an
      !! all-wet grid ⇒ ice_conc == 0.7 everywhere (interior + ghosts).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      integer :: iconc
      checks: block
         call setup_state_ice(grid, state, ncat=3)
         state%ice%part_size(:, :, 1) = 0.3_wp
         state%ice%part_size(:, :, 2) = 0.4_wp
         state%ice%m_ice(:, :, 1) = 500.0_wp
         state%ice%m_ice(:, :, 2) = 800.0_wp

         call register_default_diags(state, dt_out=1.0_wp)
         iconc = find_var(state, "ice_conc")
         call check(error, iconc > 0, "ice_conc must be registered when ice is on")
         if (allocated(error)) exit checks

         call ocean_state_enter_data(state)
         call state%diag%step(state, dt=1.0_wp, t=1.0_wp)
         call ocean_state_exit_data(state)

         call check(error, maxval(abs(state%diag%vars(iconc)%output_buffer - 0.7_wp)) < TOL, &
                    "ice_conc should equal sum(part_size(1:2)) = 0.7 everywhere")
      end block checks
      call state%destroy()
   end subroutine test_ice_conc_sum_part_size

   subroutine test_ice_conc_lumped_ncat1(error)
      !! ncat=1 legacy lumped mode: m_ice > 0 on half the columns (i<=4)
      !! ⇒ ice_conc == 1 there, 0 elsewhere.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      integer :: iconc, i, j
      real(wp) :: expected
      checks: block
         call setup_state_ice(grid, state, ncat=1)
         do j = 1, NY
            do i = 1, NX
               if (i <= 4) state%ice%m_ice(i, j, 1) = 100.0_wp
            end do
         end do

         call register_default_diags(state, dt_out=1.0_wp)
         iconc = find_var(state, "ice_conc")
         call check(error, iconc > 0, "ice_conc must be registered")
         if (allocated(error)) exit checks

         call ocean_state_enter_data(state)
         call state%diag%step(state, dt=1.0_wp, t=1.0_wp)
         call ocean_state_exit_data(state)

         do j = 1, NY
            do i = 1, NX
               expected = merge(1.0_wp, 0.0_wp, i <= 4)
               call check(error, &
                          abs(state%diag%vars(iconc)%output_buffer(i, j, 1) - expected) < TOL, &
                          "ice_conc lumped 0/1 pattern mismatch")
               if (allocated(error)) exit checks
            end do
            if (allocated(error)) exit checks
         end do
      end block checks
      call state%destroy()
   end subroutine test_ice_conc_lumped_ncat1

   subroutine test_ice_conc_supersaturated_clamp(error)
      !! The `min(1, Σ part_size)` clamp: ncat=3 with part_size(1)=0.7,
      !! part_size(2)=0.6 (Σ=1.3 > 1) ⇒ ice_conc == 1.0 everywhere.  Guards
      !! the min() — every other ITD test sums ≤ 1, so deleting the clamp
      !! would pass the rest of the suite.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      integer :: iconc
      checks: block
         call setup_state_ice(grid, state, ncat=3)
         state%ice%part_size(:, :, 1) = 0.7_wp
         state%ice%part_size(:, :, 2) = 0.6_wp

         call register_default_diags(state, dt_out=1.0_wp)
         iconc = find_var(state, "ice_conc")
         call check(error, iconc > 0, "ice_conc must be registered")
         if (allocated(error)) exit checks

         call ocean_state_enter_data(state)
         call state%diag%step(state, dt=1.0_wp, t=1.0_wp)
         call ocean_state_exit_data(state)

         call check(error, maxval(abs(state%diag%vars(iconc)%output_buffer - 1.0_wp)) < TOL, &
                    "ice_conc must clamp to 1.0 when Σ part_size = 1.3 > 1")
      end block checks
      call state%destroy()
   end subroutine test_ice_conc_supersaturated_clamp

   subroutine test_ice_thick_mass_over_rho(error)
      !! ice_thick = mice / ICE_RHO_ICE for both the ITD (ncat=3) and the
      !! legacy lumped (ncat=1) mode.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      integer :: ithick
      checks: block
         ! ncat=3: part_size(1)=0.5, m_ice(1)=905 => mice=452.5 => 0.5 m.
         call setup_state_ice(grid, state, ncat=3)
         state%ice%part_size(:, :, 1) = 0.5_wp
         state%ice%m_ice(:, :, 1) = 905.0_wp

         call register_default_diags(state, dt_out=1.0_wp)
         ithick = find_var(state, "ice_thick")
         call check(error, ithick > 0, "ice_thick must be registered")
         if (allocated(error)) exit checks

         call ocean_state_enter_data(state)
         call state%diag%step(state, dt=1.0_wp, t=1.0_wp)
         call ocean_state_exit_data(state)

         call check(error, maxval(abs(state%diag%vars(ithick)%output_buffer - 0.5_wp)) < TOL, &
                    "ice_thick (ITD) should equal 452.5/905 = 0.5 m")
      end block checks
      call state%destroy()

      checks2: block
         ! ncat=1: m_ice = 452.5 => 0.5 m directly (per-cell legacy mode).
         call setup_state_ice(grid, state, ncat=1)
         state%ice%m_ice(:, :, 1) = 452.5_wp

         call register_default_diags(state, dt_out=1.0_wp)
         ithick = find_var(state, "ice_thick")
         call check(error, ithick > 0, "ice_thick must be registered (ncat=1)")
         if (allocated(error)) exit checks2

         call ocean_state_enter_data(state)
         call state%diag%step(state, dt=1.0_wp, t=1.0_wp)
         call ocean_state_exit_data(state)

         call check(error, maxval(abs(state%diag%vars(ithick)%output_buffer - 0.5_wp)) < TOL, &
                    "ice_thick (lumped) should equal 452.5/905 = 0.5 m")
      end block checks2
      call state%destroy()
   end subroutine test_ice_thick_mass_over_rho

   subroutine test_ice_speed_face_average(error)
      !! Derived ice_speed/ice_u/ice_v.  u_ice is SPATIALLY VARYING
      !! (`u_ice(i,j)=0.05*i`, so u(i)/=u(i+1)) to pin the two-face
      !! average — a no-average `buf=u_ice(i,j)` would give 0.15 at cell
      !! (3,2), not the correct 0.175.  v_ice is uniform-NEGATIVE (-0.2)
      !! to pin the sign (a stray abs would flip it).  Asserted at the
      !! interior T-cell (3,2): with `grid%init(NX,NY,1)` (nghost=1) the
      !! physical interior is i in [2,7], j in [2,5], and both u_ice(3,2)
      !! and u_ice(4,2) are in bounds (u_ice is (nx_total+1, ny_total)).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      integer :: ispd, iu, iv, i, j
      real(wp), parameter :: EXP_U = 0.5_wp*(0.05_wp*3.0_wp + 0.05_wp*4.0_wp)  ! 0.175
      real(wp), parameter :: EXP_V = -0.2_wp
      real(wp), parameter :: EXP_SPD = sqrt(EXP_U*EXP_U + EXP_V*EXP_V)
      checks: block
         call setup_state_ice(grid, state, ncat=3)
         do j = 1, size(state%ice%u_ice, 2)
            do i = 1, size(state%ice%u_ice, 1)
               state%ice%u_ice(i, j) = 0.05_wp*real(i, wp)
            end do
         end do
         state%ice%v_ice = -0.2_wp

         call register_derived(state, "ice_speed", time_op=DIAG_OP_INSTANT, dt_out=1.0_wp)
         call register_derived(state, "ice_u", time_op=DIAG_OP_INSTANT, dt_out=1.0_wp)
         call register_derived(state, "ice_v", time_op=DIAG_OP_INSTANT, dt_out=1.0_wp)
         ispd = find_var(state, "ice_speed")
         iu = find_var(state, "ice_u")
         iv = find_var(state, "ice_v")
         call check(error, ispd > 0 .and. iu > 0 .and. iv > 0, &
                    "ice_speed/ice_u/ice_v must all register")
         if (allocated(error)) exit checks

         call ocean_state_enter_data(state)
         call state%diag%step(state, dt=1.0_wp, t=1.0_wp)
         call ocean_state_exit_data(state)

         call check(error, abs(state%diag%vars(iu)%output_buffer(3, 2, 1) - EXP_U) < TOL, &
                    "ice_u(3,2) must be the two-face average 0.5*(0.15+0.20)=0.175")
         if (allocated(error)) exit checks
         call check(error, abs(state%diag%vars(iv)%output_buffer(3, 2, 1) - EXP_V) < TOL, &
                    "ice_v(3,2) must be the signed -0.2 (no stray abs)")
         if (allocated(error)) exit checks
         call check(error, abs(state%diag%vars(ispd)%output_buffer(3, 2, 1) - EXP_SPD) < TOL, &
                    "ice_speed(3,2) must be sqrt(0.175^2 + 0.2^2)")
      end block checks
      call state%destroy()
   end subroutine test_ice_speed_face_average

   subroutine test_ice_diags_land_mask(error)
      !! A land cell (wet_T=0) reads 0 for both conc + thick; neighbours
      !! are unaffected.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      integer :: iconc, ithick
      checks: block
         call setup_state_ice(grid, state, ncat=3)
         state%ice%part_size(:, :, 1) = 0.3_wp
         state%ice%part_size(:, :, 2) = 0.4_wp
         state%ice%m_ice(:, :, 1) = 500.0_wp
         state%ice%m_ice(:, :, 2) = 800.0_wp

         ! Land-mask stamp AFTER metrics_finalize, BEFORE enter_data.
         state%metrics%wet_T(2, 2) = 0.0_wp

         call register_default_diags(state, dt_out=1.0_wp)
         iconc = find_var(state, "ice_conc")
         ithick = find_var(state, "ice_thick")
         call check(error, iconc > 0 .and. ithick > 0, "ice diags must be registered")
         if (allocated(error)) exit checks

         call ocean_state_enter_data(state)
         call state%diag%step(state, dt=1.0_wp, t=1.0_wp)
         call ocean_state_exit_data(state)

         call check(error, abs(state%diag%vars(iconc)%output_buffer(2, 2, 1)) < TOL, &
                    "land cell ice_conc should be 0")
         if (allocated(error)) exit checks
         call check(error, abs(state%diag%vars(ithick)%output_buffer(2, 2, 1)) < TOL, &
                    "land cell ice_thick should be 0")
         if (allocated(error)) exit checks
         call check(error, abs(state%diag%vars(iconc)%output_buffer(3, 2, 1) - 0.7_wp) < TOL, &
                    "neighbouring wet cell ice_conc should be unaffected (0.7)")
      end block checks
      call state%destroy()
   end subroutine test_ice_diags_land_mask

   subroutine test_ice_fill_matches_concentration_impl(error)
      !! Anti-divergence gate: a non-uniform part_size/m_ice pattern run
      !! through the diag-manager fills must match the host-side
      !! `ice_cell_concentration_impl` (the canonical two-mode gather) to
      !! FP round-off, so the fill copy cannot silently diverge from the
      !! EVP/coupler view.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      integer :: iconc, ithick, i, j, c, nx, ny, ncat
      real(wp), allocatable :: m_snow(:, :, :), mis(:, :), mice(:, :), ci(:, :)
      checks: block
         ncat = 3
         call setup_state_ice(grid, state, ncat=ncat)
         nx = grid%nx_total
         ny = grid%ny_total
         do c = 1, ncat
            do j = 1, ny
               do i = 1, nx
                  state%ice%part_size(i, j, c) = 0.1_wp*real(c, wp) + 0.01_wp*real(i, wp)
                  state%ice%m_ice(i, j, c) = 50.0_wp*real(c, wp) + real(j, wp)
               end do
            end do
         end do

         call register_default_diags(state, dt_out=1.0_wp)
         iconc = find_var(state, "ice_conc")
         ithick = find_var(state, "ice_thick")
         call check(error, iconc > 0 .and. ithick > 0, "ice diags must be registered")
         if (allocated(error)) exit checks

         call ocean_state_enter_data(state)
         call state%diag%step(state, dt=1.0_wp, t=1.0_wp)
         call ocean_state_exit_data(state)

         allocate (m_snow(nx, ny, ncat), source=0.0_wp)
         allocate (mis(nx, ny), mice(nx, ny), ci(nx, ny))
         call ice_cell_concentration_impl(state%metrics%wet_T, state%ice%part_size, &
                                          state%ice%m_ice, m_snow, mis, mice, ci, &
                                          ncat, nx, ny)

         call check(error, &
                    maxval(abs(state%diag%vars(iconc)%output_buffer(:, :, 1) - ci)) < TOL, &
                    "ice_conc fill must match ice_cell_concentration_impl's ci")
         if (allocated(error)) exit checks
         call check(error, &
                    maxval(abs(state%diag%vars(ithick)%output_buffer(:, :, 1) &
                               - mice/ICE_RHO_ICE)) < TOL, &
                    "ice_thick fill must match ice_cell_concentration_impl's mice/ICE_RHO_ICE")
      end block checks
      call state%destroy()
   end subroutine test_ice_fill_matches_concentration_impl

   subroutine test_console_ice_totals_reference(error)
      !! Pins the THIRD gather copy — the console-side `compute_ice_totals`
      !! (`rdb_ocean_console_stats`) — against a hand-computed reference, so
      !! it cannot silently drift from the fills / `ice_cell_concentration_impl`
      !! convention (the anti-divergence fill test only covers the fills copy).
      !! Reuses the supersaturated column (Σ part_size = 1.3 ⇒ ci clamps to 1)
      !! with a known mice.  `compute_ice_totals` reads its inputs on-DEVICE
      !! (!$acc reductions), so bracket the call with enter/exit_data; the out
      !! scalars are host-written by the reduction.  N_phys / areaT are derived
      !! from the grid so the goldens self-adjust if the constants change.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      integer :: nphys
      real(wp) :: area, wet_area, ci_area, hi_area
      real(wp) :: expected_wet, expected_ci, expected_hi
      real(wp), parameter :: MICE = 0.7_wp*905.0_wp + 0.6_wp*905.0_wp  ! 1176.5
      checks: block
         call setup_state_ice(grid, state, ncat=3)
         state%ice%part_size(:, :, 1) = 0.7_wp
         state%ice%part_size(:, :, 2) = 0.6_wp
         state%ice%m_ice(:, :, 1) = 905.0_wp
         state%ice%m_ice(:, :, 2) = 905.0_wp

         nphys = (grid%nx_total - 2*grid%nghost)*(grid%ny_total - 2*grid%nghost)
         area = state%metrics%areaT(2, 2)
         expected_wet = real(nphys, wp)*area
         expected_ci = 1.0_wp*real(nphys, wp)*area          ! ci clamped to 1
         expected_hi = (MICE/ICE_RHO_ICE)*real(nphys, wp)*area

         call ocean_state_enter_data(state)
         call compute_ice_totals(state%metrics%wet_T, state%metrics%areaT, &
                                 state%ice%part_size, state%ice%m_ice, 3, grid%nghost, &
                                 wet_area, ci_area, hi_area)
         call ocean_state_exit_data(state)

         call check(error, abs(wet_area - expected_wet) < TOL, &
                    "compute_ice_totals wet_area = N_phys * areaT")
         if (allocated(error)) exit checks
         call check(error, abs(ci_area - expected_ci) < TOL, &
                    "compute_ice_totals ci_area = clamped-1 conc * N_phys * areaT")
         if (allocated(error)) exit checks
         call check(error, abs(hi_area - expected_hi) < TOL, &
                    "compute_ice_totals hi_area = (mice/ICE_RHO_ICE) * N_phys * areaT")
      end block checks
      call state%destroy()
   end subroutine test_console_ice_totals_reference

   subroutine test_ice_off_no_ice_diags(error)
      !! Registry-level bit-identity: ice OFF ⇒ no ice_conc / ice_thick
      !! registered. (Membership contract only — the total var COUNT is
      !! deliberately not asserted; it couples this test to the unrelated
      !! size of the default diag set.)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      checks: block
         call setup_state(grid, state)
         call register_default_diags(state, dt_out=3600.0_wp)

         call check(error, find_var(state, "ice_conc") == 0, "ice_conc absent when ice is off")
         if (allocated(error)) exit checks
         call check(error, find_var(state, "ice_thick") == 0, "ice_thick absent when ice is off")
      end block checks
      call state%destroy()
   end subroutine test_ice_off_no_ice_diags

   subroutine test_ice_canonical_names(error)
      !! ice_conc/ice_thick are canonical; ice_speed routes to the derived
      !! catalog (not canonical).
      type(error_type), allocatable, intent(out) :: error
      checks: block
         call check(error, is_canonical_diag_name("ice_conc"), "ice_conc is canonical")
         if (allocated(error)) exit checks
         call check(error, is_canonical_diag_name("ice_thick"), "ice_thick is canonical")
         if (allocated(error)) exit checks
         call check(error,.not. is_canonical_diag_name("ice_speed"), &
                    "ice_speed is derived, not canonical")
      end block checks
   end subroutine test_ice_canonical_names

   subroutine test_ice_derived_catalog_has_velocity(error)
      !! The derived catalog LISTS ice_speed/ice_u/ice_v. (Membership only
      !! — the total catalog SIZE is deliberately not asserted; it couples
      !! this test to the unrelated count of derived diagnostics.)
      type(error_type), allocatable, intent(out) :: error
      integer :: n, i
      logical :: has_speed, has_u, has_v
      checks: block
         n = derived_catalog_size()
         has_speed = .false.
         has_u = .false.
         has_v = .false.
         do i = 1, n
            if (trim(derived_catalog_name(i)) == "ice_speed") has_speed = .true.
            if (trim(derived_catalog_name(i)) == "ice_u") has_u = .true.
            if (trim(derived_catalog_name(i)) == "ice_v") has_v = .true.
         end do
         call check(error, has_speed .and. has_u .and. has_v, &
                    "catalog must list ice_speed, ice_u, ice_v")
      end block checks
   end subroutine test_ice_derived_catalog_has_velocity

   subroutine test_ice_speed_zero_when_ice_off(error)
      !! ice_speed named on an ice-off run gets a zero field, not a crash
      !! (the is_init shim guard).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_state_t) :: state
      integer :: ispd
      checks: block
         call setup_state(grid, state)
         call register_derived(state, "ice_speed", time_op=DIAG_OP_INSTANT, dt_out=1.0_wp)
         ispd = find_var(state, "ice_speed")
         call check(error, ispd > 0, "ice_speed should register even with ice off")
         if (allocated(error)) exit checks

         call ocean_state_enter_data(state)
         call state%diag%step(state, dt=1.0_wp, t=1.0_wp)
         call ocean_state_exit_data(state)

         call check(error, maxval(abs(state%diag%vars(ispd)%output_buffer)) < TOL, &
                    "ice_speed buffer should be all zeros with ice off")
      end block checks
      call state%destroy()
   end subroutine test_ice_speed_zero_when_ice_off

end module test_ocean_ice_diags
