!! Regression test for the open-edge barotropic-η ghost fill
!! (rdb_barotropic_substep: x-pass / y-pass zero-gradient extrapolation).
!!
!! Bug (Stage-1 OBC bisection, 2026-06-18): open-edge η ghost cells were
!! never filled (only PERIODIC + the tripolar fold wrap them).  At an OPEN
!! edge the η ghost evolved under Pass-1's zeroed array-edge fluxes, seeding
!! an undamped 2Δx checkerboard null-mode in the ghost row/column.  `bebt`
!! damps only the physical-face velocity update and never reaches the ghost
!! η, so with the multilayer split re-exciting it each outer step it grew
!! ~1-step e-fold and NaN'd within ~12 steps over non-trivial bathymetry
!! (flat-bottom is stable — the topographic H-variation is the excitation;
!! 1-layer is stable — the baroclinic split is what feeds the ghost mode).
!!
!! Test (deterministic invariant): seed every per-layer velocity GHOST face
!! (edges AND ghost×ghost corners) with a large sentinel, zero the physical
!! faces, apply the open boundary once, and assert the ghosts are now
!! zero-gradient-filled (≈0) rather than retaining the sentinel.  This
!! exercises `fill_uv_layer_ghosts` — the corner velocity fill whose absence
!! is the fatal open-edge NaN.  The full nonlinear blow-up only triggers at
!! basin scale (it is config-marginal in a small clean box), so it is covered
!! end-to-end by the spoon-open validation run rather than this unit test; the
!! η-corner fill in the barotropic substep uses the identical x-then-y
!! corner-safe pattern verified here.
module test_ocean_obc_eta_ghost
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_continuity, only: continuity_t
   use rdb_coriolis_adv, only: coriolis_adv_t
   use rdb_eos, only: eos_t
   use rdb_ocean_pressure_force, only: ocean_pressure_force_t, OPGF_VARIANT_FV_MOM6
   use rdb_ocean_horizontal_viscosity, only: ocean_horizontal_viscosity_t
   use rdb_ocean_bottom_drag, only: ocean_bottom_drag_t
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t
   use rdb_ocean_vertical_advection, only: ocean_vertical_advection_t
   use rdb_ocean_hdiff_tracer, only: ocean_hdiff_tracer_t
   use rdb_ocean_vdiff, only: ocean_vdiff_t
   use rdb_ocean_vmix, only: ocean_vmix_t
   use rdb_ocean_dyn, only: ocean_dyn_t, ocean_dyn_step_split
   use rdb_ocean_obc_baroclinic, only: ocean_obc_apply_baroclinic, ocean_obc_fill_ghosts, &
                                       ocean_obc_refill_ghost_ssh
   use rdb_ocean_boundary_types, only: ocean_bc_state_t, ocean_bc_state_init, &
                                       ocean_bc_state_destroy, OBC_OPEN, OBC_WALL, OBC_PERIODIC
   implicit none
   private

   public :: collect_ocean_obc_eta_ghost_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 4

contains

   subroutine collect_ocean_obc_eta_ghost_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("obc_open_corner_ghost_fill", test_open_corner_ghost_fill), &
                  new_unittest("obc_wallopen_corner_ghost_fill", test_wallopen_corner_ghost_fill), &
                  new_unittest("obc_open_corner_tracer_fill", test_open_corner_tracer_fill), &
                  new_unittest("obc_refill_ghost_ssh_bstep", test_refill_ghost_ssh_bstep), &
                  new_unittest("obc_refill_ghost_ssh_flatb", test_refill_ghost_ssh_flatb), &
                  new_unittest("obc_refill_ghost_ssh_noop", test_refill_ghost_ssh_noop) &
                  ]
   end subroutine collect_ocean_obc_eta_ghost_tests

   ! -----------------------------------------------------------------
   ! Setup helpers (mirror test_ocean_periodic.F90)
   ! -----------------------------------------------------------------

   subroutine init_all(grid, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn)
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      type(continuity_t), intent(inout) :: ct
      type(coriolis_adv_t), intent(inout) :: cor
      type(ocean_pressure_force_t), intent(inout) :: pgf
      type(ocean_horizontal_viscosity_t), intent(inout) :: hv
      type(ocean_bottom_drag_t), intent(inout) :: bd
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_vertical_advection_t), intent(inout) :: va
      type(ocean_hdiff_tracer_t), intent(inout) :: hd
      type(ocean_vdiff_t), intent(inout) :: vd
      type(ocean_vmix_t), intent(inout) :: vmix
      type(eos_t), intent(inout) :: eos
      type(ocean_dyn_t), intent(inout) :: dyn
      ms%nz_ml = NZ
      call ms%init(grid)
      call ct%init(grid, nz_ml=NZ)
      call cor%init(grid, nz_ml=NZ)
      call pgf%init(grid, nz_ml=NZ)
      call hv%init(grid, nz_ml=NZ)
      call bd%init(grid, nz_ml=NZ)
      call ss%init(grid, nz_ml=NZ)
      call va%init(grid, nz_ml=NZ)
      call hd%init(grid, nz_ml=NZ)
      call vd%init(grid, nz_ml=NZ)
      call vmix%init(grid, nz_ml=NZ)
      call eos%init(grid)
      call dyn%init(grid, nz_ml=NZ)
   end subroutine init_all

   subroutine map_in(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
      type(hgrid_t), intent(in) :: grid
      type(ocean_metrics_t), intent(inout) :: metrics
      type(multilayer_state_t), intent(inout) :: ms
      type(continuity_t), intent(inout) :: ct
      type(coriolis_adv_t), intent(inout) :: cor
      type(ocean_pressure_force_t), intent(inout) :: pgf
      type(ocean_horizontal_viscosity_t), intent(inout) :: hv
      type(ocean_bottom_drag_t), intent(inout) :: bd
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_vertical_advection_t), intent(inout) :: va
      type(ocean_hdiff_tracer_t), intent(inout) :: hd
      type(ocean_vdiff_t), intent(inout) :: vd
      type(ocean_vmix_t), intent(inout) :: vmix
      type(ocean_dyn_t), intent(inout) :: dyn
      call make_cartesian_metrics(metrics, grid)
      !$acc enter data copyin(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
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
   end subroutine map_in

   subroutine map_out(metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
      type(ocean_metrics_t), intent(inout) :: metrics
      type(multilayer_state_t), intent(inout) :: ms
      type(continuity_t), intent(inout) :: ct
      type(coriolis_adv_t), intent(inout) :: cor
      type(ocean_pressure_force_t), intent(inout) :: pgf
      type(ocean_horizontal_viscosity_t), intent(inout) :: hv
      type(ocean_bottom_drag_t), intent(inout) :: bd
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_vertical_advection_t), intent(inout) :: va
      type(ocean_hdiff_tracer_t), intent(inout) :: hd
      type(ocean_vdiff_t), intent(inout) :: vd
      type(ocean_vmix_t), intent(inout) :: vmix
      type(ocean_dyn_t), intent(inout) :: dyn
      call destroy_cartesian_metrics(metrics)
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
      !$acc exit data delete(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
   end subroutine map_out

   subroutine destroy_all(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn)
      type(multilayer_state_t), intent(inout) :: ms
      type(continuity_t), intent(inout) :: ct
      type(coriolis_adv_t), intent(inout) :: cor
      type(ocean_pressure_force_t), intent(inout) :: pgf
      type(ocean_horizontal_viscosity_t), intent(inout) :: hv
      type(ocean_bottom_drag_t), intent(inout) :: bd
      type(ocean_surface_stress_t), intent(inout) :: ss
      type(ocean_vertical_advection_t), intent(inout) :: va
      type(ocean_hdiff_tracer_t), intent(inout) :: hd
      type(ocean_vdiff_t), intent(inout) :: vd
      type(ocean_vmix_t), intent(inout) :: vmix
      type(eos_t), intent(inout) :: eos
      type(ocean_dyn_t), intent(inout) :: dyn
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
   end subroutine destroy_all

   ! -----------------------------------------------------------------
   ! Test: open-edge per-layer VELOCITY ghost fill covers the CORNERS
   ! -----------------------------------------------------------------

   subroutine test_open_corner_ghost_fill(error)
      !! Invariant test of the fatal open-edge gap: the per-layer face
      !! velocities at the ghost×ghost CORNERS (and edge ghosts) must be
      !! zero-gradient-filled by the OBC apply.  Left unfilled, those corner
      !! velocities are read by the C-grid Coriolis/KE/vorticity stencil at the
      !! corner-adjacent PHYSICAL cell and blow the column up (the open-edge
      !! NaN; reproduced end-to-end by the spoon-open validation run).
      !!
      !! Deterministic construction: seed EVERY ghost face with a large
      !! sentinel and zero only the physical faces, then apply the OBC once.
      !! With the fill (fill_uv_layer_ghosts, corners included) the at-rest
      !! interior makes every ghost zero-gradient to ~0; without it the ghost
      !! faces — corners especially — keep the sentinel.  A whole-array max
      !! therefore discriminates fixed (≈0) from unfixed (= sentinel) cleanly,
      !! without depending on a marginal nonlinear blow-up that only triggers
      !! at basin scale.
      type(error_type), allocatable, intent(out) :: error
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
      type(eos_t) :: eos
      type(ocean_dyn_t) :: dyn
      type(ocean_bc_state_t) :: bc

      real(wp), parameter :: H0 = 1000.0_wp        ! flat total depth (m)
      real(wp), parameter :: SENTINEL = 1.0e3_wp    ! ghost-seed marker (m/s)
      real(wp), parameter :: DX = 5000.0_wp
      real(wp), parameter :: DT = 300.0_wp
      integer :: nx_p, ny_p, k, ig, i0, i1, j0, j1
      real(wp) :: max_u, max_v

      checks: block

         nx_p = 24
         ny_p = 20
         call grid%init(nx_p, ny_p, NGHOST, DX, DX)
         call init_all(grid, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn)

         ig = grid%nghost
         i0 = ig + 1; i1 = ig + grid%nx_phys     ! first / last physical CELL
         j0 = ig + 1; j1 = ig + grid%ny_phys

         ! Flat, at-rest base state.
         ms%h_layer = H0/real(NZ, wp)
         dyn%bt_work%bt_H_ref = H0
         dyn%bt_work%bt_eta = 0.0_wp
         dyn%bt_work%bt_ubt_end = 0.0_wp         ! anomaly scheme reads these
         dyn%bt_work%bt_vbt_end = 0.0_wp
         do k = 1, NZ
            ms%tracers(ms%idx_temperature)%hTr(:, :, k) = 10.0_wp*ms%h_layer(:, :, k)
            ms%tracers(ms%idx_salinity)%hTr(:, :, k) = eos%S_ref*ms%h_layer(:, :, k)
         end do

         ! Seed every face with the sentinel, then zero the PHYSICAL faces, so
         ! only the ghost faces (edges AND corners) carry it.  Physical u faces
         ! span the wall-to-wall range i0..i1+1; physical v faces j0..j1+1.
         ms%u_face_x_layer = SENTINEL
         ms%v_face_y_layer = SENTINEL
         ms%u_face_x_layer(i0:i1 + 1, j0:j1, :) = 0.0_wp
         ms%v_face_y_layer(i0:i1, j0:j1 + 1, :) = 0.0_wp

         ! All four edges OPEN.
         call ocean_bc_state_init(bc, grid, nz_ml=NZ, n_tracers=2)
         bc%west%bc_type = OBC_OPEN
         bc%east%bc_type = OBC_OPEN
         bc%south%bc_type = OBC_OPEN
         bc%north%bc_type = OBC_OPEN

         call map_in(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
         !$acc update device(dyn%bt_work%bt_H_ref, dyn%bt_work%bt_eta, &
         !$acc&              dyn%bt_work%bt_ubt_end, dyn%bt_work%bt_vbt_end, &
         !$acc&              ms%u_face_x_layer, ms%v_face_y_layer)

         ! Apply the open boundary once — sets the physical edge faces and
         ! (the fix) zero-gradient-fills the velocity ghosts, corners included.
         call ocean_obc_apply_baroclinic(grid, bc, dyn%bt_work, ms, DT)

         !$acc update self(ms%u_face_x_layer, ms%v_face_y_layer)
         max_u = maxval(abs(ms%u_face_x_layer))   ! full array — ghosts + corners
         max_v = maxval(abs(ms%v_face_y_layer))

         call map_out(metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)

         ! With the fill every ghost (incl. corners) is zero-gradient ⇒ ~0;
         ! without it the seeded sentinel survives in the ghost faces.
         call check(error, max_u < 1.0_wp, &
                    "open-edge u ghost (incl. corners) left unfilled — max|u|="// &
                    trim(adjusted_str(max_u))//" (seed was "//trim(adjusted_str(SENTINEL))//")")
         if (allocated(error)) exit checks
         call check(error, max_v < 1.0_wp, &
                    "open-edge v ghost (incl. corners) left unfilled — max|v|="// &
                    trim(adjusted_str(max_v))//" (seed was "//trim(adjusted_str(SENTINEL))//")")

      end block checks

      call ocean_bc_state_destroy(bc)
      call destroy_all(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn)
   end subroutine test_open_corner_ghost_fill

   ! -----------------------------------------------------------------
   ! Test: WALL×OPEN corner — a south-open edge must fill its velocity
   ! ghost CORNER even against a WALL side (west/east), not just open×open.
   ! -----------------------------------------------------------------
   subroutine test_wallopen_corner_ghost_fill(error)
      !! Only SOUTH is open; west/east/north are WALL.  Seed every velocity
      !! ghost with a sentinel; apply the OBC.  The south-open edge must
      !! zero-gradient-fill its ghost rows over the FULL i extent — including
      !! the wall-side corner columns — because a non-periodic (wall) x-edge
      !! is no longer excluded from the corner fill.  Without that, the SW/SE
      !! ghost corners keep the sentinel (the day-20 EAC velocity NaN).
      type(error_type), allocatable, intent(out) :: error
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
      type(eos_t) :: eos
      type(ocean_dyn_t) :: dyn
      type(ocean_bc_state_t) :: bc
      real(wp), parameter :: H0 = 1000.0_wp, SENTINEL = 1.0e3_wp, DX = 5000.0_wp, DT = 300.0_wp
      integer :: nx_p, ny_p, k, ig, i0, i1, j0, j1
      real(wp) :: max_u, max_v

      checks: block
         nx_p = 24; ny_p = 20
         call grid%init(nx_p, ny_p, NGHOST, DX, DX)
         call init_all(grid, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn)
         ig = grid%nghost
         i0 = ig + 1; i1 = ig + grid%nx_phys
         j0 = ig + 1; j1 = ig + grid%ny_phys
         ms%h_layer = H0/real(NZ, wp)
         dyn%bt_work%bt_H_ref = H0
         dyn%bt_work%bt_eta = 0.0_wp
         dyn%bt_work%bt_ubt_end = 0.0_wp
         dyn%bt_work%bt_vbt_end = 0.0_wp
         do k = 1, NZ
            ms%tracers(ms%idx_temperature)%hTr(:, :, k) = 10.0_wp*ms%h_layer(:, :, k)
            ms%tracers(ms%idx_salinity)%hTr(:, :, k) = eos%S_ref*ms%h_layer(:, :, k)
         end do
         ! Seed ONLY the south ghost rows (incl. the wall-side corner columns)
         ! with the sentinel; the y-pass source row (first physical row) stays
         ! 0, so a correct south-open fill drives the whole south ghost to 0.
         ! (Seeding the wall-ghost physical rows too would be unfair — the wall
         ! BC, absent in this isolated call, normally keeps that source sane.)
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         ms%u_face_x_layer(:, 1:ig, :) = SENTINEL
         ms%v_face_y_layer(:, 1:ig, :) = SENTINEL

         ! Only SOUTH open; the other three are WALL.
         call ocean_bc_state_init(bc, grid, nz_ml=NZ, n_tracers=2)
         bc%west%bc_type = OBC_WALL
         bc%east%bc_type = OBC_WALL
         bc%north%bc_type = OBC_WALL
         bc%south%bc_type = OBC_OPEN

         call map_in(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
         !$acc update device(dyn%bt_work%bt_H_ref, dyn%bt_work%bt_eta, &
         !$acc&              dyn%bt_work%bt_ubt_end, dyn%bt_work%bt_vbt_end, &
         !$acc&              ms%u_face_x_layer, ms%v_face_y_layer)
         call ocean_obc_apply_baroclinic(grid, bc, dyn%bt_work, ms, DT)
         !$acc update self(ms%u_face_x_layer, ms%v_face_y_layer)

         ! Check the SOUTH ghost rows over the FULL i extent (includes the
         ! wall-side corner columns).  With the wall×open fix these are
         ! zero-gradient-filled (~0); without it the corner columns keep the
         ! sentinel.
         max_u = maxval(abs(ms%u_face_x_layer(:, 1:ig, :)))
         max_v = maxval(abs(ms%v_face_y_layer(:, 1:ig, :)))
         call map_out(metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)

         call check(error, max_u < 1.0_wp, &
                    "wall×open: south-ghost u corner left unfilled — max|u|="// &
                    trim(adjusted_str(max_u)))
         if (allocated(error)) exit checks
         call check(error, max_v < 1.0_wp, &
                    "wall×open: south-ghost v corner left unfilled — max|v|="// &
                    trim(adjusted_str(max_v)))
      end block checks

      call ocean_bc_state_destroy(bc)
      call destroy_all(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn)
   end subroutine test_wallopen_corner_ghost_fill

   ! -----------------------------------------------------------------
   ! Test: tracer hTr ghost CORNERS get zero-gradient-filled at open edges.
   ! -----------------------------------------------------------------
   subroutine test_open_corner_tracer_fill(error)
      !! Seed every tracer ghost (incl. corners) with a large sentinel
      !! CONCENTRATION; apply the ghost fill.  The zero-gradient pre-fill must
      !! reset every open-edge ghost concentration — corners included — to the
      !! interior value, so the upwind per-edge fill leaves the corners
      !! zero-gradient rather than the sentinel (the inflow-driven day-2 T/S
      !! NaN).  At rest (u=v=0) every open edge is outflow ⇒ all ghosts revert
      !! to the interior concentration.
      type(error_type), allocatable, intent(out) :: error
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
      type(eos_t) :: eos
      type(ocean_dyn_t) :: dyn
      type(ocean_bc_state_t) :: bc
      real(wp), parameter :: H0 = 1000.0_wp, S_INT = 35.0_wp, SENTINEL = 999.0_wp, DX = 5000.0_wp
      integer :: nx_p, ny_p, k, ig, i0, i1, j0, j1, it
      real(wp) :: hlay, max_conc

      checks: block
         nx_p = 24; ny_p = 20
         call grid%init(nx_p, ny_p, NGHOST, DX, DX)
         call init_all(grid, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn)
         ig = grid%nghost
         i0 = ig + 1; i1 = ig + grid%nx_phys
         j0 = ig + 1; j1 = ig + grid%ny_phys
         hlay = H0/real(NZ, wp)
         ms%h_layer = hlay
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         it = ms%idx_salinity
         ! Seed the WHOLE salinity field with the sentinel concentration, then
         ! overwrite the PHYSICAL interior with the real concentration, so only
         ! the ghosts (edges + corners) carry the sentinel.
         do k = 1, NZ
            ms%tracers(it)%hTr(:, :, k) = SENTINEL*hlay
            ms%tracers(it)%hTr(i0:i1, j0:j1, k) = S_INT*hlay
            ms%tracers(ms%idx_temperature)%hTr(:, :, k) = 10.0_wp*hlay
         end do

         call ocean_bc_state_init(bc, grid, nz_ml=NZ, n_tracers=2)
         bc%west%bc_type = OBC_OPEN
         bc%east%bc_type = OBC_OPEN
         bc%south%bc_type = OBC_OPEN
         bc%north%bc_type = OBC_OPEN

         call map_in(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
         !$acc update device(ms%h_layer, ms%u_face_x_layer, ms%v_face_y_layer, &
         !$acc&              ms%tracers(it)%hTr)
         call ocean_obc_fill_ghosts(grid, bc, ms)
         !$acc update self(ms%tracers(it)%hTr, ms%h_layer)

         ! Concentration (hTr/h) over the full array — corners included.  With
         ! the corner pre-fill every ghost reverts to ~S_INT; without it the
         ! ghost corners keep the sentinel concentration.
         max_conc = maxval(ms%tracers(it)%hTr/ms%h_layer)
         call map_out(metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)

         call check(error, max_conc < 100.0_wp, &
                    "tracer ghost corner left unfilled — max concentration="// &
                    trim(adjusted_str(max_conc))//" (interior "// &
                    trim(adjusted_str(S_INT))//", seed "//trim(adjusted_str(SENTINEL))//")")
      end block checks

      call ocean_bc_state_destroy(bc)
      call destroy_all(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn)
   end subroutine test_open_corner_tracer_fill

   ! -----------------------------------------------------------------
   ! Tests: end-of-step open-edge ghost SSH refill (ocean_obc_refill_ghost_ssh)
   !
   ! Bug (open-boundary SSH blow-up, 2026-07): the slow continuity updates
   ! the GHOST h_layer over the full array with zeroed array-edge fluxes and
   ! the conservative ALE remap preserves that drift, so the diagnosed
   ! `SSH = Σh_layer − b` grew to tens of metres in the halo (worst at open
   ! corners) while the physical interior stayed healthy.  The fix rescales
   ! each open-edge ghost column at end of outer step so
   ! `Σh_ghost = b_ghost + η_interior` (zero-gradient free surface), which
   ! also absorbs any static `b_interior − b_ghost` step from formula
   ! topographies (spoon/seamount clamp their ghost b instead of
   ! constant-extrapolating like the file loader).
   ! -----------------------------------------------------------------

   subroutine test_refill_ghost_ssh_bstep(error)
      !! Formula-topo static-step + injected drift case.  Spatially-varying
      !! bathymetry where every ghost `b` differs from its first interior
      !! cell; interior seeded with a KNOWN uniform free surface η = ETA0;
      !! open-edge ghost h_layer columns then corrupted with a drift sentinel
      !! (simulating the continuity/remap ghost drift).  After
      !! `ocean_obc_refill_ghost_ssh`, EVERY open-edge ghost cell — including
      !! the open×open corner ghosts — must satisfy
      !! `Σ_k h_ghost − b_ghost == ETA0` to round-off, and the tracer ghost
      !! CONCENTRATION must be zero-gradient (== interior C0).
      type(error_type), allocatable, intent(out) :: error
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
      type(eos_t) :: eos
      type(ocean_dyn_t) :: dyn
      type(ocean_bc_state_t) :: bc
      real(wp), parameter :: ETA0 = 0.5_wp        ! known interior free surface (m)
      real(wp), parameter :: C0 = 35.0_wp         ! interior tracer concentration
      real(wp), parameter :: DRIFT_H = 3.0_wp     ! ghost drift sentinel (m / layer)
      real(wp), parameter :: DX = 5000.0_wp
      real(wp), parameter :: TOL = 1.0e-8_wp
      integer :: nx_p, ny_p, nxt, nyt, i, j, k, ig, it
      real(wp) :: ssh_err_we, ssh_err_s, conc_err, ssh, conc

      checks: block
         nx_p = 24; ny_p = 20
         call grid%init(nx_p, ny_p, NGHOST, DX, DX)
         call init_all(grid, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn)
         ig = grid%nghost
         nxt = grid%nx_total; nyt = grid%ny_total
         it = ms%idx_salinity

         ! Spatially-varying bathymetry over the FULL array (ghosts included):
         ! every ghost b differs from the first interior cell — the formula-
         ! topography static-step case (spoon clamps its ghost b).
         do j = 1, nyt
            do i = 1, nxt
               dyn%bt_work%bt_H_ref(i, j) = 500.0_wp + 10.0_wp*real(i, wp) + 7.0_wp*real(j, wp)
            end do
         end do
         ! Consistent state everywhere: Σh = b + ETA0 (uniform free surface),
         ! uniform layer split, tracer concentration C0.
         do k = 1, NZ
            ms%h_layer(:, :, k) = (dyn%bt_work%bt_H_ref + ETA0)/real(NZ, wp)
            ms%tracers(it)%hTr(:, :, k) = C0*ms%h_layer(:, :, k)
            ms%tracers(ms%idx_temperature)%hTr(:, :, k) = 10.0_wp*ms%h_layer(:, :, k)
         end do
         ! Inject the drift sentinel into the OPEN-edge ghost regions only
         ! (W/E ghost columns over the full j extent, S ghost rows over the
         ! full i extent — the exact region the end-of-step refill owns).
         ms%h_layer(1:ig, :, :) = DRIFT_H
         ms%h_layer(nxt - ig + 1:nxt, :, :) = DRIFT_H
         ms%h_layer(:, 1:ig, :) = DRIFT_H
         ms%tracers(it)%hTr(1:ig, :, :) = 999.0_wp*DRIFT_H
         ms%tracers(it)%hTr(nxt - ig + 1:nxt, :, :) = 999.0_wp*DRIFT_H
         ms%tracers(it)%hTr(:, 1:ig, :) = 999.0_wp*DRIFT_H

         call ocean_bc_state_init(bc, grid, nz_ml=NZ, n_tracers=2)
         bc%west%bc_type = OBC_OPEN
         bc%east%bc_type = OBC_OPEN
         bc%south%bc_type = OBC_OPEN
         bc%north%bc_type = OBC_WALL

         call map_in(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
         !$acc update device(dyn%bt_work%bt_H_ref, ms%h_layer, ms%tracers(it)%hTr, &
         !$acc&              ms%tracers(ms%idx_temperature)%hTr)

         call ocean_obc_refill_ghost_ssh(grid, bc, ms, dyn%bt_work%bt_H_ref)

         !$acc update self(ms%h_layer, ms%tracers(it)%hTr)

         ! Max |SSH − ETA0| over the W/E ghost columns (full j extent — the
         ! open×open SW/SE corner ghosts and the wall-side N corner ghosts
         ! that the x-pass owns are all included).
         ssh_err_we = 0.0_wp
         do j = 1, nyt
            do i = 1, ig
               ssh = sum(ms%h_layer(i, j, :)) - dyn%bt_work%bt_H_ref(i, j)
               ssh_err_we = max(ssh_err_we, abs(ssh - ETA0))
               ssh = sum(ms%h_layer(nxt - i + 1, j, :)) - dyn%bt_work%bt_H_ref(nxt - i + 1, j)
               ssh_err_we = max(ssh_err_we, abs(ssh - ETA0))
            end do
         end do
         ! Max |SSH − ETA0| over the S ghost rows (full i extent — corners).
         ssh_err_s = 0.0_wp
         do j = 1, ig
            do i = 1, nxt
               ssh = sum(ms%h_layer(i, j, :)) - dyn%bt_work%bt_H_ref(i, j)
               ssh_err_s = max(ssh_err_s, abs(ssh - ETA0))
            end do
         end do
         ! Tracer ghost concentration must be zero-gradient == C0 everywhere
         ! in the same refill-owned ghost regions.
         conc_err = 0.0_wp
         do k = 1, NZ
            do j = 1, nyt
               do i = 1, ig
                  conc = ms%tracers(it)%hTr(i, j, k)/ms%h_layer(i, j, k)
                  conc_err = max(conc_err, abs(conc - C0))
                  conc = ms%tracers(it)%hTr(nxt - i + 1, j, k)/ms%h_layer(nxt - i + 1, j, k)
                  conc_err = max(conc_err, abs(conc - C0))
               end do
            end do
            do j = 1, ig
               do i = 1, nxt
                  conc = ms%tracers(it)%hTr(i, j, k)/ms%h_layer(i, j, k)
                  conc_err = max(conc_err, abs(conc - C0))
               end do
            end do
         end do

         call map_out(metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)

         call check(error, ssh_err_we < TOL, &
                    "refill: W/E ghost SSH /= eta_interior — max err="// &
                    trim(adjusted_str(ssh_err_we)))
         if (allocated(error)) exit checks
         call check(error, ssh_err_s < TOL, &
                    "refill: S ghost SSH /= eta_interior (corners incl.) — max err="// &
                    trim(adjusted_str(ssh_err_s)))
         if (allocated(error)) exit checks
         call check(error, conc_err < TOL, &
                    "refill: ghost tracer concentration not zero-gradient — max err="// &
                    trim(adjusted_str(conc_err)))
      end block checks

      call ocean_bc_state_destroy(bc)
      call destroy_all(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn)
   end subroutine test_refill_ghost_ssh_bstep

   subroutine test_refill_ghost_ssh_flatb(error)
      !! Scale-collapses-to-1 case: with CONSTANT-EXTRAPOLATED ghost `b`
      !! (the file-bathymetry loader convention: x-columns first, then
      !! y-rows, corners from the x-filled columns) the refill's column
      !! scale is 1 and the operation reduces to a plain zero-gradient copy
      !! of the interior column — per-layer, not just in the column sum.
      !! Non-uniform layer fractions make the per-layer assertion sharp.
      type(error_type), allocatable, intent(out) :: error
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
      type(eos_t) :: eos
      type(ocean_dyn_t) :: dyn
      type(ocean_bc_state_t) :: bc
      real(wp), parameter :: ETA0 = 0.25_wp
      real(wp), parameter :: DRIFT_H = 3.0_wp
      real(wp), parameter :: DX = 5000.0_wp
      real(wp), parameter :: RTOL = 1.0e-12_wp
      integer :: nx_p, ny_p, nxt, nyt, i, j, k, ig, i0, i1, j0, j1
      real(wp) :: frac, rel_err, href, hgot

      checks: block
         nx_p = 24; ny_p = 20
         call grid%init(nx_p, ny_p, NGHOST, DX, DX)
         call init_all(grid, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn)
         ig = grid%nghost
         nxt = grid%nx_total; nyt = grid%ny_total
         i0 = ig + 1; i1 = ig + grid%nx_phys
         j0 = ig + 1; j1 = ig + grid%ny_phys

         ! Varying interior bathymetry; ghost b by CONSTANT EXTRAPOLATION
         ! (x-columns first, then y-rows — same order as the file loader,
         ! so the corner ghosts copy the x-filled columns).
         do j = j0, j1
            do i = i0, i1
               dyn%bt_work%bt_H_ref(i, j) = 800.0_wp + 15.0_wp*real(i, wp) + 5.0_wp*real(j, wp)
            end do
         end do
         do j = 1, nyt
            do i = 1, ig
               dyn%bt_work%bt_H_ref(i, j) = dyn%bt_work%bt_H_ref(i0, min(max(j, j0), j1))
               dyn%bt_work%bt_H_ref(nxt - i + 1, j) = dyn%bt_work%bt_H_ref(i1, min(max(j, j0), j1))
            end do
         end do
         do j = 1, ig
            do i = 1, nxt
               dyn%bt_work%bt_H_ref(i, j) = dyn%bt_work%bt_H_ref(i, j0)
               dyn%bt_work%bt_H_ref(i, nyt - j + 1) = dyn%bt_work%bt_H_ref(i, j1)
            end do
         end do
         ! Non-uniform layer split: fraction f_k = 2k / (NZ (NZ+1)).
         do k = 1, NZ
            frac = 2.0_wp*real(k, wp)/(real(NZ, wp)*real(NZ + 1, wp))
            ms%h_layer(:, :, k) = frac*(dyn%bt_work%bt_H_ref + ETA0)
            ms%tracers(ms%idx_salinity)%hTr(:, :, k) = 35.0_wp*ms%h_layer(:, :, k)
            ms%tracers(ms%idx_temperature)%hTr(:, :, k) = 10.0_wp*ms%h_layer(:, :, k)
         end do
         ! Drift sentinel in the open ghost regions.
         ms%h_layer(1:ig, :, :) = DRIFT_H
         ms%h_layer(nxt - ig + 1:nxt, :, :) = DRIFT_H
         ms%h_layer(:, 1:ig, :) = DRIFT_H

         call ocean_bc_state_init(bc, grid, nz_ml=NZ, n_tracers=2)
         bc%west%bc_type = OBC_OPEN
         bc%east%bc_type = OBC_OPEN
         bc%south%bc_type = OBC_OPEN
         bc%north%bc_type = OBC_WALL

         call map_in(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
         !$acc update device(dyn%bt_work%bt_H_ref, ms%h_layer, &
         !$acc&              ms%tracers(ms%idx_salinity)%hTr, &
         !$acc&              ms%tracers(ms%idx_temperature)%hTr)

         call ocean_obc_refill_ghost_ssh(grid, bc, ms, dyn%bt_work%bt_H_ref)

         !$acc update self(ms%h_layer)

         ! Per-layer zero-gradient copy: every W/E ghost layer equals the
         ! first interior cell of its row; every S ghost layer equals the
         ! first interior row of its (x-refilled) column.
         rel_err = 0.0_wp
         do k = 1, NZ
            do j = 1, nyt
               do i = 1, ig
                  href = ms%h_layer(i0, j, k)
                  hgot = ms%h_layer(i, j, k)
                  rel_err = max(rel_err, abs(hgot - href)/href)
                  href = ms%h_layer(i1, j, k)
                  hgot = ms%h_layer(nxt - i + 1, j, k)
                  rel_err = max(rel_err, abs(hgot - href)/href)
               end do
            end do
            do j = 1, ig
               do i = 1, nxt
                  href = ms%h_layer(i, j0, k)
                  hgot = ms%h_layer(i, j, k)
                  rel_err = max(rel_err, abs(hgot - href)/href)
               end do
            end do
         end do

         call map_out(metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)

         call check(error, rel_err < RTOL, &
                    "refill with extrapolated ghost b must be a plain zero-gradient "// &
                    "copy — max per-layer rel err="//trim(adjusted_str(rel_err)))
      end block checks

      call ocean_bc_state_destroy(bc)
      call destroy_all(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn)
   end subroutine test_refill_ghost_ssh_flatb

   subroutine test_refill_ghost_ssh_noop(error)
      !! No-op guarantee: with WALL and PERIODIC edges only (no open-ish
      !! edge), the refill must leave h_layer BIT-IDENTICAL — the routine
      !! returns before launching any kernel, so wall/periodic configs stay
      !! byte-unchanged.
      type(error_type), allocatable, intent(out) :: error
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
      type(eos_t) :: eos
      type(ocean_dyn_t) :: dyn
      type(ocean_bc_state_t) :: bc
      real(wp), parameter :: DX = 5000.0_wp
      integer :: nx_p, ny_p, i, j, k
      real(wp), allocatable :: h_before(:, :, :)
      logical :: identical

      checks: block
         nx_p = 24; ny_p = 20
         call grid%init(nx_p, ny_p, NGHOST, DX, DX)
         call init_all(grid, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn)

         dyn%bt_work%bt_H_ref = 1000.0_wp
         ! Deliberately RAGGED h_layer (ghosts inconsistent with b) so any
         ! accidental write would be visible.
         do k = 1, NZ
            do j = 1, grid%ny_total
               do i = 1, grid%nx_total
                  ms%h_layer(i, j, k) = 100.0_wp + real(i, wp) + 0.3_wp*real(j, wp) + 7.0_wp*real(k, wp)
               end do
            end do
            ms%tracers(ms%idx_salinity)%hTr(:, :, k) = 35.0_wp*ms%h_layer(:, :, k)
            ms%tracers(ms%idx_temperature)%hTr(:, :, k) = 10.0_wp*ms%h_layer(:, :, k)
         end do
         h_before = ms%h_layer

         ! PERIODIC x + WALL y: no open-ish edge anywhere.
         call ocean_bc_state_init(bc, grid, nz_ml=NZ, n_tracers=2)
         bc%west%bc_type = OBC_PERIODIC
         bc%east%bc_type = OBC_PERIODIC
         bc%south%bc_type = OBC_WALL
         bc%north%bc_type = OBC_WALL

         call map_in(grid, metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)
         !$acc update device(dyn%bt_work%bt_H_ref, ms%h_layer)

         call ocean_obc_refill_ghost_ssh(grid, bc, ms, dyn%bt_work%bt_H_ref)

         !$acc update self(ms%h_layer)
         identical = all(ms%h_layer == h_before)
         call map_out(metrics, ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, dyn)

         call check(error, identical, &
                    "refill must be bit-identical no-op for WALL/PERIODIC edges")
      end block checks

      call ocean_bc_state_destroy(bc)
      call destroy_all(ms, ct, cor, pgf, hv, bd, ss, va, hd, vd, vmix, eos, dyn)
   end subroutine test_refill_ghost_ssh_noop

   pure function adjusted_str(x) result(s)
      real(wp), intent(in) :: x
      character(len=32) :: s
      write (s, '(es12.4)') x
      s = adjustl(s)
   end function adjusted_str

end module test_ocean_obc_eta_ghost
