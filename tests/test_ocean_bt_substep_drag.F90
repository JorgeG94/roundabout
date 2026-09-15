!! Unit tests for the BT-substep multiplicative drag damping
!! (MOM6 `bt_rem_u` / `bt_rem_v` at `MOM_barotropic.F90:3648`).
!!
!! The damping factor is computed by `compute_bt_rem` (in
!! `rdb_barotropic_coupling`) as
!!     bt_rem_face = Htot_face / (Htot_face + r·HBBL·dt_inner)
!! and is applied INSIDE the BT substep loop as
!!     ubt_new = bt_rem · (ubt_old + dt·forces)
!! Per substep this contributes a small drag; cumulative over many
!! substeps it gives an exponential-like decay of the BT mode.
!!
!! Tests:
!!   1. `bt_rem_unity_init` — fresh bt_work has bt_rem_u/v = 1.0
!!      everywhere (so the substep multiplication is a no-op without
!!      the knob fired).
!!   2. `compute_bt_rem_formula` — given (r, HBBL, dt, h) the per-face
!!      bt_rem matches the closed form `Htot/(Htot+r·HBBL·dt)`.
!!   3. `compute_bt_rem_zero_drag` — r=0 ⇒ bt_rem = 1 everywhere
!!      (no-op).
!!   4. `compute_bt_rem_walls_unity` — wall faces (i=1, i=nx+1) keep
!!      1.0 (no damping where velocity is zeroed).
!!   5. `substep_applies_bt_rem` — drives the REAL
!!      `barotropic_substep_nonlinear_interior` (not just the coefficient
!!      formula): hand-set `bt_rem_u/v = 0.5`, one substep from a nonzero
!!      rest-column IC with zero forcing must give EXACTLY u = 0.5*u0.
!!      Closes the "both ends pass, the wire is never exercised" gap the
!!      first four cases share (PR-29 plan §2.4) — for the ORIGINAL
!!      `bt_substep_drag` knob too, not only the newer wave-drag one.
module test_ocean_bt_substep_drag
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_barotropic_workstate, only: barotropic_workstate_t
   use rdb_barotropic_coupling, only: compute_bt_rem
   use rdb_ocean_metrics, only: ocean_metrics_t
   use ocean_test_metrics, only: make_cartesian_metrics, destroy_cartesian_metrics
   use rdb_coriolis_adv, only: coriolis_adv_t
   use rdb_ocean_dyn, only: ocean_dyn_t
   use rdb_barotropic_substep, only: barotropic_substep_nonlinear_interior
   implicit none
   private

   public :: collect_bt_substep_drag_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 2

contains

   subroutine collect_bt_substep_drag_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("bt_rem_unity_init", test_bt_rem_unity_init), &
                  new_unittest("compute_bt_rem_formula", test_bt_rem_formula), &
                  new_unittest("compute_bt_rem_zero_drag", test_bt_rem_zero_drag), &
                  new_unittest("compute_bt_rem_walls_unity", test_bt_rem_walls), &
                  new_unittest("substep_applies_bt_rem", test_substep_applies_bt_rem) &
                  ]
   end subroutine collect_bt_substep_drag_tests

   subroutine build_state(grid, ms, bt_work, nz_ml)
      type(hgrid_t), intent(out) :: grid
      type(multilayer_state_t), intent(out) :: ms
      type(barotropic_workstate_t), intent(out) :: bt_work
      integer, intent(in) :: nz_ml
      call grid%init(6, 6, NGHOST, 1.0_wp, 1.0_wp)
      ms%nz_ml = nz_ml
      call ms%init(grid)
      call bt_work%init(grid, nz_ml=nz_ml)
   end subroutine build_state

   subroutine test_bt_rem_unity_init(error)
      !! Fresh bt_work allocation should leave bt_rem_u/v at exactly
      !! 1.0 (source=1.0_wp in init) ⇒ substep multiplication is a
      !! no-op without compute_bt_rem being called.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(barotropic_workstate_t) :: bt_work
      real(wp) :: dev_u, dev_v
      checks: block
         call build_state(grid, ms, bt_work, NZ)
         dev_u = maxval(abs(bt_work%bt_rem_u - 1.0_wp))
         dev_v = maxval(abs(bt_work%bt_rem_v - 1.0_wp))
         call check(error, dev_u < 1.0e-14_wp, "bt_rem_u init != 1")
         if (allocated(error)) exit checks
         call check(error, dev_v < 1.0e-14_wp, "bt_rem_v init != 1")
      end block checks
      call bt_work%destroy(); call ms%destroy()
   end subroutine test_bt_rem_unity_init

   subroutine test_bt_rem_formula(error)
      !! With (r, HBBL, dt_inner, h_layer) ⇒ bt_rem_face must equal
      !! `Htot/(Htot + r·HBBL·dt_inner)` to roundoff.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(barotropic_workstate_t) :: bt_work
      real(wp), parameter :: R_LINEAR = 2.5e-5_wp
      real(wp), parameter :: HBBL = 10.0_wp
      real(wp), parameter :: DT_INNER = 240.0_wp
      real(wp), parameter :: H_BED = 71.0_wp, H_SURF = 990.0_wp
      real(wp) :: expected, htot
      integer :: i, j
      checks: block
         call build_state(grid, ms, bt_work, NZ)
         ms%h_layer(:, :, 1) = H_BED
         ms%h_layer(:, :, 2) = H_SURF

         call compute_bt_rem(grid, bt_work, ms, R_LINEAR, HBBL, DT_INNER)

         htot = H_BED + H_SURF
         expected = htot/(htot + R_LINEAR*HBBL*DT_INNER)

         i = NGHOST + 3; j = NGHOST + 3
         call check(error, abs(bt_work%bt_rem_u(i, j) - expected) < 1.0e-12_wp, &
                    "bt_rem_u does not match Htot/(Htot+r·HBBL·dt)")
         if (allocated(error)) exit checks
         call check(error, abs(bt_work%bt_rem_v(i, j) - expected) < 1.0e-12_wp, &
                    "bt_rem_v does not match formula")

      end block checks
      call bt_work%destroy(); call ms%destroy()
   end subroutine test_bt_rem_formula

   subroutine test_bt_rem_zero_drag(error)
      !! r = 0 ⇒ bt_rem identically 1 everywhere (no-op).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(barotropic_workstate_t) :: bt_work
      real(wp) :: dev
      checks: block
         call build_state(grid, ms, bt_work, NZ)
         ms%h_layer = 500.0_wp
         call compute_bt_rem(grid, bt_work, ms, 0.0_wp, 10.0_wp, 240.0_wp)
         dev = maxval(abs(bt_work%bt_rem_u - 1.0_wp))
         call check(error, dev < 1.0e-14_wp, &
                    "r=0: bt_rem_u not exactly 1 everywhere")
      end block checks
      call bt_work%destroy(); call ms%destroy()
   end subroutine test_bt_rem_zero_drag

   subroutine test_bt_rem_walls(error)
      !! Wall faces (i=1, i=nx+1, j=1, j=ny+1) must keep bt_rem = 1
      !! (we don't damp the wall, where ubt is zeroed anyway).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(barotropic_workstate_t) :: bt_work
      integer :: nx, ny
      checks: block
         call build_state(grid, ms, bt_work, NZ)
         ms%h_layer = 500.0_wp
         call compute_bt_rem(grid, bt_work, ms, 2.5e-5_wp, 10.0_wp, 240.0_wp)
         nx = size(ms%h_layer, 1)
         ny = size(ms%h_layer, 2)
         call check(error, abs(bt_work%bt_rem_u(1, ny/2) - 1.0_wp) < 1.0e-14_wp, &
                    "west-wall bt_rem_u != 1")
         if (allocated(error)) exit checks
         call check(error, abs(bt_work%bt_rem_u(nx + 1, ny/2) - 1.0_wp) < 1.0e-14_wp, &
                    "east-wall bt_rem_u != 1")
         if (allocated(error)) exit checks
         call check(error, abs(bt_work%bt_rem_v(nx/2, 1) - 1.0_wp) < 1.0e-14_wp, &
                    "south-wall bt_rem_v != 1")
         if (allocated(error)) exit checks
         call check(error, abs(bt_work%bt_rem_v(nx/2, ny + 1) - 1.0_wp) < 1.0e-14_wp, &
                    "north-wall bt_rem_v != 1")

      end block checks
      call bt_work%destroy(); call ms%destroy()
   end subroutine test_bt_rem_walls

   subroutine test_substep_applies_bt_rem(error)
      !! Drives the REAL substep (not just compute_bt_rem's formula):
      !! hand-set bt_rem_u/v = 0.5 uniformly, one substep from a nonzero
      !! rest-column IC (eta=0, u=u0/=0) with zero forcing must give
      !! EXACTLY u = 0.5*u0.  With a spatially uniform IC the interior
      !! probe sees zero eta-gradient/vorticity/KE-gradient tendency at
      !! this first step (those all depend on fields that start uniform),
      !! so the update reduces exactly to bt_rem*(u+dt*0) -- round-off
      !! exact, no wall-reflection contamination.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(ocean_metrics_t) :: metrics
      type(ocean_dyn_t) :: dyn
      type(coriolis_adv_t) :: cor
      real(wp), parameter :: U0 = 0.1_wp, H_REF = 200.0_wp, DT_INNER = 1.0_wp
      integer :: ip, jp
      checks: block
         call grid%init(12, 8, NGHOST, 1000.0_wp, 1000.0_wp)
         call make_cartesian_metrics(metrics, grid)
         ! nz_ml=1 required: barotropic_substep_nonlinear_interior forwards
         ! bt_work%F_bt_u_fast/F_bt_v_fast as the forcing arrays, and those
         ! are only allocated when `init` is passed `nz_ml` (unlike
         ! run_fast_nonlinear's sibling helper in test_ocean_barotropic_substep,
         ! which passes explicit local fu/fv arrays instead).
         call dyn%init(grid, nz_ml=1)
         call cor%init(grid)
         cor%f_corner = 0.0_wp
         dyn%bt_work%bt_H_ref = H_REF
         dyn%bt_work%bt_ubt = U0
         dyn%bt_work%bt_rem_u = 0.5_wp
         dyn%bt_work%bt_rem_v = 0.5_wp

         !$acc enter data copyin(dyn, cor)
         call dyn%enter_data()
         call cor%enter_data()
         call barotropic_substep_nonlinear_interior(grid, metrics, dyn%bt_work, &
                                                    cor%f_corner, 1, DT_INNER)
         !$acc update self(dyn%bt_work%bt_ubt)
         call dyn%exit_data()
         call cor%exit_data()
         !$acc exit data delete(dyn, cor)

         ip = grid%nghost + 6; jp = grid%nghost + 4
         call check(error, dyn%bt_work%bt_ubt(ip, jp) == 0.5_wp*U0, &
                    "substep did not apply the hand-set bt_rem_u exactly")
      end block checks
      call dyn%destroy(); call cor%destroy()
      call destroy_cartesian_metrics(metrics)
   end subroutine test_substep_applies_bt_rem

end module test_ocean_bt_substep_drag
