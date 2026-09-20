!! Analytical tests for the ice-shelf TOP drag (`rdb_ocean_top_drag`,
!! `&ocean_tdrag_nml`) — the momentum sink at `k = nz` under an ice
!! shelf, the mirror of the bottom drag at `k = 1`.
!!
!! Nothing here asserts bit-equality between two DIFFERENT code paths;
!! every comparison carries a stated tolerance derived from the operation
!! count, because an FMA-fusing build is free to associate the same
!! algebra differently.  Bit-identity is asserted only where one path is
!! required to leave an array UNTOUCHED.
!!
!! Cases:
!!   * `quadratic_layer_only_analytic` — one step of the layer-`nz`
!!     quadratic drag reproduces `u - dt*C_d*|u|*u/h_nz` exactly;
!!     layers `k < nz` and uncovered faces are untouched.
!!   * `mirror_of_bottom_drag` — reflect the column (h and u flipped in
!!     `k`) and the distributed top drag returns the distributed bottom
!!     drag's answer, layer for layer.  This is the test that keeps the
!!     two kernels one closure rather than two.
!!   * `linear_distributed_decay_explicit` / `_implicit` — closed-form
!!     spin-down at the DISTRIBUTED rate `r*htbl/h_nz`, over N steps,
!!     both time-discretisations.
!!   * `implicit_stable_where_explicit_overshoots` — at `dt*lambda = 2.5`
!!     the explicit form amplifies and flips sign on step one; the
!!     implicit form stays same-signed and strictly decreasing.
!!   * `drag_removes_kinetic_energy` — a general multilayer flow loses
!!     KE, every face contracts toward zero, and no face overshoots.
!!   * `calving_front_face_mask` — the OR rule puts drag ON the frontal
!!     face and none one face further out.
!!   * `disabled_is_exact_no_op` — knob off ⇒ byte-identical velocities
!!     and an essentially-zero byte count.
!!   * `parse_and_gate` — `form` round-trips; a typo is `TDRAG_INVALID`.
!!   * `validate_refusals` — the configure-time matrix, including the
!!     ONE-`C_d` agreement rule against `&ocean_cavity_melt_nml
!!     cdrag_top`.
module test_ocean_top_drag
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_config, only: config_t, read_config_from_string, validate_config
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_CONFIG_VALIDATE
   use rdb_ocean_bottom_drag, only: ocean_bottom_drag_t, BDRAG_QUADRATIC, &
                                    ocean_bottom_drag_compute_tendencies
   use rdb_ocean_top_drag, only: ocean_top_drag_t, &
                                 TDRAG_LINEAR, TDRAG_QUADRATIC, TDRAG_INVALID, &
                                 ocean_top_drag_compute_tendencies, &
                                 ocean_top_drag_apply_tendencies, &
                                 top_drag_fill_face_cover_impl, &
                                 parse_tdrag_variant, tdrag_variant_is_implemented
   implicit none
   private

   public :: collect_ocean_top_drag_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 3
   real(wp), parameter :: TOL_REL = 1.0e-12_wp
      !! Relative tolerance for a closed-form comparison.  The kernels
      !! here do O(10) flops per face per step and the decay tests take
      !! at most 24 steps, so the accumulated relative round-off is
      !! O(100*eps) ~ 2e-14 in double; 1e-12 is two decades of headroom
      !! and still four decades below any physically interesting error.

contains

   subroutine collect_ocean_top_drag_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("quadratic_layer_only_analytic", test_quadratic_layer_only), &
                  new_unittest("mirror_of_bottom_drag", test_mirror_of_bottom), &
                  new_unittest("linear_distributed_decay_explicit", test_linear_decay_explicit), &
                  new_unittest("linear_distributed_decay_implicit", test_linear_decay_implicit), &
                  new_unittest("implicit_stable_where_explicit_overshoots", test_implicit_stable), &
                  new_unittest("drag_removes_kinetic_energy", test_energy_sink), &
                  new_unittest("calving_front_face_mask", test_calving_front), &
                  new_unittest("disabled_is_exact_no_op", test_disabled_no_op), &
                  new_unittest("parse_and_gate", test_parse_and_gate), &
                  new_unittest("validate_refusals", test_validate_refusals) &
                  ]
   end subroutine collect_ocean_top_drag_tests

   ! ------------------------------------------------------------------
   ! Harness
   ! ------------------------------------------------------------------

   subroutine make_grid(grid, nx_phys, ny_phys)
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nx_phys, ny_phys
      call grid%init(nx_phys, ny_phys, NGHOST, 1.0_wp, 1.0_wp)
   end subroutine make_grid

   subroutine cover_all(td)
      !! Fill the slot's face + centre cover masks with "everything is
      !! under ice".  Host-side, before `map_in`.
      type(ocean_top_drag_t), intent(inout) :: td
      td%cover_u(:, :) = 1.0_wp
      td%cover_v(:, :) = 1.0_wp
      td%cover_t(:, :) = 1.0_wp
   end subroutine cover_all

   subroutine map_in(ms, td)
      !! `mem:separate` discipline: map the state object AND the slot,
      !! and PUSH the host-set cover masks + velocities to the device.
      !! `ocean_top_drag_t%enter_data` uses `copyin` for its 2-D fields,
      !! so the masks arrive with the map; the `!$acc update device` is
      !! belt-and-braces for anything a case sets after the map.  Every
      !! directive here is an inert no-op on a host build.
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_top_drag_t), intent(inout) :: td
      !$acc enter data copyin(ms)
      call ms%enter_data()
      !$acc enter data copyin(td)
      call td%enter_data()
      !$acc update device(ms%u_face_x_layer, ms%v_face_y_layer, &
      !$acc               ms%h_layer, ms%wet_mask)
      !$acc update device(td%cover_u, td%cover_v, td%cover_t)
   end subroutine map_in

   subroutine map_out(ms, td)
      !! Pull the COMPONENT arrays the checks read — never the aggregate
      !! derived type (an aggregate D->H copy overwrites the host
      !! descriptors with device addresses).
      type(multilayer_state_t), intent(inout) :: ms
      type(ocean_top_drag_t), intent(inout) :: td
      !$acc update self(ms%u_face_x_layer, ms%v_face_y_layer)
      !$acc update self(td%du_drag%data, td%dv_drag%data, td%stress_top)
      call td%exit_data()
      !$acc exit data delete(td)
      call ms%exit_data()
      !$acc exit data delete(ms)
   end subroutine map_out

   ! ------------------------------------------------------------------
   ! (i) explicit quadratic, layer-nz only
   ! ------------------------------------------------------------------

   subroutine test_quadratic_layer_only(error)
      !! Uniform `u = U0` in the TOP layer under full cover, `v = 0`, so
      !! the face speed is exactly `|U0|` and the one-step answer is the
      !! closed form `U0 - dt*C_d*U0^2/h_nz`.  Layers `k < nz` must be
      !! untouched EXACTLY (they are never written), and so must the
      !! faces of the uncovered half of the domain.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_top_drag_t) :: td
      real(wp), parameter :: U0 = 0.5_wp
      real(wp), parameter :: CD = 2.5e-3_wp
      real(wp), parameter :: H = 20.0_wp
      real(wp), parameter :: DT = 100.0_wp
      integer :: i, j, nx, ny, i_split
      real(wp) :: expect, got, max_lower, max_open
      checks: block

         call make_grid(grid, 12, 10)
         ms%nz_ml = NZ
         call ms%init(grid)
         td%enable = .true.
         call td%init(grid, nz_ml=NZ)
         td%variant = TDRAG_QUADRATIC
         td%c_drag = CD
         td%rho0 = 1027.51_wp
         nx = grid%nx_total
         ny = grid%ny_total
         i_split = nx/2

         ms%h_layer = H
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         ms%u_face_x_layer(:, :, NZ) = U0

         ! Covered west of `i_split`, open east of it.  Cell-centred
         ! cover; the face masks are the OR of the two abutting cells,
         ! so `cover_u(i_split+1)` is covered (the frontal face) and
         ! `cover_u(i_split+2)` is not.
         td%cover_t(:, :) = 0.0_wp
         td%cover_t(1:i_split, :) = 1.0_wp
         call top_drag_fill_face_cover_impl(td%cover_u, td%cover_v, td%cover_t, nx, ny)

         call map_in(ms, td)
         call ocean_top_drag_compute_tendencies(td, ms, DT)
         call ocean_top_drag_apply_tendencies(td, ms, DT)
         call map_out(ms, td)

         expect = U0 - DT*CD*U0*U0/H
         got = ms%u_face_x_layer(i_split/2, ny/2, NZ)
         call check(error, abs(got - expect) <= TOL_REL*abs(expect), &
                    "layer-only quadratic top drag off its closed form")
         if (allocated(error)) exit checks

         ! Layers below the top are never written — exact.
         max_lower = max(maxval(abs(ms%u_face_x_layer(:, :, 1:NZ - 1))), &
                         maxval(abs(ms%v_face_y_layer(:, :, 1:NZ - 1))))
         call check(error, max_lower == 0.0_wp, &
                    "top drag leaked into layers k < nz")
         if (allocated(error)) exit checks

         ! Open faces (both cells uncovered) keep U0 exactly.
         max_open = 0.0_wp
         do j = 1, ny
            do i = i_split + 2, nx
               max_open = max(max_open, abs(ms%u_face_x_layer(i, j, NZ) - U0))
            end do
         end do
         call check(error, max_open == 0.0_wp, &
                    "top drag fired on an UNCOVERED face")
         if (allocated(error)) exit checks

         ! The stress diagnostic is the same algebra, cell-centred.
         call check(error, abs(td%stress_top(i_split/2, ny/2) - &
                               td%rho0*CD*U0*U0) <= TOL_REL*td%rho0*CD*U0*U0, &
                    "stress_top off rho0*C_d*|U|^2 under cover")
         if (allocated(error)) exit checks
         call check(error, td%stress_top(nx - 1, ny/2) == 0.0_wp, &
                    "stress_top non-zero on an uncovered cell")

      end block checks
      call td%destroy(); call ms%destroy()
   end subroutine test_quadratic_layer_only

   ! ------------------------------------------------------------------
   ! (i, part 2) mirror symmetry against the bottom drag
   ! ------------------------------------------------------------------

   subroutine test_mirror_of_bottom(error)
      !! Reflect the column about its middle — `h` and `u` flipped in
      !! `k` — and the HBBL-distributed bottom drag and the
      !! HTBL-distributed top drag must return the same per-layer
      !! tendency, layer for mirrored layer.  A non-uniform `h` is
      !! essential: it is what makes the band walk (partial last layer,
      !! `h_in/h_face` weighting) do real work rather than reduce to a
      !! single layer.
      !!
      !! Tolerance, not bit-equality: the two kernels are separate
      !! translation units and a fusing build may associate the band
      !! integral differently.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms_b, ms_t
      type(ocean_bottom_drag_t) :: bd
      type(ocean_top_drag_t) :: td
      real(wp), parameter :: CD = 2.5e-3_wp
      real(wp), parameter :: BL = 6.0_wp
      real(wp), parameter :: DT = 60.0_wp
      real(wp), parameter :: H_BOT(NZ) = [5.0_wp, 3.0_wp, 8.0_wp]
      real(wp), parameter :: U_BOT(NZ) = [0.4_wp, 0.2_wp, 0.1_wp]
      integer :: k, nx, ny, i_p, j_p
      real(wp) :: a_bot, a_top, scale
      checks: block

         call make_grid(grid, 12, 10)
         nx = grid%nx_total
         ny = grid%ny_total
         i_p = nx/2
         j_p = ny/2

         ms_b%nz_ml = NZ
         call ms_b%init(grid)
         ms_t%nz_ml = NZ
         call ms_t%init(grid)
         call bd%init(grid, nz_ml=NZ)
         bd%variant = BDRAG_QUADRATIC
         bd%c_drag = CD
         bd%hbbl = BL
         td%enable = .true.
         call td%init(grid, nz_ml=NZ)
         td%variant = TDRAG_QUADRATIC
         td%c_drag = CD
         td%htbl = BL
         call cover_all(td)

         ms_b%u_face_x_layer = 0.0_wp
         ms_b%v_face_y_layer = 0.0_wp
         ms_t%u_face_x_layer = 0.0_wp
         ms_t%v_face_y_layer = 0.0_wp
         do k = 1, NZ
            ms_b%h_layer(:, :, k) = H_BOT(k)
            ms_b%u_face_x_layer(:, :, k) = U_BOT(k)
            ! The reflection: layer k of the top case IS layer nz+1-k of
            ! the bottom case.
            ms_t%h_layer(:, :, NZ + 1 - k) = H_BOT(k)
            ms_t%u_face_x_layer(:, :, NZ + 1 - k) = U_BOT(k)
         end do

         !$acc enter data copyin(ms_b)
         call ms_b%enter_data()
         !$acc enter data copyin(bd)
         call bd%enter_data()
         !$acc update device(ms_b%u_face_x_layer, ms_b%v_face_y_layer, &
         !$acc               ms_b%h_layer, ms_b%wet_mask)
         call ocean_bottom_drag_compute_tendencies(grid, bd, ms_b, DT)
         !$acc update self(bd%du_drag%data)

         call map_in(ms_t, td)
         call ocean_top_drag_compute_tendencies(td, ms_t, DT)
         call map_out(ms_t, td)

         do k = 1, NZ
            a_bot = bd%du_drag%data(i_p, j_p, k)
            a_top = td%du_drag%data(i_p, j_p, NZ + 1 - k)
            scale = max(abs(a_bot), 1.0e-30_wp)
            call check(error, abs(a_top - a_bot) <= TOL_REL*scale, &
                       "top drag is not the mirror of the bottom drag at the "// &
                       "reflected layer")
            if (allocated(error)) exit checks
         end do
         ! ... and the band really did span more than one layer, else
         ! the mirror above would be a much weaker statement.
         call check(error, abs(bd%du_drag%data(i_p, j_p, 2)) > 0.0_wp, &
                    "the HBBL/HTBL band collapsed to a single layer — the "// &
                    "mirror test would then not exercise the band walk")

      end block checks
      call bd%exit_data()
      !$acc exit data delete(bd)
      call ms_b%exit_data()
      !$acc exit data delete(ms_b)
      call bd%destroy(); call td%destroy()
      call ms_b%destroy(); call ms_t%destroy()
   end subroutine test_mirror_of_bottom

   ! ------------------------------------------------------------------
   ! (ii) linear spin-down under a full lid — the analytic decay rate
   ! ------------------------------------------------------------------

   subroutine spin_down_case(implicit_form, u_final, rate, error)
      !! Barotropic-uniform flow under a full lid, flat bed, no rotation,
      !! LINEAR top drag distributed over `htbl`.
      !!
      !! `htbl = 4 m` inside uniform `h = 10 m` layers: the band is the
      !! top layer alone (the walk exits once the cumulative thickness
      !! reaches `htbl`) and its share of the stress is
      !! `h_in/h_face = 4/10`.  So the top layer relaxes at the
      !! DISTRIBUTED rate
      !!
      !!     lambda = r * htbl / h_nz = 0.4 * r
      !!
      !! and the exact N-step solutions are
      !!
      !!     explicit:  u_N = U0 * (1 - dt*lambda)^N
      !!     implicit:  u_N = U0 / (1 + dt*lambda)^N
      !!
      !! This routine runs the kernel and hands both back; the two tests
      !! below assert against the closed forms.
      logical, intent(in) :: implicit_form
      real(wp), intent(out) :: u_final, rate
      type(error_type), allocatable, intent(inout) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_top_drag_t) :: td
      real(wp), parameter :: U0 = 0.3_wp
      real(wp), parameter :: R = 1.0e-3_wp
      real(wp), parameter :: H = 10.0_wp
      real(wp), parameter :: HTBL = 4.0_wp
      real(wp), parameter :: DT = 50.0_wp
      integer, parameter :: N_STEPS = 24
      integer :: step, nx, ny

      u_final = 0.0_wp
      rate = R*HTBL/H

      call make_grid(grid, 12, 10)
      ms%nz_ml = NZ
      call ms%init(grid)
      td%enable = .true.
      call td%init(grid, nz_ml=NZ)
      td%variant = TDRAG_LINEAR
      td%r_linear = R
      td%htbl = HTBL
      td%implicit = implicit_form
      call cover_all(td)
      nx = grid%nx_total
      ny = grid%ny_total

      ms%h_layer = H
      ms%u_face_x_layer = U0
      ms%v_face_y_layer = 0.0_wp

      call map_in(ms, td)
      do step = 1, N_STEPS
         call ocean_top_drag_compute_tendencies(td, ms, DT)
         call ocean_top_drag_apply_tendencies(td, ms, DT)
      end do
      call map_out(ms, td)

      u_final = ms%u_face_x_layer(nx/2, ny/2, NZ)
      call check(error, abs(ms%u_face_x_layer(nx/2, ny/2, 1) - U0) == 0.0_wp, &
                 "the bed layer moved under a TOP drag")
      call td%destroy(); call ms%destroy()
   end subroutine spin_down_case

   subroutine test_linear_decay_explicit(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: U0 = 0.3_wp, DT = 50.0_wp
      integer, parameter :: N_STEPS = 24
      real(wp) :: u_final, rate, expect

      call spin_down_case(.false., u_final, rate, error)
      if (allocated(error)) return
      expect = U0*(1.0_wp - DT*rate)**N_STEPS
      call check(error, abs(u_final - expect) <= TOL_REL*abs(expect), &
                 "explicit distributed linear spin-down off (1 - dt*r*htbl/h)^N")
   end subroutine test_linear_decay_explicit

   subroutine test_linear_decay_implicit(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: U0 = 0.3_wp, DT = 50.0_wp
      integer, parameter :: N_STEPS = 24
      real(wp) :: u_final, rate, expect

      call spin_down_case(.true., u_final, rate, error)
      if (allocated(error)) return
      expect = U0/(1.0_wp + DT*rate)**N_STEPS
      call check(error, abs(u_final - expect) <= TOL_REL*abs(expect), &
                 "implicit distributed linear spin-down off U0/(1 + dt*r*htbl/h)^N")
   end subroutine test_linear_decay_implicit

   subroutine test_implicit_stable(error)
      !! At `dt*lambda = 2.5` the explicit form has amplification
      !! `1 - 2.5 = -1.5`: it flips the sign of the velocity and GROWS
      !! it, which is the failure the implicit form exists to remove.
      !! The implicit amplification is `1/3.5`, so the flow stays
      !! same-signed and strictly decreasing.  Both statements are
      !! asserted, so the test cannot pass by the two forms agreeing.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_top_drag_t) :: td
      real(wp), parameter :: U0 = 0.2_wp
      real(wp), parameter :: H = 10.0_wp
      real(wp), parameter :: DT = 100.0_wp
      real(wp), parameter :: R = 2.5e-2_wp
         !! `dt*r = 2.5` with `htbl = 0` (layer-only ⇒ lambda = r).
      integer, parameter :: N_STEPS = 6
      integer :: step, nx, ny, variant_pass
      real(wp) :: u_prev, u_now, u_expl_1
      logical :: monotone, same_sign
      checks: block

         ! --- explicit: one step is enough to show the overshoot ---
         call make_grid(grid, 8, 6)
         ms%nz_ml = NZ
         call ms%init(grid)
         td%enable = .true.
         call td%init(grid, nz_ml=NZ)
         td%variant = TDRAG_LINEAR
         td%r_linear = R
         td%implicit = .false.
         call cover_all(td)
         nx = grid%nx_total
         ny = grid%ny_total
         ms%h_layer = H
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         ms%u_face_x_layer(:, :, NZ) = U0
         call map_in(ms, td)
         call ocean_top_drag_compute_tendencies(td, ms, DT)
         call ocean_top_drag_apply_tendencies(td, ms, DT)
         call map_out(ms, td)
         u_expl_1 = ms%u_face_x_layer(nx/2, ny/2, NZ)
         call td%destroy(); call ms%destroy()

         call check(error, u_expl_1*U0 < 0.0_wp .and. abs(u_expl_1) > abs(U0), &
                    "the explicit control must actually overshoot at dt*lambda = 2.5 "// &
                    "(else the implicit test below proves nothing)")
         if (allocated(error)) exit checks

         ! --- implicit: monotone, same-signed, contracting ---
         call make_grid(grid, 8, 6)
         ms%nz_ml = NZ
         call ms%init(grid)
         td%enable = .true.
         call td%init(grid, nz_ml=NZ)
         td%variant = TDRAG_LINEAR
         td%r_linear = R
         td%implicit = .true.
         call cover_all(td)
         nx = grid%nx_total
         ny = grid%ny_total
         ms%h_layer = H
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         ms%u_face_x_layer(:, :, NZ) = U0
         call map_in(ms, td)
         monotone = .true.
         same_sign = .true.
         u_prev = U0
         do step = 1, N_STEPS
            call ocean_top_drag_compute_tendencies(td, ms, DT)
            call ocean_top_drag_apply_tendencies(td, ms, DT)
            !$acc update self(ms%u_face_x_layer)
            u_now = ms%u_face_x_layer(nx/2, ny/2, NZ)
            if (.not. (abs(u_now) < abs(u_prev))) monotone = .false.
            if (u_now*U0 <= 0.0_wp) same_sign = .false.
            u_prev = u_now
         end do
         call map_out(ms, td)
         variant_pass = td%variant

         call check(error, monotone, &
                    "implicit top drag is not strictly contracting at dt*lambda = 2.5")
         if (allocated(error)) exit checks
         call check(error, same_sign, &
                    "implicit top drag reversed the flow — it must not overshoot")
         if (allocated(error)) exit checks
         call check(error, abs(u_prev - U0/(1.0_wp + DT*R)**N_STEPS) <= &
                    TOL_REL*abs(U0/(1.0_wp + DT*R)**N_STEPS), &
                    "implicit top drag off U0/(1 + dt*r)^N")
         if (allocated(error)) exit checks
         call check(error, variant_pass == TDRAG_LINEAR, "variant tag was mutated")

      end block checks
      call td%destroy(); call ms%destroy()
   end subroutine test_implicit_stable

   ! ------------------------------------------------------------------
   ! (iii) energy: the top drag only ever REMOVES kinetic energy
   ! ------------------------------------------------------------------

   subroutine test_energy_sink(error)
      !! A general (non-uniform, sheared, two-component) flow under full
      !! cover.  After one step:
      !!   * every face has contracted toward zero — `|u_new| <= |u_old|`
      !!     and no sign flip (a drag is a contraction, not a rotation);
      !!   * the column-integrated kinetic energy `sum h*(u^2+v^2)/2`
      !!     has STRICTLY decreased.
      !! The per-face statement is the stronger of the two: an energy sum
      !! can fall while individual faces are accelerated.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_top_drag_t) :: td
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp), parameter :: CD = 2.5e-3_wp
      real(wp), parameter :: DT = 200.0_wp
      real(wp), allocatable :: u0(:, :, :), v0(:, :, :)
      real(wp) :: ke_before, ke_after
      integer :: i, j, k, nx, ny
      logical :: contracts
      checks: block

         call make_grid(grid, 14, 12)
         ms%nz_ml = NZ
         call ms%init(grid)
         td%enable = .true.
         call td%init(grid, nz_ml=NZ)
         td%variant = TDRAG_QUADRATIC
         td%c_drag = CD
         td%htbl = 12.0_wp
         td%drag_bg_vel = 0.05_wp
         td%rho0 = 1027.51_wp
         call cover_all(td)
         nx = grid%nx_total
         ny = grid%ny_total

         do k = 1, NZ
            ms%h_layer(:, :, k) = 6.0_wp + 2.0_wp*real(k, wp)
            do j = 1, ny
               do i = 1, nx + 1
                  ms%u_face_x_layer(i, j, k) = &
                     0.3_wp*sin(2.0_wp*PI*real(i + 2*k, wp)/real(nx, wp))
               end do
            end do
            do j = 1, ny + 1
               do i = 1, nx
                  ms%v_face_y_layer(i, j, k) = &
                     0.2_wp*cos(2.0_wp*PI*real(j + k, wp)/real(ny, wp))
               end do
            end do
         end do
         allocate (u0, source=ms%u_face_x_layer)
         allocate (v0, source=ms%v_face_y_layer)

         ke_before = column_ke(ms%u_face_x_layer, ms%v_face_y_layer, ms%h_layer, nx, ny)

         call map_in(ms, td)
         call ocean_top_drag_compute_tendencies(td, ms, DT)
         call ocean_top_drag_apply_tendencies(td, ms, DT)
         call map_out(ms, td)

         contracts = .true.
         do k = 1, NZ
            do j = 1, ny
               do i = 2, nx
                  if (abs(ms%u_face_x_layer(i, j, k)) > abs(u0(i, j, k))) contracts = .false.
                  if (ms%u_face_x_layer(i, j, k)*u0(i, j, k) < 0.0_wp) contracts = .false.
               end do
            end do
            do j = 2, ny
               do i = 1, nx
                  if (abs(ms%v_face_y_layer(i, j, k)) > abs(v0(i, j, k))) contracts = .false.
                  if (ms%v_face_y_layer(i, j, k)*v0(i, j, k) < 0.0_wp) contracts = .false.
               end do
            end do
         end do
         ke_after = column_ke(ms%u_face_x_layer, ms%v_face_y_layer, ms%h_layer, nx, ny)

         call check(error, contracts, &
                    "top drag accelerated or reversed a face — it is a momentum SINK")
         if (allocated(error)) exit checks
         call check(error, ke_after < ke_before, &
                    "top drag did not reduce the kinetic energy")

      end block checks
      call td%destroy(); call ms%destroy()
   end subroutine test_energy_sink

   pure function column_ke(u, v, h, nx, ny) result(ke)
      !! Layer-thickness-weighted kinetic energy over the interior
      !! faces.  Nearest-cell thickness (not a face average) — the
      !! monotonicity statement is about the velocities, and any fixed
      !! positive weighting makes it.
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: u(:, :, :), v(:, :, :), h(:, :, :)
      real(wp) :: ke
      integer :: i, j, k
      ke = 0.0_wp
      do k = 1, size(h, 3)
         do j = 1, ny
            do i = 2, nx
               ke = ke + 0.5_wp*h(i, j, k)*u(i, j, k)*u(i, j, k)
            end do
         end do
         do j = 2, ny
            do i = 1, nx
               ke = ke + 0.5_wp*h(i, j, k)*v(i, j, k)*v(i, j, k)
            end do
         end do
      end do
   end function column_ke

   ! ------------------------------------------------------------------
   ! (iv) the calving front
   ! ------------------------------------------------------------------

   subroutine test_calving_front(error)
      !! `cover_frac` steps from 1 to 0 between cells `i0` and `i0+1`.
      !! The OR rule puts the frontal face `i0+1` (between the last
      !! covered cell and the first open one) UNDER ice; face `i0+2`,
      !! whose two cells are both open, is not.  This is the documented
      !! choice, and it is asserted on both the mask and the tendency so
      !! a change of rule cannot slip through as a mask-only edit.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_top_drag_t) :: td
      real(wp), parameter :: CD = 2.5e-3_wp, U0 = 0.4_wp, DT = 60.0_wp
      integer :: nx, ny, i0, j_p
      checks: block

         call make_grid(grid, 12, 10)
         ms%nz_ml = NZ
         call ms%init(grid)
         td%enable = .true.
         call td%init(grid, nz_ml=NZ)
         td%variant = TDRAG_QUADRATIC
         td%c_drag = CD
         nx = grid%nx_total
         ny = grid%ny_total
         i0 = nx/2
         j_p = ny/2

         td%cover_t(:, :) = 0.0_wp
         td%cover_t(1:i0, :) = 1.0_wp
         call top_drag_fill_face_cover_impl(td%cover_u, td%cover_v, td%cover_t, nx, ny)

         call check(error, td%cover_u(i0, j_p) == 1.0_wp, &
                    "an interior covered face must be under ice")
         if (allocated(error)) exit checks
         call check(error, td%cover_u(i0 + 1, j_p) == 1.0_wp, &
                    "the CALVING-FRONT face (one covered cell, one open) must be "// &
                    "under ice — the OR rule")
         if (allocated(error)) exit checks
         call check(error, td%cover_u(i0 + 2, j_p) == 0.0_wp, &
                    "a face with two OPEN cells must not be under ice")
         if (allocated(error)) exit checks

         ms%h_layer = 20.0_wp
         ms%u_face_x_layer = 0.0_wp
         ms%v_face_y_layer = 0.0_wp
         ms%u_face_x_layer(:, :, NZ) = U0

         call map_in(ms, td)
         call ocean_top_drag_compute_tendencies(td, ms, DT)
         call map_out(ms, td)

         call check(error, td%du_drag%data(i0 + 1, j_p, NZ) < 0.0_wp, &
                    "no drag at the calving-front face")
         if (allocated(error)) exit checks
         call check(error, td%du_drag%data(i0 + 2, j_p, NZ) == 0.0_wp, &
                    "drag fired one face beyond the front")

      end block checks
      call td%destroy(); call ms%destroy()
   end subroutine test_calving_front

   ! ------------------------------------------------------------------
   ! (vi) knob off
   ! ------------------------------------------------------------------

   subroutine test_disabled_no_op(error)
      !! A disabled slot launches no kernel, allocates only placeholders,
      !! and leaves the velocities BYTE-identical.  `enable = .false.` is
      !! the default, so this is the statement that keeps every shipped
      !! namelist bit-identical.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_top_drag_t) :: td
      real(wp), allocatable :: u0(:, :, :)
      checks: block

         call make_grid(grid, 10, 8)
         ms%nz_ml = NZ
         call ms%init(grid)
         td%enable = .false.
         call td%init(grid, nz_ml=NZ)
         td%variant = TDRAG_QUADRATIC
         td%c_drag = 2.5e-3_wp

         ms%h_layer = 20.0_wp
         ms%u_face_x_layer = 0.7_wp
         ms%v_face_y_layer = -0.3_wp
         allocate (u0, source=ms%u_face_x_layer)

         call map_in(ms, td)
         call ocean_top_drag_compute_tendencies(td, ms, 100.0_wp)
         call ocean_top_drag_apply_tendencies(td, ms, 100.0_wp)
         call map_out(ms, td)

         call check(error, all(ms%u_face_x_layer == u0), &
                    "a disabled top-drag slot changed the velocity")
         if (allocated(error)) exit checks
         call check(error, all(ms%v_face_y_layer == -0.3_wp), &
                    "a disabled top-drag slot changed the meridional velocity")
         if (allocated(error)) exit checks
         call check(error, td%bytes() < int(1000, kind(td%bytes())), &
                    "a gated-off top-drag slot must count essentially zero bytes")

      end block checks
      call td%destroy(); call ms%destroy()
   end subroutine test_disabled_no_op

   ! ------------------------------------------------------------------
   ! (v) parsing + the configure-time matrix
   ! ------------------------------------------------------------------

   subroutine test_parse_and_gate(error)
      type(error_type), allocatable, intent(out) :: error
      checks: block
         call check(error, parse_tdrag_variant("quadratic") == TDRAG_QUADRATIC, &
                    "'quadratic' must parse to TDRAG_QUADRATIC")
         if (allocated(error)) exit checks
         call check(error, parse_tdrag_variant("linear") == TDRAG_LINEAR, &
                    "'linear' must parse to TDRAG_LINEAR")
         if (allocated(error)) exit checks
         call check(error, parse_tdrag_variant("quadratci") == TDRAG_INVALID, &
                    "a typo must be TDRAG_INVALID, never defaulted")
         if (allocated(error)) exit checks
         call check(error, tdrag_variant_is_implemented(TDRAG_QUADRATIC) .and. &
                    tdrag_variant_is_implemented(TDRAG_LINEAR), &
                    "both shipped variants must report as implemented")
         if (allocated(error)) exit checks
         call check(error,.not. tdrag_variant_is_implemented(TDRAG_INVALID), &
                    "TDRAG_INVALID must not report as implemented")
      end block checks
   end subroutine test_parse_and_gate

   subroutine test_validate_refusals(error)
      !! The configure-time matrix.  The load-bearing row is the ONE-`C_d`
      !! agreement rule: this model has a single ice-base drag
      !! coefficient, and a run that sets two different ones is told
      !! rather than silently given one of them.
      type(error_type), allocatable, intent(out) :: error

      call expect_refused("top drag without the cavity geometry group", &
                          tdrag_nml("enable = .true., cd = 2.5e-3", .false., ""), error)
      if (allocated(error)) return
      ! A typo in `form` is caught EARLIER, by the strict schema's
      ! `allowed=` set, so it never reaches `validate_config`.  Asserted
      ! at the layer that actually refuses it — `validate_config`'s own
      ! `parse_tdrag_variant` gate stays as the belt-and-braces for a
      ! programmatically built `config_t`, and `parse_and_gate` above
      ! pins that gate directly.
      call expect_parse_refused("an unrecognised top-drag form", &
                                tdrag_nml("enable = .true., form = 'quadratci', "// &
                                          "cd = 2.5e-3", .true., ""), error)
      if (allocated(error)) return
      call expect_refused("quadratic top drag with cd = 0", &
                          tdrag_nml("enable = .true.", .true., ""), error)
      if (allocated(error)) return
      call expect_refused("linear top drag with r = 0", &
                          tdrag_nml("enable = .true., form = 'linear'", .true., ""), error)
      if (allocated(error)) return
      call expect_refused("a top-drag cd that disagrees with the melt cdrag_top", &
                          tdrag_nml("enable = .true., cd = 1.0e-3", .true., &
                                    melt_group("cdrag_top = 2.5e-3")), error)
      if (allocated(error)) return

      call expect_accepted("top drag alongside the cavity geometry", &
                           tdrag_nml("enable = .true., cd = 2.5e-3", .true., ""), error)
      if (allocated(error)) return
      call expect_accepted("linear top drag with r > 0", &
                           tdrag_nml("enable = .true., form = 'linear', r = 1.0e-4", &
                                     .true., ""), error)
      if (allocated(error)) return
      call expect_accepted("top drag and melt sharing ONE cd", &
                           tdrag_nml("enable = .true., cd = 2.5e-3", .true., &
                                     melt_group("cdrag_top = 2.5e-3")), error)
      if (allocated(error)) return
      call expect_accepted("the default (top drag off)", &
                           tdrag_nml("", .true., ""), error)
   end subroutine test_validate_refusals

   function melt_group(body) result(g)
      !! `&ocean_cavity_melt_nml` in its minimal in-envelope form, plus
      !! whatever the caller wants to override.
      character(len=*), intent(in) :: body
      character(len=:), allocatable :: g
      g = "&ocean_cavity_melt_nml enable = .true., "//body//" /"
   end function melt_group

   function tdrag_nml(tdrag_body, with_geometry, extra) result(nml)
      !! One complete test namelist.  Every group appears EXACTLY ONCE —
      !! the strict schema refuses a duplicate group, and a refusal test
      !! that tripped on that would be testing the parser.
      character(len=*), intent(in) :: tdrag_body
         !! `&ocean_tdrag_nml` body without the group name or `/`; empty
         !! ⇒ the group is omitted entirely.
      logical, intent(in) :: with_geometry
         !! Append `&ocean_cavity_dyn_nml` (the top-drag prerequisite).
      character(len=*), intent(in) :: extra
         !! One extra complete group (with name and `/`), or empty.
      character(len=:), allocatable :: nml

      nml = "&sim_nml sim_type = 'ocean' /"//new_line("a")// &
            "&grid_nml nx = 8, ny = 6, nghost = 2, dx = 1000.0, dy = 1000.0 /"// &
            new_line("a")// &
            "&nonhydrostatic_nml nz_layers = 4 /"//new_line("a")// &
            "&time_nml t_end = 3600.0, dt_fixed = 60.0 /"//new_line("a")// &
            "&ocean_topo_nml max_depth = 1000.0 /"//new_line("a")// &
            "&ocean_pgf_nml form = 'fv_mom6' /"//new_line("a")// &
            "&vcoord_nml vcoord_type = 'sigma' /"//new_line("a")// &
            "&ocean_bt_nml auto_n_inner = .false., n_inner = 8 /"//new_line("a")// &
            "&ocean_diag_nml enabled = .false. /"//new_line("a")// &
            "&output_nml output_to_file = .false. /"//new_line("a")// &
            "&ocean_eos_nml tfreeze_set = 'isomip' /"//new_line("a")// &
            "&ocean_forcing_nml enable_components = .true. /"//new_line("a")
      if (with_geometry) then
         nml = nml//"&ocean_cavity_dyn_nml enable = .true., draft_config = 'flat', "// &
               "draft_depth = 300.0 /"//new_line("a")
      end if
      if (len_trim(tdrag_body) > 0) then
         nml = nml//"&ocean_tdrag_nml "//tdrag_body//" /"//new_line("a")
      end if
      if (len_trim(extra) > 0) nml = nml//extra//new_line("a")
   end function tdrag_nml

   subroutine expect_refused(what, nml, error)
      !! The namelist must PARSE and then be refused by `validate_config`.
      character(len=*), intent(in) :: what, nml
      type(error_type), allocatable, intent(inout) :: error
      type(config_t) :: cfg
      integer :: ierr

      call read_config_from_string(nml, cfg, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK, &
                 "the refusal case must PARSE cleanly (else it tests the schema, "// &
                 "not the rule): "//what)
      if (allocated(error)) return
      call validate_config(cfg, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_ERR_CONFIG_VALIDATE, &
                 "validate_config must refuse: "//what)
   end subroutine expect_refused

   subroutine expect_parse_refused(what, nml, error)
      !! The namelist must be refused by the strict SCHEMA, before
      !! `validate_config` ever sees it.
      character(len=*), intent(in) :: what, nml
      type(error_type), allocatable, intent(inout) :: error
      type(config_t) :: cfg
      integer :: ierr

      call read_config_from_string(nml, cfg, ierr=ierr)
      call check(error, ierr /= OCEAN_STATUS_OK, &
                 "the strict schema must refuse: "//what)
   end subroutine expect_parse_refused

   subroutine expect_accepted(what, nml, error)
      character(len=*), intent(in) :: what, nml
      type(error_type), allocatable, intent(inout) :: error
      type(config_t) :: cfg
      integer :: ierr
      call read_config_from_string(nml, cfg, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "must parse: "//what)
      if (allocated(error)) return
      call validate_config(cfg, ierr=ierr)
      call check(error, ierr == OCEAN_STATUS_OK, "must be accepted: "//what)
   end subroutine expect_accepted

end module test_ocean_top_drag
