!! `k_bot` — the shared index of the first LIVE layer counting UP from the
!! bed, and the bed-side consumers routed through it.
!!
!! Under `vcoord_type = "z_fixed"` every column shallower than the
!! nominal stack carries inert FILLERS (`h <= H_VANISHED`) below its
!! partial bed cell, so `k = 1` is not the bed-adjacent layer.  Before
!! `k_bot`, the bottom drag (explicit, HBBL and the implicit-fold rate),
!! the vdiff bed row and the tidal-mixing bed anchor all sat on `k = 1`.
!! The vdiff solve cuts every interface that touches a filler, so an
!! implicit drag folded into row 1 never reached water: the measured
!! day-10 kinetic energy of the 1-degree Southern Ocean was bit-identical
!! with and without it.
!!
!! Cases:
!!   1. the builder on a hand-built staircase: centre rule, the `max` face
!!      rule, agreement with the closed-face mask, the `1` fallback;
!!   2. explicit linear + quadratic bottom drag act on the face's first
!!      live layer and decelerate it at the analytic rate on a uniform-
!!      flow staircase channel, and touch nothing else;
!!   3. the implicit fold (`implicit_drag`) decelerates the first live
!!      layer at the backward-Euler rate `1/(1 + dt·r)` per step and
!!      REMOVES kinetic energy (zero before `k_bot`);
!!   4. a filler column and its sigma twin give the same live-layer answer
!!      through the implicit fold with interior viscosity and through the
!!      HBBL-distributed walk;
!!   5. geothermal heat lands in the first live layer and the budget closes;
!!   6. tidal mixing with the `e_compute` energy source: a filler column
!!      mixes exactly like its sigma twin (before, `N_bot` was read across
!!      two fillers and the source was ZERO);
!!   7. sigma: with `k_bot ≡ 1` the drag tendency is the
!!      literal `k = 1` formula (the bit-identity gate is the goldens).
module test_ocean_zfixed_k_bot
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, H_VANISHED
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_vcoord, only: ocean_vcoord_z_fixed_target_uniform, &
                               ocean_vcoord_k_bot_from_target, &
                               ocean_vcoord_k_top_from_target, &
                               ocean_vcoord_closed_face_masks
   use rdb_ocean_bottom_drag, only: ocean_bottom_drag_t, &
                                    ocean_bottom_drag_compute_tendencies, &
                                    ocean_bottom_drag_apply_tendencies, &
                                    BDRAG_LINEAR, BDRAG_QUADRATIC
   use rdb_ocean_vdiff, only: ocean_vdiff_t, vdiff_apply_momentum
   use rdb_ocean_geothermal, only: ocean_geothermal_t, ocean_geothermal_apply_tracers
   use rdb_ocean_tidal_mixing, only: ocean_tidal_mixing_t, tidal_mixing_compute
   use rdb_eos, only: EOS_VARIANT_LINEAR
   implicit none
   private

   public :: collect_ocean_zfixed_k_bot_tests

   integer, parameter :: NGHOST = 2
   integer, parameter :: NXP = 4, NYP = 3
      !! Physical extent; total = +2*NGHOST.
   integer, parameter :: NZ = 8
   real(wp), parameter :: H_REF = 400.0_wp
   real(wp), parameter :: H_NOM = H_REF/real(NZ, wp)
   real(wp), parameter :: H_FILL = 1.0e-4_wp
      !! `zstar_h_min`; strictly below `H_VANISHED`.
   integer, parameter :: NLIVE = 5, NFILL = 3
      !! Uniform-column twin: NFILL fillers under NLIVE live layers vs a
      !! NLIVE-layer sigma column.
   real(wp), parameter :: H_LIVE = 40.0_wp
   real(wp), parameter :: U0 = 0.3_wp
   real(wp), parameter :: RHO0 = 1035.0_wp

contains

   subroutine collect_ocean_zfixed_k_bot_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("k_bot_staircase_builder", test_builder), &
                  new_unittest("k_bot_no_filler_falls_back_to_1", test_fallback), &
                  new_unittest("explicit_drag_on_first_live_layer", test_explicit_drag), &
                  new_unittest("implicit_fold_decelerates_first_live_layer", test_implicit_fold), &
                  new_unittest("implicit_fold_with_viscosity_matches_sigma", test_fold_matches_sigma), &
                  new_unittest("hbbl_walk_starts_at_k_bot", test_hbbl_matches_sigma), &
                  new_unittest("geothermal_lands_in_first_live_layer", test_geothermal), &
                  new_unittest("tidal_mixing_anchors_at_k_bot", test_tidal), &
                  new_unittest("sigma_drag_is_the_literal_k1_formula", test_sigma_bitid) &
                  ]
   end subroutine collect_ocean_zfixed_k_bot_tests

   ! ------------------------------------------------------------------
   ! Builders
   ! ------------------------------------------------------------------

   subroutine staircase_target(nx, ny, tgt)
      !! Column depth decreasing eastward by 40 m per column from the full
      !! 400 m stack, so successive columns gain bed fillers.
      integer, intent(in) :: nx, ny
      real(wp), intent(out) :: tgt(nx, ny, NZ)
      real(wp) :: total_h(nx, ny), eta(nx, ny), z_top(nx, ny)
      integer :: i

      eta = 0.0_wp
      z_top = 0.0_wp
      do i = 1, nx
         total_h(i, :) = H_REF - 40.0_wp*real(i - 1, wp)
      end do
      call ocean_vcoord_z_fixed_target_uniform(tgt, total_h, eta, z_top, &
                                               nx, ny, NZ, H_NOM, H_FILL)
   end subroutine staircase_target

   subroutine build_staircase(grid, ms)
      !! A staircase channel: h = the z_fixed target, `k_bot` from the
      !! production builder, and a uniform eastward flow `U0` on every face
      !! layer that is live on BOTH sides (closed faces hold 0, as
      !! `mask_layer_velocities` leaves them).
      type(hgrid_t), intent(out) :: grid
      type(multilayer_state_t), intent(out) :: ms
      real(wp), allocatable :: tgt(:, :, :), open_u(:, :, :), open_v(:, :, :)
      integer :: nx, ny

      call grid%init(NXP, NYP, NGHOST, 1000.0_wp, 1000.0_wp)
      ms%nz_ml = NZ
      call ms%init(grid)
      nx = grid%nx_total
      ny = grid%ny_total
      allocate (tgt(nx, ny, NZ), open_u(nx + 1, ny, NZ), open_v(nx, ny + 1, NZ))
      call staircase_target(nx, ny, tgt)
      ms%h_layer = tgt
      call ocean_vcoord_k_bot_from_target(ms%k_bot, ms%k_bot_u, ms%k_bot_v, &
                                          tgt, nx, ny, NZ, H_VANISHED)
      call ocean_vcoord_closed_face_masks(open_u, open_v, tgt, nx, ny, NZ, H_VANISHED)
      ms%u_face_x_layer = U0*open_u
      ms%v_face_y_layer = 0.0_wp
   end subroutine build_staircase

   subroutine build_twin(grid, ms, with_fill)
      !! One repeated column: NLIVE live 40 m layers, under which (when
      !! `with_fill`) NFILL inert fillers sit.  T/S follow I1′: a filler
      !! carries its donor's (the first live layer's) concentration.
      type(hgrid_t), intent(out) :: grid
      type(multilayer_state_t), intent(out) :: ms
      logical, intent(in) :: with_fill
      integer :: k, kb, nzl, idx_t, idx_s
      real(wp) :: t_k

      call grid%init(NXP, NYP, NGHOST, 1000.0_wp, 1000.0_wp)
      kb = 1
      if (with_fill) kb = NFILL + 1
      nzl = NLIVE + kb - 1
      ms%nz_ml = nzl
      call ms%init(grid)
      idx_t = ms%idx_temperature
      idx_s = ms%idx_salinity
      do k = 1, nzl
         t_k = 2.0_wp + 0.5_wp*real(max(k, kb) - kb, wp)   ! stable, warmer up
         if (k < kb) then
            ms%h_layer(:, :, k) = H_FILL
         else
            ms%h_layer(:, :, k) = H_LIVE
         end if
         ms%tracers(idx_t)%hTr(:, :, k) = ms%h_layer(:, :, k)*t_k
         ms%tracers(idx_s)%hTr(:, :, k) = ms%h_layer(:, :, k)*34.5_wp
         if (k < kb) then
            ms%u_face_x_layer(:, :, k) = 0.0_wp
         else
            ms%u_face_x_layer(:, :, k) = U0*(1.0_wp + 0.1_wp*real(k - kb, wp))
         end if
      end do
      ms%v_face_y_layer = 0.0_wp
      ms%k_bot = kb
      ms%k_bot_u = kb
      ms%k_bot_v = kb
   end subroutine build_twin

   subroutine map_in(ms)
      type(multilayer_state_t), intent(inout) :: ms
      !$omp target enter data map(to: ms)
      call ms%enter_data()
   end subroutine map_in

   subroutine map_out(ms)
      type(multilayer_state_t), intent(inout) :: ms
      !$omp target update from(ms%u_face_x_layer, ms%v_face_y_layer)
      call ms%exit_data()
      !$omp target exit data map(delete: ms)
   end subroutine map_out

   ! ------------------------------------------------------------------
   ! 1. Builder
   ! ------------------------------------------------------------------

   subroutine test_builder(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 10, NY = 3
      real(wp) :: tgt(NX, NY, NZ), open_u(NX + 1, NY, NZ), open_v(NX, NY + 1, NZ)
      integer :: k_bot(NX, NY), k_bot_u(NX + 1, NY), k_bot_v(NX, NY + 1)
      integer :: k_top(NX, NY), k_top_u(NX + 1, NY), k_top_v(NX, NY + 1)
      integer :: i, j, k
      logical :: ok

      call staircase_target(NX, NY, tgt)
      call ocean_vcoord_k_bot_from_target(k_bot, k_bot_u, k_bot_v, tgt, NX, NY, NZ, H_VANISHED)
      call ocean_vcoord_k_top_from_target(k_top, k_top_u, k_top_v, tgt, NX, NY, NZ, H_VANISHED)
      call ocean_vcoord_closed_face_masks(open_u, open_v, tgt, NX, NY, NZ, H_VANISHED)

      call check(error, k_bot(1, 2) == 1 .and. k_bot(NX, 2) > 2, &
                 "the staircase really moves k_bot (1 at full depth, >2 at the shallow end)")
      if (allocated(error)) return

      ! Centre rule: first live layer from the bed, every layer below a filler.
      ok = .true.
      do j = 1, NY
         do i = 1, NX
            if (tgt(i, j, k_bot(i, j)) <= H_VANISHED) ok = .false.
            do k = 1, k_bot(i, j) - 1
               if (tgt(i, j, k) > H_VANISHED) ok = .false.
            end do
            if (k_bot(i, j) > k_top(i, j)) ok = .false.
         end do
      end do
      call check(error, ok, "k_bot is the first live layer and everything below it is a filler")
      if (allocated(error)) return

      ! Face rule: MAX of the two columns (the shallower bottom).
      ok = .true.
      do j = 1, NY
         do i = 2, NX
            if (k_bot_u(i, j) /= max(k_bot(i - 1, j), k_bot(i, j))) ok = .false.
         end do
      end do
      do j = 2, NY
         do i = 1, NX
            if (k_bot_v(i, j) /= max(k_bot(i, j - 1), k_bot(i, j))) ok = .false.
         end do
      end do
      call check(error, ok, "face twins are the max of their two columns")
      if (allocated(error)) return

      ! Agreement with the closed-face mask: closed below k_bot_u, open at it.
      ok = .true.
      do j = 1, NY
         do i = 2, NX
            do k = 1, NZ
               if (k < k_bot_u(i, j) .and. open_u(i, j, k) /= 0.0_wp) ok = .false.
            end do
            if (open_u(i, j, k_bot_u(i, j)) /= 1.0_wp) ok = .false.
         end do
      end do
      call check(error, ok, "open_u is closed below k_bot_u and open AT it")
   end subroutine test_builder

   subroutine test_fallback(error)
      !! No bed filler anywhere ⇒ `≡ 1`; a dead (all-marker) column ⇒ `1`.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 5, NY = 4
      real(wp) :: tgt(NX, NY, NZ), total_h(NX, NY), eta(NX, NY), z_top(NX, NY)
      integer :: k_bot(NX, NY), k_bot_u(NX + 1, NY), k_bot_v(NX, NY + 1)

      total_h = H_REF
      eta = 0.0_wp
      z_top = 0.0_wp
      call ocean_vcoord_z_fixed_target_uniform(tgt, total_h, eta, z_top, &
                                               NX, NY, NZ, H_NOM, H_FILL)
      call ocean_vcoord_k_bot_from_target(k_bot, k_bot_u, k_bot_v, tgt, NX, NY, NZ, H_VANISHED)
      call check(error, all(k_bot == 1) .and. all(k_bot_u == 1) .and. all(k_bot_v == 1), &
                 "full-depth column: k_bot ≡ 1 on centres and both faces")
      if (allocated(error)) return

      tgt = H_VANISHED
      call ocean_vcoord_k_bot_from_target(k_bot, k_bot_u, k_bot_v, tgt, NX, NY, NZ, H_VANISHED)
      call check(error, all(k_bot == 1) .and. all(k_bot_u == 1) .and. all(k_bot_v == 1), &
                 "a dead column falls back to 1")
   end subroutine test_fallback

   ! ------------------------------------------------------------------
   ! 2. Explicit drag on the staircase channel
   ! ------------------------------------------------------------------

   subroutine test_explicit_drag(error)
      !! Linear: `du(k_bot_u) = -r·U0` exactly on every interior face,
      !! nothing elsewhere, and after N forward-Euler steps the first live
      !! layer reads `U0·(1 - dt·r)^N`.  Quadratic: `du = -c_d·|U|·U/h_face`
      !! with `h_face` the two-column mean AT `k_bot_u`.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_bottom_drag_t) :: bd
      real(wp), parameter :: R = 1.0e-4_wp, DT = 600.0_wp, CD = 2.5e-3_wp
      integer, parameter :: NSTEP = 10
      real(wp), allocatable :: u_ic(:, :, :)
      integer :: i, j, k, kb, nx, ny
      real(wp) :: err_lin, err_q, err_other, err_decay, err_rest, hf, expect
      logical :: shallow_seen

      call build_staircase(grid, ms)
      nx = grid%nx_total
      ny = grid%ny_total
      call bd%init(grid, nz_ml=NZ)
      bd%variant = BDRAG_LINEAR
      bd%r_linear = R
      allocate (u_ic, source=ms%u_face_x_layer)

      call map_in(ms)
      !$omp target enter data map(to: bd)
      call bd%enter_data()
      call ocean_bottom_drag_compute_tendencies(grid, bd, ms, DT)
      !$omp target update from(bd%du_drag%data)

      err_lin = 0.0_wp
      err_other = 0.0_wp
      shallow_seen = .false.
      do j = 1, ny
         do i = 2, nx
            kb = ms%k_bot_u(i, j)
            if (kb > 1) shallow_seen = .true.
            err_lin = max(err_lin, abs(bd%du_drag%data(i, j, kb) + R*U0))
            do k = 1, NZ
               if (k /= kb) err_other = max(err_other, abs(bd%du_drag%data(i, j, k)))
            end do
         end do
      end do

      do i = 1, NSTEP
         call ocean_bottom_drag_compute_tendencies(grid, bd, ms, DT)
         call ocean_bottom_drag_apply_tendencies(bd, ms, DT)
      end do
      !$omp target update from(ms%u_face_x_layer)
      err_decay = 0.0_wp
      err_rest = 0.0_wp
      expect = U0*(1.0_wp - DT*R)**NSTEP
      do j = 1, ny
         do i = 2, nx
            kb = ms%k_bot_u(i, j)
            err_decay = max(err_decay, abs(ms%u_face_x_layer(i, j, kb) - expect))
            do k = 1, NZ
               if (k /= kb) err_rest = max(err_rest, &
                                           abs(ms%u_face_x_layer(i, j, k) - u_ic(i, j, k)))
            end do
         end do
      end do

      ! Quadratic, one evaluation from the post-decay state.
      bd%variant = BDRAG_QUADRATIC
      bd%c_drag = CD
      call ocean_bottom_drag_compute_tendencies(grid, bd, ms, DT)
      !$omp target update from(bd%du_drag%data)
      err_q = 0.0_wp
      do j = 1, ny
         do i = 2, nx
            kb = ms%k_bot_u(i, j)
            hf = max(0.5_wp*(ms%h_layer(i - 1, j, kb) + ms%h_layer(i, j, kb)), bd%h_min)
            expect = -CD*abs(ms%u_face_x_layer(i, j, kb))*ms%u_face_x_layer(i, j, kb)/hf
            err_q = max(err_q, abs(bd%du_drag%data(i, j, kb) - expect)/abs(expect))
         end do
      end do

      call bd%exit_data()
      !$omp target exit data map(delete: bd)
      call map_out(ms)

      call check(error, shallow_seen, "the channel has faces whose first live layer is not k=1")
      if (allocated(error)) return
      call check(error, err_lin == 0.0_wp, "linear drag on k_bot_u is exactly -r*U0")
      if (allocated(error)) return
      call check(error, err_other == 0.0_wp, "no drag tendency on any other layer (fillers included)")
      if (allocated(error)) return
      call check(error, err_decay < 1.0e-14_wp, &
                 "first live layer decays as U0*(1-dt*r)^N")
      if (allocated(error)) return
      call check(error, err_rest == 0.0_wp, "every other layer is untouched")
      if (allocated(error)) return
      call check(error, err_q < 1.0e-13_wp, &
                 "quadratic drag uses |U|/h at the face's first live layer")
      call bd%destroy()
      call ms%destroy()
   end subroutine test_explicit_drag

   ! ------------------------------------------------------------------
   ! 3. Implicit fold on the staircase channel
   ! ------------------------------------------------------------------

   subroutine test_implicit_fold(error)
      !! `&ocean_vdiff_nml implicit_drag` with the closed-face coupling cut
      !! (`zlevel_faces`, production z_fixed) and K_v = 0: the fold must
      !! decelerate the first live layer as `U0/(1 + dt·r)^N` and remove
      !! kinetic energy.  Before `k_bot` the rate sat on row 1 — a filler
      !! the cut decoupled — and the energy loss was exactly ZERO.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_bottom_drag_t) :: bd
      type(ocean_vdiff_t) :: vd
      real(wp), parameter :: R = 1.0e-4_wp, DT = 600.0_wp
      integer, parameter :: NSTEP = 10
      real(wp), allocatable :: tau_u(:, :), tau_v(:, :)
      integer :: i, j, k, kb, nx, ny, step
      real(wp) :: ke0, ke1, err_decay, err_rest, expect, hf
      real(wp), allocatable :: u_ic(:, :, :)

      call build_staircase(grid, ms)
      nx = grid%nx_total
      ny = grid%ny_total
      call bd%init(grid, nz_ml=NZ)
      bd%variant = BDRAG_LINEAR
      bd%r_linear = R
      bd%implicit_fold = .true.
      call vd%init(grid, nz_ml=NZ)
      vd%K_v_momentum = 0.0_wp
      vd%implicit_drag = .true.
      vd%zlevel_faces = .true.
      allocate (tau_u(nx + 1, ny), source=0.0_wp)
      allocate (tau_v(nx, ny + 1), source=0.0_wp)
      allocate (u_ic, source=ms%u_face_x_layer)

      ! Kinetic energy over the faces that HAVE bed fillers (`k_bot_u > 1`)
      ! only — full-depth faces were always dragged on k = 1.
      ke0 = 0.0_wp
      do k = 1, NZ
         do j = 1, ny
            do i = 2, nx
               if (ms%k_bot_u(i, j) == 1) cycle
               hf = min(ms%h_layer(i - 1, j, k), ms%h_layer(i, j, k))
               ke0 = ke0 + 0.5_wp*hf*ms%u_face_x_layer(i, j, k)**2
            end do
         end do
      end do

      call map_in(ms)
      !$omp target enter data map(to: bd, vd)
      call bd%enter_data()
      call vd%enter_data()
      !$omp target enter data map(to: tau_u, tau_v)
      do step = 1, NSTEP
         call ocean_bottom_drag_compute_tendencies(grid, bd, ms, DT)
         call vdiff_apply_momentum(grid, vd, ms, DT, tau_u=tau_u, tau_v=tau_v, &
                                   lambda_bot_u=bd%lambda_bot_u, &
                                   lambda_bot_v=bd%lambda_bot_v, rho0=RHO0)
      end do
      !$omp target exit data map(delete: tau_u, tau_v)
      call vd%exit_data()
      call bd%exit_data()
      !$omp target exit data map(delete: bd, vd)
      call map_out(ms)

      ke1 = 0.0_wp
      err_decay = 0.0_wp
      err_rest = 0.0_wp
      expect = U0/(1.0_wp + DT*R)**NSTEP
      do k = 1, NZ
         do j = 1, ny
            do i = 2, nx
               kb = ms%k_bot_u(i, j)
               hf = min(ms%h_layer(i - 1, j, k), ms%h_layer(i, j, k))
               if (kb > 1) ke1 = ke1 + 0.5_wp*hf*ms%u_face_x_layer(i, j, k)**2
               if (k == kb) then
                  err_decay = max(err_decay, abs(ms%u_face_x_layer(i, j, k) - expect))
               else
                  err_rest = max(err_rest, abs(ms%u_face_x_layer(i, j, k) - u_ic(i, j, k)))
               end if
            end do
         end do
      end do

      call check(error, ke1 < ke0 - 1.0e-6_wp*ke0, &
                 "implicit drag removes kinetic energy on the filler faces (exactly zero on k=1)")
      if (allocated(error)) return
      call check(error, err_decay < 1.0e-14_wp, &
                 "first live layer decays at the backward-Euler rate U0/(1+dt*r)^N")
      if (allocated(error)) return
      call check(error, err_rest < 1.0e-15_wp, &
                 "K_v = 0: no other layer moves")
      call vd%destroy()
      call bd%destroy()
      call ms%destroy()
   end subroutine test_implicit_fold

   ! ------------------------------------------------------------------
   ! 4. Filler column == its sigma twin
   ! ------------------------------------------------------------------

   subroutine test_fold_matches_sigma(error)
      !! Quadratic implicit fold WITH interior viscosity: the live layers of
      !! the filler column must reproduce the sigma column (same water, no
      !! fillers) — the fillers are identity rows, the bed row is `k_bot`.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: u_s(NLIVE), u_z(NLIVE)
      integer :: m
      real(wp) :: err

      call run_fold(.false., u_s)
      call run_fold(.true., u_z)
      err = 0.0_wp
      do m = 1, NLIVE
         err = max(err, abs(u_z(m) - u_s(m)))
      end do
      call check(error, abs(u_s(1) - U0) > 1.0e-4_wp, "the drag actually acted on the sigma bed")
      if (allocated(error)) return
      call check(error, err <= 1.0e-15_wp, &
                 "implicit fold + viscosity: filler column == sigma column on the live layers")
   end subroutine test_fold_matches_sigma

   subroutine run_fold(with_fill, u_live)
      logical, intent(in) :: with_fill
      real(wp), intent(out) :: u_live(NLIVE)
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_bottom_drag_t) :: bd
      type(ocean_vdiff_t) :: vd
      real(wp), parameter :: DT = 1800.0_wp, CD = 2.5e-3_wp
      real(wp), allocatable :: tau_u(:, :), tau_v(:, :)
      integer :: nx, ny, kb, m, ig, jg

      call build_twin(grid, ms, with_fill)
      nx = grid%nx_total
      ny = grid%ny_total
      kb = ms%k_bot(1, 1)
      call bd%init(grid, nz_ml=ms%nz_ml)
      bd%variant = BDRAG_QUADRATIC
      bd%c_drag = CD
      bd%implicit_fold = .true.
      call vd%init(grid, nz_ml=ms%nz_ml)
      vd%K_v_momentum = 1.0e-2_wp
      vd%implicit_drag = .true.
      vd%zlevel_faces = with_fill
      allocate (tau_u(nx + 1, ny), source=0.0_wp)
      allocate (tau_v(nx, ny + 1), source=0.0_wp)

      call map_in(ms)
      !$omp target enter data map(to: bd, vd)
      call bd%enter_data()
      call vd%enter_data()
      !$omp target enter data map(to: tau_u, tau_v)
      call ocean_bottom_drag_compute_tendencies(grid, bd, ms, DT)
      call vdiff_apply_momentum(grid, vd, ms, DT, tau_u=tau_u, tau_v=tau_v, &
                                lambda_bot_u=bd%lambda_bot_u, &
                                lambda_bot_v=bd%lambda_bot_v, rho0=RHO0)
      !$omp target exit data map(delete: tau_u, tau_v)
      call vd%exit_data()
      call bd%exit_data()
      !$omp target exit data map(delete: bd, vd)
      call map_out(ms)

      ig = NGHOST + 2
      jg = NGHOST + 2
      do m = 1, NLIVE
         u_live(m) = ms%u_face_x_layer(ig, jg, kb + m - 1)
      end do
      call vd%destroy()
      call bd%destroy()
      call ms%destroy()
   end subroutine run_fold

   subroutine test_hbbl_matches_sigma(error)
      !! HBBL-distributed quadratic drag (`hbbl = 60 m`, two live layers
      !! deep, `bed_factor = 2`): the band walk starts at `k_bot`, so the
      !! fillers take no share and the live tendencies equal the sigma
      !! column's.  Before, the walk started at k = 1 and the fillers
      !! received the bed factor and a full-rate share of the band.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: du_s(NLIVE + NFILL), du_z(NLIVE + NFILL)
      integer :: m
      real(wp) :: err, fill_max

      call run_hbbl(.false., du_s)
      call run_hbbl(.true., du_z)
      err = 0.0_wp
      do m = 1, NLIVE
         err = max(err, abs(du_z(NFILL + m) - du_s(m)))
      end do
      fill_max = maxval(abs(du_z(1:NFILL)))
      call check(error, du_s(1) < 0.0_wp .and. du_s(2) < 0.0_wp .and. du_s(3) == 0.0_wp, &
                 "the band covers exactly the two bed-most live layers")
      if (allocated(error)) return
      call check(error, fill_max == 0.0_wp, "no HBBL share on a bed filler")
      if (allocated(error)) return
      call check(error, err == 0.0_wp, "HBBL tendencies on the live layers == sigma column")
   end subroutine test_hbbl_matches_sigma

   subroutine run_hbbl(with_fill, du)
      logical, intent(in) :: with_fill
      real(wp), intent(out) :: du(NLIVE + NFILL)
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_bottom_drag_t) :: bd
      integer :: k, ig, jg

      call build_twin(grid, ms, with_fill)
      call bd%init(grid, nz_ml=ms%nz_ml)
      bd%variant = BDRAG_QUADRATIC
      bd%c_drag = 2.5e-3_wp
      bd%hbbl = 60.0_wp
      bd%bed_factor = 2.0_wp
      call map_in(ms)
      !$omp target enter data map(to: bd)
      call bd%enter_data()
      call ocean_bottom_drag_compute_tendencies(grid, bd, ms, 600.0_wp)
      !$omp target update from(bd%du_drag%data)
      call bd%exit_data()
      !$omp target exit data map(delete: bd)
      call map_out(ms)
      ig = NGHOST + 2
      jg = NGHOST + 2
      du = 0.0_wp
      do k = 1, ms%nz_ml
         du(k) = bd%du_drag%data(ig, jg, k)
      end do
      call bd%destroy()
      call ms%destroy()
   end subroutine run_hbbl

   ! ------------------------------------------------------------------
   ! 5. Geothermal
   ! ------------------------------------------------------------------

   subroutine test_geothermal(error)
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_geothermal_t) :: geo
      real(wp), parameter :: DT = 3600.0_wp, QGEO = 0.1_wp
      real(wp), allocatable :: ht0(:, :, :)
      real(wp) :: src, err_kb, err_other, budget_err
      integer :: i, j, k, kb, nx, ny

      call build_twin(grid, ms, .true.)
      nx = grid%nx_total
      ny = grid%ny_total
      kb = NFILL + 1
      call geo%init(grid)
      geo%enable = .true.
      geo%q_geo_const = QGEO
      geo%rho0 = RHO0
      allocate (ht0, source=ms%tracers(ms%idx_temperature)%hTr)
      src = DT*QGEO/(geo%rho0*geo%cp)

      !$omp target enter data map(to: ms, geo)
      call ms%enter_data()
      call ocean_geothermal_apply_tracers(grid, geo, ms, DT)
      associate (hT => ms%tracers(ms%idx_temperature)%hTr, &
                 bgeo => ms%heat_budget_geothermal)
         !$omp target update from(hT, bgeo)
      end associate
      call ms%exit_data()
      !$omp target exit data map(delete: ms, geo)

      err_kb = 0.0_wp
      err_other = 0.0_wp
      do k = 1, ms%nz_ml
         do j = 1, ny
            do i = 1, nx
               if (k == kb) then
                  err_kb = max(err_kb, abs(ms%tracers(ms%idx_temperature)%hTr(i, j, k) &
                                           - ht0(i, j, k) - src))
               else
                  err_other = max(err_other, abs(ms%tracers(ms%idx_temperature)%hTr(i, j, k) &
                                                 - ht0(i, j, k)))
               end if
            end do
         end do
      end do
      budget_err = abs(sum(ms%heat_budget_geothermal) - src*real(nx*ny, wp))/(src*real(nx*ny, wp))

      call check(error, err_kb < 1.0e-14_wp, "geothermal heat lands in the first live layer")
      if (allocated(error)) return
      call check(error, err_other == 0.0_wp, "no heat in a filler or any other layer")
      if (allocated(error)) return
      call check(error, budget_err < 1.0e-13_wp, "geothermal budget closes")
      call geo%destroy()
      call ms%destroy()
   end subroutine test_geothermal

   ! ------------------------------------------------------------------
   ! 6. Tidal mixing
   ! ------------------------------------------------------------------

   subroutine test_tidal(error)
      !! `e_compute` (Jayne & St Laurent): `E ∝ N_bot`.  On the filler
      !! column `N_bot` must be read across the two deepest LIVE layers —
      !! the sigma twin's `n2_col(1)` — so every live interface diffusivity
      !! matches the twin.  Before, it read across two fillers carrying the
      !! same donor T/S: `N_bot = 0`, no mixing at all.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: kd_s(NLIVE + NFILL + 1), kd_z(NLIVE + NFILL + 1)
      integer :: m
      real(wp) :: err

      call run_tidal(.false., kd_s)
      call run_tidal(.true., kd_z)
      err = 0.0_wp
      do m = 1, NLIVE + 1
         err = max(err, abs(kd_z(NFILL + m) - kd_s(m))/max(abs(kd_s(m)), 1.0e-30_wp))
      end do
      call check(error, maxval(kd_s) > 0.0_wp, "the sigma column mixes (E > 0)")
      if (allocated(error)) return
      call check(error, maxval(abs(kd_z(1:NFILL))) == 0.0_wp, &
                 "no tidal Kd at or below the filler interfaces")
      if (allocated(error)) return
      call check(error, err < 1.0e-12_wp, "live interfaces match the sigma column")
   end subroutine test_tidal

   subroutine run_tidal(with_fill, kd)
      logical, intent(in) :: with_fill
      real(wp), intent(out) :: kd(NLIVE + NFILL + 1)
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_tidal_mixing_t) :: tm
      integer :: k

      call build_twin(grid, ms, with_fill)
      call tm%init(grid, nz_ml=ms%nz_ml)
      tm%enable = .true.
      tm%gamma = 0.3333_wp
      tm%mu = 0.2_wp
      tm%zeta = 100.0_wp
      tm%kd_max = -1.0_wp
      tm%min_zbot = 0.0_wp
      tm%eos%variant = EOS_VARIANT_LINEAR
      tm%eos%alpha_T = 0.2_wp
      tm%eos%beta_S = RHO0*7.6e-4_wp
      tm%eos%rho0 = RHO0
      tm%rho0 = RHO0
      call tm%set_e_uniform(0.0_wp)
      tm%e_compute = .true.
      tm%kappa_itides = 6.2832e-4_wp
      tm%kappa_h2 = 1.0_wp
      tm%utide = 0.05_wp
      tm%h2_rough = 200.0_wp
      tm%frac_rough = 0.1_wp
      tm%e_max = 1.0e3_wp

      !$omp target enter data map(to: ms, tm)
      call ms%enter_data()
      call tm%enter_data()
      call tidal_mixing_compute(grid, tm, ms, 1800.0_wp)
      !$omp target update from(tm%kd_int)
      call tm%exit_data()
      call ms%exit_data()
      !$omp target exit data map(delete: ms, tm)

      kd = 0.0_wp
      do k = 1, ms%nz_ml + 1
         kd(k) = tm%kd_int(NGHOST + 2, NGHOST + 2, k)
      end do
      call tm%destroy()
      call ms%destroy()
   end subroutine run_tidal

   ! ------------------------------------------------------------------
   ! 7. Sigma bit-identity
   ! ------------------------------------------------------------------

   subroutine test_sigma_bitid(error)
      !! `k_bot` at its `1` default (every coordinate but z_fixed): the
      !! quadratic tendency IS the literal `k = 1` expression (to the last
      !! bit of a host/device FMA difference) and every other layer is
      !! exactly zero.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      type(ocean_bottom_drag_t) :: bd
      real(wp), parameter :: CD = 2.5e-3_wp, DT = 600.0_wp
      integer :: i, j, k, nx, ny
      real(wp) :: v_at_u, sp, hf, expect
      logical :: exact, others_zero

      call build_twin(grid, ms, .false.)
      nx = grid%nx_total
      ny = grid%ny_total
      do j = 1, ny + 1
         do i = 1, nx
            ms%v_face_y_layer(i, j, :) = 0.01_wp*real(i - j, wp)
         end do
      end do
      call bd%init(grid, nz_ml=ms%nz_ml)
      bd%variant = BDRAG_QUADRATIC
      bd%c_drag = CD
      call map_in(ms)
      !$omp target enter data map(to: bd)
      call bd%enter_data()
      call ocean_bottom_drag_compute_tendencies(grid, bd, ms, DT)
      !$omp target update from(bd%du_drag%data)
      call bd%exit_data()
      !$omp target exit data map(delete: bd)
      call map_out(ms)

      call check(error, all(ms%k_bot == 1) .and. all(ms%k_bot_u == 1), &
                 "k_bot defaults to 1 off z_fixed")
      if (allocated(error)) return
      exact = .true.
      others_zero = .true.
      do j = 1, ny
         do i = 2, nx
            v_at_u = 0.25_wp*(ms%v_face_y_layer(i - 1, j, 1) + ms%v_face_y_layer(i, j, 1) + &
                              ms%v_face_y_layer(i - 1, j + 1, 1) + ms%v_face_y_layer(i, j + 1, 1))
            hf = max(0.5_wp*(ms%h_layer(i - 1, j, 1) + ms%h_layer(i, j, 1)), bd%h_min)
            sp = sqrt(ms%u_face_x_layer(i, j, 1)**2 + v_at_u**2)
            expect = min(ms%wet_mask(i - 1, j), ms%wet_mask(i, j))* &
                     (-CD*sp*ms%u_face_x_layer(i, j, 1)/(hf + 0.0_wp*CD*sp))
            ! Host vs device FMA contraction may differ in the last bit;
            ! the tree-wide bit-identity gate is the golden suite.
            if (abs(bd%du_drag%data(i, j, 1) - expect) > 2.0e-15_wp*abs(expect)) exact = .false.
            do k = 2, ms%nz_ml
               if (bd%du_drag%data(i, j, k) /= 0.0_wp) others_zero = .false.
            end do
         end do
      end do
      call check(error, exact, "sigma: k=1 quadratic tendency is the literal formula (to FMA rounding)")
      if (allocated(error)) return
      call check(error, others_zero, "sigma: no tendency above k=1")
      call bd%destroy()
      call ms%destroy()
   end subroutine test_sigma_bitid

end module test_ocean_zfixed_k_bot
