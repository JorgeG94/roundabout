!! The top-side consumers routed through `k_top` (P6.4) — melt, the
!! surface-flux deposit and its budget mirror, and the ice-shelf top
!! drag.
!!
!! Every case is the SAME PHYSICAL COLUMN twice: once as a sigma-like
!! stack of `NLIVE` live layers, and once as a `z_fixed` stack with
!! `NFILL` inert fillers sitting on top of those same live layers inside
!! the ice draft.  The filler column's answer must equal the sigma
!! column's — that is the whole claim of the slice, and it is asserted
!! bit-for-bit wherever the arithmetic permits.
!!
!! ### The leak these tests fence, written down as a number
!!
!! Before `k_top`, the deposit went to `k = nz`.  On a covered column
!! under `z_fixed` that layer is `zstar_h_min = 1.0e-4 m` of water, so
!! `Q_heat = 100 W/m^2` over `dt = 900 s` deposits
!! `dt*Q/(rho_0*cp) = 2.18e-2 K m` into it — an implied concentration of
!! `2.18e-2 / 1.0e-4 = 218 degC`.  The vdiff tracer matrix decouples a
!! vanished row to the identity, so it is not diffused down; the next ALE
!! remap drains it on `h_old <= H_FLOOR`, so it is not kept; and
!! `heat_budget_surface` counted it, so the console reports a source that
!! never entered the ocean.  A per-thermo-step leak of exactly the
!! deposited flux, on every covered column.
module test_ocean_ktop_consumers
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rdb_constants, only: wp, H_VANISHED
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_surface_flux, only: ocean_surface_flux_t, &
                                     ocean_surface_flux_apply_tracers, &
                                     SEAWATER_CP
   use rdb_ocean_cavity_flux, only: cavity_far_field_impl
   use rdb_ocean_top_drag, only: top_drag_tendencies_impl, TDRAG_LINEAR
   implicit none
   private

   public :: collect_ocean_ktop_consumers_tests

   integer, parameter :: NGHOST = 1
   integer, parameter :: NLIVE = 8
      !! Live layers, identical in both stacks.
   integer, parameter :: NFILL = 3
      !! Inert filler layers the ice draft leaves on top of them.
   integer, parameter :: NZ_S = NLIVE
      !! Sigma-like reference stack.
   integer, parameter :: NZ_Z = NLIVE + NFILL
      !! z_fixed stack with the fillers.
   real(wp), parameter :: H_LIVE = 50.0_wp
   real(wp), parameter :: H_FILL = 1.0e-4_wp
      !! `zstar_h_min`; strictly below `H_VANISHED = 1.5e-4`.
   real(wp), parameter :: T0 = -1.2_wp, S0 = 34.5_wp
   real(wp), parameter :: DT_STEP = 900.0_wp
   real(wp), parameter :: RHO0 = 1028.0_wp

contains

   subroutine collect_ocean_ktop_consumers_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("far_field_sample_ignores_the_fillers", test_far_field), &
                  new_unittest("melt_deposit_lands_on_the_live_layer", test_deposit_live), &
                  new_unittest("deposit_matches_the_sigma_column", test_deposit_matches), &
                  new_unittest("surface_budget_closes_from_step_zero", test_budget_closes), &
                  new_unittest("top_drag_spin_down_matches_no_filler", test_top_drag_spindown), &
                  new_unittest("top_drag_fold_rate_is_on_the_live_row", test_fold_row) &
                  ]
   end subroutine collect_ocean_ktop_consumers_tests

   ! ------------------------------------------------------------------
   ! Column builders
   ! ------------------------------------------------------------------

   subroutine build_state(ms, grid, nz, nfill)
      !! A 3x3(+ghost) plane of one repeated column: `nfill` fillers on
      !! top of `NLIVE` 50 m layers, uniform `T0`/`S0` in the live part
      !! and EXACTLY ZERO tracer load in the fillers (which is what the
      !! remap drain leaves there, and what makes a filler's implied
      !! concentration meaningless).  `k_top` is set the way
      !! `configure_ocean_k_top` sets it.
      type(multilayer_state_t), intent(out) :: ms
      type(hgrid_t), intent(out) :: grid
      integer, intent(in) :: nz, nfill
      integer :: k, idx_t, idx_s

      call grid%init(3, 3, NGHOST, 1000.0_wp, 1000.0_wp)
      ms%nz_ml = nz
      call ms%init(grid)
      idx_t = ms%idx_temperature
      idx_s = ms%idx_salinity
      do k = 1, nz
         if (k > nz - nfill) then
            ms%h_layer(:, :, k) = H_FILL
            ms%tracers(idx_t)%hTr(:, :, k) = 0.0_wp
            ms%tracers(idx_s)%hTr(:, :, k) = 0.0_wp
         else
            ms%h_layer(:, :, k) = H_LIVE
            ms%tracers(idx_t)%hTr(:, :, k) = H_LIVE*T0
            ms%tracers(idx_s)%hTr(:, :, k) = H_LIVE*S0
         end if
      end do
      ms%k_top = nz - nfill
      ms%k_top_u = nz - nfill
      ms%k_top_v = nz - nfill
   end subroutine build_state

   subroutine apply_flux(grid, ms, q_heat, q_salt)
      !! One thermo-step surface deposit through the production path.
      type(hgrid_t), intent(in) :: grid
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: q_heat, q_salt
      type(ocean_surface_flux_t) :: sf

      call sf%init(grid)
      sf%rho0 = RHO0
      sf%cp = SEAWATER_CP
      call sf%set_surface_flux_const(q_heat, q_salt)

      !$acc enter data copyin(ms, sf)
      call ms%enter_data()
      call sf%enter_data()
      call ocean_surface_flux_apply_tracers(grid, sf, ms, DT_STEP)
      associate (hT => ms%tracers(ms%idx_temperature)%hTr, &
                 hS => ms%tracers(ms%idx_salinity)%hTr, &
                 bT => ms%heat_budget_surface, bS => ms%salt_budget_surface)
         !$acc update self(hT, hS, bT, bS)
      end associate
      call sf%exit_data()
      call ms%exit_data()
      !$acc exit data delete(ms, sf)
      call sf%destroy()
   end subroutine apply_flux

   ! ------------------------------------------------------------------
   ! Cases
   ! ------------------------------------------------------------------

   subroutine test_far_field(error)
      !! The melt far-field sampler is the ONLY part of the melt path
      !! that sees the layer stack: the three-equation solve is a pure
      !! function of `(T_far, S_far, u_far, v_far, p_b)`.  So if the
      !! sampler returns bitwise the same quadruple for the filler column
      !! as for the sigma column, the MELT RATE is bitwise the same too —
      !! which is the oracle this slice owes.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 3, NY = 3
      real(wp) :: cover(NX, NY), wet(NX, NY)
      real(wp) :: h_s(NX, NY, NZ_S), t_s(NX, NY, NZ_S), s_s(NX, NY, NZ_S)
      real(wp) :: u_s(NX + 1, NY, NZ_S), v_s(NX, NY + 1, NZ_S)
      real(wp) :: h_z(NX, NY, NZ_Z), t_z(NX, NY, NZ_Z), s_z(NX, NY, NZ_Z)
      real(wp) :: u_z(NX + 1, NY, NZ_Z), v_z(NX, NY + 1, NZ_Z)
      real(wp) :: a1(NX, NY), t1(NX, NY), s1(NX, NY), uu1(NX, NY), vv1(NX, NY)
      real(wp) :: a2(NX, NY), t2(NX, NY), s2(NX, NY), uu2(NX, NY), vv2(NX, NY)
      integer :: k

      cover = 1.0_wp
      wet = 1.0_wp
      ! Sigma stack: 8 live layers, a linear T/S profile so the
      ! thickness-weighted mean is a non-trivial number.
      do k = 1, NZ_S
         h_s(:, :, k) = H_LIVE
         t_s(:, :, k) = H_LIVE*(T0 + 0.05_wp*real(k, wp))
         s_s(:, :, k) = H_LIVE*(S0 + 0.01_wp*real(k, wp))
         u_s(:, :, k) = 0.02_wp*real(k, wp)
         v_s(:, :, k) = -0.01_wp*real(k, wp)
      end do
      ! z_fixed stack: the SAME eight layers with three fillers on top.
      do k = 1, NZ_Z
         if (k > NLIVE) then
            h_z(:, :, k) = H_FILL
            t_z(:, :, k) = 0.0_wp
            s_z(:, :, k) = 0.0_wp
            u_z(:, :, k) = 0.0_wp
            v_z(:, :, k) = 0.0_wp
         else
            h_z(:, :, k) = h_s(1, 1, k)
            t_z(:, :, k) = t_s(1, 1, k)
            s_z(:, :, k) = s_s(1, 1, k)
            u_z(:, :, k) = u_s(1, 1, k)
            v_z(:, :, k) = v_s(1, 1, k)
         end if
      end do

      call cavity_far_field_impl(NX, NY, NZ_S, 120.0_wp, cover, wet, &
                                 h_s, t_s, s_s, u_s, v_s, a1, t1, s1, uu1, vv1)
      call cavity_far_field_impl(NX, NY, NZ_Z, 120.0_wp, cover, wet, &
                                 h_z, t_z, s_z, u_z, v_z, a2, t2, s2, uu2, vv2)

      call check(error, all(a1 == 1.0_wp) .and. all(a2 == 1.0_wp), &
                 "both columns solve")
      if (allocated(error)) return
      call check(error, all(t1 == t2), "T_far is bitwise identical")
      if (allocated(error)) return
      call check(error, all(s1 == s2), "S_far is bitwise identical")
      if (allocated(error)) return
      call check(error, all(uu1 == uu2) .and. all(vv1 == vv2), &
                 "u_far/v_far are bitwise identical")
   end subroutine test_far_field

   subroutine test_deposit_live(error)
      !! The melt heat/salt deposit lands on the first LIVE layer, the
      !! fillers are untouched, and the implied concentration is a sane
      !! ocean number rather than the ~218 degC the module docstring
      !! records for a deposit into a 1.0e-4 m filler.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      real(wp), parameter :: QH = 100.0_wp, QS = 0.0_wp
      real(wp) :: t_top
      integer :: idx_t, kt, k
      logical :: ok

      call build_state(ms, grid, NZ_Z, NFILL)
      call apply_flux(grid, ms, QH, QS)
      idx_t = ms%idx_temperature
      kt = NLIVE

      ok = .true.
      do k = kt + 1, NZ_Z
         if (ms%tracers(idx_t)%hTr(2, 2, k) /= 0.0_wp) ok = .false.
         if (ms%heat_budget_surface(2, 2, k) /= 0.0_wp) ok = .false.
      end do
      call check(error, ok, "no filler layer was written")
      if (allocated(error)) return

      call check(error, ms%heat_budget_surface(2, 2, kt) > 0.0_wp, &
                 "the budget mirror is on the live row")
      if (allocated(error)) return

      t_top = ms%tracers(idx_t)%hTr(2, 2, kt)/ms%h_layer(2, 2, kt)
      call check(error, abs(t_top - T0) < 1.0e-3_wp, &
                 "and the live layer's temperature moved by a sane amount")
      if (allocated(error)) return

      ! The structural half of the old bug: had the deposit gone to nz,
      ! the ALE remap drain's `h_old > H_FLOOR` test would be FALSE
      ! there, so the tracer load would have been written to zero while
      ! the budget kept counting it.
      call check(error,.not. (ms%h_layer(2, 2, NZ_Z) > H_VANISHED), &
                 "k = nz really is a layer the remap drain would empty")
      call ms%destroy()
   end subroutine test_deposit_live

   subroutine test_deposit_matches(error)
      !! Same flux, same live column, three fillers or none: the
      !! increment into the first live layer is BITWISE the same.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: g1, g2
      type(multilayer_state_t) :: m1, m2
      real(wp), parameter :: QH = 100.0_wp, QS = 5.0e-5_wp
      real(wp) :: d_heat_1, d_heat_2, d_salt_1, d_salt_2

      call build_state(m1, g1, NZ_S, 0)
      call apply_flux(g1, m1, QH, QS)
      call build_state(m2, g2, NZ_Z, NFILL)
      call apply_flux(g2, m2, QH, QS)

      d_heat_1 = m1%tracers(m1%idx_temperature)%hTr(2, 2, NLIVE) - H_LIVE*T0
      d_heat_2 = m2%tracers(m2%idx_temperature)%hTr(2, 2, NLIVE) - H_LIVE*T0
      d_salt_1 = m1%tracers(m1%idx_salinity)%hTr(2, 2, NLIVE) - H_LIVE*S0
      d_salt_2 = m2%tracers(m2%idx_salinity)%hTr(2, 2, NLIVE) - H_LIVE*S0

      call check(error, d_heat_1 /= 0.0_wp, "the reference deposit is non-zero")
      if (allocated(error)) return
      call check(error, d_heat_1 == d_heat_2, &
                 "heat increment is bitwise identical with fillers on top")
      if (allocated(error)) return
      call check(error, d_salt_1 == d_salt_2, &
                 "salt increment is bitwise identical with fillers on top")
      call m1%destroy()
      call m2%destroy()
   end subroutine test_deposit_matches

   subroutine test_budget_closes(error)
      !! The budget principle, from step 0: the change in the column's
      !! total tracer load equals the tracked source, to round-off.  With
      !! the deposit on a filler the two would differ by the WHOLE
      !! deposit as soon as the remap ran, while the budget kept counting
      !! it — that is the leak, and this is the assertion that sees it.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(multilayer_state_t) :: ms
      real(wp), parameter :: QH = 100.0_wp, QS = 5.0e-5_wp
      real(wp) :: total_h0, total_s0, total_h1, total_s1
      real(wp) :: src_h, src_s, scale

      call build_state(ms, grid, NZ_Z, NFILL)
      total_h0 = sum(ms%tracers(ms%idx_temperature)%hTr)
      total_s0 = sum(ms%tracers(ms%idx_salinity)%hTr)
      call apply_flux(grid, ms, QH, QS)
      total_h1 = sum(ms%tracers(ms%idx_temperature)%hTr)
      total_s1 = sum(ms%tracers(ms%idx_salinity)%hTr)
      src_h = sum(ms%heat_budget_surface)
      src_s = sum(ms%salt_budget_surface)

      ! Round-off here is set by the TOTAL, not by the increment: the
      ! difference of two ~1e4 K m column totals cannot resolve better
      ! than ~1e-16 of that, whatever the source is worth.
      scale = max(abs(total_h0), abs(src_h))
      call check(error, abs((total_h1 - total_h0) - src_h) <= 1.0e-12_wp*scale, &
                 "heat: delta total = tracked source")
      if (allocated(error)) return
      scale = max(abs(total_s0), abs(src_s))
      call check(error, abs((total_s1 - total_s0) - src_s) <= 1.0e-12_wp*scale, &
                 "salt: delta total = tracked source")
      if (allocated(error)) return
      ! And the source really is the flux, not zero.
      call check(error, src_h > 0.0_wp .and. src_s > 0.0_wp, &
                 "the tracked source is non-zero")
      call ms%destroy()
   end subroutine test_budget_closes

   subroutine test_top_drag_spindown(error)
      !! LINEAR top drag, layer-only mode, under full ice cover: the
      !! first live layer decays at exactly `r` per unit time whether or
      !! not there are fillers above it, so N explicit steps reproduce
      !! `(1 - r*dt)^N` — the same closed form the no-filler column
      !! gives, bitwise.  Before `k_top` the whole stress landed on the
      !! massless filler and the live layer never felt it.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 3, NY = 3, NSTEP = 40
      real(wp), parameter :: R_LIN = 1.0e-4_wp, DTD = 600.0_wp
      real(wp) :: u_s(NX + 1, NY, NZ_S), v_s(NX, NY + 1, NZ_S)
      real(wp) :: u_z(NX + 1, NY, NZ_Z), v_z(NX, NY + 1, NZ_Z)
      real(wp) :: analytic, got_s, got_z
      integer :: n

      call spin(NZ_S, 0, NSTEP, R_LIN, DTD, u_s, v_s)
      call spin(NZ_Z, NFILL, NSTEP, R_LIN, DTD, u_z, v_z)

      analytic = 1.0_wp
      do n = 1, NSTEP
         analytic = analytic*(1.0_wp - R_LIN*DTD)
      end do
      got_s = u_s(2, 2, NLIVE)
      got_z = u_z(2, 2, NLIVE)

      call check(error, abs(got_s - analytic) <= 1.0e-13_wp, &
                 "the no-filler column reproduces the analytic decay")
      if (allocated(error)) return
      call check(error, got_z == got_s, &
                 "and the filler column is bitwise the same")
      if (allocated(error)) return
      call check(error, got_z < 0.99_wp, "the drag actually acted")
      if (allocated(error)) return
      ! The fillers themselves are untouched: no stress on a massless
      ! layer (the `H_VANISHED` band gate, not the old `<= 0` one).
      call check(error, all(u_z(2, 2, NLIVE + 1:NZ_Z) == 1.0_wp), &
                 "the fillers keep their velocity, undragged")
   end subroutine test_top_drag_spindown

   subroutine spin(nz, nfill, nstep, r_lin, dtd, u_face, v_face)
      !! `nstep` explicit top-drag steps on a uniform `u = 1` column.
      integer, intent(in) :: nz, nfill, nstep
      real(wp), intent(in) :: r_lin, dtd
      real(wp), intent(out) :: u_face(4, 3, nz), v_face(3, 4, nz)
      integer, parameter :: NX = 3, NY = 3
      real(wp) :: h(NX, NY, nz), wet(NX, NY)
      real(wp) :: cover_u(NX + 1, NY), cover_v(NX, NY + 1)
      real(wp) :: du(4, 3, nz), dv(3, 4, nz)
      real(wp) :: lam_u(NX + 1, NY), lam_v(NX, NY + 1)
      integer :: ktu(NX + 1, NY), ktv(NX, NY + 1)
      integer :: k, n

      wet = 1.0_wp
      cover_u = 1.0_wp
      cover_v = 1.0_wp
      ktu = nz - nfill
      ktv = nz - nfill
      do k = 1, nz
         if (k > nz - nfill) then
            h(:, :, k) = H_FILL
         else
            h(:, :, k) = H_LIVE
         end if
      end do
      u_face = 1.0_wp
      v_face = 0.0_wp

      do n = 1, nstep
         call top_drag_tendencies_impl(du, dv, lam_u, lam_v, &
                                       u_face, v_face, h, wet, &
                                       cover_u, cover_v, ktu, ktv, H_VANISHED, &
                                       TDRAG_LINEAR, r_lin, 0.0_wp, 1.0e-3_wp, &
                                       0.0_wp, 0.0_wp, 1.0e-3_wp, 0.0_wp, .false., &
                                       4, 3, 3, 4, NX, NY, nz)
         u_face = u_face + dtd*du
         v_face = v_face + dtd*dv
      end do
   end subroutine spin

   subroutine test_fold_row(error)
      !! The implicit-fold rate `lambda_top` is captured on the row the
      !! vdiff diagonal will add it to.  With fillers on top that row is
      !! `k_top`, not `nz` — and the rate must be the plain linear `r`
      !! there (layer-only mode, full cover, `h_in/h_face = 1`).
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 3, NY = 3
      real(wp), parameter :: R_LIN = 2.5e-4_wp
      real(wp) :: h(NX, NY, NZ_Z), wet(NX, NY)
      real(wp) :: u_face(NX + 1, NY, NZ_Z), v_face(NX, NY + 1, NZ_Z)
      real(wp) :: cover_u(NX + 1, NY), cover_v(NX, NY + 1)
      real(wp) :: du(NX + 1, NY, NZ_Z), dv(NX, NY + 1, NZ_Z)
      real(wp) :: lam_u(NX + 1, NY), lam_v(NX, NY + 1)
      integer :: ktu(NX + 1, NY), ktv(NX, NY + 1)
      integer :: k

      wet = 1.0_wp
      cover_u = 1.0_wp
      cover_v = 1.0_wp
      ktu = NLIVE
      ktv = NLIVE
      u_face = 0.5_wp
      v_face = 0.0_wp
      do k = 1, NZ_Z
         if (k > NLIVE) then
            h(:, :, k) = H_FILL
         else
            h(:, :, k) = H_LIVE
         end if
      end do

      call top_drag_tendencies_impl(du, dv, lam_u, lam_v, &
                                    u_face, v_face, h, wet, &
                                    cover_u, cover_v, ktu, ktv, H_VANISHED, &
                                    TDRAG_LINEAR, R_LIN, 0.0_wp, 1.0e-3_wp, &
                                    0.0_wp, 0.0_wp, 1.0e-3_wp, 0.0_wp, .true., &
                                    NX + 1, NY, NX, NY + 1, NX, NY, NZ_Z)

      call check(error, abs(lam_u(2, 2) - R_LIN) <= 1.0e-15_wp, &
                 "the u-face fold rate is the linear r")
      if (allocated(error)) return
      call check(error, abs(lam_v(2, 2) - R_LIN) <= 1.0e-15_wp, &
                 "the v-face fold rate is the linear r")
      if (allocated(error)) return
      call check(error, du(2, 2, NLIVE) /= 0.0_wp, &
                 "and the tendency is on the live row")
      if (allocated(error)) return
      call check(error, all(du(2, 2, NLIVE + 1:NZ_Z) == 0.0_wp), &
                 "with nothing on the fillers")
   end subroutine test_fold_row

end module test_ocean_ktop_consumers
