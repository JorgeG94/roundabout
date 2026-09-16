!! Baroclinic open-boundary kernels for the ocean dyn-core.
module rdb_ocean_obc_baroclinic
   !! Open-boundary baroclinic schemes for the ocean dyn-core.
   !!
   !!   * `ocean_obc_apply_baroclinic` — per-layer normal velocity at open
   !!     faces.  Two radiation schemes via `bc%radiation_scheme`:
   !!       0 (default, "anomaly"): Flather mean + zero-gradient baroclinic
   !!         anomaly (bit-identical to prior behaviour).
   !!       1 ("orlanski"): per-layer implicit-upwind radiation (Orlanski 1976)
   !!         with running-mean phase speed `rx` + optional inflow/outflow
   !!         nudging toward `clamped_u/v` (Marchesiello et al. 2001).
   !!     The continuity renorm at the wall still forces per-layer mass fluxes
   !!     to sum to the BT transport, so mass/eta consistency holds regardless
   !!     of scheme.  RESTART CAVEAT: `rx`/`u_prev` are NOT restart-registered;
   !!     a restart cold-starts the Orlanski scheme with rx = 0 / u_prev = 0.
   !!
   !!   * `ocean_obc_fill_ghosts` — fills h_layer and tracer hTr open-edge
   !!     ghosts (zero-gradient h, upwind-aware hTr), called before
   !!     `continuity_tracer_step_split`.
   !!
   !!   * `ocean_obc_update_reservoirs` — evolves per-edge reservoir
   !!     concentrations `tres` toward interior/external values via implicit
   !!     backward-Euler (Marchesiello et al. 2001).  Called after
   !!     `continuity_tracer_step_split` while mass_flux_* still hold the
   !!     stage's wall-face fluxes.  No-op unless `res_lscale_out/in` > 0;
   !!     when active, `ocean_obc_fill_ghosts` uses the reservoir value.
   !!
   !! GPU rules: explicit-shape dummies; do concurrent j-outer/i-inner;
   !! depth-mean via sequential k-loop in the DC body (nz small); no per-step
   !! allocation; per-tracer loops OUTSIDE the DC kernels (outer-shim); never
   !! deref bc%... inside a DC body (pass rx/u_prev/tres arrays as dummies).
   use rdb_constants, only: wp, H_VANISHED
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_barotropic_workstate, only: barotropic_workstate_t
   use rdb_ocean_boundary_types, only: ocean_bc_state_t, &
                                       OBC_WALL, OBC_OPEN, OBC_TIDAL, OBC_CHAPMAN, &
                                       OBC_NESTED, OBC_CLAMPED, OBC_PERIODIC
   implicit none
   private

   public :: ocean_obc_apply_baroclinic
   public :: ocean_obc_fill_ghosts
   public :: ocean_obc_refill_ghost_ssh
   public :: ocean_obc_update_reservoirs

   ! Minimum depth for outward-normal velocity computation in reservoir update.
   ! Matches the continuity h_min floor (1e-6 m).
   real(wp), parameter :: RES_H_MIN = 1.0e-6_wp

   ! Integer tags for the radiation scheme (mirrors bc%radiation_scheme).
   integer, parameter :: RAD_ANOMALY = 0  !! BT-mean + zero-gradient anomaly
   integer, parameter :: RAD_ORLANSKI = 1  !! Orlanski 1976 per-layer radiation

contains

   subroutine ocean_obc_apply_baroclinic(grid, bc, bt_work, ms, dt)
      !! Set per-layer normal velocity at open-ish faces after
      !! `apply_bt_correction`.  Scheme via `bc%radiation_scheme`:
      !!   RAD_ANOMALY (0, default): Flather mean + zero-gradient anomaly.
      !!   RAD_ORLANSKI (1): per-layer implicit-upwind radiation (Orlanski 1976),
      !!     running-mean phase speed `rx`; `u_prev` snapshot refreshed at END
      !!     so cold-start (u_prev=0) sees rx=0 for a quiescent IC.
      !! Optional nudging (Marchesiello et al. 2001) composes AFTER either
      !! scheme when the selected tau > 0, toward the edge `clamped_u/v`:
      !! tau_in when incoming (dhdt·dhdx ≤ 0 Orlanski; outward vel ≤ 0 anomaly),
      !! else tau_out.  CLAMPED: u_layer(wall) = clamped_u, all layers.
      !! No-op when no edge is open-ish.
      type(hgrid_t), intent(in) :: grid
      type(ocean_bc_state_t), intent(inout) :: bc
      type(barotropic_workstate_t), intent(in) :: bt_work
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: dt
         !! Outer baroclinic timestep (seconds).  Needed for nudging and Orlanski rx.

      integer :: bc_w, bc_e, bc_s, bc_n
      logical :: any_open
      integer :: i_w, i_e, j_s, j_n
      integer :: nz, nx, ny, nxt, nyt
      integer :: rad_scheme

      bc_w = bc%west%bc_type
      bc_e = bc%east%bc_type
      bc_s = bc%south%bc_type
      bc_n = bc%north%bc_type
      ! MPI-seam neutralisation (O0): a seam edge carries no physical BC.
      ! WALL is a no-op throughout this routine (verified: all dispatch is
      ! guarded by is_open_ish/is_radiating, neither of which includes
      ! OBC_WALL), so remapping the cached tag makes every per-edge block
      ! skip the seam.
      if (.not. bc%has_west) bc_w = OBC_WALL
      if (.not. bc%has_east) bc_e = OBC_WALL
      if (.not. bc%has_south) bc_s = OBC_WALL
      if (.not. bc%has_north) bc_n = OBC_WALL

      any_open = is_open_ish(bc_w) .or. is_open_ish(bc_e) .or. &
                 is_open_ish(bc_s) .or. is_open_ish(bc_n)
      if (.not. any_open) return

      nz = ms%nz_ml
      nxt = grid%nx_total
      nyt = grid%ny_total
      nx = grid%nx_phys
      ny = grid%ny_phys
      rad_scheme = bc%radiation_scheme

      ! Physical wall-face indices (u_face_x is nx_total+1 wide).
      i_w = grid%nghost + 1
      i_e = grid%nghost + nx + 1
      j_s = grid%nghost + 1
      j_n = grid%nghost + ny + 1

      if (rad_scheme == RAD_ORLANSKI) then
         ! ---- Orlanski per-layer radiation ----
         ! Zonal edges
         if (is_radiating(bc_w) .and. allocated(bc%rx_west)) then
            call apply_orlanski_west(ms%u_face_x_layer, &
                                     bc%rx_west, bc%u_prev_west, &
                                     nxt, nyt, nz, i_w, j_s, j_n, &
                                     bc%orlanski_rx_max, bc%orlanski_gamma, &
                                     bc%nudge_tau_in, bc%nudge_tau_out, &
                                     dt, bc%west%clamped_u)
         else if (bc_w == OBC_CLAMPED) then
            call apply_clamped_zonal_west(ms%u_face_x_layer, nxt, nyt, nz, i_w, j_s, j_n, &
                                          bc%west%clamped_u)
         end if
         if (is_radiating(bc_e) .and. allocated(bc%rx_east)) then
            call apply_orlanski_east(ms%u_face_x_layer, &
                                     bc%rx_east, bc%u_prev_east, &
                                     nxt, nyt, nz, i_e, j_s, j_n, &
                                     bc%orlanski_rx_max, bc%orlanski_gamma, &
                                     bc%nudge_tau_in, bc%nudge_tau_out, &
                                     dt, bc%east%clamped_u)
         else if (bc_e == OBC_CLAMPED) then
            call apply_clamped_zonal_east(ms%u_face_x_layer, nxt, nyt, nz, i_e, j_s, j_n, &
                                          bc%east%clamped_u)
         end if
         ! Meridional edges
         if (is_radiating(bc_s) .and. allocated(bc%rx_south)) then
            call apply_orlanski_south(ms%v_face_y_layer, &
                                      bc%rx_south, bc%u_prev_south, &
                                      nxt, nyt, nz, j_s, i_w, i_e - 1, &
                                      bc%orlanski_rx_max, bc%orlanski_gamma, &
                                      bc%nudge_tau_in, bc%nudge_tau_out, &
                                      dt, bc%south%clamped_v)
         else if (bc_s == OBC_CLAMPED) then
            call apply_clamped_meridional_south(ms%v_face_y_layer, nxt, nyt, nz, j_s, i_w, i_e - 1, &
                                                bc%south%clamped_v)
         end if
         if (is_radiating(bc_n) .and. allocated(bc%rx_north)) then
            call apply_orlanski_north(ms%v_face_y_layer, &
                                      bc%rx_north, bc%u_prev_north, &
                                      nxt, nyt, nz, j_n, i_w, i_e - 1, &
                                      bc%orlanski_rx_max, bc%orlanski_gamma, &
                                      bc%nudge_tau_in, bc%nudge_tau_out, &
                                      dt, bc%north%clamped_v)
         else if (bc_n == OBC_CLAMPED) then
            call apply_clamped_meridional_north(ms%v_face_y_layer, nxt, nyt, nz, j_n, i_w, i_e - 1, &
                                                bc%north%clamped_v)
         end if

         ! Refresh u_prev snapshots (first-interior-face normal velocity).
         ! Done AFTER all four edge updates so the current call's u_new is
         ! consistent with what the next call will see.
         if (allocated(bc%u_prev_west)) then
            call snapshot_u_prev_west(ms%u_face_x_layer, bc%u_prev_west, &
                                      nxt, nyt, nz, i_w, j_s, j_n)
         end if
         if (allocated(bc%u_prev_east)) then
            call snapshot_u_prev_east(ms%u_face_x_layer, bc%u_prev_east, &
                                      nxt, nyt, nz, i_e, j_s, j_n)
         end if
         if (allocated(bc%u_prev_south)) then
            call snapshot_v_prev_south(ms%v_face_y_layer, bc%u_prev_south, &
                                       nxt, nyt, nz, j_s, i_w, i_e - 1)
         end if
         if (allocated(bc%u_prev_north)) then
            call snapshot_v_prev_north(ms%v_face_y_layer, bc%u_prev_north, &
                                       nxt, nyt, nz, j_n, i_w, i_e - 1)
         end if

      else
         ! ---- Anomaly scheme (default) ----
         ! ---- Zonal open faces (west / east) ----
         if (is_open_ish(bc_w) .or. is_open_ish(bc_e)) then
            call apply_zonal_baroclinic(ms%u_face_x_layer, &
                                        ms%h_layer, &
                                        bt_work%bt_ubt_end, &
                                        nxt, nyt, nz, &
                                        i_w, i_e, j_s, j_n, &
                                        bc_w, bc_e, &
                                        bc%west%clamped_u, bc%east%clamped_u)
         end if

         ! ---- Meridional open faces (south / north) ----
         if (is_open_ish(bc_s) .or. is_open_ish(bc_n)) then
            call apply_meridional_baroclinic(ms%v_face_y_layer, &
                                             ms%h_layer, &
                                             bt_work%bt_vbt_end, &
                                             nxt, nyt, nz, &
                                             i_w, i_e - 1, j_s, j_n, &
                                             bc_s, bc_n, &
                                             bc%south%clamped_v, bc%north%clamped_v)
         end if

         ! ---- Anomaly-scheme nudging (Marchesiello et al. 2001) ----
         ! Applied when either tau > 0.  Sign convention for "incoming":
         ! anomaly scheme uses the outward velocity sign:
         !   west  inflow: u_wall(i_w) > 0  → incoming
         !   east  inflow: u_wall(i_e) < 0  → incoming
         !   south inflow: v_wall(j_s) > 0  → incoming
         !   north inflow: v_wall(j_n) < 0  → incoming
         if (bc%nudge_tau_in > 0.0_wp .or. bc%nudge_tau_out > 0.0_wp) then
            if (is_radiating(bc_w)) then
               call apply_nudge_zonal_west(ms%u_face_x_layer, &
                                           nxt, nyt, nz, i_w, j_s, j_n, &
                                           bc%nudge_tau_in, bc%nudge_tau_out, &
                                           dt, bc%west%clamped_u)
            end if
            if (is_radiating(bc_e)) then
               call apply_nudge_zonal_east(ms%u_face_x_layer, &
                                           nxt, nyt, nz, i_e, j_s, j_n, &
                                           bc%nudge_tau_in, bc%nudge_tau_out, &
                                           dt, bc%east%clamped_u)
            end if
            if (is_radiating(bc_s)) then
               call apply_nudge_meridional_south(ms%v_face_y_layer, &
                                                 nxt, nyt, nz, j_s, i_w, i_e - 1, &
                                                 bc%nudge_tau_in, bc%nudge_tau_out, &
                                                 dt, bc%south%clamped_v)
            end if
            if (is_radiating(bc_n)) then
               call apply_nudge_meridional_north(ms%v_face_y_layer, &
                                                 nxt, nyt, nz, j_n, i_w, i_e - 1, &
                                                 bc%nudge_tau_in, bc%nudge_tau_out, &
                                                 dt, bc%north%clamped_v)
            end if
         end if
      end if

      ! Zero-gradient fill of the per-layer VELOCITY ghosts incl. ghost×ghost
      ! corners: the corner-adjacent C-grid Coriolis/KE/vorticity stencil reads
      ! the diagonal corner ghost, so unfilled corner garbage blows up (open-edge
      ! NaN).  x-then-y; cross-axis range clipped to physical span unless the
      ! adjacent edge is also open, keeping WALL/PERIODIC ghosts untouched.
      call fill_uv_layer_ghosts(ms%u_face_x_layer, ms%v_face_y_layer, &
                                nxt, nyt, nz, grid%nghost, nx, ny, &
                                bc_w, bc_e, bc_s, bc_n)
   end subroutine ocean_obc_apply_baroclinic

   pure subroutine fill_uv_layer_ghosts(u_layer, v_layer, &
                                        nx_total, ny_total, nz, nghost, &
                                        nx_phys, ny_phys, bc_w, bc_e, bc_s, bc_n)
      !! Zero-gradient fill of per-layer face velocities into the OPEN-edge
      !! ghost region, corners included.  u_layer (nx_total+1, ny_total, nz);
      !! v_layer (nx_total, ny_total+1, nz).  x-pass over the full cross extent;
      !! y-pass clipped to the physical span unless the adjacent x-edge is open.
      !! Gated on open-ish tags ⇒ WALL/PERIODIC ⇒ no DC ⇒ bit-identical.
      integer, intent(in) :: nx_total, ny_total, nz, nghost, nx_phys, ny_phys
      real(wp), intent(inout) :: u_layer(nx_total + 1, ny_total, nz)
      real(wp), intent(inout) :: v_layer(nx_total, ny_total + 1, nz)
      integer, intent(in) :: bc_w, bc_e, bc_s, bc_n
      integer :: i, j, k, g, ng
      integer :: iu_lo, iu_hi, iv_lo, iv_hi
      logical :: fw, fe, fs, fn

      ng = nghost
      ! Fill set = SAME as the barotropic-substep η-ghost fill (OPEN/TIDAL/
      ! CHAPMAN/CLAMPED), deliberately NOT `is_open_ish` (excludes NESTED, which
      ! the C-grid barotropic treats as a wall) so velocity and η fills stay
      ! consistent and don't re-trip the corner NaN.
      fw = open_ghost_fill_edge(bc_w)
      fe = open_ghost_fill_edge(bc_e)
      fs = open_ghost_fill_edge(bc_s)
      fn = open_ghost_fill_edge(bc_n)

      ! ---- x-pass: west / east, over the full j extent (covers corners) ----
      if (fw) then
         do concurrent(k=1:nz, j=1:ny_total, g=1:ng)
            u_layer(g, j, k) = u_layer(ng + 1, j, k)             ! u faces 1..ng <- wall face ng+1
         end do
         do concurrent(k=1:nz, j=1:ny_total + 1, g=1:ng)
            v_layer(g, j, k) = v_layer(ng + 1, j, k)             ! v cols 1..ng <- first phys col
         end do
      end if
      if (fe) then
         do concurrent(k=1:nz, j=1:ny_total, g=1:ng)
            u_layer(ng + nx_phys + 1 + g, j, k) = u_layer(ng + nx_phys + 1, j, k)  ! east ghost faces
         end do
         do concurrent(k=1:nz, j=1:ny_total + 1, g=1:ng)
            v_layer(ng + nx_phys + g, j, k) = v_layer(ng + nx_phys, j, k)          ! east ghost cols
         end do
      end if

      ! ---- y-pass: south / north, over the corner-safe i extent ----
      ! Extend the y-pass into an x-ghost column for any NON-PERIODIC x-edge so a
      ! south/north-open edge fills its ghost CORNER even against a WALL side.  A
      ! PERIODIC x-edge is excluded — it owns its ghost via the wrap (preserves
      ! periodic+open mixed configs, e.g. Eady).
      iu_lo = ng + 1                 ! u face i-range (physical span by default)
      iu_hi = ng + nx_phys + 1
      iv_lo = ng + 1                 ! v cell i-range (physical span by default)
      iv_hi = ng + nx_phys
      if (bc_w /= OBC_PERIODIC) then
         iu_lo = 1
         iv_lo = 1
      end if
      if (bc_e /= OBC_PERIODIC) then
         iu_hi = nx_total + 1
         iv_hi = nx_total
      end if
      if (fs) then
         do concurrent(k=1:nz, g=1:ng, i=iu_lo:iu_hi)
            u_layer(i, g, k) = u_layer(i, ng + 1, k)            ! u south ghost rows <- first phys row
         end do
         do concurrent(k=1:nz, g=1:ng, i=iv_lo:iv_hi)
            v_layer(i, g, k) = v_layer(i, ng + 1, k)            ! v south ghost faces <- wall face
         end do
      end if
      if (fn) then
         do concurrent(k=1:nz, g=1:ng, i=iu_lo:iu_hi)
            u_layer(i, ng + ny_phys + g, k) = u_layer(i, ng + ny_phys, k)          ! u north ghost rows
         end do
         do concurrent(k=1:nz, g=1:ng, i=iv_lo:iv_hi)
            v_layer(i, ng + ny_phys + 1 + g, k) = v_layer(i, ng + ny_phys + 1, k)  ! v north ghost faces
         end do
      end if
   end subroutine fill_uv_layer_ghosts

   pure logical function open_ghost_fill_edge(bc_type) result(res)
      !! Edge types that receive the open zero-gradient ghost fill (velocity
      !! here, η in the barotropic substep — kept identical on purpose).
      !! OPEN / TIDAL / CHAPMAN / CLAMPED; NOT NESTED (barotropic-wall path).
      integer, intent(in) :: bc_type
      res = (bc_type == OBC_OPEN .or. bc_type == OBC_TIDAL .or. &
             bc_type == OBC_CHAPMAN .or. bc_type == OBC_CLAMPED)
   end function open_ghost_fill_edge

   pure logical function is_open_ish(bc_type) result(res)
      !! Returns .true. for edge types that need baroclinic velocity treatment.
      integer, intent(in) :: bc_type
      res = (bc_type == OBC_OPEN .or. bc_type == OBC_TIDAL .or. &
             bc_type == OBC_CHAPMAN .or. bc_type == OBC_NESTED .or. &
             bc_type == OBC_CLAMPED)
   end function is_open_ish

   pure logical function is_radiating(bc_type) result(res)
      !! Returns .true. for Flather-class edges (not CLAMPED).
      integer, intent(in) :: bc_type
      res = (bc_type == OBC_OPEN .or. bc_type == OBC_TIDAL .or. &
             bc_type == OBC_CHAPMAN .or. bc_type == OBC_NESTED)
   end function is_radiating

   ! -----------------------------------------------------------------
   ! Orlanski per-layer radiation kernels (Orlanski 1976).
   ! Per edge: 1st/2nd interior faces are stepped inward from the wall;
   !   dhdt = u_prev - u_new (1st int), dhdx = u(1st int) - u(2nd int) (toward
   !   the boundary).  rx_raw = min(dhdt/dhdx, rx_max) when dhdt·dhdx > 0
   !   (outgoing), else 0 (incoming, no radiation).
   ! Nudging (Marchesiello et al. 2001): tau = tau_in when dhdt·dhdx <= 0
   !   (incoming) else tau_out; g2 = dt/(tau+dt) (skip tau<=0);
   !   u_b = (1-g2)*u_b + g2*u_data.
   ! -----------------------------------------------------------------

   pure subroutine apply_orlanski_west(u_layer, &
                                       rx, u_prev, &
                                       nx_total, ny_total, nz, &
                                       i_w, j0, j1, &
                                       rx_max, gamma_u, &
                                       tau_in, tau_out, dt, u_data)
      !! Orlanski radiation for the west open edge (Orlanski 1976).
      !! Interior face index I = i_w+1 (1st), I-1 = i_w+2 (2nd).
      !! Outward normal = -x.  dhdx = u(i_w+1) - u(i_w+2) (westward gradient).
      !! rx updated in-place; u_prev NOT updated here (done by the caller
      !! after all edges are set).
      integer, intent(in) :: nx_total, ny_total, nz
      real(wp), intent(inout) :: u_layer(nx_total + 1, ny_total, nz)
      real(wp), intent(inout) :: rx(ny_total, nz)
      real(wp), intent(in)    :: u_prev(ny_total, nz)
      integer, intent(in) :: i_w, j0, j1
      real(wp), intent(in) :: rx_max, gamma_u, tau_in, tau_out, dt, u_data

      integer  :: j, k
      real(wp) :: u_int_1
      real(wp) :: dhdt_jk, dhdx_jk, rx_raw, rx_new, u_b
      real(wp) :: tau, g2

      do concurrent(j=j0:j1) &
         local(k, u_int_1, dhdt_jk, dhdx_jk, rx_raw, rx_new, u_b, tau, g2)

         do k = 1, nz
            u_int_1 = u_layer(i_w + 1, j, k)  ! first interior face velocity

            ! West edge: outward = -x; dhdx = u(i_w+1) - u(i_w+2).
            dhdt_jk = u_prev(j, k) - u_int_1
            dhdx_jk = u_int_1 - u_layer(i_w + 2, j, k)

            if (dhdt_jk*dhdx_jk > 0.0_wp) then
               rx_raw = min(dhdt_jk/dhdx_jk, rx_max)
            else
               rx_raw = 0.0_wp
            end if

            ! Running-mean update (gamma_u = 1 = instant, no memory).
            rx_new = (1.0_wp - gamma_u)*rx(j, k) + gamma_u*rx_raw
            rx(j, k) = rx_new

            ! Semi-implicit boundary update: u_b = (u_b_old + rx*u_int_1)/(1+rx).
            ! Anchors on the wall face's OWN per-layer value, so rx=0 leaves the
            ! baroclinic structure untouched (anchoring on BT velocity would
            ! collapse the boundary shear every quiet phase).
            u_b = (u_layer(i_w, j, k) + rx_new*u_int_1)/(1.0_wp + rx_new)

            ! Optional nudging: "incoming" = dhdt*dhdx <= 0 ⟹ tau_in.
            if (dhdt_jk*dhdx_jk <= 0.0_wp) then
               tau = tau_in
            else
               tau = tau_out
            end if
            if (tau > 0.0_wp) then
               g2 = dt/(tau + dt)
               u_b = (1.0_wp - g2)*u_b + g2*u_data
            end if

            u_layer(i_w, j, k) = u_b
         end do
      end do
   end subroutine apply_orlanski_west

   pure subroutine apply_orlanski_east(u_layer, &
                                       rx, u_prev, &
                                       nx_total, ny_total, nz, &
                                       i_e, j0, j1, &
                                       rx_max, gamma_u, &
                                       tau_in, tau_out, dt, u_data)
      !! Orlanski radiation for the east open edge (Orlanski 1976).
      !! Interior face index I = i_e-1 (1st), I-1 = i_e-2 (2nd).
      !! Outward normal = +x.  dhdx = u(i_e-1) - u(i_e-2) (eastward gradient).
      integer, intent(in) :: nx_total, ny_total, nz
      real(wp), intent(inout) :: u_layer(nx_total + 1, ny_total, nz)
      real(wp), intent(inout) :: rx(ny_total, nz)
      real(wp), intent(in)    :: u_prev(ny_total, nz)
      integer, intent(in) :: i_e, j0, j1
      real(wp), intent(in) :: rx_max, gamma_u, tau_in, tau_out, dt, u_data

      integer  :: j, k
      real(wp) :: u_int_1
      real(wp) :: dhdt_jk, dhdx_jk, rx_raw, rx_new, u_b
      real(wp) :: tau, g2

      do concurrent(j=j0:j1) &
         local(k, u_int_1, dhdt_jk, dhdx_jk, rx_raw, rx_new, u_b, tau, g2)

         do k = 1, nz
            u_int_1 = u_layer(i_e - 1, j, k)

            ! East edge: outward = +x; dhdx = u(i_e-1) - u(i_e-2).
            dhdt_jk = u_prev(j, k) - u_int_1
            dhdx_jk = u_int_1 - u_layer(i_e - 2, j, k)

            if (dhdt_jk*dhdx_jk > 0.0_wp) then
               rx_raw = min(dhdt_jk/dhdx_jk, rx_max)
            else
               rx_raw = 0.0_wp
            end if

            rx_new = (1.0_wp - gamma_u)*rx(j, k) + gamma_u*rx_raw
            rx(j, k) = rx_new

            ! Per-layer anchor — see the west kernel's comment.
            u_b = (u_layer(i_e, j, k) + rx_new*u_int_1)/(1.0_wp + rx_new)

            if (dhdt_jk*dhdx_jk <= 0.0_wp) then
               tau = tau_in
            else
               tau = tau_out
            end if
            if (tau > 0.0_wp) then
               g2 = dt/(tau + dt)
               u_b = (1.0_wp - g2)*u_b + g2*u_data
            end if

            u_layer(i_e, j, k) = u_b
         end do
      end do
   end subroutine apply_orlanski_east

   pure subroutine apply_orlanski_south(v_layer, &
                                        rx, v_prev, &
                                        nx_total, ny_total, nz, &
                                        j_s, i0, i1, &
                                        rx_max, gamma_u, &
                                        tau_in, tau_out, dt, v_data)
      !! Orlanski radiation for the south open edge (Orlanski 1976).
      !! Interior face J = j_s+1 (1st), J-1 = j_s+2 (2nd).
      !! Outward normal = -y.  dhdx = v(j_s+1) - v(j_s+2) (southward gradient).
      integer, intent(in) :: nx_total, ny_total, nz
      real(wp), intent(inout) :: v_layer(nx_total, ny_total + 1, nz)
      real(wp), intent(inout) :: rx(nx_total, nz)
      real(wp), intent(in)    :: v_prev(nx_total, nz)
      integer, intent(in) :: j_s, i0, i1
      real(wp), intent(in) :: rx_max, gamma_u, tau_in, tau_out, dt, v_data

      integer  :: i, k
      real(wp) :: v_int_1
      real(wp) :: dhdt_ik, dhdx_ik, rx_raw, rx_new, v_b
      real(wp) :: tau, g2

      do concurrent(i=i0:i1) &
         local(k, v_int_1, dhdt_ik, dhdx_ik, rx_raw, rx_new, v_b, tau, g2)

         do k = 1, nz
            v_int_1 = v_layer(i, j_s + 1, k)

            ! South edge: outward = -y; dhdx = v(j_s+1) - v(j_s+2).
            dhdt_ik = v_prev(i, k) - v_int_1
            dhdx_ik = v_int_1 - v_layer(i, j_s + 2, k)

            if (dhdt_ik*dhdx_ik > 0.0_wp) then
               rx_raw = min(dhdt_ik/dhdx_ik, rx_max)
            else
               rx_raw = 0.0_wp
            end if

            rx_new = (1.0_wp - gamma_u)*rx(i, k) + gamma_u*rx_raw
            rx(i, k) = rx_new

            ! Per-layer anchor — see the west kernel's comment.
            v_b = (v_layer(i, j_s, k) + rx_new*v_int_1)/(1.0_wp + rx_new)

            if (dhdt_ik*dhdx_ik <= 0.0_wp) then
               tau = tau_in
            else
               tau = tau_out
            end if
            if (tau > 0.0_wp) then
               g2 = dt/(tau + dt)
               v_b = (1.0_wp - g2)*v_b + g2*v_data
            end if

            v_layer(i, j_s, k) = v_b
         end do
      end do
   end subroutine apply_orlanski_south

   pure subroutine apply_orlanski_north(v_layer, &
                                        rx, v_prev, &
                                        nx_total, ny_total, nz, &
                                        j_n, i0, i1, &
                                        rx_max, gamma_u, &
                                        tau_in, tau_out, dt, v_data)
      !! Orlanski radiation for the north open edge (Orlanski 1976).
      !! Interior face J = j_n-1 (1st), J-1 = j_n-2 (2nd).
      !! Outward normal = +y.  dhdx = v(j_n-1) - v(j_n-2) (northward gradient).
      integer, intent(in) :: nx_total, ny_total, nz
      real(wp), intent(inout) :: v_layer(nx_total, ny_total + 1, nz)
      real(wp), intent(inout) :: rx(nx_total, nz)
      real(wp), intent(in)    :: v_prev(nx_total, nz)
      integer, intent(in) :: j_n, i0, i1
      real(wp), intent(in) :: rx_max, gamma_u, tau_in, tau_out, dt, v_data

      integer  :: i, k
      real(wp) :: v_int_1
      real(wp) :: dhdt_ik, dhdx_ik, rx_raw, rx_new, v_b
      real(wp) :: tau, g2

      do concurrent(i=i0:i1) &
         local(k, v_int_1, dhdt_ik, dhdx_ik, rx_raw, rx_new, v_b, tau, g2)

         do k = 1, nz
            v_int_1 = v_layer(i, j_n - 1, k)

            ! North edge: outward = +y; dhdx = v(j_n-1) - v(j_n-2).
            dhdt_ik = v_prev(i, k) - v_int_1
            dhdx_ik = v_int_1 - v_layer(i, j_n - 2, k)

            if (dhdt_ik*dhdx_ik > 0.0_wp) then
               rx_raw = min(dhdt_ik/dhdx_ik, rx_max)
            else
               rx_raw = 0.0_wp
            end if

            rx_new = (1.0_wp - gamma_u)*rx(i, k) + gamma_u*rx_raw
            rx(i, k) = rx_new

            ! Per-layer anchor — see the west kernel's comment.
            v_b = (v_layer(i, j_n, k) + rx_new*v_int_1)/(1.0_wp + rx_new)

            if (dhdt_ik*dhdx_ik <= 0.0_wp) then
               tau = tau_in
            else
               tau = tau_out
            end if
            if (tau > 0.0_wp) then
               g2 = dt/(tau + dt)
               v_b = (1.0_wp - g2)*v_b + g2*v_data
            end if

            v_layer(i, j_n, k) = v_b
         end do
      end do
   end subroutine apply_orlanski_north

   ! -----------------------------------------------------------------
   ! CLAMPED helpers for the Orlanski dispatch path (simple wrappers)
   ! -----------------------------------------------------------------

   pure subroutine apply_clamped_zonal_west(u_layer, nx_total, ny_total, nz, &
                                            i_w, j0, j1, clamp_u)
      integer, intent(in) :: nx_total, ny_total, nz, i_w, j0, j1
      real(wp), intent(inout) :: u_layer(nx_total + 1, ny_total, nz)
      real(wp), intent(in) :: clamp_u
      integer :: j, k
      do concurrent(k=1:nz, j=j0:j1)
         u_layer(i_w, j, k) = clamp_u
      end do
   end subroutine apply_clamped_zonal_west

   pure subroutine apply_clamped_zonal_east(u_layer, nx_total, ny_total, nz, &
                                            i_e, j0, j1, clamp_u)
      integer, intent(in) :: nx_total, ny_total, nz, i_e, j0, j1
      real(wp), intent(inout) :: u_layer(nx_total + 1, ny_total, nz)
      real(wp), intent(in) :: clamp_u
      integer :: j, k
      do concurrent(k=1:nz, j=j0:j1)
         u_layer(i_e, j, k) = clamp_u
      end do
   end subroutine apply_clamped_zonal_east

   pure subroutine apply_clamped_meridional_south(v_layer, nx_total, ny_total, nz, &
                                                  j_s, i0, i1, clamp_v)
      integer, intent(in) :: nx_total, ny_total, nz, j_s, i0, i1
      real(wp), intent(inout) :: v_layer(nx_total, ny_total + 1, nz)
      real(wp), intent(in) :: clamp_v
      integer :: i, k
      do concurrent(k=1:nz, i=i0:i1)
         v_layer(i, j_s, k) = clamp_v
      end do
   end subroutine apply_clamped_meridional_south

   pure subroutine apply_clamped_meridional_north(v_layer, nx_total, ny_total, nz, &
                                                  j_n, i0, i1, clamp_v)
      integer, intent(in) :: nx_total, ny_total, nz, j_n, i0, i1
      real(wp), intent(inout) :: v_layer(nx_total, ny_total + 1, nz)
      real(wp), intent(in) :: clamp_v
      integer :: i, k
      do concurrent(k=1:nz, i=i0:i1)
         v_layer(i, j_n, k) = clamp_v
      end do
   end subroutine apply_clamped_meridional_north

   ! -----------------------------------------------------------------
   ! u_prev snapshot helpers (copy first-interior-face velocity)
   ! -----------------------------------------------------------------

   pure subroutine snapshot_u_prev_west(u_layer, u_prev, nx_total, ny_total, nz, &
                                        i_w, j0, j1)
      !! Snapshot u at the first interior face (i_w+1) into u_prev_west.
      integer, intent(in) :: nx_total, ny_total, nz, i_w, j0, j1
      real(wp), intent(in)    :: u_layer(nx_total + 1, ny_total, nz)
      real(wp), intent(inout) :: u_prev(ny_total, nz)
      integer :: j, k
      do concurrent(k=1:nz, j=j0:j1)
         u_prev(j, k) = u_layer(i_w + 1, j, k)
      end do
   end subroutine snapshot_u_prev_west

   pure subroutine snapshot_u_prev_east(u_layer, u_prev, nx_total, ny_total, nz, &
                                        i_e, j0, j1)
      !! Snapshot u at the first interior face (i_e-1) into u_prev_east.
      integer, intent(in) :: nx_total, ny_total, nz, i_e, j0, j1
      real(wp), intent(in)    :: u_layer(nx_total + 1, ny_total, nz)
      real(wp), intent(inout) :: u_prev(ny_total, nz)
      integer :: j, k
      do concurrent(k=1:nz, j=j0:j1)
         u_prev(j, k) = u_layer(i_e - 1, j, k)
      end do
   end subroutine snapshot_u_prev_east

   pure subroutine snapshot_v_prev_south(v_layer, v_prev, nx_total, ny_total, nz, &
                                         j_s, i0, i1)
      !! Snapshot v at the first interior face (j_s+1) into u_prev_south.
      integer, intent(in) :: nx_total, ny_total, nz, j_s, i0, i1
      real(wp), intent(in)    :: v_layer(nx_total, ny_total + 1, nz)
      real(wp), intent(inout) :: v_prev(nx_total, nz)
      integer :: i, k
      do concurrent(k=1:nz, i=i0:i1)
         v_prev(i, k) = v_layer(i, j_s + 1, k)
      end do
   end subroutine snapshot_v_prev_south

   pure subroutine snapshot_v_prev_north(v_layer, v_prev, nx_total, ny_total, nz, &
                                         j_n, i0, i1)
      !! Snapshot v at the first interior face (j_n-1) into u_prev_north.
      integer, intent(in) :: nx_total, ny_total, nz, j_n, i0, i1
      real(wp), intent(in)    :: v_layer(nx_total, ny_total + 1, nz)
      real(wp), intent(inout) :: v_prev(nx_total, nz)
      integer :: i, k
      do concurrent(k=1:nz, i=i0:i1)
         v_prev(i, k) = v_layer(i, j_n - 1, k)
      end do
   end subroutine snapshot_v_prev_north

   ! -----------------------------------------------------------------
   ! Nudging-only helpers for the anomaly-scheme path
   ! -----------------------------------------------------------------

   pure subroutine apply_nudge_zonal_west(u_layer, nx_total, ny_total, nz, &
                                          i_w, j0, j1, tau_in, tau_out, dt, u_data)
      !! Post-anomaly nudging for the west edge.
      !! "Incoming" at west = u_wall > 0 (eastward, into domain) ⟹ tau_in.
      integer, intent(in) :: nx_total, ny_total, nz, i_w, j0, j1
      real(wp), intent(inout) :: u_layer(nx_total + 1, ny_total, nz)
      real(wp), intent(in) :: tau_in, tau_out, dt, u_data
      integer :: j, k
      real(wp) :: tau, g2, u_b
      do concurrent(j=j0:j1, k=1:nz) local(tau, g2, u_b)
         if (u_layer(i_w, j, k) > 0.0_wp) then
            tau = tau_in
         else
            tau = tau_out
         end if
         if (tau > 0.0_wp) then
            g2 = dt/(tau + dt)
            u_b = (1.0_wp - g2)*u_layer(i_w, j, k) + g2*u_data
            u_layer(i_w, j, k) = u_b
         end if
      end do
   end subroutine apply_nudge_zonal_west

   pure subroutine apply_nudge_zonal_east(u_layer, nx_total, ny_total, nz, &
                                          i_e, j0, j1, tau_in, tau_out, dt, u_data)
      !! Post-anomaly nudging for the east edge.
      !! "Incoming" at east = u_wall < 0 (westward, into domain) ⟹ tau_in.
      integer, intent(in) :: nx_total, ny_total, nz, i_e, j0, j1
      real(wp), intent(inout) :: u_layer(nx_total + 1, ny_total, nz)
      real(wp), intent(in) :: tau_in, tau_out, dt, u_data
      integer :: j, k
      real(wp) :: tau, g2, u_b
      do concurrent(j=j0:j1, k=1:nz) local(tau, g2, u_b)
         if (u_layer(i_e, j, k) < 0.0_wp) then
            tau = tau_in
         else
            tau = tau_out
         end if
         if (tau > 0.0_wp) then
            g2 = dt/(tau + dt)
            u_b = (1.0_wp - g2)*u_layer(i_e, j, k) + g2*u_data
            u_layer(i_e, j, k) = u_b
         end if
      end do
   end subroutine apply_nudge_zonal_east

   pure subroutine apply_nudge_meridional_south(v_layer, nx_total, ny_total, nz, &
                                                j_s, i0, i1, tau_in, tau_out, dt, v_data)
      !! Post-anomaly nudging for the south edge.
      !! "Incoming" at south = v_wall > 0 ⟹ tau_in.
      integer, intent(in) :: nx_total, ny_total, nz, j_s, i0, i1
      real(wp), intent(inout) :: v_layer(nx_total, ny_total + 1, nz)
      real(wp), intent(in) :: tau_in, tau_out, dt, v_data
      integer :: i, k
      real(wp) :: tau, g2, v_b
      do concurrent(i=i0:i1, k=1:nz) local(tau, g2, v_b)
         if (v_layer(i, j_s, k) > 0.0_wp) then
            tau = tau_in
         else
            tau = tau_out
         end if
         if (tau > 0.0_wp) then
            g2 = dt/(tau + dt)
            v_b = (1.0_wp - g2)*v_layer(i, j_s, k) + g2*v_data
            v_layer(i, j_s, k) = v_b
         end if
      end do
   end subroutine apply_nudge_meridional_south

   pure subroutine apply_nudge_meridional_north(v_layer, nx_total, ny_total, nz, &
                                                j_n, i0, i1, tau_in, tau_out, dt, v_data)
      !! Post-anomaly nudging for the north edge.
      !! "Incoming" at north = v_wall < 0 ⟹ tau_in.
      integer, intent(in) :: nx_total, ny_total, nz, j_n, i0, i1
      real(wp), intent(inout) :: v_layer(nx_total, ny_total + 1, nz)
      real(wp), intent(in) :: tau_in, tau_out, dt, v_data
      integer :: i, k
      real(wp) :: tau, g2, v_b
      do concurrent(i=i0:i1, k=1:nz) local(tau, g2, v_b)
         if (v_layer(i, j_n, k) < 0.0_wp) then
            tau = tau_in
         else
            tau = tau_out
         end if
         if (tau > 0.0_wp) then
            g2 = dt/(tau + dt)
            v_b = (1.0_wp - g2)*v_layer(i, j_n, k) + g2*v_data
            v_layer(i, j_n, k) = v_b
         end if
      end do
   end subroutine apply_nudge_meridional_north

   pure subroutine apply_zonal_baroclinic(u_layer, h_layer, ubt_end, &
                                          nx_total, ny_total, nz, &
                                          i_w, i_e, j0, j1, &
                                          bc_w, bc_e, clamp_u_w, clamp_u_e)
      !! Set west and east wall-face per-layer u to
      !! Flather mean + zero-gradient baroclinic anomaly (radiating edges)
      !! or uniform clamped_u (CLAMPED edge).
      !! Explicit-shape dummies; one do concurrent per edge (j outer, k inner).
      integer, intent(in) :: nx_total, ny_total, nz
      real(wp), intent(inout) :: u_layer(nx_total + 1, ny_total, nz)
         !! u_face_x_layer, shape (nx_total+1, ny_total, nz).
      real(wp), intent(in)    :: h_layer(nx_total, ny_total, nz)
      real(wp), intent(in)    :: ubt_end(nx_total + 1, ny_total)
         !! bt_ubt_end, shape (nx_total+1, ny_total).
      integer, intent(in) :: i_w, i_e   !! west / east physical wall-face indices
      integer, intent(in) :: j0, j1     !! first / last physical j-cell
      integer, intent(in) :: bc_w, bc_e
      real(wp), intent(in) :: clamp_u_w, clamp_u_e

      integer  :: j, k
      real(wp) :: h_tot, u_bar, u_int
      integer  :: i_int_w, i_int_e  ! first interior face index (one step inside)

      ! Interior face = one face step inward from the wall face.
      ! For west wall at i_w: interior face is i_w+1 (reads h_layer cell i_w).
      ! For east wall at i_e: interior face is i_e-1 (reads h_layer cell i_e-1).
      i_int_w = i_w + 1
      i_int_e = i_e - 1

      ! ---- West wall ----
      if (bc_w == OBC_CLAMPED) then
         do concurrent(k=1:nz, j=j0:j1)
            u_layer(i_w, j, k) = clamp_u_w
         end do
      else if (is_radiating(bc_w)) then
         do concurrent(j=j0:j1) local(k, h_tot, u_bar, u_int)
            ! Sequential k-loop for depth mean (nz is small, no reduction needed).
            h_tot = 0.0_wp
            u_bar = 0.0_wp
            do k = 1, nz
               ! Face i_int_w is between cell i_w and cell i_w+1.
               ! Use the upstream cell (i_w) thickness for the wall side.
               h_tot = h_tot + h_layer(i_w, j, k)
               u_bar = u_bar + u_layer(i_int_w, j, k)*h_layer(i_w, j, k)
            end do
            if (h_tot > 0.0_wp) then
               u_bar = u_bar/h_tot
            else
               u_bar = 0.0_wp
            end if
            do k = 1, nz
               u_int = u_layer(i_int_w, j, k)
               u_layer(i_w, j, k) = ubt_end(i_w, j) + (u_int - u_bar)
            end do
         end do
      end if

      ! ---- East wall ----
      if (bc_e == OBC_CLAMPED) then
         do concurrent(k=1:nz, j=j0:j1)
            u_layer(i_e, j, k) = clamp_u_e
         end do
      else if (is_radiating(bc_e)) then
         do concurrent(j=j0:j1) local(k, h_tot, u_bar, u_int)
            h_tot = 0.0_wp
            u_bar = 0.0_wp
            do k = 1, nz
               ! Use cell (i_e-1) thickness — the last interior cell.
               h_tot = h_tot + h_layer(i_e - 1, j, k)
               u_bar = u_bar + u_layer(i_int_e, j, k)*h_layer(i_e - 1, j, k)
            end do
            if (h_tot > 0.0_wp) then
               u_bar = u_bar/h_tot
            else
               u_bar = 0.0_wp
            end if
            do k = 1, nz
               u_int = u_layer(i_int_e, j, k)
               u_layer(i_e, j, k) = ubt_end(i_e, j) + (u_int - u_bar)
            end do
         end do
      end if
   end subroutine apply_zonal_baroclinic

   pure subroutine apply_meridional_baroclinic(v_layer, h_layer, vbt_end, &
                                               nx_total, ny_total, nz, &
                                               i0, i1, j_s, j_n, &
                                               bc_s, bc_n, clamp_v_s, clamp_v_n)
      !! Set south and north wall-face per-layer v.  Mirror of
      !! `apply_zonal_baroclinic` for the y-direction.
      integer, intent(in) :: nx_total, ny_total, nz
      real(wp), intent(inout) :: v_layer(nx_total, ny_total + 1, nz)
         !! v_face_y_layer, shape (nx_total, ny_total+1, nz).
      real(wp), intent(in)    :: h_layer(nx_total, ny_total, nz)
      real(wp), intent(in)    :: vbt_end(nx_total, ny_total + 1)
      integer, intent(in) :: i0, i1     !! first / last physical i-cell
      integer, intent(in) :: j_s, j_n   !! south / north physical wall-face indices
      integer, intent(in) :: bc_s, bc_n
      real(wp), intent(in) :: clamp_v_s, clamp_v_n

      integer  :: i, k
      real(wp) :: h_tot, v_bar, v_int
      integer  :: j_int_s, j_int_n

      j_int_s = j_s + 1
      j_int_n = j_n - 1

      ! ---- South wall ----
      if (bc_s == OBC_CLAMPED) then
         do concurrent(k=1:nz, i=i0:i1)
            v_layer(i, j_s, k) = clamp_v_s
         end do
      else if (is_radiating(bc_s)) then
         do concurrent(i=i0:i1) local(k, h_tot, v_bar, v_int)
            h_tot = 0.0_wp
            v_bar = 0.0_wp
            do k = 1, nz
               h_tot = h_tot + h_layer(i, j_s, k)
               v_bar = v_bar + v_layer(i, j_int_s, k)*h_layer(i, j_s, k)
            end do
            if (h_tot > 0.0_wp) then
               v_bar = v_bar/h_tot
            else
               v_bar = 0.0_wp
            end if
            do k = 1, nz
               v_int = v_layer(i, j_int_s, k)
               v_layer(i, j_s, k) = vbt_end(i, j_s) + (v_int - v_bar)
            end do
         end do
      end if

      ! ---- North wall ----
      if (bc_n == OBC_CLAMPED) then
         do concurrent(k=1:nz, i=i0:i1)
            v_layer(i, j_n, k) = clamp_v_n
         end do
      else if (is_radiating(bc_n)) then
         do concurrent(i=i0:i1) local(k, h_tot, v_bar, v_int)
            h_tot = 0.0_wp
            v_bar = 0.0_wp
            do k = 1, nz
               h_tot = h_tot + h_layer(i, j_n - 1, k)
               v_bar = v_bar + v_layer(i, j_int_n, k)*h_layer(i, j_n - 1, k)
            end do
            if (h_tot > 0.0_wp) then
               v_bar = v_bar/h_tot
            else
               v_bar = 0.0_wp
            end if
            do k = 1, nz
               v_int = v_layer(i, j_int_n, k)
               v_layer(i, j_n, k) = vbt_end(i, j_n) + (v_int - v_bar)
            end do
         end do
      end if
   end subroutine apply_meridional_baroclinic

   ! Ghost fills at open edges (h_layer + tracer hTr)

   subroutine ocean_obc_fill_ghosts(grid, bc, ms)
      !! Fill h_layer and tracer hTr ghosts at open-ish edges.
      !! h_layer: zero-gradient (copy adjacent interior column).
      !! Tracer hTr per (j,k): outflow ⇒ ghost := interior (zero-gradient);
      !! inflow ⇒ ghost hTr := clamped_tracer(it) * h_ghost.  Inflow criterion
      !! uses the per-layer wall-face velocity sign (outward-normal convention):
      !! west inflow u_wall>0, east u_wall<0, south v_wall>0, north v_wall<0.
      !! No-op when no edge is open-ish.  Per-tracer loop outside the DCs.
      type(hgrid_t), intent(in) :: grid
      type(ocean_bc_state_t), intent(in) :: bc
      type(multilayer_state_t), intent(inout) :: ms

      integer :: bc_w, bc_e, bc_s, bc_n
      logical :: any_open
      integer :: it, nz, nxt, nyt, nx, ny
      integer :: i_w, i_e, j_s, j_n

      bc_w = bc%west%bc_type
      bc_e = bc%east%bc_type
      bc_s = bc%south%bc_type
      bc_n = bc%north%bc_type
      ! MPI-seam neutralisation (O0): a seam edge carries no physical BC.
      ! WALL is a no-op throughout this routine (verified: all dispatch is
      ! guarded by is_open_ish, which excludes OBC_WALL), so remapping the
      ! cached tag makes every per-edge block skip the seam.
      if (.not. bc%has_west) bc_w = OBC_WALL
      if (.not. bc%has_east) bc_e = OBC_WALL
      if (.not. bc%has_south) bc_s = OBC_WALL
      if (.not. bc%has_north) bc_n = OBC_WALL

      any_open = is_open_ish(bc_w) .or. is_open_ish(bc_e) .or. &
                 is_open_ish(bc_s) .or. is_open_ish(bc_n)
      if (.not. any_open) return

      nz = ms%nz_ml
      nxt = grid%nx_total
      nyt = grid%ny_total
      nx = grid%nx_phys
      ny = grid%ny_phys
      i_w = grid%nghost + 1          ! west physical wall-face / first interior cell
      i_e = grid%nghost + nx         ! last interior cell (east)
      j_s = grid%nghost + 1          ! south physical wall-face / first interior cell
      j_n = grid%nghost + ny         ! last interior cell (north)

      ! ---- h_layer ghosts: zero-gradient ----
      call fill_h_layer_ghosts(ms%h_layer, nxt, nyt, nz, &
                               i_w, i_e, j_s, j_n, grid%nghost, &
                               bc_w, bc_e, bc_s, bc_n)

      ! ---- Tracer hTr ghosts: upwind-aware or reservoir-based ----
      ! Reservoir active (tres_* allocated): unconditional hTr_ghost = tres*h_ghost
      ! (Marchesiello et al. 2001).  Absent (default): sign-switch path,
      ! bit-identical to prior behaviour.  Dispatch is per-edge-per-direction.
      if (allocated(ms%tracers)) then
         do it = 1, size(ms%tracers)
            if (.not. allocated(ms%tracers(it)%hTr)) cycle

            ! Zero-gradient CONCENTRATION pre-fill over the FULL open-edge ghost
            ! region (corners included): the per-edge upwind fills below only
            ! touch the physical cross-extent, so unfilled ghost×ghost corners
            ! would accumulate garbage under inflow and blow up T/S.  Upwind
            ! fills then overwrite the physical edges, leaving corners zero-grad.
            call fill_tracer_ghosts_zerograd(ms%tracers(it)%hTr, ms%h_layer, &
                                             nxt, nyt, nz, i_w, i_e, j_s, j_n, &
                                             grid%nghost, bc_w, bc_e, bc_s, bc_n)

            ! ---- West ghost fill ----
            if (is_open_ish(bc_w)) then
               if (allocated(bc%tres_west)) then
                  call fill_ghost_west_res( &
                     ms%tracers(it)%hTr, ms%h_layer, bc%tres_west, &
                     nxt, nyt, nz, i_w, j_s, j_n, grid%nghost, it, size(bc%tres_west, 3))
               else
                  call fill_ghost_west_sign(ms%tracers(it)%hTr, ms%h_layer, &
                                            ms%u_face_x_layer, &
                                            nxt, nyt, nz, i_w, j_s, j_n, grid%nghost, &
                                            get_clamped_tracer(bc%west%clamped_tracer, it))
               end if
            end if

            ! ---- East ghost fill ----
            if (is_open_ish(bc_e)) then
               if (allocated(bc%tres_east)) then
                  call fill_ghost_east_res( &
                     ms%tracers(it)%hTr, ms%h_layer, bc%tres_east, &
                     nxt, nyt, nz, i_e, j_s, j_n, grid%nghost, it, size(bc%tres_east, 3))
               else
                  call fill_ghost_east_sign(ms%tracers(it)%hTr, ms%h_layer, &
                                            ms%u_face_x_layer, &
                                            nxt, nyt, nz, i_e, j_s, j_n, grid%nghost, &
                                            get_clamped_tracer(bc%east%clamped_tracer, it))
               end if
            end if

            ! ---- South ghost fill ----
            if (is_open_ish(bc_s)) then
               if (allocated(bc%tres_south)) then
                  call fill_ghost_south_res( &
                     ms%tracers(it)%hTr, ms%h_layer, bc%tres_south, &
                     nxt, nyt, nz, j_s, i_w, i_e, grid%nghost, it, size(bc%tres_south, 3))
               else
                  call fill_ghost_south_sign(ms%tracers(it)%hTr, ms%h_layer, &
                                             ms%v_face_y_layer, &
                                             nxt, nyt, nz, j_s, i_w, i_e, grid%nghost, &
                                             get_clamped_tracer(bc%south%clamped_tracer, it))
               end if
            end if

            ! ---- North ghost fill ----
            if (is_open_ish(bc_n)) then
               if (allocated(bc%tres_north)) then
                  call fill_ghost_north_res( &
                     ms%tracers(it)%hTr, ms%h_layer, bc%tres_north, &
                     nxt, nyt, nz, j_n, i_w, i_e, grid%nghost, it, size(bc%tres_north, 3))
               else
                  call fill_ghost_north_sign(ms%tracers(it)%hTr, ms%h_layer, &
                                             ms%v_face_y_layer, &
                                             nxt, nyt, nz, j_n, i_w, i_e, grid%nghost, &
                                             get_clamped_tracer(bc%north%clamped_tracer, it))
               end if
            end if

         end do
      end if
   end subroutine ocean_obc_fill_ghosts

   subroutine ocean_obc_refill_ghost_ssh(grid, bc, ms, bt_H_ref)
      !! Re-establish a zero-gradient free surface in the open-edge GHOST
      !! columns, called at the END of the outer step (after the ALE remap)
      !! so the diagnostic manager sees a consistent halo.
      !!
      !! Why this is needed: the ghost cells are updated by the slow
      !! continuity over the full array (`i = 1..nx_total`) using the
      !! zeroed array-edge fluxes, and the conservative ALE remap
      !! (`remap_h_ref = Σh − bt_eta`) PRESERVES that spurious ghost
      !! transport divergence.  Nothing re-imposes the open-boundary
      !! zero-gradient invariant between the remap and the diagnostic
      !! read, so `SSH = Σ_k h_layer − b` in the ghosts drifts to tens of
      !! metres (worst where the boundary bathymetry is steep and at the
      !! corners where two open edges meet), while the physical interior
      !! stays healthy.
      !!
      !! Fix: overwrite each open-edge ghost column with the nearest
      !! interior column SCALED so the ghost column total equals
      !! `bt_H_ref(ghost) + η_interior`, i.e. the free-surface anomaly is
      !! flat across the open boundary.  Since `bt_H_ref == b`, this makes
      !! `SSH_ghost = η_interior` for ANY bathymetry — including formula
      !! topographies whose ghost `b` differs from the first interior cell.
      !! For file bathymetry (ghost `b` constant-extrapolated) the scale
      !! collapses to 1, i.e. a plain zero-gradient re-copy that simply
      !! discards the accumulated ghost drift.
      !!
      !! Physics-neutral: the next outer step's stage-1 `ocean_obc_fill_ghosts`
      !! re-fills the ghost h_layer (thickness copy) BEFORE any physics
      !! kernel reads it, so this pass only affects what the diagnostics
      !! (and other halo consumers) observe — the prognostic trajectory is
      !! bit-identical.  Gated on open-ish edges ⇒ WALL/PERIODIC ⇒ no-op.
      type(hgrid_t), intent(in) :: grid
      type(ocean_bc_state_t), intent(in) :: bc
      type(multilayer_state_t), intent(inout) :: ms
      real(wp), intent(in) :: bt_H_ref(grid%nx_total, grid%ny_total)
         !! Mode-split reference column depth (== seeded bathymetry `b`).

      integer :: bc_w, bc_e, bc_s, bc_n
      logical :: any_open
      integer :: it, nz, nxt, nyt, nx, ny
      integer :: i_w, i_e, j_s, j_n

      bc_w = bc%west%bc_type
      bc_e = bc%east%bc_type
      bc_s = bc%south%bc_type
      bc_n = bc%north%bc_type

      any_open = is_open_ish(bc_w) .or. is_open_ish(bc_e) .or. &
                 is_open_ish(bc_s) .or. is_open_ish(bc_n)
      if (.not. any_open) return
      if (.not. allocated(ms%h_layer)) return

      nz = ms%nz_ml
      nxt = grid%nx_total
      nyt = grid%ny_total
      nx = grid%nx_phys
      ny = grid%ny_phys
      i_w = grid%nghost + 1          ! west physical wall-face / first interior cell
      i_e = grid%nghost + nx         ! last interior cell (east)
      j_s = grid%nghost + 1          ! south physical wall-face / first interior cell
      j_n = grid%nghost + ny         ! last interior cell (north)

      ! ---- h_layer ghosts: zero-gradient free surface (scaled to ghost b) ----
      call refill_h_ghost_scaled(ms%h_layer, bt_H_ref, nxt, nyt, nz, &
                                 i_w, i_e, j_s, j_n, grid%nghost, &
                                 bc_w, bc_e, bc_s, bc_n)

      ! ---- Tracer hTr ghosts: zero-gradient concentration onto the new h ----
      if (allocated(ms%tracers)) then
         do it = 1, size(ms%tracers)
            if (.not. allocated(ms%tracers(it)%hTr)) cycle
            call fill_tracer_ghosts_zerograd(ms%tracers(it)%hTr, ms%h_layer, &
                                             nxt, nyt, nz, i_w, i_e, j_s, j_n, &
                                             grid%nghost, bc_w, bc_e, bc_s, bc_n)
         end do
      end if
   end subroutine ocean_obc_refill_ghost_ssh

   pure subroutine refill_h_ghost_scaled(h_layer, bt_H_ref, nx_total, ny_total, nz, &
                                         i_w, i_e, j_s, j_n, nghost, &
                                         bc_w, bc_e, bc_s, bc_n)
      !! Overwrite open-edge h_layer ghost columns with the nearest interior
      !! column scaled so `Σ_k h_ghost = bt_H_ref(ghost) + η_interior`
      !! (η_interior = `Σ_k h_int − bt_H_ref(int)`).  x-pass (west/east) over
      !! the full j extent (covers corner rows); y-pass (south/north) over the
      !! full i extent, reading the already x-filled corner columns — same
      !! corner-coverage order as `fill_h_layer_ghosts`.
      integer, intent(in) :: nx_total, ny_total, nz, nghost
      real(wp), intent(inout) :: h_layer(nx_total, ny_total, nz)
      real(wp), intent(in)    :: bt_H_ref(nx_total, ny_total)
      integer, intent(in) :: i_w, i_e, j_s, j_n
      integer, intent(in) :: bc_w, bc_e, bc_s, bc_n

      integer :: i, j, k, g
      real(wp) :: sh_int, eta_ref, sc

      ! West ghosts
      if (is_open_ish(bc_w)) then
         do concurrent(j=1:ny_total, g=1:nghost) local(k, sh_int, eta_ref, sc)
            sh_int = 0.0_wp
            do k = 1, nz
               sh_int = sh_int + h_layer(i_w, j, k)
            end do
            eta_ref = sh_int - bt_H_ref(i_w, j)
            sc = (bt_H_ref(g, j) + eta_ref)/max(sh_int, H_VANISHED)
            do k = 1, nz
               h_layer(g, j, k) = h_layer(i_w, j, k)*sc
            end do
         end do
      end if
      ! East ghosts
      if (is_open_ish(bc_e)) then
         do concurrent(j=1:ny_total, g=1:nghost) local(k, sh_int, eta_ref, sc)
            sh_int = 0.0_wp
            do k = 1, nz
               sh_int = sh_int + h_layer(i_e, j, k)
            end do
            eta_ref = sh_int - bt_H_ref(i_e, j)
            sc = (bt_H_ref(nx_total - g + 1, j) + eta_ref)/max(sh_int, H_VANISHED)
            do k = 1, nz
               h_layer(nx_total - g + 1, j, k) = h_layer(i_e, j, k)*sc
            end do
         end do
      end if
      ! South ghosts (read the x-filled corner columns)
      if (is_open_ish(bc_s)) then
         do concurrent(i=1:nx_total, g=1:nghost) local(k, sh_int, eta_ref, sc)
            sh_int = 0.0_wp
            do k = 1, nz
               sh_int = sh_int + h_layer(i, j_s, k)
            end do
            eta_ref = sh_int - bt_H_ref(i, j_s)
            sc = (bt_H_ref(i, g) + eta_ref)/max(sh_int, H_VANISHED)
            do k = 1, nz
               h_layer(i, g, k) = h_layer(i, j_s, k)*sc
            end do
         end do
      end if
      ! North ghosts (read the x-filled corner columns)
      if (is_open_ish(bc_n)) then
         do concurrent(i=1:nx_total, g=1:nghost) local(k, sh_int, eta_ref, sc)
            sh_int = 0.0_wp
            do k = 1, nz
               sh_int = sh_int + h_layer(i, j_n, k)
            end do
            eta_ref = sh_int - bt_H_ref(i, j_n)
            sc = (bt_H_ref(i, ny_total - g + 1) + eta_ref)/max(sh_int, H_VANISHED)
            do k = 1, nz
               h_layer(i, ny_total - g + 1, k) = h_layer(i, j_n, k)*sc
            end do
         end do
      end if
   end subroutine refill_h_ghost_scaled

   pure function get_clamped_tracer(arr, it) result(val)
      !! Safe accessor: return clamped_tracer(it) or 0 if unallocated/out of range.
      real(wp), allocatable, intent(in) :: arr(:)
      integer, intent(in) :: it
      real(wp) :: val
      val = 0.0_wp
      if (allocated(arr)) then
         if (it <= size(arr)) val = arr(it)
      end if
   end function get_clamped_tracer

   pure subroutine fill_tracer_ghosts_zerograd(hTr, h_layer, nx_total, ny_total, nz, &
                                               i_w, i_e, j_s, j_n, nghost, &
                                               bc_w, bc_e, bc_s, bc_n)
      !! Zero-gradient CONCENTRATION fill of a tracer's hTr ghosts at open-ish
      !! edges over the FULL cross-extent (ghost×ghost corners included):
      !! hTr_ghost = (hTr_int / h_int) * h_ghost.  Runs BEFORE the upwind-aware
      !! per-edge fill, so corners keep this zero-gradient value (no corner T/S
      !! blow-up under inflow).  Same gating/corner coverage as fill_h_layer_ghosts.
      integer, intent(in) :: nx_total, ny_total, nz, nghost
      real(wp), intent(inout) :: hTr(nx_total, ny_total, nz)
      real(wp), intent(in)    :: h_layer(nx_total, ny_total, nz)
      integer, intent(in) :: i_w, i_e, j_s, j_n
      integer, intent(in) :: bc_w, bc_e, bc_s, bc_n
      integer :: i, j, k, g

      ! West / East over the full j extent (covers corner rows).
      if (is_open_ish(bc_w)) then
         do concurrent(k=1:nz, j=1:ny_total, g=1:nghost)
            hTr(g, j, k) = (hTr(i_w, j, k)/max(h_layer(i_w, j, k), H_VANISHED))*h_layer(g, j, k)
         end do
      end if
      if (is_open_ish(bc_e)) then
         do concurrent(k=1:nz, j=1:ny_total, g=1:nghost)
            hTr(nx_total - g + 1, j, k) = &
               (hTr(i_e, j, k)/max(h_layer(i_e, j, k), H_VANISHED))*h_layer(nx_total - g + 1, j, k)
         end do
      end if
      ! South / North over the full i extent (reads the x-filled corner cols).
      if (is_open_ish(bc_s)) then
         do concurrent(k=1:nz, j=1:nghost, i=1:nx_total)
            hTr(i, j, k) = (hTr(i, j_s, k)/max(h_layer(i, j_s, k), H_VANISHED))*h_layer(i, j, k)
         end do
      end if
      if (is_open_ish(bc_n)) then
         do concurrent(k=1:nz, g=1:nghost, i=1:nx_total)
            hTr(i, ny_total - g + 1, k) = &
               (hTr(i, j_n, k)/max(h_layer(i, j_n, k), H_VANISHED))*h_layer(i, ny_total - g + 1, k)
         end do
      end if
   end subroutine fill_tracer_ghosts_zerograd

   pure subroutine fill_h_layer_ghosts(h_layer, nx_total, ny_total, nz, &
                                       i_w, i_e, j_s, j_n, nghost, &
                                       bc_w, bc_e, bc_s, bc_n)
      !! Zero-gradient fill of h_layer ghosts at open-ish edges.
      integer, intent(in) :: nx_total, ny_total, nz, nghost
      real(wp), intent(inout) :: h_layer(nx_total, ny_total, nz)
      integer, intent(in) :: i_w, i_e, j_s, j_n
      integer, intent(in) :: bc_w, bc_e, bc_s, bc_n

      integer :: i, j, k, g

      ! West ghosts: copy first interior column (i_w) into ghost columns
      if (is_open_ish(bc_w)) then
         do concurrent(k=1:nz, j=1:ny_total, g=1:nghost)
            h_layer(g, j, k) = h_layer(i_w, j, k)
         end do
      end if
      ! East ghosts
      if (is_open_ish(bc_e)) then
         do concurrent(k=1:nz, j=1:ny_total, g=1:nghost)
            h_layer(nx_total - g + 1, j, k) = h_layer(i_e, j, k)
         end do
      end if
      ! South ghosts
      if (is_open_ish(bc_s)) then
         do concurrent(k=1:nz, j=1:nghost, i=1:nx_total)
            h_layer(i, j, k) = h_layer(i, j_s, k)
         end do
      end if
      ! North ghosts
      if (is_open_ish(bc_n)) then
         do concurrent(k=1:nz, g=1:nghost, i=1:nx_total)
            h_layer(i, ny_total - g + 1, k) = h_layer(i, j_n, k)
         end do
      end if
   end subroutine fill_h_layer_ghosts

   pure subroutine fill_tracer_ghosts_zonal(hTr, h_layer, u_layer, &
                                            nx_total, ny_total, nz, &
                                            i_w, i_e, j_s, j_n, nghost, &
                                            bc_w, bc_e, &
                                            clamped_tr_w, clamped_tr_e)
      !! Upwind-aware tracer ghost fill for west and east open edges.
      !! West inflow  : u(i_w,j,k) > 0  → ghost = clamped_tr * h_ghost
      !! West outflow : u(i_w,j,k) <= 0 → ghost = interior (zero-gradient)
      !! East inflow  : u(i_e,j,k) < 0
      !! East outflow : u(i_e,j,k) >= 0
      integer, intent(in) :: nx_total, ny_total, nz, nghost
      real(wp), intent(inout) :: hTr(nx_total, ny_total, nz)
      real(wp), intent(in)    :: h_layer(nx_total, ny_total, nz)
      real(wp), intent(in)    :: u_layer(nx_total + 1, ny_total, nz)
      integer, intent(in) :: i_w, i_e, j_s, j_n
      integer, intent(in) :: bc_w, bc_e
      real(wp), intent(in) :: clamped_tr_w, clamped_tr_e

      integer :: j, k, g
      real(wp) :: hTr_bc, hTr_int

      ! West ghosts
      if (is_open_ish(bc_w)) then
         do concurrent(k=1:nz, j=j_s:j_n, g=1:nghost) local(hTr_bc, hTr_int)
            ! Use wall-face velocity at u_layer(i_w, j, k).
            ! West inflow = u > 0 (eastward into domain).
            if (u_layer(i_w, j, k) > 0.0_wp) then
               ! Inflow: set ghost hTr = boundary concentration * ghost h
               hTr_bc = clamped_tr_w*h_layer(g, j, k)
               hTr(g, j, k) = hTr_bc
            else
               ! Outflow: zero-gradient — copy interior
               hTr_int = hTr(i_w, j, k)
               hTr(g, j, k) = hTr_int
            end if
         end do
      end if

      ! East ghosts
      if (is_open_ish(bc_e)) then
         do concurrent(k=1:nz, j=j_s:j_n, g=1:nghost) local(hTr_bc, hTr_int)
            ! East inflow = u < 0 (westward into domain from east).
            if (u_layer(i_e + 1, j, k) < 0.0_wp) then
               hTr_bc = clamped_tr_e*h_layer(nx_total - g + 1, j, k)
               hTr(nx_total - g + 1, j, k) = hTr_bc
            else
               hTr_int = hTr(i_e, j, k)
               hTr(nx_total - g + 1, j, k) = hTr_int
            end if
         end do
      end if
   end subroutine fill_tracer_ghosts_zonal

   pure subroutine fill_tracer_ghosts_meridional(hTr, h_layer, v_layer, &
                                                 nx_total, ny_total, nz, &
                                                 i_w, i_e, j_s, j_n, nghost, &
                                                 bc_s, bc_n, &
                                                 clamped_tr_s, clamped_tr_n)
      !! Upwind-aware tracer ghost fill for south and north open edges.
      integer, intent(in) :: nx_total, ny_total, nz, nghost
      real(wp), intent(inout) :: hTr(nx_total, ny_total, nz)
      real(wp), intent(in)    :: h_layer(nx_total, ny_total, nz)
      real(wp), intent(in)    :: v_layer(nx_total, ny_total + 1, nz)
      integer, intent(in) :: i_w, i_e, j_s, j_n
      integer, intent(in) :: bc_s, bc_n
      real(wp), intent(in) :: clamped_tr_s, clamped_tr_n

      integer :: i, k, g
      real(wp) :: hTr_bc, hTr_int

      ! South ghosts
      if (is_open_ish(bc_s)) then
         do concurrent(k=1:nz, i=i_w:i_e, g=1:nghost) local(hTr_bc, hTr_int)
            ! South inflow = v > 0 (northward into domain).
            if (v_layer(i, j_s, k) > 0.0_wp) then
               hTr_bc = clamped_tr_s*h_layer(i, g, k)
               hTr(i, g, k) = hTr_bc
            else
               hTr_int = hTr(i, j_s, k)
               hTr(i, g, k) = hTr_int
            end if
         end do
      end if

      ! North ghosts
      if (is_open_ish(bc_n)) then
         do concurrent(k=1:nz, i=i_w:i_e, g=1:nghost) local(hTr_bc, hTr_int)
            ! North inflow = v < 0 (southward into domain from north).
            if (v_layer(i, j_n + 1, k) < 0.0_wp) then
               hTr_bc = clamped_tr_n*h_layer(i, ny_total - g + 1, k)
               hTr(i, ny_total - g + 1, k) = hTr_bc
            else
               hTr_int = hTr(i, j_n, k)
               hTr(i, ny_total - g + 1, k) = hTr_int
            end if
         end do
      end if
   end subroutine fill_tracer_ghosts_meridional

   ! -----------------------------------------------------------------
   ! Per-edge helpers: reservoir path (unconditional concentration)
   ! -----------------------------------------------------------------

   pure subroutine fill_ghost_west_res(hTr, h_layer, tres_w, &
                                       nx_total, ny_total, nz, &
                                       i_w, j0, j1, nghost, it, n_tr)
      !! Reservoir-based ghost fill, west edge.
      !! hTr_ghost = tres_w(j, k, it) * h_ghost — unconditional.
      integer, intent(in) :: nx_total, ny_total, nz, nghost, it, n_tr
      real(wp), intent(inout) :: hTr(nx_total, ny_total, nz)
      real(wp), intent(in)    :: h_layer(nx_total, ny_total, nz)
      real(wp), intent(in)    :: tres_w(ny_total, nz, n_tr)
      integer, intent(in) :: i_w, j0, j1
      integer :: j, k, g
      real(wp) :: hTr_bc
      do concurrent(k=1:nz, j=j0:j1, g=1:nghost) local(hTr_bc)
         hTr_bc = tres_w(j, k, it)*h_layer(g, j, k)
         hTr(g, j, k) = hTr_bc
      end do
   end subroutine fill_ghost_west_res

   pure subroutine fill_ghost_east_res(hTr, h_layer, tres_e, &
                                       nx_total, ny_total, nz, &
                                       i_e, j0, j1, nghost, it, n_tr)
      !! Reservoir-based ghost fill, east edge.
      integer, intent(in) :: nx_total, ny_total, nz, nghost, it, n_tr
      real(wp), intent(inout) :: hTr(nx_total, ny_total, nz)
      real(wp), intent(in)    :: h_layer(nx_total, ny_total, nz)
      real(wp), intent(in)    :: tres_e(ny_total, nz, n_tr)
      integer, intent(in) :: i_e, j0, j1
      integer :: j, k, g
      real(wp) :: hTr_bc
      do concurrent(k=1:nz, j=j0:j1, g=1:nghost) local(hTr_bc)
         hTr_bc = tres_e(j, k, it)*h_layer(nx_total - g + 1, j, k)
         hTr(nx_total - g + 1, j, k) = hTr_bc
      end do
   end subroutine fill_ghost_east_res

   pure subroutine fill_ghost_south_res(hTr, h_layer, tres_s, &
                                        nx_total, ny_total, nz, &
                                        j_s, i0, i1, nghost, it, n_tr)
      !! Reservoir-based ghost fill, south edge.
      integer, intent(in) :: nx_total, ny_total, nz, nghost, it, n_tr
      real(wp), intent(inout) :: hTr(nx_total, ny_total, nz)
      real(wp), intent(in)    :: h_layer(nx_total, ny_total, nz)
      real(wp), intent(in)    :: tres_s(nx_total, nz, n_tr)
      integer, intent(in) :: j_s, i0, i1
      integer :: i, k, g
      real(wp) :: hTr_bc
      do concurrent(k=1:nz, i=i0:i1, g=1:nghost) local(hTr_bc)
         hTr_bc = tres_s(i, k, it)*h_layer(i, g, k)
         hTr(i, g, k) = hTr_bc
      end do
   end subroutine fill_ghost_south_res

   pure subroutine fill_ghost_north_res(hTr, h_layer, tres_n, &
                                        nx_total, ny_total, nz, &
                                        j_n, i0, i1, nghost, it, n_tr)
      !! Reservoir-based ghost fill, north edge.
      integer, intent(in) :: nx_total, ny_total, nz, nghost, it, n_tr
      real(wp), intent(inout) :: hTr(nx_total, ny_total, nz)
      real(wp), intent(in)    :: h_layer(nx_total, ny_total, nz)
      real(wp), intent(in)    :: tres_n(nx_total, nz, n_tr)
      integer, intent(in) :: j_n, i0, i1
      integer :: i, k, g
      real(wp) :: hTr_bc
      do concurrent(k=1:nz, i=i0:i1, g=1:nghost) local(hTr_bc)
         hTr_bc = tres_n(i, k, it)*h_layer(i, ny_total - g + 1, k)
         hTr(i, ny_total - g + 1, k) = hTr_bc
      end do
   end subroutine fill_ghost_north_res

   ! Per-edge helpers: sign-switch path (default fallback)

   pure subroutine fill_ghost_west_sign(hTr, h_layer, u_layer, &
                                        nx_total, ny_total, nz, &
                                        i_w, j0, j1, nghost, clamped_tr_w)
      !! Sign-switch ghost fill for the west edge.
      integer, intent(in) :: nx_total, ny_total, nz, nghost
      real(wp), intent(inout) :: hTr(nx_total, ny_total, nz)
      real(wp), intent(in)    :: h_layer(nx_total, ny_total, nz)
      real(wp), intent(in)    :: u_layer(nx_total + 1, ny_total, nz)
      integer, intent(in) :: i_w, j0, j1
      real(wp), intent(in) :: clamped_tr_w
      integer :: j, k, g
      real(wp) :: hTr_bc, hTr_int
      do concurrent(k=1:nz, j=j0:j1, g=1:nghost) local(hTr_bc, hTr_int)
         if (u_layer(i_w, j, k) > 0.0_wp) then
            hTr_bc = clamped_tr_w*h_layer(g, j, k)
            hTr(g, j, k) = hTr_bc
         else
            hTr_int = hTr(i_w, j, k)
            hTr(g, j, k) = hTr_int
         end if
      end do
   end subroutine fill_ghost_west_sign

   pure subroutine fill_ghost_east_sign(hTr, h_layer, u_layer, &
                                        nx_total, ny_total, nz, &
                                        i_e, j0, j1, nghost, clamped_tr_e)
      !! Sign-switch ghost fill for the east edge.
      integer, intent(in) :: nx_total, ny_total, nz, nghost
      real(wp), intent(inout) :: hTr(nx_total, ny_total, nz)
      real(wp), intent(in)    :: h_layer(nx_total, ny_total, nz)
      real(wp), intent(in)    :: u_layer(nx_total + 1, ny_total, nz)
      integer, intent(in) :: i_e, j0, j1
      real(wp), intent(in) :: clamped_tr_e
      integer :: j, k, g
      real(wp) :: hTr_bc, hTr_int
      do concurrent(k=1:nz, j=j0:j1, g=1:nghost) local(hTr_bc, hTr_int)
         if (u_layer(i_e + 1, j, k) < 0.0_wp) then
            hTr_bc = clamped_tr_e*h_layer(nx_total - g + 1, j, k)
            hTr(nx_total - g + 1, j, k) = hTr_bc
         else
            hTr_int = hTr(i_e, j, k)
            hTr(nx_total - g + 1, j, k) = hTr_int
         end if
      end do
   end subroutine fill_ghost_east_sign

   pure subroutine fill_ghost_south_sign(hTr, h_layer, v_layer, &
                                         nx_total, ny_total, nz, &
                                         j_s, i0, i1, nghost, clamped_tr_s)
      !! Sign-switch ghost fill for the south edge.
      integer, intent(in) :: nx_total, ny_total, nz, nghost
      real(wp), intent(inout) :: hTr(nx_total, ny_total, nz)
      real(wp), intent(in)    :: h_layer(nx_total, ny_total, nz)
      real(wp), intent(in)    :: v_layer(nx_total, ny_total + 1, nz)
      integer, intent(in) :: j_s, i0, i1
      real(wp), intent(in) :: clamped_tr_s
      integer :: i, k, g
      real(wp) :: hTr_bc, hTr_int
      do concurrent(k=1:nz, i=i0:i1, g=1:nghost) local(hTr_bc, hTr_int)
         if (v_layer(i, j_s, k) > 0.0_wp) then
            hTr_bc = clamped_tr_s*h_layer(i, g, k)
            hTr(i, g, k) = hTr_bc
         else
            hTr_int = hTr(i, j_s, k)
            hTr(i, g, k) = hTr_int
         end if
      end do
   end subroutine fill_ghost_south_sign

   pure subroutine fill_ghost_north_sign(hTr, h_layer, v_layer, &
                                         nx_total, ny_total, nz, &
                                         j_n, i0, i1, nghost, clamped_tr_n)
      !! Sign-switch ghost fill for the north edge.
      integer, intent(in) :: nx_total, ny_total, nz, nghost
      real(wp), intent(inout) :: hTr(nx_total, ny_total, nz)
      real(wp), intent(in)    :: h_layer(nx_total, ny_total, nz)
      real(wp), intent(in)    :: v_layer(nx_total, ny_total + 1, nz)
      integer, intent(in) :: j_n, i0, i1
      real(wp), intent(in) :: clamped_tr_n
      integer :: i, k, g
      real(wp) :: hTr_bc, hTr_int
      do concurrent(k=1:nz, i=i0:i1, g=1:nghost) local(hTr_bc, hTr_int)
         if (v_layer(i, j_n + 1, k) < 0.0_wp) then
            hTr_bc = clamped_tr_n*h_layer(i, ny_total - g + 1, k)
            hTr(i, ny_total - g + 1, k) = hTr_bc
         else
            hTr_int = hTr(i, j_n, k)
            hTr(i, ny_total - g + 1, k) = hTr_int
         end if
      end do
   end subroutine fill_ghost_north_sign

   ! Open-edge tracer reservoir update

   subroutine ocean_obc_update_reservoirs(grid, bc, ms, dt)
      !! Evolve per-edge reservoir concentrations `tres` one timestep.
      !! Called after `continuity_tracer_step_split` while
      !! `ms%mass_flux_{x,y}_layer` still hold the stage's wall-face fluxes.
      !! Per open-ish edge (with allocated `tres_*`), per (j|i, k, tracer):
      !!   u_n = sign_edge * mass_flux(wall) / max(h_int, h_min)  (outward normal)
      !! Implicit backward-Euler (Marchesiello et al. 2001):
      !!   c_out = max(0,u_n)*dt/L_out, c_in = max(0,-u_n)*dt/L_in  (0 if L==0)
      !!   tres = (tres + c_out*T_int + c_in*T_data)/(1 + c_out + c_in)
      !! Degenerate L_out==0 & u_n>0 ⇒ tres = T_int (L_in==0 & u_n<0 ⇒ T_data).
      !! T_int = hTr(int)/max(h(int),h_min); T_data = clamped_tracer(it).
      !! No-op when both length scales are zero (default path).
      type(hgrid_t), intent(in) :: grid
      type(ocean_bc_state_t), intent(inout) :: bc
      type(multilayer_state_t), intent(in) :: ms
      real(wp), intent(in) :: dt

      logical :: do_res
      integer :: nxt, nyt, nz, nx, ny
      integer :: i_w, i_e, j_s, j_n
      integer :: it, n_tr

      do_res = (bc%res_lscale_out > 0.0_wp .or. bc%res_lscale_in > 0.0_wp)
      if (.not. do_res) return

      nxt = grid%nx_total
      nyt = grid%ny_total
      nz = ms%nz_ml
      nx = grid%nx_phys
      ny = grid%ny_phys
      i_w = grid%nghost + 1
      i_e = grid%nghost + nx
      j_s = grid%nghost + 1
      j_n = grid%nghost + ny

      n_tr = bc%n_tracers
      if (n_tr <= 0) return
      if (.not. allocated(ms%tracers)) return

      ! Per-tracer outer shim: pass full explicit-shape arrays to kernels.
      do it = 1, n_tr
         if (.not. allocated(ms%tracers(it)%hTr)) cycle

         ! ---- West ----
         ! MPI-seam neutralisation (O0): dispatch here is on allocated(tres_*),
         ! not bc_type, so the remap trick cannot be used — explicit has_* guard.
         if (allocated(bc%tres_west) .and. bc%has_west) then
            call update_reservoir_zonal_west( &
               ms%tracers(it)%hTr, ms%h_layer, ms%mass_flux_x_layer, &
               bc%tres_west, &
               get_clamped_tracer(bc%west%clamped_tracer, it), &
               nxt, nyt, nz, i_w, j_s, j_n, &
               bc%res_lscale_out, bc%res_lscale_in, dt, it, n_tr)
         end if

         ! ---- East ----
         if (allocated(bc%tres_east) .and. bc%has_east) then
            call update_reservoir_zonal_east( &
               ms%tracers(it)%hTr, ms%h_layer, ms%mass_flux_x_layer, &
               bc%tres_east, &
               get_clamped_tracer(bc%east%clamped_tracer, it), &
               nxt, nyt, nz, i_e, j_s, j_n, &
               bc%res_lscale_out, bc%res_lscale_in, dt, it, n_tr)
         end if

         ! ---- South ----
         if (allocated(bc%tres_south) .and. bc%has_south) then
            call update_reservoir_meridional_south( &
               ms%tracers(it)%hTr, ms%h_layer, ms%mass_flux_y_layer, &
               bc%tres_south, &
               get_clamped_tracer(bc%south%clamped_tracer, it), &
               nxt, nyt, nz, j_s, i_w, i_e, &
               bc%res_lscale_out, bc%res_lscale_in, dt, it, n_tr)
         end if

         ! ---- North ----
         if (allocated(bc%tres_north) .and. bc%has_north) then
            call update_reservoir_meridional_north( &
               ms%tracers(it)%hTr, ms%h_layer, ms%mass_flux_y_layer, &
               bc%tres_north, &
               get_clamped_tracer(bc%north%clamped_tracer, it), &
               nxt, nyt, nz, j_n, i_w, i_e, &
               bc%res_lscale_out, bc%res_lscale_in, dt, it, n_tr)
         end if

      end do
   end subroutine ocean_obc_update_reservoirs

   pure subroutine update_reservoir_zonal_west(hTr, h_layer, mass_flux_x, &
                                               tres_w, T_data, &
                                               nx_total, ny_total, nz, &
                                               i_w, j0, j1, &
                                               L_out, L_in, dt, it, n_tr)
      !! Update reservoir for the west open edge.  Wall-face/interior cell = i_w.
      !! Outward normal = -x: u_n = -mass_flux_x(i_w,j,k)/max(h,hmin).
      integer, intent(in) :: nx_total, ny_total, nz, n_tr, it
      integer, intent(in) :: i_w, j0, j1
      real(wp), intent(in) :: hTr(nx_total, ny_total, nz)
      real(wp), intent(in) :: h_layer(nx_total, ny_total, nz)
      real(wp), intent(in) :: mass_flux_x(nx_total + 1, ny_total, nz)
      real(wp), intent(inout) :: tres_w(ny_total, nz, n_tr)
      real(wp), intent(in) :: T_data
      real(wp), intent(in) :: L_out, L_in, dt

      integer :: j, k
      real(wp) :: h_int, T_int, u_n, c_out, c_in, tres_old, denom

      do concurrent(j=j0:j1, k=1:nz) local(h_int, T_int, u_n, c_out, c_in, tres_old, denom)
         h_int = max(h_layer(i_w, j, k), RES_H_MIN)
         T_int = hTr(i_w, j, k)/h_int
         ! mass_flux_x(i_w,...) > 0 = eastward into domain (inflow at west).
         u_n = -mass_flux_x(i_w, j, k)/h_int
         tres_old = tres_w(j, k, it)
         ! Degenerate limits: L_out==0 & u_n>0 ⇒ T_int; L_in==0 & u_n<0 ⇒ T_data
         ! instantly.  Both L zero ⇒ implicit form is a no-op (tres unchanged).
         if (L_out == 0.0_wp .and. u_n > 0.0_wp) then
            tres_w(j, k, it) = T_int
         else if (L_in == 0.0_wp .and. u_n < 0.0_wp) then
            tres_w(j, k, it) = T_data
         else
            c_out = merge(max(0.0_wp, u_n)*dt/L_out, 0.0_wp, L_out > 0.0_wp)
            c_in = merge(max(0.0_wp, -u_n)*dt/L_in, 0.0_wp, L_in > 0.0_wp)
            denom = 1.0_wp + c_out + c_in
            tres_w(j, k, it) = (tres_old + c_out*T_int + c_in*T_data)/denom
         end if
      end do
   end subroutine update_reservoir_zonal_west

   pure subroutine update_reservoir_zonal_east(hTr, h_layer, mass_flux_x, &
                                               tres_e, T_data, &
                                               nx_total, ny_total, nz, &
                                               i_e, j0, j1, &
                                               L_out, L_in, dt, it, n_tr)
      !! Update reservoir for the east open edge.
      !! Wall-face index for mass flux = i_e+1  (east face of the last physical cell i_e).
      !! Interior cell = i_e.
      !! Outward normal at east = +x: u_n = +mass_flux_x(i_e+1,j,k)/max(h,hmin).
      integer, intent(in) :: nx_total, ny_total, nz, n_tr, it
      integer, intent(in) :: i_e, j0, j1
      real(wp), intent(in) :: hTr(nx_total, ny_total, nz)
      real(wp), intent(in) :: h_layer(nx_total, ny_total, nz)
      real(wp), intent(in) :: mass_flux_x(nx_total + 1, ny_total, nz)
      real(wp), intent(inout) :: tres_e(ny_total, nz, n_tr)
      real(wp), intent(in) :: T_data
      real(wp), intent(in) :: L_out, L_in, dt

      integer :: j, k
      real(wp) :: h_int, T_int, u_n, c_out, c_in, tres_old, denom

      do concurrent(j=j0:j1, k=1:nz) local(h_int, T_int, u_n, c_out, c_in, tres_old, denom)
         h_int = max(h_layer(i_e, j, k), RES_H_MIN)
         T_int = hTr(i_e, j, k)/h_int
         ! Outward normal at east = +x.
         u_n = mass_flux_x(i_e + 1, j, k)/h_int
         tres_old = tres_e(j, k, it)
         if (L_out == 0.0_wp .and. u_n > 0.0_wp) then
            tres_e(j, k, it) = T_int
         else if (L_in == 0.0_wp .and. u_n < 0.0_wp) then
            tres_e(j, k, it) = T_data
         else
            c_out = merge(max(0.0_wp, u_n)*dt/L_out, 0.0_wp, L_out > 0.0_wp)
            c_in = merge(max(0.0_wp, -u_n)*dt/L_in, 0.0_wp, L_in > 0.0_wp)
            denom = 1.0_wp + c_out + c_in
            tres_e(j, k, it) = (tres_old + c_out*T_int + c_in*T_data)/denom
         end if
      end do
   end subroutine update_reservoir_zonal_east

   pure subroutine update_reservoir_meridional_south(hTr, h_layer, mass_flux_y, &
                                                     tres_s, T_data, &
                                                     nx_total, ny_total, nz, &
                                                     j_s, i0, i1, &
                                                     L_out, L_in, dt, it, n_tr)
      !! Update reservoir for the south open edge.
      !! Wall-face index = j_s; interior cell = j_s.
      !! Outward normal = -y: u_n = -mass_flux_y(i,j_s,k)/max(h,hmin).
      integer, intent(in) :: nx_total, ny_total, nz, n_tr, it
      integer, intent(in) :: j_s, i0, i1
      real(wp), intent(in) :: hTr(nx_total, ny_total, nz)
      real(wp), intent(in) :: h_layer(nx_total, ny_total, nz)
      real(wp), intent(in) :: mass_flux_y(nx_total, ny_total + 1, nz)
      real(wp), intent(inout) :: tres_s(nx_total, nz, n_tr)
      real(wp), intent(in) :: T_data
      real(wp), intent(in) :: L_out, L_in, dt

      integer :: i, k
      real(wp) :: h_int, T_int, u_n, c_out, c_in, tres_old, denom

      do concurrent(i=i0:i1, k=1:nz) local(h_int, T_int, u_n, c_out, c_in, tres_old, denom)
         h_int = max(h_layer(i, j_s, k), RES_H_MIN)
         T_int = hTr(i, j_s, k)/h_int
         u_n = -mass_flux_y(i, j_s, k)/h_int
         tres_old = tres_s(i, k, it)
         if (L_out == 0.0_wp .and. u_n > 0.0_wp) then
            tres_s(i, k, it) = T_int
         else if (L_in == 0.0_wp .and. u_n < 0.0_wp) then
            tres_s(i, k, it) = T_data
         else
            c_out = merge(max(0.0_wp, u_n)*dt/L_out, 0.0_wp, L_out > 0.0_wp)
            c_in = merge(max(0.0_wp, -u_n)*dt/L_in, 0.0_wp, L_in > 0.0_wp)
            denom = 1.0_wp + c_out + c_in
            tres_s(i, k, it) = (tres_old + c_out*T_int + c_in*T_data)/denom
         end if
      end do
   end subroutine update_reservoir_meridional_south

   pure subroutine update_reservoir_meridional_north(hTr, h_layer, mass_flux_y, &
                                                     tres_n, T_data, &
                                                     nx_total, ny_total, nz, &
                                                     j_n, i0, i1, &
                                                     L_out, L_in, dt, it, n_tr)
      !! Update reservoir for the north open edge.
      !! Wall-face index = j_n+1; interior cell = j_n.
      !! Outward normal = +y: u_n = +mass_flux_y(i,j_n+1,k)/max(h,hmin).
      integer, intent(in) :: nx_total, ny_total, nz, n_tr, it
      integer, intent(in) :: j_n, i0, i1
      real(wp), intent(in) :: hTr(nx_total, ny_total, nz)
      real(wp), intent(in) :: h_layer(nx_total, ny_total, nz)
      real(wp), intent(in) :: mass_flux_y(nx_total, ny_total + 1, nz)
      real(wp), intent(inout) :: tres_n(nx_total, nz, n_tr)
      real(wp), intent(in) :: T_data
      real(wp), intent(in) :: L_out, L_in, dt

      integer :: i, k
      real(wp) :: h_int, T_int, u_n, c_out, c_in, tres_old, denom

      do concurrent(i=i0:i1, k=1:nz) local(h_int, T_int, u_n, c_out, c_in, tres_old, denom)
         h_int = max(h_layer(i, j_n, k), RES_H_MIN)
         T_int = hTr(i, j_n, k)/h_int
         u_n = mass_flux_y(i, j_n + 1, k)/h_int
         tres_old = tres_n(i, k, it)
         if (L_out == 0.0_wp .and. u_n > 0.0_wp) then
            tres_n(i, k, it) = T_int
         else if (L_in == 0.0_wp .and. u_n < 0.0_wp) then
            tres_n(i, k, it) = T_data
         else
            c_out = merge(max(0.0_wp, u_n)*dt/L_out, 0.0_wp, L_out > 0.0_wp)
            c_in = merge(max(0.0_wp, -u_n)*dt/L_in, 0.0_wp, L_in > 0.0_wp)
            denom = 1.0_wp + c_out + c_in
            tres_n(i, k, it) = (tres_old + c_out*T_int + c_in*T_data)/denom
         end if
      end do
   end subroutine update_reservoir_meridional_north

end module rdb_ocean_obc_baroclinic
