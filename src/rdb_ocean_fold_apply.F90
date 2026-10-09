!! State-level tripolar north-fold seam application for the ocean dyn-core.
module rdb_ocean_fold_apply
   !! Orchestration over the pure seam operators in `rdb_ocean_fold`:
   !! dereferences the state slots (outer-shim per-tracer loop) and calls
   !! the explicit-shape fold kernels.
   !!
   !! Ordering contract: periodic-x is wrapped FIRST, fold SECOND — every
   !! fold routine is called AFTER the matching periodic wrap so it reads
   !! the already cyclically-wrapped corner columns.
   !!
   !! Stagger map:
   !!   * h_layer, η, tracer hTr  → centre fold (copy, no sign flip)
   !!   * u_face_x_layer, bt_ubt  → u-face fold (negate — true vector)
   !!   * v_face_y_layer, bt_vbt  → v-face fold (negate + on-row
   !!     antisymmetric projection of the fold line, storage row
   !!     nghost+ny_phys+1 — see the `rdb_ocean_fold` header)
   !!
   !! Every routine no-ops when `bc%north_fold` is .false. ⇒ non-tripolar
   !! runs stay bit-identical.  `bc%north_fold` is RANK-LOCAL: true only on
   !! the rank that owns the physical north edge (`has_north`), so on a
   !! north-south split the other ranks leave their north ghosts to the MPI
   !! exchange that precedes every call here (exchange → periodic wrap →
   !! fold).  See the `rdb_ocean_fold` header for the decomposition limits.
   !!
   !! ## px dispatch
   !!
   !! `px = 1` (`ocean_fold_is_distributed()` false): the local kernels,
   !! textually unchanged.  `px > 1`: every routine becomes ONE collective
   !! group of the owner-routed exchange (`rdb_ocean_fold_exchange`) over
   !! the north rank row — all its ranks have `north_fold` set and reach the
   !! same call — with the same fields, signs and fold-line projection.  The
   !! optional `device_resident` only matters on that path (host-side setup
   !! calls pass `.false.`, as for the halo primitives).
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_ocean_boundary_types, only: ocean_bc_state_t
   use rdb_ocean_fold, only: fold_north_centre, fold_north_u_face, &
                             fold_north_v_face
   use rdb_ocean_fold_exchange, only: ocean_fold_is_distributed, ocean_fold_begin, &
                                      ocean_fold_pack, ocean_fold_exchange, &
                                      ocean_fold_unpack, ocean_fold_end, &
                                      FOLD_STAG_T, FOLD_STAG_U, FOLD_STAG_V
   implicit none
   private

   public :: ocean_fold_wrap_state
   public :: ocean_fold_wrap_centre_3d_state
   public :: ocean_fold_wrap_eta_2d
   public :: ocean_fold_wrap_time_means
   public :: ocean_fold_wrap_stress
   public :: ocean_fold_wrap_visc_rem
   public :: ocean_fold_wrap_centre_flat

contains

   pure integer function n_tracers(ms) result(n)
      !! Registered tracers carrying an allocated `hTr`.
      type(multilayer_state_t), intent(in) :: ms
      integer :: it
      n = 0
      if (.not. allocated(ms%tracers)) return
      do it = 1, size(ms%tracers)
         if (allocated(ms%tracers(it)%hTr)) n = n + 1
      end do
   end function n_tracers

   subroutine ocean_fold_wrap_state(grid, bc, ms, device_resident)
      !! Fold the north seam of h_layer, u/v layer faces, and every
      !! registered tracer.  Call AFTER `ocean_periodic_wrap_state`.
      !! No-op when `bc%north_fold` is .false.
      type(hgrid_t), intent(in) :: grid
      type(ocean_bc_state_t), intent(in) :: bc
      type(multilayer_state_t), intent(inout) :: ms
      logical, intent(in), optional :: device_resident
         !! px > 1 only: `.false.` for host-side (pre-`enter_data`) calls.

      integer :: it
      integer :: nx, ny, nz, nx_phys, ny_phys, nghost

      if (.not. bc%north_fold) return

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml
      nx_phys = grid%nx_phys
      ny_phys = grid%ny_phys
      nghost = grid%nghost

      if (ocean_fold_is_distributed()) then
         ! One group: h (T), u (u, −), v (v, −, + fold row), every hTr (T).
         call ocean_fold_begin((nghost + 1)*nz*(3 + n_tracers(ms)))
         call ocean_fold_pack(ms%h_layer, nx, ny, nz, FOLD_STAG_T, device_resident)
         call ocean_fold_pack(ms%u_face_x_layer, nx + 1, ny, nz, FOLD_STAG_U, device_resident)
         call ocean_fold_pack(ms%v_face_y_layer, nx, ny + 1, nz, FOLD_STAG_V, device_resident)
         if (allocated(ms%tracers)) then
            do it = 1, size(ms%tracers)
               if (.not. allocated(ms%tracers(it)%hTr)) cycle
               call ocean_fold_pack(ms%tracers(it)%hTr, nx, ny, nz, FOLD_STAG_T, device_resident)
            end do
         end if
         call ocean_fold_exchange(device_resident)
         call ocean_fold_unpack(ms%h_layer, nx, ny, nz, FOLD_STAG_T, .false., device_resident)
         call ocean_fold_unpack(ms%u_face_x_layer, nx + 1, ny, nz, FOLD_STAG_U, .true., &
                                device_resident)
         call ocean_fold_unpack(ms%v_face_y_layer, nx, ny + 1, nz, FOLD_STAG_V, .true., &
                                device_resident)
         if (allocated(ms%tracers)) then
            do it = 1, size(ms%tracers)
               if (.not. allocated(ms%tracers(it)%hTr)) cycle
               call ocean_fold_unpack(ms%tracers(it)%hTr, nx, ny, nz, FOLD_STAG_T, .false., &
                                      device_resident)
            end do
         end if
         call ocean_fold_end()
         return
      end if

      ! Centre (T): h_layer.
      call fold_north_centre(ms%h_layer, nx, ny, nz, &
                             nx_phys, ny_phys, nghost)
      ! u-face (Cu): negate.
      call fold_north_u_face(ms%u_face_x_layer, nx + 1, ny, nz, &
                             nx_phys, ny_phys, nghost)
      ! v-face (Cv): negate + on-row antisymmetric projection.
      call fold_north_v_face(ms%v_face_y_layer, nx, ny + 1, nz, &
                             nx_phys, ny_phys, nghost)

      ! Per-tracer loop OUTSIDE the DC kernels (outer-shim for array-of-DTs).
      if (allocated(ms%tracers)) then
         do it = 1, size(ms%tracers)
            if (.not. allocated(ms%tracers(it)%hTr)) cycle
            call fold_north_centre(ms%tracers(it)%hTr, nx, ny, nz, &
                                   nx_phys, ny_phys, nghost)
         end do
      end if
   end subroutine ocean_fold_wrap_state

   subroutine ocean_fold_wrap_centre_3d_state(grid, bc, ms, device_resident)
      !! Fold ONLY h_layer + tracers (centre fields) — the continuity
      !! mid-split site, which re-wraps the centre fields between the
      !! zonal and meridional Lie-split halves.  No-op when not folding.
      type(hgrid_t), intent(in) :: grid
      type(ocean_bc_state_t), intent(in) :: bc
      type(multilayer_state_t), intent(inout) :: ms
      logical, intent(in), optional :: device_resident
         !! px > 1 only: `.false.` for host-side calls.

      integer :: it
      integer :: nx, ny, nz, nx_phys, ny_phys, nghost

      if (.not. bc%north_fold) return

      nx = grid%nx_total
      ny = grid%ny_total
      nz = ms%nz_ml
      nx_phys = grid%nx_phys
      ny_phys = grid%ny_phys
      nghost = grid%nghost

      if (ocean_fold_is_distributed()) then
         ! One group: h + every hTr (all T, copy).
         call ocean_fold_begin(nghost*nz*(1 + n_tracers(ms)))
         call ocean_fold_pack(ms%h_layer, nx, ny, nz, FOLD_STAG_T, device_resident)
         if (allocated(ms%tracers)) then
            do it = 1, size(ms%tracers)
               if (.not. allocated(ms%tracers(it)%hTr)) cycle
               call ocean_fold_pack(ms%tracers(it)%hTr, nx, ny, nz, FOLD_STAG_T, device_resident)
            end do
         end if
         call ocean_fold_exchange(device_resident)
         call ocean_fold_unpack(ms%h_layer, nx, ny, nz, FOLD_STAG_T, .false., device_resident)
         if (allocated(ms%tracers)) then
            do it = 1, size(ms%tracers)
               if (.not. allocated(ms%tracers(it)%hTr)) cycle
               call ocean_fold_unpack(ms%tracers(it)%hTr, nx, ny, nz, FOLD_STAG_T, .false., &
                                      device_resident)
            end do
         end if
         call ocean_fold_end()
         return
      end if

      call fold_north_centre(ms%h_layer, nx, ny, nz, &
                             nx_phys, ny_phys, nghost)
      if (allocated(ms%tracers)) then
         do it = 1, size(ms%tracers)
            if (.not. allocated(ms%tracers(it)%hTr)) cycle
            call fold_north_centre(ms%tracers(it)%hTr, nx, ny, nz, &
                                   nx_phys, ny_phys, nghost)
         end do
      end if
   end subroutine ocean_fold_wrap_centre_3d_state

   subroutine ocean_fold_wrap_eta_2d(grid, bc, eta, device_resident)
      !! Fold a 2D cell-centred η field (driver-level SSH wrap site).
      !! No-op when not folding.
      type(hgrid_t), intent(in) :: grid
      type(ocean_bc_state_t), intent(in) :: bc
      real(wp), intent(inout) :: eta(:, :)
      logical, intent(in), optional :: device_resident
         !! px > 1 only: `.false.` for host-side calls.

      if (.not. bc%north_fold) return

      if (ocean_fold_is_distributed()) then
         call ocean_fold_begin(grid%nghost)
         call ocean_fold_pack(eta, grid%nx_total, grid%ny_total, FOLD_STAG_T, device_resident)
         call ocean_fold_exchange(device_resident)
         call ocean_fold_unpack(eta, grid%nx_total, grid%ny_total, FOLD_STAG_T, .false., &
                                device_resident)
         call ocean_fold_end()
         return
      end if

      call fold_north_centre(eta, grid%nx_total, grid%ny_total, &
                             grid%nx_phys, grid%ny_phys, grid%nghost)
   end subroutine ocean_fold_wrap_eta_2d

   subroutine ocean_fold_wrap_time_means(grid, bc, ms)
      !! Fold the `pred_corr` step time-means u_av (u, −), v_av (v, −, +
      !! fold-row projection) and h_av (T) — the stage-entry site that
      !! mirrors the prognostic fold for the Coriolis / viscosity inputs.
      !! Call after their periodic wrap.  No-op when not folding or when
      !! the means are not allocated (ssp_rk2).  Device-only: unlike the
      !! other dispatchers it takes no `device_resident` flag, because its
      !! production caller (`run_stage_split`) always runs on the mapped state.
      type(hgrid_t), intent(in) :: grid
      type(ocean_bc_state_t), intent(in) :: bc
      type(multilayer_state_t), intent(inout) :: ms

      integer :: nxu, nyu, nxv, nyv, nxh, nyh, nz

      if (.not. bc%north_fold) return
      if (.not. allocated(ms%u_av_layer)) return

      nxu = size(ms%u_av_layer, 1)
      nyu = size(ms%u_av_layer, 2)
      nxv = size(ms%v_av_layer, 1)
      nyv = size(ms%v_av_layer, 2)
      nxh = size(ms%h_av_layer, 1)
      nyh = size(ms%h_av_layer, 2)
      nz = size(ms%u_av_layer, 3)

      if (ocean_fold_is_distributed()) then
         call ocean_fold_begin(3*(grid%nghost + 1)*nz)
         call ocean_fold_pack(ms%u_av_layer, nxu, nyu, nz, FOLD_STAG_U)
         call ocean_fold_pack(ms%v_av_layer, nxv, nyv, nz, FOLD_STAG_V)
         call ocean_fold_pack(ms%h_av_layer, nxh, nyh, nz, FOLD_STAG_T)
         call ocean_fold_exchange()
         call ocean_fold_unpack(ms%u_av_layer, nxu, nyu, nz, FOLD_STAG_U, .true.)
         call ocean_fold_unpack(ms%v_av_layer, nxv, nyv, nz, FOLD_STAG_V, .true.)
         call ocean_fold_unpack(ms%h_av_layer, nxh, nyh, nz, FOLD_STAG_T, .false.)
         call ocean_fold_end()
         return
      end if

      call fold_north_u_face(ms%u_av_layer, nxu, nyu, nz, &
                             grid%nx_phys, grid%ny_phys, grid%nghost)
      call fold_north_v_face(ms%v_av_layer, nxv, nyv, nz, &
                             grid%nx_phys, grid%ny_phys, grid%nghost)
      call fold_north_centre(ms%h_av_layer, nxh, nyh, nz, &
                             grid%nx_phys, grid%ny_phys, grid%nghost)
   end subroutine ocean_fold_wrap_time_means

   subroutine ocean_fold_wrap_stress(grid, bc, tau_x, tau_y, device_resident)
      !! Fold the surface-stress pair: `tau_x` (u, −) and `tau_y` (v, − +
      !! fold-row projection) — true vector components.  Call after the
      !! pair's exchange + periodic wrap.  No-op when not folding.
      type(hgrid_t), intent(in) :: grid
      type(ocean_bc_state_t), intent(in) :: bc
      real(wp), intent(inout) :: tau_x(grid%nx_total + 1, grid%ny_total)
         !! x-face stress (nx_total+1, ny_total).
      real(wp), intent(inout) :: tau_y(grid%nx_total, grid%ny_total + 1)
         !! y-face stress (nx_total, ny_total+1).
      logical, intent(in), optional :: device_resident
         !! px > 1 only: `.false.` for host-side calls.

      integer :: nxt, nyt

      if (.not. bc%north_fold) return
      nxt = grid%nx_total
      nyt = grid%ny_total

      if (ocean_fold_is_distributed()) then
         call ocean_fold_begin(2*(grid%nghost + 1))
         call ocean_fold_pack(tau_x, nxt + 1, nyt, FOLD_STAG_U, device_resident)
         call ocean_fold_pack(tau_y, nxt, nyt + 1, FOLD_STAG_V, device_resident)
         call ocean_fold_exchange(device_resident)
         call ocean_fold_unpack(tau_x, nxt + 1, nyt, FOLD_STAG_U, .true., device_resident)
         call ocean_fold_unpack(tau_y, nxt, nyt + 1, FOLD_STAG_V, .true., device_resident)
         call ocean_fold_end()
         return
      end if

      call fold_north_u_face(tau_x, nxt + 1, nyt, grid%nx_phys, grid%ny_phys, grid%nghost)
      call fold_north_v_face(tau_y, nxt, nyt + 1, grid%nx_phys, grid%ny_phys, grid%nghost)
   end subroutine ocean_fold_wrap_stress

   subroutine ocean_fold_wrap_visc_rem(grid, bc, visc_rem_u, visc_rem_v, device_resident)
      !! Fold the viscous-remnant pair: `visc_rem_u` (u-face) and
      !! `visc_rem_v` (v-face) — PR-1's `bt_work%visc_rem_u/v` seam.
      !! UNLIKE `ocean_fold_wrap_stress` (its vector twin, tau_x/tau_y),
      !! `visc_rem` is a POSITIVE SCALAR (the fraction of a barotropic
      !! acceleration a layer still feels after one implicit-friction
      !! step, MOM6 `vertvisc_remnant`), not a flux/velocity component,
      !! so both face kernels are called
      !! with `negate=.false.`: the 180-degree fold still swaps which side
      !! of the seam the ghost value comes from, but the value itself does
      !! not change sign, and the v-face fold-line duplicate DOF is forced
      !! EQUAL (not opposite) across the seam.  Call after the pair's halo
      !! exchange + periodic wrap (MOM6's `pass_visc_rem` group pass, run
      !! after every one of the three `vertvisc_remnant` calls).  No-op
      !! when not folding.
      type(hgrid_t), intent(in) :: grid
      type(ocean_bc_state_t), intent(in) :: bc
      real(wp), intent(inout) :: visc_rem_u(:, :, :)
         !! u-face per-layer remnant, shape (nx_total+1, ny_total, nz).
      real(wp), intent(inout) :: visc_rem_v(:, :, :)
         !! v-face per-layer remnant, shape (nx_total, ny_total+1, nz).
      logical, intent(in), optional :: device_resident
         !! px > 1 only: `.false.` for host-side calls.

      integer :: nxu, nyu, nxv, nyv, nz

      if (.not. bc%north_fold) return

      nxu = size(visc_rem_u, 1)
      nyu = size(visc_rem_u, 2)
      nxv = size(visc_rem_v, 1)
      nyv = size(visc_rem_v, 2)
      nz = size(visc_rem_u, 3)

      if (ocean_fold_is_distributed()) then
         call ocean_fold_begin(2*(grid%nghost + 1)*nz)
         call ocean_fold_pack(visc_rem_u, nxu, nyu, nz, FOLD_STAG_U, device_resident)
         call ocean_fold_pack(visc_rem_v, nxv, nyv, nz, FOLD_STAG_V, device_resident)
         call ocean_fold_exchange(device_resident)
         call ocean_fold_unpack(visc_rem_u, nxu, nyu, nz, FOLD_STAG_U, .false., device_resident)
         call ocean_fold_unpack(visc_rem_v, nxv, nyv, nz, FOLD_STAG_V, .false., device_resident)
         call ocean_fold_end()
         return
      end if

      call fold_north_u_face(visc_rem_u, nxu, nyu, nz, &
                             grid%nx_phys, grid%ny_phys, grid%nghost, negate=.false.)
      call fold_north_v_face(visc_rem_v, nxv, nyv, nz, &
                             grid%nx_phys, grid%ny_phys, grid%nghost, negate=.false.)
   end subroutine ocean_fold_wrap_visc_rem

   subroutine ocean_fold_wrap_centre_flat(grid, bc, fld, nxt, nyt, nz, device_resident)
      !! Fold the north seam of an arbitrary flat cell-centred (T-stagger)
      !! SCALAR field — copy, no sign flip, exactly `fold_north_centre`'s
      !! contract.  Exists so a caller holding a sequence-associated flat
      !! view of a higher-rank array (the sea-ice category state's
      !! `ice_halo_centre_flat` pattern — `part_size`/`m_ice`/`m_snow`/
      !! `enth_ice`/`sal_ice`/`enth_snow`/`mca_ice`/`mca_snow`, and the
      !! per-cell flux diagnostics `salt_flux_diag`/`heat_flux_diag`/
      !! `sw_thru_diag`) can fold it without going through
      !! `multilayer_state_t`.  Every one of these is a per-category mass,
      !! enthalpy, salinity or fractional area — a scalar, not a vector
      !! component — so `negate` is never offered here.  Call AFTER the
      !! field's halo exchange (which also performs the single-rank
      !! periodic wrap).  No-op when `bc%north_fold` is `.false.`.
      type(hgrid_t), intent(in) :: grid
      type(ocean_bc_state_t), intent(in) :: bc
      integer, intent(in) :: nxt, nyt, nz
      real(wp), intent(inout) :: fld(nxt, nyt, nz)
      logical, intent(in), optional :: device_resident
         !! px > 1 only: `.false.` for host-side (pre-`enter_data`) calls.

      if (.not. bc%north_fold) return

      if (ocean_fold_is_distributed()) then
         call ocean_fold_begin(grid%nghost*nz)
         call ocean_fold_pack(fld, nxt, nyt, nz, FOLD_STAG_T, device_resident)
         call ocean_fold_exchange(device_resident)
         call ocean_fold_unpack(fld, nxt, nyt, nz, FOLD_STAG_T, .false., device_resident)
         call ocean_fold_end()
         return
      end if

      call fold_north_centre(fld, nxt, nyt, nz, grid%nx_phys, grid%ny_phys, grid%nghost)
   end subroutine ocean_fold_wrap_centre_flat

end module rdb_ocean_fold_apply
