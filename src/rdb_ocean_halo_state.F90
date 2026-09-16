!! State-level shim routing ML prognostic ghost exchange through the O1
!! comm primitives.
!!
!! This module provides a single call site for the outer-step ghost fill of
!! the ocean multilayer C-grid prognostic fields.  The routing follows D0
!! (MPI-agnostic solver code): on a single-rank non-periodic build the calls
!! resolve to no-ops; on a single-rank periodic build they resolve to local
!! wrap copies (handled inside `rdb_ocean_halo`); on a multi-rank build they
!! resolve to messages.  Corner ghosts are valid after the call because the
!! O1 primitives perform two-pass (E/W then N/S) exchanges internally (D2).
!!
!! The per-tracer loop is OUTSIDE any `do concurrent` region (outer-shim
!! pattern: array-of-derived-types cannot be dereferenced on-device).
module rdb_ocean_halo_state
   use rdb_ocean_halo, only: ocean_halo_centre, &
                             ocean_halo_face_x, &
                             ocean_halo_face_y, &
                             ocean_halo_is_decomposed_x, &
                             ocean_halo_is_decomposed_y
   use rdb_ocean_halo_counters, only: oh_count_ml_state, &
                                      oh_count_suppress_on, oh_count_suppress_off
   use rdb_profiler, only: profiler_start, profiler_stop
   use rdb_multilayer_state, only: multilayer_state_t
   use rdb_grid, only: hgrid_t
   use rdb_ocean_boundary_types, only: ocean_bc_state_t
   use rdb_ocean_periodic, only: ocean_periodic_wrap_face_x_2d, &
                                 ocean_periodic_wrap_face_y_2d
   use rdb_ocean_fold, only: fold_north_u_face, fold_north_v_face
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t, &
                                       ocean_surface_stress_set_derived
   implicit none
   private

   public :: ocean_halo_exchange_ml_state
   public :: ocean_seam_refresh_surface_stress

contains

   subroutine ocean_halo_exchange_ml_state(ms, device_resident)
      !! Exchange ghost cells for the four multilayer prognostic field kinds
      !! via the O1 halo primitives:
      !!
      !!   * `h_layer`          — cell-centred layer thickness (centre_3d)
      !!   * `u_face_x_layer`   — east-face layer velocity (face_x_3d)
      !!   * `v_face_y_layer`   — north-face layer velocity (face_y_3d)
      !!   * `tracers(it)%hTr`  — per-tracer thickness-weighted scalar
      !!                          (centre_3d, outer-shim loop)
      !!
      !! Unconditional (D0): single-rank + non-periodic ⇒ no-op;
      !! single-rank + periodic ⇒ local wrap; multi-rank ⇒ messages.
      !! Corner ghosts valid on return (D2 two-pass inside the primitives).
      type(multilayer_state_t), intent(inout) :: ms
         !! Multilayer C-grid state whose ghost bands are to be filled.
      logical, intent(in), optional :: device_resident
         !! Forwarded to every primitive call.  Pass .false. for init-time
         !! host-side exchanges that occur before ocean_state_enter_data.
         !! Default (.true.) is the normal device-resident path.

      integer :: it

      call profiler_start("ocean_comms_ml")
      call oh_count_ml_state()
      call oh_count_suppress_on()
      call ocean_halo_centre(ms%h_layer, ms%nz_ml, device_resident)
      call ocean_halo_face_x(ms%u_face_x_layer, ms%nz_ml, device_resident)
      call ocean_halo_face_y(ms%v_face_y_layer, ms%nz_ml, device_resident)

      ! Per-tracer loop outside DC — outer-shim pattern: array-of-DT
      ! cannot be dereferenced inside a do concurrent on-device.
      if (allocated(ms%tracers)) then
         do it = 1, size(ms%tracers)
            if (.not. allocated(ms%tracers(it)%hTr)) cycle
            call ocean_halo_centre(ms%tracers(it)%hTr, ms%nz_ml, device_resident)
         end do
      end if
      call oh_count_suppress_off()
      call profiler_stop("ocean_comms_ml")

   end subroutine ocean_halo_exchange_ml_state

   subroutine ocean_seam_refresh_surface_stress(ss, grid, bc, device_resident)
      !! Make the surface-stress pair valid in every ghost cell, then
      !! re-derive `stress_mag` from it.
      !!
      !! **Why this exists.**  `tau_x`/`tau_y` are C-grid face fields that
      !! several kernels read ONE CELL BEYOND the cell they write:
      !!
      !!   * `ocean_surfstress_derived_impl` averages `tau_x(i)`+`tau_x(i+1)`
      !!     into the cell-centred `stress_mag`, which feeds KPP/EPBL `u_*`.
      !!   * `mle_face_ustar_x/y` (Fox-Kemper / Bodner) take a 4-point
      !!     corner average reaching `tau_y(i-1, ·)` / `tau_x(·, j-1)`.
      !!
      !! At an MPI seam those reads land in ghost cells that belong to the
      !! neighbour rank, so they MUST come from an exchange.  Nothing may
      !! extrapolate them: a zero-gradient / edge-copy fill silently
      !! substitutes this rank's edge value for the neighbour's real data,
      !! which is decomposition-dependent and therefore invisible to any
      !! single-rank test.  (This mirrors the convention in MOM6, which
      !! halo-exchanges the stress pair and never extrapolates forcing.)
      !!
      !! **Order is load-bearing** and matches the prognostic-state path in
      !! `ocean_dyn_step_split`: exchange, THEN periodic wrap on any axis
      !! the exchange did not own, THEN the north fold.  The periodic
      !! kernels are skipped on a decomposed axis because the halo already
      !! filled those ghosts — running both would overwrite correct
      !! neighbour data with a local wrap (see `rdb_ocean_periodic`).
      !!
      !! `stress_mag` needs no wrap or fold of its own: it is recomputed
      !! LAST, over the full array, from a `tau` pair whose ghosts are by
      !! then already valid, so its ghosts come out right for free.
      !!
      !! Safe to call before `ocean_state_enter_data` with
      !! `device_resident = .false.` (the configure-time seed path).
      type(ocean_surface_stress_t), intent(inout) :: ss
         !! Stress slot whose `tau_x`/`tau_y` ghosts are to be filled and
         !! whose `stress_mag` is then refreshed.
      type(hgrid_t), intent(in) :: grid
      type(ocean_bc_state_t), intent(in) :: bc
         !! Supplies `periodic_x`/`periodic_y`/`north_fold`.
      logical, intent(in), optional :: device_resident
         !! Forwarded to the halo primitives; `.false.` for host-side
         !! configure-time calls.

      integer :: nxt, nyt, nxp, nyp, ng
      logical :: wrap_x, wrap_y

      if (.not. allocated(ss%tau_x) .or. .not. allocated(ss%tau_y)) return

      nxt = grid%nx_total
      nyt = grid%ny_total
      nxp = grid%nx_phys
      nyp = grid%ny_phys
      ng = grid%nghost

      call profiler_start("ocean_comms_stress")
      call oh_count_suppress_on()
      call ocean_halo_face_x(ss%tau_x, device_resident)
      call ocean_halo_face_y(ss%tau_y, device_resident)
      call oh_count_suppress_off()
      call profiler_stop("ocean_comms_stress")

      ! Local periodic wrap only on an axis the halo did NOT own.
      wrap_x = bc%periodic_x .and. (.not. ocean_halo_is_decomposed_x())
      wrap_y = bc%periodic_y .and. (.not. ocean_halo_is_decomposed_y())
      if (wrap_x .or. wrap_y) then
         call ocean_periodic_wrap_face_x_2d(ss%tau_x, nxt + 1, nyt, &
                                            nxp, nyp, ng, wrap_x, wrap_y)
         call ocean_periodic_wrap_face_y_2d(ss%tau_y, nxt, nyt + 1, &
                                            nxp, nyp, ng, wrap_x, wrap_y)
      end if

      ! Tripolar north fold.  `tau_x`/`tau_y` are TRUE VECTOR components,
      ! so the sign-flipping u/v-face variants are the correct ones (the
      ! scalar-copy duplicates in `rdb_ocean_metrics` exist precisely
      ! because those are NOT vectors).  Fold requires px == 1, which the
      ! driver already fences (D6).
      if (bc%north_fold) then
         call fold_north_u_face(ss%tau_x, nxt + 1, nyt, nxp, nyp, ng)
         call fold_north_v_face(ss%tau_y, nxt, nyt + 1, nxp, nyp, ng)
      end if

      call ocean_surface_stress_set_derived(grid, ss)
   end subroutine ocean_seam_refresh_surface_stress

end module rdb_ocean_halo_state
