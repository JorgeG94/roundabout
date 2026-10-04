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
   use rdb_ocean_fold_apply, only: ocean_fold_wrap_stress
   use rdb_ocean_surface_stress, only: ocean_surface_stress_t, &
                                       ocean_surface_stress_set_derived
   use rdb_ice_state, only: ocean_sea_ice_t
   use rdb_constants, only: wp
   implicit none
   private

   public :: ocean_halo_exchange_ml_state
   public :: ocean_seam_refresh_surface_stress
   public :: ocean_halo_exchange_ice_state
   public :: ocean_halo_exchange_ice_fluxes
   public :: ocean_halo_exchange_ice_transport

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
      ! because those are NOT vectors).  px = 1: the local kernels; px > 1:
      ! one owner-routed exchange group (`ocean_fold_wrap_stress`).
      if (bc%north_fold) call ocean_fold_wrap_stress(grid, bc, ss%tau_x, ss%tau_y, device_resident)

      call ocean_surface_stress_set_derived(grid, ss)
   end subroutine ocean_seam_refresh_surface_stress

   subroutine ocean_halo_exchange_ice_state(ice, device_resident)
      !! Make the sea-ice CATEGORY state valid in every ghost cell (X1 of
      !! the sea-ice MPI plan): `part_size`, `m_ice`, `m_snow`,
      !! `enth_ice`, `sal_ice`, `enth_snow`, one two-pass centre exchange
      !! each, all categories (and ice layers) in one message per
      !! direction.
      !!
      !! **Why.**  Nothing inside the ice step refreshes these ghosts — the
      !! column thermodynamics, ITD and (PR 4b) transport compress write
      !! PHYSICAL cells only — yet three consumers read them one cell into
      !! the halo: the EVP's category gather (`mis`/`mice`/`ci` over the
      !! full array, then strength, face mass and corner ratios), and the
      !! stress coupler's face concentration `a_u = (ci(i-1)+ci(i))/2` at
      !! the west/south-most owned face.  On one rank with a periodic axis
      !! the primitives' local wrap closes the seam the same way (the
      !! pre-existing "D7" stale-ghost note in `rdb_ice_evp`).
      !!
      !! **When.**  At the end of every thermo block (the category state
      !! changes only there) and once at cold-start configure, host-side,
      !! BEFORE `enter_data` and NEVER on a warm restart (the checkpoint
      !! carries the writer's ghosts; re-deriving them resumed a different
      !! state, `8e1931f20`).  Single-rank non-periodic: a no-op.  Requires
      !! `ocean_halo_init`; the caller skips it otherwise.
      !!
      !! The rank-4 `enth_ice`/`sal_ice`/`enth_snow` and the `0:ncat`
      !! `part_size` are contiguous and go out as one flat `nz` each,
      !! through an explicit-shape seam (`ice_halo_centre_flat`), never
      !! the aggregate `ice` (CLAUDE.md: component arrays only).
      type(ocean_sea_ice_t), intent(inout) :: ice
         !! Live sea-ice slot (`ice%is_init`); a no-op otherwise.
      logical, intent(in), optional :: device_resident
         !! Forwarded to the halo primitives; `.false.` for the host-side
         !! configure-time call.

      integer :: nxt, nyt

      if (.not. ice%is_init) return
      nxt = ice%nx_total
      nyt = ice%ny_total

      call profiler_start("ice_comms_state")
      call oh_count_suppress_on()
      call ice_halo_centre_flat(ice%part_size, nxt, nyt, ice%ncat + 1, device_resident)
      call ice_halo_centre_flat(ice%m_ice, nxt, nyt, ice%ncat, device_resident)
      call ice_halo_centre_flat(ice%m_snow, nxt, nyt, ice%ncat, device_resident)
      call ice_halo_centre_flat(ice%enth_ice, nxt, nyt, ice%ncat*ice%nk_ice, device_resident)
      call ice_halo_centre_flat(ice%sal_ice, nxt, nyt, ice%ncat*ice%nk_ice, device_resident)
      call ice_halo_centre_flat(ice%enth_snow, nxt, nyt, ice%ncat, device_resident)
      call oh_count_suppress_off()
      call profiler_stop("ice_comms_state")
   end subroutine ocean_halo_exchange_ice_state

   subroutine ocean_halo_exchange_ice_fluxes(ice, device_resident)
      !! Seam ghosts of the three per-cell ice->ocean flux diagnostics the
      !! couplers hand to the ocean — `salt_flux_diag`, `heat_flux_diag`,
      !! `sw_thru_diag` — one two-pass centre exchange each.
      !!
      !! **Why.**  The column driver, frazil uptake and snowfall share write
      !! them on PHYSICAL cells only, but the brine / heat / shortwave
      !! couplers copy them over the FULL array into `Q_salt` / `Q_heat` /
      !! `q_sw` (or their components), and the ocean's surface-flux
      !! application reads a seam ghost before its next exchange.  A stale
      !! ghost there is the neighbour's flux replaced by this tile's old
      !! one (measured: the first decomposed run diverged from the serial
      !! one in the outer step after the first thermo block).  Call after
      !! the last contributor and before the couplers.  Single-rank
      !! non-periodic: a no-op; requires `ocean_halo_init`.
      type(ocean_sea_ice_t), intent(inout) :: ice
         !! Live sea-ice slot (`ice%is_init`); a no-op otherwise.
      logical, intent(in), optional :: device_resident
         !! Forwarded to the halo primitives.

      if (.not. ice%is_init) return
      call profiler_start("ice_comms_fluxes")
      call oh_count_suppress_on()
      call ocean_halo_centre(ice%salt_flux_diag, device_resident)
      call ocean_halo_centre(ice%heat_flux_diag, device_resident)
      call ocean_halo_centre(ice%sw_thru_diag, device_resident)
      call oh_count_suppress_off()
      call profiler_stop("ice_comms_fluxes")
   end subroutine ocean_halo_exchange_ice_fluxes

   subroutine ocean_halo_exchange_ice_transport(ice, device_resident)
      !! X4 of the sea-ice MPI plan: the seam ghosts every advective
      !! substep of `ice_transport_step` reads — the cell-averaged
      !! category masses `mca_ice`/`mca_snow` (the PPM donors, 5-point
      !! stencil) and the riding intensive tracers `m_ice`, `enth_ice`,
      !! `sal_ice`, `enth_snow` (the PCM donors).  `mca_*` ghosts are
      !! zeroed by the IST->CAS conversion and the ride/mass updates leave
      !! the ghost band one substep old, so this runs at the top of EVERY
      !! substep.  On one rank with a periodic axis the primitives wrap.
      type(ocean_sea_ice_t), intent(inout) :: ice
         !! Live sea-ice slot (`ice%is_init`); a no-op otherwise.
      logical, intent(in), optional :: device_resident
         !! Forwarded to the halo primitives.

      integer :: nxt, nyt

      if (.not. ice%is_init) return
      nxt = ice%nx_total
      nyt = ice%ny_total
      call profiler_start("ice_comms_transport")
      call oh_count_suppress_on()
      call ice_halo_centre_flat(ice%mca_ice, nxt, nyt, ice%ncat, device_resident)
      call ice_halo_centre_flat(ice%mca_snow, nxt, nyt, ice%ncat, device_resident)
      call ice_halo_centre_flat(ice%m_ice, nxt, nyt, ice%ncat, device_resident)
      call ice_halo_centre_flat(ice%enth_ice, nxt, nyt, ice%ncat*ice%nk_ice, device_resident)
      call ice_halo_centre_flat(ice%sal_ice, nxt, nyt, ice%ncat*ice%nk_ice, device_resident)
      call ice_halo_centre_flat(ice%enth_snow, nxt, nyt, ice%ncat, device_resident)
      call oh_count_suppress_off()
      call profiler_stop("ice_comms_transport")
   end subroutine ocean_halo_exchange_ice_transport

   subroutine ice_halo_centre_flat(fld, nxt, nyt, nz, device_resident)
      !! Explicit-shape seam: a contiguous ice array of any rank (`0:ncat`
      !! third bound, rank-4 category x layer) is handed in by sequence
      !! association and exchanged as one `(nxt, nyt, nz)` centre field.
      !! The generic `ocean_halo_centre` resolves on the DUMMY's rank, so
      !! the rank-4 actuals cannot call it directly.
      integer, intent(in) :: nxt, nyt, nz
      real(wp), intent(inout) :: fld(nxt, nyt, nz)
      logical, intent(in), optional :: device_resident

      call ocean_halo_centre(fld, nz, device_resident)
   end subroutine ice_halo_centre_flat

end module rdb_ocean_halo_state
