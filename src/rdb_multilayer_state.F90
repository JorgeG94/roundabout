!! 3D hydrostatic multilayer state on an Arakawa C-grid.
module rdb_multilayer_state
   !! Per-layer prognostic state for the ocean dynamical core.
   !! Layout mirrors the 2D `barotropic_state_t`, lifted to (i, j, k):
   !!
   !!   * Scalars (h_layer, tracers) at cell centres, shape
   !!     (nx_total, ny_total, nz_ml).
   !!   * x-velocity / x-momentum on east faces, shape
   !!     (nx_total+1, ny_total, nz_ml).
   !!   * y-velocity / y-momentum on north faces, shape
   !!     (nx_total, ny_total+1, nz_ml).
   !!   * Per-face per-layer mass fluxes — the continuity-PPM
   !!     output that the tracer kernels consume.
   !!
   !! Lives at `ocean_state%multilayer` when `sim_type='ocean'`.
   !! Caller must set `this%nz_ml` before calling `init` — typically
   !! threaded through from `cfg%nz_layers` in `state_init_from_config`.
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_tracer, only: tracer_t, TRACER_BUDGET_NONE, TRACER_BUDGET_HEAT, TRACER_BUDGET_SALT
   use rdb_mem_report, only: arr_bytes
   use pic_logger, only: global_logger
   use rdb_error_ring, only: error_ring_push
   implicit none
   private

   public :: multilayer_state_t

   type :: multilayer_state_t
      !! Per-layer multilayer C-grid state.

      logical :: is_init = .false.
         !! True between `init` and `destroy`.  Prefer this to
         !! `allocated(...)` — tracks GPU device attachment too.

      ! Conservation-budget accumulators (host scalars; not device-mapped).
      real(wp) :: mass_out = 0.0_wp
         !! Cumulative mass (kg) that has left the domain through its open
         !! boundaries since t=0 (positive = outflow), accumulated per RK2
         !! stage from the continuity divergence (`flux_h_layer`) so the
         !! console mass `Error` closes to round-off even with open BCs.
      logical :: mass_out_tracked = .false.
         !! Set once the dyn step has accumulated `mass_out`, so the console
         !! only activates the mass budget on a path that feeds it.

      integer :: nz_ml = 0
         !! Number of multilayer levels (k=1 bed, k=nz_ml surface).

      ! ---- Cell-centred per-layer thickness ----
      real(wp), allocatable :: h_layer(:, :, :)
         !! Layer thickness at cell centres (m), shape (nx, ny, nz_ml).
      real(wp), allocatable :: h_layer0(:, :, :)
         !! RK2 save of h_layer at start of outer step.

      ! ---- Face-located layer velocities + momenta ----
      ! Same face-indexing convention as the barotropic C-grid state.
      real(wp), allocatable :: u_face_x_layer(:, :, :)
         !! x-velocity at the WEST face of cell (i,j,k), shape
         !! (nx+1, ny, nz_ml): face i sits between cells i-1 and i, so a
         !! cell's divergence reads faces (i, i+1) — the convention every
         !! consumer (continuity `flux(i+1)-flux(i)`, `metrics%wet_u`,
         !! the kappa-shear centre average) actually uses.  ("east face"
         !! here previously was a stale docstring.)
      real(wp), allocatable :: hu_face_x_layer(:, :, :)
         !! h*u at the same west-face stagger (m^2/s).
      real(wp), allocatable :: v_face_y_layer(:, :, :)
         !! y-velocity at the SOUTH face of cell (i,j,k), shape
         !! (nx, ny+1, nz_ml): face j sits between cells j-1 and j
         !! (same stagger rule as `u_face_x_layer`).
      real(wp), allocatable :: hv_face_y_layer(:, :, :)
         !! h*v at the same south-face stagger.

      ! ---- MOM6 split-RK2 time-mean fields (docs/MOM6_SPLIT_RK2_SPEC.md §1) ----
      ! MOM6 carries an INSTANTANEOUS prognostic velocity AND a step time-mean
      ! (`CS%u_av`, `CS%h_av`).  Every slow tendency (CorAdCalc, horizontal
      ! viscosity) is evaluated on the TIME-MEAN, never on the prognostic —
      ! that is where the predictor-corrector gets its second-order character
      ! without an SSP stage average.  `u_av` is produced by continuity's
      ! `u_cor` (the transport-matched velocity satisfying Sum_k u*h = uhbt)
      ! and MUST NOT be written back into the prognostic.
      real(wp), allocatable :: u_av_layer(:, :, :)
         !! Step time-mean x face velocity, shape (nx+1, ny, nz_ml). MOM6 `u_av`.
      real(wp), allocatable :: v_av_layer(:, :, :)
         !! Step time-mean y face velocity, shape (nx, ny+1, nz_ml). MOM6 `v_av`.
      real(wp), allocatable :: h_av_layer(:, :, :)
         !! Step time-mean layer thickness, shape (nx, ny, nz_ml). MOM6 `h_av`.

      ! ---- RK2 face-velocity saves ----
      real(wp), allocatable :: u_face_x_layer0(:, :, :)
      real(wp), allocatable :: v_face_y_layer0(:, :, :)

      ! ---- Per-face per-layer mass fluxes (continuity-PPM output) ----
      real(wp), allocatable :: mass_flux_x_layer(:, :, :)
      real(wp), allocatable :: mass_flux_y_layer(:, :, :)

      ! ---- Cell-centred per-layer flux divergence ----
      ! continuity-PPM writes here; the apply step reads it for the
      ! forward-Euler h_layer update.  Shape (nx, ny, nz_ml).
      real(wp), allocatable :: flux_h_layer(:, :, :)

      ! ---- Cell-centred density ----
      ! EOS writes here from T, S; the PGF kernel reads it.  Shape
      ! (nx, ny, nz_ml), matches h_layer.
      real(wp), allocatable :: rho_layer(:, :, :)

      ! ---- Vertical (cross-layer) velocity at layer interfaces ----
      ! Diagnostic in z*; reads as residual w under ALE.  Shape
      ! (nx, ny, nz_ml+1) with k=1 the bed (0) and k=nz_ml+1 the surface.
      real(wp), allocatable :: w_interface(:, :, :)

      ! ---- Land / ocean mask (surface-forcing mask) ----
      ! 2D wet-cell indicator at cell centres: 1.0 = ocean, 0.0 = land.
      ! Populated at IC time from the bathymetry threshold; consumed by
      ! the three surface-forcing kernels (`ocean_surface_stress`,
      ! `ocean_bottom_drag`, `ocean_surface_flux`) so forcing doesn't
      ! drive fictitious land-column currents.  Default 1.0 everywhere
      ! is bit-identical.  Not yet plumbed through the prognostic
      ! kernels (continuity, PGF, advection, viscosity, vdiff).
      real(wp), allocatable :: wet_mask(:, :)

      ! ---- Tracer registry ----
      ! Reuses the coastal `tracer_t` verbatim, cell-centred shape
      ! `(nx, ny, nz)` (stagger only affects velocities).  Salinity +
      ! temperature first; passive tracers append.  Special-shape
      ! kernels (EOS, surface flux) locate S/T via the named indices.
      type(tracer_t), allocatable :: tracers(:)
         !! Registered prognostic tracers.
      integer :: idx_salinity = 0
         !! Index into tracers(:) for salinity. 0 = not registered.
      integer :: idx_temperature = 0
         !! Index into tracers(:) for temperature. 0 = not registered.
      integer :: idx_age = 0
         !! Index into tracers(:) for the ideal-age tracer.  0 = not
         !! registered.  When > 0, the dyn step ages at 1 s/s and zeros
         !! the surface layer (k = nz_ml) every step.  No EOS coupling.
      integer :: idx_pseudo_salt = 0
         !! Index into tracers(:) for the pseudo-salt verification tracer.
         !! 0 = not registered.
      logical :: registry_locked = .false.
         !! Set by `enter_data`, cleared by `exit_data`.  While locked,
         !! `register_passive_tracer` REFUSES: the device map snapshots
         !! `tracers(:)` element-by-element, so a slot appended after the
         !! map has no device `hTr` and the first kernel touching it faults
         !! on the `mem:separate` build.

      ! ---- Budget contributor slots ----
      ! Per-kernel mass / heat / salt accumulators: each physics kernel
      ! writes its per-step `dt · delta` here; the budget module's
      ! `drain_contributors` spatially integrates + zeroes them at each
      ! eval cadence.
      real(wp), allocatable :: mass_budget_continuity(:, :, :)
         !! Mass change per cell per step attributed to continuity-PPM
         !! divergence (m·dt units; the budget integral over volume
         !! recovers m³).  Shape (nx, ny, nz_ml).  Zero in a closed
         !! basin (perfect telescope).
      real(wp), allocatable :: heat_budget_surface(:, :, :)
         !! hTr (K·m) change per cell per step attributed to the surface
         !! heat-flux kernel.  Populated only at k=nz_ml (the surface
         !! layer); zero elsewhere.  Sign convention: positive = source
         !! into the ocean.
      real(wp), allocatable :: salt_budget_surface(:, :, :)
         !! hTr (PSU·m) change per cell per step attributed to the
         !! surface salt-flux kernel.  Same shape + indexing convention
         !! as `heat_budget_surface`.
      real(wp), allocatable :: heat_budget_geothermal(:, :, :)
         !! hTr (K·m) change per cell per step attributed to the
         !! geothermal bottom-heat-flux kernel.  Populated only at the
         !! lowest massive layer (k=1 in the common case); zero
         !! elsewhere.  Sign convention: positive = source into the
         !! ocean from below.
      real(wp), allocatable :: heat_budget_sponge(:, :, :)
         !! hTr (K·m) change per cell per layer per step attributed to
         !! the map-driven sponge's tracer relaxation
         !! (`rdb_ocean_sponge::relax_tracer_budget_impl`).  Zero unless
         !! `&ocean_sponge_nml enable=.true., relax_tracers=.true.`.
         !! Sign convention: positive = source into the ocean (relaxing
         !! toward a warmer reference).  Not drained by the sponge.
      real(wp), allocatable :: salt_budget_sponge(:, :, :)
         !! hTr (PSU·m) change per cell per layer per step attributed to
         !! the map-driven sponge's tracer relaxation.  Same shape +
         !! indexing convention as `heat_budget_sponge`.
      real(wp), allocatable :: heat_budget_vert_adv(:, :, :)
         !! hTr (K·m) change per cell per layer per step attributed to
         !! the first-order-upwind vertical-advection kernel for
         !! temperature.  Closed-BC kernel: column sum telescopes to
         !! zero, so the spatial integral is zero to FP.
      real(wp), allocatable :: salt_budget_vert_adv(:, :, :)
         !! hTr (PSU·m) change per cell per layer per step for
         !! salinity.  Same closed-BC telescope property.
      real(wp), allocatable :: heat_budget_vdiff(:, :, :)
         !! hTr (K·m) change per cell per layer per step attributed to
         !! the backward-Euler vertical-diffusion tridiag solve for
         !! temperature.  Closed BCs (no flux through bed or surface)
         !! ⇒ column sum telescopes to zero.
      real(wp), allocatable :: salt_budget_vdiff(:, :, :)
         !! Salinity analogue of `heat_budget_vdiff`.
      real(wp), allocatable :: heat_budget_hdiff(:, :, :)
         !! hTr (K·m) change per cell per layer per step attributed to
         !! the Laplacian horizontal-diffusion kernel for temperature.
         !! Closed walls (wall faces forced to zero) ⇒ spatial integral
         !! telescopes to zero.
      real(wp), allocatable :: salt_budget_hdiff(:, :, :)
         !! Salinity analogue of `heat_budget_hdiff`.
      real(wp), allocatable :: heat_budget_horiz_adv(:, :, :)
         !! hTr (°C·m) change per cell per layer, accumulated (summed)
         !! across BOTH RK2 stages from t=0 (`+=` each stage, never zeroed
         !! mid-step; no 0.5 weight here — the console applies it),
         !! attributed to the continuity-PPM HORIZONTAL tracer advection
         !! (zonal + meridional divergence of the tracer mass flux).
         !! Interior sum telescopes to the net advective flux across the
         !! open boundaries; zero in a closed basin.  Only the fused
         !! (`dt_tracer_advect_ratio=1`) path fills it — see the console
         !! reporter for the >1 fallback.
      real(wp), allocatable :: salt_budget_horiz_adv(:, :, :)
         !! Salinity analogue (PSU·m) of `heat_budget_horiz_adv`.
      real(wp), allocatable :: mass_budget_remap(:, :, :)
         !! Per-cell `h_layer_new − h_layer_old` from the ALE remap
         !! step.  Column-conservative ⇒ Σ_k = 0 per (i, j).  Non-zero
         !! residual flags a remap conservation leak.
      real(wp), allocatable :: heat_budget_remap(:, :, :)
         !! Per-cell `hTr_new − hTr_old` from the ALE remap step for
         !! temperature.  Same column-telescope property.
      real(wp), allocatable :: salt_budget_remap(:, :, :)
         !! Salinity analogue.

   contains
      procedure, non_overridable :: init => multilayer_state_init
      procedure, non_overridable :: destroy => multilayer_state_destroy
      procedure, non_overridable :: enter_data => multilayer_state_enter_data
      procedure, non_overridable :: exit_data => multilayer_state_exit_data
      procedure, non_overridable :: bytes => multilayer_state_bytes
      procedure, non_overridable :: register_passive_tracer => &
         multilayer_register_passive_tracer
   end type multilayer_state_t

contains

   subroutine multilayer_state_init(this, grid, with_ideal_age)
      !! Allocate per-layer C-grid arrays at the grid size and the
      !! configured layer count (caller must set `this%nz_ml` first).
      !! Registers salinity + temperature with default identity strings.
      !! When `with_ideal_age` is present and true, also registers an
      !! ideal-age tracer at index 3 (see `rdb_ocean_ideal_age`).
      class(multilayer_state_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      logical, intent(in), optional :: with_ideal_age

      integer :: nx, ny, nz_ml, ntracers
      logical :: age_on

      age_on = .false.
      if (present(with_ideal_age)) age_on = with_ideal_age

      nx = grid%nx_total
      ny = grid%ny_total
      nz_ml = this%nz_ml

      ! Cell-centred per-layer arrays
      allocate (this%h_layer(nx, ny, nz_ml), source=0.0_wp)
      allocate (this%h_layer0(nx, ny, nz_ml), source=0.0_wp)
      allocate (this%h_av_layer(nx, ny, nz_ml), source=0.0_wp)

      ! East-face per-layer arrays
      allocate (this%u_face_x_layer(nx + 1, ny, nz_ml), source=0.0_wp)
      allocate (this%hu_face_x_layer(nx + 1, ny, nz_ml), source=0.0_wp)
      allocate (this%mass_flux_x_layer(nx + 1, ny, nz_ml), source=0.0_wp)
      allocate (this%u_face_x_layer0(nx + 1, ny, nz_ml), source=0.0_wp)
      allocate (this%u_av_layer(nx + 1, ny, nz_ml), source=0.0_wp)

      ! North-face per-layer arrays
      allocate (this%v_face_y_layer(nx, ny + 1, nz_ml), source=0.0_wp)
      allocate (this%hv_face_y_layer(nx, ny + 1, nz_ml), source=0.0_wp)
      allocate (this%mass_flux_y_layer(nx, ny + 1, nz_ml), source=0.0_wp)
      allocate (this%v_face_y_layer0(nx, ny + 1, nz_ml), source=0.0_wp)
      allocate (this%v_av_layer(nx, ny + 1, nz_ml), source=0.0_wp)

      ! Cell-centred per-layer flux divergence (pure workspace)
      allocate (this%flux_h_layer(nx, ny, nz_ml), source=0.0_wp)

      ! Budget contributor slots, drained at every budget eval cadence.
      allocate (this%mass_budget_continuity(nx, ny, nz_ml), source=0.0_wp)
      allocate (this%heat_budget_surface(nx, ny, nz_ml), source=0.0_wp)
      allocate (this%salt_budget_surface(nx, ny, nz_ml), source=0.0_wp)
      allocate (this%heat_budget_geothermal(nx, ny, nz_ml), source=0.0_wp)
      allocate (this%heat_budget_sponge(nx, ny, nz_ml), source=0.0_wp)
      allocate (this%salt_budget_sponge(nx, ny, nz_ml), source=0.0_wp)
      allocate (this%heat_budget_vert_adv(nx, ny, nz_ml), source=0.0_wp)
      allocate (this%salt_budget_vert_adv(nx, ny, nz_ml), source=0.0_wp)
      allocate (this%heat_budget_vdiff(nx, ny, nz_ml), source=0.0_wp)
      allocate (this%salt_budget_vdiff(nx, ny, nz_ml), source=0.0_wp)
      allocate (this%heat_budget_hdiff(nx, ny, nz_ml), source=0.0_wp)
      allocate (this%salt_budget_hdiff(nx, ny, nz_ml), source=0.0_wp)
      allocate (this%heat_budget_horiz_adv(nx, ny, nz_ml), source=0.0_wp)
      allocate (this%salt_budget_horiz_adv(nx, ny, nz_ml), source=0.0_wp)
      allocate (this%mass_budget_remap(nx, ny, nz_ml), source=0.0_wp)
      allocate (this%heat_budget_remap(nx, ny, nz_ml), source=0.0_wp)
      allocate (this%salt_budget_remap(nx, ny, nz_ml), source=0.0_wp)

      ! Cell-centred density (filled by the EOS each step)
      allocate (this%rho_layer(nx, ny, nz_ml), source=0.0_wp)

      ! Cross-layer vertical velocity at interfaces (k=1 bed .. k=nz+1 surface)
      allocate (this%w_interface(nx, ny, nz_ml + 1), source=0.0_wp)

      ! Wet-cell mask: default to all-ocean (1.0) so analytical / flat-
      ! bottom tests don't see a behaviour change.  Realistic-bathy runs
      ! overwrite this in `ocean_state_seed_from_cfg` after the bathy load.
      allocate (this%wet_mask(nx, ny), source=1.0_wp)

      ! Tracer registry: salinity at index 1, temperature at index 2,
      ! ideal-age at index 3 (optional).
      ntracers = 2
      if (age_on) ntracers = 3
      allocate (this%tracers(ntracers))
      this%idx_salinity = 1
      this%idx_temperature = 2
      call this%tracers(this%idx_salinity)%init(grid, nz_ml)
      call this%tracers(this%idx_temperature)%init(grid, nz_ml)
      this%tracers(this%idx_salinity)%name = "salinity"
      this%tracers(this%idx_salinity)%long_name = "Sea water salinity"
      this%tracers(this%idx_salinity)%units = "PSU"
      this%tracers(this%idx_salinity)%standard_name = "sea_water_salinity"
      this%tracers(this%idx_salinity)%budget_id = TRACER_BUDGET_SALT
      this%tracers(this%idx_temperature)%name = "temperature"
      this%tracers(this%idx_temperature)%long_name = "Sea water potential temperature"
      this%tracers(this%idx_temperature)%units = "degC"
      this%tracers(this%idx_temperature)%standard_name = "sea_water_potential_temperature"
      this%tracers(this%idx_temperature)%budget_id = TRACER_BUDGET_HEAT
      if (age_on) then
         this%idx_age = 3
         call this%tracers(this%idx_age)%init(grid, nz_ml)
         this%tracers(this%idx_age)%name = "age"
         this%tracers(this%idx_age)%long_name = "Ideal age of sea water"
         this%tracers(this%idx_age)%units = "s"
         this%tracers(this%idx_age)%standard_name = "age_of_sea_water"
      end if

      this%is_init = .true.
   end subroutine multilayer_state_init

   subroutine multilayer_state_destroy(this)
      class(multilayer_state_t), intent(inout) :: this
      integer :: i
      this%is_init = .false.
      if (allocated(this%tracers)) then
         do i = 1, size(this%tracers)
            call this%tracers(i)%destroy()
         end do
         deallocate (this%tracers)
      end if
      this%idx_salinity = 0
      this%idx_temperature = 0
      this%idx_age = 0
      this%idx_pseudo_salt = 0
      this%registry_locked = .false.
      if (allocated(this%h_layer)) deallocate (this%h_layer)
      if (allocated(this%h_layer0)) deallocate (this%h_layer0)
      if (allocated(this%u_face_x_layer)) deallocate (this%u_face_x_layer)
      if (allocated(this%hu_face_x_layer)) deallocate (this%hu_face_x_layer)
      if (allocated(this%v_face_y_layer)) deallocate (this%v_face_y_layer)
      if (allocated(this%hv_face_y_layer)) deallocate (this%hv_face_y_layer)
      if (allocated(this%u_face_x_layer0)) deallocate (this%u_face_x_layer0)
      if (allocated(this%u_av_layer)) deallocate (this%u_av_layer)
      if (allocated(this%v_av_layer)) deallocate (this%v_av_layer)
      if (allocated(this%h_av_layer)) deallocate (this%h_av_layer)
      if (allocated(this%v_face_y_layer0)) deallocate (this%v_face_y_layer0)
      if (allocated(this%mass_flux_x_layer)) deallocate (this%mass_flux_x_layer)
      if (allocated(this%mass_flux_y_layer)) deallocate (this%mass_flux_y_layer)
      if (allocated(this%flux_h_layer)) deallocate (this%flux_h_layer)
      if (allocated(this%rho_layer)) deallocate (this%rho_layer)
      if (allocated(this%w_interface)) deallocate (this%w_interface)
      if (allocated(this%mass_budget_continuity)) deallocate (this%mass_budget_continuity)
      if (allocated(this%heat_budget_surface)) deallocate (this%heat_budget_surface)
      if (allocated(this%salt_budget_surface)) deallocate (this%salt_budget_surface)
      if (allocated(this%heat_budget_geothermal)) deallocate (this%heat_budget_geothermal)
      if (allocated(this%heat_budget_sponge)) deallocate (this%heat_budget_sponge)
      if (allocated(this%salt_budget_sponge)) deallocate (this%salt_budget_sponge)
      if (allocated(this%heat_budget_vert_adv)) deallocate (this%heat_budget_vert_adv)
      if (allocated(this%salt_budget_vert_adv)) deallocate (this%salt_budget_vert_adv)
      if (allocated(this%heat_budget_vdiff)) deallocate (this%heat_budget_vdiff)
      if (allocated(this%salt_budget_vdiff)) deallocate (this%salt_budget_vdiff)
      if (allocated(this%heat_budget_hdiff)) deallocate (this%heat_budget_hdiff)
      if (allocated(this%salt_budget_hdiff)) deallocate (this%salt_budget_hdiff)
      if (allocated(this%heat_budget_horiz_adv)) deallocate (this%heat_budget_horiz_adv)
      if (allocated(this%salt_budget_horiz_adv)) deallocate (this%salt_budget_horiz_adv)
      if (allocated(this%mass_budget_remap)) deallocate (this%mass_budget_remap)
      if (allocated(this%heat_budget_remap)) deallocate (this%heat_budget_remap)
      if (allocated(this%salt_budget_remap)) deallocate (this%salt_budget_remap)
      if (allocated(this%wet_mask)) deallocate (this%wet_mask)
   end subroutine multilayer_state_destroy

   pure function multilayer_state_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the ocean C-grid layer slot: the
      !! layer prognostics + face transports + RK2 saves + density/vertical
      !! diagnostics, the per-tracer registry (each tracer sums its own
      !! arrays; ideal-age rides the registry when on), and the
      !! device-resident conservation-budget accumulators.
      class(multilayer_state_t), intent(in) :: this
      integer(int64) :: nbytes
      integer :: it

      nbytes = arr_bytes(this%h_layer) + arr_bytes(this%h_layer0) &
               + arr_bytes(this%u_face_x_layer) + arr_bytes(this%hu_face_x_layer) &
               + arr_bytes(this%v_face_y_layer) + arr_bytes(this%hv_face_y_layer) &
               + arr_bytes(this%u_face_x_layer0) + arr_bytes(this%v_face_y_layer0) &
               + arr_bytes(this%u_av_layer) + arr_bytes(this%v_av_layer) &
               + arr_bytes(this%h_av_layer) &
               + arr_bytes(this%mass_flux_x_layer) + arr_bytes(this%mass_flux_y_layer) &
               + arr_bytes(this%flux_h_layer) + arr_bytes(this%rho_layer) &
               + arr_bytes(this%w_interface) + arr_bytes(this%wet_mask) &
               + arr_bytes(this%mass_budget_continuity) &
               + arr_bytes(this%heat_budget_surface) + arr_bytes(this%salt_budget_surface) &
               + arr_bytes(this%heat_budget_geothermal) &
               + arr_bytes(this%heat_budget_sponge) + arr_bytes(this%salt_budget_sponge) &
               + arr_bytes(this%heat_budget_vert_adv) + arr_bytes(this%salt_budget_vert_adv) &
               + arr_bytes(this%heat_budget_vdiff) + arr_bytes(this%salt_budget_vdiff) &
               + arr_bytes(this%heat_budget_hdiff) + arr_bytes(this%salt_budget_hdiff) &
               + arr_bytes(this%heat_budget_horiz_adv) + arr_bytes(this%salt_budget_horiz_adv) &
               + arr_bytes(this%mass_budget_remap) &
               + arr_bytes(this%heat_budget_remap) + arr_bytes(this%salt_budget_remap)

      if (allocated(this%tracers)) then
         do it = 1, size(this%tracers)
            nbytes = nbytes + this%tracers(it)%bytes()
         end do
      end if
   end function multilayer_state_bytes

   subroutine multilayer_register_passive_tracer(this, grid, name, units, &
                                                 long_name, idx)
      !! Append a passive tracer (`eos_coeff = 0`, `budget_id = NONE`) to
      !! the registry, growing `tracers(:)` past the default S/T[/age] set.
      !! Returns its slot in `idx`, or `idx = 0` on refusal (registry not
      !! init'd, or locked by `enter_data`).  MUST be called after `init`
      !! and BEFORE `enter_data` — and, on the ocean path, before
      !! `ocean_bc_state_init` sizes `bc%n_tracers`.  Caller populates
      !! `hTr` once layer thicknesses exist, and may set
      !! `tracers(idx)%standard_name` / the pipeline opt-outs directly
      !! (public components).  S/T/age keep their indices.
      class(multilayer_state_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      character(len=*), intent(in) :: name
      character(len=*), intent(in) :: units
      character(len=*), intent(in) :: long_name
      integer, intent(out) :: idx

      type(tracer_t), allocatable :: tmp(:)
      integer :: n_old, it

      idx = 0

      if (.not. this%is_init) then
         call error_ring_push("multilayer_state: register_passive_tracer("// &
                              trim(name)//") refused — registry not initialised "// &
                              "(call after init())")
         call global_logger%error("multilayer_state: register_passive_tracer("// &
                                  trim(name)//") refused — registry not initialised "// &
                                  "(call after init())")
         return
      end if

      if (this%registry_locked) then
         call error_ring_push("multilayer_state: register_passive_tracer("// &
                              trim(name)//") refused — registry locked by "// &
                              "enter_data (a slot appended now would never be "// &
                              "device-mapped on the mem:separate GPU build)")
         call global_logger%error("multilayer_state: register_passive_tracer("// &
                                  trim(name)//") refused — registry locked by "// &
                                  "enter_data (a slot appended now would never be "// &
                                  "device-mapped on the mem:separate GPU build)")
         return
      end if

      n_old = size(this%tracers)
      call move_alloc(this%tracers, tmp)
      allocate (this%tracers(n_old + 1))
      ! Deep-copy the existing tracers (intrinsic assignment reallocates
      ! the allocatable hTr / hTr0 components and copies their data).
      do it = 1, n_old
         this%tracers(it) = tmp(it)
      end do

      it = n_old + 1
      call this%tracers(it)%init(grid, this%nz_ml)
      this%tracers(it)%name = trim(name)
      this%tracers(it)%long_name = trim(long_name)
      this%tracers(it)%units = trim(units)
      this%tracers(it)%standard_name = ""
      this%tracers(it)%eos_coeff = 0.0_wp
      this%tracers(it)%eos_ref = 0.0_wp
      this%tracers(it)%budget_id = TRACER_BUDGET_NONE

      idx = it
   end subroutine multilayer_register_passive_tracer

   subroutine multilayer_state_enter_data(this)
      !! Attach the C-grid multilayer allocatables to the device.  The
      !! tracer registry uses the two-step pattern: array descriptor
      !! first, then each element's hTr / hTr0 — NVHPC stdpar can't
      !! dereference `tracers(it)%hTr` from a do-concurrent body
      !! otherwise.
      class(multilayer_state_t), intent(inout) :: this
      select type (this)
      type is (multilayer_state_t)
         call multilayer_state_enter_data_impl(this)
      end select
   end subroutine multilayer_state_enter_data

   subroutine multilayer_state_enter_data_impl(this)
      type(multilayer_state_t), intent(inout) :: this
      integer :: it

      !$acc enter data copyin(this%h_layer, this%h_layer0, &
      !$acc&                  this%u_face_x_layer, this%hu_face_x_layer, &
      !$acc&                  this%v_face_y_layer, this%hv_face_y_layer, &
      !$acc&                  this%u_face_x_layer0, this%v_face_y_layer0, &
      !$acc&                  this%u_av_layer, this%v_av_layer, this%h_av_layer, &
      !$acc&                  this%w_interface)
      !$acc enter data create(this%mass_flux_x_layer, this%mass_flux_y_layer, &
      !$acc&                  this%flux_h_layer)
      !$acc enter data copyin(this%rho_layer, this%mass_budget_continuity, &
      !$acc&                  this%heat_budget_surface, this%salt_budget_surface, &
      !$acc&                  this%heat_budget_geothermal, &
      !$acc&                  this%heat_budget_sponge, this%salt_budget_sponge, &
      !$acc&                  this%heat_budget_vert_adv, this%salt_budget_vert_adv, &
      !$acc&                  this%heat_budget_vdiff, this%salt_budget_vdiff, &
      !$acc&                  this%heat_budget_hdiff, this%salt_budget_hdiff, &
      !$acc&                  this%heat_budget_horiz_adv, this%salt_budget_horiz_adv, &
      !$acc&                  this%mass_budget_remap, this%heat_budget_remap, &
      !$acc&                  this%salt_budget_remap, this%wet_mask)

      if (allocated(this%tracers)) then
         !$acc enter data copyin(this%tracers)
         do it = 1, size(this%tracers)
            !$acc enter data copyin(this%tracers(it)%hTr, &
            !$acc&                  this%tracers(it)%hTr0)
         end do
      end if

      this%registry_locked = .true.
   end subroutine multilayer_state_enter_data_impl

   subroutine multilayer_state_exit_data(this)
      !! Reverse of `enter_data`.  Copy out the prognostic fields and
      !! the tracer hTr arrays (so post-run host inspection works),
      !! drop scratch + RK saves.  Tracer registry tears down per-
      !! element first, then the array descriptor — mirror of
      !! `enter_data` order.
      class(multilayer_state_t), intent(inout) :: this
      select type (this)
      type is (multilayer_state_t)
         call multilayer_state_exit_data_impl(this)
      end select
   end subroutine multilayer_state_exit_data

   subroutine multilayer_state_exit_data_impl(this)
      type(multilayer_state_t), intent(inout) :: this
      integer :: it

      if (allocated(this%tracers)) then
         do it = 1, size(this%tracers)
            !$acc exit data copyout(this%tracers(it)%hTr)
            !$acc exit data delete(this%tracers(it)%hTr0)
         end do
         !$acc exit data delete(this%tracers)
      end if

      this%registry_locked = .false.

      !$acc exit data copyout(this%h_layer, &
      !$acc&                  this%u_face_x_layer, this%v_face_y_layer, &
      !$acc&                  this%hu_face_x_layer, this%hv_face_y_layer, &
      !$acc&                  this%w_interface)
      !$acc exit data copyout(this%rho_layer)
      !$acc exit data delete(this%h_layer0, &
      !$acc&                 this%u_av_layer, this%v_av_layer, this%h_av_layer, &
      !$acc&                 this%u_face_x_layer0, this%v_face_y_layer0, &
      !$acc&                 this%mass_flux_x_layer, this%mass_flux_y_layer, &
      !$acc&                 this%flux_h_layer, this%mass_budget_continuity, &
      !$acc&                 this%heat_budget_surface, this%salt_budget_surface, &
      !$acc&                 this%heat_budget_geothermal, &
      !$acc&                 this%heat_budget_sponge, this%salt_budget_sponge, &
      !$acc&                 this%heat_budget_vert_adv, this%salt_budget_vert_adv, &
      !$acc&                 this%heat_budget_vdiff, this%salt_budget_vdiff, &
      !$acc&                 this%heat_budget_hdiff, this%salt_budget_hdiff, &
      !$acc&                 this%heat_budget_horiz_adv, this%salt_budget_horiz_adv, &
      !$acc&                 this%mass_budget_remap, this%heat_budget_remap, &
      !$acc&                 this%salt_budget_remap, this%wet_mask)
   end subroutine multilayer_state_exit_data_impl

end module rdb_multilayer_state
