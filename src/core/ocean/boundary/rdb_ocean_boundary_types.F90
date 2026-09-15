!! Ocean open-boundary condition types.
module rdb_ocean_boundary_types
   !! Type taxonomy + per-edge BC config + the composed slot on `ocean_state_t`.
   !! Six BC types implemented end-to-end (WALL, OPEN, TIDAL, CLAMPED, SPONGE,
   !! CHAPMAN); INFLOW/DISCHARGE/NESTED tags are declared for cross-backend
   !! alignment but `error stop` if encountered. Per-edge granularity: each of
   !! the four outer edges carries one `ocean_bc_face_tag_t`, read independently
   !! by the dispatch helpers in `rdb_ocean_boundary`.
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
   use rdb_ocean_tide_astro, only: TIDE_OMEGA, TIDES_CATALOG_SIZE
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_mem_report, only: arr_bytes
   use pic_ascii, only: to_lower
   use pic_logger, only: logger => global_logger
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_SETUP
   implicit none
   private

   public :: ocean_bc_face_tag_t
   public :: ocean_bc_state_t
   public :: ocean_bc_state_init, ocean_bc_state_destroy
   public :: ocean_bc_state_enter_data, ocean_bc_state_exit_data
   public :: ocean_bc_state_set_edges
   public :: ocean_bc_state_set_topology
   public :: ocean_bc_type_from_string
   public :: ocean_bc_has_tracer_open_edge
   public :: ocean_bc_validate_periodic
   public :: ocean_bc_validate_fold
   public :: obc_match_constituent
   public :: obc_tide_nodal_fill

   real(wp), parameter, public :: OBC_TIDE_MATCH_TOL = 1.0e-4_wp
      !! Relative tolerance for matching an OBC edge constituent's angular
      !! frequency to a catalog entry (`rdb_ocean_tide_astro::TIDE_OMEGA`).
      !! The catalog ω's are well separated (min gap S2↔K2 ≈ 0.3 %), so any
      !! physical constituent resolves unambiguously; a farther-than-tol
      !! nearest match signals an unknown constituent (fail-loud at setup).

   ! ---- BC type tags ----
   ! Numbered to mirror the coastal taxonomy (rdb_boundary_types.F90) so
   ! config strings parse to the same integer on either backend.
   integer, parameter, public :: OBC_WALL = 1
      !! Closed wall (hard-zero). Default for every edge.
   integer, parameter, public :: OBC_OPEN = 2
      !! Flather radiation — gravity-wave outflow + η clamped to a reference
      !! (zero by default, supplied via `data_eta_*`).
   integer, parameter, public :: OBC_TIDAL = 3
      !! Prescribed multi-constituent η; composed into `data_eta_*`.
   integer, parameter, public :: OBC_NESTED = 4
      !! Two-way nesting (not yet implemented). Behaves like OPEN to kernels.
   integer, parameter, public :: OBC_INFLOW = 5
      !! Prescribed normal velocity + tracer. Cross-backend symmetry only.
   integer, parameter, public :: OBC_DISCHARGE = 6
      !! Prescribed volume flux. Cross-backend symmetry only.
   integer, parameter, public :: OBC_CLAMPED = 7
      !! Hard Dirichlet on η + u + v + per-tracer values, sourced from `data_*`.
   integer, parameter, public :: OBC_SPONGE = 8
      !! Relaxation band — BC kernel falls through to WALL at the outer face;
      !! the sponge kernel relaxes the interior band toward `data_*` targets.
   integer, parameter, public :: OBC_CHAPMAN = 9
      !! Orlanski radiation on η with implicit phase-speed estimation. Uses
      !! persistent `eta_old_<edge>` state across timesteps.
   integer, parameter, public :: OBC_PERIODIC = 10
      !! Ghost-wrap periodic boundary: ghost columns/rows hold copies of the
      !! opposite interior so kernels see a seamless domain. Requires
      !! `nghost >= 3` (PPM + biharmonic stencil) and must be paired
      !! (west ⟺ east; south ⟺ north). Cannot combine with OBC_SPONGE on the
      !! same edge. Keep tag in sync with coastal rdb_boundary_types.F90.
   integer, parameter, public :: OBC_TRIPOLAR_FOLD = 11
      !! Tripolar north-fold seam (Murray 1996). NORTH edge only. The fold
      !! exchange (`rdb_ocean_fold`) halo-fills the north ghost rows
      !! (reversed-i, sign-flipped for vector normals) and antisymmetrically
      !! projects the on-line v/corner row. Requires `grid_config="tripolar"`
      !! and periodic west/east (fold reads already-wrapped corner columns).
   integer, parameter, public :: OBC_INVALID = -1
      !! Sentinel returned by `ocean_bc_type_from_string` for an
      !! unrecognised edge string (PR-6 fail-loud).  A misspelled edge
      !! must NOT silently close the boundary to a wall — the user asked
      !! for a specific (often open) boundary and a wall reflects every
      !! outgoing gravity wave, a materially different model with no
      !! message.  `validate_config` rejects it before `configure_ocean_bc`
      !! ever consumes the parse result.

   ! Multi-constituent tidal forcing — width matches coastal so one
   ! TPXO/FES table feeds both backends.
   integer, parameter, public :: OBC_MAX_TIDAL_CONSTITUENTS = 8

   ! ---- Per-edge tag ----
   type :: ocean_bc_face_tag_t
      !! Configuration for one outer edge. Defaults yield a closed wall.
      integer :: bc_type = OBC_WALL

      ! Clamped (Dirichlet) reference values
      real(wp) :: clamped_eta = 0.0_wp
      real(wp) :: clamped_u = 0.0_wp
      real(wp) :: clamped_v = 0.0_wp
      real(wp), allocatable :: clamped_tracer(:)
         !! Per-tracer Dirichlet values for CLAMPED inflow.  Size
         !! `n_tracers`; only consulted on inflow faces.

      ! Tidal constituents (per-edge — different tide at each open edge)
      integer  :: n_tidal_constituents = 0
      real(wp) :: tidal_amp(OBC_MAX_TIDAL_CONSTITUENTS) = 0.0_wp
      real(wp) :: tidal_phase(OBC_MAX_TIDAL_CONSTITUENTS) = 0.0_wp
      real(wp) :: tidal_omega(OBC_MAX_TIDAL_CONSTITUENTS) = 0.0_wp
      ! Nodal/astronomical correction (capability C3).  Baked once at setup
      ! from the shared tide reference epoch when `ocean_bc_state_t%tidal_nodal`
      ! is on; the defaults (f=1, arg=0) leave the legacy static-phase OBC sum
      ! bit-identical.  Fixed-size members ⇒ ride the parent bc GPU mapping.
      real(wp) :: tidal_fnodal(OBC_MAX_TIDAL_CONSTITUENTS) = 1.0_wp
         !! 18.6-yr nodal amplitude factor `f_c` per edge constituent.
      real(wp) :: tidal_arg(OBC_MAX_TIDAL_CONSTITUENTS) = 0.0_wp
         !! Equilibrium + nodal phase `(V_c + u_c)` (rad) per edge constituent.

      ! Sponge config (only meaningful when bc_type == OBC_SPONGE).
      integer  :: sponge_width = 0
      real(wp) :: sponge_strength = 0.0_wp
      logical  :: sponge_relax_tracers = .false.
         !! When .true. the legacy band sponge relaxes tracer `hTr`
         !! (concentration held toward `clamped_tracer(:)`, mass `h_layer`
         !! left untouched) in the edge band. Default .false. is
         !! bit-identical. `h_layer` relaxation is a separate, map-driven
         !! capability gated by `&ocean_sponge_nml relax_h` (PR-23b),
         !! restricted to `VCOORD_LAGRANGIAN` — see `rdb_ocean_sponge.F90`.
   end type ocean_bc_face_tag_t

   ! ---- Composed slot on ocean_state_t ----
   type :: ocean_bc_state_t
      !! Per-state OBC bookkeeping. Composed onto `ocean_state_t`; the
      !! dispatch helpers in `rdb_ocean_boundary` consume it via `class(*)`
      !! polymorphism, keeping kernels decoupled from the full state.
      logical :: is_init = .false.

      ! Per-edge tag (the load-bearing config).
      type(ocean_bc_face_tag_t) :: west, east, south, north

      ! Derived periodic flags — cached at init from the edge tags.
      ! Kernels read these two logicals (never re-derive from tags per step).
      logical :: periodic_x = .false.
         !! True when west and east edges are both OBC_PERIODIC.
      logical :: periodic_y = .false.
         !! True when south and north edges are both OBC_PERIODIC.
      logical :: north_fold = .false.
         !! True when the north edge is OBC_TRIPOLAR_FOLD. Gates every fold
         !! exchange in the dyn loop + BT substep; default .false. is
         !! bit-identical.
      logical :: has_west = .true.
         !! False when the west edge of this subdomain is an MPI seam (a
         !! neighbouring rank owns the cells beyond it), true when it is a
         !! physical domain edge.  Set from decomp%has_west at BC configure.
         !! Default .true. => single-rank / physical-edge behaviour
         !! (bit-identical to the pre-decomp code).
      logical :: has_east = .true.
         !! False when the east edge of this subdomain is an MPI seam; true
         !! when it is a physical domain edge.  (analogous to has_west)
      logical :: has_south = .true.
         !! False when the south edge of this subdomain is an MPI seam; true
         !! when it is a physical domain edge.  (analogous to has_west)
      logical :: has_north = .true.
         !! False when the north edge of this subdomain is an MPI seam; true
         !! when it is a physical domain edge.  (analogous to has_west)
      logical :: tidal_nodal = .false.
         !! Global switch (capability C3): apply the 18.6-yr nodal factor `f_c`
         !! + equilibrium/nodal phase `(V_c + u_c)` to the OBC tidal elevation
         !! forcing. Baked into the per-edge `tidal_fnodal`/`tidal_arg` at setup
         !! from the shared `&ocean_tides_nml` reference epoch. Cached into the
         !! barotropic substep like `periodic_x`. Default .false. ⇒ bit-identical
         !! legacy static-phase sum. When .true. the phase convention flips: the
         !! Greenwich phase `tidal_phase` becomes a LAG (subtracted).

      ! Persistent Chapman state — previous-timestep η at the wet side of each
      ! open edge. Per-cell `eta_old_<edge>(:)` arrays are declared for a
      ! future per-face adaptive Orlanski; currently unallocated (scalar
      ! Chapman on edge-mean η is used instead).
      real(wp), allocatable :: eta_old_west(:), eta_old_east(:)
      real(wp), allocatable :: eta_old_south(:), eta_old_north(:)

      ! Scalar Chapman state — one persistent η per edge, used by the
      ! barotropic substep's OBC_CHAPMAN dispatch when no per-cell array is
      ! supplied.  Updated at the end of each barotropic-substep call so the
      ! next call sees the radiation history.
      real(wp) :: eta_old_chapman_w = 0.0_wp
      real(wp) :: eta_old_chapman_e = 0.0_wp
      real(wp) :: eta_old_chapman_s = 0.0_wp
      real(wp) :: eta_old_chapman_n = 0.0_wp

      ! Open-edge tracer reservoirs. Allocated only when res_lscale_out > 0
      ! or res_lscale_in > 0 on an open-ish edge. Shape per edge:
      !   tres_west/east : (ny_total, nz_ml, n_tracers)
      !   tres_south/north : (nx_total, nz_ml, n_tracers)
      real(wp), allocatable :: tres_west(:, :, :)   !! West  reservoir concentration
      real(wp), allocatable :: tres_east(:, :, :)   !! East  reservoir concentration
      real(wp), allocatable :: tres_south(:, :, :)  !! South reservoir concentration
      real(wp), allocatable :: tres_north(:, :, :)  !! North reservoir concentration
      ! Cached length-scale knobs (m, default 0 = feature disabled).
      real(wp) :: res_lscale_out = 0.0_wp
         !! Outflow reservoir length scale (m).  0 ⇒ instantaneous outflow.
      real(wp) :: res_lscale_in = 0.0_wp
         !! Inflow reservoir length scale (m).  0 ⇒ instantaneous inflow.

      ! Per-layer Orlanski radiation state. Allocated when
      ! radiation_scheme == "orlanski" AND the edge radiates.
      ! Shapes: rx_*/u_prev_* west/east (ny_total, nz_ml), south/north
      ! (nx_total, nz_ml).
      ! rx:     running-mean nondimensional phase speed (grid cells / step).
      !         Not restart-registered ⇒ restarts cold.
      ! u_prev: first-interior-face normal velocity from the previous call
      !         (seeded from u_new on first call ⇒ rx = 0 cold start).
      real(wp), allocatable :: rx_west(:, :)     !! Running-mean rx, west  edge.
      real(wp), allocatable :: rx_east(:, :)     !! Running-mean rx, east  edge.
      real(wp), allocatable :: rx_south(:, :)    !! Running-mean rx, south edge.
      real(wp), allocatable :: rx_north(:, :)    !! Running-mean rx, north edge.
      real(wp), allocatable :: u_prev_west(:, :)  !! Prev-call u at first interior face, west.
      real(wp), allocatable :: u_prev_east(:, :)  !! Prev-call u at first interior face, east.
      real(wp), allocatable :: u_prev_south(:, :)  !! Prev-call v at first interior face, south.
      real(wp), allocatable :: u_prev_north(:, :)  !! Prev-call v at first interior face, north.

      ! Cached Orlanski / nudging / Flather knobs (set from config in configure_ocean_bc).
      integer  :: radiation_scheme = 0
         !! 0 = anomaly (default); 1 = orlanski.
      real(wp) :: orlanski_rx_max = 10.0_wp
         !! Upper clamp on the nondimensional phase speed (Orlanski 1976).
      real(wp) :: orlanski_gamma = 1.0_wp
         !! Running-mean weight.  1.0 = no running mean (instant rx).
      real(wp) :: nudge_tau_in = 0.0_wp
         !! Inflow nudging timescale (s, Marchesiello et al. 2001).  0 = off.
      real(wp) :: nudge_tau_out = 0.0_wp
         !! Outflow nudging timescale (s).  0 = off.

      ! Cached Flather-form knob.
      logical  :: use_full_flather = .false.
         !! .false. = legacy (default, bit-identical); .true. = full Flather
         !! (Flather 1976 half-characteristic form with exterior velocity).
      real(wp) :: ext_u_west = 0.0_wp  !! Exterior barotropic u, west  (m/s).
      real(wp) :: ext_u_east = 0.0_wp  !! Exterior barotropic u, east  (m/s).
      real(wp) :: ext_v_south = 0.0_wp  !! Exterior barotropic v, south (m/s).
      real(wp) :: ext_v_north = 0.0_wp  !! Exterior barotropic v, north (m/s).

      ! Per-step boundary data populated by the data source before each outer
      ! step. Shapes (ny[+1], nz_ml) or (nx[+1], nz_ml) per edge. Empty when
      ! all edges are WALL (the data source sizes them only when needed).
      real(wp), allocatable :: data_u_west(:, :), data_u_east(:, :)
      real(wp), allocatable :: data_v_south(:, :), data_v_north(:, :)
      real(wp), allocatable :: data_eta_west(:), data_eta_east(:)
      real(wp), allocatable :: data_eta_south(:), data_eta_north(:)
      real(wp), allocatable :: data_tracer_west(:, :, :), data_tracer_east(:, :, :)
      real(wp), allocatable :: data_tracer_south(:, :, :), data_tracer_north(:, :, :)
         !! (ny|nx, nz_ml, n_tracers).

      ! Cached grid metadata (avoids passing the grid through every helper).
      integer :: nx_total = 0, ny_total = 0
      integer :: nx_phys = 0, ny_phys = 0
      integer :: nghost = 0
      integer :: nz_ml = 0
      integer :: n_tracers = 0
   contains
      procedure, non_overridable :: bytes => ocean_bc_state_bytes
   end type ocean_bc_state_t

contains

   subroutine ocean_bc_state_init(this, grid, nz_ml, n_tracers)
      !! Cache grid extents and derive periodic flags. Data buffers stay
      !! unallocated until a data source asks for them. Call
      !! `ocean_bc_validate_periodic` after setting per-edge tags if any edge
      !! is OBC_PERIODIC; init itself only derives the convenience flags.
      type(ocean_bc_state_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz_ml
      integer, intent(in), optional :: n_tracers

      this%nx_total = grid%nx_total
      this%ny_total = grid%ny_total
      this%nx_phys = grid%nx_phys
      this%ny_phys = grid%ny_phys
      this%nghost = grid%nghost
      this%nz_ml = nz_ml
      this%n_tracers = 0
      if (present(n_tracers)) this%n_tracers = n_tracers
      ! Reset edge tags to WALL so re-init after destroy doesn't inherit
      ! stale tags.
      this%west%bc_type = OBC_WALL
      this%east%bc_type = OBC_WALL
      this%south%bc_type = OBC_WALL
      this%north%bc_type = OBC_WALL
      ! Pre-allocate per-edge clamped_tracer arrays (zero) so callers can set
      ! Dirichlet values without allocating; only consulted on OBC_CLAMPED.
      if (this%n_tracers > 0) then
         allocate (this%west%clamped_tracer(this%n_tracers), source=0.0_wp)
         allocate (this%east%clamped_tracer(this%n_tracers), source=0.0_wp)
         allocate (this%south%clamped_tracer(this%n_tracers), source=0.0_wp)
         allocate (this%north%clamped_tracer(this%n_tracers), source=0.0_wp)
      end if
      ! Derive convenience periodic flags.  Tags may still be at the
      ! default OBC_WALL here if the caller sets them after init — the
      ! caller must call ocean_bc_validate_periodic after finalising tags.
      this%periodic_x = (this%west%bc_type == OBC_PERIODIC .and. &
                         this%east%bc_type == OBC_PERIODIC)
      this%periodic_y = (this%south%bc_type == OBC_PERIODIC .and. &
                         this%north%bc_type == OBC_PERIODIC)
      this%north_fold = (this%north%bc_type == OBC_TRIPOLAR_FOLD)
      this%is_init = .true.
   end subroutine ocean_bc_state_init

   subroutine ocean_bc_validate_periodic(this, ierr)
      !! Validate periodic pairing + ghost-width + sponge incompatibility.
      !! Call after all per-edge tags are set and after ocean_bc_state_init.
      !! Derives `periodic_x` / `periodic_y` from the final tags and
      !! stops with a diagnostic message if any rule is violated.
      !!
      !! Rules (design §1.3):
      !!   (a) west periodic ⟺ east periodic (must be paired).
      !!   (b) south periodic ⟺ north periodic (must be paired).
      !!   (c) periodic requires nghost >= 3 (PPM 5-point + biharmonic).
      !!   (d) a periodic edge cannot be paired with OBC_SPONGE on any
      !!       edge in the same axis-pair.
      type(ocean_bc_state_t), intent(inout) :: this
      integer, intent(out), optional :: ierr
         !! Non-zero on a periodic-BC pairing/ghost-width violation when
         !! present; absent behaves as today (`error stop`).

      logical :: w_per, e_per, s_per, n_per

      w_per = (this%west%bc_type == OBC_PERIODIC)
      e_per = (this%east%bc_type == OBC_PERIODIC)
      s_per = (this%south%bc_type == OBC_PERIODIC)
      n_per = (this%north%bc_type == OBC_PERIODIC)

      ! Rule (a): zonal pairing
      if (w_per .neqv. e_per) then
         call logger%error("ocean_bc_validate_periodic: west periodic requires east "// &
                           "periodic (must be paired)")
         if (present(ierr)) then
            ierr = OCEAN_STATUS_ERR_SETUP
            return
         end if

         error stop "ocean_bc_validate_periodic: west periodic requires east periodic (must be paired)"
      end if
      ! Rule (b): meridional pairing
      if (s_per .neqv. n_per) then
         call logger%error("ocean_bc_validate_periodic: south periodic requires north "// &
                           "periodic (must be paired)")
         if (present(ierr)) then
            ierr = OCEAN_STATUS_ERR_SETUP
            return
         end if

         error stop "ocean_bc_validate_periodic: south periodic requires north periodic (must be paired)"
      end if
      ! Rule (c): ghost width
      if ((w_per .or. s_per) .and. this%nghost < 3) then
         call logger%error("ocean_bc_validate_periodic: periodic BC requires nghost "// &
                           ">= 3 (PPM + biharmonic stencil depth)")
         if (present(ierr)) then
            ierr = OCEAN_STATUS_ERR_SETUP
            return
         end if

         error stop "ocean_bc_validate_periodic: periodic BC requires nghost >= 3 (PPM + biharmonic stencil depth)"
      end if
      ! Rule (d) "a periodic edge cannot be a sponge edge" is enforced
      ! structurally: each edge carries exactly one tag, and rules (a)/(b)
      ! force the partner edge of a periodic edge to be periodic too — so
      ! no edge on a periodic axis can carry OBC_SPONGE.  Cross-axis
      ! combinations (e.g. periodic-x with sponge bands at the y-walls —
      ! the reentrant-channel configuration) are deliberately allowed.

      ! Update derived flags now that tags are finalised.
      this%periodic_x = w_per
      this%periodic_y = s_per
      if (present(ierr)) ierr = OCEAN_STATUS_OK
   end subroutine ocean_bc_validate_periodic

   subroutine ocean_bc_validate_fold(this, ierr)
      !! Validate the tripolar north-fold tag.  Call after all per-edge
      !! tags are set and after `ocean_bc_state_init`.  Refreshes
      !! `north_fold` and stops with a diagnostic if any rule fails.
      !!
      !! Rules (design Appendix A):
      !!   (a) OBC_TRIPOLAR_FOLD is accepted ONLY on the NORTH edge —
      !!       west/east/south carrying it is a config error.
      !!   (b) the fold requires periodic west AND east (the fold reads
      !!       already-cyclically-wrapped corner columns).
      !!   (c) the fold requires nghost >= 3 (PPM + biharmonic stencil
      !!       depth — same as periodic).
      !! The reverse rule (tripolar grid REQUIRES the fold tag) is
      !! checked at the configure-metrics site, which knows grid_config.
      type(ocean_bc_state_t), intent(inout) :: this
      integer, intent(out), optional :: ierr
         !! Non-zero on a tripolar-fold configuration violation when
         !! present; absent behaves as today (`error stop`).

      logical :: n_fold

      n_fold = (this%north%bc_type == OBC_TRIPOLAR_FOLD)

      ! Rule (a): only the north edge may carry the fold tag.
      if (this%west%bc_type == OBC_TRIPOLAR_FOLD .or. &
          this%east%bc_type == OBC_TRIPOLAR_FOLD .or. &
          this%south%bc_type == OBC_TRIPOLAR_FOLD) then
         call logger%error("ocean_bc_validate_fold: tripolar_fold is accepted "// &
                           "only on the north edge")
         if (present(ierr)) then
            ierr = OCEAN_STATUS_ERR_SETUP
            return
         end if

         error stop "ocean_bc_validate_fold: tripolar_fold is accepted only on the north edge"
      end if

      if (n_fold) then
         ! Rule (b): periodic west+east mandatory.
         if (this%west%bc_type /= OBC_PERIODIC .or. &
             this%east%bc_type /= OBC_PERIODIC) then
            call logger%error("ocean_bc_validate_fold: north='tripolar_fold' "// &
                              "requires periodic west+east edges")
            if (present(ierr)) then
               ierr = OCEAN_STATUS_ERR_SETUP
               return
            end if

            error stop "ocean_bc_validate_fold: north='tripolar_fold' requires periodic west+east edges"
         end if
         ! Rule (c): ghost width.
         if (this%nghost < 3) then
            call logger%error("ocean_bc_validate_fold: tripolar_fold requires "// &
                              "nghost >= 3 (PPM + biharmonic stencil depth)")
            if (present(ierr)) then
               ierr = OCEAN_STATUS_ERR_SETUP
               return
            end if

            error stop "ocean_bc_validate_fold: tripolar_fold requires nghost >= 3 (PPM + biharmonic stencil depth)"
         end if
      end if

      this%north_fold = n_fold
      if (present(ierr)) ierr = OCEAN_STATUS_OK
   end subroutine ocean_bc_validate_fold

   subroutine ocean_bc_state_destroy(this)
      type(ocean_bc_state_t), intent(inout) :: this
      this%is_init = .false.
      if (allocated(this%west%clamped_tracer)) deallocate (this%west%clamped_tracer)
      if (allocated(this%east%clamped_tracer)) deallocate (this%east%clamped_tracer)
      if (allocated(this%south%clamped_tracer)) deallocate (this%south%clamped_tracer)
      if (allocated(this%north%clamped_tracer)) deallocate (this%north%clamped_tracer)
      if (allocated(this%eta_old_west)) deallocate (this%eta_old_west)
      if (allocated(this%eta_old_east)) deallocate (this%eta_old_east)
      if (allocated(this%eta_old_south)) deallocate (this%eta_old_south)
      if (allocated(this%eta_old_north)) deallocate (this%eta_old_north)
      ! Reservoir arrays (§1, v2).
      if (allocated(this%tres_west)) deallocate (this%tres_west)
      if (allocated(this%tres_east)) deallocate (this%tres_east)
      if (allocated(this%tres_south)) deallocate (this%tres_south)
      if (allocated(this%tres_north)) deallocate (this%tres_north)
      ! Orlanski radiation arrays (§2, v2).
      if (allocated(this%rx_west)) deallocate (this%rx_west)
      if (allocated(this%rx_east)) deallocate (this%rx_east)
      if (allocated(this%rx_south)) deallocate (this%rx_south)
      if (allocated(this%rx_north)) deallocate (this%rx_north)
      if (allocated(this%u_prev_west)) deallocate (this%u_prev_west)
      if (allocated(this%u_prev_east)) deallocate (this%u_prev_east)
      if (allocated(this%u_prev_south)) deallocate (this%u_prev_south)
      if (allocated(this%u_prev_north)) deallocate (this%u_prev_north)
      if (allocated(this%data_u_west)) deallocate (this%data_u_west)
      if (allocated(this%data_u_east)) deallocate (this%data_u_east)
      if (allocated(this%data_v_south)) deallocate (this%data_v_south)
      if (allocated(this%data_v_north)) deallocate (this%data_v_north)
      if (allocated(this%data_eta_west)) deallocate (this%data_eta_west)
      if (allocated(this%data_eta_east)) deallocate (this%data_eta_east)
      if (allocated(this%data_eta_south)) deallocate (this%data_eta_south)
      if (allocated(this%data_eta_north)) deallocate (this%data_eta_north)
      if (allocated(this%data_tracer_west)) deallocate (this%data_tracer_west)
      if (allocated(this%data_tracer_east)) deallocate (this%data_tracer_east)
      if (allocated(this%data_tracer_south)) deallocate (this%data_tracer_south)
      if (allocated(this%data_tracer_north)) deallocate (this%data_tracer_north)
   end subroutine ocean_bc_state_destroy

   subroutine ocean_bc_state_enter_data(this)
      !! GPU mapping for `ocean_bc_state_t`.
      !!
      !! Parent-first rule (design §1, GPU storage): `copyin(this)` MUST precede
      !! the component copies so the device descriptor for `this` is live before
      !! the component attach.  Reverse applies on exit (components first, parent
      !! last).  Missing the parent copyin causes UVM-page-fault per DC launch.
      !!
      !! The reservoir arrays (§1) and Orlanski radiation arrays (§2) are mapped;
      !! the data_* buffers remain host-only (not consumed by device kernels).
      type(ocean_bc_state_t), intent(inout) :: this
      if (.not. this%is_init) return
      !$acc enter data copyin(this)
      if (allocated(this%tres_west)) then
         !$acc enter data copyin(this%tres_west)
      end if
      if (allocated(this%tres_east)) then
         !$acc enter data copyin(this%tres_east)
      end if
      if (allocated(this%tres_south)) then
         !$acc enter data copyin(this%tres_south)
      end if
      if (allocated(this%tres_north)) then
         !$acc enter data copyin(this%tres_north)
      end if
      ! Orlanski rx and u_prev arrays (§2, v2).
      if (allocated(this%rx_west)) then
         !$acc enter data copyin(this%rx_west)
      end if
      if (allocated(this%rx_east)) then
         !$acc enter data copyin(this%rx_east)
      end if
      if (allocated(this%rx_south)) then
         !$acc enter data copyin(this%rx_south)
      end if
      if (allocated(this%rx_north)) then
         !$acc enter data copyin(this%rx_north)
      end if
      if (allocated(this%u_prev_west)) then
         !$acc enter data copyin(this%u_prev_west)
      end if
      if (allocated(this%u_prev_east)) then
         !$acc enter data copyin(this%u_prev_east)
      end if
      if (allocated(this%u_prev_south)) then
         !$acc enter data copyin(this%u_prev_south)
      end if
      if (allocated(this%u_prev_north)) then
         !$acc enter data copyin(this%u_prev_north)
      end if
   end subroutine ocean_bc_state_enter_data

   subroutine ocean_bc_state_exit_data(this)
      !! GPU unmapping — components first, parent last (reverse of enter_data).
      type(ocean_bc_state_t), intent(inout) :: this
      if (.not. this%is_init) return
      ! Orlanski arrays (§2, v2) — reverse allocation order.
      if (allocated(this%u_prev_north)) then
         !$acc exit data delete(this%u_prev_north)
      end if
      if (allocated(this%u_prev_south)) then
         !$acc exit data delete(this%u_prev_south)
      end if
      if (allocated(this%u_prev_east)) then
         !$acc exit data delete(this%u_prev_east)
      end if
      if (allocated(this%u_prev_west)) then
         !$acc exit data delete(this%u_prev_west)
      end if
      if (allocated(this%rx_north)) then
         !$acc exit data delete(this%rx_north)
      end if
      if (allocated(this%rx_south)) then
         !$acc exit data delete(this%rx_south)
      end if
      if (allocated(this%rx_east)) then
         !$acc exit data delete(this%rx_east)
      end if
      if (allocated(this%rx_west)) then
         !$acc exit data delete(this%rx_west)
      end if
      ! Reservoir arrays (§1, v2).
      if (allocated(this%tres_north)) then
         !$acc exit data delete(this%tres_north)
      end if
      if (allocated(this%tres_south)) then
         !$acc exit data delete(this%tres_south)
      end if
      if (allocated(this%tres_east)) then
         !$acc exit data delete(this%tres_east)
      end if
      if (allocated(this%tres_west)) then
         !$acc exit data delete(this%tres_west)
      end if
      !$acc exit data delete(this)
   end subroutine ocean_bc_state_exit_data

   pure subroutine ocean_bc_state_set_edges(this, has_west, has_east, has_south, has_north)
      !! Set the physical-domain-edge flags from a decomposition descriptor.
      !! Called once by the driver after configure_ocean_bc so kernels can
      !! gate wall / BC / periodic closures on physical edges (a subdomain
      !! seam is never a wall).  Default .true. keeps single-rank bit-identity.
      type(ocean_bc_state_t), intent(inout) :: this
      logical, intent(in) :: has_west
         !! True when the west edge is a physical domain edge, false at an MPI seam.
      logical, intent(in) :: has_east
         !! True when the east edge is a physical domain edge, false at an MPI seam.
      logical, intent(in) :: has_south
         !! True when the south edge is a physical domain edge, false at an MPI seam.
      logical, intent(in) :: has_north
         !! True when the north edge is a physical domain edge, false at an MPI seam.
      this%has_west = has_west
      this%has_east = has_east
      this%has_south = has_south
      this%has_north = has_north
   end subroutine ocean_bc_state_set_edges

   subroutine ocean_bc_state_set_topology(this, periodic_x, periodic_y, ierr)
      !! Pre-create GRID TOPOLOGY injection (Python runtime API plan,
      !! P2.5): force per-dimension periodicity the Oceananigans way
      !! (`docs/ocean_python_api_plan.md` S5b) — "the grid owns
      !! periodicity", not the per-edge `&ocean_bc_nml` tags. Sets
      !! `periodic_x`/`periodic_y` directly and back-fills the edge tags
      !! on every axis the caller marks periodic (both edges together, so
      !! a west/east — or south/north — mismatch is structurally
      !! unrepresentable through this entry point, unlike the namelist
      !! path which needs `ocean_bc_validate_periodic` to catch one). An
      !! axis the caller does NOT mark periodic is left untouched: its
      !! edge tags keep whatever physical BC `configure_ocean_bc` already
      !! derived from `&ocean_bc_nml` — periodicity is a GRID property,
      !! but the wall/open/clamped/... physics for a Bounded dimension
      !! stays the namelist's job.
      !!
      !! Call AFTER `configure_ocean_bc` (whose namelist-derived tags this
      !! may override) and BEFORE anything that reads `periodic_x`/`_y` —
      !! the init-time periodic ghost wrap, `ocean_halo_init`,
      !! `configure_ocean_land_mask` (`engine_setup`'s ordering).
      type(ocean_bc_state_t), intent(inout) :: this
      logical, intent(in) :: periodic_x
      logical, intent(in) :: periodic_y
      integer, intent(out), optional :: ierr
         !! Non-zero (`OCEAN_STATUS_ERR_SETUP`) when a requested periodic
         !! axis violates the `nghost >= 3` PPM/biharmonic stencil-depth
         !! requirement (`ocean_bc_validate_periodic`'s rule (c), mirrored
         !! here since this entry point bypasses that routine) when
         !! present; absent behaves as today (`error stop`).

      if ((periodic_x .or. periodic_y) .and. this%nghost < 3) then
         call logger%error("ocean_bc_state_set_topology: periodic topology requires "// &
                           "nghost >= 3 (PPM + biharmonic stencil depth)")
         if (present(ierr)) then
            ierr = OCEAN_STATUS_ERR_SETUP
            return
         end if
         error stop "ocean_bc_state_set_topology: periodic topology requires nghost >= 3"
      end if

      this%periodic_x = periodic_x
      this%periodic_y = periodic_y
      if (periodic_x) then
         this%west%bc_type = OBC_PERIODIC
         this%east%bc_type = OBC_PERIODIC
      end if
      if (periodic_y) then
         this%south%bc_type = OBC_PERIODIC
         this%north%bc_type = OBC_PERIODIC
      end if
      if (present(ierr)) ierr = OCEAN_STATUS_OK
   end subroutine ocean_bc_state_set_topology

   pure function ocean_bc_type_from_string(name) result(bc_type)
      !! Convert a config-namelist edge string → integer OBC tag.
      !! Case-INSENSITIVE (`to_lower`), so "OPEN"/"Open"/"open" all parse
      !! to `OBC_OPEN`.  An unrecognised name returns `OBC_INVALID`
      !! (PR-6 fail-loud): a typo must NOT silently close the boundary to
      !! a wall.  `validate_config` rejects `OBC_INVALID` (naming the
      !! edge) before `configure_ocean_bc` consumes any parse result, so
      !! no production caller ever sees the sentinel at a live edge.
      character(len=*), intent(in) :: name
      integer :: bc_type
      select case (to_lower(trim(name)))
      case ("wall")
         bc_type = OBC_WALL
      case ("open")
         bc_type = OBC_OPEN
      case ("tidal")
         bc_type = OBC_TIDAL
      case ("nested")
         bc_type = OBC_NESTED
      case ("inflow")
         bc_type = OBC_INFLOW
      case ("discharge")
         bc_type = OBC_DISCHARGE
      case ("clamped")
         bc_type = OBC_CLAMPED
      case ("sponge")
         bc_type = OBC_SPONGE
      case ("chapman")
         bc_type = OBC_CHAPMAN
      case ("periodic")
         bc_type = OBC_PERIODIC
      case ("tripolar_fold")
         bc_type = OBC_TRIPOLAR_FOLD
      case default
         bc_type = OBC_INVALID
      end select
   end function ocean_bc_type_from_string

   pure logical function ocean_bc_has_tracer_open_edge(bc) result(res)
      !! `.true.` if ANY outer edge is a tracer-open boundary — the set
      !! (OPEN / TIDAL / CHAPMAN / CLAMPED / NESTED) across which a lateral
      !! tracer flux (e.g. Redi neutral diffusion) can leave the domain
      !! without being mirrored into the console `out` budget.  Used to fall
      !! the closed Salt/Heat `Error` back to raw drift.  Reads the PARSED
      !! per-edge `bc_type` (not the raw config string), so it is robust to
      !! tag case and covers every open-type edge, not just literal "open".
      type(ocean_bc_state_t), intent(in) :: bc
      res = edge_is_tracer_open(bc%west%bc_type) &
            .or. edge_is_tracer_open(bc%east%bc_type) &
            .or. edge_is_tracer_open(bc%south%bc_type) &
            .or. edge_is_tracer_open(bc%north%bc_type)
   end function ocean_bc_has_tracer_open_edge

   pure logical function edge_is_tracer_open(bc_type) result(res)
      !! One-edge tracer-open test (module-internal helper).
      integer, intent(in) :: bc_type
      res = (bc_type == OBC_OPEN .or. bc_type == OBC_TIDAL &
             .or. bc_type == OBC_CHAPMAN .or. bc_type == OBC_CLAMPED &
             .or. bc_type == OBC_NESTED)
   end function edge_is_tracer_open

   pure function obc_match_constituent(omega) result(ic)
      !! Resolve an OBC edge constituent's angular frequency `omega` (rad/s)
      !! to the tide catalog index (`rdb_ocean_tide_astro::TIDE_OMEGA`) whose
      !! frequency matches within the relative tolerance `OBC_TIDE_MATCH_TOL`.
      !! Returns 0 when no catalog entry is within tolerance (unknown
      !! constituent) or when `omega <= 0` — the caller (OBC setup) converts a
      !! 0 to a fail-loud `error stop`, keeping this function `pure`.
      real(wp), intent(in) :: omega
      integer :: ic
      integer :: c, best_c
      real(wp) :: best_rel, rel
      ic = 0
      if (omega <= 0.0_wp) return
      best_c = 0
      best_rel = huge(1.0_wp)
      do c = 1, TIDES_CATALOG_SIZE
         rel = abs(omega - TIDE_OMEGA(c))/omega
         if (rel < best_rel) then
            best_rel = rel
            best_c = c
         end if
      end do
      if (best_rel <= OBC_TIDE_MATCH_TOL) ic = best_c
   end function obc_match_constituent

   pure subroutine obc_tide_nodal_fill(face, f_all, u_all, v_all, ierr)
      !! Bake the nodal/astronomical correction into one edge's per-constituent
      !! `tidal_fnodal` / `tidal_arg`.  For each of `face%n_tidal_constituents`,
      !! resolve the constituent by frequency (`obc_match_constituent`) and set
      !! `tidal_fnodal(nc) = f_all(ic)`, `tidal_arg(nc) = v_all(ic) + u_all(ic)`.
      !! `f_all` / `u_all` come from `nodal_fu`, `v_all` from
      !! `equilibrium_arguments`, all sized `TIDES_CATALOG_SIZE`.  On an
      !! unmatched constituent it leaves that entry untouched and returns
      !! `ierr = nc` (the 1-based edge slot that failed) so the caller can fail
      !! loud; `ierr = 0` on success.  `pure` — no logging / no `error stop`.
      type(ocean_bc_face_tag_t), intent(inout) :: face
      real(wp), intent(in) :: f_all(TIDES_CATALOG_SIZE)
      real(wp), intent(in) :: u_all(TIDES_CATALOG_SIZE)
      real(wp), intent(in) :: v_all(TIDES_CATALOG_SIZE)
      integer, intent(out) :: ierr
      integer :: nc, ic
      ierr = 0
      do nc = 1, face%n_tidal_constituents
         ic = obc_match_constituent(face%tidal_omega(nc))
         if (ic == 0) then
            ierr = nc
            return
         end if
         face%tidal_fnodal(nc) = f_all(ic)
         face%tidal_arg(nc) = v_all(ic) + u_all(ic)
      end do
   end subroutine obc_tide_nodal_fill

   pure function ocean_bc_state_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the boundary state slot (0 when
      !! unallocated).
      class(ocean_bc_state_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = arr_bytes(this%eta_old_west) &
               + arr_bytes(this%eta_old_east) &
               + arr_bytes(this%eta_old_south) &
               + arr_bytes(this%eta_old_north) &
               + arr_bytes(this%tres_west) &
               + arr_bytes(this%tres_east) &
               + arr_bytes(this%tres_south) &
               + arr_bytes(this%tres_north) &
               + arr_bytes(this%rx_west) &
               + arr_bytes(this%rx_east) &
               + arr_bytes(this%rx_south) &
               + arr_bytes(this%rx_north) &
               + arr_bytes(this%u_prev_west) &
               + arr_bytes(this%u_prev_east) &
               + arr_bytes(this%u_prev_south) &
               + arr_bytes(this%u_prev_north) &
               + arr_bytes(this%data_u_west) &
               + arr_bytes(this%data_u_east) &
               + arr_bytes(this%data_v_south) &
               + arr_bytes(this%data_v_north) &
               + arr_bytes(this%data_eta_west) &
               + arr_bytes(this%data_eta_east) &
               + arr_bytes(this%data_eta_south) &
               + arr_bytes(this%data_eta_north) &
               + arr_bytes(this%data_tracer_west) &
               + arr_bytes(this%data_tracer_east) &
               + arr_bytes(this%data_tracer_south) &
               + arr_bytes(this%data_tracer_north)
   end function ocean_bc_state_bytes

end module rdb_ocean_boundary_types
