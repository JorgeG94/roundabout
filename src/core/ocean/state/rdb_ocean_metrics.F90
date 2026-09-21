!! Orthogonal curvilinear horizontal metrics for the ocean dyn-core.
module rdb_ocean_metrics
   !! `ocean_metrics_t` — the full 2D metric arrays the curvilinear ocean
   !! dyn-core reads.
   !! Coordinates are GENERATORS that fill these arrays; kernels consume
   !! the metrics only and never recompute `1/dx` or `dx*dy` themselves.
   !!
   !! Storage convention (mirrors the C-grid prognostic sizing exactly,
   !! `nx = grid%nx_total`, `ny = grid%ny_total`):
   !!   * T  (cell centre)  arrays: `(nx,   ny)`     — like `h_layer`.
   !!   * Cu (east  u-face) arrays: `(nx+1, ny)`     — like `u_face_x`.
   !!   * Cv (north v-face) arrays: `(nx,   ny+1)`   — like `v_face_y`.
   !!   * Bu (NE  corner)   arrays: `(nx+1, ny+1)`   — like `f_corner`.
   !! All staggers are filled INCLUDING ghost rows/columns — the metric
   !! formulae extend naturally and unfilled ghosts are a known
   !! EOS-blowup class of bug (formula bathymetry ghost-fill gotcha).
   !!
   !! Inverses + the hvisc ratio bundle are single-sourced: computed ONCE
   !! in `metrics_finalize` from the arrays a generator wrote, via the
   !! Adcroft reciprocal (`1/x` with `0 -> 0`).  Kernels never recompute
   !! them: a round-trip `1/(1/dx)` mismatch breaks the exact telescoping
   !! that continuity relies on (D4).  `areaT` is load-bearing and
   !! `dx*dy` is dead: on supergrid / tripolar grids `areaT /= dxT*dyT`
   !! (D5), so areas are stored independently.
   !!
   !! References: MOM6 grid architecture (`MOM_dyn_horgrid` / `MOM_grid`
   !! metric vocabulary, studied 2026-06-11); Adcroft reciprocal.  This
   !! is an independent implementation — no source ported.
   use rdb_constants, only: wp
   use rdb_grid, only: hgrid_t
#ifndef RDB_NO_NETCDF
   use rdb_io_netcdf, only: nc_check, nc_open_read, nc_close, &
                            nc_get_dim_len, nc_get_varid, nc_get_var_2d
#endif
   use rdb_ocean_bipolar, only: bipolar_corner_latlon
   use rdb_ocean_fold, only: fold_north_centre, fold_north_corner
   ! NOTE (O3): the multi-rank wet_mask seam exchange is done by the CALLER
   ! (rdb_ocean_setup::configure_ocean_land_mask) on the wet_mask array BEFORE
   ! metrics_apply_land_mask runs — NOT here.  Importing the comm-layer
   ! rdb_ocean_halo into this low-level metrics leaf creates an NVFORTRAN USE
   ! cycle (rdb_config -> rdb_ocean_lateral_mix -> rdb_ocean_metrics ->
   ! rdb_ocean_halo -> ... -> back).  Keeping the exchange in the caller (which
   ! is downstream of rdb_config) avoids the inversion.
   use pic_logger, only: logger => global_logger
   use rdb_error_ring, only: fail
   use pic_strings, only: to_string
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_mem_report, only: arr_bytes
   use rdb_ocean_status, only: OCEAN_STATUS_OK, OCEAN_STATUS_ERR_IO
   implicit none
   private

   public :: ocean_metrics_t
   public :: metrics_finalize
   public :: metrics_apply_land_mask
   public :: metrics_fill_cartesian
   public :: metrics_fill_spherical
   public :: metrics_fill_from_supergrid
   public :: metrics_assemble_from_supergrid_arrays
   public :: metrics_fill_tripolar
   public :: metrics_fill_coriolis
   public :: metrics_porous_alloc
   public :: metrics_closed_faces_alloc
   public :: adcroft_recip
   public :: GRID_CONFIG_CARTESIAN, GRID_CONFIG_SPHERICAL, GRID_CONFIG_SUPERGRID
   public :: GRID_CONFIG_TRIPOLAR
   public :: CORIOLIS_SCHEME_BETA_PLANE, CORIOLIS_SCHEME_PLANETARY
   public :: parse_grid_config, parse_coriolis_scheme

   ! ---- Grid-config enum (mirrors `&ocean_grid_nml grid_config`) ----
   integer, parameter :: GRID_CONFIG_CARTESIAN = 0
      !! Uniform Cartesian: every metric constant (bit-identity gate).
   integer, parameter :: GRID_CONFIG_SPHERICAL = 1
      !! Spherical lon-lat sector (analytic-derivative form).
   integer, parameter :: GRID_CONFIG_SUPERGRID = 2
      !! MOM6 supergrid (mosaic) NetCDF reader (v1 stub).
   integer, parameter :: GRID_CONFIG_TRIPOLAR = 3
      !! Analytic TRIPOLAR (Murray 1996): lon-lat below `phi_join`,
      !! bipolar Arctic cap above.  See `metrics_fill_tripolar`.

   ! ---- Coriolis-scheme enum (mirrors `&ocean_grid_nml coriolis_scheme`) ----
   integer, parameter :: CORIOLIS_SCHEME_BETA_PLANE = 0
      !! `f = f_0 + beta*(y - y_ref)`, y from the Cartesian coordinate.
   integer, parameter :: CORIOLIS_SCHEME_PLANETARY = 1
      !! `f = 2*omega*sin(geolat)` at the respective stagger.

   real(wp), parameter :: DEG2RAD = 3.14159265358979323846_wp/180.0_wp
      !! Degrees -> radians.

   type :: ocean_metrics_t
      logical :: is_init = .false.
         !! True between `init` and `destroy`.  Guard on this, never on
         !! `allocated(...)` (host pointer only; misses GPU mapping).

      ! ---- Lengths (m) ----
      real(wp), allocatable :: dxT(:, :), dyT(:, :)
         !! Cell-centre (T) zonal/meridional grid lengths (m), `(nx,ny)`.
      real(wp), allocatable :: dxCu(:, :), dyCu(:, :)
         !! u-face (Cu) lengths (m), `(nx+1,ny)`.
      real(wp), allocatable :: dxCv(:, :), dyCv(:, :)
         !! v-face (Cv) lengths (m), `(nx,ny+1)`.
      real(wp), allocatable :: dxBu(:, :), dyBu(:, :)
         !! Corner (Bu) lengths (m), `(nx+1,ny+1)`.

      ! ---- Topography-aware face widths (m) ----
      real(wp), allocatable :: dy_cu(:, :)
         !! Open zonal width of the u-face for transport (m), `(nx+1,ny)`.
         !! v1: filled = `dyCu` but a SEPARATE array, so the transport
         !! kernels read the right name once porous/partial cells arrive.
      real(wp), allocatable :: dx_cv(:, :)
         !! Open meridional width of the v-face for transport (m),
         !! `(nx,ny+1)`.  v1: filled = `dxCv`, separate array.

      real(wp), allocatable :: dy_cu_bt(:, :)
         !! Open zonal u-face width the BAROTROPIC substep transports on
         !! (m), `(nx+1,ny)`.  ALWAYS full size and byte-equal to `dy_cu`
         !! unless porous barriers are on, in which case the per-step
         !! refresh scales it by the COLUMN-INTEGRATED open fraction
         !! `(A(eta_top) - A(eta_bed)) / (eta_top - eta_bed)`, which is
         !! identically the THICKNESS-WEIGHTED MEAN of the per-layer
         !! fractions (both are the same integral of `w` over the column,
         !! so the identity is exact, not an approximation).  Without it
         !! the barotropic solve would be porous-blind and the layer
         !! renormalisation (which drives `sum_k flux_k = uhbt`) would
         !! hand the blocked transport straight back.
         !!
         !! NOT a claim of `BT_cont` parity.  MOM6's production barotropic
         !! face area is `sum_k (dy_Cu*por_k) * h_marginal_k * visc_rem_k`
         !! — weighted by the PPM MARGINAL thickness and by `visc_rem`,
         !! neither of which appears here; its `sum_k h_k*(dy_Cu*por_k)`
         !! form is the open-boundary-segment branch only, and its
         !! `set_local_BT_cont_types` carries no `por` at all.  What this
         !! array reproduces is the telescoping identity above, applied to
         !! the plain layer thicknesses.
      real(wp), allocatable :: dx_cv_bt(:, :)
         !! v-face twin, `(nx,ny+1)`.

      ! ---- Porous barriers (Adcroft 2013; see rdb_ocean_porous) ----
      logical :: use_porous = .false.
         !! Master switch (`&ocean_porous_nml enable`).  OFF ⇒ the
         !! `por_face_area_*` arrays stay at their `(1,1,1)` placeholder
         !! size and every transport kernel takes the un-narrowed
         !! `dy_cu` / `dx_cv` branch — byte-identical to a build without
         !! porous barriers.
      integer :: porous_eta_interp = 0
         !! Interface-at-velocity-point rule, a `POROUS_ETA_*` value
         !! (`rdb_ocean_porous`).  0 = MAX (the shallower interface).
      real(wp) :: porous_mask_depth = 0.0_wp
         !! Gate HEIGHT (m, positive up, `<= 0`): faces whose mean
         !! along-face height is at or above this stay fully open
         !! (MOM6 `PORBAR_MASKING_DEPTH`, sign-flipped to a height).
      real(wp), allocatable :: por_bed(:, :)
         !! Static snapshot of the bottom topographic HEIGHT at cell
         !! centres (m, positive up — i.e. `-barotropic%b`, which is the
         !! positive-down reference depth), `(nx,ny)` when `use_porous`,
         !! `(1,1)` otherwise.
         !! The porous curve works in ABSOLUTE heights, so the recompute
         !! needs the bed on the same datum as `por_d*`; keeping a copy
         !! here makes the kernel self-contained (no barotropic-state
         !! argument threaded through the dynamics).
      real(wp), allocatable :: por_dmin_u(:, :), por_dmax_u(:, :), por_davg_u(:, :)
         !! u-face along-face deepest / shallowest / mean topographic
         !! height (m, positive up), `(nx+1,ny)` when `use_porous`,
         !! `(1,1)` otherwise.  Static — filled once at setup.
      real(wp), allocatable :: por_dmin_v(:, :), por_dmax_v(:, :), por_davg_v(:, :)
         !! v-face twins, `(nx,ny+1)` when `use_porous`, `(1,1)` otherwise.
      real(wp), allocatable :: por_face_area_u(:, :, :)
         !! u-face layer-averaged OPEN-AREA fraction (nondim, `[0,1]`),
         !! `(nx+1,ny,nz)` when `use_porous`, `(1,1,1)` otherwise.
         !! Recomputed every RK2 stage (interface-height dependent) and
         !! MULTIPLIED into `dy_cu` by the transport kernels.
      real(wp), allocatable :: por_face_area_v(:, :, :)
         !! v-face twin, `(nx,ny+1,nz)` when `use_porous`, `(1,1,1)`
         !! otherwise.

      ! ---- Partial-step z-level face closure (VCOORD_Z_FIXED) ----
      logical :: use_closed_faces = .false.
         !! Master switch (`&vcoord_nml zfixed_closed_faces`), latched by
         !! `configure_ocean_closed_faces`.  OFF ⇒ `open_u`/`open_v` stay
         !! at their `(1,1,1)` placeholder size, no kernel branch is
         !! taken, byte-identical to a build without the feature.
      real(wp), allocatable :: open_u(:, :, :)
         !! u-face per-layer 0/1 OPEN mask, `(nx+1,ny,nz)` when
         !! `use_closed_faces`, `(1,1,1)` otherwise.  1 = the layer has
         !! water on BOTH sides of the face; 0 = it is an inert `z_fixed`
         !! filler on at least one side and the face is a z-LEVEL WALL for
         !! that layer (Adcroft, Hill & Marshall 1997; Losch 2008).
         !! STATIC — built once at configure by
         !! `ocean_vcoord_closed_face_masks` from the `z_fixed` target at
         !! `η = 0`, never refreshed (the bed and the draft are static and
         !! `η` is absorbed by the first live layer).
      real(wp), allocatable :: open_v(:, :, :)
         !! v-face twin, `(nx,ny+1,nz)` when `use_closed_faces`,
         !! `(1,1,1)` otherwise.
         !!
         !! ### THE COMPOSITION RULE (stated once, here)
         !!
         !! The three face gates are INDEPENDENT and compose by
         !! multiplication — none replaces another:
         !! ```
         !! dy_eff(I,j,k) = dy_cu(I,j) · por_face_area_u(I,j,k) · open_u(I,j,k)
         !! dx_eff(i,J,k) = dx_cv(i,J) · por_face_area_v(i,J,k) · open_v(i,J,k)
         !! ```
         !! `dy_cu`/`dx_cv` carry the 2-D LAND decision (metric zeroing in
         !! `metrics_apply_land_mask`); `por_face_area_*` narrows
         !! continuously for unresolved SUBGRID sills (Adcroft 2013); and
         !! `open_*` closes per LAYER for the resolved z-level staircase.
         !! Porous barriers and closed faces are therefore NOT mutually
         !! exclusive.
         !!
         !! Every consumer applies the two 3-D factors as SEPARATE,
         !! separately host-gated, INLINE `do concurrent` passes rather
         !! than pre-composing them into a third array.  Two reasons:
         !! a composed array would have to be recomputed whenever the
         !! porous fit is refreshed (per outer step) and so could not be
         !! static; and an inert host-gated branch that never names the
         !! array costs nothing, whereas handing a state array to an
         !! external helper pessimises every `do concurrent` in the
         !! calling routine even when the branch is not taken (CLAUDE.md,
         !! measured at +4.8 % for an inert porous pass).

      ! ---- Static ice-shelf cavity geometry (P5.1; see rdb_ocean_cavity) ----
      logical :: use_cavity = .false.
         !! Master switch (`&ocean_cavity_dyn_nml enable`), latched in
         !! `ocean_state_init_from_config` BEFORE `init` so the allocation
         !! gate below can read it.  OFF ⇒ `z_draft` / `cover_frac` /
         !! `p_ice_ref` stay at their `(1,1)` placeholder size, `bt_H_ref`
         !! latches the bed as it always did, and every path is
         !! byte-identical to a build without cavities.
      real(wp), allocatable :: z_draft(:, :)
         !! Prescribed STATIC ice-base depth (m, positive DOWN, `>= 0`),
         !! `(nx_total, ny_total)` INCLUDING ghosts when `use_cavity`,
         !! `(1,1)` otherwise.  Filled by the formula setters in
         !! `rdb_ocean_cavity` immediately after the bathymetry, then
         !! carried through the SAME periodic/fold re-wrap + halo sequence
         !! `barotropic%b` gets (ordering is load-bearing: the draft must
         !! exist before the wet mask is seeded from `b - z_draft`).
         !! `z_draft = 0` is open ocean — including beyond the calving
         !! front.
      real(wp), allocatable :: cover_frac(:, :)
         !! Ice-covered area fraction (nondimensional, `[0,1]`), same
         !! shape + gating as `z_draft`.  v1 is BINARY, `merge(1, 0,
         !! z_draft > 0)`; an area-blended calving front belongs to the
         !! melt work.  Allocated alongside `z_draft` so the thermodynamic
         !! slice does not have to re-open this lifecycle.
      real(wp), allocatable :: p_ice_ref(:, :)
         !! Boussinesq-isostatic (flotation) ice load `rho_ref*GRAVITY*
         !! z_draft` (Pa, `>= 0`), same shape + gating as `z_draft`.  Built
         !! ONCE at configure from the SAME product the FV_MOM6 surface BC
         !! forms (`rho_ref*GRAVITY`), which is what makes
         !! `pa(nz+1) = rho_ref*g*eta_geo + p_ice_ref` cancel to bit-zero
         !! at rest.  Stored rather than recomputed so `GRAVITY`/`rho_ref`
         !! cannot drift between the two users.  CONSUMED as the static
         !! half of `multilayer_state_t%p_top = p_ice_ref + sf%p_surf` —
         !! seeded in `configure_ocean_cavity` and rebuilt each outer step
         !! in `ocean_dyn_step_split` whenever the psurf seam makes
         !! `sf%p_surf` live.  It is the load's route into the PRESSURE
         !! (the FV_MOM6 `pa(nz+1)` top BC and the in-situ EOS); its route
         !! into the BAROTROPIC mode is the datum `bt_H_ref = b - z_draft`
         !! and nothing else, which is why it never joins `sf%p_surf`.

      ! ---- Static land masks (real 0/1; derived in metrics_apply_land_mask) ----
      real(wp), allocatable :: wet_T(:, :)
         !! T-cell wet (1) / land (0) mask, `(nx,ny)` — the HALO-VALID
         !! working copy of `multilayer%wet_mask` (periodic-wrapped +
         !! north-folded, R5a) that `wet_u/wet_v/wet_q` are derived from.
         !! Kept device-resident so the continuity + tracer PPM
         !! reconstruction can mirror a land neighbour's thickness to the
         !! local cell (spec §14 C2 / MOM6's reflected-coast PPM).
         !! All-wet domain ⇒ `wet_T≡1` ⇒ mirror never triggers (no-op).
      real(wp), allocatable :: wet_u(:, :)
         !! u-face (Cu) open mask, `(nx+1,ny)`.  `wet_u(i,j) =
         !! wet_T(i-1,j)*wet_T(i,j)` — a u-face is open iff BOTH adjacent
         !! T-cells are wet.  `mass_flux_x(i,j)` is the west face of cell
         !! `(i,j)` (continuity divergence reads `flux(i+1)-flux(i)`), so
         !! the `i-1`/`i` pairing matches `dy_cu`'s stagger exactly.
         !! All-wet domain ⇒ `wet_u≡1` ⇒ masking is a literal no-op.
      real(wp), allocatable :: wet_v(:, :)
         !! v-face (Cv) open mask, `(nx,ny+1)`.  `wet_v(i,j) =
         !! wet_T(i,j-1)*wet_T(i,j)`.
      real(wp), allocatable :: wet_q(:, :)
         !! Corner (Bu) open mask, `(nx+1,ny+1)`.  Free-slip product of
         !! the 4 surrounding T-cells: `wet_q(i,j) =
         !! wet_T(i-1,j-1)*wet_T(i,j-1)*wet_T(i-1,j)*wet_T(i,j)`.  Consumed
         !! by the relative-vorticity / strain factor (CHUNK B).

      ! ---- Areas (m^2) — load-bearing; `dx*dy` is dead (D5) ----
      real(wp), allocatable :: areaT(:, :)
         !! T-cell area (m^2), `(nx,ny)`.
      real(wp), allocatable :: areaCu(:, :)
         !! Cu-cell area (m^2), `(nx+1,ny)`.
      real(wp), allocatable :: areaCv(:, :)
         !! Cv-cell area (m^2), `(nx,ny+1)`.
      real(wp), allocatable :: areaBu(:, :)
         !! Bu-cell area (m^2), `(nx+1,ny+1)`.

      ! ---- Stored inverses (Adcroft reciprocal; filled in finalize) ----
      real(wp), allocatable :: idxT(:, :), idyT(:, :)
         !! 1/dxT, 1/dyT (1/m), `(nx,ny)`.
      real(wp), allocatable :: idxCu(:, :), idyCu(:, :)
         !! 1/dxCu, 1/dyCu (1/m), `(nx+1,ny)`.
      real(wp), allocatable :: idxCv(:, :), idyCv(:, :)
         !! 1/dxCv, 1/dyCv (1/m), `(nx,ny+1)`.
      real(wp), allocatable :: iareaT(:, :)
         !! 1/areaT (1/m^2), `(nx,ny)`.
      real(wp), allocatable :: iareaBu(:, :)
         !! 1/areaBu (1/m^2), `(nx+1,ny+1)`.
      real(wp), allocatable :: iareaCu(:, :)
         !! 1/areaCu (1/m^2), `(nx+1,ny)`.
      real(wp), allocatable :: iareaCv(:, :)
         !! 1/areaCv (1/m^2), `(nx,ny+1)`.

      ! ---- Geography (degrees) ----
      real(wp), allocatable :: geolatT(:, :), geolonT(:, :)
         !! Latitude / longitude at T points (degrees), `(nx,ny)`.
      real(wp), allocatable :: geolatBu(:, :), geolonBu(:, :)
         !! Latitude / longitude at Bu corners (degrees), `(nx+1,ny+1)`.

      ! ---- hvisc ratio bundle (dimensionless / m; filled in finalize) ----
      real(wp), allocatable :: dy_dxT(:, :)
         !! dyT/dxT at T (dimensionless), `(nx,ny)`.  =1 on Cartesian.
      real(wp), allocatable :: dx_dyT(:, :)
         !! dxT/dyT at T (dimensionless), `(nx,ny)`.
      real(wp), allocatable :: dy_dxBu(:, :)
         !! dyBu/dxBu at Bu (dimensionless), `(nx+1,ny+1)`.
      real(wp), allocatable :: dx_dyBu(:, :)
         !! dxBu/dyBu at Bu (dimensionless), `(nx+1,ny+1)`.
      real(wp), allocatable :: dx2h(:, :)
         !! dxT^2 at T (m^2), `(nx,ny)`.
      real(wp), allocatable :: dy2h(:, :)
         !! dyT^2 at T (m^2), `(nx,ny)`.
      real(wp), allocatable :: dx2q(:, :)
         !! dxBu^2 at Bu (m^2), `(nx+1,ny+1)`.
      real(wp), allocatable :: dy2q(:, :)
         !! dyBu^2 at Bu (m^2), `(nx+1,ny+1)`.
   contains
      procedure, non_overridable :: init => ocean_metrics_init
      procedure, non_overridable :: destroy => ocean_metrics_destroy
      procedure, non_overridable :: enter_data => ocean_metrics_enter_data
      procedure, non_overridable :: exit_data => ocean_metrics_exit_data
      procedure, non_overridable :: bytes => ocean_metrics_bytes
   end type ocean_metrics_t

contains

   ! =================================================================
   ! Adcroft reciprocal
   ! =================================================================

   elemental pure function adcroft_recip(x) result(r)
      !! Adcroft reciprocal: `1/x`, but `0 -> 0` (zero-width faces give
      !! zero inverse, no NaN/Inf).  Single source for every metric
      !! inverse (D4).
      real(wp), intent(in) :: x
      real(wp) :: r
      if (x /= 0.0_wp) then
         r = 1.0_wp/x
      else
         r = 0.0_wp
      end if
   end function adcroft_recip

   ! =================================================================
   ! Lifecycle
   ! =================================================================

   subroutine ocean_metrics_init(this, grid)
      !! Allocate + zero every metric array.  Always allocates (configure
      !! runs after init, before `enter_data`); off-cost is ~24
      !! `(nx,ny)`-class arrays (~2 MB at Tasman size).
      class(ocean_metrics_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      integer :: nx, ny

      nx = grid%nx_total
      ny = grid%ny_total

      ! Lengths
      allocate (this%dxT(nx, ny), source=0.0_wp)
      allocate (this%dyT(nx, ny), source=0.0_wp)
      allocate (this%dxCu(nx + 1, ny), source=0.0_wp)
      allocate (this%dyCu(nx + 1, ny), source=0.0_wp)
      allocate (this%dxCv(nx, ny + 1), source=0.0_wp)
      allocate (this%dyCv(nx, ny + 1), source=0.0_wp)
      allocate (this%dxBu(nx + 1, ny + 1), source=0.0_wp)
      allocate (this%dyBu(nx + 1, ny + 1), source=0.0_wp)
      ! Topo face widths
      allocate (this%dy_cu(nx + 1, ny), source=0.0_wp)
      allocate (this%dx_cv(nx, ny + 1), source=0.0_wp)
      ! BT transport widths: ALWAYS full size (2D, cheap) so the
      ! barotropic substep can read them unconditionally — no branch, no
      ! placeholder.  `metrics_finalize` / `metrics_apply_land_mask` keep
      ! them a byte copy of `dy_cu` / `dx_cv` until porous barriers narrow
      ! them, which is what makes the knob-off path bit-identical.
      allocate (this%dy_cu_bt(nx + 1, ny), source=0.0_wp)
      allocate (this%dx_cv_bt(nx, ny + 1), source=0.0_wp)
      ! Porous barriers start at the (1,1)/(1,1,1) placeholder size:
      ! `metrics_porous_alloc` grows them at configure time (which runs
      ! after init and before enter_data) only when the knob is on, so a
      ! default run pays ~7 words instead of two full 3D fields.
      allocate (this%por_bed(1, 1), source=0.0_wp)
      allocate (this%por_dmin_u(1, 1), source=0.0_wp)
      allocate (this%por_dmax_u(1, 1), source=0.0_wp)
      allocate (this%por_davg_u(1, 1), source=0.0_wp)
      allocate (this%por_dmin_v(1, 1), source=0.0_wp)
      allocate (this%por_dmax_v(1, 1), source=0.0_wp)
      allocate (this%por_davg_v(1, 1), source=0.0_wp)
      allocate (this%por_face_area_u(1, 1, 1), source=1.0_wp)
      allocate (this%por_face_area_v(1, 1, 1), source=1.0_wp)
      ! z-level closed faces: same placeholder discipline as the porous
      ! arrays -- `metrics_closed_faces_alloc` grows them at configure
      ! (after init, before enter_data) only when the knob is on.  The
      ! placeholder is 1 (fully open) so an accidental read is inert, but
      ! it must NEVER reach an explicit-shape device dummy: every consumer
      ! names `open_u`/`open_v` only inside a branch guarded by
      ! `use_closed_faces`.
      allocate (this%open_u(1, 1, 1), source=1.0_wp)
      allocate (this%open_v(1, 1, 1), source=1.0_wp)
      ! Ice-shelf cavity statics.  Unlike the porous arrays (grown at
      ! configure), these are sized HERE off the `use_cavity` flag that
      ! `init_from_config` latches before `init` — the draft has to exist
      ! before `ocean_state_seed_from_cfg` seeds the wet mask and the
      ! layer thicknesses from `b - z_draft`, which is well before any
      ! `configure_ocean_*` runs.  Knob off ⇒ three `(1,1)` placeholders.
      if (this%use_cavity) then
         allocate (this%z_draft(nx, ny), source=0.0_wp)
         allocate (this%cover_frac(nx, ny), source=0.0_wp)
         allocate (this%p_ice_ref(nx, ny), source=0.0_wp)
      else
         allocate (this%z_draft(1, 1), source=0.0_wp)
         allocate (this%cover_frac(1, 1), source=0.0_wp)
         allocate (this%p_ice_ref(1, 1), source=0.0_wp)
      end if
      ! Land masks default ALL-WET (1.0): if metrics_apply_land_mask is
      ! never called (no land), the masks stay inert (×1) and the 6 face
      ! metrics are never altered — bit-identical to a no-mask build.
      allocate (this%wet_T(nx, ny), source=1.0_wp)
      allocate (this%wet_u(nx + 1, ny), source=1.0_wp)
      allocate (this%wet_v(nx, ny + 1), source=1.0_wp)
      allocate (this%wet_q(nx + 1, ny + 1), source=1.0_wp)
      ! Areas
      allocate (this%areaT(nx, ny), source=0.0_wp)
      allocate (this%areaCu(nx + 1, ny), source=0.0_wp)
      allocate (this%areaCv(nx, ny + 1), source=0.0_wp)
      allocate (this%areaBu(nx + 1, ny + 1), source=0.0_wp)
      ! Inverses
      allocate (this%idxT(nx, ny), source=0.0_wp)
      allocate (this%idyT(nx, ny), source=0.0_wp)
      allocate (this%idxCu(nx + 1, ny), source=0.0_wp)
      allocate (this%idyCu(nx + 1, ny), source=0.0_wp)
      allocate (this%idxCv(nx, ny + 1), source=0.0_wp)
      allocate (this%idyCv(nx, ny + 1), source=0.0_wp)
      allocate (this%iareaT(nx, ny), source=0.0_wp)
      allocate (this%iareaBu(nx + 1, ny + 1), source=0.0_wp)
      allocate (this%iareaCu(nx + 1, ny), source=0.0_wp)
      allocate (this%iareaCv(nx, ny + 1), source=0.0_wp)
      ! Geography
      allocate (this%geolatT(nx, ny), source=0.0_wp)
      allocate (this%geolonT(nx, ny), source=0.0_wp)
      allocate (this%geolatBu(nx + 1, ny + 1), source=0.0_wp)
      allocate (this%geolonBu(nx + 1, ny + 1), source=0.0_wp)
      ! hvisc ratio bundle
      allocate (this%dy_dxT(nx, ny), source=0.0_wp)
      allocate (this%dx_dyT(nx, ny), source=0.0_wp)
      allocate (this%dy_dxBu(nx + 1, ny + 1), source=0.0_wp)
      allocate (this%dx_dyBu(nx + 1, ny + 1), source=0.0_wp)
      allocate (this%dx2h(nx, ny), source=0.0_wp)
      allocate (this%dy2h(nx, ny), source=0.0_wp)
      allocate (this%dx2q(nx + 1, ny + 1), source=0.0_wp)
      allocate (this%dy2q(nx + 1, ny + 1), source=0.0_wp)

      this%is_init = .true.
   end subroutine ocean_metrics_init

   subroutine ocean_metrics_destroy(this)
      class(ocean_metrics_t), intent(inout) :: this
      this%is_init = .false.
      if (allocated(this%dxT)) deallocate (this%dxT)
      if (allocated(this%dyT)) deallocate (this%dyT)
      if (allocated(this%dxCu)) deallocate (this%dxCu)
      if (allocated(this%dyCu)) deallocate (this%dyCu)
      if (allocated(this%dxCv)) deallocate (this%dxCv)
      if (allocated(this%dyCv)) deallocate (this%dyCv)
      if (allocated(this%dxBu)) deallocate (this%dxBu)
      if (allocated(this%dyBu)) deallocate (this%dyBu)
      if (allocated(this%dy_cu)) deallocate (this%dy_cu)
      if (allocated(this%dx_cv)) deallocate (this%dx_cv)
      if (allocated(this%dy_cu_bt)) deallocate (this%dy_cu_bt)
      if (allocated(this%dx_cv_bt)) deallocate (this%dx_cv_bt)
      if (allocated(this%por_bed)) deallocate (this%por_bed)
      if (allocated(this%por_dmin_u)) deallocate (this%por_dmin_u)
      if (allocated(this%por_dmax_u)) deallocate (this%por_dmax_u)
      if (allocated(this%por_davg_u)) deallocate (this%por_davg_u)
      if (allocated(this%por_dmin_v)) deallocate (this%por_dmin_v)
      if (allocated(this%por_dmax_v)) deallocate (this%por_dmax_v)
      if (allocated(this%por_davg_v)) deallocate (this%por_davg_v)
      if (allocated(this%por_face_area_u)) deallocate (this%por_face_area_u)
      if (allocated(this%por_face_area_v)) deallocate (this%por_face_area_v)
      if (allocated(this%open_u)) deallocate (this%open_u)
      if (allocated(this%open_v)) deallocate (this%open_v)
      if (allocated(this%z_draft)) deallocate (this%z_draft)
      if (allocated(this%cover_frac)) deallocate (this%cover_frac)
      if (allocated(this%p_ice_ref)) deallocate (this%p_ice_ref)
      if (allocated(this%wet_T)) deallocate (this%wet_T)
      if (allocated(this%wet_u)) deallocate (this%wet_u)
      if (allocated(this%wet_v)) deallocate (this%wet_v)
      if (allocated(this%wet_q)) deallocate (this%wet_q)
      if (allocated(this%areaT)) deallocate (this%areaT)
      if (allocated(this%areaCu)) deallocate (this%areaCu)
      if (allocated(this%areaCv)) deallocate (this%areaCv)
      if (allocated(this%areaBu)) deallocate (this%areaBu)
      if (allocated(this%idxT)) deallocate (this%idxT)
      if (allocated(this%idyT)) deallocate (this%idyT)
      if (allocated(this%idxCu)) deallocate (this%idxCu)
      if (allocated(this%idyCu)) deallocate (this%idyCu)
      if (allocated(this%idxCv)) deallocate (this%idxCv)
      if (allocated(this%idyCv)) deallocate (this%idyCv)
      if (allocated(this%iareaT)) deallocate (this%iareaT)
      if (allocated(this%iareaBu)) deallocate (this%iareaBu)
      if (allocated(this%iareaCu)) deallocate (this%iareaCu)
      if (allocated(this%iareaCv)) deallocate (this%iareaCv)
      if (allocated(this%geolatT)) deallocate (this%geolatT)
      if (allocated(this%geolonT)) deallocate (this%geolonT)
      if (allocated(this%geolatBu)) deallocate (this%geolatBu)
      if (allocated(this%geolonBu)) deallocate (this%geolonBu)
      if (allocated(this%dy_dxT)) deallocate (this%dy_dxT)
      if (allocated(this%dx_dyT)) deallocate (this%dx_dyT)
      if (allocated(this%dy_dxBu)) deallocate (this%dy_dxBu)
      if (allocated(this%dx_dyBu)) deallocate (this%dx_dyBu)
      if (allocated(this%dx2h)) deallocate (this%dx2h)
      if (allocated(this%dy2h)) deallocate (this%dy2h)
      if (allocated(this%dx2q)) deallocate (this%dx2q)
      if (allocated(this%dy2q)) deallocate (this%dy2q)
   end subroutine ocean_metrics_destroy

   subroutine metrics_porous_alloc(this, grid, nz)
      !! Grow the porous-barrier arrays from their `(1,1)`/`(1,1,1)`
      !! placeholder size to full face size.  Call ONLY when
      !! `&ocean_porous_nml enable` is on, at configure time — i.e. after
      !! `init` and BEFORE `ocean_state_enter_data`, so the device map
      !! captures the final shapes (a realloc after `enter_data` would
      !! leave the device pointing at freed host memory).
      !!
      !! The open fractions start at 1 (fully open) so that a stage which
      !! reads them before the first recompute sees an inert scheme.
      type(ocean_metrics_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
         !! Number of layers (`multilayer%nz_ml`).

      integer :: nx, ny

      nx = grid%nx_total
      ny = grid%ny_total

      if (allocated(this%por_bed)) deallocate (this%por_bed)
      if (allocated(this%por_dmin_u)) deallocate (this%por_dmin_u)
      if (allocated(this%por_dmax_u)) deallocate (this%por_dmax_u)
      if (allocated(this%por_davg_u)) deallocate (this%por_davg_u)
      if (allocated(this%por_dmin_v)) deallocate (this%por_dmin_v)
      if (allocated(this%por_dmax_v)) deallocate (this%por_dmax_v)
      if (allocated(this%por_davg_v)) deallocate (this%por_davg_v)
      if (allocated(this%por_face_area_u)) deallocate (this%por_face_area_u)
      if (allocated(this%por_face_area_v)) deallocate (this%por_face_area_v)

      allocate (this%por_bed(nx, ny), source=0.0_wp)
      allocate (this%por_dmin_u(nx + 1, ny), source=0.0_wp)
      allocate (this%por_dmax_u(nx + 1, ny), source=0.0_wp)
      allocate (this%por_davg_u(nx + 1, ny), source=0.0_wp)
      allocate (this%por_dmin_v(nx, ny + 1), source=0.0_wp)
      allocate (this%por_dmax_v(nx, ny + 1), source=0.0_wp)
      allocate (this%por_davg_v(nx, ny + 1), source=0.0_wp)
      allocate (this%por_face_area_u(nx + 1, ny, nz), source=1.0_wp)
      allocate (this%por_face_area_v(nx, ny + 1, nz), source=1.0_wp)
   end subroutine metrics_porous_alloc

   subroutine metrics_closed_faces_alloc(this, grid, nz)
      !! Grow the z-level closed-face masks from their `(1,1,1)`
      !! placeholder to full face size.  Call ONLY when
      !! `&vcoord_nml zfixed_closed_faces` is on, at configure time —
      !! after `init` and BEFORE `ocean_state_enter_data`, so the device
      !! map captures the final shapes (a realloc after `enter_data`
      !! would leave the device pointing at freed host memory).
      !!
      !! Seeded fully OPEN (1) so a stage that somehow reads them before
      !! the builder runs sees an inert mask.
      type(ocean_metrics_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
         !! Number of layers (`multilayer%nz_ml`).

      integer :: nx, ny

      nx = grid%nx_total
      ny = grid%ny_total

      if (allocated(this%open_u)) deallocate (this%open_u)
      if (allocated(this%open_v)) deallocate (this%open_v)

      allocate (this%open_u(nx + 1, ny, nz), source=1.0_wp)
      allocate (this%open_v(nx, ny + 1, nz), source=1.0_wp)
   end subroutine metrics_closed_faces_alloc

   subroutine ocean_metrics_enter_data(this)
      class(ocean_metrics_t), intent(inout) :: this
      select type (this)
      type is (ocean_metrics_t)
         call ocean_metrics_enter_data_impl(this)
      end select
   end subroutine ocean_metrics_enter_data

   subroutine ocean_metrics_enter_data_impl(this)
      !! Arrays-only attach.  The parent `ocean_state_t` is mapped by the
      !! orchestrator BEFORE this runs.
      type(ocean_metrics_t), intent(inout) :: this
      !$acc enter data copyin(this%dxT, this%dyT, this%dxCu, this%dyCu)
      !$acc enter data copyin(this%dxCv, this%dyCv, this%dxBu, this%dyBu)
      !$acc enter data copyin(this%dy_cu, this%dx_cv)
      !$acc enter data copyin(this%dy_cu_bt, this%dx_cv_bt)
      !$acc enter data copyin(this%por_bed)
      !$acc enter data copyin(this%por_dmin_u, this%por_dmax_u, this%por_davg_u)
      !$acc enter data copyin(this%por_dmin_v, this%por_dmax_v, this%por_davg_v)
      !$acc enter data copyin(this%por_face_area_u, this%por_face_area_v)
      !$acc enter data copyin(this%open_u, this%open_v)
      !$acc enter data copyin(this%z_draft, this%cover_frac, this%p_ice_ref)
      !$acc enter data copyin(this%wet_T, this%wet_u, this%wet_v, this%wet_q)
      !$acc enter data copyin(this%areaT, this%areaCu, this%areaCv, this%areaBu)
      !$acc enter data copyin(this%idxT, this%idyT, this%idxCu, this%idyCu)
      !$acc enter data copyin(this%idxCv, this%idyCv)
      !$acc enter data copyin(this%iareaT, this%iareaBu, this%iareaCu, this%iareaCv)
      !$acc enter data copyin(this%geolatT, this%geolonT, this%geolatBu, this%geolonBu)
      !$acc enter data copyin(this%dy_dxT, this%dx_dyT, this%dy_dxBu, this%dx_dyBu)
      !$acc enter data copyin(this%dx2h, this%dy2h, this%dx2q, this%dy2q)
   end subroutine ocean_metrics_enter_data_impl

   subroutine ocean_metrics_exit_data(this)
      class(ocean_metrics_t), intent(inout) :: this
      select type (this)
      type is (ocean_metrics_t)
         call ocean_metrics_exit_data_impl(this)
      end select
   end subroutine ocean_metrics_exit_data

   subroutine ocean_metrics_exit_data_impl(this)
      type(ocean_metrics_t), intent(inout) :: this
      !$acc exit data delete(this%dx2h, this%dy2h, this%dx2q, this%dy2q)
      !$acc exit data delete(this%dy_dxT, this%dx_dyT, this%dy_dxBu, this%dx_dyBu)
      !$acc exit data delete(this%geolatT, this%geolonT, this%geolatBu, this%geolonBu)
      !$acc exit data delete(this%iareaT, this%iareaBu, this%iareaCu, this%iareaCv)
      !$acc exit data delete(this%idxCv, this%idyCv)
      !$acc exit data delete(this%idxT, this%idyT, this%idxCu, this%idyCu)
      !$acc exit data delete(this%areaT, this%areaCu, this%areaCv, this%areaBu)
      !$acc exit data delete(this%wet_T, this%wet_u, this%wet_v, this%wet_q)
      !$acc exit data delete(this%z_draft, this%cover_frac, this%p_ice_ref)
      !$acc exit data delete(this%open_u, this%open_v)
      !$acc exit data delete(this%por_face_area_u, this%por_face_area_v)
      !$acc exit data delete(this%por_dmin_v, this%por_dmax_v, this%por_davg_v)
      !$acc exit data delete(this%por_dmin_u, this%por_dmax_u, this%por_davg_u)
      !$acc exit data delete(this%por_bed)
      !$acc exit data delete(this%dy_cu_bt, this%dx_cv_bt)
      !$acc exit data delete(this%dy_cu, this%dx_cv)
      !$acc exit data delete(this%dxCv, this%dyCv, this%dxBu, this%dyBu)
      !$acc exit data delete(this%dxT, this%dyT, this%dxCu, this%dyCu)
   end subroutine ocean_metrics_exit_data_impl

   ! =================================================================
   ! Finalize — single-source the inverses + the hvisc ratio bundle
   ! =================================================================

   subroutine metrics_finalize(this)
      !! Compute every stored inverse + the hvisc ratio bundle ONCE from
      !! the length/area arrays a generator already wrote, via the
      !! Adcroft reciprocal (D4).  No kernel ever recomputes these.
      type(ocean_metrics_t), intent(inout) :: this

      this%idxT = adcroft_recip(this%dxT)
      this%idyT = adcroft_recip(this%dyT)
      this%idxCu = adcroft_recip(this%dxCu)
      this%idyCu = adcroft_recip(this%dyCu)
      this%idxCv = adcroft_recip(this%dxCv)
      this%idyCv = adcroft_recip(this%dyCv)
      this%iareaT = adcroft_recip(this%areaT)
      this%iareaBu = adcroft_recip(this%areaBu)
      this%iareaCu = adcroft_recip(this%areaCu)
      this%iareaCv = adcroft_recip(this%areaCv)

      ! hvisc ratio bundle (all == 1 on uniform Cartesian).
      this%dy_dxT = this%dyT*this%idxT
      this%dx_dyT = this%dxT*this%idyT
      this%dy_dxBu = this%dyBu*adcroft_recip(this%dxBu)
      this%dx_dyBu = this%dxBu*adcroft_recip(this%dyBu)
      this%dx2h = this%dxT*this%dxT
      this%dy2h = this%dyT*this%dyT
      this%dx2q = this%dxBu*this%dxBu
      this%dy2q = this%dyBu*this%dyBu

      ! BT transport widths start as an exact copy of the slow-path
      ! widths.  Re-synced after land masking (which zeroes `dy_cu` /
      ! `dx_cv`) and overwritten per step only when porous barriers are on.
      this%dy_cu_bt = this%dy_cu
      this%dx_cv_bt = this%dx_cv
   end subroutine metrics_finalize

   ! =================================================================
   ! Static land masking (CHUNK A foundation)
   ! =================================================================

   subroutine metrics_apply_land_mask(this, wet_mask, grid, &
                                      periodic_x, periodic_y, north_fold, &
                                      mask_wall_velocity, &
                                      wall_west, wall_east, wall_south, wall_north)
      !! Derive the static C-grid face / corner masks from the T-cell
      !! `wet_mask` and zero the face metrics at land faces, so every
      !! transport / gradient / circulation operator that rides those
      !! metrics couples across NO land face (MOM6 pre-masks the face
      !! LENGTHS; Adcroft & Hallberg 2006).
      !!
      !! Must run at SETUP, AFTER `metrics_finalize` (the inverses
      !! `idxCu`/`idyCv` are masked here, so they must already exist) and
      !! AFTER `wet_mask` is seeded, but BEFORE `ocean_state_enter_data`
      !! (the host edit is what the GPU copyin captures).  Plain host
      !! loops — `do concurrent` before `enter_data` would round-trip the
      !! unmapped arrays through the device per loop.
      !!
      !! Halo-aware (R5a): a working copy of `wet_mask` is first filled in
      !! the ghost columns/rows by the SAME periodic wrap / north fold the
      !! metric ghosts use, so the seam u-faces (e.g. a continent that
      !! straddles `x=0≡x=1`) mask correctly.  Wall ghosts already carry
      !! the constant-extrapolated `wet_mask` from the bathymetry fill.
      !!
      !! Masks: `wet_u(i,j) = wet_T(i-1,j)*wet_T(i,j)` (Cu),
      !! `wet_v(i,j) = wet_T(i,j-1)*wet_T(i,j)` (Cv),
      !! `wet_q(i,j) = product of the 4 T-cells around corner (i,j)` (Bu,
      !! free-slip).  Zeroed metrics (the EXACT 6 — spec §14 C3):
      !!   `dy_cu, idxCu, dxCu` at `wet_u==0`;
      !!   `dx_cv, idyCv, dyCv` at `wet_v==0`.
      !! NOT touched: `iareaT, areaT, areaCu, areaCv, iareaBu` (zeroing
      !! them would break wet-cell divergence / KE / Coriolis
      !! corner-area normalization).
      !!
      !! Bit-identity: all-wet ⇒ every `wet_*≡1` ⇒ the 6 metrics are
      !! multiplied by 1 (byte-unchanged) and the masks stay inert.
      !!
      !! Solid-wall velocity masking (`mask_wall_velocity`, opt-in): a flat
      !! all-wet channel has `wet_mask≡1` in the WALL-edge ghosts too, so
      !! `wet_v`/`wet_u` at the wall face = 1·1 = 1 (unmasked) — the wall
      !! flux is masked but the raw wall-normal velocity drifts to garbage
      !! (spurious vorticity band).  When enabled, the ghost `wm` beyond
      !! each SOLID WALL edge is zeroed (MOM6: the halo beyond a wall is
      !! land), so the derived face masks are 0 at the wall and the existing
      !! per-stage `mask_layer_velocities` clears the velocity — no bespoke
      !! velocity BC.  Only WALL edges are touched: periodic edges keep
      !! their wrapped (wet) ghosts, open/OBC edges keep the interior value.
      type(ocean_metrics_t), intent(inout) :: this
      real(wp), intent(in) :: wet_mask(:, :)
         !! T-cell wet (1) / land (0) mask, `(nx_total, ny_total)`.
      type(hgrid_t), intent(in) :: grid
      logical, intent(in) :: periodic_x, periodic_y, north_fold
         !! Boundary topology of `wet_mask`'s ghost halo (from the bc
         !! state) — selects the ghost wrap before deriving the masks.
      logical, intent(in), optional :: mask_wall_velocity
         !! Opt-in solid-wall velocity masking (default absent ⇒ .false. ⇒
         !! wall ghosts untouched ⇒ bit-identical to the legacy path).
      logical, intent(in), optional :: wall_west, wall_east, wall_south, wall_north
         !! Per-edge solid-WALL flags (an edge that is NOT periodic, NOT
         !! north-fold, NOT open/OBC).  Only consulted when
         !! `mask_wall_velocity` is .true.; each defaults to "wall" on any
         !! non-periodic / non-fold edge (the closed-default assumption).

      integer :: nx, ny, ni, nj, ng, i, j
      real(wp), allocatable :: wm(:, :)
      logical :: do_wall_mask, w_wall, e_wall, s_wall, n_wall

      nx = grid%nx_total
      ny = grid%ny_total
      ni = grid%nx_phys
      nj = grid%ny_phys
      ng = grid%nghost

      ! ---- Halo-valid working copy of wet_mask (R5a) ----
      ! The incoming `wet_mask` already carries multi-rank seam ghosts: the
      ! caller (configure_ocean_land_mask) exchanges it via ocean_halo_centre
      ! BEFORE this routine so the seam ghost columns hold the neighbour rank's
      ! real wet_T (O3 land x decomp).  Here we only add the PHYSICAL periodic /
      ! fold ghost wraps (disjoint from MPI seams).
      allocate (wm(nx, ny))
      wm = wet_mask
      if (periodic_x) call metrics_periodic_x_2d(wm, grid)
      if (periodic_y) call metrics_periodic_y_2d(wm, grid)
      if (north_fold) call fold_north_centre(wm, nx, ny, ni, nj, ng)

      ! ---- Solid-wall land fill (opt-in; MOM6 mask-in-the-update) ----
      ! Zero the ghost rows/columns beyond each SOLID WALL edge AFTER the
      ! periodic/fold wraps and BEFORE deriving the face masks, so the wall
      ! face products (wet_v/wet_u = wm(interior)*wm(ghost) = *0) vanish and
      ! `mask_layer_velocities` clears the wall-normal velocity each stage.
      ! Untouched when disabled ⇒ bit-identical.  Only WALL edges: periodic
      ! edges keep their wrapped ghosts, open/OBC edges the interior value.
      do_wall_mask = .false.
      if (present(mask_wall_velocity)) do_wall_mask = mask_wall_velocity
      if (do_wall_mask) then
         ! Default: any non-periodic, non-fold edge is a wall (closed default);
         ! callers thread the true per-edge WALL/OPEN flags to override.
         w_wall = .not. periodic_x
         e_wall = .not. periodic_x
         s_wall = .not. periodic_y
         n_wall = (.not. periodic_y) .and. (.not. north_fold)
         if (present(wall_west)) w_wall = wall_west
         if (present(wall_east)) e_wall = wall_east
         if (present(wall_south)) s_wall = wall_south
         if (present(wall_north)) n_wall = wall_north
         ! Ghost cells: west i=1..ng, east i=ng+ni+1..nx,
         !              south j=1..ng, north j=ng+nj+1..ny.
         if (w_wall) then
            do j = 1, ny
               do i = 1, ng
                  wm(i, j) = 0.0_wp
               end do
            end do
         end if
         if (e_wall) then
            do j = 1, ny
               do i = ng + ni + 1, nx
                  wm(i, j) = 0.0_wp
               end do
            end do
         end if
         if (s_wall) then
            do j = 1, ng
               do i = 1, nx
                  wm(i, j) = 0.0_wp
               end do
            end do
         end if
         if (n_wall) then
            do j = ng + nj + 1, ny
               do i = 1, nx
                  wm(i, j) = 0.0_wp
               end do
            end do
         end if
      end if

      ! ---- Store the halo-valid T-cell mask (consumed by PPM mirror-h) ----
      this%wet_T = wm

      ! ---- Derive face / corner masks (plain host loops) ----
      ! wet_u(i,j): u-face i = west face of T-cell (i,j); pairs (i-1,i).
      ! Array outer ring (i=1, nx+1) is outside the ghost band => always
      ! land, decomposition-invariant; physical/seam interface faces at
      ! nghost+1 are bathymetry-masked (wet_T product below), not edge-
      ! position-masked, so no has_west/has_east gate is needed here
      ! (O0 verified: seam face wet_T product = 1*1 = 1, mask stays open).
      do j = 1, ny
         this%wet_u(1, j) = 0.0_wp      ! west outer wall (no T-cell i=0)
         do i = 2, nx
            this%wet_u(i, j) = wm(i - 1, j)*wm(i, j)
         end do
         this%wet_u(nx + 1, j) = 0.0_wp  ! east outer wall (no T-cell nx+1)
      end do
      ! wet_v(i,j): v-face j = south face of T-cell (i,j); pairs (j-1,j).
      do i = 1, nx
         this%wet_v(i, 1) = 0.0_wp
         do j = 2, ny
            this%wet_v(i, j) = wm(i, j - 1)*wm(i, j)
         end do
         this%wet_v(i, ny + 1) = 0.0_wp
      end do
      ! wet_q(i,j): SW corner of T-cell (i,j); product of the 4 T-cells
      ! (i-1,j-1),(i,j-1),(i-1,j),(i,j).  Outer ring (i=1/nx+1, j=1/ny+1)
      ! has a missing T-neighbour ⇒ land (matches the domain-wall corner).
      do j = 1, ny + 1
         do i = 1, nx + 1
            if (i >= 2 .and. i <= nx .and. j >= 2 .and. j <= ny) then
               this%wet_q(i, j) = wm(i - 1, j - 1)*wm(i, j - 1)* &
                                  wm(i - 1, j)*wm(i, j)
            else
               this%wet_q(i, j) = 0.0_wp
            end if
         end do
      end do

      ! ---- Zero the 6 face metrics at land faces (spec §14 C3) ----
      ! u-faces (Cu): dy_cu, idxCu, dxCu.
      do j = 1, ny
         do i = 1, nx + 1
            this%dy_cu(i, j) = this%dy_cu(i, j)*this%wet_u(i, j)
            this%idxCu(i, j) = this%idxCu(i, j)*this%wet_u(i, j)
            this%dxCu(i, j) = this%dxCu(i, j)*this%wet_u(i, j)
         end do
      end do
      ! v-faces (Cv): dx_cv, idyCv, dyCv.
      do j = 1, ny + 1
         do i = 1, nx
            this%dx_cv(i, j) = this%dx_cv(i, j)*this%wet_v(i, j)
            this%idyCv(i, j) = this%idyCv(i, j)*this%wet_v(i, j)
            this%dyCv(i, j) = this%dyCv(i, j)*this%wet_v(i, j)
         end do
      end do

      ! Re-sync the BT transport widths with the freshly masked slow-path
      ! widths (this runs AFTER metrics_finalize, which set them equal).
      this%dy_cu_bt = this%dy_cu
      this%dx_cv_bt = this%dx_cv

      deallocate (wm)
   end subroutine metrics_apply_land_mask

   subroutine metrics_periodic_y_2d(arr, grid)
      !! South/north ghost rows of a T-array by periodic wrap (the y
      !! analogue of `metrics_periodic_x_2d`).  Only used by the land-mask
      !! ghost fill; the metric tripolar path wraps x only.
      real(wp), intent(inout) :: arr(:, :)
      type(hgrid_t), intent(in) :: grid
      integer :: ng, nj, i, j
      ng = grid%nghost
      nj = grid%ny_phys
      do j = 1, ng
         do i = 1, size(arr, 1)
            arr(i, j) = arr(i, j + nj)
            arr(i, ng + nj + j) = arr(i, ng + j)
         end do
      end do
   end subroutine metrics_periodic_y_2d

   ! =================================================================
   ! Generator: uniform Cartesian (bit-identity reference)
   ! =================================================================

   subroutine metrics_fill_cartesian(this, grid, dx, dy)
      !! Uniform Cartesian: every length is constant, `areaX = dx*dy`.
      !! Geography is left at zero (a Cartesian beta-plane has no lat/lon
      !! — the Coriolis fill uses the Cartesian y coordinate, D7).  Fills
      !! all ghost rows/columns (constants, trivially).  Call
      !! `metrics_finalize` afterwards.
      type(ocean_metrics_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: dx, dy

      this%dxT = dx
      this%dyT = dy
      this%dxCu = dx
      this%dyCu = dy
      this%dxCv = dx
      this%dyCv = dy
      this%dxBu = dx
      this%dyBu = dy
      this%dy_cu = dy
      this%dx_cv = dx
      this%areaT = dx*dy
      this%areaCu = dx*dy
      this%areaCv = dx*dy
      this%areaBu = dx*dy
      ! `grid` is part of the generator contract (shapes already come
      ! from `init(grid)`); cartesian metrics are pure constants so the
      ! geometry is not re-read here.
      associate (unused => grid%nx_total)
      end associate
   end subroutine metrics_fill_cartesian

   ! =================================================================
   ! Generator: spherical lon-lat sector (analytic-derivative form, D6)
   ! =================================================================

   subroutine metrics_fill_spherical(this, grid, lon_west, lat_south, &
                                     dlon_deg, dlat_deg, rad_earth)
      !! Spherical lon-lat sector.  For each stagger, geolat/geolon are
      !! evaluated at THAT point's own location; the metric lengths use
      !! the cos of that stagger's own latitude (the consistency trick
      !! that keeps the C-grid metrics compatible, D6):
      !!   dx = rad_earth * cos(lat) * dlon_rad
      !!   dy = rad_earth * dlat_rad
      !!   area = dx * dy   (analytic-derivative form, NOT great-circle).
      !!
      !! Indexing: the first INTERIOR T cell is `(1+nghost, 1+nghost)`,
      !! centred at `(lon_west + (i+i_offset_global-nghost-0.5)*dlon,
      !! lat_south + (j+j_offset_global-nghost-0.5)*dlat)`.  On an
      !! undecomposed grid the offsets are 0 and the formula reduces to
      !! the original single-rank form.  Under MPI decomposition the
      !! offsets shift the local (i,j) to the correct GLOBAL coordinate
      !! so every rank computes the right geolat/geolon.  Corners (Bu)
      !! sit half a cell up/right of their cell centre.  u-faces share
      !! the T latitude, v-faces / corners use the corner latitude.
      !! ALL ghost rows/columns are filled (the formula extends naturally
      !! past the physical sector).
      type(ocean_metrics_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: lon_west, lat_south, dlon_deg, dlat_deg, rad_earth

      integer :: i, j, nx, ny, ng
      real(wp) :: dlon_rad, dlat_rad, dy_len
      real(wp) :: lat_t, lon_t, lat_b, lon_b

      nx = grid%nx_total
      ny = grid%ny_total
      ng = grid%nghost
      dlon_rad = dlon_deg*DEG2RAD
      dlat_rad = dlat_deg*DEG2RAD
      dy_len = rad_earth*dlat_rad   ! meridional length is lat-independent

      ! ---- T points + u-faces (share the T-row latitude) ----
      do j = 1, ny
         lat_t = lat_south + (real(j + grid%j_offset_global - ng, wp) - 0.5_wp)*dlat_deg
         do i = 1, nx
            lon_t = lon_west + (real(i + grid%i_offset_global - ng, wp) - 0.5_wp)*dlon_deg
            this%geolatT(i, j) = lat_t
            this%geolonT(i, j) = lon_t
            this%dxT(i, j) = rad_earth*cos(lat_t*DEG2RAD)*dlon_rad
            this%dyT(i, j) = dy_len
            this%areaT(i, j) = this%dxT(i, j)*dy_len
         end do
      end do

      ! Cu (u-face): x-face on the T-row latitude.  Cu(i,j) sits on the
      ! west edge of T(i,j); its length uses lat_t (same row).
      !
      ! Uniform-dlon coincidence note: the correct T-to-T dxCu definition is
      ! the sum of the two half-segments straddling the face node (one from
      ! the cell to the west, one from the cell to the east).  On this analytic
      ! generator DLON is constant, so both half-segments are equal and the
      ! two-half-segment sum = R·cos(lat_face)·dlon_rad — which is exactly
      ! what this formula computes (lat_face = lat_t for Cu on the T row).
      ! For a variable-resolution supergrid the two definitions diverge; that
      ! path is corrected in `metrics_fill_from_supergrid` (see comment there).
      do j = 1, ny
         lat_t = lat_south + (real(j + grid%j_offset_global - ng, wp) - 0.5_wp)*dlat_deg
         do i = 1, nx + 1
            this%dxCu(i, j) = rad_earth*cos(lat_t*DEG2RAD)*dlon_rad
            this%dyCu(i, j) = dy_len
            this%dy_cu(i, j) = dy_len
            this%areaCu(i, j) = this%dxCu(i, j)*dy_len
         end do
      end do

      ! Cv (v-face) + Bu (corner): on the corner latitude row.
      ! dyCv uniform-dlat coincidence note: the T-to-T dyCv definition is the
      ! sum of two half-segments straddling the face node row.  On this analytic
      ! generator DLAT is constant so the sum = R·dlat_rad = dy_len — identical
      ! to what is coded here.  Variable-resolution corrected in supergrid reader.
      do j = 1, ny + 1
         lat_b = lat_south + real(j + grid%j_offset_global - ng - 1, wp)*dlat_deg
         do i = 1, nx
            this%dxCv(i, j) = rad_earth*cos(lat_b*DEG2RAD)*dlon_rad
            this%dyCv(i, j) = dy_len
            this%dx_cv(i, j) = this%dxCv(i, j)
            this%areaCv(i, j) = this%dxCv(i, j)*dy_len
         end do
      end do

      do j = 1, ny + 1
         lat_b = lat_south + real(j + grid%j_offset_global - ng - 1, wp)*dlat_deg
         do i = 1, nx + 1
            lon_b = lon_west + real(i + grid%i_offset_global - ng - 1, wp)*dlon_deg
            this%geolatBu(i, j) = lat_b
            this%geolonBu(i, j) = lon_b
            this%dxBu(i, j) = rad_earth*cos(lat_b*DEG2RAD)*dlon_rad
            this%dyBu(i, j) = dy_len
            this%areaBu(i, j) = this%dxBu(i, j)*dy_len
         end do
      end do
   end subroutine metrics_fill_spherical

   ! =================================================================
   ! Generator: MOM6 supergrid (mosaic) NetCDF reader
   ! =================================================================
   !
   ! Supergrid index convention (1-based, physical domain i∈[1,ni],
   ! j∈[1,nj]; supergrid node (1,1) = SW corner of the physical domain):
   !
   !   node (2i-1, 2j-1): SW corner of T(i,j) = Bu(i,j) stagger  (ODD/ODD)
   !   node (2i,   2j  ): T-cell centre                            (EVEN/EVEN)
   !   node (2i-1, 2j  ): u-face (Cu) midpoint of west face T(i,j) (ODD/EVEN)
   !   node (2i,   2j-1): v-face (Cv) midpoint of south face T(i,j) (EVEN/ODD)
   !
   !   Stagger       sg (s_i, s_j)        model range
   !   T(i,j)        (2i,   2j  )         i∈[1,ni], j∈[1,nj]
   !   Bu(i,j)       (2i-1, 2j-1)         i∈[1,ni+1], j∈[1,nj+1]
   !   Cu(i,j)       (2i-1, 2j  )         i∈[1,ni+1], j∈[1,nj]
   !   Cv(i,j)       (2i,   2j-1)         i∈[1,ni], j∈[1,nj+1]
   !
   ! Segment arrays:
   !   dx(m,n)  : along-i segment from node (m,n) to (m+1,n).
   !              Shape (2ni, 2nj+1), m∈[1,2ni], n∈[1,2nj+1].
   !   dy(m,n)  : along-j segment from node (m,n) to (m,n+1).
   !              Shape (2ni+1, 2nj), m∈[1,2ni+1], n∈[1,2nj].
   !   area(m,n): sub-cell area, SW corner at node (m,n).
   !              Shape (2ni, 2nj), m∈[1,2ni], n∈[1,2nj].
   !
   ! Metric sums (1-based physical i,j; ng offsets applied in code):
   !   dxT(i,j)   = dx(2i-1,2j) + dx(2i,  2j)
   !   dyT(i,j)   = dy(2i,2j-1) + dy(2i,  2j)
   !   areaT(i,j) = area(2i-1,2j-1)+area(2i,2j-1)+area(2i-1,2j)+area(2i,2j)
   !   dxBu(i,j)  = dx(2i-2,2j-1) + dx(2i-1,2j-1)  [i≥2; i=1 extrapolated]
   !   dyBu(i,j)  = dy(2i-1,2j-2) + dy(2i-1,2j-1)  [j≥2; j=1 extrapolated]
   !   dxCu(i,j)  = dx(2i-1,2j) + dx(2i,2j) = dxT(i,j)
   !   dyCu(i,j)  = dy(2i-1,2j-1) + dy(2i-1,2j)  [col 2i-1]
   !   dxCv(i,j)  = dx(2i-1,2j-1) + dx(2i,2j-1)  [row 2j-1]
   !   dyCv(i,j)  = dy(2i,2j-1) + dy(2i,2j) = dyT(i,j)
   !
   ! Ghost-row fill: the file covers the PHYSICAL domain only.
   ! Ghost rows/columns are filled by constant extrapolation of the nearest
   ! physical value.  Periodic/fold ghost metric fill is the M4 exchange job.

   subroutine metrics_fill_from_supergrid(this, grid, supergrid_file, ierr)
      !! Load an MOM6 supergrid (mosaic) NetCDF file and fill all metric
      !! arrays.  After this call the caller must invoke `metrics_finalize`
      !! to compute the inverses + hvisc ratio bundle.
      !!
      !! The file must contain variables `x`, `y` (degrees, shape
      !! `(2*ni+1, 2*nj+1)`), `dx` (m, shape `(2*ni, 2*nj+1)`), `dy` (m,
      !! shape `(2*ni+1, 2*nj)`), and `area` (m^2, shape `(2*ni, 2*nj)`),
      !! where `ni = grid%nx_phys`, `nj = grid%ny_phys`.  Dimensions must
      !! be named `nxp`/`nyp` (size 2ni+1 / 2nj+1) and `nx`/`ny` (size
      !! 2ni / 2nj).  Ghost rows are filled by constant extrapolation.
      type(ocean_metrics_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
         !! Grid metadata — supplies `nx_phys`, `ny_phys`, `nghost`.
      character(len=*), intent(in) :: supergrid_file
         !! Path to the MOM6 mosaic supergrid NetCDF file.
      integer, intent(out), optional :: ierr
         !! Non-zero on a dimension mismatch, an unreadable/missing file,
         !! or a missing NetCDF build when present; absent behaves as
         !! today (`error stop`).
#ifndef RDB_NO_NETCDF

      integer :: ncid
      integer :: ni, nj, ng
      integer :: sg_nxp, sg_nyp, sg_nx, sg_ny
      integer :: varid_x, varid_y, varid_dx, varid_dy, varid_area
      integer :: local_ierr
      real(wp), allocatable :: sg_x(:, :), sg_y(:, :)
      real(wp), allocatable :: sg_dx(:, :), sg_dy(:, :), sg_area(:, :)

      ni = grid%nx_phys
      nj = grid%ny_phys
      ng = grid%nghost

      call logger%info("Loading supergrid metrics from: "//trim(supergrid_file))
      call nc_open_read(supergrid_file, ncid, ierr=local_ierr)
      if (.not. supergrid_io_ok(local_ierr, ierr)) return

      ! ---- Validate supergrid dimensions against the model grid ----
      call nc_get_dim_len(ncid, "nxp", sg_nxp, ierr=local_ierr)
      if (.not. supergrid_io_ok(local_ierr, ierr, ncid)) return
      call nc_get_dim_len(ncid, "nyp", sg_nyp, ierr=local_ierr)
      if (.not. supergrid_io_ok(local_ierr, ierr, ncid)) return
      call nc_get_dim_len(ncid, "nx", sg_nx, ierr=local_ierr)
      if (.not. supergrid_io_ok(local_ierr, ierr, ncid)) return
      call nc_get_dim_len(ncid, "ny", sg_ny, ierr=local_ierr)
      if (.not. supergrid_io_ok(local_ierr, ierr, ncid)) return

      if (sg_nxp /= 2*ni + 1) then
         call nc_close(ncid)
         call fail("Supergrid nxp mismatch: file has "// &
                   to_string(sg_nxp)//" but expected "// &
                   to_string(2*ni + 1)//" (2*nx_phys+1)", ierr, OCEAN_STATUS_ERR_IO)
         return
      end if
      if (sg_nyp /= 2*nj + 1) then
         call nc_close(ncid)
         call fail("Supergrid nyp mismatch: file has "// &
                   to_string(sg_nyp)//" but expected "// &
                   to_string(2*nj + 1)//" (2*ny_phys+1)", ierr, OCEAN_STATUS_ERR_IO)
         return
      end if
      if (sg_nx /= 2*ni) then
         call nc_close(ncid)
         call fail("Supergrid nx mismatch: file has "// &
                   to_string(sg_nx)//" but expected "// &
                   to_string(2*ni)//" (2*nx_phys)", ierr, OCEAN_STATUS_ERR_IO)
         return
      end if
      if (sg_ny /= 2*nj) then
         call nc_close(ncid)
         call fail("Supergrid ny mismatch: file has "// &
                   to_string(sg_ny)//" but expected "// &
                   to_string(2*nj)//" (2*ny_phys)", ierr, OCEAN_STATUS_ERR_IO)
         return
      end if

      ! ---- Read supergrid arrays ----
      call nc_get_varid(ncid, "x", varid_x, ierr=local_ierr)
      if (.not. supergrid_io_ok(local_ierr, ierr, ncid)) return
      call nc_get_varid(ncid, "y", varid_y, ierr=local_ierr)
      if (.not. supergrid_io_ok(local_ierr, ierr, ncid)) return
      call nc_get_varid(ncid, "dx", varid_dx, ierr=local_ierr)
      if (.not. supergrid_io_ok(local_ierr, ierr, ncid)) return
      call nc_get_varid(ncid, "dy", varid_dy, ierr=local_ierr)
      if (.not. supergrid_io_ok(local_ierr, ierr, ncid)) return
      call nc_get_varid(ncid, "area", varid_area, ierr=local_ierr)
      if (.not. supergrid_io_ok(local_ierr, ierr, ncid)) return

      allocate (sg_x(sg_nxp, sg_nyp))
      allocate (sg_y(sg_nxp, sg_nyp))
      allocate (sg_dx(sg_nx, sg_nyp))
      allocate (sg_dy(sg_nxp, sg_ny))
      allocate (sg_area(sg_nx, sg_ny))

      call nc_get_var_2d(ncid, varid_x, sg_x, ierr=local_ierr)
      if (.not. supergrid_io_ok(local_ierr, ierr, ncid)) return
      call nc_get_var_2d(ncid, varid_y, sg_y, ierr=local_ierr)
      if (.not. supergrid_io_ok(local_ierr, ierr, ncid)) return
      call nc_get_var_2d(ncid, varid_dx, sg_dx, ierr=local_ierr)
      if (.not. supergrid_io_ok(local_ierr, ierr, ncid)) return
      call nc_get_var_2d(ncid, varid_dy, sg_dy, ierr=local_ierr)
      if (.not. supergrid_io_ok(local_ierr, ierr, ncid)) return
      call nc_get_var_2d(ncid, varid_area, sg_area, ierr=local_ierr)
      if (.not. supergrid_io_ok(local_ierr, ierr, ncid)) return
      call nc_close(ncid)

      ! Assemble model metrics from the supergrid node/segment arrays via
      ! the shared even/odd index-sum logic (also used by the tripolar
      ! analytic generator).
      call metrics_assemble_from_supergrid_arrays(this, grid, &
                                                  sg_x, sg_y, sg_dx, sg_dy, sg_area)

      deallocate (sg_x, sg_y, sg_dx, sg_dy, sg_area)
      call logger%info("Supergrid metrics loaded for "// &
                       to_string(ni)//"x"//to_string(nj)//" physical grid")
      if (present(ierr)) ierr = OCEAN_STATUS_OK
#else
      call fail("grid_config='supergrid' requires RDB_ENABLE_NETCDF=ON "// &
                "(the supergrid reader needs the NetCDF-backed rdb_io_netcdf)", ierr, OCEAN_STATUS_ERR_IO)
      return
#endif
   end subroutine metrics_fill_from_supergrid

#ifndef RDB_NO_NETCDF
   function supergrid_io_ok(local_ierr, ierr, ncid) result(ok)
      !! Translate a raw `nc_check`-style status (0 = ok) from one of the
      !! `nc_*` reader calls in `metrics_fill_from_supergrid` into the
      !! caller's `ierr` contract: `.true.` on success; on failure,
      !! returns `.false.` with `ierr = OCEAN_STATUS_ERR_IO` when `ierr`
      !! is present (closing `ncid` first, when given, so a mid-read
      !! failure does not leak the file handle), or `error stop`s with
      !! the SAME generic text `nc_check` itself would have used had the
      !! caller's `ierr` never been threaded through — this is what keeps
      !! the legacy (no `ierr`) behaviour byte-identical while unblocking
      !! the `ierr`-present return path (F1/F2 of the P0.1 review).
      integer, intent(in) :: local_ierr
      integer, intent(out), optional :: ierr
      integer, intent(in), optional :: ncid
      logical :: ok

      integer :: discard_ierr

      ok = (local_ierr == 0)
      if (ok) return

      if (present(ierr)) then
         if (present(ncid)) call nc_close(ncid, ierr=discard_ierr)
         ierr = OCEAN_STATUS_ERR_IO
         return
      end if

      error stop "NetCDF operation failed"
   end function supergrid_io_ok
#endif

   ! =================================================================
   ! Shared supergrid-array assembler (NetCDF reader + tripolar both use it)
   ! =================================================================

   subroutine metrics_assemble_from_supergrid_arrays(this, grid, &
                                                     sg_x, sg_y, sg_dx, sg_dy, sg_area)
      !! Fill all model metric arrays from an in-memory MOM6-style
      !! supergrid (2x-refined corner geography + edge segments +
      !! sub-cell areas), using the even/odd index sums.  This is the
      !! battle-tested assembly path the NetCDF reader used inline; the
      !! tripolar generator builds the supergrid analytically and feeds it
      !! here so tripolar metrics flow through identical index logic.
      !!
      !! Supergrid index convention (1-based, node (1,1) = SW corner of
      !! the physical domain; 2*ni+1 nodes in i, 2*nj+1 in j):
      !!     node (2i-1, 2j-1): SW corner of T(i,j)  -- ODD/ODD = Bu
      !!     node (2i,   2j  ): T-cell centre          -- EVEN/EVEN = T
      !!     node (2i-1, 2j  ): west  face of T(i,j)    -- ODD/EVEN = Cu
      !!     node (2i,   2j-1): south face of T(i,j)    -- EVEN/ODD = Cv
      !!   sg_dx(m,n): along-i segment node (m,n)->(m+1,n), shape (2ni,2nj+1).
      !!   sg_dy(m,n): along-j segment node (m,n)->(m,n+1), shape (2ni+1,2nj).
      !!   sg_area(m,n): sub-cell area SW at node (m,n), shape (2ni,2nj).
      !! Boundary (Cu i=1/ni+1, Cv j=1/nj+1, Bu edges) by extrapolation.
      !! Ghost rows by constant extrapolation of the nearest physical value
      !! (periodic / fold ghost fill is the M4c exchange job).
      type(ocean_metrics_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: sg_x(:, :), sg_y(:, :)
      real(wp), intent(in) :: sg_dx(:, :), sg_dy(:, :), sg_area(:, :)

      integer :: ni, nj, ng, sg_nxp, sg_nyp
      integer :: i, j, si, sj, si1, sj1

      ni = grid%nx_phys
      nj = grid%ny_phys
      ng = grid%nghost
      sg_nxp = 2*ni + 1
      sg_nyp = 2*nj + 1

      ! ---- T metrics (all physical i,j) ----
      do j = 1, nj
         sj = 2*j       ! T-centre sg j-index (even)
         sj1 = 2*j - 1  ! lower corner / SW sg j-index (odd)
         do i = 1, ni
            si = 2*i       ! T-centre sg i-index (even)
            si1 = 2*i - 1  ! left corner / SW sg i-index (odd)
            this%geolatT(ng + i, ng + j) = sg_y(si, sj)
            this%geolonT(ng + i, ng + j) = sg_x(si, sj)
            this%dxT(ng + i, ng + j) = sg_dx(si1, sj) + sg_dx(si, sj)
            this%dyT(ng + i, ng + j) = sg_dy(si, sj1) + sg_dy(si, sj)
            this%areaT(ng + i, ng + j) = sg_area(si1, sj1) + sg_area(si, sj1) + &
                                         sg_area(si1, sj) + sg_area(si, sj)
         end do
      end do

      ! ---- Cu metrics (u-face: i ∈ [1,ni+1], j ∈ [1,nj]) ----
      do j = 1, nj
         sj = 2*j
         sj1 = 2*j - 1
         do i = 1, ni + 1
            si = min(max(2*i - 1, 1), sg_nxp)  ! col of Cu face node
            if (i >= 2 .and. i <= ni) then
               this%dxCu(ng + i, ng + j) = sg_dx(2*i - 2, sj) + sg_dx(2*i - 1, sj)
            else
               this%dxCu(ng + i, ng + j) = 0.0_wp
            end if
            this%dyCu(ng + i, ng + j) = sg_dy(si, sj1) + sg_dy(si, sj)
            this%areaCu(ng + i, ng + j) = this%dxCu(ng + i, ng + j)* &
                                          this%dyCu(ng + i, ng + j)
            this%dy_cu(ng + i, ng + j) = this%dyCu(ng + i, ng + j)
         end do
      end do
      do j = 1, nj
         this%dxCu(ng + 1, ng + j) = this%dxCu(ng + 2, ng + j)
         this%areaCu(ng + 1, ng + j) = this%dxCu(ng + 1, ng + j)* &
                                       this%dyCu(ng + 1, ng + j)
         this%dxCu(ng + ni + 1, ng + j) = this%dxCu(ng + ni, ng + j)
         this%areaCu(ng + ni + 1, ng + j) = this%dxCu(ng + ni + 1, ng + j)* &
                                            this%dyCu(ng + ni + 1, ng + j)
      end do

      ! ---- Cv metrics (v-face: i ∈ [1,ni], j ∈ [1,nj+1]) ----
      do i = 1, ni
         si = 2*i
         si1 = 2*i - 1
         do j = 1, nj + 1
            sj = min(max(2*j - 1, 1), sg_nyp)  ! row of Cv face node
            this%dxCv(ng + i, ng + j) = sg_dx(si1, sj) + sg_dx(si, sj)
            if (j >= 2 .and. j <= nj) then
               this%dyCv(ng + i, ng + j) = sg_dy(si, 2*j - 2) + sg_dy(si, 2*j - 1)
            else
               this%dyCv(ng + i, ng + j) = 0.0_wp
            end if
            this%areaCv(ng + i, ng + j) = this%dxCv(ng + i, ng + j)* &
                                          this%dyCv(ng + i, ng + j)
            this%dx_cv(ng + i, ng + j) = this%dxCv(ng + i, ng + j)
         end do
      end do
      do i = 1, ni
         this%dyCv(ng + i, ng + 1) = this%dyCv(ng + i, ng + 2)
         this%areaCv(ng + i, ng + 1) = this%dxCv(ng + i, ng + 1)* &
                                       this%dyCv(ng + i, ng + 1)
         this%dx_cv(ng + i, ng + 1) = this%dxCv(ng + i, ng + 1)
         this%dyCv(ng + i, ng + nj + 1) = this%dyCv(ng + i, ng + nj)
         this%areaCv(ng + i, ng + nj + 1) = this%dxCv(ng + i, ng + nj + 1)* &
                                            this%dyCv(ng + i, ng + nj + 1)
         this%dx_cv(ng + i, ng + nj + 1) = this%dxCv(ng + i, ng + nj + 1)
      end do

      ! ---- Bu metrics (corner: i ∈ [1,ni+1], j ∈ [1,nj+1]) ----
      do j = 1, nj + 1
         sj = 2*j - 1   ! ODD: sg j-index of Bu corner
         do i = 1, ni + 1
            si = 2*i - 1  ! ODD: sg i-index of Bu corner
            this%geolatBu(ng + i, ng + j) = sg_y(si, sj)
            this%geolonBu(ng + i, ng + j) = sg_x(si, sj)
            if (i >= 2 .and. i <= ni) then
               this%dxBu(ng + i, ng + j) = sg_dx(2*i - 2, sj) + sg_dx(2*i - 1, sj)
            else
               this%dxBu(ng + i, ng + j) = 0.0_wp
            end if
            if (j >= 2 .and. j <= nj) then
               this%dyBu(ng + i, ng + j) = sg_dy(si, 2*j - 2) + sg_dy(si, 2*j - 1)
            else
               this%dyBu(ng + i, ng + j) = 0.0_wp
            end if
            this%areaBu(ng + i, ng + j) = this%dxBu(ng + i, ng + j)* &
                                          this%dyBu(ng + i, ng + j)
         end do
      end do
      do j = 1, nj + 1
         this%dxBu(ng + 1, ng + j) = this%dxBu(ng + 2, ng + j)
         this%areaBu(ng + 1, ng + j) = this%areaBu(ng + 2, ng + j)
         this%dxBu(ng + ni + 1, ng + j) = this%dxBu(ng + ni, ng + j)
         this%areaBu(ng + ni + 1, ng + j) = this%areaBu(ng + ni, ng + j)
      end do
      do i = 1, ni + 1
         this%dyBu(ng + i, ng + 1) = this%dyBu(ng + i, ng + 2)
         if (i > 1) then
            this%areaBu(ng + i, ng + 1) = this%dxBu(ng + i, ng + 1)* &
                                          this%dyBu(ng + i, ng + 1)
         end if
         this%dyBu(ng + i, ng + nj + 1) = this%dyBu(ng + i, ng + nj)
         this%areaBu(ng + i, ng + nj + 1) = this%dxBu(ng + i, ng + nj + 1)* &
                                            this%dyBu(ng + i, ng + nj + 1)
      end do
      this%areaBu(ng + 1, ng + 1) = this%areaBu(ng + 2, ng + 2)

      ! ---- Ghost extrapolation: constant copy from nearest physical cell ----
      call supergrid_ghost_fill_2d(this%dxT, grid)
      call supergrid_ghost_fill_2d(this%dyT, grid)
      call supergrid_ghost_fill_2d(this%areaT, grid)
      call supergrid_ghost_fill_2d(this%geolatT, grid)
      call supergrid_ghost_fill_2d(this%geolonT, grid)
      call supergrid_ghost_fill_cu(this%dxCu, grid)
      call supergrid_ghost_fill_cu(this%dyCu, grid)
      call supergrid_ghost_fill_cu(this%areaCu, grid)
      call supergrid_ghost_fill_cu(this%dy_cu, grid)
      call supergrid_ghost_fill_cv(this%dxCv, grid)
      call supergrid_ghost_fill_cv(this%dyCv, grid)
      call supergrid_ghost_fill_cv(this%areaCv, grid)
      call supergrid_ghost_fill_cv(this%dx_cv, grid)
      call supergrid_ghost_fill_bu(this%dxBu, grid)
      call supergrid_ghost_fill_bu(this%dyBu, grid)
      call supergrid_ghost_fill_bu(this%areaBu, grid)
      call supergrid_ghost_fill_bu(this%geolatBu, grid)
      call supergrid_ghost_fill_bu(this%geolonBu, grid)
   end subroutine metrics_assemble_from_supergrid_arrays

   ! =================================================================
   ! Generator: analytic TRIPOLAR (Murray 1996)
   ! =================================================================

   subroutine metrics_fill_tripolar(this, grid, lon_west, lat_south, &
                                    dlon_deg, dlat_deg, rad_earth, phi_join, lon_pole)
      !! Fill every metric array for an analytic TRIPOLAR grid (Murray
      !! 1996): ordinary lon-lat for cell-corner latitude <= `phi_join`,
      !! and a conformal bipolar Arctic cap above (two grid poles at
      !! `(phi_join, lon_pole)` and `(phi_join, lon_pole+180)`; see
      !! `rdb_ocean_bipolar`).  The construction generates an in-memory
      !! MOM6-style supergrid (2x-refined corner geography, great-circle
      !! edge lengths, sub-cell areas) from the analytic map, then feeds
      !! the SAME `metrics_assemble_from_supergrid_arrays` index-sum path
      !! the NetCDF reader uses — so tripolar metrics flow through the
      !! battle-tested supergrid assembly.  Call `metrics_finalize` after.
      !!
      !! Logical layout: the i-direction wraps the full 360 deg of
      !! pseudo-longitude (lon_west .. lon_west+360); the j-direction goes
      !! from lat_south up.  Corner Bu(i,j) = SW corner of T(i,j) sits at
      !! geographic lon `lon_west + (i-ng-1)*dlon` BELOW the join.  The cap
      !! starts where the corner latitude exceeds `phi_join`; inside it the
      !! row fraction `s = (lat_lonlat - phi_join)/(lat_top - phi_join)` in
      !! [0,1] drives the bipolar map (s=0 reproduces the join ring exactly
      !! -> C0 continuity; s=1 is the north-fold line).
      type(ocean_metrics_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: lon_west, lat_south, dlon_deg, dlat_deg
      real(wp), intent(in) :: rad_earth, phi_join, lon_pole

      integer :: ni, nj, sg_nxp, sg_nyp, m, n
      real(wp), allocatable :: sg_x(:, :), sg_y(:, :)
      real(wp), allocatable :: sg_dx(:, :), sg_dy(:, :), sg_area(:, :)
      real(wp) :: lat_top, dlam, dlat_sg

      ni = grid%nx_phys
      nj = grid%ny_phys
      sg_nxp = 2*ni + 1
      sg_nyp = 2*nj + 1

      ! Supergrid node spacing: half a model cell in each logical direction.
      ! The i-direction must wrap exactly 360 deg of pseudo-longitude over
      ! ni model cells, so dlam per supergrid i-step = 360/(2*ni).
      dlam = 360.0_wp/real(2*ni, wp)
      dlat_sg = dlat_deg*0.5_wp

      ! Top corner latitude of the lon-lat ladder (the cap is mapped between
      ! phi_join and lat_top).  Bu corner j runs 1..nj+1; supergrid corner
      ! row n = 2*nj+1 is the top corner = lat_south + nj*dlat.
      lat_top = lat_south + real(nj, wp)*dlat_deg

      allocate (sg_x(sg_nxp, sg_nyp), source=0.0_wp)
      allocate (sg_y(sg_nxp, sg_nyp), source=0.0_wp)
      allocate (sg_dx(2*ni, sg_nyp), source=0.0_wp)
      allocate (sg_dy(sg_nxp, 2*nj), source=0.0_wp)
      allocate (sg_area(2*ni, 2*nj), source=0.0_wp)

      ! ---- Supergrid node geography (lon-lat below join, bipolar above) ----
      ! Supergrid node (m,n), m=1..2ni+1, n=1..2nj+1.  Corner Bu(1,1) is at
      ! node (1,1) = (lon_west, lat_south).  A T-centre lies at even (m,n).
      ! lon-lat geographic at node (m,n):
      !   lon0 = lon_west + (m-1)*dlam
      !   lat0 = lat_south + (n-1)*dlat_sg
      do n = 1, sg_nyp
         do m = 1, sg_nxp
            call tripolar_node_latlon(m, n, lon_west, lat_south, dlam, dlat_sg, &
                                      phi_join, lat_top, lon_pole, &
                                      sg_y(m, n), sg_x(m, n))
         end do
      end do

      ! ---- Edge lengths by great-circle distance between adjacent nodes ----
      do n = 1, sg_nyp
         do m = 1, 2*ni
            sg_dx(m, n) = great_circle(rad_earth, sg_y(m, n), sg_x(m, n), &
                                       sg_y(m + 1, n), sg_x(m + 1, n))
         end do
      end do
      do n = 1, 2*nj
         do m = 1, sg_nxp
            sg_dy(m, n) = great_circle(rad_earth, sg_y(m, n), sg_x(m, n), &
                                       sg_y(m, n + 1), sg_x(m, n + 1))
         end do
      end do

      ! ---- Sub-cell areas: spherical quad from the four corner nodes ----
      do n = 1, 2*nj
         do m = 1, 2*ni
            sg_area(m, n) = spherical_quad_area(rad_earth, &
                                                sg_y(m, n), sg_x(m, n), &
                                                sg_y(m + 1, n), sg_x(m + 1, n), &
                                                sg_y(m + 1, n + 1), sg_x(m + 1, n + 1), &
                                                sg_y(m, n + 1), sg_x(m, n + 1))
         end do
      end do

      call metrics_assemble_from_supergrid_arrays(this, grid, &
                                                  sg_x, sg_y, sg_dx, sg_dy, sg_area)

      ! Replace the assembler's constant-extrapolation ghosts on the seam
      ! edges with the physically-correct fold (north) + periodic (east-west)
      ! values.  Every metric/geography array is a SCALAR under the fold
      ! reflection (lengths/areas invariant; geography reads the conjugate
      ! point's stored coordinate), so all fold ops use negate=.false.
      call metrics_fold_periodic_ghosts(this, grid)

      deallocate (sg_x, sg_y, sg_dx, sg_dy, sg_area)
   end subroutine metrics_fill_tripolar

   subroutine metrics_fold_periodic_ghosts(this, grid)
      !! Tripolar ghost-metric fill (M4c): periodic-x wrap of the
      !! east/west ghost columns + north-fold of the north ghost rows,
      !! for EVERY metric + geography array.  Replaces the constant
      !! extrapolation the supergrid assembler left on those edges.
      !!
      !! Composition (Appendix A): periodic-x FIRST so the fold reads the
      !! cyclically-wrapped corner columns.  All arrays fold as scalars
      !! (negate=.false.) — lengths/areas are reflection-invariant and
      !! geography reads the conjugate point's stored lat/lon.
      type(ocean_metrics_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      integer :: ng, ni, nj

      ng = grid%nghost
      ni = grid%nx_phys
      nj = grid%ny_phys

      ! ---- (1) Periodic-x wrap of east/west ghost columns ----
      call metrics_periodic_x_2d(this%dxT, grid)
      call metrics_periodic_x_2d(this%dyT, grid)
      call metrics_periodic_x_2d(this%areaT, grid)
      call metrics_periodic_x_2d(this%geolatT, grid)
      call metrics_periodic_x_2d(this%geolonT, grid)
      call metrics_periodic_x_cu(this%dxCu, grid)
      call metrics_periodic_x_cu(this%dyCu, grid)
      call metrics_periodic_x_cu(this%areaCu, grid)
      call metrics_periodic_x_cu(this%dy_cu, grid)
      call metrics_periodic_x_cv(this%dxCv, grid)
      call metrics_periodic_x_cv(this%dyCv, grid)
      call metrics_periodic_x_cv(this%areaCv, grid)
      call metrics_periodic_x_cv(this%dx_cv, grid)
      call metrics_periodic_x_bu(this%dxBu, grid)
      call metrics_periodic_x_bu(this%dyBu, grid)
      call metrics_periodic_x_bu(this%areaBu, grid)
      call metrics_periodic_x_bu(this%geolatBu, grid)
      call metrics_periodic_x_bu(this%geolonBu, grid)

      ! ---- (2) North fold of the north ghost rows (scalars: negate=.false.) ----
      ! T-stagger (centre).
      call fold_north_centre(this%dxT, grid%nx_total, grid%ny_total, ni, nj, ng)
      call fold_north_centre(this%dyT, grid%nx_total, grid%ny_total, ni, nj, ng)
      call fold_north_centre(this%areaT, grid%nx_total, grid%ny_total, ni, nj, ng)
      call fold_north_centre(this%geolatT, grid%nx_total, grid%ny_total, ni, nj, ng)
      call fold_north_centre(this%geolonT, grid%nx_total, grid%ny_total, ni, nj, ng)
      ! u-stagger (Cu) — scalar copy.  fold_north_u_face NEGATES (it is
      ! built for the vector u-component), so metrics use a local scalar-copy
      ! variant with the same (nx+1,ny) index map.
      call metrics_fold_north_cu_scalar(this%dxCu, grid)
      call metrics_fold_north_cu_scalar(this%dyCu, grid)
      call metrics_fold_north_cu_scalar(this%areaCu, grid)
      call metrics_fold_north_cu_scalar(this%dy_cu, grid)
      ! v-stagger (Cv) — scalar copy on the (nx,ny+1) extent.
      call metrics_fold_north_cv_scalar(this%dxCv, grid)
      call metrics_fold_north_cv_scalar(this%dyCv, grid)
      call metrics_fold_north_cv_scalar(this%areaCv, grid)
      call metrics_fold_north_cv_scalar(this%dx_cv, grid)
      ! corner-stagger (Bu) — scalar copy via the corner op with negate=.false.
      call fold_north_corner(this%dxBu, grid%nx_total + 1, grid%ny_total + 1, &
                             ni, nj, ng, negate=.false.)
      call fold_north_corner(this%dyBu, grid%nx_total + 1, grid%ny_total + 1, &
                             ni, nj, ng, negate=.false.)
      call fold_north_corner(this%areaBu, grid%nx_total + 1, grid%ny_total + 1, &
                             ni, nj, ng, negate=.false.)
      call fold_north_corner(this%geolatBu, grid%nx_total + 1, grid%ny_total + 1, &
                             ni, nj, ng, negate=.false.)
      call fold_north_corner(this%geolonBu, grid%nx_total + 1, grid%ny_total + 1, &
                             ni, nj, ng, negate=.false.)
   end subroutine metrics_fold_periodic_ghosts

   ! Scalar north-fold for Cu-shaped (nx+1,ny) metric arrays.  The reusable
   ! fold_north_u_face NEGATES (vector); metrics are scalars, so we
   ! replicate the same index map with a straight copy.
   subroutine metrics_fold_north_cu_scalar(arr, grid)
      real(wp), intent(inout) :: arr(:, :)
      type(hgrid_t), intent(in) :: grid
      integer :: ng, ni, nj, fsum, jsum, j_lo, i, j
      ng = grid%nghost
      ni = grid%nx_phys
      nj = grid%ny_phys
      fsum = 2*ng + ni + 2
      jsum = 2*ng + 2*nj + 1
      j_lo = ng + nj + 1
      do j = j_lo, size(arr, 2)
         do i = 1, size(arr, 1)
            arr(i, j) = arr(fsum - i, jsum - j)
         end do
      end do
   end subroutine metrics_fold_north_cu_scalar

   ! Scalar north-fold for Cv-shaped (nx,ny+1) metric arrays.  Like
   ! fold_north_v_face but a COPY (scalar) — and only the halo rows above
   ! the fold line (the on-line row j=j_fold is physical, not a ghost, for
   ! metrics, and already holds the correct value from assembly).
   subroutine metrics_fold_north_cv_scalar(arr, grid)
      real(wp), intent(inout) :: arr(:, :)
      type(hgrid_t), intent(in) :: grid
      integer :: ng, ni, nj, isum, jsum, j_fold, i, j
      ng = grid%nghost
      ni = grid%nx_phys
      nj = grid%ny_phys
      isum = 2*ng + ni + 1
      jsum = 2*ng + 2*nj
      j_fold = ng + nj
      do j = j_fold + 1, size(arr, 2)
         do i = 1, size(arr, 1)
            arr(i, j) = arr(isum - i, jsum - j)
         end do
      end do
   end subroutine metrics_fold_north_cv_scalar

   ! ---- Periodic-x ghost-column fill helpers (tripolar east/west) ----
   subroutine metrics_periodic_x_2d(arr, grid)
      !! West/east ghost columns of a T-array (nx_total,ny_total) by
      !! periodic wrap (column i <= ng ← i+ni; i > ng+ni ← i-ni).
      real(wp), intent(inout) :: arr(:, :)
      type(hgrid_t), intent(in) :: grid
      integer :: ng, ni, i, j
      ng = grid%nghost
      ni = grid%nx_phys
      do j = 1, size(arr, 2)
         do i = 1, ng
            arr(i, j) = arr(i + ni, j)
            arr(ng + ni + i, j) = arr(ng + i, j)
         end do
      end do
   end subroutine metrics_periodic_x_2d

   subroutine metrics_periodic_x_cu(arr, grid)
      !! Cu-array (nx_total+1,ny_total): faces 1..ni+1 physical at i=ng+1..ng+ni+1.
      real(wp), intent(inout) :: arr(:, :)
      type(hgrid_t), intent(in) :: grid
      integer :: ng, ni, i, j
      ng = grid%nghost
      ni = grid%nx_phys
      do j = 1, size(arr, 2)
         do i = 1, ng
            arr(i, j) = arr(i + ni, j)
            arr(ng + ni + 1 + i, j) = arr(ng + 1 + i, j)
         end do
      end do
   end subroutine metrics_periodic_x_cu

   subroutine metrics_periodic_x_cv(arr, grid)
      !! Cv-array (nx_total,ny_total+1): centre-type in x, same as T.
      real(wp), intent(inout) :: arr(:, :)
      type(hgrid_t), intent(in) :: grid
      integer :: ng, ni, i, j
      ng = grid%nghost
      ni = grid%nx_phys
      do j = 1, size(arr, 2)
         do i = 1, ng
            arr(i, j) = arr(i + ni, j)
            arr(ng + ni + i, j) = arr(ng + i, j)
         end do
      end do
   end subroutine metrics_periodic_x_cv

   subroutine metrics_periodic_x_bu(arr, grid)
      !! Bu-array (nx_total+1,ny_total+1): face-type in x, same as Cu.
      real(wp), intent(inout) :: arr(:, :)
      type(hgrid_t), intent(in) :: grid
      integer :: ng, ni, i, j
      ng = grid%nghost
      ni = grid%nx_phys
      do j = 1, size(arr, 2)
         do i = 1, ng
            arr(i, j) = arr(i + ni, j)
            arr(ng + ni + 1 + i, j) = arr(ng + 1 + i, j)
         end do
      end do
   end subroutine metrics_periodic_x_bu

   pure subroutine tripolar_node_latlon(m, n, lon_west, lat_south, dlam, dlat_sg, &
                                        phi_join, lat_top, lon_pole, lat, lon)
      !! Geographic (lat, lon) of supergrid node (m,n).  Below the join
      !! (lon-lat corner latitude <= phi_join) it is plain lon-lat; above,
      !! the bipolar cap map (s = fraction of the cap row span).
      integer, intent(in) :: m, n
      real(wp), intent(in) :: lon_west, lat_south, dlam, dlat_sg
      real(wp), intent(in) :: phi_join, lat_top, lon_pole
      real(wp), intent(out) :: lat, lon
      real(wp) :: lon0, lat0, lam, s

      lon0 = lon_west + real(m - 1, wp)*dlam
      lat0 = lat_south + real(n - 1, wp)*dlat_sg

      if (lat0 <= phi_join .or. lat_top <= phi_join) then
         ! Below the join (or a degenerate cap): plain lon-lat.  Leave the
         ! longitude UNWRAPPED so geography matches metrics_fill_spherical
         ! exactly below the join (lon = lon_west + (m-1)*dlam).
         lat = lat0
         lon = lon0
      else
         ! Cap: pseudo-longitude is the i-coordinate's geographic lon; row
         ! fraction s maps the lon-lat ladder latitude into the bipolar cap.
         lam = lon0
         s = (lat0 - phi_join)/(lat_top - phi_join)
         if (s > 1.0_wp) s = 1.0_wp
         call bipolar_corner_latlon(lam, s, phi_join, lon_pole, lat, lon)
      end if
   end subroutine tripolar_node_latlon

   pure function great_circle(r, lat1, lon1, lat2, lon2) result(d)
      !! Great-circle distance (m) between two geographic points (deg),
      !! via the haversine formula (numerically stable for short arcs).
      real(wp), intent(in) :: r, lat1, lon1, lat2, lon2
      real(wp) :: d
      real(wp) :: p1, p2, dphi, dlam, a, h
      p1 = lat1*DEG2RAD
      p2 = lat2*DEG2RAD
      dphi = (lat2 - lat1)*DEG2RAD
      dlam = (lon2 - lon1)*DEG2RAD
      ! wrap dlam into [-pi, pi] so the antimeridian doesn't blow up
      dlam = modulo(dlam + 3.14159265358979323846_wp, 2.0_wp*3.14159265358979323846_wp) &
             - 3.14159265358979323846_wp
      h = sin(0.5_wp*dphi)**2 + cos(p1)*cos(p2)*sin(0.5_wp*dlam)**2
      a = 2.0_wp*atan2(sqrt(h), sqrt(max(0.0_wp, 1.0_wp - h)))
      d = r*a
   end function great_circle

   pure function spherical_quad_area(r, lat1, lon1, lat2, lon2, &
                                     lat3, lon3, lat4, lon4) result(area)
      !! Area (m^2) of a spherical quadrilateral with the four corners
      !! (1,2,3,4 counter-clockwise) given in degrees, via L'Huilier's
      !! theorem on the two triangles (1,2,3) and (1,3,4).
      real(wp), intent(in) :: r
      real(wp), intent(in) :: lat1, lon1, lat2, lon2, lat3, lon3, lat4, lon4
      real(wp) :: area
      area = (spherical_tri_area(r, lat1, lon1, lat2, lon2, lat3, lon3) &
              + spherical_tri_area(r, lat1, lon1, lat3, lon3, lat4, lon4))
   end function spherical_quad_area

   pure function spherical_tri_area(r, lat1, lon1, lat2, lon2, lat3, lon3) result(area)
      !! Area (m^2) of a spherical triangle (corners in degrees) via the
      !! spherical-excess form of L'Huilier's theorem.  Side lengths are
      !! angular (great-circle distance / r).
      real(wp), intent(in) :: r, lat1, lon1, lat2, lon2, lat3, lon3
      real(wp) :: area
      real(wp) :: a, b, c, sps, e, t
      a = great_circle(1.0_wp, lat2, lon2, lat3, lon3)
      b = great_circle(1.0_wp, lat1, lon1, lat3, lon3)
      c = great_circle(1.0_wp, lat1, lon1, lat2, lon2)
      sps = 0.5_wp*(a + b + c)
      ! tan(E/4) = sqrt(tan(s/2) tan((s-a)/2) tan((s-b)/2) tan((s-c)/2))
      t = tan(0.5_wp*sps)*tan(0.5_wp*(sps - a))* &
          tan(0.5_wp*(sps - b))*tan(0.5_wp*(sps - c))
      if (t <= 0.0_wp) then
         area = 0.0_wp
      else
         e = 4.0_wp*atan(sqrt(t))
         area = r*r*e
      end if
   end function spherical_tri_area

   ! -----------------------------------------------------------------
   ! Ghost-fill helpers (constant extrapolation from nearest physical)
   ! -----------------------------------------------------------------

   subroutine supergrid_ghost_fill_2d(arr, grid)
      !! Fill ghost rows/columns by constant extrapolation, for a T-point
      !! array `(nx_total, ny_total)`.  Interior = `[ng+1, ng+ni]` x
      !! `[ng+1, ng+nj]`.
      real(wp), intent(inout) :: arr(:, :)
      type(hgrid_t), intent(in) :: grid
      integer :: ng, nx, ny, i, j
      ng = grid%nghost
      nx = grid%nx_phys
      ny = grid%ny_phys
      do j = 1, size(arr, 2)
         do i = 1, ng
            arr(i, j) = arr(ng + 1, j)             ! west
            arr(ng + nx + i, j) = arr(ng + nx, j)  ! east
         end do
      end do
      do j = 1, ng
         do i = 1, size(arr, 1)
            arr(i, j) = arr(i, ng + 1)            ! south
            arr(i, ng + ny + j) = arr(i, ng + ny)  ! north
         end do
      end do
   end subroutine supergrid_ghost_fill_2d

   subroutine supergrid_ghost_fill_cu(arr, grid)
      !! Ghost fill for Cu arrays `(nx_total+1, ny_total)`.
      !! Physical i-range is `[ng+1, ng+ni+1]` (ni+1 faces), j-range
      !! `[ng+1, ng+nj]`.
      real(wp), intent(inout) :: arr(:, :)
      type(hgrid_t), intent(in) :: grid
      integer :: ng, nx, ny, i, j
      ng = grid%nghost
      nx = grid%nx_phys
      ny = grid%ny_phys
      ! west/east ghost columns
      do j = 1, size(arr, 2)
         do i = 1, ng
            arr(i, j) = arr(ng + 1, j)
            arr(ng + nx + 1 + i, j) = arr(ng + nx + 1, j)
         end do
      end do
      ! south/north ghost rows
      do j = 1, ng
         do i = 1, size(arr, 1)
            arr(i, j) = arr(i, ng + 1)
            arr(i, ng + ny + j) = arr(i, ng + ny)
         end do
      end do
   end subroutine supergrid_ghost_fill_cu

   subroutine supergrid_ghost_fill_cv(arr, grid)
      !! Ghost fill for Cv arrays `(nx_total, ny_total+1)`.
      !! Physical i-range `[ng+1, ng+ni]`, j-range `[ng+1, ng+nj+1]`.
      real(wp), intent(inout) :: arr(:, :)
      type(hgrid_t), intent(in) :: grid
      integer :: ng, nx, ny, i, j
      ng = grid%nghost
      nx = grid%nx_phys
      ny = grid%ny_phys
      ! west/east ghost columns
      do j = 1, size(arr, 2)
         do i = 1, ng
            arr(i, j) = arr(ng + 1, j)
            arr(ng + nx + i, j) = arr(ng + nx, j)
         end do
      end do
      ! south/north ghost rows
      do j = 1, ng
         do i = 1, size(arr, 1)
            arr(i, j) = arr(i, ng + 1)
            arr(i, ng + ny + 1 + j) = arr(i, ng + ny + 1)
         end do
      end do
   end subroutine supergrid_ghost_fill_cv

   subroutine supergrid_ghost_fill_bu(arr, grid)
      !! Ghost fill for Bu arrays `(nx_total+1, ny_total+1)`.
      !! Physical i-range `[ng+1, ng+ni+1]`, j-range `[ng+1, ng+nj+1]`.
      real(wp), intent(inout) :: arr(:, :)
      type(hgrid_t), intent(in) :: grid
      integer :: ng, nx, ny, i, j
      ng = grid%nghost
      nx = grid%nx_phys
      ny = grid%ny_phys
      ! west/east ghost columns
      do j = 1, size(arr, 2)
         do i = 1, ng
            arr(i, j) = arr(ng + 1, j)
            arr(ng + nx + 1 + i, j) = arr(ng + nx + 1, j)
         end do
      end do
      ! south/north ghost rows
      do j = 1, ng
         do i = 1, size(arr, 1)
            arr(i, j) = arr(i, ng + 1)
            arr(i, ng + ny + 1 + j) = arr(i, ng + ny + 1)
         end do
      end do
   end subroutine supergrid_ghost_fill_bu

   ! =================================================================
   ! Coriolis fills (D7) — corner + centre, from one routine
   ! =================================================================

   subroutine metrics_fill_coriolis(this, scheme, f_0, beta, y_ref, omega, &
                                    grid, f_corner, f_centre)
      !! Fill a corner array AND a centre array with the Coriolis
      !! parameter, from one of two schemes (D7).  Does NOT touch any
      !! existing fill sites in coriolis_adv / EPBL / kappa-shear (that
      !! re-routing is M2d); this routine just exists + is tested.
      !!
      !!   `beta_plane` : f = f_0 + beta*(y - y_ref), y from the
      !!      GLOBAL CARTESIAN coordinate.  On an undecomposed grid
      !!      (`grid%j_offset_global == 0`) the corner fill is
      !!      BIT-IDENTICAL to `coriolis_adv_set_beta_plane`
      !!      (y = (j-1-nghost)*dy, raw f) and the centre fill is
      !!      BIT-IDENTICAL to the EPBL / kappa-shear `set_f_centre`
      !!      (y = (j-nghost-0.5)*dy, abs(f)).  Under MPI y-decomposition
      !!      `grid%j_offset_global` shifts the local index to the GLOBAL
      !!      row so each rank's beta-plane y is correct.
      !!   `planetary` : f = 2*omega*sin(geolat) at the respective
      !!      stagger — uses `this%geolatBu` (corner) and `this%geolatT`
      !!      (centre), so a spherical generator must have run first.
      !!
      !! `f_corner` is shaped `(nx+1,ny+1)` (like the metric corner
      !! arrays / `coriolis_adv%f_corner`); `f_centre` is `(nx,ny)`.
      type(ocean_metrics_t), intent(in) :: this
      integer, intent(in) :: scheme
      real(wp), intent(in) :: f_0, beta, y_ref, omega
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(out) :: f_corner(:, :)
         !! Coriolis at C-grid corners (1/s).
      real(wp), intent(out) :: f_centre(:, :)
         !! |Coriolis| at cell centres (1/s).

      integer :: i, j, ng
      real(wp) :: y

      ng = grid%nghost

      select case (scheme)
      case (CORIOLIS_SCHEME_PLANETARY)
         ! 2*omega*sin(geolat) at the respective stagger.
         do j = 1, size(f_corner, 2)
            do i = 1, size(f_corner, 1)
               f_corner(i, j) = 2.0_wp*omega*sin(this%geolatBu(i, j)*DEG2RAD)
            end do
         end do
         do j = 1, size(f_centre, 2)
            do i = 1, size(f_centre, 1)
               f_centre(i, j) = abs(2.0_wp*omega*sin(this%geolatT(i, j)*DEG2RAD))
            end do
         end do
      case default   ! CORIOLIS_SCHEME_BETA_PLANE
         ! Corner: BIT-IDENTICAL to coriolis_adv_set_beta_plane on an
         ! undecomposed grid (j_offset_global == 0).  Under MPI y-split
         ! the global row index is j + j_offset_global, giving the correct
         ! physical y coordinate on every rank.
         do j = 1, size(f_corner, 2)
            y = real(j + grid%j_offset_global - 1 - ng, wp)*grid%dy
            do i = 1, size(f_corner, 1)
               f_corner(i, j) = f_0 + beta*(y - y_ref)
            end do
         end do
         ! Centre: BIT-IDENTICAL to EPBL / kappa-shear set_f_centre on an
         ! undecomposed grid (j_offset_global == 0).
         do j = 1, size(f_centre, 2)
            y = (real(j + grid%j_offset_global - ng, wp) - 0.5_wp)*grid%dy
            do i = 1, size(f_centre, 1)
               f_centre(i, j) = abs(f_0 + beta*(y - y_ref))
            end do
         end do
      end select
   end subroutine metrics_fill_coriolis

   ! =================================================================
   ! Enum parsers (host-side, configure time)
   ! =================================================================

   pure function parse_grid_config(s) result(cfg_enum)
      !! Map a `&ocean_grid_nml grid_config` string onto its enum.
      !! Unknown -> cartesian (the schema enum already validates the
      !! set; this is the canonical-name dispatch).
      character(len=*), intent(in) :: s
      integer :: cfg_enum
      select case (trim(s))
      case ("spherical")
         cfg_enum = GRID_CONFIG_SPHERICAL
      case ("supergrid")
         cfg_enum = GRID_CONFIG_SUPERGRID
      case ("tripolar")
         cfg_enum = GRID_CONFIG_TRIPOLAR
      case default
         cfg_enum = GRID_CONFIG_CARTESIAN
      end select
   end function parse_grid_config

   pure function parse_coriolis_scheme(s) result(scheme_enum)
      !! Map a `&ocean_grid_nml coriolis_scheme` string onto its enum.
      character(len=*), intent(in) :: s
      integer :: scheme_enum
      select case (trim(s))
      case ("planetary")
         scheme_enum = CORIOLIS_SCHEME_PLANETARY
      case default
         scheme_enum = CORIOLIS_SCHEME_BETA_PLANE
      end select
   end function parse_coriolis_scheme

   pure function ocean_metrics_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the grid metrics slot
      !! (0 when unallocated). One arr_bytes term per array — add a
      !! term here when a new allocatable joins the type.
      class(ocean_metrics_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = arr_bytes(this%dxT) &
               + arr_bytes(this%dyT) &
               + arr_bytes(this%dxCu) &
               + arr_bytes(this%dyCu) &
               + arr_bytes(this%dxCv) &
               + arr_bytes(this%dyCv) &
               + arr_bytes(this%dxBu) &
               + arr_bytes(this%dyBu) &
               + arr_bytes(this%dy_cu) &
               + arr_bytes(this%dx_cv) &
               + arr_bytes(this%dy_cu_bt) &
               + arr_bytes(this%dx_cv_bt) &
               + arr_bytes(this%por_bed) &
               + arr_bytes(this%por_dmin_u) &
               + arr_bytes(this%por_dmax_u) &
               + arr_bytes(this%por_davg_u) &
               + arr_bytes(this%por_dmin_v) &
               + arr_bytes(this%por_dmax_v) &
               + arr_bytes(this%por_davg_v) &
               + arr_bytes(this%por_face_area_u) &
               + arr_bytes(this%por_face_area_v) &
               + arr_bytes(this%open_u) &
               + arr_bytes(this%open_v) &
               + arr_bytes(this%z_draft) &
               + arr_bytes(this%cover_frac) &
               + arr_bytes(this%p_ice_ref) &
               + arr_bytes(this%wet_T) &
               + arr_bytes(this%wet_u) &
               + arr_bytes(this%wet_v) &
               + arr_bytes(this%wet_q) &
               + arr_bytes(this%areaT) &
               + arr_bytes(this%areaCu) &
               + arr_bytes(this%areaCv) &
               + arr_bytes(this%areaBu) &
               + arr_bytes(this%idxT) &
               + arr_bytes(this%idyT) &
               + arr_bytes(this%idxCu) &
               + arr_bytes(this%idyCu) &
               + arr_bytes(this%idxCv) &
               + arr_bytes(this%idyCv) &
               + arr_bytes(this%iareaT) &
               + arr_bytes(this%iareaBu) &
               + arr_bytes(this%iareaCu) &
               + arr_bytes(this%iareaCv) &
               + arr_bytes(this%geolatT) &
               + arr_bytes(this%geolonT) &
               + arr_bytes(this%geolatBu) &
               + arr_bytes(this%geolonBu) &
               + arr_bytes(this%dy_dxT) &
               + arr_bytes(this%dx_dyT) &
               + arr_bytes(this%dy_dxBu) &
               + arr_bytes(this%dx_dyBu) &
               + arr_bytes(this%dx2h) &
               + arr_bytes(this%dy2h) &
               + arr_bytes(this%dx2q) &
               + arr_bytes(this%dy2q)
   end function ocean_metrics_bytes

end module rdb_ocean_metrics
