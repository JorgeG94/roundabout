!! Ocean vertical-coordinate + ALE remap state.
module rdb_ocean_vcoord
   !! Holds the vertical-coordinate configuration and the per-step
   !! target grid that the ALE remap step relamps `multilayer.h_layer`
   !! + every `multilayer.tracers(t)%hTr` onto.  The coastal path has
   !! `VCOORD_SIGMA`, `VCOORD_ZSIGMA`, `VCOORD_ZSTAR`, `VCOORD_ZSTAR_SIGMA`,
   !! `VCOORD_ZSTAR_FULL` plus `rdb_remap_column`; this module is the
   !! C-grid counterpart for the ocean dynamical core.
   !!
   !! Phase 5g build-out status (this commit, Layer 2):
   !!
   !!   - `target_h(nx, ny, nz_ml)` is allocated on init + bound to the
   !!     device via `enter_data` / `exit_data`.
   !!   - `compute_target_h(this, total_h, eta)` populates `target_h`
   !!     per the current `coord_type`.  Working bodies:
   !!       VCOORD_EULERIAN_Z  — H · dsig(k)  (η ignored)
   !!       VCOORD_SIGMA       — (H + η) · dsig(k)
   !!       VCOORD_ZSTAR       — same formula in the barotropic limit
   !!       VCOORD_ZSIGMA      — smoothstep blend(sigma, fixed z-levels)
   !!       VCOORD_ZSTAR_SIGMA — smoothstep blend(sigma, z*-lite)
   !!       VCOORD_ZSTAR_FULL  — per-column z_ref + vanishing-layer floors
   !!   - `build_zref_full(this, h_bed)` populates `z_ref(:, :, 0:nz_ml)`
   !!     per column from local bathymetry.  Call once at init (or any
   !!     time the bathymetry changes); the per-step `compute_target_h`
   !!     then walks the cached `z_ref` table.
   !!   - Driver call site (remap kernel invocation between slow stages
   !!     in `ocean_dyn_step_split`) is NOT wired yet —
   !!     simulation state is untouched.  External callers can invoke
   !!     `compute_target_h` / `build_zref_full` to inspect the target
   !!     grid without advancing dynamics.
   !!
   !! Collaborator hand-off — remaining Phase 5g work (Layer 3):
   !!   1. Driver wiring in `ocean_dyn_step_split`: call
   !!      `compute_target_h(total_h_2d, dyn%bt_work%bt_eta)` then invoke the
   !!      remap kernel on `multilayer.h_layer` + each `tracers(t)%hTr`
   !!      + face velocities.  Recompute `bt_eta` from the new sum of
   !!      `h_layer - bt_H_ref`.
   !!   2. Face-velocity remap adapter — MOM6 reconstructs u at centres,
   !!      remaps, then projects back to faces with a divergence-free
   !!      correction.  Donor-cell on `(h·u)_face` is the simpler
   !!      fallback.
   !!   3. `VANISHING_LAYER_TOL`-gated CWC for `VCOORD_ZSTAR_FULL`,
   !!      lifted from the coastal multilayer kernels.
   !!
   !! See `docs/ROADMAP_OCEAN.md` Phase 5g for the full scope.
#ifdef LFORTRAN_PASSING
   use rdb_constants, only: wp, REMAP_PPM, H_VANISHED, &
                            VCOORD_LAGRANGIAN, VCOORD_EULERIAN_Z, &
                            VCOORD_SIGMA, VCOORD_ZSIGMA, VCOORD_ZSTAR, &
                            VCOORD_ZSTAR_FULL, VCOORD_ZSTAR_SIGMA, &
                            VCOORD_Z_FIXED, VCOORD_RHO, VCOORD_HYCOM
#else
   use rdb_constants, only: NZ_STACK_MAX, wp, REMAP_PPM, H_VANISHED, &
                            VCOORD_LAGRANGIAN, VCOORD_EULERIAN_Z, &
                            VCOORD_SIGMA, VCOORD_ZSIGMA, VCOORD_ZSTAR, &
                            VCOORD_ZSTAR_FULL, VCOORD_ZSTAR_SIGMA, &
                            VCOORD_Z_FIXED, VCOORD_RHO, VCOORD_HYCOM
#endif
   use rdb_grid, only: hgrid_t
   use rdb_vcoord, only: STRETCH_UNIFORM, STRETCH_LOG, parse_vcoord_type
   use rdb_eos, only: eos_t, eos_density_point
   use, intrinsic :: iso_fortran_env, only: int64
   use rdb_mem_report, only: arr_bytes
   implicit none
   private
#ifdef LFORTRAN_PASSING
   integer, parameter :: NZ_STACK_MAX = 64
      !! LFortran 0.64 workaround: module-local copy of the rdb_constants value
      !! (an imported parameter used as an explicit-shape dummy bound inside a
      !! PURE call becomes an impure getter under LFortran). Keep in sync (=64).
#endif

   public :: ocean_vcoord_t
   public :: parse_ocean_vcoord_type
   public :: invert_density_targets
   public :: VCOORD_EULERIAN_Z
   public :: VCOORD_LAGRANGIAN
   public :: VCOORD_Z_FIXED
   public :: VCOORD_RHO
   public :: VCOORD_HYCOM
   public :: STRETCH_UNIFORM, STRETCH_LOG

   ! ---- RHO (isopycnal) inversion parameters ----
   integer, parameter :: NR_ITERS = 8
      !! Fixed (GPU-uniform) Newton iteration budget for the
      !! density→depth inversion.  Unrolled, no data-dependent while.
   real(wp), parameter :: NR_TOL = 1.0e-12_wp
      !! Newton convergence tolerance — tested on |delta| AFTER xi += delta.
   real(wp), parameter :: NR_OFFSET = 1.0e-6_wp
      !! Out-of-range nudge applied only when the boundary gradient ≈ 0.

   type :: ocean_vcoord_t
      logical :: is_init = .false.
         !! True between `init` and `destroy`.  Prefer this to
         !! `allocated(...)` — tracks GPU device attachment too.

      ! ---- Active coordinate ----
      integer :: coord_type = VCOORD_EULERIAN_Z
         !! Selected vertical-coordinate variant.

      ! ---- Per-layer fractional thickness ----
      ! Sums to 1.0 across the column.  For `VCOORD_SIGMA` and
      ! `VCOORD_ZSTAR` this is the target σ stencil:
      !   target_h(i,j,k) = (H(i,j) + eta(i,j)) * dsig(k)
      ! For `VCOORD_ZSTAR_FULL` dsig is a fallback used when the per-
      ! column `z_ref` table is not populated.
      real(wp), allocatable :: dsig(:)
         !! Per-layer σ-fraction.  Sums to 1.0; size `nz_ml`.

      ! ---- Global z-level reference profile ----
      ! Reference z-interfaces in metres (positive-down), `z_ref_global(0)
      ! = 0` is the surface, `z_ref_global(nz_ml)` is the deepest
      ! reference interface.  Drives the z-level branch of `VCOORD_ZSIGMA`
      ! and the z*-lite branch of `VCOORD_ZSTAR_SIGMA`.  Default at init
      ! is uniform 0..1 (normalised) — the namelist parser populates it
      ! with absolute depths when those cases are activated.
      real(wp), allocatable :: z_ref_global(:)
         !! Global reference z-interfaces (m, positive-down), shape
         !! `0:nz_ml`.  Used by ZSIGMA / ZSTAR_SIGMA / ZSTAR.

      ! ---- Per-column target thickness ----
      ! Recomputed every outer step from (H, eta) per the coord_type
      ! case.  Consumed by the remap kernel (Layer 3) to advance
      ! `multilayer.h_layer` and every `tracers(t)%hTr`.
      real(wp), allocatable :: target_h(:, :, :)
         !! Target layer thickness (m), shape `(nx, ny, nz_ml)`.

      ! ---- Per-column z* reference profile ----
      ! Anchored to local bathymetry for `VCOORD_ZSTAR_FULL`.  Indexed
      ! from `k=0` (surface) to `k=nz_ml` (bed) — opposite of the
      ! bottom-up state convention so the surface anchor is at index 0
      ! (matches the coastal convention in `vcoord_target_dz_column_zstar_full`).
      ! Populated by `build_zref_full(h_bed)`; consumed by the
      ! ZSTAR_FULL branch of `compute_target_h`.
      real(wp), allocatable :: z_ref(:, :, :)
         !! Per-column z* reference (m), shape `(nx, ny, 0:nz_ml)`.

      ! ---- Isopycnal (VCOORD_RHO) target densities ----
      ! Monotone-increasing nominal interface potential densities
      ! (kg/m³) referenced to `rho_ref_pressure`.  Indexed `0:nz_ml`:
      ! `rho_target(0)` is the lightest (surface, maps to the k=nz
      ! interface in bottom-up state); `rho_target(nz_ml)` the densest
      ! (bed, k=1 interface).  The `compute_target_h_rho` inversion
      ! places interior interfaces where the reconstructed column
      ! density equals each interior target.  Sized + populated only
      ! when `coord_type == VCOORD_RHO`; ignored otherwise.
      real(wp), allocatable :: rho_target(:)
         !! Target interface potential densities (kg/m³), shape `0:nz_ml`.

      ! ---- ALE remap workspaces ----
      ! Allocated once at init + enter_data'd to device.  Previously
      ! these were allocated per call inside `ocean_apply_ale_remap_*`
      ! which generated thousands of host-allocated buffers per day
      ! that the device-side DCs in the remap step had to implicitly
      ! transfer back and forth.  Persistent device-mapped scratch
      ! eliminates that overhead entirely.
      real(wp), allocatable :: remap_total_h(:, :)
         !! Column-total h_layer scratch, shape `(nx, ny)`.
      real(wp), allocatable :: remap_h_ref(:, :)
         !! H reference (total_h − bt_eta) scratch, shape `(nx, ny)`.
      real(wp), allocatable :: remap_h_old(:, :, :)
         !! Snapshot of `h_layer` before the remap, shape `(nx, ny, nz)`.
      real(wp), allocatable :: remap_conc_t(:, :, :)
         !! Layer-mean T concentration scratch for the `VCOORD_RHO`
         !! density inversion, shape `(nx, ny, nz)`.  Built from
         !! `hTr / remap_h_old` (vanishing-layer-guarded) once per remap.
      real(wp), allocatable :: remap_conc_s(:, :, :)
         !! Layer-mean S concentration scratch for `VCOORD_RHO`, shape
         !! `(nx, ny, nz)`.

      ! ---- Tuning knobs ----
      integer :: remap_method = REMAP_PPM
         !! ALE remap reconstruction order (REMAP_PCM/PLM/PPM/PPM_H4/PQM).
         !! PQM falls back to PPM for nz < 5 (see `remap_column_pqm`).
      real(wp) :: zstar_h_surf_target = 5.0_wp
         !! Surface-layer thickness anchor for `VCOORD_ZSTAR_FULL` (m).
      real(wp) :: zstar_h_min = 1.0e-4_wp
         !! Bed-side vanishing-layer floor (m).  **Two contracts, picked by
         !! the coordinate family — see `rdb_vcoord :: vcoord_h_min_role`.**
         !!
         !! On the GEOMETRIC families (`VCOORD_ZSTAR_FULL`, `VCOORD_Z_FIXED`)
         !! this is the thickness handed to filler layers that lie BELOW the
         !! local bed.  They hold no water; the floor exists ONLY so
         !! `target_h` is never exactly zero and no h-dividing kernel can
         !! 1/0.  They are MEANT to be classified vanished downstream, so the
         !! default sits deliberately BELOW the D4 skip/merge marker
         !! `H_VANISHED = 1.5e-4` — not by accident, and not a floor in the
         !! `angstrom_h` sense (the D4 taxonomy forbids using `H_VANISHED`
         !! as a positivity floor).  Thinner is also better physics here:
         !! each filler interface carries the full topographic slope, so the
         !! spurious rest PGF transport it drives scales WITH the floor (the
         !! same argument that took `seed_h_layer_uniform_z_impl` off
         !! `2*H_VANISHED`).  `validate_config` warns on a value above
         !! `H_VANISHED` under these families, and refuses a non-positive one.
         !!
         !! On the DENSITY families (`VCOORD_RHO`, `VCOORD_HYCOM`) the
         !! collapsed layers are real layers the inversion squeezed shut
         !! anywhere in the column; they carry tracer mass, so
         !! `compute_target_h_rho_impl` inflates them to
         !! `max(zstar_h_min, 2*H_VANISHED)` to keep them above the remap
         !! drain.  There `zstar_h_min` is additionally the pre-compaction
         !! strip threshold, so a large value is meaningful rather than wrong.
      integer  :: zstar_n_surf = 0
         !! Number of fine near-surface layers for ZSTAR_FULL.  ≤ 0 =
         !! auto-pick (max(1, nz_ml/3)).
      integer  :: zstar_stretching = STRETCH_UNIFORM
         !! Stretching mode for ZSTAR_FULL.  `STRETCH_UNIFORM` (default)
         !! or `STRETCH_LOG` for a geometric near-surface fine zone.
      real(wp) :: zsigma_depth_transition = 200.0_wp
         !! Sigma → z* transition depth (m) for `VCOORD_ZSTAR_SIGMA`.
      real(wp) :: zsigma_blend_width = 100.0_wp
         !! Smoothstep blend width (m) above the transition depth.
      real(wp) :: rho_ref_pressure = 2.0e7_wp
         !! Reference pressure (Pa, default 2e7 = 2000 dbar) for the
         !! potential density that defines the `VCOORD_RHO` coordinate.
         !! A rdb convention (not MOM6-inherited).
      real(wp) :: z_fixed_h_ref = 0.0_wp
         !! Total reference depth (m) for `VCOORD_Z_FIXED`.  Layer
         !! interfaces sit at `z = k · h_ref / nz_ml` from the surface,
         !! same as MOM6's `COORD_CONFIG = "gprime"` with `MAXIMUM_DEPTH
         !! = h_ref`.  Driver writes from `cfg%ocean%topo%max_depth` at init.
         !! When 0 (default) the `compute_target_h` Z_FIXED branch falls
         !! back to a uniform `H · dsig(k)` target so the path stays
         !! sane in tests that don't explicitly set this knob.
      real(wp) :: regrid_time_scale = 0.0_wp
         !! Grid time-filter timescale τ (s) for the ALE regrid.  After
         !! `compute_target_h` builds the new target grid, the remap step
         !! relaxes the coordinate a fraction `dt/(τ+dt)` toward that
         !! target each outer step rather than jumping to it — damping the
         !! per-step grid-motion shock that drives the σ/z* PGE
         !! (White & Adcroft 2008, the grid time-filter).  Scalar on the
         !! type, reaches the device through the existing `copyin(this)`;
         !! no new device array.  Default `0.0` ⇒ `wtd = 1` ⇒ jump to
         !! target ⇒ bit-identical to the no-filter remap.
      logical :: remap_boundary_extrap = .false.
         !! Close the ALE remap's reconstruction at the two boundary cells
         !! (`k=1`, `k=nz`) with the linear-exact one-sided edge pair
         !! instead of the PCM flatten (MOM6 `BOUNDARY_EXTRAPOLATION`).
         !!
         !! The default closure makes PLM/PPM/PPM_H4/PQM first-order in
         !! exactly the two cells adjacent to the bed and the surface, so
         !! a column whose tracer is linear in z is remapped with an O(h)
         !! error there every thermo step.  Under a terrain-following
         !! coordinate over a slope that error differs between neighbouring
         !! columns, which is a horizontal density gradient, which is a
         !! spurious pressure-gradient force — and with rotation it feeds a
         !! growing grid mode trapped in those same layers (see
         !! `docs/CAPABILITIES_AND_LIMITATIONS.md`).  Scalar on the type,
         !! reaches the device through the existing `copyin(this)`; no new
         !! device array.  Default `.false.` ⇒ bit-identical.
      logical :: remap_vel_conserve_ke = .false.
         !! Enable the KE-conserving rescale of the remapped layer
         !! velocities.  After the per-face column remap (which already
         !! conserves `u·h`, i.e. momentum), rescale the BAROCLINIC
         !! velocity anomaly per column so column KE `Σ ½ h·u²` is
         !! preserved (Adcroft & Hallberg 2006 layer-velocity remap),
         !! capped at a 1.25× rescale factor.  The barotropic/depth-mean
         !! component is never touched (mode-split consistency).  Default
         !! `.false.` ⇒ velocities unchanged ⇒ bit-identical.

      ! ---- Cached extents (for kernel loops + sanity checks) ----
      integer :: nx_total = 0
         !! Total i-extent of `target_h` (incl. halos).
      integer :: ny_total = 0
         !! Total j-extent of `target_h` (incl. halos).
      integer :: nz_ml = 0
         !! Number of active layers.
   contains
      procedure, non_overridable :: init => ocean_vcoord_init
      procedure, non_overridable :: destroy => ocean_vcoord_destroy
      procedure, non_overridable :: enter_data => ocean_vcoord_enter_data
      procedure, non_overridable :: exit_data => ocean_vcoord_exit_data
      procedure, non_overridable :: compute_target_h => ocean_vcoord_compute_target_h
      procedure, non_overridable :: compute_target_h_rho => ocean_vcoord_compute_target_h_rho
      procedure, non_overridable :: build_zref_full => ocean_vcoord_build_zref_full
      procedure, non_overridable :: bytes => ocean_vcoord_bytes
   end type ocean_vcoord_t

contains

   subroutine ocean_vcoord_init(this, grid, nz_ml)
      !! Allocate every per-column array the slot owns: `dsig`,
      !! `z_ref_global`, `target_h`, `z_ref`.  All sized once at init —
      !! grids don't resize.  Host allocations only; `enter_data` ships
      !! them to the device.
      class(ocean_vcoord_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      integer, intent(in), optional :: nz_ml
      integer :: nz_local, k
      nz_local = 1
      if (present(nz_ml)) nz_local = nz_ml
      if (nz_local < 1) nz_local = 1
      this%nx_total = grid%nx_total
      this%ny_total = grid%ny_total
      this%nz_ml = nz_local
      allocate (this%dsig(nz_local))
      do k = 1, nz_local
         this%dsig(k) = 1.0_wp/real(nz_local, wp)
      end do
      ! Default reference z-interfaces: uniform 0..1 normalised.  The
      ! ZSIGMA / ZSTAR_SIGMA branches that consume this expect absolute
      ! metre values from the namelist parser; the normalised default
      ! is only useful for VCOORD_SIGMA / VCOORD_ZSTAR (which ignore it)
      ! and for unit tests that pre-populate before running.
      allocate (this%z_ref_global(0:nz_local))
      do k = 0, nz_local
         this%z_ref_global(k) = real(k, wp)/real(nz_local, wp)
      end do
      allocate (this%target_h(grid%nx_total, grid%ny_total, nz_local), source=0.0_wp)
      ! Per-column z* reference table — sized but not populated.  Callers
      ! invoke `build_zref_full(h_bed)` once at setup to fill it.  Until
      ! then, ZSTAR_FULL's `compute_target_h` walks a column of zeros,
      ! which the formula safely degrades to an all-h_min vanishing-
      ! layer column (the dry / shallow degenerate case).
      allocate (this%z_ref(grid%nx_total, grid%ny_total, 0:nz_local), source=0.0_wp)

      ! Isopycnal target densities — sized `0:nz_ml`, populated by the
      ! setup wiring only when `coord_type == VCOORD_RHO`.  Default is a
      ! benign monotone ramp (1020..1030 kg/m³) so the slot is always
      ! valid; the setup path overwrites it from the namelist.
      allocate (this%rho_target(0:nz_local))
      do k = 0, nz_local
         this%rho_target(k) = 1020.0_wp + 10.0_wp*real(k, wp)/real(nz_local, wp)
      end do

      ! ALE remap scratch — persistent so the remap step doesn't
      ! allocate fresh host buffers per call (the DCs there were
      ! reading device-resident `ms%h_layer` and writing into freshly-
      ! allocated host arrays, forcing per-iteration H↔D transfers
      ! that dominated wallclock in production runs).
      allocate (this%remap_total_h(grid%nx_total, grid%ny_total), source=0.0_wp)
      allocate (this%remap_h_ref(grid%nx_total, grid%ny_total), source=0.0_wp)
      allocate (this%remap_h_old(grid%nx_total, grid%ny_total, nz_local), source=0.0_wp)
      allocate (this%remap_conc_t(grid%nx_total, grid%ny_total, nz_local), source=0.0_wp)
      allocate (this%remap_conc_s(grid%nx_total, grid%ny_total, nz_local), source=0.0_wp)
      this%is_init = .true.
   end subroutine ocean_vcoord_init

   subroutine ocean_vcoord_destroy(this)
      class(ocean_vcoord_t), intent(inout) :: this
      this%is_init = .false.
      if (allocated(this%dsig)) deallocate (this%dsig)
      if (allocated(this%z_ref_global)) deallocate (this%z_ref_global)
      if (allocated(this%target_h)) deallocate (this%target_h)
      if (allocated(this%z_ref)) deallocate (this%z_ref)
      if (allocated(this%rho_target)) deallocate (this%rho_target)
      if (allocated(this%remap_total_h)) deallocate (this%remap_total_h)
      if (allocated(this%remap_h_ref)) deallocate (this%remap_h_ref)
      if (allocated(this%remap_h_old)) deallocate (this%remap_h_old)
      if (allocated(this%remap_conc_t)) deallocate (this%remap_conc_t)
      if (allocated(this%remap_conc_s)) deallocate (this%remap_conc_s)
      this%nx_total = 0
      this%ny_total = 0
      this%nz_ml = 0
   end subroutine ocean_vcoord_destroy

   subroutine ocean_vcoord_enter_data(this)
      !! Map every host allocatable onto the device.  Idempotent guard
      !! via `is_init`.
      class(ocean_vcoord_t), intent(inout) :: this
      select type (this)
      type is (ocean_vcoord_t)
         call ocean_vcoord_enter_data_impl(this)
      end select
   end subroutine ocean_vcoord_enter_data

   subroutine ocean_vcoord_enter_data_impl(this)
      type(ocean_vcoord_t), intent(inout) :: this
      if (.not. this%is_init) return
      !$acc enter data copyin(this%dsig, this%z_ref_global, this%target_h, this%z_ref)
      !$acc enter data copyin(this%rho_target)
      !$acc enter data copyin(this%remap_total_h, this%remap_h_ref, this%remap_h_old)
      !$acc enter data copyin(this%remap_conc_t, this%remap_conc_s)
   end subroutine ocean_vcoord_enter_data_impl

   subroutine ocean_vcoord_exit_data(this)
      class(ocean_vcoord_t), intent(inout) :: this
      select type (this)
      type is (ocean_vcoord_t)
         call ocean_vcoord_exit_data_impl(this)
      end select
   end subroutine ocean_vcoord_exit_data

   subroutine ocean_vcoord_exit_data_impl(this)
      type(ocean_vcoord_t), intent(inout) :: this
      if (.not. this%is_init) return
      !$acc exit data delete(this%remap_conc_t, this%remap_conc_s)
      !$acc exit data delete(this%remap_h_old, this%remap_h_ref, this%remap_total_h)
      !$acc exit data delete(this%rho_target)
      !$acc exit data delete(this%z_ref, this%target_h, this%z_ref_global, this%dsig)
   end subroutine ocean_vcoord_exit_data_impl

   pure subroutine ocean_vcoord_build_zref_full(this, h_bed)
      !! Populate `z_ref(i, j, 0:nz_ml)` per column from the local
      !! bathymetry `h_bed(i, j)`.  Mirrors `zstar_full_build_column`
      !! from `src/ALE/rdb_vcoord.F90` but as a 2D loop owned by this
      !! slot — keeps the coastal helper untouched while letting the
      !! ocean path own its z_ref lifecycle.
      !!
      !! Top-down indexing inside the column: `z_ref(:, :, 0) = 0` is
      !! the surface, `z_ref(:, :, nz_ml) = h_bed(:, :)` is the bed.
      !! `compute_target_h(VCOORD_ZSTAR_FULL)` walks this table and
      !! emits ROMS-ordered `target_h(:, :, 1..nz_ml)`.
      !!
      !! Degenerate columns (`h_bed ≤ 0`) get a column of zeros — the
      !! ZSTAR_FULL branch then produces an all-h_min vanishing-layer
      !! result which the downstream dry-cell guards already handle.
      !!
      !! Call sites: once at setup, then any time bathymetry changes
      !! (which today is "never" — the ocean path doesn't move the bed).
      class(ocean_vcoord_t), intent(inout) :: this
      real(wp), intent(in) :: h_bed(:, :)
         !! Bed depth at cell centres (m, positive-down).
      integer  :: n_surf_eff, n_coarse, nz, i, j, k
      real(wp) :: h_surf_eff, h_fine, h_coarse, dz_uniform, r, base, w
      real(wp) :: r_lo, r_hi, r_mid, f_mid, target_ratio
      integer  :: it

      if (.not. this%is_init) return
      nz = this%nz_ml

      h_surf_eff = this%zstar_h_surf_target
      if (this%zstar_n_surf <= 0) then
         n_surf_eff = max(1, nz/3)
      else
         n_surf_eff = max(1, min(nz - 1, this%zstar_n_surf))
      end if
      n_coarse = nz - n_surf_eff

      ! The per-column work is column-local — different (i, j) cells
      ! don't read each other.  The do-concurrent `local()` clause keeps
      ! the scalar scratch private per thread; the bisection branch in
      ! the LOG stretching case requires a sequential inner loop so we
      ! keep it as a non-concurrent block.
      ! One-time z_ref setup runs on the HOST (plain do, not do concurrent):
      ! it writes this%z_ref, a vcoord allocatable component that is not yet
      ! device-mapped at seed time (build_zref_full runs before
      ! ocean_state_enter_data).  A device kernel writing an un-present
      ! derived-type component faults under OpenMP-target offload (the
      ! component can't be implicitly mapped); stdpar tolerates it but this
      ! is one-time init, so host is correct and costs nothing.
      do j = 1, this%ny_total
         do i = 1, this%nx_total
         if (h_bed(i, j) <= 0.0_wp) then
            ! Degenerate / dry column — leave z_ref at zero.
            do k = 0, nz
               this%z_ref(i, j, k) = 0.0_wp
            end do
         else if (h_surf_eff <= 0.0_wp .or. nz == 1) then
            ! Uniform spacing fallback.
            this%z_ref(i, j, 0) = 0.0_wp
            do k = 1, nz
               this%z_ref(i, j, k) = h_bed(i, j)*real(k, wp)/real(nz, wp)
            end do
         else
            this%z_ref(i, j, 0) = 0.0_wp
            if (this%zstar_stretching == STRETCH_LOG &
                .and. n_surf_eff >= 2 .and. n_coarse >= 1) then
               ! Geometric fine zone: layer k thickness = h_surf · r^(k-1).
               ! Bisect for r so total = h_bed.
               target_ratio = h_bed(i, j)/h_surf_eff
               r_lo = 1.000001_wp
               r_hi = 10.0_wp
               do it = 1, 60
                  r_mid = 0.5_wp*(r_lo + r_hi)
                  f_mid = (r_mid**n_surf_eff - 1.0_wp)/(r_mid - 1.0_wp) &
                          + r_mid**(n_surf_eff - 1)*real(n_coarse, wp)
                  if (f_mid > target_ratio) then
                     r_hi = r_mid
                  else
                     r_lo = r_mid
                  end if
                  if (r_hi - r_lo < 1.0e-9_wp) exit
               end do
               r = 0.5_wp*(r_lo + r_hi)
               w = 1.0_wp
               base = 0.0_wp
               do k = 1, n_surf_eff
                  base = base + h_surf_eff*w
                  this%z_ref(i, j, k) = base
                  w = w*r
               end do
               dz_uniform = h_surf_eff*r**(n_surf_eff - 1)
               do k = n_surf_eff + 1, nz
                  this%z_ref(i, j, k) = this%z_ref(i, j, k - 1) + dz_uniform
               end do
            else
               ! Uniform fine zone + uniform coarse fill.
               h_fine = h_surf_eff*real(n_surf_eff, wp)
               if (h_fine >= h_bed(i, j)) then
                  ! Column too shallow to honour h_surf for all fine layers.
                  ! Match MOM6 isopycnal NK=2: surface absorbs the water,
                  ! deeper layers vanish.  Reserve `zstar_h_min` for each
                  ! vanishing layer so target_h ≥ h_min downstream — kernels
                  ! that divide by h_layer can't see exact zero (would NaN).
                  ! Fine layers fill top-down at h_surf_eff each; bottom-most
                  ! fine grabs the residual; coarse vanishes to h_min.
                  base = 0.0_wp
                  do k = 1, n_surf_eff
                     base = base + h_surf_eff
                     ! Leave room for h_min of every layer below this one.
                     if (base > h_bed(i, j) - real(nz - k, wp)*this%zstar_h_min) then
                        base = h_bed(i, j) - real(nz - k, wp)*this%zstar_h_min
                     end if
                     this%z_ref(i, j, k) = base
                  end do
                  do k = n_surf_eff + 1, nz
                     this%z_ref(i, j, k) = this%z_ref(i, j, k - 1) + this%zstar_h_min
                  end do
               else
                  h_coarse = h_bed(i, j) - h_fine
                  do k = 1, n_surf_eff
                     this%z_ref(i, j, k) = h_surf_eff*real(k, wp)
                  end do
                  if (n_coarse > 0) then
                     dz_uniform = h_coarse/real(n_coarse, wp)
                     do k = n_surf_eff + 1, nz
                        this%z_ref(i, j, k) = h_fine + dz_uniform*real(k - n_surf_eff, wp)
                     end do
                  end if
               end if
            end if
            ! Pin the bed interface to h_bed exactly (round-off guard).
            this%z_ref(i, j, nz) = h_bed(i, j)
         end if
         end do
      end do
   end subroutine ocean_vcoord_build_zref_full

   pure subroutine ocean_vcoord_compute_target_h(this, total_h, eta)
      !! Thin polymorphic wrapper.  A type-bound procedure's passed object
      !! must be `class(...)`, but mapping a polymorphic list item into a
      !! `target` / offload region is unspecified behaviour (gfortran
      !! `-Wopenmp`; see FORTRAN_STYLE.md) — so the `do concurrent` kernels
      !! live in the `type(ocean_vcoord_t)` `_impl` and this wrapper only
      !! resolves the concrete type.  `ocean_vcoord_t` is never extended,
      !! so the dynamic type is always the declared type.  Mirrors the
      !! enter_data/exit_data split.
      class(ocean_vcoord_t), intent(inout) :: this
      real(wp), intent(in) :: total_h(:, :)  ! assumed-shape-ok: thin TBP wrapper, no kernel
      real(wp), intent(in) :: eta(:, :)      ! assumed-shape-ok: thin TBP wrapper, no kernel
      select type (this)
      type is (ocean_vcoord_t)
         call ocean_vcoord_compute_target_h_impl(this, total_h, eta)
      end select
   end subroutine ocean_vcoord_compute_target_h

   pure subroutine ocean_vcoord_compute_target_h_impl(this, total_h, eta)
      !! Populate `target_h(i,j,k)` from the column-total depth (`H`,
      !! constant per column for the ocean path; coastal would pass
      !! the bathymetry) and the free-surface anomaly `η`.
      !!
      !! Cases:
      !!
      !!   VCOORD_EULERIAN_Z — `target_h(:,:,k) = H(:,:) · dsig(k)`.
      !!     The reference grid that `rdb_ocean_vertical_advection`
      !!     pins `h_layer` to via its cancellation trick.  No η
      !!     dependence — that's the whole point of the mode.
      !!
      !!   VCOORD_SIGMA      — `target_h(:,:,k) = (H + η) · dsig(k)`.
      !!     Pure terrain-following.
      !!
      !!   VCOORD_ZSTAR      — same formula as SIGMA in this barotropic
      !!     limit; distinction (z-anchored dsig) surfaces with a
      !!     non-uniform stencil.
      !!
      !!   VCOORD_ZSIGMA / VCOORD_ZSTAR_SIGMA — smoothstep blends, see
      !!     module head comment.
      !!
      !!   VCOORD_ZSTAR_FULL — walks the per-column `z_ref(:, :, 0:nz)`
      !!     table populated by `build_zref_full`.  Surface layer absorbs
      !!     η; subsurface layers keep their reference thicknesses
      !!     when the column is at or above reference depth.  When the
      !!     column is shallower than reference (η < 0), bed-side
      !!     layers vanish to `zstar_h_min` and the surface trim makes
      !!     sum(target_h) = H exactly.  Mirrors coastal
      !!     `vcoord_target_dz_column_zstar_full`.
      !!
      !! All cases preserve `sum_k target_h(i,j,k) = H + η` (or `= H`
      !! for EULERIAN_Z).  The remap kernel will rely on this.
      type(ocean_vcoord_t), intent(inout) :: this
      ! assumed-shape-ok: cadence-bounded (once per outer ALE step); flat impl
      ! behind the thin ocean_vcoord_compute_target_h class wrapper below.
      real(wp), intent(in) :: total_h(:, :)
         !! Column-total depth H(i, j) (m).
      real(wp), intent(in) :: eta(:, :)  ! assumed-shape-ok: same reason as total_h above
         !! Free-surface anomaly η(i, j) (m).
      integer :: i, j, k, nz
      real(wp) :: column_total, alpha, x, z_top_k, z_bot_k, dz_z, dz_sum, deficit
      real(wp) :: z_ref_nz_inv
      real(wp) :: h_bed_ref, eta_loc, H_eff, z_upper, z_lower, sum_dz
      real(wp) :: h_nominal, h_min, z_below_loc, z_above_nominal_loc

      if (.not. this%is_init) return
      ! Lagrangian / isopycnal: the target IS the current h_layer —
      ! nothing to compute, the ALE remap step is a no-op (see
      ! `ocean_apply_ale_remap_step`).  Return before touching
      ! `target_h` so the caller keeps the live `h_layer`.
      if (this%coord_type == VCOORD_LAGRANGIAN) return
      nz = this%nz_ml

      ! VCOORD_LAGRANGIAN: target is whatever `h_layer` already is.  The
      ! ALE remap is a no-op for this case (see `ocean_apply_ale_remap_step`),
      ! so `target_h` is never read; we can skip the compute entirely to
      ! avoid burning a kernel launch.
      if (this%coord_type == VCOORD_LAGRANGIAN) return

      ! Bind every derived-type component the offloaded do-concurrent loops
      ! below touch to a plain associate-name.  ifx's do-concurrent ->
      ! OpenMP-target lowering ICEs when a loop body references a derived-type
      ! allocatable component directly (it can't map the parent type into the
      ! target region); associate-names lower as plain array/scalar selectors,
      ! so the outliner never sees a `this%` inside the kernel.  Pure alias --
      ! no copy, device mapping of the components is unchanged.
      associate (target_h => this%target_h, dsig => this%dsig, &
                 z_ref_global => this%z_ref_global, z_ref => this%z_ref, &
                 nx_total => this%nx_total, ny_total => this%ny_total, &
                 zsigma_depth_transition => this%zsigma_depth_transition, &
                 zsigma_blend_width => this%zsigma_blend_width, &
                 zstar_h_min => this%zstar_h_min, &
                 z_fixed_h_ref => this%z_fixed_h_ref)
         select case (this%coord_type)

         case (VCOORD_EULERIAN_Z)
            do concurrent(k=1:nz, j=1:ny_total, i=1:nx_total)
               target_h(i, j, k) = total_h(i, j)*dsig(k)
            end do

         case (VCOORD_SIGMA, VCOORD_ZSTAR)
            do concurrent(k=1:nz, j=1:ny_total, i=1:nx_total)
               column_total = total_h(i, j) + eta(i, j)
               target_h(i, j, k) = column_total*dsig(k)
            end do

         case (VCOORD_ZSIGMA)
            ! Smoothstep blend: sigma in shallow, fixed z-levels in deep.
            ! Mirrors `vcoord_target_dz_column` in `src/ALE/rdb_vcoord.F90`
            ! but built directly on the 2D (H, η) fields.  The deep-branch
            ! z-level intervals are clipped to the local column total so
            ! sum_k target_h = H + η exactly even when the column is shallower
            ! than the deepest reference interface; any residual deficit is
            ! deposited in the bed-side layer (k=1) to preserve the sum.
            do concurrent(j=1:ny_total, i=1:nx_total) &
               local(column_total, alpha, x, k, z_top_k, z_bot_k, dz_z, dz_sum, deficit)
               column_total = total_h(i, j) + eta(i, j)
               if (column_total <= zsigma_depth_transition) then
                  do k = 1, nz
                     target_h(i, j, k) = dsig(k)*column_total
                  end do
               else
                  if (zsigma_blend_width > 0.0_wp) then
                     x = (column_total - zsigma_depth_transition)/zsigma_blend_width
                     x = max(0.0_wp, min(1.0_wp, x))
                     alpha = x*x*(3.0_wp - 2.0_wp*x)
                  else
                     alpha = 1.0_wp
                  end if
                  dz_sum = 0.0_wp
                  do k = 1, nz
                     z_top_k = min(z_ref_global(nz - k), column_total)
                     z_bot_k = min(z_ref_global(nz - k + 1), column_total)
                     dz_z = max(z_bot_k - z_top_k, 0.0_wp)
                     target_h(i, j, k) = (1.0_wp - alpha)*dsig(k)*column_total &
                                         + alpha*dz_z
                     dz_sum = dz_sum + target_h(i, j, k)
                  end do
                  deficit = column_total - dz_sum
                  target_h(i, j, 1) = target_h(i, j, 1) + deficit
               end if
            end do

         case (VCOORD_ZSTAR_SIGMA)
            z_ref_nz_inv = 0.0_wp
            if (z_ref_global(nz) > 0.0_wp) then
               z_ref_nz_inv = 1.0_wp/z_ref_global(nz)
            end if
            do concurrent(j=1:ny_total, i=1:nx_total) &
               local(column_total, alpha, x, k, dz_z)
               column_total = total_h(i, j) + eta(i, j)
               if (column_total <= zsigma_depth_transition .or. z_ref_nz_inv == 0.0_wp) then
                  do k = 1, nz
                     target_h(i, j, k) = dsig(k)*column_total
                  end do
               else
                  if (zsigma_blend_width > 0.0_wp) then
                     x = (column_total - zsigma_depth_transition)/zsigma_blend_width
                     x = max(0.0_wp, min(1.0_wp, x))
                     alpha = x*x*(3.0_wp - 2.0_wp*x)
                  else
                     alpha = 1.0_wp
                  end if
                  do k = 1, nz
                     dz_z = (z_ref_global(nz - k + 1) - z_ref_global(nz - k)) &
                            *column_total*z_ref_nz_inv
                     target_h(i, j, k) = (1.0_wp - alpha)*dsig(k)*column_total &
                                         + alpha*dz_z
                  end do
               end if
            end do

         case (VCOORD_ZSTAR_FULL)
            ! Per-column z*-full: walk the cached `z_ref(i, j, 0:nz)` from
            ! `build_zref_full`.  Surface layer absorbs η when η ≥ 0;
            ! bed-side layers vanish to `zstar_h_min` and the surface gets
            ! trimmed for exact conservation when η < 0.  Mirrors
            ! `vcoord_target_dz_column_zstar_full` (coastal) but emits
            ! ROMS order (k=1 bed, k=nz surface) directly.
            do concurrent(j=1:ny_total, i=1:nx_total) &
               local(k, h_bed_ref, eta_loc, H_eff, z_upper, z_lower, sum_dz, deficit)
               h_bed_ref = z_ref(i, j, nz)
               eta_loc = (total_h(i, j) + eta(i, j)) - h_bed_ref
               H_eff = max(total_h(i, j) + eta(i, j), 0.0_wp)
               if (h_bed_ref <= 0.0_wp) then
                  ! Degenerate column: emit a single vanishing-layer stack.
                  do k = 1, nz
                     target_h(i, j, k) = zstar_h_min
                  end do
               else if (eta_loc >= 0.0_wp) then
                  ! Column at or above reference: subsurface = z_ref intervals,
                  ! surface (k=nz) gets the +η.
                  do k = 1, nz
                     ! k_top = nz - k + 1 in the top-down z_ref convention.
                     target_h(i, j, k) = max( &
                                         z_ref(i, j, nz - k + 1) - z_ref(i, j, nz - k), &
                                         0.0_wp)
                  end do
                  target_h(i, j, nz) = target_h(i, j, nz) + eta_loc
               else
                  ! Column shallower than reference (η < 0).  Walk top-down,
                  ! clip layers to H_eff, vanish below.
                  do k = 1, nz
                     z_upper = z_ref(i, j, nz - k)
                     z_lower = z_ref(i, j, nz - k + 1)
                     if (z_lower <= H_eff) then
                        target_h(i, j, k) = z_lower - z_upper
                     else if (z_upper < H_eff) then
                        target_h(i, j, k) = H_eff - z_upper
                     else
                        target_h(i, j, k) = zstar_h_min
                     end if
                  end do
                  ! Surface trim: drop the vanishing-layer overhead from the
                  ! surface to make sum = H exactly.  If the surface would
                  ! itself fall below h_min, leave it at h_min and let the
                  ! downstream dry-cell guards handle the deficit.
                  sum_dz = 0.0_wp
                  do k = 1, nz
                     sum_dz = sum_dz + target_h(i, j, k)
                  end do
                  deficit = sum_dz - H_eff
                  if (deficit > 0.0_wp) then
                     if (target_h(i, j, nz) - deficit >= zstar_h_min) then
                        target_h(i, j, nz) = target_h(i, j, nz) - deficit
                     else
                        target_h(i, j, nz) = zstar_h_min
                     end if
                  end if
               end if
            end do

         case (VCOORD_Z_FIXED)
            ! Fixed-z interfaces with vanishing layers in shallow water.
            ! Per-column algorithm mirrors `seed_h_layer_z_fixed_impl`:
            ! walk surface-down, each layer takes its nominal thickness
            ! `h_nominal = h_ref / nz_ml` if the interface below it sits
            ! inside the remaining water column; else collapses to
            ! `h_min` while the surface absorbs the residual.  When
            ! `z_fixed_h_ref = 0` (knob unset) fall back to uniform
            ! `H · dsig(k)` so tests that omit the knob still get
            ! something sensible.
            !
            ! `total_h(i,j) + eta(i,j)` is the live column total; ALE
            ! remaps toward these per-step targets every outer step,
            ! keeping the interface anchored at `z_target(k) = k · h_nominal`.
            h_nominal = 0.0_wp
            if (z_fixed_h_ref > 0.0_wp) then
               h_nominal = z_fixed_h_ref/real(nz, wp)
            end if
            h_min = zstar_h_min
            if (h_nominal <= 0.0_wp) then
               ! Uniform-sigma fallback — keeps the path active when
               ! the knob isn't set (matches sigma behaviour).
               do concurrent(k=1:nz, j=1:ny_total, i=1:nx_total)
                  target_h(i, j, k) = (total_h(i, j) + eta(i, j))*dsig(k)
               end do
            else
               do concurrent(j=1:ny_total, i=1:nx_total) &
                  local(k, z_below_loc, z_above_nominal_loc)
                  z_below_loc = total_h(i, j) + eta(i, j)
                  do k = 1, nz
                     z_above_nominal_loc = real(nz - k, wp)*h_nominal
                     if (z_above_nominal_loc > z_below_loc - h_min) then
                        ! Above the column / would be sub-h_min — vanish.
                        target_h(i, j, k) = h_min
                        z_below_loc = z_below_loc - h_min
                     else
                        ! Layer fits — take nominal thickness (or surface
                        ! residual for k=nz where z_above_nominal=0).
                        target_h(i, j, k) = z_below_loc - z_above_nominal_loc
                        z_below_loc = z_above_nominal_loc
                     end if
                  end do
               end do
            end if

         case default
            error stop "ocean_vcoord_compute_target_h: unknown coord_type."
         end select
      end associate
   end subroutine ocean_vcoord_compute_target_h_impl

   pure subroutine ocean_vcoord_compute_target_h_rho(this, total_h, eta, T, S, eos, hybrid)
      !! Thin polymorphic wrapper for the isopycnal (`VCOORD_RHO`) and
      !! hybrid z*/isopycnal (`VCOORD_HYCOM`) target-grid build.  Mirrors
      !! `ocean_vcoord_compute_target_h` (resolves the concrete type so
      !! the `do concurrent` kernel runs on `type(ocean_vcoord_t)`, never
      !! a polymorphic list item).  Separate from `compute_target_h`
      !! because the RHO/HYCOM branch needs the per-layer T/S
      !! concentrations and the device-resident EOS coefficients, which
      !! the shared `pure (total_h, eta)` TBP cannot carry.
      !!
      !! `hybrid` selects the variant: `.false.` (default) = pure
      !! `VCOORD_RHO` (bit-identical with the P2 kernel); `.true.` =
      !! `VCOORD_HYCOM` (adds the bottom-up density monotonize and the
      !! z* nominal-floor sweep inside the same column kernel).
      class(ocean_vcoord_t), intent(inout) :: this
      real(wp), intent(in) :: total_h(:, :)  ! assumed-shape-ok: thin TBP wrapper, no kernel
      real(wp), intent(in) :: eta(:, :)      ! assumed-shape-ok: thin TBP wrapper, no kernel
      real(wp), intent(in) :: T(:, :, :)     ! assumed-shape-ok: thin TBP wrapper, no kernel
      real(wp), intent(in) :: S(:, :, :)     ! assumed-shape-ok: thin TBP wrapper, no kernel
      type(eos_t), intent(in) :: eos
      logical, intent(in), optional :: hybrid
         !! Enable the HYCOM hybrid deltas (default `.false.` = pure RHO).
      logical :: hybrid_loc
      hybrid_loc = .false.
      if (present(hybrid)) hybrid_loc = hybrid
      select type (this)
      type is (ocean_vcoord_t)
         call ocean_vcoord_compute_target_h_rho_impl(this, total_h, eta, T, S, eos, hybrid_loc)
      end select
   end subroutine ocean_vcoord_compute_target_h_rho

   pure subroutine ocean_vcoord_compute_target_h_rho_impl(this, total_h, eta, T, S, eos, hybrid)
      !! Isopycnal regrid: place layer interfaces on the prescribed
      !! `rho_target(0:nz)` potential-density surfaces.  ONE
      !! `do concurrent(j,i)` over columns, each running the full
      !! density-space inversion on `NZ_STACK_MAX` fixed-size locals —
      !! no host loop, no per-call allocate.
      !!
      !! Algorithm (clean-room from Bleck 2002 / White & Adcroft 2008;
      !! MOM6 `coord_rho` is the behavioural oracle), per the spec:
      !!
      !!   0. Pre-compaction — strip source layers `h <= H_MIN`, donate
      !!      their volume to the thickest survivor (volume-conserving),
      !!      build the survivor count.  Fast path: `<= 1` survivor →
      !!      `h_new = h_old` (no inversion).
      !!   1. Layer potential densities via `eos_density_point` at
      !!      `rho_ref_pressure` on the compacted column.
      !!   2. PPM (Colella-Woodward monotone) reconstruction of the
      !!      density profile over the compacted thicknesses.
      !!   3. Per interior target, bracket-ordered inversion: light
      !!      boundary → surface; discontinuous-jump sweep; dense
      !!      boundary → bed; else fixed-8-iter Newton on `xi in [0,1]`
      !!      (convergence on `|delta| < NR_TOL` AFTER `xi += delta`;
      !!      zero-gradient `NR_OFFSET` escape at both ends; masked
      !!      fallback to the previous interface on no-bracket).
      !!   4. Monotone non-decreasing interfaces → `h_new`.
      !!   5. MOM6 min-thickness inflation, floor = `max(zstar_h_min,
      !!      H_VANISHED)` (MUST-HAVE #2: keeps RHO-collapsed layers
      !!      above the remap-drain `H_FLOOR` so the next regrid does
      !!      not zero their tracer mass); debit the single thickest
      !!      layer once.
      !!
      !! Internal working frame is TOP-DOWN (index 1 = surface, +down),
      !! matching the validated prototype and the `rho_target(0)` =
      !! lightest = surface convention.  The final assignment FLIPS to
      !! the bottom-up state (MUST-HAVE #3): `target_h(:,:,k) =
      !! h_new_td(nz - k + 1)`, so the lightest target lands at the
      !! surface (k=nz) and the densest at the bed (k=1).  Sum is
      !! conserved exactly so `sum_k target_h = H` (η is implicit in
      !! `total_h` here — the caller passes the live column total as
      !! `total_h`, see `ocean_apply_ale_remap_step`).
      !!
      !! HYCOM hybrid (`hybrid = .true.`, Bleck 2002 / MOM6
      !! `build_hycom1_column`): two deltas around the unchanged RHO
      !! inversion, both inside this same column kernel.
      !!   (1b) BOTTOM-UP density monotonize before the PPM reconstruction:
      !!        in the top-down work frame, cap each cell by the one below
      !!        it (toward the bed) — `do k=nk-1,1,-1: rhoc(k)=min(rhoc(k),
      !!        rhoc(k+1))`.  Pure RHO omits this (assumes a monotone
      !!        profile + leans on the PPM limiter); HYCOM enforces it so
      !!        the inversion always sees a monotone column.
      !!   (4b) z* NOMINAL-FLOOR sweep after the inversion, before the
      !!        monotone/inflation: walk interior+bottom interfaces from
      !!        the surface down, accumulating `dsig·total_h·stretching`
      !!        (= `dsig·(H+η)`, since `stretching = (H+η)/H`) and pushing
      !!        each interface DOWN to at least that nominal z* depth
      !!        (clamped to the column bottom).  This is a surface-side
      !!        minimum-depth floor: it protects the near-surface band
      !!        from collapse (fixed z* resolution) while leaving deep
      !!        isopycnal interfaces — already below the floor — untouched.
      !!        CRITICAL: the accumulation factor is `total_h` (= H), not
      !!        `H+η`; `dsig` sums to 1, so `dsig·H·stretching = dsig·(H+η)`.
      !!        Writing `dsig·(H+η)·stretching` over-stretches by `(H+η)/H`
      !!        — invisible at η=0, wrong with a free surface.  No
      !!        renormalize after the floor sweep (would break the floor
      !!        invariant `z(k) ≥ Σ_{j≤k} dsig·(H+η)`; MOM6 pins the bottom
      !!        interface + uses the debit-thickest inflation instead).
      !! `hybrid = .false.` (the `VCOORD_RHO` path) skips BOTH deltas and
      !! is bit-identical to the P2 kernel.
      type(ocean_vcoord_t), intent(inout) :: this
      ! assumed-shape-ok: cadence-bounded (once per outer ALE step); the per-
      ! column inversion below uses fixed-size NZ_STACK_MAX locals only.
      real(wp), intent(in) :: total_h(:, :)
         !! Column reference depth H(i, j) (m).  Caller passes the live
         !! column total (sum of h_layer) so the new grid spans it exactly.
      real(wp), intent(in) :: eta(:, :)  ! assumed-shape-ok: see total_h
         !! Free-surface anomaly η(i, j) (m).  Added to `total_h` to form
         !! the column extent the new interfaces span.
      real(wp), intent(in) :: T(:, :, :)  ! assumed-shape-ok: see total_h
         !! Layer-mean potential temperature concentration (°C), `(nx,ny,nz)`.
      real(wp), intent(in) :: S(:, :, :)  ! assumed-shape-ok: see total_h
         !! Layer-mean salinity concentration (PSU), `(nx,ny,nz)`.
      type(eos_t), intent(in) :: eos
         !! Shared device-resident EOS handle (flat POD, by value).
      logical, intent(in) :: hybrid
         !! `.true.` = HYCOM (apply the monotonize + z*-floor deltas);
         !! `.false.` = pure RHO (bit-identical with the P2 kernel).

      integer :: i, j, nz
      integer :: k, kk, nk, ns, idx_thick, src, ii
      integer :: mapping(NZ_STACK_MAX)
      real(wp) :: h_col(NZ_STACK_MAX), t_col(NZ_STACK_MAX), s_col(NZ_STACK_MAX)
      real(wp) :: hc(NZ_STACK_MAX), rhoc(NZ_STACK_MAX), rtgt(NZ_STACK_MAX)
      real(wp) :: z_new(NZ_STACK_MAX + 1), h_new(NZ_STACK_MAX)
      real(wp) :: col_extent, h_floor_eff, h_min, donate
      real(wp) :: total_need, thick_max
      real(wp) :: nominal_z, stretching, h_ref_col

      if (.not. this%is_init) return
      nz = this%nz_ml
      h_min = this%zstar_h_min
      ! Inflation floor must be STRICTLY above H_VANISHED: the remap drain
      ! (`ocean_remap_tracer_field`) gates on `h_old > H_FLOOR` (== H_VANISHED)
      ! with a strict `>`, so a layer sitting exactly at H_VANISHED has its
      ! tracer concentration zeroed when this column is fed back as `h_old`
      ! on the next regrid.  Floor at 2·H_VANISHED so inflated layers always
      ! survive the drain (closes the multi-regrid mass-loss footgun).
      h_floor_eff = max(this%zstar_h_min, 2.0_wp*H_VANISHED)

      ! One column per (j,i).  NZ_STACK_MAX fixed-size locals; no name
      ! shadows a Fortran intrinsic; cross-module pure EOS helper carries
      ! its own `!$acc routine seq`.
      ! associate the components the kernel touches to plain names — ifx's
      ! do-concurrent -> OpenMP-target lowering ICEs on a `this%<component>`
      ! reference inside the loop body (see compute_target_h_impl above).
      associate (target_h => this%target_h, dsig => this%dsig, &
                 rho_target => this%rho_target, remap_h_old => this%remap_h_old, &
                 rho_ref_pressure => this%rho_ref_pressure, &
                 nx_total => this%nx_total, ny_total => this%ny_total)
      do concurrent(j=1:ny_total, i=1:nx_total) &
         local(k, kk, ii, nk, ns, idx_thick, src, mapping, &
               h_col, t_col, s_col, hc, rhoc, rtgt, &
               z_new, h_new, col_extent, donate, &
               total_need, thick_max, nominal_z, stretching, h_ref_col)

         ! --- gather TOP-DOWN: working index 1 = surface = state k=nz ---
         ! Source thicknesses come from the remap snapshot the
         ! orchestrator placed in `remap_h_old` (the live, pre-remap
         ! `h_layer`); T/S are the layer-mean concentrations.  Flip the
         ! bottom-up state index (k=nz surface) into the top-down work
         ! frame (work index 1 = surface).
         col_extent = max(total_h(i, j) + eta(i, j), 0.0_wp)
         do k = 1, nz
            ii = nz - k + 1                 ! state (bottom-up) index
            h_col(k) = remap_h_old(i, j, ii)
            t_col(k) = T(i, j, ii)
            s_col(k) = S(i, j, ii)
         end do

         ! --- step 0: pre-compaction (strip h <= h_min, donate) ---
         nk = 0
         do k = 1, nz
            if (h_col(k) > h_min) then
               nk = nk + 1
               mapping(nk) = k
               hc(nk) = h_col(k)
            end if
         end do
         if (nk <= 1) then
            ! Fast path: <= 1 finite layer.  nz == nk_state here, so
            ! keep the source thicknesses unchanged (h_new = h_old),
            ! flipped back into the bottom-up state.
            do k = 1, nz
               target_h(i, j, k) = remap_h_old(i, j, k)
            end do
            cycle
         end if

         ! Donate the stripped volume to the thickest survivor.
         donate = col_extent
         do kk = 1, nk
            donate = donate - hc(kk)
         end do
         if (donate > 0.0_wp) then
            idx_thick = 1
            do kk = 2, nk
               if (hc(kk) > hc(idx_thick)) idx_thick = kk
            end do
            hc(idx_thick) = hc(idx_thick) + donate
         end if

         ! --- step 1: layer potential densities on the compacted column ---
         do kk = 1, nk
            src = mapping(kk)
            rhoc(kk) = eos_density_point(eos, t_col(src), s_col(src), &
                                         rho_ref_pressure)
         end do

         ! --- step 1b (HYCOM only): bottom-up density monotonize ---
         ! Work frame is top-down (index 1 = surface, nk = bed): cap each
         ! cell by the one below it sweeping bed-up so density is
         ! non-decreasing downward.  Pure RHO (hybrid=.false.) skips this
         ! and stays bit-identical to the merged RHO regrid.
         if (hybrid) then
            do kk = nk - 1, 1, -1
               rhoc(kk) = min(rhoc(kk), rhoc(kk + 1))
            end do
         end if

         ! --- steps 2-4: PPM reconstruct + invert each interior target
         !     density to an interface depth + monotone interfaces.  Shared
         !     density-space inversion (also used by the DENSITY diagnostic
         !     remap, `rdb_ocean_diag_fills`) — single source of truth for
         !     the bracket + fixed-iter-Newton solve.  Copy the interior
         !     targets into a stack array so the device call passes a whole
         !     fixed-size local (no derived-type section descriptor in the
         !     hot per-column kernel).  For HYCOM the rhoc fed in was
         !     monotonized above; the z* floor below then lifts the result.
         do kk = 1, nz - 1
            rtgt(kk) = rho_target(kk)
         end do
         call invert_density_targets(nk, hc, rhoc, nz - 1, rtgt, z_new)

         ! --- step 4b (HYCOM only): z* nominal-floor sweep ---
         ! Surface-side minimum-depth floor on the isopycnal interfaces.
         ! Walk interfaces from the surface down, accumulating the nominal
         ! z* depth and pushing any too-shallow interface DOWN to it
         ! (clamped to the column bottom); deep interfaces already below the
         ! floor are untouched.  stretching = col_extent/total_h (the
         ! SSH-following z* stretch); the accumulation is
         ! dsig*total_h*stretching (= dsig*col_extent), and dsig sums to 1.
         ! dsig is stored bottom-up (dsig(nz) = surface layer); the work
         ! layer above interface kk maps to bottom-up dsig index nz-kk+2.
         ! The floor can break monotonicity, so re-monotonize after it.
         ! Pure RHO (hybrid=.false.) skips this and is bit-identical.
         if (hybrid) then
            h_ref_col = total_h(i, j)
            if (h_ref_col > 0.0_wp) then
               stretching = col_extent/h_ref_col
            else
               stretching = 1.0_wp
            end if
            nominal_z = 0.0_wp
            do kk = 2, nz + 1
               nominal_z = nominal_z + dsig(nz - kk + 2)*h_ref_col*stretching
               if (z_new(kk) < nominal_z) z_new(kk) = nominal_z
               if (z_new(kk) > col_extent) z_new(kk) = col_extent
            end do
            do kk = 2, nz + 1
               if (z_new(kk) < z_new(kk - 1)) z_new(kk) = z_new(kk - 1)
            end do
         end if
         do kk = 1, nz
            h_new(kk) = z_new(kk + 1) - z_new(kk)
         end do

         ! --- step 5: MOM6 min-thickness inflation (floor h_floor_eff) ---
         ns = 0
         do kk = 1, nz
            if (h_new(kk) > h_floor_eff) ns = ns + 1
         end do
         if (ns == nz) then
            ! all OK
         else if (ns == 0) then
            do kk = 1, nz
               h_new(kk) = h_floor_eff
            end do
         else
            total_need = 0.0_wp
            do kk = 1, nz
               if (h_new(kk) <= h_floor_eff) then
                  total_need = total_need + (h_floor_eff - h_new(kk))
                  h_new(kk) = h_floor_eff
               end if
            end do
            ! debit the single thickest layer once
            idx_thick = 1
            thick_max = h_new(1)
            do kk = 2, nz
               if (h_new(kk) > thick_max) then
                  thick_max = h_new(kk)
                  idx_thick = kk
               end if
            end do
            h_new(idx_thick) = h_new(idx_thick) - total_need
         end if

         ! --- assignment: FLIP top-down working -> bottom-up state ---
         do k = 1, nz
            target_h(i, j, k) = h_new(nz - k + 1)
         end do
      end do
      end associate
   end subroutine ocean_vcoord_compute_target_h_rho_impl

   pure function parse_ocean_vcoord_type(name) result(code)
      !! Ocean-path wrapper around the canonical `parse_vcoord_type`
      !! in `rdb_vcoord`.  Pins the unrecognised-string fallback to
      !! `VCOORD_EULERIAN_Z` — the ocean path's "do nothing" default,
      !! distinct from the coastal path's `VCOORD_SIGMA` fallback.
      !! Kept as a thin name-preserving wrapper so the ocean-only
      !! semantic (fallback choice) is visible at the call site.
      character(len=*), intent(in) :: name
      integer :: code
      code = parse_vcoord_type(name, default_code=VCOORD_EULERIAN_Z)
   end function parse_ocean_vcoord_type

   pure subroutine invert_density_targets(nk, hc, rhoc, n_int, rho_tgt, z_new)
      !! Density-space interface inversion — the single source of truth
      !! for the RHO vcoord regrid (`ocean_vcoord_compute_target_h_rho_impl`)
      !! AND the DENSITY diagnostic remap (`rdb_ocean_diag_fills`).
      !!
      !! Given a TOP-DOWN compacted column of `nk` layers with thicknesses
      !! `hc(1:nk)` (index 1 = surface, depth positive-down) and layer-mean
      !! potential densities `rhoc(1:nk)`, plus `n_int` monotone-increasing
      !! interior target densities `rho_tgt(1:n_int)`, return the `n_int + 2`
      !! interface depths `z_new(1:n_int+2)` (top-down, 0 .. column total),
      !! monotone non-decreasing.  `z_new(1) = 0` (surface), `z_new(n_int+2)`
      !! = Σ hc (bed); interior interfaces 2..n_int+1 invert each target.
      !!
      !! Algorithm (clean-room from Bleck 2002 / White & Adcroft 2008;
      !! MOM6 `coord_rho` is the behavioural oracle):
      !!   1. PPM (Colella-Woodward monotone) reconstruction of the density
      !!      profile over `hc`.
      !!   2. Per interior target: light boundary → surface; discontinuous-
      !!      jump sweep; dense boundary → bed; else fixed-`NR_ITERS`-iter
      !!      Newton on `xi in [0,1]` (convergence on `|delta| < NR_TOL`
      !!      AFTER `xi += delta`; zero-gradient `NR_OFFSET` escape at both
      !!      ends; masked fallback to the previous interface on no-bracket).
      !!   3. Monotone non-decreasing clamp on the interfaces.
      !!
      !! Caller supplies `nk >= 2`.  The RHO regrid kernel pre-compacts
      !! vanished layers (so `nk` is the surviving count) and fast-paths
      !! `nk <= 1` upstream; the DENSITY diagnostic remap passes the full
      !! `nz` column (it assumes a non-vanished column — vanished-layer
      !! compaction for diagnostics is a deferred refinement).  Lightest
      !! target maps to the surface (index 2), densest to the bed (the
      !! surface→bed ordering the callers FLIP into the bottom-up state).
      !$acc routine seq
      integer, intent(in)  :: nk
      integer, intent(in)  :: n_int
      real(wp), intent(in)  :: hc(nk)
      real(wp), intent(in)  :: rhoc(nk)
      real(wp), intent(in)  :: rho_tgt(n_int)
      real(wp), intent(out) :: z_new(n_int + 2)
      integer  :: kk, ii, k
      real(wp) :: rhoL(NZ_STACK_MAX), rhoR(NZ_STACK_MAX)
      real(wp) :: edge(NZ_STACK_MAX + 1), z_old(NZ_STACK_MAX + 1)
      real(wp) :: tgt, lo, hi, q6, xi, fval, df, delta, grad, ww
      real(wp) :: rho_light, rho_dense, dd, six
      logical  :: placed

      ! --- step 1: PPM (Colella-Woodward) edges on the compacted column ---
      edge(1) = rhoc(1)
      edge(nk + 1) = rhoc(nk)
      do kk = 2, nk
         ww = hc(kk - 1) + hc(kk)
         ! Floor guards an uncompacted caller (the DENSITY diagnostic remap
         ! passes the raw column) where a vanished layer pair sums to ~0;
         ! the RHO regrid pre-compacts so ww is always >> the floor there
         ! (this branch is inert for it — bit-identical).
         if (ww < 1.0e-30_wp) ww = 1.0e-30_wp
         edge(kk) = (hc(kk)*rhoc(kk - 1) + hc(kk - 1)*rhoc(kk))/ww
      end do
      do kk = 1, nk
         rhoL(kk) = edge(kk)
         rhoR(kk) = edge(kk + 1)
         ! CW monotonic limiter
         if ((rhoR(kk) - rhoc(kk))*(rhoc(kk) - rhoL(kk)) <= 0.0_wp) then
            rhoL(kk) = rhoc(kk)
            rhoR(kk) = rhoc(kk)
         else
            dd = rhoR(kk) - rhoL(kk)
            six = 6.0_wp*(rhoc(kk) - 0.5_wp*(rhoL(kk) + rhoR(kk)))
            if (dd*six > dd*dd) then
               rhoL(kk) = 3.0_wp*rhoc(kk) - 2.0_wp*rhoR(kk)
            else if (dd*six < -dd*dd) then
               rhoR(kk) = 3.0_wp*rhoc(kk) - 2.0_wp*rhoL(kk)
            end if
         end if
      end do

      ! Cumulative compacted-grid interface depths (top-down, 0..H).
      z_old(1) = 0.0_wp
      do kk = 1, nk
         z_old(kk + 1) = z_old(kk) + hc(kk)
      end do
      rho_light = rhoL(1)
      rho_dense = rhoR(nk)

      ! --- step 2: invert each interior target interface (top-down) ---
      z_new(1) = 0.0_wp
      z_new(n_int + 2) = z_old(nk + 1)     ! total compacted depth (= col extent)
      do kk = 2, n_int + 1                  ! interior interfaces
         tgt = rho_tgt(kk - 1)
         if (tgt <= rho_light) then
            z_new(kk) = 0.0_wp              ! lighter than column -> surface
         else if (tgt >= rho_dense) then
            z_new(kk) = z_old(nk + 1)       ! denser than column -> bed
         else
            placed = .false.
            do ii = 1, nk
               lo = rhoL(ii)
               hi = rhoR(ii)
               ! Discontinuous jump at the TOP interface of cell ii
               ! (between cell ii-1's right edge and cell ii's left
               ! edge).  Checked FIRST and independent of whether the
               ! cells are limiter-flattened — a 2-layer column flattens
               ! both cells, and a target inside the jump must still
               ! land on the interface (MOM6 coord_rho behaviour).
               if (ii > 1) then
                  if (rhoR(ii - 1) <= tgt .and. tgt <= rhoL(ii)) then
                     z_new(kk) = z_old(ii)
                     placed = .true.
                     exit
                  end if
               end if
               if (lo == hi) then
                  ! Flat (limiter-collapsed) cell: only an exact match
                  ! places here; otherwise advance to the next cell.
                  if (abs(tgt - lo) < NR_TOL) then
                     z_new(kk) = z_old(ii)
                     placed = .true.
                     exit
                  end if
                  cycle
               end if
               if ((lo - tgt)*(hi - tgt) <= 0.0_wp) then
                  ! Newton on xi in [0,1]: ppm(xi) - tgt = 0
                  q6 = 6.0_wp*rhoc(ii) - 3.0_wp*(lo + hi)
                  xi = 0.5_wp
                  do k = 1, NR_ITERS
                     fval = lo + xi*((hi - lo) + q6*(1.0_wp - xi)) - tgt
                     df = (hi - lo) + q6*(1.0_wp - 2.0_wp*xi)
                     if (abs(df) > 1.0e-30_wp) then
                        delta = -fval/df
                     else
                        delta = 0.0_wp
                     end if
                     xi = xi + delta
                     ! clamp inside the iteration; zero-gradient nudge
                     if (xi < 0.0_wp) then
                        xi = 0.0_wp
                        grad = (hi - lo) + q6           ! d(ppm)/dxi at xi=0
                        if (abs(grad) < 1.0e-30_wp) xi = NR_OFFSET
                     else if (xi > 1.0_wp) then
                        xi = 1.0_wp
                        grad = (hi - lo) - q6           ! d(ppm)/dxi at xi=1
                        if (abs(grad) < 1.0e-30_wp) xi = 1.0_wp - NR_OFFSET
                     end if
                     if (abs(delta) < NR_TOL) exit
                  end do
                  z_new(kk) = z_old(ii) + xi*hc(ii)
                  placed = .true.
                  exit
               end if
            end do
            if (.not. placed) then
               ! masked fallback: previous interface (no FATAL in a DC)
               z_new(kk) = z_new(kk - 1)
            end if
         end if
      end do

      ! --- step 3: monotone non-decreasing interfaces ---
      do kk = 2, n_int + 2
         if (z_new(kk) < z_new(kk - 1)) z_new(kk) = z_new(kk - 1)
      end do
   end subroutine invert_density_targets

   pure function ocean_vcoord_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the vertical coordinate slot
      !! (0 when unallocated). One arr_bytes term per array — add a
      !! term here when a new allocatable joins the type.
      class(ocean_vcoord_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = arr_bytes(this%dsig) &
               + arr_bytes(this%z_ref_global) &
               + arr_bytes(this%target_h) &
               + arr_bytes(this%z_ref) &
               + arr_bytes(this%rho_target) &
               + arr_bytes(this%remap_total_h) &
               + arr_bytes(this%remap_h_ref) &
               + arr_bytes(this%remap_h_old) &
               + arr_bytes(this%remap_conc_t) &
               + arr_bytes(this%remap_conc_s)
   end function ocean_vcoord_bytes

end module rdb_ocean_vcoord
